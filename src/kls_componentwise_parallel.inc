/* Optional row-owned execution of the same accumulated-error proof. This
 * module caches only CSC traversal positions, never values or certificates. */
#include "kls_tuning.inc"
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

struct kls_componentwise_plan {
  size_t n;
  int64_t nnz;
  size_t *rowptr;
  int64_t *positions,*columns;
  int eligible,threads;
  int verdict_only; /* Output contract, never a cached answer. */
  kls_componentwise_executor executor;
  void *executor_context;
  struct {
    int samples,choice; /* 0 measuring, -1 serial rows, +1 parallel rows */
    double ratios[KLS_CERTIFICATE_PROBE_PAIRS];
  } direction[2];
};

void kls_componentwise_plan_reset(kls_componentwise_plan *plan) {
  if(!plan) return;
  free(plan->rowptr);free(plan->positions);free(plan->columns);
  *plan=(kls_componentwise_plan){0};
}

void kls_componentwise_plan_destroy(kls_componentwise_plan *plan) {
  kls_componentwise_plan_reset(plan);free(plan);
}

#if KLS_COMPONENTWISE_HAS_ACCUMULATED && defined(_OPENMP)
static int componentwise_plan_environment(void) {
  unsigned short control;
  __asm__ volatile("fnstcw %0":"=m"(control));
  return fegetround()==FE_TONEAREST && (control&0x0f00u)==0x0300u;
}

static void componentwise_plan_build(kls_componentwise_plan *plan,size_t n,
    const int64_t *p,const int64_t *rows) {
  size_t *cursor=NULL;
  plan->n=n;plan->nnz=p[n];plan->eligible=-1;
  if((uintmax_t)p[n]>SIZE_MAX/sizeof(int64_t)) return;
  const size_t count=p[n]?(size_t)p[n]:1;
  plan->rowptr=calloc(n+1,sizeof(*plan->rowptr));
  plan->positions=malloc(count*sizeof(*plan->positions));
  plan->columns=malloc(count*sizeof(*plan->columns));
  cursor=malloc(n*sizeof(*cursor));
  if(!plan->rowptr || !plan->positions || !plan->columns || !cursor) goto declined;
  for(size_t i=0;i<n;++i) cursor[i]=SIZE_MAX;
  for(size_t j=0;j<n;++j) {
    if(p[j]<0 || p[j]>p[j+1] || p[j+1]>p[n]) goto declined;
    for(int64_t k=p[j];k<p[j+1];++k) {
      if(rows[k]<0 || (uintmax_t)rows[k]>=n) goto declined;
      const size_t i=(size_t)rows[k];
      /* Coalesced magnitudes require the interval fallback for duplicates. */
      if(cursor[i]==j) goto declined;
      cursor[i]=j;++plan->rowptr[i+1];
    }
  }
  for(size_t i=0;i<n;++i) {plan->rowptr[i+1]+=plan->rowptr[i];cursor[i]=plan->rowptr[i];}
  /* Stable scatter preserves exactly the original CSC encounter order. */
  for(size_t j=0;j<n;++j) for(int64_t k=p[j];k<p[j+1];++k) {
    size_t at=cursor[rows[k]]++;
    plan->positions[at]=k;plan->columns[at]=(int64_t)j;
  }
  free(cursor);plan->eligible=1;return;
declined:
  free(cursor);kls_componentwise_plan_reset(plan);
  plan->n=n;plan->nnz=p[n];plan->eligible=-1;
}

/* This is the same per-row bound as the scatter proof above. Serial and
 * parallel execution call this same function in the same per-row order. */
static int componentwise_binary64_normal_or_zero(double value) {
  uint64_t bits;
  memcpy(&bits,&value,sizeof(bits));
  bits&=UINT64_C(0x7fffffffffffffff);
  return bits==0 || (bits>=UINT64_C(0x0010000000000000) &&
                    bits<UINT64_C(0x7ff0000000000000));
}
/* Nonnegative IEEE binary64 values only. One adjacent representable value
 * encloses an operation rounded to nearest, including subnormal results.
 * Zero is a valid lower bound for a nonnegative exact result. */
static double componentwise_positive_up(double value) {
  uint64_t bits;memcpy(&bits,&value,sizeof(bits));
  if(bits<UINT64_C(0x7ff0000000000000)) ++bits;
  memcpy(&value,&bits,sizeof(bits));return value;
}
static double componentwise_positive_down(double value) {
  uint64_t bits;memcpy(&bits,&value,sizeof(bits));
  if(bits) --bits;
  memcpy(&value,&bits,sizeof(bits));return value;
}
/* This x86-64 proof requires SSE2. Explicit scalar instructions retain one
 * binary64 rounding per operation without temporary-memory stores. The asm
 * boundaries prohibit contraction even when the including TU enables FMA. */
static inline double componentwise_binary64_mul(double left,double right) {
  __asm__ volatile("mulsd {%1,%0|%0,%1}":"+x"(left):"x"(right));return left;
}
static inline double componentwise_binary64_add(double left,double right) {
  __asm__ volatile("addsd {%1,%0|%0,%1}":"+x"(left):"x"(right));return left;
}
static inline double componentwise_binary64_sub(double left,double right) {
  __asm__ volatile("subsd {%1,%0|%0,%1}":"+x"(left):"x"(right));return left;
}
static inline double componentwise_binary64_div(double left,double right) {
  __asm__ volatile("divsd {%1,%0|%0,%1}":"+x"(left):"x"(right));return left;
}
/* Sufficient binary64 proof with an extended-precision fallback. Explicit
 * scalar instructions require separately rounded operations even with
 * contraction enabled. Exceptional ranges use the extended proof below. */
static int componentwise_binary64_row_value(const kls_componentwise_plan *plan,size_t row,
    const int64_t *p,const int64_t *rows,const double *a,const double *b,
    const double *x,int transpose,long double tolerance,double *upper) {
  /* Production callers are inside a row batch with nearest SSE rounding and
   * gradual underflow scoped by componentwise_plan_run or its callback. */
  double r=b[row],m=fabs(r);
  if(!componentwise_binary64_normal_or_zero(r) || !isfinite(x[row])) return 0;
  const size_t first=transpose?(size_t)p[row]:plan->rowptr[row];
  const size_t last=transpose?(size_t)p[row+1]:plan->rowptr[row+1];
  for(size_t k=first;k<last;++k) {
    const size_t at=transpose?k:(size_t)plan->positions[k];
    const size_t col=transpose?(size_t)rows[k]:(size_t)plan->columns[k];
    const double av=a[at],xv=x[col];
    if(!componentwise_binary64_normal_or_zero(av) ||
       !componentwise_binary64_normal_or_zero(xv)) return 0;
    const double product=componentwise_binary64_mul(av,xv);
    if(!componentwise_binary64_normal_or_zero(product) ||
       (product==0 && av!=0 && xv!=0)) return 0;
    const double next_r=componentwise_binary64_sub(r,product),
        next_m=componentwise_binary64_add(m,fabs(product));
    if(!componentwise_binary64_normal_or_zero(next_r) ||
       !componentwise_binary64_normal_or_zero(next_m)) return 0;
    r=next_r;m=next_m;
  }
  if(m==0) {*upper=0;return r==0;}
  /* Allow 2u per operation: a just-subnormal exact product can round up to
   * DBL_MIN. For at most 2k operations, t=4ku and gamma/(1-gamma)=t/(1-2t).
   * c=8(k+1)u<=1/4 ensures t<=1/8 and c dominates that error coefficient.
   * c<=1/4 also bounds k+1 by 2^48, so its integer conversion and the
   * power-of-two coefficient are exact in binary64. Outward binary64 bounds
   * enclose the same proof; exceptional ranges retain the extended fallback.
   * Explicit scalar instructions prevent contraction and excess intermediate
   * precision. */
  const double c=4*((double)(last-first)+1)*DBL_EPSILON;
  if(!(c<=0.25)) return 0;
  const double rounded_error=componentwise_binary64_mul(c,m);
  const double error=componentwise_positive_up(rounded_error);
  const double rounded_numerator=componentwise_binary64_add(fabs(r),error);
  const double numerator=componentwise_positive_up(rounded_numerator);
  if(!(m>error) || !isfinite(numerator)) return 0;
  const double rounded_denominator=componentwise_binary64_sub(m,error);
  const double denominator=componentwise_positive_down(rounded_denominator);
  const double tolerance64=componentwise_positive_down((double)tolerance);
  if(!(denominator>0)) return 0;
  const double rounded_threshold=componentwise_binary64_mul(tolerance64,denominator);
  if(numerator>componentwise_positive_down(rounded_threshold)) return 0;
  const double rounded_ratio=componentwise_binary64_div(numerator,denominator);
  const double ratio=componentwise_positive_up(rounded_ratio);
  if(!isfinite(ratio) || ratio>tolerance64) return 0;
  *upper=ratio;return 1;
}

/* Preserve the scalar test/row interface; complete binary64 batches convert
 * their maximum to extended precision only once. */
static int componentwise_binary64_row(const kls_componentwise_plan *plan,size_t row,
    const int64_t *p,const int64_t *rows,const double *a,const double *b,
    const double *x,int transpose,long double tolerance,long double *upper) {
  double bound;
  if(!componentwise_binary64_row_value(plan,row,p,rows,a,b,x,transpose,tolerance,&bound)) return 0;
  *upper=(long double)bound;return 1;
}

static int componentwise_plan_row(const kls_componentwise_plan *plan,size_t row,
    const int64_t *p,const int64_t *rows,const double *a,const double *b,
    const double *x,int transpose,int binary64_permitted,
    long double tolerance,long double *upper) {
  if(binary64_permitted &&
     componentwise_binary64_row(plan,row,p,rows,a,b,x,transpose,tolerance,upper)) return 1;
  volatile long double bj=(long double)b[row],xi=(long double)x[row];
  if(!isfinite(bj) || !isfinite(xi)) return 0;
  long double r=bj,magnitude=fabsl(bj);
  size_t first=transpose?(size_t)p[row]:plan->rowptr[row];
  size_t last=transpose?(size_t)p[row+1]:plan->rowptr[row+1];
  for(size_t k=first;k<last;++k) {
    const size_t at=transpose?k:(size_t)plan->positions[k];
    const size_t col=transpose?(size_t)rows[k]:(size_t)plan->columns[k];
    volatile long double av=(long double)a[at],xv=(long double)x[col];
    if(!isfinite(av)) return 0;
    const long double product=av*xv;
    r-=product;magnitude+=fabsl(product);
  }
  if(!isfinite(r) || !isfinite(magnitude)) return 0;
  if(magnitude==0) {*upper=0;return r==0;}
  const long double c=4*((long double)(last-first)+1)*LDBL_EPSILON;
  if(!(c<=0.25L)) return 0;
  const long double error=interval_outward(c*magnitude,INFINITY);
  const long double numerator=interval_outward(fabsl(r)+error,INFINITY);
  const long double denominator=interval_outward(magnitude-error,-INFINITY);
  if(!(denominator>0) || !isfinite(numerator) ||
     numerator>interval_outward(tolerance*denominator,-INFINITY)) return 0;
  *upper=interval_outward(numerator/denominator,INFINITY);
  return 1;
}

/* Taking the maximum commutes with exact binary64-to-extended conversion.
 * Keep binary64 bounds in their native format through a complete row range;
 * retain the extended maximum for any row requiring that proof. */
static int componentwise_plan_range_original(const kls_componentwise_plan *plan,
    size_t first,size_t last,const int64_t *p,const int64_t *rows,
    const double *a,const double *b,const double *x,int transpose,
    int binary64_permitted,long double tolerance,long double *upper) {
  double binary64_maximum=0;
  long double extended_maximum=0;
  for(size_t row=first;row<last;++row) {
    double binary64_bound;
    if(binary64_permitted && componentwise_binary64_row_value(
        plan,row,p,rows,a,b,x,transpose,tolerance,&binary64_bound)) {
      if(binary64_bound>binary64_maximum) binary64_maximum=binary64_bound;
    } else {
      long double extended_bound;
      if(!componentwise_plan_row(plan,row,p,rows,a,b,x,transpose,0,tolerance,&extended_bound)) return 0;
      if(extended_bound>extended_maximum) extended_maximum=extended_bound;
    }
  }
  const long double binary64_extended=(long double)binary64_maximum;
  *upper=extended_maximum>binary64_extended?extended_maximum:binary64_extended;
  return 1;
}

typedef struct { double value; int ready; } componentwise_tolerance64_cache;
#include "kls_componentwise_verdict_witness.inc"

static int componentwise_plan_range(const kls_componentwise_plan *plan,
    size_t first,size_t last,const int64_t *p,const int64_t *rows,
    const double *a,const double *b,const double *x,int transpose,
    int binary64_permitted,long double tolerance,long double *upper) {
  /* Degree+1 must fit the exact signed conversion in the witness. The plan
   * has already validated every coordinate and rejected duplicates. */
  if(plan->verdict_only && (uintmax_t)plan->nnz<=UINT64_C(0x0000ffffffffffff))
    return componentwise_plan_range_verdict_exact_domain(plan,first,last,p,rows,
        a,b,x,transpose,binary64_permitted,tolerance,upper);
  return componentwise_plan_range_original(plan,first,last,p,rows,a,b,x,
      transpose,binary64_permitted,tolerance,upper);
}

typedef struct {
  long double upper;
  int valid,environment;
} componentwise_worker_result;
typedef struct {
  const kls_componentwise_plan *plan;
  const int64_t *p,*rows;
  const double *a,*b,*x;
  int transpose,threads;
  long double tolerance;
  componentwise_worker_result *results;
} componentwise_worker_job;

static void componentwise_plan_worker(void *opaque,int tid) {
  componentwise_worker_job *job=opaque;
  if(tid<0 || tid>=job->threads) return;
  componentwise_worker_result *out=job->results+tid;
  out->environment=componentwise_plan_environment();
  if(!out->environment) return;
  /* Scope gradual underflow to certificate arithmetic, once per worker.
   * Restore every MXCSR bit, including the caller's exception flags. */
  const unsigned saved_csr=_mm_getcsr();
  const int binary64_permitted=!(saved_csr&(3u<<13));
  _mm_setcsr(saved_csr&~((1u<<15)|(1u<<6)));
  out->valid=1;
  /* Overflow-safe contiguous partition; each row keeps its original order. */
  size_t q=job->plan->n/(size_t)job->threads,r=job->plan->n%(size_t)job->threads;
  size_t first=q*(size_t)tid+(r<(size_t)tid?r:(size_t)tid);
  size_t last=first+q+((size_t)tid<r);
  out->valid=componentwise_plan_range(job->plan,first,last,job->p,job->rows,
      job->a,job->b,job->x,job->transpose,binary64_permitted,job->tolerance,&out->upper);
  _mm_setcsr(saved_csr);
}

static int componentwise_plan_run(const kls_componentwise_plan *plan,int threads,
    const int64_t *p,const int64_t *rows,const double *a,const double *b,
    const double *x,int transpose,long double *upper) {
  int valid=1,environment=1;
  long double maximum=0;
  *upper=INFINITY;
  const long double tolerance=interval_outward(
      1.0L/KLS_COMPONENTWISE_TOLERANCE_DENOMINATOR,-INFINITY);
  if(threads>1 && plan->executor) {
    if((size_t)threads>SIZE_MAX/sizeof(componentwise_worker_result)) return -2;
    componentwise_worker_result *results=calloc((size_t)threads,sizeof(*results));
    if(!results) return -2;
    componentwise_worker_job job={plan,p,rows,a,b,x,transpose,threads,tolerance,results};
    const int completed=plan->executor(plan->executor_context,componentwise_plan_worker,&job,threads);
    for(int tid=0;tid<threads;++tid) {
      environment &= results[tid].environment;
      valid &= results[tid].valid;
      if(results[tid].upper>maximum) maximum=results[tid].upper;
    }
    free(results);
    if(!completed || !environment) return -2;
    if(valid) *upper=maximum;
    return valid;
  }
  if(threads==1) {
    if(!componentwise_plan_environment()) return -2;
    const unsigned saved_csr=_mm_getcsr();
    const int binary64_permitted=!(saved_csr&(3u<<13));
    _mm_setcsr(saved_csr&~((1u<<15)|(1u<<6)));
    valid=componentwise_plan_range(plan,0,plan->n,p,rows,a,b,x,transpose,
        binary64_permitted,tolerance,&maximum);
    _mm_setcsr(saved_csr);
  } else {
#pragma omp parallel num_threads(threads) reduction(&:valid,environment) reduction(max:maximum)
    {
      const int supported=componentwise_plan_environment();
      environment&=supported;
      const unsigned saved_csr=_mm_getcsr();
      const int binary64_permitted=!(saved_csr&(3u<<13));
      if(supported) _mm_setcsr(saved_csr&~((1u<<15)|(1u<<6)));
      /* Use the actual runtime team size, including dynamic team reduction.
       * Independent contiguous ranges can stop at rejection without leaving
       * an OpenMP worksharing loop or skipping rows of an accepted proof. */
      const size_t count=(size_t)omp_get_num_threads(),tid=(size_t)omp_get_thread_num();
      const size_t q=plan->n/count,r=plan->n%count;
      const size_t first=q*tid+(r<tid?r:tid),last=first+q+(tid<r);
      if(!supported) valid=0;
      else valid&=componentwise_plan_range(plan,first,last,p,rows,a,b,x,transpose,
          binary64_permitted,tolerance,&maximum);
      _mm_setcsr(saved_csr);
    }
  }
  if(!environment) return -2;
  if(valid) *upper=maximum;
  return valid;
}

static int componentwise_plan_select(kls_componentwise_plan *plan,int threads,
    const int64_t *p,const int64_t *rows,const double *a,const double *b,
    const double *x,int transpose,long double *upper) {
  if(plan->threads!=threads) {
    plan->threads=threads;
    for(int d=0;d<2;++d) {plan->direction[d].samples=0;plan->direction[d].choice=0;}
  }
  const int choice=plan->direction[transpose].choice;
  if(choice) return componentwise_plan_run(plan,choice<0?1:threads,p,rows,a,b,x,transpose,upper);
  const int sample=plan->direction[transpose].samples;
  int result[2];long double bounds[2];double seconds[2];
  /* Clock conversion and timing ratios use SSE arithmetic too. Keep their
   * exception flags out of the caller state, including an early rejection. */
  const unsigned timing_csr=_mm_getcsr();
  for(int trial=0;trial<2;++trial) {
    const int side=(sample+trial)%2;
    const double start=omp_get_wtime();
    result[side]=componentwise_plan_run(plan,side?threads:1,p,rows,a,b,x,transpose,bounds+side);
    seconds[side]=omp_get_wtime()-start;
    /* A rejected first trial already requires the full certificate. Avoid a
     * second incomplete proof and calibrate only on accepted full traversals. */
    if(trial==0 && result[side]==0) {_mm_setcsr(timing_csr);return 0;}
  }
  if(result[0]<0 || result[1]<0 || result[0]!=result[1] ||
     (result[0]==1 && bounds[0]!=bounds[1])) {
    /* Never let a timing trial override disagreement or an unsupported worker. */
    plan->eligible=-1;_mm_setcsr(timing_csr);return 0;
  }
  if(!(seconds[0]>0) || !(seconds[1]>0) || !isfinite(seconds[0]) || !isfinite(seconds[1]))
    plan->direction[transpose].choice=-1;
  else {
    plan->direction[transpose].ratios[sample]=seconds[1]/seconds[0];
    if(++plan->direction[transpose].samples==KLS_CERTIFICATE_PROBE_PAIRS) {
      double sorted[KLS_CERTIFICATE_PROBE_PAIRS];
      for(int i=0;i<KLS_CERTIFICATE_PROBE_PAIRS;++i) {
        sorted[i]=plan->direction[transpose].ratios[i];
        for(int j=i;j>0 && sorted[j]<sorted[j-1];--j) {
          double v=sorted[j];sorted[j]=sorted[j-1];sorted[j-1]=v;
        }
      }
      plan->direction[transpose].choice=sorted[KLS_CERTIFICATE_PROBE_PAIRS/2]<1?1:-1;
    }
  }
  if(result[0]==1) *upper=bounds[0];
  _mm_setcsr(timing_csr);
  return result[0];
}
#include "kls_componentwise_interval_rows.inc"
#endif

static int componentwise_certify_executor_body(kls_componentwise_plan **owner,int threads,
    size_t n,const int64_t *p,const int64_t *rows,const double *a,
    const double *b,const double *x,int transpose,long double *upper,
    kls_componentwise_executor executor,void *context,double **correction,int verdict_only) {
  if(upper) *upper=INFINITY;
  int accumulated_rejected=0;
#if KLS_COMPONENTWISE_HAS_ACCUMULATED && defined(_OPENMP)
  if(owner && threads>=1 && n && p && b && x && upper &&
     (transpose==0 || transpose==1) && n<=SIZE_MAX/(5*sizeof(long double)) &&
     n<=SIZE_MAX/(2*sizeof(size_t)) && p[0]==0 && p[n]>=0 &&
     (!p[n] || (a && rows)) && !omp_in_parallel() && componentwise_plan_environment()) {
    if(!*owner) *owner=calloc(1,sizeof(**owner));
    if(*owner) {
      kls_componentwise_plan *plan=*owner;
      if(plan->eligible && (plan->n!=n || plan->nnz!=p[n])) kls_componentwise_plan_reset(plan);
      if(!plan->eligible) componentwise_plan_build(plan,n,p,rows);
      if(plan->eligible>0) {
        if(plan->verdict_only!=verdict_only) {
          plan->verdict_only=verdict_only;
          plan->threads=0; /* Re-measure after changing the output contract. */
        }
        if(plan->executor!=executor || plan->executor_context!=context) {
          plan->executor=executor;plan->executor_context=context;
          plan->threads=0; /* Reset both timing choices on executor changes. */
        }
        /* An existing numerical worker team owns the CPU budget. In that
         * case the caller requests serial rows without starting OpenMP trials.
         * Otherwise even discarded trials can leave competing workers awake. */
        int accepted=threads==1
          ?componentwise_plan_run(plan,1,p,rows,a,b,x,transpose,upper)
          :componentwise_plan_select(plan,threads,p,rows,a,b,x,transpose,upper);
        if(accepted==1) return 1;
        if(accepted<0) plan->eligible=-1;
        /* A supported row rejection has already run the extended accumulated
         * proof for that row. Stable scatter preserves its CSC entry order,
         * degree and separately rounded arithmetic; repeating the whole CSC
         * accumulated proof would reject the same row. Go directly to the
         * full interval certificate, which still owns correction capture.
         * A declined executor/environment or disagreeing trials never supply
         * this fact and retain the ordinary stateless path. */
        accumulated_rejected=accepted==0 && plan->eligible>0;
        if(accumulated_rejected && correction) {
          const int interval=componentwise_interval_rows(plan,p,rows,a,b,x,
              transpose,upper,correction);
          if(interval!=-2) return interval;
        }
      }
    }
  }
#else
  (void)owner;(void)threads;(void)executor;(void)context;(void)verdict_only;
#endif
  return componentwise_certify_body(n,p,rows,a,b,x,transpose,NULL,upper,correction,
                                   accumulated_rejected);
}

int kls_componentwise_certify_executor(kls_componentwise_plan **owner,int threads,
    size_t n,const int64_t *p,const int64_t *rows,const double *a,
    const double *b,const double *x,int transpose,long double *upper,
    kls_componentwise_executor executor,void *context) {
  return componentwise_certify_executor_body(owner,threads,n,p,rows,a,b,x,
      transpose,upper,executor,context,NULL,0);
}

int kls_componentwise_certify_correction(kls_componentwise_plan **owner,int threads,
    size_t n,const int64_t *p,const int64_t *rows,const double *a,
    const double *b,const double *x,int transpose,long double *upper,
    kls_componentwise_executor executor,void *context,double **correction) {
  return componentwise_certify_executor_body(owner,threads,n,p,rows,a,b,x,
      transpose,upper,executor,context,correction,0);
}

/* On acceptance only the verdict is consumed. On rejection the existing
 * interval path still supplies the bound and correction RHS, unchanged. */
int kls_componentwise_certify_verdict_correction(kls_componentwise_plan **owner,int threads,
    size_t n,const int64_t *p,const int64_t *rows,const double *a,
    const double *b,const double *x,int transpose,long double *failure_upper,
    kls_componentwise_executor executor,void *context,double **correction) {
  return componentwise_certify_executor_body(owner,threads,n,p,rows,a,b,x,
      transpose,failure_upper,executor,context,correction,1);
}

int kls_componentwise_certify_planned(kls_componentwise_plan **owner,int threads,
    size_t n,const int64_t *p,const int64_t *rows,const double *a,
    const double *b,const double *x,int transpose,long double *upper) {
  return kls_componentwise_certify_executor(owner,threads,n,p,rows,a,b,x,
                                          transpose,upper,NULL,NULL);
}
