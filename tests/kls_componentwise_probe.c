/* Test-only access: no diagnostic entry points in the production library. */
#include <stdio.h>
#include "../src/kls_componentwise.c"

int kls_test_fast_available(void) { return KLS_COMPONENTWISE_HAS_ACCUMULATED; }

int kls_test_binary64_classification(void) {
#if KLS_COMPONENTWISE_HAS_ACCUMULATED && defined(_OPENMP)
  const unsigned saved=_mm_getcsr();
  _mm_setcsr(saved&~((1u<<15)|(1u<<6)));
  const uint64_t fractions[]={0,1,UINT64_C(0x8000000000000),UINT64_C(0xfffffffffffff)};
  int valid=1;
  for(unsigned sign=0;sign<2;++sign) for(unsigned exponent=0;exponent<2048;++exponent)
    for(size_t f=0;f<sizeof(fractions)/sizeof(*fractions);++f) {
      uint64_t bits=((uint64_t)sign<<63)|((uint64_t)exponent<<52)|fractions[f];
      double value;memcpy(&value,&bits,sizeof(value));
      int expected=isfinite(value) && (value==0 || isnormal(value));
      valid&=componentwise_binary64_normal_or_zero(value)==expected;
    }
  _mm_setcsr(saved);return valid;
#else
  return -1;
#endif
}

#if KLS_COMPONENTWISE_HAS_ACCUMULATED && defined(_OPENMP)
static int test_scoped_executor(void *context,void (*worker)(void *,int),
    void *job,int threads) {
  (void)context;
  for(int tid=0;tid<threads;++tid) {
    unsigned before=_mm_getcsr();
    worker(job,tid);
    if(_mm_getcsr()!=before) return 0;
  }
  return 1;
}
#endif

int kls_test_binary64_scoped_environment(int mode,int reject,int transpose) {
#if KLS_COMPONENTWISE_HAS_ACCUMULATED && defined(_OPENMP)
  const unsigned original=_mm_getcsr();
  int64_t p[]={0,1},rows[]={0};double a[]={1},b[]={reject?2:1},x[]={1};
  kls_componentwise_plan plan={0};long double upper;
  componentwise_plan_build(&plan,1,p,rows);
  if(mode==3) plan.executor=test_scoped_executor;
  const unsigned scoped=original|(1u<<15)|(1u<<6);
  _mm_setcsr(scoped);
  int result=componentwise_plan_run(&plan,mode==1?1:2,p,rows,a,b,x,transpose,&upper);
  int restored=_mm_getcsr()==scoped;
  _mm_setcsr(original);
  kls_componentwise_plan_reset(&plan);
  /* The double proof gives a larger valid bound than the extended proof:
   * check that acceptance actually traversed the new arithmetic path. */
  return restored && result==!reject && (reject || upper>DBL_EPSILON);
#else
  (void)mode;(void)reject;(void)transpose;return -1;
#endif
}

/* Save and restore the complete environment: these changes are test-local. */
int kls_test_binary64_environment(int mode) {
#if KLS_COMPONENTWISE_HAS_ACCUMULATED && defined(_OPENMP)
  fenv_t saved;
  if(fegetenv(&saved)) return -1;
  unsigned csr=_mm_getcsr();
  if(mode==1) _mm_setcsr(csr|(1u<<15));
  if(mode==2) _mm_setcsr(csr|(1u<<6));
  if(mode==3) fesetround(FE_UPWARD);
  if(mode==4) fesetround(FE_DOWNWARD);
  if(mode==5) fesetround(FE_TOWARDZERO);
  int64_t p[]={0,1},rows[]={0};double a[]={1},b[]={1},x[]={1};
  kls_componentwise_plan plan={0};long double upper;
  componentwise_plan_build(&plan,1,p,rows);
  int result=componentwise_plan_environment() && componentwise_binary64_row(
      &plan,0,p,rows,a,b,x,0,1.0L/100000000,&upper);
  kls_componentwise_plan_reset(&plan);
  fesetenv(&saved);_mm_setcsr(csr);
  return result;
#else
  (void)mode;return -1;
#endif
}

int kls_test_binary64(size_t n,const int64_t *p,const int64_t *rows,
    const double *a,const double *b,const double *x,int transpose,
    double *residual,long double *upper) {
  (void)residual;*upper=INFINITY;
#if KLS_COMPONENTWISE_HAS_ACCUMULATED && defined(_OPENMP)
  kls_componentwise_plan plan={0};
  if(!componentwise_plan_environment()) return 0;
  componentwise_plan_build(&plan,n,p,rows);
  int result=plan.eligible>0;
  long double maximum=0,tolerance=interval_outward(
      1.0L/KLS_COMPONENTWISE_TOLERANCE_DENOMINATOR,-INFINITY);
  for(size_t i=0;result && i<n;++i) {
    long double bound;
    result=componentwise_binary64_row(&plan,i,p,rows,a,b,x,transpose,tolerance,&bound);
    if(result && bound>maximum) maximum=bound;
  }
  if(result) *upper=maximum;
  kls_componentwise_plan_reset(&plan);return result;
#else
  (void)n;(void)p;(void)rows;(void)a;(void)b;(void)x;(void)transpose;
  return 0;
#endif
}

int kls_test_fast(size_t n, const int64_t *p, const int64_t *rows,
    const double *a, const double *b, const double *x, int transpose,
    double *residual, long double *upper) {
  (void)residual;
  *upper=INFINITY;
  /* Fixtures supply valid CSC and the default arithmetic environment. */
  return componentwise_accumulated_certificate(n,p,rows,a,b,x,transpose,upper);
}

int kls_test_bound(const long double *upper, char *out, size_t capacity) {
  /* Hexadecimal text retains every significand bit, unlike ctypes.value. */
  return snprintf(out,capacity,"%.*La",(LDBL_MANT_DIG+3)/4,*upper);
}

int kls_test_parallel_available(void) {
#if KLS_COMPONENTWISE_HAS_ACCUMULATED && defined(_OPENMP)
  return 1;
#else
  return 0;
#endif
}

int kls_test_parallel(size_t n,const int64_t *p,const int64_t *rows,
    const double *a,const double *b,const double *x,int transpose,
    double *residual,long double *upper,int threads) {
  (void)residual;*upper=INFINITY;
#if KLS_COMPONENTWISE_HAS_ACCUMULATED && defined(_OPENMP)
  kls_componentwise_plan plan={0};
  if(!n) {*upper=0;return 1;}
  componentwise_plan_build(&plan,n,p,rows);
  int result=plan.eligible>0?componentwise_plan_run(&plan,threads,p,rows,a,b,x,transpose,upper):0;
  kls_componentwise_plan_reset(&plan);
  return result;
#else
  (void)n;(void)p;(void)rows;(void)a;(void)b;(void)x;(void)transpose;(void)threads;
  return 0;
#endif
}
