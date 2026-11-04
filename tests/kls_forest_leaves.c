/* Current-input differential coverage; no policy switches in the installed API. */
#define _GNU_SOURCE
#define _POSIX_C_SOURCE 200809L
#include <stdlib.h>
#include <stdatomic.h>
static atomic_int fail_next,fail_next_calloc;
static unsigned tail_mask;
static void *leaf_malloc(size_t);
static void *leaf_calloc(size_t,size_t);
#undef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#define malloc leaf_malloc
#define calloc leaf_calloc
static _Thread_local unsigned long long remaining_plan_test_calls;
#define KLS_FOREST_REMAINING_TEST_EVENT() (++remaining_plan_test_calls)
static _Thread_local unsigned long long staged_leaf_test_calls;
#define KLS_FOREST_STAGE_TEST_EVENT() (++staged_leaf_test_calls)
#include "../src/kls.c"
#undef malloc
#undef calloc
#undef NDEBUG
#include <assert.h>
static void *leaf_malloc(size_t n) {
  return atomic_exchange_explicit(&fail_next,0,memory_order_relaxed)?NULL:malloc(n);
}
static void *leaf_calloc(size_t n,size_t bytes) {
  return atomic_exchange_explicit(&fail_next_calloc,0,memory_order_relaxed)?NULL:calloc(n,bytes);
}
#if defined(__x86_64__) && FLT_EVAL_METHOD==0
static void same_numeric(kls_solver *a,kls_solver *b) {
  const size_t bytes=(size_t)a->n*sizeof(double);
  assert(!memcmp(a->numeric->Udiag,b->numeric->Udiag,bytes));
  if(a->numeric->Rs) assert(!memcmp(a->numeric->Rs,b->numeric->Rs,bytes));
  for(UF_long k=0;k<a->n;++k) for(int u=0;u<2;++u) {
    UF_long *ai,*bi,an,bn;double *av,*bv;
    kls_klu_get_pointer(a->numeric->LUbx[0],u?a->numeric->Uip:a->numeric->Lip,
      u?a->numeric->Ulen:a->numeric->Llen,k,&ai,&av,&an);
    kls_klu_get_pointer(b->numeric->LUbx[0],u?b->numeric->Uip:b->numeric->Lip,
      u?b->numeric->Ulen:b->numeric->Llen,k,&bi,&bv,&bn);
    assert(an==bn);
    if(an) {assert(!memcmp(ai,bi,(size_t)an*sizeof(*ai)));assert(!memcmp(av,bv,(size_t)an*sizeof(*av)));}
  }
}
static kls_solver *create(int n,int64_t *p,int64_t *ix,double *a,
    int csr,int transpose,int scale) {
  kls_options o;kls_default_options(&o);o.backend=KLS_BACKEND_SERIAL;
  o.threads=1;o.use_btf=0;o.ordering=KLS_ORDERING_NATURAL;o.static_pivoting=0;
  o.orientation=transpose?KLS_ORIENTATION_TRANSPOSE:KLS_ORIENTATION_NORMAL;o.scale=scale;
  kls_solver *s=NULL;assert(kls_create(&s)==KLS_OK);
  assert((csr?kls_analyze_csr:kls_analyze_csc)(s,KLS_INDEX_INT64,n,p,ix,0,&o)==KLS_OK);
  assert(kls_factor(s,a)==KLS_OK);
  s->compact_map32_probe_active=1;
  assert(kls_build_refactor_map(s)&&kls_i32_solve_ready(s));
  s->compact_map32_probe_active=0;
  free(s->fused_refactor_solve_work);
  s->fused_refactor_solve_work=calloc((size_t)2*n,sizeof(double));
  assert(s->fused_refactor_solve_work);s->fused_refactor_solve_work_n=n;
  return s;
}
static void answer(int n,const double *x,const double *truth) {
  for(int i=0;i<n;++i) assert(isfinite(x[i])&&fabs(x[i]-truth[i])<1e-10);
}
static void rhs(int n,const int64_t *p,const int64_t *ix,const double *a,
    int csr,int transpose,const double *truth,double *b) {
  memset(b,0,(size_t)n*sizeof(*b));
  for(int j=0;j<n;++j) for(int64_t z=p[j];z<p[j+1];++z) {
    int row=csr?j:(int)ix[z],col=csr?(int)ix[z]:j;
    if(transpose) b[col]+=a[z]*truth[row];else b[row]+=a[z]*truth[col];
  }
}
static void shadow_checks(kls_solver *s,double *values,const double *b) {
  const kls_forest_leaf_plan *p=s->forest_leaves;const UF_long count=p->count,n=s->n;
  const int32_t *offset=p->indices+5*count;
  const double *scale=s->numeric->Rs?s->fused_refactor_solve_work+n:NULL;
  double *diagonal=malloc((size_t)n*sizeof(*diagonal));
  double *forward=malloc((size_t)n*sizeof(*forward));
  double *lower=malloc((size_t)count*sizeof(*lower));assert(diagonal&&forward&&lower);
  memcpy(diagonal,s->numeric->Udiag,(size_t)n*sizeof(*diagonal));
  memcpy(forward,s->fused_refactor_solve_work,(size_t)n*sizeof(*forward));
  for(UF_long i=0;i<count;++i) lower[i]=((double *)s->numeric->LUbx[0])[offset[i]];
  fenv_t environment;assert(fegetenv(&environment)==0);const unsigned csr=_mm_getcsr();
  feclearexcept(FE_ALL_EXCEPT);
  assert(kls_forest_leaf_validate_scalar(p,values,b,scale,0));
  const int scalar_flags=fetestexcept(FE_ALL_EXCEPT);
  kls_forest_leaf_publish_scalar(s,p,values,b,scale,s->fused_refactor_solve_work,s->i32solve_udiag_recip,0);
#if KLS_HAVE_AVX512_KERNELS
  if(kls_avx512_scatter_enabled()) {
    feclearexcept(FE_ALL_EXCEPT);
    assert(kls_forest_leaf_validate_avx512(p,values,b,scale));
    assert(scalar_flags==fetestexcept(FE_ALL_EXCEPT));
    const UF_long first=kls_forest_leaf_publish_avx512(s,p,values,b,scale,
      s->fused_refactor_solve_work,s->i32solve_udiag_recip);
    kls_forest_leaf_publish_scalar(s,p,values,b,scale,s->fused_refactor_solve_work,s->i32solve_udiag_recip,first);
  }
#else
  (void)scalar_flags;
#endif
  assert(!memcmp(diagonal,s->numeric->Udiag,(size_t)n*sizeof(*diagonal)));
  assert(!memcmp(forward,s->fused_refactor_solve_work,(size_t)n*sizeof(*forward)));
  for(UF_long i=0;i<count;++i) assert(!memcmp(lower+i,(double *)s->numeric->LUbx[0]+offset[i],sizeof(*lower)));
  const size_t saved_worksize=s->numeric->worksize;
  s->numeric->worksize=(size_t)n*sizeof(double);
  const unsigned long long staged_before=staged_leaf_test_calls;
  assert(kls_forest_leaf_begin(s,values,b,scale,s->fused_refactor_solve_work,s->i32solve_udiag_recip));
  assert(staged_leaf_test_calls==staged_before);
  s->numeric->worksize=saved_worksize;
  assert(!memcmp(diagonal,s->numeric->Udiag,(size_t)n*sizeof(*diagonal)));
  assert(!memcmp(forward,s->fused_refactor_solve_work,(size_t)n*sizeof(*forward)));
  for(UF_long i=0;i<count;++i)assert(!memcmp(lower+i,(double *)s->numeric->LUbx[0]+offset[i],sizeof(*lower)));
  const int32_t source=p->indices[3*count-1];const double old=values[source];
  const double bad[]={0,DBL_TRUE_MIN,DBL_MAX,INFINITY,NAN};
  for(size_t t=0;t<sizeof(bad)/sizeof(*bad);++t) {
    values[source]=bad[t];const unsigned before=_mm_getcsr();
    assert(!kls_forest_leaf_begin(s,values,b,scale,s->fused_refactor_solve_work,s->i32solve_udiag_recip));
    const unsigned after=_mm_getcsr();
    assert(after==before);
    assert(!memcmp(diagonal,s->numeric->Udiag,(size_t)n*sizeof(*diagonal)));
    assert(!memcmp(forward,s->fused_refactor_solve_work,(size_t)n*sizeof(*forward)));
    for(UF_long i=0;i<count;++i) assert(!memcmp(lower+i,(double *)s->numeric->LUbx[0]+offset[i],sizeof(*lower)));
  }
  values[source]=old;
  _mm_setcsr((csr&~0x3fu)&~(1u<<7));const unsigned unmasked=_mm_getcsr();
  assert(!kls_forest_leaf_begin(s,values,b,scale,s->fused_refactor_solve_work,NULL));
  assert(_mm_getcsr()==unmasked);
  _mm_setcsr((csr&~0x6000u)|0x2000u);
  const unsigned mismatched=_mm_getcsr();
  assert(fegetround()==FE_TONEAREST);
  assert(!kls_forest_leaf_begin(s,values,b,scale,s->fused_refactor_solve_work,NULL));
  assert(_mm_getcsr()==mismatched);
  assert(fesetenv(&environment)==0);_mm_setcsr(csr);
  free(diagonal);free(forward);free(lower);
}
#define KLS_FOREST_RESUME_TEST 1
#include "forest_unscaled_checks.inc"

static void exercise(int groups,int csr,int transpose,int scale) {
  const int n=4*groups+2,capacity=3*n-2;
  int *parent=malloc((size_t)n*sizeof(*parent));
  int64_t *p=malloc((size_t)(n+1)*sizeof(*p)),*ix=malloc((size_t)capacity*sizeof(*ix));
  double *a=malloc((size_t)capacity*sizeof(*a)),*initial=malloc((size_t)capacity*sizeof(*initial));
  double *b=malloc((size_t)n*sizeof(*b)),*truth=malloc((size_t)n*sizeof(*truth));
  double *x=malloc((size_t)n*sizeof(*x)),*y=malloc((size_t)n*sizeof(*y));
  assert(parent&&p&&ix&&a&&initial&&b&&truth&&x&&y);
  for(int g=0;g<groups;++g) {
    for(int j=0;j<3;++j) parent[4*g+j]=4*g+3;
    parent[4*g+3]=g+1<groups?4*g+7:n-1;
  }
  parent[n-2]=n-1;parent[n-1]=-1;
  const int isolated=groups==18&&csr&&transpose&&scale==2;
  if(isolated) parent[n-2]=-1;
  int nz=0;
  for(int j=0;j<n;++j) {
    p[j]=nz;truth[j]=1+.03125*(j%7-3);
    for(int i=0;i<n;++i) if(i==j||parent[i]==j||parent[j]==i) {
      ix[nz]=i;int row=csr?j:i,col=csr?i:j;
      a[nz]=initial[nz]=row==col?16+.125*(row%7):-.03125*(1+(row+2*col)%3);++nz;
    }
  }
  p[n]=nz;assert(nz==capacity-2*isolated);
  if(groups==17&&csr&&transpose&&scale==2)forest_unscaled_public_checks(n,p,ix,a,nz);
  kls_solver *r=create(n,p,ix,a,csr,transpose,scale),*f=create(n,p,ix,a,csr,transpose,scale);
  r->forest_map_state=-1;f->options.backend=KLS_BACKEND_AUTO;
  f->private_caller_depth=1;f->solve_contract_probe=2;
  for(int epoch=0;epoch<4;++epoch) {
    for(int z=0;z<nz;++z) a[z]=initial[z]*(1+.0009765625*(epoch+1)*(1+z%3));
    for(int j=0;j<n;++j) truth[j]+=.015625;
    rhs(n,p,ix,a,csr,r->orientation==KLS_ORIENTATION_NORMAL,truth,b);
    double *rv=NULL,*fv=NULL;
    assert(prepare_numeric_values(r,a,&rv)==KLS_OK&&prepare_numeric_values(f,a,&fv)==KLS_OK);
    r->fused_refactor_solve_rhs=f->fused_refactor_solve_rhs=b;
    r->transpose_cycle.active=f->transpose_cycle.active=1;
    feclearexcept(FE_ALL_EXCEPT);assert(kls_lean_single_block_map32_refactor(r,rv,1)==1);
    int flags=fetestexcept(FE_ALL_EXCEPT);
    feclearexcept(FE_ALL_EXCEPT);assert(kls_lean_single_block_map32_refactor(f,fv,1)==1);
    assert(flags==fetestexcept(FE_ALL_EXCEPT));same_numeric(r,f);
    assert(f->forest_leaf_state==1&&f->forest_leaf_ready);
    tail_mask|=1u<<(f->forest_leaves->count%8);
    assert(f->forest_leaf_numerics==(uint64_t)epoch+1);
    assert(!memcmp(r->fused_refactor_solve_work,f->fused_refactor_solve_work,
      (size_t)n*(scale>0?2:1)*sizeof(double)));
    if(epoch==0) shadow_checks(f,fv,b);
    assert(kls_compact_transpose_solve(r,b,x)&&kls_compact_transpose_solve(f,b,y));
    assert(f->forest_leaf_solves==(uint64_t)epoch+1&&!f->forest_leaf_ready);
    answer(n,x,truth);answer(n,y,truth);
    /* A later RHS must execute the ordinary complete triangular solve. */
    for(int j=0;j<n;++j) b[j]=-b[j];
    assert(kls_compact_transpose_solve(f,b,y));
    for(int j=0;j<n;++j) assert(fabs(y[j]+truth[j])<1e-10);
    assert(f->forest_leaf_solves==(uint64_t)epoch+1);
  }
  /* Optional plan failure is sticky and changes no numerical value. */
  kls_free_forest_map(f);assert(kls_prepare_forest_map(f));
  memcpy(x,f->numeric->Udiag,(size_t)n*sizeof(*x));
  atomic_store_explicit(&fail_next,1,memory_order_relaxed);
  assert(!kls_prepare_forest_leaves(f)&&f->forest_leaf_state<0);
  assert(!atomic_load_explicit(&fail_next,memory_order_relaxed));
  assert(!memcmp(x,f->numeric->Udiag,(size_t)n*sizeof(*x)));
  kls_free_forest_map(f);assert(kls_prepare_forest_map(f));
  atomic_store_explicit(&fail_next_calloc,1,memory_order_relaxed);
  assert(!kls_prepare_forest_leaves(f)&&f->forest_leaf_state<0);
  assert(!atomic_load_explicit(&fail_next_calloc,memory_order_relaxed));
  assert(!memcmp(x,f->numeric->Udiag,(size_t)n*sizeof(*x)));
  kls_free_forest_map(f);assert(kls_prepare_forest_map(f)&&kls_prepare_forest_leaves(f));
  /* Unsupported rounding uses the unmodified natural-order incumbent. */
  assert(fesetround(FE_UPWARD)==0);
  double *fv=NULL,*rv=NULL;assert(prepare_numeric_values(f,a,&fv)==KLS_OK);
  assert(prepare_numeric_values(r,a,&rv)==KLS_OK);
  f->fused_refactor_solve_rhs=r->fused_refactor_solve_rhs=b;
  assert(kls_lean_single_block_map32_refactor(r,rv,1)==1);
  uint64_t before=f->forest_leaf_numerics;
  assert(kls_lean_single_block_map32_refactor(f,fv,1)==1);
  assert(f->forest_leaf_numerics==before&&!f->forest_leaf_ready);same_numeric(r,f);
  assert(fesetround(FE_TONEAREST)==0);
  if(transpose&&csr) {
    /* A genuinely public combined call still captures CURRENT A and RHS,
       runs the full caller certificate, and publishes only successful output.
       Settle the already validated compact owner independently of timing. */
    f->private_caller_depth=0;f->transpose_cycle.active=0;
    f->factor_preps_deferred=f->prestatic_deferred=f->rowmatch_deferred=0;
    f->metis_race_deferred=f->auto_scale_deferred=f->tight_pivot_deferred=0;
    f->lean_choice=-1;f->lean_probe_arm=0;f->unchanged_refactor_state=-1;
    f->transpose_cycle.choice=1;f->transpose_cycle.cycles=32;
    /* The serial constructor made exactly the native factors validated
       above. Construct the historical tag needed by this private fixture;
       this is not evidence of public policy or real-circuit activation. */
    kls_set_last_factor_path(f,KLS_FACTOR_PATH_KLU_FIRST);
    kls_set_last_refactor_path(f,KLS_REFACTOR_PATH_KLU);
    for(int crew=1;crew<=8;crew+=7) {
      f->options.threads=crew;
      assert(kls_compact_transpose_cycle_capable(f));
      for(int epoch=0;epoch<3;++epoch) {
        for(int z=0;z<nz;++z) a[z]=initial[z]*(1+.0009765625*(epoch+crew+7)*(1+z%3));
        for(int j=0;j<n;++j) truth[j]+=.015625;
        rhs(n,p,ix,a,csr,0,truth,b);memcpy(x,b,(size_t)n*sizeof(*x));
        const uint64_t numeric=f->forest_leaf_numerics,solves=f->forest_leaf_solves;
        double *destination=epoch==0&&crew==1?a:x;
        assert(kls_refactor_solve(f,a,1,epoch&1?x:b,0,destination,0)==KLS_OK);
        assert(f->forest_leaf_numerics==numeric+1&&f->forest_leaf_solves==solves+1);
        assert(!f->forest_leaf_ready&&f->private_original_valid);answer(n,destination,truth);
      }
    }
    double *multi=malloc((size_t)2*n*sizeof(*multi)),*mx=malloc((size_t)2*n*sizeof(*mx));
    assert(multi&&mx);
    for(int j=0;j<n;++j) {multi[j]=b[j];multi[n+j]=2*b[j];}
    const uint64_t before_multi=f->forest_leaf_numerics;
    assert(kls_refactor_solve(f,a,2,multi,0,mx,0)==KLS_OK);
    assert(f->forest_leaf_numerics==before_multi&&!f->forest_leaf_ready);
    answer(n,mx,truth);
    for(int j=0;j<n;++j) assert(fabs(mx[n+j]-2*truth[j])<1e-10);
    rhs(n,p,ix,a,csr,1,truth,b);
    assert(kls_solve_transpose(f,1,b,0,x,0)==KLS_OK);answer(n,x,truth);
    free(multi);free(mx);f->private_caller_depth=1;
  }
  /* Restore matching complete numerics before checking partial singular
     publication. The public lifecycle above deliberately changed A again. */
  assert(prepare_numeric_values(r,a,&rv)==KLS_OK&&prepare_numeric_values(f,a,&fv)==KLS_OK);
  r->fused_refactor_solve_rhs=f->fused_refactor_solve_rhs=b;
  assert(kls_lean_single_block_map32_refactor(r,rv,1)==1);
  assert(kls_lean_single_block_map32_refactor(f,fv,1)==1);same_numeric(r,f);
  if(groups==17&&csr&&transpose&&scale==2) {
    /* A singular spine is a partial numeric failure, not optional decline.
       Leaf inputs remain unchanged, so even the incumbent's partial values
       agree in this fixture. No prepared result may survive the failure. */
    const UF_long k=(UF_long)f->forest_leaves->indices[6*f->forest_leaves->count];
    for(UF_long j=r->refactor_col_ptr[k];j<r->refactor_col_ptr[k+1];++j)
      rv[r->refactor_input_pos32[j]]=0;
    for(UF_long j=f->refactor_col_ptr[k];j<f->refactor_col_ptr[k+1];++j)
      fv[f->refactor_input_pos32[j]]=0;
    const uint64_t previous=f->forest_leaf_numerics;
    assert(kls_lean_single_block_map32_refactor(r,rv,1)==0);
    assert(kls_lean_single_block_map32_refactor(f,fv,1)==0);
    assert(!f->forest_leaf_ready&&f->forest_leaf_numerics==previous);
    assert(r->common.status==f->common.status&&r->common.numerical_rank==f->common.numerical_rank);
    same_numeric(r,f);
  }
  /* A zero leaf pivot declines before publication; the singular incumbent's
     rank, factor bytes and saved forward values remain unchanged. */
  memset(a,0,(size_t)nz*sizeof(*a));
  assert(prepare_numeric_values(r,a,&rv)==KLS_OK&&prepare_numeric_values(f,a,&fv)==KLS_OK);
  feclearexcept(FE_ALL_EXCEPT);
  const int ro=kls_lean_single_block_map32_refactor(r,rv,1);
  const int singular_flags=fetestexcept(FE_ALL_EXCEPT);
  feclearexcept(FE_ALL_EXCEPT);
  const int fo=kls_lean_single_block_map32_refactor(f,fv,1);
  assert(singular_flags==fetestexcept(FE_ALL_EXCEPT));
  assert(ro==fo&&r->common.status==f->common.status);
  assert(r->common.numerical_rank==f->common.numerical_rank&&!f->forest_leaf_ready);
  same_numeric(r,f);
  kls_invalidate_i32_solve(f);assert(!f->forest_leaves&&!f->forest_leaf_state&&!f->forest_leaf_ready);
  f->private_caller_depth=0;
  for(int j=0;j<n;++j) x[j]=-123;
  assert(kls_refactor_solve(f,a,1,b,0,x,0)!=KLS_OK);
  for(int j=0;j<n;++j) assert(x[j]==-123);
  assert(!f->forest_leaf_ready);
  kls_destroy(r);kls_destroy(f);free(parent);free(p);free(ix);free(a);free(initial);
  free(b);free(truth);free(x);free(y);
}
#if KLS_HAVE_AVX512_KERNELS
static void reciprocal_tails(void) {
  const double values[]={1,2,3,0x1p-1022,0x1p1023,DBL_TRUE_MIN,
    DBL_MAX,INFINITY,NAN,0,-0.0,-3};
  double scale[71],reference[71],actual[71];
  unsigned original;__asm__ volatile("stmxcsr %0":"=m"(original)::"memory");
  for(int mode=0;mode<4;++mode) for(int shift=0;shift<12;++shift)
    for(int n=64;n<=71;++n) {
      const unsigned control=(original&~0x3fu)|(mode&1?1u<<15:0)|(mode&2?1u<<6:0);
      for(int k=0;k<n;++k) scale[k]=values[(k+shift)%12];
      __asm__ volatile("ldmxcsr %0"::"m"(control):"memory");
      for(int k=0;k<n;++k) {volatile double operand=scale[k];reference[k]=1.0/operand;}
      unsigned expected;__asm__ volatile("stmxcsr %0":"=m"(expected)::"memory");
      __asm__ volatile("ldmxcsr %0"::"m"(control):"memory");
      kls_forest_scale_reciprocals(n,scale,actual);
      unsigned observed;__asm__ volatile("stmxcsr %0":"=m"(observed)::"memory");
      assert(expected==observed);
      assert(!memcmp(reference,actual,(size_t)n*sizeof(double)));
    }
  __asm__ volatile("ldmxcsr %0"::"m"(original):"memory");
}
#endif
int main(void) {
#if !KLS_HAVE_AVX512_KERNELS
  return 0;
#else
  if(!kls_avx512_scatter_enabled()) return 0;
#endif
#if KLS_HAVE_AVX512_KERNELS
  reciprocal_tails();
#endif
  assert(!setenv("KLS_EXPERIMENT_FOREST_STAGE_LEAVES","1",1));
  assert(!setenv("KLS_EXPERIMENT_FOREST_REMAINING_PLAN","1",1));
  atomic_init(&fail_next,0);atomic_init(&fail_next_calloc,0);
  for(int groups=17;groups<=24;++groups) for(int csr=0;csr<2;++csr)
    for(int transpose=0;transpose<2;++transpose) for(int scale=0;scale<3;++scale)
      exercise(groups,csr,transpose,scale?scale:-1);
  assert(remaining_plan_test_calls>0);
  fprintf(stderr,"FOREST_REMAINING_PLAN_CALLS %llu\n",remaining_plan_test_calls);
  assert(staged_leaf_test_calls>0);
  fprintf(stderr,"FOREST_STAGED_LEAF_CALLS %llu\n",staged_leaf_test_calls);
  assert(tail_mask==255);
  const unsigned csr=_mm_getcsr();
  _mm_setcsr(csr|(1u<<15)|(1u<<6));
  exercise(17,1,1,2);
  _mm_setcsr(csr);
  return 0;
}

#else
int main(void) {return 0;}
#endif
