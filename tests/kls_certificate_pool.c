#include "../src/kls.c"
#undef NDEBUG
#include <assert.h>
#include <fenv.h>

static void count_worker(void *opaque,int tid) {
  int *hits=opaque;assert(tid>=0 && tid<8);++hits[tid];
}
#if defined(_OPENMP) && defined(__SSE__) && defined(__x86_64__) && \
    FLT_RADIX == 2 && DBL_MANT_DIG == 53 && LDBL_MANT_DIG == 64 && \
    LDBL_MAX_EXP == 16384
typedef struct {fenv_t saved[8];int restore,sse_only;} certificate_team_environment;
static void change_certificate_worker_environment(void *opaque,int tid) {
  certificate_team_environment *env=opaque;
  assert(tid>=0 && tid<8);
  if(tid==0) return; /* The public caller retains its original environment. */
  if(env->restore) assert(!fesetenv(env->saved+tid));
  else {
    assert(!fegetenv(env->saved+tid));
    if(env->sse_only) {
      /* Leave x87 at nearest while independently changing SSE rounding. */
      _mm_setcsr((_mm_getcsr()&~(3u<<13))|(2u<<13));
    } else assert(!fesetround(FE_DOWNWARD));
  }
}
static void exercise_dedicated_environment(kls_solver *s,int threads,
    const int64_t *p,const int64_t *rows,const double *a,int sse_only) {
  certificate_team_environment env={.sse_only=sse_only};
  const int caller_rounding=fegetround();
  assert(private_componentwise_execute_dedicated(s,
      change_certificate_worker_environment,&env,threads));
  kls_componentwise_plan *plan=NULL;
  double b[]={a[0],4,5,7},bt[]={a[0]+1,3,5,7},x[]={1,1,1,1};
  for(int tr=0;tr<2;++tr) {
    const unsigned long before=s->componentwise_team->generation;
    long double upper,reference;
    assert(kls_componentwise_certify(4,p,rows,a,tr?bt:b,x,tr,NULL,&reference)==1);
    /* Unsupported x87 rounding declines the worker proof and uses the
     * stateless certificate. Independent SSE rounding uses the supported
     * extended-precision row proof, whose conservative bound can differ. */
    assert(kls_componentwise_certify_executor(&plan,threads,4,p,rows,a,tr?bt:b,x,tr,
             &upper,private_componentwise_execute_dedicated,s)==1);
    if(!sse_only) assert(upper==reference);
    else assert(isfinite(upper) && upper>=0 &&
                upper<=1.0L/KLS_COMPONENTWISE_TOLERANCE_DENOMINATOR);
    if(!tr) assert(s->componentwise_team->generation!=before);
    assert(fegetround()==caller_rounding);
    assert(!s->componentwise_work && !s->componentwise_work_context);
  }
  env.restore=1;
  assert(private_componentwise_execute_dedicated(s,
      change_certificate_worker_environment,&env,threads));
  kls_componentwise_plan_destroy(plan);
}
#endif
static void exercise(int family,int threads) {
  int64_t p[]={0,2,3,4,5},rows[]={0,1,1,2,3};
  double a[]={2,1,3,5,7},b[]={2,4,5,7},bt[]={3,3,5,7},x[]={1,1,1,1};
  kls_options options;kls_default_options(&options);options.threads=1;
  kls_solver *s=NULL;assert(kls_create(&s)==KLS_OK);
  assert(kls_set_accuracy_policy(s,KLS_ACCURACY_COMPONENTWISE_BACKWARD_ERROR)==KLS_OK);
  assert(kls_analyze_csc(s,KLS_INDEX_INT64,4,p,rows,0,&options)==KLS_OK);
  assert(kls_factor(s,a)==KLS_OK);
  destroy_refactor_pool(s);destroy_egraph_refactor_pool(s);
  if(family) assert(ensure_egraph_refactor_pool(s,threads));
  else assert(ensure_refactor_pool(s,threads));
  assert(private_componentwise_pool_threads(s)==threads);
  int hits[8]={0};
  /* A stale numeric mode must never be entered by a certificate dispatch. */
  if(family) s->egraph_pool->shared.solve_permute_mode=1;
  for(int repeat=0;repeat<20;++repeat)
    assert(private_componentwise_execute(s,count_worker,hits,threads)==1);
  for(int tid=0;tid<threads;++tid) assert(hits[tid]==20);
  assert(!s->componentwise_work && !s->componentwise_work_context);
  assert(!private_componentwise_execute(s,count_worker,hits,threads+1));
  if(family) s->egraph_pool->shared.solve_permute_mode=0;
  kls_componentwise_plan *plan=NULL,*serial_plan=NULL;
  for(int tr=0;tr<2;++tr) for(int repeat=0;repeat<8;++repeat) {
    long double upper,reference;
    assert(kls_componentwise_certify_executor(&plan,threads,4,p,rows,a,tr?bt:b,x,tr,
             &upper,private_componentwise_execute,s)==1);
    assert(kls_componentwise_certify(4,p,rows,a,tr?bt:b,x,tr,NULL,&reference)==1);
    /* These integer fixtures have exact zero residual. Different valid proof
     * kernels need not report identical conservative upper bounds. Compare
     * worker dispatch with the same row proof executed serially instead. */
    assert(upper>=0 && upper<=1.0L/100000000.0L);
    assert(kls_componentwise_certify_planned(&serial_plan,1,4,p,rows,a,tr?bt:b,x,tr,
             &reference)==1);
    assert(upper==reference);
  }
  kls_componentwise_plan_destroy(plan);
  kls_componentwise_plan_destroy(serial_plan);
  /* Exercise public certification, including strided multiple right sides.
   * The tiny matrix would not create a team itself; the fixture supplies one
   * of each family so dispatch coverage is independent of tuning thresholds. */
  s->options.threads=threads;
  unsigned long before=family?atomic_load(&s->egraph_pool->generation)
                             :atomic_load(&s->refactor_pool->generation);
  for(int tr=0;tr<2;++tr) for(int repeat=0;repeat<8;++repeat) {
    double rhs[12],answer[12];
    for(int i=0;i<12;++i) rhs[i]=answer[i]=-123;
    for(int j=0;j<2;++j) for(int i=0;i<4;++i)
      rhs[j*6+i]=(j+1)*(tr?bt[i]:b[i]);
    assert((tr?kls_solve_transpose(s,2,rhs,6,answer,6)
              :kls_solve(s,2,rhs,6,answer,6))==KLS_OK);
    for(int j=0;j<2;++j) {
      for(int i=0;i<4;++i) assert(answer[j*6+i]==j+1);
      assert(answer[j*6+4]==-123 && answer[j*6+5]==-123);
    }
    assert(!s->componentwise_work && !s->componentwise_work_context);
  }
  unsigned long after=family?atomic_load(&s->egraph_pool->generation)
                            :atomic_load(&s->refactor_pool->generation);
#if defined(_OPENMP) && defined(__SSE__) && defined(__x86_64__) && \
    FLT_RADIX == 2 && DBL_MANT_DIG == 53 && LDBL_MANT_DIG == 64 && \
    LDBL_MAX_EXP == 16384
  assert(after>before); /* At least the executor calibration actually ran. */
  assert(s->componentwise_plan);
#else
  assert(after==before); /* The stateless fallback must not dispatch workers. */
  assert(!s->componentwise_plan);
#endif
  /* Numeric work after certificate jobs must not see polluted dispatch state. */
  a[0]=4;b[0]=4;
  assert(kls_refactor(s,a)==KLS_OK);
  assert(kls_solve(s,1,b,0,x,0)==KLS_OK);
  for(int i=0;i<4;++i) assert(x[i]==1);
  kls_destroy(s);
}
static void exercise_dedicated(int family,int threads) {
  int64_t p[]={0,2,3,4,5},rows[]={0,1,1,2,3};
  double a[]={2,1,3,5,7},x[4];
  kls_options options;kls_default_options(&options);
  kls_solver *s=NULL;assert(kls_create(&s)==KLS_OK);
  assert(kls_analyze_csc(s,KLS_INDEX_INT64,4,p,rows,0,&options)==KLS_OK);
  assert(kls_factor(s,a)==KLS_OK);
  s->options.threads=threads;
  /* A declined optional team must still certify with serial partitions. */
  s->componentwise_team_declined=1;
  double initial[]={2,4,5,7};
  assert(kls_solve(s,1,initial,0,x,0)==KLS_OK);
  for(int i=0;i<4;++i) assert(x[i]==1);
  assert(!s->componentwise_team);
  assert(!s->componentwise_work && !s->componentwise_work_context);
  s->componentwise_team_declined=0;
  kls_componentwise_plan_reset(s->componentwise_plan);
  for(int generation=0;generation<4;++generation) {
    a[0]=2+generation;assert(kls_refactor(s,a)==KLS_OK);
    /* Supply only a certificate team, independently of numerical tuning. */
    destroy_refactor_pool(s);destroy_egraph_refactor_pool(s);
    for(int tr=0;tr<2;++tr) {
      double b[]={a[0],4,5,7},bt[]={a[0]+1,3,5,7};
      assert((tr?kls_solve_transpose(s,1,bt,0,x,0):kls_solve(s,1,b,0,x,0))==KLS_OK);
      for(int i=0;i<4;++i) assert(x[i]==1);
      assert(!s->componentwise_work && !s->componentwise_work_context);
    }
  }
#if defined(_OPENMP) && defined(__SSE__) && defined(__x86_64__) && \
    FLT_RADIX == 2 && DBL_MANT_DIG == 53 && LDBL_MANT_DIG == 64 && \
    LDBL_MAX_EXP == 16384
  assert(s->componentwise_team && s->componentwise_team->count==threads);
  assert(!s->componentwise_team->pending && !s->componentwise_team->work);
  exercise_dedicated_environment(s,threads,p,rows,a,0);
  exercise_dedicated_environment(s,threads,p,rows,a,1);
#else
  assert(!s->componentwise_team);
#endif
  /* If a numerical team becomes available, it must take ownership without
   * retaining another certificate crew or its cached executor choice. */
  if(family) assert(ensure_egraph_refactor_pool(s,threads));
  else assert(ensure_refactor_pool(s,threads));
  double b[]={a[0],4,5,7};
  assert(!s->componentwise_team); /* Release before starting numerical workers. */
  assert(kls_solve(s,1,b,0,x,0)==KLS_OK);
  for(int i=0;i<4;++i) assert(x[i]==1);
  assert(!s->componentwise_team);
  assert(!s->componentwise_work && !s->componentwise_work_context);
  assert(kls_analyze_csc(s,KLS_INDEX_INT64,4,p,rows,0,&options)==KLS_OK);
  assert(!s->componentwise_team && !s->componentwise_team_declined);
  kls_destroy(s);
}
static void *independent(void *opaque) {
  int family=*(int *)opaque;exercise(family,2);exercise_dedicated(family,2);return NULL;
}
static void exercise_serial(void) {
  int64_t p[]={0,2,3,4,5},rows[]={0,1,1,2,3};
  double a[]={2,1,3,5,7};
  kls_options options;kls_default_options(&options);
  kls_solver *s=NULL;assert(kls_create(&s)==KLS_OK);
  assert(kls_analyze_csc(s,KLS_INDEX_INT64,4,p,rows,0,&options)==KLS_OK);
  for(int generation=0;generation<8;++generation) {
    a[0]=2+generation;
    assert((generation?kls_refactor(s,a):kls_factor(s,a))==KLS_OK);
    for(int tr=0;tr<2;++tr) {
      double b[]={a[0],4,5,7},bt[]={a[0]+1,3,5,7},x[4];
      assert((tr?kls_solve_transpose(s,1,bt,0,x,0):kls_solve(s,1,b,0,x,0))==KLS_OK);
      for(int i=0;i<4;++i) assert(x[i]==1);
      assert(private_componentwise_pool_threads(s)==0);
#if defined(_OPENMP) && defined(__SSE__) && defined(__x86_64__) && \
    FLT_RADIX == 2 && DBL_MANT_DIG == 53 && LDBL_MANT_DIG == 64 && \
    LDBL_MAX_EXP == 16384
      assert(s->componentwise_plan);
#else
      assert(!s->componentwise_plan);
#endif
      assert(!s->componentwise_team);
    }
  }
  kls_destroy(s);
}
int main(void) {
  exercise_serial();
  for(int family=0;family<2;++family) {
    exercise(family,2);exercise(family,3);exercise(family,8);
    exercise_dedicated(family,2);exercise_dedicated(family,3);exercise_dedicated(family,8);
  }
  pthread_t workers[2];int family[]={0,1};
  for(int i=0;i<2;++i) assert(!pthread_create(workers+i,NULL,independent,family+i));
  for(int i=0;i<2;++i) assert(!pthread_join(workers[i],NULL));
  puts("PASS certificate reuse: both pools, team sizes, repeated dispatch, numeric continuation, concurrent owners");
}
