#define _GNU_SOURCE
#define _POSIX_C_SOURCE 200809L
#include <time.h>
#include <pthread.h>
#include <stdatomic.h>
#undef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
static atomic_int expire_worker;
static pthread_t caller_thread;
static int fixture_clock_gettime(clockid_t clock,struct timespec *ts);
#define clock_gettime fixture_clock_gettime
#include "../src/kls.c"
#undef clock_gettime
#undef NDEBUG
#include <assert.h>

static int fixture_clock_gettime(clockid_t clock,struct timespec *ts) {
  int result=clock_gettime(clock,ts);
  if(!result && atomic_load_explicit(&expire_worker,memory_order_relaxed) &&
      !pthread_equal(caller_thread,pthread_self())) ts->tv_sec+=1000;
  return result;
}

static void exercise(int pipeline,int csr) {
  enum {N=256,NZ=3*N};int64_t p[N+1],rows[NZ];double a[NZ];
  double truth[N],rhs[N],answer[N],transpose_rhs[N];int nz=0;
  for(int j=0;j<N;++j) {
    p[j]=nz;
    for(int i=0;i<N;++i) {
      int linked=pipeline?((i/16==j/16)&&(i%16==15 || j%16==15))
        :((i&~1)==(j&~1));
      if(i==j || linked) {rows[nz]=i;a[nz++]=i==j?4:((i+2*j)%3+1)*0.015625;}
    }
    truth[j]=1+(j%7)*0.0625;
  }
  p[N]=nz;assert(nz<=NZ);
  kls_options o;kls_default_options(&o);o.threads=1;o.use_btf=0;
  o.scale=0;o.static_pivoting=0;o.ordering=KLS_ORDERING_NATURAL;
  o.orientation=KLS_ORIENTATION_NORMAL;o.expected_refactorizations=o.expected_solves=0;
  kls_solver *s=NULL;assert(kls_create(&s)==KLS_OK);
  assert((csr?kls_analyze_csr:kls_analyze_csc)(s,KLS_INDEX_INT64,N,p,rows,0,&o)==KLS_OK);
  assert(kls_factor(s,a)==KLS_OK);
  s->options.threads=8;s->tuning.egraph_min_size=0;
  s->tuning.egraph_min_flops_per_thread=1e12;
  s->lean_choice=-1;s->lean_probe_arm=0;s->solve_contract_probe=2;
  s->factor_preps_deferred=0;s->rowmatch_deferred=0;
  s->unchanged_refactor_state=-1;
  assert(kls_build_refactor_map(s));
  kls_set_last_refactor_path(s,KLS_REFACTOR_PATH_MAPPED);
  assert(kls_observed_parallel_capable(s));
  assert(!kls_build_refactor_schedule(s));
  for(int i=0;i<N;++i) rhs[i]=transpose_rhs[i]=0;
  for(int j=0;j<N;++j) for(int k=p[j];k<p[j+1];++k) {
    rhs[rows[k]]+=a[k]*truth[j];transpose_rhs[j]+=a[k]*truth[rows[k]];
  }
  /* Deterministic activation: seed observed costs only in this fixture.
     No assertion depends on the noisy measured winner. */
  s->observed_parallel.cycles=200;s->observed_parallel.total=200;
  for(int i=0;i<8;++i) s->observed_parallel.history[i]=1;
  a[0]+=0.125;rhs[0]+=0.125*truth[0];transpose_rhs[0]+=0.125*truth[0];
  assert(kls_refactor_solve(s,a,1,csr?transpose_rhs:rhs,0,answer,0)==KLS_OK);
  assert(s->observed_parallel.spent>0);
  assert(s->observed_parallel.candidate[0]>0);
  assert(!s->observed_parallel.active && !s->observed_parallel.deadline);
  for(int j=0;j<N;++j) assert(fabs(answer[j]-truth[j])<1e-11);
  assert(s->stats.last_refactor_path==KLS_REFACTOR_PATH_MAPPED);
  /* Explicit hints and loaded profiles retain their existing policy. */
  s->options.expected_refactorizations=1;assert(!kls_observed_parallel_capable(s));
  s->options.expected_refactorizations=0;s->stats.tuning_profile_active=1;
  assert(!kls_observed_parallel_capable(s));s->stats.tuning_profile_active=0;
  s->rowmatch_deferred=1;assert(!kls_observed_parallel_capable(s));
  s->rowmatch_deferred=0;
  s->accuracy_policy=KLS_ACCURACY_APPLICATION_MANAGED;
  assert(!kls_observed_parallel_capable(s));
  s->accuracy_policy=KLS_ACCURACY_COMPONENTWISE_BACKWARD_ERROR;
  /* Invalidation resets the one-attempt policy without relying on graph free. */
  unsigned epoch=s->observed_parallel_epoch;kls_invalidate_i32_solve(s);
  assert(epoch!=s->observed_parallel_epoch && s->observed_parallel.choice==0);
  s->observed_parallel.active=1;
  assert(kls_build_refactor_schedule(s));assert(kls_egraph_refactor_is_eligible(s));
  assert((s->refactor_cluster_level_count<s->refactor_level_count)==pipeline);
  s->observed_parallel.deadline=kls_now_seconds()-1;
  free_refactor_schedule(s);assert(!kls_build_refactor_schedule(s));
  assert(!s->refactor_level_ptr);
  s->observed_parallel.deadline=0;assert(kls_build_refactor_schedule(s));
  /* The caller's graph preparation succeeds; workers deterministically see
     expiry during real cluster/pipeline dispatch. All barriers must quiesce. */
  s->observed_parallel.deadline=kls_now_seconds()+10;
  atomic_store_explicit(&expire_worker,1,memory_order_relaxed);
  double *values=NULL;assert(prepare_numeric_values(s,a,&values)==KLS_OK);
  int stopped=kls_egraph_mapped_refactor(s,values,0);
  atomic_store_explicit(&expire_worker,0,memory_order_relaxed);
  assert(stopped==-1 && s->egraph_pool);
  assert(atomic_load_explicit(&s->egraph_pool->shared.trial_deadline_hit,memory_order_acquire));
  assert(!s->egraph_pool->shared.trial_deadline_seconds && s->egraph_worker_scratch_dirty);
  assert(!atomic_load_explicit(&s->egraph_pool->active_workers,memory_order_acquire));
  s->observed_parallel=(kls_observed_parallel_state){.choice=-1};
  /* Refactor and certify the real current matrix after partial publication. */
  assert(kls_refactor_solve(s,a,1,csr?transpose_rhs:rhs,0,answer,0)==KLS_OK);
  assert(kls_solve_transpose(s,1,csr?rhs:transpose_rhs,0,answer,0)==KLS_OK);
  for(int j=0;j<N;++j) assert(fabs(answer[j]-truth[j])<1e-11);
  memcpy(answer,csr?transpose_rhs:rhs,sizeof(answer));
  assert(kls_refactor_solve(s,a,1,answer,0,answer,0)==KLS_OK);
  for(int j=0;j<N;++j) assert(fabs(answer[j]-truth[j])<1e-11);
  double multi_rhs[2*N],multi_answer[2*N];
  for(int j=0;j<N;++j) {
    multi_rhs[j]=(csr?transpose_rhs:rhs)[j];multi_rhs[N+j]=2*multi_rhs[j];
  }
  assert(kls_solve(s,2,multi_rhs,0,multi_answer,0)==KLS_OK);
  for(int j=0;j<N;++j) {
    assert(fabs(multi_answer[j]-truth[j])<1e-11);
    assert(fabs(multi_answer[N+j]-2*truth[j])<1e-11);
  }
  /* A joined first-job deadline defers once; it cannot permanently reject
     a fast representation before later changed-input funding arrives. */
  s->observed_parallel=(kls_observed_parallel_state){.cycles=200,.total=200};
  for(int i=0;i<8;++i)s->observed_parallel.history[i]=1;
  kls_set_last_refactor_path(s,KLS_REFACTOR_PATH_MAPPED);
  assert(kls_observed_parallel_capable(s));
  a[0]+=.125;rhs[0]+=.125*truth[0];transpose_rhs[0]+=.125*truth[0];
  atomic_store_explicit(&expire_worker,1,memory_order_relaxed);
  assert(kls_refactor_solve(s,a,1,csr?transpose_rhs:rhs,0,answer,0)==KLS_OK);
  atomic_store_explicit(&expire_worker,0,memory_order_relaxed);
  assert(s->observed_parallel.choice==0 && s->observed_parallel.retry_pending &&
    s->observed_parallel.deadline_retries==1 && s->observed_parallel.spent>0);
  assert(s->egraph_pool && !atomic_load_explicit(&s->egraph_pool->active_workers,memory_order_acquire));
  const double saved_floor=s->tuning.egraph_min_flops_per_thread;
  s->tuning.egraph_min_flops_per_thread=0;
  assert(!kls_egraph_refactor_is_eligible(s)); /* Cached plan cannot bypass funding. */
  s->observed_parallel.active=1;assert(kls_egraph_refactor_is_eligible(s));
  s->observed_parallel.active=0;s->tuning.egraph_min_flops_per_thread=saved_floor;
  for(int j=0;j<N;++j)assert(fabs(answer[j]-truth[j])<1e-11);
  a[0]+=.125;rhs[0]+=.125*truth[0];transpose_rhs[0]+=.125*truth[0];
  assert(kls_refactor_solve(s,a,1,csr?transpose_rhs:rhs,0,answer,0)==KLS_OK);
  assert(s->observed_parallel.deadline_retries==1 && s->observed_parallel.candidate[0]>0);
  assert(!s->observed_parallel.active && !s->observed_parallel.deadline);
  for(int j=0;j<N;++j)assert(fabs(answer[j]-truth[j])<1e-11);
  /* A failed speculative epoch cannot turn an old factor/answer into
     success, and a failed public combined call must not modify output. */
  s->observed_parallel=(kls_observed_parallel_state){.cycles=200,.total=200};
  for(int i=0;i<8;++i) s->observed_parallel.history[i]=1;
  kls_set_last_refactor_path(s,KLS_REFACTOR_PATH_MAPPED);
  assert(kls_observed_parallel_capable(s));
  for(int k=0;k<nz;++k) a[k]=0;
  for(int j=0;j<N;++j) answer[j]=-123.5;
  atomic_store_explicit(&expire_worker,1,memory_order_relaxed);
  assert(kls_refactor_solve(s,a,1,rhs,0,answer,0)!=KLS_OK);
  atomic_store_explicit(&expire_worker,0,memory_order_relaxed);
  assert(s->observed_parallel.choice<0 && !s->observed_parallel.retry_pending);
  for(int j=0;j<N;++j) assert(answer[j]==-123.5);
  assert(!s->observed_parallel.active && !s->observed_parallel.deadline);
  kls_destroy(s);
}

int main(void) {
  caller_thread=pthread_self();atomic_init(&expire_worker,0);
  for(int pipeline=0;pipeline<2;++pipeline) for(int csr=0;csr<2;++csr)
    exercise(pipeline,csr);
  return 0;
}
