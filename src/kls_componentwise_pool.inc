/* Reuse a complete numerical team, without entering its numerical dispatcher
 * or constructing another runtime. Same-instance public calls are serialized
 * by the caller, as required for all existing solver operations. */
static int private_componentwise_pool_threads(const kls_solver *s) {
  const int ordinary=s->refactor_pool && s->refactor_pool->created_count>0;
  const int egraph=s->egraph_pool && s->egraph_pool->created_count>0;
  if(ordinary==egraph) return 0; /* No unambiguous existing team. */
  if(ordinary) {
    const kls_refactor_pool *pool=s->refactor_pool;
    return pool->created_count==pool->thread_count?pool->thread_count:0;
  }
  const kls_egraph_refactor_pool *pool=s->egraph_pool;
  return pool->created_count==pool->thread_count-1?pool->thread_count:0;
}

static int private_componentwise_execute(void *opaque,
    void (*work)(void *,int),void *job,int threads) {
  kls_solver *s=opaque;
  if(!s || !work || threads<=1 || s->componentwise_work ||
     private_componentwise_pool_threads(s)!=threads) return 0;
  if(s->refactor_pool && s->refactor_pool->created_count>0) {
    kls_refactor_pool *pool=s->refactor_pool;
    kls_parallel_refactor_shared *shared=&pool->shared;
    pthread_mutex_lock(&shared->lock);
    if(atomic_load_explicit(&pool->active_workers,memory_order_acquire)!=0 ||
       atomic_load_explicit(&pool->shutdown,memory_order_acquire)) {
      pthread_mutex_unlock(&shared->lock);return 0;
    }
    kls_solver *previous=shared->solver;
    shared->solver=s;s->componentwise_work=work;s->componentwise_work_context=job;
    atomic_store_explicit(&pool->active_workers,threads,memory_order_relaxed);
    atomic_fetch_add_explicit(&pool->generation,1ul,memory_order_release);
    pthread_cond_broadcast(&pool->work_cond);
    while(atomic_load_explicit(&pool->active_workers,memory_order_acquire)>0)
      pthread_cond_wait(&pool->done_cond,&shared->lock);
    s->componentwise_work=NULL;s->componentwise_work_context=NULL;
    shared->solver=previous;
    pthread_mutex_unlock(&shared->lock);
  } else {
    kls_egraph_refactor_pool *pool=s->egraph_pool;
    kls_egraph_refactor_shared *shared=&pool->shared;
    pthread_mutex_lock(&shared->lock);
    if(atomic_load_explicit(&pool->active_workers,memory_order_acquire)!=0 ||
       atomic_load_explicit(&pool->shutdown,memory_order_acquire)) {
      pthread_mutex_unlock(&shared->lock);return 0;
    }
    kls_solver *previous=shared->solver;
    const int old_threads=shared->thread_count;
    shared->solver=s;shared->thread_count=threads;
    s->componentwise_work=work;s->componentwise_work_context=job;
    kls_egraph_pool_dispatch_and_spin_wait(pool,shared,threads);
    s->componentwise_work=NULL;s->componentwise_work_context=NULL;
    shared->solver=previous;shared->thread_count=old_threads;
    pthread_mutex_unlock(&shared->lock);
  }
  return 1;
}

/* Called lazily by a supported row certificate's parallel trial. Unsupported
 * environments and formats never construct this team. Creation failure uses
 * serial partitions of the same proof and lets the ordinary timing choice
 * select serial rows; it must not turn an optional accelerator into failure. */
static int private_componentwise_execute_dedicated(void *opaque,
    void (*work)(void *,int),void *job,int threads) {
  kls_solver *s=opaque;
  if(!s || !work || threads<=1 || threads!=s->options.threads ||
     s->componentwise_work ||
     (s->refactor_pool && s->refactor_pool->created_count>0) ||
     (s->egraph_pool && s->egraph_pool->created_count>0)) return 0;
  if(s->componentwise_team && s->componentwise_team->count!=threads) {
    kls_componentwise_team_destroy(s->componentwise_team);
    s->componentwise_team=NULL;s->componentwise_team_declined=0;
  }
  if(!s->componentwise_team && !s->componentwise_team_declined) {
    s->componentwise_team=kls_componentwise_team_create(threads);
    if(!s->componentwise_team) s->componentwise_team_declined=1;
  }
  s->componentwise_work=work;s->componentwise_work_context=job;
  int completed=1;
  if(s->componentwise_team)
    completed=kls_componentwise_team_execute(s->componentwise_team,work,job,threads);
  else
    for(int tid=0;tid<threads;++tid) work(job,tid);
  s->componentwise_work=NULL;s->componentwise_work_context=NULL;
  return completed;
}
