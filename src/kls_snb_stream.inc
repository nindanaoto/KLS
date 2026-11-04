/* Owned crew callback. The caller joins the complete crew before publishing
 * factors or using this stack-owned context again. */
typedef struct kls_snb_stream_job {
  kls_solver *solver;
  const double *values;
  fenv_t environment;
  int private_mode;
  _Atomic int64_t cursor;
  _Atomic int environment_failed;
  int flags[64];
  unsigned mxflags[64];
} kls_snb_stream_job;

static void kls_snb_stream_worker(void *opaque,int tid) {
  kls_snb_stream_job *job=opaque;
  kls_solver *s=job->solver;
  kls_snb_state *st=s->snb;
  fenv_t saved;
  if(tid<0 || tid>=st->nworkers || tid>=64 || fegetenv(&saved)) {
    atomic_store_explicit(&job->environment_failed,1,memory_order_release);
    atomic_store_explicit(&st->failed,1,memory_order_release);return;
  }
  int environment_failure=0;
#ifdef KLS_SNB_STREAM_TEST_ENV_FAILURE
  environment_failure=KLS_SNB_STREAM_TEST_ENV_FAILURE(tid);
#endif
  if(environment_failure || fesetenv(&job->environment)) {
    atomic_store_explicit(&job->environment_failed,1,memory_order_release);
    atomic_store_explicit(&st->failed,1,memory_order_release);
    (void)fesetenv(&saved);return;
  }
  const unsigned saved_column_mxflags=kls_snb_column_mxflags;
  const int saved_column_flags=kls_snb_column_flags;
  kls_snb_column_mxflags=0;kls_snb_column_flags=0;
  const int saved_validation=kls_snb_owned_masked_environment;
  kls_snb_owned_masked_environment=1;
  if(st->transactional) {
#if defined(__x86_64__) || defined(__i386__)
    __asm__ volatile("fnclex":::"memory");
    _mm_setcsr(_mm_getcsr()&~63u);
#endif
  }
  if(job->private_mode) {
    for(int64_t at=st->private_ptr[tid];at<st->private_ptr[tid+1];++at) {
      if(atomic_load_explicit(&st->failed,memory_order_acquire))break;
      const int32_t node=st->private_nodes[at];
      kls_snb_process_target(s,st,&st->blocks[0],node,job->values,
        s->numeric->Rs,(int)s->common.scale,st->gemm_w[tid],0);
      if(!st->transactional||!atomic_load_explicit(&st->failed,memory_order_acquire))
        atomic_store_explicit(&st->blocks[0].published[node],
          st->blocks[0].snodes[node].width,memory_order_release);
    }
  }
  const int64_t end=job->private_mode?st->private_pipeline_count:
    st->slevel_ptr[st->slevel_count];
  while(!atomic_load_explicit(&st->failed,memory_order_acquire)) {
    const int64_t task=atomic_fetch_add_explicit(&job->cursor,1,memory_order_relaxed);
    if(task>=end)break;
    const int32_t block=job->private_mode?0:st->slevel_block[task];
    const int32_t node=job->private_mode?st->private_pipeline[task]:st->slevel_snode[task];
    kls_snb_process_target(s,st,&st->blocks[block],node,job->values,
      s->numeric->Rs,(int)s->common.scale,st->gemm_w[tid],1);
  }
  job->flags[tid]=fetestexcept(FE_ALL_EXCEPT)|kls_snb_column_flags;
#if defined(__x86_64__) || defined(__i386__)
  unsigned mx;
  __asm__ volatile("stmxcsr %0":"=m"(mx)::"memory");
  job->mxflags[tid]=(mx&0x3fu)|kls_snb_column_mxflags;
#endif
  kls_snb_column_mxflags=saved_column_mxflags;
  kls_snb_column_flags=saved_column_flags;
  kls_snb_owned_masked_environment=saved_validation;
  if(fesetenv(&saved)) {
    atomic_store_explicit(&job->environment_failed,1,memory_order_release);
    atomic_store_explicit(&st->failed,1,memory_order_release);
  }
}

/* -1 declines before dispatch; -2 declines joined speculative arithmetic;
 * 0 reports failure; 1 is complete. A -2 native payload is not authoritative:
 * the caller must recompute from current A through the ordinary engine. */
static int kls_snb_stream_refactor(kls_solver *s,const double *values,int threads) {
  if(!s || !s->snb || !s->numeric || !values || threads<2 || threads>64 ||
     s->componentwise_work || s->refactor_pool || s->common.scale>0 ||
     s->numeric->Rs || s->refactor_l_packed_valid || s->refactor_u_packed_valid)
    return -1;
#if (defined(__x86_64__) || defined(__i386__)) && FLT_EVAL_METHOD==0
  unsigned mx;unsigned short cw;
  __asm__ volatile("stmxcsr %0":"=m"(mx)::"memory");
  __asm__ volatile("fnstcw %0":"=m"(cw)::"memory");
  if((mx&0x1f80u)!=0x1f80u || (cw&0x3fu)!=0x3fu ||
     fegetround()!=FE_TONEAREST)return -1;
#else
  return -1;
#endif
  kls_snb_state *st=s->snb;
#if defined(__x86_64__) || defined(__i386__)
  if(st->transactional&&(mx&0x6000u))return -1; /* nearest; exceptional speculative arithmetic rolls back. */
#endif
  atomic_store_explicit(&st->declined,0,memory_order_relaxed);
  if(st->transactional) {
    for(UF_long p=0;p<s->nnz;++p) {
      uint64_t bits;memcpy(&bits,values+p,sizeof(bits));bits&=UINT64_C(0x7fffffffffffffff);
      if(bits>=UINT64_C(0x7ff0000000000000)) {
        atomic_store_explicit(&st->declined,1,memory_order_release);return -2;
      }
    }
  }
  if(threads>st->nworkers || !st->slevel_ptr || st->slevel_count<=0)return -1;
  kls_egraph_refactor_pool *pool=s->egraph_pool;
  if(!pool || pool->thread_count!=threads || pool->created_count!=threads-1)return -1;
  kls_egraph_refactor_shared *shared=&pool->shared;
  if(pthread_mutex_trylock(&shared->lock))return -1;
  if(atomic_load_explicit(&pool->active_workers,memory_order_acquire)) {
    pthread_mutex_unlock(&shared->lock);return -1;
  }
  kls_snb_stream_job job={0};job.solver=s;job.values=values;
  job.private_mode=st->private_threads==threads && st->private_owner &&
    st->private_nodes && st->private_pipeline;
  atomic_init(&job.cursor,0);
  atomic_init(&job.environment_failed,0);
  if(fegetenv(&job.environment)) {pthread_mutex_unlock(&shared->lock);return -1;}
  for(UF_long b=0;b<st->nblocks;++b) {
    if(st->blocks[b].snode_count && !st->blocks[b].published) {
      pthread_mutex_unlock(&shared->lock);return -1;
    }
    for(int64_t k=0;k<st->blocks[b].snode_count;++k)
      atomic_store_explicit(&st->blocks[b].published[k],0,memory_order_relaxed);
  }
  int result=0;
  if(kls_snb_refactor_prelude(s,values)) {
    shared->solver=s;s->componentwise_work_context=&job;
    s->componentwise_work=kls_snb_stream_worker;
    kls_egraph_pool_dispatch_and_wait(pool,shared,threads);
    s->componentwise_work=NULL;s->componentwise_work_context=NULL;
    ++st->stream_jobs;
    if(st->transactional&&atomic_load_explicit(&st->declined,memory_order_acquire)&&
       !atomic_load_explicit(&job.environment_failed,memory_order_acquire)) {
      /* All workers have restored their entry environment and joined. Do
       * not merge any discarded arithmetic flags into the caller. Native
       * LU indices/counts and native workspaces have not been changed. */
      pthread_mutex_unlock(&shared->lock);return -2;
    }
    result=atomic_load_explicit(&st->failed,memory_order_acquire)?0:
      (int)kls_snb_refactor_epilogue(s);
    int flags=0;unsigned mxflags=0;
    for(int t=0;t<threads;++t){flags|=job.flags[t];mxflags|=job.mxflags[t];}
    if(flags)feraiseexcept(flags);
#if defined(__x86_64__) || defined(__i386__)
    unsigned mx;__asm__ volatile("stmxcsr %0":"=m"(mx)::"memory");
    mx|=mxflags;__asm__ volatile("ldmxcsr %0"::"m"(mx):"memory");
#endif
  }
  pthread_mutex_unlock(&shared->lock);return result;
}
