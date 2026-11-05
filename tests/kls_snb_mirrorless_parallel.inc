/* Four independent dense components deliberately stay in one BTF block.
 * The private schedule gives their roots to four distinct workers. Two
 * barriers force all roots into repair and all conversions to be ready
 * concurrently; each conversion must use a different owned scratch buffer. */
#ifdef KLS_SNB_HAVE_ROW_PANEL
static _Atomic int parallel_underflow_enabled;
static _Atomic uint64_t parallel_underflow_mask;
static _Atomic unsigned parallel_underflow_arrived,parallel_repair_arrived,parallel_barrier_timeout;
static double *parallel_scratch[4],*parallel_panels[4];
static void parallel_wait(_Atomic unsigned *arrived) {
  unsigned spins=0;
  while(atomic_load_explicit(arrived,memory_order_acquire)<4) {
    if(++spins==100000000u) {
      atomic_store_explicit(&parallel_barrier_timeout,1,memory_order_relaxed);break;
    }
    kls_cpu_relax();
  }
}
static int parallel_root(kls_solver *s,int64_t node) {
  int64_t start=s->snb->blocks[0].snodes[node].start;
  if(start<64||start>=384||(start-64)%96)return -1;
  return (int)((start-64)/96);
}
static void mirrorless_parallel_after_writeback(kls_solver *s,int64_t node) {
  if(!atomic_load_explicit(&parallel_underflow_enabled,memory_order_acquire))return;
  int root=parallel_root(s,node);if(root<0)return;
  const uint64_t bit=UINT64_C(1)<<(unsigned)root;
  if(atomic_fetch_or_explicit(&parallel_underflow_mask,bit,memory_order_acq_rel)&bit)return;
  atomic_fetch_add_explicit(&parallel_underflow_arrived,1,memory_order_acq_rel);
  parallel_wait(&parallel_underflow_arrived);assert(!feraiseexcept(FE_UNDERFLOW));
}
static void mirrorless_parallel_repair_begin(kls_solver *s,int64_t node,double *scratch) {
  if(!atomic_load_explicit(&parallel_underflow_enabled,memory_order_acquire))return;
  int root=parallel_root(s,node);assert(root>=0);
  parallel_scratch[root]=scratch;
  parallel_panels[root]=s->snb->blocks[0].panels+s->snb->blocks[0].snodes[node].panel;
  atomic_fetch_add_explicit(&parallel_repair_arrived,1,memory_order_acq_rel);
  parallel_wait(&parallel_repair_arrived);
}
static double parallel_truth(int column,int phase) {return 1.+.125*((column+phase)%7);}
static void parallel_rhs(const double *a,double *b,int phase) {
  memset(b,0,384*sizeof(double));
  for(int c=0;c<384;++c)for(int r=0;r<96;++r)
    b[(c/96)*96+r]+=a[c*96+r]*parallel_truth(c,phase);
}
static void parallel_matrix(double *a,double *b,int phase) {
  for(int c=0;c<384;++c)for(int r=0;r<96;++r)
    a[c*96+r]=(c%96)==r?64.+.125*phase:.03125*(1.+.015625*phase);
  parallel_rhs(a,b,phase);
}
static void parallel_reset(void) {
  atomic_store_explicit(&parallel_underflow_mask,0,memory_order_relaxed);
  atomic_store_explicit(&parallel_underflow_arrived,0,memory_order_relaxed);
  atomic_store_explicit(&parallel_repair_arrived,0,memory_order_relaxed);
  atomic_store_explicit(&parallel_barrier_timeout,0,memory_order_relaxed);
  memset(parallel_scratch,0,sizeof parallel_scratch);memset(parallel_panels,0,sizeof parallel_panels);
}
static void mirrorless_parallel_checks(void) {
  enum{N=384,NZ=N*96};
  int64_t *ap=malloc((N+1)*sizeof(*ap)),*ai=malloc(NZ*sizeof(*ai));
  double *a=malloc(NZ*sizeof(*a)),b[N],x[N];assert(ap&&ai&&a);
  for(int c=0;c<N;++c) {ap[c]=c*96;for(int r=0;r<96;++r)ai[c*96+r]=(c/96)*96+r;}
  ap[N]=NZ;parallel_matrix(a,b,0);
  atomic_store_explicit(&simulated_late_slow,0,memory_order_relaxed);
  atomic_store_explicit(&underflow_width,0,memory_order_relaxed);
  atomic_store_explicit(&force_transaction_decline,0,memory_order_relaxed);
  atomic_store_explicit(&repair_overflow,0,memory_order_relaxed);
  kls_options o;kls_default_options(&o);o.threads=4;o.use_btf=0;o.scale=0;
  o.ordering=KLS_ORDERING_NATURAL;o.orientation=KLS_ORIENTATION_TRANSPOSE;
  kls_solver *s=NULL;parallel_partition=1;
  assert(kls_create(&s)==KLS_OK);
  assert(kls_analyze_csc(s,KLS_INDEX_INT64,N,ap,ai,0,&o)==KLS_OK);
  assert(kls_factor(s,a)==KLS_OK);
  if(s->snb_auto.choice!=1)fprintf(stderr,"PARALLEL_AUTO_NOT_ARMED n=%lld nnz=%lld flops=%.17g choice=%d recipe=%d blocks=%lld observed=%d native=%d transpose=%d\n",
    (long long)s->n,(long long)s->nnz,s->common.flops,s->snb_auto.choice,s->snb_recipe,
    (long long)s->symbolic->nblocks,s->observed_parallel.choice,s->native_i16_policy.choice,s->transpose_cycle.choice);
  assert(s->snb_auto.choice==1);
  s->snb_auto.factor_credit=1.;auto_test_warm_pool(s);
  for(int epoch=1;epoch<=8;++epoch) {
    parallel_matrix(a,b,epoch);assert(kls_refactor_solve(s,a,1,b,0,x,0)==KLS_OK);
  }
  assert(s->snb_auto.choice==3&&s->snb_auto_candidate);
  kls_snb_state *st=s->snb_auto_candidate;assert(st->row_layout&&st->nworkers==4);
  assert(st->private_threads==4&&st->private_owner&&st->private_nodes&&st->private_pipeline);
  unsigned owners=0,roots=0;
  for(int64_t j=0;j<st->blocks[0].snode_count;++j) {
    kls_snb_snode *sn=st->blocks[0].snodes+j;
    if(sn->start%96==0) {
      assert(st->private_owner[j]>=0&&st->private_owner[j]<4);
      owners|=1u<<(unsigned)st->private_owner[j];++roots;
    }
  }
  assert(test_use_dead_upper_images(st)==8);
  for(int64_t j=0;j<st->blocks[0].snode_count;++j) {
    const kls_snb_snode *sn=st->blocks[0].snodes+j;
    if(sn->start%96==64)assert(sn->uw==2*(16+sn->below)&&sn->column_panel==st->blocks[0].panels+sn->panel);
  }
  assert(roots==4&&owners==15&&kls_snb_row_inplace_storage_fits(st));
  size_t capacity=N+(size_t)s->numeric->LUsize[0]+(size_t)s->numeric->Offp[N];
  double *expected=malloc(capacity*sizeof(double)),*actual=malloc(capacity*sizeof(double));
  assert(expected&&actual);
  for(int epoch=9;epoch<=11;++epoch) {
    parallel_matrix(a,b,epoch);
    if(epoch<11) {
      assert(kls_refactor_solve(s,a,1,b,0,x,0)==KLS_OK);
      assert(s->stats.last_refactor_path==KLS_REFACTOR_PATH_SNB);
    }
    assert(transaction_capture(s,expected)==capacity);
    parallel_rhs(a,b,epoch+3);for(int r=0;r<N;++r)x[r]=733.;
    s->unchanged_refactor_state=-1;s->refactor_input_snapshot_valid=0;
    parallel_reset();
    const unsigned before=atomic_load_explicit(&row_inplace_repairs,memory_order_relaxed);
    const unsigned packed_before=atomic_load_explicit(&packed_row_calls,memory_order_relaxed);
    if(epoch==11)st->gemm_scratch_bytes=0;
    assert(!feclearexcept(FE_ALL_EXCEPT));
    atomic_store_explicit(&parallel_underflow_enabled,1,memory_order_release);
    assert(kls_refactor_solve(s,a,1,b,0,x,0)==KLS_OK);
    atomic_store_explicit(&parallel_underflow_enabled,0,memory_order_release);
    const unsigned csr=_mm_getcsr();const int flags=fetestexcept(FE_ALL_EXCEPT);
    assert(!fetestexcept(FE_UNDERFLOW|FE_INVALID|FE_OVERFLOW|FE_DIVBYZERO));
    assert(atomic_load_explicit(&parallel_underflow_mask,memory_order_acquire)==15);
    assert(atomic_load_explicit(&parallel_underflow_arrived,memory_order_acquire)==4);
    assert(!atomic_load_explicit(&parallel_barrier_timeout,memory_order_relaxed));
    assert(!atomic_load_explicit(&s->egraph_pool->active_workers,memory_order_acquire));
    assert(!s->componentwise_work&&!s->componentwise_work_context);
    for(int r=0;r<N;++r)assert(fabs(x[r]-parallel_truth(r,epoch+3))<1e-11);
    assert(transaction_capture(s,actual)==capacity);
    if(epoch<11) {
      assert(s->stats.last_refactor_path==KLS_REFACTOR_PATH_SNB);
      assert(atomic_load_explicit(&parallel_repair_arrived,memory_order_acquire)==4);
      assert(atomic_load_explicit(&row_inplace_repairs,memory_order_relaxed)==before+4);
      assert(atomic_load_explicit(&packed_row_calls,memory_order_relaxed)==packed_before+16);
      assert(!memcmp(expected,actual,capacity*sizeof(double)));
      for(int i=0;i<4;++i)for(int j=0;j<i;++j) {
        assert(parallel_scratch[i]&&parallel_scratch[i]!=parallel_scratch[j]);
        assert(parallel_panels[i]&&parallel_panels[i]!=parallel_panels[j]);
      }
    } else {
      assert(!s->snb_auto_candidate&&s->snb_auto.choice<0);
      assert(s->stats.last_refactor_path!=KLS_REFACTOR_PATH_SNB);
      assert(atomic_load_explicit(&row_inplace_repairs,memory_order_relaxed)==before);
      assert(memcmp(expected,actual,capacity*sizeof(double))); /* Current A changed. */
      s->unchanged_refactor_state=-1;s->refactor_input_snapshot_valid=0;
      parallel_rhs(a,b,epoch+5);assert(!feclearexcept(FE_ALL_EXCEPT));
      assert(kls_refactor_solve(s,a,1,b,0,x,0)==KLS_OK);
      assert(s->stats.last_refactor_path!=KLS_REFACTOR_PATH_UNCHANGED);
      assert(transaction_capture(s,expected)==capacity);
      assert(!memcmp(expected,actual,capacity*sizeof(double)));
      assert(_mm_getcsr()==csr&&fetestexcept(FE_ALL_EXCEPT)==flags);
      for(int r=0;r<N;++r)assert(fabs(x[r]-parallel_truth(r,epoch+5))<1e-11);
    }
  }
  free(expected);free(actual);kls_destroy(s);parallel_partition=0;free(ap);free(ai);free(a);
  puts("MIRRORLESS_PARALLEL_ORACLES workers=4 repairs=8 packed_producer_edges=32 alias_repairs=8 exact_dead_upper_boundary=1 distinct_scratch=1 fresh_rhs=1 joined_scratch_decline=1 current_A_native_bits_flags=1");
}
#else
static void mirrorless_parallel_after_writeback(kls_solver *s,int64_t node) {(void)s;(void)node;}
static void mirrorless_parallel_repair_begin(kls_solver *s,int64_t node,double *scratch) {(void)s;(void)node;(void)scratch;}
static void mirrorless_parallel_checks(void) {}
#endif
