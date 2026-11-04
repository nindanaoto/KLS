#define _GNU_SOURCE
#define _POSIX_C_SOURCE 200809L
#define KLS_SNB_AUTO_TEST_SECONDS(arm,actual) ((arm)?(atomic_load_explicit(&simulated_late_slow,memory_order_relaxed)?.004:.001):.002)
#include <stdatomic.h>
#include <fenv.h>
#include "../include/kls/kls.h"
static void *row_test_build(kls_solver *,double,size_t);
#define KLS_SNB_AUTO_TEST_BUILD(s,deadline,limit) row_test_build(s,deadline,limit)
static _Atomic unsigned row_edges,row_getrfs,row_small_edges;
#define KLS_RPANEL_TEST_EDGE() ((void)atomic_fetch_add_explicit(&row_edges,1,memory_order_relaxed))
#define KLS_RPANEL_TEST_GETRF() ((void)atomic_fetch_add_explicit(&row_getrfs,1,memory_order_relaxed))
#define KLS_RPANEL_TEST_SMALL_EDGE() ((void)atomic_fetch_add_explicit(&row_small_edges,1,memory_order_relaxed))
static _Atomic int force_transaction_decline,transaction_native_changed,simulated_late_slow;
static _Atomic unsigned transaction_declines;
static _Atomic int underflow_width;
static _Atomic unsigned underflow_injections,row_underflow_repairs;
static _Atomic int repair_overflow;
#define KLS_SNB_ROW_UNDERFLOW_REPAIR() do { \
 atomic_fetch_add_explicit(&row_underflow_repairs,1,memory_order_relaxed); \
 if(atomic_exchange_explicit(&repair_overflow,0,memory_order_relaxed))feraiseexcept(FE_OVERFLOW); \
} while(0)
static _Atomic unsigned transaction_writebacks;
static double transaction_old_diagonal[129];
#define KLS_SNB_TRANSACTION_AFTER_WRITEBACK(s,node) do { \
  atomic_fetch_add_explicit(&transaction_writebacks,1,memory_order_relaxed); \
  int uf_width=atomic_load_explicit(&underflow_width,memory_order_relaxed); \
  if(uf_width&&(s)->snb->blocks[0].snodes[node].width==uf_width&& \
     atomic_exchange_explicit(&underflow_width,0,memory_order_relaxed)) { \
    atomic_fetch_add_explicit(&underflow_injections,1,memory_order_relaxed); \
    feraiseexcept(FE_UNDERFLOW); \
  } \
  if((!(s)->snb->row_layout||(s)->snb->blocks[0].snodes[node].width==16) && \
     atomic_exchange_explicit(&force_transaction_decline,0,memory_order_relaxed)) { \
    size_t k=(size_t)(s)->snb->blocks[0].snodes[node].start; \
    atomic_store_explicit(&transaction_native_changed, \
      memcmp((double *)(s)->numeric->Udiag+k,transaction_old_diagonal+k,sizeof(double))!=0,memory_order_relaxed); \
    atomic_fetch_add_explicit(&transaction_declines,1,memory_order_relaxed); \
    feraiseexcept(FE_OVERFLOW); \
  } \
} while(0)
#include "../src/kls.c"
#undef NDEBUG
#include <assert.h>
static void *row_test_build(kls_solver *s,double deadline,size_t limit) {
 assert(s->symbolic->nblocks==1);int64_t *starts=malloc((size_t)s->n*sizeof(*starts));assert(starts);
 const int widths[]={3,16,2,16,1,16};int64_t count=0,col=0;
 while(col<(int64_t)s->n){starts[count]=col;col+=widths[count%6];++count;}
 kls_snb_partition p={starts,&count};void *st=kls_snb_build(s,&p,deadline,limit);free(starts);return st;
}
static size_t transaction_capture(kls_solver *s,double *out) {
 size_t at=0;trilinos_klu_l_numeric *n=s->numeric;
 memcpy(out,n->Udiag,(size_t)s->n*sizeof(double));at+=s->n;
 for(UF_long b=0;b<s->symbolic->nblocks;++b) {
  memcpy(out+at,n->LUbx[b],(size_t)n->LUsize[b]*sizeof(double));at+=(size_t)n->LUsize[b];
 }
 memcpy(out+at,n->Offx,(size_t)n->Offp[s->n]*sizeof(double));return at+(size_t)n->Offp[s->n];
}
static void exercise(int transpose,int transactional) {
 unsigned before_getrfs=atomic_load_explicit(&row_getrfs,memory_order_relaxed);
 unsigned before_edges=atomic_load_explicit(&row_edges,memory_order_relaxed);
 unsigned before_small=atomic_load_explicit(&row_small_edges,memory_order_relaxed);
 enum{N=129};int64_t ap[N+1],ai[N*N];double a[N*N],b[N],x[N];
 for(int c=0;c<N;++c){ap[c]=c*N;for(int r=0;r<N;++r){ai[c*N+r]=r;a[c*N+r]=r==c?N+1.:.03125;}}ap[N]=N*N;
 kls_options o;kls_default_options(&o);o.threads=2;o.use_btf=0;o.scale=0;
 o.ordering=KLS_ORDERING_NATURAL;o.orientation=transpose?KLS_ORIENTATION_TRANSPOSE:KLS_ORIENTATION_NORMAL;
 kls_solver *s=NULL;assert(kls_create(&s)==KLS_OK);
 assert(kls_analyze_csc(s,KLS_INDEX_INT64,N,ap,ai,0,&o)==KLS_OK);
 assert(kls_factor(s,a)==KLS_OK);assert(s->snb_auto.choice==1);
 /* These are funded arithmetic/ownership tests, not wall-clock admission
    tests. Instrumentation and scheduler delays must not skip their kernels;
    the policy test separately checks fixed setup-credit boundaries. */
 s->snb_auto.factor_credit=1.;
 const kls_snb_auto_policy prepared=s->snb_auto;
 for(int other=0;other<3;++other) {
  int *active=other==0?&s->transpose_cycle.active:
    other==1?&s->observed_parallel.active:&s->native_i16_policy.active;
  assert(!*active);*active=1;s->snb_auto.choice=0;
  double elapsed=0;assert(!kls_snb_auto_factor_prepare(s,&elapsed));
  assert(s->snb_auto.choice==0);*active=0;s->snb_auto=prepared;
 }
 assert(!kls_compact_transpose_cycle_capable(s));
 kls_observed_parallel_state observed=s->observed_parallel;
 (void)kls_observed_parallel_before(s,NULL,1,NULL,NULL);
 assert(!memcmp(&observed,&s->observed_parallel,sizeof observed));

 for(int epoch=1;epoch<=12;++epoch) {
  memset(b,0,sizeof b);
  for(int c=0;c<N;++c)for(int r=0;r<N;++r){
   a[c*N+r]=r==c?N+1.+.125*epoch:.03125*(1+.0625*epoch);
   b[r]+=a[c*N+r]*(1.+.125*((c+epoch)%7));
  }
  assert(kls_refactor_solve(s,a,1,b,0,x,0)==KLS_OK);
  for(int c=0;c<N;++c)assert(fabs(x[c]-(1.+.125*((c+epoch)%7)))<1e-11);
  if(epoch<=8) {
   if(s->snb_auto.phase!=epoch)fprintf(stderr,"AUTO_TEST activation epoch=%d phase=%d choice=%d spent=%.9g limit=%.9g transpose=%d transactional=%d\n",epoch,s->snb_auto.phase,s->snb_auto.choice,s->snb_auto.spent,s->snb_auto.limit,transpose,transactional);
   assert(s->snb_auto.phase==epoch);
  }
 }
 assert(s->snb_auto.choice==3 && s->snb_auto_candidate);
 assert(s->snb_auto_candidate->stream_jobs>=4);
 if(!transpose)assert(s->native_i16_normal_runs>=4);
#ifdef KLS_SNB_HAVE_ROW_PANEL
 assert(s->snb_auto_candidate->row_layout&&s->snb_auto_candidate->transactional);
 assert(atomic_load_explicit(&row_getrfs,memory_order_relaxed)>before_getrfs);
 assert(atomic_load_explicit(&row_edges,memory_order_relaxed)>before_edges);
 assert(atomic_load_explicit(&row_small_edges,memory_order_relaxed)>before_small);
 unsigned mirrors=0; kls_snb_block *mirror_block=s->snb_auto_candidate->blocks;
 for(int64_t j=0;j<mirror_block->snode_count;++j) {
  const kls_snb_snode *node=mirror_block->snodes+j;if(!node->column_panel)continue;
  assert(node->width==16);++mirrors;
  for(int c=0;c<16;++c)for(int r=0;r<node->height;++r)
   assert(!memcmp(node->column_panel+(size_t)c*node->height+r,mirror_block->panels+node->panel+kls_rpanel_at(node->height,c,r),8));
 }
 assert(mirrors>0);
 assert(kls_snb_enable_row_layout(s->snb_auto_candidate,0,SIZE_MAX)); /* Idempotent. */
#else
 (void)before_getrfs;(void)before_edges;(void)before_small;
#endif
 if(transactional==3) {
  const double original_limit=s->snb_auto.limit;
  s->snb_auto.factor_credit+=100.; /* A recheck must not mint new credit. */
  unsigned monitor_count=s->snb_auto.retained_count,monitor_cursor=s->snb_auto.retained_cursor;
  atomic_store_explicit(&simulated_late_slow,1,memory_order_relaxed);
  double other_rhs[N];for(int r=0;r<N;++r)other_rhs[r]=2*b[r];
  assert(kls_refactor_solve(s,a,1,other_rhs,0,x,0)==KLS_OK);
  for(int c=0;c<N;++c)assert(fabs(x[c]-2*(1.+.125*((c+12)%7)))<2e-11);
  assert(!s->snb_auto.rechecked&&s->snb_auto.retained_count==monitor_count&&s->snb_auto.retained_cursor==monitor_cursor);
  for(int epoch=13;epoch<=48;++epoch) {
   memset(b,0,sizeof b);
   for(int c=0;c<N;++c)for(int r=0;r<N;++r) {
    a[c*N+r]=r==c?N+1.+.125*epoch:.03125*(1+.0625*epoch);
    b[r]+=a[c*N+r]*(1.+.125*((c+epoch)%7));
   }
   assert(kls_refactor_solve(s,a,1,b,0,x,0)==KLS_OK);
   for(int c=0;c<N;++c)assert(fabs(x[c]-(1.+.125*((c+epoch)%7)))<1e-11);
   assert(s->snb_auto.limit==original_limit);
   assert(!s->componentwise_work&&!s->componentwise_work_context);
   if(s->egraph_pool)assert(!atomic_load_explicit(&s->egraph_pool->active_workers,memory_order_acquire));
  }
  assert(s->snb_auto.rechecked&&s->snb_auto.phase==8&&s->snb_auto.choice<0);
  assert(s->snb_auto.count[0]==4&&s->snb_auto.count[1]==4&&!s->snb_auto_candidate);
  atomic_store_explicit(&simulated_late_slow,0,memory_order_relaxed);
  kls_destroy(s);puts("PASS late certified ABBA recheck, unchanged-A exclusion, fixed credit and fresh RHS epochs");return;
 }
 if(transactional) {
  kls_snb_state *candidate=s->snb_auto_candidate,*selected=s->snb;
  candidate->transactional=1;s->snb=candidate;
  double *prepared_values=NULL;assert(prepare_numeric_values(s,a,&prepared_values)==KLS_OK);
  double *bad=malloc(N*N*sizeof(double));size_t units=(size_t)s->numeric->LUsize[0];
  double *old_lu=malloc(units*sizeof(double)),old_work[N];assert(bad&&old_lu);
  memcpy(bad,prepared_values,N*N*sizeof(double));
  const uint64_t nan_bits=UINT64_C(0x7ff8000000000042);memcpy(bad,&nan_bits,8);
  memcpy(old_lu,s->numeric->LUbx[0],units*sizeof(double));
  memcpy(old_work,s->numeric->Xwork,sizeof old_work);
  memcpy(transaction_old_diagonal,s->numeric->Udiag,sizeof transaction_old_diagonal);
  fenv_t environment;assert(!fegetenv(&environment));assert(!feclearexcept(FE_ALL_EXCEPT));
  assert(kls_snb_stream_refactor(s,bad,2)==-2);
  assert(!fetestexcept(FE_ALL_EXCEPT));
  assert(!memcmp(old_lu,s->numeric->LUbx[0],units*sizeof(double)));
  assert(!memcmp(old_work,s->numeric->Xwork,sizeof old_work));
  assert(!memcmp(transaction_old_diagonal,s->numeric->Udiag,sizeof transaction_old_diagonal));
  assert(!atomic_load_explicit(&s->egraph_pool->active_workers,memory_order_acquire));
  assert(!s->componentwise_work&&!s->componentwise_work_context);
  unsigned control=_mm_getcsr();const unsigned unsupported[]={0x2000u};
  for(unsigned i=0;i<sizeof(unsupported)/sizeof(*unsupported);++i) {
   unsigned selected=(control&~0xe040u)|unsupported[i];_mm_setcsr(selected);
   assert(kls_snb_stream_refactor(s,prepared_values,2)==-1);
   assert(_mm_getcsr()==selected);
   assert(!memcmp(old_lu,s->numeric->LUbx[0],units*sizeof(double)));
  }
  _mm_setcsr(control);
  assert(kls_snb_refactor(s,prepared_values)==0);
  assert(kls_snb_refactor_parallel(s,prepared_values,2)==0);
  assert(!memcmp(old_lu,s->numeric->LUbx[0],units*sizeof(double)));
#ifdef KLS_SNB_HAVE_ROW_PANEL
  for(int which=0;which<3;++which) {
   assert(!feclearexcept(FE_ALL_EXCEPT));
   unsigned injections=atomic_load_explicit(&underflow_injections,memory_order_relaxed);
   unsigned repairs=atomic_load_explicit(&row_underflow_repairs,memory_order_relaxed);
   atomic_store_explicit(&repair_overflow,which==2,memory_order_relaxed);
   atomic_store_explicit(&underflow_width,which?16:3,memory_order_relaxed);
   assert(kls_snb_stream_refactor(s,prepared_values,2)==(which==2?-2:1));
   assert(!atomic_load_explicit(&repair_overflow,memory_order_relaxed));
   assert(atomic_load_explicit(&row_underflow_repairs,memory_order_relaxed)==repairs+(unsigned)(which!=0));
   assert(atomic_load_explicit(&underflow_injections,memory_order_relaxed)==injections+1);
   assert(!atomic_load_explicit(&underflow_width,memory_order_relaxed));
   assert(!!fetestexcept(FE_UNDERFLOW)==!which);
   assert(!fetestexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW));
   assert(!memcmp(old_lu,s->numeric->LUbx[0],units*sizeof(double)));
   assert(!memcmp(transaction_old_diagonal,s->numeric->Udiag,sizeof transaction_old_diagonal));
   assert(!atomic_load_explicit(&s->egraph_pool->active_workers,memory_order_acquire));
   assert(!s->componentwise_work&&!s->componentwise_work_context);
  }
#endif
  assert(!fesetenv(&environment));free(bad);free(old_lu);s->snb=selected;
  if(transactional==2) {
   unsigned before=atomic_load_explicit(&transaction_writebacks,memory_order_relaxed);
#ifdef KLS_SNB_HAVE_ROW_PANEL
   unsigned row_before=atomic_load_explicit(&row_getrfs,memory_order_relaxed);
#endif
   for(int c=0;c<N;++c)for(int r=0;r<N;++r)a[c*N+r]=r==c&&r!=3?1.:0.;
   for(int r=0;r<N;++r){b[r]=1.;x[r]=733.;}
   assert(kls_refactor_solve(s,a,1,b,0,x,0)!=KLS_OK);
   for(int r=0;r<N;++r)assert(x[r]==733.);
   assert(atomic_load_explicit(&transaction_writebacks,memory_order_relaxed)>before);
#ifdef KLS_SNB_HAVE_ROW_PANEL
   assert(atomic_load_explicit(&row_getrfs,memory_order_relaxed)>row_before);
#endif
   assert(!s->snb_auto_candidate&&s->snb_auto.choice<0);
   assert(!atomic_load_explicit(&s->egraph_pool->active_workers,memory_order_acquire));
   assert(!s->componentwise_work&&!s->componentwise_work_context);
   kls_destroy(s);puts("PASS singular row-panel decline and unchanged failed-output staging");return;
  }
  atomic_store_explicit(&transaction_declines,0,memory_order_relaxed);
  atomic_store_explicit(&transaction_native_changed,0,memory_order_relaxed);
  atomic_store_explicit(&force_transaction_decline,1,memory_order_relaxed);
  size_t capacity=N+(size_t)s->numeric->LUsize[0]+(size_t)s->numeric->Offp[N];
  double *expected_native=malloc(capacity*sizeof(double)),*actual_native=malloc(capacity*sizeof(double));
  assert(expected_native&&actual_native);
  for(int epoch=13;epoch<=14;++epoch) {
   memset(b,0,sizeof b);
   for(int c=0;c<N;++c)for(int r=0;r<N;++r) {
    a[c*N+r]=r==c?N+1.+.125*epoch:.03125*(1+.0625*epoch);
    b[r]+=a[c*N+r]*(1.+.125*((c+epoch)%7));
   }
   assert(!feclearexcept(FE_ALL_EXCEPT));
   assert(kls_refactor_solve(s,a,1,b,0,x,0)==KLS_OK);
   unsigned flags=_mm_getcsr();unsigned short x87;__asm__ volatile("fnstsw %0":"=am"(x87));
   for(int c=0;c<N;++c)assert(fabs(x[c]-(1.+.125*((c+epoch)%7)))<1e-11);
   assert(!s->snb_auto_candidate&&s->snb_auto.choice<0);
   assert(s->stats.last_refactor_path!=KLS_REFACTOR_PATH_SNB);
   assert(!atomic_load_explicit(&s->egraph_pool->active_workers,memory_order_acquire));
   assert(!s->componentwise_work&&!s->componentwise_work_context);
   assert(!fetestexcept(FE_INVALID|FE_OVERFLOW|FE_UNDERFLOW|FE_DIVBYZERO));
   assert((_mm_getcsr()&0x1fu)==0);
   if(epoch==13) {
    assert(transaction_capture(s,expected_native)==capacity);
    /* Recompute the exact same current operator through the now-unopposed
       incumbent, explicitly disabling unchanged-input reuse for this test. */
    s->unchanged_refactor_state=-1;s->refactor_input_snapshot_valid=0;
    assert(!feclearexcept(FE_ALL_EXCEPT));
    assert(kls_refactor_solve(s,a,1,b,0,x,0)==KLS_OK);
    unsigned incumbent_flags=_mm_getcsr();unsigned short incumbent_x87;
    __asm__ volatile("fnstsw %0":"=am"(incumbent_x87));
    assert(s->stats.last_refactor_path!=KLS_REFACTOR_PATH_UNCHANGED);
    assert(transaction_capture(s,actual_native)==capacity);
    assert(!memcmp(expected_native,actual_native,capacity*sizeof(double)));
    assert(flags==incumbent_flags&&(x87&63u)==(incumbent_x87&63u));
    for(int c=0;c<N;++c)assert(fabs(x[c]-(1.+.125*((c+epoch)%7)))<1e-11);
   }
  }
  assert(atomic_load_explicit(&transaction_declines,memory_order_relaxed)==1);
  assert(atomic_load_explicit(&transaction_native_changed,memory_order_relaxed));
  free(expected_native);free(actual_native);
  kls_destroy(s);puts("PASS joined speculative decline, partial native writes, same-call recomputation and fresh next epoch");return;
 }
 /* Cache refusal preserves an ordinary plan byte for byte. */
 kls_snb_state *cache_control=kls_snb_build(s,NULL,0,SIZE_MAX);assert(cache_control);
 size_t bytes_before=cache_control->allocation_bytes;
 kls_snb_block *cb=cache_control->blocks;
 int64_t *scatter_before=malloc((size_t)s->nnz*sizeof(int64_t));assert(scatter_before);
 for(int64_t i=0;i<cache_control->col_blk_ptr[s->n];++i)scatter_before[i]=cache_control->blk[i].dst;
 assert(kls_snb_try_row_layout(s,cache_control,0,0));
 assert(!cache_control->row_layout&&!cache_control->transactional&&cache_control->allocation_bytes==bytes_before);
 for(int64_t j=0;j<cb->snode_count;++j)assert(!cb->snodes[j].column_panel);
 for(int64_t i=0;i<cache_control->col_blk_ptr[s->n];++i)assert(scatter_before[i]==cache_control->blk[i].dst);
 free(scatter_before);kls_snb_destroy(cache_control);
 /* An independent schedule must survive retirement of the candidate. */
 kls_snb_state *ordinary=kls_snb_build(s,NULL,0,SIZE_MAX);assert(ordinary);
 assert(ordinary!=s->snb_auto_candidate);s->snb=ordinary;
 s->private_caller_depth=1;
 (void)kls_snb_auto_before(s,a,1,kls_now_seconds());
 s->private_caller_depth=0;
 assert(s->snb==ordinary && !s->snb_auto_candidate && s->snb_auto.choice<0);
 kls_destroy(s);puts("PASS AUTO stream balanced certified calls and current RHS");
}

static void row_cache_checks(void) {
 size_t bytes=10;assert(kls_snb_add_resident_bytes(&bytes,3,4)&&bytes==22);
 assert(!kls_snb_add_resident_bytes(&bytes,UINT64_MAX,8)&&bytes==22);
 bytes=SIZE_MAX;assert(kls_snb_add_resident_bytes(&bytes,0,8)&&bytes==SIZE_MAX);
 assert(!kls_snb_add_resident_bytes(&bytes,1,1)&&bytes==SIZE_MAX);

 kls_snb_snode node={0};node.width=16;node.height=17;
 kls_snb_block block={0};block.snodes=&node;block.snode_count=1;
 kls_snb_state st={0};st.blocks=&block;st.nblocks=1;st.allocation_bytes=100;
 const size_t total=100+200+17*16*sizeof(double);
 assert(kls_snb_row_cache_fits(&st,200,total));
 assert(!kls_snb_row_cache_fits(&st,200,total-1));
 assert(!kls_snb_row_cache_fits(&st,SIZE_MAX,total));
 st.allocation_bytes=SIZE_MAX;assert(!kls_snb_row_cache_fits(&st,200,SIZE_MAX));
 st.allocation_bytes=100;node.height=0;assert(!kls_snb_row_cache_fits(&st,200,total));
 puts("PASS combined native and mirrored panel cache capacity boundaries");
}
static void flag_scope_checks(void) {
 fenv_t saved;assert(!feholdexcept(&saved));
 unsigned old_mx=kls_snb_column_mxflags;int old_flags=kls_snb_column_flags;
 kls_snb_column_mxflags=0;kls_snb_column_flags=0;
 assert(!feraiseexcept(FE_UNDERFLOW));_mm_setcsr(_mm_getcsr()|0x30u);
 kls_snb_row_transaction_begin();
 assert(!kls_snb_transaction_bad_flags());
 assert(kls_snb_column_flags&FE_UNDERFLOW);assert((kls_snb_column_mxflags&0x30u)==0x30u);
 assert(!feraiseexcept(FE_UNDERFLOW));assert(kls_snb_transaction_bad_flags()&0x10u);
 kls_snb_column_mxflags=old_mx;kls_snb_column_flags=old_flags;assert(!fesetenv(&saved));
 puts("PASS incumbent underflow preservation and fresh speculative row flag frame");
}
int main(void){row_cache_checks();flag_scope_checks();unsigned original=_mm_getcsr(); for(int ftz=0;ftz<4;++ftz){_mm_setcsr((original&~0x8040u)|((ftz&1?0x8000u:0)|(ftz&2?0x40u:0)));exercise(0,0);exercise(1,0);exercise(0,1);exercise(1,1);exercise(0,2);exercise(1,2);exercise(0,3);exercise(1,3);} _mm_setcsr(original);return 0;}
