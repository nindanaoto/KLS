/* Optional schedule construction must not install, invalidate, or overwrite
   the incumbent, including on a rejected byte budget or malformed partition. */
#define _GNU_SOURCE
#define _POSIX_C_SOURCE 200809L
#include <stdlib.h>
static int allocation_failure=-1,allocation_calls=0;
static int fail_allocation(void) {
  ++allocation_calls;
  if(allocation_failure<0)return 0;
  if(allocation_failure--==0){allocation_failure=-1;return 1;}
  return 0;
}
static void *checked_malloc(size_t n) {return fail_allocation()?NULL:malloc(n);}
static void *checked_calloc(size_t n,size_t size) {return fail_allocation()?NULL:calloc(n,size);}
#define malloc checked_malloc
#define calloc checked_calloc
#include "../src/kls.c"
#undef malloc
#undef calloc

static int check_build(int n,int transpose,int csr,int alias,int reject) {
  int *ptr=calloc((size_t)n+1,sizeof(int));
  int *idx=malloc((size_t)n*n*sizeof(int));
  double *a=malloc((size_t)n*n*sizeof(double));
  double *b=calloc((size_t)n,sizeof(double)),*x=malloc((size_t)n*sizeof(double));
  if(!ptr||!idx||!a||!b||!x)return 1;
  int nz=0;
  for(int c=0;c<n;++c) {
    ptr[c]=nz;
    for(int r=0;r<n;++r) if(abs(r-c)<=8 || (71*c+31*r)%97==0) {
      idx[nz]=r;a[nz]=r==c?4.0:0.001*((r+c)%7-3);b[csr?c:r]+=a[nz++];
    }
  }
  ptr[n]=nz;
  kls_options opts;kls_default_options(&opts);
  opts.threads=1;opts.backend=KLS_BACKEND_SERIAL;
  opts.ordering=KLS_ORDERING_NATURAL;opts.use_btf=0;
  opts.orientation=transpose?KLS_ORIENTATION_TRANSPOSE:KLS_ORIENTATION_NORMAL;
  kls_solver *s=NULL;
  if(kls_create(&s)||(csr?kls_analyze_csr:kls_analyze_csc)(s,KLS_INDEX_INT32,n,ptr,idx,0,&opts)||
     kls_factor(s,a))return 2;
  if(s->snb==NULL) {
    s->snb_cost.enabled=1;
    if(kls_refactor_solve(s,a,1,b,0,x,0)||s->snb_group||s->snb_group_disabled)return 29;
    memset(&s->snb_cost,0,sizeof(s->snb_cost));
  }
  s->snode_prepared=1;s->snb_declined=0;
  if(!kls_parallel_lu_sort(s)||!kls_snb_prepare(s))return 3;
  kls_snb_state *incumbent=s->snb;
  const int decision=s->snb_decision;
  kls_snb_partition p={malloc((size_t)n*sizeof(int64_t)),calloc((size_t)s->symbolic->nblocks,sizeof(int64_t))};
  if(!p.starts||!p.counts)return 4;
  for(UF_long block=0;block<incumbent->nblocks;++block) {
    kls_snb_block *blk=&incumbent->blocks[block];
    p.counts[block]=blk->snode_count;
    for(int64_t j=0;j<blk->snode_count;++j)p.starts[blk->k1+j]=blk->snodes[j].start;
  }
  allocation_calls=0;
  kls_snb_state *copy=kls_snb_build(s,&p,0.0,incumbent->allocation_bytes);
  const int build_allocations=allocation_calls;
  if(!copy || copy->allocation_bytes!=incumbent->allocation_bytes ||
     s->snb!=incumbent || s->snb_decision!=decision)return 5;
  kls_snb_destroy(copy);
  for(int failure=0;failure<build_allocations;++failure) {
    allocation_failure=failure;
    copy=kls_snb_build(s,&p,0.0,SIZE_MAX);
    allocation_failure=-1;
    /* Optional parallel scheduling metadata may fall back to the serial
       schedule; mandatory allocation failures reject the candidate. Neither
       outcome may mutate the incumbent or leave partial panels live. */
    if(s->snb!=incumbent||s->snb_decision!=decision)return 25;
    if(copy) {
      if(copy->allocation_bytes>incumbent->allocation_bytes)return 26;
      for(UF_long block=0;block<copy->nblocks;++block) {
        const kls_snb_block *blk=&copy->blocks[block];
        if(blk->snode_count && (!blk->panels||!blk->writeback||!blk->col2snode))return 27;
      }
      kls_snb_destroy(copy);
    }
  }
  /* Compare incremental proposal pricing against an independently built
     schedule, including consumers outside each changed range. */
  s->snb_cost.enabled=1;
  for(int i=0;i<3;++i)s->snb_cost.narrow[i]=1e-9;
  for(int i=0;i<6;++i)s->snb_cost.dense[i]=1e-9;
  for(int i=0;i<4;++i)s->snb_cost.panel[i]=1e-9;
  const double base_price=kls_snb_group_cost(s);
  double delta=0.0;
  for(UF_long block=0;block<incumbent->nblocks;++block) {
    const kls_snb_block *blk=&incumbent->blocks[block];
    double before=0.0,after=0.0;const UF_long k=blk->k1;
    kls_snb_group_model unions={0};
    unions.n=(int)blk->nk;unions.lu=s->numeric->LUbx[block];
    unions.lip=s->numeric->Lip+k;unions.llen=s->numeric->Llen+k;
    unions.mark=calloc((size_t)blk->nk,sizeof(int));
    unions.pmark=calloc((size_t)blk->nk,sizeof(int));
    if(!unions.mark||!unions.pmark)return 45;
    /* Every legal width/offset, including overlapping ranges and stamp wrap,
       must match independent per-range unions exactly. */
    for(int c=0;c<unions.n;++c)for(int w=1;w<=32&&c+w<=unions.n;++w) {
      int left[33],right[33];
      if(c==0&&w==32)unions.stamp=INT_MAX-1;
      if(!kls_snb_model_split_below(&unions,c,w,left,right))return 46;
      for(int split=1;split<w;++split)
        if(left[split]!=kls_snb_model_below(&unions,c,split)||
           right[split]!=kls_snb_model_below(&unions,c+split,w-split))return 47;
    }
    unions.deadline=kls_now_seconds()-1.0;
    int expired_left[33],expired_right[33];
    if(kls_snb_model_split_below(&unions,0,2,expired_left,expired_right)||!unions.expired)return 48;
    free(unions.mark);free(unions.pmark);
    p.counts[block]=kls_snb_model_regroup(s->numeric->LUbx[block],
      s->numeric->Lip+k,s->numeric->Llen+k,s->numeric->Uip+k,s->numeric->Ulen+k,
      blk->nk,p.starts+k,blk->snode_count,&s->snb_cost,0.0,&before,&after);
    delta+=after-before;
  }
  copy=kls_snb_build(s,&p,0.0,SIZE_MAX);
  if(!copy)return 23;
  s->snb=copy;
  double candidate_price=kls_snb_group_cost(s);
  s->snb=incumbent;kls_snb_destroy(copy);
  if(fabs(candidate_price-base_price-delta)>1e-10*base_price)return 24;
  memset(&s->snb_cost,0,sizeof(s->snb_cost));
  for(UF_long block=0;block<incumbent->nblocks;++block) {
    const kls_snb_block *blk=&incumbent->blocks[block];
    p.counts[block]=blk->snode_count;
    for(int64_t j=0;j<blk->snode_count;++j)p.starts[blk->k1+j]=blk->snodes[j].start;
  }
  if(kls_snb_build(s,&p,0.0,incumbent->allocation_bytes-1)!=NULL ||
     kls_snb_build(s,&p,kls_now_seconds()-1.0,SIZE_MAX)!=NULL)return 6;
  p.starts[0]=1;
  if(kls_snb_build(s,&p,0.0,SIZE_MAX)!=NULL || s->snb!=incumbent ||
     s->snb_decision!=decision || s->snb_declined)return 7;
  s->snb_decision=1;
  for(int i=0;i<3;++i) {
    a[0]+=0.001;b[0]+=0.001;
    if(kls_refactor_solve(s,a,1,b,0,x,0))return 8;
    for(int r=0;r<n;++r)if(!isfinite(x[r])||fabs(x[r]-1.0)>1e-10)return 9;
  }
  if(s->snb_group!=NULL)return 10; /* No profile: no adaptive state. */
  if(n==129) {
    p.starts[0]=0;
    s->snb_cost.enabled=1;
    for(int i=0;i<3;++i)s->snb_cost.narrow[i]=1e-9;
    for(int i=0;i<6;++i)s->snb_cost.dense[i]=1e-9;
    for(int i=0;i<4;++i)s->snb_cost.panel[i]=1e-9;
    s->snb_group=calloc(1,sizeof(*s->snb_group));
    if(!s->snb_group)return 11;
    kls_snb_group_cycle *g=s->snb_group;
    g->started=1;g->cycles=16;g->samples=KLS_SNB_GROUP_SAMPLES-1;
    g->incumbent_total=1000.0;g->experiment_limit=20.0;
    g->other=kls_snb_build(s,&p,0.0,SIZE_MAX);
    if(!g->other)return 12;
    kls_snb_state *candidate=g->other;
    /* Deterministic verdict fixture; numerical cycles themselves remain real. */
    for(int i=0;i<KLS_SNB_GROUP_SAMPLES;++i) {
      g->history[i]=1.0;g->incumbent[i]=reject?1e-9:1.0;
      g->challenger[i]=reject?1.0:1e-9;
    }
    s->unchanged_refactor_state=-1;s->direct_klu_choice=-1;
    s->lean_choice=-1;s->lean_probe_arm=0;s->snb_decision=1;
    s->factor_preps_deferred=0;s->numeric_is_predicted=0;
    kls_set_last_refactor_path(s,KLS_REFACTOR_PATH_SNB);
    if(!kls_snb_group_capable(s))return 13;
    const int old_row_audit=s->row_reaudit_state,old_lean_audit=s->lean_reaudit_state;
    const int old_publication=s->row_publish_experiment;
    s->row_accept_pending_side=1;
    if(kls_snb_group_capable(s))return 30;
    s->row_accept_pending_side=0;s->row_reaudit_state=4;
    if(kls_snb_group_capable(s))return 31;
    s->row_reaudit_state=old_row_audit;s->lean_reaudit_state=3;
    if(kls_snb_group_capable(s))return 32;
    s->lean_reaudit_state=old_lean_audit;s->row_publish_experiment=2;
    if(kls_snb_group_capable(s))return 33;
    s->row_publish_experiment=old_publication;
    /* Real short reuse sequence: an explicit horizon must prevent paying
       for an unaffordable experiment, across formats/orientations/aliasing. */
    kls_snb_group_cycle saved_cycle=*g;
    const int64_t saved_hint=s->options.expected_refactorizations;
    memset(g,0,sizeof(*g));s->options.expected_refactorizations=99;
    for(int step=0;step<99;++step) {
      a[0]+=0.00001;b[0]+=0.00001;
      if(alias)memcpy(x,b,(size_t)n*sizeof(double));
      if(kls_refactor_solve(s,a,1,alias?x:b,0,x,0))return 45;
      for(int r=0;r<n;++r)if(!isfinite(x[r])||fabs(x[r]-1.0)>1e-10)return 46;
      if(g->started||g->other||g->proposal.starts)return 47;
    }
    *g=saved_cycle;s->options.expected_refactorizations=saved_hint;
    a[0]+=0.001;b[0]+=0.001;
    memcpy(x,b,(size_t)n*sizeof(double));
    if(reject==2) {
      /* Singular input forces speculative numeric recovery. Its failure
         must not surface as an invalid-argument error or retain a dangling
         challenger. The analysis remains owned until normal destruction. */
      memset(a,0,(size_t)nz*sizeof(double));
      int status=kls_refactor_solve(s,a,1,alias?x:b,0,x,0);
      if(status==KLS_OK || status==KLS_ERR_INVALID_ARGUMENT || g->active || g->other)return 28;
      goto done;
    }
    if(kls_refactor_solve(s,a,1,alias?x:b,0,x,0))return 14;
    if(s->snb!=(reject?incumbent:candidate)||!s->snb_group_disabled||g->other)return 15;
    for(int r=0;r<n;++r)if(!isfinite(x[r])||fabs(x[r]-1.0)>1e-10)return 16;
    /* Ordinary APIs must consume the same native factor after promotion. */
    if(kls_solve(s,1,b,0,x,0))return 17;
    for(int r=0;r<n;++r)if(fabs(x[r]-1.0)>1e-10)return 18;
    /* An invalidating free owns both schedules, including the detached one. */
    s->snb_group_disabled=0;
    g->other=kls_snb_build(s,&p,0.0,SIZE_MAX);
    if(!g->other)return 19;
    const unsigned epoch=s->snb_group_epoch;
    kls_snb_free(s);
    if(s->snb||g->other||!s->snb_group_disabled||s->snb_group_epoch==epoch)return 20;
    if((csr?kls_analyze_csr:kls_analyze_csc)(s,KLS_INDEX_INT32,n,ptr,idx,0,&opts)||
       s->snb_group||s->snb_group_disabled||kls_factor(s,a)||kls_solve(s,1,b,0,x,0))return 21;
    for(int r=0;r<n;++r)if(fabs(x[r]-1.0)>1e-10)return 22;
  }
done:
  free(p.starts);free(p.counts);kls_destroy(s);
  free(ptr);free(idx);free(a);free(b);free(x);
  return 0;
}
int main(void) {
  /* Budget tests use no matrix identity, synthetic model, or CPU coefficient.
     A short hinted workload cannot fund the complete search/build/trial envelope;
     an unhinted workload earns its allowance from completed real work. */
  kls_solver budget_solver={0};kls_snb_group_cycle budget_cycle={0};
  kls_snb_state budget_state={0};budget_state.preparation_seconds=0.001;
  budget_solver.snb=&budget_state;
  budget_solver.snb_group=&budget_cycle;
  budget_solver.options.expected_refactorizations=99;
  double budget=kls_snb_group_allowance(&budget_solver,0.002);
  if(fabs(budget-0.00396)>1e-12 ||
     budget>=kls_snb_group_envelope(&budget_solver,0.002))return 40;
  budget_solver.options.expected_refactorizations=0;
  if(kls_snb_group_allowance(&budget_solver,0.002)!=0.0)return 41;
  budget_cycle.incumbent_total=2.0;
  if(fabs(kls_snb_group_allowance(&budget_solver,0.002)-0.04)>1e-12)return 42;
  if(kls_snb_group_allowance(&budget_solver,0.002)<=
     kls_snb_group_envelope(&budget_solver,0.002))return 44;
  if(kls_snb_group_allowance(&budget_solver,NAN)!=0.0)return 43;
  const int sizes[]={17,32,33,64,129};
  for(unsigned i=0;i<sizeof(sizes)/sizeof(sizes[0]);++i) {
    int status=check_build(sizes[i],0,0,0,0);
    if(status){fprintf(stderr,"SNB build n=%d failed at %d\n",sizes[i],status);return status;}
  }
  for(int transpose=0;transpose<2;++transpose)for(int csr=0;csr<2;++csr)
    for(int alias=0;alias<2;++alias)for(int reject=0;reject<3;++reject) {
      int status=check_build(129,transpose,csr,alias,reject);
      if(status){fprintf(stderr,"SNB cycle transpose=%d csr=%d alias=%d reject=%d failed at %d\n",
        transpose,csr,alias,reject,status);return status;}
    }
  return 0;
}
