/* Optional schedule construction must not install, invalidate, or overwrite
   the incumbent, including on a rejected byte budget or malformed partition. */
#define _GNU_SOURCE
#define _POSIX_C_SOURCE 200809L
#include <stdlib.h>
static int allocation_failure=-1,allocation_calls=0,alternate_ordering;
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

/* Exercise the largest short offset and the first geometry requiring the
 * wide map, including row layout, column repair, scaling, and arena guards. */
static int check_scatter_boundaries(void) {
  kls_solver *solver=calloc(1,sizeof(*solver));kls_snb_state *st=calloc(1,sizeof(*st));
  kls_snb_block *block=calloc(1,sizeof(*block));kls_snb_snode *node=calloc(1,sizeof(*node));
  if(!solver||!st||!block||!node)return 70;
  solver->symbolic=calloc(1,sizeof(*solver->symbolic));if(!solver->symbolic)return 70;
  UF_long q[16]={0},cp[2]={0,1},ri[1]={0};double az[1]={6.},rs[1]={2.};
  solver->symbolic->Q=q;solver->col_ptr=cp;solver->row_idx=ri;
  st->nblocks=1;st->blocks=block;block->snode_count=1;block->snodes=node;node->width=16;
  if(sizeof(*node)!=64)return 71; /* column_ld uses existing descriptor padding. */
  const unsigned columns[]={0,7,8,15};
  const int heights[]={65536,65537,4096,4097},widths[]={1,1,16,16};
  for(int shape=0;shape<4;++shape) {
    const int H=heights[shape],width=widths[shape],wide=shape&1;
    node->height=H;node->width=width;node->below=H-width;
    if(kls_snb_short_scatter_fits(st)!=!wide)return 72;
    const size_t count=(size_t)width*H;double *arena=calloc(count+2,8),*repair=calloc(count+2,8);
    if(!arena||!repair)return 70;
    block->panels=arena+1;
    const uint64_t guard=UINT64_C(0x123456789abcdef0);
    memcpy(arena,&guard,8);memcpy(arena+count+1,&guard,8);
    memcpy(repair,&guard,8);memcpy(repair+count+1,&guard,8);
    uint16_t short_at=0;uint32_t wide_at=0;
    st->scatter_dst16=wide?NULL:&short_at;st->scatter_dst=wide?&wide_at:NULL;
    for(int row=0;row<(width==16?2:1);++row)for(unsigned ci=0;ci<(width==16?4u:1u);++ci)for(int scale=0;scale<=2;scale+=2) {
      const size_t c=columns[ci],r=(size_t)H-1;
      const size_t at=row?(c/8)*8*(size_t)H+r*8+c%8:c*(size_t)H+r;
      st->row_layout=row;wide_at=(uint32_t)at;short_at=(uint16_t)at;
      const double expected=scale?3.:6.;
      kls_snb_scatter_input_column(solver,st,block,node,(UF_long)c,az,rs,scale,NULL,NULL);
      if(block->panels[at]!=expected)return 73;
      block->panels[at]=0.;
      if(row) {
        kls_snb_scatter_input_column(solver,st,block,node,(UF_long)c,az,rs,scale,NULL,repair+1);
        if(repair[1+c*(size_t)H+r]!=expected)return 74;
        repair[1+c*(size_t)H+r]=0.;
      }
      if(memcmp(arena,&guard,8)||memcmp(arena+count+1,&guard,8)||
          memcmp(repair,&guard,8)||memcmp(repair+count+1,&guard,8))return 75;
    }
    free(arena);free(repair);
  }
  free(solver->symbolic);free(solver);free(st);free(block);free(node);
  puts("SHORT_SCATTER_BOUNDARIES short_panel=65536 wide_panel=65537 row_panel=65536 row_wide_panel=65552 row_column=1 repair=1 scaled=1 guards=1 descriptor_bytes=64");
  return 0;
}
static int check_compact_image_policy(void) {
#ifdef KLS_SNB_HAVE_ROW_PANEL
  kls_snb_state *st=calloc(1,sizeof(*st));kls_snb_block *block=calloc(1,sizeof(*block));
  kls_snb_snode *nodes=calloc(3,sizeof(*nodes));kls_snb_edge *edges=calloc(4,sizeof(*edges));
  if(!st||!block||!nodes||!edges)return 76;
  st->nblocks=1;st->blocks=block;st->wnarrow=8;st->allocation_bytes=100;
  block->snodes=nodes;block->snode_count=3;block->edges=edges;
  for(int j=0;j<2;++j){nodes[j].width=16;nodes[j].height=32;nodes[j].below=16;}
  nodes[2].width=1;nodes[2].edge_count=4;
  for(int j=0;j<4;++j){edges[j].producer=j?1:0;edges[j].sub_count=16;}
  const size_t bytes=32*16*sizeof(double);
  if(!kls_snb_row_compact_images(st,0,100+bytes-1)||st->allocation_bytes!=100||
      nodes[0].column_panel||nodes[1].column_panel||nodes[0].column_ld||nodes[1].column_ld)return 77;
  allocation_failure=0;
  if(!kls_snb_row_compact_images(st,0,100+bytes)||allocation_failure!=-1||st->allocation_bytes!=100||
      nodes[0].column_panel||nodes[1].column_panel||nodes[0].column_ld||nodes[1].column_ld)return 78;
  if(!kls_snb_row_compact_images(st,0,100+bytes)||st->allocation_bytes!=100+bytes||
      nodes[0].column_panel||!nodes[1].column_panel||nodes[0].column_ld||nodes[1].column_ld!=32)return 79;
  if(!kls_snb_row_compact_images(st,0,100+bytes)||st->allocation_bytes!=100+bytes||nodes[0].column_panel)return 80;
  free(nodes[1].column_panel);free(edges);free(nodes);free(block);free(st);
  puts("COMPACT_IMAGE_POLICY reuse_rank=1 exact_limit=1 one_short=1 allocation_failure=1 repeat=1 ledger_exact=1");
#endif
  return 0;
}
static int check_preferred_storage(void) {
#ifdef KLS_SNB_HAVE_ROW_PANEL
  for(int failure=0;failure<3;++failure) {
    kls_snb_state *st=calloc(1,sizeof(*st));kls_snb_block *block=calloc(1,sizeof(*block));
    kls_snb_snode *node=calloc(1,sizeof(*node));uint32_t *ranges=calloc(17,sizeof(*ranges));
    uint32_t *writeback=calloc(17,sizeof(*writeback));
    if(!st||!block||!node||!ranges||!writeback)return 81;
    st->blocks=block;st->nblocks=1;st->col_blk_ptr=ranges;
    block->snodes=node;block->snode_count=1;block->nk=16;block->wb_col_off=writeback;
    node->width=16;node->height=32;node->below=16;
    const size_t original=sizeof(*st)+sizeof(*block)+sizeof(*node)+34*sizeof(uint32_t);
    st->allocation_bytes=original;
    const size_t metadata=kls_snb_row_preference_bytes(st);
    allocation_failure=failure<2?failure:-1;
    const int enabled=kls_snb_enable_row_layout_mode(NULL,st,0,SIZE_MAX,1);
    allocation_failure=-1;
    if(enabled!=(failure==2)||st->row_layout!=(failure==2))return 82;
    if(!!st->row_column_preferred!=(failure>0)||!!node->column_panel!=(failure==2))return 83;
    const size_t expected=original+(failure>0?metadata:0)+(failure==2?32*16*sizeof(double):0);
    if(st->allocation_bytes!=expected)return 84;
    for(int j=0;j<17;++j)if(ranges[j]||writeback[j])return 85;
    if(failure==2) {
      if(kls_snb_prefers_column(st,0))return 86;
      kls_snb_prefer_column(st,0);kls_snb_prefer_column(st,0);
      if(!kls_snb_prefers_column(st,0)||st->allocation_bytes!=expected)return 87;
      kls_snb_retry_row(st,0);kls_snb_retry_row(st,0);
      if(kls_snb_prefers_column(st,0)||st->allocation_bytes!=expected)return 88;
    }
    kls_snb_destroy(st);
  }
  puts("PREFERRED_STORAGE_ORACLES metadata_failure=1 mirror_failure=1 zeroed=1 ledger_exact=1 destruction=1");
#endif
  return 0;
}
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
  opts.ordering=alternate_ordering?KLS_ORDERING_AMD:KLS_ORDERING_NATURAL;opts.use_btf=0;
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
  if(alternate_ordering) {
    int changed=0;for(UF_long k=0;k<s->n;++k)changed|=s->symbolic->Q[k]!=k;
    if(!changed)return 59;
  }
  if(sizeof(kls_snb_entry)!=16 || sizeof(*incumbent->col_off_ptr)!=4 ||
      sizeof(*incumbent->col_blk_ptr)!=4)return 51;
  if(incumbent->scatter_dst16||incumbent->scatter_dst) {
    if(incumbent->off||incumbent->blk||incumbent->col_off_ptr||incumbent->col_blk_ptr||
        incumbent->nblocks!=1||sizeof(*incumbent->blocks[0].wb_col_off)!=4)return 52;
    for(UF_long at=0;at<s->nnz;++at)
      if(incumbent->scatter_dst&&incumbent->scatter_dst[at]>=(uint64_t)incumbent->blocks[0].panel_doubles)return 53;
    if(incumbent->scatter_dst16) {
      if(sizeof(*incumbent->scatter_dst16)!=2||incumbent->scatter_dst)return 62;
      const kls_snb_block *block=incumbent->blocks;
      for(UF_long c=0;c<s->n;++c) {
        const kls_snb_snode *node=block->snodes+block->col2snode[c];
        const UF_long oldcol=s->symbolic->Q[c];
        for(UF_long at=s->col_ptr[oldcol];at<s->col_ptr[oldcol+1];++at)
          if(incumbent->scatter_dst16[at]>=(uint64_t)node->height*node->width)return 63;
      }
    }
  } else {
    for(UF_long c=0;c<s->n;++c) {
      if(incumbent->col_off_ptr[c]>incumbent->col_off_ptr[c+1] ||
          incumbent->col_blk_ptr[c]>incumbent->col_blk_ptr[c+1])return 52;
    }
    if((uint64_t)incumbent->col_off_ptr[s->n]+incumbent->col_blk_ptr[s->n]!=(uint64_t)s->nnz)return 53;
    for(uint64_t at=0;at<incumbent->col_off_ptr[s->n];++at)
      if(incumbent->off[at].src>=(uint64_t)s->nnz||incumbent->off[at].oldrow>=(uint64_t)s->n)return 54;
    for(uint64_t at=0;at<incumbent->col_blk_ptr[s->n];++at)
      if(incumbent->blk[at].src>=(uint64_t)s->nnz||incumbent->blk[at].oldrow>=(uint64_t)s->n)return 55;
  }
  if(sizeof(UF_long)>4) {
    UF_long saved=s->n;s->n=(UF_long)((uint64_t)UINT32_MAX+1);
    allocation_calls=0;kls_snb_state *bad=kls_snb_build(s,NULL,0,SIZE_MAX);s->n=saved;
    if(bad||allocation_calls||s->snb!=incumbent||s->snb_decision!=decision)return 56;
    saved=s->nnz;s->nnz=(UF_long)((uint64_t)UINT32_MAX+1);
    allocation_calls=0;bad=kls_snb_build(s,NULL,0,SIZE_MAX);s->nnz=saved;
    if(bad||allocation_calls||s->snb!=incumbent||s->snb_decision!=decision)return 57;
  }

  kls_snb_partition p={malloc((size_t)n*sizeof(int64_t)),calloc((size_t)s->symbolic->nblocks,sizeof(int64_t))};
  if(!p.starts||!p.counts)return 4;
  for(UF_long block=0;block<incumbent->nblocks;++block) {
    kls_snb_block *blk=&incumbent->blocks[block];
    /* Distinct producer below rows must map injectively into the target. */
    for(int64_t j=0;j<blk->snode_count;++j) {
      const kls_snb_snode *target=&blk->snodes[j];
      for(int32_t e=0;e<target->edge_count;++e) {
        const kls_snb_edge *edge=&blk->edges[target->edges_off+e];
        const kls_snb_snode *producer=&blk->snodes[edge->producer];
        const int32_t *dest=blk->edge_dsts+edge->dst_off;
        for(int32_t r=0;r<producer->below;++r) if(dest[r]>=0) {
          if(dest[r]>=target->height)return 49;
          for(int32_t prior=0;prior<r;++prior)if(dest[prior]==dest[r])return 50;
        }
      }
    }
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
static int check_general_btf(int transpose,int scale) {
  enum{N=12};int ap[N+1],ai[N*N],nz=0;double a[N*N],b[N],x[N];
  for(int c=0;c<N;++c) {
    ap[c]=nz;
    for(int r=0;r<N;++r)if(r/6==c/6||(r<6&&c>=6&&r==c%6)) {
      ai[nz]=r;a[nz++]=r==c?8.:.03125;
    }
  }
  ap[N]=nz;kls_options o;kls_default_options(&o);o.backend=KLS_BACKEND_SERIAL;
  o.threads=1;o.use_btf=1;o.scale=scale;o.ordering=KLS_ORDERING_NATURAL;
  o.orientation=transpose?KLS_ORIENTATION_TRANSPOSE:KLS_ORIENTATION_NORMAL;
  kls_solver *s=NULL;if(kls_create(&s)||kls_analyze_csc(s,KLS_INDEX_INT32,N,ap,ai,0,&o)||kls_factor(s,a))return 60;
  s->snode_prepared=1;s->snb_declined=0;
  if(!kls_parallel_lu_sort(s)||!kls_snb_prepare(s)||s->symbolic->nblocks<2||
      s->snb->scatter_dst16||s->snb->scatter_dst||!s->snb->off||!s->snb->blk||!s->snb->col_off_ptr||!s->snb->col_blk_ptr)return 61;
  if(!s->snb->col_off_ptr[N])return 62;
  s->snb_decision=1;
  for(int epoch=1;epoch<=3;++epoch) {
    memset(b,0,sizeof b);
    for(int c=0;c<N;++c)for(int p=ap[c];p<ap[c+1];++p) {
      if(ai[p]==c)a[p]+=.125;
      b[ai[p]]+=a[p]*(1.+.125*((c+epoch)%7));
    }
    if(kls_refactor_solve(s,a,1,b,0,x,0))return 63;
    if(s->stats.last_refactor_path!=KLS_REFACTOR_PATH_SNB)return 65;
    for(int r=0;r<N;++r)if(!isfinite(x[r])||fabs(x[r]-(1.+.125*((r+epoch)%7)))>1e-11)return 64;
  }
  kls_destroy(s);return 0;
}
int main(void) {
  setvbuf(stdout,NULL,_IOLBF,0);
  {int boundary=check_scatter_boundaries();if(boundary)return boundary;}
  {int policy=check_compact_image_policy();if(policy)return policy;}
  {int storage=check_preferred_storage();if(storage)return storage;}
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
  alternate_ordering=1;
  for(int transpose=0;transpose<2;++transpose)for(int csr=0;csr<2;++csr) {
    int status=check_build(129,transpose,csr,0,0);
    if(status){fprintf(stderr,"SNB compact AMD CSC/CSR failed at %d\n",status);return status;}
  }
  alternate_ordering=0;
  for(int transpose=0;transpose<2;++transpose)for(int scale=0;scale<=2;scale+=2) {
    int status=check_general_btf(transpose,scale);
    if(status){fprintf(stderr,"SNB general BTF transpose=%d scale=%d failed at %d\n",transpose,scale,status);return status;}
  }
  puts("CSC_SCATTER_ORACLES compact_map=1 amd_permutation=1 csc_csr=1 scaled=1 btf_offdiag_fallback=1 current_rhs=1");
  return 0;
}
