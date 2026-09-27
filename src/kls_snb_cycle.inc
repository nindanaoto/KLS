/* Optional single-RHS grouping trial. An internal sample always runs the
   incumbent last; the caller never receives the speculative answer. */
static int kls_snb_group_capable(const kls_solver *s) {
  return s->snb_cost.enabled && !s->snb_group_disabled &&
    s->options.threads == 1 && s->snb != NULL && s->snb_decision > 0 &&
    s->stats.last_refactor_path == KLS_REFACTOR_PATH_SNB &&
    s->unchanged_refactor_state < 0 &&
    !s->factor_preps_deferred && !s->prestatic_deferred &&
    !s->rowmatch_deferred && !s->metis_race_deferred &&
    !s->auto_scale_deferred && !s->tight_pivot_deferred &&
    !s->numeric_is_predicted && !s->numeric_needs_refinement &&
    !s->row_accept_pending_side && !s->lean_reaudit_pending_side &&
    (s->row_reaudit_state==0 || s->row_reaudit_state==7) &&
    (s->lean_reaudit_state==0 || s->lean_reaudit_state==5) &&
    (s->row_publish_experiment==0 || s->row_publish_experiment==3) &&
    !s->transpose_cycle.active &&
    s->row_perm == NULL && s->user_col_perm == NULL &&
    s->row_scale == NULL && s->col_scale == NULL;
}

static double kls_snb_group_median(const double *values) {
  double v[KLS_SNB_GROUP_SAMPLES];
  memcpy(v,values,sizeof(v));
  for(int i=1;i<KLS_SNB_GROUP_SAMPLES;++i) {
    double x=v[i];int j=i;
    while(j && v[j-1]>x){v[j]=v[j-1];--j;}v[j]=x;
  }
  return (v[KLS_SNB_GROUP_SAMPLES/2-1]+v[KLS_SNB_GROUP_SAMPLES/2])*0.5;
}

static void kls_snb_group_reject(kls_solver *s) {
  kls_snb_group_cycle *g=s->snb_group;
  s->snb_group_disabled=1;
  if(g == NULL)return;
  kls_snb_destroy(g->other);g->other=NULL;
  free(g->proposal.starts);free(g->proposal.counts);
  g->proposal=(kls_snb_partition){0};
}

static double kls_snb_group_cost(const kls_solver *s) {
  const kls_snb_cost_model *m=&s->snb_cost;
  double total=0.0;
  for(UF_long b=0;b<s->snb->nblocks;++b) {
    const kls_snb_block *blk=&s->snb->blocks[b];
    for(int64_t k=0;k<blk->snode_count;++k) {
      const kls_snb_snode *sn=&blk->snodes[k];
      const int w=sn->width,ph=sn->below;
      double updates=0.0;
      for(int i=0;i<w;++i)updates+=(double)(w-i-1)*(w+ph-i-1);
      total+=m->panel[0]+m->panel[1]*((double)w*ph+(double)w*(w-1)/2)+
        m->panel[2]*updates+m->panel[3]*(double)w*sn->height;
      for(int e=0;e<sn->edge_count;++e) {
        const kls_snb_edge *edge=&blk->edges[sn->edges_off+e];
        total+=kls_snb_model_edge(m,w,edge->sub_count,blk->snodes[edge->producer].below);
      }
    }
  }
  return total;
}

/* The allowance does not depend on the proposed partition's predicted win.
   Hints estimate total work; absent hints, only already-observed work funds
   exploration. Snapshot the allowance when search starts, rather than
   spending now against a budget that might become available later. */
static double kls_snb_group_allowance(const kls_solver *s,double cycle) {
  const double hinted=s->options.expected_refactorizations>0?
    (double)s->options.expected_refactorizations*cycle:0.0;
  const double work=fmax(s->snb_group->incumbent_total,hinted);
  return isfinite(work)&&cycle>0.0&&isfinite(cycle)?
    KLS_SNB_GROUP_MAX_OVERHEAD_FRACTION*work:0.0;
}

static double kls_snb_group_envelope(const kls_solver *s,double cycle) {
  const double search=KLS_SNB_GROUP_SAMPLES*cycle;
  const double samples=KLS_SNB_GROUP_SAMPLES*cycle;
  /* Match the existing builder's cold-preparation reserve: incumbent
     preparation plus the search envelope, in addition to search itself. */
  const double construction=s->snb->preparation_seconds+search;
  return search+construction+samples;
}

static int kls_snb_group_propose(kls_solver *s,double cycle) {
  kls_snb_group_cycle *g=s->snb_group;
  const double start=kls_now_seconds();
  const double deadline=start+KLS_SNB_GROUP_SAMPLES*cycle;
  const double base=kls_snb_group_cost(s);
  if(!isfinite(base)||base<=0.0)return 0;
  /* Bound extra proposal/reverse-incidence storage using already-owned
     schedule and native value storage, not a matrix-size exception. */
  const long double workspace=(long double)s->n*
    (12*sizeof(int)+sizeof(uint64_t)+sizeof(int64_t))+
    (long double)s->snb->nblocks*sizeof(int64_t)+sizeof(int)+
    (long double)s->numeric->unz*sizeof(int);
  const long double allowed=(long double)s->snb->allocation_bytes+
    ((long double)s->numeric->lnz+s->numeric->unz)*sizeof(double)+
    (long double)s->n*(2*sizeof(double)+sizeof(long double));
  if(workspace>allowed)return 0;
  g->proposal.starts=malloc((size_t)s->n*sizeof(int64_t));
  g->proposal.counts=calloc((size_t)s->snb->nblocks,sizeof(int64_t));
  if(!g->proposal.starts||!g->proposal.counts)return 0;
  for(UF_long b=0;b<s->snb->nblocks;++b) {
    const kls_snb_block *blk=&s->snb->blocks[b];
    g->proposal.counts[b]=blk->snode_count;
    for(int64_t k=0;k<blk->snode_count;++k)
      g->proposal.starts[blk->k1+k]=blk->snodes[k].start;
  }
  double saving=0.0;
  for(UF_long b=0;b<s->snb->nblocks && kls_now_seconds()<deadline;++b) {
    const kls_snb_block *blk=&s->snb->blocks[b];
    if(blk->nk<2 || !blk->snode_count)continue;
    const UF_long k=blk->k1;
    double before=0.0,after=0.0;
    g->proposal.counts[b]=kls_snb_model_regroup(s->numeric->LUbx[b],
      s->numeric->Lip+k,s->numeric->Llen+k,s->numeric->Uip+k,s->numeric->Ulen+k,
      blk->nk,g->proposal.starts+k,blk->snode_count,&s->snb_cost,deadline,&before,&after);
    saving+=before-after;
  }
  g->overhead+=kls_now_seconds()-start;
  g->predicted_fraction=saving/base;
  return isfinite(g->predicted_fraction) && g->predicted_fraction>0.0 &&
    g->predicted_fraction<1.0;
}

/* The admission guard excludes external row/column transforms. Numeric Rs
   scaling and either internal orientation are supported: the stored CSC
   pattern and input permutation still describe the original public matrix. */
static int kls_snb_group_residual(const kls_solver *s,const double *values,
    const double *b,const double *x,long double *r) {
  /* input_to_csc maps public positions TO internal positions, not vice
     versa. SNB has just prepared this unscaled CSC values array. */
  const double *a_csc=s->input_to_csc?s->values:values;
  if(!a_csc)return 0;
  long double rr=0.0L,bb=0.0L;
  for(UF_long i=0;i<s->n;++i){r[i]=b[i];bb+=(long double)b[i]*b[i];}
  for(UF_long col=0;col<s->n;++col)for(UF_long p=s->col_ptr[col];p<s->col_ptr[col+1];++p) {
    const UF_long row=s->row_idx[p];
    const double a=a_csc[p];
    if(s->orientation==KLS_ORIENTATION_TRANSPOSE)r[col]-=(long double)a*x[row];
    else r[row]-=(long double)a*x[col];
  }
  for(UF_long i=0;i<s->n;++i)rr+=r[i]*r[i];
  return isfinite(rr) && isfinite(bb) && rr<=1e-16L*(bb>0.0L?bb:1.0L);
}

static int private_refactor_solve_body(kls_solver *s,const double *values,int64_t nrhs,
    const double *b,int64_t ldb,double *x,int64_t ldx) {
  /* Unsupported call shapes use the unchanged API, including its validation. */
  if(!s||!values||!b||!x||nrhs!=1||(ldb && ldb<(int64_t)s->n)||
     (ldx && ldx<(int64_t)s->n)||!s->snb_cost.enabled||
     s->snb_group_disabled||s->options.threads!=1||
     (s->snb_group && s->snb_group->active))
    return kls_refactor_solve_impl(s,values,nrhs,b,ldb,x,ldx);

  /* Ineligible engines must not pay even the sampling-clock overhead.
     A later ordinary call may qualify SNB; begin observing on that next
     call. If an existing experiment loses eligibility, cancel it once. */
  if(!kls_snb_group_capable(s)) {
    if(s->snb_group)kls_snb_group_reject(s);
    return kls_refactor_solve_impl(s,values,nrhs,b,ldb,x,ldx);
  }
  const double call_start=kls_now_seconds();
  const unsigned epoch=s->snb_group_epoch;
  if(kls_snb_group_capable(s) && s->snb_group==NULL) {
    s->snb_group=calloc(1,sizeof(*s->snb_group));
    if(s->snb_group==NULL)s->snb_group_disabled=1;
  }
  kls_snb_group_cycle *g=s->snb_group;
  const double prior_overhead=g?g->overhead:0.0;
  double extra=0.0,candidate_seconds=0.0;
  int sampled=0,candidate_ok=0;
  if(g && kls_snb_group_capable(s) && g->cycles>=KLS_SNB_GROUP_SAMPLES) {
    const double cycle=kls_snb_group_median(g->history);
    const double horizon=fmax((double)g->cycles,
      (double)s->options.expected_refactorizations-g->cycles);
    const double allowance=kls_snb_group_allowance(s,cycle);
    const double trial_reserve=KLS_SNB_GROUP_SAMPLES*cycle;
    if(!g->started && allowance>kls_snb_group_envelope(s,cycle) &&
        (g->cycles>=KLS_REPEATED_UPDATE_MIN_REFACTORIZATIONS ||
        horizon>=KLS_REPEATED_UPDATE_MIN_REFACTORIZATIONS)) {
      g->started=1;
      g->experiment_limit=allowance;
      if(!kls_snb_group_propose(s,cycle))kls_snb_group_reject(s);
    }
    if(!s->snb_group_disabled && g->proposal.starts && !g->other) {
      const double budget=fmin(g->experiment_limit-g->overhead-trial_reserve,
        horizon*g->predicted_fraction*cycle/
        KLS_COMPACT_CYCLE_PREPARATION_PAYBACK-g->overhead-trial_reserve);
      /* Starting at the first incumbent-build-time crossover makes the sole
         attempt vulnerable to small cold-build variations. Reserve another
         already-measured preparation envelope, rather than fitting a time
         or matrix-size margin to one workload. This still waits entirely
         within the same prospective payback allowance. */
      if(budget>s->snb->preparation_seconds+g->overhead) {
        const double build_start=kls_now_seconds();
        g->other=kls_snb_build(s,&g->proposal,build_start+budget,s->snb->allocation_bytes);
        g->overhead+=kls_now_seconds()-build_start;
        free(g->proposal.starts);free(g->proposal.counts);
        g->proposal=(kls_snb_partition){0};
        if(!g->other)kls_snb_group_reject(s);
      }
    }
    if(g->other && g->samples<KLS_SNB_GROUP_SAMPLES &&
        g->overhead+(KLS_SNB_GROUP_SAMPLES-g->samples)*cycle<g->experiment_limit) {
      double *trial_x=malloc((size_t)s->n*sizeof(double));
      long double *residual=malloc((size_t)s->n*sizeof(long double));
      if(!trial_x||!residual){free(trial_x);free(residual);kls_snb_group_reject(s);}
      else {
        kls_snb_state *incumbent=s->snb;
        s->snb=g->other;g->other=incumbent;g->active=1;
        const double saved_adaptive=s->adaptive_refactor_seconds;
        const double sample_start=kls_now_seconds();
        int status=kls_refactor_solve_impl(s,values,1,b,ldb,trial_x,0);
        candidate_seconds=kls_now_seconds()-sample_start;
        double trial_seconds=candidate_seconds;
        for(int i=0;i<g->samples;++i)trial_seconds+=g->challenger[i];
        candidate_ok=status==KLS_OK && epoch==s->snb_group_epoch &&
          trial_seconds<=KLS_SNB_GROUP_SAMPLES*cycle &&
          s->stats.last_refactor_path==KLS_REFACTOR_PATH_SNB &&
          kls_snb_group_residual(s,values,b,trial_x,residual);
        if(epoch==s->snb_group_epoch) {
          g->other=s->snb;s->snb=incumbent;s->snb_decision=1;
        }
        g->active=0;s->adaptive_refactor_seconds=saved_adaptive;
        g->overhead+=kls_now_seconds()-sample_start;
        free(trial_x);free(residual);sampled=1;
        if(!candidate_ok)kls_snb_group_reject(s);
      }
    }
  }
  extra=kls_now_seconds()-call_start;
  if(g && g->started)g->overhead=prior_overhead+extra;
  if(g && g->started && g->overhead>g->experiment_limit)
    kls_snb_group_reject(s);
  /* Recompute the complete incumbent cycle, including cache invalidation and
     normal recovery. Speculation must never leave a candidate answer live. */
  const double incumbent_start=kls_now_seconds();
  int status;
  if(epoch!=s->snb_group_epoch && s->numeric==NULL && s->symbolic!=NULL) {
    /* Failed speculative recovery may have discarded the native factor as
       well as both schedules. Rebuild through the normal factor/solve path;
       refactor alone would incorrectly turn this into INVALID_ARGUMENT. */
    status=kls_factor(s,values);
    s->stats.refactor_seconds=kls_now_seconds()-incumbent_start;
    s->stats.solve_seconds=0.0;
    if(status==KLS_OK)status=kls_solve(s,nrhs,b,ldb,x,ldx);
  } else {
    status=kls_refactor_solve_impl(s,values,nrhs,b,ldb,x,ldx);
  }
  const double incumbent_seconds=kls_now_seconds()-incumbent_start;
  const double bookkeeping_start=kls_now_seconds();
  if(g && status==KLS_OK && epoch==s->snb_group_epoch && kls_snb_group_capable(s)) {
    g->history[g->cycles%KLS_SNB_GROUP_SAMPLES]=incumbent_seconds;
    g->incumbent_total+=incumbent_seconds;
    if(g->cycles<UINT_MAX)++g->cycles;
    if(sampled && candidate_ok && !s->snb_group_disabled) {
      g->incumbent[g->samples]=incumbent_seconds;
      g->challenger[g->samples++]=candidate_seconds;
      if(g->samples==KLS_SNB_GROUP_SAMPLES) {
        const double base=kls_snb_group_median(g->incumbent);
        const double saving=base-kls_snb_group_median(g->challenger);
        if(!(saving>=KLS_COMPACT_CYCLE_MIN_SECONDS_SAVED &&
             saving>KLS_COMPACT_CYCLE_MIN_RELATIVE_SAVING*base))kls_snb_group_reject(s);
      }
    }
    if(!s->snb_group_disabled && g->other && g->samples==KLS_SNB_GROUP_SAMPLES) {
      const double saving=kls_snb_group_median(g->incumbent)-kls_snb_group_median(g->challenger);
      const double horizon=fmax((double)g->cycles,
        (double)s->options.expected_refactorizations-g->cycles);
      if(horizon*saving>KLS_COMPACT_CYCLE_PREPARATION_PAYBACK*
          (g->overhead+kls_now_seconds()-bookkeeping_start)) {
        kls_snb_destroy(s->snb);s->snb=g->other;g->other=NULL;
        s->snb_group_disabled=1; /* Accepted: no further experiment this analysis. */
      }
    }
  } else if(g && g->started)kls_snb_group_reject(s);
  if(g && g->started)g->overhead+=kls_now_seconds()-bookkeeping_start;
  s->stats.refactor_seconds+=kls_now_seconds()-call_start-incumbent_seconds;
  return status;
}
