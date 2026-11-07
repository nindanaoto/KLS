/* Symbolic proposal only. Unknown METIS flop estimates are not made positive
   by invention: a missing estimate requires a material predicted-fill saving.
   Actual factors/answers still follow the ordinary native/certificate paths. */
#ifdef __linux__
static uint64_t kls_execution_domain_llc_bytes(void);
#endif
static int kls_small_btf_metis_improves(double old_fill,double old_work,
    double fill,double work) {
  return isfinite(old_fill) && old_fill>0 && old_fill<=UINT16_MAX &&
    isfinite(old_work) && old_work>0 && old_work<100000 &&
    isfinite(fill) && fill>0 && fill<=.80*old_fill &&
    isfinite(work) && (work==-1.0 || (work>0 && work<=old_work));
}
static void kls_auto_small_btf_candidate(kls_pattern_candidate *c,const kls_options *o) {
  if(!c || !o || o->ordering!=KLS_ORDERING_AUTO || o->backend!=KLS_BACKEND_AUTO ||
     o->scale!=KLS_SCALE_AUTO || o->threads<1 || o->threads>64 ||
     (o->tuning_profile_path && !kls_policy_model_only(c->policy)) ||
     c->orientation!=KLS_ORIENTATION_TRANSPOSE || c->selected_ordering!=KLS_ORDERING_AMD ||
     c->small_btf_fallback.symbolic || c->block_order_perm || !c->symbolic ||
     c->n<=1 || c->n>4096 || !c->nnz || c->nnz>UINT16_MAX ||
     !c->symbolic->do_btf || c->symbolic->nblocks<=1 || c->symbolic->maxblock<=1 ||
     c->symbolic->structural_rank!=c->n ||
     kls_prestatic_ordering_ctx || kls_matched_low_pair_ordering_ctx ||
     kls_value_matched_ordering_ctx || kls_block_order_proposal ||
     !isfinite(c->score) || c->score<=0 || c->score>UINT16_MAX ||
     !isfinite(c->symbolic->est_flops) || c->symbolic->est_flops<=0 ||
     c->symbolic->est_flops>=100000)return;
  /* Includes simultaneous compact maps, both symbolic owners, native index
     mirrors, and the bounded input deck. No numerical result is retained. */
  const uint64_t bytes=512+256*(uint64_t)c->n+24*(uint64_t)c->nnz+
    16*(uint64_t)ceil(c->score);
#ifdef __linux__
  const uint64_t llc=kls_execution_domain_llc_bytes();
#else
  const uint64_t llc=0;
#endif
  if(!llc || bytes>llc/32)return;
  kls_options trial_options=*o;trial_options.ordering=KLS_ORDERING_METIS;
  trial_options.use_btf=1;
  trilinos_klu_l_symbolic *trial=NULL;trilinos_klu_l_common common;
  kls_separator_analysis separator={0};const double start=kls_now_seconds();
  const int status=analyze_with_ordering(c->n,c->col_ptr,c->row_idx,&trial_options,
    KLS_ORDERING_METIS,&trial,&common,&separator);
  /* Both symbolics own their internal permutations. External coordinate
     rewrites are excluded; public matrix coordinates must stay identical. */
  const int block_rewrite=kls_block_order_proposal!=NULL;
  free(kls_block_order_proposal);kls_block_order_proposal=NULL;
  const double score=status==KLS_OK && trial?symbolic_score(trial):DBL_MAX;
  int accepted=status==KLS_OK && trial && !block_rewrite &&
    trial->do_btf &&
    /* Given-order KLU leaves rank unknown. The identical input pattern has
       already been proved structurally full rank by the retained AMD BTF. */
    (trial->structural_rank==c->n || trial->structural_rank==KLS_KLU_EMPTY) &&
    trial->nblocks==c->symbolic->nblocks && trial->maxblock==c->symbolic->maxblock &&
    kls_small_btf_metis_improves(c->score,c->symbolic->est_flops,score,trial->est_flops);

  UF_long *pattern=NULL;
#ifndef KLS_SMALL_BTF_PATTERN_MALLOC
#define KLS_SMALL_BTF_PATTERN_MALLOC malloc
#endif
  if(accepted) {
    pattern=KLS_SMALL_BTF_PATTERN_MALLOC(((size_t)c->n+1+(size_t)c->nnz)*sizeof(*pattern));
    if(!pattern)accepted=0;
  }
#ifdef KLS_SMALL_BTF_AUTO_TEST_EVENT
  KLS_SMALL_BTF_AUTO_TEST_EVENT(c,trial,score,accepted);
#endif
  if(accepted) {
    memcpy(pattern,c->col_ptr,((size_t)c->n+1)*sizeof(*pattern));
    memcpy(pattern+c->n+1,c->row_idx,(size_t)c->nnz*sizeof(*pattern));
    kls_small_btf_fallback *f=&c->small_btf_fallback;
    f->owner=trial;f->symbolic=c->symbolic;f->common=c->common;
    f->ordering=c->selected_ordering;f->orientation=c->orientation;
    f->col_ptr=c->col_ptr;f->row_idx=c->row_idx;f->n=c->n;f->nnz=c->nnz;
    f->pattern=pattern;
    kls_separator_analysis_move(&f->separator,&c->separator);
    c->symbolic=trial;trial=NULL;c->common=common;c->score=score;
    c->selected_ordering=KLS_ORDERING_METIS;c->generic_amf3_span_variant_selected=0;
    kls_separator_analysis_move(&c->separator,&separator);
  }
  if(trial)trilinos_klu_l_free_symbolic(&trial,&common);
  kls_separator_analysis_clear(&separator);
  (void)start; /* Public analyze timing includes the entire proposal. */
}

/* Consume a lease only after a real failed public numeric operation. */
static int kls_small_btf_restore(kls_solver *s) {
  if(!kls_small_btf_owner_matches(s))return 0;
  kls_small_btf_fallback *f=&s->small_btf_fallback;
  if(!f->pattern || memcmp(f->pattern,s->col_ptr,((size_t)s->n+1)*sizeof(UF_long)) ||
     memcmp(f->pattern+s->n+1,s->row_idx,(size_t)s->nnz*sizeof(UF_long))) {
    kls_small_btf_fallback_clear(f);return 0;
  }
  kls_small_btf_fallback saved=*f;memset(f,0,sizeof(*f));
  kls_metis_race_abandon(s);free_numeric(s);
  trilinos_klu_l_free_symbolic(&s->symbolic,&s->common);
  kls_separator_analysis_clear(&s->separator);
  s->symbolic=saved.symbolic;s->common=saved.common;
  kls_separator_analysis_move(&s->separator,&saved.separator);
  s->stats.selected_ordering=saved.ordering;
  s->stats.nblocks=(int64_t)s->symbolic->nblocks;
  s->stats.max_block=(int64_t)s->symbolic->maxblock;
  s->stats.structural_rank=(int64_t)s->symbolic->structural_rank;
  s->stats.estimated_flops=s->symbolic->est_flops;
  s->generic_nd_portfolio_selected=0;s->generic_amf3_span_variant_selected=0;
  kls_invalidate_factor_etree_stats(s);free(saved.pattern);
#ifdef KLS_SMALL_BTF_RESTORE_TEST_EVENT
  KLS_SMALL_BTF_RESTORE_TEST_EVENT(s);
#endif
  return 1;
}
static int kls_small_btf_numeric_failure(int status) {
  return status==KLS_ERR_FACTOR_FAILED || status==KLS_ERR_REFACTOR_FAILED ||
    status==KLS_ERR_SINGULAR || status==KLS_ERR_OUT_OF_MEMORY;
}
