/* Complete-call consultation. Speculative success is never a caller answer. */
static int kls_observed_parallel_capable(const kls_solver *s) {
  return s && s->symbolic && s->numeric &&
    s->accuracy_policy==KLS_ACCURACY_COMPONENTWISE_BACKWARD_ERROR &&
    s->options.backend==KLS_BACKEND_AUTO && s->options.threads>1 &&
    s->options.expected_refactorizations<=0 && s->options.expected_solves<=0 &&
    !s->stats.tuning_profile_active &&
    s->n>=(UF_long)s->tuning.egraph_min_size &&
    s->common.flops>0 && s->common.flops<
      s->tuning.egraph_min_flops_per_thread*s->options.threads &&
    s->common.scale<=0 && s->numeric->Rs==NULL &&
    s->row_perm==NULL && s->user_col_perm==NULL &&
    s->row_scale==NULL && s->col_scale==NULL &&
    !s->dense_tail_cols && !s->numeric_is_predicted &&
    !s->numeric_needs_refinement && !s->tight_tol_refine &&
    !s->pivot_nudge_count && !s->common.kls_perturb_count &&
    !s->in_solve_refinement && !s->private_scaled_recovery &&
    !s->factor_preps_deferred && !s->prestatic_deferred &&
    !s->rowmatch_deferred && !s->metis_race_deferred &&
    !s->auto_scale_deferred && !s->tight_pivot_deferred &&
    !s->pts_build_deferred && !s->generic_nd_portfolio_selected &&
    s->unchanged_refactor_state<0 &&
    !s->snb && s->lean_choice<0 && s->lean_probe_arm<=0 &&
    !s->row_refactor_values_ready && !s->row_refactor_values_dirty &&
    !s->row_accept_pending_side && !s->lean_reaudit_pending_side &&
    !s->refactor_l_packed_valid && !s->refactor_u_packed_valid &&
    s->solve_contract_probe>=1 && s->solve_contract_probe<=2 &&
    s->transpose_cycle.choice<=0 && !s->transpose_cycle.active &&
    s->private_original_valid &&
    private_original_operator_input_order(s->private_original) &&
    !getenv("KLS_EGRAPH_REFACTOR_FLOOR") &&
    !getenv("KLS_DISABLE_EGRAPH_REFACTOR") &&
    !getenv("KLS_DISABLE_OBSERVED_PARALLEL") &&
    !kls_row_refactor_env_enabled();
}

typedef struct {
  unsigned epoch;
  int observed,changed,sampled,candidate_ok;
  double candidate_seconds,extra,incumbent_start;
} kls_observed_parallel_call;

static void kls_observed_parallel_restore_native(kls_solver *s) {
  kls_sync_authoritative_packed_l_values(s);
  kls_sync_authoritative_packed_u_values(s);
  s->refactor_l_packed_valid=0;s->refactor_u_packed_valid=0;
}

__attribute__((noinline))
static kls_observed_parallel_call kls_observed_parallel_before(kls_solver *s,
    const double *values,int64_t nrhs,const double *b,double *x) {
  kls_observed_parallel_call call={0};
  if(s && (s->native_i16_policy.active || s->snb_auto.choice>0))return call;
  kls_observed_parallel_state *g=&s->observed_parallel;
  if(nrhs!=1 || !kls_observed_parallel_capable(s)) {
    if(g->choice>0 || g->retry_pending) {
      g->choice=-1;g->retry_pending=0;kls_observed_parallel_restore_native(s);
      destroy_egraph_refactor_pool(s);free_refactor_schedule(s);
    }
    return call;
  }
  if(g->choice==2 || g->choice<0 ||
      s->stats.last_refactor_path!=KLS_REFACTOR_PATH_MAPPED) return call;
  call.epoch=s->observed_parallel_epoch;call.observed=1;
  /* The snapshot is used solely to classify observed workload. It never
     proves factor validity, a solve, nonsingularity or accuracy. Most changed
     inputs differ early; exact repeats contribute no exploration funding. */
  call.changed=s->private_original->nnz==(size_t)s->nnz &&
    s->private_original->values && memcmp(values,s->private_original->values,
                                        (size_t)s->nnz*sizeof(double))!=0;
  if(!call.changed || !kls_observed_parallel_start(g)) return call;
  /* A process-level engine override is an explicit policy. Collection may
     precede the funding threshold, but no competing consultation starts. */
  extern char **environ;
  for(char **entry=environ;entry && *entry;++entry)
    if(strncmp(*entry,"KLS_",4)==0) {g->choice=-1;return call;}
  if(fegetround()!=FE_TONEAREST) {g->choice=-1;return call;}
  const double start=kls_now_seconds();
  const double prior_spent=g->spent;
  const double budget=kls_observed_parallel_preparation_budget(g);
  const int retained_retry=g->retry_pending && s->egraph_pool && s->refactor_level_ptr;
  g->active=1;g->deadline=start+(retained_retry?fmax(budget,g->cycle):budget);
  /* The complete plan and joined pool survived a deadline, not a layout
     change. The actual refactor still revalidates every map/cache under its
     numeric deadline; only duplicate preparation is avoided here. */
  int ready=(retained_retry || (budget>0 && kls_build_refactor_schedule(s))) &&
    kls_egraph_refactor_is_eligible(s);
  g->spent+=kls_now_seconds()-start;
  const double remaining=g->limit-g->spent;
  if(ready && remaining>=g->cycle &&
      call.epoch==s->observed_parallel_epoch) {
    const double trial_start=kls_now_seconds();
    g->deadline=trial_start+fmin(g->cycle,remaining);
    int status=kls_refactor_solve_impl(s,values,1,b,0,x,0);
    status=private_componentwise_finish(s,0,1,b,x,0,status);
    call.candidate_seconds=kls_now_seconds()-trial_start;
    call.sampled=1;
    call.candidate_ok=status==KLS_OK &&
      call.epoch==s->observed_parallel_epoch &&
      s->stats.last_refactor_path==KLS_REFACTOR_PATH_EGRAPH &&
      !s->numeric_is_predicted && !s->numeric_needs_refinement &&
      !s->private_scaled_recovery && !s->pivot_nudge_count &&
      !s->common.kls_perturb_count;
    const int deadline=s->egraph_pool &&
      atomic_load_explicit(&s->egraph_pool->shared.trial_deadline_hit,memory_order_acquire) &&
      !atomic_load_explicit(&s->egraph_pool->active_workers,memory_order_acquire);
    kls_observed_parallel_restore_native(s);
    if(status==KLS_OK && !call.candidate_ok && deadline && s->refactor_level_ptr &&
       call.epoch==s->observed_parallel_epoch && kls_observed_parallel_capable(s))
      (void)kls_observed_parallel_defer_deadline(g);
  } else g->choice=-1;
  g->deadline=0;g->active=2;
  call.extra=kls_now_seconds()-start;
  /* Full extra work, including setup, numeric, snapshot, solve and proof,
     is charged even on failed or timed-out trials. */
  if(call.epoch==s->observed_parallel_epoch) {
    g->spent=prior_spent+call.extra;
  }
  call.incumbent_start=kls_now_seconds();
  return call;
}

__attribute__((noinline))
static void kls_observed_parallel_after(kls_solver *s,
    kls_observed_parallel_call call,int status) {
  kls_observed_parallel_state *g=&s->observed_parallel;
  const int stable=call.epoch==s->observed_parallel_epoch;
  if(call.incumbent_start) {
    const double incumbent=kls_now_seconds()-call.incumbent_start;
    if(stable && call.sampled)
      (void)kls_observed_parallel_record(g,call.candidate_ok && status==KLS_OK &&
          s->stats.last_refactor_path==KLS_REFACTOR_PATH_MAPPED,
          call.candidate_seconds,incumbent);
    else if(g->choice==1) g->choice=-1;
    if(!stable || status!=KLS_OK ||
       (g->retry_pending && s->stats.last_refactor_path!=KLS_REFACTOR_PATH_MAPPED))g->choice=-1;
    g->active=0;
    s->stats.refactor_seconds+=call.extra;
    if(g->choice<0) {
      g->retry_pending=0;
      const double cleanup_start=kls_now_seconds();
      destroy_egraph_refactor_pool(s);
      free_refactor_schedule(s);
      const double cleanup_seconds=kls_now_seconds()-cleanup_start;
      g->spent+=cleanup_seconds;s->stats.refactor_seconds+=cleanup_seconds;
    }
  }
  if(call.observed && stable && !call.incumbent_start) {
    const int clean=status==KLS_OK && kls_observed_parallel_capable(s) &&
      s->stats.last_refactor_path==KLS_REFACTOR_PATH_MAPPED;
    kls_observed_parallel_observe(g,call.changed,clean,
        s->stats.refactor_seconds+s->stats.solve_seconds);
    if(!clean) kls_observed_parallel_invalidate(s);
  }
}

/* Keep optional consultation and its recovery out of the public entry point.
   Inlining speculative numeric work here needlessly enlarges every call's
   instruction footprint, including small serial solves. */
__attribute__((noinline))
static int kls_observed_parallel_refactor_solve(kls_solver *s,
    const double *values,int64_t nrhs,double *b,double *x) {
  const kls_observed_parallel_call call=
    kls_observed_parallel_before(s,values,nrhs,b,x);
  int status;
  if(call.incumbent_start && !s->numeric && s->symbolic) {
    status=kls_factor(s,values);
    if(status==KLS_OK) status=kls_solve(s,nrhs,b,0,x,0);
  } else status=private_refactor_solve_body(s,values,nrhs,b,0,x,0);
  status=private_caller_finish(s,0,nrhs,b,x,0,status);
  kls_observed_parallel_after(s,call,status);
  return status;
}
