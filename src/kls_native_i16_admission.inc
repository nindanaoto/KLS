/* One paid consultation per stable native layout. No speculative answer is
 * the caller answer: each trial is followed by a fresh original native call. */
typedef struct kls_native_i16_call {
  unsigned epoch;
  int outer,observed,changed,sampled,candidate_ok;
  uint64_t serial_before,runs_before,incumbent_serial,incumbent_runs;
  double candidate_seconds,extra,incumbent_start;
} kls_native_i16_call;

static int kls_native_i16_group_busy(const kls_solver *s) {
  const kls_snb_group_cycle *g=s->snb_group;
  return g && !s->snb_group_disabled &&
    !(g->started && g->samples>=KLS_SNB_GROUP_SAMPLES && !g->active && !g->proposal.starts);
}
static int kls_native_i16_capable(const kls_solver *s) {
  return s && s->numeric && s->symbolic && s->symbolic->nblocks==1 &&
    s->n>1 && s->n<=UINT16_MAX && s->private_caller_depth==1 &&
    s->native_i16_combined_depth==1 && s->native_i16_last_native &&
    s->options.backend==KLS_BACKEND_AUTO &&
    s->accuracy_policy==KLS_ACCURACY_COMPONENTWISE_BACKWARD_ERROR &&
    s->orientation==KLS_ORIENTATION_TRANSPOSE && s->unchanged_refactor_state<0 &&
    !kls_kernel_tuning_overridden(s) && !s->options.tuning_profile_path &&
    s->options.expected_refactorizations<=0 && s->options.expected_solves<=0 &&
    s->stats.last_refactor_path!=KLS_REFACTOR_PATH_NONE &&
    s->stats.last_refactor_path!=KLS_REFACTOR_PATH_UNCHANGED &&
    !s->row_perm && !s->user_col_perm && !s->row_scale && !s->col_scale &&
    !s->dense_tail_cols && !s->padded_run_of && !s->numeric_is_predicted &&
    !s->numeric_needs_refinement && !s->tight_tol_refine &&
    !s->pivot_nudge_count && !s->common.kls_perturb_count &&
    !s->in_solve_refinement && !s->private_scaled_recovery &&
    !s->row_refactor_values_ready && !s->row_refactor_values_dirty &&
    !s->lean_probe_arm && !s->row_accept_pending_side && !s->lean_reaudit_pending_side &&
    !s->padded_pending && !s->eg_pair_pending &&
    !s->refactor_l_packed_valid && !s->refactor_u_packed_valid &&
    !kls_authoritative_packed_lu_enabled(s) &&
    !s->observed_parallel.active && s->observed_parallel.choice!=1 &&
    s->observed_parallel.choice!=2 &&
    (s->stats.last_refactor_path!=KLS_REFACTOR_PATH_MAPPED ||
      s->observed_parallel.choice<0 || !kls_observed_parallel_capable(s)) &&
    s->transpose_cycle.choice<=0 && !s->transpose_cycle.active &&
    !kls_native_i16_group_busy(s) &&
    s->private_original_valid && private_original_operator_input_order(s->private_original) &&
    s->private_original->nnz==(size_t)s->nnz && (uintmax_t)s->nnz<=SIZE_MAX/sizeof(double) &&
    kls_native_i16_environment() && fegetround()==FE_TONEAREST;
}

static kls_native_i16_call kls_native_i16_before(kls_solver *s,
    const double *values,int64_t nrhs,const double *b,double *x) {
  kls_native_i16_call call={0};
  call.outer=s->native_i16_combined_depth==1;
  call.epoch=s->native_i16_epoch;call.serial_before=s->native_i16_serial_calls;call.runs_before=s->native_i16_runs;
  if(!call.outer)return call;
  kls_observed_parallel_state *g=&s->native_i16_policy;
  if(g->active || g->choice<0)return call;
  if(nrhs!=1 || !kls_native_i16_capable(s)) {
    if(g->choice>0 || s->native_i16) {
      g->choice=-1;const double start=kls_now_seconds();kls_native_i16_drop(s);
      call.extra=kls_now_seconds()-start;g->spent+=call.extra;
    }
    return call;
  }
  if(g->choice==2)return call;
  call.observed=1;
  call.changed=memcmp(values,s->private_original->values,(size_t)s->nnz*sizeof(double))!=0;
  if(!call.changed || !kls_observed_parallel_start(g))return call;
  const double start=kls_now_seconds(),prior_spent=g->spent;
  extern char **environ;
  for(char **entry=environ;entry && *entry;++entry)
    if(strncmp(*entry,"KLS_",4)==0){g->choice=-1;goto incumbent;}
  g->active=1;
  g->deadline=start+kls_observed_parallel_preparation_budget(g);
  if(!s->native_i16) {
    const uint64_t llc=kls_execution_domain_llc_bytes();
    uint64_t now=0,deadline=0;
    const double left=g->deadline-kls_now_seconds();
    if(left>0 && left<(double)(UINT64_MAX/1000000000u) && kls_native_i16_clock(&now)) {
      const uint64_t span=(uint64_t)(left*1e9);
      if(span && span<=UINT64_MAX-now)deadline=now+span;
    }
    if(deadline && llc && llc/16<=SIZE_MAX)
      s->native_i16=kls_native_i16_build(s,deadline,(size_t)(llc/16));
  }
  g->spent+=kls_now_seconds()-start;
  const double remaining=g->limit-g->spent;
  if(s->native_i16 && call.epoch==s->native_i16_epoch && remaining>=g->cycle) {
    const uint64_t runs=s->native_i16_runs;
    const double trial_start=kls_now_seconds();g->deadline=trial_start+fmin(g->cycle,remaining);
    int status=kls_refactor_solve_impl(s,values,1,b,0,x,0);
    status=private_componentwise_finish(s,0,1,b,x,0,status);
    call.candidate_seconds=kls_now_seconds()-trial_start;call.sampled=1;
    call.candidate_ok=status==KLS_OK && call.epoch==s->native_i16_epoch &&
      s->native_i16 && s->native_i16_runs!=runs &&
      !s->numeric_is_predicted && !s->numeric_needs_refinement &&
      !s->private_scaled_recovery && !s->pivot_nudge_count && !s->common.kls_perturb_count;
  } else if(call.epoch==s->native_i16_epoch)g->choice=-1;
incumbent:
  if(call.epoch==s->native_i16_epoch){g->deadline=0;g->active=2;}
  call.extra=kls_now_seconds()-start;
  if(call.epoch==s->native_i16_epoch)g->spent=prior_spent+call.extra;
  call.incumbent_serial=s->native_i16_serial_calls;call.incumbent_runs=s->native_i16_runs;
  call.incumbent_start=kls_now_seconds();
  return call;
}

static void kls_native_i16_after(kls_solver *s,kls_native_i16_call call,
    int status,int other_experiment) {
  if(!call.outer)return;
  kls_observed_parallel_state *g=&s->native_i16_policy;
  const int stable=call.epoch==s->native_i16_epoch;
  const uint64_t serial=call.incumbent_start?call.incumbent_serial:call.serial_before;
  const uint64_t runs=call.incumbent_start?call.incumbent_runs:call.runs_before;
  s->native_i16_last_native=s->native_i16_serial_calls!=serial || s->native_i16_runs!=runs;
  if(call.incumbent_start) {
    const double incumbent=kls_now_seconds()-call.incumbent_start;
    if(stable) {
      if(call.sampled)
        (void)kls_observed_parallel_record(g,call.candidate_ok && status==KLS_OK &&
          s->native_i16_serial_calls!=call.incumbent_serial && !other_experiment,
          call.candidate_seconds,incumbent);
      else if(g->choice==1)g->choice=-1;
      g->active=0;
      if(g->choice<0) {
        const double start=kls_now_seconds();kls_native_i16_drop(s);
        const double cleanup=kls_now_seconds()-start;g->spent+=cleanup;call.extra+=cleanup;
      }
    }
  } else if(call.observed && stable) {
    const int clean=status==KLS_OK && !other_experiment && kls_native_i16_capable(s);
    kls_observed_parallel_observe(g,call.changed,clean,s->stats.refactor_seconds+s->stats.solve_seconds);
    if(!clean && g->choice==1) {
      g->choice=-1;const double start=kls_now_seconds();kls_native_i16_drop(s);
      const double cleanup=kls_now_seconds()-start;g->spent+=cleanup;call.extra+=cleanup;
    }
  }
  /* Both failed experiments and preflight cleanup remain caller-visible work. */
  s->stats.refactor_seconds+=call.extra;
}
