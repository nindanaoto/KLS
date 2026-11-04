/* Complete certified production calls own this consultation. No matrix or RHS
 * is retained as an answer; the ordinary public frame publishes every epoch. */
typedef struct {
  unsigned epoch;
  int observed,changed,arm;
  double start,extra;
  uint64_t jobs,runs,forward,packed;
} kls_snb_auto_call;

static int kls_snb_auto_capable(const kls_solver *s) {
  return s && s->numeric && s->numeric==s->snb_auto_numeric &&
    (!s->snb || s->snb==s->snb_auto_candidate) &&
    s->symbolic && s->symbolic->nblocks==1 &&
    (s->orientation==KLS_ORIENTATION_TRANSPOSE || s->orientation==KLS_ORIENTATION_NORMAL) &&
    s->options.backend==KLS_BACKEND_AUTO && s->options.threads>1 && s->options.threads<=64 && s->n>1 && s->n<=UINT16_MAX &&
    s->accuracy_policy==KLS_ACCURACY_COMPONENTWISE_BACKWARD_ERROR &&
    s->common.scale<=0 && !s->numeric->Rs && !s->row_perm && !s->user_col_perm &&
    !s->row_scale && !s->col_scale && !s->numeric_is_predicted &&
    !s->numeric_needs_refinement && !s->tight_tol_refine && !s->pivot_nudge_count &&
    !s->common.kls_perturb_count && !s->row_refactor_values_ready &&
    !s->row_refactor_values_dirty && !s->row_accept_pending_side &&
    !s->lean_reaudit_pending_side && !s->lean_probe_arm && !s->private_scaled_recovery &&
    !s->transpose_cycle.active && !s->observed_parallel.active && !s->native_i16_policy.active &&
    !s->componentwise_work && !s->refactor_pool &&
    kls_native_i16_environment() && fegetround()==FE_TONEAREST;
}
static void kls_snb_auto_drop(kls_solver *s) {
  const double start=kls_now_seconds();
  kls_snb_state *candidate=s->snb_auto_candidate;
  if(candidate) {
    if(s->snb==candidate){s->snb=NULL;s->snb_decision=0;}
    kls_snb_destroy(candidate);s->snb_auto_candidate=NULL;
    kls_native_i16_drop(s);
  }
  s->snb_recipe=0;s->snb_auto.choice=-1;
  const double extra=kls_now_seconds()-start;
  s->snb_auto.spent+=extra;s->stats.refactor_seconds+=extra;
}
static int kls_snb_auto_build(kls_solver *s) {
  kls_snb_auto_policy *p=&s->snb_auto;
  const double start=kls_now_seconds(),left=p->limit-p->spent;
  if(!(left>0) || !isfinite(left) || !s->snode_prepared)return 0;
  const double deadline=start+left;
  const long double fallback_bytes=s->auto_amf_fallback?s->auto_amf_fallback_common.memusage:0;
  if(fallback_bytes>=s->snb_auto_memory_limit)return 0;
  const size_t limit=s->snb_auto_memory_limit-(size_t)fallback_bytes;
#ifdef KLS_SNB_AUTO_TEST_BUILD
  kls_snb_state *st=KLS_SNB_AUTO_TEST_BUILD(s,deadline,limit);
#else
  kls_snb_state *st=kls_snb_build(s,NULL,deadline,limit);
#endif
  s->snb_auto_candidate=st;
  int ok=st!=NULL && kls_snb_build_private(st,deadline,limit) &&
    kls_snb_try_row_layout(s,st,deadline,limit) &&
    kls_now_seconds()<deadline &&
    ensure_egraph_refactor_pool(s,s->options.threads)!=NULL;
  if(ok && s->orientation==KLS_ORIENTATION_NORMAL) {
    uint64_t now=0,until=0;const double remaining=deadline-kls_now_seconds();
    if(remaining>0 && remaining<(double)(UINT64_MAX/1000000000u) && kls_native_i16_clock(&now)) {
      const uint64_t span=(uint64_t)(remaining*1e9);
      if(span && span<=UINT64_MAX-now)until=now+span;
    }
    kls_native_i16_drop(s);
    if(until && st->allocation_bytes<limit)
      s->native_i16=kls_native_i16_build(s,until,limit-st->allocation_bytes);
    ok=s->native_i16!=NULL;
  }
  p->spent+=kls_now_seconds()-start;
  return ok && p->spent<=p->limit && kls_now_seconds()<deadline;
}
static __attribute__((noinline)) kls_snb_auto_call kls_snb_auto_before(kls_solver *s,const double *values,
    int64_t nrhs,double start) {
  kls_snb_auto_call call={0};
  if(s->snb_auto.choice<=0 || s->private_caller_depth!=1)return call;
  call.epoch=s->snb_auto_epoch;call.start=start;
  if(nrhs!=1 || !kls_snb_auto_capable(s)) {
    const double cleanup=kls_now_seconds();kls_snb_auto_drop(s);
    call.extra=kls_now_seconds()-cleanup;return call;
  }
  if(s->snb_auto.choice!=3 || !s->snb_auto.rechecked) {
    if(!s->private_original_valid || !s->private_original ||
       !private_original_operator_input_order(s->private_original) ||
       s->private_original->nnz!=(size_t)s->nnz || !s->private_original->values) {
      const double cleanup=kls_now_seconds();kls_snb_auto_drop(s);
      call.extra=kls_now_seconds()-cleanup;return call;
    }
    call.observed=1;
    call.changed=memcmp(values,s->private_original->values,(size_t)s->nnz*sizeof(double))!=0;
    if(!call.changed)return call;
  }
  call.arm=kls_snb_auto_arm(&s->snb_auto);
  const double setup=kls_now_seconds(),prior_spent=s->snb_auto.spent;
  if(call.arm) {
    if(!s->snb_auto_candidate && !kls_snb_auto_build(s)) {
      kls_snb_auto_drop(s);call.observed=0;call.extra=kls_now_seconds()-setup;
      s->snb_auto.spent=prior_spent+call.extra;return call;
    }
    s->snb=s->snb_auto_candidate;s->snb_decision=1;s->snb_declined=0;
    s->snb_recipe=0;
    call.jobs=s->snb->stream_jobs;call.runs=s->native_i16_normal_runs;
  } else if(s->snb_auto_candidate && s->snb==s->snb_auto_candidate) {
    s->snb=NULL;s->snb_decision=0;s->snb_recipe=0;
  }
  call.extra=kls_now_seconds()-setup;
  if(s->snb_auto.choice!=3)s->snb_auto.spent=prior_spent+call.extra;
  return call;
}
static __attribute__((noinline)) void kls_snb_auto_after(kls_solver *s,kls_snb_auto_call call,int status) {
  s->stats.refactor_seconds+=call.extra;
  if(call.epoch!=s->snb_auto_epoch || s->snb_auto.choice<=0)return;
  if(status!=KLS_OK || !kls_snb_auto_capable(s)) {kls_snb_auto_drop(s);return;}
  const int clean=!call.arm || (!call.changed &&
      s->stats.last_refactor_path==KLS_REFACTOR_PATH_UNCHANGED) ||
    (s->snb==s->snb_auto_candidate && s->stats.last_refactor_path==KLS_REFACTOR_PATH_SNB &&
     s->snb->stream_jobs>call.jobs && (s->orientation!=KLS_ORIENTATION_NORMAL ||
       s->native_i16_normal_runs>call.runs));
  /* A retained executor that falls back is no longer the measured engine.
     Exact unchanged-input reuse remains valid and is never a new sample. */
  if(!clean){kls_snb_auto_drop(s);return;}
  if(!call.observed || !call.changed)return;
  double seconds=kls_now_seconds()-call.start-call.extra;
#ifdef KLS_SNB_AUTO_TEST_SECONDS
  seconds=KLS_SNB_AUTO_TEST_SECONDS(call.arm,seconds);
#endif
  if(s->snb_auto.choice==3) {
    (void)kls_snb_auto_recheck(&s->snb_auto,seconds);return;
  }
  kls_snb_auto_sample(&s->snb_auto,1,clean,seconds);
  if(s->snb_auto.choice<0)kls_snb_auto_drop(s);
  else if(s->snb_auto.choice==3) {
    /* The final sample is incumbent. Publish only the decision; the next
       changed call executes fresh values in the retained representation. */
    fprintf(stderr,"KLS_AUTO_SNB retained=1 incumbent=%.9g candidate=%.9g spent=%.9g row_layout=%d\n",
      kls_snb_auto_median(s->snb_auto.seconds[0]),
      kls_snb_auto_median(s->snb_auto.seconds[1]),s->snb_auto.spent,
      s->snb_auto_candidate?s->snb_auto_candidate->row_layout:0);
  }
}
