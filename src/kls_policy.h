#ifndef KLS_POLICY_H
#define KLS_POLICY_H

#include <float.h>
#include <math.h>
#include <stdatomic.h>
#include <stdint.h>

#include "kls_tuning_profile.h"
#include "kls_tuning.inc"

/* Private, per-analysis accounting. Workers share reservations, never solver
 * numerics or observations from another solver. No matrix identity is an input. */
typedef struct kls_policy_context {
  atomic_flag lock;
  kls_policy_cost_model cost;
  double limit, reserved, spent;
  unsigned admitted, declined;
  int funded, kernel_overrides;
} kls_policy_context;

static inline void kls_policy_init(kls_policy_context *p,
    const kls_policy_cost_model *cost) {
  atomic_flag_clear(&p->lock);
  p->cost = cost && cost->enabled ? *cost : (kls_policy_cost_model){
    .seconds = {KLS_POLICY_SCAN_SECONDS, KLS_POLICY_PAIR_SECONDS,
      KLS_POLICY_SORT_SECONDS, KLS_POLICY_FLOP_SECONDS,
      KLS_POLICY_ENTRY_SECONDS, KLS_POLICY_SOLVE_ENTRY_SECONDS,
      KLS_POLICY_COPY_ENTRY_SECONDS, KLS_POLICY_DISPATCH_SECONDS,
      KLS_POLICY_FLOP_SECONDS, KLS_POLICY_ENTRY_SECONDS,
      KLS_POLICY_DISPATCH_SECONDS, 0.0}};
  p->limit=p->reserved=p->spent=0;
  p->admitted=p->declined=0;p->funded=p->kernel_overrides=0;
}

/* Cost-only profiles preserve every default kernel threshold. Legacy and
 * explicit kernel profiles retain their original conservative admission. */
static inline int kls_policy_model_only(const kls_policy_context *p) {
  return p && p->cost.enabled && !p->kernel_overrides;
}

static inline void kls_policy_lock(kls_policy_context *p) {
  while(atomic_flag_test_and_set_explicit(&p->lock,memory_order_acquire)) {}
}
static inline void kls_policy_unlock(kls_policy_context *p) {
  atomic_flag_clear_explicit(&p->lock,memory_order_release);
}

static inline double kls_policy_numeric_seconds(const kls_policy_context *p,
    double input, double fill, double flops, double blocks,
    double refactors, double solves) {
  if(!p || !isfinite(input) || input<0 || !isfinite(fill) || fill<=0 ||
     fill>=DBL_MAX/4 || !isfinite(flops) || flops<0 || !isfinite(blocks) || blocks<0 ||
     !isfinite(refactors) || refactors<0 || !isfinite(solves) || solves<0)
    return NAN;
  const double *c=p->cost.seconds;
  const double steady=c[KLS_POLICY_FLOP]*flops+c[KLS_POLICY_ENTRY]*fill+
    c[KLS_POLICY_DISPATCH]*blocks;
  const double factor=c[KLS_POLICY_FACTOR_FLOP]*flops+
    c[KLS_POLICY_FACTOR_ENTRY]*fill+c[KLS_POLICY_FACTOR_DISPATCH]*blocks+
    c[KLS_POLICY_COPY_ENTRY]*input;
  /* Include value preparation and original-operator certification. These
   * remain costs even when a sparse numeric contains almost no arithmetic. */
  const double refactor=steady+c[KLS_POLICY_COPY_ENTRY]*input;
  const double solve=c[KLS_POLICY_SOLVE_ENTRY]*fill+
    c[KLS_POLICY_SCAN]*input+c[KLS_POLICY_DISPATCH]*blocks;
  /* First-changed preparation is paid once, never multiplied by the whole
   * horizon. Initial factor allocation is not a steady refactor cost. */
  const double first=refactors>0?c[KLS_POLICY_FIRST_REFACTOR_ENTRY]*fill:0;
  const double total=factor+first+refactors*refactor+solves*solve;
  return isfinite(total)&&total>0?total:NAN;
}

static inline double kls_policy_graph_seconds(const kls_policy_context *p,
    double vertices,double entries,double pairs) {
  if(!p || !isfinite(vertices) || vertices<0 || !isfinite(entries) || entries<0 ||
     !isfinite(pairs) || pairs<0)return NAN;
  const double total=p->cost.seconds[KLS_POLICY_SCAN]*(vertices+entries)+
    p->cost.seconds[KLS_POLICY_PAIR]*pairs;
  return isfinite(total)&&total>0?total:NAN;
}

static inline void kls_policy_fund(kls_policy_context *p,double lifecycle) {
  if(!p || !isfinite(lifecycle) || lifecycle<=0)return;
  kls_policy_lock(p);
  if(!p->funded) {
    p->limit=KLS_POLICY_EXPLORATION_FRACTION*lifecycle;p->funded=1;
  }
  kls_policy_unlock(p);
}

static inline int kls_policy_reserve(kls_policy_context *p,double estimate) {
  if(!p)return 1; /* Private arithmetic/reference callers have no lifecycle. */
  kls_policy_lock(p);
  const int ok=p->funded && isfinite(estimate) && estimate>0 &&
    estimate<=p->limit-p->spent-p->reserved;
  if(ok){p->reserved+=estimate;++p->admitted;}else ++p->declined;
  kls_policy_unlock(p);return ok;
}

static inline void kls_policy_charge(kls_policy_context *p,double reserve,
    double seconds) {
  if(!p)return;
  kls_policy_lock(p);
  if(isfinite(reserve)&&reserve>0)p->reserved=fmax(0,p->reserved-reserve);
  /* Invalid clock evidence cannot create credit or authorize more trials. */
  p->spent=isfinite(seconds)&&seconds>=0 && isfinite(p->spent+seconds)
    ?p->spent+seconds:INFINITY;
  kls_policy_unlock(p);
}

#endif
