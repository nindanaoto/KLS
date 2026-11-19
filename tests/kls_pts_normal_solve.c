/* The subtree partition carries the normal-orientation solve (what every
   column-ordered caller gets) in a timed trial against whichever serial
   engine is the incumbent.  The pool kernel must agree with that engine, the
   trial must settle after six alternating solves without a wrong answer, a
   partition whose shared top is too wide for the untimed gates must still be
   published when its predicted cost leaves room, and a factor that gets no
   partition must keep answering exactly as before. */
#include "../src/kls.c"
#undef NDEBUG
#include <assert.h>

enum { CHAINS = 32, LINKMAX = 512, BAND = 5, TOPMAX = 640,
       NMAX = CHAINS * LINKMAX + TOPMAX,
       NZMAX = (2 * BAND + 4) * NMAX + TOPMAX * TOPMAX + TOPMAX * BAND };

static int64_t cp[NMAX + 1], ri[NZMAX];
static double av[NZMAX], truth[NMAX], rhs[NMAX], direct[NMAX], routed[NMAX];

static int nz;
static void put(int row, double value) { ri[nz] = row; av[nz++] = value; }

/* Column-ordered matrix: 32 unsymmetric banded chains, each tied both ways
   to one node of a dense top of the given size. */
static int build(int link, int top) {
  const int body = CHAINS * link, n = body + top;
  nz = 0;
  for (int j = 0; j < n; ++j) {
    cp[j] = nz;
    if (j < body) {
      const int c = j / link, k = j % link;
      for (int d = -BAND; d <= BAND; ++d) {
        if (k + d < 0 || k + d >= link) continue;
        put(j + d, d == 0 ? 8.0 + 0.01 * (j % 5)
                 : d < 0 ? -0.25 - 0.001 * ((j + d) % 7)
                         : -0.5 + 0.002 * ((j + 3 * d) % 3));
      }
      if (k + BAND >= link) put(body + c, 0.125);
    } else {
      const int c = j - body;
      if (c < CHAINS)
        for (int d = 1; d <= BAND; ++d) put(c * link + link - d, -0.375 / d);
      for (int t = 0; t < top; ++t)
        put(body + t, t == c ? 0.5 * top + 6.0
                             : 0.03125 * (1 + (t + c) % 3));
    }
  }
  cp[n] = nz;
  assert(nz <= NZMAX);
  for (int i = 0; i < n; ++i) truth[i] = 1.0 + 0.0625 * (i % 11);
  for (int i = 0; i < n; ++i) rhs[i] = 0.0;
  for (int j = 0; j < n; ++j)
    for (int64_t p = cp[j]; p < cp[j + 1]; ++p) rhs[ri[p]] += av[p] * truth[j];
  return n;
}

static void close_to_truth(const double *x, int n) {
  for (int i = 0; i < n; ++i) assert(fabs(x[i] - truth[i]) < 1e-9);
}

static kls_solver *prepare(int n) {
  kls_options o;
  kls_default_options(&o);
  o.threads = 8;
  o.use_btf = 0;
  o.scale = 0;
  o.static_pivoting = 0;
  o.ordering = KLS_ORDERING_NATURAL;
  kls_solver *s = NULL;
  assert(kls_create(&s) == KLS_OK);
  assert(kls_analyze_csc(s, KLS_INDEX_INT64, n, cp, ri, 0, &o) == KLS_OK);
  assert(kls_factor(s, av) == KLS_OK);
  assert(s->orientation == KLS_ORIENTATION_NORMAL);
  assert(s->row_perm == NULL && s->user_col_perm == NULL &&
         s->row_scale == NULL && s->col_scale == NULL);
  assert(s->symbolic->nblocks == 1u);
  /* The plan needs row-sorted factor columns, as the first refactor of a
     repeated-update run leaves them. */
  assert(kls_parallel_lu_sort(s));
  return s;
}

/* One compact solve with the subtree arm named by the caller. */
static int trial_solve(kls_solver *s, const double *in, double *out) {
  s->pts_trial_solve = 1;
  assert(kls_i32_solve(s, in, out, NULL, NULL, NULL, 0) == 1u);
  const int carried = s->pts_trial_solve == 2;
  s->pts_trial_solve = 0;
  return carried;
}

static void exercise(int link, int top, int wide_top) {
  const int n = build(link, top);
  kls_solver *s = prepare(n);
  assert(kls_i32_solve_ready(s));
  assert(s->pts != NULL && s->pts->streams_ordered && s->pts->nthreads >= 2);
  if (wide_top) {
    /* Published for the timed trial only: no untimed solve, no refactor. */
    assert(n < 16384 && s->numeric->lnz > (UF_long)UINT16_MAX);
    assert(!s->pts->solve_ok && !s->pts->refactor_ok);
    assert(5 * s->pts->ntop >= 4 * top);
  }
  /* No pool yet: the plan alone never admits the solve. */
  assert(!kls_pts_normal_solve_candidate(s));
  assert(ensure_egraph_refactor_pool(s, s->pts->nthreads) != NULL);
  assert(kls_pts_normal_solve_candidate(s));

  /* Kernel against the incumbent route on the first factor's values. */
  assert(trial_solve(s, rhs, direct));
  assert(kls_solve(s, 1, rhs, 0, routed, 0) == KLS_OK);
  close_to_truth(direct, n);
  close_to_truth(routed, n);
  assert(s->pts->nsolve_samples[0] + s->pts->nsolve_samples[1] == 0);

  /* Changed values through the public combined call, then again. */
  for (int round = 0; round < 3; ++round) {
    for (int j = 0; j < n; ++j) {
      for (int64_t p = cp[j]; p < cp[j + 1]; ++p) {
        if (ri[p] == j) {
          av[p] += 0.125;
          rhs[j] += 0.125 * truth[j];
        }
      }
    }
    assert(kls_refactor_solve(s, av, 1, rhs, 0, routed, 0) == KLS_OK);
    close_to_truth(routed, n);
    if (!kls_i32_solve_ready(s) || !kls_pts_normal_solve_candidate(s) ||
        s->row_refactor_values_ready || s->row_refactor_values_dirty) {
      continue;
    }
    assert(trial_solve(s, rhs, direct));
    close_to_truth(direct, n);
  }
  assert(kls_publish_row_refactor_values(s));
  assert(kls_i32_solve_ready(s) && kls_pts_normal_solve_candidate(s));

  /* A busy pool does not change the answer. */
  atomic_store_explicit(&s->egraph_pool->active_workers, 1,
                        memory_order_release);
  (void)trial_solve(s, rhs, direct);
  close_to_truth(direct, n);
  atomic_store_explicit(&s->egraph_pool->active_workers, 0,
                        memory_order_release);

  /* The routed trial only runs behind a pool refactor and never while the
     complete-call consultation is comparing engines. */
  memset(s->pts->nsolve_samples, 0, sizeof(s->pts->nsolve_samples));
  s->pts->nsolve_decision = 0;
  kls_set_last_refactor_path(s, KLS_REFACTOR_PATH_EGRAPH);
  s->observed_parallel.choice = 1;
  assert(kls_solve(s, 1, rhs, 0, routed, 0) == KLS_OK);
  assert(s->pts->nsolve_samples[0] + s->pts->nsolve_samples[1] == 0);
  s->observed_parallel.choice = 2;
  kls_set_last_refactor_path(s, KLS_REFACTOR_PATH_MAPPED);
  assert(kls_solve(s, 1, rhs, 0, routed, 0) == KLS_OK);
  assert(s->pts->nsolve_samples[0] + s->pts->nsolve_samples[1] == 0);
  kls_set_last_refactor_path(s, KLS_REFACTOR_PATH_EGRAPH);
  for (int call = 0; call < 6; ++call) {
    assert(s->pts->nsolve_decision == 0);
    assert(kls_solve(s, 1, rhs, 0, routed, 0) == KLS_OK);
    assert(s->pts->nsolve_samples[0] + s->pts->nsolve_samples[1] == call + 1);
    close_to_truth(routed, n);
  }
  assert(s->pts->nsolve_samples[0] == 3 && s->pts->nsolve_samples[1] == 3);
  assert(s->pts->nsolve_decision != 0);
  /* Either verdict keeps answering correctly, in place as well. */
  for (int verdict = -1; verdict <= 1; verdict += 2) {
    s->pts->nsolve_decision = verdict;
    memcpy(routed, rhs, (size_t)n * sizeof(*rhs));
    assert(kls_solve(s, 1, routed, 0, routed, 0) == KLS_OK);
    close_to_truth(routed, n);
  }
  kls_destroy(s);
}

/* A top so wide that no partition is worth a trial: the small vector keeps
   the native index walk it had before the planner was consulted. */
static void no_partition(int link, int top) {
  const int n = build(link, top);
  kls_solver *s = prepare(n);
  assert(n < 16384 && s->numeric->lnz > (UF_long)UINT16_MAX);
  assert(!kls_i32_solve_ready(s));
  assert(s->pts == NULL && s->i32solve_state == -1);
  assert(kls_solve(s, 1, rhs, 0, routed, 0) == KLS_OK);
  close_to_truth(routed, n);
  kls_destroy(s);
}

int main(void) {
  exercise(512, 32, 0);
  exercise(256, 416, 1);
  no_partition(256, 640);
  return 0;
}
