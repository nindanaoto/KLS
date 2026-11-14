/* The subtree partition carries the transpose solve (the orientation every
   row-ordered caller uses): the pool kernel must agree with the incumbent
   solve on single- and multi-block factors, and the timed trial must settle
   after six alternating solves without ever returning a wrong answer. */
#include "../src/kls.c"
#undef NDEBUG
#include <assert.h>

enum { CHAINS = 32, LINK = 512, BAND = 5, TOP = CHAINS,
       MAIN = CHAINS * LINK + TOP, EXTRA = 6, NMAX = MAIN + EXTRA,
       NZMAX = (2 * BAND + 4) * NMAX + TOP * TOP + TOP * BAND };

static int64_t rp[NMAX + 1], ci[NZMAX];
static double av[NZMAX], truth[NMAX], rhs[NMAX], direct[NMAX], routed[NMAX];

static int nz;
static void put(int col, double value) { ci[nz] = col; av[nz++] = value; }

/* Row-ordered matrix: 32 unsymmetric banded chains (wide enough that the
   factor outgrows the 16-bit solve cache), each tied both
   ways to one node of a dense top; optional extra unknowns coupled one way
   only, so BTF splits them into their own blocks with off-diagonal entries. */
static int build(int extra) {
  const int n = MAIN + (extra ? EXTRA : 0);
  nz = 0;
  for (int i = 0; i < n; ++i) {
    rp[i] = nz;
    if (i < CHAINS * LINK) {
      const int c = i / LINK, k = i % LINK;
      for (int d = -BAND; d <= BAND; ++d) {
        if (k + d < 0 || k + d >= LINK) continue;
        put(i + d, d == 0 ? 8.0 + 0.01 * (i % 5)
                 : d < 0 ? -0.25 - 0.001 * ((i + d) % 7)
                         : -0.5 + 0.002 * ((i + 3 * d) % 3));
      }
      if (k + BAND >= LINK) put(CHAINS * LINK + c, 0.125);
      if (extra && k == 3) put(MAIN + 3 + c % 3, 0.0625);
    } else if (i < MAIN) {
      const int c = i - CHAINS * LINK;
      for (int d = 1; d <= BAND; ++d) put(c * LINK + LINK - d, -0.375 / d);
      for (int t = 0; t < TOP; ++t)
        put(CHAINS * LINK + t, t == c ? 6.0 : 0.03125 * (1 + (t + c) % 3));
    } else {
      const int e = i - MAIN;
      if (e < 3) put(e * LINK + 7, 0.1875);
      put(i, 3.0 + 0.5 * e);
    }
  }
  rp[n] = nz;
  assert(nz <= NZMAX);
  for (int i = 0; i < n; ++i) truth[i] = 1.0 + 0.0625 * (i % 11);
  for (int i = 0; i < n; ++i) {
    rhs[i] = 0.0;
    for (int64_t p = rp[i]; p < rp[i + 1]; ++p) rhs[i] += av[p] * truth[ci[p]];
  }
  return n;
}

static void compare(const kls_solver *s, int n) {
  (void)s;
  for (int i = 0; i < n; ++i) {
    assert(fabs(direct[i] - truth[i]) < 1e-10);
    assert(fabs(routed[i] - truth[i]) < 1e-10);
    assert(fabs(direct[i] - routed[i]) <= 1e-12 * fabs(routed[i]));
  }
}

static void exercise(int extra) {
  const int n = build(extra);
  kls_options o;
  kls_default_options(&o);
  o.threads = 8;
  o.use_btf = extra;
  o.scale = 0;
  o.static_pivoting = 0;
  o.ordering = KLS_ORDERING_NATURAL;
  kls_solver *s = NULL;
  assert(kls_create(&s) == KLS_OK);
  assert(kls_analyze_csr(s, KLS_INDEX_INT64, n, rp, ci, 0, &o) == KLS_OK);
  assert(kls_factor(s, av) == KLS_OK);
  assert(s->orientation == KLS_ORIENTATION_TRANSPOSE);
  assert(s->row_perm == NULL && s->user_col_perm == NULL &&
         s->row_scale == NULL && s->col_scale == NULL);
  assert(extra ? s->symbolic->nblocks > 1u : s->symbolic->nblocks == 1u);

  /* The plan needs row-sorted factor columns, as the first refactor of a
     repeated-update run leaves them. */
  assert(kls_parallel_lu_sort(s));
  /* No pool yet: the plan alone never admits the solve. */
  assert(kls_i32_solve_ready(s));
  assert(s->pts != NULL && s->pts->solve_ok && s->pts->nthreads >= 2);
  assert(!kls_pts_transpose_solve_candidate(s));
  assert(ensure_egraph_refactor_pool(s, s->pts->nthreads) != NULL);
  assert(kls_pts_transpose_solve_candidate(s));

  /* Kernel against the incumbent route on the first factor's values. */
  memcpy(direct, rhs, (size_t)n * sizeof(*rhs));
  assert(kls_pts_transpose_solve(s, direct) == 1u);
  assert(kls_solve(s, 1, rhs, 0, routed, 0) == KLS_OK);
  compare(s, n);
  assert(s->pts->tsolve_samples[0] + s->pts->tsolve_samples[1] == 0);

  /* Changed values through the public combined call, then again. */
  for (int round = 0; round < 3; ++round) {
    for (int i = 0; i < n; ++i) {
      for (int64_t p = rp[i]; p < rp[i + 1]; ++p) {
        if (ci[p] == i) {
          av[p] += 0.125;
          rhs[i] += 0.125 * truth[i];
        }
      }
    }
    assert(kls_refactor_solve(s, av, 1, rhs, 0, routed, 0) == KLS_OK);
    assert(kls_i32_solve_ready(s) && kls_pts_transpose_solve_candidate(s));
    memcpy(direct, rhs, (size_t)n * sizeof(*rhs));
    assert(kls_pts_transpose_solve(s, direct) == 1u);
    compare(s, n);
  }

  /* A busy pool declines without touching the caller's vector. */
  atomic_store_explicit(&s->egraph_pool->active_workers, 1,
                        memory_order_release);
  memcpy(direct, rhs, (size_t)n * sizeof(*rhs));
  assert(kls_pts_transpose_solve(s, direct) == 0u);
  assert(memcmp(direct, rhs, (size_t)n * sizeof(*rhs)) == 0);
  atomic_store_explicit(&s->egraph_pool->active_workers, 0,
                        memory_order_release);

  /* The routed trial only runs behind a pool refactor and never while the
     complete-call consultation is comparing engines. */
  kls_set_last_refactor_path(s, KLS_REFACTOR_PATH_EGRAPH);
  s->observed_parallel.choice = 1;
  assert(kls_solve(s, 1, rhs, 0, routed, 0) == KLS_OK);
  assert(s->pts->tsolve_samples[0] + s->pts->tsolve_samples[1] == 0);
  s->observed_parallel.choice = 2;
  kls_set_last_refactor_path(s, KLS_REFACTOR_PATH_MAPPED);
  assert(kls_solve(s, 1, rhs, 0, routed, 0) == KLS_OK);
  assert(s->pts->tsolve_samples[0] + s->pts->tsolve_samples[1] == 0);
  kls_set_last_refactor_path(s, KLS_REFACTOR_PATH_EGRAPH);
  for (int call = 0; call < 6; ++call) {
    assert(s->pts->tsolve_decision == 0);
    assert(kls_solve(s, 1, rhs, 0, routed, 0) == KLS_OK);
    assert(s->pts->tsolve_samples[0] + s->pts->tsolve_samples[1] == call + 1);
    for (int i = 0; i < n; ++i) assert(fabs(routed[i] - truth[i]) < 1e-10);
  }
  assert(s->pts->tsolve_samples[0] == 3 && s->pts->tsolve_samples[1] == 3);
  assert(s->pts->tsolve_decision != 0);
  /* Either verdict keeps answering correctly, in place as well. */
  for (int verdict = -1; verdict <= 1; verdict += 2) {
    s->pts->tsolve_decision = verdict;
    memcpy(routed, rhs, (size_t)n * sizeof(*rhs));
    assert(kls_solve(s, 1, routed, 0, routed, 0) == KLS_OK);
    for (int i = 0; i < n; ++i) assert(fabs(routed[i] - truth[i]) < 1e-10);
  }
  kls_destroy(s);
}

int main(void) {
  exercise(0);
  exercise(1);
  return 0;
}
