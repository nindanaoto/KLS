/* Exercise a real numerical dispatch immediately after thread creation.
 * Private admission floors keep this small fixture eligible without changing
 * production policy. Run with TSan to cover the startup publication boundary.
 */
#include "../src/kls.c"
#undef NDEBUG
#include <assert.h>

static void exercise(int crew) {
  enum { N = 128, NZ = 256 };
  int64_t ptr[N + 1], rows[NZ];
  double a[NZ], truth[N], rhs[N], transpose_rhs[N], answer[N];
  for (int j = 0; j < N; ++j) {
    ptr[j] = 2 * j;
    rows[2 * j] = j & ~1;
    rows[2 * j + 1] = (j & ~1) + 1;
    a[2 * j] = j % 2 ? 0.125 : 4.0;
    a[2 * j + 1] = j % 2 ? 4.0 : 0.25;
    truth[j] = 1.0 + 0.25 * (j % 7 - 3);
  }
  ptr[N] = NZ;
  for (int j = 0; j < N; ++j) {
    rhs[j] = 4.0 * truth[j] + (j % 2 ? 0.25 : 0.125) * truth[j ^ 1];
    transpose_rhs[j] =
      4.0 * truth[j] + (j % 2 ? 0.125 : 0.25) * truth[j ^ 1];
  }
  kls_options options;
  kls_default_options(&options);
  options.threads = 1;
  options.use_btf = 0;
  options.scale = 0;
  options.ordering = KLS_ORDERING_NATURAL;
  options.orientation = KLS_ORIENTATION_NORMAL;
  options.expected_refactorizations = options.expected_solves = 0;
  kls_solver *solver = NULL;
  assert(kls_create(&solver) == KLS_OK);
  assert(kls_set_accuracy_policy(solver,
    KLS_ACCURACY_COMPONENTWISE_BACKWARD_ERROR) == KLS_OK);
  assert(kls_analyze_csc(solver, KLS_INDEX_INT64, N, ptr, rows, 0,
                        &options) == KLS_OK);
  assert(kls_factor(solver, a) == KLS_OK);
  solver->options.threads = crew;
  solver->tuning.egraph_min_size = 0;
  solver->tuning.egraph_min_flops_per_thread = 0;
  for (int repeat = 0; repeat < 8; ++repeat) {
    double *numeric_values = NULL;
    assert(prepare_numeric_values(solver, a, &numeric_values) == KLS_OK);
    destroy_egraph_refactor_pool(solver);
    assert(kls_egraph_mapped_refactor(solver, numeric_values, 1) == 1);
    assert(solver->egraph_pool->created_count == crew - 1);
    assert(kls_solve(solver, 1, rhs, 0, answer, 0) == KLS_OK);
    for (int j = 0; j < N; ++j) assert(fabs(answer[j] - truth[j]) < 1e-12);
    assert(kls_solve_transpose(solver, 1, transpose_rhs, 0, answer, 0) == KLS_OK);
    for (int j = 0; j < N; ++j) assert(fabs(answer[j] - truth[j]) < 1e-12);
  }
  kls_destroy(solver);
}

int main(void) {
  exercise(2);
  exercise(8);
  puts("PASS cold EGraph numerical startup: crews 2/8, normal/transpose, repeated destruction");
  return 0;
}
