/* Exercise a real numerical dispatch immediately after thread creation.
 * Private admission floors keep this small fixture eligible without changing
 * production policy. Run with TSan to cover the startup publication boundary.
 */
#include "../src/kls.c"
#undef NDEBUG
#include <assert.h>

typedef struct {
  kls_egraph_refactor_pool *pool;
  atomic_int stop;
  atomic_ulong changes;
} budget_writer;

static void *update_idle_budget(void *argument) {
  budget_writer *writer = argument;
  unsigned next = 12345;
  while (!atomic_load_explicit(&writer->stop, memory_order_acquire)) {
    /* Assignment also supplies the deliberately racy unfixed control. */
    writer->pool->worker_spin_iters = next;
    next = next == 12345 ? 23456 : 12345;
    atomic_fetch_add_explicit(&writer->changes, 1ul, memory_order_relaxed);
  }
  return NULL;
}

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
  budget_writer writer = {.pool = solver->egraph_pool};
  atomic_init(&writer.stop, 0);
  atomic_init(&writer.changes, 0ul);
  pthread_t thread;
  assert(pthread_create(&thread, NULL, update_idle_budget, &writer) == 0);
  while (!atomic_load_explicit(&writer.changes, memory_order_relaxed))
    kls_cpu_relax();
  for (int repeat = 0; repeat < 64; ++repeat) {
    double *numeric_values = NULL;
    assert(prepare_numeric_values(solver, a, &numeric_values) == KLS_OK);
    assert(kls_egraph_mapped_refactor(solver, numeric_values, 1) == 1);
  }
  atomic_store_explicit(&writer.stop, 1, memory_order_release);
  assert(pthread_join(thread, NULL) == 0);
  assert(atomic_load_explicit(&writer.changes, memory_order_relaxed) > 0);
  assert(kls_solve(solver, 1, rhs, 0, answer, 0) == KLS_OK);
  for (int j = 0; j < N; ++j) assert(fabs(answer[j] - truth[j]) < 1e-12);
  assert(kls_solve_transpose(solver, 1, transpose_rhs, 0, answer, 0) == KLS_OK);
  for (int j = 0; j < N; ++j) assert(fabs(answer[j] - truth[j]) < 1e-12);
  kls_destroy(solver);
}

int main(void) {
  exercise(2);
  exercise(8);
  puts("PASS cold EGraph numerical startup: crews 2/8, normal/transpose, repeated destruction");
  puts("PASS concurrent EGraph idle-budget updates: crews 2/8, 128 numeric refactors");
  return 0;
}
