/* Exercise the production structural planner, without timing or numeric
   dependencies. Publication-order dispatch is covered separately. */
#include "../src/kls.c"

#define CHECK(c) do { if (!(c)) { \
  fprintf(stderr, "batch plan contract failed at line %d\n", __LINE__); \
  return 1; } } while (0)

static int check_padded_prefix(void) {
  enum { N = 80, WIDTH = 5, TAIL = 64 };
  UF_long run_of[N] = {0}, start[1], length[1] = {WIDTH};
  UF_long union_ptr[2] = {0, TAIL}, rows[TAIL], panel_ptr[1] = {0};
  UF_long ui[WIDTH] = {0, 1, 2, 3, 4};
  double panel[WIDTH * TAIL], x[N], expected[N], ux[WIDTH], eu[WIDTH];
  atomic_uint published[N];
  kls_solver solver = {0};
  solver.n = N;
  solver.padded_run_of = run_of;
  solver.padded_run_start = start;
  solver.padded_run_len = length;
  solver.padded_union_ptr = union_ptr;
  solver.padded_union_rows = rows;
  solver.padded_panel_ptr = panel_ptr;
  solver.padded_panel_values = panel;
  kls_egraph_refactor_shared shared = {0};
  shared.solver = &solver;
  shared.pipeline_done = published;
  shared.pipeline_generation = 1;
  atomic_init(&shared.stop, 0);
  for (int r = 0; r < TAIL; ++r) rows[r] = (UF_long)r + 1;
  for (int i = 0; i < WIDTH * TAIL; ++i) panel[i] = 0.001 * (i % 11);
  for (UF_long offset = 0; offset <= 8; offset += 8) {
    memset(run_of, 0, sizeof(run_of));
    start[0] = offset;
    for (int i = 0; i < N; ++i) atomic_init(&published[i], 0);
    for (int i = 0; i < WIDTH; ++i) run_of[offset + i] = 1;
    for (int i = 0; i < 3; ++i) atomic_store(&published[offset + i], 1);
    CHECK(kls_padded_batch_length(&solver, offset, 0, ui, 3, 0) == 3);
    CHECK(kls_padded_batch_length(&solver, offset, 0, ui, 2, 0) == 0);
    CHECK(kls_padded_batch_length(&solver, offset, 0, ui, 0, 0) == 0);
    ui[2] = 4;
    CHECK(kls_padded_batch_length(&solver, offset, 0, ui, 3, 0) == 0);
    ui[2] = 2;
    for (int i = 0; i < N; ++i) x[i] = expected[i] = i * 0.017;
    memset(ux, 0, sizeof(ux));
    memset(eu, 0, sizeof(eu));
    /* Consumer 3 and later run columns are unpublished. Waiting for the
       whole run would deadlock here; only producers 0..2 are required. */
    const kls_cached_batch_result result = kls_egraph_padded_consume(
      &shared, offset, 3, 0, ui, ux, 3, 0, x, 1);
    CHECK(!result.pending && result.consumed == 3);
    CHECK(kls_padded_run_consume(&solver, offset, 0, ui, eu, 3, 0, expected) == 3);
    CHECK(memcmp(x, expected, sizeof(x)) == 0);
    CHECK(memcmp(ux, eu, sizeof(ux)) == 0);
    atomic_store(&published[offset + 1], 0);
    atomic_store(&shared.stop, 1);
    const kls_cached_batch_result stopped = kls_egraph_padded_consume(
      &shared, offset, 3, 0, ui, ux, 3, 0, x, 1);
    CHECK(stopped.pending && stopped.consumed == 0);
    CHECK(memcmp(x, expected, sizeof(x)) == 0);
    CHECK(memcmp(ux, eu, sizeof(ux)) == 0);
    atomic_store(&shared.stop, 0);
  }
  return 0;
}

static int check_fusion(void) {
  enum { N = 160, C = 4 };
  kls_solver solver = {0};
  solver.n = N;
  solver.tuning.snode_min_batch = 2;
  solver.tuning.snode_min_batch_work = 16;
  kls_egraph_refactor_shared shared = {0};
  shared.solver = &solver;
  atomic_init(&shared.stop, 0);
  atomic_uint published[3];
  for (int c = 0; c < 3; ++c) atomic_init(&published[c], 1u);
  shared.pipeline_done = published;
  shared.pipeline_generation = 1u;
  UF_long u[3] = {0, 1, 2}, runs[N] = {3, 3, 3};
  int32_t u32[3] = {0, 1, 2};
  UF_long rows[3][N], lengths[3];
  UF_long *indices[3] = {rows[0], rows[1], rows[2]};
  double values[3][N], *lv[3] = {values[0], values[1], values[2]};
  double initial[C][N], x[C][N], expected[C][N];
  double ux[C][3], eu[C][3];
  const UF_long *ui[C] = {u, u, u, u};
  const int32_t *ui32[C] = {NULL, u32, NULL, u32};
  UF_long lens[C] = {3, 3, 3, 3}, positions[C] = {0};
  double *xv[C] = {x[0], x[1], x[2], x[3]};
  double *uv[C] = {ux[0], ux[1], ux[2], ux[3]};
  const int tails[] = {8, 65, 145};
  for (unsigned trial = 0; trial < 60; ++trial) {
    const int tail = tails[trial % 3];
    for (int c = 0; c < 3; ++c) {
      lengths[c] = (UF_long)(tail + 2 - c);
      for (UF_long r = 0; r < lengths[c]; ++r) {
        rows[c][r] = (UF_long)c + r + 1;
        values[c][r] = ((int)((r * 17 + trial * 13 + c) % 29) - 14) * 0.013;
      }
    }
    for (int c = 0; c < C; ++c) {
      for (int r = 0; r < N; ++r) initial[c][r] = (r - 21 + c) * 0.017;
    }
    memcpy(expected, initial, sizeof(initial));
    memset(eu, 0, sizeof(eu));
    for (int c = 0; c < C; ++c) {
      CHECK(kls_snode_batch_consume_cached(indices, lv, lengths, ui[c],
        ui32[c], eu[c], 3, 0, expected[c], 0, 3, runs, &shared, 1).consumed == 3);
    }
    memcpy(x, initial, sizeof(x));
    memset(ux, 0, sizeof(ux));
    CHECK(kls_snode_batch_consume_cached_pair(indices, lv, lengths,
      u, NULL, ux[0], 3, 0, x[0], u, u32, ux[1], 3, 0, x[1],
      0, 3, runs, &shared) == 3);
    CHECK(memcmp(x, expected, 2 * sizeof(x[0])) == 0);
    CHECK(memcmp(ux, eu, 2 * sizeof(ux[0])) == 0);
    for (int nc = 2; nc <= C; ++nc) {
      memcpy(x, initial, sizeof(x));
      memset(ux, 0, sizeof(ux));
      CHECK(kls_snode_batch_consume_cached_multi(indices, lv, lengths,
        nc, ui, ui32, uv, lens, positions, xv, 0, 3, runs, &shared) == 3);
      CHECK(memcmp(x, expected, (size_t)nc * sizeof(x[0])) == 0);
      CHECK(memcmp(ux, eu, (size_t)nc * sizeof(ux[0])) == 0);
    }
    /* A shorter partner must not silently change the other consumer's
       grouping. Rejection must leave both SPAs and U outputs untouched. */
    memcpy(x, initial, sizeof(x));
    memset(ux, 0, sizeof(ux));
    double saved_u[C][3];
    memcpy(saved_u, ux, sizeof(ux));
    CHECK(kls_snode_batch_consume_cached_pair(indices, lv, lengths,
      u, NULL, ux[0], 3, 0, x[0], u, u32, ux[1], 2, 0, x[1],
      0, 3, runs, &shared) == 0);
    lens[1] = 2;
    CHECK(kls_snode_batch_consume_cached_multi(indices, lv, lengths,
      C, ui, ui32, uv, lens, positions, xv, 0, 3, runs, &shared) == 0);
    lens[1] = 3;
    CHECK(memcmp(x, initial, sizeof(x)) == 0);
    CHECK(memcmp(ux, saved_u, sizeof(ux)) == 0);
    /* Combining consumers cannot promote individually ineligible work. */
    solver.tuning.snode_min_batch_work = (unsigned)(3 * tail + 1);
    CHECK(kls_snode_batch_consume_cached_pair(indices, lv, lengths,
      u, NULL, ux[0], 3, 0, x[0], u, u32, ux[1], 3, 0, x[1],
      0, 3, runs, &shared) == 0);
    CHECK(kls_snode_batch_consume_cached_multi(indices, lv, lengths,
      C, ui, ui32, uv, lens, positions, xv, 0, 3, runs, &shared) == 0);
    CHECK(memcmp(x, initial, sizeof(x)) == 0);
    CHECK(memcmp(ux, saved_u, sizeof(ux)) == 0);
    solver.tuning.snode_min_batch_work = 16;
    atomic_store(&published[1], 0u);
    const kls_cached_batch_result waiting = kls_snode_batch_try_cached(
      indices, lv, lengths, u, NULL, ux[0], 3, 0, x[0], 0, 3,
      runs, &shared, 1);
    CHECK(waiting.pending && waiting.consumed == 0u);
    atomic_store(&shared.stop, 1);
    const kls_cached_batch_result interrupted = kls_snode_batch_consume_cached(
      indices, lv, lengths, u, NULL, ux[0], 3, 0, x[0], 0, 3,
      runs, &shared, 1);
    CHECK(interrupted.pending && interrupted.consumed == 0u);
    atomic_store(&shared.stop, 0);
    const kls_cached_batch_result absent = kls_snode_batch_try_cached(
      indices, lv, lengths, u, NULL, ux[0], 1, 0, x[0], 0, 3,
      runs, &shared, 1);
    CHECK(!absent.pending && absent.consumed == 0u);
    CHECK(kls_snode_batch_consume_cached_pair(indices, lv, lengths,
      u, NULL, ux[0], 3, 0, x[0], u, u32, ux[1], 3, 0, x[1],
      0, 3, runs, &shared) == 0);
    CHECK(kls_snode_batch_consume_cached_multi(indices, lv, lengths,
      C, ui, ui32, uv, lens, positions, xv, 0, 3, runs, &shared) == 0);
    CHECK(memcmp(x, initial, sizeof(x)) == 0);
    CHECK(memcmp(ux, saved_u, sizeof(ux)) == 0);
    atomic_store(&published[1], 1u);
    const kls_cached_batch_result ready = kls_snode_batch_try_cached(
      indices, lv, lengths, u, NULL, ux[0], 3, 0, x[0], 0, 3,
      runs, &shared, 1);
    CHECK(!ready.pending && ready.consumed == 3u);
    CHECK(memcmp(x[0], expected[0], sizeof(x[0])) == 0);
    CHECK(memcmp(ux[0], eu[0], sizeof(ux[0])) == 0);
  }
  return 0;
}

int main(void) {
  enum { N = 128 };
  kls_solver solver = {0};
  solver.tuning.snode_min_batch = 2;
  solver.tuning.snode_min_batch_work = 16;
  UF_long ui[N], runs[N], lengths[N], row = 0;
  int32_t ui32[N];
  UF_long *indices[N];
  for (UF_long i = 0; i < N; ++i) {
    ui[i] = i;
    ui32[i] = (int32_t)i;
    runs[i] = N;
    lengths[i] = 32;
    indices[i] = &row;
  }
  /* Global run and pointer indices; block-relative U and length indices. */
  for (UF_long offset = 0; offset < 8; ++offset) {
    for (UF_long end = 0; end <= 70; ++end) {
      for (UF_long up = 0; up <= end; ++up) {
        for (UF_long limit = up; limit <= end + 1; ++limit) {
          UF_long expected = end - up;
          if (expected > KLS_SNODE_MAX_BATCH) expected = KLS_SNODE_MAX_BATCH;
          if (expected < 2 || up + expected > limit) expected = 0;
          const UF_long a = kls_snode_cached_batch_length(
            &solver, indices, lengths, ui, NULL, end, up, offset, limit, runs);
          const UF_long b = kls_snode_cached_batch_length(
            &solver, indices, lengths, NULL, ui32, end, up, offset, limit, runs);
          CHECK(a == expected && b == expected);
        }
      }
    }
  }
  /* A U gap ends the prefix, even when the producer run continues. */
  ui[3] = 5; ui32[3] = 5;
  CHECK(kls_snode_cached_batch_length(&solver, indices, lengths, ui, NULL,
    8, 0, 0, 8, runs) == 3);
  CHECK(kls_snode_cached_batch_length(&solver, indices, lengths, NULL, ui32,
    8, 0, 0, 8, runs) == 3);
  /* Unequal consumer prefixes stay unequal: no common-prefix clipping. */
  CHECK(kls_snode_cached_batch_length(&solver, indices, lengths, ui, NULL,
    2, 0, 0, 8, runs) == 2);
  lengths[2] = 5;
  CHECK(kls_snode_cached_batch_length(&solver, indices, lengths, ui, NULL,
    8, 0, 0, 8, runs) == 0);
  lengths[2] = 6;
  CHECK(kls_snode_cached_batch_length(&solver, indices, lengths, ui, NULL,
    8, 0, 0, 8, runs) == 3);
  indices[2] = NULL;
  CHECK(kls_snode_cached_batch_length(&solver, indices, lengths, ui, NULL,
    8, 0, 0, 8, runs) == 0);
  indices[2] = &row;
  runs[0] = 1;
  CHECK(kls_snode_cached_batch_length(&solver, indices, lengths, ui, NULL,
    8, 0, 0, 8, runs) == 0);
  runs[0] = 2;
  CHECK(kls_snode_cached_batch_length(&solver, indices, lengths, ui, NULL,
    8, 0, 0, 8, runs) == 2);
  CHECK(check_fusion() == 0);
  CHECK(check_padded_prefix() == 0);
  puts("EGraph cached batch structural and fusion contracts passed");
  return 0;
}
