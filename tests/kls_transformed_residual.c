#define _POSIX_C_SOURCE 200809L
#include "kls/kls.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* A row-permuted, independently row/column-scaled cyclic band matrix.
   Exercise both CSR column widths with changing values and RHS vectors;
   transpose and multiple RHS retain their ordinary fallback paths. */
static int run_case(int32_t n, int scale) {
  int32_t *ap = malloc(((size_t)n + 1) * sizeof(*ap));
  int32_t *ai = malloc(3 * (size_t)n * sizeof(*ai));
  double *ax = malloc(3 * (size_t)n * sizeof(*ax));
  double *b = malloc(2 * (size_t)n * sizeof(*b));
  double *x = malloc(2 * (size_t)n * sizeof(*x));
  double *expected = malloc(2 * (size_t)n * sizeof(*expected));
  kls_solver *solver = NULL;
  int ok = ap && ai && ax && b && x && expected;
  kls_options options;
  kls_default_options(&options);
  options.backend = KLS_BACKEND_KLS;
  options.ordering = KLS_ORDERING_AMD;
  options.orientation = KLS_ORIENTATION_NORMAL;
  options.threads = 4;
  options.scale = scale;
  options.expected_refactorizations = 99;
  options.expected_solves = 100;
  if (!ok) goto done;
  for (int32_t col = 0; col < n; ++col) {
    ap[col] = 3 * col;
    /* Unsorted input intentionally also exercises CSC normalization. */
    for (int k = 0; k < 3; ++k) ai[3 * col + k] = (col + k) % n;
  }
  ap[n] = 3 * n;
  ok = kls_create(&solver) == KLS_OK &&
    kls_analyze_csc(solver, KLS_INDEX_INT32, n, ap, ai, 0, &options) == KLS_OK;
  for (int update = 0; ok && update < 6; ++update) {
    for (int32_t col = 0; col < n; ++col) {
      for (int k = 0; k < 3; ++k) {
        int32_t p = 3 * col + k;
        ax[p] = (k == 1 ? 8.0 : 0.125) *
          (1.0 + 0.1 * (ai[p] % 17)) * (1.0 + 0.1 * (col % 13)) *
          (1.0 + 0.001 * update * (1 + p % 11));
      }
    }
    ok = (update ? kls_refactor(solver, ax) : kls_factor(solver, ax)) == KLS_OK;
    for (int mode = 0; ok && mode < 3; ++mode) {
      const int transpose = mode == 2;
      const int nrhs = mode == 1 ? 2 : 1;
      memset(b, 0, 2 * (size_t)n * sizeof(*b));
      for (int j = 0; j < nrhs; ++j) {
        for (int32_t i = 0; i < n; ++i)
          expected[(size_t)j * n + i] = i % 7 == 0 ? 0.0 :
            0.5 + 0.03125 * ((i + j + update) % 19);
        for (int32_t col = 0; col < n; ++col) {
          for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
            int32_t dst = transpose ? col : ai[p];
            int32_t src = transpose ? ai[p] : col;
            b[(size_t)j * n + dst] += ax[p] * expected[(size_t)j * n + src];
          }
        }
      }
      ok = (transpose ? kls_solve_transpose(solver, nrhs, b, n, x, n)
                      : kls_solve(solver, nrhs, b, n, x, n)) == KLS_OK;
      for (size_t i = 0; ok && i < (size_t)n * nrhs; ++i)
        if (!isfinite(x[i]) || fabs(x[i] - expected[i]) > 1e-9) ok = 0;
    }
  }
done:
  if (!ok) fprintf(stderr, "transformed residual failed: n=%d scale=%d\n", n, scale);
  kls_destroy(solver);
  free(ap); free(ai); free(ax); free(b); free(x); free(expected);
  return ok;
}

int main(void) {
  /* These existing diagnostic controls make coverage deterministic. They
     affect only this standalone test process, never benchmark timings. */
  if (setenv("KLS_FORCE_STATIC_MATCH", "1", 1) != 0 ||
      setenv("KLS_FORCE_GENERIC_PARALLEL_CONTRACT_RESIDUAL", "1", 1) != 0 ||
      setenv("KLS_ENABLE_SOLVE_REFINEMENT", "1", 1) != 0) return EXIT_FAILURE;
  if (!run_case(4096, 0) || !run_case(4096, 2) || !run_case(65536, 2) ||
      setenv("KLS_FORCE_GENERIC_SPRAL_MATCH", "1", 1) != 0) return EXIT_FAILURE;
  return run_case(4096, 2) ? EXIT_SUCCESS : EXIT_FAILURE;
}
