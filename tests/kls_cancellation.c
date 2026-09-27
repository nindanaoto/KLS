/* Nearly floating resistor-chain systems: certify the answer, not the
 * rounding error of a cancellation-prone residual SpMV. Also change the
 * conductances after certification to exercise numeric-epoch invalidation. */
#include "kls/kls.h"
#include <math.h>
#include <stdio.h>

static int check_scaled_chain(int exponent, int zero_initial_rhs) {
  const int ptr[] = {0, 2, 5, 7};
  const int idx[] = {0, 1, 0, 1, 2, 1, 2};
  double a[] = {0.88624328757768067, -0.88624328757768067,
                -0.88624328757768067, 1.7724865851557305,
                -0.88624328757768067, -0.88624328757768067,
                0.88624328757768067};
  for (int p = 0; p < 7; ++p) a[p] = ldexp(a[p], exponent);
  kls_options options;
  kls_default_options(&options);
  options.threads = 1;
  options.orientation = KLS_ORIENTATION_TRANSPOSE;
  kls_solver *solver = NULL;
  if (kls_create(&solver) != KLS_OK) return 1;
  /* This fixture checks strict RHS-relative recovery and its condition probe. */
  if (kls_set_accuracy_policy(solver,KLS_ACCURACY_STRICT_RHS_L2)!=KLS_OK) {
    kls_destroy(solver);return 1;
  }
  int status = kls_analyze_csr(solver, KLS_INDEX_INT32, 3, ptr, idx, 0,
                             &options);
  for (int epoch = 0; status == KLS_OK && epoch < 4; ++epoch) {
    double b[] = {0, -7.0002584799185852e-10, 0}, x[3];
    if (zero_initial_rhs) b[1] = 0.0;
    b[1] = ldexp(b[1], exponent);
    if (epoch) {
      a[3] += ldexp(epoch == 1 ? 0.01 : 1e-8, exponent);
      for (int row = 0; row < 3; ++row) {
        long double rhs = 0;
        for (int p = ptr[row]; p < ptr[row + 1]; ++p)
          rhs += (long double)a[p] * -0.07L;
        b[row] = (double)rhs;
      }
    }
    if (!epoch) {
      status = kls_factor(solver, a);
      if (status == KLS_OK) status = kls_solve(solver, 1, b, 0, x, 0);
    } else {
      status = kls_refactor_solve(solver, a, 1, b, 0, x, 0);
    }
    /* A retained pivot-family policy must not cache a previous RHS answer.
       Exercise a second RHS before advancing to the next numeric epoch. */
    for (int rhs = 0; status == KLS_OK && rhs < 2; ++rhs) {
      if (rhs) {
        for (int row = 0; row < 3; ++row) b[row] = -b[row];
        status = kls_solve(solver, 1, b, 0, x, 0);
        if (status != KLS_OK) break;
      }
      long double r2 = 0, b2 = 0;
      for (int row = 0; row < 3; ++row) {
        long double r = b[row];
        for (int p = ptr[row]; p < ptr[row + 1]; ++p)
          r -= (long double)a[p] * (long double)x[idx[p]];
        r2 += r * r;
        b2 += (long double)b[row] * b[row];
      }
      if (!isfinite(r2) || !(r2 <= 25e-18L * b2)) status = -1;
    }
    if (status != KLS_OK)
      fprintf(stderr, "cancellation epoch %d failed: %d\n", epoch, status);
    if (epoch == 1 && status == KLS_OK) {
      kls_stats stats = {0};
      stats.struct_size = sizeof(stats);
      if (kls_get_stats(solver, &stats) != KLS_OK || !(stats.rcond > 1e-6)) {
        fprintf(stderr, "condition estimate retained from prior epoch\n");
        status = -1;
      }
    }
  }
  kls_destroy(solver);
  return status == KLS_OK ? 0 : 1;
}

int main(void) {
  int failed = check_scaled_chain(0, 1);
  failed |= check_scaled_chain(0, 0);
  failed |= check_scaled_chain(20, 0);
  failed |= check_scaled_chain(-20, 0);
  return failed;
}
