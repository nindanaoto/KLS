/* Private-kernel coverage must not depend on a noisy runtime timing verdict.
   Including the implementation here adds no test switches to the library. */
#include "../src/kls.c"

static int check_answer(int n, const int *ptr, const int *col, const double *a,
                        const double *b, const double *x, const double *truth) {
  long double r2 = 0, b2 = 0;
  for (int row = 0; row < n; ++row) {
    if (!isfinite(x[row]) || fabs(x[row] - truth[row]) > 1e-10) return 0;
    long double residual = b[row];
    for (int p = ptr[row]; p < ptr[row + 1]; ++p)
      residual -= (long double)a[p] * x[col[p]];
    r2 += residual * residual;
    b2 += (long double)b[row] * b[row];
  }
  return isfinite(r2) && r2 <= 1e-24L * b2;
}

static int check_case(int n, int tree, int scale, int exponent, int public_api) {
  int *ptr = malloc(((size_t)n + 1) * sizeof(*ptr));
  const size_t capacity = (tree == 2 ? (size_t)n : 4u) * (size_t)n;
  int *col = malloc(capacity * sizeof(*col));
  double *a = malloc(capacity * sizeof(*a));
  double *b = malloc((size_t)n * sizeof(*b));
  double *x = malloc((size_t)n * sizeof(*x));
  double *truth = malloc((size_t)n * sizeof(*truth));
  if (!ptr || !col || !a || !b || !x || !truth) return 1;
  int nz = 0;
  // Multiplication by two is a bijection for these odd sizes: test a
  // relabelled chain/tree rather than relying on natural node numbering.
  for (int row = 0; row < n; ++row) {
    ptr[row] = nz;
    if (tree == 2) {
      for (int c = 0; c < n; ++c) col[nz++] = c;
      continue;
    }
    const int node = (row + (row % 2 ? n : 0)) / 2;
    col[nz++] = row;
    if (node > 0) col[nz++] = (2 * (tree ? (node - 1) / 2 : node - 1)) % n;
    const int first = tree ? 2 * node + 1 : node + 1;
    if (first < n) col[nz++] = (2 * first) % n;
    if (tree && first + 1 < n) col[nz++] = (2 * (first + 1)) % n;
    for (int i = ptr[row] + 1; i < nz; ++i)
      for (int j = i; j > ptr[row] && col[j] < col[j - 1]; --j) {
        int swap = col[j]; col[j] = col[j - 1]; col[j - 1] = swap;
      }
  }
  ptr[n] = nz;
  kls_options options;
  kls_default_options(&options);
  options.threads = 1;
  options.backend = public_api ? KLS_BACKEND_AUTO : KLS_BACKEND_SERIAL;
  options.orientation = KLS_ORIENTATION_TRANSPOSE;
  options.ordering = KLS_ORDERING_AMD;
  options.use_btf = 0;
  options.scale = scale;
  options.static_pivoting = 0;
  kls_solver *s = NULL;
  int status = kls_create(&s);
  if (status == KLS_OK)
    status = kls_analyze_csr(s, KLS_INDEX_INT32, n, ptr, col, 0, &options);
  for (int epoch = 0; status == KLS_OK && epoch < (public_api ? 64 : 8); ++epoch) {
    for (int row = 0; row < n; ++row) {
      truth[row] = 0.5 + (row % 17) * 0.03125 + epoch * 0.0625;
      for (int p = ptr[row]; p < ptr[row + 1]; ++p)
        a[p] = ldexp(col[p] == row ? (tree == 2 ? n + 4.0 : 4.0) + epoch * 0.125
                                  : -0.125 - (row % 5) * 0.03125, exponent);
    }
    for (int row = 0; row < n; ++row) {
      long double value = 0;
      for (int p = ptr[row]; p < ptr[row + 1]; ++p)
        value += (long double)a[p] * truth[col[p]];
      b[row] = (double)value;
    }
    if (public_api) {
      /* Exercise the public timing coordinator without requiring a particular
         timing verdict. Refactorization must reset any settled trial state;
         alternating aliases and fresh RHSs must not reuse an earlier answer. */
      if (epoch == 0 || epoch == 32) status = kls_factor(s, a);
      if (status != KLS_OK) break;
      if (epoch % 2) memcpy(x, b, (size_t)n * sizeof(*x));
      status = kls_refactor_solve(s, a, 1, epoch % 2 ? x : b, 0, x, 0);
      if (status == KLS_OK && !check_answer(n, ptr, col, a, b, x, truth)) status = -1;
      for (int row = 0; row < n; ++row) { b[row] = -b[row]; truth[row] = -truth[row]; }
      if (status == KLS_OK) status = kls_solve(s, 1, b, 0, x, 0);
      if (status == KLS_OK && !check_answer(n, ptr, col, a, b, x, truth)) status = -1;
      continue;
    }
    if (epoch == 0) {
      status = kls_factor(s, a);
      if (status != KLS_OK) break;
      s->compact_map32_probe_active = 1;
      if (!kls_build_refactor_map(s) || !kls_i32_solve_ready(s)) status = -1;
      s->compact_map32_probe_active = 0;
      s->fused_refactor_solve_work = malloc(2u * (size_t)n * sizeof(double));
      s->fused_refactor_solve_work_n = n;
      if (!s->fused_refactor_solve_work) status = -1;
    }
    if (status != KLS_OK) break;
    s->transpose_cycle.active = 1;
    s->fused_refactor_solve_rhs = b;
    s->compact_map32_probe_active = 1;
    if (kls_lean_single_block_map32_refactor(s, a, 1) != 1) status = -1;
    s->compact_map32_probe_active = 0;
    if (status != KLS_OK) break;
    if (!s->fused_refactor_solve_computed ||
        !kls_compact_transpose_solve(s, b, x) ||
        !check_answer(n, ptr, col, a, b, x, truth)) status = -1;
    // A second RHS must execute U^T again, not reuse the fused partial result.
    for (int row = 0; row < n; ++row) { b[row] = -b[row]; truth[row] = -truth[row]; }
    memcpy(x, b, (size_t)n * sizeof(*x));
    if (!kls_compact_transpose_solve(s, x, x) ||
        !check_answer(n, ptr, col, a, b, x, truth)) status = -1;
    // Standalone refactor then solve, including dirty KLU solve workspace.
    s->compact_map32_probe_active = 1;
    if (kls_lean_single_block_map32_refactor(s, a, 0) != 1) status = -1;
    s->compact_map32_probe_active = 0;
    if (!kls_compact_transpose_solve(s, b, x) ||
        !check_answer(n, ptr, col, a, b, x, truth)) status = -1;
    /* Native indices use no optional mirrors. Exercise the same scaled,
       permuted factors with dirty scratch and an aliased RHS/output. */
    memcpy(x, b, (size_t)n * sizeof(*x));
    if (!kls_native_transpose_solve(s, x) ||
        !check_answer(n, ptr, col, a, b, x, truth)) status = -1;
  }
  if (status != KLS_OK)
    fprintf(stderr, "compact transpose failed n=%d tree=%d scale=%d exponent=%d status=%d\n",
            n, tree, scale, exponent, status);
  kls_destroy(s);
  free(ptr); free(col); free(a); free(b); free(x); free(truth);
  return status != KLS_OK;
}

int main(void) {
  const int sizes[] = {17, 257, 4095, 4097, 65537};
  int failed = 0;
  for (unsigned i = 0; i < sizeof(sizes) / sizeof(sizes[0]); ++i)
    for (int tree = 0; tree < 2; ++tree)
      for (int scale = 0; scale <= 2; ++scale)
        failed |= check_case(sizes[i], tree, scale ? scale : -1, tree ? -20 : 20, 0);
  for (int tree = 0; tree < 2; ++tree)
    for (int scale = 1; scale <= 2; ++scale)
      failed |= check_case(4097, tree, scale, tree ? -20 : 20, 1);
  /* Dense factors exercise the independent accumulators and their tails. */
  for (int n = 7; n <= 33; ++n)
    for (int scale = -1; scale <= 2; ++scale)
      failed |= check_case(n, 2, scale, n % 2 ? -20 : 20, 0);
  return failed;
}
