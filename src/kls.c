#define _POSIX_C_SOURCE 200809L

#include "kls/kls.h"

#include "trilinos_klu_decl.h"

#include <errno.h>
#include <float.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef enum kls_input_format {
  KLS_INPUT_NONE = 0,
  KLS_INPUT_CSC = 1,
  KLS_INPUT_CSR = 2
} kls_input_format;

struct kls_solver {
  UF_long n;
  UF_long nnz;
  UF_long *col_ptr;
  UF_long *row_idx;
  UF_long *input_to_csc;
  double *values;
  double *solve_work;
  size_t solve_work_capacity;
  kls_input_format input_format;
  kls_options options;
  kls_stats stats;
  trilinos_klu_l_common common;
  trilinos_klu_l_symbolic *symbolic;
  trilinos_klu_l_numeric *numeric;
};

static double kls_now_seconds(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static int64_t read_index(kls_index_type type, const void *data, int64_t i) {
  if (type == KLS_INDEX_INT32) {
    const int32_t *p = (const int32_t *)data;
    return (int64_t)p[i];
  }
  const int64_t *p = (const int64_t *)data;
  return p[i];
}

static void free_symbolic(kls_solver *solver) {
  if (solver->symbolic != NULL) {
    trilinos_klu_l_free_symbolic(&solver->symbolic, &solver->common);
    solver->symbolic = NULL;
  }
}

static void free_numeric(kls_solver *solver) {
  if (solver->numeric != NULL) {
    trilinos_klu_l_free_numeric(&solver->numeric, &solver->common);
    solver->numeric = NULL;
  }
}

static void clear_matrix(kls_solver *solver) {
  free_numeric(solver);
  free_symbolic(solver);
  free(solver->col_ptr);
  free(solver->row_idx);
  free(solver->input_to_csc);
  free(solver->values);
  free(solver->solve_work);
  solver->col_ptr = NULL;
  solver->row_idx = NULL;
  solver->input_to_csc = NULL;
  solver->values = NULL;
  solver->solve_work = NULL;
  solver->solve_work_capacity = 0;
  solver->n = 0;
  solver->nnz = 0;
  solver->input_format = KLS_INPUT_NONE;
  memset(&solver->stats, 0, sizeof(solver->stats));
  solver->stats.struct_size = sizeof(solver->stats);
}

static int apply_options_to_common(trilinos_klu_l_common *common, const kls_options *options) {
  if (!trilinos_klu_l_defaults(common)) {
    return KLS_ERR_ANALYZE_FAILED;
  }
  common->btf = options->use_btf ? 1 : 0;
  common->scale = options->scale;
  common->tol = options->pivot_tolerance;
  if (options->memory_growth > 0.0) {
    common->memgrow = options->memory_growth;
  }
  common->halt_if_singular = options->halt_if_singular ? 1 : 0;
  return KLS_OK;
}

static int analyze_with_ordering(kls_solver *solver,
                                 kls_ordering ordering,
                                 trilinos_klu_l_symbolic **symbolic_out,
                                 trilinos_klu_l_common *common_out) {
  trilinos_klu_l_common common;
  int status = apply_options_to_common(&common, &solver->options);
  if (status != KLS_OK) {
    return status;
  }

  trilinos_klu_l_symbolic *symbolic = NULL;
  if (ordering == KLS_ORDERING_NATURAL) {
    symbolic = trilinos_klu_l_analyze_given(solver->n, solver->col_ptr,
                                            solver->row_idx, NULL, NULL,
                                            &common);
  } else {
    common.ordering = (ordering == KLS_ORDERING_COLAMD) ? 1 : 0;
    symbolic = trilinos_klu_l_analyze(solver->n, solver->col_ptr,
                                      solver->row_idx, &common);
  }

  if (symbolic == NULL || common.status < 0) {
    if (symbolic != NULL) {
      trilinos_klu_l_free_symbolic(&symbolic, &common);
    }
    return KLS_ERR_ANALYZE_FAILED;
  }

  *symbolic_out = symbolic;
  *common_out = common;
  return KLS_OK;
}

static double symbolic_score(const trilinos_klu_l_symbolic *symbolic) {
  const double estimated_fill = symbolic->lnz + symbolic->unz;
  if (estimated_fill > 0.0) {
    return estimated_fill;
  }
  if (symbolic->est_flops > 0.0) {
    return symbolic->est_flops;
  }
  return DBL_MAX;
}

static int choose_symbolic(kls_solver *solver) {
  if (solver->options.ordering != KLS_ORDERING_AUTO) {
    int status = analyze_with_ordering(solver, solver->options.ordering,
                                       &solver->symbolic, &solver->common);
    if (status == KLS_OK) {
      solver->stats.selected_ordering = solver->options.ordering;
    }
    return status;
  }

  const kls_ordering candidates[] = {KLS_ORDERING_AMD, KLS_ORDERING_COLAMD};
  trilinos_klu_l_symbolic *best_symbolic = NULL;
  trilinos_klu_l_common best_common;
  kls_ordering best_ordering = KLS_ORDERING_AMD;
  double best_score = 0.0;
  int any_ok = 0;

  for (size_t i = 0; i < sizeof(candidates) / sizeof(candidates[0]); ++i) {
    trilinos_klu_l_symbolic *candidate_symbolic = NULL;
    trilinos_klu_l_common candidate_common;
    int status = analyze_with_ordering(solver, candidates[i],
                                       &candidate_symbolic, &candidate_common);
    if (status != KLS_OK) {
      continue;
    }
    const double score = symbolic_score(candidate_symbolic);
    if (!any_ok || score < best_score) {
      if (best_symbolic != NULL) {
        trilinos_klu_l_free_symbolic(&best_symbolic, &best_common);
      }
      best_symbolic = candidate_symbolic;
      best_common = candidate_common;
      best_ordering = candidates[i];
      best_score = score;
      any_ok = 1;
    } else {
      trilinos_klu_l_free_symbolic(&candidate_symbolic, &candidate_common);
    }
    /* Stop once an ordering provides a real fill estimate; unknown estimates
       should not displace a known-good circuit ordering. */
    if (best_ordering == candidates[i] && best_score < DBL_MAX) {
      break;
    }
  }

  if (!any_ok) {
    return KLS_ERR_ANALYZE_FAILED;
  }
  solver->symbolic = best_symbolic;
  solver->common = best_common;
  solver->stats.selected_ordering = best_ordering;
  return KLS_OK;
}

static int validate_index_base(int index_base) {
  return index_base == 0 || index_base == 1;
}

static int validate_n(int64_t n) {
  return n > 0 && n < INT64_MAX / 2 && n <= (int64_t)(UF_long_max / 2);
}

static int validate_options(const kls_options *options) {
  if (options->ordering < KLS_ORDERING_AUTO || options->ordering > KLS_ORDERING_NATURAL) {
    return 0;
  }
  if (options->scale < -1 || options->scale > 2) {
    return 0;
  }
  if (options->pivot_tolerance < 0.0) {
    return 0;
  }
  if (options->memory_growth <= 0.0) {
    return 0;
  }
  return 1;
}

static int allocate_structure(kls_solver *solver, int64_t n, int64_t nnz, int need_map) {
  solver->col_ptr = (UF_long *)calloc((size_t)n + 1u, sizeof(UF_long));
  solver->row_idx = (UF_long *)calloc((size_t)nnz, sizeof(UF_long));
  if (need_map) {
    solver->input_to_csc = (UF_long *)calloc((size_t)nnz, sizeof(UF_long));
    solver->values = (double *)calloc((size_t)nnz, sizeof(double));
  }
  if (solver->col_ptr == NULL || solver->row_idx == NULL ||
      (need_map && (solver->values == NULL || solver->input_to_csc == NULL))) {
    return KLS_ERR_OUT_OF_MEMORY;
  }
  solver->n = (UF_long)n;
  solver->nnz = (UF_long)nnz;
  return KLS_OK;
}

static void fill_symbolic_stats(kls_solver *solver, double elapsed) {
  solver->stats.struct_size = sizeof(solver->stats);
  solver->stats.n = (int64_t)solver->n;
  solver->stats.nnz = (int64_t)solver->nnz;
  solver->stats.analysis_seconds = elapsed;
  if (solver->symbolic != NULL) {
    solver->stats.last_kernel_status = (int)solver->common.status;
    solver->stats.nblocks = (int64_t)solver->symbolic->nblocks;
    solver->stats.max_block = (int64_t)solver->symbolic->maxblock;
    solver->stats.structural_rank = (int64_t)solver->symbolic->structural_rank;
    solver->stats.estimated_flops = solver->symbolic->est_flops;
    solver->stats.nnz_l = (int64_t)solver->symbolic->lnz;
    solver->stats.nnz_u = (int64_t)solver->symbolic->unz;
  }
  solver->stats.memory_bytes = solver->common.memusage;
  solver->stats.memory_peak_bytes = solver->common.mempeak;
}

static void fill_numeric_stats(kls_solver *solver) {
  solver->stats.last_kernel_status = (int)solver->common.status;
  solver->stats.numerical_rank = (int64_t)solver->common.numerical_rank;
  solver->stats.singular_col = (int64_t)solver->common.singular_col;
  solver->stats.offdiag_pivots = (int64_t)solver->common.noffdiag;
  solver->stats.reallocations = (int64_t)solver->common.nrealloc;
  solver->stats.factor_flops = solver->common.flops;
  solver->stats.rcond = solver->common.rcond;
  solver->stats.condest = solver->common.condest;
  solver->stats.rgrowth = solver->common.rgrowth;
  solver->stats.memory_bytes = solver->common.memusage;
  solver->stats.memory_peak_bytes = solver->common.mempeak;
  if (solver->numeric != NULL) {
    solver->stats.nnz_l = (int64_t)solver->numeric->lnz;
    solver->stats.nnz_u = (int64_t)solver->numeric->unz;
  }
}

void kls_default_options(kls_options *options) {
  if (options == NULL) {
    return;
  }
  memset(options, 0, sizeof(*options));
  options->struct_size = sizeof(*options);
  options->threads = 1;
  options->ordering = KLS_ORDERING_AUTO;
  options->use_btf = 1;
  options->scale = 2;
  options->pivot_tolerance = 0.001;
  options->memory_growth = 1.5;
  options->halt_if_singular = 1;
}

int kls_create(kls_solver **solver_out) {
  if (solver_out == NULL) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  kls_solver *solver = (kls_solver *)calloc(1, sizeof(*solver));
  if (solver == NULL) {
    return KLS_ERR_OUT_OF_MEMORY;
  }
  kls_default_options(&solver->options);
  solver->stats.struct_size = sizeof(solver->stats);
  if (!trilinos_klu_l_defaults(&solver->common)) {
    free(solver);
    return KLS_ERR_ANALYZE_FAILED;
  }
  *solver_out = solver;
  return KLS_OK;
}

void kls_destroy(kls_solver *solver) {
  if (solver == NULL) {
    return;
  }
  clear_matrix(solver);
  free(solver);
}

int kls_analyze_csc(kls_solver *solver,
                    kls_index_type index_type,
                    int64_t n,
                    const void *col_ptr,
                    const void *row_idx,
                    int index_base,
                    const kls_options *options) {
  if (solver == NULL || col_ptr == NULL || row_idx == NULL ||
      !validate_index_base(index_base) || !validate_n(n)) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  if (index_type != KLS_INDEX_INT32 && index_type != KLS_INDEX_INT64) {
    return KLS_ERR_INVALID_ARGUMENT;
  }

  kls_options normalized;
  kls_default_options(&normalized);
  if (options != NULL) {
    if (options->struct_size < sizeof(kls_options)) {
      return KLS_ERR_INVALID_ARGUMENT;
    }
    normalized = *options;
  }
  if (!validate_options(&normalized)) {
    return KLS_ERR_INVALID_ARGUMENT;
  }

  const int64_t base = index_base;
  const int64_t nnz = read_index(index_type, col_ptr, n) - base;
  if (nnz < 0) {
    return KLS_ERR_INVALID_ARGUMENT;
  }

  clear_matrix(solver);
  solver->options = normalized;
  int status = allocate_structure(solver, n, nnz, 0);
  if (status != KLS_OK) {
    clear_matrix(solver);
    return status;
  }

  int64_t prev = 0;
  for (int64_t j = 0; j <= n; ++j) {
    const int64_t ptr = read_index(index_type, col_ptr, j) - base;
    if (ptr < prev || ptr < 0 || ptr > nnz) {
      clear_matrix(solver);
      return KLS_ERR_INVALID_ARGUMENT;
    }
    solver->col_ptr[j] = (UF_long)ptr;
    prev = ptr;
  }
  for (int64_t p = 0; p < nnz; ++p) {
    const int64_t row = read_index(index_type, row_idx, p) - base;
    if (row < 0 || row >= n) {
      clear_matrix(solver);
      return KLS_ERR_INVALID_ARGUMENT;
    }
    solver->row_idx[p] = (UF_long)row;
  }

  solver->input_format = KLS_INPUT_CSC;
  const double start = kls_now_seconds();
  status = choose_symbolic(solver);
  const double elapsed = kls_now_seconds() - start;
  if (status != KLS_OK) {
    clear_matrix(solver);
    return status;
  }
  fill_symbolic_stats(solver, elapsed);
  return KLS_OK;
}

int kls_analyze_csr(kls_solver *solver,
                    kls_index_type index_type,
                    int64_t n,
                    const void *row_ptr,
                    const void *col_idx,
                    int index_base,
                    const kls_options *options) {
  if (solver == NULL || row_ptr == NULL || col_idx == NULL ||
      !validate_index_base(index_base) || !validate_n(n)) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  if (index_type != KLS_INDEX_INT32 && index_type != KLS_INDEX_INT64) {
    return KLS_ERR_INVALID_ARGUMENT;
  }

  kls_options normalized;
  kls_default_options(&normalized);
  if (options != NULL) {
    if (options->struct_size < sizeof(kls_options)) {
      return KLS_ERR_INVALID_ARGUMENT;
    }
    normalized = *options;
  }
  if (!validate_options(&normalized)) {
    return KLS_ERR_INVALID_ARGUMENT;
  }

  const int64_t base = index_base;
  const int64_t nnz = read_index(index_type, row_ptr, n) - base;
  if (nnz < 0) {
    return KLS_ERR_INVALID_ARGUMENT;
  }

  clear_matrix(solver);
  solver->options = normalized;
  int status = allocate_structure(solver, n, nnz, 1);
  if (status != KLS_OK) {
    clear_matrix(solver);
    return status;
  }

  int64_t prev = 0;
  for (int64_t i = 0; i <= n; ++i) {
    const int64_t ptr = read_index(index_type, row_ptr, i) - base;
    if (ptr < prev || ptr < 0 || ptr > nnz) {
      clear_matrix(solver);
      return KLS_ERR_INVALID_ARGUMENT;
    }
    prev = ptr;
  }

  for (int64_t p = 0; p < nnz; ++p) {
    const int64_t col = read_index(index_type, col_idx, p) - base;
    if (col < 0 || col >= n) {
      clear_matrix(solver);
      return KLS_ERR_INVALID_ARGUMENT;
    }
    solver->col_ptr[col + 1]++;
  }
  for (int64_t j = 0; j < n; ++j) {
    solver->col_ptr[j + 1] += solver->col_ptr[j];
  }

  UF_long *next = (UF_long *)malloc(((size_t)n) * sizeof(UF_long));
  if (next == NULL) {
    clear_matrix(solver);
    return KLS_ERR_OUT_OF_MEMORY;
  }
  memcpy(next, solver->col_ptr, ((size_t)n) * sizeof(UF_long));

  for (int64_t row = 0; row < n; ++row) {
    const int64_t begin = read_index(index_type, row_ptr, row) - base;
    const int64_t end = read_index(index_type, row_ptr, row + 1) - base;
    for (int64_t p = begin; p < end; ++p) {
      const int64_t col = read_index(index_type, col_idx, p) - base;
      const UF_long dst = next[col]++;
      solver->row_idx[dst] = (UF_long)row;
      solver->input_to_csc[p] = dst;
    }
  }
  free(next);

  solver->input_format = KLS_INPUT_CSR;
  const double start = kls_now_seconds();
  status = choose_symbolic(solver);
  const double elapsed = kls_now_seconds() - start;
  if (status != KLS_OK) {
    clear_matrix(solver);
    return status;
  }
  fill_symbolic_stats(solver, elapsed);
  return KLS_OK;
}

static int prepare_numeric_values(kls_solver *solver, const double *values, double **values_out) {
  if (solver == NULL || values == NULL || values_out == NULL) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  if (solver->input_format == KLS_INPUT_CSC) {
    *values_out = (double *)values;
    return KLS_OK;
  }
  if (solver->input_format == KLS_INPUT_CSR) {
    if (solver->values == NULL || solver->input_to_csc == NULL) {
      return KLS_ERR_INVALID_ARGUMENT;
    }
    for (int64_t p = 0; p < solver->nnz; ++p) {
      solver->values[solver->input_to_csc[p]] = values[p];
    }
    *values_out = solver->values;
    return KLS_OK;
  }
  return KLS_ERR_INVALID_ARGUMENT;
}

static int ensure_solve_work(kls_solver *solver, int64_t total) {
  if (total <= 0 || (uint64_t)total > SIZE_MAX / sizeof(double)) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  const size_t required = (size_t)total;
  if (solver->solve_work_capacity >= required) {
    return KLS_OK;
  }
  double *next = (double *)realloc(solver->solve_work, required * sizeof(double));
  if (next == NULL) {
    return KLS_ERR_OUT_OF_MEMORY;
  }
  solver->solve_work = next;
  solver->solve_work_capacity = required;
  return KLS_OK;
}

int kls_factor(kls_solver *solver, const double *values) {
  if (solver == NULL || solver->symbolic == NULL || values == NULL) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  double *numeric_values = NULL;
  int status = prepare_numeric_values(solver, values, &numeric_values);
  if (status != KLS_OK) {
    return status;
  }
  free_numeric(solver);
  const double start = kls_now_seconds();
  solver->numeric = trilinos_klu_l_factor(solver->col_ptr, solver->row_idx,
                                          numeric_values, solver->symbolic,
                                          &solver->common);
  solver->stats.factor_seconds = kls_now_seconds() - start;
  if (solver->numeric == NULL || solver->common.status < 0) {
    fill_numeric_stats(solver);
    return solver->common.status == TRILINOS_KLU_SINGULAR ? KLS_ERR_SINGULAR
                                                          : KLS_ERR_FACTOR_FAILED;
  }
  (void)trilinos_klu_l_flops(solver->symbolic, solver->numeric, &solver->common);
  (void)trilinos_klu_l_rcond(solver->symbolic, solver->numeric, &solver->common);
  fill_numeric_stats(solver);
  if (solver->common.status == TRILINOS_KLU_SINGULAR) {
    return KLS_ERR_SINGULAR;
  }
  return KLS_OK;
}

int kls_refactor(kls_solver *solver, const double *values) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL || values == NULL) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  double *numeric_values = NULL;
  int status = prepare_numeric_values(solver, values, &numeric_values);
  if (status != KLS_OK) {
    return status;
  }
  const double start = kls_now_seconds();
  const UF_long ok = trilinos_klu_l_refactor(solver->col_ptr, solver->row_idx,
                                             numeric_values, solver->symbolic,
                                             solver->numeric, &solver->common);
  solver->stats.refactor_seconds = kls_now_seconds() - start;
  (void)trilinos_klu_l_flops(solver->symbolic, solver->numeric, &solver->common);
  (void)trilinos_klu_l_rcond(solver->symbolic, solver->numeric, &solver->common);
  fill_numeric_stats(solver);
  if (!ok || solver->common.status < 0) {
    return solver->common.status == TRILINOS_KLU_SINGULAR ? KLS_ERR_SINGULAR
                                                          : KLS_ERR_REFACTOR_FAILED;
  }
  if (solver->common.status == TRILINOS_KLU_SINGULAR) {
    return KLS_ERR_SINGULAR;
  }
  return KLS_OK;
}

static int solve_impl(kls_solver *solver,
                      int transpose,
                      int64_t nrhs,
                      const double *b,
                      int64_t ldb,
                      double *x,
                      int64_t ldx) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      b == NULL || x == NULL || nrhs <= 0) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  if (ldb == 0) {
    ldb = solver->n;
  }
  if (ldx == 0) {
    ldx = solver->n;
  }
  if (ldb < solver->n || ldx < solver->n) {
    return KLS_ERR_INVALID_ARGUMENT;
  }

  if (nrhs > INT64_MAX / (int64_t)solver->n) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  const int64_t total = (int64_t)solver->n * nrhs;
  int status = ensure_solve_work(solver, total);
  if (status != KLS_OK) {
    return status;
  }
  double *work = solver->solve_work;
  for (int64_t rhs = 0; rhs < nrhs; ++rhs) {
    memcpy(work + rhs * solver->n, b + rhs * ldb, (size_t)solver->n * sizeof(double));
  }

  const double start = kls_now_seconds();
  const UF_long ok = transpose
    ? trilinos_klu_l_tsolve(solver->symbolic, solver->numeric, solver->n,
                            (UF_long)nrhs, work, &solver->common)
    : trilinos_klu_l_solve(solver->symbolic, solver->numeric, solver->n,
                           (UF_long)nrhs, work, &solver->common);
  solver->stats.solve_seconds = kls_now_seconds() - start;
  fill_numeric_stats(solver);

  if (!ok || solver->common.status < 0) {
    return KLS_ERR_SOLVE_FAILED;
  }

  for (int64_t rhs = 0; rhs < nrhs; ++rhs) {
    memcpy(x + rhs * ldx, work + rhs * solver->n, (size_t)solver->n * sizeof(double));
  }
  return KLS_OK;
}

int kls_solve(kls_solver *solver,
              int64_t nrhs,
              const double *b,
              int64_t ldb,
              double *x,
              int64_t ldx) {
  return solve_impl(solver, 0, nrhs, b, ldb, x, ldx);
}

int kls_solve_transpose(kls_solver *solver,
                        int64_t nrhs,
                        const double *b,
                        int64_t ldb,
                        double *x,
                        int64_t ldx) {
  return solve_impl(solver, 1, nrhs, b, ldb, x, ldx);
}

int kls_get_stats(const kls_solver *solver, kls_stats *stats) {
  if (solver == NULL || stats == NULL || stats->struct_size < sizeof(kls_stats)) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  *stats = solver->stats;
  return KLS_OK;
}

const char *kls_status_string(int status) {
  switch (status) {
    case KLS_OK: return "ok";
    case KLS_ERR_INVALID_ARGUMENT: return "invalid argument";
    case KLS_ERR_OUT_OF_MEMORY: return "out of memory";
    case KLS_ERR_ANALYZE_FAILED: return "analysis failed";
    case KLS_ERR_FACTOR_FAILED: return "factorization failed";
    case KLS_ERR_REFACTOR_FAILED: return "refactorization failed";
    case KLS_ERR_SOLVE_FAILED: return "solve failed";
    case KLS_ERR_SINGULAR: return "singular matrix";
    case KLS_ERR_UNSUPPORTED: return "unsupported";
    default: return "unknown status";
  }
}

const char *kls_ordering_name(kls_ordering ordering) {
  switch (ordering) {
    case KLS_ORDERING_AUTO: return "auto";
    case KLS_ORDERING_AMD: return "amd";
    case KLS_ORDERING_COLAMD: return "colamd";
    case KLS_ORDERING_NATURAL: return "natural";
    default: return "unknown";
  }
}
