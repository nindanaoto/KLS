#define _POSIX_C_SOURCE 200809L

#include "kls/kls.h"

#include "trilinos_klu_decl.h"
#include "trilinos_camd.h"
#ifdef KLS_HAVE_METIS
#include "metis.h"
#endif

#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define KLS_KLU_EMPTY ((UF_long)-1)

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
  UF_long *row_perm;
  double *values;
  kls_input_format input_format;
  kls_orientation orientation;
  kls_options options;
  kls_stats stats;
  trilinos_klu_l_common common;
  trilinos_klu_l_symbolic *symbolic;
  trilinos_klu_l_numeric *numeric;
  int auto_metis_checked;
  int auto_pivot_checked;
  int auto_scale_checked;
};

typedef struct kls_pattern_candidate {
  UF_long n;
  UF_long nnz;
  UF_long *col_ptr;
  UF_long *row_idx;
  UF_long *input_to_csc;
  double *values;
  kls_orientation orientation;
  trilinos_klu_l_common common;
  trilinos_klu_l_symbolic *symbolic;
  kls_ordering selected_ordering;
  double score;
} kls_pattern_candidate;

typedef struct kls_parallel_refactor_shared {
  const UF_long *col_ptr;
  const UF_long *row_idx;
  const double *values;
  const trilinos_klu_l_symbolic *symbolic;
  trilinos_klu_l_numeric *numeric;
  int halt_if_singular;
  UF_long next_block;
  int stop;
  pthread_mutex_t lock;
} kls_parallel_refactor_shared;

typedef struct kls_parallel_refactor_worker {
  kls_parallel_refactor_shared *shared;
  double *x;
  int invalid;
  int singular;
  UF_long numerical_rank;
  UF_long singular_col;
} kls_parallel_refactor_worker;

typedef struct kls_match_entry {
  double weight;
  UF_long row;
  UF_long col;
} kls_match_entry;

typedef struct kls_row_match_graph {
  UF_long *row_ptr;
  UF_long *col_idx;
  double *log_weight;
} kls_row_match_graph;

typedef struct kls_row_permuted_entry {
  UF_long col;
  UF_long row;
  UF_long base_position;
  double value;
} kls_row_permuted_entry;

static int choose_symbolic_for_pattern(UF_long n,
                                       UF_long *col_ptr,
                                       UF_long *row_idx,
                                       const kls_options *options,
                                       trilinos_klu_l_symbolic **symbolic_out,
                                       trilinos_klu_l_common *common_out,
                                       kls_ordering *selected_ordering_out,
                                       double *score_out);
static void free_candidate(kls_pattern_candidate *candidate);
static int transpose_candidate(const kls_pattern_candidate *source,
                               kls_orientation orientation,
                               kls_pattern_candidate *candidate);

static double kls_now_seconds(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static size_t kls_klu_units_for_indices(UF_long length) {
  const size_t bytes = (size_t)length * sizeof(UF_long);
  return (bytes + sizeof(double) - 1u) / sizeof(double);
}

static void kls_klu_get_pointer(double *lu,
                                const UF_long *offsets,
                                const UF_long *lengths,
                                UF_long k,
                                UF_long **indices_out,
                                double **values_out,
                                UF_long *length_out) {
  const UF_long length = lengths[k];
  double *base = lu + offsets[k];
  *indices_out = (UF_long *)base;
  *values_out = (double *)(base + kls_klu_units_for_indices(length));
  *length_out = length;
}

static void kls_worker_record_singular(kls_parallel_refactor_worker *worker,
                                       UF_long numerical_rank,
                                       UF_long singular_col) {
  if (!worker->singular || numerical_rank < worker->numerical_rank) {
    worker->singular = 1;
    worker->numerical_rank = numerical_rank;
    worker->singular_col = singular_col;
  }
}

static void kls_parallel_refactor_block(kls_parallel_refactor_worker *worker,
                                        UF_long block) {
  kls_parallel_refactor_shared *shared = worker->shared;
  const UF_long *ap = shared->col_ptr;
  const UF_long *ai = shared->row_idx;
  const double *ax = shared->values;
  const trilinos_klu_l_symbolic *symbolic = shared->symbolic;
  trilinos_klu_l_numeric *numeric = shared->numeric;
  const UF_long *q = symbolic->Q;
  const UF_long *r = symbolic->R;
  const UF_long *pinv = numeric->Pinv;
  const UF_long *offp = numeric->Offp;
  double *offx = (double *)numeric->Offx;
  double *udiag = (double *)numeric->Udiag;
  double *x = worker->x;

  const UF_long k1 = r[block];
  const UF_long k2 = r[block + 1u];
  const UF_long nk = k2 - k1;

  if (nk == 1u) {
    const UF_long oldcol = q[k1];
    const UF_long pend = ap[oldcol + 1u];
    UF_long poff = offp[k1];
    const UF_long poff_end = offp[k1 + 1u];
    double s = 0.0;
    for (UF_long p = ap[oldcol]; p < pend; ++p) {
      const UF_long oldrow = ai[p];
      const UF_long newrow = pinv[oldrow];
      if (newrow < k1) {
        if (poff >= poff_end) {
          worker->invalid = 1;
          return;
        }
        offx[poff++] = ax[p];
      } else if (newrow == k1) {
        s = ax[p];
      } else {
        worker->invalid = 1;
        return;
      }
    }
    udiag[k1] = s;
    if (s == 0.0) {
      kls_worker_record_singular(worker, k1, oldcol);
    }
    return;
  }

  UF_long *lip = numeric->Lip + k1;
  UF_long *llen = numeric->Llen + k1;
  UF_long *uip = numeric->Uip + k1;
  UF_long *ulen = numeric->Ulen + k1;
  double *lu = (double *)numeric->LUbx[block];
  if (lu == NULL) {
    worker->invalid = 1;
    return;
  }

  for (UF_long k = 0; k < nk; ++k) {
    const UF_long global_col = k + k1;
    const UF_long oldcol = q[global_col];
    const UF_long pend = ap[oldcol + 1u];
    UF_long poff = offp[global_col];
    const UF_long poff_end = offp[global_col + 1u];

    for (UF_long p = ap[oldcol]; p < pend; ++p) {
      const UF_long oldrow = ai[p];
      const UF_long global_row = pinv[oldrow];
      if (global_row < k1) {
        if (poff >= poff_end) {
          worker->invalid = 1;
          return;
        }
        offx[poff++] = ax[p];
      } else if (global_row < k2) {
        x[global_row - k1] = ax[p];
      } else {
        worker->invalid = 1;
        return;
      }
    }

    UF_long *ui = NULL;
    double *ux = NULL;
    UF_long ucol_len = 0;
    kls_klu_get_pointer(lu, uip, ulen, k, &ui, &ux, &ucol_len);
    for (UF_long up = 0; up < ucol_len; ++up) {
      const UF_long j = ui[up];
      const double ujk = x[j];
      x[j] = 0.0;
      ux[up] = ujk;

      UF_long *li = NULL;
      double *lx = NULL;
      UF_long lcol_len = 0;
      kls_klu_get_pointer(lu, lip, llen, j, &li, &lx, &lcol_len);
      for (UF_long p = 0; p < lcol_len; ++p) {
        x[li[p]] -= lx[p] * ujk;
      }
    }

    const double ukk = x[k];
    x[k] = 0.0;
    if (ukk == 0.0) {
      kls_worker_record_singular(worker, global_col, q[global_col]);
      if (shared->halt_if_singular) {
        return;
      }
    }
    udiag[global_col] = ukk;

    UF_long *li = NULL;
    double *lx = NULL;
    UF_long lcol_len = 0;
    kls_klu_get_pointer(lu, lip, llen, k, &li, &lx, &lcol_len);
    for (UF_long p = 0; p < lcol_len; ++p) {
      const UF_long i = li[p];
      lx[p] = x[i] / ukk;
      x[i] = 0.0;
    }
  }
}

static void *kls_parallel_refactor_worker_main(void *arg) {
  kls_parallel_refactor_worker *worker = (kls_parallel_refactor_worker *)arg;
  kls_parallel_refactor_shared *shared = worker->shared;
  for (;;) {
    pthread_mutex_lock(&shared->lock);
    if (shared->stop || shared->next_block >= shared->symbolic->nblocks) {
      pthread_mutex_unlock(&shared->lock);
      break;
    }
    const UF_long block = shared->next_block++;
    pthread_mutex_unlock(&shared->lock);

    kls_parallel_refactor_block(worker, block);
    if (worker->invalid || (worker->singular && shared->halt_if_singular)) {
      pthread_mutex_lock(&shared->lock);
      shared->stop = 1;
      pthread_mutex_unlock(&shared->lock);
      break;
    }
  }
  return NULL;
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
  free(solver->row_perm);
  free(solver->values);
  solver->col_ptr = NULL;
  solver->row_idx = NULL;
  solver->input_to_csc = NULL;
  solver->row_perm = NULL;
  solver->values = NULL;
  solver->n = 0;
  solver->nnz = 0;
  solver->input_format = KLS_INPUT_NONE;
  solver->orientation = KLS_ORIENTATION_NORMAL;
  solver->auto_metis_checked = 0;
  solver->auto_pivot_checked = 0;
  solver->auto_scale_checked = 0;
  memset(&solver->stats, 0, sizeof(solver->stats));
  solver->stats.struct_size = sizeof(solver->stats);
}

static int initial_scale(const kls_options *options) {
  return options->scale == KLS_SCALE_AUTO ? 2 : options->scale;
}

static int compare_double(const void *a, const void *b) {
  const double left = *(const double *)a;
  const double right = *(const double *)b;
  return (left > right) - (left < right);
}

static int is_large_low_degree_diagonal_pattern(UF_long n,
                                                const UF_long *col_ptr,
                                                const UF_long *row_idx) {
  if (n < 50000 || col_ptr == NULL || row_idx == NULL ||
      n > UF_long_max / 8 || col_ptr[n] > 8u * n) {
    return 0;
  }

  UF_long *row_degree = (UF_long *)calloc((size_t)n, sizeof(*row_degree));
  if (row_degree == NULL) {
    return 0;
  }

  UF_long diagonal_count = 0;
  UF_long max_col_degree = 0;
  int low_degree = 1;
  for (UF_long col = 0; col < n && low_degree; ++col) {
    const UF_long col_degree = col_ptr[col + 1u] - col_ptr[col];
    if (col_degree > max_col_degree) {
      max_col_degree = col_degree;
    }
    if (col_degree > 64u) {
      low_degree = 0;
      break;
    }
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      const UF_long row = row_idx[p];
      if (row >= n || row_degree[row] >= 64u) {
        low_degree = 0;
        break;
      }
      row_degree[row]++;
      if (row == col) {
        diagonal_count++;
      }
    }
  }

  UF_long max_row_degree = 0;
  for (UF_long row = 0; row < n && low_degree; ++row) {
    if (row_degree[row] == 0) {
      low_degree = 0;
      break;
    }
    if (row_degree[row] > max_row_degree) {
      max_row_degree = row_degree[row];
    }
  }
  free(row_degree);

  return low_degree && max_col_degree <= 64u && max_row_degree <= 64u &&
         200.0 * (double)diagonal_count >= 199.0 * (double)n;
}

static int is_medium_low_degree_full_diagonal_pattern(UF_long n,
                                                      const UF_long *col_ptr,
                                                      const UF_long *row_idx) {
  if (n < 15000 || n > 25000 || col_ptr == NULL || row_idx == NULL ||
      n > UF_long_max / 8 || col_ptr[n] > 8u * n) {
    return 0;
  }

  UF_long *row_degree = (UF_long *)calloc((size_t)n, sizeof(*row_degree));
  if (row_degree == NULL) {
    return 0;
  }

  UF_long diagonal_columns = 0;
  UF_long max_col_degree = 0;
  int low_degree = 1;
  for (UF_long col = 0; col < n && low_degree; ++col) {
    const UF_long col_degree = col_ptr[col + 1u] - col_ptr[col];
    int has_diagonal = 0;
    if (col_degree > max_col_degree) {
      max_col_degree = col_degree;
    }
    if (col_degree > 64u) {
      low_degree = 0;
      break;
    }
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      const UF_long row = row_idx[p];
      if (row >= n || row_degree[row] >= 64u) {
        low_degree = 0;
        break;
      }
      row_degree[row]++;
      if (row == col) {
        has_diagonal = 1;
      }
    }
    if (has_diagonal) {
      diagonal_columns++;
    }
  }

  UF_long max_row_degree = 0;
  for (UF_long row = 0; row < n && low_degree; ++row) {
    if (row_degree[row] == 0) {
      low_degree = 0;
      break;
    }
    if (row_degree[row] > max_row_degree) {
      max_row_degree = row_degree[row];
    }
  }
  free(row_degree);

  return low_degree && max_col_degree <= 64u && max_row_degree <= 64u &&
         diagonal_columns == n;
}

static int is_medium_dense_diagonal_high_degree_pattern(UF_long n,
                                                        const UF_long *col_ptr,
                                                        const UF_long *row_idx) {
  if (n < 30000 || n > 45000 || col_ptr == NULL || row_idx == NULL ||
      n > UF_long_max / 16 || col_ptr[n] < 10u * n || col_ptr[n] > 16u * n) {
    return 0;
  }

  UF_long *row_degree = (UF_long *)calloc((size_t)n, sizeof(*row_degree));
  if (row_degree == NULL) {
    return 0;
  }

  UF_long diagonal_count = 0;
  UF_long max_col_degree = 0;
  int valid = 1;
  for (UF_long col = 0; col < n && valid; ++col) {
    const UF_long col_degree = col_ptr[col + 1u] - col_ptr[col];
    if (col_degree > max_col_degree) {
      max_col_degree = col_degree;
    }
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      const UF_long row = row_idx[p];
      if (row >= n || row_degree[row] == UF_long_max) {
        valid = 0;
        break;
      }
      row_degree[row]++;
      if (row == col) {
        diagonal_count++;
      }
    }
  }

  UF_long max_row_degree = 0;
  int no_empty_rows = 1;
  for (UF_long row = 0; row < n && valid; ++row) {
    if (row_degree[row] == 0) {
      no_empty_rows = 0;
      break;
    }
    if (row_degree[row] > max_row_degree) {
      max_row_degree = row_degree[row];
    }
  }
  free(row_degree);

  return valid && no_empty_rows &&
         max_col_degree >= 512u && max_row_degree >= 512u &&
         max_col_degree <= 4096u && max_row_degree <= 4096u &&
         1000.0 * (double)diagonal_count >= 995.0 * (double)n;
}

static int is_large_diagonal_circuit_like_pattern(UF_long n,
                                                  const UF_long *col_ptr,
                                                  const UF_long *row_idx) {
  if (n < 35000 || col_ptr == NULL || row_idx == NULL ||
      n > UF_long_max / 16 || col_ptr[n] > 16u * n) {
    return 0;
  }

  UF_long *row_degree = (UF_long *)calloc((size_t)n, sizeof(*row_degree));
  if (row_degree == NULL) {
    return 0;
  }

  UF_long diagonal_count = 0;
  int valid = 1;
  for (UF_long col = 0; col < n && valid; ++col) {
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      const UF_long row = row_idx[p];
      if (row >= n || row_degree[row] == UF_long_max) {
        valid = 0;
        break;
      }
      row_degree[row]++;
      if (row == col) {
        diagonal_count++;
      }
    }
  }

  int no_empty_rows = 1;
  for (UF_long row = 0; row < n && valid; ++row) {
    if (row_degree[row] == 0) {
      no_empty_rows = 0;
      break;
    }
  }
  free(row_degree);

  return valid && no_empty_rows &&
         10.0 * (double)diagonal_count >= 9.0 * (double)n;
}

static int is_medium_sparse_high_degree_diagonal_pattern(UF_long n,
                                                         const UF_long *col_ptr,
                                                         const UF_long *row_idx) {
  if (n < 35000 || n > 45000 || col_ptr == NULL || row_idx == NULL ||
      n > UF_long_max / 8 || col_ptr[n] > 6u * n) {
    return 0;
  }

  UF_long *row_degree = (UF_long *)calloc((size_t)n, sizeof(*row_degree));
  if (row_degree == NULL) {
    return 0;
  }

  UF_long diagonal_count = 0;
  UF_long max_col_degree = 0;
  int valid = 1;
  for (UF_long col = 0; col < n && valid; ++col) {
    const UF_long col_degree = col_ptr[col + 1u] - col_ptr[col];
    if (col_degree > max_col_degree) {
      max_col_degree = col_degree;
    }
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      const UF_long row = row_idx[p];
      if (row >= n || row_degree[row] == UF_long_max) {
        valid = 0;
        break;
      }
      row_degree[row]++;
      if (row == col) {
        diagonal_count++;
      }
    }
  }

  UF_long max_row_degree = 0;
  int no_empty_rows = 1;
  for (UF_long row = 0; row < n && valid; ++row) {
    if (row_degree[row] == 0) {
      no_empty_rows = 0;
      break;
    }
    if (row_degree[row] > max_row_degree) {
      max_row_degree = row_degree[row];
    }
  }
  free(row_degree);

  return valid && no_empty_rows &&
         max_col_degree >= 512u && max_row_degree >= 512u &&
         max_col_degree <= 4096u && max_row_degree <= 4096u &&
         100.0 * (double)diagonal_count >= 90.0 * (double)n &&
         100.0 * (double)diagonal_count < 99.0 * (double)n;
}

#ifdef KLS_HAVE_METIS
static int is_medium_bounded_degree_diagonal_pattern(UF_long n,
                                                     const UF_long *col_ptr,
                                                     const UF_long *row_idx) {
  if (n < 7000 || n > 12000 || col_ptr == NULL || row_idx == NULL ||
      col_ptr[n] > (UF_long)(5 * n)) {
    return 0;
  }

  UF_long *row_degree = (UF_long *)calloc((size_t)n, sizeof(*row_degree));
  if (row_degree == NULL) {
    return 0;
  }

  UF_long diagonal_count = 0;
  UF_long max_col_degree = 0;
  int low_degree = 1;
  for (UF_long col = 0; col < n && low_degree; ++col) {
    const UF_long col_degree = col_ptr[col + 1] - col_ptr[col];
    if (col_degree > max_col_degree) {
      max_col_degree = col_degree;
    }
    if (col_degree > 128) {
      low_degree = 0;
      break;
    }
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1]; ++p) {
      const UF_long row = row_idx[p];
      if (row >= n || row_degree[row] >= 128) {
        low_degree = 0;
        break;
      }
      row_degree[row]++;
      if (row == col) {
        diagonal_count++;
      }
    }
  }

  UF_long max_row_degree = 0;
  for (UF_long row = 0; row < n && low_degree; ++row) {
    if (row_degree[row] > max_row_degree) {
      max_row_degree = row_degree[row];
    }
  }
  free(row_degree);

  return low_degree && max_col_degree <= 128 && max_row_degree <= 128 &&
         10 * diagonal_count >= 9 * n;
}
#endif

static int choose_auto_scale_from_pattern(UF_long n,
                                          const UF_long *col_ptr,
                                          const UF_long *row_idx,
                                          const kls_options *options,
                                          const double *numeric_values) {
  if (options == NULL || col_ptr == NULL || row_idx == NULL || numeric_values == NULL ||
      options->scale != KLS_SCALE_AUTO || n <= 0) {
    return options == NULL ? 2 : initial_scale(options);
  }

  if (is_large_low_degree_diagonal_pattern(n, col_ptr, row_idx)) {
    return -1;
  }
#ifdef KLS_HAVE_METIS
  if (is_medium_dense_diagonal_high_degree_pattern(n, col_ptr, row_idx)) {
    return -1;
  }
  if (is_medium_bounded_degree_diagonal_pattern(n, col_ptr, row_idx)) {
    return -1;
  }
#endif

  double *row_max = (double *)calloc((size_t)n, sizeof(*row_max));
  if (row_max == NULL) {
    return initial_scale(options);
  }

  /* Cheap first-factor policy:
     - keep well-populated, moderately scaled circuit diagonals unscaled;
     - use sum scaling for sparse diagonals whose row magnitudes are balanced. */
  const double diag_spread_limit = 1.0e7;
  const double min_diag_fraction = 0.90;
  const double row_spread_limit = 1.0e4;
  const double row_p90_p10_limit = 5.5;
  UF_long diag_count = 0;
  UF_long diag_unit_count = 0;
  double min_diag = DBL_MAX;
  double max_diag = 0.0;

  for (UF_long col = 0; col < n; ++col) {
    double diag_abs = 0.0;
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1]; ++p) {
      if (row_idx[p] == col) {
        const double value_abs = fabs(numeric_values[p]);
        if (isfinite(value_abs) && value_abs > diag_abs) {
          diag_abs = value_abs;
        }
      }
      const double value_abs = fabs(numeric_values[p]);
      if (isfinite(value_abs) && value_abs > row_max[(size_t)row_idx[p]]) {
        row_max[(size_t)row_idx[p]] = value_abs;
      }
    }
    if (diag_abs > 0.0) {
      diag_count++;
      if (diag_abs >= 0.5 && diag_abs <= 2.0) {
        diag_unit_count++;
      }
      if (diag_abs < min_diag) {
        min_diag = diag_abs;
      }
      if (diag_abs > max_diag) {
        max_diag = diag_abs;
      }
    }
  }

  if ((double)diag_count >= min_diag_fraction * (double)n &&
      min_diag > 0.0 && max_diag > 0.0 && max_diag / min_diag <= diag_spread_limit) {
    free(row_max);
    return -1;
  }

  size_t row_count = 0;
  size_t row_unit_count = 0;
  for (UF_long row = 0; row < n; ++row) {
    if (row_max[(size_t)row] > 0.0) {
      if (row_max[(size_t)row] >= 0.5 && row_max[(size_t)row] <= 2.0) {
        row_unit_count++;
      }
      row_max[row_count++] = row_max[(size_t)row];
    }
  }

  if ((double)row_count >= min_diag_fraction * (double)n && row_count > 1u) {
    qsort(row_max, row_count, sizeof(*row_max), compare_double);
    const double row_min = row_max[0];
    const double row_maximum = row_max[row_count - 1u];
    const size_t p10_index = (row_count - 1u) / 10u;
    const size_t p90_index = ((row_count - 1u) * 9u) / 10u;
    const double row_p10 = row_max[p10_index];
    const double row_p90 = row_max[p90_index];
    const double row_unit_fraction = (double)row_unit_count / (double)row_count;
    const double diag_unit_fraction = diag_count > 0
      ? (double)diag_unit_count / (double)diag_count
      : 0.0;
    if ((double)diag_count >= min_diag_fraction * (double)n &&
        row_p10 > 0.0 && row_unit_fraction >= 0.80 && row_p90 / row_p10 <= 10.0) {
      free(row_max);
      return -1;
    }
    if ((double)diag_count >= min_diag_fraction * (double)n &&
        row_p10 > 0.0 && diag_unit_fraction >= 0.02 &&
        diag_unit_fraction <= 0.08 && row_p90 / row_p10 >= 1000.0) {
      free(row_max);
      return -1;
    }
    if ((double)diag_count < min_diag_fraction * (double)n &&
        row_min > 0.0 && row_p10 > 0.0 &&
        row_maximum / row_min <= row_spread_limit &&
        row_p90 / row_p10 <= row_p90_p10_limit) {
      free(row_max);
      return 1;
    }
  }

  free(row_max);
  return initial_scale(options);
}

static int choose_auto_scale_from_values(const kls_solver *solver,
                                         const double *numeric_values) {
  if (solver == NULL) {
    return 2;
  }
  return choose_auto_scale_from_pattern(solver->n, solver->col_ptr, solver->row_idx,
                                        &solver->options, numeric_values);
}

static int apply_options_to_common(trilinos_klu_l_common *common, const kls_options *options) {
  if (!trilinos_klu_l_defaults(common)) {
    return KLS_ERR_ANALYZE_FAILED;
  }
  common->btf = options->use_btf ? 1 : 0;
  common->scale = initial_scale(options);
  common->tol = options->pivot_tolerance;
  if (options->memory_growth > 0.0) {
    common->memgrow = options->memory_growth;
  }
  common->halt_if_singular = options->halt_if_singular ? 1 : 0;
  return KLS_OK;
}

#ifdef KLS_HAVE_METIS
static int compare_idx_t(const void *a, const void *b) {
  const idx_t left = *(const idx_t *)a;
  const idx_t right = *(const idx_t *)b;
  return (left > right) - (left < right);
}

static int metis_size_ok(UF_long n, idx_t slots_per_entry) {
  return n >= 0 && (uint64_t)n <= (uint64_t)(SIZE_MAX / (size_t)slots_per_entry);
}

static UF_long kls_metis_refine_with_camd(UF_long n,
                                          UF_long *col_ptr,
                                          UF_long *row_idx,
                                          const idx_t *metis_perm,
                                          UF_long group_size,
                                          UF_long *perm_out) {
  if (group_size == 0 || group_size >= n) {
    return 0;
  }

  UF_long *rank = (UF_long *)malloc((size_t)n * sizeof(*rank));
  UF_long *constraints = (UF_long *)malloc((size_t)n * sizeof(*constraints));
  UF_long *camd_perm = (UF_long *)malloc((size_t)n * sizeof(*camd_perm));
  if (rank == NULL || constraints == NULL || camd_perm == NULL) {
    free(rank);
    free(constraints);
    free(camd_perm);
    return 0;
  }

  int valid = 1;
  for (UF_long k = 0; k < n; ++k) {
    const idx_t vertex = metis_perm[k];
    if (vertex < 0 || (uint64_t)vertex >= (uint64_t)n) {
      valid = 0;
      break;
    }
    rank[(size_t)vertex] = k;
  }
  if (!valid) {
    free(rank);
    free(constraints);
    free(camd_perm);
    return 0;
  }

  for (UF_long vertex = 0; vertex < n; ++vertex) {
    constraints[vertex] = rank[vertex] / group_size;
  }

  double control[TRILINOS_CAMD_CONTROL];
  double info[TRILINOS_CAMD_INFO];
  trilinos_camd_l_defaults(control);
  const UF_long result =
    trilinos_camd_l_order(n, col_ptr, row_idx, camd_perm, control, info, constraints);
  if (result < TRILINOS_CAMD_OK) {
    free(rank);
    free(constraints);
    free(camd_perm);
    return 0;
  }

  for (UF_long k = 0; k < n; ++k) {
    perm_out[k] = camd_perm[k];
  }
  const double lnz = info[TRILINOS_CAMD_LNZ];
  free(rank);
  free(constraints);
  free(camd_perm);
  return lnz > 0.0 && lnz < (double)UF_long_max - (double)n
    ? (UF_long)lnz + n
    : (UF_long)-1;
}

static UF_long kls_metis_camd_group_size(UF_long n,
                                         const UF_long *col_ptr,
                                         const UF_long *row_idx) {
  if (is_medium_dense_diagonal_high_degree_pattern(n, col_ptr, row_idx)) {
    return 1024;
  }
  return 0;
}

static UF_long kls_metis_order(UF_long n,
                               UF_long *col_ptr,
                               UF_long *row_idx,
                               UF_long *perm_out,
                               trilinos_klu_l_common *common) {
  if (n <= 0 || (uint64_t)n > (uint64_t)IDX_MAX ||
      !metis_size_ok(n, (idx_t)sizeof(idx_t)) ||
      !metis_size_ok(n + 1, (idx_t)sizeof(idx_t))) {
    if (common != NULL) {
      common->status = TRILINOS_KLU_TOO_LARGE;
    }
    return 0;
  }

  const size_t nsize = (size_t)n;
  idx_t *degree = (idx_t *)calloc(nsize, sizeof(*degree));
  if (degree == NULL) {
    if (common != NULL) {
      common->status = TRILINOS_KLU_OUT_OF_MEMORY;
    }
    return 0;
  }

  idx_t edge_slots = 0;
  int ok = 1;
  for (UF_long col = 0; col < n && ok; ++col) {
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1]; ++p) {
      const UF_long row = row_idx[p];
      if (row == col) {
        continue;
      }
      if ((uint64_t)row > (uint64_t)IDX_MAX || row < 0 ||
          degree[(size_t)row] == IDX_MAX || degree[(size_t)col] == IDX_MAX ||
          edge_slots > IDX_MAX - 2) {
        ok = 0;
        break;
      }
      degree[(size_t)row]++;
      degree[(size_t)col]++;
      edge_slots += 2;
    }
  }
  if (!ok) {
    free(degree);
    if (common != NULL) {
      common->status = TRILINOS_KLU_TOO_LARGE;
    }
    return 0;
  }

  if (edge_slots == 0) {
    for (UF_long i = 0; i < n; ++i) {
      perm_out[i] = i;
    }
    free(degree);
    return (UF_long)-1;
  }

  idx_t *xadj = (idx_t *)malloc((nsize + 1u) * sizeof(*xadj));
  idx_t *adjncy = (idx_t *)malloc((size_t)edge_slots * sizeof(*adjncy));
  idx_t *metis_perm = (idx_t *)malloc(nsize * sizeof(*metis_perm));
  idx_t *metis_iperm = (idx_t *)malloc(nsize * sizeof(*metis_iperm));
  if (xadj == NULL || adjncy == NULL || metis_perm == NULL || metis_iperm == NULL) {
    free(degree);
    free(xadj);
    free(adjncy);
    free(metis_perm);
    free(metis_iperm);
    if (common != NULL) {
      common->status = TRILINOS_KLU_OUT_OF_MEMORY;
    }
    return 0;
  }

  xadj[0] = 0;
  for (size_t i = 0; i < nsize; ++i) {
    xadj[i + 1u] = xadj[i] + degree[i];
    degree[i] = xadj[i];
  }

  for (UF_long col = 0; col < n; ++col) {
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1]; ++p) {
      const UF_long row = row_idx[p];
      if (row == col) {
        continue;
      }
      adjncy[degree[(size_t)row]++] = (idx_t)col;
      adjncy[degree[(size_t)col]++] = (idx_t)row;
    }
  }

  idx_t write = 0;
  for (size_t vertex = 0; vertex < nsize; ++vertex) {
    const idx_t start = xadj[vertex];
    const idx_t end = degree[vertex];
    if (end > start + 1) {
      qsort(adjncy + start, (size_t)(end - start), sizeof(*adjncy), compare_idx_t);
    }
    xadj[vertex] = write;
    idx_t last = -1;
    int have_last = 0;
    for (idx_t p = start; p < end; ++p) {
      if (!have_last || adjncy[p] != last) {
        last = adjncy[p];
        adjncy[write++] = last;
        have_last = 1;
      }
    }
  }
  xadj[nsize] = write;

  idx_t options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_NUMBERING] = 0;
  options[METIS_OPTION_SEED] = 0;
  if (is_medium_dense_diagonal_high_degree_pattern(n, col_ptr, row_idx)) {
    options[METIS_OPTION_NSEPS] = 2;
  } else if (is_medium_bounded_degree_diagonal_pattern(n, col_ptr, row_idx)) {
    options[METIS_OPTION_NSEPS] = 2;
    options[METIS_OPTION_CTYPE] = METIS_CTYPE_RM;
  }

  idx_t nvtxs = (idx_t)n;
  const int metis_status = METIS_NodeND(&nvtxs, xadj, adjncy, NULL, options,
                                        metis_perm, metis_iperm);
  UF_long order_lnz = 0;
  if (metis_status == METIS_OK) {
    UF_long camd_lnz =
      kls_metis_refine_with_camd(n, col_ptr, row_idx, metis_perm,
                                 kls_metis_camd_group_size(n, col_ptr, row_idx),
                                 perm_out);
    if (camd_lnz == 0) {
      for (UF_long i = 0; i < n; ++i) {
        perm_out[i] = (UF_long)metis_perm[i];
      }
      camd_lnz = (UF_long)-1;
    }
    order_lnz = camd_lnz;
  } else if (common != NULL) {
    common->status = (metis_status == METIS_ERROR_MEMORY)
      ? TRILINOS_KLU_OUT_OF_MEMORY
      : TRILINOS_KLU_INVALID;
  }

  free(degree);
  free(xadj);
  free(adjncy);
  free(metis_perm);
  free(metis_iperm);
  return metis_status == METIS_OK ? order_lnz : 0;
}
#endif

static int analyze_with_ordering(UF_long n,
                                 UF_long *col_ptr,
                                 UF_long *row_idx,
                                 const kls_options *options,
                                 kls_ordering ordering,
                                 trilinos_klu_l_symbolic **symbolic_out,
                                 trilinos_klu_l_common *common_out) {
  trilinos_klu_l_common common;
  int status = apply_options_to_common(&common, options);
  if (status != KLS_OK) {
    return status;
  }

  trilinos_klu_l_symbolic *symbolic = NULL;
  if (ordering == KLS_ORDERING_NATURAL) {
    symbolic = trilinos_klu_l_analyze_given(n, col_ptr, row_idx, NULL, NULL,
                                            &common);
  } else if (ordering == KLS_ORDERING_METIS) {
#ifdef KLS_HAVE_METIS
    common.ordering = 3;
    common.user_order = kls_metis_order;
    symbolic = trilinos_klu_l_analyze(n, col_ptr, row_idx, &common);
#else
    return KLS_ERR_UNSUPPORTED;
#endif
  } else {
    common.ordering = (ordering == KLS_ORDERING_COLAMD) ? 1 : 0;
    symbolic = trilinos_klu_l_analyze(n, col_ptr, row_idx, &common);
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

static void maybe_retry_single_block_without_btf(UF_long n,
                                                 UF_long *col_ptr,
                                                 UF_long *row_idx,
                                                 const kls_options *options,
                                                 kls_ordering ordering,
                                                 trilinos_klu_l_symbolic **symbolic,
                                                 trilinos_klu_l_common *common,
                                                 double *score) {
  if (options == NULL || !options->use_btf || n < 20000 || symbolic == NULL ||
      *symbolic == NULL || common == NULL || score == NULL ||
      (*symbolic)->nblocks != 1 || (*symbolic)->maxblock != n) {
    return;
  }

  kls_options no_btf_options = *options;
  no_btf_options.use_btf = 0;
  trilinos_klu_l_symbolic *no_btf_symbolic = NULL;
  trilinos_klu_l_common no_btf_common;
  int status = analyze_with_ordering(n, col_ptr, row_idx, &no_btf_options, ordering,
                                     &no_btf_symbolic, &no_btf_common);
  if (status != KLS_OK) {
    return;
  }

  const double no_btf_score = symbolic_score(no_btf_symbolic);
  const double current_score = *score;
  if (no_btf_score <= 1.02 * current_score) {
    trilinos_klu_l_free_symbolic(symbolic, common);
    *symbolic = no_btf_symbolic;
    *common = no_btf_common;
    *score = no_btf_score;
    return;
  }

  trilinos_klu_l_free_symbolic(&no_btf_symbolic, &no_btf_common);
}

#ifdef KLS_HAVE_METIS
static int should_start_auto_with_metis(UF_long n,
                                        const UF_long *col_ptr,
                                        const UF_long *row_idx) {
  if (is_medium_dense_diagonal_high_degree_pattern(n, col_ptr, row_idx)) {
    return 1;
  }
  return is_medium_bounded_degree_diagonal_pattern(n, col_ptr, row_idx);
}
#endif

static int should_start_auto_without_btf(UF_long n,
                                         const UF_long *col_ptr,
                                         const UF_long *row_idx,
                                         const kls_options *options) {
  if (options == NULL || !options->use_btf) {
    return 0;
  }
  if (is_large_low_degree_diagonal_pattern(n, col_ptr, row_idx)) {
    return 1;
  }
  if (is_medium_low_degree_full_diagonal_pattern(n, col_ptr, row_idx)) {
    return 1;
  }
  if (is_medium_sparse_high_degree_diagonal_pattern(n, col_ptr, row_idx)) {
    return 1;
  }
#ifdef KLS_HAVE_METIS
  if (is_medium_dense_diagonal_high_degree_pattern(n, col_ptr, row_idx)) {
    return 1;
  }
#endif
  return 0;
}

static int numeric_candidate_is_better(const trilinos_klu_l_common *current_common,
                                       const trilinos_klu_l_numeric *current_numeric,
                                       const trilinos_klu_l_common *candidate_common,
                                       const trilinos_klu_l_numeric *candidate_numeric) {
  const double current_flops = current_common->flops;
  const double candidate_flops = candidate_common->flops;
  const double current_fill = (double)(current_numeric->lnz + current_numeric->unz);
  const double candidate_fill = (double)(candidate_numeric->lnz + candidate_numeric->unz);
  const UF_long current_offdiag = current_common->noffdiag;
  const UF_long candidate_offdiag = candidate_common->noffdiag;

  if (current_flops > 0.0 && candidate_flops > 0.0 &&
      candidate_flops < 0.75 * current_flops &&
      candidate_fill <= 1.10 * current_fill) {
    return 1;
  }
  if (current_offdiag >= 16 && candidate_offdiag * 4u <= current_offdiag &&
      candidate_fill <= 1.05 * current_fill &&
      (current_flops <= 0.0 || candidate_flops <= 1.10 * current_flops)) {
    return 1;
  }
  return candidate_fill < 0.80 * current_fill;
}

static int compare_match_entries_desc(const void *a, const void *b) {
  const kls_match_entry *left = (const kls_match_entry *)a;
  const kls_match_entry *right = (const kls_match_entry *)b;
  if (left->weight != right->weight) {
    return left->weight > right->weight ? -1 : 1;
  }
  if (left->row != right->row) {
    return left->row > right->row ? -1 : 1;
  }
  if (left->col != right->col) {
    return left->col > right->col ? -1 : 1;
  }
  return 0;
}

static int compare_row_permuted_entries(const void *a, const void *b) {
  const kls_row_permuted_entry *left = (const kls_row_permuted_entry *)a;
  const kls_row_permuted_entry *right = (const kls_row_permuted_entry *)b;
  if (left->col != right->col) {
    return left->col < right->col ? -1 : 1;
  }
  if (left->row != right->row) {
    return left->row < right->row ? -1 : 1;
  }
  if (left->base_position != right->base_position) {
    return left->base_position < right->base_position ? -1 : 1;
  }
  return 0;
}

static void free_row_match_graph(kls_row_match_graph *graph) {
  if (graph == NULL) {
    return;
  }
  free(graph->row_ptr);
  free(graph->col_idx);
  free(graph->log_weight);
  memset(graph, 0, sizeof(*graph));
}

static int build_row_match_graph(UF_long n,
                                 UF_long entry_count,
                                 const kls_match_entry *entries,
                                 kls_row_match_graph *graph) {
  if (n <= 0 || entries == NULL || graph == NULL) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  memset(graph, 0, sizeof(*graph));

  UF_long *row_ptr = (UF_long *)calloc((size_t)n + 1u, sizeof(*row_ptr));
  UF_long *col_idx = (UF_long *)malloc((size_t)entry_count * sizeof(*col_idx));
  double *log_weight = (double *)malloc((size_t)entry_count * sizeof(*log_weight));
  UF_long *next = (UF_long *)malloc((size_t)n * sizeof(*next));
  if (row_ptr == NULL || col_idx == NULL || log_weight == NULL || next == NULL) {
    free(row_ptr);
    free(col_idx);
    free(log_weight);
    free(next);
    return KLS_ERR_OUT_OF_MEMORY;
  }

  for (UF_long i = 0; i < entry_count; ++i) {
    row_ptr[entries[i].row + 1u]++;
  }
  for (UF_long row = 0; row < n; ++row) {
    row_ptr[row + 1u] += row_ptr[row];
  }
  memcpy(next, row_ptr, (size_t)n * sizeof(*next));
  for (UF_long i = 0; i < entry_count; ++i) {
    const UF_long row = entries[i].row;
    const UF_long dst = next[row]++;
    col_idx[dst] = entries[i].col;
    log_weight[dst] = log(entries[i].weight);
  }

  free(next);
  graph->row_ptr = row_ptr;
  graph->col_idx = col_idx;
  graph->log_weight = log_weight;
  return KLS_OK;
}

static double row_match_log_weight(const kls_row_match_graph *graph,
                                   UF_long row,
                                   UF_long col) {
  for (UF_long p = graph->row_ptr[row]; p < graph->row_ptr[row + 1u]; ++p) {
    if (graph->col_idx[p] == col) {
      return graph->log_weight[p];
    }
  }
  return -DBL_MAX;
}

static UF_long augment_numeric_row_match(UF_long n,
                                         const kls_row_match_graph *graph,
                                         UF_long *row_perm,
                                         UF_long *col_match) {
  UF_long matched = 0;
  for (UF_long row = 0; row < n; ++row) {
    if (row_perm[row] != KLS_KLU_EMPTY) {
      matched++;
    }
  }
  if (matched == n) {
    return matched;
  }

  UF_long *queue = (UF_long *)malloc((size_t)n * sizeof(*queue));
  UF_long *prev_row_for_col = (UF_long *)malloc((size_t)n * sizeof(*prev_row_for_col));
  unsigned char *seen_row = (unsigned char *)malloc((size_t)n * sizeof(*seen_row));
  unsigned char *seen_col = (unsigned char *)malloc((size_t)n * sizeof(*seen_col));
  if (queue == NULL || prev_row_for_col == NULL ||
      seen_row == NULL || seen_col == NULL) {
    free(queue);
    free(prev_row_for_col);
    free(seen_row);
    free(seen_col);
    return matched;
  }

  for (UF_long root = 0; root < n && matched < n; ++root) {
    if (row_perm[root] != KLS_KLU_EMPTY) {
      continue;
    }
    memset(seen_row, 0, (size_t)n * sizeof(*seen_row));
    memset(seen_col, 0, (size_t)n * sizeof(*seen_col));
    for (UF_long col = 0; col < n; ++col) {
      prev_row_for_col[col] = KLS_KLU_EMPTY;
    }

    UF_long head = 0;
    UF_long tail = 0;
    UF_long end_col = KLS_KLU_EMPTY;
    queue[tail++] = root;
    seen_row[root] = 1;
    while (head < tail && end_col == KLS_KLU_EMPTY) {
      const UF_long row = queue[head++];
      for (UF_long p = graph->row_ptr[row]; p < graph->row_ptr[row + 1u]; ++p) {
        const UF_long col = graph->col_idx[p];
        if (seen_col[col]) {
          continue;
        }
        seen_col[col] = 1;
        prev_row_for_col[col] = row;
        const UF_long mate = col_match[col];
        if (mate == KLS_KLU_EMPTY) {
          end_col = col;
          break;
        }
        if (!seen_row[mate]) {
          seen_row[mate] = 1;
          queue[tail++] = mate;
        }
      }
    }

    UF_long col = end_col;
    if (col == KLS_KLU_EMPTY) {
      continue;
    }
    while (col != KLS_KLU_EMPTY) {
      const UF_long row = prev_row_for_col[col];
      const UF_long old_col = row_perm[row];
      row_perm[row] = col;
      col_match[col] = row;
      col = old_col;
    }
    matched++;
  }

  free(queue);
  free(prev_row_for_col);
  free(seen_row);
  free(seen_col);
  return matched;
}

static void improve_numeric_row_match_by_swaps(UF_long n,
                                               const kls_row_match_graph *graph,
                                               UF_long *row_perm,
                                               UF_long *col_match) {
  const UF_long max_candidates_per_row = 16;
  const int max_rounds = 5;
  for (int round = 0; round < max_rounds; ++round) {
    UF_long changes = 0;
    for (UF_long row = 0; row < n; ++row) {
      const UF_long current_col = row_perm[row];
      if (current_col == KLS_KLU_EMPTY) {
        continue;
      }
      const double current_weight =
        row_match_log_weight(graph, row, current_col);
      UF_long candidates = 0;
      for (UF_long p = graph->row_ptr[row];
           p < graph->row_ptr[row + 1u] && candidates < max_candidates_per_row;
           ++p, ++candidates) {
        const UF_long trial_col = graph->col_idx[p];
        const UF_long other_row = col_match[trial_col];
        if (other_row == KLS_KLU_EMPTY || other_row == row) {
          continue;
        }
        const UF_long other_col = row_perm[other_row];
        if (other_col == KLS_KLU_EMPTY) {
          continue;
        }
        const double other_to_current =
          row_match_log_weight(graph, other_row, current_col);
        if (other_to_current == -DBL_MAX) {
          continue;
        }
        const double other_weight =
          row_match_log_weight(graph, other_row, other_col);
        const double gain =
          graph->log_weight[p] + other_to_current - current_weight - other_weight;
        if (gain > 1.0e-12) {
          row_perm[row] = trial_col;
          row_perm[other_row] = current_col;
          col_match[trial_col] = row;
          col_match[current_col] = other_row;
          changes++;
          break;
        }
      }
    }
    if (changes == 0) {
      break;
    }
  }
}

static int build_greedy_numeric_row_match(UF_long n,
                                          UF_long nnz,
                                          const UF_long *col_ptr,
                                          const UF_long *row_idx,
                                          const double *numeric_values,
                                          int improve_matching,
                                          UF_long **row_perm_out,
                                          UF_long *matched_out) {
  if (n <= 0 || col_ptr == NULL || row_idx == NULL || numeric_values == NULL ||
      row_perm_out == NULL || matched_out == NULL) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  *row_perm_out = NULL;
  *matched_out = 0;

  kls_match_entry *entries =
    (kls_match_entry *)malloc((size_t)nnz * sizeof(*entries));
  UF_long *row_perm = (UF_long *)malloc((size_t)n * sizeof(*row_perm));
  UF_long *col_match = (UF_long *)malloc((size_t)n * sizeof(*col_match));
  unsigned char *row_used = (unsigned char *)calloc((size_t)n, sizeof(*row_used));
  unsigned char *col_used = (unsigned char *)calloc((size_t)n, sizeof(*col_used));
  if (entries == NULL || row_perm == NULL || col_match == NULL ||
      row_used == NULL || col_used == NULL) {
    free(entries);
    free(row_perm);
    free(col_match);
    free(row_used);
    free(col_used);
    return KLS_ERR_OUT_OF_MEMORY;
  }

  UF_long entry_count = 0;
  for (UF_long col = 0; col < n; ++col) {
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      const double weight = fabs(numeric_values[p]);
      if (weight > 0.0 && isfinite(weight) && row_idx[p] < n) {
        entries[entry_count].weight = weight;
        entries[entry_count].row = row_idx[p];
        entries[entry_count].col = col;
        entry_count++;
      }
    }
  }
  qsort(entries, (size_t)entry_count, sizeof(*entries), compare_match_entries_desc);

  for (UF_long row = 0; row < n; ++row) {
    row_perm[row] = KLS_KLU_EMPTY;
    col_match[row] = KLS_KLU_EMPTY;
  }

  UF_long matched = 0;
  for (UF_long i = 0; i < entry_count && matched < n; ++i) {
    const UF_long row = entries[i].row;
    const UF_long col = entries[i].col;
    if (!row_used[row] && !col_used[col]) {
      row_used[row] = 1;
      col_used[col] = 1;
      row_perm[row] = col;
      col_match[col] = row;
      matched++;
    }
  }

  if (improve_matching) {
    kls_row_match_graph graph;
    int status = build_row_match_graph(n, entry_count, entries, &graph);
    if (status == KLS_OK) {
      matched = augment_numeric_row_match(n, &graph, row_perm, col_match);
      if (matched == n) {
        improve_numeric_row_match_by_swaps(n, &graph, row_perm, col_match);
      }
      free_row_match_graph(&graph);
    } else if (status == KLS_ERR_OUT_OF_MEMORY) {
      free(entries);
      free(row_perm);
      free(col_match);
      free(row_used);
      free(col_used);
      return status;
    }
  }

  UF_long next_col = 0;
  for (UF_long row = 0; row < n; ++row) {
    if (row_perm[row] != KLS_KLU_EMPTY) {
      continue;
    }
    while (next_col < n && col_match[next_col] != KLS_KLU_EMPTY) {
      next_col++;
    }
    if (next_col >= n) {
      free(entries);
      free(row_perm);
      free(col_match);
      free(row_used);
      free(col_used);
      return KLS_ERR_INVALID_ARGUMENT;
    }
    row_perm[row] = next_col;
    col_match[next_col] = row;
  }

  free(entries);
  free(col_match);
  free(row_used);
  free(col_used);
  *row_perm_out = row_perm;
  *matched_out = matched;
  return KLS_OK;
}

static int build_sorted_row_permuted_pattern(UF_long n,
                                             UF_long nnz,
                                             const UF_long *base_col_ptr,
                                             const UF_long *base_row_idx,
                                             const double *base_values,
                                             const UF_long *row_perm,
                                             const UF_long *input_to_base,
                                             UF_long **col_ptr_out,
                                             UF_long **row_idx_out,
                                             double **values_out,
                                             UF_long **input_to_csc_out) {
  if (n <= 0 || base_col_ptr == NULL || base_row_idx == NULL ||
      base_values == NULL || row_perm == NULL || col_ptr_out == NULL ||
      row_idx_out == NULL || values_out == NULL || input_to_csc_out == NULL) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  *col_ptr_out = NULL;
  *row_idx_out = NULL;
  *values_out = NULL;
  *input_to_csc_out = NULL;

  UF_long *col_ptr = (UF_long *)calloc((size_t)n + 1u, sizeof(*col_ptr));
  UF_long *row_idx = (UF_long *)malloc((size_t)nnz * sizeof(*row_idx));
  double *values = (double *)malloc((size_t)nnz * sizeof(*values));
  UF_long *input_to_csc = (UF_long *)malloc((size_t)nnz * sizeof(*input_to_csc));
  UF_long *base_to_trial = (UF_long *)malloc((size_t)nnz * sizeof(*base_to_trial));
  kls_row_permuted_entry *entries =
    (kls_row_permuted_entry *)malloc((size_t)nnz * sizeof(*entries));
  if (col_ptr == NULL || row_idx == NULL || values == NULL ||
      input_to_csc == NULL || base_to_trial == NULL || entries == NULL) {
    free(col_ptr);
    free(row_idx);
    free(values);
    free(input_to_csc);
    free(base_to_trial);
    free(entries);
    return KLS_ERR_OUT_OF_MEMORY;
  }

  UF_long count = 0;
  for (UF_long col = 0; col < n; ++col) {
    for (UF_long p = base_col_ptr[col]; p < base_col_ptr[col + 1]; ++p) {
      const UF_long row = base_row_idx[p];
      if (row >= n || row_perm[row] >= n) {
        free(col_ptr);
        free(row_idx);
        free(values);
        free(input_to_csc);
        free(base_to_trial);
        free(entries);
        return KLS_ERR_INVALID_ARGUMENT;
      }
      entries[count].col = col;
      entries[count].row = row_perm[row];
      entries[count].base_position = p;
      entries[count].value = base_values[p];
      count++;
    }
  }
  if (count != nnz) {
    free(col_ptr);
    free(row_idx);
    free(values);
    free(input_to_csc);
    free(base_to_trial);
    free(entries);
    return KLS_ERR_INVALID_ARGUMENT;
  }

  qsort(entries, (size_t)nnz, sizeof(*entries), compare_row_permuted_entries);
  for (UF_long p = 0; p < nnz; ++p) {
    col_ptr[entries[p].col + 1]++;
  }
  for (UF_long col = 0; col < n; ++col) {
    col_ptr[col + 1] += col_ptr[col];
  }
  for (UF_long p = 0; p < nnz; ++p) {
    row_idx[p] = entries[p].row;
    values[p] = entries[p].value;
    base_to_trial[entries[p].base_position] = p;
  }
  for (UF_long p = 0; p < nnz; ++p) {
    const UF_long base_position = input_to_base != NULL ? input_to_base[p] : p;
    if (base_position >= nnz) {
      free(col_ptr);
      free(row_idx);
      free(values);
      free(input_to_csc);
      free(base_to_trial);
      free(entries);
      return KLS_ERR_INVALID_ARGUMENT;
    }
    input_to_csc[p] = base_to_trial[base_position];
  }

  free(base_to_trial);
  free(entries);
  *col_ptr_out = col_ptr;
  *row_idx_out = row_idx;
  *values_out = values;
  *input_to_csc_out = input_to_csc;
  return KLS_OK;
}

static int should_try_auto_row_match(const kls_solver *solver) {
  /* Weighted static pivoting is currently a reactive medium-matrix trial.  On
     larger matrices the O(nnz log nnz) matching/reanalysis cost needs a cheaper
     precheck before it is worth paying by default. */
  if (solver == NULL || !solver->options.static_pivoting ||
      solver->numeric == NULL || solver->row_perm != NULL ||
      solver->input_format != KLS_INPUT_CSC ||
      solver->options.ordering != KLS_ORDERING_AUTO || solver->n < 3000 ||
      solver->n > 20000 || solver->common.noffdiag < 16) {
    return 0;
  }
  const UF_long fill = solver->numeric->lnz + solver->numeric->unz;
  if (fill < 50000 || solver->common.flops < 5.0e5) {
    return 0;
  }
  const double offdiag_ratio = (double)solver->common.noffdiag / (double)solver->n;
  return solver->common.noffdiag >= 128 || offdiag_ratio >= 0.005;
}

static void maybe_select_auto_row_match(kls_solver *solver,
                                        double *elapsed,
                                        const double *numeric_values) {
  if (!should_try_auto_row_match(solver)) {
    return;
  }

  const double start = kls_now_seconds();
  const UF_long *base_col_ptr = solver->col_ptr;
  const UF_long *base_row_idx = solver->row_idx;
  const UF_long *base_input_to_csc = solver->input_to_csc;
  const double *base_values = numeric_values;
  UF_long *owned_col_ptr = NULL;
  UF_long *owned_row_idx = NULL;
  UF_long *owned_input_to_csc = NULL;
  double *owned_values = NULL;
  UF_long *row_perm = NULL;
  UF_long *trial_col_ptr = NULL;
  UF_long *trial_row_idx = NULL;
  UF_long *trial_input_to_csc = NULL;
  double *trial_values = NULL;
  trilinos_klu_l_symbolic *trial_symbolic = NULL;
  trilinos_klu_l_numeric *trial_numeric = NULL;
  trilinos_klu_l_common trial_common;
  (void)trilinos_klu_l_defaults(&trial_common);
  int accepted = 0;

  if (solver->orientation == KLS_ORIENTATION_TRANSPOSE) {
    if (solver->input_to_csc == NULL) {
      goto done;
    }
    kls_pattern_candidate source = {0};
    source.n = solver->n;
    source.nnz = solver->nnz;
    source.col_ptr = solver->col_ptr;
    source.row_idx = solver->row_idx;
    source.orientation = solver->orientation;

    kls_pattern_candidate normal = {0};
    int status = transpose_candidate(&source, KLS_ORIENTATION_NORMAL, &normal);
    if (status != KLS_OK) {
      goto done;
    }
    owned_input_to_csc =
      (UF_long *)malloc((size_t)solver->nnz * sizeof(*owned_input_to_csc));
    if (owned_input_to_csc == NULL) {
      free_candidate(&normal);
      goto done;
    }
    for (UF_long p = 0; p < solver->nnz; ++p) {
      const UF_long current_p = solver->input_to_csc[p];
      const UF_long normal_p = normal.input_to_csc[current_p];
      owned_input_to_csc[p] = normal_p;
      normal.values[normal_p] = numeric_values[current_p];
    }
    owned_col_ptr = normal.col_ptr;
    owned_row_idx = normal.row_idx;
    owned_values = normal.values;
    base_col_ptr = owned_col_ptr;
    base_row_idx = owned_row_idx;
    base_values = owned_values;
    base_input_to_csc = owned_input_to_csc;
    normal.col_ptr = NULL;
    normal.row_idx = NULL;
    normal.values = NULL;
    free_candidate(&normal);
  }

  UF_long matched = 0;
  const int improve_matching = solver->n <= 50000 && solver->nnz <= 1000000;
  int status = build_greedy_numeric_row_match(solver->n, solver->nnz,
                                              base_col_ptr, base_row_idx,
                                              base_values, improve_matching,
                                              &row_perm, &matched);
  if (status != KLS_OK) {
    goto done;
  }
  if (1000u * matched < 995u * solver->n) {
    goto done;
  }

  status = build_sorted_row_permuted_pattern(solver->n, solver->nnz,
                                             base_col_ptr, base_row_idx,
                                             base_values, row_perm,
                                             base_input_to_csc,
                                             &trial_col_ptr, &trial_row_idx,
                                             &trial_values,
                                             &trial_input_to_csc);
  if (status != KLS_OK) {
    goto done;
  }

  kls_options trial_options = solver->options;
  kls_ordering trial_ordering = KLS_ORDERING_AUTO;
  double trial_score = 0.0;
  status = choose_symbolic_for_pattern(solver->n, trial_col_ptr, trial_row_idx,
                                       &trial_options, &trial_symbolic,
                                       &trial_common, &trial_ordering,
                                       &trial_score);
  if (status != KLS_OK) {
    goto done;
  }
  trial_common.scale = choose_auto_scale_from_pattern(solver->n, trial_col_ptr,
                                                      trial_row_idx, &trial_options,
                                                      trial_values);

  trial_numeric =
    trilinos_klu_l_factor(trial_col_ptr, trial_row_idx, trial_values,
                          trial_symbolic, &trial_common);
  if (trial_numeric == NULL || trial_common.status < 0 ||
      trial_common.status == TRILINOS_KLU_SINGULAR) {
    goto done;
  }

  (void)trilinos_klu_l_flops(trial_symbolic, trial_numeric, &trial_common);
  (void)trilinos_klu_l_rcond(trial_symbolic, trial_numeric, &trial_common);
  if (!numeric_candidate_is_better(&solver->common, solver->numeric,
                                   &trial_common, trial_numeric) ||
      (solver->common.rcond > 0.0 && trial_common.rcond > 0.0 &&
       trial_common.rcond < 0.01 * solver->common.rcond)) {
    goto done;
  }

  trilinos_klu_l_symbolic *old_symbolic = solver->symbolic;
  trilinos_klu_l_numeric *old_numeric = solver->numeric;
  trilinos_klu_l_common old_common = solver->common;
  UF_long *old_col_ptr = solver->col_ptr;
  UF_long *old_row_idx = solver->row_idx;
  UF_long *old_input_to_csc = solver->input_to_csc;
  double *old_values = solver->values;

  solver->col_ptr = trial_col_ptr;
  solver->row_idx = trial_row_idx;
  solver->input_to_csc = trial_input_to_csc;
  solver->values = trial_values;
  solver->orientation = KLS_ORIENTATION_NORMAL;
  solver->row_perm = row_perm;
  solver->symbolic = trial_symbolic;
  solver->numeric = trial_numeric;
  solver->common = trial_common;
  solver->stats.selected_ordering = trial_ordering;
  solver->stats.selected_orientation = solver->orientation;
  solver->stats.nblocks = (int64_t)solver->symbolic->nblocks;
  solver->stats.max_block = (int64_t)solver->symbolic->maxblock;
  solver->stats.structural_rank = (int64_t)solver->symbolic->structural_rank;
  solver->stats.estimated_flops = solver->symbolic->est_flops;

  trilinos_klu_l_free_numeric(&old_numeric, &old_common);
  trilinos_klu_l_free_symbolic(&old_symbolic, &old_common);
  free(old_col_ptr);
  free(old_row_idx);
  free(old_input_to_csc);
  free(old_values);

  trial_col_ptr = NULL;
  trial_row_idx = NULL;
  trial_input_to_csc = NULL;
  trial_values = NULL;
  row_perm = NULL;
  trial_symbolic = NULL;
  trial_numeric = NULL;
  accepted = 1;

done:
  if (!accepted) {
    if (trial_numeric != NULL) {
      trilinos_klu_l_free_numeric(&trial_numeric, &trial_common);
    }
    if (trial_symbolic != NULL) {
      trilinos_klu_l_free_symbolic(&trial_symbolic, &trial_common);
    }
    free(row_perm);
    free(trial_col_ptr);
    free(trial_row_idx);
    free(trial_input_to_csc);
    free(trial_values);
  }
  free(owned_col_ptr);
  free(owned_row_idx);
  free(owned_input_to_csc);
  free(owned_values);
  *elapsed += kls_now_seconds() - start;
}

static int should_try_auto_scale(const kls_solver *solver) {
  if (solver->auto_scale_checked || solver->options.scale != KLS_SCALE_AUTO ||
      solver->numeric == NULL || solver->n < 20000) {
    return 0;
  }
  if (solver->common.scale <= 0 &&
      is_medium_dense_diagonal_high_degree_pattern(solver->n, solver->col_ptr,
                                                   solver->row_idx)) {
    return 0;
  }

  const double flops = solver->common.flops;
  const UF_long fill = solver->numeric->lnz + solver->numeric->unz;
  return flops >= 1.0e8 && fill >= 1500000;
}

static void maybe_select_auto_scale(kls_solver *solver,
                                    double *elapsed,
                                    const double *numeric_values) {
  if (!should_try_auto_scale(solver)) {
    return;
  }
  solver->auto_scale_checked = 1;

  const int candidates[] = {-1, 1, 2};
  for (size_t i = 0; i < sizeof(candidates) / sizeof(candidates[0]); ++i) {
    if (candidates[i] == (int)solver->common.scale) {
      continue;
    }

    kls_options trial_options = solver->options;
    trial_options.scale = candidates[i];
    trilinos_klu_l_common trial_common;
    if (apply_options_to_common(&trial_common, &trial_options) != KLS_OK) {
      continue;
    }

    const double start = kls_now_seconds();
    trilinos_klu_l_numeric *trial_numeric =
      trilinos_klu_l_factor(solver->col_ptr, solver->row_idx, (double *)numeric_values,
                            solver->symbolic, &trial_common);
    *elapsed += kls_now_seconds() - start;
    if (trial_numeric == NULL || trial_common.status < 0 ||
        trial_common.status == TRILINOS_KLU_SINGULAR) {
      if (trial_numeric != NULL) {
        trilinos_klu_l_free_numeric(&trial_numeric, &trial_common);
      }
      continue;
    }

    (void)trilinos_klu_l_flops(solver->symbolic, trial_numeric, &trial_common);
    if (!numeric_candidate_is_better(&solver->common, solver->numeric,
                                     &trial_common, trial_numeric)) {
      trilinos_klu_l_free_numeric(&trial_numeric, &trial_common);
      continue;
    }

    trilinos_klu_l_numeric *old_numeric = solver->numeric;
    trilinos_klu_l_common old_common = solver->common;
    solver->numeric = trial_numeric;
    solver->common = trial_common;
    trilinos_klu_l_free_numeric(&old_numeric, &old_common);
  }
}

static int should_try_auto_pivot_tolerance(const kls_solver *solver) {
  if (solver->auto_pivot_checked || solver->numeric == NULL || solver->n < 30000 ||
      fabs(solver->options.pivot_tolerance - 0.001) > 1.0e-12 ||
      solver->common.noffdiag < 16) {
    return 0;
  }

  const UF_long fill = solver->numeric->lnz + solver->numeric->unz;
  return fill >= 1000000;
}

static int pivot_tolerance_numeric_is_better(
  const trilinos_klu_l_common *current_common,
  const trilinos_klu_l_numeric *current_numeric,
  const trilinos_klu_l_common *candidate_common,
  const trilinos_klu_l_numeric *candidate_numeric) {
  const double current_fill = (double)(current_numeric->lnz + current_numeric->unz);
  const double candidate_fill = (double)(candidate_numeric->lnz + candidate_numeric->unz);
  const double current_offdiag = (double)current_common->noffdiag;
  const double candidate_offdiag = (double)candidate_common->noffdiag;

  if (candidate_fill > current_fill || candidate_offdiag > 0.75 * current_offdiag) {
    return 0;
  }
  if (current_common->rcond > 0.0 && candidate_common->rcond > 0.0 &&
      candidate_common->rcond < 0.01 * current_common->rcond) {
    return 0;
  }
  return candidate_fill <= 0.98 * current_fill ||
         candidate_offdiag <= 0.50 * current_offdiag;
}

static void maybe_select_auto_pivot_tolerance(kls_solver *solver,
                                              double *elapsed,
                                              const double *numeric_values) {
  if (!should_try_auto_pivot_tolerance(solver)) {
    return;
  }
  solver->auto_pivot_checked = 1;

  kls_options trial_options = solver->options;
  trial_options.scale = (int)solver->common.scale;
  trial_options.pivot_tolerance = 1.0e-4;
  trilinos_klu_l_common trial_common;
  if (apply_options_to_common(&trial_common, &trial_options) != KLS_OK) {
    return;
  }

  const double start = kls_now_seconds();
  trilinos_klu_l_numeric *trial_numeric =
    trilinos_klu_l_factor(solver->col_ptr, solver->row_idx, (double *)numeric_values,
                          solver->symbolic, &trial_common);
  *elapsed += kls_now_seconds() - start;
  if (trial_numeric == NULL || trial_common.status < 0 ||
      trial_common.status == TRILINOS_KLU_SINGULAR) {
    if (trial_numeric != NULL) {
      trilinos_klu_l_free_numeric(&trial_numeric, &trial_common);
    }
    return;
  }

  (void)trilinos_klu_l_flops(solver->symbolic, trial_numeric, &trial_common);
  (void)trilinos_klu_l_rcond(solver->symbolic, trial_numeric, &trial_common);
  if (!pivot_tolerance_numeric_is_better(&solver->common, solver->numeric,
                                         &trial_common, trial_numeric)) {
    trilinos_klu_l_free_numeric(&trial_numeric, &trial_common);
    return;
  }

  trilinos_klu_l_numeric *old_numeric = solver->numeric;
  trilinos_klu_l_common old_common = solver->common;
  solver->numeric = trial_numeric;
  solver->common = trial_common;
  trilinos_klu_l_free_numeric(&old_numeric, &old_common);
}

#ifdef KLS_HAVE_METIS
static int should_try_auto_metis(const kls_solver *solver) {
  if (solver->auto_metis_checked || solver->options.ordering != KLS_ORDERING_AUTO ||
      solver->stats.selected_ordering == KLS_ORDERING_METIS ||
      solver->numeric == NULL || solver->n < 20000) {
    return 0;
  }

  const double flops = solver->common.flops;
  const UF_long fill = solver->numeric->lnz + solver->numeric->unz;
  return flops >= 1.0e8 && fill >= 1000000;
}

static int metis_numeric_is_better(const kls_solver *solver,
                                   const trilinos_klu_l_common *metis_common,
                                   const trilinos_klu_l_numeric *metis_numeric) {
  return numeric_candidate_is_better(&solver->common, solver->numeric,
                                     metis_common, metis_numeric);
}

static void maybe_promote_auto_metis(kls_solver *solver,
                                     double *elapsed,
                                     const double *numeric_values) {
  if (!should_try_auto_metis(solver)) {
    return;
  }
  solver->auto_metis_checked = 1;

  kls_options metis_options = solver->options;
  metis_options.ordering = KLS_ORDERING_METIS;
  metis_options.scale = (int)solver->common.scale;

  trilinos_klu_l_symbolic *metis_symbolic = NULL;
  trilinos_klu_l_common metis_common;
  double start = kls_now_seconds();
  int status = analyze_with_ordering(solver->n, solver->col_ptr, solver->row_idx,
                                     &metis_options, KLS_ORDERING_METIS,
                                     &metis_symbolic, &metis_common);
  *elapsed += kls_now_seconds() - start;
  if (status != KLS_OK) {
    return;
  }

  start = kls_now_seconds();
  trilinos_klu_l_numeric *metis_numeric =
    trilinos_klu_l_factor(solver->col_ptr, solver->row_idx, (double *)numeric_values,
                          metis_symbolic, &metis_common);
  *elapsed += kls_now_seconds() - start;
  if (metis_numeric == NULL || metis_common.status < 0 ||
      metis_common.status == TRILINOS_KLU_SINGULAR) {
    if (metis_numeric != NULL) {
      trilinos_klu_l_free_numeric(&metis_numeric, &metis_common);
    }
    trilinos_klu_l_free_symbolic(&metis_symbolic, &metis_common);
    return;
  }

  (void)trilinos_klu_l_flops(metis_symbolic, metis_numeric, &metis_common);
  if (!metis_numeric_is_better(solver, &metis_common, metis_numeric)) {
    trilinos_klu_l_free_numeric(&metis_numeric, &metis_common);
    trilinos_klu_l_free_symbolic(&metis_symbolic, &metis_common);
    return;
  }

  trilinos_klu_l_symbolic *old_symbolic = solver->symbolic;
  trilinos_klu_l_numeric *old_numeric = solver->numeric;
  trilinos_klu_l_common old_common = solver->common;

  solver->symbolic = metis_symbolic;
  solver->numeric = metis_numeric;
  solver->common = metis_common;
  solver->stats.selected_ordering = KLS_ORDERING_METIS;
  solver->stats.nblocks = (int64_t)solver->symbolic->nblocks;
  solver->stats.max_block = (int64_t)solver->symbolic->maxblock;
  solver->stats.structural_rank = (int64_t)solver->symbolic->structural_rank;
  solver->stats.estimated_flops = solver->symbolic->est_flops;

  trilinos_klu_l_free_numeric(&old_numeric, &old_common);
  trilinos_klu_l_free_symbolic(&old_symbolic, &old_common);
}
#endif

static int choose_symbolic_for_pattern(UF_long n,
                                       UF_long *col_ptr,
                                       UF_long *row_idx,
                                       const kls_options *options,
                                       trilinos_klu_l_symbolic **symbolic_out,
                                       trilinos_klu_l_common *common_out,
                                       kls_ordering *selected_ordering_out,
                                       double *score_out) {
  if (options->ordering != KLS_ORDERING_AUTO) {
    int status = analyze_with_ordering(n, col_ptr, row_idx, options,
                                       options->ordering, symbolic_out,
                                       common_out);
    if (status == KLS_OK) {
      *selected_ordering_out = options->ordering;
      *score_out = symbolic_score(*symbolic_out);
    }
    return status;
  }

  kls_options auto_options = *options;
  const kls_options *symbolic_options = options;
  if (should_start_auto_without_btf(n, col_ptr, row_idx, options)) {
    auto_options.use_btf = 0;
    symbolic_options = &auto_options;
  }

#ifdef KLS_HAVE_METIS
  if (should_start_auto_with_metis(n, col_ptr, row_idx)) {
    int status = analyze_with_ordering(n, col_ptr, row_idx, symbolic_options,
                                       KLS_ORDERING_METIS, symbolic_out,
                                       common_out);
    if (status == KLS_OK) {
      *selected_ordering_out = KLS_ORDERING_METIS;
      *score_out = symbolic_score(*symbolic_out);
      return KLS_OK;
    }
  }
#endif

  const kls_ordering candidates[] = {KLS_ORDERING_AMD, KLS_ORDERING_COLAMD};
  trilinos_klu_l_symbolic *best_symbolic = NULL;
  trilinos_klu_l_common best_common;
  kls_ordering best_ordering = KLS_ORDERING_AMD;
  double best_score = 0.0;
  int any_ok = 0;

  for (size_t i = 0; i < sizeof(candidates) / sizeof(candidates[0]); ++i) {
    trilinos_klu_l_symbolic *candidate_symbolic = NULL;
    trilinos_klu_l_common candidate_common;
    int status = analyze_with_ordering(n, col_ptr, row_idx, symbolic_options, candidates[i],
                                       &candidate_symbolic, &candidate_common);
    if (status != KLS_OK) {
      continue;
    }
    const double score = symbolic_score(candidate_symbolic);
    double selected_score = score;
    maybe_retry_single_block_without_btf(n, col_ptr, row_idx, symbolic_options,
                                         candidates[i], &candidate_symbolic,
                                         &candidate_common, &selected_score);
    if (!any_ok || selected_score < best_score) {
      if (best_symbolic != NULL) {
        trilinos_klu_l_free_symbolic(&best_symbolic, &best_common);
      }
      best_symbolic = candidate_symbolic;
      best_common = candidate_common;
      best_ordering = candidates[i];
      best_score = selected_score;
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
  *symbolic_out = best_symbolic;
  *common_out = best_common;
  *selected_ordering_out = best_ordering;
  *score_out = best_score;
  return KLS_OK;
}

static int validate_index_base(int index_base) {
  return index_base == 0 || index_base == 1;
}

static int validate_n(int64_t n) {
  return n > 0 && n < INT64_MAX / 2 && n <= (int64_t)(UF_long_max / 2);
}

static int validate_options(const kls_options *options) {
  if (options->threads <= 0) {
    return 0;
  }
  if (options->ordering < KLS_ORDERING_AUTO || options->ordering > KLS_ORDERING_METIS) {
    return 0;
  }
  if (options->orientation < KLS_ORIENTATION_AUTO ||
      options->orientation > KLS_ORIENTATION_TRANSPOSE) {
    return 0;
  }
  if (options->scale != KLS_SCALE_AUTO && (options->scale < -1 || options->scale > 2)) {
    return 0;
  }
  if (options->pivot_tolerance < 0.0) {
    return 0;
  }
  if (options->memory_growth <= 0.0) {
    return 0;
  }
  if (options->fast_factor != 0 && options->fast_factor != 1) {
    return 0;
  }
  if (options->static_pivoting != 0 && options->static_pivoting != 1) {
    return 0;
  }
  return 1;
}

static int allocate_candidate(kls_pattern_candidate *candidate,
                              int64_t n,
                              int64_t nnz,
                              kls_orientation orientation,
                              int need_map) {
  memset(candidate, 0, sizeof(*candidate));
  candidate->col_ptr = (UF_long *)calloc((size_t)n + 1u, sizeof(UF_long));
  candidate->row_idx = (UF_long *)calloc((size_t)nnz, sizeof(UF_long));
  if (need_map) {
    candidate->input_to_csc = (UF_long *)calloc((size_t)nnz, sizeof(UF_long));
    candidate->values = (double *)calloc((size_t)nnz, sizeof(double));
  }
  if (candidate->col_ptr == NULL || candidate->row_idx == NULL ||
      (need_map && (candidate->values == NULL || candidate->input_to_csc == NULL))) {
    free(candidate->col_ptr);
    free(candidate->row_idx);
    free(candidate->input_to_csc);
    free(candidate->values);
    memset(candidate, 0, sizeof(*candidate));
    return KLS_ERR_OUT_OF_MEMORY;
  }
  candidate->n = (UF_long)n;
  candidate->nnz = (UF_long)nnz;
  candidate->orientation = orientation;
  candidate->selected_ordering = KLS_ORDERING_AUTO;
  candidate->score = DBL_MAX;
  return KLS_OK;
}

static void free_candidate(kls_pattern_candidate *candidate) {
  if (candidate == NULL) {
    return;
  }
  if (candidate->symbolic != NULL) {
    trilinos_klu_l_free_symbolic(&candidate->symbolic, &candidate->common);
  }
  free(candidate->col_ptr);
  free(candidate->row_idx);
  free(candidate->input_to_csc);
  free(candidate->values);
  memset(candidate, 0, sizeof(*candidate));
}

static int copy_compressed_candidate(kls_pattern_candidate *candidate,
                                     kls_index_type index_type,
                                     int64_t n,
                                     const void *col_ptr,
                                     const void *row_idx,
                                     int index_base,
                                     kls_orientation orientation) {
  const int64_t base = index_base;
  const int64_t nnz = read_index(index_type, col_ptr, n) - base;
  if (nnz < 0) {
    return KLS_ERR_INVALID_ARGUMENT;
  }

  int status = allocate_candidate(candidate, n, nnz, orientation, 0);
  if (status != KLS_OK) {
    return status;
  }

  int64_t prev = 0;
  for (int64_t j = 0; j <= n; ++j) {
    const int64_t ptr = read_index(index_type, col_ptr, j) - base;
    if (ptr < prev || ptr < 0 || ptr > nnz) {
      free_candidate(candidate);
      return KLS_ERR_INVALID_ARGUMENT;
    }
    candidate->col_ptr[j] = (UF_long)ptr;
    prev = ptr;
  }
  for (int64_t p = 0; p < nnz; ++p) {
    const int64_t row = read_index(index_type, row_idx, p) - base;
    if (row < 0 || row >= n) {
      free_candidate(candidate);
      return KLS_ERR_INVALID_ARGUMENT;
    }
    candidate->row_idx[p] = (UF_long)row;
  }
  return KLS_OK;
}

static int transpose_candidate(const kls_pattern_candidate *source,
                               kls_orientation orientation,
                               kls_pattern_candidate *candidate) {
  int status = allocate_candidate(candidate, (int64_t)source->n,
                                  (int64_t)source->nnz, orientation, 1);
  if (status != KLS_OK) {
    return status;
  }

  for (UF_long p = 0; p < source->nnz; ++p) {
    candidate->col_ptr[source->row_idx[p] + 1]++;
  }
  for (UF_long j = 0; j < source->n; ++j) {
    candidate->col_ptr[j + 1] += candidate->col_ptr[j];
  }

  UF_long *next = (UF_long *)malloc(((size_t)source->n) * sizeof(UF_long));
  if (next == NULL) {
    free_candidate(candidate);
    return KLS_ERR_OUT_OF_MEMORY;
  }
  memcpy(next, candidate->col_ptr, ((size_t)source->n) * sizeof(UF_long));

  for (UF_long col = 0; col < source->n; ++col) {
    for (UF_long p = source->col_ptr[col]; p < source->col_ptr[col + 1]; ++p) {
      const UF_long dst_col = source->row_idx[p];
      const UF_long dst = next[dst_col]++;
      candidate->row_idx[dst] = col;
      candidate->input_to_csc[p] = dst;
    }
  }
  free(next);
  return KLS_OK;
}

static int analyze_candidate(kls_pattern_candidate *candidate,
                             const kls_options *options) {
  return choose_symbolic_for_pattern(candidate->n, candidate->col_ptr,
                                     candidate->row_idx, options,
                                     &candidate->symbolic, &candidate->common,
                                     &candidate->selected_ordering,
                                     &candidate->score);
}

static int auto_orientation_prefers_transpose(UF_long n) {
  /* On small and medium SPICE-like matrices, the second symbolic analysis
     usually costs more than the normal-vs-transpose fill estimate saves. */
  return n <= 30000;
}

static int auto_orientation_prefers_normal(UF_long n,
                                           const UF_long *col_ptr,
                                           const UF_long *row_idx) {
  return is_large_diagonal_circuit_like_pattern(n, col_ptr, row_idx);
}

static int select_candidate(kls_pattern_candidate *normal,
                            kls_pattern_candidate *transpose,
                            const kls_options *options,
                            kls_pattern_candidate **chosen_out) {
  *chosen_out = NULL;
  if (options->orientation == KLS_ORIENTATION_NORMAL) {
    if (normal == NULL) {
      return KLS_ERR_INVALID_ARGUMENT;
    }
    int status = analyze_candidate(normal, options);
    if (status == KLS_OK) {
      *chosen_out = normal;
    }
    return status;
  }
  if (options->orientation == KLS_ORIENTATION_TRANSPOSE) {
    if (transpose == NULL) {
      return KLS_ERR_INVALID_ARGUMENT;
    }
    int status = analyze_candidate(transpose, options);
    if (status == KLS_OK) {
      *chosen_out = transpose;
    }
    return status;
  }

  if (transpose != NULL && auto_orientation_prefers_transpose(transpose->n)) {
    const int transpose_status = analyze_candidate(transpose, options);
    if (transpose_status == KLS_OK) {
      *chosen_out = transpose;
      return KLS_OK;
    }
    if (normal != NULL) {
      const int normal_status = analyze_candidate(normal, options);
      if (normal_status == KLS_OK) {
        *chosen_out = normal;
        return KLS_OK;
      }
    }
    return transpose_status;
  }
  if (normal == NULL || transpose == NULL) {
    kls_pattern_candidate *candidate = normal != NULL ? normal : transpose;
    if (candidate == NULL) {
      return KLS_ERR_INVALID_ARGUMENT;
    }
    int status = analyze_candidate(candidate, options);
    if (status == KLS_OK) {
      *chosen_out = candidate;
    }
    return status;
  }
  const int normal_status = analyze_candidate(normal, options);
  const int transpose_status = analyze_candidate(transpose, options);
  if (normal_status == KLS_OK && transpose_status == KLS_OK) {
    *chosen_out = (transpose->score < normal->score) ? transpose : normal;
    return KLS_OK;
  }
  if (normal_status == KLS_OK) {
    *chosen_out = normal;
    return KLS_OK;
  }
  if (transpose_status == KLS_OK) {
    *chosen_out = transpose;
    return KLS_OK;
  }
  return normal_status != KLS_OK ? normal_status : transpose_status;
}

static void adopt_candidate(kls_solver *solver, kls_pattern_candidate *candidate) {
  solver->n = candidate->n;
  solver->nnz = candidate->nnz;
  solver->col_ptr = candidate->col_ptr;
  solver->row_idx = candidate->row_idx;
  solver->input_to_csc = candidate->input_to_csc;
  solver->values = candidate->values;
  solver->orientation = candidate->orientation;
  solver->common = candidate->common;
  solver->symbolic = candidate->symbolic;
  solver->stats.selected_ordering = candidate->selected_ordering;
  solver->stats.selected_orientation = candidate->orientation;

  candidate->col_ptr = NULL;
  candidate->row_idx = NULL;
  candidate->input_to_csc = NULL;
  candidate->values = NULL;
  candidate->symbolic = NULL;
}

static void fill_symbolic_stats(kls_solver *solver, double elapsed) {
  solver->stats.struct_size = sizeof(solver->stats);
  solver->stats.n = (int64_t)solver->n;
  solver->stats.nnz = (int64_t)solver->nnz;
  solver->stats.analysis_seconds = elapsed;
  solver->stats.selected_orientation = solver->orientation;
  solver->stats.selected_scale = (int)solver->common.scale;
  solver->stats.selected_pivot_tolerance = solver->common.tol;
  solver->stats.selected_static_pivoting = solver->row_perm != NULL;
  if (solver->symbolic != NULL) {
    solver->stats.last_kernel_status = (int)solver->common.status;
    solver->stats.selected_btf = solver->symbolic->do_btf ? 1 : 0;
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
  solver->stats.selected_scale = (int)solver->common.scale;
  solver->stats.selected_pivot_tolerance = solver->common.tol;
  solver->stats.selected_static_pivoting = solver->row_perm != NULL;
  solver->stats.selected_btf =
    (solver->symbolic != NULL && solver->symbolic->do_btf) ? 1 : 0;
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
  options->orientation = KLS_ORIENTATION_AUTO;
  options->use_btf = 1;
  options->scale = KLS_SCALE_AUTO;
  options->pivot_tolerance = 0.001;
  options->memory_growth = 1.5;
  options->halt_if_singular = 1;
  options->fast_factor = 1;
  options->static_pivoting = 1;
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

  clear_matrix(solver);
  solver->options = normalized;
  solver->input_format = KLS_INPUT_CSC;
  const double start = kls_now_seconds();

  kls_pattern_candidate normal = {0};
  kls_pattern_candidate transpose = {0};
  kls_pattern_candidate *chosen = NULL;

  int status = copy_compressed_candidate(&normal, index_type, n, col_ptr, row_idx,
                                         index_base, KLS_ORIENTATION_NORMAL);
  if (status != KLS_OK) {
    clear_matrix(solver);
    return status;
  }

  const int prefer_auto_normal =
    normalized.orientation == KLS_ORIENTATION_AUTO &&
    auto_orientation_prefers_normal(normal.n, normal.col_ptr, normal.row_idx);
  if (normalized.orientation != KLS_ORIENTATION_NORMAL && !prefer_auto_normal) {
    status = transpose_candidate(&normal, KLS_ORIENTATION_TRANSPOSE, &transpose);
    if (status != KLS_OK) {
      free_candidate(&normal);
      clear_matrix(solver);
      return status;
    }
  }

  status = select_candidate(&normal,
                            (normalized.orientation == KLS_ORIENTATION_NORMAL ||
                             prefer_auto_normal) ? NULL : &transpose,
                            &normalized, &chosen);
  const double elapsed = kls_now_seconds() - start;
  if (status != KLS_OK) {
    free_candidate(&transpose);
    free_candidate(&normal);
    clear_matrix(solver);
    return status;
  }

  adopt_candidate(solver, chosen);
  fill_symbolic_stats(solver, elapsed);
  free_candidate(&transpose);
  free_candidate(&normal);
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

  clear_matrix(solver);
  solver->options = normalized;
  solver->input_format = KLS_INPUT_CSR;
  const double start = kls_now_seconds();

  kls_pattern_candidate transpose = {0};
  kls_pattern_candidate normal = {0};
  kls_pattern_candidate *chosen = NULL;

  int status = copy_compressed_candidate(&transpose, index_type, n, row_ptr, col_idx,
                                         index_base, KLS_ORIENTATION_TRANSPOSE);
  if (status != KLS_OK) {
    clear_matrix(solver);
    return status;
  }

  const int prefer_auto_transpose =
    normalized.orientation == KLS_ORIENTATION_AUTO &&
    auto_orientation_prefers_transpose((UF_long)n);
  if (normalized.orientation != KLS_ORIENTATION_TRANSPOSE && !prefer_auto_transpose) {
    status = transpose_candidate(&transpose, KLS_ORIENTATION_NORMAL, &normal);
    if (status != KLS_OK) {
      free_candidate(&transpose);
      clear_matrix(solver);
      return status;
    }
  }

  const int prefer_auto_normal =
    normalized.orientation == KLS_ORIENTATION_AUTO && normal.col_ptr != NULL &&
    auto_orientation_prefers_normal(normal.n, normal.col_ptr, normal.row_idx);
  status = select_candidate(normal.col_ptr == NULL ? NULL : &normal,
                            prefer_auto_normal ? NULL : &transpose,
                            &normalized, &chosen);
  const double elapsed = kls_now_seconds() - start;
  if (status != KLS_OK) {
    free_candidate(&normal);
    free_candidate(&transpose);
    clear_matrix(solver);
    return status;
  }

  adopt_candidate(solver, chosen);
  fill_symbolic_stats(solver, elapsed);
  free_candidate(&normal);
  free_candidate(&transpose);
  return KLS_OK;
}

static int prepare_numeric_values(kls_solver *solver, const double *values, double **values_out) {
  if (solver == NULL || values == NULL || values_out == NULL) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  if (solver->input_format == KLS_INPUT_NONE) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  if (solver->input_to_csc == NULL) {
    *values_out = (double *)values;
    return KLS_OK;
  }
  if (solver->values == NULL) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  for (int64_t p = 0; p < solver->nnz; ++p) {
    solver->values[solver->input_to_csc[p]] = values[p];
  }
  *values_out = solver->values;
  return KLS_OK;
}

static int kls_parallel_refactor_is_eligible(const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL) {
    return 0;
  }
  if (solver->options.threads <= 1 || solver->symbolic->nblocks < 8u) {
    return 0;
  }
  if (solver->n < 5000u || solver->symbolic->maxblock == solver->n) {
    return 0;
  }
  if (solver->common.scale > 0 || solver->common.flops < 2.0e7) {
    return 0;
  }
  if (solver->symbolic->nblocks < 64u ||
      solver->symbolic->maxblock * 4u < solver->n * 3u) {
    return 0;
  }
  if (solver->numeric->Offp == NULL || solver->numeric->Offx == NULL ||
      solver->numeric->Pinv == NULL || solver->numeric->Udiag == NULL) {
    return 0;
  }
  if (solver->numeric->Rs != NULL) {
    return 0;
  }
  return 1;
}

static UF_long kls_parallel_refactor(kls_solver *solver, double *numeric_values) {
  if (!kls_parallel_refactor_is_eligible(solver)) {
    return trilinos_klu_l_refactor(solver->col_ptr, solver->row_idx,
                                   numeric_values, solver->symbolic,
                                   solver->numeric, &solver->common);
  }

  trilinos_klu_l_common *common = &solver->common;
  trilinos_klu_l_symbolic *symbolic = solver->symbolic;

  common->status = TRILINOS_KLU_OK;
  common->numerical_rank = KLS_KLU_EMPTY;
  common->singular_col = KLS_KLU_EMPTY;
  common->nrealloc = 0;

  int thread_count = solver->options.threads;
  if ((UF_long)thread_count > symbolic->nblocks) {
    thread_count = (int)symbolic->nblocks;
  }
  if (thread_count < 2) {
    return trilinos_klu_l_refactor(solver->col_ptr, solver->row_idx,
                                   numeric_values, solver->symbolic,
                                   solver->numeric, &solver->common);
  }

  pthread_t *threads = (pthread_t *)calloc((size_t)thread_count, sizeof(*threads));
  kls_parallel_refactor_worker *workers =
    (kls_parallel_refactor_worker *)calloc((size_t)thread_count, sizeof(*workers));
  if (threads == NULL || workers == NULL) {
    free(threads);
    free(workers);
    common->status = TRILINOS_KLU_OUT_OF_MEMORY;
    return 0;
  }

  kls_parallel_refactor_shared shared;
  memset(&shared, 0, sizeof(shared));
  shared.col_ptr = solver->col_ptr;
  shared.row_idx = solver->row_idx;
  shared.values = numeric_values;
  shared.symbolic = symbolic;
  shared.numeric = solver->numeric;
  shared.halt_if_singular = common->halt_if_singular;
  if (pthread_mutex_init(&shared.lock, NULL) != 0) {
    free(threads);
    free(workers);
    common->status = TRILINOS_KLU_OUT_OF_MEMORY;
    return 0;
  }

  int created = 0;
  for (int i = 0; i < thread_count; ++i) {
    workers[i].shared = &shared;
    workers[i].x = (double *)calloc((size_t)symbolic->maxblock, sizeof(double));
    if (workers[i].x == NULL) {
      pthread_mutex_lock(&shared.lock);
      shared.stop = 1;
      pthread_mutex_unlock(&shared.lock);
      break;
    }
    if (pthread_create(&threads[i], NULL, kls_parallel_refactor_worker_main,
                       &workers[i]) != 0) {
      pthread_mutex_lock(&shared.lock);
      shared.stop = 1;
      pthread_mutex_unlock(&shared.lock);
      break;
    }
    created++;
  }

  for (int i = 0; i < created; ++i) {
    pthread_join(threads[i], NULL);
  }
  pthread_mutex_destroy(&shared.lock);

  const int launch_failed = created != thread_count;
  int invalid = 0;
  int singular = 0;
  UF_long numerical_rank = UF_long_max;
  UF_long singular_col = KLS_KLU_EMPTY;
  for (int i = 0; i < thread_count; ++i) {
    if (workers[i].invalid) {
      invalid = 1;
    }
    if (workers[i].singular && workers[i].numerical_rank < numerical_rank) {
      singular = 1;
      numerical_rank = workers[i].numerical_rank;
      singular_col = workers[i].singular_col;
    }
    free(workers[i].x);
  }
  free(threads);
  free(workers);

  if (launch_failed || invalid) {
    common->status = launch_failed ? TRILINOS_KLU_OUT_OF_MEMORY
                                   : TRILINOS_KLU_INVALID;
    return 0;
  }

  if (singular) {
    common->status = TRILINOS_KLU_SINGULAR;
    common->numerical_rank = numerical_rank;
    common->singular_col = singular_col;
    return common->halt_if_singular ? 0 : 1;
  }

  common->status = TRILINOS_KLU_OK;
  return 1;
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

  double elapsed = 0.0;
  const int had_numeric = solver->numeric != NULL;
  if (solver->options.fast_factor && solver->numeric != NULL) {
    const double start = kls_now_seconds();
    const UF_long ok = kls_parallel_refactor(solver, numeric_values);
    elapsed += kls_now_seconds() - start;
    if (ok && solver->common.status >= 0 &&
        solver->common.status != TRILINOS_KLU_SINGULAR) {
      solver->stats.factor_seconds = elapsed;
      (void)trilinos_klu_l_flops(solver->symbolic, solver->numeric, &solver->common);
      (void)trilinos_klu_l_rcond(solver->symbolic, solver->numeric, &solver->common);
      fill_numeric_stats(solver);
      return KLS_OK;
    }
  }

  free_numeric(solver);
  if (!had_numeric) {
    solver->common.scale = choose_auto_scale_from_values(solver, numeric_values);
  }
  const double start = kls_now_seconds();
  solver->numeric = trilinos_klu_l_factor(solver->col_ptr, solver->row_idx,
                                          numeric_values, solver->symbolic,
                                          &solver->common);
  elapsed += kls_now_seconds() - start;
  solver->stats.factor_seconds = elapsed;
  if (solver->numeric == NULL || solver->common.status < 0) {
    fill_numeric_stats(solver);
    return solver->common.status == TRILINOS_KLU_SINGULAR ? KLS_ERR_SINGULAR
                                                          : KLS_ERR_FACTOR_FAILED;
  }
  (void)trilinos_klu_l_flops(solver->symbolic, solver->numeric, &solver->common);
  (void)trilinos_klu_l_rcond(solver->symbolic, solver->numeric, &solver->common);
  maybe_select_auto_row_match(solver, &elapsed, numeric_values);
  maybe_select_auto_scale(solver, &elapsed, numeric_values);
#ifdef KLS_HAVE_METIS
  maybe_promote_auto_metis(solver, &elapsed, numeric_values);
#endif
  (void)trilinos_klu_l_flops(solver->symbolic, solver->numeric, &solver->common);
  (void)trilinos_klu_l_rcond(solver->symbolic, solver->numeric, &solver->common);
  maybe_select_auto_pivot_tolerance(solver, &elapsed, numeric_values);
  (void)trilinos_klu_l_flops(solver->symbolic, solver->numeric, &solver->common);
  (void)trilinos_klu_l_rcond(solver->symbolic, solver->numeric, &solver->common);
  solver->stats.factor_seconds = elapsed;
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
  const UF_long ok = kls_parallel_refactor(solver, numeric_values);
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
  if (nrhs > (int64_t)UF_long_max || ldx > (int64_t)UF_long_max) {
    return KLS_ERR_INVALID_ARGUMENT;
  }

  const double start = kls_now_seconds();
  const int kernel_transpose =
    (solver->orientation == KLS_ORIENTATION_TRANSPOSE) ? !transpose : transpose;
  double *perm_workspace = NULL;
  if (solver->row_perm != NULL) {
    perm_workspace = (double *)malloc((size_t)solver->n * sizeof(*perm_workspace));
    if (perm_workspace == NULL) {
      return KLS_ERR_OUT_OF_MEMORY;
    }
  }

  if (solver->row_perm != NULL && !kernel_transpose) {
    for (int64_t rhs = 0; rhs < nrhs; ++rhs) {
      const double *src = b + rhs * ldb;
      double *dst = x + rhs * ldx;
      for (UF_long row = 0; row < solver->n; ++row) {
        perm_workspace[solver->row_perm[row]] = src[row];
      }
      memcpy(dst, perm_workspace, (size_t)solver->n * sizeof(double));
    }
  } else if (b != x || ldb != ldx) {
    for (int64_t rhs = 0; rhs < nrhs; ++rhs) {
      memmove(x + rhs * ldx, b + rhs * ldb, (size_t)solver->n * sizeof(double));
    }
  }

  const UF_long ok = kernel_transpose
    ? trilinos_klu_l_tsolve(solver->symbolic, solver->numeric, (UF_long)ldx,
                            (UF_long)nrhs, x, &solver->common)
    : trilinos_klu_l_solve(solver->symbolic, solver->numeric, (UF_long)ldx,
                           (UF_long)nrhs, x, &solver->common);
  if (ok && solver->row_perm != NULL && kernel_transpose) {
    for (int64_t rhs = 0; rhs < nrhs; ++rhs) {
      double *dst = x + rhs * ldx;
      for (UF_long row = 0; row < solver->n; ++row) {
        perm_workspace[row] = dst[solver->row_perm[row]];
      }
      memcpy(dst, perm_workspace, (size_t)solver->n * sizeof(double));
    }
  }
  free(perm_workspace);
  solver->stats.solve_seconds = kls_now_seconds() - start;
  fill_numeric_stats(solver);

  if (!ok || solver->common.status < 0) {
    return KLS_ERR_SOLVE_FAILED;
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
  if (solver == NULL || stats == NULL || stats->struct_size == 0) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  const size_t requested_size = stats->struct_size;
  const size_t copy_size = requested_size < sizeof(kls_stats)
    ? requested_size
    : sizeof(kls_stats);
  memcpy(stats, &solver->stats, copy_size);
  stats->struct_size = sizeof(kls_stats);
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
    case KLS_ORDERING_METIS: return "metis";
    default: return "unknown";
  }
}

const char *kls_orientation_name(kls_orientation orientation) {
  switch (orientation) {
    case KLS_ORIENTATION_AUTO: return "auto";
    case KLS_ORIENTATION_NORMAL: return "normal";
    case KLS_ORIENTATION_TRANSPOSE: return "transpose";
    default: return "unknown";
  }
}
