#define _POSIX_C_SOURCE 200809L

#include "kls/kls.h"

#include "trilinos_klu_decl.h"
#ifdef KLS_HAVE_METIS
#include "metis.h"
#endif

#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
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
  kls_input_format input_format;
  kls_orientation orientation;
  kls_options options;
  kls_stats stats;
  trilinos_klu_l_common common;
  trilinos_klu_l_symbolic *symbolic;
  trilinos_klu_l_numeric *numeric;
  int auto_metis_checked;
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
  solver->col_ptr = NULL;
  solver->row_idx = NULL;
  solver->input_to_csc = NULL;
  solver->values = NULL;
  solver->n = 0;
  solver->nnz = 0;
  solver->input_format = KLS_INPUT_NONE;
  solver->orientation = KLS_ORIENTATION_NORMAL;
  solver->auto_metis_checked = 0;
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

static int choose_auto_scale_from_values(const kls_solver *solver,
                                         const double *numeric_values) {
  if (solver == NULL || numeric_values == NULL ||
      solver->options.scale != KLS_SCALE_AUTO || solver->n <= 0) {
    return initial_scale(&solver->options);
  }

  double *row_max = (double *)calloc((size_t)solver->n, sizeof(*row_max));
  if (row_max == NULL) {
    return initial_scale(&solver->options);
  }

  /* Cheap first-factor policy:
     - keep well-populated, moderately scaled circuit diagonals unscaled;
     - use sum scaling for sparse diagonals whose row magnitudes are balanced. */
  const double diag_spread_limit = 1.0e7;
  const double min_diag_fraction = 0.90;
  const double row_spread_limit = 1.0e4;
  const double row_p90_p10_limit = 5.5;
  UF_long diag_count = 0;
  double min_diag = DBL_MAX;
  double max_diag = 0.0;

  for (UF_long col = 0; col < solver->n; ++col) {
    double diag_abs = 0.0;
    for (UF_long p = solver->col_ptr[col]; p < solver->col_ptr[col + 1]; ++p) {
      if (solver->row_idx[p] == col) {
        const double value_abs = fabs(numeric_values[p]);
        if (isfinite(value_abs) && value_abs > diag_abs) {
          diag_abs = value_abs;
        }
      }
      const double value_abs = fabs(numeric_values[p]);
      if (isfinite(value_abs) && value_abs > row_max[(size_t)solver->row_idx[p]]) {
        row_max[(size_t)solver->row_idx[p]] = value_abs;
      }
    }
    if (diag_abs > 0.0) {
      diag_count++;
      if (diag_abs < min_diag) {
        min_diag = diag_abs;
      }
      if (diag_abs > max_diag) {
        max_diag = diag_abs;
      }
    }
  }

  if ((double)diag_count >= min_diag_fraction * (double)solver->n &&
      min_diag > 0.0 && max_diag > 0.0 && max_diag / min_diag <= diag_spread_limit) {
    free(row_max);
    return 0;
  }

  size_t row_count = 0;
  for (UF_long row = 0; row < solver->n; ++row) {
    if (row_max[(size_t)row] > 0.0) {
      row_max[row_count++] = row_max[(size_t)row];
    }
  }

  if ((double)diag_count < min_diag_fraction * (double)solver->n &&
      (double)row_count >= min_diag_fraction * (double)solver->n && row_count > 1u) {
    qsort(row_max, row_count, sizeof(*row_max), compare_double);
    const double row_min = row_max[0];
    const double row_maximum = row_max[row_count - 1u];
    const size_t p10_index = (row_count - 1u) / 10u;
    const size_t p90_index = ((row_count - 1u) * 9u) / 10u;
    const double row_p10 = row_max[p10_index];
    const double row_p90 = row_max[p90_index];
    if (row_min > 0.0 && row_p10 > 0.0 &&
        row_maximum / row_min <= row_spread_limit &&
        row_p90 / row_p10 <= row_p90_p10_limit) {
      free(row_max);
      return 1;
    }
  }

  free(row_max);
  return initial_scale(&solver->options);
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

  idx_t nvtxs = (idx_t)n;
  const int metis_status = METIS_NodeND(&nvtxs, xadj, adjncy, NULL, options,
                                        metis_perm, metis_iperm);
  if (metis_status == METIS_OK) {
    for (UF_long i = 0; i < n; ++i) {
      perm_out[i] = (UF_long)metis_perm[i];
    }
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
  return metis_status == METIS_OK ? (UF_long)-1 : 0;
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

static int numeric_candidate_is_better(const trilinos_klu_l_common *current_common,
                                       const trilinos_klu_l_numeric *current_numeric,
                                       const trilinos_klu_l_common *candidate_common,
                                       const trilinos_klu_l_numeric *candidate_numeric) {
  const double current_flops = current_common->flops;
  const double candidate_flops = candidate_common->flops;
  const double current_fill = (double)(current_numeric->lnz + current_numeric->unz);
  const double candidate_fill = (double)(candidate_numeric->lnz + candidate_numeric->unz);

  if (current_flops > 0.0 && candidate_flops > 0.0 &&
      candidate_flops < 0.75 * current_flops &&
      candidate_fill <= 1.10 * current_fill) {
    return 1;
  }
  return candidate_fill < 0.80 * current_fill;
}

static int should_try_auto_scale(const kls_solver *solver) {
  if (solver->auto_scale_checked || solver->options.scale != KLS_SCALE_AUTO ||
      solver->numeric == NULL || solver->n < 20000) {
    return 0;
  }

  const double flops = solver->common.flops;
  const UF_long fill = solver->numeric->lnz + solver->numeric->unz;
  return flops >= 5.0e6 && fill >= 300000;
}

static void maybe_select_auto_scale(kls_solver *solver,
                                    double *elapsed,
                                    const double *numeric_values) {
  if (!should_try_auto_scale(solver)) {
    return;
  }
  solver->auto_scale_checked = 1;

  const int candidates[] = {0, 1, 2};
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

  const kls_ordering candidates[] = {KLS_ORDERING_AMD, KLS_ORDERING_COLAMD};
  trilinos_klu_l_symbolic *best_symbolic = NULL;
  trilinos_klu_l_common best_common;
  kls_ordering best_ordering = KLS_ORDERING_AMD;
  double best_score = 0.0;
  int any_ok = 0;

  for (size_t i = 0; i < sizeof(candidates) / sizeof(candidates[0]); ++i) {
    trilinos_klu_l_symbolic *candidate_symbolic = NULL;
    trilinos_klu_l_common candidate_common;
    int status = analyze_with_ordering(n, col_ptr, row_idx, options, candidates[i],
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

  if (normal == NULL || transpose == NULL) {
    return KLS_ERR_INVALID_ARGUMENT;
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
  solver->stats.selected_scale = (int)solver->common.scale;
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

  if (normalized.orientation != KLS_ORIENTATION_NORMAL) {
    status = transpose_candidate(&normal, KLS_ORIENTATION_TRANSPOSE, &transpose);
    if (status != KLS_OK) {
      free_candidate(&normal);
      clear_matrix(solver);
      return status;
    }
  }

  status = select_candidate(&normal,
                            normalized.orientation == KLS_ORIENTATION_NORMAL ? NULL : &transpose,
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

  if (normalized.orientation != KLS_ORIENTATION_TRANSPOSE) {
    status = transpose_candidate(&transpose, KLS_ORIENTATION_NORMAL, &normal);
    if (status != KLS_OK) {
      free_candidate(&transpose);
      clear_matrix(solver);
      return status;
    }
  }

  status = select_candidate(normal.col_ptr == NULL ? NULL : &normal, &transpose,
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
    const UF_long ok = trilinos_klu_l_refactor(solver->col_ptr, solver->row_idx,
                                               numeric_values, solver->symbolic,
                                               solver->numeric, &solver->common);
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
  maybe_select_auto_scale(solver, &elapsed, numeric_values);
#ifdef KLS_HAVE_METIS
  maybe_promote_auto_metis(solver, &elapsed, numeric_values);
#endif
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
  if (nrhs > (int64_t)UF_long_max || ldx > (int64_t)UF_long_max) {
    return KLS_ERR_INVALID_ARGUMENT;
  }

  const double start = kls_now_seconds();
  if (b != x || ldb != ldx) {
    for (int64_t rhs = 0; rhs < nrhs; ++rhs) {
      memmove(x + rhs * ldx, b + rhs * ldb, (size_t)solver->n * sizeof(double));
    }
  }
  const int kernel_transpose =
    (solver->orientation == KLS_ORIENTATION_TRANSPOSE) ? !transpose : transpose;
  const UF_long ok = kernel_transpose
    ? trilinos_klu_l_tsolve(solver->symbolic, solver->numeric, (UF_long)ldx,
                            (UF_long)nrhs, x, &solver->common)
    : trilinos_klu_l_solve(solver->symbolic, solver->numeric, (UF_long)ldx,
                           (UF_long)nrhs, x, &solver->common);
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
