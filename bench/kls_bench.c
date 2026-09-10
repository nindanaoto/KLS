#define _GNU_SOURCE
#define _POSIX_C_SOURCE 200809L

#include "kls/kls.h"
#include "bench_value_sequence.h"

#include <ctype.h>
#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

typedef struct triplet {
  int64_t row;
  int64_t col;
  double value;
} triplet;

typedef struct matrix {
  int64_t n;
  int64_t nnz;
  int64_t *col_ptr;
  int64_t *row_idx;
  double *values;
} matrix;

typedef enum bench_index_mode {
  BENCH_INDEX_AUTO = 0,
  BENCH_INDEX_INT32,
  BENCH_INDEX_INT64
} bench_index_mode;

typedef struct bench_index_view {
  kls_index_type type;
  int bytes;
  const void *col_ptr;
  const void *row_idx;
  int32_t *col_ptr32;
  int32_t *row_idx32;
} bench_index_view;


static int bench_env_enabled(const char *name) {
  const char *value = getenv(name);
  return value != NULL && value[0] != '\0' && strcmp(value, "0") != 0 &&
         strcasecmp(value, "false") != 0 && strcasecmp(value, "off") != 0 &&
         strcasecmp(value, "no") != 0;
}

static void matrix_free(matrix *a) {
  if (a == NULL) return;
  free(a->col_ptr);
  free(a->row_idx);
  free(a->values);
  memset(a, 0, sizeof(*a));
}

static void bench_index_view_free(bench_index_view *view) {
  if (view == NULL) return;
  free(view->col_ptr32);
  free(view->row_idx32);
  memset(view, 0, sizeof(*view));
}

static int cmp_triplet(const void *lhs, const void *rhs) {
  const triplet *a = (const triplet *)lhs;
  const triplet *b = (const triplet *)rhs;
  if (a->col != b->col) return (a->col < b->col) ? -1 : 1;
  if (a->row != b->row) return (a->row < b->row) ? -1 : 1;
  return 0;
}

static int starts_with_ci(const char *s, const char *prefix) {
  while (*prefix != '\0') {
    if (tolower((unsigned char)*s) != tolower((unsigned char)*prefix)) return 0;
    ++s;
    ++prefix;
  }
  return 1;
}

static int contains_token_ci(const char *s, const char *needle) {
  const size_t n = strlen(needle);
  for (; *s != '\0'; ++s) {
    if (strncasecmp(s, needle, n) == 0) return 1;
  }
  return 0;
}

static int push_triplet(triplet **items, int64_t *count, int64_t *capacity,
                        int64_t row, int64_t col, double value) {
  if (*count == *capacity) {
    const int64_t next_capacity = (*capacity == 0) ? 1024 : (*capacity * 2);
    triplet *next = (triplet *)realloc(*items, (size_t)next_capacity * sizeof(triplet));
    if (next == NULL) return 0;
    *items = next;
    *capacity = next_capacity;
  }
  (*items)[*count].row = row;
  (*items)[*count].col = col;
  (*items)[*count].value = value;
  ++(*count);
  return 1;
}

static int read_matrix_market(const char *path, matrix *out) {
  FILE *fp = fopen(path, "r");
  if (fp == NULL) {
    fprintf(stderr, "cannot open %s: %s\n", path, strerror(errno));
    return 0;
  }

  char line[4096];
  if (fgets(line, sizeof(line), fp) == NULL ||
      !starts_with_ci(line, "%%MatrixMarket matrix coordinate")) {
    fprintf(stderr, "%s is not a coordinate MatrixMarket file\n", path);
    fclose(fp);
    return 0;
  }

  const int is_pattern = contains_token_ci(line, "pattern");
  const int is_complex = contains_token_ci(line, "complex");
  const int is_symmetric = contains_token_ci(line, "symmetric");
  const int is_skew = contains_token_ci(line, "skew-symmetric");
  if (is_complex) {
    fprintf(stderr, "%s uses complex values; this baseline supports real matrices only\n", path);
    fclose(fp);
    return 0;
  }

  do {
    if (fgets(line, sizeof(line), fp) == NULL) {
      fprintf(stderr, "%s has no size line\n", path);
      fclose(fp);
      return 0;
    }
  } while (line[0] == '%');

  int64_t rows = 0, cols = 0, entries = 0;
  if (sscanf(line, "%" SCNd64 " %" SCNd64 " %" SCNd64, &rows, &cols, &entries) != 3 ||
      rows <= 0 || cols <= 0 || rows != cols || entries < 0) {
    fprintf(stderr, "%s must be a square sparse matrix\n", path);
    fclose(fp);
    return 0;
  }

  triplet *items = NULL;
  int64_t count = 0;
  int64_t capacity = 0;
  for (int64_t k = 0; k < entries; ++k) {
    if (fgets(line, sizeof(line), fp) == NULL) {
      fprintf(stderr, "%s ended before all entries were read\n", path);
      free(items);
      fclose(fp);
      return 0;
    }
    int64_t row = 0, col = 0;
    double value = 1.0;
    int fields = 0;
    if (is_pattern) {
      fields = sscanf(line, "%" SCNd64 " %" SCNd64, &row, &col);
    } else {
      fields = sscanf(line, "%" SCNd64 " %" SCNd64 " %lf", &row, &col, &value);
    }
    if (fields < 2 || row <= 0 || row > rows || col <= 0 || col > cols) {
      fprintf(stderr, "%s has an invalid entry near %" PRId64 "\n", path, k + 1);
      free(items);
      fclose(fp);
      return 0;
    }
    --row;
    --col;
    if (!push_triplet(&items, &count, &capacity, row, col, value)) {
      free(items);
      fclose(fp);
      return 0;
    }
    if ((is_symmetric || is_skew) && row != col) {
      const double mirrored = is_skew ? -value : value;
      if (!push_triplet(&items, &count, &capacity, col, row, mirrored)) {
        free(items);
        fclose(fp);
        return 0;
      }
    }
  }
  fclose(fp);

  qsort(items, (size_t)count, sizeof(triplet), cmp_triplet);
  int64_t unique = 0;
  for (int64_t i = 0; i < count; ) {
    int64_t j = i + 1;
    double sum = items[i].value;
    while (j < count && items[j].row == items[i].row && items[j].col == items[i].col) {
      sum += items[j].value;
      ++j;
    }
    if (sum != 0.0) {
      items[unique] = items[i];
      items[unique].value = sum;
      ++unique;
    }
    i = j;
  }

  out->n = rows;
  out->nnz = unique;
  out->col_ptr = (int64_t *)calloc((size_t)rows + 1u, sizeof(int64_t));
  out->row_idx = (int64_t *)calloc((size_t)unique, sizeof(int64_t));
  out->values = (double *)calloc((size_t)unique, sizeof(double));
  if (out->col_ptr == NULL || out->row_idx == NULL || out->values == NULL) {
    free(items);
    matrix_free(out);
    return 0;
  }

  for (int64_t p = 0; p < unique; ++p) {
    out->col_ptr[items[p].col + 1]++;
  }
  for (int64_t col = 0; col < rows; ++col) {
    out->col_ptr[col + 1] += out->col_ptr[col];
  }
  int64_t *next = (int64_t *)malloc((size_t)rows * sizeof(int64_t));
  if (next == NULL) {
    free(items);
    matrix_free(out);
    return 0;
  }
  memcpy(next, out->col_ptr, (size_t)rows * sizeof(int64_t));
  for (int64_t p = 0; p < unique; ++p) {
    const int64_t dst = next[items[p].col]++;
    out->row_idx[dst] = items[p].row;
    out->values[dst] = items[p].value;
  }
  free(next);
  free(items);
  return 1;
}

static void matvec_values(const matrix *a, const double *values,
                          const double *x, double *y) {
  memset(y, 0, (size_t)a->n * sizeof(double));
  for (int64_t col = 0; col < a->n; ++col) {
    const double xj = x[col];
    for (int64_t p = a->col_ptr[col]; p < a->col_ptr[col + 1]; ++p) {
      y[a->row_idx[p]] += values[p] * xj;
    }
  }
}

static void make_refactor_values(const matrix *a, const double *base_values,
                                 bench_refactor_value_mode mode,
                                 uint64_t generation, double amplitude,
                                 double *values_out) {
  for (int64_t col = 0; col < a->n; ++col) {
    for (int64_t p = a->col_ptr[col]; p < a->col_ptr[col + 1]; ++p) {
      values_out[p] = bench_generated_refactor_value(
        mode, base_values[p], (uint64_t)a->row_idx[p], (uint64_t)col,
        (uint64_t)a->n, generation, amplitude);
    }
  }
}

static double residual_norm_values(const matrix *a, const double *values,
                                   const double *x, const double *b,
                                   double *relative_out) {
  double *ax = (double *)calloc((size_t)a->n, sizeof(double));
  if (ax == NULL) return INFINITY;
  matvec_values(a, values, x, ax);
  double r2 = 0.0;
  double b2 = 0.0;
  for (int64_t i = 0; i < a->n; ++i) {
    const double r = ax[i] - b[i];
    r2 += r * r;
    b2 += b[i] * b[i];
  }
  free(ax);
  const double rn = sqrt(r2);
  *relative_out = (b2 > 0.0) ? rn / sqrt(b2) : rn;
  return rn;
}

static kls_ordering parse_ordering(const char *s) {
  if (strcmp(s, "amd") == 0) return KLS_ORDERING_AMD;
  if (strcmp(s, "colamd") == 0) return KLS_ORDERING_COLAMD;
  if (strcmp(s, "natural") == 0) return KLS_ORDERING_NATURAL;
  if (strcmp(s, "metis") == 0) return KLS_ORDERING_METIS;
  if (strcmp(s, "scotch") == 0) return KLS_ORDERING_SCOTCH;
  if (strcmp(s, "amf") == 0) return KLS_ORDERING_AMF;
  if (strcmp(s, "ammf") == 0 || strcmp(s, "amf2") == 0) {
    return KLS_ORDERING_AMMF;
  }
  if (strcmp(s, "amf3") == 0) return KLS_ORDERING_AMF3;
  return KLS_ORDERING_AUTO;
}

static kls_orientation parse_orientation(const char *s) {
  if (strcmp(s, "normal") == 0) return KLS_ORIENTATION_NORMAL;
  if (strcmp(s, "transpose") == 0) return KLS_ORIENTATION_TRANSPOSE;
  return KLS_ORIENTATION_AUTO;
}

static int parse_backend(const char *s, kls_backend *backend_out) {
  if (strcmp(s, "auto") == 0) {
    *backend_out = KLS_BACKEND_AUTO;
    return 1;
  }
  if (strcmp(s, "kls") == 0) {
    *backend_out = KLS_BACKEND_KLS;
    return 1;
  }
  if (strcmp(s, "serial") == 0) {
    *backend_out = KLS_BACKEND_SERIAL;
    return 1;
  }
  return 0;
}

static int parse_scale(const char *s, int *scale_out) {
  if (strcmp(s, "auto") == 0) {
    *scale_out = KLS_SCALE_AUTO;
    return 1;
  }
  char *end = NULL;
  const long value = strtol(s, &end, 10);
  if (end == s || *end != '\0' || value < -1 || value > 2) {
    return 0;
  }
  *scale_out = (int)value;
  return 1;
}

static int parse_nonnegative_double(const char *s, double *value_out) {
  char *end = NULL;
  errno = 0;
  const double value = strtod(s, &end);
  if (errno != 0 || end == s || *end != '\0' || !isfinite(value) || value < 0.0) {
    return 0;
  }
  *value_out = value;
  return 1;
}

static int parse_index_mode(const char *s, bench_index_mode *mode_out) {
  if (strcmp(s, "auto") == 0) {
    *mode_out = BENCH_INDEX_AUTO;
    return 1;
  }
  if (strcmp(s, "32") == 0 || strcmp(s, "int32") == 0) {
    *mode_out = BENCH_INDEX_INT32;
    return 1;
  }
  if (strcmp(s, "64") == 0 || strcmp(s, "int64") == 0) {
    *mode_out = BENCH_INDEX_INT64;
    return 1;
  }
  return 0;
}

static const char *index_mode_name(bench_index_mode mode) {
  switch (mode) {
    case BENCH_INDEX_AUTO: return "auto";
    case BENCH_INDEX_INT32: return "32";
    case BENCH_INDEX_INT64: return "64";
    default: return "unknown";
  }
}

static int matrix_fits_int32_indices(const matrix *a) {
  return a != NULL &&
         a->n >= 0 && a->n <= (int64_t)INT32_MAX &&
         a->nnz >= 0 && a->nnz <= (int64_t)INT32_MAX;
}

static int prepare_index_view(const matrix *a,
                              bench_index_mode mode,
                              bench_index_view *view) {
  if (a == NULL || view == NULL || a->col_ptr == NULL || a->row_idx == NULL) {
    return 0;
  }
  memset(view, 0, sizeof(*view));

  const int fits_int32 = matrix_fits_int32_indices(a);
  const int use_int32 =
    mode == BENCH_INDEX_INT32 || (mode == BENCH_INDEX_AUTO && fits_int32);
  if (!use_int32) {
    if (mode == BENCH_INDEX_INT32) {
      fprintf(stderr,
              "matrix order/nnz exceed 32-bit CSC index range; use --input-index auto or 64\n");
      return 0;
    }
    view->type = KLS_INDEX_INT64;
    view->bytes = 8;
    view->col_ptr = a->col_ptr;
    view->row_idx = a->row_idx;
    return 1;
  }

  if (!fits_int32 ||
      a->n + 1 > (int64_t)(SIZE_MAX / sizeof(*view->col_ptr32)) ||
      a->nnz > (int64_t)(SIZE_MAX / sizeof(*view->row_idx32))) {
    fprintf(stderr,
            "matrix order/nnz exceed 32-bit CSC index range; use --input-index auto or 64\n");
    return 0;
  }
  view->col_ptr32 =
    (int32_t *)malloc((size_t)(a->n + 1) * sizeof(*view->col_ptr32));
  view->row_idx32 = a->nnz > 0
    ? (int32_t *)malloc((size_t)a->nnz * sizeof(*view->row_idx32))
    : NULL;
  if (view->col_ptr32 == NULL || (a->nnz > 0 && view->row_idx32 == NULL)) {
    bench_index_view_free(view);
    return 0;
  }
  for (int64_t i = 0; i <= a->n; ++i) {
    if (a->col_ptr[i] < 0 || a->col_ptr[i] > (int64_t)INT32_MAX) {
      fprintf(stderr,
              "matrix column pointer exceeds 32-bit CSC index range; use --input-index auto or 64\n");
      bench_index_view_free(view);
      return 0;
    }
    view->col_ptr32[i] = (int32_t)a->col_ptr[i];
  }
  for (int64_t p = 0; p < a->nnz; ++p) {
    if (a->row_idx[p] < 0 || a->row_idx[p] > (int64_t)INT32_MAX) {
      fprintf(stderr,
              "matrix row index exceeds 32-bit CSC index range; use --input-index auto or 64\n");
      bench_index_view_free(view);
      return 0;
    }
    view->row_idx32[p] = (int32_t)a->row_idx[p];
  }
  view->type = KLS_INDEX_INT32;
  view->bytes = 4;
  view->col_ptr = view->col_ptr32;
  view->row_idx = view->row_idx32;
  return 1;
}

static int valid_row_refactor_control(const char *s) {
  return strcmp(s, "env") == 0 ||
         strcmp(s, "off") == 0 ||
         strcmp(s, "refactor") == 0 ||
         strcmp(s, "checked") == 0 ||
         strcmp(s, "all") == 0;
}

static int apply_row_refactor_control(const char *s) {
  if (strcmp(s, "env") == 0) {
    return 1;
  }
  const int enable_refactor =
    strcmp(s, "refactor") == 0 || strcmp(s, "all") == 0;
  const int enable_checked =
    strcmp(s, "checked") == 0 || strcmp(s, "all") == 0;
  if (setenv("KLS_ENABLE_ROW_REFACTOR",
             enable_refactor ? "1" : "0", 1) != 0) {
    return 0;
  }
  if (setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR",
             enable_checked ? "1" : "0", 1) != 0) {
    return 0;
  }
  return 1;
}

static const char *scale_name(int scale) {
  switch (scale) {
    case KLS_SCALE_AUTO: return "auto";
    case -1: return "-1";
    case 0: return "0";
    case 1: return "1";
    case 2: return "2";
    default: return "unknown";
  }
}

static void usage(const char *argv0) {
  fprintf(stderr,
          "Usage: %s <matrix.mtx> [--lifecycle-systems N] [--repeat N] [--factor-repeat N] [--refactor-repeat N] [--expected-refactors N] [--expected-solves N] [--refactor-values unchanged|rank-preserving|entrywise|localized-entrywise] [--refactor-value-amplitude A] [--threads N] [--tuning-profile FILE] [--backend auto|kls|serial] [--ordering auto|amd|colamd|natural|metis|scotch|amf|ammf|amf3] [--orientation auto|normal|transpose] [--scale auto|-1|0|1|2] [--input-index auto|32|64] [--pivot-tol T] [--row-refactor env|off|refactor|checked|all] [--no-btf] [--no-fast-factor] [--no-static-pivoting] [--no-transpose-solve] [--analyze-only] [--json]\n",
          argv0);
}

int main(int argc, char **argv) {
  if (argc < 2) {
    usage(argv[0]);
    return EXIT_FAILURE;
  }
  const char *path = argv[1];
  int repeat = 5;
  int factor_repeat = -1;
  int refactor_repeat = 5;
  int transpose_solve = 1;
  int json = 0;
  int analyze_only = 0;
  int lifecycle_systems = 0;
  int expected_refactors_explicit = 0;
  int expected_solves_explicit = 0;
  bench_refactor_value_mode refactor_value_mode =
    BENCH_REFACTOR_VALUES_UNCHANGED;
  double refactor_value_amplitude = 1.0e-3;
  const char *row_refactor_control = "env";
  bench_index_mode input_index_mode = BENCH_INDEX_AUTO;
  kls_options options;
  kls_default_options(&options);
  /* The reported spice_cycle_seconds metric is an H100 lifecycle even when
     fewer samples are requested to estimate its steady terms. */
  options.expected_refactorizations = 99;
  options.expected_solves = 100;

  for (int i = 2; i < argc; ++i) {
    if (strcmp(argv[i], "--json") == 0) {
      json = 1;
    } else if (strcmp(argv[i], "--analyze-only") == 0) {
      analyze_only = 1;
    } else if (strcmp(argv[i], "--lifecycle-systems") == 0 &&
               i + 1 < argc) {
      lifecycle_systems = atoi(argv[++i]);
    } else if (strcmp(argv[i], "--repeat") == 0 && i + 1 < argc) {
      repeat = atoi(argv[++i]);
    } else if (strcmp(argv[i], "--factor-repeat") == 0 && i + 1 < argc) {
      factor_repeat = atoi(argv[++i]);
    } else if (strcmp(argv[i], "--refactor-repeat") == 0 && i + 1 < argc) {
      refactor_repeat = atoi(argv[++i]);
    } else if (strcmp(argv[i], "--expected-refactors") == 0 &&
               i + 1 < argc) {
      options.expected_refactorizations = atoll(argv[++i]);
      expected_refactors_explicit = 1;
    } else if (strcmp(argv[i], "--expected-solves") == 0 &&
               i + 1 < argc) {
      options.expected_solves = atoll(argv[++i]);
      expected_solves_explicit = 1;
    } else if (strcmp(argv[i], "--refactor-values") == 0 &&
               i + 1 < argc) {
      if (!bench_parse_refactor_value_mode(argv[++i], &refactor_value_mode)) {
        usage(argv[0]);
        return EXIT_FAILURE;
      }
    } else if (strcmp(argv[i], "--refactor-value-amplitude") == 0 &&
               i + 1 < argc) {
      if (!parse_nonnegative_double(argv[++i], &refactor_value_amplitude)) {
        usage(argv[0]);
        return EXIT_FAILURE;
      }
    } else if (strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
      options.threads = atoi(argv[++i]);
    } else if (strcmp(argv[i], "--tuning-profile") == 0 && i + 1 < argc) {
      options.tuning_profile_path = argv[++i];
    } else if (strcmp(argv[i], "--backend") == 0 && i + 1 < argc) {
      if (!parse_backend(argv[++i], &options.backend)) {
        usage(argv[0]);
        return EXIT_FAILURE;
      }
    } else if (strcmp(argv[i], "--ordering") == 0 && i + 1 < argc) {
      options.ordering = parse_ordering(argv[++i]);
    } else if (strcmp(argv[i], "--orientation") == 0 && i + 1 < argc) {
      options.orientation = parse_orientation(argv[++i]);
    } else if (strcmp(argv[i], "--scale") == 0 && i + 1 < argc) {
      if (!parse_scale(argv[++i], &options.scale)) {
        usage(argv[0]);
        return EXIT_FAILURE;
      }
    } else if (strcmp(argv[i], "--input-index") == 0 && i + 1 < argc) {
      if (!parse_index_mode(argv[++i], &input_index_mode)) {
        usage(argv[0]);
        return EXIT_FAILURE;
      }
    } else if (strcmp(argv[i], "--pivot-tol") == 0 && i + 1 < argc) {
      if (!parse_nonnegative_double(argv[++i], &options.pivot_tolerance)) {
        usage(argv[0]);
        return EXIT_FAILURE;
      }
    } else if (strcmp(argv[i], "--row-refactor") == 0 && i + 1 < argc) {
      row_refactor_control = argv[++i];
      if (!valid_row_refactor_control(row_refactor_control)) {
        usage(argv[0]);
        return EXIT_FAILURE;
      }
    } else if (strcmp(argv[i], "--no-btf") == 0) {
      options.use_btf = 0;
    } else if (strcmp(argv[i], "--no-fast-factor") == 0) {
      options.fast_factor = 0;
    } else if (strcmp(argv[i], "--no-static-pivoting") == 0) {
      options.static_pivoting = 0;
    } else if (strcmp(argv[i], "--no-transpose-solve") == 0) {
      transpose_solve = 0;
    } else {
      usage(argv[0]);
      return EXIT_FAILURE;
    }
  }
  if (lifecycle_systems > 0) {
    repeat = 1;
    factor_repeat = 0;
    refactor_repeat = lifecycle_systems - 1;
    transpose_solve = 0;
    if (!expected_refactors_explicit) {
      options.expected_refactorizations = lifecycle_systems - 1;
    }
    if (!expected_solves_explicit) {
      options.expected_solves = lifecycle_systems;
    }
  }
  if (factor_repeat < 0) factor_repeat = repeat;
  if (repeat <= 0 || factor_repeat < 0 || refactor_repeat < 0 ||
      lifecycle_systems < 0 ||
      options.threads <= 0 || refactor_value_amplitude >= 1.0 ||
      (refactor_value_mode != BENCH_REFACTOR_VALUES_UNCHANGED &&
       refactor_value_amplitude <= 0.0)) {
    usage(argv[0]);
    return EXIT_FAILURE;
  }
  if (!apply_row_refactor_control(row_refactor_control)) {
    perror("apply row-refactor control");
    return EXIT_FAILURE;
  }

  matrix a = {0};
  if (!read_matrix_market(path, &a)) {
    return EXIT_FAILURE;
  }
  bench_index_view input_index = {0};
  if (!prepare_index_view(&a, input_index_mode, &input_index)) {
    matrix_free(&a);
    return EXIT_FAILURE;
  }

  if (analyze_only) {
    kls_solver *solver = NULL;
    int status = kls_create(&solver);
    if (status == KLS_OK) {
      status = kls_analyze_csc(solver, input_index.type, a.n,
                               input_index.col_ptr, input_index.row_idx, 0,
                               &options);
    }
    if (status != KLS_OK) {
      fprintf(stderr, "KLS analyze failed: %s (%d)\n",
              kls_status_string(status), status);
      kls_destroy(solver);
      bench_index_view_free(&input_index);
      matrix_free(&a);
      return EXIT_FAILURE;
    }
    kls_stats stats;
    stats.struct_size = sizeof(stats);
    kls_get_stats(solver, &stats);
    if (json) {
      printf("{\"matrix\":\"%s\",\"n\":%" PRId64 ",\"nnz\":%" PRId64
             ",\"threads\":%d"
             ",\"backend\":\"%s\""
             ",\"requested_input_index\":\"%s\""
             ",\"input_index_bytes\":%d"
             ",\"internal_index_bytes\":%d"
             ",\"requested_orientation\":\"%s\",\"orientation\":\"%s\""
             ",\"ordering\":\"%s\",\"requested_scale\":\"%s\""
             ",\"row_refactor_control\":\"%s\""
             ",\"requested_btf\":%s,\"btf\":%s"
             ",\"analysis_seconds\":%.9g"
             ",\"nblocks\":%" PRId64 ",\"max_block\":%" PRId64
             ",\"structural_rank\":%" PRId64
             ",\"nnz_l\":%" PRId64 ",\"nnz_u\":%" PRId64
             ",\"estimated_flops\":%.9g"
             ",\"build_has_metis\":%s,\"build_has_scotch\":%s"
             ",\"build_has_spral_scaling\":%s"
             ",\"tuning_profile_active\":%s"
             ",\"tuning_profile_id\":\"%016" PRIx64 "\""
             ",\"tuning_profile_field_count\":%d"
             ",\"analyze_only\":true,\"status\":%d}\n",
             path, a.n, a.nnz, options.threads,
             kls_backend_name(options.backend),
             index_mode_name(input_index_mode),
             input_index.bytes,
             stats.internal_index_bytes,
             kls_orientation_name(options.orientation),
             kls_orientation_name(stats.selected_orientation),
             kls_ordering_name(stats.selected_ordering),
             scale_name(options.scale),
             row_refactor_control,
             options.use_btf ? "true" : "false",
             stats.selected_btf ? "true" : "false",
             stats.analysis_seconds, stats.nblocks, stats.max_block,
             stats.structural_rank,
             stats.nnz_l, stats.nnz_u,
             stats.estimated_flops,
             stats.build_has_metis ? "true" : "false",
             stats.build_has_scotch ? "true" : "false",
             stats.build_has_spral_scaling ? "true" : "false",
             stats.tuning_profile_active ? "true" : "false",
             stats.tuning_profile_id, stats.tuning_profile_field_count, status);
    } else {
      printf("matrix: %s\n", path);
      printf("n: %" PRId64 ", nnz: %" PRId64 "\n", a.n, a.nnz);
      printf("input index bytes: %d (requested %s)\n",
             input_index.bytes, index_mode_name(input_index_mode));
      printf("requested orientation: %s\n",
             kls_orientation_name(options.orientation));
      printf("selected orientation: %s\n",
             kls_orientation_name(stats.selected_orientation));
      printf("ordering: %s\n", kls_ordering_name(stats.selected_ordering));
      printf("requested scale: %s\n", scale_name(options.scale));
      printf("row refactor control: %s\n", row_refactor_control);
      printf("requested btf: %s\n", options.use_btf ? "on" : "off");
      printf("selected btf: %s\n", stats.selected_btf ? "on" : "off");
      printf("analysis: %.6f s\n", stats.analysis_seconds);
      printf("blocks: %" PRId64 ", max block: %" PRId64 "\n",
             stats.nblocks, stats.max_block);
      printf("structural rank: %" PRId64 "\n", stats.structural_rank);
      printf("estimated nnz(L): %" PRId64 ", nnz(U): %" PRId64 "\n",
             stats.nnz_l, stats.nnz_u);
      printf("estimated flops: %.6e\n", stats.estimated_flops);
    }
    kls_destroy(solver);
    bench_index_view_free(&input_index);
    matrix_free(&a);
    return EXIT_SUCCESS;
  }

  const double *run_values = a.values;

  double *x_true = (double *)malloc((size_t)a.n * sizeof(double));
  double *b = (double *)calloc((size_t)a.n, sizeof(double));
  double *x = (double *)calloc((size_t)a.n, sizeof(double));
  if (x_true == NULL || b == NULL || x == NULL) {
    bench_index_view_free(&input_index);
    matrix_free(&a);
    free(x_true);
    free(b);
    free(x);
    return EXIT_FAILURE;
  }
  for (int64_t i = 0; i < a.n; ++i) {
    x_true[i] = 1.0 + (double)(i % 17) * 0.01;
  }
  matvec_values(&a, run_values, x_true, b);

  kls_solver *solver = NULL;
  int status = kls_create(&solver);
  if (status == KLS_OK) {
    status = kls_analyze_csc(solver, input_index.type, a.n,
                             input_index.col_ptr, input_index.row_idx, 0,
                             &options);
  }
  if (status == KLS_OK) {
    status = kls_factor(solver, run_values);
  }
  if (status != KLS_OK) {
    fprintf(stderr, "KLS setup failed: %s (%d)\n", kls_status_string(status), status);
    kls_destroy(solver);
    bench_index_view_free(&input_index);
    matrix_free(&a);
    free(x_true);
    free(b);
    free(x);
    return EXIT_FAILURE;
  }

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  kls_get_stats(solver, &stats);
  const double initial_factor_seconds = stats.factor_seconds;
  const kls_factor_path initial_factor_path = stats.last_factor_path;

  double *generated_values = NULL;
  const double *current_values = run_values;
  if (refactor_value_mode != BENCH_REFACTOR_VALUES_UNCHANGED) {
    generated_values =
      (double *)malloc((size_t)a.nnz * sizeof(*generated_values));
    if (generated_values == NULL) {
      kls_destroy(solver);
      bench_index_view_free(&input_index);
      matrix_free(&a);
      free(x_true);
      free(b);
      free(x);
      return EXIT_FAILURE;
    }
  }

  double factor_total = 0.0;
  double refactor_total = 0.0;
  double refactor_solve_total = 0.0;
  double solve_total = 0.0;
  double tsolve_total = 0.0;
  double refactor_max_relative_residual = 0.0;
  const int verify_each_refactor =
      bench_env_enabled("KLS_BENCH_VERIFY_EACH_REFACTOR");

  /* Direct lifecycle mode measures the initial system in its actual place:
     analyze + factor + solve, followed by changed-value refactor/solve pairs.
     It does not substitute a final-state solve or project sampled averages. */
  if (lifecycle_systems > 0) {
    status = kls_solve(solver, 1, b, 0, x, 0);
    if (status == KLS_OK) {
      kls_get_stats(solver, &stats);
      solve_total = stats.solve_seconds;
      if (verify_each_refactor) {
        double relative = 0.0;
        (void)residual_norm_values(&a, current_values, x, b, &relative);
        refactor_max_relative_residual = isfinite(relative) ? relative : INFINITY;
      }
    }
  }

  for (int i = 0; i < factor_repeat; ++i) {
    status = kls_factor(solver, run_values);
    if (status != KLS_OK) break;
    kls_get_stats(solver, &stats);
    factor_total += stats.factor_seconds;
  }
  double refactor_first = 0.0;
  double refactor_solve_first = 0.0;
  for (int i = 0; i < refactor_repeat && status == KLS_OK; ++i) {
    if (generated_values != NULL) {
      make_refactor_values(&a, run_values, refactor_value_mode,
                           (uint64_t)i + 1u,
                           refactor_value_amplitude, generated_values);
      current_values = generated_values;
      matvec_values(&a, current_values, x_true, b);
    }
    status = kls_refactor_solve(solver, current_values, 1, b, 0, x, 0);
    if (status != KLS_OK) break;
    kls_get_stats(solver, &stats);
    refactor_total += stats.refactor_seconds;
    if (i == 0) {
      /* the first refactorization carries the deferred engine preps;
         reporting it separately lets short suites reconstruct the
         99-iteration SPICE cycle exactly:
         cycle = shot + (first + solve) + 98*(steady + solve) */
      refactor_first = stats.refactor_seconds;
    }
    /* SPICE-shaped: every refactor is followed by a solve.  Keep these
       solve samples separate from refactor_total so H100 can charge the
       actual changed-numeric pairs, including any solve-side recovery,
       while the repeated final-state loop below supplies the initial-state
       solve sample. */
    if (status == KLS_OK) {
      kls_get_stats(solver, &stats);
      refactor_solve_total += stats.solve_seconds;
      if (i == 0) {
        refactor_solve_first = stats.solve_seconds;
      }
    }
    if (status == KLS_OK && verify_each_refactor) {
      double relative = 0.0;
      (void)residual_norm_values(&a, current_values, x, b, &relative);
      if (!isfinite(relative)) {
        refactor_max_relative_residual = INFINITY;
      } else if (relative > refactor_max_relative_residual) {
        refactor_max_relative_residual = relative;
      }
    }
  }
  const int standalone_solve_repeats = lifecycle_systems > 0 ? 0 : repeat;
  for (int i = 0; i < standalone_solve_repeats && status == KLS_OK; ++i) {
    status = kls_solve(solver, 1, b, 0, x, 0);
    if (status != KLS_OK) break;
    kls_get_stats(solver, &stats);
    solve_total += stats.solve_seconds;
  }
  if (!transpose_solve) {
    tsolve_total = -1.0 * (double)repeat;
  } else {
    for (int i = 0; i < repeat && status == KLS_OK; ++i) {
      status = kls_solve_transpose(solver, 1, b, 0, x, 0);
      if (status == KLS_ERR_UNSUPPORTED) {
        status = KLS_OK;
        tsolve_total = -1.0 * (double)repeat;
        break;
      }
      if (status != KLS_OK) break;
      kls_get_stats(solver, &stats);
      tsolve_total += stats.solve_seconds;
    }
  }
  if (status != KLS_OK) {
    fprintf(stderr, "KLS benchmark failed: %s (%d)\n", kls_status_string(status), status);
    kls_destroy(solver);
    bench_index_view_free(&input_index);
    matrix_free(&a);
    free(generated_values);
    free(x_true);
    free(b);
    free(x);
    return EXIT_FAILURE;
  }
  if (verify_each_refactor) {
    fprintf(stderr, "KLS refactor max relative residual: %.17g\n",
            refactor_max_relative_residual);
  }

  status = kls_solve(solver, 1, b, 0, x, 0);
  if (status != KLS_OK) {
    fprintf(stderr, "KLS final solve failed: %s (%d)\n", kls_status_string(status), status);
    kls_destroy(solver);
    bench_index_view_free(&input_index);
    matrix_free(&a);
    free(generated_values);
    free(x_true);
    free(b);
    free(x);
    return EXIT_FAILURE;
  }
  kls_get_stats(solver, &stats);

  double rel_residual = 0.0;
  const double residual =
      residual_norm_values(&a, current_values, x, b, &rel_residual);
  const double factor_avg = factor_repeat > 0
    ? factor_total / (double)factor_repeat : initial_factor_seconds;
  const double refactor_avg = refactor_repeat > 0 ? refactor_total / (double)refactor_repeat : 0.0;
  const double refactor_steady_avg = refactor_repeat > 1
    ? (refactor_total - refactor_first) / (double)(refactor_repeat - 1)
    : refactor_first;
  const double solve_avg = solve_total / (double)repeat;
  const double refactor_solve_avg = refactor_repeat > 0
    ? refactor_solve_total / (double)refactor_repeat : solve_avg;
  const double refactor_solve_first_effective = refactor_repeat > 0
    ? refactor_solve_first : solve_avg;
  const double refactor_solve_steady_avg = refactor_repeat > 1
    ? (refactor_solve_total - refactor_solve_first) /
        (double)(refactor_repeat - 1)
    : refactor_solve_first_effective;
  const double spice_cycle_seconds =
    stats.analysis_seconds + initial_factor_seconds + solve_avg +
    refactor_first + refactor_solve_first_effective +
    98.0 * (refactor_steady_avg + refactor_solve_steady_avg);
  const double measured_lifecycle_seconds = lifecycle_systems > 0
    ? stats.analysis_seconds + initial_factor_seconds + solve_total +
        refactor_total + refactor_solve_total
    : -1.0;
  const double tsolve_avg = tsolve_total / (double)repeat;

  if (json) {
    printf("{\"matrix\":\"%s\",\"n\":%" PRId64 ",\"nnz\":%" PRId64
           ",\"threads\":%d,\"repeat\":%d,\"factor_repeat\":%d"
           ",\"refactor_repeat\":%d"
           ",\"lifecycle_mode\":\"%s\",\"lifecycle_systems\":%d"
           ",\"measured_lifecycle_seconds\":%.9g"
           ",\"refactor_value_mode\":\"%s\""
           ",\"refactor_value_amplitude\":%.9g"
           ",\"backend\":\"%s\",\"requested_input_index\":\"%s\""
           ",\"input_index_bytes\":%d,\"internal_index_bytes\":%d"
           ",\"requested_orientation\":\"%s\",\"orientation\":\"%s\""
           ",\"ordering\":\"%s\",\"requested_scale\":\"%s\",\"scale\":%d"
           ",\"requested_btf\":%s,\"btf\":%s"
           ",\"initial_factor_path\":\"%s\",\"last_factor_path\":\"%s\""
           ",\"last_refactor_path\":\"%s\""
           ",\"analysis_seconds\":%.9g,\"initial_factor_seconds\":%.9g"
           ",\"factor_seconds_avg\":%.9g"
           ",\"refactor_first_seconds\":%.9g"
           ",\"refactor_steady_seconds_avg\":%.9g"
           ",\"refactor_seconds_avg\":%.9g"
           ",\"refactor_solve_first_seconds\":%.9g"
           ",\"refactor_solve_steady_seconds_avg\":%.9g"
           ",\"refactor_solve_seconds_avg\":%.9g"
           ",\"solve_seconds_avg\":%.9g,\"transpose_solve_seconds_avg\":%.9g"
           ",\"spice_cycle_seconds\":%.9g"
           ",\"residual_l2\":%.9g,\"relative_residual_l2\":%.9g"
           ",\"verify_each_refactor\":%s"
           ",\"refactor_max_relative_residual\":%.9g"
           ",\"nnz_l\":%" PRId64 ",\"nnz_u\":%" PRId64
           ",\"estimated_flops\":%.9g,\"factor_flops\":%.9g"
           ",\"rcond\":%.9g,\"rgrowth\":%.9g"
           ",\"memory_bytes\":%zu,\"memory_peak_bytes\":%zu"
           ",\"build_has_metis\":%s,\"build_has_scotch\":%s"
           ",\"build_has_spral_scaling\":%s"
           ",\"tuning_profile_active\":%s"
           ",\"tuning_profile_id\":\"%016" PRIx64 "\""
           ",\"tuning_profile_field_count\":%d,\"status\":%d}\n",
           path, a.n, a.nnz, options.threads, repeat, factor_repeat,
           refactor_repeat, lifecycle_systems > 0 ? "direct" : "projected",
           lifecycle_systems, measured_lifecycle_seconds,
           bench_refactor_value_mode_name(refactor_value_mode),
           refactor_value_mode != BENCH_REFACTOR_VALUES_UNCHANGED
             ? refactor_value_amplitude : 0.0,
           kls_backend_name(options.backend), index_mode_name(input_index_mode),
           input_index.bytes, stats.internal_index_bytes,
           kls_orientation_name(options.orientation),
           kls_orientation_name(stats.selected_orientation),
           kls_ordering_name(stats.selected_ordering), scale_name(options.scale),
           stats.selected_scale, options.use_btf ? "true" : "false",
           stats.selected_btf ? "true" : "false",
           kls_factor_path_name(initial_factor_path),
           kls_factor_path_name(stats.last_factor_path),
           kls_refactor_path_name(stats.last_refactor_path),
           stats.analysis_seconds, initial_factor_seconds, factor_avg,
           refactor_first, refactor_steady_avg, refactor_avg,
           refactor_solve_first_effective, refactor_solve_steady_avg,
           refactor_solve_avg, solve_avg, tsolve_avg, spice_cycle_seconds,
           residual, rel_residual,
           verify_each_refactor ? "true" : "false",
           refactor_max_relative_residual, stats.nnz_l, stats.nnz_u,
           stats.estimated_flops, stats.factor_flops, stats.rcond, stats.rgrowth,
           stats.memory_bytes, stats.memory_peak_bytes,
           stats.build_has_metis ? "true" : "false",
           stats.build_has_scotch ? "true" : "false",
           stats.build_has_spral_scaling ? "true" : "false",
           stats.tuning_profile_active ? "true" : "false",
           stats.tuning_profile_id, stats.tuning_profile_field_count, status);
  } else {
    printf("matrix: %s\n", path);
    printf("n: %" PRId64 ", nnz: %" PRId64 ", threads: %d\n",
           a.n, a.nnz, options.threads);
    printf("backend/orientation/ordering/scale: %s/%s/%s/%d\n",
           kls_backend_name(options.backend),
           kls_orientation_name(stats.selected_orientation),
           kls_ordering_name(stats.selected_ordering), stats.selected_scale);
    printf("factor paths: initial %s, final %s; refactor path: %s\n",
           kls_factor_path_name(initial_factor_path),
           kls_factor_path_name(stats.last_factor_path),
           kls_refactor_path_name(stats.last_refactor_path));
    printf("analysis %.6g, initial factor %.6g, refactor %.6g, solve %.6g\n",
           stats.analysis_seconds, initial_factor_seconds, refactor_avg,
           solve_avg);
    printf("SPICE cycle %.6g, relative residual %.6g\n",
           spice_cycle_seconds, rel_residual);
    printf("nnz(L/U): %" PRId64 "/%" PRId64
           ", factor flops %.6g, memory %zu bytes (peak %zu)\n",
           stats.nnz_l, stats.nnz_u, stats.factor_flops,
           stats.memory_bytes, stats.memory_peak_bytes);
  }

  kls_destroy(solver);
  bench_index_view_free(&input_index);
  matrix_free(&a);
  free(generated_values);
  free(x_true);
  free(b);
  free(x);
  return EXIT_SUCCESS;
}
