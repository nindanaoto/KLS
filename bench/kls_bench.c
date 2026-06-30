#define _POSIX_C_SOURCE 200809L

#include "kls/kls.h"

#include <ctype.h>
#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#if defined(__has_include)
#if __has_include(<valgrind/callgrind.h>)
#include <valgrind/callgrind.h>
#define KLS_BENCH_HAVE_CALLGRIND 1
#endif
#endif

#ifndef KLS_BENCH_HAVE_CALLGRIND
#define CALLGRIND_START_INSTRUMENTATION ((void)0)
#define CALLGRIND_STOP_INSTRUMENTATION ((void)0)
#define CALLGRIND_ZERO_STATS ((void)0)
#define CALLGRIND_DUMP_STATS ((void)0)
#endif

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
  return KLS_ORDERING_AUTO;
}

static kls_orientation parse_orientation(const char *s) {
  if (strcmp(s, "normal") == 0) return KLS_ORIENTATION_NORMAL;
  if (strcmp(s, "transpose") == 0) return KLS_ORIENTATION_TRANSPOSE;
  return KLS_ORIENTATION_AUTO;
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

static int parse_int64_arg(const char *s, int64_t *value_out) {
  char *end = NULL;
  errno = 0;
  const long long value = strtoll(s, &end, 10);
  if (errno != 0 || end == s || *end != '\0') {
    return 0;
  }
  *value_out = (int64_t)value;
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

static int valid_kls_first_factor_control(const char *s) {
  return strcmp(s, "env") == 0 ||
         strcmp(s, "off") == 0 ||
         strcmp(s, "on") == 0;
}

static int valid_row_solve_control(const char *s) {
  return strcmp(s, "env") == 0 ||
         strcmp(s, "off") == 0 ||
         strcmp(s, "on") == 0;
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

static int apply_kls_first_factor_control(const char *s) {
  if (strcmp(s, "env") == 0) {
    return 1;
  }
  if (setenv("KLS_ENABLE_KLS_FIRST_FACTOR",
             strcmp(s, "on") == 0 ? "1" : "0", 1) != 0) {
    return 0;
  }
  return 1;
}

static int apply_row_solve_control(const char *s) {
  if (strcmp(s, "env") == 0) {
    return 1;
  }
  if (setenv("KLS_ENABLE_ROW_SOLVE_FROM_NUMERIC",
             strcmp(s, "on") == 0 ? "1" : "0", 1) != 0) {
    return 0;
  }
  return 1;
}

static double *make_stressed_values(const matrix *a,
                                    double diagonal_scale,
                                    int64_t diagonal_column,
                                    int64_t *entries_out) {
  if (a == NULL || a->values == NULL || entries_out == NULL ||
      !isfinite(diagonal_scale) || diagonal_scale < 0.0 ||
      diagonal_column < -1 || diagonal_column >= a->n) {
    return NULL;
  }
  *entries_out = 0;
  double *values = (double *)malloc((size_t)a->nnz * sizeof(*values));
  if (values == NULL) {
    return NULL;
  }
  memcpy(values, a->values, (size_t)a->nnz * sizeof(*values));
  for (int64_t col = 0; col < a->n; ++col) {
    if (diagonal_column >= 0 && col != diagonal_column) {
      continue;
    }
    for (int64_t p = a->col_ptr[col]; p < a->col_ptr[col + 1]; ++p) {
      if (a->row_idx[p] == col) {
        values[p] *= diagonal_scale;
        ++(*entries_out);
      }
    }
  }
  if (*entries_out == 0) {
    free(values);
    return NULL;
  }
  return values;
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
          "Usage: %s <matrix.mtx> [--repeat N] [--refactor-repeat N] [--threads N] [--ordering auto|amd|colamd|natural|metis|scotch] [--orientation auto|normal|transpose] [--scale auto|-1|0|1|2] [--input-index auto|32|64] [--pivot-tol T] [--row-refactor env|off|refactor|checked|all] [--kls-first-factor env|off|on] [--row-solve env|off|on] [--stress-diagonal-scale S] [--stress-diagonal-column C] [--no-btf] [--no-fast-factor] [--no-static-pivoting] [--analyze-only] [--json]\n",
          argv0);
}

int main(int argc, char **argv) {
  if (argc < 2) {
    usage(argv[0]);
    return EXIT_FAILURE;
  }
  const char *path = argv[1];
  int repeat = 5;
  int refactor_repeat = 5;
  int json = 0;
  int analyze_only = 0;
  double stress_diagonal_scale = 1.0;
  int64_t stress_diagonal_column = -1;
  const char *row_refactor_control = "env";
  const char *kls_first_factor_control = "env";
  const char *row_solve_control = "env";
  bench_index_mode input_index_mode = BENCH_INDEX_AUTO;
  kls_options options;
  kls_default_options(&options);

  for (int i = 2; i < argc; ++i) {
    if (strcmp(argv[i], "--json") == 0) {
      json = 1;
    } else if (strcmp(argv[i], "--analyze-only") == 0) {
      analyze_only = 1;
    } else if (strcmp(argv[i], "--repeat") == 0 && i + 1 < argc) {
      repeat = atoi(argv[++i]);
    } else if (strcmp(argv[i], "--refactor-repeat") == 0 && i + 1 < argc) {
      refactor_repeat = atoi(argv[++i]);
    } else if (strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
      options.threads = atoi(argv[++i]);
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
    } else if (strcmp(argv[i], "--kls-first-factor") == 0 &&
               i + 1 < argc) {
      kls_first_factor_control = argv[++i];
      if (!valid_kls_first_factor_control(kls_first_factor_control)) {
        usage(argv[0]);
        return EXIT_FAILURE;
      }
    } else if (strcmp(argv[i], "--row-solve") == 0 && i + 1 < argc) {
      row_solve_control = argv[++i];
      if (!valid_row_solve_control(row_solve_control)) {
        usage(argv[0]);
        return EXIT_FAILURE;
      }
    } else if (strcmp(argv[i], "--stress-diagonal-scale") == 0 &&
               i + 1 < argc) {
      if (!parse_nonnegative_double(argv[++i], &stress_diagonal_scale)) {
        usage(argv[0]);
        return EXIT_FAILURE;
      }
    } else if (strcmp(argv[i], "--stress-diagonal-column") == 0 &&
               i + 1 < argc) {
      if (!parse_int64_arg(argv[++i], &stress_diagonal_column) ||
          stress_diagonal_column < -1) {
        usage(argv[0]);
        return EXIT_FAILURE;
      }
    } else if (strcmp(argv[i], "--no-btf") == 0) {
      options.use_btf = 0;
    } else if (strcmp(argv[i], "--no-fast-factor") == 0) {
      options.fast_factor = 0;
    } else if (strcmp(argv[i], "--no-static-pivoting") == 0) {
      options.static_pivoting = 0;
    } else {
      usage(argv[0]);
      return EXIT_FAILURE;
    }
  }
  if (repeat <= 0 || refactor_repeat < 0 || options.threads <= 0) {
    usage(argv[0]);
    return EXIT_FAILURE;
  }
  if (!apply_row_refactor_control(row_refactor_control)) {
    perror("apply row-refactor control");
    return EXIT_FAILURE;
  }
  if (!apply_kls_first_factor_control(kls_first_factor_control)) {
    perror("apply kls-first-factor control");
    return EXIT_FAILURE;
  }
  if (!apply_row_solve_control(row_solve_control)) {
    perror("apply row-solve control");
    return EXIT_FAILURE;
  }

  matrix a = {0};
  if (!read_matrix_market(path, &a)) {
    return EXIT_FAILURE;
  }
  if (stress_diagonal_column >= a.n) {
    fprintf(stderr, "stress diagonal column %" PRId64
            " is outside matrix order %" PRId64 "\n",
            stress_diagonal_column, a.n);
    matrix_free(&a);
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
             ",\"requested_input_index\":\"%s\""
             ",\"input_index_bytes\":%d"
             ",\"internal_index_bytes\":%d"
             ",\"requested_orientation\":\"%s\",\"orientation\":\"%s\""
             ",\"ordering\":\"%s\",\"requested_scale\":\"%s\""
             ",\"row_refactor_control\":\"%s\""
             ",\"kls_first_factor_control\":\"%s\""
             ",\"row_solve_control\":\"%s\""
             ",\"requested_btf\":%s,\"btf\":%s"
             ",\"analysis_seconds\":%.9g"
             ",\"nblocks\":%" PRId64 ",\"max_block\":%" PRId64
             ",\"structural_rank\":%" PRId64
             ",\"factor_etree_block_start\":%" PRId64
             ",\"factor_etree_block_size\":%" PRId64
             ",\"factor_etree_levels\":%" PRId64
             ",\"factor_etree_max_width\":%" PRId64
             ",\"factor_etree_edges\":%" PRId64
             ",\"factor_etree_root_columns\":%" PRId64
             ",\"factor_etree_leaf_columns\":%" PRId64
             ",\"factor_etree_max_fanout\":%" PRId64
             ",\"separator_analyzed_rows\":%" PRId64
             ",\"separator_global_begin\":%" PRId64
             ",\"separator_global_end\":%" PRId64
             ",\"separator_thread_count\":%" PRId64
             ",\"separator_component_count\":%" PRId64
             ",\"separator_private_components\":%" PRId64
             ",\"separator_pipeline_components\":%" PRId64
             ",\"separator_private_rows\":%" PRId64
             ",\"separator_pipeline_rows\":%" PRId64
             ",\"separator_private_max_rows\":%" PRId64
             ",\"separator_pipeline_max_rows\":%" PRId64
             ",\"nnz_l\":%" PRId64 ",\"nnz_u\":%" PRId64
             ",\"estimated_flops\":%.9g"
             ",\"parallel_model_r1\":%.9g"
             ",\"parallel_model_r2\":%.9g"
             ",\"parallel_model_recommends_parallel\":%d"
             ",\"parallel_task_flow_threads\":%" PRId64
             ",\"parallel_task_flow_dependencies\":%" PRId64
             ",\"parallel_task_flow_work\":%.9g"
             ",\"parallel_task_flow_finish_time\":%.9g"
             ",\"parallel_task_flow_speedup\":%.9g"
             ",\"parallel_task_flow_recommends_parallel\":%d"
             ",\"analyze_only\":true}\n",
             path, a.n, a.nnz, options.threads,
             index_mode_name(input_index_mode),
             input_index.bytes,
             stats.internal_index_bytes,
             kls_orientation_name(options.orientation),
             kls_orientation_name(stats.selected_orientation),
             kls_ordering_name(stats.selected_ordering),
             scale_name(options.scale),
             row_refactor_control,
             kls_first_factor_control,
             row_solve_control,
             options.use_btf ? "true" : "false",
             stats.selected_btf ? "true" : "false",
             stats.analysis_seconds, stats.nblocks, stats.max_block,
             stats.structural_rank,
             stats.factor_etree_block_start,
             stats.factor_etree_block_size,
             stats.factor_etree_levels,
             stats.factor_etree_max_width,
             stats.factor_etree_edges,
             stats.factor_etree_root_columns,
             stats.factor_etree_leaf_columns,
             stats.factor_etree_max_fanout,
             stats.separator_analyzed_rows,
             stats.separator_global_begin,
             stats.separator_global_end,
             stats.separator_thread_count,
             stats.separator_component_count,
             stats.separator_private_components,
             stats.separator_pipeline_components,
             stats.separator_private_rows,
             stats.separator_pipeline_rows,
             stats.separator_private_max_rows,
             stats.separator_pipeline_max_rows,
             stats.nnz_l, stats.nnz_u,
             stats.estimated_flops,
             stats.parallel_model_r1,
             stats.parallel_model_r2,
             stats.parallel_model_recommends_parallel,
             stats.parallel_task_flow_threads,
             stats.parallel_task_flow_dependencies,
             stats.parallel_task_flow_work,
             stats.parallel_task_flow_finish_time,
             stats.parallel_task_flow_speedup,
             stats.parallel_task_flow_recommends_parallel);
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
      printf("KLS first-factor control: %s\n", kls_first_factor_control);
      printf("row solve control: %s\n", row_solve_control);
      printf("requested btf: %s\n", options.use_btf ? "on" : "off");
      printf("selected btf: %s\n", stats.selected_btf ? "on" : "off");
      printf("analysis: %.6f s\n", stats.analysis_seconds);
      printf("blocks: %" PRId64 ", max block: %" PRId64 "\n",
             stats.nblocks, stats.max_block);
      printf("structural rank: %" PRId64 "\n", stats.structural_rank);
      printf("factor ETree largest block: start %" PRId64 ", size %" PRId64
             ", levels %" PRId64 ", max width %" PRId64 ", edges %" PRId64
             ", roots %" PRId64 ", leaves %" PRId64
             ", max fanout %" PRId64 "\n",
             stats.factor_etree_block_start,
             stats.factor_etree_block_size,
             stats.factor_etree_levels,
             stats.factor_etree_max_width,
             stats.factor_etree_edges,
             stats.factor_etree_root_columns,
             stats.factor_etree_leaf_columns,
             stats.factor_etree_max_fanout);
      printf("separator queues: rows %" PRId64 ", global [%" PRId64
             ", %" PRId64 "), threads %" PRId64
             ", components %" PRId64 " (private %" PRId64
             ", pipeline %" PRId64 "), row split %" PRId64 "/%" PRId64
             ", max private/pipeline %" PRId64 "/%" PRId64 "\n",
             stats.separator_analyzed_rows,
             stats.separator_global_begin,
             stats.separator_global_end,
             stats.separator_thread_count,
             stats.separator_component_count,
             stats.separator_private_components,
             stats.separator_pipeline_components,
             stats.separator_private_rows,
             stats.separator_pipeline_rows,
             stats.separator_private_max_rows,
             stats.separator_pipeline_max_rows);
      printf("estimated nnz(L): %" PRId64 ", nnz(U): %" PRId64 "\n",
             stats.nnz_l, stats.nnz_u);
      printf("estimated flops: %.6e\n", stats.estimated_flops);
      printf("parallel model R1/R2: %.6g / %.6g, recommends parallel: %s\n",
             stats.parallel_model_r1,
             stats.parallel_model_r2,
             stats.parallel_model_recommends_parallel ? "yes" : "no");
      printf("NICSLU task-flow model: threads %" PRId64
             ", deps %" PRId64 ", work %.6g, finish %.6g"
             ", speedup %.6g, recommends parallel: %s\n",
             stats.parallel_task_flow_threads,
             stats.parallel_task_flow_dependencies,
             stats.parallel_task_flow_work,
             stats.parallel_task_flow_finish_time,
             stats.parallel_task_flow_speedup,
             stats.parallel_task_flow_recommends_parallel ? "yes" : "no");
    }
    kls_destroy(solver);
    bench_index_view_free(&input_index);
    matrix_free(&a);
    return EXIT_SUCCESS;
  }

  const int stress_requested =
      (stress_diagonal_scale != 1.0 || stress_diagonal_column >= 0);
  int64_t stress_entries = 0;
  double *stressed_values = NULL;
  const double *run_values = a.values;
  if (stress_requested) {
    stressed_values = make_stressed_values(&a, stress_diagonal_scale,
                                           stress_diagonal_column,
                                           &stress_entries);
    if (stressed_values == NULL) {
      fprintf(stderr, "stress diagonal selection matched no diagonal entries\n");
      bench_index_view_free(&input_index);
      matrix_free(&a);
      return EXIT_FAILURE;
    }
    run_values = stressed_values;
  }

  double *x_true = (double *)malloc((size_t)a.n * sizeof(double));
  double *b = (double *)calloc((size_t)a.n, sizeof(double));
  double *x = (double *)calloc((size_t)a.n, sizeof(double));
  if (x_true == NULL || b == NULL || x == NULL) {
    bench_index_view_free(&input_index);
    matrix_free(&a);
    free(stressed_values);
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
    status = kls_factor(solver, a.values);
  }
  if (status != KLS_OK) {
    fprintf(stderr, "KLS setup failed: %s (%d)\n", kls_status_string(status), status);
    kls_destroy(solver);
    bench_index_view_free(&input_index);
    matrix_free(&a);
    free(stressed_values);
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

  double factor_total = 0.0;
  double refactor_total = 0.0;
  double solve_total = 0.0;
  double tsolve_total = 0.0;
  const int callgrind_refactor =
      bench_env_enabled("KLS_BENCH_CALLGRIND_REFACTOR");

  for (int i = 0; i < repeat; ++i) {
    status = kls_factor(solver, run_values);
    if (status != KLS_OK) break;
    kls_get_stats(solver, &stats);
    factor_total += stats.factor_seconds;
  }
  if (callgrind_refactor) {
    CALLGRIND_START_INSTRUMENTATION;
    CALLGRIND_ZERO_STATS;
  }
  for (int i = 0; i < refactor_repeat && status == KLS_OK; ++i) {
    status = kls_refactor(solver, run_values);
    if (status != KLS_OK) break;
    kls_get_stats(solver, &stats);
    refactor_total += stats.refactor_seconds;
  }
  if (callgrind_refactor) {
    CALLGRIND_DUMP_STATS;
    CALLGRIND_STOP_INSTRUMENTATION;
  }
  for (int i = 0; i < repeat && status == KLS_OK; ++i) {
    status = kls_solve(solver, 1, b, 0, x, 0);
    if (status != KLS_OK) break;
    kls_get_stats(solver, &stats);
    solve_total += stats.solve_seconds;
  }
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
  if (status != KLS_OK) {
    fprintf(stderr, "KLS benchmark failed: %s (%d)\n", kls_status_string(status), status);
    kls_destroy(solver);
    bench_index_view_free(&input_index);
    matrix_free(&a);
    free(stressed_values);
    free(x_true);
    free(b);
    free(x);
    return EXIT_FAILURE;
  }

  status = kls_solve(solver, 1, b, 0, x, 0);
  if (status != KLS_OK) {
    fprintf(stderr, "KLS final solve failed: %s (%d)\n", kls_status_string(status), status);
    kls_destroy(solver);
    bench_index_view_free(&input_index);
    matrix_free(&a);
    free(stressed_values);
    free(x_true);
    free(b);
    free(x);
    return EXIT_FAILURE;
  }
  kls_get_stats(solver, &stats);

  double rel_residual = 0.0;
  const double residual =
      residual_norm_values(&a, run_values, x, b, &rel_residual);
  const double factor_avg = factor_total / (double)repeat;
  const double refactor_avg = refactor_repeat > 0 ? refactor_total / (double)refactor_repeat : 0.0;
  const double solve_avg = solve_total / (double)repeat;
  const double tsolve_avg = tsolve_total / (double)repeat;

  if (json) {
    printf("{\"matrix\":\"%s\",\"n\":%" PRId64 ",\"nnz\":%" PRId64
           ",\"threads\":%d"
           ",\"build_has_metis\":%s"
           ",\"build_has_scotch\":%s"
           ",\"build_has_spral_scaling\":%s"
           ",\"build_has_cblas\":%s"
           ",\"requested_input_index\":\"%s\""
           ",\"input_index_bytes\":%d"
           ",\"internal_index_bytes\":%d"
           ",\"requested_orientation\":\"%s\",\"orientation\":\"%s\""
           ",\"ordering\":\"%s\",\"requested_scale\":\"%s\",\"scale\":%d"
           ",\"pivot_tolerance\":%.9g,\"selected_pivot_tolerance\":%.9g"
           ",\"row_refactor_control\":\"%s\""
           ",\"kls_first_factor_control\":\"%s\""
           ",\"row_solve_control\":\"%s\""
           ",\"stress_diagonal_scale\":%.9g"
           ",\"stress_diagonal_column\":%" PRId64
           ",\"stress_diagonal_entries\":%" PRId64
           ",\"requested_btf\":%s,\"btf\":%s,\"fast_factor\":%s"
           ",\"static_pivoting\":%s,\"selected_static_pivoting\":%s"
           ",\"selected_exact_matching\":%s"
           ",\"selected_exact_matching_scaling\":%s"
           ",\"selected_spral_matching\":%s"
           ",\"analysis_seconds\":%.9g"
           ",\"initial_factor_seconds\":%.9g"
           ",\"initial_factor_path\":\"%s\""
           ",\"last_factor_path\":\"%s\""
           ",\"last_refactor_path\":\"%s\""
           ",\"kls_tail_last_mapped_columns\":%" PRId64
           ",\"kls_tail_mapped_column_count\":%" PRId64
           ",\"kls_first_last_row_uplooking_columns\":%" PRId64
           ",\"kls_first_row_uplooking_column_count\":%" PRId64
           ",\"kls_first_last_dominant_btf_pipeline\":%d"
           ",\"kls_first_dominant_btf_pipeline_count\":%" PRId64
           ",\"kls_first_last_dominant_btf_pipeline_block\":%" PRId64
           ",\"kls_first_last_dominant_btf_pipeline_rows\":%" PRId64
           ",\"kls_first_last_dominant_btf_pipeline_has_separator\":%d"
           ",\"kls_first_dominant_btf_pipeline_without_separator_count\":%" PRId64
           ",\"kls_first_last_row_refactor_seeded_rows\":%" PRId64
           ",\"kls_first_row_refactor_seeded_row_count\":%" PRId64
           ",\"kls_first_last_dynamic_column_pivots\":%" PRId64
           ",\"kls_first_dynamic_column_pivot_count\":%" PRId64
           ",\"kls_first_last_separator_dynamic_column_pivots\":%" PRId64
           ",\"kls_first_separator_dynamic_column_pivot_count\":%" PRId64
           ",\"kls_first_last_separator_extent_dynamic_column_pivots\":%" PRId64
           ",\"kls_first_separator_extent_dynamic_column_pivot_count\":%" PRId64
           ",\"kls_first_last_separator_dynamic_column_fallbacks\":%" PRId64
           ",\"kls_first_separator_dynamic_column_fallback_count\":%" PRId64
           ",\"kls_first_last_separator_dynamic_column_rejects\":%" PRId64
           ",\"kls_first_separator_dynamic_column_reject_count\":%" PRId64
           ",\"kls_first_last_separator_queue\":%d"
           ",\"kls_first_separator_queue_run_count\":%" PRId64
           ",\"kls_first_last_separator_queue_private_components\":%" PRId64
           ",\"kls_first_last_separator_queue_pipeline_components\":%" PRId64
           ",\"kls_first_last_separator_queue_private_rows\":%" PRId64
           ",\"kls_first_last_separator_queue_pipeline_rows\":%" PRId64
           ",\"kls_first_last_separator_queue_nonempty_threads\":%" PRId64
           ",\"kls_first_last_separator_queue_max_thread_rows\":%" PRId64
           ",\"kls_first_last_separator_queue_min_thread_work\":%.9g"
           ",\"kls_first_last_separator_queue_max_thread_work\":%.9g"
           ",\"kls_first_last_separator_queue_partitioned\":%d"
           ",\"kls_first_separator_queue_partitioned_count\":%" PRId64
           ",\"kls_first_last_separator_queue_split_components\":%" PRId64,
           path, a.n, a.nnz, options.threads,
           stats.build_has_metis ? "true" : "false",
           stats.build_has_scotch ? "true" : "false",
           stats.build_has_spral_scaling ? "true" : "false",
           stats.build_has_cblas ? "true" : "false",
           index_mode_name(input_index_mode),
           input_index.bytes,
           stats.internal_index_bytes,
           kls_orientation_name(options.orientation),
           kls_orientation_name(stats.selected_orientation),
           kls_ordering_name(stats.selected_ordering),
           scale_name(options.scale), stats.selected_scale,
           options.pivot_tolerance, stats.selected_pivot_tolerance,
           row_refactor_control,
           kls_first_factor_control,
           row_solve_control,
           stress_requested ? stress_diagonal_scale : 1.0,
           stress_requested ? stress_diagonal_column : -1,
           stress_entries,
           options.use_btf ? "true" : "false",
           stats.selected_btf ? "true" : "false",
           options.fast_factor ? "true" : "false",
           options.static_pivoting ? "true" : "false",
           stats.selected_static_pivoting ? "true" : "false",
           stats.selected_exact_matching ? "true" : "false",
           stats.selected_exact_matching_scaling ? "true" : "false",
           stats.selected_spral_matching ? "true" : "false",
           stats.analysis_seconds, initial_factor_seconds,
           kls_factor_path_name(initial_factor_path),
           kls_factor_path_name(stats.last_factor_path),
           kls_refactor_path_name(stats.last_refactor_path),
           stats.kls_tail_last_mapped_columns,
           stats.kls_tail_mapped_column_count,
           stats.kls_first_last_row_uplooking_columns,
           stats.kls_first_row_uplooking_column_count,
           stats.kls_first_last_dominant_btf_pipeline,
           stats.kls_first_dominant_btf_pipeline_count,
           stats.kls_first_last_dominant_btf_pipeline_block,
           stats.kls_first_last_dominant_btf_pipeline_rows,
           stats.kls_first_last_dominant_btf_pipeline_has_separator,
           stats.kls_first_dominant_btf_pipeline_without_separator_count,
           stats.kls_first_last_row_refactor_seeded_rows,
           stats.kls_first_row_refactor_seeded_row_count,
           stats.kls_first_last_dynamic_column_pivots,
           stats.kls_first_dynamic_column_pivot_count,
           stats.kls_first_last_separator_dynamic_column_pivots,
           stats.kls_first_separator_dynamic_column_pivot_count,
           stats.kls_first_last_separator_extent_dynamic_column_pivots,
           stats.kls_first_separator_extent_dynamic_column_pivot_count,
           stats.kls_first_last_separator_dynamic_column_fallbacks,
           stats.kls_first_separator_dynamic_column_fallback_count,
           stats.kls_first_last_separator_dynamic_column_rejects,
           stats.kls_first_separator_dynamic_column_reject_count,
           stats.kls_first_last_separator_queue,
           stats.kls_first_separator_queue_run_count,
           stats.kls_first_last_separator_queue_private_components,
           stats.kls_first_last_separator_queue_pipeline_components,
           stats.kls_first_last_separator_queue_private_rows,
           stats.kls_first_last_separator_queue_pipeline_rows,
           stats.kls_first_last_separator_queue_nonempty_threads,
           stats.kls_first_last_separator_queue_max_thread_rows,
           stats.kls_first_last_separator_queue_min_thread_work,
           stats.kls_first_last_separator_queue_max_thread_work,
           stats.kls_first_last_separator_queue_partitioned,
           stats.kls_first_separator_queue_partitioned_count,
           stats.kls_first_last_separator_queue_split_components);
    printf(",\"kls_first_last_separator_queue_executed\":%d"
           ",\"kls_first_separator_queue_executed_run_count\":%" PRId64
           ",\"kls_first_last_separator_queue_executed_private_rows\":%" PRId64
           ",\"kls_first_last_separator_queue_executed_pipeline_rows\":%" PRId64
           ",\"kls_first_last_separator_queue_parallel_private\":%d"
           ",\"kls_first_separator_queue_parallel_private_run_count\":%" PRId64
           ",\"kls_first_last_separator_queue_parallel_private_rows\":%" PRId64
           ",\"kls_first_last_separator_queue_parallel_private_threads\":%" PRId64
           ",\"kls_first_last_parallel_btf_blocks\":%" PRId64
           ",\"kls_first_parallel_btf_block_count\":%" PRId64
           ",\"kls_first_auto_skipped_scaled_single_block\":%d"
           ",\"kls_first_auto_skipped_scaled_single_block_count\":%" PRId64
           ",\"factor_seconds_avg\":%.9g,\"refactor_seconds_avg\":%.9g"
           ",\"solve_seconds_avg\":%.9g,\"transpose_solve_seconds_avg\":%.9g"
           ",\"residual_l2\":%.9g,\"relative_residual_l2\":%.9g"
           ",\"nblocks\":%" PRId64 ",\"max_block\":%" PRId64
           ",\"structural_rank\":%" PRId64 ",\"numerical_rank\":%" PRId64
           ",\"offdiag_pivots\":%" PRId64 ",\"reallocations\":%" PRId64
           ",\"fast_rejected_pivot\":%" PRId64
           ",\"fast_rejected_pivot_col\":%" PRId64
           ",\"fast_rejected_row\":%" PRId64
           ",\"fast_rejected_multiplier_abs\":%.9g"
           ",\"fast_rejected_pivot_abs\":%.9g"
           ",\"fast_rejected_candidate_abs\":%.9g"
           ",\"fast_rejected_tail_candidate_row\":%" PRId64
           ",\"fast_rejected_tail_candidate_abs\":%.9g"
           ",\"fast_rejected_tail_candidate_count\":%" PRId64
           ",\"fast_rejected_tail_candidate_position\":%" PRId64
           ",\"fast_rejected_tail_repair_ready\":%d"
           ",\"fast_repaired_pivot_row\":%" PRId64
           ",\"fast_repaired_pivot_matches_tail_candidate\":%d"
           ",\"fast_repaired_first_changed_pivot\":%" PRId64
           ",\"fast_repaired_prefix_changed_pivots\":%" PRId64
           ",\"fast_repaired_suffix_changed_pivots\":%" PRId64
           ",\"fast_repaired_tail_restart_ready\":%d"
           ",\"fast_repaired_block_work\":%.9g"
           ",\"fast_repaired_tail_restart_columns\":%" PRId64
           ",\"fast_repaired_tail_restart_work\":%.9g"
           ",\"fast_repaired_tail_restart_saved_work\":%.9g"
           ",\"fast_repaired_tail_restart_overcompute_columns\":%" PRId64
           ",\"fast_repaired_tail_restart_overcompute_work\":%.9g"
           ",\"fast_repaired_tail_restart_skipped_columns\":%" PRId64
           ",\"fast_repaired_tail_restart_skipped_work\":%.9g"
           ",\"fast_repaired_tail_restart_exact_mask\":%d"
           ",\"fast_repaired_tail_restart_etree_mask\":%d"
           ",\"fast_block_restarts\":%d"
           ",\"fast_kls_block_restarts\":%d"
           ",\"fast_kls_rebuild_restarts\":%d"
           ",\"fast_tail_restarts\":%d"
           ",\"fast_repaired_last_offdiag_suffix_refresh\":%d"
           ",\"fast_repaired_offdiag_suffix_refresh_count\":%" PRId64
           ",\"fast_repaired_offdiag_full_refresh_count\":%" PRId64
           ",\"fast_repaired_parallel_tail_blocks\":%" PRId64,
           stats.kls_first_last_separator_queue_executed,
           stats.kls_first_separator_queue_executed_run_count,
           stats.kls_first_last_separator_queue_executed_private_rows,
           stats.kls_first_last_separator_queue_executed_pipeline_rows,
           stats.kls_first_last_separator_queue_parallel_private,
           stats.kls_first_separator_queue_parallel_private_run_count,
           stats.kls_first_last_separator_queue_parallel_private_rows,
           stats.kls_first_last_separator_queue_parallel_private_threads,
           stats.kls_first_last_parallel_btf_blocks,
           stats.kls_first_parallel_btf_block_count,
           stats.kls_first_auto_skipped_scaled_single_block,
           stats.kls_first_auto_skipped_scaled_single_block_count,
           factor_avg, refactor_avg, solve_avg, tsolve_avg,
           residual, rel_residual, stats.nblocks, stats.max_block,
           stats.structural_rank, stats.numerical_rank,
           stats.offdiag_pivots, stats.reallocations,
           stats.fast_rejected_pivot, stats.fast_rejected_pivot_col,
           stats.fast_rejected_row,
           stats.fast_rejected_multiplier_abs,
           stats.fast_rejected_pivot_abs,
           stats.fast_rejected_candidate_abs,
           stats.fast_rejected_tail_candidate_row,
           stats.fast_rejected_tail_candidate_abs,
           stats.fast_rejected_tail_candidate_count,
           stats.fast_rejected_tail_candidate_position,
           stats.fast_rejected_tail_repair_ready,
           stats.fast_repaired_pivot_row,
           stats.fast_repaired_pivot_matches_tail_candidate,
           stats.fast_repaired_first_changed_pivot,
           stats.fast_repaired_prefix_changed_pivots,
           stats.fast_repaired_suffix_changed_pivots,
           stats.fast_repaired_tail_restart_ready,
           stats.fast_repaired_block_work,
           stats.fast_repaired_tail_restart_columns,
           stats.fast_repaired_tail_restart_work,
           stats.fast_repaired_tail_restart_saved_work,
           stats.fast_repaired_tail_restart_overcompute_columns,
           stats.fast_repaired_tail_restart_overcompute_work,
           stats.fast_repaired_tail_restart_skipped_columns,
           stats.fast_repaired_tail_restart_skipped_work,
           stats.fast_repaired_tail_restart_exact_mask,
           stats.fast_repaired_tail_restart_etree_mask,
           stats.fast_block_restarts,
           stats.fast_kls_block_restarts,
           stats.fast_kls_rebuild_restarts,
           stats.fast_tail_restarts,
           stats.fast_repaired_last_offdiag_suffix_refresh,
           stats.fast_repaired_offdiag_suffix_refresh_count,
           stats.fast_repaired_offdiag_full_refresh_count,
           stats.fast_repaired_parallel_tail_blocks);
    printf(",\"fast_kls_block_restart_last_row_pipeline\":%d"
           ",\"fast_kls_block_restart_row_pipeline_count\":%" PRId64
           ",\"fast_kls_block_restart_last_row_pipeline_rows\":%" PRId64
           ",\"fast_kls_block_restart_last_row_pipeline_threads\":%" PRId64
           ",\"fast_kls_block_restart_last_row_pipeline_prefix_rows\":%" PRId64
           ",\"fast_kls_block_restart_last_row_pipeline_suffix_rows\":%" PRId64
           ",\"fast_kls_block_restart_last_row_pipeline_gap_rows\":%" PRId64
           ",\"fast_kls_block_restart_last_row_pipeline_etree_tail\":%d"
           ",\"fast_kls_block_restart_row_pipeline_etree_tail_count\":%" PRId64
           ",\"fast_kls_block_restart_last_row_pipeline_etree_tail_rows\":%" PRId64
           ",\"fast_kls_block_restart_last_row_pipeline_etree_tail_gap_rows\":%" PRId64
           ",\"fast_kls_block_restart_last_row_pipeline_etree_tail_exact_mask\":%d"
           ",\"fast_kls_block_restart_last_row_pipeline_etree_ready\":%d"
           ",\"fast_kls_block_restart_row_pipeline_etree_ready_count\":%" PRId64
           ",\"fast_kls_block_restart_last_row_pipeline_etree_ready_rows\":%" PRId64
           ",\"fast_kls_block_restart_last_row_pipeline_etree_ready_threads\":%" PRId64
           ",\"fast_kls_block_restart_last_row_pipeline_etree_prefactor\":%d"
           ",\"fast_kls_block_restart_row_pipeline_etree_prefactor_count\":%" PRId64
           ",\"fast_kls_block_restart_last_row_pipeline_etree_prefactor_rows\":%" PRId64
           ",\"fast_kls_block_restart_last_row_pipeline_etree_prefactor_threads\":%" PRId64
           ",\"fast_kls_block_restart_last_row_pipeline_etree_prefactor_wait_rows\":%" PRId64
           ",\"fast_kls_block_restart_last_row_pipeline_etree_prefactor_wait_deps\":%" PRId64
           ",\"fast_kls_block_restart_last_row_pipeline_separator_tail_scope\":%d"
           ",\"fast_kls_block_restart_row_pipeline_separator_tail_scope_count\":%" PRId64
           ",\"fast_kls_block_restart_last_row_pipeline_separator_tail_scope_rows\":%" PRId64
           ",\"fast_kls_block_restart_last_row_pipeline_separator_queue\":%d"
           ",\"fast_kls_block_restart_row_pipeline_separator_queue_count\":%" PRId64
           ",\"fast_kls_block_restart_last_row_pipeline_separator_private_rows\":%" PRId64
           ",\"fast_kls_block_restart_last_row_pipeline_separator_pipeline_rows\":%" PRId64
           ",\"fast_kls_block_restart_last_row_pipeline_separator_private_threads\":%" PRId64
           ",\"fast_kls_block_restart_last_row_pipeline_separator_partitioned\":%d"
           ",\"fast_kls_block_restart_last_row_pipeline_separator_split_components\":%" PRId64
           ",\"fast_kls_block_restart_last_row_pipeline_pivot_tail_rows\":%" PRId64
           ",\"fast_kls_block_restart_last_row_pipeline_pivot_restarts\":%" PRId64
           ",\"fast_kls_block_restart_last_row_pipeline_supernode_update_groups\":%" PRId64
           ",\"fast_kls_block_restart_last_row_pipeline_supernode_update_rows\":%" PRId64
           ",\"fast_kls_block_restart_last_row_pipeline_supernode_panel_update_groups\":%" PRId64
           ",\"fast_kls_block_restart_last_row_pipeline_supernode_panel_update_rows\":%" PRId64,
           stats.fast_kls_block_restart_last_row_pipeline,
           stats.fast_kls_block_restart_row_pipeline_count,
           stats.fast_kls_block_restart_last_row_pipeline_rows,
           stats.fast_kls_block_restart_last_row_pipeline_threads,
           stats.fast_kls_block_restart_last_row_pipeline_prefix_rows,
           stats.fast_kls_block_restart_last_row_pipeline_suffix_rows,
           stats.fast_kls_block_restart_last_row_pipeline_gap_rows,
           stats.fast_kls_block_restart_last_row_pipeline_etree_tail,
           stats.fast_kls_block_restart_row_pipeline_etree_tail_count,
           stats.fast_kls_block_restart_last_row_pipeline_etree_tail_rows,
           stats.fast_kls_block_restart_last_row_pipeline_etree_tail_gap_rows,
           stats.fast_kls_block_restart_last_row_pipeline_etree_tail_exact_mask,
           stats.fast_kls_block_restart_last_row_pipeline_etree_ready,
           stats.fast_kls_block_restart_row_pipeline_etree_ready_count,
           stats.fast_kls_block_restart_last_row_pipeline_etree_ready_rows,
           stats.fast_kls_block_restart_last_row_pipeline_etree_ready_threads,
           stats.fast_kls_block_restart_last_row_pipeline_etree_prefactor,
           stats.fast_kls_block_restart_row_pipeline_etree_prefactor_count,
           stats.fast_kls_block_restart_last_row_pipeline_etree_prefactor_rows,
           stats.fast_kls_block_restart_last_row_pipeline_etree_prefactor_threads,
           stats.fast_kls_block_restart_last_row_pipeline_etree_prefactor_wait_rows,
           stats.fast_kls_block_restart_last_row_pipeline_etree_prefactor_wait_deps,
           stats.fast_kls_block_restart_last_row_pipeline_separator_tail_scope,
           stats.fast_kls_block_restart_row_pipeline_separator_tail_scope_count,
           stats.fast_kls_block_restart_last_row_pipeline_separator_tail_scope_rows,
           stats.fast_kls_block_restart_last_row_pipeline_separator_queue,
           stats.fast_kls_block_restart_row_pipeline_separator_queue_count,
           stats.fast_kls_block_restart_last_row_pipeline_separator_private_rows,
           stats.fast_kls_block_restart_last_row_pipeline_separator_pipeline_rows,
           stats.fast_kls_block_restart_last_row_pipeline_separator_private_threads,
           stats.fast_kls_block_restart_last_row_pipeline_separator_partitioned,
           stats.fast_kls_block_restart_last_row_pipeline_separator_split_components,
           stats.fast_kls_block_restart_last_row_pipeline_pivot_tail_rows,
           stats.fast_kls_block_restart_last_row_pipeline_pivot_restarts,
           stats.fast_kls_block_restart_last_row_pipeline_supernode_update_groups,
           stats.fast_kls_block_restart_last_row_pipeline_supernode_update_rows,
           stats.fast_kls_block_restart_last_row_pipeline_supernode_panel_update_groups,
           stats.fast_kls_block_restart_last_row_pipeline_supernode_panel_update_rows);
    printf(",\"kls_first_last_row_pipeline\":%d"
           ",\"kls_first_row_pipeline_run_count\":%" PRId64
           ",\"kls_first_last_row_pipeline_rows\":%" PRId64
           ",\"kls_first_last_row_pipeline_threads\":%" PRId64
           ",\"kls_first_last_row_pipeline_partial\":%d"
           ",\"kls_first_row_pipeline_partial_run_count\":%" PRId64
           ",\"kls_first_last_row_pipeline_partial_rows\":%" PRId64
           ",\"kls_first_last_row_pipeline_partial_threads\":%" PRId64
           ",\"kls_first_last_row_pipeline_pivot_tail\":%d"
           ",\"kls_first_row_pipeline_pivot_tail_run_count\":%" PRId64
           ",\"kls_first_last_row_pipeline_pivot_tail_rows\":%" PRId64
           ",\"kls_first_last_row_pipeline_pivot_restarts\":%" PRId64
           ",\"kls_first_row_pipeline_pivot_restart_count\":%" PRId64
           ",\"kls_first_last_row_pipeline_pivot_serial_rows\":%" PRId64
           ",\"kls_first_last_row_pipeline_prefix_panel_rebuild\":%d"
           ",\"kls_first_row_pipeline_prefix_panel_rebuild_count\":%" PRId64
           ",\"kls_first_last_row_pipeline_prefix_panel_rebuild_rows\":%" PRId64
           ",\"kls_first_active_rank_pivot_reset_count\":%" PRId64
           ",\"kls_first_active_rank_pivot_reset_rows\":%" PRId64
           ",\"kls_first_active_rank_pivot_panel_rebuild_count\":%" PRId64
           ",\"kls_first_active_rank_pivot_panel_rebuild_rows\":%" PRId64
           ",\"kls_first_row_panel_cache_build_count\":%" PRId64
           ",\"kls_first_row_panel_cache_build_panels\":%" PRId64
           ",\"kls_first_row_panel_cache_build_entries\":%" PRId64
           ",\"kls_first_row_panel_cache_append_count\":%" PRId64
           ",\"kls_first_row_panel_cache_append_panels\":%" PRId64
           ",\"kls_first_row_panel_cache_append_entries\":%" PRId64,
           stats.kls_first_last_row_pipeline,
           stats.kls_first_row_pipeline_run_count,
           stats.kls_first_last_row_pipeline_rows,
           stats.kls_first_last_row_pipeline_threads,
           stats.kls_first_last_row_pipeline_partial,
           stats.kls_first_row_pipeline_partial_run_count,
           stats.kls_first_last_row_pipeline_partial_rows,
           stats.kls_first_last_row_pipeline_partial_threads,
           stats.kls_first_last_row_pipeline_pivot_tail,
           stats.kls_first_row_pipeline_pivot_tail_run_count,
           stats.kls_first_last_row_pipeline_pivot_tail_rows,
           stats.kls_first_last_row_pipeline_pivot_restarts,
           stats.kls_first_row_pipeline_pivot_restart_count,
           stats.kls_first_last_row_pipeline_pivot_serial_rows,
           stats.kls_first_last_row_pipeline_prefix_panel_rebuild,
           stats.kls_first_row_pipeline_prefix_panel_rebuild_count,
           stats.kls_first_last_row_pipeline_prefix_panel_rebuild_rows,
           stats.kls_first_active_rank_pivot_reset_count,
           stats.kls_first_active_rank_pivot_reset_rows,
           stats.kls_first_active_rank_pivot_panel_rebuild_count,
           stats.kls_first_active_rank_pivot_panel_rebuild_rows,
           stats.kls_first_row_panel_cache_build_count,
           stats.kls_first_row_panel_cache_build_panels,
           stats.kls_first_row_panel_cache_build_entries,
           stats.kls_first_row_panel_cache_append_count,
           stats.kls_first_row_panel_cache_append_panels,
           stats.kls_first_row_panel_cache_append_entries);
    printf(",\"kls_first_last_row_supernode_update\":%d"
           ",\"kls_first_row_supernode_update_run_count\":%" PRId64
           ",\"kls_first_row_supernode_update_groups\":%" PRId64
           ",\"kls_first_row_supernode_update_rows\":%" PRId64
           ",\"kls_first_last_row_supernode_update_groups\":%" PRId64
           ",\"kls_first_last_row_supernode_update_rows\":%" PRId64
           ",\"kls_first_last_row_supernode_panel_update\":%d"
           ",\"kls_first_row_supernode_panel_update_run_count\":%" PRId64
           ",\"kls_first_row_supernode_panel_update_groups\":%" PRId64
           ",\"kls_first_row_supernode_panel_update_rows\":%" PRId64
           ",\"kls_first_last_row_supernode_panel_update_groups\":%" PRId64
           ",\"kls_first_last_row_supernode_panel_update_rows\":%" PRId64,
           stats.kls_first_last_row_supernode_update,
           stats.kls_first_row_supernode_update_run_count,
           stats.kls_first_row_supernode_update_groups,
           stats.kls_first_row_supernode_update_rows,
           stats.kls_first_last_row_supernode_update_groups,
           stats.kls_first_last_row_supernode_update_rows,
           stats.kls_first_last_row_supernode_panel_update,
           stats.kls_first_row_supernode_panel_update_run_count,
           stats.kls_first_row_supernode_panel_update_groups,
           stats.kls_first_row_supernode_panel_update_rows,
           stats.kls_first_last_row_supernode_panel_update_groups,
           stats.kls_first_last_row_supernode_panel_update_rows);
    printf(",\"kls_first_last_separator_queue_parallel_pipeline\":%d"
           ",\"kls_first_separator_queue_parallel_pipeline_run_count\":%" PRId64
           ",\"kls_first_last_separator_queue_parallel_pipeline_rows\":%" PRId64
           ",\"kls_first_last_separator_queue_parallel_pipeline_threads\":%" PRId64
           ",\"kls_first_last_separator_queue_pipeline_partial\":%d"
           ",\"kls_first_separator_queue_pipeline_partial_run_count\":%" PRId64
           ",\"kls_first_last_separator_queue_pipeline_partial_rows\":%" PRId64
           ",\"kls_first_last_separator_queue_pipeline_partial_threads\":%" PRId64
           ",\"kls_first_last_separator_queue_pipeline_wait_partial\":%d"
           ",\"kls_first_separator_queue_pipeline_wait_partial_run_count\":%" PRId64
           ",\"kls_first_last_separator_queue_pipeline_wait_partial_rows\":%" PRId64
           ",\"kls_first_last_separator_queue_pipeline_wait_partial_deps\":%" PRId64
           ",\"kls_first_last_separator_queue_pipeline_supernode_update\":%d"
           ",\"kls_first_separator_queue_pipeline_supernode_update_run_count\":%" PRId64
           ",\"kls_first_last_separator_queue_pipeline_supernode_update_groups\":%" PRId64
           ",\"kls_first_last_separator_queue_pipeline_supernode_update_rows\":%" PRId64
           ",\"kls_first_last_separator_queue_pipeline_supernode_panel_update\":%d"
           ",\"kls_first_separator_queue_pipeline_supernode_panel_update_run_count\":%" PRId64
           ",\"kls_first_last_separator_queue_pipeline_supernode_panel_update_groups\":%" PRId64
           ",\"kls_first_last_separator_queue_pipeline_supernode_panel_update_rows\":%" PRId64
           ",\"kls_first_last_separator_queue_pipeline_pivot_tail\":%d"
           ",\"kls_first_separator_queue_pipeline_pivot_tail_run_count\":%" PRId64
           ",\"kls_first_last_separator_queue_pipeline_pivot_tail_rows\":%" PRId64
           ",\"kls_first_last_separator_queue_pipeline_pivot_restarts\":%" PRId64
           ",\"kls_first_separator_queue_pipeline_pivot_restart_count\":%" PRId64
           ",\"kls_first_last_separator_queue_pipeline_pivot_serial_rows\":%" PRId64
           ",\"kls_first_last_separator_queue_pipeline_prefix_panel_rebuild\":%d"
           ",\"kls_first_separator_queue_pipeline_prefix_panel_rebuild_count\":%" PRId64
           ",\"kls_first_last_separator_queue_pipeline_prefix_panel_rebuild_rows\":%" PRId64,
           stats.kls_first_last_separator_queue_parallel_pipeline,
           stats.kls_first_separator_queue_parallel_pipeline_run_count,
           stats.kls_first_last_separator_queue_parallel_pipeline_rows,
           stats.kls_first_last_separator_queue_parallel_pipeline_threads,
           stats.kls_first_last_separator_queue_pipeline_partial,
           stats.kls_first_separator_queue_pipeline_partial_run_count,
           stats.kls_first_last_separator_queue_pipeline_partial_rows,
           stats.kls_first_last_separator_queue_pipeline_partial_threads,
           stats.kls_first_last_separator_queue_pipeline_wait_partial,
           stats.kls_first_separator_queue_pipeline_wait_partial_run_count,
           stats.kls_first_last_separator_queue_pipeline_wait_partial_rows,
           stats.kls_first_last_separator_queue_pipeline_wait_partial_deps,
           stats.kls_first_last_separator_queue_pipeline_supernode_update,
           stats.kls_first_separator_queue_pipeline_supernode_update_run_count,
           stats.kls_first_last_separator_queue_pipeline_supernode_update_groups,
           stats.kls_first_last_separator_queue_pipeline_supernode_update_rows,
           stats.kls_first_last_separator_queue_pipeline_supernode_panel_update,
           stats.kls_first_separator_queue_pipeline_supernode_panel_update_run_count,
           stats.kls_first_last_separator_queue_pipeline_supernode_panel_update_groups,
           stats.kls_first_last_separator_queue_pipeline_supernode_panel_update_rows,
           stats.kls_first_last_separator_queue_pipeline_pivot_tail,
           stats.kls_first_separator_queue_pipeline_pivot_tail_run_count,
           stats.kls_first_last_separator_queue_pipeline_pivot_tail_rows,
           stats.kls_first_last_separator_queue_pipeline_pivot_restarts,
           stats.kls_first_separator_queue_pipeline_pivot_restart_count,
           stats.kls_first_last_separator_queue_pipeline_pivot_serial_rows,
           stats.kls_first_last_separator_queue_pipeline_prefix_panel_rebuild,
           stats.kls_first_separator_queue_pipeline_prefix_panel_rebuild_count,
           stats.kls_first_last_separator_queue_pipeline_prefix_panel_rebuild_rows);
    printf(",\"fast_rejected_block_start\":%" PRId64
           ",\"fast_rejected_block_size\":%" PRId64
           ",\"fast_rejected_suffix_columns\":%" PRId64
           ",\"fast_rejected_descendant_columns\":%" PRId64
           ",\"fast_rejected_descendant_work\":%.9g"
           ",\"fast_rejected_row_tail_columns\":%" PRId64
           ",\"fast_rejected_row_tail_work\":%.9g"
           ",\"fast_rejected_group_tail_groups\":%" PRId64
           ",\"fast_rejected_group_tail_rows\":%" PRId64
           ",\"fast_rejected_group_tail_work\":%.9g"
           ",\"fast_rejected_etree_columns\":%" PRId64
           ",\"fast_rejected_etree_work\":%.9g"
           ",\"fast_rejected_pivoting_tail_columns\":%" PRId64
           ",\"fast_rejected_pivoting_tail_work\":%.9g"
           ",\"fast_rejected_pivoting_tail_first\":%" PRId64
           ",\"fast_rejected_pivoting_tail_last\":%" PRId64
           ",\"fast_rejected_pivoting_tail_contains_reject\":%d"
           ",\"fast_rejected_pivoting_tail_topological\":%d"
           ",\"fast_rejected_pivoting_tail_seed_columns\":%" PRId64
           ",\"fast_rejected_pivoting_tail_row_seed_columns\":%" PRId64
           ",\"fast_rejected_pivoting_tail_block_seed_columns\":%" PRId64
           ",\"fast_rejected_pivoting_tail_contiguous\":%d"
           ",\"fast_rejected_pivoting_tail_suffix_exact\":%d"
           ",\"fast_rejected_pivoting_tail_gap_columns\":%" PRId64
           ",\"fast_rejected_pivoting_tail_suffix_overcompute_columns\":%" PRId64
           ",\"fast_rejected_pivoting_tail_suffix_overcompute_work\":%.9g"
           ",\"fast_rejected_pivoting_tail_etree_edges\":%" PRId64
           ",\"fast_rejected_pivoting_tail_etree_roots\":%" PRId64
           ",\"fast_rejected_pivoting_tail_etree_leaves\":%" PRId64
           ",\"fast_rejected_pivoting_tail_etree_max_fanout\":%" PRId64
           ",\"fast_rejected_pivoting_tail_etree_levels\":%" PRId64
           ",\"fast_rejected_pivoting_tail_etree_max_width\":%" PRId64
           ",\"fast_rejected_refresh_state\":%d"
           ",\"fast_factor_fail_reason\":%d"
           ",\"fast_factor_fail_status\":%d"
           ",\"fast_rejected_prefix_refresh_columns\":%" PRId64
           ",\"fast_rejected_prefix_refresh_count\":%" PRId64,
           stats.fast_rejected_block_start,
           stats.fast_rejected_block_size,
           stats.fast_rejected_suffix_columns,
           stats.fast_rejected_descendant_columns,
           stats.fast_rejected_descendant_work,
           stats.fast_rejected_row_tail_columns,
           stats.fast_rejected_row_tail_work,
           stats.fast_rejected_group_tail_groups,
           stats.fast_rejected_group_tail_rows,
           stats.fast_rejected_group_tail_work,
           stats.fast_rejected_etree_columns,
           stats.fast_rejected_etree_work,
           stats.fast_rejected_pivoting_tail_columns,
           stats.fast_rejected_pivoting_tail_work,
           stats.fast_rejected_pivoting_tail_first,
           stats.fast_rejected_pivoting_tail_last,
           stats.fast_rejected_pivoting_tail_contains_reject,
           stats.fast_rejected_pivoting_tail_topological,
           stats.fast_rejected_pivoting_tail_seed_columns,
           stats.fast_rejected_pivoting_tail_row_seed_columns,
           stats.fast_rejected_pivoting_tail_block_seed_columns,
           stats.fast_rejected_pivoting_tail_contiguous,
           stats.fast_rejected_pivoting_tail_suffix_exact,
           stats.fast_rejected_pivoting_tail_gap_columns,
           stats.fast_rejected_pivoting_tail_suffix_overcompute_columns,
           stats.fast_rejected_pivoting_tail_suffix_overcompute_work,
           stats.fast_rejected_pivoting_tail_etree_edges,
           stats.fast_rejected_pivoting_tail_etree_roots,
           stats.fast_rejected_pivoting_tail_etree_leaves,
           stats.fast_rejected_pivoting_tail_etree_max_fanout,
           stats.fast_rejected_pivoting_tail_etree_levels,
           stats.fast_rejected_pivoting_tail_etree_max_width,
           stats.fast_rejected_refresh_state,
           stats.fast_factor_fail_reason,
           stats.fast_factor_fail_status,
           stats.fast_rejected_prefix_refresh_columns,
           stats.fast_rejected_prefix_refresh_count);
    printf(",\"factor_etree_block_start\":%" PRId64
           ",\"factor_etree_block_size\":%" PRId64
           ",\"factor_etree_levels\":%" PRId64
           ",\"factor_etree_max_width\":%" PRId64
           ",\"factor_etree_edges\":%" PRId64
           ",\"factor_etree_root_columns\":%" PRId64
           ",\"factor_etree_leaf_columns\":%" PRId64
           ",\"factor_etree_max_fanout\":%" PRId64
           ",\"separator_analyzed_rows\":%" PRId64
           ",\"separator_global_begin\":%" PRId64
           ",\"separator_global_end\":%" PRId64
           ",\"separator_thread_count\":%" PRId64
           ",\"separator_component_count\":%" PRId64
           ",\"separator_private_components\":%" PRId64
           ",\"separator_pipeline_components\":%" PRId64
           ",\"separator_private_rows\":%" PRId64
           ",\"separator_pipeline_rows\":%" PRId64
           ",\"separator_private_max_rows\":%" PRId64
           ",\"separator_pipeline_max_rows\":%" PRId64
           ",\"refactor_dependency_levels\":%" PRId64
           ",\"refactor_dependency_max_width\":%" PRId64
           ",\"refactor_dependency_edges\":%" PRId64
           ",\"refactor_dependency_root_columns\":%" PRId64
           ",\"refactor_dependency_leaf_columns\":%" PRId64
           ",\"refactor_dependency_max_fanout\":%" PRId64
           ",\"refactor_dependency_max_column_work\":%.9g"
           ",\"refactor_dependency_pipeline_max_column_work\":%.9g"
           ",\"refactor_supernode_candidate_count\":%" PRId64
           ",\"refactor_supernode_candidate_rows\":%" PRId64
           ",\"refactor_supernode_candidate_max_width\":%" PRId64
           ",\"refactor_supernode_candidate_dense_entries\":%.9g"
           ",\"refactor_supernode_candidate_trailing_entries\":%.9g"
           ",\"refactor_supernode_consumer_run_count\":%" PRId64
           ",\"refactor_supernode_consumer_run_rows\":%" PRId64
           ",\"refactor_supernode_consumer_run_max_width\":%" PRId64
           ",\"refactor_supernode_consumer_suffix_count\":%" PRId64
           ",\"refactor_supernode_consumer_l_entries\":%.9g"
           ",\"refactor_supernode_consumer_internal_entries\":%.9g"
           ",\"refactor_u_supernode_pattern_count\":%" PRId64
           ",\"refactor_u_supernode_pattern_rows\":%" PRId64
           ",\"refactor_u_supernode_pattern_max_width\":%" PRId64
           ",\"refactor_u_supernode_pattern_right_entries\":%" PRId64
           ",\"refactor_u_supernode_pattern_internal_entries\":%.9g"
           ",\"refactor_u_supernode_value_dense_entries\":%" PRId64
           ",\"refactor_u_supernode_value_right_entries\":%" PRId64
           ",\"refactor_last_u_supernode_value_dense_writes\":%" PRId64
           ",\"refactor_last_u_supernode_value_right_writes\":%" PRId64
           ",\"refactor_u_supernode_value_dense_write_count\":%" PRId64
           ",\"refactor_u_supernode_value_right_write_count\":%" PRId64
           ",\"refactor_u_supernode_l_panel_count\":%" PRId64
           ",\"refactor_u_supernode_l_dense_entries\":%" PRId64
           ",\"refactor_u_supernode_l_trailing_entries\":%" PRId64
           ",\"refactor_u_supernode_l_prune_count\":%" PRId64
           ",\"refactor_u_supernode_l_pruned_panels\":%" PRId64
           ",\"refactor_u_supernode_l_pruned_dense_entries\":%" PRId64
           ",\"refactor_u_supernode_l_pruned_trailing_entries\":%" PRId64
           ",\"refactor_last_u_supernode_l_update_runs\":%" PRId64
           ",\"refactor_last_u_supernode_l_update_rows\":%" PRId64
           ",\"refactor_last_u_supernode_l_update_entries\":%" PRId64
           ",\"refactor_last_u_supernode_l_probe_attempts\":%" PRId64
           ",\"refactor_last_u_supernode_l_panel_misses\":%" PRId64
           ",\"refactor_last_u_supernode_l_short_rejects\":%" PRId64
           ",\"refactor_last_u_supernode_l_stream_rejects\":%" PRId64
           ",\"refactor_last_u_supernode_l_work_rejects\":%" PRId64
           ",\"refactor_u_supernode_l_exec_disabled\":%d"
           ",\"refactor_u_supernode_l_exec_disable_count\":%" PRId64
           ",\"refactor_u_supernode_l_update_run_count\":%" PRId64
           ",\"refactor_u_supernode_l_update_rows\":%" PRId64
           ",\"refactor_u_supernode_l_update_entries\":%" PRId64
           ",\"refactor_supernode_panel_count\":%" PRId64
           ",\"refactor_supernode_panel_used_count\":%" PRId64
           ",\"refactor_supernode_consumer_plan_cached_panel_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_cached_panel_rows\":%" PRId64
           ",\"refactor_supernode_consumer_plan_strict_cached_panel_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_strict_cached_panel_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_strict_cached_panel_trailing_entries\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_built\":%d"
           ",\"refactor_supernode_consumer_plan_group_l_storage_limited\":%d"
           ",\"refactor_supernode_consumer_plan_group_l_panel_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_rows\":%" PRId64
           ",\"refactor_supernode_consumer_plan_group_l_run_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_exec_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_exec_run_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_batch_candidate_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_batch_candidate_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_batch_candidate_run_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_batch_candidate_entries\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_batch_candidate_max_runs\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_batch_candidate_advance_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_batch_candidate_advance_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_batch_candidate_advance_run_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_batch_candidate_advance_work\":%.9g"
           ",\"refactor_supernode_consumer_plan_group_l_batch_candidate_advance_max_work\":%.9g"
           ",\"refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_run_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_entries\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_advance_work\":%.9g"
           ",\"refactor_supernode_consumer_plan_group_l_batch_candidate_max_payoff_ratio\":%.9g"
           ",\"refactor_supernode_consumer_plan_group_l_dense_entries\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_trailing_entries\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_bytes\":%" PRId64
           ",\"refactor_supernode_consumer_plan_group_l_state_built\":%d"
           ",\"refactor_supernode_consumer_plan_group_l_state_storage_limited\":%d"
           ",\"refactor_supernode_consumer_plan_group_l_state_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_state_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_state_max_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_state_bytes\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_state_focus_enabled\":%d"
           ",\"refactor_supernode_consumer_plan_group_l_state_candidate_group_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_state_candidate_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_state_candidate_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_state_selected_group_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_state_selected_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_state_selected_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_state_selected_bytes\":%"
           PRId64
           ",\"refactor_last_supernode_consumer_plan_group_l_dense_writes\":%"
           PRId64
           ",\"refactor_last_supernode_consumer_plan_group_l_trailing_writes\":%"
           PRId64
           ",\"refactor_last_supernode_consumer_plan_group_l_invalidations\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_dense_write_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_trailing_write_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_invalidation_count\":%"
           PRId64
           ",\"refactor_last_supernode_consumer_plan_group_l_update_runs\":%"
           PRId64
           ",\"refactor_last_supernode_consumer_plan_group_l_update_rows\":%"
           PRId64
           ",\"refactor_last_supernode_consumer_plan_group_l_update_entries\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_update_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_update_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_update_entries\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_deferred_columns\":%" PRId64
           ",\"refactor_supernode_consumer_plan_deferred_entries\":%" PRId64
           ",\"refactor_supernode_consumer_plan_deferred_unique_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_batch_deferred_columns\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_batch_deferred_entries\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_batch_deferred_unique_rows\":%"
           PRId64
           ",\"row_refactor_group_count\":%" PRId64
           ",\"row_refactor_group_single_count\":%" PRId64
           ",\"row_refactor_group_batch_count\":%" PRId64
           ",\"row_refactor_group_batch_rows\":%" PRId64
           ",\"row_refactor_group_batch_max_width\":%" PRId64
           ",\"row_refactor_group_batch_width_le_4_count\":%" PRId64
           ",\"row_refactor_group_batch_width_le_8_count\":%" PRId64
           ",\"row_refactor_group_scalar_candidate_count\":%" PRId64
           ",\"row_refactor_group_scalar_candidate_rows\":%" PRId64
           ",\"row_refactor_group_scalar_short_count\":%" PRId64
           ",\"row_refactor_group_scalar_short_rows\":%" PRId64
           ",\"row_refactor_group_scalar_stop_level_mismatch_count\":%" PRId64
           ",\"row_refactor_group_scalar_stop_internal_dep_count\":%" PRId64
           ",\"row_refactor_group_scalar_stop_next_segment_count\":%" PRId64
           ",\"row_refactor_group_scalar_stop_max_width_count\":%" PRId64
           ",\"row_refactor_group_scalar_stop_matrix_end_count\":%" PRId64
           ",\"row_refactor_group_generic_count\":%" PRId64
           ",\"row_refactor_group_generic_rows\":%" PRId64
           ",\"row_refactor_group_generic_max_width\":%" PRId64
           ",\"row_refactor_group_dense_count\":%" PRId64
           ",\"row_refactor_group_dense_rows\":%" PRId64
           ",\"row_refactor_group_dense_max_width\":%" PRId64
           ",\"row_refactor_group_single_work\":%.9g"
           ",\"row_refactor_group_batch_work\":%.9g"
           ",\"row_refactor_group_generic_work\":%.9g"
           ",\"row_refactor_group_dense_work\":%.9g"
           ",\"row_refactor_group_level_count\":%" PRId64
           ",\"row_refactor_group_level_max_width\":%" PRId64
           ",\"row_refactor_group_cluster_levels\":%" PRId64
           ",\"row_refactor_group_pipeline_groups\":%" PRId64
           ",\"row_refactor_group_pipeline_rows\":%" PRId64
           ",\"row_refactor_group_pipeline_work\":%.9g"
           ",\"row_refactor_total_group_work\":%.9g"
           ",\"row_refactor_auto_enabled\":%d"
           ",\"row_refactor_auto_values_ready\":%d"
           ",\"row_refactor_auto_work_allowed\":%d"
           ",\"row_refactor_auto_should_run\":%d"
           ",\"row_refactor_auto_lower_bound_work\":%.9g"
           ",\"row_refactor_auto_lower_bound_rejected\":%d"
           ",\"row_refactor_auto_pattern_build_failed\":%d"
           ",\"row_refactor_auto_value_copy_failed\":%d"
           ",\"row_refactor_auto_model_recommended\":%d"
           ",\"row_refactor_auto_model_attempted\":%d"
           ",\"row_refactor_auto_model_accepted\":%d"
           ",\"row_refactor_group_dependency_edges\":%" PRId64
           ",\"row_refactor_group_root_count\":%" PRId64
           ",\"row_refactor_group_leaf_count\":%" PRId64
           ",\"row_refactor_group_max_fanout\":%" PRId64
           ",\"row_refactor_last_run\":%d"
           ",\"row_refactor_last_checked\":%d"
           ",\"row_refactor_last_parallel\":%d"
           ",\"row_refactor_last_ready_queue\":%d"
           ",\"row_refactor_run_count\":%" PRId64
           ",\"row_refactor_checked_run_count\":%" PRId64
           ",\"row_refactor_parallel_run_count\":%" PRId64
           ",\"row_refactor_ready_queue_run_count\":%" PRId64
           ",\"row_refactor_ready_queue_group_count\":%" PRId64
           ",\"row_refactor_last_done_bitmap\":%d"
           ",\"row_refactor_done_bitmap_run_count\":%" PRId64
           ",\"row_refactor_last_prefactor\":%d"
           ",\"row_refactor_last_prefactor_rows\":%" PRId64
           ",\"row_refactor_last_prefactor_deps\":%" PRId64
           ",\"row_refactor_prefactor_run_count\":%" PRId64
           ",\"row_refactor_prefactor_rows\":%" PRId64
           ",\"row_refactor_prefactor_deps\":%" PRId64
           ",\"row_refactor_last_prefactor_supernode\":%d"
           ",\"row_refactor_last_prefactor_supernode_rows\":%" PRId64
           ",\"row_refactor_last_prefactor_supernode_deps\":%" PRId64
           ",\"row_refactor_prefactor_supernode_run_count\":%" PRId64
           ",\"row_refactor_prefactor_supernode_rows\":%" PRId64
           ",\"row_refactor_prefactor_supernode_deps\":%" PRId64
           ",\"row_refactor_input_cleanup_rows\":%" PRId64
           ",\"row_refactor_input_cleanup_entries\":%" PRId64
           ",\"row_refactor_last_defer_value_scatter\":%d"
           ",\"row_refactor_defer_value_scatter_run_count\":%" PRId64
           ",\"row_refactor_values_dirty\":%d"
           ",\"row_refactor_last_lazy_value_scatter\":%d"
           ",\"row_refactor_lazy_value_scatter_run_count\":%" PRId64
           ",\"row_refactor_last_row_solve\":%d"
           ",\"row_refactor_row_solve_run_count\":%" PRId64,
           stats.factor_etree_block_start,
           stats.factor_etree_block_size,
           stats.factor_etree_levels,
           stats.factor_etree_max_width,
           stats.factor_etree_edges,
           stats.factor_etree_root_columns,
           stats.factor_etree_leaf_columns,
           stats.factor_etree_max_fanout,
           stats.separator_analyzed_rows,
           stats.separator_global_begin,
           stats.separator_global_end,
           stats.separator_thread_count,
           stats.separator_component_count,
           stats.separator_private_components,
           stats.separator_pipeline_components,
           stats.separator_private_rows,
           stats.separator_pipeline_rows,
           stats.separator_private_max_rows,
           stats.separator_pipeline_max_rows,
           stats.refactor_dependency_levels,
           stats.refactor_dependency_max_width,
           stats.refactor_dependency_edges,
           stats.refactor_dependency_root_columns,
           stats.refactor_dependency_leaf_columns,
           stats.refactor_dependency_max_fanout,
           stats.refactor_dependency_max_column_work,
           stats.refactor_dependency_pipeline_max_column_work,
           stats.refactor_supernode_candidate_count,
           stats.refactor_supernode_candidate_rows,
           stats.refactor_supernode_candidate_max_width,
           stats.refactor_supernode_candidate_dense_entries,
           stats.refactor_supernode_candidate_trailing_entries,
           stats.refactor_supernode_consumer_run_count,
           stats.refactor_supernode_consumer_run_rows,
           stats.refactor_supernode_consumer_run_max_width,
           stats.refactor_supernode_consumer_suffix_count,
           stats.refactor_supernode_consumer_l_entries,
           stats.refactor_supernode_consumer_internal_entries,
           stats.refactor_u_supernode_pattern_count,
           stats.refactor_u_supernode_pattern_rows,
           stats.refactor_u_supernode_pattern_max_width,
           stats.refactor_u_supernode_pattern_right_entries,
           stats.refactor_u_supernode_pattern_internal_entries,
           stats.refactor_u_supernode_value_dense_entries,
           stats.refactor_u_supernode_value_right_entries,
           stats.refactor_last_u_supernode_value_dense_writes,
           stats.refactor_last_u_supernode_value_right_writes,
           stats.refactor_u_supernode_value_dense_write_count,
           stats.refactor_u_supernode_value_right_write_count,
           stats.refactor_u_supernode_l_panel_count,
           stats.refactor_u_supernode_l_dense_entries,
           stats.refactor_u_supernode_l_trailing_entries,
           stats.refactor_u_supernode_l_prune_count,
           stats.refactor_u_supernode_l_pruned_panels,
           stats.refactor_u_supernode_l_pruned_dense_entries,
           stats.refactor_u_supernode_l_pruned_trailing_entries,
           stats.refactor_last_u_supernode_l_update_runs,
           stats.refactor_last_u_supernode_l_update_rows,
           stats.refactor_last_u_supernode_l_update_entries,
           stats.refactor_last_u_supernode_l_probe_attempts,
           stats.refactor_last_u_supernode_l_panel_misses,
           stats.refactor_last_u_supernode_l_short_rejects,
           stats.refactor_last_u_supernode_l_stream_rejects,
           stats.refactor_last_u_supernode_l_work_rejects,
           stats.refactor_u_supernode_l_exec_disabled,
           stats.refactor_u_supernode_l_exec_disable_count,
           stats.refactor_u_supernode_l_update_run_count,
           stats.refactor_u_supernode_l_update_rows,
           stats.refactor_u_supernode_l_update_entries,
           stats.refactor_supernode_panel_count,
           stats.refactor_supernode_panel_used_count,
           stats.refactor_supernode_consumer_plan_cached_panel_count,
           stats.refactor_supernode_consumer_plan_cached_panel_rows,
           stats.refactor_supernode_consumer_plan_strict_cached_panel_count,
           stats.refactor_supernode_consumer_plan_strict_cached_panel_rows,
           stats.refactor_supernode_consumer_plan_strict_cached_panel_trailing_entries,
           stats.refactor_supernode_consumer_plan_group_l_built,
           stats.refactor_supernode_consumer_plan_group_l_storage_limited,
           stats.refactor_supernode_consumer_plan_group_l_panel_count,
           stats.refactor_supernode_consumer_plan_group_l_run_count,
           stats.refactor_supernode_consumer_plan_group_l_rows,
           stats.refactor_supernode_consumer_plan_group_l_run_rows,
           stats.refactor_supernode_consumer_plan_group_l_exec_run_count,
           stats.refactor_supernode_consumer_plan_group_l_exec_run_rows,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_count,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_run_count,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_run_rows,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_entries,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_max_runs,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_advance_count,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_advance_run_count,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_advance_run_rows,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_advance_work,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_advance_max_work,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_count,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_run_count,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_run_rows,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_entries,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_advance_work,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_max_payoff_ratio,
           stats.refactor_supernode_consumer_plan_group_l_dense_entries,
           stats.refactor_supernode_consumer_plan_group_l_trailing_entries,
           stats.refactor_supernode_consumer_plan_group_l_bytes,
           stats.refactor_supernode_consumer_plan_group_l_state_built,
           stats.refactor_supernode_consumer_plan_group_l_state_storage_limited,
           stats.refactor_supernode_consumer_plan_group_l_state_run_count,
           stats.refactor_supernode_consumer_plan_group_l_state_rows,
           stats.refactor_supernode_consumer_plan_group_l_state_max_rows,
           stats.refactor_supernode_consumer_plan_group_l_state_bytes,
           stats.refactor_supernode_consumer_plan_group_l_state_focus_enabled,
           stats
             .refactor_supernode_consumer_plan_group_l_state_candidate_group_count,
           stats
             .refactor_supernode_consumer_plan_group_l_state_candidate_run_count,
           stats.refactor_supernode_consumer_plan_group_l_state_candidate_rows,
           stats
             .refactor_supernode_consumer_plan_group_l_state_selected_group_count,
           stats
             .refactor_supernode_consumer_plan_group_l_state_selected_run_count,
           stats.refactor_supernode_consumer_plan_group_l_state_selected_rows,
           stats.refactor_supernode_consumer_plan_group_l_state_selected_bytes,
           stats.refactor_last_supernode_consumer_plan_group_l_dense_writes,
           stats.refactor_last_supernode_consumer_plan_group_l_trailing_writes,
           stats.refactor_last_supernode_consumer_plan_group_l_invalidations,
           stats.refactor_supernode_consumer_plan_group_l_dense_write_count,
           stats.refactor_supernode_consumer_plan_group_l_trailing_write_count,
           stats.refactor_supernode_consumer_plan_group_l_invalidation_count,
           stats.refactor_last_supernode_consumer_plan_group_l_update_runs,
           stats.refactor_last_supernode_consumer_plan_group_l_update_rows,
           stats.refactor_last_supernode_consumer_plan_group_l_update_entries,
           stats.refactor_supernode_consumer_plan_group_l_update_run_count,
           stats.refactor_supernode_consumer_plan_group_l_update_rows,
           stats.refactor_supernode_consumer_plan_group_l_update_entries,
           stats.refactor_supernode_consumer_plan_deferred_columns,
           stats.refactor_supernode_consumer_plan_deferred_entries,
           stats.refactor_supernode_consumer_plan_deferred_unique_rows,
           stats.refactor_supernode_consumer_plan_batch_deferred_columns,
           stats.refactor_supernode_consumer_plan_batch_deferred_entries,
           stats.refactor_supernode_consumer_plan_batch_deferred_unique_rows,
           stats.row_refactor_group_count,
           stats.row_refactor_group_single_count,
           stats.row_refactor_group_batch_count,
           stats.row_refactor_group_batch_rows,
           stats.row_refactor_group_batch_max_width,
           stats.row_refactor_group_batch_width_le_4_count,
           stats.row_refactor_group_batch_width_le_8_count,
           stats.row_refactor_group_scalar_candidate_count,
           stats.row_refactor_group_scalar_candidate_rows,
           stats.row_refactor_group_scalar_short_count,
           stats.row_refactor_group_scalar_short_rows,
           stats.row_refactor_group_scalar_stop_level_mismatch_count,
           stats.row_refactor_group_scalar_stop_internal_dep_count,
           stats.row_refactor_group_scalar_stop_next_segment_count,
           stats.row_refactor_group_scalar_stop_max_width_count,
           stats.row_refactor_group_scalar_stop_matrix_end_count,
           stats.row_refactor_group_generic_count,
           stats.row_refactor_group_generic_rows,
           stats.row_refactor_group_generic_max_width,
           stats.row_refactor_group_dense_count,
           stats.row_refactor_group_dense_rows,
           stats.row_refactor_group_dense_max_width,
           stats.row_refactor_group_single_work,
           stats.row_refactor_group_batch_work,
           stats.row_refactor_group_generic_work,
           stats.row_refactor_group_dense_work,
           stats.row_refactor_group_level_count,
           stats.row_refactor_group_level_max_width,
           stats.row_refactor_group_cluster_levels,
           stats.row_refactor_group_pipeline_groups,
           stats.row_refactor_group_pipeline_rows,
           stats.row_refactor_group_pipeline_work,
           stats.row_refactor_total_group_work,
           stats.row_refactor_auto_enabled,
           stats.row_refactor_auto_values_ready,
           stats.row_refactor_auto_work_allowed,
           stats.row_refactor_auto_should_run,
           stats.row_refactor_auto_lower_bound_work,
           stats.row_refactor_auto_lower_bound_rejected,
           stats.row_refactor_auto_pattern_build_failed,
           stats.row_refactor_auto_value_copy_failed,
           stats.row_refactor_auto_model_recommended,
           stats.row_refactor_auto_model_attempted,
           stats.row_refactor_auto_model_accepted,
           stats.row_refactor_group_dependency_edges,
           stats.row_refactor_group_root_count,
           stats.row_refactor_group_leaf_count,
           stats.row_refactor_group_max_fanout,
           stats.row_refactor_last_run,
           stats.row_refactor_last_checked,
           stats.row_refactor_last_parallel,
           stats.row_refactor_last_ready_queue,
           stats.row_refactor_run_count,
           stats.row_refactor_checked_run_count,
           stats.row_refactor_parallel_run_count,
           stats.row_refactor_ready_queue_run_count,
           stats.row_refactor_ready_queue_group_count,
           stats.row_refactor_last_done_bitmap,
           stats.row_refactor_done_bitmap_run_count,
           stats.row_refactor_last_prefactor,
           stats.row_refactor_last_prefactor_rows,
           stats.row_refactor_last_prefactor_deps,
           stats.row_refactor_prefactor_run_count,
           stats.row_refactor_prefactor_rows,
           stats.row_refactor_prefactor_deps,
           stats.row_refactor_last_prefactor_supernode,
           stats.row_refactor_last_prefactor_supernode_rows,
           stats.row_refactor_last_prefactor_supernode_deps,
           stats.row_refactor_prefactor_supernode_run_count,
           stats.row_refactor_prefactor_supernode_rows,
           stats.row_refactor_prefactor_supernode_deps,
           stats.row_refactor_input_cleanup_rows,
           stats.row_refactor_input_cleanup_entries,
           stats.row_refactor_last_defer_value_scatter,
           stats.row_refactor_defer_value_scatter_run_count,
           stats.row_refactor_values_dirty,
           stats.row_refactor_last_lazy_value_scatter,
           stats.row_refactor_lazy_value_scatter_run_count,
           stats.row_refactor_last_row_solve,
           stats.row_refactor_row_solve_run_count);
    printf(",\"refactor_supernode_consumer_panel_count\":%" PRId64
           ",\"refactor_supernode_consumer_reused_panel_count\":%" PRId64
           ",\"refactor_supernode_consumer_reused_run_count\":%" PRId64
           ",\"refactor_supernode_consumer_reused_run_rows\":%" PRId64
           ",\"refactor_supernode_consumer_panel_max_runs\":%" PRId64
           ",\"refactor_supernode_consumer_panel_max_rows\":%" PRId64
           ",\"refactor_supernode_consumer_reused_l_entries\":%.9g"
           ",\"refactor_supernode_consumer_reused_internal_entries\":%.9g",
           stats.refactor_supernode_consumer_panel_count,
           stats.refactor_supernode_consumer_reused_panel_count,
           stats.refactor_supernode_consumer_reused_run_count,
           stats.refactor_supernode_consumer_reused_run_rows,
           stats.refactor_supernode_consumer_panel_max_runs,
           stats.refactor_supernode_consumer_panel_max_rows,
           stats.refactor_supernode_consumer_reused_l_entries,
           stats.refactor_supernode_consumer_reused_internal_entries);
    printf(",\"refactor_supernode_consumer_plan_panel_count\":%" PRId64
           ",\"refactor_supernode_consumer_plan_reused_panel_count\":%" PRId64
           ",\"refactor_supernode_consumer_plan_run_count\":%" PRId64
           ",\"refactor_supernode_consumer_plan_run_rows\":%" PRId64
           ",\"refactor_supernode_consumer_plan_positioned_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_positioned_run_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_max_panel_runs\":%" PRId64
           ",\"refactor_supernode_consumer_plan_max_panel_rows\":%" PRId64
           ",\"refactor_supernode_consumer_plan_l_entries\":%.9g"
           ",\"refactor_supernode_consumer_plan_internal_entries\":%.9g"
           ",\"refactor_supernode_consumer_plan_bytes\":%" PRId64
           ",\"refactor_supernode_consumer_plan_small_run_count\":%" PRId64
           ",\"refactor_supernode_consumer_plan_small_run_rows\":%" PRId64
           ",\"refactor_supernode_consumer_plan_batch_panel_count\":%" PRId64
           ",\"refactor_supernode_consumer_plan_batch_run_count\":%" PRId64
           ",\"refactor_supernode_consumer_plan_batch_run_rows\":%" PRId64
           ",\"refactor_supernode_consumer_plan_batch_small_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_batch_small_run_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_column_count\":%" PRId64
           ",\"refactor_supernode_consumer_plan_max_column_runs\":%" PRId64
           ",\"refactor_supernode_consumer_plan_max_column_rows\":%" PRId64
           ",\"refactor_supernode_consumer_plan_column_batch_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_column_batch_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_column_batch_run_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_column_batch_small_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_column_batch_small_run_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_column_shape_batch_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_column_shape_batch_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_column_shape_batch_run_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_column_shape_batch_small_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_column_shape_batch_small_run_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_column_shape_batch_max_runs\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_shape_batch_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_shape_batch_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_shape_batch_run_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_shape_batch_small_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_shape_batch_small_run_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_shape_batch_max_runs\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_built\":%d"
           ",\"refactor_supernode_consumer_plan_group_l_storage_limited\":%d"
           ",\"refactor_supernode_consumer_plan_group_l_panel_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_rows\":%" PRId64
           ",\"refactor_supernode_consumer_plan_group_l_run_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_exec_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_exec_run_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_batch_candidate_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_batch_candidate_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_batch_candidate_run_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_batch_candidate_entries\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_batch_candidate_max_runs\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_batch_candidate_advance_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_batch_candidate_advance_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_batch_candidate_advance_run_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_batch_candidate_advance_work\":%.9g"
           ",\"refactor_supernode_consumer_plan_group_l_batch_candidate_advance_max_work\":%.9g"
           ",\"refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_run_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_entries\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_advance_work\":%.9g"
           ",\"refactor_supernode_consumer_plan_group_l_batch_candidate_max_payoff_ratio\":%.9g"
           ",\"refactor_supernode_consumer_plan_group_l_dense_entries\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_trailing_entries\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_bytes\":%" PRId64
           ",\"refactor_supernode_consumer_plan_group_l_state_built\":%d"
           ",\"refactor_supernode_consumer_plan_group_l_state_storage_limited\":%d"
           ",\"refactor_supernode_consumer_plan_group_l_state_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_state_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_state_max_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_state_bytes\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_state_focus_enabled\":%d"
           ",\"refactor_supernode_consumer_plan_group_l_state_candidate_group_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_state_candidate_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_state_candidate_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_state_selected_group_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_state_selected_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_state_selected_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_state_selected_bytes\":%"
           PRId64
           ",\"refactor_last_supernode_consumer_plan_group_l_dense_writes\":%"
           PRId64
           ",\"refactor_last_supernode_consumer_plan_group_l_trailing_writes\":%"
           PRId64
           ",\"refactor_last_supernode_consumer_plan_group_l_invalidations\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_dense_write_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_trailing_write_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_invalidation_count\":%"
           PRId64
           ",\"refactor_last_supernode_consumer_plan_group_l_update_runs\":%"
           PRId64
           ",\"refactor_last_supernode_consumer_plan_group_l_update_rows\":%"
           PRId64
           ",\"refactor_last_supernode_consumer_plan_group_l_update_entries\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_update_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_update_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_group_l_update_entries\":%"
           PRId64
           ",\"refactor_last_supernode_consumer_plan_attempts\":%" PRId64
           ",\"refactor_last_supernode_consumer_plan_hits\":%" PRId64
           ",\"refactor_last_supernode_consumer_plan_applied\":%" PRId64
           ",\"refactor_last_supernode_consumer_plan_rows\":%" PRId64
           ",\"refactor_last_supernode_consumer_plan_entries\":%" PRId64
           ",\"refactor_supernode_consumer_plan_attempt_count\":%" PRId64
           ",\"refactor_supernode_consumer_plan_hit_count\":%" PRId64
           ",\"refactor_supernode_consumer_plan_apply_count\":%" PRId64
           ",\"refactor_supernode_consumer_plan_apply_rows\":%" PRId64
           ",\"refactor_supernode_consumer_plan_apply_entries\":%" PRId64
           ",\"refactor_supernode_consumer_plan_exec_disabled\":%d"
           ",\"refactor_supernode_consumer_plan_exec_disable_count\":%" PRId64
           ",\"refactor_last_supernode_consumer_plan_claimed_columns\":%"
           PRId64
           ",\"refactor_last_supernode_consumer_plan_claim_skips\":%" PRId64
           ",\"refactor_last_supernode_consumer_plan_claim_waits\":%" PRId64
           ",\"refactor_supernode_consumer_plan_claimed_columns\":%" PRId64
           ",\"refactor_supernode_consumer_plan_claim_skip_count\":%" PRId64
           ",\"refactor_supernode_consumer_plan_claim_wait_count\":%" PRId64,
           stats.refactor_supernode_consumer_plan_panel_count,
           stats.refactor_supernode_consumer_plan_reused_panel_count,
           stats.refactor_supernode_consumer_plan_run_count,
           stats.refactor_supernode_consumer_plan_run_rows,
           stats.refactor_supernode_consumer_plan_positioned_run_count,
           stats.refactor_supernode_consumer_plan_positioned_run_rows,
           stats.refactor_supernode_consumer_plan_max_panel_runs,
           stats.refactor_supernode_consumer_plan_max_panel_rows,
           stats.refactor_supernode_consumer_plan_l_entries,
           stats.refactor_supernode_consumer_plan_internal_entries,
           stats.refactor_supernode_consumer_plan_bytes,
           stats.refactor_supernode_consumer_plan_small_run_count,
           stats.refactor_supernode_consumer_plan_small_run_rows,
           stats.refactor_supernode_consumer_plan_batch_panel_count,
           stats.refactor_supernode_consumer_plan_batch_run_count,
           stats.refactor_supernode_consumer_plan_batch_run_rows,
           stats.refactor_supernode_consumer_plan_batch_small_run_count,
           stats.refactor_supernode_consumer_plan_batch_small_run_rows,
           stats.refactor_supernode_consumer_plan_column_count,
           stats.refactor_supernode_consumer_plan_max_column_runs,
           stats.refactor_supernode_consumer_plan_max_column_rows,
           stats.refactor_supernode_consumer_plan_column_batch_count,
           stats.refactor_supernode_consumer_plan_column_batch_run_count,
           stats.refactor_supernode_consumer_plan_column_batch_run_rows,
           stats.refactor_supernode_consumer_plan_column_batch_small_run_count,
           stats.refactor_supernode_consumer_plan_column_batch_small_run_rows,
           stats.refactor_supernode_consumer_plan_column_shape_batch_count,
           stats.refactor_supernode_consumer_plan_column_shape_batch_run_count,
           stats.refactor_supernode_consumer_plan_column_shape_batch_run_rows,
           stats.refactor_supernode_consumer_plan_column_shape_batch_small_run_count,
           stats.refactor_supernode_consumer_plan_column_shape_batch_small_run_rows,
           stats.refactor_supernode_consumer_plan_column_shape_batch_max_runs,
           stats.refactor_supernode_consumer_plan_shape_batch_count,
           stats.refactor_supernode_consumer_plan_shape_batch_run_count,
           stats.refactor_supernode_consumer_plan_shape_batch_run_rows,
           stats.refactor_supernode_consumer_plan_shape_batch_small_run_count,
           stats.refactor_supernode_consumer_plan_shape_batch_small_run_rows,
           stats.refactor_supernode_consumer_plan_shape_batch_max_runs,
           stats.refactor_supernode_consumer_plan_group_l_built,
           stats.refactor_supernode_consumer_plan_group_l_storage_limited,
           stats.refactor_supernode_consumer_plan_group_l_panel_count,
           stats.refactor_supernode_consumer_plan_group_l_run_count,
           stats.refactor_supernode_consumer_plan_group_l_rows,
           stats.refactor_supernode_consumer_plan_group_l_run_rows,
           stats.refactor_supernode_consumer_plan_group_l_exec_run_count,
           stats.refactor_supernode_consumer_plan_group_l_exec_run_rows,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_count,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_run_count,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_run_rows,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_entries,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_max_runs,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_advance_count,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_advance_run_count,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_advance_run_rows,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_advance_work,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_advance_max_work,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_count,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_run_count,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_run_rows,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_entries,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_advance_work,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_max_payoff_ratio,
           stats.refactor_supernode_consumer_plan_group_l_dense_entries,
           stats.refactor_supernode_consumer_plan_group_l_trailing_entries,
           stats.refactor_supernode_consumer_plan_group_l_bytes,
           stats.refactor_supernode_consumer_plan_group_l_state_built,
           stats.refactor_supernode_consumer_plan_group_l_state_storage_limited,
           stats.refactor_supernode_consumer_plan_group_l_state_run_count,
           stats.refactor_supernode_consumer_plan_group_l_state_rows,
           stats.refactor_supernode_consumer_plan_group_l_state_max_rows,
           stats.refactor_supernode_consumer_plan_group_l_state_bytes,
           stats.refactor_supernode_consumer_plan_group_l_state_focus_enabled,
           stats
             .refactor_supernode_consumer_plan_group_l_state_candidate_group_count,
           stats
             .refactor_supernode_consumer_plan_group_l_state_candidate_run_count,
           stats.refactor_supernode_consumer_plan_group_l_state_candidate_rows,
           stats
             .refactor_supernode_consumer_plan_group_l_state_selected_group_count,
           stats
             .refactor_supernode_consumer_plan_group_l_state_selected_run_count,
           stats.refactor_supernode_consumer_plan_group_l_state_selected_rows,
           stats.refactor_supernode_consumer_plan_group_l_state_selected_bytes,
           stats.refactor_last_supernode_consumer_plan_group_l_dense_writes,
           stats.refactor_last_supernode_consumer_plan_group_l_trailing_writes,
           stats.refactor_last_supernode_consumer_plan_group_l_invalidations,
           stats.refactor_supernode_consumer_plan_group_l_dense_write_count,
           stats.refactor_supernode_consumer_plan_group_l_trailing_write_count,
           stats.refactor_supernode_consumer_plan_group_l_invalidation_count,
           stats.refactor_last_supernode_consumer_plan_group_l_update_runs,
           stats.refactor_last_supernode_consumer_plan_group_l_update_rows,
           stats.refactor_last_supernode_consumer_plan_group_l_update_entries,
           stats.refactor_supernode_consumer_plan_group_l_update_run_count,
           stats.refactor_supernode_consumer_plan_group_l_update_rows,
           stats.refactor_supernode_consumer_plan_group_l_update_entries,
           stats.refactor_last_supernode_consumer_plan_attempts,
           stats.refactor_last_supernode_consumer_plan_hits,
           stats.refactor_last_supernode_consumer_plan_applied,
           stats.refactor_last_supernode_consumer_plan_rows,
           stats.refactor_last_supernode_consumer_plan_entries,
           stats.refactor_supernode_consumer_plan_attempt_count,
           stats.refactor_supernode_consumer_plan_hit_count,
           stats.refactor_supernode_consumer_plan_apply_count,
           stats.refactor_supernode_consumer_plan_apply_rows,
           stats.refactor_supernode_consumer_plan_apply_entries,
           stats.refactor_supernode_consumer_plan_exec_disabled,
           stats.refactor_supernode_consumer_plan_exec_disable_count,
           stats.refactor_last_supernode_consumer_plan_claimed_columns,
           stats.refactor_last_supernode_consumer_plan_claim_skips,
           stats.refactor_last_supernode_consumer_plan_claim_waits,
           stats.refactor_supernode_consumer_plan_claimed_columns,
           stats.refactor_supernode_consumer_plan_claim_skip_count,
           stats.refactor_supernode_consumer_plan_claim_wait_count);
    printf(",\"refactor_supernode_consumer_plan_shape_targets_built\":%d"
           ",\"refactor_supernode_consumer_plan_shape_target_group_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_shape_target_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_shape_target_entries\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_shape_target_rows\":%" PRId64
           ",\"refactor_supernode_consumer_plan_shape_target_max_rows\":%"
           PRId64,
           stats.refactor_supernode_consumer_plan_shape_targets_built,
           stats.refactor_supernode_consumer_plan_shape_target_group_count,
           stats.refactor_supernode_consumer_plan_shape_target_run_count,
           stats.refactor_supernode_consumer_plan_shape_target_entries,
           stats.refactor_supernode_consumer_plan_shape_target_rows,
           stats.refactor_supernode_consumer_plan_shape_target_max_rows);
    printf(",\"refactor_supernode_consumer_plan_first_dep_shape_batch_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_first_dep_shape_batch_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_first_dep_shape_batch_run_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_first_dep_shape_batch_small_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_first_dep_shape_batch_small_run_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_first_dep_shape_batch_max_runs\":%"
           PRId64,
           stats.refactor_supernode_consumer_plan_first_dep_shape_batch_count,
           stats
             .refactor_supernode_consumer_plan_first_dep_shape_batch_run_count,
           stats
             .refactor_supernode_consumer_plan_first_dep_shape_batch_run_rows,
           stats
             .refactor_supernode_consumer_plan_first_dep_shape_batch_small_run_count,
           stats
             .refactor_supernode_consumer_plan_first_dep_shape_batch_small_run_rows,
           stats
             .refactor_supernode_consumer_plan_first_dep_shape_batch_max_runs);
    printf(",\"refactor_supernode_consumer_plan_shape_batch_advance_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_shape_batch_advance_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_shape_batch_advance_run_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_shape_batch_advance_dep_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_shape_batch_advance_work\":%.9g"
           ",\"refactor_supernode_consumer_plan_shape_batch_advance_max_deps\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_shape_batch_advance_max_work\":%.9g",
           stats.refactor_supernode_consumer_plan_shape_batch_advance_count,
           stats
             .refactor_supernode_consumer_plan_shape_batch_advance_run_count,
           stats
             .refactor_supernode_consumer_plan_shape_batch_advance_run_rows,
           stats
             .refactor_supernode_consumer_plan_shape_batch_advance_dep_count,
           stats.refactor_supernode_consumer_plan_shape_batch_advance_work,
           stats
             .refactor_supernode_consumer_plan_shape_batch_advance_max_deps,
           stats.refactor_supernode_consumer_plan_shape_batch_advance_max_work);
    printf(",\"refactor_supernode_consumer_plan_prefix_advance_batch_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_prefix_advance_batch_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_prefix_advance_batch_run_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_prefix_advance_batch_dep_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_prefix_advance_batch_work\":%.9g"
           ",\"refactor_supernode_consumer_plan_prefix_advance_batch_max_runs\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_prefix_advance_batch_max_deps\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_prefix_advance_batch_max_work\":%.9g",
           stats.refactor_supernode_consumer_plan_prefix_advance_batch_count,
           stats
             .refactor_supernode_consumer_plan_prefix_advance_batch_run_count,
           stats
             .refactor_supernode_consumer_plan_prefix_advance_batch_run_rows,
           stats
             .refactor_supernode_consumer_plan_prefix_advance_batch_dep_count,
           stats.refactor_supernode_consumer_plan_prefix_advance_batch_work,
           stats
             .refactor_supernode_consumer_plan_prefix_advance_batch_max_runs,
           stats
             .refactor_supernode_consumer_plan_prefix_advance_batch_max_deps,
           stats
             .refactor_supernode_consumer_plan_prefix_advance_batch_max_work);
    printf(",\"refactor_supernode_consumer_plan_shape_bounded_advance_dep_limit\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_shape_bounded_advance_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_shape_bounded_advance_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_shape_bounded_advance_run_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_shape_bounded_advance_dep_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_shape_bounded_advance_update_entries\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_shape_bounded_advance_work\":%.9g"
           ",\"refactor_supernode_consumer_plan_shape_bounded_advance_max_run_deps\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_shape_bounded_advance_max_run_work\":%.9g"
           ",\"refactor_supernode_consumer_plan_shape_bounded_advance_payoff_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_shape_bounded_advance_payoff_run_count\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_shape_bounded_advance_payoff_run_rows\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_shape_bounded_advance_payoff_update_entries\":%"
           PRId64
           ",\"refactor_supernode_consumer_plan_shape_bounded_advance_payoff_advance_work\":%.9g"
           ",\"refactor_supernode_consumer_plan_shape_bounded_advance_max_payoff_ratio\":%.9g",
           stats
             .refactor_supernode_consumer_plan_shape_bounded_advance_dep_limit,
           stats.refactor_supernode_consumer_plan_shape_bounded_advance_count,
           stats
             .refactor_supernode_consumer_plan_shape_bounded_advance_run_count,
           stats
             .refactor_supernode_consumer_plan_shape_bounded_advance_run_rows,
           stats
             .refactor_supernode_consumer_plan_shape_bounded_advance_dep_count,
           stats
             .refactor_supernode_consumer_plan_shape_bounded_advance_update_entries,
           stats.refactor_supernode_consumer_plan_shape_bounded_advance_work,
           stats
             .refactor_supernode_consumer_plan_shape_bounded_advance_max_run_deps,
           stats
             .refactor_supernode_consumer_plan_shape_bounded_advance_max_run_work,
           stats
             .refactor_supernode_consumer_plan_shape_bounded_advance_payoff_count,
           stats
             .refactor_supernode_consumer_plan_shape_bounded_advance_payoff_run_count,
           stats
             .refactor_supernode_consumer_plan_shape_bounded_advance_payoff_run_rows,
           stats
             .refactor_supernode_consumer_plan_shape_bounded_advance_payoff_update_entries,
           stats
             .refactor_supernode_consumer_plan_shape_bounded_advance_payoff_advance_work,
           stats
             .refactor_supernode_consumer_plan_shape_bounded_advance_max_payoff_ratio);
    printf(",\"refactor_supernode_algorithm5_large_panel_count\":%" PRId64
           ",\"refactor_supernode_algorithm5_large_panel_rows\":%" PRId64
           ",\"refactor_supernode_algorithm5_large_panel_prefix_rows\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_large_panel_max_width\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_candidate_run_count\":%" PRId64
           ",\"refactor_supernode_algorithm5_candidate_run_rows\":%" PRId64
           ",\"refactor_supernode_algorithm5_prefix_run_count\":%" PRId64
           ",\"refactor_supernode_algorithm5_prefix_run_rows\":%" PRId64
           ",\"refactor_supernode_algorithm5_prefix_update_work\":%.9g"
           ",\"refactor_supernode_algorithm5_crossing_run_count\":%" PRId64
           ",\"refactor_supernode_algorithm5_crossing_run_rows\":%" PRId64,
           stats.refactor_supernode_algorithm5_large_panel_count,
           stats.refactor_supernode_algorithm5_large_panel_rows,
           stats.refactor_supernode_algorithm5_large_panel_prefix_rows,
           stats.refactor_supernode_algorithm5_large_panel_max_width,
           stats.refactor_supernode_algorithm5_candidate_run_count,
           stats.refactor_supernode_algorithm5_candidate_run_rows,
           stats.refactor_supernode_algorithm5_prefix_run_count,
           stats.refactor_supernode_algorithm5_prefix_run_rows,
           stats.refactor_supernode_algorithm5_prefix_update_work,
           stats.refactor_supernode_algorithm5_crossing_run_count,
           stats.refactor_supernode_algorithm5_crossing_run_rows);
    printf(",\"refactor_supernode_algorithm5_prefix_advance_run_count\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_prefix_advance_run_rows\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_prefix_advance_dep_count\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_prefix_advance_work\":%.9g"
           ",\"refactor_supernode_algorithm5_prefix_payoff_run_count\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_prefix_payoff_run_rows\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_prefix_payoff_update_work\":%.9g"
           ",\"refactor_supernode_algorithm5_prefix_payoff_advance_work\":%.9g"
           ",\"refactor_last_egraph_algorithm5_prefactor_columns\":%" PRId64
           ",\"refactor_last_egraph_algorithm5_prefactor_deps\":%" PRId64
           ",\"refactor_egraph_algorithm5_prefactor_column_count\":%" PRId64
           ",\"refactor_egraph_algorithm5_prefactor_dep_count\":%" PRId64
           ",\"refactor_supernode_algorithm5_prefix_panel_count\":%" PRId64
           ",\"refactor_supernode_algorithm5_prefix_panel_run_count\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_prefix_panel_run_rows\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_prefix_panel_advance_dep_count\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_prefix_panel_update_work\":%.9g"
           ",\"refactor_supernode_algorithm5_prefix_panel_advance_work\":%.9g"
           ",\"refactor_supernode_algorithm5_prefix_panel_max_runs\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_prefix_panel_payoff_count\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_prefix_panel_payoff_run_count\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_prefix_panel_payoff_run_rows\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_prefix_panel_payoff_update_work\":%.9g"
           ",\"refactor_supernode_algorithm5_prefix_panel_payoff_advance_work\":%.9g"
           ",\"refactor_supernode_algorithm5_prefix_panel_payoff_subset_count\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_prefix_panel_payoff_subset_run_count\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_prefix_panel_payoff_subset_run_rows\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_prefix_panel_payoff_subset_update_work\":%.9g"
           ",\"refactor_supernode_algorithm5_prefix_panel_payoff_subset_advance_work\":%.9g"
           ",\"refactor_supernode_algorithm5_prefix_panel_payoff_subset_max_runs\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_prefix_panel_payoff_subset_current_count\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_prefix_panel_payoff_subset_multi_current_count\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_prefix_panel_payoff_subset_max_currents\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_count\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_prefix_rows\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_max_run_rows\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_positioned_runs\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_current_total\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_multi_current_count\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_max_currents\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_workspace_rows\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_max_workspace_rows\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_advance_deps\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_max_advance_deps\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_zero_advance_runs\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_advance_slots\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_max_advance_slots\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_max_run_advance_slots\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_target_entries\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_max_target_entries\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_max_run_target_entries\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_target_slots\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_max_target_slots\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_max_run_target_slots\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_pattern_width\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_max_pattern_width\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_current_state_rows\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_current_state_max_rows\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_current_state_span_rows\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_current_state_max_span_rows\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_runtime_workspace_rows\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_runtime_target_slots\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_runtime_current_count\":%"
           PRId64,
           stats.refactor_supernode_algorithm5_prefix_advance_run_count,
           stats.refactor_supernode_algorithm5_prefix_advance_run_rows,
           stats.refactor_supernode_algorithm5_prefix_advance_dep_count,
           stats.refactor_supernode_algorithm5_prefix_advance_work,
           stats.refactor_supernode_algorithm5_prefix_payoff_run_count,
           stats.refactor_supernode_algorithm5_prefix_payoff_run_rows,
           stats.refactor_supernode_algorithm5_prefix_payoff_update_work,
           stats.refactor_supernode_algorithm5_prefix_payoff_advance_work,
           stats.refactor_last_egraph_algorithm5_prefactor_columns,
           stats.refactor_last_egraph_algorithm5_prefactor_deps,
           stats.refactor_egraph_algorithm5_prefactor_column_count,
           stats.refactor_egraph_algorithm5_prefactor_dep_count,
           stats.refactor_supernode_algorithm5_prefix_panel_count,
           stats.refactor_supernode_algorithm5_prefix_panel_run_count,
           stats.refactor_supernode_algorithm5_prefix_panel_run_rows,
           stats
             .refactor_supernode_algorithm5_prefix_panel_advance_dep_count,
           stats.refactor_supernode_algorithm5_prefix_panel_update_work,
           stats.refactor_supernode_algorithm5_prefix_panel_advance_work,
           stats.refactor_supernode_algorithm5_prefix_panel_max_runs,
           stats.refactor_supernode_algorithm5_prefix_panel_payoff_count,
           stats.refactor_supernode_algorithm5_prefix_panel_payoff_run_count,
           stats.refactor_supernode_algorithm5_prefix_panel_payoff_run_rows,
           stats
             .refactor_supernode_algorithm5_prefix_panel_payoff_update_work,
           stats
             .refactor_supernode_algorithm5_prefix_panel_payoff_advance_work,
           stats
             .refactor_supernode_algorithm5_prefix_panel_payoff_subset_count,
           stats
             .refactor_supernode_algorithm5_prefix_panel_payoff_subset_run_count,
           stats
             .refactor_supernode_algorithm5_prefix_panel_payoff_subset_run_rows,
           stats
             .refactor_supernode_algorithm5_prefix_panel_payoff_subset_update_work,
           stats
             .refactor_supernode_algorithm5_prefix_panel_payoff_subset_advance_work,
           stats
             .refactor_supernode_algorithm5_prefix_panel_payoff_subset_max_runs,
           stats
             .refactor_supernode_algorithm5_prefix_panel_payoff_subset_current_count,
           stats
             .refactor_supernode_algorithm5_prefix_panel_payoff_subset_multi_current_count,
           stats
             .refactor_supernode_algorithm5_prefix_panel_payoff_subset_max_currents,
           stats.refactor_supernode_algorithm5_payoff_group_count,
           stats.refactor_supernode_algorithm5_payoff_group_prefix_rows,
           stats.refactor_supernode_algorithm5_payoff_group_max_run_rows,
           stats.refactor_supernode_algorithm5_payoff_group_positioned_runs,
           stats.refactor_supernode_algorithm5_payoff_group_current_total,
           stats
             .refactor_supernode_algorithm5_payoff_group_multi_current_count,
           stats.refactor_supernode_algorithm5_payoff_group_max_currents,
           stats.refactor_supernode_algorithm5_payoff_group_workspace_rows,
           stats
             .refactor_supernode_algorithm5_payoff_group_max_workspace_rows,
           stats.refactor_supernode_algorithm5_payoff_group_advance_deps,
           stats
             .refactor_supernode_algorithm5_payoff_group_max_advance_deps,
           stats
             .refactor_supernode_algorithm5_payoff_group_zero_advance_runs,
           stats.refactor_supernode_algorithm5_payoff_group_advance_slots,
           stats
             .refactor_supernode_algorithm5_payoff_group_max_advance_slots,
           stats
             .refactor_supernode_algorithm5_payoff_group_max_run_advance_slots,
           stats.refactor_supernode_algorithm5_payoff_group_target_entries,
           stats
             .refactor_supernode_algorithm5_payoff_group_max_target_entries,
           stats
             .refactor_supernode_algorithm5_payoff_group_max_run_target_entries,
           stats.refactor_supernode_algorithm5_payoff_group_target_slots,
           stats
             .refactor_supernode_algorithm5_payoff_group_max_target_slots,
           stats
             .refactor_supernode_algorithm5_payoff_group_max_run_target_slots,
           stats.refactor_supernode_algorithm5_payoff_group_pattern_width,
           stats
             .refactor_supernode_algorithm5_payoff_group_max_pattern_width,
           stats.refactor_supernode_algorithm5_payoff_current_state_rows,
           stats
             .refactor_supernode_algorithm5_payoff_current_state_max_rows,
           stats
             .refactor_supernode_algorithm5_payoff_current_state_span_rows,
           stats
             .refactor_supernode_algorithm5_payoff_current_state_max_span_rows,
           stats.refactor_supernode_algorithm5_payoff_runtime_workspace_rows,
           stats.refactor_supernode_algorithm5_payoff_runtime_target_slots,
           stats.refactor_supernode_algorithm5_payoff_runtime_current_count);
    printf(",\"refactor_supernode_algorithm5_payoff_group_advance_unique_deps\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_advance_duplicate_deps\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_advance_shared_deps\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_advance_max_dep_fanout\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_advance_duplicate_update_entries\":%"
           PRId64,
           stats
             .refactor_supernode_algorithm5_payoff_group_advance_unique_deps,
           stats
             .refactor_supernode_algorithm5_payoff_group_advance_duplicate_deps,
           stats
             .refactor_supernode_algorithm5_payoff_group_advance_shared_deps,
           stats
             .refactor_supernode_algorithm5_payoff_group_advance_max_dep_fanout,
           stats
             .refactor_supernode_algorithm5_payoff_group_advance_duplicate_update_entries);
    printf(",\"refactor_supernode_algorithm5_payoff_group_suffix_deps\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_max_suffix_deps\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_suffix_update_entries\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_max_run_suffix_update_entries\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_suffix_unique_deps\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_suffix_duplicate_deps\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_suffix_shared_deps\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_suffix_max_dep_fanout\":%"
           PRId64
           ",\"refactor_supernode_algorithm5_payoff_group_suffix_duplicate_update_entries\":%"
           PRId64,
           stats.refactor_supernode_algorithm5_payoff_group_suffix_deps,
           stats
             .refactor_supernode_algorithm5_payoff_group_max_suffix_deps,
           stats
             .refactor_supernode_algorithm5_payoff_group_suffix_update_entries,
           stats
             .refactor_supernode_algorithm5_payoff_group_max_run_suffix_update_entries,
           stats
             .refactor_supernode_algorithm5_payoff_group_suffix_unique_deps,
           stats
             .refactor_supernode_algorithm5_payoff_group_suffix_duplicate_deps,
           stats
             .refactor_supernode_algorithm5_payoff_group_suffix_shared_deps,
           stats
             .refactor_supernode_algorithm5_payoff_group_suffix_max_dep_fanout,
           stats
             .refactor_supernode_algorithm5_payoff_group_suffix_duplicate_update_entries);
    printf(",\"refactor_last_supernode_algorithm5_payoff_slot_accum_runs\":%"
           PRId64
           ",\"refactor_last_supernode_algorithm5_payoff_slot_accum_rows\":%"
           PRId64
           ",\"refactor_last_supernode_algorithm5_payoff_slot_accum_target_entries\":%"
           PRId64
           ",\"refactor_last_supernode_algorithm5_payoff_slot_accum_target_slots\":%"
           PRId64
           ",\"refactor_last_supernode_algorithm5_payoff_prefix_prep_runs\":%"
           PRId64
           ",\"refactor_last_supernode_algorithm5_payoff_prefix_prep_rows\":%"
           PRId64
           ",\"refactor_last_supernode_algorithm5_payoff_prefix_prep_target_entries\":%"
           PRId64
           ",\"refactor_last_supernode_algorithm5_payoff_prefix_prep_target_slots\":%"
           PRId64
           ",\"refactor_last_supernode_algorithm5_payoff_advance_seed_runs\":%"
           PRId64
           ",\"refactor_last_supernode_algorithm5_payoff_advance_seed_deps\":%"
           PRId64
           ",\"refactor_last_supernode_algorithm5_payoff_advance_seed_slots\":%"
           PRId64
           ",\"refactor_last_supernode_algorithm5_payoff_current_state_seed_runs\":%"
           PRId64
           ",\"refactor_last_supernode_algorithm5_payoff_current_state_seed_deps\":%"
           PRId64
           ",\"refactor_last_supernode_algorithm5_payoff_current_state_seed_rows\":%"
           PRId64,
           stats.refactor_last_supernode_algorithm5_payoff_slot_accum_runs,
           stats.refactor_last_supernode_algorithm5_payoff_slot_accum_rows,
           stats
             .refactor_last_supernode_algorithm5_payoff_slot_accum_target_entries,
           stats
             .refactor_last_supernode_algorithm5_payoff_slot_accum_target_slots,
           stats.refactor_last_supernode_algorithm5_payoff_prefix_prep_runs,
           stats.refactor_last_supernode_algorithm5_payoff_prefix_prep_rows,
           stats
             .refactor_last_supernode_algorithm5_payoff_prefix_prep_target_entries,
           stats
             .refactor_last_supernode_algorithm5_payoff_prefix_prep_target_slots,
           stats.refactor_last_supernode_algorithm5_payoff_advance_seed_runs,
           stats.refactor_last_supernode_algorithm5_payoff_advance_seed_deps,
           stats.refactor_last_supernode_algorithm5_payoff_advance_seed_slots,
           stats
             .refactor_last_supernode_algorithm5_payoff_current_state_seed_runs,
           stats
             .refactor_last_supernode_algorithm5_payoff_current_state_seed_deps,
           stats
             .refactor_last_supernode_algorithm5_payoff_current_state_seed_rows);
    printf(",\"refactor_last_supernode_algorithm5_payoff_final_trigger_batches\":%"
           PRId64
           ",\"refactor_last_supernode_algorithm5_payoff_final_trigger_multi_batches\":%"
           PRId64
           ",\"refactor_last_supernode_algorithm5_payoff_final_trigger_claims\":%"
           PRId64
           ",\"refactor_last_supernode_algorithm5_payoff_final_trigger_suffix_deps\":%"
           PRId64,
           stats
             .refactor_last_supernode_algorithm5_payoff_final_trigger_batches,
           stats
             .refactor_last_supernode_algorithm5_payoff_final_trigger_multi_batches,
           stats
             .refactor_last_supernode_algorithm5_payoff_final_trigger_claims,
           stats
             .refactor_last_supernode_algorithm5_payoff_final_trigger_suffix_deps);
    printf(",\"refactor_last_supernode_algorithm5_payoff_suffix_advance_slots\":%"
           PRId64
           ",\"refactor_last_supernode_algorithm5_payoff_suffix_advance_deps\":%"
           PRId64
           ",\"refactor_last_supernode_algorithm5_payoff_suffix_advance_updates\":%"
           PRId64
           ",\"refactor_last_supernode_algorithm5_payoff_suffix_advance_finished\":%"
           PRId64,
           stats.refactor_last_supernode_algorithm5_payoff_suffix_advance_slots,
           stats.refactor_last_supernode_algorithm5_payoff_suffix_advance_deps,
           stats
             .refactor_last_supernode_algorithm5_payoff_suffix_advance_updates,
           stats
             .refactor_last_supernode_algorithm5_payoff_suffix_advance_finished);
    printf(",\"refactor_l_pattern_columns\":%" PRId64
           ",\"refactor_l_pattern_entries\":%" PRId64
           ",\"refactor_l_adjacent_run_count\":%" PRId64
           ",\"refactor_l_adjacent_run_entries\":%" PRId64
           ",\"refactor_l_adjacent_run_max_len\":%" PRId64
           ",\"refactor_l_contiguous_suffix_columns\":%" PRId64
           ",\"refactor_l_contiguous_suffix_entries\":%" PRId64
           ",\"refactor_l_contiguous_suffix_max_len\":%" PRId64
           ",\"refactor_map_index32_enabled\":%d"
           ",\"refactor_map_index32_entries\":%" PRId64
           ",\"refactor_l_index32_enabled\":%d"
           ",\"refactor_l_index32_entries\":%" PRId64
           ",\"refactor_u_index32_enabled\":%d"
           ",\"refactor_u_index32_entries\":%" PRId64,
           stats.refactor_l_pattern_columns,
           stats.refactor_l_pattern_entries,
           stats.refactor_l_adjacent_run_count,
           stats.refactor_l_adjacent_run_entries,
           stats.refactor_l_adjacent_run_max_len,
           stats.refactor_l_contiguous_suffix_columns,
           stats.refactor_l_contiguous_suffix_entries,
           stats.refactor_l_contiguous_suffix_max_len,
           stats.refactor_map_index32_enabled,
           stats.refactor_map_index32_entries,
           stats.refactor_l_index32_enabled,
           stats.refactor_l_index32_entries,
           stats.refactor_u_index32_enabled,
           stats.refactor_u_index32_entries);
    printf(",\"row_solve_parallel_run_count\":%" PRId64
           ",\"row_solve_parallel_l_slice_runs\":%" PRId64
           ",\"row_solve_parallel_u_slice_runs\":%" PRId64
           ",\"row_solve_parallel_l_sparse_level_runs\":%" PRId64
           ",\"row_solve_parallel_u_sparse_level_runs\":%" PRId64
           ",\"row_solve_thread_count\":%" PRId64
           ",\"row_solve_l_thread_max_rect_entries\":%" PRId64
           ",\"row_solve_u_thread_max_rect_entries\":%" PRId64
           ",\"row_solve_partition_ready\":%d"
           ",\"row_solve_partition_slices\":%" PRId64
           ",\"row_solve_l_sparse_level_count\":%" PRId64
           ",\"row_solve_l_sparse_cluster_levels\":%" PRId64
           ",\"row_solve_l_sparse_level_max_width\":%" PRId64
           ",\"row_solve_l_dense_tail_start\":%" PRId64
           ",\"row_solve_l_dense_tail_rows\":%" PRId64
           ",\"row_solve_l_dense_tail_entries\":%" PRId64
           ",\"row_solve_l_slice_max_entries\":%" PRId64
           ",\"row_solve_l_segmented_rows\":%" PRId64
           ",\"row_solve_l_rect_entries\":%" PRId64
           ",\"row_solve_l_tri_entries\":%" PRId64
           ",\"row_solve_u_sparse_level_count\":%" PRId64
           ",\"row_solve_u_sparse_cluster_levels\":%" PRId64
           ",\"row_solve_u_sparse_level_max_width\":%" PRId64
           ",\"row_solve_u_dense_tail_start\":%" PRId64
           ",\"row_solve_u_dense_tail_rows\":%" PRId64
           ",\"row_solve_u_dense_tail_entries\":%" PRId64
           ",\"row_solve_u_slice_max_entries\":%" PRId64
           ",\"row_solve_u_segmented_rows\":%" PRId64
           ",\"row_solve_u_rect_entries\":%" PRId64
           ",\"row_solve_u_tri_entries\":%" PRId64,
           stats.row_solve_parallel_run_count,
           stats.row_solve_parallel_l_slice_runs,
           stats.row_solve_parallel_u_slice_runs,
           stats.row_solve_parallel_l_sparse_level_runs,
           stats.row_solve_parallel_u_sparse_level_runs,
           stats.row_solve_thread_count,
           stats.row_solve_l_thread_max_rect_entries,
           stats.row_solve_u_thread_max_rect_entries,
           stats.row_solve_partition_ready,
           stats.row_solve_partition_slices,
           stats.row_solve_l_sparse_level_count,
           stats.row_solve_l_sparse_cluster_levels,
           stats.row_solve_l_sparse_level_max_width,
           stats.row_solve_l_dense_tail_start,
           stats.row_solve_l_dense_tail_rows,
           stats.row_solve_l_dense_tail_entries,
           stats.row_solve_l_slice_max_entries,
           stats.row_solve_l_segmented_rows,
           stats.row_solve_l_rect_entries,
           stats.row_solve_l_tri_entries,
           stats.row_solve_u_sparse_level_count,
           stats.row_solve_u_sparse_cluster_levels,
           stats.row_solve_u_sparse_level_max_width,
           stats.row_solve_u_dense_tail_start,
           stats.row_solve_u_dense_tail_rows,
           stats.row_solve_u_dense_tail_entries,
           stats.row_solve_u_slice_max_entries,
           stats.row_solve_u_segmented_rows,
           stats.row_solve_u_rect_entries,
           stats.row_solve_u_tri_entries);
    printf(",\"row_refactor_last_work_ready_queue\":%d"
           ",\"row_refactor_work_ready_queue_run_count\":%" PRId64
           ",\"row_refactor_ready_queue_workspace_groups\":%" PRId64
           ",\"row_refactor_last_partial_supernode_pipeline\":%d"
           ",\"row_refactor_last_partial_supernode_pipeline_groups\":%" PRId64
           ",\"row_refactor_last_partial_supernode_pipeline_rows\":%" PRId64
           ",\"row_refactor_partial_supernode_pipeline_run_count\":%" PRId64
           ",\"row_refactor_last_local_ready_groups\":%" PRId64
           ",\"row_refactor_local_ready_group_count\":%" PRId64
           ",\"row_refactor_last_private_ready_groups\":%" PRId64
           ",\"row_refactor_private_ready_group_count\":%" PRId64
           ",\"row_refactor_last_separator_private_queue\":%d"
           ",\"row_refactor_separator_private_queue_run_count\":%" PRId64
           ",\"row_refactor_last_separator_private_components\":%" PRId64
           ",\"row_refactor_separator_private_component_count\":%" PRId64
           ",\"row_refactor_last_separator_flop_queue\":%d"
           ",\"row_refactor_separator_flop_queue_run_count\":%" PRId64
           ",\"row_refactor_last_separator_flop_ordered_private\":%d"
           ",\"row_refactor_separator_flop_ordered_private_run_count\":%" PRId64
           ",\"row_refactor_last_separator_flop_components\":%" PRId64
           ",\"row_refactor_separator_flop_component_count\":%" PRId64
           ",\"row_refactor_last_separator_flop_private_groups\":%" PRId64
           ",\"row_refactor_last_separator_flop_pipeline_groups\":%" PRId64
           ",\"row_refactor_last_separator_flop_closure_groups\":%" PRId64
           ",\"row_refactor_separator_flop_private_group_count\":%" PRId64
           ",\"row_refactor_separator_flop_pipeline_group_count\":%" PRId64
           ",\"row_refactor_separator_flop_closure_group_count\":%" PRId64
           ",\"row_refactor_segment_count\":%" PRId64
           ",\"row_refactor_segment_rows\":%" PRId64
           ",\"row_refactor_segment_max_width\":%" PRId64
           ",\"row_refactor_segment_dense_entries\":%.9g"
           ",\"row_refactor_segment_trailing_entries\":%.9g"
           ",\"row_refactor_dense_segment_count\":%" PRId64
           ",\"row_refactor_dense_segment_rows\":%" PRId64
           ",\"row_refactor_dense_segment_max_width\":%" PRId64
           ",\"row_refactor_dense_segment_dense_entries\":%.9g"
           ",\"row_refactor_dense_segment_trailing_entries\":%.9g"
           ",\"row_refactor_last_compact_dense_panel\":%d"
           ",\"row_refactor_compact_dense_panel_count\":%" PRId64
           ",\"row_refactor_compact_dense_panel_eligible_count\":%" PRId64
           ",\"row_refactor_compact_dense_panel_eligible_rows\":%" PRId64
           ",\"row_refactor_compact_dense_panel_update_work\":%.9g"
           ",\"row_refactor_compact_dense_panel_entries\":%.9g"
           ",\"row_refactor_dense_producer_run_count\":%" PRId64
           ",\"row_refactor_dense_producer_run_rows\":%" PRId64
           ",\"row_refactor_dense_producer_run_dep_rows\":%" PRId64
           ",\"row_refactor_dense_producer_run_max_per_row\":%" PRId64
           ",\"row_refactor_dense_producer_full_suffix_run_count\":%" PRId64
           ",\"row_refactor_dense_producer_full_suffix_rows\":%" PRId64
           ",\"row_refactor_dense_producer_multi_run_rows\":%" PRId64
           ",\"row_refactor_dense_producer_fragmented_rows\":%" PRId64
           ",\"row_refactor_dense_producer_target_count\":%" PRId64
           ",\"row_refactor_dense_producer_target_none_count\":%" PRId64
           ",\"row_refactor_dense_producer_target_external_count\":%" PRId64
           ",\"row_refactor_dense_producer_target_dense_count\":%" PRId64
           ",\"row_refactor_dense_producer_target_pivot_count\":%" PRId64
           ",\"row_refactor_dense_producer_target_trailing_count\":%" PRId64
           ",\"row_refactor_compact_dense_panel_persistent_groups\":%" PRId64
           ",\"row_refactor_compact_dense_panel_persistent_entries\":%" PRId64
           ",\"row_refactor_last_compact_dense_panel_persistent\":%d"
           ",\"row_refactor_compact_dense_panel_persistent_run_count\":%" PRId64
           ",\"row_refactor_last_compact_supernode_update\":%d"
           ",\"row_refactor_compact_supernode_update_count\":%" PRId64
           ",\"row_refactor_compact_supernode_update_rows\":%" PRId64
           ",\"row_refactor_compact_supernode_update_entries\":%" PRId64
           ",\"row_refactor_last_compact_supernode_partial_update\":%d"
           ",\"row_refactor_compact_supernode_partial_update_count\":%" PRId64
           ",\"row_refactor_compact_supernode_partial_update_rows\":%" PRId64
           ",\"row_refactor_compact_supernode_partial_update_entries\":%" PRId64
           ",\"row_refactor_last_compact_supernode_gemv\":%d"
           ",\"row_refactor_compact_supernode_gemv_count\":%" PRId64
           ",\"row_refactor_compact_supernode_gemv_rows\":%" PRId64
           ",\"row_refactor_compact_supernode_gemv_entries\":%" PRId64
           ",\"row_refactor_last_compact_supernode_trsv\":%d"
           ",\"row_refactor_compact_supernode_trsv_count\":%" PRId64
           ",\"row_refactor_compact_supernode_trsv_rows\":%" PRId64
           ",\"row_refactor_compact_supernode_trsv_entries\":%" PRId64
           ",\"row_refactor_last_compact_supernode_batch\":%d"
           ",\"row_refactor_compact_supernode_batch_count\":%" PRId64
           ",\"row_refactor_compact_supernode_batch_rows\":%" PRId64
           ",\"row_refactor_compact_supernode_batch_dep_rows\":%" PRId64
           ",\"row_refactor_compact_supernode_batch_entries\":%" PRId64,
           stats.row_refactor_last_work_ready_queue,
           stats.row_refactor_work_ready_queue_run_count,
           stats.row_refactor_ready_queue_workspace_groups,
           stats.row_refactor_last_partial_supernode_pipeline,
           stats.row_refactor_last_partial_supernode_pipeline_groups,
           stats.row_refactor_last_partial_supernode_pipeline_rows,
           stats.row_refactor_partial_supernode_pipeline_run_count,
           stats.row_refactor_last_local_ready_groups,
           stats.row_refactor_local_ready_group_count,
           stats.row_refactor_last_private_ready_groups,
           stats.row_refactor_private_ready_group_count,
           stats.row_refactor_last_separator_private_queue,
           stats.row_refactor_separator_private_queue_run_count,
           stats.row_refactor_last_separator_private_components,
           stats.row_refactor_separator_private_component_count,
           stats.row_refactor_last_separator_flop_queue,
           stats.row_refactor_separator_flop_queue_run_count,
           stats.row_refactor_last_separator_flop_ordered_private,
           stats.row_refactor_separator_flop_ordered_private_run_count,
           stats.row_refactor_last_separator_flop_components,
           stats.row_refactor_separator_flop_component_count,
           stats.row_refactor_last_separator_flop_private_groups,
           stats.row_refactor_last_separator_flop_pipeline_groups,
           stats.row_refactor_last_separator_flop_closure_groups,
           stats.row_refactor_separator_flop_private_group_count,
           stats.row_refactor_separator_flop_pipeline_group_count,
           stats.row_refactor_separator_flop_closure_group_count,
           stats.row_refactor_segment_count,
           stats.row_refactor_segment_rows,
           stats.row_refactor_segment_max_width,
           stats.row_refactor_segment_dense_entries,
           stats.row_refactor_segment_trailing_entries,
           stats.row_refactor_dense_segment_count,
           stats.row_refactor_dense_segment_rows,
           stats.row_refactor_dense_segment_max_width,
           stats.row_refactor_dense_segment_dense_entries,
           stats.row_refactor_dense_segment_trailing_entries,
           stats.row_refactor_last_compact_dense_panel,
           stats.row_refactor_compact_dense_panel_count,
           stats.row_refactor_compact_dense_panel_eligible_count,
           stats.row_refactor_compact_dense_panel_eligible_rows,
           stats.row_refactor_compact_dense_panel_update_work,
           stats.row_refactor_compact_dense_panel_entries,
           stats.row_refactor_dense_producer_run_count,
           stats.row_refactor_dense_producer_run_rows,
           stats.row_refactor_dense_producer_run_dep_rows,
           stats.row_refactor_dense_producer_run_max_per_row,
           stats.row_refactor_dense_producer_full_suffix_run_count,
           stats.row_refactor_dense_producer_full_suffix_rows,
           stats.row_refactor_dense_producer_multi_run_rows,
           stats.row_refactor_dense_producer_fragmented_rows,
           stats.row_refactor_dense_producer_target_count,
           stats.row_refactor_dense_producer_target_none_count,
           stats.row_refactor_dense_producer_target_external_count,
           stats.row_refactor_dense_producer_target_dense_count,
           stats.row_refactor_dense_producer_target_pivot_count,
           stats.row_refactor_dense_producer_target_trailing_count,
           stats.row_refactor_compact_dense_panel_persistent_groups,
           stats.row_refactor_compact_dense_panel_persistent_entries,
           stats.row_refactor_last_compact_dense_panel_persistent,
           stats.row_refactor_compact_dense_panel_persistent_run_count,
           stats.row_refactor_last_compact_supernode_update,
           stats.row_refactor_compact_supernode_update_count,
           stats.row_refactor_compact_supernode_update_rows,
           stats.row_refactor_compact_supernode_update_entries,
           stats.row_refactor_last_compact_supernode_partial_update,
           stats.row_refactor_compact_supernode_partial_update_count,
           stats.row_refactor_compact_supernode_partial_update_rows,
           stats.row_refactor_compact_supernode_partial_update_entries,
           stats.row_refactor_last_compact_supernode_gemv,
           stats.row_refactor_compact_supernode_gemv_count,
           stats.row_refactor_compact_supernode_gemv_rows,
           stats.row_refactor_compact_supernode_gemv_entries,
           stats.row_refactor_last_compact_supernode_trsv,
           stats.row_refactor_compact_supernode_trsv_count,
           stats.row_refactor_compact_supernode_trsv_rows,
           stats.row_refactor_compact_supernode_trsv_entries,
           stats.row_refactor_last_compact_supernode_batch,
           stats.row_refactor_compact_supernode_batch_count,
           stats.row_refactor_compact_supernode_batch_rows,
           stats.row_refactor_compact_supernode_batch_dep_rows,
           stats.row_refactor_compact_supernode_batch_entries);
    printf(",\"row_refactor_last_separator_flop_private_threads\":%" PRId64
           ",\"row_refactor_last_separator_flop_private_min_groups\":%" PRId64
           ",\"row_refactor_last_separator_flop_private_max_groups\":%" PRId64
           ",\"row_refactor_last_separator_flop_private_min_work\":%.9g"
           ",\"row_refactor_last_separator_flop_private_max_work\":%.9g",
           stats.row_refactor_last_separator_flop_private_threads,
           stats.row_refactor_last_separator_flop_private_min_groups,
           stats.row_refactor_last_separator_flop_private_max_groups,
           stats.row_refactor_last_separator_flop_private_min_work,
           stats.row_refactor_last_separator_flop_private_max_work);
    printf(",\"row_refactor_last_compact_dense_panel_blocked\":%d"
           ",\"row_refactor_compact_dense_panel_blocked_run_count\":%" PRId64
           ",\"row_refactor_compact_dense_panel_blocked_rows\":%" PRId64
           ",\"row_refactor_compact_dense_panel_blocked_entries\":%" PRId64,
           stats.row_refactor_last_compact_dense_panel_blocked,
           stats.row_refactor_compact_dense_panel_blocked_run_count,
           stats.row_refactor_compact_dense_panel_blocked_rows,
           stats.row_refactor_compact_dense_panel_blocked_entries);
    printf(",\"row_refactor_native_row_panel_enabled\":%d"
           ",\"row_refactor_last_native_row_panel\":%d"
           ",\"row_refactor_native_row_panel_auto_disabled\":%d"
           ",\"row_refactor_native_row_panel_count\":%" PRId64
           ",\"row_refactor_native_row_panel_rows\":%" PRId64
           ",\"row_refactor_native_row_panel_entries\":%" PRId64
           ",\"row_refactor_native_row_panel_blocked_count\":%" PRId64
           ",\"row_refactor_native_row_panel_blocked_rows\":%" PRId64
           ",\"row_refactor_native_row_panel_blocked_entries\":%" PRId64
           ",\"row_refactor_native_row_panel_fallback_count\":%" PRId64
           ",\"row_refactor_native_row_panel_checked_reject_count\":%" PRId64
           ",\"row_refactor_native_row_panel_auto_disable_count\":%" PRId64,
           stats.row_refactor_native_row_panel_enabled,
           stats.row_refactor_last_native_row_panel,
           stats.row_refactor_native_row_panel_auto_disabled,
           stats.row_refactor_native_row_panel_count,
           stats.row_refactor_native_row_panel_rows,
           stats.row_refactor_native_row_panel_entries,
           stats.row_refactor_native_row_panel_blocked_count,
           stats.row_refactor_native_row_panel_blocked_rows,
           stats.row_refactor_native_row_panel_blocked_entries,
           stats.row_refactor_native_row_panel_fallback_count,
           stats.row_refactor_native_row_panel_checked_reject_count,
           stats.row_refactor_native_row_panel_auto_disable_count);
    printf(",\"row_refactor_last_compact_dense_panel_direct_input_rows\":%" PRId64
           ",\"row_refactor_compact_dense_panel_direct_input_rows\":%" PRId64,
           stats.row_refactor_last_compact_dense_panel_direct_input_rows,
           stats.row_refactor_compact_dense_panel_direct_input_rows);
    printf(",\"row_refactor_last_compact_panel_solve_values\":%" PRId64
           ",\"row_refactor_compact_panel_solve_values\":%" PRId64,
           stats.row_refactor_last_compact_panel_solve_values,
           stats.row_refactor_compact_panel_solve_values);
    printf(",\"row_refactor_last_compact_panel_group_solve_rows\":%" PRId64
           ",\"row_refactor_compact_panel_group_solve_rows\":%" PRId64
           ",\"row_refactor_last_compact_panel_group_solve_entries\":%" PRId64
           ",\"row_refactor_compact_panel_group_solve_entries\":%" PRId64,
           stats.row_refactor_last_compact_panel_group_solve_rows,
           stats.row_refactor_compact_panel_group_solve_rows,
           stats.row_refactor_last_compact_panel_group_solve_entries,
           stats.row_refactor_compact_panel_group_solve_entries);
    printf(",\"row_refactor_last_compact_panel_scalar_update_rows\":%" PRId64
           ",\"row_refactor_compact_panel_scalar_update_rows\":%" PRId64
           ",\"row_refactor_last_compact_panel_scalar_update_entries\":%" PRId64
           ",\"row_refactor_compact_panel_scalar_update_entries\":%" PRId64,
           stats.row_refactor_last_compact_panel_scalar_update_rows,
           stats.row_refactor_compact_panel_scalar_update_rows,
           stats.row_refactor_last_compact_panel_scalar_update_entries,
           stats.row_refactor_compact_panel_scalar_update_entries);
    printf(",\"row_refactor_last_dense_segment_direct_input_rows\":%" PRId64
           ",\"row_refactor_dense_segment_direct_input_rows\":%" PRId64,
           stats.row_refactor_last_dense_segment_direct_input_rows,
           stats.row_refactor_dense_segment_direct_input_rows);
    printf(",\"row_refactor_last_sparse_segment_direct_input_rows\":%" PRId64
           ",\"row_refactor_sparse_segment_direct_input_rows\":%" PRId64,
           stats.row_refactor_last_sparse_segment_direct_input_rows,
           stats.row_refactor_sparse_segment_direct_input_rows);
    printf(",\"row_refactor_last_batch_direct_input_rows\":%" PRId64
           ",\"row_refactor_batch_direct_input_rows\":%" PRId64,
           stats.row_refactor_last_batch_direct_input_rows,
           stats.row_refactor_batch_direct_input_rows);
    printf(",\"row_refactor_segment_input_target_rows\":%" PRId64
           ",\"row_refactor_segment_input_target_entries\":%" PRId64
           ",\"row_refactor_last_segment_target_input_rows\":%" PRId64
           ",\"row_refactor_segment_target_input_rows\":%" PRId64
           ",\"row_refactor_segment_input_cleanup_rows\":%" PRId64
           ",\"row_refactor_segment_input_cleanup_entries\":%" PRId64
           ",\"row_refactor_last_segment_target_cleanup_rows\":%" PRId64
           ",\"row_refactor_segment_target_cleanup_rows\":%" PRId64
           ",\"row_refactor_last_segment_target_cleanup_entries\":%" PRId64
           ",\"row_refactor_segment_target_cleanup_entries\":%" PRId64,
           stats.row_refactor_segment_input_target_rows,
           stats.row_refactor_segment_input_target_entries,
           stats.row_refactor_last_segment_target_input_rows,
           stats.row_refactor_segment_target_input_rows,
           stats.row_refactor_segment_input_cleanup_rows,
           stats.row_refactor_segment_input_cleanup_entries,
           stats.row_refactor_last_segment_target_cleanup_rows,
           stats.row_refactor_segment_target_cleanup_rows,
           stats.row_refactor_last_segment_target_cleanup_entries,
           stats.row_refactor_segment_target_cleanup_entries);
    printf(",\"row_refactor_compact_supernode_batch_pattern_count\":%" PRId64
           ",\"row_refactor_compact_supernode_batch_pattern_rows\":%" PRId64
           ",\"row_refactor_compact_supernode_batch_candidate_count\":%" PRId64
           ",\"row_refactor_compact_supernode_batch_candidate_rows\":%" PRId64
           ",\"row_refactor_compact_supernode_batch_candidate_dep_rows\":%" PRId64
           ",\"row_refactor_compact_supernode_batch_rejected_work_count\":%" PRId64
           ",\"refactor_dependency_cluster_levels\":%" PRId64
           ",\"refactor_dependency_pipeline_columns\":%" PRId64
           ",\"refactor_dependency_work\":%.9g"
           ",\"refactor_dependency_pipeline_work\":%.9g"
           ",\"refactor_stream_dependency_entries\":%.9g"
           ",\"refactor_stream_pivot_entries\":%.9g"
           ",\"refactor_stream_output_entries\":%.9g"
           ",\"refactor_last_supernode_pipeline_tasks\":%" PRId64
           ",\"refactor_last_supernode_pipeline_columns\":%" PRId64
           ",\"refactor_supernode_pipeline_task_count\":%" PRId64
           ",\"refactor_supernode_pipeline_column_count\":%" PRId64
           ",\"refactor_last_supernode_update_runs\":%" PRId64
           ",\"refactor_last_supernode_update_rows\":%" PRId64
           ",\"refactor_last_supernode_update_entries\":%" PRId64
           ",\"refactor_supernode_update_run_count\":%" PRId64
           ",\"refactor_supernode_update_rows\":%" PRId64
           ",\"refactor_supernode_update_entries\":%" PRId64
           ",\"refactor_last_supernode_cblas_update_runs\":%" PRId64
           ",\"refactor_last_supernode_cblas_update_rows\":%" PRId64
           ",\"refactor_last_supernode_cblas_update_entries\":%" PRId64
           ",\"refactor_supernode_cblas_update_run_count\":%" PRId64
           ",\"refactor_supernode_cblas_update_rows\":%" PRId64
           ",\"refactor_supernode_cblas_update_entries\":%" PRId64
           ",\"refactor_last_supernode_blocked_update_runs\":%" PRId64
           ",\"refactor_last_supernode_blocked_update_rows\":%" PRId64
           ",\"refactor_last_supernode_blocked_update_entries\":%" PRId64
           ",\"refactor_supernode_blocked_update_run_count\":%" PRId64
           ",\"refactor_supernode_blocked_update_rows\":%" PRId64
           ",\"refactor_supernode_blocked_update_entries\":%" PRId64
           ",\"refactor_last_supernode_cached_probe_attempts\":%" PRId64
           ",\"refactor_last_supernode_cached_probe_panel_hits\":%" PRId64
           ",\"refactor_last_supernode_cached_probe_contiguous\":%" PRId64
           ",\"refactor_last_supernode_cached_probe_allowed\":%" PRId64
           ",\"refactor_last_supernode_cached_probe_allowed_rows\":%" PRId64
           ",\"refactor_last_supernode_cached_probe_applied\":%" PRId64
           ",\"refactor_last_supernode_cached_probe_applied_rows\":%" PRId64
           ",\"refactor_supernode_cached_probe_attempt_count\":%" PRId64
           ",\"refactor_supernode_cached_probe_panel_hits\":%" PRId64
           ",\"refactor_supernode_cached_probe_contiguous\":%" PRId64
           ",\"refactor_supernode_cached_probe_allowed\":%" PRId64
           ",\"refactor_supernode_cached_probe_allowed_rows\":%" PRId64
           ",\"refactor_supernode_cached_probe_applied\":%" PRId64
           ",\"refactor_supernode_cached_probe_applied_rows\":%" PRId64
           ",\"refactor_supernode_cached_probe_disabled\":%d"
           ",\"refactor_supernode_cached_probe_disable_count\":%" PRId64
           ",\"refactor_last_btf_scalar_run_candidates\":%" PRId64
           ",\"refactor_last_btf_scalar_run_rows\":%" PRId64
           ",\"refactor_last_btf_scalar_run_entries\":%" PRId64
           ",\"refactor_last_btf_scalar_run_max_rows\":%" PRId64
           ",\"refactor_btf_scalar_run_candidate_count\":%" PRId64
           ",\"refactor_btf_scalar_run_rows\":%" PRId64
           ",\"refactor_btf_scalar_run_entries\":%" PRId64
           ",\"refactor_btf_scalar_run_max_rows\":%" PRId64
           ",\"refactor_last_btf_scalar_run_exec_runs\":%" PRId64
           ",\"refactor_last_btf_scalar_run_exec_rows\":%" PRId64
           ",\"refactor_last_btf_scalar_run_exec_entries\":%" PRId64
           ",\"refactor_last_btf_scalar_run_exec_max_rows\":%" PRId64
           ",\"refactor_btf_scalar_run_exec_count\":%" PRId64
           ",\"refactor_btf_scalar_run_exec_rows\":%" PRId64
           ",\"refactor_btf_scalar_run_exec_entries\":%" PRId64
           ",\"refactor_btf_scalar_run_exec_max_rows\":%" PRId64
           ",\"refactor_supernode_update_disabled\":%d"
           ",\"refactor_supernode_update_disable_count\":%" PRId64
           ",\"refactor_last_ready_queue_columns\":%" PRId64
           ",\"refactor_ready_queue_run_count\":%" PRId64
           ",\"nnz_l\":%" PRId64 ",\"nnz_u\":%" PRId64
           ",\"estimated_flops\":%.9g,\"factor_flops\":%.9g"
           ",\"parallel_model_r1\":%.9g"
           ",\"parallel_model_r2\":%.9g"
           ",\"parallel_model_recommends_parallel\":%d"
           ",\"parallel_task_flow_threads\":%" PRId64
           ",\"parallel_task_flow_dependencies\":%" PRId64
           ",\"parallel_task_flow_work\":%.9g"
           ",\"parallel_task_flow_finish_time\":%.9g"
           ",\"parallel_task_flow_speedup\":%.9g"
           ",\"parallel_task_flow_recommends_parallel\":%d"
           ",\"rcond\":%.9g,\"rgrowth\":%.9g"
           ",\"memory_bytes\":%zu,\"memory_peak_bytes\":%zu",
           stats.row_refactor_compact_supernode_batch_pattern_count,
           stats.row_refactor_compact_supernode_batch_pattern_rows,
           stats.row_refactor_compact_supernode_batch_candidate_count,
           stats.row_refactor_compact_supernode_batch_candidate_rows,
           stats.row_refactor_compact_supernode_batch_candidate_dep_rows,
           stats.row_refactor_compact_supernode_batch_rejected_work_count,
           stats.refactor_dependency_cluster_levels,
           stats.refactor_dependency_pipeline_columns,
           stats.refactor_dependency_work,
           stats.refactor_dependency_pipeline_work,
           stats.refactor_stream_dependency_entries,
           stats.refactor_stream_pivot_entries,
           stats.refactor_stream_output_entries,
           stats.refactor_last_supernode_pipeline_tasks,
           stats.refactor_last_supernode_pipeline_columns,
           stats.refactor_supernode_pipeline_task_count,
           stats.refactor_supernode_pipeline_column_count,
           stats.refactor_last_supernode_update_runs,
           stats.refactor_last_supernode_update_rows,
           stats.refactor_last_supernode_update_entries,
           stats.refactor_supernode_update_run_count,
           stats.refactor_supernode_update_rows,
           stats.refactor_supernode_update_entries,
           stats.refactor_last_supernode_cblas_update_runs,
           stats.refactor_last_supernode_cblas_update_rows,
           stats.refactor_last_supernode_cblas_update_entries,
           stats.refactor_supernode_cblas_update_run_count,
           stats.refactor_supernode_cblas_update_rows,
           stats.refactor_supernode_cblas_update_entries,
           stats.refactor_last_supernode_blocked_update_runs,
           stats.refactor_last_supernode_blocked_update_rows,
           stats.refactor_last_supernode_blocked_update_entries,
           stats.refactor_supernode_blocked_update_run_count,
           stats.refactor_supernode_blocked_update_rows,
           stats.refactor_supernode_blocked_update_entries,
           stats.refactor_last_supernode_cached_probe_attempts,
           stats.refactor_last_supernode_cached_probe_panel_hits,
           stats.refactor_last_supernode_cached_probe_contiguous,
           stats.refactor_last_supernode_cached_probe_allowed,
           stats.refactor_last_supernode_cached_probe_allowed_rows,
           stats.refactor_last_supernode_cached_probe_applied,
           stats.refactor_last_supernode_cached_probe_applied_rows,
           stats.refactor_supernode_cached_probe_attempt_count,
           stats.refactor_supernode_cached_probe_panel_hits,
           stats.refactor_supernode_cached_probe_contiguous,
           stats.refactor_supernode_cached_probe_allowed,
           stats.refactor_supernode_cached_probe_allowed_rows,
           stats.refactor_supernode_cached_probe_applied,
           stats.refactor_supernode_cached_probe_applied_rows,
           stats.refactor_supernode_cached_probe_disabled,
           stats.refactor_supernode_cached_probe_disable_count,
           stats.refactor_last_btf_scalar_run_candidates,
           stats.refactor_last_btf_scalar_run_rows,
           stats.refactor_last_btf_scalar_run_entries,
           stats.refactor_last_btf_scalar_run_max_rows,
           stats.refactor_btf_scalar_run_candidate_count,
           stats.refactor_btf_scalar_run_rows,
           stats.refactor_btf_scalar_run_entries,
           stats.refactor_btf_scalar_run_max_rows,
           stats.refactor_last_btf_scalar_run_exec_runs,
           stats.refactor_last_btf_scalar_run_exec_rows,
           stats.refactor_last_btf_scalar_run_exec_entries,
           stats.refactor_last_btf_scalar_run_exec_max_rows,
           stats.refactor_btf_scalar_run_exec_count,
           stats.refactor_btf_scalar_run_exec_rows,
           stats.refactor_btf_scalar_run_exec_entries,
           stats.refactor_btf_scalar_run_exec_max_rows,
           stats.refactor_supernode_update_disabled,
           stats.refactor_supernode_update_disable_count,
           stats.refactor_last_ready_queue_columns,
           stats.refactor_ready_queue_run_count,
           stats.nnz_l, stats.nnz_u,
           stats.estimated_flops, stats.factor_flops,
           stats.parallel_model_r1,
           stats.parallel_model_r2,
           stats.parallel_model_recommends_parallel,
           stats.parallel_task_flow_threads,
           stats.parallel_task_flow_dependencies,
           stats.parallel_task_flow_work,
           stats.parallel_task_flow_finish_time,
           stats.parallel_task_flow_speedup,
           stats.parallel_task_flow_recommends_parallel,
           stats.rcond, stats.rgrowth,
           stats.memory_bytes, stats.memory_peak_bytes);
    printf(",\"refactor_btf_scalar_run_group_built\":%d"
           ",\"refactor_btf_scalar_run_group_count\":%" PRId64
           ",\"refactor_btf_scalar_run_group_current_total\":%" PRId64
           ",\"refactor_btf_scalar_run_group_multi_count\":%" PRId64
           ",\"refactor_btf_scalar_run_group_multi_current_total\":%" PRId64
           ",\"refactor_btf_scalar_run_group_rows\":%" PRId64
           ",\"refactor_btf_scalar_run_group_reused_rows\":%" PRId64
           ",\"refactor_btf_scalar_run_group_entries\":%" PRId64
           ",\"refactor_btf_scalar_run_group_reused_entries\":%" PRId64
           ",\"refactor_btf_scalar_run_group_max_currents\":%" PRId64
           ",\"refactor_btf_scalar_run_group_max_rows\":%" PRId64
           ",\"refactor_btf_scalar_run_group_state_rows_total\":%" PRId64
           ",\"refactor_btf_scalar_run_group_state_max_rows\":%" PRId64,
           stats.refactor_btf_scalar_run_group_built,
           stats.refactor_btf_scalar_run_group_count,
           stats.refactor_btf_scalar_run_group_current_total,
           stats.refactor_btf_scalar_run_group_multi_count,
           stats.refactor_btf_scalar_run_group_multi_current_total,
           stats.refactor_btf_scalar_run_group_rows,
           stats.refactor_btf_scalar_run_group_reused_rows,
           stats.refactor_btf_scalar_run_group_entries,
           stats.refactor_btf_scalar_run_group_reused_entries,
           stats.refactor_btf_scalar_run_group_max_currents,
           stats.refactor_btf_scalar_run_group_max_rows,
           stats.refactor_btf_scalar_run_group_state_rows_total,
           stats.refactor_btf_scalar_run_group_state_max_rows);
    printf(",\"refactor_last_btf_scalar_run_group_waits\":%" PRId64
           ",\"refactor_last_btf_scalar_run_group_wait_rows\":%" PRId64
           ",\"refactor_last_btf_scalar_run_group_wait_entries\":%" PRId64
           ",\"refactor_last_btf_scalar_run_group_overlap_waits\":%" PRId64
           ",\"refactor_last_btf_scalar_run_group_overlap_rows\":%" PRId64
           ",\"refactor_last_btf_scalar_run_group_overlap_entries\":%" PRId64
           ",\"refactor_last_btf_scalar_run_group_max_live\":%" PRId64
           ",\"refactor_btf_scalar_run_group_wait_count\":%" PRId64
           ",\"refactor_btf_scalar_run_group_wait_rows\":%" PRId64
           ",\"refactor_btf_scalar_run_group_wait_entries\":%" PRId64
           ",\"refactor_btf_scalar_run_group_overlap_count\":%" PRId64
           ",\"refactor_btf_scalar_run_group_overlap_rows\":%" PRId64
           ",\"refactor_btf_scalar_run_group_overlap_entries\":%" PRId64
           ",\"refactor_btf_scalar_run_group_max_live\":%" PRId64
           ",\"refactor_last_btf_scalar_run_group_claim_surface_triggers\":%" PRId64
           ",\"refactor_last_btf_scalar_run_group_claim_surface_groups\":%" PRId64
           ",\"refactor_last_btf_scalar_run_group_claim_surface_currents\":%" PRId64
           ",\"refactor_last_btf_scalar_run_group_claim_triggers\":%" PRId64
           ",\"refactor_last_btf_scalar_run_group_claim_groups\":%" PRId64
           ",\"refactor_last_btf_scalar_run_group_claim_currents\":%" PRId64
           ",\"refactor_last_btf_scalar_run_group_prefix_ready_currents\":%" PRId64
           ",\"refactor_last_btf_scalar_run_group_prefix_ready_rows\":%" PRId64
           ",\"refactor_last_btf_scalar_run_group_prefix_ready_entries\":%" PRId64
           ",\"refactor_last_btf_scalar_run_group_run_ready_currents\":%" PRId64
           ",\"refactor_last_btf_scalar_run_group_run_ready_rows\":%" PRId64
           ",\"refactor_last_btf_scalar_run_group_run_ready_entries\":%" PRId64
           ",\"refactor_last_btf_scalar_run_group_wake_armed_currents\":%" PRId64
           ",\"refactor_last_btf_scalar_run_group_wake_armed_rows\":%" PRId64
           ",\"refactor_last_btf_scalar_run_group_wake_armed_entries\":%" PRId64
           ",\"refactor_last_btf_scalar_run_group_wake_ready_currents\":%" PRId64
           ",\"refactor_last_btf_scalar_run_group_wake_ready_rows\":%" PRId64
           ",\"refactor_last_btf_scalar_run_group_wake_ready_entries\":%" PRId64
           ",\"refactor_last_btf_scalar_run_group_state_materialized_currents\":%" PRId64
           ",\"refactor_last_btf_scalar_run_group_state_materialized_rows\":%" PRId64
           ",\"refactor_last_btf_scalar_run_group_state_materialized_prefix_deps\":%" PRId64
           ",\"refactor_last_btf_scalar_run_group_state_advanced_currents\":%" PRId64
           ",\"refactor_last_btf_scalar_run_group_state_advanced_rows\":%" PRId64
           ",\"refactor_last_btf_scalar_run_group_state_advanced_entries\":%" PRId64
           ",\"refactor_last_btf_scalar_run_group_state_rejects\":%" PRId64
           ",\"refactor_btf_scalar_run_group_claim_surface_trigger_count\":%" PRId64
           ",\"refactor_btf_scalar_run_group_claim_surface_group_count\":%" PRId64
           ",\"refactor_btf_scalar_run_group_claim_surface_current_count\":%" PRId64
           ",\"refactor_btf_scalar_run_group_claim_trigger_count\":%" PRId64
           ",\"refactor_btf_scalar_run_group_claim_group_count\":%" PRId64
           ",\"refactor_btf_scalar_run_group_claim_current_count\":%" PRId64
           ",\"refactor_btf_scalar_run_group_prefix_ready_current_count\":%" PRId64
           ",\"refactor_btf_scalar_run_group_prefix_ready_rows\":%" PRId64
           ",\"refactor_btf_scalar_run_group_prefix_ready_entries\":%" PRId64
           ",\"refactor_btf_scalar_run_group_run_ready_current_count\":%" PRId64
           ",\"refactor_btf_scalar_run_group_run_ready_rows\":%" PRId64
           ",\"refactor_btf_scalar_run_group_run_ready_entries\":%" PRId64
           ",\"refactor_btf_scalar_run_group_wake_armed_current_count\":%" PRId64
           ",\"refactor_btf_scalar_run_group_wake_armed_rows\":%" PRId64
           ",\"refactor_btf_scalar_run_group_wake_armed_entries\":%" PRId64
           ",\"refactor_btf_scalar_run_group_wake_ready_current_count\":%" PRId64
           ",\"refactor_btf_scalar_run_group_wake_ready_rows\":%" PRId64
           ",\"refactor_btf_scalar_run_group_wake_ready_entries\":%" PRId64
           ",\"refactor_btf_scalar_run_group_state_materialized_current_count\":%" PRId64
           ",\"refactor_btf_scalar_run_group_state_materialized_rows\":%" PRId64
           ",\"refactor_btf_scalar_run_group_state_materialized_prefix_deps\":%" PRId64
           ",\"refactor_btf_scalar_run_group_state_advanced_current_count\":%" PRId64
           ",\"refactor_btf_scalar_run_group_state_advanced_rows\":%" PRId64
           ",\"refactor_btf_scalar_run_group_state_advanced_entries\":%" PRId64
           ",\"refactor_btf_scalar_run_group_state_reject_count\":%" PRId64,
           stats.refactor_last_btf_scalar_run_group_waits,
           stats.refactor_last_btf_scalar_run_group_wait_rows,
           stats.refactor_last_btf_scalar_run_group_wait_entries,
           stats.refactor_last_btf_scalar_run_group_overlap_waits,
           stats.refactor_last_btf_scalar_run_group_overlap_rows,
           stats.refactor_last_btf_scalar_run_group_overlap_entries,
           stats.refactor_last_btf_scalar_run_group_max_live,
           stats.refactor_btf_scalar_run_group_wait_count,
           stats.refactor_btf_scalar_run_group_wait_rows,
           stats.refactor_btf_scalar_run_group_wait_entries,
           stats.refactor_btf_scalar_run_group_overlap_count,
           stats.refactor_btf_scalar_run_group_overlap_rows,
           stats.refactor_btf_scalar_run_group_overlap_entries,
           stats.refactor_btf_scalar_run_group_max_live,
           stats.refactor_last_btf_scalar_run_group_claim_surface_triggers,
           stats.refactor_last_btf_scalar_run_group_claim_surface_groups,
           stats.refactor_last_btf_scalar_run_group_claim_surface_currents,
           stats.refactor_last_btf_scalar_run_group_claim_triggers,
           stats.refactor_last_btf_scalar_run_group_claim_groups,
           stats.refactor_last_btf_scalar_run_group_claim_currents,
           stats.refactor_last_btf_scalar_run_group_prefix_ready_currents,
           stats.refactor_last_btf_scalar_run_group_prefix_ready_rows,
           stats.refactor_last_btf_scalar_run_group_prefix_ready_entries,
           stats.refactor_last_btf_scalar_run_group_run_ready_currents,
           stats.refactor_last_btf_scalar_run_group_run_ready_rows,
           stats.refactor_last_btf_scalar_run_group_run_ready_entries,
           stats.refactor_last_btf_scalar_run_group_wake_armed_currents,
           stats.refactor_last_btf_scalar_run_group_wake_armed_rows,
           stats.refactor_last_btf_scalar_run_group_wake_armed_entries,
           stats.refactor_last_btf_scalar_run_group_wake_ready_currents,
           stats.refactor_last_btf_scalar_run_group_wake_ready_rows,
           stats.refactor_last_btf_scalar_run_group_wake_ready_entries,
           stats.refactor_last_btf_scalar_run_group_state_materialized_currents,
           stats.refactor_last_btf_scalar_run_group_state_materialized_rows,
           stats.refactor_last_btf_scalar_run_group_state_materialized_prefix_deps,
           stats.refactor_last_btf_scalar_run_group_state_advanced_currents,
           stats.refactor_last_btf_scalar_run_group_state_advanced_rows,
           stats.refactor_last_btf_scalar_run_group_state_advanced_entries,
           stats.refactor_last_btf_scalar_run_group_state_rejects,
           stats.refactor_btf_scalar_run_group_claim_surface_trigger_count,
           stats.refactor_btf_scalar_run_group_claim_surface_group_count,
           stats.refactor_btf_scalar_run_group_claim_surface_current_count,
           stats.refactor_btf_scalar_run_group_claim_trigger_count,
           stats.refactor_btf_scalar_run_group_claim_group_count,
           stats.refactor_btf_scalar_run_group_claim_current_count,
           stats.refactor_btf_scalar_run_group_prefix_ready_current_count,
           stats.refactor_btf_scalar_run_group_prefix_ready_rows,
           stats.refactor_btf_scalar_run_group_prefix_ready_entries,
           stats.refactor_btf_scalar_run_group_run_ready_current_count,
           stats.refactor_btf_scalar_run_group_run_ready_rows,
           stats.refactor_btf_scalar_run_group_run_ready_entries,
           stats.refactor_btf_scalar_run_group_wake_armed_current_count,
           stats.refactor_btf_scalar_run_group_wake_armed_rows,
           stats.refactor_btf_scalar_run_group_wake_armed_entries,
           stats.refactor_btf_scalar_run_group_wake_ready_current_count,
           stats.refactor_btf_scalar_run_group_wake_ready_rows,
           stats.refactor_btf_scalar_run_group_wake_ready_entries,
           stats.refactor_btf_scalar_run_group_state_materialized_current_count,
           stats.refactor_btf_scalar_run_group_state_materialized_rows,
           stats.refactor_btf_scalar_run_group_state_materialized_prefix_deps,
           stats.refactor_btf_scalar_run_group_state_advanced_current_count,
           stats.refactor_btf_scalar_run_group_state_advanced_rows,
           stats.refactor_btf_scalar_run_group_state_advanced_entries,
           stats.refactor_btf_scalar_run_group_state_reject_count);
    printf(",\"refactor_last_supernode_cached_probe_shape_rejects\":%" PRId64
           ",\"refactor_last_supernode_cached_probe_shape_reject_rows\":%" PRId64
           ",\"refactor_last_supernode_cached_probe_stream_rejects\":%" PRId64
           ",\"refactor_last_supernode_cached_probe_stream_reject_rows\":%" PRId64
           ",\"refactor_last_supernode_cached_probe_work_rejects\":%" PRId64
           ",\"refactor_last_supernode_cached_probe_work_reject_rows\":%" PRId64
           ",\"refactor_last_supernode_cached_probe_workspace_rejects\":%" PRId64
           ",\"refactor_last_supernode_cached_probe_workspace_reject_rows\":%" PRId64
           ",\"refactor_supernode_cached_probe_shape_rejects\":%" PRId64
           ",\"refactor_supernode_cached_probe_shape_reject_rows\":%" PRId64
           ",\"refactor_supernode_cached_probe_stream_rejects\":%" PRId64
           ",\"refactor_supernode_cached_probe_stream_reject_rows\":%" PRId64
           ",\"refactor_supernode_cached_probe_work_rejects\":%" PRId64
           ",\"refactor_supernode_cached_probe_work_reject_rows\":%" PRId64
           ",\"refactor_supernode_cached_probe_workspace_rejects\":%" PRId64
           ",\"refactor_supernode_cached_probe_workspace_reject_rows\":%" PRId64
           "}\n",
           stats.refactor_last_supernode_cached_probe_shape_rejects,
           stats.refactor_last_supernode_cached_probe_shape_reject_rows,
           stats.refactor_last_supernode_cached_probe_stream_rejects,
           stats.refactor_last_supernode_cached_probe_stream_reject_rows,
           stats.refactor_last_supernode_cached_probe_work_rejects,
           stats.refactor_last_supernode_cached_probe_work_reject_rows,
           stats.refactor_last_supernode_cached_probe_workspace_rejects,
           stats.refactor_last_supernode_cached_probe_workspace_reject_rows,
           stats.refactor_supernode_cached_probe_shape_rejects,
           stats.refactor_supernode_cached_probe_shape_reject_rows,
           stats.refactor_supernode_cached_probe_stream_rejects,
           stats.refactor_supernode_cached_probe_stream_reject_rows,
           stats.refactor_supernode_cached_probe_work_rejects,
           stats.refactor_supernode_cached_probe_work_reject_rows,
           stats.refactor_supernode_cached_probe_workspace_rejects,
           stats.refactor_supernode_cached_probe_workspace_reject_rows);
  } else {
    printf("matrix: %s\n", path);
    printf("n: %" PRId64 ", nnz: %" PRId64 "\n", a.n, a.nnz);
    printf("threads: %d\n", options.threads);
    printf("build features: METIS %s, SCOTCH %s, SPRAL scaling %s, CBLAS %s\n",
           stats.build_has_metis ? "on" : "off",
           stats.build_has_scotch ? "on" : "off",
           stats.build_has_spral_scaling ? "on" : "off",
           stats.build_has_cblas ? "on" : "off");
    printf("input index bytes: %d (requested %s)\n",
           input_index.bytes, index_mode_name(input_index_mode));
    printf("requested orientation: %s\n", kls_orientation_name(options.orientation));
    printf("selected orientation: %s\n", kls_orientation_name(stats.selected_orientation));
    printf("ordering: %s\n", kls_ordering_name(stats.selected_ordering));
    printf("requested scale: %s\n", scale_name(options.scale));
    printf("selected scale: %d\n", stats.selected_scale);
    printf("requested pivot tolerance: %.6g\n", options.pivot_tolerance);
    printf("selected pivot tolerance: %.6g\n", stats.selected_pivot_tolerance);
    printf("row refactor control: %s\n", row_refactor_control);
    printf("KLS first-factor control: %s\n", kls_first_factor_control);
    printf("row solve control: %s\n", row_solve_control);
    if (stress_requested) {
      printf("stress diagonal scale: %.6g, column: %" PRId64
             ", entries: %" PRId64 "\n",
             stress_diagonal_scale, stress_diagonal_column, stress_entries);
    }
    printf("requested btf: %s\n", options.use_btf ? "on" : "off");
    printf("selected btf: %s\n", stats.selected_btf ? "on" : "off");
    printf("fast factor: %s\n", options.fast_factor ? "on" : "off");
    printf("static pivoting: %s\n", options.static_pivoting ? "on" : "off");
    printf("selected static pivoting: %s\n",
           stats.selected_static_pivoting ? "on" : "off");
    printf("selected exact matching: %s\n",
           stats.selected_exact_matching ? "on" : "off");
    printf("selected exact matching scaling: %s\n",
           stats.selected_exact_matching_scaling ? "on" : "off");
    printf("selected SPRAL matching: %s\n",
           stats.selected_spral_matching ? "on" : "off");
    printf("analysis: %.6f s\n", stats.analysis_seconds);
    printf("initial factor: %.6f s\n", initial_factor_seconds);
    printf("initial factor path: %s\n",
           kls_factor_path_name(initial_factor_path));
    printf("last factor path: %s\n",
           kls_factor_path_name(stats.last_factor_path));
    printf("last refactor path: %s\n",
           kls_refactor_path_name(stats.last_refactor_path));
    printf("KLS tail mapped columns: last %" PRId64
           ", total %" PRId64 "\n",
           stats.kls_tail_last_mapped_columns,
           stats.kls_tail_mapped_column_count);
    printf("KLS row-up-looking first columns: last %" PRId64
           ", total %" PRId64 "\n",
           stats.kls_first_last_row_uplooking_columns,
           stats.kls_first_row_uplooking_column_count);
    printf("KLS row-up-looking dominant BTF pipeline: last %d"
           ", count %" PRId64 ", block %" PRId64 ", rows %" PRId64
           ", separator %d, no-separator count %" PRId64 "\n",
           stats.kls_first_last_dominant_btf_pipeline,
           stats.kls_first_dominant_btf_pipeline_count,
           stats.kls_first_last_dominant_btf_pipeline_block,
           stats.kls_first_last_dominant_btf_pipeline_rows,
           stats.kls_first_last_dominant_btf_pipeline_has_separator,
           stats.kls_first_dominant_btf_pipeline_without_separator_count);
    printf("KLS first row-refactor seed rows: last %" PRId64
           ", total %" PRId64 "\n",
           stats.kls_first_last_row_refactor_seeded_rows,
           stats.kls_first_row_refactor_seeded_row_count);
    printf("KLS row-up-looking natural pipeline: last %d"
           ", runs %" PRId64 ", rows %" PRId64 ", threads %" PRId64
           ", partial last %d, partial runs %" PRId64
           ", partial rows %" PRId64 ", partial threads %" PRId64
           ", pivot tail %d/%" PRId64 "/%" PRId64
           ", restarts %" PRId64 "/%" PRId64
           ", serial rows %" PRId64
           ", prefix rebuild %d/%" PRId64 "/%" PRId64 "\n",
           stats.kls_first_last_row_pipeline,
           stats.kls_first_row_pipeline_run_count,
           stats.kls_first_last_row_pipeline_rows,
           stats.kls_first_last_row_pipeline_threads,
           stats.kls_first_last_row_pipeline_partial,
           stats.kls_first_row_pipeline_partial_run_count,
           stats.kls_first_last_row_pipeline_partial_rows,
           stats.kls_first_last_row_pipeline_partial_threads,
           stats.kls_first_last_row_pipeline_pivot_tail,
           stats.kls_first_row_pipeline_pivot_tail_run_count,
           stats.kls_first_last_row_pipeline_pivot_tail_rows,
           stats.kls_first_last_row_pipeline_pivot_restarts,
           stats.kls_first_row_pipeline_pivot_restart_count,
           stats.kls_first_last_row_pipeline_pivot_serial_rows,
           stats.kls_first_last_row_pipeline_prefix_panel_rebuild,
           stats.kls_first_row_pipeline_prefix_panel_rebuild_count,
           stats.kls_first_last_row_pipeline_prefix_panel_rebuild_rows);
    printf("KLS row-up-looking active-rank pivot resets: resets %" PRId64
           ", reset rows %" PRId64 ", panel rebuilds %" PRId64
           ", rebuild rows %" PRId64 "\n",
           stats.kls_first_active_rank_pivot_reset_count,
           stats.kls_first_active_rank_pivot_reset_rows,
           stats.kls_first_active_rank_pivot_panel_rebuild_count,
           stats.kls_first_active_rank_pivot_panel_rebuild_rows);
    printf("KLS row-up-looking panel cache: builds %" PRId64
           ", build panels %" PRId64 ", build entries %" PRId64
           ", appends %" PRId64 ", append panels %" PRId64
           ", append entries %" PRId64 "\n",
           stats.kls_first_row_panel_cache_build_count,
           stats.kls_first_row_panel_cache_build_panels,
           stats.kls_first_row_panel_cache_build_entries,
           stats.kls_first_row_panel_cache_append_count,
           stats.kls_first_row_panel_cache_append_panels,
           stats.kls_first_row_panel_cache_append_entries);
    printf("KLS row-up-looking supernode updates: last %d"
           ", runs %" PRId64 ", total groups %" PRId64
           ", total rows %" PRId64 ", last groups %" PRId64
           ", last rows %" PRId64 "\n",
           stats.kls_first_last_row_supernode_update,
           stats.kls_first_row_supernode_update_run_count,
           stats.kls_first_row_supernode_update_groups,
           stats.kls_first_row_supernode_update_rows,
           stats.kls_first_last_row_supernode_update_groups,
           stats.kls_first_last_row_supernode_update_rows);
    printf("KLS row-up-looking supernode panels: last %d"
           ", runs %" PRId64 ", total groups %" PRId64
           ", total rows %" PRId64 ", last groups %" PRId64
           ", last rows %" PRId64 "\n",
           stats.kls_first_last_row_supernode_panel_update,
           stats.kls_first_row_supernode_panel_update_run_count,
           stats.kls_first_row_supernode_panel_update_groups,
           stats.kls_first_row_supernode_panel_update_rows,
           stats.kls_first_last_row_supernode_panel_update_groups,
           stats.kls_first_last_row_supernode_panel_update_rows);
    printf("KLS row-up-looking dynamic column pivots: last %" PRId64
           ", total %" PRId64 "\n",
           stats.kls_first_last_dynamic_column_pivots,
           stats.kls_first_dynamic_column_pivot_count);
    printf("KLS row-up-looking separator dynamic pivots: last %" PRId64
           ", total %" PRId64 ", extent last/total %" PRId64 "/%" PRId64
           ", fallbacks last/total %" PRId64 "/%" PRId64
           ", rejects last/total %" PRId64 "/%" PRId64 "\n",
           stats.kls_first_last_separator_dynamic_column_pivots,
           stats.kls_first_separator_dynamic_column_pivot_count,
           stats.kls_first_last_separator_extent_dynamic_column_pivots,
           stats.kls_first_separator_extent_dynamic_column_pivot_count,
           stats.kls_first_last_separator_dynamic_column_fallbacks,
           stats.kls_first_separator_dynamic_column_fallback_count,
           stats.kls_first_last_separator_dynamic_column_rejects,
           stats.kls_first_separator_dynamic_column_reject_count);
    printf("KLS row-up-looking separator queue: last %d, runs %" PRId64
           ", components private/pipeline %" PRId64 "/%" PRId64
           ", rows private/pipeline %" PRId64 "/%" PRId64
           ", threads/max_rows %" PRId64 "/%" PRId64
           ", work min/max %.9g/%.9g"
           ", partitioned %d/%" PRId64 ", split components %" PRId64 "\n",
           stats.kls_first_last_separator_queue,
           stats.kls_first_separator_queue_run_count,
           stats.kls_first_last_separator_queue_private_components,
           stats.kls_first_last_separator_queue_pipeline_components,
           stats.kls_first_last_separator_queue_private_rows,
           stats.kls_first_last_separator_queue_pipeline_rows,
           stats.kls_first_last_separator_queue_nonempty_threads,
           stats.kls_first_last_separator_queue_max_thread_rows,
           stats.kls_first_last_separator_queue_min_thread_work,
           stats.kls_first_last_separator_queue_max_thread_work,
           stats.kls_first_last_separator_queue_partitioned,
           stats.kls_first_separator_queue_partitioned_count,
           stats.kls_first_last_separator_queue_split_components);
    printf("KLS row-up-looking separator queue executed: last %d"
           ", runs %" PRId64 ", rows private/pipeline %" PRId64 "/%" PRId64
           "\n",
           stats.kls_first_last_separator_queue_executed,
           stats.kls_first_separator_queue_executed_run_count,
           stats.kls_first_last_separator_queue_executed_private_rows,
           stats.kls_first_last_separator_queue_executed_pipeline_rows);
    printf("KLS row-up-looking separator private parallel: last %d"
           ", runs %" PRId64 ", rows %" PRId64 ", threads %" PRId64 "\n",
           stats.kls_first_last_separator_queue_parallel_private,
           stats.kls_first_separator_queue_parallel_private_run_count,
           stats.kls_first_last_separator_queue_parallel_private_rows,
           stats.kls_first_last_separator_queue_parallel_private_threads);
    printf("KLS row-up-looking separator pipeline parallel: last %d"
           ", runs %" PRId64 ", rows %" PRId64 ", threads %" PRId64 "\n",
           stats.kls_first_last_separator_queue_parallel_pipeline,
           stats.kls_first_separator_queue_parallel_pipeline_run_count,
           stats.kls_first_last_separator_queue_parallel_pipeline_rows,
           stats.kls_first_last_separator_queue_parallel_pipeline_threads);
    printf("KLS row-up-looking separator pipeline partial: last %d"
           ", runs %" PRId64 ", rows %" PRId64 ", threads %" PRId64 "\n",
           stats.kls_first_last_separator_queue_pipeline_partial,
           stats.kls_first_separator_queue_pipeline_partial_run_count,
           stats.kls_first_last_separator_queue_pipeline_partial_rows,
           stats.kls_first_last_separator_queue_pipeline_partial_threads);
    printf("KLS row-up-looking separator pipeline wait partial: last %d"
           ", runs %" PRId64 ", rows %" PRId64 ", deps %" PRId64 "\n",
           stats.kls_first_last_separator_queue_pipeline_wait_partial,
           stats.kls_first_separator_queue_pipeline_wait_partial_run_count,
           stats.kls_first_last_separator_queue_pipeline_wait_partial_rows,
           stats.kls_first_last_separator_queue_pipeline_wait_partial_deps);
    printf("KLS row-up-looking separator pipeline supernode update: last %d"
           ", runs %" PRId64 ", groups %" PRId64 ", rows %" PRId64 "\n",
           stats.kls_first_last_separator_queue_pipeline_supernode_update,
           stats.kls_first_separator_queue_pipeline_supernode_update_run_count,
           stats.kls_first_last_separator_queue_pipeline_supernode_update_groups,
           stats.kls_first_last_separator_queue_pipeline_supernode_update_rows);
    printf("KLS row-up-looking separator pipeline supernode panels: last %d"
           ", runs %" PRId64 ", groups %" PRId64 ", rows %" PRId64 "\n",
           stats.kls_first_last_separator_queue_pipeline_supernode_panel_update,
           stats.kls_first_separator_queue_pipeline_supernode_panel_update_run_count,
           stats.kls_first_last_separator_queue_pipeline_supernode_panel_update_groups,
           stats.kls_first_last_separator_queue_pipeline_supernode_panel_update_rows);
    printf("KLS row-up-looking separator pipeline pivot tail: last %d"
           ", runs %" PRId64 ", rows %" PRId64
           ", restarts %" PRId64 ", serial rows %" PRId64 "\n",
           stats.kls_first_last_separator_queue_pipeline_pivot_tail,
           stats.kls_first_separator_queue_pipeline_pivot_tail_run_count,
           stats.kls_first_last_separator_queue_pipeline_pivot_tail_rows,
           stats.kls_first_last_separator_queue_pipeline_pivot_restarts,
           stats.kls_first_last_separator_queue_pipeline_pivot_serial_rows);
    printf("KLS row-up-looking separator pipeline prefix panel rebuilds:"
           " last %d, count %" PRId64 ", rows %" PRId64 "\n",
           stats.kls_first_last_separator_queue_pipeline_prefix_panel_rebuild,
           stats.kls_first_separator_queue_pipeline_prefix_panel_rebuild_count,
           stats.kls_first_last_separator_queue_pipeline_prefix_panel_rebuild_rows);
    printf("KLS row-up-looking parallel BTF blocks: last %" PRId64
           ", total %" PRId64 "\n",
           stats.kls_first_last_parallel_btf_blocks,
           stats.kls_first_parallel_btf_block_count);
    printf("KLS first auto skipped scaled single-block: last %d, total %" PRId64
           "\n",
           stats.kls_first_auto_skipped_scaled_single_block,
           stats.kls_first_auto_skipped_scaled_single_block_count);
    printf("factor avg: %.6f s\n", factor_avg);
    printf("refactor avg: %.6f s\n", refactor_avg);
    printf("solve avg: %.6f s\n", solve_avg);
    printf("transpose solve avg: %.6f s\n", tsolve_avg);
    printf("residual: %.6e, relative: %.6e\n", residual, rel_residual);
    printf("blocks: %" PRId64 ", max block: %" PRId64 "\n",
           stats.nblocks, stats.max_block);
    printf("structural rank: %" PRId64 ", numerical rank: %" PRId64 "\n",
           stats.structural_rank, stats.numerical_rank);
    printf("factor ETree largest block: start %" PRId64 ", size %" PRId64
           ", levels %" PRId64 ", max width %" PRId64 ", edges %" PRId64
           ", roots %" PRId64 ", leaves %" PRId64
           ", max fanout %" PRId64 "\n",
           stats.factor_etree_block_start,
           stats.factor_etree_block_size,
           stats.factor_etree_levels,
           stats.factor_etree_max_width,
           stats.factor_etree_edges,
           stats.factor_etree_root_columns,
           stats.factor_etree_leaf_columns,
           stats.factor_etree_max_fanout);
    printf("separator queues: rows %" PRId64 ", global [%" PRId64
           ", %" PRId64 "), threads %" PRId64
           ", components %" PRId64 " (private %" PRId64
           ", pipeline %" PRId64 "), row split %" PRId64 "/%" PRId64
           ", max private/pipeline %" PRId64 "/%" PRId64 "\n",
           stats.separator_analyzed_rows,
           stats.separator_global_begin,
           stats.separator_global_end,
           stats.separator_thread_count,
           stats.separator_component_count,
           stats.separator_private_components,
           stats.separator_pipeline_components,
           stats.separator_private_rows,
           stats.separator_pipeline_rows,
           stats.separator_private_max_rows,
           stats.separator_pipeline_max_rows);
    printf("off-diagonal pivots: %" PRId64 ", reallocations: %" PRId64 "\n",
           stats.offdiag_pivots, stats.reallocations);
    printf("fast rejected pivot: %" PRId64 ", original column: %" PRId64
           ", row: %" PRId64 ", |L|: %.6g, |pivot|: %.6g"
           ", |candidate|: %.6g, best tail row: %" PRId64
           ", best tail |candidate|: %.6g, tail candidates: %" PRId64 "\n",
           stats.fast_rejected_pivot, stats.fast_rejected_pivot_col,
           stats.fast_rejected_row, stats.fast_rejected_multiplier_abs,
           stats.fast_rejected_pivot_abs,
           stats.fast_rejected_candidate_abs,
           stats.fast_rejected_tail_candidate_row,
           stats.fast_rejected_tail_candidate_abs,
           stats.fast_rejected_tail_candidate_count);
    printf("fast rejected row-tail repair: ready %d, candidate position %" PRId64
           "\n",
           stats.fast_rejected_tail_repair_ready,
           stats.fast_rejected_tail_candidate_position);
    printf("fast repaired pivot: row %" PRId64 ", matches tail candidate %d"
           ", first changed pivot %" PRId64 ", prefix changes %" PRId64
           ", suffix changes %" PRId64 ", tail restart ready %d"
           ", block work %.6g, tail columns %" PRId64
           ", tail work %.6g, saved work %.6g"
           ", overcompute columns %" PRId64 ", overcompute work %.6g"
           ", skipped columns %" PRId64 ", skipped work %.6g"
           ", exact mask %d, etree mask %d\n",
           stats.fast_repaired_pivot_row,
           stats.fast_repaired_pivot_matches_tail_candidate,
           stats.fast_repaired_first_changed_pivot,
           stats.fast_repaired_prefix_changed_pivots,
           stats.fast_repaired_suffix_changed_pivots,
           stats.fast_repaired_tail_restart_ready,
           stats.fast_repaired_block_work,
           stats.fast_repaired_tail_restart_columns,
           stats.fast_repaired_tail_restart_work,
           stats.fast_repaired_tail_restart_saved_work,
           stats.fast_repaired_tail_restart_overcompute_columns,
           stats.fast_repaired_tail_restart_overcompute_work,
           stats.fast_repaired_tail_restart_skipped_columns,
           stats.fast_repaired_tail_restart_skipped_work,
           stats.fast_repaired_tail_restart_exact_mask,
           stats.fast_repaired_tail_restart_etree_mask);
    printf("fast block restarts: %d, KLS block restarts: %d"
           ", KLS rebuild restarts: %d"
           ", tail restarts: %d"
           ", KLS block repair row pipeline %d/%" PRId64
           " rows %" PRId64 ", active threads %" PRId64
           ", prefix rows %" PRId64 ", suffix rows %" PRId64
           ", gap rows %" PRId64
           ", etree tail %d/%" PRId64 " rows %" PRId64
           ", etree gaps %" PRId64 ", etree exact mask %d"
           ", etree ready %d/%" PRId64 " rows %" PRId64
           ", ready threads %" PRId64
           ", etree prefactor %d/%" PRId64 " rows %" PRId64
           ", threads %" PRId64 ", wait rows/deps %" PRId64 "/%" PRId64
           ", separator tail scope %d/%" PRId64 " rows %" PRId64
           ", separator queue %d/%" PRId64 " private/pipeline rows %" PRId64
           "/%" PRId64 ", private threads %" PRId64
           ", partitioned %d split components %" PRId64
           ", pivot-tail rows %" PRId64 ", pivot restarts %" PRId64
           ", supernode groups %" PRId64 ", supernode rows %" PRId64
           ", panel groups %" PRId64 ", panel rows %" PRId64
           ", offdiag suffix refresh last %d, suffix refreshes %" PRId64
           ", full refreshes %" PRId64
           ", parallel tail blocks %" PRId64 "\n",
           stats.fast_block_restarts, stats.fast_kls_block_restarts,
           stats.fast_kls_rebuild_restarts,
           stats.fast_tail_restarts,
           stats.fast_kls_block_restart_last_row_pipeline,
           stats.fast_kls_block_restart_row_pipeline_count,
           stats.fast_kls_block_restart_last_row_pipeline_rows,
           stats.fast_kls_block_restart_last_row_pipeline_threads,
           stats.fast_kls_block_restart_last_row_pipeline_prefix_rows,
           stats.fast_kls_block_restart_last_row_pipeline_suffix_rows,
           stats.fast_kls_block_restart_last_row_pipeline_gap_rows,
           stats.fast_kls_block_restart_last_row_pipeline_etree_tail,
           stats.fast_kls_block_restart_row_pipeline_etree_tail_count,
           stats.fast_kls_block_restart_last_row_pipeline_etree_tail_rows,
           stats.fast_kls_block_restart_last_row_pipeline_etree_tail_gap_rows,
           stats.fast_kls_block_restart_last_row_pipeline_etree_tail_exact_mask,
           stats.fast_kls_block_restart_last_row_pipeline_etree_ready,
           stats.fast_kls_block_restart_row_pipeline_etree_ready_count,
           stats.fast_kls_block_restart_last_row_pipeline_etree_ready_rows,
           stats.fast_kls_block_restart_last_row_pipeline_etree_ready_threads,
           stats.fast_kls_block_restart_last_row_pipeline_etree_prefactor,
           stats.fast_kls_block_restart_row_pipeline_etree_prefactor_count,
           stats.fast_kls_block_restart_last_row_pipeline_etree_prefactor_rows,
           stats.fast_kls_block_restart_last_row_pipeline_etree_prefactor_threads,
           stats.fast_kls_block_restart_last_row_pipeline_etree_prefactor_wait_rows,
           stats.fast_kls_block_restart_last_row_pipeline_etree_prefactor_wait_deps,
           stats.fast_kls_block_restart_last_row_pipeline_separator_tail_scope,
           stats.fast_kls_block_restart_row_pipeline_separator_tail_scope_count,
           stats.fast_kls_block_restart_last_row_pipeline_separator_tail_scope_rows,
           stats.fast_kls_block_restart_last_row_pipeline_separator_queue,
           stats.fast_kls_block_restart_row_pipeline_separator_queue_count,
           stats.fast_kls_block_restart_last_row_pipeline_separator_private_rows,
           stats.fast_kls_block_restart_last_row_pipeline_separator_pipeline_rows,
           stats.fast_kls_block_restart_last_row_pipeline_separator_private_threads,
           stats.fast_kls_block_restart_last_row_pipeline_separator_partitioned,
           stats.fast_kls_block_restart_last_row_pipeline_separator_split_components,
           stats.fast_kls_block_restart_last_row_pipeline_pivot_tail_rows,
           stats.fast_kls_block_restart_last_row_pipeline_pivot_restarts,
           stats.fast_kls_block_restart_last_row_pipeline_supernode_update_groups,
           stats.fast_kls_block_restart_last_row_pipeline_supernode_update_rows,
           stats.fast_kls_block_restart_last_row_pipeline_supernode_panel_update_groups,
           stats.fast_kls_block_restart_last_row_pipeline_supernode_panel_update_rows,
           stats.fast_repaired_last_offdiag_suffix_refresh,
           stats.fast_repaired_offdiag_suffix_refresh_count,
           stats.fast_repaired_offdiag_full_refresh_count,
           stats.fast_repaired_parallel_tail_blocks);
    printf("fast rejected block: start %" PRId64 ", size %" PRId64
           ", suffix %" PRId64 ", descendants %" PRId64
           ", descendant work %.6g, row tail %" PRId64
           ", row-tail work %.6g, group tail %" PRId64
           ", group rows %" PRId64 ", group work %.6g, etree %" PRId64
           ", etree work %.6g, pivoting tail %" PRId64
           ", pivoting-tail work %.6g, first %" PRId64
           ", last %" PRId64 ", contains reject %d, topological %d"
           ", seed %" PRId64 ", row seed %" PRId64
           ", block seed %" PRId64
           ", contiguous %d, suffix exact %d"
           ", gaps %" PRId64 ", suffix overcompute columns %" PRId64
           ", suffix overcompute work %.6g, etree edges %" PRId64
           ", etree roots %" PRId64 ", etree leaves %" PRId64
           ", etree max fanout %" PRId64 ", etree levels %" PRId64
           ", etree max width %" PRId64 ", refresh state %d"
           ", prefix refresh columns %" PRId64
           ", prefix refreshes %" PRId64 "\n",
           stats.fast_rejected_block_start,
           stats.fast_rejected_block_size,
           stats.fast_rejected_suffix_columns,
           stats.fast_rejected_descendant_columns,
           stats.fast_rejected_descendant_work,
           stats.fast_rejected_row_tail_columns,
           stats.fast_rejected_row_tail_work,
           stats.fast_rejected_group_tail_groups,
           stats.fast_rejected_group_tail_rows,
           stats.fast_rejected_group_tail_work,
           stats.fast_rejected_etree_columns,
           stats.fast_rejected_etree_work,
           stats.fast_rejected_pivoting_tail_columns,
           stats.fast_rejected_pivoting_tail_work,
           stats.fast_rejected_pivoting_tail_first,
           stats.fast_rejected_pivoting_tail_last,
           stats.fast_rejected_pivoting_tail_contains_reject,
           stats.fast_rejected_pivoting_tail_topological,
           stats.fast_rejected_pivoting_tail_seed_columns,
           stats.fast_rejected_pivoting_tail_row_seed_columns,
           stats.fast_rejected_pivoting_tail_block_seed_columns,
           stats.fast_rejected_pivoting_tail_contiguous,
           stats.fast_rejected_pivoting_tail_suffix_exact,
           stats.fast_rejected_pivoting_tail_gap_columns,
           stats.fast_rejected_pivoting_tail_suffix_overcompute_columns,
           stats.fast_rejected_pivoting_tail_suffix_overcompute_work,
           stats.fast_rejected_pivoting_tail_etree_edges,
           stats.fast_rejected_pivoting_tail_etree_roots,
           stats.fast_rejected_pivoting_tail_etree_leaves,
           stats.fast_rejected_pivoting_tail_etree_max_fanout,
           stats.fast_rejected_pivoting_tail_etree_levels,
           stats.fast_rejected_pivoting_tail_etree_max_width,
           stats.fast_rejected_refresh_state,
           stats.fast_rejected_prefix_refresh_columns,
           stats.fast_rejected_prefix_refresh_count);
    printf("refactor dependency levels: %" PRId64
           ", max width: %" PRId64 ", edges: %" PRId64 "\n",
           stats.refactor_dependency_levels,
           stats.refactor_dependency_max_width,
           stats.refactor_dependency_edges);
    printf("refactor dependency roots: %" PRId64
           ", leaves: %" PRId64 ", max fanout: %" PRId64 "\n",
           stats.refactor_dependency_root_columns,
           stats.refactor_dependency_leaf_columns,
           stats.refactor_dependency_max_fanout);
    printf("refactor dependency max column work: %.6g"
           ", pipeline max column work: %.6g\n",
           stats.refactor_dependency_max_column_work,
           stats.refactor_dependency_pipeline_max_column_work);
    printf("refactor supernode candidates: %" PRId64
           ", rows: %" PRId64 ", max width: %" PRId64
           ", dense entries: %.6g, trailing entries: %.6g"
           ", consumer runs: %" PRId64 ", consumer rows: %" PRId64
           ", consumer max width: %" PRId64
           ", suffix runs: %" PRId64
           ", consumer L entries: %.6g"
           ", consumer internal entries: %.6g"
           ", consumer panels: %" PRId64 "/%" PRId64
           ", reused runs: %" PRId64 "/%" PRId64
           ", max panel reuse: %" PRId64 "/%" PRId64
           ", reused entries: %.6g/%.6g"
           ", U patterns: %" PRId64 "/%" PRId64
           ", U pattern right entries: %" PRId64
           ", U pattern internal entries: %.6g"
           ", U value entries: %" PRId64 "/%" PRId64
           ", U value writes: %" PRId64 "/%" PRId64
           ", ragged L panels: %" PRId64
           ", ragged L entries: %" PRId64 "/%" PRId64
           ", ragged L pruned: %" PRId64 "/%" PRId64
           ", ragged L updates: %" PRId64 "/%" PRId64 "/%" PRId64
           ", ragged L misses: %" PRId64 "/%" PRId64 "/%" PRId64
           "/%" PRId64 "/%" PRId64
           ", ragged L disabled: %d/%" PRId64
           ", panels: %" PRId64 ", used panels: %" PRId64
           ", plan cached panels: %" PRId64 "/%" PRId64
           ", strict plan cached panels: %" PRId64 "/%" PRId64
           "/%" PRId64
           ", plan deferred output: %" PRId64 "/%" PRId64 "/%" PRId64
           ", batch deferred output: %" PRId64 "/%" PRId64 "/%" PRId64
           ", cached probe disabled: %d/%" PRId64 "\n",
           stats.refactor_supernode_candidate_count,
           stats.refactor_supernode_candidate_rows,
           stats.refactor_supernode_candidate_max_width,
           stats.refactor_supernode_candidate_dense_entries,
           stats.refactor_supernode_candidate_trailing_entries,
           stats.refactor_supernode_consumer_run_count,
           stats.refactor_supernode_consumer_run_rows,
           stats.refactor_supernode_consumer_run_max_width,
           stats.refactor_supernode_consumer_suffix_count,
           stats.refactor_supernode_consumer_l_entries,
           stats.refactor_supernode_consumer_internal_entries,
           stats.refactor_supernode_consumer_panel_count,
           stats.refactor_supernode_consumer_reused_panel_count,
           stats.refactor_supernode_consumer_reused_run_count,
           stats.refactor_supernode_consumer_reused_run_rows,
           stats.refactor_supernode_consumer_panel_max_runs,
           stats.refactor_supernode_consumer_panel_max_rows,
           stats.refactor_supernode_consumer_reused_l_entries,
           stats.refactor_supernode_consumer_reused_internal_entries,
           stats.refactor_u_supernode_pattern_count,
           stats.refactor_u_supernode_pattern_rows,
           stats.refactor_u_supernode_pattern_right_entries,
           stats.refactor_u_supernode_pattern_internal_entries,
           stats.refactor_u_supernode_value_dense_entries,
           stats.refactor_u_supernode_value_right_entries,
           stats.refactor_last_u_supernode_value_dense_writes,
           stats.refactor_last_u_supernode_value_right_writes,
           stats.refactor_u_supernode_l_panel_count,
           stats.refactor_u_supernode_l_dense_entries,
           stats.refactor_u_supernode_l_trailing_entries,
           stats.refactor_u_supernode_l_prune_count,
           stats.refactor_u_supernode_l_pruned_panels,
           stats.refactor_last_u_supernode_l_update_runs,
           stats.refactor_last_u_supernode_l_update_rows,
           stats.refactor_last_u_supernode_l_update_entries,
           stats.refactor_last_u_supernode_l_probe_attempts,
           stats.refactor_last_u_supernode_l_panel_misses,
           stats.refactor_last_u_supernode_l_short_rejects,
           stats.refactor_last_u_supernode_l_stream_rejects,
           stats.refactor_last_u_supernode_l_work_rejects,
           stats.refactor_u_supernode_l_exec_disabled,
           stats.refactor_u_supernode_l_exec_disable_count,
           stats.refactor_supernode_panel_count,
           stats.refactor_supernode_panel_used_count,
           stats.refactor_supernode_consumer_plan_cached_panel_count,
           stats.refactor_supernode_consumer_plan_cached_panel_rows,
           stats.refactor_supernode_consumer_plan_strict_cached_panel_count,
           stats.refactor_supernode_consumer_plan_strict_cached_panel_rows,
           stats.refactor_supernode_consumer_plan_strict_cached_panel_trailing_entries,
           stats.refactor_supernode_consumer_plan_deferred_columns,
           stats.refactor_supernode_consumer_plan_deferred_entries,
           stats.refactor_supernode_consumer_plan_deferred_unique_rows,
           stats.refactor_supernode_consumer_plan_batch_deferred_columns,
           stats.refactor_supernode_consumer_plan_batch_deferred_entries,
           stats.refactor_supernode_consumer_plan_batch_deferred_unique_rows,
           stats.refactor_supernode_cached_probe_disabled,
           stats.refactor_supernode_cached_probe_disable_count);
    printf("refactor BTF scalar producer runs: last %" PRId64
           "/%" PRId64 "/%" PRId64 ", max rows %" PRId64
           ", cumulative %" PRId64 "/%" PRId64 "/%" PRId64
           ", max rows %" PRId64
           ", exec %" PRId64 "/%" PRId64 "/%" PRId64
           ", exec max rows %" PRId64
           ", cumulative exec %" PRId64 "/%" PRId64 "/%" PRId64
           ", exec max rows %" PRId64
           ", groups built %d"
           ", groups %" PRId64 "/%" PRId64
           ", multi %" PRId64 "/%" PRId64
           ", rows %" PRId64 " reused %" PRId64
           ", entries %" PRId64 " reused %" PRId64
           ", max currents %" PRId64 ", max rows %" PRId64
           ", waits %" PRId64 "/%" PRId64 "/%" PRId64
           ", overlaps %" PRId64 "/%" PRId64 "/%" PRId64
           ", max live %" PRId64
           ", cumulative waits %" PRId64 "/%" PRId64 "/%" PRId64
           ", cumulative overlaps %" PRId64 "/%" PRId64 "/%" PRId64
           ", cumulative max live %" PRId64
           ", claim surface %" PRId64 "/%" PRId64 "/%" PRId64
           ", claims %" PRId64 "/%" PRId64 "/%" PRId64
           ", prefix ready %" PRId64 "/%" PRId64 "/%" PRId64
           ", run ready %" PRId64 "/%" PRId64 "/%" PRId64
           ", wake armed %" PRId64 "/%" PRId64 "/%" PRId64
           ", wake ready %" PRId64 "/%" PRId64 "/%" PRId64
           ", cumulative claim surface %" PRId64 "/%" PRId64 "/%" PRId64
           ", cumulative claims %" PRId64 "/%" PRId64 "/%" PRId64
           ", cumulative prefix ready %" PRId64 "/%" PRId64 "/%" PRId64
           ", cumulative run ready %" PRId64 "/%" PRId64 "/%" PRId64
           ", cumulative wake armed %" PRId64 "/%" PRId64 "/%" PRId64
           ", cumulative wake ready %" PRId64 "/%" PRId64 "/%" PRId64 "\n",
           stats.refactor_last_btf_scalar_run_candidates,
           stats.refactor_last_btf_scalar_run_rows,
           stats.refactor_last_btf_scalar_run_entries,
           stats.refactor_last_btf_scalar_run_max_rows,
           stats.refactor_btf_scalar_run_candidate_count,
           stats.refactor_btf_scalar_run_rows,
           stats.refactor_btf_scalar_run_entries,
           stats.refactor_btf_scalar_run_max_rows,
           stats.refactor_last_btf_scalar_run_exec_runs,
           stats.refactor_last_btf_scalar_run_exec_rows,
           stats.refactor_last_btf_scalar_run_exec_entries,
           stats.refactor_last_btf_scalar_run_exec_max_rows,
           stats.refactor_btf_scalar_run_exec_count,
           stats.refactor_btf_scalar_run_exec_rows,
           stats.refactor_btf_scalar_run_exec_entries,
           stats.refactor_btf_scalar_run_exec_max_rows,
           stats.refactor_btf_scalar_run_group_built,
           stats.refactor_btf_scalar_run_group_count,
           stats.refactor_btf_scalar_run_group_current_total,
           stats.refactor_btf_scalar_run_group_multi_count,
           stats.refactor_btf_scalar_run_group_multi_current_total,
           stats.refactor_btf_scalar_run_group_rows,
           stats.refactor_btf_scalar_run_group_reused_rows,
           stats.refactor_btf_scalar_run_group_entries,
           stats.refactor_btf_scalar_run_group_reused_entries,
           stats.refactor_btf_scalar_run_group_max_currents,
           stats.refactor_btf_scalar_run_group_max_rows,
           stats.refactor_last_btf_scalar_run_group_waits,
           stats.refactor_last_btf_scalar_run_group_wait_rows,
           stats.refactor_last_btf_scalar_run_group_wait_entries,
           stats.refactor_last_btf_scalar_run_group_overlap_waits,
           stats.refactor_last_btf_scalar_run_group_overlap_rows,
           stats.refactor_last_btf_scalar_run_group_overlap_entries,
           stats.refactor_last_btf_scalar_run_group_max_live,
           stats.refactor_btf_scalar_run_group_wait_count,
           stats.refactor_btf_scalar_run_group_wait_rows,
           stats.refactor_btf_scalar_run_group_wait_entries,
           stats.refactor_btf_scalar_run_group_overlap_count,
           stats.refactor_btf_scalar_run_group_overlap_rows,
           stats.refactor_btf_scalar_run_group_overlap_entries,
           stats.refactor_btf_scalar_run_group_max_live,
           stats.refactor_last_btf_scalar_run_group_claim_surface_triggers,
           stats.refactor_last_btf_scalar_run_group_claim_surface_groups,
           stats.refactor_last_btf_scalar_run_group_claim_surface_currents,
           stats.refactor_last_btf_scalar_run_group_claim_triggers,
           stats.refactor_last_btf_scalar_run_group_claim_groups,
           stats.refactor_last_btf_scalar_run_group_claim_currents,
           stats.refactor_last_btf_scalar_run_group_prefix_ready_currents,
           stats.refactor_last_btf_scalar_run_group_prefix_ready_rows,
           stats.refactor_last_btf_scalar_run_group_prefix_ready_entries,
           stats.refactor_last_btf_scalar_run_group_run_ready_currents,
           stats.refactor_last_btf_scalar_run_group_run_ready_rows,
           stats.refactor_last_btf_scalar_run_group_run_ready_entries,
           stats.refactor_last_btf_scalar_run_group_wake_armed_currents,
           stats.refactor_last_btf_scalar_run_group_wake_armed_rows,
           stats.refactor_last_btf_scalar_run_group_wake_armed_entries,
           stats.refactor_last_btf_scalar_run_group_wake_ready_currents,
           stats.refactor_last_btf_scalar_run_group_wake_ready_rows,
           stats.refactor_last_btf_scalar_run_group_wake_ready_entries,
           stats.refactor_btf_scalar_run_group_claim_surface_trigger_count,
           stats.refactor_btf_scalar_run_group_claim_surface_group_count,
           stats.refactor_btf_scalar_run_group_claim_surface_current_count,
           stats.refactor_btf_scalar_run_group_claim_trigger_count,
           stats.refactor_btf_scalar_run_group_claim_group_count,
           stats.refactor_btf_scalar_run_group_claim_current_count,
           stats.refactor_btf_scalar_run_group_prefix_ready_current_count,
           stats.refactor_btf_scalar_run_group_prefix_ready_rows,
           stats.refactor_btf_scalar_run_group_prefix_ready_entries,
           stats.refactor_btf_scalar_run_group_run_ready_current_count,
           stats.refactor_btf_scalar_run_group_run_ready_rows,
           stats.refactor_btf_scalar_run_group_run_ready_entries,
           stats.refactor_btf_scalar_run_group_wake_armed_current_count,
           stats.refactor_btf_scalar_run_group_wake_armed_rows,
           stats.refactor_btf_scalar_run_group_wake_armed_entries,
           stats.refactor_btf_scalar_run_group_wake_ready_current_count,
           stats.refactor_btf_scalar_run_group_wake_ready_rows,
           stats.refactor_btf_scalar_run_group_wake_ready_entries);
    printf("refactor supernode consumer plan: panels %" PRId64 "/%" PRId64
           ", runs %" PRId64 ", rows %" PRId64
           ", positioned %" PRId64 "/%" PRId64
           ", max panel %" PRId64 "/%" PRId64
           ", entries %.6g/%.6g, bytes %" PRId64
           ", small %" PRId64 "/%" PRId64
           ", batch %" PRId64 "/%" PRId64 "/%" PRId64
           ", batch small %" PRId64 "/%" PRId64
           ", columns %" PRId64 ", max column %" PRId64 "/%" PRId64
           ", column batch %" PRId64 "/%" PRId64 "/%" PRId64
           ", column batch small %" PRId64 "/%" PRId64
           ", column exact-shape batch %" PRId64 "/%" PRId64 "/%" PRId64
           ", column exact-shape small %" PRId64 "/%" PRId64
           ", column exact-shape max %" PRId64
           ", exact-shape batch %" PRId64 "/%" PRId64 "/%" PRId64
           ", exact-shape small %" PRId64 "/%" PRId64
           ", exact-shape max %" PRId64
           ", first-dep exact-shape %" PRId64 "/%" PRId64 "/%" PRId64
           ", first-dep small %" PRId64 "/%" PRId64
           ", first-dep max %" PRId64
           ", advance exact-shape %" PRId64 "/%" PRId64 "/%" PRId64
           ", advance deps/work %" PRId64 "/%.9g"
           ", advance max %" PRId64 "/%.9g"
           ", prefix-advance batch %" PRId64 "/%" PRId64 "/%" PRId64
           ", prefix-advance deps/work %" PRId64 "/%.9g"
           ", prefix-advance max %" PRId64 "/%" PRId64 "/%.9g"
           ", bounded advance limit %" PRId64
           " groups %" PRId64 "/%" PRId64 "/%" PRId64
           " deps/update/work %" PRId64 "/%" PRId64 "/%.9g"
           " max %" PRId64 "/%.9g"
           " payoff %" PRId64 "/%" PRId64 "/%" PRId64
           "/%" PRId64 " advance %.9g max ratio %.9g"
           ", alg5 split panels %" PRId64 "/%" PRId64 "/%" PRId64
           " max %" PRId64
           ", alg5 split runs %" PRId64 "/%" PRId64
           ", prefix %" PRId64 "/%" PRId64 " work %.9g"
           ", crossing %" PRId64 "/%" PRId64
           ", alg5 advance %" PRId64 "/%" PRId64 "/%" PRId64
           " work %.9g"
           ", alg5 payoff %" PRId64 "/%" PRId64 " work %.9g/%.9g"
           ", alg5 prefactor %" PRId64 "/%" PRId64
           " cumulative %" PRId64 "/%" PRId64
           ", alg5 panels %" PRId64 "/%" PRId64 "/%" PRId64
           " deps %" PRId64 " work %.9g/%.9g max %" PRId64
           ", alg5 panel payoff %" PRId64 "/%" PRId64 "/%" PRId64
           " work %.9g/%.9g"
           ", alg5 subset payoff %" PRId64 "/%" PRId64 "/%" PRId64
           " work %.9g/%.9g max %" PRId64
           " currents %" PRId64 "/%" PRId64 "/%" PRId64
           " groups %" PRId64 "/%" PRId64 "/%" PRId64 "/%" PRId64
           " currents %" PRId64 "/%" PRId64 "/%" PRId64
           " workspace %" PRId64 "/%" PRId64
           " advance %" PRId64 "/%" PRId64 "/%" PRId64
           " advance share %" PRId64 "/%" PRId64 "/%" PRId64 "/%" PRId64
           "/%" PRId64
           " suffix %" PRId64 "/%" PRId64 "/%" PRId64 "/%" PRId64
           " suffix share %" PRId64 "/%" PRId64 "/%" PRId64 "/%" PRId64
           "/%" PRId64
           " advance slots %" PRId64 "/%" PRId64 "/%" PRId64
           " target %" PRId64 "/%" PRId64 "/%" PRId64
           " target slots %" PRId64 "/%" PRId64 "/%" PRId64
           " pattern %" PRId64 "/%" PRId64
           " state rows %" PRId64 "/%" PRId64
           " state span %" PRId64 "/%" PRId64
           " runtime %" PRId64 "/%" PRId64 "/%" PRId64
           " slot accum %" PRId64 "/%" PRId64 "/%" PRId64 "/%" PRId64
           " prefix prep %" PRId64 "/%" PRId64 "/%" PRId64 "/%" PRId64
           " current-state seed %" PRId64 "/%" PRId64 "/%" PRId64
           " final trigger %" PRId64 "/%" PRId64 "/%" PRId64 "/%" PRId64
           " suffix advance %" PRId64 "/%" PRId64 "/%" PRId64 "/%" PRId64
           ", group L %d/%d %" PRId64 "/%" PRId64 "/%" PRId64
           "/%" PRId64 "/%" PRId64 "/%" PRId64 "/%" PRId64
           "/%" PRId64 "/%" PRId64
          ", group L state %d/%d focus %d stored %" PRId64 "/%" PRId64
          "/%" PRId64 "/%" PRId64
          " candidate %" PRId64 "/%" PRId64 "/%" PRId64
          " selected %" PRId64 "/%" PRId64 "/%" PRId64 "/%" PRId64
           ", group L batch candidates %" PRId64 "/%" PRId64 "/%" PRId64
           "/%" PRId64 "/%" PRId64
           ", group L batch advance %" PRId64 "/%" PRId64 "/%" PRId64
           " work %.9g max %.9g"
           ", group L batch payoff %" PRId64 "/%" PRId64 "/%" PRId64
           "/%" PRId64 " advance %.9g max ratio %.9g"
           ", group L writes %" PRId64 "/%" PRId64 "/%" PRId64
           " cumulative %" PRId64 "/%" PRId64 "/%" PRId64
           ", group L updates %" PRId64 "/%" PRId64 "/%" PRId64
           " cumulative %" PRId64 "/%" PRId64 "/%" PRId64 "\n",
           stats.refactor_supernode_consumer_plan_panel_count,
           stats.refactor_supernode_consumer_plan_reused_panel_count,
           stats.refactor_supernode_consumer_plan_run_count,
           stats.refactor_supernode_consumer_plan_run_rows,
           stats.refactor_supernode_consumer_plan_positioned_run_count,
           stats.refactor_supernode_consumer_plan_positioned_run_rows,
           stats.refactor_supernode_consumer_plan_max_panel_runs,
           stats.refactor_supernode_consumer_plan_max_panel_rows,
           stats.refactor_supernode_consumer_plan_l_entries,
           stats.refactor_supernode_consumer_plan_internal_entries,
           stats.refactor_supernode_consumer_plan_bytes,
           stats.refactor_supernode_consumer_plan_small_run_count,
           stats.refactor_supernode_consumer_plan_small_run_rows,
           stats.refactor_supernode_consumer_plan_batch_panel_count,
           stats.refactor_supernode_consumer_plan_batch_run_count,
           stats.refactor_supernode_consumer_plan_batch_run_rows,
           stats.refactor_supernode_consumer_plan_batch_small_run_count,
           stats.refactor_supernode_consumer_plan_batch_small_run_rows,
           stats.refactor_supernode_consumer_plan_column_count,
           stats.refactor_supernode_consumer_plan_max_column_runs,
           stats.refactor_supernode_consumer_plan_max_column_rows,
           stats.refactor_supernode_consumer_plan_column_batch_count,
           stats.refactor_supernode_consumer_plan_column_batch_run_count,
           stats.refactor_supernode_consumer_plan_column_batch_run_rows,
           stats.refactor_supernode_consumer_plan_column_batch_small_run_count,
           stats.refactor_supernode_consumer_plan_column_batch_small_run_rows,
           stats.refactor_supernode_consumer_plan_column_shape_batch_count,
           stats.refactor_supernode_consumer_plan_column_shape_batch_run_count,
           stats.refactor_supernode_consumer_plan_column_shape_batch_run_rows,
           stats.refactor_supernode_consumer_plan_column_shape_batch_small_run_count,
           stats.refactor_supernode_consumer_plan_column_shape_batch_small_run_rows,
           stats.refactor_supernode_consumer_plan_column_shape_batch_max_runs,
           stats.refactor_supernode_consumer_plan_shape_batch_count,
           stats.refactor_supernode_consumer_plan_shape_batch_run_count,
           stats.refactor_supernode_consumer_plan_shape_batch_run_rows,
           stats.refactor_supernode_consumer_plan_shape_batch_small_run_count,
           stats.refactor_supernode_consumer_plan_shape_batch_small_run_rows,
           stats.refactor_supernode_consumer_plan_shape_batch_max_runs,
           stats.refactor_supernode_consumer_plan_first_dep_shape_batch_count,
           stats.refactor_supernode_consumer_plan_first_dep_shape_batch_run_count,
           stats.refactor_supernode_consumer_plan_first_dep_shape_batch_run_rows,
           stats.refactor_supernode_consumer_plan_first_dep_shape_batch_small_run_count,
           stats.refactor_supernode_consumer_plan_first_dep_shape_batch_small_run_rows,
           stats.refactor_supernode_consumer_plan_first_dep_shape_batch_max_runs,
           stats.refactor_supernode_consumer_plan_shape_batch_advance_count,
           stats.refactor_supernode_consumer_plan_shape_batch_advance_run_count,
           stats.refactor_supernode_consumer_plan_shape_batch_advance_run_rows,
           stats.refactor_supernode_consumer_plan_shape_batch_advance_dep_count,
           stats.refactor_supernode_consumer_plan_shape_batch_advance_work,
           stats.refactor_supernode_consumer_plan_shape_batch_advance_max_deps,
           stats.refactor_supernode_consumer_plan_shape_batch_advance_max_work,
           stats.refactor_supernode_consumer_plan_prefix_advance_batch_count,
           stats
             .refactor_supernode_consumer_plan_prefix_advance_batch_run_count,
           stats
             .refactor_supernode_consumer_plan_prefix_advance_batch_run_rows,
           stats
             .refactor_supernode_consumer_plan_prefix_advance_batch_dep_count,
           stats.refactor_supernode_consumer_plan_prefix_advance_batch_work,
           stats
             .refactor_supernode_consumer_plan_prefix_advance_batch_max_runs,
           stats
             .refactor_supernode_consumer_plan_prefix_advance_batch_max_deps,
           stats
             .refactor_supernode_consumer_plan_prefix_advance_batch_max_work,
           stats
             .refactor_supernode_consumer_plan_shape_bounded_advance_dep_limit,
           stats.refactor_supernode_consumer_plan_shape_bounded_advance_count,
           stats
             .refactor_supernode_consumer_plan_shape_bounded_advance_run_count,
           stats
             .refactor_supernode_consumer_plan_shape_bounded_advance_run_rows,
           stats
             .refactor_supernode_consumer_plan_shape_bounded_advance_dep_count,
           stats
             .refactor_supernode_consumer_plan_shape_bounded_advance_update_entries,
           stats.refactor_supernode_consumer_plan_shape_bounded_advance_work,
           stats
             .refactor_supernode_consumer_plan_shape_bounded_advance_max_run_deps,
           stats
             .refactor_supernode_consumer_plan_shape_bounded_advance_max_run_work,
           stats
             .refactor_supernode_consumer_plan_shape_bounded_advance_payoff_count,
           stats
             .refactor_supernode_consumer_plan_shape_bounded_advance_payoff_run_count,
           stats
             .refactor_supernode_consumer_plan_shape_bounded_advance_payoff_run_rows,
           stats
             .refactor_supernode_consumer_plan_shape_bounded_advance_payoff_update_entries,
           stats
             .refactor_supernode_consumer_plan_shape_bounded_advance_payoff_advance_work,
           stats
             .refactor_supernode_consumer_plan_shape_bounded_advance_max_payoff_ratio,
           stats.refactor_supernode_algorithm5_large_panel_count,
           stats.refactor_supernode_algorithm5_large_panel_rows,
           stats.refactor_supernode_algorithm5_large_panel_prefix_rows,
           stats.refactor_supernode_algorithm5_large_panel_max_width,
           stats.refactor_supernode_algorithm5_candidate_run_count,
           stats.refactor_supernode_algorithm5_candidate_run_rows,
           stats.refactor_supernode_algorithm5_prefix_run_count,
           stats.refactor_supernode_algorithm5_prefix_run_rows,
           stats.refactor_supernode_algorithm5_prefix_update_work,
           stats.refactor_supernode_algorithm5_crossing_run_count,
           stats.refactor_supernode_algorithm5_crossing_run_rows,
           stats.refactor_supernode_algorithm5_prefix_advance_run_count,
           stats.refactor_supernode_algorithm5_prefix_advance_run_rows,
           stats.refactor_supernode_algorithm5_prefix_advance_dep_count,
           stats.refactor_supernode_algorithm5_prefix_advance_work,
           stats.refactor_supernode_algorithm5_prefix_payoff_run_count,
           stats.refactor_supernode_algorithm5_prefix_payoff_run_rows,
           stats.refactor_supernode_algorithm5_prefix_payoff_update_work,
           stats.refactor_supernode_algorithm5_prefix_payoff_advance_work,
           stats.refactor_last_egraph_algorithm5_prefactor_columns,
           stats.refactor_last_egraph_algorithm5_prefactor_deps,
           stats.refactor_egraph_algorithm5_prefactor_column_count,
           stats.refactor_egraph_algorithm5_prefactor_dep_count,
           stats.refactor_supernode_algorithm5_prefix_panel_count,
           stats.refactor_supernode_algorithm5_prefix_panel_run_count,
           stats.refactor_supernode_algorithm5_prefix_panel_run_rows,
           stats
             .refactor_supernode_algorithm5_prefix_panel_advance_dep_count,
           stats.refactor_supernode_algorithm5_prefix_panel_update_work,
           stats.refactor_supernode_algorithm5_prefix_panel_advance_work,
           stats.refactor_supernode_algorithm5_prefix_panel_max_runs,
           stats.refactor_supernode_algorithm5_prefix_panel_payoff_count,
           stats.refactor_supernode_algorithm5_prefix_panel_payoff_run_count,
           stats.refactor_supernode_algorithm5_prefix_panel_payoff_run_rows,
           stats
             .refactor_supernode_algorithm5_prefix_panel_payoff_update_work,
           stats
             .refactor_supernode_algorithm5_prefix_panel_payoff_advance_work,
           stats
             .refactor_supernode_algorithm5_prefix_panel_payoff_subset_count,
           stats
             .refactor_supernode_algorithm5_prefix_panel_payoff_subset_run_count,
           stats
             .refactor_supernode_algorithm5_prefix_panel_payoff_subset_run_rows,
           stats
             .refactor_supernode_algorithm5_prefix_panel_payoff_subset_update_work,
           stats
             .refactor_supernode_algorithm5_prefix_panel_payoff_subset_advance_work,
           stats
             .refactor_supernode_algorithm5_prefix_panel_payoff_subset_max_runs,
           stats
             .refactor_supernode_algorithm5_prefix_panel_payoff_subset_current_count,
           stats
             .refactor_supernode_algorithm5_prefix_panel_payoff_subset_multi_current_count,
           stats
             .refactor_supernode_algorithm5_prefix_panel_payoff_subset_max_currents,
           stats.refactor_supernode_algorithm5_payoff_group_count,
           stats.refactor_supernode_algorithm5_payoff_group_prefix_rows,
           stats.refactor_supernode_algorithm5_payoff_group_max_run_rows,
           stats.refactor_supernode_algorithm5_payoff_group_positioned_runs,
           stats.refactor_supernode_algorithm5_payoff_group_current_total,
           stats
             .refactor_supernode_algorithm5_payoff_group_multi_current_count,
           stats.refactor_supernode_algorithm5_payoff_group_max_currents,
           stats.refactor_supernode_algorithm5_payoff_group_workspace_rows,
           stats
             .refactor_supernode_algorithm5_payoff_group_max_workspace_rows,
           stats.refactor_supernode_algorithm5_payoff_group_advance_deps,
           stats
             .refactor_supernode_algorithm5_payoff_group_max_advance_deps,
           stats
             .refactor_supernode_algorithm5_payoff_group_zero_advance_runs,
           stats
             .refactor_supernode_algorithm5_payoff_group_advance_unique_deps,
           stats
             .refactor_supernode_algorithm5_payoff_group_advance_duplicate_deps,
           stats
             .refactor_supernode_algorithm5_payoff_group_advance_shared_deps,
           stats
             .refactor_supernode_algorithm5_payoff_group_advance_max_dep_fanout,
           stats
             .refactor_supernode_algorithm5_payoff_group_advance_duplicate_update_entries,
           stats.refactor_supernode_algorithm5_payoff_group_suffix_deps,
           stats
             .refactor_supernode_algorithm5_payoff_group_max_suffix_deps,
           stats
             .refactor_supernode_algorithm5_payoff_group_suffix_update_entries,
           stats
             .refactor_supernode_algorithm5_payoff_group_max_run_suffix_update_entries,
           stats
             .refactor_supernode_algorithm5_payoff_group_suffix_unique_deps,
           stats
             .refactor_supernode_algorithm5_payoff_group_suffix_duplicate_deps,
           stats
             .refactor_supernode_algorithm5_payoff_group_suffix_shared_deps,
           stats
             .refactor_supernode_algorithm5_payoff_group_suffix_max_dep_fanout,
           stats
             .refactor_supernode_algorithm5_payoff_group_suffix_duplicate_update_entries,
           stats.refactor_supernode_algorithm5_payoff_group_advance_slots,
           stats
             .refactor_supernode_algorithm5_payoff_group_max_advance_slots,
           stats
             .refactor_supernode_algorithm5_payoff_group_max_run_advance_slots,
           stats.refactor_supernode_algorithm5_payoff_group_target_entries,
           stats
             .refactor_supernode_algorithm5_payoff_group_max_target_entries,
           stats
             .refactor_supernode_algorithm5_payoff_group_max_run_target_entries,
           stats.refactor_supernode_algorithm5_payoff_group_target_slots,
           stats
             .refactor_supernode_algorithm5_payoff_group_max_target_slots,
           stats
             .refactor_supernode_algorithm5_payoff_group_max_run_target_slots,
           stats.refactor_supernode_algorithm5_payoff_group_pattern_width,
           stats
             .refactor_supernode_algorithm5_payoff_group_max_pattern_width,
           stats.refactor_supernode_algorithm5_payoff_current_state_rows,
           stats
             .refactor_supernode_algorithm5_payoff_current_state_max_rows,
           stats
             .refactor_supernode_algorithm5_payoff_current_state_span_rows,
           stats
             .refactor_supernode_algorithm5_payoff_current_state_max_span_rows,
           stats.refactor_supernode_algorithm5_payoff_runtime_workspace_rows,
           stats.refactor_supernode_algorithm5_payoff_runtime_target_slots,
           stats.refactor_supernode_algorithm5_payoff_runtime_current_count,
           stats.refactor_last_supernode_algorithm5_payoff_slot_accum_runs,
           stats.refactor_last_supernode_algorithm5_payoff_slot_accum_rows,
           stats
             .refactor_last_supernode_algorithm5_payoff_slot_accum_target_entries,
           stats
             .refactor_last_supernode_algorithm5_payoff_slot_accum_target_slots,
           stats.refactor_last_supernode_algorithm5_payoff_prefix_prep_runs,
           stats.refactor_last_supernode_algorithm5_payoff_prefix_prep_rows,
           stats
             .refactor_last_supernode_algorithm5_payoff_prefix_prep_target_entries,
           stats
             .refactor_last_supernode_algorithm5_payoff_prefix_prep_target_slots,
           stats
             .refactor_last_supernode_algorithm5_payoff_current_state_seed_runs,
           stats
             .refactor_last_supernode_algorithm5_payoff_current_state_seed_deps,
           stats
             .refactor_last_supernode_algorithm5_payoff_current_state_seed_rows,
           stats
             .refactor_last_supernode_algorithm5_payoff_final_trigger_batches,
           stats
             .refactor_last_supernode_algorithm5_payoff_final_trigger_multi_batches,
           stats
             .refactor_last_supernode_algorithm5_payoff_final_trigger_claims,
           stats
             .refactor_last_supernode_algorithm5_payoff_final_trigger_suffix_deps,
           stats.refactor_last_supernode_algorithm5_payoff_suffix_advance_slots,
           stats.refactor_last_supernode_algorithm5_payoff_suffix_advance_deps,
           stats
             .refactor_last_supernode_algorithm5_payoff_suffix_advance_updates,
           stats
             .refactor_last_supernode_algorithm5_payoff_suffix_advance_finished,
           stats.refactor_supernode_consumer_plan_group_l_built,
           stats.refactor_supernode_consumer_plan_group_l_storage_limited,
           stats.refactor_supernode_consumer_plan_group_l_panel_count,
           stats.refactor_supernode_consumer_plan_group_l_run_count,
           stats.refactor_supernode_consumer_plan_group_l_rows,
           stats.refactor_supernode_consumer_plan_group_l_run_rows,
           stats.refactor_supernode_consumer_plan_group_l_exec_run_count,
           stats.refactor_supernode_consumer_plan_group_l_exec_run_rows,
           stats.refactor_supernode_consumer_plan_group_l_dense_entries,
           stats.refactor_supernode_consumer_plan_group_l_trailing_entries,
           stats.refactor_supernode_consumer_plan_group_l_bytes,
           stats.refactor_supernode_consumer_plan_group_l_state_built,
           stats.refactor_supernode_consumer_plan_group_l_state_storage_limited,
           stats.refactor_supernode_consumer_plan_group_l_state_focus_enabled,
           stats.refactor_supernode_consumer_plan_group_l_state_run_count,
           stats.refactor_supernode_consumer_plan_group_l_state_rows,
           stats.refactor_supernode_consumer_plan_group_l_state_max_rows,
           stats.refactor_supernode_consumer_plan_group_l_state_bytes,
           stats
             .refactor_supernode_consumer_plan_group_l_state_candidate_group_count,
           stats
             .refactor_supernode_consumer_plan_group_l_state_candidate_run_count,
           stats.refactor_supernode_consumer_plan_group_l_state_candidate_rows,
           stats
             .refactor_supernode_consumer_plan_group_l_state_selected_group_count,
           stats
             .refactor_supernode_consumer_plan_group_l_state_selected_run_count,
           stats.refactor_supernode_consumer_plan_group_l_state_selected_rows,
           stats.refactor_supernode_consumer_plan_group_l_state_selected_bytes,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_count,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_run_count,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_run_rows,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_entries,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_max_runs,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_advance_count,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_advance_run_count,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_advance_run_rows,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_advance_work,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_advance_max_work,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_count,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_run_count,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_run_rows,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_entries,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_advance_work,
           stats.refactor_supernode_consumer_plan_group_l_batch_candidate_max_payoff_ratio,
           stats.refactor_last_supernode_consumer_plan_group_l_dense_writes,
           stats.refactor_last_supernode_consumer_plan_group_l_trailing_writes,
           stats.refactor_last_supernode_consumer_plan_group_l_invalidations,
           stats.refactor_supernode_consumer_plan_group_l_dense_write_count,
           stats.refactor_supernode_consumer_plan_group_l_trailing_write_count,
           stats.refactor_supernode_consumer_plan_group_l_invalidation_count,
           stats.refactor_last_supernode_consumer_plan_group_l_update_runs,
           stats.refactor_last_supernode_consumer_plan_group_l_update_rows,
           stats.refactor_last_supernode_consumer_plan_group_l_update_entries,
           stats.refactor_supernode_consumer_plan_group_l_update_run_count,
           stats.refactor_supernode_consumer_plan_group_l_update_rows,
           stats.refactor_supernode_consumer_plan_group_l_update_entries);
    printf("refactor supernode consumer plan exec: attempts %" PRId64
           ", hits %" PRId64 ", applied %" PRId64
           ", rows %" PRId64 ", entries %" PRId64
           ", cumulative %" PRId64 "/%" PRId64 "/%" PRId64
           "/%" PRId64 "/%" PRId64
           ", disabled %d/%" PRId64
           ", claims %" PRId64 "/%" PRId64 "/%" PRId64
           ", cumulative claims %" PRId64 "/%" PRId64 "/%" PRId64 "\n",
           stats.refactor_last_supernode_consumer_plan_attempts,
           stats.refactor_last_supernode_consumer_plan_hits,
           stats.refactor_last_supernode_consumer_plan_applied,
           stats.refactor_last_supernode_consumer_plan_rows,
           stats.refactor_last_supernode_consumer_plan_entries,
           stats.refactor_supernode_consumer_plan_attempt_count,
           stats.refactor_supernode_consumer_plan_hit_count,
           stats.refactor_supernode_consumer_plan_apply_count,
           stats.refactor_supernode_consumer_plan_apply_rows,
           stats.refactor_supernode_consumer_plan_apply_entries,
           stats.refactor_supernode_consumer_plan_exec_disabled,
           stats.refactor_supernode_consumer_plan_exec_disable_count,
           stats.refactor_last_supernode_consumer_plan_claimed_columns,
           stats.refactor_last_supernode_consumer_plan_claim_skips,
           stats.refactor_last_supernode_consumer_plan_claim_waits,
           stats.refactor_supernode_consumer_plan_claimed_columns,
           stats.refactor_supernode_consumer_plan_claim_skip_count,
           stats.refactor_supernode_consumer_plan_claim_wait_count);
    printf("refactor cached probe rejects: shape %" PRId64 "/%" PRId64
           ", stream %" PRId64 "/%" PRId64
           ", work %" PRId64 "/%" PRId64
           ", workspace %" PRId64 "/%" PRId64
           ", cumulative %" PRId64 "/%" PRId64 "/%" PRId64 "/%" PRId64
           "\n",
           stats.refactor_last_supernode_cached_probe_shape_rejects,
           stats.refactor_last_supernode_cached_probe_shape_reject_rows,
           stats.refactor_last_supernode_cached_probe_stream_rejects,
           stats.refactor_last_supernode_cached_probe_stream_reject_rows,
           stats.refactor_last_supernode_cached_probe_work_rejects,
           stats.refactor_last_supernode_cached_probe_work_reject_rows,
           stats.refactor_last_supernode_cached_probe_workspace_rejects,
           stats.refactor_last_supernode_cached_probe_workspace_reject_rows,
           stats.refactor_supernode_cached_probe_shape_rejects,
           stats.refactor_supernode_cached_probe_stream_rejects,
           stats.refactor_supernode_cached_probe_work_rejects,
           stats.refactor_supernode_cached_probe_workspace_rejects);
    printf("refactor supernode update disabled: %d, disable count %" PRId64
           "\n",
           stats.refactor_supernode_update_disabled,
           stats.refactor_supernode_update_disable_count);
    printf("refactor L pattern: columns %" PRId64
           ", entries %" PRId64 ", adjacent runs %" PRId64
           ", adjacent entries %" PRId64 ", max run %" PRId64
           ", suffix columns %" PRId64 ", suffix entries %" PRId64
           ", max suffix %" PRId64
           ", map index32 %d/%" PRId64
           ", L index32 %d/%" PRId64
           ", U index32 %d/%" PRId64 "\n",
           stats.refactor_l_pattern_columns,
           stats.refactor_l_pattern_entries,
           stats.refactor_l_adjacent_run_count,
           stats.refactor_l_adjacent_run_entries,
           stats.refactor_l_adjacent_run_max_len,
           stats.refactor_l_contiguous_suffix_columns,
           stats.refactor_l_contiguous_suffix_entries,
           stats.refactor_l_contiguous_suffix_max_len,
           stats.refactor_map_index32_enabled,
           stats.refactor_map_index32_entries,
           stats.refactor_l_index32_enabled,
           stats.refactor_l_index32_entries,
           stats.refactor_u_index32_enabled,
           stats.refactor_u_index32_entries);
    printf("refactor supernode pipeline tasks: %" PRId64
           ", columns: %" PRId64 ", cumulative tasks: %" PRId64
           ", cumulative columns: %" PRId64 "\n",
           stats.refactor_last_supernode_pipeline_tasks,
           stats.refactor_last_supernode_pipeline_columns,
           stats.refactor_supernode_pipeline_task_count,
           stats.refactor_supernode_pipeline_column_count);
    printf("refactor supernode numeric updates: %" PRId64
           ", rows: %" PRId64 ", entries: %" PRId64
           ", cumulative updates: %" PRId64
           ", rows: %" PRId64 ", entries: %" PRId64 "\n",
           stats.refactor_last_supernode_update_runs,
           stats.refactor_last_supernode_update_rows,
           stats.refactor_last_supernode_update_entries,
           stats.refactor_supernode_update_run_count,
           stats.refactor_supernode_update_rows,
           stats.refactor_supernode_update_entries);
    printf("refactor supernode CBLAS updates: %" PRId64
           ", rows: %" PRId64 ", entries: %" PRId64
           ", cumulative updates: %" PRId64
           ", rows: %" PRId64 ", entries: %" PRId64 "\n",
           stats.refactor_last_supernode_cblas_update_runs,
           stats.refactor_last_supernode_cblas_update_rows,
           stats.refactor_last_supernode_cblas_update_entries,
           stats.refactor_supernode_cblas_update_run_count,
           stats.refactor_supernode_cblas_update_rows,
           stats.refactor_supernode_cblas_update_entries);
    printf("refactor supernode blocked updates: %" PRId64
           ", rows: %" PRId64 ", entries: %" PRId64
           ", cumulative updates: %" PRId64
           ", rows: %" PRId64 ", entries: %" PRId64 "\n",
           stats.refactor_last_supernode_blocked_update_runs,
           stats.refactor_last_supernode_blocked_update_rows,
           stats.refactor_last_supernode_blocked_update_entries,
           stats.refactor_supernode_blocked_update_run_count,
           stats.refactor_supernode_blocked_update_rows,
           stats.refactor_supernode_blocked_update_entries);
    printf("refactor ready queue columns: %" PRId64
           ", runs: %" PRId64 "\n",
           stats.refactor_last_ready_queue_columns,
           stats.refactor_ready_queue_run_count);
    printf("row refactor groups: %" PRId64
           " (single %" PRId64 ", batch %" PRId64 "/%" PRId64
           " rows, batch max/le4/le8 %" PRId64 "/%" PRId64 "/%" PRId64
           ", scalar candidates %" PRId64 "/%" PRId64
           " rows, short %" PRId64 "/%" PRId64
           " rows, stops level/internal/segment/max/end "
           "%" PRId64 "/%" PRId64 "/%" PRId64 "/%" PRId64 "/%" PRId64
           ", generic %" PRId64 "/%" PRId64 " rows max %" PRId64
           ", dense %" PRId64 "/%" PRId64 " rows max %" PRId64 ")"
           ", work single/batch/generic/dense %.6g/%.6g/%.6g/%.6g"
           ", levels: %" PRId64 ", max level width: %" PRId64
           ", cluster levels: %" PRId64 ", pipeline groups: %" PRId64
           ", pipeline rows: %" PRId64 ", pipeline work: %.6g"
           ", total work: %.6g"
           ", auto enabled: %d, auto values ready: %d"
           ", auto work allowed: %d, auto should run: %d"
           ", auto lower bound work/rejected: %.6g/%d"
           ", auto pattern/value failures: %d/%d"
           ", auto model rec/attempt/accept: %d/%d/%d"
           ", edges: %" PRId64 ", roots: %" PRId64
           ", leaves: %" PRId64 ", max fanout: %" PRId64
           ", last run: %d, last checked: %d, last parallel: %d"
           ", last ready queue: %d"
           ", runs: %" PRId64 ", checked runs: %" PRId64
           ", parallel runs: %" PRId64
           ", ready queue runs: %" PRId64
           ", ready queue groups: %" PRId64
           ", last done bitmap: %d"
           ", done bitmap runs: %" PRId64
           ", prefactor: %d rows/deps %" PRId64 "/%" PRId64
           ", prefactor totals: %" PRId64 " runs, rows/deps %" PRId64
           "/%" PRId64
           ", prefactor supernode: %d rows/deps %" PRId64 "/%" PRId64
           ", prefactor supernode totals: %" PRId64
           " runs, rows/deps %" PRId64 "/%" PRId64
           ", input cleanup rows: %" PRId64
           ", input cleanup entries: %" PRId64
           ", last defer scatter: %d"
           ", defer scatter runs: %" PRId64
           ", values dirty: %d"
           ", last lazy scatter: %d"
           ", lazy scatter runs: %" PRId64
           ", last work queue: %d"
           ", work queue runs: %" PRId64
           ", queue workspace groups: %" PRId64
           ", partial supernode pipeline: %d"
           ", partial supernode groups/rows: %" PRId64 "/%" PRId64
           ", partial supernode runs: %" PRId64
           ", local ready groups: %" PRId64 "/%" PRId64
           ", private ready groups: %" PRId64 "/%" PRId64
           ", separator-private queue: %d, components %" PRId64 "/%" PRId64
           "\n",
           stats.row_refactor_group_count,
           stats.row_refactor_group_single_count,
           stats.row_refactor_group_batch_count,
           stats.row_refactor_group_batch_rows,
           stats.row_refactor_group_batch_max_width,
           stats.row_refactor_group_batch_width_le_4_count,
           stats.row_refactor_group_batch_width_le_8_count,
           stats.row_refactor_group_scalar_candidate_count,
           stats.row_refactor_group_scalar_candidate_rows,
           stats.row_refactor_group_scalar_short_count,
           stats.row_refactor_group_scalar_short_rows,
           stats.row_refactor_group_scalar_stop_level_mismatch_count,
           stats.row_refactor_group_scalar_stop_internal_dep_count,
           stats.row_refactor_group_scalar_stop_next_segment_count,
           stats.row_refactor_group_scalar_stop_max_width_count,
           stats.row_refactor_group_scalar_stop_matrix_end_count,
           stats.row_refactor_group_generic_count,
           stats.row_refactor_group_generic_rows,
           stats.row_refactor_group_generic_max_width,
           stats.row_refactor_group_dense_count,
           stats.row_refactor_group_dense_rows,
           stats.row_refactor_group_dense_max_width,
           stats.row_refactor_group_single_work,
           stats.row_refactor_group_batch_work,
           stats.row_refactor_group_generic_work,
           stats.row_refactor_group_dense_work,
           stats.row_refactor_group_level_count,
           stats.row_refactor_group_level_max_width,
           stats.row_refactor_group_cluster_levels,
           stats.row_refactor_group_pipeline_groups,
           stats.row_refactor_group_pipeline_rows,
           stats.row_refactor_group_pipeline_work,
           stats.row_refactor_total_group_work,
           stats.row_refactor_auto_enabled,
           stats.row_refactor_auto_values_ready,
           stats.row_refactor_auto_work_allowed,
           stats.row_refactor_auto_should_run,
           stats.row_refactor_auto_lower_bound_work,
           stats.row_refactor_auto_lower_bound_rejected,
           stats.row_refactor_auto_pattern_build_failed,
           stats.row_refactor_auto_value_copy_failed,
           stats.row_refactor_auto_model_recommended,
           stats.row_refactor_auto_model_attempted,
           stats.row_refactor_auto_model_accepted,
           stats.row_refactor_group_dependency_edges,
           stats.row_refactor_group_root_count,
           stats.row_refactor_group_leaf_count,
           stats.row_refactor_group_max_fanout,
           stats.row_refactor_last_run,
           stats.row_refactor_last_checked,
           stats.row_refactor_last_parallel,
           stats.row_refactor_last_ready_queue,
           stats.row_refactor_run_count,
           stats.row_refactor_checked_run_count,
           stats.row_refactor_parallel_run_count,
           stats.row_refactor_ready_queue_run_count,
           stats.row_refactor_ready_queue_group_count,
           stats.row_refactor_last_done_bitmap,
           stats.row_refactor_done_bitmap_run_count,
           stats.row_refactor_last_prefactor,
           stats.row_refactor_last_prefactor_rows,
           stats.row_refactor_last_prefactor_deps,
           stats.row_refactor_prefactor_run_count,
           stats.row_refactor_prefactor_rows,
           stats.row_refactor_prefactor_deps,
           stats.row_refactor_last_prefactor_supernode,
           stats.row_refactor_last_prefactor_supernode_rows,
           stats.row_refactor_last_prefactor_supernode_deps,
           stats.row_refactor_prefactor_supernode_run_count,
           stats.row_refactor_prefactor_supernode_rows,
           stats.row_refactor_prefactor_supernode_deps,
           stats.row_refactor_input_cleanup_rows,
           stats.row_refactor_input_cleanup_entries,
           stats.row_refactor_last_defer_value_scatter,
           stats.row_refactor_defer_value_scatter_run_count,
           stats.row_refactor_values_dirty,
           stats.row_refactor_last_lazy_value_scatter,
           stats.row_refactor_lazy_value_scatter_run_count,
           stats.row_refactor_last_work_ready_queue,
           stats.row_refactor_work_ready_queue_run_count,
           stats.row_refactor_ready_queue_workspace_groups,
           stats.row_refactor_last_partial_supernode_pipeline,
           stats.row_refactor_last_partial_supernode_pipeline_groups,
           stats.row_refactor_last_partial_supernode_pipeline_rows,
           stats.row_refactor_partial_supernode_pipeline_run_count,
           stats.row_refactor_last_local_ready_groups,
           stats.row_refactor_local_ready_group_count,
           stats.row_refactor_last_private_ready_groups,
           stats.row_refactor_private_ready_group_count,
           stats.row_refactor_last_separator_private_queue,
           stats.row_refactor_last_separator_private_components,
           stats.row_refactor_separator_private_component_count);
    printf("row refactor separator FLOP queue: %d"
           ", ordered private %d/%" PRId64
           ", components %" PRId64 "/%" PRId64
           ", private groups %" PRId64 "/%" PRId64
           ", pipeline groups %" PRId64 "/%" PRId64
           ", closure groups %" PRId64 "/%" PRId64
           ", private threads %" PRId64
           ", private groups min/max %" PRId64 "/%" PRId64
           ", private work min/max %.6g/%.6g\n",
           stats.row_refactor_last_separator_flop_queue,
           stats.row_refactor_last_separator_flop_ordered_private,
           stats.row_refactor_separator_flop_ordered_private_run_count,
           stats.row_refactor_last_separator_flop_components,
           stats.row_refactor_separator_flop_component_count,
           stats.row_refactor_last_separator_flop_private_groups,
           stats.row_refactor_separator_flop_private_group_count,
           stats.row_refactor_last_separator_flop_pipeline_groups,
           stats.row_refactor_separator_flop_pipeline_group_count,
           stats.row_refactor_last_separator_flop_closure_groups,
           stats.row_refactor_separator_flop_closure_group_count,
           stats.row_refactor_last_separator_flop_private_threads,
           stats.row_refactor_last_separator_flop_private_min_groups,
           stats.row_refactor_last_separator_flop_private_max_groups,
           stats.row_refactor_last_separator_flop_private_min_work,
           stats.row_refactor_last_separator_flop_private_max_work);
    printf("row solve parallel: runs %" PRId64
           ", L slice runs %" PRId64 ", U slice runs %" PRId64
           ", L sparse level runs %" PRId64
           ", U sparse level runs %" PRId64
           ", threads %" PRId64
           ", max thread rect entries L/U %" PRId64 "/%" PRId64 "\n",
           stats.row_solve_parallel_run_count,
           stats.row_solve_parallel_l_slice_runs,
           stats.row_solve_parallel_u_slice_runs,
           stats.row_solve_parallel_l_sparse_level_runs,
           stats.row_solve_parallel_u_sparse_level_runs,
           stats.row_solve_thread_count,
           stats.row_solve_l_thread_max_rect_entries,
           stats.row_solve_u_thread_max_rect_entries);
    printf("row solve partition: ready %d, slices %" PRId64
           ", L sparse levels %" PRId64 ", cluster levels %" PRId64
           ", max width %" PRId64
           ", L dense tail start %" PRId64 ", rows %" PRId64
           ", entries %" PRId64 ", max slice entries %" PRId64
           ", segmented rows %" PRId64 ", rect/tri entries %" PRId64 "/%" PRId64
           ", U sparse levels %" PRId64 ", cluster levels %" PRId64
           ", max width %" PRId64
           ", U dense tail start %" PRId64 ", rows %" PRId64
           ", entries %" PRId64 ", max slice entries %" PRId64
           ", segmented rows %" PRId64 ", rect/tri entries %" PRId64 "/%" PRId64
           "\n",
           stats.row_solve_partition_ready,
           stats.row_solve_partition_slices,
           stats.row_solve_l_sparse_level_count,
           stats.row_solve_l_sparse_cluster_levels,
           stats.row_solve_l_sparse_level_max_width,
           stats.row_solve_l_dense_tail_start,
           stats.row_solve_l_dense_tail_rows,
           stats.row_solve_l_dense_tail_entries,
           stats.row_solve_l_slice_max_entries,
           stats.row_solve_l_segmented_rows,
           stats.row_solve_l_rect_entries,
           stats.row_solve_l_tri_entries,
           stats.row_solve_u_sparse_level_count,
           stats.row_solve_u_sparse_cluster_levels,
           stats.row_solve_u_sparse_level_max_width,
           stats.row_solve_u_dense_tail_start,
           stats.row_solve_u_dense_tail_rows,
           stats.row_solve_u_dense_tail_entries,
           stats.row_solve_u_slice_max_entries,
           stats.row_solve_u_segmented_rows,
           stats.row_solve_u_rect_entries,
           stats.row_solve_u_tri_entries);
    printf("row refactor segments: %" PRId64
           ", rows: %" PRId64 ", max width: %" PRId64
           ", dense entries: %.6g, trailing entries: %.6g\n",
           stats.row_refactor_segment_count,
           stats.row_refactor_segment_rows,
           stats.row_refactor_segment_max_width,
           stats.row_refactor_segment_dense_entries,
           stats.row_refactor_segment_trailing_entries);
    printf("row refactor dense segments: %" PRId64
           ", rows: %" PRId64 ", max width: %" PRId64
           ", dense entries: %.6g, trailing entries: %.6g"
           ", compact panel: %d/%" PRId64
           ", eligible: %" PRId64 "/%" PRId64
           ", work: %.6g, entries: %.6g"
           ", producer runs: %" PRId64 "/%" PRId64 "/%" PRId64
           " max %" PRId64
           ", full suffix: %" PRId64 "/%" PRId64
           ", multi/fragmented rows: %" PRId64 "/%" PRId64
           ", persistent: %" PRId64 " groups/%" PRId64
           " entries, used: %d/%" PRId64
           ", blocked: %d/%" PRId64
           " rows/entries %" PRId64 "/%" PRId64
           ", native panel: %d/%d/%" PRId64
           " rows/entries %" PRId64 "/%" PRId64
           ", auto disabled: %d/%" PRId64
           ", native blocked: %" PRId64
           " rows/entries %" PRId64 "/%" PRId64
           ", native fallback/reject: %" PRId64 "/%" PRId64
           ", direct input rows: %" PRId64 "/%" PRId64
           ", compact solve values: %" PRId64 "/%" PRId64
           ", compact group solves: %" PRId64 "/%" PRId64
           " rows, %" PRId64 "/%" PRId64 " entries"
           ", compact scalar updates: %" PRId64 "/%" PRId64
           " rows, %" PRId64 "/%" PRId64 " entries"
           ", native direct input rows: %" PRId64 "/%" PRId64
           ", sparse direct input rows: %" PRId64 "/%" PRId64
           ", batch direct input rows: %" PRId64 "/%" PRId64
           ", sparse input targets: %" PRId64 "/%" PRId64
           ", target direct rows: %" PRId64 "/%" PRId64
           ", target cleanup plan: %" PRId64 "/%" PRId64
           ", target cleanup rows: %" PRId64 "/%" PRId64
           ", target cleanup entries: %" PRId64 "/%" PRId64
           ", supernode updates: %d/%" PRId64
           " rows/entries %" PRId64 "/%" PRId64
           ", partial: %d/%" PRId64
           " rows/entries %" PRId64 "/%" PRId64
           ", gemv: %d/%" PRId64
           " rows/entries %" PRId64 "/%" PRId64
           ", trsv: %d/%" PRId64
           " rows/entries %" PRId64 "/%" PRId64
           ", batch: %d/%" PRId64
           " rows/dep_rows/entries %" PRId64 "/%" PRId64 "/%" PRId64
           ", batch patterns: %" PRId64 "/%" PRId64
           ", candidates: %" PRId64 "/%" PRId64 "/%" PRId64
           ", rejected work: %" PRId64 "\n",
           stats.row_refactor_dense_segment_count,
           stats.row_refactor_dense_segment_rows,
           stats.row_refactor_dense_segment_max_width,
           stats.row_refactor_dense_segment_dense_entries,
           stats.row_refactor_dense_segment_trailing_entries,
           stats.row_refactor_last_compact_dense_panel,
           stats.row_refactor_compact_dense_panel_count,
           stats.row_refactor_compact_dense_panel_eligible_count,
           stats.row_refactor_compact_dense_panel_eligible_rows,
           stats.row_refactor_compact_dense_panel_update_work,
           stats.row_refactor_compact_dense_panel_entries,
           stats.row_refactor_dense_producer_run_count,
           stats.row_refactor_dense_producer_run_rows,
           stats.row_refactor_dense_producer_run_dep_rows,
           stats.row_refactor_dense_producer_run_max_per_row,
           stats.row_refactor_dense_producer_full_suffix_run_count,
           stats.row_refactor_dense_producer_full_suffix_rows,
           stats.row_refactor_dense_producer_multi_run_rows,
           stats.row_refactor_dense_producer_fragmented_rows,
           stats.row_refactor_compact_dense_panel_persistent_groups,
           stats.row_refactor_compact_dense_panel_persistent_entries,
           stats.row_refactor_last_compact_dense_panel_persistent,
           stats.row_refactor_compact_dense_panel_persistent_run_count,
           stats.row_refactor_last_compact_dense_panel_blocked,
           stats.row_refactor_compact_dense_panel_blocked_run_count,
           stats.row_refactor_compact_dense_panel_blocked_rows,
           stats.row_refactor_compact_dense_panel_blocked_entries,
           stats.row_refactor_native_row_panel_enabled,
           stats.row_refactor_last_native_row_panel,
           stats.row_refactor_native_row_panel_count,
           stats.row_refactor_native_row_panel_rows,
           stats.row_refactor_native_row_panel_entries,
           stats.row_refactor_native_row_panel_auto_disabled,
           stats.row_refactor_native_row_panel_auto_disable_count,
           stats.row_refactor_native_row_panel_blocked_count,
           stats.row_refactor_native_row_panel_blocked_rows,
           stats.row_refactor_native_row_panel_blocked_entries,
           stats.row_refactor_native_row_panel_fallback_count,
           stats.row_refactor_native_row_panel_checked_reject_count,
           stats.row_refactor_last_compact_dense_panel_direct_input_rows,
           stats.row_refactor_compact_dense_panel_direct_input_rows,
           stats.row_refactor_last_compact_panel_solve_values,
           stats.row_refactor_compact_panel_solve_values,
           stats.row_refactor_last_compact_panel_group_solve_rows,
           stats.row_refactor_compact_panel_group_solve_rows,
           stats.row_refactor_last_compact_panel_group_solve_entries,
           stats.row_refactor_compact_panel_group_solve_entries,
           stats.row_refactor_last_compact_panel_scalar_update_rows,
           stats.row_refactor_compact_panel_scalar_update_rows,
           stats.row_refactor_last_compact_panel_scalar_update_entries,
           stats.row_refactor_compact_panel_scalar_update_entries,
           stats.row_refactor_last_dense_segment_direct_input_rows,
           stats.row_refactor_dense_segment_direct_input_rows,
           stats.row_refactor_last_sparse_segment_direct_input_rows,
           stats.row_refactor_sparse_segment_direct_input_rows,
           stats.row_refactor_last_batch_direct_input_rows,
           stats.row_refactor_batch_direct_input_rows,
           stats.row_refactor_segment_input_target_rows,
           stats.row_refactor_segment_input_target_entries,
           stats.row_refactor_last_segment_target_input_rows,
           stats.row_refactor_segment_target_input_rows,
           stats.row_refactor_segment_input_cleanup_rows,
           stats.row_refactor_segment_input_cleanup_entries,
           stats.row_refactor_last_segment_target_cleanup_rows,
           stats.row_refactor_segment_target_cleanup_rows,
           stats.row_refactor_last_segment_target_cleanup_entries,
           stats.row_refactor_segment_target_cleanup_entries,
           stats.row_refactor_last_compact_supernode_update,
           stats.row_refactor_compact_supernode_update_count,
           stats.row_refactor_compact_supernode_update_rows,
           stats.row_refactor_compact_supernode_update_entries,
           stats.row_refactor_last_compact_supernode_partial_update,
           stats.row_refactor_compact_supernode_partial_update_count,
           stats.row_refactor_compact_supernode_partial_update_rows,
           stats.row_refactor_compact_supernode_partial_update_entries,
           stats.row_refactor_last_compact_supernode_gemv,
           stats.row_refactor_compact_supernode_gemv_count,
           stats.row_refactor_compact_supernode_gemv_rows,
           stats.row_refactor_compact_supernode_gemv_entries,
           stats.row_refactor_last_compact_supernode_trsv,
           stats.row_refactor_compact_supernode_trsv_count,
           stats.row_refactor_compact_supernode_trsv_rows,
           stats.row_refactor_compact_supernode_trsv_entries,
           stats.row_refactor_last_compact_supernode_batch,
           stats.row_refactor_compact_supernode_batch_count,
           stats.row_refactor_compact_supernode_batch_rows,
           stats.row_refactor_compact_supernode_batch_dep_rows,
           stats.row_refactor_compact_supernode_batch_entries,
           stats.row_refactor_compact_supernode_batch_pattern_count,
           stats.row_refactor_compact_supernode_batch_pattern_rows,
           stats.row_refactor_compact_supernode_batch_candidate_count,
           stats.row_refactor_compact_supernode_batch_candidate_rows,
           stats.row_refactor_compact_supernode_batch_candidate_dep_rows,
           stats.row_refactor_compact_supernode_batch_rejected_work_count);
    printf("refactor dependency cluster levels: %" PRId64
           ", pipeline columns: %" PRId64 "\n",
           stats.refactor_dependency_cluster_levels,
           stats.refactor_dependency_pipeline_columns);
    printf("refactor dependency work: %.6g, pipeline work: %.6g\n",
           stats.refactor_dependency_work,
           stats.refactor_dependency_pipeline_work);
    printf("refactor stream split: dependency %.6g, pivot %.6g, output %.6g\n",
           stats.refactor_stream_dependency_entries,
           stats.refactor_stream_pivot_entries,
           stats.refactor_stream_output_entries);
    printf("nnz(L): %" PRId64 ", nnz(U): %" PRId64 "\n", stats.nnz_l, stats.nnz_u);
    printf("estimated flops: %.6e, factor flops: %.6e\n", stats.estimated_flops, stats.factor_flops);
    printf("parallel model R1/R2: %.6g / %.6g, recommends parallel: %s\n",
           stats.parallel_model_r1,
           stats.parallel_model_r2,
           stats.parallel_model_recommends_parallel ? "yes" : "no");
    printf("NICSLU task-flow model: threads %" PRId64
           ", deps %" PRId64 ", work %.6g, finish %.6g"
           ", speedup %.6g, recommends parallel: %s\n",
           stats.parallel_task_flow_threads,
           stats.parallel_task_flow_dependencies,
           stats.parallel_task_flow_work,
           stats.parallel_task_flow_finish_time,
           stats.parallel_task_flow_speedup,
           stats.parallel_task_flow_recommends_parallel ? "yes" : "no");
    printf("rcond: %.6e, rgrowth: %.6e\n", stats.rcond, stats.rgrowth);
    printf("memory: %zu bytes, peak: %zu bytes\n", stats.memory_bytes, stats.memory_peak_bytes);
  }

  kls_destroy(solver);
  bench_index_view_free(&input_index);
  matrix_free(&a);
  free(stressed_values);
  free(x_true);
  free(b);
  free(x);
  return EXIT_SUCCESS;
}
