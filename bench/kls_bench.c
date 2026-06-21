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

static void matrix_free(matrix *a) {
  if (a == NULL) return;
  free(a->col_ptr);
  free(a->row_idx);
  free(a->values);
  memset(a, 0, sizeof(*a));
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

static void matvec(const matrix *a, const double *x, double *y) {
  memset(y, 0, (size_t)a->n * sizeof(double));
  for (int64_t col = 0; col < a->n; ++col) {
    const double xj = x[col];
    for (int64_t p = a->col_ptr[col]; p < a->col_ptr[col + 1]; ++p) {
      y[a->row_idx[p]] += a->values[p] * xj;
    }
  }
}

static double residual_norm(const matrix *a, const double *x, const double *b, double *relative_out) {
  double *ax = (double *)calloc((size_t)a->n, sizeof(double));
  if (ax == NULL) return INFINITY;
  matvec(a, x, ax);
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
          "Usage: %s <matrix.mtx> [--repeat N] [--refactor-repeat N] [--ordering auto|amd|colamd|natural|metis] [--orientation auto|normal|transpose] [--scale auto|-1|0|1|2] [--no-btf] [--no-fast-factor] [--json]\n",
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
  kls_options options;
  kls_default_options(&options);

  for (int i = 2; i < argc; ++i) {
    if (strcmp(argv[i], "--json") == 0) {
      json = 1;
    } else if (strcmp(argv[i], "--repeat") == 0 && i + 1 < argc) {
      repeat = atoi(argv[++i]);
    } else if (strcmp(argv[i], "--refactor-repeat") == 0 && i + 1 < argc) {
      refactor_repeat = atoi(argv[++i]);
    } else if (strcmp(argv[i], "--ordering") == 0 && i + 1 < argc) {
      options.ordering = parse_ordering(argv[++i]);
    } else if (strcmp(argv[i], "--orientation") == 0 && i + 1 < argc) {
      options.orientation = parse_orientation(argv[++i]);
    } else if (strcmp(argv[i], "--scale") == 0 && i + 1 < argc) {
      if (!parse_scale(argv[++i], &options.scale)) {
        usage(argv[0]);
        return EXIT_FAILURE;
      }
    } else if (strcmp(argv[i], "--no-btf") == 0) {
      options.use_btf = 0;
    } else if (strcmp(argv[i], "--no-fast-factor") == 0) {
      options.fast_factor = 0;
    } else {
      usage(argv[0]);
      return EXIT_FAILURE;
    }
  }
  if (repeat <= 0 || refactor_repeat < 0) {
    usage(argv[0]);
    return EXIT_FAILURE;
  }

  matrix a = {0};
  if (!read_matrix_market(path, &a)) {
    return EXIT_FAILURE;
  }

  double *x_true = (double *)malloc((size_t)a.n * sizeof(double));
  double *b = (double *)calloc((size_t)a.n, sizeof(double));
  double *x = (double *)calloc((size_t)a.n, sizeof(double));
  if (x_true == NULL || b == NULL || x == NULL) {
    matrix_free(&a);
    free(x_true);
    free(b);
    free(x);
    return EXIT_FAILURE;
  }
  for (int64_t i = 0; i < a.n; ++i) {
    x_true[i] = 1.0 + (double)(i % 17) * 0.01;
  }
  matvec(&a, x_true, b);

  kls_solver *solver = NULL;
  int status = kls_create(&solver);
  if (status == KLS_OK) {
    status = kls_analyze_csc(solver, KLS_INDEX_INT64, a.n,
                             a.col_ptr, a.row_idx, 0,
                             &options);
  }
  if (status == KLS_OK) {
    status = kls_factor(solver, a.values);
  }
  if (status != KLS_OK) {
    fprintf(stderr, "KLS setup failed: %s (%d)\n", kls_status_string(status), status);
    kls_destroy(solver);
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

  double factor_total = 0.0;
  double refactor_total = 0.0;
  double solve_total = 0.0;
  double tsolve_total = 0.0;

  for (int i = 0; i < repeat; ++i) {
    status = kls_factor(solver, a.values);
    if (status != KLS_OK) break;
    kls_get_stats(solver, &stats);
    factor_total += stats.factor_seconds;
  }
  for (int i = 0; i < refactor_repeat && status == KLS_OK; ++i) {
    status = kls_refactor(solver, a.values);
    if (status != KLS_OK) break;
    kls_get_stats(solver, &stats);
    refactor_total += stats.refactor_seconds;
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
    matrix_free(&a);
    free(x_true);
    free(b);
    free(x);
    return EXIT_FAILURE;
  }

  status = kls_solve(solver, 1, b, 0, x, 0);
  if (status != KLS_OK) {
    fprintf(stderr, "KLS final solve failed: %s (%d)\n", kls_status_string(status), status);
    kls_destroy(solver);
    matrix_free(&a);
    free(x_true);
    free(b);
    free(x);
    return EXIT_FAILURE;
  }
  kls_get_stats(solver, &stats);

  double rel_residual = 0.0;
  const double residual = residual_norm(&a, x, b, &rel_residual);
  const double factor_avg = factor_total / (double)repeat;
  const double refactor_avg = refactor_repeat > 0 ? refactor_total / (double)refactor_repeat : 0.0;
  const double solve_avg = solve_total / (double)repeat;
  const double tsolve_avg = tsolve_total / (double)repeat;

  if (json) {
    printf("{\"matrix\":\"%s\",\"n\":%" PRId64 ",\"nnz\":%" PRId64
           ",\"requested_orientation\":\"%s\",\"orientation\":\"%s\""
           ",\"ordering\":\"%s\",\"requested_scale\":\"%s\",\"scale\":%d"
           ",\"btf\":%s,\"fast_factor\":%s"
           ",\"analysis_seconds\":%.9g"
           ",\"initial_factor_seconds\":%.9g"
           ",\"factor_seconds_avg\":%.9g,\"refactor_seconds_avg\":%.9g"
           ",\"solve_seconds_avg\":%.9g,\"transpose_solve_seconds_avg\":%.9g"
           ",\"residual_l2\":%.9g,\"relative_residual_l2\":%.9g"
           ",\"nnz_l\":%" PRId64 ",\"nnz_u\":%" PRId64
           ",\"estimated_flops\":%.9g,\"factor_flops\":%.9g"
           ",\"memory_bytes\":%zu,\"memory_peak_bytes\":%zu}\n",
           path, a.n, a.nnz, kls_orientation_name(options.orientation),
           kls_orientation_name(stats.selected_orientation),
           kls_ordering_name(stats.selected_ordering),
           scale_name(options.scale), stats.selected_scale,
           options.use_btf ? "true" : "false",
           options.fast_factor ? "true" : "false",
           stats.analysis_seconds, initial_factor_seconds,
           factor_avg, refactor_avg, solve_avg, tsolve_avg,
           residual, rel_residual, stats.nnz_l, stats.nnz_u,
           stats.estimated_flops, stats.factor_flops,
           stats.memory_bytes, stats.memory_peak_bytes);
  } else {
    printf("matrix: %s\n", path);
    printf("n: %" PRId64 ", nnz: %" PRId64 "\n", a.n, a.nnz);
    printf("requested orientation: %s\n", kls_orientation_name(options.orientation));
    printf("selected orientation: %s\n", kls_orientation_name(stats.selected_orientation));
    printf("ordering: %s\n", kls_ordering_name(stats.selected_ordering));
    printf("requested scale: %s\n", scale_name(options.scale));
    printf("selected scale: %d\n", stats.selected_scale);
    printf("btf: %s\n", options.use_btf ? "on" : "off");
    printf("fast factor: %s\n", options.fast_factor ? "on" : "off");
    printf("analysis: %.6f s\n", stats.analysis_seconds);
    printf("initial factor: %.6f s\n", initial_factor_seconds);
    printf("factor avg: %.6f s\n", factor_avg);
    printf("refactor avg: %.6f s\n", refactor_avg);
    printf("solve avg: %.6f s\n", solve_avg);
    printf("transpose solve avg: %.6f s\n", tsolve_avg);
    printf("residual: %.6e, relative: %.6e\n", residual, rel_residual);
    printf("nnz(L): %" PRId64 ", nnz(U): %" PRId64 "\n", stats.nnz_l, stats.nnz_u);
    printf("estimated flops: %.6e, factor flops: %.6e\n", stats.estimated_flops, stats.factor_flops);
    printf("memory: %zu bytes, peak: %zu bytes\n", stats.memory_bytes, stats.memory_peak_bytes);
  }

  kls_destroy(solver);
  matrix_free(&a);
  free(x_true);
  free(b);
  free(x);
  return EXIT_SUCCESS;
}
