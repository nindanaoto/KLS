#define _POSIX_C_SOURCE 200809L

#include "klu.h"

#include <ctype.h>
#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <time.h>

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

typedef struct run_stats {
  double analysis_seconds;
  double initial_factor_seconds;
  double factor_seconds_avg;
  double refactor_seconds_avg;
  double solve_seconds_avg;
  double spice_cycle_seconds;
  int status;
  int nblocks;
  int64_t nnz_l;
  int64_t nnz_u;
} run_stats;

static double now_seconds(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1.0e-9 * (double)ts.tv_nsec;
}

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
    if (tolower((unsigned char)*s) !=
        tolower((unsigned char)*prefix)) {
      return 0;
    }
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
    triplet *next =
      (triplet *)realloc(*items, (size_t)next_capacity * sizeof(*next));
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
    fprintf(stderr, "%s uses complex values; this tool supports real matrices only\n",
            path);
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

  int64_t rows = 0;
  int64_t cols = 0;
  int64_t entries = 0;
  if (sscanf(line, "%" SCNd64 " %" SCNd64 " %" SCNd64, &rows, &cols,
             &entries) != 3 ||
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
    int64_t row = 0;
    int64_t col = 0;
    double value = 1.0;
    const int fields = is_pattern
                         ? sscanf(line, "%" SCNd64 " %" SCNd64, &row, &col)
                         : sscanf(line, "%" SCNd64 " %" SCNd64 " %lf", &row,
                                  &col, &value);
    if (fields < 2 || row <= 0 || row > rows || col <= 0 || col > cols) {
      fprintf(stderr, "%s has an invalid entry near %" PRId64 "\n", path,
              k + 1);
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

  qsort(items, (size_t)count, sizeof(*items), cmp_triplet);
  int64_t unique = 0;
  for (int64_t i = 0; i < count;) {
    int64_t j = i + 1;
    double sum = items[i].value;
    while (j < count && items[j].row == items[i].row &&
           items[j].col == items[i].col) {
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
  out->col_ptr = (int64_t *)calloc((size_t)rows + 1u, sizeof(*out->col_ptr));
  out->row_idx = (int64_t *)calloc((size_t)unique, sizeof(*out->row_idx));
  out->values = (double *)calloc((size_t)unique, sizeof(*out->values));
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
  int64_t *next = (int64_t *)malloc((size_t)rows * sizeof(*next));
  if (next == NULL) {
    free(items);
    matrix_free(out);
    return 0;
  }
  memcpy(next, out->col_ptr, (size_t)rows * sizeof(*next));
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
  memset(y, 0, (size_t)a->n * sizeof(*y));
  for (int64_t col = 0; col < a->n; ++col) {
    const double xj = x[col];
    for (int64_t p = a->col_ptr[col]; p < a->col_ptr[col + 1]; ++p) {
      y[a->row_idx[p]] += a->values[p] * xj;
    }
  }
}

static int copy_int_arrays(const matrix *a, int **ap_out, int **ai_out) {
  if (a->n > INT_MAX || a->nnz > INT_MAX) return 0;
  int *ap = (int *)malloc(((size_t)a->n + 1u) * sizeof(*ap));
  int *ai = (int *)malloc((size_t)a->nnz * sizeof(*ai));
  if (ap == NULL || ai == NULL) {
    free(ap);
    free(ai);
    return 0;
  }
  for (int64_t i = 0; i <= a->n; ++i) ap[i] = (int)a->col_ptr[i];
  for (int64_t i = 0; i < a->nnz; ++i) ai[i] = (int)a->row_idx[i];
  *ap_out = ap;
  *ai_out = ai;
  return 1;
}

static int copy_long_arrays(const matrix *a, int64_t **ap_out,
                            int64_t **ai_out) {
  int64_t *ap = (int64_t *)malloc(((size_t)a->n + 1u) * sizeof(*ap));
  int64_t *ai = (int64_t *)malloc((size_t)a->nnz * sizeof(*ai));
  if (ap == NULL || ai == NULL) {
    free(ap);
    free(ai);
    return 0;
  }
  for (int64_t i = 0; i <= a->n; ++i) ap[i] = (int64_t)a->col_ptr[i];
  for (int64_t i = 0; i < a->nnz; ++i) ai[i] = (int64_t)a->row_idx[i];
  *ap_out = ap;
  *ai_out = ai;
  return 1;
}

static double cycle_seconds(const run_stats *s) {
  return s->analysis_seconds + s->initial_factor_seconds + s->solve_seconds_avg +
         99.0 * (s->refactor_seconds_avg + s->solve_seconds_avg);
}

static int run_klu32(const matrix *a, const double *rhs, int repeat,
                     int refactor_repeat, int ordering, int btf, int scale,
                     run_stats *out) {
  int *ap = NULL;
  int *ai = NULL;
  if (!copy_int_arrays(a, &ap, &ai)) return 0;

  klu_common common;
  klu_defaults(&common);
  common.ordering = ordering;
  common.btf = btf;
  common.scale = scale;

  const double analysis_start = now_seconds();
  klu_symbolic *symbolic =
    klu_analyze((int)a->n, ap, ai, &common);
  out->analysis_seconds = now_seconds() - analysis_start;
  if (symbolic == NULL || common.status < 0) {
    out->status = common.status;
    free(ap);
    free(ai);
    return 0;
  }
  out->nblocks = symbolic->nblocks;

  const double factor_start = now_seconds();
  klu_numeric *numeric =
    klu_factor(ap, ai, a->values, symbolic, &common);
  out->initial_factor_seconds = now_seconds() - factor_start;
  if (numeric == NULL || common.status < 0) {
    out->status = common.status;
    klu_free_symbolic(&symbolic, &common);
    free(ap);
    free(ai);
    return 0;
  }
  out->nnz_l = numeric->lnz;
  out->nnz_u = numeric->unz;

  double total = 0.0;
  for (int i = 0; i < repeat; ++i) {
    klu_free_numeric(&numeric, &common);
    const double start = now_seconds();
    numeric = klu_factor(ap, ai, a->values, symbolic, &common);
    total += now_seconds() - start;
    if (numeric == NULL || common.status < 0) break;
  }
  out->factor_seconds_avg = total / (double)repeat;
  if (numeric == NULL || common.status < 0) {
    out->status = common.status;
    klu_free_symbolic(&symbolic, &common);
    free(ap);
    free(ai);
    return 0;
  }

  total = 0.0;
  for (int i = 0; i < refactor_repeat; ++i) {
    const double start = now_seconds();
    const int ok = klu_refactor(ap, ai, a->values, symbolic, numeric,
                                         &common);
    total += now_seconds() - start;
    if (!ok || common.status < 0) break;
  }
  out->refactor_seconds_avg =
    refactor_repeat > 0 ? total / (double)refactor_repeat : 0.0;
  if (common.status < 0) {
    out->status = common.status;
    klu_free_numeric(&numeric, &common);
    klu_free_symbolic(&symbolic, &common);
    free(ap);
    free(ai);
    return 0;
  }

  double *work = (double *)malloc((size_t)a->n * sizeof(*work));
  if (work == NULL) {
    out->status = KLU_OUT_OF_MEMORY;
    klu_free_numeric(&numeric, &common);
    klu_free_symbolic(&symbolic, &common);
    free(ap);
    free(ai);
    return 0;
  }
  total = 0.0;
  for (int i = 0; i < repeat; ++i) {
    memcpy(work, rhs, (size_t)a->n * sizeof(*work));
    const double start = now_seconds();
    const int ok = klu_solve(symbolic, numeric, (int)a->n, 1, work,
                                      &common);
    total += now_seconds() - start;
    if (!ok || common.status < 0) break;
  }
  out->solve_seconds_avg = total / (double)repeat;
  out->spice_cycle_seconds = cycle_seconds(out);
  out->status = common.status;

  free(work);
  klu_free_numeric(&numeric, &common);
  klu_free_symbolic(&symbolic, &common);
  free(ap);
  free(ai);
  return out->status >= 0;
}

static int run_klu64(const matrix *a, const double *rhs, int repeat,
                     int refactor_repeat, int ordering, int btf, int scale,
                     run_stats *out) {
  int64_t *ap = NULL;
  int64_t *ai = NULL;
  if (!copy_long_arrays(a, &ap, &ai)) return 0;

  klu_l_common common;
  klu_l_defaults(&common);
  common.ordering = (int64_t)ordering;
  common.btf = (int64_t)btf;
  common.scale = (int64_t)scale;

  const double analysis_start = now_seconds();
  klu_l_symbolic *symbolic =
    klu_l_analyze((int64_t)a->n, ap, ai, &common);
  out->analysis_seconds = now_seconds() - analysis_start;
  if (symbolic == NULL || common.status < 0) {
    out->status = (int)common.status;
    free(ap);
    free(ai);
    return 0;
  }
  out->nblocks = (int)symbolic->nblocks;

  const double factor_start = now_seconds();
  klu_l_numeric *numeric =
    klu_l_factor(ap, ai, a->values, symbolic, &common);
  out->initial_factor_seconds = now_seconds() - factor_start;
  if (numeric == NULL || common.status < 0) {
    out->status = (int)common.status;
    klu_l_free_symbolic(&symbolic, &common);
    free(ap);
    free(ai);
    return 0;
  }
  out->nnz_l = (int64_t)numeric->lnz;
  out->nnz_u = (int64_t)numeric->unz;

  double total = 0.0;
  for (int i = 0; i < repeat; ++i) {
    klu_l_free_numeric(&numeric, &common);
    const double start = now_seconds();
    numeric = klu_l_factor(ap, ai, a->values, symbolic, &common);
    total += now_seconds() - start;
    if (numeric == NULL || common.status < 0) break;
  }
  out->factor_seconds_avg = total / (double)repeat;
  if (numeric == NULL || common.status < 0) {
    out->status = (int)common.status;
    klu_l_free_symbolic(&symbolic, &common);
    free(ap);
    free(ai);
    return 0;
  }

  total = 0.0;
  for (int i = 0; i < refactor_repeat; ++i) {
    const double start = now_seconds();
    const int64_t ok =
      klu_l_refactor(ap, ai, a->values, symbolic, numeric, &common);
    total += now_seconds() - start;
    if (!ok || common.status < 0) break;
  }
  out->refactor_seconds_avg =
    refactor_repeat > 0 ? total / (double)refactor_repeat : 0.0;
  if (common.status < 0) {
    out->status = (int)common.status;
    klu_l_free_numeric(&numeric, &common);
    klu_l_free_symbolic(&symbolic, &common);
    free(ap);
    free(ai);
    return 0;
  }

  double *work = (double *)malloc((size_t)a->n * sizeof(*work));
  if (work == NULL) {
    out->status = KLU_OUT_OF_MEMORY;
    klu_l_free_numeric(&numeric, &common);
    klu_l_free_symbolic(&symbolic, &common);
    free(ap);
    free(ai);
    return 0;
  }
  total = 0.0;
  for (int i = 0; i < repeat; ++i) {
    memcpy(work, rhs, (size_t)a->n * sizeof(*work));
    const double start = now_seconds();
    const int64_t ok =
      klu_l_solve(symbolic, numeric, (int64_t)a->n, 1u, work,
                           &common);
    total += now_seconds() - start;
    if (!ok || common.status < 0) break;
  }
  out->solve_seconds_avg = total / (double)repeat;
  out->spice_cycle_seconds = cycle_seconds(out);
  out->status = (int)common.status;

  free(work);
  klu_l_free_numeric(&numeric, &common);
  klu_l_free_symbolic(&symbolic, &common);
  free(ap);
  free(ai);
  return out->status >= 0;
}

static int parse_ordering(const char *name) {
  if (strcmp(name, "amd") == 0) return 0;
  if (strcmp(name, "colamd") == 0) return 1;
  if (strcmp(name, "natural") == 0 || strcmp(name, "given") == 0) return 2;
  return -1;
}

static void print_stats_json(const char *name, const run_stats *s) {
  printf("\"%s\":{\"status\":%d,\"analysis_seconds\":%.9g,"
         "\"initial_factor_seconds\":%.9g,\"factor_seconds_avg\":%.9g,"
         "\"refactor_seconds_avg\":%.9g,\"solve_seconds_avg\":%.9g,"
         "\"spice_cycle_seconds\":%.9g,\"nblocks\":%d,"
         "\"nnz_l\":%" PRId64 ",\"nnz_u\":%" PRId64 "}",
         name, s->status, s->analysis_seconds, s->initial_factor_seconds,
         s->factor_seconds_avg, s->refactor_seconds_avg, s->solve_seconds_avg,
         s->spice_cycle_seconds, s->nblocks, s->nnz_l, s->nnz_u);
}

int main(int argc, char **argv) {
  if (argc < 2) {
    fprintf(stderr,
            "Usage: %s <matrix.mtx> [--repeat N] [--refactor-repeat N] "
            "[--ordering amd|colamd|natural] [--scale -1|0|1|2] "
            "[--no-btf] [--json]\n",
            argv[0]);
    return EXIT_FAILURE;
  }

  const char *path = argv[1];
  int repeat = 5;
  int refactor_repeat = 5;
  int ordering = 0;
  int btf = 1;
  int scale = 2;
  int json = 0;
  for (int i = 2; i < argc; ++i) {
    if (strcmp(argv[i], "--repeat") == 0 && i + 1 < argc) {
      repeat = atoi(argv[++i]);
    } else if (strcmp(argv[i], "--refactor-repeat") == 0 && i + 1 < argc) {
      refactor_repeat = atoi(argv[++i]);
    } else if (strcmp(argv[i], "--ordering") == 0 && i + 1 < argc) {
      ordering = parse_ordering(argv[++i]);
      if (ordering < 0) {
        fprintf(stderr, "unknown ordering\n");
        return EXIT_FAILURE;
      }
    } else if (strcmp(argv[i], "--scale") == 0 && i + 1 < argc) {
      scale = atoi(argv[++i]);
    } else if (strcmp(argv[i], "--no-btf") == 0) {
      btf = 0;
    } else if (strcmp(argv[i], "--json") == 0) {
      json = 1;
    } else {
      fprintf(stderr, "unknown argument: %s\n", argv[i]);
      return EXIT_FAILURE;
    }
  }
  if (repeat <= 0 || refactor_repeat < 0 ||
      (scale != -1 && scale != 0 && scale != 1 && scale != 2)) {
    fprintf(stderr, "invalid repeat or scale argument\n");
    return EXIT_FAILURE;
  }

  matrix a;
  memset(&a, 0, sizeof(a));
  if (!read_matrix_market(path, &a)) return EXIT_FAILURE;

  double *x_true = (double *)malloc((size_t)a.n * sizeof(*x_true));
  double *rhs = (double *)malloc((size_t)a.n * sizeof(*rhs));
  if (x_true == NULL || rhs == NULL) {
    matrix_free(&a);
    free(x_true);
    free(rhs);
    return EXIT_FAILURE;
  }
  for (int64_t i = 0; i < a.n; ++i) {
    x_true[i] = 1.0 + (double)((i % 17) - 8) * 0.01;
  }
  matvec(&a, x_true, rhs);

  run_stats s32;
  run_stats s64;
  memset(&s32, 0, sizeof(s32));
  memset(&s64, 0, sizeof(s64));
  const int ok32 = run_klu32(&a, rhs, repeat, refactor_repeat, ordering, btf,
                             scale, &s32);
  const int ok64 = run_klu64(&a, rhs, repeat, refactor_repeat, ordering, btf,
                             scale, &s64);

  if (json) {
    printf("{\"matrix\":\"%s\",\"n\":%" PRId64 ",\"nnz\":%" PRId64
           ",\"repeat\":%d,\"refactor_repeat\":%d,\"ordering\":%d,"
           "\"btf\":%s,\"scale\":%d,",
           path, a.n, a.nnz, repeat, refactor_repeat, ordering,
           btf ? "true" : "false", scale);
    print_stats_json("klu32", &s32);
    printf(",");
    print_stats_json("klu64", &s64);
    if (ok32 && ok64) {
      printf(",\"ratio32_over_64\":{\"analysis\":%.9g,"
             "\"initial_factor\":%.9g,\"factor\":%.9g,"
             "\"refactor\":%.9g,\"solve\":%.9g,\"spice_cycle\":%.9g}",
             s32.analysis_seconds / s64.analysis_seconds,
             s32.initial_factor_seconds / s64.initial_factor_seconds,
             s32.factor_seconds_avg / s64.factor_seconds_avg,
             s32.refactor_seconds_avg / s64.refactor_seconds_avg,
             s32.solve_seconds_avg / s64.solve_seconds_avg,
             s32.spice_cycle_seconds / s64.spice_cycle_seconds);
    }
    printf("}\n");
  } else {
    printf("KLU32: analysis %.6f factor %.6f refactor %.6f solve %.6f cycle %.6f\n",
           s32.analysis_seconds, s32.initial_factor_seconds,
           s32.refactor_seconds_avg, s32.solve_seconds_avg,
           s32.spice_cycle_seconds);
    printf("KLU64: analysis %.6f factor %.6f refactor %.6f solve %.6f cycle %.6f\n",
           s64.analysis_seconds, s64.initial_factor_seconds,
           s64.refactor_seconds_avg, s64.solve_seconds_avg,
           s64.spice_cycle_seconds);
  }

  matrix_free(&a);
  free(x_true);
  free(rhs);
  return (ok32 && ok64) ? EXIT_SUCCESS : EXIT_FAILURE;
}
