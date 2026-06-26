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
          "Usage: %s <matrix.mtx> [--repeat N] [--refactor-repeat N] [--threads N] [--ordering auto|amd|colamd|natural|metis|scotch] [--orientation auto|normal|transpose] [--scale auto|-1|0|1|2] [--pivot-tol T] [--row-refactor env|off|refactor|checked|all] [--kls-first-factor env|off|on] [--row-solve env|off|on] [--stress-diagonal-scale S] [--stress-diagonal-column C] [--no-btf] [--no-fast-factor] [--no-static-pivoting] [--analyze-only] [--json]\n",
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

  if (analyze_only) {
    kls_solver *solver = NULL;
    int status = kls_create(&solver);
    if (status == KLS_OK) {
      status = kls_analyze_csc(solver, KLS_INDEX_INT64, a.n,
                               a.col_ptr, a.row_idx, 0,
                               &options);
    }
    if (status != KLS_OK) {
      fprintf(stderr, "KLS analyze failed: %s (%d)\n",
              kls_status_string(status), status);
      kls_destroy(solver);
      matrix_free(&a);
      return EXIT_FAILURE;
    }
    kls_stats stats;
    stats.struct_size = sizeof(stats);
    kls_get_stats(solver, &stats);
    if (json) {
      printf("{\"matrix\":\"%s\",\"n\":%" PRId64 ",\"nnz\":%" PRId64
             ",\"threads\":%d"
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
             ",\"analyze_only\":true}\n",
             path, a.n, a.nnz, options.threads,
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
             stats.parallel_model_recommends_parallel);
    } else {
      printf("matrix: %s\n", path);
      printf("n: %" PRId64 ", nnz: %" PRId64 "\n", a.n, a.nnz);
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
    }
    kls_destroy(solver);
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
      matrix_free(&a);
      return EXIT_FAILURE;
    }
    run_values = stressed_values;
  }

  double *x_true = (double *)malloc((size_t)a.n * sizeof(double));
  double *b = (double *)calloc((size_t)a.n, sizeof(double));
  double *x = (double *)calloc((size_t)a.n, sizeof(double));
  if (x_true == NULL || b == NULL || x == NULL) {
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

  for (int i = 0; i < repeat; ++i) {
    status = kls_factor(solver, run_values);
    if (status != KLS_OK) break;
    kls_get_stats(solver, &stats);
    factor_total += stats.factor_seconds;
  }
  for (int i = 0; i < refactor_repeat && status == KLS_OK; ++i) {
    status = kls_refactor(solver, run_values);
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
           ",\"selected_spral_matching\":%s"
           ",\"analysis_seconds\":%.9g"
           ",\"initial_factor_seconds\":%.9g"
           ",\"initial_factor_path\":\"%s\""
           ",\"last_factor_path\":\"%s\""
           ",\"kls_tail_last_mapped_columns\":%" PRId64
           ",\"kls_tail_mapped_column_count\":%" PRId64
           ",\"kls_first_last_row_uplooking_columns\":%" PRId64
           ",\"kls_first_row_uplooking_column_count\":%" PRId64
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
           ",\"kls_first_last_separator_queue_executed\":%d"
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
           ",\"fast_block_restarts\":%d"
           ",\"fast_kls_block_restarts\":%d"
           ",\"fast_tail_restarts\":%d"
           ",\"fast_repaired_last_offdiag_suffix_refresh\":%d"
           ",\"fast_repaired_offdiag_suffix_refresh_count\":%" PRId64
           ",\"fast_repaired_offdiag_full_refresh_count\":%" PRId64
           ",\"fast_repaired_parallel_tail_blocks\":%" PRId64,
           path, a.n, a.nnz, options.threads,
           stats.build_has_metis ? "true" : "false",
           stats.build_has_scotch ? "true" : "false",
           stats.build_has_spral_scaling ? "true" : "false",
           stats.build_has_cblas ? "true" : "false",
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
           stats.selected_spral_matching ? "true" : "false",
           stats.analysis_seconds, initial_factor_seconds,
           kls_factor_path_name(initial_factor_path),
           kls_factor_path_name(stats.last_factor_path),
           stats.kls_tail_last_mapped_columns,
           stats.kls_tail_mapped_column_count,
           stats.kls_first_last_row_uplooking_columns,
           stats.kls_first_row_uplooking_column_count,
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
           stats.fast_block_restarts,
           stats.fast_kls_block_restarts,
           stats.fast_tail_restarts,
           stats.fast_repaired_last_offdiag_suffix_refresh,
           stats.fast_repaired_offdiag_suffix_refresh_count,
           stats.fast_repaired_offdiag_full_refresh_count,
           stats.fast_repaired_parallel_tail_blocks);
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
           ",\"kls_first_last_separator_queue_pipeline_pivot_tail\":%d"
           ",\"kls_first_separator_queue_pipeline_pivot_tail_run_count\":%" PRId64
           ",\"kls_first_last_separator_queue_pipeline_pivot_tail_rows\":%" PRId64
           ",\"kls_first_last_separator_queue_pipeline_pivot_restarts\":%" PRId64
           ",\"kls_first_separator_queue_pipeline_pivot_restart_count\":%" PRId64
           ",\"kls_first_last_separator_queue_pipeline_pivot_serial_rows\":%" PRId64,
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
           stats.kls_first_last_separator_queue_pipeline_pivot_tail,
           stats.kls_first_separator_queue_pipeline_pivot_tail_run_count,
           stats.kls_first_last_separator_queue_pipeline_pivot_tail_rows,
           stats.kls_first_last_separator_queue_pipeline_pivot_restarts,
           stats.kls_first_separator_queue_pipeline_pivot_restart_count,
           stats.kls_first_last_separator_queue_pipeline_pivot_serial_rows);
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
           ",\"fast_rejected_pivoting_tail_contiguous\":%d"
           ",\"fast_rejected_pivoting_tail_suffix_exact\":%d"
           ",\"fast_rejected_pivoting_tail_gap_columns\":%" PRId64
           ",\"fast_rejected_pivoting_tail_suffix_overcompute_columns\":%" PRId64
           ",\"fast_rejected_pivoting_tail_suffix_overcompute_work\":%.9g"
           ",\"fast_rejected_refresh_state\":%d"
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
           stats.fast_rejected_pivoting_tail_contiguous,
           stats.fast_rejected_pivoting_tail_suffix_exact,
           stats.fast_rejected_pivoting_tail_gap_columns,
           stats.fast_rejected_pivoting_tail_suffix_overcompute_columns,
           stats.fast_rejected_pivoting_tail_suffix_overcompute_work,
           stats.fast_rejected_refresh_state,
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
           ",\"row_refactor_group_count\":%" PRId64
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
           stats.row_refactor_group_count,
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
           stats.row_refactor_input_cleanup_rows,
           stats.row_refactor_input_cleanup_entries,
           stats.row_refactor_last_defer_value_scatter,
           stats.row_refactor_defer_value_scatter_run_count,
           stats.row_refactor_values_dirty,
           stats.row_refactor_last_lazy_value_scatter,
           stats.row_refactor_lazy_value_scatter_run_count,
           stats.row_refactor_last_row_solve,
           stats.row_refactor_row_solve_run_count);
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
           ",\"refactor_last_ready_queue_columns\":%" PRId64
           ",\"refactor_ready_queue_run_count\":%" PRId64
           ",\"nnz_l\":%" PRId64 ",\"nnz_u\":%" PRId64
           ",\"estimated_flops\":%.9g,\"factor_flops\":%.9g"
           ",\"parallel_model_r1\":%.9g"
           ",\"parallel_model_r2\":%.9g"
           ",\"parallel_model_recommends_parallel\":%d"
           ",\"rcond\":%.9g,\"rgrowth\":%.9g"
           ",\"memory_bytes\":%zu,\"memory_peak_bytes\":%zu}\n",
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
           stats.refactor_last_ready_queue_columns,
           stats.refactor_ready_queue_run_count,
           stats.nnz_l, stats.nnz_u,
           stats.estimated_flops, stats.factor_flops,
           stats.parallel_model_r1,
           stats.parallel_model_r2,
           stats.parallel_model_recommends_parallel,
           stats.rcond, stats.rgrowth,
           stats.memory_bytes, stats.memory_peak_bytes);
  } else {
    printf("matrix: %s\n", path);
    printf("n: %" PRId64 ", nnz: %" PRId64 "\n", a.n, a.nnz);
    printf("threads: %d\n", options.threads);
    printf("build features: METIS %s, SCOTCH %s, SPRAL scaling %s, CBLAS %s\n",
           stats.build_has_metis ? "on" : "off",
           stats.build_has_scotch ? "on" : "off",
           stats.build_has_spral_scaling ? "on" : "off",
           stats.build_has_cblas ? "on" : "off");
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
    printf("selected SPRAL matching: %s\n",
           stats.selected_spral_matching ? "on" : "off");
    printf("analysis: %.6f s\n", stats.analysis_seconds);
    printf("initial factor: %.6f s\n", initial_factor_seconds);
    printf("initial factor path: %s\n",
           kls_factor_path_name(initial_factor_path));
    printf("last factor path: %s\n",
           kls_factor_path_name(stats.last_factor_path));
    printf("KLS tail mapped columns: last %" PRId64
           ", total %" PRId64 "\n",
           stats.kls_tail_last_mapped_columns,
           stats.kls_tail_mapped_column_count);
    printf("KLS row-up-looking first columns: last %" PRId64
           ", total %" PRId64 "\n",
           stats.kls_first_last_row_uplooking_columns,
           stats.kls_first_row_uplooking_column_count);
    printf("KLS first row-refactor seed rows: last %" PRId64
           ", total %" PRId64 "\n",
           stats.kls_first_last_row_refactor_seeded_rows,
           stats.kls_first_row_refactor_seeded_row_count);
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
           ", threads/max_rows %" PRId64 "/%" PRId64 "\n",
           stats.kls_first_last_separator_queue,
           stats.kls_first_separator_queue_run_count,
           stats.kls_first_last_separator_queue_private_components,
           stats.kls_first_last_separator_queue_pipeline_components,
           stats.kls_first_last_separator_queue_private_rows,
           stats.kls_first_last_separator_queue_pipeline_rows,
           stats.kls_first_last_separator_queue_nonempty_threads,
           stats.kls_first_last_separator_queue_max_thread_rows);
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
    printf("KLS row-up-looking separator pipeline pivot tail: last %d"
           ", runs %" PRId64 ", rows %" PRId64
           ", restarts %" PRId64 ", serial rows %" PRId64 "\n",
           stats.kls_first_last_separator_queue_pipeline_pivot_tail,
           stats.kls_first_separator_queue_pipeline_pivot_tail_run_count,
           stats.kls_first_last_separator_queue_pipeline_pivot_tail_rows,
           stats.kls_first_last_separator_queue_pipeline_pivot_restarts,
           stats.kls_first_last_separator_queue_pipeline_pivot_serial_rows);
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
           ", exact mask %d\n",
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
           stats.fast_repaired_tail_restart_exact_mask);
    printf("fast block restarts: %d, KLS block restarts: %d"
           ", tail restarts: %d"
           ", offdiag suffix refresh last %d, suffix refreshes %" PRId64
           ", full refreshes %" PRId64
           ", parallel tail blocks %" PRId64 "\n",
           stats.fast_block_restarts, stats.fast_kls_block_restarts,
           stats.fast_tail_restarts,
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
           ", contiguous %d, suffix exact %d"
           ", gaps %" PRId64 ", suffix overcompute columns %" PRId64
           ", suffix overcompute work %.6g, refresh state %d"
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
           stats.fast_rejected_pivoting_tail_contiguous,
           stats.fast_rejected_pivoting_tail_suffix_exact,
           stats.fast_rejected_pivoting_tail_gap_columns,
           stats.fast_rejected_pivoting_tail_suffix_overcompute_columns,
           stats.fast_rejected_pivoting_tail_suffix_overcompute_work,
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
           ", dense entries: %.6g, trailing entries: %.6g\n",
           stats.refactor_supernode_candidate_count,
           stats.refactor_supernode_candidate_rows,
           stats.refactor_supernode_candidate_max_width,
           stats.refactor_supernode_candidate_dense_entries,
           stats.refactor_supernode_candidate_trailing_entries);
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
    printf("refactor ready queue columns: %" PRId64
           ", runs: %" PRId64 "\n",
           stats.refactor_last_ready_queue_columns,
           stats.refactor_ready_queue_run_count);
    printf("row refactor groups: %" PRId64
           ", levels: %" PRId64 ", max level width: %" PRId64
           ", cluster levels: %" PRId64 ", pipeline groups: %" PRId64
           ", pipeline rows: %" PRId64 ", pipeline work: %.6g"
           ", total work: %.6g"
           ", auto enabled: %d, auto values ready: %d"
           ", auto work allowed: %d, auto should run: %d"
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
           ", components %" PRId64 "/%" PRId64
           ", private groups %" PRId64 "/%" PRId64
           ", pipeline groups %" PRId64 "/%" PRId64
           ", closure groups %" PRId64 "/%" PRId64 "\n",
           stats.row_refactor_last_separator_flop_queue,
           stats.row_refactor_last_separator_flop_components,
           stats.row_refactor_separator_flop_component_count,
           stats.row_refactor_last_separator_flop_private_groups,
           stats.row_refactor_separator_flop_private_group_count,
           stats.row_refactor_last_separator_flop_pipeline_groups,
           stats.row_refactor_separator_flop_pipeline_group_count,
           stats.row_refactor_last_separator_flop_closure_groups,
           stats.row_refactor_separator_flop_closure_group_count);
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
    printf("nnz(L): %" PRId64 ", nnz(U): %" PRId64 "\n", stats.nnz_l, stats.nnz_u);
    printf("estimated flops: %.6e, factor flops: %.6e\n", stats.estimated_flops, stats.factor_flops);
    printf("parallel model R1/R2: %.6g / %.6g, recommends parallel: %s\n",
           stats.parallel_model_r1,
           stats.parallel_model_r2,
           stats.parallel_model_recommends_parallel ? "yes" : "no");
    printf("rcond: %.6e, rgrowth: %.6e\n", stats.rcond, stats.rgrowth);
    printf("memory: %zu bytes, peak: %zu bytes\n", stats.memory_bytes, stats.memory_peak_bytes);
  }

  kls_destroy(solver);
  matrix_free(&a);
  free(stressed_values);
  free(x_true);
  free(b);
  free(x);
  return EXIT_SUCCESS;
}
