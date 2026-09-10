#define _GNU_SOURCE
#include "kls/kls.h"
#include "kls_tuning_profile.h"

#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <sched.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

typedef struct synthetic_matrix {
  int64_t n, nnz;
  int64_t *ptr, *rows;
  double *values, *base, *rhs, *x;
} synthetic_matrix;

static double now_seconds(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static uint64_t next_random(uint64_t *state) {
  uint64_t x = *state;
  x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
  *state = x;
  return x * UINT64_C(2685821657736338717);
}

/* A deterministic diagonally dominant stencil exercises sparse, supernodal,
 * BTF, and dense-tail crossovers without embedding any benchmark matrix. */
static int make_stencil(synthetic_matrix *a, int64_t n, int bandwidth,
                        uint64_t seed) {
  memset(a, 0, sizeof(*a));
  a->n = n;
  a->ptr = calloc((size_t)n + 1u, sizeof(*a->ptr));
  if (a->ptr == NULL) return 0;
  for (int64_t col = 0; col < n; ++col) {
    int64_t lo = col > bandwidth ? col - bandwidth : 0;
    int64_t hi = col + bandwidth + 1 < n ? col + bandwidth + 1 : n;
    a->ptr[col + 1] = a->ptr[col] + hi - lo;
  }
  a->nnz = a->ptr[n];
  a->rows = malloc((size_t)a->nnz * sizeof(*a->rows));
  a->values = malloc((size_t)a->nnz * sizeof(*a->values));
  a->base = malloc((size_t)a->nnz * sizeof(*a->base));
  a->rhs = malloc((size_t)n * sizeof(*a->rhs));
  a->x = malloc((size_t)n * sizeof(*a->x));
  if (a->rows == NULL || a->values == NULL || a->base == NULL ||
      a->rhs == NULL || a->x == NULL) return 0;
  for (int64_t i = 0; i < n; ++i) a->rhs[i] = 0.0;
  for (int64_t col = 0; col < n; ++col) {
    int64_t lo = col > bandwidth ? col - bandwidth : 0;
    int64_t hi = col + bandwidth + 1 < n ? col + bandwidth + 1 : n;
    double off_sum = 0.0;
    for (int64_t row = lo, p = a->ptr[col]; row < hi; ++row, ++p) {
      a->rows[p] = row;
      if (row == col) continue;
      double v = -0.01 * (1.0 + (double)(next_random(&seed) & 255u) / 1024.0);
      a->base[p] = v;
      off_sum += fabs(v);
      a->rhs[row] += v;
    }
    for (int64_t p = a->ptr[col]; p < a->ptr[col + 1]; ++p) {
      if (a->rows[p] == col) {
        a->base[p] = 1.0 + 2.0 * off_sum;
        a->rhs[col] += a->base[p];
        break;
      }
    }
  }
  memcpy(a->values, a->base, (size_t)a->nnz * sizeof(*a->values));
  return 1;
}

static void free_matrix(synthetic_matrix *a) {
  free(a->ptr); free(a->rows); free(a->values); free(a->base);
  free(a->rhs); free(a->x); memset(a, 0, sizeof(*a));
}

static double median(double *values, int count) {
  for (int i = 1; i < count; ++i) {
    double value = values[i];
    int j = i;
    while (j > 0 && values[j - 1] > value) {
      values[j] = values[j - 1]; --j;
    }
    values[j] = value;
  }
  return values[count / 2];
}

static double relative_residual(const synthetic_matrix *a) {
  double *r = malloc((size_t)a->n * sizeof(*r));
  if (r == NULL) return HUGE_VAL;
  memcpy(r, a->rhs, (size_t)a->n * sizeof(*r));
  double norm_a = 0.0, norm_x = 0.0, norm_b = 0.0, norm_r = 0.0;
  for (int64_t col = 0; col < a->n; ++col) {
    double col_sum = 0.0;
    for (int64_t p = a->ptr[col]; p < a->ptr[col + 1]; ++p) {
      r[a->rows[p]] -= a->values[p] * a->x[col];
      col_sum += fabs(a->values[p]);
    }
    if (col_sum > norm_a) norm_a = col_sum;
    if (fabs(a->x[col]) > norm_x) norm_x = fabs(a->x[col]);
    if (fabs(a->rhs[col]) > norm_b) norm_b = fabs(a->rhs[col]);
  }
  for (int64_t row = 0; row < a->n; ++row)
    if (fabs(r[row]) > norm_r) norm_r = fabs(r[row]);
  free(r);
  const double scale = norm_a * norm_x + norm_b;
  return scale > 0.0 ? norm_r / scale : norm_r;
}

static double run_case(const synthetic_matrix *source, const char *profile,
                       int threads, uint64_t seed, int quick) {
  double samples[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  const int sample_count = quick ? 3 : 5;
  for (int sample = 0; sample < sample_count; ++sample) {
    synthetic_matrix a = *source;
    memcpy(a.values, a.base, (size_t)a.nnz * sizeof(*a.values));
    kls_options options; kls_default_options(&options);
    options.threads = threads;
    options.ordering = KLS_ORDERING_NATURAL;
    options.expected_refactorizations = 99;
    options.expected_solves = 100;
    options.tuning_profile_path = profile;
    kls_solver *solver = NULL;
    double start = now_seconds();
    int status = kls_create(&solver);
    if (status == KLS_OK)
      status = kls_analyze_csc(solver, KLS_INDEX_INT64, a.n, a.ptr, a.rows,
                               0, &options);
    if (status == KLS_OK) status = kls_factor(solver, a.values);
    const int cycles = quick ? 4 : 20;
    for (int cycle = 0; status == KLS_OK && cycle < cycles; ++cycle) {
      for (int64_t p = 0; p < a.nnz; ++p) {
        const double jitter = 1e-4 *
          ((double)(next_random(&seed) & 1023u) / 1023.0 - 0.5);
        a.values[p] = a.base[p] * (1.0 + jitter);
      }
      status = kls_refactor(solver, a.values);
      if (status == KLS_OK)
        status = kls_solve(solver, 1, a.rhs, a.n, a.x, a.n);
    }
    const double elapsed = now_seconds() - start;
    if (status == KLS_OK && relative_residual(&a) > 1e-8)
      status = KLS_ERR_SOLVE_FAILED;
    samples[sample] = status == KLS_OK ? elapsed : HUGE_VAL;
    kls_destroy(solver);
  }
  return median(samples, sample_count);
}

static int write_candidate(const char *path, int threads,
                           const kls_tuning_values *values, uint64_t seed) {
  kls_tuning_metadata metadata;
  return kls_tuning_write_host_profile(path, threads, values, seed, 1,
                                       &metadata) == KLS_OK;
}

static double score_profile(synthetic_matrix *cases, size_t count,
                            const char *profile, int threads, uint64_t seed,
                            int quick, double *per_case) {
  double log_sum = 0.0;
  for (size_t i = 0; i < count; ++i) {
    per_case[i] = run_case(&cases[i], profile, threads,
                           seed + UINT64_C(0x9e3779b97f4a7c15) * i, quick);
    if (!isfinite(per_case[i]) || per_case[i] <= 0.0) return HUGE_VAL;
    log_sum += log(per_case[i]);
  }
  return exp(log_sum / (double)count);
}

typedef enum field_kind { FIELD_DOUBLE, FIELD_UNSIGNED, FIELD_U64 } field_kind;
typedef struct search_field {
  const char *name; field_kind kind; size_t offset;
} search_field;
#define SD(member) {#member, FIELD_DOUBLE, offsetof(kls_tuning_values, member)}
#define SU(member) {#member, FIELD_UNSIGNED, offsetof(kls_tuning_values, member)}
#define SQ(member) {#member, FIELD_U64, offsetof(kls_tuning_values, member)}
static const search_field search_fields[] = {
  SD(row_dense_min_work), SD(row_compact_panel_min_work),
  SD(row_batch_snode_min_work), SU(dense_help_tail_slice),
  SQ(solve_dense_tail_min_nnz), SQ(solve_parallel_rect_min_nnz),
  SQ(solve_min_entries_per_sync), SD(egraph_min_flops_per_thread),
  SQ(egraph_min_size), SU(egraph_worker_spin_iters),
  SQ(snb_coop_min_doubles)
};

static void scale_field(kls_tuning_values *v, const search_field *f,
                        double scale) {
  unsigned char *p = (unsigned char *)v + f->offset;
  if (f->kind == FIELD_DOUBLE) *(double *)p *= scale;
  else if (f->kind == FIELD_UNSIGNED) {
    unsigned x = (unsigned)((double)*(unsigned *)p * scale);
    *(unsigned *)p = x > 0 ? x : 1;
  } else {
    uint64_t x = (uint64_t)((double)*(uint64_t *)p * scale);
    *(uint64_t *)p = x > 0 ? x : 1;
  }
}

static int pin_cpus(const char *list, int threads) {
  cpu_set_t available, selected;
  if (sched_getaffinity(0, sizeof(available), &available) != 0) return 0;
  CPU_ZERO(&selected);
  if (list != NULL) {
    const char *p = list; int count = 0;
    while (*p != '\0') {
      char *end = NULL; long cpu = strtol(p, &end, 10);
      if (end == p || cpu < 0 || cpu >= CPU_SETSIZE || !CPU_ISSET(cpu, &available)) return 0;
      CPU_SET((int)cpu, &selected); count++;
      p = end; if (*p == ',') p++; else if (*p != '\0') return 0;
    }
    if (count != threads) return 0;
  } else {
    int count = 0;
    for (int cpu = 0; cpu < CPU_SETSIZE && count < threads; ++cpu)
      if (CPU_ISSET(cpu, &available)) { CPU_SET(cpu, &selected); count++; }
    if (count != threads) return 0;
  }
  return sched_setaffinity(0, sizeof(selected), &selected) == 0;
}

static void usage(const char *name) {
  fprintf(stderr, "Usage: %s --threads N --output FILE [--cpus C,C,...] [--seed N] [--force] [--self-test]\n", name);
}

int main(int argc, char **argv) {
#if !defined(__linux__) || !defined(__x86_64__)
  (void)argc; (void)argv;
  fprintf(stderr, "kls-autotune currently supports Linux x86-64 only\n");
  return 2;
#else
  int threads = 0, force = 0, self_test = 0;
  const char *output = NULL, *cpus = NULL;
  uint64_t seed = UINT64_C(0x4b4c5354554e4531);
  for (int i = 1; i < argc; ++i) {
    if (!strcmp(argv[i], "--threads") && i + 1 < argc) threads = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--output") && i + 1 < argc) output = argv[++i];
    else if (!strcmp(argv[i], "--cpus") && i + 1 < argc) cpus = argv[++i];
    else if (!strcmp(argv[i], "--seed") && i + 1 < argc) seed = strtoull(argv[++i], NULL, 0);
    else if (!strcmp(argv[i], "--force")) force = 1;
    else if (!strcmp(argv[i], "--self-test")) self_test = 1;
    else { usage(argv[0]); return 2; }
  }
  if (threads <= 0 || output == NULL || !pin_cpus(cpus, threads)) {
    usage(argv[0]); return 2;
  }
  if (!force && access(output, F_OK) == 0) {
    fprintf(stderr, "%s exists; pass --force to replace it\n", output);
    return 2;
  }
  const int64_t full_sizes[] = {4096, 8192, 16384, 32768, 49152, 65536};
  const int full_widths[] = {1, 2, 4, 8, 16, 24};
  const int64_t self_sizes[] = {96, 128, 160};
  const int self_widths[] = {1, 4, 12};
  const size_t case_count = self_test ? 3u : 6u;
  synthetic_matrix cases[6]; memset(cases, 0, sizeof(cases));
  for (size_t i = 0; i < case_count; ++i) {
    const int64_t size = self_test ? self_sizes[i] : full_sizes[i];
    const int width = self_test ? self_widths[i] : full_widths[i];
    if (!make_stencil(&cases[i], size, width, seed + i)) {
      fprintf(stderr, "failed to create synthetic calibration case\n");
      return 1;
    }
  }
  char candidate_path[4096];
  snprintf(candidate_path, sizeof(candidate_path), "%s.candidate.%ld", output, (long)getpid());
  kls_tuning_values best; kls_tuning_defaults(&best);
  double best_cases[6], trial_cases[6];
  double best_score = score_profile(cases, case_count, NULL, threads, seed,
                                    self_test, best_cases);
  if (!isfinite(best_score)) { fprintf(stderr, "portable baseline failed\n"); return 1; }
  fprintf(stderr, "portable synthetic lifecycle geomean %.6g s\n", best_score);
  if (!self_test) {
    for (int pass = 0; pass < 2; ++pass) {
      for (size_t f = 0; f < sizeof(search_fields) / sizeof(search_fields[0]); ++f) {
        kls_tuning_values accepted = best;
        double accepted_score = best_score;
        for (int side = 0; side < 2; ++side) {
          kls_tuning_values trial = best;
          scale_field(&trial, &search_fields[f], side == 0 ? 0.5 : 2.0);
          if (!write_candidate(candidate_path, threads, &trial, seed)) continue;
          double score = score_profile(cases, case_count, candidate_path, threads,
                                       seed, 0, trial_cases);
          int guarded = isfinite(score) && score <= 0.99 * best_score;
          for (size_t c = 0; c < case_count && guarded; ++c)
            if (trial_cases[c] > 1.025 * best_cases[c]) guarded = 0;
          if (guarded && score < accepted_score) {
            accepted = trial; accepted_score = score;
          }
        }
        if (accepted_score < best_score) {
          fprintf(stderr, "accepted %-38s %.3f%%\n", search_fields[f].name,
                  100.0 * (accepted_score / best_score - 1.0));
          best = accepted; best_score = accepted_score;
          if (write_candidate(candidate_path, threads, &best, seed))
            (void)score_profile(cases, case_count, candidate_path, threads,
                                seed, 0, best_cases);
        }
      }
    }
  }
  unlink(candidate_path);
  kls_tuning_metadata metadata;
  int status = kls_tuning_write_host_profile(output, threads, &best, seed,
                                              force, &metadata);
  for (size_t i = 0; i < case_count; ++i) free_matrix(&cases[i]);
  if (status != KLS_OK) {
    fprintf(stderr, "cannot write profile: %s\n", kls_status_string(status));
    return 1;
  }
  printf("profile=%s id=%016" PRIx64 " tuned_fields=%d threads=%d score=%.9g\n",
         output, metadata.profile_id, metadata.field_count, threads, best_score);
  return 0;
#endif
}
