#define _POSIX_C_SOURCE 200809L

#include "kls/kls.h"

#include <math.h>
#include <inttypes.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int require_ok(int status, const char *what) {
  if (status != KLS_OK) {
    fprintf(stderr, "%s failed: %s (%d)\n", what, kls_status_string(status), status);
    return 0;
  }
  return 1;
}

static int close_enough(double a, double b) {
  return fabs(a - b) < 1e-10;
}

static int require_pivoting_tail_plan(const kls_stats *stats,
                                      const char *what) {
  if (stats == NULL || stats->fast_rejected_pivot < 0 ||
      stats->fast_rejected_block_start < 0 ||
      stats->fast_rejected_block_size <= 0 ||
      stats->fast_rejected_pivoting_tail_columns <= 0 ||
      stats->fast_rejected_pivoting_tail_first <
        stats->fast_rejected_block_start ||
      stats->fast_rejected_pivoting_tail_first >
        stats->fast_rejected_pivot ||
      stats->fast_rejected_pivoting_tail_last <
        stats->fast_rejected_pivot ||
      stats->fast_rejected_pivoting_tail_last >=
        stats->fast_rejected_block_start +
          stats->fast_rejected_block_size ||
      !stats->fast_rejected_pivoting_tail_contains_reject ||
      !stats->fast_rejected_pivoting_tail_topological) {
    fprintf(stderr,
            "unexpected pivoting tail plan for %s: pivot=%" PRId64
            ", block=[%" PRId64 ",%" PRId64 "), cols=%" PRId64
            ", first=%" PRId64 ", last=%" PRId64
            ", contains=%d, topo=%d\n",
            what,
            stats != NULL ? stats->fast_rejected_pivot : -1,
            stats != NULL ? stats->fast_rejected_block_start : -1,
            stats != NULL
              ? stats->fast_rejected_block_start +
                  stats->fast_rejected_block_size
              : -1,
            stats != NULL ? stats->fast_rejected_pivoting_tail_columns : 0,
            stats != NULL ? stats->fast_rejected_pivoting_tail_first : -1,
            stats != NULL ? stats->fast_rejected_pivoting_tail_last : -1,
            stats != NULL
              ? stats->fast_rejected_pivoting_tail_contains_reject
              : 0,
            stats != NULL
              ? stats->fast_rejected_pivoting_tail_topological
              : 0);
    return 0;
  }
  return 1;
}

static int test_csc(void) {
  const int32_t ap[] = {0, 2, 5, 7};
  const int32_t ai[] = {0, 1, 0, 1, 2, 1, 2};
  const double ax[] = {4.0, 1.0, 1.0, 3.0, 1.0, 1.0, 2.0};
  const double b[] = {6.0, 10.0, 8.0};
  double x[3] = {0.0, 0.0, 0.0};

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_AUTO;

  if (!require_ok(kls_create(&solver), "create")) return 0;
  if (!require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 3, ap, ai, 0, &options), "analyze csc")) return 0;
  if (!require_ok(kls_factor(solver, ax), "factor")) return 0;
  if (!require_ok(kls_solve(solver, 1, b, 0, x, 0), "solve")) return 0;

  const int ok = close_enough(x[0], 1.0) && close_enough(x[1], 2.0) && close_enough(x[2], 3.0);
  if (!ok) {
    fprintf(stderr, "unexpected csc solution: %.17g %.17g %.17g\n", x[0], x[1], x[2]);
  }

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (!require_ok(kls_get_stats(solver, &stats), "stats")) return 0;
  if (stats.selected_ordering != KLS_ORDERING_AMD && stats.selected_ordering != KLS_ORDERING_COLAMD) {
    fprintf(stderr, "unexpected selected ordering: %s\n", kls_ordering_name(stats.selected_ordering));
    return 0;
  }
  if (stats.selected_scale != -1) {
    fprintf(stderr, "unexpected auto-selected scale: %d\n", stats.selected_scale);
    return 0;
  }

  kls_destroy(solver);
  return ok;
}

static int test_csr_and_refactor(void) {
  const int64_t rp[] = {0, 2, 5, 7};
  const int64_t ci[] = {0, 1, 0, 1, 2, 1, 2};
  const double ax[] = {4.0, 1.0, 1.0, 3.0, 1.0, 1.0, 2.0};
  const double b[] = {6.0, 10.0, 8.0};
  double x[3] = {0.0, 0.0, 0.0};

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_NATURAL;
  options.use_btf = 0;

  if (!require_ok(kls_create(&solver), "create")) return 0;
  if (!require_ok(kls_analyze_csr(solver, KLS_INDEX_INT64, 3, rp, ci, 0, &options), "analyze csr")) return 0;
  if (!require_ok(kls_factor(solver, ax), "factor csr")) return 0;
  if (!require_ok(kls_refactor(solver, ax), "refactor csr")) return 0;
  if (!require_ok(kls_solve(solver, 1, b, 0, x, 0), "solve csr")) return 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (!require_ok(kls_get_stats(solver, &stats), "stats csr")) return 0;
  if (stats.selected_btf != 0) {
    fprintf(stderr, "unexpected selected btf for no-btf solve: %d\n", stats.selected_btf);
    return 0;
  }

  const int ok = close_enough(x[0], 1.0) && close_enough(x[1], 2.0) && close_enough(x[2], 3.0);
  if (!ok) {
    fprintf(stderr, "unexpected csr solution: %.17g %.17g %.17g\n", x[0], x[1], x[2]);
  }
  kls_destroy(solver);
  return ok;
}

static int test_sparse_diagonal_auto_scale(void) {
  const int32_t ap[] = {0, 1, 3, 4};
  const int32_t ai[] = {1, 0, 1, 2};
  const double ax[] = {1.0, 1.0, 2.0, 3.0};
  const double b[] = {2.0, 5.0, 9.0};
  double x[3] = {0.0, 0.0, 0.0};

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_NATURAL;

  if (!require_ok(kls_create(&solver), "create")) return 0;
  if (!require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 3, ap, ai, 0, &options),
                  "analyze sparse diagonal")) return 0;
  if (!require_ok(kls_factor(solver, ax), "factor sparse diagonal")) return 0;
  if (!require_ok(kls_solve(solver, 1, b, 0, x, 0), "solve sparse diagonal")) return 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (!require_ok(kls_get_stats(solver, &stats), "stats sparse diagonal")) return 0;
  if (stats.selected_scale != 1) {
    fprintf(stderr, "unexpected sparse-diagonal auto scale: %d\n", stats.selected_scale);
    return 0;
  }

  const int ok = close_enough(x[0], 1.0) && close_enough(x[1], 2.0) && close_enough(x[2], 3.0);
  if (!ok) {
    fprintf(stderr, "unexpected sparse-diagonal solution: %.17g %.17g %.17g\n",
            x[0], x[1], x[2]);
  }
  kls_destroy(solver);
  return ok;
}

static int test_forced_transpose_orientation(void) {
  const int32_t ap[] = {0, 2, 5, 8};
  const int32_t ai[] = {0, 1, 0, 1, 2, 0, 1, 2};
  const double ax[] = {4.0, 1.0, 2.0, 3.0, 5.0, 7.0, 1.0, 2.0};
  const double b[] = {29.0, 10.0, 16.0};
  const double bt[] = {7.0, 3.5, 14.0};
  double x[3] = {0.0, 0.0, 0.0};
  double xt[3] = {0.0, 0.0, 0.0};

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_NATURAL;
  options.orientation = KLS_ORIENTATION_TRANSPOSE;

  if (!require_ok(kls_create(&solver), "create")) return 0;
  if (!require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 3, ap, ai, 0, &options),
                  "analyze csc transpose")) return 0;
  if (!require_ok(kls_factor(solver, ax), "factor transpose")) return 0;
  if (!require_ok(kls_solve(solver, 1, b, 0, x, 0), "solve transpose-oriented")) return 0;
  if (!require_ok(kls_solve_transpose(solver, 1, bt, 0, xt, 0), "transpose solve")) return 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (!require_ok(kls_get_stats(solver, &stats), "stats transpose")) return 0;
  if (stats.selected_orientation != KLS_ORIENTATION_TRANSPOSE) {
    fprintf(stderr, "unexpected selected orientation: %s\n",
            kls_orientation_name(stats.selected_orientation));
    kls_destroy(solver);
    return 0;
  }

  const int ok =
    close_enough(x[0], 1.0) && close_enough(x[1], 2.0) && close_enough(x[2], 3.0) &&
    close_enough(xt[0], 2.0) && close_enough(xt[1], -1.0) && close_enough(xt[2], 0.5);
  if (!ok) {
    fprintf(stderr,
            "unexpected transpose-oriented solutions: x=(%.17g %.17g %.17g), xt=(%.17g %.17g %.17g)\n",
            x[0], x[1], x[2], xt[0], xt[1], xt[2]);
  }
  kls_destroy(solver);
  return ok;
}

static int test_csr_forced_transpose_orientation(void) {
  const int64_t rp[] = {0, 3, 6, 8};
  const int64_t ci[] = {0, 1, 2, 0, 1, 2, 1, 2};
  const double ax[] = {4.0, 2.0, 7.0, 1.0, 3.0, 1.0, 5.0, 2.0};
  const double b[] = {29.0, 10.0, 16.0};
  double x[3] = {0.0, 0.0, 0.0};

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_NATURAL;
  options.orientation = KLS_ORIENTATION_TRANSPOSE;

  if (!require_ok(kls_create(&solver), "create")) return 0;
  if (!require_ok(kls_analyze_csr(solver, KLS_INDEX_INT64, 3, rp, ci, 0, &options),
                  "analyze csr transpose")) return 0;
  if (!require_ok(kls_factor(solver, ax), "factor csr transpose")) return 0;
  if (!require_ok(kls_refactor(solver, ax), "refactor csr transpose")) return 0;
  if (!require_ok(kls_solve(solver, 1, b, 0, x, 0), "solve csr transpose-oriented")) return 0;

  const int ok = close_enough(x[0], 1.0) && close_enough(x[1], 2.0) && close_enough(x[2], 3.0);
  if (!ok) {
    fprintf(stderr, "unexpected csr transpose-oriented solution: %.17g %.17g %.17g\n",
            x[0], x[1], x[2]);
  }
  kls_destroy(solver);
  return ok;
}

static int test_solve_strides_and_in_place(void) {
  const int32_t ap[] = {0, 2, 5, 7};
  const int32_t ai[] = {0, 1, 0, 1, 2, 1, 2};
  const double ax[] = {4.0, 1.0, 1.0, 3.0, 1.0, 1.0, 2.0};
  const double b_strided[] = {6.0, 10.0, 8.0, -1.0, 11.0, 15.0, 11.0, -2.0};
  double x_strided[] = {0.0, 0.0, 0.0, 77.0, 0.0, 0.0, 0.0, 88.0};
  double in_place[] = {6.0, 10.0, 8.0};

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_NATURAL;

  if (!require_ok(kls_create(&solver), "create")) return 0;
  if (!require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 3, ap, ai, 0, &options),
                  "analyze strides")) return 0;
  if (!require_ok(kls_factor(solver, ax), "factor strides")) return 0;
  if (!require_ok(kls_solve(solver, 2, b_strided, 4, x_strided, 4),
                  "solve strided")) return 0;
  if (!require_ok(kls_solve(solver, 1, in_place, 0, in_place, 0),
                  "solve in-place")) return 0;

  const int ok =
    close_enough(x_strided[0], 1.0) && close_enough(x_strided[1], 2.0) &&
    close_enough(x_strided[2], 3.0) && close_enough(x_strided[3], 77.0) &&
    close_enough(x_strided[4], 2.0) && close_enough(x_strided[5], 3.0) &&
    close_enough(x_strided[6], 4.0) && close_enough(x_strided[7], 88.0) &&
    close_enough(in_place[0], 1.0) && close_enough(in_place[1], 2.0) &&
    close_enough(in_place[2], 3.0);
  if (!ok) {
    fprintf(stderr,
            "unexpected strided/in-place solution: xs=(%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g), xi=(%.17g %.17g %.17g)\n",
            x_strided[0], x_strided[1], x_strided[2], x_strided[3],
            x_strided[4], x_strided[5], x_strided[6], x_strided[7],
            in_place[0], in_place[1], in_place[2]);
  }
  kls_destroy(solver);
  return ok;
}

static int test_fast_factor_pivot_check_fallback(void) {
  const int32_t ap[] = {0, 2, 4};
  const int32_t ai[] = {0, 1, 0, 1};
  const double ax0[] = {2.0, 1.0, 1.0, 2.0};
  const double ax1[] = {1.0e-12, 1.0, 1.0, 2.0};
  const double b[] = {2.000000000001, 5.0};
  double x[2] = {0.0, 0.0};

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_NATURAL;
  options.use_btf = 0;
  options.scale = -1;
  options.pivot_tolerance = 0.001;
  options.static_pivoting = 0;

  int ok = 1;
  if (!require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 2, ap, ai, 0,
                                        &options),
                        "analyze pivot-check fallback")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax0), "factor pivot-check base")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax1), "factor pivot-check fallback")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve pivot-check fallback")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats), "stats pivot-check fallback")) {
    ok = 0;
  }
  if (ok && stats.offdiag_pivots < 1) {
    fprintf(stderr, "fast factor did not fall back to pivoting factorization\n");
    ok = 0;
  }
  if (ok && (stats.fast_rejected_pivot != 0 || stats.fast_rejected_pivot_col != 0)) {
    fprintf(stderr,
            "unexpected fast rejected pivot: pivot=%" PRId64 ", col=%" PRId64 "\n",
            stats.fast_rejected_pivot, stats.fast_rejected_pivot_col);
    ok = 0;
  }
  if (ok && stats.fast_block_restarts != 1) {
    fprintf(stderr, "unexpected fast block restarts: %d\n",
            stats.fast_block_restarts);
    ok = 0;
  }
  if (ok && stats.fast_rejected_refresh_state !=
              KLS_FAST_REJECT_REFRESH_PREFIX) {
    fprintf(stderr, "unexpected fast reject refresh state: %d\n",
            stats.fast_rejected_refresh_state);
    ok = 0;
  }
  if (ok && (stats.fast_rejected_block_start != 0 ||
             stats.fast_rejected_block_size != 2 ||
             stats.fast_rejected_suffix_columns != 2 ||
             stats.fast_rejected_descendant_columns < 1 ||
             stats.fast_rejected_descendant_columns >
               stats.fast_rejected_suffix_columns ||
             stats.fast_rejected_etree_columns < 1 ||
             stats.fast_rejected_etree_columns >
               stats.fast_rejected_suffix_columns ||
             stats.fast_rejected_pivoting_tail_columns < 1 ||
             stats.fast_rejected_pivoting_tail_columns >
               stats.fast_rejected_suffix_columns)) {
    fprintf(stderr,
            "unexpected fast reject tail stats: start=%" PRId64
            ", size=%" PRId64 ", suffix=%" PRId64
            ", descendants=%" PRId64 ", etree=%" PRId64
            ", pivoting_tail=%" PRId64 "\n",
            stats.fast_rejected_block_start,
            stats.fast_rejected_block_size,
            stats.fast_rejected_suffix_columns,
            stats.fast_rejected_descendant_columns,
            stats.fast_rejected_etree_columns,
            stats.fast_rejected_pivoting_tail_columns);
    ok = 0;
  }
  if (ok && !require_pivoting_tail_plan(&stats, "pivot-check fallback")) {
    ok = 0;
  }
  if (ok && (!close_enough(x[0], 1.0) || !close_enough(x[1], 2.0))) {
    fprintf(stderr, "unexpected pivot-check fallback solution: %.17g %.17g\n",
            x[0], x[1]);
    ok = 0;
  }

  kls_destroy(solver);
  return ok;
}

static int test_checked_row_fast_factor_block_restart(void) {
  const int32_t ap[] = {0, 2, 4};
  const int32_t ai[] = {0, 1, 0, 1};
  const double ax0[] = {2.0, 1.0, 1.0, 2.0};
  const double ax1[] = {1.0e-12, 1.0, 1.0, 2.0};
  const double b[] = {2.000000000001, 5.0};
  double x[2] = {0.0, 0.0};

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_NATURAL;
  options.use_btf = 0;
  options.scale = -1;
  options.pivot_tolerance = 0.001;
  options.static_pivoting = 0;
  const char *saved_env_value = getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  char *saved_env = saved_env_value != NULL ? strdup(saved_env_value) : NULL;
  const int had_saved_env = saved_env_value != NULL;

  int ok = 1;
  if (had_saved_env && saved_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_CHECKED_ROW_REFACTOR\n");
    ok = 0;
  }
  if (!require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 2, ap, ai, 0,
                                        &options),
                        "analyze checked-row block restart")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor checked-row base")) ok = 0;
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(solver, ax1),
                        "factor checked-row block restart")) ok = 0;
  if (had_saved_env && saved_env != NULL) {
    if (setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", saved_env, 1) != 0) {
      perror("restore KLS_ENABLE_CHECKED_ROW_REFACTOR");
      ok = 0;
    }
  } else if (!had_saved_env) {
    if (unsetenv("KLS_ENABLE_CHECKED_ROW_REFACTOR") != 0) {
      perror("unsetenv KLS_ENABLE_CHECKED_ROW_REFACTOR");
      ok = 0;
    }
  }
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve checked-row block restart")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats checked-row block restart")) {
    ok = 0;
  }
  if (ok && stats.offdiag_pivots < 1) {
    fprintf(stderr, "checked-row restart did not pivot the repaired block\n");
    ok = 0;
  }
  if (ok && (stats.fast_rejected_pivot != 0 ||
             stats.fast_rejected_pivot_col != 0 ||
             stats.fast_rejected_row != 1 ||
             stats.fast_rejected_multiplier_abs <= 1.0e6 ||
             stats.fast_rejected_pivot_abs >= 1.0e-9 ||
             stats.fast_rejected_candidate_abs <= 0.5 ||
             stats.fast_rejected_tail_candidate_row != 1 ||
             stats.fast_rejected_tail_candidate_abs <= 0.5 ||
             stats.fast_rejected_tail_candidate_count < 1 ||
             stats.fast_rejected_tail_candidate_position < 0 ||
             stats.fast_rejected_tail_repair_ready != 1 ||
             stats.fast_repaired_pivot_row != 1 ||
             stats.fast_repaired_pivot_matches_tail_candidate != 1 ||
             stats.fast_repaired_first_changed_pivot != 0 ||
             stats.fast_repaired_prefix_changed_pivots != 0 ||
             stats.fast_repaired_suffix_changed_pivots < 1 ||
             stats.fast_repaired_tail_restart_ready != 0 ||
             stats.fast_repaired_block_work <= 0.0 ||
             stats.fast_repaired_tail_restart_columns != 0 ||
             stats.fast_repaired_tail_restart_work != 0.0 ||
             stats.fast_block_restarts != 1 ||
             stats.fast_tail_restarts != 0 ||
             stats.fast_rejected_refresh_state !=
               KLS_FAST_REJECT_REFRESH_PREFIX)) {
    fprintf(stderr,
            "unexpected checked-row reject stats: pivot=%" PRId64
            ", col=%" PRId64 ", row=%" PRId64 ", |L|=%.6g"
            ", |pivot|=%.6g, |candidate|=%.6g"
            ", tail_row=%" PRId64 ", tail_candidate=%.6g"
            ", tail_count=%" PRId64 ", tail_pos=%" PRId64
            ", tail_ready=%d, repaired_row=%" PRId64
            ", repaired_match=%d, first_changed=%" PRId64
            ", prefix_changed=%" PRId64 ", suffix_changed=%" PRId64
            ", tail_restart_ready=%d, block_work=%.6g"
            ", tail_cols=%" PRId64 ", tail_work=%.6g"
            ", saved_work=%.6g, restarts=%d, tail_restarts=%d"
            ", refresh=%d\n",
            stats.fast_rejected_pivot, stats.fast_rejected_pivot_col,
            stats.fast_rejected_row, stats.fast_rejected_multiplier_abs,
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
            stats.fast_block_restarts, stats.fast_tail_restarts,
            stats.fast_rejected_refresh_state);
    ok = 0;
  }
  if (ok && (stats.fast_rejected_block_start != 0 ||
             stats.fast_rejected_block_size != 2 ||
             stats.fast_rejected_suffix_columns != 2 ||
             stats.fast_rejected_descendant_columns < 1 ||
             stats.fast_rejected_descendant_columns >
               stats.fast_rejected_suffix_columns ||
             stats.fast_rejected_row_tail_columns < 1 ||
             stats.fast_rejected_row_tail_columns >
               stats.fast_rejected_suffix_columns ||
             stats.fast_rejected_row_tail_work <= 0.0 ||
             stats.fast_rejected_etree_columns < 1 ||
             stats.fast_rejected_etree_columns >
               stats.fast_rejected_suffix_columns ||
             stats.fast_rejected_pivoting_tail_columns < 1 ||
             stats.fast_rejected_pivoting_tail_columns >
               stats.fast_rejected_suffix_columns)) {
    fprintf(stderr,
            "unexpected checked-row tail stats: start=%" PRId64
            ", size=%" PRId64 ", suffix=%" PRId64
            ", descendants=%" PRId64 ", row_tail=%" PRId64
            ", etree=%" PRId64 ", pivoting_tail=%" PRId64 "\n",
            stats.fast_rejected_block_start,
            stats.fast_rejected_block_size,
            stats.fast_rejected_suffix_columns,
            stats.fast_rejected_descendant_columns,
            stats.fast_rejected_row_tail_columns,
            stats.fast_rejected_etree_columns,
            stats.fast_rejected_pivoting_tail_columns);
    ok = 0;
  }
  if (ok && !require_pivoting_tail_plan(&stats, "checked-row block restart")) {
    ok = 0;
  }
  if (ok && (!close_enough(x[0], 1.0) || !close_enough(x[1], 2.0))) {
    fprintf(stderr,
            "unexpected checked-row restart solution: %.17g %.17g\n",
            x[0], x[1]);
    ok = 0;
  }

  kls_destroy(solver);
  free(saved_env);
  return ok;
}

static int test_parallel_checked_row_fast_factor_block_restart(void) {
  const int32_t ap[] = {0, 1, 3, 5, 6};
  const int32_t ai[] = {0, 1, 2, 1, 2, 3};
  const double ax0[] = {2.0, 2.0, 1.0, 1.0, 2.0, 4.0};
  const double ax1[] = {2.0, 1.0e-12, 1.0, 1.0, 2.0, 4.0};
  const double b[] = {2.0, 3.000000000002, 8.0, 16.0};
  double x[4] = {0.0, 0.0, 0.0, 0.0};

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.threads = 2;
  options.ordering = KLS_ORDERING_NATURAL;
  options.use_btf = 0;
  options.scale = -1;
  options.pivot_tolerance = 0.001;
  options.static_pivoting = 0;
  const char *saved_env_value = getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  char *saved_env = saved_env_value != NULL ? strdup(saved_env_value) : NULL;
  const int had_saved_env = saved_env_value != NULL;

  int ok = 1;
  if (had_saved_env && saved_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_CHECKED_ROW_REFACTOR\n");
    ok = 0;
  }
  if (!require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 4, ap, ai, 0,
                                        &options),
                        "analyze parallel checked-row restart")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor parallel checked-row base")) ok = 0;
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(solver, ax1),
                        "factor parallel checked-row repair")) ok = 0;
  if (had_saved_env && saved_env != NULL) {
    if (setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", saved_env, 1) != 0) {
      perror("restore KLS_ENABLE_CHECKED_ROW_REFACTOR");
      ok = 0;
    }
  } else if (!had_saved_env) {
    if (unsetenv("KLS_ENABLE_CHECKED_ROW_REFACTOR") != 0) {
      perror("unsetenv KLS_ENABLE_CHECKED_ROW_REFACTOR");
      ok = 0;
    }
  }
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve parallel checked-row repair")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats parallel checked-row restart")) {
    ok = 0;
  }
  if (ok && (stats.fast_rejected_pivot != 1 ||
             stats.fast_rejected_pivot_col != 1 ||
             stats.fast_rejected_row != 2 ||
             stats.fast_rejected_refresh_state !=
               KLS_FAST_REJECT_REFRESH_PREFIX ||
             stats.fast_rejected_tail_repair_ready != 1 ||
             stats.fast_block_restarts != 1 ||
             stats.fast_tail_restarts != 1 ||
             stats.fast_rejected_pivoting_tail_columns < 1 ||
             stats.fast_rejected_pivoting_tail_topological != 1)) {
    fprintf(stderr,
            "unexpected parallel checked-row stats: pivot=%" PRId64
            ", col=%" PRId64 ", row=%" PRId64 ", refresh=%d"
            ", tail_ready=%d, restarts=%d, tail_restarts=%d"
            ", tail_cols=%" PRId64
            ", tail_topo=%d\n",
            stats.fast_rejected_pivot,
            stats.fast_rejected_pivot_col,
            stats.fast_rejected_row,
            stats.fast_rejected_refresh_state,
            stats.fast_rejected_tail_repair_ready,
            stats.fast_block_restarts,
            stats.fast_tail_restarts,
            stats.fast_rejected_pivoting_tail_columns,
            stats.fast_rejected_pivoting_tail_topological);
    ok = 0;
  }
  if (ok && !require_pivoting_tail_plan(&stats,
                                        "parallel checked-row restart")) {
    ok = 0;
  }
  if (ok && (!close_enough(x[0], 1.0) || !close_enough(x[1], 2.0) ||
             !close_enough(x[2], 3.0) || !close_enough(x[3], 4.0))) {
    fprintf(stderr,
            "unexpected parallel checked-row solution: %.17g %.17g %.17g %.17g\n",
            x[0], x[1], x[2], x[3]);
    ok = 0;
  }

  kls_destroy(solver);
  free(saved_env);
  return ok;
}

static int test_fast_factor_tail_prefix_state_validation(void) {
  const int32_t ap[] = {0, 1, 3, 5};
  const int32_t ai[] = {0, 1, 2, 1, 2};
  const double ax0[] = {2.0, 2.0, 1.0, 1.0, 2.0};
  const double ax1[] = {2.0, 1.0e-12, 1.0, 1.0, 2.0};
  const double b[] = {2.0, 3.000000000002, 8.0};
  double x[3] = {0.0, 0.0, 0.0};

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_NATURAL;
  options.use_btf = 0;
  options.scale = -1;
  options.pivot_tolerance = 0.001;
  const char *saved_env_value = getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  char *saved_env = saved_env_value != NULL ? strdup(saved_env_value) : NULL;
  const int had_saved_env = saved_env_value != NULL;

  int ok = 1;
  if (had_saved_env && saved_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_CHECKED_ROW_REFACTOR\n");
    ok = 0;
  }
  if (!require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 3, ap, ai, 0,
                                        &options),
                        "analyze tail-prefix validation")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor tail-prefix base")) ok = 0;
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(solver, ax1),
                        "factor tail-prefix repair")) ok = 0;
  if (had_saved_env && saved_env != NULL) {
    if (setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", saved_env, 1) != 0) {
      perror("restore KLS_ENABLE_CHECKED_ROW_REFACTOR");
      ok = 0;
    }
  } else if (!had_saved_env) {
    if (unsetenv("KLS_ENABLE_CHECKED_ROW_REFACTOR") != 0) {
      perror("unsetenv KLS_ENABLE_CHECKED_ROW_REFACTOR");
      ok = 0;
    }
  }
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve tail-prefix repair")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats tail-prefix validation")) {
    ok = 0;
  }
  if (ok && (stats.fast_rejected_pivot != 1 ||
             stats.fast_rejected_block_start != 0 ||
             stats.fast_rejected_block_size != 3 ||
             stats.fast_rejected_suffix_columns != 2 ||
             stats.fast_repaired_prefix_changed_pivots != 0 ||
             stats.fast_repaired_suffix_changed_pivots < 1 ||
             stats.fast_repaired_tail_restart_ready != 1 ||
             stats.fast_repaired_tail_restart_columns < 1 ||
             stats.fast_repaired_tail_restart_work <= 0.0 ||
             stats.fast_block_restarts != 1 ||
             stats.fast_tail_restarts != 1)) {
    fprintf(stderr,
            "unexpected tail-prefix validation stats: pivot=%" PRId64
            ", start=%" PRId64 ", size=%" PRId64 ", suffix=%" PRId64
            ", prefix_changed=%" PRId64 ", suffix_changed=%" PRId64
            ", tail_candidate=%" PRId64 ", tail_candidate_count=%" PRId64
            ", tail_repair_ready=%d, repaired_row=%" PRId64
            ", repaired_match=%d"
            ", tail_ready=%d, tail_cols=%" PRId64 ", tail_work=%.6g"
            ", pivoting_tail_cols=%" PRId64 ", pivoting_tail_work=%.6g"
            ", restarts=%d, tail_restarts=%d\n",
            stats.fast_rejected_pivot,
            stats.fast_rejected_block_start,
            stats.fast_rejected_block_size,
            stats.fast_rejected_suffix_columns,
            stats.fast_repaired_prefix_changed_pivots,
            stats.fast_repaired_suffix_changed_pivots,
            stats.fast_rejected_tail_candidate_row,
            stats.fast_rejected_tail_candidate_count,
            stats.fast_rejected_tail_repair_ready,
            stats.fast_repaired_pivot_row,
            stats.fast_repaired_pivot_matches_tail_candidate,
            stats.fast_repaired_tail_restart_ready,
            stats.fast_repaired_tail_restart_columns,
            stats.fast_repaired_tail_restart_work,
            stats.fast_rejected_pivoting_tail_columns,
            stats.fast_rejected_pivoting_tail_work,
            stats.fast_block_restarts,
            stats.fast_tail_restarts);
    ok = 0;
  }
  if (ok && !require_pivoting_tail_plan(&stats, "tail-prefix validation")) {
    ok = 0;
  }
  if (ok && (!close_enough(x[0], 1.0) || !close_enough(x[1], 2.0) ||
             !close_enough(x[2], 3.0))) {
    fprintf(stderr,
            "unexpected tail-prefix solution: %.17g %.17g %.17g\n",
            x[0], x[1], x[2]);
    ok = 0;
  }

  kls_destroy(solver);
  free(saved_env);
  return ok;
}

static int test_mapped_fast_factor_prefix_tail_restart(void) {
  const int32_t ap[] = {0, 1, 3, 5};
  const int32_t ai[] = {0, 1, 2, 1, 2};
  const double ax0[] = {2.0, 2.0, 1.0, 1.0, 2.0};
  const double ax1[] = {2.0, 1.0e-12, 1.0, 1.0, 2.0};
  const double b[] = {2.0, 3.000000000002, 8.0};
  double x[3] = {0.0, 0.0, 0.0};

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.threads = 1;
  options.ordering = KLS_ORDERING_NATURAL;
  options.use_btf = 0;
  options.scale = -1;
  options.pivot_tolerance = 0.001;
  options.static_pivoting = 0;
  const char *saved_env_value = getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  char *saved_env = saved_env_value != NULL ? strdup(saved_env_value) : NULL;
  const int had_saved_env = saved_env_value != NULL;

  int ok = 1;
  if (had_saved_env && saved_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_CHECKED_ROW_REFACTOR\n");
    ok = 0;
  }
  if (!require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 3, ap, ai, 0,
                                        &options),
                        "analyze mapped prefix-tail restart")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor mapped prefix-tail base")) ok = 0;
  if (ok && unsetenv("KLS_ENABLE_CHECKED_ROW_REFACTOR") != 0) {
    perror("unsetenv KLS_ENABLE_CHECKED_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(solver, ax1),
                        "factor mapped prefix-tail repair")) ok = 0;
  if (had_saved_env && saved_env != NULL) {
    if (setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", saved_env, 1) != 0) {
      perror("restore KLS_ENABLE_CHECKED_ROW_REFACTOR");
      ok = 0;
    }
  } else if (!had_saved_env) {
    if (unsetenv("KLS_ENABLE_CHECKED_ROW_REFACTOR") != 0) {
      perror("unsetenv KLS_ENABLE_CHECKED_ROW_REFACTOR");
      ok = 0;
    }
  }
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve mapped prefix-tail repair")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats mapped prefix-tail restart")) {
    ok = 0;
  }
  if (ok && (stats.fast_rejected_pivot != 1 ||
             stats.fast_rejected_block_start != 0 ||
             stats.fast_rejected_block_size != 3 ||
             stats.fast_rejected_refresh_state !=
               KLS_FAST_REJECT_REFRESH_PREFIX ||
             stats.fast_repaired_prefix_changed_pivots != 0 ||
             stats.fast_repaired_suffix_changed_pivots < 1 ||
             stats.fast_repaired_tail_restart_ready != 1 ||
             stats.fast_repaired_tail_restart_columns < 1 ||
             stats.fast_repaired_tail_restart_work <= 0.0 ||
             stats.fast_block_restarts != 1 ||
             stats.fast_tail_restarts != 1)) {
    fprintf(stderr,
            "unexpected mapped prefix-tail restart stats: pivot=%" PRId64
            ", start=%" PRId64 ", size=%" PRId64 ", refresh=%d"
            ", tail_candidate=%" PRId64 ", tail_candidate_count=%" PRId64
            ", tail_repair_ready=%d, repaired_match=%d"
            ", prefix_changed=%" PRId64 ", suffix_changed=%" PRId64
            ", tail_ready=%d, tail_cols=%" PRId64 ", tail_work=%.6g"
            ", block_restarts=%d, tail_restarts=%d\n",
            stats.fast_rejected_pivot,
            stats.fast_rejected_block_start,
            stats.fast_rejected_block_size,
            stats.fast_rejected_refresh_state,
            stats.fast_rejected_tail_candidate_row,
            stats.fast_rejected_tail_candidate_count,
            stats.fast_rejected_tail_repair_ready,
            stats.fast_repaired_pivot_matches_tail_candidate,
            stats.fast_repaired_prefix_changed_pivots,
            stats.fast_repaired_suffix_changed_pivots,
            stats.fast_repaired_tail_restart_ready,
            stats.fast_repaired_tail_restart_columns,
            stats.fast_repaired_tail_restart_work,
            stats.fast_block_restarts,
            stats.fast_tail_restarts);
    ok = 0;
  }
  if (ok && !require_pivoting_tail_plan(&stats,
                                        "mapped prefix-tail restart")) {
    ok = 0;
  }
  if (ok && (!close_enough(x[0], 1.0) || !close_enough(x[1], 2.0) ||
             !close_enough(x[2], 3.0))) {
    fprintf(stderr,
            "unexpected mapped prefix-tail solution: %.17g %.17g %.17g\n",
            x[0], x[1], x[2]);
    ok = 0;
  }

  kls_destroy(solver);
  free(saved_env);
  return ok;
}

static int test_scaled_fast_factor_block_restart(void) {
  const int32_t ap[] = {0, 2, 4};
  const int32_t ai[] = {0, 1, 0, 1};
  const double ax0[] = {2.0, 1.0, 1.0, 2.0};
  const double ax1[] = {1.0e-12, 1.0, 1.0, 2.0};
  const double b[] = {2.000000000001, 5.0};
  double x[2] = {0.0, 0.0};

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_NATURAL;
  options.use_btf = 0;
  options.scale = 2;
  options.pivot_tolerance = 0.001;

  int ok = 1;
  if (!require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 2, ap, ai, 0,
                                        &options),
                        "analyze scaled pivot-check fallback")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor scaled pivot-check base")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax1),
                        "factor scaled pivot-check block restart")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve scaled pivot-check block restart")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats scaled pivot-check block restart")) {
    ok = 0;
  }
  if (ok && stats.offdiag_pivots < 1) {
    fprintf(stderr, "scaled block restart did not pivot the rejected block\n");
    ok = 0;
  }
  if (ok && (stats.fast_rejected_pivot != 0 || stats.fast_rejected_pivot_col != 0)) {
    fprintf(stderr,
            "unexpected scaled fast rejected pivot: pivot=%" PRId64
            ", col=%" PRId64 "\n",
            stats.fast_rejected_pivot, stats.fast_rejected_pivot_col);
    ok = 0;
  }
  if (ok && stats.fast_block_restarts != 1) {
    fprintf(stderr, "scaled path restart count was %d\n",
            stats.fast_block_restarts);
    ok = 0;
  }
  if (ok && stats.fast_rejected_refresh_state !=
              KLS_FAST_REJECT_REFRESH_PREFIX) {
    fprintf(stderr, "unexpected scaled reject refresh state: %d\n",
            stats.fast_rejected_refresh_state);
    ok = 0;
  }
  if (ok && (stats.fast_rejected_block_start != 0 ||
             stats.fast_rejected_block_size != 2 ||
             stats.fast_rejected_suffix_columns != 2 ||
             stats.fast_rejected_descendant_columns < 1 ||
             stats.fast_rejected_descendant_columns >
               stats.fast_rejected_suffix_columns ||
             stats.fast_rejected_etree_columns < 1 ||
             stats.fast_rejected_etree_columns >
               stats.fast_rejected_suffix_columns ||
             stats.fast_rejected_pivoting_tail_columns < 1 ||
             stats.fast_rejected_pivoting_tail_columns >
               stats.fast_rejected_suffix_columns)) {
    fprintf(stderr,
            "unexpected scaled reject tail stats: start=%" PRId64
            ", size=%" PRId64 ", suffix=%" PRId64
            ", descendants=%" PRId64 ", etree=%" PRId64
            ", pivoting_tail=%" PRId64 "\n",
            stats.fast_rejected_block_start,
            stats.fast_rejected_block_size,
            stats.fast_rejected_suffix_columns,
            stats.fast_rejected_descendant_columns,
            stats.fast_rejected_etree_columns,
            stats.fast_rejected_pivoting_tail_columns);
    ok = 0;
  }
  if (ok && !require_pivoting_tail_plan(&stats, "scaled block restart")) {
    ok = 0;
  }
  if (ok && (!close_enough(x[0], 1.0) || !close_enough(x[1], 2.0))) {
    fprintf(stderr, "unexpected scaled pivot-check restart solution: %.17g %.17g\n",
            x[0], x[1]);
    ok = 0;
  }

  kls_destroy(solver);
  return ok;
}

static int test_scaled_fast_factor_prefix_tail_restart(void) {
  const int32_t ap[] = {0, 1, 3, 5};
  const int32_t ai[] = {0, 1, 2, 1, 2};
  const double ax0[] = {2.0, 2.0, 1.0, 1.0, 2.0};
  const double ax1[] = {2.0, 1.0e-12, 1.0, 1.0, 2.0};
  const double b[] = {2.0, 3.000000000002, 8.0};
  double x[3] = {0.0, 0.0, 0.0};

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_NATURAL;
  options.use_btf = 0;
  options.scale = 2;
  options.pivot_tolerance = 0.001;

  int ok = 1;
  if (!require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 3, ap, ai, 0,
                                        &options),
                        "analyze scaled prefix-tail restart")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor scaled prefix-tail base")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax1),
                        "factor scaled prefix-tail repair")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve scaled prefix-tail repair")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats scaled prefix-tail restart")) {
    ok = 0;
  }
  if (ok && (stats.fast_rejected_pivot != 1 ||
             stats.fast_rejected_block_start != 0 ||
             stats.fast_rejected_block_size != 3 ||
             stats.fast_rejected_refresh_state !=
               KLS_FAST_REJECT_REFRESH_PREFIX ||
             stats.fast_repaired_prefix_changed_pivots != 0 ||
             stats.fast_repaired_suffix_changed_pivots < 1 ||
             stats.fast_repaired_tail_restart_ready != 1 ||
             stats.fast_repaired_tail_restart_columns < 1 ||
             stats.fast_repaired_tail_restart_work <= 0.0 ||
             stats.fast_block_restarts != 1 ||
             stats.fast_tail_restarts != 1)) {
    fprintf(stderr,
            "unexpected scaled prefix-tail stats: pivot=%" PRId64
            ", start=%" PRId64 ", size=%" PRId64 ", refresh=%d"
            ", prefix_changed=%" PRId64 ", suffix_changed=%" PRId64
            ", tail_ready=%d, tail_cols=%" PRId64 ", tail_work=%.6g"
            ", block_restarts=%d, tail_restarts=%d\n",
            stats.fast_rejected_pivot,
            stats.fast_rejected_block_start,
            stats.fast_rejected_block_size,
            stats.fast_rejected_refresh_state,
            stats.fast_repaired_prefix_changed_pivots,
            stats.fast_repaired_suffix_changed_pivots,
            stats.fast_repaired_tail_restart_ready,
            stats.fast_repaired_tail_restart_columns,
            stats.fast_repaired_tail_restart_work,
            stats.fast_block_restarts,
            stats.fast_tail_restarts);
    ok = 0;
  }
  if (ok && !require_pivoting_tail_plan(&stats,
                                        "scaled prefix-tail restart")) {
    ok = 0;
  }
  if (ok && (!close_enough(x[0], 1.0) || !close_enough(x[1], 2.0) ||
             !close_enough(x[2], 3.0))) {
    fprintf(stderr,
            "unexpected scaled prefix-tail solution: %.17g %.17g %.17g\n",
            x[0], x[1], x[2]);
    ok = 0;
  }

  kls_destroy(solver);
  return ok;
}

static int test_scaled_btf_fast_factor_tail_continuation(void) {
  const int32_t ap[] = {0, 2, 4, 7, 9};
  const int32_t ai[] = {0, 1, 0, 1, 0, 2, 3, 2, 3};
  const double ax0[] = {2.0, 1.0, 1.0, 2.0,
                        0.5, 2.0, 1.0, 1.0, 2.0};
  const double ax1[] = {2.0, 1.0, 1.0, 2.0,
                        0.5, 1.0e-12, 1.0, 1.0, 2.0};
  const double b[] = {5.5, 5.0, 4.000000000003, 11.0};
  double x[4] = {0.0, 0.0, 0.0, 0.0};

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_NATURAL;
  options.scale = 2;
  options.pivot_tolerance = 0.001;

  int ok = 1;
  if (!require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 4, ap, ai, 0,
                                        &options),
                        "analyze scaled btf tail continuation")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor scaled btf tail base")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax1),
                        "factor scaled btf tail continuation")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve scaled btf tail continuation")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats scaled btf tail continuation")) {
    ok = 0;
  }
  if (ok && stats.nblocks < 2) {
    fprintf(stderr,
            "scaled btf tail test did not form multiple blocks: %" PRId64 "\n",
            stats.nblocks);
    ok = 0;
  }
  if (ok && (stats.fast_rejected_pivot < stats.fast_rejected_block_start ||
             stats.fast_rejected_pivot >=
               stats.fast_rejected_block_start +
                 stats.fast_rejected_block_size ||
             stats.fast_rejected_block_start +
               stats.fast_rejected_block_size >= stats.n ||
             stats.fast_rejected_refresh_state !=
               KLS_FAST_REJECT_REFRESH_PREFIX ||
             stats.fast_block_restarts != 1 ||
             stats.fast_tail_restarts != 0)) {
    fprintf(stderr,
            "unexpected scaled btf tail stats: pivot=%" PRId64
            ", col=%" PRId64 ", block=[%" PRId64 ",%" PRId64 ")"
            ", n=%" PRId64 ", refresh=%d, block_restarts=%d"
            ", tail_restarts=%d\n",
            stats.fast_rejected_pivot,
            stats.fast_rejected_pivot_col,
            stats.fast_rejected_block_start,
            stats.fast_rejected_block_start + stats.fast_rejected_block_size,
            stats.n,
            stats.fast_rejected_refresh_state,
            stats.fast_block_restarts,
            stats.fast_tail_restarts);
    ok = 0;
  }
  if (ok && (stats.fast_rejected_block_start < 0 ||
             stats.fast_rejected_block_start +
               stats.fast_rejected_block_size >= stats.n ||
             stats.fast_rejected_block_size != 2 ||
             stats.fast_rejected_suffix_columns != 2 ||
             stats.fast_rejected_pivoting_tail_columns < 1 ||
             stats.fast_rejected_pivoting_tail_columns >
               stats.fast_rejected_suffix_columns)) {
    fprintf(stderr,
            "unexpected scaled btf tail plan stats: start=%" PRId64
            ", size=%" PRId64 ", suffix=%" PRId64
            ", pivoting_tail=%" PRId64 "\n",
            stats.fast_rejected_block_start,
            stats.fast_rejected_block_size,
            stats.fast_rejected_suffix_columns,
            stats.fast_rejected_pivoting_tail_columns);
    ok = 0;
  }
  if (ok && !require_pivoting_tail_plan(&stats,
                                        "scaled btf tail continuation")) {
    ok = 0;
  }
  if (ok && (!close_enough(x[0], 1.0) || !close_enough(x[1], 2.0) ||
             !close_enough(x[2], 3.0) || !close_enough(x[3], 4.0))) {
    fprintf(stderr,
            "unexpected scaled btf tail solution: %.17g %.17g %.17g %.17g\n",
            x[0], x[1], x[2], x[3]);
    ok = 0;
  }

  kls_destroy(solver);
  return ok;
}

static int test_btf_fast_factor_block_restart(void) {
  const int32_t ap[] = {0, 2, 4, 7, 9};
  const int32_t ai[] = {0, 1, 0, 1, 0, 2, 3, 2, 3};
  const double ax0[] = {2.0, 1.0, 1.0, 2.0,
                        0.5, 2.0, 1.0, 1.0, 2.0};
  const double ax1[] = {2.0, 1.0, 1.0, 2.0,
                        0.5, 1.0e-12, 1.0, 1.0, 2.0};
  const double b[] = {5.5, 5.0, 4.000000000003, 11.0};
  double x[4] = {0.0, 0.0, 0.0, 0.0};

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_NATURAL;
  options.scale = -1;
  options.pivot_tolerance = 0.001;

  int ok = 1;
  if (!require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 4, ap, ai, 0,
                                        &options),
                        "analyze btf block restart")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor btf block restart base")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax1),
                        "factor btf block restart repair")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve btf block restart")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats btf block restart")) {
    ok = 0;
  }
  if (ok && stats.nblocks < 2) {
    fprintf(stderr, "btf restart test did not form multiple blocks: %" PRId64 "\n",
            stats.nblocks);
    ok = 0;
  }
  if (ok && stats.fast_block_restarts < 1) {
    fprintf(stderr, "btf restart test did not use block restart\n");
    ok = 0;
  }
  if (ok && stats.fast_rejected_refresh_state !=
              KLS_FAST_REJECT_REFRESH_PREFIX) {
    fprintf(stderr, "unexpected btf reject refresh state: %d\n",
            stats.fast_rejected_refresh_state);
    ok = 0;
  }
  if (ok && (stats.fast_rejected_block_start < 0 ||
             stats.fast_rejected_pivot < stats.fast_rejected_block_start ||
             stats.fast_rejected_pivot >=
               stats.fast_rejected_block_start +
                 stats.fast_rejected_block_size ||
             stats.fast_rejected_block_size != 2 ||
             stats.fast_rejected_suffix_columns != 2 ||
             stats.fast_rejected_descendant_columns < 1 ||
             stats.fast_rejected_descendant_columns >
               stats.fast_rejected_suffix_columns ||
             stats.fast_rejected_etree_columns < 1 ||
             stats.fast_rejected_etree_columns >
               stats.fast_rejected_suffix_columns ||
             stats.fast_rejected_pivoting_tail_columns < 1 ||
             stats.fast_rejected_pivoting_tail_columns >
               stats.fast_rejected_suffix_columns)) {
    fprintf(stderr,
            "unexpected btf reject tail stats: start=%" PRId64
            ", size=%" PRId64 ", suffix=%" PRId64
            ", descendants=%" PRId64 ", etree=%" PRId64
            ", pivoting_tail=%" PRId64 "\n",
            stats.fast_rejected_block_start,
            stats.fast_rejected_block_size,
            stats.fast_rejected_suffix_columns,
            stats.fast_rejected_descendant_columns,
            stats.fast_rejected_etree_columns,
            stats.fast_rejected_pivoting_tail_columns);
    ok = 0;
  }
  if (ok && !require_pivoting_tail_plan(&stats, "btf block restart")) {
    ok = 0;
  }
  if (ok && (!close_enough(x[0], 1.0) || !close_enough(x[1], 2.0) ||
             !close_enough(x[2], 3.0) || !close_enough(x[3], 4.0))) {
    fprintf(stderr,
            "unexpected btf block restart solution: %.17g %.17g %.17g %.17g\n",
            x[0], x[1], x[2], x[3]);
    ok = 0;
  }

  kls_destroy(solver);
  return ok;
}

static int test_btf_prefix_tail_restart_with_offblock(void) {
  const int32_t ap[] = {0, 1, 5, 9, 13};
  const int32_t ai[] = {
    0,
    0, 1, 2, 3,
    0, 1, 2, 3,
    0, 1, 2, 3
  };
  const double ax0[] = {
    4.0,
    0.1, 8.0, 0.01, 0.02,
    0.2, 0.01, 8.0, 0.03,
    0.3, 0.02, 0.03, 8.0
  };
  const double ax1[] = {
    4.0,
    0.1, 8.0, 0.01, 0.02,
    1.2, 0.01, 1.0e-12, 2.0,
    -0.4, 0.02, 0.03, 8.0
  };
  const double expected[] = {1.0, 2.0, 3.0, 4.0};
  double b[4] = {0.0, 0.0, 0.0, 0.0};
  double x[4] = {0.0, 0.0, 0.0, 0.0};
  for (int32_t col = 0; col < 4; ++col) {
    for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
      b[ai[p]] += ax1[p] * expected[col];
    }
  }

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_NATURAL;
  options.scale = -1;
  options.pivot_tolerance = 0.001;
  options.static_pivoting = 0;

  int ok = 1;
  if (!require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 4, ap, ai, 0,
                                        &options),
                        "analyze btf offblock tail")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor btf offblock tail base")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax1),
                        "factor btf offblock tail repair")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve btf offblock tail")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats btf offblock tail")) {
    ok = 0;
  }
  if (ok && (stats.nblocks < 2 ||
             stats.fast_rejected_block_start < 0 ||
             stats.fast_rejected_pivot <= stats.fast_rejected_block_start ||
             stats.fast_rejected_pivot >=
               stats.fast_rejected_block_start +
                 stats.fast_rejected_block_size ||
             stats.fast_rejected_refresh_state !=
               KLS_FAST_REJECT_REFRESH_PREFIX ||
             stats.fast_block_restarts != 1 ||
             stats.fast_tail_restarts != 1 ||
             stats.fast_repaired_tail_restart_ready != 1 ||
             stats.fast_repaired_tail_restart_columns < 1 ||
             stats.fast_repaired_tail_restart_work <= 0.0)) {
    fprintf(stderr,
            "unexpected btf offblock tail stats: nblocks=%" PRId64
            ", pivot=%" PRId64 ", start=%" PRId64 ", size=%" PRId64
            ", refresh=%d, restarts=%d, tail_restarts=%d"
            ", tail_ready=%d, tail_cols=%" PRId64 ", tail_work=%.6g\n",
            stats.nblocks,
            stats.fast_rejected_pivot,
            stats.fast_rejected_block_start,
            stats.fast_rejected_block_size,
            stats.fast_rejected_refresh_state,
            stats.fast_block_restarts,
            stats.fast_tail_restarts,
            stats.fast_repaired_tail_restart_ready,
            stats.fast_repaired_tail_restart_columns,
            stats.fast_repaired_tail_restart_work);
    ok = 0;
  }
  if (ok && !require_pivoting_tail_plan(&stats, "btf offblock tail")) {
    ok = 0;
  }
  double max_solution_error = 0.0;
  for (int32_t i = 0; i < 4; ++i) {
    const double err = fabs(x[i] - expected[i]);
    if (err > max_solution_error) {
      max_solution_error = err;
    }
  }
  if (ok && max_solution_error > 1.0e-8) {
    fprintf(stderr,
            "unexpected btf offblock tail solution: %.17g %.17g %.17g %.17g"
            ", max_err=%.17g\n",
            x[0], x[1], x[2], x[3], max_solution_error);
    ok = 0;
  }

  kls_destroy(solver);
  return ok;
}

static int test_fast_factor_restart_after_prior_pivot(void) {
  const int32_t ap[] = {0, 3, 6, 9};
  const int32_t ai[] = {0, 1, 2, 0, 1, 2, 0, 1, 2};
  const double ax0[] = {
    5.8113162339081946e-4, 4.6413198306709501e-2,
    1.5278393799834244e-12, 8.4541227476277784e5,
    6.7324793225895394e-2, 3.0341635684366865e-3,
    5.7525546837879347e1, 1.0034706005124141e4,
    7.5697097316485318e-9
  };
  const double ax1[] = {
    4.2033293705135975e3, 7.3244855455671168e-9,
    1.8481368028583987e-3, 3.7564831178168104e-12,
    9.9130486746531957e-12, 4.0619351468812578e-11,
    4.7736904569941925e2, 8.9161242869151965e-10,
    1.5076094143845339e2
  };
  const double b[] = {
    5.6354365076118629e3, 1.0019148928990983e-8,
    4.5228467245224425e2
  };
  double x[3] = {0.0, 0.0, 0.0};

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_NATURAL;
  options.use_btf = 0;
  options.scale = -1;
  options.pivot_tolerance = 0.001;

  int ok = 1;
  if (!require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 3, ap, ai, 0,
                                        &options),
                        "analyze prior-pivot restart")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor prior-pivot base")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax1),
                        "factor prior-pivot repair")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve prior-pivot repair")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats prior-pivot restart")) {
    ok = 0;
  }
  if (ok && stats.fast_block_restarts != 1) {
    fprintf(stderr,
            "prior-pivot restart count was %d, rejected=%" PRId64
            ", col=%" PRId64 ", offdiag=%" PRId64 "\n",
            stats.fast_block_restarts, stats.fast_rejected_pivot,
            stats.fast_rejected_pivot_col, stats.offdiag_pivots);
    ok = 0;
  }
  if (ok && (stats.fast_rejected_pivot != 0 ||
             stats.fast_rejected_pivot_col != 0)) {
    fprintf(stderr,
            "unexpected prior-pivot rejected pivot: pivot=%" PRId64
            ", col=%" PRId64 "\n",
            stats.fast_rejected_pivot, stats.fast_rejected_pivot_col);
    ok = 0;
  }
  if (ok && (stats.fast_rejected_block_start != 0 ||
             stats.fast_rejected_block_size != 3 ||
             stats.fast_rejected_suffix_columns != 3 ||
             stats.fast_rejected_descendant_columns < 1 ||
             stats.fast_rejected_descendant_columns >
               stats.fast_rejected_suffix_columns ||
             stats.fast_rejected_etree_columns < 1 ||
             stats.fast_rejected_etree_columns >
               stats.fast_rejected_suffix_columns ||
             stats.fast_rejected_pivoting_tail_columns < 1 ||
             stats.fast_rejected_pivoting_tail_columns >
               stats.fast_rejected_suffix_columns)) {
    fprintf(stderr,
            "unexpected prior-pivot tail stats: start=%" PRId64
            ", size=%" PRId64 ", suffix=%" PRId64
            ", descendants=%" PRId64 ", etree=%" PRId64
            ", pivoting_tail=%" PRId64 "\n",
            stats.fast_rejected_block_start,
            stats.fast_rejected_block_size,
            stats.fast_rejected_suffix_columns,
            stats.fast_rejected_descendant_columns,
            stats.fast_rejected_etree_columns,
            stats.fast_rejected_pivoting_tail_columns);
    ok = 0;
  }
  if (ok && !require_pivoting_tail_plan(&stats, "prior-pivot restart")) {
    ok = 0;
  }
  if (ok && (!close_enough(x[0], 1.0) || !close_enough(x[1], 2.0) ||
             !close_enough(x[2], 3.0))) {
    fprintf(stderr, "unexpected prior-pivot restart solution: %.17g %.17g %.17g\n",
            x[0], x[1], x[2]);
    ok = 0;
  }

  kls_destroy(solver);
  return ok;
}

static int test_checked_row_dense_prefix_scatter_tail_restart(void) {
  const int32_t n = 48;
  const size_t nnz = (size_t)n * (size_t)n;
  int32_t *ap = (int32_t *)malloc(((size_t)n + 1u) * sizeof(*ap));
  int32_t *ai = (int32_t *)malloc(nnz * sizeof(*ai));
  double *ax0 = (double *)malloc(nnz * sizeof(*ax0));
  double *ax1 = (double *)malloc(nnz * sizeof(*ax1));
  double *b = (double *)calloc((size_t)n, sizeof(*b));
  double *x = (double *)calloc((size_t)n, sizeof(*x));
  double *expected = (double *)malloc((size_t)n * sizeof(*expected));
  if (ap == NULL || ai == NULL || ax0 == NULL || ax1 == NULL ||
      b == NULL || x == NULL || expected == NULL) {
    free(ap);
    free(ai);
    free(ax0);
    free(ax1);
    free(b);
    free(x);
    free(expected);
    return 0;
  }

  for (int32_t col = 0; col <= n; ++col) {
    ap[col] = col * n;
  }
  for (int32_t col = 0; col < n; ++col) {
    expected[col] = 1.0 + 0.125 * (double)(col % 7);
    for (int32_t row = 0; row < n; ++row) {
      const size_t p = (size_t)col * (size_t)n + (size_t)row;
      ai[p] = row;
      ax0[p] = row == col
        ? 10.0 + 0.01 * (double)col
        : 0.001 * (1.0 + (double)((row + 3 * col) % 7));
      ax1[p] = ax0[p];
    }
  }
  ax1[0] = 20.0;
  ax1[(size_t)10 * (size_t)n + 11u] = 2.0e3;
  for (int32_t col = 0; col < n; ++col) {
    for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
      b[ai[p]] += ax1[p] * expected[col];
    }
  }

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.threads = 1;
  options.ordering = KLS_ORDERING_NATURAL;
  options.orientation = KLS_ORIENTATION_NORMAL;
  options.use_btf = 0;
  options.scale = -1;
  options.pivot_tolerance = 0.01;
  options.static_pivoting = 0;

  const char *saved_env_value = getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  char *saved_env = saved_env_value != NULL ? strdup(saved_env_value) : NULL;
  const int had_saved_env = saved_env_value != NULL;

  int ok = 1;
  if (had_saved_env && saved_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_CHECKED_ROW_REFACTOR\n");
    ok = 0;
  }
  if (!require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, n, ap, ai, 0,
                                        &options),
                        "analyze dense checked-row prefix")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor dense checked-row base")) ok = 0;
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(solver, ax1),
                        "factor dense checked-row repair")) ok = 0;
  if (had_saved_env && saved_env != NULL) {
    if (setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", saved_env, 1) != 0) {
      perror("restore KLS_ENABLE_CHECKED_ROW_REFACTOR");
      ok = 0;
    }
  } else if (!had_saved_env) {
    if (unsetenv("KLS_ENABLE_CHECKED_ROW_REFACTOR") != 0) {
      perror("unsetenv KLS_ENABLE_CHECKED_ROW_REFACTOR");
      ok = 0;
    }
  }
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve dense checked-row prefix")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats dense checked-row prefix")) {
    ok = 0;
  }
  if (ok && (stats.fast_rejected_pivot != 10 ||
             stats.fast_rejected_refresh_state !=
               KLS_FAST_REJECT_REFRESH_PREFIX ||
             stats.fast_block_restarts != 1 ||
             stats.fast_tail_restarts != 1 ||
             stats.fast_repaired_tail_restart_ready != 1)) {
    fprintf(stderr,
            "unexpected dense checked-row prefix stats: pivot=%" PRId64
            ", refresh=%d, block_restarts=%d, tail_restarts=%d"
            ", tail_ready=%d\n",
            stats.fast_rejected_pivot,
            stats.fast_rejected_refresh_state,
            stats.fast_block_restarts,
            stats.fast_tail_restarts,
            stats.fast_repaired_tail_restart_ready);
    ok = 0;
  }
  double max_solution_error = 0.0;
  for (int32_t i = 0; i < n; ++i) {
    const double err = fabs(x[i] - expected[i]);
    if (err > max_solution_error) {
      max_solution_error = err;
    }
  }
  double max_residual = 0.0;
  double max_rhs = 0.0;
  for (int32_t row = 0; row < n; ++row) {
    double residual = -b[row];
    for (int32_t col = 0; col < n; ++col) {
      residual += ax1[(size_t)col * (size_t)n + (size_t)row] * x[col];
    }
    const double residual_abs = fabs(residual);
    if (residual_abs > max_residual) {
      max_residual = residual_abs;
    }
    const double rhs_abs = fabs(b[row]);
    if (rhs_abs > max_rhs) {
      max_rhs = rhs_abs;
    }
  }
  const double relative_residual =
    max_residual / (max_rhs > 0.0 ? max_rhs : 1.0);
  if (ok && (max_solution_error > 2.0e-1 ||
             relative_residual > 1.0e-5)) {
    fprintf(stderr,
            "unexpected dense checked-row accuracy: max_x_err=%.17g"
            ", rel_resid=%.17g\n",
            max_solution_error, relative_residual);
    ok = 0;
  }

  kls_destroy(solver);
  free(saved_env);
  free(ap);
  free(ai);
  free(ax0);
  free(ax1);
  free(b);
  free(x);
  free(expected);
  return ok;
}

static int test_pre_static_pivoting(void) {
  const int32_t n = 3000;
  int32_t *ap = (int32_t *)malloc(((size_t)n + 1u) * sizeof(*ap));
  int32_t *ai = (int32_t *)malloc((size_t)n * sizeof(*ai));
  double *ax = (double *)malloc((size_t)n * sizeof(*ax));
  double *b = (double *)malloc((size_t)n * sizeof(*b));
  double *x = (double *)calloc((size_t)n, sizeof(*x));
  double *expected = (double *)malloc((size_t)n * sizeof(*expected));
  if (ap == NULL || ai == NULL || ax == NULL || b == NULL ||
      x == NULL || expected == NULL) {
    free(ap);
    free(ai);
    free(ax);
    free(b);
    free(x);
    free(expected);
    return 0;
  }

  for (int32_t col = 0; col <= n; ++col) {
    ap[col] = col;
  }
  for (int32_t col = 0; col < n; ++col) {
    ai[col] = (col + n - 1) % n;
    ax[col] = 1.0;
    expected[col] = 1.0 + (double)(col % 17);
  }
  for (int32_t row = 0; row < n; ++row) {
    b[row] = expected[(row + 1) % n];
  }

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_AUTO;

  int ok = 1;
  if (!require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, n, ap, ai, 0,
                                        &options),
                        "analyze pre-static pivot")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax), "factor pre-static pivot")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve pre-static pivot")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats), "stats pre-static pivot")) {
    ok = 0;
  }
  if (ok && !stats.selected_static_pivoting) {
    fprintf(stderr, "pre-static pivoting was not selected\n");
    ok = 0;
  }
  if (ok && !stats.selected_exact_matching) {
    fprintf(stderr, "pre-static pivoting did not use exact matching\n");
    ok = 0;
  }
  if (ok && stats.selected_spral_matching) {
    fprintf(stderr, "small pre-static pivoting unexpectedly used SPRAL\n");
    ok = 0;
  }
  for (int32_t i = 0; ok && i < n; ++i) {
    if (!close_enough(x[i], expected[i])) {
      fprintf(stderr, "unexpected pre-static solution at %d: %.17g != %.17g\n",
              (int)i, x[i], expected[i]);
      ok = 0;
    }
  }

  kls_destroy(solver);
  free(ap);
  free(ai);
  free(ax);
  free(b);
  free(x);
  free(expected);
  return ok;
}

static int test_pre_static_pivoting_with_scaling(void) {
  const int32_t n = 3000;
  int32_t *ap = (int32_t *)malloc(((size_t)n + 1u) * sizeof(*ap));
  int32_t *ai = (int32_t *)malloc((size_t)n * sizeof(*ai));
  double *ax = (double *)malloc((size_t)n * sizeof(*ax));
  double *b = (double *)calloc((size_t)n, sizeof(*b));
  double *bt = (double *)calloc((size_t)n, sizeof(*bt));
  double *x = (double *)calloc((size_t)n, sizeof(*x));
  double *xt = (double *)calloc((size_t)n, sizeof(*xt));
  double *expected = (double *)malloc((size_t)n * sizeof(*expected));
  double *expected_t = (double *)malloc((size_t)n * sizeof(*expected_t));
  if (ap == NULL || ai == NULL || ax == NULL || b == NULL || bt == NULL ||
      x == NULL || xt == NULL || expected == NULL || expected_t == NULL) {
    free(ap);
    free(ai);
    free(ax);
    free(b);
    free(bt);
    free(x);
    free(xt);
    free(expected);
    free(expected_t);
    return 0;
  }

  for (int32_t col = 0; col <= n; ++col) {
    ap[col] = col;
  }
  for (int32_t col = 0; col < n; ++col) {
    const int32_t row = (col + n - 1) % n;
    const double magnitude = (col % 2 == 0) ? 1.0e-6 : 1.0e6;
    ai[col] = row;
    ax[col] = magnitude;
    expected[col] = 1.0 + (double)(col % 19);
    expected_t[row] = 2.0 + (double)(row % 23);
    b[row] = magnitude * expected[col];
    bt[col] = magnitude * expected_t[row];
  }

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_AUTO;

  int ok = 1;
  if (!require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, n, ap, ai, 0,
                                        &options),
                        "analyze scaled pre-static pivot")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax),
                        "factor scaled pre-static pivot")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve scaled pre-static pivot")) ok = 0;
  if (ok && !require_ok(kls_solve_transpose(solver, 1, bt, 0, xt, 0),
                        "transpose solve scaled pre-static pivot")) ok = 0;
  if (ok) {
    memset(x, 0, (size_t)n * sizeof(*x));
    memset(xt, 0, (size_t)n * sizeof(*xt));
  }
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "repeat solve scaled pre-static pivot")) ok = 0;
  if (ok && !require_ok(kls_solve_transpose(solver, 1, bt, 0, xt, 0),
                        "repeat transpose solve scaled pre-static pivot")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats scaled pre-static pivot")) {
    ok = 0;
  }
  if (ok && !stats.selected_static_pivoting) {
    fprintf(stderr, "scaled pre-static pivoting was not selected\n");
    ok = 0;
  }
  if (ok && !stats.selected_exact_matching) {
    fprintf(stderr, "scaled pre-static pivoting did not use exact matching\n");
    ok = 0;
  }
  if (ok && stats.selected_spral_matching) {
    fprintf(stderr, "small scaled pre-static pivoting unexpectedly used SPRAL\n");
    ok = 0;
  }
  if (ok && stats.selected_scale != -1) {
    fprintf(stderr, "matching equilibration did not switch to no-scale mode: %d\n",
            stats.selected_scale);
    ok = 0;
  }
  for (int32_t i = 0; ok && i < n; ++i) {
    if (!close_enough(x[i], expected[i])) {
      fprintf(stderr, "unexpected scaled pre-static solution at %d: %.17g != %.17g\n",
              (int)i, x[i], expected[i]);
      ok = 0;
    }
    if (ok && !close_enough(xt[i], expected_t[i])) {
      fprintf(stderr, "unexpected scaled pre-static transpose solution at %d: %.17g != %.17g\n",
              (int)i, xt[i], expected_t[i]);
      ok = 0;
    }
  }

  kls_destroy(solver);
  free(ap);
  free(ai);
  free(ax);
  free(b);
  free(bt);
  free(x);
  free(xt);
  free(expected);
  free(expected_t);
  return ok;
}

static int test_parallel_row_refactor_pipeline_scope(void) {
  const int32_t n = 12;
  int32_t ap[13];
  int32_t ai[49];
  double ax0[49];
  double ax1[49];
  double expected[12];
  double b[12] = {0.0};
  double x[12] = {0.0};

  int32_t p = 0;
  for (int32_t col = 0; col < n; ++col) {
    ap[col] = p;
    expected[col] = 1.0 + 0.25 * (double)col;
    const int32_t row_block = col < 4 ? 0 : (col < 8 ? 4 : 8);
    for (int32_t row = row_block; row < row_block + 4; ++row) {
      ai[p] = row;
      ax0[p] = row == col
        ? 12.0 + (double)col
        : 0.05 * (double)(1 + ((row + 2 * col) % 5));
      ax1[p] = ax0[p] + (row == col ? 0.125 : 0.01);
      p++;
    }
    if (col == 0) {
      ai[p] = 4;
      ax0[p] = 0.035;
      ax1[p] = 0.041;
      p++;
    }
  }
  ap[n] = p;
  for (int32_t col = 0; col < n; ++col) {
    for (int32_t q = ap[col]; q < ap[col + 1]; ++q) {
      b[ai[q]] += ax1[q] * expected[col];
    }
  }

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.threads = 4;
  options.orientation = KLS_ORIENTATION_NORMAL;
  options.ordering = KLS_ORDERING_NATURAL;
  options.use_btf = 0;
  options.scale = -1;
  options.static_pivoting = 0;

  const char *saved_env_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_env = saved_env_value != NULL ? strdup(saved_env_value) : NULL;
  const int had_saved_env = saved_env_value != NULL;
  const char *saved_checked_env_value =
    getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  char *saved_checked_env = saved_checked_env_value != NULL
    ? strdup(saved_checked_env_value) : NULL;
  const int had_saved_checked_env = saved_checked_env_value != NULL;

  int ok = 1;
  if (had_saved_env && saved_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_ROW_REFACTOR\n");
    ok = 0;
  }
  if (had_saved_checked_env && saved_checked_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_CHECKED_ROW_REFACTOR\n");
    ok = 0;
  }
  if (!require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, n, ap, ai, 0,
                                        &options),
                        "analyze row-pipeline refactor")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor row-pipeline base")) ok = 0;
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(solver, ax1),
                        "checked row-pipeline factor")) ok = 0;
  kls_stats checked_stats;
  checked_stats.struct_size = sizeof(checked_stats);
  if (ok && !require_ok(kls_get_stats(solver, &checked_stats),
                        "stats checked row-pipeline factor")) {
    ok = 0;
  }
  const int checked_expect_defer =
    checked_stats.row_refactor_dense_segment_count > 0;
  if (ok && (checked_stats.row_refactor_last_run != 1 ||
             checked_stats.row_refactor_last_checked != 1 ||
             checked_stats.row_refactor_last_parallel != 1 ||
             checked_stats.row_refactor_last_ready_queue != 1 ||
             checked_stats.row_refactor_last_done_bitmap != 1 ||
             checked_stats.row_refactor_last_defer_value_scatter !=
               checked_expect_defer ||
             checked_stats.row_refactor_last_work_ready_queue != 1 ||
             checked_stats.row_refactor_run_count != 1 ||
             checked_stats.row_refactor_checked_run_count != 1 ||
             checked_stats.row_refactor_parallel_run_count != 1 ||
             checked_stats.row_refactor_ready_queue_run_count != 1 ||
             checked_stats.row_refactor_ready_queue_group_count < 3 ||
             checked_stats.row_refactor_done_bitmap_run_count != 1 ||
             checked_stats.row_refactor_defer_value_scatter_run_count !=
               checked_expect_defer ||
             checked_stats.row_refactor_work_ready_queue_run_count != 1)) {
    fprintf(stderr,
            "unexpected checked row stats: last=%d/%d/%d/%d/%d/%d/%d"
            ", runs=%" PRId64 "/%" PRId64 "/%" PRId64
            ", ready=%" PRId64 "/%" PRId64
            ", done_bitmap=%" PRId64
            ", defer_scatter=%" PRId64
            ", work_queue=%" PRId64 "\n",
            checked_stats.row_refactor_last_run,
            checked_stats.row_refactor_last_checked,
            checked_stats.row_refactor_last_parallel,
            checked_stats.row_refactor_last_ready_queue,
            checked_stats.row_refactor_last_done_bitmap,
            checked_stats.row_refactor_last_defer_value_scatter,
            checked_stats.row_refactor_last_work_ready_queue,
            checked_stats.row_refactor_run_count,
            checked_stats.row_refactor_checked_run_count,
            checked_stats.row_refactor_parallel_run_count,
            checked_stats.row_refactor_ready_queue_run_count,
            checked_stats.row_refactor_ready_queue_group_count,
            checked_stats.row_refactor_done_bitmap_run_count,
            checked_stats.row_refactor_defer_value_scatter_run_count,
            checked_stats.row_refactor_work_ready_queue_run_count);
    ok = 0;
  }
  if (ok && !require_ok(kls_refactor(solver, ax1),
                        "checked-only default refactor")) ok = 0;
  kls_stats checked_only_refactor_stats;
  checked_only_refactor_stats.struct_size =
    sizeof(checked_only_refactor_stats);
  if (ok && !require_ok(kls_get_stats(solver, &checked_only_refactor_stats),
                        "stats checked-only default refactor")) {
    ok = 0;
  }
  if (ok && (checked_only_refactor_stats.row_refactor_last_run != 0 ||
             checked_only_refactor_stats.row_refactor_last_ready_queue != 0 ||
             checked_only_refactor_stats.row_refactor_last_done_bitmap != 0 ||
             checked_only_refactor_stats.row_refactor_last_defer_value_scatter != 0 ||
             checked_only_refactor_stats.row_refactor_last_work_ready_queue != 0 ||
             checked_only_refactor_stats.row_refactor_run_count != 1 ||
             checked_only_refactor_stats.row_refactor_checked_run_count != 1)) {
    fprintf(stderr,
            "checked-only refactor incorrectly used row path: last=%d/%d/%d/%d/%d"
            ", runs=%" PRId64 ", checked=%" PRId64 "\n",
            checked_only_refactor_stats.row_refactor_last_run,
            checked_only_refactor_stats.row_refactor_last_ready_queue,
            checked_only_refactor_stats.row_refactor_last_done_bitmap,
            checked_only_refactor_stats.row_refactor_last_defer_value_scatter,
            checked_only_refactor_stats.row_refactor_last_work_ready_queue,
            checked_only_refactor_stats.row_refactor_run_count,
            checked_only_refactor_stats.row_refactor_checked_run_count);
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && !require_ok(kls_refactor(solver, ax1),
                        "row-pipeline refactor")) ok = 0;
  if (had_saved_env && saved_env != NULL) {
    if (setenv("KLS_ENABLE_ROW_REFACTOR", saved_env, 1) != 0) {
      perror("restore KLS_ENABLE_ROW_REFACTOR");
      ok = 0;
    }
  } else if (!had_saved_env) {
    if (unsetenv("KLS_ENABLE_ROW_REFACTOR") != 0) {
      perror("unsetenv KLS_ENABLE_ROW_REFACTOR");
      ok = 0;
    }
  }
  if (had_saved_checked_env && saved_checked_env != NULL) {
    if (setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR",
               saved_checked_env, 1) != 0) {
      perror("restore KLS_ENABLE_CHECKED_ROW_REFACTOR");
      ok = 0;
    }
  } else if (!had_saved_checked_env) {
    if (unsetenv("KLS_ENABLE_CHECKED_ROW_REFACTOR") != 0) {
      perror("unsetenv KLS_ENABLE_CHECKED_ROW_REFACTOR");
      ok = 0;
    }
  }
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve row-pipeline refactor")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats row-pipeline refactor")) {
    ok = 0;
  }
  const int expect_defer =
    stats.row_refactor_dense_segment_count > 0;
  if (ok && (stats.row_refactor_group_count < 2 ||
             stats.row_refactor_group_level_count < 1 ||
             stats.row_refactor_group_cluster_levels != 0 ||
             stats.row_refactor_group_pipeline_groups < 2 ||
             stats.row_refactor_group_pipeline_rows != n ||
             stats.row_refactor_group_pipeline_work <= 0.0 ||
             stats.row_refactor_group_dependency_edges <= 0 ||
             stats.row_refactor_group_root_count <= 0 ||
             stats.row_refactor_group_root_count >
               stats.row_refactor_group_count ||
             stats.row_refactor_group_leaf_count <= 0 ||
             stats.row_refactor_group_leaf_count >
               stats.row_refactor_group_count ||
             (stats.row_refactor_group_dependency_edges == 0 &&
              stats.row_refactor_group_max_fanout != 0) ||
             (stats.row_refactor_group_dependency_edges > 0 &&
              stats.row_refactor_group_max_fanout <= 0) ||
             stats.row_refactor_last_run != 1 ||
             stats.row_refactor_last_checked != 0 ||
             stats.row_refactor_last_parallel != 1 ||
             stats.row_refactor_last_ready_queue != 1 ||
             stats.row_refactor_last_done_bitmap != 0 ||
             stats.row_refactor_last_defer_value_scatter != 0 ||
             stats.row_refactor_last_work_ready_queue != 1 ||
             stats.row_refactor_run_count != 2 ||
             stats.row_refactor_checked_run_count != 1 ||
             stats.row_refactor_parallel_run_count != 2 ||
             stats.row_refactor_ready_queue_run_count != 2 ||
             stats.row_refactor_ready_queue_group_count <
               stats.row_refactor_group_pipeline_groups ||
             stats.row_refactor_done_bitmap_run_count != 1 ||
             stats.row_refactor_input_cleanup_rows != 0 ||
             stats.row_refactor_input_cleanup_entries != 0 ||
             stats.row_refactor_defer_value_scatter_run_count !=
               expect_defer ||
             stats.row_refactor_work_ready_queue_run_count != 2)) {
    fprintf(stderr,
            "unexpected row pipeline stats: groups=%" PRId64
            ", levels=%" PRId64 ", cluster=%" PRId64
            ", pipe_groups=%" PRId64 ", pipe_rows=%" PRId64
            ", pipe_work=%.6g, edges=%" PRId64
            ", roots=%" PRId64 ", leaves=%" PRId64
            ", max_fanout=%" PRId64 ", last=%d/%d/%d/%d/%d/%d/%d"
            ", runs=%" PRId64 "/%" PRId64 "/%" PRId64
            ", ready=%" PRId64 "/%" PRId64
            ", done_bitmap=%" PRId64
            ", cleanup=%" PRId64 "/%" PRId64
            ", defer_scatter=%" PRId64
            ", work_queue=%" PRId64 "\n",
            stats.row_refactor_group_count,
            stats.row_refactor_group_level_count,
            stats.row_refactor_group_cluster_levels,
            stats.row_refactor_group_pipeline_groups,
            stats.row_refactor_group_pipeline_rows,
            stats.row_refactor_group_pipeline_work,
            stats.row_refactor_group_dependency_edges,
            stats.row_refactor_group_root_count,
            stats.row_refactor_group_leaf_count,
            stats.row_refactor_group_max_fanout,
            stats.row_refactor_last_run,
            stats.row_refactor_last_checked,
            stats.row_refactor_last_parallel,
            stats.row_refactor_last_ready_queue,
            stats.row_refactor_last_done_bitmap,
            stats.row_refactor_last_defer_value_scatter,
            stats.row_refactor_last_work_ready_queue,
            stats.row_refactor_run_count,
            stats.row_refactor_checked_run_count,
            stats.row_refactor_parallel_run_count,
            stats.row_refactor_ready_queue_run_count,
            stats.row_refactor_ready_queue_group_count,
            stats.row_refactor_done_bitmap_run_count,
            stats.row_refactor_input_cleanup_rows,
            stats.row_refactor_input_cleanup_entries,
            stats.row_refactor_defer_value_scatter_run_count,
            stats.row_refactor_work_ready_queue_run_count);
    ok = 0;
  }
  const size_t legacy_stats_size =
    offsetof(kls_stats, row_refactor_group_cluster_levels);
  if (ok && legacy_stats_size <
              offsetof(kls_stats, row_refactor_dense_segment_trailing_entries) +
                sizeof(stats.row_refactor_dense_segment_trailing_entries)) {
    fprintf(stderr, "row pipeline stats are not append-only\n");
    ok = 0;
  }
  kls_stats legacy_stats;
  memset(&legacy_stats, 0, sizeof(legacy_stats));
  legacy_stats.struct_size = legacy_stats_size;
  if (ok && !require_ok(kls_get_stats(solver, &legacy_stats),
                        "legacy-size stats row-pipeline refactor")) {
    ok = 0;
  }
  if (ok && (legacy_stats.row_refactor_group_count !=
               stats.row_refactor_group_count ||
             legacy_stats.row_refactor_group_level_count !=
               stats.row_refactor_group_level_count ||
             legacy_stats.row_refactor_group_cluster_levels != 0)) {
    fprintf(stderr, "legacy-size row pipeline stats copy shifted fields\n");
    ok = 0;
  }
  for (int32_t i = 0; ok && i < n; ++i) {
    if (!close_enough(x[i], expected[i])) {
      fprintf(stderr, "unexpected row-pipeline solution at %d: %.17g != %.17g\n",
              (int)i, x[i], expected[i]);
      ok = 0;
    }
  }

  kls_destroy(solver);
  free(saved_env);
  free(saved_checked_env);
  return ok;
}

int main(void) {
  if (!test_csc()) {
    return EXIT_FAILURE;
  }
  if (!test_csr_and_refactor()) {
    return EXIT_FAILURE;
  }
  if (!test_sparse_diagonal_auto_scale()) {
    return EXIT_FAILURE;
  }
  if (!test_forced_transpose_orientation()) {
    return EXIT_FAILURE;
  }
  if (!test_csr_forced_transpose_orientation()) {
    return EXIT_FAILURE;
  }
  if (!test_solve_strides_and_in_place()) {
    return EXIT_FAILURE;
  }
  if (!test_fast_factor_pivot_check_fallback()) {
    return EXIT_FAILURE;
  }
  if (!test_checked_row_fast_factor_block_restart()) {
    return EXIT_FAILURE;
  }
  if (!test_parallel_checked_row_fast_factor_block_restart()) {
    return EXIT_FAILURE;
  }
  if (!test_fast_factor_tail_prefix_state_validation()) {
    return EXIT_FAILURE;
  }
  if (!test_mapped_fast_factor_prefix_tail_restart()) {
    return EXIT_FAILURE;
  }
  if (!test_scaled_fast_factor_block_restart()) {
    return EXIT_FAILURE;
  }
  if (!test_scaled_fast_factor_prefix_tail_restart()) {
    return EXIT_FAILURE;
  }
  if (!test_scaled_btf_fast_factor_tail_continuation()) {
    return EXIT_FAILURE;
  }
  if (!test_btf_fast_factor_block_restart()) {
    return EXIT_FAILURE;
  }
  if (!test_btf_prefix_tail_restart_with_offblock()) {
    return EXIT_FAILURE;
  }
  if (!test_fast_factor_restart_after_prior_pivot()) {
    return EXIT_FAILURE;
  }
  if (!test_checked_row_dense_prefix_scatter_tail_restart()) {
    return EXIT_FAILURE;
  }
  if (!test_pre_static_pivoting()) {
    return EXIT_FAILURE;
  }
  if (!test_pre_static_pivoting_with_scaling()) {
    return EXIT_FAILURE;
  }
  if (!test_parallel_row_refactor_pipeline_scope()) {
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
