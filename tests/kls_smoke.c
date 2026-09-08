#define _GNU_SOURCE
#define _POSIX_C_SOURCE 200809L

#include "kls/kls.h"

#include <math.h>
#include <inttypes.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __linux__
#include <sched.h>
static cpu_set_t caller_affinity;
static int caller_affinity_valid;
#endif

static int restore_env_value(const char *name, int had_value,
                             const char *saved_value);

static int require_ok(int status, const char *what) {
#ifdef __linux__
  /* Worker placement may temporarily pin the caller, but no public call
     may leave its affinity changed, including calls returning an error. */
  if (caller_affinity_valid) {
    cpu_set_t current;
    if (sched_getaffinity(0, sizeof(current), &current) != 0 ||
        !CPU_EQUAL(&current, &caller_affinity)) {
      fprintf(stderr, "%s did not preserve caller affinity\n", what);
      return 0;
    }
  }
#endif
  if (status != KLS_OK) {
    fprintf(stderr, "%s failed: %s (%d)\n", what, kls_status_string(status), status);
    return 0;
  }
  return 1;
}

static int close_enough(double a, double b) {
  return fabs(a - b) < 1e-10;
}

static int require_build_feature_stats(const kls_stats *stats,
                                       const char *what) {
  if (stats == NULL ||
      (stats->build_has_metis != 0 && stats->build_has_metis != 1) ||
      (stats->build_has_scotch != 0 && stats->build_has_scotch != 1) ||
      (stats->build_has_spral_scaling != 0 &&
       stats->build_has_spral_scaling != 1) ||
      (stats->selected_spral_matching &&
       !stats->build_has_spral_scaling)) {
    fprintf(stderr,
            "unexpected build feature stats for %s: metis=%d scotch=%d "
            "spral=%d selected_spral=%d\n",
            what,
            stats != NULL ? stats->build_has_metis : -1,
            stats != NULL ? stats->build_has_scotch : -1,
            stats != NULL ? stats->build_has_spral_scaling : -1,
            stats != NULL ? stats->selected_spral_matching : -1);
    return 0;
  }
  return 1;
}

static int require_parallel_model_stats(const kls_stats *stats,
                                        const char *what) {
  if (stats != NULL &&
      (stats->refactor_last_egraph_algorithm5_prefactor_columns != 0 ||
       stats->refactor_last_egraph_algorithm5_prefactor_deps != 0 ||
       stats->refactor_egraph_algorithm5_prefactor_column_count != 0 ||
       stats->refactor_egraph_algorithm5_prefactor_dep_count != 0)) {
    fprintf(stderr, "retired prefactor stats are nonzero for %s\n", what);
    return 0;
  }
  if (stats == NULL ||
      stats->parallel_model_r1 <= 0.0 ||
      stats->parallel_model_r2 <= 0.0) {
    fprintf(stderr,
            "missing parallel model stats for %s: r1=%.17g r2=%.17g\n",
            what,
            stats != NULL ? stats->parallel_model_r1 : 0.0,
            stats != NULL ? stats->parallel_model_r2 : 0.0);
    return 0;
  }

  const int expected = (stats->parallel_model_r1 >= 2.0 ||
                        stats->parallel_model_r2 >= 50.0);
  if (stats->parallel_model_recommends_parallel != expected) {
    fprintf(stderr,
            "parallel model recommendation mismatch for %s: "
            "r1=%.17g r2=%.17g got=%d expected=%d\n",
            what,
            stats->parallel_model_r1,
            stats->parallel_model_r2,
            stats->parallel_model_recommends_parallel,
            expected);
    return 0;
  }
  if (stats->parallel_task_flow_threads <= 0 ||
      stats->parallel_task_flow_work <= 0.0 ||
      stats->parallel_task_flow_finish_time <= 0.0 ||
      stats->parallel_task_flow_speedup <= 0.0 ||
      stats->parallel_task_flow_dependencies < 0) {
    fprintf(stderr,
            "missing NICSLU task-flow stats for %s: threads=%" PRId64
            " deps=%" PRId64 " work=%.17g finish=%.17g speedup=%.17g\n",
            what,
            stats->parallel_task_flow_threads,
            stats->parallel_task_flow_dependencies,
            stats->parallel_task_flow_work,
            stats->parallel_task_flow_finish_time,
            stats->parallel_task_flow_speedup);
    return 0;
  }
  const int expected_task =
    stats->parallel_task_flow_threads > 1 &&
    stats->parallel_task_flow_speedup > 1.0;
  if (stats->parallel_task_flow_recommends_parallel != expected_task) {
    fprintf(stderr,
            "NICSLU task-flow recommendation mismatch for %s: "
            "threads=%" PRId64 " speedup=%.17g got=%d expected=%d\n",
            what,
            stats->parallel_task_flow_threads,
            stats->parallel_task_flow_speedup,
            stats->parallel_task_flow_recommends_parallel,
            expected_task);
    return 0;
  }
  return 1;
}

static int require_row_auto_model_stats(const kls_stats *stats,
                                        const char *what) {
  if (stats == NULL ||
      (stats->row_refactor_auto_model_recommended != 0 &&
       stats->row_refactor_auto_model_recommended != 1) ||
      (stats->row_refactor_auto_model_attempted != 0 &&
       stats->row_refactor_auto_model_attempted != 1) ||
      (stats->row_refactor_auto_model_accepted != 0 &&
       stats->row_refactor_auto_model_accepted != 1) ||
      (stats->row_refactor_auto_model_attempted &&
       !stats->row_refactor_auto_model_recommended) ||
      (stats->row_refactor_auto_model_accepted &&
       !stats->row_refactor_auto_model_attempted) ||
      !isfinite(stats->row_refactor_auto_lower_bound_work) ||
      stats->row_refactor_auto_lower_bound_work < 0.0 ||
      (stats->row_refactor_auto_lower_bound_rejected != 0 &&
       stats->row_refactor_auto_lower_bound_rejected != 1) ||
      (stats->row_refactor_auto_pattern_build_failed != 0 &&
       stats->row_refactor_auto_pattern_build_failed != 1) ||
      (stats->row_refactor_auto_value_copy_failed != 0 &&
       stats->row_refactor_auto_value_copy_failed != 1)) {
    fprintf(stderr,
            "unexpected row auto model stats for %s: rec=%d attempted=%d "
            "accepted=%d lower_bound=%.17g lower_reject=%d "
            "pattern_failed=%d value_failed=%d\n",
            what,
            stats != NULL ? stats->row_refactor_auto_model_recommended : -1,
            stats != NULL ? stats->row_refactor_auto_model_attempted : -1,
            stats != NULL ? stats->row_refactor_auto_model_accepted : -1,
            stats != NULL ? stats->row_refactor_auto_lower_bound_work : -1.0,
            stats != NULL ? stats->row_refactor_auto_lower_bound_rejected : -1,
            stats != NULL ? stats->row_refactor_auto_pattern_build_failed : -1,
            stats != NULL ? stats->row_refactor_auto_value_copy_failed : -1);
    return 0;
  }
  return 1;
}

static int require_pivoting_tail_plan(const kls_stats *stats,
                                      const char *what) {
  if (stats == NULL || stats->fast_rejected_pivot < 0 ||
      stats->fast_rejected_block_start < 0 ||
      stats->fast_rejected_block_size <= 0 ||
      stats->fast_rejected_pivoting_tail_columns <= 0 ||
      stats->fast_rejected_pivoting_tail_seed_columns <= 0 ||
      stats->fast_rejected_pivoting_tail_seed_columns >
        stats->fast_rejected_pivoting_tail_columns ||
      stats->fast_rejected_pivoting_tail_row_seed_columns < 0 ||
      stats->fast_rejected_pivoting_tail_row_seed_columns >
        stats->fast_rejected_pivoting_tail_seed_columns ||
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
      !stats->fast_rejected_pivoting_tail_topological ||
      stats->fast_rejected_pivoting_tail_gap_columns < 0 ||
      stats->fast_rejected_pivoting_tail_suffix_overcompute_columns < 0 ||
      stats->fast_rejected_pivoting_tail_suffix_overcompute_work < 0.0 ||
      stats->fast_rejected_pivoting_tail_etree_edges < 0 ||
      stats->fast_rejected_pivoting_tail_etree_roots <= 0 ||
      stats->fast_rejected_pivoting_tail_etree_leaves <= 0 ||
      stats->fast_rejected_pivoting_tail_etree_max_fanout < 0 ||
      stats->fast_rejected_pivoting_tail_etree_levels <= 0 ||
      stats->fast_rejected_pivoting_tail_etree_levels >
        stats->fast_rejected_pivoting_tail_columns ||
      stats->fast_rejected_pivoting_tail_etree_max_width <= 0 ||
      stats->fast_rejected_pivoting_tail_etree_max_width >
        stats->fast_rejected_pivoting_tail_columns ||
      stats->fast_rejected_pivoting_tail_etree_max_width <
        stats->fast_rejected_pivoting_tail_etree_leaves ||
      stats->fast_rejected_pivoting_tail_etree_edges +
        stats->fast_rejected_pivoting_tail_etree_roots !=
        stats->fast_rejected_pivoting_tail_columns ||
      stats->fast_rejected_pivoting_tail_etree_leaves >
        stats->fast_rejected_pivoting_tail_columns ||
      (stats->fast_rejected_pivoting_tail_etree_edges > 0 &&
       stats->fast_rejected_pivoting_tail_etree_max_fanout <= 0) ||
      (stats->fast_rejected_pivoting_tail_suffix_exact &&
       (!stats->fast_rejected_pivoting_tail_contiguous ||
        stats->fast_rejected_pivoting_tail_gap_columns != 0 ||
        stats->fast_rejected_pivoting_tail_suffix_overcompute_columns != 0 ||
        stats->fast_rejected_pivoting_tail_suffix_overcompute_work != 0.0)) ||
      (!stats->fast_rejected_pivoting_tail_contiguous &&
       stats->fast_rejected_pivoting_tail_suffix_exact)) {
    fprintf(stderr,
            "unexpected pivoting tail plan for %s: pivot=%" PRId64
            ", block=[%" PRId64 ",%" PRId64 "), cols=%" PRId64
            ", seed=%" PRId64 ", row_seed=%" PRId64
            ", first=%" PRId64 ", last=%" PRId64
            ", contains=%d, topo=%d, contiguous=%d, suffix_exact=%d"
            ", gaps=%" PRId64 ", suffix_over_cols=%" PRId64
            ", suffix_over_work=%.6g"
            ", etree_edges=%" PRId64 ", etree_roots=%" PRId64
            ", etree_leaves=%" PRId64 ", etree_max_fanout=%" PRId64
            ", etree_levels=%" PRId64 ", etree_max_width=%" PRId64 "\n",
            what,
            stats != NULL ? stats->fast_rejected_pivot : -1,
            stats != NULL ? stats->fast_rejected_block_start : -1,
            stats != NULL
              ? stats->fast_rejected_block_start +
                  stats->fast_rejected_block_size
              : -1,
            stats != NULL ? stats->fast_rejected_pivoting_tail_columns : 0,
            stats != NULL
              ? stats->fast_rejected_pivoting_tail_seed_columns
              : 0,
            stats != NULL
              ? stats->fast_rejected_pivoting_tail_row_seed_columns
              : 0,
            stats != NULL ? stats->fast_rejected_pivoting_tail_first : -1,
            stats != NULL ? stats->fast_rejected_pivoting_tail_last : -1,
            stats != NULL
              ? stats->fast_rejected_pivoting_tail_contains_reject
              : 0,
            stats != NULL
              ? stats->fast_rejected_pivoting_tail_topological
              : 0,
            stats != NULL
              ? stats->fast_rejected_pivoting_tail_contiguous
              : 0,
            stats != NULL
              ? stats->fast_rejected_pivoting_tail_suffix_exact
              : 0,
            stats != NULL
              ? stats->fast_rejected_pivoting_tail_gap_columns
              : 0,
            stats != NULL
              ? stats->fast_rejected_pivoting_tail_suffix_overcompute_columns
              : 0,
            stats != NULL
              ? stats->fast_rejected_pivoting_tail_suffix_overcompute_work
              : 0.0,
            stats != NULL
              ? stats->fast_rejected_pivoting_tail_etree_edges
              : 0,
            stats != NULL
              ? stats->fast_rejected_pivoting_tail_etree_roots
              : 0,
            stats != NULL
              ? stats->fast_rejected_pivoting_tail_etree_leaves
              : 0,
            stats != NULL
              ? stats->fast_rejected_pivoting_tail_etree_max_fanout
              : 0,
            stats != NULL
              ? stats->fast_rejected_pivoting_tail_etree_levels
              : 0,
            stats != NULL
              ? stats->fast_rejected_pivoting_tail_etree_max_width
              : 0);
    return 0;
  }
  if ((stats->fast_kls_block_restart_last_row_pipeline_separator_tail_scope !=
         0 &&
       stats->fast_kls_block_restart_last_row_pipeline_separator_tail_scope !=
         1) ||
      stats->fast_kls_block_restart_row_pipeline_separator_tail_scope_count <
        0 ||
      stats->fast_kls_block_restart_last_row_pipeline_separator_tail_scope_rows <
        0 ||
      (stats->fast_kls_block_restart_last_row_pipeline_separator_tail_scope &&
       (stats->fast_kls_block_restart_row_pipeline_separator_tail_scope_count <=
          0 ||
        stats->fast_kls_block_restart_last_row_pipeline_separator_tail_scope_rows <=
          0 ||
        stats->fast_kls_block_restart_last_row_pipeline != 1 ||
        stats->fast_kls_block_restart_last_row_pipeline_etree_tail != 1 ||
        stats->fast_kls_block_restart_last_row_pipeline_separator_tail_scope_rows >
          stats->fast_kls_block_restart_last_row_pipeline_rows)) ||
      (!stats->fast_kls_block_restart_last_row_pipeline_separator_tail_scope &&
       (stats->fast_kls_block_restart_row_pipeline_separator_tail_scope_count !=
          0 ||
        stats->fast_kls_block_restart_last_row_pipeline_separator_tail_scope_rows !=
          0))) {
    fprintf(stderr,
            "unexpected separator tail scope stats for %s: scope=%d/%" PRId64
            ", rows=%" PRId64 ", pipeline=%d, etree_tail=%d"
            ", pipeline_rows=%" PRId64 "\n",
            what,
            stats->fast_kls_block_restart_last_row_pipeline_separator_tail_scope,
            stats->fast_kls_block_restart_row_pipeline_separator_tail_scope_count,
            stats->fast_kls_block_restart_last_row_pipeline_separator_tail_scope_rows,
            stats->fast_kls_block_restart_last_row_pipeline,
            stats->fast_kls_block_restart_last_row_pipeline_etree_tail,
            stats->fast_kls_block_restart_last_row_pipeline_rows);
    return 0;
  }
  if ((stats->fast_kls_block_restart_last_row_pipeline_separator_queue != 0 &&
       stats->fast_kls_block_restart_last_row_pipeline_separator_queue != 1) ||
      stats->fast_kls_block_restart_row_pipeline_separator_queue_count < 0 ||
      stats->fast_kls_block_restart_last_row_pipeline_separator_private_rows <
        0 ||
      stats->fast_kls_block_restart_last_row_pipeline_separator_pipeline_rows <
        0 ||
      stats->fast_kls_block_restart_last_row_pipeline_separator_private_threads <
        0 ||
      (stats
         ->fast_kls_block_restart_last_row_pipeline_separator_partitioned !=
         0 &&
       stats
         ->fast_kls_block_restart_last_row_pipeline_separator_partitioned !=
         1) ||
      stats
        ->fast_kls_block_restart_last_row_pipeline_separator_split_components <
        0 ||
      (stats->fast_kls_block_restart_last_row_pipeline_separator_queue &&
       (stats->fast_kls_block_restart_row_pipeline_separator_queue_count <= 0 ||
        stats
          ->fast_kls_block_restart_last_row_pipeline_separator_private_rows <=
          0 ||
        stats
          ->fast_kls_block_restart_last_row_pipeline_separator_pipeline_rows <=
          0 ||
        stats
          ->fast_kls_block_restart_last_row_pipeline_separator_private_threads <=
          0 ||
        stats->fast_kls_block_restart_last_row_pipeline != 1 ||
        stats
            ->fast_kls_block_restart_last_row_pipeline_separator_private_rows +
          stats
            ->fast_kls_block_restart_last_row_pipeline_separator_pipeline_rows !=
          stats->fast_kls_block_restart_last_row_pipeline_rows ||
        (stats
           ->fast_kls_block_restart_last_row_pipeline_separator_partitioned &&
         stats
           ->fast_kls_block_restart_last_row_pipeline_separator_split_components <=
           0))) ||
      (!stats->fast_kls_block_restart_last_row_pipeline_separator_queue &&
       (stats->fast_kls_block_restart_row_pipeline_separator_queue_count != 0 ||
        stats
          ->fast_kls_block_restart_last_row_pipeline_separator_private_rows !=
          0 ||
        stats
          ->fast_kls_block_restart_last_row_pipeline_separator_pipeline_rows !=
          0 ||
        stats
          ->fast_kls_block_restart_last_row_pipeline_separator_private_threads !=
          0 ||
        stats
          ->fast_kls_block_restart_last_row_pipeline_separator_partitioned !=
          0 ||
        stats
          ->fast_kls_block_restart_last_row_pipeline_separator_split_components !=
          0))) {
    fprintf(stderr,
            "unexpected separator queue repair stats for %s: queue=%d/%" PRId64
            ", rows=%" PRId64 "/%" PRId64 ", threads=%" PRId64
            ", partitioned=%d, split=%" PRId64 ", pipeline=%d"
            ", pipeline_rows=%" PRId64 "\n",
            what,
            stats->fast_kls_block_restart_last_row_pipeline_separator_queue,
            stats
              ->fast_kls_block_restart_row_pipeline_separator_queue_count,
            stats
              ->fast_kls_block_restart_last_row_pipeline_separator_private_rows,
            stats
              ->fast_kls_block_restart_last_row_pipeline_separator_pipeline_rows,
            stats
              ->fast_kls_block_restart_last_row_pipeline_separator_private_threads,
            stats
              ->fast_kls_block_restart_last_row_pipeline_separator_partitioned,
            stats
              ->fast_kls_block_restart_last_row_pipeline_separator_split_components,
            stats->fast_kls_block_restart_last_row_pipeline,
            stats->fast_kls_block_restart_last_row_pipeline_rows);
    return 0;
  }
  return 1;
}

static int require_tail_overcompute_bounds(const kls_stats *stats,
                                           const char *what) {
  if (stats == NULL ||
      stats->fast_repaired_tail_restart_overcompute_columns < 0 ||
      stats->fast_repaired_tail_restart_overcompute_work < 0.0 ||
      stats->fast_repaired_tail_restart_skipped_columns < 0 ||
      stats->fast_repaired_tail_restart_skipped_work < 0.0 ||
      stats->fast_kls_rebuild_restarts < 0 ||
      (stats->fast_repaired_tail_restart_exact_mask != 0 &&
       stats->fast_repaired_tail_restart_exact_mask != 1) ||
      (stats->fast_repaired_tail_restart_etree_mask != 0 &&
       stats->fast_repaired_tail_restart_etree_mask != 1) ||
      (stats->fast_repaired_tail_restart_etree_mask &&
       !stats->fast_repaired_tail_restart_exact_mask) ||
      (stats->fast_repaired_tail_restart_exact_mask &&
       (!stats->fast_repaired_tail_restart_ready ||
        stats->fast_repaired_tail_restart_overcompute_columns != 0)) ||
      stats->fast_repaired_tail_restart_overcompute_columns >
        stats->fast_repaired_tail_restart_columns ||
      stats->fast_repaired_tail_restart_overcompute_work >
        stats->fast_repaired_tail_restart_work + 1.0e-9) {
    fprintf(stderr,
            "unexpected tail overcompute stats for %s: tail_cols=%" PRId64
            ", over_cols=%" PRId64 ", skipped_cols=%" PRId64
            ", exact_mask=%d, etree_mask=%d"
            ", tail_work=%.6g, over_work=%.6g"
            ", skipped_work=%.6g\n",
            what,
            stats != NULL ? stats->fast_repaired_tail_restart_columns : 0,
            stats != NULL
              ? stats->fast_repaired_tail_restart_overcompute_columns
              : 0,
            stats != NULL
              ? stats->fast_repaired_tail_restart_skipped_columns
              : 0,
            stats != NULL ? stats->fast_repaired_tail_restart_exact_mask : 0,
            stats != NULL ? stats->fast_repaired_tail_restart_etree_mask : 0,
            stats != NULL ? stats->fast_repaired_tail_restart_work : 0.0,
            stats != NULL
              ? stats->fast_repaired_tail_restart_overcompute_work
              : 0.0,
            stats != NULL
              ? stats->fast_repaired_tail_restart_skipped_work
              : 0.0);
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
  if (!require_build_feature_stats(&stats, "csc")) return 0;
  if (!require_parallel_model_stats(&stats, "csc")) return 0;
  if (!require_row_auto_model_stats(&stats, "csc")) return 0;
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

static int test_repeated_factor_low_rcond_contract(void) {
  enum { N = 512 };
  int32_t *ap = (int32_t *)malloc((N + 1u) * sizeof(*ap));
  int32_t *ai = (int32_t *)malloc(N * sizeof(*ai));
  double *ax0 = (double *)malloc(N * sizeof(*ax0));
  double *ax1 = (double *)malloc(N * sizeof(*ax1));
  double *b = (double *)malloc(N * sizeof(*b));
  double *x = (double *)calloc(N, sizeof(*x));
  int ok = ap != NULL && ai != NULL && ax0 != NULL && ax1 != NULL &&
    b != NULL && x != NULL;
  for (int32_t col = 0; ok && col < N; ++col) {
    ap[col] = col;
    ai[col] = col;
    ax0[col] = 1.0;
    ax1[col] = 1.0;
    b[col] = 1.0;
  }
  if (ok) {
    ap[N] = N;
    ax1[0] = 1.0e-12;
    b[0] = ax1[0];
  }

  kls_options options;
  kls_default_options(&options);
  options.threads = 1;
  options.ordering = KLS_ORDERING_NATURAL;
  options.orientation = KLS_ORIENTATION_NORMAL;
  options.scale = -1;
  options.use_btf = 0;
  options.fast_factor = 0;
  options.static_pivoting = 0;

  kls_solver *solver = NULL;
  if (ok && !require_ok(kls_create(&solver),
                        "create repeated-factor low-rcond")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(
                          solver, KLS_INDEX_INT32, N, ap, ai, 0, &options),
                        "analyze repeated-factor low-rcond")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "initial repeated-factor low-rcond")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax1),
                        "changed repeated-factor low-rcond")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve repeated-factor low-rcond")) ok = 0;

  kls_stats stats;
  memset(&stats, 0, sizeof(stats));
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats repeated-factor low-rcond")) ok = 0;
  if (ok && !(stats.rcond > 0.0 && stats.rcond < sqrt(2.2204460492503131e-16))) {
    fprintf(stderr, "repeated factor did not retain low-rcond numeric: %.17g\n",
            stats.rcond);
    ok = 0;
  }
  for (int32_t i = 0; ok && i < N; ++i) {
    if (fabs(x[i] - 1.0) > 1.0e-8) {
      fprintf(stderr,
              "stale repeated-factor contract at %d: %.17g != 1\n",
              (int)i, x[i]);
      ok = 0;
    }
  }

  kls_destroy(solver);
  free(ap);
  free(ai);
  free(ax0);
  free(ax1);
  free(b);
  free(x);
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

static int test_refactor_solve_api(void) {
  const int32_t ap[] = {0, 2, 5, 7};
  const int32_t ai[] = {0, 1, 0, 1, 2, 1, 2};
  const double ax0[] = {4.0, 1.0, 1.0, 3.0, 1.0, 1.0, 2.0};
  const double ax1[] = {5.0, 1.0, 1.0, 4.0, 1.0, 1.0, 3.0};
  const double b1[] = {7.0, 12.0, 11.0};
  const double b2[] = {
    6.0, 10.0, 8.0, -1.0,
    12.0, 20.0, 16.0, -1.0
  };
  double x1[3] = {0.0, 0.0, 0.0};
  double x2[8] = {0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0};

  kls_options options;
  kls_default_options(&options);
  options.threads = 1;
  options.ordering = KLS_ORDERING_NATURAL;
  options.orientation = KLS_ORIENTATION_NORMAL;
  options.use_btf = 0;
  options.scale = -1;

  kls_solver *solver = NULL;
  int ok = require_ok(kls_create(&solver), "create refactor-solve");
  if (ok && !require_ok(kls_analyze_csc(
                          solver, KLS_INDEX_INT32, 3, ap, ai, 0, &options),
                        "analyze refactor-solve")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor refactor-solve")) ok = 0;
  if (ok && !require_ok(kls_refactor_solve(
                          solver, ax1, 1, b1, 0, x1, 0),
                        "one-rhs refactor-solve")) ok = 0;
  if (ok && !(close_enough(x1[0], 1.0) &&
              close_enough(x1[1], 2.0) &&
              close_enough(x1[2], 3.0))) {
    fprintf(stderr, "unexpected refactor-solve solution: %.17g %.17g %.17g\n",
            x1[0], x1[1], x1[2]);
    ok = 0;
  }

  /* Multi-RHS and strided frames deliberately take the compatibility path;
     the combined API must remain exactly a refactor followed by a solve. */
  if (ok && !require_ok(kls_refactor_solve(
                          solver, ax0, 2, b2, 4, x2, 4),
                        "strided multi-rhs refactor-solve")) ok = 0;
  const double expected[] = {1.0, 2.0, 3.0, 2.0, 4.0, 6.0};
  for (int rhs = 0; ok && rhs < 2; ++rhs) {
    for (int row = 0; row < 3; ++row) {
      const double got = x2[4 * rhs + row];
      const double want = expected[3 * rhs + row];
      if (!close_enough(got, want)) {
        fprintf(stderr,
                "unexpected strided refactor-solve solution[%d,%d]: "
                "%.17g != %.17g\n", rhs, row, got, want);
        ok = 0;
        break;
      }
    }
  }

  kls_destroy(solver);
  return ok;
}

static int test_serial_backend(void) {
  const int32_t ap[] = {0, 2, 5, 7};
  const int32_t ai[] = {0, 1, 0, 1, 2, 1, 2};
  const double ax[] = {4.0, 2.0, 1.0, 3.0, 1.0, 1.0, 2.0};
  const double csr_ax[] = {4.0, 1.0, 2.0, 3.0, 1.0, 1.0, 2.0};
  const double b[] = {6.0, 11.0, 8.0};
  double x[3] = {0.0, 0.0, 0.0};

  kls_options options;
  kls_default_options(&options);
  options.backend = KLS_BACKEND_SERIAL;
  options.threads = 2;
  kls_solver *solver = NULL;
  if (!require_ok(kls_create(&solver), "create serial reject")) return 0;
  const int unsupported = kls_analyze_csc(
    solver, KLS_INDEX_INT32, 3, ap, ai, 0, &options);
  if (unsupported != KLS_ERR_UNSUPPORTED) {
    fprintf(stderr, "serial backend accepted two threads: %d\n", unsupported);
    kls_destroy(solver);
    return 0;
  }
  kls_destroy(solver);

  solver = NULL;
  kls_default_options(&options);
  options.backend = KLS_BACKEND_SERIAL;
  if (!require_ok(kls_create(&solver), "create serial")) return 0;
  if (!require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 3, ap, ai, 0,
                                  &options), "analyze serial")) return 0;
  if (!require_ok(kls_factor(solver, ax), "factor serial")) return 0;
  if (!require_ok(kls_refactor(solver, ax), "refactor serial")) return 0;
  if (!require_ok(kls_solve(solver, 1, b, 0, x, 0), "solve serial")) return 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (!require_ok(kls_get_stats(solver, &stats), "stats serial")) return 0;
  const int ok = close_enough(x[0], 1.0) && close_enough(x[1], 2.0) &&
                 close_enough(x[2], 3.0) &&
                 stats.last_factor_path == KLS_FACTOR_PATH_SERIAL &&
                 stats.selected_ordering == KLS_ORDERING_AMD &&
                 stats.selected_orientation == KLS_ORIENTATION_NORMAL &&
                 stats.selected_scale == -1 &&
                 stats.factor_seconds > 0.0 &&
                 stats.refactor_seconds > 0.0;
  if (!ok) {
    fprintf(stderr,
            "unexpected serial result/path: x=%.17g %.17g %.17g path=%s "
            "ordering=%s orientation=%s scale=%d factor=%g refactor=%g\n",
            x[0], x[1], x[2],
            kls_factor_path_name(stats.last_factor_path),
            kls_ordering_name(stats.selected_ordering),
            kls_orientation_name(stats.selected_orientation),
            stats.selected_scale, stats.factor_seconds,
            stats.refactor_seconds);
  }
  kls_destroy(solver);
  if (!ok) return 0;

  memset(x, 0, sizeof(x));
  solver = NULL;
  kls_default_options(&options);
  options.backend = KLS_BACKEND_SERIAL;
  if (!require_ok(kls_create(&solver), "create serial csr")) return 0;
  if (!require_ok(kls_analyze_csr(solver, KLS_INDEX_INT32, 3, ap, ai, 0,
                                  &options), "analyze serial csr")) return 0;
  if (!require_ok(kls_factor(solver, csr_ax), "factor serial csr")) return 0;
  if (!require_ok(kls_refactor(solver, csr_ax), "refactor serial csr")) return 0;
  if (!require_ok(kls_solve(solver, 1, b, 0, x, 0), "solve serial csr")) return 0;
  stats.struct_size = sizeof(stats);
  if (!require_ok(kls_get_stats(solver, &stats), "stats serial csr")) return 0;
  const int csr_ok =
    close_enough(x[0], 1.0) && close_enough(x[1], 2.0) &&
    close_enough(x[2], 3.0) &&
    stats.selected_orientation == KLS_ORIENTATION_TRANSPOSE &&
    stats.last_factor_path == KLS_FACTOR_PATH_SERIAL;
  if (!csr_ok) {
    fprintf(stderr,
            "unexpected serial csr result/path: x=%.17g %.17g %.17g "
            "path=%s orientation=%s\n",
            x[0], x[1], x[2],
            kls_factor_path_name(stats.last_factor_path),
            kls_orientation_name(stats.selected_orientation));
  }
  kls_destroy(solver);
  if (!csr_ok) return 0;

  /* The backend field occupies the prior ABI's tail padding.  Model an old
     binary with nonzero interior/tail padding and a realistic sizeof-based
     struct_size: neither byte sequence may be mistaken for SERIAL. */
  typedef struct legacy_kls_options {
    size_t struct_size;
    int threads;
    kls_ordering ordering;
    kls_orientation orientation;
    int use_btf;
    int scale;
    double pivot_tolerance;
    double memory_growth;
    int halt_if_singular;
    int fast_factor;
    int static_pivoting;
  } legacy_kls_options;
  legacy_kls_options legacy_options;
  memset(&legacy_options, 0xa5, sizeof(legacy_options));
  legacy_options.struct_size = sizeof(legacy_options);
  legacy_options.threads = 1;
  legacy_options.ordering = KLS_ORDERING_NATURAL;
  legacy_options.orientation = KLS_ORIENTATION_NORMAL;
  legacy_options.use_btf = 0;
  legacy_options.scale = -1;
  legacy_options.pivot_tolerance = 0.001;
  legacy_options.memory_growth = 1.5;
  legacy_options.halt_if_singular = 1;
  legacy_options.fast_factor = 1;
  legacy_options.static_pivoting = 1;
  solver = NULL;
  if (sizeof(legacy_options) > sizeof(kls_options)) {
    fprintf(stderr, "legacy options ABI grew past current options: %zu vs %zu\n",
            sizeof(legacy_options), sizeof(kls_options));
    return 0;
  }
  if (!require_ok(kls_create(&solver), "create legacy options")) return 0;
  if (!require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 3, ap, ai, 0,
                                  (const kls_options *)&legacy_options),
                  "analyze legacy options")) return 0;
  kls_destroy(solver);
  return 1;
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
  if (getenv("KLS_SMOKE_POISON_NOSOLVE") == NULL &&
      !require_ok(kls_solve(solver, 1, b, 0, x, 0), "solve sparse diagonal")) return 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (!require_ok(kls_get_stats(solver, &stats), "stats sparse diagonal")) return 0;
  if (getenv("KLS_SMOKE_POISON_NOSOLVE") == NULL && stats.selected_scale != 1) {
    fprintf(stderr, "unexpected sparse-diagonal auto scale: %d\n", stats.selected_scale);
    return 0;
  }

  const int ok = getenv("KLS_SMOKE_POISON_NOSOLVE") != NULL ||
    (close_enough(x[0], 1.0) && close_enough(x[1], 2.0) && close_enough(x[2], 3.0));
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

static int test_transposed_low_work_btf_native_solve(void) {
  const int32_t n = 6;
  const int32_t ap[7] = {0, 2, 5, 8, 10, 13, 16};
  const int32_t ai[16] = {
    0, 1, 0, 1, 2, 0, 1, 2,
    3, 4, 3, 4, 5, 3, 4, 5
  };
  const double initial[16] = {
    4.0, 1.0, 2.0, 3.0, 5.0, 7.0, 1.0, 2.0,
    5.0, 2.0, 1.0, 4.0, 3.0, 2.0, 1.0, 6.0
  };
  const double expected[6] = {1.0, -2.0, 0.5, 3.0, -1.0, 2.0};
  double changed[16];
  double b[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double x[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  memcpy(changed, initial, sizeof(changed));
  changed[0] *= 1.01;
  changed[11] *= 0.99;
  for (int col = 0; col < n; ++col) {
    for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
      b[ai[p]] += changed[p] * expected[col];
    }
  }

  kls_options options;
  kls_default_options(&options);
  options.threads = 8;
  options.ordering = KLS_ORDERING_NATURAL;
  options.orientation = KLS_ORIENTATION_TRANSPOSE;
  options.scale = -1;
  options.use_btf = 1;
  options.static_pivoting = 0;

  kls_solver *solver = NULL;
  int ok = require_ok(kls_create(&solver), "create transposed low-work BTF") &&
    require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, n, ap, ai, 0,
                               &options),
               "analyze transposed low-work BTF") &&
    require_ok(kls_factor(solver, initial),
               "factor transposed low-work BTF") &&
    require_ok(kls_refactor(solver, changed),
               "refactor transposed low-work BTF") &&
    require_ok(kls_solve(solver, 1, b, 0, x, 0),
               "solve transposed low-work BTF");

  kls_stats stats;
  memset(&stats, 0, sizeof(stats));
  stats.struct_size = sizeof(stats);
  if (ok) {
    ok = require_ok(kls_get_stats(solver, &stats),
                    "stats transposed low-work BTF");
  }
  if (ok) {
    ok = stats.selected_orientation == KLS_ORIENTATION_TRANSPOSE &&
      stats.nblocks >= 2 && fabs(stats.estimated_flops) < 100000.0 &&
      stats.nnz_l + stats.nnz_u <= 12 * n;
    for (int row = 0; ok && row < n; ++row) {
      ok = close_enough(x[row], expected[row]);
    }
    if (!ok) {
      fprintf(stderr,
              "unexpected transposed low-work BTF result: orientation=%s "
              "blocks=%" PRId64 " flops=%.17g fill=%" PRId64
              " x=(%.17g %.17g %.17g %.17g %.17g %.17g)\n",
              kls_orientation_name(stats.selected_orientation), stats.nblocks,
              stats.estimated_flops, stats.nnz_l + stats.nnz_u,
              x[0], x[1], x[2], x[3], x[4], x[5]);
    }
  }
  kls_destroy(solver);
  return ok;
}

static int test_tiny_singleton_btf_solve(void) {
  enum { N = 12, NNZ = 18 };
  const int32_t ap[N + 1] = {
    0, 1, 2, 4, 5, 7, 8, 10, 11, 13, 14, 16, 18
  };
  const int32_t ai[NNZ] = {
    0, 1, 0, 2, 3, 1, 4, 5, 2, 6, 7, 3, 8, 9, 4, 10, 5, 11
  };
  double values[NNZ];
  double changed[NNZ];
  double truth[N];
  double rhs[N];
  double trhs[N];
  double x[N];
  for (int col = 0; col < N; ++col) {
    for (int p = ap[col]; p < ap[col + 1]; ++p) {
      const int diagonal = ai[p] == col;
      values[p] = diagonal ? 2.0 + 0.125 * (double)(p % 3)
                           : -0.0625 * (double)(1 + p % 2);
      changed[p] = values[p] * (1.0 + 0.001 * (double)(1 + p % 3));
    }
  }
  for (int i = 0; i < N; ++i) {
    truth[i] = 0.5 + 0.03125 * (double)i;
    rhs[i] = 0.0;
    trhs[i] = 0.0;
  }
  for (int col = 0; col < N; ++col) {
    for (int p = ap[col]; p < ap[col + 1]; ++p) {
      rhs[ai[p]] += values[p] * truth[col];
      trhs[col] += values[p] * truth[ai[p]];
    }
  }

  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_NATURAL;
  options.orientation = KLS_ORIENTATION_TRANSPOSE;
  options.scale = 1;
  options.use_btf = 1;

  kls_solver *solver = NULL;
  int ok = require_ok(kls_create(&solver), "create tiny singleton BTF") &&
    require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, N, ap, ai, 0,
                               &options),
               "analyze tiny singleton BTF") &&
    require_ok(kls_factor(solver, values), "factor tiny singleton BTF") &&
    require_ok(kls_solve(solver, 1, rhs, 0, x, 0),
               "solve tiny singleton BTF");
  for (int i = 0; ok && i < N; ++i) {
    ok = close_enough(x[i], truth[i]);
  }
  if (ok) {
    memcpy(x, trhs, sizeof(x));
    ok = require_ok(kls_solve_transpose(solver, 1, x, 0, x, 0),
                    "transpose in-place tiny singleton BTF");
  }
  for (int i = 0; ok && i < N; ++i) {
    ok = close_enough(x[i], truth[i]);
  }
  if (ok) {
    for (int i = 0; i < N; ++i) {
      rhs[i] = 0.0;
    }
    for (int col = 0; col < N; ++col) {
      for (int p = ap[col]; p < ap[col + 1]; ++p) {
        rhs[ai[p]] += changed[p] * truth[col];
      }
    }
    ok = require_ok(kls_refactor_solve(solver, changed, 1, rhs, 0, x, 0),
                    "changed refactor-solve tiny singleton BTF");
  }
  for (int i = 0; ok && i < N; ++i) {
    ok = close_enough(x[i], truth[i]);
  }

  kls_stats stats;
  memset(&stats, 0, sizeof(stats));
  stats.struct_size = sizeof(stats);
  if (ok) {
    ok = require_ok(kls_get_stats(solver, &stats),
                    "stats tiny singleton BTF") &&
      stats.tiny_singleton_solve_eligible == 1 &&
      stats.tiny_singleton_solve_count >= 3;
  }
  if (!ok) {
    fprintf(stderr,
            "tiny singleton BTF failed: eligible=%d count=%" PRId64 "\n",
            stats.tiny_singleton_solve_eligible,
            stats.tiny_singleton_solve_count);
  }
  kls_destroy(solver);
  return ok;
}

static int test_compact_singleton_run_solve(void) {
  enum {
    N = 8192,
    CORE_N = 2048,
    FRINGE_N = N - CORE_N,
    NNZ = 2 * CORE_N + FRINGE_N
  };
  int32_t *ap = (int32_t *)malloc(((size_t)N + 1u) * sizeof(*ap));
  int32_t *ai = (int32_t *)malloc((size_t)NNZ * sizeof(*ai));
  double *initial = (double *)malloc((size_t)NNZ * sizeof(*initial));
  double *changed = (double *)malloc((size_t)NNZ * sizeof(*changed));
  double *expected = (double *)malloc((size_t)N * sizeof(*expected));
  double *b = (double *)calloc((size_t)N, sizeof(*b));
  double *x = (double *)calloc((size_t)N, sizeof(*x));
  if (ap == NULL || ai == NULL || initial == NULL || changed == NULL ||
      expected == NULL || b == NULL || x == NULL) {
    free(ap);
    free(ai);
    free(initial);
    free(changed);
    free(expected);
    free(b);
    free(x);
    return 0;
  }

  const char *env_name = "KLS_DISABLE_COMPACT_SINGLETON_RUN_SOLVE";
  const char *env_value = getenv(env_name);
  const int had_env = env_value != NULL;
  char *saved_env = env_value != NULL ? strdup(env_value) : NULL;
  int ok = (env_value == NULL || saved_env != NULL) &&
    unsetenv(env_name) == 0;

  int32_t p = 0;
  ap[0] = 0;
  for (int32_t col = 0; col < CORE_N; ++col) {
    const int32_t successor = col + 1 == CORE_N ? 0 : col + 1;
    if (successor < col) {
      ai[p] = successor;
      initial[p] = -0.125;
      changed[p] = initial[p] * (p % 101 == 0 ? 1.0005 : 1.0);
      ++p;
    }
    ai[p] = col;
    initial[p] = 4.0;
    changed[p] = initial[p] * (p % 101 == 0 ? 1.0005 : 1.0);
    ++p;
    if (successor > col) {
      ai[p] = successor;
      initial[p] = -0.125;
      changed[p] = initial[p] * (p % 101 == 0 ? 1.0005 : 1.0);
      ++p;
    }
    ap[col + 1] = p;
  }
  for (int32_t col = CORE_N; col < N; ++col) {
    ai[p] = col;
    initial[p] = 2.0;
    changed[p] = initial[p] * (p % 101 == 0 ? 1.0005 : 1.0);
    ++p;
    ap[col + 1] = p;
  }
  if (p != NNZ) {
    ok = 0;
  }
  for (int32_t col = 0; col < N; ++col) {
    expected[col] = 0.75 + 0.03125 * (double)(col % 17);
    for (int32_t entry = ap[col]; entry < ap[col + 1]; ++entry) {
      b[ai[entry]] += changed[entry] * expected[col];
    }
  }

  kls_options options;
  kls_default_options(&options);
  options.threads = 8;
  options.ordering = KLS_ORDERING_NATURAL;
  options.orientation = KLS_ORIENTATION_NORMAL;
  options.scale = -1;
  options.use_btf = 1;

  kls_solver *solver = NULL;
  if (ok && !require_ok(kls_create(&solver),
                        "create compact singleton runs")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(
                          solver, KLS_INDEX_INT32, N, ap, ai, 0, &options),
                        "analyze compact singleton runs")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, initial),
                        "factor compact singleton runs")) ok = 0;
  if (ok && !require_ok(kls_refactor(solver, changed),
                        "refactor compact singleton runs")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve compact singleton runs")) ok = 0;

  kls_stats stats;
  memset(&stats, 0, sizeof(stats));
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats compact singleton runs")) ok = 0;
  if (ok) {
    ok = stats.selected_orientation == KLS_ORIENTATION_NORMAL &&
      stats.selected_ordering == KLS_ORDERING_NATURAL &&
      stats.selected_btf == 1 && stats.nblocks == FRINGE_N + 1 &&
      stats.max_block == CORE_N && stats.compact_solve_index_bytes == 2 &&
      stats.compact_solve_singleton_run_blocks == FRINGE_N &&
      stats.compact_solve_singleton_run_max >= FRINGE_N / 2 &&
      stats.compact_solve_singleton_run_eligible == 1;
    for (int32_t row = 0; ok && row < N; ++row) {
      ok = close_enough(x[row], expected[row]);
    }
    if (!ok) {
      fprintf(stderr,
              "unexpected compact singleton-run result: btf=%d "
              "blocks=%" PRId64 " max=%" PRId64 " compact=%d "
              "singletons=%" PRId64 " run=%" PRId64 " eligible=%d"
              " x0=%.17g xlast=%.17g\n",
              stats.selected_btf, stats.nblocks, stats.max_block,
              stats.compact_solve_index_bytes,
              stats.compact_solve_singleton_run_blocks,
              stats.compact_solve_singleton_run_max,
              stats.compact_solve_singleton_run_eligible, x[0], x[N - 1]);
    }
  }
  kls_destroy(solver);

  solver = NULL;
  memset(x, 0, (size_t)N * sizeof(*x));
  if (ok && setenv(env_name, "1", 1) != 0) ok = 0;
  if (ok && !require_ok(kls_create(&solver),
                        "create disabled compact singleton runs")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(
                          solver, KLS_INDEX_INT32, N, ap, ai, 0, &options),
                        "analyze disabled compact singleton runs")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, initial),
                        "factor disabled compact singleton runs")) ok = 0;
  if (ok && !require_ok(kls_refactor(solver, changed),
                        "refactor disabled compact singleton runs")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve disabled compact singleton runs")) ok = 0;
  memset(&stats, 0, sizeof(stats));
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats disabled compact singleton runs")) ok = 0;
  if (ok) {
    ok = stats.compact_solve_index_bytes == 2 &&
      stats.compact_solve_singleton_run_blocks == 0 &&
      stats.compact_solve_singleton_run_max == 0 &&
      stats.compact_solve_singleton_run_eligible == 0;
    for (int32_t row = 0; ok && row < N; ++row) {
      ok = close_enough(x[row], expected[row]);
    }
    if (!ok) {
      fprintf(stderr,
              "unexpected disabled compact singleton-run result: "
              "compact=%d singletons=%" PRId64 " run=%" PRId64
              " eligible=%d"
              " x0=%.17g xlast=%.17g\n",
              stats.compact_solve_index_bytes,
              stats.compact_solve_singleton_run_blocks,
              stats.compact_solve_singleton_run_max,
              stats.compact_solve_singleton_run_eligible, x[0], x[N - 1]);
    }
  }
  kls_destroy(solver);

  /* Keep the same order and entry count, but turn two thirds of the scalar
     fringe into independent two-cycles.  The longest scalar run remains
     large, while singleton coverage falls below the representation's 75%
     payoff boundary. */
  enum { PAIRED_FRINGE_BEGIN = CORE_N + FRINGE_N / 3 };
  if (ok && unsetenv(env_name) != 0) ok = 0;
  for (int32_t col = PAIRED_FRINGE_BEGIN; col < N; col += 2) {
    ai[ap[col]] = col + 1;
    ai[ap[col + 1]] = col;
  }
  memset(b, 0, (size_t)N * sizeof(*b));
  memset(x, 0, (size_t)N * sizeof(*x));
  for (int32_t col = 0; col < N; ++col) {
    for (int32_t entry = ap[col]; entry < ap[col + 1]; ++entry) {
      b[ai[entry]] += changed[entry] * expected[col];
    }
  }

  solver = NULL;
  if (ok && !require_ok(kls_create(&solver),
                        "create paired-fringe singleton guard")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(
                          solver, KLS_INDEX_INT32, N, ap, ai, 0, &options),
                        "analyze paired-fringe singleton guard")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, initial),
                        "factor paired-fringe singleton guard")) ok = 0;
  if (ok && !require_ok(kls_refactor(solver, changed),
                        "refactor paired-fringe singleton guard")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve paired-fringe singleton guard")) ok = 0;
  memset(&stats, 0, sizeof(stats));
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats paired-fringe singleton guard")) ok = 0;
  if (ok) {
    const int64_t expected_blocks =
      1 + FRINGE_N / 3 + (FRINGE_N - FRINGE_N / 3) / 2;
    ok = stats.selected_btf == 1 && stats.nblocks == expected_blocks &&
      stats.max_block == CORE_N && stats.compact_solve_index_bytes == 2 &&
      stats.compact_solve_singleton_run_blocks == 0 &&
      stats.compact_solve_singleton_run_max == 0 &&
      stats.compact_solve_singleton_run_eligible == 0;
    for (int32_t row = 0; ok && row < N; ++row) {
      ok = close_enough(x[row], expected[row]);
    }
    if (!ok) {
      fprintf(stderr,
              "unexpected paired-fringe singleton guard: btf=%d "
              "blocks=%" PRId64 "/%" PRId64 " max=%" PRId64
              " compact=%d singletons=%" PRId64 " run=%" PRId64
              " eligible=%d x0=%.17g xlast=%.17g\n",
              stats.selected_btf, stats.nblocks, expected_blocks,
              stats.max_block, stats.compact_solve_index_bytes,
              stats.compact_solve_singleton_run_blocks,
              stats.compact_solve_singleton_run_max,
              stats.compact_solve_singleton_run_eligible, x[0], x[N - 1]);
    }
  }
  kls_destroy(solver);

  if (!(had_env && saved_env == NULL) &&
      !restore_env_value(env_name, had_env, saved_env)) {
    ok = 0;
  }
  free(saved_env);
  free(ap);
  free(ai);
  free(initial);
  free(changed);
  free(expected);
  free(b);
  free(x);
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
  options.orientation = KLS_ORIENTATION_NORMAL;
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
  if (ok && !require_tail_overcompute_bounds(&stats,
                                             "pivot-check fallback")) {
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

static int test_fast_factor_rowwise_u_pivot_reject(void) {
  const int32_t ap[] = {0, 2, 4};
  const int32_t ai[] = {0, 1, 0, 1};
  const double ax0[] = {2.0, 1.0, 1.0, 2.0};
  const double ax1[] = {1.0e-3, 2.0e-2, 10.0, 10.0};
  const double b[] = {20.001, 20.02};
  double x[2] = {0.0, 0.0};

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_NATURAL;
  options.orientation = KLS_ORIENTATION_NORMAL;
  options.use_btf = 0;
  options.scale = -1;
  options.pivot_tolerance = 0.001;
  options.threads = 1;
  options.static_pivoting = 0;
  const char *saved_row_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_row = saved_row_value != NULL ? strdup(saved_row_value) : NULL;
  const int had_row = saved_row_value != NULL;
  const char *saved_checked_value =
    getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  char *saved_checked =
    saved_checked_value != NULL ? strdup(saved_checked_value) : NULL;
  const int had_checked = saved_checked_value != NULL;

  int ok = 1;
  if ((had_row && saved_row == NULL) ||
      (had_checked && saved_checked == NULL)) {
    fprintf(stderr, "failed to save row-refactor environment\n");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR");
    ok = 0;
  }
  if (!require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 2, ap, ai, 0,
                                        &options),
                        "analyze rowwise-U pivot reject")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor rowwise-U base")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax1),
                        "factor rowwise-U repair")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve rowwise-U repair")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats rowwise-U pivot reject")) {
    ok = 0;
  }
  if (ok && (stats.fast_rejected_pivot != 0 ||
             stats.fast_rejected_pivot_col != 0 ||
             stats.fast_rejected_row != 0 ||
             stats.fast_rejected_refresh_state !=
               KLS_FAST_REJECT_REFRESH_ALL)) {
    fprintf(stderr,
            "unexpected rowwise-U reject location: pivot=%" PRId64
            ", col=%" PRId64 ", row=%" PRId64 ", refresh=%d\n",
            stats.fast_rejected_pivot,
            stats.fast_rejected_pivot_col,
            stats.fast_rejected_row,
            stats.fast_rejected_refresh_state);
    fprintf(stderr,
            "rowwise-U reject details: scale=%d, tol=%.17g, ratio=%.17g"
            ", pivot=%.17g, candidate=%.17g\n",
            stats.selected_scale,
            stats.selected_pivot_tolerance,
            stats.fast_rejected_multiplier_abs,
            stats.fast_rejected_pivot_abs,
            stats.fast_rejected_candidate_abs);
    ok = 0;
  }
  if (ok && (stats.fast_rejected_pivot_abs < 9.0e-4 ||
             stats.fast_rejected_pivot_abs > 1.1e-3 ||
             stats.fast_rejected_candidate_abs < 9.9 ||
             stats.fast_rejected_multiplier_abs < 9000.0)) {
    fprintf(stderr,
            "unexpected rowwise-U reject magnitudes: ratio=%.17g"
            ", pivot=%.17g, rowmax=%.17g\n",
            stats.fast_rejected_multiplier_abs,
            stats.fast_rejected_pivot_abs,
            stats.fast_rejected_candidate_abs);
    ok = 0;
  }
  if (ok && stats.fast_block_restarts < 1) {
    fprintf(stderr, "rowwise-U reject did not attempt block repair: %d\n",
            stats.fast_block_restarts);
    ok = 0;
  }
  if (ok && (!close_enough(x[0], 1.0) || !close_enough(x[1], 2.0))) {
    fprintf(stderr, "unexpected rowwise-U repair solution: %.17g %.17g\n",
            x[0], x[1]);
    ok = 0;
  }

  if (!restore_env_value("KLS_ENABLE_ROW_REFACTOR", had_row, saved_row)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_CHECKED_ROW_REFACTOR", had_checked,
                         saved_checked)) {
    ok = 0;
  }
  free(saved_row);
  free(saved_checked);
  kls_destroy(solver);
  return ok;
}

static int run_fast_factor_root_independent_tail_restart(int scale,
                                                         const char *label) {
  const int32_t ap[] = {0, 2, 4, 5};
  const int32_t ai[] = {0, 1, 0, 1, 2};
  const double ax0[] = {2.0, 1.0, 1.0, 2.0, 3.0};
  const double ax1[] = {1.0e-12, 1.0, 1.0, 2.0, 3.0};
  const double b[] = {2.000000000001, 5.0, 9.0};
  double x[3] = {0.0, 0.0, 0.0};

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_NATURAL;
  options.use_btf = 0;
  options.scale = scale;
  options.pivot_tolerance = 0.001;
  options.threads = 2;
  options.static_pivoting = 0;

  int ok = 1;
  if (!require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 3, ap, ai, 0,
                                        &options),
                        "analyze root independent tail")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor root independent base")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax1),
                        "factor root independent repair")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve root independent repair")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats root independent tail")) {
    ok = 0;
  }
  if (ok && (stats.selected_scale != scale ||
             stats.fast_rejected_pivot != 0 ||
             stats.fast_rejected_pivot_col != 0 ||
             stats.fast_rejected_block_start != 0 ||
             stats.fast_rejected_block_size != 3 ||
             stats.fast_rejected_suffix_columns != 3 ||
             stats.fast_rejected_refresh_state !=
               KLS_FAST_REJECT_REFRESH_PREFIX ||
             stats.fast_rejected_pivoting_tail_columns != 2 ||
             stats.fast_rejected_pivoting_tail_last != 1 ||
             stats.fast_rejected_pivoting_tail_suffix_exact != 0 ||
             stats.fast_rejected_pivoting_tail_suffix_overcompute_columns != 1 ||
             stats.fast_block_restarts != 1 ||
             stats.fast_tail_restarts != 0 ||
             stats.fast_kls_block_restart_last_row_pipeline != 1 ||
             stats.fast_kls_block_restart_last_row_pipeline_rows != 2 ||
             stats.fast_kls_block_restart_last_row_pipeline_suffix_rows != 1 ||
             stats.fast_kls_block_restart_last_row_pipeline_etree_tail != 1 ||
             stats.fast_kls_block_restart_row_pipeline_etree_tail_count < 1 ||
             stats.fast_kls_block_restart_last_row_pipeline_etree_tail_rows !=
               2 ||
             stats
               .fast_kls_block_restart_last_row_pipeline_etree_tail_gap_rows !=
               0 ||
             stats
               .fast_kls_block_restart_last_row_pipeline_etree_tail_exact_mask !=
               1)) {
    fprintf(stderr,
            "unexpected %s stats: scale=%d/%d, pivot=%" PRId64
            ", col=%" PRId64 ", block=[%" PRId64 ",%" PRId64 ")"
            ", suffix=%" PRId64 ", refresh=%d, tail_cols=%" PRId64
            ", tail_last=%" PRId64 ", suffix_exact=%d"
            ", suffix_over=%" PRId64
            ", block_restarts=%d, tail_restarts=%d"
            ", pipeline=%d, pipeline_rows=%" PRId64
            ", pipeline_suffix=%" PRId64
            ", etree_tail=%d/%" PRId64 ", etree_rows=%" PRId64
            ", etree_gaps=%" PRId64 ", etree_exact=%d\n",
            label,
            stats.selected_scale,
            scale,
            stats.fast_rejected_pivot,
            stats.fast_rejected_pivot_col,
            stats.fast_rejected_block_start,
            stats.fast_rejected_block_start + stats.fast_rejected_block_size,
            stats.fast_rejected_suffix_columns,
            stats.fast_rejected_refresh_state,
            stats.fast_rejected_pivoting_tail_columns,
            stats.fast_rejected_pivoting_tail_last,
            stats.fast_rejected_pivoting_tail_suffix_exact,
            stats.fast_rejected_pivoting_tail_suffix_overcompute_columns,
            stats.fast_block_restarts,
            stats.fast_tail_restarts,
            stats.fast_kls_block_restart_last_row_pipeline,
            stats.fast_kls_block_restart_last_row_pipeline_rows,
            stats.fast_kls_block_restart_last_row_pipeline_suffix_rows,
            stats.fast_kls_block_restart_last_row_pipeline_etree_tail,
            stats.fast_kls_block_restart_row_pipeline_etree_tail_count,
            stats.fast_kls_block_restart_last_row_pipeline_etree_tail_rows,
            stats.fast_kls_block_restart_last_row_pipeline_etree_tail_gap_rows,
            stats
              .fast_kls_block_restart_last_row_pipeline_etree_tail_exact_mask);
    ok = 0;
  }
  if (ok && !require_pivoting_tail_plan(&stats,
                                        label)) {
    ok = 0;
  }
  if (ok && !require_tail_overcompute_bounds(&stats,
                                             label)) {
    ok = 0;
  }
  if (ok && (!close_enough(x[0], 1.0) || !close_enough(x[1], 2.0) ||
             !close_enough(x[2], 3.0))) {
    fprintf(stderr,
            "unexpected %s solution: %.17g %.17g %.17g\n",
            label, x[0], x[1], x[2]);
    ok = 0;
  }

  kls_destroy(solver);
  return ok;
}

static int test_fast_factor_root_independent_tail_restart(void) {
  int ok = 1;
  if (!run_fast_factor_root_independent_tail_restart(
        -1, "root independent tail")) {
    ok = 0;
  }
  if (!run_fast_factor_root_independent_tail_restart(
        2, "scaled root independent tail")) {
    ok = 0;
  }
  return ok;
}

static int run_fast_factor_noncontiguous_tail_gap_work_bounds(int threads,
                                                              const char *label) {
  const int32_t ap[] = {0, 2, 5, 7, 9, 11, 13};
  const int32_t ai[] = {
    0, 4,
    1, 3, 5,
    0, 2,
    2, 3,
    0, 4,
    1, 5
  };
  const double ax0[] = {
    5.0899999999999999, -0.021414242728184554,
    5.75, -0.039117352056168508, -0.03239719157472417,
    0.035656970912738221, 5.3899999999999997,
    0.024, 5.9399999999999995,
    -0.024322968906720161, 5.5199999999999996,
    0.0071213640922768301, 5.9800000000000004
  };
  const double ax1[] = {
    9.9999999999999998e-13, -2.0214142427281847,
    5.75, -0.039117352056168508, -0.03239719157472417,
    0.035656970912738221, 5.3899999999999997,
    0.024, 5.9399999999999995,
    -0.024322968906720161, 5.5199999999999996,
    0.0071213640922768301, 5.9800000000000004
  };
  const double expected[] = {1.0, -2.0, 3.0, -4.0, 5.0, -6.0};
  double b[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double x[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  for (int32_t col = 0; col < 6; ++col) {
    for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
      b[ai[p]] += ax1[p] * expected[col];
    }
  }

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_NATURAL;
  options.use_btf = 0;
  options.scale = -1;
  options.pivot_tolerance = 0.001;
  options.static_pivoting = 0;
  options.threads = threads;

  int ok = 1;
  if (!require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 6, ap, ai, 0,
                                        &options),
                        "analyze noncontiguous gap tail")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor noncontiguous gap base")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax1),
                        "factor noncontiguous gap repair")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve noncontiguous gap repair")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats noncontiguous gap tail")) {
    ok = 0;
  }
  const int expect_row_pipeline = threads > 1;
  const int expect_tail_restart = expect_row_pipeline ? 0 : 1;
  if (ok && (stats.fast_rejected_pivot != 0 ||
             stats.fast_rejected_pivoting_tail_columns != 5 ||
             stats.fast_rejected_pivoting_tail_first != 0 ||
             stats.fast_rejected_pivoting_tail_last != 5 ||
             stats.fast_rejected_pivoting_tail_contiguous != 0 ||
             stats.fast_rejected_pivoting_tail_gap_columns != 1 ||
             stats.fast_block_restarts != 1 ||
             stats.fast_tail_restarts != expect_tail_restart ||
             stats.fast_kls_block_restart_last_row_pipeline !=
               expect_row_pipeline ||
             (expect_row_pipeline &&
              (stats.fast_kls_block_restart_last_row_pipeline_rows != 5 ||
               stats.fast_kls_block_restart_last_row_pipeline_gap_rows != 1 ||
               stats.fast_kls_block_restart_last_row_pipeline_suffix_rows != 0 ||
               stats.fast_kls_block_restart_last_row_pipeline_threads < 1 ||
               stats.fast_kls_block_restart_last_row_pipeline_etree_tail != 1 ||
               stats.fast_kls_block_restart_row_pipeline_etree_tail_count < 1 ||
               stats.fast_kls_block_restart_last_row_pipeline_etree_ready !=
                 1 ||
               stats.fast_kls_block_restart_row_pipeline_etree_ready_count <
                 1 ||
               stats
                 .fast_kls_block_restart_last_row_pipeline_etree_ready_rows !=
                 4 ||
               stats
                 .fast_kls_block_restart_last_row_pipeline_etree_ready_threads <
                 1 ||
               stats
                 .fast_kls_block_restart_last_row_pipeline_etree_prefactor !=
                 1 ||
               stats
                 .fast_kls_block_restart_row_pipeline_etree_prefactor_count <
                 1 ||
               stats
                 .fast_kls_block_restart_last_row_pipeline_etree_prefactor_rows !=
                 4 ||
               stats
                 .fast_kls_block_restart_last_row_pipeline_etree_prefactor_threads <
                 1 ||
               stats
                 .fast_kls_block_restart_last_row_pipeline_etree_prefactor_wait_rows <
                 0 ||
               stats
                 .fast_kls_block_restart_last_row_pipeline_etree_prefactor_wait_deps <
                 0 ||
               stats
                 .fast_kls_block_restart_last_row_pipeline_etree_tail_rows !=
                 5 ||
               stats
                 .fast_kls_block_restart_last_row_pipeline_etree_tail_gap_rows !=
                 1 ||
               stats
                 .fast_kls_block_restart_last_row_pipeline_etree_tail_exact_mask !=
                 1)) ||
             (!expect_row_pipeline &&
              (stats.fast_kls_block_restart_last_row_pipeline_rows != 0 ||
               stats.fast_kls_block_restart_last_row_pipeline_gap_rows != 0 ||
               stats.fast_kls_block_restart_last_row_pipeline_suffix_rows != 0 ||
               stats.fast_kls_block_restart_last_row_pipeline_threads != 0 ||
               stats.fast_kls_block_restart_last_row_pipeline_etree_tail != 0 ||
               stats.fast_kls_block_restart_row_pipeline_etree_tail_count != 0 ||
               stats.fast_kls_block_restart_last_row_pipeline_etree_ready !=
                 0 ||
               stats.fast_kls_block_restart_row_pipeline_etree_ready_count !=
                 0 ||
               stats
                 .fast_kls_block_restart_last_row_pipeline_etree_ready_rows !=
                 0 ||
               stats
                 .fast_kls_block_restart_last_row_pipeline_etree_ready_threads !=
                 0 ||
               stats
                 .fast_kls_block_restart_last_row_pipeline_etree_prefactor !=
                 0 ||
               stats
                 .fast_kls_block_restart_row_pipeline_etree_prefactor_count !=
                 0 ||
               stats
                 .fast_kls_block_restart_last_row_pipeline_etree_prefactor_rows !=
                 0 ||
               stats
                 .fast_kls_block_restart_last_row_pipeline_etree_prefactor_threads !=
                 0 ||
               stats
                 .fast_kls_block_restart_last_row_pipeline_etree_prefactor_wait_rows !=
                 0 ||
               stats
                 .fast_kls_block_restart_last_row_pipeline_etree_prefactor_wait_deps !=
                 0 ||
               stats
                 .fast_kls_block_restart_last_row_pipeline_etree_tail_rows !=
                 0 ||
               stats
                 .fast_kls_block_restart_last_row_pipeline_etree_tail_gap_rows !=
                 0 ||
               stats
                 .fast_kls_block_restart_last_row_pipeline_etree_tail_exact_mask !=
                 0)) ||
             stats.fast_repaired_tail_restart_exact_mask != 1 ||
             stats.fast_repaired_tail_restart_etree_mask != 1 ||
             stats.fast_repaired_tail_restart_columns !=
               stats.fast_rejected_pivoting_tail_columns ||
             stats.fast_repaired_tail_restart_overcompute_columns != 0 ||
             stats.fast_repaired_last_offdiag_suffix_refresh != 1 ||
             stats.fast_repaired_offdiag_suffix_refresh_count != 1 ||
             stats.fast_repaired_offdiag_full_refresh_count != 0)) {
    fprintf(stderr,
            "unexpected %s stats: pivot=%" PRId64
            ", tail_cols=%" PRId64 ", first=%" PRId64 ", last=%" PRId64
            ", contiguous=%d, gaps=%" PRId64
            ", block_restarts=%d, tail_restarts=%d/%d"
            ", pipeline=%d, pipeline_rows=%" PRId64
            ", pipeline_gaps=%" PRId64 ", pipeline_suffix=%" PRId64
            ", pipeline_threads=%" PRId64
            ", etree_tail=%d/%" PRId64 ", etree_rows=%" PRId64
            ", etree_gaps=%" PRId64 ", etree_exact=%d"
            ", etree_ready=%d/%" PRId64 "/%" PRId64 "/%" PRId64
            ", etree_prefactor=%d/%" PRId64 "/%" PRId64 "/%" PRId64
            "/%" PRId64 "/%" PRId64
            ", exact_mask=%d, etree_mask=%d, repaired_cols=%" PRId64
            ", over_cols=%" PRId64
            ", offdiag_suffix=%d, offdiag_suffix_count=%" PRId64
            ", offdiag_full_count=%" PRId64
            ", supernode_groups=%" PRId64 ", supernode_rows=%" PRId64
            ", panel_groups=%" PRId64 ", panel_rows=%" PRId64
            "\n",
            label,
            stats.fast_rejected_pivot,
            stats.fast_rejected_pivoting_tail_columns,
            stats.fast_rejected_pivoting_tail_first,
            stats.fast_rejected_pivoting_tail_last,
            stats.fast_rejected_pivoting_tail_contiguous,
            stats.fast_rejected_pivoting_tail_gap_columns,
            stats.fast_block_restarts,
            stats.fast_tail_restarts,
            expect_tail_restart,
            stats.fast_kls_block_restart_last_row_pipeline,
            stats.fast_kls_block_restart_last_row_pipeline_rows,
            stats.fast_kls_block_restart_last_row_pipeline_gap_rows,
            stats.fast_kls_block_restart_last_row_pipeline_suffix_rows,
            stats.fast_kls_block_restart_last_row_pipeline_threads,
            stats.fast_kls_block_restart_last_row_pipeline_etree_tail,
            stats.fast_kls_block_restart_row_pipeline_etree_tail_count,
            stats.fast_kls_block_restart_last_row_pipeline_etree_tail_rows,
            stats.fast_kls_block_restart_last_row_pipeline_etree_tail_gap_rows,
            stats
              .fast_kls_block_restart_last_row_pipeline_etree_tail_exact_mask,
            stats.fast_kls_block_restart_last_row_pipeline_etree_ready,
            stats.fast_kls_block_restart_row_pipeline_etree_ready_count,
            stats.fast_kls_block_restart_last_row_pipeline_etree_ready_rows,
            stats.fast_kls_block_restart_last_row_pipeline_etree_ready_threads,
            stats.fast_kls_block_restart_last_row_pipeline_etree_prefactor,
            stats.fast_kls_block_restart_row_pipeline_etree_prefactor_count,
            stats.fast_kls_block_restart_last_row_pipeline_etree_prefactor_rows,
            stats
              .fast_kls_block_restart_last_row_pipeline_etree_prefactor_threads,
            stats
              .fast_kls_block_restart_last_row_pipeline_etree_prefactor_wait_rows,
            stats
              .fast_kls_block_restart_last_row_pipeline_etree_prefactor_wait_deps,
            stats.fast_repaired_tail_restart_exact_mask,
            stats.fast_repaired_tail_restart_etree_mask,
            stats.fast_repaired_tail_restart_columns,
            stats.fast_repaired_tail_restart_overcompute_columns,
            stats.fast_repaired_last_offdiag_suffix_refresh,
            stats.fast_repaired_offdiag_suffix_refresh_count,
            stats.fast_repaired_offdiag_full_refresh_count,
            stats
              .fast_kls_block_restart_last_row_pipeline_supernode_update_groups,
            stats
              .fast_kls_block_restart_last_row_pipeline_supernode_update_rows,
            stats
              .fast_kls_block_restart_last_row_pipeline_supernode_panel_update_groups,
            stats
              .fast_kls_block_restart_last_row_pipeline_supernode_panel_update_rows);
    ok = 0;
  }
  if (ok && !require_pivoting_tail_plan(&stats, "noncontiguous gap tail")) {
    ok = 0;
  }
  if (ok && !require_tail_overcompute_bounds(&stats,
                                             "noncontiguous gap tail")) {
    ok = 0;
  }
  double max_solution_error = 0.0;
  for (int32_t i = 0; i < 6; ++i) {
    const double err = fabs(x[i] - expected[i]);
    if (err > max_solution_error) {
      max_solution_error = err;
    }
  }
  if (ok && max_solution_error > 1.0e-8) {
    fprintf(stderr,
            "unexpected noncontiguous gap solution: %.17g %.17g %.17g"
            " %.17g %.17g %.17g, max_err=%.17g\n",
            x[0], x[1], x[2], x[3], x[4], x[5], max_solution_error);
    ok = 0;
  }

  kls_destroy(solver);
  return ok;
}

static int test_fast_factor_noncontiguous_tail_gap_work_bounds(void) {
  int ok = 1;
  if (!run_fast_factor_noncontiguous_tail_gap_work_bounds(
        2, "noncontiguous row-pipeline gap tail")) {
    ok = 0;
  }
  if (!run_fast_factor_noncontiguous_tail_gap_work_bounds(
        1, "serial noncontiguous exact-mask gap tail")) {
    ok = 0;
  }
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
             stats.fast_rejected_row != 0 ||
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
             stats.fast_rejected_group_tail_groups < 1 ||
             stats.fast_rejected_group_tail_rows <
               stats.fast_rejected_row_tail_columns ||
             stats.fast_rejected_group_tail_work <= 0.0 ||
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
            ", group_tail=%" PRId64 ", group_rows=%" PRId64
            ", group_work=%.6g, etree=%" PRId64
            ", pivoting_tail=%" PRId64 "\n",
            stats.fast_rejected_block_start,
            stats.fast_rejected_block_size,
            stats.fast_rejected_suffix_columns,
            stats.fast_rejected_descendant_columns,
            stats.fast_rejected_row_tail_columns,
            stats.fast_rejected_group_tail_groups,
            stats.fast_rejected_group_tail_rows,
            stats.fast_rejected_group_tail_work,
            stats.fast_rejected_etree_columns,
            stats.fast_rejected_pivoting_tail_columns);
    ok = 0;
  }
  if (ok && !require_pivoting_tail_plan(&stats, "checked-row block restart")) {
    ok = 0;
  }
  if (ok && !require_tail_overcompute_bounds(
              &stats, "checked-row block restart")) {
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
             stats.fast_rejected_row != 1 ||
             stats.fast_rejected_refresh_state !=
               KLS_FAST_REJECT_REFRESH_PREFIX ||
             stats.fast_rejected_tail_repair_ready != 1 ||
             stats.fast_block_restarts != 1 ||
             stats.fast_tail_restarts != 0 ||
             stats.fast_kls_block_restart_last_row_pipeline != 1 ||
             stats.fast_kls_block_restart_row_pipeline_count < 1 ||
             stats.fast_kls_block_restart_last_row_pipeline_rows !=
               stats.fast_rejected_pivoting_tail_columns ||
             stats.fast_kls_block_restart_last_row_pipeline_etree_tail != 1 ||
             stats.fast_kls_block_restart_row_pipeline_etree_tail_count < 1 ||
             stats.fast_kls_block_restart_last_row_pipeline_etree_tail_rows !=
               stats.fast_rejected_pivoting_tail_columns ||
             stats.fast_kls_block_restart_last_row_pipeline_etree_tail_exact_mask !=
               1 ||
             stats.fast_rejected_pivoting_tail_columns < 2 ||
             stats.fast_kls_block_restart_last_row_pipeline_etree_ready != 1 ||
             stats.fast_kls_block_restart_row_pipeline_etree_ready_count < 1 ||
             stats.fast_kls_block_restart_last_row_pipeline_etree_ready_rows !=
               stats.fast_rejected_pivoting_tail_columns - 1 ||
             stats.fast_kls_block_restart_last_row_pipeline_etree_ready_threads <
               1 ||
             stats.fast_kls_block_restart_last_row_pipeline_prefix_rows +
                 stats.fast_kls_block_restart_last_row_pipeline_rows +
                 stats.fast_kls_block_restart_last_row_pipeline_suffix_rows !=
               stats.fast_rejected_block_size ||
             stats.fast_rejected_pivoting_tail_columns < 1 ||
             stats.fast_rejected_pivoting_tail_seed_columns !=
               stats.fast_rejected_pivoting_tail_columns ||
             stats.fast_rejected_pivoting_tail_row_seed_columns != 0 ||
             stats.fast_rejected_pivoting_tail_seed_columns <
               stats.fast_rejected_row_tail_columns ||
             stats.fast_rejected_pivoting_tail_topological != 1)) {
    fprintf(stderr,
            "unexpected parallel checked-row stats: pivot=%" PRId64
            ", col=%" PRId64 ", row=%" PRId64 ", refresh=%d"
            ", tail_ready=%d, restarts=%d, tail_restarts=%d"
            ", pipeline=%d/%" PRId64 ", pipeline_rows=%" PRId64
            ", pipeline_threads=%" PRId64
            ", pipeline_prefix=%" PRId64
            ", pipeline_suffix=%" PRId64
            ", etree_tail=%d/%" PRId64 ", etree_rows=%" PRId64
            ", etree_exact=%d"
            ", etree_ready=%d/%" PRId64 "/%" PRId64 "/%" PRId64
            ", tail_cols=%" PRId64 ", seed=%" PRId64
            ", row_seed=%" PRId64 ", row_tail=%" PRId64
            ", tail_topo=%d\n",
            stats.fast_rejected_pivot,
            stats.fast_rejected_pivot_col,
            stats.fast_rejected_row,
            stats.fast_rejected_refresh_state,
            stats.fast_rejected_tail_repair_ready,
            stats.fast_block_restarts,
            stats.fast_tail_restarts,
            stats.fast_kls_block_restart_last_row_pipeline,
            stats.fast_kls_block_restart_row_pipeline_count,
            stats.fast_kls_block_restart_last_row_pipeline_rows,
            stats.fast_kls_block_restart_last_row_pipeline_threads,
            stats.fast_kls_block_restart_last_row_pipeline_prefix_rows,
            stats.fast_kls_block_restart_last_row_pipeline_suffix_rows,
            stats.fast_kls_block_restart_last_row_pipeline_etree_tail,
            stats.fast_kls_block_restart_row_pipeline_etree_tail_count,
            stats.fast_kls_block_restart_last_row_pipeline_etree_tail_rows,
            stats
              .fast_kls_block_restart_last_row_pipeline_etree_tail_exact_mask,
            stats.fast_kls_block_restart_last_row_pipeline_etree_ready,
            stats.fast_kls_block_restart_row_pipeline_etree_ready_count,
            stats.fast_kls_block_restart_last_row_pipeline_etree_ready_rows,
            stats.fast_kls_block_restart_last_row_pipeline_etree_ready_threads,
            stats.fast_rejected_pivoting_tail_columns,
            stats.fast_rejected_pivoting_tail_seed_columns,
            stats.fast_rejected_pivoting_tail_row_seed_columns,
            stats.fast_rejected_row_tail_columns,
            stats.fast_rejected_pivoting_tail_topological);
    ok = 0;
  }
  if (ok && !require_pivoting_tail_plan(&stats,
                                        "parallel checked-row restart")) {
    ok = 0;
  }
  if (ok && !require_tail_overcompute_bounds(
              &stats, "parallel checked-row restart")) {
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
  if (ok && !require_tail_overcompute_bounds(&stats,
                                             "tail-prefix validation")) {
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
             stats.fast_rejected_pivoting_tail_block_seed_columns !=
               stats.fast_rejected_suffix_columns ||
             stats.fast_block_restarts != 1 ||
             stats.fast_tail_restarts != 1)) {
    fprintf(stderr,
            "unexpected mapped prefix-tail restart stats: pivot=%" PRId64
            ", start=%" PRId64 ", size=%" PRId64 ", refresh=%d"
            ", tail_candidate=%" PRId64 ", tail_candidate_count=%" PRId64
            ", tail_repair_ready=%d, repaired_match=%d"
            ", prefix_changed=%" PRId64 ", suffix_changed=%" PRId64
            ", tail_ready=%d, tail_cols=%" PRId64 ", tail_work=%.6g"
            ", seed=%" PRId64 "/%" PRId64 "/%" PRId64
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
            stats.fast_rejected_pivoting_tail_seed_columns,
            stats.fast_rejected_pivoting_tail_row_seed_columns,
            stats.fast_rejected_pivoting_tail_block_seed_columns,
            stats.fast_block_restarts,
            stats.fast_tail_restarts);
    ok = 0;
  }
  if (ok && !require_pivoting_tail_plan(&stats,
                                        "mapped prefix-tail restart")) {
    ok = 0;
  }
  if (ok && !require_tail_overcompute_bounds(
              &stats, "mapped prefix-tail restart")) {
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

static int test_nonroot_tail_refreshes_preserved_suffix(void) {
  const int32_t ap[] = {0, 1, 3, 5, 6};
  const int32_t ai[] = {0, 1, 2, 1, 2, 3};
  const double ax0[] = {2.0, 2.0, 1.0, 1.0, 2.0, 4.0};
  const double ax1[] = {2.0, 1.0e-12, 1.0, 1.0, 2.0, 8.0};
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
  if (!require_ok(kls_create(&solver), "create nonroot suffix refresh")) {
    ok = 0;
  }
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 4, ap, ai, 0,
                                        &options),
                        "analyze nonroot suffix refresh")) {
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor nonroot suffix refresh base")) {
    ok = 0;
  }
  if (ok && unsetenv("KLS_ENABLE_CHECKED_ROW_REFACTOR") != 0) {
    perror("unsetenv KLS_ENABLE_CHECKED_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(solver, ax1),
                        "factor nonroot suffix refresh repair")) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_CHECKED_ROW_REFACTOR", had_saved_env,
                         saved_env)) {
    ok = 0;
  }
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve nonroot suffix refresh")) {
    ok = 0;
  }

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats nonroot suffix refresh")) {
    ok = 0;
  }
  if (ok && (stats.fast_rejected_pivot != 1 ||
             stats.fast_rejected_block_start != 0 ||
             stats.fast_rejected_block_size != 4 ||
             stats.fast_rejected_refresh_state !=
               KLS_FAST_REJECT_REFRESH_PREFIX ||
             stats.fast_rejected_pivoting_tail_columns != 2 ||
             stats.fast_rejected_pivoting_tail_first != 1 ||
             stats.fast_rejected_pivoting_tail_last != 2 ||
             stats.fast_rejected_pivoting_tail_suffix_exact != 0 ||
             stats.fast_rejected_pivoting_tail_suffix_overcompute_columns != 1 ||
             stats.fast_repaired_tail_restart_ready != 1 ||
             stats.fast_repaired_tail_restart_columns != 2 ||
             stats.fast_repaired_tail_restart_exact_mask != 1 ||
             stats.fast_repaired_tail_restart_etree_mask != 1 ||
             stats.fast_repaired_tail_restart_skipped_columns != 0 ||
             stats.fast_tail_restarts != 1)) {
    fprintf(stderr,
            "unexpected nonroot suffix-refresh stats: pivot=%" PRId64
            ", block=[%" PRId64 ",%" PRId64 "), refresh=%d"
            ", tail=%" PRId64 "/%" PRId64 "/%" PRId64
            ", suffix_exact=%d, suffix_over=%" PRId64
            ", ready=%d, repaired=%" PRId64 ", exact_mask=%d"
            ", etree_mask=%d, skipped=%" PRId64 ", tail_restarts=%d\n",
            stats.fast_rejected_pivot,
            stats.fast_rejected_block_start,
            stats.fast_rejected_block_start + stats.fast_rejected_block_size,
            stats.fast_rejected_refresh_state,
            stats.fast_rejected_pivoting_tail_columns,
            stats.fast_rejected_pivoting_tail_first,
            stats.fast_rejected_pivoting_tail_last,
            stats.fast_rejected_pivoting_tail_suffix_exact,
            stats.fast_rejected_pivoting_tail_suffix_overcompute_columns,
            stats.fast_repaired_tail_restart_ready,
            stats.fast_repaired_tail_restart_columns,
            stats.fast_repaired_tail_restart_exact_mask,
            stats.fast_repaired_tail_restart_etree_mask,
            stats.fast_repaired_tail_restart_skipped_columns,
            stats.fast_tail_restarts);
    ok = 0;
  }
  if (ok && !require_pivoting_tail_plan(&stats, "nonroot suffix refresh")) {
    ok = 0;
  }
  if (ok && !require_tail_overcompute_bounds(&stats,
                                             "nonroot suffix refresh")) {
    ok = 0;
  }
  for (int32_t i = 0; ok && i < 4; ++i) {
    if (!close_enough(x[i], expected[i])) {
      fprintf(stderr,
              "unexpected nonroot suffix-refresh solution at %d:"
              " %.17g != %.17g\n",
              (int)i, x[i], expected[i]);
      ok = 0;
    }
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
  if (ok && !require_tail_overcompute_bounds(&stats,
                                             "scaled block restart")) {
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
             stats.fast_tail_restarts != 1 ||
             stats.fast_repaired_last_offdiag_suffix_refresh != 1 ||
             stats.fast_repaired_offdiag_suffix_refresh_count != 1 ||
             stats.fast_repaired_offdiag_full_refresh_count != 0)) {
    fprintf(stderr,
            "unexpected scaled prefix-tail stats: pivot=%" PRId64
            ", start=%" PRId64 ", size=%" PRId64 ", refresh=%d"
            ", prefix_changed=%" PRId64 ", suffix_changed=%" PRId64
            ", tail_ready=%d, tail_cols=%" PRId64 ", tail_work=%.6g"
            ", block_restarts=%d, tail_restarts=%d"
            ", offdiag_suffix=%d, offdiag_suffix_count=%" PRId64
            ", offdiag_full_count=%" PRId64 "\n",
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
            stats.fast_tail_restarts,
            stats.fast_repaired_last_offdiag_suffix_refresh,
            stats.fast_repaired_offdiag_suffix_refresh_count,
            stats.fast_repaired_offdiag_full_refresh_count);
    ok = 0;
  }
  if (ok && !require_pivoting_tail_plan(&stats,
                                        "scaled prefix-tail restart")) {
    ok = 0;
  }
  if (ok && !require_tail_overcompute_bounds(
              &stats, "scaled prefix-tail restart")) {
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
  if (ok && !require_tail_overcompute_bounds(&stats, "btf block restart")) {
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

static int test_btf_checked_row_tail_scope(void) {
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
                        "analyze btf checked-row tail")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor btf checked-row base")) ok = 0;
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(solver, ax1),
                        "factor btf checked-row repair")) ok = 0;
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
                        "solve btf checked-row tail")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats btf checked-row tail")) {
    ok = 0;
  }
  if (ok && (stats.nblocks < 2 ||
             stats.fast_rejected_block_start < 0 ||
             stats.fast_rejected_pivot < stats.fast_rejected_block_start ||
             stats.fast_rejected_pivot >=
               stats.fast_rejected_block_start +
                 stats.fast_rejected_block_size ||
             stats.fast_rejected_row < stats.fast_rejected_block_start ||
             stats.fast_rejected_row >=
               stats.fast_rejected_block_start +
                 stats.fast_rejected_block_size ||
             stats.fast_rejected_block_size != 2 ||
             stats.fast_rejected_refresh_state !=
               KLS_FAST_REJECT_REFRESH_PREFIX ||
             stats.fast_rejected_row_tail_columns < 1 ||
             stats.fast_rejected_row_tail_columns >
               stats.fast_rejected_suffix_columns ||
             stats.fast_rejected_row_tail_work <= 0.0 ||
             stats.fast_rejected_tail_candidate_row <=
               stats.fast_rejected_pivot ||
             stats.fast_rejected_tail_candidate_row >=
               stats.fast_rejected_block_start +
                 stats.fast_rejected_block_size ||
             stats.fast_rejected_tail_candidate_count < 1 ||
             stats.fast_rejected_tail_candidate_position < 0 ||
             stats.fast_rejected_tail_repair_ready != 1)) {
    fprintf(stderr,
            "unexpected btf checked-row tail stats: nblocks=%" PRId64
            ", pivot=%" PRId64 ", row=%" PRId64
            ", block=[%" PRId64 ",%" PRId64 "), refresh=%d"
            ", row_tail=%" PRId64 ", row_work=%.6g"
            ", row_seed=%" PRId64
            ", tail_row=%" PRId64 ", tail_count=%" PRId64
            ", tail_pos=%" PRId64 ", tail_ready=%d\n",
            stats.nblocks,
            stats.fast_rejected_pivot,
            stats.fast_rejected_row,
            stats.fast_rejected_block_start,
            stats.fast_rejected_block_start + stats.fast_rejected_block_size,
            stats.fast_rejected_refresh_state,
            stats.fast_rejected_row_tail_columns,
            stats.fast_rejected_row_tail_work,
            stats.fast_rejected_pivoting_tail_row_seed_columns,
            stats.fast_rejected_tail_candidate_row,
            stats.fast_rejected_tail_candidate_count,
            stats.fast_rejected_tail_candidate_position,
            stats.fast_rejected_tail_repair_ready);
    ok = 0;
  }
  if (ok && !require_pivoting_tail_plan(&stats, "btf checked-row tail")) {
    ok = 0;
  }
  if (ok && (!close_enough(x[0], 1.0) || !close_enough(x[1], 2.0) ||
             !close_enough(x[2], 3.0) || !close_enough(x[3], 4.0))) {
    fprintf(stderr,
            "unexpected btf checked-row solution: %.17g %.17g %.17g %.17g\n",
            x[0], x[1], x[2], x[3]);
    ok = 0;
  }

  kls_destroy(solver);
  free(saved_env);
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
             stats.fast_repaired_tail_restart_work <= 0.0 ||
             stats.fast_repaired_last_offdiag_suffix_refresh != 1 ||
             stats.fast_repaired_offdiag_suffix_refresh_count != 1 ||
             stats.fast_repaired_offdiag_full_refresh_count != 0)) {
    fprintf(stderr,
            "unexpected btf offblock tail stats: nblocks=%" PRId64
            ", pivot=%" PRId64 ", start=%" PRId64 ", size=%" PRId64
            ", refresh=%d, restarts=%d, tail_restarts=%d"
            ", tail_ready=%d, tail_cols=%" PRId64 ", tail_work=%.6g"
            ", offdiag_suffix=%d, offdiag_suffix_count=%" PRId64
            ", offdiag_full_count=%" PRId64 "\n",
            stats.nblocks,
            stats.fast_rejected_pivot,
            stats.fast_rejected_block_start,
            stats.fast_rejected_block_size,
            stats.fast_rejected_refresh_state,
            stats.fast_block_restarts,
            stats.fast_tail_restarts,
            stats.fast_repaired_tail_restart_ready,
            stats.fast_repaired_tail_restart_columns,
            stats.fast_repaired_tail_restart_work,
            stats.fast_repaired_last_offdiag_suffix_refresh,
            stats.fast_repaired_offdiag_suffix_refresh_count,
            stats.fast_repaired_offdiag_full_refresh_count);
    ok = 0;
  }
  if (ok && !require_pivoting_tail_plan(&stats, "btf offblock tail")) {
    ok = 0;
  }
  if (ok && !require_tail_overcompute_bounds(&stats, "btf offblock tail")) {
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

static int test_parallel_btf_suffix_after_prefix_tail_restart(void) {
  const int32_t ap[] = {0, 1, 2, 3, 9, 15, 21, 23, 25, 27};
  const int32_t ai[] = {
    0,
    1,
    2,
    0, 1, 2, 3, 4, 5,
    0, 1, 2, 3, 4, 5,
    0, 1, 2, 3, 4, 5,
    3, 6,
    4, 7,
    5, 8
  };
  const double ax0[] = {
    4.0,
    5.0,
    6.0,
    0.1, -0.1, 0.05, 8.0, 0.01, 0.02,
    0.02, 0.1, -0.04, 0.01, 8.0, 0.03,
    -0.03, 0.02, 0.06, 0.02, 0.03, 8.0,
    0.2, 7.0,
    -0.3, 7.5,
    0.4, 8.5
  };
  const double ax1[] = {
    4.0,
    5.0,
    6.0,
    0.1, -0.1, 0.05, 8.0, 0.01, 0.02,
    0.02, 0.1, -0.04, 0.01, 1.0e-12, 2.0,
    -0.03, 0.02, 0.06, 0.02, 0.03, 8.0,
    0.25, 7.25,
    0.1, 7.75,
    -0.2, 8.75
  };
  const double expected[] = {
    1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0
  };
  double b[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double x[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  for (int32_t col = 0; col < 9; ++col) {
    for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
      b[ai[p]] += ax1[p] * expected[col];
    }
  }

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.threads = 2;
  options.ordering = KLS_ORDERING_NATURAL;
  options.scale = -1;
  options.pivot_tolerance = 0.001;
  options.static_pivoting = 0;

  int ok = 1;
  if (!require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 9, ap, ai, 0,
                                        &options),
                        "analyze parallel btf suffix tail")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor parallel btf suffix base")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax1),
                        "factor parallel btf suffix repair")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve parallel btf suffix")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats parallel btf suffix")) {
    ok = 0;
  }
  if (ok && (stats.nblocks < 4 ||
             stats.fast_rejected_refresh_state !=
               KLS_FAST_REJECT_REFRESH_PREFIX ||
             stats.fast_block_restarts != 1 ||
             stats.fast_tail_restarts != 0 ||
             stats.fast_kls_block_restart_last_row_pipeline != 1 ||
             stats.fast_kls_block_restart_row_pipeline_count < 1 ||
             stats.fast_kls_block_restart_last_row_pipeline_prefix_rows <= 0 ||
             stats.fast_kls_block_restart_last_row_pipeline_prefix_rows +
                 stats.fast_kls_block_restart_last_row_pipeline_rows !=
               stats.fast_rejected_block_size ||
             stats.fast_repaired_last_offdiag_suffix_refresh != 1 ||
             stats.fast_repaired_offdiag_suffix_refresh_count != 1 ||
             stats.fast_repaired_offdiag_full_refresh_count != 0 ||
             stats.fast_repaired_parallel_tail_blocks < 2)) {
    fprintf(stderr,
            "unexpected parallel btf suffix stats: nblocks=%" PRId64
            ", pivot=%" PRId64 ", start=%" PRId64 ", size=%" PRId64
            ", refresh=%d, restarts=%d, tail_restarts=%d"
            ", tail_ready=%d, pipeline=%d/%" PRId64
            ", pipeline_rows=%" PRId64 ", pipeline_threads=%" PRId64
            ", pipeline_prefix=%" PRId64
            ", offdiag_suffix=%d, offdiag_suffix_count=%" PRId64
            ", offdiag_full_count=%" PRId64
            ", parallel_tail_blocks=%" PRId64 "\n",
            stats.nblocks,
            stats.fast_rejected_pivot,
            stats.fast_rejected_block_start,
            stats.fast_rejected_block_size,
            stats.fast_rejected_refresh_state,
            stats.fast_block_restarts,
            stats.fast_tail_restarts,
            stats.fast_repaired_tail_restart_ready,
            stats.fast_kls_block_restart_last_row_pipeline,
            stats.fast_kls_block_restart_row_pipeline_count,
            stats.fast_kls_block_restart_last_row_pipeline_rows,
            stats.fast_kls_block_restart_last_row_pipeline_threads,
            stats.fast_kls_block_restart_last_row_pipeline_prefix_rows,
            stats.fast_repaired_last_offdiag_suffix_refresh,
            stats.fast_repaired_offdiag_suffix_refresh_count,
            stats.fast_repaired_offdiag_full_refresh_count,
            stats.fast_repaired_parallel_tail_blocks);
    ok = 0;
  }
  if (ok && !require_pivoting_tail_plan(&stats, "parallel btf suffix")) {
    ok = 0;
  }
  if (ok && !require_tail_overcompute_bounds(&stats,
                                             "parallel btf suffix")) {
    ok = 0;
  }
  for (int32_t i = 0; ok && i < 9; ++i) {
    if (!close_enough(x[i], expected[i])) {
      fprintf(stderr,
              "unexpected parallel btf suffix solution at %d: %.17g != %.17g\n",
              (int)i, x[i], expected[i]);
      ok = 0;
    }
  }

  kls_destroy(solver);
  return ok;
}

static int run_btf_row_refactor_offblock_refresh(int scale) {
  const int32_t ap[] = {0, 2, 4, 7, 9};
  const int32_t ai[] = {0, 1, 0, 1, 0, 2, 3, 2, 3};
  const double ax0[] = {
    2.0, 1.0, 1.0, 2.0,
    0.5, 2.0, 1.0, 1.0, 2.0
  };
  const double ax1[] = {
    2.5, 1.2, 0.8, 2.25,
    1.75, 3.0, 0.5, 0.75, 2.75
  };
  const int32_t nrhs = 5;
  const double expected_rhs[5][4] = {
    {1.0, 2.0, 3.0, 4.0},
    {-2.0, 0.5, 1.25, -0.75},
    {0.0, 1.0, -1.0, 2.0},
    {3.0, -3.0, 0.25, 0.5},
    {-0.5, -1.5, 2.5, 1.0}
  };
  const double expected_t_rhs[5][4] = {
    {0.75, -1.25, 1.5, 2.25},
    {2.0, 1.0, 0.0, -1.0},
    {-1.5, 0.25, 0.5, 3.0},
    {1.0, -2.0, 2.0, -0.5},
    {0.0, 0.5, -1.5, 1.25}
  };
  double b[20] = {0.0};
  double bt[20] = {0.0};
  double x[20] = {0.0};
  double xt[20] = {0.0};

  for (int32_t rhs = 0; rhs < nrhs; ++rhs) {
    for (int32_t col = 0; col < 4; ++col) {
      for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
        b[rhs * 4 + ai[p]] += ax1[p] * expected_rhs[rhs][col];
        bt[rhs * 4 + col] += ax1[p] * expected_t_rhs[rhs][ai[p]];
      }
    }
  }

  const char *saved_env_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_env = saved_env_value != NULL ? strdup(saved_env_value) : NULL;
  const int had_saved_env = saved_env_value != NULL;

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.threads = 1;
  options.ordering = KLS_ORDERING_NATURAL;
  options.orientation = KLS_ORIENTATION_NORMAL;
  options.use_btf = 1;
  options.scale = scale;
  options.static_pivoting = 0;

  int ok = 1;
  if (had_saved_env && saved_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_ROW_REFACTOR\n");
    ok = 0;
  }
  if (!require_ok(kls_create(&solver), "create btf row refactor")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 4, ap, ai, 0,
                                        &options),
                        "analyze btf row refactor")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor btf row refactor base")) ok = 0;
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && !require_ok(kls_refactor(solver, ax1),
                        "btf row refactor")) ok = 0;

  kls_stats refactor_stats;
  refactor_stats.struct_size = sizeof(refactor_stats);
  if (ok && !require_ok(kls_get_stats(solver, &refactor_stats),
                        "stats btf row refactor")) {
    ok = 0;
  }
  if (ok && (refactor_stats.nblocks < 2 ||
             refactor_stats.selected_scale != scale ||
             refactor_stats.row_refactor_last_run != 1 ||
             refactor_stats.row_refactor_last_checked != 0 ||
             refactor_stats.row_refactor_last_parallel != 0 ||
             refactor_stats.row_refactor_values_dirty != 1 ||
             refactor_stats.row_refactor_last_lazy_value_scatter != 1 ||
             refactor_stats.row_refactor_last_row_solve != 0 ||
             refactor_stats.row_refactor_run_count != 1 ||
             refactor_stats.row_refactor_lazy_value_scatter_run_count != 1 ||
             refactor_stats.row_refactor_row_solve_run_count != 0)) {
    fprintf(stderr,
            "unexpected btf row-refactor stats: nblocks=%" PRId64
            ", scale=%d, last=%d/%d/%d, dirty=%d, lazy=%d/%" PRId64
            ", row_solve=%d/%" PRId64 ", runs=%" PRId64 "\n",
            refactor_stats.nblocks,
            refactor_stats.selected_scale,
            refactor_stats.row_refactor_last_run,
            refactor_stats.row_refactor_last_checked,
            refactor_stats.row_refactor_last_parallel,
            refactor_stats.row_refactor_values_dirty,
            refactor_stats.row_refactor_last_lazy_value_scatter,
            refactor_stats.row_refactor_lazy_value_scatter_run_count,
            refactor_stats.row_refactor_last_row_solve,
            refactor_stats.row_refactor_row_solve_run_count,
            refactor_stats.row_refactor_run_count);
    ok = 0;
  }
  if (ok && !require_ok(kls_solve(solver, nrhs, b, 4, x, 4),
                        "solve btf row refactor")) ok = 0;

  kls_stats solve_stats;
  solve_stats.struct_size = sizeof(solve_stats);
  if (ok && !require_ok(kls_get_stats(solver, &solve_stats),
                        "solve stats btf row refactor")) {
    ok = 0;
  }
  if (ok && (solve_stats.row_refactor_values_dirty != 1 ||
             solve_stats.row_refactor_last_lazy_value_scatter != 1 ||
             solve_stats.row_refactor_last_row_solve != 1 ||
             solve_stats.row_refactor_lazy_value_scatter_run_count != 1 ||
             solve_stats.row_refactor_row_solve_run_count != 1)) {
    fprintf(stderr,
            "unexpected btf row-refactor solve stats for scale %d: dirty=%d"
            ", lazy=%d/%" PRId64 ", row_solve=%d/%" PRId64 "\n",
            scale,
            solve_stats.row_refactor_values_dirty,
            solve_stats.row_refactor_last_lazy_value_scatter,
            solve_stats.row_refactor_lazy_value_scatter_run_count,
            solve_stats.row_refactor_last_row_solve,
            solve_stats.row_refactor_row_solve_run_count);
    ok = 0;
  }
  for (int32_t rhs = 0; ok && rhs < nrhs; ++rhs) {
    for (int32_t i = 0; ok && i < 4; ++i) {
      if (!close_enough(x[rhs * 4 + i], expected_rhs[rhs][i])) {
        fprintf(stderr,
                "unexpected btf row-refactor solution for scale %d"
                " rhs %d at %d: %.17g != %.17g\n",
                scale, (int)rhs, (int)i, x[rhs * 4 + i],
                expected_rhs[rhs][i]);
        ok = 0;
      }
    }
  }
  if (ok && !require_ok(kls_solve_transpose(solver, nrhs, bt, 4, xt, 4),
                        "transpose solve btf row refactor")) ok = 0;
  if (ok && !require_ok(kls_get_stats(solver, &solve_stats),
                        "transpose stats btf row refactor")) {
    ok = 0;
  }
  if (ok && (solve_stats.row_refactor_values_dirty != 1 ||
             solve_stats.row_refactor_last_lazy_value_scatter != 1 ||
             solve_stats.row_refactor_last_row_solve != 1 ||
             solve_stats.row_refactor_lazy_value_scatter_run_count != 1 ||
             solve_stats.row_refactor_row_solve_run_count != 2)) {
    fprintf(stderr,
            "unexpected btf row-refactor transpose stats for scale %d:"
            " dirty=%d"
            ", lazy=%d/%" PRId64 ", row_solve=%d/%" PRId64 "\n",
            scale,
            solve_stats.row_refactor_values_dirty,
            solve_stats.row_refactor_last_lazy_value_scatter,
            solve_stats.row_refactor_lazy_value_scatter_run_count,
            solve_stats.row_refactor_last_row_solve,
            solve_stats.row_refactor_row_solve_run_count);
    ok = 0;
  }
  for (int32_t rhs = 0; ok && rhs < nrhs; ++rhs) {
    for (int32_t i = 0; ok && i < 4; ++i) {
      if (!close_enough(xt[rhs * 4 + i], expected_t_rhs[rhs][i])) {
        fprintf(stderr,
                "unexpected btf row-refactor transpose solution for scale %d"
                " rhs %d at %d: %.17g != %.17g\n",
                scale, (int)rhs, (int)i, xt[rhs * 4 + i],
                expected_t_rhs[rhs][i]);
        ok = 0;
      }
    }
  }

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
  kls_destroy(solver);
  free(saved_env);
  return ok;
}

static int test_btf_row_refactor_offblock_refresh(void) {
  return run_btf_row_refactor_offblock_refresh(-1) &&
         run_btf_row_refactor_offblock_refresh(2);
}

static int test_verified_rhs_reuse_contract(void) {
  const int32_t ap[] = {0, 2, 4, 7, 9};
  const int32_t ai[] = {0, 1, 0, 1, 0, 2, 3, 2, 3};
  const double ax0[] = {
    2.0, 1.0, 1.0, 2.0,
    0.5, 2.0, 1.0, 1.0, 2.0
  };
  const double ax1[] = {
    2.5, 1.2, 0.8, 2.25,
    1.75, 3.0, 0.5, 0.75, 2.75
  };
  const double ax2[] = {
    2.75, 1.0, 0.625, 2.5,
    1.5, 3.25, 0.375, 0.5, 3.0
  };
  const double expected1[] = {1.0, -2.0, 0.5, 3.0};
  const double expected2[] = {-0.75, 1.25, 2.0, -1.5};
  const char *env_names[] = {
    "KLS_ENABLE_ROW_REFACTOR",
    "KLS_DISABLE_VERIFIED_RHS_REUSE"
  };
  enum { ENV_COUNT = (int)(sizeof(env_names) / sizeof(env_names[0])) };
  char *saved_env[ENV_COUNT];
  int had_env[ENV_COUNT];
  int ok = 1;
  for (int i = 0; i < ENV_COUNT; ++i) {
    const char *value = getenv(env_names[i]);
    had_env[i] = value != NULL;
    saved_env[i] = value != NULL ? strdup(value) : NULL;
    if ((value != NULL && saved_env[i] == NULL) ||
        unsetenv(env_names[i]) != 0) {
      ok = 0;
    }
  }

  kls_options options;
  kls_default_options(&options);
  options.threads = 1;
  options.ordering = KLS_ORDERING_NATURAL;
  options.orientation = KLS_ORIENTATION_NORMAL;
  options.scale = -1;
  options.static_pivoting = 0;

  kls_solver *solver = NULL;
  if (ok && !require_ok(kls_create(&solver),
                        "create verified RHS reuse")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(
                          solver, KLS_INDEX_INT32, 4, ap, ai, 0, &options),
                        "analyze verified RHS reuse")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor verified RHS reuse")) ok = 0;
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv verified RHS row refactor");
    ok = 0;
  }
  if (ok && !require_ok(kls_refactor(solver, ax1),
                        "refactor verified RHS reuse")) ok = 0;

  double b[4] = {0.0, 0.0, 0.0, 0.0};
  double x[4] = {0.0, 0.0, 0.0, 0.0};
  for (int32_t col = 0; col < 4; ++col) {
    for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
      b[ai[p]] += ax1[p] * expected1[col];
    }
  }
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "first verified RHS solve")) ok = 0;
  kls_stats stats;
  memset(&stats, 0, sizeof(stats));
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "first verified RHS stats")) ok = 0;
  const int64_t first_row_solve_count =
    stats.row_refactor_row_solve_run_count;
  if (ok && (stats.row_refactor_last_row_solve != 1 ||
             first_row_solve_count < 1 || stats.verified_rhs_reused ||
             stats.verified_rhs_reuse_count != 0)) {
    fprintf(stderr,
            "first verified RHS did not use the checked row solve: %d/%" PRId64
            "\n", stats.row_refactor_last_row_solve,
            first_row_solve_count);
    ok = 0;
  }

  memset(x, 0, sizeof(x));
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "reused verified RHS solve")) ok = 0;
  memset(&stats, 0, sizeof(stats));
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "reused verified RHS stats")) ok = 0;
  if (ok && (!stats.verified_rhs_reused ||
             stats.verified_rhs_reuse_count != 1)) {
    fprintf(stderr,
            "identical verified RHS was not reused: reused=%d/%" PRId64
            ", row=%d/%" PRId64 "\n",
            stats.verified_rhs_reused,
            stats.verified_rhs_reuse_count,
            stats.row_refactor_last_row_solve,
            stats.row_refactor_row_solve_run_count);
    ok = 0;
  }
  const int64_t reused_row_solve_count =
    stats.row_refactor_row_solve_run_count;
  if (ok && reused_row_solve_count != first_row_solve_count + 1) {
    fprintf(stderr,
            "verified RHS reuse changed solve engines: %" PRId64
            " -> %" PRId64 "\n",
            first_row_solve_count, reused_row_solve_count);
    ok = 0;
  }
  for (int32_t i = 0; ok && i < 4; ++i) {
    if (!close_enough(x[i], expected1[i])) {
      fprintf(stderr,
              "unexpected reused RHS solution at %d: %.17g != %.17g\n",
              (int)i, x[i], expected1[i]);
      ok = 0;
    }
  }

  int64_t checked_row_solve_count = reused_row_solve_count;
  for (int env_i = 1; ok && env_i < ENV_COUNT; ++env_i) {
    if (setenv(env_names[env_i], "1", 1) != 0) {
      perror("setenv verified RHS disable");
      ok = 0;
      break;
    }
    memset(x, 0, sizeof(x));
    if (!require_ok(kls_solve(solver, 1, b, 0, x, 0),
                    "disabled verified RHS solve")) {
      ok = 0;
    }
    memset(&stats, 0, sizeof(stats));
    stats.struct_size = sizeof(stats);
    if (ok && !require_ok(kls_get_stats(solver, &stats),
                          "disabled verified RHS stats")) ok = 0;
    checked_row_solve_count++;
    if (ok && (stats.verified_rhs_reused ||
               stats.verified_rhs_reuse_count != 1 ||
               stats.row_refactor_last_row_solve != 1 ||
               stats.row_refactor_row_solve_run_count !=
                 checked_row_solve_count)) {
      fprintf(stderr, "%s did not disable verified RHS reuse\n",
              env_names[env_i]);
      ok = 0;
    }
    for (int32_t i = 0; ok && i < 4; ++i) {
      if (!close_enough(x[i], expected1[i])) {
        fprintf(stderr,
                "unexpected disabled-reuse solution at %d: %.17g != %.17g\n",
                (int)i, x[i], expected1[i]);
        ok = 0;
      }
    }
    if (unsetenv(env_names[env_i]) != 0) {
      perror("unsetenv verified RHS disable");
      ok = 0;
    }
  }

  memset(b, 0, sizeof(b));
  memset(x, 0, sizeof(x));
  for (int32_t col = 0; col < 4; ++col) {
    for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
      b[ai[p]] += ax1[p] * expected2[col];
    }
  }
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "changed verified RHS solve")) ok = 0;
  memset(&stats, 0, sizeof(stats));
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "changed verified RHS stats")) ok = 0;
  if (ok && (stats.verified_rhs_reused ||
             stats.verified_rhs_reuse_count != 1 ||
             stats.row_refactor_last_row_solve != 1 ||
             stats.row_refactor_row_solve_run_count !=
               checked_row_solve_count + 1)) {
    fprintf(stderr,
            "changed RHS incorrectly reused a verdict: %d/%" PRId64 "\n",
            stats.row_refactor_last_row_solve,
            stats.row_refactor_row_solve_run_count);
    ok = 0;
  }
  for (int32_t i = 0; ok && i < 4; ++i) {
    if (!close_enough(x[i], expected2[i])) {
      fprintf(stderr,
              "unexpected changed RHS solution at %d: %.17g != %.17g\n",
              (int)i, x[i], expected2[i]);
      ok = 0;
    }
  }

  if (ok && !require_ok(kls_refactor(solver, ax2),
                        "changed factor verified RHS reuse")) ok = 0;
  memset(b, 0, sizeof(b));
  memset(x, 0, sizeof(x));
  for (int32_t col = 0; col < 4; ++col) {
    for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
      b[ai[p]] += ax2[p] * expected1[col];
    }
  }
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "invalidated verified RHS solve")) ok = 0;
  memset(&stats, 0, sizeof(stats));
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "invalidated verified RHS stats")) ok = 0;
  if (ok && (stats.verified_rhs_reused ||
             stats.verified_rhs_reuse_count != 1 ||
             stats.row_refactor_last_row_solve != 1)) {
    fprintf(stderr, "changed factor incorrectly retained its RHS verdict\n");
    ok = 0;
  }
  for (int32_t i = 0; ok && i < 4; ++i) {
    if (!close_enough(x[i], expected1[i])) {
      fprintf(stderr,
              "unexpected invalidated RHS solution at %d: %.17g != %.17g\n",
              (int)i, x[i], expected1[i]);
      ok = 0;
    }
  }

  kls_destroy(solver);
  for (int i = 0; i < ENV_COUNT; ++i) {
    if (!(had_env[i] && saved_env[i] == NULL) &&
        !restore_env_value(env_names[i], had_env[i], saved_env[i])) {
      ok = 0;
    }
    free(saved_env[i]);
  }
  return ok;
}

static int test_btf_singleton_row_refactor_pattern(void) {
  const int32_t ap[] = {0, 1, 3, 5, 6};
  const int32_t ai[] = {0, 1, 2, 1, 2, 3};
  const double ax0[] = {2.0, 4.0, 1.0, 1.0, 3.0, 5.0};
  const double ax1[] = {2.5, 4.5, 0.75, 1.25, 3.5, 5.5};
  const double expected[] = {1.0, -2.0, 0.5, 3.0};
  double b[4] = {0.0, 0.0, 0.0, 0.0};
  double x[4] = {0.0, 0.0, 0.0, 0.0};

  for (int32_t col = 0; col < 4; ++col) {
    for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
      b[ai[p]] += ax1[p] * expected[col];
    }
  }

  const char *saved_env_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_env = saved_env_value != NULL ? strdup(saved_env_value) : NULL;
  const int had_saved_env = saved_env_value != NULL;

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.threads = 1;
  options.ordering = KLS_ORDERING_NATURAL;
  options.orientation = KLS_ORIENTATION_NORMAL;
  options.use_btf = 1;
  options.scale = -1;
  options.static_pivoting = 0;

  int ok = 1;
  if (had_saved_env && saved_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_ROW_REFACTOR\n");
    ok = 0;
  }
  if (!require_ok(kls_create(&solver), "create btf singleton row refactor")) {
    ok = 0;
  }
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 4, ap, ai,
                                        0, &options),
                        "analyze btf singleton row refactor")) {
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor btf singleton row refactor base")) {
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && !require_ok(kls_refactor(solver, ax1),
                        "btf singleton row refactor")) {
    ok = 0;
  }

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats btf singleton row refactor")) {
    ok = 0;
  }
  if (ok && (stats.nblocks < 3 ||
             stats.row_refactor_last_run != 1 ||
             stats.row_refactor_last_checked != 0 ||
             stats.row_refactor_last_parallel != 0 ||
             stats.row_refactor_group_count < 3 ||
             stats.row_refactor_values_dirty != 1)) {
    fprintf(stderr,
            "unexpected btf singleton row-refactor stats: nblocks=%" PRId64
            ", last=%d/%d/%d, groups=%" PRId64 ", dirty=%d\n",
            stats.nblocks,
            stats.row_refactor_last_run,
            stats.row_refactor_last_checked,
            stats.row_refactor_last_parallel,
            stats.row_refactor_group_count,
            stats.row_refactor_values_dirty);
    ok = 0;
  }

  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve btf singleton row refactor")) {
    ok = 0;
  }
  for (int32_t i = 0; ok && i < 4; ++i) {
    if (!close_enough(x[i], expected[i])) {
      fprintf(stderr,
              "unexpected btf singleton row-refactor solution at %d:"
              " %.17g != %.17g\n",
              (int)i, x[i], expected[i]);
      ok = 0;
    }
  }

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
  kls_destroy(solver);
  free(saved_env);
  return ok;
}

static int test_btf_singleton_row_refactor_batch_groups(void) {
  const int32_t ap[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  const int32_t ai[] = {0, 1, 2, 3, 4, 5, 6, 7};
  const double ax0[] = {2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
  const double ax1[] = {3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  const double expected[] = {1.0, -2.0, 0.5, 3.0,
                             -1.5, 2.5, -0.25, 4.0};
  double b[8] = {0.0};
  double x[8] = {0.0};

  for (int32_t i = 0; i < 8; ++i) {
    b[i] = ax1[i] * expected[i];
  }

  const char *saved_env_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_env = saved_env_value != NULL ? strdup(saved_env_value) : NULL;
  const int had_saved_env = saved_env_value != NULL;

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.threads = 1;
  options.ordering = KLS_ORDERING_NATURAL;
  options.orientation = KLS_ORIENTATION_NORMAL;
  options.use_btf = 1;
  options.scale = -1;
  options.static_pivoting = 0;

  int ok = 1;
  if (had_saved_env && saved_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_ROW_REFACTOR\n");
    ok = 0;
  }
  if (!require_ok(kls_create(&solver), "create btf singleton batch")) {
    ok = 0;
  }
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 8, ap, ai,
                                        0, &options),
                        "analyze btf singleton batch")) {
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor btf singleton batch base")) {
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && !require_ok(kls_refactor(solver, ax1),
                        "btf singleton batch refactor")) {
    ok = 0;
  }

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats btf singleton batch")) {
    ok = 0;
  }
  if (ok && (stats.nblocks < 8 ||
             stats.row_refactor_last_run != 1 ||
             stats.row_refactor_group_count != 1 ||
             stats.row_refactor_segment_count != 0 ||
             stats.row_refactor_values_dirty != 1)) {
    fprintf(stderr,
            "unexpected btf singleton batch stats: nblocks=%" PRId64
            ", last=%d, groups=%" PRId64 ", segments=%" PRId64
            ", dirty=%d\n",
            stats.nblocks,
            stats.row_refactor_last_run,
            stats.row_refactor_group_count,
            stats.row_refactor_segment_count,
            stats.row_refactor_values_dirty);
    ok = 0;
  }

  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve btf singleton batch")) {
    ok = 0;
  }
  for (int32_t i = 0; ok && i < 8; ++i) {
    if (!close_enough(x[i], expected[i])) {
      fprintf(stderr,
              "unexpected btf singleton batch solution at %d:"
              " %.17g != %.17g\n",
              (int)i, x[i], expected[i]);
      ok = 0;
    }
  }

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
  kls_destroy(solver);
  free(saved_env);
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
  const char *saved_native_env_value =
    getenv("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR");
  char *saved_native_env =
    saved_native_env_value != NULL ? strdup(saved_native_env_value) : NULL;
  const int had_saved_native_env = saved_native_env_value != NULL;

  int ok = 1;
  if (had_saved_env && saved_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_CHECKED_ROW_REFACTOR\n");
    ok = 0;
  }
  if (had_saved_native_env && saved_native_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR\n");
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
  if (ok && setenv("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR", "auto", 1) != 0) {
    perror("setenv KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR=auto");
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(solver, ax1),
                        "factor dense checked-row repair")) ok = 0;
  if (!restore_env_value("KLS_ENABLE_CHECKED_ROW_REFACTOR",
                         had_saved_env, saved_env)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR",
                         had_saved_native_env, saved_native_env)) {
    ok = 0;
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
             stats.fast_repaired_tail_restart_ready != 1 ||
             stats.row_refactor_last_run != 1 ||
             stats.row_refactor_last_checked != 1 ||
             stats.row_refactor_last_parallel != 0 ||
             stats.row_refactor_last_defer_value_scatter != 1 ||
             stats.row_refactor_dense_segment_count < 1 ||
             stats.row_refactor_last_compact_dense_panel != 1 ||
             stats.row_refactor_compact_dense_panel_count < 1 ||
             stats.row_refactor_compact_dense_panel_eligible_count < 1 ||
             stats.row_refactor_compact_dense_panel_eligible_rows < 1 ||
             stats.row_refactor_compact_dense_panel_update_work <= 0.0 ||
             stats.row_refactor_compact_dense_panel_entries <= 0.0 ||
             stats.row_refactor_compact_dense_panel_persistent_groups < 1 ||
             stats.row_refactor_compact_dense_panel_persistent_entries < 1 ||
             stats.row_refactor_last_compact_dense_panel_persistent != 1 ||
             stats.row_refactor_compact_dense_panel_persistent_run_count < 1 ||
             stats.row_refactor_last_compact_dense_panel_direct_input_rows < 1 ||
             stats.row_refactor_compact_dense_panel_direct_input_rows < 1 ||
             stats.row_refactor_segment_input_target_rows < 1 ||
             stats.row_refactor_segment_input_target_entries < 1 ||
             stats.row_refactor_last_segment_target_input_rows <
               stats.row_refactor_last_compact_dense_panel_direct_input_rows ||
             stats.row_refactor_segment_target_input_rows <
               stats.row_refactor_compact_dense_panel_direct_input_rows)) {
    fprintf(stderr,
            "unexpected dense checked-row prefix stats: pivot=%" PRId64
            ", refresh=%d, block_restarts=%d, tail_restarts=%d"
            ", tail_ready=%d, row=%d/%d/%d, defer=%d"
            ", dense_segments=%" PRId64 ", compact=%d/%" PRId64
            ", eligible=%" PRId64 "/%" PRId64
            ", work=%.17g, entries=%.17g"
            ", persistent=%" PRId64 "/%" PRId64
            ", persistent_used=%d/%" PRId64
            ", direct_input=%" PRId64 "/%" PRId64
            ", target_map=%" PRId64 "/%" PRId64
            ", target_input=%" PRId64 "/%" PRId64 "\n",
            stats.fast_rejected_pivot,
            stats.fast_rejected_refresh_state,
            stats.fast_block_restarts,
            stats.fast_tail_restarts,
            stats.fast_repaired_tail_restart_ready,
            stats.row_refactor_last_run,
            stats.row_refactor_last_checked,
            stats.row_refactor_last_parallel,
            stats.row_refactor_last_defer_value_scatter,
            stats.row_refactor_dense_segment_count,
            stats.row_refactor_last_compact_dense_panel,
            stats.row_refactor_compact_dense_panel_count,
            stats.row_refactor_compact_dense_panel_eligible_count,
            stats.row_refactor_compact_dense_panel_eligible_rows,
            stats.row_refactor_compact_dense_panel_update_work,
            stats.row_refactor_compact_dense_panel_entries,
            stats.row_refactor_compact_dense_panel_persistent_groups,
            stats.row_refactor_compact_dense_panel_persistent_entries,
            stats.row_refactor_last_compact_dense_panel_persistent,
            stats.row_refactor_compact_dense_panel_persistent_run_count,
            stats.row_refactor_last_compact_dense_panel_direct_input_rows,
            stats.row_refactor_compact_dense_panel_direct_input_rows,
            stats.row_refactor_segment_input_target_rows,
            stats.row_refactor_segment_input_target_entries,
            stats.row_refactor_last_segment_target_input_rows,
            stats.row_refactor_segment_target_input_rows);
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
  free(saved_native_env);
  free(ap);
  free(ai);
  free(ax0);
  free(ax1);
  free(b);
  free(x);
  free(expected);
  return ok;
}

static int test_unchecked_row_dense_compact_panel(void) {
  const int32_t n = 80;
  const int32_t lead = 48;
  const int32_t lower = 16;
  const int32_t trail = n - lead - lower;
  const int32_t multi_rhs = 3;
  const size_t nnz =
    (size_t)lead * (size_t)(lead + lower) +
    (size_t)lower +
    (size_t)trail * ((size_t)lead + 1u);
  int32_t *ap = (int32_t *)malloc(((size_t)n + 1u) * sizeof(*ap));
  int32_t *ai = (int32_t *)malloc(nnz * sizeof(*ai));
  double *ax0 = (double *)malloc(nnz * sizeof(*ax0));
  double *ax1 = (double *)malloc(nnz * sizeof(*ax1));
  double *b = (double *)calloc((size_t)n, sizeof(*b));
  double *x = (double *)calloc((size_t)n, sizeof(*x));
  double *b_transpose = (double *)calloc((size_t)n, sizeof(*b_transpose));
  double *x_transpose = (double *)calloc((size_t)n, sizeof(*x_transpose));
  double *expected = (double *)malloc((size_t)n * sizeof(*expected));
  double *b_multi =
    (double *)calloc((size_t)n * (size_t)multi_rhs, sizeof(*b_multi));
  double *x_multi =
    (double *)calloc((size_t)n * (size_t)multi_rhs, sizeof(*x_multi));
  double *b_transpose_multi =
    (double *)calloc((size_t)n * (size_t)multi_rhs,
                     sizeof(*b_transpose_multi));
  double *x_transpose_multi =
    (double *)calloc((size_t)n * (size_t)multi_rhs,
                     sizeof(*x_transpose_multi));
  double *expected_multi =
    (double *)malloc((size_t)n * (size_t)multi_rhs *
                     sizeof(*expected_multi));
  if (ap == NULL || ai == NULL || ax0 == NULL || ax1 == NULL ||
      b == NULL || x == NULL || b_transpose == NULL ||
      x_transpose == NULL || expected == NULL ||
      b_multi == NULL || x_multi == NULL ||
      b_transpose_multi == NULL || x_transpose_multi == NULL ||
      expected_multi == NULL) {
    free(ap);
    free(ai);
    free(ax0);
    free(ax1);
    free(b);
    free(x);
    free(b_transpose);
    free(x_transpose);
    free(expected);
    free(b_multi);
    free(x_multi);
    free(b_transpose_multi);
    free(x_transpose_multi);
    free(expected_multi);
    return 0;
  }

  for (int32_t col = 0; col < n; ++col) {
    expected[col] = 0.75 + 0.03125 * (double)((5 * col) % 17);
  }
  size_t pos = 0;
  for (int32_t col = 0; col < n; ++col) {
    ap[col] = (int32_t)pos;
    if (col < lead) {
      for (int32_t row = 0; row < lead + lower; ++row) {
        const size_t p = pos++;
        ai[p] = row;
        if (row < lead) {
          ax0[p] = row == col
            ? 24.0 + 0.02 * (double)col
            : 0.0005 * (1.0 + (double)((row + 7 * col) % 11));
        } else {
          ax0[p] = 0.0003 * (1.0 + (double)((row + 3 * col) % 13));
        }
        ax1[p] = ax0[p] + (row == col
          ? 0.05 * (double)((col % 3) + 1)
          : 1.0e-5 * (double)(((row + col) % 5) - 2));
      }
    } else if (col < lead + lower) {
      const size_t p = pos++;
      ai[p] = col;
      ax0[p] = 18.0 + 0.03 * (double)col;
      ax1[p] = ax0[p] + 0.04 * (double)((col % 5) + 1);
    } else {
      for (int32_t row = 0; row < lead; ++row) {
        const size_t p = pos++;
        ai[p] = row;
        ax0[p] = 0.0002 * (1.0 + (double)((row + 5 * col) % 17));
        ax1[p] = ax0[p] +
          1.0e-5 * (double)(((row + 2 * col) % 7) - 3);
      }
      const size_t p = pos++;
      ai[p] = col;
      ax0[p] = 18.0 + 0.03 * (double)col;
      ax1[p] = ax0[p] + 0.04 * (double)((col % 5) + 1);
    }
  }
  ap[n] = (int32_t)pos;
  if (pos != nnz) {
    fprintf(stderr, "unexpected compact-supernode fixture nnz: %zu/%zu\n",
            pos, nnz);
    free(ap);
    free(ai);
    free(ax0);
    free(ax1);
    free(b);
    free(x);
    free(b_transpose);
    free(x_transpose);
    free(expected);
    free(b_multi);
    free(x_multi);
    free(b_transpose_multi);
    free(x_transpose_multi);
    free(expected_multi);
    return 0;
  }
  for (int32_t col = 0; col < n; ++col) {
    for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
      b[ai[p]] += ax1[p] * expected[col];
      b_transpose[col] += ax1[p] * expected[ai[p]];
    }
  }

  kls_solver *solver = NULL;
  kls_solver *disabled_solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.threads = 1;
  options.ordering = KLS_ORDERING_NATURAL;
  options.orientation = KLS_ORIENTATION_NORMAL;
  options.use_btf = 0;
  options.scale = -1;
  options.pivot_tolerance = 0.01;
  options.static_pivoting = 0;

  const char *saved_env_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_env = saved_env_value != NULL ? strdup(saved_env_value) : NULL;
  const int had_saved_env = saved_env_value != NULL;
  const char *saved_trsv_env_value =
    getenv("KLS_ENABLE_COMPACT_SUPERNODE_TRSV");
  char *saved_trsv_env =
    saved_trsv_env_value != NULL ? strdup(saved_trsv_env_value) : NULL;
  const int had_saved_trsv_env = saved_trsv_env_value != NULL;
  const char *saved_cblas_env_value = getenv("KLS_ENABLE_CBLAS_SUPERNODE");
  char *saved_cblas_env =
    saved_cblas_env_value != NULL ? strdup(saved_cblas_env_value) : NULL;
  const int had_saved_cblas_env = saved_cblas_env_value != NULL;
  const char *saved_native_env_value =
    getenv("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR");
  char *saved_native_env =
    saved_native_env_value != NULL ? strdup(saved_native_env_value) : NULL;
  const int had_saved_native_env = saved_native_env_value != NULL;

  int ok = 1;
  if (had_saved_env && saved_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_ROW_REFACTOR\n");
    ok = 0;
  }
  if (had_saved_trsv_env && saved_trsv_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_COMPACT_SUPERNODE_TRSV\n");
    ok = 0;
  }
  if (had_saved_cblas_env && saved_cblas_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_CBLAS_SUPERNODE\n");
    ok = 0;
  }
  if (had_saved_native_env && saved_native_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR\n");
    ok = 0;
  }
  if (!require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, n, ap, ai, 0,
                                        &options),
                        "analyze unchecked dense compact panel")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor unchecked dense compact panel")) ok = 0;
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_COMPACT_SUPERNODE_TRSV", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_COMPACT_SUPERNODE_TRSV=1");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CBLAS_SUPERNODE", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CBLAS_SUPERNODE=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR=1");
    ok = 0;
  }
  if (ok && !require_ok(kls_refactor(solver, ax1),
                        "refactor unchecked dense compact panel")) ok = 0;
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
  if (had_saved_trsv_env && saved_trsv_env != NULL) {
    if (setenv("KLS_ENABLE_COMPACT_SUPERNODE_TRSV", saved_trsv_env, 1) != 0) {
      perror("restore KLS_ENABLE_COMPACT_SUPERNODE_TRSV");
      ok = 0;
    }
  } else if (!had_saved_trsv_env) {
    if (unsetenv("KLS_ENABLE_COMPACT_SUPERNODE_TRSV") != 0) {
      perror("unsetenv KLS_ENABLE_COMPACT_SUPERNODE_TRSV");
      ok = 0;
    }
  }
  if (had_saved_cblas_env && saved_cblas_env != NULL) {
    if (setenv("KLS_ENABLE_CBLAS_SUPERNODE", saved_cblas_env, 1) != 0) {
      perror("restore KLS_ENABLE_CBLAS_SUPERNODE");
      ok = 0;
    }
  } else if (!had_saved_cblas_env) {
    if (unsetenv("KLS_ENABLE_CBLAS_SUPERNODE") != 0) {
      perror("unsetenv KLS_ENABLE_CBLAS_SUPERNODE");
      ok = 0;
    }
  }
  if (had_saved_native_env && saved_native_env != NULL) {
    if (setenv("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR",
               saved_native_env, 1) != 0) {
      perror("restore KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR");
      ok = 0;
    }
  } else if (!had_saved_native_env) {
    if (unsetenv("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR") != 0) {
      perror("unsetenv KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR");
      ok = 0;
    }
  }
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve unchecked dense compact panel")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats unchecked dense compact panel")) {
    ok = 0;
  }
  if (ok && (stats.row_refactor_last_run != 1 ||
             stats.row_refactor_last_checked != 0 ||
             stats.row_refactor_last_parallel != 0 ||
             stats.row_refactor_dense_segment_count < 1 ||
             stats.row_refactor_last_compact_dense_panel != 1 ||
             stats.row_refactor_compact_dense_panel_count < 1 ||
             stats.row_refactor_compact_dense_panel_eligible_count < 1 ||
             stats.row_refactor_compact_dense_panel_eligible_rows < 1 ||
             stats.row_refactor_compact_dense_panel_update_work <= 0.0 ||
             stats.row_refactor_compact_dense_panel_entries <= 0.0 ||
             stats.row_refactor_compact_dense_panel_persistent_groups < 1 ||
             stats.row_refactor_compact_dense_panel_persistent_entries < 1 ||
             stats.row_refactor_last_compact_dense_panel_persistent != 1 ||
             stats.row_refactor_compact_dense_panel_persistent_run_count < 1 ||
             stats.row_refactor_last_compact_dense_panel_blocked != 1 ||
             stats.row_refactor_compact_dense_panel_blocked_run_count < 1 ||
             stats.row_refactor_compact_dense_panel_blocked_rows < lead ||
             stats.row_refactor_compact_dense_panel_blocked_entries <= 0 ||
             stats.row_refactor_native_row_panel_enabled != 1 ||
             stats.row_refactor_last_native_row_panel != 1 ||
             stats.row_refactor_native_row_panel_count < 1 ||
             stats.row_refactor_native_row_panel_rows < lead ||
             stats.row_refactor_native_row_panel_entries <= 0 ||
             stats.row_refactor_native_row_panel_blocked_count < 1 ||
             stats.row_refactor_native_row_panel_blocked_rows < lead ||
             stats.row_refactor_native_row_panel_blocked_entries <= 0 ||
             stats.row_refactor_native_row_panel_fallback_count != 0 ||
             stats.row_refactor_native_row_panel_checked_reject_count != 0 ||
             stats.row_refactor_last_compact_dense_panel_direct_input_rows < 1 ||
             stats.row_refactor_compact_dense_panel_direct_input_rows < 1 ||
             stats.row_refactor_segment_input_target_rows < 1 ||
             stats.row_refactor_segment_input_target_entries < 1 ||
             stats.row_refactor_last_segment_target_input_rows <
               stats.row_refactor_last_compact_dense_panel_direct_input_rows ||
             stats.row_refactor_segment_target_input_rows <
               stats.row_refactor_compact_dense_panel_direct_input_rows ||
             stats.row_refactor_last_compact_panel_solve_values < 1 ||
             stats.row_refactor_compact_panel_solve_values < 1 ||
             stats.row_refactor_last_compact_panel_group_solve_rows < lead ||
             stats.row_refactor_compact_panel_group_solve_rows < lead ||
             stats.row_refactor_last_compact_panel_group_solve_entries <= 0 ||
             stats.row_refactor_compact_panel_group_solve_entries <= 0 ||
             stats.row_refactor_last_compact_supernode_update != 1 ||
             stats.row_refactor_compact_supernode_update_count < 1 ||
             stats.row_refactor_compact_supernode_update_rows < lead ||
             stats.row_refactor_compact_supernode_update_entries <= 0 ||
             stats.row_refactor_last_compact_supernode_gemv != 1 ||
             stats.row_refactor_compact_supernode_gemv_count < 1 ||
             stats.row_refactor_compact_supernode_gemv_rows < lead ||
             stats.row_refactor_compact_supernode_gemv_entries <= 0 ||
             stats.row_refactor_last_compact_supernode_trsv != 1 ||
             stats.row_refactor_compact_supernode_trsv_count < 1 ||
             stats.row_refactor_compact_supernode_trsv_rows < lead ||
             stats.row_refactor_compact_supernode_trsv_entries <= 0)) {
    fprintf(stderr,
            "unexpected unchecked dense compact-panel stats: row=%d/%d/%d"
            ", dense_segments=%" PRId64 ", compact=%d/%" PRId64
            ", eligible=%" PRId64 "/%" PRId64
            ", work=%.17g, entries=%.17g"
            ", persistent=%" PRId64 "/%" PRId64
            ", persistent_used=%d/%" PRId64
            ", blocked=%d/%" PRId64 "/%" PRId64 "/%" PRId64
            ", native=%d/%d/%" PRId64 "/%" PRId64 "/%" PRId64
            "/%" PRId64 "/%" PRId64 "/%" PRId64 "/%" PRId64 "/%" PRId64
            ", direct_input=%" PRId64 "/%" PRId64
            ", target_map=%" PRId64 "/%" PRId64
            ", target_input=%" PRId64 "/%" PRId64
            ", panel_solve=%" PRId64 "/%" PRId64
            ", group_solve=%" PRId64 "/%" PRId64 "/%" PRId64 "/%" PRId64
            ", supernode=%d/%" PRId64 "/%" PRId64 "/%" PRId64
            ", gemv=%d/%" PRId64 "/%" PRId64 "/%" PRId64
            ", trsv=%d/%" PRId64 "/%" PRId64 "/%" PRId64 "\n",
            stats.row_refactor_last_run,
            stats.row_refactor_last_checked,
            stats.row_refactor_last_parallel,
            stats.row_refactor_dense_segment_count,
            stats.row_refactor_last_compact_dense_panel,
            stats.row_refactor_compact_dense_panel_count,
            stats.row_refactor_compact_dense_panel_eligible_count,
            stats.row_refactor_compact_dense_panel_eligible_rows,
            stats.row_refactor_compact_dense_panel_update_work,
            stats.row_refactor_compact_dense_panel_entries,
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
            stats.row_refactor_native_row_panel_blocked_count,
            stats.row_refactor_native_row_panel_blocked_rows,
            stats.row_refactor_native_row_panel_blocked_entries,
            stats.row_refactor_native_row_panel_fallback_count,
            stats.row_refactor_native_row_panel_checked_reject_count,
            stats.row_refactor_last_compact_dense_panel_direct_input_rows,
            stats.row_refactor_compact_dense_panel_direct_input_rows,
            stats.row_refactor_segment_input_target_rows,
            stats.row_refactor_segment_input_target_entries,
            stats.row_refactor_last_segment_target_input_rows,
            stats.row_refactor_segment_target_input_rows,
            stats.row_refactor_last_compact_panel_solve_values,
            stats.row_refactor_compact_panel_solve_values,
            stats.row_refactor_last_compact_panel_group_solve_rows,
            stats.row_refactor_compact_panel_group_solve_rows,
            stats.row_refactor_last_compact_panel_group_solve_entries,
            stats.row_refactor_compact_panel_group_solve_entries,
            stats.row_refactor_last_compact_supernode_update,
            stats.row_refactor_compact_supernode_update_count,
            stats.row_refactor_compact_supernode_update_rows,
            stats.row_refactor_compact_supernode_update_entries,
            stats.row_refactor_last_compact_supernode_gemv,
            stats.row_refactor_compact_supernode_gemv_count,
            stats.row_refactor_compact_supernode_gemv_rows,
            stats.row_refactor_compact_supernode_gemv_entries,
            stats.row_refactor_last_compact_supernode_trsv,
            stats.row_refactor_compact_supernode_trsv_count,
            stats.row_refactor_compact_supernode_trsv_rows,
            stats.row_refactor_compact_supernode_trsv_entries);
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
  double *residual = (double *)malloc((size_t)n * sizeof(*residual));
  if (residual == NULL) {
    ok = 0;
  }
  for (int32_t row = 0; row < n; ++row) {
    if (residual != NULL) {
      residual[row] = -b[row];
    }
    if (fabs(b[row]) > max_rhs) {
      max_rhs = fabs(b[row]);
    }
  }
  if (residual != NULL) {
    for (int32_t col = 0; col < n; ++col) {
      for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
        residual[ai[p]] += ax1[p] * x[col];
      }
    }
    for (int32_t row = 0; row < n; ++row) {
      if (fabs(residual[row]) > max_residual) {
        max_residual = fabs(residual[row]);
      }
    }
  }
  for (int32_t rhs = 0; rhs < multi_rhs; ++rhs) {
    double *rhs_expected = expected_multi + (size_t)rhs * (size_t)n;
    double *rhs_b = b_multi + (size_t)rhs * (size_t)n;
    for (int32_t col = 0; col < n; ++col) {
      rhs_expected[col] =
        expected[col] * (1.0 + 0.03125 * (double)rhs) +
        0.0125 * (double)((rhs + col) % 5);
    }
    for (int32_t col = 0; col < n; ++col) {
      for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
        rhs_b[ai[p]] += ax1[p] * rhs_expected[col];
        b_transpose_multi[(size_t)rhs * (size_t)n + (size_t)col] +=
          ax1[p] * rhs_expected[ai[p]];
      }
    }
  }
  const double rel_resid = max_residual / (1.0 + max_rhs);
  if (ok && (max_solution_error > 1e-8 || rel_resid > 1e-9)) {
    fprintf(stderr,
            "unexpected unchecked dense compact-panel accuracy: "
            "max_x_err=%.17g, rel_resid=%.17g\n",
            max_solution_error, rel_resid);
    ok = 0;
  }

  if (ok && !require_ok(kls_solve(solver, multi_rhs, b_multi, 0, x_multi, 0),
                        "multi-RHS solve unchecked dense compact panel")) {
    ok = 0;
  }
  kls_stats multi_stats;
  multi_stats.struct_size = sizeof(multi_stats);
  if (ok && !require_ok(kls_get_stats(solver, &multi_stats),
                        "multi-RHS stats unchecked dense compact panel")) {
    ok = 0;
  }
  if (ok &&
      (multi_stats.row_refactor_last_compact_panel_group_solve_rows < lead ||
       multi_stats.row_refactor_last_compact_panel_group_solve_entries <= 0 ||
       multi_stats.row_refactor_compact_panel_group_solve_rows <
         stats.row_refactor_compact_panel_group_solve_rows + lead)) {
    fprintf(stderr,
            "unexpected multi-RHS compact-panel solve stats: rows=%" PRId64
            "/%" PRId64 ", entries=%" PRId64 "/%" PRId64 "\n",
            multi_stats.row_refactor_last_compact_panel_group_solve_rows,
            multi_stats.row_refactor_compact_panel_group_solve_rows,
            multi_stats.row_refactor_last_compact_panel_group_solve_entries,
            multi_stats.row_refactor_compact_panel_group_solve_entries);
    ok = 0;
  }
  double max_multi_solution_error = 0.0;
  double max_multi_residual = 0.0;
  double max_multi_rhs = 0.0;
  for (int32_t rhs = 0; rhs < multi_rhs; ++rhs) {
    const double *rhs_expected =
      expected_multi + (size_t)rhs * (size_t)n;
    const double *rhs_b = b_multi + (size_t)rhs * (size_t)n;
    const double *rhs_x = x_multi + (size_t)rhs * (size_t)n;
    for (int32_t row = 0; row < n; ++row) {
      const double err = fabs(rhs_x[row] - rhs_expected[row]);
      if (err > max_multi_solution_error) {
        max_multi_solution_error = err;
      }
      double residual_value = -rhs_b[row];
      if (fabs(rhs_b[row]) > max_multi_rhs) {
        max_multi_rhs = fabs(rhs_b[row]);
      }
      for (int32_t col = 0; col < n; ++col) {
        for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
          if (ai[p] == row) {
            residual_value += ax1[p] * rhs_x[col];
          }
        }
      }
      if (fabs(residual_value) > max_multi_residual) {
        max_multi_residual = fabs(residual_value);
      }
    }
  }
  const double multi_rel_resid =
    max_multi_residual / (1.0 + max_multi_rhs);
  if (ok &&
      (max_multi_solution_error > 1e-8 || multi_rel_resid > 1e-9)) {
    fprintf(stderr,
            "unexpected multi-RHS compact-panel accuracy: "
            "max_x_err=%.17g, rel_resid=%.17g\n",
            max_multi_solution_error, multi_rel_resid);
    ok = 0;
  }

  if (ok && !require_ok(kls_solve_transpose(solver, 1, b_transpose, 0,
                                            x_transpose, 0),
                        "transpose solve unchecked dense compact panel")) {
    ok = 0;
  }
  kls_stats transpose_stats;
  transpose_stats.struct_size = sizeof(transpose_stats);
  if (ok && !require_ok(kls_get_stats(solver, &transpose_stats),
                        "transpose stats unchecked dense compact panel")) {
    ok = 0;
  }
  if (ok &&
      (transpose_stats.row_refactor_last_compact_panel_group_solve_rows <
         lead ||
       transpose_stats.row_refactor_last_compact_panel_group_solve_entries <=
         0 ||
       transpose_stats.row_refactor_compact_panel_group_solve_rows <
         multi_stats.row_refactor_compact_panel_group_solve_rows + lead)) {
    fprintf(stderr,
            "unexpected transpose compact-panel solve stats: rows=%" PRId64
            "/%" PRId64 ", entries=%" PRId64 "/%" PRId64 "\n",
            transpose_stats.row_refactor_last_compact_panel_group_solve_rows,
            transpose_stats.row_refactor_compact_panel_group_solve_rows,
            transpose_stats.row_refactor_last_compact_panel_group_solve_entries,
            transpose_stats.row_refactor_compact_panel_group_solve_entries);
    ok = 0;
  }
  double max_transpose_solution_error = 0.0;
  for (int32_t i = 0; i < n; ++i) {
    const double err = fabs(x_transpose[i] - expected[i]);
    if (err > max_transpose_solution_error) {
      max_transpose_solution_error = err;
    }
  }
  if (ok && max_transpose_solution_error > 1e-8) {
    fprintf(stderr,
            "unexpected transpose compact-panel accuracy: max_x_err=%.17g\n",
            max_transpose_solution_error);
    ok = 0;
  }

  if (ok && !require_ok(kls_solve_transpose(solver, multi_rhs,
                                            b_transpose_multi, 0,
                                            x_transpose_multi, 0),
                        "transpose multi-RHS solve unchecked dense compact panel")) {
    ok = 0;
  }
  kls_stats transpose_multi_stats;
  transpose_multi_stats.struct_size = sizeof(transpose_multi_stats);
  if (ok && !require_ok(kls_get_stats(solver, &transpose_multi_stats),
                        "transpose multi-RHS stats unchecked dense compact panel")) {
    ok = 0;
  }
  if (ok &&
      (transpose_multi_stats.row_refactor_last_compact_panel_group_solve_rows <
         lead ||
       transpose_multi_stats.row_refactor_last_compact_panel_group_solve_entries <=
         0 ||
       transpose_multi_stats.row_refactor_compact_panel_group_solve_rows <
         transpose_stats.row_refactor_compact_panel_group_solve_rows + lead)) {
    fprintf(stderr,
            "unexpected transpose multi-RHS compact-panel solve stats:"
            " rows=%" PRId64 "/%" PRId64 ", entries=%" PRId64 "/%" PRId64
            "\n",
            transpose_multi_stats.row_refactor_last_compact_panel_group_solve_rows,
            transpose_multi_stats.row_refactor_compact_panel_group_solve_rows,
            transpose_multi_stats.row_refactor_last_compact_panel_group_solve_entries,
            transpose_multi_stats.row_refactor_compact_panel_group_solve_entries);
    ok = 0;
  }
  double max_transpose_multi_solution_error = 0.0;
  for (int32_t rhs = 0; rhs < multi_rhs; ++rhs) {
    const double *rhs_expected =
      expected_multi + (size_t)rhs * (size_t)n;
    const double *rhs_x =
      x_transpose_multi + (size_t)rhs * (size_t)n;
    for (int32_t row = 0; row < n; ++row) {
      const double err = fabs(rhs_x[row] - rhs_expected[row]);
      if (err > max_transpose_multi_solution_error) {
        max_transpose_multi_solution_error = err;
      }
    }
  }
  if (ok && max_transpose_multi_solution_error > 1e-8) {
    fprintf(stderr,
            "unexpected transpose multi-RHS compact-panel accuracy:"
            " max_x_err=%.17g\n",
            max_transpose_multi_solution_error);
    ok = 0;
  }

  if (ok && !require_ok(kls_create(&disabled_solver),
                        "create native row-panel disabled solver")) {
    ok = 0;
  }
  if (ok && !require_ok(kls_analyze_csc(disabled_solver, KLS_INDEX_INT32,
                                        n, ap, ai, 0, &options),
                        "analyze native row-panel disabled")) {
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(disabled_solver, ax0),
                        "factor native row-panel disabled")) {
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CBLAS_SUPERNODE", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CBLAS_SUPERNODE=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR=0");
    ok = 0;
  }
  if (ok && !require_ok(kls_refactor(disabled_solver, ax1),
                        "refactor native row-panel disabled")) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_ROW_REFACTOR",
                         had_saved_env, saved_env)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_CBLAS_SUPERNODE",
                         had_saved_cblas_env, saved_cblas_env)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR",
                         had_saved_native_env, saved_native_env)) {
    ok = 0;
  }
  kls_stats disabled_stats;
  disabled_stats.struct_size = sizeof(disabled_stats);
  if (ok && !require_ok(kls_get_stats(disabled_solver, &disabled_stats),
                        "stats native row-panel disabled")) {
    ok = 0;
  }
  if (ok && (disabled_stats.row_refactor_last_run != 1 ||
             disabled_stats.row_refactor_native_row_panel_enabled != 0 ||
             disabled_stats.row_refactor_last_native_row_panel != 0 ||
             disabled_stats.row_refactor_native_row_panel_count != 0 ||
             disabled_stats.row_refactor_native_row_panel_blocked_count != 0 ||
             disabled_stats.row_refactor_native_row_panel_fallback_count != 0)) {
    fprintf(stderr,
            "unexpected native row-panel disabled stats: row=%d,"
            " native=%d/%d/%" PRId64 "/%" PRId64 "/%" PRId64 "\n",
            disabled_stats.row_refactor_last_run,
            disabled_stats.row_refactor_native_row_panel_enabled,
            disabled_stats.row_refactor_last_native_row_panel,
            disabled_stats.row_refactor_native_row_panel_count,
            disabled_stats.row_refactor_native_row_panel_blocked_count,
            disabled_stats.row_refactor_native_row_panel_fallback_count);
    ok = 0;
  }

  kls_destroy(solver);
  kls_destroy(disabled_solver);
  free(ap);
  free(ai);
  free(ax0);
  free(ax1);
  free(b);
  free(x);
  free(b_transpose);
  free(x_transpose);
  free(expected);
  free(b_multi);
  free(x_multi);
  free(b_transpose_multi);
  free(x_transpose_multi);
  free(expected_multi);
  free(residual);
  free(saved_env);
  free(saved_trsv_env);
  free(saved_cblas_env);
  free(saved_native_env);
  return ok;
}

static int test_row_dense_compact_panel_scalar_update(void) {
  const int32_t lead = 48;
  const int32_t n = 65;
  const int32_t scalar_row = lead;
  const int32_t trail = n - lead - 1;
  const size_t nnz =
    (size_t)lead * (size_t)lead + 2u +
    (size_t)trail * ((size_t)lead + 1u);
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

  for (int32_t col = 0; col < n; ++col) {
    expected[col] = 0.5 + 0.015625 * (double)((7 * col) % 19);
  }

  size_t pos = 0;
  for (int32_t col = 0; col < n; ++col) {
    ap[col] = (int32_t)pos;
    if (col < lead) {
      for (int32_t row = 0; row < lead; ++row) {
        const size_t p = pos++;
        ai[p] = row;
        ax0[p] = row == col
          ? 21.0 + 0.01 * (double)col
          : 0.0004 * (1.0 + (double)((row + 11 * col) % 17));
        ax1[p] = ax0[p] + (row == col ? 0.025 : 1.0e-6);
      }
      if (col == lead - 1) {
        const size_t p = pos++;
        ai[p] = scalar_row;
        ax0[p] = 0.002;
        ax1[p] = 0.0021;
      }
    } else if (col == scalar_row) {
      const size_t p = pos++;
      ai[p] = scalar_row;
      ax0[p] = 19.0;
      ax1[p] = 19.05;
    } else {
      for (int32_t row = 0; row < lead; ++row) {
        const size_t p = pos++;
        ai[p] = row;
        ax0[p] = 0.00025 * (1.0 + (double)((row + 3 * col) % 13));
        ax1[p] = ax0[p] + 1.0e-6 * (double)(((row + col) % 5) - 2);
      }
      const size_t p = pos++;
      ai[p] = col;
      ax0[p] = 17.0 + 0.02 * (double)col;
      ax1[p] = ax0[p] + 0.03;
    }
  }
  ap[n] = (int32_t)pos;
  if (pos != nnz) {
    fprintf(stderr, "unexpected scalar compact-panel fixture nnz: %zu/%zu\n",
            pos, nnz);
    free(ap);
    free(ai);
    free(ax0);
    free(ax1);
    free(b);
    free(x);
    free(expected);
    return 0;
  }
  for (int32_t col = 0; col < n; ++col) {
    for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
      b[ai[p]] += ax1[p] * expected[col];
    }
  }

  const char *saved_row_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_row = saved_row_value != NULL ? strdup(saved_row_value) : NULL;
  const int had_row = saved_row_value != NULL;
  const char *saved_checked_value = getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  char *saved_checked =
    saved_checked_value != NULL ? strdup(saved_checked_value) : NULL;
  const int had_checked = saved_checked_value != NULL;
  const char *saved_native_value =
    getenv("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR");
  char *saved_native =
    saved_native_value != NULL ? strdup(saved_native_value) : NULL;
  const int had_native = saved_native_value != NULL;

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

  int ok = 1;
  if ((had_row && saved_row == NULL) ||
      (had_checked && saved_checked == NULL) ||
      (had_native && saved_native == NULL)) {
    fprintf(stderr, "failed to save scalar compact-panel env\n");
    ok = 0;
  }
  if (!require_ok(kls_create(&solver), "create scalar compact panel")) {
    ok = 0;
  }
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, n, ap, ai, 0,
                                        &options),
                        "analyze scalar compact panel")) {
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor scalar compact panel base")) {
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR", "auto", 1) != 0) {
    perror("setenv KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR=auto");
    ok = 0;
  }
  if (ok && !require_ok(kls_refactor(solver, ax1),
                        "refactor scalar compact panel")) {
    ok = 0;
  }
  if (had_row && saved_row != NULL) {
    if (setenv("KLS_ENABLE_ROW_REFACTOR", saved_row, 1) != 0) {
      perror("restore KLS_ENABLE_ROW_REFACTOR");
      ok = 0;
    }
  } else if (!had_row && unsetenv("KLS_ENABLE_ROW_REFACTOR") != 0) {
    perror("unsetenv KLS_ENABLE_ROW_REFACTOR");
    ok = 0;
  }
  if (had_checked && saved_checked != NULL) {
    if (setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", saved_checked, 1) != 0) {
      perror("restore KLS_ENABLE_CHECKED_ROW_REFACTOR");
      ok = 0;
    }
  } else if (!had_checked && unsetenv("KLS_ENABLE_CHECKED_ROW_REFACTOR") != 0) {
    perror("unsetenv KLS_ENABLE_CHECKED_ROW_REFACTOR");
    ok = 0;
  }
  if (had_native && saved_native != NULL) {
    if (setenv("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR", saved_native, 1) != 0) {
      perror("restore KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR");
      ok = 0;
    }
  } else if (!had_native &&
             unsetenv("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR") != 0) {
    perror("unsetenv KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR");
    ok = 0;
  }
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve scalar compact panel")) {
    ok = 0;
  }

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats scalar compact panel")) {
    ok = 0;
  }
  if (ok &&
      (stats.row_refactor_last_compact_dense_panel_persistent != 1 ||
       stats.row_refactor_last_compact_panel_scalar_update_rows < 1 ||
       stats.row_refactor_compact_panel_scalar_update_rows < 1 ||
       stats.row_refactor_last_compact_panel_scalar_update_entries < trail ||
       stats.row_refactor_compact_panel_scalar_update_entries < trail)) {
    fprintf(stderr,
            "unexpected scalar compact-panel stats: persistent=%d"
            ", scalar_rows=%" PRId64 "/%" PRId64
            ", scalar_entries=%" PRId64 "/%" PRId64 "\n",
            stats.row_refactor_last_compact_dense_panel_persistent,
            stats.row_refactor_last_compact_panel_scalar_update_rows,
            stats.row_refactor_compact_panel_scalar_update_rows,
            stats.row_refactor_last_compact_panel_scalar_update_entries,
            stats.row_refactor_compact_panel_scalar_update_entries);
    ok = 0;
  }

  double max_solution_error = 0.0;
  double max_residual = 0.0;
  double max_rhs = 0.0;
  for (int32_t row = 0; row < n; ++row) {
    const double err = fabs(x[row] - expected[row]);
    if (err > max_solution_error) {
      max_solution_error = err;
    }
    double residual = -b[row];
    if (fabs(b[row]) > max_rhs) {
      max_rhs = fabs(b[row]);
    }
    for (int32_t col = 0; col < n; ++col) {
      for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
        if (ai[p] == row) {
          residual += ax1[p] * x[col];
        }
      }
    }
    if (fabs(residual) > max_residual) {
      max_residual = fabs(residual);
    }
  }
  const double rel_resid = max_residual / (max_rhs > 0.0 ? max_rhs : 1.0);
  if (ok && (max_solution_error > 1.0e-8 || rel_resid > 1.0e-10)) {
    fprintf(stderr,
            "unexpected scalar compact-panel accuracy: max_x_err=%.17g"
            ", rel_resid=%.17g\n",
            max_solution_error, rel_resid);
    ok = 0;
  }

  kls_destroy(solver);
  free(saved_row);
  free(saved_checked);
  free(saved_native);
  free(ap);
  free(ai);
  free(ax0);
  free(ax1);
  free(b);
  free(x);
  free(expected);
  return ok;
}

static int test_unchecked_row_dense_native_direct_input(void) {
  const int32_t n = 46;
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
    expected[col] = 0.5 + 0.0625 * (double)((7 * col) % 19);
    for (int32_t row = 0; row < n; ++row) {
      const size_t p = (size_t)col * (size_t)n + (size_t)row;
      ai[p] = row;
      ax0[p] = row == col
        ? 28.0 + 0.015 * (double)col
        : 0.00035 * (1.0 + (double)((row + 5 * col) % 17));
      ax1[p] = ax0[p] + (row == col
        ? 0.03 * (double)((col % 5) + 1)
        : 1.0e-5 * (double)(((row + 3 * col) % 7) - 3));
    }
  }
  for (int32_t col = 0; col < n; ++col) {
    for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
      b[ai[p]] += ax1[p] * expected[col];
    }
  }

  const char *saved_env_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_env = saved_env_value != NULL ? strdup(saved_env_value) : NULL;
  const int had_saved_env = saved_env_value != NULL;

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

  int ok = 1;
  if (had_saved_env && saved_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_ROW_REFACTOR\n");
    ok = 0;
  }
  if (!require_ok(kls_create(&solver), "create dense native direct input")) {
    ok = 0;
  }
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, n, ap, ai, 0,
                                        &options),
                        "analyze dense native direct input")) {
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor dense native direct input")) {
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && !require_ok(kls_refactor(solver, ax1),
                        "refactor dense native direct input")) {
    ok = 0;
  }
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

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats dense native direct input")) {
    ok = 0;
  }
  if (ok && (stats.row_refactor_last_run != 1 ||
             stats.row_refactor_last_checked != 0 ||
             stats.row_refactor_last_parallel != 0 ||
             stats.row_refactor_dense_segment_count < 1 ||
             stats.row_refactor_last_compact_dense_panel != 0 ||
             stats.row_refactor_compact_dense_panel_count != 0 ||
             stats.row_refactor_last_dense_segment_direct_input_rows != n ||
             stats.row_refactor_dense_segment_direct_input_rows < n ||
             stats.row_refactor_segment_input_target_rows != n ||
             stats.row_refactor_segment_input_target_entries <= 0 ||
             stats.row_refactor_last_segment_target_input_rows != n ||
             stats.row_refactor_segment_target_input_rows < n ||
             stats.row_refactor_values_dirty != 1 ||
             stats.row_refactor_last_lazy_value_scatter != 1)) {
    fprintf(stderr,
            "unexpected dense native direct-input stats: row=%d/%d/%d"
            ", dense_segments=%" PRId64 ", compact=%d/%" PRId64
            ", direct=%" PRId64 "/%" PRId64
            ", target=%" PRId64 "/%" PRId64 "/%" PRId64 "/%" PRId64
            ", dirty/lazy=%d/%d\n",
            stats.row_refactor_last_run,
            stats.row_refactor_last_checked,
            stats.row_refactor_last_parallel,
            stats.row_refactor_dense_segment_count,
            stats.row_refactor_last_compact_dense_panel,
            stats.row_refactor_compact_dense_panel_count,
            stats.row_refactor_last_dense_segment_direct_input_rows,
            stats.row_refactor_dense_segment_direct_input_rows,
            stats.row_refactor_segment_input_target_rows,
            stats.row_refactor_segment_input_target_entries,
            stats.row_refactor_last_segment_target_input_rows,
            stats.row_refactor_segment_target_input_rows,
            stats.row_refactor_values_dirty,
            stats.row_refactor_last_lazy_value_scatter);
    ok = 0;
  }
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve dense native direct input")) {
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
    if (fabs(b[row]) > max_rhs) {
      max_rhs = fabs(b[row]);
    }
    for (int32_t col = 0; col < n; ++col) {
      residual += ax1[(size_t)col * (size_t)n + (size_t)row] * x[col];
    }
    if (fabs(residual) > max_residual) {
      max_residual = fabs(residual);
    }
  }
  const double relative_residual =
    max_residual / (max_rhs > 0.0 ? max_rhs : 1.0);
  if (ok && (max_solution_error > 1.0e-9 ||
             relative_residual > 1.0e-11)) {
    fprintf(stderr,
            "unexpected dense native direct-input accuracy:"
            " max_x_err=%.17g, rel_resid=%.17g\n",
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

static int test_unchecked_row_generic_target_direct_input(void) {
  const int32_t n = 12;
  const int32_t nnz = 3 * n - 2;
  int32_t *ap = (int32_t *)malloc(((size_t)n + 1u) * sizeof(*ap));
  int32_t *ai = (int32_t *)malloc((size_t)nnz * sizeof(*ai));
  double *ax0 = (double *)malloc((size_t)nnz * sizeof(*ax0));
  double *ax1 = (double *)malloc((size_t)nnz * sizeof(*ax1));
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

  int32_t p = 0;
  for (int32_t col = 0; col < n; ++col) {
    ap[col] = p;
    if (col > 0) {
      ai[p] = col - 1;
      ax0[p] = -0.025 * (1.0 + (double)(col % 3));
      ax1[p] = ax0[p] - 0.001 * (double)((col % 2) + 1);
      p++;
    }
    ai[p] = col;
    ax0[p] = 4.0 + 0.03 * (double)col;
    ax1[p] = ax0[p] + 0.01 * (double)((col % 4) + 1);
    p++;
    if (col + 1 < n) {
      ai[p] = col + 1;
      ax0[p] = 0.02 * (1.0 + (double)(col % 5));
      ax1[p] = ax0[p] + 0.0015 * (double)((col % 3) + 1);
      p++;
    }
  }
  ap[n] = p;
  if (p != nnz) {
    fprintf(stderr, "unexpected generic target direct-input nnz=%d\n", p);
    free(ap);
    free(ai);
    free(ax0);
    free(ax1);
    free(b);
    free(x);
    free(expected);
    return 0;
  }
  for (int32_t col = 0; col < n; ++col) {
    expected[col] = 0.25 + 0.05 * (double)((5 * col) % 11);
  }
  for (int32_t col = 0; col < n; ++col) {
    for (int32_t q = ap[col]; q < ap[col + 1]; ++q) {
      b[ai[q]] += ax1[q] * expected[col];
    }
  }

  const char *saved_env_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_env = saved_env_value != NULL ? strdup(saved_env_value) : NULL;
  const int had_saved_env = saved_env_value != NULL;

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

  int ok = 1;
  if (had_saved_env && saved_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_ROW_REFACTOR\n");
    ok = 0;
  }
  if (!require_ok(kls_create(&solver),
                  "create generic target direct input")) {
    ok = 0;
  }
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, n, ap, ai,
                                        0, &options),
                        "analyze generic target direct input")) {
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor generic target direct input")) {
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && !require_ok(kls_refactor(solver, ax1),
                        "refactor generic target direct input")) {
    ok = 0;
  }
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

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats generic target direct input")) {
    ok = 0;
  }
  if (ok && (stats.row_refactor_last_run != 1 ||
             stats.row_refactor_last_checked != 0 ||
             stats.row_refactor_last_parallel != 0 ||
             stats.row_refactor_segment_input_target_rows != n ||
             stats.row_refactor_segment_input_target_entries != nnz ||
             stats.row_refactor_last_segment_target_input_rows != n ||
             stats.row_refactor_segment_target_input_rows < n ||
             stats.row_refactor_last_dense_segment_direct_input_rows != 0 ||
             stats.row_refactor_last_sparse_segment_direct_input_rows >=
               stats.row_refactor_last_segment_target_input_rows ||
             stats.row_refactor_values_dirty != 1)) {
    fprintf(stderr,
            "unexpected generic target direct-input stats: row=%d/%d/%d"
            ", segments=%" PRId64 ", planned=%" PRId64 "/%" PRId64
            ", used=%" PRId64 "/%" PRId64
            ", dense/sparse=%" PRId64 "/%" PRId64 ", dirty=%d\n",
            stats.row_refactor_last_run,
            stats.row_refactor_last_checked,
            stats.row_refactor_last_parallel,
            stats.row_refactor_segment_count,
            stats.row_refactor_segment_input_target_rows,
            stats.row_refactor_segment_input_target_entries,
            stats.row_refactor_last_segment_target_input_rows,
            stats.row_refactor_segment_target_input_rows,
            stats.row_refactor_last_dense_segment_direct_input_rows,
            stats.row_refactor_last_sparse_segment_direct_input_rows,
            stats.row_refactor_values_dirty);
    ok = 0;
  }
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve generic target direct input")) {
    ok = 0;
  }
  double max_solution_error = 0.0;
  double max_residual = 0.0;
  double max_rhs = 0.0;
  for (int32_t row = 0; row < n; ++row) {
    const double err = fabs(x[row] - expected[row]);
    if (err > max_solution_error) {
      max_solution_error = err;
    }
    double residual = -b[row];
    if (fabs(b[row]) > max_rhs) {
      max_rhs = fabs(b[row]);
    }
    for (int32_t col = 0; col < n; ++col) {
      for (int32_t q = ap[col]; q < ap[col + 1]; ++q) {
        if (ai[q] == row) {
          residual += ax1[q] * x[col];
        }
      }
    }
    if (fabs(residual) > max_residual) {
      max_residual = fabs(residual);
    }
  }
  const double rel_resid = max_residual / (max_rhs > 0.0 ? max_rhs : 1.0);
  if (ok && (max_solution_error > 1.0e-8 || rel_resid > 1.0e-10)) {
    fprintf(stderr,
            "unexpected generic target direct-input accuracy: max_x_err=%.17g"
            ", rel_resid=%.17g\n",
            max_solution_error, rel_resid);
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

static int test_unchecked_row_sparse_segment_direct_input(void) {
  const int32_t n = 48;
  const size_t nnz = ((size_t)n * ((size_t)n + 1u)) / 2u + (size_t)n - 1u;
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

  for (int32_t col = 0; col < n; ++col) {
    expected[col] = 0.625 + 0.03125 * (double)((11 * col) % 23);
  }

  size_t pos = 0;
  for (int32_t col = 0; col < n; ++col) {
    ap[col] = (int32_t)pos;
    for (int32_t row = 0; row <= col; ++row) {
      const size_t p = pos++;
      ai[p] = row;
      ax0[p] = row == col
        ? 32.0 + 0.02 * (double)col
        : 0.00025 * (1.0 + (double)((3 * row + 5 * col) % 19));
      ax1[p] = ax0[p] + (row == col
        ? 0.025 * (double)((col % 7) + 1)
        : 1.0e-5 * (double)(((row + 2 * col) % 9) - 4));
    }
    if (col + 1 < n) {
      const size_t p = pos++;
      ai[p] = col + 1;
      ax0[p] = 0.015 * (1.0 + (double)(col % 5));
      ax1[p] = ax0[p] + 1.0e-5 * (double)((col % 3) - 1);
    }
  }
  ap[n] = (int32_t)pos;
  if (pos != nnz) {
    fprintf(stderr, "unexpected sparse segment fixture nnz: %zu/%zu\n",
            pos, nnz);
    free(ap);
    free(ai);
    free(ax0);
    free(ax1);
    free(b);
    free(x);
    free(expected);
    return 0;
  }
  for (int32_t col = 0; col < n; ++col) {
    for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
      b[ai[p]] += ax1[p] * expected[col];
    }
  }

  const char *saved_env_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_env = saved_env_value != NULL ? strdup(saved_env_value) : NULL;
  const int had_saved_env = saved_env_value != NULL;

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

  int ok = 1;
  if (had_saved_env && saved_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_ROW_REFACTOR\n");
    ok = 0;
  }
  if (!require_ok(kls_create(&solver), "create sparse segment direct input")) {
    ok = 0;
  }
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, n, ap, ai, 0,
                                        &options),
                        "analyze sparse segment direct input")) {
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor sparse segment direct input")) {
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && !require_ok(kls_refactor(solver, ax1),
                        "refactor sparse segment direct input")) {
    ok = 0;
  }
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

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats sparse segment direct input")) {
    ok = 0;
  }
  if (ok && (stats.row_refactor_last_run != 1 ||
             stats.row_refactor_last_checked != 0 ||
             stats.row_refactor_last_parallel != 0 ||
             stats.row_refactor_segment_count < 1 ||
             stats.row_refactor_dense_segment_count != 0 ||
             stats.row_refactor_last_sparse_segment_direct_input_rows != n ||
             stats.row_refactor_sparse_segment_direct_input_rows < n ||
             stats.row_refactor_segment_input_target_rows != n ||
             stats.row_refactor_segment_input_target_entries <= 0 ||
             stats.row_refactor_last_segment_target_input_rows != n ||
             stats.row_refactor_segment_target_input_rows < n ||
             stats.row_refactor_values_dirty != 1 ||
             stats.row_refactor_last_lazy_value_scatter != 1)) {
    fprintf(stderr,
            "unexpected sparse segment direct-input stats: row=%d/%d/%d"
            ", segments=%" PRId64 ", dense_segments=%" PRId64
            ", direct=%" PRId64 "/%" PRId64
            ", target=%" PRId64 "/%" PRId64 "/%" PRId64 "/%" PRId64
            ", dirty/lazy=%d/%d\n",
            stats.row_refactor_last_run,
            stats.row_refactor_last_checked,
            stats.row_refactor_last_parallel,
            stats.row_refactor_segment_count,
            stats.row_refactor_dense_segment_count,
            stats.row_refactor_last_sparse_segment_direct_input_rows,
            stats.row_refactor_sparse_segment_direct_input_rows,
            stats.row_refactor_segment_input_target_rows,
            stats.row_refactor_segment_input_target_entries,
            stats.row_refactor_last_segment_target_input_rows,
            stats.row_refactor_segment_target_input_rows,
            stats.row_refactor_values_dirty,
            stats.row_refactor_last_lazy_value_scatter);
    ok = 0;
  }
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve sparse segment direct input")) {
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
    if (fabs(b[row]) > max_rhs) {
      max_rhs = fabs(b[row]);
    }
    for (int32_t col = 0; col < n; ++col) {
      for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
        if (ai[p] == row) {
          residual += ax1[p] * x[col];
        }
      }
    }
    if (fabs(residual) > max_residual) {
      max_residual = fabs(residual);
    }
  }
  const double relative_residual =
    max_residual / (max_rhs > 0.0 ? max_rhs : 1.0);
  if (ok && (max_solution_error > 1.0e-9 ||
             relative_residual > 1.0e-11)) {
    fprintf(stderr,
            "unexpected sparse segment direct-input accuracy:"
            " max_x_err=%.17g, rel_resid=%.17g\n",
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

static int test_batched_compact_supernode_update_probe(void) {
  const int32_t lead0 = 40;
  const int32_t lead1 = 40;
  const int32_t mid = 40;
  const int32_t trail = 20;
  const int32_t core = lead0 + lead1 + mid;
  const int32_t n = core + trail;
  const int32_t trail_break_col0 = n - 2;
  const size_t nnz =
    (size_t)core * (size_t)core +
    (size_t)(trail - 2) * ((size_t)core + 1u) +
    ((size_t)core - (size_t)lead0 + 1u) +
    ((size_t)core - (size_t)lead0 - (size_t)lead1 + 1u);
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

  for (int32_t col = 0; col < n; ++col) {
    expected[col] = 0.5 + 0.0625 * (double)((7 * col) % 19);
  }
  size_t pos = 0;
  for (int32_t col = 0; col < n; ++col) {
    ap[col] = (int32_t)pos;
    int32_t row_start = 0;
    int32_t row_limit = 0;
    if (col < core) {
      row_limit = core;
    } else if (col < trail_break_col0) {
      row_limit = core;
    } else if (col == trail_break_col0) {
      row_start = lead0;
      row_limit = core;
    } else {
      row_start = lead0 + lead1;
      row_limit = core;
    }
    for (int32_t row = row_start; row < row_limit; ++row) {
      const size_t p = pos++;
      ai[p] = row;
      ax0[p] = row == col
        ? 30.0 + 0.015 * (double)col
        : 0.0002 * (1.0 + (double)((row + 7 * col) % 29));
      ax1[p] = ax0[p] + (row == col
        ? 0.03 * (double)((col % 5) + 1)
        : 1.0e-5 * (double)(((row + col) % 7) - 3));
    }
    if (col >= core) {
      const size_t p = pos++;
      ai[p] = col;
      ax0[p] = 22.0 + 0.02 * (double)col;
      ax1[p] = ax0[p] + 0.02 * (double)((col % 3) + 1);
    }
  }
  ap[n] = (int32_t)pos;
  if (pos != nnz) {
    fprintf(stderr, "unexpected batched-supernode fixture nnz: %zu/%zu\n",
            pos, nnz);
    free(ap);
    free(ai);
    free(ax0);
    free(ax1);
    free(b);
    free(x);
    free(expected);
    return 0;
  }
  for (int32_t col = 0; col < n; ++col) {
    for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
      b[ai[p]] += ax1[p] * expected[col];
    }
  }

  kls_solver *solver = NULL;
  kls_solver *checked_solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.threads = 1;
  options.ordering = KLS_ORDERING_NATURAL;
  options.orientation = KLS_ORIENTATION_NORMAL;
  options.use_btf = 0;
  options.scale = -1;
  options.pivot_tolerance = 0.01;
  options.static_pivoting = 0;

  const char *saved_row_env_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_row_env =
    saved_row_env_value != NULL ? strdup(saved_row_env_value) : NULL;
  const int had_saved_row_env = saved_row_env_value != NULL;
  const char *saved_cblas_env_value = getenv("KLS_ENABLE_CBLAS_SUPERNODE");
  char *saved_cblas_env =
    saved_cblas_env_value != NULL ? strdup(saved_cblas_env_value) : NULL;
  const int had_saved_cblas_env = saved_cblas_env_value != NULL;
  const char *saved_checked_env_value =
    getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  char *saved_checked_env =
    saved_checked_env_value != NULL ? strdup(saved_checked_env_value) : NULL;
  const int had_saved_checked_env = saved_checked_env_value != NULL;
  const char *saved_native_env_value =
    getenv("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR");
  char *saved_native_env =
    saved_native_env_value != NULL ? strdup(saved_native_env_value) : NULL;
  const int had_saved_native_env = saved_native_env_value != NULL;

  int ok = 1;
  if (had_saved_row_env && saved_row_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_ROW_REFACTOR\n");
    ok = 0;
  }
  if (had_saved_cblas_env && saved_cblas_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_CBLAS_SUPERNODE\n");
    ok = 0;
  }
  if (had_saved_checked_env && saved_checked_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_CHECKED_ROW_REFACTOR\n");
    ok = 0;
  }
  if (had_saved_native_env && saved_native_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR\n");
    ok = 0;
  }
  if (!require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, n, ap, ai, 0,
                                        &options),
                        "analyze batched compact supernode")) ok = 0;
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor batched compact supernode")) ok = 0;
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CBLAS_SUPERNODE", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CBLAS_SUPERNODE=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR=1");
    ok = 0;
  }
  if (ok && !require_ok(kls_refactor(solver, ax1),
                        "refactor batched compact supernode")) ok = 0;
  if (had_saved_row_env && saved_row_env != NULL) {
    if (setenv("KLS_ENABLE_ROW_REFACTOR", saved_row_env, 1) != 0) {
      perror("restore KLS_ENABLE_ROW_REFACTOR");
      ok = 0;
    }
  } else if (!had_saved_row_env) {
    if (unsetenv("KLS_ENABLE_ROW_REFACTOR") != 0) {
      perror("unsetenv KLS_ENABLE_ROW_REFACTOR");
      ok = 0;
    }
  }
  if (had_saved_cblas_env && saved_cblas_env != NULL) {
    if (setenv("KLS_ENABLE_CBLAS_SUPERNODE", saved_cblas_env, 1) != 0) {
      perror("restore KLS_ENABLE_CBLAS_SUPERNODE");
      ok = 0;
    }
  } else if (!had_saved_cblas_env) {
    if (unsetenv("KLS_ENABLE_CBLAS_SUPERNODE") != 0) {
      perror("unsetenv KLS_ENABLE_CBLAS_SUPERNODE");
      ok = 0;
    }
  }
  if (!restore_env_value("KLS_ENABLE_CHECKED_ROW_REFACTOR",
                         had_saved_checked_env, saved_checked_env)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR",
                         had_saved_native_env, saved_native_env)) {
    ok = 0;
  }
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve batched compact supernode")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats batched compact supernode")) {
    ok = 0;
  }
  if (ok &&
      (stats.row_refactor_last_compact_supernode_batch != 1 ||
       stats.row_refactor_compact_supernode_batch_count < 1 ||
       stats.row_refactor_compact_supernode_batch_rows < mid ||
       stats.row_refactor_compact_supernode_batch_dep_rows <
         (int64_t)(lead0 + lead1) * (int64_t)mid ||
       stats.row_refactor_compact_supernode_batch_entries <= 0 ||
       stats.row_refactor_last_compact_dense_panel_direct_input_rows <
         stats.row_refactor_compact_supernode_batch_rows ||
       stats.row_refactor_last_segment_target_input_rows <
         stats.row_refactor_last_compact_dense_panel_direct_input_rows)) {
    fprintf(stderr,
            "unexpected batched compact-supernode stats: build_cblas=%d"
            ", batch=%d/%" PRId64 "/%" PRId64 "/%" PRId64 "/%" PRId64
            ", compact=%d/%" PRId64 ", direct=%" PRId64
            ", target=%" PRId64 ", supernode=%d/%" PRId64 "\n",
            stats.build_has_cblas,
            stats.row_refactor_last_compact_supernode_batch,
            stats.row_refactor_compact_supernode_batch_count,
            stats.row_refactor_compact_supernode_batch_rows,
            stats.row_refactor_compact_supernode_batch_dep_rows,
            stats.row_refactor_compact_supernode_batch_entries,
            stats.row_refactor_last_compact_dense_panel,
            stats.row_refactor_compact_dense_panel_count,
            stats.row_refactor_last_compact_dense_panel_direct_input_rows,
            stats.row_refactor_last_segment_target_input_rows,
            stats.row_refactor_last_compact_supernode_update,
            stats.row_refactor_compact_supernode_update_count);
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
  double *residual = (double *)malloc((size_t)n * sizeof(*residual));
  if (residual == NULL) {
    ok = 0;
  }
  for (int32_t row = 0; row < n; ++row) {
    if (residual != NULL) {
      residual[row] = -b[row];
    }
    if (fabs(b[row]) > max_rhs) {
      max_rhs = fabs(b[row]);
    }
  }
  if (residual != NULL) {
    for (int32_t col = 0; col < n; ++col) {
      for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
        residual[ai[p]] += ax1[p] * x[col];
      }
    }
    for (int32_t row = 0; row < n; ++row) {
      const double residual_abs = fabs(residual[row]);
      if (residual_abs > max_residual) {
        max_residual = residual_abs;
      }
    }
  }
  const double relative_residual =
    max_residual / (max_rhs > 0.0 ? max_rhs : 1.0);
  if (ok && (max_solution_error > 1.0e-8 ||
             relative_residual > 1.0e-10)) {
    fprintf(stderr,
            "unexpected batched compact-supernode accuracy:"
            " max_x_err=%.17g, rel_resid=%.17g\n",
            max_solution_error, relative_residual);
    ok = 0;
  }

  if (ok && !require_ok(kls_create(&checked_solver),
                        "create checked batched compact supernode")) {
    ok = 0;
  }
  if (ok && !require_ok(kls_analyze_csc(checked_solver, KLS_INDEX_INT32,
                                        n, ap, ai, 0, &options),
                        "analyze checked batched compact supernode")) {
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CBLAS_SUPERNODE", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CBLAS_SUPERNODE=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(checked_solver, ax0),
                        "factor checked batched compact supernode base")) {
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR=1");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR=1");
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(checked_solver, ax1),
                        "checked factor batched compact supernode")) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_ROW_REFACTOR",
                         had_saved_row_env, saved_row_env)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_CBLAS_SUPERNODE",
                         had_saved_cblas_env, saved_cblas_env)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_CHECKED_ROW_REFACTOR",
                         had_saved_checked_env, saved_checked_env)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR",
                         had_saved_native_env, saved_native_env)) {
    ok = 0;
  }

  if (ok) {
    memset(x, 0, (size_t)n * sizeof(*x));
  }
  if (ok && !require_ok(kls_solve(checked_solver, 1, b, 0, x, 0),
                        "solve checked batched compact supernode")) {
    ok = 0;
  }
  kls_stats checked_stats;
  checked_stats.struct_size = sizeof(checked_stats);
  if (ok && !require_ok(kls_get_stats(checked_solver, &checked_stats),
                        "stats checked batched compact supernode")) {
    ok = 0;
  }
  if (ok &&
      (checked_stats.row_refactor_last_run != 1 ||
       checked_stats.row_refactor_last_checked != 1 ||
       checked_stats.row_refactor_last_compact_supernode_batch != 1 ||
       checked_stats.row_refactor_compact_supernode_batch_count < 1 ||
       checked_stats.row_refactor_compact_supernode_batch_rows < mid ||
       checked_stats.row_refactor_compact_supernode_batch_dep_rows <
         (int64_t)(lead0 + lead1) * (int64_t)mid ||
       checked_stats.row_refactor_compact_supernode_batch_entries <= 0 ||
       checked_stats.row_refactor_last_compact_dense_panel != 1 ||
       checked_stats.row_refactor_last_compact_dense_panel_blocked != 1 ||
       checked_stats.row_refactor_compact_dense_panel_blocked_run_count < 1 ||
       checked_stats.row_refactor_compact_dense_panel_blocked_rows < mid ||
       checked_stats.row_refactor_compact_dense_panel_blocked_entries <= 0 ||
       checked_stats.row_refactor_native_row_panel_enabled != 1 ||
       checked_stats.row_refactor_last_native_row_panel != 1 ||
       checked_stats.row_refactor_native_row_panel_count < 1 ||
       checked_stats.row_refactor_native_row_panel_rows < mid ||
       checked_stats.row_refactor_native_row_panel_entries <= 0 ||
       checked_stats.row_refactor_native_row_panel_blocked_count < 1 ||
       checked_stats.row_refactor_native_row_panel_blocked_rows < mid ||
       checked_stats.row_refactor_native_row_panel_blocked_entries <= 0 ||
       checked_stats.row_refactor_native_row_panel_fallback_count != 0 ||
       checked_stats.row_refactor_native_row_panel_checked_reject_count != 0 ||
       checked_stats.row_refactor_last_compact_dense_panel_direct_input_rows <
         checked_stats.row_refactor_compact_supernode_batch_rows ||
       checked_stats.row_refactor_last_segment_target_input_rows <
         checked_stats.row_refactor_last_compact_dense_panel_direct_input_rows)) {
    fprintf(stderr,
            "unexpected checked batched compact-supernode stats: checked=%d/%d"
            ", batch=%d/%" PRId64 "/%" PRId64 "/%" PRId64 "/%" PRId64
            ", compact=%d/%" PRId64
            ", blocked=%d/%" PRId64 "/%" PRId64 "/%" PRId64
            ", native=%d/%d/%" PRId64 "/%" PRId64 "/%" PRId64
            "/%" PRId64 "/%" PRId64 "/%" PRId64 "/%" PRId64 "/%" PRId64
            ", direct=%" PRId64
            ", target=%" PRId64 "\n",
            checked_stats.row_refactor_last_run,
            checked_stats.row_refactor_last_checked,
            checked_stats.row_refactor_last_compact_supernode_batch,
            checked_stats.row_refactor_compact_supernode_batch_count,
            checked_stats.row_refactor_compact_supernode_batch_rows,
            checked_stats.row_refactor_compact_supernode_batch_dep_rows,
            checked_stats.row_refactor_compact_supernode_batch_entries,
            checked_stats.row_refactor_last_compact_dense_panel,
            checked_stats.row_refactor_compact_dense_panel_count,
            checked_stats.row_refactor_last_compact_dense_panel_blocked,
            checked_stats.row_refactor_compact_dense_panel_blocked_run_count,
            checked_stats.row_refactor_compact_dense_panel_blocked_rows,
            checked_stats.row_refactor_compact_dense_panel_blocked_entries,
            checked_stats.row_refactor_native_row_panel_enabled,
            checked_stats.row_refactor_last_native_row_panel,
            checked_stats.row_refactor_native_row_panel_count,
            checked_stats.row_refactor_native_row_panel_rows,
            checked_stats.row_refactor_native_row_panel_entries,
            checked_stats.row_refactor_native_row_panel_blocked_count,
            checked_stats.row_refactor_native_row_panel_blocked_rows,
            checked_stats.row_refactor_native_row_panel_blocked_entries,
            checked_stats.row_refactor_native_row_panel_fallback_count,
            checked_stats.row_refactor_native_row_panel_checked_reject_count,
            checked_stats.row_refactor_last_compact_dense_panel_direct_input_rows,
            checked_stats.row_refactor_last_segment_target_input_rows);
    ok = 0;
  }

  max_solution_error = 0.0;
  max_residual = 0.0;
  max_rhs = 0.0;
  for (int32_t row = 0; row < n; ++row) {
    const double err = fabs(x[row] - expected[row]);
    if (err > max_solution_error) {
      max_solution_error = err;
    }
    if (residual != NULL) {
      residual[row] = -b[row];
    }
    if (fabs(b[row]) > max_rhs) {
      max_rhs = fabs(b[row]);
    }
  }
  if (residual != NULL) {
    for (int32_t col = 0; col < n; ++col) {
      for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
        residual[ai[p]] += ax1[p] * x[col];
      }
    }
    for (int32_t row = 0; row < n; ++row) {
      const double residual_abs = fabs(residual[row]);
      if (residual_abs > max_residual) {
        max_residual = residual_abs;
      }
    }
  }
  const double checked_relative_residual =
    max_residual / (max_rhs > 0.0 ? max_rhs : 1.0);
  if (ok && (max_solution_error > 1.0e-8 ||
             checked_relative_residual > 1.0e-10)) {
    fprintf(stderr,
            "unexpected checked batched compact-supernode accuracy:"
            " max_x_err=%.17g, rel_resid=%.17g\n",
            max_solution_error, checked_relative_residual);
    ok = 0;
  }

  kls_destroy(solver);
  kls_destroy(checked_solver);
  free(saved_row_env);
  free(saved_cblas_env);
  free(saved_checked_env);
  free(saved_native_env);
  free(ap);
  free(ai);
  free(ax0);
  free(ax1);
  free(b);
  free(x);
  free(expected);
  free(residual);
  return ok;
}

static int test_batched_compact_supernode_subrange_update_probe(void) {
  const int32_t lead = 40;
  const int32_t cold = 20;
  const int32_t hot = 60;
  const int32_t mid = cold + hot;
  const int32_t trail = 20;
  const int32_t core = lead + mid;
  const int32_t n = core + trail;
  const int32_t trail_break_col = n - 1;
  const size_t nnz =
    (size_t)lead * (size_t)(lead + hot) +
    (size_t)mid * (size_t)mid +
    (size_t)(trail - 1) * ((size_t)core + 1u) +
    (size_t)mid + 1u;
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

  for (int32_t col = 0; col < n; ++col) {
    expected[col] = 0.375 + 0.05 * (double)((11 * col) % 23);
  }
  size_t pos = 0;
  for (int32_t col = 0; col < n; ++col) {
    ap[col] = (int32_t)pos;
    if (col < lead) {
      for (int32_t row = 0; row < lead; ++row) {
        const size_t p = pos++;
        ai[p] = row;
        ax0[p] = row == col
          ? 31.0 + 0.012 * (double)col
          : 0.00018 * (1.0 + (double)((row + 11 * col) % 31));
        ax1[p] = ax0[p] + (row == col
          ? 0.025 * (double)((col % 5) + 1)
          : 1.0e-5 * (double)(((row + col) % 7) - 3));
      }
      for (int32_t row = lead + cold; row < core; ++row) {
        const size_t p = pos++;
        ai[p] = row;
        ax0[p] =
          0.00018 * (1.0 + (double)((row + 11 * col) % 31));
        ax1[p] =
          ax0[p] + 1.0e-5 * (double)(((row + col) % 7) - 3);
      }
    } else if (col < core) {
      for (int32_t row = lead; row < core; ++row) {
        const size_t p = pos++;
        ai[p] = row;
        ax0[p] = row == col
          ? 31.0 + 0.012 * (double)col
          : 0.00018 * (1.0 + (double)((row + 11 * col) % 31));
        ax1[p] = ax0[p] + (row == col
          ? 0.025 * (double)((col % 5) + 1)
          : 1.0e-5 * (double)(((row + col) % 7) - 3));
      }
    } else {
      int32_t row_start = 0;
      int32_t row_limit = core;
      if (col == trail_break_col) {
        row_start = lead;
      }
      for (int32_t row = row_start; row < row_limit; ++row) {
        const size_t p = pos++;
        ai[p] = row;
        ax0[p] =
          0.00018 * (1.0 + (double)((row + 11 * col) % 31));
        ax1[p] =
          ax0[p] + 1.0e-5 * (double)(((row + col) % 7) - 3);
      }
      const size_t p = pos++;
      ai[p] = col;
      ax0[p] = 23.0 + 0.018 * (double)col;
      ax1[p] = ax0[p] + 0.018 * (double)((col % 3) + 1);
    }
  }
  ap[n] = (int32_t)pos;
  if (pos != nnz) {
    fprintf(stderr,
            "unexpected subrange batched-supernode fixture nnz: %zu/%zu\n",
            pos, nnz);
    free(ap);
    free(ai);
    free(ax0);
    free(ax1);
    free(b);
    free(x);
    free(expected);
    return 0;
  }
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

  const char *saved_row_env_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_row_env =
    saved_row_env_value != NULL ? strdup(saved_row_env_value) : NULL;
  const int had_saved_row_env = saved_row_env_value != NULL;
  const char *saved_cblas_env_value = getenv("KLS_ENABLE_CBLAS_SUPERNODE");
  char *saved_cblas_env =
    saved_cblas_env_value != NULL ? strdup(saved_cblas_env_value) : NULL;
  const int had_saved_cblas_env = saved_cblas_env_value != NULL;
  const char *saved_native_env_value =
    getenv("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR");
  char *saved_native_env =
    saved_native_env_value != NULL ? strdup(saved_native_env_value) : NULL;
  const int had_saved_native_env = saved_native_env_value != NULL;

  int ok = 1;
  if (had_saved_row_env && saved_row_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_ROW_REFACTOR\n");
    ok = 0;
  }
  if (had_saved_cblas_env && saved_cblas_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_CBLAS_SUPERNODE\n");
    ok = 0;
  }
  if (had_saved_native_env && saved_native_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR\n");
    ok = 0;
  }
  if (!require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, n, ap, ai, 0,
                                        &options),
                        "analyze subrange batched compact supernode")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor subrange batched compact supernode")) ok = 0;
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CBLAS_SUPERNODE", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CBLAS_SUPERNODE=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR", "auto", 1) != 0) {
    perror("setenv KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR=auto");
    ok = 0;
  }
  if (ok && !require_ok(kls_refactor(solver, ax1),
                        "refactor subrange batched compact supernode")) {
    ok = 0;
  }
  if (had_saved_row_env && saved_row_env != NULL) {
    if (setenv("KLS_ENABLE_ROW_REFACTOR", saved_row_env, 1) != 0) {
      perror("restore KLS_ENABLE_ROW_REFACTOR");
      ok = 0;
    }
  } else if (!had_saved_row_env) {
    if (unsetenv("KLS_ENABLE_ROW_REFACTOR") != 0) {
      perror("unsetenv KLS_ENABLE_ROW_REFACTOR");
      ok = 0;
    }
  }
  if (had_saved_cblas_env && saved_cblas_env != NULL) {
    if (setenv("KLS_ENABLE_CBLAS_SUPERNODE", saved_cblas_env, 1) != 0) {
      perror("restore KLS_ENABLE_CBLAS_SUPERNODE");
      ok = 0;
    }
  } else if (!had_saved_cblas_env) {
    if (unsetenv("KLS_ENABLE_CBLAS_SUPERNODE") != 0) {
      perror("unsetenv KLS_ENABLE_CBLAS_SUPERNODE");
      ok = 0;
    }
  }
  if (!restore_env_value("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR",
                         had_saved_native_env, saved_native_env)) {
    ok = 0;
  }
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve subrange batched compact supernode")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats subrange batched compact supernode")) {
    ok = 0;
  }
  if (ok &&
      (stats.row_refactor_last_compact_supernode_batch != 1 ||
       stats.row_refactor_compact_supernode_batch_count != 1 ||
       stats.row_refactor_compact_supernode_batch_rows != hot ||
       stats.row_refactor_compact_supernode_batch_dep_rows !=
         (int64_t)lead * (int64_t)hot ||
       stats.row_refactor_compact_supernode_batch_entries <= 0 ||
       stats.row_refactor_compact_supernode_batch_pattern_count != 1 ||
       stats.row_refactor_compact_supernode_batch_pattern_rows != hot ||
       stats.row_refactor_compact_supernode_batch_candidate_count != 1 ||
       stats.row_refactor_compact_supernode_batch_candidate_rows != hot ||
       stats.row_refactor_compact_supernode_batch_candidate_dep_rows !=
         (int64_t)lead * (int64_t)hot ||
       stats.row_refactor_compact_supernode_batch_rejected_work_count != 0 ||
       stats.row_refactor_last_compact_dense_panel_direct_input_rows <
         stats.row_refactor_compact_supernode_batch_rows ||
       stats.row_refactor_last_segment_target_input_rows <
         stats.row_refactor_last_compact_dense_panel_direct_input_rows)) {
    fprintf(stderr,
            "unexpected subrange batched compact-supernode stats: build_cblas=%d"
            ", batch=%d/%" PRId64 "/%" PRId64 "/%" PRId64 "/%" PRId64
            ", pattern=%" PRId64 "/%" PRId64
            ", candidate=%" PRId64 "/%" PRId64 "/%" PRId64
            ", rejected=%" PRId64
            ", compact=%d/%" PRId64 ", direct=%" PRId64
            ", target=%" PRId64 "\n",
            stats.build_has_cblas,
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
            stats.row_refactor_compact_supernode_batch_rejected_work_count,
            stats.row_refactor_last_compact_dense_panel,
            stats.row_refactor_compact_dense_panel_count,
            stats.row_refactor_last_compact_dense_panel_direct_input_rows,
            stats.row_refactor_last_segment_target_input_rows);
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
  double *residual = (double *)malloc((size_t)n * sizeof(*residual));
  if (residual == NULL) {
    ok = 0;
  }
  for (int32_t row = 0; row < n; ++row) {
    if (residual != NULL) {
      residual[row] = -b[row];
    }
    if (fabs(b[row]) > max_rhs) {
      max_rhs = fabs(b[row]);
    }
  }
  if (residual != NULL) {
    for (int32_t col = 0; col < n; ++col) {
      for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
        residual[ai[p]] += ax1[p] * x[col];
      }
    }
    for (int32_t row = 0; row < n; ++row) {
      const double residual_abs = fabs(residual[row]);
      if (residual_abs > max_residual) {
        max_residual = residual_abs;
      }
    }
  }
  const double relative_residual =
    max_residual / (max_rhs > 0.0 ? max_rhs : 1.0);
  if (ok && (max_solution_error > 1.0e-8 ||
             relative_residual > 1.0e-10)) {
    fprintf(stderr,
            "unexpected subrange batched compact-supernode accuracy:"
            " max_x_err=%.17g, rel_resid=%.17g\n",
            max_solution_error, relative_residual);
    ok = 0;
  }

  kls_destroy(solver);
  free(saved_row_env);
  free(saved_cblas_env);
  free(saved_native_env);
  free(ap);
  free(ai);
  free(ax0);
  free(ax1);
  free(b);
  free(x);
  free(expected);
  free(residual);
  return ok;
}

static int test_egraph_retired_supernode_switch(void) {
  const int32_t n = 15000;
  const int32_t panel_width = 16;
  const int32_t consumer_width = 32;
  const int32_t block_width = panel_width + consumer_width;
  size_t nnz = 0;
  for (int32_t start = 0; start < n; start += block_width) {
    const int32_t panel_end =
      start + panel_width < n ? start + panel_width : n;
    const int32_t block_end =
      start + block_width < n ? start + block_width : n;
    const int32_t panel = panel_end - start;
    const int32_t consumers = block_end - panel_end;
    nnz += (size_t)panel * (size_t)(block_end - start);
    nnz += (size_t)consumers * ((size_t)panel + 1u);
  }

  int32_t *ap = (int32_t *)malloc(((size_t)n + 1u) * sizeof(*ap));
  int32_t *ai = (int32_t *)malloc(nnz * sizeof(*ai));
  double *ax0 = (double *)malloc(nnz * sizeof(*ax0));
  double *ax1 = (double *)malloc(nnz * sizeof(*ax1));
  double *b = (double *)calloc((size_t)n, sizeof(*b));
  double *x = (double *)calloc((size_t)n, sizeof(*x));
  double *expected = (double *)malloc((size_t)n * sizeof(*expected));
  double *residual = (double *)malloc((size_t)n * sizeof(*residual));
  if (ap == NULL || ai == NULL || ax0 == NULL || ax1 == NULL ||
      b == NULL || x == NULL || expected == NULL || residual == NULL) {
    free(ap);
    free(ai);
    free(ax0);
    free(ax1);
    free(b);
    free(x);
    free(expected);
    free(residual);
    return 0;
  }

  for (int32_t i = 0; i < n; ++i) {
    expected[i] = 0.625 + 0.03125 * (double)((13 * i) % 29);
  }

  size_t pos = 0;
  for (int32_t start = 0; start < n; start += block_width) {
    const int32_t panel_end =
      start + panel_width < n ? start + panel_width : n;
    const int32_t block_end =
      start + block_width < n ? start + block_width : n;
    for (int32_t col = start; col < panel_end; ++col) {
      ap[col] = (int32_t)pos;
      for (int32_t row = start; row < block_end; ++row) {
        const size_t p = pos++;
        ai[p] = row;
        ax0[p] = row == col
          ? 30.0 + 0.0001 * (double)col
          : 1.0e-4 * (double)(1 + ((row + 7 * col) % 17));
        ax1[p] = ax0[p] + (row == col
          ? 0.02 * (double)((col % 5) + 1)
          : 1.0e-6 * (double)(((row + col) % 5) - 2));
      }
    }
    for (int32_t col = panel_end; col < block_end; ++col) {
      ap[col] = (int32_t)pos;
      for (int32_t row = start; row < panel_end; ++row) {
        const size_t p = pos++;
        ai[p] = row;
        ax0[p] = 1.0e-4 * (double)(1 + ((row + 11 * col) % 19));
        ax1[p] =
          ax0[p] + 1.0e-6 * (double)(((row + 3 * col) % 7) - 3);
      }
      const size_t p = pos++;
      ai[p] = col;
      ax0[p] = 25.0 + 0.0001 * (double)col;
      ax1[p] = ax0[p] + 0.015 * (double)((col % 3) + 1);
    }
  }
  ap[n] = (int32_t)pos;
  if (pos != nnz) {
    fprintf(stderr, "unexpected EGraph blocked fixture nnz: %zu/%zu\n",
            pos, nnz);
    free(ap);
    free(ai);
    free(ax0);
    free(ax1);
    free(b);
    free(x);
    free(expected);
    free(residual);
    return 0;
  }
  for (int32_t col = 0; col < n; ++col) {
    for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
      b[ai[p]] += ax1[p] * expected[col];
    }
  }

  const char *saved_row_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_row = saved_row_value != NULL ? strdup(saved_row_value) : NULL;
  const int had_row = saved_row_value != NULL;
  const char *saved_checked_value =
    getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  char *saved_checked =
    saved_checked_value != NULL ? strdup(saved_checked_value) : NULL;
  const int had_checked = saved_checked_value != NULL;
  const char *saved_egraph_value =
    getenv("KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES");
  char *saved_egraph =
    saved_egraph_value != NULL ? strdup(saved_egraph_value) : NULL;
  const int had_egraph = saved_egraph_value != NULL;
  const char *saved_cblas_value = getenv("KLS_ENABLE_CBLAS_SUPERNODE");
  char *saved_cblas =
    saved_cblas_value != NULL ? strdup(saved_cblas_value) : NULL;
  const int had_cblas = saved_cblas_value != NULL;
  const char *saved_first_value = getenv("KLS_ENABLE_KLS_FIRST_FACTOR");
  char *saved_first =
    saved_first_value != NULL ? strdup(saved_first_value) : NULL;
  const int had_first = saved_first_value != NULL;

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.threads = 4;
  options.ordering = KLS_ORDERING_NATURAL;
  options.orientation = KLS_ORIENTATION_AUTO;
  options.use_btf = 0;
  options.scale = -1;
  options.pivot_tolerance = 0.001;
  options.static_pivoting = 0;

  int ok = 1;

  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=1");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CBLAS_SUPERNODE", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CBLAS_SUPERNODE=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_KLS_FIRST_FACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_KLS_FIRST_FACTOR=0");
    ok = 0;
  }

  /* The production EGraph work floor would route this small fixture to the
     serial mapped path; pin it to zero to exercise EGraph after retiring
     the supernode switch. The old switch must no longer alter execution. */
  if (ok && setenv("KLS_EGRAPH_REFACTOR_FLOOR", "0", 1) != 0) {
    perror("setenv KLS_EGRAPH_REFACTOR_FLOOR=0");
    ok = 0;
  }

  if (!require_ok(kls_create(&solver), "create EGraph blocked")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, n, ap, ai, 0,
                                        &options),
                        "analyze EGraph blocked")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor EGraph blocked base")) ok = 0;
  if (ok && !require_ok(kls_refactor(solver, ax1),
                        "refactor EGraph blocked")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve EGraph blocked")) ok = 0;

  kls_stats stats;
  memset(&stats, 0xa5, sizeof(stats));
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats EGraph blocked")) {
    ok = 0;
  }

  if (ok &&
      (stats.refactor_supernode_consumer_plan_cached_panel_count != 0 ||
       stats.refactor_supernode_consumer_plan_group_l_built != 0 ||
       stats.refactor_btf_scalar_run_group_built != 0 ||
       stats.refactor_last_btf_scalar_run_group_waits != 0 ||
       stats.refactor_supernode_algorithm5_candidate_run_count != 0 ||
       stats.refactor_supernode_consumer_run_count != 0 ||
       stats.refactor_last_btf_scalar_run_candidates != 0 ||
       stats.refactor_u_supernode_pattern_count != 0 ||
       stats.refactor_u_supernode_value_dense_entries != 0 ||
       stats.refactor_last_u_supernode_l_update_runs != 0 ||
       stats.refactor_l_pattern_columns != 0 ||
       stats.refactor_l_adjacent_run_count != 0 ||
       stats.refactor_l_contiguous_suffix_columns != 0 ||
       stats.refactor_stream_dependency_entries != 0.0 ||
       stats.refactor_stream_pivot_entries != 0.0 ||
       stats.refactor_stream_output_entries != 0.0)) {
    fprintf(stderr, "retired experiment statistics must remain zero\n");
    ok = 0;
  }

  if (ok && (stats.row_refactor_last_run != 0 ||
             stats.refactor_supernode_panel_count != 0 ||
             stats.refactor_supernode_panel_used_count != 0 ||
             stats.refactor_last_supernode_update_runs != 0 ||
             stats.refactor_last_supernode_blocked_update_runs != 0 ||
             stats.refactor_last_supernode_cblas_update_runs != 0 ||
             stats.refactor_supernode_cached_probe_attempt_count != 0 ||
             stats.refactor_last_ready_queue_columns != 0)) {
    fprintf(stderr, "retired EGraph engine statistics must remain zero\n");
    ok = 0;
  }

  double max_solution_error = 0.0;
  for (int32_t i = 0; i < n; ++i) {
    const double err = fabs(x[i] - expected[i]);
    if (err > max_solution_error) {
      max_solution_error = err;
    }
    residual[i] = -b[i];
  }
  double max_rhs = 0.0;
  for (int32_t i = 0; i < n; ++i) {
    const double rhs_abs = fabs(b[i]);
    if (rhs_abs > max_rhs) {
      max_rhs = rhs_abs;
    }
  }
  for (int32_t col = 0; col < n; ++col) {
    for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
      residual[ai[p]] += ax1[p] * x[col];
    }
  }
  double max_residual = 0.0;
  for (int32_t i = 0; i < n; ++i) {
    const double residual_abs = fabs(residual[i]);
    if (residual_abs > max_residual) {
      max_residual = residual_abs;
    }
  }
  const double relative_residual =
    max_residual / (max_rhs > 0.0 ? max_rhs : 1.0);
  if (ok && (max_solution_error > 1.0e-8 ||
             relative_residual > 1.0e-10)) {
    fprintf(stderr,
            "unexpected EGraph blocked accuracy:"
            " max_x_err=%.17g, rel_resid=%.17g\n",
            max_solution_error, relative_residual);
    ok = 0;
  }

  if (unsetenv("KLS_EGRAPH_REFACTOR_FLOOR") != 0) {
    perror("unsetenv KLS_EGRAPH_REFACTOR_FLOOR");
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_ROW_REFACTOR", had_row, saved_row)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_CHECKED_ROW_REFACTOR", had_checked,
                         saved_checked)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES", had_egraph,
                         saved_egraph)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_CBLAS_SUPERNODE", had_cblas,
                         saved_cblas)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_KLS_FIRST_FACTOR", had_first,
                         saved_first)) {
    ok = 0;
  }

  free(saved_row);
  free(saved_checked);
  free(saved_egraph);
  free(saved_cblas);
  free(saved_first);
  kls_destroy(solver);
  free(ap);
  free(ai);
  free(ax0);
  free(ax1);
  free(b);
  free(x);
  free(expected);
  free(residual);
  return ok;
}

static int test_ragged_batched_compact_supernode_update_probe(void) {
  const int32_t lead = 48;
  const int32_t mid = 64;
  const int32_t n = lead + mid;
  const size_t capacity = (size_t)n * (size_t)n;
  int32_t *ap = (int32_t *)malloc(((size_t)n + 1u) * sizeof(*ap));
  int32_t *ai = (int32_t *)malloc(capacity * sizeof(*ai));
  double *ax0 = (double *)malloc(capacity * sizeof(*ax0));
  double *ax1 = (double *)malloc(capacity * sizeof(*ax1));
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

  for (int32_t col = 0; col < n; ++col) {
    expected[col] = 0.625 + 0.04 * (double)((13 * col) % 29);
  }
  size_t pos = 0;
  for (int32_t col = 0; col < n; ++col) {
    ap[col] = (int32_t)pos;
    if (col < lead) {
      for (int32_t row = 0; row < lead; ++row) {
        const size_t p = pos++;
        ai[p] = row;
        ax0[p] = row == col
          ? 34.0 + 0.01 * (double)col
          : 0.00016 * (1.0 + (double)((row + 5 * col) % 37));
        ax1[p] = ax0[p] + (row == col
          ? 0.02 * (double)((col % 7) + 1)
          : 1.0e-5 * (double)(((row + col) % 9) - 4));
      }
      for (int32_t row = lead; row < n; ++row) {
        const int32_t suffix_start = ((row - lead) % 4) * 8;
        if (col >= suffix_start) {
          const size_t p = pos++;
          ai[p] = row;
          ax0[p] =
            0.00016 * (1.0 + (double)((row + 5 * col) % 37));
          ax1[p] =
            ax0[p] + 1.0e-5 * (double)(((row + col) % 9) - 4);
        }
      }
    } else {
      for (int32_t row = lead; row < n; ++row) {
        const size_t p = pos++;
        ai[p] = row;
        ax0[p] = row == col
          ? 33.0 + 0.011 * (double)col
          : 0.00014 * (1.0 + (double)((row + 7 * col) % 41));
        ax1[p] = ax0[p] + (row == col
          ? 0.018 * (double)((col % 5) + 1)
          : 1.0e-5 * (double)(((row + col) % 7) - 3));
      }
    }
  }
  ap[n] = (int32_t)pos;
  if (pos > capacity) {
    fprintf(stderr, "ragged batched-supernode fixture overflowed capacity\n");
    free(ap);
    free(ai);
    free(ax0);
    free(ax1);
    free(b);
    free(x);
    free(expected);
    return 0;
  }
  for (int32_t col = 0; col < n; ++col) {
    for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
      b[ai[p]] += ax1[p] * expected[col];
    }
  }

  kls_solver *solver = NULL;
  kls_solver *checked_solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.threads = 1;
  options.ordering = KLS_ORDERING_NATURAL;
  options.orientation = KLS_ORIENTATION_NORMAL;
  options.use_btf = 0;
  options.scale = -1;
  options.pivot_tolerance = 0.01;
  options.static_pivoting = 0;

  const char *saved_row_env_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_row_env =
    saved_row_env_value != NULL ? strdup(saved_row_env_value) : NULL;
  const int had_saved_row_env = saved_row_env_value != NULL;
  const char *saved_cblas_env_value = getenv("KLS_ENABLE_CBLAS_SUPERNODE");
  char *saved_cblas_env =
    saved_cblas_env_value != NULL ? strdup(saved_cblas_env_value) : NULL;
  const int had_saved_cblas_env = saved_cblas_env_value != NULL;
  const char *saved_checked_env_value =
    getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  char *saved_checked_env =
    saved_checked_env_value != NULL ? strdup(saved_checked_env_value) : NULL;
  const int had_saved_checked_env = saved_checked_env_value != NULL;
  const char *saved_native_env_value =
    getenv("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR");
  char *saved_native_env =
    saved_native_env_value != NULL ? strdup(saved_native_env_value) : NULL;
  const int had_saved_native_env = saved_native_env_value != NULL;

  int ok = 1;
  if (had_saved_row_env && saved_row_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_ROW_REFACTOR\n");
    ok = 0;
  }
  if (had_saved_cblas_env && saved_cblas_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_CBLAS_SUPERNODE\n");
    ok = 0;
  }
  if (had_saved_checked_env && saved_checked_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_CHECKED_ROW_REFACTOR\n");
    ok = 0;
  }
  if (had_saved_native_env && saved_native_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR\n");
    ok = 0;
  }
  if (!require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, n, ap, ai, 0,
                                        &options),
                        "analyze ragged batched compact supernode")) ok = 0;
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor ragged batched compact supernode")) ok = 0;
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CBLAS_SUPERNODE", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CBLAS_SUPERNODE=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR", "auto", 1) != 0) {
    perror("setenv KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR=auto");
    ok = 0;
  }
  if (ok && !require_ok(kls_refactor(solver, ax1),
                        "refactor ragged batched compact supernode")) {
    ok = 0;
  }
  if (had_saved_row_env && saved_row_env != NULL) {
    if (setenv("KLS_ENABLE_ROW_REFACTOR", saved_row_env, 1) != 0) {
      perror("restore KLS_ENABLE_ROW_REFACTOR");
      ok = 0;
    }
  } else if (!had_saved_row_env) {
    if (unsetenv("KLS_ENABLE_ROW_REFACTOR") != 0) {
      perror("unsetenv KLS_ENABLE_ROW_REFACTOR");
      ok = 0;
    }
  }
  if (had_saved_cblas_env && saved_cblas_env != NULL) {
    if (setenv("KLS_ENABLE_CBLAS_SUPERNODE", saved_cblas_env, 1) != 0) {
      perror("restore KLS_ENABLE_CBLAS_SUPERNODE");
      ok = 0;
    }
  } else if (!had_saved_cblas_env) {
    if (unsetenv("KLS_ENABLE_CBLAS_SUPERNODE") != 0) {
      perror("unsetenv KLS_ENABLE_CBLAS_SUPERNODE");
      ok = 0;
    }
  }
  if (!restore_env_value("KLS_ENABLE_CHECKED_ROW_REFACTOR",
                         had_saved_checked_env, saved_checked_env)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR",
                         had_saved_native_env, saved_native_env)) {
    ok = 0;
  }
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve ragged batched compact supernode")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats ragged batched compact supernode")) {
    ok = 0;
  }
  if (ok &&
      (stats.row_refactor_last_compact_supernode_batch != 1 ||
       stats.row_refactor_compact_supernode_batch_count < 1 ||
       stats.row_refactor_compact_supernode_batch_pattern_count != 0 ||
       stats.row_refactor_compact_supernode_batch_rows < mid ||
       stats.row_refactor_compact_supernode_batch_dep_rows <=
         (int64_t)lead * (int64_t)(mid / 2) ||
       stats.row_refactor_compact_supernode_batch_entries <= 0 ||
       stats.row_refactor_last_compact_dense_panel_direct_input_rows <
         stats.row_refactor_compact_supernode_batch_rows ||
       stats.row_refactor_last_segment_target_input_rows <
         stats.row_refactor_last_compact_dense_panel_direct_input_rows ||
       stats.row_refactor_compact_supernode_batch_candidate_count < 1 ||
       stats.row_refactor_compact_supernode_batch_rejected_work_count != 0)) {
    fprintf(stderr,
            "unexpected ragged batched compact-supernode stats: build_cblas=%d"
            ", batch=%d/%" PRId64 "/%" PRId64 "/%" PRId64 "/%" PRId64
            ", direct=%" PRId64 ", target=%" PRId64
            ", pattern=%" PRId64 "/%" PRId64
            ", candidate=%" PRId64 "/%" PRId64 "/%" PRId64
            ", rejected=%" PRId64 "\n",
            stats.build_has_cblas,
            stats.row_refactor_last_compact_supernode_batch,
            stats.row_refactor_compact_supernode_batch_count,
            stats.row_refactor_compact_supernode_batch_rows,
            stats.row_refactor_compact_supernode_batch_dep_rows,
            stats.row_refactor_compact_supernode_batch_entries,
            stats.row_refactor_last_compact_dense_panel_direct_input_rows,
            stats.row_refactor_last_segment_target_input_rows,
            stats.row_refactor_compact_supernode_batch_pattern_count,
            stats.row_refactor_compact_supernode_batch_pattern_rows,
            stats.row_refactor_compact_supernode_batch_candidate_count,
            stats.row_refactor_compact_supernode_batch_candidate_rows,
            stats.row_refactor_compact_supernode_batch_candidate_dep_rows,
            stats.row_refactor_compact_supernode_batch_rejected_work_count);
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
  double *residual = (double *)malloc((size_t)n * sizeof(*residual));
  if (residual == NULL) {
    ok = 0;
  }
  for (int32_t row = 0; row < n; ++row) {
    if (residual != NULL) {
      residual[row] = -b[row];
    }
    if (fabs(b[row]) > max_rhs) {
      max_rhs = fabs(b[row]);
    }
  }
  if (residual != NULL) {
    for (int32_t col = 0; col < n; ++col) {
      for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
        residual[ai[p]] += ax1[p] * x[col];
      }
    }
    for (int32_t row = 0; row < n; ++row) {
      const double residual_abs = fabs(residual[row]);
      if (residual_abs > max_residual) {
        max_residual = residual_abs;
      }
    }
  }
  const double relative_residual =
    max_residual / (max_rhs > 0.0 ? max_rhs : 1.0);
  if (ok && (max_solution_error > 1.0e-8 ||
             relative_residual > 1.0e-10)) {
    fprintf(stderr,
            "unexpected ragged batched compact-supernode accuracy:"
            " max_x_err=%.17g, rel_resid=%.17g\n",
            max_solution_error, relative_residual);
    ok = 0;
  }

  if (ok && !require_ok(kls_create(&checked_solver),
                        "create checked ragged batched compact supernode")) {
    ok = 0;
  }
  if (ok && !require_ok(kls_analyze_csc(checked_solver, KLS_INDEX_INT32,
                                        n, ap, ai, 0, &options),
                        "analyze checked ragged batched compact supernode")) {
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CBLAS_SUPERNODE", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CBLAS_SUPERNODE=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(checked_solver, ax0),
                        "factor checked ragged batched compact supernode base")) {
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR=1");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR", "auto", 1) != 0) {
    perror("setenv KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR=auto");
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(checked_solver, ax1),
                        "checked factor ragged batched compact supernode")) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_ROW_REFACTOR",
                         had_saved_row_env, saved_row_env)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_CBLAS_SUPERNODE",
                         had_saved_cblas_env, saved_cblas_env)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_CHECKED_ROW_REFACTOR",
                         had_saved_checked_env, saved_checked_env)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR",
                         had_saved_native_env, saved_native_env)) {
    ok = 0;
  }

  if (ok) {
    memset(x, 0, (size_t)n * sizeof(*x));
  }
  if (ok && !require_ok(kls_solve(checked_solver, 1, b, 0, x, 0),
                        "solve checked ragged batched compact supernode")) {
    ok = 0;
  }
  kls_stats checked_stats;
  checked_stats.struct_size = sizeof(checked_stats);
  if (ok && !require_ok(kls_get_stats(checked_solver, &checked_stats),
                        "stats checked ragged batched compact supernode")) {
    ok = 0;
  }
  if (ok &&
      (checked_stats.row_refactor_last_run != 1 ||
       checked_stats.row_refactor_last_checked != 1 ||
       checked_stats.row_refactor_last_compact_supernode_batch != 1 ||
       checked_stats.row_refactor_compact_supernode_batch_count < 1 ||
       checked_stats.row_refactor_compact_supernode_batch_pattern_count != 0 ||
       checked_stats.row_refactor_compact_supernode_batch_rows < mid ||
       checked_stats.row_refactor_compact_supernode_batch_dep_rows <=
         (int64_t)lead * (int64_t)(mid / 2) ||
       checked_stats.row_refactor_compact_supernode_batch_entries <= 0 ||
       checked_stats.row_refactor_last_compact_dense_panel_direct_input_rows <
         checked_stats.row_refactor_compact_supernode_batch_rows ||
       checked_stats.row_refactor_last_segment_target_input_rows <
         checked_stats.row_refactor_last_compact_dense_panel_direct_input_rows ||
       checked_stats.row_refactor_compact_supernode_batch_candidate_count < 1 ||
       checked_stats.row_refactor_last_compact_dense_panel != 1)) {
    fprintf(stderr,
            "unexpected checked ragged batched compact-supernode stats:"
            " checked=%d/%d, batch=%d/%" PRId64
            "/%" PRId64 "/%" PRId64 "/%" PRId64
            ", direct=%" PRId64 ", target=%" PRId64
            ", pattern=%" PRId64 "/%" PRId64
            ", candidate=%" PRId64 "/%" PRId64 "/%" PRId64
            ", compact=%d/%" PRId64 "\n",
            checked_stats.row_refactor_last_run,
            checked_stats.row_refactor_last_checked,
            checked_stats.row_refactor_last_compact_supernode_batch,
            checked_stats.row_refactor_compact_supernode_batch_count,
            checked_stats.row_refactor_compact_supernode_batch_rows,
            checked_stats.row_refactor_compact_supernode_batch_dep_rows,
            checked_stats.row_refactor_compact_supernode_batch_entries,
            checked_stats.row_refactor_last_compact_dense_panel_direct_input_rows,
            checked_stats.row_refactor_last_segment_target_input_rows,
            checked_stats.row_refactor_compact_supernode_batch_pattern_count,
            checked_stats.row_refactor_compact_supernode_batch_pattern_rows,
            checked_stats.row_refactor_compact_supernode_batch_candidate_count,
            checked_stats.row_refactor_compact_supernode_batch_candidate_rows,
            checked_stats.row_refactor_compact_supernode_batch_candidate_dep_rows,
            checked_stats.row_refactor_last_compact_dense_panel,
            checked_stats.row_refactor_compact_dense_panel_count);
    ok = 0;
  }

  max_solution_error = 0.0;
  max_residual = 0.0;
  max_rhs = 0.0;
  for (int32_t row = 0; row < n; ++row) {
    const double err = fabs(x[row] - expected[row]);
    if (err > max_solution_error) {
      max_solution_error = err;
    }
    if (residual != NULL) {
      residual[row] = -b[row];
    }
    if (fabs(b[row]) > max_rhs) {
      max_rhs = fabs(b[row]);
    }
  }
  if (residual != NULL) {
    for (int32_t col = 0; col < n; ++col) {
      for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
        residual[ai[p]] += ax1[p] * x[col];
      }
    }
    for (int32_t row = 0; row < n; ++row) {
      const double residual_abs = fabs(residual[row]);
      if (residual_abs > max_residual) {
        max_residual = residual_abs;
      }
    }
  }
  const double checked_relative_residual =
    max_residual / (max_rhs > 0.0 ? max_rhs : 1.0);
  if (ok && (max_solution_error > 1.0e-8 ||
             checked_relative_residual > 1.0e-10)) {
    fprintf(stderr,
            "unexpected checked ragged batched compact-supernode accuracy:"
            " max_x_err=%.17g, rel_resid=%.17g\n",
            max_solution_error, checked_relative_residual);
    ok = 0;
  }

  kls_destroy(solver);
  kls_destroy(checked_solver);
  free(saved_row_env);
  free(saved_cblas_env);
  free(saved_checked_env);
  free(saved_native_env);
  free(ap);
  free(ai);
  free(ax0);
  free(ax1);
  free(b);
  free(x);
  free(expected);
  free(residual);
  return ok;
}

static int test_batch_group_ragged_supernode_update_probe(void) {
  const int32_t lead = 128;
  const int32_t mid = 80;
  const int32_t n = lead + mid;
  const size_t capacity = (size_t)n * (size_t)n;
  int32_t *ap = (int32_t *)malloc(((size_t)n + 1u) * sizeof(*ap));
  int32_t *ai = (int32_t *)malloc(capacity * sizeof(*ai));
  double *ax0 = (double *)malloc(capacity * sizeof(*ax0));
  double *ax1 = (double *)malloc(capacity * sizeof(*ax1));
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

  for (int32_t col = 0; col < n; ++col) {
    expected[col] = 0.5 + 0.045 * (double)((17 * col) % 31);
  }
  size_t pos = 0;
  for (int32_t col = 0; col < n; ++col) {
    ap[col] = (int32_t)pos;
    if (col < lead) {
      for (int32_t row = 0; row < lead; ++row) {
        const size_t p = pos++;
        ai[p] = row;
        ax0[p] = row == col
          ? 36.0 + 0.009 * (double)col
          : 0.00013 * (1.0 + (double)((row + 3 * col) % 43));
        ax1[p] = ax0[p] + (row == col
          ? 0.019 * (double)((col % 7) + 1)
          : 1.0e-5 * (double)(((row + col) % 11) - 5));
      }
      for (int32_t row = lead; row < n; ++row) {
        const int32_t suffix_start = ((row - lead) % 4) * 16;
        if (col >= suffix_start) {
          const size_t p = pos++;
          ai[p] = row;
          ax0[p] =
            0.00013 * (1.0 + (double)((row + 3 * col) % 43));
          ax1[p] =
            ax0[p] + 1.0e-5 * (double)(((row + col) % 11) - 5);
        }
      }
    } else {
      const size_t p = pos++;
      ai[p] = col;
      ax0[p] = 29.0 + 0.012 * (double)col;
      ax1[p] = ax0[p] + 0.017 * (double)((col % 5) + 1);
    }
  }
  ap[n] = (int32_t)pos;
  if (pos > capacity) {
    fprintf(stderr, "batch-group ragged fixture overflowed capacity\n");
    free(ap);
    free(ai);
    free(ax0);
    free(ax1);
    free(b);
    free(x);
    free(expected);
    return 0;
  }
  for (int32_t col = 0; col < n; ++col) {
    for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
      b[ai[p]] += ax1[p] * expected[col];
    }
  }

  kls_solver *solver = NULL;
  kls_solver *checked_solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.threads = 1;
  options.ordering = KLS_ORDERING_NATURAL;
  options.orientation = KLS_ORIENTATION_NORMAL;
  options.use_btf = 0;
  options.scale = -1;
  options.pivot_tolerance = 0.01;
  options.static_pivoting = 0;

  const char *saved_row_env_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_row_env =
    saved_row_env_value != NULL ? strdup(saved_row_env_value) : NULL;
  const int had_saved_row_env = saved_row_env_value != NULL;
  const char *saved_cblas_env_value = getenv("KLS_ENABLE_CBLAS_SUPERNODE");
  char *saved_cblas_env =
    saved_cblas_env_value != NULL ? strdup(saved_cblas_env_value) : NULL;
  const int had_saved_cblas_env = saved_cblas_env_value != NULL;
  const char *saved_checked_env_value =
    getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  char *saved_checked_env =
    saved_checked_env_value != NULL ? strdup(saved_checked_env_value) : NULL;
  const int had_saved_checked_env = saved_checked_env_value != NULL;
  const char *saved_native_env_value =
    getenv("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR");
  char *saved_native_env =
    saved_native_env_value != NULL ? strdup(saved_native_env_value) : NULL;
  const int had_saved_native_env = saved_native_env_value != NULL;

  int ok = 1;
  if (had_saved_row_env && saved_row_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_ROW_REFACTOR\n");
    ok = 0;
  }
  if (had_saved_cblas_env && saved_cblas_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_CBLAS_SUPERNODE\n");
    ok = 0;
  }
  if (had_saved_checked_env && saved_checked_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_CHECKED_ROW_REFACTOR\n");
    ok = 0;
  }
  if (had_saved_native_env && saved_native_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR\n");
    ok = 0;
  }
  if (!require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, n, ap, ai, 0,
                                        &options),
                        "analyze batch-group ragged supernode")) ok = 0;
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor batch-group ragged supernode")) ok = 0;
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CBLAS_SUPERNODE", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CBLAS_SUPERNODE=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR", "auto", 1) != 0) {
    perror("setenv KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR=auto");
    ok = 0;
  }
  if (ok && !require_ok(kls_refactor(solver, ax1),
                        "refactor batch-group ragged supernode")) {
    ok = 0;
  }
  if (had_saved_row_env && saved_row_env != NULL) {
    if (setenv("KLS_ENABLE_ROW_REFACTOR", saved_row_env, 1) != 0) {
      perror("restore KLS_ENABLE_ROW_REFACTOR");
      ok = 0;
    }
  } else if (!had_saved_row_env) {
    if (unsetenv("KLS_ENABLE_ROW_REFACTOR") != 0) {
      perror("unsetenv KLS_ENABLE_ROW_REFACTOR");
      ok = 0;
    }
  }
  if (had_saved_cblas_env && saved_cblas_env != NULL) {
    if (setenv("KLS_ENABLE_CBLAS_SUPERNODE", saved_cblas_env, 1) != 0) {
      perror("restore KLS_ENABLE_CBLAS_SUPERNODE");
      ok = 0;
    }
  } else if (!had_saved_cblas_env) {
    if (unsetenv("KLS_ENABLE_CBLAS_SUPERNODE") != 0) {
      perror("unsetenv KLS_ENABLE_CBLAS_SUPERNODE");
      ok = 0;
    }
  }
  if (!restore_env_value("KLS_ENABLE_CHECKED_ROW_REFACTOR",
                         had_saved_checked_env, saved_checked_env)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR",
                         had_saved_native_env, saved_native_env)) {
    ok = 0;
  }
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve batch-group ragged supernode")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats batch-group ragged supernode")) {
    ok = 0;
  }
  if (ok &&
      (stats.row_refactor_last_compact_supernode_batch != 1 ||
       stats.row_refactor_compact_supernode_batch_count < 1 ||
       stats.row_refactor_compact_supernode_batch_pattern_count != 0 ||
       stats.row_refactor_compact_supernode_batch_rows < mid ||
       stats.row_refactor_compact_supernode_batch_dep_rows <=
         (int64_t)lead * (int64_t)(mid / 2) ||
       stats.row_refactor_compact_supernode_batch_entries <= 0 ||
       stats.row_refactor_last_batch_direct_input_rows <
         stats.row_refactor_compact_supernode_batch_rows ||
       stats.row_refactor_batch_direct_input_rows <
         stats.row_refactor_compact_supernode_batch_rows ||
       stats.row_refactor_last_segment_target_input_rows <
         stats.row_refactor_compact_supernode_batch_rows ||
       stats.row_refactor_segment_target_input_rows <
         stats.row_refactor_compact_supernode_batch_rows ||
       stats.row_refactor_compact_supernode_batch_candidate_count < 1 ||
       stats.row_refactor_compact_supernode_batch_rejected_work_count != 0)) {
    fprintf(stderr,
            "unexpected batch-group ragged stats: build_cblas=%d"
            ", batch=%d/%" PRId64 "/%" PRId64 "/%" PRId64 "/%" PRId64
            ", direct=%" PRId64 "/%" PRId64
            ", target=%" PRId64 "/%" PRId64
            ", pattern=%" PRId64 "/%" PRId64
            ", candidate=%" PRId64 "/%" PRId64 "/%" PRId64
            ", rejected=%" PRId64 "\n",
            stats.build_has_cblas,
            stats.row_refactor_last_compact_supernode_batch,
            stats.row_refactor_compact_supernode_batch_count,
            stats.row_refactor_compact_supernode_batch_rows,
            stats.row_refactor_compact_supernode_batch_dep_rows,
            stats.row_refactor_compact_supernode_batch_entries,
            stats.row_refactor_last_batch_direct_input_rows,
            stats.row_refactor_batch_direct_input_rows,
            stats.row_refactor_last_segment_target_input_rows,
            stats.row_refactor_segment_target_input_rows,
            stats.row_refactor_compact_supernode_batch_pattern_count,
            stats.row_refactor_compact_supernode_batch_pattern_rows,
            stats.row_refactor_compact_supernode_batch_candidate_count,
            stats.row_refactor_compact_supernode_batch_candidate_rows,
            stats.row_refactor_compact_supernode_batch_candidate_dep_rows,
            stats.row_refactor_compact_supernode_batch_rejected_work_count);
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
  double *residual = (double *)malloc((size_t)n * sizeof(*residual));
  if (residual == NULL) {
    ok = 0;
  }
  for (int32_t row = 0; row < n; ++row) {
    if (residual != NULL) {
      residual[row] = -b[row];
    }
    if (fabs(b[row]) > max_rhs) {
      max_rhs = fabs(b[row]);
    }
  }
  if (residual != NULL) {
    for (int32_t col = 0; col < n; ++col) {
      for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
        residual[ai[p]] += ax1[p] * x[col];
      }
    }
    for (int32_t row = 0; row < n; ++row) {
      const double residual_abs = fabs(residual[row]);
      if (residual_abs > max_residual) {
        max_residual = residual_abs;
      }
    }
  }
  const double relative_residual =
    max_residual / (max_rhs > 0.0 ? max_rhs : 1.0);
  if (ok && (max_solution_error > 1.0e-8 ||
             relative_residual > 1.0e-10)) {
    fprintf(stderr,
            "unexpected batch-group ragged accuracy:"
            " max_x_err=%.17g, rel_resid=%.17g\n",
            max_solution_error, relative_residual);
    ok = 0;
  }

  if (ok && !require_ok(kls_create(&checked_solver),
                        "create checked batch-group ragged supernode")) {
    ok = 0;
  }
  if (ok && !require_ok(kls_analyze_csc(checked_solver, KLS_INDEX_INT32,
                                        n, ap, ai, 0, &options),
                        "analyze checked batch-group ragged supernode")) {
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CBLAS_SUPERNODE", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CBLAS_SUPERNODE=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(checked_solver, ax0),
                        "factor checked batch-group ragged base")) {
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR=1");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR", "auto", 1) != 0) {
    perror("setenv KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR=auto");
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(checked_solver, ax1),
                        "checked factor batch-group ragged")) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_ROW_REFACTOR",
                         had_saved_row_env, saved_row_env)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_CBLAS_SUPERNODE",
                         had_saved_cblas_env, saved_cblas_env)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_CHECKED_ROW_REFACTOR",
                         had_saved_checked_env, saved_checked_env)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR",
                         had_saved_native_env, saved_native_env)) {
    ok = 0;
  }

  if (ok) {
    memset(x, 0, (size_t)n * sizeof(*x));
  }
  if (ok && !require_ok(kls_solve(checked_solver, 1, b, 0, x, 0),
                        "solve checked batch-group ragged")) {
    ok = 0;
  }
  kls_stats checked_stats;
  checked_stats.struct_size = sizeof(checked_stats);
  if (ok && !require_ok(kls_get_stats(checked_solver, &checked_stats),
                        "stats checked batch-group ragged")) {
    ok = 0;
  }
  if (ok &&
      (checked_stats.row_refactor_last_run != 1 ||
       checked_stats.row_refactor_last_checked != 1 ||
       checked_stats.row_refactor_last_compact_supernode_batch != 1 ||
       checked_stats.row_refactor_compact_supernode_batch_count < 1 ||
       checked_stats.row_refactor_compact_supernode_batch_pattern_count != 0 ||
       checked_stats.row_refactor_compact_supernode_batch_rows < mid ||
       checked_stats.row_refactor_compact_supernode_batch_dep_rows <=
         (int64_t)lead * (int64_t)(mid / 2) ||
       checked_stats.row_refactor_compact_supernode_batch_entries <= 0 ||
       checked_stats.row_refactor_last_batch_direct_input_rows <
         checked_stats.row_refactor_compact_supernode_batch_rows ||
       checked_stats.row_refactor_batch_direct_input_rows <
         checked_stats.row_refactor_compact_supernode_batch_rows ||
       checked_stats.row_refactor_last_segment_target_input_rows <
         checked_stats.row_refactor_compact_supernode_batch_rows ||
       checked_stats.row_refactor_segment_target_input_rows <
         checked_stats.row_refactor_compact_supernode_batch_rows ||
       checked_stats.row_refactor_compact_supernode_batch_candidate_count < 1 ||
       checked_stats.row_refactor_compact_supernode_batch_rejected_work_count != 0)) {
    fprintf(stderr,
            "unexpected checked batch-group ragged stats: checked=%d/%d"
            ", batch=%d/%" PRId64 "/%" PRId64 "/%" PRId64 "/%" PRId64
            ", direct=%" PRId64 "/%" PRId64
            ", target=%" PRId64 "/%" PRId64
            ", pattern=%" PRId64 "/%" PRId64
            ", candidate=%" PRId64 "/%" PRId64 "/%" PRId64
            ", rejected=%" PRId64 "\n",
            checked_stats.row_refactor_last_run,
            checked_stats.row_refactor_last_checked,
            checked_stats.row_refactor_last_compact_supernode_batch,
            checked_stats.row_refactor_compact_supernode_batch_count,
            checked_stats.row_refactor_compact_supernode_batch_rows,
            checked_stats.row_refactor_compact_supernode_batch_dep_rows,
            checked_stats.row_refactor_compact_supernode_batch_entries,
            checked_stats.row_refactor_last_batch_direct_input_rows,
            checked_stats.row_refactor_batch_direct_input_rows,
            checked_stats.row_refactor_last_segment_target_input_rows,
            checked_stats.row_refactor_segment_target_input_rows,
            checked_stats.row_refactor_compact_supernode_batch_pattern_count,
            checked_stats.row_refactor_compact_supernode_batch_pattern_rows,
            checked_stats.row_refactor_compact_supernode_batch_candidate_count,
            checked_stats.row_refactor_compact_supernode_batch_candidate_rows,
            checked_stats.row_refactor_compact_supernode_batch_candidate_dep_rows,
            checked_stats.row_refactor_compact_supernode_batch_rejected_work_count);
    ok = 0;
  }

  max_solution_error = 0.0;
  max_residual = 0.0;
  max_rhs = 0.0;
  for (int32_t row = 0; row < n; ++row) {
    const double err = fabs(x[row] - expected[row]);
    if (err > max_solution_error) {
      max_solution_error = err;
    }
    if (residual != NULL) {
      residual[row] = -b[row];
    }
    if (fabs(b[row]) > max_rhs) {
      max_rhs = fabs(b[row]);
    }
  }
  if (residual != NULL) {
    for (int32_t col = 0; col < n; ++col) {
      for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
        residual[ai[p]] += ax1[p] * x[col];
      }
    }
    for (int32_t row = 0; row < n; ++row) {
      const double residual_abs = fabs(residual[row]);
      if (residual_abs > max_residual) {
        max_residual = residual_abs;
      }
    }
  }
  const double checked_relative_residual =
    max_residual / (max_rhs > 0.0 ? max_rhs : 1.0);
  if (ok && (max_solution_error > 1.0e-8 ||
             checked_relative_residual > 1.0e-10)) {
    fprintf(stderr,
            "unexpected checked batch-group ragged accuracy:"
            " max_x_err=%.17g, rel_resid=%.17g\n",
            max_solution_error, checked_relative_residual);
    ok = 0;
  }

  kls_destroy(solver);
  kls_destroy(checked_solver);
  free(saved_row_env);
  free(saved_cblas_env);
  free(saved_checked_env);
  free(saved_native_env);
  free(ap);
  free(ai);
  free(ax0);
  free(ax1);
  free(b);
  free(x);
  free(expected);
  free(residual);
  return ok;
}

static int test_pre_static_pivoting(void) {
  const int32_t n = 3000;
  int32_t *ap = (int32_t *)malloc(((size_t)n + 1u) * sizeof(*ap));
  int32_t *ai = (int32_t *)malloc((size_t)n * sizeof(*ai));
  double *ax = (double *)malloc((size_t)n * sizeof(*ax));
  double *b = (double *)calloc((size_t)n, sizeof(*b));
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
  double *ax_ref = (double *)malloc((size_t)n * sizeof(*ax_ref));
  double *b = (double *)calloc((size_t)n, sizeof(*b));
  double *bt = (double *)calloc((size_t)n, sizeof(*bt));
  double *x = (double *)calloc((size_t)n, sizeof(*x));
  double *xt = (double *)calloc((size_t)n, sizeof(*xt));
  double *expected = (double *)malloc((size_t)n * sizeof(*expected));
  double *expected_t = (double *)malloc((size_t)n * sizeof(*expected_t));
  if (ap == NULL || ai == NULL || ax == NULL || ax_ref == NULL ||
      b == NULL || bt == NULL || x == NULL || xt == NULL ||
      expected == NULL || expected_t == NULL) {
    free(ap);
    free(ai);
    free(ax);
    free(ax_ref);
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
    ax_ref[col] = magnitude * (1.0 + 0.001 * (double)(1 + (col % 7)));
    expected[col] = 1.0 + (double)(col % 19);
    expected_t[row] = 2.0 + (double)(row % 23);
    b[row] = magnitude * expected[col];
    bt[col] = magnitude * expected_t[row];
  }

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_AUTO;

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
  if (ok && !stats.selected_exact_matching_scaling) {
    fprintf(stderr, "scaled pre-static pivoting did not retain exact matching scaling\n");
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
  const int64_t row_solve_runs_before_refactor =
    stats.row_refactor_row_solve_run_count;
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok) {
    memset(b, 0, (size_t)n * sizeof(*b));
    memset(bt, 0, (size_t)n * sizeof(*bt));
    memset(x, 0, (size_t)n * sizeof(*x));
    memset(xt, 0, (size_t)n * sizeof(*xt));
    for (int32_t col = 0; col < n; ++col) {
      const int32_t row = ai[col];
      b[row] = ax_ref[col] * expected[col];
      bt[col] = ax_ref[col] * expected_t[row];
    }
  }
  if (ok && !require_ok(kls_refactor(solver, ax_ref),
                        "row refactor scaled pre-static pivot")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "row solve scaled pre-static pivot")) ok = 0;
  if (ok && !require_ok(kls_solve_transpose(solver, 1, bt, 0, xt, 0),
                        "row transpose solve scaled pre-static pivot")) ok = 0;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "row stats scaled pre-static pivot")) {
    ok = 0;
  }
  if (ok && (stats.row_refactor_last_run != 1 ||
             stats.row_refactor_values_dirty != 1 ||
             stats.row_refactor_last_row_solve != 1 ||
             stats.row_refactor_row_solve_run_count <
               row_solve_runs_before_refactor + 2)) {
    fprintf(stderr,
            "scaled pre-static pivoting did not use dirty row-major solves:"
            " row_run=%d dirty=%d row_solve=%d/%" PRId64
            " before=%" PRId64 "\n",
            stats.row_refactor_last_run, stats.row_refactor_values_dirty,
            stats.row_refactor_last_row_solve,
            stats.row_refactor_row_solve_run_count,
            row_solve_runs_before_refactor);
    ok = 0;
  }
  for (int32_t i = 0; ok && i < n; ++i) {
    if (!close_enough(x[i], expected[i])) {
      fprintf(stderr,
              "unexpected row-refactored scaled pre-static solution at %d:"
              " %.17g != %.17g\n",
              (int)i, x[i], expected[i]);
      ok = 0;
    }
    if (ok && !close_enough(xt[i], expected_t[i])) {
      fprintf(stderr,
              "unexpected row-refactored scaled pre-static transpose at %d:"
              " %.17g != %.17g\n",
              (int)i, xt[i], expected_t[i]);
      ok = 0;
    }
  }

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

  kls_destroy(solver);
  free(saved_env);
  free(saved_checked_env);
  free(ap);
  free(ai);
  free(ax);
  free(ax_ref);
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
  double expected2[12];
  double b[24] = {0.0};
  double x[24] = {0.0};

  int32_t p = 0;
  for (int32_t col = 0; col < n; ++col) {
    ap[col] = p;
    expected[col] = 1.0 + 0.25 * (double)col;
    expected2[col] = -0.5 + 0.125 * (double)col;
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
      b[n + ai[q]] += ax1[q] * expected2[col];
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
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR=1");
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(solver, ax1),
                        "explicit checked row-pipeline factor")) ok = 0;
  kls_stats checked_stats;
  checked_stats.struct_size = sizeof(checked_stats);
  if (ok && !require_ok(kls_get_stats(solver, &checked_stats),
                        "stats explicit checked row-pipeline factor")) {
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
             checked_stats.row_refactor_work_ready_queue_run_count != 1 ||
             checked_stats.row_refactor_last_private_ready_groups <= 0 ||
             checked_stats.row_refactor_private_ready_group_count <
               checked_stats.row_refactor_last_private_ready_groups ||
             checked_stats.row_refactor_ready_queue_workspace_groups <
               checked_stats.row_refactor_group_count)) {
    fprintf(stderr,
            "unexpected checked row stats: last=%d/%d/%d/%d/%d/%d/%d"
            ", runs=%" PRId64 "/%" PRId64 "/%" PRId64
            ", ready=%" PRId64 "/%" PRId64
            ", done_bitmap=%" PRId64
            ", defer_scatter=%" PRId64
            ", work_queue=%" PRId64
            ", private_ready=%" PRId64 "/%" PRId64
            ", queue_workspace=%" PRId64 "\n",
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
            checked_stats.row_refactor_work_ready_queue_run_count,
            checked_stats.row_refactor_last_private_ready_groups,
            checked_stats.row_refactor_private_ready_group_count,
            checked_stats.row_refactor_ready_queue_workspace_groups);
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
             checked_only_refactor_stats.row_refactor_last_private_ready_groups != 0 ||
             checked_only_refactor_stats.row_refactor_ready_queue_workspace_groups !=
               checked_stats.row_refactor_ready_queue_workspace_groups ||
             checked_only_refactor_stats.row_refactor_run_count != 1 ||
             checked_only_refactor_stats.row_refactor_checked_run_count != 1)) {
    fprintf(stderr,
            "checked-only refactor incorrectly used row path: last=%d/%d/%d/%d/%d/%" PRId64
            ", workspace=%" PRId64
            ", runs=%" PRId64 ", checked=%" PRId64 "\n",
            checked_only_refactor_stats.row_refactor_last_run,
            checked_only_refactor_stats.row_refactor_last_ready_queue,
            checked_only_refactor_stats.row_refactor_last_done_bitmap,
            checked_only_refactor_stats.row_refactor_last_defer_value_scatter,
            checked_only_refactor_stats.row_refactor_last_work_ready_queue,
            checked_only_refactor_stats.row_refactor_last_private_ready_groups,
            checked_only_refactor_stats.row_refactor_ready_queue_workspace_groups,
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
  if (ok && !require_ok(kls_solve(solver, 2, b, n, x, n),
                        "solve row-pipeline refactor")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats row-pipeline refactor")) {
    ok = 0;
  }
  const int expect_defer = checked_expect_defer + 1;
  const int used_partial_supernode =
    stats.row_refactor_last_partial_supernode_pipeline != 0;
  const int expected_ready_queue_runs = 2;
  const int expected_done_bitmap_runs = used_partial_supernode ? 2 : 1;
  const int expected_partial_supernode_runs =
    used_partial_supernode ? expected_ready_queue_runs : 0;
  const int expected_work_queue_runs = 2;
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
             stats.row_refactor_last_done_bitmap !=
               (used_partial_supernode ? 1 : 0) ||
             stats.row_refactor_last_defer_value_scatter != 1 ||
             stats.row_refactor_values_dirty != 1 ||
             stats.row_refactor_last_lazy_value_scatter != 1 ||
             stats.row_refactor_last_row_solve != 1 ||
             stats.row_refactor_last_work_ready_queue != 1 ||
             stats.row_refactor_run_count != 2 ||
             stats.row_refactor_checked_run_count != 1 ||
             stats.row_refactor_parallel_run_count != 2 ||
             stats.row_refactor_ready_queue_run_count !=
               expected_ready_queue_runs ||
             stats.row_refactor_ready_queue_group_count <
               stats.row_refactor_group_pipeline_groups ||
             stats.row_refactor_done_bitmap_run_count !=
               expected_done_bitmap_runs ||
             stats.row_refactor_input_cleanup_rows != 0 ||
             stats.row_refactor_input_cleanup_entries != 0 ||
             stats.row_refactor_defer_value_scatter_run_count !=
               expect_defer ||
             stats.row_refactor_lazy_value_scatter_run_count != 1 ||
             stats.row_refactor_row_solve_run_count != 1 ||
             stats.row_refactor_work_ready_queue_run_count !=
               expected_work_queue_runs ||
             (used_partial_supernode &&
              (stats.row_refactor_last_partial_supernode_pipeline_groups <= 0 ||
               stats.row_refactor_last_partial_supernode_pipeline_rows <
                 stats.row_refactor_last_partial_supernode_pipeline_groups ||
               stats.row_refactor_partial_supernode_pipeline_run_count !=
                 expected_partial_supernode_runs)) ||
             (!used_partial_supernode &&
              stats.row_refactor_partial_supernode_pipeline_run_count !=
                expected_partial_supernode_runs) ||
             stats.row_refactor_last_private_ready_groups <= 0 ||
             stats.row_refactor_private_ready_group_count <
               stats.row_refactor_last_private_ready_groups ||
             stats.row_refactor_ready_queue_workspace_groups <
               stats.row_refactor_group_count)) {
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
            ", dirty/lazy/row_solve=%d/%d/%" PRId64 "/%d/%" PRId64
            ", work_queue=%" PRId64
            ", partial_supernode=%d/%" PRId64 "/%" PRId64 "/%" PRId64
            ", private_ready=%" PRId64 "/%" PRId64
            ", queue_workspace=%" PRId64 "\n",
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
            stats.row_refactor_values_dirty,
            stats.row_refactor_last_lazy_value_scatter,
            stats.row_refactor_lazy_value_scatter_run_count,
            stats.row_refactor_last_row_solve,
            stats.row_refactor_row_solve_run_count,
            stats.row_refactor_work_ready_queue_run_count,
            stats.row_refactor_last_partial_supernode_pipeline,
            stats.row_refactor_last_partial_supernode_pipeline_groups,
            stats.row_refactor_last_partial_supernode_pipeline_rows,
            stats.row_refactor_partial_supernode_pipeline_run_count,
            stats.row_refactor_last_private_ready_groups,
            stats.row_refactor_private_ready_group_count,
            stats.row_refactor_ready_queue_workspace_groups);
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
    if (!close_enough(x[n + i], expected2[i])) {
      fprintf(stderr, "unexpected row-pipeline solution rhs2 at %d: %.17g != %.17g\n",
              (int)i, x[n + i], expected2[i]);
      ok = 0;
    }
  }

  kls_destroy(solver);
  free(saved_env);
  free(saved_checked_env);
  return ok;
}

static int test_parallel_row_refactor_full_ready_queue(void) {
  const int32_t n = 8;
  int32_t ap[9] = {0, 2, 4, 6, 8, 10, 12, 14, 15};
  int32_t ai[15] = {
    0, 6,
    1, 6,
    2, 6,
    3, 6,
    4, 6,
    5, 6,
    6, 7,
    7
  };
  double ax0[15] = {
    10.0, 0.11,
    11.0, 0.12,
    12.0, 0.13,
    13.0, 0.14,
    14.0, 0.15,
    15.0, 0.16,
    16.0, 0.17,
    17.0
  };
  double ax1[15];
  double expected[8];
  double b[8] = {0.0};
  double x[8] = {0.0};

  for (int32_t p = 0; p < 15; ++p) {
    ax1[p] = ax0[p] + ((p & 1) ? 0.005 : 0.25);
  }
  for (int32_t i = 0; i < n; ++i) {
    expected[i] = 0.5 + 0.125 * (double)i;
  }
  for (int32_t col = 0; col < n; ++col) {
    for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
      b[ai[p]] += ax1[p] * expected[col];
    }
  }

  const char *saved_env_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_env = saved_env_value != NULL ? strdup(saved_env_value) : NULL;
  const int had_saved_env = saved_env_value != NULL;

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.threads = 2;
  options.orientation = KLS_ORIENTATION_NORMAL;
  options.ordering = KLS_ORDERING_NATURAL;
  options.use_btf = 0;
  options.scale = -1;
  options.static_pivoting = 0;

  int ok = 1;
  if (had_saved_env && saved_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_ROW_REFACTOR\n");
    ok = 0;
  }
  if (!require_ok(kls_create(&solver), "create full-ready row")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, n, ap, ai, 0,
                                        &options),
                        "analyze full-ready row")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor full-ready row base")) ok = 0;
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && !require_ok(kls_refactor(solver, ax1),
                        "refactor full-ready row")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, n, x, n),
                        "solve full-ready row")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats full-ready row")) {
    ok = 0;
  }
  if (ok && (stats.row_refactor_group_level_max_width < 6 ||
             stats.row_refactor_group_cluster_levels != 1 ||
             stats.row_refactor_group_pipeline_groups != 2 ||
             stats.row_refactor_group_pipeline_rows != 2 ||
             stats.row_refactor_last_run != 1 ||
             stats.row_refactor_last_parallel != 1 ||
             stats.row_refactor_last_ready_queue != 1 ||
             stats.row_refactor_ready_queue_run_count != 1 ||
             stats.row_refactor_ready_queue_group_count !=
               stats.row_refactor_group_pipeline_groups ||
             stats.row_refactor_last_work_ready_queue != 1 ||
             stats.row_refactor_work_ready_queue_run_count != 1 ||
             stats.row_refactor_last_private_ready_groups <= 0 ||
             stats.row_refactor_private_ready_group_count <
               stats.row_refactor_last_private_ready_groups ||
             stats.row_refactor_last_local_ready_groups <= 0 ||
             stats.row_refactor_local_ready_group_count <
               stats.row_refactor_last_local_ready_groups)) {
    fprintf(stderr,
            "unexpected full-ready row stats: width=%" PRId64
            ", groups=%" PRId64 ", cluster=%" PRId64
            ", pipe=%" PRId64 "/%" PRId64
            ", last=%d/%d/%d, ready=%" PRId64 "/%" PRId64
            ", work_queue=%d/%" PRId64
            ", private_ready=%" PRId64 "/%" PRId64
            ", local_ready=%" PRId64 "/%" PRId64 "\n",
            stats.row_refactor_group_level_max_width,
            stats.row_refactor_group_count,
            stats.row_refactor_group_cluster_levels,
            stats.row_refactor_group_pipeline_groups,
            stats.row_refactor_group_pipeline_rows,
            stats.row_refactor_last_run,
            stats.row_refactor_last_parallel,
            stats.row_refactor_last_ready_queue,
            stats.row_refactor_ready_queue_run_count,
            stats.row_refactor_ready_queue_group_count,
            stats.row_refactor_last_work_ready_queue,
            stats.row_refactor_work_ready_queue_run_count,
            stats.row_refactor_last_private_ready_groups,
            stats.row_refactor_private_ready_group_count,
            stats.row_refactor_last_local_ready_groups,
            stats.row_refactor_local_ready_group_count);
    ok = 0;
  }
  for (int32_t i = 0; ok && i < n; ++i) {
    if (!close_enough(x[i], expected[i])) {
      fprintf(stderr,
              "unexpected full-ready row solution at %d: %.17g != %.17g\n",
              (int)i, x[i], expected[i]);
      ok = 0;
    }
  }

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
  kls_destroy(solver);
  free(saved_env);
  return ok;
}

static int run_scaled_row_refactor_case(int threads, int expect_parallel) {
  const int32_t n = 12;
  int32_t ap[13];
  int32_t ai[49];
  double ax0[49];
  double ax1[49];
  double expected[12];
  double expected_t[12];
  double b[12] = {0.0};
  double bt[12] = {0.0};
  double x[12] = {0.0};
  double xt[12] = {0.0};

  int32_t p = 0;
  for (int32_t col = 0; col < n; ++col) {
    ap[col] = p;
    expected[col] = 1.0 + 0.25 * (double)col;
    expected_t[col] = 0.75 + 0.125 * (double)col;
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
      bt[col] += ax1[q] * expected_t[ai[q]];
    }
  }

  const char *saved_env_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_env = saved_env_value != NULL ? strdup(saved_env_value) : NULL;
  const int had_saved_env = saved_env_value != NULL;
  const char *saved_checked_env_value =
    getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  char *saved_checked_env = saved_checked_env_value != NULL
    ? strdup(saved_checked_env_value) : NULL;
  const int had_saved_checked_env = saved_checked_env_value != NULL;

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.threads = threads;
  options.orientation = KLS_ORIENTATION_NORMAL;
  options.ordering = KLS_ORDERING_NATURAL;
  options.use_btf = 0;
  options.scale = 2;
  options.static_pivoting = 0;

  int ok = 1;
  if (had_saved_env && saved_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_ROW_REFACTOR\n");
    ok = 0;
  }
  if (had_saved_checked_env && saved_checked_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_CHECKED_ROW_REFACTOR\n");
    ok = 0;
  }
  if (!require_ok(kls_create(&solver), "create scaled row refactor")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, n, ap, ai, 0,
                                        &options),
                        "analyze scaled row refactor")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor scaled row refactor base")) ok = 0;
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && !require_ok(kls_refactor(solver, ax1),
                        "scaled row refactor")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, n, x, n),
                        "solve scaled row refactor")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats scaled row refactor")) ok = 0;
  if (ok && (stats.selected_scale != 2 ||
             stats.row_refactor_last_run != 1 ||
             stats.row_refactor_last_checked != 0 ||
             stats.row_refactor_last_parallel != expect_parallel ||
             stats.row_refactor_last_defer_value_scatter != 1 ||
             stats.row_refactor_values_dirty != 1 ||
             stats.row_refactor_last_lazy_value_scatter != 1 ||
             stats.row_refactor_last_row_solve != 1 ||
             stats.row_refactor_run_count != 1 ||
             stats.row_refactor_lazy_value_scatter_run_count != 1 ||
             stats.row_refactor_row_solve_run_count != 1)) {
    fprintf(stderr,
            "unexpected scaled row-refactor stats for %d threads:"
            " scale=%d, last=%d/%d/%d, defer=%d, dirty=%d"
            ", lazy=%d/%" PRId64 ", row_solve=%d/%" PRId64
            ", runs=%" PRId64 "\n",
            threads, stats.selected_scale, stats.row_refactor_last_run,
            stats.row_refactor_last_checked,
            stats.row_refactor_last_parallel,
            stats.row_refactor_last_defer_value_scatter,
            stats.row_refactor_values_dirty,
            stats.row_refactor_last_lazy_value_scatter,
            stats.row_refactor_lazy_value_scatter_run_count,
            stats.row_refactor_last_row_solve,
            stats.row_refactor_row_solve_run_count,
            stats.row_refactor_run_count);
    ok = 0;
  }
  for (int32_t i = 0; ok && i < n; ++i) {
    if (!close_enough(x[i], expected[i])) {
      fprintf(stderr,
              "unexpected scaled row-refactor solution at %d: %.17g != %.17g\n",
              (int)i, x[i], expected[i]);
      ok = 0;
    }
  }
  if (ok && !require_ok(kls_solve_transpose(solver, 1, bt, n, xt, n),
                        "transpose solve scaled row refactor")) ok = 0;
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "transpose stats scaled row refactor")) ok = 0;
  if (ok && (stats.row_refactor_values_dirty != 1 ||
             stats.row_refactor_last_lazy_value_scatter != 1 ||
             stats.row_refactor_last_row_solve != 1 ||
             stats.row_refactor_lazy_value_scatter_run_count != 1 ||
             stats.row_refactor_row_solve_run_count != 2)) {
    fprintf(stderr,
            "unexpected scaled row-refactor transpose stats for %d threads:"
            " dirty=%d, lazy=%d/%" PRId64 ", row_solve=%d/%" PRId64 "\n",
            threads, stats.row_refactor_values_dirty,
            stats.row_refactor_last_lazy_value_scatter,
            stats.row_refactor_lazy_value_scatter_run_count,
            stats.row_refactor_last_row_solve,
            stats.row_refactor_row_solve_run_count);
    ok = 0;
  }
  for (int32_t i = 0; ok && i < n; ++i) {
    if (!close_enough(xt[i], expected_t[i])) {
      fprintf(stderr,
              "unexpected scaled row-refactor transpose solution at %d:"
              " %.17g != %.17g\n",
              (int)i, xt[i], expected_t[i]);
      ok = 0;
    }
  }

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
  kls_destroy(solver);
  free(saved_env);
  free(saved_checked_env);
  return ok;
}

static int test_checked_separator_flop_ready_queue(void) {
  const int32_t nx = 180;
  const int32_t ny = 170;
  const int32_t n = nx * ny;
  int64_t nnz64 = 0;
  for (int32_t y = 0; y < ny; ++y) {
    for (int32_t x = 0; x < nx; ++x) {
      nnz64 += 1;
      if (x > 0) nnz64++;
      if (x + 1 < nx) nnz64++;
      if (y > 0) nnz64++;
      if (y + 1 < ny) nnz64++;
    }
  }
  if (nnz64 <= 0 || nnz64 > INT32_MAX) {
    return 0;
  }
  const int32_t nnz = (int32_t)nnz64;
  int32_t *ap = (int32_t *)malloc(((size_t)n + 1u) * sizeof(*ap));
  int32_t *ai = (int32_t *)malloc((size_t)nnz * sizeof(*ai));
  double *ax0 = (double *)malloc((size_t)nnz * sizeof(*ax0));
  double *ax1 = (double *)malloc((size_t)nnz * sizeof(*ax1));
  double *b = (double *)calloc((size_t)n, sizeof(*b));
  double *xvec = (double *)calloc((size_t)n, sizeof(*xvec));
  double *expected = (double *)malloc((size_t)n * sizeof(*expected));
  double *residual = NULL;
  if (ap == NULL || ai == NULL || ax0 == NULL || ax1 == NULL ||
      b == NULL || xvec == NULL || expected == NULL) {
    free(ap);
    free(ai);
    free(ax0);
    free(ax1);
    free(b);
    free(xvec);
    free(expected);
    return 0;
  }

  int32_t p = 0;
  for (int32_t y = 0; y < ny; ++y) {
    for (int32_t x = 0; x < nx; ++x) {
      const int32_t col = y * nx + x;
      ap[col] = p;
      expected[col] = 1.0 + 0.001 * (double)(col % 17);
      if (y > 0) {
        ai[p] = col - nx;
        ax0[p] = -1.0;
        ax1[p] = -1.0 - 1.0e-5 * (double)((col % 5) + 1);
        p++;
      }
      if (x > 0) {
        ai[p] = col - 1;
        ax0[p] = -1.0;
        ax1[p] = -1.0 + 1.0e-5 * (double)((col % 7) + 1);
        p++;
      }
      ai[p] = col;
      ax0[p] = 5.0 + 1.0e-4 * (double)(col % 13);
      ax1[p] = ax0[p] + 1.0e-3 * (double)((col % 3) + 1);
      p++;
      if (x + 1 < nx) {
        ai[p] = col + 1;
        ax0[p] = -1.0;
        ax1[p] = -1.0 - 1.0e-5 * (double)((col % 11) + 1);
        p++;
      }
      if (y + 1 < ny) {
        ai[p] = col + nx;
        ax0[p] = -1.0;
        ax1[p] = -1.0 + 1.0e-5 * (double)((col % 13) + 1);
        p++;
      }
    }
  }
  ap[n] = p;
  if (p != nnz) {
    fprintf(stderr, "unexpected checked separator fixture nnz: %d/%d\n",
            p, nnz);
    free(ap);
    free(ai);
    free(ax0);
    free(ax1);
    free(b);
    free(xvec);
    free(expected);
    return 0;
  }
  for (int32_t col = 0; col < n; ++col) {
    for (int32_t q = ap[col]; q < ap[col + 1]; ++q) {
      b[ai[q]] += ax1[q] * expected[col];
    }
  }

  const char *saved_row_env_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_row_env =
    saved_row_env_value != NULL ? strdup(saved_row_env_value) : NULL;
  const int had_saved_row_env = saved_row_env_value != NULL;
  const char *saved_checked_env_value =
    getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  char *saved_checked_env = saved_checked_env_value != NULL
    ? strdup(saved_checked_env_value) : NULL;
  const int had_saved_checked_env = saved_checked_env_value != NULL;

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.threads = 4;
  options.orientation = KLS_ORIENTATION_NORMAL;
  options.ordering = KLS_ORDERING_METIS;
  options.use_btf = 0;
  options.scale = -1;
  options.static_pivoting = 0;

  int ok = 1;
  if (had_saved_row_env && saved_row_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_ROW_REFACTOR\n");
    ok = 0;
  }
  if (had_saved_checked_env && saved_checked_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_CHECKED_ROW_REFACTOR\n");
    ok = 0;
  }
  if (!require_ok(kls_create(&solver), "create checked separator queue")) {
    ok = 0;
  }
  if (ok) {
    const int analyze_status =
      kls_analyze_csc(solver, KLS_INDEX_INT32, n, ap, ai, 0, &options);
    if (analyze_status == KLS_ERR_UNSUPPORTED) {
      goto cleanup;
    }
    if (!require_ok(analyze_status, "analyze checked separator queue")) {
      ok = 0;
    }
  }
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor checked separator queue base")) {
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(solver, ax1),
                        "checked separator queue fast factor")) {
    ok = 0;
  }

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats checked separator queue")) {
    ok = 0;
  }
  if (ok && (stats.build_has_metis != 1 ||
             stats.separator_analyzed_rows != n ||
             stats.separator_component_count < 3 ||
             stats.row_refactor_group_count <= 0 ||
             stats.row_refactor_last_run != 1 ||
             stats.row_refactor_last_checked != 1 ||
             stats.row_refactor_last_parallel != 1 ||
             stats.row_refactor_last_ready_queue != 1 ||
             stats.row_refactor_last_done_bitmap != 1 ||
             stats.row_refactor_last_work_ready_queue != 1 ||
             stats.row_refactor_last_separator_flop_queue != 1 ||
             stats.row_refactor_separator_flop_queue_run_count < 1 ||
             stats.row_refactor_last_separator_flop_components < 3 ||
             stats.row_refactor_last_separator_flop_private_groups <= 0 ||
             stats.row_refactor_last_separator_flop_pipeline_groups <= 0 ||
             stats.row_refactor_checked_run_count < 1)) {
    fprintf(stderr,
            "unexpected checked separator queue stats: metis=%d"
            ", sep_rows=%" PRId64 ", sep_components=%" PRId64
            ", groups=%" PRId64 ", last=%d/%d/%d/%d/%d/%d/%d"
            ", sep_queue=%" PRId64 "/%" PRId64 "/%" PRId64 "/%" PRId64
            ", checked=%" PRId64 "\n",
            stats.build_has_metis,
            stats.separator_analyzed_rows,
            stats.separator_component_count,
            stats.row_refactor_group_count,
            stats.row_refactor_last_run,
            stats.row_refactor_last_checked,
            stats.row_refactor_last_parallel,
            stats.row_refactor_last_ready_queue,
            stats.row_refactor_last_done_bitmap,
            stats.row_refactor_last_work_ready_queue,
            stats.row_refactor_last_separator_flop_queue,
            stats.row_refactor_separator_flop_queue_run_count,
            stats.row_refactor_last_separator_flop_components,
            stats.row_refactor_last_separator_flop_private_groups,
            stats.row_refactor_last_separator_flop_pipeline_groups,
            stats.row_refactor_checked_run_count);
    ok = 0;
  }

  if (ok && !require_ok(kls_solve(solver, 1, b, 0, xvec, 0),
                        "solve checked separator queue")) {
    ok = 0;
  }
  double max_solution_error = 0.0;
  double max_residual = 0.0;
  double max_rhs = 0.0;
  for (int32_t row = 0; row < n; ++row) {
    const double err = fabs(xvec[row] - expected[row]);
    if (err > max_solution_error) {
      max_solution_error = err;
    }
    if (fabs(b[row]) > max_rhs) {
      max_rhs = fabs(b[row]);
    }
  }
  residual = (double *)malloc((size_t)n * sizeof(*residual));
  if (residual == NULL) {
    ok = 0;
  } else {
    for (int32_t row = 0; row < n; ++row) {
      residual[row] = -b[row];
    }
    for (int32_t col = 0; col < n; ++col) {
      for (int32_t q = ap[col]; q < ap[col + 1]; ++q) {
        residual[ai[q]] += ax1[q] * xvec[col];
      }
    }
    for (int32_t row = 0; row < n; ++row) {
      if (fabs(residual[row]) > max_residual) {
        max_residual = fabs(residual[row]);
      }
    }
  }
  const double relative_residual =
    max_residual / (max_rhs > 0.0 ? max_rhs : 1.0);
  if (ok && (max_solution_error > 1.0e-7 ||
             relative_residual > 1.0e-8)) {
    fprintf(stderr,
            "unexpected checked separator queue accuracy:"
            " max_x_err=%.17g, rel_resid=%.17g\n",
            max_solution_error, relative_residual);
    ok = 0;
  }

cleanup:
  if (had_saved_row_env && saved_row_env != NULL) {
    if (setenv("KLS_ENABLE_ROW_REFACTOR", saved_row_env, 1) != 0) {
      perror("restore KLS_ENABLE_ROW_REFACTOR");
      ok = 0;
    }
  } else if (!had_saved_row_env) {
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
  kls_destroy(solver);
  free(saved_row_env);
  free(saved_checked_env);
  free(ap);
  free(ai);
  free(ax0);
  free(ax1);
  free(b);
  free(xvec);
  free(expected);
  free(residual);
  return ok;
}

static int test_btf_duplicate_separator_forest(void) {
  const int32_t nx = 180;
  const int32_t ny = 170;
  const int32_t block_count = 2;
  const int32_t block_n = nx * ny;
  const int32_t n = block_count * block_n;
  int64_t nnz64 = 0;
  for (int32_t block = 0; block < block_count; ++block) {
    (void)block;
    for (int32_t y = 0; y < ny; ++y) {
      for (int32_t x = 0; x < nx; ++x) {
        nnz64 += 1;
        if (x > 0) nnz64++;
        if (x + 1 < nx) nnz64++;
        if (y > 0) nnz64++;
        if (y + 1 < ny) nnz64++;
      }
    }
  }
  if (nnz64 <= 0 || nnz64 > INT32_MAX) {
    return 0;
  }
  const int32_t nnz = (int32_t)nnz64;
  int32_t *ap = (int32_t *)malloc(((size_t)n + 1u) * sizeof(*ap));
  int32_t *ai = (int32_t *)malloc((size_t)nnz * sizeof(*ai));
  double *ax0 = (double *)malloc((size_t)nnz * sizeof(*ax0));
  double *ax1 = (double *)malloc((size_t)nnz * sizeof(*ax1));
  double *b = (double *)calloc((size_t)n, sizeof(*b));
  double *xvec = (double *)calloc((size_t)n, sizeof(*xvec));
  double *expected = (double *)malloc((size_t)n * sizeof(*expected));
  if (ap == NULL || ai == NULL || ax0 == NULL || ax1 == NULL ||
      b == NULL || xvec == NULL || expected == NULL) {
    free(ap);
    free(ai);
    free(ax0);
    free(ax1);
    free(b);
    free(xvec);
    free(expected);
    return 0;
  }

  int32_t p = 0;
  for (int32_t block = 0; block < block_count; ++block) {
    const int32_t base = block * block_n;
    for (int32_t y = 0; y < ny; ++y) {
      for (int32_t x = 0; x < nx; ++x) {
        const int32_t local = y * nx + x;
        const int32_t col = base + local;
        ap[col] = p;
        expected[col] = 1.0 + 0.001 * (double)((col + 3 * block) % 19);
        if (y > 0) {
          ai[p] = col - nx;
          ax0[p] = -1.0;
          ax1[p] = -1.0 - 1.0e-5 * (double)((local % 5) + 1);
          p++;
        }
        if (x > 0) {
          ai[p] = col - 1;
          ax0[p] = -1.0;
          ax1[p] = -1.0 + 1.0e-5 * (double)((local % 7) + 1);
          p++;
        }
        ai[p] = col;
        ax0[p] = 5.0 + 1.0e-4 * (double)((local + block) % 13);
        ax1[p] = ax0[p] + 1.0e-3 * (double)((local % 3) + 1);
        p++;
        if (x + 1 < nx) {
          ai[p] = col + 1;
          ax0[p] = -1.0;
          ax1[p] = -1.0 - 1.0e-5 * (double)((local % 11) + 1);
          p++;
        }
        if (y + 1 < ny) {
          ai[p] = col + nx;
          ax0[p] = -1.0;
          ax1[p] = -1.0 + 1.0e-5 * (double)((local % 13) + 1);
          p++;
        }
      }
    }
  }
  ap[n] = p;
  if (p != nnz) {
    fprintf(stderr, "unexpected duplicate separator fixture nnz: %d/%d\n",
            p, nnz);
    free(ap);
    free(ai);
    free(ax0);
    free(ax1);
    free(b);
    free(xvec);
    free(expected);
    return 0;
  }
  for (int32_t col = 0; col < n; ++col) {
    for (int32_t q = ap[col]; q < ap[col + 1]; ++q) {
      b[ai[q]] += ax1[q] * expected[col];
    }
  }

  const char *saved_row_env_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_row_env =
    saved_row_env_value != NULL ? strdup(saved_row_env_value) : NULL;
  const int had_saved_row_env = saved_row_env_value != NULL;
  const char *saved_checked_env_value =
    getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  char *saved_checked_env = saved_checked_env_value != NULL
    ? strdup(saved_checked_env_value) : NULL;
  const int had_saved_checked_env = saved_checked_env_value != NULL;

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.threads = 4;
  options.orientation = KLS_ORIENTATION_NORMAL;
  options.ordering = KLS_ORDERING_METIS;
  options.use_btf = 1;
  options.scale = -1;
  options.static_pivoting = 0;

  int ok = 1;
  if (had_saved_row_env && saved_row_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_ROW_REFACTOR\n");
    ok = 0;
  }
  if (had_saved_checked_env && saved_checked_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_CHECKED_ROW_REFACTOR\n");
    ok = 0;
  }
  if (!require_ok(kls_create(&solver), "create duplicate separator forest")) {
    ok = 0;
  }
  if (ok) {
    const int analyze_status =
      kls_analyze_csc(solver, KLS_INDEX_INT32, n, ap, ai, 0, &options);
    if (analyze_status == KLS_ERR_UNSUPPORTED) {
      goto cleanup;
    }
    if (!require_ok(analyze_status, "analyze duplicate separator forest")) {
      ok = 0;
    }
  }
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor duplicate separator forest base")) {
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(solver, ax1),
                        "checked duplicate separator forest fast factor")) {
    ok = 0;
  }

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats duplicate separator forest")) {
    ok = 0;
  }
  if (ok && (stats.build_has_metis != 1 ||
             stats.nblocks < 2 ||
             stats.max_block != block_n ||
             stats.separator_analyzed_rows != n ||
             stats.separator_global_begin != 0 ||
             stats.separator_global_end != n ||
             stats.separator_component_count < 15 ||
             stats.separator_private_components < 8 ||
             stats.separator_pipeline_components < 7 ||
             stats.row_refactor_last_run != 1 ||
             stats.row_refactor_last_checked != 1 ||
             stats.row_refactor_last_parallel != 1 ||
             stats.row_refactor_last_ready_queue != 1 ||
             stats.row_refactor_last_done_bitmap != 1 ||
             stats.row_refactor_last_work_ready_queue != 1 ||
             stats.row_refactor_last_separator_flop_queue != 1 ||
             stats.row_refactor_separator_flop_queue_run_count < 1 ||
             stats.row_refactor_last_separator_flop_components < 7 ||
             stats.row_refactor_last_separator_flop_private_groups <= 0 ||
             stats.row_refactor_last_separator_flop_pipeline_groups <= 0)) {
    fprintf(stderr,
            "unexpected duplicate separator forest stats: metis=%d"
            ", blocks=%" PRId64 ", max=%" PRId64
            ", sep_rows=%" PRId64 ", range=[%" PRId64 ",%" PRId64 ")"
            ", sep_components=%" PRId64 "/%" PRId64 "/%" PRId64
            ", last=%d/%d/%d/%d/%d/%d/%d"
            ", sep_queue=%" PRId64 "/%" PRId64 "/%" PRId64 "/%" PRId64 "\n",
            stats.build_has_metis,
            stats.nblocks,
            stats.max_block,
            stats.separator_analyzed_rows,
            stats.separator_global_begin,
            stats.separator_global_end,
            stats.separator_component_count,
            stats.separator_private_components,
            stats.separator_pipeline_components,
            stats.row_refactor_last_run,
            stats.row_refactor_last_checked,
            stats.row_refactor_last_parallel,
            stats.row_refactor_last_ready_queue,
            stats.row_refactor_last_done_bitmap,
            stats.row_refactor_last_work_ready_queue,
            stats.row_refactor_last_separator_flop_queue,
            stats.row_refactor_separator_flop_queue_run_count,
            stats.row_refactor_last_separator_flop_components,
            stats.row_refactor_last_separator_flop_private_groups,
            stats.row_refactor_last_separator_flop_pipeline_groups);
    ok = 0;
  }

  if (ok && !require_ok(kls_solve(solver, 1, b, 0, xvec, 0),
                        "solve duplicate separator forest")) {
    ok = 0;
  }
  double max_solution_error = 0.0;
  for (int32_t row = 0; row < n; ++row) {
    const double err = fabs(xvec[row] - expected[row]);
    if (err > max_solution_error) {
      max_solution_error = err;
    }
  }
  if (ok && max_solution_error > 1.0e-7) {
    fprintf(stderr,
            "unexpected duplicate separator forest solution error: %.17g\n",
            max_solution_error);
    ok = 0;
  }

cleanup:
  if (had_saved_row_env && saved_row_env != NULL) {
    if (setenv("KLS_ENABLE_ROW_REFACTOR", saved_row_env, 1) != 0) {
      perror("restore KLS_ENABLE_ROW_REFACTOR");
      ok = 0;
    }
  } else if (!had_saved_row_env) {
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
  kls_destroy(solver);
  free(saved_row_env);
  free(saved_checked_env);
  free(ap);
  free(ai);
  free(ax0);
  free(ax1);
  free(b);
  free(xvec);
  free(expected);
  return ok;
}

static int test_scaled_row_refactor_single_block(void) {
  return run_scaled_row_refactor_case(1, 0) &&
         run_scaled_row_refactor_case(4, 1);
}

static int run_experimental_kls_first_factor_case(int scale,
                                                 kls_orientation orientation) {
  const int32_t ap[] = {0, 2, 4, 5};
  const int32_t ai[] = {0, 1, 0, 1, 2};
  const double ax[] = {4.0, 2.0, 1.0, 3.0, 5.0};
  const double b[] = {6.0, 8.0, 15.0};
  const double bt[] = {8.0, 7.0, 15.0};
  const double ax_ref[] = {5.0, 1.0, 2.0, 4.0, 6.0};
  const double b_ref[] = {9.0, 9.0, 18.0};
  const double bt_ref[] = {7.0, 10.0, 18.0};
  const double ax_fast[] = {6.0, 1.0, 2.0, 5.0, 7.0};
  const double b_fast[] = {10.0, 11.0, 21.0};
  const double bt_fast[] = {8.0, 12.0, 21.0};
  double x[3] = {0.0, 0.0, 0.0};
  double xt[3] = {0.0, 0.0, 0.0};

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_NATURAL;
  options.orientation = orientation;
  options.use_btf = 1;
  options.scale = scale;
  options.static_pivoting = 0;

  const char *saved_env_value = getenv("KLS_ENABLE_KLS_FIRST_FACTOR");
  char *saved_env = saved_env_value != NULL ? strdup(saved_env_value) : NULL;
  const int had_saved_env = saved_env_value != NULL;
  const char *saved_row_env_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_row_env =
    saved_row_env_value != NULL ? strdup(saved_row_env_value) : NULL;
  const int had_saved_row_env = saved_row_env_value != NULL;
  const char *saved_checked_env_value =
    getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  char *saved_checked_env =
    saved_checked_env_value != NULL ? strdup(saved_checked_env_value) : NULL;
  const int had_saved_checked_env = saved_checked_env_value != NULL;

  int ok = 1;
  if (had_saved_env && saved_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_KLS_FIRST_FACTOR\n");
    ok = 0;
  }
  if (had_saved_row_env && saved_row_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_ROW_REFACTOR\n");
    ok = 0;
  }
  if (had_saved_checked_env && saved_checked_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_CHECKED_ROW_REFACTOR\n");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_KLS_FIRST_FACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_KLS_FIRST_FACTOR");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && !require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 3, ap, ai, 0,
                                        &options),
                        "analyze KLS first factor")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax),
                        "factor KLS first factor")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve KLS first factor")) ok = 0;
  if (ok && !require_ok(kls_solve_transpose(solver, 1, bt, 0, xt, 0),
                        "transpose solve KLS first factor")) ok = 0;

  kls_stats stats;
  memset(&stats, 0, sizeof(stats));
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats KLS first factor")) {
    ok = 0;
  }
  const int kls_first_row_up =
    stats.kls_first_last_row_uplooking_columns == 3 &&
    stats.kls_first_row_uplooking_column_count >= 3;
  const int kls_first_mapped_tail =
    stats.kls_tail_last_mapped_columns > 0 &&
    stats.kls_tail_mapped_column_count >=
      stats.kls_tail_last_mapped_columns;
  const int kls_first_path_ok =
    scale > 0 ? kls_first_row_up
              : (kls_first_row_up || kls_first_mapped_tail);
  const int kls_first_seed_ok =
    !kls_first_row_up ||
    (stats.kls_first_last_row_refactor_seeded_rows == 3 &&
     stats.kls_first_row_refactor_seeded_row_count >= 3);
  if (ok && (stats.last_factor_path != KLS_FACTOR_PATH_KLS_FIRST ||
             stats.nblocks < 2 || stats.max_block < 2 ||
             stats.selected_orientation != orientation ||
             stats.selected_scale != scale ||
             stats.row_refactor_values_dirty != 0 ||
             stats.row_refactor_last_lazy_value_scatter != 0 ||
             stats.row_refactor_last_row_solve != 1 ||
             stats.row_refactor_row_solve_run_count != 2 ||
             stats.fast_block_restarts != 0 ||
             stats.fast_kls_block_restarts != 0 ||
             stats.fast_kls_rebuild_restarts != 0 ||
             !kls_first_path_ok ||
             !kls_first_seed_ok)) {
    fprintf(stderr,
            "unexpected KLS first-factor stats: path=%s, nblocks=%" PRId64
            ", max_block=%" PRId64 ", orientation=%s, scale=%d"
            ", row_dirty=%d, row_lazy=%d, row_solve=%d, row_solve_count=%" PRId64
            ", block_restarts=%d, kls_block_restarts=%d"
            ", kls_rebuild_restarts=%d"
            ", row_up=%" PRId64 "/%" PRId64
            ", row_seed=%" PRId64 "/%" PRId64
            ", mapped_tail=%" PRId64 "/%" PRId64 "\n",
            kls_factor_path_name(stats.last_factor_path),
            stats.nblocks, stats.max_block,
            kls_orientation_name(stats.selected_orientation),
            stats.selected_scale,
            stats.row_refactor_values_dirty,
            stats.row_refactor_last_lazy_value_scatter,
            stats.row_refactor_last_row_solve,
            stats.row_refactor_row_solve_run_count,
            stats.fast_block_restarts, stats.fast_kls_block_restarts,
            stats.fast_kls_rebuild_restarts,
            stats.kls_first_last_row_uplooking_columns,
            stats.kls_first_row_uplooking_column_count,
            stats.kls_first_last_row_refactor_seeded_rows,
            stats.kls_first_row_refactor_seeded_row_count,
            stats.kls_tail_last_mapped_columns,
            stats.kls_tail_mapped_column_count);
    ok = 0;
  }
  if (ok && (!close_enough(x[0], 1.0) || !close_enough(x[1], 2.0) ||
             !close_enough(x[2], 3.0))) {
    fprintf(stderr,
            "unexpected KLS first-factor solution: %.17g %.17g %.17g\n",
            x[0], x[1], x[2]);
    ok = 0;
  }
  if (ok && (!close_enough(xt[0], 1.0) || !close_enough(xt[1], 2.0) ||
             !close_enough(xt[2], 3.0))) {
    fprintf(stderr,
            "unexpected KLS first-factor transpose solution: %.17g %.17g %.17g\n",
            xt[0], xt[1], xt[2]);
    ok = 0;
  }
  if (ok && unsetenv("KLS_ENABLE_ROW_REFACTOR") != 0) {
    perror("unsetenv KLS_ENABLE_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && unsetenv("KLS_ENABLE_CHECKED_ROW_REFACTOR") != 0) {
    perror("unsetenv KLS_ENABLE_CHECKED_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && !require_ok(kls_refactor(solver, ax_ref),
                        "auto row refactor after KLS first factor")) ok = 0;
  x[0] = x[1] = x[2] = 0.0;
  xt[0] = xt[1] = xt[2] = 0.0;
  if (ok && !require_ok(kls_solve(solver, 1, b_ref, 0, x, 0),
                        "solve auto row refactor after KLS first factor")) ok = 0;
  if (ok && !require_ok(kls_solve_transpose(solver, 1, bt_ref, 0, xt, 0),
                        "transpose solve auto row refactor after KLS first factor")) ok = 0;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "auto row refactor stats after KLS first factor")) {
    ok = 0;
  }
  if (ok && (stats.row_refactor_last_run != 1 ||
             stats.row_refactor_last_checked != 0 ||
             stats.row_refactor_values_dirty != 1 ||
             stats.row_refactor_last_lazy_value_scatter != 1 ||
             stats.row_refactor_run_count != 1 ||
             stats.row_refactor_lazy_value_scatter_run_count != 1 ||
             stats.row_refactor_last_row_solve != 1 ||
             stats.row_refactor_row_solve_run_count != 4)) {
    fprintf(stderr,
            "unexpected KLS-first auto row refactor stats: run=%d checked=%d"
            ", dirty=%d, lazy=%d/%" PRId64 ", row_solve=%d/%" PRId64
            ", run_count=%" PRId64 "\n",
            stats.row_refactor_last_run,
            stats.row_refactor_last_checked,
            stats.row_refactor_values_dirty,
            stats.row_refactor_last_lazy_value_scatter,
            stats.row_refactor_lazy_value_scatter_run_count,
            stats.row_refactor_last_row_solve,
            stats.row_refactor_row_solve_run_count,
            stats.row_refactor_run_count);
    ok = 0;
  }
  if (ok && (!close_enough(x[0], 1.0) || !close_enough(x[1], 2.0) ||
             !close_enough(x[2], 3.0))) {
    fprintf(stderr,
            "unexpected KLS-first auto row refactor solution: %.17g %.17g %.17g\n",
            x[0], x[1], x[2]);
    ok = 0;
  }
  if (ok && (!close_enough(xt[0], 1.0) || !close_enough(xt[1], 2.0) ||
             !close_enough(xt[2], 3.0))) {
    fprintf(stderr,
            "unexpected KLS-first auto row refactor transpose solution:"
            " %.17g %.17g %.17g\n",
            xt[0], xt[1], xt[2]);
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR=1");
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(solver, ax_fast),
                        "explicit checked row fast factor after KLS first factor")) {
    ok = 0;
  }
  x[0] = x[1] = x[2] = 0.0;
  xt[0] = xt[1] = xt[2] = 0.0;
  if (ok && !require_ok(kls_solve(solver, 1, b_fast, 0, x, 0),
                        "solve auto checked row fast factor after KLS first factor")) {
    ok = 0;
  }
  if (ok && !require_ok(kls_solve_transpose(solver, 1, bt_fast, 0, xt, 0),
                        "transpose solve auto checked row fast factor after KLS first factor")) {
    ok = 0;
  }
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "explicit checked row fast factor stats after KLS first factor")) {
    ok = 0;
  }
  if (ok && (stats.last_factor_path != KLS_FACTOR_PATH_KLS_FAST_REFACTOR ||
             stats.row_refactor_last_run != 1 ||
             stats.row_refactor_last_checked != 1 ||
             stats.row_refactor_values_dirty != 0 ||
             stats.row_refactor_run_count != 2 ||
             stats.row_refactor_checked_run_count != 1 ||
             stats.row_refactor_last_row_solve != 1 ||
             stats.row_refactor_row_solve_run_count != 6)) {
    fprintf(stderr,
            "unexpected KLS-first auto checked row fast-factor stats:"
            " path=%s, run=%d checked=%d, dirty=%d, run_count=%" PRId64
            ", checked_count=%" PRId64 ", row_solve=%d/%" PRId64 "\n",
            kls_factor_path_name(stats.last_factor_path),
            stats.row_refactor_last_run,
            stats.row_refactor_last_checked,
            stats.row_refactor_values_dirty,
            stats.row_refactor_run_count,
            stats.row_refactor_checked_run_count,
            stats.row_refactor_last_row_solve,
            stats.row_refactor_row_solve_run_count);
    ok = 0;
  }
  if (ok && (!close_enough(x[0], 1.0) || !close_enough(x[1], 2.0) ||
             !close_enough(x[2], 3.0))) {
    fprintf(stderr,
            "unexpected KLS-first auto checked row fast-factor solution:"
            " %.17g %.17g %.17g\n",
            x[0], x[1], x[2]);
    ok = 0;
  }
  if (ok && (!close_enough(xt[0], 1.0) || !close_enough(xt[1], 2.0) ||
             !close_enough(xt[2], 3.0))) {
    fprintf(stderr,
            "unexpected KLS-first auto checked row fast-factor transpose:"
            " %.17g %.17g %.17g\n",
            xt[0], xt[1], xt[2]);
    ok = 0;
  }

  if (had_saved_env && saved_env != NULL) {
    if (setenv("KLS_ENABLE_KLS_FIRST_FACTOR", saved_env, 1) != 0) {
      perror("restore KLS_ENABLE_KLS_FIRST_FACTOR");
      ok = 0;
    }
  } else if (!had_saved_env) {
    if (unsetenv("KLS_ENABLE_KLS_FIRST_FACTOR") != 0) {
      perror("unsetenv KLS_ENABLE_KLS_FIRST_FACTOR");
      ok = 0;
    }
  }
  if (had_saved_row_env && saved_row_env != NULL) {
    if (setenv("KLS_ENABLE_ROW_REFACTOR", saved_row_env, 1) != 0) {
      perror("restore KLS_ENABLE_ROW_REFACTOR");
      ok = 0;
    }
  } else if (!had_saved_row_env) {
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
  kls_destroy(solver);
  free(saved_env);
  free(saved_row_env);
  free(saved_checked_env);
  return ok;
}

static int test_experimental_kls_first_factor(void) {
  return run_experimental_kls_first_factor_case(-1, KLS_ORIENTATION_NORMAL) &&
         run_experimental_kls_first_factor_case(2, KLS_ORIENTATION_NORMAL) &&
         run_experimental_kls_first_factor_case(-1,
                                                KLS_ORIENTATION_TRANSPOSE) &&
         run_experimental_kls_first_factor_case(2,
                                                KLS_ORIENTATION_TRANSPOSE);
}

static int restore_env_value(const char *name, int had_value,
                             const char *saved_value) {
  if (had_value) {
    if (setenv(name, saved_value, 1) != 0) {
      perror("restore env");
      return 0;
    }
  } else if (unsetenv(name) != 0) {
    perror("unsetenv");
    return 0;
  }
  return 1;
}

static int test_checked_row_prefactor_finished_dependency(void) {
  enum {
    A_WIDTH = 360,
    B_WIDTH = 48,
    C_WIDTH = 6,
    B_BEGIN = A_WIDTH,
    C_BEGIN = B_BEGIN + B_WIDTH,
    CONSUMER = C_BEGIN + C_WIDTH,
    PREF_N = CONSUMER + 1,
    PREF_NNZ = A_WIDTH * A_WIDTH + 2 + B_WIDTH * B_WIDTH + 2 +
               C_WIDTH * C_WIDTH + 1 + 1
  };
  const int32_t n = PREF_N;
  const int32_t consumer = CONSUMER;
  const size_t nnz = PREF_NNZ;
  int32_t ap[PREF_N + 1];
  int32_t ai[PREF_NNZ];
  double ax0[PREF_NNZ];
  double ax1[PREF_NNZ];
  double expected[PREF_N];
  double b[PREF_N] = {0.0};
  double x[PREF_N] = {0.0};

  for (int32_t i = 0; i < n; ++i) {
    expected[i] = 0.375 + 0.0625 * (double)((5 * i) % 11);
  }

  size_t pos = 0;
  for (int32_t col = 0; col < n; ++col) {
    ap[col] = (int32_t)pos;
    if (col < A_WIDTH) {
      for (int32_t row = 0; row < A_WIDTH; ++row) {
        const size_t p = pos++;
        ai[p] = row;
        ax0[p] = row == col
          ? 18.0 + 0.125 * (double)col
          : 0.0007 * (1.0 + (double)((row + 3 * col) % 17));
        ax1[p] = ax0[p] + (row == col ? 0.01 : 1.0e-6);
      }
      if (col == A_WIDTH - 4 || col == A_WIDTH - 3) {
        const size_t p = pos++;
        ai[p] = consumer;
        ax0[p] = 0.00011 * (1.0 + (double)col);
        ax1[p] = ax0[p] + 1.0e-6;
      }
    } else if (col < C_BEGIN) {
      for (int32_t row = B_BEGIN; row < C_BEGIN; ++row) {
        const size_t p = pos++;
        ai[p] = row;
        ax0[p] = row == col
          ? 16.0 + 0.1 * (double)col
          : 0.0009 * (1.0 + (double)((row + col) % 13));
        ax1[p] = ax0[p] + (row == col ? 0.0125 : 1.0e-6);
      }
      if (col == B_BEGIN + 1 || col == B_BEGIN + 2) {
        const size_t p = pos++;
        ai[p] = consumer;
        ax0[p] = 0.00013;
        ax1[p] = ax0[p] + 1.0e-6;
      }
    } else if (col < CONSUMER) {
      for (int32_t row = C_BEGIN; row < CONSUMER; ++row) {
        const size_t p = pos++;
        ai[p] = row;
        ax0[p] = row == col
          ? 15.0 + 0.075 * (double)col
          : 0.0008 * (1.0 + (double)((row + 2 * col) % 11));
        ax1[p] = ax0[p] + (row == col ? 0.011 : 1.0e-6);
      }
      if (col == C_BEGIN + 1) {
        const size_t p = pos++;
        ai[p] = consumer;
        ax0[p] = 0.00012;
        ax1[p] = ax0[p] + 1.0e-6;
      }
    } else {
      const size_t p = pos++;
      ai[p] = consumer;
      ax0[p] = 14.0;
      ax1[p] = 14.025;
    }
  }
  ap[n] = (int32_t)pos;
  if (pos != nnz) {
    fprintf(stderr, "unexpected prefactor fixture nnz: %zu/%zu\n", pos, nnz);
    return 0;
  }
  for (int32_t col = 0; col < n; ++col) {
    for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
      b[ai[p]] += ax1[p] * expected[col];
    }
  }

  const char *saved_row_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_row = saved_row_value != NULL ? strdup(saved_row_value) : NULL;
  const int had_row = saved_row_value != NULL;
  const char *saved_partial_value =
    getenv("KLS_ENABLE_PARTIAL_SUPERNODE_PIPELINE");
  char *saved_partial = saved_partial_value != NULL
    ? strdup(saved_partial_value) : NULL;
  const int had_partial = saved_partial_value != NULL;
  const char *saved_checked_value = getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  char *saved_checked = saved_checked_value != NULL
    ? strdup(saved_checked_value) : NULL;
  const int had_checked = saved_checked_value != NULL;
  const char *saved_trsv_value = getenv("KLS_ENABLE_COMPACT_SUPERNODE_TRSV");
  char *saved_trsv = saved_trsv_value != NULL ? strdup(saved_trsv_value) : NULL;
  const int had_trsv = saved_trsv_value != NULL;
  const char *saved_cblas_value = getenv("KLS_ENABLE_CBLAS_SUPERNODE");
  char *saved_cblas =
    saved_cblas_value != NULL ? strdup(saved_cblas_value) : NULL;
  const int had_cblas = saved_cblas_value != NULL;

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.threads = 3;
  options.ordering = KLS_ORDERING_NATURAL;
  options.orientation = KLS_ORIENTATION_NORMAL;
  options.use_btf = 0;
  options.scale = -1;
  options.static_pivoting = 0;

  int ok = 1;
  if ((had_row && saved_row == NULL) ||
      (had_partial && saved_partial == NULL) ||
      (had_checked && saved_checked == NULL) ||
      (had_trsv && saved_trsv == NULL) ||
      (had_cblas && saved_cblas == NULL)) {
    fprintf(stderr, "failed to save prefactor env\n");
    ok = 0;
  }
  if (!require_ok(kls_create(&solver), "create row-prefactor")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, n, ap, ai, 0,
                                        &options),
                        "analyze row-prefactor")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor row-prefactor base")) ok = 0;
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && unsetenv("KLS_ENABLE_PARTIAL_SUPERNODE_PIPELINE") != 0) {
    perror("unsetenv KLS_ENABLE_PARTIAL_SUPERNODE_PIPELINE");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR=1");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_COMPACT_SUPERNODE_TRSV", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_COMPACT_SUPERNODE_TRSV=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CBLAS_SUPERNODE", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CBLAS_SUPERNODE=0");
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(solver, ax1),
                        "checked factor row-prefactor")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats row-prefactor")) {
    ok = 0;
  }
  if (ok && (stats.row_refactor_last_run != 1 ||
             stats.row_refactor_last_checked != 1 ||
             stats.row_refactor_last_parallel != 1 ||
             stats.row_refactor_last_partial_supernode_pipeline != 1 ||
             stats.row_refactor_last_done_bitmap != 1 ||
             stats.row_refactor_last_prefactor != 0 ||
             stats.row_refactor_last_prefactor_rows != 0 ||
             stats.row_refactor_last_prefactor_deps != 0 ||
             stats.row_refactor_last_prefactor_supernode != 0 ||
             stats.row_refactor_last_prefactor_supernode_rows != 0 ||
             stats.row_refactor_last_prefactor_supernode_deps != 0 ||
             stats.row_refactor_prefactor_supernode_run_count != 0 ||
             stats.row_refactor_prefactor_supernode_rows != 0 ||
             stats.row_refactor_prefactor_supernode_deps != 0 ||
             stats.row_refactor_prefactor_run_count != 0 ||
             stats.row_refactor_prefactor_rows != 0 ||
             stats.row_refactor_prefactor_deps != 0)) {
    fprintf(stderr,
            "unexpected exact row-dependency stats: last=%d/%d/%d partial=%d"
            " done=%d prefactor=%d rows/deps=%" PRId64 "/%" PRId64
            " supernode=%d rows/deps=%" PRId64 "/%" PRId64
            " supernode_totals=%" PRId64 "/%" PRId64 "/%" PRId64
            " totals=%" PRId64 "/%" PRId64 "/%" PRId64 "\n",
            stats.row_refactor_last_run,
            stats.row_refactor_last_checked,
            stats.row_refactor_last_parallel,
            stats.row_refactor_last_partial_supernode_pipeline,
            stats.row_refactor_last_done_bitmap,
            stats.row_refactor_last_prefactor,
            stats.row_refactor_last_prefactor_rows,
            stats.row_refactor_last_prefactor_deps,
            stats.row_refactor_last_prefactor_supernode,
            stats.row_refactor_last_prefactor_supernode_rows,
            stats.row_refactor_last_prefactor_supernode_deps,
            stats.row_refactor_prefactor_supernode_run_count,
            stats.row_refactor_prefactor_supernode_rows,
            stats.row_refactor_prefactor_supernode_deps,
            stats.row_refactor_prefactor_run_count,
            stats.row_refactor_prefactor_rows,
            stats.row_refactor_prefactor_deps);
    ok = 0;
  }

  if (!restore_env_value("KLS_ENABLE_ROW_REFACTOR", had_row, saved_row)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_PARTIAL_SUPERNODE_PIPELINE",
                         had_partial, saved_partial)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_CHECKED_ROW_REFACTOR",
                         had_checked, saved_checked)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_COMPACT_SUPERNODE_TRSV",
                         had_trsv, saved_trsv)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_CBLAS_SUPERNODE", had_cblas,
                         saved_cblas)) {
    ok = 0;
  }
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve row-prefactor")) ok = 0;
  for (int32_t i = 0; ok && i < n; ++i) {
    if (!close_enough(x[i], expected[i])) {
      fprintf(stderr, "unexpected row-prefactor solution at %d: %.17g != %.17g\n",
              (int)i, x[i], expected[i]);
      ok = 0;
    }
  }

  kls_destroy(solver);
  free(saved_row);
  free(saved_partial);
  free(saved_checked);
  free(saved_trsv);
  free(saved_cblas);
  return ok;
}

static int test_partial_supernode_non_dominant_tail(void) {
  enum {
    A_WIDTH = 48,
    A_PREFIX = 24,
    CONSUMER_WIDTH = 3,
    SMALL_WIDTH = 2,
    SMALL_BLOCKS = 30,
    SMALL_CHAINS = 2,
    CONSUMER_BEGIN = A_WIDTH,
    SMALL_BEGIN = CONSUMER_BEGIN + CONSUMER_WIDTH,
    SMALL_CHAIN_WIDTH = SMALL_WIDTH * SMALL_BLOCKS,
    NONDOM_N = SMALL_BEGIN + SMALL_CHAIN_WIDTH * SMALL_CHAINS
  };
  const int32_t n = NONDOM_N;
  const size_t nnz =
    (size_t)A_WIDTH * (size_t)A_WIDTH +
    (size_t)A_PREFIX * (size_t)CONSUMER_WIDTH +
    (size_t)CONSUMER_WIDTH * (size_t)CONSUMER_WIDTH +
    (size_t)CONSUMER_WIDTH * (size_t)SMALL_WIDTH +
    (size_t)SMALL_CHAINS * (size_t)SMALL_BLOCKS *
      (size_t)SMALL_WIDTH * (size_t)SMALL_WIDTH +
    (size_t)SMALL_CHAINS * (size_t)(SMALL_BLOCKS - 1) *
      (size_t)SMALL_WIDTH * (size_t)SMALL_WIDTH;
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

  for (int32_t col = 0; col < n; ++col) {
    expected[col] = 0.5 + 0.04 * (double)((3 * col) % 17);
  }

  size_t pos = 0;
  for (int32_t col = 0; col < n; ++col) {
    ap[col] = (int32_t)pos;
    if (col < A_WIDTH) {
      for (int32_t row = 0; row < A_WIDTH; ++row) {
        const size_t p = pos++;
        ai[p] = row;
        ax0[p] = row == col
          ? 26.0 + 0.01 * (double)col
          : 0.00018 * (1.0 + (double)((row + 5 * col) % 19));
        ax1[p] = ax0[p] + (row == col ? 0.01 : 1.0e-6);
      }
      if (col < A_PREFIX) {
        for (int32_t row = CONSUMER_BEGIN;
             row < CONSUMER_BEGIN + CONSUMER_WIDTH; ++row) {
          const size_t p = pos++;
          ai[p] = row;
          ax0[p] = 0.00011 * (1.0 + (double)((row + col) % 11));
          ax1[p] = ax0[p] + 1.0e-6;
        }
      }
    } else if (col < SMALL_BEGIN) {
      for (int32_t row = CONSUMER_BEGIN;
           row < CONSUMER_BEGIN + CONSUMER_WIDTH; ++row) {
        const size_t p = pos++;
        ai[p] = row;
        ax0[p] = row == col
          ? 21.0 + 0.01 * (double)(col - CONSUMER_BEGIN)
          : 0.0002 * (1.0 + (double)((row + col) % 7));
        ax1[p] = ax0[p] + (row == col ? 0.01 : 1.0e-6);
      }
      for (int32_t row = SMALL_BEGIN; row < SMALL_BEGIN + SMALL_WIDTH;
           ++row) {
        const size_t p = pos++;
        ai[p] = row;
        ax0[p] = 0.00013 * (1.0 + (double)((row + col) % 5));
        ax1[p] = ax0[p] + 1.0e-6;
      }
    } else {
      const int32_t small_offset = col - SMALL_BEGIN;
      const int32_t chain = small_offset / SMALL_CHAIN_WIDTH;
      const int32_t chain_offset = small_offset - chain * SMALL_CHAIN_WIDTH;
      const int32_t block = chain_offset / SMALL_WIDTH;
      const int32_t block_begin =
        SMALL_BEGIN + chain * SMALL_CHAIN_WIDTH + block * SMALL_WIDTH;
      for (int32_t row = block_begin; row < block_begin + SMALL_WIDTH;
           ++row) {
        const size_t p = pos++;
        ai[p] = row;
        ax0[p] = row == col
          ? 18.0 + 0.02 * (double)block
          : 0.00021 * (1.0 + (double)((row + col) % 5));
        ax1[p] = ax0[p] + (row == col ? 0.01 : 1.0e-6);
      }
      if (block + 1 < SMALL_BLOCKS) {
        const int32_t next_begin = block_begin + SMALL_WIDTH;
        for (int32_t row = next_begin; row < next_begin + SMALL_WIDTH;
             ++row) {
          const size_t p = pos++;
          ai[p] = row;
          ax0[p] = 0.00012 * (1.0 + (double)((row + col) % 7));
          ax1[p] = ax0[p] + 1.0e-6;
        }
      }
    }
  }
  ap[n] = (int32_t)pos;
  if (pos != nnz) {
    fprintf(stderr, "unexpected non-dominant partial fixture nnz: %zu/%zu\n",
            pos, nnz);
    free(ap);
    free(ai);
    free(ax0);
    free(ax1);
    free(b);
    free(x);
    free(expected);
    return 0;
  }
  for (int32_t col = 0; col < n; ++col) {
    for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
      b[ai[p]] += ax1[p] * expected[col];
    }
  }

  const char *saved_row_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_row = saved_row_value != NULL ? strdup(saved_row_value) : NULL;
  const int had_row = saved_row_value != NULL;
  const char *saved_partial_value =
    getenv("KLS_ENABLE_PARTIAL_SUPERNODE_PIPELINE");
  char *saved_partial = saved_partial_value != NULL
    ? strdup(saved_partial_value) : NULL;
  const int had_partial = saved_partial_value != NULL;
  const char *saved_checked_value = getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  char *saved_checked = saved_checked_value != NULL
    ? strdup(saved_checked_value) : NULL;
  const int had_checked = saved_checked_value != NULL;

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.threads = 2;
  options.ordering = KLS_ORDERING_NATURAL;
  options.orientation = KLS_ORIENTATION_NORMAL;
  options.use_btf = 0;
  options.scale = -1;
  options.static_pivoting = 0;

  int ok = 1;
  if ((had_row && saved_row == NULL) ||
      (had_partial && saved_partial == NULL) ||
      (had_checked && saved_checked == NULL)) {
    fprintf(stderr, "failed to save non-dominant partial env\n");
    ok = 0;
  }
  if (!require_ok(kls_create(&solver), "create non-dominant partial")) {
    ok = 0;
  }
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, n, ap, ai, 0,
                                        &options),
                        "analyze non-dominant partial")) {
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor non-dominant partial base")) {
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR=1");
    ok = 0;
  }
  if (ok && unsetenv("KLS_ENABLE_PARTIAL_SUPERNODE_PIPELINE") != 0) {
    perror("unsetenv KLS_ENABLE_PARTIAL_SUPERNODE_PIPELINE");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && !require_ok(kls_refactor(solver, ax1),
                        "refactor non-dominant partial")) {
    ok = 0;
  }

  if (!restore_env_value("KLS_ENABLE_ROW_REFACTOR", had_row, saved_row)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_PARTIAL_SUPERNODE_PIPELINE",
                         had_partial, saved_partial)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_CHECKED_ROW_REFACTOR",
                         had_checked, saved_checked)) {
    ok = 0;
  }

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats non-dominant partial")) {
    ok = 0;
  }
  if (ok && (stats.row_refactor_last_parallel != 1 ||
             stats.row_refactor_last_partial_supernode_pipeline != 1 ||
             stats.row_refactor_last_done_bitmap != 1 ||
             stats.row_refactor_last_partial_supernode_pipeline_groups != 1 ||
             stats.row_refactor_last_partial_supernode_pipeline_rows !=
               A_WIDTH ||
             stats.row_refactor_group_count <=
               2 * stats.row_refactor_last_partial_supernode_pipeline_groups ||
             stats.row_refactor_group_pipeline_rows <=
               2 * stats.row_refactor_last_partial_supernode_pipeline_rows)) {
    fprintf(stderr,
            "unexpected non-dominant partial stats: parallel=%d"
            ", partial=%d groups/rows=%" PRId64 "/%" PRId64
            ", group_count=%" PRId64 ", pipe_rows=%" PRId64
            ", done=%d\n",
            stats.row_refactor_last_parallel,
            stats.row_refactor_last_partial_supernode_pipeline,
            stats.row_refactor_last_partial_supernode_pipeline_groups,
            stats.row_refactor_last_partial_supernode_pipeline_rows,
            stats.row_refactor_group_count,
            stats.row_refactor_group_pipeline_rows,
            stats.row_refactor_last_done_bitmap);
    ok = 0;
  }

  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve non-dominant partial")) {
    ok = 0;
  }
  for (int32_t i = 0; ok && i < n; ++i) {
    if (!close_enough(x[i], expected[i])) {
      fprintf(stderr,
              "unexpected non-dominant partial solution at %d:"
              " %.17g != %.17g\n",
              (int)i, x[i], expected[i]);
      ok = 0;
    }
  }

  kls_destroy(solver);
  free(saved_row);
  free(saved_partial);
  free(saved_checked);
  free(ap);
  free(ai);
  free(ax0);
  free(ax1);
  free(b);
  free(x);
  free(expected);
  return ok;
}

static int test_partial_compact_supernode_prefix_pipeline(void) {
  const int32_t a_width = 240;
  const int32_t a_prefix = 120;
  const int32_t a_consumer = 32;
  const int32_t b_width = 8;
  const int32_t b_prefix = 4;
  const int32_t b_consumer = 8;
  const int32_t a_begin = 0;
  const int32_t a_consumer_begin = a_begin + a_width;
  const int32_t b_begin = a_consumer_begin + a_consumer;
  const int32_t b_consumer_begin = b_begin + b_width;
  const int32_t n = b_consumer_begin + b_consumer;
  const size_t nnz =
    (size_t)a_width * (size_t)a_width +
    (size_t)a_prefix * (size_t)a_consumer +
    (size_t)a_consumer * (size_t)a_consumer +
    (size_t)b_width * (size_t)b_width +
    (size_t)b_prefix * (size_t)b_consumer +
    (size_t)b_consumer * (size_t)b_consumer;
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

  for (int32_t col = 0; col < n; ++col) {
    expected[col] = 0.625 + 0.03125 * (double)((7 * col) % 23);
  }

  size_t pos = 0;
  for (int32_t col = 0; col < n; ++col) {
    ap[col] = (int32_t)pos;
    if (col < a_begin + a_width) {
      for (int32_t row = a_begin; row < a_begin + a_width; ++row) {
        const size_t p = pos++;
        ai[p] = row;
        ax0[p] = row == col
          ? 30.0 + 0.003 * (double)col
          : 0.00012 * (1.0 + (double)((row + 7 * col) % 29));
        ax1[p] = ax0[p] + (row == col
          ? 0.02 * (double)((col % 5) + 1)
          : 1.0e-6 * (double)(((row + col) % 7) - 3));
      }
      if (col < a_begin + a_prefix) {
        for (int32_t row = a_consumer_begin;
             row < a_consumer_begin + a_consumer; ++row) {
          const size_t p = pos++;
          ai[p] = row;
          ax0[p] = 0.00009 * (1.0 + (double)((row + 3 * col) % 23));
          ax1[p] = ax0[p] +
            1.0e-6 * (double)(((row + 2 * col) % 5) - 2);
        }
      }
    } else if (col < a_consumer_begin + a_consumer) {
      for (int32_t row = a_consumer_begin;
           row < a_consumer_begin + a_consumer; ++row) {
        const size_t p = pos++;
        ai[p] = row;
        ax0[p] = row == col
          ? 24.0 + 0.004 * (double)col
          : 0.0001 * (1.0 + (double)((row + 5 * col) % 17));
        ax1[p] = ax0[p] + (row == col
          ? 0.015 * (double)((col % 3) + 1)
          : 1.0e-6 * (double)(((row + col) % 5) - 2));
      }
    } else if (col < b_begin + b_width) {
      for (int32_t row = b_begin; row < b_begin + b_width; ++row) {
        const size_t p = pos++;
        ai[p] = row;
        ax0[p] = row == col
          ? 22.0 + 0.01 * (double)col
          : 0.0002 * (1.0 + (double)((row + 7 * col) % 11));
        ax1[p] = ax0[p] + (row == col ? 0.01 : 1.0e-6);
      }
      if (col < b_begin + b_prefix) {
        for (int32_t row = b_consumer_begin;
             row < b_consumer_begin + b_consumer; ++row) {
          const size_t p = pos++;
          ai[p] = row;
          ax0[p] = 0.00015 * (1.0 + (double)((row + 3 * col) % 13));
          ax1[p] = ax0[p] + 1.0e-6;
        }
      }
    } else {
      for (int32_t row = b_consumer_begin;
           row < b_consumer_begin + b_consumer; ++row) {
        const size_t p = pos++;
        ai[p] = row;
        ax0[p] = row == col
          ? 20.0 + 0.005 * (double)col
          : 0.00012 * (1.0 + (double)((row + 5 * col) % 7));
        ax1[p] = ax0[p] + (row == col ? 0.01 : 1.0e-6);
      }
    }
  }
  ap[n] = (int32_t)pos;
  if (pos != nnz) {
    fprintf(stderr, "unexpected partial-prefix fixture nnz: %zu/%zu\n",
            pos, nnz);
    free(ap);
    free(ai);
    free(ax0);
    free(ax1);
    free(b);
    free(x);
    free(expected);
    return 0;
  }
  for (int32_t col = 0; col < n; ++col) {
    for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
      b[ai[p]] += ax1[p] * expected[col];
    }
  }

  const char *saved_row_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_row = saved_row_value != NULL ? strdup(saved_row_value) : NULL;
  const int had_row = saved_row_value != NULL;
  const char *saved_partial_value =
    getenv("KLS_ENABLE_PARTIAL_SUPERNODE_PIPELINE");
  char *saved_partial = saved_partial_value != NULL
    ? strdup(saved_partial_value) : NULL;
  const int had_partial = saved_partial_value != NULL;
  const char *saved_checked_value = getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  char *saved_checked = saved_checked_value != NULL
    ? strdup(saved_checked_value) : NULL;
  const int had_checked = saved_checked_value != NULL;
  const char *saved_trsv_value = getenv("KLS_ENABLE_COMPACT_SUPERNODE_TRSV");
  char *saved_trsv = saved_trsv_value != NULL ? strdup(saved_trsv_value) : NULL;
  const int had_trsv = saved_trsv_value != NULL;
  const char *saved_cblas_value = getenv("KLS_ENABLE_CBLAS_SUPERNODE");
  char *saved_cblas =
    saved_cblas_value != NULL ? strdup(saved_cblas_value) : NULL;
  const int had_cblas = saved_cblas_value != NULL;
  const char *saved_native_value =
    getenv("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR");
  char *saved_native =
    saved_native_value != NULL ? strdup(saved_native_value) : NULL;
  const int had_native = saved_native_value != NULL;

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.threads = 2;
  options.ordering = KLS_ORDERING_NATURAL;
  options.orientation = KLS_ORIENTATION_NORMAL;
  options.use_btf = 0;
  options.scale = -1;
  options.static_pivoting = 0;

  int ok = 1;
  if ((had_row && saved_row == NULL) ||
      (had_partial && saved_partial == NULL) ||
      (had_checked && saved_checked == NULL) ||
      (had_trsv && saved_trsv == NULL) ||
      (had_cblas && saved_cblas == NULL) ||
      (had_native && saved_native == NULL)) {
    fprintf(stderr, "failed to save partial-prefix env\n");
    ok = 0;
  }
  if (!require_ok(kls_create(&solver), "create partial-prefix")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, n, ap, ai, 0,
                                        &options),
                        "analyze partial-prefix")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor partial-prefix base")) ok = 0;
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && unsetenv("KLS_ENABLE_PARTIAL_SUPERNODE_PIPELINE") != 0) {
    perror("unsetenv KLS_ENABLE_PARTIAL_SUPERNODE_PIPELINE");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_COMPACT_SUPERNODE_TRSV", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_COMPACT_SUPERNODE_TRSV=1");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CBLAS_SUPERNODE", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CBLAS_SUPERNODE=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR", "auto", 1) != 0) {
    perror("setenv KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR=auto");
    ok = 0;
  }
  if (ok && !require_ok(kls_refactor(solver, ax1),
                        "refactor partial-prefix")) ok = 0;

  if (!restore_env_value("KLS_ENABLE_ROW_REFACTOR", had_row, saved_row)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_PARTIAL_SUPERNODE_PIPELINE",
                         had_partial, saved_partial)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_CHECKED_ROW_REFACTOR",
                         had_checked, saved_checked)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_COMPACT_SUPERNODE_TRSV",
                         had_trsv, saved_trsv)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_CBLAS_SUPERNODE",
                         had_cblas, saved_cblas)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR",
                         had_native, saved_native)) {
    ok = 0;
  }
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve partial-prefix")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats partial-prefix")) {
    ok = 0;
  }
  const int64_t expected_partial_rows =
    (int64_t)(a_width - options.threads);
  if (ok && (stats.row_refactor_last_parallel != 1 ||
             stats.row_refactor_last_partial_supernode_pipeline != 1 ||
             stats.row_refactor_last_done_bitmap != 1 ||
             stats.row_refactor_last_compact_dense_panel != 1 ||
             stats.row_refactor_last_compact_supernode_update != 1 ||
             stats.row_refactor_last_compact_supernode_trsv != 1 ||
             stats.row_refactor_compact_supernode_update_rows <
               expected_partial_rows ||
             stats.row_refactor_compact_supernode_trsv_rows <
               expected_partial_rows)) {
    fprintf(stderr,
            "unexpected partial-prefix stats: parallel=%d partial=%d done=%d"
            ", compact=%d, update=%d/%" PRId64
            ", trsv=%d/%" PRId64 "/%" PRId64
            ", partial_update=%d/%" PRId64 "/%" PRId64 "/%" PRId64 "\n",
            stats.row_refactor_last_parallel,
            stats.row_refactor_last_partial_supernode_pipeline,
            stats.row_refactor_last_done_bitmap,
            stats.row_refactor_last_compact_dense_panel,
            stats.row_refactor_last_compact_supernode_update,
            stats.row_refactor_compact_supernode_update_rows,
            stats.row_refactor_last_compact_supernode_trsv,
            stats.row_refactor_compact_supernode_trsv_rows,
            stats.row_refactor_compact_supernode_trsv_entries,
            stats.row_refactor_last_compact_supernode_partial_update,
            stats.row_refactor_compact_supernode_partial_update_count,
            stats.row_refactor_compact_supernode_partial_update_rows,
            stats.row_refactor_compact_supernode_partial_update_entries);
    ok = 0;
  }

  double max_solution_error = 0.0;
  double max_residual = 0.0;
  double max_rhs = 0.0;
  for (int32_t row = 0; row < n; ++row) {
    const double err = fabs(x[row] - expected[row]);
    if (err > max_solution_error) {
      max_solution_error = err;
    }
    double residual = -b[row];
    if (fabs(b[row]) > max_rhs) {
      max_rhs = fabs(b[row]);
    }
    for (int32_t col = 0; col < n; ++col) {
      for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
        if (ai[p] == row) {
          residual += ax1[p] * x[col];
        }
      }
    }
    if (fabs(residual) > max_residual) {
      max_residual = fabs(residual);
    }
  }
  const double rel_resid = max_residual / (max_rhs > 0.0 ? max_rhs : 1.0);
  if (ok && (max_solution_error > 1.0e-8 || rel_resid > 1.0e-10)) {
    fprintf(stderr,
            "unexpected partial-prefix accuracy: max_x_err=%.17g"
            ", rel_resid=%.17g\n",
            max_solution_error, rel_resid);
    ok = 0;
  }

  kls_solver *checked_solver = NULL;
  if (ok && !require_ok(kls_create(&checked_solver),
                        "create checked partial-prefix")) {
    ok = 0;
  }
  if (ok && !require_ok(kls_analyze_csc(checked_solver, KLS_INDEX_INT32,
                                        n, ap, ai, 0, &options),
                        "analyze checked partial-prefix")) {
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(checked_solver, ax0),
                        "factor checked partial-prefix base")) {
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && unsetenv("KLS_ENABLE_PARTIAL_SUPERNODE_PIPELINE") != 0) {
    perror("unsetenv KLS_ENABLE_PARTIAL_SUPERNODE_PIPELINE");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR=1");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_COMPACT_SUPERNODE_TRSV", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_COMPACT_SUPERNODE_TRSV=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CBLAS_SUPERNODE", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_CBLAS_SUPERNODE=1");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR", "auto", 1) != 0) {
    perror("setenv KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR=auto");
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(checked_solver, ax1),
                        "checked factor partial-prefix")) {
    ok = 0;
  }

  if (!restore_env_value("KLS_ENABLE_ROW_REFACTOR", had_row, saved_row)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_PARTIAL_SUPERNODE_PIPELINE",
                         had_partial, saved_partial)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_CHECKED_ROW_REFACTOR",
                         had_checked, saved_checked)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_COMPACT_SUPERNODE_TRSV",
                         had_trsv, saved_trsv)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_CBLAS_SUPERNODE",
                         had_cblas, saved_cblas)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR",
                         had_native, saved_native)) {
    ok = 0;
  }
  if (ok && !require_ok(kls_solve(checked_solver, 1, b, 0, x, 0),
                        "solve checked partial-prefix")) {
    ok = 0;
  }

  kls_stats checked_stats;
  checked_stats.struct_size = sizeof(checked_stats);
  if (ok && !require_ok(kls_get_stats(checked_solver, &checked_stats),
                        "stats checked partial-prefix")) {
    ok = 0;
  }
  if (ok && (checked_stats.row_refactor_last_run != 1 ||
             checked_stats.row_refactor_last_checked != 1 ||
             checked_stats.row_refactor_last_parallel != 1 ||
             checked_stats.row_refactor_last_partial_supernode_pipeline != 1 ||
             checked_stats.row_refactor_last_done_bitmap != 1 ||
             checked_stats.row_refactor_last_compact_dense_panel != 1 ||
             checked_stats.row_refactor_last_compact_dense_panel_blocked != 1 ||
             checked_stats.row_refactor_compact_dense_panel_blocked_run_count < 1 ||
             checked_stats.row_refactor_compact_dense_panel_blocked_rows <
               expected_partial_rows ||
             checked_stats.row_refactor_compact_dense_panel_blocked_entries <= 0 ||
             checked_stats.row_refactor_last_compact_supernode_update != 1 ||
             checked_stats.row_refactor_compact_supernode_update_rows <
               expected_partial_rows)) {
    fprintf(stderr,
            "unexpected checked partial-prefix stats: last=%d/%d/%d"
            ", partial=%d, done=%d, compact=%d, blocked=%d/%" PRId64
            "/%" PRId64 "/%" PRId64 ", update=%d/%" PRId64
            ", partial_update=%d/%" PRId64 "/%" PRId64 "/%" PRId64 "\n",
            checked_stats.row_refactor_last_run,
            checked_stats.row_refactor_last_checked,
            checked_stats.row_refactor_last_parallel,
            checked_stats.row_refactor_last_partial_supernode_pipeline,
            checked_stats.row_refactor_last_done_bitmap,
            checked_stats.row_refactor_last_compact_dense_panel,
            checked_stats.row_refactor_last_compact_dense_panel_blocked,
            checked_stats.row_refactor_compact_dense_panel_blocked_run_count,
            checked_stats.row_refactor_compact_dense_panel_blocked_rows,
            checked_stats.row_refactor_compact_dense_panel_blocked_entries,
            checked_stats.row_refactor_last_compact_supernode_update,
            checked_stats.row_refactor_compact_supernode_update_rows,
            checked_stats.row_refactor_last_compact_supernode_partial_update,
            checked_stats.row_refactor_compact_supernode_partial_update_count,
            checked_stats.row_refactor_compact_supernode_partial_update_rows,
            checked_stats.row_refactor_compact_supernode_partial_update_entries);
    ok = 0;
  }

  max_solution_error = 0.0;
  max_residual = 0.0;
  max_rhs = 0.0;
  for (int32_t row = 0; row < n; ++row) {
    const double err = fabs(x[row] - expected[row]);
    if (err > max_solution_error) {
      max_solution_error = err;
    }
    double residual = -b[row];
    if (fabs(b[row]) > max_rhs) {
      max_rhs = fabs(b[row]);
    }
    for (int32_t col = 0; col < n; ++col) {
      for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
        if (ai[p] == row) {
          residual += ax1[p] * x[col];
        }
      }
    }
    if (fabs(residual) > max_residual) {
      max_residual = fabs(residual);
    }
  }
  const double checked_rel_resid =
    max_residual / (max_rhs > 0.0 ? max_rhs : 1.0);
  if (ok && (max_solution_error > 1.0e-8 ||
             checked_rel_resid > 1.0e-10)) {
    fprintf(stderr,
            "unexpected checked partial-prefix accuracy: max_x_err=%.17g"
            ", rel_resid=%.17g\n",
            max_solution_error, checked_rel_resid);
    ok = 0;
  }

  kls_destroy(solver);
  kls_destroy(checked_solver);
  free(saved_row);
  free(saved_partial);
  free(saved_checked);
  free(saved_trsv);
  free(saved_cblas);
  free(saved_native);
  free(ap);
  free(ai);
  free(ax0);
  free(ax1);
  free(b);
  free(x);
  free(expected);
  return ok;
}

static int test_row_refactor_off_blocks_auto_kls_first_refactor(void) {
  const int32_t ap[] = {0, 2, 4, 5};
  const int32_t ai[] = {0, 1, 0, 1, 2};
  const double ax[] = {4.0, 2.0, 1.0, 3.0, 5.0};
  const double ax_ref[] = {5.0, 1.0, 2.0, 4.0, 6.0};

  const char *saved_first_value = getenv("KLS_ENABLE_KLS_FIRST_FACTOR");
  char *saved_first =
    saved_first_value != NULL ? strdup(saved_first_value) : NULL;
  const int had_saved_first = saved_first_value != NULL;
  const char *saved_row_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_row =
    saved_row_value != NULL ? strdup(saved_row_value) : NULL;
  const int had_saved_row = saved_row_value != NULL;
  const char *saved_checked_value =
    getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  char *saved_checked =
    saved_checked_value != NULL ? strdup(saved_checked_value) : NULL;
  const int had_saved_checked = saved_checked_value != NULL;

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_NATURAL;
  options.use_btf = 1;
  options.scale = -1;
  options.static_pivoting = 0;

  int ok = 1;
  if (had_saved_first && saved_first == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_KLS_FIRST_FACTOR\n");
    ok = 0;
  }
  if (had_saved_row && saved_row == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_ROW_REFACTOR\n");
    ok = 0;
  }
  if (had_saved_checked && saved_checked == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_CHECKED_ROW_REFACTOR\n");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_KLS_FIRST_FACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_KLS_FIRST_FACTOR");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR=0");
    ok = 0;
  }

  if (ok && !require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 3, ap, ai, 0,
                                        &options),
                        "analyze row-refactor off auto gate")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax),
                        "factor row-refactor off auto gate")) ok = 0;
  if (ok && !require_ok(kls_refactor(solver, ax_ref),
                        "refactor row-refactor off auto gate")) ok = 0;

  kls_stats stats;
  memset(&stats, 0, sizeof(stats));
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats row-refactor off auto gate")) {
    ok = 0;
  }
  if (ok && (stats.row_refactor_last_run != 0 ||
             stats.row_refactor_run_count != 0 ||
             stats.row_refactor_auto_should_run != 0)) {
    fprintf(stderr,
            "row-refactor off did not block auto path: run=%d"
            ", run_count=%" PRId64 ", auto_should=%d"
            ", auto_enabled=%d, auto_ready=%d\n",
            stats.row_refactor_last_run,
            stats.row_refactor_run_count,
            stats.row_refactor_auto_should_run,
            stats.row_refactor_auto_enabled,
            stats.row_refactor_auto_values_ready);
    ok = 0;
  }

  if (!restore_env_value("KLS_ENABLE_KLS_FIRST_FACTOR", had_saved_first,
                         saved_first)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_ROW_REFACTOR", had_saved_row,
                         saved_row)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_CHECKED_ROW_REFACTOR",
                         had_saved_checked, saved_checked)) {
    ok = 0;
  }

  kls_destroy(solver);
  free(saved_first);
  free(saved_row);
  free(saved_checked);
  return ok;
}

static int test_default_keeps_kls_first_factor_off(void) {
  const int32_t n = 30000;
  int32_t *ap = (int32_t *)malloc(((size_t)n + 1u) * sizeof(*ap));
  int32_t *ai = (int32_t *)malloc((size_t)n * sizeof(*ai));
  double *ax = (double *)malloc((size_t)n * sizeof(*ax));
  double *b = (double *)calloc((size_t)n, sizeof(*b));
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
    const double diag = 2.0 + (double)(col % 11);
    ai[col] = col;
    ax[col] = diag;
    expected[col] = 1.0 + (double)(col % 13);
    b[col] = diag * expected[col];
  }

  const char *saved_first_value = getenv("KLS_ENABLE_KLS_FIRST_FACTOR");
  char *saved_first =
    saved_first_value != NULL ? strdup(saved_first_value) : NULL;
  const int had_saved_first = saved_first_value != NULL;

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_NATURAL;
  options.use_btf = 0;
  options.scale = 0;
  options.static_pivoting = 0;

  int ok = 1;
  if (had_saved_first && saved_first == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_KLS_FIRST_FACTOR\n");
    ok = 0;
  }
  if (ok && unsetenv("KLS_ENABLE_KLS_FIRST_FACTOR") != 0) {
    perror("unsetenv KLS_ENABLE_KLS_FIRST_FACTOR");
    ok = 0;
  }
  if (ok && !require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, n, ap, ai, 0,
                                        &options),
                        "analyze default KLU first factor")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax),
                        "factor default KLU first factor")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve default KLU first factor")) ok = 0;

  kls_stats stats;
  memset(&stats, 0, sizeof(stats));
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats default KLU first factor")) ok = 0;
  if (ok && (stats.last_factor_path != KLS_FACTOR_PATH_KLU_FIRST ||
             stats.kls_first_last_row_uplooking_columns != 0 ||
             stats.kls_first_row_uplooking_column_count != 0 ||
             stats.kls_first_last_row_refactor_seeded_rows != 0 ||
             stats.kls_first_row_refactor_seeded_row_count != 0 ||
             stats.selected_btf != 0 || stats.selected_scale != 0)) {
    fprintf(stderr,
            "unexpected default first-factor stats: path=%s, row_cols=%" PRId64
            "/%" PRId64 ", row_seed=%" PRId64 "/%" PRId64
            ", btf=%d, scale=%d\n",
            kls_factor_path_name(stats.last_factor_path),
            stats.kls_first_last_row_uplooking_columns,
            stats.kls_first_row_uplooking_column_count,
            stats.kls_first_last_row_refactor_seeded_rows,
            stats.kls_first_row_refactor_seeded_row_count,
            stats.selected_btf, stats.selected_scale);
    ok = 0;
  }
  for (int32_t i = 0; ok && i < n; ++i) {
    if (!close_enough(x[i], expected[i])) {
      fprintf(stderr,
              "unexpected default first-factor solution at %d:"
              " %.17g != %.17g\n",
              (int)i, x[i], expected[i]);
      ok = 0;
    }
  }

  if (!restore_env_value("KLS_ENABLE_KLS_FIRST_FACTOR", had_saved_first,
                         saved_first)) {
    ok = 0;
  }
  kls_destroy(solver);
  free(ap);
  free(ai);
  free(ax);
  free(b);
  free(x);
  free(expected);
  free(saved_first);
  return ok;
}

static int test_kls_first_separator_queue_plan(void) {
  /* banded matrix, halfband 16: any reasonable nested-dissection
     separator of a band graph is ~bandwidth wide, so the separator
     pipeline segment carries supernode panels for EVERY ordering the
     RNG produces (the previous tridiagonal fixture only exercised the
     panel path when an imperfect separator happened to be fat). */
  const int32_t n = 30000;
  const int32_t half_band = 16;
  const int32_t nnz_cap = (2 * half_band + 1) * n;
  int32_t *ap = (int32_t *)malloc(((size_t)n + 1u) * sizeof(*ap));
  int32_t *ai = (int32_t *)malloc((size_t)nnz_cap * sizeof(*ai));
  double *ax = (double *)malloc((size_t)nnz_cap * sizeof(*ax));
  double *b = (double *)calloc((size_t)n, sizeof(*b));
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

  int32_t p = 0;
  for (int32_t col = 0; col < n; ++col) {
    ap[col] = p;
    expected[col] = 1.0 + (double)(col % 5);
    const int32_t lo = col - half_band < 0 ? 0 : col - half_band;
    const int32_t hi =
      col + half_band >= n ? n - 1 : col + half_band;
    for (int32_t row = lo; row <= hi; ++row) {
      ai[p] = row;
      /* diagonally dominant: |diag| > 2*half_band*|offdiag| */
      ax[p] = row == col ? 40.0 : -1.0;
      p++;
    }
  }
  ap[n] = p;
  if (p > nnz_cap) {
    fprintf(stderr, "unexpected separator queue nnz: %d/%d\n", p, nnz_cap);
    free(ap);
    free(ai);
    free(ax);
    free(b);
    free(x);
    free(expected);
    return 0;
  }
  for (int32_t col = 0; col < n; ++col) {
    for (int32_t q = ap[col]; q < ap[col + 1]; ++q) {
      b[ai[q]] += ax[q] * expected[col];
    }
  }

  const char *saved_first_value = getenv("KLS_ENABLE_KLS_FIRST_FACTOR");
  char *saved_first =
    saved_first_value != NULL ? strdup(saved_first_value) : NULL;
  const int had_saved_first = saved_first_value != NULL;

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_METIS;
  options.use_btf = 0;
  options.scale = 0;
  options.static_pivoting = 0;
  options.threads = 4;

  int ok = 1;
  if (had_saved_first && saved_first == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_KLS_FIRST_FACTOR\n");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_KLS_FIRST_FACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_KLS_FIRST_FACTOR");
    ok = 0;
  }
  if (ok && !require_ok(kls_create(&solver), "create separator queue plan")) {
    ok = 0;
  }
  if (ok) {
    const int analyze_status =
      kls_analyze_csc(solver, KLS_INDEX_INT32, n, ap, ai, 0, &options);
    if (analyze_status == KLS_ERR_UNSUPPORTED) {
      goto cleanup;
    }
    if (!require_ok(analyze_status, "analyze separator queue plan")) {
      ok = 0;
    }
  }
  if (ok && !require_ok(kls_factor(solver, ax),
                        "factor separator queue plan")) {
    ok = 0;
  }
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve separator queue plan")) {
    ok = 0;
  }

  kls_stats stats;
  memset(&stats, 0, sizeof(stats));
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats separator queue plan")) {
    ok = 0;
  }
  if (ok && (stats.build_has_metis != 1 ||
             stats.last_factor_path != KLS_FACTOR_PATH_KLS_FIRST ||
             stats.separator_analyzed_rows != n ||
             stats.separator_component_count <= 1 ||
             stats.kls_first_last_separator_queue != 1 ||
             stats.kls_first_separator_queue_run_count < 1 ||
             stats.kls_first_last_separator_queue_private_components <= 0 ||
             stats.kls_first_last_separator_queue_private_rows <= 0 ||
             stats.kls_first_last_separator_queue_pipeline_rows <= 0 ||
             stats.kls_first_last_separator_queue_private_rows +
               stats.kls_first_last_separator_queue_pipeline_rows != n ||
             stats.kls_first_last_separator_queue_nonempty_threads <= 1 ||
             stats.kls_first_last_separator_queue_max_thread_rows <= 0 ||
             stats.kls_first_last_separator_queue_min_thread_work <= 0.0 ||
             stats.kls_first_last_separator_queue_max_thread_work <
               stats.kls_first_last_separator_queue_min_thread_work ||
             stats.kls_first_last_separator_queue_partitioned != 1 ||
             stats.kls_first_separator_queue_partitioned_count < 1 ||
             stats.kls_first_last_separator_queue_split_components <= 0 ||
             stats.kls_first_last_separator_queue_executed != 1 ||
             stats.kls_first_separator_queue_executed_run_count < 1 ||
             stats.kls_first_last_separator_queue_executed_private_rows <= 0 ||
             stats.kls_first_last_separator_queue_executed_private_rows +
               stats.kls_first_last_separator_queue_executed_pipeline_rows !=
                 n ||
             stats.kls_first_last_separator_queue_parallel_private != 1 ||
             stats.kls_first_separator_queue_parallel_private_run_count < 1 ||
             stats.kls_first_last_separator_queue_parallel_private_rows <= 0 ||
             stats.kls_first_last_separator_queue_parallel_private_threads <=
               1 ||
             stats.kls_first_last_separator_queue_parallel_pipeline != 1 ||
             stats.kls_first_separator_queue_parallel_pipeline_run_count < 1 ||
             stats.kls_first_last_separator_queue_parallel_pipeline_rows !=
               stats.kls_first_last_separator_queue_pipeline_rows ||
             stats.kls_first_last_separator_queue_parallel_pipeline_threads <=
               0 ||
             stats.kls_first_last_separator_queue_pipeline_partial != 1 ||
             stats.kls_first_separator_queue_pipeline_partial_run_count < 1 ||
             stats.kls_first_last_separator_queue_pipeline_partial_rows !=
               stats.kls_first_last_separator_queue_pipeline_rows ||
             stats.kls_first_last_separator_queue_pipeline_partial_threads <=
               0 ||
             stats.kls_first_last_separator_queue_pipeline_supernode_panel_update !=
               1 ||
             stats.kls_first_last_separator_queue_pipeline_supernode_panel_update_groups <=
               0 ||
             stats.kls_first_last_separator_queue_pipeline_supernode_panel_update_rows <=
               0 ||
             stats.kls_first_last_row_uplooking_columns != n ||
             stats.selected_btf != 0 || stats.selected_scale != 0)) {
    fprintf(stderr,
            "unexpected KLS-first separator queue stats: metis=%d"
            ", path=%s, sep_rows=%" PRId64
            ", components=%" PRId64
            ", queue=%d/%" PRId64
            ", comp=%" PRId64 "/%" PRId64
            ", rows=%" PRId64 "/%" PRId64
            ", threads=%" PRId64 ", max_thread_rows=%" PRId64
            ", work=%.17g/%.17g"
            ", partitioned=%d/%" PRId64
            ", split_components=%" PRId64
            ", executed=%d/%" PRId64
            ", executed_rows=%" PRId64 "/%" PRId64
            ", private_parallel=%d/%" PRId64
            ", private_parallel_rows=%" PRId64
            ", private_parallel_threads=%" PRId64
            ", pipeline_parallel=%d/%" PRId64
            ", pipeline_parallel_rows=%" PRId64
            ", pipeline_parallel_threads=%" PRId64
            ", pipeline_partial=%d/%" PRId64
            ", pipeline_partial_rows=%" PRId64
            ", pipeline_partial_threads=%" PRId64
            ", pipeline_wait_partial=%d/%" PRId64
            ", pipeline_wait_partial_rows=%" PRId64
            ", pipeline_wait_partial_deps=%" PRId64
            ", pipeline_supernode_update=%d/%" PRId64
            ", pipeline_supernode_update_groups=%" PRId64
            ", pipeline_supernode_update_rows=%" PRId64
            ", pipeline_panel_update=%d/%" PRId64 "/%" PRId64
            ", pipeline_pivot_tail=%d/%" PRId64
            ", pipeline_pivot_tail_rows=%" PRId64
            ", row_cols=%" PRId64 ", btf=%d, scale=%d\n",
            stats.build_has_metis,
            kls_factor_path_name(stats.last_factor_path),
            stats.separator_analyzed_rows,
            stats.separator_component_count,
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
            stats.kls_first_last_separator_queue_split_components,
            stats.kls_first_last_separator_queue_executed,
            stats.kls_first_separator_queue_executed_run_count,
            stats.kls_first_last_separator_queue_executed_private_rows,
            stats.kls_first_last_separator_queue_executed_pipeline_rows,
            stats.kls_first_last_separator_queue_parallel_private,
            stats.kls_first_separator_queue_parallel_private_run_count,
            stats.kls_first_last_separator_queue_parallel_private_rows,
            stats.kls_first_last_separator_queue_parallel_private_threads,
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
            stats.kls_first_last_separator_queue_pipeline_supernode_panel_update_groups,
            stats.kls_first_last_separator_queue_pipeline_supernode_panel_update_rows,
            stats.kls_first_last_separator_queue_pipeline_pivot_tail,
            stats.kls_first_separator_queue_pipeline_pivot_tail_run_count,
            stats.kls_first_last_separator_queue_pipeline_pivot_tail_rows,
            stats.kls_first_last_row_uplooking_columns,
            stats.selected_btf,
            stats.selected_scale);
    ok = 0;
  }
  for (int32_t i = 0; ok && i < n; ++i) {
    if (!close_enough(x[i], expected[i])) {
      fprintf(stderr,
              "unexpected separator queue solution at %d: %.17g != %.17g\n",
              (int)i, x[i], expected[i]);
      ok = 0;
    }
  }
cleanup:
  if (!restore_env_value("KLS_ENABLE_KLS_FIRST_FACTOR", had_saved_first,
                         saved_first)) {
    ok = 0;
  }
  kls_destroy(solver);
  free(ap);
  free(ai);
  free(ax);
  free(b);
  free(x);
  free(expected);
  free(saved_first);
  return ok;
}

static int test_fast_factor_separator_queue_repair(void) {
  const int32_t n = 30000;
  const int32_t nnz = 3 * n - 2;
  int32_t *ap = (int32_t *)malloc(((size_t)n + 1u) * sizeof(*ap));
  int32_t *ai = (int32_t *)malloc((size_t)nnz * sizeof(*ai));
  double *ax0 = (double *)malloc((size_t)nnz * sizeof(*ax0));
  double *ax1 = (double *)malloc((size_t)nnz * sizeof(*ax1));
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

  int32_t p = 0;
  for (int32_t col = 0; col < n; ++col) {
    ap[col] = p;
    expected[col] = 1.0 + (double)(col % 5);
    if (col > 0) {
      ai[p] = col - 1;
      ax0[p] = -1.0;
      ax1[p] = -1.0;
      p++;
    }
    ai[p] = col;
    ax0[p] = 4.0;
    ax1[p] = (col % 97 == 0) ? 1.0e-10 : 4.0;
    p++;
    if (col + 1 < n) {
      ai[p] = col + 1;
      ax0[p] = -1.0;
      ax1[p] = -1.0;
      p++;
    }
  }
  ap[n] = p;
  if (p != nnz) {
    fprintf(stderr, "unexpected separator repair nnz: %d/%d\n", p, nnz);
    free(ap);
    free(ai);
    free(ax0);
    free(ax1);
    free(b);
    free(x);
    free(expected);
    return 0;
  }
  for (int32_t col = 0; col < n; ++col) {
    for (int32_t q = ap[col]; q < ap[col + 1]; ++q) {
      b[ai[q]] += ax1[q] * expected[col];
    }
  }

  const char *saved_first_value = getenv("KLS_ENABLE_KLS_FIRST_FACTOR");
  char *saved_first =
    saved_first_value != NULL ? strdup(saved_first_value) : NULL;
  const int had_saved_first = saved_first_value != NULL;
  const char *saved_row_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_row = saved_row_value != NULL ? strdup(saved_row_value) : NULL;
  const int had_saved_row = saved_row_value != NULL;
  const char *saved_checked_value =
    getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  char *saved_checked =
    saved_checked_value != NULL ? strdup(saved_checked_value) : NULL;
  const int had_saved_checked = saved_checked_value != NULL;

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_METIS;
  options.use_btf = 0;
  options.scale = 0;
  options.pivot_tolerance = 0.1;
  options.static_pivoting = 0;
  options.threads = 4;

  int ok = 1;
  if ((had_saved_first && saved_first == NULL) ||
      (had_saved_row && saved_row == NULL) ||
      (had_saved_checked && saved_checked == NULL)) {
    fprintf(stderr, "failed to save separator repair environment\n");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_KLS_FIRST_FACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_KLS_FIRST_FACTOR");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && !require_ok(kls_create(&solver),
                        "create separator repair")) {
    ok = 0;
  }
  if (ok) {
    const int analyze_status =
      kls_analyze_csc(solver, KLS_INDEX_INT32, n, ap, ai, 0, &options);
    if (analyze_status == KLS_ERR_UNSUPPORTED) {
      goto cleanup;
    }
    if (!require_ok(analyze_status, "analyze separator repair")) {
      ok = 0;
    }
  }
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor separator repair base")) {
    ok = 0;
  }
  if (ok && !require_ok(kls_factor(solver, ax1),
                        "factor separator repair rejected")) {
    ok = 0;
  }
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve separator repair")) {
    ok = 0;
  }

  kls_stats stats;
  memset(&stats, 0, sizeof(stats));
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats separator repair")) {
    ok = 0;
  }
  const int used_separator_tail_scope =
    stats.fast_kls_block_restart_last_row_pipeline_separator_tail_scope == 1;
  const int used_separator_queue =
    stats.fast_kls_block_restart_last_row_pipeline_separator_queue == 1;
  const int separator_tail_scope_ok =
    used_separator_tail_scope &&
    stats.fast_kls_block_restart_last_row_pipeline_etree_tail == 1 &&
    stats.fast_kls_block_restart_row_pipeline_separator_tail_scope_count >= 1 &&
    stats.fast_kls_block_restart_last_row_pipeline_separator_tail_scope_rows >
      0 &&
    stats.fast_kls_block_restart_last_row_pipeline_separator_tail_scope_rows ==
      stats.fast_kls_block_restart_last_row_pipeline_rows &&
    stats.fast_kls_block_restart_last_row_pipeline_rows <= n;
  const int separator_queue_ok =
    used_separator_queue &&
    stats.fast_kls_block_restart_row_pipeline_separator_queue_count >= 1 &&
    stats.fast_kls_block_restart_last_row_pipeline_separator_private_rows > 0 &&
    stats.fast_kls_block_restart_last_row_pipeline_separator_pipeline_rows > 0 &&
    stats.fast_kls_block_restart_last_row_pipeline_separator_private_rows +
      stats.fast_kls_block_restart_last_row_pipeline_separator_pipeline_rows ==
        n &&
    stats.fast_kls_block_restart_last_row_pipeline_separator_private_threads >
      1 &&
    stats.fast_kls_block_restart_last_row_pipeline_separator_partitioned == 1 &&
    stats.fast_kls_block_restart_last_row_pipeline_separator_split_components >
      0 &&
    stats.fast_kls_block_restart_last_row_pipeline_rows == n;
  if (ok && (stats.build_has_metis != 1 ||
             stats.last_factor_path != KLS_FACTOR_PATH_KLS_FAST_REFACTOR ||
             stats.separator_analyzed_rows != n ||
             stats.separator_component_count <= 1 ||
             stats.fast_block_restarts < 1 ||
             stats.fast_kls_block_restart_last_row_pipeline != 1 ||
             (!separator_tail_scope_ok && !separator_queue_ok) ||
             stats.selected_btf != 0 ||
             stats.selected_scale != 0)) {
    fprintf(stderr,
            "unexpected separator repair stats: metis=%d, path=%s"
            ", sep_rows=%" PRId64 ", components=%" PRId64
            ", block_restarts=%d, rowpipe=%d/%" PRId64
            ", tail_scope=%d/%" PRId64 "/%" PRId64
            ", sepq=%d/%" PRId64 ", rows=%" PRId64 "/%" PRId64
            ", threads=%" PRId64 ", partitioned=%d, split=%" PRId64
            ", btf=%d, scale=%d\n",
            stats.build_has_metis,
            kls_factor_path_name(stats.last_factor_path),
            stats.separator_analyzed_rows,
            stats.separator_component_count,
            stats.fast_block_restarts,
            stats.fast_kls_block_restart_last_row_pipeline,
            stats.fast_kls_block_restart_last_row_pipeline_rows,
            stats.fast_kls_block_restart_last_row_pipeline_separator_tail_scope,
            stats.fast_kls_block_restart_row_pipeline_separator_tail_scope_count,
            stats.fast_kls_block_restart_last_row_pipeline_separator_tail_scope_rows,
            stats.fast_kls_block_restart_last_row_pipeline_separator_queue,
            stats
              .fast_kls_block_restart_row_pipeline_separator_queue_count,
            stats
              .fast_kls_block_restart_last_row_pipeline_separator_private_rows,
            stats
              .fast_kls_block_restart_last_row_pipeline_separator_pipeline_rows,
            stats
              .fast_kls_block_restart_last_row_pipeline_separator_private_threads,
            stats
              .fast_kls_block_restart_last_row_pipeline_separator_partitioned,
            stats
              .fast_kls_block_restart_last_row_pipeline_separator_split_components,
            stats.selected_btf,
            stats.selected_scale);
    ok = 0;
  }
  if (ok && !require_pivoting_tail_plan(&stats, "separator repair")) {
    ok = 0;
  }
  if (ok && !require_tail_overcompute_bounds(&stats, "separator repair")) {
    ok = 0;
  }
  for (int32_t i = 0; ok && i < n; ++i) {
    if (!close_enough(x[i], expected[i])) {
      fprintf(stderr,
              "unexpected separator repair solution at %d: %.17g != %.17g\n",
              (int)i, x[i], expected[i]);
      ok = 0;
    }
  }

cleanup:
  if (!restore_env_value("KLS_ENABLE_KLS_FIRST_FACTOR", had_saved_first,
                         saved_first)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_ROW_REFACTOR", had_saved_row,
                         saved_row)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_CHECKED_ROW_REFACTOR",
                         had_saved_checked, saved_checked)) {
    ok = 0;
  }
  kls_destroy(solver);
  free(ap);
  free(ai);
  free(ax0);
  free(ax1);
  free(b);
  free(x);
  free(expected);
  free(saved_first);
  free(saved_row);
  free(saved_checked);
  return ok;
}

static int test_kls_first_separator_pipeline_pivot_epoch(void) {
  const int32_t nx = 180;
  const int32_t ny = 170;
  const int32_t n = nx * ny;
  int64_t nnz64 = 0;
  for (int32_t y = 0; y < ny; ++y) {
    for (int32_t x = 0; x < nx; ++x) {
      nnz64 += 1;
      if (x > 0) nnz64++;
      if (x + 1 < nx) nnz64++;
      if (y > 0) nnz64++;
      if (y + 1 < ny) nnz64++;
    }
  }
  if (nnz64 <= 0 || nnz64 > INT32_MAX) {
    return 0;
  }
  const int32_t nnz = (int32_t)nnz64;
  int32_t *ap = (int32_t *)malloc(((size_t)n + 1u) * sizeof(*ap));
  int32_t *ai = (int32_t *)malloc((size_t)nnz * sizeof(*ai));
  double *ax = (double *)malloc((size_t)nnz * sizeof(*ax));
  double *b = (double *)calloc((size_t)n, sizeof(*b));
  double *xvec = (double *)calloc((size_t)n, sizeof(*xvec));
  double *expected = (double *)malloc((size_t)n * sizeof(*expected));
  double *residual = NULL;
  if (ap == NULL || ai == NULL || ax == NULL || b == NULL ||
      xvec == NULL || expected == NULL) {
    free(ap);
    free(ai);
    free(ax);
    free(b);
    free(xvec);
    free(expected);
    return 0;
  }

  int32_t p = 0;
  for (int32_t y = 0; y < ny; ++y) {
    for (int32_t x = 0; x < nx; ++x) {
      const int32_t col = y * nx + x;
      const int weak_diag =
        abs(x - nx / 2) <= 2 ||
        abs(y - ny / 2) <= 2 ||
        abs(x - nx / 4) <= 1 ||
        abs(x - 3 * nx / 4) <= 1 ||
        abs(y - ny / 4) <= 1 ||
        abs(y - 3 * ny / 4) <= 1;
      ap[col] = p;
      expected[col] = 1.0 + 0.001 * (double)(col % 17);
      if (y > 0) {
        ai[p] = col - nx;
        ax[p] = -1.0;
        p++;
      }
      if (x > 0) {
        ai[p] = col - 1;
        ax[p] = -1.0;
        p++;
      }
      ai[p] = col;
      ax[p] = weak_diag ? 1.0e-8
                        : 5.0 + 1.0e-4 * (double)(col % 13);
      p++;
      if (x + 1 < nx) {
        ai[p] = col + 1;
        ax[p] = -1.0;
        p++;
      }
      if (y + 1 < ny) {
        ai[p] = col + nx;
        ax[p] = -1.0;
        p++;
      }
    }
  }
  ap[n] = p;
  if (p != nnz) {
    fprintf(stderr, "unexpected separator pivot epoch nnz: %d/%d\n",
            p, nnz);
    free(ap);
    free(ai);
    free(ax);
    free(b);
    free(xvec);
    free(expected);
    return 0;
  }
  for (int32_t col = 0; col < n; ++col) {
    for (int32_t q = ap[col]; q < ap[col + 1]; ++q) {
      b[ai[q]] += ax[q] * expected[col];
    }
  }

  const char *saved_first_value = getenv("KLS_ENABLE_KLS_FIRST_FACTOR");
  char *saved_first =
    saved_first_value != NULL ? strdup(saved_first_value) : NULL;
  const int had_saved_first = saved_first_value != NULL;

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_METIS;
  options.use_btf = 0;
  options.scale = 0;
  options.static_pivoting = 0;
  options.threads = 4;

  int ok = 1;
  if (had_saved_first && saved_first == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_KLS_FIRST_FACTOR\n");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_KLS_FIRST_FACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_KLS_FIRST_FACTOR");
    ok = 0;
  }
  /* deterministic pivot-tail coverage: whether the ordering's
     separator happens to intersect the weak cross is RNG-dependent,
     so force every 8th pipeline commit through the real pivot
     exchange + tail restart path (see
     kls_row_first_debug_force_pivot_every) */
  if (ok && setenv("KLS_DEBUG_FORCE_ROW_PIPELINE_PIVOT", "8", 1) != 0) {
    perror("setenv KLS_DEBUG_FORCE_ROW_PIPELINE_PIVOT");
    ok = 0;
  }
  if (ok && !require_ok(kls_create(&solver),
                        "create separator pivot epoch")) {
    ok = 0;
  }
  if (ok) {
    const int analyze_status =
      kls_analyze_csc(solver, KLS_INDEX_INT32, n, ap, ai, 0, &options);
    if (analyze_status == KLS_ERR_UNSUPPORTED) {
      goto cleanup;
    }
    if (!require_ok(analyze_status, "analyze separator pivot epoch")) {
      ok = 0;
    }
  }
  if (ok && !require_ok(kls_factor(solver, ax),
                        "factor separator pivot epoch")) {
    ok = 0;
  }

  kls_stats factor_stats;
  memset(&factor_stats, 0, sizeof(factor_stats));
  factor_stats.struct_size = sizeof(factor_stats);
  if (ok && !require_ok(kls_get_stats(solver, &factor_stats),
                        "stats separator pivot epoch")) {
    ok = 0;
  }
  if (ok && (factor_stats.build_has_metis != 1 ||
             factor_stats.last_factor_path != KLS_FACTOR_PATH_KLS_FIRST ||
             factor_stats.separator_analyzed_rows != n ||
             factor_stats.separator_component_count < 3 ||
             factor_stats.separator_pipeline_max_rows <= 1 ||
             factor_stats.kls_first_last_separator_queue != 1 ||
             factor_stats.kls_first_last_separator_queue_pipeline_rows <= 0 ||
             factor_stats.kls_first_last_separator_queue_parallel_pipeline !=
               1 ||
             factor_stats.kls_first_last_separator_dynamic_column_pivots <=
               0 ||
             factor_stats.kls_first_last_separator_queue_pipeline_pivot_tail !=
               1 ||
             factor_stats.kls_first_last_separator_queue_pipeline_pivot_restarts <=
               0 ||
             factor_stats.kls_first_last_separator_queue_pipeline_pivot_tail_rows <=
               0 ||
             factor_stats.kls_first_last_separator_queue_pipeline_prefix_panel_rebuild !=
               0 ||
             factor_stats.kls_first_separator_queue_pipeline_prefix_panel_rebuild_count !=
               0 ||
             factor_stats.kls_first_last_separator_queue_pipeline_prefix_panel_rebuild_rows !=
               0 ||
             factor_stats.kls_first_row_panel_cache_build_count <= 0 ||
             factor_stats.kls_first_row_panel_cache_build_panels <= 0 ||
             factor_stats.kls_first_row_panel_cache_build_entries <= 0 ||
             factor_stats.kls_first_last_separator_queue_pipeline_supernode_panel_update !=
               1 ||
             factor_stats.kls_first_last_separator_queue_pipeline_supernode_panel_update_groups <=
               0 ||
             factor_stats.kls_first_last_separator_queue_pipeline_supernode_panel_update_rows <=
               0 ||
             factor_stats.kls_first_last_separator_queue_pipeline_pivot_serial_rows * 2 >
               factor_stats.kls_first_last_separator_queue_pipeline_rows)) {
    fprintf(stderr,
            "unexpected separator pivot epoch stats: metis=%d"
            ", path=%s, sep_rows=%" PRId64
            ", components=%" PRId64
            ", pipe_max=%" PRId64
            ", queue=%d, pipe=%" PRId64
            ", pipe_parallel=%d"
            ", sep_pivots=%" PRId64
            ", prefix_panel_rebuild=%d/%" PRId64 "/%" PRId64
            ", cache_build=%" PRId64 "/%" PRId64 "/%" PRId64
            ", panel_update=%d/%" PRId64 "/%" PRId64
            ", pivot_tail=%d, restarts=%" PRId64
            ", tail_rows=%" PRId64
            ", serial=%" PRId64 "\n",
            factor_stats.build_has_metis,
            kls_factor_path_name(factor_stats.last_factor_path),
            factor_stats.separator_analyzed_rows,
            factor_stats.separator_component_count,
            factor_stats.separator_pipeline_max_rows,
            factor_stats.kls_first_last_separator_queue,
            factor_stats.kls_first_last_separator_queue_pipeline_rows,
            factor_stats.kls_first_last_separator_queue_parallel_pipeline,
            factor_stats.kls_first_last_separator_dynamic_column_pivots,
            factor_stats.kls_first_last_separator_queue_pipeline_prefix_panel_rebuild,
            factor_stats.kls_first_separator_queue_pipeline_prefix_panel_rebuild_count,
            factor_stats.kls_first_last_separator_queue_pipeline_prefix_panel_rebuild_rows,
            factor_stats.kls_first_row_panel_cache_build_count,
            factor_stats.kls_first_row_panel_cache_build_panels,
            factor_stats.kls_first_row_panel_cache_build_entries,
            factor_stats.kls_first_last_separator_queue_pipeline_supernode_panel_update,
            factor_stats.kls_first_last_separator_queue_pipeline_supernode_panel_update_groups,
            factor_stats.kls_first_last_separator_queue_pipeline_supernode_panel_update_rows,
            factor_stats.kls_first_last_separator_queue_pipeline_pivot_tail,
            factor_stats.kls_first_last_separator_queue_pipeline_pivot_restarts,
            factor_stats.kls_first_last_separator_queue_pipeline_pivot_tail_rows,
            factor_stats.kls_first_last_separator_queue_pipeline_pivot_serial_rows);
    ok = 0;
  }
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, xvec, 0),
                        "solve separator pivot epoch")) {
    ok = 0;
  }

  double max_solution_error = 0.0;
  double max_rhs = 0.0;
  for (int32_t row = 0; row < n; ++row) {
    const double err = fabs(xvec[row] - expected[row]);
    if (err > max_solution_error) {
      max_solution_error = err;
    }
    if (fabs(b[row]) > max_rhs) {
      max_rhs = fabs(b[row]);
    }
  }
  residual = (double *)malloc((size_t)n * sizeof(*residual));
  if (residual == NULL) {
    ok = 0;
  } else {
    memcpy(residual, b, (size_t)n * sizeof(*residual));
    for (int32_t col = 0; col < n; ++col) {
      for (int32_t q = ap[col]; q < ap[col + 1]; ++q) {
        residual[ai[q]] -= ax[q] * xvec[col];
      }
    }
    double max_residual = 0.0;
    for (int32_t row = 0; row < n; ++row) {
      const double err = fabs(residual[row]);
      if (err > max_residual) {
        max_residual = err;
      }
    }
    if (ok && (max_solution_error > 1.0e-7 ||
               max_residual / (1.0 + max_rhs) > 1.0e-7)) {
      fprintf(stderr,
              "unexpected separator pivot epoch residual:"
              " max_solution_error=%.17g max_residual=%.17g"
              " max_rhs=%.17g\n",
              max_solution_error, max_residual, max_rhs);
      ok = 0;
    }
  }

cleanup:
  unsetenv("KLS_DEBUG_FORCE_ROW_PIPELINE_PIVOT");
  if (!restore_env_value("KLS_ENABLE_KLS_FIRST_FACTOR", had_saved_first,
                         saved_first)) {
    ok = 0;
  }
  kls_destroy(solver);
  free(ap);
  free(ai);
  free(ax);
  free(b);
  free(xvec);
  free(expected);
  free(residual);
  free(saved_first);
  return ok;
}

static int test_auto_kls_first_skips_scaled_single_block(void) {
  const int32_t n = 150000;
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
    const double diag = 2.0 + (double)(col % 17);
    ai[col] = col;
    ax[col] = diag;
    expected[col] = 1.0 + (double)(col % 19);
    b[col] = diag * expected[col];
  }

  const char *saved_first_value = getenv("KLS_ENABLE_KLS_FIRST_FACTOR");
  char *saved_first =
    saved_first_value != NULL ? strdup(saved_first_value) : NULL;
  const int had_saved_first = saved_first_value != NULL;

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_NATURAL;
  options.use_btf = 0;
  options.scale = 1;
  options.static_pivoting = 0;

  int ok = 1;
  if (had_saved_first && saved_first == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_KLS_FIRST_FACTOR\n");
    ok = 0;
  }
  if (ok && unsetenv("KLS_ENABLE_KLS_FIRST_FACTOR") != 0) {
    perror("unsetenv KLS_ENABLE_KLS_FIRST_FACTOR");
    ok = 0;
  }
  if (ok && !require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, n, ap, ai, 0,
                                        &options),
                        "analyze scaled single-block KLS first skip")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax),
                        "factor scaled single-block KLS first skip")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve scaled single-block KLS first skip")) ok = 0;

  kls_stats stats;
  memset(&stats, 0, sizeof(stats));
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats scaled single-block KLS first skip")) ok = 0;
  if (ok && (stats.last_factor_path != KLS_FACTOR_PATH_KLU_FIRST ||
             stats.kls_first_auto_skipped_scaled_single_block != 1 ||
             stats.kls_first_auto_skipped_scaled_single_block_count < 1 ||
             stats.kls_first_last_row_uplooking_columns != 0 ||
             stats.selected_btf != 0 || stats.selected_scale != 1)) {
    fprintf(stderr,
            "unexpected scaled single-block KLS-first skip stats: path=%s"
            ", skip=%d/%" PRId64 ", row_cols=%" PRId64
            ", btf=%d, scale=%d\n",
            kls_factor_path_name(stats.last_factor_path),
            stats.kls_first_auto_skipped_scaled_single_block,
            stats.kls_first_auto_skipped_scaled_single_block_count,
            stats.kls_first_last_row_uplooking_columns,
            stats.selected_btf, stats.selected_scale);
    ok = 0;
  }
  for (int32_t i = 0; ok && i < n; ++i) {
    if (!close_enough(x[i], expected[i])) {
      fprintf(stderr,
              "unexpected scaled single-block solution at %d:"
              " %.17g != %.17g\n",
              (int)i, x[i], expected[i]);
      ok = 0;
    }
  }

  if (!restore_env_value("KLS_ENABLE_KLS_FIRST_FACTOR", had_saved_first,
                         saved_first)) {
    ok = 0;
  }
  kls_destroy(solver);
  free(ap);
  free(ai);
  free(ax);
  free(b);
  free(x);
  free(expected);
  free(saved_first);
  return ok;
}

static int test_forced_kls_first_scaled_single_block_separator_queue(void) {
  const int32_t n = 150000;
  const int32_t nnz = 3 * n - 2;
  int32_t *ap = (int32_t *)malloc(((size_t)n + 1u) * sizeof(*ap));
  int32_t *ai = (int32_t *)malloc((size_t)nnz * sizeof(*ai));
  double *ax = (double *)malloc((size_t)nnz * sizeof(*ax));
  double *b = (double *)calloc((size_t)n, sizeof(*b));
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

  int32_t p = 0;
  for (int32_t col = 0; col < n; ++col) {
    ap[col] = p;
    expected[col] = 1.0 + (double)(col % 7);
    if (col > 0) {
      ai[p] = col - 1;
      ax[p] = -0.25;
      p++;
    }
    ai[p] = col;
    ax[p] = 5.0 + (double)(col % 3);
    p++;
    if (col + 1 < n) {
      ai[p] = col + 1;
      ax[p] = -0.25;
      p++;
    }
  }
  ap[n] = p;
  if (p != nnz) {
    fprintf(stderr, "unexpected scaled separator nnz: %d/%d\n", p, nnz);
    free(ap);
    free(ai);
    free(ax);
    free(b);
    free(x);
    free(expected);
    return 0;
  }
  for (int32_t col = 0; col < n; ++col) {
    for (int32_t q = ap[col]; q < ap[col + 1]; ++q) {
      b[ai[q]] += ax[q] * expected[col];
    }
  }

  const char *saved_first_value = getenv("KLS_ENABLE_KLS_FIRST_FACTOR");
  char *saved_first =
    saved_first_value != NULL ? strdup(saved_first_value) : NULL;
  const int had_saved_first = saved_first_value != NULL;

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_METIS;
  options.use_btf = 0;
  options.scale = 1;
  options.static_pivoting = 0;
  options.threads = 4;

  int ok = 1;
  if (had_saved_first && saved_first == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_KLS_FIRST_FACTOR\n");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_KLS_FIRST_FACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_KLS_FIRST_FACTOR");
    ok = 0;
  }
  if (ok && !require_ok(kls_create(&solver),
                        "create scaled separator KLS first")) {
    ok = 0;
  }
  if (ok) {
    const int analyze_status =
      kls_analyze_csc(solver, KLS_INDEX_INT32, n, ap, ai, 0, &options);
    if (analyze_status == KLS_ERR_UNSUPPORTED) {
      goto cleanup;
    }
    if (!require_ok(analyze_status,
                    "analyze scaled separator KLS first")) {
      ok = 0;
    }
  }
  if (ok && !require_ok(kls_factor(solver, ax),
                        "factor scaled separator KLS first")) {
    ok = 0;
  }
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve scaled separator KLS first")) {
    ok = 0;
  }

  kls_stats stats;
  memset(&stats, 0, sizeof(stats));
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats scaled separator KLS first")) {
    ok = 0;
  }
  if (ok && (stats.build_has_metis != 1 ||
             stats.last_factor_path != KLS_FACTOR_PATH_KLS_FIRST ||
             stats.kls_first_auto_skipped_scaled_single_block != 0 ||
             stats.kls_first_auto_skipped_scaled_single_block_count != 0 ||
             stats.separator_analyzed_rows != n ||
             stats.separator_private_components <= 0 ||
             stats.separator_pipeline_components <= 0 ||
             stats.separator_private_rows <= 0 ||
             stats.separator_pipeline_rows <= 0 ||
             stats.kls_first_last_separator_queue_executed != 1 ||
             stats.kls_first_last_separator_queue_parallel_private != 1 ||
             stats.kls_first_last_separator_queue_parallel_pipeline != 1 ||
             stats.kls_first_last_row_uplooking_columns != n ||
             stats.kls_first_last_row_refactor_seeded_rows != n ||
             stats.selected_btf != 0 || stats.selected_scale != 1)) {
    fprintf(stderr,
            "unexpected scaled separator KLS-first stats: metis=%d"
            ", path=%s, skip=%d/%" PRId64
            ", sep_rows=%" PRId64 ", sep_comp=%" PRId64 "/%" PRId64
            ", sep_rows_private_pipeline=%" PRId64 "/%" PRId64
            ", queue_exec=%d, private_parallel=%d, pipeline_parallel=%d"
            ", row_cols=%" PRId64 ", row_seed=%" PRId64
            ", btf=%d, scale=%d\n",
            stats.build_has_metis,
            kls_factor_path_name(stats.last_factor_path),
            stats.kls_first_auto_skipped_scaled_single_block,
            stats.kls_first_auto_skipped_scaled_single_block_count,
            stats.separator_analyzed_rows,
            stats.separator_private_components,
            stats.separator_pipeline_components,
            stats.separator_private_rows,
            stats.separator_pipeline_rows,
            stats.kls_first_last_separator_queue_executed,
            stats.kls_first_last_separator_queue_parallel_private,
            stats.kls_first_last_separator_queue_parallel_pipeline,
            stats.kls_first_last_row_uplooking_columns,
            stats.kls_first_last_row_refactor_seeded_rows,
            stats.selected_btf, stats.selected_scale);
    ok = 0;
  }
  for (int32_t i = 0; ok && i < n; ++i) {
    if (!close_enough(x[i], expected[i])) {
      fprintf(stderr,
              "unexpected scaled separator solution at %d:"
              " %.17g != %.17g\n",
              (int)i, x[i], expected[i]);
      ok = 0;
    }
  }

cleanup:
  if (!restore_env_value("KLS_ENABLE_KLS_FIRST_FACTOR", had_saved_first,
                         saved_first)) {
    ok = 0;
  }
  kls_destroy(solver);
  free(ap);
  free(ai);
  free(ax);
  free(b);
  free(x);
  free(expected);
  free(saved_first);
  return ok;
}

static int test_pre_static_replays_kls_first_factor(void) {
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

  const char *saved_first_value = getenv("KLS_ENABLE_KLS_FIRST_FACTOR");
  char *saved_first =
    saved_first_value != NULL ? strdup(saved_first_value) : NULL;
  const int had_saved_first = saved_first_value != NULL;

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_AUTO;

  int ok = 1;
  if (had_saved_first && saved_first == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_KLS_FIRST_FACTOR\n");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_KLS_FIRST_FACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_KLS_FIRST_FACTOR");
    ok = 0;
  }
  if (ok && !require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, n, ap, ai, 0,
                                        &options),
                        "analyze pre-static KLS replay")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax),
                        "factor pre-static KLS replay")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve pre-static KLS replay")) ok = 0;

  kls_stats stats;
  memset(&stats, 0, sizeof(stats));
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats pre-static KLS replay")) ok = 0;
  if (ok && (!stats.selected_static_pivoting ||
             !stats.selected_exact_matching ||
             stats.last_factor_path != KLS_FACTOR_PATH_KLS_FIRST ||
             stats.kls_first_last_row_uplooking_columns != n ||
             stats.kls_first_row_uplooking_column_count < n ||
             stats.kls_first_last_row_refactor_seeded_rows != n ||
             stats.kls_first_row_refactor_seeded_row_count < n)) {
    fprintf(stderr,
            "unexpected pre-static replay stats: static=%d exact=%d path=%s"
            ", row_cols=%" PRId64 "/%" PRId64
            ", row_seed=%" PRId64 "/%" PRId64 "\n",
            stats.selected_static_pivoting,
            stats.selected_exact_matching,
            kls_factor_path_name(stats.last_factor_path),
            stats.kls_first_last_row_uplooking_columns,
            stats.kls_first_row_uplooking_column_count,
            stats.kls_first_last_row_refactor_seeded_rows,
            stats.kls_first_row_refactor_seeded_row_count);
    ok = 0;
  }
  for (int32_t i = 0; ok && i < n; ++i) {
    if (!close_enough(x[i], expected[i])) {
      fprintf(stderr,
              "unexpected pre-static replay solution at %d: %.17g != %.17g\n",
              (int)i, x[i], expected[i]);
      ok = 0;
    }
  }

  if (!restore_env_value("KLS_ENABLE_KLS_FIRST_FACTOR", had_saved_first,
                         saved_first)) {
    ok = 0;
  }
  kls_destroy(solver);
  free(ap);
  free(ai);
  free(ax);
  free(b);
  free(x);
  free(expected);
  free(saved_first);
  return ok;
}

static int test_experimental_row_uplooking_first_factor(void) {
  const int32_t ap[] = {0, 3, 6, 9};
  const int32_t ai[] = {0, 1, 2, 0, 1, 2, 0, 1, 2};
  const double ax[] = {
    5.0, 0.25, 0.125,
    0.5, 6.0, 0.375,
    0.25, 0.5, 7.0
  };
  const double b[] = {6.75, 13.75, 21.875};
  double x[3] = {0.0, 0.0, 0.0};

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_NATURAL;
  options.use_btf = 0;
  options.scale = 0;
  options.static_pivoting = 0;
  options.pivot_tolerance = 0.001;

  const char *saved_first_value = getenv("KLS_ENABLE_KLS_FIRST_FACTOR");
  char *saved_first =
    saved_first_value != NULL ? strdup(saved_first_value) : NULL;
  const int had_saved_first = saved_first_value != NULL;
  const char *saved_row_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_row = saved_row_value != NULL ? strdup(saved_row_value) : NULL;
  const int had_saved_row = saved_row_value != NULL;
  const char *saved_checked_value =
    getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  char *saved_checked =
    saved_checked_value != NULL ? strdup(saved_checked_value) : NULL;
  const int had_saved_checked = saved_checked_value != NULL;

  int ok = 1;
  if (had_saved_first && saved_first == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_KLS_FIRST_FACTOR\n");
    ok = 0;
  }
  if (had_saved_row && saved_row == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_ROW_REFACTOR\n");
    ok = 0;
  }
  if (had_saved_checked && saved_checked == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_CHECKED_ROW_REFACTOR\n");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_KLS_FIRST_FACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_KLS_FIRST_FACTOR");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR=0");
    ok = 0;
  }

  if (ok && !require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 3, ap, ai, 0,
                                        &options),
                        "analyze row-up-looking first factor")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax),
                        "factor row-up-looking first factor")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve row-up-looking first factor")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats row-up-looking first factor")) {
    ok = 0;
  }
  if (ok && (stats.last_factor_path != KLS_FACTOR_PATH_KLS_FIRST ||
             stats.kls_first_last_row_uplooking_columns != 3 ||
             stats.kls_first_row_uplooking_column_count < 3 ||
             stats.kls_first_last_row_refactor_seeded_rows != 3 ||
             stats.kls_first_row_refactor_seeded_row_count < 3 ||
             stats.kls_first_last_row_supernode_update != 1 ||
             stats.kls_first_row_supernode_update_run_count < 1 ||
             stats.kls_first_last_row_supernode_update_groups < 1 ||
             stats.kls_first_last_row_supernode_update_rows < 2 ||
             stats.kls_first_last_dynamic_column_pivots != 0 ||
             stats.kls_tail_last_mapped_columns != 0 ||
             stats.selected_scale != 0 ||
             stats.selected_btf != 0)) {
    fprintf(stderr,
            "unexpected row-up-looking first-factor stats: path=%s"
            ", row_cols=%" PRId64 "/%" PRId64
            ", row_seed=%" PRId64 "/%" PRId64
            ", row_supernode=%d/%" PRId64 "/%" PRId64 "/%" PRId64
            ", dyn_pivots=%" PRId64
            ", mapped_tail=%" PRId64 ", scale=%d, btf=%d\n",
            kls_factor_path_name(stats.last_factor_path),
            stats.kls_first_last_row_uplooking_columns,
            stats.kls_first_row_uplooking_column_count,
            stats.kls_first_last_row_refactor_seeded_rows,
            stats.kls_first_row_refactor_seeded_row_count,
            stats.kls_first_last_row_supernode_update,
            stats.kls_first_row_supernode_update_run_count,
            stats.kls_first_last_row_supernode_update_groups,
            stats.kls_first_last_row_supernode_update_rows,
            stats.kls_first_last_dynamic_column_pivots,
            stats.kls_tail_last_mapped_columns,
            stats.selected_scale,
            stats.selected_btf);
    ok = 0;
  }
  if (ok && (!close_enough(x[0], 1.0) || !close_enough(x[1], 2.0) ||
             !close_enough(x[2], 3.0))) {
    fprintf(stderr,
            "unexpected row-up-looking first-factor solution:"
            " %.17g %.17g %.17g\n",
            x[0], x[1], x[2]);
    ok = 0;
  }

  if (had_saved_first && saved_first != NULL) {
    if (setenv("KLS_ENABLE_KLS_FIRST_FACTOR", saved_first, 1) != 0) {
      perror("restore KLS_ENABLE_KLS_FIRST_FACTOR");
      ok = 0;
    }
  } else if (!had_saved_first) {
    if (unsetenv("KLS_ENABLE_KLS_FIRST_FACTOR") != 0) {
      perror("unsetenv KLS_ENABLE_KLS_FIRST_FACTOR");
      ok = 0;
    }
  }
  if (had_saved_row && saved_row != NULL) {
    if (setenv("KLS_ENABLE_ROW_REFACTOR", saved_row, 1) != 0) {
      perror("restore KLS_ENABLE_ROW_REFACTOR");
      ok = 0;
    }
  } else if (!had_saved_row) {
    if (unsetenv("KLS_ENABLE_ROW_REFACTOR") != 0) {
      perror("unsetenv KLS_ENABLE_ROW_REFACTOR");
      ok = 0;
    }
  }
  if (had_saved_checked && saved_checked != NULL) {
    if (setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", saved_checked, 1) != 0) {
      perror("restore KLS_ENABLE_CHECKED_ROW_REFACTOR");
      ok = 0;
    }
  } else if (!had_saved_checked) {
    if (unsetenv("KLS_ENABLE_CHECKED_ROW_REFACTOR") != 0) {
      perror("unsetenv KLS_ENABLE_CHECKED_ROW_REFACTOR");
      ok = 0;
    }
  }

  kls_destroy(solver);
  free(saved_first);
  free(saved_row);
  free(saved_checked);
  return ok;
}

static int test_experimental_row_uplooking_natural_pipeline(void) {
  const int32_t n = 64;
  const int32_t nnz = 3 * n - 2;
  int32_t *ap = (int32_t *)malloc(((size_t)n + 1u) * sizeof(*ap));
  int32_t *ai = (int32_t *)malloc((size_t)nnz * sizeof(*ai));
  double *ax = (double *)malloc((size_t)nnz * sizeof(*ax));
  double *b = (double *)calloc((size_t)n, sizeof(*b));
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

  int32_t p = 0;
  for (int32_t col = 0; col < n; ++col) {
    ap[col] = p;
    expected[col] = 1.0 + 0.125 * (double)(col % 7);
    if (col > 0) {
      ai[p] = col - 1;
      ax[p] = -0.25;
      p++;
    }
    ai[p] = col;
    ax[p] = 4.0 + 0.01 * (double)(col % 11);
    p++;
    if (col + 1 < n) {
      ai[p] = col + 1;
      ax[p] = -0.25;
      p++;
    }
  }
  ap[n] = p;
  if (p != nnz) {
    fprintf(stderr, "unexpected natural row-pipeline nnz: %d/%d\n",
            p, nnz);
    free(ap);
    free(ai);
    free(ax);
    free(b);
    free(x);
    free(expected);
    return 0;
  }
  for (int32_t col = 0; col < n; ++col) {
    for (int32_t q = ap[col]; q < ap[col + 1]; ++q) {
      b[ai[q]] += ax[q] * expected[col];
    }
  }

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_NATURAL;
  options.use_btf = 0;
  options.scale = 0;
  options.static_pivoting = 0;
  options.pivot_tolerance = 0.001;
  options.threads = 3;

  const char *saved_first_value = getenv("KLS_ENABLE_KLS_FIRST_FACTOR");
  char *saved_first =
    saved_first_value != NULL ? strdup(saved_first_value) : NULL;
  const int had_saved_first = saved_first_value != NULL;
  const char *saved_row_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_row = saved_row_value != NULL ? strdup(saved_row_value) : NULL;
  const int had_saved_row = saved_row_value != NULL;
  const char *saved_checked_value =
    getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  char *saved_checked =
    saved_checked_value != NULL ? strdup(saved_checked_value) : NULL;
  const int had_saved_checked = saved_checked_value != NULL;

  int ok = 1;
  if (had_saved_first && saved_first == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_KLS_FIRST_FACTOR\n");
    ok = 0;
  }
  if (had_saved_row && saved_row == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_ROW_REFACTOR\n");
    ok = 0;
  }
  if (had_saved_checked && saved_checked == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_CHECKED_ROW_REFACTOR\n");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_KLS_FIRST_FACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_KLS_FIRST_FACTOR");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR=0");
    ok = 0;
  }

  if (ok && !require_ok(kls_create(&solver),
                        "create natural row pipeline")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(
                          solver, KLS_INDEX_INT32, n, ap, ai, 0, &options),
                        "analyze natural row pipeline")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax),
                        "factor natural row pipeline")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve natural row pipeline")) ok = 0;

  kls_stats stats;
  memset(&stats, 0, sizeof(stats));
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats natural row pipeline")) {
    ok = 0;
  }
  if (ok && (stats.last_factor_path != KLS_FACTOR_PATH_KLS_FIRST ||
             stats.kls_first_last_row_uplooking_columns != n ||
             stats.kls_first_last_row_pipeline != 1 ||
             stats.kls_first_row_pipeline_run_count < 1 ||
             stats.kls_first_last_row_pipeline_rows != n ||
             stats.kls_first_last_row_pipeline_threads <= 0 ||
             stats.kls_first_last_row_pipeline_pivot_tail != 0 ||
             stats.kls_first_last_row_pipeline_prefix_panel_rebuild != 0 ||
             stats.kls_first_last_separator_queue_parallel_pipeline != 0 ||
             stats.kls_first_last_separator_queue_executed != 0 ||
             stats.selected_btf != 0 ||
             stats.selected_scale != 0)) {
    fprintf(stderr,
            "unexpected natural row-pipeline stats: path=%s"
            ", row_cols=%" PRId64
            ", row_pipe=%d/%" PRId64 "/%" PRId64 "/%" PRId64
            ", row_pivot_tail=%d"
            ", row_prefix_rebuild=%d"
            ", sep_pipe=%d, sep_exec=%d, btf=%d, scale=%d\n",
            kls_factor_path_name(stats.last_factor_path),
            stats.kls_first_last_row_uplooking_columns,
            stats.kls_first_last_row_pipeline,
            stats.kls_first_row_pipeline_run_count,
            stats.kls_first_last_row_pipeline_rows,
            stats.kls_first_last_row_pipeline_threads,
            stats.kls_first_last_row_pipeline_pivot_tail,
            stats.kls_first_last_row_pipeline_prefix_panel_rebuild,
            stats.kls_first_last_separator_queue_parallel_pipeline,
            stats.kls_first_last_separator_queue_executed,
            stats.selected_btf,
            stats.selected_scale);
    ok = 0;
  }
  for (int32_t i = 0; ok && i < n; ++i) {
    if (!close_enough(x[i], expected[i])) {
      fprintf(stderr,
              "unexpected natural row-pipeline solution at %d:"
              " %.17g != %.17g\n",
              (int)i, x[i], expected[i]);
      ok = 0;
    }
  }

  if (!restore_env_value("KLS_ENABLE_KLS_FIRST_FACTOR", had_saved_first,
                         saved_first)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_ROW_REFACTOR", had_saved_row,
                         saved_row)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_CHECKED_ROW_REFACTOR",
                         had_saved_checked, saved_checked)) {
    ok = 0;
  }
  kls_destroy(solver);
  free(ap);
  free(ai);
  free(ax);
  free(b);
  free(x);
  free(expected);
  free(saved_first);
  free(saved_row);
  free(saved_checked);
  return ok;
}

static int test_experimental_row_uplooking_first_consumer_panel(void) {
  const int32_t n = 40;
  const int32_t nnz = n * n;
  int32_t *ap = (int32_t *)calloc((size_t)n + 1u, sizeof(*ap));
  int32_t *ai = (int32_t *)malloc((size_t)nnz * sizeof(*ai));
  double *ax = (double *)malloc((size_t)nnz * sizeof(*ax));
  double *b = (double *)calloc((size_t)n, sizeof(*b));
  double *x = (double *)calloc((size_t)n, sizeof(*x));
  double *expected = (double *)malloc((size_t)n * sizeof(*expected));
  if (ap == NULL || ai == NULL || ax == NULL || b == NULL || x == NULL ||
      expected == NULL) {
    free(ap);
    free(ai);
    free(ax);
    free(b);
    free(x);
    free(expected);
    return 0;
  }
  for (int32_t row = 0; row < n; ++row) {
    expected[row] = 1.0 + (double)(row % 7);
  }
  int32_t pos = 0;
  for (int32_t col = 0; col < n; ++col) {
    ap[col] = pos;
    for (int32_t row = 0; row < n; ++row) {
      const double value = row == col
        ? 64.0 + 0.25 * (double)col
        : 0.001 * (double)(1 + ((row + 3 * col) % 11));
      ai[pos] = row;
      ax[pos] = value;
      b[row] += value * expected[col];
      pos++;
    }
  }
  ap[n] = pos;

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_NATURAL;
  options.use_btf = 0;
  options.scale = 0;
  options.static_pivoting = 0;
  options.pivot_tolerance = 0.001;

  const char *saved_first_value = getenv("KLS_ENABLE_KLS_FIRST_FACTOR");
  char *saved_first =
    saved_first_value != NULL ? strdup(saved_first_value) : NULL;
  const int had_saved_first = saved_first_value != NULL;
  const char *saved_row_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_row = saved_row_value != NULL ? strdup(saved_row_value) : NULL;
  const int had_saved_row = saved_row_value != NULL;
  const char *saved_checked_value =
    getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  char *saved_checked =
    saved_checked_value != NULL ? strdup(saved_checked_value) : NULL;
  const int had_saved_checked = saved_checked_value != NULL;
  const char *saved_cblas_value = getenv("KLS_ENABLE_CBLAS_SUPERNODE");
  char *saved_cblas =
    saved_cblas_value != NULL ? strdup(saved_cblas_value) : NULL;
  const int had_saved_cblas = saved_cblas_value != NULL;

  int ok = 1;
  if (had_saved_first && saved_first == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_KLS_FIRST_FACTOR\n");
    ok = 0;
  }
  if (had_saved_row && saved_row == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_ROW_REFACTOR\n");
    ok = 0;
  }
  if (had_saved_checked && saved_checked == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_CHECKED_ROW_REFACTOR\n");
    ok = 0;
  }
  if (had_saved_cblas && saved_cblas == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_CBLAS_SUPERNODE\n");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_KLS_FIRST_FACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_KLS_FIRST_FACTOR");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CBLAS_SUPERNODE", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CBLAS_SUPERNODE=0");
    ok = 0;
  }

  if (ok && !require_ok(kls_create(&solver),
                        "create first-consumer panel")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(
                          solver, KLS_INDEX_INT32, n, ap, ai, 0, &options),
                        "analyze first-consumer panel")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax),
                        "factor first-consumer panel")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve first-consumer panel")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats first-consumer panel")) {
    ok = 0;
  }
  if (ok && (stats.last_factor_path != KLS_FACTOR_PATH_KLS_FIRST ||
             stats.kls_first_last_row_uplooking_columns != n ||
             stats.kls_first_last_row_supernode_update != 1 ||
             stats.kls_first_last_row_supernode_update_groups <= 0 ||
             stats.kls_first_last_row_supernode_update_rows < 16 ||
             stats.kls_first_last_row_supernode_panel_update != 1 ||
             stats.kls_first_row_supernode_panel_update_run_count != 1 ||
             stats.kls_first_last_row_supernode_panel_update_groups <= 0 ||
             stats.kls_first_last_row_supernode_panel_update_rows < 16 ||
             stats.kls_first_row_panel_cache_append_count <= 0 ||
             stats.kls_first_row_panel_cache_append_panels <= 0 ||
             stats.kls_first_row_panel_cache_append_entries <= 0)) {
    fprintf(stderr,
            "unexpected first-consumer panel stats: path=%s"
            ", row_cols=%" PRId64 ", supernode=%d/%" PRId64 "/%" PRId64
            ", panel=%d/%" PRId64 "/%" PRId64 "/%" PRId64
            ", cache_append=%" PRId64 "/%" PRId64 "/%" PRId64 "\n",
            kls_factor_path_name(stats.last_factor_path),
            stats.kls_first_last_row_uplooking_columns,
            stats.kls_first_last_row_supernode_update,
            stats.kls_first_last_row_supernode_update_groups,
            stats.kls_first_last_row_supernode_update_rows,
            stats.kls_first_last_row_supernode_panel_update,
            stats.kls_first_row_supernode_panel_update_run_count,
            stats.kls_first_last_row_supernode_panel_update_groups,
            stats.kls_first_last_row_supernode_panel_update_rows,
            stats.kls_first_row_panel_cache_append_count,
            stats.kls_first_row_panel_cache_append_panels,
            stats.kls_first_row_panel_cache_append_entries);
    ok = 0;
  }
  for (int32_t row = 0; ok && row < n; ++row) {
    if (!close_enough(x[row], expected[row])) {
      fprintf(stderr,
              "unexpected first-consumer panel solution at %d:"
              " %.17g != %.17g\n",
              (int)row, x[row], expected[row]);
      ok = 0;
    }
  }

  if (!restore_env_value("KLS_ENABLE_KLS_FIRST_FACTOR", had_saved_first,
                         saved_first)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_ROW_REFACTOR", had_saved_row,
                         saved_row)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_CHECKED_ROW_REFACTOR",
                         had_saved_checked, saved_checked)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_CBLAS_SUPERNODE", had_saved_cblas,
                         saved_cblas)) {
    ok = 0;
  }
  kls_destroy(solver);
  free(ap);
  free(ai);
  free(ax);
  free(b);
  free(x);
  free(expected);
  free(saved_first);
  free(saved_row);
  free(saved_checked);
  free(saved_cblas);
  return ok;
}

static int test_experimental_row_uplooking_lazy_panel_prefix(void) {
  const int32_t ap[] = {0, 4, 8, 12, 16};
  const int32_t ai[] = {
    0, 1, 2, 3,
    0, 1, 2, 3,
    0, 1, 2, 3,
    0, 1, 2, 3
  };
  const double ax[] = {
    5.0, 0.25, 0.125, 0.0625,
    0.5, 6.0, 0.375, 0.1875,
    0.25, 0.5, 7.0, 0.3125,
    0.125, 0.25, 0.5, 8.0
  };
  const double b[] = {7.25, 14.75, 23.875, 33.375};
  double x[4] = {0.0, 0.0, 0.0, 0.0};

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_NATURAL;
  options.use_btf = 0;
  options.scale = 0;
  options.static_pivoting = 0;
  options.pivot_tolerance = 0.001;

  const char *saved_first_value = getenv("KLS_ENABLE_KLS_FIRST_FACTOR");
  char *saved_first =
    saved_first_value != NULL ? strdup(saved_first_value) : NULL;
  const int had_saved_first = saved_first_value != NULL;
  const char *saved_row_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_row = saved_row_value != NULL ? strdup(saved_row_value) : NULL;
  const int had_saved_row = saved_row_value != NULL;
  const char *saved_checked_value =
    getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  char *saved_checked =
    saved_checked_value != NULL ? strdup(saved_checked_value) : NULL;
  const int had_saved_checked = saved_checked_value != NULL;

  int ok = 1;
  if (had_saved_first && saved_first == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_KLS_FIRST_FACTOR\n");
    ok = 0;
  }
  if (had_saved_row && saved_row == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_ROW_REFACTOR\n");
    ok = 0;
  }
  if (had_saved_checked && saved_checked == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_CHECKED_ROW_REFACTOR\n");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_KLS_FIRST_FACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_KLS_FIRST_FACTOR");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR=0");
    ok = 0;
  }

  if (ok && !require_ok(kls_create(&solver),
                        "create row-up lazy panel prefix")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(
                          solver, KLS_INDEX_INT32, 4, ap, ai, 0, &options),
                        "analyze row-up lazy panel prefix")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax),
                        "factor row-up lazy panel prefix")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve row-up lazy panel prefix")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats row-up lazy panel prefix")) {
    ok = 0;
  }
  if (ok && (stats.last_factor_path != KLS_FACTOR_PATH_KLS_FIRST ||
             stats.kls_first_last_row_uplooking_columns != 4 ||
             stats.kls_first_last_row_supernode_update != 1 ||
             stats.kls_first_last_row_supernode_update_rows < 4 ||
             stats.kls_first_last_row_supernode_panel_update != 0 ||
             stats.kls_first_row_supernode_panel_update_run_count != 0 ||
             stats.kls_first_last_row_supernode_panel_update_groups != 0 ||
             stats.kls_first_last_row_supernode_panel_update_rows != 0 ||
             stats.kls_first_row_panel_cache_append_count != 0 ||
             stats.kls_first_row_panel_cache_append_panels != 0 ||
             stats.kls_first_row_panel_cache_append_entries != 0)) {
    fprintf(stderr,
            "unexpected row-up lazy panel stats: path=%s, row_cols=%" PRId64
            ", row_supernode=%d/%" PRId64 ", panel=%d/%" PRId64
            "/%" PRId64 "/%" PRId64
            ", cache_append=%" PRId64 "/%" PRId64 "/%" PRId64 "\n",
            kls_factor_path_name(stats.last_factor_path),
            stats.kls_first_last_row_uplooking_columns,
            stats.kls_first_last_row_supernode_update,
            stats.kls_first_last_row_supernode_update_rows,
            stats.kls_first_last_row_supernode_panel_update,
            stats.kls_first_row_supernode_panel_update_run_count,
            stats.kls_first_last_row_supernode_panel_update_groups,
            stats.kls_first_last_row_supernode_panel_update_rows,
            stats.kls_first_row_panel_cache_append_count,
            stats.kls_first_row_panel_cache_append_panels,
            stats.kls_first_row_panel_cache_append_entries);
    ok = 0;
  }
  if (ok && (!close_enough(x[0], 1.0) || !close_enough(x[1], 2.0) ||
             !close_enough(x[2], 3.0) || !close_enough(x[3], 4.0))) {
    fprintf(stderr,
            "unexpected row-up lazy panel solution:"
            " %.17g %.17g %.17g %.17g\n",
            x[0], x[1], x[2], x[3]);
    ok = 0;
  }

  if (!restore_env_value("KLS_ENABLE_KLS_FIRST_FACTOR", had_saved_first,
                         saved_first)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_ROW_REFACTOR", had_saved_row,
                         saved_row)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_CHECKED_ROW_REFACTOR",
                         had_saved_checked, saved_checked)) {
    ok = 0;
  }
  kls_destroy(solver);
  free(saved_first);
  free(saved_row);
  free(saved_checked);
  return ok;
}

static int test_experimental_row_uplooking_dynamic_column_pivot(void) {
  const int32_t ap[] = {0, 2, 4, 5};
  const int32_t ai[] = {0, 1, 0, 1, 2};
  const double ax[] = {0.001, 2.0, 10.0, 1.0, 3.0};
  const double b[] = {20.001, 4.0, 9.0};
  double x[3] = {0.0, 0.0, 0.0};

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_NATURAL;
  options.use_btf = 0;
  options.scale = 0;
  options.static_pivoting = 0;
  options.pivot_tolerance = 0.1;
  options.threads = 2;

  const char *saved_first_value = getenv("KLS_ENABLE_KLS_FIRST_FACTOR");
  char *saved_first =
    saved_first_value != NULL ? strdup(saved_first_value) : NULL;
  const int had_saved_first = saved_first_value != NULL;
  const char *saved_row_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_row = saved_row_value != NULL ? strdup(saved_row_value) : NULL;
  const int had_saved_row = saved_row_value != NULL;
  const char *saved_checked_value =
    getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  char *saved_checked =
    saved_checked_value != NULL ? strdup(saved_checked_value) : NULL;
  const int had_saved_checked = saved_checked_value != NULL;

  int ok = 1;
  if (had_saved_first && saved_first == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_KLS_FIRST_FACTOR\n");
    ok = 0;
  }
  if (had_saved_row && saved_row == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_ROW_REFACTOR\n");
    ok = 0;
  }
  if (had_saved_checked && saved_checked == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_CHECKED_ROW_REFACTOR\n");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_KLS_FIRST_FACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_KLS_FIRST_FACTOR");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR=0");
    ok = 0;
  }

  if (ok && !require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 3, ap, ai, 0,
                                        &options),
                        "analyze row-up-looking dynamic pivot")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax),
                        "factor row-up-looking dynamic pivot")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve row-up-looking dynamic pivot")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats row-up-looking dynamic pivot")) {
    ok = 0;
  }
  if (ok && (stats.last_factor_path != KLS_FACTOR_PATH_KLS_FIRST ||
             stats.kls_first_last_row_uplooking_columns != 3 ||
             stats.kls_first_last_row_refactor_seeded_rows != 3 ||
             stats.kls_first_row_refactor_seeded_row_count < 3 ||
             stats.kls_first_last_dynamic_column_pivots != 1 ||
             stats.kls_first_dynamic_column_pivot_count < 1 ||
             stats.kls_first_last_row_pipeline != 1 ||
             stats.kls_first_last_row_pipeline_pivot_tail != 1 ||
             stats.kls_first_row_pipeline_pivot_tail_run_count < 1 ||
             stats.kls_first_last_row_pipeline_pivot_tail_rows <= 0 ||
             stats.kls_first_last_row_pipeline_pivot_restarts <= 0 ||
             stats.kls_first_last_row_pipeline_pivot_serial_rows != 0 ||
             stats.kls_first_last_row_pipeline_prefix_panel_rebuild != 0 ||
             stats.kls_first_row_pipeline_prefix_panel_rebuild_count != 0 ||
             stats.kls_first_last_row_pipeline_prefix_panel_rebuild_rows !=
               0 ||
             stats.kls_first_last_separator_queue_parallel_pipeline != 0 ||
             stats.kls_first_last_separator_queue_pipeline_prefix_panel_rebuild !=
               0 ||
             stats.selected_scale != 0 ||
             stats.selected_btf != 0)) {
    fprintf(stderr,
            "unexpected row-up-looking dynamic pivot stats: path=%s"
            ", row_cols=%" PRId64 ", row_seed=%" PRId64 "/%" PRId64
            ", dyn_pivots=%" PRId64 "/%" PRId64
            ", row_pipe=%d"
            ", row_pivot_tail=%d/%" PRId64 "/%" PRId64
            ", restarts=%" PRId64
            ", serial=%" PRId64
            ", prefix_rebuild=%d/%" PRId64 "/%" PRId64
            ", sep_pipe=%d/%d"
            ", scale=%d, btf=%d\n",
            kls_factor_path_name(stats.last_factor_path),
            stats.kls_first_last_row_uplooking_columns,
            stats.kls_first_last_row_refactor_seeded_rows,
            stats.kls_first_row_refactor_seeded_row_count,
            stats.kls_first_last_dynamic_column_pivots,
            stats.kls_first_dynamic_column_pivot_count,
            stats.kls_first_last_row_pipeline,
            stats.kls_first_last_row_pipeline_pivot_tail,
            stats.kls_first_row_pipeline_pivot_tail_run_count,
            stats.kls_first_last_row_pipeline_pivot_tail_rows,
            stats.kls_first_last_row_pipeline_pivot_restarts,
            stats.kls_first_last_row_pipeline_pivot_serial_rows,
            stats.kls_first_last_row_pipeline_prefix_panel_rebuild,
            stats.kls_first_row_pipeline_prefix_panel_rebuild_count,
            stats.kls_first_last_row_pipeline_prefix_panel_rebuild_rows,
            stats.kls_first_last_separator_queue_parallel_pipeline,
            stats.kls_first_last_separator_queue_pipeline_prefix_panel_rebuild,
            stats.selected_scale,
            stats.selected_btf);
    ok = 0;
  }
  if (ok && (!close_enough(x[0], 1.0) || !close_enough(x[1], 2.0) ||
             !close_enough(x[2], 3.0))) {
    fprintf(stderr,
            "unexpected row-up-looking dynamic pivot solution:"
            " %.17g %.17g %.17g\n",
            x[0], x[1], x[2]);
    ok = 0;
  }

  if (had_saved_first && saved_first != NULL) {
    if (setenv("KLS_ENABLE_KLS_FIRST_FACTOR", saved_first, 1) != 0) {
      perror("restore KLS_ENABLE_KLS_FIRST_FACTOR");
      ok = 0;
    }
  } else if (!had_saved_first) {
    if (unsetenv("KLS_ENABLE_KLS_FIRST_FACTOR") != 0) {
      perror("unsetenv KLS_ENABLE_KLS_FIRST_FACTOR");
      ok = 0;
    }
  }
  if (had_saved_row && saved_row != NULL) {
    if (setenv("KLS_ENABLE_ROW_REFACTOR", saved_row, 1) != 0) {
      perror("restore KLS_ENABLE_ROW_REFACTOR");
      ok = 0;
    }
  } else if (!had_saved_row) {
    if (unsetenv("KLS_ENABLE_ROW_REFACTOR") != 0) {
      perror("unsetenv KLS_ENABLE_ROW_REFACTOR");
      ok = 0;
    }
  }
  if (had_saved_checked && saved_checked != NULL) {
    if (setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", saved_checked, 1) != 0) {
      perror("restore KLS_ENABLE_CHECKED_ROW_REFACTOR");
      ok = 0;
    }
  } else if (!had_saved_checked) {
    if (unsetenv("KLS_ENABLE_CHECKED_ROW_REFACTOR") != 0) {
      perror("unsetenv KLS_ENABLE_CHECKED_ROW_REFACTOR");
      ok = 0;
    }
  }

  kls_destroy(solver);
  free(saved_first);
  free(saved_row);
  free(saved_checked);
  return ok;
}

static int test_experimental_row_uplooking_btf_blocks(void) {
  const int32_t ap[] = {0, 2, 4, 7, 10};
  const int32_t ai[] = {0, 1, 0, 1, 0, 2, 3, 1, 2, 3};
  const double ax[] = {
    5.0, 0.5,
    1.0, 4.0,
    0.25, 0.001, 2.0,
    0.5, 10.0, 1.0
  };
  const double b[] = {7.75, 10.5, 40.003, 10.0};
  double x[4] = {0.0, 0.0, 0.0, 0.0};

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_NATURAL;
  options.use_btf = 1;
  options.scale = 0;
  options.static_pivoting = 0;
  options.pivot_tolerance = 0.1;
  options.threads = 2;

  const char *saved_first_value = getenv("KLS_ENABLE_KLS_FIRST_FACTOR");
  char *saved_first =
    saved_first_value != NULL ? strdup(saved_first_value) : NULL;
  const int had_saved_first = saved_first_value != NULL;
  const char *saved_row_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_row = saved_row_value != NULL ? strdup(saved_row_value) : NULL;
  const int had_saved_row = saved_row_value != NULL;
  const char *saved_checked_value =
    getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  char *saved_checked =
    saved_checked_value != NULL ? strdup(saved_checked_value) : NULL;
  const int had_saved_checked = saved_checked_value != NULL;

  int ok = 1;
  if (had_saved_first && saved_first == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_KLS_FIRST_FACTOR\n");
    ok = 0;
  }
  if (had_saved_row && saved_row == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_ROW_REFACTOR\n");
    ok = 0;
  }
  if (had_saved_checked && saved_checked == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_CHECKED_ROW_REFACTOR\n");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_KLS_FIRST_FACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_KLS_FIRST_FACTOR");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR=0");
    ok = 0;
  }

  if (ok && !require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 4, ap, ai, 0,
                                        &options),
                        "analyze row-up-looking BTF blocks")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax),
                        "factor row-up-looking BTF blocks")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve row-up-looking BTF blocks")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats row-up-looking BTF blocks")) {
    ok = 0;
  }
  if (ok && (stats.last_factor_path != KLS_FACTOR_PATH_KLS_FIRST ||
             stats.nblocks != 2 ||
             stats.kls_first_last_row_uplooking_columns != 4 ||
             stats.kls_first_last_row_refactor_seeded_rows != 4 ||
             stats.kls_first_row_refactor_seeded_row_count < 4 ||
             stats.kls_first_last_dynamic_column_pivots != 1 ||
             stats.kls_first_last_parallel_btf_blocks != 2 ||
             stats.kls_first_parallel_btf_block_count < 2 ||
             stats.selected_scale != 0 ||
             stats.selected_btf != 1)) {
    fprintf(stderr,
            "unexpected row-up-looking BTF stats: path=%s"
            ", nblocks=%" PRId64 ", row_cols=%" PRId64
            ", row_seed=%" PRId64 "/%" PRId64
            ", dyn_pivots=%" PRId64
            ", parallel_btf=%" PRId64 "/%" PRId64
            ", scale=%d, btf=%d\n",
            kls_factor_path_name(stats.last_factor_path),
            stats.nblocks,
            stats.kls_first_last_row_uplooking_columns,
            stats.kls_first_last_row_refactor_seeded_rows,
            stats.kls_first_row_refactor_seeded_row_count,
            stats.kls_first_last_dynamic_column_pivots,
            stats.kls_first_last_parallel_btf_blocks,
            stats.kls_first_parallel_btf_block_count,
            stats.selected_scale,
            stats.selected_btf);
    ok = 0;
  }
  if (ok && (!close_enough(x[0], 1.0) || !close_enough(x[1], 2.0) ||
             !close_enough(x[2], 3.0) || !close_enough(x[3], 4.0))) {
    fprintf(stderr,
            "unexpected row-up-looking BTF solution:"
            " %.17g %.17g %.17g %.17g\n",
            x[0], x[1], x[2], x[3]);
    ok = 0;
  }

  if (had_saved_first && saved_first != NULL) {
    if (setenv("KLS_ENABLE_KLS_FIRST_FACTOR", saved_first, 1) != 0) {
      perror("restore KLS_ENABLE_KLS_FIRST_FACTOR");
      ok = 0;
    }
  } else if (!had_saved_first) {
    if (unsetenv("KLS_ENABLE_KLS_FIRST_FACTOR") != 0) {
      perror("unsetenv KLS_ENABLE_KLS_FIRST_FACTOR");
      ok = 0;
    }
  }
  if (had_saved_row && saved_row != NULL) {
    if (setenv("KLS_ENABLE_ROW_REFACTOR", saved_row, 1) != 0) {
      perror("restore KLS_ENABLE_ROW_REFACTOR");
      ok = 0;
    }
  } else if (!had_saved_row) {
    if (unsetenv("KLS_ENABLE_ROW_REFACTOR") != 0) {
      perror("unsetenv KLS_ENABLE_ROW_REFACTOR");
      ok = 0;
    }
  }
  if (had_saved_checked && saved_checked != NULL) {
    if (setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", saved_checked, 1) != 0) {
      perror("restore KLS_ENABLE_CHECKED_ROW_REFACTOR");
      ok = 0;
    }
  } else if (!had_saved_checked) {
    if (unsetenv("KLS_ENABLE_CHECKED_ROW_REFACTOR") != 0) {
      perror("unsetenv KLS_ENABLE_CHECKED_ROW_REFACTOR");
      ok = 0;
    }
  }

  kls_destroy(solver);
  free(saved_first);
  free(saved_row);
  free(saved_checked);
  return ok;
}

static int test_kls_first_parallel_pivoted_btf_fallback(void) {
  const int32_t ap[] = {0, 2, 4, 6, 8};
  const int32_t ai[] = {0, 1, 0, 1, 2, 3, 2, 3};
  const double ax[] = {1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 2.0};

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_NATURAL;
  options.use_btf = 1;
  options.scale = -1;
  options.static_pivoting = 0;
  options.halt_if_singular = 0;
  options.threads = 2;

  const char *saved_first_value = getenv("KLS_ENABLE_KLS_FIRST_FACTOR");
  char *saved_first =
    saved_first_value != NULL ? strdup(saved_first_value) : NULL;
  const int had_saved_first = saved_first_value != NULL;
  const char *saved_row_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_row = saved_row_value != NULL ? strdup(saved_row_value) : NULL;
  const int had_saved_row = saved_row_value != NULL;
  const char *saved_checked_value =
    getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  char *saved_checked =
    saved_checked_value != NULL ? strdup(saved_checked_value) : NULL;
  const int had_saved_checked = saved_checked_value != NULL;

  int ok = 1;
  if (had_saved_first && saved_first == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_KLS_FIRST_FACTOR\n");
    ok = 0;
  }
  if (had_saved_row && saved_row == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_ROW_REFACTOR\n");
    ok = 0;
  }
  if (had_saved_checked && saved_checked == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_CHECKED_ROW_REFACTOR\n");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_KLS_FIRST_FACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_KLS_FIRST_FACTOR");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR=0");
    ok = 0;
  }

  if (ok && !require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 4, ap, ai, 0,
                                        &options),
                        "analyze pivoted BTF fallback")) ok = 0;
  if (ok) {
    const int status = kls_factor(solver, ax);
    if (status != KLS_ERR_SINGULAR) {
      fprintf(stderr, "expected singular pivoted BTF fallback, got %d\n",
              status);
      ok = 0;
    }
  }

  kls_stats stats;
  memset(&stats, 0, sizeof(stats));
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats pivoted BTF fallback")) {
    ok = 0;
  }
  if (ok && (stats.last_factor_path != KLS_FACTOR_PATH_KLS_FIRST ||
             stats.nblocks != 2 ||
             stats.kls_first_last_row_uplooking_columns != 0 ||
             stats.kls_first_last_parallel_btf_blocks != 2 ||
             stats.kls_first_parallel_btf_block_count < 2 ||
             stats.selected_scale != -1 ||
             stats.selected_btf != 1)) {
    fprintf(stderr,
            "unexpected pivoted BTF fallback stats: path=%s"
            ", nblocks=%" PRId64 ", row_cols=%" PRId64
            ", parallel_btf=%" PRId64 "/%" PRId64
            ", scale=%d, btf=%d\n",
            kls_factor_path_name(stats.last_factor_path),
            stats.nblocks,
            stats.kls_first_last_row_uplooking_columns,
            stats.kls_first_last_parallel_btf_blocks,
            stats.kls_first_parallel_btf_block_count,
            stats.selected_scale,
            stats.selected_btf);
    ok = 0;
  }

  if (had_saved_first && saved_first != NULL) {
    if (setenv("KLS_ENABLE_KLS_FIRST_FACTOR", saved_first, 1) != 0) {
      perror("restore KLS_ENABLE_KLS_FIRST_FACTOR");
      ok = 0;
    }
  } else if (!had_saved_first) {
    if (unsetenv("KLS_ENABLE_KLS_FIRST_FACTOR") != 0) {
      perror("unsetenv KLS_ENABLE_KLS_FIRST_FACTOR");
      ok = 0;
    }
  }
  if (had_saved_row && saved_row != NULL) {
    if (setenv("KLS_ENABLE_ROW_REFACTOR", saved_row, 1) != 0) {
      perror("restore KLS_ENABLE_ROW_REFACTOR");
      ok = 0;
    }
  } else if (!had_saved_row) {
    if (unsetenv("KLS_ENABLE_ROW_REFACTOR") != 0) {
      perror("unsetenv KLS_ENABLE_ROW_REFACTOR");
      ok = 0;
    }
  }
  if (had_saved_checked && saved_checked != NULL) {
    if (setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", saved_checked, 1) != 0) {
      perror("restore KLS_ENABLE_CHECKED_ROW_REFACTOR");
      ok = 0;
    }
  } else if (!had_saved_checked) {
    if (unsetenv("KLS_ENABLE_CHECKED_ROW_REFACTOR") != 0) {
      perror("unsetenv KLS_ENABLE_CHECKED_ROW_REFACTOR");
      ok = 0;
    }
  }

  kls_destroy(solver);
  free(saved_first);
  free(saved_row);
  free(saved_checked);
  return ok;
}

static int test_structurally_full_rank_singular_rejection(void) {
  enum {
    n = 15066,
    lower_n = 15064,
    extra_columns = 1938,
    nnz = 62198
  };
  int32_t *ap = (int32_t *)malloc((size_t)(n + 1) * sizeof(*ap));
  int32_t *ai = (int32_t *)malloc((size_t)nnz * sizeof(*ai));
  double *ax = (double *)malloc((size_t)nnz * sizeof(*ax));

  const char *saved_first_value = getenv("KLS_ENABLE_KLS_FIRST_FACTOR");
  char *saved_first = saved_first_value != NULL
    ? strdup(saved_first_value) : NULL;
  const int had_first = saved_first_value != NULL;

  int ok = 1;
  if (ap == NULL || ai == NULL || ax == NULL) {
    fprintf(stderr, "failed to allocate singular fingerprint regression\n");
    ok = 0;
  }
  if ((had_first && saved_first == NULL)) {
    fprintf(stderr, "failed to save singular fingerprint environment\n");
    ok = 0;
  }

  int32_t position = 0;
  for (int32_t col = 0; ok && col < n; ++col) {
    ap[col] = position;
    int32_t rows[5];
    double values[5];
    int32_t count = 0;
    if (col < 2) {
      rows[count] = 0;
      values[count++] = 1.0;
      rows[count] = 1;
      values[count++] = 1.0;
    } else {
      const int32_t local = col - 2;
      rows[count] = col;
      values[count++] = 4.0;
      rows[count] = 2 + (local + 1) % lower_n;
      values[count++] = 0.01;
      rows[count] = 2 + (local + 37) % lower_n;
      values[count++] = -0.02;
      rows[count] = 2 + (local + 701) % lower_n;
      values[count++] = 0.015;
      if (local < extra_columns) {
        rows[count] = 2 + (local + 1901) % lower_n;
        values[count++] = 0.005;
      }
    }
    for (int32_t i = 1; i < count; ++i) {
      const int32_t row = rows[i];
      const double value = values[i];
      int32_t j = i;
      while (j > 0 && rows[j - 1] > row) {
        rows[j] = rows[j - 1];
        values[j] = values[j - 1];
        --j;
      }
      rows[j] = row;
      values[j] = value;
    }
    if (position > nnz - count) {
      fprintf(stderr, "singular fingerprint generator overflow\n");
      ok = 0;
      break;
    }
    for (int32_t i = 0; i < count; ++i) {
      ai[position] = rows[i];
      ax[position++] = values[i];
    }
  }
  if (ok) {
    ap[n] = position;
    if (position != nnz) {
      fprintf(stderr,
              "singular fingerprint generator produced %d entries, expected %d\n",
              position, nnz);
      ok = 0;
    }
  }

  if (ok && setenv("KLS_ENABLE_KLS_FIRST_FACTOR", "0", 1) != 0) {
    perror("configure singular fingerprint regression");
    ok = 0;
  }

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.threads = 8;
  if (ok && !require_ok(kls_create(&solver),
                        "create singular fingerprint regression")) {
    ok = 0;
  }
  if (ok && !require_ok(kls_analyze_csc(
                          solver, KLS_INDEX_INT32, n, ap, ai, 0, &options),
                        "analyze singular fingerprint regression")) {
    ok = 0;
  }
  if (ok) {
    const int status = kls_factor(solver, ax);
    if (status != KLS_ERR_SINGULAR) {
      fprintf(stderr,
              "benchmark-shaped singular input was implicitly completed: %d\n",
              status);
      ok = 0;
    }
  }

  kls_destroy(solver);
  if (!restore_env_value("KLS_ENABLE_KLS_FIRST_FACTOR", had_first,
                         saved_first != NULL ? saved_first : "")) {
    ok = 0;
  }
  free(ap);
  free(ai);
  free(ax);
  free(saved_first);
  return ok;
}

static int test_parallel_btf_row_supernode_update(void) {
  const int32_t ap[] = {0, 3, 6, 9, 10};
  const int32_t ai[] = {
    0, 1, 2,
    0, 1, 2,
    0, 1, 2,
    3
  };
  const double ax[] = {
    5.0, 0.25, 0.125,
    0.5, 6.0, 0.375,
    0.25, 0.5, 7.0,
    8.0
  };
  const double b[] = {6.75, 13.75, 21.875, 32.0};
  double x[4] = {0.0, 0.0, 0.0, 0.0};

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_NATURAL;
  options.use_btf = 1;
  options.scale = 0;
  options.static_pivoting = 0;
  options.pivot_tolerance = 0.001;
  options.threads = 2;

  const char *saved_first_value = getenv("KLS_ENABLE_KLS_FIRST_FACTOR");
  char *saved_first =
    saved_first_value != NULL ? strdup(saved_first_value) : NULL;
  const int had_saved_first = saved_first_value != NULL;
  const char *saved_row_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_row = saved_row_value != NULL ? strdup(saved_row_value) : NULL;
  const int had_saved_row = saved_row_value != NULL;
  const char *saved_checked_value =
    getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  char *saved_checked =
    saved_checked_value != NULL ? strdup(saved_checked_value) : NULL;
  const int had_saved_checked = saved_checked_value != NULL;

  int ok = 1;
  if (had_saved_first && saved_first == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_KLS_FIRST_FACTOR\n");
    ok = 0;
  }
  if (had_saved_row && saved_row == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_ROW_REFACTOR\n");
    ok = 0;
  }
  if (had_saved_checked && saved_checked == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_CHECKED_ROW_REFACTOR\n");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_KLS_FIRST_FACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_KLS_FIRST_FACTOR");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR=0");
    ok = 0;
  }

  if (ok && !require_ok(kls_create(&solver),
                        "create parallel BTF row supernode")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(
                          solver, KLS_INDEX_INT32, 4, ap, ai, 0, &options),
                        "analyze parallel BTF row supernode")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax),
                        "factor parallel BTF row supernode")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve parallel BTF row supernode")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats parallel BTF row supernode")) {
    ok = 0;
  }
  if (ok && (stats.last_factor_path != KLS_FACTOR_PATH_KLS_FIRST ||
             stats.nblocks != 2 ||
             stats.kls_first_last_parallel_btf_blocks != 2 ||
             stats.kls_first_last_row_supernode_update != 1 ||
             stats.kls_first_row_supernode_update_run_count < 1 ||
             stats.kls_first_last_row_supernode_update_groups < 1 ||
             stats.kls_first_last_row_supernode_update_rows < 2 ||
             stats.selected_btf != 1 || stats.selected_scale != 0)) {
    fprintf(stderr,
            "unexpected parallel BTF row-supernode stats: path=%s"
            ", nblocks=%" PRId64 ", parallel_btf=%" PRId64
            ", row_supernode=%d/%" PRId64 "/%" PRId64 "/%" PRId64
            ", btf=%d, scale=%d\n",
            kls_factor_path_name(stats.last_factor_path),
            stats.nblocks,
            stats.kls_first_last_parallel_btf_blocks,
            stats.kls_first_last_row_supernode_update,
            stats.kls_first_row_supernode_update_run_count,
            stats.kls_first_last_row_supernode_update_groups,
            stats.kls_first_last_row_supernode_update_rows,
            stats.selected_btf,
            stats.selected_scale);
    ok = 0;
  }
  if (ok && (!close_enough(x[0], 1.0) || !close_enough(x[1], 2.0) ||
             !close_enough(x[2], 3.0) || !close_enough(x[3], 4.0))) {
    fprintf(stderr,
            "unexpected parallel BTF row-supernode solution:"
            " %.17g %.17g %.17g %.17g\n",
            x[0], x[1], x[2], x[3]);
    ok = 0;
  }

  if (!restore_env_value("KLS_ENABLE_KLS_FIRST_FACTOR", had_saved_first,
                         saved_first)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_ROW_REFACTOR", had_saved_row,
                         saved_row)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_CHECKED_ROW_REFACTOR",
                         had_saved_checked, saved_checked)) {
    ok = 0;
  }
  kls_destroy(solver);
  free(saved_first);
  free(saved_row);
  free(saved_checked);
  return ok;
}

static int test_row_solve_from_numeric_after_klu_first(void) {
  const int32_t ap[] = {0, 2, 5, 7};
  const int32_t ai[] = {0, 1, 0, 1, 2, 1, 2};
  const double ax[] = {4.0, 1.0, 1.0, 3.0, 1.0, 1.0, 2.0};
  const double b[] = {6.0, 10.0, 8.0};
  double x[3] = {0.0, 0.0, 0.0};
  double xt[3] = {0.0, 0.0, 0.0};

  kls_solver *solver = NULL;
  kls_options options;
  kls_default_options(&options);
  options.ordering = KLS_ORDERING_NATURAL;
  options.use_btf = 0;
  options.scale = -1;
  options.static_pivoting = 0;

  const char *saved_row_solve_env_value =
    getenv("KLS_ENABLE_ROW_SOLVE_FROM_NUMERIC");
  char *saved_row_solve_env =
    saved_row_solve_env_value != NULL ? strdup(saved_row_solve_env_value)
                                      : NULL;
  const int had_saved_row_solve_env = saved_row_solve_env_value != NULL;
  const char *saved_first_env_value = getenv("KLS_ENABLE_KLS_FIRST_FACTOR");
  char *saved_first_env =
    saved_first_env_value != NULL ? strdup(saved_first_env_value) : NULL;
  const int had_saved_first_env = saved_first_env_value != NULL;
  const char *saved_row_env_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_row_env =
    saved_row_env_value != NULL ? strdup(saved_row_env_value) : NULL;
  const int had_saved_row_env = saved_row_env_value != NULL;
  const char *saved_checked_env_value =
    getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  char *saved_checked_env =
    saved_checked_env_value != NULL ? strdup(saved_checked_env_value) : NULL;
  const int had_saved_checked_env = saved_checked_env_value != NULL;

  int ok = 1;
  if (had_saved_row_solve_env && saved_row_solve_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_ROW_SOLVE_FROM_NUMERIC\n");
    ok = 0;
  }
  if (had_saved_first_env && saved_first_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_KLS_FIRST_FACTOR\n");
    ok = 0;
  }
  if (had_saved_row_env && saved_row_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_ROW_REFACTOR\n");
    ok = 0;
  }
  if (had_saved_checked_env && saved_checked_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_CHECKED_ROW_REFACTOR\n");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_ROW_SOLVE_FROM_NUMERIC", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_SOLVE_FROM_NUMERIC");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_KLS_FIRST_FACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_KLS_FIRST_FACTOR=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR=0");
    ok = 0;
  }

  if (!require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 3, ap, ai, 0,
                                        &options),
                        "analyze row solve seed")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax),
                        "factor row solve seed")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats row solve seed factor")) {
    ok = 0;
  }
  if (ok && (stats.last_factor_path != KLS_FACTOR_PATH_KLU_FIRST ||
             stats.row_refactor_auto_enabled != 0 ||
             stats.row_refactor_auto_values_ready != 0 ||
             stats.row_refactor_values_dirty != 0 ||
             stats.row_refactor_group_count != 0 ||
             stats.row_refactor_segment_count != 0 ||
             stats.row_solve_partition_ready != 0 ||
             stats.row_solve_partition_slices != 0 ||
             stats.row_solve_l_slice_max_entries != 0 ||
             stats.row_solve_u_slice_max_entries != 0 ||
             stats.row_solve_l_segmented_rows != 0 ||
             stats.row_solve_l_rect_entries != 0 ||
             stats.row_solve_l_tri_entries != 0 ||
             stats.row_solve_u_segmented_rows != 0 ||
             stats.row_solve_u_rect_entries != 0 ||
             stats.row_solve_u_tri_entries != 0)) {
    fprintf(stderr,
            "unexpected row-solve seed factor stats: path=%s"
            ", auto=%d, ready=%d, dirty=%d, groups=%" PRId64
            ", segments=%" PRId64 ", solve_partition=%d/%" PRId64
            ", max_slices=%" PRId64 "/%" PRId64
            ", l_seg=%" PRId64 "/%" PRId64 "/%" PRId64
            ", u_seg=%" PRId64 "/%" PRId64 "/%" PRId64 "\n",
            kls_factor_path_name(stats.last_factor_path),
            stats.row_refactor_auto_enabled,
            stats.row_refactor_auto_values_ready,
            stats.row_refactor_values_dirty,
            stats.row_refactor_group_count,
            stats.row_refactor_segment_count,
            stats.row_solve_partition_ready,
            stats.row_solve_partition_slices,
            stats.row_solve_l_slice_max_entries,
            stats.row_solve_u_slice_max_entries,
            stats.row_solve_l_segmented_rows,
            stats.row_solve_l_rect_entries,
            stats.row_solve_l_tri_entries,
            stats.row_solve_u_segmented_rows,
            stats.row_solve_u_rect_entries,
            stats.row_solve_u_tri_entries);
    ok = 0;
  }

  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve row solve seed")) ok = 0;
  if (ok && !require_ok(kls_solve_transpose(solver, 1, b, 0, xt, 0),
                        "transpose solve row solve seed")) ok = 0;
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats row solve seed solves")) {
    ok = 0;
  }
  if (ok && (stats.row_refactor_last_row_solve != 0 ||
             stats.row_refactor_row_solve_run_count != 0 ||
             stats.row_solve_parallel_run_count != 0 ||
             stats.row_solve_parallel_l_slice_runs != 0 ||
             stats.row_solve_parallel_u_slice_runs != 0 ||
             stats.row_solve_parallel_l_sparse_level_runs != 0 ||
             stats.row_solve_parallel_u_sparse_level_runs != 0 ||
             stats.row_solve_thread_count != 0 ||
             stats.row_solve_l_thread_max_rect_entries != 0 ||
             stats.row_solve_u_thread_max_rect_entries != 0 ||
             stats.row_refactor_values_dirty != 0)) {
    fprintf(stderr,
            "unexpected row-solve seed solve stats: row_solve=%d/%" PRId64
            ", parallel=%" PRId64 "/%" PRId64 "/%" PRId64
            "/%" PRId64 "/%" PRId64
            ", thread_balance=%" PRId64 "/%" PRId64 "/%" PRId64
            ", dirty=%d\n",
            stats.row_refactor_last_row_solve,
            stats.row_refactor_row_solve_run_count,
            stats.row_solve_parallel_run_count,
            stats.row_solve_parallel_l_slice_runs,
            stats.row_solve_parallel_u_slice_runs,
            stats.row_solve_parallel_l_sparse_level_runs,
            stats.row_solve_parallel_u_sparse_level_runs,
            stats.row_solve_thread_count,
            stats.row_solve_l_thread_max_rect_entries,
            stats.row_solve_u_thread_max_rect_entries,
            stats.row_refactor_values_dirty);
    ok = 0;
  }
  if (ok && (!close_enough(x[0], 1.0) || !close_enough(x[1], 2.0) ||
             !close_enough(x[2], 3.0) ||
             !close_enough(xt[0], 1.0) || !close_enough(xt[1], 2.0) ||
             !close_enough(xt[2], 3.0))) {
    fprintf(stderr,
            "unexpected row-solve seed solutions:"
            " x=(%.17g %.17g %.17g), xt=(%.17g %.17g %.17g)\n",
            x[0], x[1], x[2], xt[0], xt[1], xt[2]);
    ok = 0;
  }

  if (had_saved_row_solve_env && saved_row_solve_env != NULL) {
    if (setenv("KLS_ENABLE_ROW_SOLVE_FROM_NUMERIC",
               saved_row_solve_env, 1) != 0) {
      perror("restore KLS_ENABLE_ROW_SOLVE_FROM_NUMERIC");
      ok = 0;
    }
  } else if (!had_saved_row_solve_env) {
    if (unsetenv("KLS_ENABLE_ROW_SOLVE_FROM_NUMERIC") != 0) {
      perror("unsetenv KLS_ENABLE_ROW_SOLVE_FROM_NUMERIC");
      ok = 0;
    }
  }
  if (had_saved_first_env && saved_first_env != NULL) {
    if (setenv("KLS_ENABLE_KLS_FIRST_FACTOR", saved_first_env, 1) != 0) {
      perror("restore KLS_ENABLE_KLS_FIRST_FACTOR");
      ok = 0;
    }
  } else if (!had_saved_first_env) {
    if (unsetenv("KLS_ENABLE_KLS_FIRST_FACTOR") != 0) {
      perror("unsetenv KLS_ENABLE_KLS_FIRST_FACTOR");
      ok = 0;
    }
  }
  if (had_saved_row_env && saved_row_env != NULL) {
    if (setenv("KLS_ENABLE_ROW_REFACTOR", saved_row_env, 1) != 0) {
      perror("restore KLS_ENABLE_ROW_REFACTOR");
      ok = 0;
    }
  } else if (!had_saved_row_env) {
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

  kls_destroy(solver);
  free(saved_row_solve_env);
  free(saved_first_env);
  free(saved_row_env);
  free(saved_checked_env);
  return ok;
}

static int test_transpose_row_solve_from_numeric_parallel(int repeated) {
  const int32_t n = 4000;
  const int64_t nnz64 = ((int64_t)n * ((int64_t)n + 1)) / 2;
  int ok = 1;
  int32_t *ap = NULL;
  int32_t *ai = NULL;
  double *ax = NULL;
  double *b = NULL;
  double *x = NULL;
  double *expected = NULL;
  kls_solver *solver = NULL;

  const char *saved_row_solve_value =
    getenv("KLS_ENABLE_ROW_SOLVE_FROM_NUMERIC");
  char *saved_row_solve =
    saved_row_solve_value != NULL ? strdup(saved_row_solve_value) : NULL;
  const int had_row_solve = saved_row_solve_value != NULL;
  const char *saved_first_value = getenv("KLS_ENABLE_KLS_FIRST_FACTOR");
  char *saved_first =
    saved_first_value != NULL ? strdup(saved_first_value) : NULL;
  const int had_first = saved_first_value != NULL;
  const char *saved_row_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_row =
    saved_row_value != NULL ? strdup(saved_row_value) : NULL;
  const int had_row = saved_row_value != NULL;
  const char *saved_checked_value = getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  char *saved_checked =
    saved_checked_value != NULL ? strdup(saved_checked_value) : NULL;
  const int had_checked = saved_checked_value != NULL;

  if (nnz64 > INT32_MAX) {
    fprintf(stderr, "transpose row-solve fixture too large\n");
    ok = 0;
  }
  if ((had_row_solve && saved_row_solve == NULL) ||
      (had_first && saved_first == NULL) ||
      (had_row && saved_row == NULL) ||
      (had_checked && saved_checked == NULL)) {
    fprintf(stderr, "failed to save row-solve environment\n");
    ok = 0;
  }

  if (ok) {
    ap = (int32_t *)malloc(((size_t)n + 1u) * sizeof(*ap));
    ai = (int32_t *)malloc((size_t)nnz64 * sizeof(*ai));
    ax = (double *)malloc((size_t)nnz64 * sizeof(*ax));
    b = (double *)calloc((size_t)n, sizeof(*b));
    x = (double *)calloc((size_t)n, sizeof(*x));
    expected = (double *)malloc((size_t)n * sizeof(*expected));
    if (ap == NULL || ai == NULL || ax == NULL || b == NULL ||
        x == NULL || expected == NULL) {
      fprintf(stderr, "failed to allocate transpose row-solve fixture\n");
      ok = 0;
    }
  }

  if (ok) {
    for (int32_t i = 0; i < n; ++i) {
      expected[i] = 1.0 + 0.001 * (double)(i % 23);
    }
    int32_t pos = 0;
    for (int32_t col = 0; col < n; ++col) {
      ap[col] = pos;
      for (int32_t row = 0; row <= col; ++row) {
        const double value =
          row == col ? 2.0
                     : 1.0e-4 * (double)(1 + ((row + col) % 7));
        ai[pos] = row;
        ax[pos] = value;
        b[col] += value * expected[row];
        pos++;
      }
    }
    ap[n] = pos;
  }

  if (ok && setenv("KLS_ENABLE_ROW_SOLVE_FROM_NUMERIC", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_SOLVE_FROM_NUMERIC");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_KLS_FIRST_FACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_KLS_FIRST_FACTOR");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR");
    ok = 0;
  }

  kls_options options;
  kls_default_options(&options);
  options.threads = 4;
  options.ordering = KLS_ORDERING_NATURAL;
  options.use_btf = 0;
  options.scale = 0;
  options.static_pivoting = 0;

  /* Exercise both ordinary scheduling and repeated-workload worker pinning.
     require_ok checks caller-affinity restoration after each public call. */
  options.expected_refactorizations = repeated ? 100 : 0;
  options.expected_solves = repeated ? 100 : 0;

  if (ok && !require_ok(kls_create(&solver),
                        "create transpose parallel row solve")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, n,
                                        ap, ai, 0, &options),
                        "analyze transpose parallel row solve")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax),
                        "factor transpose parallel row solve")) ok = 0;
  /* engine/solve preps are deferred to the first refactorization; the
     seeded partition is part of that contract */
  if (ok && !require_ok(kls_refactor(solver, ax),
                        "refactor transpose parallel row solve")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "factor stats transpose parallel row solve")) ok = 0;
  if (ok && stats.row_solve_partition_ready != 1) {
    fprintf(stderr,
            "transpose row-solve factor did not seed partition: ready=%d\n",
            stats.row_solve_partition_ready);
    ok = 0;
  }

  if (ok && !require_ok(kls_solve_transpose(solver, 1, b, 0, x, 0),
                        "transpose parallel row solve")) ok = 0;
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "transpose parallel row solve stats")) ok = 0;
  if (ok && (stats.row_refactor_last_row_solve != 1 ||
             stats.row_refactor_row_solve_run_count != 1 ||
             stats.row_solve_parallel_run_count != 1 ||
             (stats.row_solve_parallel_l_slice_runs <= 0 &&
              stats.row_solve_parallel_l_sparse_level_runs <= 0 &&
              stats.row_solve_parallel_u_slice_runs <= 0 &&
              stats.row_solve_parallel_u_sparse_level_runs <= 0))) {
    fprintf(stderr,
            "unexpected transpose parallel row-solve stats:"
            " row_solve=%d/%" PRId64 ", parallel=%" PRId64
            ", l_runs=%" PRId64 "/%" PRId64
            ", u_runs=%" PRId64 "/%" PRId64 "\n",
            stats.row_refactor_last_row_solve,
            stats.row_refactor_row_solve_run_count,
            stats.row_solve_parallel_run_count,
            stats.row_solve_parallel_l_slice_runs,
            stats.row_solve_parallel_l_sparse_level_runs,
            stats.row_solve_parallel_u_slice_runs,
            stats.row_solve_parallel_u_sparse_level_runs);
    ok = 0;
  }
  for (int32_t i = 0; ok && i < n; ++i) {
    if (fabs(x[i] - expected[i]) > 1.0e-8) {
      fprintf(stderr,
              "unexpected transpose parallel row-solve solution at %d:"
              " %.17g != %.17g\n",
              (int)i, x[i], expected[i]);
      ok = 0;
    }
  }

  if (!restore_env_value("KLS_ENABLE_ROW_SOLVE_FROM_NUMERIC",
                         had_row_solve, saved_row_solve)) ok = 0;
  if (!restore_env_value("KLS_ENABLE_KLS_FIRST_FACTOR",
                         had_first, saved_first)) ok = 0;
  if (!restore_env_value("KLS_ENABLE_ROW_REFACTOR", had_row, saved_row)) {
    ok = 0;
  }
  if (!restore_env_value("KLS_ENABLE_CHECKED_ROW_REFACTOR",
                         had_checked, saved_checked)) ok = 0;
  kls_destroy(solver);
  free(ap);
  free(ai);
  free(ax);
  free(b);
  free(x);
  free(expected);
  free(saved_row_solve);
  free(saved_first);
  free(saved_row);
  free(saved_checked);
  return ok;
}

static int test_kls_first_reseeds_after_pivot_repair(void) {
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
  options.scale = 0;
  options.pivot_tolerance = 0.001;
  options.static_pivoting = 0;

  const char *saved_env_value = getenv("KLS_ENABLE_KLS_FIRST_FACTOR");
  char *saved_env = saved_env_value != NULL ? strdup(saved_env_value) : NULL;
  const int had_saved_env = saved_env_value != NULL;
  const char *saved_row_env_value = getenv("KLS_ENABLE_ROW_REFACTOR");
  char *saved_row_env =
    saved_row_env_value != NULL ? strdup(saved_row_env_value) : NULL;
  const int had_saved_row_env = saved_row_env_value != NULL;
  const char *saved_checked_env_value =
    getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  char *saved_checked_env =
    saved_checked_env_value != NULL ? strdup(saved_checked_env_value) : NULL;
  const int had_saved_checked_env = saved_checked_env_value != NULL;

  int ok = 1;
  if (had_saved_env && saved_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_KLS_FIRST_FACTOR\n");
    ok = 0;
  }
  if (had_saved_row_env && saved_row_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_ROW_REFACTOR\n");
    ok = 0;
  }
  if (had_saved_checked_env && saved_checked_env == NULL) {
    fprintf(stderr, "failed to save KLS_ENABLE_CHECKED_ROW_REFACTOR\n");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_KLS_FIRST_FACTOR", "1", 1) != 0) {
    perror("setenv KLS_ENABLE_KLS_FIRST_FACTOR");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_ROW_REFACTOR=0");
    ok = 0;
  }
  if (ok && setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "0", 1) != 0) {
    perror("setenv KLS_ENABLE_CHECKED_ROW_REFACTOR=0");
    ok = 0;
  }

  if (!require_ok(kls_create(&solver), "create")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 3, ap, ai, 0,
                                        &options),
                        "analyze KLS-first repair reseed")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax0),
                        "factor KLS-first repair reseed base")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax1),
                        "factor KLS-first repair reseed pivot")) ok = 0;
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve KLS-first repair reseed")) ok = 0;

  kls_stats stats;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats KLS-first repair reseed")) {
    ok = 0;
  }
  if (ok && (stats.last_factor_path != KLS_FACTOR_PATH_KLS_FAST_REFACTOR ||
             stats.fast_block_restarts != 1 ||
             stats.fast_kls_block_restarts != 1 ||
             stats.fast_tail_restarts != 1 ||
             stats.row_refactor_values_dirty != 0 ||
             stats.row_refactor_last_row_solve != 1 ||
             stats.row_refactor_row_solve_run_count != 1)) {
    fprintf(stderr,
            "unexpected KLS-first repair reseed stats: path=%s"
            ", restarts=%d, kls_restarts=%d, tail_restarts=%d"
            ", dirty=%d, row_solve=%d/%" PRId64 "\n",
            kls_factor_path_name(stats.last_factor_path),
            stats.fast_block_restarts,
            stats.fast_kls_block_restarts,
            stats.fast_tail_restarts,
            stats.row_refactor_values_dirty,
            stats.row_refactor_last_row_solve,
            stats.row_refactor_row_solve_run_count);
    ok = 0;
  }
  if (ok && !require_pivoting_tail_plan(&stats,
                                        "KLS-first repair reseed")) {
    ok = 0;
  }
  if (ok && (!close_enough(x[0], 1.0) || !close_enough(x[1], 2.0) ||
             !close_enough(x[2], 3.0))) {
    fprintf(stderr,
            "unexpected KLS-first repair reseed solution:"
            " %.17g %.17g %.17g\n",
            x[0], x[1], x[2]);
    ok = 0;
  }

  if (had_saved_env && saved_env != NULL) {
    if (setenv("KLS_ENABLE_KLS_FIRST_FACTOR", saved_env, 1) != 0) {
      perror("restore KLS_ENABLE_KLS_FIRST_FACTOR");
      ok = 0;
    }
  } else if (!had_saved_env) {
    if (unsetenv("KLS_ENABLE_KLS_FIRST_FACTOR") != 0) {
      perror("unsetenv KLS_ENABLE_KLS_FIRST_FACTOR");
      ok = 0;
    }
  }
  if (had_saved_row_env && saved_row_env != NULL) {
    if (setenv("KLS_ENABLE_ROW_REFACTOR", saved_row_env, 1) != 0) {
      perror("restore KLS_ENABLE_ROW_REFACTOR");
      ok = 0;
    }
  } else if (!had_saved_row_env) {
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

  kls_destroy(solver);
  free(saved_env);
  free(saved_row_env);
  free(saved_checked_env);
  return ok;
}

/* supernodal panel factor (phase 1b): random R x W panels, verify
   P.A = L.U to 1e-11 with the row-map permutation */
extern long KLS_SN_PANEL_FACTOR(double *A, long R, long W, long *rowids,
                                const long *diagrows, double tol);

static int run_sn_panel_factor_test(void) {
  unsigned int seed = 12345u;
  int trial;
  for (trial = 0; trial < 20; ++trial) {
    long R = 4 + (long)(rand_r(&seed) % 60);
    long W = 1 + (long)(rand_r(&seed) % (R < 32 ? R : 32));
    long i, j, k;
    double *A = (double *)malloc((size_t)(R * W) * sizeof(double));
    double *orig = (double *)malloc((size_t)(R * W) * sizeof(double));
    long *rowids = (long *)malloc((size_t)R * sizeof(long));
    long *origids = (long *)malloc((size_t)R * sizeof(long));
    if (A == NULL || orig == NULL || rowids == NULL || origids == NULL) {
      free(A); free(orig); free(rowids); free(origids);
      return 0;
    }
    for (i = 0; i < R; ++i) {
      rowids[i] = origids[i] = i;
      for (j = 0; j < W; ++j) {
        A[i * W + j] =
          ((double)(rand_r(&seed) % 2000) - 1000.0) / 97.0;
      }
    }
    memcpy(orig, A, (size_t)(R * W) * sizeof(double));
    if (KLS_SN_PANEL_FACTOR(A, R, W, rowids, NULL, 0.001) != W) {
      /* exact zero column is possible but vanishingly unlikely here */
      free(A); free(orig); free(rowids); free(origids);
      continue;
    }
    /* verify: for each factored position, orig[rowids[i]][j] ==
       sum_k L[i][k] * U[k][j] (L unit lower incl. multipliers below,
       U upper W x W) */
    for (i = 0; i < R; ++i) {
      for (j = 0; j < W; ++j) {
        double acc = 0.0;
        for (k = 0; k < W && k <= j; ++k) {
          double l;
          if (i < W && k > i) break;
          l = (i == k) ? 1.0 : (k < i || i >= W ? A[i * W + k] : 0.0);
          if (k > i && i < W) l = 0.0;
          acc += l * A[k * W + j];
        }
        {
          double ref = orig[rowids[i] * W + j];
          double d = acc - ref;
          if (d < 0) d = -d;
          double scale = ref < 0 ? -ref : ref;
          if (scale < 1.0) scale = 1.0;
          if (d / scale > 1e-11) {
            fprintf(stderr,
                    "sn panel mismatch trial %d (%ld x %ld) i=%ld j=%ld "
                    "acc=%g ref=%g\n", trial, R, W, i, j, acc, ref);
            free(A); free(orig); free(rowids); free(origids);
            return 0;
          }
        }
      }
    }
    free(A); free(orig); free(rowids); free(origids);
  }
  return 1;
}

/* A solve immediately after factorization builds the compact i32 row-index
   streams.  Large factors defer their engine preparations until the first
   changed refactorization, and that preparation sorts the packed LU columns.
   The sort must invalidate and rebuild the streams before the next solve. */
static int test_deferred_sort_rebuilds_compact_solve(void) {
  enum { block_size = 64, block_count = 256, n = block_size * block_count };
  const size_t nnz = (size_t)n * block_size;
  const char *env_names[] = {
    "KLS_ENABLE_KLS_FIRST_FACTOR",
    "KLS_ENABLE_ROW_REFACTOR",
    "KLS_ENABLE_CHECKED_ROW_REFACTOR",
    "KLS_ENABLE_ROW_SOLVE_FROM_NUMERIC",
    "KLS_SYNC_FACTOR_PREPS",
    "KLS_DISABLE_I32_SOLVE",
    "KLS_DISABLE_SNB_REFACTOR",
    "KLS_SNB_WMAX"
  };
  enum { env_count = (int)(sizeof(env_names) / sizeof(env_names[0])) };
  char *saved[env_count];
  int had[env_count];
  int32_t *ap = NULL;
  int32_t *ai = NULL;
  double *ax = NULL;
  double *b = NULL;
  double *x = NULL;
  double *expected = NULL;
  kls_solver *solver = NULL;
  int ok = 1;

  memset(saved, 0, sizeof(saved));
  memset(had, 0, sizeof(had));
  for (int e = 0; e < env_count; ++e) {
    const char *value = getenv(env_names[e]);
    had[e] = value != NULL;
    saved[e] = value != NULL ? strdup(value) : NULL;
    if (had[e] && saved[e] == NULL) {
      ok = 0;
    }
  }
  if (ok &&
      (setenv("KLS_ENABLE_KLS_FIRST_FACTOR", "0", 1) != 0 ||
       setenv("KLS_ENABLE_ROW_REFACTOR", "0", 1) != 0 ||
       setenv("KLS_ENABLE_CHECKED_ROW_REFACTOR", "0", 1) != 0 ||
       setenv("KLS_ENABLE_ROW_SOLVE_FROM_NUMERIC", "0", 1) != 0 ||
       unsetenv("KLS_SYNC_FACTOR_PREPS") != 0 ||
       unsetenv("KLS_DISABLE_I32_SOLVE") != 0 ||
       unsetenv("KLS_DISABLE_SNB_REFACTOR") != 0 ||
       setenv("KLS_SNB_WMAX", "4", 1) != 0)) {
    perror("configure deferred-sort compact-solve test");
    ok = 0;
  }

  if (ok) {
    ap = (int32_t *)malloc(((size_t)n + 1u) * sizeof(*ap));
    ai = (int32_t *)malloc(nnz * sizeof(*ai));
    ax = (double *)malloc(nnz * sizeof(*ax));
    b = (double *)calloc((size_t)n, sizeof(*b));
    x = (double *)calloc((size_t)n, sizeof(*x));
    expected = (double *)malloc((size_t)n * sizeof(*expected));
    if (ap == NULL || ai == NULL || ax == NULL || b == NULL || x == NULL ||
        expected == NULL) {
      fprintf(stderr, "failed to allocate deferred-sort fixture\n");
      ok = 0;
    }
  }

  if (ok) {
    size_t pos = 0u;
    for (int32_t col = 0; col < n; ++col) {
      const int32_t block = col / block_size;
      const int32_t local_col = col % block_size;
      const int32_t first = block * block_size;
      ap[col] = (int32_t)pos;
      expected[col] = 1.0 + 0.001 * (double)(col % 29);
      /* Four-column SNB panels split each dense block across 16 producers,
         exercising first/interior/last producer and subcolumn lookups in
         both packed-U writeback and input scatter construction.
         Deliberately permute input rows.  KLU accepts unsorted CSC input and
         its packed factors therefore give the deferred sort real work. */
      for (int32_t q = 0; q < block_size; ++q) {
        const int32_t local_row = (17 * q + 7 * local_col) % block_size;
        const int32_t row = first + local_row;
        ai[pos] = row;
        ax[pos] = local_row == local_col
                    ? 8.0
                    : 0.002 * (double)(1 + ((13 * local_row +
                                              11 * local_col) % 11));
        ++pos;
      }
    }
    ap[n] = (int32_t)pos;
    if (pos != nnz) {
      fprintf(stderr, "deferred-sort fixture nnz mismatch\n");
      ok = 0;
    }
  }

  kls_options options;
  kls_default_options(&options);
  options.backend = KLS_BACKEND_AUTO;
  options.threads = 1;
  options.ordering = KLS_ORDERING_NATURAL;
  options.orientation = KLS_ORIENTATION_NORMAL;
  options.scale = -1;
  options.use_btf = 0;
  options.static_pivoting = 0;

  if (ok && !require_ok(kls_create(&solver),
                        "create deferred-sort compact solve")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, n, ap, ai,
                                        0, &options),
                        "analyze deferred-sort compact solve")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax),
                        "factor deferred-sort compact solve")) ok = 0;

  if (ok) {
    for (int32_t col = 0; col < n; ++col) {
      for (int32_t p2 = ap[col]; p2 < ap[col + 1]; ++p2) {
        b[ai[p2]] += ax[p2] * expected[col];
      }
    }
    /* This is the cache-building solve that precedes deferred preparation. */
    if (!require_ok(kls_solve(solver, 1, b, 0, x, 0),
                    "pre-sort compact solve")) {
      ok = 0;
    }
  }
  for (int32_t i = 0; ok && i < n; ++i) {
    if (fabs(x[i] - expected[i]) > 1.0e-9) {
      fprintf(stderr, "pre-sort compact solve mismatch at %d: %.17g\n",
              (int)i, x[i]);
      ok = 0;
    }
  }

  if (ok) {
    memset(b, 0, (size_t)n * sizeof(*b));
    memset(x, 0, (size_t)n * sizeof(*x));
    for (size_t p2 = 0; p2 < nnz; ++p2) {
      const int delta = (int)(p2 % 17u) - 8;
      ax[p2] *= 1.0 + 1.0e-4 * (double)delta;
    }
    for (int32_t col = 0; col < n; ++col) {
      for (int32_t p2 = ap[col]; p2 < ap[col + 1]; ++p2) {
        b[ai[p2]] += ax[p2] * expected[col];
      }
    }
    if (!require_ok(kls_refactor(solver, ax),
                    "deferred sorted refactor")) {
      ok = 0;
    }
    if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                          "post-sort compact solve")) {
      ok = 0;
    }
  }
  for (int32_t i = 0; ok && i < n; ++i) {
    if (fabs(x[i] - expected[i]) > 1.0e-8) {
      fprintf(stderr, "post-sort compact solve mismatch at %d: %.17g\n",
              (int)i, x[i]);
      ok = 0;
    }
  }

  if (ok) {
    memset(b, 0, (size_t)n * sizeof(*b));
    for (int32_t col = 0; col < n; ++col) {
      for (int32_t p2 = ap[col]; p2 < ap[col + 1]; ++p2) {
        b[col] += ax[p2] * expected[ai[p2]];
      }
    }
    if (!require_ok(kls_solve_transpose(solver, 1, b, 0, x, 0),
                    "post-sort transpose solve")) ok = 0;
    for (int32_t i = 0; ok && i < n; ++i) {
      if (fabs(x[i] - expected[i]) > 1.0e-8) {
        fprintf(stderr, "post-sort transpose mismatch at %d\n", (int)i);
        ok = 0;
      }
    }
  }

  kls_destroy(solver);
  free(ap);
  free(ai);
  free(ax);
  free(b);
  free(x);
  free(expected);
  for (int e = 0; e < env_count; ++e) {
    if (!restore_env_value(env_names[e], had[e],
                           saved[e] != NULL ? saved[e] : "")) {
      ok = 0;
    }
    free(saved[e]);
  }
  return ok;
}

static int test_exact_unchanged_refactor_reuse(void) {
  const int32_t n = 3;
  const int32_t ap[4] = {0, 2, 5, 7};
  const int32_t ai[7] = {0, 1, 0, 1, 2, 1, 2};
  double ax[7] = {4.0, 2.0, 1.0, 5.0, 3.0, 1.0, 6.0};
  const double expected[3] = {1.0, 2.0, 3.0};
  double b[3] = {0.0, 0.0, 0.0};
  double x[3] = {0.0, 0.0, 0.0};
  const char *saved_value = getenv("KLS_DISABLE_UNCHANGED_REFACTOR");
  char *saved = saved_value != NULL ? strdup(saved_value) : NULL;
  const int had_saved = saved_value != NULL;
  kls_solver *solver = NULL;
  kls_options options;
  kls_stats stats;
  int ok = 1;

  if (had_saved && saved == NULL) {
    return 0;
  }
  if (unsetenv("KLS_DISABLE_UNCHANGED_REFACTOR") != 0) {
    perror("unsetenv KLS_DISABLE_UNCHANGED_REFACTOR");
    free(saved);
    return 0;
  }
  kls_default_options(&options);
  options.threads = 1;
  options.ordering = KLS_ORDERING_NATURAL;
  options.orientation = KLS_ORIENTATION_NORMAL;
  options.scale = -1;
  options.use_btf = 0;
  options.static_pivoting = 0;

  if (!require_ok(kls_create(&solver), "create unchanged-refactor")) ok = 0;
  if (ok && !require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, n, ap, ai,
                                        0, &options),
                        "analyze unchanged-refactor")) ok = 0;
  if (ok && !require_ok(kls_factor(solver, ax),
                        "factor unchanged-refactor")) ok = 0;
  if (ok && !require_ok(kls_refactor(solver, ax),
                        "reuse unchanged factor")) ok = 0;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats unchanged factor")) ok = 0;
  if (ok && stats.last_refactor_path != KLS_REFACTOR_PATH_UNCHANGED) {
    fprintf(stderr, "unchanged refactor did not reuse factor: %s\n",
            kls_refactor_path_name(stats.last_refactor_path));
    ok = 0;
  }

  /* Change the final entry in place so the exact comparison must inspect the
     whole array, then verify that the changed matrix is actually factored. */
  ax[6] = 7.0;
  for (int col = 0; col < n; ++col) {
    for (int32_t p = ap[col]; p < ap[col + 1]; ++p) {
      b[ai[p]] += ax[p] * expected[col];
    }
  }
  if (ok && !require_ok(kls_refactor(solver, ax),
                        "changed after unchanged refactor")) ok = 0;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats changed after unchanged")) ok = 0;
  if (ok && stats.last_refactor_path == KLS_REFACTOR_PATH_UNCHANGED) {
    fprintf(stderr, "changed values incorrectly reused factor\n");
    ok = 0;
  }
  if (ok && !require_ok(kls_solve(solver, 1, b, 0, x, 0),
                        "solve changed after unchanged")) ok = 0;
  for (int i = 0; ok && i < n; ++i) {
    if (!close_enough(x[i], expected[i])) {
      fprintf(stderr,
              "changed-after-unchanged solution mismatch at %d: %.17g\n",
              i, x[i]);
      ok = 0;
    }
  }

  /* Once repetition was observed, a changed step refreshes the snapshot and
     the following identical step is eligible again. */
  if (ok && !require_ok(kls_refactor(solver, ax),
                        "reuse changed snapshot")) ok = 0;
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "stats reused changed snapshot")) ok = 0;
  if (ok && stats.last_refactor_path != KLS_REFACTOR_PATH_UNCHANGED) {
    fprintf(stderr, "refreshed unchanged snapshot was not reused: %s\n",
            kls_refactor_path_name(stats.last_refactor_path));
    ok = 0;
  }

  kls_destroy(solver);
  if (!restore_env_value("KLS_DISABLE_UNCHANGED_REFACTOR", had_saved,
                         saved != NULL ? saved : "")) {
    ok = 0;
  }
  free(saved);
  return ok;
}

/* Checker-mode coverage for the experimental partial-BTF path.  The fast
   solver and a second solver forced through the ordinary refactor see each
   generated numeric and RHS; their complete solutions must agree. */

static int removed_shape_stats_are_zero(const kls_stats *stats) {
  return stats != NULL &&
    stats->moderate_fragmented_policy_eligible == 0 &&
    stats->compact_amf_two_block_policy_eligible == 0 &&
    stats->dense_reciprocal_hub_policy_eligible == 0 &&
    stats->symmetric_scalar_fringe_policy_eligible == 0 &&
    stats->pivoted_high_work_single_block_policy_eligible == 0 &&
    stats->low_work_many_fringe_btf_pts_policy_eligible == 0 &&
    stats->high_work_tiny_scalar_fringe_policy_eligible == 0 &&
    stats->low_work_hubbed_scalar_fringe_pts_policy_eligible == 0 &&
    stats->nearly_missing_diagonal_early_match_selected == 0 &&
    stats->low_work_tiny_block_btf_policy_eligible == 0 &&
    stats->low_work_tiny_block_btf_symbolic_eligible == 0 &&
    stats->scaled_fragmented_compact_row_policy_eligible == 0 &&
    stats->compact_missing_diagonal_match_candidate == 0 &&
    stats->compact_missing_diagonal_match_selected == 0 &&
    stats->compact_missing_diagonal_factor_eligible == 0 &&
    stats->symmetric_partial_diagonal_match_candidate == 0 &&
    stats->symmetric_partial_diagonal_match_selected == 0 &&
    stats->symmetric_partial_diagonal_factor_eligible == 0 &&
    stats->symmetric_partial_diagonal_low_work_eligible == 0 &&
    stats->sparse_symmetric_fragmented_metis_symbolic_eligible == 0 &&
    stats->sparse_symmetric_fragmented_metis_policy_eligible == 0 &&
    stats->sparse_spiked_predicted_candidate == 0 &&
    stats->sparse_spiked_predicted_factor_eligible == 0 &&
    stats->sparse_spiked_predicted_clustered_eligible == 0 &&
    stats->dense_fragmented_scaled_row_factor_eligible == 0 &&
    stats->sparse_full_diagonal_metis_row_candidate == 0 &&
    stats->sparse_full_diagonal_metis_row_symbolic_eligible == 0 &&
    stats->sparse_full_diagonal_metis_row_factor_eligible == 0 &&
    stats->giant_symmetric_scalar_fringe_metis_row_candidate == 0 &&
    stats->giant_symmetric_scalar_fringe_metis_row_symbolic_eligible == 0 &&
    stats->giant_symmetric_scalar_fringe_metis_row_factor_eligible == 0 &&
    stats->hybrid_huge_single_egraph_factor_eligible == 0 &&
    stats->bounded_degree_retained_preconditioner_candidate == 0 &&
    stats->bounded_degree_retained_preconditioner_symbolic_eligible == 0 &&
    stats->bounded_degree_retained_preconditioner_factor_eligible == 0 &&
    stats->bounded_degree_retained_preconditioner_reuse_count == 0 &&
    stats->asymmetric_bounded_degree_direct_metis_candidate == 0 &&
    stats->asymmetric_bounded_degree_direct_metis_tuning_class == 0 &&
    stats->asymmetric_bounded_degree_direct_metis_symbolic_eligible == 0 &&
    stats->asymmetric_bounded_degree_direct_metis_factor_eligible == 0 &&
    stats->near_symmetric_mega_hub_amd_candidate == 0 &&
    stats->near_symmetric_mega_hub_amd_symbolic_eligible == 0 &&
    stats->near_symmetric_mega_hub_amd_factor_eligible == 0 &&
    stats->giant_dominant_hub_metis_dense_tail_candidate == 0 &&
    stats->giant_dominant_hub_metis_dense_tail_symbolic_eligible == 0 &&
    stats->giant_dominant_hub_metis_dense_tail_factor_eligible == 0 &&
    stats->high_work_tiny_fringe_btf_pts_factor_eligible == 0 &&
    stats->medium_spike_minfill_candidate == 0 &&
    stats->medium_spike_minfill_symbolic_eligible == 0 &&
    stats->sparse_broad_column_amf_no_btf_candidate == 0 &&
    stats->sparse_broad_column_amf_no_btf_symbolic_eligible == 0 &&
    stats->bounded_degree_amf_no_btf_candidate == 0 &&
    stats->bounded_degree_amf_no_btf_symbolic_eligible == 0;
}

static int test_scaled_update_refactor(void) {
  const int32_t ap[] = {0, 3, 6, 9};
  const int32_t ai[] = {0, 1, 2, 0, 1, 2, 0, 1, 2};
  const double base[] = {4, 2, 1, 1, 5, 2, 2, 1, 6};
  const double expected[] = {1, -2, 0.5, -0.5, 1, 3};
  for (int mode = 0; mode < 4; ++mode) {
    kls_solver *solver = NULL;
    kls_options options;
    kls_default_options(&options);
    options.backend = mode < 2 ? KLS_BACKEND_SERIAL : KLS_BACKEND_KLS;
    options.threads = mode < 2 ? 1 : 4;
    options.orientation = mode % 2 ? KLS_ORIENTATION_TRANSPOSE
                                   : KLS_ORIENTATION_NORMAL;
    int ok = require_ok(kls_create(&solver), "create scaled-update solver");
    if (ok) ok = require_ok(kls_analyze_csc(solver, KLS_INDEX_INT32, 3,
                                           ap, ai, 0, &options),
                            "analyze scaled-update solver");
    if (ok) ok = require_ok(kls_factor(solver, base), "factor scaled-update base");
    for (int update = 0; ok && update < 4; ++update) {
      double ax[9];
      for (int col = 0; col < 3; ++col) {
        for (int p = ap[col]; p < ap[col + 1]; ++p) {
          ax[p] = base[p] * (1.0 + 0.25 * (update + ai[p])) *
                  (1.0 + 0.125 * (update + col));
        }
      }
      /* Alternate separable scaling and ordinary entrywise changes. */
      if (update % 2) ax[1] += 0.125;
      ok = require_ok(kls_refactor(solver, ax), "refactor scaled update");
      for (int transpose = 0; ok && transpose < 2; ++transpose) {
        double rhs[6] = {0};
        for (int j = 0; j < 2; ++j) {
          for (int col = 0; col < 3; ++col) {
            for (int p = ap[col]; p < ap[col + 1]; ++p) {
              const int dst = transpose ? col : ai[p];
              const int src = transpose ? ai[p] : col;
              rhs[3 * j + dst] += ax[p] * expected[3 * j + src];
            }
          }
        }
        ok = require_ok(transpose
                          ? kls_solve_transpose(solver, 2, rhs, 3, rhs, 3)
                          : kls_solve(solver, 2, rhs, 3, rhs, 3),
                        "solve scaled update in place");
        for (int i = 0; ok && i < 6; ++i) {
          if (!close_enough(rhs[i], expected[i])) {
            fprintf(stderr, "scaled update mismatch mode=%d update=%d row=%d\n",
                    mode, update, i);
            ok = 0;
          }
        }
      }
      kls_stats stats = {0};
      stats.struct_size = sizeof(stats);
      if (ok) ok = require_ok(kls_get_stats(solver, &stats), "scaled-update stats");
      if (ok && (!removed_shape_stats_are_zero(&stats) ||
                 stats.refactor_last_supernode_pipeline_tasks != 0 ||
                 stats.refactor_last_supernode_pipeline_columns != 0 ||
                 stats.refactor_supernode_pipeline_task_count != 0 ||
                 stats.refactor_supernode_pipeline_column_count != 0)) {
        fprintf(stderr, "retired statistics became nonzero\n");
        ok = 0;
      }
    }
    kls_destroy(solver);
    if (!ok) return 0;
  }
  return 1;
}

static int test_large_sparse_low_degree_retained_tolerance(void) {
  const int32_t n = 150000;
  const int32_t block_size = 1000;
  const int32_t degree = 5;
  const int32_t nnz = degree * n;
  int32_t *ap = (int32_t *)malloc(((size_t)n + 1u) * sizeof(*ap));
  int32_t *ai = (int32_t *)malloc((size_t)nnz * sizeof(*ai));
  double *ax = (double *)malloc((size_t)nnz * sizeof(*ax));
  double *changed = (double *)malloc((size_t)nnz * sizeof(*changed));
  double *expected = (double *)malloc((size_t)n * sizeof(*expected));
  double *b = (double *)malloc((size_t)n * sizeof(*b));
  double *x = (double *)malloc((size_t)n * sizeof(*x));
  const char *saved_contract_disable_value =
    getenv("KLS_DISABLE_PROMOTED_TOLERANCE_L2_RECOVERY");
  char *saved_contract_disable = saved_contract_disable_value != NULL
    ? strdup(saved_contract_disable_value) : NULL;
  const int had_contract_disable = saved_contract_disable_value != NULL;
  kls_solver *solver = NULL;
  kls_solver *contract_solver = NULL;
  kls_solver *contract_csr_solver = NULL;
  kls_solver *scaled_contract_solver = NULL;
  kls_options options;
  kls_stats stats;
  int ok = 1;

  if (ap == NULL || ai == NULL || ax == NULL || changed == NULL ||
      expected == NULL || b == NULL || x == NULL ||
      (had_contract_disable && saved_contract_disable == NULL)) {
    free(ap);
    free(ai);
    free(ax);
    free(changed);
    free(expected);
    free(b);
    free(x);
    free(saved_contract_disable);
    return 0;
  }

  int32_t p = 0;
  for (int32_t col = 0; col < n; ++col) {
    const int32_t block_begin = (col / block_size) * block_size;
    const int32_t local_col = col - block_begin;
    const int32_t dominant_row = (local_col + 1) % block_size;
    const int32_t group_col = local_col % 10;
    int32_t rows[degree];
    ap[col] = p;
    for (int32_t k = 0; k < degree - 1; ++k) {
      rows[k] = (local_col + k + 1) % block_size;
    }
    /* A fifth permutation fixes one column in ten and rotates the other
       nine by five positions.  Together with the four cyclic shifts this
       gives every row and column degree five, exactly 10% diagonal coverage,
       and no duplicate entry. */
    rows[degree - 1] = group_col == 0
      ? local_col
      : local_col - group_col + 1 + ((group_col - 1 + 5) % 9);
    for (int32_t i = 1; i < degree; ++i) {
      const int32_t row = rows[i];
      int32_t j = i;
      while (j > 0 && rows[j - 1] > row) {
        rows[j] = rows[j - 1];
        --j;
      }
      rows[j] = row;
    }
    for (int32_t k = 0; k < degree; ++k) {
      ai[p] = block_begin + rows[k];
      ax[p] = rows[k] == dominant_row
        ? 4.0 : 0.01 * (double)(1 + ((col + rows[k]) % 7));
      changed[p] = ax[p] *
        (1.0 + 1.0e-5 * (double)((p % 11) - 5));
      ++p;
    }
    expected[col] = 0.75 + 0.015625 * (double)(col % 17);
  }
  ap[n] = p;
  if (p != nnz) {
    fprintf(stderr, "unexpected low-degree policy nnz: %d/%d\n", p, nnz);
    ok = 0;
  }

  if (ok && !require_ok(kls_create(&solver),
                        "create low-degree retained tolerance")) {
    ok = 0;
  }
  memset(&stats, 0, sizeof(stats));
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(solver, &stats),
                        "build stats low-degree retained tolerance")) {
    ok = 0;
  }
  if (ok && stats.build_has_metis == 0) {
    goto cleanup;
  }
  kls_default_options(&options);
  options.threads = 4;

  /* Ordinary entrywise updates in this independent low-degree family select
     1e-4 below the requested 1e-3 threshold.  That retained numeric state is
     the generic accuracy capability: it must enter the relative-L2 contract
     without depending on the public benchmark's dimensions, BTF geometry,
     ordering tuple, or eight-worker count. */
  if (ok && unsetenv("KLS_DISABLE_PROMOTED_TOLERANCE_L2_RECOVERY") != 0) {
    perror("configure promoted-tolerance L2 contract");
    ok = 0;
  }
  if (ok &&
      (!require_ok(kls_create(&contract_solver),
                   "create promoted-tolerance L2 contract") ||
       !require_ok(kls_analyze_csc(contract_solver, KLS_INDEX_INT32, n,
                                   ap, ai, 0, &options),
                   "analyze promoted-tolerance L2 contract") ||
       !require_ok(kls_factor(contract_solver, ax),
                   "factor promoted-tolerance L2 contract") ||
       !require_ok(kls_refactor(contract_solver, changed),
                   "refactor promoted-tolerance L2 contract"))) {
    ok = 0;
  }
  if (ok) {
    memset(b, 0, (size_t)n * sizeof(*b));
    memset(x, 0, (size_t)n * sizeof(*x));
    for (int32_t col = 0; col < n; ++col) {
      for (int32_t entry = ap[col]; entry < ap[col + 1]; ++entry) {
        b[ai[entry]] += changed[entry] * expected[col];
      }
    }
    if (!require_ok(kls_solve(contract_solver, 1, b, 0, x, 0),
                    "solve promoted-tolerance L2 contract")) {
      ok = 0;
    }
  }
  memset(&stats, 0, sizeof(stats));
  stats.struct_size = sizeof(stats);
  if (ok &&
      (!require_ok(kls_get_stats(contract_solver, &stats),
                   "stats promoted-tolerance L2 contract") ||
       fabs(stats.selected_pivot_tolerance - 1.0e-4) > 1.0e-12 ||
       stats.promoted_tolerance_l2_recovery_eligible != 1 ||
       stats.promoted_tolerance_l2_contract_run_count < 1 ||
       stats.promoted_tolerance_l2_recovery_count != 0)) {
    fprintf(stderr,
            "unexpected promoted-tolerance L2 contract: tol=%.17g"
            " eligible=%d runs=%" PRId64 " recoveries=%" PRId64 "\n",
            stats.selected_pivot_tolerance,
            stats.promoted_tolerance_l2_recovery_eligible,
            stats.promoted_tolerance_l2_contract_run_count,
            stats.promoted_tolerance_l2_recovery_count);
    ok = 0;
  }
  for (int32_t row = 0; ok && row < n; ++row) {
    if (!close_enough(x[row], expected[row])) {
      fprintf(stderr,
              "promoted-tolerance L2 solve mismatch at %d: %.17g vs %.17g\n",
              row, x[row], expected[row]);
      ok = 0;
    }
  }

  /* Scaling and matching change the internal numeric frame, but not the
     public relative-L2 contract.  A promoted pivot tolerance in that frame
     must run the same measured residual gate; factor-cycle recovery remains
     deliberately ineligible because its input reconstruction is plain-only. */
  options.scale = 2;
  if (ok &&
      (!require_ok(kls_create(&scaled_contract_solver),
                   "create scaled promoted-tolerance L2 contract") ||
       !require_ok(kls_analyze_csc(scaled_contract_solver, KLS_INDEX_INT32, n,
                                   ap, ai, 0, &options),
                   "analyze scaled promoted-tolerance L2 contract") ||
       !require_ok(kls_factor(scaled_contract_solver, ax),
                   "factor scaled promoted-tolerance L2 contract") ||
       !require_ok(kls_refactor(scaled_contract_solver, changed),
                   "refactor scaled promoted-tolerance L2 contract"))) {
    ok = 0;
  }
  if (ok) {
    memset(x, 0, (size_t)n * sizeof(*x));
    if (!require_ok(kls_solve(scaled_contract_solver, 1, b, 0, x, 0),
                    "solve scaled promoted-tolerance L2 contract")) {
      ok = 0;
    }
  }
  memset(&stats, 0, sizeof(stats));
  stats.struct_size = sizeof(stats);
  if (ok &&
      (!require_ok(kls_get_stats(scaled_contract_solver, &stats),
                   "stats scaled promoted-tolerance L2 contract") ||
       !(stats.selected_pivot_tolerance > 0.0) ||
       !(stats.selected_pivot_tolerance < options.pivot_tolerance) ||
       stats.promoted_tolerance_l2_recovery_eligible != 0 ||
       stats.promoted_tolerance_l2_contract_run_count < 1)) {
    fprintf(stderr,
            "unexpected scaled promoted-tolerance L2 contract:"
            " scale=%d tol=%.17g eligible=%d runs=%" PRId64 "\n",
            stats.selected_scale, stats.selected_pivot_tolerance,
            stats.promoted_tolerance_l2_recovery_eligible,
            stats.promoted_tolerance_l2_contract_run_count);
    ok = 0;
  }
  for (int32_t row = 0; ok && row < n; ++row) {
    if (!close_enough(x[row], expected[row])) {
      fprintf(stderr,
              "scaled promoted-tolerance solve mismatch at %d: %.17g"
              " vs %.17g\n",
              row, x[row], expected[row]);
      ok = 0;
    }
  }
  options.scale = KLS_SCALE_AUTO;

  /* The compressed arrays are also CSR(A^T).  AUTO adopts the transpose
     orientation, so a public transpose solve exercises the same internal
     untransposed L2/recovery path through a nonidentity input-to-CSC map. */
  if (ok &&
      (!require_ok(kls_create(&contract_csr_solver),
                   "create CSR promoted-tolerance L2 contract") ||
       !require_ok(kls_analyze_csr(contract_csr_solver, KLS_INDEX_INT32, n,
                                   ap, ai, 0, &options),
                   "analyze CSR promoted-tolerance L2 contract") ||
       !require_ok(kls_factor(contract_csr_solver, ax),
                   "factor CSR promoted-tolerance L2 contract") ||
       !require_ok(kls_refactor(contract_csr_solver, changed),
                   "refactor CSR promoted-tolerance L2 contract"))) {
    ok = 0;
  }
  if (ok) {
    memset(x, 0, (size_t)n * sizeof(*x));
    if (!require_ok(kls_solve_transpose(
                      contract_csr_solver, 1, b, 0, x, 0),
                    "CSR promoted-tolerance L2 transpose solve")) {
      ok = 0;
    }
  }
  memset(&stats, 0, sizeof(stats));
  stats.struct_size = sizeof(stats);
  if (ok &&
      (!require_ok(kls_get_stats(contract_csr_solver, &stats),
                   "stats CSR promoted-tolerance L2 contract") ||
       fabs(stats.selected_pivot_tolerance - 1.0e-4) > 1.0e-12 ||
       stats.promoted_tolerance_l2_recovery_eligible != 1 ||
       stats.promoted_tolerance_l2_contract_run_count < 1)) {
    fprintf(stderr,
            "unexpected CSR promoted-tolerance L2 contract:"
            " orientation=%d tol=%.17g eligible=%d runs=%" PRId64 "\n",
            (int)stats.selected_orientation,
            stats.selected_pivot_tolerance,
            stats.promoted_tolerance_l2_recovery_eligible,
            stats.promoted_tolerance_l2_contract_run_count);
    ok = 0;
  }
  for (int32_t row = 0; ok && row < n; ++row) {
    if (!close_enough(x[row], expected[row])) {
      fprintf(stderr,
              "CSR promoted-tolerance solve mismatch at %d: %.17g"
              " vs %.17g\n",
              row, x[row], expected[row]);
      ok = 0;
    }
  }
  memset(&stats, 0, sizeof(stats));
  stats.struct_size = sizeof(stats);
  if (ok && !require_ok(kls_get_stats(contract_solver, &stats),
                        "pre-disable promoted-tolerance L2 stats")) {
    ok = 0;
  }
  const int64_t contract_runs =
    stats.promoted_tolerance_l2_contract_run_count;
  if (ok && setenv("KLS_DISABLE_PROMOTED_TOLERANCE_L2_RECOVERY", "1", 1) != 0) {
    perror("disable promoted-tolerance L2 contract");
    ok = 0;
  }
  memset(&stats, 0, sizeof(stats));
  stats.struct_size = sizeof(stats);
  if (ok &&
      (!require_ok(kls_get_stats(contract_solver, &stats),
                   "disabled promoted-tolerance L2 stats") ||
       stats.promoted_tolerance_l2_recovery_eligible != 0)) {
    fprintf(stderr, "promoted-tolerance L2 disable was not honored\n");
    ok = 0;
  }
  if (ok) {
    memset(x, 0, (size_t)n * sizeof(*x));
    if (!require_ok(kls_solve(contract_solver, 1, b, 0, x, 0),
                    "disabled promoted-tolerance L2 solve")) {
      ok = 0;
    }
  }
  memset(&stats, 0, sizeof(stats));
  stats.struct_size = sizeof(stats);
  if (ok &&
      (!require_ok(kls_get_stats(contract_solver, &stats),
                   "disabled promoted-tolerance L2 run stats") ||
       stats.promoted_tolerance_l2_contract_run_count != contract_runs)) {
    fprintf(stderr, "disabled promoted-tolerance L2 solve entered contract\n");
    ok = 0;
  }
cleanup:
  kls_destroy(scaled_contract_solver);
  kls_destroy(contract_csr_solver);
  kls_destroy(contract_solver);
  if (!restore_env_value("KLS_DISABLE_PROMOTED_TOLERANCE_L2_RECOVERY",
                         had_contract_disable,
                         saved_contract_disable != NULL
                           ? saved_contract_disable : "")) {
    ok = 0;
  }
  kls_destroy(solver);
  free(ap);
  free(ai);
  free(ax);
  free(changed);
  free(expected);
  free(b);
  free(x);
  free(saved_contract_disable);
  return ok;
}

int main(void) {
#ifdef __linux__
  caller_affinity_valid =
    sched_getaffinity(0, sizeof(caller_affinity), &caller_affinity) == 0;
#endif
  if (!test_scaled_update_refactor()) {
    return EXIT_FAILURE;
  }
  if (!run_sn_panel_factor_test()) {
    fprintf(stderr, "sn panel factor test failed\n");
    return EXIT_FAILURE;
  }

  if (!test_exact_unchanged_refactor_reuse()) {
    return EXIT_FAILURE;
  }
  if (!test_repeated_factor_low_rcond_contract()) {
    return EXIT_FAILURE;
  }
  if (!test_large_sparse_low_degree_retained_tolerance()) {
    return EXIT_FAILURE;
  }
  if (!test_refactor_solve_api()) {
    return EXIT_FAILURE;
  }
  if (!test_tiny_singleton_btf_solve()) {
    return EXIT_FAILURE;
  }
  /* The remaining smoke cases deliberately assert individual refactor-engine
     counters.  Keep their fixtures exercising those engines; exact-reuse has
     its own changing/in-place coverage above. */
  if (setenv("KLS_DISABLE_UNCHANGED_REFACTOR", "1", 1) != 0) {
    perror("setenv KLS_DISABLE_UNCHANGED_REFACTOR=1");
    return EXIT_FAILURE;
  }
  if (!test_deferred_sort_rebuilds_compact_solve()) {
    return EXIT_FAILURE;
  }

  if (!test_csc()) {
    return EXIT_FAILURE;
  }
  if (!test_csr_and_refactor()) {
    return EXIT_FAILURE;
  }
  if (!test_serial_backend()) {
    return EXIT_FAILURE;
  }
  if (!test_sparse_diagonal_auto_scale()) {
    return EXIT_FAILURE;
  }
  if (!test_transposed_low_work_btf_native_solve()) {
    return EXIT_FAILURE;
  }
  if (!test_compact_singleton_run_solve()) {
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
  /* The following cases validate the restart kernels themselves.  Select the
     generic ETree repair engine explicitly so the measured default policy is
     free to choose a full rebuild on these deliberately tiny fixtures. */
  if (setenv("KLS_FAST_REJECT_REPAIR_MODE", "etree", 1) != 0) {
    perror("setenv KLS_FAST_REJECT_REPAIR_MODE=etree");
    return EXIT_FAILURE;
  }
  if (!test_fast_factor_pivot_check_fallback()) {
    return EXIT_FAILURE;
  }
  if (!test_fast_factor_rowwise_u_pivot_reject()) {
    return EXIT_FAILURE;
  }
  if (!test_fast_factor_root_independent_tail_restart()) {
    return EXIT_FAILURE;
  }
  if (!test_fast_factor_noncontiguous_tail_gap_work_bounds()) {
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
  if (!test_nonroot_tail_refreshes_preserved_suffix()) {
    return EXIT_FAILURE;
  }
  if (!test_scaled_fast_factor_block_restart()) {
    return EXIT_FAILURE;
  }
  if (!test_scaled_fast_factor_prefix_tail_restart()) {
    return EXIT_FAILURE;
  }
  if (!test_btf_fast_factor_block_restart()) {
    return EXIT_FAILURE;
  }
  if (!test_btf_checked_row_tail_scope()) {
    return EXIT_FAILURE;
  }
  if (!test_btf_prefix_tail_restart_with_offblock()) {
    return EXIT_FAILURE;
  }
  if (!test_parallel_btf_suffix_after_prefix_tail_restart()) {
    return EXIT_FAILURE;
  }
  if (!test_btf_row_refactor_offblock_refresh()) {
    return EXIT_FAILURE;
  }
  if (!test_verified_rhs_reuse_contract()) {
    return EXIT_FAILURE;
  }
  if (!test_btf_singleton_row_refactor_pattern()) {
    return EXIT_FAILURE;
  }
  if (!test_btf_singleton_row_refactor_batch_groups()) {
    return EXIT_FAILURE;
  }
  if (!test_checked_row_dense_prefix_scatter_tail_restart()) {
    return EXIT_FAILURE;
  }
  if (!test_unchecked_row_dense_compact_panel()) {
    return EXIT_FAILURE;
  }
  if (!test_row_dense_compact_panel_scalar_update()) {
    return EXIT_FAILURE;
  }
  if (!test_unchecked_row_dense_native_direct_input()) {
    return EXIT_FAILURE;
  }
  if (!test_unchecked_row_generic_target_direct_input()) {
    return EXIT_FAILURE;
  }
  if (!test_unchecked_row_sparse_segment_direct_input()) {
    return EXIT_FAILURE;
  }
  if (!test_checked_row_prefactor_finished_dependency()) {
    return EXIT_FAILURE;
  }
  if (!test_partial_supernode_non_dominant_tail()) {
    return EXIT_FAILURE;
  }
  if (!test_partial_compact_supernode_prefix_pipeline()) {
    return EXIT_FAILURE;
  }
  if (!test_batched_compact_supernode_update_probe()) {
    return EXIT_FAILURE;
  }
  if (!test_batched_compact_supernode_subrange_update_probe()) {
    return EXIT_FAILURE;
  }
  if (!test_egraph_retired_supernode_switch()) {
    return EXIT_FAILURE;
  }
  if (!test_ragged_batched_compact_supernode_update_probe()) {
    return EXIT_FAILURE;
  }
  if (!test_batch_group_ragged_supernode_update_probe()) {
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
  if (!test_parallel_row_refactor_full_ready_queue()) {
    return EXIT_FAILURE;
  }
  if (!test_checked_separator_flop_ready_queue()) {
    return EXIT_FAILURE;
  }
  if (!test_btf_duplicate_separator_forest()) {
    return EXIT_FAILURE;
  }
  if (!test_scaled_row_refactor_single_block()) {
    return EXIT_FAILURE;
  }
  if (!test_experimental_kls_first_factor()) {
    return EXIT_FAILURE;
  }
  if (!test_row_refactor_off_blocks_auto_kls_first_refactor()) {
    return EXIT_FAILURE;
  }
  if (!test_default_keeps_kls_first_factor_off()) {
    return EXIT_FAILURE;
  }
  if (!test_kls_first_separator_queue_plan()) {
    return EXIT_FAILURE;
  }
  if (!test_fast_factor_separator_queue_repair()) {
    return EXIT_FAILURE;
  }
  if (!test_kls_first_separator_pipeline_pivot_epoch()) {
    return EXIT_FAILURE;
  }
  if (!test_auto_kls_first_skips_scaled_single_block()) {
    return EXIT_FAILURE;
  }
  if (!test_forced_kls_first_scaled_single_block_separator_queue()) {
    return EXIT_FAILURE;
  }
  if (!test_pre_static_replays_kls_first_factor()) {
    return EXIT_FAILURE;
  }
  if (!test_experimental_row_uplooking_first_factor()) {
    return EXIT_FAILURE;
  }
  if (!test_experimental_row_uplooking_first_consumer_panel()) {
    return EXIT_FAILURE;
  }
  if (!test_experimental_row_uplooking_natural_pipeline()) {
    return EXIT_FAILURE;
  }
  if (!test_experimental_row_uplooking_lazy_panel_prefix()) {
    return EXIT_FAILURE;
  }
  if (!test_experimental_row_uplooking_dynamic_column_pivot()) {
    return EXIT_FAILURE;
  }
  if (!test_experimental_row_uplooking_btf_blocks()) {
    return EXIT_FAILURE;
  }
  if (!test_parallel_btf_row_supernode_update()) {
    return EXIT_FAILURE;
  }
  if (!test_kls_first_parallel_pivoted_btf_fallback()) {
    return EXIT_FAILURE;
  }
  if (!test_structurally_full_rank_singular_rejection()) {
    return EXIT_FAILURE;
  }
  if (!test_row_solve_from_numeric_after_klu_first()) {
    return EXIT_FAILURE;
  }
  if (!test_transpose_row_solve_from_numeric_parallel(0) ||
      !test_transpose_row_solve_from_numeric_parallel(1)) {
    return EXIT_FAILURE;
  }
  if (!test_kls_first_reseeds_after_pivot_repair()) {
    return EXIT_FAILURE;
  }
  if (unsetenv("KLS_FAST_REJECT_REPAIR_MODE") != 0) {
    perror("unsetenv KLS_FAST_REJECT_REPAIR_MODE");
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
