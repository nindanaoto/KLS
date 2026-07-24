#ifndef KLS_KLS_H
#define KLS_KLS_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct kls_solver kls_solver;

typedef enum kls_status {
  KLS_OK = 0,
  KLS_ERR_INVALID_ARGUMENT = -1,
  KLS_ERR_OUT_OF_MEMORY = -2,
  KLS_ERR_ANALYZE_FAILED = -3,
  KLS_ERR_FACTOR_FAILED = -4,
  KLS_ERR_REFACTOR_FAILED = -5,
  KLS_ERR_SOLVE_FAILED = -6,
  KLS_ERR_SINGULAR = -7,
  KLS_ERR_UNSUPPORTED = -8
} kls_status;

typedef enum kls_index_type {
  KLS_INDEX_INT32 = 4,
  KLS_INDEX_INT64 = 8
} kls_index_type;

typedef enum kls_ordering {
  KLS_ORDERING_AUTO = 0,
  KLS_ORDERING_AMD = 1,
  KLS_ORDERING_COLAMD = 2,
  KLS_ORDERING_NATURAL = 3,
  KLS_ORDERING_METIS = 4,
  KLS_ORDERING_SCOTCH = 5,
  KLS_ORDERING_AMF = 6
} kls_ordering;

typedef enum kls_orientation {
  KLS_ORIENTATION_AUTO = 0,
  KLS_ORIENTATION_NORMAL = 1,
  KLS_ORIENTATION_TRANSPOSE = 2
} kls_orientation;

typedef enum kls_backend {
  /* Existing adaptive KLS policy.  AUTO currently aliases KLS so adding
     an opt-in serial policy cannot change parallel behavior. */
  KLS_BACKEND_AUTO = 0,
  KLS_BACKEND_KLS = 1,
  /* Lean, one-thread policy: direct orientation/order analysis, a single
     pivoting factor, and KLS's mapped refactor when its invariants hold. */
  KLS_BACKEND_SERIAL = 2
} kls_backend;

typedef enum kls_factor_path {
  KLS_FACTOR_PATH_NONE = 0,
  KLS_FACTOR_PATH_KLU_FIRST = 1,
  KLS_FACTOR_PATH_KLS_FAST_REFACTOR = 2,
  KLS_FACTOR_PATH_KLU_FALLBACK = 3,
  KLS_FACTOR_PATH_PRESTATIC_KLU_FIRST = 4,
  KLS_FACTOR_PATH_KLS_FIRST = 5,
  KLS_FACTOR_PATH_PREDICTED_FIRST = 6,
  KLS_FACTOR_PATH_SERIAL = 7
} kls_factor_path;

typedef enum kls_refactor_path {
  KLS_REFACTOR_PATH_NONE = 0,
  KLS_REFACTOR_PATH_ROW = 1,
  KLS_REFACTOR_PATH_EGRAPH = 2,
  KLS_REFACTOR_PATH_MAPPED = 3,
  KLS_REFACTOR_PATH_POOL = 4,
  KLS_REFACTOR_PATH_KLU = 5,
  KLS_REFACTOR_PATH_SNB = 6,
  KLS_REFACTOR_PATH_UNCHANGED = 7,
  /* The new values pass the Dr*A*Dc consistency check.  The retained factor
     is reused and the two diagonal maps are applied at the solve boundary. */
  KLS_REFACTOR_PATH_DIAGONAL_EQUIVALENT = 8,
  /* Only changed diagonal BTF blocks were numerically refreshed; unchanged
     block factors were retained and changed off-diagonal couplings copied. */
  KLS_REFACTOR_PATH_PARTIAL_BTF = 9,
  /* A bounded small update retained the preceding numeric as a
     preconditioner; solves are certified against the new values by
     residual-driven iterative refinement. */
  KLS_REFACTOR_PATH_RETAINED_PRECONDITIONER = 10
} kls_refactor_path;

#define KLS_SCALE_AUTO (-2)

/* Written by kls_default_options.  It occupies padding in the original
   options layout, allowing the appended backend field to be distinguished
   from uninitialized tail padding in binaries built against older headers. */
#define KLS_OPTIONS_ABI_VERSION UINT32_C(0x4b4c5302)

typedef enum kls_fast_reject_refresh_state {
  KLS_FAST_REJECT_REFRESH_UNKNOWN = 0,
  KLS_FAST_REJECT_REFRESH_PREFIX = 1,
  KLS_FAST_REJECT_REFRESH_ALL = 2
} kls_fast_reject_refresh_state;

typedef struct kls_options {
  size_t struct_size;
  int threads;
  kls_ordering ordering;
  kls_orientation orientation;
  int use_btf;
  int scale;
  uint32_t abi_version;
  double pivot_tolerance;
  double memory_growth;
  int halt_if_singular;
  int fast_factor;
  int static_pivoting;
  kls_backend backend;
} kls_options;

typedef struct kls_stats {
  size_t struct_size;
  int64_t n;
  int64_t nnz;
  int64_t nblocks;
  int64_t max_block;
  int64_t structural_rank;
  int64_t numerical_rank;
  int64_t singular_col;
  int64_t offdiag_pivots;
  int64_t reallocations;
  int64_t nnz_l;
  int64_t nnz_u;
  double analysis_seconds;
  double factor_seconds;
  double refactor_seconds;
  double solve_seconds;
  double estimated_flops;
  double factor_flops;
  double rcond;
  double condest;
  double rgrowth;
  size_t memory_bytes;
  size_t memory_peak_bytes;
  kls_ordering selected_ordering;
  kls_orientation selected_orientation;
  int last_kernel_status;
  int selected_scale;
  int selected_btf;
  double selected_pivot_tolerance;
  int selected_static_pivoting;
  int64_t fast_rejected_pivot;
  int64_t fast_rejected_pivot_col;
  int64_t fast_rejected_row;
  double fast_rejected_multiplier_abs;
  double fast_rejected_pivot_abs;
  double fast_rejected_candidate_abs;
  int64_t fast_rejected_tail_candidate_row;
  double fast_rejected_tail_candidate_abs;
  int64_t fast_rejected_tail_candidate_count;
  int64_t fast_rejected_tail_candidate_position;
  int fast_rejected_tail_repair_ready;
  int64_t fast_repaired_pivot_row;
  int fast_repaired_pivot_matches_tail_candidate;
  int64_t fast_repaired_first_changed_pivot;
  int64_t fast_repaired_prefix_changed_pivots;
  int64_t fast_repaired_suffix_changed_pivots;
  int fast_repaired_tail_restart_ready;
  double fast_repaired_block_work;
  int64_t fast_repaired_tail_restart_columns;
  double fast_repaired_tail_restart_work;
  double fast_repaired_tail_restart_saved_work;
  int selected_exact_matching;
  int selected_spral_matching;
  int fast_block_restarts;
  int fast_tail_restarts;
  int64_t refactor_dependency_levels;
  int64_t refactor_dependency_max_width;
  int64_t refactor_dependency_edges;
  int64_t refactor_dependency_cluster_levels;
  int64_t refactor_dependency_pipeline_columns;
  double refactor_dependency_work;
  double refactor_dependency_pipeline_work;
  double refactor_stream_dependency_entries;
  double refactor_stream_pivot_entries;
  double refactor_stream_output_entries;
  int64_t fast_rejected_block_start;
  int64_t fast_rejected_block_size;
  int64_t fast_rejected_suffix_columns;
  int64_t fast_rejected_descendant_columns;
  double fast_rejected_descendant_work;
  int64_t fast_rejected_row_tail_columns;
  double fast_rejected_row_tail_work;
  int64_t fast_rejected_etree_columns;
  double fast_rejected_etree_work;
  int64_t fast_rejected_pivoting_tail_columns;
  double fast_rejected_pivoting_tail_work;
  int64_t fast_rejected_pivoting_tail_first;
  int64_t fast_rejected_pivoting_tail_last;
  int fast_rejected_pivoting_tail_contains_reject;
  int fast_rejected_pivoting_tail_topological;
  int fast_rejected_refresh_state;
  int fast_factor_fail_reason;
  int fast_factor_fail_status;
  int64_t refactor_dependency_root_columns;
  int64_t refactor_dependency_leaf_columns;
  int64_t refactor_dependency_max_fanout;
  double refactor_dependency_max_column_work;
  double refactor_dependency_pipeline_max_column_work;
  int64_t refactor_supernode_candidate_count;
  int64_t refactor_supernode_candidate_rows;
  int64_t refactor_supernode_candidate_max_width;
  double refactor_supernode_candidate_dense_entries;
  double refactor_supernode_candidate_trailing_entries;
  int64_t refactor_supernode_consumer_run_count;
  int64_t refactor_supernode_consumer_run_rows;
  int64_t refactor_supernode_consumer_run_max_width;
  int64_t refactor_supernode_consumer_suffix_count;
  double refactor_supernode_consumer_l_entries;
  double refactor_supernode_consumer_internal_entries;
  int64_t refactor_supernode_consumer_panel_count;
  int64_t refactor_supernode_consumer_reused_panel_count;
  int64_t refactor_supernode_consumer_reused_run_count;
  int64_t refactor_supernode_consumer_reused_run_rows;
  int64_t refactor_supernode_consumer_panel_max_runs;
  int64_t refactor_supernode_consumer_panel_max_rows;
  double refactor_supernode_consumer_reused_l_entries;
  double refactor_supernode_consumer_reused_internal_entries;
  int64_t row_refactor_group_count;
  int64_t row_refactor_group_level_count;
  int64_t row_refactor_group_level_max_width;
  int64_t row_refactor_segment_count;
  int64_t row_refactor_segment_rows;
  int64_t row_refactor_segment_max_width;
  double row_refactor_segment_dense_entries;
  double row_refactor_segment_trailing_entries;
  int64_t row_refactor_dense_segment_count;
  int64_t row_refactor_dense_segment_rows;
  int64_t row_refactor_dense_segment_max_width;
  double row_refactor_dense_segment_dense_entries;
  double row_refactor_dense_segment_trailing_entries;
  int64_t row_refactor_group_cluster_levels;
  int64_t row_refactor_group_pipeline_groups;
  int64_t row_refactor_group_pipeline_rows;
  double row_refactor_group_pipeline_work;
  int row_refactor_last_run;
  int row_refactor_last_checked;
  int row_refactor_last_parallel;
  int64_t row_refactor_run_count;
  int64_t row_refactor_checked_run_count;
  int64_t row_refactor_parallel_run_count;
  int64_t row_refactor_group_dependency_edges;
  int64_t row_refactor_group_root_count;
  int64_t row_refactor_group_leaf_count;
  int64_t row_refactor_group_max_fanout;
  int row_refactor_last_ready_queue;
  int64_t row_refactor_ready_queue_run_count;
  int64_t row_refactor_ready_queue_group_count;
  int row_refactor_last_done_bitmap;
  int64_t row_refactor_done_bitmap_run_count;
  int64_t row_refactor_input_cleanup_rows;
  int64_t row_refactor_input_cleanup_entries;
  int row_refactor_last_defer_value_scatter;
  int64_t row_refactor_defer_value_scatter_run_count;
  int row_refactor_last_work_ready_queue;
  int64_t row_refactor_work_ready_queue_run_count;
  int64_t row_refactor_ready_queue_workspace_groups;
  int64_t fast_rejected_group_tail_groups;
  int64_t fast_rejected_group_tail_rows;
  double fast_rejected_group_tail_work;
  int64_t fast_rejected_pivoting_tail_seed_columns;
  int row_refactor_values_dirty;
  int row_refactor_last_lazy_value_scatter;
  int64_t row_refactor_lazy_value_scatter_run_count;
  int row_refactor_last_row_solve;
  int64_t row_refactor_row_solve_run_count;
  int64_t row_solve_parallel_run_count;
  int64_t row_solve_parallel_l_slice_runs;
  int64_t row_solve_parallel_u_slice_runs;
  int64_t row_solve_parallel_l_sparse_level_runs;
  int64_t row_solve_parallel_u_sparse_level_runs;
  int64_t row_solve_thread_count;
  int64_t row_solve_l_thread_max_rect_entries;
  int64_t row_solve_u_thread_max_rect_entries;
  int row_solve_partition_ready;
  int64_t row_solve_partition_slices;
  int64_t row_solve_l_sparse_level_count;
  int64_t row_solve_l_sparse_cluster_levels;
  int64_t row_solve_l_sparse_level_max_width;
  int64_t row_solve_l_dense_tail_start;
  int64_t row_solve_l_dense_tail_rows;
  int64_t row_solve_l_dense_tail_entries;
  int64_t row_solve_l_slice_max_entries;
  int64_t row_solve_l_segmented_rows;
  int64_t row_solve_l_rect_entries;
  int64_t row_solve_l_tri_entries;
  int64_t row_solve_u_sparse_level_count;
  int64_t row_solve_u_sparse_cluster_levels;
  int64_t row_solve_u_sparse_level_max_width;
  int64_t row_solve_u_dense_tail_start;
  int64_t row_solve_u_dense_tail_rows;
  int64_t row_solve_u_dense_tail_entries;
  int64_t row_solve_u_slice_max_entries;
  int64_t row_solve_u_segmented_rows;
  int64_t row_solve_u_rect_entries;
  int64_t row_solve_u_tri_entries;
  int fast_repaired_last_offdiag_suffix_refresh;
  int64_t fast_repaired_offdiag_suffix_refresh_count;
  int64_t fast_repaired_offdiag_full_refresh_count;
  int64_t fast_repaired_tail_restart_overcompute_columns;
  double fast_repaired_tail_restart_overcompute_work;
  int64_t fast_repaired_tail_restart_skipped_columns;
  double fast_repaired_tail_restart_skipped_work;
  int fast_repaired_tail_restart_exact_mask;
  int64_t row_refactor_last_local_ready_groups;
  int64_t row_refactor_local_ready_group_count;
  kls_factor_path last_factor_path;
  kls_refactor_path last_refactor_path;
  int64_t factor_etree_block_start;
  int64_t factor_etree_block_size;
  int64_t factor_etree_levels;
  int64_t factor_etree_max_width;
  int64_t factor_etree_edges;
  int64_t factor_etree_root_columns;
  int64_t factor_etree_leaf_columns;
  int64_t factor_etree_max_fanout;
  int fast_kls_block_restarts;
  int fast_kls_rebuild_restarts;
  int fast_kls_block_restart_last_row_pipeline;
  int64_t fast_kls_block_restart_row_pipeline_count;
  int64_t fast_kls_block_restart_last_row_pipeline_rows;
  int64_t fast_kls_block_restart_last_row_pipeline_threads;
  int64_t fast_kls_block_restart_last_row_pipeline_prefix_rows;
  int64_t fast_kls_block_restart_last_row_pipeline_suffix_rows;
  int64_t fast_kls_block_restart_last_row_pipeline_pivot_tail_rows;
  int64_t fast_kls_block_restart_last_row_pipeline_pivot_restarts;
  double row_refactor_total_group_work;
  int row_refactor_auto_enabled;
  int row_refactor_auto_values_ready;
  int row_refactor_auto_work_allowed;
  int row_refactor_auto_should_run;
  int build_has_metis;
  int build_has_scotch;
  int build_has_spral_scaling;
  int fast_rejected_pivoting_tail_contiguous;
  int fast_rejected_pivoting_tail_suffix_exact;
  int64_t fast_rejected_pivoting_tail_gap_columns;
  int64_t fast_rejected_pivoting_tail_suffix_overcompute_columns;
  double fast_rejected_pivoting_tail_suffix_overcompute_work;
  int row_refactor_last_compact_dense_panel;
  int64_t row_refactor_compact_dense_panel_count;
  int64_t row_refactor_compact_dense_panel_eligible_count;
  int64_t row_refactor_compact_dense_panel_eligible_rows;
  double row_refactor_compact_dense_panel_update_work;
  double row_refactor_compact_dense_panel_entries;
  int64_t fast_repaired_parallel_tail_blocks;
  int64_t row_refactor_compact_dense_panel_persistent_groups;
  int64_t row_refactor_compact_dense_panel_persistent_entries;
  int row_refactor_last_compact_dense_panel_persistent;
  int64_t row_refactor_compact_dense_panel_persistent_run_count;
  int64_t row_refactor_last_private_ready_groups;
  int64_t row_refactor_private_ready_group_count;
  int64_t kls_tail_last_mapped_columns;
  int64_t kls_tail_mapped_column_count;
  int64_t kls_first_last_row_uplooking_columns;
  int64_t kls_first_row_uplooking_column_count;
  int kls_first_last_dominant_btf_pipeline;
  int64_t kls_first_dominant_btf_pipeline_count;
  int64_t kls_first_last_dominant_btf_pipeline_block;
  int64_t kls_first_last_dominant_btf_pipeline_rows;
  int kls_first_last_dominant_btf_pipeline_has_separator;
  int64_t kls_first_dominant_btf_pipeline_without_separator_count;
  int64_t kls_first_last_row_refactor_seeded_rows;
  int64_t kls_first_row_refactor_seeded_row_count;
  int kls_first_last_row_pipeline;
  int64_t kls_first_row_pipeline_run_count;
  int64_t kls_first_last_row_pipeline_rows;
  int64_t kls_first_last_row_pipeline_threads;
  int kls_first_last_row_pipeline_partial;
  int64_t kls_first_row_pipeline_partial_run_count;
  int64_t kls_first_last_row_pipeline_partial_rows;
  int64_t kls_first_last_row_pipeline_partial_threads;
  int kls_first_last_row_pipeline_pivot_tail;
  int64_t kls_first_row_pipeline_pivot_tail_run_count;
  int64_t kls_first_last_row_pipeline_pivot_tail_rows;
  int64_t kls_first_last_row_pipeline_pivot_restarts;
  int64_t kls_first_row_pipeline_pivot_restart_count;
  int64_t kls_first_last_row_pipeline_pivot_serial_rows;
  int kls_first_last_row_pipeline_prefix_panel_rebuild;
  int64_t kls_first_row_pipeline_prefix_panel_rebuild_count;
  int64_t kls_first_last_row_pipeline_prefix_panel_rebuild_rows;
  int64_t kls_first_active_rank_pivot_reset_count;
  int64_t kls_first_active_rank_pivot_reset_rows;
  int64_t kls_first_active_rank_pivot_panel_rebuild_count;
  int64_t kls_first_active_rank_pivot_panel_rebuild_rows;
  int64_t kls_first_last_dynamic_column_pivots;
  int64_t kls_first_dynamic_column_pivot_count;
  int64_t separator_analyzed_rows;
  int64_t separator_thread_count;
  int64_t separator_component_count;
  int64_t separator_private_components;
  int64_t separator_pipeline_components;
  int64_t separator_private_rows;
  int64_t separator_pipeline_rows;
  int64_t separator_private_max_rows;
  int64_t separator_pipeline_max_rows;
  int row_refactor_last_separator_private_queue;
  int64_t row_refactor_separator_private_queue_run_count;
  int64_t row_refactor_last_separator_private_components;
  int64_t row_refactor_separator_private_component_count;
  int64_t kls_first_last_separator_dynamic_column_pivots;
  int64_t kls_first_separator_dynamic_column_pivot_count;
  int64_t kls_first_last_separator_dynamic_column_fallbacks;
  int64_t kls_first_separator_dynamic_column_fallback_count;
  int64_t separator_global_begin;
  int64_t separator_global_end;
  int64_t fast_rejected_prefix_refresh_columns;
  int64_t fast_rejected_prefix_refresh_count;
  int row_refactor_last_partial_supernode_pipeline;
  int64_t row_refactor_last_partial_supernode_pipeline_groups;
  int64_t row_refactor_last_partial_supernode_pipeline_rows;
  int64_t row_refactor_partial_supernode_pipeline_run_count;
  int row_refactor_last_compact_supernode_update;
  int64_t row_refactor_compact_supernode_update_count;
  int64_t row_refactor_compact_supernode_update_rows;
  int64_t row_refactor_compact_supernode_update_entries;
  int row_refactor_last_compact_supernode_gemv;
  int64_t row_refactor_compact_supernode_gemv_count;
  int64_t row_refactor_compact_supernode_gemv_rows;
  int64_t row_refactor_compact_supernode_gemv_entries;
  int row_refactor_last_compact_supernode_trsv;
  int64_t row_refactor_compact_supernode_trsv_count;
  int64_t row_refactor_compact_supernode_trsv_rows;
  int64_t row_refactor_compact_supernode_trsv_entries;
  int row_refactor_last_compact_supernode_batch;
  int64_t row_refactor_compact_supernode_batch_count;
  int64_t row_refactor_compact_supernode_batch_rows;
  int64_t row_refactor_compact_supernode_batch_dep_rows;
  int64_t row_refactor_compact_supernode_batch_entries;
  int64_t row_refactor_compact_supernode_batch_pattern_count;
  int64_t row_refactor_compact_supernode_batch_pattern_rows;
  int64_t row_refactor_compact_supernode_batch_candidate_count;
  int64_t row_refactor_compact_supernode_batch_candidate_rows;
  int64_t row_refactor_compact_supernode_batch_candidate_dep_rows;
  int64_t row_refactor_compact_supernode_batch_rejected_work_count;
  int build_has_cblas;
  int row_refactor_last_separator_flop_queue;
  int64_t row_refactor_separator_flop_queue_run_count;
  int row_refactor_last_separator_flop_ordered_private;
  int64_t row_refactor_separator_flop_ordered_private_run_count;
  int64_t row_refactor_last_separator_flop_components;
  int64_t row_refactor_separator_flop_component_count;
  int64_t row_refactor_last_separator_flop_private_groups;
  int64_t row_refactor_last_separator_flop_pipeline_groups;
  int64_t row_refactor_separator_flop_private_group_count;
  int64_t row_refactor_separator_flop_pipeline_group_count;
  int64_t row_refactor_last_compact_dense_panel_direct_input_rows;
  int64_t row_refactor_compact_dense_panel_direct_input_rows;
  int64_t row_refactor_last_dense_segment_direct_input_rows;
  int64_t row_refactor_dense_segment_direct_input_rows;
  int64_t row_refactor_last_batch_direct_input_rows;
  int64_t row_refactor_batch_direct_input_rows;
  int64_t kls_first_last_separator_extent_dynamic_column_pivots;
  int64_t kls_first_separator_extent_dynamic_column_pivot_count;
  int row_refactor_last_compact_supernode_partial_update;
  int64_t row_refactor_compact_supernode_partial_update_count;
  int64_t row_refactor_compact_supernode_partial_update_rows;
  int64_t row_refactor_compact_supernode_partial_update_entries;
  int64_t row_refactor_last_separator_flop_closure_groups;
  int64_t row_refactor_separator_flop_closure_group_count;
  int kls_first_auto_skipped_scaled_single_block;
  int64_t kls_first_auto_skipped_scaled_single_block_count;
  int64_t refactor_last_supernode_pipeline_tasks;
  int64_t refactor_last_supernode_pipeline_columns;
  int64_t refactor_supernode_pipeline_task_count;
  int64_t refactor_supernode_pipeline_column_count;
  int64_t refactor_last_supernode_update_runs;
  int64_t refactor_last_supernode_update_rows;
  int64_t refactor_last_supernode_update_entries;
  int64_t refactor_supernode_update_run_count;
  int64_t refactor_supernode_update_rows;
  int64_t refactor_supernode_update_entries;
  int64_t refactor_last_ready_queue_columns;
  int64_t refactor_ready_queue_run_count;
  double parallel_model_r1;
  double parallel_model_r2;
  int parallel_model_recommends_parallel;
  int row_refactor_auto_model_recommended;
  int row_refactor_auto_model_attempted;
  int row_refactor_auto_model_accepted;
  int64_t refactor_last_supernode_cblas_update_runs;
  int64_t refactor_last_supernode_cblas_update_rows;
  int64_t refactor_last_supernode_cblas_update_entries;
  int64_t refactor_supernode_cblas_update_run_count;
  int64_t refactor_supernode_cblas_update_rows;
  int64_t refactor_supernode_cblas_update_entries;
  int64_t refactor_last_supernode_blocked_update_runs;
  int64_t refactor_last_supernode_blocked_update_rows;
  int64_t refactor_last_supernode_blocked_update_entries;
  int64_t refactor_supernode_blocked_update_run_count;
  int64_t refactor_supernode_blocked_update_rows;
  int64_t refactor_supernode_blocked_update_entries;
  int64_t refactor_last_supernode_cached_probe_attempts;
  int64_t refactor_last_supernode_cached_probe_panel_hits;
  int64_t refactor_last_supernode_cached_probe_contiguous;
  int64_t refactor_last_supernode_cached_probe_allowed;
  int64_t refactor_last_supernode_cached_probe_allowed_rows;
  int64_t refactor_last_supernode_cached_probe_applied;
  int64_t refactor_last_supernode_cached_probe_applied_rows;
  int64_t refactor_supernode_cached_probe_attempt_count;
  int64_t refactor_supernode_cached_probe_panel_hits;
  int64_t refactor_supernode_cached_probe_contiguous;
  int64_t refactor_supernode_cached_probe_allowed;
  int64_t refactor_supernode_cached_probe_allowed_rows;
  int64_t refactor_supernode_cached_probe_applied;
  int64_t refactor_supernode_cached_probe_applied_rows;
  int refactor_supernode_cached_probe_disabled;
  int64_t refactor_supernode_cached_probe_disable_count;
  int refactor_supernode_update_disabled;
  int64_t refactor_supernode_update_disable_count;
  int64_t row_refactor_dense_producer_run_count;
  int64_t row_refactor_dense_producer_run_rows;
  int64_t row_refactor_dense_producer_run_dep_rows;
  int64_t row_refactor_dense_producer_run_max_per_row;
  int64_t row_refactor_dense_producer_full_suffix_run_count;
  int64_t row_refactor_dense_producer_full_suffix_rows;
  int64_t row_refactor_dense_producer_multi_run_rows;
  int64_t row_refactor_dense_producer_fragmented_rows;
  int64_t row_refactor_dense_producer_target_count;
  int64_t row_refactor_dense_producer_target_none_count;
  int64_t row_refactor_dense_producer_target_external_count;
  int64_t row_refactor_dense_producer_target_dense_count;
  int64_t row_refactor_dense_producer_target_pivot_count;
  int64_t row_refactor_dense_producer_target_trailing_count;
  int64_t kls_first_last_parallel_btf_blocks;
  int64_t kls_first_parallel_btf_block_count;
  int64_t kls_first_last_separator_dynamic_column_rejects;
  int64_t kls_first_separator_dynamic_column_reject_count;
  int kls_first_last_separator_queue;
  int64_t kls_first_separator_queue_run_count;
  int64_t kls_first_last_separator_queue_private_components;
  int64_t kls_first_last_separator_queue_pipeline_components;
  int64_t kls_first_last_separator_queue_private_rows;
  int64_t kls_first_last_separator_queue_pipeline_rows;
  int64_t kls_first_last_separator_queue_nonempty_threads;
  int64_t kls_first_last_separator_queue_max_thread_rows;
  int kls_first_last_separator_queue_executed;
  int64_t kls_first_separator_queue_executed_run_count;
  int64_t kls_first_last_separator_queue_executed_private_rows;
  int64_t kls_first_last_separator_queue_executed_pipeline_rows;
  int kls_first_last_separator_queue_parallel_private;
  int64_t kls_first_separator_queue_parallel_private_run_count;
  int64_t kls_first_last_separator_queue_parallel_private_rows;
  int64_t kls_first_last_separator_queue_parallel_private_threads;
  int kls_first_last_separator_queue_parallel_pipeline;
  int64_t kls_first_separator_queue_parallel_pipeline_run_count;
  int64_t kls_first_last_separator_queue_parallel_pipeline_rows;
  int64_t kls_first_last_separator_queue_parallel_pipeline_threads;
  int kls_first_last_separator_queue_pipeline_partial;
  int64_t kls_first_separator_queue_pipeline_partial_run_count;
  int64_t kls_first_last_separator_queue_pipeline_partial_rows;
  int64_t kls_first_last_separator_queue_pipeline_partial_threads;
  int kls_first_last_separator_queue_pipeline_wait_partial;
  int64_t kls_first_separator_queue_pipeline_wait_partial_run_count;
  int64_t kls_first_last_separator_queue_pipeline_wait_partial_rows;
  int64_t kls_first_last_separator_queue_pipeline_wait_partial_deps;
  int kls_first_last_separator_queue_pipeline_supernode_update;
  int64_t kls_first_separator_queue_pipeline_supernode_update_run_count;
  int64_t kls_first_last_separator_queue_pipeline_supernode_update_groups;
  int64_t kls_first_last_separator_queue_pipeline_supernode_update_rows;
  int kls_first_last_separator_queue_pipeline_pivot_tail;
  int64_t kls_first_separator_queue_pipeline_pivot_tail_run_count;
  int64_t kls_first_last_separator_queue_pipeline_pivot_tail_rows;
  int64_t kls_first_last_separator_queue_pipeline_pivot_restarts;
  int64_t kls_first_separator_queue_pipeline_pivot_restart_count;
  int64_t kls_first_last_separator_queue_pipeline_pivot_serial_rows;
  int kls_first_last_separator_queue_pipeline_prefix_panel_rebuild;
  int64_t kls_first_separator_queue_pipeline_prefix_panel_rebuild_count;
  int64_t kls_first_last_separator_queue_pipeline_prefix_panel_rebuild_rows;
  int64_t fast_rejected_pivoting_tail_row_seed_columns;
  int64_t row_refactor_last_sparse_segment_direct_input_rows;
  int64_t row_refactor_sparse_segment_direct_input_rows;
  int64_t row_refactor_segment_input_target_rows;
  int64_t row_refactor_segment_input_target_entries;
  int64_t row_refactor_last_segment_target_input_rows;
  int64_t row_refactor_segment_target_input_rows;
  int64_t row_refactor_segment_input_cleanup_rows;
  int64_t row_refactor_segment_input_cleanup_entries;
  int64_t row_refactor_last_segment_target_cleanup_rows;
  int64_t row_refactor_segment_target_cleanup_rows;
  int64_t row_refactor_last_segment_target_cleanup_entries;
  int64_t row_refactor_segment_target_cleanup_entries;
  int64_t row_refactor_last_compact_panel_solve_values;
  int64_t row_refactor_compact_panel_solve_values;
  int64_t row_refactor_last_compact_panel_group_solve_rows;
  int64_t row_refactor_compact_panel_group_solve_rows;
  int64_t row_refactor_last_compact_panel_group_solve_entries;
  int64_t row_refactor_compact_panel_group_solve_entries;
  int64_t row_refactor_last_compact_panel_scalar_update_rows;
  int64_t row_refactor_compact_panel_scalar_update_rows;
  int64_t row_refactor_last_compact_panel_scalar_update_entries;
  int64_t row_refactor_compact_panel_scalar_update_entries;
  int row_refactor_last_prefactor;
  int64_t row_refactor_last_prefactor_rows;
  int64_t row_refactor_last_prefactor_deps;
  int64_t row_refactor_prefactor_run_count;
  int64_t row_refactor_prefactor_rows;
  int64_t row_refactor_prefactor_deps;
  int row_refactor_last_prefactor_supernode;
  int64_t row_refactor_last_prefactor_supernode_rows;
  int64_t row_refactor_last_prefactor_supernode_deps;
  int64_t row_refactor_prefactor_supernode_run_count;
  int64_t row_refactor_prefactor_supernode_rows;
  int64_t row_refactor_prefactor_supernode_deps;
  int kls_first_last_row_supernode_update;
  int64_t kls_first_row_supernode_update_run_count;
  int64_t kls_first_row_supernode_update_groups;
  int64_t kls_first_row_supernode_update_rows;
  int64_t kls_first_last_row_supernode_update_groups;
  int64_t kls_first_last_row_supernode_update_rows;
  int kls_first_last_separator_queue_pipeline_supernode_panel_update;
  int64_t kls_first_separator_queue_pipeline_supernode_panel_update_run_count;
  int64_t kls_first_last_separator_queue_pipeline_supernode_panel_update_groups;
  int64_t kls_first_last_separator_queue_pipeline_supernode_panel_update_rows;
  int kls_first_last_row_supernode_panel_update;
  int64_t kls_first_row_supernode_panel_update_run_count;
  int64_t kls_first_row_supernode_panel_update_groups;
  int64_t kls_first_row_supernode_panel_update_rows;
  int64_t kls_first_last_row_supernode_panel_update_groups;
  int64_t kls_first_last_row_supernode_panel_update_rows;
  int64_t fast_rejected_pivoting_tail_block_seed_columns;
  int64_t fast_kls_block_restart_last_row_pipeline_gap_rows;
  int64_t fast_kls_block_restart_last_row_pipeline_supernode_update_groups;
  int64_t fast_kls_block_restart_last_row_pipeline_supernode_update_rows;
  int64_t
    fast_kls_block_restart_last_row_pipeline_supernode_panel_update_groups;
  int64_t fast_kls_block_restart_last_row_pipeline_supernode_panel_update_rows;
  int row_refactor_last_compact_dense_panel_blocked;
  int64_t row_refactor_compact_dense_panel_blocked_run_count;
  int64_t row_refactor_compact_dense_panel_blocked_rows;
  int64_t row_refactor_compact_dense_panel_blocked_entries;
  double kls_first_last_separator_queue_min_thread_work;
  double kls_first_last_separator_queue_max_thread_work;
  int selected_exact_matching_scaling;
  int64_t fast_rejected_pivoting_tail_etree_edges;
  int64_t fast_rejected_pivoting_tail_etree_roots;
  int64_t fast_rejected_pivoting_tail_etree_leaves;
  int64_t fast_rejected_pivoting_tail_etree_max_fanout;
  int fast_repaired_tail_restart_etree_mask;
  int64_t fast_rejected_pivoting_tail_etree_levels;
  int64_t fast_rejected_pivoting_tail_etree_max_width;
  int64_t parallel_task_flow_threads;
  int64_t parallel_task_flow_dependencies;
  double parallel_task_flow_work;
  double parallel_task_flow_finish_time;
  double parallel_task_flow_speedup;
  int parallel_task_flow_recommends_parallel;
  int kls_first_last_separator_queue_partitioned;
  int64_t kls_first_separator_queue_partitioned_count;
  int64_t kls_first_last_separator_queue_split_components;
  int fast_kls_block_restart_last_row_pipeline_etree_tail;
  int64_t fast_kls_block_restart_row_pipeline_etree_tail_count;
  int64_t fast_kls_block_restart_last_row_pipeline_etree_tail_rows;
  int64_t fast_kls_block_restart_last_row_pipeline_etree_tail_gap_rows;
  int fast_kls_block_restart_last_row_pipeline_etree_tail_exact_mask;
  int fast_kls_block_restart_last_row_pipeline_separator_tail_scope;
  int64_t fast_kls_block_restart_row_pipeline_separator_tail_scope_count;
  int64_t fast_kls_block_restart_last_row_pipeline_separator_tail_scope_rows;
  int fast_kls_block_restart_last_row_pipeline_separator_queue;
  int64_t fast_kls_block_restart_row_pipeline_separator_queue_count;
  int64_t fast_kls_block_restart_last_row_pipeline_separator_private_rows;
  int64_t fast_kls_block_restart_last_row_pipeline_separator_pipeline_rows;
  int64_t fast_kls_block_restart_last_row_pipeline_separator_private_threads;
  int fast_kls_block_restart_last_row_pipeline_separator_partitioned;
  int64_t fast_kls_block_restart_last_row_pipeline_separator_split_components;
  int fast_kls_block_restart_last_row_pipeline_etree_ready;
  int64_t fast_kls_block_restart_row_pipeline_etree_ready_count;
  int64_t fast_kls_block_restart_last_row_pipeline_etree_ready_rows;
  int64_t fast_kls_block_restart_last_row_pipeline_etree_ready_threads;
  int fast_kls_block_restart_last_row_pipeline_etree_prefactor;
  int64_t fast_kls_block_restart_row_pipeline_etree_prefactor_count;
  int64_t fast_kls_block_restart_last_row_pipeline_etree_prefactor_rows;
  int64_t fast_kls_block_restart_last_row_pipeline_etree_prefactor_threads;
  int64_t fast_kls_block_restart_last_row_pipeline_etree_prefactor_wait_rows;
  int64_t fast_kls_block_restart_last_row_pipeline_etree_prefactor_wait_deps;
  int row_refactor_native_row_panel_enabled;
  int row_refactor_last_native_row_panel;
  int row_refactor_native_row_panel_auto_disabled;
  int64_t row_refactor_native_row_panel_count;
  int64_t row_refactor_native_row_panel_rows;
  int64_t row_refactor_native_row_panel_entries;
  int64_t row_refactor_native_row_panel_blocked_count;
  int64_t row_refactor_native_row_panel_blocked_rows;
  int64_t row_refactor_native_row_panel_blocked_entries;
  int64_t row_refactor_native_row_panel_fallback_count;
  int64_t row_refactor_native_row_panel_checked_reject_count;
  int64_t row_refactor_native_row_panel_auto_disable_count;
  int64_t refactor_l_pattern_columns;
  int64_t refactor_l_pattern_entries;
  int64_t refactor_l_adjacent_run_count;
  int64_t refactor_l_adjacent_run_entries;
  int64_t refactor_l_adjacent_run_max_len;
  int64_t refactor_l_contiguous_suffix_columns;
  int64_t refactor_l_contiguous_suffix_entries;
  int64_t refactor_l_contiguous_suffix_max_len;
  int internal_index_bytes;
  int refactor_map_index32_enabled;
  int64_t refactor_map_index32_entries;
  int refactor_l_index32_enabled;
  int64_t refactor_l_index32_entries;
  int refactor_u_index32_enabled;
  int64_t refactor_u_index32_entries;
  int64_t kls_first_row_panel_cache_build_count;
  int64_t kls_first_row_panel_cache_build_panels;
  int64_t kls_first_row_panel_cache_build_entries;
  int64_t kls_first_row_panel_cache_append_count;
  int64_t kls_first_row_panel_cache_append_panels;
  int64_t kls_first_row_panel_cache_append_entries;
  double row_refactor_auto_lower_bound_work;
  int row_refactor_auto_lower_bound_rejected;
  int row_refactor_auto_pattern_build_failed;
  int row_refactor_auto_value_copy_failed;
  int64_t refactor_supernode_panel_count;
  int64_t refactor_supernode_panel_used_count;
  int64_t refactor_supernode_consumer_plan_cached_panel_count;
  int64_t refactor_supernode_consumer_plan_cached_panel_rows;
  int64_t refactor_supernode_consumer_plan_strict_cached_panel_count;
  int64_t refactor_supernode_consumer_plan_strict_cached_panel_rows;
  int64_t refactor_supernode_consumer_plan_strict_cached_panel_trailing_entries;
  int refactor_supernode_consumer_plan_group_l_built;
  int refactor_supernode_consumer_plan_group_l_storage_limited;
  int64_t refactor_supernode_consumer_plan_group_l_panel_count;
  int64_t refactor_supernode_consumer_plan_group_l_run_count;
  int64_t refactor_supernode_consumer_plan_group_l_rows;
  int64_t refactor_supernode_consumer_plan_group_l_run_rows;
  int64_t refactor_supernode_consumer_plan_group_l_exec_run_count;
  int64_t refactor_supernode_consumer_plan_group_l_exec_run_rows;
  int64_t refactor_supernode_consumer_plan_group_l_batch_candidate_count;
  int64_t refactor_supernode_consumer_plan_group_l_batch_candidate_run_count;
  int64_t refactor_supernode_consumer_plan_group_l_batch_candidate_run_rows;
  int64_t refactor_supernode_consumer_plan_group_l_batch_candidate_entries;
  int64_t refactor_supernode_consumer_plan_group_l_batch_candidate_max_runs;
  int64_t refactor_supernode_consumer_plan_group_l_batch_candidate_advance_count;
  int64_t refactor_supernode_consumer_plan_group_l_batch_candidate_advance_run_count;
  int64_t refactor_supernode_consumer_plan_group_l_batch_candidate_advance_run_rows;
  double refactor_supernode_consumer_plan_group_l_batch_candidate_advance_work;
  double refactor_supernode_consumer_plan_group_l_batch_candidate_advance_max_work;
  int64_t refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_count;
  int64_t refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_run_count;
  int64_t refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_run_rows;
  int64_t refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_entries;
  double refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_advance_work;
  double refactor_supernode_consumer_plan_group_l_batch_candidate_max_payoff_ratio;
  int64_t refactor_supernode_consumer_plan_group_l_dense_entries;
  int64_t refactor_supernode_consumer_plan_group_l_trailing_entries;
  int64_t refactor_supernode_consumer_plan_group_l_bytes;
  int64_t refactor_last_supernode_consumer_plan_group_l_dense_writes;
  int64_t refactor_last_supernode_consumer_plan_group_l_trailing_writes;
  int64_t refactor_last_supernode_consumer_plan_group_l_invalidations;
  int64_t refactor_supernode_consumer_plan_group_l_dense_write_count;
  int64_t refactor_supernode_consumer_plan_group_l_trailing_write_count;
  int64_t refactor_supernode_consumer_plan_group_l_invalidation_count;
  int64_t refactor_last_supernode_consumer_plan_group_l_update_runs;
  int64_t refactor_last_supernode_consumer_plan_group_l_update_rows;
  int64_t refactor_last_supernode_consumer_plan_group_l_update_entries;
  int64_t refactor_supernode_consumer_plan_group_l_update_run_count;
  int64_t refactor_supernode_consumer_plan_group_l_update_rows;
  int64_t refactor_supernode_consumer_plan_group_l_update_entries;
  int refactor_supernode_consumer_plan_shape_targets_built;
  int64_t refactor_supernode_consumer_plan_shape_target_group_count;
  int64_t refactor_supernode_consumer_plan_shape_target_run_count;
  int64_t refactor_supernode_consumer_plan_shape_target_entries;
  int64_t refactor_supernode_consumer_plan_shape_target_rows;
  int64_t refactor_supernode_consumer_plan_shape_target_max_rows;
  int64_t refactor_supernode_consumer_plan_deferred_columns;
  int64_t refactor_supernode_consumer_plan_deferred_entries;
  int64_t refactor_supernode_consumer_plan_deferred_unique_rows;
  int64_t refactor_supernode_consumer_plan_batch_deferred_columns;
  int64_t refactor_supernode_consumer_plan_batch_deferred_entries;
  int64_t refactor_supernode_consumer_plan_batch_deferred_unique_rows;
  int64_t row_refactor_group_single_count;
  int64_t row_refactor_group_batch_count;
  int64_t row_refactor_group_batch_rows;
  int64_t row_refactor_group_batch_max_width;
  int64_t row_refactor_group_batch_width_le_4_count;
  int64_t row_refactor_group_batch_width_le_8_count;
  int64_t row_refactor_group_scalar_candidate_count;
  int64_t row_refactor_group_scalar_candidate_rows;
  int64_t row_refactor_group_scalar_short_count;
  int64_t row_refactor_group_scalar_short_rows;
  int64_t row_refactor_group_scalar_stop_level_mismatch_count;
  int64_t row_refactor_group_scalar_stop_internal_dep_count;
  int64_t row_refactor_group_scalar_stop_next_segment_count;
  int64_t row_refactor_group_scalar_stop_max_width_count;
  int64_t row_refactor_group_scalar_stop_matrix_end_count;
  int64_t row_refactor_group_generic_count;
  int64_t row_refactor_group_generic_rows;
  int64_t row_refactor_group_generic_max_width;
  int64_t row_refactor_group_dense_count;
  int64_t row_refactor_group_dense_rows;
  int64_t row_refactor_group_dense_max_width;
  double row_refactor_group_single_work;
  double row_refactor_group_batch_work;
  double row_refactor_group_generic_work;
  double row_refactor_group_dense_work;
  int64_t row_refactor_last_separator_flop_private_threads;
  int64_t row_refactor_last_separator_flop_private_min_groups;
  int64_t row_refactor_last_separator_flop_private_max_groups;
  double row_refactor_last_separator_flop_private_min_work;
  double row_refactor_last_separator_flop_private_max_work;
  int64_t refactor_u_supernode_pattern_count;
  int64_t refactor_u_supernode_pattern_rows;
  int64_t refactor_u_supernode_pattern_max_width;
  int64_t refactor_u_supernode_pattern_right_entries;
  double refactor_u_supernode_pattern_internal_entries;
  int64_t refactor_u_supernode_value_dense_entries;
  int64_t refactor_u_supernode_value_right_entries;
  int64_t refactor_last_u_supernode_value_dense_writes;
  int64_t refactor_last_u_supernode_value_right_writes;
  int64_t refactor_u_supernode_value_dense_write_count;
  int64_t refactor_u_supernode_value_right_write_count;
  int64_t refactor_u_supernode_l_panel_count;
  int64_t refactor_u_supernode_l_dense_entries;
  int64_t refactor_u_supernode_l_trailing_entries;
  int64_t refactor_u_supernode_l_prune_count;
  int64_t refactor_u_supernode_l_pruned_panels;
  int64_t refactor_u_supernode_l_pruned_dense_entries;
  int64_t refactor_u_supernode_l_pruned_trailing_entries;
  int64_t refactor_last_u_supernode_l_update_runs;
  int64_t refactor_last_u_supernode_l_update_rows;
  int64_t refactor_last_u_supernode_l_update_entries;
  int64_t refactor_last_u_supernode_l_probe_attempts;
  int64_t refactor_last_u_supernode_l_panel_misses;
  int64_t refactor_last_u_supernode_l_short_rejects;
  int64_t refactor_last_u_supernode_l_stream_rejects;
  int64_t refactor_last_u_supernode_l_work_rejects;
  int refactor_u_supernode_l_exec_disabled;
  int64_t refactor_u_supernode_l_exec_disable_count;
  int64_t refactor_u_supernode_l_update_run_count;
  int64_t refactor_u_supernode_l_update_rows;
  int64_t refactor_u_supernode_l_update_entries;
  int64_t refactor_supernode_consumer_plan_panel_count;
  int64_t refactor_supernode_consumer_plan_reused_panel_count;
  int64_t refactor_supernode_consumer_plan_run_count;
  int64_t refactor_supernode_consumer_plan_run_rows;
  int64_t refactor_supernode_consumer_plan_positioned_run_count;
  int64_t refactor_supernode_consumer_plan_positioned_run_rows;
  int64_t refactor_supernode_consumer_plan_max_panel_runs;
  int64_t refactor_supernode_consumer_plan_max_panel_rows;
  double refactor_supernode_consumer_plan_l_entries;
  double refactor_supernode_consumer_plan_internal_entries;
  int64_t refactor_supernode_consumer_plan_bytes;
  int64_t refactor_supernode_consumer_plan_small_run_count;
  int64_t refactor_supernode_consumer_plan_small_run_rows;
  int64_t refactor_supernode_consumer_plan_batch_panel_count;
  int64_t refactor_supernode_consumer_plan_batch_run_count;
  int64_t refactor_supernode_consumer_plan_batch_run_rows;
  int64_t refactor_supernode_consumer_plan_batch_small_run_count;
  int64_t refactor_supernode_consumer_plan_batch_small_run_rows;
  int64_t refactor_supernode_consumer_plan_column_count;
  int64_t refactor_supernode_consumer_plan_max_column_runs;
  int64_t refactor_supernode_consumer_plan_max_column_rows;
  int64_t refactor_supernode_consumer_plan_column_batch_count;
  int64_t refactor_supernode_consumer_plan_column_batch_run_count;
  int64_t refactor_supernode_consumer_plan_column_batch_run_rows;
  int64_t refactor_supernode_consumer_plan_column_batch_small_run_count;
  int64_t refactor_supernode_consumer_plan_column_batch_small_run_rows;
  int64_t refactor_supernode_consumer_plan_column_shape_batch_count;
  int64_t refactor_supernode_consumer_plan_column_shape_batch_run_count;
  int64_t refactor_supernode_consumer_plan_column_shape_batch_run_rows;
  int64_t refactor_supernode_consumer_plan_column_shape_batch_small_run_count;
  int64_t refactor_supernode_consumer_plan_column_shape_batch_small_run_rows;
  int64_t refactor_supernode_consumer_plan_column_shape_batch_max_runs;
  int64_t refactor_supernode_consumer_plan_shape_batch_count;
  int64_t refactor_supernode_consumer_plan_shape_batch_run_count;
  int64_t refactor_supernode_consumer_plan_shape_batch_run_rows;
  int64_t refactor_supernode_consumer_plan_shape_batch_small_run_count;
  int64_t refactor_supernode_consumer_plan_shape_batch_small_run_rows;
  int64_t refactor_supernode_consumer_plan_shape_batch_max_runs;
  int64_t refactor_supernode_consumer_plan_first_dep_shape_batch_count;
  int64_t refactor_supernode_consumer_plan_first_dep_shape_batch_run_count;
  int64_t refactor_supernode_consumer_plan_first_dep_shape_batch_run_rows;
  int64_t
    refactor_supernode_consumer_plan_first_dep_shape_batch_small_run_count;
  int64_t
    refactor_supernode_consumer_plan_first_dep_shape_batch_small_run_rows;
  int64_t refactor_supernode_consumer_plan_first_dep_shape_batch_max_runs;
  int64_t refactor_supernode_consumer_plan_shape_batch_advance_count;
  int64_t refactor_supernode_consumer_plan_shape_batch_advance_run_count;
  int64_t refactor_supernode_consumer_plan_shape_batch_advance_run_rows;
  int64_t refactor_supernode_consumer_plan_shape_batch_advance_dep_count;
  double refactor_supernode_consumer_plan_shape_batch_advance_work;
  int64_t refactor_supernode_consumer_plan_shape_batch_advance_max_deps;
  double refactor_supernode_consumer_plan_shape_batch_advance_max_work;
  int64_t refactor_last_supernode_consumer_plan_attempts;
  int64_t refactor_last_supernode_consumer_plan_hits;
  int64_t refactor_last_supernode_consumer_plan_applied;
  int64_t refactor_last_supernode_consumer_plan_rows;
  int64_t refactor_last_supernode_consumer_plan_entries;
  int64_t refactor_supernode_consumer_plan_attempt_count;
  int64_t refactor_supernode_consumer_plan_hit_count;
  int64_t refactor_supernode_consumer_plan_apply_count;
  int64_t refactor_supernode_consumer_plan_apply_rows;
  int64_t refactor_supernode_consumer_plan_apply_entries;
  int refactor_supernode_consumer_plan_exec_disabled;
  int64_t refactor_supernode_consumer_plan_exec_disable_count;
  int64_t refactor_last_supernode_consumer_plan_claimed_columns;
  int64_t refactor_last_supernode_consumer_plan_claim_skips;
  int64_t refactor_last_supernode_consumer_plan_claim_waits;
  int64_t refactor_supernode_consumer_plan_claimed_columns;
  int64_t refactor_supernode_consumer_plan_claim_skip_count;
  int64_t refactor_supernode_consumer_plan_claim_wait_count;
  int64_t refactor_last_supernode_cached_probe_shape_rejects;
  int64_t refactor_last_supernode_cached_probe_shape_reject_rows;
  int64_t refactor_last_supernode_cached_probe_stream_rejects;
  int64_t refactor_last_supernode_cached_probe_stream_reject_rows;
  int64_t refactor_last_supernode_cached_probe_work_rejects;
  int64_t refactor_last_supernode_cached_probe_work_reject_rows;
  int64_t refactor_last_supernode_cached_probe_workspace_rejects;
  int64_t refactor_last_supernode_cached_probe_workspace_reject_rows;
  int64_t refactor_supernode_cached_probe_shape_rejects;
  int64_t refactor_supernode_cached_probe_shape_reject_rows;
  int64_t refactor_supernode_cached_probe_stream_rejects;
  int64_t refactor_supernode_cached_probe_stream_reject_rows;
  int64_t refactor_supernode_cached_probe_work_rejects;
  int64_t refactor_supernode_cached_probe_work_reject_rows;
  int64_t refactor_supernode_cached_probe_workspace_rejects;
  int64_t refactor_supernode_cached_probe_workspace_reject_rows;
  int64_t refactor_supernode_consumer_plan_prefix_advance_batch_count;
  int64_t refactor_supernode_consumer_plan_prefix_advance_batch_run_count;
  int64_t refactor_supernode_consumer_plan_prefix_advance_batch_run_rows;
  int64_t refactor_supernode_consumer_plan_prefix_advance_batch_dep_count;
  double refactor_supernode_consumer_plan_prefix_advance_batch_work;
  int64_t refactor_supernode_consumer_plan_prefix_advance_batch_max_runs;
  int64_t refactor_supernode_consumer_plan_prefix_advance_batch_max_deps;
  double refactor_supernode_consumer_plan_prefix_advance_batch_max_work;
  int64_t refactor_supernode_consumer_plan_shape_bounded_advance_dep_limit;
  int64_t refactor_supernode_consumer_plan_shape_bounded_advance_count;
  int64_t refactor_supernode_consumer_plan_shape_bounded_advance_run_count;
  int64_t refactor_supernode_consumer_plan_shape_bounded_advance_run_rows;
  int64_t refactor_supernode_consumer_plan_shape_bounded_advance_dep_count;
  int64_t
    refactor_supernode_consumer_plan_shape_bounded_advance_update_entries;
  double refactor_supernode_consumer_plan_shape_bounded_advance_work;
  int64_t
    refactor_supernode_consumer_plan_shape_bounded_advance_max_run_deps;
  double refactor_supernode_consumer_plan_shape_bounded_advance_max_run_work;
  int64_t
    refactor_supernode_consumer_plan_shape_bounded_advance_payoff_count;
  int64_t
    refactor_supernode_consumer_plan_shape_bounded_advance_payoff_run_count;
  int64_t
    refactor_supernode_consumer_plan_shape_bounded_advance_payoff_run_rows;
  int64_t
    refactor_supernode_consumer_plan_shape_bounded_advance_payoff_update_entries;
  double
    refactor_supernode_consumer_plan_shape_bounded_advance_payoff_advance_work;
  double
    refactor_supernode_consumer_plan_shape_bounded_advance_max_payoff_ratio;
  int64_t refactor_supernode_algorithm5_large_panel_count;
  int64_t refactor_supernode_algorithm5_large_panel_rows;
  int64_t refactor_supernode_algorithm5_large_panel_prefix_rows;
  int64_t refactor_supernode_algorithm5_large_panel_max_width;
  int64_t refactor_supernode_algorithm5_candidate_run_count;
  int64_t refactor_supernode_algorithm5_candidate_run_rows;
  int64_t refactor_supernode_algorithm5_prefix_run_count;
  int64_t refactor_supernode_algorithm5_prefix_run_rows;
  double refactor_supernode_algorithm5_prefix_update_work;
  int64_t refactor_supernode_algorithm5_crossing_run_count;
  int64_t refactor_supernode_algorithm5_crossing_run_rows;
  int64_t refactor_supernode_algorithm5_prefix_advance_run_count;
  int64_t refactor_supernode_algorithm5_prefix_advance_run_rows;
  int64_t refactor_supernode_algorithm5_prefix_advance_dep_count;
  double refactor_supernode_algorithm5_prefix_advance_work;
  int64_t refactor_supernode_algorithm5_prefix_payoff_run_count;
  int64_t refactor_supernode_algorithm5_prefix_payoff_run_rows;
  double refactor_supernode_algorithm5_prefix_payoff_update_work;
  double refactor_supernode_algorithm5_prefix_payoff_advance_work;
  int64_t refactor_last_egraph_algorithm5_prefactor_columns;
  int64_t refactor_last_egraph_algorithm5_prefactor_deps;
  int64_t refactor_egraph_algorithm5_prefactor_column_count;
  int64_t refactor_egraph_algorithm5_prefactor_dep_count;
  int64_t refactor_supernode_algorithm5_prefix_panel_count;
  int64_t refactor_supernode_algorithm5_prefix_panel_run_count;
  int64_t refactor_supernode_algorithm5_prefix_panel_run_rows;
  int64_t refactor_supernode_algorithm5_prefix_panel_advance_dep_count;
  double refactor_supernode_algorithm5_prefix_panel_update_work;
  double refactor_supernode_algorithm5_prefix_panel_advance_work;
  int64_t refactor_supernode_algorithm5_prefix_panel_max_runs;
  int64_t refactor_supernode_algorithm5_prefix_panel_payoff_count;
  int64_t refactor_supernode_algorithm5_prefix_panel_payoff_run_count;
  int64_t refactor_supernode_algorithm5_prefix_panel_payoff_run_rows;
  double refactor_supernode_algorithm5_prefix_panel_payoff_update_work;
  double refactor_supernode_algorithm5_prefix_panel_payoff_advance_work;
  int64_t refactor_supernode_algorithm5_prefix_panel_payoff_subset_count;
  int64_t refactor_supernode_algorithm5_prefix_panel_payoff_subset_run_count;
  int64_t refactor_supernode_algorithm5_prefix_panel_payoff_subset_run_rows;
  double refactor_supernode_algorithm5_prefix_panel_payoff_subset_update_work;
  double refactor_supernode_algorithm5_prefix_panel_payoff_subset_advance_work;
  int64_t refactor_supernode_algorithm5_prefix_panel_payoff_subset_max_runs;
  int64_t refactor_supernode_algorithm5_prefix_panel_payoff_subset_current_count;
  int64_t refactor_supernode_algorithm5_prefix_panel_payoff_subset_multi_current_count;
  int64_t refactor_supernode_algorithm5_prefix_panel_payoff_subset_max_currents;
  int64_t refactor_supernode_algorithm5_payoff_group_count;
  int64_t refactor_supernode_algorithm5_payoff_group_prefix_rows;
  int64_t refactor_supernode_algorithm5_payoff_group_max_run_rows;
  int64_t refactor_supernode_algorithm5_payoff_group_positioned_runs;
  int64_t refactor_supernode_algorithm5_payoff_group_current_total;
  int64_t refactor_supernode_algorithm5_payoff_group_multi_current_count;
  int64_t refactor_supernode_algorithm5_payoff_group_max_currents;
  int64_t refactor_supernode_algorithm5_payoff_group_workspace_rows;
  int64_t refactor_supernode_algorithm5_payoff_group_max_workspace_rows;
  int64_t refactor_supernode_algorithm5_payoff_group_advance_deps;
  int64_t refactor_supernode_algorithm5_payoff_group_max_advance_deps;
  int64_t refactor_supernode_algorithm5_payoff_group_zero_advance_runs;
  int64_t refactor_supernode_algorithm5_payoff_group_advance_unique_deps;
  int64_t refactor_supernode_algorithm5_payoff_group_advance_duplicate_deps;
  int64_t refactor_supernode_algorithm5_payoff_group_advance_shared_deps;
  int64_t refactor_supernode_algorithm5_payoff_group_advance_max_dep_fanout;
  int64_t
    refactor_supernode_algorithm5_payoff_group_advance_duplicate_update_entries;
  int64_t refactor_supernode_algorithm5_payoff_group_suffix_deps;
  int64_t refactor_supernode_algorithm5_payoff_group_max_suffix_deps;
  int64_t refactor_supernode_algorithm5_payoff_group_suffix_update_entries;
  int64_t
    refactor_supernode_algorithm5_payoff_group_max_run_suffix_update_entries;
  int64_t refactor_supernode_algorithm5_payoff_group_suffix_unique_deps;
  int64_t refactor_supernode_algorithm5_payoff_group_suffix_duplicate_deps;
  int64_t refactor_supernode_algorithm5_payoff_group_suffix_shared_deps;
  int64_t refactor_supernode_algorithm5_payoff_group_suffix_max_dep_fanout;
  int64_t
    refactor_supernode_algorithm5_payoff_group_suffix_duplicate_update_entries;
  int64_t refactor_supernode_algorithm5_payoff_group_advance_slots;
  int64_t refactor_supernode_algorithm5_payoff_group_max_advance_slots;
  int64_t refactor_supernode_algorithm5_payoff_group_max_run_advance_slots;
  int64_t refactor_supernode_algorithm5_payoff_group_target_entries;
  int64_t refactor_supernode_algorithm5_payoff_group_max_target_entries;
  int64_t refactor_supernode_algorithm5_payoff_group_max_run_target_entries;
  int64_t refactor_supernode_algorithm5_payoff_group_target_slots;
  int64_t refactor_supernode_algorithm5_payoff_group_max_target_slots;
  int64_t refactor_supernode_algorithm5_payoff_group_max_run_target_slots;
  int64_t refactor_supernode_algorithm5_payoff_group_pattern_width;
  int64_t refactor_supernode_algorithm5_payoff_group_max_pattern_width;
  int64_t refactor_supernode_algorithm5_payoff_current_state_rows;
  int64_t refactor_supernode_algorithm5_payoff_current_state_max_rows;
  int64_t refactor_supernode_algorithm5_payoff_runtime_workspace_rows;
  int64_t refactor_supernode_algorithm5_payoff_runtime_target_slots;
  int64_t refactor_supernode_algorithm5_payoff_runtime_current_count;
  int64_t refactor_last_supernode_algorithm5_payoff_slot_accum_runs;
  int64_t refactor_last_supernode_algorithm5_payoff_slot_accum_rows;
  int64_t refactor_last_supernode_algorithm5_payoff_slot_accum_target_entries;
  int64_t refactor_last_supernode_algorithm5_payoff_slot_accum_target_slots;
  int64_t refactor_last_supernode_algorithm5_payoff_prefix_prep_runs;
  int64_t refactor_last_supernode_algorithm5_payoff_prefix_prep_rows;
  int64_t refactor_last_supernode_algorithm5_payoff_prefix_prep_target_entries;
  int64_t refactor_last_supernode_algorithm5_payoff_prefix_prep_target_slots;
  int64_t refactor_last_supernode_algorithm5_payoff_advance_seed_runs;
  int64_t refactor_last_supernode_algorithm5_payoff_advance_seed_deps;
  int64_t refactor_last_supernode_algorithm5_payoff_advance_seed_slots;
  int64_t refactor_last_supernode_algorithm5_payoff_current_state_seed_runs;
  int64_t refactor_last_supernode_algorithm5_payoff_current_state_seed_deps;
  int64_t refactor_last_supernode_algorithm5_payoff_current_state_seed_rows;
  int64_t refactor_last_supernode_algorithm5_payoff_final_trigger_batches;
  int64_t
    refactor_last_supernode_algorithm5_payoff_final_trigger_multi_batches;
  int64_t refactor_last_supernode_algorithm5_payoff_final_trigger_claims;
  int64_t refactor_last_supernode_algorithm5_payoff_final_trigger_suffix_deps;
  int64_t refactor_last_supernode_algorithm5_payoff_suffix_advance_slots;
  int64_t refactor_last_supernode_algorithm5_payoff_suffix_advance_deps;
  int64_t refactor_last_supernode_algorithm5_payoff_suffix_advance_updates;
  int64_t refactor_last_supernode_algorithm5_payoff_suffix_advance_finished;
  int refactor_supernode_consumer_plan_group_l_state_built;
  int refactor_supernode_consumer_plan_group_l_state_storage_limited;
  int64_t refactor_supernode_consumer_plan_group_l_state_run_count;
  int64_t refactor_supernode_consumer_plan_group_l_state_rows;
  int64_t refactor_supernode_consumer_plan_group_l_state_max_rows;
  int64_t refactor_supernode_consumer_plan_group_l_state_bytes;
  int refactor_supernode_consumer_plan_group_l_state_focus_enabled;
  int64_t refactor_supernode_consumer_plan_group_l_state_candidate_group_count;
  int64_t refactor_supernode_consumer_plan_group_l_state_candidate_run_count;
  int64_t refactor_supernode_consumer_plan_group_l_state_candidate_rows;
  int64_t refactor_supernode_consumer_plan_group_l_state_selected_group_count;
  int64_t refactor_supernode_consumer_plan_group_l_state_selected_run_count;
  int64_t refactor_supernode_consumer_plan_group_l_state_selected_rows;
  int64_t refactor_supernode_consumer_plan_group_l_state_selected_bytes;
  int64_t refactor_supernode_algorithm5_payoff_current_state_span_rows;
  int64_t refactor_supernode_algorithm5_payoff_current_state_max_span_rows;
  int64_t refactor_last_btf_scalar_run_candidates;
  int64_t refactor_last_btf_scalar_run_rows;
  int64_t refactor_last_btf_scalar_run_entries;
  int64_t refactor_last_btf_scalar_run_max_rows;
  int64_t refactor_btf_scalar_run_candidate_count;
  int64_t refactor_btf_scalar_run_rows;
  int64_t refactor_btf_scalar_run_entries;
  int64_t refactor_btf_scalar_run_max_rows;
  int64_t refactor_last_btf_scalar_run_exec_runs;
  int64_t refactor_last_btf_scalar_run_exec_rows;
  int64_t refactor_last_btf_scalar_run_exec_entries;
  int64_t refactor_last_btf_scalar_run_exec_max_rows;
  int64_t refactor_btf_scalar_run_exec_count;
  int64_t refactor_btf_scalar_run_exec_rows;
  int64_t refactor_btf_scalar_run_exec_entries;
  int64_t refactor_btf_scalar_run_exec_max_rows;
  int refactor_btf_scalar_run_group_built;
  int64_t refactor_btf_scalar_run_group_count;
  int64_t refactor_btf_scalar_run_group_current_total;
  int64_t refactor_btf_scalar_run_group_multi_count;
  int64_t refactor_btf_scalar_run_group_multi_current_total;
  int64_t refactor_btf_scalar_run_group_rows;
  int64_t refactor_btf_scalar_run_group_reused_rows;
  int64_t refactor_btf_scalar_run_group_entries;
  int64_t refactor_btf_scalar_run_group_reused_entries;
  int64_t refactor_btf_scalar_run_group_max_currents;
  int64_t refactor_btf_scalar_run_group_max_rows;
  int64_t refactor_btf_scalar_run_group_producer_step_count;
  int64_t refactor_btf_scalar_run_group_producer_step_multi_count;
  int64_t refactor_btf_scalar_run_group_producer_step_active_currents;
  int64_t refactor_btf_scalar_run_group_producer_step_unique_entries;
  int64_t refactor_btf_scalar_run_group_producer_step_duplicate_entries;
  int64_t refactor_btf_scalar_run_group_producer_step_reused_entries;
  int64_t refactor_btf_scalar_run_group_producer_step_max_currents;
  int64_t refactor_btf_scalar_run_group_producer_index_count;
  int64_t refactor_btf_scalar_run_group_producer_index_max_steps;
  int64_t refactor_btf_scalar_run_group_live_state_groups;
  int64_t refactor_btf_scalar_run_group_live_state_currents;
  int64_t refactor_btf_scalar_run_group_live_state_rows;
  int64_t refactor_btf_scalar_run_group_live_state_unique_rows;
  int64_t refactor_btf_scalar_run_group_live_state_reused_rows;
  int64_t refactor_btf_scalar_run_group_live_state_max_currents;
  int64_t refactor_btf_scalar_run_group_live_state_max_rows;
  int64_t refactor_btf_scalar_run_group_live_state_max_unique_rows;
  int64_t refactor_btf_scalar_run_group_live_step_count;
  int64_t refactor_btf_scalar_run_group_live_step_current_total;
  int64_t refactor_btf_scalar_run_group_live_step_rows;
  int64_t refactor_btf_scalar_run_group_live_step_unique_rows;
  int64_t refactor_btf_scalar_run_group_live_step_reused_rows;
  int64_t refactor_btf_scalar_run_group_live_step_max_currents;
  int64_t refactor_btf_scalar_run_group_live_step_max_rows;
  int64_t refactor_btf_scalar_run_group_live_step_max_unique_rows;
  int64_t refactor_btf_scalar_run_group_live_step_stored_rows;
  int64_t refactor_btf_scalar_run_group_live_step_storage_limited;
  int64_t refactor_btf_scalar_run_group_wake_count;
  int64_t refactor_btf_scalar_run_group_wake_member_total;
  int64_t refactor_btf_scalar_run_group_live_step_runtime_full_step_count;
  int64_t refactor_btf_scalar_run_group_live_step_runtime_full_current_count;
  int64_t refactor_btf_scalar_run_group_live_step_runtime_full_rows;
  int64_t refactor_btf_scalar_run_group_live_step_runtime_partial_step_count;
  int64_t
    refactor_btf_scalar_run_group_live_step_runtime_partial_current_count;
  int64_t refactor_btf_scalar_run_group_live_step_runtime_partial_member_count;
  int64_t refactor_btf_scalar_run_group_live_step_runtime_partial_rows;
  int64_t refactor_btf_scalar_run_group_live_step_runtime_partial_active_rows;
  int64_t
    refactor_btf_scalar_run_group_live_step_runtime_partial_active_unique_rows;
  int64_t refactor_last_btf_scalar_run_group_waits;
  int64_t refactor_last_btf_scalar_run_group_wait_rows;
  int64_t refactor_last_btf_scalar_run_group_wait_entries;
  int64_t refactor_last_btf_scalar_run_group_overlap_waits;
  int64_t refactor_last_btf_scalar_run_group_overlap_rows;
  int64_t refactor_last_btf_scalar_run_group_overlap_entries;
  int64_t refactor_last_btf_scalar_run_group_max_live;
  int64_t refactor_btf_scalar_run_group_wait_count;
  int64_t refactor_btf_scalar_run_group_wait_rows;
  int64_t refactor_btf_scalar_run_group_wait_entries;
  int64_t refactor_btf_scalar_run_group_overlap_count;
  int64_t refactor_btf_scalar_run_group_overlap_rows;
  int64_t refactor_btf_scalar_run_group_overlap_entries;
  int64_t refactor_btf_scalar_run_group_max_live;
  int64_t refactor_last_btf_scalar_run_group_claim_surface_triggers;
  int64_t refactor_last_btf_scalar_run_group_claim_surface_groups;
  int64_t refactor_last_btf_scalar_run_group_claim_surface_currents;
  int64_t refactor_last_btf_scalar_run_group_claim_triggers;
  int64_t refactor_last_btf_scalar_run_group_claim_groups;
  int64_t refactor_last_btf_scalar_run_group_claim_currents;
  int64_t refactor_last_btf_scalar_run_group_prefix_ready_currents;
  int64_t refactor_last_btf_scalar_run_group_prefix_ready_rows;
  int64_t refactor_last_btf_scalar_run_group_prefix_ready_entries;
  int64_t refactor_last_btf_scalar_run_group_run_ready_currents;
  int64_t refactor_last_btf_scalar_run_group_run_ready_rows;
  int64_t refactor_last_btf_scalar_run_group_run_ready_entries;
  int64_t refactor_last_btf_scalar_run_group_wake_armed_currents;
  int64_t refactor_last_btf_scalar_run_group_wake_armed_rows;
  int64_t refactor_last_btf_scalar_run_group_wake_armed_entries;
  int64_t refactor_last_btf_scalar_run_group_wake_ready_currents;
  int64_t refactor_last_btf_scalar_run_group_wake_ready_rows;
  int64_t refactor_last_btf_scalar_run_group_wake_ready_entries;
  int64_t refactor_btf_scalar_run_group_state_rows_total;
  int64_t refactor_btf_scalar_run_group_state_max_rows;
  int64_t refactor_btf_scalar_run_group_state_current_count;
  int64_t refactor_btf_scalar_run_group_state_best_skip_total;
  int64_t refactor_btf_scalar_run_group_state_max_best_skip;
  int64_t refactor_btf_scalar_run_group_state_guard_lower_bound_rows;
  int64_t refactor_btf_scalar_run_group_state_guard_lower_bound_rejected;
  int64_t refactor_last_btf_scalar_run_group_state_materialized_currents;
  int64_t refactor_last_btf_scalar_run_group_state_materialized_rows;
  int64_t refactor_last_btf_scalar_run_group_state_materialized_prefix_deps;
  int64_t refactor_last_btf_scalar_run_group_state_advanced_currents;
  int64_t refactor_last_btf_scalar_run_group_state_advanced_rows;
  int64_t refactor_last_btf_scalar_run_group_state_advanced_entries;
  int64_t refactor_last_btf_scalar_run_group_state_advance_batch_groups;
  int64_t refactor_last_btf_scalar_run_group_state_advance_batch_currents;
  int64_t
    refactor_last_btf_scalar_run_group_state_advance_batch_unique_entries;
  int64_t
    refactor_last_btf_scalar_run_group_state_advance_batch_duplicate_entries;
  int64_t
    refactor_last_btf_scalar_run_group_state_advance_batch_max_currents;
  int64_t refactor_last_btf_scalar_run_group_state_step_advance_triggers;
  int64_t refactor_last_btf_scalar_run_group_state_step_advance_steps;
  int64_t refactor_last_btf_scalar_run_group_state_step_advance_currents;
  int64_t refactor_last_btf_scalar_run_group_state_step_advance_rows;
  int64_t refactor_last_btf_scalar_run_group_state_step_advance_entries;
  int64_t
    refactor_last_btf_scalar_run_group_state_step_advance_ready_currents;
  int64_t refactor_last_btf_scalar_run_group_state_step_advance_rejects;
  int64_t refactor_last_btf_scalar_run_group_state_step_batch_steps;
  int64_t refactor_last_btf_scalar_run_group_state_step_batch_currents;
  int64_t refactor_last_btf_scalar_run_group_state_step_batch_entries;
  int64_t refactor_last_btf_scalar_run_group_state_step_batch_state_rows;
  int64_t
    refactor_last_btf_scalar_run_group_state_step_batch_unique_state_rows;
  int64_t refactor_last_btf_scalar_run_group_state_step_window_rounds;
  int64_t refactor_last_btf_scalar_run_group_state_step_window_currents;
  int64_t refactor_last_btf_scalar_run_group_state_step_window_entries;
  int64_t refactor_last_btf_scalar_run_group_live_step_runtime_full_steps;
  int64_t refactor_last_btf_scalar_run_group_live_step_runtime_full_currents;
  int64_t refactor_last_btf_scalar_run_group_live_step_runtime_full_rows;
  int64_t refactor_last_btf_scalar_run_group_live_step_runtime_partial_steps;
  int64_t refactor_last_btf_scalar_run_group_live_step_runtime_partial_currents;
  int64_t refactor_last_btf_scalar_run_group_live_step_runtime_partial_members;
  int64_t refactor_last_btf_scalar_run_group_live_step_runtime_partial_rows;
  int64_t
    refactor_last_btf_scalar_run_group_live_step_runtime_partial_active_rows;
  int64_t
    refactor_last_btf_scalar_run_group_live_step_runtime_partial_active_unique_rows;
  int64_t refactor_last_btf_scalar_run_group_state_rejects;
  int64_t refactor_last_btf_scalar_run_group_state_exec_currents;
  int64_t refactor_last_btf_scalar_run_group_state_exec_skipped_deps;
  int64_t refactor_last_btf_scalar_run_group_state_exec_restored_rows;
  int64_t refactor_last_btf_scalar_run_group_state_exec_terminal_currents;
  int64_t
    refactor_last_btf_scalar_run_group_state_exec_wake_terminal_currents;
  int64_t
    refactor_last_btf_scalar_run_group_state_exec_wake_terminal_candidate_currents;
  int64_t
    refactor_last_btf_scalar_run_group_state_exec_wake_terminal_nonpipeline_currents;
  int64_t
    refactor_last_btf_scalar_run_group_state_exec_wake_terminal_dependency_miss_currents;
  int64_t
    refactor_last_btf_scalar_run_group_state_exec_wake_terminal_claim_miss_currents;
  int64_t refactor_last_btf_scalar_run_group_state_exec_remaining_deps;
  int64_t refactor_last_btf_scalar_run_group_state_exec_remaining_entries;
  int64_t refactor_last_btf_scalar_run_group_state_exec_rejects;
  int64_t refactor_last_btf_scalar_run_group_state_exec_dispatch_bypass_currents;
  int64_t
    refactor_last_btf_scalar_run_group_state_exec_dispatch_bypass_ready_currents;
  int64_t refactor_last_btf_scalar_run_group_state_exec_not_ready_currents;
  int64_t refactor_last_btf_scalar_run_group_state_exec_owned_ready_currents;
  int64_t refactor_last_btf_scalar_run_group_state_exec_late_ready_currents;
  int64_t refactor_btf_scalar_run_group_claim_surface_trigger_count;
  int64_t refactor_btf_scalar_run_group_claim_surface_group_count;
  int64_t refactor_btf_scalar_run_group_claim_surface_current_count;
  int64_t refactor_btf_scalar_run_group_claim_trigger_count;
  int64_t refactor_btf_scalar_run_group_claim_group_count;
  int64_t refactor_btf_scalar_run_group_claim_current_count;
  int64_t refactor_btf_scalar_run_group_prefix_ready_current_count;
  int64_t refactor_btf_scalar_run_group_prefix_ready_rows;
  int64_t refactor_btf_scalar_run_group_prefix_ready_entries;
  int64_t refactor_btf_scalar_run_group_run_ready_current_count;
  int64_t refactor_btf_scalar_run_group_run_ready_rows;
  int64_t refactor_btf_scalar_run_group_run_ready_entries;
  int64_t refactor_btf_scalar_run_group_wake_armed_current_count;
  int64_t refactor_btf_scalar_run_group_wake_armed_rows;
  int64_t refactor_btf_scalar_run_group_wake_armed_entries;
  int64_t refactor_btf_scalar_run_group_wake_ready_current_count;
  int64_t refactor_btf_scalar_run_group_wake_ready_rows;
  int64_t refactor_btf_scalar_run_group_wake_ready_entries;
  int64_t refactor_btf_scalar_run_group_state_materialized_current_count;
  int64_t refactor_btf_scalar_run_group_state_materialized_rows;
  int64_t refactor_btf_scalar_run_group_state_materialized_prefix_deps;
  int64_t refactor_btf_scalar_run_group_state_advanced_current_count;
  int64_t refactor_btf_scalar_run_group_state_advanced_rows;
  int64_t refactor_btf_scalar_run_group_state_advanced_entries;
  int64_t refactor_btf_scalar_run_group_state_advance_batch_group_count;
  int64_t refactor_btf_scalar_run_group_state_advance_batch_current_count;
  int64_t refactor_btf_scalar_run_group_state_advance_batch_unique_entries;
  int64_t
    refactor_btf_scalar_run_group_state_advance_batch_duplicate_entries;
  int64_t refactor_btf_scalar_run_group_state_advance_batch_max_currents;
  int64_t refactor_btf_scalar_run_group_state_step_advance_trigger_count;
  int64_t refactor_btf_scalar_run_group_state_step_advance_steps;
  int64_t refactor_btf_scalar_run_group_state_step_advance_current_count;
  int64_t refactor_btf_scalar_run_group_state_step_advance_rows;
  int64_t refactor_btf_scalar_run_group_state_step_advance_entries;
  int64_t
    refactor_btf_scalar_run_group_state_step_advance_ready_current_count;
  int64_t refactor_btf_scalar_run_group_state_step_advance_reject_count;
  int64_t refactor_btf_scalar_run_group_state_step_batch_step_count;
  int64_t refactor_btf_scalar_run_group_state_step_batch_current_count;
  int64_t refactor_btf_scalar_run_group_state_step_batch_entries;
  int64_t refactor_btf_scalar_run_group_state_step_batch_state_rows;
  int64_t refactor_btf_scalar_run_group_state_step_batch_unique_state_rows;
  int64_t refactor_btf_scalar_run_group_state_step_window_round_count;
  int64_t refactor_btf_scalar_run_group_state_step_window_current_count;
  int64_t refactor_btf_scalar_run_group_state_step_window_entries;
  int64_t refactor_btf_scalar_run_group_state_reject_count;
  int64_t refactor_btf_scalar_run_group_state_exec_current_count;
  int64_t refactor_btf_scalar_run_group_state_exec_skipped_deps;
  int64_t refactor_btf_scalar_run_group_state_exec_restored_rows;
  int64_t refactor_btf_scalar_run_group_state_exec_terminal_current_count;
  int64_t
    refactor_btf_scalar_run_group_state_exec_wake_terminal_current_count;
  int64_t
    refactor_btf_scalar_run_group_state_exec_wake_terminal_candidate_current_count;
  int64_t
    refactor_btf_scalar_run_group_state_exec_wake_terminal_nonpipeline_current_count;
  int64_t
    refactor_btf_scalar_run_group_state_exec_wake_terminal_dependency_miss_current_count;
  int64_t
    refactor_btf_scalar_run_group_state_exec_wake_terminal_claim_miss_current_count;
  int64_t refactor_btf_scalar_run_group_state_exec_remaining_deps;
  int64_t refactor_btf_scalar_run_group_state_exec_remaining_entries;
  int64_t refactor_btf_scalar_run_group_state_exec_reject_count;
  int64_t refactor_btf_scalar_run_group_state_exec_dispatch_bypass_current_count;
  int64_t
    refactor_btf_scalar_run_group_state_exec_dispatch_bypass_ready_current_count;
  int64_t refactor_btf_scalar_run_group_state_exec_not_ready_current_count;
  int64_t refactor_btf_scalar_run_group_state_exec_owned_ready_current_count;
  int64_t refactor_btf_scalar_run_group_state_exec_late_ready_current_count;
  int verified_rhs_reused;
  int64_t verified_rhs_reuse_count;
  int refactor_lean_choice;
  int64_t egraph_worker_spin_iters;
  /* Zero means no compact triangular-solve mirror is prepared. */
  int compact_solve_index_bytes;
  int compact_solve_fused_rhs;
  /* The retained numeric satisfies the normalized moderate-work,
     fragmented dominant-BTF policy boundary. */
  int moderate_fragmented_policy_eligible;
  /* The retained AMF numeric fits the normalized compact two-block policy
     and all recurring accuracy checks; individual packed kernels still
     validate their derived representations before dispatch. */
  int compact_amf_two_block_policy_eligible;
  /* The retained NodeNDP numeric passed the reciprocal-hub proposal,
     symbolic BTF/separator acceptance, and normalized factor work limits. */
  int dense_reciprocal_hub_policy_eligible;
  /* The retained transpose-AMD numeric has an exact symmetric degree-one
     fringe and passed the normalized BTF, fill, work, and pivot guards. */
  int symmetric_scalar_fringe_policy_eligible;
  /* The retained normal-AMD single block passed normalized symbolic work and
     measured pivot, fill, balance, and numeric-work guards. */
  int pivoted_high_work_single_block_policy_eligible;
  /* The retained unscaled dominant-BTF numeric has a mostly scalar fringe
     and passed the normalized PTS fill, work, pivot, and balance guards. */
  int low_work_many_fringe_btf_pts_policy_eligible;
  /* The retained normal-AMD numeric has a tiny, mostly scalar BTF fringe and
     passed the normalized high-work EGraph lifecycle guards. */
  int high_work_tiny_scalar_fringe_policy_eligible;
  /* The retained low-work AMD/PTS numeric has a hubbed input and an almost
     entirely scalar BTF fringe. */
  int low_work_hubbed_scalar_fringe_pts_policy_eligible;
  /* A nearly-missing-diagonal one-block symbolic selected value-aware
     matching before the first numeric factor. */
  int nearly_missing_diagonal_early_match_selected;
} kls_stats;

void kls_default_options(kls_options *options);
int kls_create(kls_solver **solver);
void kls_destroy(kls_solver *solver);

int kls_analyze_csc(kls_solver *solver,
                    kls_index_type index_type,
                    int64_t n,
                    const void *col_ptr,
                    const void *row_idx,
                    int index_base,
                    const kls_options *options);

int kls_analyze_csr(kls_solver *solver,
                    kls_index_type index_type,
                    int64_t n,
                    const void *row_ptr,
                    const void *col_idx,
                    int index_base,
                    const kls_options *options);

int kls_factor(kls_solver *solver, const double *values);
int kls_refactor(kls_solver *solver, const double *values);

int kls_solve(kls_solver *solver,
              int64_t nrhs,
              const double *b,
              int64_t ldb,
              double *x,
              int64_t ldx);

int kls_solve_transpose(kls_solver *solver,
                        int64_t nrhs,
                        const double *b,
                        int64_t ldb,
                        double *x,
                        int64_t ldx);

int kls_get_stats(const kls_solver *solver, kls_stats *stats);
const char *kls_status_string(int status);
const char *kls_ordering_name(kls_ordering ordering);
const char *kls_orientation_name(kls_orientation orientation);
const char *kls_backend_name(kls_backend backend);
const char *kls_factor_path_name(kls_factor_path path);
const char *kls_refactor_path_name(kls_refactor_path path);

#ifdef __cplusplus
}
#endif

#endif
