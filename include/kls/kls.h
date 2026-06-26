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
  KLS_ORDERING_SCOTCH = 5
} kls_ordering;

typedef enum kls_orientation {
  KLS_ORIENTATION_AUTO = 0,
  KLS_ORIENTATION_NORMAL = 1,
  KLS_ORIENTATION_TRANSPOSE = 2
} kls_orientation;

typedef enum kls_factor_path {
  KLS_FACTOR_PATH_NONE = 0,
  KLS_FACTOR_PATH_KLU_FIRST = 1,
  KLS_FACTOR_PATH_KLS_FAST_REFACTOR = 2,
  KLS_FACTOR_PATH_KLU_FALLBACK = 3,
  KLS_FACTOR_PATH_PRESTATIC_KLU_FIRST = 4,
  KLS_FACTOR_PATH_KLS_FIRST = 5
} kls_factor_path;

#define KLS_SCALE_AUTO (-2)

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
  double pivot_tolerance;
  double memory_growth;
  int halt_if_singular;
  int fast_factor;
  int static_pivoting;
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
  int64_t kls_first_last_row_supernode_update_groups;
  int64_t kls_first_last_row_supernode_update_rows;
  int kls_first_last_separator_queue_pipeline_supernode_panel_update;
  int64_t kls_first_separator_queue_pipeline_supernode_panel_update_run_count;
  int64_t kls_first_last_separator_queue_pipeline_supernode_panel_update_groups;
  int64_t kls_first_last_separator_queue_pipeline_supernode_panel_update_rows;
  int kls_first_last_row_supernode_panel_update;
  int64_t kls_first_row_supernode_panel_update_run_count;
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
const char *kls_factor_path_name(kls_factor_path path);

#ifdef __cplusplus
}
#endif

#endif
