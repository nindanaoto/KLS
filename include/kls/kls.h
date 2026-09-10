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
  KLS_ORDERING_AMF = 6,
  /* Deterministic approximate-minimum-fill variants used by AUTO's
     lifecycle tournament.  AMF retains the historical adaptive/default
     mode; AMMF amortizes local fill by supervariable size and AMF3 uses
     the tighter clique-deficiency bound. */
  KLS_ORDERING_AMMF = 7,
  KLS_ORDERING_AMF3 = 8
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
  KLS_REFACTOR_PATH_UNCHANGED = 7
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
  /* Expected calls after the initial factorization.  They are workload
     hints, not correctness controls: zero preserves the one-shot policy.
     Repeated-workload AUTO may use them to include value-remapping and
     triangular-solve traffic in representation selection. */
  int64_t expected_refactorizations;
  int64_t expected_solves;
  /* Set to zero when an outer integration owns timing and the overhead of
     timing each tiny singleton solve would be part of the measured work. */
  int record_tiny_solve_timing;
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
  /* Retired stream dependency/pivot/output diagnostics: always zero. */
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
  /* Retired supernode_consumer diagnostics are ABI-reserved, always zero. */
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
  int row_refactor_last_compact_supernode_partial_update;
  int64_t row_refactor_compact_supernode_partial_update_count;
  int64_t row_refactor_compact_supernode_partial_update_rows;
  int64_t row_refactor_compact_supernode_partial_update_entries;
  int64_t row_refactor_last_separator_flop_closure_groups;
  int64_t row_refactor_separator_flop_closure_group_count;
  /* Retired pipeline counters; retained for ABI compatibility, always zero. */
  /* Retired optional EGraph update/ready-queue counters: always zero. */
  double parallel_model_r1;
  double parallel_model_r2;
  int parallel_model_recommends_parallel;
  int row_refactor_auto_model_recommended;
  int row_refactor_auto_model_attempted;
  int row_refactor_auto_model_accepted;
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
  /* Retired L-pattern/adjacent-run/contiguous-suffix diagnostic fields are
     ABI-reserved, always zero. Sorted-L eligibility statistics remain active. */
  int internal_index_bytes;
  int refactor_map_index32_enabled;
  int64_t refactor_map_index32_entries;
  int refactor_l_index32_enabled;
  int64_t refactor_l_index32_entries;
  int refactor_u_index32_enabled;
  int64_t refactor_u_index32_entries;
  double row_refactor_auto_lower_bound_work;
  int row_refactor_auto_lower_bound_rejected;
  int row_refactor_auto_pattern_build_failed;
  int row_refactor_auto_value_copy_failed;
  /* Retired optional EGraph panel-cache counters: always zero. */
  /* Retired consumer_plan fields throughout this struct are ABI-reserved,
     always zero. */
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
  /* Retired u_supernode experiment fields throughout this struct are
     ABI-reserved, always zero. Ordinary supernode update statistics remain. */
  int64_t
    refactor_supernode_consumer_plan_first_dep_shape_batch_small_run_count;
  int64_t
    refactor_supernode_consumer_plan_first_dep_shape_batch_small_run_rows;
  int64_t
    refactor_supernode_consumer_plan_shape_bounded_advance_update_entries;
  int64_t
    refactor_supernode_consumer_plan_shape_bounded_advance_max_run_deps;
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
  /* Reserved ABI fields: retired EGraph prefactor experiment; always zero. */
  /* Retired Algorithm5 payoff statistics: retained for ABI compatibility;
     these fields through refactor_last_supernode_algorithm5_payoff_* are zero. */
  int64_t
    refactor_supernode_algorithm5_payoff_group_advance_duplicate_update_entries;
  int64_t
    refactor_supernode_algorithm5_payoff_group_max_run_suffix_update_entries;
  int64_t
    refactor_supernode_algorithm5_payoff_group_suffix_duplicate_update_entries;
  int64_t
    refactor_last_supernode_algorithm5_payoff_final_trigger_multi_batches;
  /* Retired surface-tracker counters (next eight fields): always zero.
     The following scalar_run_exec counters remain active. */
  int64_t refactor_last_btf_scalar_run_exec_runs;
  int64_t refactor_last_btf_scalar_run_exec_rows;
  int64_t refactor_last_btf_scalar_run_exec_entries;
  int64_t refactor_last_btf_scalar_run_exec_max_rows;
  int64_t refactor_btf_scalar_run_exec_count;
  int64_t refactor_btf_scalar_run_exec_rows;
  int64_t refactor_btf_scalar_run_exec_entries;
  int64_t refactor_btf_scalar_run_exec_max_rows;
  /* Retired btf_scalar_run_group fields are ABI-reserved, always zero.
     The ordinary BTF scalar-run executor remains supported. */
  int64_t
    refactor_btf_scalar_run_group_live_step_runtime_partial_current_count;
  int64_t
    refactor_btf_scalar_run_group_live_step_runtime_partial_active_unique_rows;
  int64_t
    refactor_last_btf_scalar_run_group_state_advance_batch_unique_entries;
  int64_t
    refactor_last_btf_scalar_run_group_state_advance_batch_duplicate_entries;
  int64_t
    refactor_last_btf_scalar_run_group_state_advance_batch_max_currents;
  int64_t
    refactor_last_btf_scalar_run_group_state_step_advance_ready_currents;
  int64_t
    refactor_last_btf_scalar_run_group_state_step_batch_unique_state_rows;
  int64_t
    refactor_last_btf_scalar_run_group_live_step_runtime_partial_active_rows;
  int64_t
    refactor_last_btf_scalar_run_group_live_step_runtime_partial_active_unique_rows;
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
  int64_t
    refactor_last_btf_scalar_run_group_state_exec_dispatch_bypass_ready_currents;
  int64_t
    refactor_btf_scalar_run_group_state_advance_batch_duplicate_entries;
  int64_t
    refactor_btf_scalar_run_group_state_step_advance_ready_current_count;
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
  int64_t
    refactor_btf_scalar_run_group_state_exec_dispatch_bypass_ready_current_count;
  int verified_rhs_reused;
  int64_t verified_rhs_reuse_count;
  int refactor_lean_choice;
  int64_t egraph_worker_spin_iters;
  /* Zero means no compact triangular-solve mirror is prepared. */
  int compact_solve_index_bytes;
  int compact_solve_fused_rhs;
  /* The retained numeric satisfies the normalized moderate-work,
     fragmented dominant-BTF policy boundary. */
  /* The retained AMF numeric fits the normalized compact two-block policy
     and all recurring accuracy checks; individual packed kernels still
     validate their derived representations before dispatch. */
  /* The retained NodeNDP numeric passed the reciprocal-hub proposal,
     symbolic BTF/separator acceptance, and normalized factor work limits. */
  /* The retained transpose-AMD numeric has an exact symmetric degree-one
     fringe and passed the normalized BTF, fill, work, and pivot guards. */
  /* The retained normal-AMD single block passed normalized symbolic work and
     measured pivot, fill, balance, and numeric-work guards. */
  /* The retained unscaled dominant-BTF numeric has a mostly scalar fringe
     and passed the normalized PTS fill, work, pivot, and balance guards. */
  /* The retained normal-AMD numeric has a tiny, mostly scalar BTF fringe and
     passed the normalized high-work EGraph lifecycle guards. */
  /* The retained low-work AMD/PTS numeric has a hubbed input and an almost
     entirely scalar BTF fringe. */
  /* A nearly-missing-diagonal one-block symbolic selected value-aware
     matching before the first numeric factor. */
  /* The retained factor has low work spread across many genuinely tiny BTF
     blocks and uses the bounded tiny-block lifecycle. */
  /* The selected symbolic has low estimated work spread across many tiny BTF
     blocks, independent of the later numeric pivot verdict. */
  /* The retained scaled fragmented factor fits the compact paired-row
     representation and its bounded-work lifecycle. */
  /* Singleton BTF blocks represented by the compact solve run cache, and
     the longest consecutive cached run.  Both are zero until a solve builds
     and adopts that optional representation. */
  int64_t compact_solve_singleton_run_blocks;
  int64_t compact_solve_singleton_run_max;
  /* The selected BTF geometry can amortize the compact singleton-run
     representation if a compact solve mirror is later built. */
  int compact_solve_singleton_run_eligible;
  /* The original AUTO input fits the compact bounded-degree,
     nearly-missing-diagonal matching representation. */
  /* The compact match covered the input and its accepted factor replaced the
     analyze-time proposal. */
  int compact_missing_diagonal_match_selected;
  /* The accepted match also produced the fragmented, low-work numeric needed
     by the compact row/direct-value lifecycle. */
  /* The original AUTO input is a bounded-degree symmetric graph with a
     material structural-diagonal defect suitable for value-aware matching. */
  /* A value-aware match for that proposal was adopted. */
  /* The adopted match produced the bounded, balanced single-block numeric
     required by the direct-value/PTS lifecycle. */
  /* The same structural family produced a smaller AMD/BTF symbolic for which
     the lightweight matched PTS lifecycle is preferable. */
  /* The retained NodeNDP symbolic passed the sparse symmetric scalar-fringe
     input proposal and normalized BTF/fill/separator acceptance. */
  /* That retained symbolic also produced an unscaled, no-pivot numeric inside
     the measured balanced fill/work regime. */
  /* The original normal AUTO input is a nearly diagonal sparse graph with a
     moderate spike, suitable for overlapped matching and NodeND. */
  /* Matching produced an accepted one-block predicted METIS numeric inside
     the normalized fill/work and separator bounds. */
  /* That factor also has the nearly all-private moderate-work separator for
     which the retained cluster schedule and relaxed-consume floors apply. */
  /* The retained scaled BTF factor has normalized dense-input,
     fragmentation, fill, work, and pivot bounds for direct row updates. */
  /* The original AUTO input is a bounded-degree, full-diagonal graph in a
     resource band where tuned METIS refinement can amortize row metadata. */
  /* Its retained one-block METIS symbolic has balanced intermediate fill and
     a mostly private separator inside the measured setup-cost bounds. */
  /* The retained fixed-pivot factor also preserves the normalized fill/work
     regime required by the recurring cooperative row lifecycle. */
  /* The original AUTO input is a giant, exactly structurally symmetric
     sparse graph with an almost-full diagonal, a scalar missing-diagonal
     fringe, and a bounded hub population. */
  /* Its raced METIS/BTF symbolic proved a giant core, mostly scalar fringe,
     bounded fill/work, and a complete mostly-private separator. */
  /* The retained unscaled fixed-pivot factor also passed normalized numeric
     fill/work, balance, and pivot-repair guards for recurring row updates. */
  /* The retained sparse, one-block METIS factor and its separator have the
     balanced fill/work and nearly all-private geometry used by the hybrid
     clustered-prefix/dependency-pipeline EGraph schedule. */
  /* A plain-frame factor selected a weaker-than-requested pivot threshold and
     retains the current matrix values needed to enforce a relative-L2 solve
     contract.  This is a numeric-state capability, not an input-shape tag. */
  int promoted_tolerance_l2_recovery_eligible;
  /* Number of public solves that entered that relative-L2 contract. */
  int64_t promoted_tolerance_l2_contract_run_count;
  /* Number of residual-stagnation episodes that rebuilt the numeric at a
     more conservative, request-relative pivot threshold. */
  int64_t promoted_tolerance_l2_recovery_count;
  /* A sparse, nearly diagonal-free, bounded in/out-degree input proposed the
     high-work retained-preconditioner lifecycle. */
  /* The directly measured transpose/AMD symbolic has a dominant SCC and
     enough fill/work to make retained-factor updates economically relevant. */
  /* The installed numeric preserved the symbolic regime and passed the
     measured factor-work and representation checks. */
  /* Number of changed-value refactors served by the retained numeric after
     the entrywise update bound accepted the new values. */
  /* The original AUTO input is a nearly full-diagonal asymmetric graph with
     bounded in/out degree and one nearly spanning raw SCC. */
  /* SCC fragmentation selected the coarse-fringe/eight-leaf class (1) or
     the thin-fringe/fourteen-leaf class (2). */
  /* Direct METIS/BTF analysis preserved that SCC class and proved the
     normalized dominant-core and separator economics. */
  /* The installed unscaled fixed-pivot numeric also passed normalized
     fill/work, balance, and pivot-repair guards. */
  /* The original AUTO input is a dense, nearly degree-balanced mega-hub
     graph with an almost-full diagonal and very small in/out-degree skew. */
  /* Its direct AMD/BTF symbolic proved full rank, dominant-core geometry,
     and normalized fill/work economics. */
  /* The installed max-row-scaled numeric also passed normalized pivot,
     fill, work, and representation guards. */
  /* The original AUTO input is a giant almost-full-diagonal graph with a
     dominant, degree-balanced in/out hub. */
  /* Its ordinary AUTO symbolic selected one unscaled METIS block and proved
     normalized fill plus a nearly all-private separator. */
  /* The installed full-rank pipelined numeric also proved the dense-tail,
     pivot, fill, work, and residual-refinement representation. */
  /* The retained full-rank min-fill factor has balanced high work and a
     1/512--1/256 mixed-component BTF fringe suitable for direct PTS solves. */
  /* The generic PTS builder proved a large elimination forest, bounded
     serial top, and balanced worker partition suitable for direct adoption. */
  int verified_large_pts_solve_policy_eligible;
  /* A medium full-diagonal input with one macroscopic column spike proposed
     the guarded AMMF/BTF ordering. */
  /* The selected AMMF symbolic also proved the dominant-core, fill, and work
     contract used by the recurring lifecycle. */
  /* A large full-diagonal sparse input with a bounded broad column proposed
     the one-block AMF/no-BTF ordering. */
  /* The selected one-block AMF symbolic passed its balanced fill/work budget. */
  /* A large almost-full-diagonal bounded-degree input proposed the guarded
     one-block AMF/no-BTF ordering. */
  /* The selected one-block AMF symbolic passed its bounded-work contract. */
  /* The current numeric admits the compact one-RHS all-singleton BTF solve. */
  int tiny_singleton_solve_eligible;
  /* Number of public solves dispatched through that compact kernel. */
  int64_t tiny_singleton_solve_count;
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

/* Refactor and immediately solve the changed numeric.  Unsupported fusion
   cases retain exactly the same semantics by falling back to kls_refactor
   followed by kls_solve. */
int kls_refactor_solve(kls_solver *solver,
                       const double *values,
                       int64_t nrhs,
                       const double *b,
                       int64_t ldb,
                       double *x,
                       int64_t ldx);

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
