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

#ifdef __cplusplus
}
#endif

#endif
