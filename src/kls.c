#define _POSIX_C_SOURCE 200809L

#include "kls/kls.h"

#include <stdio.h>

#define DLONG 1
#include "trilinos_klu_internal.h"
#undef DLONG
#include "trilinos_klu_decl.h"
#include "trilinos_camd.h"
#ifdef KLS_HAVE_METIS
#include "metis.h"
#endif
#ifdef KLS_HAVE_SCOTCH
#include "scotch.h"
#endif
#ifdef KLS_HAVE_SPRAL_SCALING
#include "spral_scaling.h"
#endif

#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define KLS_KLU_EMPTY ((UF_long)-1)
#define KLS_ROW_REFACTOR_DENSE_MIN_WORK 1024.0

typedef struct kls_refactor_pool kls_refactor_pool;
typedef struct kls_egraph_refactor_pool kls_egraph_refactor_pool;

typedef enum kls_input_format {
  KLS_INPUT_NONE = 0,
  KLS_INPUT_CSC = 1,
  KLS_INPUT_CSR = 2
} kls_input_format;

typedef enum kls_row_refactor_group_kind {
  KLS_ROW_REFACTOR_GROUP_SINGLE = 0,
  KLS_ROW_REFACTOR_GROUP_GENERIC = 1,
  KLS_ROW_REFACTOR_GROUP_DENSE = 2
} kls_row_refactor_group_kind;

struct kls_solver {
  UF_long n;
  UF_long nnz;
  UF_long *col_ptr;
  UF_long *row_idx;
  UF_long *input_to_csc;
  UF_long *row_perm;
  double *row_scale;
  double *col_scale;
  double *values;
  double *solve_perm_workspace;
  UF_long solve_perm_workspace_n;
  kls_egraph_refactor_pool *egraph_pool;
  UF_long *refactor_col_ptr;
  UF_long *refactor_row_idx;
  UF_long *refactor_input_pos;
  UF_long *refactor_block_start;
  UF_long *refactor_col_block;
  UF_long **refactor_l_indices;
  double **refactor_l_values;
  UF_long **refactor_u_indices;
  double **refactor_u_values;
  UF_long refactor_lu_pointer_count;
  UF_long *row_refactor_l_ptr;
  UF_long *row_refactor_l_cols;
  double **row_refactor_l_values;
  double *row_refactor_l_row_values;
  UF_long *row_refactor_u_ptr;
  UF_long *row_refactor_u_cols;
  double **row_refactor_u_values;
  double *row_refactor_u_row_values;
  UF_long *row_refactor_input_ptr;
  UF_long *row_refactor_input_cols;
  UF_long *row_refactor_input_pos;
  UF_long *row_refactor_successor_ptr;
  UF_long *row_refactor_successor_rows;
  UF_long *row_refactor_tail_rows;
  unsigned int *row_refactor_tail_marks;
  UF_long *row_refactor_level_ptr;
  UF_long *row_refactor_level_rows;
  UF_long *row_refactor_group_ptr;
  UF_long *row_refactor_group_dep_ptr;
  UF_long *row_refactor_group_dep_rows;
  UF_long *row_refactor_group_level_ptr;
  UF_long *row_refactor_level_groups;
  UF_long *row_refactor_group_level_thread_ptr;
  int row_refactor_group_level_thread_count;
  UF_long *row_refactor_l_internal_ptr;
  UF_long *row_refactor_etree_parent;
  UF_long *row_refactor_group_trailing_len;
  unsigned char *row_refactor_group_dense;
  unsigned char *row_refactor_group_kind;
  UF_long row_refactor_pattern_n;
  UF_long row_refactor_group_count;
  UF_long row_refactor_level_count;
  UF_long row_refactor_level_max_width;
  UF_long row_refactor_cluster_level_count;
  UF_long row_refactor_pipeline_group_count;
  UF_long row_refactor_pipeline_row_count;
  double row_refactor_pipeline_work;
  int row_refactor_last_run;
  int row_refactor_last_checked;
  int row_refactor_last_parallel;
  UF_long row_refactor_run_count;
  UF_long row_refactor_checked_run_count;
  UF_long row_refactor_parallel_run_count;
  UF_long row_refactor_segment_count;
  UF_long row_refactor_segment_rows;
  UF_long row_refactor_segment_max_width;
  double row_refactor_segment_dense_entries;
  double row_refactor_segment_trailing_entries;
  UF_long row_refactor_dense_segment_count;
  UF_long row_refactor_dense_segment_rows;
  UF_long row_refactor_dense_segment_max_width;
  double row_refactor_dense_segment_dense_entries;
  double row_refactor_dense_segment_trailing_entries;
  unsigned int row_refactor_tail_mark;
  UF_long row_refactor_tail_count;
  UF_long *refactor_level_ptr;
  UF_long *refactor_level_cols;
  UF_long *refactor_level_thread_ptr;
  int refactor_level_thread_count;
  UF_long refactor_level_count;
  UF_long refactor_level_max_width;
  UF_long refactor_dependency_edges;
  UF_long refactor_dependency_root_columns;
  UF_long refactor_dependency_leaf_columns;
  UF_long refactor_dependency_max_fanout;
  double refactor_dependency_max_column_work;
  double refactor_dependency_pipeline_max_column_work;
  UF_long refactor_supernode_candidate_count;
  UF_long refactor_supernode_candidate_rows;
  UF_long refactor_supernode_candidate_max_width;
  double refactor_supernode_candidate_dense_entries;
  double refactor_supernode_candidate_trailing_entries;
  UF_long refactor_cluster_level_count;
  UF_long refactor_pipeline_column_count;
  double refactor_dependency_work;
  double refactor_pipeline_work;
  double **egraph_worker_scratch;
  UF_long egraph_worker_scratch_size;
  int egraph_worker_scratch_count;
  int egraph_worker_scratch_dirty;
  atomic_uint *egraph_pipeline_done;
  UF_long egraph_pipeline_done_size;
  unsigned int egraph_pipeline_generation;
  kls_input_format input_format;
  kls_orientation orientation;
  kls_options options;
  kls_stats stats;
  trilinos_klu_l_common common;
  trilinos_klu_l_symbolic *symbolic;
  trilinos_klu_l_numeric *numeric;
  kls_refactor_pool *refactor_pool;
  int auto_metis_checked;
  int auto_pivot_checked;
  int auto_scale_checked;
  int exact_matching_selected;
  int spral_matching_selected;
  int fast_block_restarts;
  int fast_tail_restarts;
  int fast_reject_refresh_state;
  UF_long *fast_reject_tail_cols;
  unsigned int *fast_reject_tail_marks;
  UF_long fast_reject_tail_capacity;
  unsigned int fast_reject_tail_mark;
  UF_long fast_reject_tail_count;
};

typedef struct kls_pattern_candidate {
  UF_long n;
  UF_long nnz;
  UF_long *col_ptr;
  UF_long *row_idx;
  UF_long *input_to_csc;
  double *values;
  kls_orientation orientation;
  trilinos_klu_l_common common;
  trilinos_klu_l_symbolic *symbolic;
  kls_ordering selected_ordering;
  double score;
} kls_pattern_candidate;

typedef struct kls_parallel_refactor_shared {
  UF_long n;
  UF_long nnz;
  const UF_long *col_ptr;
  const UF_long *row_idx;
  const UF_long *map_col_ptr;
  const UF_long *map_row_idx;
  const UF_long *map_input_pos;
  const UF_long *map_block_start;
  const double *values;
  const trilinos_klu_l_symbolic *symbolic;
  trilinos_klu_l_numeric *numeric;
  const double *rs;
  int scale;
  int halt_if_singular;
  int check_pivots;
  double pivot_tolerance;
  UF_long next_block;
  UF_long block_chunk;
  unsigned char *block_done;
  int stop;
  pthread_mutex_t lock;
} kls_parallel_refactor_shared;

typedef struct kls_parallel_refactor_worker {
  kls_parallel_refactor_shared *shared;
  kls_refactor_pool *pool;
  double *x;
  int invalid;
  int pivot_rejected;
  int singular;
  UF_long rejected_pivot;
  UF_long rejected_pivot_col;
  UF_long rejected_row;
  double rejected_multiplier_abs;
  double rejected_pivot_abs;
  double rejected_candidate_abs;
  UF_long numerical_rank;
  UF_long singular_col;
} kls_parallel_refactor_worker;

typedef enum kls_egraph_refactor_kernel {
  KLS_EGRAPH_REFACTOR_KERNEL_GENERIC = 0,
  KLS_EGRAPH_REFACTOR_KERNEL_SINGLE_UNSCALED = 1,
  KLS_EGRAPH_REFACTOR_KERNEL_BTF_UNSCALED = 2
} kls_egraph_refactor_kernel;

typedef struct kls_egraph_refactor_shared {
  kls_solver *solver;
  const double *values;
  const double *rs;
  int check_pivots;
  int scale;
  kls_egraph_refactor_kernel kernel;
  int thread_count;
  atomic_int stop;
  int invalid;
  int pivot_rejected;
  int singular;
  UF_long rejected_pivot;
  UF_long rejected_pivot_col;
  UF_long rejected_row;
  double rejected_multiplier_abs;
  double rejected_pivot_abs;
  double rejected_candidate_abs;
  UF_long numerical_rank;
  UF_long singular_col;
  pthread_mutex_t lock;
  atomic_uint *pipeline_done;
  unsigned int pipeline_generation;
  atomic_ulong next_pipeline_pos;
  UF_long pipeline_pos_end;
  UF_long cluster_level_count;
  int pipeline_natural_order;
  int row_refactor_mode;
  int row_refactor_defer_value_scatter;
  pthread_barrier_t barrier;
} kls_egraph_refactor_shared;

typedef struct kls_egraph_refactor_worker {
  kls_egraph_refactor_shared *shared;
  kls_egraph_refactor_pool *pool;
  int tid;
  double *x;
  double *segment_panel;
  UF_long segment_panel_size;
} kls_egraph_refactor_worker;

struct kls_egraph_refactor_pool {
  kls_egraph_refactor_shared shared;
  pthread_cond_t work_cond;
  pthread_cond_t done_cond;
  pthread_t *threads;
  kls_egraph_refactor_worker *workers;
  unsigned long generation;
  int thread_count;
  int created_count;
  int active_workers;
  int shutdown;
  int conds_initialized;
  int lock_initialized;
  int barrier_initialized;
};

struct kls_refactor_pool {
  kls_parallel_refactor_shared shared;
  pthread_cond_t work_cond;
  pthread_cond_t done_cond;
  pthread_t *threads;
  kls_parallel_refactor_worker *workers;
  unsigned char *block_done;
  UF_long block_done_capacity;
  UF_long maxblock;
  unsigned long generation;
  int thread_count;
  int created_count;
  int active_workers;
  int shutdown;
  int scratch_dirty;
  int conds_initialized;
  int lock_initialized;
};

typedef struct kls_match_entry {
  double weight;
  UF_long row;
  UF_long col;
} kls_match_entry;

static int kls_build_refactor_schedule(kls_solver *solver);
static UF_long kls_block_for_pivot(const kls_solver *solver, UF_long pivot);
static void kls_fill_fast_reject_tail_stats(kls_solver *solver,
                                            UF_long rejected_pivot);
static void kls_fill_fast_reject_observed_tail_candidate(kls_solver *solver);
static void destroy_egraph_refactor_pool(kls_solver *solver);
static void kls_egraph_refactor_record_invalid(
  kls_egraph_refactor_shared *shared);
static void kls_egraph_refactor_record_reject(
  kls_egraph_refactor_shared *shared,
  UF_long rejected_pivot,
  UF_long rejected_pivot_col,
  UF_long rejected_row,
  double rejected_multiplier_abs,
  double rejected_pivot_abs,
  double rejected_candidate_abs);
static void kls_egraph_refactor_record_singular(
  kls_egraph_refactor_shared *shared,
  UF_long numerical_rank,
  UF_long singular_col);
static int kls_egraph_refactor_should_stop(
  kls_egraph_refactor_shared *shared);
static void kls_egraph_refactor_mark_done(
  kls_egraph_refactor_shared *shared,
  UF_long col);
static int kls_egraph_refactor_wait_done(
  kls_egraph_refactor_shared *shared,
  UF_long col);
static int kls_egraph_refreshed_prefix(
  const kls_egraph_refactor_shared *shared,
  UF_long rejected_pivot);
static kls_egraph_refactor_pool *ensure_egraph_refactor_pool(
  kls_solver *solver,
  int thread_count);

#define KLS_MATCH_PATH_MAX_DEPTH 4u

typedef struct kls_weighted_path_search {
  const struct kls_row_match_graph *graph;
  const UF_long *row_perm;
  const UF_long *col_match;
  const double *current_log_weight;
  UF_long max_candidates_per_row;
  UF_long max_depth;
  UF_long path_rows[KLS_MATCH_PATH_MAX_DEPTH];
  UF_long path_cols[KLS_MATCH_PATH_MAX_DEPTH];
  double path_weights[KLS_MATCH_PATH_MAX_DEPTH];
  UF_long best_rows[KLS_MATCH_PATH_MAX_DEPTH];
  UF_long best_cols[KLS_MATCH_PATH_MAX_DEPTH];
  double best_weights[KLS_MATCH_PATH_MAX_DEPTH];
  UF_long best_len;
  double best_gain;
} kls_weighted_path_search;

typedef struct kls_row_match_graph {
  UF_long *row_ptr;
  UF_long *col_idx;
  double *log_weight;
} kls_row_match_graph;

typedef struct kls_heap_item {
  size_t node;
  double key;
} kls_heap_item;

typedef struct kls_min_heap {
  kls_heap_item *items;
  size_t size;
  size_t capacity;
} kls_min_heap;

typedef struct kls_row_permuted_entry {
  UF_long col;
  UF_long row;
  UF_long base_position;
  double value;
} kls_row_permuted_entry;

static int choose_symbolic_for_pattern(UF_long n,
                                       UF_long *col_ptr,
                                       UF_long *row_idx,
                                       const kls_options *options,
                                       trilinos_klu_l_symbolic **symbolic_out,
                                       trilinos_klu_l_common *common_out,
                                       kls_ordering *selected_ordering_out,
                                       double *score_out);
static void free_candidate(kls_pattern_candidate *candidate);
static int transpose_candidate(const kls_pattern_candidate *source,
                               kls_orientation orientation,
                               kls_pattern_candidate *candidate);

static double kls_now_seconds(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static size_t kls_klu_units_for_indices(UF_long length) {
  const size_t bytes = (size_t)length * sizeof(UF_long);
  return (bytes + sizeof(double) - 1u) / sizeof(double);
}

static void kls_klu_get_pointer(double *lu,
                                const UF_long *offsets,
                                const UF_long *lengths,
                                UF_long k,
                                UF_long **indices_out,
                                double **values_out,
                                UF_long *length_out) {
  const UF_long length = lengths[k];
  double *base = lu + offsets[k];
  *indices_out = (UF_long *)base;
  *values_out = (double *)(base + kls_klu_units_for_indices(length));
  *length_out = length;
}

static inline void kls_scatter_subtract(double *restrict x,
                                        const UF_long *restrict rows,
                                        const double *restrict values,
                                        UF_long length,
                                        double scale) {
  UF_long p = 0;
  for (; p + 3u < length; p += 4u) {
    x[rows[p]] -= values[p] * scale;
    x[rows[p + 1u]] -= values[p + 1u] * scale;
    x[rows[p + 2u]] -= values[p + 2u] * scale;
    x[rows[p + 3u]] -= values[p + 3u] * scale;
  }
  for (; p < length; ++p) {
    x[rows[p]] -= values[p] * scale;
  }
}

static inline uint64_t kls_mix_u64(uint64_t x) {
  x ^= x >> 33u;
  x *= UINT64_C(0xff51afd7ed558ccd);
  x ^= x >> 33u;
  x *= UINT64_C(0xc4ceb9fe1a85ec53);
  x ^= x >> 33u;
  return x;
}

static inline uint64_t kls_supernode_hash1(UF_long col) {
  return kls_mix_u64((uint64_t)col + UINT64_C(0x9e3779b97f4a7c15));
}

static inline uint64_t kls_supernode_hash2(UF_long col) {
  return kls_mix_u64((uint64_t)col ^ UINT64_C(0xd1b54a32d192ed03));
}

static void kls_record_supernode_candidate(
  kls_solver *solver,
  UF_long start,
  UF_long width,
  const UF_long *successor_counts) {
  if (solver == NULL || successor_counts == NULL || width <= 1u) {
    return;
  }
  solver->refactor_supernode_candidate_count++;
  solver->refactor_supernode_candidate_rows += width;
  if (width > solver->refactor_supernode_candidate_max_width) {
    solver->refactor_supernode_candidate_max_width = width;
  }
  solver->refactor_supernode_candidate_dense_entries +=
    (double)width * (double)(width - 1u) * 0.5;
  solver->refactor_supernode_candidate_trailing_entries +=
    (double)width * (double)successor_counts[start + width - 1u];
}

static void free_refactor_map(kls_solver *solver) {
  if (solver == NULL) {
    return;
  }
  free(solver->refactor_col_ptr);
  free(solver->refactor_row_idx);
  free(solver->refactor_input_pos);
  free(solver->refactor_block_start);
  free(solver->refactor_col_block);
  solver->refactor_col_ptr = NULL;
  solver->refactor_row_idx = NULL;
  solver->refactor_input_pos = NULL;
  solver->refactor_block_start = NULL;
  solver->refactor_col_block = NULL;
}

static void free_refactor_lu_pointer_cache(kls_solver *solver) {
  if (solver == NULL) {
    return;
  }
  free(solver->refactor_l_indices);
  free(solver->refactor_l_values);
  free(solver->refactor_u_indices);
  free(solver->refactor_u_values);
  solver->refactor_l_indices = NULL;
  solver->refactor_l_values = NULL;
  solver->refactor_u_indices = NULL;
  solver->refactor_u_values = NULL;
  solver->refactor_lu_pointer_count = 0;
}

static void free_fast_reject_tail_plan(kls_solver *solver) {
  if (solver == NULL) {
    return;
  }
  free(solver->fast_reject_tail_cols);
  free(solver->fast_reject_tail_marks);
  solver->fast_reject_tail_cols = NULL;
  solver->fast_reject_tail_marks = NULL;
  solver->fast_reject_tail_capacity = 0;
  solver->fast_reject_tail_mark = 0u;
  solver->fast_reject_tail_count = 0;
}

static void free_row_refactor_pattern(kls_solver *solver) {
  if (solver == NULL) {
    return;
  }
  free(solver->row_refactor_l_ptr);
  free(solver->row_refactor_l_cols);
  free(solver->row_refactor_l_values);
  free(solver->row_refactor_l_row_values);
  free(solver->row_refactor_u_ptr);
  free(solver->row_refactor_u_cols);
  free(solver->row_refactor_u_values);
  free(solver->row_refactor_u_row_values);
  free(solver->row_refactor_input_ptr);
  free(solver->row_refactor_input_cols);
  free(solver->row_refactor_input_pos);
  free(solver->row_refactor_successor_ptr);
  free(solver->row_refactor_successor_rows);
  free(solver->row_refactor_tail_rows);
  free(solver->row_refactor_tail_marks);
  free(solver->row_refactor_level_ptr);
  free(solver->row_refactor_level_rows);
  free(solver->row_refactor_group_ptr);
  free(solver->row_refactor_group_dep_ptr);
  free(solver->row_refactor_group_dep_rows);
  free(solver->row_refactor_group_level_ptr);
  free(solver->row_refactor_level_groups);
  free(solver->row_refactor_group_level_thread_ptr);
  free(solver->row_refactor_l_internal_ptr);
  free(solver->row_refactor_etree_parent);
  free(solver->row_refactor_group_trailing_len);
  free(solver->row_refactor_group_dense);
  free(solver->row_refactor_group_kind);
  solver->row_refactor_l_ptr = NULL;
  solver->row_refactor_l_cols = NULL;
  solver->row_refactor_l_values = NULL;
  solver->row_refactor_l_row_values = NULL;
  solver->row_refactor_u_ptr = NULL;
  solver->row_refactor_u_cols = NULL;
  solver->row_refactor_u_values = NULL;
  solver->row_refactor_u_row_values = NULL;
  solver->row_refactor_input_ptr = NULL;
  solver->row_refactor_input_cols = NULL;
  solver->row_refactor_input_pos = NULL;
  solver->row_refactor_successor_ptr = NULL;
  solver->row_refactor_successor_rows = NULL;
  solver->row_refactor_tail_rows = NULL;
  solver->row_refactor_tail_marks = NULL;
  solver->row_refactor_level_ptr = NULL;
  solver->row_refactor_level_rows = NULL;
  solver->row_refactor_group_ptr = NULL;
  solver->row_refactor_group_dep_ptr = NULL;
  solver->row_refactor_group_dep_rows = NULL;
  solver->row_refactor_group_level_ptr = NULL;
  solver->row_refactor_level_groups = NULL;
  solver->row_refactor_group_level_thread_ptr = NULL;
  solver->row_refactor_group_level_thread_count = 0;
  solver->row_refactor_l_internal_ptr = NULL;
  solver->row_refactor_etree_parent = NULL;
  solver->row_refactor_group_trailing_len = NULL;
  solver->row_refactor_group_dense = NULL;
  solver->row_refactor_group_kind = NULL;
  solver->row_refactor_pattern_n = 0;
  solver->row_refactor_group_count = 0;
  solver->row_refactor_level_count = 0;
  solver->row_refactor_level_max_width = 0;
  solver->row_refactor_cluster_level_count = 0;
  solver->row_refactor_pipeline_group_count = 0;
  solver->row_refactor_pipeline_row_count = 0;
  solver->row_refactor_pipeline_work = 0.0;
  solver->row_refactor_last_run = 0;
  solver->row_refactor_last_checked = 0;
  solver->row_refactor_last_parallel = 0;
  solver->row_refactor_run_count = 0;
  solver->row_refactor_checked_run_count = 0;
  solver->row_refactor_parallel_run_count = 0;
  solver->row_refactor_segment_count = 0;
  solver->row_refactor_segment_rows = 0;
  solver->row_refactor_segment_max_width = 0;
  solver->row_refactor_segment_dense_entries = 0.0;
  solver->row_refactor_segment_trailing_entries = 0.0;
  solver->row_refactor_dense_segment_count = 0;
  solver->row_refactor_dense_segment_rows = 0;
  solver->row_refactor_dense_segment_max_width = 0;
  solver->row_refactor_dense_segment_dense_entries = 0.0;
  solver->row_refactor_dense_segment_trailing_entries = 0.0;
  solver->row_refactor_tail_mark = 0u;
  solver->row_refactor_tail_count = 0;
}

static void free_refactor_schedule(kls_solver *solver) {
  if (solver == NULL) {
    return;
  }
  free(solver->refactor_level_ptr);
  free(solver->refactor_level_cols);
  free(solver->refactor_level_thread_ptr);
  solver->refactor_level_ptr = NULL;
  solver->refactor_level_cols = NULL;
  solver->refactor_level_thread_ptr = NULL;
  solver->refactor_level_thread_count = 0;
  solver->refactor_level_count = 0;
  solver->refactor_level_max_width = 0;
  solver->refactor_dependency_edges = 0;
  solver->refactor_dependency_root_columns = 0;
  solver->refactor_dependency_leaf_columns = 0;
  solver->refactor_dependency_max_fanout = 0;
  solver->refactor_dependency_max_column_work = 0.0;
  solver->refactor_dependency_pipeline_max_column_work = 0.0;
  solver->refactor_supernode_candidate_count = 0;
  solver->refactor_supernode_candidate_rows = 0;
  solver->refactor_supernode_candidate_max_width = 0;
  solver->refactor_supernode_candidate_dense_entries = 0.0;
  solver->refactor_supernode_candidate_trailing_entries = 0.0;
  solver->refactor_cluster_level_count = 0;
  solver->refactor_pipeline_column_count = 0;
  solver->refactor_dependency_work = 0.0;
  solver->refactor_pipeline_work = 0.0;
}

static void free_egraph_worker_scratch(kls_solver *solver) {
  if (solver == NULL) {
    return;
  }
  if (solver->egraph_worker_scratch != NULL) {
    for (int i = 0; i < solver->egraph_worker_scratch_count; ++i) {
      free(solver->egraph_worker_scratch[i]);
    }
  }
  free(solver->egraph_worker_scratch);
  solver->egraph_worker_scratch = NULL;
  solver->egraph_worker_scratch_size = 0;
  solver->egraph_worker_scratch_count = 0;
  solver->egraph_worker_scratch_dirty = 0;
}

static double **ensure_egraph_worker_scratch(kls_solver *solver,
                                             int thread_count,
                                             UF_long scratch_size) {
  if (solver == NULL || thread_count <= 0 || scratch_size == 0u) {
    return NULL;
  }
  if (solver->egraph_worker_scratch != NULL &&
      solver->egraph_worker_scratch_count == thread_count &&
      solver->egraph_worker_scratch_size >= scratch_size) {
    if (solver->egraph_worker_scratch_dirty) {
      for (int i = 0; i < thread_count; ++i) {
        memset(solver->egraph_worker_scratch[i], 0,
               (size_t)solver->egraph_worker_scratch_size *
                 sizeof(*solver->egraph_worker_scratch[i]));
      }
      solver->egraph_worker_scratch_dirty = 0;
    }
    return solver->egraph_worker_scratch;
  }

  free_egraph_worker_scratch(solver);
  double **scratch =
    (double **)calloc((size_t)thread_count, sizeof(*scratch));
  if (scratch == NULL) {
    return NULL;
  }
  int allocated = 0;
  for (; allocated < thread_count; ++allocated) {
    scratch[allocated] =
      (double *)calloc((size_t)scratch_size, sizeof(*scratch[allocated]));
    if (scratch[allocated] == NULL) {
      break;
    }
  }
  if (allocated != thread_count) {
    for (int i = 0; i < allocated; ++i) {
      free(scratch[i]);
    }
    free(scratch);
    return NULL;
  }
  solver->egraph_worker_scratch = scratch;
  solver->egraph_worker_scratch_count = thread_count;
  solver->egraph_worker_scratch_size = scratch_size;
  solver->egraph_worker_scratch_dirty = 0;
  return scratch;
}

static void free_egraph_pipeline_done(kls_solver *solver) {
  if (solver == NULL) {
    return;
  }
  free(solver->egraph_pipeline_done);
  solver->egraph_pipeline_done = NULL;
  solver->egraph_pipeline_done_size = 0;
  solver->egraph_pipeline_generation = 0;
}

static atomic_uint *ensure_egraph_pipeline_done(
  kls_solver *solver,
  unsigned int *generation_out) {
  if (generation_out != NULL) {
    *generation_out = 0;
  }
  if (solver == NULL || solver->n == 0u) {
    return NULL;
  }
  if (solver->egraph_pipeline_done == NULL ||
      solver->egraph_pipeline_done_size != solver->n) {
    free_egraph_pipeline_done(solver);
    solver->egraph_pipeline_done =
      (atomic_uint *)malloc((size_t)solver->n *
                            sizeof(*solver->egraph_pipeline_done));
    if (solver->egraph_pipeline_done == NULL) {
      return NULL;
    }
    solver->egraph_pipeline_done_size = solver->n;
    for (UF_long k = 0; k < solver->n; ++k) {
      atomic_init(&solver->egraph_pipeline_done[k], 0u);
    }
  } else if (solver->egraph_pipeline_generation == UINT_MAX) {
    for (UF_long k = 0; k < solver->egraph_pipeline_done_size; ++k) {
      atomic_store_explicit(&solver->egraph_pipeline_done[k], 0u,
                            memory_order_relaxed);
    }
    solver->egraph_pipeline_generation = 0;
  }

  solver->egraph_pipeline_generation++;
  if (generation_out != NULL) {
    *generation_out = solver->egraph_pipeline_generation;
  }
  return solver->egraph_pipeline_done;
}

static void kls_worker_record_singular(kls_parallel_refactor_worker *worker,
                                       UF_long numerical_rank,
                                       UF_long singular_col) {
  if (!worker->singular || numerical_rank < worker->numerical_rank) {
    worker->singular = 1;
    worker->numerical_rank = numerical_rank;
    worker->singular_col = singular_col;
  }
}

static void kls_clear_fast_reject_stats(kls_solver *solver) {
  if (solver == NULL) {
    return;
  }
  solver->stats.fast_rejected_pivot = -1;
  solver->stats.fast_rejected_pivot_col = -1;
  solver->stats.fast_rejected_row = -1;
  solver->stats.fast_rejected_multiplier_abs = -1.0;
  solver->stats.fast_rejected_pivot_abs = -1.0;
  solver->stats.fast_rejected_candidate_abs = -1.0;
  solver->stats.fast_rejected_tail_candidate_row = -1;
  solver->stats.fast_rejected_tail_candidate_abs = -1.0;
  solver->stats.fast_rejected_tail_candidate_count = 0;
  solver->stats.fast_rejected_tail_candidate_position = -1;
  solver->stats.fast_rejected_tail_repair_ready = 0;
  solver->stats.fast_repaired_pivot_row = -1;
  solver->stats.fast_repaired_pivot_matches_tail_candidate = 0;
  solver->stats.fast_repaired_first_changed_pivot = -1;
  solver->stats.fast_repaired_prefix_changed_pivots = 0;
  solver->stats.fast_repaired_suffix_changed_pivots = 0;
  solver->stats.fast_repaired_tail_restart_ready = 0;
  solver->stats.fast_repaired_block_work = 0.0;
  solver->stats.fast_repaired_tail_restart_columns = 0;
  solver->stats.fast_repaired_tail_restart_work = 0.0;
  solver->stats.fast_repaired_tail_restart_saved_work = 0.0;
  solver->stats.fast_block_restarts = 0;
  solver->stats.fast_tail_restarts = 0;
  solver->stats.fast_rejected_block_start = -1;
  solver->stats.fast_rejected_block_size = 0;
  solver->stats.fast_rejected_suffix_columns = 0;
  solver->stats.fast_rejected_descendant_columns = 0;
  solver->stats.fast_rejected_descendant_work = 0.0;
  solver->stats.fast_rejected_row_tail_columns = 0;
  solver->stats.fast_rejected_row_tail_work = 0.0;
  solver->stats.fast_rejected_etree_columns = 0;
  solver->stats.fast_rejected_etree_work = 0.0;
  solver->stats.fast_rejected_pivoting_tail_columns = 0;
  solver->stats.fast_rejected_pivoting_tail_work = 0.0;
  solver->stats.fast_rejected_pivoting_tail_first = -1;
  solver->stats.fast_rejected_pivoting_tail_last = -1;
  solver->stats.fast_rejected_pivoting_tail_contains_reject = 0;
  solver->stats.fast_rejected_pivoting_tail_topological = 0;
  solver->stats.fast_rejected_refresh_state =
    KLS_FAST_REJECT_REFRESH_UNKNOWN;
  solver->fast_block_restarts = 0;
  solver->fast_tail_restarts = 0;
  solver->fast_reject_refresh_state = KLS_FAST_REJECT_REFRESH_UNKNOWN;
  solver->fast_reject_tail_count = 0;
}

static void kls_clear_row_refactor_last_stats(kls_solver *solver) {
  if (solver == NULL) {
    return;
  }
  solver->row_refactor_last_run = 0;
  solver->row_refactor_last_checked = 0;
  solver->row_refactor_last_parallel = 0;
  solver->stats.row_refactor_last_run = 0;
  solver->stats.row_refactor_last_checked = 0;
  solver->stats.row_refactor_last_parallel = 0;
}

static void kls_record_row_refactor_run(kls_solver *solver,
                                        int check_pivots,
                                        int parallel) {
  if (solver == NULL) {
    return;
  }
  solver->row_refactor_last_run = 1;
  solver->row_refactor_last_checked = check_pivots ? 1 : 0;
  solver->row_refactor_last_parallel = parallel ? 1 : 0;
  solver->row_refactor_run_count++;
  if (check_pivots) {
    solver->row_refactor_checked_run_count++;
  }
  if (parallel) {
    solver->row_refactor_parallel_run_count++;
  }
}

static void kls_record_fast_reject_detail(kls_solver *solver,
                                          UF_long rejected_pivot,
                                          UF_long rejected_pivot_col,
                                          UF_long rejected_row,
                                          double rejected_multiplier_abs,
                                          double rejected_pivot_abs,
                                          double rejected_candidate_abs) {
  if (solver == NULL || rejected_pivot == KLS_KLU_EMPTY) {
    return;
  }
  solver->stats.fast_rejected_pivot = (int64_t)rejected_pivot;
  solver->stats.fast_rejected_pivot_col =
    rejected_pivot_col == KLS_KLU_EMPTY ? -1 : (int64_t)rejected_pivot_col;
  solver->stats.fast_rejected_row =
    rejected_row == KLS_KLU_EMPTY ? -1 : (int64_t)rejected_row;
  solver->stats.fast_rejected_multiplier_abs =
    rejected_multiplier_abs >= 0.0 ? rejected_multiplier_abs : -1.0;
  solver->stats.fast_rejected_pivot_abs =
    rejected_pivot_abs >= 0.0 ? rejected_pivot_abs : -1.0;
  solver->stats.fast_rejected_candidate_abs =
    rejected_candidate_abs >= 0.0 ? rejected_candidate_abs : -1.0;
  solver->stats.fast_rejected_tail_candidate_row = -1;
  solver->stats.fast_rejected_tail_candidate_abs = -1.0;
  solver->stats.fast_rejected_tail_candidate_count = 0;
  solver->stats.fast_rejected_tail_candidate_position = -1;
  solver->stats.fast_rejected_tail_repair_ready = 0;
  solver->stats.fast_repaired_pivot_row = -1;
  solver->stats.fast_repaired_pivot_matches_tail_candidate = 0;
  solver->stats.fast_repaired_first_changed_pivot = -1;
  solver->stats.fast_repaired_prefix_changed_pivots = 0;
  solver->stats.fast_repaired_suffix_changed_pivots = 0;
  solver->stats.fast_repaired_tail_restart_ready = 0;
  solver->stats.fast_repaired_block_work = 0.0;
  solver->stats.fast_repaired_tail_restart_columns = 0;
  solver->stats.fast_repaired_tail_restart_work = 0.0;
  solver->stats.fast_repaired_tail_restart_saved_work = 0.0;
  solver->stats.fast_rejected_refresh_state =
    solver->fast_reject_refresh_state;
  kls_fill_fast_reject_tail_stats(solver, rejected_pivot);
  kls_fill_fast_reject_observed_tail_candidate(solver);
}

static void kls_record_fast_reject(kls_solver *solver,
                                   UF_long rejected_pivot,
                                   UF_long rejected_pivot_col) {
  kls_record_fast_reject_detail(solver, rejected_pivot, rejected_pivot_col,
                                KLS_KLU_EMPTY, -1.0, -1.0, -1.0);
}

static int kls_parallel_refactor_value(const kls_parallel_refactor_shared *shared,
                                       UF_long oldrow,
                                       double raw_value,
                                       double *value_out) {
  if (shared->scale <= 0) {
    *value_out = raw_value;
    return 1;
  }
  if (oldrow >= shared->n || shared->rs == NULL || shared->rs[oldrow] == 0.0) {
    return 0;
  }
  *value_out = raw_value / shared->rs[oldrow];
  return 1;
}

static int kls_parallel_refactor_mapped_value(
  const kls_parallel_refactor_shared *shared,
  UF_long input_pos,
  double *value_out) {
  if (input_pos >= shared->nnz) {
    return 0;
  }
  if (shared->scale <= 0) {
    *value_out = shared->values[input_pos];
    return 1;
  }
  return kls_parallel_refactor_value(shared, shared->row_idx[input_pos],
                                     shared->values[input_pos], value_out);
}

static int kls_checked_refactor_best_reject_candidate(
  const UF_long *rows,
  UF_long row_count,
  const double *x,
  double pivot,
  double tolerance,
  UF_long row_base,
  UF_long *rejected_row_out,
  UF_long *rejected_local_row_out,
  double *rejected_multiplier_abs_out,
  double *rejected_pivot_abs_out,
  double *rejected_candidate_abs_out) {
  if (rows == NULL || x == NULL || rejected_row_out == NULL ||
      rejected_local_row_out == NULL || rejected_multiplier_abs_out == NULL ||
      rejected_pivot_abs_out == NULL ||
      rejected_candidate_abs_out == NULL) {
    return 0;
  }

  const double pivot_abs = fabs(pivot);
  UF_long best_local_row = KLS_KLU_EMPTY;
  double best_candidate_abs = -1.0;
  double best_multiplier_abs = -1.0;
  for (UF_long p = 0; p < row_count; ++p) {
    const UF_long local_row = rows[p];
    const double candidate_abs = fabs(x[local_row]);
    const double multiplier_abs = fabs(x[local_row] / pivot);
    if (isfinite(multiplier_abs) &&
        multiplier_abs * tolerance <= 1.0 + 1.0e-12) {
      continue;
    }
    if (best_local_row == KLS_KLU_EMPTY ||
        (!isnan(candidate_abs) &&
         (isnan(best_candidate_abs) ||
          candidate_abs > best_candidate_abs))) {
      best_local_row = local_row;
      best_candidate_abs = candidate_abs;
      best_multiplier_abs = multiplier_abs;
    }
  }
  if (best_local_row == KLS_KLU_EMPTY) {
    return 0;
  }

  *rejected_row_out = row_base + best_local_row;
  *rejected_local_row_out = best_local_row;
  *rejected_multiplier_abs_out = best_multiplier_abs;
  *rejected_pivot_abs_out = pivot_abs;
  *rejected_candidate_abs_out = best_candidate_abs;
  return 1;
}

static void kls_parallel_refactor_block(kls_parallel_refactor_worker *worker,
                                        UF_long block) {
  kls_parallel_refactor_shared *shared = worker->shared;
  const UF_long *ap = shared->col_ptr;
  const UF_long *ai = shared->row_idx;
  const UF_long *map_col_ptr = shared->map_col_ptr;
  const UF_long *map_row_idx = shared->map_row_idx;
  const UF_long *map_input_pos = shared->map_input_pos;
  const UF_long *map_block_start = shared->map_block_start;
  const double *ax = shared->values;
  const trilinos_klu_l_symbolic *symbolic = shared->symbolic;
  trilinos_klu_l_numeric *numeric = shared->numeric;
  const UF_long *q = symbolic->Q;
  const UF_long *r = symbolic->R;
  const UF_long *pinv = numeric->Pinv;
  const UF_long *offp = numeric->Offp;
  double *offx = (double *)numeric->Offx;
  double *udiag = (double *)numeric->Udiag;
  double *x = worker->x;

  const UF_long k1 = r[block];
  const UF_long k2 = r[block + 1u];
  const UF_long nk = k2 - k1;

  if (nk == 1u) {
    const UF_long oldcol = q[k1];
    UF_long poff = offp[k1];
    const UF_long poff_end = offp[k1 + 1u];
    double s = 0.0;
    if (map_col_ptr != NULL && map_block_start != NULL) {
      for (UF_long p = map_col_ptr[k1]; p < map_block_start[k1]; ++p) {
        if (poff >= poff_end) {
          worker->invalid = 1;
          return;
        }
        const UF_long input_pos = map_input_pos[p];
        double value = 0.0;
        if (!kls_parallel_refactor_mapped_value(shared, input_pos, &value)) {
          worker->invalid = 1;
          return;
        }
        offx[poff++] = value;
      }
      for (UF_long p = map_block_start[k1]; p < map_col_ptr[k1 + 1u]; ++p) {
        if (map_row_idx[p] != k1) {
          worker->invalid = 1;
          return;
        }
        const UF_long input_pos = map_input_pos[p];
        if (!kls_parallel_refactor_mapped_value(shared, input_pos, &s)) {
          worker->invalid = 1;
          return;
        }
      }
    } else {
      const UF_long pend = ap[oldcol + 1u];
      for (UF_long p = ap[oldcol]; p < pend; ++p) {
        const UF_long oldrow = ai[p];
        double value = 0.0;
        if (!kls_parallel_refactor_value(shared, oldrow, ax[p], &value)) {
          worker->invalid = 1;
          return;
        }
        const UF_long newrow = pinv[oldrow];
        if (newrow < k1) {
          if (poff >= poff_end) {
            worker->invalid = 1;
            return;
          }
          offx[poff++] = value;
        } else if (newrow == k1) {
          s = value;
        } else {
          worker->invalid = 1;
          return;
        }
      }
    }
    udiag[k1] = s;
    if (s == 0.0) {
      kls_worker_record_singular(worker, k1, oldcol);
    }
    return;
  }

  UF_long *lip = numeric->Lip + k1;
  UF_long *llen = numeric->Llen + k1;
  UF_long *uip = numeric->Uip + k1;
  UF_long *ulen = numeric->Ulen + k1;
  double *lu = (double *)numeric->LUbx[block];
  if (lu == NULL) {
    worker->invalid = 1;
    return;
  }

  for (UF_long k = 0; k < nk; ++k) {
    const UF_long global_col = k + k1;
    const UF_long oldcol = q[global_col];
    UF_long poff = offp[global_col];
    const UF_long poff_end = offp[global_col + 1u];

    if (map_col_ptr != NULL && map_block_start != NULL) {
      for (UF_long p = map_col_ptr[global_col];
           p < map_block_start[global_col]; ++p) {
        if (poff >= poff_end) {
          worker->invalid = 1;
          return;
        }
        const UF_long input_pos = map_input_pos[p];
        double value = 0.0;
        if (!kls_parallel_refactor_mapped_value(shared, input_pos, &value)) {
          worker->invalid = 1;
          return;
        }
        offx[poff++] = value;
      }
      for (UF_long p = map_block_start[global_col];
           p < map_col_ptr[global_col + 1u]; ++p) {
        const UF_long global_row = map_row_idx[p];
        if (global_row < k1 || global_row >= k2) {
          worker->invalid = 1;
          return;
        }
        const UF_long input_pos = map_input_pos[p];
        double value = 0.0;
        if (!kls_parallel_refactor_mapped_value(shared, input_pos, &value)) {
          worker->invalid = 1;
          return;
        }
        x[global_row - k1] = value;
      }
    } else {
      const UF_long pend = ap[oldcol + 1u];
      for (UF_long p = ap[oldcol]; p < pend; ++p) {
        const UF_long oldrow = ai[p];
        double value = 0.0;
        if (!kls_parallel_refactor_value(shared, oldrow, ax[p], &value)) {
          worker->invalid = 1;
          return;
        }
        const UF_long global_row = pinv[oldrow];
        if (global_row < k1) {
          if (poff >= poff_end) {
            worker->invalid = 1;
            return;
          }
          offx[poff++] = value;
        } else if (global_row < k2) {
          x[global_row - k1] = value;
        } else {
          worker->invalid = 1;
          return;
        }
      }
    }

    UF_long *ui = NULL;
    double *ux = NULL;
    UF_long ucol_len = 0;
    kls_klu_get_pointer(lu, uip, ulen, k, &ui, &ux, &ucol_len);
    for (UF_long up = 0; up < ucol_len; ++up) {
      const UF_long j = ui[up];
      const double ujk = x[j];
      x[j] = 0.0;
      ux[up] = ujk;

      UF_long *li = NULL;
      double *lx = NULL;
      UF_long lcol_len = 0;
      kls_klu_get_pointer(lu, lip, llen, j, &li, &lx, &lcol_len);
      kls_scatter_subtract(x, li, lx, lcol_len, ujk);
    }

    const double ukk = x[k];
    x[k] = 0.0;
    if (ukk == 0.0) {
      kls_worker_record_singular(worker, global_col, q[global_col]);
      if (shared->halt_if_singular) {
        return;
      }
    }
    udiag[global_col] = ukk;

    UF_long *li = NULL;
    double *lx = NULL;
    UF_long lcol_len = 0;
    kls_klu_get_pointer(lu, lip, llen, k, &li, &lx, &lcol_len);
    UF_long rejected_row = KLS_KLU_EMPTY;
    UF_long rejected_local_row = KLS_KLU_EMPTY;
    double rejected_multiplier_abs = -1.0;
    double rejected_pivot_abs = -1.0;
    double rejected_candidate_abs = -1.0;
    if (shared->check_pivots &&
        kls_checked_refactor_best_reject_candidate(
          li, lcol_len, x, ukk, shared->pivot_tolerance, k1,
          &rejected_row, &rejected_local_row, &rejected_multiplier_abs,
          &rejected_pivot_abs, &rejected_candidate_abs)) {
      worker->pivot_rejected = 1;
      worker->rejected_pivot = global_col;
      worker->rejected_pivot_col = q[global_col];
      worker->rejected_row = rejected_row;
      worker->rejected_multiplier_abs = rejected_multiplier_abs;
      worker->rejected_pivot_abs = rejected_pivot_abs;
      worker->rejected_candidate_abs = rejected_candidate_abs;
      x[rejected_local_row] = 0.0;
      return;
    }
    for (UF_long p = 0; p < lcol_len; ++p) {
      const UF_long i = li[p];
      const double lij = x[i] / ukk;
      lx[p] = lij;
      x[i] = 0.0;
    }
  }
}

static void *kls_refactor_pool_worker_main(void *arg) {
  kls_parallel_refactor_worker *worker = (kls_parallel_refactor_worker *)arg;
  kls_refactor_pool *pool = worker->pool;
  kls_parallel_refactor_shared *shared = worker->shared;
  unsigned long seen_generation = 0;

  pthread_mutex_lock(&shared->lock);
  for (;;) {
    while (!pool->shutdown && pool->generation == seen_generation) {
      pthread_cond_wait(&pool->work_cond, &shared->lock);
    }
    if (pool->shutdown) {
      pthread_mutex_unlock(&shared->lock);
      return NULL;
    }
    seen_generation = pool->generation;
    pthread_mutex_unlock(&shared->lock);

    for (;;) {
      pthread_mutex_lock(&shared->lock);
      if (shared->stop || shared->next_block >= shared->symbolic->nblocks) {
        pthread_mutex_unlock(&shared->lock);
        break;
      }
      const UF_long begin = shared->next_block;
      UF_long end = begin + shared->block_chunk;
      if (end < begin || end > shared->symbolic->nblocks) {
        end = shared->symbolic->nblocks;
      }
      shared->next_block = end;
      pthread_mutex_unlock(&shared->lock);

      for (UF_long block = begin; block < end; ++block) {
        kls_parallel_refactor_block(worker, block);
        if (worker->invalid || worker->pivot_rejected ||
            (worker->singular && shared->halt_if_singular)) {
          pthread_mutex_lock(&shared->lock);
          shared->stop = 1;
          pthread_mutex_unlock(&shared->lock);
          break;
        }
        if (shared->block_done != NULL) {
          shared->block_done[block] = 1u;
        }
      }
      if (worker->invalid || worker->pivot_rejected ||
          (worker->singular && shared->halt_if_singular)) {
        break;
      }
    }

    pthread_mutex_lock(&shared->lock);
    pool->active_workers--;
    if (pool->active_workers == 0) {
      pthread_cond_signal(&pool->done_cond);
    }
  }
}

static void destroy_refactor_pool(kls_solver *solver) {
  if (solver == NULL || solver->refactor_pool == NULL) {
    return;
  }
  kls_refactor_pool *pool = solver->refactor_pool;
  if (pool->lock_initialized) {
    pthread_mutex_lock(&pool->shared.lock);
    pool->shutdown = 1;
    if (pool->conds_initialized) {
      pthread_cond_broadcast(&pool->work_cond);
    }
    pthread_mutex_unlock(&pool->shared.lock);
  }
  if (pool->threads != NULL) {
    for (int i = 0; i < pool->created_count; ++i) {
      pthread_join(pool->threads[i], NULL);
    }
  }
  if (pool->conds_initialized) {
    pthread_cond_destroy(&pool->work_cond);
    pthread_cond_destroy(&pool->done_cond);
  }
  if (pool->lock_initialized) {
    pthread_mutex_destroy(&pool->shared.lock);
  }
  if (pool->workers != NULL) {
    for (int i = 0; i < pool->thread_count; ++i) {
      free(pool->workers[i].x);
    }
  }
  free(pool->threads);
  free(pool->workers);
  free(pool->block_done);
  free(pool);
  solver->refactor_pool = NULL;
}

static int ensure_refactor_pool(kls_solver *solver, int thread_count) {
  if (solver == NULL || solver->symbolic == NULL || thread_count < 2) {
    return 0;
  }
  if (solver->refactor_pool != NULL &&
      (solver->refactor_pool->thread_count != thread_count ||
       solver->refactor_pool->maxblock < solver->symbolic->maxblock)) {
    destroy_refactor_pool(solver);
  }
  if (solver->refactor_pool != NULL) {
    return 1;
  }

  kls_refactor_pool *pool = (kls_refactor_pool *)calloc(1, sizeof(*pool));
  if (pool == NULL) {
    return 0;
  }
  pool->thread_count = thread_count;
  pool->maxblock = solver->symbolic->maxblock;
  pool->threads = (pthread_t *)calloc((size_t)thread_count, sizeof(*pool->threads));
  pool->workers =
    (kls_parallel_refactor_worker *)calloc((size_t)thread_count, sizeof(*pool->workers));
  if (pool->threads == NULL || pool->workers == NULL) {
    free(pool->threads);
    free(pool->workers);
    free(pool);
    return 0;
  }
  if (pthread_mutex_init(&pool->shared.lock, NULL) != 0) {
    free(pool->threads);
    free(pool->workers);
    free(pool);
    return 0;
  }
  pool->lock_initialized = 1;
  if (pthread_cond_init(&pool->work_cond, NULL) != 0) {
    pthread_mutex_destroy(&pool->shared.lock);
    free(pool->threads);
    free(pool->workers);
    free(pool);
    return 0;
  }
  if (pthread_cond_init(&pool->done_cond, NULL) != 0) {
    pthread_cond_destroy(&pool->work_cond);
    pthread_mutex_destroy(&pool->shared.lock);
    free(pool->threads);
    free(pool->workers);
    free(pool);
    return 0;
  }
  pool->conds_initialized = 1;

  for (int i = 0; i < thread_count; ++i) {
    pool->workers[i].shared = &pool->shared;
    pool->workers[i].pool = pool;
    pool->workers[i].x =
      (double *)calloc((size_t)solver->symbolic->maxblock, sizeof(double));
    if (pool->workers[i].x == NULL ||
        pthread_create(&pool->threads[i], NULL, kls_refactor_pool_worker_main,
                       &pool->workers[i]) != 0) {
      break;
    }
    pool->created_count++;
  }
  if (pool->created_count != thread_count) {
    solver->refactor_pool = pool;
    destroy_refactor_pool(solver);
    return 0;
  }
  solver->refactor_pool = pool;
  return 1;
}

static UF_long kls_refactor_pool_block_chunk(const kls_solver *solver,
                                             int thread_count) {
  if (solver == NULL || solver->symbolic == NULL || thread_count <= 0) {
    return 1u;
  }
  const UF_long nblocks = solver->symbolic->nblocks;
  if (nblocks < 4096u ||
      (double)solver->symbolic->maxblock > 0.50 * (double)solver->n) {
    return 1u;
  }
  const UF_long target_ranges = (UF_long)thread_count * 4096u;
  UF_long chunk = target_ranges > 0u ? nblocks / target_ranges : 1u;
  if (chunk < 1u) {
    chunk = 1u;
  }
  if (chunk > 64u) {
    chunk = 64u;
  }
  return chunk;
}

static int kls_refactor_pool_rejected_prefix_current(
  const kls_solver *solver,
  const unsigned char *block_done,
  UF_long rejected_pivot) {
  if (solver == NULL || solver->symbolic == NULL ||
      solver->symbolic->R == NULL || block_done == NULL ||
      rejected_pivot == KLS_KLU_EMPTY) {
    return 0;
  }
  const UF_long block = kls_block_for_pivot(solver, rejected_pivot);
  if (block == KLS_KLU_EMPTY) {
    return 0;
  }
  for (UF_long b = 0; b < block; ++b) {
    if (!block_done[b]) {
      return 0;
    }
  }
  return 1;
}

static int run_refactor_pool(kls_solver *solver,
                             double *numeric_values,
                             int thread_count,
                             int check_pivots,
                             int *invalid_out,
                             int *pivot_rejected_out,
                             UF_long *rejected_pivot_out,
                             UF_long *rejected_pivot_col_out,
                             UF_long *rejected_row_out,
                             double *rejected_multiplier_abs_out,
                             double *rejected_pivot_abs_out,
                             double *rejected_candidate_abs_out,
                             int *singular_out,
                             UF_long *numerical_rank_out,
                             UF_long *singular_col_out,
                             int *prefix_current_out) {
  if (!ensure_refactor_pool(solver, thread_count)) {
    return 0;
  }
  kls_refactor_pool *pool = solver->refactor_pool;
  kls_parallel_refactor_shared *shared = &pool->shared;
  if (prefix_current_out != NULL) {
    *prefix_current_out = 0;
  }

  pthread_mutex_lock(&shared->lock);
  if (pool->active_workers != 0) {
    pthread_mutex_unlock(&shared->lock);
    return 0;
  }
  if (pool->block_done_capacity < solver->symbolic->nblocks) {
    unsigned char *done =
      (unsigned char *)realloc(pool->block_done,
                               (size_t)solver->symbolic->nblocks *
                                 sizeof(*pool->block_done));
    if (done == NULL) {
      pthread_mutex_unlock(&shared->lock);
      return 0;
    }
    pool->block_done = done;
    pool->block_done_capacity = solver->symbolic->nblocks;
  }
  if (pool->block_done != NULL && solver->symbolic->nblocks > 0u) {
    memset(pool->block_done, 0,
           (size_t)solver->symbolic->nblocks * sizeof(*pool->block_done));
  }

  shared->col_ptr = solver->col_ptr;
  shared->row_idx = solver->row_idx;
  shared->nnz = solver->nnz;
  shared->map_col_ptr = solver->refactor_col_ptr;
  shared->map_row_idx = solver->refactor_row_idx;
  shared->map_input_pos = solver->refactor_input_pos;
  shared->map_block_start = solver->refactor_block_start;
  shared->values = numeric_values;
  shared->symbolic = solver->symbolic;
  shared->numeric = solver->numeric;
  shared->n = solver->n;
  shared->rs = solver->numeric->Rs;
  shared->scale = (int)solver->common.scale;
  shared->halt_if_singular = solver->common.halt_if_singular;
  shared->check_pivots = check_pivots;
  shared->pivot_tolerance = solver->common.tol;
  shared->next_block = 0;
  shared->block_chunk = kls_refactor_pool_block_chunk(solver, thread_count);
  shared->block_done = pool->block_done;
  shared->stop = 0;
  pool->active_workers = thread_count;

  for (int i = 0; i < thread_count; ++i) {
    kls_parallel_refactor_worker *worker = &pool->workers[i];
    worker->invalid = 0;
    worker->pivot_rejected = 0;
    worker->singular = 0;
    worker->rejected_pivot = KLS_KLU_EMPTY;
    worker->rejected_pivot_col = KLS_KLU_EMPTY;
    worker->rejected_row = KLS_KLU_EMPTY;
    worker->rejected_multiplier_abs = -1.0;
    worker->rejected_pivot_abs = -1.0;
    worker->rejected_candidate_abs = -1.0;
    worker->numerical_rank = UF_long_max;
    worker->singular_col = KLS_KLU_EMPTY;
    /* Successful KLU-style refactors clear touched X slots as they go. */
    if (pool->scratch_dirty) {
      memset(worker->x, 0, (size_t)solver->symbolic->maxblock * sizeof(*worker->x));
    }
  }
  pool->scratch_dirty = 0;

  pool->generation++;
  pthread_cond_broadcast(&pool->work_cond);
  while (pool->active_workers > 0) {
    pthread_cond_wait(&pool->done_cond, &shared->lock);
  }
  pthread_mutex_unlock(&shared->lock);

  int invalid = 0;
  int pivot_rejected = 0;
  UF_long rejected_pivot = KLS_KLU_EMPTY;
  UF_long rejected_pivot_col = KLS_KLU_EMPTY;
  UF_long rejected_row = KLS_KLU_EMPTY;
  double rejected_multiplier_abs = -1.0;
  double rejected_pivot_abs = -1.0;
  double rejected_candidate_abs = -1.0;
  int singular = 0;
  UF_long numerical_rank = UF_long_max;
  UF_long singular_col = KLS_KLU_EMPTY;
  for (int i = 0; i < thread_count; ++i) {
    const kls_parallel_refactor_worker *worker = &pool->workers[i];
    if (worker->invalid) {
      invalid = 1;
    }
    if (worker->pivot_rejected) {
      pivot_rejected = 1;
      if (worker->rejected_pivot != KLS_KLU_EMPTY &&
          (rejected_pivot == KLS_KLU_EMPTY ||
           worker->rejected_pivot < rejected_pivot)) {
        rejected_pivot = worker->rejected_pivot;
        rejected_pivot_col = worker->rejected_pivot_col;
        rejected_row = worker->rejected_row;
        rejected_multiplier_abs = worker->rejected_multiplier_abs;
        rejected_pivot_abs = worker->rejected_pivot_abs;
        rejected_candidate_abs = worker->rejected_candidate_abs;
      }
    }
    if (worker->singular && worker->numerical_rank < numerical_rank) {
      singular = 1;
      numerical_rank = worker->numerical_rank;
      singular_col = worker->singular_col;
    }
  }

  *invalid_out = invalid;
  *pivot_rejected_out = pivot_rejected;
  if (rejected_pivot_out != NULL) {
    *rejected_pivot_out = rejected_pivot;
  }
  if (rejected_pivot_col_out != NULL) {
    *rejected_pivot_col_out = rejected_pivot_col;
  }
  if (rejected_row_out != NULL) {
    *rejected_row_out = rejected_row;
  }
  if (rejected_multiplier_abs_out != NULL) {
    *rejected_multiplier_abs_out = rejected_multiplier_abs;
  }
  if (rejected_pivot_abs_out != NULL) {
    *rejected_pivot_abs_out = rejected_pivot_abs;
  }
  if (rejected_candidate_abs_out != NULL) {
    *rejected_candidate_abs_out = rejected_candidate_abs;
  }
  *singular_out = singular;
  *numerical_rank_out = numerical_rank;
  *singular_col_out = singular_col;
  if (prefix_current_out != NULL && pivot_rejected && !invalid) {
    *prefix_current_out =
      kls_refactor_pool_rejected_prefix_current(solver, pool->block_done,
                                                rejected_pivot);
  }
  pool->scratch_dirty = invalid || pivot_rejected ||
                        (singular && solver->common.halt_if_singular);
  return 1;
}

static int64_t read_index(kls_index_type type, const void *data, int64_t i) {
  if (type == KLS_INDEX_INT32) {
    const int32_t *p = (const int32_t *)data;
    return (int64_t)p[i];
  }
  const int64_t *p = (const int64_t *)data;
  return p[i];
}

static void free_symbolic(kls_solver *solver) {
  if (solver->symbolic != NULL) {
    trilinos_klu_l_free_symbolic(&solver->symbolic, &solver->common);
    solver->symbolic = NULL;
  }
}

static void free_numeric(kls_solver *solver) {
  destroy_egraph_refactor_pool(solver);
  free_egraph_worker_scratch(solver);
  free_egraph_pipeline_done(solver);
  free_row_refactor_pattern(solver);
  free_fast_reject_tail_plan(solver);
  free_refactor_lu_pointer_cache(solver);
  if (solver->numeric != NULL) {
    free_refactor_map(solver);
    free_refactor_schedule(solver);
    trilinos_klu_l_free_numeric(&solver->numeric, &solver->common);
    solver->numeric = NULL;
  }
}

static void clear_matrix(kls_solver *solver) {
  destroy_refactor_pool(solver);
  free_numeric(solver);
  free_egraph_worker_scratch(solver);
  free_egraph_pipeline_done(solver);
  free_symbolic(solver);
  free(solver->col_ptr);
  free(solver->row_idx);
  free(solver->input_to_csc);
  free(solver->row_perm);
  free(solver->row_scale);
  free(solver->col_scale);
  free(solver->values);
  free(solver->solve_perm_workspace);
  free_refactor_map(solver);
  free_refactor_schedule(solver);
  solver->col_ptr = NULL;
  solver->row_idx = NULL;
  solver->input_to_csc = NULL;
  solver->row_perm = NULL;
  solver->row_scale = NULL;
  solver->col_scale = NULL;
  solver->values = NULL;
  solver->solve_perm_workspace = NULL;
  solver->solve_perm_workspace_n = 0;
  solver->n = 0;
  solver->nnz = 0;
  solver->input_format = KLS_INPUT_NONE;
  solver->orientation = KLS_ORIENTATION_NORMAL;
  solver->auto_metis_checked = 0;
  solver->auto_pivot_checked = 0;
  solver->auto_scale_checked = 0;
  solver->exact_matching_selected = 0;
  solver->spral_matching_selected = 0;
  solver->fast_block_restarts = 0;
  solver->fast_tail_restarts = 0;
  memset(&solver->stats, 0, sizeof(solver->stats));
  solver->stats.struct_size = sizeof(solver->stats);
  kls_clear_fast_reject_stats(solver);
}

static int initial_scale(const kls_options *options) {
  return options->scale == KLS_SCALE_AUTO ? 2 : options->scale;
}

static int compare_double(const void *a, const void *b) {
  const double left = *(const double *)a;
  const double right = *(const double *)b;
  return (left > right) - (left < right);
}

static int is_large_low_degree_diagonal_pattern(UF_long n,
                                                const UF_long *col_ptr,
                                                const UF_long *row_idx) {
  if (n < 50000 || col_ptr == NULL || row_idx == NULL ||
      n > UF_long_max / 8 || col_ptr[n] > 8u * n) {
    return 0;
  }

  UF_long *row_degree = (UF_long *)calloc((size_t)n, sizeof(*row_degree));
  if (row_degree == NULL) {
    return 0;
  }

  UF_long diagonal_count = 0;
  UF_long max_col_degree = 0;
  int low_degree = 1;
  for (UF_long col = 0; col < n && low_degree; ++col) {
    const UF_long col_degree = col_ptr[col + 1u] - col_ptr[col];
    if (col_degree > max_col_degree) {
      max_col_degree = col_degree;
    }
    if (col_degree > 64u) {
      low_degree = 0;
      break;
    }
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      const UF_long row = row_idx[p];
      if (row >= n || row_degree[row] >= 64u) {
        low_degree = 0;
        break;
      }
      row_degree[row]++;
      if (row == col) {
        diagonal_count++;
      }
    }
  }

  UF_long max_row_degree = 0;
  for (UF_long row = 0; row < n && low_degree; ++row) {
    if (row_degree[row] == 0) {
      low_degree = 0;
      break;
    }
    if (row_degree[row] > max_row_degree) {
      max_row_degree = row_degree[row];
    }
  }
  free(row_degree);

  return low_degree && max_col_degree <= 64u && max_row_degree <= 64u &&
         200.0 * (double)diagonal_count >= 199.0 * (double)n;
}

static int is_medium_low_degree_full_diagonal_pattern(UF_long n,
                                                      const UF_long *col_ptr,
                                                      const UF_long *row_idx) {
  if (n < 15000 || n > 25000 || col_ptr == NULL || row_idx == NULL ||
      n > UF_long_max / 8 || col_ptr[n] > 8u * n) {
    return 0;
  }

  UF_long *row_degree = (UF_long *)calloc((size_t)n, sizeof(*row_degree));
  if (row_degree == NULL) {
    return 0;
  }

  UF_long diagonal_columns = 0;
  UF_long max_col_degree = 0;
  int low_degree = 1;
  for (UF_long col = 0; col < n && low_degree; ++col) {
    const UF_long col_degree = col_ptr[col + 1u] - col_ptr[col];
    int has_diagonal = 0;
    if (col_degree > max_col_degree) {
      max_col_degree = col_degree;
    }
    if (col_degree > 64u) {
      low_degree = 0;
      break;
    }
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      const UF_long row = row_idx[p];
      if (row >= n || row_degree[row] >= 64u) {
        low_degree = 0;
        break;
      }
      row_degree[row]++;
      if (row == col) {
        has_diagonal = 1;
      }
    }
    if (has_diagonal) {
      diagonal_columns++;
    }
  }

  UF_long max_row_degree = 0;
  for (UF_long row = 0; row < n && low_degree; ++row) {
    if (row_degree[row] == 0) {
      low_degree = 0;
      break;
    }
    if (row_degree[row] > max_row_degree) {
      max_row_degree = row_degree[row];
    }
  }
  free(row_degree);

  return low_degree && max_col_degree <= 64u && max_row_degree <= 64u &&
         diagonal_columns == n;
}

static int is_medium_dense_diagonal_high_degree_pattern(UF_long n,
                                                        const UF_long *col_ptr,
                                                        const UF_long *row_idx) {
  if (n < 30000 || n > 45000 || col_ptr == NULL || row_idx == NULL ||
      n > UF_long_max / 16 || col_ptr[n] < 10u * n || col_ptr[n] > 16u * n) {
    return 0;
  }

  UF_long *row_degree = (UF_long *)calloc((size_t)n, sizeof(*row_degree));
  if (row_degree == NULL) {
    return 0;
  }

  UF_long diagonal_count = 0;
  UF_long max_col_degree = 0;
  int valid = 1;
  for (UF_long col = 0; col < n && valid; ++col) {
    const UF_long col_degree = col_ptr[col + 1u] - col_ptr[col];
    if (col_degree > max_col_degree) {
      max_col_degree = col_degree;
    }
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      const UF_long row = row_idx[p];
      if (row >= n || row_degree[row] == UF_long_max) {
        valid = 0;
        break;
      }
      row_degree[row]++;
      if (row == col) {
        diagonal_count++;
      }
    }
  }

  UF_long max_row_degree = 0;
  int no_empty_rows = 1;
  for (UF_long row = 0; row < n && valid; ++row) {
    if (row_degree[row] == 0) {
      no_empty_rows = 0;
      break;
    }
    if (row_degree[row] > max_row_degree) {
      max_row_degree = row_degree[row];
    }
  }
  free(row_degree);

  return valid && no_empty_rows &&
         max_col_degree >= 512u && max_row_degree >= 512u &&
         max_col_degree <= 4096u && max_row_degree <= 4096u &&
         1000.0 * (double)diagonal_count >= 995.0 * (double)n;
}

#ifdef KLS_HAVE_METIS
static int is_medium_spiked_low_diagonal_pattern(UF_long n,
                                                 const UF_long *col_ptr,
                                                 const UF_long *row_idx) {
  if (n < 50000 || n > 125000 || col_ptr == NULL || row_idx == NULL ||
      n > UF_long_max / 32 || col_ptr[n] < 20u * n || col_ptr[n] > 32u * n) {
    return 0;
  }

  UF_long *row_degree = (UF_long *)calloc((size_t)n, sizeof(*row_degree));
  if (row_degree == NULL) {
    return 0;
  }

  UF_long diagonal_count = 0;
  UF_long max_col_degree = 0;
  int valid = 1;
  for (UF_long col = 0; col < n && valid; ++col) {
    const UF_long col_degree = col_ptr[col + 1u] - col_ptr[col];
    if (col_degree > max_col_degree) {
      max_col_degree = col_degree;
    }
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      const UF_long row = row_idx[p];
      if (row >= n || row_degree[row] == UF_long_max) {
        valid = 0;
        break;
      }
      row_degree[row]++;
      if (row == col) {
        diagonal_count++;
      }
    }
  }

  UF_long max_row_degree = 0;
  int no_empty_rows = 1;
  for (UF_long row = 0; row < n && valid; ++row) {
    if (row_degree[row] == 0) {
      no_empty_rows = 0;
      break;
    }
    if (row_degree[row] > max_row_degree) {
      max_row_degree = row_degree[row];
    }
  }
  free(row_degree);

  return valid && no_empty_rows &&
         100.0 * (double)diagonal_count >= 20.0 * (double)n &&
         100.0 * (double)diagonal_count <= 35.0 * (double)n &&
         10.0 * (double)max_col_degree >= 4.0 * (double)n &&
         10.0 * (double)max_row_degree >= 4.0 * (double)n &&
         10.0 * (double)max_col_degree <= 6.0 * (double)n &&
         10.0 * (double)max_row_degree <= 6.0 * (double)n;
}

static int is_small_spiked_low_diagonal_pattern(UF_long n,
                                                const UF_long *col_ptr,
                                                const UF_long *row_idx) {
  if (n < 2000 || n > 20000 || col_ptr == NULL || row_idx == NULL ||
      n > UF_long_max / 14 || col_ptr[n] < 8u * n || col_ptr[n] > 14u * n) {
    return 0;
  }

  UF_long *row_degree = (UF_long *)calloc((size_t)n, sizeof(*row_degree));
  if (row_degree == NULL) {
    return 0;
  }

  UF_long diagonal_count = 0;
  UF_long max_col_degree = 0;
  int valid = 1;
  for (UF_long col = 0; col < n && valid; ++col) {
    const UF_long col_degree = col_ptr[col + 1u] - col_ptr[col];
    if (col_degree > max_col_degree) {
      max_col_degree = col_degree;
    }
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      const UF_long row = row_idx[p];
      if (row >= n || row_degree[row] == UF_long_max) {
        valid = 0;
        break;
      }
      row_degree[row]++;
      if (row == col) {
        diagonal_count++;
      }
    }
  }

  UF_long max_row_degree = 0;
  int no_empty_rows = 1;
  for (UF_long row = 0; row < n && valid; ++row) {
    if (row_degree[row] == 0) {
      no_empty_rows = 0;
      break;
    }
    if (row_degree[row] > max_row_degree) {
      max_row_degree = row_degree[row];
    }
  }
  free(row_degree);

  return valid && no_empty_rows &&
         100.0 * (double)diagonal_count >= 20.0 * (double)n &&
         100.0 * (double)diagonal_count <= 35.0 * (double)n &&
         10.0 * (double)max_col_degree >= 4.0 * (double)n &&
         10.0 * (double)max_row_degree >= 4.0 * (double)n &&
         10.0 * (double)max_col_degree <= 6.0 * (double)n &&
         10.0 * (double)max_row_degree <= 6.0 * (double)n;
}
#endif

static int is_large_diagonal_circuit_like_pattern(UF_long n,
                                                  const UF_long *col_ptr,
                                                  const UF_long *row_idx) {
  if (n < 35000 || col_ptr == NULL || row_idx == NULL ||
      n > UF_long_max / 16 || col_ptr[n] > 16u * n) {
    return 0;
  }

  UF_long *row_degree = (UF_long *)calloc((size_t)n, sizeof(*row_degree));
  if (row_degree == NULL) {
    return 0;
  }

  UF_long diagonal_count = 0;
  int valid = 1;
  for (UF_long col = 0; col < n && valid; ++col) {
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      const UF_long row = row_idx[p];
      if (row >= n || row_degree[row] == UF_long_max) {
        valid = 0;
        break;
      }
      row_degree[row]++;
      if (row == col) {
        diagonal_count++;
      }
    }
  }

  int no_empty_rows = 1;
  for (UF_long row = 0; row < n && valid; ++row) {
    if (row_degree[row] == 0) {
      no_empty_rows = 0;
      break;
    }
  }
  free(row_degree);

  return valid && no_empty_rows &&
         10.0 * (double)diagonal_count >= 9.0 * (double)n;
}

static int is_medium_sparse_high_degree_diagonal_pattern(UF_long n,
                                                         const UF_long *col_ptr,
                                                         const UF_long *row_idx) {
  if (n < 35000 || n > 45000 || col_ptr == NULL || row_idx == NULL ||
      n > UF_long_max / 8 || col_ptr[n] > 6u * n) {
    return 0;
  }

  UF_long *row_degree = (UF_long *)calloc((size_t)n, sizeof(*row_degree));
  if (row_degree == NULL) {
    return 0;
  }

  UF_long diagonal_count = 0;
  UF_long max_col_degree = 0;
  int valid = 1;
  for (UF_long col = 0; col < n && valid; ++col) {
    const UF_long col_degree = col_ptr[col + 1u] - col_ptr[col];
    if (col_degree > max_col_degree) {
      max_col_degree = col_degree;
    }
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      const UF_long row = row_idx[p];
      if (row >= n || row_degree[row] == UF_long_max) {
        valid = 0;
        break;
      }
      row_degree[row]++;
      if (row == col) {
        diagonal_count++;
      }
    }
  }

  UF_long max_row_degree = 0;
  int no_empty_rows = 1;
  for (UF_long row = 0; row < n && valid; ++row) {
    if (row_degree[row] == 0) {
      no_empty_rows = 0;
      break;
    }
    if (row_degree[row] > max_row_degree) {
      max_row_degree = row_degree[row];
    }
  }
  free(row_degree);

  return valid && no_empty_rows &&
         max_col_degree >= 512u && max_row_degree >= 512u &&
         max_col_degree <= 4096u && max_row_degree <= 4096u &&
         100.0 * (double)diagonal_count >= 90.0 * (double)n &&
         100.0 * (double)diagonal_count < 99.0 * (double)n;
}

#ifdef KLS_HAVE_METIS
static int is_large_very_low_degree_full_diagonal_pattern(UF_long n,
                                                          const UF_long *col_ptr,
                                                          const UF_long *row_idx) {
  if (n < 100000 || col_ptr == NULL || row_idx == NULL ||
      n > UF_long_max / 5 || col_ptr[n] > 5u * n) {
    return 0;
  }

  UF_long *row_degree = (UF_long *)calloc((size_t)n, sizeof(*row_degree));
  if (row_degree == NULL) {
    return 0;
  }

  UF_long diagonal_count = 0;
  UF_long max_col_degree = 0;
  int very_low_degree = 1;
  for (UF_long col = 0; col < n && very_low_degree; ++col) {
    const UF_long col_degree = col_ptr[col + 1] - col_ptr[col];
    if (col_degree > max_col_degree) {
      max_col_degree = col_degree;
    }
    if (col_degree > 8) {
      very_low_degree = 0;
      break;
    }
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1]; ++p) {
      const UF_long row = row_idx[p];
      if (row >= n || row_degree[row] >= 8) {
        very_low_degree = 0;
        break;
      }
      row_degree[row]++;
      if (row == col) {
        diagonal_count++;
      }
    }
  }

  UF_long max_row_degree = 0;
  for (UF_long row = 0; row < n && very_low_degree; ++row) {
    if (row_degree[row] == 0) {
      very_low_degree = 0;
      break;
    }
    if (row_degree[row] > max_row_degree) {
      max_row_degree = row_degree[row];
    }
  }
  free(row_degree);

  return very_low_degree && max_col_degree <= 8 && max_row_degree <= 8 &&
         diagonal_count == n;
}

static int is_large_diagonal_metis_start_pattern(UF_long n,
                                                 const UF_long *col_ptr,
                                                 const UF_long *row_idx) {
  if (n < 90000 || col_ptr == NULL || row_idx == NULL ||
      n > UF_long_max / 10) {
    return 0;
  }

  UF_long *row_degree = (UF_long *)calloc((size_t)n, sizeof(*row_degree));
  if (row_degree == NULL) {
    return 0;
  }

  UF_long diagonal_count = 0;
  UF_long max_col_degree = 0;
  int valid = 1;
  for (UF_long col = 0; col < n && valid; ++col) {
    const UF_long col_degree = col_ptr[col + 1] - col_ptr[col];
    if (col_degree > max_col_degree) {
      max_col_degree = col_degree;
    }
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1]; ++p) {
      const UF_long row = row_idx[p];
      if (row >= n || row_degree[row] == UF_long_max) {
        valid = 0;
        break;
      }
      row_degree[row]++;
      if (row == col) {
        diagonal_count++;
      }
    }
  }

  UF_long max_row_degree = 0;
  for (UF_long row = 0; row < n && valid; ++row) {
    if (row_degree[row] == 0) {
      valid = 0;
      break;
    }
    if (row_degree[row] > max_row_degree) {
      max_row_degree = row_degree[row];
    }
  }
  free(row_degree);

  if (!valid) {
    return 0;
  }

  const int full_sparse =
    diagonal_count == n && col_ptr[n] <= 6u * n &&
    max_col_degree <= 512u && max_row_degree <= 512u &&
    (max_col_degree > 64u || max_row_degree > 64u);
  const int dense_spike_density =
    col_ptr[n] >= 8u * n || (n >= 250000 && col_ptr[n] >= 6u * n);
  const int dense_spike =
    1000.0 * (double)diagonal_count >= 999.0 * (double)n &&
    dense_spike_density &&
    col_ptr[n] <= 10u * n &&
    (2.0 * (double)max_col_degree >= (double)n ||
     2.0 * (double)max_row_degree >= (double)n);
  return full_sparse || dense_spike;
}

static int is_large_nearly_diagonal_spiked_metis_pattern(
  UF_long n,
  const UF_long *col_ptr,
  const UF_long *row_idx) {
  if (n < 200000 || col_ptr == NULL || row_idx == NULL ||
      n > UF_long_max / 12 || col_ptr[n] < 6u * n ||
      col_ptr[n] > 12u * n) {
    return 0;
  }

  UF_long *row_degree = (UF_long *)calloc((size_t)n, sizeof(*row_degree));
  if (row_degree == NULL) {
    return 0;
  }

  UF_long diagonal_count = 0;
  UF_long max_col_degree = 0;
  int valid = 1;
  for (UF_long col = 0; col < n && valid; ++col) {
    const UF_long col_degree = col_ptr[col + 1u] - col_ptr[col];
    if (col_degree > max_col_degree) {
      max_col_degree = col_degree;
    }
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      const UF_long row = row_idx[p];
      if (row >= n || row_degree[row] == UF_long_max) {
        valid = 0;
        break;
      }
      row_degree[row]++;
      if (row == col) {
        diagonal_count++;
      }
    }
  }

  UF_long max_row_degree = 0;
  int no_empty_rows = 1;
  for (UF_long row = 0; row < n && valid; ++row) {
    if (row_degree[row] == 0) {
      no_empty_rows = 0;
      break;
    }
    if (row_degree[row] > max_row_degree) {
      max_row_degree = row_degree[row];
    }
  }
  free(row_degree);
  if (!valid || !no_empty_rows) {
    return 0;
  }

  const double nnz = (double)col_ptr[n];
  const double avg_degree = nnz / (double)n;
  const double max_degree =
    (double)(max_col_degree > max_row_degree ? max_col_degree : max_row_degree);
  const int sparse_spike =
    1000.0 * (double)diagonal_count >= 995.0 * (double)n &&
    avg_degree >= 6.0 && avg_degree <= 8.0 &&
    max_degree >= 0.05 * (double)n &&
    max_degree <= 0.25 * (double)n;
  const int dense_spike =
    1000.0 * (double)diagonal_count >= 980.0 * (double)n &&
    avg_degree >= 8.0 && avg_degree <= 12.0 &&
    max_degree >= 0.50 * (double)n;
  return sparse_spike || dense_spike;
}

static int is_large_sparse_diagonal_low_degree_pattern(UF_long n,
                                                       const UF_long *col_ptr,
                                                       const UF_long *row_idx) {
  if (n < 150000 || n > 250000 || col_ptr == NULL || row_idx == NULL ||
      n > UF_long_max / 8 || col_ptr[n] < 5u * n || col_ptr[n] > 8u * n) {
    return 0;
  }

  UF_long *row_degree = (UF_long *)calloc((size_t)n, sizeof(*row_degree));
  if (row_degree == NULL) {
    return 0;
  }

  UF_long diagonal_count = 0;
  UF_long max_col_degree = 0;
  int valid = 1;
  for (UF_long col = 0; col < n && valid; ++col) {
    const UF_long col_degree = col_ptr[col + 1u] - col_ptr[col];
    if (col_degree > max_col_degree) {
      max_col_degree = col_degree;
    }
    if (col_degree > 64u) {
      valid = 0;
      break;
    }
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      const UF_long row = row_idx[p];
      if (row >= n || row_degree[row] >= 64u) {
        valid = 0;
        break;
      }
      row_degree[row]++;
      if (row == col) {
        diagonal_count++;
      }
    }
  }

  UF_long max_row_degree = 0;
  int no_empty_rows = 1;
  for (UF_long row = 0; row < n && valid; ++row) {
    if (row_degree[row] == 0) {
      no_empty_rows = 0;
      break;
    }
    if (row_degree[row] > max_row_degree) {
      max_row_degree = row_degree[row];
    }
  }
  free(row_degree);

  return valid && no_empty_rows && max_col_degree <= 64u &&
         max_row_degree <= 64u &&
         100.0 * (double)diagonal_count >= 5.0 * (double)n &&
         100.0 * (double)diagonal_count <= 20.0 * (double)n;
}

static int is_medium_bounded_degree_diagonal_pattern(UF_long n,
                                                     const UF_long *col_ptr,
                                                     const UF_long *row_idx) {
  if (n < 7000 || n > 12000 || col_ptr == NULL || row_idx == NULL ||
      col_ptr[n] > (UF_long)(5 * n)) {
    return 0;
  }

  UF_long *row_degree = (UF_long *)calloc((size_t)n, sizeof(*row_degree));
  if (row_degree == NULL) {
    return 0;
  }

  UF_long diagonal_count = 0;
  UF_long max_col_degree = 0;
  int low_degree = 1;
  for (UF_long col = 0; col < n && low_degree; ++col) {
    const UF_long col_degree = col_ptr[col + 1] - col_ptr[col];
    if (col_degree > max_col_degree) {
      max_col_degree = col_degree;
    }
    if (col_degree > 128) {
      low_degree = 0;
      break;
    }
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1]; ++p) {
      const UF_long row = row_idx[p];
      if (row >= n || row_degree[row] >= 128) {
        low_degree = 0;
        break;
      }
      row_degree[row]++;
      if (row == col) {
        diagonal_count++;
      }
    }
  }

  UF_long max_row_degree = 0;
  for (UF_long row = 0; row < n && low_degree; ++row) {
    if (row_degree[row] > max_row_degree) {
      max_row_degree = row_degree[row];
    }
  }
  free(row_degree);

  return low_degree && max_col_degree <= 128 && max_row_degree <= 128 &&
         1000.0 * (double)diagonal_count >= 990.0 * (double)n;
}
#endif

static int symbolic_is_low_work_dominant_btf(
  UF_long n,
  const trilinos_klu_l_symbolic *symbolic);

static int choose_auto_scale_from_pattern(UF_long n,
                                          const UF_long *col_ptr,
                                          const UF_long *row_idx,
                                          const kls_options *options,
                                          const double *numeric_values) {
  if (options == NULL || col_ptr == NULL || row_idx == NULL || numeric_values == NULL ||
      options->scale != KLS_SCALE_AUTO || n <= 0) {
    return options == NULL ? 2 : initial_scale(options);
  }

  if (is_large_low_degree_diagonal_pattern(n, col_ptr, row_idx)) {
    return -1;
  }
#ifdef KLS_HAVE_METIS
  if (is_large_diagonal_metis_start_pattern(n, col_ptr, row_idx)) {
    return -1;
  }
  if (is_large_sparse_diagonal_low_degree_pattern(n, col_ptr, row_idx)) {
    return 1;
  }
  if (is_medium_dense_diagonal_high_degree_pattern(n, col_ptr, row_idx)) {
    return -1;
  }
  if (is_medium_spiked_low_diagonal_pattern(n, col_ptr, row_idx)) {
    return 1;
  }
  if (is_small_spiked_low_diagonal_pattern(n, col_ptr, row_idx)) {
    return 0;
  }
  if (is_medium_bounded_degree_diagonal_pattern(n, col_ptr, row_idx)) {
    return -1;
  }
#endif

  double *row_max = (double *)calloc((size_t)n, sizeof(*row_max));
  double *diag_max = (double *)calloc((size_t)n, sizeof(*diag_max));
  if (row_max == NULL || diag_max == NULL) {
    free(row_max);
    free(diag_max);
    return initial_scale(options);
  }

  /* Cheap first-factor policy:
     - keep well-populated, moderately scaled circuit diagonals unscaled;
     - use sum scaling for sparse diagonals whose row magnitudes are balanced. */
  const double diag_spread_limit = 1.0e7;
  const double min_diag_fraction = 0.90;
  const double row_spread_limit = 1.0e4;
  const double row_p90_p10_limit = 5.5;
  UF_long diag_count = 0;
  UF_long diag_unit_count = 0;
  double min_diag = DBL_MAX;
  double max_diag = 0.0;

  for (UF_long col = 0; col < n; ++col) {
    double diag_abs = 0.0;
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1]; ++p) {
      if (row_idx[p] == col) {
        const double value_abs = fabs(numeric_values[p]);
        if (isfinite(value_abs) && value_abs > diag_abs) {
          diag_abs = value_abs;
        }
        if (isfinite(value_abs) && value_abs > diag_max[(size_t)col]) {
          diag_max[(size_t)col] = value_abs;
        }
      }
      const double value_abs = fabs(numeric_values[p]);
      if (isfinite(value_abs) && value_abs > row_max[(size_t)row_idx[p]]) {
        row_max[(size_t)row_idx[p]] = value_abs;
      }
    }
    if (diag_abs > 0.0) {
      diag_count++;
      if (diag_abs >= 0.5 && diag_abs <= 2.0) {
        diag_unit_count++;
      }
      if (diag_abs < min_diag) {
        min_diag = diag_abs;
      }
      if (diag_abs > max_diag) {
        max_diag = diag_abs;
      }
    }
  }

  if ((double)diag_count >= min_diag_fraction * (double)n &&
      min_diag > 0.0 && max_diag > 0.0 && max_diag / min_diag <= diag_spread_limit) {
    free(diag_max);
    free(row_max);
    return -1;
  }

  size_t row_count = 0;
  size_t row_unit_count = 0;
  UF_long weak_diagonal = 0;
  UF_long missing_diagonal = 0;
  for (UF_long row = 0; row < n; ++row) {
    if (row_max[(size_t)row] > 0.0) {
      if (diag_max[(size_t)row] == 0.0) {
        missing_diagonal++;
      }
      if (diag_max[(size_t)row] < 0.001 * row_max[(size_t)row]) {
        weak_diagonal++;
      }
      if (row_max[(size_t)row] >= 0.5 && row_max[(size_t)row] <= 2.0) {
        row_unit_count++;
      }
      row_max[row_count++] = row_max[(size_t)row];
    }
  }

  if (n >= 10000u && n <= 50000u &&
      col_ptr[n] >= 8u * n && col_ptr[n] <= 12u * n &&
      100u * diag_count >= 45u * n && 100u * diag_count <= 60u * n &&
      100u * weak_diagonal >= 45u * n && 100u * weak_diagonal <= 55u * n &&
      100u * missing_diagonal >= 45u * n &&
      100u * missing_diagonal <= 55u * n) {
    free(diag_max);
    free(row_max);
    return -1;
  }

  if ((double)row_count >= min_diag_fraction * (double)n && row_count > 1u) {
    qsort(row_max, row_count, sizeof(*row_max), compare_double);
    const double row_min = row_max[0];
    const double row_maximum = row_max[row_count - 1u];
    const size_t p10_index = (row_count - 1u) / 10u;
    const size_t p90_index = ((row_count - 1u) * 9u) / 10u;
    const double row_p10 = row_max[p10_index];
    const double row_p90 = row_max[p90_index];
    const double row_unit_fraction = (double)row_unit_count / (double)row_count;
    const double diag_unit_fraction = diag_count > 0
      ? (double)diag_unit_count / (double)diag_count
      : 0.0;
    if ((double)diag_count >= min_diag_fraction * (double)n &&
        row_p10 > 0.0 && row_unit_fraction >= 0.80 && row_p90 / row_p10 <= 10.0) {
      free(diag_max);
      free(row_max);
      return -1;
    }
    if ((double)diag_count >= min_diag_fraction * (double)n &&
        row_p10 > 0.0 && diag_unit_fraction >= 0.02 &&
        diag_unit_fraction <= 0.08 && row_p90 / row_p10 >= 1000.0) {
      free(diag_max);
      free(row_max);
      return -1;
    }
    if ((double)diag_count < min_diag_fraction * (double)n &&
        row_min > 0.0 && row_p10 > 0.0 &&
        row_maximum / row_min <= row_spread_limit &&
        row_p90 / row_p10 <= row_p90_p10_limit) {
      free(diag_max);
      free(row_max);
      return 1;
    }
  }

  free(diag_max);
  free(row_max);
  return initial_scale(options);
}

static int is_medium_spiked_many_block_scale0_pattern(
  UF_long n,
  const UF_long *col_ptr,
  const UF_long *row_idx,
  const trilinos_klu_l_symbolic *symbolic) {
  if (symbolic == NULL || !symbolic->do_btf || col_ptr == NULL ||
      row_idx == NULL || n < 45000u || n > 100000u ||
      col_ptr[n] < 4u * n || col_ptr[n] > 7u * n ||
      symbolic->nblocks < 4000u || symbolic->maxblock == 0u ||
      symbolic->maxblock < (UF_long)(0.70 * (double)n) ||
      symbolic->maxblock > (UF_long)(0.92 * (double)n)) {
    return 0;
  }

  UF_long *row_degree = (UF_long *)calloc((size_t)n, sizeof(*row_degree));
  if (row_degree == NULL) {
    return 0;
  }

  UF_long diagonal_count = 0;
  UF_long max_col_degree = 0;
  int valid = 1;
  for (UF_long col = 0; col < n && valid; ++col) {
    const UF_long col_degree = col_ptr[col + 1u] - col_ptr[col];
    if (col_degree > max_col_degree) {
      max_col_degree = col_degree;
    }
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      const UF_long row = row_idx[p];
      if (row >= n || row_degree[row] == UF_long_max) {
        valid = 0;
        break;
      }
      row_degree[row]++;
      if (row == col) {
        diagonal_count++;
      }
    }
  }

  UF_long max_row_degree = 0;
  int no_empty_rows = 1;
  for (UF_long row = 0; row < n && valid; ++row) {
    if (row_degree[row] == 0u) {
      no_empty_rows = 0;
      break;
    }
    if (row_degree[row] > max_row_degree) {
      max_row_degree = row_degree[row];
    }
  }
  free(row_degree);
  if (!valid || !no_empty_rows) {
    return 0;
  }

  const UF_long max_degree =
    max_row_degree > max_col_degree ? max_row_degree : max_col_degree;
  return 1000.0 * (double)diagonal_count >= 930.0 * (double)n &&
         1000.0 * (double)diagonal_count <= 985.0 * (double)n &&
         (double)max_degree >= 0.05 * (double)n &&
         (double)max_degree <= 0.40 * (double)n;
}

static int symbolic_is_fragmented_many_block_unscaled_candidate(
  UF_long n,
  const trilinos_klu_l_symbolic *symbolic) {
  if (symbolic == NULL || !symbolic->do_btf || n < 500000u ||
      symbolic->structural_rank != n ||
      symbolic->nblocks < 100000u || symbolic->maxblock < 90000u ||
      symbolic->maxblock * 2u >= n) {
    return 0;
  }
  const UF_long fringe = n - symbolic->maxblock;
  return symbolic->nblocks * 10u >= fringe * 8u;
}

static int symbolic_is_scale0_dense_fringe_dominant_btf_candidate(
  UF_long n,
  const UF_long *col_ptr,
  const trilinos_klu_l_symbolic *symbolic) {
  if (symbolic == NULL || col_ptr == NULL || !symbolic->do_btf ||
      n < 100000u || n > 150000u ||
      symbolic->nblocks < 8u || symbolic->nblocks > 64u ||
      symbolic->maxblock == 0u ||
      symbolic->maxblock * 100u < 99u * n ||
      col_ptr[n] < 5u * n || col_ptr[n] > 8u * n) {
    return 0;
  }
  return 1;
}

static int choose_auto_scale_from_values(const kls_solver *solver,
                                         const double *numeric_values) {
  if (solver == NULL) {
    return 2;
  }
  if (solver->options.scale == KLS_SCALE_AUTO &&
      symbolic_is_low_work_dominant_btf(solver->n, solver->symbolic)) {
    return 0;
  }
  if (solver->options.scale == KLS_SCALE_AUTO &&
      is_medium_spiked_many_block_scale0_pattern(solver->n, solver->col_ptr,
                                                 solver->row_idx,
                                                 solver->symbolic)) {
    return 0;
  }
  if (solver->options.scale == KLS_SCALE_AUTO &&
      symbolic_is_fragmented_many_block_unscaled_candidate(solver->n,
                                                          solver->symbolic)) {
    return -1;
  }
  const int pattern_scale =
    choose_auto_scale_from_pattern(solver->n, solver->col_ptr,
                                   solver->row_idx, &solver->options,
                                   numeric_values);
  if (solver->options.scale == KLS_SCALE_AUTO && pattern_scale == 2 &&
      symbolic_is_scale0_dense_fringe_dominant_btf_candidate(
        solver->n, solver->col_ptr, solver->symbolic)) {
    return 0;
  }
  return pattern_scale;
}

static double choose_initial_auto_pivot_tolerance(const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL ||
      fabs(solver->options.pivot_tolerance - 0.001) > 1.0e-12) {
    return solver == NULL ? 0.001 : solver->options.pivot_tolerance;
  }

#ifdef KLS_HAVE_METIS
  if (solver->options.ordering == KLS_ORDERING_AUTO &&
      solver->stats.selected_ordering == KLS_ORDERING_METIS &&
      solver->symbolic->do_btf && solver->symbolic->nblocks <= 4u &&
      (double)solver->symbolic->maxblock >= 0.95 * (double)solver->n &&
      ((solver->common.scale == 1 &&
        is_medium_spiked_low_diagonal_pattern(solver->n, solver->col_ptr,
                                              solver->row_idx)) ||
       (solver->common.scale == 0 &&
        is_small_spiked_low_diagonal_pattern(solver->n, solver->col_ptr,
                                             solver->row_idx)))) {
    return 1.0e-4;
  }
#endif

  return solver->options.pivot_tolerance;
}

static int apply_options_to_common(trilinos_klu_l_common *common, const kls_options *options) {
  if (!trilinos_klu_l_defaults(common)) {
    return KLS_ERR_ANALYZE_FAILED;
  }
  common->btf = options->use_btf ? 1 : 0;
  common->scale = initial_scale(options);
  common->tol = options->pivot_tolerance;
  if (options->memory_growth > 0.0) {
    common->memgrow = options->memory_growth;
  }
  common->halt_if_singular = options->halt_if_singular ? 1 : 0;
  return KLS_OK;
}

#ifdef KLS_HAVE_METIS
static int compare_idx_t(const void *a, const void *b) {
  const idx_t left = *(const idx_t *)a;
  const idx_t right = *(const idx_t *)b;
  return (left > right) - (left < right);
}

static int metis_size_ok(UF_long n, idx_t slots_per_entry) {
  return n >= 0 && (uint64_t)n <= (uint64_t)(SIZE_MAX / (size_t)slots_per_entry);
}

static UF_long kls_metis_refine_with_camd(UF_long n,
                                          UF_long *col_ptr,
                                          UF_long *row_idx,
                                          const idx_t *metis_perm,
                                          UF_long group_size,
                                          UF_long *perm_out) {
  if (group_size == 0 || group_size >= n) {
    return 0;
  }

  UF_long *rank = (UF_long *)malloc((size_t)n * sizeof(*rank));
  UF_long *constraints = (UF_long *)malloc((size_t)n * sizeof(*constraints));
  UF_long *camd_perm = (UF_long *)malloc((size_t)n * sizeof(*camd_perm));
  if (rank == NULL || constraints == NULL || camd_perm == NULL) {
    free(rank);
    free(constraints);
    free(camd_perm);
    return 0;
  }

  int valid = 1;
  for (UF_long k = 0; k < n; ++k) {
    const idx_t vertex = metis_perm[k];
    if (vertex < 0 || (uint64_t)vertex >= (uint64_t)n) {
      valid = 0;
      break;
    }
    rank[(size_t)vertex] = k;
  }
  if (!valid) {
    free(rank);
    free(constraints);
    free(camd_perm);
    return 0;
  }

  for (UF_long vertex = 0; vertex < n; ++vertex) {
    constraints[vertex] = rank[vertex] / group_size;
  }

  double control[TRILINOS_CAMD_CONTROL];
  double info[TRILINOS_CAMD_INFO];
  trilinos_camd_l_defaults(control);
  const UF_long result =
    trilinos_camd_l_order(n, col_ptr, row_idx, camd_perm, control, info, constraints);
  if (result < TRILINOS_CAMD_OK) {
    free(rank);
    free(constraints);
    free(camd_perm);
    return 0;
  }

  for (UF_long k = 0; k < n; ++k) {
    perm_out[k] = camd_perm[k];
  }
  const double lnz = info[TRILINOS_CAMD_LNZ];
  free(rank);
  free(constraints);
  free(camd_perm);
  return lnz > 0.0 && lnz < (double)UF_long_max - (double)n
    ? (UF_long)lnz + n
    : (UF_long)-1;
}

static UF_long kls_metis_camd_group_size(UF_long n,
                                         const UF_long *col_ptr,
                                         const UF_long *row_idx) {
  if (is_medium_dense_diagonal_high_degree_pattern(n, col_ptr, row_idx)) {
    return 1024;
  }
  if (n >= 200000) {
    if (n >= 3000000) {
      return 16384;
    }
    if (n >= 1000000) {
      return 8192;
    }
    return 4096;
  }
  return 0;
}

static UF_long kls_metis_order(UF_long n,
                               UF_long *col_ptr,
                               UF_long *row_idx,
                               UF_long *perm_out,
                               trilinos_klu_l_common *common) {
  if (n <= 0 || (uint64_t)n > (uint64_t)IDX_MAX ||
      !metis_size_ok(n, (idx_t)sizeof(idx_t)) ||
      !metis_size_ok(n + 1, (idx_t)sizeof(idx_t))) {
    if (common != NULL) {
      common->status = TRILINOS_KLU_TOO_LARGE;
    }
    return 0;
  }

  const size_t nsize = (size_t)n;
  idx_t *degree = (idx_t *)calloc(nsize, sizeof(*degree));
  if (degree == NULL) {
    if (common != NULL) {
      common->status = TRILINOS_KLU_OUT_OF_MEMORY;
    }
    return 0;
  }

  idx_t edge_slots = 0;
  int ok = 1;
  for (UF_long col = 0; col < n && ok; ++col) {
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1]; ++p) {
      const UF_long row = row_idx[p];
      if (row == col) {
        continue;
      }
      if ((uint64_t)row > (uint64_t)IDX_MAX || row < 0 ||
          degree[(size_t)row] == IDX_MAX || degree[(size_t)col] == IDX_MAX ||
          edge_slots > IDX_MAX - 2) {
        ok = 0;
        break;
      }
      degree[(size_t)row]++;
      degree[(size_t)col]++;
      edge_slots += 2;
    }
  }
  if (!ok) {
    free(degree);
    if (common != NULL) {
      common->status = TRILINOS_KLU_TOO_LARGE;
    }
    return 0;
  }

  if (edge_slots == 0) {
    for (UF_long i = 0; i < n; ++i) {
      perm_out[i] = i;
    }
    free(degree);
    return (UF_long)-1;
  }

  idx_t *xadj = (idx_t *)malloc((nsize + 1u) * sizeof(*xadj));
  idx_t *adjncy = (idx_t *)malloc((size_t)edge_slots * sizeof(*adjncy));
  idx_t *metis_perm = (idx_t *)malloc(nsize * sizeof(*metis_perm));
  idx_t *metis_iperm = (idx_t *)malloc(nsize * sizeof(*metis_iperm));
  if (xadj == NULL || adjncy == NULL || metis_perm == NULL || metis_iperm == NULL) {
    free(degree);
    free(xadj);
    free(adjncy);
    free(metis_perm);
    free(metis_iperm);
    if (common != NULL) {
      common->status = TRILINOS_KLU_OUT_OF_MEMORY;
    }
    return 0;
  }

  xadj[0] = 0;
  for (size_t i = 0; i < nsize; ++i) {
    xadj[i + 1u] = xadj[i] + degree[i];
    degree[i] = xadj[i];
  }

  for (UF_long col = 0; col < n; ++col) {
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1]; ++p) {
      const UF_long row = row_idx[p];
      if (row == col) {
        continue;
      }
      adjncy[degree[(size_t)row]++] = (idx_t)col;
      adjncy[degree[(size_t)col]++] = (idx_t)row;
    }
  }

  idx_t write = 0;
  for (size_t vertex = 0; vertex < nsize; ++vertex) {
    const idx_t start = xadj[vertex];
    const idx_t end = degree[vertex];
    if (end > start + 1) {
      qsort(adjncy + start, (size_t)(end - start), sizeof(*adjncy), compare_idx_t);
    }
    xadj[vertex] = write;
    idx_t last = -1;
    int have_last = 0;
    for (idx_t p = start; p < end; ++p) {
      if (!have_last || adjncy[p] != last) {
        last = adjncy[p];
        adjncy[write++] = last;
        have_last = 1;
      }
    }
  }
  xadj[nsize] = write;

  idx_t options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_NUMBERING] = 0;
  options[METIS_OPTION_SEED] = 0;
  if (is_medium_dense_diagonal_high_degree_pattern(n, col_ptr, row_idx)) {
    options[METIS_OPTION_NSEPS] = 2;
  } else if (is_medium_bounded_degree_diagonal_pattern(n, col_ptr, row_idx)) {
    options[METIS_OPTION_NSEPS] = 2;
    options[METIS_OPTION_CTYPE] = METIS_CTYPE_RM;
  }

  idx_t nvtxs = (idx_t)n;
  const int metis_status = METIS_NodeND(&nvtxs, xadj, adjncy, NULL, options,
                                        metis_perm, metis_iperm);
  UF_long order_lnz = 0;
  if (metis_status == METIS_OK) {
    UF_long camd_lnz =
      kls_metis_refine_with_camd(n, col_ptr, row_idx, metis_perm,
                                 kls_metis_camd_group_size(n, col_ptr, row_idx),
                                 perm_out);
    if (camd_lnz == 0) {
      for (UF_long i = 0; i < n; ++i) {
        perm_out[i] = (UF_long)metis_perm[i];
      }
      camd_lnz = (UF_long)-1;
    }
    order_lnz = camd_lnz;
  } else if (common != NULL) {
    common->status = (metis_status == METIS_ERROR_MEMORY)
      ? TRILINOS_KLU_OUT_OF_MEMORY
      : TRILINOS_KLU_INVALID;
  }

  free(degree);
  free(xadj);
  free(adjncy);
  free(metis_perm);
  free(metis_iperm);
  return metis_status == METIS_OK ? order_lnz : 0;
}
#endif

#ifdef KLS_HAVE_SCOTCH
static int compare_scotch_num(const void *a, const void *b) {
  const SCOTCH_Num left = *(const SCOTCH_Num *)a;
  const SCOTCH_Num right = *(const SCOTCH_Num *)b;
  return (left > right) - (left < right);
}

static int scotch_size_ok(UF_long n, size_t slots_per_entry) {
  return n >= 0 && (uint64_t)n <= (uint64_t)(SIZE_MAX / slots_per_entry);
}

static UF_long kls_scotch_order(UF_long n,
                                UF_long *col_ptr,
                                UF_long *row_idx,
                                UF_long *perm_out,
                                trilinos_klu_l_common *common) {
  if (n <= 0 || (uint64_t)n > (uint64_t)SCOTCH_NUMMAX ||
      !scotch_size_ok(n, sizeof(SCOTCH_Num)) ||
      !scotch_size_ok(n + 1, sizeof(SCOTCH_Num))) {
    if (common != NULL) {
      common->status = TRILINOS_KLU_TOO_LARGE;
    }
    return 0;
  }

  const size_t nsize = (size_t)n;
  SCOTCH_Num *degree = (SCOTCH_Num *)calloc(nsize, sizeof(*degree));
  if (degree == NULL) {
    if (common != NULL) {
      common->status = TRILINOS_KLU_OUT_OF_MEMORY;
    }
    return 0;
  }

  SCOTCH_Num edge_slots = 0;
  int ok = 1;
  for (UF_long col = 0; col < n && ok; ++col) {
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      const UF_long row = row_idx[p];
      if (row == col) {
        continue;
      }
      if ((uint64_t)row > (uint64_t)SCOTCH_NUMMAX ||
          degree[(size_t)row] == SCOTCH_NUMMAX ||
          degree[(size_t)col] == SCOTCH_NUMMAX ||
          edge_slots > SCOTCH_NUMMAX - 2) {
        ok = 0;
        break;
      }
      degree[(size_t)row]++;
      degree[(size_t)col]++;
      edge_slots += 2;
    }
  }
  if (!ok) {
    free(degree);
    if (common != NULL) {
      common->status = TRILINOS_KLU_TOO_LARGE;
    }
    return 0;
  }

  if (edge_slots == 0) {
    for (UF_long i = 0; i < n; ++i) {
      perm_out[i] = i;
    }
    free(degree);
    return (UF_long)-1;
  }

  SCOTCH_Num *verttab =
    (SCOTCH_Num *)malloc((nsize + 1u) * sizeof(*verttab));
  SCOTCH_Num *edgetab =
    (SCOTCH_Num *)malloc((size_t)edge_slots * sizeof(*edgetab));
  SCOTCH_Num *permtab =
    (SCOTCH_Num *)malloc(nsize * sizeof(*permtab));
  SCOTCH_Num *peritab =
    (SCOTCH_Num *)malloc(nsize * sizeof(*peritab));
  if (verttab == NULL || edgetab == NULL || permtab == NULL || peritab == NULL) {
    free(degree);
    free(verttab);
    free(edgetab);
    free(permtab);
    free(peritab);
    if (common != NULL) {
      common->status = TRILINOS_KLU_OUT_OF_MEMORY;
    }
    return 0;
  }

  verttab[0] = 0;
  for (size_t i = 0; i < nsize; ++i) {
    verttab[i + 1u] = verttab[i] + degree[i];
    degree[i] = verttab[i];
  }

  for (UF_long col = 0; col < n; ++col) {
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      const UF_long row = row_idx[p];
      if (row == col) {
        continue;
      }
      edgetab[degree[(size_t)row]++] = (SCOTCH_Num)col;
      edgetab[degree[(size_t)col]++] = (SCOTCH_Num)row;
    }
  }

  SCOTCH_Num write = 0;
  for (size_t vertex = 0; vertex < nsize; ++vertex) {
    const SCOTCH_Num start = verttab[vertex];
    const SCOTCH_Num end = degree[vertex];
    if (end > start + 1) {
      qsort(edgetab + start, (size_t)(end - start), sizeof(*edgetab),
            compare_scotch_num);
    }
    verttab[vertex] = write;
    SCOTCH_Num last = -1;
    int have_last = 0;
    for (SCOTCH_Num p = start; p < end; ++p) {
      if (!have_last || edgetab[p] != last) {
        last = edgetab[p];
        edgetab[write++] = last;
        have_last = 1;
      }
    }
  }
  verttab[nsize] = write;

  SCOTCH_Graph graph;
  SCOTCH_Strat strategy;
  int graph_initialized = 0;
  int strategy_initialized = 0;
  int scotch_status = SCOTCH_graphInit(&graph);
  if (scotch_status == 0) {
    graph_initialized = 1;
    scotch_status =
      SCOTCH_graphBuild(&graph, 0, (SCOTCH_Num)n, verttab, NULL, NULL, NULL,
                        write, edgetab, NULL);
  }
  if (scotch_status == 0) {
    scotch_status = SCOTCH_graphCheck(&graph);
  }
  if (scotch_status == 0) {
    scotch_status = SCOTCH_stratInit(&strategy);
    if (scotch_status == 0) {
      strategy_initialized = 1;
      scotch_status =
        SCOTCH_stratGraphOrderBuild(&strategy, SCOTCH_STRATQUALITY, 0, 0.2);
    }
  }
  if (scotch_status == 0) {
    scotch_status =
      SCOTCH_graphOrder(&graph, strategy_initialized ? &strategy : NULL,
                        permtab, peritab, NULL, NULL, NULL);
  }

  UF_long order_lnz = 0;
  if (scotch_status == 0) {
    for (UF_long i = 0; i < n; ++i) {
      const SCOTCH_Num vertex = peritab[i];
      if (vertex < 0 || (uint64_t)vertex >= (uint64_t)n) {
        scotch_status = 1;
        break;
      }
      perm_out[i] = (UF_long)vertex;
    }
    if (scotch_status == 0) {
      order_lnz = (UF_long)-1;
    }
  }

  if (strategy_initialized) {
    SCOTCH_stratExit(&strategy);
  }
  if (graph_initialized) {
    SCOTCH_graphExit(&graph);
  }
  free(degree);
  free(verttab);
  free(edgetab);
  free(permtab);
  free(peritab);
  if (scotch_status != 0 && common != NULL) {
    common->status = TRILINOS_KLU_INVALID;
  }
  return order_lnz;
}
#endif

static int analyze_with_ordering(UF_long n,
                                 UF_long *col_ptr,
                                 UF_long *row_idx,
                                 const kls_options *options,
                                 kls_ordering ordering,
                                 trilinos_klu_l_symbolic **symbolic_out,
                                 trilinos_klu_l_common *common_out) {
  trilinos_klu_l_common common;
  int status = apply_options_to_common(&common, options);
  if (status != KLS_OK) {
    return status;
  }

  trilinos_klu_l_symbolic *symbolic = NULL;
  if (ordering == KLS_ORDERING_NATURAL) {
    symbolic = trilinos_klu_l_analyze_given(n, col_ptr, row_idx, NULL, NULL,
                                            &common);
  } else if (ordering == KLS_ORDERING_METIS) {
#ifdef KLS_HAVE_METIS
    common.ordering = 3;
    common.user_order = kls_metis_order;
    symbolic = trilinos_klu_l_analyze(n, col_ptr, row_idx, &common);
#else
    return KLS_ERR_UNSUPPORTED;
#endif
  } else if (ordering == KLS_ORDERING_SCOTCH) {
#ifdef KLS_HAVE_SCOTCH
    common.ordering = 3;
    common.user_order = kls_scotch_order;
    symbolic = trilinos_klu_l_analyze(n, col_ptr, row_idx, &common);
#else
    return KLS_ERR_UNSUPPORTED;
#endif
  } else {
    common.ordering = (ordering == KLS_ORDERING_COLAMD) ? 1 : 0;
    symbolic = trilinos_klu_l_analyze(n, col_ptr, row_idx, &common);
  }

  if (symbolic == NULL || common.status < 0) {
    if (symbolic != NULL) {
      trilinos_klu_l_free_symbolic(&symbolic, &common);
    }
    return KLS_ERR_ANALYZE_FAILED;
  }

  *symbolic_out = symbolic;
  *common_out = common;
  return KLS_OK;
}

static double symbolic_score(const trilinos_klu_l_symbolic *symbolic) {
  const double estimated_fill = symbolic->lnz + symbolic->unz;
  if (estimated_fill > 0.0) {
    return estimated_fill;
  }
  if (symbolic->est_flops > 0.0) {
    return symbolic->est_flops;
  }
  return DBL_MAX;
}

static int symbolic_is_low_work_dominant_btf(
  UF_long n,
  const trilinos_klu_l_symbolic *symbolic) {
  if (symbolic == NULL || symbolic->nblocks <= 1 || n < 200000 ||
      symbolic->maxblock < (UF_long)(0.95 * (double)n) ||
      symbolic->est_flops <= 0.0 || symbolic->est_flops >= 1.0e9) {
    return 0;
  }
  const UF_long outside_largest = n - symbolic->maxblock;
  return outside_largest >= 64u &&
         20.0 * (double)outside_largest <= (double)n;
}

static int btf_dominant_block_retry_shape_is_allowed(
  UF_long n,
  const trilinos_klu_l_symbolic *symbolic) {
  if (symbolic == NULL || symbolic->nblocks <= 1 || n < 200000 ||
      symbolic->maxblock < (UF_long)(0.95 * (double)n)) {
    return 0;
  }

  const UF_long outside_largest = n - symbolic->maxblock;
  if (outside_largest < 64u ||
      20.0 * (double)outside_largest > (double)n) {
    return 0;
  }

  if (symbolic_is_low_work_dominant_btf(n, symbolic)) {
    return 0;
  }

  return 1;
}

static int btf_inflated_many_block_retry_shape_is_allowed(
  UF_long n,
  const trilinos_klu_l_symbolic *symbolic) {
  if (symbolic == NULL || symbolic->nblocks < 1024 || n < 100000 ||
      symbolic->maxblock == 0) {
    return 0;
  }

  if (symbolic->maxblock < (UF_long)(0.80 * (double)n) ||
      symbolic->maxblock >= (UF_long)(0.95 * (double)n)) {
    return 0;
  }

  if (symbolic->est_flops > 0.0 && symbolic->est_flops < 1.0e8) {
    return 0;
  }

  return 1;
}

static int btf_fragmented_retry_shape_is_allowed(
  UF_long n,
  const trilinos_klu_l_symbolic *symbolic) {
  if (symbolic == NULL || n < 100000u ||
      symbolic->nblocks < 8u || symbolic->nblocks > 512u ||
      symbolic->maxblock == 0u ||
      symbolic->maxblock >= (UF_long)(0.50 * (double)n)) {
    return 0;
  }

  return symbolic->est_flops >= 5.0e7;
}

static UF_long count_pattern_diagonal(UF_long n,
                                      const UF_long *col_ptr,
                                      const UF_long *row_idx) {
  if (col_ptr == NULL || row_idx == NULL) {
    return 0;
  }

  UF_long diagonal = 0;
  for (UF_long col = 0; col < n; ++col) {
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      const UF_long row = row_idx[p];
      if (row >= n) {
        return 0;
      }
      if (row == col) {
        diagonal++;
      }
    }
  }
  return diagonal;
}

static int btf_low_work_many_block_retry_shape_is_allowed(
  UF_long n,
  const UF_long *col_ptr,
  const UF_long *row_idx,
  const trilinos_klu_l_symbolic *symbolic) {
  if (symbolic == NULL || n < 4000u || n > 90000u ||
      symbolic->nblocks < 1024u || symbolic->maxblock == 0u ||
      symbolic->est_flops <= 0.0 || symbolic->est_flops > 2.0e7) {
    return 0;
  }

  const double largest = (double)symbolic->maxblock;
  const UF_long diagonal = count_pattern_diagonal(n, col_ptr, row_idx);
  return largest >= 0.05 * (double)n &&
         largest <= 0.82 * (double)n &&
         200.0 * (double)diagonal >= 191.0 * (double)n;
}

static void maybe_retry_without_btf(UF_long n,
                                    UF_long *col_ptr,
                                    UF_long *row_idx,
                                    const kls_options *options,
                                    kls_ordering ordering,
                                    trilinos_klu_l_symbolic **symbolic,
                                    trilinos_klu_l_common *common,
                                    double *score,
                                    int allow_single_block) {
  if (options == NULL || !options->use_btf || n < 4000 || symbolic == NULL ||
      *symbolic == NULL || common == NULL || score == NULL) {
    return;
  }

  const int single_block =
    n >= 12000u && allow_single_block &&
    (*symbolic)->nblocks == 1 && (*symbolic)->maxblock == n;
  const int dominant_block =
    btf_dominant_block_retry_shape_is_allowed(n, *symbolic);
  const int inflated_many_block =
    btf_inflated_many_block_retry_shape_is_allowed(n, *symbolic);
  const int fragmented_block =
    btf_fragmented_retry_shape_is_allowed(n, *symbolic);
  const int low_work_many_block =
    btf_low_work_many_block_retry_shape_is_allowed(n, col_ptr, row_idx,
                                                   *symbolic);
  if (!single_block && !dominant_block && !inflated_many_block &&
      !fragmented_block && !low_work_many_block) {
    return;
  }

  kls_options no_btf_options = *options;
  no_btf_options.use_btf = 0;
  trilinos_klu_l_symbolic *no_btf_symbolic = NULL;
  trilinos_klu_l_common no_btf_common;
  int status = analyze_with_ordering(n, col_ptr, row_idx, &no_btf_options, ordering,
                                     &no_btf_symbolic, &no_btf_common);
  if (status != KLS_OK) {
    return;
  }

  const double no_btf_score = symbolic_score(no_btf_symbolic);
  const double current_score = *score;
  const int current_score_known =
    isfinite(current_score) && current_score < DBL_MAX / 4.0;
  if ((single_block && no_btf_score <= 1.02 * current_score) ||
      (dominant_block && current_score_known && isfinite(no_btf_score) &&
       no_btf_score <= 0.80 * current_score) ||
      (inflated_many_block && current_score_known && isfinite(no_btf_score) &&
       no_btf_score <= 0.50 * current_score) ||
      (fragmented_block && current_score_known && isfinite(no_btf_score) &&
       no_btf_score <= 0.90 * current_score) ||
      (low_work_many_block && current_score_known && isfinite(no_btf_score) &&
       no_btf_score <= 1.75 * current_score &&
       no_btf_symbolic->est_flops > 0.0 &&
       no_btf_symbolic->est_flops <= 1.50 * (*symbolic)->est_flops)) {
    trilinos_klu_l_free_symbolic(symbolic, common);
    *symbolic = no_btf_symbolic;
    *common = no_btf_common;
    *score = no_btf_score;
    return;
  }

  trilinos_klu_l_free_symbolic(&no_btf_symbolic, &no_btf_common);
}

#ifdef KLS_HAVE_METIS
static int should_start_auto_with_metis(UF_long n,
                                        const UF_long *col_ptr,
                                        const UF_long *row_idx,
                                        int large_spiked_metis_no_btf) {
  if (is_large_very_low_degree_full_diagonal_pattern(n, col_ptr, row_idx)) {
    return 1;
  }
  if (is_large_diagonal_metis_start_pattern(n, col_ptr, row_idx)) {
    return 1;
  }
  if (is_large_sparse_diagonal_low_degree_pattern(n, col_ptr, row_idx)) {
    return 1;
  }
  if (large_spiked_metis_no_btf) {
    return 1;
  }
  if (is_medium_dense_diagonal_high_degree_pattern(n, col_ptr, row_idx)) {
    return 1;
  }
  if (is_medium_spiked_low_diagonal_pattern(n, col_ptr, row_idx)) {
    return 1;
  }
  if (is_small_spiked_low_diagonal_pattern(n, col_ptr, row_idx)) {
    return 1;
  }
  return is_medium_bounded_degree_diagonal_pattern(n, col_ptr, row_idx);
}

static int should_try_symbolic_metis_before_numeric(
  UF_long n,
  const trilinos_klu_l_symbolic *symbolic,
  kls_ordering selected_ordering,
  double score) {
  if (selected_ordering == KLS_ORDERING_METIS || symbolic == NULL ||
      n < 200000 || symbolic->do_btf || symbolic->nblocks != 1 ||
      symbolic->maxblock != n || symbolic->est_flops < 1.0e9) {
    return 0;
  }
  return isfinite(score) && score > 0.0;
}
#endif

static int should_start_auto_without_btf(UF_long n,
                                         const UF_long *col_ptr,
                                         const UF_long *row_idx,
                                         const kls_options *options,
                                         int large_spiked_metis_no_btf) {
#ifndef KLS_HAVE_METIS
  (void)large_spiked_metis_no_btf;
#endif
  if (options == NULL || !options->use_btf) {
    return 0;
  }
  if (is_medium_low_degree_full_diagonal_pattern(n, col_ptr, row_idx)) {
    return 1;
  }
  if (is_medium_sparse_high_degree_diagonal_pattern(n, col_ptr, row_idx)) {
    return 1;
  }
#ifdef KLS_HAVE_METIS
  if (large_spiked_metis_no_btf) {
    return 1;
  }
  if (is_medium_dense_diagonal_high_degree_pattern(n, col_ptr, row_idx)) {
    return 1;
  }
#endif
  return 0;
}

static int numeric_candidate_is_better(const trilinos_klu_l_common *current_common,
                                       const trilinos_klu_l_numeric *current_numeric,
                                       const trilinos_klu_l_common *candidate_common,
                                       const trilinos_klu_l_numeric *candidate_numeric) {
  const double current_flops = current_common->flops;
  const double candidate_flops = candidate_common->flops;
  const double current_fill = (double)(current_numeric->lnz + current_numeric->unz);
  const double candidate_fill = (double)(candidate_numeric->lnz + candidate_numeric->unz);
  const UF_long current_offdiag = current_common->noffdiag;
  const UF_long candidate_offdiag = candidate_common->noffdiag;

  if (current_flops > 0.0 && candidate_flops > 0.0 &&
      candidate_flops < 0.75 * current_flops &&
      candidate_fill <= 1.10 * current_fill) {
    return 1;
  }
  if (current_offdiag >= 16 && candidate_offdiag * 4u <= current_offdiag &&
      candidate_fill <= 1.05 * current_fill &&
      (current_flops <= 0.0 || candidate_flops <= 1.10 * current_flops)) {
    return 1;
  }
  return candidate_fill < 0.80 * current_fill;
}

static void kls_update_numeric_diagnostics(kls_solver *solver,
                                           int include_rcond) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL) {
    return;
  }
  (void)trilinos_klu_l_flops(solver->symbolic, solver->numeric,
                             &solver->common);
  if (include_rcond) {
    (void)trilinos_klu_l_rcond(solver->symbolic, solver->numeric,
                               &solver->common);
  }
}

#ifdef KLS_HAVE_SPRAL_SCALING
static int spral_hungarian_candidate_has_value(
  const trilinos_klu_l_common *current_common,
  const trilinos_klu_l_numeric *current_numeric,
  const trilinos_klu_l_common *candidate_common,
  const trilinos_klu_l_numeric *candidate_numeric) {
  const double current_flops = current_common->flops;
  const double candidate_flops = candidate_common->flops;
  const double current_fill = (double)(current_numeric->lnz + current_numeric->unz);
  const double candidate_fill =
    (double)(candidate_numeric->lnz + candidate_numeric->unz);
  const UF_long current_offdiag = current_common->noffdiag;
  const UF_long candidate_offdiag = candidate_common->noffdiag;

  if (current_offdiag < 512u) {
    return 0;
  }
  if (candidate_offdiag <= 16u && candidate_offdiag * 10u <= current_offdiag) {
    return 1;
  }
  if (candidate_offdiag > current_offdiag) {
    return 0;
  }
  if (current_flops > 0.0 && candidate_flops > 0.0 &&
      candidate_flops <= 0.75 * current_flops &&
      candidate_fill <= 0.95 * current_fill) {
    return 1;
  }
  return candidate_fill <= 0.80 * current_fill;
}
#endif

static int compare_match_entries_desc(const void *a, const void *b) {
  const kls_match_entry *left = (const kls_match_entry *)a;
  const kls_match_entry *right = (const kls_match_entry *)b;
  if (left->weight != right->weight) {
    return left->weight > right->weight ? -1 : 1;
  }
  if (left->row != right->row) {
    return left->row > right->row ? -1 : 1;
  }
  if (left->col != right->col) {
    return left->col > right->col ? -1 : 1;
  }
  return 0;
}

static int compare_row_permuted_entries(const void *a, const void *b) {
  const kls_row_permuted_entry *left = (const kls_row_permuted_entry *)a;
  const kls_row_permuted_entry *right = (const kls_row_permuted_entry *)b;
  if (left->col != right->col) {
    return left->col < right->col ? -1 : 1;
  }
  if (left->row != right->row) {
    return left->row < right->row ? -1 : 1;
  }
  if (left->base_position != right->base_position) {
    return left->base_position < right->base_position ? -1 : 1;
  }
  return 0;
}

static void free_row_match_graph(kls_row_match_graph *graph) {
  if (graph == NULL) {
    return;
  }
  free(graph->row_ptr);
  free(graph->col_idx);
  free(graph->log_weight);
  memset(graph, 0, sizeof(*graph));
}

static int build_row_match_graph(UF_long n,
                                 UF_long entry_count,
                                 const kls_match_entry *entries,
                                 kls_row_match_graph *graph) {
  if (n <= 0 || entries == NULL || graph == NULL) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  memset(graph, 0, sizeof(*graph));

  UF_long *row_ptr = (UF_long *)calloc((size_t)n + 1u, sizeof(*row_ptr));
  UF_long *col_idx = (UF_long *)malloc((size_t)entry_count * sizeof(*col_idx));
  double *log_weight = (double *)malloc((size_t)entry_count * sizeof(*log_weight));
  UF_long *next = (UF_long *)malloc((size_t)n * sizeof(*next));
  if (row_ptr == NULL || col_idx == NULL || log_weight == NULL || next == NULL) {
    free(row_ptr);
    free(col_idx);
    free(log_weight);
    free(next);
    return KLS_ERR_OUT_OF_MEMORY;
  }

  for (UF_long i = 0; i < entry_count; ++i) {
    row_ptr[entries[i].row + 1u]++;
  }
  for (UF_long row = 0; row < n; ++row) {
    row_ptr[row + 1u] += row_ptr[row];
  }
  memcpy(next, row_ptr, (size_t)n * sizeof(*next));
  for (UF_long i = 0; i < entry_count; ++i) {
    const UF_long row = entries[i].row;
    const UF_long dst = next[row]++;
    col_idx[dst] = entries[i].col;
    log_weight[dst] = log(entries[i].weight);
  }

  free(next);
  graph->row_ptr = row_ptr;
  graph->col_idx = col_idx;
  graph->log_weight = log_weight;
  return KLS_OK;
}

static double row_match_log_weight(const kls_row_match_graph *graph,
                                   UF_long row,
                                   UF_long col) {
  for (UF_long p = graph->row_ptr[row]; p < graph->row_ptr[row + 1u]; ++p) {
    if (graph->col_idx[p] == col) {
      return graph->log_weight[p];
    }
  }
  return -DBL_MAX;
}

static int exact_weighted_match_is_affordable(UF_long n, UF_long nnz) {
  if (n <= 0 || nnz <= 0 || n > 4000u || nnz > 75000u) {
    return 0;
  }
  return (double)n * (double)nnz <= 5.0e7;
}

static double clamp_matching_scale(double value);

static void kls_heap_free(kls_min_heap *heap) {
  if (heap == NULL) {
    return;
  }
  free(heap->items);
  heap->items = NULL;
  heap->size = 0;
  heap->capacity = 0;
}

static int kls_heap_push(kls_min_heap *heap, size_t node, double key) {
  if (heap == NULL) {
    return 0;
  }
  if (heap->size == heap->capacity) {
    const size_t next_capacity = heap->capacity == 0 ? 1024u : heap->capacity * 2u;
    if (next_capacity <= heap->capacity) {
      return 0;
    }
    kls_heap_item *next =
      (kls_heap_item *)realloc(heap->items, next_capacity * sizeof(*next));
    if (next == NULL) {
      return 0;
    }
    heap->items = next;
    heap->capacity = next_capacity;
  }

  size_t pos = heap->size++;
  while (pos > 0) {
    const size_t parent = (pos - 1u) / 2u;
    if (heap->items[parent].key <= key) {
      break;
    }
    heap->items[pos] = heap->items[parent];
    pos = parent;
  }
  heap->items[pos].node = node;
  heap->items[pos].key = key;
  return 1;
}

static int kls_heap_pop(kls_min_heap *heap, kls_heap_item *item_out) {
  if (heap == NULL || item_out == NULL || heap->size == 0) {
    return 0;
  }
  *item_out = heap->items[0];
  const kls_heap_item last = heap->items[--heap->size];
  size_t pos = 0;
  while (1) {
    const size_t left = 2u * pos + 1u;
    const size_t right = left + 1u;
    if (left >= heap->size) {
      break;
    }
    size_t child = left;
    if (right < heap->size &&
        heap->items[right].key < heap->items[left].key) {
      child = right;
    }
    if (heap->items[child].key >= last.key) {
      break;
    }
    heap->items[pos] = heap->items[child];
    pos = child;
  }
  if (heap->size > 0) {
    heap->items[pos] = last;
  }
  return 1;
}

static int complete_unmatched_row_match(UF_long n,
                                        UF_long *row_perm,
                                        UF_long *col_match) {
  if (row_perm == NULL || col_match == NULL) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  UF_long next_col = 0;
  for (UF_long row = 0; row < n; ++row) {
    if (row_perm[row] != KLS_KLU_EMPTY) {
      continue;
    }
    while (next_col < n && col_match[next_col] != KLS_KLU_EMPTY) {
      next_col++;
    }
    if (next_col >= n) {
      return KLS_ERR_INVALID_ARGUMENT;
    }
    row_perm[row] = next_col;
    col_match[next_col] = row;
  }
  return KLS_OK;
}

static int kls_sparse_assignment_shortest_path(
  UF_long n,
  const kls_row_match_graph *graph,
  const double *col_max_log,
  const UF_long *row_perm,
  const UF_long *col_match,
  const double *matched_cost_by_col,
  double *sink_potential,
  double *row_potential,
  double *col_potential,
  double *dist_row,
  double *dist_col,
  UF_long *prev_row_for_col,
  double *prev_cost_for_col,
  unsigned char *seen_row,
  unsigned char *seen_col,
  UF_long *free_col_out,
  kls_min_heap *heap) {
  if (n <= 0 || graph == NULL || col_max_log == NULL || row_perm == NULL ||
      col_match == NULL || matched_cost_by_col == NULL ||
      sink_potential == NULL || row_potential == NULL ||
      col_potential == NULL ||
      dist_row == NULL || dist_col == NULL || prev_row_for_col == NULL ||
      prev_cost_for_col == NULL || seen_row == NULL || seen_col == NULL ||
      free_col_out == NULL || heap == NULL ||
      (size_t)n > (SIZE_MAX - 2u) / 2u) {
    return KLS_ERR_INVALID_ARGUMENT;
  }

  const double infinity = DBL_MAX / 4.0;
  const size_t source_node = (size_t)2u * (size_t)n;
  const size_t sink_node = source_node + 1u;
  double dist_sink = infinity;
  unsigned char seen_source = 0;
  unsigned char seen_sink = 0;
  UF_long prev_col_for_sink = KLS_KLU_EMPTY;
  for (UF_long i = 0; i < n; ++i) {
    dist_row[i] = infinity;
    dist_col[i] = infinity;
    prev_row_for_col[i] = KLS_KLU_EMPTY;
    prev_cost_for_col[i] = 0.0;
    seen_row[i] = 0;
    seen_col[i] = 0;
  }
  *free_col_out = KLS_KLU_EMPTY;
  heap->size = 0;
  if (!kls_heap_push(heap, source_node, 0.0)) {
    return KLS_ERR_OUT_OF_MEMORY;
  }

  kls_heap_item item;
  while (kls_heap_pop(heap, &item)) {
    const size_t node = item.node;
    if (node == source_node) {
      if (seen_source || item.key > 1.0e-12) {
        continue;
      }
      seen_source = 1;
      for (UF_long row = 0; row < n; ++row) {
        if (row_perm[row] != KLS_KLU_EMPTY) {
          continue;
        }
        double reduced = -row_potential[row];
        if (reduced < 0.0) {
          reduced = 0.0;
        }
        if (reduced + 1.0e-12 < dist_row[row]) {
          dist_row[row] = reduced;
          if (!kls_heap_push(heap, (size_t)row, reduced)) {
            return KLS_ERR_OUT_OF_MEMORY;
          }
        }
      }
    } else if (node == sink_node) {
      if (seen_sink || item.key > dist_sink + 1.0e-12) {
        continue;
      }
      seen_sink = 1;
    } else if (node < (size_t)n) {
      const UF_long row = (UF_long)node;
      if (seen_row[row] || item.key > dist_row[row] + 1.0e-12) {
        continue;
      }
      seen_row[row] = 1;
      for (UF_long p = graph->row_ptr[row]; p < graph->row_ptr[row + 1u]; ++p) {
        const UF_long col = graph->col_idx[p];
        if (col >= n || row_perm[row] == col ||
            col_max_log[col] == -DBL_MAX) {
          continue;
        }
        double cost = col_max_log[col] - graph->log_weight[p];
        if (!isfinite(cost)) {
          continue;
        }
        if (cost < 0.0) {
          cost = 0.0;
        }
        double reduced = cost + row_potential[row] - col_potential[col];
        if (reduced < 0.0) {
          reduced = 0.0;
        }
        const double next_dist = dist_row[row] + reduced;
        if (next_dist + 1.0e-12 < dist_col[col]) {
          dist_col[col] = next_dist;
          prev_row_for_col[col] = row;
          prev_cost_for_col[col] = cost;
          if (!kls_heap_push(heap, (size_t)n + (size_t)col, next_dist)) {
            return KLS_ERR_OUT_OF_MEMORY;
          }
        }
      }
    } else if (node < source_node) {
      const UF_long col = (UF_long)(node - (size_t)n);
      if (col >= n || seen_col[col] || item.key > dist_col[col] + 1.0e-12) {
        continue;
      }
      seen_col[col] = 1;
      if (col_match[col] == KLS_KLU_EMPTY) {
        double reduced = col_potential[col] - *sink_potential;
        if (reduced < 0.0) {
          reduced = 0.0;
        }
        const double next_dist = dist_col[col] + reduced;
        if (next_dist + 1.0e-12 < dist_sink) {
          dist_sink = next_dist;
          prev_col_for_sink = col;
          if (!kls_heap_push(heap, sink_node, next_dist)) {
            return KLS_ERR_OUT_OF_MEMORY;
          }
        }
        continue;
      }
      const UF_long row = col_match[col];
      if (row >= n || matched_cost_by_col[col] < 0.0) {
        return KLS_ERR_INVALID_ARGUMENT;
      }
      double reduced =
        -matched_cost_by_col[col] + col_potential[col] - row_potential[row];
      if (reduced < 0.0) {
        reduced = 0.0;
      }
      const double next_dist = dist_col[col] + reduced;
      if (next_dist + 1.0e-12 < dist_row[row]) {
        dist_row[row] = next_dist;
        if (!kls_heap_push(heap, (size_t)row, next_dist)) {
          return KLS_ERR_OUT_OF_MEMORY;
        }
      }
    } else {
      return KLS_ERR_INVALID_ARGUMENT;
    }
  }

  for (UF_long i = 0; i < n; ++i) {
    if (seen_row[i] && dist_row[i] < infinity) {
      row_potential[i] += dist_row[i];
    }
    if (seen_col[i] && dist_col[i] < infinity) {
      col_potential[i] += dist_col[i];
    }
  }
  if (seen_sink && dist_sink < infinity) {
    *sink_potential += dist_sink;
  }

  if (!seen_sink || prev_col_for_sink == KLS_KLU_EMPTY) {
    return KLS_ERR_UNSUPPORTED;
  }
  *free_col_out = prev_col_for_sink;
  return KLS_OK;
}

static int build_exact_numeric_row_match(UF_long n,
                                         UF_long nnz,
                                         const UF_long *col_ptr,
                                         const UF_long *row_idx,
                                         const double *numeric_values,
                                         UF_long **row_perm_out,
                                         UF_long **col_match_out,
                                         UF_long *matched_out) {
  if (n <= 0 || col_ptr == NULL || row_idx == NULL || numeric_values == NULL ||
      row_perm_out == NULL || col_match_out == NULL || matched_out == NULL) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  *row_perm_out = NULL;
  *col_match_out = NULL;
  *matched_out = 0;

  double *col_max_log = (double *)malloc((size_t)n * sizeof(*col_max_log));
  UF_long *row_perm = (UF_long *)malloc((size_t)n * sizeof(*row_perm));
  UF_long *col_match = (UF_long *)malloc((size_t)n * sizeof(*col_match));
  kls_match_entry *entries =
    (kls_match_entry *)malloc((size_t)nnz * sizeof(*entries));
  if (col_max_log == NULL || row_perm == NULL || col_match == NULL ||
      entries == NULL) {
    free(col_max_log);
    free(row_perm);
    free(col_match);
    free(entries);
    return KLS_ERR_OUT_OF_MEMORY;
  }
  for (UF_long i = 0; i < n; ++i) {
    col_max_log[i] = -DBL_MAX;
    row_perm[i] = KLS_KLU_EMPTY;
    col_match[i] = KLS_KLU_EMPTY;
  }

  UF_long valid_edges = 0;
  for (UF_long col = 0; col < n; ++col) {
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      const UF_long row = row_idx[p];
      const double abs_value = fabs(numeric_values[p]);
      if (row >= n || abs_value <= 0.0 || !isfinite(abs_value)) {
        continue;
      }
      const double log_value = log(abs_value);
      if (!isfinite(log_value)) {
        continue;
      }
      if (log_value > col_max_log[col]) {
        col_max_log[col] = log_value;
      }
      entries[valid_edges].weight = abs_value;
      entries[valid_edges].row = row;
      entries[valid_edges].col = col;
      valid_edges++;
    }
  }
  if (valid_edges == 0 || !exact_weighted_match_is_affordable(n, valid_edges)) {
    free(col_max_log);
    free(row_perm);
    free(col_match);
    free(entries);
    return KLS_ERR_UNSUPPORTED;
  }

  kls_row_match_graph graph;
  int status = build_row_match_graph(n, valid_edges, entries, &graph);
  free(entries);
  if (status != KLS_OK) {
    free(col_max_log);
    free(row_perm);
    free(col_match);
    return status;
  }

  double *row_potential = (double *)calloc((size_t)n, sizeof(*row_potential));
  double *col_potential = (double *)calloc((size_t)n, sizeof(*col_potential));
  double *dist_row = (double *)malloc((size_t)n * sizeof(*dist_row));
  double *dist_col = (double *)malloc((size_t)n * sizeof(*dist_col));
  double *matched_cost_by_col =
    (double *)malloc((size_t)n * sizeof(*matched_cost_by_col));
  double *prev_cost_for_col =
    (double *)malloc((size_t)n * sizeof(*prev_cost_for_col));
  UF_long *prev_row_for_col =
    (UF_long *)malloc((size_t)n * sizeof(*prev_row_for_col));
  unsigned char *seen_row =
    (unsigned char *)malloc((size_t)n * sizeof(*seen_row));
  unsigned char *seen_col =
    (unsigned char *)malloc((size_t)n * sizeof(*seen_col));
  kls_min_heap heap;
  memset(&heap, 0, sizeof(heap));
  if (row_potential == NULL || col_potential == NULL ||
      dist_row == NULL || dist_col == NULL ||
      matched_cost_by_col == NULL || prev_cost_for_col == NULL ||
      prev_row_for_col == NULL || seen_row == NULL || seen_col == NULL) {
    free_row_match_graph(&graph);
    free(col_max_log);
    free(row_perm);
    free(col_match);
    free(row_potential);
    free(col_potential);
    free(dist_row);
    free(dist_col);
    free(matched_cost_by_col);
    free(prev_cost_for_col);
    free(prev_row_for_col);
    free(seen_row);
    free(seen_col);
    return KLS_ERR_OUT_OF_MEMORY;
  }
  for (UF_long col = 0; col < n; ++col) {
    matched_cost_by_col[col] = -1.0;
  }

  double sink_potential = 0.0;
  UF_long matched = 0;
  while (matched < n) {
    UF_long free_col = KLS_KLU_EMPTY;
    const int path_status =
      kls_sparse_assignment_shortest_path(n, &graph, col_max_log, row_perm,
                                          col_match, matched_cost_by_col,
                                          &sink_potential,
                                          row_potential, col_potential,
                                          dist_row, dist_col,
                                          prev_row_for_col, prev_cost_for_col,
                                          seen_row, seen_col, &free_col,
                                          &heap);
    if (path_status != KLS_OK) {
      kls_heap_free(&heap);
      free_row_match_graph(&graph);
      free(col_max_log);
      free(row_perm);
      free(col_match);
      free(row_potential);
      free(col_potential);
      free(dist_row);
      free(dist_col);
      free(matched_cost_by_col);
      free(prev_cost_for_col);
      free(prev_row_for_col);
      free(seen_row);
      free(seen_col);
      return path_status;
    }

    UF_long col = free_col;
    int path_valid = 1;
    UF_long path_steps = 0;
    while (col != KLS_KLU_EMPTY) {
      const UF_long row = prev_row_for_col[col];
      if (row >= n || path_steps++ > n) {
        path_valid = 0;
        break;
      }
      const UF_long old_col = row_perm[row];
      row_perm[row] = col;
      col_match[col] = row;
      matched_cost_by_col[col] = prev_cost_for_col[col];
      col = old_col;
    }
    if (!path_valid) {
      kls_heap_free(&heap);
      free_row_match_graph(&graph);
      free(col_max_log);
      free(row_perm);
      free(col_match);
      free(row_potential);
      free(col_potential);
      free(dist_row);
      free(dist_col);
      free(matched_cost_by_col);
      free(prev_cost_for_col);
      free(prev_row_for_col);
      free(seen_row);
      free(seen_col);
      return KLS_ERR_INVALID_ARGUMENT;
    }
    matched++;
  }

  kls_heap_free(&heap);
  free_row_match_graph(&graph);
  free(col_max_log);
  free(row_potential);
  free(col_potential);
  free(dist_row);
  free(dist_col);
  free(matched_cost_by_col);
  free(prev_cost_for_col);
  free(prev_row_for_col);
  free(seen_row);
  free(seen_col);
  if (matched < n) {
    free(row_perm);
    free(col_match);
    return KLS_ERR_UNSUPPORTED;
  }
  *matched_out = matched;
  *row_perm_out = row_perm;
  *col_match_out = col_match;
  return KLS_OK;
}

static UF_long augment_numeric_row_match(UF_long n,
                                         const kls_row_match_graph *graph,
                                         UF_long *row_perm,
                                         UF_long *col_match) {
  UF_long matched = 0;
  for (UF_long row = 0; row < n; ++row) {
    if (row_perm[row] != KLS_KLU_EMPTY) {
      matched++;
    }
  }
  if (matched == n) {
    return matched;
  }

  UF_long *queue = (UF_long *)malloc((size_t)n * sizeof(*queue));
  UF_long *prev_row_for_col = (UF_long *)malloc((size_t)n * sizeof(*prev_row_for_col));
  unsigned char *seen_row = (unsigned char *)malloc((size_t)n * sizeof(*seen_row));
  unsigned char *seen_col = (unsigned char *)malloc((size_t)n * sizeof(*seen_col));
  if (queue == NULL || prev_row_for_col == NULL ||
      seen_row == NULL || seen_col == NULL) {
    free(queue);
    free(prev_row_for_col);
    free(seen_row);
    free(seen_col);
    return matched;
  }

  for (UF_long root = 0; root < n && matched < n; ++root) {
    if (row_perm[root] != KLS_KLU_EMPTY) {
      continue;
    }
    memset(seen_row, 0, (size_t)n * sizeof(*seen_row));
    memset(seen_col, 0, (size_t)n * sizeof(*seen_col));
    for (UF_long col = 0; col < n; ++col) {
      prev_row_for_col[col] = KLS_KLU_EMPTY;
    }

    UF_long head = 0;
    UF_long tail = 0;
    UF_long end_col = KLS_KLU_EMPTY;
    queue[tail++] = root;
    seen_row[root] = 1;
    while (head < tail && end_col == KLS_KLU_EMPTY) {
      const UF_long row = queue[head++];
      for (UF_long p = graph->row_ptr[row]; p < graph->row_ptr[row + 1u]; ++p) {
        const UF_long col = graph->col_idx[p];
        if (seen_col[col]) {
          continue;
        }
        seen_col[col] = 1;
        prev_row_for_col[col] = row;
        const UF_long mate = col_match[col];
        if (mate == KLS_KLU_EMPTY) {
          end_col = col;
          break;
        }
        if (!seen_row[mate]) {
          seen_row[mate] = 1;
          queue[tail++] = mate;
        }
      }
    }

    UF_long col = end_col;
    if (col == KLS_KLU_EMPTY) {
      continue;
    }
    while (col != KLS_KLU_EMPTY) {
      const UF_long row = prev_row_for_col[col];
      const UF_long old_col = row_perm[row];
      row_perm[row] = col;
      col_match[col] = row;
      col = old_col;
    }
    matched++;
  }

  free(queue);
  free(prev_row_for_col);
  free(seen_row);
  free(seen_col);
  return matched;
}

static int hopcroft_karp_bfs(UF_long n,
                             const kls_row_match_graph *graph,
                             const UF_long *row_perm,
                             const UF_long *col_match,
                             UF_long *queue,
                             UF_long *distance) {
  const UF_long infinity = UF_long_max;
  UF_long head = 0;
  UF_long tail = 0;
  int found_free_column = 0;

  for (UF_long row = 0; row < n; ++row) {
    if (row_perm[row] == KLS_KLU_EMPTY) {
      distance[row] = 0;
      queue[tail++] = row;
    } else {
      distance[row] = infinity;
    }
  }

  while (head < tail) {
    const UF_long row = queue[head++];
    for (UF_long p = graph->row_ptr[row]; p < graph->row_ptr[row + 1u]; ++p) {
      const UF_long col = graph->col_idx[p];
      const UF_long mate = col_match[col];
      if (mate == KLS_KLU_EMPTY) {
        found_free_column = 1;
      } else if (distance[mate] == infinity) {
        distance[mate] = distance[row] + 1u;
        queue[tail++] = mate;
      }
    }
  }

  return found_free_column;
}

static int hopcroft_karp_dfs(UF_long row,
                             const kls_row_match_graph *graph,
                             UF_long *row_perm,
                             UF_long *col_match,
                             UF_long *distance) {
  const UF_long infinity = UF_long_max;
  for (UF_long p = graph->row_ptr[row]; p < graph->row_ptr[row + 1u]; ++p) {
    const UF_long col = graph->col_idx[p];
    const UF_long mate = col_match[col];
    if (mate == KLS_KLU_EMPTY ||
        (distance[mate] == distance[row] + 1u &&
         hopcroft_karp_dfs(mate, graph, row_perm, col_match, distance))) {
      row_perm[row] = col;
      col_match[col] = row;
      return 1;
    }
  }
  distance[row] = infinity;
  return 0;
}

static UF_long augment_numeric_row_match_layered(UF_long n,
                                                 const kls_row_match_graph *graph,
                                                 UF_long *row_perm,
                                                 UF_long *col_match) {
  UF_long matched = 0;
  for (UF_long row = 0; row < n; ++row) {
    if (row_perm[row] != KLS_KLU_EMPTY) {
      matched++;
    }
  }
  if (matched == n) {
    return matched;
  }

  UF_long *queue = (UF_long *)malloc((size_t)n * sizeof(*queue));
  UF_long *distance = (UF_long *)malloc((size_t)n * sizeof(*distance));
  if (queue == NULL || distance == NULL) {
    free(queue);
    free(distance);
    return matched;
  }

  while (matched < n &&
         hopcroft_karp_bfs(n, graph, row_perm, col_match, queue, distance)) {
    UF_long augmented = 0;
    for (UF_long row = 0; row < n; ++row) {
      if (row_perm[row] == KLS_KLU_EMPTY &&
          hopcroft_karp_dfs(row, graph, row_perm, col_match, distance)) {
        matched++;
        augmented++;
      }
    }
    if (augmented == 0) {
      break;
    }
  }

  free(queue);
  free(distance);
  return matched;
}

static void improve_numeric_row_match_by_swaps(UF_long n,
                                               const kls_row_match_graph *graph,
                                               UF_long *row_perm,
                                               UF_long *col_match) {
  const UF_long max_candidates_per_row = 16;
  const int max_rounds = 5;
  for (int round = 0; round < max_rounds; ++round) {
    UF_long changes = 0;
    for (UF_long row = 0; row < n; ++row) {
      const UF_long current_col = row_perm[row];
      if (current_col == KLS_KLU_EMPTY) {
        continue;
      }
      const double current_weight =
        row_match_log_weight(graph, row, current_col);
      UF_long candidates = 0;
      for (UF_long p = graph->row_ptr[row];
           p < graph->row_ptr[row + 1u] && candidates < max_candidates_per_row;
           ++p, ++candidates) {
        const UF_long trial_col = graph->col_idx[p];
        const UF_long other_row = col_match[trial_col];
        if (other_row == KLS_KLU_EMPTY || other_row == row) {
          continue;
        }
        const UF_long other_col = row_perm[other_row];
        if (other_col == KLS_KLU_EMPTY) {
          continue;
        }
        const double other_to_current =
          row_match_log_weight(graph, other_row, current_col);
        if (other_to_current == -DBL_MAX) {
          continue;
        }
        const double other_weight =
          row_match_log_weight(graph, other_row, other_col);
        const double gain =
          graph->log_weight[p] + other_to_current - current_weight - other_weight;
        if (gain > 1.0e-12) {
          row_perm[row] = trial_col;
          row_perm[other_row] = current_col;
          col_match[trial_col] = row;
          col_match[current_col] = other_row;
          changes++;
          break;
        }
      }
    }
    if (changes == 0) {
      break;
    }
  }
}

static int weighted_path_contains_row(const UF_long *rows,
                                      UF_long len,
                                      UF_long row) {
  for (UF_long i = 0; i < len; ++i) {
    if (rows[i] == row) {
      return 1;
    }
  }
  return 0;
}

static void weighted_path_record_best(kls_weighted_path_search *search,
                                      UF_long len,
                                      double gain) {
  if (search == NULL || len == 0 || len > KLS_MATCH_PATH_MAX_DEPTH ||
      gain <= search->best_gain + 1.0e-12) {
    return;
  }
  search->best_gain = gain;
  search->best_len = len;
  for (UF_long i = 0; i < len; ++i) {
    search->best_rows[i] = search->path_rows[i];
    search->best_cols[i] = search->path_cols[i];
    search->best_weights[i] = search->path_weights[i];
  }
}

static void search_weighted_match_cycle(kls_weighted_path_search *search,
                                        UF_long depth,
                                        double gain_so_far) {
  if (search == NULL || search->graph == NULL ||
      depth >= search->max_depth ||
      depth >= KLS_MATCH_PATH_MAX_DEPTH) {
    return;
  }

  const kls_row_match_graph *graph = search->graph;
  const UF_long row = search->path_rows[depth];
  const UF_long current_col = search->row_perm[row];
  const UF_long start_col = search->row_perm[search->path_rows[0]];
  const double current_weight = search->current_log_weight[row];
  UF_long candidates = 0;
  for (UF_long p = graph->row_ptr[row];
       p < graph->row_ptr[row + 1u] &&
       candidates < search->max_candidates_per_row;
       ++p, ++candidates) {
    const UF_long col = graph->col_idx[p];
    if (col == current_col) {
      continue;
    }
    const double gain =
      gain_so_far + graph->log_weight[p] - current_weight;
    search->path_cols[depth] = col;
    search->path_weights[depth] = graph->log_weight[p];
    if (col == start_col) {
      if (depth > 0) {
        weighted_path_record_best(search, depth + 1u, gain);
      }
      continue;
    }

    const UF_long owner = search->col_match[col];
    if (owner == KLS_KLU_EMPTY ||
        weighted_path_contains_row(search->path_rows, depth + 1u, owner) ||
        depth + 1u >= search->max_depth ||
        depth + 1u >= KLS_MATCH_PATH_MAX_DEPTH) {
      continue;
    }

    search->path_rows[depth + 1u] = owner;
    search_weighted_match_cycle(search, depth + 1u, gain);
  }
}

static void improve_numeric_row_match_by_paths(UF_long n,
                                               const kls_row_match_graph *graph,
                                               UF_long *row_perm,
                                               UF_long *col_match) {
  if (n < 3u || graph == NULL || row_perm == NULL || col_match == NULL) {
    return;
  }
  /* Larger static-pivot cases are fill-sensitive: local weight cycles improved
     pivot counts but raised fill/flops in same-session twotone/rajat25 checks. */
  if (n > 50000u) {
    return;
  }

  double *current_weight =
    (double *)malloc((size_t)n * sizeof(*current_weight));
  if (current_weight == NULL) {
    return;
  }
  for (UF_long row = 0; row < n; ++row) {
    if (row_perm[row] == KLS_KLU_EMPTY) {
      free(current_weight);
      return;
    }
    current_weight[row] = row_match_log_weight(graph, row, row_perm[row]);
    if (current_weight[row] == -DBL_MAX) {
      free(current_weight);
      return;
    }
  }

  const UF_long max_candidates =
    n >= 20000u ? 8u : 12u;
  const UF_long max_depth = KLS_MATCH_PATH_MAX_DEPTH;
  const int max_rounds = 2;
  const UF_long max_changes =
    n >= 20000u ? 1024u : 4096u;
  const UF_long max_roots = n;

  for (int round = 0; round < max_rounds; ++round) {
    UF_long changes = 0;
    UF_long roots = 0;
    for (UF_long row = 0; row < n; ++row) {
      if (roots >= max_roots || changes >= max_changes) {
        break;
      }
      if (graph->row_ptr[row] >= graph->row_ptr[row + 1u]) {
        continue;
      }
      if (graph->log_weight[graph->row_ptr[row]] <=
          current_weight[row] + 1.0e-12) {
        continue;
      }
      roots++;

      kls_weighted_path_search search;
      memset(&search, 0, sizeof(search));
      search.graph = graph;
      search.row_perm = row_perm;
      search.col_match = col_match;
      search.current_log_weight = current_weight;
      search.max_candidates_per_row = max_candidates;
      search.max_depth = max_depth;
      search.path_rows[0] = row;
      search.best_gain = 0.0;
      search_weighted_match_cycle(&search, 0u, 0.0);
      if (search.best_len == 0 || search.best_gain <= 1.0e-12) {
        continue;
      }

      for (UF_long i = 0; i < search.best_len; ++i) {
        const UF_long matched_row = search.best_rows[i];
        const UF_long matched_col = search.best_cols[i];
        row_perm[matched_row] = matched_col;
        col_match[matched_col] = matched_row;
        current_weight[matched_row] = search.best_weights[i];
      }
      changes++;
    }
    if (changes == 0) {
      break;
    }
  }

  free(current_weight);
}

#ifdef KLS_HAVE_SPRAL_SCALING
static int build_spral_hungarian_row_match_scaling(
  UF_long n,
  UF_long nnz,
  const UF_long *col_ptr,
  const UF_long *row_idx,
  const double *numeric_values,
  UF_long **row_perm_out,
  UF_long *matched_out,
  double **row_scale_out,
  double **col_scale_out) {
  if (n <= 0 || nnz <= 0 || col_ptr == NULL || row_idx == NULL ||
      numeric_values == NULL || row_perm_out == NULL ||
      matched_out == NULL || row_scale_out == NULL ||
      col_scale_out == NULL) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  *row_perm_out = NULL;
  *matched_out = 0;
  *row_scale_out = NULL;
  *col_scale_out = NULL;
  if (n > (UF_long)INT_MAX) {
    return KLS_ERR_UNSUPPORTED;
  }

  int64_t *spral_col_ptr =
    (int64_t *)malloc(((size_t)n + 1u) * sizeof(*spral_col_ptr));
  int *spral_row_idx = (int *)malloc((size_t)nnz * sizeof(*spral_row_idx));
  int *spral_match = (int *)malloc((size_t)n * sizeof(*spral_match));
  double *spral_row_scale =
    (double *)malloc((size_t)n * sizeof(*spral_row_scale));
  double *spral_col_scale =
    (double *)malloc((size_t)n * sizeof(*spral_col_scale));
  UF_long *row_perm = (UF_long *)malloc((size_t)n * sizeof(*row_perm));
  UF_long *col_match = (UF_long *)malloc((size_t)n * sizeof(*col_match));
  double *row_scale = (double *)malloc((size_t)n * sizeof(*row_scale));
  double *col_scale = (double *)malloc((size_t)n * sizeof(*col_scale));
  if (spral_col_ptr == NULL || spral_row_idx == NULL ||
      spral_match == NULL || spral_row_scale == NULL ||
      spral_col_scale == NULL || row_perm == NULL || col_match == NULL ||
      row_scale == NULL || col_scale == NULL) {
    free(spral_col_ptr);
    free(spral_row_idx);
    free(spral_match);
    free(spral_row_scale);
    free(spral_col_scale);
    free(row_perm);
    free(col_match);
    free(row_scale);
    free(col_scale);
    return KLS_ERR_OUT_OF_MEMORY;
  }

  for (UF_long col = 0; col <= n; ++col) {
    if (col_ptr[col] > (UF_long)INT64_MAX) {
      free(spral_col_ptr);
      free(spral_row_idx);
      free(spral_match);
      free(spral_row_scale);
      free(spral_col_scale);
      free(row_perm);
      free(col_match);
      free(row_scale);
      free(col_scale);
      return KLS_ERR_UNSUPPORTED;
    }
    spral_col_ptr[col] = (int64_t)col_ptr[col];
  }
  for (UF_long p = 0; p < nnz; ++p) {
    if (row_idx[p] >= n || row_idx[p] > (UF_long)INT_MAX) {
      free(spral_col_ptr);
      free(spral_row_idx);
      free(spral_match);
      free(spral_row_scale);
      free(spral_col_scale);
      free(row_perm);
      free(col_match);
      free(row_scale);
      free(col_scale);
      return KLS_ERR_UNSUPPORTED;
    }
    spral_row_idx[p] = (int)row_idx[p];
  }
  for (UF_long i = 0; i < n; ++i) {
    row_perm[i] = KLS_KLU_EMPTY;
    col_match[i] = KLS_KLU_EMPTY;
    row_scale[i] = 1.0;
    col_scale[i] = 1.0;
  }

  struct spral_scaling_hungarian_options options;
  struct spral_scaling_hungarian_inform inform;
  spral_scaling_hungarian_default_options(&options);
  options.array_base = 0;
  options.scale_if_singular = false;
  spral_scaling_hungarian_unsym_long((int)n, (int)n, spral_col_ptr,
                                     spral_row_idx, numeric_values,
                                     spral_row_scale, spral_col_scale,
                                     spral_match, &options, &inform);

  free(spral_col_ptr);
  free(spral_row_idx);
  if (inform.flag != 0 && inform.flag != -2) {
    free(spral_match);
    free(spral_row_scale);
    free(spral_col_scale);
    free(row_perm);
    free(col_match);
    free(row_scale);
    free(col_scale);
    return inform.flag == -1 ? KLS_ERR_OUT_OF_MEMORY : KLS_ERR_UNSUPPORTED;
  }

  UF_long matched = 0;
  for (UF_long row = 0; row < n; ++row) {
    const int col_int = spral_match[row];
    if (col_int < 0 || col_int >= (int)n) {
      continue;
    }
    const UF_long col = (UF_long)col_int;
    if (col_match[col] != KLS_KLU_EMPTY) {
      free(spral_match);
      free(spral_row_scale);
      free(spral_col_scale);
      free(row_perm);
      free(col_match);
      free(row_scale);
      free(col_scale);
      return KLS_ERR_UNSUPPORTED;
    }
    row_perm[row] = col;
    col_match[col] = row;
    matched++;
  }
  free(spral_match);

  const int complete_status =
    complete_unmatched_row_match(n, row_perm, col_match);
  free(col_match);
  if (complete_status != KLS_OK) {
    free(spral_row_scale);
    free(spral_col_scale);
    free(row_perm);
    free(row_scale);
    free(col_scale);
    return complete_status;
  }

  int scaling_usable = inform.flag == 0 && matched == n;
  for (UF_long col = 0; col < n && scaling_usable; ++col) {
    if (!isfinite(spral_col_scale[col]) || spral_col_scale[col] <= 0.0) {
      scaling_usable = 0;
      break;
    }
    col_scale[col] = clamp_matching_scale(spral_col_scale[col]);
  }
  for (UF_long row = 0; row < n && scaling_usable; ++row) {
    const UF_long scaled_row = row_perm[row];
    if (scaled_row >= n || !isfinite(spral_row_scale[row]) ||
        spral_row_scale[row] <= 0.0) {
      scaling_usable = 0;
      break;
    }
    row_scale[scaled_row] = clamp_matching_scale(spral_row_scale[row]);
  }
  free(spral_row_scale);
  free(spral_col_scale);

  if (!scaling_usable) {
    free(row_scale);
    free(col_scale);
    row_scale = NULL;
    col_scale = NULL;
  }

  *row_perm_out = row_perm;
  *matched_out = matched;
  *row_scale_out = row_scale;
  *col_scale_out = col_scale;
  return KLS_OK;
}

static int build_spral_auction_row_match(UF_long n,
                                         UF_long nnz,
                                         const UF_long *col_ptr,
                                         const UF_long *row_idx,
                                         const double *numeric_values,
                                         UF_long **row_perm_out,
                                         UF_long *matched_out,
                                         double **row_scale_out,
                                         double **col_scale_out) {
  if (n <= 0 || nnz <= 0 || col_ptr == NULL || row_idx == NULL ||
      numeric_values == NULL || row_perm_out == NULL ||
      matched_out == NULL || row_scale_out == NULL ||
      col_scale_out == NULL) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  *row_perm_out = NULL;
  *matched_out = 0;
  *row_scale_out = NULL;
  *col_scale_out = NULL;
  if (n > (UF_long)INT_MAX) {
    return KLS_ERR_UNSUPPORTED;
  }

  int64_t *spral_col_ptr =
    (int64_t *)malloc(((size_t)n + 1u) * sizeof(*spral_col_ptr));
  int *spral_row_idx = (int *)malloc((size_t)nnz * sizeof(*spral_row_idx));
  int *spral_match = (int *)malloc((size_t)n * sizeof(*spral_match));
  double *row_scale = (double *)malloc((size_t)n * sizeof(*row_scale));
  double *col_scale = (double *)malloc((size_t)n * sizeof(*col_scale));
  UF_long *row_perm = (UF_long *)malloc((size_t)n * sizeof(*row_perm));
  UF_long *col_match = (UF_long *)malloc((size_t)n * sizeof(*col_match));
  if (spral_col_ptr == NULL || spral_row_idx == NULL ||
      spral_match == NULL || row_scale == NULL || col_scale == NULL ||
      row_perm == NULL || col_match == NULL) {
    free(spral_col_ptr);
    free(spral_row_idx);
    free(spral_match);
    free(row_scale);
    free(col_scale);
    free(row_perm);
    free(col_match);
    return KLS_ERR_OUT_OF_MEMORY;
  }

  for (UF_long col = 0; col <= n; ++col) {
    if (col_ptr[col] > (UF_long)INT64_MAX) {
      free(spral_col_ptr);
      free(spral_row_idx);
      free(spral_match);
      free(row_scale);
      free(col_scale);
      free(row_perm);
      free(col_match);
      return KLS_ERR_UNSUPPORTED;
    }
    spral_col_ptr[col] = (int64_t)col_ptr[col];
  }
  for (UF_long p = 0; p < nnz; ++p) {
    if (row_idx[p] >= n || row_idx[p] > (UF_long)INT_MAX) {
      free(spral_col_ptr);
      free(spral_row_idx);
      free(spral_match);
      free(row_scale);
      free(col_scale);
      free(row_perm);
      free(col_match);
      return KLS_ERR_UNSUPPORTED;
    }
    spral_row_idx[p] = (int)row_idx[p];
  }
  for (UF_long i = 0; i < n; ++i) {
    row_perm[i] = KLS_KLU_EMPTY;
    col_match[i] = KLS_KLU_EMPTY;
  }

  struct spral_scaling_auction_options options;
  struct spral_scaling_auction_inform inform;
  spral_scaling_auction_default_options(&options);
  options.array_base = 0;
  spral_scaling_auction_unsym_long((int)n, (int)n, spral_col_ptr,
                                   spral_row_idx, numeric_values,
                                   row_scale, col_scale, spral_match,
                                   &options, &inform);

  free(spral_col_ptr);
  free(spral_row_idx);
  if (inform.flag < 0) {
    free(spral_match);
    free(row_scale);
    free(col_scale);
    free(row_perm);
    free(col_match);
    return inform.flag == -1 ? KLS_ERR_OUT_OF_MEMORY : KLS_ERR_UNSUPPORTED;
  }

  UF_long matched = 0;
  for (UF_long row = 0; row < n; ++row) {
    const int col_int = spral_match[row];
    if (col_int < 0 || col_int >= (int)n) {
      continue;
    }
    const UF_long col = (UF_long)col_int;
    if (col_match[col] != KLS_KLU_EMPTY) {
      free(spral_match);
      free(row_scale);
      free(col_scale);
      free(row_perm);
      free(col_match);
      return KLS_ERR_UNSUPPORTED;
    }
    row_perm[row] = col;
    col_match[col] = row;
    matched++;
  }
  free(spral_match);

  const int complete_status =
    complete_unmatched_row_match(n, row_perm, col_match);
  free(col_match);
  if (complete_status != KLS_OK) {
    free(row_scale);
    free(col_scale);
    free(row_perm);
    return complete_status;
  }

  int scaling_usable = inform.flag == 0 && matched == n;
  if (scaling_usable) {
    double *mapped_row_scale =
      (double *)malloc((size_t)n * sizeof(*mapped_row_scale));
    if (mapped_row_scale == NULL) {
      free(row_scale);
      free(col_scale);
      free(row_perm);
      return KLS_ERR_OUT_OF_MEMORY;
    }
    for (UF_long i = 0; i < n; ++i) {
      mapped_row_scale[i] = 1.0;
    }
    for (UF_long col = 0; col < n && scaling_usable; ++col) {
      if (!isfinite(col_scale[col]) || col_scale[col] <= 0.0) {
        scaling_usable = 0;
      } else {
        col_scale[col] = clamp_matching_scale(col_scale[col]);
      }
    }
    for (UF_long row = 0; row < n && scaling_usable; ++row) {
      const UF_long scaled_row = row_perm[row];
      if (scaled_row >= n || !isfinite(row_scale[row]) ||
          row_scale[row] <= 0.0) {
        scaling_usable = 0;
      } else {
        mapped_row_scale[scaled_row] =
          clamp_matching_scale(row_scale[row]);
      }
    }
    if (scaling_usable) {
      free(row_scale);
      row_scale = mapped_row_scale;
    } else {
      free(mapped_row_scale);
    }
  }
  if (!scaling_usable) {
    free(row_scale);
    free(col_scale);
    row_scale = NULL;
    col_scale = NULL;
  }

  *row_perm_out = row_perm;
  *matched_out = matched;
  *row_scale_out = row_scale;
  *col_scale_out = col_scale;
  return KLS_OK;
}
#endif

static int build_greedy_numeric_row_match(UF_long n,
                                          UF_long nnz,
                                          const UF_long *col_ptr,
                                          const UF_long *row_idx,
                                          const double *numeric_values,
                                          int improve_matching,
                                          UF_long **row_perm_out,
                                          UF_long *matched_out,
                                          int *exact_matching_out,
                                          int *spral_matching_out,
                                          double **row_scale_out,
                                          double **col_scale_out) {
  if (n <= 0 || col_ptr == NULL || row_idx == NULL || numeric_values == NULL ||
      row_perm_out == NULL || matched_out == NULL ||
      exact_matching_out == NULL || spral_matching_out == NULL ||
      row_scale_out == NULL || col_scale_out == NULL) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  *row_perm_out = NULL;
  *matched_out = 0;
  *exact_matching_out = 0;
  *spral_matching_out = 0;
  *row_scale_out = NULL;
  *col_scale_out = NULL;

  if (improve_matching && exact_weighted_match_is_affordable(n, nnz)) {
    UF_long *exact_row_perm = NULL;
    UF_long *exact_col_match = NULL;
    UF_long exact_matched = 0;
    const int exact_status =
      build_exact_numeric_row_match(n, nnz, col_ptr, row_idx, numeric_values,
                                    &exact_row_perm, &exact_col_match,
                                    &exact_matched);
    if (exact_status == KLS_OK) {
      const int complete_status =
        complete_unmatched_row_match(n, exact_row_perm, exact_col_match);
      free(exact_col_match);
      if (complete_status != KLS_OK) {
        free(exact_row_perm);
        return complete_status;
      }
      *row_perm_out = exact_row_perm;
      *matched_out = exact_matched;
      *exact_matching_out = 1;
      return KLS_OK;
    }
    if (exact_status == KLS_ERR_OUT_OF_MEMORY) {
      return exact_status;
    }
  }

  kls_match_entry *entries =
    (kls_match_entry *)malloc((size_t)nnz * sizeof(*entries));
  UF_long *row_perm = (UF_long *)malloc((size_t)n * sizeof(*row_perm));
  UF_long *col_match = (UF_long *)malloc((size_t)n * sizeof(*col_match));
  unsigned char *row_used = (unsigned char *)calloc((size_t)n, sizeof(*row_used));
  unsigned char *col_used = (unsigned char *)calloc((size_t)n, sizeof(*col_used));
  if (entries == NULL || row_perm == NULL || col_match == NULL ||
      row_used == NULL || col_used == NULL) {
    free(entries);
    free(row_perm);
    free(col_match);
    free(row_used);
    free(col_used);
    return KLS_ERR_OUT_OF_MEMORY;
  }

  UF_long entry_count = 0;
  for (UF_long col = 0; col < n; ++col) {
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      const double weight = fabs(numeric_values[p]);
      if (weight > 0.0 && isfinite(weight) && row_idx[p] < n) {
        entries[entry_count].weight = weight;
        entries[entry_count].row = row_idx[p];
        entries[entry_count].col = col;
        entry_count++;
      }
    }
  }
  qsort(entries, (size_t)entry_count, sizeof(*entries), compare_match_entries_desc);

  for (UF_long row = 0; row < n; ++row) {
    row_perm[row] = KLS_KLU_EMPTY;
    col_match[row] = KLS_KLU_EMPTY;
  }

  UF_long matched = 0;
  for (UF_long i = 0; i < entry_count && matched < n; ++i) {
    const UF_long row = entries[i].row;
    const UF_long col = entries[i].col;
    if (!row_used[row] && !col_used[col]) {
      row_used[row] = 1;
      col_used[col] = 1;
      row_perm[row] = col;
      col_match[col] = row;
      matched++;
    }
  }

  if (improve_matching) {
    kls_row_match_graph graph;
    int status = build_row_match_graph(n, entry_count, entries, &graph);
    if (status == KLS_OK) {
      if (n >= 50000u && n - matched >= 1024u) {
        matched = augment_numeric_row_match_layered(n, &graph, row_perm,
                                                    col_match);
      } else {
        matched = augment_numeric_row_match(n, &graph, row_perm, col_match);
      }
      if (matched == n) {
        improve_numeric_row_match_by_swaps(n, &graph, row_perm, col_match);
        improve_numeric_row_match_by_paths(n, &graph, row_perm, col_match);
      }
#ifdef KLS_HAVE_SPRAL_SCALING
      if (matched < n && n >= 4000u && n <= 150000u && nnz <= 1500000u) {
        UF_long *spral_row_perm = NULL;
        UF_long spral_matched = 0;
        double *spral_row_scale = NULL;
        double *spral_col_scale = NULL;
        const int spral_status =
          build_spral_hungarian_row_match_scaling(n, nnz, col_ptr, row_idx,
                                                  numeric_values,
                                                  &spral_row_perm,
                                                  &spral_matched,
                                                  &spral_row_scale,
                                                  &spral_col_scale);
        if (spral_status == KLS_OK) {
          const int better_cardinality =
            spral_matched > matched && matched < n;
          if (better_cardinality) {
            free(row_perm);
            row_perm = spral_row_perm;
            matched = spral_matched;
            spral_row_perm = NULL;
            for (UF_long col = 0; col < n; ++col) {
              col_match[col] = KLS_KLU_EMPTY;
            }
            for (UF_long row = 0; row < n; ++row) {
              const UF_long col = row_perm[row];
              if (col >= n || col_match[col] != KLS_KLU_EMPTY) {
                free(entries);
                free(row_perm);
                free(col_match);
                free(row_used);
                free(col_used);
                free(spral_row_scale);
                free(spral_col_scale);
                return KLS_ERR_INVALID_ARGUMENT;
              }
              col_match[col] = row;
            }
            if (spral_row_scale != NULL && spral_col_scale != NULL) {
              *row_scale_out = spral_row_scale;
              *col_scale_out = spral_col_scale;
              spral_row_scale = NULL;
              spral_col_scale = NULL;
            }
            *exact_matching_out = 1;
            *spral_matching_out = 1;
          }
          free(spral_row_perm);
          free(spral_row_scale);
          free(spral_col_scale);
        } else if (spral_status == KLS_ERR_OUT_OF_MEMORY) {
          free(entries);
          free(row_perm);
          free(col_match);
          free(row_used);
          free(col_used);
          return spral_status;
        }
      }
#endif
      free_row_match_graph(&graph);
    } else if (status == KLS_ERR_OUT_OF_MEMORY) {
      free(entries);
      free(row_perm);
      free(col_match);
      free(row_used);
      free(col_used);
      return status;
    }
  }

#ifdef KLS_HAVE_SPRAL_SCALING
  if (improve_matching && n >= 200000u && n - matched >= 4096u &&
      1000.0 * (double)matched < 995.0 * (double)n) {
    UF_long *spral_row_perm = NULL;
    UF_long spral_matched = 0;
    const int spral_status =
      build_spral_auction_row_match(n, nnz, col_ptr, row_idx, numeric_values,
                                    &spral_row_perm, &spral_matched,
                                    row_scale_out, col_scale_out);
    if (spral_status == KLS_OK) {
      if (spral_matched > matched) {
        free(row_perm);
        row_perm = spral_row_perm;
        matched = spral_matched;
        *spral_matching_out = 1;
        for (UF_long col = 0; col < n; ++col) {
          col_match[col] = KLS_KLU_EMPTY;
        }
        for (UF_long row = 0; row < n; ++row) {
          const UF_long col = row_perm[row];
          if (col >= n || col_match[col] != KLS_KLU_EMPTY) {
            free(entries);
            free(row_perm);
            free(col_match);
            free(row_used);
            free(col_used);
            return KLS_ERR_INVALID_ARGUMENT;
          }
          col_match[col] = row;
        }
      } else {
        free(spral_row_perm);
        free(*row_scale_out);
        free(*col_scale_out);
        *row_scale_out = NULL;
        *col_scale_out = NULL;
      }
    } else if (spral_status == KLS_ERR_OUT_OF_MEMORY) {
      free(entries);
      free(row_perm);
      free(col_match);
      free(row_used);
      free(col_used);
      return spral_status;
    }
  }
#endif

  const int complete_status =
    complete_unmatched_row_match(n, row_perm, col_match);
  if (complete_status != KLS_OK) {
    free(entries);
    free(row_perm);
    free(col_match);
    free(row_used);
    free(col_used);
    return complete_status;
  }

  free(entries);
  free(col_match);
  free(row_used);
  free(col_used);
  *row_perm_out = row_perm;
  *matched_out = matched;
  return KLS_OK;
}

static int build_sorted_row_permuted_pattern(UF_long n,
                                             UF_long nnz,
                                             const UF_long *base_col_ptr,
                                             const UF_long *base_row_idx,
                                             const double *base_values,
                                             const UF_long *row_perm,
                                             const UF_long *input_to_base,
                                             UF_long **col_ptr_out,
                                             UF_long **row_idx_out,
                                             double **values_out,
                                             UF_long **input_to_csc_out) {
  if (n <= 0 || base_col_ptr == NULL || base_row_idx == NULL ||
      base_values == NULL || row_perm == NULL || col_ptr_out == NULL ||
      row_idx_out == NULL || values_out == NULL || input_to_csc_out == NULL) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  *col_ptr_out = NULL;
  *row_idx_out = NULL;
  *values_out = NULL;
  *input_to_csc_out = NULL;

  UF_long *col_ptr = (UF_long *)calloc((size_t)n + 1u, sizeof(*col_ptr));
  UF_long *row_idx = (UF_long *)malloc((size_t)nnz * sizeof(*row_idx));
  double *values = (double *)malloc((size_t)nnz * sizeof(*values));
  UF_long *input_to_csc = (UF_long *)malloc((size_t)nnz * sizeof(*input_to_csc));
  UF_long *base_to_trial = (UF_long *)malloc((size_t)nnz * sizeof(*base_to_trial));
  kls_row_permuted_entry *entries =
    (kls_row_permuted_entry *)malloc((size_t)nnz * sizeof(*entries));
  if (col_ptr == NULL || row_idx == NULL || values == NULL ||
      input_to_csc == NULL || base_to_trial == NULL || entries == NULL) {
    free(col_ptr);
    free(row_idx);
    free(values);
    free(input_to_csc);
    free(base_to_trial);
    free(entries);
    return KLS_ERR_OUT_OF_MEMORY;
  }

  UF_long count = 0;
  for (UF_long col = 0; col < n; ++col) {
    for (UF_long p = base_col_ptr[col]; p < base_col_ptr[col + 1]; ++p) {
      const UF_long row = base_row_idx[p];
      if (row >= n || row_perm[row] >= n) {
        free(col_ptr);
        free(row_idx);
        free(values);
        free(input_to_csc);
        free(base_to_trial);
        free(entries);
        return KLS_ERR_INVALID_ARGUMENT;
      }
      entries[count].col = col;
      entries[count].row = row_perm[row];
      entries[count].base_position = p;
      entries[count].value = base_values[p];
      count++;
    }
  }
  if (count != nnz) {
    free(col_ptr);
    free(row_idx);
    free(values);
    free(input_to_csc);
    free(base_to_trial);
    free(entries);
    return KLS_ERR_INVALID_ARGUMENT;
  }

  qsort(entries, (size_t)nnz, sizeof(*entries), compare_row_permuted_entries);
  for (UF_long p = 0; p < nnz; ++p) {
    col_ptr[entries[p].col + 1]++;
  }
  for (UF_long col = 0; col < n; ++col) {
    col_ptr[col + 1] += col_ptr[col];
  }
  for (UF_long p = 0; p < nnz; ++p) {
    row_idx[p] = entries[p].row;
    values[p] = entries[p].value;
    base_to_trial[entries[p].base_position] = p;
  }
  for (UF_long p = 0; p < nnz; ++p) {
    const UF_long base_position = input_to_base != NULL ? input_to_base[p] : p;
    if (base_position >= nnz) {
      free(col_ptr);
      free(row_idx);
      free(values);
      free(input_to_csc);
      free(base_to_trial);
      free(entries);
      return KLS_ERR_INVALID_ARGUMENT;
    }
    input_to_csc[p] = base_to_trial[base_position];
  }

  free(base_to_trial);
  free(entries);
  *col_ptr_out = col_ptr;
  *row_idx_out = row_idx;
  *values_out = values;
  *input_to_csc_out = input_to_csc;
  return KLS_OK;
}

static double clamp_matching_scale(double value) {
  const double min_scale = 1.0e-12;
  const double max_scale = 1.0e12;
  if (!isfinite(value) || value <= 0.0) {
    return 1.0;
  }
  if (value < min_scale) {
    return min_scale;
  }
  if (value > max_scale) {
    return max_scale;
  }
  return value;
}

static int build_matching_dual_scaling(UF_long n,
                                       const UF_long *base_col_ptr,
                                       const UF_long *base_row_idx,
                                       const double *base_values,
                                       const UF_long *row_perm,
                                       const double *diag,
                                       double *row_scale,
                                       double *col_scale) {
  if (n <= 0 || base_col_ptr == NULL || base_row_idx == NULL ||
      base_values == NULL || row_perm == NULL || diag == NULL ||
      row_scale == NULL || col_scale == NULL) {
    return 0;
  }

  double *matched_log = (double *)malloc((size_t)n * sizeof(*matched_log));
  double *col_potential = (double *)calloc((size_t)n, sizeof(*col_potential));
  if (matched_log == NULL || col_potential == NULL) {
    free(matched_log);
    free(col_potential);
    return 0;
  }

  for (UF_long i = 0; i < n; ++i) {
    if (diag[i] <= 0.0 || !isfinite(diag[i])) {
      free(matched_log);
      free(col_potential);
      return 0;
    }
    matched_log[i] = log(diag[i]);
  }

  const int max_rounds = n >= 200000 ? 32 : 64;
  const double update_tol = 1.0e-10;
  int converged = 0;
  /* MC64-style dual potentials: x[col] is -log(col_scale[col]).  If the
     matched diagonal is normalized to one, every nonmatched entry gives a
     lower-bound constraint on x[col]. */
  for (int round = 0; round < max_rounds; ++round) {
    UF_long changes = 0;
    for (UF_long col = 0; col < n; ++col) {
      for (UF_long p = base_col_ptr[col]; p < base_col_ptr[col + 1u]; ++p) {
        const UF_long row = base_row_idx[p];
        if (row >= n || row_perm[row] >= n) {
          continue;
        }
        const double abs_value = fabs(base_values[p]);
        if (abs_value <= 0.0 || !isfinite(abs_value)) {
          continue;
        }
        const UF_long matched_col = row_perm[row];
        const double lower_bound =
          col_potential[matched_col] + log(abs_value) - matched_log[matched_col];
        if (lower_bound > col_potential[col] + update_tol) {
          col_potential[col] = lower_bound;
          changes++;
        }
      }
    }

    double min_potential = DBL_MAX;
    double max_potential = -DBL_MAX;
    for (UF_long i = 0; i < n; ++i) {
      if (col_potential[i] < min_potential) {
        min_potential = col_potential[i];
      }
      if (col_potential[i] > max_potential) {
        max_potential = col_potential[i];
      }
    }
    if (!isfinite(min_potential) || !isfinite(max_potential) ||
        max_potential - min_potential > 120.0) {
      free(matched_log);
      free(col_potential);
      return 0;
    }
    if (min_potential != 0.0 && isfinite(min_potential)) {
      for (UF_long i = 0; i < n; ++i) {
        col_potential[i] -= min_potential;
      }
    }
    if (changes == 0) {
      converged = 1;
      break;
    }
  }

  double max_excess = -DBL_MAX;
  for (UF_long col = 0; col < n; ++col) {
    for (UF_long p = base_col_ptr[col]; p < base_col_ptr[col + 1u]; ++p) {
      const UF_long row = base_row_idx[p];
      if (row >= n || row_perm[row] >= n) {
        continue;
      }
      const double abs_value = fabs(base_values[p]);
      if (abs_value <= 0.0 || !isfinite(abs_value)) {
        continue;
      }
      const UF_long matched_col = row_perm[row];
      const double excess = log(abs_value) - matched_log[matched_col] -
                            col_potential[col] + col_potential[matched_col];
      if (excess > max_excess) {
        max_excess = excess;
      }
    }
  }
  if (!converged && max_excess > log(4.0)) {
    free(matched_log);
    free(col_potential);
    return 0;
  }

  double row_log_sum = 0.0;
  double col_log_sum = 0.0;
  for (UF_long i = 0; i < n; ++i) {
    row_log_sum += col_potential[i] - matched_log[i];
    col_log_sum += -col_potential[i];
  }
  const double shift =
    0.5 * (col_log_sum - row_log_sum) / (double)n;
  const double max_log_scale = log(1.0e12);
  for (UF_long i = 0; i < n; ++i) {
    const double row_log_scale = col_potential[i] - matched_log[i] + shift;
    const double col_log_scale = -col_potential[i] - shift;
    if (!isfinite(row_log_scale) || !isfinite(col_log_scale) ||
        fabs(row_log_scale) > max_log_scale ||
        fabs(col_log_scale) > max_log_scale) {
      free(matched_log);
      free(col_potential);
      return 0;
    }
    row_scale[i] = clamp_matching_scale(exp(row_log_scale));
    col_scale[i] = clamp_matching_scale(exp(col_log_scale));
  }

  free(matched_log);
  free(col_potential);
  return 1;
}

static int build_matching_equilibration(UF_long n,
                                        const UF_long *base_col_ptr,
                                        const UF_long *base_row_idx,
                                        const double *base_values,
                                        const UF_long *row_perm,
                                        double **row_scale_out,
                                        double **col_scale_out) {
  if (n <= 0 || base_col_ptr == NULL || base_row_idx == NULL ||
      base_values == NULL || row_perm == NULL ||
      row_scale_out == NULL || col_scale_out == NULL) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  *row_scale_out = NULL;
  *col_scale_out = NULL;

  double *diag = (double *)calloc((size_t)n, sizeof(*diag));
  double *row_scale = (double *)malloc((size_t)n * sizeof(*row_scale));
  double *col_scale = (double *)malloc((size_t)n * sizeof(*col_scale));
  double *row_max = (double *)malloc((size_t)n * sizeof(*row_max));
  if (diag == NULL || row_scale == NULL || col_scale == NULL || row_max == NULL) {
    free(diag);
    free(row_scale);
    free(col_scale);
    free(row_max);
    return KLS_ERR_OUT_OF_MEMORY;
  }

  for (UF_long i = 0; i < n; ++i) {
    row_scale[i] = 1.0;
    col_scale[i] = 1.0;
  }

  for (UF_long col = 0; col < n; ++col) {
    for (UF_long p = base_col_ptr[col]; p < base_col_ptr[col + 1u]; ++p) {
      const UF_long row = base_row_idx[p];
      if (row < n && row_perm[row] == col) {
        const double value_abs = fabs(base_values[p]);
        if (isfinite(value_abs) && value_abs > diag[col]) {
          diag[col] = value_abs;
        }
      }
    }
  }

  double min_diag = DBL_MAX;
  double max_diag = 0.0;
  for (UF_long i = 0; i < n; ++i) {
    if (diag[i] <= 0.0 || !isfinite(diag[i])) {
      free(diag);
      free(row_scale);
      free(col_scale);
      free(row_max);
      return KLS_OK;
    }
    if (diag[i] < min_diag) {
      min_diag = diag[i];
    }
    if (diag[i] > max_diag) {
      max_diag = diag[i];
    }
  }

  if (min_diag <= 0.0 || max_diag / min_diag < 1.0e4) {
    free(diag);
    free(row_scale);
    free(col_scale);
    free(row_max);
    return KLS_OK;
  }

  if (!build_matching_dual_scaling(n, base_col_ptr, base_row_idx, base_values,
                                   row_perm, diag, row_scale, col_scale)) {
    for (UF_long i = 0; i < n; ++i) {
      const double scale = clamp_matching_scale(exp(-0.5 * log(diag[i])));
      row_scale[i] = scale;
      col_scale[i] = scale;
    }

    for (int round = 0; round < 6; ++round) {
      for (UF_long i = 0; i < n; ++i) {
        row_max[i] = 0.0;
      }
      for (UF_long col = 0; col < n; ++col) {
        const double cs = col_scale[col];
        for (UF_long p = base_col_ptr[col]; p < base_col_ptr[col + 1u]; ++p) {
          const UF_long row = base_row_idx[p];
          if (row >= n || row_perm[row] >= n) {
            continue;
          }
          const UF_long scaled_row = row_perm[row];
          const double scaled = fabs(base_values[p]) * row_scale[scaled_row] * cs;
          if (isfinite(scaled) && scaled > row_max[scaled_row]) {
            row_max[scaled_row] = scaled;
          }
        }
      }

      UF_long changed = 0;
      for (UF_long i = 0; i < n; ++i) {
        if (row_max[i] > 4.0) {
          const double factor = sqrt(row_max[i]);
          row_scale[i] = clamp_matching_scale(row_scale[i] / factor);
          col_scale[i] = clamp_matching_scale(col_scale[i] * factor);
          changed++;
        }
      }
      if (changed == 0) {
        break;
      }
    }
  }

  free(diag);
  free(row_max);
  *row_scale_out = row_scale;
  *col_scale_out = col_scale;
  return KLS_OK;
}

static void apply_value_scaling(UF_long n,
                                const UF_long *col_ptr,
                                const UF_long *row_idx,
                                const double *row_scale,
                                const double *col_scale,
                                double *values) {
  if (n <= 0 || col_ptr == NULL || row_idx == NULL || values == NULL ||
      (row_scale == NULL && col_scale == NULL)) {
    return;
  }
  for (UF_long col = 0; col < n; ++col) {
    const double cs = col_scale != NULL ? col_scale[col] : 1.0;
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      const UF_long row = row_idx[p];
      const double rs = (row_scale != NULL && row < n) ? row_scale[row] : 1.0;
      values[p] *= rs * cs;
    }
  }
}

static int matching_equilibration_is_better(
  const trilinos_klu_l_common *base_common,
  const trilinos_klu_l_numeric *base_numeric,
  const trilinos_klu_l_common *scaled_common,
  const trilinos_klu_l_numeric *scaled_numeric) {
  const double base_fill = (double)(base_numeric->lnz + base_numeric->unz);
  const double scaled_fill = (double)(scaled_numeric->lnz + scaled_numeric->unz);
  const double base_flops = base_common->flops;
  const double scaled_flops = scaled_common->flops;

  if (scaled_common->noffdiag > base_common->noffdiag ||
      scaled_fill > 1.02 * base_fill ||
      (base_flops > 0.0 && scaled_flops > 1.05 * base_flops)) {
    return 0;
  }
  if (base_common->rcond <= 0.0 && scaled_common->rcond > 0.0) {
    return 1;
  }
  if (base_common->rcond > 0.0 && scaled_common->rcond > 10.0 * base_common->rcond) {
    return 1;
  }
  if (base_common->scale > 0 &&
      (base_common->rcond <= 0.0 ||
       scaled_common->rcond >= 0.10 * base_common->rcond)) {
    return 1;
  }
  return scaled_fill < 0.98 * base_fill ||
         scaled_common->noffdiag < base_common->noffdiag;
}

static void maybe_use_matching_equilibration(
  UF_long n,
  UF_long nnz,
  UF_long *col_ptr,
  UF_long *row_idx,
  const kls_options *options,
  trilinos_klu_l_symbolic *symbolic,
  double **values_io,
  double **row_scale_io,
  double **col_scale_io,
  trilinos_klu_l_numeric **numeric_io,
  trilinos_klu_l_common *common_io) {
  if (n <= 0 || nnz <= 0 || col_ptr == NULL || row_idx == NULL ||
      options == NULL || symbolic == NULL || values_io == NULL ||
      row_scale_io == NULL || col_scale_io == NULL || numeric_io == NULL ||
      common_io == NULL || *values_io == NULL || *row_scale_io == NULL ||
      *col_scale_io == NULL || *numeric_io == NULL) {
    return;
  }

  double *scaled_values = (double *)malloc((size_t)nnz * sizeof(*scaled_values));
  if (scaled_values == NULL) {
    free(*row_scale_io);
    free(*col_scale_io);
    *row_scale_io = NULL;
    *col_scale_io = NULL;
    return;
  }
  memcpy(scaled_values, *values_io, (size_t)nnz * sizeof(*scaled_values));
  apply_value_scaling(n, col_ptr, row_idx, *row_scale_io, *col_scale_io,
                      scaled_values);

  kls_options scaled_options = *options;
  scaled_options.scale = -1;
  trilinos_klu_l_common scaled_common;
  if (apply_options_to_common(&scaled_common, &scaled_options) != KLS_OK) {
    free(scaled_values);
    free(*row_scale_io);
    free(*col_scale_io);
    *row_scale_io = NULL;
    *col_scale_io = NULL;
    return;
  }

  trilinos_klu_l_numeric *scaled_numeric =
    trilinos_klu_l_factor(col_ptr, row_idx, scaled_values, symbolic,
                          &scaled_common);
  if (scaled_numeric == NULL || scaled_common.status < 0 ||
      scaled_common.status == TRILINOS_KLU_SINGULAR) {
    if (scaled_numeric != NULL) {
      trilinos_klu_l_free_numeric(&scaled_numeric, &scaled_common);
    }
    free(scaled_values);
    free(*row_scale_io);
    free(*col_scale_io);
    *row_scale_io = NULL;
    *col_scale_io = NULL;
    return;
  }

  (void)trilinos_klu_l_flops(symbolic, scaled_numeric, &scaled_common);
  (void)trilinos_klu_l_rcond(symbolic, scaled_numeric, &scaled_common);
  if (!matching_equilibration_is_better(common_io, *numeric_io,
                                        &scaled_common, scaled_numeric)) {
    trilinos_klu_l_free_numeric(&scaled_numeric, &scaled_common);
    free(scaled_values);
    free(*row_scale_io);
    free(*col_scale_io);
    *row_scale_io = NULL;
    *col_scale_io = NULL;
    return;
  }

  trilinos_klu_l_numeric *old_numeric = *numeric_io;
  trilinos_klu_l_common old_common = *common_io;
  trilinos_klu_l_free_numeric(&old_numeric, &old_common);
  free(*values_io);
  *values_io = scaled_values;
  *numeric_io = scaled_numeric;
  *common_io = scaled_common;
}

#ifdef KLS_HAVE_SPRAL_SCALING
static int maybe_accept_spral_hungarian_numeric_trial(
  kls_solver *solver,
  const UF_long *base_col_ptr,
  const UF_long *base_row_idx,
  const double *base_values,
  const UF_long *base_input_to_csc) {
  if (solver == NULL || solver->numeric == NULL || solver->symbolic == NULL ||
      base_col_ptr == NULL || base_row_idx == NULL || base_values == NULL) {
    return 0;
  }

  UF_long *row_perm = NULL;
  UF_long matched = 0;
  double *row_scale = NULL;
  double *col_scale = NULL;
  UF_long *trial_col_ptr = NULL;
  UF_long *trial_row_idx = NULL;
  UF_long *trial_input_to_csc = NULL;
  double *trial_values = NULL;
  trilinos_klu_l_symbolic *trial_symbolic = NULL;
  trilinos_klu_l_numeric *trial_numeric = NULL;
  trilinos_klu_l_common trial_common;
  (void)trilinos_klu_l_defaults(&trial_common);
  int accepted = 0;

  int status =
    build_spral_hungarian_row_match_scaling(solver->n, solver->nnz,
                                            base_col_ptr, base_row_idx,
                                            base_values, &row_perm, &matched,
                                            &row_scale, &col_scale);
  if (status != KLS_OK || matched != solver->n ||
      row_scale == NULL || col_scale == NULL) {
    goto done;
  }

  status = build_sorted_row_permuted_pattern(solver->n, solver->nnz,
                                             base_col_ptr, base_row_idx,
                                             base_values, row_perm,
                                             base_input_to_csc,
                                             &trial_col_ptr, &trial_row_idx,
                                             &trial_values,
                                             &trial_input_to_csc);
  if (status != KLS_OK) {
    goto done;
  }

  kls_options trial_options = solver->options;
  kls_ordering trial_ordering = KLS_ORDERING_AUTO;
  double trial_score = 0.0;
  status = choose_symbolic_for_pattern(solver->n, trial_col_ptr, trial_row_idx,
                                       &trial_options, &trial_symbolic,
                                       &trial_common, &trial_ordering,
                                       &trial_score);
  if (status != KLS_OK) {
    goto done;
  }
  trial_common.scale = choose_auto_scale_from_pattern(solver->n, trial_col_ptr,
                                                      trial_row_idx,
                                                      &trial_options,
                                                      trial_values);
  trial_common.tol = solver->common.tol;

  trial_numeric =
    trilinos_klu_l_factor(trial_col_ptr, trial_row_idx, trial_values,
                          trial_symbolic, &trial_common);
  if (trial_numeric == NULL || trial_common.status < 0 ||
      trial_common.status == TRILINOS_KLU_SINGULAR) {
    goto done;
  }

  (void)trilinos_klu_l_flops(trial_symbolic, trial_numeric, &trial_common);
  (void)trilinos_klu_l_rcond(trial_symbolic, trial_numeric, &trial_common);
  maybe_use_matching_equilibration(solver->n, solver->nnz, trial_col_ptr,
                                   trial_row_idx, &solver->options,
                                   trial_symbolic, &trial_values,
                                   &row_scale, &col_scale, &trial_numeric,
                                   &trial_common);
  if (!numeric_candidate_is_better(&solver->common, solver->numeric,
                                   &trial_common, trial_numeric) ||
      !spral_hungarian_candidate_has_value(&solver->common, solver->numeric,
                                           &trial_common, trial_numeric) ||
      (solver->common.rcond > 0.0 && trial_common.rcond > 0.0 &&
       trial_common.rcond < 0.01 * solver->common.rcond)) {
    goto done;
  }

  trilinos_klu_l_symbolic *old_symbolic = solver->symbolic;
  trilinos_klu_l_numeric *old_numeric = solver->numeric;
  trilinos_klu_l_common old_common = solver->common;
  UF_long *old_col_ptr = solver->col_ptr;
  UF_long *old_row_idx = solver->row_idx;
  UF_long *old_input_to_csc = solver->input_to_csc;
  UF_long *old_row_perm = solver->row_perm;
  double *old_row_scale = solver->row_scale;
  double *old_col_scale = solver->col_scale;
  double *old_values = solver->values;

  solver->col_ptr = trial_col_ptr;
  solver->row_idx = trial_row_idx;
  solver->input_to_csc = trial_input_to_csc;
  solver->row_perm = row_perm;
  solver->row_scale = row_scale;
  solver->col_scale = col_scale;
  solver->values = trial_values;
  solver->orientation = KLS_ORIENTATION_NORMAL;
  solver->exact_matching_selected = 1;
  solver->spral_matching_selected = 1;
  solver->symbolic = trial_symbolic;
  solver->numeric = trial_numeric;
  solver->common = trial_common;
  solver->stats.selected_ordering = trial_ordering;
  solver->stats.selected_orientation = solver->orientation;
  solver->stats.nblocks = (int64_t)solver->symbolic->nblocks;
  solver->stats.max_block = (int64_t)solver->symbolic->maxblock;
  solver->stats.structural_rank = (int64_t)solver->symbolic->structural_rank;
  solver->stats.estimated_flops = solver->symbolic->est_flops;

  trilinos_klu_l_free_numeric(&old_numeric, &old_common);
  trilinos_klu_l_free_symbolic(&old_symbolic, &old_common);
  free(old_col_ptr);
  free(old_row_idx);
  free(old_input_to_csc);
  free(old_row_perm);
  free(old_row_scale);
  free(old_col_scale);
  free(old_values);

  trial_col_ptr = NULL;
  trial_row_idx = NULL;
  trial_input_to_csc = NULL;
  row_perm = NULL;
  row_scale = NULL;
  col_scale = NULL;
  trial_values = NULL;
  trial_symbolic = NULL;
  trial_numeric = NULL;
  accepted = 1;

done:
  if (!accepted) {
    if (trial_numeric != NULL) {
      trilinos_klu_l_free_numeric(&trial_numeric, &trial_common);
    }
    if (trial_symbolic != NULL) {
      trilinos_klu_l_free_symbolic(&trial_symbolic, &trial_common);
    }
    free(row_perm);
    free(row_scale);
    free(col_scale);
    free(trial_col_ptr);
    free(trial_row_idx);
    free(trial_input_to_csc);
    free(trial_values);
  }
  return accepted;
}

static int should_try_spral_hungarian_numeric_trial(
  const kls_solver *solver) {
  if (solver == NULL || !solver->options.static_pivoting ||
      solver->numeric == NULL || solver->symbolic == NULL ||
      solver->row_perm != NULL || solver->input_format != KLS_INPUT_CSC ||
      solver->options.ordering != KLS_ORDERING_AUTO ||
      solver->n < 20000u || solver->n > 750000u ||
      solver->nnz > 8000000u || solver->common.noffdiag < 16u) {
    return 0;
  }
  const UF_long fill = solver->numeric->lnz + solver->numeric->unz;
  const double offdiag_ratio =
    (double)solver->common.noffdiag / (double)solver->n;
  if (solver->common.noffdiag < 1024u && offdiag_ratio < 0.02) {
    return 0;
  }
  if (solver->common.flops >= 2.0e8 && fill >= 2000000u) {
    return 1;
  }
  return solver->common.noffdiag >= 1024u || offdiag_ratio >= 0.02;
}

static int maybe_select_spral_hungarian_row_match(
  kls_solver *solver,
  double *elapsed,
  const double *numeric_values) {
  if (elapsed == NULL || numeric_values == NULL ||
      !should_try_spral_hungarian_numeric_trial(solver)) {
    return 0;
  }

  const double start = kls_now_seconds();
  const UF_long *base_col_ptr = solver->col_ptr;
  const UF_long *base_row_idx = solver->row_idx;
  const UF_long *base_input_to_csc = solver->input_to_csc;
  const double *base_values = numeric_values;
  UF_long *owned_col_ptr = NULL;
  UF_long *owned_row_idx = NULL;
  UF_long *owned_input_to_csc = NULL;
  double *owned_values = NULL;
  int accepted = 0;

  if (solver->orientation == KLS_ORIENTATION_TRANSPOSE) {
    if (solver->input_to_csc == NULL) {
      goto done;
    }
    kls_pattern_candidate source = {0};
    source.n = solver->n;
    source.nnz = solver->nnz;
    source.col_ptr = solver->col_ptr;
    source.row_idx = solver->row_idx;
    source.orientation = solver->orientation;

    kls_pattern_candidate normal = {0};
    int status = transpose_candidate(&source, KLS_ORIENTATION_NORMAL, &normal);
    if (status != KLS_OK) {
      goto done;
    }
    owned_input_to_csc =
      (UF_long *)malloc((size_t)solver->nnz * sizeof(*owned_input_to_csc));
    if (owned_input_to_csc == NULL) {
      free_candidate(&normal);
      goto done;
    }
    for (UF_long p = 0; p < solver->nnz; ++p) {
      const UF_long current_p = solver->input_to_csc[p];
      const UF_long normal_p = normal.input_to_csc[current_p];
      owned_input_to_csc[p] = normal_p;
      normal.values[normal_p] = numeric_values[current_p];
    }
    owned_col_ptr = normal.col_ptr;
    owned_row_idx = normal.row_idx;
    owned_values = normal.values;
    base_col_ptr = owned_col_ptr;
    base_row_idx = owned_row_idx;
    base_values = owned_values;
    base_input_to_csc = owned_input_to_csc;
    normal.col_ptr = NULL;
    normal.row_idx = NULL;
    normal.values = NULL;
    free_candidate(&normal);
  }

  accepted =
    maybe_accept_spral_hungarian_numeric_trial(solver, base_col_ptr,
                                               base_row_idx, base_values,
                                               base_input_to_csc);

done:
  free(owned_col_ptr);
  free(owned_row_idx);
  free(owned_input_to_csc);
  free(owned_values);
  *elapsed += kls_now_seconds() - start;
  return accepted;
}
#endif

static UF_long count_weak_diagonal_rows(UF_long n,
                                        const UF_long *col_ptr,
                                        const UF_long *row_idx,
                                        const double *values,
                                        double tolerance,
                                        UF_long *missing_diagonal_out) {
  double *row_max = (double *)calloc((size_t)n, sizeof(*row_max));
  double *diag = (double *)calloc((size_t)n, sizeof(*diag));
  if (missing_diagonal_out != NULL) {
    *missing_diagonal_out = 0;
  }
  if (row_max == NULL || diag == NULL) {
    free(row_max);
    free(diag);
    return 0;
  }

  for (UF_long col = 0; col < n; ++col) {
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1]; ++p) {
      const UF_long row = row_idx[p];
      if (row >= n) {
        continue;
      }
      const double abs_value = fabs(values[p]);
      if (abs_value > row_max[row]) {
        row_max[row] = abs_value;
      }
      if (row == col && abs_value > diag[row]) {
        diag[row] = abs_value;
      }
    }
  }

  UF_long weak = 0;
  UF_long missing = 0;
  for (UF_long row = 0; row < n; ++row) {
    if (diag[row] == 0.0) {
      missing++;
    }
    if (row_max[row] > 0.0 && diag[row] < tolerance * row_max[row]) {
      weak++;
    }
  }
  if (missing_diagonal_out != NULL) {
    *missing_diagonal_out = missing;
  }
  free(row_max);
  free(diag);
  return weak;
}

static int static_match_prefers_unscaled(UF_long n,
                                         UF_long nnz,
                                         UF_long weak_diagonal,
                                         UF_long missing_diagonal) {
  /* Majority-missing medium circuit blocks can lose sparsity from matching
     equilibration; keep the static row permutation but leave values unscaled. */
  return n >= 80000u && n <= 150000u && nnz <= 1500000u &&
         missing_diagonal * 2u >= n &&
         weak_diagonal * 2u >= n;
}

static int reactive_static_match_setup_is_unlikely_to_pay(
  const kls_solver *solver,
  UF_long weak_diagonal,
  UF_long missing_diagonal) {
  if (solver == NULL || solver->symbolic == NULL || solver->n > 20000u ||
      solver->common.flops <= 0.0) {
    return 0;
  }

  const UF_long n = solver->n;
  const int severe_weak_diagonal =
    weak_diagonal * 10u >= 9u * n || missing_diagonal * 10u >= 9u * n;
  if (severe_weak_diagonal) {
    return 0;
  }

  if (solver->symbolic->nblocks == 1u && solver->symbolic->maxblock == n &&
      solver->common.flops < 3.0e6 && weak_diagonal * 5u < 3u * n) {
    return 1;
  }

  if (n <= 12000u && solver->symbolic->nblocks >= 512u &&
      solver->common.flops < 1.0e6 && weak_diagonal * 20u < n &&
      missing_diagonal * 20u < n && solver->common.noffdiag < 512u) {
    return 1;
  }

  if (solver->stats.selected_ordering == KLS_ORDERING_METIS &&
      solver->symbolic->nblocks >= 512u &&
      solver->common.flops < 2.0e6 &&
      weak_diagonal * 5u < n &&
      solver->common.noffdiag < 256u) {
    return 1;
  }

  return 0;
}

#ifdef KLS_HAVE_METIS
static int pre_static_metis_refinement_is_worthwhile(
  UF_long n,
  UF_long nnz,
  UF_long weak_diagonal,
  const trilinos_klu_l_symbolic *symbolic,
  const trilinos_klu_l_common *common,
  kls_ordering selected_ordering) {
  if (selected_ordering == KLS_ORDERING_METIS || symbolic == NULL ||
      common == NULL || n < 50000u || n > 150000u ||
      nnz > 1500000u || weak_diagonal * 100u < n ||
      !symbolic->do_btf || symbolic->nblocks < 1024u ||
      symbolic->maxblock == 0u ||
      symbolic->maxblock < (UF_long)(0.90 * (double)n) ||
      symbolic->maxblock >= (UF_long)(0.99 * (double)n) ||
      common->flops < 1.0e8) {
    return 0;
  }
  return 1;
}

static int maybe_refine_pre_static_with_metis(
  UF_long n,
  UF_long nnz,
  UF_long weak_diagonal,
  UF_long *trial_col_ptr,
  UF_long *trial_row_idx,
  double *trial_values,
  const kls_options *trial_options,
  trilinos_klu_l_symbolic **trial_symbolic_io,
  trilinos_klu_l_numeric **trial_numeric_io,
  trilinos_klu_l_common *trial_common_io,
  kls_ordering *trial_ordering_io,
  double *trial_score_io) {
  if (trial_col_ptr == NULL || trial_row_idx == NULL ||
      trial_values == NULL || trial_options == NULL ||
      trial_symbolic_io == NULL || *trial_symbolic_io == NULL ||
      trial_numeric_io == NULL || *trial_numeric_io == NULL ||
      trial_common_io == NULL || trial_ordering_io == NULL ||
      trial_score_io == NULL ||
      !pre_static_metis_refinement_is_worthwhile(
        n, nnz, weak_diagonal, *trial_symbolic_io, trial_common_io,
        *trial_ordering_io)) {
    return 0;
  }

  kls_options metis_options = *trial_options;
  metis_options.ordering = KLS_ORDERING_METIS;
  metis_options.scale = (int)trial_common_io->scale;
  metis_options.pivot_tolerance = trial_common_io->tol;

  trilinos_klu_l_symbolic *metis_symbolic = NULL;
  trilinos_klu_l_common metis_common;
  int status = analyze_with_ordering(n, trial_col_ptr, trial_row_idx,
                                     &metis_options, KLS_ORDERING_METIS,
                                     &metis_symbolic, &metis_common);
  if (status != KLS_OK) {
    return 0;
  }

  trilinos_klu_l_numeric *metis_numeric =
    trilinos_klu_l_factor(trial_col_ptr, trial_row_idx, trial_values,
                          metis_symbolic, &metis_common);
  if (metis_numeric == NULL || metis_common.status < 0 ||
      metis_common.status == TRILINOS_KLU_SINGULAR) {
    if (metis_numeric != NULL) {
      trilinos_klu_l_free_numeric(&metis_numeric, &metis_common);
    }
    trilinos_klu_l_free_symbolic(&metis_symbolic, &metis_common);
    return 0;
  }

  (void)trilinos_klu_l_flops(metis_symbolic, metis_numeric, &metis_common);
  (void)trilinos_klu_l_rcond(metis_symbolic, metis_numeric, &metis_common);
  if (!numeric_candidate_is_better(trial_common_io, *trial_numeric_io,
                                   &metis_common, metis_numeric) ||
      (trial_common_io->rcond > 0.0 && metis_common.rcond > 0.0 &&
       metis_common.rcond < 0.01 * trial_common_io->rcond)) {
    trilinos_klu_l_free_numeric(&metis_numeric, &metis_common);
    trilinos_klu_l_free_symbolic(&metis_symbolic, &metis_common);
    return 0;
  }

  trilinos_klu_l_free_numeric(trial_numeric_io, trial_common_io);
  trilinos_klu_l_free_symbolic(trial_symbolic_io, trial_common_io);
  *trial_symbolic_io = metis_symbolic;
  *trial_numeric_io = metis_numeric;
  *trial_common_io = metis_common;
  *trial_ordering_io = KLS_ORDERING_METIS;
  *trial_score_io = symbolic_score(metis_symbolic);
  return 1;
}
#endif

static int should_try_auto_row_match(const kls_solver *solver,
                                     const double *numeric_values) {
  /* Weighted static pivoting is currently a reactive medium-matrix trial.  On
     larger matrices the O(nnz log nnz) matching/reanalysis cost needs a cheaper
     precheck before it is worth paying by default. */
  if (solver == NULL || !solver->options.static_pivoting ||
      solver->numeric == NULL || solver->row_perm != NULL ||
      solver->input_format != KLS_INPUT_CSC ||
      solver->options.ordering != KLS_ORDERING_AUTO || solver->n < 3000 ||
      solver->n > 20000 || solver->common.noffdiag < 16 ||
      numeric_values == NULL) {
    return 0;
  }
#ifdef KLS_HAVE_METIS
  if (is_small_spiked_low_diagonal_pattern(solver->n, solver->col_ptr,
                                           solver->row_idx)) {
    return 0;
  }
#endif
  UF_long missing_diagonal = 0;
  const UF_long weak_diagonal =
    count_weak_diagonal_rows(solver->n, solver->col_ptr, solver->row_idx,
                             numeric_values, solver->common.tol,
                             &missing_diagonal);
  if (reactive_static_match_setup_is_unlikely_to_pay(
        solver, weak_diagonal, missing_diagonal)) {
    return 0;
  }
  const UF_long fill = solver->numeric->lnz + solver->numeric->unz;
  if (fill < 50000 || solver->common.flops < 5.0e5) {
    return 0;
  }
  const double offdiag_ratio = (double)solver->common.noffdiag / (double)solver->n;
  return solver->common.noffdiag >= 128 || offdiag_ratio >= 0.005;
}

static int maybe_select_auto_row_match(kls_solver *solver,
                                       double *elapsed,
                                       const double *numeric_values) {
  if (!should_try_auto_row_match(solver, numeric_values)) {
    return 0;
  }

  const double start = kls_now_seconds();
  const UF_long *base_col_ptr = solver->col_ptr;
  const UF_long *base_row_idx = solver->row_idx;
  const UF_long *base_input_to_csc = solver->input_to_csc;
  const double *base_values = numeric_values;
  UF_long *owned_col_ptr = NULL;
  UF_long *owned_row_idx = NULL;
  UF_long *owned_input_to_csc = NULL;
  double *owned_values = NULL;
  UF_long *row_perm = NULL;
  UF_long *trial_col_ptr = NULL;
  UF_long *trial_row_idx = NULL;
  UF_long *trial_input_to_csc = NULL;
  double *trial_values = NULL;
  double *trial_row_scale = NULL;
  double *trial_col_scale = NULL;
  double *match_row_scale = NULL;
  double *match_col_scale = NULL;
  trilinos_klu_l_symbolic *trial_symbolic = NULL;
  trilinos_klu_l_numeric *trial_numeric = NULL;
  trilinos_klu_l_common trial_common;
  (void)trilinos_klu_l_defaults(&trial_common);
  int accepted = 0;

  if (solver->orientation == KLS_ORIENTATION_TRANSPOSE) {
    if (solver->input_to_csc == NULL) {
      goto done;
    }
    kls_pattern_candidate source = {0};
    source.n = solver->n;
    source.nnz = solver->nnz;
    source.col_ptr = solver->col_ptr;
    source.row_idx = solver->row_idx;
    source.orientation = solver->orientation;

    kls_pattern_candidate normal = {0};
    int status = transpose_candidate(&source, KLS_ORIENTATION_NORMAL, &normal);
    if (status != KLS_OK) {
      goto done;
    }
    owned_input_to_csc =
      (UF_long *)malloc((size_t)solver->nnz * sizeof(*owned_input_to_csc));
    if (owned_input_to_csc == NULL) {
      free_candidate(&normal);
      goto done;
    }
    for (UF_long p = 0; p < solver->nnz; ++p) {
      const UF_long current_p = solver->input_to_csc[p];
      const UF_long normal_p = normal.input_to_csc[current_p];
      owned_input_to_csc[p] = normal_p;
      normal.values[normal_p] = numeric_values[current_p];
    }
    owned_col_ptr = normal.col_ptr;
    owned_row_idx = normal.row_idx;
    owned_values = normal.values;
    base_col_ptr = owned_col_ptr;
    base_row_idx = owned_row_idx;
    base_values = owned_values;
    base_input_to_csc = owned_input_to_csc;
    normal.col_ptr = NULL;
    normal.row_idx = NULL;
    normal.values = NULL;
    free_candidate(&normal);
  }

  UF_long matched = 0;
  int exact_matching = 0;
  int spral_matching = 0;
  const int improve_matching = solver->n <= 50000 && solver->nnz <= 1000000;
  int status = build_greedy_numeric_row_match(solver->n, solver->nnz,
                                              base_col_ptr, base_row_idx,
                                              base_values, improve_matching,
                                              &row_perm, &matched,
                                              &exact_matching,
                                              &spral_matching,
                                              &match_row_scale,
                                              &match_col_scale);
  if (status != KLS_OK) {
    goto done;
  }
  if (1000u * matched < 995u * solver->n) {
    goto done;
  }

  status = build_sorted_row_permuted_pattern(solver->n, solver->nnz,
                                             base_col_ptr, base_row_idx,
                                             base_values, row_perm,
                                             base_input_to_csc,
                                             &trial_col_ptr, &trial_row_idx,
                                             &trial_values,
                                             &trial_input_to_csc);
  if (status != KLS_OK) {
    goto done;
  }

  kls_options trial_options = solver->options;
  if (trial_options.scale == KLS_SCALE_AUTO) {
    if (match_row_scale != NULL && match_col_scale != NULL) {
      trial_row_scale = match_row_scale;
      trial_col_scale = match_col_scale;
      match_row_scale = NULL;
      match_col_scale = NULL;
    } else {
      status = build_matching_equilibration(solver->n, base_col_ptr,
                                            base_row_idx, base_values,
                                            row_perm, &trial_row_scale,
                                            &trial_col_scale);
      if (status != KLS_OK) {
        goto done;
      }
    }
  }
  kls_ordering trial_ordering = KLS_ORDERING_AUTO;
  double trial_score = 0.0;
  status = choose_symbolic_for_pattern(solver->n, trial_col_ptr, trial_row_idx,
                                       &trial_options, &trial_symbolic,
                                       &trial_common, &trial_ordering,
                                       &trial_score);
  if (status != KLS_OK) {
    goto done;
  }
  trial_common.scale = choose_auto_scale_from_pattern(solver->n, trial_col_ptr,
                                                      trial_row_idx,
                                                      &trial_options,
                                                      trial_values);

  trial_numeric =
    trilinos_klu_l_factor(trial_col_ptr, trial_row_idx, trial_values,
                          trial_symbolic, &trial_common);
  if (trial_numeric == NULL || trial_common.status < 0 ||
      trial_common.status == TRILINOS_KLU_SINGULAR) {
    goto done;
  }

  (void)trilinos_klu_l_flops(trial_symbolic, trial_numeric, &trial_common);
  (void)trilinos_klu_l_rcond(trial_symbolic, trial_numeric, &trial_common);
  maybe_use_matching_equilibration(solver->n, solver->nnz, trial_col_ptr,
                                   trial_row_idx, &solver->options,
                                   trial_symbolic, &trial_values,
                                   &trial_row_scale, &trial_col_scale,
                                   &trial_numeric, &trial_common);
  if (!numeric_candidate_is_better(&solver->common, solver->numeric,
                                   &trial_common, trial_numeric) ||
      (solver->common.rcond > 0.0 && trial_common.rcond > 0.0 &&
       trial_common.rcond < 0.01 * solver->common.rcond)) {
    goto done;
  }

  trilinos_klu_l_symbolic *old_symbolic = solver->symbolic;
  trilinos_klu_l_numeric *old_numeric = solver->numeric;
  trilinos_klu_l_common old_common = solver->common;
  UF_long *old_col_ptr = solver->col_ptr;
  UF_long *old_row_idx = solver->row_idx;
  UF_long *old_input_to_csc = solver->input_to_csc;
  double *old_row_scale = solver->row_scale;
  double *old_col_scale = solver->col_scale;
  double *old_values = solver->values;

  solver->col_ptr = trial_col_ptr;
  solver->row_idx = trial_row_idx;
  solver->input_to_csc = trial_input_to_csc;
  solver->row_scale = trial_row_scale;
  solver->col_scale = trial_col_scale;
  solver->values = trial_values;
  solver->orientation = KLS_ORIENTATION_NORMAL;
  solver->row_perm = row_perm;
  solver->exact_matching_selected = exact_matching;
  solver->spral_matching_selected = spral_matching;
  solver->symbolic = trial_symbolic;
  solver->numeric = trial_numeric;
  solver->common = trial_common;
  solver->stats.selected_ordering = trial_ordering;
  solver->stats.selected_orientation = solver->orientation;
  solver->stats.nblocks = (int64_t)solver->symbolic->nblocks;
  solver->stats.max_block = (int64_t)solver->symbolic->maxblock;
  solver->stats.structural_rank = (int64_t)solver->symbolic->structural_rank;
  solver->stats.estimated_flops = solver->symbolic->est_flops;

  trilinos_klu_l_free_numeric(&old_numeric, &old_common);
  trilinos_klu_l_free_symbolic(&old_symbolic, &old_common);
  free(old_col_ptr);
  free(old_row_idx);
  free(old_input_to_csc);
  free(old_row_scale);
  free(old_col_scale);
  free(old_values);

  trial_col_ptr = NULL;
  trial_row_idx = NULL;
  trial_input_to_csc = NULL;
  trial_values = NULL;
  trial_row_scale = NULL;
  trial_col_scale = NULL;
  row_perm = NULL;
  trial_symbolic = NULL;
  trial_numeric = NULL;
  free(match_row_scale);
  free(match_col_scale);
  match_row_scale = NULL;
  match_col_scale = NULL;
  accepted = 1;

done:
  if (!accepted) {
    if (trial_numeric != NULL) {
      trilinos_klu_l_free_numeric(&trial_numeric, &trial_common);
    }
    if (trial_symbolic != NULL) {
      trilinos_klu_l_free_symbolic(&trial_symbolic, &trial_common);
    }
    free(row_perm);
    free(trial_col_ptr);
    free(trial_row_idx);
    free(trial_input_to_csc);
    free(trial_values);
    free(trial_row_scale);
    free(trial_col_scale);
    free(match_row_scale);
    free(match_col_scale);
  }
  free(owned_col_ptr);
  free(owned_row_idx);
  free(owned_input_to_csc);
  free(owned_values);
  *elapsed += kls_now_seconds() - start;
  return accepted;
}

static void maybe_select_pre_static_row_match(kls_solver *solver,
                                              double *elapsed,
                                              const double *numeric_values) {
  if (solver == NULL || !solver->options.static_pivoting ||
      solver->numeric != NULL || solver->row_perm != NULL ||
      solver->input_format != KLS_INPUT_CSC ||
      solver->options.ordering != KLS_ORDERING_AUTO ||
      solver->n < 3000) {
    return;
  }
#ifdef KLS_HAVE_METIS
  if (is_small_spiked_low_diagonal_pattern(solver->n, solver->col_ptr,
                                           solver->row_idx)) {
    return;
  }
#endif

  const int small_candidate = solver->n <= 20000u;
  const int medium_weak_candidate =
    solver->n <= 150000u && solver->nnz <= 1500000u;
#ifdef KLS_HAVE_SPRAL_SCALING
  const int large_spral_candidate =
    solver->n > 150000u && solver->n <= 750000u &&
    solver->nnz <= 8000000u;
#else
  const int large_spral_candidate = 0;
#endif
  if (!small_candidate && !medium_weak_candidate &&
      !large_spral_candidate) {
    return;
  }

  const double start = kls_now_seconds();
  const UF_long *base_col_ptr = solver->col_ptr;
  const UF_long *base_row_idx = solver->row_idx;
  const UF_long *base_input_to_csc = solver->input_to_csc;
  const double *base_values = numeric_values;
  UF_long *owned_col_ptr = NULL;
  UF_long *owned_row_idx = NULL;
  UF_long *owned_input_to_csc = NULL;
  double *owned_values = NULL;
  UF_long *row_perm = NULL;
  UF_long *trial_col_ptr = NULL;
  UF_long *trial_row_idx = NULL;
  UF_long *trial_input_to_csc = NULL;
  double *trial_values = NULL;
  double *trial_row_scale = NULL;
  double *trial_col_scale = NULL;
  double *match_row_scale = NULL;
  double *match_col_scale = NULL;
  trilinos_klu_l_symbolic *trial_symbolic = NULL;
  trilinos_klu_l_numeric *trial_numeric = NULL;
  trilinos_klu_l_common trial_common;
  (void)trilinos_klu_l_defaults(&trial_common);
  int accepted = 0;

  if (solver->orientation == KLS_ORIENTATION_TRANSPOSE) {
    if (solver->input_to_csc == NULL) {
      goto done;
    }
    kls_pattern_candidate source = {0};
    source.n = solver->n;
    source.nnz = solver->nnz;
    source.col_ptr = solver->col_ptr;
    source.row_idx = solver->row_idx;
    source.orientation = solver->orientation;

    kls_pattern_candidate normal = {0};
    int status = transpose_candidate(&source, KLS_ORIENTATION_NORMAL, &normal);
    if (status != KLS_OK) {
      goto done;
    }
    owned_input_to_csc =
      (UF_long *)malloc((size_t)solver->nnz * sizeof(*owned_input_to_csc));
    if (owned_input_to_csc == NULL) {
      free_candidate(&normal);
      goto done;
    }
    for (UF_long p = 0; p < solver->nnz; ++p) {
      const UF_long current_p = solver->input_to_csc[p];
      const UF_long normal_p = normal.input_to_csc[current_p];
      owned_input_to_csc[p] = normal_p;
      normal.values[normal_p] = numeric_values[current_p];
    }
    owned_col_ptr = normal.col_ptr;
    owned_row_idx = normal.row_idx;
    owned_values = normal.values;
    base_col_ptr = owned_col_ptr;
    base_row_idx = owned_row_idx;
    base_values = owned_values;
    base_input_to_csc = owned_input_to_csc;
    normal.col_ptr = NULL;
    normal.row_idx = NULL;
    normal.values = NULL;
    free_candidate(&normal);
  }

  UF_long missing_diagonal = 0;
  const UF_long weak =
    count_weak_diagonal_rows(solver->n, base_col_ptr, base_row_idx,
                             base_values, solver->common.tol,
                             &missing_diagonal);
  const int prefer_unscaled_static_match =
    static_match_prefers_unscaled(solver->n, solver->nnz, weak,
                                  missing_diagonal);
#ifdef KLS_HAVE_SPRAL_SCALING
  int use_large_spral_match = 0;
#endif
  if (!small_candidate) {
    const int majority_weak = weak * 2u >= solver->n && weak >= 5000;
    const int partial_weak =
      solver->n >= 80000 && solver->n <= 100000 &&
      solver->nnz >= 5u * solver->n && solver->nnz <= 1000000 &&
      missing_diagonal * 100u <= 3u * solver->n &&
      weak >= 2000 && weak * 100u >= solver->n;
    int use_medium_gate = 1;
#ifdef KLS_HAVE_SPRAL_SCALING
    if (!medium_weak_candidate) {
      const int large_spral_weak =
        large_spral_candidate &&
        solver->symbolic != NULL && solver->symbolic->do_btf &&
        solver->symbolic->maxblock * 4u >= solver->n * 3u &&
        ((weak >= 32768u && weak * 100u >= solver->n) ||
         (missing_diagonal >= 4096u &&
          missing_diagonal * 200u >= solver->n));
      if (!large_spral_weak) {
        goto done;
      }
      use_large_spral_match = 1;
      use_medium_gate = 0;
    }
#endif
    if (use_medium_gate && !majority_weak && !partial_weak) {
      goto done;
    }
  } else if (weak * 2u < solver->n) {
    goto done;
  }

  UF_long matched = 0;
  int exact_matching = 0;
  int spral_matching = 0;
#ifdef KLS_HAVE_SPRAL_SCALING
  int status = KLS_OK;
  if (use_large_spral_match) {
    status = build_spral_auction_row_match(
      solver->n, solver->nnz, base_col_ptr, base_row_idx, base_values,
      &row_perm, &matched, &match_row_scale, &match_col_scale);
    exact_matching = 0;
    spral_matching = status == KLS_OK;
  } else
#else
  int status = KLS_OK;
#endif
  {
    const int improve_matching =
      small_candidate || medium_weak_candidate;
    status = build_greedy_numeric_row_match(solver->n, solver->nnz,
                                            base_col_ptr, base_row_idx,
                                            base_values, improve_matching,
                                            &row_perm,
                                            &matched,
                                            &exact_matching,
                                            &spral_matching,
                                            &match_row_scale,
                                            &match_col_scale);
  }
  if (status != KLS_OK || 1000u * matched < 995u * solver->n) {
    goto done;
  }

  status = build_sorted_row_permuted_pattern(solver->n, solver->nnz,
                                             base_col_ptr, base_row_idx,
                                             base_values, row_perm,
                                             base_input_to_csc,
                                             &trial_col_ptr, &trial_row_idx,
                                             &trial_values,
                                             &trial_input_to_csc);
  if (status != KLS_OK) {
    goto done;
  }

  kls_options trial_options = solver->options;
  if (trial_options.scale == KLS_SCALE_AUTO && !prefer_unscaled_static_match) {
    if (match_row_scale != NULL && match_col_scale != NULL) {
      trial_row_scale = match_row_scale;
      trial_col_scale = match_col_scale;
      match_row_scale = NULL;
      match_col_scale = NULL;
    } else {
      status = build_matching_equilibration(solver->n, base_col_ptr,
                                            base_row_idx, base_values,
                                            row_perm, &trial_row_scale,
                                            &trial_col_scale);
      if (status != KLS_OK) {
        goto done;
      }
    }
  }
  kls_ordering trial_ordering = KLS_ORDERING_AUTO;
  double trial_score = 0.0;
  status = choose_symbolic_for_pattern(solver->n, trial_col_ptr, trial_row_idx,
                                       &trial_options, &trial_symbolic,
                                       &trial_common, &trial_ordering,
                                       &trial_score);
  if (status != KLS_OK) {
    goto done;
  }
  trial_common.scale = choose_auto_scale_from_pattern(solver->n, trial_col_ptr,
                                                      trial_row_idx,
                                                      &trial_options,
                                                      trial_values);
  if (trial_options.scale == KLS_SCALE_AUTO && prefer_unscaled_static_match) {
    trial_common.scale = -1;
  }

  trial_numeric =
    trilinos_klu_l_factor(trial_col_ptr, trial_row_idx, trial_values,
                          trial_symbolic, &trial_common);
  if (trial_numeric == NULL || trial_common.status < 0 ||
      trial_common.status == TRILINOS_KLU_SINGULAR) {
    goto done;
  }

  (void)trilinos_klu_l_flops(trial_symbolic, trial_numeric, &trial_common);
  (void)trilinos_klu_l_rcond(trial_symbolic, trial_numeric, &trial_common);
#ifdef KLS_HAVE_METIS
  (void)maybe_refine_pre_static_with_metis(
    solver->n, solver->nnz, weak, trial_col_ptr, trial_row_idx, trial_values,
    &trial_options, &trial_symbolic, &trial_numeric, &trial_common,
    &trial_ordering, &trial_score);
#endif
  maybe_use_matching_equilibration(solver->n, solver->nnz, trial_col_ptr,
                                   trial_row_idx, &solver->options,
                                   trial_symbolic, &trial_values,
                                   &trial_row_scale, &trial_col_scale,
                                   &trial_numeric, &trial_common);
  if (trial_common.noffdiag > weak / 20u + 16u ||
      trial_common.rcond <= 0.0) {
    goto done;
  }

  trilinos_klu_l_symbolic *old_symbolic = solver->symbolic;
  trilinos_klu_l_common old_common = solver->common;
  free(solver->col_ptr);
  free(solver->row_idx);
  free(solver->input_to_csc);
  free(solver->row_scale);
  free(solver->col_scale);
  free(solver->values);

  solver->col_ptr = trial_col_ptr;
  solver->row_idx = trial_row_idx;
  solver->input_to_csc = trial_input_to_csc;
  solver->row_scale = trial_row_scale;
  solver->col_scale = trial_col_scale;
  solver->values = trial_values;
  solver->orientation = KLS_ORIENTATION_NORMAL;
  solver->row_perm = row_perm;
  solver->exact_matching_selected = exact_matching;
  solver->spral_matching_selected = spral_matching;
  solver->symbolic = trial_symbolic;
  solver->numeric = trial_numeric;
  solver->common = trial_common;
  solver->stats.selected_ordering = trial_ordering;
  solver->stats.selected_orientation = solver->orientation;
  solver->stats.nblocks = (int64_t)solver->symbolic->nblocks;
  solver->stats.max_block = (int64_t)solver->symbolic->maxblock;
  solver->stats.structural_rank = (int64_t)solver->symbolic->structural_rank;
  solver->stats.estimated_flops = solver->symbolic->est_flops;
  trilinos_klu_l_free_symbolic(&old_symbolic, &old_common);

  trial_col_ptr = NULL;
  trial_row_idx = NULL;
  trial_input_to_csc = NULL;
  trial_values = NULL;
  trial_row_scale = NULL;
  trial_col_scale = NULL;
  row_perm = NULL;
  trial_symbolic = NULL;
  trial_numeric = NULL;
  free(match_row_scale);
  free(match_col_scale);
  match_row_scale = NULL;
  match_col_scale = NULL;
  accepted = 1;

done:
  if (!accepted) {
    if (trial_numeric != NULL) {
      trilinos_klu_l_free_numeric(&trial_numeric, &trial_common);
    }
    if (trial_symbolic != NULL) {
      trilinos_klu_l_free_symbolic(&trial_symbolic, &trial_common);
    }
    free(row_perm);
    free(trial_col_ptr);
    free(trial_row_idx);
    free(trial_input_to_csc);
    free(trial_values);
    free(trial_row_scale);
    free(trial_col_scale);
    free(match_row_scale);
    free(match_col_scale);
  }
  free(owned_col_ptr);
  free(owned_row_idx);
  free(owned_input_to_csc);
  free(owned_values);
  *elapsed += kls_now_seconds() - start;
}

static int should_try_auto_scale(const kls_solver *solver) {
  if (solver->auto_scale_checked || solver->options.scale != KLS_SCALE_AUTO ||
      solver->numeric == NULL || solver->n < 20000) {
    return 0;
  }
  if (solver->common.scale <= 0 &&
      is_medium_dense_diagonal_high_degree_pattern(solver->n, solver->col_ptr,
                                                   solver->row_idx)) {
    return 0;
  }
#ifdef KLS_HAVE_METIS
  if (is_large_very_low_degree_full_diagonal_pattern(solver->n, solver->col_ptr,
                                                     solver->row_idx) ||
      is_large_sparse_diagonal_low_degree_pattern(solver->n, solver->col_ptr,
                                                  solver->row_idx)) {
    return 0;
  }
#endif

#ifdef KLS_HAVE_METIS
  if (solver->common.scale == 1 &&
      is_medium_spiked_low_diagonal_pattern(solver->n, solver->col_ptr,
                                            solver->row_idx)) {
    return 0;
  }
#endif

  const double flops = solver->common.flops;
  const UF_long fill = solver->numeric->lnz + solver->numeric->unz;
  if (solver->options.ordering == KLS_ORDERING_AUTO &&
      solver->stats.selected_ordering == KLS_ORDERING_METIS &&
      solver->symbolic != NULL && !solver->symbolic->do_btf &&
      solver->symbolic->nblocks == 1u && solver->common.scale == 2 &&
      solver->n >= 200000u && flops >= 1.0e9 && fill >= 10000000u) {
    return 0;
  }
  return flops >= 1.0e8 && fill >= 1500000;
}

static int maybe_select_auto_scale(kls_solver *solver,
                                   double *elapsed,
                                   const double *numeric_values) {
  if (!should_try_auto_scale(solver)) {
    return 0;
  }
  solver->auto_scale_checked = 1;

  int accepted = 0;
  const int candidates[] = {-1, 1, 2};
  for (size_t i = 0; i < sizeof(candidates) / sizeof(candidates[0]); ++i) {
    if (candidates[i] == (int)solver->common.scale) {
      continue;
    }

    kls_options trial_options = solver->options;
    trial_options.scale = candidates[i];
    trilinos_klu_l_common trial_common;
    if (apply_options_to_common(&trial_common, &trial_options) != KLS_OK) {
      continue;
    }

    const double start = kls_now_seconds();
    trilinos_klu_l_numeric *trial_numeric =
      trilinos_klu_l_factor(solver->col_ptr, solver->row_idx, (double *)numeric_values,
                            solver->symbolic, &trial_common);
    *elapsed += kls_now_seconds() - start;
    if (trial_numeric == NULL || trial_common.status < 0 ||
        trial_common.status == TRILINOS_KLU_SINGULAR) {
      if (trial_numeric != NULL) {
        trilinos_klu_l_free_numeric(&trial_numeric, &trial_common);
      }
      continue;
    }

    (void)trilinos_klu_l_flops(solver->symbolic, trial_numeric, &trial_common);
    if (!numeric_candidate_is_better(&solver->common, solver->numeric,
                                     &trial_common, trial_numeric)) {
      trilinos_klu_l_free_numeric(&trial_numeric, &trial_common);
      continue;
    }

    trilinos_klu_l_numeric *old_numeric = solver->numeric;
    trilinos_klu_l_common old_common = solver->common;
    solver->numeric = trial_numeric;
    solver->common = trial_common;
    trilinos_klu_l_free_numeric(&old_numeric, &old_common);
    accepted = 1;
  }
  return accepted;
}

static int should_try_auto_pivot_tolerance(const kls_solver *solver) {
  if (solver->auto_pivot_checked || solver->numeric == NULL || solver->n < 30000 ||
      fabs(solver->options.pivot_tolerance - 0.001) > 1.0e-12 ||
      fabs(solver->common.tol - 0.001) > 1.0e-12 ||
      solver->common.noffdiag < 16) {
    return 0;
  }

  const UF_long fill = solver->numeric->lnz + solver->numeric->unz;
  return fill >= 1000000;
}

static int pivot_tolerance_numeric_is_better(
  const trilinos_klu_l_common *current_common,
  const trilinos_klu_l_numeric *current_numeric,
  const trilinos_klu_l_common *candidate_common,
  const trilinos_klu_l_numeric *candidate_numeric) {
  const double current_fill = (double)(current_numeric->lnz + current_numeric->unz);
  const double candidate_fill = (double)(candidate_numeric->lnz + candidate_numeric->unz);
  const double current_offdiag = (double)current_common->noffdiag;
  const double candidate_offdiag = (double)candidate_common->noffdiag;

  if (candidate_fill > current_fill || candidate_offdiag > 0.75 * current_offdiag) {
    return 0;
  }
  if (current_common->rcond > 0.0 && candidate_common->rcond > 0.0 &&
      candidate_common->rcond < 0.01 * current_common->rcond) {
    return 0;
  }
  return candidate_fill <= 0.98 * current_fill ||
         candidate_offdiag <= 0.50 * current_offdiag;
}

static int maybe_select_auto_pivot_tolerance(kls_solver *solver,
                                             double *elapsed,
                                             const double *numeric_values) {
  if (!should_try_auto_pivot_tolerance(solver)) {
    return 0;
  }
  solver->auto_pivot_checked = 1;

  kls_options trial_options = solver->options;
  trial_options.scale = (int)solver->common.scale;
  trial_options.pivot_tolerance = 1.0e-4;
  trilinos_klu_l_common trial_common;
  if (apply_options_to_common(&trial_common, &trial_options) != KLS_OK) {
    return 0;
  }

  const double start = kls_now_seconds();
  trilinos_klu_l_numeric *trial_numeric =
    trilinos_klu_l_factor(solver->col_ptr, solver->row_idx, (double *)numeric_values,
                          solver->symbolic, &trial_common);
  *elapsed += kls_now_seconds() - start;
  if (trial_numeric == NULL || trial_common.status < 0 ||
      trial_common.status == TRILINOS_KLU_SINGULAR) {
    if (trial_numeric != NULL) {
      trilinos_klu_l_free_numeric(&trial_numeric, &trial_common);
    }
    return 0;
  }

  (void)trilinos_klu_l_flops(solver->symbolic, trial_numeric, &trial_common);
  (void)trilinos_klu_l_rcond(solver->symbolic, trial_numeric, &trial_common);
  if (!pivot_tolerance_numeric_is_better(&solver->common, solver->numeric,
                                         &trial_common, trial_numeric)) {
    trilinos_klu_l_free_numeric(&trial_numeric, &trial_common);
    return 0;
  }

  trilinos_klu_l_numeric *old_numeric = solver->numeric;
  trilinos_klu_l_common old_common = solver->common;
  solver->numeric = trial_numeric;
  solver->common = trial_common;
  trilinos_klu_l_free_numeric(&old_numeric, &old_common);
  return 1;
}

#ifdef KLS_HAVE_METIS
static int should_try_auto_metis(const kls_solver *solver) {
  if (solver->auto_metis_checked || solver->options.ordering != KLS_ORDERING_AUTO ||
      solver->stats.selected_ordering == KLS_ORDERING_METIS ||
      solver->numeric == NULL) {
    return 0;
  }

  const double flops = solver->common.flops;
  const UF_long fill = solver->numeric->lnz + solver->numeric->unz;
  if (solver->symbolic != NULL && solver->symbolic->do_btf &&
      solver->symbolic->nblocks <= 4 &&
      (double)solver->symbolic->maxblock >= 0.90 * (double)solver->n &&
      solver->common.noffdiag >= 128 &&
      flops >= 1.0e7 && fill >= 150000) {
    return 1;
  }
  if (solver->n < 20000) {
    return 0;
  }
  return flops >= 1.0e8 && fill >= 1000000;
}

static int metis_numeric_is_better(const kls_solver *solver,
                                   const trilinos_klu_l_common *metis_common,
                                   const trilinos_klu_l_numeric *metis_numeric) {
  return numeric_candidate_is_better(&solver->common, solver->numeric,
                                     metis_common, metis_numeric);
}

static int maybe_promote_auto_metis(kls_solver *solver,
                                    double *elapsed,
                                    const double *numeric_values) {
  if (!should_try_auto_metis(solver)) {
    return 0;
  }
  solver->auto_metis_checked = 1;

  kls_options metis_options = solver->options;
  metis_options.ordering = KLS_ORDERING_METIS;
  metis_options.scale = (int)solver->common.scale;
  if (solver->symbolic != NULL) {
    metis_options.use_btf = solver->symbolic->do_btf ? 1 : 0;
  }

  trilinos_klu_l_symbolic *metis_symbolic = NULL;
  trilinos_klu_l_common metis_common;
  double start = kls_now_seconds();
  int status = analyze_with_ordering(solver->n, solver->col_ptr, solver->row_idx,
                                     &metis_options, KLS_ORDERING_METIS,
                                     &metis_symbolic, &metis_common);
  *elapsed += kls_now_seconds() - start;
  if (status != KLS_OK) {
    return 0;
  }

  start = kls_now_seconds();
  trilinos_klu_l_numeric *metis_numeric =
    trilinos_klu_l_factor(solver->col_ptr, solver->row_idx, (double *)numeric_values,
                          metis_symbolic, &metis_common);
  *elapsed += kls_now_seconds() - start;
  if (metis_numeric == NULL || metis_common.status < 0 ||
      metis_common.status == TRILINOS_KLU_SINGULAR) {
    if (metis_numeric != NULL) {
      trilinos_klu_l_free_numeric(&metis_numeric, &metis_common);
    }
    trilinos_klu_l_free_symbolic(&metis_symbolic, &metis_common);
    return 0;
  }

  (void)trilinos_klu_l_flops(metis_symbolic, metis_numeric, &metis_common);
  if (!metis_numeric_is_better(solver, &metis_common, metis_numeric)) {
    trilinos_klu_l_free_numeric(&metis_numeric, &metis_common);
    trilinos_klu_l_free_symbolic(&metis_symbolic, &metis_common);
    return 0;
  }

  trilinos_klu_l_symbolic *old_symbolic = solver->symbolic;
  trilinos_klu_l_numeric *old_numeric = solver->numeric;
  trilinos_klu_l_common old_common = solver->common;

  solver->symbolic = metis_symbolic;
  solver->numeric = metis_numeric;
  solver->common = metis_common;
  solver->stats.selected_ordering = KLS_ORDERING_METIS;
  solver->stats.nblocks = (int64_t)solver->symbolic->nblocks;
  solver->stats.max_block = (int64_t)solver->symbolic->maxblock;
  solver->stats.structural_rank = (int64_t)solver->symbolic->structural_rank;
  solver->stats.estimated_flops = solver->symbolic->est_flops;

  trilinos_klu_l_free_numeric(&old_numeric, &old_common);
  trilinos_klu_l_free_symbolic(&old_symbolic, &old_common);
  return 1;
}
#endif

static int choose_symbolic_for_pattern(UF_long n,
                                       UF_long *col_ptr,
                                       UF_long *row_idx,
                                       const kls_options *options,
                                       trilinos_klu_l_symbolic **symbolic_out,
                                       trilinos_klu_l_common *common_out,
                                       kls_ordering *selected_ordering_out,
                                       double *score_out) {
  if (options->ordering != KLS_ORDERING_AUTO) {
    int status = analyze_with_ordering(n, col_ptr, row_idx, options,
                                       options->ordering, symbolic_out,
                                       common_out);
    if (status == KLS_OK) {
      *selected_ordering_out = options->ordering;
      *score_out = symbolic_score(*symbolic_out);
    }
    return status;
  }

#ifdef KLS_HAVE_METIS
  if (is_large_very_low_degree_full_diagonal_pattern(n, col_ptr, row_idx)) {
    int status = analyze_with_ordering(n, col_ptr, row_idx, options,
                                       KLS_ORDERING_METIS, symbolic_out,
                                       common_out);
    if (status == KLS_OK) {
      *selected_ordering_out = KLS_ORDERING_METIS;
      double selected_score = symbolic_score(*symbolic_out);
      maybe_retry_without_btf(n, col_ptr, row_idx, options,
                              KLS_ORDERING_METIS, symbolic_out, common_out,
                              &selected_score, 0);
      *score_out = selected_score;
      return KLS_OK;
    }
  }
#endif

  kls_options auto_options = *options;
  const kls_options *symbolic_options = options;
#ifdef KLS_HAVE_METIS
  const int large_spiked_metis_no_btf =
    is_large_nearly_diagonal_spiked_metis_pattern(n, col_ptr, row_idx);
#else
  const int large_spiked_metis_no_btf = 0;
#endif
  if (should_start_auto_without_btf(n, col_ptr, row_idx, options,
                                    large_spiked_metis_no_btf)) {
    auto_options.use_btf = 0;
    symbolic_options = &auto_options;
  }

#ifdef KLS_HAVE_METIS
  if (should_start_auto_with_metis(n, col_ptr, row_idx,
                                   large_spiked_metis_no_btf)) {
    int status = analyze_with_ordering(n, col_ptr, row_idx, symbolic_options,
                                       KLS_ORDERING_METIS, symbolic_out,
                                       common_out);
    if (status == KLS_OK) {
      *selected_ordering_out = KLS_ORDERING_METIS;
      double selected_score = symbolic_score(*symbolic_out);
      maybe_retry_without_btf(n, col_ptr, row_idx, symbolic_options,
                              KLS_ORDERING_METIS, symbolic_out, common_out,
                              &selected_score, 0);
      *score_out = selected_score;
      return KLS_OK;
    }
  }
#endif

  const kls_ordering candidates[] = {KLS_ORDERING_AMD, KLS_ORDERING_COLAMD};
  trilinos_klu_l_symbolic *best_symbolic = NULL;
  trilinos_klu_l_common best_common;
  kls_ordering best_ordering = KLS_ORDERING_AMD;
  double best_score = 0.0;
  int any_ok = 0;

  for (size_t i = 0; i < sizeof(candidates) / sizeof(candidates[0]); ++i) {
    trilinos_klu_l_symbolic *candidate_symbolic = NULL;
    trilinos_klu_l_common candidate_common;
    int status = analyze_with_ordering(n, col_ptr, row_idx, symbolic_options, candidates[i],
                                       &candidate_symbolic, &candidate_common);
    if (status != KLS_OK) {
      continue;
    }
    const double score = symbolic_score(candidate_symbolic);
    double selected_score = score;
    maybe_retry_without_btf(n, col_ptr, row_idx, symbolic_options,
                            candidates[i], &candidate_symbolic,
                            &candidate_common, &selected_score, 1);
    if (!any_ok || selected_score < best_score) {
      if (best_symbolic != NULL) {
        trilinos_klu_l_free_symbolic(&best_symbolic, &best_common);
      }
      best_symbolic = candidate_symbolic;
      best_common = candidate_common;
      best_ordering = candidates[i];
      best_score = selected_score;
      any_ok = 1;
    } else {
      trilinos_klu_l_free_symbolic(&candidate_symbolic, &candidate_common);
    }
    /* Stop once an ordering provides a real fill estimate; unknown estimates
       should not displace a known-good circuit ordering. */
    if (best_ordering == candidates[i] && best_score < DBL_MAX) {
      break;
    }
  }

  if (!any_ok) {
    return KLS_ERR_ANALYZE_FAILED;
  }

#ifdef KLS_HAVE_METIS
  if (should_try_symbolic_metis_before_numeric(n, best_symbolic,
                                               best_ordering, best_score)) {
    kls_options metis_options = *symbolic_options;
    metis_options.use_btf = best_symbolic->do_btf ? 1 : 0;
    trilinos_klu_l_symbolic *metis_symbolic = NULL;
    trilinos_klu_l_common metis_common;
    int status = analyze_with_ordering(n, col_ptr, row_idx, &metis_options,
                                       KLS_ORDERING_METIS, &metis_symbolic,
                                       &metis_common);
    if (status == KLS_OK) {
      double metis_score = symbolic_score(metis_symbolic);
      maybe_retry_without_btf(n, col_ptr, row_idx, &metis_options,
                              KLS_ORDERING_METIS, &metis_symbolic,
                              &metis_common, &metis_score, 0);
      if (isfinite(metis_score) && metis_score <= 0.90 * best_score) {
        trilinos_klu_l_free_symbolic(&best_symbolic, &best_common);
        best_symbolic = metis_symbolic;
        best_common = metis_common;
        best_ordering = KLS_ORDERING_METIS;
        best_score = metis_score;
      } else {
        trilinos_klu_l_free_symbolic(&metis_symbolic, &metis_common);
      }
    }
  }
#endif

  *symbolic_out = best_symbolic;
  *common_out = best_common;
  *selected_ordering_out = best_ordering;
  *score_out = best_score;
  return KLS_OK;
}

static int validate_index_base(int index_base) {
  return index_base == 0 || index_base == 1;
}

static int validate_n(int64_t n) {
  return n > 0 && n < INT64_MAX / 2 && n <= (int64_t)(UF_long_max / 2);
}

static int validate_options(const kls_options *options) {
  if (options->threads <= 0) {
    return 0;
  }
  if (options->ordering < KLS_ORDERING_AUTO || options->ordering > KLS_ORDERING_SCOTCH) {
    return 0;
  }
  if (options->orientation < KLS_ORIENTATION_AUTO ||
      options->orientation > KLS_ORIENTATION_TRANSPOSE) {
    return 0;
  }
  if (options->scale != KLS_SCALE_AUTO && (options->scale < -1 || options->scale > 2)) {
    return 0;
  }
  if (options->pivot_tolerance < 0.0) {
    return 0;
  }
  if (options->memory_growth <= 0.0) {
    return 0;
  }
  if (options->fast_factor != 0 && options->fast_factor != 1) {
    return 0;
  }
  if (options->static_pivoting != 0 && options->static_pivoting != 1) {
    return 0;
  }
  return 1;
}

static int allocate_candidate(kls_pattern_candidate *candidate,
                              int64_t n,
                              int64_t nnz,
                              kls_orientation orientation,
                              int need_map) {
  memset(candidate, 0, sizeof(*candidate));
  candidate->col_ptr = (UF_long *)calloc((size_t)n + 1u, sizeof(UF_long));
  candidate->row_idx = (UF_long *)calloc((size_t)nnz, sizeof(UF_long));
  if (need_map) {
    candidate->input_to_csc = (UF_long *)calloc((size_t)nnz, sizeof(UF_long));
    candidate->values = (double *)calloc((size_t)nnz, sizeof(double));
  }
  if (candidate->col_ptr == NULL || candidate->row_idx == NULL ||
      (need_map && (candidate->values == NULL || candidate->input_to_csc == NULL))) {
    free(candidate->col_ptr);
    free(candidate->row_idx);
    free(candidate->input_to_csc);
    free(candidate->values);
    memset(candidate, 0, sizeof(*candidate));
    return KLS_ERR_OUT_OF_MEMORY;
  }
  candidate->n = (UF_long)n;
  candidate->nnz = (UF_long)nnz;
  candidate->orientation = orientation;
  candidate->selected_ordering = KLS_ORDERING_AUTO;
  candidate->score = DBL_MAX;
  return KLS_OK;
}

static void free_candidate(kls_pattern_candidate *candidate) {
  if (candidate == NULL) {
    return;
  }
  if (candidate->symbolic != NULL) {
    trilinos_klu_l_free_symbolic(&candidate->symbolic, &candidate->common);
  }
  free(candidate->col_ptr);
  free(candidate->row_idx);
  free(candidate->input_to_csc);
  free(candidate->values);
  memset(candidate, 0, sizeof(*candidate));
}

static int copy_compressed_candidate(kls_pattern_candidate *candidate,
                                     kls_index_type index_type,
                                     int64_t n,
                                     const void *col_ptr,
                                     const void *row_idx,
                                     int index_base,
                                     kls_orientation orientation) {
  const int64_t base = index_base;
  const int64_t nnz = read_index(index_type, col_ptr, n) - base;
  if (nnz < 0) {
    return KLS_ERR_INVALID_ARGUMENT;
  }

  int status = allocate_candidate(candidate, n, nnz, orientation, 0);
  if (status != KLS_OK) {
    return status;
  }

  int64_t prev = 0;
  for (int64_t j = 0; j <= n; ++j) {
    const int64_t ptr = read_index(index_type, col_ptr, j) - base;
    if (ptr < prev || ptr < 0 || ptr > nnz) {
      free_candidate(candidate);
      return KLS_ERR_INVALID_ARGUMENT;
    }
    candidate->col_ptr[j] = (UF_long)ptr;
    prev = ptr;
  }
  for (int64_t p = 0; p < nnz; ++p) {
    const int64_t row = read_index(index_type, row_idx, p) - base;
    if (row < 0 || row >= n) {
      free_candidate(candidate);
      return KLS_ERR_INVALID_ARGUMENT;
    }
    candidate->row_idx[p] = (UF_long)row;
  }
  return KLS_OK;
}

static int transpose_candidate(const kls_pattern_candidate *source,
                               kls_orientation orientation,
                               kls_pattern_candidate *candidate) {
  int status = allocate_candidate(candidate, (int64_t)source->n,
                                  (int64_t)source->nnz, orientation, 1);
  if (status != KLS_OK) {
    return status;
  }

  for (UF_long p = 0; p < source->nnz; ++p) {
    candidate->col_ptr[source->row_idx[p] + 1]++;
  }
  for (UF_long j = 0; j < source->n; ++j) {
    candidate->col_ptr[j + 1] += candidate->col_ptr[j];
  }

  UF_long *next = (UF_long *)malloc(((size_t)source->n) * sizeof(UF_long));
  if (next == NULL) {
    free_candidate(candidate);
    return KLS_ERR_OUT_OF_MEMORY;
  }
  memcpy(next, candidate->col_ptr, ((size_t)source->n) * sizeof(UF_long));

  for (UF_long col = 0; col < source->n; ++col) {
    for (UF_long p = source->col_ptr[col]; p < source->col_ptr[col + 1]; ++p) {
      const UF_long dst_col = source->row_idx[p];
      const UF_long dst = next[dst_col]++;
      candidate->row_idx[dst] = col;
      candidate->input_to_csc[p] = dst;
    }
  }
  free(next);
  return KLS_OK;
}

static int analyze_candidate(kls_pattern_candidate *candidate,
                             const kls_options *options) {
  return choose_symbolic_for_pattern(candidate->n, candidate->col_ptr,
                                     candidate->row_idx, options,
                                     &candidate->symbolic, &candidate->common,
                                     &candidate->selected_ordering,
                                     &candidate->score);
}

static int auto_orientation_prefers_transpose(UF_long n) {
  /* On small and medium SPICE-like matrices, the second symbolic analysis
     usually costs more than the normal-vs-transpose fill estimate saves. */
  return n <= 30000;
}

static int auto_orientation_prefers_normal(UF_long n,
                                           const UF_long *col_ptr,
                                           const UF_long *row_idx) {
#ifdef KLS_HAVE_METIS
  if (is_large_sparse_diagonal_low_degree_pattern(n, col_ptr, row_idx)) {
    return 1;
  }
#endif
  return is_large_diagonal_circuit_like_pattern(n, col_ptr, row_idx);
}

static int select_candidate(kls_pattern_candidate *normal,
                            kls_pattern_candidate *transpose,
                            const kls_options *options,
                            kls_pattern_candidate **chosen_out) {
  *chosen_out = NULL;
  if (options->orientation == KLS_ORIENTATION_NORMAL) {
    if (normal == NULL) {
      return KLS_ERR_INVALID_ARGUMENT;
    }
    int status = analyze_candidate(normal, options);
    if (status == KLS_OK) {
      *chosen_out = normal;
    }
    return status;
  }
  if (options->orientation == KLS_ORIENTATION_TRANSPOSE) {
    if (transpose == NULL) {
      return KLS_ERR_INVALID_ARGUMENT;
    }
    int status = analyze_candidate(transpose, options);
    if (status == KLS_OK) {
      *chosen_out = transpose;
    }
    return status;
  }

  if (transpose != NULL && auto_orientation_prefers_transpose(transpose->n)) {
    const int transpose_status = analyze_candidate(transpose, options);
    if (transpose_status == KLS_OK) {
      *chosen_out = transpose;
      return KLS_OK;
    }
    if (normal != NULL) {
      const int normal_status = analyze_candidate(normal, options);
      if (normal_status == KLS_OK) {
        *chosen_out = normal;
        return KLS_OK;
      }
    }
    return transpose_status;
  }
  if (normal == NULL || transpose == NULL) {
    kls_pattern_candidate *candidate = normal != NULL ? normal : transpose;
    if (candidate == NULL) {
      return KLS_ERR_INVALID_ARGUMENT;
    }
    int status = analyze_candidate(candidate, options);
    if (status == KLS_OK) {
      *chosen_out = candidate;
    }
    return status;
  }
  const int normal_status = analyze_candidate(normal, options);
  const int transpose_status = analyze_candidate(transpose, options);
  if (normal_status == KLS_OK && transpose_status == KLS_OK) {
    *chosen_out = (transpose->score < normal->score) ? transpose : normal;
    return KLS_OK;
  }
  if (normal_status == KLS_OK) {
    *chosen_out = normal;
    return KLS_OK;
  }
  if (transpose_status == KLS_OK) {
    *chosen_out = transpose;
    return KLS_OK;
  }
  return normal_status != KLS_OK ? normal_status : transpose_status;
}

static void adopt_candidate(kls_solver *solver, kls_pattern_candidate *candidate) {
  solver->n = candidate->n;
  solver->nnz = candidate->nnz;
  solver->col_ptr = candidate->col_ptr;
  solver->row_idx = candidate->row_idx;
  solver->input_to_csc = candidate->input_to_csc;
  solver->values = candidate->values;
  solver->orientation = candidate->orientation;
  solver->common = candidate->common;
  solver->symbolic = candidate->symbolic;
  solver->stats.selected_ordering = candidate->selected_ordering;
  solver->stats.selected_orientation = candidate->orientation;

  candidate->col_ptr = NULL;
  candidate->row_idx = NULL;
  candidate->input_to_csc = NULL;
  candidate->values = NULL;
  candidate->symbolic = NULL;
}

static void fill_symbolic_stats(kls_solver *solver, double elapsed) {
  solver->stats.struct_size = sizeof(solver->stats);
  solver->stats.n = (int64_t)solver->n;
  solver->stats.nnz = (int64_t)solver->nnz;
  solver->stats.analysis_seconds = elapsed;
  solver->stats.selected_orientation = solver->orientation;
  solver->stats.selected_scale = (int)solver->common.scale;
  solver->stats.selected_pivot_tolerance = solver->common.tol;
  solver->stats.selected_static_pivoting = solver->row_perm != NULL;
  solver->stats.selected_exact_matching = solver->exact_matching_selected;
  solver->stats.selected_spral_matching = solver->spral_matching_selected;
  solver->stats.fast_block_restarts = solver->fast_block_restarts;
  solver->stats.fast_tail_restarts = solver->fast_tail_restarts;
  solver->stats.fast_rejected_refresh_state =
    solver->fast_reject_refresh_state;
  if (solver->symbolic != NULL) {
    solver->stats.last_kernel_status = (int)solver->common.status;
    solver->stats.selected_btf = solver->symbolic->do_btf ? 1 : 0;
    solver->stats.nblocks = (int64_t)solver->symbolic->nblocks;
    solver->stats.max_block = (int64_t)solver->symbolic->maxblock;
    solver->stats.structural_rank = (int64_t)solver->symbolic->structural_rank;
    solver->stats.estimated_flops = solver->symbolic->est_flops;
    solver->stats.nnz_l = (int64_t)solver->symbolic->lnz;
    solver->stats.nnz_u = (int64_t)solver->symbolic->unz;
  }
  solver->stats.memory_bytes = solver->common.memusage;
  solver->stats.memory_peak_bytes = solver->common.mempeak;
}

static void fill_numeric_stats(kls_solver *solver) {
  solver->stats.last_kernel_status = (int)solver->common.status;
  solver->stats.selected_scale = (int)solver->common.scale;
  solver->stats.selected_pivot_tolerance = solver->common.tol;
  solver->stats.selected_static_pivoting = solver->row_perm != NULL;
  solver->stats.selected_exact_matching = solver->exact_matching_selected;
  solver->stats.selected_spral_matching = solver->spral_matching_selected;
  solver->stats.fast_block_restarts = solver->fast_block_restarts;
  solver->stats.fast_tail_restarts = solver->fast_tail_restarts;
  solver->stats.fast_rejected_refresh_state =
    solver->fast_reject_refresh_state;
  solver->stats.selected_btf =
    (solver->symbolic != NULL && solver->symbolic->do_btf) ? 1 : 0;
  solver->stats.numerical_rank = (int64_t)solver->common.numerical_rank;
  solver->stats.singular_col = (int64_t)solver->common.singular_col;
  solver->stats.offdiag_pivots = (int64_t)solver->common.noffdiag;
  solver->stats.reallocations = (int64_t)solver->common.nrealloc;
  solver->stats.factor_flops = solver->common.flops;
  solver->stats.rcond = solver->common.rcond;
  solver->stats.condest = solver->common.condest;
  solver->stats.rgrowth = solver->common.rgrowth;
  solver->stats.memory_bytes = solver->common.memusage;
  solver->stats.memory_peak_bytes = solver->common.mempeak;
  if (solver->numeric != NULL) {
    solver->stats.nnz_l = (int64_t)solver->numeric->lnz;
    solver->stats.nnz_u = (int64_t)solver->numeric->unz;
  }
  solver->stats.refactor_dependency_levels =
    (int64_t)solver->refactor_level_count;
  solver->stats.refactor_dependency_max_width =
    (int64_t)solver->refactor_level_max_width;
  solver->stats.refactor_dependency_edges =
    (int64_t)solver->refactor_dependency_edges;
  solver->stats.refactor_dependency_root_columns =
    (int64_t)solver->refactor_dependency_root_columns;
  solver->stats.refactor_dependency_leaf_columns =
    (int64_t)solver->refactor_dependency_leaf_columns;
  solver->stats.refactor_dependency_max_fanout =
    (int64_t)solver->refactor_dependency_max_fanout;
  solver->stats.refactor_dependency_max_column_work =
    solver->refactor_dependency_max_column_work;
  solver->stats.refactor_dependency_pipeline_max_column_work =
    solver->refactor_dependency_pipeline_max_column_work;
  solver->stats.refactor_supernode_candidate_count =
    (int64_t)solver->refactor_supernode_candidate_count;
  solver->stats.refactor_supernode_candidate_rows =
    (int64_t)solver->refactor_supernode_candidate_rows;
  solver->stats.refactor_supernode_candidate_max_width =
    (int64_t)solver->refactor_supernode_candidate_max_width;
  solver->stats.refactor_supernode_candidate_dense_entries =
    solver->refactor_supernode_candidate_dense_entries;
  solver->stats.refactor_supernode_candidate_trailing_entries =
    solver->refactor_supernode_candidate_trailing_entries;
  solver->stats.row_refactor_group_count =
    (int64_t)solver->row_refactor_group_count;
  solver->stats.row_refactor_group_level_count =
    (int64_t)solver->row_refactor_level_count;
  solver->stats.row_refactor_group_level_max_width =
    (int64_t)solver->row_refactor_level_max_width;
  solver->stats.row_refactor_group_cluster_levels =
    (int64_t)solver->row_refactor_cluster_level_count;
  solver->stats.row_refactor_group_pipeline_groups =
    (int64_t)solver->row_refactor_pipeline_group_count;
  solver->stats.row_refactor_group_pipeline_rows =
    (int64_t)solver->row_refactor_pipeline_row_count;
  solver->stats.row_refactor_group_pipeline_work =
    solver->row_refactor_pipeline_work;
  solver->stats.row_refactor_last_run = solver->row_refactor_last_run;
  solver->stats.row_refactor_last_checked = solver->row_refactor_last_checked;
  solver->stats.row_refactor_last_parallel = solver->row_refactor_last_parallel;
  solver->stats.row_refactor_run_count =
    (int64_t)solver->row_refactor_run_count;
  solver->stats.row_refactor_checked_run_count =
    (int64_t)solver->row_refactor_checked_run_count;
  solver->stats.row_refactor_parallel_run_count =
    (int64_t)solver->row_refactor_parallel_run_count;
  solver->stats.row_refactor_segment_count =
    (int64_t)solver->row_refactor_segment_count;
  solver->stats.row_refactor_segment_rows =
    (int64_t)solver->row_refactor_segment_rows;
  solver->stats.row_refactor_segment_max_width =
    (int64_t)solver->row_refactor_segment_max_width;
  solver->stats.row_refactor_segment_dense_entries =
    solver->row_refactor_segment_dense_entries;
  solver->stats.row_refactor_segment_trailing_entries =
    solver->row_refactor_segment_trailing_entries;
  solver->stats.row_refactor_dense_segment_count =
    (int64_t)solver->row_refactor_dense_segment_count;
  solver->stats.row_refactor_dense_segment_rows =
    (int64_t)solver->row_refactor_dense_segment_rows;
  solver->stats.row_refactor_dense_segment_max_width =
    (int64_t)solver->row_refactor_dense_segment_max_width;
  solver->stats.row_refactor_dense_segment_dense_entries =
    solver->row_refactor_dense_segment_dense_entries;
  solver->stats.row_refactor_dense_segment_trailing_entries =
    solver->row_refactor_dense_segment_trailing_entries;
  solver->stats.refactor_dependency_cluster_levels =
    (int64_t)solver->refactor_cluster_level_count;
  solver->stats.refactor_dependency_pipeline_columns =
    (int64_t)solver->refactor_pipeline_column_count;
  solver->stats.refactor_dependency_work =
    solver->refactor_dependency_work;
  solver->stats.refactor_dependency_pipeline_work =
    solver->refactor_pipeline_work;
}

static double *ensure_solve_perm_workspace(kls_solver *solver) {
  if (solver == NULL || solver->n == 0) {
    return NULL;
  }
  if (solver->solve_perm_workspace != NULL &&
      solver->solve_perm_workspace_n == solver->n) {
    return solver->solve_perm_workspace;
  }
  free(solver->solve_perm_workspace);
  solver->solve_perm_workspace = NULL;
  solver->solve_perm_workspace_n = 0;
  solver->solve_perm_workspace =
    (double *)malloc((size_t)solver->n * sizeof(*solver->solve_perm_workspace));
  if (solver->solve_perm_workspace == NULL) {
    return NULL;
  }
  solver->solve_perm_workspace_n = solver->n;
  return solver->solve_perm_workspace;
}

void kls_default_options(kls_options *options) {
  if (options == NULL) {
    return;
  }
  memset(options, 0, sizeof(*options));
  options->struct_size = sizeof(*options);
  options->threads = 1;
  options->ordering = KLS_ORDERING_AUTO;
  options->orientation = KLS_ORIENTATION_AUTO;
  options->use_btf = 1;
  options->scale = KLS_SCALE_AUTO;
  options->pivot_tolerance = 0.001;
  options->memory_growth = 1.5;
  options->halt_if_singular = 1;
  options->fast_factor = 1;
  options->static_pivoting = 1;
}

int kls_create(kls_solver **solver_out) {
  if (solver_out == NULL) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  kls_solver *solver = (kls_solver *)calloc(1, sizeof(*solver));
  if (solver == NULL) {
    return KLS_ERR_OUT_OF_MEMORY;
  }
  kls_default_options(&solver->options);
  solver->stats.struct_size = sizeof(solver->stats);
  kls_clear_fast_reject_stats(solver);
  if (!trilinos_klu_l_defaults(&solver->common)) {
    free(solver);
    return KLS_ERR_ANALYZE_FAILED;
  }
  *solver_out = solver;
  return KLS_OK;
}

void kls_destroy(kls_solver *solver) {
  if (solver == NULL) {
    return;
  }
  clear_matrix(solver);
  free(solver);
}

int kls_analyze_csc(kls_solver *solver,
                    kls_index_type index_type,
                    int64_t n,
                    const void *col_ptr,
                    const void *row_idx,
                    int index_base,
                    const kls_options *options) {
  if (solver == NULL || col_ptr == NULL || row_idx == NULL ||
      !validate_index_base(index_base) || !validate_n(n)) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  if (index_type != KLS_INDEX_INT32 && index_type != KLS_INDEX_INT64) {
    return KLS_ERR_INVALID_ARGUMENT;
  }

  kls_options normalized;
  kls_default_options(&normalized);
  if (options != NULL) {
    if (options->struct_size < sizeof(kls_options)) {
      return KLS_ERR_INVALID_ARGUMENT;
    }
    normalized = *options;
  }
  if (!validate_options(&normalized)) {
    return KLS_ERR_INVALID_ARGUMENT;
  }

  clear_matrix(solver);
  solver->options = normalized;
  solver->input_format = KLS_INPUT_CSC;
  const double start = kls_now_seconds();

  kls_pattern_candidate normal = {0};
  kls_pattern_candidate transpose = {0};
  kls_pattern_candidate *chosen = NULL;

  int status = copy_compressed_candidate(&normal, index_type, n, col_ptr, row_idx,
                                         index_base, KLS_ORIENTATION_NORMAL);
  if (status != KLS_OK) {
    clear_matrix(solver);
    return status;
  }

  const int prefer_auto_normal =
    normalized.orientation == KLS_ORIENTATION_AUTO &&
    auto_orientation_prefers_normal(normal.n, normal.col_ptr, normal.row_idx);
  if (normalized.orientation != KLS_ORIENTATION_NORMAL && !prefer_auto_normal) {
    status = transpose_candidate(&normal, KLS_ORIENTATION_TRANSPOSE, &transpose);
    if (status != KLS_OK) {
      free_candidate(&normal);
      clear_matrix(solver);
      return status;
    }
  }

  status = select_candidate(&normal,
                            (normalized.orientation == KLS_ORIENTATION_NORMAL ||
                             prefer_auto_normal) ? NULL : &transpose,
                            &normalized, &chosen);
  const double elapsed = kls_now_seconds() - start;
  if (status != KLS_OK) {
    free_candidate(&transpose);
    free_candidate(&normal);
    clear_matrix(solver);
    return status;
  }

  adopt_candidate(solver, chosen);
  fill_symbolic_stats(solver, elapsed);
  free_candidate(&transpose);
  free_candidate(&normal);
  return KLS_OK;
}

int kls_analyze_csr(kls_solver *solver,
                    kls_index_type index_type,
                    int64_t n,
                    const void *row_ptr,
                    const void *col_idx,
                    int index_base,
                    const kls_options *options) {
  if (solver == NULL || row_ptr == NULL || col_idx == NULL ||
      !validate_index_base(index_base) || !validate_n(n)) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  if (index_type != KLS_INDEX_INT32 && index_type != KLS_INDEX_INT64) {
    return KLS_ERR_INVALID_ARGUMENT;
  }

  kls_options normalized;
  kls_default_options(&normalized);
  if (options != NULL) {
    if (options->struct_size < sizeof(kls_options)) {
      return KLS_ERR_INVALID_ARGUMENT;
    }
    normalized = *options;
  }
  if (!validate_options(&normalized)) {
    return KLS_ERR_INVALID_ARGUMENT;
  }

  clear_matrix(solver);
  solver->options = normalized;
  solver->input_format = KLS_INPUT_CSR;
  const double start = kls_now_seconds();

  kls_pattern_candidate transpose = {0};
  kls_pattern_candidate normal = {0};
  kls_pattern_candidate *chosen = NULL;

  int status = copy_compressed_candidate(&transpose, index_type, n, row_ptr, col_idx,
                                         index_base, KLS_ORIENTATION_TRANSPOSE);
  if (status != KLS_OK) {
    clear_matrix(solver);
    return status;
  }

  const int prefer_auto_transpose =
    normalized.orientation == KLS_ORIENTATION_AUTO &&
    auto_orientation_prefers_transpose((UF_long)n);
  if (normalized.orientation != KLS_ORIENTATION_TRANSPOSE && !prefer_auto_transpose) {
    status = transpose_candidate(&transpose, KLS_ORIENTATION_NORMAL, &normal);
    if (status != KLS_OK) {
      free_candidate(&transpose);
      clear_matrix(solver);
      return status;
    }
  }

  const int prefer_auto_normal =
    normalized.orientation == KLS_ORIENTATION_AUTO && normal.col_ptr != NULL &&
    auto_orientation_prefers_normal(normal.n, normal.col_ptr, normal.row_idx);
  status = select_candidate(normal.col_ptr == NULL ? NULL : &normal,
                            prefer_auto_normal ? NULL : &transpose,
                            &normalized, &chosen);
  const double elapsed = kls_now_seconds() - start;
  if (status != KLS_OK) {
    free_candidate(&normal);
    free_candidate(&transpose);
    clear_matrix(solver);
    return status;
  }

  adopt_candidate(solver, chosen);
  fill_symbolic_stats(solver, elapsed);
  free_candidate(&normal);
  free_candidate(&transpose);
  return KLS_OK;
}

static int prepare_numeric_values(kls_solver *solver, const double *values, double **values_out) {
  if (solver == NULL || values == NULL || values_out == NULL) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  if (solver->input_format == KLS_INPUT_NONE) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  if (solver->input_to_csc == NULL &&
      solver->row_scale == NULL && solver->col_scale == NULL) {
    *values_out = (double *)values;
    return KLS_OK;
  }
  if (solver->values == NULL) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  if (solver->input_to_csc == NULL) {
    memcpy(solver->values, values, (size_t)solver->nnz * sizeof(*solver->values));
  } else {
    for (int64_t p = 0; p < solver->nnz; ++p) {
      solver->values[solver->input_to_csc[p]] = values[p];
    }
  }
  apply_value_scaling(solver->n, solver->col_ptr, solver->row_idx,
                      solver->row_scale, solver->col_scale, solver->values);
  *values_out = solver->values;
  return KLS_OK;
}

static int kls_parallel_refactor_is_eligible(const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL) {
    return 0;
  }
  if (solver->options.threads <= 1 || solver->symbolic->nblocks < 8u) {
    return 0;
  }
  if (solver->n < 5000u || solver->symbolic->maxblock == solver->n) {
    return 0;
  }
  if (solver->common.flops < 2.0e7) {
    return 0;
  }
  if (solver->symbolic->nblocks < 64u ||
      solver->symbolic->maxblock * 4u < solver->n * 3u) {
    return 0;
  }
  if (solver->symbolic->nblocks >= 1024u &&
      solver->symbolic->maxblock * 20u >= solver->n * 19u &&
      solver->symbolic->maxblock < 90000u &&
      solver->common.flops < 1.0e9) {
    return 0;
  }
  if (solver->numeric->Offp == NULL || solver->numeric->Offx == NULL ||
      solver->numeric->Pinv == NULL || solver->numeric->Udiag == NULL) {
    return 0;
  }
  if (solver->common.scale > 0) {
    const UF_long independent_rows = solver->n - solver->symbolic->maxblock;
    /* Scaled refactor pays O(n) scale recomputation and Rs permutation before
       any threaded block work, so require enough non-dominant BTF rows. */
    if (independent_rows < 2048u ||
        (double)independent_rows < 0.02 * (double)solver->n) {
      return 0;
    }
    return solver->numeric->Rs != NULL && solver->numeric->Pnum != NULL &&
           solver->numeric->Xwork != NULL;
  }
  if (solver->numeric->Rs != NULL) {
    return 0;
  }
  return 1;
}

static int kls_refactor_pool_map_is_worthwhile(const kls_solver *solver) {
  if (!kls_parallel_refactor_is_eligible(solver) ||
      solver->symbolic == NULL || solver->symbolic->nblocks <= 1u) {
    return 0;
  }
  if ((double)solver->symbolic->maxblock < 0.75 * (double)solver->n) {
    return 0;
  }
  if (solver->symbolic->nblocks > 20000u || solver->nnz > 3000000u) {
    return 0;
  }
  return 1;
}

static int kls_scaled_serial_mapped_btf_is_worthwhile(
  const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL || solver->n == 0u ||
      solver->common.scale <= 0 || solver->symbolic->nblocks < 1024u ||
      solver->symbolic->nblocks == 1u) {
    return 0;
  }
  const double coverage =
    (double)solver->symbolic->maxblock / (double)solver->n;
  return coverage >= 0.95 &&
         solver->symbolic->maxblock < 90000u &&
         solver->common.flops >= 1.0e8;
}

static int kls_parallel_refactor_permute_scale(kls_solver *solver) {
  trilinos_klu_l_numeric *numeric = solver->numeric;
  double *rs = numeric->Rs;
  double *xwork = (double *)numeric->Xwork;
  if (rs == NULL || xwork == NULL || numeric->Pnum == NULL) {
    return 0;
  }
  for (UF_long k = 0; k < solver->n; ++k) {
    const UF_long row = numeric->Pnum[k];
    if (row >= solver->n) {
      return 0;
    }
    xwork[k] = rs[row];
  }
  for (UF_long k = 0; k < solver->n; ++k) {
    rs[k] = xwork[k];
  }
  return 1;
}

static int kls_numeric_pivots_pass_threshold_from_block(
  const kls_solver *solver,
  UF_long start_block,
  UF_long *rejected_pivot_out,
  UF_long *rejected_pivot_col_out) {
  if (rejected_pivot_out != NULL) {
    *rejected_pivot_out = KLS_KLU_EMPTY;
  }
  if (rejected_pivot_col_out != NULL) {
    *rejected_pivot_col_out = KLS_KLU_EMPTY;
  }
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL) {
    return 0;
  }
  const double tol = solver->common.tol;
  if (tol <= DBL_MIN) {
    return 1;
  }
  const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  const trilinos_klu_l_numeric *numeric = solver->numeric;
  if (numeric->Lip == NULL || numeric->Llen == NULL || numeric->LUbx == NULL) {
    return 0;
  }
  if (start_block > symbolic->nblocks) {
    return 0;
  }

  for (UF_long block = start_block; block < symbolic->nblocks; ++block) {
    const UF_long k1 = symbolic->R[block];
    const UF_long k2 = symbolic->R[block + 1u];
    const UF_long nk = k2 - k1;
    if (nk <= 1u) {
      continue;
    }
    double *lu = (double *)numeric->LUbx[block];
    if (lu == NULL) {
      return 0;
    }
    const UF_long *lip = numeric->Lip + k1;
    const UF_long *llen = numeric->Llen + k1;
    for (UF_long k = 0; k < nk; ++k) {
      UF_long *li = NULL;
      double *lx = NULL;
      UF_long lcol_len = 0;
      kls_klu_get_pointer(lu, lip, llen, k, &li, &lx, &lcol_len);
      (void)li;
      for (UF_long p = 0; p < lcol_len; ++p) {
        const double value_abs = fabs(lx[p]);
        if (!isfinite(value_abs) || value_abs * tol > 1.0 + 1.0e-12) {
          const UF_long global_col = k1 + k;
          if (rejected_pivot_out != NULL) {
            *rejected_pivot_out = global_col;
          }
          if (rejected_pivot_col_out != NULL) {
            *rejected_pivot_col_out =
              symbolic->Q != NULL ? symbolic->Q[global_col] : KLS_KLU_EMPTY;
          }
          return 0;
        }
      }
    }
  }
  return 1;
}

static int kls_numeric_pivots_pass_threshold(const kls_solver *solver,
                                             UF_long *rejected_pivot_out,
                                             UF_long *rejected_pivot_col_out) {
  return kls_numeric_pivots_pass_threshold_from_block(
    solver, 0u, rejected_pivot_out, rejected_pivot_col_out);
}

static UF_long kls_block_for_pivot(const kls_solver *solver, UF_long pivot) {
  if (solver == NULL || solver->symbolic == NULL ||
      solver->symbolic->R == NULL || pivot >= solver->n) {
    return KLS_KLU_EMPTY;
  }
  UF_long lo = 0;
  UF_long hi = solver->symbolic->nblocks;
  while (lo + 1u < hi) {
    const UF_long mid = lo + (hi - lo) / 2u;
    if (solver->symbolic->R[mid] <= pivot) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  return pivot >= solver->symbolic->R[lo] &&
         pivot < solver->symbolic->R[lo + 1u] ? lo : KLS_KLU_EMPTY;
}

static int kls_fast_reject_column_work(const kls_solver *solver,
                                       UF_long block,
                                       UF_long k1,
                                       UF_long nk,
                                       UF_long local_col,
                                       double *work_out) {
  if (solver == NULL || solver->numeric == NULL || work_out == NULL ||
      block >= solver->numeric->nblocks || local_col >= nk ||
      solver->numeric->Uip == NULL || solver->numeric->Ulen == NULL ||
      solver->numeric->Llen == NULL || solver->numeric->LUbx == NULL ||
      solver->numeric->LUbx[block] == NULL) {
    return 0;
  }
  double *lu = (double *)solver->numeric->LUbx[block];
  const UF_long *uip = solver->numeric->Uip + k1;
  const UF_long *ulen = solver->numeric->Ulen + k1;
  UF_long *ui = NULL;
  double *ux = NULL;
  UF_long ucol_len = 0;
  kls_klu_get_pointer(lu, uip, ulen, local_col, &ui, &ux, &ucol_len);
  (void)ux;

  double work = 1.0;
  for (UF_long p = 0; p < ucol_len; ++p) {
    const UF_long dep = ui[p];
    if (dep >= local_col || dep >= nk) {
      return 0;
    }
    work += 1.0 + (double)solver->numeric->Llen[k1 + dep];
  }
  *work_out = work;
  return 1;
}

static int kls_fast_reject_block_work(const kls_solver *solver,
                                      UF_long block,
                                      UF_long k1,
                                      UF_long nk,
                                      double *work_out) {
  if (solver == NULL || work_out == NULL) {
    return 0;
  }
  double work = 0.0;
  for (UF_long k = 0; k < nk; ++k) {
    double column_work = 0.0;
    if (!kls_fast_reject_column_work(solver, block, k1, nk, k,
                                     &column_work)) {
      return 0;
    }
    work += column_work;
  }
  *work_out = work;
  return 1;
}

static int kls_fast_reject_suffix_work(const kls_solver *solver,
                                       UF_long block,
                                       UF_long k1,
                                       UF_long nk,
                                       UF_long local_reject,
                                       UF_long *columns_out,
                                       double *work_out) {
  if (solver == NULL || columns_out == NULL || work_out == NULL ||
      local_reject >= nk) {
    return 0;
  }
  double work = 0.0;
  for (UF_long k = local_reject; k < nk; ++k) {
    double column_work = 0.0;
    if (!kls_fast_reject_column_work(solver, block, k1, nk, k,
                                     &column_work)) {
      return 0;
    }
    work += column_work;
  }
  *columns_out = nk - local_reject;
  *work_out = work;
  return 1;
}

static int kls_build_ordered_block_etree(const kls_solver *solver,
                                         UF_long k1,
                                         UF_long k2,
                                         UF_long *parent) {
  if (solver == NULL || solver->symbolic == NULL || parent == NULL ||
      solver->col_ptr == NULL || solver->row_idx == NULL ||
      solver->symbolic->P == NULL || solver->symbolic->Q == NULL ||
      k1 > k2 || k2 > solver->n) {
    return 0;
  }
  const UF_long nk = k2 - k1;
  UF_long *pinv = (UF_long *)malloc((size_t)solver->n * sizeof(*pinv));
  UF_long *ancestor = (UF_long *)malloc((size_t)nk * sizeof(*ancestor));
  UF_long *prev = (UF_long *)malloc((size_t)nk * sizeof(*prev));
  if (pinv == NULL || ancestor == NULL || prev == NULL) {
    free(pinv);
    free(ancestor);
    free(prev);
    return 0;
  }

  for (UF_long k = 0; k < solver->n; ++k) {
    pinv[k] = KLS_KLU_EMPTY;
  }
  for (UF_long k = 0; k < solver->n; ++k) {
    const UF_long row = solver->symbolic->P[k];
    if (row >= solver->n) {
      free(pinv);
      free(ancestor);
      free(prev);
      return 0;
    }
    pinv[row] = k;
  }
  for (UF_long k = 0; k < nk; ++k) {
    parent[k] = KLS_KLU_EMPTY;
    ancestor[k] = KLS_KLU_EMPTY;
    prev[k] = KLS_KLU_EMPTY;
  }

  for (UF_long local_col = 0; local_col < nk; ++local_col) {
    const UF_long old_col = solver->symbolic->Q[k1 + local_col];
    if (old_col >= solver->n) {
      free(pinv);
      free(ancestor);
      free(prev);
      return 0;
    }
    for (UF_long p = solver->col_ptr[old_col];
         p < solver->col_ptr[old_col + 1u]; ++p) {
      const UF_long old_row = solver->row_idx[p];
      if (old_row >= solver->n) {
        free(pinv);
        free(ancestor);
        free(prev);
        return 0;
      }
      const UF_long ordered_row = pinv[old_row];
      if (ordered_row < k1 || ordered_row >= k2) {
        continue;
      }
      const UF_long local_row = ordered_row - k1;
      UF_long j = prev[local_row];
      while (j != KLS_KLU_EMPTY && j < local_col) {
        const UF_long next = ancestor[j];
        ancestor[j] = local_col;
        if (parent[j] == KLS_KLU_EMPTY) {
          parent[j] = local_col;
        }
        j = next;
      }
      prev[local_row] = local_col;
    }
  }

  free(pinv);
  free(ancestor);
  free(prev);
  return 1;
}

static int kls_ensure_row_refactor_etree_parent(kls_solver *solver) {
  if (solver == NULL || solver->row_refactor_pattern_n != solver->n ||
      solver->symbolic == NULL || solver->symbolic->nblocks != 1u) {
    return 0;
  }
  if (solver->row_refactor_etree_parent != NULL || solver->n == 0u) {
    return 1;
  }
  UF_long *parent =
    (UF_long *)malloc((size_t)solver->n * sizeof(*parent));
  if (parent == NULL) {
    return 0;
  }
  if (!kls_build_ordered_block_etree(solver, 0u, solver->n, parent)) {
    free(parent);
    return 0;
  }
  solver->row_refactor_etree_parent = parent;
  return 1;
}

static int kls_get_ordered_block_etree_parent(kls_solver *solver,
                                              UF_long block,
                                              UF_long k1,
                                              UF_long k2,
                                              const UF_long **parent_out,
                                              UF_long **owned_parent_out) {
  if (parent_out != NULL) {
    *parent_out = NULL;
  }
  if (owned_parent_out != NULL) {
    *owned_parent_out = NULL;
  }
  if (solver == NULL || solver->symbolic == NULL || parent_out == NULL ||
      owned_parent_out == NULL || k1 > k2 || k2 > solver->n) {
    return 0;
  }

  if (block == 0u && k1 == 0u && k2 == solver->n &&
      kls_ensure_row_refactor_etree_parent(solver)) {
    *parent_out = solver->row_refactor_etree_parent;
    return *parent_out != NULL || solver->n == 0u;
  }

  const UF_long nk = k2 - k1;
  UF_long *owned_parent =
    (UF_long *)malloc((size_t)nk * sizeof(*owned_parent));
  if (owned_parent == NULL) {
    return 0;
  }
  if (!kls_build_ordered_block_etree(solver, k1, k2, owned_parent)) {
    free(owned_parent);
    return 0;
  }
  *parent_out = owned_parent;
  *owned_parent_out = owned_parent;
  return 1;
}

static int kls_ensure_fast_reject_tail_plan_storage(kls_solver *solver,
                                                    UF_long capacity) {
  if (solver == NULL) {
    return 0;
  }
  if (capacity == 0u) {
    return 1;
  }
  if (solver->fast_reject_tail_capacity >= capacity &&
      solver->fast_reject_tail_cols != NULL &&
      solver->fast_reject_tail_marks != NULL) {
    return 1;
  }
  if (capacity > (UF_long)(SIZE_MAX / sizeof(*solver->fast_reject_tail_cols)) ||
      capacity > (UF_long)(SIZE_MAX / sizeof(*solver->fast_reject_tail_marks))) {
    return 0;
  }

  UF_long *cols =
    (UF_long *)malloc((size_t)capacity *
                      sizeof(*solver->fast_reject_tail_cols));
  unsigned int *marks =
    (unsigned int *)calloc((size_t)capacity,
                           sizeof(*solver->fast_reject_tail_marks));
  if (cols == NULL || marks == NULL) {
    free(cols);
    free(marks);
    return 0;
  }
  free(solver->fast_reject_tail_cols);
  free(solver->fast_reject_tail_marks);
  solver->fast_reject_tail_cols = cols;
  solver->fast_reject_tail_marks = marks;
  solver->fast_reject_tail_capacity = capacity;
  solver->fast_reject_tail_mark = 0u;
  solver->fast_reject_tail_count = 0;
  return 1;
}

static int kls_mark_fast_reject_tail_col(kls_solver *solver,
                                         UF_long nk,
                                         UF_long local_col,
                                         unsigned int mark,
                                         UF_long *tail_count) {
  if (solver == NULL || local_col >= nk || tail_count == NULL ||
      solver->fast_reject_tail_cols == NULL ||
      solver->fast_reject_tail_marks == NULL) {
    return 0;
  }
  if (solver->fast_reject_tail_marks[local_col] == mark) {
    return 1;
  }
  solver->fast_reject_tail_marks[local_col] = mark;
  solver->fast_reject_tail_cols[(*tail_count)++] = local_col;
  return 1;
}

static void kls_record_fast_reject_pivoting_tail_plan(
  kls_solver *solver,
  UF_long k1,
  UF_long nk,
  UF_long local_reject,
  const UF_long *parent,
  unsigned int mark,
  UF_long columns,
  double work) {
  if (solver == NULL) {
    return;
  }
  solver->stats.fast_rejected_pivoting_tail_columns = (int64_t)columns;
  solver->stats.fast_rejected_pivoting_tail_work = work;
  solver->stats.fast_rejected_pivoting_tail_first = -1;
  solver->stats.fast_rejected_pivoting_tail_last = -1;
  solver->stats.fast_rejected_pivoting_tail_contains_reject = 0;
  solver->stats.fast_rejected_pivoting_tail_topological = 0;
  if (columns == 0u || solver->fast_reject_tail_cols == NULL ||
      solver->fast_reject_tail_marks == NULL || parent == NULL ||
      local_reject >= nk) {
    return;
  }

  UF_long previous = KLS_KLU_EMPTY;
  int sorted = 1;
  int contains_reject = 0;
  int topological = 1;
  for (UF_long p = 0; p < columns; ++p) {
    const UF_long global_col = solver->fast_reject_tail_cols[p];
    if (global_col < k1 || global_col >= k1 + nk) {
      sorted = 0;
      topological = 0;
      continue;
    }
    const UF_long local_col = global_col - k1;
    if (p == 0u) {
      solver->stats.fast_rejected_pivoting_tail_first =
        (int64_t)global_col;
    } else if (previous >= global_col) {
      sorted = 0;
      topological = 0;
    }
    previous = global_col;
    if (local_col == local_reject) {
      contains_reject = 1;
    }
    const UF_long next = parent[local_col];
    if (next == KLS_KLU_EMPTY || next == local_col) {
      continue;
    }
    if (next >= nk) {
      topological = 0;
      continue;
    }
    if (solver->fast_reject_tail_marks[next] == mark && next <= local_col) {
      topological = 0;
    }
  }
  solver->stats.fast_rejected_pivoting_tail_last = (int64_t)previous;
  solver->stats.fast_rejected_pivoting_tail_contains_reject =
    contains_reject;
  solver->stats.fast_rejected_pivoting_tail_topological =
    sorted && contains_reject && topological;
}

static int kls_build_fast_reject_pivoting_tail_plan(
  kls_solver *solver,
  UF_long block,
  UF_long k1,
  UF_long k2,
  UF_long local_reject,
  const UF_long **tail_cols_out,
  UF_long *tail_count_out,
  double *tail_work_out) {
  if (tail_cols_out != NULL) {
    *tail_cols_out = NULL;
  }
  if (tail_count_out != NULL) {
    *tail_count_out = 0;
  }
  if (tail_work_out != NULL) {
    *tail_work_out = 0.0;
  }
  if (solver == NULL || tail_cols_out == NULL || tail_count_out == NULL ||
      tail_work_out == NULL || local_reject >= k2 - k1) {
    return 0;
  }

  const UF_long nk = k2 - k1;
  if (!kls_ensure_fast_reject_tail_plan_storage(solver, nk)) {
    return 0;
  }
  if (solver->fast_reject_tail_mark == UINT_MAX) {
    memset(solver->fast_reject_tail_marks, 0,
           (size_t)solver->fast_reject_tail_capacity *
             sizeof(*solver->fast_reject_tail_marks));
    solver->fast_reject_tail_mark = 0u;
  }

  const UF_long *parent = NULL;
  UF_long *owned_parent = NULL;
  if (!kls_get_ordered_block_etree_parent(solver, block, k1, k2, &parent,
                                          &owned_parent) ||
      parent == NULL) {
    free(owned_parent);
    return 0;
  }

  const unsigned int mark = ++solver->fast_reject_tail_mark;
  UF_long queued = 0;
  if (solver->fast_reject_refresh_state == KLS_FAST_REJECT_REFRESH_PREFIX) {
    for (UF_long local_col = local_reject; local_col < nk; ++local_col) {
      if (!kls_mark_fast_reject_tail_col(solver, nk, local_col, mark,
                                         &queued)) {
        free(owned_parent);
        return 0;
      }
    }
  } else if (!kls_mark_fast_reject_tail_col(solver, nk, local_reject, mark,
                                            &queued)) {
    free(owned_parent);
    return 0;
  }

  UF_long head = 0;
  while (head < queued) {
    const UF_long local_col = solver->fast_reject_tail_cols[head++];
    const UF_long next = parent[local_col];
    if (next == KLS_KLU_EMPTY || next == local_col) {
      continue;
    }
    if (next >= nk ||
        !kls_mark_fast_reject_tail_col(solver, nk, next, mark, &queued)) {
      free(owned_parent);
      return 0;
    }
  }

  UF_long columns = 0;
  double work = 0.0;
  for (UF_long local_col = 0; local_col < nk; ++local_col) {
    if (solver->fast_reject_tail_marks[local_col] != mark) {
      continue;
    }
    double column_work = 0.0;
    if (!kls_fast_reject_column_work(solver, block, k1, nk, local_col,
                                     &column_work)) {
      free(owned_parent);
      return 0;
    }
    solver->fast_reject_tail_cols[columns++] = k1 + local_col;
    work += column_work;
  }

  solver->fast_reject_tail_count = columns;
  kls_record_fast_reject_pivoting_tail_plan(
    solver, k1, nk, local_reject, parent, mark, columns, work);
  *tail_cols_out = solver->fast_reject_tail_cols;
  *tail_count_out = columns;
  *tail_work_out = work;
  free(owned_parent);
  return 1;
}

static double kls_row_refactor_tail_row_work(const kls_solver *solver,
                                             UF_long row) {
  if (solver == NULL || row >= solver->n ||
      solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_u_ptr == NULL ||
      solver->row_refactor_input_ptr == NULL ||
      solver->row_refactor_l_cols == NULL) {
    return 1.0;
  }
  double work = 1.0;
  work += (double)(solver->row_refactor_input_ptr[row + 1u] -
                   solver->row_refactor_input_ptr[row]);
  work += (double)(solver->row_refactor_u_ptr[row + 1u] -
                   solver->row_refactor_u_ptr[row]);
  for (UF_long p = solver->row_refactor_l_ptr[row];
       p < solver->row_refactor_l_ptr[row + 1u]; ++p) {
    const UF_long dep = solver->row_refactor_l_cols[p];
    if (dep < solver->n) {
      work += 1.0 + (double)(solver->row_refactor_u_ptr[dep + 1u] -
                             solver->row_refactor_u_ptr[dep]);
    }
  }
  return work;
}

static int kls_build_row_refactor_tail(kls_solver *solver,
                                       UF_long root,
                                       const UF_long **rows_out,
                                       UF_long *count_out) {
  if (rows_out != NULL) {
    *rows_out = NULL;
  }
  if (count_out != NULL) {
    *count_out = 0;
  }
  if (solver == NULL || root >= solver->n ||
      solver->row_refactor_pattern_n != solver->n ||
      solver->row_refactor_successor_ptr == NULL ||
      (solver->row_refactor_successor_ptr[solver->n] > 0u &&
       solver->row_refactor_successor_rows == NULL) ||
      rows_out == NULL || count_out == NULL) {
    return 0;
  }

  const UF_long n = solver->n;
  if (solver->row_refactor_tail_rows == NULL) {
    solver->row_refactor_tail_rows =
      (UF_long *)malloc((size_t)n * sizeof(*solver->row_refactor_tail_rows));
  }
  if (solver->row_refactor_tail_marks == NULL) {
    solver->row_refactor_tail_marks =
      (unsigned int *)calloc((size_t)n,
                             sizeof(*solver->row_refactor_tail_marks));
  }
  if (solver->row_refactor_tail_rows == NULL ||
      solver->row_refactor_tail_marks == NULL) {
    return 0;
  }

  if (solver->row_refactor_tail_mark == UINT_MAX) {
    memset(solver->row_refactor_tail_marks, 0,
           (size_t)n * sizeof(*solver->row_refactor_tail_marks));
    solver->row_refactor_tail_mark = 0u;
  }
  const unsigned int mark = ++solver->row_refactor_tail_mark;
  unsigned int *marks = solver->row_refactor_tail_marks;
  UF_long *queue = solver->row_refactor_tail_rows;

  UF_long head = 0;
  UF_long tail = 0;
  marks[root] = mark;
  queue[tail++] = root;
  while (head < tail) {
    const UF_long row = queue[head++];
    for (UF_long p = solver->row_refactor_successor_ptr[row];
         p < solver->row_refactor_successor_ptr[row + 1u]; ++p) {
      const UF_long successor = solver->row_refactor_successor_rows[p];
      if (successor >= n) {
        solver->row_refactor_tail_count = 0;
        return 0;
      }
      if (marks[successor] != mark) {
        marks[successor] = mark;
        queue[tail++] = successor;
      }
    }
  }

  UF_long count = 0;
  for (UF_long row = root; row < n; ++row) {
    if (marks[row] == mark) {
      solver->row_refactor_tail_rows[count++] = row;
    }
  }
  solver->row_refactor_tail_count = count;
  *rows_out = solver->row_refactor_tail_rows;
  *count_out = count;
  return 1;
}

static void kls_fill_fast_reject_row_tail_stats(kls_solver *solver,
                                                UF_long block,
                                                UF_long k1,
                                                UF_long k2,
                                                UF_long local_reject) {
  if (solver == NULL || local_reject >= k2 - k1) {
    return;
  }
  solver->row_refactor_tail_count = 0;

  if (block == 0u && k1 == 0u && k2 == solver->n) {
    const UF_long *tail_rows = NULL;
    UF_long columns = 0;
    if (kls_build_row_refactor_tail(solver, local_reject,
                                    &tail_rows, &columns)) {
      double work = 0.0;
      for (UF_long p = 0; p < columns; ++p) {
        work += kls_row_refactor_tail_row_work(solver, tail_rows[p]);
      }

      solver->stats.fast_rejected_row_tail_columns = (int64_t)columns;
      solver->stats.fast_rejected_row_tail_work = work;
      return;
    }
  }

  if (solver->numeric == NULL || block >= solver->numeric->nblocks ||
      solver->numeric->LUbx == NULL || solver->numeric->LUbx[block] == NULL ||
      solver->numeric->Lip == NULL || solver->numeric->Llen == NULL) {
    return;
  }

  const UF_long nk = k2 - k1;
  if (nk == 0u || k2 > solver->n) {
    return;
  }
  if (solver->row_refactor_tail_rows == NULL) {
    solver->row_refactor_tail_rows =
      (UF_long *)malloc((size_t)solver->n *
                        sizeof(*solver->row_refactor_tail_rows));
  }
  if (solver->row_refactor_tail_marks == NULL) {
    solver->row_refactor_tail_marks =
      (unsigned int *)calloc((size_t)solver->n,
                             sizeof(*solver->row_refactor_tail_marks));
  }
  if (solver->row_refactor_tail_rows == NULL ||
      solver->row_refactor_tail_marks == NULL) {
    return;
  }
  if (solver->row_refactor_tail_mark == UINT_MAX) {
    memset(solver->row_refactor_tail_marks, 0,
           (size_t)solver->n * sizeof(*solver->row_refactor_tail_marks));
    solver->row_refactor_tail_mark = 0u;
  }

  const unsigned int mark = ++solver->row_refactor_tail_mark;
  unsigned int *marks = solver->row_refactor_tail_marks;
  UF_long *queue = solver->row_refactor_tail_rows;
  double *lu = (double *)solver->numeric->LUbx[block];
  const UF_long *lip = solver->numeric->Lip + k1;
  const UF_long *llen = solver->numeric->Llen + k1;

  UF_long head = 0;
  UF_long tail = 0;
  marks[k1 + local_reject] = mark;
  queue[tail++] = local_reject;
  while (head < tail) {
    const UF_long local_col = queue[head++];
    UF_long *li = NULL;
    double *lx = NULL;
    UF_long lcol_len = 0;
    kls_klu_get_pointer(lu, lip, llen, local_col, &li, &lx, &lcol_len);
    (void)lx;
    for (UF_long p = 0; p < lcol_len; ++p) {
      const UF_long successor = li[p];
      if (successor <= local_col || successor >= nk) {
        return;
      }
      const UF_long global_successor = k1 + successor;
      if (marks[global_successor] != mark) {
        marks[global_successor] = mark;
        queue[tail++] = successor;
      }
    }
  }

  UF_long columns = 0;
  double work = 0.0;
  for (UF_long local_col = local_reject; local_col < nk; ++local_col) {
    if (marks[k1 + local_col] != mark) {
      continue;
    }
    double column_work = 0.0;
    if (!kls_fast_reject_column_work(solver, block, k1, nk, local_col,
                                     &column_work)) {
      return;
    }
    solver->row_refactor_tail_rows[columns] = k1 + local_col;
    columns++;
    work += column_work;
  }
  solver->row_refactor_tail_count = columns;
  solver->stats.fast_rejected_row_tail_columns = (int64_t)columns;
  solver->stats.fast_rejected_row_tail_work = work;
}

static void kls_fill_fast_reject_observed_tail_candidate(
  kls_solver *solver) {
  if (solver == NULL ||
      solver->stats.fast_rejected_tail_candidate_row >= 0 ||
      solver->stats.fast_rejected_row < 0 ||
      solver->stats.fast_rejected_candidate_abs < 0.0 ||
      solver->stats.fast_rejected_pivot_abs < 0.0 ||
      solver->row_refactor_tail_rows == NULL ||
      solver->row_refactor_tail_count == 0u) {
    return;
  }

  const UF_long row = (UF_long)solver->stats.fast_rejected_row;
  UF_long position = KLS_KLU_EMPTY;
  for (UF_long pos = 0; pos < solver->row_refactor_tail_count; ++pos) {
    if (solver->row_refactor_tail_rows[pos] == row) {
      position = pos;
      break;
    }
  }
  if (position == KLS_KLU_EMPTY) {
    return;
  }

  const double candidate_abs = solver->stats.fast_rejected_candidate_abs;
  const double pivot_abs = solver->stats.fast_rejected_pivot_abs;
  const double tol = solver->common.tol;
  const int tolerance_valid =
    tol > DBL_MIN && candidate_abs * tol > pivot_abs;

  solver->stats.fast_rejected_tail_candidate_row = (int64_t)row;
  solver->stats.fast_rejected_tail_candidate_abs = candidate_abs;
  solver->stats.fast_rejected_tail_candidate_count =
    tolerance_valid ? 1 : 0;
  solver->stats.fast_rejected_tail_candidate_position = (int64_t)position;
  solver->stats.fast_rejected_tail_repair_ready =
    solver->fast_reject_refresh_state == KLS_FAST_REJECT_REFRESH_PREFIX &&
    tolerance_valid;
}

static int kls_row_refactor_touch(UF_long col,
                                  unsigned char *touched,
                                  UF_long *touched_cols,
                                  UF_long *touched_count,
                                  UF_long n) {
  if (col >= n || touched == NULL || touched_cols == NULL ||
      touched_count == NULL) {
    return 0;
  }
  if (!touched[col]) {
    touched[col] = 1u;
    touched_cols[(*touched_count)++] = col;
  }
  return 1;
}

static int kls_row_refactor_prefix_candidate_value(
  const kls_solver *solver,
  const double *numeric_values,
  UF_long row,
  UF_long col,
  double *x,
  unsigned char *touched,
  UF_long *touched_cols,
  double *value_out) {
  if (solver == NULL || numeric_values == NULL || row >= solver->n ||
      col >= solver->n || x == NULL || touched == NULL ||
      touched_cols == NULL || value_out == NULL ||
      solver->row_refactor_input_ptr == NULL ||
      solver->row_refactor_input_cols == NULL ||
      solver->row_refactor_input_pos == NULL ||
      solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_l_cols == NULL ||
      solver->row_refactor_u_ptr == NULL ||
      solver->row_refactor_u_cols == NULL ||
      solver->row_refactor_u_row_values == NULL ||
      solver->numeric == NULL || solver->numeric->Udiag == NULL) {
    return 0;
  }

  UF_long touched_count = 0;
  int ok = 1;
  for (UF_long p = solver->row_refactor_input_ptr[row];
       ok && p < solver->row_refactor_input_ptr[row + 1u]; ++p) {
    const UF_long input_col = solver->row_refactor_input_cols[p];
    if (!kls_row_refactor_touch(input_col, touched, touched_cols,
                                &touched_count, solver->n)) {
      ok = 0;
      break;
    }
    x[input_col] = numeric_values[solver->row_refactor_input_pos[p]];
  }

  const double *udiag = (const double *)solver->numeric->Udiag;
  for (UF_long p = solver->row_refactor_l_ptr[row];
       ok && p < solver->row_refactor_l_ptr[row + 1u]; ++p) {
    const UF_long dep = solver->row_refactor_l_cols[p];
    if (dep >= col) {
      break;
    }
    if (dep >= solver->n || udiag[dep] == 0.0) {
      ok = 0;
      break;
    }
    if (!kls_row_refactor_touch(dep, touched, touched_cols,
                                &touched_count, solver->n)) {
      ok = 0;
      break;
    }
    const double lij = x[dep] / udiag[dep];
    x[dep] = 0.0;
    const UF_long u_begin = solver->row_refactor_u_ptr[dep];
    const UF_long u_end = solver->row_refactor_u_ptr[dep + 1u];
    const UF_long *u_cols = solver->row_refactor_u_cols + u_begin;
    const double *u_values =
      solver->row_refactor_u_row_values + u_begin;
    for (UF_long offset = 0; offset < u_end - u_begin; ++offset) {
      const UF_long target = u_cols[offset];
      if (target > col) {
        continue;
      }
      if (!kls_row_refactor_touch(target, touched, touched_cols,
                                  &touched_count, solver->n)) {
        ok = 0;
        break;
      }
      x[target] -= lij * u_values[offset];
    }
  }

  if (ok) {
    if (!kls_row_refactor_touch(col, touched, touched_cols,
                                &touched_count, solver->n)) {
      ok = 0;
    } else {
      *value_out = x[col];
    }
  }

  for (UF_long p = 0; p < touched_count; ++p) {
    const UF_long touched_col = touched_cols[p];
    x[touched_col] = 0.0;
    touched[touched_col] = 0u;
  }
  return ok;
}

static void kls_fill_fast_reject_row_tail_candidate_stats(
  kls_solver *solver,
  const double *numeric_values,
  UF_long rejected_pivot,
  double *x) {
  if (solver == NULL || numeric_values == NULL || x == NULL ||
      rejected_pivot >= solver->n ||
      solver->row_refactor_tail_rows == NULL ||
      solver->row_refactor_tail_count == 0u ||
      solver->stats.fast_rejected_pivot_abs < 0.0) {
    return;
  }

  unsigned char *touched =
    (unsigned char *)calloc((size_t)solver->n, sizeof(*touched));
  UF_long *touched_cols =
    (UF_long *)malloc((size_t)solver->n * sizeof(*touched_cols));
  if (touched == NULL || touched_cols == NULL) {
    free(touched);
    free(touched_cols);
    return;
  }

  UF_long best_row = KLS_KLU_EMPTY;
  UF_long best_pos = KLS_KLU_EMPTY;
  double best_abs = solver->stats.fast_rejected_pivot_abs;
  UF_long candidate_count = 0;
  const double pivot_abs = solver->stats.fast_rejected_pivot_abs;
  const double tol = solver->common.tol;
  for (UF_long pos = 0; pos < solver->row_refactor_tail_count; ++pos) {
    const UF_long row = solver->row_refactor_tail_rows[pos];
    if (row <= rejected_pivot || row >= solver->n) {
      continue;
    }
    double value = 0.0;
    if (!kls_row_refactor_prefix_candidate_value(
          solver, numeric_values, row, rejected_pivot, x, touched,
          touched_cols, &value)) {
      continue;
    }
    const double value_abs = fabs(value);
    if (tol > DBL_MIN && value_abs * tol > pivot_abs) {
      candidate_count++;
    }
    if (best_row == KLS_KLU_EMPTY || value_abs > best_abs) {
      best_row = row;
      best_pos = pos;
      best_abs = value_abs;
    }
  }

  if (best_row != KLS_KLU_EMPTY) {
    solver->stats.fast_rejected_tail_candidate_row = (int64_t)best_row;
    solver->stats.fast_rejected_tail_candidate_abs = best_abs;
    solver->stats.fast_rejected_tail_candidate_count =
      (int64_t)candidate_count;
    solver->stats.fast_rejected_tail_candidate_position =
      (int64_t)best_pos;
    solver->stats.fast_rejected_tail_repair_ready =
      solver->fast_reject_refresh_state == KLS_FAST_REJECT_REFRESH_PREFIX &&
      candidate_count > 0u && tol > DBL_MIN &&
      best_abs * tol > pivot_abs;
  }
  free(touched);
  free(touched_cols);
}

static void kls_fill_fast_reject_etree_tail_stats(kls_solver *solver,
                                                  UF_long block,
                                                  UF_long k1,
                                                  UF_long k2,
                                                  UF_long local_reject) {
  if (solver == NULL || local_reject >= k2 - k1) {
    return;
  }
  const UF_long nk = k2 - k1;
  const UF_long *parent = NULL;
  UF_long *owned_parent = NULL;
  if (!kls_get_ordered_block_etree_parent(solver, block, k1, k2, &parent,
                                          &owned_parent) ||
      parent == NULL) {
    free(owned_parent);
    return;
  }

  UF_long columns = 0;
  double work = 0.0;
  UF_long col = local_reject;
  while (col != KLS_KLU_EMPTY && col < nk) {
    double column_work = 0.0;
    if (!kls_fast_reject_column_work(solver, block, k1, nk, col,
                                     &column_work)) {
      free(owned_parent);
      return;
    }
    columns++;
    work += column_work;
    const UF_long next = parent[col];
    if (next == col) {
      break;
    }
    col = next;
  }

  solver->stats.fast_rejected_etree_columns = (int64_t)columns;
  solver->stats.fast_rejected_etree_work = work;
  free(owned_parent);

  const UF_long *tail_cols = NULL;
  UF_long tail_columns = 0;
  double tail_work = 0.0;
  if (kls_build_fast_reject_pivoting_tail_plan(
        solver, block, k1, k2, local_reject, &tail_cols, &tail_columns,
        &tail_work)) {
    (void)tail_cols;
    solver->stats.fast_rejected_pivoting_tail_columns =
      (int64_t)tail_columns;
    solver->stats.fast_rejected_pivoting_tail_work = tail_work;
  }
}

static void kls_fill_fast_reject_tail_stats(kls_solver *solver,
                                            UF_long rejected_pivot) {
  if (solver == NULL || solver->symbolic == NULL ||
      solver->numeric == NULL || rejected_pivot == KLS_KLU_EMPTY) {
    return;
  }
  const UF_long block = kls_block_for_pivot(solver, rejected_pivot);
  if (block == KLS_KLU_EMPTY) {
    return;
  }
  const UF_long k1 = solver->symbolic->R[block];
  const UF_long k2 = solver->symbolic->R[block + 1u];
  if (rejected_pivot < k1 || rejected_pivot >= k2) {
    return;
  }
  const UF_long nk = k2 - k1;
  const UF_long local_reject = rejected_pivot - k1;
  solver->stats.fast_rejected_block_start = (int64_t)k1;
  solver->stats.fast_rejected_block_size = (int64_t)nk;
  solver->stats.fast_rejected_suffix_columns =
    (int64_t)(k2 - rejected_pivot);
  solver->stats.fast_rejected_descendant_columns = 0;
  solver->stats.fast_rejected_descendant_work = 0.0;
  solver->stats.fast_rejected_row_tail_columns = 0;
  solver->stats.fast_rejected_row_tail_work = 0.0;
  solver->stats.fast_rejected_etree_columns = 0;
  solver->stats.fast_rejected_etree_work = 0.0;
  solver->stats.fast_rejected_pivoting_tail_columns = 0;
  solver->stats.fast_rejected_pivoting_tail_work = 0.0;
  solver->stats.fast_rejected_pivoting_tail_first = -1;
  solver->stats.fast_rejected_pivoting_tail_last = -1;
  solver->stats.fast_rejected_pivoting_tail_contains_reject = 0;
  solver->stats.fast_rejected_pivoting_tail_topological = 0;
  solver->fast_reject_tail_count = 0;

  kls_fill_fast_reject_row_tail_stats(solver, block, k1, k2,
                                      local_reject);
  kls_fill_fast_reject_etree_tail_stats(solver, block, k1, k2,
                                        local_reject);

  if (nk == 0u || solver->numeric->Uip == NULL ||
      solver->numeric->Ulen == NULL || solver->numeric->Llen == NULL ||
      solver->numeric->LUbx == NULL || solver->numeric->LUbx[block] == NULL) {
    return;
  }

  unsigned char *affected =
    (unsigned char *)calloc((size_t)nk, sizeof(*affected));
  if (affected == NULL) {
    return;
  }

  double *lu = (double *)solver->numeric->LUbx[block];
  const UF_long *uip = solver->numeric->Uip + k1;
  const UF_long *ulen = solver->numeric->Ulen + k1;
  UF_long descendant_columns = 0;
  double descendant_work = 0.0;
  int valid = 1;

  for (UF_long k = local_reject; k < nk; ++k) {
    UF_long *ui = NULL;
    double *ux = NULL;
    UF_long ucol_len = 0;
    kls_klu_get_pointer(lu, uip, ulen, k, &ui, &ux, &ucol_len);
    (void)ux;

    int is_affected = (k == local_reject);
    double column_work = 0.0;
    if (!kls_fast_reject_column_work(solver, block, k1, nk, k,
                                     &column_work)) {
      valid = 0;
      break;
    }
    for (UF_long p = 0; p < ucol_len; ++p) {
      const UF_long dep = ui[p];
      if (dep >= k || dep >= nk) {
        valid = 0;
        break;
      }
      if (affected[dep]) {
        is_affected = 1;
      }
    }
    if (!valid) {
      break;
    }
    if (is_affected) {
      affected[k] = 1u;
      descendant_columns++;
      descendant_work += column_work;
    }
  }

  free(affected);
  if (valid) {
    solver->stats.fast_rejected_descendant_columns =
      (int64_t)descendant_columns;
    solver->stats.fast_rejected_descendant_work = descendant_work;
  }
}

static int kls_fast_repair_covers_remaining_columns(
  const kls_solver *solver,
  UF_long rejected_pivot) {
  if (solver == NULL || solver->symbolic == NULL ||
      solver->symbolic->nblocks == 0u ||
      rejected_pivot == KLS_KLU_EMPTY) {
    return 0;
  }
  if (solver->symbolic->nblocks == 1u ||
      solver->fast_reject_refresh_state == KLS_FAST_REJECT_REFRESH_ALL) {
    return 1;
  }
  if (solver->fast_reject_refresh_state != KLS_FAST_REJECT_REFRESH_PREFIX) {
    return 0;
  }
  const UF_long block = kls_block_for_pivot(solver, rejected_pivot);
  return block != KLS_KLU_EMPTY && block + 1u == solver->symbolic->nblocks;
}

static int kls_recompute_offdiag_from_pinv(kls_solver *solver,
                                           const double *numeric_values) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      numeric_values == NULL || solver->symbolic->Q == NULL ||
      solver->symbolic->R == NULL || solver->numeric->Pinv == NULL ||
      solver->numeric->Offp == NULL || solver->numeric->Offi == NULL ||
      solver->numeric->Offx == NULL ||
      (solver->common.scale > 0 && solver->numeric->Rs == NULL) ||
      (solver->common.scale <= 0 && solver->numeric->Rs != NULL)) {
    return 0;
  }

  const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  trilinos_klu_l_numeric *numeric = solver->numeric;
  UF_long poff = 0;
  UF_long block = 0;
  for (UF_long k = 0; k < solver->n; ++k) {
    while (block + 1u < symbolic->nblocks &&
           symbolic->R[block + 1u] <= k) {
      block++;
    }
    const UF_long k1 = symbolic->R[block];
    const UF_long k2 = symbolic->R[block + 1u];
    const UF_long oldcol = symbolic->Q[k];
    if (oldcol >= solver->n) {
      return 0;
    }
    numeric->Offp[k] = poff;
    for (UF_long p = solver->col_ptr[oldcol];
         p < solver->col_ptr[oldcol + 1u]; ++p) {
      const UF_long oldrow = solver->row_idx[p];
      if (oldrow >= solver->n) {
        return 0;
      }
      const UF_long row = numeric->Pinv[oldrow];
      if (row >= solver->n) {
        return 0;
      }
      if (row < k1) {
        if (poff >= numeric->nzoff) {
          return 0;
        }
        double value = numeric_values[p];
        if (solver->common.scale > 0) {
          const double rs = numeric->Rs[oldrow];
          if (rs == 0.0) {
            return 0;
          }
          value /= rs;
        }
        numeric->Offi[poff] = row;
        ((double *)numeric->Offx)[poff] = value;
        poff++;
      } else if (row >= k2) {
        return 0;
      }
    }
  }
  numeric->Offp[solver->n] = poff;
  return poff == numeric->nzoff;
}

static int kls_recompute_offdiag_suffix_from_pinv(kls_solver *solver,
                                                  const double *numeric_values,
                                                  UF_long col_begin) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      numeric_values == NULL || solver->symbolic->Q == NULL ||
      solver->symbolic->R == NULL || solver->numeric->Pinv == NULL ||
      solver->numeric->Offp == NULL || solver->numeric->Offi == NULL ||
      solver->numeric->Offx == NULL ||
      col_begin > solver->n ||
      (solver->common.scale > 0 && solver->numeric->Rs == NULL) ||
      (solver->common.scale <= 0 && solver->numeric->Rs != NULL)) {
    return 0;
  }

  const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  trilinos_klu_l_numeric *numeric = solver->numeric;
  if (col_begin == solver->n) {
    return numeric->Offp[col_begin] == numeric->nzoff;
  }
  UF_long poff = numeric->Offp[col_begin];
  if (poff > numeric->nzoff) {
    return 0;
  }

  UF_long block = 0;
  while (block + 1u < symbolic->nblocks &&
         symbolic->R[block + 1u] <= col_begin) {
    block++;
  }

  for (UF_long k = col_begin; k < solver->n; ++k) {
    while (block + 1u < symbolic->nblocks &&
           symbolic->R[block + 1u] <= k) {
      block++;
    }
    const UF_long k1 = symbolic->R[block];
    const UF_long k2 = symbolic->R[block + 1u];
    const UF_long oldcol = symbolic->Q[k];
    if (oldcol >= solver->n) {
      return 0;
    }
    numeric->Offp[k] = poff;
    for (UF_long p = solver->col_ptr[oldcol];
         p < solver->col_ptr[oldcol + 1u]; ++p) {
      const UF_long oldrow = solver->row_idx[p];
      if (oldrow >= solver->n) {
        return 0;
      }
      const UF_long row = numeric->Pinv[oldrow];
      if (row >= solver->n) {
        return 0;
      }
      if (row < k1) {
        if (poff >= numeric->nzoff) {
          return 0;
        }
        double value = numeric_values[p];
        if (solver->common.scale > 0) {
          const double rs = numeric->Rs[oldrow];
          if (rs == 0.0) {
            return 0;
          }
          value /= rs;
        }
        numeric->Offi[poff] = row;
        ((double *)numeric->Offx)[poff] = value;
        poff++;
      } else if (row >= k2) {
        return 0;
      }
    }
  }
  numeric->Offp[solver->n] = poff;
  return poff == numeric->nzoff;
}

static int kls_rebuild_numeric_pinv(kls_solver *solver) {
  if (solver == NULL || solver->numeric == NULL ||
      solver->numeric->Pnum == NULL || solver->numeric->Pinv == NULL) {
    return 0;
  }
  for (UF_long k = 0; k < solver->n; ++k) {
    solver->numeric->Pinv[k] = KLS_KLU_EMPTY;
  }
  for (UF_long k = 0; k < solver->n; ++k) {
    const UF_long row = solver->numeric->Pnum[k];
    if (row >= solver->n) {
      return 0;
    }
    solver->numeric->Pinv[row] = k;
  }
  for (UF_long k = 0; k < solver->n; ++k) {
    if (solver->numeric->Pinv[k] == KLS_KLU_EMPTY) {
      return 0;
    }
  }
  return 1;
}

static int kls_klu_column_span_end(const UF_long *offsets,
                                   const UF_long *lengths,
                                   UF_long k,
                                   size_t lusize,
                                   size_t *end_out) {
  if (offsets == NULL || lengths == NULL || end_out == NULL ||
      offsets[k] < 0 || lengths[k] < 0) {
    return 0;
  }
  const size_t offset = (size_t)offsets[k];
  const UF_long length = lengths[k];
  if (length > (UF_long)(SIZE_MAX / sizeof(UF_long))) {
    return 0;
  }
  const size_t index_units = kls_klu_units_for_indices(length);
  if (offset > SIZE_MAX - index_units) {
    return 0;
  }
  size_t end = offset + index_units;
  if ((size_t)length > SIZE_MAX - end) {
    return 0;
  }
  end += (size_t)length;
  if (end > lusize) {
    return 0;
  }
  *end_out = end;
  return 1;
}

static int kls_unfinalize_l_row(const UF_long *final_pblock,
                                UF_long nk,
                                UF_long final_order,
                                UF_long *local_row_out) {
  if (final_pblock == NULL || local_row_out == NULL ||
      final_order < 0 || final_order >= nk) {
    return 0;
  }
  const UF_long local_row = final_pblock[final_order];
  if (local_row < 0 || local_row >= nk) {
    return 0;
  }
  *local_row_out = local_row;
  return 1;
}

static int kls_prefix_l_partition_matches(const UF_long *final_pblock,
                                          UF_long nk,
                                          const UF_long *li,
                                          UF_long llen,
                                          UF_long pivoted_head,
                                          const UF_long *live_pinv) {
  if (final_pblock == NULL || li == NULL || live_pinv == NULL ||
      pivoted_head > llen) {
    return 0;
  }
  for (UF_long p = 0; p < llen; ++p) {
    UF_long local_row = KLS_KLU_EMPTY;
    if (!kls_unfinalize_l_row(final_pblock, nk, li[p], &local_row)) {
      return 0;
    }
    const int pivotal = live_pinv[local_row] >= 0;
    if ((p < pivoted_head) != pivotal) {
      return 0;
    }
  }
  return 1;
}

static int kls_replay_prefix_prune(const kls_solver *solver,
                                   UF_long block,
                                   UF_long k1,
                                   UF_long nk,
                                   UF_long k,
                                   UF_long pivrow,
                                   const UF_long *final_pblock,
                                   const UF_long *live_pinv,
                                   UF_long *lpend) {
  if (solver == NULL || solver->numeric == NULL ||
      solver->numeric->LUbx == NULL || solver->numeric->LUbx[block] == NULL ||
      solver->numeric->Llen == NULL || solver->numeric->Ulen == NULL ||
      solver->numeric->Lip == NULL || solver->numeric->Uip == NULL ||
      final_pblock == NULL || live_pinv == NULL || lpend == NULL ||
      k >= nk || pivrow >= nk) {
    return 0;
  }

  double *lu = (double *)solver->numeric->LUbx[block];
  const UF_long *lip = solver->numeric->Lip + k1;
  const UF_long *llen = solver->numeric->Llen + k1;
  const UF_long *uip = solver->numeric->Uip + k1;
  const UF_long *ulen = solver->numeric->Ulen + k1;

  UF_long *ui = NULL;
  double *ux = NULL;
  UF_long ucol_len = 0;
  kls_klu_get_pointer(lu, uip, ulen, k, &ui, &ux, &ucol_len);
  (void)ux;
  for (UF_long p = 0; p < ucol_len; ++p) {
    const UF_long j = ui[p];
    if (j < 0 || j >= k || j >= nk) {
      return 0;
    }
    if (lpend[j] != KLS_KLU_EMPTY) {
      continue;
    }

    UF_long *li = NULL;
    double *lx = NULL;
    UF_long lcol_len = 0;
    kls_klu_get_pointer(lu, lip, llen, j, &li, &lx, &lcol_len);
    (void)lx;
    int found_pivrow = 0;
    UF_long pivotal_count = 0;
    for (UF_long p2 = 0; p2 < lcol_len; ++p2) {
      UF_long local_row = KLS_KLU_EMPTY;
      if (!kls_unfinalize_l_row(final_pblock, nk, li[p2], &local_row)) {
        return 0;
      }
      if (local_row == pivrow) {
        found_pivrow = 1;
      }
      if (live_pinv[local_row] >= 0) {
        pivotal_count++;
      }
    }
    if (found_pivrow) {
      if (!kls_prefix_l_partition_matches(final_pblock, nk, li, lcol_len,
                                          pivotal_count, live_pinv)) {
        return 0;
      }
      lpend[j] = pivotal_count;
    }
  }
  return 1;
}

static int kls_reconstruct_block_live_prefix_state(
  const kls_solver *solver,
  UF_long block,
  UF_long k1,
  UF_long nk,
  UF_long prefix_cols,
  const UF_long *final_pblock,
  UF_long *live_p,
  UF_long *live_pinv,
  UF_long *final_pinv,
  UF_long *lpend,
  size_t *prefix_lup_out,
  UF_long *prefix_lnz_out,
  UF_long *prefix_unz_out) {
  if (prefix_lup_out != NULL) {
    *prefix_lup_out = 0u;
  }
  if (prefix_lnz_out != NULL) {
    *prefix_lnz_out = 0;
  }
  if (prefix_unz_out != NULL) {
    *prefix_unz_out = 0;
  }
  if (solver == NULL || solver->numeric == NULL ||
      solver->numeric->LUbx == NULL || solver->numeric->LUsize == NULL ||
      solver->numeric->Llen == NULL || solver->numeric->Ulen == NULL ||
      solver->numeric->Lip == NULL || solver->numeric->Uip == NULL ||
      final_pblock == NULL || live_p == NULL || live_pinv == NULL ||
      final_pinv == NULL || lpend == NULL || prefix_cols > nk ||
      k1 > solver->n || nk > solver->n - k1 ||
      block >= solver->numeric->nblocks ||
      solver->numeric->LUbx[block] == NULL) {
    return 0;
  }

  const size_t lusize = solver->numeric->LUsize[block];
  if (lusize == 0u) {
    return 0;
  }
  for (UF_long k = 0; k < nk; ++k) {
    live_p[k] = k;
    live_pinv[k] = FLIP(k);
    final_pinv[k] = KLS_KLU_EMPTY;
    lpend[k] = KLS_KLU_EMPTY;
  }
  for (UF_long k = 0; k < nk; ++k) {
    const UF_long row = final_pblock[k];
    if (row < 0 || row >= nk || final_pinv[row] != KLS_KLU_EMPTY) {
      return 0;
    }
    final_pinv[row] = k;
  }

  double *lu = (double *)solver->numeric->LUbx[block];
  const UF_long *lip = solver->numeric->Lip + k1;
  const UF_long *llen = solver->numeric->Llen + k1;
  const UF_long *uip = solver->numeric->Uip + k1;
  const UF_long *ulen = solver->numeric->Ulen + k1;

  size_t prefix_lup = 0u;
  UF_long prefix_lnz = 0;
  UF_long prefix_unz = 0;
  for (UF_long k = 0; k < prefix_cols; ++k) {
    size_t l_end = 0u;
    size_t u_end = 0u;
    if (!kls_klu_column_span_end(lip, llen, k, lusize, &l_end) ||
        !kls_klu_column_span_end(uip, ulen, k, lusize, &u_end) ||
        uip[k] < 0 || (size_t)uip[k] != l_end ||
        (k + 1u < nk && (lip[k + 1u] < 0 ||
                         (size_t)lip[k + 1u] != u_end))) {
      return 0;
    }

    UF_long *ui = NULL;
    double *ux = NULL;
    UF_long ucol_len = 0;
    kls_klu_get_pointer(lu, uip, ulen, k, &ui, &ux, &ucol_len);
    (void)ux;
    for (UF_long p = 0; p < ucol_len; ++p) {
      if (ui[p] < 0 || ui[p] >= k) {
        return 0;
      }
    }

    const UF_long pivrow = final_pblock[k];
    const UF_long diagrow = live_p[k];
    if (pivrow < 0 || pivrow >= nk || diagrow < 0 || diagrow >= nk ||
        live_pinv[pivrow] >= 0) {
      return 0;
    }
    if (pivrow != diagrow && live_pinv[diagrow] < 0) {
      const UF_long kbar = FLIP(live_pinv[pivrow]);
      if (kbar < 0 || kbar >= nk) {
        return 0;
      }
      live_p[kbar] = diagrow;
      live_pinv[diagrow] = FLIP(kbar);
    }
    live_p[k] = pivrow;
    live_pinv[pivrow] = k;

    UF_long *li = NULL;
    double *lx = NULL;
    UF_long lcol_len = 0;
    kls_klu_get_pointer(lu, lip, llen, k, &li, &lx, &lcol_len);
    (void)lx;
    for (UF_long p = 0; p < lcol_len; ++p) {
      const UF_long final_order = li[p];
      if (final_order <= k || final_order >= nk) {
        return 0;
      }
      UF_long local_row = KLS_KLU_EMPTY;
      if (!kls_unfinalize_l_row(final_pblock, nk, final_order,
                                &local_row) ||
          live_pinv[local_row] >= 0) {
        return 0;
      }
    }

    if (!kls_replay_prefix_prune(solver, block, k1, nk, k, pivrow,
                                 final_pblock, live_pinv, lpend)) {
      return 0;
    }
    prefix_lup = u_end;
    prefix_lnz += llen[k] + 1u;
    prefix_unz += ulen[k] + 1u;
  }

  if (prefix_lup_out != NULL) {
    *prefix_lup_out = prefix_lup;
  }
  if (prefix_lnz_out != NULL) {
    *prefix_lnz_out = prefix_lnz;
  }
  if (prefix_unz_out != NULL) {
    *prefix_unz_out = prefix_unz;
  }
  return 1;
}

static int kls_repaired_tail_prefix_state_valid(const kls_solver *solver,
                                                UF_long block,
                                                UF_long k1,
                                                UF_long nk,
                                                UF_long prefix_cols,
                                                const UF_long *final_pblock) {
  if (prefix_cols == 0u) {
    return 1;
  }
  if (solver == NULL || final_pblock == NULL || prefix_cols > nk ||
      nk > (UF_long)(SIZE_MAX / (4u * sizeof(UF_long)))) {
    return 0;
  }
  UF_long *scratch =
    (UF_long *)malloc(4u * (size_t)nk * sizeof(*scratch));
  if (scratch == NULL) {
    return 0;
  }
  UF_long *live_p = scratch;
  UF_long *live_pinv = scratch + nk;
  UF_long *final_pinv = scratch + 2u * nk;
  UF_long *lpend = scratch + 3u * nk;
  size_t prefix_lup = 0u;
  UF_long prefix_lnz = 0;
  UF_long prefix_unz = 0;
  const int ok = kls_reconstruct_block_live_prefix_state(
    solver, block, k1, nk, prefix_cols, final_pblock, live_p, live_pinv,
    final_pinv, lpend, &prefix_lup, &prefix_lnz, &prefix_unz);
  (void)prefix_lup;
  (void)prefix_lnz;
  (void)prefix_unz;
  free(scratch);
  return ok;
}

static Int kls_tail_dfs(Int j,
                        Int k,
                        Int pinv[],
                        Int llen[],
                        Int lip[],
                        Int stack[],
                        Int flag[],
                        Int lpend[],
                        Int top,
                        Unit lu[],
                        Int *lik,
                        Int *plength,
                        Int ap_pos[]) {
  Int l_length = *plength;
  Int head = 0;
  stack[0] = j;

  while (head >= 0) {
    j = stack[head];
    const Int jnew = pinv[j];
    if (jnew < 0 || jnew >= k) {
      return -1;
    }

    if (flag[j] != k) {
      flag[j] = k;
      ap_pos[head] =
        (lpend[jnew] == TRILINOS_KLU_EMPTY) ? llen[jnew] : lpend[jnew];
    }

    Int *li = (Int *)(lu + lip[jnew]);
    Int pos = 0;
    for (pos = --ap_pos[head]; pos >= 0; --pos) {
      const Int i = li[pos];
      if (flag[i] == k) {
        continue;
      }
      if (pinv[i] >= 0) {
        ap_pos[head] = pos;
        stack[++head] = i;
        break;
      }
      flag[i] = k;
      lik[l_length++] = i;
    }

    if (pos == -1) {
      head--;
      stack[--top] = j;
    }
  }

  *plength = l_length;
  return top;
}

static Int kls_tail_lsolve_symbolic(Int n,
                                    Int k,
                                    Int ap[],
                                    Int ai[],
                                    Int q[],
                                    Int pinv[],
                                    Int stack[],
                                    Int flag[],
                                    Int lpend[],
                                    Int ap_pos[],
                                    Unit lu[],
                                    Int lup,
                                    Int llen[],
                                    Int lip[],
                                    Int k1,
                                    Int psinv[]) {
  Int top = n;
  Int l_length = 0;
  Int *lik = (Int *)(lu + lup);
  const Int kglobal = k + k1;
  const Int oldcol = q[kglobal];

  for (Int p = ap[oldcol]; p < ap[oldcol + 1]; ++p) {
    const Int i = psinv[ai[p]] - k1;
    if (i < 0) {
      continue;
    }
    if (i >= n) {
      return -1;
    }
    if (flag[i] == k) {
      continue;
    }
    if (pinv[i] >= 0) {
      top = kls_tail_dfs(i, k, pinv, llen, lip, stack, flag, lpend,
                         top, lu, lik, &l_length, ap_pos);
      if (top < 0) {
        return -1;
      }
    } else {
      flag[i] = k;
      lik[l_length++] = i;
    }
  }

  llen[k] = l_length;
  return top;
}

static int kls_tail_construct_column(Int k,
                                     Int ap[],
                                     Int ai[],
                                     Entry ax[],
                                     Int q[],
                                     Entry x[],
                                     Int n,
                                     Int k1,
                                     Int psinv[],
                                     double rs[],
                                     Int scale,
                                     Int offp[],
                                     Int offi[],
                                     Entry offx[],
                                     Int nzoff) {
  const Int kglobal = k + k1;
  const Int oldcol = q[kglobal];
  Int poff = offp[kglobal];
  if (poff < 0 || poff > nzoff) {
    return 0;
  }

  for (Int p = ap[oldcol]; p < ap[oldcol + 1]; ++p) {
    const Int oldrow = ai[p];
    const Int i = psinv[oldrow] - k1;
    Entry aik = ax[p];
    if (scale > 0) {
      if (rs == NULL || rs[oldrow] == 0.0) {
        return 0;
      }
      SCALE_DIV(aik, rs[oldrow]);
    }
    if (i < 0) {
      if (poff >= nzoff) {
        return 0;
      }
      if (offi != NULL) {
        offi[poff] = oldrow;
      }
      if (offx != NULL) {
        offx[poff] = aik;
      }
      poff++;
    } else if (i >= n) {
      return 0;
    } else {
      x[i] = aik;
    }
  }
  offp[kglobal + 1] = poff;
  return 1;
}

static void kls_tail_lsolve_numeric(Int pinv[],
                                    Unit *lu,
                                    Int stack[],
                                    Int lip[],
                                    Int top,
                                    Int n,
                                    Int llen[],
                                    Entry x[]) {
  for (Int s = top; s < n; ++s) {
    const Int j = stack[s];
    const Int jnew = pinv[j];
    Entry xj = x[j];
    Entry *lx = NULL;
    Int *li = NULL;
    Int len = 0;
    GET_POINTER(lu, lip, llen, li, lx, jnew, len);
    for (Int p = 0; p < len; ++p) {
      MULT_SUB(x[li[p]], lx[p], xj);
    }
  }
}

static Int kls_tail_lpivot(Int diagrow,
                           Int *p_pivrow,
                           Entry *p_pivot,
                           double *p_abs_pivot,
                           double tol,
                           Entry x[],
                           Unit *lu,
                           Int lip[],
                           Int llen[],
                           Int k,
                           Int n,
                           Int pinv[],
                           Int *p_firstrow,
                           TRILINOS_KLU_common *common) {
  Int pivrow = TRILINOS_KLU_EMPTY;
  if (llen[k] == 0) {
    if (common->halt_if_singular) {
      return FALSE;
    }
    for (Int firstrow = *p_firstrow; firstrow < n; ++firstrow) {
      if (pinv[firstrow] < 0) {
        pivrow = firstrow;
        *p_firstrow = firstrow;
        break;
      }
    }
    if (pivrow < 0 || pivrow >= n) {
      return FALSE;
    }
    *p_pivrow = pivrow;
    *p_pivot = 0.0;
    *p_abs_pivot = 0.0;
    return FALSE;
  }

  Int *li = NULL;
  Entry *lx = NULL;
  Int len = 0;
  Int ppivrow = TRILINOS_KLU_EMPTY;
  Int pdiag = TRILINOS_KLU_EMPTY;
  double abs_pivot = -1.0;

  llen[k]--;
  GET_POINTER(lu, lip, llen, li, lx, k, len);
  const Int last_row_index = li[len];

  for (Int p = 0; p < len; ++p) {
    const Int i = li[p];
    const Entry xi = x[i];
    x[i] = 0.0;
    lx[p] = xi;
    const double xabs = fabs(xi);
    if (i == diagrow) {
      pdiag = p;
    }
    if (xabs > abs_pivot) {
      abs_pivot = xabs;
      ppivrow = p;
    }
  }

  const double last_abs = fabs(x[last_row_index]);
  if (last_abs > abs_pivot) {
    abs_pivot = last_abs;
    ppivrow = TRILINOS_KLU_EMPTY;
  }
  if (last_row_index == diagrow) {
    if (last_abs >= tol * abs_pivot) {
      abs_pivot = last_abs;
      ppivrow = TRILINOS_KLU_EMPTY;
    }
  } else if (pdiag != TRILINOS_KLU_EMPTY) {
    const double diag_abs = fabs(lx[pdiag]);
    if (diag_abs >= tol * abs_pivot) {
      abs_pivot = diag_abs;
      ppivrow = pdiag;
    }
  }

  Entry pivot = 0.0;
  if (ppivrow != TRILINOS_KLU_EMPTY) {
    pivrow = li[ppivrow];
    pivot = lx[ppivrow];
    li[ppivrow] = last_row_index;
    lx[ppivrow] = x[last_row_index];
  } else {
    pivrow = last_row_index;
    pivot = x[last_row_index];
  }
  x[last_row_index] = 0.0;

  *p_pivrow = pivrow;
  *p_pivot = pivot;
  *p_abs_pivot = abs_pivot;
  if (pivrow < 0 || pivrow >= n) {
    return FALSE;
  }
  if (pivot == 0.0 && common->halt_if_singular) {
    return FALSE;
  }

  for (Int p = 0; p < llen[k]; ++p) {
    DIV(lx[p], lx[p], pivot);
  }
  return TRUE;
}

static void kls_tail_prune(Int lpend[],
                           Int pinv[],
                           Int k,
                           Int pivrow,
                           Unit *lu,
                           Int uip[],
                           Int lip[],
                           Int ulen[],
                           Int llen[]) {
  Int *ui = NULL;
  Entry *ux = NULL;
  Int ucol_len = 0;
  GET_POINTER(lu, uip, ulen, ui, ux, k, ucol_len);
  (void)ux;
  for (Int p = 0; p < ucol_len; ++p) {
    const Int j = ui[p];
    if (j < 0 || j >= k || lpend[j] != TRILINOS_KLU_EMPTY) {
      continue;
    }
    Int *li = NULL;
    Entry *lx = NULL;
    Int lcol_len = 0;
    GET_POINTER(lu, lip, llen, li, lx, j, lcol_len);
    for (Int p2 = 0; p2 < lcol_len; ++p2) {
      if (pivrow != li[p2]) {
        continue;
      }
      Int phead = 0;
      Int ptail = llen[j];
      while (phead < ptail) {
        const Int i = li[phead];
        if (pinv[i] >= 0) {
          phead++;
        } else {
          ptail--;
          li[phead] = li[ptail];
          li[ptail] = i;
          const Entry tmp = lx[phead];
          lx[phead] = lx[ptail];
          lx[ptail] = tmp;
        }
      }
      lpend[j] = ptail;
      break;
    }
  }
}

static int kls_copy_live_prefix_lu(Unit *new_lu,
                                   size_t lusize,
                                   const Unit *old_lu,
                                   const UF_long *old_lip,
                                   const UF_long *old_llen,
                                   const UF_long *old_uip,
                                   const UF_long *old_ulen,
                                   const Entry *old_udiag,
                                   const UF_long *final_pblock,
                                   UF_long nk,
                                   UF_long prefix_cols,
                                   UF_long *lip,
                                   UF_long *llen,
                                   UF_long *uip,
                                   UF_long *ulen,
                                   Entry *udiag,
                                   size_t *lup_out) {
  if (new_lu == NULL || old_lu == NULL || final_pblock == NULL ||
      lup_out == NULL) {
    return 0;
  }
  size_t lup = 0u;
  for (UF_long k = 0; k < prefix_cols; ++k) {
    const UF_long l_len = old_llen[k];
    const UF_long u_len = old_ulen[k];
    const size_t l_units = kls_klu_units_for_indices(l_len);
    const size_t u_units = kls_klu_units_for_indices(u_len);
    if (lup > lusize || l_units > lusize - lup ||
        l_len > (UF_long)(lusize - lup - l_units)) {
      return 0;
    }
    lip[k] = (UF_long)lup;
    llen[k] = l_len;
    UF_long *dst_li = (UF_long *)(new_lu + lup);
    Entry *dst_lx = (Entry *)(new_lu + lup + l_units);
    Int *src_li = NULL;
    Entry *src_lx = NULL;
    Int src_len = 0;
    GET_POINTER((Unit *)old_lu, old_lip, old_llen, src_li, src_lx, k,
                src_len);
    if ((UF_long)src_len != l_len) {
      return 0;
    }
    for (UF_long p = 0; p < l_len; ++p) {
      const UF_long final_order = src_li[p];
      if (final_order >= nk || final_pblock[final_order] >= nk) {
        return 0;
      }
      dst_li[p] = final_pblock[final_order];
      dst_lx[p] = src_lx[p];
    }
    lup += l_units + (size_t)l_len;

    if (lup > lusize || u_units > lusize - lup ||
        u_len > (UF_long)(lusize - lup - u_units)) {
      return 0;
    }
    uip[k] = (UF_long)lup;
    ulen[k] = u_len;
    UF_long *dst_ui = (UF_long *)(new_lu + lup);
    Entry *dst_ux = (Entry *)(new_lu + lup + u_units);
    Int *src_ui = NULL;
    Entry *src_ux = NULL;
    src_len = 0;
    GET_POINTER((Unit *)old_lu, old_uip, old_ulen, src_ui, src_ux, k,
                src_len);
    if ((UF_long)src_len != u_len) {
      return 0;
    }
    for (UF_long p = 0; p < u_len; ++p) {
      if (src_ui[p] >= (Int)k) {
        return 0;
      }
      dst_ui[p] = (UF_long)src_ui[p];
      dst_ux[p] = src_ux[p];
    }
    lup += u_units + (size_t)u_len;
    udiag[k] = old_udiag[k];
  }
  *lup_out = lup;
  return 1;
}

static int kls_try_pivot_tail_restart_rejected_block(
  kls_solver *solver,
  double *numeric_values,
  UF_long block,
  UF_long k1,
  UF_long nk,
  UF_long local_reject,
  const UF_long *old_pblock,
  const UF_long *psinv,
  Unit **new_lu_out,
  size_t *new_size_out,
  UF_long *lnz_block_out,
  UF_long *unz_block_out,
  UF_long *pblock) {
  if (new_lu_out != NULL) {
    *new_lu_out = NULL;
  }
  if (new_size_out != NULL) {
    *new_size_out = 0u;
  }
  if (lnz_block_out != NULL) {
    *lnz_block_out = 0;
  }
  if (unz_block_out != NULL) {
    *unz_block_out = 0;
  }
  const int reusable_prefix_state =
    solver != NULL &&
    (solver->fast_reject_refresh_state == KLS_FAST_REJECT_REFRESH_PREFIX ||
     solver->fast_reject_refresh_state == KLS_FAST_REJECT_REFRESH_ALL);
  const int scaled = solver != NULL && solver->common.scale > 0;
  if (solver == NULL || solver->numeric == NULL ||
      solver->symbolic == NULL || numeric_values == NULL ||
      old_pblock == NULL || psinv == NULL || new_lu_out == NULL ||
      new_size_out == NULL || lnz_block_out == NULL ||
      unz_block_out == NULL || pblock == NULL ||
      !reusable_prefix_state ||
      local_reject == 0u || local_reject >= nk ||
      block >= solver->numeric->nblocks ||
      solver->numeric->LUbx == NULL ||
      solver->numeric->LUbx[block] == NULL ||
      solver->numeric->LUsize == NULL ||
      solver->numeric->LUsize[block] == 0u ||
      solver->numeric->Xwork == NULL ||
      solver->numeric->Llen == NULL || solver->numeric->Ulen == NULL ||
      solver->numeric->Lip == NULL || solver->numeric->Uip == NULL ||
      solver->numeric->Udiag == NULL || solver->numeric->Offp == NULL ||
      solver->numeric->Offi == NULL || solver->numeric->Offx == NULL) {
    return 0;
  }
  if ((scaled && (solver->numeric->Rs == NULL ||
                  solver->fast_reject_refresh_state !=
                    KLS_FAST_REJECT_REFRESH_PREFIX)) ||
      (!scaled && solver->numeric->Rs != NULL)) {
    return 0;
  }

  if (nk > (UF_long)(SIZE_MAX / (7u * sizeof(UF_long)))) {
    return 0;
  }
  UF_long *scratch =
    (UF_long *)malloc(7u * (size_t)nk * sizeof(*scratch));
  if (scratch == NULL) {
    return 0;
  }
  UF_long *live_p = scratch;
  UF_long *live_pinv = scratch + nk;
  UF_long *final_pinv = scratch + 2u * nk;
  UF_long *lpend = scratch + 3u * nk;
  UF_long *stack = scratch + 4u * nk;
  UF_long *flag = scratch + 5u * nk;
  UF_long *ap_pos = scratch + 6u * nk;

  size_t prefix_lup = 0u;
  UF_long prefix_lnz = 0;
  UF_long prefix_unz = 0;
  if (!kls_reconstruct_block_live_prefix_state(
        solver, block, k1, nk, local_reject, old_pblock, live_p,
        live_pinv, final_pinv, lpend, &prefix_lup, &prefix_lnz,
        &prefix_unz)) {
    free(scratch);
    return 0;
  }
  (void)final_pinv;

  const size_t old_lusize = solver->numeric->LUsize[block];
  Unit *new_lu =
    (Unit *)TRILINOS_KLU_malloc(old_lusize, sizeof(Unit), &solver->common);
  if (new_lu == NULL) {
    free(scratch);
    return 0;
  }

  UF_long *offp = (UF_long *)malloc(((size_t)solver->n + 1u) * sizeof(*offp));
  if (offp == NULL) {
    free(offp);
    (void)TRILINOS_KLU_free(new_lu, old_lusize, sizeof(Unit),
                            &solver->common);
    free(scratch);
    return 0;
  }
  memcpy(offp, solver->numeric->Offp,
         ((size_t)solver->n + 1u) * sizeof(*offp));

  Unit *old_lu = (Unit *)solver->numeric->LUbx[block];
  UF_long *lip = solver->numeric->Lip + k1;
  UF_long *llen = solver->numeric->Llen + k1;
  UF_long *uip = solver->numeric->Uip + k1;
  UF_long *ulen = solver->numeric->Ulen + k1;
  Entry *udiag = ((Entry *)solver->numeric->Udiag) + k1;

  size_t lup = 0u;
  if (!kls_copy_live_prefix_lu(
        new_lu, old_lusize, old_lu, lip, llen, uip, ulen, udiag,
        old_pblock, nk, local_reject, lip, llen, uip, ulen, udiag,
        &lup) ||
      lup != prefix_lup) {
    free(offp);
    (void)TRILINOS_KLU_free(new_lu, old_lusize, sizeof(Unit),
                            &solver->common);
    free(scratch);
    return 0;
  }

  for (UF_long k = 0; k < nk; ++k) {
    pblock[k] = live_p[k];
    flag[k] = KLS_KLU_EMPTY;
    ((Entry *)solver->numeric->Xwork)[k] = 0.0;
  }

  size_t lusize = old_lusize;
  UF_long lnz = prefix_lnz;
  UF_long unz = prefix_unz;
  Int firstrow = 0;
  while (firstrow < (Int)nk && live_pinv[firstrow] >= 0) {
    firstrow++;
  }
  const double tol = solver->common.tol;
  const double memgrow = solver->common.memgrow;
  for (UF_long kk = local_reject; kk < nk; ++kk) {
    const Int k = (Int)kk;
    const double nunits =
      DUNITS(Int, (Int)nk - k) + DUNITS(Int, k) +
      DUNITS(Entry, (Int)nk - k) + DUNITS(Entry, k);
    const double xsize = (double)lup + nunits;
    if (xsize > (double)lusize) {
      if (INT_OVERFLOW(memgrow * (double)lusize + 4.0 * (double)nk + 1.0)) {
        solver->common.status = TRILINOS_KLU_TOO_LARGE;
        goto fail;
      }
      size_t grown = (size_t)(memgrow * (double)lusize) + 2u * (size_t)nk + 1u;
      Unit *grown_lu =
        (Unit *)TRILINOS_KLU_realloc(grown, lusize, sizeof(Unit), new_lu,
                                     &solver->common);
      if (solver->common.status == TRILINOS_KLU_OUT_OF_MEMORY ||
          grown_lu == NULL) {
        goto fail;
      }
      new_lu = grown_lu;
      lusize = grown;
      solver->common.nrealloc++;
    }

    lip[k] = (UF_long)lup;
    const Int top =
      kls_tail_lsolve_symbolic((Int)nk, k, solver->col_ptr, solver->row_idx,
                               solver->symbolic->Q, live_pinv, stack, flag,
                               lpend, ap_pos, new_lu, (Int)lup, llen, lip,
                               (Int)k1, (Int *)psinv);
    if (top < 0) {
      solver->common.status = TRILINOS_KLU_INVALID;
      goto fail;
    }
    /* Trial off-block row/value entries are discarded: accepted repairs
       rebuild Offi/Offx from the final Pinv before publishing stats. */
    if (!kls_tail_construct_column(
          k, solver->col_ptr, solver->row_idx, numeric_values,
          solver->symbolic->Q, (Entry *)solver->numeric->Xwork, (Int)nk,
          (Int)k1, (Int *)psinv,
          scaled ? solver->numeric->Rs : NULL,
          scaled ? (Int)solver->common.scale : 0,
          offp, NULL, NULL,
          (Int)solver->numeric->nzoff)) {
      solver->common.status = TRILINOS_KLU_INVALID;
      goto fail;
    }
    kls_tail_lsolve_numeric(live_pinv, new_lu, stack, lip, top, (Int)nk,
                            llen, (Entry *)solver->numeric->Xwork);

    const Int diagrow = pblock[k];
    Int pivrow = TRILINOS_KLU_EMPTY;
    Entry pivot = 0.0;
    double abs_pivot = 0.0;
    if (!kls_tail_lpivot(diagrow, &pivrow, &pivot, &abs_pivot, tol,
                         (Entry *)solver->numeric->Xwork, new_lu, lip, llen,
                         k, (Int)nk, live_pinv, &firstrow,
                         &solver->common)) {
      solver->common.status = TRILINOS_KLU_SINGULAR;
      if (solver->common.numerical_rank == TRILINOS_KLU_EMPTY) {
        solver->common.numerical_rank = k + (Int)k1;
        solver->common.singular_col = solver->symbolic->Q[k + (Int)k1];
      }
      if (solver->common.halt_if_singular) {
        goto fail;
      }
    }
    if (pivrow < 0 || pivrow >= (Int)nk || live_pinv[pivrow] >= 0) {
      solver->common.status = TRILINOS_KLU_INVALID;
      goto fail;
    }

    uip[k] = lip[k] + (UF_long)UNITS(Int, llen[k]) +
             (UF_long)UNITS(Entry, llen[k]);
    lup += (size_t)UNITS(Int, llen[k]) + (size_t)UNITS(Entry, llen[k]);
    ulen[k] = (UF_long)((Int)nk - top);
    Int *ui = NULL;
    Entry *ux = NULL;
    Int u_len = 0;
    GET_POINTER(new_lu, uip, ulen, ui, ux, k, u_len);
    for (Int p = top, i = 0; p < (Int)nk; ++p, ++i) {
      const Int j = stack[p];
      ui[i] = live_pinv[j];
      ux[i] = ((Entry *)solver->numeric->Xwork)[j];
      ((Entry *)solver->numeric->Xwork)[j] = 0.0;
    }
    lup += (size_t)UNITS(Int, ulen[k]) + (size_t)UNITS(Entry, ulen[k]);
    udiag[k] = pivot;

    if (pivrow != diagrow) {
      solver->common.noffdiag++;
      if (live_pinv[diagrow] < 0) {
        const Int kbar = FLIP(live_pinv[pivrow]);
        if (kbar < 0 || kbar >= (Int)nk) {
          solver->common.status = TRILINOS_KLU_INVALID;
          goto fail;
        }
        pblock[kbar] = diagrow;
        live_pinv[diagrow] = FLIP(kbar);
      }
    }
    pblock[k] = pivrow;
    live_pinv[pivrow] = k;

    kls_tail_prune(lpend, live_pinv, k, pivrow, new_lu, uip, lip, ulen,
                   llen);
    lnz += llen[k] + 1u;
    unz += ulen[k] + 1u;
  }

  for (UF_long k = 0; k < nk; ++k) {
    Int *li = NULL;
    Entry *lx = NULL;
    Int l_len = 0;
    GET_POINTER(new_lu, lip, llen, li, lx, (Int)k, l_len);
    (void)lx;
    for (Int p = 0; p < l_len; ++p) {
      const Int row = li[p];
      if (row < 0 || row >= (Int)nk || live_pinv[row] < 0) {
        solver->common.status = TRILINOS_KLU_INVALID;
        goto fail;
      }
      li[p] = live_pinv[row];
    }
  }

  {
    Unit *shrunk =
      (Unit *)TRILINOS_KLU_realloc(lup, lusize, sizeof(Unit), new_lu,
                                   &solver->common);
    if (shrunk != NULL) {
      new_lu = shrunk;
      lusize = lup;
    }
  }
  free(offp);
  free(scratch);
  *new_lu_out = new_lu;
  *new_size_out = lusize;
  *lnz_block_out = lnz;
  *unz_block_out = unz;
  return 1;

fail:
  for (UF_long k = 0; k < nk; ++k) {
    ((Entry *)solver->numeric->Xwork)[k] = 0.0;
  }
  free(offp);
  (void)TRILINOS_KLU_free(new_lu, lusize, sizeof(Unit), &solver->common);
  free(scratch);
  if (new_lu_out != NULL) {
    *new_lu_out = NULL;
  }
  if (new_size_out != NULL) {
    *new_size_out = 0u;
  }
  return 0;
}

static void kls_record_fast_repaired_block_stats(kls_solver *solver,
                                                 UF_long block,
                                                 UF_long k1,
                                                 UF_long nk,
                                                 UF_long rejected_pivot,
                                                 const UF_long *old_pblock,
                                                 const UF_long *pblock) {
  if (solver == NULL || old_pblock == NULL || pblock == NULL ||
      rejected_pivot < k1 || rejected_pivot >= k1 + nk) {
    return;
  }

  const UF_long local_reject = rejected_pivot - k1;
  UF_long first_changed = KLS_KLU_EMPTY;
  UF_long prefix_changed = 0;
  UF_long suffix_changed = 0;
  int valid = 1;
  for (UF_long k = 0; k < nk; ++k) {
    const UF_long old_local_row = old_pblock[k];
    const UF_long local_row = pblock[k];
    if (old_local_row >= nk || local_row >= nk) {
      valid = 0;
      break;
    }
    if (local_row != old_local_row) {
      if (first_changed == KLS_KLU_EMPTY) {
        first_changed = k1 + k;
      }
      if (k < local_reject) {
        prefix_changed++;
      } else {
        suffix_changed++;
      }
    }
  }
  if (!valid) {
    return;
  }

  const int64_t repaired_row = (int64_t)(k1 + pblock[local_reject]);
  solver->stats.fast_repaired_pivot_row = repaired_row;
  solver->stats.fast_repaired_pivot_matches_tail_candidate =
    solver->stats.fast_rejected_tail_candidate_row >= 0 &&
    solver->stats.fast_rejected_tail_candidate_row == repaired_row;
  solver->stats.fast_repaired_first_changed_pivot =
    first_changed == KLS_KLU_EMPTY ? -1 : (int64_t)first_changed;
  solver->stats.fast_repaired_prefix_changed_pivots =
    (int64_t)prefix_changed;
  solver->stats.fast_repaired_suffix_changed_pivots =
    (int64_t)suffix_changed;
  const int pivoting_tail_ready =
    solver->stats.fast_rejected_pivoting_tail_columns > 0 &&
    solver->stats.fast_rejected_pivoting_tail_work > 0.0 &&
    solver->stats.fast_rejected_pivoting_tail_contains_reject &&
    solver->stats.fast_rejected_pivoting_tail_topological;
  const int serial_tail_state =
    solver->fast_reject_refresh_state == KLS_FAST_REJECT_REFRESH_PREFIX ||
    solver->fast_reject_refresh_state == KLS_FAST_REJECT_REFRESH_ALL;
  solver->stats.fast_repaired_tail_restart_ready =
    serial_tail_state &&
    local_reject > 0u &&
    prefix_changed == 0u &&
    (first_changed == KLS_KLU_EMPTY || first_changed >= rejected_pivot) &&
    suffix_changed > 0u &&
    pivoting_tail_ready;
  if (solver->stats.fast_repaired_tail_restart_ready &&
      !kls_repaired_tail_prefix_state_valid(solver, block, k1, nk,
                                            local_reject, pblock)) {
    solver->stats.fast_repaired_tail_restart_ready = 0;
  }
  double block_work = 0.0;
  if (kls_fast_reject_block_work(solver, block, k1, nk, &block_work)) {
    solver->stats.fast_repaired_block_work = block_work;
    if (solver->stats.fast_repaired_tail_restart_ready) {
      UF_long tail_columns = 0;
      double tail_work = 0.0;
      if (kls_fast_reject_suffix_work(solver, block, k1, nk, local_reject,
                                      &tail_columns, &tail_work)) {
        solver->stats.fast_repaired_tail_restart_columns =
          (int64_t)tail_columns;
        solver->stats.fast_repaired_tail_restart_work = tail_work;
        solver->stats.fast_repaired_tail_restart_saved_work =
          block_work > tail_work ? block_work - tail_work : 0.0;
      } else {
        solver->stats.fast_repaired_tail_restart_ready = 0;
      }
    }
  }
}

static int kls_pivot_restart_rejected_block(kls_solver *solver,
                                            double *numeric_values,
                                            UF_long rejected_pivot) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      numeric_values == NULL || rejected_pivot == KLS_KLU_EMPTY ||
      solver->symbolic->R == NULL || solver->symbolic->P == NULL ||
      solver->symbolic->Q == NULL ||
      solver->symbolic->Lnz == NULL || solver->numeric->LUbx == NULL ||
      solver->numeric->LUsize == NULL || solver->numeric->Udiag == NULL ||
      solver->numeric->Llen == NULL || solver->numeric->Ulen == NULL ||
      solver->numeric->Lip == NULL || solver->numeric->Uip == NULL ||
      solver->numeric->Pnum == NULL || solver->numeric->Pinv == NULL ||
      solver->numeric->Xwork == NULL || solver->numeric->Iwork == NULL) {
    return 0;
  }
  const int scaled = solver->common.scale > 0;
  if ((scaled && solver->numeric->Rs == NULL) ||
      (!scaled && solver->numeric->Rs != NULL)) {
    return 0;
  }
  if (scaled &&
      solver->fast_reject_refresh_state == KLS_FAST_REJECT_REFRESH_ALL) {
    return 0;
  }

  const UF_long block = kls_block_for_pivot(solver, rejected_pivot);
  if (block == KLS_KLU_EMPTY) {
    return 0;
  }
  const UF_long k1 = solver->symbolic->R[block];
  const UF_long k2 = solver->symbolic->R[block + 1u];
  const UF_long nk = k2 - k1;
  if (nk <= 1u || nk > solver->symbolic->maxblock) {
    return 0;
  }
  if (rejected_pivot < k1 || rejected_pivot >= k2) {
    return 0;
  }
  const UF_long local_reject = rejected_pivot - k1;

  UF_long *psinv = solver->numeric->Pinv;
  for (UF_long k = 0; k < solver->n; ++k) {
    const UF_long row = solver->symbolic->P[k];
    if (row >= solver->n) {
      (void)kls_rebuild_numeric_pinv(solver);
      return 0;
    }
    psinv[row] = k;
  }

  UF_long old_lnz_block = 0;
  UF_long old_unz_block = 0;
  for (UF_long k = k1; k < k2; ++k) {
    old_lnz_block += solver->numeric->Llen[k] + 1u;
    old_unz_block += solver->numeric->Ulen[k] + 1u;
  }

  UF_long *old_pblock =
    (UF_long *)malloc((size_t)nk * sizeof(*old_pblock));
  if (old_pblock == NULL) {
    (void)kls_rebuild_numeric_pinv(solver);
    return 0;
  }
  for (UF_long k = 0; k < nk; ++k) {
    const UF_long old_row = solver->numeric->Pnum[k1 + k];
    if (old_row >= solver->n || psinv[old_row] < k1 ||
        psinv[old_row] >= k2) {
      free(old_pblock);
      (void)kls_rebuild_numeric_pinv(solver);
      return 0;
    }
    old_pblock[k] = psinv[old_row] - k1;
  }

  double lsize = 0.0;
  if (solver->symbolic->Lnz[block] < 0.0) {
    lsize = -(solver->common.initmem);
  } else {
    lsize = solver->common.initmem_amd * solver->symbolic->Lnz[block] +
            (double)nk;
  }

  Unit *new_lu = NULL;
  UF_long lnz_block = 0;
  UF_long unz_block = 0;
  UF_long *pblock =
    solver->numeric->Iwork + 5u * (size_t)solver->symbolic->maxblock;
  const int old_status = (int)solver->common.status;
  const UF_long old_numerical_rank = (UF_long)solver->common.numerical_rank;
  const UF_long old_singular_col = (UF_long)solver->common.singular_col;
  const UF_long old_noffdiag = (UF_long)solver->common.noffdiag;
  solver->common.status = TRILINOS_KLU_OK;
  solver->common.numerical_rank = KLS_KLU_EMPTY;
  solver->common.singular_col = KLS_KLU_EMPTY;

  size_t new_size = 0u;
  int tail_restart_used =
    kls_try_pivot_tail_restart_rejected_block(
      solver, numeric_values, block, k1, nk, local_reject, old_pblock,
      psinv, &new_lu, &new_size, &lnz_block, &unz_block, pblock);
  if (!tail_restart_used) {
    solver->common.status = TRILINOS_KLU_OK;
    solver->common.numerical_rank = KLS_KLU_EMPTY;
    solver->common.singular_col = KLS_KLU_EMPTY;
    solver->common.noffdiag = old_noffdiag;
    new_size =
      TRILINOS_KLU_kernel_factor(nk, solver->col_ptr, solver->row_idx,
                                 numeric_values, solver->symbolic->Q,
                                 lsize, &new_lu,
                                 ((double *)solver->numeric->Udiag) + k1,
                                 solver->numeric->Llen + k1,
                                 solver->numeric->Ulen + k1,
                                 solver->numeric->Lip + k1,
                                 solver->numeric->Uip + k1,
                                 pblock, &lnz_block, &unz_block,
                                 (double *)solver->numeric->Xwork,
                                 solver->numeric->Iwork, k1,
                                 psinv, scaled ? solver->numeric->Rs : NULL,
                                 solver->numeric->Offp,
                                 solver->numeric->Offi,
                                 (double *)solver->numeric->Offx,
                                 &solver->common);
  }
  if (new_size == 0 || new_lu == NULL || solver->common.status < 0 ||
      (solver->common.status == TRILINOS_KLU_SINGULAR &&
       solver->common.halt_if_singular)) {
    if (new_lu != NULL) {
      (void)TRILINOS_KLU_free(new_lu, new_size, sizeof(Unit), &solver->common);
    }
    free(old_pblock);
    solver->common.status = old_status;
    solver->common.numerical_rank = old_numerical_rank;
    solver->common.singular_col = old_singular_col;
    solver->common.noffdiag = old_noffdiag;
    (void)kls_rebuild_numeric_pinv(solver);
    return 0;
  }

  (void)TRILINOS_KLU_free(solver->numeric->LUbx[block],
                          solver->numeric->LUsize[block],
                          sizeof(Unit), &solver->common);
  solver->numeric->LUbx[block] = new_lu;
  solver->numeric->LUsize[block] = new_size;

  for (UF_long k = 0; k < nk; ++k) {
    const UF_long local_row = pblock[k];
    if (local_row >= nk) {
      free(old_pblock);
      solver->common.status = old_status;
      solver->common.numerical_rank = old_numerical_rank;
      solver->common.singular_col = old_singular_col;
      solver->common.noffdiag = old_noffdiag;
      (void)kls_rebuild_numeric_pinv(solver);
      return 0;
    }
    solver->numeric->Pnum[k1 + k] = solver->symbolic->P[k1 + local_row];
  }
  kls_record_fast_repaired_block_stats(solver, block, k1, nk, rejected_pivot,
                                       old_pblock, pblock);
  free(old_pblock);

  if (!kls_rebuild_numeric_pinv(solver)) {
    solver->common.status = old_status;
    solver->common.numerical_rank = old_numerical_rank;
    solver->common.singular_col = old_singular_col;
    solver->common.noffdiag = old_noffdiag;
    return 0;
  }

  int offdiag_refreshed = 0;
  if (tail_restart_used && !scaled) {
    offdiag_refreshed =
      kls_recompute_offdiag_suffix_from_pinv(solver, numeric_values,
                                             rejected_pivot);
  }
  if (!offdiag_refreshed) {
    offdiag_refreshed = kls_recompute_offdiag_from_pinv(solver, numeric_values);
  }
  if (!offdiag_refreshed) {
    solver->common.status = old_status;
    solver->common.numerical_rank = old_numerical_rank;
    solver->common.singular_col = old_singular_col;
    solver->common.noffdiag = old_noffdiag;
    return 0;
  }
  if (scaled && !kls_parallel_refactor_permute_scale(solver)) {
    solver->common.status = old_status;
    solver->common.numerical_rank = old_numerical_rank;
    solver->common.singular_col = old_singular_col;
    solver->common.noffdiag = old_noffdiag;
    return 0;
  }

  solver->numeric->lnz =
    solver->numeric->lnz - old_lnz_block + lnz_block;
  solver->numeric->unz =
    solver->numeric->unz - old_unz_block + unz_block;
  solver->numeric->max_lnz_block =
    solver->numeric->max_lnz_block < lnz_block
      ? lnz_block : solver->numeric->max_lnz_block;
  solver->numeric->max_unz_block =
    solver->numeric->max_unz_block < unz_block
      ? unz_block : solver->numeric->max_unz_block;
  if (solver->symbolic->Lnz[block] < 0.0) {
    solver->symbolic->Lnz[block] =
      (double)(lnz_block > unz_block ? lnz_block : unz_block);
  }
  UF_long offdiag = 0;
  for (UF_long k = 0; k < solver->n; ++k) {
    if (solver->numeric->Pnum[k] != solver->symbolic->P[k]) {
      offdiag++;
    }
  }
  solver->common.noffdiag = offdiag;
  solver->fast_block_restarts++;
  if (tail_restart_used) {
    solver->fast_tail_restarts++;
  }
  solver->stats.fast_block_restarts = solver->fast_block_restarts;
  solver->stats.fast_tail_restarts = solver->fast_tail_restarts;
  free_refactor_map(solver);
  free_refactor_schedule(solver);
  free_row_refactor_pattern(solver);
  free_refactor_lu_pointer_cache(solver);
  solver->common.status = TRILINOS_KLU_OK;
  solver->common.numerical_rank = KLS_KLU_EMPTY;
  solver->common.singular_col = KLS_KLU_EMPTY;
  return 1;
}

static int kls_refactor_map_is_eligible(const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      solver->col_ptr == NULL || solver->row_idx == NULL) {
    return 0;
  }
  if (solver->symbolic->R == NULL || solver->symbolic->Q == NULL ||
      solver->symbolic->R[0] != 0 ||
      solver->symbolic->R[solver->symbolic->nblocks] != solver->n) {
    return 0;
  }
  return solver->numeric->Pinv != NULL;
}

static int kls_build_refactor_map(kls_solver *solver) {
  if (!kls_refactor_map_is_eligible(solver)) {
    return 0;
  }
  if (solver->refactor_col_ptr != NULL && solver->refactor_row_idx != NULL &&
      solver->refactor_input_pos != NULL &&
      (solver->symbolic->nblocks == 1u ||
       (solver->refactor_block_start != NULL &&
        solver->refactor_col_block != NULL))) {
    return 1;
  }

  free_refactor_map(solver);
  UF_long *col_ptr =
    (UF_long *)malloc(((size_t)solver->n + 1u) * sizeof(*col_ptr));
  UF_long *row_idx =
    (UF_long *)malloc((size_t)solver->nnz * sizeof(*row_idx));
  UF_long *input_pos =
    (UF_long *)malloc((size_t)solver->nnz * sizeof(*input_pos));
  if (col_ptr == NULL || row_idx == NULL || input_pos == NULL) {
    free(col_ptr);
    free(row_idx);
    free(input_pos);
    return 0;
  }

  if (solver->symbolic->nblocks == 1u) {
    const UF_long *q = solver->symbolic->Q;
    const UF_long *pinv = solver->numeric->Pinv;
    UF_long count = 0;
    for (UF_long k = 0; k < solver->n; ++k) {
      col_ptr[k] = count;
      const UF_long oldcol = q[k];
      if (oldcol >= solver->n) {
        free(col_ptr);
        free(row_idx);
        free(input_pos);
        return 0;
      }
      for (UF_long p = solver->col_ptr[oldcol]; p < solver->col_ptr[oldcol + 1u]; ++p) {
        const UF_long oldrow = solver->row_idx[p];
        if (oldrow >= solver->n || count >= solver->nnz) {
          free(col_ptr);
          free(row_idx);
          free(input_pos);
          return 0;
        }
        const UF_long row = pinv[oldrow];
        if (row >= solver->n) {
          free(col_ptr);
          free(row_idx);
          free(input_pos);
          return 0;
        }
        row_idx[count] = row;
        input_pos[count] = p;
        count++;
      }
    }
    col_ptr[solver->n] = count;
    if (count != solver->nnz) {
      free(col_ptr);
      free(row_idx);
      free(input_pos);
      return 0;
    }
    solver->refactor_col_ptr = col_ptr;
    solver->refactor_row_idx = row_idx;
    solver->refactor_input_pos = input_pos;
    return 1;
  }

  UF_long *block_start =
    (UF_long *)malloc((size_t)solver->n * sizeof(*block_start));
  UF_long *col_block =
    (UF_long *)malloc((size_t)solver->n * sizeof(*col_block));
  if (block_start == NULL || col_block == NULL) {
    free(col_ptr);
    free(row_idx);
    free(input_pos);
    free(block_start);
    free(col_block);
    return 0;
  }

  const UF_long *q = solver->symbolic->Q;
  const UF_long *pinv = solver->numeric->Pinv;
  UF_long count = 0;
  UF_long block = 0;
  for (UF_long k = 0; k < solver->n; ++k) {
    while (block + 1u < solver->symbolic->nblocks &&
           solver->symbolic->R[block + 1u] <= k) {
      block++;
    }
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    col_ptr[k] = count;
    col_block[k] = block;
    const UF_long oldcol = q[k];
    if (oldcol >= solver->n) {
      free(col_ptr);
      free(row_idx);
      free(input_pos);
      free(block_start);
      free(col_block);
      return 0;
    }
    UF_long off_count = 0;
    UF_long block_count = 0;
    for (UF_long p = solver->col_ptr[oldcol]; p < solver->col_ptr[oldcol + 1u]; ++p) {
      const UF_long oldrow = solver->row_idx[p];
      if (oldrow >= solver->n) {
        free(col_ptr);
        free(row_idx);
        free(input_pos);
        free(block_start);
        free(col_block);
        return 0;
      }
      const UF_long row = pinv[oldrow];
      if (row >= solver->n) {
        free(col_ptr);
        free(row_idx);
        free(input_pos);
        free(block_start);
        free(col_block);
        return 0;
      }
      if (row < k1) {
        off_count++;
      } else if (row < k2) {
        block_count++;
      } else {
        free(col_ptr);
        free(row_idx);
        free(input_pos);
        free(block_start);
        free(col_block);
        return 0;
      }
    }
    if (solver->nnz - count < off_count + block_count) {
      free(col_ptr);
      free(row_idx);
      free(input_pos);
      free(block_start);
      free(col_block);
      return 0;
    }
    block_start[k] = count + off_count;
    UF_long off_pos = count;
    UF_long block_pos = block_start[k];
    for (UF_long p = solver->col_ptr[oldcol]; p < solver->col_ptr[oldcol + 1u]; ++p) {
      const UF_long row = pinv[solver->row_idx[p]];
      UF_long dst = 0;
      if (row < k1) {
        dst = off_pos++;
      } else {
        dst = block_pos++;
      }
      row_idx[dst] = row;
      input_pos[dst] = p;
    }
    count += off_count + block_count;
  }
  col_ptr[solver->n] = count;
  if (count != solver->nnz) {
    free(col_ptr);
    free(row_idx);
    free(input_pos);
    free(block_start);
    free(col_block);
    return 0;
  }

  solver->refactor_col_ptr = col_ptr;
  solver->refactor_row_idx = row_idx;
  solver->refactor_input_pos = input_pos;
  solver->refactor_block_start = block_start;
  solver->refactor_col_block = col_block;
  return 1;
}

static int kls_build_refactor_lu_pointer_cache(kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      solver->symbolic->R == NULL || solver->numeric->LUbx == NULL ||
      solver->numeric->Lip == NULL || solver->numeric->Llen == NULL ||
      solver->numeric->Uip == NULL || solver->numeric->Ulen == NULL ||
      solver->n == 0u) {
    return 0;
  }
  if (solver->refactor_lu_pointer_count == solver->n &&
      solver->refactor_l_indices != NULL &&
      solver->refactor_l_values != NULL &&
      solver->refactor_u_indices != NULL &&
      solver->refactor_u_values != NULL) {
    return 1;
  }

  free_refactor_lu_pointer_cache(solver);
  UF_long **l_indices =
    (UF_long **)calloc((size_t)solver->n, sizeof(*l_indices));
  double **l_values =
    (double **)calloc((size_t)solver->n, sizeof(*l_values));
  UF_long **u_indices =
    (UF_long **)calloc((size_t)solver->n, sizeof(*u_indices));
  double **u_values =
    (double **)calloc((size_t)solver->n, sizeof(*u_values));
  if (l_indices == NULL || l_values == NULL ||
      u_indices == NULL || u_values == NULL) {
    free(l_indices);
    free(l_values);
    free(u_indices);
    free(u_values);
    return 0;
  }

  for (UF_long block = 0; block < solver->symbolic->nblocks; ++block) {
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    const UF_long nk = k2 - k1;
    if (nk <= 1u) {
      continue;
    }
    double *lu = (double *)solver->numeric->LUbx[block];
    if (lu == NULL) {
      free(l_indices);
      free(l_values);
      free(u_indices);
      free(u_values);
      return 0;
    }
    const UF_long *lip = solver->numeric->Lip + k1;
    const UF_long *llen = solver->numeric->Llen + k1;
    const UF_long *uip = solver->numeric->Uip + k1;
    const UF_long *ulen = solver->numeric->Ulen + k1;
    for (UF_long k = 0; k < nk; ++k) {
      UF_long len = 0;
      kls_klu_get_pointer(lu, lip, llen, k,
                          &l_indices[k1 + k], &l_values[k1 + k], &len);
      kls_klu_get_pointer(lu, uip, ulen, k,
                          &u_indices[k1 + k], &u_values[k1 + k], &len);
    }
  }

  solver->refactor_l_indices = l_indices;
  solver->refactor_l_values = l_values;
  solver->refactor_u_indices = u_indices;
  solver->refactor_u_values = u_values;
  solver->refactor_lu_pointer_count = solver->n;
  return 1;
}

static int kls_row_refactor_pattern_is_eligible(const kls_solver *solver) {
  return solver != NULL && solver->symbolic != NULL &&
         solver->numeric != NULL && solver->symbolic->nblocks == 1u &&
         solver->common.scale <= 0 && solver->numeric->Rs == NULL &&
         solver->numeric->Udiag != NULL && solver->numeric->Xwork != NULL &&
         solver->refactor_lu_pointer_count == solver->n &&
         solver->refactor_col_ptr != NULL &&
         solver->refactor_row_idx != NULL &&
         solver->refactor_input_pos != NULL &&
         solver->refactor_l_indices != NULL &&
         solver->refactor_l_values != NULL &&
         solver->refactor_u_indices != NULL &&
         solver->refactor_u_values != NULL;
}

static int kls_build_row_refactor_pattern(kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      solver->common.scale > 0 || solver->numeric->Rs != NULL ||
      solver->symbolic->nblocks != 1u) {
    return 0;
  }
  if (solver->row_refactor_pattern_n == solver->n &&
      solver->row_refactor_l_ptr != NULL &&
      (solver->row_refactor_l_ptr[solver->n] == 0u ||
       solver->row_refactor_l_row_values != NULL) &&
      solver->row_refactor_u_ptr != NULL &&
      (solver->row_refactor_u_ptr[solver->n] == 0u ||
       solver->row_refactor_u_row_values != NULL) &&
      solver->row_refactor_input_ptr != NULL &&
      solver->row_refactor_successor_ptr != NULL &&
      (solver->row_refactor_successor_ptr[solver->n] == 0u ||
       solver->row_refactor_successor_rows != NULL) &&
      solver->row_refactor_group_dep_ptr != NULL &&
      (solver->row_refactor_group_dep_ptr[solver->row_refactor_group_count] == 0u ||
       solver->row_refactor_group_dep_rows != NULL)) {
    return 1;
  }
  if (!kls_build_refactor_map(solver) ||
      !kls_build_refactor_lu_pointer_cache(solver) ||
      !kls_row_refactor_pattern_is_eligible(solver)) {
    return 0;
  }

  free_row_refactor_pattern(solver);
  const UF_long n = solver->n;
  UF_long *l_ptr = (UF_long *)calloc((size_t)n + 1u, sizeof(*l_ptr));
  UF_long *u_ptr = (UF_long *)calloc((size_t)n + 1u, sizeof(*u_ptr));
  UF_long *input_ptr =
    (UF_long *)calloc((size_t)n + 1u, sizeof(*input_ptr));
  UF_long *next = (UF_long *)malloc((size_t)n * sizeof(*next));
  if (l_ptr == NULL || u_ptr == NULL || input_ptr == NULL ||
      (n > 0u && next == NULL)) {
    free(l_ptr);
    free(u_ptr);
    free(input_ptr);
    free(next);
    return 0;
  }

  for (UF_long j = 0; j < n; ++j) {
    UF_long *li = solver->refactor_l_indices[j];
    const UF_long len = solver->numeric->Llen[j];
    if (len > 0u && li == NULL) {
      free(l_ptr);
      free(u_ptr);
      free(input_ptr);
      free(next);
      return 0;
    }
    for (UF_long p = 0; p < len; ++p) {
      const UF_long row = li[p];
      if (row <= j || row >= n) {
        free(l_ptr);
        free(u_ptr);
        free(input_ptr);
        free(next);
        return 0;
      }
      l_ptr[row + 1u]++;
    }
  }
  for (UF_long k = 0; k < n; ++k) {
    UF_long *ui = solver->refactor_u_indices[k];
    const UF_long len = solver->numeric->Ulen[k];
    if (len > 0u && ui == NULL) {
      free(l_ptr);
      free(u_ptr);
      free(input_ptr);
      free(next);
      return 0;
    }
    for (UF_long p = 0; p < len; ++p) {
      const UF_long row = ui[p];
      if (row >= k || row >= n) {
        free(l_ptr);
        free(u_ptr);
        free(input_ptr);
        free(next);
        return 0;
      }
      u_ptr[row + 1u]++;
    }
  }
  for (UF_long k = 0; k < n; ++k) {
    for (UF_long p = solver->refactor_col_ptr[k];
         p < solver->refactor_col_ptr[k + 1u]; ++p) {
      const UF_long row = solver->refactor_row_idx[p];
      if (row >= n || solver->refactor_input_pos[p] >= solver->nnz) {
        free(l_ptr);
        free(u_ptr);
        free(input_ptr);
        free(next);
        return 0;
      }
      input_ptr[row + 1u]++;
    }
  }

  for (UF_long i = 0; i < n; ++i) {
    l_ptr[i + 1u] += l_ptr[i];
    u_ptr[i + 1u] += u_ptr[i];
    input_ptr[i + 1u] += input_ptr[i];
  }
  const UF_long l_nnz = l_ptr[n];
  const UF_long u_nnz = u_ptr[n];
  const UF_long input_nnz = input_ptr[n];
  UF_long *l_cols = l_nnz > 0u
    ? (UF_long *)malloc((size_t)l_nnz * sizeof(*l_cols)) : NULL;
  double **l_values = l_nnz > 0u
    ? (double **)malloc((size_t)l_nnz * sizeof(*l_values)) : NULL;
  double *l_row_values = l_nnz > 0u
    ? (double *)malloc((size_t)l_nnz * sizeof(*l_row_values)) : NULL;
  UF_long *u_cols = u_nnz > 0u
    ? (UF_long *)malloc((size_t)u_nnz * sizeof(*u_cols)) : NULL;
  double **u_values = u_nnz > 0u
    ? (double **)malloc((size_t)u_nnz * sizeof(*u_values)) : NULL;
  double *u_row_values = u_nnz > 0u
    ? (double *)malloc((size_t)u_nnz * sizeof(*u_row_values)) : NULL;
  UF_long *input_cols = input_nnz > 0u
    ? (UF_long *)malloc((size_t)input_nnz * sizeof(*input_cols)) : NULL;
  UF_long *input_pos = input_nnz > 0u
    ? (UF_long *)malloc((size_t)input_nnz * sizeof(*input_pos)) : NULL;
  if ((l_nnz > 0u &&
       (l_cols == NULL || l_values == NULL || l_row_values == NULL)) ||
      (u_nnz > 0u &&
       (u_cols == NULL || u_values == NULL || u_row_values == NULL)) ||
      (input_nnz > 0u && (input_cols == NULL || input_pos == NULL))) {
    free(l_ptr);
    free(u_ptr);
    free(input_ptr);
    free(next);
    free(l_cols);
    free(l_values);
    free(l_row_values);
    free(u_cols);
    free(u_values);
    free(u_row_values);
    free(input_cols);
    free(input_pos);
    return 0;
  }

  memcpy(next, l_ptr, (size_t)n * sizeof(*next));
  for (UF_long j = 0; j < n; ++j) {
    UF_long *li = solver->refactor_l_indices[j];
    double *lx = solver->refactor_l_values[j];
    const UF_long len = solver->numeric->Llen[j];
    for (UF_long p = 0; p < len; ++p) {
      const UF_long row = li[p];
      const UF_long dst = next[row]++;
      l_cols[dst] = j;
      l_values[dst] = lx + p;
    }
  }
  memcpy(next, u_ptr, (size_t)n * sizeof(*next));
  for (UF_long k = 0; k < n; ++k) {
    UF_long *ui = solver->refactor_u_indices[k];
    double *ux = solver->refactor_u_values[k];
    const UF_long len = solver->numeric->Ulen[k];
    for (UF_long p = 0; p < len; ++p) {
      const UF_long row = ui[p];
      const UF_long dst = next[row]++;
      u_cols[dst] = k;
      u_values[dst] = ux + p;
    }
  }
  memcpy(next, input_ptr, (size_t)n * sizeof(*next));
  for (UF_long k = 0; k < n; ++k) {
    for (UF_long p = solver->refactor_col_ptr[k];
         p < solver->refactor_col_ptr[k + 1u]; ++p) {
      const UF_long row = solver->refactor_row_idx[p];
      const UF_long dst = next[row]++;
      input_cols[dst] = k;
      input_pos[dst] = solver->refactor_input_pos[p];
    }
  }

  UF_long *successor_ptr =
    (UF_long *)calloc((size_t)n + 1u, sizeof(*successor_ptr));
  UF_long *successor_rows = l_nnz > 0u
    ? (UF_long *)malloc((size_t)l_nnz * sizeof(*successor_rows)) : NULL;
  if (successor_ptr == NULL || (l_nnz > 0u && successor_rows == NULL)) {
    free(l_ptr);
    free(u_ptr);
    free(input_ptr);
    free(next);
    free(l_cols);
    free(l_values);
    free(l_row_values);
    free(u_cols);
    free(u_values);
    free(u_row_values);
    free(input_cols);
    free(input_pos);
    free(successor_ptr);
    free(successor_rows);
    return 0;
  }
  for (UF_long row = 0; row < n; ++row) {
    for (UF_long p = l_ptr[row]; p < l_ptr[row + 1u]; ++p) {
      const UF_long dep = l_cols[p];
      if (dep >= row || dep >= n) {
        free(l_ptr);
        free(u_ptr);
        free(input_ptr);
        free(next);
        free(l_cols);
        free(l_values);
        free(l_row_values);
        free(u_cols);
        free(u_values);
        free(u_row_values);
        free(input_cols);
        free(input_pos);
        free(successor_ptr);
        free(successor_rows);
        return 0;
      }
      successor_ptr[dep + 1u]++;
    }
  }
  for (UF_long row = 0; row < n; ++row) {
    successor_ptr[row + 1u] += successor_ptr[row];
  }
  memcpy(next, successor_ptr, (size_t)n * sizeof(*next));
  for (UF_long row = 0; row < n; ++row) {
    for (UF_long p = l_ptr[row]; p < l_ptr[row + 1u]; ++p) {
      const UF_long dep = l_cols[p];
      successor_rows[next[dep]++] = row;
    }
  }
  free(next);

  solver->row_refactor_l_ptr = l_ptr;
  solver->row_refactor_l_cols = l_cols;
  solver->row_refactor_l_values = l_values;
  solver->row_refactor_l_row_values = l_row_values;
  solver->row_refactor_u_ptr = u_ptr;
  solver->row_refactor_u_cols = u_cols;
  solver->row_refactor_u_values = u_values;
  solver->row_refactor_u_row_values = u_row_values;
  solver->row_refactor_input_ptr = input_ptr;
  solver->row_refactor_input_cols = input_cols;
  solver->row_refactor_input_pos = input_pos;
  solver->row_refactor_successor_ptr = successor_ptr;
  solver->row_refactor_successor_rows = successor_rows;
  solver->row_refactor_pattern_n = n;

  UF_long *levels = (UF_long *)calloc((size_t)n, sizeof(*levels));
  if (n > 0u && levels == NULL) {
    free_row_refactor_pattern(solver);
    return 0;
  }
  UF_long max_level = 0;
  for (UF_long i = 0; i < n; ++i) {
    UF_long level = 0;
    for (UF_long p = l_ptr[i]; p < l_ptr[i + 1u]; ++p) {
      const UF_long dep = l_cols[p];
      if (dep >= i) {
        free(levels);
        free_row_refactor_pattern(solver);
        return 0;
      }
      const UF_long dep_level = levels[dep] + 1u;
      if (dep_level > level) {
        level = dep_level;
      }
    }
    levels[i] = level;
    if (level > max_level) {
      max_level = level;
    }
  }
  const UF_long level_count = n > 0u ? max_level + 1u : 0u;
  UF_long *level_ptr =
    (UF_long *)calloc((size_t)level_count + 1u, sizeof(*level_ptr));
  UF_long *level_rows = n > 0u
    ? (UF_long *)malloc((size_t)n * sizeof(*level_rows)) : NULL;
  if (level_ptr == NULL || (n > 0u && level_rows == NULL)) {
    free(levels);
    free(level_ptr);
    free(level_rows);
    free_row_refactor_pattern(solver);
    return 0;
  }
  for (UF_long i = 0; i < n; ++i) {
    level_ptr[levels[i] + 1u]++;
  }
  for (UF_long level = 0; level < level_count; ++level) {
    level_ptr[level + 1u] += level_ptr[level];
  }
  UF_long *write_pos = level_count > 0u
    ? (UF_long *)malloc((size_t)level_count * sizeof(*write_pos)) : NULL;
  if (level_count > 0u && write_pos == NULL) {
    free(levels);
    free(level_ptr);
    free(level_rows);
    free_row_refactor_pattern(solver);
    return 0;
  }
  memcpy(write_pos, level_ptr, (size_t)level_count * sizeof(*write_pos));
  for (UF_long i = 0; i < n; ++i) {
    level_rows[write_pos[levels[i]]++] = i;
  }
  free(write_pos);
  free(levels);
  UF_long max_width = 0;
  for (UF_long level = 0; level < level_count; ++level) {
    const UF_long width = level_ptr[level + 1u] - level_ptr[level];
    if (width > max_width) {
      max_width = width;
    }
  }
  solver->row_refactor_level_ptr = level_ptr;
  solver->row_refactor_level_rows = level_rows;
  solver->row_refactor_level_count = level_count;
  solver->row_refactor_level_max_width = max_width;

  UF_long *row_group = n > 0u
    ? (UF_long *)malloc((size_t)n * sizeof(*row_group)) : NULL;
  UF_long *group_ptr = n > 0u
    ? (UF_long *)malloc(((size_t)n + 1u) * sizeof(*group_ptr)) : NULL;
  if ((n > 0u && row_group == NULL) || (n > 0u && group_ptr == NULL)) {
    free(row_group);
    free(group_ptr);
    free_row_refactor_pattern(solver);
    return 0;
  }
  UF_long group_count = 0;
  UF_long row = 0;
  while (row < n) {
    group_ptr[group_count] = row;
    UF_long end = row + 1u;
    while (end < n) {
      const UF_long prev = end - 1u;
      const UF_long prev_begin = u_ptr[prev];
      const UF_long prev_end = u_ptr[prev + 1u];
      const UF_long curr_begin = u_ptr[end];
      const UF_long curr_end = u_ptr[end + 1u];
      const UF_long prev_len = prev_end - prev_begin;
      const UF_long curr_len = curr_end - curr_begin;
      int extends = 0;
      if (prev_len == curr_len + 1u &&
          prev_len > 0u &&
          u_cols[prev_begin] == end) {
        if (curr_len == 0u ||
            memcmp(u_cols + prev_begin + 1u,
                   u_cols + curr_begin,
                   (size_t)curr_len * sizeof(*u_cols)) == 0) {
          extends = 1;
        }
      }
      if (!extends) {
        break;
      }
      end++;
    }
    for (UF_long i = row; i < end; ++i) {
      row_group[i] = group_count;
    }
    group_count++;
    row = end;
  }
  group_ptr[group_count] = n;

  UF_long *l_internal_ptr = n > 0u
    ? (UF_long *)calloc((size_t)n, sizeof(*l_internal_ptr)) : NULL;
  UF_long *group_trailing_len = group_count > 0u
    ? (UF_long *)calloc((size_t)group_count, sizeof(*group_trailing_len))
    : NULL;
  unsigned char *group_dense = group_count > 0u
    ? (unsigned char *)calloc((size_t)group_count, sizeof(*group_dense))
    : NULL;
  unsigned char *group_kind = group_count > 0u
    ? (unsigned char *)calloc((size_t)group_count, sizeof(*group_kind))
    : NULL;
  if ((n > 0u && l_internal_ptr == NULL) ||
      (group_count > 0u &&
       (group_trailing_len == NULL || group_dense == NULL ||
        group_kind == NULL))) {
    free(row_group);
    free(group_ptr);
    free(l_internal_ptr);
    free(group_trailing_len);
    free(group_dense);
    free(group_kind);
    free_row_refactor_pattern(solver);
    return 0;
  }

  for (UF_long g = 0; g < group_count; ++g) {
    group_kind[g] = KLS_ROW_REFACTOR_GROUP_SINGLE;
  }

  UF_long *group_dep_ptr =
    (UF_long *)calloc((size_t)group_count + 1u, sizeof(*group_dep_ptr));
  UF_long *group_dep_rows = NULL;
  UF_long *dep_marks = n > 0u
    ? (UF_long *)malloc((size_t)n * sizeof(*dep_marks)) : NULL;
  if (group_dep_ptr == NULL || (n > 0u && dep_marks == NULL)) {
    free(row_group);
    free(group_ptr);
    free(l_internal_ptr);
    free(group_trailing_len);
    free(group_dense);
    free(group_kind);
    free(group_dep_ptr);
    free(dep_marks);
    free_row_refactor_pattern(solver);
    return 0;
  }
  for (UF_long i = 0; i < n; ++i) {
    dep_marks[i] = KLS_KLU_EMPTY;
  }
  for (UF_long g = 0; g < group_count; ++g) {
    for (UF_long i = group_ptr[g]; i < group_ptr[g + 1u]; ++i) {
      for (UF_long p = l_ptr[i]; p < l_ptr[i + 1u]; ++p) {
        const UF_long dep = l_cols[p];
        if (dep >= n || row_group[dep] == g || dep_marks[dep] == g) {
          continue;
        }
        dep_marks[dep] = g;
        group_dep_ptr[g + 1u]++;
      }
    }
  }
  for (UF_long g = 0; g < group_count; ++g) {
    group_dep_ptr[g + 1u] += group_dep_ptr[g];
  }
  const UF_long group_dep_count = group_dep_ptr[group_count];
  group_dep_rows = group_dep_count > 0u
    ? (UF_long *)malloc((size_t)group_dep_count * sizeof(*group_dep_rows))
    : NULL;
  UF_long *group_dep_write = group_count > 0u
    ? (UF_long *)malloc((size_t)group_count * sizeof(*group_dep_write))
    : NULL;
  if ((group_dep_count > 0u && group_dep_rows == NULL) ||
      (group_count > 0u && group_dep_write == NULL)) {
    free(row_group);
    free(group_ptr);
    free(l_internal_ptr);
    free(group_trailing_len);
    free(group_dense);
    free(group_kind);
    free(group_dep_ptr);
    free(group_dep_rows);
    free(group_dep_write);
    free(dep_marks);
    free_row_refactor_pattern(solver);
    return 0;
  }
  for (UF_long i = 0; i < n; ++i) {
    dep_marks[i] = KLS_KLU_EMPTY;
  }
  if (group_count > 0u) {
    memcpy(group_dep_write, group_dep_ptr,
           (size_t)group_count * sizeof(*group_dep_write));
  }
  for (UF_long g = 0; g < group_count; ++g) {
    for (UF_long i = group_ptr[g]; i < group_ptr[g + 1u]; ++i) {
      for (UF_long p = l_ptr[i]; p < l_ptr[i + 1u]; ++p) {
        const UF_long dep = l_cols[p];
        if (dep >= n || row_group[dep] == g || dep_marks[dep] == g) {
          continue;
        }
        dep_marks[dep] = g;
        group_dep_rows[group_dep_write[g]++] = dep;
      }
    }
  }
  free(group_dep_write);
  free(dep_marks);

  UF_long segment_count = 0;
  UF_long segment_rows = 0;
  UF_long segment_max_width = 0;
  double segment_dense_entries = 0.0;
  double segment_trailing_entries = 0.0;
  UF_long dense_segment_count = 0;
  UF_long dense_segment_rows = 0;
  UF_long dense_segment_max_width = 0;
  double dense_segment_dense_entries = 0.0;
  double dense_segment_trailing_entries = 0.0;
  for (UF_long g = 0; g < group_count; ++g) {
    const UF_long row_begin = group_ptr[g];
    const UF_long row_end = group_ptr[g + 1u];
    const UF_long width = row_end - row_begin;
    if (width <= 1u) {
      continue;
    }
    group_kind[g] = KLS_ROW_REFACTOR_GROUP_GENERIC;
    const UF_long trailing_len =
      u_ptr[row_end] - u_ptr[row_end - 1u];
    const double dense_entries =
      (double)width * (double)(width - 1u) * 0.5;
    const double trailing_entries =
      (double)width * (double)trailing_len;
    group_trailing_len[g] = trailing_len;
    segment_count++;
    segment_rows += width;
    if (width > segment_max_width) {
      segment_max_width = width;
    }
    segment_dense_entries += dense_entries;
    segment_trailing_entries += trailing_entries;

    if (dense_entries + trailing_entries <
        KLS_ROW_REFACTOR_DENSE_MIN_WORK) {
      continue;
    }

    int dense = 1;
    const UF_long trailing_begin = u_ptr[row_end - 1u];
    const UF_long *trailing_cols = u_cols + trailing_begin;
    for (UF_long i = row_begin; dense && i < row_end; ++i) {
      UF_long p = l_ptr[i];
      const UF_long l_end = l_ptr[i + 1u];
      while (p < l_end && l_cols[p] < row_begin) {
        p++;
      }
      l_internal_ptr[i] = p;
      for (UF_long dep = row_begin; dep < i; ++dep, ++p) {
        if (p >= l_end || l_cols[p] != dep) {
          dense = 0;
          break;
        }
      }
      if (dense && p != l_end) {
        dense = 0;
      }

      const UF_long dense_len = row_end - i - 1u;
      const UF_long u_begin = u_ptr[i];
      const UF_long u_end = u_ptr[i + 1u];
      if (dense && u_end - u_begin != dense_len + trailing_len) {
        dense = 0;
      }
      for (UF_long offset = 0; dense && offset < dense_len; ++offset) {
        if (u_cols[u_begin + offset] != i + 1u + offset) {
          dense = 0;
        }
      }
      for (UF_long offset = 0; dense && offset < trailing_len; ++offset) {
        if (u_cols[u_begin + dense_len + offset] != trailing_cols[offset]) {
          dense = 0;
        }
      }
    }
    if (!dense) {
      continue;
    }
    group_dense[g] = 1u;
    group_kind[g] = KLS_ROW_REFACTOR_GROUP_DENSE;
    dense_segment_count++;
    dense_segment_rows += width;
    if (width > dense_segment_max_width) {
      dense_segment_max_width = width;
    }
    dense_segment_dense_entries += dense_entries;
    dense_segment_trailing_entries += trailing_entries;
  }

  UF_long *group_levels = group_count > 0u
    ? (UF_long *)calloc((size_t)group_count, sizeof(*group_levels)) : NULL;
  if (group_count > 0u && group_levels == NULL) {
    free(row_group);
    free(group_ptr);
    free(l_internal_ptr);
    free(group_trailing_len);
    free(group_dense);
    free(group_kind);
    free(group_dep_ptr);
    free(group_dep_rows);
    free_row_refactor_pattern(solver);
    return 0;
  }
  UF_long group_max_level = 0;
  for (UF_long g = 0; g < group_count; ++g) {
    UF_long group_level = 0;
    for (UF_long i = group_ptr[g]; i < group_ptr[g + 1u]; ++i) {
      for (UF_long p = l_ptr[i]; p < l_ptr[i + 1u]; ++p) {
        const UF_long dep_group = row_group[l_cols[p]];
        if (dep_group == g) {
          continue;
        }
        if (dep_group > g) {
          free(row_group);
          free(group_ptr);
          free(group_levels);
          free(l_internal_ptr);
          free(group_trailing_len);
          free(group_dense);
          free(group_kind);
          free(group_dep_ptr);
          free(group_dep_rows);
          free_row_refactor_pattern(solver);
          return 0;
        }
        const UF_long dep_level = group_levels[dep_group] + 1u;
        if (dep_level > group_level) {
          group_level = dep_level;
        }
      }
    }
    group_levels[g] = group_level;
    if (group_level > group_max_level) {
      group_max_level = group_level;
    }
  }
  const UF_long group_level_count =
    group_count > 0u ? group_max_level + 1u : 0u;
  UF_long *group_level_ptr =
    (UF_long *)calloc((size_t)group_level_count + 1u,
                      sizeof(*group_level_ptr));
  UF_long *level_groups = group_count > 0u
    ? (UF_long *)malloc((size_t)group_count * sizeof(*level_groups)) : NULL;
  if (group_level_ptr == NULL ||
      (group_count > 0u && level_groups == NULL)) {
    free(row_group);
    free(group_ptr);
    free(group_levels);
    free(l_internal_ptr);
    free(group_trailing_len);
    free(group_dense);
    free(group_kind);
    free(group_dep_ptr);
    free(group_dep_rows);
    free(group_level_ptr);
    free(level_groups);
    free_row_refactor_pattern(solver);
    return 0;
  }
  for (UF_long g = 0; g < group_count; ++g) {
    group_level_ptr[group_levels[g] + 1u]++;
  }
  for (UF_long level = 0; level < group_level_count; ++level) {
    group_level_ptr[level + 1u] += group_level_ptr[level];
  }
  UF_long *group_write_pos = group_level_count > 0u
    ? (UF_long *)malloc((size_t)group_level_count *
                        sizeof(*group_write_pos)) : NULL;
  if (group_level_count > 0u && group_write_pos == NULL) {
    free(row_group);
    free(group_ptr);
    free(group_levels);
    free(l_internal_ptr);
    free(group_trailing_len);
    free(group_dense);
    free(group_kind);
    free(group_dep_ptr);
    free(group_dep_rows);
    free(group_level_ptr);
    free(level_groups);
    free_row_refactor_pattern(solver);
    return 0;
  }
  memcpy(group_write_pos, group_level_ptr,
         (size_t)group_level_count * sizeof(*group_write_pos));
  for (UF_long g = 0; g < group_count; ++g) {
    level_groups[group_write_pos[group_levels[g]]++] = g;
  }
  UF_long group_max_width = 0;
  for (UF_long level = 0; level < group_level_count; ++level) {
    const UF_long width = group_level_ptr[level + 1u] - group_level_ptr[level];
    if (width > group_max_width) {
      group_max_width = width;
    }
  }
  free(group_write_pos);
  free(group_levels);
  free(row_group);

  solver->row_refactor_group_ptr = group_ptr;
  solver->row_refactor_group_dep_ptr = group_dep_ptr;
  solver->row_refactor_group_dep_rows = group_dep_rows;
  solver->row_refactor_group_count = group_count;
  solver->row_refactor_group_level_ptr = group_level_ptr;
  solver->row_refactor_level_groups = level_groups;
  solver->row_refactor_l_internal_ptr = l_internal_ptr;
  solver->row_refactor_group_trailing_len = group_trailing_len;
  solver->row_refactor_group_dense = group_dense;
  solver->row_refactor_group_kind = group_kind;
  solver->row_refactor_level_count = group_level_count;
  solver->row_refactor_level_max_width = group_max_width;
  solver->row_refactor_cluster_level_count = group_level_count;
  solver->row_refactor_pipeline_group_count = 0;
  solver->row_refactor_pipeline_row_count = 0;
  solver->row_refactor_pipeline_work = 0.0;
  solver->row_refactor_segment_count = segment_count;
  solver->row_refactor_segment_rows = segment_rows;
  solver->row_refactor_segment_max_width = segment_max_width;
  solver->row_refactor_segment_dense_entries = segment_dense_entries;
  solver->row_refactor_segment_trailing_entries = segment_trailing_entries;
  solver->row_refactor_dense_segment_count = dense_segment_count;
  solver->row_refactor_dense_segment_rows = dense_segment_rows;
  solver->row_refactor_dense_segment_max_width = dense_segment_max_width;
  solver->row_refactor_dense_segment_dense_entries =
    dense_segment_dense_entries;
  solver->row_refactor_dense_segment_trailing_entries =
    dense_segment_trailing_entries;
  return 1;
}

static double kls_row_refactor_group_work(const kls_solver *solver,
                                          UF_long group) {
  if (solver == NULL || group >= solver->row_refactor_group_count ||
      solver->row_refactor_group_ptr == NULL ||
      solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_u_ptr == NULL ||
      solver->row_refactor_input_ptr == NULL ||
      solver->row_refactor_l_cols == NULL) {
    return 1.0;
  }
  const UF_long row_begin = solver->row_refactor_group_ptr[group];
  const UF_long row_end = solver->row_refactor_group_ptr[group + 1u];
  if (row_begin >= row_end || row_end > solver->n) {
    return 1.0;
  }

  double work = 1.0;
  for (UF_long row = row_begin; row < row_end; ++row) {
    work += (double)(solver->row_refactor_input_ptr[row + 1u] -
                     solver->row_refactor_input_ptr[row]);
    work += (double)(solver->row_refactor_u_ptr[row + 1u] -
                     solver->row_refactor_u_ptr[row]);
    for (UF_long p = solver->row_refactor_l_ptr[row];
         p < solver->row_refactor_l_ptr[row + 1u]; ++p) {
      const UF_long dep = solver->row_refactor_l_cols[p];
      if (dep < solver->n) {
        work += 1.0 + (double)(solver->row_refactor_u_ptr[dep + 1u] -
                               solver->row_refactor_u_ptr[dep]);
      }
    }
  }

  if (solver->row_refactor_group_dense != NULL &&
      solver->row_refactor_group_dense[group]) {
    const double width = (double)(row_end - row_begin);
    const double trailing = solver->row_refactor_group_trailing_len != NULL
      ? (double)solver->row_refactor_group_trailing_len[group] : 0.0;
    work += width * width * width / 3.0;
    work += width * width * trailing * 0.5;
  }
  return work;
}

static UF_long kls_row_refactor_choose_cluster_levels(
  const kls_solver *solver,
  int thread_count) {
  if (solver == NULL || thread_count <= 1 ||
      solver->row_refactor_group_level_ptr == NULL ||
      solver->row_refactor_level_count == 0u) {
    return solver != NULL ? solver->row_refactor_level_count : 0u;
  }

  const double width_limit = 2.0 * (double)thread_count;
  for (UF_long level = 0; level < solver->row_refactor_level_count; ++level) {
    const UF_long width =
      solver->row_refactor_group_level_ptr[level + 1u] -
      solver->row_refactor_group_level_ptr[level];
    if ((double)width < width_limit) {
      return level;
    }
  }
  return solver->row_refactor_level_count;
}

static void kls_row_refactor_record_pipeline_scope(kls_solver *solver,
                                                   UF_long cluster_levels) {
  if (solver == NULL) {
    return;
  }
  if (cluster_levels > solver->row_refactor_level_count) {
    cluster_levels = solver->row_refactor_level_count;
  }
  solver->row_refactor_cluster_level_count = cluster_levels;
  solver->row_refactor_pipeline_group_count = 0;
  solver->row_refactor_pipeline_row_count = 0;
  solver->row_refactor_pipeline_work = 0.0;
  if (solver->row_refactor_group_level_ptr == NULL ||
      solver->row_refactor_level_groups == NULL ||
      solver->row_refactor_group_ptr == NULL ||
      cluster_levels >= solver->row_refactor_level_count) {
    return;
  }

  const UF_long begin = solver->row_refactor_group_level_ptr[cluster_levels];
  for (UF_long pos = begin; pos < solver->row_refactor_group_count; ++pos) {
    const UF_long group = solver->row_refactor_level_groups[pos];
    if (group >= solver->row_refactor_group_count) {
      continue;
    }
    const UF_long row_begin = solver->row_refactor_group_ptr[group];
    const UF_long row_end = solver->row_refactor_group_ptr[group + 1u];
    if (row_begin <= row_end && row_end <= solver->n) {
      solver->row_refactor_pipeline_row_count += row_end - row_begin;
    }
    solver->row_refactor_pipeline_group_count++;
    solver->row_refactor_pipeline_work +=
      kls_row_refactor_group_work(solver, group);
  }
}

static int kls_build_row_refactor_group_thread_slices(kls_solver *solver,
                                                      int thread_count) {
  if (solver == NULL || thread_count <= 1 ||
      solver->row_refactor_group_level_ptr == NULL ||
      solver->row_refactor_level_groups == NULL ||
      solver->row_refactor_level_count == 0u) {
    return 1;
  }
  if (solver->row_refactor_group_level_thread_ptr != NULL &&
      solver->row_refactor_group_level_thread_count == thread_count) {
    return 1;
  }

  free(solver->row_refactor_group_level_thread_ptr);
  solver->row_refactor_group_level_thread_ptr = NULL;
  solver->row_refactor_group_level_thread_count = 0;

  const size_t levels = (size_t)solver->row_refactor_level_count;
  const size_t stride = (size_t)thread_count + 1u;
  if (levels > SIZE_MAX / stride ||
      levels * stride > SIZE_MAX / sizeof(UF_long)) {
    return 0;
  }
  UF_long *thread_ptr =
    (UF_long *)malloc(levels * stride * sizeof(*thread_ptr));
  if (thread_ptr == NULL) {
    return 0;
  }

  for (UF_long level = 0; level < solver->row_refactor_level_count; ++level) {
    const UF_long begin = solver->row_refactor_group_level_ptr[level];
    const UF_long end = solver->row_refactor_group_level_ptr[level + 1u];
    UF_long *parts = thread_ptr + (size_t)level * stride;
    parts[0] = begin;
    parts[thread_count] = end;

    double level_work = 0.0;
    for (UF_long pos = begin; pos < end; ++pos) {
      const UF_long group = solver->row_refactor_level_groups[pos];
      level_work += kls_row_refactor_group_work(solver, group);
    }

    UF_long pos = begin;
    double prefix_work = 0.0;
    for (int t = 1; t < thread_count; ++t) {
      const double target =
        level_work * (double)t / (double)thread_count;
      while (pos < end && prefix_work < target) {
        const UF_long group = solver->row_refactor_level_groups[pos];
        prefix_work += kls_row_refactor_group_work(solver, group);
        pos++;
      }
      parts[t] = pos;
    }
  }

  solver->row_refactor_group_level_thread_ptr = thread_ptr;
  solver->row_refactor_group_level_thread_count = thread_count;
  return 1;
}

static int kls_row_refactor_env_enabled(void) {
  const char *value = getenv("KLS_ENABLE_ROW_REFACTOR");
  return value != NULL && value[0] != '\0' &&
         !(value[0] == '0' && value[1] == '\0');
}

static int kls_checked_row_refactor_env_enabled(void) {
  const char *value = getenv("KLS_ENABLE_CHECKED_ROW_REFACTOR");
  return value != NULL && value[0] != '\0' &&
         !(value[0] == '0' && value[1] == '\0');
}

static int kls_row_refactor_should_defer_value_scatter(const kls_solver *solver) {
  return solver != NULL && solver->row_refactor_dense_segment_count > 0u;
}

static int kls_scatter_row_refactor_l_prefix_values(kls_solver *solver,
                                                    UF_long end_row) {
  if (solver == NULL || solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_l_values == NULL) {
    return 0;
  }
  if (end_row > solver->n) {
    return 0;
  }
  const UF_long l_nnz = solver->row_refactor_l_ptr[end_row];
  if (l_nnz == 0u) {
    return 1;
  }
  if (solver->row_refactor_l_row_values == NULL) {
    return 0;
  }
  for (UF_long p = 0; p < l_nnz; ++p) {
    *solver->row_refactor_l_values[p] =
      solver->row_refactor_l_row_values[p];
  }
  return 1;
}

static int kls_scatter_row_refactor_u_prefix_values(kls_solver *solver,
                                                    UF_long end_row) {
  if (solver == NULL || solver->row_refactor_u_ptr == NULL ||
      solver->row_refactor_u_values == NULL) {
    return 0;
  }
  if (end_row > solver->n) {
    return 0;
  }
  const UF_long u_nnz = solver->row_refactor_u_ptr[end_row];
  if (u_nnz == 0u) {
    return 1;
  }
  if (solver->row_refactor_u_row_values == NULL) {
    return 0;
  }
  for (UF_long p = 0; p < u_nnz; ++p) {
    *solver->row_refactor_u_values[p] =
      solver->row_refactor_u_row_values[p];
  }
  return 1;
}

static int kls_scatter_row_refactor_prefix_values(kls_solver *solver,
                                                  UF_long end_row) {
  return kls_scatter_row_refactor_l_prefix_values(solver, end_row) &&
         kls_scatter_row_refactor_u_prefix_values(solver, end_row);
}

static int kls_scatter_row_refactor_l_values(kls_solver *solver) {
  return solver != NULL &&
         kls_scatter_row_refactor_l_prefix_values(solver, solver->n);
}

static int kls_scatter_row_refactor_u_values(kls_solver *solver) {
  return solver != NULL &&
         kls_scatter_row_refactor_u_prefix_values(solver, solver->n);
}

static int kls_row_refactor_multiplier_rejects(double lij, double tol) {
  if (tol <= DBL_MIN) {
    return 0;
  }
  const double lij_abs = fabs(lij);
  return !isfinite(lij_abs) || lij_abs * tol > 1.0 + 1.0e-12;
}

static int kls_single_block_row_refactor(kls_solver *solver,
                                         double *numeric_values,
                                         int check_pivots) {
  if (numeric_values == NULL || solver == NULL || solver->symbolic == NULL ||
      solver->numeric == NULL || solver->symbolic->nblocks != 1u ||
      solver->common.scale > 0 || solver->numeric->Rs != NULL ||
      solver->numeric->Udiag == NULL || solver->numeric->Xwork == NULL ||
      !kls_build_row_refactor_pattern(solver)) {
    return -1;
  }

  trilinos_klu_l_common *common = &solver->common;
  const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  trilinos_klu_l_numeric *numeric = solver->numeric;
  double *x = (double *)numeric->Xwork;
  double *udiag = (double *)numeric->Udiag;
  const UF_long n = solver->n;

  memset(x, 0, (size_t)n * sizeof(*x));
  common->status = TRILINOS_KLU_OK;
  common->numerical_rank = KLS_KLU_EMPTY;
  common->singular_col = KLS_KLU_EMPTY;
  common->nrealloc = 0;
  const int defer_value_scatter =
    kls_row_refactor_should_defer_value_scatter(solver);
  kls_row_refactor_record_pipeline_scope(solver,
                                         solver->row_refactor_level_count);
  kls_record_row_refactor_run(solver, check_pivots, 0);

  for (UF_long i = 0; i < n; ++i) {
    for (UF_long p = solver->row_refactor_input_ptr[i];
         p < solver->row_refactor_input_ptr[i + 1u]; ++p) {
      const UF_long col = solver->row_refactor_input_cols[p];
      x[col] = numeric_values[solver->row_refactor_input_pos[p]];
    }

    for (UF_long p = solver->row_refactor_l_ptr[i];
         p < solver->row_refactor_l_ptr[i + 1u]; ++p) {
      const UF_long dep = solver->row_refactor_l_cols[p];
      const double lij = x[dep] / udiag[dep];
      const double lij_abs = check_pivots ? fabs(lij) : 0.0;
      if (check_pivots &&
          kls_row_refactor_multiplier_rejects(lij, common->tol)) {
        const double candidate_abs = fabs(x[dep]);
        x[dep] = 0.0;
        memset(x, 0, (size_t)n * sizeof(*x));
        if (defer_value_scatter &&
            !kls_scatter_row_refactor_prefix_values(solver, dep)) {
          common->status = TRILINOS_KLU_INVALID;
          return 0;
        }
        solver->fast_reject_refresh_state = KLS_FAST_REJECT_REFRESH_PREFIX;
        kls_record_fast_reject_detail(
          solver, dep,
          symbolic->Q != NULL ? symbolic->Q[dep] : KLS_KLU_EMPTY,
          i, lij_abs, fabs(udiag[dep]), candidate_abs);
        kls_fill_fast_reject_row_tail_candidate_stats(
          solver, numeric_values, dep, x);
        common->status = TRILINOS_KLU_OK;
        return 0;
      }
      if (defer_value_scatter) {
        solver->row_refactor_l_row_values[p] = lij;
      } else {
        *solver->row_refactor_l_values[p] = lij;
      }
      x[dep] = 0.0;
      const UF_long u_begin = solver->row_refactor_u_ptr[dep];
      const UF_long u_end = solver->row_refactor_u_ptr[dep + 1u];
      const UF_long *u_cols = solver->row_refactor_u_cols + u_begin;
      const double *u_values =
        solver->row_refactor_u_row_values + u_begin;
      for (UF_long offset = 0; offset < u_end - u_begin; ++offset) {
        x[u_cols[offset]] -= lij * u_values[offset];
      }
    }

    const double pivot = x[i];
    x[i] = 0.0;
    if (pivot == 0.0) {
      common->status = TRILINOS_KLU_SINGULAR;
      if (common->numerical_rank == KLS_KLU_EMPTY) {
        common->numerical_rank = i;
        common->singular_col = symbolic->Q[i];
      }
      if (common->halt_if_singular) {
        memset(x, 0, (size_t)n * sizeof(*x));
        return 0;
      }
    }
    udiag[i] = pivot;

    const UF_long u_begin = solver->row_refactor_u_ptr[i];
    const UF_long u_end = solver->row_refactor_u_ptr[i + 1u];
    double *row_u_values = solver->row_refactor_u_row_values + u_begin;
    if (defer_value_scatter) {
      for (UF_long p = u_begin; p < u_end; ++p) {
        const UF_long offset = p - u_begin;
        const UF_long col = solver->row_refactor_u_cols[p];
        row_u_values[offset] = x[col];
        x[col] = 0.0;
      }
    } else {
      double **klu_u_values = solver->row_refactor_u_values + u_begin;
      for (UF_long p = u_begin; p < u_end; ++p) {
        const UF_long offset = p - u_begin;
        const UF_long col = solver->row_refactor_u_cols[p];
        const double value = x[col];
        row_u_values[offset] = value;
        *klu_u_values[offset] = value;
        x[col] = 0.0;
      }
    }
    for (UF_long p = solver->row_refactor_input_ptr[i];
         p < solver->row_refactor_input_ptr[i + 1u]; ++p) {
      x[solver->row_refactor_input_cols[p]] = 0.0;
    }
  }
  if (defer_value_scatter &&
      (!kls_scatter_row_refactor_l_values(solver) ||
       !kls_scatter_row_refactor_u_values(solver))) {
    common->status = TRILINOS_KLU_INVALID;
    return 0;
  }
  return 1;
}

static int kls_parallel_row_refactor_rejects_multiplier(
  kls_egraph_refactor_worker *worker,
  UF_long row,
  UF_long dep,
  double candidate,
  double lij) {
  if (worker == NULL || worker->shared == NULL ||
      !worker->shared->check_pivots) {
    return 0;
  }
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  if (solver == NULL || solver->symbolic == NULL ||
      solver->symbolic->Q == NULL || solver->numeric == NULL ||
      solver->numeric->Udiag == NULL || row >= solver->n ||
      dep >= solver->n ||
      !kls_row_refactor_multiplier_rejects(lij, solver->common.tol)) {
    return 0;
  }
  const double *udiag = (const double *)solver->numeric->Udiag;
  kls_egraph_refactor_record_reject(
    shared, dep, solver->symbolic->Q[dep], row, fabs(lij),
    fabs(udiag[dep]), fabs(candidate));
  return 1;
}

static int kls_parallel_row_refactor_process_row(
  kls_egraph_refactor_worker *worker,
  UF_long row,
  int wait_for_dependencies) {
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  trilinos_klu_l_numeric *numeric = solver->numeric;
  const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  double *x = worker->x;
  double *udiag = (double *)numeric->Udiag;

  for (UF_long p = solver->row_refactor_input_ptr[row];
       p < solver->row_refactor_input_ptr[row + 1u]; ++p) {
    const UF_long col = solver->row_refactor_input_cols[p];
    x[col] = shared->values[solver->row_refactor_input_pos[p]];
  }

  for (UF_long p = solver->row_refactor_l_ptr[row];
       p < solver->row_refactor_l_ptr[row + 1u]; ++p) {
    const UF_long dep = solver->row_refactor_l_cols[p];
    if (wait_for_dependencies &&
        !kls_egraph_refactor_wait_done(shared, dep)) {
      x[dep] = 0.0;
      return 0;
    }
    const double candidate = x[dep];
    const double lij = candidate / udiag[dep];
    if (kls_parallel_row_refactor_rejects_multiplier(
          worker, row, dep, candidate, lij)) {
      x[dep] = 0.0;
      return 0;
    }
    if (shared->row_refactor_defer_value_scatter) {
      solver->row_refactor_l_row_values[p] = lij;
    } else {
      *solver->row_refactor_l_values[p] = lij;
    }
    x[dep] = 0.0;
    const UF_long u_begin = solver->row_refactor_u_ptr[dep];
    const UF_long u_end = solver->row_refactor_u_ptr[dep + 1u];
    const UF_long *u_cols = solver->row_refactor_u_cols + u_begin;
    const double *u_values = solver->row_refactor_u_row_values + u_begin;
    for (UF_long offset = 0; offset < u_end - u_begin; ++offset) {
      x[u_cols[offset]] -= lij * u_values[offset];
    }
  }

  const double pivot = x[row];
  x[row] = 0.0;
  if (pivot == 0.0) {
    kls_egraph_refactor_record_singular(shared, row, symbolic->Q[row]);
    if (solver->common.halt_if_singular) {
      memset(x, 0, (size_t)solver->n * sizeof(*x));
      return 0;
    }
  }
  udiag[row] = pivot;

  const UF_long u_begin = solver->row_refactor_u_ptr[row];
  const UF_long u_end = solver->row_refactor_u_ptr[row + 1u];
  double *row_u_values = solver->row_refactor_u_row_values + u_begin;
  if (shared->row_refactor_defer_value_scatter) {
    for (UF_long p = u_begin; p < u_end; ++p) {
      const UF_long offset = p - u_begin;
      const UF_long col = solver->row_refactor_u_cols[p];
      row_u_values[offset] = x[col];
      x[col] = 0.0;
    }
  } else {
    double **klu_u_values = solver->row_refactor_u_values + u_begin;
    for (UF_long p = u_begin; p < u_end; ++p) {
      const UF_long offset = p - u_begin;
      const UF_long col = solver->row_refactor_u_cols[p];
      const double value = x[col];
      row_u_values[offset] = value;
      *klu_u_values[offset] = value;
      x[col] = 0.0;
    }
  }
  for (UF_long p = solver->row_refactor_input_ptr[row];
       p < solver->row_refactor_input_ptr[row + 1u]; ++p) {
    x[solver->row_refactor_input_cols[p]] = 0.0;
  }
  return 1;
}

static double *kls_egraph_worker_segment_panel(
  kls_egraph_refactor_worker *worker,
  UF_long entry_count) {
  if (entry_count == 0u) {
    return NULL;
  }
  if (worker == NULL ||
      entry_count > (UF_long)(SIZE_MAX / sizeof(*worker->segment_panel))) {
    return NULL;
  }
  if (worker->segment_panel != NULL &&
      worker->segment_panel_size >= entry_count) {
    return worker->segment_panel;
  }

  double *panel =
    (double *)realloc(worker->segment_panel,
                      (size_t)entry_count * sizeof(*panel));
  if (panel == NULL) {
    return NULL;
  }
  worker->segment_panel = panel;
  worker->segment_panel_size = entry_count;
  return panel;
}

static int kls_parallel_row_refactor_process_dense_group_native(
  kls_egraph_refactor_worker *worker,
  UF_long group,
  UF_long row_begin,
  UF_long row_end,
  UF_long trailing_len,
  const UF_long *trailing_cols,
  int wait_for_dependencies) {
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  if (!shared->row_refactor_defer_value_scatter ||
      solver->row_refactor_group_dense == NULL ||
      solver->row_refactor_group_trailing_len == NULL ||
      solver->row_refactor_l_internal_ptr == NULL ||
      solver->row_refactor_l_row_values == NULL ||
      solver->row_refactor_u_row_values == NULL ||
      !solver->row_refactor_group_dense[group] ||
      solver->row_refactor_group_trailing_len[group] != trailing_len) {
    return -1;
  }

  trilinos_klu_l_numeric *numeric = solver->numeric;
  const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  double *x = worker->x;
  double *udiag = (double *)numeric->Udiag;

  for (UF_long row = row_begin; row < row_end; ++row) {
    for (UF_long p = solver->row_refactor_input_ptr[row];
         p < solver->row_refactor_input_ptr[row + 1u]; ++p) {
      const UF_long col = solver->row_refactor_input_cols[p];
      x[col] = shared->values[solver->row_refactor_input_pos[p]];
    }

    UF_long lp = solver->row_refactor_l_ptr[row];
    const UF_long l_internal = solver->row_refactor_l_internal_ptr[row];
    for (; lp < l_internal; ++lp) {
      const UF_long dep = solver->row_refactor_l_cols[lp];
      if (wait_for_dependencies &&
          !kls_egraph_refactor_wait_done(shared, dep)) {
        x[dep] = 0.0;
        return 0;
      }
      const double candidate = x[dep];
      const double lij = candidate / udiag[dep];
      if (kls_parallel_row_refactor_rejects_multiplier(
            worker, row, dep, candidate, lij)) {
        x[dep] = 0.0;
        return 0;
      }
      solver->row_refactor_l_row_values[lp] = lij;
      x[dep] = 0.0;

      const UF_long u_begin = solver->row_refactor_u_ptr[dep];
      const UF_long u_end = solver->row_refactor_u_ptr[dep + 1u];
      const UF_long *u_cols = solver->row_refactor_u_cols + u_begin;
      const double *u_values =
        solver->row_refactor_u_row_values + u_begin;
      for (UF_long offset = 0; offset < u_end - u_begin; ++offset) {
        x[u_cols[offset]] -= lij * u_values[offset];
      }
    }

    const UF_long l_end = solver->row_refactor_l_ptr[row + 1u];
    double *row_l_values = solver->row_refactor_l_row_values + l_internal;
    for (UF_long dep = row_begin; dep < row; ++dep, ++lp) {
      if (lp >= l_end || solver->row_refactor_l_cols[lp] != dep) {
        kls_egraph_refactor_record_invalid(shared);
        return 0;
      }
      row_l_values[dep - row_begin] = x[dep];
      x[dep] = 0.0;
    }
    if (lp != l_end) {
      kls_egraph_refactor_record_invalid(shared);
      return 0;
    }

    udiag[row] = x[row];
    x[row] = 0.0;

    const UF_long dense_len = row_end - row - 1u;
    const UF_long u_begin = solver->row_refactor_u_ptr[row];
    const UF_long u_end = solver->row_refactor_u_ptr[row + 1u];
    if (u_end - u_begin != dense_len + trailing_len) {
      kls_egraph_refactor_record_invalid(shared);
      return 0;
    }
    double *row_u_values = solver->row_refactor_u_row_values + u_begin;
    for (UF_long offset = 0; offset < dense_len; ++offset) {
      const UF_long col = row + 1u + offset;
      row_u_values[offset] = x[col];
      x[col] = 0.0;
    }
    for (UF_long offset = 0; offset < trailing_len; ++offset) {
      const UF_long col = trailing_cols[offset];
      row_u_values[dense_len + offset] = x[col];
      x[col] = 0.0;
    }
    for (UF_long p = solver->row_refactor_input_ptr[row];
         p < solver->row_refactor_input_ptr[row + 1u]; ++p) {
      x[solver->row_refactor_input_cols[p]] = 0.0;
    }
  }

  for (UF_long row = row_begin; row < row_end; ++row) {
    const UF_long l_internal = solver->row_refactor_l_internal_ptr[row];
    double *row_l_values = solver->row_refactor_l_row_values + l_internal;
    const UF_long row_dense_len = row_end - row - 1u;
    const UF_long row_u_begin = solver->row_refactor_u_ptr[row];
    double *row_u_values =
      solver->row_refactor_u_row_values + row_u_begin;
    double *row_trailing_values = row_u_values + row_dense_len;

    for (UF_long dep = row_begin; dep < row; ++dep) {
      const UF_long local_dep = dep - row_begin;
      const double candidate = row_l_values[local_dep];
      const double lij = candidate / udiag[dep];
      if (kls_parallel_row_refactor_rejects_multiplier(
            worker, row, dep, candidate, lij)) {
        return 0;
      }
      row_l_values[local_dep] = lij;

      const UF_long dep_dense_len = row_end - dep - 1u;
      const UF_long dep_u_begin = solver->row_refactor_u_ptr[dep];
      const UF_long dep_u_end = solver->row_refactor_u_ptr[dep + 1u];
      if (dep_u_end - dep_u_begin != dep_dense_len + trailing_len) {
        kls_egraph_refactor_record_invalid(shared);
        return 0;
      }
      const double *dep_u_values =
        solver->row_refactor_u_row_values + dep_u_begin;
      const double *dep_trailing_values = dep_u_values + dep_dense_len;

      for (UF_long target = dep + 1u; target < row; ++target) {
        row_l_values[target - row_begin] -=
          lij * dep_u_values[target - dep - 1u];
      }
      udiag[row] -= lij * dep_u_values[row - dep - 1u];
      for (UF_long target = row + 1u; target < row_end; ++target) {
        row_u_values[target - row - 1u] -=
          lij * dep_u_values[target - dep - 1u];
      }
      for (UF_long offset = 0; offset < trailing_len; ++offset) {
        row_trailing_values[offset] -= lij * dep_trailing_values[offset];
      }
    }

    const double pivot = udiag[row];
    if (pivot == 0.0) {
      kls_egraph_refactor_record_singular(shared, row, symbolic->Q[row]);
      if (solver->common.halt_if_singular) {
        memset(x, 0, (size_t)solver->n * sizeof(*x));
        return 0;
      }
    }
  }
  return 1;
}

static int kls_parallel_row_refactor_process_dense_group(
  kls_egraph_refactor_worker *worker,
  UF_long group,
  UF_long row_begin,
  UF_long row_end,
  UF_long trailing_len,
  const UF_long *trailing_cols,
  int wait_for_dependencies) {
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  if (solver->row_refactor_group_dense == NULL ||
      solver->row_refactor_group_trailing_len == NULL ||
      solver->row_refactor_l_internal_ptr == NULL ||
      !solver->row_refactor_group_dense[group] ||
      solver->row_refactor_group_trailing_len[group] != trailing_len) {
    return -1;
  }

  const int native_status =
    kls_parallel_row_refactor_process_dense_group_native(
      worker, group, row_begin, row_end, trailing_len, trailing_cols,
      wait_for_dependencies);
  if (native_status >= 0) {
    return native_status;
  }

  trilinos_klu_l_numeric *numeric = solver->numeric;
  const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  double *x = worker->x;
  double *udiag = (double *)numeric->Udiag;
  const UF_long width = row_end - row_begin;
  const UF_long max_panel_entries =
    (UF_long)(SIZE_MAX / sizeof(*worker->segment_panel));
  if (width > max_panel_entries / width) {
    return -1;
  }
  const UF_long dense_panel_entries = width * width;
  UF_long trailing_panel_entries = 0;
  if (trailing_len > 0u) {
    if (width > (max_panel_entries - dense_panel_entries) / trailing_len) {
      return -1;
    }
    trailing_panel_entries = width * trailing_len;
  }
  double *dense_panel =
    kls_egraph_worker_segment_panel(
      worker, dense_panel_entries + trailing_panel_entries);
  if (dense_panel == NULL) {
    return -1;
  }
  double *trailing_panel = NULL;
  if (trailing_len > 0u) {
    trailing_panel = dense_panel + dense_panel_entries;
  }

  for (UF_long row = row_begin; row < row_end; ++row) {
    const UF_long local_row = row - row_begin;
    double *row_dense_panel = dense_panel + local_row * width;
    for (UF_long p = solver->row_refactor_input_ptr[row];
         p < solver->row_refactor_input_ptr[row + 1u]; ++p) {
      const UF_long col = solver->row_refactor_input_cols[p];
      x[col] = shared->values[solver->row_refactor_input_pos[p]];
    }

    UF_long lp = solver->row_refactor_l_ptr[row];
    const UF_long l_internal = solver->row_refactor_l_internal_ptr[row];
    for (; lp < l_internal; ++lp) {
      const UF_long dep = solver->row_refactor_l_cols[lp];
      if (wait_for_dependencies &&
          !kls_egraph_refactor_wait_done(shared, dep)) {
        x[dep] = 0.0;
        return 0;
      }
      const double candidate = x[dep];
      const double lij = candidate / udiag[dep];
      if (kls_parallel_row_refactor_rejects_multiplier(
            worker, row, dep, candidate, lij)) {
        x[dep] = 0.0;
        return 0;
      }
      if (shared->row_refactor_defer_value_scatter) {
        solver->row_refactor_l_row_values[lp] = lij;
      } else {
        *solver->row_refactor_l_values[lp] = lij;
      }
      x[dep] = 0.0;
      const UF_long u_begin = solver->row_refactor_u_ptr[dep];
      const UF_long u_end = solver->row_refactor_u_ptr[dep + 1u];
      const UF_long *u_cols = solver->row_refactor_u_cols + u_begin;
      const double *u_values =
        solver->row_refactor_u_row_values + u_begin;
      for (UF_long offset = 0; offset < u_end - u_begin; ++offset) {
        x[u_cols[offset]] -= lij * u_values[offset];
      }
    }

    for (UF_long dep = row_begin; dep < row; ++dep, ++lp) {
      row_dense_panel[dep - row_begin] = x[dep];
      x[dep] = 0.0;
    }

    udiag[row] = x[row];
    x[row] = 0.0;

    const UF_long dense_len = row_end - row - 1u;
    double *row_panel = trailing_panel != NULL
      ? trailing_panel + local_row * trailing_len : NULL;
    for (UF_long offset = 0; offset < dense_len; ++offset) {
      const UF_long col = row + 1u + offset;
      row_dense_panel[local_row + 1u + offset] = x[col];
      x[col] = 0.0;
    }
    for (UF_long offset = 0; offset < trailing_len; ++offset) {
      const UF_long col = trailing_cols[offset];
      row_panel[offset] = x[col];
      x[col] = 0.0;
    }
    for (UF_long p = solver->row_refactor_input_ptr[row];
         p < solver->row_refactor_input_ptr[row + 1u]; ++p) {
      x[solver->row_refactor_input_cols[p]] = 0.0;
    }
  }

  for (UF_long row = row_begin; row < row_end; ++row) {
    const UF_long local_row = row - row_begin;
    double *row_dense_panel = dense_panel + local_row * width;

    for (UF_long dep = row_begin; dep < row; ++dep) {
      const UF_long local_dep = dep - row_begin;
      const double candidate = row_dense_panel[local_dep];
      const double lij = candidate / udiag[dep];
      if (kls_parallel_row_refactor_rejects_multiplier(
            worker, row, dep, candidate, lij)) {
        return 0;
      }
      row_dense_panel[local_dep] = lij;

      double *dep_dense_panel = dense_panel + local_dep * width;

      for (UF_long target = local_dep + 1u; target < local_row; ++target) {
        row_dense_panel[target] -= lij * dep_dense_panel[target];
      }
      udiag[row] -= lij * dep_dense_panel[local_row];
      for (UF_long target = local_row + 1u; target < width; ++target) {
        row_dense_panel[target] -= lij * dep_dense_panel[target];
      }
    }

    const double pivot = udiag[row];
    if (pivot == 0.0) {
      kls_egraph_refactor_record_singular(shared, row, symbolic->Q[row]);
      if (solver->common.halt_if_singular) {
        memset(x, 0, (size_t)solver->n * sizeof(*x));
        return 0;
      }
    }
  }

  /* Keep the dense segment and shared trailing panel in compact row-major
     scratch so this path can move to persistent row/segment storage without
     changing semantics. */
  for (UF_long row = row_begin; row < row_end; ++row) {
    const UF_long local_row = row - row_begin;
    double *row_dense_panel = dense_panel + local_row * width;
    double *row_panel = trailing_panel != NULL
      ? trailing_panel + local_row * trailing_len : NULL;

    if (trailing_len > 0u) {
      for (UF_long dep = row_begin; dep < row; ++dep) {
        const double lij = row_dense_panel[dep - row_begin];
        double *dep_panel =
          trailing_panel + (dep - row_begin) * trailing_len;

        for (UF_long offset = 0; offset < trailing_len; ++offset) {
          row_panel[offset] -= lij * dep_panel[offset];
        }
      }
    }

    const UF_long l_dense_begin =
      solver->row_refactor_l_internal_ptr[row];
    if (shared->row_refactor_defer_value_scatter) {
      double *row_l_values =
        solver->row_refactor_l_row_values + l_dense_begin;
      for (UF_long dep = row_begin; dep < row; ++dep) {
        row_l_values[dep - row_begin] =
          row_dense_panel[dep - row_begin];
      }
    } else {
      double **row_l_values =
        solver->row_refactor_l_values + l_dense_begin;
      for (UF_long dep = row_begin; dep < row; ++dep) {
        *row_l_values[dep - row_begin] =
          row_dense_panel[dep - row_begin];
      }
    }

    const UF_long row_dense_len = row_end - row - 1u;
    const UF_long row_u_begin = solver->row_refactor_u_ptr[row];
    double *row_u_row_values =
      solver->row_refactor_u_row_values + row_u_begin;
    if (shared->row_refactor_defer_value_scatter) {
      for (UF_long offset = 0; offset < row_dense_len; ++offset) {
        row_u_row_values[offset] =
          row_dense_panel[local_row + 1u + offset];
      }
      if (trailing_len > 0u) {
        for (UF_long offset = 0; offset < trailing_len; ++offset) {
          row_u_row_values[row_dense_len + offset] = row_panel[offset];
        }
      }
    } else {
      double **row_u_values = solver->row_refactor_u_values + row_u_begin;
      for (UF_long offset = 0; offset < row_dense_len; ++offset) {
        const double value = row_dense_panel[local_row + 1u + offset];
        row_u_row_values[offset] = value;
        *row_u_values[offset] = value;
      }
      if (trailing_len > 0u) {
        for (UF_long offset = 0; offset < trailing_len; ++offset) {
          const double value = row_panel[offset];
          row_u_row_values[row_dense_len + offset] = value;
          *row_u_values[row_dense_len + offset] = value;
        }
      }
    }
  }
  return 1;
}

static kls_row_refactor_group_kind kls_row_refactor_group_kind_for(
  const kls_solver *solver,
  UF_long group,
  UF_long width) {
  if (width <= 1u) {
    return KLS_ROW_REFACTOR_GROUP_SINGLE;
  }
  if (solver != NULL && solver->row_refactor_group_kind != NULL &&
      group < solver->row_refactor_group_count) {
    const unsigned char kind = solver->row_refactor_group_kind[group];
    if (kind == KLS_ROW_REFACTOR_GROUP_DENSE) {
      return KLS_ROW_REFACTOR_GROUP_DENSE;
    }
    if (kind == KLS_ROW_REFACTOR_GROUP_GENERIC) {
      return KLS_ROW_REFACTOR_GROUP_GENERIC;
    }
  }
  if (solver != NULL && solver->row_refactor_group_dense != NULL &&
      group < solver->row_refactor_group_count &&
      solver->row_refactor_group_dense[group]) {
    return KLS_ROW_REFACTOR_GROUP_DENSE;
  }
  return KLS_ROW_REFACTOR_GROUP_GENERIC;
}

static int kls_parallel_row_refactor_process_group(
  kls_egraph_refactor_worker *worker,
  UF_long group,
  int wait_for_dependencies) {
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  if (group >= solver->row_refactor_group_count ||
      solver->row_refactor_group_ptr == NULL) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }
  const UF_long row_begin = solver->row_refactor_group_ptr[group];
  const UF_long row_end = solver->row_refactor_group_ptr[group + 1u];
  if (row_begin >= row_end || row_end > solver->n) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }
  const UF_long width = row_end - row_begin;
  const kls_row_refactor_group_kind group_kind =
    kls_row_refactor_group_kind_for(solver, group, width);
  if (group_kind == KLS_ROW_REFACTOR_GROUP_SINGLE) {
    return kls_parallel_row_refactor_process_row(worker, row_begin,
                                                wait_for_dependencies);
  }

  trilinos_klu_l_numeric *numeric = solver->numeric;
  const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  double *x = worker->x;
  double *udiag = (double *)numeric->Udiag;
  const UF_long trailing_begin = solver->row_refactor_u_ptr[row_end - 1u];
  const UF_long trailing_len =
    solver->row_refactor_group_trailing_len != NULL
      ? solver->row_refactor_group_trailing_len[group]
      : solver->row_refactor_u_ptr[row_end] - trailing_begin;
  const UF_long *trailing_cols =
    solver->row_refactor_u_cols + trailing_begin;

  if (group_kind == KLS_ROW_REFACTOR_GROUP_DENSE) {
    const int dense_status =
      kls_parallel_row_refactor_process_dense_group(
        worker, group, row_begin, row_end, trailing_len, trailing_cols,
        wait_for_dependencies);
    if (dense_status >= 0) {
      return dense_status;
    }
  }

  for (UF_long row = row_begin; row < row_end; ++row) {
    for (UF_long p = solver->row_refactor_input_ptr[row];
         p < solver->row_refactor_input_ptr[row + 1u]; ++p) {
      const UF_long col = solver->row_refactor_input_cols[p];
      x[col] = shared->values[solver->row_refactor_input_pos[p]];
    }

    for (UF_long p = solver->row_refactor_l_ptr[row];
         p < solver->row_refactor_l_ptr[row + 1u]; ++p) {
      const UF_long dep = solver->row_refactor_l_cols[p];
      if (wait_for_dependencies && dep < row_begin &&
          !kls_egraph_refactor_wait_done(shared, dep)) {
        x[dep] = 0.0;
        return 0;
      }
      const double candidate = x[dep];
      const double lij = candidate / udiag[dep];
      if (kls_parallel_row_refactor_rejects_multiplier(
            worker, row, dep, candidate, lij)) {
        x[dep] = 0.0;
        return 0;
      }
      if (shared->row_refactor_defer_value_scatter) {
        solver->row_refactor_l_row_values[p] = lij;
      } else {
        *solver->row_refactor_l_values[p] = lij;
      }
      x[dep] = 0.0;

      if (dep >= row_begin) {
        const UF_long dense_len = row_end - dep - 1u;
        const UF_long u_begin = solver->row_refactor_u_ptr[dep];
        const UF_long u_end = solver->row_refactor_u_ptr[dep + 1u];
        if (u_end - u_begin != dense_len + trailing_len) {
          kls_egraph_refactor_record_invalid(shared);
          return 0;
        }
        const double *u_values =
          solver->row_refactor_u_row_values + u_begin;
        for (UF_long offset = 0; offset < dense_len; ++offset) {
          x[dep + 1u + offset] -= lij * u_values[offset];
        }
        for (UF_long offset = 0; offset < trailing_len; ++offset) {
          x[trailing_cols[offset]] -=
            lij * u_values[dense_len + offset];
        }
      } else {
        const UF_long u_begin = solver->row_refactor_u_ptr[dep];
        const UF_long u_end = solver->row_refactor_u_ptr[dep + 1u];
        const UF_long *u_cols = solver->row_refactor_u_cols + u_begin;
        const double *u_values =
          solver->row_refactor_u_row_values + u_begin;
        for (UF_long offset = 0; offset < u_end - u_begin; ++offset) {
          x[u_cols[offset]] -= lij * u_values[offset];
        }
      }
    }

    const double pivot = x[row];
    x[row] = 0.0;
    if (pivot == 0.0) {
      kls_egraph_refactor_record_singular(shared, row, symbolic->Q[row]);
      if (solver->common.halt_if_singular) {
        memset(x, 0, (size_t)solver->n * sizeof(*x));
        return 0;
      }
    }
    udiag[row] = pivot;

    const UF_long dense_len = row_end - row - 1u;
    const UF_long u_begin = solver->row_refactor_u_ptr[row];
    const UF_long u_end = solver->row_refactor_u_ptr[row + 1u];
    if (u_end - u_begin != dense_len + trailing_len) {
      kls_egraph_refactor_record_invalid(shared);
      return 0;
    }
    double *row_u_values = solver->row_refactor_u_row_values + u_begin;
    if (shared->row_refactor_defer_value_scatter) {
      for (UF_long offset = 0; offset < dense_len; ++offset) {
        const UF_long col = row + 1u + offset;
        row_u_values[offset] = x[col];
        x[col] = 0.0;
      }
      for (UF_long offset = 0; offset < trailing_len; ++offset) {
        const UF_long col = trailing_cols[offset];
        row_u_values[dense_len + offset] = x[col];
        x[col] = 0.0;
      }
    } else {
      double **klu_u_values = solver->row_refactor_u_values + u_begin;
      for (UF_long offset = 0; offset < dense_len; ++offset) {
        const UF_long col = row + 1u + offset;
        const double value = x[col];
        row_u_values[offset] = value;
        *klu_u_values[offset] = value;
        x[col] = 0.0;
      }
      for (UF_long offset = 0; offset < trailing_len; ++offset) {
        const UF_long col = trailing_cols[offset];
        const double value = x[col];
        row_u_values[dense_len + offset] = value;
        *klu_u_values[dense_len + offset] = value;
        x[col] = 0.0;
      }
    }
    for (UF_long p = solver->row_refactor_input_ptr[row];
         p < solver->row_refactor_input_ptr[row + 1u]; ++p) {
      x[solver->row_refactor_input_cols[p]] = 0.0;
    }
  }
  return 1;
}

static void kls_parallel_row_refactor_mark_group_done(
  kls_egraph_refactor_shared *shared,
  const kls_solver *solver,
  UF_long group) {
  if (shared == NULL || solver == NULL ||
      shared->pipeline_done == NULL ||
      solver->row_refactor_group_ptr == NULL ||
      group >= solver->row_refactor_group_count) {
    return;
  }
  const UF_long row_begin = solver->row_refactor_group_ptr[group];
  const UF_long row_end = solver->row_refactor_group_ptr[group + 1u];
  if (row_begin >= row_end || row_end > solver->n) {
    return;
  }
  for (UF_long row = row_begin; row < row_end; ++row) {
    kls_egraph_refactor_mark_done(shared, row);
  }
}

static int kls_parallel_row_refactor_wait_group_dependencies(
  kls_egraph_refactor_worker *worker,
  UF_long group) {
  if (worker == NULL || worker->shared == NULL) {
    return 0;
  }
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  if (shared->pipeline_done == NULL) {
    return 1;
  }
  if (solver == NULL || solver->row_refactor_group_dep_ptr == NULL ||
      group >= solver->row_refactor_group_count) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }
  const UF_long begin = solver->row_refactor_group_dep_ptr[group];
  const UF_long end = solver->row_refactor_group_dep_ptr[group + 1u];
  if (begin == end) {
    return 1;
  }
  if (solver->row_refactor_group_dep_rows == NULL) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }
  for (UF_long pos = begin; pos < end; ++pos) {
    if (!kls_egraph_refactor_wait_done(
          shared, solver->row_refactor_group_dep_rows[pos])) {
      return 0;
    }
  }
  return 1;
}

static int kls_row_refactor_group_prefers_external_prewait(
  const kls_solver *solver,
  UF_long group) {
  if (solver == NULL || solver->row_refactor_group_ptr == NULL ||
      solver->row_refactor_group_dep_ptr == NULL ||
      group >= solver->row_refactor_group_count) {
    return 0;
  }
  if (solver->row_refactor_dense_segment_count > 0u) {
    return 0;
  }
  if (solver->row_refactor_group_dep_ptr[group] ==
      solver->row_refactor_group_dep_ptr[group + 1u]) {
    return 0;
  }
  const UF_long row_begin = solver->row_refactor_group_ptr[group];
  const UF_long row_end = solver->row_refactor_group_ptr[group + 1u];
  if (row_begin >= row_end || row_end > solver->n) {
    return 0;
  }
  const UF_long width = row_end - row_begin;
  return kls_row_refactor_group_kind_for(solver, group, width) !=
         KLS_ROW_REFACTOR_GROUP_DENSE;
}

static void kls_row_refactor_worker_run(kls_egraph_refactor_worker *worker) {
  if (worker == NULL || worker->shared == NULL) {
    return;
  }
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  if (solver == NULL) {
    kls_egraph_refactor_record_invalid(shared);
    return;
  }

  if (solver->row_refactor_group_level_ptr != NULL &&
      solver->row_refactor_level_groups != NULL &&
      solver->row_refactor_group_ptr != NULL) {
    UF_long cluster_levels = shared->cluster_level_count;
    if (cluster_levels > solver->row_refactor_level_count) {
      cluster_levels = solver->row_refactor_level_count;
    }
    for (UF_long level = 0; level < cluster_levels; ++level) {
      if (!kls_egraph_refactor_should_stop(shared)) {
        const UF_long begin = solver->row_refactor_group_level_ptr[level];
        const UF_long end = solver->row_refactor_group_level_ptr[level + 1u];
        if (solver->row_refactor_group_level_thread_ptr != NULL &&
            solver->row_refactor_group_level_thread_count ==
              shared->thread_count) {
          const UF_long *parts =
            solver->row_refactor_group_level_thread_ptr +
            level * (UF_long)(shared->thread_count + 1);
          for (UF_long pos = parts[worker->tid];
               pos < parts[worker->tid + 1]; ++pos) {
            const UF_long group = solver->row_refactor_level_groups[pos];
            if (!kls_parallel_row_refactor_process_group(worker, group, 0)) {
              break;
            }
            kls_parallel_row_refactor_mark_group_done(shared, solver, group);
          }
        } else {
          for (UF_long pos = begin + (UF_long)worker->tid;
               pos < end; pos += (UF_long)shared->thread_count) {
            const UF_long group = solver->row_refactor_level_groups[pos];
            if (!kls_parallel_row_refactor_process_group(worker, group, 0)) {
              break;
            }
            kls_parallel_row_refactor_mark_group_done(shared, solver, group);
          }
        }
      }

      (void)pthread_barrier_wait(&shared->barrier);
      if (kls_egraph_refactor_should_stop(shared)) {
        break;
      }
    }
    if (shared->pipeline_done != NULL &&
        cluster_levels < solver->row_refactor_level_count &&
        !kls_egraph_refactor_should_stop(shared)) {
      for (;;) {
        if (kls_egraph_refactor_should_stop(shared)) {
          break;
        }
        const UF_long pos =
          (UF_long)atomic_fetch_add_explicit(&shared->next_pipeline_pos, 1ul,
                                             memory_order_relaxed);
        if (pos >= shared->pipeline_pos_end) {
          break;
        }
        const UF_long group = solver->row_refactor_level_groups[pos];
        const int prewait =
          kls_row_refactor_group_prefers_external_prewait(solver, group);
        if ((prewait &&
             !kls_parallel_row_refactor_wait_group_dependencies(worker,
                                                                group)) ||
            !kls_parallel_row_refactor_process_group(worker, group,
                                                     prewait ? 0 : 1)) {
          break;
        }
        kls_parallel_row_refactor_mark_group_done(shared, solver, group);
      }
    }
    return;
  }

  if (solver->row_refactor_level_ptr == NULL ||
      solver->row_refactor_level_rows == NULL) {
    kls_egraph_refactor_record_invalid(shared);
    return;
  }

  for (UF_long level = 0; level < solver->row_refactor_level_count; ++level) {
    if (!kls_egraph_refactor_should_stop(shared)) {
      const UF_long begin = solver->row_refactor_level_ptr[level];
      const UF_long end = solver->row_refactor_level_ptr[level + 1u];
      for (UF_long pos = begin + (UF_long)worker->tid;
           pos < end; pos += (UF_long)shared->thread_count) {
        const UF_long row = solver->row_refactor_level_rows[pos];
        if (row >= solver->n ||
            !kls_parallel_row_refactor_process_row(worker, row, 0)) {
          break;
        }
        kls_egraph_refactor_mark_done(shared, row);
      }
    }

    (void)pthread_barrier_wait(&shared->barrier);
    if (kls_egraph_refactor_should_stop(shared)) {
      break;
    }
  }
}

static int kls_single_block_parallel_row_refactor(kls_solver *solver,
                                                 double *numeric_values,
                                                 int check_pivots) {
  if (numeric_values == NULL || solver == NULL || solver->symbolic == NULL ||
      solver->numeric == NULL || solver->symbolic->nblocks != 1u ||
      solver->common.scale > 0 || solver->numeric->Rs != NULL ||
      solver->numeric->Udiag == NULL || solver->options.threads <= 1 ||
      !kls_build_row_refactor_pattern(solver) ||
      solver->row_refactor_level_max_width < 2u) {
    return -1;
  }

  int thread_count = solver->options.threads;
  if ((UF_long)thread_count > solver->row_refactor_level_max_width) {
    thread_count = (int)solver->row_refactor_level_max_width;
  }
  if (thread_count < 2) {
    return -1;
  }
  if (!kls_build_row_refactor_group_thread_slices(solver, thread_count)) {
    return -1;
  }
  UF_long cluster_levels =
    kls_row_refactor_choose_cluster_levels(solver, thread_count);
  int use_group_pipeline = 0;
  UF_long pipeline_group_begin = 0;
  if (solver->row_refactor_group_level_ptr != NULL &&
      solver->row_refactor_level_groups != NULL &&
      solver->row_refactor_group_ptr != NULL &&
      cluster_levels < solver->row_refactor_level_count) {
    pipeline_group_begin = solver->row_refactor_group_level_ptr[cluster_levels];
    use_group_pipeline = pipeline_group_begin < solver->row_refactor_group_count;
  }
  if (!use_group_pipeline) {
    cluster_levels = solver->row_refactor_level_count;
  }
  kls_row_refactor_record_pipeline_scope(solver, cluster_levels);
  atomic_uint *pipeline_done = NULL;
  unsigned int pipeline_generation = 0;
  if (check_pivots || use_group_pipeline) {
    pipeline_done =
      ensure_egraph_pipeline_done(solver, &pipeline_generation);
    if (pipeline_done == NULL) {
      return -1;
    }
  }

  double **scratch =
    ensure_egraph_worker_scratch(solver, thread_count, solver->n);
  if (scratch == NULL) {
    return -1;
  }
  kls_egraph_refactor_pool *pool =
    ensure_egraph_refactor_pool(solver, thread_count);
  if (pool == NULL) {
    return -1;
  }

  trilinos_klu_l_common *common = &solver->common;
  common->status = TRILINOS_KLU_OK;
  common->numerical_rank = KLS_KLU_EMPTY;
  common->singular_col = KLS_KLU_EMPTY;
  common->nrealloc = 0;

  kls_egraph_refactor_shared *shared = &pool->shared;
  pthread_mutex_lock(&shared->lock);
  if (pool->active_workers != 0) {
    pthread_mutex_unlock(&shared->lock);
    return -1;
  }

  shared->solver = solver;
  shared->values = numeric_values;
  shared->rs = NULL;
  shared->check_pivots = check_pivots;
  shared->scale = 0;
  shared->kernel = KLS_EGRAPH_REFACTOR_KERNEL_GENERIC;
  shared->thread_count = thread_count;
  shared->row_refactor_mode = 1;
  shared->row_refactor_defer_value_scatter =
    kls_row_refactor_should_defer_value_scatter(solver);
  atomic_store_explicit(&shared->stop, 0, memory_order_release);
  shared->invalid = 0;
  shared->pivot_rejected = 0;
  shared->singular = 0;
  shared->rejected_pivot = KLS_KLU_EMPTY;
  shared->rejected_pivot_col = KLS_KLU_EMPTY;
  shared->rejected_row = KLS_KLU_EMPTY;
  shared->rejected_multiplier_abs = -1.0;
  shared->rejected_pivot_abs = -1.0;
  shared->rejected_candidate_abs = -1.0;
  shared->numerical_rank = UF_long_max;
  shared->singular_col = KLS_KLU_EMPTY;
  shared->pipeline_done = pipeline_done;
  shared->pipeline_generation = pipeline_generation;
  shared->pipeline_natural_order = 0;
  atomic_store_explicit(&shared->next_pipeline_pos,
                        (unsigned long)(use_group_pipeline
                                          ? pipeline_group_begin : 0u),
                        memory_order_release);
  shared->pipeline_pos_end =
    use_group_pipeline ? solver->row_refactor_group_count : 0u;
  shared->cluster_level_count = cluster_levels;
  kls_record_row_refactor_run(solver, check_pivots, 1);

  for (int i = 0; i < thread_count; ++i) {
    pool->workers[i].shared = shared;
    pool->workers[i].x = scratch[i];
  }

  pool->active_workers = thread_count;
  pool->generation++;
  pthread_cond_broadcast(&pool->work_cond);
  while (pool->active_workers > 0) {
    pthread_cond_wait(&pool->done_cond, &shared->lock);
  }
  shared->row_refactor_mode = 0;
  pthread_mutex_unlock(&shared->lock);

  if (shared->invalid || shared->pivot_rejected ||
      (shared->singular && common->halt_if_singular)) {
    solver->egraph_worker_scratch_dirty = 1;
  }
  if (shared->invalid) {
    common->status = TRILINOS_KLU_INVALID;
    return 0;
  }
  if (shared->pivot_rejected) {
    const int prefix_ready =
      kls_egraph_refreshed_prefix(shared, shared->rejected_pivot);
    if (prefix_ready && shared->row_refactor_defer_value_scatter &&
        !kls_scatter_row_refactor_prefix_values(solver,
                                                shared->rejected_pivot)) {
      common->status = TRILINOS_KLU_INVALID;
      return 0;
    }
    solver->fast_reject_refresh_state =
      prefix_ready ? KLS_FAST_REJECT_REFRESH_PREFIX
                   : KLS_FAST_REJECT_REFRESH_UNKNOWN;
    kls_record_fast_reject_detail(solver, shared->rejected_pivot,
                                  shared->rejected_pivot_col,
                                  shared->rejected_row,
                                  shared->rejected_multiplier_abs,
                                  shared->rejected_pivot_abs,
                                  shared->rejected_candidate_abs);
    common->status = TRILINOS_KLU_OK;
    return 0;
  }
  if (shared->singular) {
    common->status = TRILINOS_KLU_SINGULAR;
    common->numerical_rank = shared->numerical_rank;
    common->singular_col = shared->singular_col;
    if (common->halt_if_singular) {
      return 0;
    }
  } else {
    common->status = TRILINOS_KLU_OK;
  }
  if (shared->row_refactor_defer_value_scatter &&
      (!kls_scatter_row_refactor_l_values(solver) ||
       !kls_scatter_row_refactor_u_values(solver))) {
    common->status = TRILINOS_KLU_INVALID;
    return 0;
  }
  return 1;
}

static int kls_single_block_mapped_refactor(kls_solver *solver,
                                            double *numeric_values,
                                            int check_pivots) {
  if (numeric_values == NULL || solver == NULL || solver->symbolic == NULL ||
      solver->numeric == NULL || solver->symbolic->nblocks != 1u ||
      solver->numeric->Udiag == NULL || solver->numeric->Lip == NULL ||
      solver->numeric->Llen == NULL || solver->numeric->Uip == NULL ||
      solver->numeric->Ulen == NULL || solver->numeric->LUbx == NULL ||
      solver->numeric->LUbx[0] == NULL || solver->numeric->Xwork == NULL ||
      !kls_build_refactor_map(solver)) {
    return -1;
  }

  trilinos_klu_l_common *common = &solver->common;
  trilinos_klu_l_numeric *numeric = solver->numeric;
  const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  double *x = (double *)numeric->Xwork;
  double *udiag = (double *)numeric->Udiag;
  UF_long *lip = numeric->Lip;
  UF_long *llen = numeric->Llen;
  UF_long *uip = numeric->Uip;
  UF_long *ulen = numeric->Ulen;
  double *lu = (double *)numeric->LUbx[0];
  const UF_long *q = symbolic->Q;

  common->status = TRILINOS_KLU_OK;
  common->numerical_rank = KLS_KLU_EMPTY;
  common->singular_col = KLS_KLU_EMPTY;
  common->nrealloc = 0;

  for (UF_long k = 0; k < solver->n; ++k) {
    for (UF_long p = solver->refactor_col_ptr[k];
         p < solver->refactor_col_ptr[k + 1u]; ++p) {
      x[solver->refactor_row_idx[p]] = numeric_values[solver->refactor_input_pos[p]];
    }

    UF_long *ui = NULL;
    double *ux = NULL;
    UF_long ucol_len = 0;
    kls_klu_get_pointer(lu, uip, ulen, k, &ui, &ux, &ucol_len);
    for (UF_long up = 0; up < ucol_len; ++up) {
      const UF_long j = ui[up];
      const double ujk = x[j];
      x[j] = 0.0;
      ux[up] = ujk;

      UF_long *li = NULL;
      double *lx = NULL;
      UF_long lcol_len = 0;
      kls_klu_get_pointer(lu, lip, llen, j, &li, &lx, &lcol_len);
      kls_scatter_subtract(x, li, lx, lcol_len, ujk);
    }

    const double ukk = x[k];
    x[k] = 0.0;
    if (ukk == 0.0) {
      common->status = TRILINOS_KLU_SINGULAR;
      if (common->numerical_rank == KLS_KLU_EMPTY) {
        common->numerical_rank = k;
        common->singular_col = q[k];
      }
      if (common->halt_if_singular) {
        memset(x, 0, (size_t)solver->n * sizeof(*x));
        return 0;
      }
    }
    udiag[k] = ukk;

    UF_long *li = NULL;
    double *lx = NULL;
    UF_long lcol_len = 0;
    kls_klu_get_pointer(lu, lip, llen, k, &li, &lx, &lcol_len);
    UF_long rejected_row = KLS_KLU_EMPTY;
    UF_long rejected_local_row = KLS_KLU_EMPTY;
    double rejected_multiplier_abs = -1.0;
    double rejected_pivot_abs = -1.0;
    double rejected_candidate_abs = -1.0;
    if (check_pivots &&
        kls_checked_refactor_best_reject_candidate(
          li, lcol_len, x, ukk, common->tol, 0u, &rejected_row,
          &rejected_local_row, &rejected_multiplier_abs,
          &rejected_pivot_abs, &rejected_candidate_abs)) {
      memset(x, 0, (size_t)solver->n * sizeof(*x));
      solver->fast_reject_refresh_state =
        KLS_FAST_REJECT_REFRESH_PREFIX;
      kls_record_fast_reject_detail(solver, k, q[k], rejected_row,
                                    rejected_multiplier_abs,
                                    rejected_pivot_abs,
                                    rejected_candidate_abs);
      common->status = TRILINOS_KLU_OK;
      return 0;
    }
    for (UF_long p = 0; p < lcol_len; ++p) {
      const UF_long i = li[p];
      const double lij = x[i] / ukk;
      lx[p] = lij;
      x[i] = 0.0;
    }
  }

  return 1;
}

static void kls_egraph_refactor_record_invalid(
  kls_egraph_refactor_shared *shared) {
  pthread_mutex_lock(&shared->lock);
  shared->invalid = 1;
  atomic_store_explicit(&shared->stop, 1, memory_order_release);
  pthread_mutex_unlock(&shared->lock);
}

static void kls_egraph_refactor_record_reject(
  kls_egraph_refactor_shared *shared,
  UF_long rejected_pivot,
  UF_long rejected_pivot_col,
  UF_long rejected_row,
  double rejected_multiplier_abs,
  double rejected_pivot_abs,
  double rejected_candidate_abs) {
  pthread_mutex_lock(&shared->lock);
  shared->pivot_rejected = 1;
  if (rejected_pivot != KLS_KLU_EMPTY &&
      (shared->rejected_pivot == KLS_KLU_EMPTY ||
       rejected_pivot < shared->rejected_pivot)) {
    shared->rejected_pivot = rejected_pivot;
    shared->rejected_pivot_col = rejected_pivot_col;
    shared->rejected_row = rejected_row;
    shared->rejected_multiplier_abs = rejected_multiplier_abs;
    shared->rejected_pivot_abs = rejected_pivot_abs;
    shared->rejected_candidate_abs = rejected_candidate_abs;
  }
  atomic_store_explicit(&shared->stop, 1, memory_order_release);
  pthread_mutex_unlock(&shared->lock);
}

static void kls_egraph_refactor_record_singular(
  kls_egraph_refactor_shared *shared,
  UF_long numerical_rank,
  UF_long singular_col) {
  pthread_mutex_lock(&shared->lock);
  if (!shared->singular || numerical_rank < shared->numerical_rank) {
    shared->singular = 1;
    shared->numerical_rank = numerical_rank;
    shared->singular_col = singular_col;
  }
  if (shared->solver->common.halt_if_singular) {
    atomic_store_explicit(&shared->stop, 1, memory_order_release);
  }
  pthread_mutex_unlock(&shared->lock);
}

static int kls_egraph_refactor_should_stop(
  kls_egraph_refactor_shared *shared) {
  return atomic_load_explicit(&shared->stop, memory_order_acquire) != 0;
}

static int kls_egraph_refactor_value(kls_egraph_refactor_shared *shared,
                                     UF_long input_pos,
                                     double *value_out) {
  kls_solver *solver = shared->solver;
  if (input_pos >= solver->nnz) {
    return 0;
  }
  double value = shared->values[input_pos];
  if (shared->scale > 0) {
    const UF_long oldrow = solver->row_idx[input_pos];
    if (oldrow >= solver->n || shared->rs == NULL ||
        shared->rs[oldrow] == 0.0) {
      return 0;
    }
    value /= shared->rs[oldrow];
  }
  *value_out = value;
  return 1;
}

static void kls_egraph_refactor_mark_done(
  kls_egraph_refactor_shared *shared,
  UF_long col) {
  if (shared->pipeline_done != NULL) {
    atomic_store_explicit(&shared->pipeline_done[col],
                          shared->pipeline_generation,
                          memory_order_release);
  }
}

static int kls_egraph_refactor_wait_done(
  kls_egraph_refactor_shared *shared,
  UF_long col) {
  if (shared->pipeline_done == NULL) {
    return 1;
  }
  const unsigned int generation = shared->pipeline_generation;
  unsigned spin = 0;
  while (atomic_load_explicit(&shared->pipeline_done[col],
                              memory_order_acquire) != generation) {
    if ((spin++ & 1023u) == 0u &&
        kls_egraph_refactor_should_stop(shared)) {
      return 0;
    }
  }
  return 1;
}

static int kls_egraph_refreshed_prefix(
  const kls_egraph_refactor_shared *shared,
  UF_long rejected_pivot) {
  if (shared == NULL || shared->pipeline_done == NULL ||
      rejected_pivot == KLS_KLU_EMPTY ||
      shared->solver == NULL || rejected_pivot >= shared->solver->n) {
    return 0;
  }
  const unsigned int generation = shared->pipeline_generation;
  if (generation == 0u) {
    return 0;
  }
  for (UF_long col = 0; col < rejected_pivot; ++col) {
    if (atomic_load_explicit(&shared->pipeline_done[col],
                             memory_order_acquire) != generation) {
      return 0;
    }
  }
  return 1;
}

static int kls_egraph_refactor_single_unscaled_column(
  kls_egraph_refactor_worker *worker,
  UF_long k,
  int wait_for_dependencies) {
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  trilinos_klu_l_numeric *numeric = solver->numeric;
  const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  double *x = worker->x;
  double *udiag = (double *)numeric->Udiag;
  UF_long **l_indices = solver->refactor_l_indices;
  double **l_values = solver->refactor_l_values;
  UF_long **u_indices = solver->refactor_u_indices;
  double **u_values = solver->refactor_u_values;
  if (solver->refactor_lu_pointer_count != solver->n ||
      l_indices == NULL || l_values == NULL ||
      u_indices == NULL || u_values == NULL) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }

  /* The EGraph dispatcher validates the map, LU arrays, and U topological
     order once before launching workers; keep this hot kernel branch-light. */
  for (UF_long p = solver->refactor_col_ptr[k];
       p < solver->refactor_col_ptr[k + 1u]; ++p) {
    x[solver->refactor_row_idx[p]] =
      shared->values[solver->refactor_input_pos[p]];
  }

  UF_long *ui = u_indices[k];
  double *ux = u_values[k];
  UF_long ucol_len = numeric->Ulen[k];
  for (UF_long up = 0; up < ucol_len; ++up) {
    const UF_long j = ui[up];
    if (wait_for_dependencies &&
        !kls_egraph_refactor_wait_done(shared, j)) {
      return 0;
    }
    const double ujk = x[j];
    x[j] = 0.0;
    ux[up] = ujk;

    UF_long *li = l_indices[j];
    double *lx = l_values[j];
    UF_long lcol_len = numeric->Llen[j];
    kls_scatter_subtract(x, li, lx, lcol_len, ujk);
  }

  const double ukk = x[k];
  x[k] = 0.0;
  if (ukk == 0.0) {
    kls_egraph_refactor_record_singular(shared, k, symbolic->Q[k]);
    if (solver->common.halt_if_singular) {
      return 0;
    }
  }
  udiag[k] = ukk;

  UF_long *li = l_indices[k];
  double *lx = l_values[k];
  UF_long lcol_len = numeric->Llen[k];
  UF_long rejected_row = KLS_KLU_EMPTY;
  UF_long rejected_local_row = KLS_KLU_EMPTY;
  double rejected_multiplier_abs = -1.0;
  double rejected_pivot_abs = -1.0;
  double rejected_candidate_abs = -1.0;
  if (shared->check_pivots &&
      kls_checked_refactor_best_reject_candidate(
        li, lcol_len, x, ukk, solver->common.tol, 0u, &rejected_row,
        &rejected_local_row, &rejected_multiplier_abs, &rejected_pivot_abs,
        &rejected_candidate_abs)) {
    x[rejected_local_row] = 0.0;
    kls_egraph_refactor_record_reject(shared, k, symbolic->Q[k],
                                      rejected_row,
                                      rejected_multiplier_abs,
                                      rejected_pivot_abs,
                                      rejected_candidate_abs);
    return 0;
  }
  for (UF_long p = 0; p < lcol_len; ++p) {
    const UF_long i = li[p];
    const double lij = x[i] / ukk;
    lx[p] = lij;
    x[i] = 0.0;
  }
  return 1;
}

static int kls_egraph_refactor_btf_unscaled_column(
  kls_egraph_refactor_worker *worker,
  UF_long k,
  int wait_for_dependencies) {
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  trilinos_klu_l_numeric *numeric = solver->numeric;
  const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  double *x = worker->x;
  double *udiag = (double *)numeric->Udiag;
  UF_long **l_indices = solver->refactor_l_indices;
  double **l_values = solver->refactor_l_values;
  UF_long **u_indices = solver->refactor_u_indices;
  double **u_values = solver->refactor_u_values;
  if (solver->refactor_lu_pointer_count != solver->n ||
      l_indices == NULL || l_values == NULL ||
      u_indices == NULL || u_values == NULL) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }

  /* The dispatcher validates the map, block metadata, Offx arrays, and
     topological U pattern before launching workers; keep this BTF hot path
     aligned with the branch-light single-block kernel. */
  const UF_long block = solver->refactor_col_block[k];
  const UF_long k1 = symbolic->R[block];
  const UF_long k2 = symbolic->R[block + 1u];
  const UF_long nk = k2 - k1;
  const UF_long local_k = k - k1;

  UF_long poff = numeric->Offp[k];
  const UF_long poff_end = numeric->Offp[k + 1u];
  double *offx = (double *)numeric->Offx;
  for (UF_long p = solver->refactor_col_ptr[k];
       p < solver->refactor_block_start[k]; ++p) {
    if (poff >= poff_end) {
      kls_egraph_refactor_record_invalid(shared);
      return 0;
    }
    offx[poff++] = shared->values[solver->refactor_input_pos[p]];
  }

  if (nk == 1u) {
    double pivot = 0.0;
    for (UF_long p = solver->refactor_block_start[k];
         p < solver->refactor_col_ptr[k + 1u]; ++p) {
      pivot = shared->values[solver->refactor_input_pos[p]];
    }
    udiag[k] = pivot;
    if (pivot == 0.0) {
      kls_egraph_refactor_record_singular(shared, k, symbolic->Q[k]);
      if (solver->common.halt_if_singular) {
        return 0;
      }
    }
    return 1;
  }

  for (UF_long p = solver->refactor_block_start[k];
       p < solver->refactor_col_ptr[k + 1u]; ++p) {
    const UF_long global_row = solver->refactor_row_idx[p];
    x[global_row - k1] = shared->values[solver->refactor_input_pos[p]];
  }

  UF_long *ui = u_indices[k];
  double *ux = u_values[k];
  UF_long ucol_len = numeric->Ulen[k];
  for (UF_long up = 0; up < ucol_len; ++up) {
    const UF_long j = ui[up];
    if (wait_for_dependencies &&
        !kls_egraph_refactor_wait_done(shared, k1 + j)) {
      return 0;
    }
    const double ujk = x[j];
    x[j] = 0.0;
    ux[up] = ujk;

    UF_long *li = l_indices[k1 + j];
    double *lx = l_values[k1 + j];
    UF_long lcol_len = numeric->Llen[k1 + j];
    kls_scatter_subtract(x, li, lx, lcol_len, ujk);
  }

  const double ukk = x[local_k];
  x[local_k] = 0.0;
  if (ukk == 0.0) {
    kls_egraph_refactor_record_singular(shared, k, symbolic->Q[k]);
    if (solver->common.halt_if_singular) {
      return 0;
    }
  }
  udiag[k] = ukk;

  UF_long *li = l_indices[k];
  double *lx = l_values[k];
  UF_long lcol_len = numeric->Llen[k];
  UF_long rejected_row = KLS_KLU_EMPTY;
  UF_long rejected_local_row = KLS_KLU_EMPTY;
  double rejected_multiplier_abs = -1.0;
  double rejected_pivot_abs = -1.0;
  double rejected_candidate_abs = -1.0;
  if (shared->check_pivots &&
      kls_checked_refactor_best_reject_candidate(
        li, lcol_len, x, ukk, solver->common.tol, k1, &rejected_row,
        &rejected_local_row, &rejected_multiplier_abs, &rejected_pivot_abs,
        &rejected_candidate_abs)) {
    x[rejected_local_row] = 0.0;
    kls_egraph_refactor_record_reject(shared, k, symbolic->Q[k],
                                      rejected_row,
                                      rejected_multiplier_abs,
                                      rejected_pivot_abs,
                                      rejected_candidate_abs);
    return 0;
  }
  for (UF_long p = 0; p < lcol_len; ++p) {
    const UF_long i = li[p];
    const double lij = x[i] / ukk;
    lx[p] = lij;
    x[i] = 0.0;
  }
  return 1;
}

static int kls_egraph_refactor_column(kls_egraph_refactor_worker *worker,
                                      UF_long k,
                                      int wait_for_dependencies) {
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  trilinos_klu_l_numeric *numeric = solver->numeric;
  const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  double *x = worker->x;
  double *udiag = (double *)numeric->Udiag;
  UF_long **l_indices = solver->refactor_l_indices;
  double **l_values = solver->refactor_l_values;
  UF_long **u_indices = solver->refactor_u_indices;
  double **u_values = solver->refactor_u_values;
  if (symbolic != NULL && shared->scale <= 0) {
    if (symbolic->nblocks == 1u) {
      return kls_egraph_refactor_single_unscaled_column(
        worker, k, wait_for_dependencies);
    }
    if (symbolic->nblocks > 1u && symbolic->maxblock >= 30000u) {
      return kls_egraph_refactor_btf_unscaled_column(
        worker, k, wait_for_dependencies);
    }
  }
  if (solver->refactor_lu_pointer_count != solver->n ||
      l_indices == NULL || l_values == NULL ||
      u_indices == NULL || u_values == NULL) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }
  if (k >= solver->n || symbolic->R == NULL || symbolic->Q == NULL) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }
  UF_long block = KLS_KLU_EMPTY;
  if (symbolic->nblocks == 1u) {
    block = 0u;
  } else if (solver->refactor_col_block != NULL) {
    block = solver->refactor_col_block[k];
  } else {
    block = kls_block_for_pivot(solver, k);
  }
  if (block == KLS_KLU_EMPTY || block >= symbolic->nblocks) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }

  const UF_long k1 = symbolic->R[block];
  const UF_long k2 = symbolic->R[block + 1u];
  const UF_long nk = k2 - k1;
  const UF_long local_k = k - k1;
  const int btf_block = symbolic->nblocks != 1u;

  if (local_k >= nk ||
      solver->refactor_col_ptr == NULL ||
      solver->refactor_row_idx == NULL ||
      solver->refactor_input_pos == NULL ||
      (btf_block && (solver->refactor_block_start == NULL ||
                     solver->refactor_col_block == NULL))) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }

  if (nk == 1u) {
    if (!btf_block || numeric->Offp == NULL || numeric->Offx == NULL) {
      kls_egraph_refactor_record_invalid(shared);
      return 0;
    }
    UF_long poff = numeric->Offp[k];
    const UF_long poff_end = numeric->Offp[k + 1u];
    double *offx = (double *)numeric->Offx;
    double pivot = 0.0;
    for (UF_long p = solver->refactor_col_ptr[k];
         p < solver->refactor_block_start[k]; ++p) {
      if (poff >= poff_end) {
        kls_egraph_refactor_record_invalid(shared);
        return 0;
      }
      double value = 0.0;
      if (!kls_egraph_refactor_value(shared, solver->refactor_input_pos[p],
                                     &value)) {
        kls_egraph_refactor_record_invalid(shared);
        return 0;
      }
      offx[poff++] = value;
    }
    for (UF_long p = solver->refactor_block_start[k];
         p < solver->refactor_col_ptr[k + 1u]; ++p) {
      if (solver->refactor_row_idx[p] != k) {
        kls_egraph_refactor_record_invalid(shared);
        return 0;
      }
      if (!kls_egraph_refactor_value(shared, solver->refactor_input_pos[p],
                                     &pivot)) {
        kls_egraph_refactor_record_invalid(shared);
        return 0;
      }
    }
    udiag[k] = pivot;
    if (pivot == 0.0) {
      kls_egraph_refactor_record_singular(shared, k, symbolic->Q[k]);
      if (solver->common.halt_if_singular) {
        return 0;
      }
    }
    return 1;
  }

  UF_long *llen = numeric->Llen + k1;
  UF_long *ulen = numeric->Ulen + k1;
  double *lu = (double *)numeric->LUbx[block];
  if (lu == NULL) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }

  if (btf_block) {
    if (numeric->Offp == NULL || numeric->Offx == NULL) {
      kls_egraph_refactor_record_invalid(shared);
      return 0;
    }
    UF_long poff = numeric->Offp[k];
    const UF_long poff_end = numeric->Offp[k + 1u];
    double *offx = (double *)numeric->Offx;
    for (UF_long p = solver->refactor_col_ptr[k];
         p < solver->refactor_block_start[k]; ++p) {
      if (poff >= poff_end) {
        kls_egraph_refactor_record_invalid(shared);
        return 0;
      }
      double value = 0.0;
      if (!kls_egraph_refactor_value(shared, solver->refactor_input_pos[p],
                                     &value)) {
        kls_egraph_refactor_record_invalid(shared);
        return 0;
      }
      offx[poff++] = value;
    }
    for (UF_long p = solver->refactor_block_start[k];
         p < solver->refactor_col_ptr[k + 1u]; ++p) {
      const UF_long global_row = solver->refactor_row_idx[p];
      if (global_row < k1 || global_row >= k2) {
        kls_egraph_refactor_record_invalid(shared);
        return 0;
      }
      double value = 0.0;
      if (!kls_egraph_refactor_value(shared, solver->refactor_input_pos[p],
                                     &value)) {
        kls_egraph_refactor_record_invalid(shared);
        return 0;
      }
      x[global_row - k1] = value;
    }
  } else {
    for (UF_long p = solver->refactor_col_ptr[k];
         p < solver->refactor_col_ptr[k + 1u]; ++p) {
      const UF_long row = solver->refactor_row_idx[p];
      const UF_long input_pos = solver->refactor_input_pos[p];
      if (row >= solver->n || input_pos >= solver->nnz) {
        kls_egraph_refactor_record_invalid(shared);
        return 0;
      }
      double value = shared->values[input_pos];
      if (shared->scale > 0) {
        const UF_long oldrow = solver->row_idx[input_pos];
        if (oldrow >= solver->n || shared->rs == NULL ||
            shared->rs[oldrow] == 0.0) {
          kls_egraph_refactor_record_invalid(shared);
          return 0;
        }
        value /= shared->rs[oldrow];
      }
      x[row] = value;
    }
  }

  UF_long *ui = u_indices[k];
  double *ux = u_values[k];
  UF_long ucol_len = ulen[local_k];
  for (UF_long up = 0; up < ucol_len; ++up) {
    const UF_long j = ui[up];
    if (j >= local_k) {
      kls_egraph_refactor_record_invalid(shared);
      return 0;
    }
    if (wait_for_dependencies &&
        !kls_egraph_refactor_wait_done(shared, k1 + j)) {
      return 0;
    }
    const double ujk = x[j];
    x[j] = 0.0;
    ux[up] = ujk;

    UF_long *li = l_indices[k1 + j];
    double *lx = l_values[k1 + j];
    UF_long lcol_len = llen[j];
    kls_scatter_subtract(x, li, lx, lcol_len, ujk);
  }

  const double ukk = x[local_k];
  x[local_k] = 0.0;
  if (ukk == 0.0) {
    kls_egraph_refactor_record_singular(shared, k, symbolic->Q[k]);
    if (solver->common.halt_if_singular) {
      return 0;
    }
  }
  udiag[k] = ukk;

  UF_long *li = l_indices[k];
  double *lx = l_values[k];
  UF_long lcol_len = llen[local_k];
  UF_long rejected_row = KLS_KLU_EMPTY;
  UF_long rejected_local_row = KLS_KLU_EMPTY;
  double rejected_multiplier_abs = -1.0;
  double rejected_pivot_abs = -1.0;
  double rejected_candidate_abs = -1.0;
  if (shared->check_pivots &&
      kls_checked_refactor_best_reject_candidate(
        li, lcol_len, x, ukk, solver->common.tol, k1, &rejected_row,
        &rejected_local_row, &rejected_multiplier_abs, &rejected_pivot_abs,
        &rejected_candidate_abs)) {
    x[rejected_local_row] = 0.0;
    kls_egraph_refactor_record_reject(shared, k, symbolic->Q[k],
                                      rejected_row,
                                      rejected_multiplier_abs,
                                      rejected_pivot_abs,
                                      rejected_candidate_abs);
    return 0;
  }
  for (UF_long p = 0; p < lcol_len; ++p) {
    const UF_long i = li[p];
    const double lij = x[i] / ukk;
    lx[p] = lij;
    x[i] = 0.0;
  }
  return 1;
}

static int kls_egraph_refactor_dispatch_column(
  kls_egraph_refactor_worker *worker,
  UF_long k,
  int wait_for_dependencies) {
  if (worker == NULL || worker->shared == NULL) {
    return 0;
  }
  switch (worker->shared->kernel) {
    case KLS_EGRAPH_REFACTOR_KERNEL_SINGLE_UNSCALED:
      return kls_egraph_refactor_single_unscaled_column(
        worker, k, wait_for_dependencies);
    case KLS_EGRAPH_REFACTOR_KERNEL_BTF_UNSCALED:
      return kls_egraph_refactor_btf_unscaled_column(
        worker, k, wait_for_dependencies);
    case KLS_EGRAPH_REFACTOR_KERNEL_GENERIC:
    default:
      return kls_egraph_refactor_column(worker, k, wait_for_dependencies);
  }
}

static void kls_egraph_refactor_worker_run(kls_egraph_refactor_worker *worker) {
  if (worker == NULL || worker->shared == NULL) {
    return;
  }
  kls_egraph_refactor_shared *shared = worker->shared;
  const kls_solver *solver = shared->solver;
  if (solver == NULL) {
    return;
  }

  UF_long cluster_levels = shared->cluster_level_count;
  if (cluster_levels > solver->refactor_level_count) {
    cluster_levels = solver->refactor_level_count;
  }

  for (UF_long level = 0; level < cluster_levels; ++level) {
    const int stopped = kls_egraph_refactor_should_stop(shared);

    if (!stopped) {
      const UF_long begin = solver->refactor_level_ptr[level];
      const UF_long end = solver->refactor_level_ptr[level + 1u];
      if (solver->refactor_level_thread_ptr != NULL &&
          solver->refactor_level_thread_count == shared->thread_count) {
        const UF_long *parts =
          solver->refactor_level_thread_ptr +
          level * (UF_long)(shared->thread_count + 1);
        for (UF_long pos = parts[worker->tid];
             pos < parts[worker->tid + 1]; ++pos) {
          const UF_long col = solver->refactor_level_cols[pos];
          if (!kls_egraph_refactor_dispatch_column(worker, col, 0)) {
            break;
          }
          kls_egraph_refactor_mark_done(shared, col);
        }
      } else {
        for (UF_long pos = begin + (UF_long)worker->tid;
             pos < end; pos += (UF_long)shared->thread_count) {
          const UF_long col = solver->refactor_level_cols[pos];
          if (!kls_egraph_refactor_dispatch_column(worker, col, 0)) {
            break;
          }
          kls_egraph_refactor_mark_done(shared, col);
        }
      }
    }

    (void)pthread_barrier_wait(&shared->barrier);

    if (kls_egraph_refactor_should_stop(shared)) {
      break;
    }
  }

  if (shared->pipeline_done != NULL &&
      cluster_levels < solver->refactor_level_count &&
      !kls_egraph_refactor_should_stop(shared)) {
    if (shared->pipeline_natural_order) {
      for (;;) {
        if (kls_egraph_refactor_should_stop(shared)) {
          break;
        }
        const UF_long col =
          (UF_long)atomic_fetch_add_explicit(&shared->next_pipeline_pos, 1ul,
                                             memory_order_relaxed);
        if (col >= shared->pipeline_pos_end) {
          break;
        }
        if (!kls_egraph_refactor_dispatch_column(worker, col, 1)) {
          break;
        }
        kls_egraph_refactor_mark_done(shared, col);
      }
    } else if (solver->refactor_level_cols != NULL) {
      for (;;) {
        if (kls_egraph_refactor_should_stop(shared)) {
          break;
        }
        const UF_long pos =
          (UF_long)atomic_fetch_add_explicit(&shared->next_pipeline_pos, 1ul,
                                             memory_order_relaxed);
        if (pos >= shared->pipeline_pos_end) {
          break;
        }
        const UF_long col = solver->refactor_level_cols[pos];
        if (!kls_egraph_refactor_dispatch_column(worker, col, 1)) {
          break;
        }
        kls_egraph_refactor_mark_done(shared, col);
      }
    } else {
      kls_egraph_refactor_record_invalid(shared);
    }
  }
}

static void *kls_egraph_refactor_pool_worker_main(void *arg) {
  kls_egraph_refactor_worker *worker = (kls_egraph_refactor_worker *)arg;
  kls_egraph_refactor_pool *pool = worker->pool;
  kls_egraph_refactor_shared *shared = worker->shared;
  unsigned long seen_generation = 0;

  pthread_mutex_lock(&shared->lock);
  for (;;) {
    while (!pool->shutdown && pool->generation == seen_generation) {
      pthread_cond_wait(&pool->work_cond, &shared->lock);
    }
    if (pool->shutdown) {
      pthread_mutex_unlock(&shared->lock);
      return NULL;
    }
    seen_generation = pool->generation;
    pthread_mutex_unlock(&shared->lock);

    if (shared->row_refactor_mode) {
      kls_row_refactor_worker_run(worker);
    } else {
      kls_egraph_refactor_worker_run(worker);
    }

    pthread_mutex_lock(&shared->lock);
    pool->active_workers--;
    if (pool->active_workers == 0) {
      pthread_cond_signal(&pool->done_cond);
    }
  }
}

static void destroy_egraph_refactor_pool(kls_solver *solver) {
  if (solver == NULL || solver->egraph_pool == NULL) {
    return;
  }
  kls_egraph_refactor_pool *pool = solver->egraph_pool;
  if (pool->lock_initialized) {
    pthread_mutex_lock(&pool->shared.lock);
    pool->shutdown = 1;
    if (pool->conds_initialized) {
      pthread_cond_broadcast(&pool->work_cond);
    }
    pthread_mutex_unlock(&pool->shared.lock);
  }
  if (pool->threads != NULL) {
    for (int i = 0; i < pool->created_count; ++i) {
      pthread_join(pool->threads[i], NULL);
    }
  }
  if (pool->barrier_initialized) {
    pthread_barrier_destroy(&pool->shared.barrier);
  }
  if (pool->conds_initialized) {
    pthread_cond_destroy(&pool->work_cond);
    pthread_cond_destroy(&pool->done_cond);
  }
  if (pool->lock_initialized) {
    pthread_mutex_destroy(&pool->shared.lock);
  }
  if (pool->workers != NULL) {
    for (int i = 0; i < pool->thread_count; ++i) {
      free(pool->workers[i].segment_panel);
    }
  }
  free(pool->threads);
  free(pool->workers);
  free(pool);
  solver->egraph_pool = NULL;
}

static kls_egraph_refactor_pool *ensure_egraph_refactor_pool(
  kls_solver *solver,
  int thread_count) {
  if (solver == NULL || thread_count < 2) {
    return NULL;
  }
  if (solver->egraph_pool != NULL &&
      solver->egraph_pool->thread_count != thread_count) {
    destroy_egraph_refactor_pool(solver);
  }
  if (solver->egraph_pool != NULL) {
    return solver->egraph_pool;
  }

  kls_egraph_refactor_pool *pool =
    (kls_egraph_refactor_pool *)calloc(1, sizeof(*pool));
  if (pool == NULL) {
    return NULL;
  }
  pool->thread_count = thread_count;
  pool->threads =
    (pthread_t *)calloc((size_t)thread_count, sizeof(*pool->threads));
  pool->workers =
    (kls_egraph_refactor_worker *)calloc((size_t)thread_count,
                                         sizeof(*pool->workers));
  if (pool->threads == NULL || pool->workers == NULL) {
    free(pool->threads);
    free(pool->workers);
    free(pool);
    return NULL;
  }
  if (pthread_mutex_init(&pool->shared.lock, NULL) != 0) {
    free(pool->threads);
    free(pool->workers);
    free(pool);
    return NULL;
  }
  pool->lock_initialized = 1;
  if (pthread_cond_init(&pool->work_cond, NULL) != 0) {
    pthread_mutex_destroy(&pool->shared.lock);
    free(pool->threads);
    free(pool->workers);
    free(pool);
    return NULL;
  }
  if (pthread_cond_init(&pool->done_cond, NULL) != 0) {
    pthread_cond_destroy(&pool->work_cond);
    pthread_mutex_destroy(&pool->shared.lock);
    free(pool->threads);
    free(pool->workers);
    free(pool);
    return NULL;
  }
  pool->conds_initialized = 1;
  if (pthread_barrier_init(&pool->shared.barrier, NULL,
                           (unsigned)thread_count) != 0) {
    pthread_cond_destroy(&pool->done_cond);
    pthread_cond_destroy(&pool->work_cond);
    pthread_mutex_destroy(&pool->shared.lock);
    free(pool->threads);
    free(pool->workers);
    free(pool);
    return NULL;
  }
  pool->barrier_initialized = 1;
  atomic_init(&pool->shared.stop, 0);
  atomic_init(&pool->shared.next_pipeline_pos, 0ul);

  for (int i = 0; i < thread_count; ++i) {
    pool->workers[i].shared = &pool->shared;
    pool->workers[i].pool = pool;
    pool->workers[i].tid = i;
    if (pthread_create(&pool->threads[i], NULL,
                       kls_egraph_refactor_pool_worker_main,
                       &pool->workers[i]) != 0) {
      break;
    }
    pool->created_count++;
  }
  if (pool->created_count != thread_count) {
    solver->egraph_pool = pool;
    destroy_egraph_refactor_pool(solver);
    return NULL;
  }
  solver->egraph_pool = pool;
  return pool;
}

static int kls_egraph_medium_heavy_dominant_btf_shape(
  const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL ||
      solver->symbolic->nblocks <= 1u || solver->n == 0u) {
    return 0;
  }
  const double coverage =
    (double)solver->symbolic->maxblock / (double)solver->n;
  return coverage >= 0.85 && coverage < 0.95 &&
         solver->symbolic->maxblock >= 30000u &&
         solver->symbolic->nblocks <= 5000u &&
         solver->common.flops >= 1.5e8;
}

static int kls_egraph_all_pipeline_dominant_btf_shape(
  const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL ||
      solver->symbolic->nblocks <= 1u || solver->n == 0u) {
    return 0;
  }
  const double coverage =
    (double)solver->symbolic->maxblock / (double)solver->n;
  /* Thousands of tiny fringe blocks around one 95%+ block do not amortize
     BTF block scheduling, but the exact EGraph is wide enough to run without
     cluster barriers by waiting only on actual U-pattern predecessors. */
  return coverage >= 0.95 &&
         solver->symbolic->nblocks >= 1024u &&
         solver->symbolic->maxblock < 90000u &&
         solver->common.flops >= 1.0e8;
}

static int kls_egraph_all_pipeline_huge_single_shape(
  const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL ||
      solver->symbolic->nblocks != 1u || solver->numeric == NULL ||
      solver->n < 100000u) {
    return 0;
  }
  if (solver->common.scale > 0) {
    if (solver->numeric->Rs == NULL || solver->numeric->Pnum == NULL ||
        solver->numeric->Xwork == NULL) {
      return 0;
    }
  } else if (solver->numeric->Rs != NULL) {
    return 0;
  }
  return solver->common.flops >= 1.0e9;
}

static int kls_egraph_non_dominant_many_block_shape(
  const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      solver->symbolic->nblocks < 100000u || solver->n == 0u ||
      solver->common.scale > 0 || solver->numeric->Rs != NULL) {
    return 0;
  }
  const double coverage =
    (double)solver->symbolic->maxblock / (double)solver->n;
  /* Extremely fragmented BTF can hide one expensive diagonal block behind
     hundreds of thousands of singleton blocks.  The BTF worker pool has too
     little off-block work to help, but the retained exact EGraph exposes
     enough intra-block independence to run the refactor in parallel. */
  return coverage < 0.50 &&
         solver->symbolic->maxblock >= 90000u &&
         solver->common.flops >= 5.0e8;
}

static int kls_egraph_compact_dominant_btf_shape(const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL ||
      solver->symbolic->nblocks < 8u || solver->n == 0u ||
      solver->common.scale > 0) {
    return 0;
  }
  const double coverage =
    (double)solver->symbolic->maxblock / (double)solver->n;
  return coverage >= 0.95 &&
         solver->symbolic->nblocks <= 512u &&
         solver->symbolic->maxblock >= 10000u &&
         solver->symbolic->maxblock < 30000u &&
         solver->common.flops >= 2.0e7;
}

static int kls_egraph_small_compact_dominant_btf_shape(
  const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL ||
      solver->symbolic->nblocks < 8u ||
      solver->symbolic->nblocks > 512u ||
      solver->n == 0u || solver->common.scale > 0) {
    return 0;
  }
  const double coverage =
    (double)solver->symbolic->maxblock / (double)solver->n;
  return coverage >= 0.90 &&
         solver->symbolic->maxblock >= 3000u &&
         solver->symbolic->maxblock < 10000u &&
         solver->common.flops >= 2.0e6;
}

static int kls_egraph_low_work_dominant_btf_shape(
  const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      solver->symbolic->nblocks <= 1u ||
      solver->symbolic->nblocks > 5000u ||
      solver->n < 30000u || solver->n > 120000u ||
      solver->common.scale > 0 || solver->common.noffdiag != 0u) {
    return 0;
  }
  const double coverage =
    (double)solver->symbolic->maxblock / (double)solver->n;
  return coverage >= 0.95 &&
         solver->symbolic->maxblock >= 60000u &&
         solver->symbolic->maxblock < 90000u &&
         solver->common.flops >= 1.0e7 &&
         solver->common.flops < 3.0e7;
}

static int kls_egraph_scaled_medium_dominant_btf_shape(
  const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      solver->symbolic->nblocks < 8u || solver->symbolic->nblocks > 512u ||
      solver->common.scale <= 0 || solver->n == 0u) {
    return 0;
  }
  const double coverage =
    (double)solver->symbolic->maxblock / (double)solver->n;
  return coverage >= 0.95 &&
         solver->symbolic->maxblock >= 30000u &&
         solver->symbolic->maxblock < 60000u &&
         solver->common.flops >= 3.0e7 &&
         solver->numeric->lnz + solver->numeric->unz >= 1000000u;
}

static int kls_egraph_moderate_single_block_shape(const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      solver->symbolic->nblocks != 1u || solver->common.scale > 0 ||
      solver->n < 30000u || solver->n >= 100000u) {
    return 0;
  }
  return solver->common.flops >= 5.0e7 &&
         solver->numeric->lnz + solver->numeric->unz >= 1000000u;
}

static int kls_egraph_low_work_single_block_shape(const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      solver->symbolic->nblocks != 1u || solver->common.scale > 0 ||
      solver->n < 15000u || solver->n > 250000u ||
      solver->common.noffdiag > 8u) {
    return 0;
  }
  return solver->common.flops >= 4.0e6 &&
         solver->common.flops < 5.0e7 &&
         solver->numeric->lnz + solver->numeric->unz >= 100000u;
}

static int kls_egraph_dominant_btf_shape(const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL ||
      solver->symbolic->nblocks <= 1u || solver->n == 0u) {
    return 0;
  }
  const double coverage =
    (double)solver->symbolic->maxblock / (double)solver->n;
  if (coverage >= 0.95 &&
      (solver->symbolic->maxblock >= 90000u ||
       solver->common.flops >= 5.0e9)) {
    return 1;
  }
  if (kls_egraph_medium_heavy_dominant_btf_shape(solver)) {
    return 1;
  }
  if (kls_egraph_compact_dominant_btf_shape(solver)) {
    return 1;
  }
  if (kls_egraph_small_compact_dominant_btf_shape(solver)) {
    return 1;
  }
  if (kls_egraph_low_work_dominant_btf_shape(solver)) {
    return 1;
  }
  if (kls_egraph_scaled_medium_dominant_btf_shape(solver)) {
    return 1;
  }
  return coverage >= 0.85 &&
         solver->symbolic->maxblock >= 100000u &&
         solver->symbolic->nblocks <= 20000u &&
         solver->common.flops >= 2.0e9;
}

static int kls_egraph_compact_large_dominant_btf_shape(
  const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL ||
      solver->symbolic->nblocks < 128u ||
      solver->symbolic->nblocks > 512u ||
      solver->common.scale > 0 || solver->n == 0u) {
    return 0;
  }
  const double coverage =
    (double)solver->symbolic->maxblock / (double)solver->n;
  return coverage >= 0.95 &&
         solver->symbolic->maxblock >= 90000u &&
         solver->common.flops >= 5.0e8;
}

static UF_long kls_egraph_refactor_size_floor(const kls_solver *solver) {
  if (kls_egraph_all_pipeline_huge_single_shape(solver)) {
    return 100000u;
  }
  if (kls_egraph_non_dominant_many_block_shape(solver)) {
    return 90000u;
  }
  if (kls_egraph_all_pipeline_dominant_btf_shape(solver)) {
    return 30000u;
  }
  if (kls_egraph_medium_heavy_dominant_btf_shape(solver)) {
    return 30000u;
  }
  if (kls_egraph_compact_dominant_btf_shape(solver)) {
    return 10000u;
  }
  if (kls_egraph_small_compact_dominant_btf_shape(solver)) {
    return 3000u;
  }
  if (kls_egraph_low_work_dominant_btf_shape(solver)) {
    return 30000u;
  }
  if (kls_egraph_scaled_medium_dominant_btf_shape(solver)) {
    return 30000u;
  }
  if (kls_egraph_moderate_single_block_shape(solver)) {
    return 30000u;
  }
  if (kls_egraph_low_work_single_block_shape(solver)) {
    return 15000u;
  }
  return kls_egraph_dominant_btf_shape(solver) ? 50000u : 100000u;
}

static int kls_egraph_refactor_is_eligible(const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      solver->options.threads <= 1 ||
      solver->refactor_level_ptr == NULL ||
      solver->refactor_level_count == 0u ||
      solver->refactor_level_max_width < (UF_long)(4 * solver->options.threads)) {
    return 0;
  }
  const int single_block = solver->symbolic->nblocks == 1u;
  const int dominant_btf = kls_egraph_dominant_btf_shape(solver);
  const int all_pipeline_btf =
    kls_egraph_all_pipeline_dominant_btf_shape(solver);
  const int all_pipeline_single =
    kls_egraph_all_pipeline_huge_single_shape(solver);
  if (solver->refactor_level_cols == NULL && !all_pipeline_single) {
    return 0;
  }
  const int non_dominant_many_block =
    kls_egraph_non_dominant_many_block_shape(solver);
  const int medium_heavy_btf =
    kls_egraph_medium_heavy_dominant_btf_shape(solver);
  const int scaled_medium_btf =
    kls_egraph_scaled_medium_dominant_btf_shape(solver);
  const int small_compact_btf =
    kls_egraph_small_compact_dominant_btf_shape(solver);
  const int low_work_btf =
    kls_egraph_low_work_dominant_btf_shape(solver);
  const int moderate_single =
    kls_egraph_moderate_single_block_shape(solver);
  const int low_work_single =
    kls_egraph_low_work_single_block_shape(solver);
  if (solver->n < kls_egraph_refactor_size_floor(solver)) {
    return 0;
  }
  if (!single_block && !dominant_btf && !all_pipeline_btf &&
      !non_dominant_many_block) {
    return 0;
  }
  const int low_work_dominant_btf =
    dominant_btf && !medium_heavy_btf && !all_pipeline_btf &&
    solver->common.flops < 1.0e8;
  const double min_dependency_work =
    all_pipeline_single ? 1.0e9 :
    non_dominant_many_block ? 2.0e8 :
    (dominant_btf || all_pipeline_btf)
      ? (medium_heavy_btf || all_pipeline_btf ? 8.0e7 :
         scaled_medium_btf ? 1.5e7 :
         small_compact_btf ? 5.0e5 :
         low_work_btf ? 1.0e6 :
         low_work_dominant_btf ? 1.0e7 : 1.0e8)
      : (moderate_single ? 2.0e7 :
         low_work_single ? 1.0e6 : 1.5e8);
  if (solver->refactor_dependency_work < min_dependency_work) {
    return 0;
  }
  if (!single_block &&
      (solver->refactor_block_start == NULL ||
       solver->refactor_col_block == NULL ||
       solver->numeric->Offp == NULL || solver->numeric->Offx == NULL)) {
    return 0;
  }
  if (solver->common.scale > 0) {
    if (solver->numeric->Rs == NULL || solver->numeric->Pnum == NULL ||
        solver->numeric->Xwork == NULL) {
      return 0;
    }
  } else if (solver->numeric->Rs != NULL) {
    return 0;
  }
  return solver->numeric->Udiag != NULL && solver->numeric->Lip != NULL &&
         solver->numeric->Llen != NULL && solver->numeric->Uip != NULL &&
         solver->numeric->Ulen != NULL && solver->numeric->LUbx != NULL &&
         (single_block ? solver->numeric->LUbx[0] != NULL : 1) &&
         kls_refactor_map_is_eligible(solver);
}

static kls_egraph_refactor_kernel kls_egraph_refactor_kernel_for(
  const kls_solver *solver,
  int scale) {
  if (solver == NULL || solver->symbolic == NULL || scale > 0) {
    return KLS_EGRAPH_REFACTOR_KERNEL_GENERIC;
  }
  if (solver->symbolic->nblocks == 1u) {
    return KLS_EGRAPH_REFACTOR_KERNEL_SINGLE_UNSCALED;
  }
  if (solver->symbolic->nblocks > 1u &&
      solver->symbolic->maxblock >= 30000u) {
    return KLS_EGRAPH_REFACTOR_KERNEL_BTF_UNSCALED;
  }
  return KLS_EGRAPH_REFACTOR_KERNEL_GENERIC;
}

static int kls_egraph_mapped_refactor(kls_solver *solver,
                                      double *numeric_values,
                                      int check_pivots) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      numeric_values == NULL || solver->options.threads <= 1 ||
      solver->n < kls_egraph_refactor_size_floor(solver)) {
    return -1;
  }
  if (solver->common.scale <= 0 && solver->numeric->Rs != NULL) {
    return -1;
  }
  if (!kls_build_refactor_schedule(solver) ||
      !kls_build_refactor_map(solver) ||
      !kls_build_refactor_lu_pointer_cache(solver) ||
      !kls_egraph_refactor_is_eligible(solver)) {
    return -1;
  }
  trilinos_klu_l_common *common = &solver->common;
  if (common->scale > 0 &&
      !trilinos_klu_l_scale((UF_long)common->scale, solver->n,
                            solver->col_ptr, solver->row_idx,
                            numeric_values, solver->numeric->Rs, NULL,
                            common)) {
    return 0;
  }

  int thread_count = solver->options.threads;
  if ((UF_long)thread_count > solver->refactor_level_max_width) {
    thread_count = (int)solver->refactor_level_max_width;
  }
  if (thread_count < 2) {
    return -1;
  }
  const int all_pipeline =
    kls_egraph_all_pipeline_dominant_btf_shape(solver) ||
    kls_egraph_all_pipeline_huge_single_shape(solver);
  const int natural_pipeline =
    kls_egraph_all_pipeline_huge_single_shape(solver);
  const UF_long cluster_level_count =
    all_pipeline ? 0u : solver->refactor_cluster_level_count;
  atomic_uint *pipeline_done = NULL;
  unsigned int pipeline_generation = 0;
  if (all_pipeline ||
      (solver->refactor_cluster_level_count < solver->refactor_level_count &&
       solver->refactor_pipeline_column_count >= (UF_long)(2 * thread_count) &&
       solver->refactor_pipeline_work >=
         0.10 * solver->refactor_dependency_work)) {
    pipeline_done =
      ensure_egraph_pipeline_done(solver, &pipeline_generation);
  }
  if (natural_pipeline && pipeline_done == NULL) {
    return -1;
  }

  const int single_block = solver->symbolic->nblocks == 1u;
  const UF_long scratch_size =
    single_block ? solver->n : solver->symbolic->maxblock;
  double **scratch =
    ensure_egraph_worker_scratch(solver, thread_count, scratch_size);
  if (scratch == NULL) {
    return -1;
  }

  kls_egraph_refactor_pool *pool =
    ensure_egraph_refactor_pool(solver, thread_count);
  if (pool == NULL) {
    return -1;
  }
  kls_egraph_refactor_shared *shared = &pool->shared;

  pthread_mutex_lock(&shared->lock);
  if (pool->active_workers != 0) {
    pthread_mutex_unlock(&shared->lock);
    return -1;
  }

  shared->solver = solver;
  shared->values = numeric_values;
  shared->rs = solver->numeric->Rs;
  shared->check_pivots = check_pivots;
  shared->scale = (int)common->scale;
  shared->kernel =
    kls_egraph_refactor_kernel_for(solver, shared->scale);
  shared->thread_count = thread_count;
  shared->row_refactor_mode = 0;
  shared->row_refactor_defer_value_scatter = 0;
  atomic_store_explicit(&shared->stop, 0, memory_order_release);
  shared->invalid = 0;
  shared->pivot_rejected = 0;
  shared->singular = 0;
  shared->rejected_pivot = KLS_KLU_EMPTY;
  shared->rejected_pivot_col = KLS_KLU_EMPTY;
  shared->rejected_row = KLS_KLU_EMPTY;
  shared->rejected_multiplier_abs = -1.0;
  shared->rejected_pivot_abs = -1.0;
  shared->rejected_candidate_abs = -1.0;
  shared->numerical_rank = UF_long_max;
  shared->singular_col = KLS_KLU_EMPTY;
  shared->pipeline_done = pipeline_done;
  shared->pipeline_generation = pipeline_generation;
  shared->pipeline_natural_order = natural_pipeline ? 1 : 0;
  atomic_store_explicit(
    &shared->next_pipeline_pos,
    (unsigned long)(pipeline_done != NULL
                      ? (all_pipeline
                           ? 0u
                           : solver->refactor_level_ptr[cluster_level_count])
                      : solver->n),
    memory_order_release);
  shared->pipeline_pos_end = solver->n;
  shared->cluster_level_count =
    pipeline_done != NULL ? cluster_level_count : solver->refactor_level_count;

  for (int i = 0; i < thread_count; ++i) {
    pool->workers[i].shared = shared;
    pool->workers[i].x = scratch[i];
  }

  pool->active_workers = thread_count;
  pool->generation++;
  pthread_cond_broadcast(&pool->work_cond);
  while (pool->active_workers > 0) {
    pthread_cond_wait(&pool->done_cond, &shared->lock);
  }
  pthread_mutex_unlock(&shared->lock);

  if (shared->invalid || shared->pivot_rejected ||
      (shared->singular && common->halt_if_singular)) {
    solver->egraph_worker_scratch_dirty = 1;
  }

  if (shared->invalid) {
    common->status = TRILINOS_KLU_INVALID;
    return 0;
  }
  if (shared->pivot_rejected) {
    solver->fast_reject_refresh_state =
      kls_egraph_refreshed_prefix(shared, shared->rejected_pivot)
        ? KLS_FAST_REJECT_REFRESH_PREFIX
        : KLS_FAST_REJECT_REFRESH_UNKNOWN;
    kls_record_fast_reject_detail(solver, shared->rejected_pivot,
                                  shared->rejected_pivot_col,
                                  shared->rejected_row,
                                  shared->rejected_multiplier_abs,
                                  shared->rejected_pivot_abs,
                                  shared->rejected_candidate_abs);
    common->status = TRILINOS_KLU_OK;
    return 0;
  }
  if (shared->singular) {
    common->status = TRILINOS_KLU_SINGULAR;
    common->numerical_rank = shared->numerical_rank;
    common->singular_col = shared->singular_col;
    if (common->halt_if_singular) {
      return 0;
    }
    if (common->scale > 0 && !kls_parallel_refactor_permute_scale(solver)) {
      common->status = TRILINOS_KLU_INVALID;
      return 0;
    }
    return 1;
  }
  common->status = TRILINOS_KLU_OK;
  common->numerical_rank = KLS_KLU_EMPTY;
  common->singular_col = KLS_KLU_EMPTY;
  common->nrealloc = 0;
  if (common->scale > 0 && !kls_parallel_refactor_permute_scale(solver)) {
    common->status = TRILINOS_KLU_INVALID;
    return 0;
  }
  return 1;
}

static int kls_mapped_refactor(kls_solver *solver,
                               double *numeric_values,
                               int check_pivots) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      numeric_values == NULL) {
    return -1;
  }
  const int scaled = solver->common.scale > 0;
  if (scaled) {
    if (solver->symbolic->nblocks == 1u ||
        solver->numeric->Rs == NULL || solver->numeric->Pnum == NULL ||
        solver->numeric->Xwork == NULL) {
      return -1;
    }
    if (!kls_scaled_serial_mapped_btf_is_worthwhile(solver)) {
      return -1;
    }
  } else if (solver->numeric->Rs != NULL) {
    return -1;
  }
  if (!kls_build_refactor_map(solver)) {
    return -1;
  }
  if (solver->symbolic->nblocks == 1u) {
    if (check_pivots && kls_checked_row_refactor_env_enabled()) {
      if (solver->options.threads > 1) {
        const int parallel_row_status =
          kls_single_block_parallel_row_refactor(solver, numeric_values, 1);
        if (parallel_row_status >= 0) {
          return parallel_row_status;
        }
      }
      const int row_status =
        kls_single_block_row_refactor(solver, numeric_values, 1);
      if (row_status >= 0) {
        return row_status;
      }
    } else if (!check_pivots && kls_row_refactor_env_enabled()) {
      if (solver->options.threads > 1) {
        const int parallel_row_status =
          kls_single_block_parallel_row_refactor(solver, numeric_values, 0);
        if (parallel_row_status >= 0) {
          return parallel_row_status;
        }
      }
      const int row_status =
        kls_single_block_row_refactor(solver, numeric_values, 0);
      if (row_status >= 0) {
        return row_status;
      }
    }
    return kls_single_block_mapped_refactor(solver, numeric_values, check_pivots);
  }
  if (solver->numeric->Udiag == NULL || solver->numeric->Offp == NULL ||
      solver->numeric->Offx == NULL || solver->numeric->Lip == NULL ||
      solver->numeric->Llen == NULL || solver->numeric->Uip == NULL ||
      solver->numeric->Ulen == NULL || solver->numeric->LUbx == NULL ||
      solver->numeric->Xwork == NULL) {
    return -1;
  }

  trilinos_klu_l_common *common = &solver->common;
  if (scaled &&
      !trilinos_klu_l_scale((UF_long)common->scale, solver->n,
                            solver->col_ptr, solver->row_idx,
                            numeric_values, solver->numeric->Rs, NULL,
                            common)) {
    return 0;
  }

  common->status = TRILINOS_KLU_OK;
  common->numerical_rank = KLS_KLU_EMPTY;
  common->singular_col = KLS_KLU_EMPTY;
  common->nrealloc = 0;

  kls_parallel_refactor_shared shared;
  memset(&shared, 0, sizeof(shared));
  shared.col_ptr = solver->col_ptr;
  shared.row_idx = solver->row_idx;
  shared.nnz = solver->nnz;
  shared.map_col_ptr = solver->refactor_col_ptr;
  shared.map_row_idx = solver->refactor_row_idx;
  shared.map_input_pos = solver->refactor_input_pos;
  shared.map_block_start = solver->refactor_block_start;
  shared.values = numeric_values;
  shared.symbolic = solver->symbolic;
  shared.numeric = solver->numeric;
  shared.rs = scaled ? solver->numeric->Rs : NULL;
  shared.n = solver->n;
  shared.scale = (int)common->scale;
  shared.halt_if_singular = common->halt_if_singular;
  shared.check_pivots = check_pivots;
  shared.pivot_tolerance = common->tol;

  kls_parallel_refactor_worker worker;
  memset(&worker, 0, sizeof(worker));
  worker.shared = &shared;
  worker.x = (double *)solver->numeric->Xwork;

  for (UF_long block = 0; block < solver->symbolic->nblocks; ++block) {
    worker.rejected_pivot = KLS_KLU_EMPTY;
    worker.rejected_pivot_col = KLS_KLU_EMPTY;
    worker.rejected_row = KLS_KLU_EMPTY;
    worker.rejected_multiplier_abs = -1.0;
    worker.rejected_pivot_abs = -1.0;
    worker.rejected_candidate_abs = -1.0;
    kls_parallel_refactor_block(&worker, block);
    if (worker.invalid || worker.pivot_rejected ||
        (worker.singular && common->halt_if_singular)) {
      memset(worker.x, 0, (size_t)solver->symbolic->maxblock * sizeof(*worker.x));
      break;
    }
  }

  if (worker.invalid) {
    common->status = TRILINOS_KLU_INVALID;
    return 0;
  }
  if (worker.pivot_rejected) {
    solver->fast_reject_refresh_state = KLS_FAST_REJECT_REFRESH_PREFIX;
    kls_record_fast_reject_detail(solver, worker.rejected_pivot,
                                  worker.rejected_pivot_col,
                                  worker.rejected_row,
                                  worker.rejected_multiplier_abs,
                                  worker.rejected_pivot_abs,
                                  worker.rejected_candidate_abs);
    common->status = TRILINOS_KLU_OK;
    return 0;
  }
  if (worker.singular) {
    common->status = TRILINOS_KLU_SINGULAR;
    common->numerical_rank = worker.numerical_rank;
    common->singular_col = worker.singular_col;
    if (common->halt_if_singular) {
      return 0;
    }
  }
  if (scaled && !kls_parallel_refactor_permute_scale(solver)) {
    common->status = TRILINOS_KLU_INVALID;
    return 0;
  }
  if (!worker.singular) {
    common->status = TRILINOS_KLU_OK;
  }
  return 1;
}

static int kls_serial_refactor_tail_from_block(kls_solver *solver,
                                               double *numeric_values,
                                               UF_long start_block,
                                               int check_pivots) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      numeric_values == NULL || solver->common.scale > 0 ||
      solver->col_ptr == NULL || solver->row_idx == NULL ||
      solver->symbolic->R == NULL || solver->symbolic->Q == NULL ||
      solver->numeric->Rs != NULL ||
      solver->numeric->Udiag == NULL || solver->numeric->Offp == NULL ||
      solver->numeric->Offx == NULL || solver->numeric->Lip == NULL ||
      solver->numeric->Llen == NULL || solver->numeric->Uip == NULL ||
      solver->numeric->Ulen == NULL || solver->numeric->LUbx == NULL ||
      solver->numeric->Xwork == NULL || solver->numeric->Pinv == NULL ||
      start_block > solver->symbolic->nblocks) {
    return -1;
  }
  if (start_block == solver->symbolic->nblocks) {
    return 1;
  }

  kls_parallel_refactor_shared shared;
  memset(&shared, 0, sizeof(shared));
  shared.col_ptr = solver->col_ptr;
  shared.row_idx = solver->row_idx;
  shared.nnz = solver->nnz;
  shared.values = numeric_values;
  shared.symbolic = solver->symbolic;
  shared.numeric = solver->numeric;
  shared.rs = NULL;
  shared.n = solver->n;
  shared.scale = 0;
  shared.halt_if_singular = solver->common.halt_if_singular;
  shared.check_pivots = check_pivots;
  shared.pivot_tolerance = solver->common.tol;

  kls_parallel_refactor_worker worker;
  memset(&worker, 0, sizeof(worker));
  worker.shared = &shared;
  worker.x = (double *)solver->numeric->Xwork;
  worker.rejected_pivot = KLS_KLU_EMPTY;
  worker.rejected_pivot_col = KLS_KLU_EMPTY;
  worker.rejected_row = KLS_KLU_EMPTY;
  worker.rejected_multiplier_abs = -1.0;
  worker.rejected_pivot_abs = -1.0;
  worker.rejected_candidate_abs = -1.0;
  worker.numerical_rank = UF_long_max;
  worker.singular_col = KLS_KLU_EMPTY;
  memset(worker.x, 0, (size_t)solver->symbolic->maxblock * sizeof(*worker.x));

  trilinos_klu_l_common *common = &solver->common;
  common->status = TRILINOS_KLU_OK;
  common->numerical_rank = KLS_KLU_EMPTY;
  common->singular_col = KLS_KLU_EMPTY;
  common->nrealloc = 0;

  for (UF_long block = start_block; block < solver->symbolic->nblocks; ++block) {
    kls_parallel_refactor_block(&worker, block);
    if (worker.invalid || worker.pivot_rejected ||
        (worker.singular && common->halt_if_singular)) {
      memset(worker.x, 0, (size_t)solver->symbolic->maxblock * sizeof(*worker.x));
      break;
    }
  }

  if (worker.invalid) {
    common->status = TRILINOS_KLU_INVALID;
    return 0;
  }
  if (worker.pivot_rejected) {
    solver->fast_reject_refresh_state = KLS_FAST_REJECT_REFRESH_PREFIX;
    kls_record_fast_reject_detail(solver, worker.rejected_pivot,
                                  worker.rejected_pivot_col,
                                  worker.rejected_row,
                                  worker.rejected_multiplier_abs,
                                  worker.rejected_pivot_abs,
                                  worker.rejected_candidate_abs);
    common->status = TRILINOS_KLU_OK;
    return 0;
  }
  if (worker.singular) {
    common->status = TRILINOS_KLU_SINGULAR;
    common->numerical_rank = worker.numerical_rank;
    common->singular_col = worker.singular_col;
    if (common->halt_if_singular) {
      return 0;
    }
  }
  if (!worker.singular) {
    common->status = TRILINOS_KLU_OK;
  }
  return 1;
}

static int kls_serial_checked_scaled_refactor_from_block(
  kls_solver *solver,
  double *numeric_values,
  UF_long start_block) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      numeric_values == NULL || solver->common.scale <= 0 ||
      solver->numeric->Rs == NULL || solver->numeric->Pnum == NULL ||
      solver->numeric->Pinv == NULL || solver->numeric->Udiag == NULL ||
      solver->numeric->Offp == NULL || solver->numeric->Offx == NULL ||
      solver->numeric->Lip == NULL || solver->numeric->Llen == NULL ||
      solver->numeric->Uip == NULL || solver->numeric->Ulen == NULL ||
      solver->numeric->LUbx == NULL || solver->numeric->Xwork == NULL ||
      start_block > solver->symbolic->nblocks) {
    return -1;
  }

  trilinos_klu_l_common *common = &solver->common;
  if (!trilinos_klu_l_scale((UF_long)common->scale, solver->n,
                            solver->col_ptr, solver->row_idx,
                            numeric_values, solver->numeric->Rs, NULL,
                            common)) {
    return 0;
  }

  common->status = TRILINOS_KLU_OK;
  common->numerical_rank = KLS_KLU_EMPTY;
  common->singular_col = KLS_KLU_EMPTY;
  common->nrealloc = 0;

  kls_parallel_refactor_shared shared;
  memset(&shared, 0, sizeof(shared));
  shared.col_ptr = solver->col_ptr;
  shared.row_idx = solver->row_idx;
  shared.nnz = solver->nnz;
  shared.values = numeric_values;
  shared.symbolic = solver->symbolic;
  shared.numeric = solver->numeric;
  shared.rs = solver->numeric->Rs;
  shared.n = solver->n;
  shared.scale = (int)common->scale;
  shared.halt_if_singular = common->halt_if_singular;
  shared.check_pivots = 1;
  shared.pivot_tolerance = common->tol;

  kls_parallel_refactor_worker worker;
  memset(&worker, 0, sizeof(worker));
  worker.shared = &shared;
  worker.x = (double *)solver->numeric->Xwork;
  worker.rejected_pivot = KLS_KLU_EMPTY;
  worker.rejected_pivot_col = KLS_KLU_EMPTY;
  worker.rejected_row = KLS_KLU_EMPTY;
  worker.rejected_multiplier_abs = -1.0;
  worker.rejected_pivot_abs = -1.0;
  worker.rejected_candidate_abs = -1.0;
  worker.numerical_rank = UF_long_max;
  worker.singular_col = KLS_KLU_EMPTY;
  memset(worker.x, 0, (size_t)solver->symbolic->maxblock * sizeof(*worker.x));

  for (UF_long block = start_block; block < solver->symbolic->nblocks; ++block) {
    kls_parallel_refactor_block(&worker, block);
    if (worker.invalid || worker.pivot_rejected ||
        (worker.singular && common->halt_if_singular)) {
      memset(worker.x, 0, (size_t)solver->symbolic->maxblock * sizeof(*worker.x));
      break;
    }
  }

  if (worker.invalid) {
    common->status = TRILINOS_KLU_INVALID;
    return 0;
  }
  if (worker.pivot_rejected) {
    solver->fast_reject_refresh_state = KLS_FAST_REJECT_REFRESH_PREFIX;
    kls_record_fast_reject_detail(solver, worker.rejected_pivot,
                                  worker.rejected_pivot_col,
                                  worker.rejected_row,
                                  worker.rejected_multiplier_abs,
                                  worker.rejected_pivot_abs,
                                  worker.rejected_candidate_abs);
    common->status = TRILINOS_KLU_OK;
    return 0;
  }
  if (worker.singular) {
    common->status = TRILINOS_KLU_SINGULAR;
    common->numerical_rank = worker.numerical_rank;
    common->singular_col = worker.singular_col;
    if (common->halt_if_singular) {
      return 0;
    }
  }
  if (!kls_parallel_refactor_permute_scale(solver)) {
    common->status = TRILINOS_KLU_INVALID;
    return 0;
  }
  if (!worker.singular) {
    common->status = TRILINOS_KLU_OK;
  }
  return 1;
}

static int kls_serial_checked_scaled_refactor(kls_solver *solver,
                                              double *numeric_values) {
  return kls_serial_checked_scaled_refactor_from_block(solver, numeric_values,
                                                       0u);
}

static void maybe_prepare_refactor_map(kls_solver *solver,
                                       double *elapsed) {
  const int pool_map = kls_refactor_pool_map_is_worthwhile(solver);
  if (solver == NULL || elapsed == NULL ||
      solver->refactor_col_ptr != NULL ||
      solver->numeric == NULL ||
      (!pool_map && (solver->common.scale > 0 ||
                     solver->numeric->Rs != NULL ||
                     kls_parallel_refactor_is_eligible(solver))) ||
      !kls_refactor_map_is_eligible(solver)) {
    return;
  }
  const double start = kls_now_seconds();
  (void)kls_build_refactor_map(solver);
  *elapsed += kls_now_seconds() - start;
}

static int kls_refactor_schedule_is_eligible(const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      solver->options.threads <= 1 || solver->n == 0) {
    return 0;
  }
  if (solver->numeric->Uip == NULL || solver->numeric->Ulen == NULL ||
      solver->numeric->LUbx == NULL || solver->symbolic->R == NULL) {
    return 0;
  }
  const int single_block = solver->symbolic->nblocks == 1u;
  const int dominant_btf = kls_egraph_dominant_btf_shape(solver);
  if (kls_egraph_all_pipeline_huge_single_shape(solver)) {
    return solver->common.flops >= 1.0e9;
  }
  if (kls_egraph_non_dominant_many_block_shape(solver)) {
    return solver->common.flops >= 5.0e8;
  }
  if (kls_egraph_all_pipeline_dominant_btf_shape(solver)) {
    return solver->common.flops >= 1.0e8;
  }
  if (!single_block && !dominant_btf) {
    return 0;
  }
  if (kls_egraph_medium_heavy_dominant_btf_shape(solver)) {
    return solver->common.flops >= 1.5e8;
  }
  if (kls_egraph_scaled_medium_dominant_btf_shape(solver)) {
    return 1;
  }
  if (kls_egraph_small_compact_dominant_btf_shape(solver)) {
    return 1;
  }
  if (kls_egraph_low_work_dominant_btf_shape(solver)) {
    return 1;
  }
  if (dominant_btf) {
    return solver->common.flops >= 2.0e7;
  }
  if (kls_egraph_moderate_single_block_shape(solver)) {
    return 1;
  }
  if (kls_egraph_low_work_single_block_shape(solver)) {
    return 1;
  }
  if (solver->common.flops >= 3.0e8) {
    return 1;
  }
  return (solver->numeric->unz + solver->numeric->lnz) >= 3000000u;
}

static int kls_build_refactor_schedule(kls_solver *solver) {
  if (!kls_refactor_schedule_is_eligible(solver)) {
    return 0;
  }
  const int all_pipeline =
    kls_egraph_all_pipeline_dominant_btf_shape(solver) ||
    kls_egraph_all_pipeline_huge_single_shape(solver);
  const int natural_pipeline =
    kls_egraph_all_pipeline_huge_single_shape(solver);
  if (solver->refactor_level_ptr != NULL &&
      (solver->refactor_level_cols != NULL || natural_pipeline)) {
    return 1;
  }

  free_refactor_schedule(solver);
  UF_long *levels = NULL;
  UF_long *successor_counts = NULL;
  UF_long *first_successor = NULL;
  uint64_t *row_hash1 = NULL;
  uint64_t *row_hash2 = NULL;
  double *column_work = NULL;
  UF_long *counts = NULL;
  UF_long *level_ptr = NULL;
  UF_long *level_cols = NULL;
  UF_long *level_thread_ptr = NULL;

#define KLS_FREE_REFACTOR_SCHEDULE_TEMP() \
  do { \
    free(levels); \
    free(successor_counts); \
    free(first_successor); \
    free(row_hash1); \
    free(row_hash2); \
    free(column_work); \
    free(counts); \
    free(level_ptr); \
    free(level_cols); \
    free(level_thread_ptr); \
  } while (0)

  levels = (UF_long *)calloc((size_t)solver->n, sizeof(*levels));
  successor_counts =
    (UF_long *)calloc((size_t)solver->n, sizeof(*successor_counts));
  first_successor =
    (UF_long *)malloc((size_t)solver->n * sizeof(*first_successor));
  row_hash1 = (uint64_t *)calloc((size_t)solver->n, sizeof(*row_hash1));
  row_hash2 = (uint64_t *)calloc((size_t)solver->n, sizeof(*row_hash2));
  column_work =
    (double *)calloc((size_t)solver->n, sizeof(*column_work));
  if (levels == NULL || successor_counts == NULL ||
      first_successor == NULL || row_hash1 == NULL ||
      row_hash2 == NULL || column_work == NULL) {
    KLS_FREE_REFACTOR_SCHEDULE_TEMP();
    return 0;
  }
  for (UF_long k = 0; k < solver->n; ++k) {
    first_successor[k] = KLS_KLU_EMPTY;
  }

  UF_long max_level = 0;
  UF_long edges = 0;
  UF_long root_columns = 0;
  double max_column_work = 0.0;
  double total_work = 0.0;
  const int weight_singleton_blocks =
    kls_egraph_non_dominant_many_block_shape(solver);
  for (UF_long block = 0; block < solver->symbolic->nblocks; ++block) {
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    const UF_long nk = k2 - k1;
    if (nk <= 1u) {
      if (nk == 1u && weight_singleton_blocks) {
        column_work[k1] = 1.0;
        if (max_column_work < 1.0) {
          max_column_work = 1.0;
        }
      }
      if (nk == 1u) {
        root_columns++;
      }
      continue;
    }
    double *lu = (double *)solver->numeric->LUbx[block];
    if (lu == NULL) {
      KLS_FREE_REFACTOR_SCHEDULE_TEMP();
      return 0;
    }
    const UF_long *uip = solver->numeric->Uip + k1;
    const UF_long *ulen = solver->numeric->Ulen + k1;
    for (UF_long k = 0; k < nk; ++k) {
      UF_long *ui = NULL;
      double *ux = NULL;
      UF_long ucol_len = 0;
      kls_klu_get_pointer(lu, uip, ulen, k, &ui, &ux, &ucol_len);
      (void)ux;
      UF_long level = 0;
      double work = 1.0;
      for (UF_long p = 0; p < ucol_len; ++p) {
        const UF_long dep = ui[p];
        if (dep >= k) {
          KLS_FREE_REFACTOR_SCHEDULE_TEMP();
          return 0;
        }
        const UF_long global_dep = k1 + dep;
        const UF_long global_col = k1 + k;
        successor_counts[global_dep]++;
        if (first_successor[global_dep] == KLS_KLU_EMPTY) {
          first_successor[global_dep] = global_col;
        }
        row_hash1[global_dep] += kls_supernode_hash1(global_col);
        row_hash2[global_dep] += kls_supernode_hash2(global_col);
        const UF_long dep_level = levels[global_dep] + 1u;
        if (dep_level > level) {
          level = dep_level;
        }
        work += 1.0 + (double)solver->numeric->Llen[k1 + dep];
        edges++;
      }
      if (ucol_len == 0u) {
        root_columns++;
      }
      levels[k1 + k] = level;
      column_work[k1 + k] = work;
      if (work > max_column_work) {
        max_column_work = work;
      }
      total_work += work;
      if (level > max_level) {
        max_level = level;
      }
    }
  }

  UF_long leaf_columns = 0;
  UF_long max_fanout = 0;
  for (UF_long k = 0; k < solver->n; ++k) {
    const UF_long fanout = successor_counts[k];
    if (fanout == 0u) {
      leaf_columns++;
    }
    if (fanout > max_fanout) {
      max_fanout = fanout;
    }
  }

  UF_long level_count = max_level + 1u;
  counts = (UF_long *)calloc((size_t)level_count + 1u, sizeof(*counts));
  level_ptr =
    (UF_long *)malloc(((size_t)level_count + 1u) * sizeof(*level_ptr));
  const int thread_count = solver->options.threads;
  level_cols = natural_pipeline ? NULL
                                : (UF_long *)malloc((size_t)solver->n *
                                                    sizeof(*level_cols));
  if (thread_count > 1 && !all_pipeline) {
    level_thread_ptr =
      (UF_long *)malloc((size_t)level_count *
                        ((size_t)thread_count + 1u) *
                        sizeof(*level_thread_ptr));
  }
  if (counts == NULL || level_ptr == NULL ||
      (!natural_pipeline && level_cols == NULL)) {
    KLS_FREE_REFACTOR_SCHEDULE_TEMP();
    return 0;
  }
  if (thread_count > 1 && !all_pipeline && level_thread_ptr == NULL) {
    KLS_FREE_REFACTOR_SCHEDULE_TEMP();
    return 0;
  }

  for (UF_long k = 0; k < solver->n; ++k) {
    if (levels[k] >= level_count) {
      KLS_FREE_REFACTOR_SCHEDULE_TEMP();
      return 0;
    }
    counts[levels[k] + 1u]++;
  }
  for (UF_long level = 0; level < level_count; ++level) {
    counts[level + 1u] += counts[level];
    level_ptr[level] = counts[level];
  }
  level_ptr[level_count] = counts[level_count];
  UF_long max_width = 0;
  UF_long cluster_levels = level_count;
  /* CKTSO uses a width threshold of alpha * threads with alpha=2.  Compact
     dominant-BTF schedules around one very large block retain enough tail work
     to benefit from a slightly earlier pipeline split in KLS's exact EGraph. */
  const double cluster_width_limit =
    (kls_egraph_compact_large_dominant_btf_shape(solver) ? 4.0 : 2.0) *
    (double)solver->options.threads;
  for (UF_long level = 0; level < level_count; ++level) {
    const UF_long width = level_ptr[level + 1u] - level_ptr[level];
    if (width > max_width) {
      max_width = width;
    }
    if (!all_pipeline &&
        cluster_levels == level_count &&
        (double)width < cluster_width_limit) {
      cluster_levels = level;
    }
  }
  if (all_pipeline) {
    cluster_levels = 0u;
  }

  if (level_cols != NULL) {
    UF_long *next = counts;
    for (UF_long k = 0; k < solver->n; ++k) {
      const UF_long level = levels[k];
      const UF_long dst = next[level]++;
      if (dst >= solver->n) {
        KLS_FREE_REFACTOR_SCHEDULE_TEMP();
        return 0;
      }
      level_cols[dst] = k;
    }
  }

  if (level_thread_ptr != NULL) {
    for (UF_long level = 0; level < level_count; ++level) {
      const UF_long begin = level_ptr[level];
      const UF_long end = level_ptr[level + 1u];
      UF_long *parts =
        level_thread_ptr + level * (UF_long)(thread_count + 1);
      parts[0] = begin;
      parts[thread_count] = end;
      /* Cluster levels keep barriers, but split contiguous slices by the
       * existing no-pivot column-work estimate instead of column count. */
      double level_work = 0.0;
      for (UF_long pos = begin; pos < end; ++pos) {
        const UF_long col = level_cols[pos];
        if (col >= solver->n) {
          KLS_FREE_REFACTOR_SCHEDULE_TEMP();
          return 0;
        }
        level_work += column_work[col];
      }
      UF_long pos = begin;
      double prefix_work = 0.0;
      for (int t = 1; t < thread_count; ++t) {
        const double target =
          level_work * (double)t / (double)thread_count;
        while (pos < end && prefix_work < target) {
          prefix_work += column_work[level_cols[pos]];
          pos++;
        }
        parts[t] = pos;
      }
    }
  }

  UF_long pipeline_columns = 0;
  double pipeline_work = 0.0;
  double pipeline_max_column_work = 0.0;
  if (cluster_levels < level_count) {
    for (UF_long pos = level_ptr[cluster_levels]; pos < solver->n; ++pos) {
      const UF_long col = natural_pipeline ? pos : level_cols[pos];
      if (!natural_pipeline && col >= solver->n) {
        KLS_FREE_REFACTOR_SCHEDULE_TEMP();
        return 0;
      }
      pipeline_columns++;
      const double work = column_work[col];
      pipeline_work += work;
      if (work > pipeline_max_column_work) {
        pipeline_max_column_work = work;
      }
    }
  }

  for (UF_long block = 0; block < solver->symbolic->nblocks; ++block) {
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    if (k2 - k1 <= 1u) {
      continue;
    }
    UF_long start = k1;
    UF_long width = 1u;
    for (UF_long k = k1; k + 1u < k2; ++k) {
      const UF_long next = k + 1u;
      const uint64_t next_hash1 = kls_supernode_hash1(next);
      const uint64_t next_hash2 = kls_supernode_hash2(next);
      const int extends =
        first_successor[k] == next &&
        successor_counts[k] > 0u &&
        successor_counts[k] - 1u == successor_counts[next] &&
        row_hash1[k] - next_hash1 == row_hash1[next] &&
        row_hash2[k] - next_hash2 == row_hash2[next];
      if (extends) {
        width++;
      } else {
        kls_record_supernode_candidate(solver, start, width,
                                       successor_counts);
        start = next;
        width = 1u;
      }
    }
    kls_record_supernode_candidate(solver, start, width, successor_counts);
  }

  free(levels);
  free(successor_counts);
  free(first_successor);
  free(row_hash1);
  free(row_hash2);
  free(column_work);
  free(counts);
#undef KLS_FREE_REFACTOR_SCHEDULE_TEMP
  solver->refactor_level_ptr = level_ptr;
  solver->refactor_level_cols = level_cols;
  solver->refactor_level_thread_ptr = level_thread_ptr;
  solver->refactor_level_thread_count = thread_count;
  solver->refactor_level_count = level_count;
  solver->refactor_level_max_width = max_width;
  solver->refactor_dependency_edges = edges;
  solver->refactor_dependency_root_columns = root_columns;
  solver->refactor_dependency_leaf_columns = leaf_columns;
  solver->refactor_dependency_max_fanout = max_fanout;
  solver->refactor_dependency_max_column_work = max_column_work;
  solver->refactor_dependency_pipeline_max_column_work =
    pipeline_max_column_work;
  solver->refactor_cluster_level_count = cluster_levels;
  solver->refactor_pipeline_column_count = pipeline_columns;
  solver->refactor_dependency_work = total_work;
  solver->refactor_pipeline_work = pipeline_work;
  return 1;
}

static void maybe_prepare_refactor_schedule(kls_solver *solver,
                                            double *elapsed) {
  if (solver == NULL || elapsed == NULL ||
      solver->refactor_level_ptr != NULL ||
      !kls_refactor_schedule_is_eligible(solver)) {
    return;
  }
  const double start = kls_now_seconds();
  (void)kls_build_refactor_schedule(solver);
  *elapsed += kls_now_seconds() - start;
}

static UF_long kls_parallel_refactor(kls_solver *solver,
                                     double *numeric_values,
                                     int check_pivots) {
  solver->fast_reject_refresh_state = KLS_FAST_REJECT_REFRESH_UNKNOWN;
  if (check_pivots && kls_checked_row_refactor_env_enabled()) {
    if (solver->options.threads > 1) {
      const int parallel_row_status =
        kls_single_block_parallel_row_refactor(solver, numeric_values, 1);
      if (parallel_row_status >= 0) {
        return (UF_long)parallel_row_status;
      }
    }
    const int row_status =
      kls_single_block_row_refactor(solver, numeric_values, 1);
    if (row_status >= 0) {
      return (UF_long)row_status;
    }
  } else if (!check_pivots && kls_row_refactor_env_enabled()) {
    if (solver->options.threads > 1) {
      const int parallel_row_status =
        kls_single_block_parallel_row_refactor(solver, numeric_values, 0);
      if (parallel_row_status >= 0) {
        return (UF_long)parallel_row_status;
      }
    }
    const int row_status =
      kls_single_block_row_refactor(solver, numeric_values, 0);
    if (row_status >= 0) {
      return (UF_long)row_status;
    }
  }
  const int egraph =
    kls_egraph_mapped_refactor(solver, numeric_values, check_pivots);
  if (egraph >= 0) {
    return (UF_long)egraph;
  }

  if (!kls_parallel_refactor_is_eligible(solver)) {
    const int mapped = kls_mapped_refactor(solver, numeric_values, check_pivots);
    if (mapped >= 0) {
      return (UF_long)mapped;
    }
    if (check_pivots && solver->common.scale > 0) {
      const int checked = kls_serial_checked_scaled_refactor(solver,
                                                            numeric_values);
      if (checked >= 0) {
        return (UF_long)checked;
      }
    }
    const UF_long ok =
      trilinos_klu_l_refactor(solver->col_ptr, solver->row_idx,
                              numeric_values, solver->symbolic,
                              solver->numeric, &solver->common);
    if (ok && check_pivots && solver->common.status >= 0 &&
        solver->common.status != TRILINOS_KLU_SINGULAR) {
      UF_long rejected_pivot = KLS_KLU_EMPTY;
      UF_long rejected_pivot_col = KLS_KLU_EMPTY;
      if (!kls_numeric_pivots_pass_threshold(solver, &rejected_pivot,
                                             &rejected_pivot_col)) {
        solver->fast_reject_refresh_state = KLS_FAST_REJECT_REFRESH_ALL;
        kls_record_fast_reject(solver, rejected_pivot, rejected_pivot_col);
        return 0;
      }
    }
    return ok;
  }

  trilinos_klu_l_common *common = &solver->common;
  trilinos_klu_l_symbolic *symbolic = solver->symbolic;

  common->status = TRILINOS_KLU_OK;
  common->numerical_rank = KLS_KLU_EMPTY;
  common->singular_col = KLS_KLU_EMPTY;
  common->nrealloc = 0;
  if (common->scale > 0 &&
      !trilinos_klu_l_scale((UF_long)common->scale, solver->n, solver->col_ptr,
                            solver->row_idx, numeric_values, solver->numeric->Rs,
                            NULL, common)) {
    return 0;
  }

  int thread_count = solver->options.threads;
  if ((UF_long)thread_count > symbolic->nblocks) {
    thread_count = (int)symbolic->nblocks;
  }
  if (thread_count < 2) {
    const UF_long ok =
      trilinos_klu_l_refactor(solver->col_ptr, solver->row_idx,
                              numeric_values, solver->symbolic,
                              solver->numeric, &solver->common);
    if (ok && check_pivots && solver->common.status >= 0 &&
        solver->common.status != TRILINOS_KLU_SINGULAR) {
      UF_long rejected_pivot = KLS_KLU_EMPTY;
      UF_long rejected_pivot_col = KLS_KLU_EMPTY;
      if (!kls_numeric_pivots_pass_threshold(solver, &rejected_pivot,
                                             &rejected_pivot_col)) {
        solver->fast_reject_refresh_state = KLS_FAST_REJECT_REFRESH_ALL;
        kls_record_fast_reject(solver, rejected_pivot, rejected_pivot_col);
        return 0;
      }
    }
    return ok;
  }

  if (kls_refactor_pool_map_is_worthwhile(solver) &&
      solver->refactor_col_ptr == NULL) {
    (void)kls_build_refactor_map(solver);
  }

  int invalid = 0;
  int pivot_rejected = 0;
  UF_long rejected_pivot = KLS_KLU_EMPTY;
  UF_long rejected_pivot_col = KLS_KLU_EMPTY;
  UF_long rejected_row = KLS_KLU_EMPTY;
  double rejected_multiplier_abs = -1.0;
  double rejected_pivot_abs = -1.0;
  double rejected_candidate_abs = -1.0;
  int singular = 0;
  UF_long numerical_rank = UF_long_max;
  UF_long singular_col = KLS_KLU_EMPTY;
  int pool_prefix_current = 0;
  if (!run_refactor_pool(solver, numeric_values, thread_count, check_pivots,
                         &invalid, &pivot_rejected,
                         &rejected_pivot, &rejected_pivot_col,
                         &rejected_row, &rejected_multiplier_abs,
                         &rejected_pivot_abs, &rejected_candidate_abs,
                         &singular,
                         &numerical_rank, &singular_col,
                         &pool_prefix_current)) {
    common->status = TRILINOS_KLU_OUT_OF_MEMORY;
    return 0;
  }

  if (invalid) {
    common->status = TRILINOS_KLU_INVALID;
    return 0;
  }
  if (pivot_rejected) {
    solver->fast_reject_refresh_state =
      pool_prefix_current ? KLS_FAST_REJECT_REFRESH_PREFIX
                          : KLS_FAST_REJECT_REFRESH_UNKNOWN;
    kls_record_fast_reject_detail(solver, rejected_pivot,
                                  rejected_pivot_col, rejected_row,
                                  rejected_multiplier_abs,
                                  rejected_pivot_abs,
                                  rejected_candidate_abs);
    common->status = TRILINOS_KLU_OK;
    return 0;
  }

  if (singular) {
    common->status = TRILINOS_KLU_SINGULAR;
    common->numerical_rank = numerical_rank;
    common->singular_col = singular_col;
    if (common->halt_if_singular) {
      return 0;
    }
  }

  if (common->scale > 0 && !kls_parallel_refactor_permute_scale(solver)) {
    common->status = TRILINOS_KLU_INVALID;
    return 0;
  }

  if (!singular) {
    common->status = TRILINOS_KLU_OK;
  }
  return 1;
}

static UF_long kls_fast_factor_with_block_restarts(kls_solver *solver,
                                                   double *numeric_values) {
  const int max_restarts = 4;
  for (int attempt = 0; attempt <= max_restarts; ++attempt) {
    const UF_long ok = kls_parallel_refactor(solver, numeric_values, 1);
    if (ok || solver->common.status < 0 ||
        solver->common.status == TRILINOS_KLU_SINGULAR) {
      return ok;
    }
    if (solver->stats.fast_rejected_pivot < 0 ||
        solver->numeric == NULL || attempt == max_restarts) {
      return 0;
    }
    const UF_long rejected_pivot = (UF_long)solver->stats.fast_rejected_pivot;
    const int repair_covers_remaining =
      kls_fast_repair_covers_remaining_columns(solver, rejected_pivot);
    if (solver->common.scale > 0) {
      if (solver->fast_reject_refresh_state == KLS_FAST_REJECT_REFRESH_ALL ||
          (!repair_covers_remaining &&
           solver->fast_reject_refresh_state !=
             KLS_FAST_REJECT_REFRESH_PREFIX)) {
        return 0;
      }
    }
    if (!kls_pivot_restart_rejected_block(
          solver, numeric_values, rejected_pivot)) {
      return 0;
    }
    const UF_long repaired_block = kls_block_for_pivot(solver, rejected_pivot);
    if (repaired_block == KLS_KLU_EMPTY) {
      return 0;
    }
    if (repair_covers_remaining) {
      UF_long next_rejected_pivot = KLS_KLU_EMPTY;
      UF_long rejected_pivot_col = KLS_KLU_EMPTY;
      if (!kls_numeric_pivots_pass_threshold_from_block(
            solver, repaired_block, &next_rejected_pivot,
            &rejected_pivot_col)) {
        kls_record_fast_reject(solver, next_rejected_pivot,
                               rejected_pivot_col);
        return 0;
      }
      return 1;
    }
    if (solver->fast_reject_refresh_state == KLS_FAST_REJECT_REFRESH_PREFIX) {
      const int tail_ok =
        solver->common.scale > 0
          ? kls_serial_checked_scaled_refactor_from_block(
              solver, numeric_values, repaired_block + 1u)
          : kls_serial_refactor_tail_from_block(solver, numeric_values,
                                                repaired_block + 1u, 1);
      if (tail_ok < 0 || solver->common.status < 0 ||
          solver->common.status == TRILINOS_KLU_SINGULAR) {
        return 0;
      }
      if (!tail_ok) {
        if (solver->stats.fast_rejected_pivot >= 0) {
          continue;
        }
        return 0;
      }
      UF_long next_rejected_pivot = KLS_KLU_EMPTY;
      UF_long rejected_pivot_col = KLS_KLU_EMPTY;
      if (!kls_numeric_pivots_pass_threshold_from_block(
            solver, repaired_block, &next_rejected_pivot,
            &rejected_pivot_col)) {
        kls_record_fast_reject(solver, next_rejected_pivot,
                               rejected_pivot_col);
        continue;
      }
      return 1;
    }
  }
  return 0;
}

int kls_factor(kls_solver *solver, const double *values) {
  if (solver == NULL || solver->symbolic == NULL || values == NULL) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  kls_clear_fast_reject_stats(solver);
  kls_clear_row_refactor_last_stats(solver);
  double *numeric_values = NULL;
  int status = prepare_numeric_values(solver, values, &numeric_values);
  if (status != KLS_OK) {
    return status;
  }

  double elapsed = 0.0;
  const int had_numeric = solver->numeric != NULL;
  if (!had_numeric) {
    maybe_select_pre_static_row_match(solver, &elapsed, numeric_values);
    if (solver->numeric != NULL) {
      maybe_prepare_refactor_map(solver, &elapsed);
      maybe_prepare_refactor_schedule(solver, &elapsed);
      solver->stats.factor_seconds = elapsed;
      fill_numeric_stats(solver);
      return solver->common.status == TRILINOS_KLU_SINGULAR ? KLS_ERR_SINGULAR
                                                            : KLS_OK;
    }
  }
  if (solver->options.fast_factor && solver->numeric != NULL) {
    const double start = kls_now_seconds();
    const UF_long ok = kls_fast_factor_with_block_restarts(solver,
                                                           numeric_values);
    elapsed += kls_now_seconds() - start;
    if (ok && solver->common.status >= 0 &&
        solver->common.status != TRILINOS_KLU_SINGULAR) {
      solver->stats.factor_seconds = elapsed;
      kls_update_numeric_diagnostics(solver, 1);
      maybe_prepare_refactor_map(solver, &elapsed);
      maybe_prepare_refactor_schedule(solver, &elapsed);
      solver->stats.factor_seconds = elapsed;
      fill_numeric_stats(solver);
      return KLS_OK;
    }
  }

  free_numeric(solver);
  if (!had_numeric) {
    solver->common.scale = choose_auto_scale_from_values(solver, numeric_values);
    solver->common.tol = choose_initial_auto_pivot_tolerance(solver);
  }
  const double start = kls_now_seconds();
  solver->numeric = trilinos_klu_l_factor(solver->col_ptr, solver->row_idx,
                                          numeric_values, solver->symbolic,
                                          &solver->common);
  elapsed += kls_now_seconds() - start;
  solver->stats.factor_seconds = elapsed;
  if (solver->numeric == NULL || solver->common.status < 0) {
    fill_numeric_stats(solver);
    return solver->common.status == TRILINOS_KLU_SINGULAR ? KLS_ERR_SINGULAR
                                                          : KLS_ERR_FACTOR_FAILED;
  }
  kls_update_numeric_diagnostics(solver, 1);
  int diagnostics_have_flops = 1;
  int diagnostics_have_rcond = 1;
  if (maybe_select_auto_row_match(solver, &elapsed, numeric_values)) {
    numeric_values = solver->values != NULL ? solver->values : numeric_values;
    diagnostics_have_flops = 1;
    diagnostics_have_rcond = 1;
  }
  if (maybe_select_auto_scale(solver, &elapsed, numeric_values)) {
    diagnostics_have_flops = 1;
    diagnostics_have_rcond = 0;
  }
#ifdef KLS_HAVE_METIS
  if (maybe_promote_auto_metis(solver, &elapsed, numeric_values)) {
    diagnostics_have_flops = 1;
    diagnostics_have_rcond = 0;
  }
#endif
  if (!diagnostics_have_flops || !diagnostics_have_rcond) {
    kls_update_numeric_diagnostics(solver, 1);
    diagnostics_have_flops = 1;
    diagnostics_have_rcond = 1;
  }
  if (maybe_select_auto_pivot_tolerance(solver, &elapsed, numeric_values)) {
    diagnostics_have_flops = 1;
    diagnostics_have_rcond = 1;
  }
#ifdef KLS_HAVE_SPRAL_SCALING
  if (maybe_select_spral_hungarian_row_match(solver, &elapsed,
                                             numeric_values)) {
    numeric_values = solver->values != NULL ? solver->values : numeric_values;
    diagnostics_have_flops = 1;
    diagnostics_have_rcond = 1;
  }
#endif
  if (!diagnostics_have_flops || !diagnostics_have_rcond) {
    kls_update_numeric_diagnostics(solver, 1);
  }
  maybe_prepare_refactor_map(solver, &elapsed);
  maybe_prepare_refactor_schedule(solver, &elapsed);
  solver->stats.factor_seconds = elapsed;
  fill_numeric_stats(solver);
  if (solver->common.status == TRILINOS_KLU_SINGULAR) {
    return KLS_ERR_SINGULAR;
  }
  return KLS_OK;
}

int kls_refactor(kls_solver *solver, const double *values) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL || values == NULL) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  kls_clear_fast_reject_stats(solver);
  kls_clear_row_refactor_last_stats(solver);
  double *numeric_values = NULL;
  int status = prepare_numeric_values(solver, values, &numeric_values);
  if (status != KLS_OK) {
    return status;
  }
  const double start = kls_now_seconds();
  const UF_long ok = kls_parallel_refactor(solver, numeric_values, 0);
  solver->stats.refactor_seconds = kls_now_seconds() - start;
  fill_numeric_stats(solver);
  if (!ok || solver->common.status < 0) {
    return solver->common.status == TRILINOS_KLU_SINGULAR ? KLS_ERR_SINGULAR
                                                          : KLS_ERR_REFACTOR_FAILED;
  }
  if (solver->common.status == TRILINOS_KLU_SINGULAR) {
    return KLS_ERR_SINGULAR;
  }
  return KLS_OK;
}

static int solve_impl(kls_solver *solver,
                      int transpose,
                      int64_t nrhs,
                      const double *b,
                      int64_t ldb,
                      double *x,
                      int64_t ldx) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      b == NULL || x == NULL || nrhs <= 0) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  if (ldb == 0) {
    ldb = solver->n;
  }
  if (ldx == 0) {
    ldx = solver->n;
  }
  if (ldb < solver->n || ldx < solver->n) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  if (nrhs > (int64_t)UF_long_max || ldx > (int64_t)UF_long_max) {
    return KLS_ERR_INVALID_ARGUMENT;
  }

  const double start = kls_now_seconds();
  const int kernel_transpose =
    (solver->orientation == KLS_ORIENTATION_TRANSPOSE) ? !transpose : transpose;
  const int has_row_scale = solver->row_scale != NULL;
  const int has_col_scale = solver->col_scale != NULL;
  double *perm_workspace = solver->row_perm != NULL
    ? ensure_solve_perm_workspace(solver)
    : NULL;
  if (solver->row_perm != NULL && perm_workspace == NULL) {
    return KLS_ERR_OUT_OF_MEMORY;
  }

  if (solver->row_perm != NULL && !kernel_transpose) {
    for (int64_t rhs = 0; rhs < nrhs; ++rhs) {
      const double *src = b + rhs * ldb;
      double *dst = x + rhs * ldx;
      for (UF_long row = 0; row < solver->n; ++row) {
        const UF_long scaled_row = solver->row_perm[row];
        const double rs = has_row_scale ? solver->row_scale[scaled_row] : 1.0;
        perm_workspace[scaled_row] = src[row] * rs;
      }
      memcpy(dst, perm_workspace, (size_t)solver->n * sizeof(double));
    }
  } else if (!kernel_transpose && has_row_scale) {
    for (int64_t rhs = 0; rhs < nrhs; ++rhs) {
      const double *src = b + rhs * ldb;
      double *dst = x + rhs * ldx;
      for (UF_long row = 0; row < solver->n; ++row) {
        dst[row] = src[row] * solver->row_scale[row];
      }
    }
  } else if (kernel_transpose && has_col_scale) {
    for (int64_t rhs = 0; rhs < nrhs; ++rhs) {
      const double *src = b + rhs * ldb;
      double *dst = x + rhs * ldx;
      for (UF_long row = 0; row < solver->n; ++row) {
        dst[row] = src[row] * solver->col_scale[row];
      }
    }
  } else if (b != x || ldb != ldx) {
    for (int64_t rhs = 0; rhs < nrhs; ++rhs) {
      memmove(x + rhs * ldx, b + rhs * ldb, (size_t)solver->n * sizeof(double));
    }
  }

  const UF_long ok = kernel_transpose
    ? trilinos_klu_l_tsolve(solver->symbolic, solver->numeric, (UF_long)ldx,
                            (UF_long)nrhs, x, &solver->common)
    : trilinos_klu_l_solve(solver->symbolic, solver->numeric, (UF_long)ldx,
                           (UF_long)nrhs, x, &solver->common);
  if (ok && !kernel_transpose && has_col_scale) {
    for (int64_t rhs = 0; rhs < nrhs; ++rhs) {
      double *dst = x + rhs * ldx;
      for (UF_long row = 0; row < solver->n; ++row) {
        dst[row] *= solver->col_scale[row];
      }
    }
  }
  if (ok && solver->row_perm != NULL && kernel_transpose) {
    for (int64_t rhs = 0; rhs < nrhs; ++rhs) {
      double *dst = x + rhs * ldx;
      for (UF_long row = 0; row < solver->n; ++row) {
        const UF_long scaled_row = solver->row_perm[row];
        const double rs = has_row_scale ? solver->row_scale[scaled_row] : 1.0;
        perm_workspace[row] = dst[scaled_row] * rs;
      }
      memcpy(dst, perm_workspace, (size_t)solver->n * sizeof(double));
    }
  } else if (ok && kernel_transpose && has_row_scale) {
    for (int64_t rhs = 0; rhs < nrhs; ++rhs) {
      double *dst = x + rhs * ldx;
      for (UF_long row = 0; row < solver->n; ++row) {
        dst[row] *= solver->row_scale[row];
      }
    }
  }
  solver->stats.solve_seconds = kls_now_seconds() - start;
  solver->stats.last_kernel_status = (int)solver->common.status;
  solver->stats.memory_bytes = solver->common.memusage;
  solver->stats.memory_peak_bytes = solver->common.mempeak;

  if (!ok || solver->common.status < 0) {
    return KLS_ERR_SOLVE_FAILED;
  }
  return KLS_OK;
}

int kls_solve(kls_solver *solver,
              int64_t nrhs,
              const double *b,
              int64_t ldb,
              double *x,
              int64_t ldx) {
  return solve_impl(solver, 0, nrhs, b, ldb, x, ldx);
}

int kls_solve_transpose(kls_solver *solver,
                        int64_t nrhs,
                        const double *b,
                        int64_t ldb,
                        double *x,
                        int64_t ldx) {
  return solve_impl(solver, 1, nrhs, b, ldb, x, ldx);
}

int kls_get_stats(const kls_solver *solver, kls_stats *stats) {
  if (solver == NULL || stats == NULL || stats->struct_size == 0) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  const size_t requested_size = stats->struct_size;
  const size_t copy_size = requested_size < sizeof(kls_stats)
    ? requested_size
    : sizeof(kls_stats);
  memcpy(stats, &solver->stats, copy_size);
  stats->struct_size = sizeof(kls_stats);
  return KLS_OK;
}

const char *kls_status_string(int status) {
  switch (status) {
    case KLS_OK: return "ok";
    case KLS_ERR_INVALID_ARGUMENT: return "invalid argument";
    case KLS_ERR_OUT_OF_MEMORY: return "out of memory";
    case KLS_ERR_ANALYZE_FAILED: return "analysis failed";
    case KLS_ERR_FACTOR_FAILED: return "factorization failed";
    case KLS_ERR_REFACTOR_FAILED: return "refactorization failed";
    case KLS_ERR_SOLVE_FAILED: return "solve failed";
    case KLS_ERR_SINGULAR: return "singular matrix";
    case KLS_ERR_UNSUPPORTED: return "unsupported";
    default: return "unknown status";
  }
}

const char *kls_ordering_name(kls_ordering ordering) {
  switch (ordering) {
    case KLS_ORDERING_AUTO: return "auto";
    case KLS_ORDERING_AMD: return "amd";
    case KLS_ORDERING_COLAMD: return "colamd";
    case KLS_ORDERING_NATURAL: return "natural";
    case KLS_ORDERING_METIS: return "metis";
    case KLS_ORDERING_SCOTCH: return "scotch";
    default: return "unknown";
  }
}

const char *kls_orientation_name(kls_orientation orientation) {
  switch (orientation) {
    case KLS_ORIENTATION_AUTO: return "auto";
    case KLS_ORIENTATION_NORMAL: return "normal";
    case KLS_ORIENTATION_TRANSPOSE: return "transpose";
    default: return "unknown";
  }
}
