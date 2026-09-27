#ifndef KLS_TUNING_PROFILE_H
#define KLS_TUNING_PROFILE_H

#include <stdint.h>

/* Costs are seconds per generated-kernel operation, never CPU defaults.
   Kept separate from hot tuning fields to preserve solver field offsets. */
typedef struct kls_snb_cost_model {
  int enabled;
  double narrow[3], dense[6], panel[4];
} kls_snb_cost_model;

typedef struct kls_tuning_values {
  unsigned row_batch_min_rows, row_batch_max_rows, row_small_sort_max;
  double row_dense_min_work, row_compact_panel_min_work;
  double row_compact_panel_min_work_per_entry, row_native_panel_min_work;
  unsigned row_blocked_trailing_min_rows, row_blocked_trailing_min_cols;
  double row_batch_snode_min_work, row_batch_snode_min_work_per_entry;
  unsigned dense_help_tail_slice, btf_scalar_min_rows, btf_scalar_max_rows;
  uint64_t solve_dense_tail_min_nnz;
  double solve_dense_tail_min_fraction;
  unsigned solve_trapezoid_slices;
  uint64_t solve_parallel_rect_min_nnz, solve_min_entries_per_sync;
  unsigned snode_min_batch, snode_min_batch_work;
  double egraph_min_flops_per_thread;
  uint64_t egraph_min_size;
  unsigned egraph_pool_spin_iters, egraph_worker_spin_iters;
  uint64_t snb_coop_min_doubles;
} kls_tuning_values;

typedef struct kls_tuning_metadata {
  uint64_t profile_id;
  int field_count;
  int threads;
  kls_snb_cost_model snb_cost;
} kls_tuning_metadata;

void kls_tuning_defaults(kls_tuning_values *values);
const char *kls_tuning_build_fingerprint(void);
int kls_tuning_write_host_profile(const char *path, int threads,
                                  const kls_tuning_values *values,
                                  unsigned long long seed, int force,
                                  kls_tuning_metadata *metadata);
int kls_tuning_load_host_profile(const char *path, int threads,
                                 kls_tuning_values *values,
                                 kls_tuning_metadata *metadata);
int kls_tuning_write_host_profile_model(const char *path, int threads,
                                  const kls_tuning_values *values,
                                  const kls_snb_cost_model *model,
                                  unsigned long long seed, int force,
                                  kls_tuning_metadata *metadata);

#endif
