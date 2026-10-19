/* Optional diagnostic link object: --wrap=kls_destroy. Only existing public
 * statistics are read, once per externally destroyed solver. This reports
 * retained structure/final counters, not a per-call trace or speed evidence.
 * Compile against the exact header matching the linked solver archive. */
#define _POSIX_C_SOURCE 200809L
#include "kls/kls.h"
#include <errno.h>
#include <fenv.h>
#include <stdatomic.h>
#include <stdio.h>

_Static_assert(ATOMIC_INT_LOCK_FREE == 2, "Diagnostic IDs require lock-free atomics");
static atomic_uint owner_snapshot_id;
extern void __real_kls_destroy(kls_solver *);

void __wrap_kls_destroy(kls_solver *solver) {
  if (!solver) { __real_kls_destroy(solver); return; }
  const int caller_errno = errno;
  fenv_t saved;
  const int held = feholdexcept(&saved) == 0;
  kls_stats stats = {0};
  stats.struct_size = offsetof(kls_stats, compact_solve_index_bytes) +
                      sizeof(stats.compact_solve_index_bytes);
  const int status = held ? kls_get_stats(solver, &stats) : KLS_ERR_UNSUPPORTED;
  const unsigned id = atomic_fetch_add_explicit(&owner_snapshot_id, 1,
                                               memory_order_relaxed);
  fprintf(stderr,
      "owner_snapshot id=%u env_supported=%d status=%d n=%lld nnz=%lld nblocks=%lld max_block=%lld nnz_l=%lld nnz_u=%lld offdiag_pivots=%lld selected_ordering=%lld selected_orientation=%lld selected_scale=%lld selected_btf=%lld last_refactor_path=%lld compact_solve_index_bytes=%lld row_refactor_run_count=%lld row_refactor_checked_run_count=%lld row_refactor_parallel_run_count=%lld row_refactor_group_count=%lld row_refactor_group_level_count=%lld row_refactor_group_level_max_width=%lld row_refactor_group_dependency_edges=%lld row_refactor_group_root_count=%lld row_refactor_group_leaf_count=%lld row_refactor_group_max_fanout=%lld row_refactor_compact_supernode_update_count=%lld row_refactor_row_solve_run_count=%lld row_refactor_values_dirty=%lld\n",
      id, held, status,
      (long long)stats.n,
      (long long)stats.nnz,
      (long long)stats.nblocks,
      (long long)stats.max_block,
      (long long)stats.nnz_l,
      (long long)stats.nnz_u,
      (long long)stats.offdiag_pivots,
      (long long)stats.selected_ordering,
      (long long)stats.selected_orientation,
      (long long)stats.selected_scale,
      (long long)stats.selected_btf,
      (long long)stats.last_refactor_path,
      (long long)stats.compact_solve_index_bytes,
      (long long)stats.row_refactor_run_count,
      (long long)stats.row_refactor_checked_run_count,
      (long long)stats.row_refactor_parallel_run_count,
      (long long)stats.row_refactor_group_count,
      (long long)stats.row_refactor_group_level_count,
      (long long)stats.row_refactor_group_level_max_width,
      (long long)stats.row_refactor_group_dependency_edges,
      (long long)stats.row_refactor_group_root_count,
      (long long)stats.row_refactor_group_leaf_count,
      (long long)stats.row_refactor_group_max_fanout,
      (long long)stats.row_refactor_compact_supernode_update_count,
      (long long)stats.row_refactor_row_solve_run_count,
      (long long)stats.row_refactor_values_dirty);
  if (held && fesetenv(&saved) != 0)
    fputs("owner_snapshot_restore_error\n", stderr);
  errno = caller_errno;
  __real_kls_destroy(solver);
}
