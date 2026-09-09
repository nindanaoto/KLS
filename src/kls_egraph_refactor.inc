/* Elimination-graph scheduling, dependency execution, and mapped refactor kernels.
 * Included by kls_numeric_engines.inc in the original definition order.
 */

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
  if (shared->solver->zero_pivot_collect != NULL) {
    const long slot = atomic_fetch_add_explicit(
      &shared->solver->zero_pivot_collect_count, 1, memory_order_relaxed);
    if (slot < shared->solver->zero_pivot_collect_cap) {
      shared->solver->zero_pivot_collect[slot] = numerical_rank;
    }
  }
}

static int kls_egraph_refactor_should_stop(
  kls_egraph_refactor_shared *shared) {
  if (atomic_load_explicit(&shared->stop, memory_order_acquire) != 0) {
    return 1;
  }
  if (shared->trial_deadline_seconds != 0.0 &&
      kls_now_seconds() > shared->trial_deadline_seconds) {
    atomic_store_explicit(&shared->trial_deadline_hit, 1,
                          memory_order_release);
    atomic_store_explicit(&shared->stop, 1, memory_order_release);
    return 1;
  }
  return 0;
}

static int kls_egraph_refactor_value(kls_egraph_refactor_shared *shared,
                                     UF_long input_pos,
                                     double *value_out) {
  return shared != NULL &&
         kls_refactor_input_value(shared->solver, shared->values, shared->rs,
                                  shared->scale, input_pos, value_out);
}

static inline UF_long kls_refactor_map_row_at(const kls_solver *solver,
                                              UF_long pos) {
  return solver->refactor_row_idx32 != NULL
    ? (UF_long)solver->refactor_row_idx32[pos]
    : solver->refactor_row_idx[pos];
}

static inline UF_long kls_refactor_map_input_at(const kls_solver *solver,
                                                UF_long pos) {
  return solver->refactor_input_pos32 != NULL
    ? (UF_long)solver->refactor_input_pos32[pos]
    : solver->refactor_input_pos[pos];
}

static inline void kls_egraph_scatter_unscaled_input(
  const kls_solver *solver,
  const double *restrict values,
  double *restrict x,
  UF_long begin,
  UF_long end) {
  if (solver->refactor_row_idx32 != NULL &&
      solver->refactor_input_pos32 != NULL) {
    const int32_t *restrict row_idx32 = solver->refactor_row_idx32;
    const int32_t *restrict input_pos32 = solver->refactor_input_pos32;
    for (UF_long p = begin; p < end; ++p) {
      x[(UF_long)row_idx32[p]] = values[(UF_long)input_pos32[p]];
    }
    return;
  }
  for (UF_long p = begin; p < end; ++p) {
    x[solver->refactor_row_idx[p]] = values[solver->refactor_input_pos[p]];
  }
}

static inline double kls_egraph_unscaled_input_value_at(
  const kls_solver *solver,
  const double *restrict values,
  UF_long pos) {
  return values[solver->refactor_input_pos32 != NULL
                  ? (UF_long)solver->refactor_input_pos32[pos]
                  : solver->refactor_input_pos[pos]];
}

static inline int kls_egraph_copy_unscaled_offblock_input(
  const kls_solver *solver,
  const double *restrict values,
  double *restrict out,
  UF_long *out_pos_io,
  UF_long out_end,
  UF_long begin,
  UF_long end) {
  UF_long out_pos = *out_pos_io;
  /* The symbolic map partitions each input column into one off-block prefix
     and one diagonal-block suffix, while Offp records the same prefix in the
     retained numeric.  Validate that invariant once before copying instead
     of branching on the output bound for every gathered value. */
  if (begin > end || out_pos > out_end ||
      end - begin != out_end - out_pos) {
    return 0;
  }
  if (solver->refactor_input_pos32 != NULL) {
    const int32_t *restrict input_pos32 = solver->refactor_input_pos32;
    for (UF_long p = begin; p < end; ++p) {
      out[out_pos++] = values[(UF_long)input_pos32[p]];
    }
  } else {
    const UF_long *restrict input_pos = solver->refactor_input_pos;
    for (UF_long p = begin; p < end; ++p) {
      out[out_pos++] = values[input_pos[p]];
    }
  }
  *out_pos_io = out_pos;
  return 1;
}

static inline void kls_egraph_scatter_btf_unscaled_input(
  const kls_solver *solver,
  const double *restrict values,
  double *restrict x,
  UF_long k1,
  UF_long begin,
  UF_long end) {
  if (solver->refactor_row_idx32 != NULL &&
      solver->refactor_input_pos32 != NULL) {
    const int32_t *restrict row_idx32 = solver->refactor_row_idx32;
    const int32_t *restrict input_pos32 = solver->refactor_input_pos32;
    for (UF_long p = begin; p < end; ++p) {
      x[(UF_long)row_idx32[p] - k1] = values[(UF_long)input_pos32[p]];
    }
    return;
  }
  for (UF_long p = begin; p < end; ++p) {
    x[solver->refactor_row_idx[p] - k1] = values[solver->refactor_input_pos[p]];
  }
}

static inline int kls_egraph_input_value_from_pos(
  const kls_solver *solver,
  const double *restrict values,
  const double *restrict rs,
  int scale,
  UF_long input_pos,
  double *value_out) {
  if (solver == NULL || values == NULL || value_out == NULL ||
      input_pos >= solver->nnz) {
    return 0;
  }
  double value = values[input_pos];
  if (scale > 0) {
    if (solver->row_idx == NULL || rs == NULL) {
      return 0;
    }
    const UF_long oldrow = solver->row_idx[input_pos];
    if (oldrow >= solver->n || rs[oldrow] == 0.0) {
      return 0;
    }
    value /= rs[oldrow];
  }
  *value_out = value;
  return 1;
}

static inline int kls_egraph_copy_refactor_values(
  const kls_solver *solver,
  const double *restrict values,
  const double *restrict rs,
  int scale,
  double *restrict out,
  UF_long *out_pos_io,
  UF_long out_end,
  UF_long begin,
  UF_long end) {
  UF_long out_pos = *out_pos_io;
  if (solver->refactor_input_pos32 != NULL) {
    const int32_t *restrict input_pos32 = solver->refactor_input_pos32;
    if (scale > 0 && rs != NULL &&
        solver->refactor_input_oldrow32 != NULL) {
      const int32_t *restrict oldrow32 = solver->refactor_input_oldrow32;
      for (UF_long p = begin; p < end; ++p) {
        const UF_long oldrow = (UF_long)oldrow32[p];
        if (out_pos >= out_end || oldrow >= solver->n || rs[oldrow] == 0.0) {
          return 0;
        }
        out[out_pos++] = values[(UF_long)input_pos32[p]] / rs[oldrow];
      }
      *out_pos_io = out_pos;
      return 1;
    }
    for (UF_long p = begin; p < end; ++p) {
      if (out_pos >= out_end) {
        return 0;
      }
      double value = 0.0;
      if (!kls_egraph_input_value_from_pos(
            solver, values, rs, scale, (UF_long)input_pos32[p], &value)) {
        return 0;
      }
      out[out_pos++] = value;
    }
  } else {
    const UF_long *restrict input_pos = solver->refactor_input_pos;
    for (UF_long p = begin; p < end; ++p) {
      if (out_pos >= out_end) {
        return 0;
      }
      double value = 0.0;
      if (!kls_egraph_input_value_from_pos(
            solver, values, rs, scale, input_pos[p], &value)) {
        return 0;
      }
      out[out_pos++] = value;
    }
  }
  *out_pos_io = out_pos;
  return 1;
}

static inline int kls_egraph_scatter_refactor_input(
  const kls_solver *solver,
  const double *restrict values,
  const double *restrict rs,
  int scale,
  double *restrict x,
  UF_long row_base,
  UF_long row_limit,
  UF_long begin,
  UF_long end) {
  if (row_base > row_limit) {
    return 0;
  }
  if (solver->refactor_row_idx32 != NULL &&
      solver->refactor_input_pos32 != NULL) {
    const int32_t *restrict row_idx32 = solver->refactor_row_idx32;
    const int32_t *restrict input_pos32 = solver->refactor_input_pos32;
    if (scale > 0 && rs != NULL &&
        solver->refactor_input_oldrow32 != NULL) {
      const int32_t *restrict oldrow32 = solver->refactor_input_oldrow32;
      for (UF_long p = begin; p < end; ++p) {
        const UF_long row = (UF_long)row_idx32[p];
        const UF_long oldrow = (UF_long)oldrow32[p];
        if (row < row_base || row >= row_limit || oldrow >= solver->n ||
            rs[oldrow] == 0.0) {
          return 0;
        }
        x[row - row_base] = values[(UF_long)input_pos32[p]] / rs[oldrow];
      }
      return 1;
    }
    for (UF_long p = begin; p < end; ++p) {
      const UF_long row = (UF_long)row_idx32[p];
      if (row < row_base || row >= row_limit) {
        return 0;
      }
      double value = 0.0;
      if (!kls_egraph_input_value_from_pos(
            solver, values, rs, scale, (UF_long)input_pos32[p], &value)) {
        return 0;
      }
      x[row - row_base] = value;
    }
  } else {
    const UF_long *restrict row_idx = solver->refactor_row_idx;
    const UF_long *restrict input_pos = solver->refactor_input_pos;
    for (UF_long p = begin; p < end; ++p) {
      const UF_long row = row_idx[p];
      if (row < row_base || row >= row_limit) {
        return 0;
      }
      double value = 0.0;
      if (!kls_egraph_input_value_from_pos(
            solver, values, rs, scale, input_pos[p], &value)) {
        return 0;
      }
      x[row - row_base] = value;
    }
  }
  return 1;
}

/* egraph-side padded panel refresh: values[] is the packed column in
   the numeric's li order (the padded slot maps were built from that
   same order). Single-block kernels only - local row frame == global. */
static inline void kls_padded_refresh_column(const kls_solver *solver,
                                             UF_long column,
                                             const double *values,
                                             UF_long length) {
  if (solver == NULL || solver->padded_run_of == NULL) {
    return;
  }
  /* refresh runs on EVERY refactor while panels exist - a skipped
     refresh leaves stale panels that a later consume would read as
     wrong values under changing-value SPICE (invisible to the
     identical-values bench); padded_active gates only the consume */
  const UF_long run1 = solver->padded_run_of[column];
  if (run1 == 0u) {
    return;
  }
  const UF_long run = run1 - 1u;
  const UF_long ulen_run =
    solver->padded_union_ptr[run + 1u] - solver->padded_union_ptr[run];
  double *panel_row = solver->padded_panel_values +
    solver->padded_panel_ptr[run] +
    (column - solver->padded_run_start[run]) * ulen_run;
  const UF_long *slots =
    solver->padded_slots + solver->padded_slot_ptr[column];
  for (UF_long p = 0; p < length; ++p) {
    panel_row[slots[p]] = values[p];
  }
}

static inline void kls_egraph_store_l_column_from_workspace(
  const kls_solver *solver,
  double *restrict x,
  UF_long column,
  const UF_long *restrict rows,
  double *restrict values,
  UF_long length,
  double pivot) {
  const int32_t *rows32 =
    solver != NULL && solver->refactor_l_indices32 != NULL
      ? solver->refactor_l_indices32[column]
      : NULL;
  if (rows32 != NULL || length == 0u) {
#if KLS_HAVE_AVX512_KERNELS
    /* Extraction has more arithmetic to amortize an indexed vector access
       than the ordinary sparse RMW kernel: eight scalar divisions become
       one vector divide while the gathered accumulator lanes are also
       cleared by one scatter.  Long realized columns therefore cross over
       well before the update-only scatter's conservative 256-entry gate.
       Keep short columns scalar and select solely from the available ISA
       and actual stream length. */
    if (length >= 64u && kls_avx512_scatter_enabled()) {
      kls_store_l_i32_avx512(x, rows32, values, length, pivot);
    } else
#endif
    {
      for (UF_long p = 0; p < length; ++p) {
        const UF_long i = (UF_long)rows32[p];
        values[p] = x[i] / pivot;
        x[i] = 0.0;
      }
    }
    if (solver != NULL && solver->refactor_l_sorted_values != NULL &&
        solver->refactor_l_sorted_pos32 != NULL) {
      double *sorted_values = solver->refactor_l_sorted_values[column];
      const int32_t *sorted_pos = solver->refactor_l_sorted_pos32[column];
      if (sorted_values != NULL && sorted_pos != NULL) {
        for (UF_long p = 0u; p < length; ++p) {
          sorted_values[p] = values[(UF_long)sorted_pos[p]];
        }
      }
    }
    kls_padded_refresh_column(solver, column, values, length);
    return;
  }
  for (UF_long p = 0; p < length; ++p) {
    const UF_long i = rows[p];
    values[p] = x[i] / pivot;
    x[i] = 0.0;
  }
  if (solver != NULL && solver->refactor_l_sorted_values != NULL &&
      solver->refactor_l_sorted_pos32 != NULL) {
    double *sorted_values = solver->refactor_l_sorted_values[column];
    const int32_t *sorted_pos = solver->refactor_l_sorted_pos32[column];
    if (sorted_values != NULL && sorted_pos != NULL) {
      for (UF_long p = 0u; p < length; ++p) {
        sorted_values[p] = values[(UF_long)sorted_pos[p]];
      }
    }
  }
  kls_padded_refresh_column(solver, column, values, length);
}

static KLS_ALWAYS_INLINE UF_long kls_egraph_llen_at(
  const kls_solver *solver,
  UF_long column) {
  return solver->i32solve_llen != NULL
    ? (UF_long)solver->i32solve_llen[column]
    : solver->numeric->Llen[column];
}

static KLS_ALWAYS_INLINE UF_long kls_egraph_ulen_at(
  const kls_solver *solver,
  UF_long column) {
  return solver->i32solve_ulen != NULL
    ? (UF_long)solver->i32solve_ulen[column]
    : solver->numeric->Ulen[column];
}

static int kls_egraph_refactor_mark_done_once(
  kls_egraph_refactor_shared *shared,
  UF_long col) {
  if (shared == NULL || shared->pipeline_done == NULL) {
    return 0;
  }
  const unsigned int generation = shared->pipeline_generation;
  if (atomic_load_explicit(&shared->pipeline_done[col],
                           memory_order_acquire) == generation) {
    return 0;
  }
  atomic_store_explicit(&shared->pipeline_done[col], generation,
                        memory_order_release);
  return 1;
}

static void kls_egraph_refactor_mark_done(
  kls_egraph_refactor_shared *shared,
  UF_long col) {
  (void)kls_egraph_refactor_mark_done_once(shared, col);
}

static void kls_row_refactor_mark_row_done(
  kls_egraph_refactor_worker *worker,
  UF_long row) {
  if (worker == NULL || worker->shared == NULL) {
    return;
  }
  kls_egraph_refactor_shared *shared = worker->shared;
  const int first_mark = kls_egraph_refactor_mark_done_once(shared, row);
  if (!first_mark || !shared->row_pipeline_row_dep_ready_queue) {
    return;
  }
  kls_solver *solver = shared->solver;
  if (solver == NULL ||
      solver->row_refactor_row_group == NULL ||
      solver->row_refactor_successor_ptr == NULL ||
      solver->row_refactor_successor_rows == NULL ||
      shared->row_pipeline_tail_groups == NULL ||
      shared->row_pipeline_remaining_preds == NULL ||
      row >= solver->n) {
    kls_egraph_refactor_record_invalid(shared);
    return;
  }
  const UF_long dep_group = solver->row_refactor_row_group[row];
  if (dep_group >= solver->row_refactor_group_count ||
      !shared->row_pipeline_tail_groups[dep_group]) {
    return;
  }
  const UF_long begin = solver->row_refactor_successor_ptr[row];
  const UF_long end = solver->row_refactor_successor_ptr[row + 1u];
  if (begin > end) {
    kls_egraph_refactor_record_invalid(shared);
    return;
  }
  for (UF_long pos = begin; pos < end; ++pos) {
    const UF_long successor_row = solver->row_refactor_successor_rows[pos];
    if (successor_row >= solver->n) {
      kls_egraph_refactor_record_invalid(shared);
      return;
    }
    const UF_long successor_group =
      solver->row_refactor_row_group[successor_row];
    if (successor_group >= solver->row_refactor_group_count) {
      kls_egraph_refactor_record_invalid(shared);
      return;
    }
    if (successor_group == dep_group ||
        !shared->row_pipeline_tail_groups[successor_group]) {
      continue;
    }
    const int successor_private =
      shared->row_pipeline_private_group_mask != NULL &&
      shared->row_pipeline_private_group_mask[successor_group];
    if (successor_private) {
      continue;
    }
    unsigned long current =
      atomic_load_explicit(
        &shared->row_pipeline_remaining_preds[successor_group],
        memory_order_acquire);
    for (;;) {
      if (current == ULONG_MAX) {
        break;
      }
      if (current == 0ul) {
        kls_egraph_refactor_record_invalid(shared);
        return;
      }
      const unsigned long next = current - 1ul;
      if (atomic_compare_exchange_weak_explicit(
            &shared->row_pipeline_remaining_preds[successor_group],
            &current, next, memory_order_acq_rel, memory_order_acquire)) {
        if (next == 0ul &&
            !kls_row_refactor_enqueue_ready_group(shared, successor_group)) {
          return;
        }
        break;
      }
    }
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
    if ((spin & 63u) == 1u &&
        kls_dense_help_try(shared, kls_dense_help_row_floor, 0)) {
      continue;
    }
    kls_cpu_relax();
    if ((spin++ & 1023u) == 0u &&
        kls_egraph_refactor_should_stop(shared)) {
      return 0;
    }
  }
  return 1;
}

static void kls_atomic_update_max_ulong(atomic_ulong *target,
                                        UF_long value) {
  if (target == NULL) {
    return;
  }
  const unsigned long capped = (unsigned long)value;
  unsigned long observed =
    atomic_load_explicit(target, memory_order_relaxed);
  while (observed < capped &&
         !atomic_compare_exchange_weak_explicit(
           target, &observed, capped, memory_order_relaxed,
           memory_order_relaxed)) {
  }
}

static void kls_egraph_record_btf_scalar_run_exec(
  kls_egraph_refactor_shared *shared,
  UF_long rows,
  UF_long entries) {
  if (shared == NULL || rows <= 1u) {
    return;
  }
  atomic_fetch_add_explicit(&shared->btf_scalar_run_exec_runs, 1ul,
                            memory_order_relaxed);
  atomic_fetch_add_explicit(&shared->btf_scalar_run_exec_rows,
                            (unsigned long)rows, memory_order_relaxed);
  atomic_fetch_add_explicit(&shared->btf_scalar_run_exec_entries,
                            (unsigned long)entries, memory_order_relaxed);
  kls_atomic_update_max_ulong(&shared->btf_scalar_run_exec_max_rows, rows);
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

static int kls_egraph_refactor_dependency_done_now(
  const kls_egraph_refactor_shared *shared,
  UF_long col) {
  if (shared == NULL || shared->pipeline_done == NULL ||
      shared->pipeline_generation == 0u ||
      shared->solver == NULL || col >= shared->solver->n) {
    return 0;
  }
  return atomic_load_explicit(&shared->pipeline_done[col],
                              memory_order_acquire) ==
         shared->pipeline_generation;
}

/* Cached-pointer variant of kls_snode_batch_consume for the EGraph column
 * kernels, which read L columns through the retained l_indices/l_values
 * caches and can walk 32-bit U index mirrors.  In pipeline mode a batch is
 * only taken when every producer in it is already published; otherwise the
 * scalar path performs its usual blocking wait.  Returns the number of
 * producers consumed (0 = no batch at position up). */
/* Keep the large panel temporaries out of the eligibility probe.  Most call
 * sites decline because the first dependency is not in a supernode run; a
 * single combined function made every one of those probes reserve the full
 * ~2.7 KiB AVX-512 update frame before inspecting run_end. */
#if defined(__GNUC__) && defined(__x86_64__) && !defined(__clang__)
/* GCC cannot instantiate the explicit AVX-512 FMA source below as its
   x86-64-v3 clone.  The default clone remains the fallback and native v4
   builds dispatch to the tuned clone. */
__attribute__((noinline, target_clones("default", "arch=x86-64-v4")))
#endif
static UF_long kls_snode_batch_consume_cached_apply(
  const kls_solver *solver,
  double *const *l_values,
  double *ux,
  UF_long up,
  double *restrict x,
  UF_long k1,
  UF_long j,
  UF_long t,
  const UF_long *tli,
  const int32_t *tli32,
  UF_long tlen) {
  double xs[KLS_SNODE_MAX_BATCH];
  const double *lx_arr[KLS_SNODE_MAX_BATCH];
  for (UF_long i = 0; i < t; ++i) {
    lx_arr[i] = l_values[k1 + j + i];
    if (lx_arr[i] == NULL) {
      return 0;
    }
    xs[i] = x[j + i];
  }
  for (UF_long i = 0; i < t; ++i) {
    x[j + i] = 0.0;
  }
  for (UF_long i = 0; i < t; ++i) {
    const double u = xs[i];
    ux[up + i] = u;
    const double *lxi = lx_arr[i];
    for (UF_long r0 = 0; r0 + i + 1u < t; ++r0) {
      xs[i + 1u + r0] -= lxi[r0] * u;
    }
  }
  const double *tlx = lx_arr[t - 1u];
  UF_long p0 = 0;
  while (p0 < tlen) {
    const UF_long remaining = tlen - p0;
#if defined(__AVX512F__)
    const UF_long pc =
                       remaining >= 144u && solver != NULL &&
                           solver->snode_wide_tail
                         ? 144u
                       : remaining >= 128u && solver != NULL &&
                           solver->snode_wide_tail
                         ? 128u
                       : remaining >= 64u ? 64u
                       : remaining > KLS_SNODE_TAIL_CHUNK &&
                           solver != NULL &&
                           solver->snode_wide_tail
                         ? remaining
                       : remaining >= KLS_SNODE_TAIL_CHUNK
                         ? KLS_SNODE_TAIL_CHUNK : remaining;
#else
    const UF_long pc = remaining < KLS_SNODE_TAIL_CHUNK
                         ? remaining : KLS_SNODE_TAIL_CHUNK;
#endif
    double acc[144];
#if defined(__AVX512F__)
    if (pc == 144u) {
      __m512d u = _mm512_set1_pd(xs[t - 1u]);
#define KLS_SNODE_INIT_ACC(I) \
      __m512d a##I = _mm512_mul_pd( \
        _mm512_loadu_pd(tlx + p0 + 8u * (I)), u)
      KLS_SNODE_INIT_ACC(0);
      KLS_SNODE_INIT_ACC(1);
      KLS_SNODE_INIT_ACC(2);
      KLS_SNODE_INIT_ACC(3);
      KLS_SNODE_INIT_ACC(4);
      KLS_SNODE_INIT_ACC(5);
      KLS_SNODE_INIT_ACC(6);
      KLS_SNODE_INIT_ACC(7);
      KLS_SNODE_INIT_ACC(8);
      KLS_SNODE_INIT_ACC(9);
      KLS_SNODE_INIT_ACC(10);
      KLS_SNODE_INIT_ACC(11);
      KLS_SNODE_INIT_ACC(12);
      KLS_SNODE_INIT_ACC(13);
      KLS_SNODE_INIT_ACC(14);
      KLS_SNODE_INIT_ACC(15);
      KLS_SNODE_INIT_ACC(16);
      KLS_SNODE_INIT_ACC(17);
#undef KLS_SNODE_INIT_ACC
      for (UF_long i = 0; i + 1u < t; ++i) {
        const double *src = lx_arr[i] + (t - 1u - i) + p0;
        u = _mm512_set1_pd(xs[i]);
#define KLS_SNODE_UPDATE_ACC(I) \
        a##I = _mm512_fmadd_pd( \
          _mm512_loadu_pd(src + 8u * (I)), u, a##I)
        KLS_SNODE_UPDATE_ACC(0);
        KLS_SNODE_UPDATE_ACC(1);
        KLS_SNODE_UPDATE_ACC(2);
        KLS_SNODE_UPDATE_ACC(3);
        KLS_SNODE_UPDATE_ACC(4);
        KLS_SNODE_UPDATE_ACC(5);
        KLS_SNODE_UPDATE_ACC(6);
        KLS_SNODE_UPDATE_ACC(7);
        KLS_SNODE_UPDATE_ACC(8);
        KLS_SNODE_UPDATE_ACC(9);
        KLS_SNODE_UPDATE_ACC(10);
        KLS_SNODE_UPDATE_ACC(11);
        KLS_SNODE_UPDATE_ACC(12);
        KLS_SNODE_UPDATE_ACC(13);
        KLS_SNODE_UPDATE_ACC(14);
        KLS_SNODE_UPDATE_ACC(15);
        KLS_SNODE_UPDATE_ACC(16);
        KLS_SNODE_UPDATE_ACC(17);
#undef KLS_SNODE_UPDATE_ACC
      }
#define KLS_SNODE_STORE_ACC(I) \
      _mm512_storeu_pd(acc + 8u * (I), a##I)
      KLS_SNODE_STORE_ACC(0);
      KLS_SNODE_STORE_ACC(1);
      KLS_SNODE_STORE_ACC(2);
      KLS_SNODE_STORE_ACC(3);
      KLS_SNODE_STORE_ACC(4);
      KLS_SNODE_STORE_ACC(5);
      KLS_SNODE_STORE_ACC(6);
      KLS_SNODE_STORE_ACC(7);
      KLS_SNODE_STORE_ACC(8);
      KLS_SNODE_STORE_ACC(9);
      KLS_SNODE_STORE_ACC(10);
      KLS_SNODE_STORE_ACC(11);
      KLS_SNODE_STORE_ACC(12);
      KLS_SNODE_STORE_ACC(13);
      KLS_SNODE_STORE_ACC(14);
      KLS_SNODE_STORE_ACC(15);
      KLS_SNODE_STORE_ACC(16);
      KLS_SNODE_STORE_ACC(17);
#undef KLS_SNODE_STORE_ACC
    } else if (pc == 128u) {
      /* A 185-entry average shared tail otherwise scans every producer four
         times (64+64+32+remainder).  Sixteen accumulators fit in the
         32-register AVX-512 file, so a 128-entry block retains the arithmetic
         order while removing one producer-pointer/broadcast scan. */
#define KLS_SNODE_INIT_ACC(I) \
      __m512d a##I = _mm512_mul_pd( \
        _mm512_loadu_pd(tlx + p0 + 8u * (I)), u)
      __m512d u = _mm512_set1_pd(xs[t - 1u]);
      KLS_SNODE_INIT_ACC(0);
      KLS_SNODE_INIT_ACC(1);
      KLS_SNODE_INIT_ACC(2);
      KLS_SNODE_INIT_ACC(3);
      KLS_SNODE_INIT_ACC(4);
      KLS_SNODE_INIT_ACC(5);
      KLS_SNODE_INIT_ACC(6);
      KLS_SNODE_INIT_ACC(7);
      KLS_SNODE_INIT_ACC(8);
      KLS_SNODE_INIT_ACC(9);
      KLS_SNODE_INIT_ACC(10);
      KLS_SNODE_INIT_ACC(11);
      KLS_SNODE_INIT_ACC(12);
      KLS_SNODE_INIT_ACC(13);
      KLS_SNODE_INIT_ACC(14);
      KLS_SNODE_INIT_ACC(15);
#undef KLS_SNODE_INIT_ACC
      for (UF_long i = 0; i + 1u < t; ++i) {
        const double *src = lx_arr[i] + (t - 1u - i) + p0;
        u = _mm512_set1_pd(xs[i]);
#define KLS_SNODE_UPDATE_ACC(I) \
        a##I = _mm512_fmadd_pd( \
          _mm512_loadu_pd(src + 8u * (I)), u, a##I)
        KLS_SNODE_UPDATE_ACC(0);
        KLS_SNODE_UPDATE_ACC(1);
        KLS_SNODE_UPDATE_ACC(2);
        KLS_SNODE_UPDATE_ACC(3);
        KLS_SNODE_UPDATE_ACC(4);
        KLS_SNODE_UPDATE_ACC(5);
        KLS_SNODE_UPDATE_ACC(6);
        KLS_SNODE_UPDATE_ACC(7);
        KLS_SNODE_UPDATE_ACC(8);
        KLS_SNODE_UPDATE_ACC(9);
        KLS_SNODE_UPDATE_ACC(10);
        KLS_SNODE_UPDATE_ACC(11);
        KLS_SNODE_UPDATE_ACC(12);
        KLS_SNODE_UPDATE_ACC(13);
        KLS_SNODE_UPDATE_ACC(14);
        KLS_SNODE_UPDATE_ACC(15);
#undef KLS_SNODE_UPDATE_ACC
      }
#define KLS_SNODE_STORE_ACC(I) \
      _mm512_storeu_pd(acc + 8u * (I), a##I)
      KLS_SNODE_STORE_ACC(0);
      KLS_SNODE_STORE_ACC(1);
      KLS_SNODE_STORE_ACC(2);
      KLS_SNODE_STORE_ACC(3);
      KLS_SNODE_STORE_ACC(4);
      KLS_SNODE_STORE_ACC(5);
      KLS_SNODE_STORE_ACC(6);
      KLS_SNODE_STORE_ACC(7);
      KLS_SNODE_STORE_ACC(8);
      KLS_SNODE_STORE_ACC(9);
      KLS_SNODE_STORE_ACC(10);
      KLS_SNODE_STORE_ACC(11);
      KLS_SNODE_STORE_ACC(12);
      KLS_SNODE_STORE_ACC(13);
      KLS_SNODE_STORE_ACC(14);
      KLS_SNODE_STORE_ACC(15);
#undef KLS_SNODE_STORE_ACC
    } else if (pc == 64u) {
      /* Eight accumulators still fit comfortably in the v4 register file.
         Consuming two old chunks at once halves producer-pointer loads and
         broadcasts on long shared tails. */
      __m512d u = _mm512_set1_pd(xs[t - 1u]);
      __m512d a0 = _mm512_mul_pd(_mm512_loadu_pd(tlx + p0), u);
      __m512d a1 = _mm512_mul_pd(_mm512_loadu_pd(tlx + p0 + 8u), u);
      __m512d a2 = _mm512_mul_pd(_mm512_loadu_pd(tlx + p0 + 16u), u);
      __m512d a3 = _mm512_mul_pd(_mm512_loadu_pd(tlx + p0 + 24u), u);
      __m512d a4 = _mm512_mul_pd(_mm512_loadu_pd(tlx + p0 + 32u), u);
      __m512d a5 = _mm512_mul_pd(_mm512_loadu_pd(tlx + p0 + 40u), u);
      __m512d a6 = _mm512_mul_pd(_mm512_loadu_pd(tlx + p0 + 48u), u);
      __m512d a7 = _mm512_mul_pd(_mm512_loadu_pd(tlx + p0 + 56u), u);
      for (UF_long i = 0; i + 1u < t; ++i) {
        const double *src = lx_arr[i] + (t - 1u - i) + p0;
        u = _mm512_set1_pd(xs[i]);
        a0 = _mm512_fmadd_pd(_mm512_loadu_pd(src), u, a0);
        a1 = _mm512_fmadd_pd(_mm512_loadu_pd(src + 8u), u, a1);
        a2 = _mm512_fmadd_pd(_mm512_loadu_pd(src + 16u), u, a2);
        a3 = _mm512_fmadd_pd(_mm512_loadu_pd(src + 24u), u, a3);
        a4 = _mm512_fmadd_pd(_mm512_loadu_pd(src + 32u), u, a4);
        a5 = _mm512_fmadd_pd(_mm512_loadu_pd(src + 40u), u, a5);
        a6 = _mm512_fmadd_pd(_mm512_loadu_pd(src + 48u), u, a6);
        a7 = _mm512_fmadd_pd(_mm512_loadu_pd(src + 56u), u, a7);
      }
      _mm512_storeu_pd(acc, a0);
      _mm512_storeu_pd(acc + 8u, a1);
      _mm512_storeu_pd(acc + 16u, a2);
      _mm512_storeu_pd(acc + 24u, a3);
      _mm512_storeu_pd(acc + 32u, a4);
      _mm512_storeu_pd(acc + 40u, a5);
      _mm512_storeu_pd(acc + 48u, a6);
      _mm512_storeu_pd(acc + 56u, a7);
    } else if (pc >= 24u && pc < KLS_SNODE_TAIL_CHUNK) {
      const __mmask8 tail_mask =
        (__mmask8)((1u << (unsigned)(pc - 24u)) - 1u);
      __m512d u = _mm512_set1_pd(xs[t - 1u]);
      __m512d a0 = _mm512_mul_pd(_mm512_loadu_pd(tlx + p0), u);
      __m512d a1 = _mm512_mul_pd(_mm512_loadu_pd(tlx + p0 + 8u), u);
      __m512d a2 = _mm512_mul_pd(_mm512_loadu_pd(tlx + p0 + 16u), u);
      __m512d a3 = _mm512_mul_pd(
        _mm512_maskz_loadu_pd(tail_mask, tlx + p0 + 24u), u);
      for (UF_long i = 0; i + 1u < t; ++i) {
        const double *src = lx_arr[i] + (t - 1u - i) + p0;
        u = _mm512_set1_pd(xs[i]);
        a0 = _mm512_fmadd_pd(_mm512_loadu_pd(src), u, a0);
        a1 = _mm512_fmadd_pd(_mm512_loadu_pd(src + 8u), u, a1);
        a2 = _mm512_fmadd_pd(_mm512_loadu_pd(src + 16u), u, a2);
        a3 = _mm512_fmadd_pd(
          _mm512_maskz_loadu_pd(tail_mask, src + 24u), u, a3);
      }
      _mm512_storeu_pd(acc, a0);
      _mm512_storeu_pd(acc + 8u, a1);
      _mm512_storeu_pd(acc + 16u, a2);
      _mm512_mask_storeu_pd(acc + 24u, tail_mask, a3);
    } else if (pc >= 16u && pc < 24u) {
      const __mmask8 tail_mask =
        (__mmask8)((1u << (unsigned)(pc - 16u)) - 1u);
      __m512d u = _mm512_set1_pd(xs[t - 1u]);
      __m512d a0 = _mm512_mul_pd(_mm512_loadu_pd(tlx + p0), u);
      __m512d a1 = _mm512_mul_pd(_mm512_loadu_pd(tlx + p0 + 8u), u);
      __m512d a2 = _mm512_mul_pd(
        _mm512_maskz_loadu_pd(tail_mask, tlx + p0 + 16u), u);
      for (UF_long i = 0; i + 1u < t; ++i) {
        const double *src = lx_arr[i] + (t - 1u - i) + p0;
        u = _mm512_set1_pd(xs[i]);
        a0 = _mm512_fmadd_pd(_mm512_loadu_pd(src), u, a0);
        a1 = _mm512_fmadd_pd(_mm512_loadu_pd(src + 8u), u, a1);
        a2 = _mm512_fmadd_pd(
          _mm512_maskz_loadu_pd(tail_mask, src + 16u), u, a2);
      }
      _mm512_storeu_pd(acc, a0);
      _mm512_storeu_pd(acc + 8u, a1);
      _mm512_mask_storeu_pd(acc + 16u, tail_mask, a2);
    } else if (pc >= 8u && pc < 16u) {
      const __mmask8 tail_mask =
        (__mmask8)((1u << (unsigned)(pc - 8u)) - 1u);
      __m512d u = _mm512_set1_pd(xs[t - 1u]);
      __m512d a0 = _mm512_mul_pd(_mm512_loadu_pd(tlx + p0), u);
      __m512d a1 = _mm512_mul_pd(
        _mm512_maskz_loadu_pd(tail_mask, tlx + p0 + 8u), u);
      for (UF_long i = 0; i + 1u < t; ++i) {
        const double *src = lx_arr[i] + (t - 1u - i) + p0;
        u = _mm512_set1_pd(xs[i]);
        a0 = _mm512_fmadd_pd(_mm512_loadu_pd(src), u, a0);
        a1 = _mm512_fmadd_pd(
          _mm512_maskz_loadu_pd(tail_mask, src + 8u), u, a1);
      }
      _mm512_storeu_pd(acc, a0);
      _mm512_mask_storeu_pd(acc + 8u, tail_mask, a1);
    } else if (pc > KLS_SNODE_TAIL_CHUNK && pc <= 40u) {
      /* A 144-entry main tile commonly leaves only 33--40 entries.  Keep
         that remainder in five zmm accumulators and mask the final vector,
         avoiding a second producer/broadcast pass for the last 1--8 rows. */
      const __mmask8 tail_mask =
        (__mmask8)((1u << (unsigned)(pc - 32u)) - 1u);
      __m512d u = _mm512_set1_pd(xs[t - 1u]);
      __m512d a0 = _mm512_mul_pd(_mm512_loadu_pd(tlx + p0), u);
      __m512d a1 = _mm512_mul_pd(_mm512_loadu_pd(tlx + p0 + 8u), u);
      __m512d a2 = _mm512_mul_pd(_mm512_loadu_pd(tlx + p0 + 16u), u);
      __m512d a3 = _mm512_mul_pd(_mm512_loadu_pd(tlx + p0 + 24u), u);
      __m512d a4 = _mm512_mul_pd(
        _mm512_maskz_loadu_pd(tail_mask, tlx + p0 + 32u), u);
      for (UF_long i = 0; i + 1u < t; ++i) {
        const double *src = lx_arr[i] + (t - 1u - i) + p0;
        u = _mm512_set1_pd(xs[i]);
        a0 = _mm512_fmadd_pd(_mm512_loadu_pd(src), u, a0);
        a1 = _mm512_fmadd_pd(_mm512_loadu_pd(src + 8u), u, a1);
        a2 = _mm512_fmadd_pd(_mm512_loadu_pd(src + 16u), u, a2);
        a3 = _mm512_fmadd_pd(_mm512_loadu_pd(src + 24u), u, a3);
        a4 = _mm512_fmadd_pd(
          _mm512_maskz_loadu_pd(tail_mask, src + 32u), u, a4);
      }
      _mm512_storeu_pd(acc, a0);
      _mm512_storeu_pd(acc + 8u, a1);
      _mm512_storeu_pd(acc + 16u, a2);
      _mm512_storeu_pd(acc + 24u, a3);
      _mm512_mask_storeu_pd(acc + 32u, tail_mask, a4);
    } else if (pc == KLS_SNODE_TAIL_CHUNK) {
      /* GCC's generic dynamic-trip loop reloads and stores acc[] for every
         producer, even in the x86-64-v4 clone.  A full tail chunk is exactly
         four vectors, so keep those accumulators live across the producer
         dimension and materialize acc[] only once for the sparse scatter. */
      __m512d u = _mm512_set1_pd(xs[t - 1u]);
      __m512d a0 = _mm512_mul_pd(_mm512_loadu_pd(tlx + p0), u);
      __m512d a1 = _mm512_mul_pd(_mm512_loadu_pd(tlx + p0 + 8u), u);
      __m512d a2 = _mm512_mul_pd(_mm512_loadu_pd(tlx + p0 + 16u), u);
      __m512d a3 = _mm512_mul_pd(_mm512_loadu_pd(tlx + p0 + 24u), u);
      for (UF_long i = 0; i + 1u < t; ++i) {
        const double *src = lx_arr[i] + (t - 1u - i) + p0;
        u = _mm512_set1_pd(xs[i]);
        a0 = _mm512_fmadd_pd(_mm512_loadu_pd(src), u, a0);
        a1 = _mm512_fmadd_pd(_mm512_loadu_pd(src + 8u), u, a1);
        a2 = _mm512_fmadd_pd(_mm512_loadu_pd(src + 16u), u, a2);
        a3 = _mm512_fmadd_pd(_mm512_loadu_pd(src + 24u), u, a3);
      }
      _mm512_storeu_pd(acc, a0);
      _mm512_storeu_pd(acc + 8u, a1);
      _mm512_storeu_pd(acc + 16u, a2);
      _mm512_storeu_pd(acc + 24u, a3);
    } else
#endif
    {
    {
      const double *src = tlx + p0;
      const double u = xs[t - 1u];
      for (UF_long p = 0; p < pc; ++p) {
        acc[p] = src[p] * u;
      }
    }
    for (UF_long i = 0; i + 1u < t; ++i) {
      const double *src = lx_arr[i] + (t - 1u - i) + p0;
      const double u = xs[i];
      for (UF_long p = 0; p < pc; ++p) {
        acc[p] += src[p] * u;
      }
    }
    }
    if (tli32 != NULL) {
      /* The persistent 32-bit row mirror halves index traffic in the
         batch's final indexed RMW and shares the same AVX-512 crossover
         machinery as ordinary L scatters.  This scatter is the remaining
         non-contiguous phase after the dense tail panel has been fused. */
      kls_scatter_subtract_refactor_i32(
        solver, x, tli32 + p0, acc, pc, 1.0);
    } else {
      for (UF_long p = 0; p < pc; ++p) {
        x[tli[p0 + p]] -= acc[p];
      }
    }
    p0 += pc;
  }
  return t;
}

static KLS_ALWAYS_INLINE UF_long kls_snode_batch_consume_cached(
  UF_long *const *l_indices,
  double *const *l_values,
  const UF_long *llen_local,
  const UF_long *ui,
  const int32_t *ui32,
  double *ux,
  UF_long ucol_len,
  UF_long up,
  double *restrict x,
  UF_long k1,
  UF_long producer_limit,
  const UF_long *snode_run_end,
  const kls_egraph_refactor_shared *shared,
  int wait_for_dependencies) {
  const UF_long j = ui32 != NULL ? (UF_long)ui32[up] : ui[up];
  const UF_long run_end = snode_run_end[k1 + j];
  if (run_end <= k1 + j + 1u) {
    return 0;
  }
  UF_long tmax = run_end - (k1 + j);
  if (tmax > KLS_SNODE_MAX_BATCH) {
    tmax = KLS_SNODE_MAX_BATCH;
  }
  if (tmax > ucol_len - up) {
    tmax = ucol_len - up;
  }
  UF_long t = 1;
  if (ui32 != NULL) {
    while (t < tmax && (UF_long)ui32[up + t] == j + t) {
      t++;
    }
  } else {
    while (t < tmax && ui[up + t] == j + t) {
      t++;
    }
  }
  if (t < KLS_SNODE_MIN_BATCH || j + t > producer_limit) {
    return 0;
  }
  if (wait_for_dependencies) {
    for (UF_long i = 0; i < t; ++i) {
      if (!kls_egraph_refactor_dependency_done_now(shared, k1 + j + i)) {
        return 0;
      }
    }
  }
  const UF_long *tli = l_indices[k1 + j + t - 1u];
  const int32_t *tli32 =
    shared != NULL && shared->solver != NULL &&
        shared->solver->refactor_l_indices32 != NULL
      ? shared->solver->refactor_l_indices32[k1 + j + t - 1u]
      : NULL;
  const UF_long tlen = llen_local[j + t - 1u];
  if (tli == NULL || t * tlen < KLS_SNODE_MIN_BATCH_WORK) {
    /* Short shared tails lose to the 32-bit-index scalar path; only pay the
       panel staging when the batched update amortizes it. */
    return 0;
  }
  return kls_snode_batch_consume_cached_apply(
    shared != NULL ? shared->solver : NULL,
    l_values, ux, up, x, k1, j, t, tli, tli32, tlen);
}

/* Fused two-consumer batch consume: one pass over the producer run's tail
   panel feeds two consumer SPAs (the panel stream is the per-core wall;
   sharing it halves the dominant refactor traffic).  Both consumers must
   hold the run contiguously at their current U positions; the run's shared
   tail pattern (last member) is producer-side structure and thus consumer
   independent.  Per-consumer arithmetic order matches the single kernel
   exactly, so paired columns are bit-identical to single dispatch.
   Returns producers consumed from BOTH positions, or 0 to fall back. */
static UF_long kls_snode_batch_consume_cached_pair(
  UF_long *const *l_indices,
  double *const *l_values,
  const UF_long *llen_local,
  const UF_long *uia,
  const int32_t *uia32,
  double *uxa,
  UF_long ulen_a,
  UF_long upa,
  double *restrict xa,
  const UF_long *uib,
  const int32_t *uib32,
  double *uxb,
  UF_long ulen_b,
  UF_long upb,
  double *restrict xb,
  UF_long k1,
  UF_long producer_limit,
  const UF_long *snode_run_end,
  const kls_egraph_refactor_shared *shared) {
  const UF_long j = uia32 != NULL ? (UF_long)uia32[upa] : uia[upa];
  const UF_long run_end = snode_run_end[k1 + j];
  if (run_end <= k1 + j + 1u) {
    return 0;
  }
  UF_long tmax = run_end - (k1 + j);
  if (tmax > KLS_SNODE_MAX_BATCH) {
    tmax = KLS_SNODE_MAX_BATCH;
  }
  if (tmax > ulen_a - upa) {
    tmax = ulen_a - upa;
  }
  if (tmax > ulen_b - upb) {
    tmax = ulen_b - upb;
  }
  UF_long t = 1;
  while (t < tmax &&
         (uia32 != NULL ? (UF_long)uia32[upa + t] : uia[upa + t]) == j + t &&
         (uib32 != NULL ? (UF_long)uib32[upb + t] : uib[upb + t]) == j + t) {
    t++;
  }
  if (t < KLS_SNODE_MIN_BATCH || j + t > producer_limit) {
    return 0;
  }
  {
    for (UF_long i = 0; i < t; ++i) {
      if (!kls_egraph_refactor_dependency_done_now(shared, k1 + j + i)) {
        return 0;
      }
    }
  }
  const UF_long *tli = l_indices[k1 + j + t - 1u];
  const UF_long tlen = llen_local[j + t - 1u];
  if (tli == NULL || 2u * t * tlen < KLS_SNODE_MIN_BATCH_WORK) {
    return 0;
  }
  double xsa[KLS_SNODE_MAX_BATCH];
  double xsb[KLS_SNODE_MAX_BATCH];
  const double *lx_arr[KLS_SNODE_MAX_BATCH];
  for (UF_long i = 0; i < t; ++i) {
    lx_arr[i] = l_values[k1 + j + i];
    if (lx_arr[i] == NULL) {
      return 0;
    }
    xsa[i] = xa[j + i];
    xsb[i] = xb[j + i];
  }
  for (UF_long i = 0; i < t; ++i) {
    xa[j + i] = 0.0;
    xb[j + i] = 0.0;
  }
  for (UF_long i = 0; i < t; ++i) {
    const double ua = xsa[i];
    const double ub = xsb[i];
    uxa[upa + i] = ua;
    uxb[upb + i] = ub;
    const double *lxi = lx_arr[i];
    for (UF_long r0 = 0; r0 + i + 1u < t; ++r0) {
      xsa[i + 1u + r0] -= lxi[r0] * ua;
      xsb[i + 1u + r0] -= lxi[r0] * ub;
    }
  }
  const double *tlx = lx_arr[t - 1u];
  for (UF_long p0 = 0; p0 < tlen; p0 += KLS_SNODE_TAIL_CHUNK) {
    const UF_long pc = tlen - p0 < KLS_SNODE_TAIL_CHUNK
                         ? tlen - p0
                         : KLS_SNODE_TAIL_CHUNK;
    double acca[KLS_SNODE_TAIL_CHUNK];
    double accb[KLS_SNODE_TAIL_CHUNK];
    {
      const double *src = tlx + p0;
      const double ua = xsa[t - 1u];
      const double ub = xsb[t - 1u];
      for (UF_long p = 0; p < pc; ++p) {
        acca[p] = src[p] * ua;
        accb[p] = src[p] * ub;
      }
    }
    for (UF_long i = 0; i + 1u < t; ++i) {
      const double *src = lx_arr[i] + (t - 1u - i) + p0;
      const double ua = xsa[i];
      const double ub = xsb[i];
      for (UF_long p = 0; p < pc; ++p) {
        acca[p] += src[p] * ua;
        accb[p] += src[p] * ub;
      }
    }
    for (UF_long p = 0; p < pc; ++p) {
      const UF_long row = tli[p0 + p];
      xa[row] -= acca[p];
      xb[row] -= accb[p];
    }
  }
  return t;
}

#define KLS_SNODE_MAX_FUSE 4

/* Fused k-consumer batch consume (k <= KLS_SNODE_MAX_FUSE): one pass
   over the producer run's tail panel feeds k consumer SPAs, dividing
   the dominant panel-stream traffic by k. All consumers' cursors must
   point at the run start (caller-checked) and hold it contiguously;
   per-consumer arithmetic order matches the single kernel exactly, so
   fused columns stay bit-identical to single dispatch. Returns
   producers consumed from EVERY cursor, or 0 to fall back. */
static UF_long kls_snode_batch_consume_cached_multi(
  UF_long *const *l_indices,
  double *const *l_values,
  const UF_long *llen_local,
  int nc,
  const UF_long *const *ui,
  const int32_t *const *ui32,
  double *const *ux,
  const UF_long *ulen,
  const UF_long *up,
  double *const *xv,
  UF_long k1,
  UF_long producer_limit,
  const UF_long *snode_run_end,
  const kls_egraph_refactor_shared *shared) {
  const UF_long j =
    ui32[0] != NULL ? (UF_long)ui32[0][up[0]] : ui[0][up[0]];
  const UF_long run_end = snode_run_end[k1 + j];
  if (run_end <= k1 + j + 1u) {
    return 0;
  }
  UF_long tmax = run_end - (k1 + j);
  if (tmax > KLS_SNODE_MAX_BATCH) {
    tmax = KLS_SNODE_MAX_BATCH;
  }
  for (int c = 0; c < nc; ++c) {
    if (tmax > ulen[c] - up[c]) {
      tmax = ulen[c] - up[c];
    }
  }
  UF_long t = 1;
  for (; t < tmax; ++t) {
    int aligned = 1;
    for (int c = 0; c < nc; ++c) {
      const UF_long jc = ui32[c] != NULL
        ? (UF_long)ui32[c][up[c] + t] : ui[c][up[c] + t];
      if (jc != j + t) {
        aligned = 0;
        break;
      }
    }
    if (!aligned) {
      break;
    }
  }
  if (t < KLS_SNODE_MIN_BATCH || j + t > producer_limit) {
    return 0;
  }
  {
    for (UF_long i = 0; i < t; ++i) {
      if (!kls_egraph_refactor_dependency_done_now(shared, k1 + j + i)) {
        return 0;
      }
    }
  }
  const UF_long *tli = l_indices[k1 + j + t - 1u];
  const UF_long tlen = llen_local[j + t - 1u];
  if (tli == NULL || (UF_long)nc * t * tlen < KLS_SNODE_MIN_BATCH_WORK) {
    return 0;
  }
  double xs[KLS_SNODE_MAX_FUSE][KLS_SNODE_MAX_BATCH];
  const double *lx_arr[KLS_SNODE_MAX_BATCH];
  for (UF_long i = 0; i < t; ++i) {
    lx_arr[i] = l_values[k1 + j + i];
    if (lx_arr[i] == NULL) {
      return 0;
    }
    for (int c = 0; c < nc; ++c) {
      xs[c][i] = xv[c][j + i];
    }
  }
  for (UF_long i = 0; i < t; ++i) {
    for (int c = 0; c < nc; ++c) {
      xv[c][j + i] = 0.0;
    }
  }
  for (UF_long i = 0; i < t; ++i) {
    const double *lxi = lx_arr[i];
    for (int c = 0; c < nc; ++c) {
      const double u = xs[c][i];
      ux[c][up[c] + i] = u;
      for (UF_long r0 = 0; r0 + i + 1u < t; ++r0) {
        xs[c][i + 1u + r0] -= lxi[r0] * u;
      }
    }
  }
  const double *tlx = lx_arr[t - 1u];
  for (UF_long p0 = 0; p0 < tlen; p0 += KLS_SNODE_TAIL_CHUNK) {
    const UF_long pc = tlen - p0 < KLS_SNODE_TAIL_CHUNK
                         ? tlen - p0
                         : KLS_SNODE_TAIL_CHUNK;
    double acc[KLS_SNODE_MAX_FUSE][KLS_SNODE_TAIL_CHUNK];
    {
      const double *src = tlx + p0;
      for (int c = 0; c < nc; ++c) {
        const double u = xs[c][t - 1u];
        for (UF_long p = 0; p < pc; ++p) {
          acc[c][p] = src[p] * u;
        }
      }
    }
    for (UF_long i = 0; i + 1u < t; ++i) {
      const double *src = lx_arr[i] + (t - 1u - i) + p0;
      for (int c = 0; c < nc; ++c) {
        const double u = xs[c][i];
        for (UF_long p = 0; p < pc; ++p) {
          acc[c][p] += src[p] * u;
        }
      }
    }
    for (UF_long p = 0; p < pc; ++p) {
      const UF_long row = tli[p0 + p];
      for (int c = 0; c < nc; ++c) {
        xv[c][row] -= acc[c][p];
      }
    }
  }
  return t;
}

static void kls_egraph_refactor_apply_btf_scalar_dep(
  kls_egraph_refactor_shared *shared,
  UF_long k1,
  UF_long dep_local,
  UF_long dep_pos,
  double *ux,
  UF_long **l_indices,
  int32_t **l_indices32,
  double **l_values,
  const UF_long *llen,
  double *x) {
  const UF_long dep_global = k1 + dep_local;
  const double ujk = x[dep_local];
  x[dep_local] = 0.0;
  ux[dep_pos] = ujk;

  if (ujk != 0.0) {
    UF_long *li = l_indices[dep_global];
    const int32_t *li32 =
      l_indices32 != NULL ? l_indices32[dep_global] : NULL;
    double *lx = l_values[dep_global];
    UF_long lcol_len = llen[dep_local];
    if (li32 != NULL || lcol_len == 0u) {
      kls_scatter_subtract_refactor_i32(
        shared->solver, x, li32, lx, lcol_len, ujk);
    } else {
      kls_scatter_subtract(x, li, lx, lcol_len, ujk);
    }
  }
}

static UF_long kls_egraph_btf_scalar_producer_run_len(
  UF_long start_pos,
  UF_long ucol_len,
  const UF_long *ui,
  const int32_t *ui32,
  UF_long local_k) {
  if (start_pos >= ucol_len || (ui == NULL && ui32 == NULL)) {
    return 0u;
  }
  const UF_long dep0 =
    ui32 != NULL ? (UF_long)ui32[start_pos] : ui[start_pos];
  if (dep0 >= local_k) {
    return 0u;
  }
  UF_long rows = 1u;
  while (start_pos + rows < ucol_len &&
         rows < KLS_BTF_SCALAR_RUN_EXEC_MAX_ROWS) {
    if (dep0 > UF_long_max - rows) {
      break;
    }
    const UF_long dep =
      ui32 != NULL ? (UF_long)ui32[start_pos + rows] : ui[start_pos + rows];
    if (dep != dep0 + rows || dep >= local_k) {
      break;
    }
    rows++;
  }
  return rows;
}

static int kls_egraph_refactor_try_btf_scalar_producer_run(
  kls_egraph_refactor_worker *worker,
  UF_long k1,
  UF_long local_k,
  UF_long *up,
  UF_long ucol_len,
  const UF_long *ui,
  const int32_t *ui32,
  double *ux,
  UF_long **l_indices,
  int32_t **l_indices32,
  double **l_values,
  const UF_long *llen,
  double *x,
  int wait_for_dependencies) {
  if (worker == NULL || worker->shared == NULL || up == NULL ||
      !worker->shared->btf_scalar_run_exec || *up >= ucol_len ||
      (ui == NULL && ui32 == NULL) || ux == NULL || l_indices == NULL ||
      l_values == NULL || llen == NULL || x == NULL) {
    return 0;
  }
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  const int use_sorted_mirror =
    solver != NULL && solver->refactor_l_sorted_enabled &&
    solver->refactor_l_sorted_indices32 != NULL &&
    solver->refactor_l_sorted_values != NULL;
  if (solver == NULL ||
      (!use_sorted_mirror && !solver->refactor_l_indices_sorted)) {
    return 0;
  }
  const UF_long start_pos = *up;
  const UF_long rows =
    kls_egraph_btf_scalar_producer_run_len(
      start_pos, ucol_len, ui, ui32, local_k);
  if (rows < KLS_BTF_SCALAR_RUN_EXEC_MIN_ROWS) {
    return 0;
  }
  const UF_long dep0 =
    ui32 != NULL ? (UF_long)ui32[start_pos] : ui[start_pos];
  if (dep0 > UF_long_max - rows) {
    return 0;
  }
  const UF_long run_end = dep0 + rows;

  double *workspace = kls_egraph_worker_supernode_workspace(worker, rows);
  if (workspace == NULL) {
    return 0;
  }

  UF_long entries = 0u;
  for (UF_long local = 0u; local < rows; ++local) {
    const UF_long dep = dep0 + local;
    if (wait_for_dependencies &&
        !kls_egraph_refactor_wait_done(shared, k1 + dep)) {
      return -1;
    }
    const UF_long lcol_len = llen[dep];
    if (lcol_len > UF_long_max - entries) {
      entries = UF_long_max;
    } else {
      entries += lcol_len;
    }
    const UF_long dep_global = k1 + dep;
    const int32_t *li32 =
      use_sorted_mirror ? solver->refactor_l_sorted_indices32[dep_global]
                        : (l_indices32 != NULL ? l_indices32[dep_global]
                                                : NULL);
    const double *lx_check =
      use_sorted_mirror ? solver->refactor_l_sorted_values[dep_global]
                        : l_values[dep_global];
    if (lcol_len > 0u &&
        (lx_check == NULL || (use_sorted_mirror && li32 == NULL))) {
      return 0;
    }
    if (li32 == NULL && lcol_len != 0u) {
      const UF_long *li = l_indices[dep_global];
      if (li == NULL) {
        return 0;
      }
    }
  }

  for (UF_long local = 0u; local < rows; ++local) {
    const UF_long dep = dep0 + local;
    workspace[local] = x[dep];
    x[dep] = 0.0;
  }

  for (UF_long local = 0u; local < rows; ++local) {
    const UF_long dep = dep0 + local;
    const UF_long dep_global = k1 + dep;
    const double ujk = workspace[local];
    ux[start_pos + local] = ujk;

    if (ujk == 0.0) {
      continue;
    }
    const UF_long lcol_len = llen[dep];
    double *lx = use_sorted_mirror
                   ? solver->refactor_l_sorted_values[dep_global]
                   : l_values[dep_global];
    const int32_t *li32 =
      use_sorted_mirror ? solver->refactor_l_sorted_indices32[dep_global]
                        : (l_indices32 != NULL ? l_indices32[dep_global]
                                                : NULL);
    if (li32 != NULL || lcol_len == 0u) {
      UF_long p = 0u;
      for (; p < lcol_len; ++p) {
        const UF_long row = (UF_long)li32[p];
        if (row >= run_end) {
          break;
        }
        workspace[row - dep0] -= ujk * lx[p];
      }
      if (p < lcol_len) {
        kls_scatter_subtract_refactor_i32(
          solver, x, li32 + p, lx + p, lcol_len - p, ujk);
      }
    } else {
      UF_long *li = l_indices[dep_global];
      UF_long p = 0u;
      for (; p < lcol_len; ++p) {
        const UF_long row = li[p];
        if (row >= run_end) {
          break;
        }
        workspace[row - dep0] -= ujk * lx[p];
      }
      if (p < lcol_len) {
        kls_scatter_subtract(x, li + p, lx + p, lcol_len - p, ujk);
      }
    }
  }

  kls_egraph_record_btf_scalar_run_exec(shared, rows, entries);
  *up = start_pos + rows;
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
  kls_egraph_scatter_unscaled_input(solver, shared->values, x,
                                    solver->refactor_col_ptr[k],
                                    solver->refactor_col_ptr[k + 1u]);

  UF_long *ui = u_indices[k];
  const int32_t *ui32 =
    solver->refactor_u_indices32 != NULL
      ? solver->refactor_u_indices32[k] : NULL;
  double *ux = u_values[k];
  UF_long ucol_len = kls_egraph_ulen_at(solver, k);
  UF_long up = 0;

  const int snode_batches_allowed =
    shared->snode_run_end != NULL;

  while (up < ucol_len) {

    if (snode_batches_allowed) {
      const UF_long consumed = kls_snode_batch_consume_cached(
        l_indices, l_values, numeric->Llen, ui, ui32, ux, ucol_len, up, x,
        0u, k, shared->snode_run_end, shared, wait_for_dependencies);
      if (consumed != 0u) {
        up += consumed;
        continue;
      }
    }
    if (snode_batches_allowed && solver->padded_run_of != NULL &&
        solver->padded_active && symbolic->nblocks == 1u) {
      const UF_long j0 = ui32 != NULL ? (UF_long)ui32[up] : ui[up];
      const UF_long r1 = solver->padded_run_of[j0];
      int run_ready = r1 != 0u;
      if (run_ready && wait_for_dependencies) {
        const UF_long re0 = solver->padded_run_start[r1 - 1u] +
                            solver->padded_run_len[r1 - 1u];
        for (UF_long c = j0; c < re0; ++c) {
          if (!kls_egraph_refactor_dependency_done_now(shared, c)) {
            run_ready = 0;
            break;
          }
        }
      }
      if (run_ready) {
        const UF_long consumed = kls_padded_run_consume(
          solver, 0u, j0, ui, ux, ucol_len, up, x);
        if (consumed != 0u) {
          up += consumed;
          continue;
        }
      }
    }
    const UF_long j = ui32 != NULL ? (UF_long)ui32[up] : ui[up];
    if (wait_for_dependencies &&
        !kls_egraph_refactor_wait_done(shared, j)) {
      return 0;
    }
    const double ujk = x[j];
    x[j] = 0.0;
    ux[up] = ujk;

    if (ujk != 0.0) {
      UF_long *li = l_indices[j];
      double *lx = l_values[j];
      UF_long lcol_len = kls_egraph_llen_at(solver, j);
      kls_scatter_subtract_refactor_l(solver, x, j, li, lx, lcol_len, ujk);
    }
    up++;
  }

  const double ukk = x[k];
  x[k] = 0.0;
  if (ukk == 0.0) {
    kls_egraph_refactor_record_singular(shared, k, symbolic->Q[k]);
    if (solver->common.halt_if_singular) {
      return 0;
    }
  }
  kls_egraph_store_udiag(shared, udiag, k, ukk);

  UF_long *li = l_indices[k];
  double *lx = l_values[k];
  UF_long lcol_len = kls_egraph_llen_at(solver, k);
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
  kls_egraph_store_l_column_from_workspace(solver, x, k, li, lx, lcol_len,
                                           ukk);

  return 1;
}

/* Once the realized single-block executor has declined every representation-
   changing numeric surface, keep its scalar dependency walk free of the
   corresponding per-dependency tests.  The capability gate is computed once
   per refactor dispatch; this kernel retains the cached strict-supernode
   consumer and the exact scalar arithmetic/storage order of the general
   single-block path. */
static int kls_egraph_refactor_single_unscaled_plain_column(
  kls_egraph_refactor_worker *worker,
  UF_long k,
  int wait_for_dependencies) {
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  trilinos_klu_l_numeric *numeric = solver->numeric;
  const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  double *restrict x = worker->x;
  double *restrict udiag = (double *)numeric->Udiag;
  UF_long **restrict l_indices = solver->refactor_l_indices;
  int32_t **restrict l_indices32 = solver->refactor_l_indices32;
  double **restrict l_values = solver->refactor_l_values;
  UF_long **restrict u_indices = solver->refactor_u_indices;
  double **restrict u_values = solver->refactor_u_values;
  if (solver->refactor_lu_pointer_count != solver->n ||
      l_indices == NULL || l_values == NULL ||
      u_indices == NULL || u_values == NULL) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }

  kls_egraph_scatter_unscaled_input(solver, shared->values, x,
                                    solver->refactor_col_ptr[k],
                                    solver->refactor_col_ptr[k + 1u]);

  UF_long *restrict ui = u_indices[k];
  const int32_t *restrict ui32 =
    solver->refactor_u_indices32 != NULL
      ? solver->refactor_u_indices32[k] : NULL;
  double *restrict ux = u_values[k];
  const UF_long ucol_len = kls_egraph_ulen_at(solver, k);
  UF_long up = 0u;
  while (up < ucol_len) {
    const UF_long consumed = kls_snode_batch_consume_cached(
      l_indices, l_values, numeric->Llen, ui, ui32, ux, ucol_len, up, x,
      0u, k, shared->snode_run_end, shared, wait_for_dependencies);
    if (consumed != 0u) {
      up += consumed;
      continue;
    }

    const UF_long j = ui32 != NULL ? (UF_long)ui32[up] : ui[up];
    if (wait_for_dependencies && shared->pipeline_done != NULL &&
        atomic_load_explicit(&shared->pipeline_done[j],
                             memory_order_acquire) !=
          shared->pipeline_generation &&
        !kls_egraph_refactor_wait_done(shared, j)) {
      return 0;
    }
    const double ujk = x[j];
    x[j] = 0.0;
    ux[up] = ujk;
    if (ujk != 0.0) {
      const UF_long lcol_len = kls_egraph_llen_at(solver, j);
      const int32_t *restrict li32 =
        l_indices32 != NULL ? l_indices32[j] : NULL;
      if (li32 != NULL || lcol_len == 0u) {
        kls_scatter_subtract_refactor_i32(
          solver, x, li32, l_values[j], lcol_len, ujk);
      } else {
        kls_scatter_subtract(x, l_indices[j], l_values[j], lcol_len, ujk);
      }
    }
    up++;
  }

  const double ukk = x[k];
  x[k] = 0.0;
  if (ukk == 0.0) {
    kls_egraph_refactor_record_singular(shared, k, symbolic->Q[k]);
    if (solver->common.halt_if_singular) {
      return 0;
    }
  }
  kls_egraph_store_udiag(shared, udiag, k, ukk);
  kls_egraph_store_l_column_from_workspace(
    solver, x, k, l_indices[k], l_values[k],
    kls_egraph_llen_at(solver, k), ukk);
  return 1;
}

/* Two READY pipeline columns (readiness implies mutual independence)
   processed with one merged dependency walk so producer runs shared by
   both U patterns stream their tail panels once.  Plain configuration
   only; kls_egraph_pair_dispatch_allowed gates the call. */
static int kls_egraph_refactor_single_unscaled_column_pair(
  kls_egraph_refactor_worker *worker,
  UF_long ka,
  UF_long kb) {
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  trilinos_klu_l_numeric *numeric = solver->numeric;
  const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  double *xa = worker->x;
  double *xb = worker->pair_x;
  double *udiag = (double *)numeric->Udiag;
  UF_long **l_indices = solver->refactor_l_indices;
  double **l_values = solver->refactor_l_values;
  UF_long **u_indices = solver->refactor_u_indices;
  double **u_values = solver->refactor_u_values;
  if (solver->refactor_lu_pointer_count != solver->n ||
      l_indices == NULL || l_values == NULL ||
      u_indices == NULL || u_values == NULL || xb == NULL) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }
  kls_egraph_scatter_unscaled_input(solver, shared->values, xa,
                                    solver->refactor_col_ptr[ka],
                                    solver->refactor_col_ptr[ka + 1u]);
  kls_egraph_scatter_unscaled_input(solver, shared->values, xb,
                                    solver->refactor_col_ptr[kb],
                                    solver->refactor_col_ptr[kb + 1u]);

  UF_long *uia = u_indices[ka];
  UF_long *uib = u_indices[kb];
  const int32_t *uia32 = solver->refactor_u_indices32 != NULL
    ? solver->refactor_u_indices32[ka] : NULL;
  const int32_t *uib32 = solver->refactor_u_indices32 != NULL
    ? solver->refactor_u_indices32[kb] : NULL;
  double *uxa = u_values[ka];
  double *uxb = u_values[kb];
  const UF_long ulen_a = numeric->Ulen[ka];
  const UF_long ulen_b = numeric->Ulen[kb];
  UF_long upa = 0;
  UF_long upb = 0;
  const UF_long producer_limit = ka < kb ? ka : kb;

  /* Cooperative walk: leased pipeline columns are not necessarily ready,
     and one side of the pair may depend (transitively) on the other, so
     no step may block.  Each side advances only through done_now
     dependencies; a completed side finalizes and is marked done
     IMMEDIATELY so chains through it (including the partner's) keep
     moving.  Liveness follows the single-column argument: the earliest
     not-done column's dependencies are all done, so its owner always
     progresses. */
  int a_done = 0;
  int b_done = 0;
  unsigned spin = 0;
  while (!a_done || !b_done) {
    if (kls_egraph_refactor_should_stop(shared)) {
      return 0;
    }
    int progress = 0;
    if (upa < ulen_a && upb < ulen_b) {
      const UF_long ja = uia32 != NULL ? (UF_long)uia32[upa] : uia[upa];
      const UF_long jb = uib32 != NULL ? (UF_long)uib32[upb] : uib[upb];
      if (ja == jb) {
        const UF_long consumed = kls_snode_batch_consume_cached_pair(
          l_indices, l_values, numeric->Llen, uia, uia32, uxa, ulen_a, upa,
          xa, uib, uib32, uxb, ulen_b, upb, xb, 0u, producer_limit,
          shared->snode_run_end, shared);
        if (consumed != 0u) {
          upa += consumed;
          upb += consumed;
          progress = 1;
        }
      }
    }
    if (!progress && upa < ulen_a) {
      const UF_long consumed = kls_snode_batch_consume_cached(
        l_indices, l_values, numeric->Llen, uia, uia32, uxa, ulen_a, upa,
        xa, 0u, ka, shared->snode_run_end, shared, 1);
      if (consumed != 0u) {
        upa += consumed;
        progress = 1;
      } else {
        const UF_long ja = uia32 != NULL ? (UF_long)uia32[upa] : uia[upa];
        if (kls_egraph_refactor_dependency_done_now(shared, ja)) {
          const double ujk = xa[ja];
          xa[ja] = 0.0;
          uxa[upa] = ujk;
          if (ujk != 0.0) {
            kls_scatter_subtract_refactor_l(solver, xa, ja, l_indices[ja],
                                            l_values[ja], numeric->Llen[ja],
                                            ujk);
          }
          upa++;
          progress = 1;
        }
      }
    }
    if (!progress && upb < ulen_b) {
      const UF_long consumed = kls_snode_batch_consume_cached(
        l_indices, l_values, numeric->Llen, uib, uib32, uxb, ulen_b, upb,
        xb, 0u, kb, shared->snode_run_end, shared, 1);
      if (consumed != 0u) {
        upb += consumed;
        progress = 1;
      } else {
        const UF_long jb = uib32 != NULL ? (UF_long)uib32[upb] : uib[upb];
        if (kls_egraph_refactor_dependency_done_now(shared, jb)) {
          const double ujk = xb[jb];
          xb[jb] = 0.0;
          uxb[upb] = ujk;
          if (ujk != 0.0) {
            kls_scatter_subtract_refactor_l(solver, xb, jb, l_indices[jb],
                                            l_values[jb], numeric->Llen[jb],
                                            ujk);
          }
          upb++;
          progress = 1;
        }
      }
    }
    if (!a_done && upa >= ulen_a) {
      const double ukka = xa[ka];
      xa[ka] = 0.0;
      if (ukka == 0.0) {
        kls_egraph_refactor_record_singular(shared, ka, symbolic->Q[ka]);
        if (solver->common.halt_if_singular) {
          return 0;
        }
      }
      udiag[ka] = ukka;
      kls_egraph_store_l_column_from_workspace(solver, xa, ka,
                                               l_indices[ka], l_values[ka],
                                               numeric->Llen[ka], ukka);

      kls_egraph_refactor_mark_done(shared, ka);
      a_done = 1;
      progress = 1;
    }
    if (!b_done && upb >= ulen_b) {
      const double ukkb = xb[kb];
      xb[kb] = 0.0;
      if (ukkb == 0.0) {
        kls_egraph_refactor_record_singular(shared, kb, symbolic->Q[kb]);
        if (solver->common.halt_if_singular) {
          return 0;
        }
      }
      udiag[kb] = ukkb;
      kls_egraph_store_l_column_from_workspace(solver, xb, kb,
                                               l_indices[kb], l_values[kb],
                                               numeric->Llen[kb], ukkb);

      kls_egraph_refactor_mark_done(shared, kb);
      b_done = 1;
      progress = 1;
    }
    if (!progress) {
      kls_cpu_relax();
      kls_egraph_pipeline_pause(&spin);
    }
  }
  return 1;
}

/* k-way fused column walk (k <= KLS_SNODE_MAX_FUSE), the pair walk
   generalized over per-side arrays. Same cooperative contract: leased
   pipeline columns are not necessarily ready and may depend on each
   other, so no step blocks; every side advances only through done_now
   dependencies, aligned sides consume runs fused (panel read once for
   k SPAs), and a completed side finalizes and marks done IMMEDIATELY
   so chains through it keep moving. */
static int kls_egraph_refactor_single_unscaled_column_multi(
  kls_egraph_refactor_worker *worker,
  const UF_long *cols,
  int nc) {
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  trilinos_klu_l_numeric *numeric = solver->numeric;
  const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  double *udiag = (double *)numeric->Udiag;
  UF_long **l_indices = solver->refactor_l_indices;
  double **l_values = solver->refactor_l_values;
  UF_long **u_indices = solver->refactor_u_indices;
  double **u_values = solver->refactor_u_values;
  double *xv[KLS_SNODE_MAX_FUSE];
  if (solver->refactor_lu_pointer_count != solver->n ||
      l_indices == NULL || l_values == NULL ||
      u_indices == NULL || u_values == NULL || nc < 2 ||
      nc > KLS_SNODE_MAX_FUSE || worker->pair_x == NULL ||
      (nc > 2 && worker->fuse_x == NULL)) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }
  xv[0] = worker->x;
  xv[1] = worker->pair_x;
  for (int c = 2; c < nc; ++c) {
    xv[c] = worker->fuse_x + (size_t)(c - 2) * (size_t)solver->n;
  }
  const UF_long *ui[KLS_SNODE_MAX_FUSE];
  const int32_t *ui32[KLS_SNODE_MAX_FUSE];
  double *ux[KLS_SNODE_MAX_FUSE];
  UF_long ulen[KLS_SNODE_MAX_FUSE];
  UF_long up[KLS_SNODE_MAX_FUSE];
  int done[KLS_SNODE_MAX_FUSE];
  for (int c = 0; c < nc; ++c) {
    const UF_long k = cols[c];
    kls_egraph_scatter_unscaled_input(solver, shared->values, xv[c],
                                      solver->refactor_col_ptr[k],
                                      solver->refactor_col_ptr[k + 1u]);
    ui[c] = u_indices[k];
    ui32[c] = solver->refactor_u_indices32 != NULL
      ? solver->refactor_u_indices32[k] : NULL;
    ux[c] = u_values[k];
    ulen[c] = numeric->Ulen[k];
    up[c] = 0;
    done[c] = 0;
  }
  int remaining = nc;
  unsigned spin = 0;
  while (remaining > 0) {
    if (kls_egraph_refactor_should_stop(shared)) {
      return 0;
    }
    int progress = 0;
    const int subset_fusion = solver->eg_subset_choice > 0 ||
      solver->eg_subset_pending == 2;
    if (subset_fusion) {
      /* Fuse every aligned subset in the deck.  Requiring all four consumers
         to share the same current producer made one divergent U stream force
         four independent reads of an otherwise common L panel. */
      unsigned int grouped = 0u;
      for (int seed = 0; seed < nc; ++seed) {
        if ((grouped & (1u << seed)) != 0u || done[seed] ||
            up[seed] >= ulen[seed]) {
          continue;
        }
      const UF_long j0 = ui32[seed] != NULL
        ? (UF_long)ui32[seed][up[seed]] : ui[seed][up[seed]];
      const UF_long *group_ui[KLS_SNODE_MAX_FUSE];
      const int32_t *group_ui32[KLS_SNODE_MAX_FUSE];
      double *group_ux[KLS_SNODE_MAX_FUSE];
      double *group_xv[KLS_SNODE_MAX_FUSE];
      UF_long group_ulen[KLS_SNODE_MAX_FUSE];
      UF_long group_up[KLS_SNODE_MAX_FUSE];
      int members[KLS_SNODE_MAX_FUSE];
      int group_count = 0;
      UF_long group_producer_limit = solver->n;
      for (int c = seed; c < nc; ++c) {
        if ((grouped & (1u << c)) != 0u || done[c] || up[c] >= ulen[c]) {
          continue;
        }
        const UF_long jc = ui32[c] != NULL
          ? (UF_long)ui32[c][up[c]] : ui[c][up[c]];
        if (jc != j0) {
          continue;
        }
        members[group_count] = c;
        group_ui[group_count] = ui[c];
        group_ui32[group_count] = ui32[c];
        group_ux[group_count] = ux[c];
        group_xv[group_count] = xv[c];
        group_ulen[group_count] = ulen[c];
        group_up[group_count] = up[c];
        if (cols[c] < group_producer_limit) {
          group_producer_limit = cols[c];
        }
        grouped |= 1u << c;
        group_count++;
      }
      if (group_count >= 2) {
        const UF_long consumed = kls_snode_batch_consume_cached_multi(
          l_indices, l_values, numeric->Llen, group_count, group_ui,
          group_ui32, group_ux, group_ulen, group_up, group_xv, 0u,
          group_producer_limit, shared->snode_run_end, shared);
        if (consumed != 0u) {
          for (int g = 0; g < group_count; ++g) {
            up[members[g]] += consumed;
          }
          progress = 1;
        }
      }
    }
    } else {
      int all_alive = 1;
      UF_long all_producer_limit = solver->n;
      for (int c = 0; c < nc; ++c) {
        if (done[c] || up[c] >= ulen[c]) {
          all_alive = 0;
          break;
        }
        if (cols[c] < all_producer_limit) {
          all_producer_limit = cols[c];
        }
      }
      if (all_alive) {
        const UF_long j0 = ui32[0] != NULL
          ? (UF_long)ui32[0][up[0]] : ui[0][up[0]];
        int aligned = 1;
        for (int c = 1; c < nc; ++c) {
          const UF_long jc = ui32[c] != NULL
            ? (UF_long)ui32[c][up[c]] : ui[c][up[c]];
          if (jc != j0) {
            aligned = 0;
            break;
          }
        }
        if (aligned) {
          const UF_long consumed = kls_snode_batch_consume_cached_multi(
            l_indices, l_values, numeric->Llen, nc, ui, ui32, ux, ulen,
            up, xv, 0u, all_producer_limit, shared->snode_run_end, shared);
          if (consumed != 0u) {
            for (int c = 0; c < nc; ++c) {
              up[c] += consumed;
            }
            progress = 1;
          }
        }
      }
    }
    /* When the current producer positions are not common, advance every
       ready consumer in this deck before polling the dependency frontier
       again.  The former first-progress break revisited the same four-side
       deck up to four times and repeated the stop/alignment scans even
       though the SPAs are independent.  Each side still advances strictly
       in U order and only through published dependencies. */
    for (int c = 0; c < nc; ++c) {
      if (done[c] || up[c] >= ulen[c]) {
        continue;
      }
      const UF_long consumed = kls_snode_batch_consume_cached(
        l_indices, l_values, numeric->Llen, ui[c], ui32[c], ux[c],
        ulen[c], up[c], xv[c], 0u, cols[c], shared->snode_run_end,
        shared, 1);
      if (consumed != 0u) {
        up[c] += consumed;
        progress = 1;
        continue;
      }
      const UF_long jc = ui32[c] != NULL
        ? (UF_long)ui32[c][up[c]] : ui[c][up[c]];
      if (kls_egraph_refactor_dependency_done_now(shared, jc)) {
        const double ujk = xv[c][jc];
        xv[c][jc] = 0.0;
        ux[c][up[c]] = ujk;
        if (ujk != 0.0) {
          kls_scatter_subtract_refactor_l(solver, xv[c], jc,
                                          l_indices[jc], l_values[jc],
                                          numeric->Llen[jc], ujk);
        }
        up[c]++;
        progress = 1;
      }
    }
    for (int c = 0; c < nc; ++c) {
      if (!done[c] && up[c] >= ulen[c]) {
        const UF_long k = cols[c];
        const double ukk = xv[c][k];
        xv[c][k] = 0.0;
        if (ukk == 0.0) {
          kls_egraph_refactor_record_singular(shared, k, symbolic->Q[k]);
          if (solver->common.halt_if_singular) {
            return 0;
          }
        }
        udiag[k] = ukk;
        kls_egraph_store_l_column_from_workspace(solver, xv[c], k,
                                                 l_indices[k],
                                                 l_values[k],
                                                 numeric->Llen[k], ukk);

        kls_egraph_refactor_mark_done(shared, k);
        done[c] = 1;
        remaining--;
        progress = 1;
      }
    }
    if (!progress) {
      kls_cpu_relax();
      kls_egraph_pipeline_pause(&spin);
    }
  }
  return 1;
}

static int kls_egraph_refactor_single_scaled_column(
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
      u_indices == NULL || u_values == NULL ||
      solver->row_idx == NULL || shared->rs == NULL) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }

  if (!kls_egraph_scatter_refactor_input(
        solver, shared->values, shared->rs, shared->scale, x, 0u, solver->n,
        solver->refactor_col_ptr[k], solver->refactor_col_ptr[k + 1u])) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }

  UF_long *ui = u_indices[k];
  const int32_t *ui32 =
    solver->refactor_u_indices32 != NULL
      ? solver->refactor_u_indices32[k] : NULL;
  double *ux = u_values[k];
  UF_long ucol_len = numeric->Ulen[k];
  UF_long up = 0;

  const int snode_batches_allowed =
    shared->snode_run_end != NULL;

  while (up < ucol_len) {

    if (snode_batches_allowed) {
      const UF_long consumed = kls_snode_batch_consume_cached(
        l_indices, l_values, numeric->Llen, ui, ui32, ux, ucol_len, up, x,
        0u, k, shared->snode_run_end, shared, wait_for_dependencies);
      if (consumed != 0u) {
        up += consumed;
        continue;
      }
    }
    if (snode_batches_allowed && solver->padded_run_of != NULL &&
        solver->padded_active && symbolic->nblocks == 1u) {
      const UF_long j0 = ui32 != NULL ? (UF_long)ui32[up] : ui[up];
      const UF_long r1 = solver->padded_run_of[j0];
      int run_ready = r1 != 0u;
      if (run_ready && wait_for_dependencies) {
        const UF_long re0 = solver->padded_run_start[r1 - 1u] +
                            solver->padded_run_len[r1 - 1u];
        for (UF_long c = j0; c < re0; ++c) {
          if (!kls_egraph_refactor_dependency_done_now(shared, c)) {
            run_ready = 0;
            break;
          }
        }
      }
      if (run_ready) {
        const UF_long consumed = kls_padded_run_consume(
          solver, 0u, j0, ui, ux, ucol_len, up, x);
        if (consumed != 0u) {
          up += consumed;
          continue;
        }
      }
    }
    const UF_long j = ui32 != NULL ? (UF_long)ui32[up] : ui[up];
    if (wait_for_dependencies &&
        !kls_egraph_refactor_wait_done(shared, j)) {
      return 0;
    }
    const double ujk = x[j];
    x[j] = 0.0;
    ux[up] = ujk;

    if (ujk != 0.0) {
      UF_long *li = l_indices[j];
      double *lx = l_values[j];
      UF_long lcol_len = numeric->Llen[j];
      kls_scatter_subtract_refactor_l(solver, x, j, li, lx, lcol_len, ujk);
    }
    up++;
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
  kls_egraph_store_l_column_from_workspace(solver, x, k, li, lx, lcol_len,
                                           ukk);

  return 1;
}

/* Cluster barriers publish all predecessors. This narrow no-wait kernel
   retains cached-supernode batches and the general kernel's scalar order,
   without testing the checked, supernode-update, or scalar-run modes. */
static int kls_egraph_btf_plain_cluster_allowed(
  const kls_egraph_refactor_shared *shared) {
  return shared != NULL && !shared->check_pivots &&
         !shared->btf_scalar_run_exec;
}

static int kls_egraph_refactor_btf_unscaled_plain_cluster_column(
  kls_egraph_refactor_worker *worker,
  UF_long k,
  int wait_for_dependencies) {
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  trilinos_klu_l_numeric *numeric = solver->numeric;
  const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  double *restrict x = worker->x;
  double *restrict udiag = (double *)numeric->Udiag;
  UF_long **restrict l_indices = solver->refactor_l_indices;
  int32_t **restrict l_indices32 = solver->refactor_l_indices32;
  double **restrict l_values = solver->refactor_l_values;
  UF_long **restrict u_indices = solver->refactor_u_indices;
  double **restrict u_values = solver->refactor_u_values;
  if (solver->refactor_lu_pointer_count != solver->n ||
      l_indices == NULL || l_values == NULL ||
      u_indices == NULL || u_values == NULL) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }

  const UF_long block = solver->refactor_col_block[k];
  const UF_long k1 = symbolic->R[block];
  const UF_long k2 = symbolic->R[block + 1u];
  const UF_long nk = k2 - k1;
  const UF_long local_k = k - k1;
  UF_long poff = numeric->Offp[k];
  const UF_long poff_end = numeric->Offp[k + 1u];
  if (!kls_egraph_copy_unscaled_offblock_input(
        solver, shared->values, (double *)numeric->Offx, &poff, poff_end,
        solver->refactor_col_ptr[k], solver->refactor_block_start[k])) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }

  if (nk == 1u) {
    double pivot = 0.0;
    for (UF_long p = solver->refactor_block_start[k];
         p < solver->refactor_col_ptr[k + 1u]; ++p) {
      pivot = kls_egraph_unscaled_input_value_at(solver, shared->values, p);
    }
    kls_egraph_store_udiag(shared, udiag, k, pivot);
    if (pivot == 0.0) {
      kls_egraph_refactor_record_singular(shared, k, symbolic->Q[k]);
      if (solver->common.halt_if_singular) {
        return 0;
      }
    }
    return 1;
  }

  UF_long *restrict ui = u_indices[k];
  const int32_t *restrict ui32 =
    solver->refactor_u_indices32 != NULL
      ? solver->refactor_u_indices32[k] : NULL;
  double *restrict ux = u_values[k];
  const UF_long ucol_len = numeric->Ulen[k];
  const UF_long *restrict llen = numeric->Llen + k1;
  const int snode_batches_allowed =
    shared->snode_run_end != NULL;

  kls_egraph_scatter_btf_unscaled_input(
    solver, shared->values, x, k1, solver->refactor_block_start[k],
    solver->refactor_col_ptr[k + 1u]);

  UF_long up = 0u;
  while (up < ucol_len) {
    if (snode_batches_allowed) {
      const UF_long consumed = kls_snode_batch_consume_cached(
        l_indices, l_values, llen, ui, ui32, ux, ucol_len, up, x, k1,
        local_k, shared->snode_run_end, shared, wait_for_dependencies);
      if (consumed != 0u) {
        up += consumed;
        continue;
      }
    }
    const UF_long j = ui32 != NULL ? (UF_long)ui32[up] : ui[up];
    if (wait_for_dependencies && shared->pipeline_done != NULL &&
        atomic_load_explicit(&shared->pipeline_done[k1 + j],
                             memory_order_acquire) !=
          shared->pipeline_generation) {
      /* The general wait helper also polls stop/help state.  Keep it for an
         actually late producer; the overwhelmingly common published case
         needs only this inlined acquire instead of a non-inlined call. */
      if (!kls_egraph_refactor_wait_done(shared, k1 + j)) {
        return 0;
      }
    }
    const double ujk = x[j];
    x[j] = 0.0;
    ux[up] = ujk;
    if (ujk != 0.0) {
      const UF_long global_j = k1 + j;
      const UF_long lcol_len = llen[j];
      const int32_t *restrict li32 =
        l_indices32 != NULL ? l_indices32[global_j] : NULL;
      if (li32 != NULL || lcol_len == 0u) {
        kls_scatter_subtract_refactor_i32(
          solver, x, li32, l_values[global_j], lcol_len, ujk);
      } else {
        kls_scatter_subtract(x, l_indices[global_j], l_values[global_j],
                             lcol_len, ujk);
      }
    }
    up++;
  }

  const double ukk = x[local_k];
  x[local_k] = 0.0;
  if (ukk == 0.0) {
    kls_egraph_refactor_record_singular(shared, k, symbolic->Q[k]);
    if (solver->common.halt_if_singular) {
      return 0;
    }
  }
  kls_egraph_store_udiag(shared, udiag, k, ukk);
  UF_long *restrict li = l_indices[k];
  double *restrict lx = l_values[k];
  const UF_long lcol_len = llen[local_k];
  kls_egraph_store_l_column_from_workspace(solver, x, k, li, lx, lcol_len,
                                           ukk);
  return 1;
}

__attribute__((aligned(32)))
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
  int32_t **l_indices32 = solver->refactor_l_indices32;
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
  if (!kls_egraph_copy_unscaled_offblock_input(
        solver, shared->values, offx, &poff, poff_end,
        solver->refactor_col_ptr[k], solver->refactor_block_start[k])) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }

  if (nk == 1u) {
    double pivot = 0.0;
    for (UF_long p = solver->refactor_block_start[k];
         p < solver->refactor_col_ptr[k + 1u]; ++p) {
      pivot = kls_egraph_unscaled_input_value_at(solver, shared->values, p);
    }
    kls_egraph_store_udiag(shared, udiag, k, pivot);
    if (pivot == 0.0) {
      kls_egraph_refactor_record_singular(shared, k, symbolic->Q[k]);
      if (solver->common.halt_if_singular) {
        return 0;
      }
    }
    return 1;
  }

  UF_long *ui = u_indices[k];
  const int32_t *ui32 =
    solver->refactor_u_indices32 != NULL
      ? solver->refactor_u_indices32[k] : NULL;
  double *ux = u_values[k];
  UF_long ucol_len = numeric->Ulen[k];

  const int snode_batches_allowed =
    shared->snode_run_end != NULL;
  const UF_long *llen = numeric->Llen + k1;

  UF_long up = 0;
  kls_egraph_scatter_btf_unscaled_input(
    solver, shared->values, x, k1, solver->refactor_block_start[k],
    solver->refactor_col_ptr[k + 1u]);

  /* Cluster columns already have their predecessors published; keep that
     CKTSO-style cluster loop free of the per-entry wait check used below. */
  if (wait_for_dependencies) {
    while (up < ucol_len) {
      if (snode_batches_allowed) {
        const UF_long consumed = kls_snode_batch_consume_cached(
          l_indices, l_values, llen, ui, ui32, ux, ucol_len, up, x, k1,
          local_k, shared->snode_run_end, shared, 1);
        if (consumed != 0u) {
          up += consumed;
          continue;
        }
      }
      if (snode_batches_allowed && solver->padded_run_of != NULL &&
          solver->padded_active) {
        const UF_long j0 = ui32 != NULL ? (UF_long)ui32[up] : ui[up];
        const UF_long r1 = solver->padded_run_of[k1 + j0];
        int run_ready = r1 != 0u;
        if (run_ready) {
          const UF_long re0 = solver->padded_run_start[r1 - 1u] +
                              solver->padded_run_len[r1 - 1u];
          for (UF_long c = k1 + j0; c < re0; ++c) {
            if (!kls_egraph_refactor_dependency_done_now(shared, c)) {
              run_ready = 0;
              break;
            }
          }
        }
        if (run_ready) {
          const UF_long consumed = kls_padded_run_consume(
            solver, k1, j0, ui, ux, ucol_len, up, x);
          if (consumed != 0u) {
            up += consumed;
            continue;
          }
        }
      }
      const UF_long j = ui32 != NULL ? (UF_long)ui32[up] : ui[up];
      const int run_status =
        kls_egraph_refactor_try_btf_scalar_producer_run(
          worker, k1, local_k, &up, ucol_len, ui, ui32, ux, l_indices,
          l_indices32, l_values, llen, x, 1);
      if (run_status < 0) {
        return 0;
      }
      if (run_status > 0) {
        continue;
      }

      if (!kls_egraph_refactor_wait_done(shared, k1 + j)) {

        return 0;
      }

      kls_egraph_refactor_apply_btf_scalar_dep(
        shared, k1, j, up, ux, l_indices, l_indices32, l_values, llen, x);
      up++;
    }
  } else {
    while (up < ucol_len) {
      if (snode_batches_allowed) {
        const UF_long consumed = kls_snode_batch_consume_cached(
          l_indices, l_values, llen, ui, ui32, ux, ucol_len, up, x, k1,
          local_k, shared->snode_run_end, shared, 0);
        if (consumed != 0u) {
          up += consumed;
          continue;
        }
      }
      const int run_status =
        kls_egraph_refactor_try_btf_scalar_producer_run(
          worker, k1, local_k, &up, ucol_len, ui, ui32, ux, l_indices,
          l_indices32, l_values, llen, x, 0);
      if (run_status < 0) {
        return 0;
      }
      if (run_status > 0) {
        continue;
      }
      const UF_long j = ui32 != NULL ? (UF_long)ui32[up] : ui[up];

      kls_egraph_refactor_apply_btf_scalar_dep(
        shared, k1, j, up, ux, l_indices, l_indices32, l_values, llen, x);
      up++;
    }
  }

  const double ukk = x[local_k];
  x[local_k] = 0.0;
  if (ukk == 0.0) {
    kls_egraph_refactor_record_singular(shared, k, symbolic->Q[k]);
    if (solver->common.halt_if_singular) {
      return 0;
    }
  }
  kls_egraph_store_udiag(shared, udiag, k, ukk);

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
  kls_egraph_store_l_column_from_workspace(solver, x, k, li, lx, lcol_len,
                                           ukk);

  return 1;
}

/* One L-column stream feeding two consumer SPAs: the scalar-dep twin of
   the fused batch (2 fmas per streamed value + one shared index load). */
static KLS_ALWAYS_INLINE void kls_scatter_subtract_pair_i32(
  double *restrict xa,
  double *restrict xb,
  const int32_t *restrict rows,
  const double *restrict values,
  UF_long length,
  double ua,
  double ub) {
  for (UF_long p = 0; p < length; ++p) {
    const UF_long row = (UF_long)rows[p];
    const double v = values[p];
    xa[row] -= ua * v;
    xb[row] -= ub * v;
  }
}

static KLS_ALWAYS_INLINE void kls_scatter_subtract_pair(
  double *restrict xa,
  double *restrict xb,
  const UF_long *restrict rows,
  const double *restrict values,
  UF_long length,
  double ua,
  double ub) {
  for (UF_long p = 0; p < length; ++p) {
    const UF_long row = rows[p];
    const double v = values[p];
    xa[row] -= ua * v;
    xb[row] -= ub * v;
  }
}

/* BTF twin of the single-unscaled pair: two leased same-block columns
   walked cooperatively so shared producer runs stream once.  The group
   state machinery is informed via the existing dispatch-bypass branch;
   the fp64 scalar-run exec is skipped (it either blocks or assumes done
   dependencies) - declined batch positions go through done_now-gated
   scalars.  Caller guarantees the plain flag surface and same block. */
static int kls_egraph_refactor_btf_unscaled_column_pair(
  kls_egraph_refactor_worker *worker,
  UF_long ka,
  UF_long kb) {
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  trilinos_klu_l_numeric *numeric = solver->numeric;
  const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  double *xa = worker->x;
  double *xb = worker->pair_x;
  double *udiag = (double *)numeric->Udiag;
  UF_long **l_indices = solver->refactor_l_indices;
  int32_t **l_indices32 = solver->refactor_l_indices32;
  double **l_values = solver->refactor_l_values;
  UF_long **u_indices = solver->refactor_u_indices;
  double **u_values = solver->refactor_u_values;
  if (solver->refactor_lu_pointer_count != solver->n ||
      l_indices == NULL || l_values == NULL ||
      u_indices == NULL || u_values == NULL || xb == NULL) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }
  const UF_long block = solver->refactor_col_block[ka];
  const UF_long k1 = symbolic->R[block];
  const UF_long local_a = ka - k1;
  const UF_long local_b = kb - k1;
  const UF_long *llen = numeric->Llen + k1;

  UF_long poff = numeric->Offp[ka];
  double *offx = (double *)numeric->Offx;
  if (!kls_egraph_copy_unscaled_offblock_input(
        solver, shared->values, offx, &poff, numeric->Offp[ka + 1u],
        solver->refactor_col_ptr[ka], solver->refactor_block_start[ka])) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }
  poff = numeric->Offp[kb];
  if (!kls_egraph_copy_unscaled_offblock_input(
        solver, shared->values, offx, &poff, numeric->Offp[kb + 1u],
        solver->refactor_col_ptr[kb], solver->refactor_block_start[kb])) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }

  kls_egraph_scatter_btf_unscaled_input(
    solver, shared->values, xa, k1, solver->refactor_block_start[ka],
    solver->refactor_col_ptr[ka + 1u]);
  kls_egraph_scatter_btf_unscaled_input(
    solver, shared->values, xb, k1, solver->refactor_block_start[kb],
    solver->refactor_col_ptr[kb + 1u]);

  UF_long *uia = u_indices[ka];
  UF_long *uib = u_indices[kb];
  const int32_t *uia32 = solver->refactor_u_indices32 != NULL
    ? solver->refactor_u_indices32[ka] : NULL;
  const int32_t *uib32 = solver->refactor_u_indices32 != NULL
    ? solver->refactor_u_indices32[kb] : NULL;
  double *uxa = u_values[ka];
  double *uxb = u_values[kb];
  const UF_long ulen_a = numeric->Ulen[ka];
  const UF_long ulen_b = numeric->Ulen[kb];
  UF_long upa = 0;
  UF_long upb = 0;
  const UF_long producer_limit = local_a < local_b ? local_a : local_b;

  int a_done = 0;
  int b_done = 0;
  unsigned spin = 0;
  while (!a_done || !b_done) {
    if (kls_egraph_refactor_should_stop(shared)) {
      return 0;
    }
    int progress = 0;
    if (upa < ulen_a && upb < ulen_b) {
      const UF_long ja = uia32 != NULL ? (UF_long)uia32[upa] : uia[upa];
      const UF_long jb = uib32 != NULL ? (UF_long)uib32[upb] : uib[upb];
      if (ja == jb) {
        const UF_long consumed = kls_snode_batch_consume_cached_pair(
          l_indices, l_values, llen, uia, uia32, uxa, ulen_a, upa, xa,
          uib, uib32, uxb, ulen_b, upb, xb, k1, producer_limit,
          shared->snode_run_end, shared);
        if (consumed != 0u) {
          upa += consumed;
          upb += consumed;
          progress = 1;
        } else if (kls_egraph_refactor_dependency_done_now(shared,
                                                           k1 + ja)) {
          const double ua = xa[ja];
          const double ub = xb[ja];
          xa[ja] = 0.0;
          xb[ja] = 0.0;
          uxa[upa] = ua;
          uxb[upb] = ub;
          if (ua != 0.0 || ub != 0.0) {
            const UF_long dep_global = k1 + ja;
            const int32_t *li32 =
              l_indices32 != NULL ? l_indices32[dep_global] : NULL;
            const double *lx = l_values[dep_global];
            const UF_long lcol_len = llen[ja];
            if (li32 != NULL || lcol_len == 0u) {
              kls_scatter_subtract_pair_i32(xa, xb, li32, lx, lcol_len, ua,
                                            ub);
            } else {
              kls_scatter_subtract_pair(xa, xb, l_indices[dep_global], lx,
                                        lcol_len, ua, ub);
            }
          }
          upa++;
          upb++;
          progress = 1;
        }
      }
    }
    if (!progress && upa < ulen_a) {
      const UF_long consumed = kls_snode_batch_consume_cached(
        l_indices, l_values, llen, uia, uia32, uxa, ulen_a, upa, xa, k1,
        local_a, shared->snode_run_end, shared, 1);
      if (consumed != 0u) {
        upa += consumed;
        progress = 1;
      } else {
        const UF_long ja = uia32 != NULL ? (UF_long)uia32[upa] : uia[upa];
        if (kls_egraph_refactor_dependency_done_now(shared, k1 + ja)) {
          kls_egraph_refactor_apply_btf_scalar_dep(
            shared, k1, ja, upa, uxa, l_indices, l_indices32, l_values, llen,
            xa);
          upa++;
          progress = 1;
        }
      }
    }
    if (!progress && upb < ulen_b) {
      const UF_long consumed = kls_snode_batch_consume_cached(
        l_indices, l_values, llen, uib, uib32, uxb, ulen_b, upb, xb, k1,
        local_b, shared->snode_run_end, shared, 1);
      if (consumed != 0u) {
        upb += consumed;
        progress = 1;
      } else {
        const UF_long jb = uib32 != NULL ? (UF_long)uib32[upb] : uib[upb];
        if (kls_egraph_refactor_dependency_done_now(shared, k1 + jb)) {
          kls_egraph_refactor_apply_btf_scalar_dep(
            shared, k1, jb, upb, uxb, l_indices, l_indices32, l_values, llen,
            xb);
          upb++;
          progress = 1;
        }
      }
    }
    if (!a_done && upa >= ulen_a) {
      const double ukka = xa[local_a];
      xa[local_a] = 0.0;
      if (ukka == 0.0) {
        kls_egraph_refactor_record_singular(shared, ka, symbolic->Q[ka]);
        if (solver->common.halt_if_singular) {
          return 0;
        }
      }
      kls_egraph_store_udiag(shared, udiag, ka, ukka);
      kls_egraph_store_l_column_from_workspace(solver, xa, ka,
                                               l_indices[ka], l_values[ka],
                                               llen[local_a], ukka);
      kls_egraph_refactor_mark_done(shared, ka);
      a_done = 1;
      progress = 1;
    }
    if (!b_done && upb >= ulen_b) {
      const double ukkb = xb[local_b];
      xb[local_b] = 0.0;
      if (ukkb == 0.0) {
        kls_egraph_refactor_record_singular(shared, kb, symbolic->Q[kb]);
        if (solver->common.halt_if_singular) {
          return 0;
        }
      }
      kls_egraph_store_udiag(shared, udiag, kb, ukkb);
      kls_egraph_store_l_column_from_workspace(solver, xb, kb,
                                               l_indices[kb], l_values[kb],
                                               llen[local_b], ukkb);
      kls_egraph_refactor_mark_done(shared, kb);
      b_done = 1;
      progress = 1;
    }
    if (!progress) {
      kls_cpu_relax();
      kls_egraph_pipeline_pause(&spin);
    }
  }
  return 1;
}

/* k-way counterpart of the BTF pair walk.  All columns belong to one BTF
   block, so their SPAs use the same local row frame.  Ready/pipeline leases
   are still only ownership claims: a side advances solely through published
   dependencies and publishes immediately on completion.  Consequently this
   is the same dependency-safe algorithm as the pair kernel, while a common
   producer's L stream can feed up to four independent consumer SPAs. */
static int kls_egraph_refactor_btf_unscaled_column_multi(
  kls_egraph_refactor_worker *worker,
  const UF_long *cols,
  int nc) {
  if (worker == NULL || worker->shared == NULL || cols == NULL) {
    return 0;
  }
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL) {
    return 0;
  }
  trilinos_klu_l_numeric *numeric = solver->numeric;
  const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  double *udiag = (double *)numeric->Udiag;
  UF_long **l_indices = solver->refactor_l_indices;
  int32_t **l_indices32 = solver->refactor_l_indices32;
  double **l_values = solver->refactor_l_values;
  UF_long **u_indices = solver->refactor_u_indices;
  double **u_values = solver->refactor_u_values;
  if (solver->refactor_lu_pointer_count != solver->n ||
      solver->refactor_col_block == NULL ||
      solver->refactor_col_ptr == NULL ||
      solver->refactor_block_start == NULL || numeric->Offp == NULL ||
      numeric->Offx == NULL || udiag == NULL || l_indices == NULL ||
      l_values == NULL || u_indices == NULL || u_values == NULL || nc < 2 ||
      nc > KLS_SNODE_MAX_FUSE || worker->pair_x == NULL ||
      (nc > 2 && worker->fuse_x == NULL)) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }

  const UF_long block = solver->refactor_col_block[cols[0]];
  if (block >= symbolic->nblocks) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }
  const UF_long k1 = symbolic->R[block];
  const UF_long k2 = symbolic->R[block + 1u];
  if (k2 <= k1 + 1u || k2 > solver->n) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }
  const UF_long *llen = numeric->Llen + k1;
  double *xv[KLS_SNODE_MAX_FUSE];
  xv[0] = worker->x;
  xv[1] = worker->pair_x;
  for (int c = 2; c < nc; ++c) {
    xv[c] = worker->fuse_x + (size_t)(c - 2) * (size_t)solver->n;
  }

  const UF_long *ui[KLS_SNODE_MAX_FUSE];
  const int32_t *ui32[KLS_SNODE_MAX_FUSE];
  double *ux[KLS_SNODE_MAX_FUSE];
  UF_long ulen[KLS_SNODE_MAX_FUSE];
  UF_long up[KLS_SNODE_MAX_FUSE];
  UF_long local_col[KLS_SNODE_MAX_FUSE];
  int done[KLS_SNODE_MAX_FUSE];
  UF_long producer_limit = k2 - k1;
  double *offx = (double *)numeric->Offx;
  for (int c = 0; c < nc; ++c) {
    const UF_long k = cols[c];
    if (k < k1 || k >= k2 || solver->refactor_col_block[k] != block) {
      kls_egraph_refactor_record_invalid(shared);
      return 0;
    }
    local_col[c] = k - k1;
    if (local_col[c] < producer_limit) {
      producer_limit = local_col[c];
    }
    UF_long poff = numeric->Offp[k];
    if (!kls_egraph_copy_unscaled_offblock_input(
          solver, shared->values, offx, &poff, numeric->Offp[k + 1u],
          solver->refactor_col_ptr[k], solver->refactor_block_start[k])) {
      kls_egraph_refactor_record_invalid(shared);
      return 0;
    }

    kls_egraph_scatter_btf_unscaled_input(
      solver, shared->values, xv[c], k1, solver->refactor_block_start[k],
      solver->refactor_col_ptr[k + 1u]);
    ui[c] = u_indices[k];
    ui32[c] = solver->refactor_u_indices32 != NULL
      ? solver->refactor_u_indices32[k] : NULL;
    ux[c] = u_values[k];
    ulen[c] = numeric->Ulen[k];
    up[c] = 0u;
    done[c] = 0;
  }

  int remaining = nc;
  unsigned spin = 0u;
  while (remaining > 0) {
    if (kls_egraph_refactor_should_stop(shared)) {
      return 0;
    }
    int progress = 0;
    int all_alive = 1;
    for (int c = 0; c < nc; ++c) {
      if (done[c] || up[c] >= ulen[c]) {
        all_alive = 0;
        break;
      }
    }
    if (all_alive) {
      const UF_long j0 = ui32[0] != NULL
        ? (UF_long)ui32[0][up[0]] : ui[0][up[0]];
      int aligned = 1;
      for (int c = 1; c < nc; ++c) {
        const UF_long jc = ui32[c] != NULL
          ? (UF_long)ui32[c][up[c]] : ui[c][up[c]];
        if (jc != j0) {
          aligned = 0;
          break;
        }
      }
      if (aligned) {
        const UF_long consumed = kls_snode_batch_consume_cached_multi(
          l_indices, l_values, llen, nc, ui, ui32, ux, ulen, up, xv, k1,
          producer_limit, shared->snode_run_end, shared);
        if (consumed != 0u) {
          for (int c = 0; c < nc; ++c) {
            up[c] += consumed;
          }
          progress = 1;
        } else if (kls_egraph_refactor_dependency_done_now(shared,
                                                            k1 + j0)) {
          /* Scalar aligned producers are worth fusing too: read each L row
             and value once, then update every independent SPA. */
          double u[KLS_SNODE_MAX_FUSE];
          int any_nonzero = 0;
          for (int c = 0; c < nc; ++c) {
            u[c] = xv[c][j0];
            xv[c][j0] = 0.0;
            ux[c][up[c]] = u[c];
            any_nonzero |= u[c] != 0.0;
          }
          if (any_nonzero) {
            const UF_long dep_global = k1 + j0;
            const UF_long length = llen[j0];
            const int32_t *rows32 =
              l_indices32 != NULL ? l_indices32[dep_global] : NULL;
            const double *values = l_values[dep_global];
            if (rows32 != NULL || length == 0u) {
              for (UF_long p = 0u; p < length; ++p) {
                const UF_long row = (UF_long)rows32[p];
                const double value = values[p];
                for (int c = 0; c < nc; ++c) {
                  xv[c][row] -= u[c] * value;
                }
              }
            } else {
              const UF_long *rows = l_indices[dep_global];
              for (UF_long p = 0u; p < length; ++p) {
                const UF_long row = rows[p];
                const double value = values[p];
                for (int c = 0; c < nc; ++c) {
                  xv[c][row] -= u[c] * value;
                }
              }
            }
          }
          for (int c = 0; c < nc; ++c) {
            up[c]++;
          }
          progress = 1;
        }
      }
    }

    /* A mismatched deck still contains independent ready sides.  Advance
       all of them before polling the shared frontier again; each cursor
       remains ordered and the done bitmap still gates every dependency. */
    for (int c = 0; c < nc; ++c) {
      if (done[c] || up[c] >= ulen[c]) {
        continue;
      }
      const UF_long consumed = kls_snode_batch_consume_cached(
        l_indices, l_values, llen, ui[c], ui32[c], ux[c], ulen[c], up[c],
        xv[c], k1, local_col[c], shared->snode_run_end, shared, 1);
      if (consumed != 0u) {
        up[c] += consumed;
        progress = 1;
        continue;
      }
      const UF_long j = ui32[c] != NULL
        ? (UF_long)ui32[c][up[c]] : ui[c][up[c]];
      if (kls_egraph_refactor_dependency_done_now(shared, k1 + j)) {
        kls_egraph_refactor_apply_btf_scalar_dep(
          shared, k1, j, up[c], ux[c], l_indices, l_indices32, l_values,
          llen, xv[c]);
        up[c]++;
        progress = 1;
      }
    }

    for (int c = 0; c < nc; ++c) {
      if (!done[c] && up[c] >= ulen[c]) {
        const UF_long k = cols[c];
        const UF_long local_k = local_col[c];
        const double pivot = xv[c][local_k];
        xv[c][local_k] = 0.0;
        if (pivot == 0.0) {
          kls_egraph_refactor_record_singular(shared, k, symbolic->Q[k]);
          if (solver->common.halt_if_singular) {
            return 0;
          }
        }
        kls_egraph_store_udiag(shared, udiag, k, pivot);
        kls_egraph_store_l_column_from_workspace(
          solver, xv[c], k, l_indices[k], l_values[k], llen[local_k], pivot);
        kls_egraph_refactor_mark_done(shared, k);
        done[c] = 1;
        remaining--;
        progress = 1;
      }
    }
    if (!progress) {
      kls_cpu_relax();
      kls_egraph_pipeline_pause(&spin);
    }
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
    if (!kls_egraph_copy_refactor_values(
          solver, shared->values, shared->rs, shared->scale, offx, &poff,
          poff_end, solver->refactor_col_ptr[k],
          solver->refactor_block_start[k])) {
      kls_egraph_refactor_record_invalid(shared);
      return 0;
    }
    for (UF_long p = solver->refactor_block_start[k];
         p < solver->refactor_col_ptr[k + 1u]; ++p) {
      if (kls_refactor_map_row_at(solver, p) != k) {
        kls_egraph_refactor_record_invalid(shared);
        return 0;
      }
      if (!kls_egraph_refactor_value(
            shared, kls_refactor_map_input_at(solver, p),
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
    if (!kls_egraph_copy_refactor_values(
          solver, shared->values, shared->rs, shared->scale, offx, &poff,
          poff_end, solver->refactor_col_ptr[k],
          solver->refactor_block_start[k])) {
      kls_egraph_refactor_record_invalid(shared);
      return 0;
    }
    if (!kls_egraph_scatter_refactor_input(
          solver, shared->values, shared->rs, shared->scale, x, k1, k2,
          solver->refactor_block_start[k],
          solver->refactor_col_ptr[k + 1u])) {
      kls_egraph_refactor_record_invalid(shared);
      return 0;
    }
  } else {
    if (!kls_egraph_scatter_refactor_input(
          solver, shared->values, shared->rs, shared->scale, x, 0u, solver->n,
          solver->refactor_col_ptr[k], solver->refactor_col_ptr[k + 1u])) {
      kls_egraph_refactor_record_invalid(shared);
      return 0;
    }
  }

  UF_long *ui = u_indices[k];
  const int32_t *ui32 =
    solver->refactor_u_indices32 != NULL
      ? solver->refactor_u_indices32[k] : NULL;
  double *ux = u_values[k];
  UF_long ucol_len = ulen[local_k];
  UF_long up = 0u;

  const int snode_batches_allowed =
    shared->snode_run_end != NULL;

  while (up < ucol_len) {

    if (snode_batches_allowed) {
      const UF_long consumed = kls_snode_batch_consume_cached(
        l_indices, l_values, llen, ui, ui32, ux, ucol_len, up, x, k1,
        local_k, shared->snode_run_end, shared, wait_for_dependencies);
      if (consumed != 0u) {
        up += consumed;
        continue;
      }
    }
    const UF_long j = ui32 != NULL ? (UF_long)ui32[up] : ui[up];
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

    if (ujk != 0.0) {
      UF_long *li = l_indices[k1 + j];
      double *lx = l_values[k1 + j];
      UF_long lcol_len = llen[j];
      kls_scatter_subtract_refactor_l(solver, x, k1 + j, li, lx, lcol_len,
                                      ujk);
    }
    up++;
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
  kls_egraph_store_l_column_from_workspace(solver, x, k, li, lx, lcol_len,
                                           ukk);

  return 1;
}

static int kls_egraph_refactor_dispatch_column(
  kls_egraph_refactor_worker *worker,
  UF_long k,
  int wait_for_dependencies) {
  if (worker == NULL || worker->shared == NULL) {
    return 0;
  }
  {
    const kls_solver *solver = worker->shared->solver;
    if (solver != NULL && solver->dense_tail_cols > 0 &&
        solver->symbolic != NULL) {
      const UF_long tail_end =
        solver->symbolic->R[solver->dense_tail_block + 1];
      if (k >= tail_end - solver->dense_tail_cols && k < tail_end) {
        /* pipe-emitted dense-tail column: the BLAS3 tail refresh
           rebuilds its values after this refactorization; completing
           it here releases dependents (only later tail columns, which
           no-op the same way) without the 24s scalar scatter walk */
        return 1;
      }
    }
  }
  switch (worker->shared->kernel) {
    case KLS_EGRAPH_REFACTOR_KERNEL_SINGLE_UNSCALED:
      if (worker->shared->single_unscaled_plain) {
        return kls_egraph_refactor_single_unscaled_plain_column(
          worker, k, wait_for_dependencies);
      }
      return kls_egraph_refactor_single_unscaled_column(
        worker, k, wait_for_dependencies);
    case KLS_EGRAPH_REFACTOR_KERNEL_SINGLE_SCALED:
      return kls_egraph_refactor_single_scaled_column(
        worker, k, wait_for_dependencies);
    case KLS_EGRAPH_REFACTOR_KERNEL_BTF_UNSCALED_DEFER_TERMINAL:
      if (k == worker->shared->deferred_pts_terminal_col) return 1;
      /* fall through */
    case KLS_EGRAPH_REFACTOR_KERNEL_BTF_UNSCALED:
      if (worker->shared->btf_unscaled_plain &&
          (!wait_for_dependencies ||
           !worker->shared->solver->padded_active)) {
        return kls_egraph_refactor_btf_unscaled_plain_cluster_column(worker,
          k, wait_for_dependencies);
      }
      return kls_egraph_refactor_btf_unscaled_column(
        worker, k, wait_for_dependencies);
    case KLS_EGRAPH_REFACTOR_KERNEL_GENERIC:
    default:
      return kls_egraph_refactor_column(worker, k, wait_for_dependencies);
  }
}

/* CKTSO tail repair can omit independent nodes outside the ETree-descendant
   tail.  Make those omitted columns current before the pivoted tail copies
   them as preserved suffix/gap columns. */
static int kls_refresh_pivot_tail_preserved_block_columns(
  kls_solver *solver,
  double *numeric_values,
  UF_long block,
  UF_long k1,
  UF_long nk,
  const unsigned char *tail_mask) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      numeric_values == NULL || tail_mask == NULL ||
      solver->symbolic->R == NULL ||
      block >= solver->symbolic->nblocks ||
      k1 != solver->symbolic->R[block] ||
      k1 + nk != solver->symbolic->R[block + 1u] ||
      k1 > solver->n || nk > solver->n - k1 ||
      (solver->common.scale > 0 && solver->numeric->Rs == NULL) ||
      (solver->common.scale <= 0 && solver->numeric->Rs != NULL) ||
      solver->numeric->Xwork == NULL ||
      solver->numeric->Ulen == NULL ||
      !kls_build_refactor_map(solver) ||
      !kls_build_refactor_lu_pointer_cache(solver)) {
    return 0;
  }
  if (solver->refactor_lu_pointer_count != solver->n ||
      solver->refactor_u_indices == NULL ||
      (solver->symbolic->nblocks > 1u &&
       (solver->refactor_block_start == NULL ||
        solver->refactor_col_block == NULL ||
        solver->numeric->Offp == NULL ||
        solver->numeric->Offx == NULL))) {
    return 0;
  }

  UF_long refresh_columns = 0;
  for (UF_long k = 0; k < nk; ++k) {
    if (tail_mask[k]) {
      continue;
    }
    const UF_long global_k = k1 + k;
    const UF_long *ui = solver->refactor_u_indices[global_k];
    const UF_long ulen = solver->numeric->Ulen[global_k];
    if (ui == NULL && ulen > 0u) {
      return 0;
    }
    for (UF_long p = 0; p < ulen; ++p) {
      const UF_long dep = ui[p];
      if (dep >= nk || tail_mask[dep]) {
        return 0;
      }
    }
    refresh_columns++;
  }
  if (refresh_columns == 0u) {
    return 1;
  }

  const int scaled = solver->common.scale > 0;
  kls_egraph_refactor_shared shared;
  memset(&shared, 0, sizeof(shared));
  if (pthread_mutex_init(&shared.lock, NULL) != 0) {
    return 0;
  }
  shared.solver = solver;
  shared.values = numeric_values;
  shared.rs = scaled ? solver->numeric->Rs : NULL;
  shared.check_pivots = 0;
  shared.scale = scaled ? (int)solver->common.scale : 0;
  shared.kernel = kls_egraph_refactor_kernel_for(solver, shared.scale);
  shared.thread_count = 1;
  shared.rejected_pivot = KLS_KLU_EMPTY;
  shared.rejected_pivot_col = KLS_KLU_EMPTY;
  shared.rejected_row = KLS_KLU_EMPTY;
  shared.rejected_multiplier_abs = -1.0;
  shared.rejected_pivot_abs = -1.0;
  shared.rejected_candidate_abs = -1.0;
  shared.numerical_rank = UF_long_max;
  shared.singular_col = KLS_KLU_EMPTY;
  shared.snode_run_end = kls_refactor_snode_run_end(solver);
  atomic_init(&shared.stop, 0);

  kls_egraph_refactor_worker worker;
  memset(&worker, 0, sizeof(worker));
  worker.shared = &shared;
  worker.x = (double *)solver->numeric->Xwork;

  int ok = 1;
  for (UF_long k = 0; k < nk; ++k) {
    if (tail_mask[k]) {
      continue;
    }
    memset(worker.x, 0, (size_t)nk * sizeof(*worker.x));
    if (!kls_egraph_refactor_dispatch_column(&worker, k1 + k, 0) ||
        shared.invalid || shared.pivot_rejected || shared.singular) {
      ok = 0;
      break;
    }
  }
  memset(worker.x, 0, (size_t)nk * sizeof(*worker.x));
  free(worker.segment_panel);
  free(worker.supernode_workspace);
  free(worker.index_workspace);
  free(worker.byte_workspace);
  free(worker.object_workspace);
  (void)pthread_mutex_destroy(&shared.lock);
  return ok;
}

static UF_long kls_restore_pivot_tail_plan(
  kls_solver *solver,
  UF_long block,
  UF_long k1,
  UF_long nk,
  UF_long local_reject,
  int refresh_state) {
  if (solver == NULL) {
    return nk;
  }
  const int saved_state = solver->fast_reject_refresh_state;
  solver->fast_reject_refresh_state = refresh_state;
  const UF_long *tail_cols = NULL;
  UF_long tail_count = 0;
  double tail_work = 0.0;
  (void)kls_build_fast_reject_pivoting_tail_plan(
    solver, block, k1, k1 + nk, local_reject, &tail_cols, &tail_count,
    &tail_work);
  solver->fast_reject_refresh_state = saved_state;
  return nk;
}

static UF_long kls_prepare_root_pivot_tail_independent_refresh(
  kls_solver *solver,
  double *numeric_values,
  UF_long block,
  UF_long k1,
  UF_long nk,
  UF_long local_reject) {
  if (solver == NULL || solver->symbolic == NULL || numeric_values == NULL ||
      local_reject != 0u || nk == 0u ||
      block >= solver->symbolic->nblocks ||
      solver->symbolic->R == NULL ||
      k1 != solver->symbolic->R[block] ||
      k1 + nk != solver->symbolic->R[block + 1u] ||
      solver->fast_reject_refresh_state != KLS_FAST_REJECT_REFRESH_PREFIX) {
    return nk;
  }

  const int saved_state = solver->fast_reject_refresh_state;
  solver->fast_reject_refresh_state = KLS_FAST_REJECT_REFRESH_UNKNOWN;
  const UF_long *tail_cols = NULL;
  UF_long tail_count = 0;
  double tail_work = 0.0;
  if (!kls_build_fast_reject_pivoting_tail_plan(
        solver, block, k1, k1 + nk, local_reject, &tail_cols, &tail_count,
        &tail_work)) {
    solver->fast_reject_refresh_state = saved_state;
    return nk;
  }
  solver->fast_reject_refresh_state = saved_state;

  const UF_long restart_end =
    kls_preferred_pivot_tail_restart_end(solver, k1, nk, local_reject);
  if (restart_end >= nk) {
    return nk;
  }
  const size_t nk_size = (size_t)nk;
  if ((UF_long)nk_size != nk) {
    return kls_restore_pivot_tail_plan(
      solver, block, k1, nk, local_reject, saved_state);
  }
  unsigned char *tail_mask =
    (unsigned char *)malloc(nk_size * sizeof(*tail_mask));
  if (tail_mask == NULL) {
    return kls_restore_pivot_tail_plan(
      solver, block, k1, nk, local_reject, saved_state);
  }
  const int mask_ok =
    kls_build_pivot_tail_restart_mask(
      solver, k1, nk, local_reject, restart_end, tail_mask);
  const int refresh_ok =
    mask_ok &&
    kls_refresh_pivot_tail_preserved_block_columns(
      solver, numeric_values, block, k1, nk, tail_mask);
  free(tail_mask);
  if (!refresh_ok) {
    return kls_restore_pivot_tail_plan(
      solver, block, k1, nk, local_reject, saved_state);
  }
  return restart_end;
}

static void kls_egraph_refactor_refresh_missing_prefix(
  kls_egraph_refactor_shared *shared,
  double *x,
  UF_long scratch_size) {
  if (shared == NULL || shared->solver == NULL ||
      shared->pipeline_done == NULL || x == NULL || scratch_size == 0u ||
      shared->pipeline_generation == 0u ||
      shared->rejected_pivot == KLS_KLU_EMPTY) {
    return;
  }
  kls_solver *solver = shared->solver;
  if (shared->rejected_pivot > solver->n) {
    return;
  }

  kls_egraph_refactor_worker worker;
  memset(&worker, 0, sizeof(worker));
  worker.shared = shared;
  worker.x = x;

  UF_long refreshed = 0;
  UF_long rejected_pivot = shared->rejected_pivot;
  for (UF_long col = 0; col < rejected_pivot; ++col) {
    if (atomic_load_explicit(&shared->pipeline_done[col],
                             memory_order_acquire) ==
        shared->pipeline_generation) {
      continue;
    }
    memset(x, 0, (size_t)scratch_size * sizeof(*x));
    if (!kls_egraph_refactor_dispatch_column(&worker, col, 0)) {
      break;
    }
    kls_egraph_refactor_mark_done(shared, col);
    refreshed++;
    if (shared->rejected_pivot != KLS_KLU_EMPTY &&
        shared->rejected_pivot < rejected_pivot) {
      rejected_pivot = shared->rejected_pivot;
    }
  }
  memset(x, 0, (size_t)scratch_size * sizeof(*x));
  free(worker.segment_panel);
  free(worker.supernode_workspace);
  free(worker.index_workspace);
  free(worker.byte_workspace);
  free(worker.object_workspace);

  if (refreshed > 0u) {
    solver->fast_rejected_prefix_refresh_columns += refreshed;
    solver->fast_rejected_prefix_refresh_count++;
    solver->stats.fast_rejected_prefix_refresh_columns =
      (int64_t)solver->fast_rejected_prefix_refresh_columns;
    solver->stats.fast_rejected_prefix_refresh_count =
      (int64_t)solver->fast_rejected_prefix_refresh_count;
  }
}

static int kls_egraph_separator_private_requested(const kls_solver *solver) {
  return solver != NULL && solver->symbolic != NULL &&
    /* A separator plan is an executor capability, not an ordering verdict.
       A dominant block exposes its private domains directly even when the
       factor arrived through a one-block, predicted, or preselected
       ordering.  Let the existing timed ordinary/private trial arbitrate
       that equivalent schedule; fragmented factors retain the ordinary
       executor. */
    solver->symbolic->maxblock >= 3u * solver->n / 4u &&
    solver->separator.n == solver->n &&
    solver->separator.thread_count >= 2u &&
    solver->separator.private_component_count >=
      solver->separator.thread_count &&
    solver->separator.private_rows <= solver->n &&
    solver->separator.pipeline_rows <= solver->n &&
    solver->separator.private_rows + solver->separator.pipeline_rows ==
      solver->n;
}

/* The separator forest proves that a private-domain dispatch is legal, not
   that it is faster than the ordinary dependency-level partition for this
   realized factor.  Compare those equivalent executors on the caller's real
   update stream.  The checked first refactor warms the ordinary path; the
   private arm then receives one unscored first-touch pass before a like-for-
   like sample from each arm. */
static int kls_egraph_separator_private_dispatch(kls_solver *solver,
                                                  int check_pivots) {
  if (solver == NULL) {
    return 0;
  }
  if (!kls_egraph_separator_private_requested(solver)) {
    solver->eg_separator_choice = -1;
    solver->eg_separator_pending = 0;
    return 0;
  }
  if (kls_repeated_update_workload(&solver->options) &&
      solver->refactor_separator_private_thread_ptr != NULL &&
      solver->refactor_separator_private_thread_count >= 2 &&
      solver->refactor_separator_private_column_count >=
        (UF_long)(2 * solver->refactor_separator_private_thread_count) &&
      solver->refactor_separator_cluster_tail_column_count == 0u &&
      solver->refactor_separator_private_unsafe_component_count == 0u &&
      solver->refactor_separator_private_component_count >=
        (UF_long)(2 * solver->refactor_separator_private_thread_count)) {
    const int threads = solver->refactor_separator_private_thread_count;
    UF_long min_columns = UF_long_max;
    UF_long max_columns = 0u;
    for (int tid = 0; tid < threads; ++tid) {
      const UF_long columns =
        solver->refactor_separator_private_thread_ptr[tid + 1] -
        solver->refactor_separator_private_thread_ptr[tid];
      if (columns < min_columns) min_columns = columns;
      if (columns > max_columns) max_columns = columns;
    }
    if (min_columns > 0u &&
        max_columns - min_columns <= min_columns) {
      /* Every scheduled cluster column belongs to a certified private
         component, the complete caller team owns nonempty and balanced
         slices, and no dependency tail remains to serialize publication.
         The ordinary level executor can expose no additional concurrency;
         settle the private layout before charging a complete alternative
         numeric update.  Partial, unsafe, or imbalanced forests retain the
         adjacent timing tournament below.  This is a realized schedule
         capability, independent of matrix dimensions or identity. */
      solver->eg_separator_choice = 1;
      solver->eg_separator_pending = 0;
      return 1;
    }
  }
  if (kls_repeated_update_workload(&solver->options) &&
      solver->refactor_separator_private_thread_ptr != NULL &&
      solver->refactor_separator_private_thread_count >= 4 &&
      solver->refactor_separator_private_component_count >=
        (UF_long)(2 * solver->refactor_separator_private_thread_count) &&
      solver->refactor_separator_private_unsafe_component_count <=
        solver->refactor_separator_private_component_count / 16u &&
      solver->refactor_separator_cluster_tail_column_count <=
        solver->refactor_separator_private_column_count / 2u) {
    const int threads = solver->refactor_separator_private_thread_count;
    const UF_long average =
      solver->refactor_separator_private_column_count / (UF_long)threads;
    const UF_long substantial_floor = average / 4u;
    int substantial_threads = 0;
    UF_long min_substantial = UF_long_max;
    UF_long max_substantial = 0u;
    for (int tid = 0; tid < threads; ++tid) {
      const UF_long columns =
        solver->refactor_separator_private_thread_ptr[tid + 1] -
        solver->refactor_separator_private_thread_ptr[tid];
      if (columns < substantial_floor) {
        continue;
      }
      substantial_threads++;
      if (columns < min_substantial) min_substantial = columns;
      if (columns > max_substantial) max_substantial = columns;
    }
    if (substantial_threads >= (3 * threads + 3) / 4 &&
        min_substantial > 0u &&
        max_substantial - min_substantial <= min_substantial) {
      /* A small number of threads may own only tiny separator leaves, while
         the remaining team owns balanced, dependency-closed domains covering
         at least two thirds of the cluster prefix.  Level scheduling cannot
         add useful concurrency inside those whole private components; it
         merely reintroduces a barrier at every dependency level.  Settle the
         certified private prefix and retain the ordinary parallel level walk
         for the bounded tail after its single join.  Coverage, balance, and
         unsafe-component limits are realized schedule properties rather than
         matrix identity or dimension thresholds. */
      solver->eg_separator_choice = 1;
      solver->eg_separator_pending = 0;
      return 1;
    }
  }
  if (check_pivots) {
    solver->eg_separator_pending = 0;
    return 0;
  }
  if (solver->eg_separator_choice != 0) {
    solver->eg_separator_pending = 0;
    return solver->eg_separator_choice > 0;
  }

  int side;
  if (solver->eg_separator_samples[1] == 0) {
    side = 1; /* unscored private-domain warm-up */
  } else if (solver->eg_separator_samples[0] == 0) {
    side = 0; /* ordinary level slices */
  } else {
    side = 1; /* scored private-domain sample */
  }
  solver->eg_separator_pending = side + 1;
  return side == 1;
}

static void kls_egraph_pipeline_pause(unsigned *spin) {
  if (spin == NULL) {
    return;
  }
  if (((*spin)++ & 4095u) == 0u) {
    struct timespec pause_time;
    pause_time.tv_sec = 0;
    pause_time.tv_nsec = 1000;
    (void)nanosleep(&pause_time, NULL);
  }
}

/* Fused dispatch supports the plain unscaled single-block and BTF kernels;
   every representation-changing or checked mode falls back to singles. */
static int kls_egraph_pair_dispatch_allowed(
  const kls_egraph_refactor_shared *shared) {
  /* Adopted by the per-matrix timed probe or armed for the probe refactor. */
  return shared->solver != NULL &&
         (shared->solver->eg_pair_choice > 0 ||
          shared->solver->eg_pair_pending) &&
         (shared->kernel == KLS_EGRAPH_REFACTOR_KERNEL_SINGLE_UNSCALED ||
          kls_egraph_kernel_is_btf_unscaled(shared->kernel)) &&
         !shared->check_pivots &&
         shared->snode_run_end != NULL;
}

/* Route a leased/popped column pair to the kernel-matched pair walk.
   Returns 1 on success, 0 on failure (caller records invalid), and -1
   when this pair cannot be fused (caller dispatches both singly). */
static int kls_egraph_refactor_dispatch_pair(
  kls_egraph_refactor_worker *worker,
  UF_long col,
  UF_long col2) {
  kls_egraph_refactor_shared *shared = worker->shared;
  if (shared->kernel == KLS_EGRAPH_REFACTOR_KERNEL_SINGLE_UNSCALED) {
    return kls_egraph_refactor_single_unscaled_column_pair(worker, col,
                                                           col2);
  }
  if (kls_egraph_kernel_is_btf_unscaled(shared->kernel)) {
    const kls_solver *solver = shared->solver;
    const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
    if (solver->refactor_col_block == NULL ||
        solver->refactor_col_block[col] !=
          solver->refactor_col_block[col2]) {
      return -1;
    }
    const UF_long block = solver->refactor_col_block[col];
    if (symbolic->R[block + 1u] - symbolic->R[block] <= 1u) {
      return -1;
    }
    return kls_egraph_refactor_btf_unscaled_column_pair(worker, col, col2);
  }
  return -1;
}

static KLS_ALWAYS_INLINE void kls_egraph_cluster_barrier_wait(
  kls_egraph_refactor_worker *worker) {
  kls_egraph_refactor_shared *shared = worker != NULL ? worker->shared : NULL;
  if (shared != NULL) {
    if (shared->use_spin_cluster_barrier) {
      const unsigned int generation = atomic_load_explicit(
        &shared->cluster_barrier_generation, memory_order_acquire);
      const unsigned int target = generation + 1u;
      atomic_store_explicit(&worker->cluster_barrier_arrival, target,
                            memory_order_release);
      if (worker->tid == 0) {
        unsigned int spins = 0u;
        for (int tid = 1; tid < shared->thread_count; ++tid) {
          while (atomic_load_explicit(
                   &worker->pool->workers[tid].cluster_barrier_arrival,
                   memory_order_acquire) != target) {
            kls_cpu_relax();
            if (++spins == 4096u) {
              sched_yield();
              spins = 0u;
            }
          }
        }
        atomic_store_explicit(&shared->cluster_barrier_generation,
                              target, memory_order_release);
        return;
      }
      unsigned int spins = 0u;
      while (atomic_load_explicit(&shared->cluster_barrier_generation,
                                  memory_order_acquire) == generation) {
        kls_cpu_relax();
        if (++spins == 4096u) {
          sched_yield();
          spins = 0u;
        }
      }
      return;
    }
    (void)pthread_barrier_wait(&shared->barrier);
  }
}

static UF_long kls_egraph_scale_row_boundary(
  const UF_long *row_ptr,
  UF_long n,
  UF_long nnz,
  int part,
  int parts) {
  if (part <= 0) {
    return 0u;
  }
  if (part >= parts) {
    return n;
  }
  const UF_long total = n + nnz;
  const UF_long target =
    (total / (UF_long)parts) * (UF_long)part +
    ((total % (UF_long)parts) * (UF_long)part) / (UF_long)parts;
  UF_long low = 0u;
  UF_long high = n;
  while (low < high) {
    const UF_long middle = low + (high - low) / 2u;
    if (row_ptr[middle] + middle < target) {
      low = middle + 1u;
    } else {
      high = middle;
    }
  }
  return low;
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

  if (shared->parallel_scale_rows) {
    const UF_long *restrict row_ptr = solver->refactor_scale_row_ptr;
    const int32_t *restrict input_pos =
      solver->refactor_scale_input_pos32;
    const double *restrict values = shared->parallel_scale_input_values;
    double *restrict rs = (double *)shared->rs;
    if (row_ptr == NULL || input_pos == NULL || values == NULL || rs == NULL) {
      if (worker->tid == 0) {
        kls_egraph_refactor_record_invalid(shared);
      }
      kls_egraph_cluster_barrier_wait(worker);
      return;
    }
    const UF_long row_begin = kls_egraph_scale_row_boundary(
      row_ptr, solver->n, solver->nnz, worker->tid,
      shared->thread_count);
    const UF_long row_end = kls_egraph_scale_row_boundary(
      row_ptr, solver->n, solver->nnz, worker->tid + 1,
      shared->thread_count);
    if (shared->parallel_scale_method == 1) {
      for (UF_long row = row_begin; row < row_end; ++row) {
        double sum = 0.0;
        for (UF_long p = row_ptr[row]; p < row_ptr[row + 1u]; ++p) {
          sum += fabs(values[(UF_long)input_pos[p]]);
        }
        rs[row] = sum != 0.0 ? sum : 1.0;
      }
    } else {
      for (UF_long row = row_begin; row < row_end; ++row) {
        double maximum = 0.0;
        for (UF_long p = row_ptr[row]; p < row_ptr[row + 1u]; ++p) {
          maximum = fmax(maximum,
                         fabs(values[(UF_long)input_pos[p]]));
        }
        rs[row] = maximum != 0.0 ? maximum : 1.0;
      }
    }
    kls_egraph_cluster_barrier_wait(worker);
  }

  UF_long cluster_levels = shared->cluster_level_count;
  if (cluster_levels > solver->refactor_level_count) {
    cluster_levels = solver->refactor_level_count;
  }
  /* The common CKTSO-style cluster phase has neither speculative column
     claims nor BTF run-group wakeups.  Its level barrier already establishes
     uniqueness and ordering, so avoid three general-pipeline consultations
     per column (claim lookup, load-before-publish, and run-group lookup).
     A release publication is still required for the later pipeline levels. */
  atomic_uint *const direct_cluster_done =
    (shared->pipeline_done != NULL)
      ? shared->pipeline_done : NULL;
  const unsigned int direct_cluster_generation = shared->pipeline_generation;
  const int use_alpha4_separator =
    shared->separator_private &&
    cluster_levels == solver->refactor_cluster_level_count_alpha4 &&
    solver->refactor_separator_private_cols_alpha4 != NULL &&
    solver->refactor_separator_private_thread_ptr_alpha4 != NULL;
  const UF_long *separator_private_cols = use_alpha4_separator
    ? solver->refactor_separator_private_cols_alpha4
    : solver->refactor_separator_private_cols;
  const UF_long *separator_private_thread_ptr = use_alpha4_separator
    ? solver->refactor_separator_private_thread_ptr_alpha4
    : solver->refactor_separator_private_thread_ptr;
  const UF_long *separator_cluster_tail_cols = use_alpha4_separator
    ? solver->refactor_separator_cluster_tail_cols_alpha4
    : solver->refactor_separator_cluster_tail_cols;
  const UF_long *separator_cluster_tail_level_ptr = use_alpha4_separator
    ? solver->refactor_separator_cluster_tail_level_ptr_alpha4
    : solver->refactor_separator_cluster_tail_level_ptr;
  const UF_long *separator_cluster_tail_level_thread_ptr =
    use_alpha4_separator
      ? solver->refactor_separator_cluster_tail_level_thread_ptr_alpha4
      : solver->refactor_separator_cluster_tail_level_thread_ptr;
  const UF_long separator_private_column_count = use_alpha4_separator
    ? solver->refactor_separator_private_column_count_alpha4
    : solver->refactor_separator_private_column_count;
  const UF_long separator_cluster_tail_column_count = use_alpha4_separator
    ? solver->refactor_separator_cluster_tail_column_count_alpha4
    : solver->refactor_separator_cluster_tail_column_count;
  const int use_separator_private =
    shared->separator_private &&
    direct_cluster_done != NULL &&
    separator_private_cols != NULL &&
    separator_private_thread_ptr != NULL &&
    solver->refactor_separator_private_thread_count == shared->thread_count &&
    (use_alpha4_separator ||
     solver->refactor_separator_private_cluster_level_count == cluster_levels) &&
    (separator_cluster_tail_column_count == 0u ||
     separator_cluster_tail_cols != NULL);

  if (use_separator_private) {
    const UF_long begin =
      separator_private_thread_ptr[worker->tid];
    const UF_long end =
      separator_private_thread_ptr[worker->tid + 1];
    if (begin <= end &&
        end <= separator_private_column_count &&
        !kls_egraph_refactor_should_stop(shared)) {
      for (UF_long pos = begin; pos < end; ++pos) {
        const UF_long col = separator_private_cols[pos];
        if (col >= solver->n ||
            !kls_egraph_refactor_dispatch_column(worker, col, 0)) {
          if (!kls_egraph_refactor_should_stop(shared)) {
            kls_egraph_refactor_record_invalid(shared);
          }
          break;
        }
        if (!shared->cluster_done_premarked) {
          atomic_store_explicit(&direct_cluster_done[col],
                                direct_cluster_generation,
                                memory_order_release);
        }
      }
    } else if (begin > end ||
               end > separator_private_column_count) {
      kls_egraph_refactor_record_invalid(shared);
    }

    /* Separator/rejected columns can depend on any completed private domain.
       Private queues have joined above.  The retained tail is grouped by its
       exact numeric dependency level, so a level's columns can run in
       parallel while the level barriers preserve every dependency.  Older
       schedules without the level index retain the conservative serial
       natural-order fallback. */
    kls_egraph_cluster_barrier_wait(worker);
    const UF_long *tail_level_ptr =
      separator_cluster_tail_level_ptr;
    const UF_long *tail_level_thread_ptr =
      separator_cluster_tail_level_thread_ptr;
    if (tail_level_ptr != NULL) {
      for (UF_long level = 0u; level < cluster_levels; ++level) {
        const UF_long begin = tail_level_ptr[level];
        const UF_long end = tail_level_ptr[level + 1u];
        if (begin > end ||
            end > separator_cluster_tail_column_count) {
          kls_egraph_refactor_record_invalid(shared);
        } else if (!kls_egraph_refactor_should_stop(shared)) {
          UF_long local_begin = begin + (UF_long)worker->tid;
          UF_long local_end = end;
          if (tail_level_thread_ptr != NULL) {
            const UF_long *parts = tail_level_thread_ptr +
              level * (UF_long)(shared->thread_count + 1);
            local_begin = parts[worker->tid];
            local_end = parts[worker->tid + 1];
            if (local_begin < begin || local_begin > local_end ||
                local_end > end) {
              kls_egraph_refactor_record_invalid(shared);
              local_begin = local_end;
            }
          }
          const UF_long stride = tail_level_thread_ptr != NULL
            ? 1u : (UF_long)shared->thread_count;
          for (UF_long pos = local_begin; pos < local_end; pos += stride) {
            const UF_long col = separator_cluster_tail_cols[pos];
            if (col >= solver->n ||
                !kls_egraph_refactor_dispatch_column(worker, col, 0)) {
              if (!kls_egraph_refactor_should_stop(shared)) {
                kls_egraph_refactor_record_invalid(shared);
              }
              break;
            }
            if (!shared->cluster_done_premarked) {
              atomic_store_explicit(&direct_cluster_done[col],
                                    direct_cluster_generation,
                                    memory_order_release);
            }
          }
        }
        kls_egraph_cluster_barrier_wait(worker);
      }
    } else {
      if (worker->tid == 0 && !kls_egraph_refactor_should_stop(shared)) {
        for (UF_long pos = 0u;
             pos < separator_cluster_tail_column_count;
             ++pos) {
          const UF_long col = separator_cluster_tail_cols[pos];
          if (col >= solver->n ||
              !kls_egraph_refactor_dispatch_column(worker, col, 0)) {
            if (!kls_egraph_refactor_should_stop(shared)) {
              kls_egraph_refactor_record_invalid(shared);
            }
            break;
          }
          if (!shared->cluster_done_premarked) {
            atomic_store_explicit(&direct_cluster_done[col],
                                  direct_cluster_generation,
                                  memory_order_release);
          }
        }
      }
      kls_egraph_cluster_barrier_wait(worker);
    }
  } else {
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
            if (direct_cluster_done != NULL) {
              if (!shared->cluster_done_premarked) {
                atomic_store_explicit(&direct_cluster_done[col],
                                      direct_cluster_generation,
                                      memory_order_release);
              }
            } else {
              kls_egraph_refactor_mark_done(shared, col);

            }
          }
        } else {
          for (UF_long pos = begin + (UF_long)worker->tid;
               pos < end; pos += (UF_long)shared->thread_count) {
            const UF_long col = solver->refactor_level_cols[pos];
            if (!kls_egraph_refactor_dispatch_column(worker, col, 0)) {
              break;
            }
            if (direct_cluster_done != NULL) {
              if (!shared->cluster_done_premarked) {
                atomic_store_explicit(&direct_cluster_done[col],
                                      direct_cluster_generation,
                                      memory_order_release);
              }
            } else {
              kls_egraph_refactor_mark_done(shared, col);

            }
          }
        }
      }

      kls_egraph_cluster_barrier_wait(worker);
    }
  }

  if (shared->pipeline_done != NULL &&
      cluster_levels < solver->refactor_level_count &&
      !kls_egraph_refactor_should_stop(shared)) {
    if (solver->refactor_level_cols != NULL) {
      const int pair_allowed = kls_egraph_pair_dispatch_allowed(shared);
      if (pair_allowed && worker->pair_x == NULL) {
        worker->pair_x =
          (double *)calloc((size_t)solver->n, sizeof(double));
      }
      const int quad_allowed = pair_allowed &&
        (shared->kernel == KLS_EGRAPH_REFACTOR_KERNEL_SINGLE_UNSCALED ||
         kls_egraph_kernel_is_btf_unscaled(shared->kernel)) &&
        (solver->eg_pair_choice == 2 || solver->eg_pair_pending == 2);
      if (quad_allowed && worker->fuse_x == NULL) {
        worker->fuse_x = (double *)calloc(
          (size_t)solver->n * (KLS_SNODE_MAX_FUSE - 2u), sizeof(double));
      }
      const int try_pair = pair_allowed && worker->pair_x != NULL;
      const unsigned long fuse_width =
        quad_allowed && worker->fuse_x != NULL ? KLS_SNODE_MAX_FUSE
        : try_pair ? 2ul : 1ul;
      for (;;) {
        if (kls_egraph_refactor_should_stop(shared)) {
          break;
        }

        if (fuse_width > 2ul) {
          /* Single-block factors retain two separately claimed pairs:
             concurrent workers' interleaving decorrelates dependencies
             (four adjacent positions serialized rajat25).  BTF multi-fusion
             instead needs columns from one block and benefits from common
             producer streams, so claim one adjacent four-column deck.  The
             existing timed quad probe remains the authority on whether that
             equivalent schedule is retained for a particular numeric. */
          const int adjacent_btf_deck =
            kls_egraph_kernel_is_btf_unscaled(shared->kernel);
          const UF_long qpos = (UF_long)atomic_fetch_add_explicit(
            &shared->next_pipeline_pos, adjacent_btf_deck ? 4ul : 2ul,
            memory_order_relaxed);
          if (qpos >= shared->pipeline_pos_end) {
            break;
          }
          const UF_long qpos2 = adjacent_btf_deck
            ? qpos + 2u
            : (UF_long)atomic_fetch_add_explicit(
                &shared->next_pipeline_pos, 2ul, memory_order_relaxed);
          UF_long qp[KLS_SNODE_MAX_FUSE];
          UF_long nqp = 0;
          for (UF_long w = 0; w < 2u; ++w) {
            if (qpos + w < shared->pipeline_pos_end) {
              qp[nqp++] = qpos + w;
            }
          }
          for (UF_long w = 0; w < 2u; ++w) {
            if (qpos2 + w < shared->pipeline_pos_end) {
              qp[nqp++] = qpos2 + w;
            }
          }
          UF_long got[KLS_SNODE_MAX_FUSE];
          int ng = 0;
          for (UF_long w = 0; w < nqp; ++w) {
            got[ng++] = solver->refactor_level_cols[qp[w]];
          }
          int ok2 = 1;
          if (ng >= 2) {
            if (shared->kernel == KLS_EGRAPH_REFACTOR_KERNEL_SINGLE_UNSCALED) {
              ok2 = kls_egraph_refactor_single_unscaled_column_multi(
                worker, got, ng);
            } else {
              const UF_long block = solver->refactor_col_block[got[0]];
              int same_nontrivial_block = block < solver->symbolic->nblocks &&
                solver->symbolic->R[block + 1u] -
                    solver->symbolic->R[block] > 1u;
              for (int g = 1; g < ng && same_nontrivial_block; ++g) {
                same_nontrivial_block =
                  solver->refactor_col_block[got[g]] == block;
              }
              if (same_nontrivial_block) {
                ok2 = kls_egraph_refactor_btf_unscaled_column_multi(
                  worker, got, ng);
              } else {
                /* A deck can straddle a BTF boundary.  Preserve ownership
                   and correctness with ordinary dispatch; later decks in a
                   common block remain eligible for measured fusion. */
                for (int g = 0; g < ng && ok2; ++g) {
                  ok2 = kls_egraph_refactor_dispatch_column(
                    worker, got[g], 1);
                  if (ok2) {
                    kls_egraph_refactor_mark_done(shared, got[g]);
                  }
                }
              }
            }
          } else if (ng == 1) {
            ok2 = kls_egraph_refactor_dispatch_column(worker, got[0], 1);
            if (ok2) {
              kls_egraph_refactor_mark_done(shared, got[0]);
            }
          }
          if (!ok2) {
            if (!kls_egraph_refactor_should_stop(shared)) {
              kls_egraph_refactor_record_invalid(shared);
            }
            break;
          }

          continue;
        }
        const UF_long pos =
          (UF_long)atomic_fetch_add_explicit(&shared->next_pipeline_pos,
                                             try_pair ? 2ul : 1ul,
                                             memory_order_relaxed);
        if (pos >= shared->pipeline_pos_end) {
          break;
        }
        const UF_long col = solver->refactor_level_cols[pos];
        if (try_pair && pos + 1u < shared->pipeline_pos_end) {
          const UF_long col2 = solver->refactor_level_cols[pos + 1u];
          const int pair_status =
            kls_egraph_refactor_dispatch_pair(worker, col, col2);
          if (pair_status == 0) {
            if (!kls_egraph_refactor_should_stop(shared)) {
              kls_egraph_refactor_record_invalid(shared);
            }
            break;
          }
          if (pair_status < 0) {
            if (!kls_egraph_refactor_dispatch_column(worker, col, 1)) {
              if (!kls_egraph_refactor_should_stop(shared)) {
                kls_egraph_refactor_record_invalid(shared);
              }
              break;
            }
            kls_egraph_refactor_mark_done(shared, col);
            if (!kls_egraph_refactor_dispatch_column(worker, col2, 1)) {
              if (!kls_egraph_refactor_should_stop(shared)) {
                kls_egraph_refactor_record_invalid(shared);
              }
              break;
            }
            kls_egraph_refactor_mark_done(shared, col2);
          }

          continue;
        }
        if (!kls_egraph_refactor_dispatch_column(worker, col, 1)) {
          if (!kls_egraph_refactor_should_stop(shared)) {
            kls_egraph_refactor_record_invalid(shared);
          }
          break;
        }
        kls_egraph_refactor_mark_done(shared, col);

      }
    } else {
      kls_egraph_refactor_record_invalid(shared);
    }
  }
}

static KLS_ALWAYS_INLINE double kls_row_solve_view_value(
  const kls_solver *solver,
  const kls_row_solve_factor_view *view,
  UF_long p) {
  const UF_long source_pos =
    view != NULL && view->source_pos != NULL ? view->source_pos[p] : p;
  return view != NULL && view->source_upper
    ? kls_row_refactor_u_solve_value(solver, source_pos)
    : kls_row_refactor_l_solve_value(solver, source_pos);
}

static int kls_row_solve_view_has_values(
  const kls_solver *solver,
  const kls_row_solve_factor_view *view) {
  if (solver == NULL || view == NULL || view->ptr == NULL) {
    return 0;
  }
  const UF_long nnz = view->ptr[solver->n];
  if (nnz == 0u) {
    return 1;
  }
  return view->cols != NULL &&
         (view->values != NULL ||
          kls_row_refactor_solve_has_values(solver, view->source_upper));
}

static void kls_row_solve_view_row(const kls_solver *solver,
                                   const kls_row_solve_factor_view *view,
                                   double *work,
  UF_long row) {
  double value = work[row];
  if (view->values != NULL) {
    for (UF_long p = view->ptr[row]; p < view->ptr[row + 1u]; ++p) {
      value -= view->values[p] * work[view->cols[p]];
    }
  } else {
    for (UF_long p = view->ptr[row]; p < view->ptr[row + 1u]; ++p) {
      value -= kls_row_solve_view_value(solver, view, p) *
               work[view->cols[p]];
    }
  }
  if (view->diagonal) {
    const double *udiag = (const double *)solver->numeric->Udiag;
    value /= udiag[row];
  }
  work[row] = value;
}

static void kls_row_solve_sparse_worker_run(
  kls_egraph_refactor_worker *worker) {
  if (worker == NULL || worker->shared == NULL) {
    return;
  }
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  const kls_row_solve_factor_view *view = shared->row_solve_view;
  if (solver == NULL || view == NULL || shared->row_solve_work == NULL ||
      shared->thread_count <= 0 ||
      view->thread_count != shared->thread_count) {
    return;
  }

  const UF_long *level_ptr = view->sparse_level_ptr;
  const UF_long *level_rows = view->sparse_level_rows;
  const UF_long level_count = view->sparse_level_count;
  UF_long cluster_levels = view->sparse_cluster_levels;
  const UF_long sparse_end = view->dense_tail_start;
  if (level_ptr == NULL || level_rows == NULL || level_count == 0u ||
      cluster_levels == 0u || sparse_end > solver->n ||
      !kls_row_solve_view_has_values(solver, view) ||
      (view->diagonal && solver->numeric->Udiag == NULL)) {
    return;
  }
  if (cluster_levels > level_count) {
    cluster_levels = level_count;
  }

  double *work = shared->row_solve_work;
  for (UF_long level = 0u; level < cluster_levels; ++level) {
    const UF_long begin = level_ptr[level];
    const UF_long end = level_ptr[level + 1u];
    for (UF_long pos = begin + (UF_long)worker->tid;
         pos < end;
         pos += (UF_long)shared->thread_count) {
      const UF_long row = level_rows[pos];
      kls_row_solve_view_row(solver, view, work, row);
    }
    (void)pthread_barrier_wait(&shared->barrier);
  }

  if (worker->tid == 0) {
    for (UF_long level = cluster_levels; level < level_count; ++level) {
      for (UF_long pos = level_ptr[level]; pos < level_ptr[level + 1u]; ++pos) {
        const UF_long row = level_rows[pos];
        kls_row_solve_view_row(solver, view, work, row);
      }
    }
  }
}


static void kls_compact_amf_two_block_residual_worker_run(
  kls_egraph_refactor_worker *worker) {
  if (worker == NULL || worker->shared == NULL) {
    return;
  }
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  if (solver == NULL || shared->values == NULL || shared->rs == NULL ||
      shared->row_solve_residual_x == NULL ||
      shared->row_solve_work == NULL || shared->thread_count <= 0 ||
      worker->tid < 0 || worker->tid >= shared->thread_count ||
      solver->solve_refine_csr_ptr16 == NULL ||
      solver->solve_refine_csr_col_pos32 == NULL) {
    return;
  }
  const UF_long begin =
    (UF_long)solver->solve_refine_csr_row_bound16[worker->tid];
  const UF_long end =
    (UF_long)solver->solve_refine_csr_row_bound16[worker->tid + 1];
  const uint16_t *restrict ptr = solver->solve_refine_csr_ptr16;
  const uint32_t *restrict col_pos = solver->solve_refine_csr_col_pos32;
  const double *restrict a = shared->values;
  const double *restrict b = shared->rs;
  const double *restrict x = shared->row_solve_residual_x;
  double *restrict residual = shared->row_solve_work;
  double *restrict results = shared->contract_rgrowth_results;
  double bmax = 0.0;
  double bnorm2 = 0.0;
  double rmax = 0.0;
  double rnorm2 = 0.0;
  for (UF_long row = begin; row < end; ++row) {
    double value = b[row];
    for (UF_long p = (UF_long)ptr[row];
         p < (UF_long)ptr[row + 1u]; ++p) {
      const uint32_t packed = col_pos[p];
      value = fma(-a[packed >> 16u],
                  x[packed & UINT32_C(0xffff)], value);
    }
    residual[row] = value;
    if (results != NULL) {
      const double ab = fabs(b[row]);
      const double ar = fabs(value);
      bmax = bmax < ab ? ab : bmax;
      bnorm2 += b[row] * b[row];
      rmax = rmax < ar ? ar : rmax;
      rnorm2 += value * value;
    }
  }
  if (results != NULL) {
    results[4 * worker->tid] = bmax;
    results[4 * worker->tid + 1] = bnorm2;
    results[4 * worker->tid + 2] = rmax;
    results[4 * worker->tid + 3] = rnorm2;
  }
}

/* Bound the rounding error of the compact CSR residual evaluation without
   recomputing the residual in extended precision.  For a row with d fused
   multiply-adds, gamma_d * (|b| + sum |a*x|) bounds the dot-product error;
   divide by (1-gamma_d) to also cover the rounded absolute-sum accumulator.
   This pass is dispatched only for a raw residual already inside the public
   contract but outside the ordinary conservative fast margin. */
static void kls_compact_residual_error_bound_worker_run(
  kls_egraph_refactor_worker *worker) {
  if (worker == NULL || worker->shared == NULL) {
    return;
  }
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  double *results = shared->contract_rgrowth_results;
  if (solver == NULL || shared->values == NULL || shared->rs == NULL ||
      shared->row_solve_residual_x == NULL || results == NULL ||
      shared->thread_count != 5 || worker->tid < 0 || worker->tid >= 5 ||
      solver->solve_refine_csr_ptr16 == NULL ||
      solver->solve_refine_csr_col_pos32 == NULL) {
    return;
  }
  const UF_long begin =
    (UF_long)solver->solve_refine_csr_row_bound16[worker->tid];
  const UF_long end =
    (UF_long)solver->solve_refine_csr_row_bound16[worker->tid + 1];
  const uint16_t *restrict ptr = solver->solve_refine_csr_ptr16;
  const uint32_t *restrict col_pos = solver->solve_refine_csr_col_pos32;
  const double *restrict a = shared->values;
  const double *restrict b = shared->rs;
  const double *restrict x = shared->row_solve_residual_x;
  double error_norm2 = 0.0;
  for (UF_long row = begin; row < end; ++row) {
    double magnitude_sum = fabs(b[row]);
    const UF_long row_begin = (UF_long)ptr[row];
    const UF_long row_end = (UF_long)ptr[row + 1u];
    for (UF_long p = row_begin; p < row_end; ++p) {
      const uint32_t packed = col_pos[p];
      magnitude_sum += fabs(a[packed >> 16u]) *
        fabs(x[packed & UINT32_C(0xffff)]);
    }
    const double operations = (double)(row_end - row_begin + 1u);
    const double gamma = operations * DBL_EPSILON;
    if (!(gamma < 0.5) || !isfinite(magnitude_sum)) {
      error_norm2 = HUGE_VAL;
      break;
    }
    const double error = (gamma / (1.0 - 2.0 * gamma)) * magnitude_sum;
    error_norm2 += error * error;
  }
  results[worker->tid] = error_norm2;
}

static void kls_csr32_residual_worker_run(
  kls_egraph_refactor_worker *worker) {
  if (worker == NULL || worker->shared == NULL) {
    return;
  }
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  if (solver == NULL || shared->values == NULL || shared->rs == NULL ||
      shared->row_solve_residual_x == NULL ||
      shared->row_solve_work == NULL || shared->thread_count < 2 ||
      shared->thread_count > 8 ||
      worker->tid < 0 || worker->tid >= shared->thread_count ||
      solver->solve_refine_csr_ptr32 == NULL ||
      solver->solve_refine_csr_pos32 == NULL ||
      (solver->solve_refine_csr_col16 == NULL &&
       solver->solve_refine_csr_col32 == NULL) ||
      solver->solve_refine_csr32_threads != shared->thread_count) {
    return;
  }
  const UF_long begin =
    (UF_long)solver->solve_refine_csr_row_bound32[worker->tid];
  const UF_long end =
    (UF_long)solver->solve_refine_csr_row_bound32[worker->tid + 1];
  const uint32_t *restrict ptr = solver->solve_refine_csr_ptr32;
  const uint32_t *restrict pos = solver->solve_refine_csr_pos32;
  const uint16_t *restrict cols16 = solver->solve_refine_csr_col16;
  const uint32_t *restrict cols32 = solver->solve_refine_csr_col32;
  const double *restrict a = shared->values;
  const double *restrict b = shared->rs;
  const double *restrict x = shared->row_solve_residual_x;
  double *restrict residual = shared->row_solve_work;
  double *restrict results = shared->contract_rgrowth_results;
  double bmax = 0.0;
  double bnorm2 = 0.0;
  double rmax = 0.0;
  double rnorm2 = 0.0;
  const UF_long *restrict rinv = solver->solve_refine_rinv;
  const UF_long *restrict cmap = solver->user_col_perm;
  const double *restrict rs_inv = solver->row_scale != NULL
    ? solver->solve_refine_rs_inv : NULL;
  const double *restrict cs = solver->col_scale;
  const int transformed = rinv != NULL || cmap != NULL ||
    rs_inv != NULL || cs != NULL;
  /* Share source while retaining compile-time column widths in both loops. */
#define KLS_TRANSFORMED_CSR_ROWS(columns) \
    for (UF_long row = begin; row < end; ++row) { \
      const UF_long public_row = rinv != NULL ? rinv[row] : row; \
      const double row_recip = rs_inv != NULL ? rs_inv[row] : 1.0; \
      double value = b[public_row]; \
      for (UF_long p = (UF_long)ptr[row]; \
           p < (UF_long)ptr[row + 1u]; ++p) { \
        const UF_long col = (UF_long)(columns)[p]; \
        double xv = x[cmap != NULL ? cmap[col] : col]; \
        if (cs != NULL) { \
          xv /= cs[col]; \
        } \
        value = fma(-a[pos[p]] * row_recip, xv, value); \
      } \
      residual[public_row] = value; \
      if (results != NULL) { \
        const double ab = fabs(b[public_row]); \
        const double ar = fabs(value); \
        bmax = bmax < ab ? ab : bmax; \
        bnorm2 += b[public_row] * b[public_row]; \
        rmax = rmax < ar ? ar : rmax; \
        rnorm2 += value * value; \
      } \
    }
  if (cols16 != NULL) {
    KLS_TRANSFORMED_CSR_ROWS(cols16)
  } else if (transformed) {
    KLS_TRANSFORMED_CSR_ROWS(cols32)
#undef KLS_TRANSFORMED_CSR_ROWS
  } else {
    for (UF_long row = begin; row < end; ++row) {
      double value = b[row];
      for (UF_long p = (UF_long)ptr[row];
           p < (UF_long)ptr[row + 1u]; ++p) {
        value = fma(-a[pos[p]], x[(UF_long)cols32[p]], value);
      }
      residual[row] = value;
      if (results != NULL) {
        const double ab = fabs(b[row]);
        const double ar = fabs(value);
        bmax = bmax < ab ? ab : bmax;
        bnorm2 += b[row] * b[row];
        rmax = rmax < ar ? ar : rmax;
        rnorm2 += value * value;
      }
    }
  }
  if (results != NULL) {
    results[4 * worker->tid] = bmax;
    results[4 * worker->tid + 1] = bnorm2;
    results[4 * worker->tid + 2] = rmax;
    results[4 * worker->tid + 3] = rnorm2;
  }
}

/* Form max and squared-L2 statistics for the RHS and residual in one
   contiguous vector pass.  This is intentionally independent of the sparse
   pattern: it only divides uniform O(n) work across an already-retained
   numeric worker pool. */
static void kls_plain_contract_vector_stats_worker_run(
  kls_egraph_refactor_worker *worker) {
  if (worker == NULL || worker->shared == NULL) {
    return;
  }
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  const int tid = worker->tid;
  const int threads = shared->thread_count;
  double *results = shared->contract_rgrowth_results;
  if (solver == NULL || shared->values == NULL || shared->rs == NULL ||
      results == NULL || threads < 2 || tid < 0 || tid >= threads) {
    return;
  }
  const UF_long begin =
    (solver->n * (UF_long)tid) / (UF_long)threads;
  const UF_long end =
    (solver->n * (UF_long)(tid + 1)) / (UF_long)threads;
  const double *restrict b = shared->values;
  const double *restrict residual = shared->rs;
  double bmax = 0.0;
  double bnorm2 = 0.0;
  double rmax = 0.0;
  double rnorm2 = 0.0;
  for (UF_long i = begin; i < end; ++i) {
    const double bv = b[i];
    const double rv = residual[i];
    const double ab = fabs(bv);
    const double ar = fabs(rv);
    bmax = bmax < ab ? ab : bmax;
    rmax = rmax < ar ? ar : rmax;
    bnorm2 += bv * bv;
    rnorm2 += rv * rv;
  }
  results[4 * tid] = bmax;
  results[4 * tid + 1] = bnorm2;
  results[4 * tid + 2] = rmax;
  results[4 * tid + 3] = rnorm2;
}

static void kls_row_solve_worker_run(kls_egraph_refactor_worker *worker) {
  if (worker == NULL || worker->shared == NULL) {
    return;
  }
  kls_egraph_refactor_shared *shared = worker->shared;
  if (shared->row_solve_mode == 4) {
    kls_compact_amf_two_block_residual_worker_run(worker);
    return;
  }
  if (shared->row_solve_mode == 7) {
    kls_compact_residual_error_bound_worker_run(worker);
    return;
  }
  if (shared->row_solve_mode == 5) {
    kls_csr32_residual_worker_run(worker);
    return;
  }
  if (shared->row_solve_mode == 6) {
    kls_plain_contract_vector_stats_worker_run(worker);
    return;
  }
  if (shared->row_solve_mode == 2) {
    kls_row_solve_sparse_worker_run(worker);
    return;
  }
  kls_solver *solver = shared->solver;
  const kls_row_solve_factor_view *view = shared->row_solve_view;
  if (solver == NULL || view == NULL || shared->row_solve_work == NULL ||
      shared->thread_count <= 0) {
    return;
  }

  const int upper = view->upper;
  const UF_long *bounds = view->slice_bounds;
  const UF_long *ptr = view->ptr;
  const UF_long *cols = view->cols;
  const UF_long *split = view->segment_split;
  const UF_long *thread_bounds = view->thread_bounds;
  if (bounds == NULL || ptr == NULL || cols == NULL ||
      !kls_row_solve_view_has_values(solver, view) ||
      split == NULL || thread_bounds == NULL ||
      view->thread_count != shared->thread_count ||
      (view->diagonal && solver->numeric->Udiag == NULL)) {
    return;
  }

  double *work = shared->row_solve_work;
  const double *values = view->values;

  for (UF_long pass = 0u; pass < KLS_ROW_SOLVE_TRAPEZOID_SLICES; ++pass) {
    const UF_long slice = upper
      ? KLS_ROW_SOLVE_TRAPEZOID_SLICES - 1u - pass
      : pass;
    const UF_long begin = bounds[slice];
    const UF_long end = bounds[slice + 1u];
    const UF_long *parts =
      thread_bounds + slice * ((UF_long)shared->thread_count + 1u);
    const UF_long worker_begin = parts[worker->tid];
    const UF_long worker_end = parts[worker->tid + 1];

    for (UF_long row = worker_begin; row < worker_end; ++row) {
      const UF_long row_begin = ptr[row];
      const UF_long row_end = ptr[row + 1u];
      const UF_long cut = split[row];
      const UF_long p_begin = upper ? cut : row_begin;
      const UF_long p_end = upper ? row_end : cut;
      double value = work[row];
      if (values != NULL) {
        for (UF_long p = p_begin; p < p_end; ++p) {
          value -= values[p] * work[cols[p]];
        }
      } else {
        for (UF_long p = p_begin; p < p_end; ++p) {
          value -= kls_row_solve_view_value(solver, view, p) *
                   work[cols[p]];
        }
      }
      work[row] = value;
    }

    (void)pthread_barrier_wait(&shared->barrier);

    if (worker->tid == 0) {
      if (upper) {
        for (UF_long remaining = end; remaining > begin; --remaining) {
          const UF_long row = remaining - 1u;
          double value = work[row];
          const UF_long cut = split[row];
          if (values != NULL) {
            for (UF_long p = ptr[row]; p < cut; ++p) {
              value -= values[p] * work[cols[p]];
            }
          } else {
            for (UF_long p = ptr[row]; p < cut; ++p) {
              value -= kls_row_solve_view_value(solver, view, p) *
                       work[cols[p]];
            }
          }
          if (view->diagonal) {
            const double *udiag = (const double *)solver->numeric->Udiag;
            value /= udiag[row];
          }
          work[row] = value;
        }
      } else {
        for (UF_long row = begin; row < end; ++row) {
          double value = work[row];
          const UF_long cut = split[row];
          if (values != NULL) {
            for (UF_long p = cut; p < ptr[row + 1u]; ++p) {
              value -= values[p] * work[cols[p]];
            }
          } else {
            for (UF_long p = cut; p < ptr[row + 1u]; ++p) {
              value -= kls_row_solve_view_value(solver, view, p) *
                       work[cols[p]];
            }
          }
          if (view->diagonal) {
            const double *udiag = (const double *)solver->numeric->Udiag;
            value /= udiag[row];
          }
          work[row] = value;
        }
      }
    }

    (void)pthread_barrier_wait(&shared->barrier);
  }
}

static KLS_ALWAYS_INLINE void kls_egraph_refactor_pool_run_worker(
  kls_egraph_refactor_worker *worker) {
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_dense_help_self = worker;
  if (shared->solve_permute_mode) {
    const UF_long begin =
      (shared->solve_permute_count * (UF_long)worker->tid) /
      (UF_long)shared->thread_count;
    const UF_long end =
      (shared->solve_permute_count * (UF_long)(worker->tid + 1)) /
      (UF_long)shared->thread_count;
    const double *restrict input = shared->solve_permute_input;
    double *restrict output = shared->solve_permute_output;
    const uint32_t *restrict map = shared->solve_permute_map;
    if (shared->solve_permute_scatter) {
      for (UF_long p = begin; p < end; ++p) {
        output[map[p]] = input[p];
      }
    } else {
      for (UF_long p = begin; p < end; ++p) {
        output[p] = input[map[p]];
      }
    }
  } else if (shared->value_prep_mode) {
    const UF_long begin =
      (shared->value_prep_nnz * (UF_long)worker->tid) /
      (UF_long)shared->thread_count;
    const UF_long end =
      (shared->value_prep_nnz * (UF_long)(worker->tid + 1)) /
      (UF_long)shared->thread_count;
    const double *restrict input = shared->value_prep_input;
    double *restrict output = shared->value_prep_output;
    const UF_long *restrict input_pos =
      shared->value_prep_input_to_csc;
    for (UF_long p = begin; p < end; ++p) {
      output[p] = input[input_pos[p]];
    }
  } else if (shared->workspace_touch_mode) {
    double *restrict values = shared->workspace_touch_values;
    const size_t stride = (size_t)shared->thread_count * 512u;
    for (size_t p = (size_t)worker->tid * 512u;
         p < shared->workspace_touch_count; p += stride) {
      values[p] = 0.0;
    }
  } else if (shared->pts_refactor_mode) {
    kls_pts_refactor_pool_worker_run(worker);
  } else if (shared->pts_solve_mode) {
    kls_pts_pool_worker_run(worker);
  } else if (shared->lean_pattern_mode) {
    kls_lean_pattern_worker_run(worker);
  } else if (shared->lean_refactor_mode) {
    kls_lean_parallel_worker_run(worker);
  } else if (shared->row_publish_mode) {
    kls_row_publish_worker_run(worker);
  } else if (shared->row_solve_mode) {
    kls_row_solve_worker_run(worker);
  } else if (shared->row_refactor_mode) {
    kls_row_refactor_worker_run(worker);
  } else {
    kls_egraph_refactor_worker_run(worker);
    if (shared->parallel_scale_permute) {
      kls_solver *solver = shared->solver;
      trilinos_klu_l_numeric *numeric = solver->numeric;
      double *restrict rs = (double *)numeric->Rs;
      double *restrict permuted = solver->refactor_scale_permute_values;
      const UF_long *restrict pnum = numeric->Pnum;
      pthread_barrier_wait(&shared->barrier);
      if (!shared->invalid && !shared->pivot_rejected &&
          !(shared->singular && solver->common.halt_if_singular)) {
        const UF_long begin =
          (solver->n / (UF_long)shared->thread_count) *
            (UF_long)worker->tid +
          ((solver->n % (UF_long)shared->thread_count) *
            (UF_long)worker->tid) / (UF_long)shared->thread_count;
        const UF_long end =
          (solver->n / (UF_long)shared->thread_count) *
            (UF_long)(worker->tid + 1) +
          ((solver->n % (UF_long)shared->thread_count) *
            (UF_long)(worker->tid + 1)) / (UF_long)shared->thread_count;
        int invalid = rs == NULL || permuted == NULL || pnum == NULL;
        for (UF_long k = begin; !invalid && k < end; ++k) {
          const UF_long row = pnum[k];
          if (row >= solver->n) {
            invalid = 1;
            break;
          }
          permuted[k] = rs[row];
        }
        if (invalid) {
          kls_egraph_refactor_record_invalid(shared);
        }
        pthread_barrier_wait(&shared->barrier);
        if (!shared->invalid) {
          memcpy(rs + begin, permuted + begin,
                 (size_t)(end - begin) * sizeof(*rs));
        }
        /* No third rendezvous is needed: every worker next performs only
           disjoint optional maintenance, and the pool completion join keeps
           the caller from observing Rs until all slices have returned. */
      }
    }
  }
  if (shared->refine_copy_input != NULL &&
      shared->refine_copy_output != NULL) {
    const UF_long begin =
      (shared->refine_copy_nnz * (UF_long)worker->tid) /
      (UF_long)shared->thread_count;
    const UF_long end =
      (shared->refine_copy_nnz * (UF_long)(worker->tid + 1)) /
      (UF_long)shared->thread_count;
    if (end > begin) {
      memcpy(shared->refine_copy_output + begin,
             shared->refine_copy_input + begin,
             (size_t)(end - begin) * sizeof(*shared->refine_copy_output));
    }
  }
  kls_dense_help_self = NULL;
}

#ifdef __linux__
typedef struct kls_affinity_scope {
  cpu_set_t saved;
  int active;
} kls_affinity_scope;

static int kls_read_sysfs_integer(const char *path, int *value_out) {
  FILE *file = fopen(path, "r");
  if (file == NULL) {
    return 0;
  }
  int value = -1;
  const int ok = fscanf(file, "%d", &value) == 1 && value >= 0;
  fclose(file);
  if (ok) {
    *value_out = value;
  }
  return ok;
}

static int kls_read_sysfs_cache_size(const char *path,
                                     uint64_t *bytes_out) {
  FILE *file = fopen(path, "r");
  if (file == NULL) {
    return 0;
  }
  unsigned long long amount = 0u;
  char unit = '\0';
  const int fields = fscanf(file, "%llu%c", &amount, &unit);
  fclose(file);
  if (fields < 1) {
    return 0;
  }
  uint64_t multiplier = 1u;
  if (unit == 'K' || unit == 'k') {
    multiplier = UINT64_C(1024);
  } else if (unit == 'M' || unit == 'm') {
    multiplier = UINT64_C(1024) * UINT64_C(1024);
  }
  if ((uint64_t)amount > UINT64_MAX / multiplier) {
    return 0;
  }
  *bytes_out = (uint64_t)amount * multiplier;
  return 1;
}

typedef struct kls_topology_cache {
  int package[CPU_SETSIZE];
  int llc_id[CPU_SETSIZE];
  int core_id[CPU_SETSIZE];
  uint64_t llc_bytes[CPU_SETSIZE];
  uint64_t max_frequency_khz[CPU_SETSIZE];
  unsigned char valid[CPU_SETSIZE];
} kls_topology_cache;

static kls_topology_cache kls_cached_topology;
static pthread_once_t kls_topology_once = PTHREAD_ONCE_INIT;

static void kls_initialize_topology_cache(void) {
  memset(&kls_cached_topology, 0, sizeof(kls_cached_topology));
  long configured = sysconf(_SC_NPROCESSORS_CONF);
  if (configured < 1) {
    return;
  }
  if (configured > CPU_SETSIZE) {
    configured = CPU_SETSIZE;
  }
  for (int cpu = 0; cpu < (int)configured; ++cpu) {
    char path[160];
    snprintf(path, sizeof(path),
             "/sys/devices/system/cpu/cpu%d/topology/physical_package_id",
             cpu);
    if (!kls_read_sysfs_integer(path,
                                &kls_cached_topology.package[cpu])) {
      continue;
    }
    snprintf(path, sizeof(path),
             "/sys/devices/system/cpu/cpu%d/topology/core_id", cpu);
    if (!kls_read_sysfs_integer(path,
                                &kls_cached_topology.core_id[cpu])) {
      continue;
    }
    snprintf(path, sizeof(path),
             "/sys/devices/system/cpu/cpu%d/cache/index3/id", cpu);
    if (!kls_read_sysfs_integer(path,
                                &kls_cached_topology.llc_id[cpu])) {
      continue;
    }
    snprintf(path, sizeof(path),
             "/sys/devices/system/cpu/cpu%d/cache/index3/size", cpu);
    if (!kls_read_sysfs_cache_size(
          path, &kls_cached_topology.llc_bytes[cpu])) {
      continue;
    }
    snprintf(path, sizeof(path),
             "/sys/devices/system/cpu/cpu%d/cpufreq/cpuinfo_max_freq", cpu);
    int max_frequency_khz = 0;
    if (kls_read_sysfs_integer(path, &max_frequency_khz)) {
      kls_cached_topology.max_frequency_khz[cpu] =
        (uint64_t)max_frequency_khz;
    }
    kls_cached_topology.valid[cpu] = 1u;
  }
}

static void kls_warm_topology_cache(void) {
  (void)pthread_once(&kls_topology_once, kls_initialize_topology_cache);
}

/* Return the cache capacity of a caller-constrained execution domain.  libc's
   LEVEL3_CACHE_SIZE describes the CPU on which the query happens (and on
   asymmetric-cache systems may instead expose a machine-wide value), so it
   cannot distinguish a taskset confined to another LLC.  When every allowed
   CPU belongs to one known LLC, that domain is authoritative.  With a broad
   or unavailable affinity mask, retain the conservative libc fallback. */
static uint64_t kls_execution_domain_llc_bytes(void) {
  kls_warm_topology_cache();
  cpu_set_t allowed;
  if (sched_getaffinity(0, sizeof(allowed), &allowed) == 0) {
    int representative = -1;
    int constrained = 1;
    for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
      if (!CPU_ISSET(cpu, &allowed)) continue;
      if (!kls_cached_topology.valid[cpu]) {
        constrained = 0;
        break;
      }
      if (representative < 0) {
        representative = cpu;
      } else if (kls_cached_topology.package[cpu] !=
                   kls_cached_topology.package[representative] ||
                 kls_cached_topology.llc_id[cpu] !=
                   kls_cached_topology.llc_id[representative]) {
        constrained = 0;
        break;
      }
    }
    if (constrained && representative >= 0) {
      return kls_cached_topology.llc_bytes[representative];
    }
  }
#if defined(_SC_LEVEL3_CACHE_SIZE)
  const long llc = sysconf(_SC_LEVEL3_CACHE_SIZE);
  return llc > 0 ? (uint64_t)llc : 0u;
#else
  return 0u;
#endif
}

/* Select one logical CPU from each physical core in a single last-level
   cache domain.  This is a hardware/workload decision, not a matrix policy:
   repeated sparse factors exchange dependency state and repeatedly revisit
   their retained numeric, so crossing LLC domains is pure coherence and
   cache-capacity cost when one domain already has enough physical cores.
   Respect the caller's cpuset and decline when topology is unavailable or
   the requested crew does not fit without SMT. */
static int kls_compact_llc_affinity_plan(const kls_solver *solver,
                                         int thread_count,
                                         int *cpus_out,
                                         cpu_set_t *allowed_out) {
  if (solver == NULL || cpus_out == NULL || allowed_out == NULL ||
      thread_count < 2 || thread_count > CPU_SETSIZE) {
    return 0;
  }
  if (!kls_repeated_update_workload(&solver->options)) {
    return 0;
  }

  kls_warm_topology_cache();
  if (sched_getaffinity(0, sizeof(*allowed_out), allowed_out) != 0) {
    return 0;
  }
  const int *package = kls_cached_topology.package;
  const int *llc_id = kls_cached_topology.llc_id;
  const int *core_id = kls_cached_topology.core_id;
  const uint64_t *llc_bytes = kls_cached_topology.llc_bytes;
  const uint64_t *max_frequency_khz =
    kls_cached_topology.max_frequency_khz;
  const unsigned char *valid = kls_cached_topology.valid;
  int allowed_count = 0;
  for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
    if (!CPU_ISSET(cpu, allowed_out)) {
      continue;
    }
    allowed_count++;
  }
  if (allowed_count < thread_count) {
    return 0;
  }

  /* A shared cpuset bounds placement but does not assign one worker per
     core. Even within one LLC, retain the explicit plan: delayed workers
     can otherwise stall dependency consumers that spin awaiting them. */

  /* Linux exposes a portable per-CPU maximum policy frequency even when it
     does not export Intel's hybrid core-type CPUID leaf through sysfs.  On a
     hybrid package, the highest allowed tier is the P-core tier.  Prefer it
     only when an entire physical-core crew fits there; otherwise preserve the
     existing LLC selection and the caller's requested parallelism. */
  uint64_t performance_frequency = 0u;
  for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
    if (CPU_ISSET(cpu, allowed_out) && valid[cpu] &&
        max_frequency_khz[cpu] > performance_frequency) {
      performance_frequency = max_frequency_khz[cpu];
    }
  }

  const int current_cpu = sched_getcpu();
  int best_cpu = -1;
  int best_core_count = 0;
  uint64_t best_bytes = 0u;
  uint64_t smallest_eligible_bytes = UINT64_MAX;
  int eligible_group_count = 0;
  int best_contains_current = 0;
  int selected_performance_only = 0;
  for (int preference_pass = 0; preference_pass < 2; ++preference_pass) {
    const int performance_only =
      preference_pass == 0 && performance_frequency > 0u;
    best_cpu = -1;
    best_core_count = 0;
    best_bytes = 0u;
    smallest_eligible_bytes = UINT64_MAX;
    eligible_group_count = 0;
    best_contains_current = 0;
    for (int representative = 0; representative < CPU_SETSIZE;
         ++representative) {
    if (!valid[representative]) {
      continue;
    }
    int already_seen = 0;
    for (int prior = 0; prior < representative; ++prior) {
      if (valid[prior] && package[prior] == package[representative] &&
          llc_id[prior] == llc_id[representative]) {
        already_seen = 1;
        break;
      }
    }
    if (already_seen) {
      continue;
    }
    int unique_cores[CPU_SETSIZE];
    int core_count = 0;
    int contains_current = 0;
    for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
      if (!CPU_ISSET(cpu, allowed_out) || !valid[cpu] ||
          package[cpu] != package[representative] ||
          llc_id[cpu] != llc_id[representative] ||
          (performance_only &&
           max_frequency_khz[cpu] != performance_frequency)) {
        continue;
      }
      contains_current |= cpu == current_cpu;
      int core_seen = 0;
      for (int i = 0; i < core_count; ++i) {
        if (unique_cores[i] == core_id[cpu]) {
          core_seen = 1;
          break;
        }
      }
      if (!core_seen) {
        unique_cores[core_count++] = core_id[cpu];
      }
    }
    if (core_count < thread_count) {
      continue;
    }
    eligible_group_count++;
    if (llc_bytes[representative] < smallest_eligible_bytes) {
      smallest_eligible_bytes = llc_bytes[representative];
    }
    if (best_cpu < 0 || llc_bytes[representative] > best_bytes ||
        (llc_bytes[representative] == best_bytes &&
         contains_current > best_contains_current) ||
        (llc_bytes[representative] == best_bytes &&
         contains_current == best_contains_current &&
         core_count < best_core_count)) {
      best_cpu = representative;
      best_core_count = core_count;
      best_bytes = llc_bytes[representative];
      best_contains_current = contains_current;
    }
  }
    if (best_cpu >= 0 && smallest_eligible_bytes != UINT64_MAX &&
        (performance_only ||
         (eligible_group_count >= 2 && best_bytes > smallest_eligible_bytes))) {
      selected_performance_only = performance_only;
      break;
    }
    if (!performance_only) {
      break;
    }
  }
  if (best_cpu < 0 || smallest_eligible_bytes == UINT64_MAX ||
      (!selected_performance_only &&
       (eligible_group_count < 2 || best_bytes <= smallest_eligible_bytes))) {
    return 0;
  }
  /* A repeated sparse pool exchanges completion words and producer columns
     on every numeric update.  When the complete requested crew fits in one
     physical LLC, keep it there rather than letting the scheduler spread a
     latency-sensitive generation across asymmetric cache domains. */
  int selected_cores[CPU_SETSIZE];
  int selected = 0;
  if (current_cpu >= 0 && current_cpu < CPU_SETSIZE &&
      valid[current_cpu] && package[current_cpu] == package[best_cpu] &&
      llc_id[current_cpu] == llc_id[best_cpu] &&
      (!selected_performance_only ||
       max_frequency_khz[current_cpu] == performance_frequency)) {
    selected_cores[selected] = core_id[current_cpu];
    cpus_out[selected++] = current_cpu;
  }
  for (int cpu = 0; cpu < CPU_SETSIZE && selected < thread_count; ++cpu) {
    if (!CPU_ISSET(cpu, allowed_out) || !valid[cpu] ||
        package[cpu] != package[best_cpu] ||
        llc_id[cpu] != llc_id[best_cpu] ||
        (selected_performance_only &&
         max_frequency_khz[cpu] != performance_frequency)) {
      continue;
    }
    int core_seen = 0;
    for (int i = 0; i < selected; ++i) {
      if (selected_cores[i] == core_id[cpu]) {
        core_seen = 1;
        break;
      }
    }
    if (!core_seen) {
      selected_cores[selected] = core_id[cpu];
      cpus_out[selected++] = cpu;
    }
  }
  if (selected != thread_count) {
    return 0;
  }
  return 1;
}

static kls_affinity_scope kls_affinity_scope_enter(int cpu) {
  kls_affinity_scope scope;
  memset(&scope, 0, sizeof(scope));
  if (cpu < 0 || cpu >= CPU_SETSIZE ||
      sched_getaffinity(0, sizeof(scope.saved), &scope.saved) != 0) {
    return scope;
  }
  if (CPU_COUNT(&scope.saved) == 1 && CPU_ISSET(cpu, &scope.saved)) {
    return scope;
  }
  cpu_set_t selected;
  CPU_ZERO(&selected);
  CPU_SET(cpu, &selected);
  if (sched_setaffinity(0, sizeof(selected), &selected) == 0) {
    scope.active = 1;
  }
  return scope;
}

static void kls_affinity_scope_leave(kls_affinity_scope *scope) {
  if (scope != NULL && scope->active) {
    (void)sched_setaffinity(0, sizeof(scope->saved), &scope->saved);
    scope->active = 0;
  }
}
#else
typedef struct kls_affinity_scope { int active; } kls_affinity_scope;
static void kls_warm_topology_cache(void) {}
static kls_affinity_scope kls_affinity_scope_enter(int cpu) {
  (void)cpu;
  kls_affinity_scope scope = {0};
  return scope;
}
static void kls_affinity_scope_leave(kls_affinity_scope *scope) {
  (void)scope;
}
#endif

static void *kls_egraph_refactor_pool_worker_main(void *arg) {
  kls_egraph_refactor_worker *worker = (kls_egraph_refactor_worker *)arg;
  kls_egraph_refactor_pool *pool = worker->pool;
  kls_egraph_refactor_shared *shared = worker->shared;
  unsigned long seen_generation = 0;

  for (;;) {
    unsigned long generation;
    unsigned spin = 0;
    for (;;) {
      if (atomic_load_explicit(&pool->shutdown, memory_order_acquire)) {
        return NULL;
      }
      generation =
        atomic_load_explicit(&pool->generation, memory_order_acquire);
      if (generation != seen_generation) {
        break;
      }
      if (pool->busy_wait && spin < pool->worker_spin_iters) {
        spin++;
        kls_cpu_relax();
        continue;
      }
      /* Repeated SPICE refactors arrive back to back, so a short spin
         usually catches the next generation without a futex round trip;
         fall back to the condition variable when the solver goes idle. */
      pthread_mutex_lock(&shared->lock);
      while (!atomic_load_explicit(&pool->shutdown, memory_order_acquire) &&
             atomic_load_explicit(&pool->generation,
                                  memory_order_acquire) ==
               seen_generation) {
        pthread_cond_wait(&pool->work_cond, &shared->lock);
      }
      pthread_mutex_unlock(&shared->lock);
      spin = 0;
    }
    seen_generation = generation;

    kls_egraph_refactor_pool_run_worker(worker);

    if (atomic_load_explicit(&pool->spin_completion,
                             memory_order_acquire)) {
      atomic_store_explicit(&worker->completed_generation, generation,
                            memory_order_release);
    } else if (atomic_fetch_sub_explicit(&pool->active_workers, 1,
                                         memory_order_acq_rel) == 1) {
      pthread_mutex_lock(&shared->lock);
      pthread_cond_signal(&pool->done_cond);
      pthread_mutex_unlock(&shared->lock);
    }
  }
}

/* Publish a new pool generation, run worker zero on the calling thread, and
   wait for the N-1 background workers. This keeps an N-thread request at N
   runnable threads while making the caller useful instead of spinning or
   sleeping beside N workers. Expects shared->lock held on entry and returns
   with it held. */
static void kls_egraph_pool_dispatch_and_wait(
  kls_egraph_refactor_pool *pool,
  kls_egraph_refactor_shared *shared,
  int thread_count) {
  kls_affinity_scope affinity =
    kls_affinity_scope_enter(pool->caller_affinity_cpu);
  atomic_store_explicit(&pool->active_workers, thread_count - 1,
                        memory_order_relaxed);
  atomic_store_explicit(&pool->spin_completion, 0, memory_order_relaxed);
  atomic_fetch_add_explicit(&pool->generation, 1ul, memory_order_release);
  pthread_cond_broadcast(&pool->work_cond);
  pthread_mutex_unlock(&shared->lock);
  kls_egraph_refactor_pool_run_worker(&pool->workers[0]);
  pthread_mutex_lock(&shared->lock);
  while (atomic_load_explicit(&pool->active_workers,
                              memory_order_acquire) > 0) {
    pthread_cond_wait(&pool->done_cond, &shared->lock);
  }
  kls_affinity_scope_leave(&affinity);
}

static void kls_egraph_pool_dispatch_and_spin_wait(
  kls_egraph_refactor_pool *pool,
  kls_egraph_refactor_shared *shared,
  int thread_count) {
  kls_affinity_scope affinity =
    kls_affinity_scope_enter(pool->caller_affinity_cpu);
  atomic_store_explicit(&pool->active_workers, thread_count - 1,
                        memory_order_relaxed);
  atomic_store_explicit(&pool->spin_completion, 1, memory_order_relaxed);
  const unsigned long generation =
    atomic_fetch_add_explicit(&pool->generation, 1ul,
                              memory_order_release) + 1ul;
  pthread_cond_broadcast(&pool->work_cond);
  pthread_mutex_unlock(&shared->lock);
  kls_egraph_refactor_pool_run_worker(&pool->workers[0]);
  for (int tid = 1; tid < thread_count; ++tid) {
    while (atomic_load_explicit(
             &pool->workers[tid].completed_generation,
             memory_order_acquire) != generation) {
      kls_cpu_relax();
    }
  }
  atomic_store_explicit(&pool->active_workers, 0, memory_order_release);
  atomic_store_explicit(&pool->spin_completion, 0, memory_order_release);
  pthread_mutex_lock(&shared->lock);
  kls_affinity_scope_leave(&affinity);
}

static void destroy_egraph_refactor_pool(kls_solver *solver) {
  if (solver == NULL || solver->egraph_pool == NULL) {
    return;
  }
  kls_egraph_refactor_pool *pool = solver->egraph_pool;
  if (pool->lock_initialized) {
    pthread_mutex_lock(&pool->shared.lock);
    atomic_store_explicit(&pool->shutdown, 1, memory_order_release);
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
      free(pool->workers[i].help_x);
      free(pool->workers[i].pair_x);
      free(pool->workers[i].fuse_x);
      free(pool->workers[i].segment_panel);
      free(pool->workers[i].supernode_workspace);
      free(pool->workers[i].index_workspace);
      free(pool->workers[i].byte_workspace);
      free(pool->workers[i].object_workspace);
      free(pool->workers[i].row_target_stamp_workspace);
      free(pool->workers[i].row_target_pos_workspace);
      free(pool->workers[i].row_target_kind_workspace);
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
  atomic_init(&pool->generation, 0ul);
  atomic_init(&pool->active_workers, 0);
  atomic_init(&pool->spin_completion, 0);
  atomic_init(&pool->shutdown, 0);
  {
    const char *busy = getenv("KLS_DISABLE_EGRAPH_BUSY_WAIT");
    pool->busy_wait = !(busy != NULL && busy[0] == '1');
    /* These retained numerics return after every update/solve pair.  Keep
       their crew live through that bounded gap so the next numeric phase
       avoids seven condition-variable wakeups. */
    unsigned initial_spin = (kls_symmetric_partial_diagonal_match_factor_cycle(solver)
          ? 100000u : KLS_EGRAPH_WORKER_SPIN_ITERS);
    const char *spin_iters = NULL;
    {
      spin_iters = getenv("KLS_EGRAPH_WORKER_SPIN_OVERRIDE");
    }
    if (spin_iters != NULL) {
      const long requested = strtol(spin_iters, NULL, 10);
      if (requested >= 0 && requested <= 1000000) {
        initial_spin = (unsigned)requested;
      }
    }
    pool->worker_spin_iters = initial_spin;
  }
  atomic_init(&pool->shared.stop, 0);
  atomic_init(&pool->shared.cluster_barrier_generation, 0u);
  atomic_init(&pool->shared.next_pipeline_pos, 0ul);

  atomic_init(&pool->shared.btf_scalar_run_exec_runs, 0ul);
  atomic_init(&pool->shared.btf_scalar_run_exec_rows, 0ul);
  atomic_init(&pool->shared.btf_scalar_run_exec_entries, 0ul);
  atomic_init(&pool->shared.btf_scalar_run_exec_max_rows, 0ul);

  atomic_init(&pool->shared.row_pipeline_private_pos, 0ul);
  atomic_init(&pool->shared.row_pipeline_local_ready_groups, 0ul);
  atomic_init(&pool->shared.row_prefactor_rows, 0ul);
  atomic_init(&pool->shared.row_prefactor_deps, 0ul);
  atomic_init(&pool->shared.row_prefactor_supernode_rows, 0ul);
  atomic_init(&pool->shared.row_prefactor_supernode_deps, 0ul);
  atomic_init(&pool->shared.row_pipeline_ready_head, 0ul);
  atomic_init(&pool->shared.row_pipeline_ready_tail, 0ul);
  atomic_init(&pool->shared.row_pipeline_completed_groups, 0ul);
  pool->shared.row_pipeline_private_group_mask = NULL;
  pool->shared.row_pipeline_ordered_private = 0;
  pool->caller_affinity_cpu = -1;

#ifdef __linux__
  int affinity_cpus[CPU_SETSIZE];
  cpu_set_t affinity_allowed;
  CPU_ZERO(&affinity_allowed);
  int compact_affinity = kls_compact_llc_affinity_plan(
    solver, thread_count, affinity_cpus, &affinity_allowed);
  int affinity_applied = 0;
#else
  int compact_affinity = 0;
#endif

  for (int i = 0; i < thread_count; ++i) {
    pool->workers[i].shared = &pool->shared;
    pool->workers[i].pool = pool;
    atomic_init(&pool->workers[i].completed_generation, 0ul);
    atomic_init(&pool->workers[i].cluster_barrier_arrival, 0u);
    pool->workers[i].tid = i;
  }
  for (int i = 1; i < thread_count; ++i) {
    if (pthread_create(&pool->threads[i - 1], NULL,
                       kls_egraph_refactor_pool_worker_main,
                       &pool->workers[i]) != 0) {
      break;
    }
    pool->created_count++;
#ifdef __linux__
    if (compact_affinity) {
      cpu_set_t selected;
      CPU_ZERO(&selected);
      CPU_SET(affinity_cpus[i], &selected);
      if (pthread_setaffinity_np(pool->threads[i - 1], sizeof(selected),
                                 &selected) != 0) {
        compact_affinity = 0;
      } else {
        affinity_applied++;
      }
    }
#endif
  }
  if (pool->created_count != thread_count - 1) {
    solver->egraph_pool = pool;
    destroy_egraph_refactor_pool(solver);
    return NULL;
  }
#ifdef __linux__
  if (compact_affinity) {
    pool->caller_affinity_cpu = affinity_cpus[0];
  } else if (affinity_applied > 0 && CPU_COUNT(&affinity_allowed) > 0) {
    /* A partial placement is worse than the scheduler's coherent policy.
       Restore every worker to the caller's original cpuset if any affinity
       request failed. */
    for (int i = 1; i < thread_count; ++i) {
      (void)pthread_setaffinity_np(pool->threads[i - 1],
                                   sizeof(affinity_allowed),
                                   &affinity_allowed);
    }
  }
#else

#endif
  solver->egraph_pool = pool;
  return pool;
}

static int kls_parallel_prepare_permuted_values(
  kls_solver *solver,
  const double *values) {
  if (solver == NULL || values == NULL || solver->values == NULL ||
      solver->input_to_csc == NULL || solver->row_scale != NULL ||
      solver->col_scale != NULL || solver->options.threads < 2 ||
      solver->nnz < 200000u ||
      !kls_build_prepared_value_input_pos(solver)) {
    return 0;
  }

  /* Do not create a crew just for this bandwidth pass.  Once EGraph has
     established its persistent pool, reuse it on every later numeric. */
  kls_egraph_refactor_pool *pool = solver->egraph_pool;
  if (pool == NULL || pool->thread_count < 2 ||
      pool->created_count != pool->thread_count - 1) {
    return 0;
  }
  kls_egraph_refactor_shared *shared = &pool->shared;
  pthread_mutex_lock(&shared->lock);
  if (atomic_load_explicit(&pool->active_workers,
                           memory_order_acquire) != 0) {
    pthread_mutex_unlock(&shared->lock);
    return 0;
  }
  shared->solver = solver;
  shared->thread_count = pool->thread_count;
  shared->value_prep_input = values;
  shared->value_prep_output = solver->values;
  shared->value_prep_input_to_csc = solver->prepared_value_input_pos;
  shared->value_prep_nnz = solver->nnz;
  shared->value_prep_mode = 1;
  for (int tid = 0; tid < pool->thread_count; ++tid) {
    pool->workers[tid].shared = shared;
  }
  kls_egraph_pool_dispatch_and_spin_wait(pool, shared, pool->thread_count);
  shared->value_prep_mode = 0;
  shared->value_prep_input = NULL;
  shared->value_prep_output = NULL;
  shared->value_prep_input_to_csc = NULL;
  shared->value_prep_nnz = 0u;
  pthread_mutex_unlock(&shared->lock);
  return 1;
}

static int kls_run_parallel_row_solve_factor(kls_solver *solver,
                                             kls_egraph_refactor_pool *pool,
                                             int thread_count,
                                             const kls_row_solve_factor_view *view,
                                             double *work) {
  if (solver == NULL || pool == NULL || thread_count < 2 ||
      view == NULL || work == NULL) {
    return 0;
  }
  const UF_long *bounds = view->slice_bounds;
  const UF_long *ptr = view->ptr;
  const UF_long *cols = view->cols;
  const UF_long *split = view->segment_split;
  const UF_long *thread_bounds = view->thread_bounds;
  if (bounds == NULL || ptr == NULL || cols == NULL ||
      !kls_row_solve_view_has_values(solver, view) ||
      split == NULL || thread_bounds == NULL ||
      view->thread_count != thread_count ||
      (view->diagonal && solver->numeric->Udiag == NULL)) {
    return 0;
  }
  for (UF_long slice = 0u; slice < KLS_ROW_SOLVE_TRAPEZOID_SLICES; ++slice) {
    const UF_long begin = bounds[slice];
    const UF_long end = bounds[slice + 1u];
    if (begin > end || end > solver->n) {
      return 0;
    }
    for (UF_long row = begin; row < end; ++row) {
      if (ptr[row] > split[row] || split[row] > ptr[row + 1u]) {
        return 0;
      }
    }
    const UF_long *parts =
      thread_bounds + slice * ((UF_long)thread_count + 1u);
    if (parts[0] != begin || parts[thread_count] != end) {
      return 0;
    }
    for (int worker = 0; worker < thread_count; ++worker) {
      if (parts[worker] > parts[worker + 1] ||
          parts[worker + 1] > end) {
        return 0;
      }
    }
  }

  kls_egraph_refactor_shared *shared = &pool->shared;
  pthread_mutex_lock(&shared->lock);
  if (pool->active_workers != 0) {
    pthread_mutex_unlock(&shared->lock);
    return 0;
  }

  shared->solver = solver;
  shared->thread_count = thread_count;
  shared->row_refactor_mode = 0;
  shared->row_solve_mode = 1;
  shared->row_solve_upper = view->upper ? 1 : 0;
  shared->row_solve_view = view;
  shared->row_solve_work = work;
  shared->row_refactor_defer_value_scatter = 0;
  shared->row_refactor_lazy_value_scatter = 0;
  shared->row_refactor_native_row_panel_active = 0;
  shared->row_pipeline_ready_queue = 0;
  shared->row_pipeline_row_dep_ready_queue = 0;
  shared->row_pipeline_ready_groups = NULL;
  shared->row_pipeline_ready_slots = NULL;
  shared->row_pipeline_remaining_preds = NULL;
  shared->row_pipeline_tail_groups = NULL;
  shared->row_pipeline_tail_count = 0;
  shared->row_pipeline_private_groups = NULL;
  shared->row_pipeline_private_thread_ptr = NULL;
  shared->row_pipeline_private_group_mask = NULL;
  shared->row_pipeline_private_count = 0;
  shared->row_pipeline_ordered_private = 0;
  shared->pipeline_done = NULL;
  shared->pipeline_generation = 0;
  shared->pipeline_pos_end = 0;
  shared->cluster_level_count = 0;
  shared->invalid = 0;
  shared->pivot_rejected = 0;
  shared->singular = 0;
  atomic_store_explicit(&shared->stop, 0, memory_order_release);

  for (int i = 0; i < thread_count; ++i) {
    pool->workers[i].shared = shared;
  }

  kls_egraph_pool_dispatch_and_wait(pool, shared, thread_count);
  const int ok = !shared->invalid &&
                 !kls_egraph_refactor_should_stop(shared);
  shared->row_solve_mode = 0;
  shared->row_solve_view = NULL;
  shared->row_solve_work = NULL;
  pthread_mutex_unlock(&shared->lock);
  return ok;
}


static int kls_run_parallel_row_solve_sparse_block(
  kls_solver *solver,
  kls_egraph_refactor_pool *pool,
  int thread_count,
  const kls_row_solve_factor_view *view,
  double *work,
  UF_long *level_runs_out) {
  KLS_SET_OPTIONAL_OUTPUT(level_runs_out, 0);
  if (solver == NULL || pool == NULL || thread_count < 2 ||
      view == NULL || work == NULL) {
    return 0;
  }
  const UF_long *level_ptr = view->sparse_level_ptr;
  const UF_long *level_rows = view->sparse_level_rows;
  const UF_long level_count = view->sparse_level_count;
  UF_long cluster_levels = view->sparse_cluster_levels;
  const UF_long sparse_end = view->dense_tail_start;
  const UF_long *ptr = view->ptr;
  const UF_long *cols = view->cols;
  if (level_ptr == NULL || level_rows == NULL || level_count == 0u ||
      cluster_levels == 0u || sparse_end == 0u || sparse_end > solver->n ||
      ptr == NULL || cols == NULL ||
      !kls_row_solve_view_has_values(solver, view) ||
      view->thread_count != thread_count ||
      (view->diagonal && solver->numeric->Udiag == NULL)) {
    return 0;
  }
  if (cluster_levels > level_count) {
    cluster_levels = level_count;
  }
  if (level_ptr[0] != 0u || level_ptr[level_count] != sparse_end) {
    return 0;
  }
  for (UF_long level = 0u; level < level_count; ++level) {
    if (level_ptr[level] > level_ptr[level + 1u] ||
        level_ptr[level + 1u] > sparse_end) {
      return 0;
    }
  }
  for (UF_long pos = 0u; pos < sparse_end; ++pos) {
    if (level_rows[pos] >= sparse_end) {
      return 0;
    }
  }

  kls_egraph_refactor_shared *shared = &pool->shared;
  pthread_mutex_lock(&shared->lock);
  if (pool->active_workers != 0) {
    pthread_mutex_unlock(&shared->lock);
    return 0;
  }

  shared->solver = solver;
  shared->thread_count = thread_count;
  shared->row_refactor_mode = 0;
  shared->row_solve_mode = 2;
  shared->row_refactor_native_row_panel_active = 0;
  shared->row_solve_upper = view->upper ? 1 : 0;
  shared->row_solve_view = view;
  shared->row_solve_work = work;
  shared->row_refactor_defer_value_scatter = 0;
  shared->row_refactor_lazy_value_scatter = 0;
  shared->row_pipeline_ready_queue = 0;
  shared->row_pipeline_row_dep_ready_queue = 0;
  shared->row_pipeline_ready_groups = NULL;
  shared->row_pipeline_ready_slots = NULL;
  shared->row_pipeline_remaining_preds = NULL;
  shared->row_pipeline_tail_groups = NULL;
  shared->row_pipeline_tail_count = 0;
  shared->row_pipeline_private_groups = NULL;
  shared->row_pipeline_private_thread_ptr = NULL;
  shared->row_pipeline_private_group_mask = NULL;
  shared->row_pipeline_private_count = 0;
  shared->row_pipeline_ordered_private = 0;
  shared->pipeline_done = NULL;
  shared->pipeline_generation = 0;
  shared->pipeline_pos_end = 0;
  shared->cluster_level_count = 0;
  shared->invalid = 0;
  shared->pivot_rejected = 0;
  shared->singular = 0;
  atomic_store_explicit(&shared->stop, 0, memory_order_release);

  for (int i = 0; i < thread_count; ++i) {
    pool->workers[i].shared = shared;
  }

  kls_egraph_pool_dispatch_and_wait(pool, shared, thread_count);
  const int ok = !shared->invalid &&
                 !kls_egraph_refactor_should_stop(shared);
  shared->row_solve_mode = 0;
  shared->row_solve_view = NULL;
  shared->row_solve_work = NULL;
  pthread_mutex_unlock(&shared->lock);
  if (ok && level_runs_out != NULL) {
    *level_runs_out = cluster_levels;
  }
  return ok;
}

static void kls_row_solve_forward_rows(const kls_solver *solver,
                                       const kls_row_solve_factor_view *view,
                                       double *work,
                                       UF_long begin,
                                       UF_long end) {
  for (UF_long row = begin; row < end; ++row) {
    kls_row_solve_view_row(solver, view, work, row);
  }
}

static void kls_row_solve_backward_rows(const kls_solver *solver,
                                        const kls_row_solve_factor_view *view,
                                        double *work,
                                        UF_long begin,
                                        UF_long end) {
  for (UF_long remaining = end; remaining > begin; --remaining) {
    const UF_long row = remaining - 1u;
    kls_row_solve_view_row(solver, view, work, row);
  }
}

static int kls_row_solve_parallel_view(kls_solver *solver,
                                       kls_egraph_refactor_pool *pool,
                                       int thread_count,
                                       const kls_row_solve_factor_view *view,
                                       double *work,
                                       UF_long *slice_runs_out,
                                       UF_long *sparse_level_runs_out) {
  KLS_SET_OPTIONAL_OUTPUT(slice_runs_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(sparse_level_runs_out, 0);
  if (solver == NULL || view == NULL || work == NULL ||
      view->ptr == NULL ||
      !kls_row_solve_view_has_values(solver, view) ||
      (view->diagonal && solver->numeric->Udiag == NULL)) {
    return 0;
  }
  const UF_long start = view->dense_tail_start;
  if (start > solver->n) {
    return 0;
  }
  UF_long sparse_level_runs = 0;
  if (!view->upper) {
    if (start > 0u && pool != NULL &&
        kls_run_parallel_row_solve_sparse_block(
          solver, pool, thread_count, view, work, &sparse_level_runs)) {
      if (sparse_level_runs_out != NULL) {
        *sparse_level_runs_out = sparse_level_runs;
      }
    } else {
      kls_row_solve_forward_rows(solver, view, work, 0u, start);
    }

    if (pool == NULL ||
        view->slice_bounds == NULL ||
        view->segment_split == NULL ||
        view->rect_entries < KLS_ROW_SOLVE_PARALLEL_RECT_MIN_NNZ ||
        start >= solver->n) {
      kls_row_solve_forward_rows(solver, view, work, start, solver->n);
      return 1;
    }

    if (view->slice_bounds[0] != start ||
        !kls_run_parallel_row_solve_factor(
          solver, pool, thread_count, view, work)) {
      return 0;
    }
    if (slice_runs_out != NULL) {
      *slice_runs_out = KLS_ROW_SOLVE_TRAPEZOID_SLICES;
    }
    return 1;
  }

  if (pool == NULL ||
      view->slice_bounds == NULL ||
      view->segment_split == NULL ||
      view->rect_entries < KLS_ROW_SOLVE_PARALLEL_RECT_MIN_NNZ ||
      start >= solver->n) {
    kls_row_solve_backward_rows(solver, view, work, start, solver->n);
    if (start > 0u && pool != NULL &&
        kls_run_parallel_row_solve_sparse_block(
          solver, pool, thread_count, view, work, &sparse_level_runs)) {
      if (sparse_level_runs_out != NULL) {
        *sparse_level_runs_out = sparse_level_runs;
      }
    } else {
      kls_row_solve_backward_rows(solver, view, work, 0u, start);
    }
    return 1;
  }

  if (view->slice_bounds[0] != start ||
      !kls_run_parallel_row_solve_factor(
        solver, pool, thread_count, view, work)) {
    return 0;
  }
  if (start > 0u &&
      kls_run_parallel_row_solve_sparse_block(
        solver, pool, thread_count, view, work, &sparse_level_runs)) {
    if (sparse_level_runs_out != NULL) {
      *sparse_level_runs_out = sparse_level_runs;
    }
  } else {
    kls_row_solve_backward_rows(solver, view, work, 0u, start);
  }
  if (slice_runs_out != NULL) {
    *slice_runs_out = KLS_ROW_SOLVE_TRAPEZOID_SLICES;
  }
  return 1;
}

static void kls_make_normal_row_solve_view(
  const kls_solver *solver,
  int upper,
  kls_row_solve_factor_view *view) {
  if (view == NULL) {
    return;
  }
  memset(view, 0, sizeof(*view));
  if (solver == NULL) {
    return;
  }
  view->ptr = upper ? solver->row_refactor_u_ptr : solver->row_refactor_l_ptr;
  view->cols =
    upper ? solver->row_refactor_u_cols : solver->row_refactor_l_cols;
  if (!kls_row_refactor_solve_uses_direct_values(solver) &&
      solver->row_refactor_values_ready) {
    view->values = upper ? solver->row_refactor_u_row_values
                         : solver->row_refactor_l_row_values;
  }
  view->slice_bounds =
    upper ? solver->row_solve_u_slice_bounds
          : solver->row_solve_l_slice_bounds;
  view->segment_split =
    upper ? solver->row_solve_u_segment_split
          : solver->row_solve_l_segment_split;
  view->thread_bounds =
    upper ? solver->row_solve_u_thread_bounds
          : solver->row_solve_l_thread_bounds;
  view->sparse_level_ptr =
    upper ? solver->row_solve_u_sparse_level_ptr
          : solver->row_solve_l_sparse_level_ptr;
  view->sparse_level_rows =
    upper ? solver->row_solve_u_sparse_level_rows
          : solver->row_solve_l_sparse_level_rows;
  view->sparse_level_count =
    upper ? solver->row_solve_u_sparse_level_count
          : solver->row_solve_l_sparse_level_count;
  view->sparse_cluster_levels =
    upper ? solver->row_solve_u_sparse_cluster_levels
          : solver->row_solve_l_sparse_cluster_levels;
  view->dense_tail_start =
    upper ? solver->row_solve_u_dense_tail_start
          : solver->row_solve_l_dense_tail_start;
  view->rect_entries =
    upper ? solver->row_solve_u_rect_entries
          : solver->row_solve_l_rect_entries;
  view->thread_count = solver->row_solve_thread_count;
  view->source_upper = upper ? 1 : 0;
  view->upper = upper ? 1 : 0;
  view->diagonal = upper ? 1 : 0;
}

static void kls_make_transpose_row_solve_view(
  const kls_row_solve_transpose_plan *plan,
  kls_row_solve_factor_view *view) {
  if (view == NULL) {
    return;
  }
  memset(view, 0, sizeof(*view));
  if (plan == NULL) {
    return;
  }
  view->ptr = plan->ptr;
  view->cols = plan->cols;
  view->source_pos = plan->source_pos;
  view->slice_bounds = plan->slice_bounds;
  view->segment_split = plan->segment_split;
  view->thread_bounds = plan->thread_bounds;
  view->sparse_level_ptr = plan->sparse_level_ptr;
  view->sparse_level_rows = plan->sparse_level_rows;
  view->sparse_level_count = plan->sparse_level_count;
  view->sparse_cluster_levels = plan->sparse_cluster_levels;
  view->dense_tail_start = plan->dense_tail_start;
  view->rect_entries = plan->rect_entries;
  view->thread_count = plan->thread_count;
  view->source_upper = plan->source_upper;
  view->upper = plan->upper;
  view->diagonal = plan->diagonal;
}


static int kls_try_parallel_row_solve_one_rhs(kls_solver *solver, double *x) {
  if (solver == NULL || x == NULL || solver->symbolic == NULL ||
      solver->numeric == NULL || solver->options.threads <= 1 ||
      solver->symbolic->nblocks != 1u ||
      solver->row_solve_partition_ready == 0 ||
      (solver->row_solve_l_segment_split == NULL &&
       solver->row_solve_u_segment_split == NULL &&
       solver->row_solve_l_sparse_level_ptr == NULL &&
       solver->row_solve_u_sparse_level_ptr == NULL) ||
      solver->numeric->Xwork == NULL || solver->numeric->Pnum == NULL ||
      solver->symbolic->Q == NULL || solver->numeric->Udiag == NULL) {
    return 0;
  }

  int thread_count = solver->options.threads;
  const UF_long rectangular_entries =
    (solver->row_solve_l_segment_split != NULL
       ? solver->row_solve_l_rect_entries : 0u) +
    (solver->row_solve_u_segment_split != NULL
       ? solver->row_solve_u_rect_entries : 0u);
  const UF_long sparse_entries =
    (solver->row_solve_l_sparse_level_ptr != NULL
       ? solver->row_refactor_l_ptr[solver->row_solve_l_dense_tail_start]
       : 0u) +
    (solver->row_solve_u_sparse_level_ptr != NULL
       ? solver->row_refactor_u_ptr[solver->row_solve_u_dense_tail_start]
       : 0u);
  if (rectangular_entries + sparse_entries <
      KLS_ROW_SOLVE_PARALLEL_RECT_MIN_NNZ) {
    return 0;
  }

  UF_long max_rows = 0;
  if (solver->row_solve_l_segment_split != NULL &&
      solver->row_solve_l_dense_tail_rows > max_rows) {
    max_rows = solver->row_solve_l_dense_tail_rows;
  }
  if (solver->row_solve_u_segment_split != NULL &&
      solver->row_solve_u_dense_tail_rows > max_rows) {
    max_rows = solver->row_solve_u_dense_tail_rows;
  }
  if (solver->row_solve_l_sparse_level_ptr != NULL &&
      solver->row_solve_l_dense_tail_start > max_rows) {
    max_rows = solver->row_solve_l_dense_tail_start;
  }
  if (solver->row_solve_u_sparse_level_ptr != NULL &&
      solver->row_solve_u_dense_tail_start > max_rows) {
    max_rows = solver->row_solve_u_dense_tail_start;
  }
  if ((UF_long)thread_count > max_rows) {
    thread_count = (int)max_rows;
  }
  if (thread_count < 2) {
    return 0;
  }
  const UF_long l_sparse_cluster_levels =
    kls_row_solve_sparse_cluster_levels(
      solver->row_solve_l_sparse_level_ptr,
      solver->row_solve_l_sparse_level_count, thread_count);
  const UF_long u_sparse_cluster_levels =
    kls_row_solve_sparse_cluster_levels(
      solver->row_solve_u_sparse_level_ptr,
      solver->row_solve_u_sparse_level_count, thread_count);
  const UF_long active_l_rect_entries =
    solver->row_solve_l_segment_split != NULL &&
        solver->row_solve_l_rect_entries >=
          KLS_ROW_SOLVE_PARALLEL_RECT_MIN_NNZ
      ? solver->row_solve_l_rect_entries
      : 0u;
  const UF_long active_u_rect_entries =
    solver->row_solve_u_segment_split != NULL &&
        solver->row_solve_u_rect_entries >=
          KLS_ROW_SOLVE_PARALLEL_RECT_MIN_NNZ
      ? solver->row_solve_u_rect_entries
      : 0u;
  const UF_long active_rect_entries =
    active_l_rect_entries + active_u_rect_entries;
  const UF_long active_sparse_entries =
    (l_sparse_cluster_levels > 0u
       ? solver->row_refactor_l_ptr[solver->row_solve_l_dense_tail_start]
       : 0u) +
    (u_sparse_cluster_levels > 0u
       ? solver->row_refactor_u_ptr[solver->row_solve_u_dense_tail_start]
       : 0u);
  const UF_long active_parallel_entries =
    active_rect_entries + active_sparse_entries;
  const UF_long total_solve_entries =
    solver->row_refactor_l_ptr[solver->n] +
    solver->row_refactor_u_ptr[solver->n];
  const UF_long sync_count =
    (active_l_rect_entries > 0u
       ? 2u * KLS_ROW_SOLVE_TRAPEZOID_SLICES
       : 0u) +
    (active_u_rect_entries > 0u
       ? 2u * KLS_ROW_SOLVE_TRAPEZOID_SLICES
       : 0u) +
    l_sparse_cluster_levels +
    u_sparse_cluster_levels;
  if (active_parallel_entries == 0u ||
      active_parallel_entries * 5u < total_solve_entries * 2u ||
      sync_count == 0u ||
      active_parallel_entries / sync_count <
        KLS_ROW_SOLVE_PARALLEL_MIN_ENTRIES_PER_SYNC) {
    return 0;
  }
  if (!kls_build_row_solve_thread_bounds(solver, thread_count)) {
    return 0;
  }
  kls_row_solve_factor_view l_view;
  kls_row_solve_factor_view u_view;
  kls_make_normal_row_solve_view(solver, 0, &l_view);
  kls_make_normal_row_solve_view(solver, 1, &u_view);

  kls_egraph_refactor_pool *pool =
    ensure_egraph_refactor_pool(solver, thread_count);
  if (pool == NULL) {
    return 0;
  }

  const UF_long n = solver->n;
  const UF_long *pnum = solver->numeric->Pnum;
  const UF_long *q = solver->symbolic->Q;
  const double *rs = solver->numeric->Rs;
  double *work = (double *)solver->numeric->Xwork;
  for (UF_long k = 0; k < n; ++k) {
    const double rhs_value = x[pnum[k]];
    work[k] = rs != NULL ? rhs_value / rs[k] : rhs_value;
  }

  UF_long l_slice_runs = 0;
  UF_long u_slice_runs = 0;
  UF_long l_sparse_level_runs = 0;
  UF_long u_sparse_level_runs = 0;
  if (!kls_row_solve_parallel_view(solver, pool, thread_count, &l_view, work,
                                   &l_slice_runs, &l_sparse_level_runs) ||
      !kls_row_solve_parallel_view(solver, pool, thread_count, &u_view, work,
                                   &u_slice_runs, &u_sparse_level_runs)) {
    return 0;
  }

  for (UF_long k = 0; k < n; ++k) {
    x[q[k]] = work[k];
  }
  if (l_slice_runs > 0u || u_slice_runs > 0u ||
      l_sparse_level_runs > 0u || u_sparse_level_runs > 0u) {
    kls_record_row_solve_parallel_run(solver, l_slice_runs, u_slice_runs,
                                      l_sparse_level_runs,
                                      u_sparse_level_runs);
  }
  return l_slice_runs > 0u || u_slice_runs > 0u ||
         l_sparse_level_runs > 0u || u_sparse_level_runs > 0u;
}

static KLS_ALWAYS_INLINE int kls_packed_factor_column_scatter(
  double *restrict work,
  const UF_long *restrict indices,
  const double *restrict values,
  UF_long length,
  UF_long bound,
  double alpha) {
  UF_long p = 0u;
  for (; p + 3u < length; p += 4u) {
    if (indices[p] >= bound || indices[p + 1u] >= bound ||
        indices[p + 2u] >= bound || indices[p + 3u] >= bound) {
      return 0;
    }
    work[indices[p]] = fma(-values[p], alpha, work[indices[p]]);
    work[indices[p + 1u]] =
      fma(-values[p + 1u], alpha, work[indices[p + 1u]]);
    work[indices[p + 2u]] =
      fma(-values[p + 2u], alpha, work[indices[p + 2u]]);
    work[indices[p + 3u]] =
      fma(-values[p + 3u], alpha, work[indices[p + 3u]]);
  }
  for (; p < length; ++p) {
    if (indices[p] >= bound) {
      return 0;
    }
    work[indices[p]] = fma(-values[p], alpha, work[indices[p]]);
  }
  return 1;
}

/* The fused refactor already formed the unit-lower result.  Finish the one
   admitted block directly from KLU's authoritative packed U columns. */
static int kls_fused_packed_upper_solve(kls_solver *solver,
                                         double *work,
                                         double *out,
                                         const UF_long *qperm) {
  if (solver == NULL || solver->symbolic == NULL ||
      solver->numeric == NULL || work == NULL ||
      solver->symbolic->nblocks != 1u ||
      solver->numeric->Uip == NULL || solver->numeric->Ulen == NULL ||
      solver->numeric->LUbx == NULL ||
      solver->numeric->LUbx[0] == NULL ||
      solver->numeric->Udiag == NULL || out == NULL || qperm == NULL) {
    return 0;
  }
  const UF_long n = solver->n;
  const UF_long *restrict uip = solver->numeric->Uip;
  const UF_long *restrict ulen = solver->numeric->Ulen;
  const double *restrict udiag =
    (const double *)solver->numeric->Udiag;
  double *restrict lu = (double *)solver->numeric->LUbx[0];

  for (UF_long remaining = n; remaining > 0u; --remaining) {
    const UF_long row = remaining - 1u;
    if (udiag[row] == 0.0 || !isfinite(udiag[row])) {
      return 0;
    }
    const double xk = work[row] / udiag[row];
    work[row] = xk;
    if (qperm[row] >= n) return 0;
    out[qperm[row]] = xk;
    if (xk == 0.0) {
      continue;
    }
    UF_long *indices = NULL;
    double *column_values = NULL;
    UF_long length = 0u;
    kls_klu_get_pointer(lu, uip, ulen, row,
                        &indices, &column_values, &length);
    if (!kls_packed_factor_column_scatter(
          work, indices, column_values, length, row, xk)) {
      return 0;
    }
  }
  return 1;
}

static KLS_ALWAYS_INLINE double kls_fused_i16_row_dot(
  const double *restrict values,
  const uint16_t *restrict cols,
  const double *restrict work,
  UF_long begin,
  UF_long end,
  double value) {
  const UF_long len = end - begin;
  switch (len) {
    case 0u:
      return value;
    case 1u:
      return fma(-values[begin], work[cols[begin]], value);
    case 2u:
      value = fma(-values[begin], work[cols[begin]], value);
      return fma(-values[begin + 1u], work[cols[begin + 1u]], value);
    case 3u:
      value = fma(-values[begin], work[cols[begin]], value);
      value = fma(-values[begin + 1u], work[cols[begin + 1u]], value);
      return fma(-values[begin + 2u], work[cols[begin + 2u]], value);
    default:
      break;
  }
  UF_long p = begin;
  if (len >= 4u) {
    double sum0 = value;
    double sum1 = 0.0;
    double sum2 = 0.0;
    double sum3 = 0.0;
    for (; p + 3u < end; p += 4u) {
      sum0 = fma(-values[p], work[cols[p]], sum0);
      sum1 = fma(-values[p + 1u], work[cols[p + 1u]], sum1);
      sum2 = fma(-values[p + 2u], work[cols[p + 2u]], sum2);
      sum3 = fma(-values[p + 3u], work[cols[p + 3u]], sum3);
    }
    value = (sum0 + sum1) + (sum2 + sum3);
  }
  for (; p < end; ++p) {
    value = fma(-values[p], work[cols[p]], value);
  }
  return value;
}

static int kls_fused_row_values_upper_solve(kls_solver *solver,
                                             double *work,
                                             double *out,
                                             const UF_long *qperm) {
  if (solver == NULL || solver->symbolic == NULL ||
      solver->numeric == NULL || work == NULL ||
      solver->symbolic->nblocks != 1u ||
      solver->row_refactor_u_ptr == NULL ||
      solver->row_refactor_u_cols16 == NULL ||
      solver->row_refactor_u_row_values == NULL ||
      solver->numeric->Udiag == NULL || out == NULL || qperm == NULL) {
    return 0;
  }
  const UF_long *restrict uptr = solver->row_refactor_u_ptr;
  const uint16_t *restrict ucols = solver->row_refactor_u_cols16;
  const double *restrict uvalues = solver->row_refactor_u_row_values;
  const double *restrict udiag =
    (const double *)solver->numeric->Udiag;
  for (UF_long remaining = solver->n; remaining > 0u; --remaining) {
    const UF_long row = remaining - 1u;
    if (udiag[row] == 0.0 || !isfinite(udiag[row])) {
      return 0;
    }
    const double value = kls_fused_i16_row_dot(
      uvalues, ucols, work, uptr[row], uptr[row + 1u], work[row]);
    work[row] = value / udiag[row];
    if (qperm[row] >= solver->n) return 0;
    out[qperm[row]] = work[row];
  }
  return 1;
}

/* Finish a normal solve whose terminal-block unit-lower sweep was formed while
   the changed L rows were published.  This helper deliberately enters the
   ordinary solve dispatcher (rather than being a separate unchecked solve):
   refinement and residual contracts after the raw kernel remain unchanged. */
static int kls_try_fused_refactor_upper_solve(kls_solver *solver,
                                               double *out) {
  if (solver == NULL || out == NULL ||
      !solver->fused_refactor_solve_ready ||
      solver->fused_refactor_solve_work == NULL ||
      solver->fused_refactor_solve_work_n != solver->n ||
      solver->symbolic == NULL || solver->numeric == NULL ||
      solver->symbolic->nblocks != 1u || solver->symbolic->R == NULL ||
      solver->symbolic->Q == NULL ||
      solver->numeric->Udiag == NULL ||
      solver->row_refactor_u_ptr == NULL ||
      solver->row_refactor_u_cols == NULL ||
      solver->row_refactor_u_values == NULL ||
      (solver->fused_refactor_solve_row_values &&
       solver->row_refactor_u_row_values == NULL)) {
    return 0;
  }

  /* A prepared generation is single-use even if a defensive metadata check
     below fails; a later unrelated solve must never consume an older RHS. */
  solver->fused_refactor_solve_ready = 0;
  const UF_long n = solver->n;
  const UF_long *restrict r = solver->symbolic->R;
  if (r[0] != 0u || r[1] != n) {
    return 0;
  }
  const int use_row_values = solver->fused_refactor_solve_row_values;
  const UF_long *restrict q = solver->symbolic->Q;
  double *restrict work = solver->fused_refactor_solve_work;

  const int upper_ok = use_row_values
    ? kls_fused_row_values_upper_solve(solver, work, out, q)
    : kls_fused_packed_upper_solve(solver, work, out, q);
  if (!upper_ok) {
    return 0;
  }
  return 1;
}

static int kls_try_parallel_row_solve_transpose_one_rhs(kls_solver *solver,
                                                        double *x) {
  if (solver == NULL || x == NULL || solver->symbolic == NULL ||
      solver->numeric == NULL || solver->options.threads <= 1 ||
      solver->symbolic->nblocks != 1u ||
      solver->row_solve_partition_ready == 0 ||
      solver->numeric->Xwork == NULL || solver->numeric->Pnum == NULL ||
      solver->symbolic->Q == NULL || solver->numeric->Udiag == NULL) {
    return 0;
  }

  if (solver->row_solve_ut_plan.serial_only ||
      solver->row_solve_lt_plan.serial_only) {
    /* A dirty lean solve may have built only the immutable transpose graph.
       Materialize the scheduling metadata if a later factor asks for the
       parallel row solver. */
    kls_build_row_solve_transpose_plans(solver);
  }

  const kls_row_solve_transpose_plan *ut = &solver->row_solve_ut_plan;
  const kls_row_solve_transpose_plan *lt = &solver->row_solve_lt_plan;
  const int have_ut = ut->ptr != NULL &&
                      (ut->segment_split != NULL ||
                       ut->sparse_level_ptr != NULL);
  const int have_lt = lt->ptr != NULL &&
                      (lt->segment_split != NULL ||
                       lt->sparse_level_ptr != NULL);
  if (!have_ut && !have_lt) {
    return 0;
  }

  int thread_count = solver->options.threads;
  const UF_long rectangular_entries =
    (have_ut && ut->segment_split != NULL ? ut->rect_entries : 0u) +
    (have_lt && lt->segment_split != NULL ? lt->rect_entries : 0u);
  const UF_long sparse_entries =
    (have_ut && ut->sparse_level_ptr != NULL && ut->ptr != NULL
       ? ut->ptr[ut->dense_tail_start] : 0u) +
    (have_lt && lt->sparse_level_ptr != NULL && lt->ptr != NULL
       ? lt->ptr[lt->dense_tail_start] : 0u);
  if (rectangular_entries + sparse_entries <
      KLS_ROW_SOLVE_PARALLEL_RECT_MIN_NNZ) {
    return 0;
  }

  UF_long max_rows = 0;
  if (have_ut && ut->segment_split != NULL &&
      ut->dense_tail_rows > max_rows) {
    max_rows = ut->dense_tail_rows;
  }
  if (have_lt && lt->segment_split != NULL &&
      lt->dense_tail_rows > max_rows) {
    max_rows = lt->dense_tail_rows;
  }
  if (have_ut && ut->sparse_level_ptr != NULL &&
      ut->dense_tail_start > max_rows) {
    max_rows = ut->dense_tail_start;
  }
  if (have_lt && lt->sparse_level_ptr != NULL &&
      lt->dense_tail_start > max_rows) {
    max_rows = lt->dense_tail_start;
  }
  if ((UF_long)thread_count > max_rows) {
    thread_count = (int)max_rows;
  }
  if (thread_count < 2) {
    return 0;
  }

  const UF_long ut_sparse_cluster_levels =
    kls_row_solve_sparse_cluster_levels(ut->sparse_level_ptr,
                                        ut->sparse_level_count,
                                        thread_count);
  const UF_long lt_sparse_cluster_levels =
    kls_row_solve_sparse_cluster_levels(lt->sparse_level_ptr,
                                        lt->sparse_level_count,
                                        thread_count);
  const UF_long active_ut_rect_entries =
    have_ut && ut->segment_split != NULL &&
        ut->rect_entries >= KLS_ROW_SOLVE_PARALLEL_RECT_MIN_NNZ
      ? ut->rect_entries : 0u;
  const UF_long active_lt_rect_entries =
    have_lt && lt->segment_split != NULL &&
        lt->rect_entries >= KLS_ROW_SOLVE_PARALLEL_RECT_MIN_NNZ
      ? lt->rect_entries : 0u;
  const UF_long active_parallel_entries =
    active_ut_rect_entries + active_lt_rect_entries +
    (ut_sparse_cluster_levels > 0u && ut->ptr != NULL
       ? ut->ptr[ut->dense_tail_start] : 0u) +
    (lt_sparse_cluster_levels > 0u && lt->ptr != NULL
       ? lt->ptr[lt->dense_tail_start] : 0u);
  const UF_long total_solve_entries =
    solver->row_refactor_l_ptr[solver->n] +
    solver->row_refactor_u_ptr[solver->n];
  const UF_long sync_count =
    (active_ut_rect_entries > 0u
       ? 2u * KLS_ROW_SOLVE_TRAPEZOID_SLICES
       : 0u) +
    (active_lt_rect_entries > 0u
       ? 2u * KLS_ROW_SOLVE_TRAPEZOID_SLICES
       : 0u) +
    ut_sparse_cluster_levels +
    lt_sparse_cluster_levels;
  if (active_parallel_entries == 0u ||
      active_parallel_entries * 5u < total_solve_entries * 2u ||
      sync_count == 0u ||
      active_parallel_entries / sync_count <
        KLS_ROW_SOLVE_PARALLEL_MIN_ENTRIES_PER_SYNC) {
    return 0;
  }
  if (!kls_build_row_solve_transpose_thread_bounds(solver, thread_count)) {
    return 0;
  }

  kls_row_solve_factor_view ut_view;
  kls_row_solve_factor_view lt_view;
  kls_make_transpose_row_solve_view(&solver->row_solve_ut_plan, &ut_view);
  kls_make_transpose_row_solve_view(&solver->row_solve_lt_plan, &lt_view);

  kls_egraph_refactor_pool *pool =
    ensure_egraph_refactor_pool(solver, thread_count);
  if (pool == NULL) {
    return 0;
  }

  const UF_long n = solver->n;
  const UF_long *pnum = solver->numeric->Pnum;
  const UF_long *q = solver->symbolic->Q;
  const double *rs = solver->numeric->Rs;
  double *work = (double *)solver->numeric->Xwork;
  for (UF_long k = 0; k < n; ++k) {
    work[k] = x[q[k]];
  }

  UF_long ut_slice_runs = 0;
  UF_long lt_slice_runs = 0;
  UF_long ut_sparse_level_runs = 0;
  UF_long lt_sparse_level_runs = 0;
  if (!kls_row_solve_parallel_view(solver, pool, thread_count, &ut_view,
                                   work, &ut_slice_runs,
                                   &ut_sparse_level_runs) ||
      !kls_row_solve_parallel_view(solver, pool, thread_count, &lt_view,
                                   work, &lt_slice_runs,
                                   &lt_sparse_level_runs)) {
    return 0;
  }

  for (UF_long k = 0; k < n; ++k) {
    const double result = work[k];
    x[pnum[k]] = rs != NULL ? result / rs[k] : result;
  }
  if (ut_slice_runs > 0u || lt_slice_runs > 0u ||
      ut_sparse_level_runs > 0u || lt_sparse_level_runs > 0u) {
    kls_record_row_solve_parallel_run(solver, lt_slice_runs, ut_slice_runs,
                                      lt_sparse_level_runs,
                                      ut_sparse_level_runs);
  }
  return ut_slice_runs > 0u || lt_slice_runs > 0u ||
         ut_sparse_level_runs > 0u || lt_sparse_level_runs > 0u;
}

static const UF_long *kls_refactor_snode_run_end(
  const kls_solver *solver) {
  if (solver == NULL || solver->snode_run_end == NULL) {
    return NULL;
  }

  /* Batch only when the realized run census exposes at least one useful row
     per worker.  This replaces the former n/flop/single-block window with a
     direct capability and amortization check. */
  return solver->refactor_supernode_candidate_rows >=
           (UF_long)(2 * solver->options.threads)
    ? solver->snode_run_end : NULL;
}

static int kls_egraph_refactor_storage_is_eligible(
  const kls_solver *solver,
  int single_block) {
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

/* A complete separator-private prefix can still leave a small number of
   very expensive separator columns in the dependency pipeline.  In that
   case prefix width is already settled by the private plan, but adjacent
   pipeline-column fusion remains an independent equivalent executor choice.
   Consult it only for a recurring numeric whose realized tail owns at least
   half of the modeled update work, contains enough columns to amortize the
   two measured probe updates, and averages enough dependency work per column
   for shared producer streaming to dominate deck bookkeeping.  This is a
   retained-task-graph capability, independent of matrix identity or
   dimensions. */
static int kls_egraph_pipeline_fusion_trial_worthwhile(
  const kls_solver *solver) {
  return solver != NULL &&
    kls_repeated_update_workload(&solver->options) &&
    solver->options.threads >= 2 &&
    solver->refactor_dependency_work > 0.0 &&
    solver->refactor_pipeline_work >=
      0.5 * solver->refactor_dependency_work &&
    solver->refactor_pipeline_column_count >=
      (UF_long)(4 * solver->options.threads) &&
    solver->refactor_pipeline_work >=
      1.0e6 * (double)solver->refactor_pipeline_column_count;
}

static int kls_egraph_refactor_is_eligible(const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      getenv("KLS_DISABLE_EGRAPH_REFACTOR") != NULL ||
      solver->options.threads <= 1 ||
      solver->refactor_level_ptr == NULL ||
      solver->refactor_level_count == 0u ||
      solver->refactor_level_max_width < (UF_long)(4 * solver->options.threads)) {
    return 0;
  }
  if (solver->dense_tail_cols > 0) {
    /* pipe-emitted dense-tail numerics (see the schedule gate): the
       tail no-ops in the dispatch, the BLAS3 refresh rebuilds it */
    return solver->common.flops >= 5.0e8;
  }

  const int single_block = solver->symbolic->nblocks == 1u;
  {
    /* The schedule has already measured the fixed-pivot dependency work and
       its widest runnable level.  Require enough total work to amortize one
       dispatch per worker and enough independent work relative to the most
       expensive column.  These are machine-resource tests, not input-shape
       identities; later timing trials remain free to reject the engine. */
    const double min_dependency_work =
      KLS_EGRAPH_REFACTOR_MIN_FLOPS_PER_THREAD *
      (double)solver->options.threads;
    /* One heavy dependency column limits ideal crew-wide speedup, but it
       does not serialize the remaining schedule.  Admit factors with enough
       independent work to keep at least half the requested crew useful;
       the measured fusion trial still selects dispatch granularity.
       This is a retained task-graph resource test, independent of input
       dimensions or matrix identity. */
    const double independence_workers =
      0.5 * (double)solver->options.threads;
    /* The size floor excludes tiny problems dominated by thread dispatch;
       the work checks decide whether schedule construction can amortize. */
    if (solver->n < KLS_EGRAPH_REFACTOR_MIN_SIZE ||
        solver->refactor_level_cols == NULL ||
        solver->refactor_dependency_work < min_dependency_work ||
        (solver->refactor_dependency_max_column_work > 0.0 &&
         solver->refactor_dependency_work <
           independence_workers *
             solver->refactor_dependency_max_column_work)) {
      return 0;
    }
    return kls_egraph_refactor_storage_is_eligible(solver, single_block);
  }
}

static kls_egraph_refactor_kernel kls_egraph_refactor_kernel_for(
  const kls_solver *solver,
  int scale) {
  if (solver == NULL || solver->symbolic == NULL) {
    return KLS_EGRAPH_REFACTOR_KERNEL_GENERIC;
  }
  if (scale > 0) {
    return solver->symbolic->nblocks == 1u
      ? KLS_EGRAPH_REFACTOR_KERNEL_SINGLE_SCALED
      : KLS_EGRAPH_REFACTOR_KERNEL_GENERIC;
  }
  if (solver->symbolic->nblocks == 1u) {
    return KLS_EGRAPH_REFACTOR_KERNEL_SINGLE_UNSCALED;
  }
  if (solver->symbolic->nblocks > 1u) {
    return KLS_EGRAPH_REFACTOR_KERNEL_BTF_UNSCALED;
  }
  return KLS_EGRAPH_REFACTOR_KERNEL_GENERIC;
}

/* Keep the caller's full team. Once the separator layout settles, record
   a warmed baseline for the independent pair/quad pipeline-fusion trial. */
static int kls_egraph_steady_thread_count(kls_solver *solver,
                                          int thread_count,
                                          int check_pivots) {
  solver->eg_tt_pending = 0;
  if (check_pivots) {
    return thread_count;
  }
  solver->eg_tt_choice = thread_count;
  const int separator_settled = solver->eg_separator_choice != 0 &&
    !solver->eg_separator_pending;
  const int pipeline_fusion_trial = separator_settled &&
    thread_count >= 4 &&
    kls_egraph_pipeline_fusion_trial_worthwhile(solver);
  if (pipeline_fusion_trial) {
    if (solver->eg_tt_samples < 2) {
      solver->eg_tt_pending = 1;
    } else if (solver->eg_pair_choice == 0 &&
               !solver->eg_pair_pending) {
      solver->eg_pair_pending = 1;
    }
  } else if (separator_settled) {
    solver->eg_pair_choice = -1;
  }
  return thread_count;
}

static void kls_egraph_release_unused_quad_scratch(kls_solver *solver) {
  if (solver == NULL || solver->egraph_pool == NULL) {
    return;
  }
  kls_egraph_refactor_pool *pool = solver->egraph_pool;
  for (int tid = 0; tid < pool->thread_count; ++tid) {
    free(pool->workers[tid].fuse_x);
    pool->workers[tid].fuse_x = NULL;
  }
}

/* Compare two equivalent barrier/pipeline cuts on the retained dependency
   graph after the other EGraph executor choices have settled.  The alpha-3
   cut exposes a little more dependency-driven tail work and wins on many
   irregular factors, while the measured alpha-2 incumbent protects graphs
   whose wider barrier levels are cheaper. */
static UF_long kls_egraph_cluster_level_dispatch(kls_solver *solver,
                                                  int check_pivots) {
  if (solver == NULL) {
    return 0u;
  }
  const UF_long incumbent = solver->refactor_cluster_level_count;
  const UF_long alpha3 = solver->refactor_cluster_level_count_alpha3;
  const UF_long alpha4 = solver->refactor_cluster_level_count_alpha4;
  const int private_alpha4 =
    solver->eg_separator_choice > 0 &&
    solver->refactor_separator_private_cols != NULL &&
    solver->refactor_separator_private_thread_ptr != NULL &&
    solver->refactor_separator_private_cols_alpha4 != NULL &&
    solver->refactor_separator_private_thread_ptr_alpha4 != NULL;
  solver->eg_cluster_pending = 0;
  if (solver->eg_separator_choice == 0) {
    return incumbent;
  }
  if (check_pivots ||
      !kls_repeated_update_workload(&solver->options) ||
      (incumbent == alpha3 && incumbent == alpha4) ||
      alpha3 > incumbent || alpha4 > incumbent ||
      (solver->eg_separator_choice > 0 && !private_alpha4)) {
    solver->eg_cluster_choice = -1;
    return incumbent;
  }
  if (solver->eg_cluster_choice != 0) {
    return solver->eg_cluster_choice == 2 ? alpha4
      : solver->eg_cluster_choice == 1 ? alpha3 : incumbent;
  }
  /* Hold width, fusion, stream, and separator layout fixed so the samples
     differ only in where level barriers give way to dependency waits. */
  if (solver->eg_tt_choice <= 0 || solver->eg_pair_choice == 0 ||
      solver->eg_subset_choice == 0 || solver->eg_subset_pending ||
      solver->eg_stream_choice == 0 || solver->eg_premark_choice != 0 ||
      solver->eg_tt_pending || solver->eg_pair_pending ||
      solver->eg_stream_pending || solver->eg_separator_pending ||
      solver->eg_premark_pending) {
    return incumbent;
  }
  int side = 0;
  const int candidate_begin = private_alpha4 ? 2 : 1;
  const int candidate_end = private_alpha4 ? 3 : 2;
  for (int candidate = candidate_begin; candidate < candidate_end;
       ++candidate) {
    if (solver->eg_cluster_samples[candidate] <
        solver->eg_cluster_samples[side]) {
      side = candidate;
    }
  }
  solver->eg_cluster_pending = side + 1;
  return side == 2 ? alpha4 : side == 1 ? alpha3 : incumbent;
}

static int kls_egraph_cluster_premark_dispatch(
  kls_solver *solver,
  const kls_egraph_refactor_shared *shared,
  int check_pivots) {
  if (solver == NULL || shared == NULL) {
    return 0;
  }
  solver->eg_premark_pending = 0;
  if (check_pivots) {
    return 0;
  }
  if (!kls_repeated_update_workload(&solver->options)) {
    solver->eg_premark_choice = -1;
    return 0;
  }
  if (solver->eg_premark_choice != 0) {
    return solver->eg_premark_choice > 0;
  }

  /* Consult this equivalent publication layout only after the other EGraph
     executor trials have settled.  Otherwise their width, fusion, scatter,
     or separator changes would be charged to the completion bitmap. */
  if (solver->eg_separator_choice == 0 || solver->eg_tt_choice <= 0 ||
      solver->eg_pair_choice == 0 || solver->eg_stream_choice == 0 ||
      solver->eg_subset_choice == 0 || solver->eg_subset_pending ||
      solver->eg_cluster_choice == 0 ||
      solver->eg_separator_pending || solver->eg_tt_pending ||
      solver->eg_pair_pending || solver->eg_stream_pending ||
      solver->eg_cluster_pending) {
    return 0;
  }

  if (shared->pipeline_done == NULL || shared->pipeline_generation == 0u ||
      shared->cluster_level_count == 0u ||
      shared->cluster_level_count > solver->refactor_level_count ||
      solver->refactor_level_ptr == NULL ||
      solver->refactor_level_cols == NULL) {
    solver->eg_premark_choice = -1;
    return 0;
  }

  /* The immediately preceding stream-kernel consultation measured the exact
     retained width/fusion/kernel combination with ordinary publication.  Use
     that already-paid sample as the incumbent instead of charging the caller
     two duplicate ordinary refactors; keep two independent samples for the
     new arm before it can be adopted. */
  if (solver->eg_premark_samples[0] == 0 &&
      solver->eg_premark_samples[1] == 0) {
    const int stream_side = solver->eg_stream_choice < 0
      ? 0 : solver->eg_stream_choice;
    const double ordinary = solver->eg_stream_min[stream_side];
    if (ordinary > 0.0) {
      solver->eg_premark_samples[0] = 2;
      solver->eg_premark_min[0] = ordinary;
    }
  }
  /* If no equivalent predecessor sample exists, alternate two warm samples
     per arm.  Both layouts execute the same columns and synchronize at the
     same final barrier; only generation-store placement differs. */
  const int side =
    solver->eg_premark_samples[0] <= solver->eg_premark_samples[1] ? 0 : 1;
  solver->eg_premark_pending = side + 1;
  return side == 1;
}

static void kls_egraph_thread_trial_record(kls_solver *solver,
                                           double seconds) {
  if (solver->eg_cluster_pending) {
    const int side = solver->eg_cluster_pending - 1;
    solver->eg_cluster_pending = 0;
    if (solver->eg_cluster_choice != 0 || side < 0 || side > 2) {
      return;
    }
    if (solver->stats.last_refactor_path != KLS_REFACTOR_PATH_EGRAPH) {
      /* An outer path tournament may run PTS after preparing this EGraph
         arm.  That elapsed time says nothing about either EGraph cut; leave
         the verdict open so the first retained EGraph update can restart a
         like-for-like comparison. */
      return;
    }
    solver->eg_cluster_samples[side]++;
    solver->eg_cluster_sum[side] += seconds;
    const int private_alpha4 =
      solver->eg_separator_choice > 0 &&
      solver->refactor_separator_private_cols_alpha4 != NULL &&
      solver->refactor_separator_private_thread_ptr_alpha4 != NULL;
    if (private_alpha4 && solver->eg_cluster_samples[0] == 1 &&
        solver->eg_cluster_samples[2] == 1) {
      const double alpha2_first = solver->eg_cluster_sum[0];
      const double alpha4_first = solver->eg_cluster_sum[2];
      if (alpha2_first > 0.0 && alpha4_first > 0.0 &&
          alpha4_first >= 0.98 * alpha2_first) {
        solver->eg_cluster_choice = -1;
        kls_free_refactor_separator_alpha4_plan(solver);
        return;
      }
    }
    const int enough_samples = private_alpha4
      ? solver->eg_cluster_samples[0] >= 2 &&
        solver->eg_cluster_samples[2] >= 2
      : solver->eg_cluster_samples[0] >= 3 &&
        solver->eg_cluster_samples[1] >= 3;
    if (enough_samples) {
      const double alpha2 = solver->eg_cluster_sum[0] /
        (double)solver->eg_cluster_samples[0];
      const double alpha3 = solver->eg_cluster_samples[1] > 0
        ? solver->eg_cluster_sum[1] /
            (double)solver->eg_cluster_samples[1]
        : 0.0;
      const double alpha4 = solver->eg_cluster_samples[2] > 0
        ? solver->eg_cluster_sum[2] /
            (double)solver->eg_cluster_samples[2]
        : 0.0;
      solver->eg_cluster_choice = -1;
      double best = alpha2;
      if (!private_alpha4 && alpha2 > 0.0 && alpha3 > 0.0 &&
          alpha3 < 0.998 * best) {
        solver->eg_cluster_choice = 1;
        best = alpha3;
      }
      if (private_alpha4 && alpha2 > 0.0 && alpha4 > 0.0 &&
          alpha4 < 0.99 * best) {
        solver->eg_cluster_choice = 2;
      }
      if (private_alpha4 && solver->eg_cluster_choice != 2) {
        kls_free_refactor_separator_alpha4_plan(solver);
      }
    }
    return;
  }
  if (solver->eg_premark_pending) {
    const int side = solver->eg_premark_pending - 1;
    solver->eg_premark_pending = 0;
    if (solver->eg_premark_choice != 0 || side < 0 || side > 1) {
      return;
    }
    if (solver->stats.last_refactor_path != KLS_REFACTOR_PATH_EGRAPH) {
      solver->eg_premark_choice = -1;
      return;
    }
    solver->eg_premark_samples[side]++;
    if (solver->eg_premark_min[side] <= 0.0 ||
        seconds < solver->eg_premark_min[side]) {
      solver->eg_premark_min[side] = seconds;
    }
    if (solver->eg_premark_samples[0] >= 2 &&
        solver->eg_premark_samples[1] >= 2) {
      const double ordinary = solver->eg_premark_min[0];
      const double premarked = solver->eg_premark_min[1];
      solver->eg_premark_choice =
        ordinary > 0.0 && premarked > 0.0 &&
            premarked < 0.995 * ordinary
          ? 1 : -1;
    }
    return;
  }
  if (solver->eg_separator_pending) {
    const int side = solver->eg_separator_pending - 1;
    solver->eg_separator_pending = 0;
    if (solver->eg_separator_choice != 0 || side < 0 || side > 1) {
      return;
    }
    if (solver->stats.last_refactor_path != KLS_REFACTOR_PATH_EGRAPH) {
      /* The EGraph executor can be prepared speculatively before an outer
         PTS/EGraph path trial chooses which implementation supplied this
         update.  Do not turn a sample from the other implementation into a
         permanent rejection of separator-private EGraph execution. */
      return;
    }
    solver->eg_separator_samples[side]++;
    if (side == 1 && solver->eg_separator_samples[side] == 1) {
      return;
    }
    if (solver->eg_separator_min[side] <= 0.0 ||
        seconds < solver->eg_separator_min[side]) {
      solver->eg_separator_min[side] = seconds;
    }
    if (solver->eg_separator_samples[0] >= 1 &&
        solver->eg_separator_samples[1] >= 2) {
      const double ordinary = solver->eg_separator_min[0];
      const double private_domains = solver->eg_separator_min[1];
      solver->eg_separator_choice =
        ordinary > 0.0 && private_domains > 0.0 &&
            private_domains < 0.98 * ordinary
          ? 1 : -1;
    }
    return;
  }
  if (solver->eg_subset_pending) {
    const int side = solver->eg_subset_pending - 1;
    solver->eg_subset_pending = 0;
    if (solver->eg_subset_choice != 0 || side != 1 ||
        solver->stats.last_refactor_path != KLS_REFACTOR_PATH_EGRAPH) {
      return;
    }
    solver->eg_subset_samples[side]++;
    if (solver->eg_subset_min[side] <= 0.0 ||
        seconds < solver->eg_subset_min[side]) {
      solver->eg_subset_min[side] = seconds;
    }
    if (solver->eg_subset_samples[side] == 1 &&
        solver->eg_subset_min[0] > 0.0 &&
        solver->eg_subset_min[1] >= 0.995 * solver->eg_subset_min[0]) {
      solver->eg_subset_choice = -1;
      return;
    }
    if (solver->eg_subset_samples[side] >= 2) {
      const double ordinary = solver->eg_subset_min[0];
      const double subsets = solver->eg_subset_min[1];
      solver->eg_subset_choice =
        ordinary > 0.0 && subsets > 0.0 && subsets < 0.995 * ordinary
          ? 1 : -1;
    }
    return;
  }
  if (solver->eg_stream_pending) {
    const int probe = solver->eg_stream_pending;
    solver->eg_stream_pending = 0;
    if (solver->eg_stream_choice != 0) {
      return;
    }
    if (solver->stats.last_refactor_path != KLS_REFACTOR_PATH_EGRAPH) {
      return;
    }
    if (probe < 1 || probe > 3) {
      return;
    }
    solver->eg_stream_min[probe - 1] = seconds;
    if (probe == 1) {
      return;
    }
    if (probe == 2 && solver->refactor_l_index32_enabled == 0) {
      /* Without a compact indexed stream the two wide-tail arms are the
         same implementation; do not spend a numeric generation measuring
         an exact duplicate. */
      solver->eg_stream_min[2] = seconds;
    } else if (probe < 3) {
      return;
    }
    const double standard = solver->eg_stream_min[0];
    const double scalar_wide = solver->eg_stream_min[1];
    const double vector_wide = solver->eg_stream_min[2];
    /* These equivalent arms isolate tail tiling from indexed scatter.  The
       old bundled comparison could not retain wide producer tiles on a
       factor whose irregular destination stream favors vector scatter. */
    solver->eg_stream_choice = -1;
    double best = standard;
    if (standard > 0.0 && scalar_wide > 0.0 &&
        scalar_wide < 0.99 * standard) {
      solver->eg_stream_choice = 1;
      best = scalar_wide;
    }
    if (standard > 0.0 && vector_wide > 0.0 &&
        vector_wide < 0.99 * standard &&
        (solver->eg_stream_choice < 0 || vector_wide < best)) {
      solver->eg_stream_choice = 2;
    }
    return;
  }
  if (solver->eg_pair_pending) {
    const int probe = solver->eg_pair_pending;
    solver->eg_pair_pending = 0;
    if (solver->eg_pair_choice == 0) {
      if (solver->stats.last_refactor_path != KLS_REFACTOR_PATH_EGRAPH) {
        solver->eg_pair_choice = -1;
        return;
      }
      solver->eg_fuse_min[probe - 1] = seconds;
      if (probe == 1) {
        const double none = solver->eg_tt_min;
        const double pairm = solver->eg_fuse_min[0];
        const double dominated_pair_ratio =
          kls_egraph_pipeline_fusion_trial_worthwhile(solver) ? 1.01 :
          1.10;
        if (none > 0.0 && pairm > dominated_pair_ratio * none) {
          /* Pair fusion and quad fusion share the same cross-column
             synchronization and extra SPAs.  For the tail-only consultation
             even a measured 1% pair loss is enough to avoid spending another
             recurring update on a wider form of the same mechanism; the
             broader legacy tournament keeps its conservative noise margin.
             Settle from realized timings and release the extra scratch. */
          solver->eg_pair_choice = -1;
          kls_egraph_release_unused_quad_scratch(solver);
          return;
        }
        /* Quad probes next.  Both plain single-block and plain BTF kernels
           have dependency-safe k-way implementations; other kernels retain
           pair dispatch, and the strict better-than-pair requirement below
           prevents a noise draw from adopting a nonexistent wider path. */
        solver->eg_pair_pending = 2;
        return;
      }
      const double none = solver->eg_tt_min;
      const double pairm = solver->eg_fuse_min[0];
      const double quadm = solver->eg_fuse_min[1];
      const double fusion_margin =
        kls_egraph_pipeline_fusion_trial_worthwhile(solver) ? 0.99 : 0.95;
      int choice = -1;
      if (none > 0.0 && pairm > 0.0 && pairm < fusion_margin * none) {
        choice = 1;
      }
      if (none > 0.0 && quadm > 0.0 && quadm < fusion_margin * none &&
          (pairm <= 0.0 || quadm < fusion_margin * pairm)) {
        choice = 2;
      }
      solver->eg_pair_choice = choice;
      if (choice != 2) {
        /* Wider fusion is a timed option, not a permanent memory tax.  The
           pool is quiescent when the completed refactor reaches this
           recorder, so rejected extra SPAs can be released immediately. */
        kls_egraph_release_unused_quad_scratch(solver);
      }
    }
    return;
  }
  const int pending = solver->eg_tt_pending;
  solver->eg_tt_pending = 0;
  if (!pending || solver->eg_tt_choice <= 0 ||
      solver->stats.last_refactor_path != KLS_REFACTOR_PATH_EGRAPH) {
    return;
  }
  solver->eg_tt_samples++;
  if (solver->eg_tt_samples == 1) {
    /* Do not score the first-touch warm-up as the fusion baseline. */
    return;
  }
  if (solver->eg_tt_min <= 0.0 || seconds < solver->eg_tt_min) {
    solver->eg_tt_min = seconds;
  }
}

/* Select between two equivalent implementations of the realized indexed
   update streams.  This consultation starts only after dispatch width and
   fusion have settled, so its samples differ solely in scatter/tail kernels.
   Applicability comes from the retained representation (32-bit L streams and
   actual supernode runs), and a strict measured margin protects numerics that
   do not benefit from the longer scalar/panel loops. */
static void kls_egraph_select_stream_kernels(kls_solver *solver,
                                             int check_pivots) {
  solver->eg_stream_pending = 0;
  solver->scalar_refactor_scatter = 0;
  solver->snode_wide_tail = 0;

  if (check_pivots) {
    return;
  }

  const int indexed_scatter_applicable =
    solver->refactor_l_index32_enabled != 0;
  const int wide_tail_applicable =
    kls_refactor_snode_run_end(solver) != NULL;
  if (!indexed_scatter_applicable && !wide_tail_applicable) {
    solver->eg_stream_choice = -1;
    return;
  }
  /* Do not confound the implementation comparison with an outstanding width
     or fusion arm. */
  if (solver->eg_separator_choice == 0 || solver->eg_tt_choice <= 0 ||
      solver->eg_pair_choice == 0 ||
      solver->eg_subset_choice == 0 || solver->eg_subset_pending ||
      solver->eg_tt_pending || solver->eg_pair_pending) {
    return;
  }

  int use_scalar_scatter = solver->eg_stream_choice == 1;
  int use_wide_tail = solver->eg_stream_choice > 0;
  if (solver->eg_stream_choice == 0) {
    const int probe = solver->eg_stream_min[0] <= 0.0 ? 1
      : solver->eg_stream_min[1] <= 0.0 ? 2 : 3;
    solver->eg_stream_pending = probe;
    use_scalar_scatter = probe == 2;
    use_wide_tail = probe >= 2;
  }
  if (use_wide_tail) {
    solver->scalar_refactor_scatter =
      use_scalar_scatter && indexed_scatter_applicable;
    solver->snode_wide_tail = wide_tail_applicable;
  }
}

static void kls_egraph_select_subset_fusion(kls_solver *solver,
                                             int check_pivots) {
  solver->eg_subset_pending = 0;
  if (check_pivots || !kls_repeated_update_workload(&solver->options) ||
      solver->eg_pair_choice != 2) {
    if (solver->eg_pair_choice != 0) {
      solver->eg_subset_choice = -1;
    }
    return;
  }
  if (solver->eg_subset_choice != 0 || solver->eg_pair_pending ||
      solver->eg_tt_pending || solver->eg_separator_pending) {
    return;
  }
  if (solver->eg_subset_samples[0] == 0 && solver->eg_fuse_min[1] > 0.0) {
    solver->eg_subset_samples[0] = 1;
    solver->eg_subset_min[0] = solver->eg_fuse_min[1];
  }
  if (solver->eg_subset_samples[1] < 2) {
    solver->eg_subset_pending = 2;
  }
}

static int kls_egraph_mapped_refactor(kls_solver *solver,
                                      double *numeric_values,
                                      int check_pivots) {
  double **native_l_values = NULL;
  double **native_u_values = NULL;
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      numeric_values == NULL || solver->options.threads <= 1 ||
      solver->n < KLS_EGRAPH_REFACTOR_MIN_SIZE) {
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
  if (kls_pts_refactor_ready(solver)) {
    (void)kls_build_refactor_user_input_pos32(solver);
  }

  int thread_count = solver->options.threads;
  if ((UF_long)thread_count > solver->refactor_level_max_width) {
    thread_count = (int)solver->refactor_level_max_width;
  }
  if (thread_count < 2) {
    return -1;
  }
  const int separator_private =
    kls_egraph_separator_private_dispatch(solver, check_pivots);
  thread_count =
    kls_egraph_steady_thread_count(solver, thread_count, check_pivots);
  kls_egraph_select_subset_fusion(solver, check_pivots);
  kls_egraph_select_stream_kernels(solver, check_pivots);

  const int btf_scalar_run_exec_requested =
    kls_refactor_btf_scalar_run_exec_enabled(solver);

  trilinos_klu_l_common *common = &solver->common;
  const int parallel_scale_requested =
    kls_repeated_update_workload(&solver->options) &&
    solver->n >= (UF_long)(4096 * thread_count) &&
    solver->nnz >= (UF_long)(32768 * thread_count);
  const int parallel_scale_rows =
    common->scale > 0 &&
    parallel_scale_requested &&
    kls_build_refactor_scale_rows(solver);
  const int parallel_permute_requested =
    kls_repeated_update_workload(&solver->options) &&
    solver->n >= (UF_long)(4096 * thread_count);
  int parallel_scale_permute =
    common->scale > 0 && parallel_permute_requested &&
    solver->n <= (UF_long)(SIZE_MAX / sizeof(double));
  if (parallel_scale_permute &&
      solver->refactor_scale_permute_values == NULL) {
    solver->refactor_scale_permute_values = (double *)malloc(
      (size_t)solver->n * sizeof(*solver->refactor_scale_permute_values));
    if (solver->refactor_scale_permute_values == NULL) {
      parallel_scale_permute = 0;
    }
  }
  if (common->scale > 0 &&
      !parallel_scale_rows &&
      !trilinos_klu_l_scale((UF_long)common->scale, solver->n,
                            solver->col_ptr, solver->row_idx,
                            numeric_values, solver->numeric->Rs, NULL,
                            common)) {
    return 0;
  }
  const kls_egraph_refactor_kernel selected_kernel =
    kls_egraph_refactor_kernel_for(solver, (int)common->scale);

  const int build_separator_private =
    kls_egraph_separator_private_requested(solver);
  const UF_long cluster_level_count =
    kls_egraph_cluster_level_dispatch(solver, check_pivots);
  atomic_uint *pipeline_done = NULL;
  unsigned int pipeline_generation = 0;

  if (solver->refactor_cluster_level_count < solver->refactor_level_count &&
       solver->refactor_pipeline_column_count >= (UF_long)(2 * thread_count) &&
       solver->refactor_pipeline_work >=
         0.10 * solver->refactor_dependency_work) {
    pipeline_done =
      ensure_egraph_pipeline_done(solver, &pipeline_generation);
  }
  if (build_separator_private &&
      solver->refactor_separator_private_cols != NULL &&
      solver->refactor_separator_private_thread_ptr != NULL &&
      solver->refactor_separator_private_thread_count == thread_count &&
      pipeline_done == NULL) {
    pipeline_done =
      ensure_egraph_pipeline_done(solver, &pipeline_generation);
  }
  if (check_pivots && pipeline_done == NULL) {
    pipeline_done =
      ensure_egraph_pipeline_done(solver, &pipeline_generation);
  }
  if (check_pivots && pipeline_done == NULL) {
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
  if (
      getenv("KLS_EGRAPH_WORKER_SPIN_OVERRIDE") == NULL &&
      !check_pivots && kls_repeated_update_workload(&solver->options) &&
      kls_pts_solve_available(solver) &&
      pool->worker_spin_iters < 200000u) {
    /* A retained EGraph crew is reused by the following solve/value phases
       and by the next numeric generation.  The generic pool default is
       deliberately short for isolated calls, but it expires before these
       back-to-back repeated-update phases and turns every generation into
       N-1 condition-variable wakeups.  Once the realized executor and the
       caller's repeated lifecycle establish reuse, keep the crew live over
       that bounded gap only when the verified parallel triangular solve will
       immediately dispatch the same pool again.  This is an implementation-
       capability decision: serial solves plus row, lean, and direct numeric
       executors retain their independently selected idle policies, and the
       public override remains authoritative. */
    pool->worker_spin_iters = 200000u;
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
  shared->parallel_scale_rows = parallel_scale_rows;
  shared->parallel_scale_method = (int)common->scale;
  shared->parallel_scale_permute = parallel_scale_permute;
  shared->parallel_scale_input_values = numeric_values;
  shared->kernel = selected_kernel;
  shared->thread_count = thread_count;
  /* A resident spin barrier wins when synchronization is frequent relative
     to the remaining dependency-pipeline work.  On long pipeline tails,
     allowing pthread_barrier_wait to stagger/sleep early arrivals avoids
     cache and bandwidth contention.  Work is the ordering-independent
     symbolic estimate built with the dependency graph, not a matrix- or
     machine-specific timing sample. */
  shared->use_spin_cluster_barrier =
    cluster_level_count > 0u &&
    solver->refactor_pipeline_work <=
      1500000.0 * (double)cluster_level_count;
  /* Generations are local to one numeric dispatch.  The crew is idle while
     shared->lock is held here, so clear every private arrival before worker
     zero starts comparing the new generation against those slots. */
  for (int tid = 0; tid < thread_count; ++tid) {
    atomic_store_explicit(&pool->workers[tid].cluster_barrier_arrival, 0u,
                          memory_order_relaxed);
  }
  atomic_store_explicit(&shared->cluster_barrier_generation, 0u,
                        memory_order_relaxed);
  shared->snode_run_end = kls_refactor_snode_run_end(solver);
  shared->deferred_pts_terminal_col = KLS_KLU_EMPTY;
  if (!check_pivots &&
      selected_kernel == KLS_EGRAPH_REFACTOR_KERNEL_BTF_UNSCALED &&
      shared->btf_unscaled_plain &&
      kls_repeated_update_workload(&solver->options) &&
      solver->refactor_u_indices != NULL) {
    shared->deferred_pts_terminal_col =
      kls_pts_egraph_terminal_candidate(solver);
    if (shared->deferred_pts_terminal_col != KLS_KLU_EMPTY) {
      shared->kernel =
        KLS_EGRAPH_REFACTOR_KERNEL_BTF_UNSCALED_DEFER_TERMINAL;
    }
  }
  shared->row_refactor_mode = 0;
  shared->row_solve_mode = 0;
  shared->row_solve_view = NULL;
  shared->row_solve_work = NULL;
  shared->row_refactor_native_row_panel_active = 0;
  shared->row_refactor_defer_value_scatter = 0;
  shared->row_pipeline_ready_queue = 0;
  shared->row_pipeline_row_dep_ready_queue = 0;
  shared->row_pipeline_ready_groups = NULL;
  shared->row_pipeline_ready_slots = NULL;
  shared->row_pipeline_remaining_preds = NULL;
  shared->row_pipeline_tail_groups = NULL;
  shared->row_pipeline_tail_count = 0;
  shared->row_pipeline_private_groups = NULL;
  shared->row_pipeline_private_thread_ptr = NULL;
  shared->row_pipeline_private_group_mask = NULL;
  shared->row_pipeline_private_count = 0;
  shared->row_pipeline_ordered_private = 0;
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
  shared->cluster_level_count =
    pipeline_done != NULL ? cluster_level_count : solver->refactor_level_count;
  shared->separator_private = separator_private;
  shared->btf_scalar_run_exec =
    (btf_scalar_run_exec_requested &&
     kls_egraph_kernel_is_btf_unscaled(shared->kernel))
      ? 1 : 0;

  shared->single_unscaled_plain =
    selected_kernel == KLS_EGRAPH_REFACTOR_KERNEL_SINGLE_UNSCALED &&
    !shared->check_pivots && !solver->padded_active &&
    shared->snode_run_end != NULL &&
    solver->eg_pair_choice < 0 && !solver->eg_pair_pending &&
    !solver->eg_stream_pending;
  /* All inputs to the plain-BTF capability test are immutable for one
     refactor dispatch.  Cache it beside the corresponding single-block
     verdict instead of replaying the long advanced-surface predicate for
     every numeric column. */
  shared->btf_unscaled_plain =
    selected_kernel == KLS_EGRAPH_REFACTOR_KERNEL_BTF_UNSCALED &&
    kls_egraph_btf_plain_cluster_allowed(shared);

  atomic_store_explicit(&shared->btf_scalar_run_exec_runs, 0ul,
                        memory_order_release);
  atomic_store_explicit(&shared->btf_scalar_run_exec_rows, 0ul,
                        memory_order_release);
  atomic_store_explicit(&shared->btf_scalar_run_exec_entries, 0ul,
                        memory_order_release);
  atomic_store_explicit(&shared->btf_scalar_run_exec_max_rows, 0ul,
                        memory_order_release);

  shared->cluster_done_premarked = 0;
  if ((kls_egraph_cluster_premark_dispatch(solver, shared, check_pivots) &&
      !check_pivots && pipeline_done != NULL && pipeline_generation != 0u) && cluster_level_count > 0u &&
      cluster_level_count <= solver->refactor_level_count &&
      solver->refactor_level_ptr != NULL &&
      solver->refactor_level_cols != NULL) {
    /* Every worker crosses the last cluster-level barrier before it can
       enter the dependency-waiting pipeline.  Seed the generation bitmap in
       one pre-dispatch pass instead of interleaving its stores with the hot
       numeric columns; the barrier remains the publication boundary for the
       numeric values themselves. */
    const UF_long cluster_end =
      solver->refactor_level_ptr[cluster_level_count];
    if (cluster_end <= solver->n) {
      for (UF_long pos = 0u; pos < cluster_end; ++pos) {
        const UF_long col = solver->refactor_level_cols[pos];
        if (col >= solver->n) {
          shared->invalid = 1;
          break;
        }
        atomic_store_explicit(&pipeline_done[col], pipeline_generation,
                              memory_order_relaxed);
      }
      shared->cluster_done_premarked = !shared->invalid;
    }
  }
  atomic_store_explicit(
    &shared->next_pipeline_pos,
    (unsigned long)(pipeline_done != NULL
                      ? solver->refactor_level_ptr[cluster_level_count]
                      : solver->n),
    memory_order_release);
  shared->pipeline_pos_end = solver->n;
  shared->cluster_level_count =
    pipeline_done != NULL ? cluster_level_count : solver->refactor_level_count;
  shared->udiag_recip =
    selected_kernel == KLS_EGRAPH_REFACTOR_KERNEL_BTF_UNSCALED &&
      solver->dense_tail_cols == 0u && solver->n >= 262144u &&
      solver->i32solve_state > 0 &&
      solver->i32solve_udiag_recip != NULL &&
      shared->btf_unscaled_plain
    ? solver->i32solve_udiag_recip : NULL;
  const int parallel_refine_copy =
    !solver->parallel_refine_values_copied &&
    solver->solve_refine_values != NULL &&
    solver->solve_refine_values != numeric_values &&
    solver->nnz >= (UF_long)(16384 * thread_count);
  shared->refine_copy_input =
    parallel_refine_copy ? numeric_values : NULL;
  shared->refine_copy_output =
    parallel_refine_copy ? solver->solve_refine_values : NULL;
  shared->refine_copy_nnz = parallel_refine_copy ? solver->nnz : 0u;

  if (kls_authoritative_packed_lu_enabled(solver) &&
      kls_ensure_authoritative_packed_l_values(solver) &&
      kls_ensure_authoritative_packed_u_values(solver)) {
    native_l_values = solver->refactor_l_values;
    native_u_values = solver->refactor_u_values;
    solver->refactor_l_values = solver->refactor_l_packed_values;
    solver->refactor_u_values = solver->refactor_u_packed_values;
  }

  for (int i = 0; i < thread_count; ++i) {
    pool->workers[i].shared = shared;
    pool->workers[i].x = scratch[i];

  }

  kls_egraph_pool_dispatch_and_wait(pool, shared, thread_count);
  if (native_l_values != NULL) {
    solver->refactor_l_values = native_l_values;
    solver->refactor_u_values = native_u_values;
    solver->refactor_l_packed_valid = 1;
    solver->refactor_u_packed_valid = 1;
    if (shared->invalid || shared->pivot_rejected || shared->singular) {
      kls_sync_authoritative_packed_l_values(solver);
      kls_sync_authoritative_packed_u_values(solver);
      solver->refactor_l_packed_valid = 0;
      solver->refactor_u_packed_valid = 0;
    }
  }
  if (!shared->invalid && !shared->pivot_rejected && !shared->singular &&
      shared->deferred_pts_terminal_col != KLS_KLU_EMPTY) {
    const UF_long deferred = shared->deferred_pts_terminal_col;
    pthread_mutex_unlock(&shared->lock);
    const int pts_ok = kls_pts_recompute_egraph_terminal(
      solver, numeric_values, scratch[0]);
    pthread_mutex_lock(&shared->lock);
    if (!pts_ok) {
      /* Eligibility is conservative, but allocation/pool availability is a
         runtime property.  Preserve exact behavior by executing the original
         leaf kernel after all of its dependencies have completed. */
      shared->deferred_pts_terminal_col = KLS_KLU_EMPTY;
      memset(scratch[0], 0,
             (size_t)solver->symbolic->maxblock * sizeof(*scratch[0]));
      double **saved_l_values = solver->refactor_l_values;
      double **saved_u_values = solver->refactor_u_values;
      if (solver->refactor_l_packed_valid &&
          solver->refactor_u_packed_valid) {
        solver->refactor_l_values = solver->refactor_l_packed_values;
        solver->refactor_u_values = solver->refactor_u_packed_values;
      }
      if (!kls_egraph_refactor_dispatch_column(&pool->workers[0], deferred,
                                                0)) {
        shared->invalid = 1;
      }
      solver->refactor_l_values = saved_l_values;
      solver->refactor_u_values = saved_u_values;
    } else if (shared->udiag_recip != NULL) {
      const double pivot = ((double *)solver->numeric->Udiag)[deferred];
      shared->udiag_recip[deferred] = 1.0 / pivot;
    }
  }
  if (!shared->invalid && !shared->pivot_rejected && !shared->singular &&
      parallel_refine_copy) {
    solver->parallel_refine_values_copied = 1;
  }
  shared->refine_copy_input = NULL;
  shared->refine_copy_output = NULL;
  shared->refine_copy_nnz = 0u;
  if (!shared->invalid && !shared->pivot_rejected && !shared->singular &&
      shared->udiag_recip != NULL) {
    /* Every BTF column published its reciprocal alongside Udiag while the
       owning numeric worker already had both values hot. */
    solver->i32solve_udiag_recip_fresh = 1;
  }
  shared->udiag_recip = NULL;

  shared->separator_private = 0;

  shared->btf_scalar_run_exec = 0;

  pthread_mutex_unlock(&shared->lock);

  const UF_long btf_scalar_run_exec_runs =
    (UF_long)atomic_load_explicit(&shared->btf_scalar_run_exec_runs,
                                  memory_order_acquire);
  const UF_long btf_scalar_run_exec_rows =
    (UF_long)atomic_load_explicit(&shared->btf_scalar_run_exec_rows,
                                  memory_order_acquire);
  const UF_long btf_scalar_run_exec_entries =
    (UF_long)atomic_load_explicit(&shared->btf_scalar_run_exec_entries,
                                  memory_order_acquire);
  const UF_long btf_scalar_run_exec_max_rows =
    (UF_long)atomic_load_explicit(&shared->btf_scalar_run_exec_max_rows,
                                  memory_order_acquire);

  solver->refactor_last_btf_scalar_run_exec_runs =
    btf_scalar_run_exec_runs;
  solver->refactor_last_btf_scalar_run_exec_rows = btf_scalar_run_exec_rows;
  solver->refactor_last_btf_scalar_run_exec_entries =
    btf_scalar_run_exec_entries;
  solver->refactor_last_btf_scalar_run_exec_max_rows =
    btf_scalar_run_exec_max_rows;
  solver->refactor_btf_scalar_run_exec_count += btf_scalar_run_exec_runs;
  solver->refactor_btf_scalar_run_exec_rows += btf_scalar_run_exec_rows;
  solver->refactor_btf_scalar_run_exec_entries += btf_scalar_run_exec_entries;
  if (btf_scalar_run_exec_max_rows >
      solver->refactor_btf_scalar_run_exec_max_rows) {
    solver->refactor_btf_scalar_run_exec_max_rows =
      btf_scalar_run_exec_max_rows;
  }

  if (shared->invalid || shared->pivot_rejected ||
      (shared->singular && common->halt_if_singular)) {
    solver->egraph_worker_scratch_dirty = 1;
  }

  int unfinished_tail_seed_recorded = 0;
  if (!shared->invalid && shared->pivot_rejected &&
      pipeline_done != NULL) {
    kls_record_fast_reject_unfinished_tail_seed(
      solver, shared, shared->rejected_pivot);
    unfinished_tail_seed_recorded =
      solver->fast_reject_tail_seed_valid ? 1 : 0;
  }
  if (!shared->invalid && shared->pivot_rejected &&
      pipeline_done != NULL && scratch[0] != NULL) {
    kls_egraph_refactor_refresh_missing_prefix(shared, scratch[0],
                                               scratch_size);
  }

  if (shared->invalid) {
    kls_record_fast_factor_failure(
      solver, KLS_FAST_FACTOR_FAIL_EGRAPH_INVALID, TRILINOS_KLU_INVALID);
    common->status = TRILINOS_KLU_INVALID;
    return 0;
  }
  if (shared->pivot_rejected) {
    const int prefix_ready =
      kls_egraph_refreshed_prefix(shared, shared->rejected_pivot);
    solver->fast_reject_refresh_state =
      prefix_ready ? KLS_FAST_REJECT_REFRESH_PREFIX
                   : KLS_FAST_REJECT_REFRESH_UNKNOWN;
    if (prefix_ready) {
      const UF_long seed_block =
        kls_block_for_pivot(solver, shared->rejected_pivot);
      if (!unfinished_tail_seed_recorded ||
          seed_block == KLS_KLU_EMPTY ||
          solver->fast_reject_tail_seed_block != seed_block) {
        kls_record_fast_reject_unfinished_tail_seed(
          solver, shared, shared->rejected_pivot);
      }
    } else {
      kls_clear_fast_reject_tail_seed(solver);
    }
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
    if (common->scale > 0 && !parallel_scale_permute &&
        !kls_parallel_refactor_permute_scale(solver)) {
      common->status = TRILINOS_KLU_INVALID;
      return 0;
    }
    return 1;
  }
  common->status = TRILINOS_KLU_OK;
  common->numerical_rank = KLS_KLU_EMPTY;
  common->singular_col = KLS_KLU_EMPTY;
  common->nrealloc = 0;
  if (common->scale > 0 && !parallel_scale_permute &&
      !kls_parallel_refactor_permute_scale(solver)) {
    common->status = TRILINOS_KLU_INVALID;
    return 0;
  }
  solver->refactor_lu_pointer_egraph_certified = 1;
  return 1;
}

#ifdef KLS_HAVE_CBLAS
/* No-pivot right-looking blocked LU of the t x t column-major block D.
   The refactorization contract keeps the first factor's pivot order, so
   the trailing dense block refreshes values-only (dgetrf would repivot).
   On return D's strictly-lower part holds the unit-L multipliers and its
   upper part (with diagonal) holds U — the layout the tail columns'
   stored patterns want back. */
static int kls_dense_tail_nopivot_lu(double *D, UF_long t) {
  const UF_long nb = 128;
  for (UF_long j0 = 0; j0 < t; j0 += nb) {
    const UF_long jb = t - j0 < nb ? t - j0 : nb;
    for (UF_long j = j0; j < j0 + jb; ++j) {
      const double piv = D[j * t + j];
      if (piv == 0.0) {
        return 0;
      }
      const double pinv = 1.0 / piv;
      for (UF_long i = j + 1; i < j0 + jb; ++i) {
        D[j * t + i] *= pinv;
      }
      for (UF_long k = j + 1; k < j0 + jb; ++k) {
        const double u = D[k * t + j];
        if (u == 0.0) {
          continue;
        }
        for (UF_long i = j + 1; i < j0 + jb; ++i) {
          D[k * t + i] -= D[j * t + i] * u;
        }
      }
    }
    const UF_long rest = t - j0 - jb;
    if (rest > 0) {
      cblas_dtrsm(CblasColMajor, CblasRight, CblasUpper, CblasNoTrans,
                  CblasNonUnit, (int)rest, (int)jb, 1.0,
                  &D[j0 * t + j0], (int)t,
                  &D[j0 * t + j0 + jb], (int)t);
      cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans,
                  CblasUnit, (int)jb, (int)rest, 1.0,
                  &D[j0 * t + j0], (int)t,
                  &D[(j0 + jb) * t + j0], (int)t);
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  (int)rest, (int)rest, (int)jb, -1.0,
                  &D[j0 * t + j0 + jb], (int)t,
                  &D[(j0 + jb) * t + j0], (int)t, 1.0,
                  &D[(j0 + jb) * t + j0 + jb], (int)t);
    }
  }
  return 1;
}

typedef struct {
  kls_solver *solver;
  const double *numeric_values;
  double *D;
  double *X;      /* this worker's private nk-sized workspace */
  UF_long c0;     /* first tail column (block-local), cyclic by stride */
  UF_long stride;
  UF_long k1;
  UF_long nk;
  UF_long s;
  UF_long t;
  int write_u;    /* production: store prefix-U values into the factor */
  int ok;
} kls_dense_tail_gather_job;

static void *kls_dense_tail_gather_main(void *arg) {
  kls_dense_tail_gather_job *g = (kls_dense_tail_gather_job *)arg;
  kls_solver *solver = g->solver;
  trilinos_klu_l_numeric *num = solver->numeric;
  const trilinos_klu_l_symbolic *sym = solver->symbolic;
  double *lu = (double *)num->LUbx[solver->dense_tail_block];
  UF_long *lip = num->Lip + g->k1;
  UF_long *llen = num->Llen + g->k1;
  UF_long *uip = num->Uip + g->k1;
  UF_long *ulen_arr = num->Ulen + g->k1;
  const double *rs = num->Rs;
  const UF_long *pinv = num->Pinv;
  const UF_long *q = sym->Q;
  double *X = g->X;
  const UF_long s = g->s;
  const UF_long nk = g->nk;
  const UF_long t = g->t;
  g->ok = 0;
  for (UF_long c = g->c0; c < nk; c += g->stride) {
    const UF_long oldcol = q[g->k1 + c];
    for (UF_long p = solver->col_ptr[oldcol];
         p < solver->col_ptr[oldcol + 1]; ++p) {
      const UF_long oldrow = solver->row_idx[p];
      const UF_long newrow = pinv[oldrow] - g->k1;
      if (newrow < 0 || newrow >= nk) {
        continue; /* off-block entry: the vendor walk owns Offx */
      }
      X[newrow] = rs != NULL ? g->numeric_values[p] / rs[g->k1 + newrow]
                             : g->numeric_values[p];
    }
    UF_long *ui;
    double *ux;
    UF_long ul;
    kls_klu_get_pointer(lu, uip, ulen_arr, c, &ui, &ux, &ul);
    for (UF_long up = 0; up < ul; ++up) {
      const UF_long j = ui[up];
      if (j >= s) {
        break; /* ascending patterns: tail entries come from the LU */
      }
      const double ujk = X[j];
      X[j] = 0.0;
      if (g->write_u) {
        ux[up] = ujk;
      }
      if (ujk == 0.0) {
        continue;
      }
      UF_long *li;
      double *lx;
      UF_long ll;
      kls_klu_get_pointer(lu, lip, llen, j, &li, &lx, &ll);
      for (UF_long lp = 0; lp < ll; ++lp) {
        X[li[lp]] -= lx[lp] * ujk;
      }
    }
    for (UF_long r = s; r < nk; ++r) {
      g->D[(size_t)(c - s) * t + (r - s)] = X[r];
      X[r] = 0.0;
    }
  }
  g->ok = 1;
  return NULL;
}

extern void openblas_set_num_threads(int) __attribute__((weak));

/* Values-only refresh of the pipe-emitted dense tail: parallel gather of
   the updated trailing block (prefix-consumption per tail column against
   the just-refreshed prefix), one no-pivot blocked LU, and a writeback
   through the stored dense patterns.  Replaces the 24s scalar scatter
   walk of ss1's 16.8M-entry tail; the vendor refactor ran with the tail
   columns skipped (off-diagonal walk preserved). */
static int kls_dense_tail_refresh(kls_solver *solver,
                                  const double *numeric_values) {
  const trilinos_klu_l_symbolic *sym = solver->symbolic;
  trilinos_klu_l_numeric *num = solver->numeric;
  const UF_long b = solver->dense_tail_block;
  const UF_long k1 = sym->R[b];
  const UF_long nk = sym->R[b + 1] - k1;
  const UF_long t = solver->dense_tail_cols;
  const UF_long s = nk - t;
  if (s <= 0 || num->LUbx == NULL || num->LUbx[b] == NULL) {
    return 0;
  }
  double *D = (double *)malloc((size_t)t * (size_t)t * sizeof(*D));
  if (D == NULL) {
    return 0;
  }
  memset(D, 0, (size_t)t * (size_t)t * sizeof(*D));
  int nthreads = (int)solver->options.threads;
  if (nthreads > 8) {
    nthreads = 8;
  }
  if (nthreads < 1) {
    nthreads = 1;
  }
  kls_dense_tail_gather_job jobs[8];
  pthread_t tids[8];
  int spawned = 0;
  int ok = 1;
  for (int w = 0; w < nthreads; ++w) {
    jobs[w].solver = solver;
    jobs[w].numeric_values = numeric_values;
    jobs[w].D = D;
    jobs[w].X = (double *)calloc((size_t)nk, sizeof(double));
    jobs[w].c0 = s + (UF_long)w;
    jobs[w].stride = (UF_long)nthreads;
    jobs[w].k1 = k1;
    jobs[w].nk = nk;
    jobs[w].s = s;
    jobs[w].t = t;
    jobs[w].write_u = 1;
    jobs[w].ok = 0;
    if (jobs[w].X == NULL) {
      ok = 0;
      break;
    }
  }
  if (ok) {
    for (int w = 1; w < nthreads; ++w) {
      if (pthread_create(&tids[w], NULL, kls_dense_tail_gather_main,
                         &jobs[w]) != 0) {
        ok = 0;
        break;
      }
      spawned = w;
    }
    if (ok) {
      kls_dense_tail_gather_main(&jobs[0]);
    }
    for (int w = 1; w <= spawned; ++w) {
      pthread_join(tids[w], NULL);
    }
    for (int w = 0; w < nthreads && ok; ++w) {
      ok = ok && jobs[w].ok;
    }
  }
  for (int w = 0; w < nthreads; ++w) {
    free(jobs[w].X);
  }
  if (ok) {
    if (openblas_set_num_threads != NULL) {
      openblas_set_num_threads(nthreads);
    }
    ok = kls_dense_tail_nopivot_lu(D, t);
    if (openblas_set_num_threads != NULL) {
      openblas_set_num_threads(1);
    }
  }
  if (ok) {
    double *lu = (double *)num->LUbx[b];
    UF_long *lip = num->Lip + k1;
    UF_long *llen = num->Llen + k1;
    UF_long *uip = num->Uip + k1;
    UF_long *ulen_arr = num->Ulen + k1;
    double *udiag = ((double *)num->Udiag) + k1;
    for (UF_long c = s; c < nk; ++c) {
      UF_long *ui;
      double *ux;
      UF_long ul;
      kls_klu_get_pointer(lu, uip, ulen_arr, c, &ui, &ux, &ul);
      for (UF_long up = 0; up < ul; ++up) {
        const UF_long j = ui[up];
        if (j < s) {
          continue; /* prefix values were stored during the gather */
        }
        ux[up] = D[(size_t)(c - s) * t + (j - s)];
      }
      udiag[c] = D[(size_t)(c - s) * t + (c - s)];
      UF_long *li;
      double *lx;
      UF_long ll;
      kls_klu_get_pointer(lu, lip, llen, c, &li, &lx, &ll);
      for (UF_long lp = 0; lp < ll; ++lp) {
        const UF_long r = li[lp];
        if (r < s) {
          continue;
        }
        lx[lp] = D[(size_t)(c - s) * t + (r - s)];
      }
    }
  }
  free(D);
  return ok;
}

/* Every serial klu refactorization of a dense-tail numeric goes through
   here: the vendor walk runs with the tail columns skipped (their
   off-diagonal entries still consumed), then the BLAS3 refresh rebuilds
   the trailing block values-only.  Failure falls back to one full
   scalar walk for a consistent factor. */
static UF_long kls_serial_klu_refactor_dense_tail(kls_solver *solver,
                                                  double *numeric_values,
                                                  int check_pivots) {
  const int dense_tail_route =
    solver->dense_tail_cols > 0 && !check_pivots;
  if (dense_tail_route) {
    kls_klu_refactor_tail_skip = (long)solver->dense_tail_cols;
    kls_klu_refactor_tail_block = (long)solver->dense_tail_block;
  }
  UF_long serial_ok =
    trilinos_klu_l_refactor(solver->col_ptr, solver->row_idx,
                            numeric_values, solver->symbolic,
                            solver->numeric, &solver->common);
  if (dense_tail_route) {
    kls_klu_refactor_tail_skip = 0;
    kls_klu_refactor_tail_block = -1;
    if (serial_ok && solver->common.status >= 0 &&
        !kls_dense_tail_refresh(solver, numeric_values)) {
      /* the BLAS3 refresh failed (zero pivot / allocation): restore a
         consistent factor with one full scalar walk */
      serial_ok =
        trilinos_klu_l_refactor(solver->col_ptr, solver->row_idx,
                                numeric_values, solver->symbolic,
                                solver->numeric, &solver->common);
    }
  }
  return serial_ok;
}

#else /* !KLS_HAVE_CBLAS */
static UF_long kls_serial_klu_refactor_dense_tail(kls_solver *solver,
                                                  double *numeric_values,
                                                  int check_pivots) {
  (void)check_pivots;
  return trilinos_klu_l_refactor(solver->col_ptr, solver->row_idx,
                                 numeric_values, solver->symbolic,
                                 solver->numeric, &solver->common);
}
#endif /* KLS_HAVE_CBLAS */

static int kls_serial_refactor_tail_from_block(kls_solver *solver,
                                               double *numeric_values,
                                               UF_long start_block,
                                               int check_pivots);

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
  const int moderate_btf_trial = !check_pivots && !scaled &&
    kls_moderate_btf_map32_trial_capable(solver);
  if (kls_pts_refactor_ready(solver) ||
      kls_low_work_btf_map32_policy_enabled(solver) ||
      kls_symmetric_partial_diagonal_match_factor_cycle(solver)) {
    (void)kls_build_refactor_user_input_pos32(solver);
    if (kls_low_work_btf_map32_policy_enabled(solver)) {
      /* The first changed-value pass uses the ordinary prepared frame, but
         build the direct Offx map while its column map is already hot.  Once
         the numeric contract is established, later updates can consume the
         caller's entrywise values without another permutation copy. */
      (void)kls_build_lean_btf_off_map(solver, 1);
    }
  }
  if (!check_pivots && !scaled &&
      kls_low_work_single_block_policy_enabled(solver)) {
    const int lean_single =
      kls_lean_single_block_map32_refactor(solver, numeric_values);
    if (lean_single >= 0) {
      return lean_single;
    }
  }
  if (!check_pivots && !scaled) {
    if (moderate_btf_trial && solver->moderate_btf_lean_choice == 0) {
      /* Build the compact Offx deck before timing either executor.  Charging
         its one-time construction to the lean arm made the old single-sample
         comparison reject a kernel that was consistently faster after the
         retained map became warm.  The setup remains charged to the public
         refactor call; it is merely excluded from the steady executor
         comparison used to settle the recurring representation. */
      if (!kls_build_lean_btf_off_map(solver, 0)) {
        solver->moderate_btf_lean_choice = -1;
      } else {
        /* Compare one complete execution of each exact representation on the
           same changed values.  The former warm-pair portfolio performed
           four numerics here and often spent more than an H100 win could
           recover.  A compact arm now needs a decisive single-pair margin;
           a close result merely admits the stronger row representation.
           Losers are restored once, while either admitted arm is already a
           valid current numeric and needs no redundant mapped refresh. */
        const double trial_start = kls_now_seconds();
        const double mapped_start = kls_now_seconds();
        int mapped = kls_serial_refactor_tail_from_block(
          solver, numeric_values, 0u, 0);
        const double mapped_seconds =
          kls_now_seconds() - mapped_start;
        const double lean_start = kls_now_seconds();
        const int lean = mapped > 0
          ? kls_lean_btf_map32_refactor(solver, numeric_values) : mapped;
        const double lean_seconds = kls_now_seconds() - lean_start;
        const double trial_seconds = kls_now_seconds() - trial_start;
        solver->moderate_btf_mapped_seconds = mapped_seconds;
        solver->moderate_btf_compact_seconds = lean_seconds;
        solver->moderate_btf_trial_seconds = trial_seconds;
        if (mapped < 0 || lean < 0) {
          solver->moderate_btf_lean_choice = -1;
        } else if (mapped == 0 || lean == 0) {
          return 0;
        } else if (lean_seconds < 0.90 * mapped_seconds &&
                   mapped_seconds - lean_seconds >= 2.0e-6 &&
                   (mapped_seconds - lean_seconds) *
                       (double)(solver->options.expected_refactorizations - 1)
                     > 1.25 * trial_seconds) {
          solver->moderate_btf_lean_choice = 1;
          return lean;
        } else {
          /* A compact stream that is competitive with the retained mapped
             walk is useful evidence even when it cannot repay adoption on
             its own.  Admit the existing row-layout tournament, which times
             and validates the stronger representation before selecting it.
             The compact execution is the current valid factor, so an
             admitted row trial returns it directly. */
          solver->moderate_btf_lean_choice =
            lean_seconds <= 1.10 * mapped_seconds ? 2 : -1;
          if (solver->moderate_btf_lean_choice == 2) return lean;
          return kls_serial_refactor_tail_from_block(
            solver, numeric_values, 0u, 0);
        }
      }
    }
    if (moderate_btf_trial && solver->moderate_btf_lean_choice == 1) {
      const int lean = kls_lean_btf_map32_refactor(solver, numeric_values);
      if (lean >= 0) {
        return lean;
      }
      solver->moderate_btf_lean_choice = -1;
    }
    if (kls_low_work_btf_map32_policy_enabled(solver)) {
      const int lean =
        kls_lean_btf_map32_refactor(solver, numeric_values);
      if (lean >= 0) {
        return lean;
      }
    }
  }
  if (solver->symbolic->nblocks == 1u) {
    if (check_pivots && kls_checked_row_refactor_env_enabled()) {
      if (solver->options.threads > 1) {
        const int parallel_row_status =
          kls_threaded_row_refactor_numeric(solver, numeric_values, 1);
        if (parallel_row_status >= 0) {
          return parallel_row_status;
        }
      }
      const int row_status =
        kls_serial_row_refactor_numeric(solver, numeric_values, 1);
      if (row_status >= 0) {
        return row_status;
      }
    } else if (!check_pivots) {
      /* The snb engine's per-matrix acceptance must see refactor calls
         even where the row engine adopted: its timed comparison uses
         the wrapper-recorded incumbent (whichever engine that is) and
         adopts only decisive wins (pre2: snb 3.2s vs row 8.1s). */
      UF_long snb_ok = 0;
      if (kls_snb_try_refactor(solver, numeric_values, 0, &snb_ok)) {
        return (int)snb_ok;
      }
      if (kls_row_refactor_acceptance_wants_row(solver)) {
      if (solver->options.threads > 1) {
        const int parallel_row_status =
          kls_threaded_row_refactor_numeric(solver, numeric_values, 0);
        if (parallel_row_status >= 0) {
          kls_set_last_refactor_path(solver, KLS_REFACTOR_PATH_ROW);
          return parallel_row_status;
        }
        /* structural decline: the attempt's schedule prep alone costs
           ~1.4s on massively-blocked numerics (ASIC_680ks: 583,771
           groups, 69% of a 20-refactor run's CPU re-prepping before
           declining every time) - settle on the column engine for
           this numeric instead of re-attempting */
        solver->row_accept_decision = -1;
      }
      if (kls_row_refactor_env_enabled()) {
        const int row_status =
          kls_serial_row_refactor_numeric(solver, numeric_values, 0);
        if (row_status >= 0) {
          kls_set_last_refactor_path(solver, KLS_REFACTOR_PATH_ROW);
          return row_status;
        }
      } else if (solver->row_accept_decision == 0) {
        /* threaded trial infeasible: settle on the column engine */
        solver->row_accept_decision = -1;
      }
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
  shared.solver = solver;
  shared.native_short_l = kls_mapped_native_short_l_enabled(solver);
  shared.col_ptr = solver->col_ptr;
  shared.row_idx = solver->row_idx;
  shared.nnz = solver->nnz;
  shared.map_col_ptr = solver->refactor_col_ptr;
  shared.map_row_idx = solver->refactor_row_idx;
  shared.map_input_pos = solver->refactor_input_pos;
  shared.map_block_start = solver->refactor_block_start;
  shared.snode_run_end = kls_refactor_snode_run_end(solver);
  shared.padded_src = solver->padded_run_of != NULL ? solver : NULL;
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
    kls_record_fast_factor_failure(
      solver, KLS_FAST_FACTOR_FAIL_MAPPED_INVALID, TRILINOS_KLU_INVALID);
    common->status = TRILINOS_KLU_INVALID;
    return 0;
  }
  if (worker.pivot_rejected) {
    solver->fast_reject_refresh_state = KLS_FAST_REJECT_REFRESH_PREFIX;
    kls_record_fast_reject_block_suffix_tail_seed(solver,
                                                 worker.rejected_pivot);
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
  shared.solver = solver;
  shared.native_short_l = kls_mapped_native_short_l_enabled(solver);
  shared.col_ptr = solver->col_ptr;
  shared.row_idx = solver->row_idx;
  shared.nnz = solver->nnz;
  shared.map_col_ptr = solver->refactor_col_ptr;
  shared.map_row_idx = solver->refactor_row_idx;
  shared.map_input_pos = solver->refactor_input_pos;
  shared.map_block_start = solver->refactor_block_start;
  shared.direct_user_values = solver->refactor_direct_user_values_active;
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
    kls_record_fast_factor_failure(
      solver, KLS_FAST_FACTOR_FAIL_MAPPED_INVALID, TRILINOS_KLU_INVALID);
    common->status = TRILINOS_KLU_INVALID;
    return 0;
  }
  if (worker.pivot_rejected) {
    solver->fast_reject_refresh_state = KLS_FAST_REJECT_REFRESH_PREFIX;
    kls_record_fast_reject_block_suffix_tail_seed(solver,
                                                 worker.rejected_pivot);
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

static int kls_parallel_refactor_tail_from_block(kls_solver *solver,
                                                 double *numeric_values,
                                                 UF_long start_block) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      numeric_values == NULL ||
      solver->col_ptr == NULL || solver->row_idx == NULL ||
      solver->symbolic->R == NULL || solver->symbolic->Q == NULL ||
      (solver->common.scale > 0 && solver->numeric->Rs == NULL) ||
      (solver->common.scale <= 0 && solver->numeric->Rs != NULL) ||
      solver->numeric->Udiag == NULL || solver->numeric->Offp == NULL ||
      solver->numeric->Offx == NULL || solver->numeric->Lip == NULL ||
      solver->numeric->Llen == NULL || solver->numeric->Uip == NULL ||
      solver->numeric->Ulen == NULL || solver->numeric->LUbx == NULL ||
      solver->numeric->Pinv == NULL ||
      start_block > solver->symbolic->nblocks) {
    return -1;
  }

  const UF_long nblocks = solver->symbolic->nblocks;
  if (start_block == nblocks) {
    return 1;
  }
  const int scaled = solver->common.scale > 0;
  const UF_long remaining_blocks = nblocks - start_block;
  int thread_count = solver->options.threads;
  if ((UF_long)thread_count > remaining_blocks) {
    thread_count = (int)remaining_blocks;
  }
  if (thread_count < 2) {
    return -1;
  }
  if (kls_refactor_pool_map_is_worthwhile(solver) &&
      solver->refactor_col_ptr == NULL) {
    (void)kls_build_refactor_map(solver);
  }

  trilinos_klu_l_common *common = &solver->common;
  common->status = TRILINOS_KLU_OK;
  common->numerical_rank = KLS_KLU_EMPTY;
  common->singular_col = KLS_KLU_EMPTY;
  common->nrealloc = 0;
  if (scaled &&
      !trilinos_klu_l_scale((UF_long)common->scale, solver->n,
                            solver->col_ptr, solver->row_idx,
                            numeric_values, solver->numeric->Rs, NULL,
                            common)) {
    return 0;
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
  if (!run_refactor_pool(solver, numeric_values, thread_count, start_block,
                         1, &invalid, &pivot_rejected,
                         &rejected_pivot, &rejected_pivot_col,
                         &rejected_row, &rejected_multiplier_abs,
                         &rejected_pivot_abs, &rejected_candidate_abs,
                         &singular, &numerical_rank, &singular_col,
                         &pool_prefix_current)) {
    kls_record_fast_factor_failure(
      solver, KLS_FAST_FACTOR_FAIL_POOL_INVALID,
      TRILINOS_KLU_OUT_OF_MEMORY);
    common->status = TRILINOS_KLU_OUT_OF_MEMORY;
    return 0;
  }

  solver->fast_repaired_parallel_tail_blocks += remaining_blocks;
  solver->stats.fast_repaired_parallel_tail_blocks =
    (int64_t)solver->fast_repaired_parallel_tail_blocks;

  if (invalid) {
    kls_record_fast_factor_failure(
      solver, KLS_FAST_FACTOR_FAIL_POOL_INVALID, TRILINOS_KLU_INVALID);
    common->status = TRILINOS_KLU_INVALID;
    return 0;
  }
  if (pivot_rejected) {
    solver->fast_reject_refresh_state =
      pool_prefix_current ? KLS_FAST_REJECT_REFRESH_PREFIX
                          : KLS_FAST_REJECT_REFRESH_UNKNOWN;
    if (pool_prefix_current) {
      kls_record_fast_reject_block_suffix_tail_seed(solver, rejected_pivot);
    } else {
      kls_clear_fast_reject_tail_seed(solver);
    }
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
  if (scaled && !kls_parallel_refactor_permute_scale(solver)) {
    common->status = TRILINOS_KLU_INVALID;
    return 0;
  }
  if (!singular) {
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
  shared.solver = solver;
  shared.native_short_l = kls_mapped_native_short_l_enabled(solver);
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
    kls_record_fast_factor_failure(
      solver, KLS_FAST_FACTOR_FAIL_MAPPED_INVALID, TRILINOS_KLU_INVALID);
    common->status = TRILINOS_KLU_INVALID;
    return 0;
  }
  if (worker.pivot_rejected) {
    solver->fast_reject_refresh_state = KLS_FAST_REJECT_REFRESH_PREFIX;
    kls_record_fast_reject_block_suffix_tail_seed(solver,
                                                 worker.rejected_pivot);
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

/* ====================================================================== */
/* Supernodal dense-panel refactor (task #12 port, refactor path only).   */
/*                                                                        */
/* Prototype: bench/supernodal_bench.c.  Amalgamated supernodes over the  */
/* sorted numeric; column-major dense panels [exact U rows | diag block | */
/* below-row union]; per-edge exact producer-column sublists and          */
/* precomputed scatter destinations; per-edge UB gather + restricted      */
/* unit-lower TRSM + GEMM into a row-major temp + row scatter; fused      */
/* no-pivot getrf.  Values are written back into the packed KLU numeric   */
/* through per-block position lists so solves and every downstream        */
/* consumer see an ordinary refactored numeric.  Serial engine; adopted   */
/* per matrix only after beating the incumbent serial path on a timed     */
/* first call (permanent fallback otherwise).                             */
/* ====================================================================== */

typedef struct kls_snb_edge {
  int32_t producer;
  int32_t u_start;
  int32_t sub_count;
  int64_t sub_off;
  int64_t dst_off;
} kls_snb_edge;

typedef struct kls_snb_snode {
  int64_t start;
  int32_t width;
  int32_t below;
  int64_t panel;
  int32_t height;
  int32_t uw;
  int64_t rows_off;
  int64_t edges_off;
  int32_t edge_count;
} kls_snb_snode;

typedef struct kls_snb_block {
  UF_long k1, nk;
  int64_t snode_count;
  kls_snb_snode *snodes;
  int32_t *col2snode;
  int32_t *below_rows;
  kls_snb_edge *edges;
  int32_t *subcols;
  int32_t *edge_dsts;
  double *panels;
  int64_t panel_doubles;
  int32_t *writeback;   /* per LU entry (block walk order): panel index */
  int64_t *wb_col_off;  /* per local column: offset into writeback[] */
  _Atomic int32_t *published; /* per snode: getrf columns finalized */
} kls_snb_block;

typedef struct kls_snb_entry {
  int64_t dst;
  int64_t src;
  int64_t oldrow;
} kls_snb_entry;

typedef struct kls_snb_state {
  kls_snb_block *blocks;
  UF_long nblocks;
  kls_snb_entry *off;      /* off-diagonal scatter entries, Offx order */
  kls_snb_entry *blk;      /* in-panel scatter entries */
  int64_t *col_off_ptr;    /* per global col: range in off[] */
  int64_t *col_blk_ptr;    /* per global col: range in blk[] */
  double **ub_w;           /* per-worker scratch */
  double **gemm_w;
  int nworkers;
  double *rs_tmp;          /* n scratch for the Rs permute */
  int wnarrow;
  int wmax;
  /* snode-level schedule for the EGraph integration (rock 1): levels of
   * the per-block snode DAGs, merged across independent blocks */
  int64_t slevel_count;
  int64_t *slevel_ptr;     /* level -> range in slevel_tasks */
  int32_t *slevel_block;   /* task -> block id */
  int32_t *slevel_snode;   /* task -> snode id within block */
  int64_t seg_count;       /* pipelined segments between coop levels */
  int64_t *seg_ptr;        /* segment -> range in the task list */
  unsigned char *seg_coop; /* segment is a single cooperative snode */
  _Atomic int64_t singular_rank; /* min singular pivot seen this refactor */
  _Atomic int failed;            /* halt_if_singular tripped */
  double coop_work_share; /* fraction in cooperative fat-chain snodes */
} kls_snb_state;

#define KLS_SNB_WMAX 32
#define KLS_SNB_ZETA 0.2
#define KLS_SNB_MAX_ENTRIES 150000000
/* Work threshold above which a lone task on a level is executed by the
 * whole crew (column-sliced edges, row-sliced GEMM/scatter and getrf). */
#define KLS_SNB_COOP_MIN_DOUBLES 32768
