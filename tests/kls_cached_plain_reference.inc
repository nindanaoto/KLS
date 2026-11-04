static void reference_plain_block(kls_parallel_refactor_worker *worker,
                                        UF_long block) {
  kls_parallel_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  const UF_long *ap = shared->col_ptr;
  const UF_long *ai = shared->row_idx;
  const UF_long *map_col_ptr = shared->map_col_ptr;
  const UF_long *map_row_idx = shared->map_row_idx;
  const UF_long *map_input_pos = shared->map_input_pos;
  const UF_long *map_block_start = shared->map_block_start;
  const int32_t *map_row_idx32 = solver->refactor_row_idx32;
  const int32_t *map_input_pos32 = shared->direct_user_values
    ? solver->refactor_input_user_pos32
    : solver->refactor_input_pos32;
  const int32_t *map_internal_pos32 =
    solver->refactor_input_pos32;
  const double *prepared_scale = solver->prepared_value_scale;
  double *owned_values = solver->values;
  const int fast_direct_scaled_map32 = shared->direct_user_values &&
    prepared_scale != NULL && map_row_idx32 != NULL &&
    map_input_pos32 != NULL && map_internal_pos32 != NULL &&
    owned_values != NULL;
  const int fast_unscaled_map32 = shared->scale <= 0 &&
    !fast_direct_scaled_map32 &&
    map_row_idx32 != NULL && map_input_pos32 != NULL;
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

  /* Short sparse L columns are latency-bound rather than bandwidth-bound.
     Their duplicate i32 scatter pays sign-extension/address-generation
     overhead without moving enough indices to repay it.  Retain the compact
     stream for longer factors and use the packed native rows for this
     representation-defined short-L regime. */
  const int native_short_l = shared->native_short_l;
  const UF_long k1 = r[block];
  const UF_long k2 = r[block + 1u];
  const UF_long nk = k2 - k1;

  if (nk == 1u) {
    const UF_long oldcol = q[k1];
    UF_long poff = offp[k1];
    const UF_long poff_end = offp[k1 + 1u];
    double s = 0.0;
    if (map_col_ptr != NULL && map_block_start != NULL) {
      if (fast_direct_scaled_map32) {
        for (UF_long p = map_col_ptr[k1]; p < map_block_start[k1]; ++p) {
          const UF_long internal = (UF_long)map_internal_pos32[p];
          const double value =
            ax[(UF_long)map_input_pos32[p]] * prepared_scale[internal];
          owned_values[internal] = value;
          offx[poff++] = value;
        }
        for (UF_long p = map_block_start[k1];
             p < map_col_ptr[k1 + 1u]; ++p) {
          const UF_long internal = (UF_long)map_internal_pos32[p];
          const double value =
            ax[(UF_long)map_input_pos32[p]] * prepared_scale[internal];
          owned_values[internal] = value;
          s = value;
        }
      } else if (fast_unscaled_map32) {
        for (UF_long p = map_col_ptr[k1]; p < map_block_start[k1]; ++p) {
          offx[poff++] = ax[(UF_long)map_input_pos32[p]];
        }
        for (UF_long p = map_block_start[k1]; p < map_col_ptr[k1 + 1u]; ++p) {
          s = ax[(UF_long)map_input_pos32[p]];
        }
      } else {
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
      if (fast_direct_scaled_map32) {
        for (UF_long p = map_col_ptr[global_col];
             p < map_block_start[global_col]; ++p) {
          const UF_long internal = (UF_long)map_internal_pos32[p];
          const double value =
            ax[(UF_long)map_input_pos32[p]] * prepared_scale[internal];
          owned_values[internal] = value;
          offx[poff++] = value;
        }
        for (UF_long p = map_block_start[global_col];
             p < map_col_ptr[global_col + 1u]; ++p) {
          const UF_long internal = (UF_long)map_internal_pos32[p];
          const double value =
            ax[(UF_long)map_input_pos32[p]] * prepared_scale[internal];
          owned_values[internal] = value;
          x[(UF_long)map_row_idx32[p] - k1] = value;
        }
      } else if (fast_unscaled_map32) {
        for (UF_long p = map_col_ptr[global_col];
             p < map_block_start[global_col]; ++p) {
          offx[poff++] = ax[(UF_long)map_input_pos32[p]];
        }
        for (UF_long p = map_block_start[global_col];
             p < map_col_ptr[global_col + 1u]; ++p) {
          x[(UF_long)map_row_idx32[p] - k1] =
            ax[(UF_long)map_input_pos32[p]];
        }
      } else {
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
    const UF_long *snode_run_end = shared->snode_run_end;
    UF_long up = 0;
    while (up < ucol_len) {
      const UF_long j = ui[up];
      if (snode_run_end != NULL) {
        const UF_long consumed = kls_snode_batch_consume(
          lu, lip, llen, ui, ux, ucol_len, up, x, k1, snode_run_end);
        if (consumed != 0u) {
          up += consumed;
          continue;
        }
      }
      if (shared->padded_src != NULL && shared->padded_src->padded_active) {
        const UF_long consumed = kls_padded_run_consume(
          shared->padded_src, k1, j, ui, ux, ucol_len, up, x);
        if (consumed != 0u) {
          up += consumed;
          continue;
        }
      }
      const double ujk = x[j];
      x[j] = 0.0;
      ux[up] = ujk;

      if (ujk != 0.0) {
        const UF_long dep_global = k1 + j;
        UF_long *li = NULL;
        double *lx = NULL;
        UF_long lcol_len = 0;
        kls_klu_get_pointer(lu, lip, llen, j, &li, &lx, &lcol_len);
        if (lcol_len <= 3u) {
          if (lcol_len > 0u) {
            x[li[0]] -= lx[0] * ujk;
          }
          if (lcol_len > 1u) {
            x[li[1]] -= lx[1] * ujk;
          }
          if (lcol_len > 2u) {
            x[li[2]] -= lx[2] * ujk;
          }
        } else if (native_short_l || solver->refactor_l_indices32 == NULL) {
          kls_scatter_subtract(x, li, lx, lcol_len, ujk);
        } else {
          kls_scatter_subtract_refactor_l(solver, x, dep_global, li, lx,
                                          lcol_len, ujk);
        }
      }
      up++;
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
    if (lcol_len == 1u) {
      lx[0] = x[li[0]] / ukk;
      x[li[0]] = 0.0;
    } else if (lcol_len > 1u) {
      const double ukk_recip = 1.0 / ukk;
      for (UF_long p = 0; p < lcol_len; ++p) {
        const UF_long i = li[p];
        lx[p] = x[i] * ukk_recip;
        x[i] = 0.0;
      }
    }
    if (shared->padded_src != NULL) {
      const kls_solver *ps = shared->padded_src;
      const UF_long run1 = ps->padded_run_of[global_col];
      if (run1 != 0u) {
        const UF_long run = run1 - 1u;
        const UF_long ulen_run =
          ps->padded_union_ptr[run + 1u] - ps->padded_union_ptr[run];
        double *panel_row = ps->padded_panel_values +
          ps->padded_panel_ptr[run] +
          (global_col - ps->padded_run_start[run]) * ulen_run;
        const UF_long *slots =
          ps->padded_slots + ps->padded_slot_ptr[global_col];
        for (UF_long p = 0; p < lcol_len; ++p) {
          panel_row[slots[p]] = lx[p];
        }
      }
    }
  }
}

