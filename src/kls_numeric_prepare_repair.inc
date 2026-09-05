/* Value preparation, fast rejection, repair, and retained numeric metadata.
 * Included by kls_numeric_engines.inc in the original definition order.
 */

static int kls_build_prepared_value_scale(kls_solver *solver) {
  if (solver == NULL || solver->col_ptr == NULL || solver->row_idx == NULL ||
      (solver->row_scale == NULL && solver->col_scale == NULL)) {
    return 0;
  }
  if (solver->prepared_value_scale != NULL || solver->nnz == 0u) {
    return 1;
  }
  if (solver->nnz > (UF_long)(SIZE_MAX / sizeof(double))) {
    return 0;
  }
  double *scale =
    (double *)malloc((size_t)solver->nnz * sizeof(*scale));
  if (scale == NULL) {
    return 0;
  }
  for (UF_long col = 0u; col < solver->n; ++col) {
    const double cs =
      solver->col_scale != NULL ? solver->col_scale[col] : 1.0;
    for (UF_long p = solver->col_ptr[col];
         p < solver->col_ptr[col + 1u]; ++p) {
      const UF_long row = solver->row_idx[p];
      if (row >= solver->n || p >= solver->nnz) {
        free(scale);
        return 0;
      }
      const double rs =
        solver->row_scale != NULL ? solver->row_scale[row] : 1.0;
      scale[p] = rs * cs;
    }
  }
  solver->prepared_value_scale = scale;
  return 1;
}

static int kls_build_prepared_value_input_pos(kls_solver *solver) {
  if (solver == NULL || solver->input_to_csc == NULL) {
    return 0;
  }
  if (solver->prepared_value_input_pos != NULL || solver->nnz == 0u) {
    return 1;
  }
  if (solver->nnz > (UF_long)(SIZE_MAX / sizeof(UF_long))) {
    return 0;
  }
  UF_long *input_pos =
    (UF_long *)malloc((size_t)solver->nnz * sizeof(*input_pos));
  if (input_pos == NULL) {
    return 0;
  }
  for (UF_long target = 0u; target < solver->nnz; ++target) {
    input_pos[target] = KLS_KLU_EMPTY;
  }
  for (UF_long p = 0u; p < solver->nnz; ++p) {
    const UF_long target = solver->input_to_csc[p];
    if (target >= solver->nnz || input_pos[target] != KLS_KLU_EMPTY) {
      free(input_pos);
      return 0;
    }
    input_pos[target] = p;
  }
  solver->prepared_value_input_pos = input_pos;
  return 1;
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
  const int generic_fused_value_prep =
    kls_repeated_update_workload(&solver->options) &&
    getenv("KLS_DISABLE_GENERIC_FUSED_VALUE_PREP") == NULL;
  if ((getenv("KLS_ENABLE_FUSED_VALUE_PREP") != NULL ||
       generic_fused_value_prep) &&
      (solver->row_scale != NULL || solver->col_scale != NULL) &&
      kls_build_prepared_value_scale(solver)) {
    const double *restrict scale = solver->prepared_value_scale;
    double *restrict prepared = solver->values;
    if (solver->input_to_csc == NULL) {
      for (UF_long p = 0u; p < solver->nnz; ++p) {
        prepared[p] = values[p] * scale[p];
      }
    } else if ((getenv("KLS_ENABLE_GATHER_FUSED_VALUE_PREP") != NULL ||
                generic_fused_value_prep) &&
               kls_build_prepared_value_input_pos(solver)) {
      const UF_long *restrict input_pos =
        solver->prepared_value_input_pos;
      for (UF_long target = 0u; target < solver->nnz; ++target) {
        prepared[target] = values[input_pos[target]] * scale[target];
      }
    } else {
      for (UF_long p = 0u; p < solver->nnz; ++p) {
        const UF_long target = solver->input_to_csc[p];
        if (target >= solver->nnz) {
          return KLS_ERR_INVALID_ARGUMENT;
        }
        prepared[target] = values[p] * scale[target];
      }
    }
    *values_out = prepared;
    return KLS_OK;
  }
  if (kls_parallel_prepare_permuted_values(solver, values)) {
    *values_out = solver->values;
    return KLS_OK;
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

/* A refactorization is a pure function of the fixed pattern and its numeric
   input.  Several real workloads (and the paper harnesses) submit occasional
   or long runs of byte-identical values.  Reusing the current factor is exact
   in that case and avoids touching the substantially larger filled factors.

   The first changed refactor after a factor disables the cache, so workloads
   that update values every step pay only one comparison and no per-step
   snapshot copy.  Once an unchanged step has been observed, keep tracking the
   input so mixed changed/unchanged workloads can benefit as well. */
static void kls_arm_unchanged_refactor_cache(kls_solver *solver,
                                             const double *values) {
  if (solver == NULL || values == NULL) {
    return;
  }
  solver->refactor_input_snapshot_valid = 0;
  solver->unchanged_refactor_state = 0;
  if (getenv("KLS_DISABLE_UNCHANGED_REFACTOR") != NULL ||
      solver->nnz > (UF_long)(SIZE_MAX / sizeof(double))) {
    free(solver->refactor_input_snapshot);
    solver->refactor_input_snapshot = NULL;
    solver->unchanged_refactor_state = -1;
    return;
  }
  if (solver->nnz > 0u && solver->refactor_input_snapshot == NULL) {
    solver->refactor_input_snapshot =
      (double *)malloc((size_t)solver->nnz * sizeof(double));
  }
  if (solver->nnz > 0u && solver->refactor_input_snapshot == NULL) {
    solver->unchanged_refactor_state = -1;
    return;
  }
  if (solver->nnz > 0u) {
    memcpy(solver->refactor_input_snapshot, values,
           (size_t)solver->nnz * sizeof(double));
  }
  solver->refactor_input_snapshot_valid = 1;
  solver->unchanged_refactor_state = 1;
}

static int kls_refactor_input_is_unchanged(kls_solver *solver,
                                           const double *values) {
  if (solver == NULL || values == NULL ||
      solver->unchanged_refactor_state <= 0 ||
      !solver->refactor_input_snapshot_valid) {
    return 0;
  }
  const int unchanged = solver->nnz == 0u ||
    memcmp(solver->refactor_input_snapshot, values,
           (size_t)solver->nnz * sizeof(double)) == 0;
  if (unchanged) {
    solver->unchanged_refactor_state = 2;
    return 1;
  }

  solver->refactor_input_snapshot_valid = 0;
  if (solver->unchanged_refactor_state == 1) {
    /* Continuously changing workload: do not add a snapshot memcpy to every
       otherwise-cheap refactorization.  A later explicit factor re-arms it. */
    free(solver->refactor_input_snapshot);
    solver->refactor_input_snapshot = NULL;
    solver->unchanged_refactor_state = -1;
  }
  return 0;
}

static void kls_refresh_unchanged_refactor_cache(kls_solver *solver,
                                                 const double *values) {
  if (solver == NULL || values == NULL ||
      solver->unchanged_refactor_state != 2 ||
      (solver->nnz > 0u && solver->refactor_input_snapshot == NULL)) {
    return;
  }
  if (solver->nnz > 0u) {
    memcpy(solver->refactor_input_snapshot, values,
           (size_t)solver->nnz * sizeof(double));
  }
  solver->refactor_input_snapshot_valid = 1;
}

/* KLU scaling is a separate O(nnz) metadata pass: it only derives Rs from
   Ax, then the numeric refactor reads Rs while still rebuilding every L/U
   value.  Long runs of an identical scaled matrix therefore need not repeat
   that scan merely because the full-factor unchanged shortcut is disabled.

   This cache is intentionally one-way for each numeric.  The first mismatch
   frees the snapshot and permanently restores the ordinary scale pass, so a
   transient or continuously changing workload pays one comparison and never
   adds a per-refactor snapshot copy. */
static void kls_reset_lean_scale_input_cache(kls_solver *solver) {
  if (solver == NULL) {
    return;
  }
  free(solver->lean_scale_input_snapshot);
  free(solver->lean_scale_rs_snapshot);
  solver->lean_scale_input_snapshot = NULL;
  solver->lean_scale_rs_snapshot = NULL;
  solver->lean_scale_input_state = 0;
}

static int kls_lean_scale_input_is_unchanged(kls_solver *solver,
                                             const double *values) {
  if (solver == NULL || values == NULL ||
      solver->lean_scale_input_state != 1 ||
      (solver->nnz > 0u && solver->lean_scale_input_snapshot == NULL) ||
      solver->numeric == NULL || solver->numeric->Rs == NULL ||
      (solver->n > 0u && solver->lean_scale_rs_snapshot == NULL)) {
    return 0;
  }
  if (solver->nnz == 0u ||
      memcmp(solver->lean_scale_input_snapshot, values,
             (size_t)solver->nnz * sizeof(double)) == 0) {
    if (solver->n > 0u) {
      /* The refactor epilogue permutes Numeric->Rs into pivot order for
         solve.  Restore the cached input-row order expected by the next
         numeric walk before declaring this scale pass reusable. */
      memcpy(solver->numeric->Rs, solver->lean_scale_rs_snapshot,
             (size_t)solver->n * sizeof(double));
    }
    return 1;
  }
  free(solver->lean_scale_input_snapshot);
  free(solver->lean_scale_rs_snapshot);
  solver->lean_scale_input_snapshot = NULL;
  solver->lean_scale_rs_snapshot = NULL;
  solver->lean_scale_input_state = -1;
  return 0;
}

static void kls_arm_lean_scale_input_cache(kls_solver *solver,
                                           const double *values) {
  if (solver == NULL || values == NULL ||
      solver->lean_scale_input_state != 0 || solver->numeric == NULL ||
      solver->numeric->Rs == NULL) {
    return;
  }
  if (getenv("KLS_DISABLE_LEAN_SCALE_CACHE") != NULL ||
      solver->nnz > (UF_long)(SIZE_MAX / sizeof(double)) ||
      solver->n > (UF_long)(SIZE_MAX / sizeof(double))) {
    solver->lean_scale_input_state = -1;
    return;
  }
  if (solver->nnz > 0u) {
    solver->lean_scale_input_snapshot =
      (double *)malloc((size_t)solver->nnz * sizeof(double));
    if (solver->lean_scale_input_snapshot == NULL) {
      solver->lean_scale_input_state = -1;
      return;
    }
  }
  if (solver->n > 0u) {
    solver->lean_scale_rs_snapshot =
      (double *)malloc((size_t)solver->n * sizeof(double));
    if (solver->lean_scale_rs_snapshot == NULL) {
      free(solver->lean_scale_input_snapshot);
      solver->lean_scale_input_snapshot = NULL;
      solver->lean_scale_input_state = -1;
      return;
    }
  }
  if (solver->nnz > 0u) {
    memcpy(solver->lean_scale_input_snapshot, values,
           (size_t)solver->nnz * sizeof(double));
  }
  if (solver->n > 0u) {
    memcpy(solver->lean_scale_rs_snapshot, solver->numeric->Rs,
           (size_t)solver->n * sizeof(double));
  }
  solver->lean_scale_input_state = 1;
}

/* The retained row pattern is already a validated permutation of every CSC
   entry.  Derive KLU's input-row scales from that stream with one contiguous
   output per row instead of scattering nnz max/sum updates through Rs while
   walking columns.  Pnum maps each numeric row back to the input row whose
   scale KLU's solve contract expects. */
static int kls_scale_from_row_refactor_stream(kls_solver *solver,
                                              const double *values) {
  if (solver == NULL || values == NULL || solver->numeric == NULL ||
      solver->numeric->Rs == NULL || solver->numeric->Pnum == NULL ||
      solver->row_refactor_input_ptr == NULL ||
      solver->row_refactor_input_pos == NULL ||
      solver->row_refactor_input_ptr[solver->n] != solver->nnz ||
      solver->common.scale <= 0) {
    return 0;
  }
  const UF_long *restrict ptr = solver->row_refactor_input_ptr;
  const UF_long *restrict pos = solver->row_refactor_input_pos;
  const UF_long *restrict pnum = solver->numeric->Pnum;
  double *restrict rs = solver->numeric->Rs;
  const int sum_scale = solver->common.scale == 1;
  for (UF_long row = 0u; row < solver->n; ++row) {
    double scale = 0.0;
    for (UF_long p = ptr[row]; p < ptr[row + 1u]; ++p) {
      const double magnitude = fabs(values[pos[p]]);
      if (sum_scale) {
        scale += magnitude;
      } else if (magnitude > scale) {
        scale = magnitude;
      }
    }
    rs[pnum[row]] = scale != 0.0 ? scale : 1.0;
  }
  return 1;
}

static int kls_parallel_refactor_is_eligible(const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL) {
    return 0;
  }
  {
    /* The block pool can only expose concurrency between BTF components.
       Base admission on the amount of real work and the number of runnable
       components, not on a matrix-family coverage window.  Near-single-block
       factors are left to the dependency-graph engine, which can parallelize
       inside their dominant block. */
    const UF_long min_blocks =
      (UF_long)(2 * solver->options.threads) > 8u
        ? (UF_long)(2 * solver->options.threads) : 8u;
    const UF_long independent_rows =
      solver->n - solver->symbolic->maxblock;
    if (solver->options.threads <= 1 || solver->n < 5000u ||
        solver->symbolic->nblocks < min_blocks ||
        solver->symbolic->maxblock == solver->n ||
        independent_rows < (UF_long)(64 * solver->options.threads) ||
        solver->common.flops < 2.0e7) {
      return 0;
    }
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
  if (!kls_parallel_refactor_is_eligible(solver)) {
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
      solver->numeric->LUsize == NULL ||
      solver->numeric->LUbx[block] == NULL) {
    return 0;
  }
  double *lu = (double *)solver->numeric->LUbx[block];
  const size_t lusize = solver->numeric->LUsize[block];
  const UF_long *uip = solver->numeric->Uip + k1;
  const UF_long *ulen = solver->numeric->Ulen + k1;
  if (uip[local_col] == KLS_KLU_EMPTY ||
      (size_t)uip[local_col] > lusize) {
    return 0;
  }
  const size_t index_units = kls_klu_units_for_indices(ulen[local_col]);
  if (index_units > lusize - (size_t)uip[local_col] ||
      ulen[local_col] >
        (UF_long)(lusize - (size_t)uip[local_col] - index_units)) {
    return 0;
  }
  UF_long *ui = NULL;
  double *ux = NULL;
  UF_long ucol_len = 0;
  kls_klu_get_pointer(lu, uip, ulen, local_col, &ui, &ux, &ucol_len);
  (void)ux;

  double work = 1.0;
  for (UF_long p = 0; p < ucol_len; ++p) {
    const UF_long dep = ui[p];
    if (dep < 0 || dep >= local_col || dep >= nk) {
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

static int kls_fast_reject_range_work(const kls_solver *solver,
                                      UF_long block,
                                      UF_long k1,
                                      UF_long nk,
                                      UF_long local_begin,
                                      UF_long local_end,
                                      UF_long *columns_out,
                                      double *work_out) {
  if (solver == NULL || columns_out == NULL || work_out == NULL ||
      local_begin >= local_end || local_end > nk) {
    return 0;
  }
  double work = 0.0;
  for (UF_long k = local_begin; k < local_end; ++k) {
    double column_work = 0.0;
    if (!kls_fast_reject_column_work(solver, block, k1, nk, k,
                                     &column_work)) {
      return 0;
    }
    work += column_work;
  }
  *columns_out = local_end - local_begin;
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
  return kls_fast_reject_range_work(solver, block, k1, nk, local_reject,
                                    nk, columns_out, work_out);
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
  KLS_SET_OPTIONAL_OUTPUT(parent_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(owned_parent_out, NULL);
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

static void kls_update_factor_etree_stats(kls_solver *solver) {
  if (solver == NULL || solver->factor_etree_stats_valid) {
    return;
  }
  kls_invalidate_factor_etree_stats(solver);
  if (solver->symbolic == NULL || solver->symbolic->R == NULL ||
      solver->symbolic->nblocks == 0u || solver->n == 0u) {
    solver->factor_etree_stats_valid = 1;
    return;
  }

  UF_long best_block = KLS_KLU_EMPTY;
  UF_long best_k1 = 0;
  UF_long best_k2 = 0;
  UF_long best_nk = 0;
  for (UF_long block = 0; block < solver->symbolic->nblocks; ++block) {
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    if (k1 > k2 || k2 > solver->n) {
      solver->factor_etree_stats_valid = 1;
      return;
    }
    if (k2 - k1 > best_nk) {
      best_block = block;
      best_k1 = k1;
      best_k2 = k2;
      best_nk = k2 - k1;
    }
  }

  if (best_block == KLS_KLU_EMPTY || best_nk == 0u) {
    solver->factor_etree_stats_valid = 1;
    return;
  }
  UF_long *parent = (UF_long *)malloc((size_t)best_nk * sizeof(*parent));
  UF_long *level = (UF_long *)calloc((size_t)best_nk, sizeof(*level));
  UF_long *width = (UF_long *)calloc((size_t)best_nk + 1u, sizeof(*width));
  UF_long *fanout = (UF_long *)calloc((size_t)best_nk, sizeof(*fanout));
  if (parent == NULL || level == NULL || width == NULL || fanout == NULL) {
    free(parent);
    free(level);
    free(width);
    free(fanout);
    solver->factor_etree_stats_valid = 1;
    return;
  }
  if (!kls_build_ordered_block_etree(solver, best_k1, best_k2, parent)) {
    free(parent);
    free(level);
    free(width);
    free(fanout);
    solver->factor_etree_stats_valid = 1;
    return;
  }

  UF_long edges = 0;
  UF_long roots = 0;
  UF_long max_fanout = 0;
  for (UF_long col = 0; col < best_nk; ++col) {
    level[col] = 1u;
    if (parent[col] == KLS_KLU_EMPTY) {
      roots++;
      continue;
    }
    if (parent[col] >= best_nk || parent[col] <= col) {
      free(parent);
      free(level);
      free(width);
      free(fanout);
      solver->factor_etree_stats_valid = 1;
      return;
    }
    edges++;
    fanout[parent[col]]++;
  }
  UF_long leaves = 0;
  for (UF_long col = 0; col < best_nk; ++col) {
    if (fanout[col] == 0u) {
      leaves++;
    }
    if (fanout[col] > max_fanout) {
      max_fanout = fanout[col];
    }
    const UF_long next = parent[col];
    if (next != KLS_KLU_EMPTY && level[next] < level[col] + 1u) {
      level[next] = level[col] + 1u;
    }
  }

  UF_long level_count = 0;
  UF_long max_width = 0;
  for (UF_long col = 0; col < best_nk; ++col) {
    if (level[col] > level_count) {
      level_count = level[col];
    }
    if (level[col] <= best_nk) {
      const UF_long count = ++width[level[col]];
      if (count > max_width) {
        max_width = count;
      }
    }
  }

  solver->stats.factor_etree_block_start = (int64_t)best_k1;
  solver->stats.factor_etree_block_size = (int64_t)best_nk;
  solver->stats.factor_etree_levels = (int64_t)level_count;
  solver->stats.factor_etree_max_width = (int64_t)max_width;
  solver->stats.factor_etree_edges = (int64_t)edges;
  solver->stats.factor_etree_root_columns = (int64_t)roots;
  solver->stats.factor_etree_leaf_columns = (int64_t)leaves;
  solver->stats.factor_etree_max_fanout = (int64_t)max_fanout;
  solver->factor_etree_stats_valid = 1;
  free(parent);
  free(level);
  free(width);
  free(fanout);
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
      solver->fast_reject_tail_seed_cols != NULL &&
      solver->fast_reject_tail_parent != NULL &&
      solver->fast_reject_tail_child_count != NULL &&
      solver->fast_reject_tail_level != NULL &&
      solver->fast_reject_tail_marks != NULL) {
    return 1;
  }
  if (capacity > (UF_long)(SIZE_MAX / sizeof(*solver->fast_reject_tail_cols)) ||
      capacity >
        (UF_long)(SIZE_MAX / sizeof(*solver->fast_reject_tail_seed_cols)) ||
      capacity >
        (UF_long)(SIZE_MAX / sizeof(*solver->fast_reject_tail_parent)) ||
      capacity >
        (UF_long)(SIZE_MAX / sizeof(*solver->fast_reject_tail_child_count)) ||
      capacity >
        (UF_long)(SIZE_MAX / sizeof(*solver->fast_reject_tail_level)) ||
      capacity > (UF_long)(SIZE_MAX / sizeof(*solver->fast_reject_tail_marks))) {
    return 0;
  }

  UF_long *cols =
    (UF_long *)malloc((size_t)capacity *
                      sizeof(*solver->fast_reject_tail_cols));
  UF_long *seed_cols =
    (UF_long *)malloc((size_t)capacity *
                      sizeof(*solver->fast_reject_tail_seed_cols));
  UF_long *tail_parent =
    (UF_long *)malloc((size_t)capacity *
                      sizeof(*solver->fast_reject_tail_parent));
  UF_long *child_count =
    (UF_long *)malloc((size_t)capacity *
                      sizeof(*solver->fast_reject_tail_child_count));
  UF_long *level =
    (UF_long *)malloc((size_t)capacity *
                      sizeof(*solver->fast_reject_tail_level));
  unsigned int *marks =
    (unsigned int *)calloc((size_t)capacity,
                           sizeof(*solver->fast_reject_tail_marks));
  if (cols == NULL || seed_cols == NULL || tail_parent == NULL ||
      child_count == NULL || level == NULL || marks == NULL) {
    free(cols);
    free(seed_cols);
    free(tail_parent);
    free(child_count);
    free(level);
    free(marks);
    return 0;
  }
  /* the tail-repair consumers index tail_parent by rows the per-column
     scan may never have visited: fill the whole array with EMPTY (-1,
     the same value the per-column init writes) so unvisited slots are
     inert instead of heap garbage (measured: the 3x3 smoke transpose
     solve flipped with heap layout through a garbage in-bounds parent) */
  memset(tail_parent, 0xff,
         (size_t)capacity * sizeof(*solver->fast_reject_tail_parent));
  memset(cols, 0xff,
         (size_t)capacity * sizeof(*solver->fast_reject_tail_cols));
  memset(seed_cols, 0xff,
         (size_t)capacity * sizeof(*solver->fast_reject_tail_seed_cols));
  memset(child_count, 0,
         (size_t)capacity * sizeof(*solver->fast_reject_tail_child_count));
  memset(level, 0xff,
         (size_t)capacity * sizeof(*solver->fast_reject_tail_level));
  free(solver->fast_reject_tail_cols);
  free(solver->fast_reject_tail_seed_cols);
  free(solver->fast_reject_tail_parent);
  free(solver->fast_reject_tail_child_count);
  free(solver->fast_reject_tail_level);
  free(solver->fast_reject_tail_marks);
  solver->fast_reject_tail_cols = cols;
  solver->fast_reject_tail_seed_cols = seed_cols;
  solver->fast_reject_tail_parent = tail_parent;
  solver->fast_reject_tail_child_count = child_count;
  solver->fast_reject_tail_level = level;
  solver->fast_reject_tail_marks = marks;
  solver->fast_reject_tail_capacity = capacity;
  solver->fast_reject_tail_mark = 0u;
  solver->fast_reject_tail_plan_mark = 0u;
  solver->fast_reject_tail_plan_block = KLS_KLU_EMPTY;
  solver->fast_reject_tail_plan_k1 = KLS_KLU_EMPTY;
  solver->fast_reject_tail_plan_nk = 0;
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

static void kls_clear_fast_reject_tail_seed(kls_solver *solver) {
  if (solver == NULL) {
    return;
  }
  solver->fast_reject_tail_seed_block = KLS_KLU_EMPTY;
  solver->fast_reject_tail_seed_count = 0;
  solver->fast_reject_tail_seed_valid = 0;
  solver->fast_reject_tail_seed_block_suffix = 0;
  solver->stats.fast_rejected_pivoting_tail_seed_columns = 0;
  solver->stats.fast_rejected_pivoting_tail_row_seed_columns = 0;
  solver->stats.fast_rejected_pivoting_tail_block_seed_columns = 0;
}

static void kls_record_fast_reject_unfinished_tail_seed(
  kls_solver *solver,
  const kls_egraph_refactor_shared *shared,
  UF_long rejected_pivot) {
  kls_clear_fast_reject_tail_seed(solver);
  if (solver == NULL || shared == NULL || shared->pipeline_done == NULL ||
      shared->pipeline_generation == 0u ||
      solver->symbolic == NULL || solver->symbolic->R == NULL ||
      rejected_pivot == KLS_KLU_EMPTY) {
    return;
  }
  const UF_long block = kls_block_for_pivot(solver, rejected_pivot);
  if (block == KLS_KLU_EMPTY || block >= solver->symbolic->nblocks) {
    return;
  }
  const UF_long k1 = solver->symbolic->R[block];
  const UF_long k2 = solver->symbolic->R[block + 1u];
  if (rejected_pivot < k1 || rejected_pivot >= k2) {
    return;
  }
  const UF_long nk = k2 - k1;
  if (!kls_ensure_fast_reject_tail_plan_storage(solver, nk)) {
    return;
  }

  const unsigned int generation = shared->pipeline_generation;
  UF_long seed_count = 0;
  int contains_reject = 0;
  for (UF_long col = k1; col < k2; ++col) {
    if (atomic_load_explicit(&shared->pipeline_done[col],
                             memory_order_acquire) == generation) {
      continue;
    }
    const UF_long local_col = col - k1;
    solver->fast_reject_tail_seed_cols[seed_count++] = local_col;
    if (col == rejected_pivot) {
      contains_reject = 1;
    }
  }
  if (!contains_reject) {
    solver->fast_reject_tail_seed_cols[seed_count++] = rejected_pivot - k1;
  }

  solver->fast_reject_tail_seed_block = block;
  solver->fast_reject_tail_seed_count = seed_count;
  solver->fast_reject_tail_seed_valid = 1;
  solver->stats.fast_rejected_pivoting_tail_seed_columns =
    (int64_t)seed_count;
}

static void kls_record_fast_reject_block_suffix_tail_seed(
  kls_solver *solver,
  UF_long rejected_pivot) {
  kls_clear_fast_reject_tail_seed(solver);
  if (solver == NULL || solver->symbolic == NULL ||
      solver->symbolic->R == NULL || rejected_pivot == KLS_KLU_EMPTY) {
    return;
  }
  const UF_long block = kls_block_for_pivot(solver, rejected_pivot);
  if (block == KLS_KLU_EMPTY || block >= solver->symbolic->nblocks) {
    return;
  }
  const UF_long k1 = solver->symbolic->R[block];
  const UF_long k2 = solver->symbolic->R[block + 1u];
  if (rejected_pivot < k1 || rejected_pivot >= k2) {
    return;
  }
  const UF_long nk = k2 - k1;
  if (!kls_ensure_fast_reject_tail_plan_storage(solver, nk)) {
    return;
  }

  UF_long seed_count = 0;
  for (UF_long col = rejected_pivot; col < k2; ++col) {
    solver->fast_reject_tail_seed_cols[seed_count++] = col - k1;
  }
  solver->fast_reject_tail_seed_block = block;
  solver->fast_reject_tail_seed_count = seed_count;
  solver->fast_reject_tail_seed_valid = seed_count > 0u ? 1 : 0;
  solver->fast_reject_tail_seed_block_suffix =
    solver->fast_reject_tail_seed_valid;
  solver->stats.fast_rejected_pivoting_tail_seed_columns =
    (int64_t)seed_count;
  solver->stats.fast_rejected_pivoting_tail_block_seed_columns =
    (int64_t)seed_count;
}

static int kls_mark_fast_reject_row_tail_seed(kls_solver *solver,
                                              UF_long k1,
                                              UF_long nk,
                                              UF_long local_reject,
                                              unsigned int mark,
                                              UF_long *queued_out) {
  if (solver == NULL || queued_out == NULL || local_reject >= nk ||
      solver->fast_reject_refresh_state !=
        KLS_FAST_REJECT_REFRESH_PREFIX ||
      solver->stats.fast_rejected_row < 0 ||
      solver->row_refactor_tail_rows == NULL ||
      solver->row_refactor_tail_count == 0u) {
    return 0;
  }

  const UF_long reject_global = k1 + local_reject;
  const UF_long block_end = k1 + nk;
  int contains_reject = 0;
  for (UF_long p = 0; p < solver->row_refactor_tail_count; ++p) {
    const UF_long row = solver->row_refactor_tail_rows[p];
    if (row < reject_global || row >= block_end) {
      return 0;
    }
    if (row == reject_global) {
      contains_reject = 1;
    }
  }
  if (!contains_reject) {
    return 0;
  }

  const UF_long queued_before = *queued_out;
  for (UF_long p = 0; p < solver->row_refactor_tail_count; ++p) {
    const UF_long local_col = solver->row_refactor_tail_rows[p] - k1;
    if (!kls_mark_fast_reject_tail_col(solver, nk, local_col, mark,
                                       queued_out)) {
      return 0;
    }
  }
  if (*queued_out == queued_before) {
    return 0;
  }
  const UF_long seed_columns = *queued_out - queued_before;
  solver->stats.fast_rejected_pivoting_tail_seed_columns =
    (int64_t)seed_columns;
  solver->stats.fast_rejected_pivoting_tail_row_seed_columns =
    (int64_t)seed_columns;
  return 1;
}

static void kls_record_fast_reject_pivoting_tail_plan(
  kls_solver *solver,
  UF_long block,
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
  solver->stats.fast_rejected_pivoting_tail_contiguous = 0;
  solver->stats.fast_rejected_pivoting_tail_suffix_exact = 0;
  solver->stats.fast_rejected_pivoting_tail_gap_columns = 0;
  solver->stats.fast_rejected_pivoting_tail_suffix_overcompute_columns = 0;
  solver->stats.fast_rejected_pivoting_tail_suffix_overcompute_work = 0.0;
  solver->stats.fast_rejected_pivoting_tail_etree_edges = 0;
  solver->stats.fast_rejected_pivoting_tail_etree_roots = 0;
  solver->stats.fast_rejected_pivoting_tail_etree_leaves = 0;
  solver->stats.fast_rejected_pivoting_tail_etree_max_fanout = 0;
  solver->stats.fast_rejected_pivoting_tail_etree_levels = 0;
  solver->stats.fast_rejected_pivoting_tail_etree_max_width = 0;
  solver->fast_reject_tail_plan_mark = 0u;
  solver->fast_reject_tail_plan_block = KLS_KLU_EMPTY;
  solver->fast_reject_tail_plan_k1 = KLS_KLU_EMPTY;
  solver->fast_reject_tail_plan_nk = 0;
  if (columns == 0u || solver->fast_reject_tail_cols == NULL ||
      solver->fast_reject_tail_parent == NULL ||
      solver->fast_reject_tail_child_count == NULL ||
      solver->fast_reject_tail_level == NULL ||
      solver->fast_reject_tail_marks == NULL || parent == NULL ||
      local_reject >= nk) {
    return;
  }
  for (UF_long local_col = 0; local_col < nk; ++local_col) {
    solver->fast_reject_tail_parent[local_col] = KLS_KLU_EMPTY;
    solver->fast_reject_tail_child_count[local_col] = 0;
    solver->fast_reject_tail_level[local_col] = 0;
  }

  UF_long previous = KLS_KLU_EMPTY;
  UF_long first_global = KLS_KLU_EMPTY;
  UF_long gap_columns = 0;
  int sorted = 1;
  int contiguous = 1;
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
      first_global = global_col;
      solver->stats.fast_rejected_pivoting_tail_first =
        (int64_t)global_col;
    } else if (previous >= global_col) {
      sorted = 0;
      topological = 0;
      contiguous = 0;
    } else if (global_col != previous + 1u) {
      contiguous = 0;
      gap_columns += global_col - previous - 1u;
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
  solver->stats.fast_rejected_pivoting_tail_contiguous =
    sorted && contiguous && contains_reject;
  solver->stats.fast_rejected_pivoting_tail_gap_columns =
    sorted ? (int64_t)gap_columns : 0;
  const UF_long suffix_columns = nk - local_reject;
  solver->stats.fast_rejected_pivoting_tail_suffix_exact =
    sorted && contiguous && contains_reject &&
    first_global == k1 + local_reject &&
    previous == k1 + nk - 1u &&
    columns == suffix_columns;
  UF_long measured_suffix_columns = 0;
  double suffix_work = 0.0;
  if (kls_fast_reject_suffix_work(solver, block, k1, nk, local_reject,
                                  &measured_suffix_columns, &suffix_work)) {
    solver->stats.fast_rejected_pivoting_tail_suffix_overcompute_columns =
      measured_suffix_columns > columns
        ? (int64_t)(measured_suffix_columns - columns)
        : 0;
    solver->stats.fast_rejected_pivoting_tail_suffix_overcompute_work =
      suffix_work > work ? suffix_work - work : 0.0;
  }
  if (solver->stats.fast_rejected_pivoting_tail_topological) {
    UF_long edge_count = 0;
    UF_long root_count = 0;
    for (UF_long local_col = 0; local_col < nk; ++local_col) {
      if (solver->fast_reject_tail_marks[local_col] != mark) {
        continue;
      }
      const UF_long next = parent[local_col];
      if (next != KLS_KLU_EMPTY && next != local_col && next < nk &&
          solver->fast_reject_tail_marks[next] == mark) {
        solver->fast_reject_tail_parent[local_col] = next;
        edge_count++;
        solver->fast_reject_tail_child_count[next]++;
      } else {
        root_count++;
      }
    }

    UF_long leaf_count = 0;
    UF_long max_fanout = 0;
    UF_long max_level = 0;
    for (UF_long local_col = 0; local_col < nk; ++local_col) {
      if (solver->fast_reject_tail_marks[local_col] != mark) {
        continue;
      }
      if (solver->fast_reject_tail_child_count[local_col] == 0u) {
        leaf_count++;
      }
      if (solver->fast_reject_tail_child_count[local_col] > max_fanout) {
        max_fanout = solver->fast_reject_tail_child_count[local_col];
      }
      const UF_long next = solver->fast_reject_tail_parent[local_col];
      if (next != KLS_KLU_EMPTY && next < nk &&
          solver->fast_reject_tail_marks[next] == mark &&
          solver->fast_reject_tail_level[next] <=
            solver->fast_reject_tail_level[local_col]) {
        solver->fast_reject_tail_level[next] =
          solver->fast_reject_tail_level[local_col] + 1u;
      }
      if (solver->fast_reject_tail_level[local_col] > max_level) {
        max_level = solver->fast_reject_tail_level[local_col];
      }
    }
    const UF_long level_count = columns > 0u ? max_level + 1u : 0u;
    UF_long max_width = 0;
    UF_long *level_width = level_count > 0u
      ? (UF_long *)calloc((size_t)level_count, sizeof(*level_width)) : NULL;
    if (level_count == 0u || level_width != NULL) {
      for (UF_long local_col = 0; local_col < nk; ++local_col) {
        if (solver->fast_reject_tail_marks[local_col] != mark) {
          continue;
        }
        const UF_long level = solver->fast_reject_tail_level[local_col];
        if (level >= level_count) {
          max_width = 0;
          break;
        }
        level_width[level]++;
        if (level_width[level] > max_width) {
          max_width = level_width[level];
        }
      }
      solver->stats.fast_rejected_pivoting_tail_etree_edges =
        (int64_t)edge_count;
      solver->stats.fast_rejected_pivoting_tail_etree_roots =
        (int64_t)root_count;
      solver->stats.fast_rejected_pivoting_tail_etree_leaves =
        (int64_t)leaf_count;
      solver->stats.fast_rejected_pivoting_tail_etree_max_fanout =
        (int64_t)max_fanout;
      solver->stats.fast_rejected_pivoting_tail_etree_levels =
        (int64_t)level_count;
      solver->stats.fast_rejected_pivoting_tail_etree_max_width =
        (int64_t)max_width;
      solver->fast_reject_tail_plan_mark = mark;
      solver->fast_reject_tail_plan_block = block;
      solver->fast_reject_tail_plan_k1 = k1;
      solver->fast_reject_tail_plan_nk = nk;
    }
    free(level_width);
  }
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
  KLS_SET_OPTIONAL_OUTPUT(tail_cols_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(tail_count_out, 0);
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

  int use_unfinished_seed =
    solver->fast_reject_refresh_state == KLS_FAST_REJECT_REFRESH_PREFIX &&
    solver->fast_reject_tail_seed_valid &&
    solver->fast_reject_tail_seed_block == block &&
    solver->fast_reject_tail_seed_count > 0u &&
    solver->fast_reject_tail_seed_cols != NULL;
  const int64_t available_block_seed_columns =
    use_unfinished_seed && solver->fast_reject_tail_seed_block_suffix
      ? (int64_t)solver->fast_reject_tail_seed_count : 0;
  const int use_block_suffix_seed =
    use_unfinished_seed && solver->fast_reject_tail_seed_block_suffix;
  enum {
    KLS_FAST_REJECT_TAIL_SEED_SUFFIX = 0,
    KLS_FAST_REJECT_TAIL_SEED_UNFINISHED = 1,
    KLS_FAST_REJECT_TAIL_SEED_ROW = 2
  };
  int seed_mode =
    solver->fast_reject_refresh_state == KLS_FAST_REJECT_REFRESH_PREFIX
      ? (use_block_suffix_seed
           ? KLS_FAST_REJECT_TAIL_SEED_ROW
           : (use_unfinished_seed
                ? KLS_FAST_REJECT_TAIL_SEED_UNFINISHED
                : KLS_FAST_REJECT_TAIL_SEED_ROW))
      : KLS_FAST_REJECT_TAIL_SEED_SUFFIX;
  solver->stats.fast_rejected_pivoting_tail_row_seed_columns = 0;
  solver->stats.fast_rejected_pivoting_tail_block_seed_columns =
    available_block_seed_columns;
  for (;;) {
    if (solver->fast_reject_tail_mark == UINT_MAX) {
      memset(solver->fast_reject_tail_marks, 0,
             (size_t)solver->fast_reject_tail_capacity *
               sizeof(*solver->fast_reject_tail_marks));
      solver->fast_reject_tail_mark = 0u;
    }
    const unsigned int mark = ++solver->fast_reject_tail_mark;
    UF_long queued = 0;
    const int used_seed_mode = seed_mode;
    if (solver->fast_reject_refresh_state ==
        KLS_FAST_REJECT_REFRESH_PREFIX) {
      if (seed_mode == KLS_FAST_REJECT_TAIL_SEED_ROW) {
        if (!kls_mark_fast_reject_row_tail_seed(
              solver, k1, nk, local_reject, mark, &queued)) {
          solver->stats.fast_rejected_pivoting_tail_row_seed_columns = 0;
          seed_mode = use_unfinished_seed
            ? KLS_FAST_REJECT_TAIL_SEED_UNFINISHED
            : KLS_FAST_REJECT_TAIL_SEED_SUFFIX;
          continue;
        }
      } else if (seed_mode == KLS_FAST_REJECT_TAIL_SEED_UNFINISHED &&
                 use_unfinished_seed) {
        const UF_long seed_count = solver->fast_reject_tail_seed_count;
        for (UF_long p = 0; p < seed_count; ++p) {
          const UF_long local_col = solver->fast_reject_tail_seed_cols[p];
          if (local_col >= nk ||
              !kls_mark_fast_reject_tail_col(solver, nk, local_col, mark,
                                             &queued)) {
            free(owned_parent);
            return 0;
          }
        }
        solver->stats.fast_rejected_pivoting_tail_seed_columns =
          (int64_t)seed_count;
        solver->stats.fast_rejected_pivoting_tail_row_seed_columns = 0;
        solver->stats.fast_rejected_pivoting_tail_block_seed_columns =
          available_block_seed_columns;
      } else {
        for (UF_long local_col = local_reject; local_col < nk; ++local_col) {
          if (!kls_mark_fast_reject_tail_col(solver, nk, local_col, mark,
                                             &queued)) {
            free(owned_parent);
            return 0;
          }
        }
        solver->stats.fast_rejected_pivoting_tail_seed_columns =
          (int64_t)queued;
        solver->stats.fast_rejected_pivoting_tail_row_seed_columns = 0;
        solver->stats.fast_rejected_pivoting_tail_block_seed_columns =
          available_block_seed_columns;
      }
    } else if (!kls_mark_fast_reject_tail_col(solver, nk, local_reject, mark,
                                              &queued)) {
      free(owned_parent);
      return 0;
    } else {
      solver->stats.fast_rejected_pivoting_tail_seed_columns = 1;
      solver->stats.fast_rejected_pivoting_tail_row_seed_columns = 0;
      solver->stats.fast_rejected_pivoting_tail_block_seed_columns =
        available_block_seed_columns;
    }
    if (queued == 0u) {
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
      solver, block, k1, nk, local_reject, parent, mark, columns, work);
    if (used_seed_mode == KLS_FAST_REJECT_TAIL_SEED_ROW &&
        !solver->stats.fast_rejected_pivoting_tail_topological) {
      solver->stats.fast_rejected_pivoting_tail_row_seed_columns = 0;
      solver->stats.fast_rejected_pivoting_tail_block_seed_columns =
        available_block_seed_columns;
      seed_mode = use_unfinished_seed
        ? KLS_FAST_REJECT_TAIL_SEED_UNFINISHED
        : KLS_FAST_REJECT_TAIL_SEED_SUFFIX;
      continue;
    }
    if (used_seed_mode == KLS_FAST_REJECT_TAIL_SEED_UNFINISHED &&
        !solver->stats.fast_rejected_pivoting_tail_topological) {
      seed_mode = KLS_FAST_REJECT_TAIL_SEED_SUFFIX;
      solver->stats.fast_rejected_pivoting_tail_block_seed_columns =
        available_block_seed_columns;
      continue;
    }
    *tail_cols_out = solver->fast_reject_tail_cols;
    *tail_count_out = columns;
    *tail_work_out = work;
    free(owned_parent);
    return 1;
  }
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
  KLS_SET_OPTIONAL_OUTPUT(rows_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(count_out, 0);
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

static void kls_fill_fast_reject_group_tail_stats(kls_solver *solver,
                                                  UF_long rejected_pivot) {
  if (solver == NULL || rejected_pivot >= solver->n ||
      solver->row_refactor_pattern_n != solver->n ||
      solver->row_refactor_row_group == NULL ||
      solver->row_refactor_group_ptr == NULL ||
      solver->row_refactor_group_successor_ptr == NULL ||
      (solver->row_refactor_group_dependency_edges > 0u &&
       solver->row_refactor_group_successor_groups == NULL) ||
      solver->row_refactor_group_count <= 0) {
    return;
  }

  const UF_long group_count = solver->row_refactor_group_count;
  const UF_long root_group = solver->row_refactor_row_group[rejected_pivot];
  if (root_group < 0 || root_group >= group_count) {
    return;
  }
  const size_t group_count_size = (size_t)group_count;
  if ((UF_long)group_count_size != group_count ||
      group_count_size > SIZE_MAX / sizeof(UF_long)) {
    return;
  }
  UF_long *queue =
    (UF_long *)malloc(group_count_size * sizeof(*queue));
  unsigned char *marks =
    (unsigned char *)calloc(group_count_size, sizeof(*marks));
  if (queue == NULL || marks == NULL) {
    free(queue);
    free(marks);
    return;
  }

  UF_long head = 0;
  UF_long tail = 0;
  marks[root_group] = 1u;
  queue[tail++] = root_group;
  while (head < tail) {
    const UF_long group = queue[head++];
    const UF_long begin = solver->row_refactor_group_successor_ptr[group];
    const UF_long end = solver->row_refactor_group_successor_ptr[group + 1u];
    if (begin > end ||
        end > solver->row_refactor_group_dependency_edges) {
      free(queue);
      free(marks);
      return;
    }
    for (UF_long p = begin; p < end; ++p) {
      const UF_long successor =
        solver->row_refactor_group_successor_groups[p];
      if (successor < 0 || successor >= group_count) {
        free(queue);
        free(marks);
        return;
      }
      if (!marks[successor]) {
        marks[successor] = 1u;
        queue[tail++] = successor;
      }
    }
  }

  UF_long groups = 0;
  UF_long rows = 0;
  double work = 0.0;
  for (UF_long group = root_group; group < group_count; ++group) {
    if (!marks[group]) {
      continue;
    }
    const UF_long row_begin = solver->row_refactor_group_ptr[group];
    const UF_long row_end = solver->row_refactor_group_ptr[group + 1u];
    if (row_begin > row_end || row_end > solver->n) {
      free(queue);
      free(marks);
      return;
    }
    groups++;
    rows += row_end - row_begin;
    work += kls_row_refactor_group_work(solver, group);
  }
  solver->stats.fast_rejected_group_tail_groups = (int64_t)groups;
  solver->stats.fast_rejected_group_tail_rows = (int64_t)rows;
  solver->stats.fast_rejected_group_tail_work = work;
  free(queue);
  free(marks);
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

  const UF_long global_reject = k1 + local_reject;
  const UF_long *tail_rows = NULL;
  UF_long row_tail_columns = 0;
  if (kls_build_row_refactor_tail(solver, global_reject,
                                  &tail_rows, &row_tail_columns)) {
    int block_local_tail = 1;
    for (UF_long p = 0; p < row_tail_columns; ++p) {
      if (tail_rows[p] < k1 || tail_rows[p] >= k2) {
        block_local_tail = 0;
        break;
      }
    }
    if (block_local_tail) {
      double work = 0.0;
      for (UF_long p = 0; p < row_tail_columns; ++p) {
        work += kls_row_refactor_tail_row_work(solver, tail_rows[p]);
      }

      solver->stats.fast_rejected_row_tail_columns =
        (int64_t)row_tail_columns;
      solver->stats.fast_rejected_row_tail_work = work;
      return;
    }
    solver->row_refactor_tail_count = 0;
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

static int kls_refactor_input_value(const kls_solver *solver,
                                    const double *numeric_values,
                                    const double *rs,
                                    int scale,
                                    UF_long input_pos,
                                    double *value_out) {
  if (solver == NULL || numeric_values == NULL || value_out == NULL ||
      solver->row_idx == NULL || input_pos >= solver->nnz) {
    return 0;
  }
  double value = numeric_values[input_pos];
  if (scale > 0) {
    const UF_long oldrow = solver->row_idx[input_pos];
    if (oldrow >= solver->n || rs == NULL || rs[oldrow] == 0.0) {
      return 0;
    }
    value /= rs[oldrow];
  }
  *value_out = value;
  return 1;
}

static int kls_refresh_row_refactor_offdiag_values(kls_solver *solver,
                                                   const double *numeric_values) {
  if (solver == NULL || solver->symbolic == NULL ||
      solver->symbolic->nblocks <= 1u) {
    return 1;
  }
  if (solver->lean_btf_off_input_pos != NULL &&
      numeric_values != NULL && solver->numeric != NULL &&
      solver->numeric->Offp != NULL && solver->numeric->Offx != NULL &&
      ((solver->common.scale > 0 && solver->numeric->Rs != NULL) ||
       (solver->common.scale <= 0 && solver->numeric->Rs == NULL))) {
    const UF_long offcount = solver->numeric->Offp[solver->n];
    double *restrict offx = (double *)solver->numeric->Offx;
    const uint32_t *restrict input = solver->lean_btf_off_input_pos;
    for (UF_long p = 0u; p < offcount; ++p) {
      if ((UF_long)input[p] >= solver->nnz) {
        return 0;
      }
      double value = 0.0;
      if (!kls_refactor_input_value(
            solver, numeric_values,
            solver->common.scale > 0 ? solver->numeric->Rs : NULL,
            (int)solver->common.scale, (UF_long)input[p], &value)) {
        return 0;
      }
      offx[p] = value;
    }
    return 1;
  }
  if (solver->input_to_csc == NULL && solver->row_perm == NULL &&
      solver->compact_match_offdiag_user_pos != NULL &&
      numeric_values != NULL && solver->numeric != NULL &&
      solver->numeric->Offp != NULL && solver->numeric->Offx != NULL &&
      solver->common.scale <= 0 && solver->numeric->Rs == NULL) {
    const UF_long offcount = solver->numeric->Offp[solver->n];
    double *restrict offx = (double *)solver->numeric->Offx;
    const uint16_t *restrict input =
      solver->compact_match_offdiag_user_pos;
    for (UF_long p = 0u; p < offcount; ++p) {
      if ((UF_long)input[p] >= solver->nnz) {
        return 0;
      }
      offx[p] = numeric_values[(UF_long)input[p]];
    }
    return 1;
  }
  if (numeric_values == NULL || solver->numeric == NULL ||
      solver->refactor_col_ptr == NULL ||
      solver->refactor_block_start == NULL ||
      solver->refactor_input_pos == NULL ||
      (solver->lean_user_values_active &&
       solver->compact_match_offdiag_user_pos == NULL) ||
      solver->numeric->Offp == NULL || solver->numeric->Offx == NULL ||
      (solver->common.scale > 0 && solver->numeric->Rs == NULL) ||
      (solver->common.scale <= 0 && solver->numeric->Rs != NULL)) {
    return 0;
  }

  double *offx = (double *)solver->numeric->Offx;
  for (UF_long k = 0; k < solver->n; ++k) {
    if (solver->refactor_block_start[k] < solver->refactor_col_ptr[k] ||
        solver->refactor_block_start[k] > solver->refactor_col_ptr[k + 1u]) {
      return 0;
    }
    UF_long poff = solver->numeric->Offp[k];
    const UF_long poff_end = solver->numeric->Offp[k + 1u];
    for (UF_long p = solver->refactor_col_ptr[k];
         p < solver->refactor_block_start[k]; ++p) {
      if (poff >= poff_end) {
        return 0;
      }
      double value = 0.0;
      if (solver->lean_user_values_active) {
        const UF_long input_pos =
          solver->compact_match_offdiag_user_pos[poff];
        if (input_pos >= solver->nnz) {
          return 0;
        }
        value = numeric_values[input_pos];
      } else {
        if (!kls_refactor_input_value(
              solver, numeric_values,
              solver->common.scale > 0 ? solver->numeric->Rs : NULL,
              (int)solver->common.scale, solver->refactor_input_pos[p],
              &value)) {
          return 0;
        }
      }
      offx[poff++] = value;
    }
    if (poff != poff_end) {
      return 0;
    }
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
    const UF_long input_pos = solver->row_refactor_input_pos[p];
    if (!kls_row_refactor_touch(input_col, touched, touched_cols,
                                &touched_count, solver->n)) {
      ok = 0;
      break;
    }
    if (!kls_refactor_input_value(
          solver, numeric_values,
          solver->common.scale > 0 ? solver->numeric->Rs : NULL,
          (int)solver->common.scale, input_pos, x + input_col)) {
      ok = 0;
      break;
    }
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
  solver->stats.fast_rejected_group_tail_groups = 0;
  solver->stats.fast_rejected_group_tail_rows = 0;
  solver->stats.fast_rejected_group_tail_work = 0.0;
  solver->stats.fast_rejected_etree_columns = 0;
  solver->stats.fast_rejected_etree_work = 0.0;
  solver->stats.fast_rejected_pivoting_tail_columns = 0;
  solver->stats.fast_rejected_pivoting_tail_work = 0.0;
  solver->stats.fast_rejected_pivoting_tail_first = -1;
  solver->stats.fast_rejected_pivoting_tail_last = -1;
  solver->stats.fast_rejected_pivoting_tail_contains_reject = 0;
  solver->stats.fast_rejected_pivoting_tail_topological = 0;
  solver->stats.fast_rejected_pivoting_tail_seed_columns =
    solver->fast_reject_tail_seed_valid
      ? (int64_t)solver->fast_reject_tail_seed_count : 0;
  solver->stats.fast_rejected_pivoting_tail_row_seed_columns = 0;
  solver->stats.fast_rejected_pivoting_tail_block_seed_columns =
    solver->fast_reject_tail_seed_valid &&
        solver->fast_reject_tail_seed_block_suffix
      ? (int64_t)solver->fast_reject_tail_seed_count : 0;
  solver->stats.fast_rejected_pivoting_tail_contiguous = 0;
  solver->stats.fast_rejected_pivoting_tail_suffix_exact = 0;
  solver->stats.fast_rejected_pivoting_tail_gap_columns = 0;
  solver->stats.fast_rejected_pivoting_tail_suffix_overcompute_columns = 0;
  solver->stats.fast_rejected_pivoting_tail_suffix_overcompute_work = 0.0;
  solver->stats.fast_rejected_pivoting_tail_etree_edges = 0;
  solver->stats.fast_rejected_pivoting_tail_etree_roots = 0;
  solver->stats.fast_rejected_pivoting_tail_etree_leaves = 0;
  solver->stats.fast_rejected_pivoting_tail_etree_max_fanout = 0;
  solver->stats.fast_rejected_pivoting_tail_etree_levels = 0;
  solver->stats.fast_rejected_pivoting_tail_etree_max_width = 0;
  solver->fast_reject_tail_count = 0;
  solver->fast_reject_tail_plan_mark = 0u;
  solver->fast_reject_tail_plan_block = KLS_KLU_EMPTY;
  solver->fast_reject_tail_plan_k1 = KLS_KLU_EMPTY;
  solver->fast_reject_tail_plan_nk = 0;

  kls_fill_fast_reject_row_tail_stats(solver, block, k1, k2,
                                      local_reject);
  kls_fill_fast_reject_group_tail_stats(solver, rejected_pivot);
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
  KLS_SET_OPTIONAL_OUTPUT(prefix_lup_out, 0u);
  KLS_SET_OPTIONAL_OUTPUT(prefix_lnz_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(prefix_unz_out, 0);
  if (solver == NULL || solver->numeric == NULL ||
      solver->numeric->LUbx == NULL || solver->numeric->LUsize == NULL ||
      solver->numeric->Llen == NULL || solver->numeric->Ulen == NULL ||
      solver->numeric->Lip == NULL || solver->numeric->Uip == NULL ||
      final_pblock == NULL || live_p == NULL || live_pinv == NULL ||
      final_pinv == NULL || lpend == NULL || prefix_cols > nk ||
      k1 > solver->n || nk > solver->n - k1 ||
      block >= solver->numeric->nblocks) {
    return 0;
  }

  if (prefix_cols > 0u &&
      (solver->numeric->LUbx[block] == NULL ||
       solver->numeric->LUsize[block] == 0u)) {
    return 0;
  }
  const size_t lusize =
    prefix_cols > 0u ? solver->numeric->LUsize[block] : 0u;
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

  double *lu = prefix_cols > 0u
    ? (double *)solver->numeric->LUbx[block] : NULL;
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

static int kls_tail_mapped_column_range(const kls_solver *solver,
                                        UF_long global_col,
                                        UF_long *begin_out,
                                        UF_long *end_out) {
  if (solver == NULL || begin_out == NULL || end_out == NULL ||
      solver->refactor_col_ptr == NULL ||
      solver->refactor_row_idx == NULL ||
      solver->refactor_input_pos == NULL ||
      solver->symbolic == NULL ||
      global_col >= solver->n ||
      solver->refactor_col_ptr[global_col] >
        solver->refactor_col_ptr[global_col + 1u] ||
      solver->refactor_col_ptr[global_col + 1u] > solver->nnz) {
    return 0;
  }

  const UF_long col_begin = solver->refactor_col_ptr[global_col];
  const UF_long col_end = solver->refactor_col_ptr[global_col + 1u];
  UF_long block_begin = col_begin;
  if (solver->symbolic->nblocks > 1u) {
    if (solver->refactor_block_start == NULL ||
        solver->refactor_block_start[global_col] < col_begin ||
        solver->refactor_block_start[global_col] > col_end) {
      return 0;
    }
    block_begin = solver->refactor_block_start[global_col];
  }

  *begin_out = block_begin;
  *end_out = col_end;
  return 1;
}

static int kls_tail_mapped_block_is_valid(const kls_solver *solver,
                                          UF_long k1,
                                          UF_long nk,
                                          const UF_long *psinv) {
  if (solver == NULL || k1 > solver->n || nk > solver->n - k1 ||
      solver->row_idx == NULL || psinv == NULL) {
    return 0;
  }
  for (UF_long kk = 0; kk < nk; ++kk) {
    const UF_long global_col = k1 + kk;
    UF_long begin = 0;
    UF_long end = 0;
    if (!kls_tail_mapped_column_range(solver, global_col, &begin, &end)) {
      return 0;
    }
    for (UF_long p = begin; p < end; ++p) {
      const UF_long global_row = solver->refactor_row_idx[p];
      const UF_long input_pos = solver->refactor_input_pos[p];
      if (global_row < k1 || global_row >= k1 + nk ||
          input_pos >= solver->nnz ||
          solver->row_idx[input_pos] >= solver->n ||
          psinv[solver->row_idx[input_pos]] != global_row) {
        return 0;
      }
    }
  }
  return 1;
}

static Int kls_tail_lsolve_symbolic_mapped(const kls_solver *solver,
                                           Int n,
                                           Int k,
                                           Int pinv[],
                                           Int stack[],
                                           Int flag[],
                                           Int lpend[],
                                           Int ap_pos[],
                                           Unit lu[],
                                           Int lup,
                                           Int llen[],
                                           Int lip[],
                                           Int k1) {
  Int top = n;
  Int l_length = 0;
  Int *lik = (Int *)(lu + lup);
  const UF_long global_col = (UF_long)(k + k1);
  UF_long begin = 0;
  UF_long end = 0;
  if (!kls_tail_mapped_column_range(solver, global_col, &begin, &end)) {
    return -1;
  }

  for (UF_long p = begin; p < end; ++p) {
    const UF_long global_row = solver->refactor_row_idx[p];
    if (global_row < (UF_long)k1 ||
        global_row >= (UF_long)k1 + (UF_long)n) {
      return -1;
    }
    const Int i = (Int)(global_row - (UF_long)k1);
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
  Int poff = 0;
  if (offp != NULL) {
    poff = offp[kglobal];
    if (poff < 0 || poff > nzoff) {
      return 0;
    }
  } else if (offi != NULL || offx != NULL) {
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
      if (offp != NULL) {
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
      }
    } else if (i >= n) {
      return 0;
    } else {
      x[i] = aik;
    }
  }
  if (offp != NULL) {
    offp[kglobal + 1] = poff;
  }
  return 1;
}

static int kls_tail_construct_column_mapped(const kls_solver *solver,
                                            const double *numeric_values,
                                            Int k,
                                            Entry x[],
                                            Int n,
                                            Int k1,
                                            double rs[],
                                            Int scale) {
  if (solver == NULL || numeric_values == NULL || x == NULL ||
      solver->row_idx == NULL) {
    return 0;
  }
  const UF_long global_col = (UF_long)(k + k1);
  UF_long begin = 0;
  UF_long end = 0;
  if (!kls_tail_mapped_column_range(solver, global_col, &begin, &end)) {
    return 0;
  }

  for (UF_long p = begin; p < end; ++p) {
    const UF_long global_row = solver->refactor_row_idx[p];
    const UF_long input_pos = solver->refactor_input_pos[p];
    if (global_row < (UF_long)k1 ||
        global_row >= (UF_long)k1 + (UF_long)n ||
        input_pos >= solver->nnz) {
      return 0;
    }
    const UF_long oldrow = solver->row_idx[input_pos];
    if (oldrow >= solver->n) {
      return 0;
    }
    Entry aik = numeric_values[input_pos];
    if (scale > 0) {
      if (rs == NULL || rs[oldrow] == 0.0) {
        return 0;
      }
      SCALE_DIV(aik, rs[oldrow]);
    }
    x[global_row - (UF_long)k1] = aik;
  }
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
  if (new_lu == NULL || final_pblock == NULL || lup_out == NULL ||
      (prefix_cols > 0u && old_lu == NULL)) {
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

static int kls_copy_preserved_lu_columns(Unit **new_lu_io,
                                         size_t *lusize_io,
                                         const Unit *old_lu,
                                         size_t old_lusize,
                                         const UF_long *old_lip,
                                         const UF_long *old_llen,
                                         const UF_long *old_uip,
                                         const UF_long *old_ulen,
                                         const Entry *old_udiag,
                                         const UF_long *old_pblock,
                                         const UF_long *final_pinv,
                                         UF_long nk,
                                         UF_long begin,
                                         UF_long end,
                                         UF_long *lip,
                                         UF_long *llen,
                                         UF_long *uip,
                                         UF_long *ulen,
                                         Entry *udiag,
                                         size_t *lup_io,
                                         UF_long *lnz_io,
                                         UF_long *unz_io,
                                         TRILINOS_KLU_common *common) {
  if (new_lu_io == NULL || *new_lu_io == NULL || lusize_io == NULL ||
      old_lu == NULL || old_lip == NULL || old_llen == NULL ||
      old_uip == NULL || old_ulen == NULL || old_udiag == NULL ||
      old_pblock == NULL || final_pinv == NULL || lip == NULL ||
      llen == NULL || uip == NULL || ulen == NULL || udiag == NULL ||
      lup_io == NULL || lnz_io == NULL || unz_io == NULL ||
      common == NULL || begin > end || end > nk) {
    return 0;
  }

  Unit *new_lu = *new_lu_io;
  size_t lusize = *lusize_io;
  size_t lup = *lup_io;
  for (UF_long k = begin; k < end; ++k) {
    const UF_long l_len = old_llen[k];
    const UF_long u_len = old_ulen[k];
    const UF_long src_lip = old_lip[k];
    const UF_long src_uip = old_uip[k];
    const size_t l_units = kls_klu_units_for_indices(l_len);
    const size_t u_units = kls_klu_units_for_indices(u_len);
    if ((size_t)src_lip > old_lusize ||
        l_units > old_lusize - (size_t)src_lip ||
        l_len > (UF_long)(old_lusize - (size_t)src_lip - l_units) ||
        (size_t)src_uip > old_lusize ||
        u_units > old_lusize - (size_t)src_uip ||
        u_len > (UF_long)(old_lusize - (size_t)src_uip - u_units)) {
      return 0;
    }
    const size_t need =
      l_units + (size_t)l_len + u_units + (size_t)u_len;
    if (need > SIZE_MAX - lup) {
      return 0;
    }
    if (lup + need > lusize) {
      const double grown_double =
        common->memgrow * (double)lusize + 2.0 * (double)nk +
        (double)need + 1.0;
      if (INT_OVERFLOW(grown_double)) {
        common->status = TRILINOS_KLU_TOO_LARGE;
        return 0;
      }
      size_t grown = (size_t)grown_double;
      if (grown < lup + need) {
        grown = lup + need;
      }
      Unit *grown_lu =
        (Unit *)TRILINOS_KLU_realloc(grown, lusize, sizeof(Unit), new_lu,
                                     common);
      if (common->status == TRILINOS_KLU_OUT_OF_MEMORY ||
          grown_lu == NULL) {
        return 0;
      }
      new_lu = grown_lu;
      lusize = grown;
      common->nrealloc++;
    }

    lip[k] = (UF_long)lup;
    llen[k] = l_len;
    UF_long *dst_li = (UF_long *)(new_lu + lup);
    Entry *dst_lx = (Entry *)(new_lu + lup + l_units);
    Int *src_li = (Int *)((Unit *)old_lu + src_lip);
    Entry *src_lx = (Entry *)((Unit *)old_lu + src_lip + l_units);
    for (UF_long p = 0; p < l_len; ++p) {
      const UF_long old_final_order = src_li[p];
      if (old_final_order >= nk) {
        return 0;
      }
      const UF_long local_row = old_pblock[old_final_order];
      if (local_row >= nk || final_pinv[local_row] <= k ||
          final_pinv[local_row] >= nk) {
        return 0;
      }
      dst_li[p] = local_row;
      dst_lx[p] = src_lx[p];
    }
    lup += l_units + (size_t)l_len;

    uip[k] = (UF_long)lup;
    ulen[k] = u_len;
    UF_long *dst_ui = (UF_long *)(new_lu + lup);
    Entry *dst_ux = (Entry *)(new_lu + lup + u_units);
    Int *src_ui = (Int *)((Unit *)old_lu + src_uip);
    Entry *src_ux = (Entry *)((Unit *)old_lu + src_uip + u_units);
    for (UF_long p = 0; p < u_len; ++p) {
      if (src_ui[p] >= (Int)k) {
        return 0;
      }
      dst_ui[p] = (UF_long)src_ui[p];
      dst_ux[p] = src_ux[p];
    }
    lup += u_units + (size_t)u_len;
    udiag[k] = old_udiag[k];
    *lnz_io += l_len + 1u;
    *unz_io += u_len + 1u;
  }

  *new_lu_io = new_lu;
  *lusize_io = lusize;
  *lup_io = lup;
  return 1;
}

static int kls_copy_finished_gap_lu_column(
  Unit **new_lu_io,
  size_t *lusize_io,
  const Unit *old_lu,
  size_t old_lusize,
  const UF_long *old_lip,
  const UF_long *old_llen,
  const UF_long *old_uip,
  const UF_long *old_ulen,
  const Entry *old_udiag,
  const UF_long *old_pblock,
  const unsigned char *restart_tail_mask,
  const UF_long *live_pinv,
  UF_long nk,
  UF_long k,
  UF_long *lip,
  UF_long *llen,
  UF_long *uip,
  UF_long *ulen,
  Entry *udiag,
  size_t *lup_io,
  UF_long *lnz_io,
  UF_long *unz_io,
  TRILINOS_KLU_common *common) {
  if (new_lu_io == NULL || *new_lu_io == NULL || lusize_io == NULL ||
      old_lu == NULL || old_lip == NULL || old_llen == NULL ||
      old_uip == NULL || old_ulen == NULL || old_udiag == NULL ||
      old_pblock == NULL || live_pinv == NULL || lip == NULL ||
      llen == NULL || uip == NULL || ulen == NULL || udiag == NULL ||
      lup_io == NULL || lnz_io == NULL || unz_io == NULL ||
      common == NULL || k >= nk || old_pblock[k] >= nk ||
      live_pinv[old_pblock[k]] >= 0) {
    return 0;
  }

  const UF_long l_len = old_llen[k];
  const UF_long u_len = old_ulen[k];
  const UF_long src_lip = old_lip[k];
  const UF_long src_uip = old_uip[k];
  const size_t l_units = kls_klu_units_for_indices(l_len);
  const size_t u_units = kls_klu_units_for_indices(u_len);
  if ((size_t)src_lip > old_lusize ||
      l_units > old_lusize - (size_t)src_lip ||
      l_len > (UF_long)(old_lusize - (size_t)src_lip - l_units) ||
      (size_t)src_uip > old_lusize ||
      u_units > old_lusize - (size_t)src_uip ||
      u_len > (UF_long)(old_lusize - (size_t)src_uip - u_units)) {
    return 0;
  }

  const Int *src_ui = (const Int *)((const Unit *)old_lu + src_uip);
  const Entry *src_ux =
    (const Entry *)((const Unit *)old_lu + src_uip + u_units);
  for (UF_long p = 0; p < u_len; ++p) {
    const Int dep = src_ui[p];
    if (dep < 0 || (UF_long)dep >= k) {
      return 0;
    }
    if (restart_tail_mask != NULL && restart_tail_mask[dep]) {
      return 0;
    }
  }

  const Int *src_li = (const Int *)((const Unit *)old_lu + src_lip);
  const Entry *src_lx =
    (const Entry *)((const Unit *)old_lu + src_lip + l_units);
  for (UF_long p = 0; p < l_len; ++p) {
    const Int final_order = src_li[p];
    if (final_order <= (Int)k || final_order < 0 ||
        (UF_long)final_order >= nk) {
      return 0;
    }
    const UF_long local_row = old_pblock[final_order];
    if (local_row >= nk || live_pinv[local_row] >= 0) {
      return 0;
    }
  }

  const size_t need =
    l_units + (size_t)l_len + u_units + (size_t)u_len;
  size_t lup = *lup_io;
  if (need > SIZE_MAX - lup) {
    return 0;
  }
  Unit *new_lu = *new_lu_io;
  size_t lusize = *lusize_io;
  if (lup + need > lusize) {
    const double grown_double =
      common->memgrow * (double)lusize + 2.0 * (double)nk +
      (double)need + 1.0;
    if (INT_OVERFLOW(grown_double)) {
      common->status = TRILINOS_KLU_TOO_LARGE;
      return 0;
    }
    size_t grown = (size_t)grown_double;
    if (grown < lup + need) {
      grown = lup + need;
    }
    Unit *grown_lu =
      (Unit *)TRILINOS_KLU_realloc(grown, lusize, sizeof(Unit), new_lu,
                                   common);
    if (common->status == TRILINOS_KLU_OUT_OF_MEMORY ||
        grown_lu == NULL) {
      return 0;
    }
    new_lu = grown_lu;
    lusize = grown;
    common->nrealloc++;
  }

  lip[k] = (UF_long)lup;
  llen[k] = l_len;
  Int *dst_li = (Int *)(new_lu + lup);
  Entry *dst_lx = (Entry *)(new_lu + lup + l_units);
  for (UF_long p = 0; p < l_len; ++p) {
    const UF_long local_row = old_pblock[src_li[p]];
    dst_li[p] = (Int)local_row;
    dst_lx[p] = src_lx[p];
  }
  lup += l_units + (size_t)l_len;

  uip[k] = (UF_long)lup;
  ulen[k] = u_len;
  Int *dst_ui = (Int *)(new_lu + lup);
  Entry *dst_ux = (Entry *)(new_lu + lup + u_units);
  for (UF_long p = 0; p < u_len; ++p) {
    dst_ui[p] = src_ui[p];
    dst_ux[p] = src_ux[p];
  }
  lup += u_units + (size_t)u_len;
  udiag[k] = old_udiag[k];

  *new_lu_io = new_lu;
  *lusize_io = lusize;
  *lup_io = lup;
  *lnz_io += l_len + 1u;
  *unz_io += u_len + 1u;
  return 1;
}

static int kls_finished_gap_lu_column_depends_on_tail(
  const Unit *old_lu,
  size_t old_lusize,
  const UF_long *old_lip,
  const UF_long *old_llen,
  const UF_long *old_uip,
  const UF_long *old_ulen,
  const UF_long *old_pblock,
  const unsigned char *restart_tail_mask,
  const UF_long *live_pinv,
  UF_long nk,
  UF_long k) {
  if (old_lu == NULL || old_lip == NULL || old_llen == NULL ||
      old_uip == NULL || old_ulen == NULL || old_pblock == NULL ||
      restart_tail_mask == NULL || live_pinv == NULL ||
      k >= nk || old_pblock[k] >= nk) {
    return 0;
  }
  const UF_long l_len = old_llen[k];
  const UF_long u_len = old_ulen[k];
  const UF_long src_lip = old_lip[k];
  const UF_long src_uip = old_uip[k];
  const size_t l_units = kls_klu_units_for_indices(l_len);
  const size_t u_units = kls_klu_units_for_indices(u_len);
  if ((size_t)src_lip > old_lusize ||
      l_units > old_lusize - (size_t)src_lip ||
      l_len > (UF_long)(old_lusize - (size_t)src_lip - l_units) ||
      (size_t)src_uip > old_lusize ||
      u_units > old_lusize - (size_t)src_uip ||
      u_len > (UF_long)(old_lusize - (size_t)src_uip - u_units)) {
    return 0;
  }

  int depends = live_pinv[old_pblock[k]] >= 0;
  const Int *src_ui = (const Int *)((const Unit *)old_lu + src_uip);
  for (UF_long p = 0; p < u_len; ++p) {
    const Int dep = src_ui[p];
    if (dep < 0 || (UF_long)dep >= k) {
      return 0;
    }
    if (restart_tail_mask[dep]) {
      depends = 1;
    }
  }

  const Int *src_li = (const Int *)((const Unit *)old_lu + src_lip);
  for (UF_long p = 0; p < l_len; ++p) {
    const Int final_order = src_li[p];
    if (final_order <= (Int)k || final_order < 0 ||
        (UF_long)final_order >= nk) {
      return 0;
    }
    const UF_long local_row = old_pblock[final_order];
    if (local_row >= nk) {
      return 0;
    }
    if (live_pinv[local_row] >= 0) {
      depends = 1;
    }
  }
  return depends;
}

static UF_long kls_preferred_pivot_tail_restart_end(const kls_solver *solver,
                                                    UF_long k1,
                                                    UF_long nk,
                                                    UF_long local_reject) {
  if (solver == NULL || local_reject >= nk ||
      solver->stats.fast_rejected_pivoting_tail_suffix_exact ||
      !solver->stats.fast_rejected_pivoting_tail_topological ||
      solver->stats.fast_rejected_pivoting_tail_columns <= 0 ||
      solver->stats.fast_rejected_pivoting_tail_first !=
        (int64_t)(k1 + local_reject) ||
      solver->stats.fast_rejected_pivoting_tail_last <
        (int64_t)(k1 + local_reject) ||
      solver->stats.fast_rejected_pivoting_tail_last >=
        (int64_t)(k1 + nk - 1u)) {
    return nk;
  }
  const UF_long restart_end =
    (UF_long)solver->stats.fast_rejected_pivoting_tail_last - k1 + 1u;
  return restart_end > local_reject && restart_end < nk ? restart_end : nk;
}

static int kls_build_pivot_tail_restart_mask_range(
  const kls_solver *solver,
  UF_long k1,
  UF_long nk,
  UF_long required_local_col,
  UF_long restart_begin,
  UF_long restart_end,
  unsigned char *mask) {
  if (solver == NULL || mask == NULL || required_local_col >= nk ||
      restart_begin > required_local_col || restart_end <= required_local_col ||
      restart_end > nk ||
      solver->fast_reject_tail_cols == NULL ||
      solver->fast_reject_tail_count == 0u) {
    return 0;
  }
  memset(mask, 0, (size_t)nk * sizeof(*mask));
  for (UF_long p = 0; p < solver->fast_reject_tail_count; ++p) {
    const UF_long global_col = solver->fast_reject_tail_cols[p];
    if (global_col < k1 || global_col >= k1 + nk) {
      return 0;
    }
    const UF_long local_col = global_col - k1;
    if (local_col < restart_begin || local_col >= restart_end) {
      return 0;
    }
    mask[local_col] = 1u;
  }
  return mask[required_local_col] && mask[restart_begin] &&
         mask[restart_end - 1u];
}

static int kls_build_pivot_tail_restart_mask(const kls_solver *solver,
                                             UF_long k1,
                                             UF_long nk,
                                             UF_long local_reject,
                                             UF_long restart_end,
                                             unsigned char *mask) {
  return kls_build_pivot_tail_restart_mask_range(
    solver, k1, nk, local_reject, local_reject, restart_end, mask);
}

static int kls_pivot_tail_restart_mask_matches_etree_forest(
  kls_solver *solver,
  UF_long block,
  UF_long k1,
  UF_long nk,
  UF_long local_reject,
  UF_long restart_begin,
  UF_long restart_end,
  const unsigned char *mask) {
  if (solver == NULL || mask == NULL || local_reject >= nk ||
      restart_begin > local_reject || restart_end <= local_reject ||
      restart_end > nk ||
      solver->fast_reject_tail_cols == NULL ||
      solver->fast_reject_tail_parent == NULL ||
      solver->fast_reject_tail_child_count == NULL ||
      solver->fast_reject_tail_level == NULL ||
      solver->fast_reject_tail_count == 0u ||
      solver->fast_reject_tail_plan_mark == 0u ||
      solver->fast_reject_tail_plan_block != block ||
      solver->fast_reject_tail_plan_k1 != k1 ||
      solver->fast_reject_tail_plan_nk != nk ||
      solver->stats.fast_rejected_pivoting_tail_columns <= 0 ||
      !solver->stats.fast_rejected_pivoting_tail_contains_reject ||
      !solver->stats.fast_rejected_pivoting_tail_topological) {
    return 0;
  }
  const UF_long expected_columns =
    (UF_long)solver->stats.fast_rejected_pivoting_tail_columns;
  if (expected_columns != solver->fast_reject_tail_count) {
    return 0;
  }
  for (UF_long p = 0; p < solver->fast_reject_tail_count; ++p) {
    const UF_long global_col = solver->fast_reject_tail_cols[p];
    if (global_col < k1 || global_col >= k1 + nk) {
      return 0;
    }
    const UF_long local_col = global_col - k1;
    if (local_col < restart_begin || local_col >= restart_end ||
        mask[local_col] == 0u) {
      return 0;
    }
  }

  UF_long columns = 0;
  UF_long edges = 0;
  UF_long roots = 0;
  UF_long leaves = 0;
  UF_long max_fanout = 0;
  UF_long max_level = 0;
  int contains_reject = 0;
  int ok = 1;
  for (UF_long local_col = 0; local_col < nk; ++local_col) {
    if (mask[local_col] == 0u) {
      continue;
    }
    if (local_col < restart_begin || local_col >= restart_end) {
      ok = 0;
      break;
    }
    columns++;
    if (local_col == local_reject) {
      contains_reject = 1;
    }
    const UF_long next = solver->fast_reject_tail_parent[local_col];
    if (next != KLS_KLU_EMPTY && next < nk &&
        mask[next] != 0u) {
      if (next <= local_col) {
        ok = 0;
        break;
      }
      edges++;
    } else {
      roots++;
    }
    if (solver->fast_reject_tail_child_count[local_col] == 0u) {
      leaves++;
    }
    if (solver->fast_reject_tail_child_count[local_col] > max_fanout) {
      max_fanout = solver->fast_reject_tail_child_count[local_col];
    }
    if (solver->fast_reject_tail_level[local_col] > max_level) {
      max_level = solver->fast_reject_tail_level[local_col];
    }
  }

  UF_long max_width = 0;
  const UF_long level_count = ok && columns > 0u ? max_level + 1u : 0u;
  UF_long *level_width = ok && level_count > 0u
    ? (UF_long *)calloc((size_t)level_count, sizeof(*level_width)) : NULL;
  if (ok && level_count > 0u && level_width == NULL) {
    ok = 0;
  }
  if (ok) {
    for (UF_long local_col = 0; local_col < nk; ++local_col) {
      if (mask[local_col] == 0u) {
        continue;
      }
      const UF_long level = solver->fast_reject_tail_level[local_col];
      if (level >= level_count) {
        ok = 0;
        break;
      }
      level_width[level]++;
      if (level_width[level] > max_width) {
        max_width = level_width[level];
      }
    }
    ok = contains_reject &&
         columns == expected_columns &&
         edges + roots == columns &&
         edges == (UF_long)
           solver->stats.fast_rejected_pivoting_tail_etree_edges &&
         roots == (UF_long)
           solver->stats.fast_rejected_pivoting_tail_etree_roots &&
         leaves == (UF_long)
           solver->stats.fast_rejected_pivoting_tail_etree_leaves &&
         max_fanout == (UF_long)
           solver->stats.fast_rejected_pivoting_tail_etree_max_fanout &&
         level_count == (UF_long)
           solver->stats.fast_rejected_pivoting_tail_etree_levels &&
         max_width == (UF_long)
           solver->stats.fast_rejected_pivoting_tail_etree_max_width;
  }
  free(level_width);
  return ok;
}

static int kls_build_pivot_tail_restart_order_range(
  const kls_solver *solver,
  UF_long k1,
  UF_long nk,
  UF_long required_local_col,
  UF_long restart_begin,
  UF_long restart_end,
  const unsigned char *mask,
  UF_long *order,
  UF_long *order_count_out) {
  KLS_SET_OPTIONAL_OUTPUT(order_count_out, 0);
  if (solver == NULL || mask == NULL || order == NULL ||
      order_count_out == NULL || required_local_col >= nk ||
      restart_begin > required_local_col ||
      restart_end <= required_local_col || restart_end > nk ||
      solver->fast_reject_tail_cols == NULL ||
      solver->fast_reject_tail_parent == NULL ||
      solver->fast_reject_tail_count == 0u ||
      solver->stats.fast_rejected_pivoting_tail_columns <= 0 ||
      !solver->stats.fast_rejected_pivoting_tail_contains_reject ||
      !solver->stats.fast_rejected_pivoting_tail_topological) {
    return 0;
  }
  const UF_long expected =
    (UF_long)solver->stats.fast_rejected_pivoting_tail_columns;
  if (expected == 0u || expected > nk ||
      expected != solver->fast_reject_tail_count) {
    return 0;
  }
  UF_long *rank = (UF_long *)malloc((size_t)nk * sizeof(*rank));
  if (rank == NULL) {
    return 0;
  }
  for (UF_long k = 0; k < nk; ++k) {
    rank[k] = KLS_KLU_EMPTY;
  }

  UF_long count = 0;
  UF_long previous = KLS_KLU_EMPTY;
  int contains_required = 0;
  int ok = 1;
  for (UF_long p = 0; p < solver->fast_reject_tail_count; ++p) {
    const UF_long global_col = solver->fast_reject_tail_cols[p];
    if (global_col < k1 || global_col >= k1 + nk) {
      ok = 0;
      break;
    }
    const UF_long local_col = global_col - k1;
    if (local_col < restart_begin || local_col >= restart_end ||
        mask[local_col] == 0u || rank[local_col] != KLS_KLU_EMPTY ||
        (previous != KLS_KLU_EMPTY && local_col <= previous) ||
        count >= expected) {
      ok = 0;
      break;
    }
    order[count] = local_col;
    rank[local_col] = count;
    count++;
    previous = local_col;
    if (local_col == required_local_col) {
      contains_required = 1;
    }
  }
  if (ok && (!contains_required || count != expected)) {
    ok = 0;
  }
  for (UF_long local_col = restart_begin; ok && local_col < restart_end;
       ++local_col) {
    if (mask[local_col] != 0u &&
        rank[local_col] == KLS_KLU_EMPTY) {
      ok = 0;
      break;
    }
  }
  for (UF_long p = 0; ok && p < count; ++p) {
    const UF_long local_col = order[p];
    const UF_long next = solver->fast_reject_tail_parent[local_col];
    if (next == KLS_KLU_EMPTY || next == local_col ||
        next >= nk || mask[next] == 0u) {
      continue;
    }
    if (rank[next] == KLS_KLU_EMPTY || rank[next] <= p) {
      ok = 0;
      break;
    }
  }
  free(rank);
  if (!ok) {
    return 0;
  }
  *order_count_out = count;
  return 1;
}

static int kls_pivot_tail_plan_local_range(
  const kls_solver *solver,
  UF_long k1,
  UF_long nk,
  UF_long local_reject,
  UF_long *begin_out,
  UF_long *end_out,
  UF_long *columns_out,
  int *contiguous_out,
  int *suffix_exact_out) {
  KLS_SET_OPTIONAL_OUTPUT(begin_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(end_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(columns_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(contiguous_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(suffix_exact_out, 0);
  if (solver == NULL || local_reject >= nk ||
      solver->stats.fast_rejected_pivoting_tail_columns <= 0 ||
      !solver->stats.fast_rejected_pivoting_tail_contains_reject ||
      !solver->stats.fast_rejected_pivoting_tail_topological ||
      solver->stats.fast_rejected_pivoting_tail_first < (int64_t)k1 ||
      solver->stats.fast_rejected_pivoting_tail_last <
        solver->stats.fast_rejected_pivoting_tail_first ||
      solver->stats.fast_rejected_pivoting_tail_last >=
        (int64_t)(k1 + nk) ||
      solver->stats.fast_rejected_pivoting_tail_first >
        (int64_t)(k1 + local_reject) ||
      solver->stats.fast_rejected_pivoting_tail_last <
        (int64_t)(k1 + local_reject)) {
    return 0;
  }
  const UF_long begin =
    (UF_long)solver->stats.fast_rejected_pivoting_tail_first - k1;
  const UF_long end =
    (UF_long)solver->stats.fast_rejected_pivoting_tail_last - k1 + 1u;
  const UF_long columns =
    (UF_long)solver->stats.fast_rejected_pivoting_tail_columns;
  if (begin > local_reject || end <= local_reject || end > nk ||
      columns == 0u || columns > end - begin) {
    return 0;
  }
  if (begin_out != NULL) {
    *begin_out = begin;
  }
  if (end_out != NULL) {
    *end_out = end;
  }
  if (columns_out != NULL) {
    *columns_out = columns;
  }
  if (contiguous_out != NULL) {
    *contiguous_out =
      solver->stats.fast_rejected_pivoting_tail_contiguous ? 1 : 0;
  }
  if (suffix_exact_out != NULL) {
    *suffix_exact_out =
      solver->stats.fast_rejected_pivoting_tail_suffix_exact ? 1 : 0;
  }
  return 1;
}

static int kls_try_pivot_tail_restart_rejected_block(
  kls_solver *solver,
  double *numeric_values,
  UF_long block,
  UF_long k1,
  UF_long nk,
  UF_long local_reject,
  UF_long restart_end,
  const unsigned char *restart_tail_mask,
  const UF_long *old_pblock,
  const UF_long *psinv,
  Unit **new_lu_out,
  size_t *new_size_out,
  UF_long *lnz_block_out,
  UF_long *unz_block_out,
  UF_long *pblock,
  UF_long *skipped_columns_out,
  double *skipped_work_out,
  int *exact_mask_used_out) {
  KLS_SET_OPTIONAL_OUTPUT(new_lu_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(new_size_out, 0u);
  KLS_SET_OPTIONAL_OUTPUT(lnz_block_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(unz_block_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(skipped_columns_out, 0);
  if (skipped_work_out != NULL) {
    *skipped_work_out = 0.0;
  }
  KLS_SET_OPTIONAL_OUTPUT(exact_mask_used_out, 0);
  if (restart_end == KLS_KLU_EMPTY) {
    restart_end = nk;
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
      skipped_columns_out == NULL || skipped_work_out == NULL ||
      exact_mask_used_out == NULL ||
      (local_reject > 0u && !reusable_prefix_state) ||
      local_reject >= nk || restart_end <= local_reject ||
      restart_end > nk ||
      block >= solver->numeric->nblocks ||
      solver->numeric->LUbx == NULL ||
      solver->numeric->LUsize == NULL ||
      solver->numeric->Xwork == NULL ||
      solver->numeric->Llen == NULL || solver->numeric->Ulen == NULL ||
      solver->numeric->Lip == NULL || solver->numeric->Uip == NULL ||
      solver->numeric->Udiag == NULL || solver->numeric->Offp == NULL ||
      solver->numeric->Offi == NULL || solver->numeric->Offx == NULL) {
    return 0;
  }
  if ((scaled && solver->numeric->Rs == NULL) ||
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

  size_t old_lusize = solver->numeric->LUsize[block];
  Unit *old_lu = (Unit *)solver->numeric->LUbx[block];
  if (old_lu == NULL || old_lusize == 0u) {
    if (local_reject > 0u) {
      free(scratch);
      return 0;
    }
    const double lsize =
      solver->symbolic->Lnz[block] < 0.0
        ? -(solver->common.initmem)
        : solver->common.initmem_amd * solver->symbolic->Lnz[block] +
            (double)nk;
    old_lusize = kls_initial_block_lusize(solver, k1, nk, lsize);
    if (old_lusize == 0u) {
      free(scratch);
      return 0;
    }
  }
  Unit *new_lu =
    (Unit *)TRILINOS_KLU_malloc(old_lusize, sizeof(Unit), &solver->common);
  if (new_lu == NULL) {
    free(scratch);
    return 0;
  }

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
  const int preserve_suffix = restart_end < nk;
  const int use_restart_mask = restart_tail_mask != NULL;
  int promoted_gap_columns = 0;
  unsigned char *active_tail_mask = NULL;
  if (use_restart_mask) {
    active_tail_mask = (unsigned char *)malloc((size_t)nk *
                                               sizeof(*active_tail_mask));
    if (active_tail_mask == NULL) {
      goto fail;
    }
    memcpy(active_tail_mask, restart_tail_mask,
           (size_t)nk * sizeof(*active_tail_mask));
  }
  const int use_mapped_tail =
    kls_tail_mapped_block_is_valid(solver, k1, nk, psinv);
  UF_long mapped_columns = 0;
  UF_long skipped_columns = 0;
  double skipped_work = 0.0;
  for (UF_long kk = local_reject; kk < restart_end; ++kk) {
    const Int k = (Int)kk;
    int lock_gap_pivot =
      use_restart_mask && active_tail_mask != NULL &&
      !active_tail_mask[kk];
    if (lock_gap_pivot && pblock[kk] != old_pblock[kk]) {
      solver->common.status = TRILINOS_KLU_INVALID;
      goto fail;
    }
    if (lock_gap_pivot) {
      if (kls_copy_finished_gap_lu_column(
            &new_lu, &lusize, old_lu, old_lusize, lip, llen, uip, ulen,
            udiag, old_pblock, active_tail_mask, live_pinv, nk, kk, lip,
            llen, uip, ulen, udiag, &lup, &lnz, &unz,
            &solver->common)) {
        const UF_long pivrow = old_pblock[kk];
        pblock[kk] = pivrow;
        live_pinv[pivrow] = kk;
        kls_tail_prune(lpend, live_pinv, (Int)kk, (Int)pivrow, new_lu, uip,
                       lip, ulen, llen);
        skipped_columns++;
        double column_work = 0.0;
        if (kls_fast_reject_column_work(solver, block, k1, nk, kk,
                                        &column_work)) {
          skipped_work += column_work;
        }
        continue;
      }
      if (solver->common.status < TRILINOS_KLU_OK) {
        goto fail;
      }
      if (kls_finished_gap_lu_column_depends_on_tail(
            old_lu, old_lusize, lip, llen, uip, ulen, old_pblock,
            active_tail_mask, live_pinv, nk, kk)) {
        active_tail_mask[kk] = 1u;
        promoted_gap_columns++;
        lock_gap_pivot = 0;
      } else {
        /* Non-dependency copy failures are structural or mapping failures, so
           keep them on the wider fallback instead of speculating here. */
        goto fail;
      }
    }
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
      use_mapped_tail
        ? kls_tail_lsolve_symbolic_mapped(
            solver, (Int)nk, k, live_pinv, stack, flag, lpend, ap_pos,
            new_lu, (Int)lup, llen, lip, (Int)k1)
        : kls_tail_lsolve_symbolic(
            (Int)nk, k, solver->col_ptr, solver->row_idx,
            solver->symbolic->Q, live_pinv, stack, flag, lpend, ap_pos,
            new_lu, (Int)lup, llen, lip, (Int)k1, (Int *)psinv);
    if (top < 0) {
      solver->common.status = TRILINOS_KLU_INVALID;
      goto fail;
    }
    /* Trial off-block row/value entries are discarded: accepted repairs
       rebuild Offi/Offx from the final Pinv before publishing stats. */
    const int constructed =
      use_mapped_tail
        ? kls_tail_construct_column_mapped(
            solver, numeric_values, k, (Entry *)solver->numeric->Xwork,
            (Int)nk, (Int)k1, scaled ? solver->numeric->Rs : NULL,
            scaled ? (Int)solver->common.scale : 0)
        : kls_tail_construct_column(
            k, solver->col_ptr, solver->row_idx, numeric_values,
            solver->symbolic->Q, (Entry *)solver->numeric->Xwork, (Int)nk,
            (Int)k1, (Int *)psinv, scaled ? solver->numeric->Rs : NULL,
            scaled ? (Int)solver->common.scale : 0,
            NULL, NULL, NULL, 0);
    if (!constructed) {
      solver->common.status = TRILINOS_KLU_INVALID;
      goto fail;
    }
    if (use_mapped_tail) {
      mapped_columns++;
    }
    kls_tail_lsolve_numeric(live_pinv, new_lu, stack, lip, top, (Int)nk,
                            llen, (Entry *)solver->numeric->Xwork);

    const Int diagrow = pblock[k];
    Int pivrow = TRILINOS_KLU_EMPTY;
    Entry pivot = 0.0;
    double abs_pivot = 0.0;
    const double pivot_tol = lock_gap_pivot ? 0.0 : tol;
    if (!kls_tail_lpivot(diagrow, &pivrow, &pivot, &abs_pivot, pivot_tol,
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
    if (lock_gap_pivot && (UF_long)pivrow != old_pblock[kk]) {
      solver->common.status = TRILINOS_KLU_INVALID;
      goto fail;
    }
    if (preserve_suffix) {
      const UF_long old_pivot_pos = final_pinv[pivrow];
      if (old_pivot_pos >= restart_end) {
        solver->common.status = TRILINOS_KLU_INVALID;
        goto fail;
      }
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
        if (preserve_suffix && (UF_long)kbar >= restart_end) {
          solver->common.status = TRILINOS_KLU_INVALID;
          goto fail;
        }
        pblock[kbar] = diagrow;
        live_pinv[diagrow] = FLIP(kbar);
      }
    }
    pblock[k] = pivrow;
    live_pinv[pivrow] = k;
    if (lock_gap_pivot && pblock[kk] != old_pblock[kk]) {
      solver->common.status = TRILINOS_KLU_INVALID;
      goto fail;
    }

    kls_tail_prune(lpend, live_pinv, k, pivrow, new_lu, uip, lip, ulen,
                   llen);
    lnz += llen[k] + 1u;
    unz += ulen[k] + 1u;
  }

  if (preserve_suffix) {
    for (UF_long k = restart_end; k < nk; ++k) {
      pblock[k] = old_pblock[k];
    }
  }
  for (UF_long k = 0; k < nk; ++k) {
    live_pinv[k] = KLS_KLU_EMPTY;
  }
  for (UF_long k = 0; k < nk; ++k) {
    const UF_long local_row = pblock[k];
    if (local_row >= nk || live_pinv[local_row] != KLS_KLU_EMPTY) {
      solver->common.status = TRILINOS_KLU_INVALID;
      goto fail;
    }
    live_pinv[local_row] = k;
  }
  if (preserve_suffix &&
      !kls_copy_preserved_lu_columns(
        &new_lu, &lusize, old_lu, old_lusize, lip, llen, uip, ulen, udiag,
        old_pblock, live_pinv, nk, restart_end, nk, lip, llen, uip, ulen,
        udiag, &lup, &lnz, &unz, &solver->common)) {
    solver->common.status = TRILINOS_KLU_INVALID;
    goto fail;
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
  free(scratch);
  free(active_tail_mask);
  *new_lu_out = new_lu;
  *new_size_out = lusize;
  *lnz_block_out = lnz;
  *unz_block_out = unz;
  *skipped_columns_out = skipped_columns;
  *skipped_work_out = skipped_work;
  *exact_mask_used_out = use_restart_mask && promoted_gap_columns == 0;
  solver->kls_tail_last_mapped_columns += mapped_columns;
  solver->kls_tail_mapped_column_count += mapped_columns;
  return 1;

fail:
  for (UF_long k = 0; k < nk; ++k) {
    ((Entry *)solver->numeric->Xwork)[k] = 0.0;
  }
  (void)TRILINOS_KLU_free(new_lu, lusize, sizeof(Unit), &solver->common);
  free(scratch);
  free(active_tail_mask);
  KLS_SET_OPTIONAL_OUTPUT(new_lu_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(new_size_out, 0u);
  KLS_SET_OPTIONAL_OUTPUT(skipped_columns_out, 0);
  if (skipped_work_out != NULL) {
    *skipped_work_out = 0.0;
  }
  KLS_SET_OPTIONAL_OUTPUT(exact_mask_used_out, 0);
  return 0;
}

static void kls_record_fast_repaired_block_stats(kls_solver *solver,
                                                 UF_long block,
                                                 UF_long k1,
                                                 UF_long nk,
                                                 UF_long rejected_pivot,
                                                 UF_long repaired_tail_begin,
                                                 UF_long repaired_tail_end,
                                                 UF_long skipped_columns,
                                                 double skipped_work,
                                                 int exact_mask_used,
                                                 int etree_mask_used,
                                                 const UF_long *old_pblock,
                                                 const UF_long *pblock) {
  if (solver == NULL || old_pblock == NULL || pblock == NULL ||
      rejected_pivot < k1 || rejected_pivot >= k1 + nk) {
    return;
  }

  const UF_long local_reject = rejected_pivot - k1;
  if (repaired_tail_begin > local_reject || repaired_tail_begin >= nk) {
    repaired_tail_begin = local_reject;
  }
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
  if (repaired_tail_end <= repaired_tail_begin || repaired_tail_end > nk) {
    repaired_tail_end = nk;
  }
  if (skipped_columns > nk) {
    skipped_columns = 0;
    skipped_work = 0.0;
  }
  if (skipped_work < 0.0) {
    skipped_work = 0.0;
  }
  solver->stats.fast_repaired_tail_restart_exact_mask = 0;
  solver->stats.fast_repaired_tail_restart_etree_mask = 0;

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
  const UF_long repaired_tail_global_begin = k1 + repaired_tail_begin;
  solver->stats.fast_repaired_tail_restart_ready =
    serial_tail_state &&
    (repaired_tail_begin > 0u || repaired_tail_end < nk ||
     skipped_columns > 0u) &&
    (first_changed == KLS_KLU_EMPTY ||
     first_changed >= repaired_tail_global_begin) &&
    pivoting_tail_ready;
  if (solver->stats.fast_repaired_tail_restart_ready &&
      !kls_repaired_tail_prefix_state_valid(solver, block, k1, nk,
                                            repaired_tail_begin, pblock)) {
    solver->stats.fast_repaired_tail_restart_ready = 0;
  }
  double block_work = 0.0;
  if (kls_fast_reject_block_work(solver, block, k1, nk, &block_work)) {
    solver->stats.fast_repaired_block_work = block_work;
    if (solver->stats.fast_repaired_tail_restart_ready) {
      UF_long tail_columns = 0;
      double tail_work = 0.0;
      if (kls_fast_reject_range_work(
            solver, block, k1, nk, repaired_tail_begin, repaired_tail_end,
            &tail_columns, &tail_work)) {
        const UF_long actual_tail_columns =
          skipped_columns < tail_columns ? tail_columns - skipped_columns
                                         : 0u;
        const double actual_tail_work =
          skipped_work < tail_work ? tail_work - skipped_work : 0.0;
        solver->stats.fast_repaired_tail_restart_skipped_columns =
          (int64_t)(tail_columns - actual_tail_columns);
        solver->stats.fast_repaired_tail_restart_skipped_work =
          tail_work - actual_tail_work;
        solver->stats.fast_repaired_tail_restart_columns =
          (int64_t)actual_tail_columns;
        solver->stats.fast_repaired_tail_restart_work = actual_tail_work;
        solver->stats.fast_repaired_tail_restart_saved_work =
          block_work > actual_tail_work ? block_work - actual_tail_work : 0.0;
        const int64_t pivoting_tail_columns =
          solver->stats.fast_rejected_pivoting_tail_columns;
        const double pivoting_tail_work =
          solver->stats.fast_rejected_pivoting_tail_work;
        solver->stats.fast_repaired_tail_restart_overcompute_columns =
          pivoting_tail_columns > 0 &&
              actual_tail_columns > (UF_long)pivoting_tail_columns
            ? (int64_t)(actual_tail_columns - (UF_long)pivoting_tail_columns)
            : 0;
        solver->stats.fast_repaired_tail_restart_overcompute_work =
          pivoting_tail_work > 0.0 && actual_tail_work > pivoting_tail_work
            ? actual_tail_work - pivoting_tail_work
            : 0.0;
        solver->stats.fast_repaired_tail_restart_exact_mask =
          exact_mask_used &&
          pivoting_tail_columns > 0 &&
          actual_tail_columns == (UF_long)pivoting_tail_columns &&
          solver->stats.fast_repaired_tail_restart_overcompute_columns == 0;
        solver->stats.fast_repaired_tail_restart_etree_mask =
          etree_mask_used &&
          solver->stats.fast_repaired_tail_restart_exact_mask;
      } else {
        solver->stats.fast_repaired_tail_restart_ready = 0;
      }
    }
  }
}

static int kls_try_fast_reject_kls_first_rebuild(kls_solver *solver,
                                                 double *numeric_values) {
  if (solver == NULL || numeric_values == NULL) {
    return 0;
  }
  if (!kls_should_try_first_factor_recovery(solver)) {
    return 0;
  }
  if (!kls_try_rebuild_current_numeric_with_kls_first_mode(
        solver, numeric_values, NULL, 1) ||
      solver->common.status < TRILINOS_KLU_OK ||
      solver->common.status == TRILINOS_KLU_SINGULAR) {
    return 0;
  }
  solver->fast_kls_rebuild_restarts++;
  solver->stats.fast_kls_rebuild_restarts =
    solver->fast_kls_rebuild_restarts;
  return 1;
}

typedef enum {
  KLS_FAST_REJECT_REPAIR_DEFAULT = 0,
  KLS_FAST_REJECT_REPAIR_ETREE = 1,
  KLS_FAST_REJECT_REPAIR_SERIAL = 2,
  KLS_FAST_REJECT_REPAIR_BLOCK = 3,
  KLS_FAST_REJECT_REPAIR_FULL = 4
} kls_fast_reject_repair_mode;

/* Diagnostic control consulted only after a checked factor has rejected a
   pivot.  It isolates the bounded ETree row pipeline from the serial exact
   tail, whole rejected-block KLU, and whole-numeric KLU fallbacks without
   putting an environment lookup on the successful fixed-pattern path. */
static kls_fast_reject_repair_mode kls_fast_reject_repair_mode_from_env(void) {
  const char *value = getenv("KLS_FAST_REJECT_REPAIR_MODE");
  if (value == NULL || value[0] == '\0') {
    return KLS_FAST_REJECT_REPAIR_DEFAULT;
  }
  if (strcmp(value, "etree") == 0) {
    return KLS_FAST_REJECT_REPAIR_ETREE;
  }
  if (strcmp(value, "serial") == 0) {
    return KLS_FAST_REJECT_REPAIR_SERIAL;
  }
  if (strcmp(value, "block") == 0) {
    return KLS_FAST_REJECT_REPAIR_BLOCK;
  }
  if (strcmp(value, "full") == 0) {
    return KLS_FAST_REJECT_REPAIR_FULL;
  }
  return KLS_FAST_REJECT_REPAIR_DEFAULT;
}

static int kls_pivot_restart_rejected_block(kls_solver *solver,
                                            double *numeric_values,
                                            UF_long rejected_pivot,
                                            kls_fast_reject_repair_mode
                                              repair_mode) {
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
  const UF_long old_nrealloc = solver->common.nrealloc;
  const int scaled_all_refresh =
    scaled &&
    solver->fast_reject_refresh_state == KLS_FAST_REJECT_REFRESH_ALL;
  /* Full KLU refactor leaves Rs in pivot order; repair kernels need input-row
     order and the accepted path permutes it back after Pnum is installed. */
  if (scaled_all_refresh &&
      !trilinos_klu_l_scale((UF_long)solver->common.scale, solver->n,
                            solver->col_ptr, solver->row_idx,
                            numeric_values, solver->numeric->Rs, NULL,
                            &solver->common)) {
    free(old_pblock);
    solver->common.status = old_status;
    solver->common.numerical_rank = old_numerical_rank;
    solver->common.singular_col = old_singular_col;
    solver->common.noffdiag = old_noffdiag;
    (void)kls_rebuild_numeric_pinv(solver);
    return 0;
  }
  solver->common.status = TRILINOS_KLU_OK;
  solver->common.numerical_rank = KLS_KLU_EMPTY;
  solver->common.singular_col = KLS_KLU_EMPTY;

  size_t new_size = 0u;
  UF_long repaired_tail_end = nk;
  UF_long repaired_tail_begin = local_reject;
  UF_long repaired_tail_skipped_columns = 0;
  double repaired_tail_skipped_work = 0.0;
  int repaired_tail_exact_mask = 0;
  int repaired_tail_etree_mask = 0;
  UF_long tail_plan_begin = 0;
  UF_long tail_plan_end = 0;
  int tail_plan_suffix_exact = 0;
  const int have_retained_tail_plan =
    kls_pivot_tail_plan_local_range(
      solver, k1, nk, local_reject, &tail_plan_begin, &tail_plan_end,
      NULL, NULL, &tail_plan_suffix_exact);
  UF_long preferred_restart_end =
    kls_preferred_pivot_tail_restart_end(solver, k1, nk, local_reject);
  const int can_try_masked_tail =
    have_retained_tail_plan && !tail_plan_suffix_exact;
  if (preferred_restart_end == nk && !can_try_masked_tail) {
    preferred_restart_end =
      kls_prepare_root_pivot_tail_independent_refresh(
        solver, numeric_values, block, k1, nk, local_reject);
  }
  int kls_block_restart_used = 0;
  int kls_row_first_block_rebuild_tried = 0;
  const int allow_partial_repair =
    repair_mode != KLS_FAST_REJECT_REPAIR_BLOCK;
  const int allow_row_pipeline =
    repair_mode != KLS_FAST_REJECT_REPAIR_SERIAL;
  const int allow_kls_first_rebuild =
    repair_mode != KLS_FAST_REJECT_REPAIR_SERIAL;
  const int prefer_row_pipeline_block_repair =
    allow_partial_repair && allow_row_pipeline &&
    solver->options.threads > 1 && nk > 1u;
  if (getenv("KLS_TRACE_FAST_REJECT_REPAIR") != NULL) {
    fprintf(stderr,
            "KLS rejected-block repair dispatch: mode=%d threads=%d nk=%ld"
            " row_pipeline=%d masked_tail=%d restart_end=%ld\n",
            (int)repair_mode, solver->options.threads, (long)nk,
            prefer_row_pipeline_block_repair, can_try_masked_tail,
            (long)preferred_restart_end);
  }
  if (prefer_row_pipeline_block_repair) {
    kls_row_first_block_rebuild_tried = 1;
    if (kls_try_row_first_rebuild_rejected_block(
          solver, numeric_values, block, k1, nk, rejected_pivot,
          old_pblock, old_lnz_block, old_unz_block, pblock)) {
      free(old_pblock);
      return 1;
    }
  }
  if (allow_partial_repair &&
      (preferred_restart_end < nk || can_try_masked_tail)) {
    unsigned char *tail_mask = NULL;
    UF_long *saved_ptrs = NULL;
    Entry *saved_udiag = NULL;
    UF_long restart_begin = local_reject;
    UF_long restart_end = preferred_restart_end;
    if (can_try_masked_tail) {
      restart_begin = tail_plan_begin;
      restart_end = tail_plan_end;
    }
    const size_t nk_size = (size_t)nk;
    if ((UF_long)nk_size == nk) {
      tail_mask = (unsigned char *)malloc(nk_size * sizeof(*tail_mask));
    }
    if ((UF_long)nk_size == nk &&
        nk_size <= SIZE_MAX / (4u * sizeof(*saved_ptrs))) {
      saved_ptrs =
        (UF_long *)malloc(4u * nk_size * sizeof(*saved_ptrs));
      saved_udiag = (Entry *)malloc(nk_size * sizeof(*saved_udiag));
    }
    const int main_mask_ok =
      tail_mask != NULL &&
      kls_build_pivot_tail_restart_mask_range(
        solver, k1, nk, local_reject, restart_begin, restart_end,
        tail_mask);
    const int main_mask_etree_ok =
      main_mask_ok &&
      kls_pivot_tail_restart_mask_matches_etree_forest(
        solver, block, k1, nk, local_reject, restart_begin, restart_end,
        tail_mask);
    if (main_mask_etree_ok && saved_ptrs != NULL && saved_udiag != NULL &&
        (solver->fast_reject_refresh_state !=
           KLS_FAST_REJECT_REFRESH_PREFIX ||
         kls_refresh_pivot_tail_preserved_block_columns(
           solver, numeric_values, block, k1, nk, tail_mask))) {
      UF_long *saved_lip = saved_ptrs;
      UF_long *saved_llen = saved_ptrs + nk;
      UF_long *saved_uip = saved_ptrs + 2u * nk;
      UF_long *saved_ulen = saved_ptrs + 3u * nk;
      memcpy(saved_lip, solver->numeric->Lip + k1,
             (size_t)nk * sizeof(*saved_lip));
      memcpy(saved_llen, solver->numeric->Llen + k1,
             (size_t)nk * sizeof(*saved_llen));
      memcpy(saved_uip, solver->numeric->Uip + k1,
             (size_t)nk * sizeof(*saved_uip));
      memcpy(saved_ulen, solver->numeric->Ulen + k1,
             (size_t)nk * sizeof(*saved_ulen));
      memcpy(saved_udiag, ((Entry *)solver->numeric->Udiag) + k1,
             (size_t)nk * sizeof(*saved_udiag));
      kls_block_restart_used =
        kls_try_pivot_tail_restart_rejected_block(
          solver, numeric_values, block, k1, nk, restart_begin,
          restart_end, tail_mask, old_pblock, psinv, &new_lu,
          &new_size, &lnz_block, &unz_block, pblock,
          &repaired_tail_skipped_columns, &repaired_tail_skipped_work,
          &repaired_tail_exact_mask);
      if (kls_block_restart_used) {
        repaired_tail_begin = restart_begin;
        repaired_tail_end = restart_end;
        repaired_tail_etree_mask = 1;
      }
      if (!kls_block_restart_used) {
        repaired_tail_skipped_columns = 0;
        repaired_tail_skipped_work = 0.0;
        repaired_tail_exact_mask = 0;
        repaired_tail_etree_mask = 0;
        solver->common.nrealloc = old_nrealloc;
        memcpy(solver->numeric->Lip + k1, saved_lip,
               (size_t)nk * sizeof(*saved_lip));
        memcpy(solver->numeric->Llen + k1, saved_llen,
               (size_t)nk * sizeof(*saved_llen));
        memcpy(solver->numeric->Uip + k1, saved_uip,
               (size_t)nk * sizeof(*saved_uip));
        memcpy(solver->numeric->Ulen + k1, saved_ulen,
               (size_t)nk * sizeof(*saved_ulen));
        memcpy(((Entry *)solver->numeric->Udiag) + k1, saved_udiag,
               (size_t)nk * sizeof(*saved_udiag));
      }
    }
    free(tail_mask);
    free(saved_ptrs);
    free(saved_udiag);
  }
  if (allow_partial_repair && !kls_block_restart_used &&
      allow_row_pipeline && can_try_masked_tail &&
      solver->stats.fast_rejected_pivoting_tail_topological) {
    kls_row_first_block_rebuild_tried = 1;
    if (kls_try_row_first_rebuild_rejected_block(
          solver, numeric_values, block, k1, nk, rejected_pivot,
          old_pblock, old_lnz_block, old_unz_block, pblock)) {
      free(old_pblock);
      return 1;
    }
  }
  if (allow_partial_repair && !kls_block_restart_used &&
      allow_kls_first_rebuild && can_try_masked_tail &&
      solver->stats.fast_rejected_pivoting_tail_topological &&
      kls_try_fast_reject_kls_first_rebuild(solver, numeric_values)) {
    free(old_pblock);
    return 1;
  }
  if (allow_partial_repair && !kls_block_restart_used) {
    UF_long *saved_ptrs = NULL;
    Entry *saved_udiag = NULL;
    const size_t nk_size = (size_t)nk;
    if ((UF_long)nk_size == nk &&
        nk_size <= SIZE_MAX / (4u * sizeof(*saved_ptrs))) {
      saved_ptrs =
        (UF_long *)malloc(4u * nk_size * sizeof(*saved_ptrs));
      saved_udiag = (Entry *)malloc(nk_size * sizeof(*saved_udiag));
    }
    if (saved_ptrs != NULL && saved_udiag != NULL) {
      UF_long *saved_lip = saved_ptrs;
      UF_long *saved_llen = saved_ptrs + nk;
      UF_long *saved_uip = saved_ptrs + 2u * nk;
      UF_long *saved_ulen = saved_ptrs + 3u * nk;
      memcpy(saved_lip, solver->numeric->Lip + k1,
             (size_t)nk * sizeof(*saved_lip));
      memcpy(saved_llen, solver->numeric->Llen + k1,
             (size_t)nk * sizeof(*saved_llen));
      memcpy(saved_uip, solver->numeric->Uip + k1,
             (size_t)nk * sizeof(*saved_uip));
      memcpy(saved_ulen, solver->numeric->Ulen + k1,
             (size_t)nk * sizeof(*saved_ulen));
      memcpy(saved_udiag, ((Entry *)solver->numeric->Udiag) + k1,
             (size_t)nk * sizeof(*saved_udiag));
    }
    solver->common.status = TRILINOS_KLU_OK;
    solver->common.numerical_rank = KLS_KLU_EMPTY;
    solver->common.singular_col = KLS_KLU_EMPTY;
    solver->common.noffdiag = old_noffdiag;
    kls_block_restart_used =
      kls_try_pivot_tail_restart_rejected_block(
        solver, numeric_values, block, k1, nk, local_reject, nk,
        NULL, old_pblock, psinv, &new_lu, &new_size, &lnz_block,
        &unz_block, pblock, &repaired_tail_skipped_columns,
        &repaired_tail_skipped_work, &repaired_tail_exact_mask);
    if (kls_block_restart_used) {
      repaired_tail_begin = local_reject;
    }
    repaired_tail_end = nk;
    if (!kls_block_restart_used && saved_ptrs != NULL &&
        saved_udiag != NULL) {
      UF_long *saved_lip = saved_ptrs;
      UF_long *saved_llen = saved_ptrs + nk;
      UF_long *saved_uip = saved_ptrs + 2u * nk;
      UF_long *saved_ulen = saved_ptrs + 3u * nk;
      solver->common.status = TRILINOS_KLU_OK;
      solver->common.numerical_rank = KLS_KLU_EMPTY;
      solver->common.singular_col = KLS_KLU_EMPTY;
      solver->common.noffdiag = old_noffdiag;
      solver->common.nrealloc = old_nrealloc;
      memcpy(solver->numeric->Lip + k1, saved_lip,
             (size_t)nk * sizeof(*saved_lip));
      memcpy(solver->numeric->Llen + k1, saved_llen,
             (size_t)nk * sizeof(*saved_llen));
      memcpy(solver->numeric->Uip + k1, saved_uip,
             (size_t)nk * sizeof(*saved_uip));
      memcpy(solver->numeric->Ulen + k1, saved_ulen,
             (size_t)nk * sizeof(*saved_ulen));
      memcpy(((Entry *)solver->numeric->Udiag) + k1, saved_udiag,
             (size_t)nk * sizeof(*saved_udiag));
      if (allow_row_pipeline && !kls_row_first_block_rebuild_tried) {
        kls_row_first_block_rebuild_tried = 1;
        if (kls_try_row_first_rebuild_rejected_block(
              solver, numeric_values, block, k1, nk, rejected_pivot,
              old_pblock, old_lnz_block, old_unz_block, pblock)) {
          free(saved_ptrs);
          free(saved_udiag);
          free(old_pblock);
          return 1;
        }
      }
      if (allow_kls_first_rebuild &&
          kls_try_fast_reject_kls_first_rebuild(solver, numeric_values)) {
        free(saved_ptrs);
        free(saved_udiag);
        free(old_pblock);
        return 1;
      }
    }
    free(saved_ptrs);
    free(saved_udiag);
  }
  const int tail_restart_used =
    kls_block_restart_used &&
    (repaired_tail_begin > 0u || repaired_tail_end < nk ||
     repaired_tail_skipped_columns > 0u);
  if (!kls_block_restart_used) {
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
    if (scaled_all_refresh) {
      (void)kls_parallel_refactor_permute_scale(solver);
    }
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
      if (scaled_all_refresh) {
        (void)kls_parallel_refactor_permute_scale(solver);
      }
      (void)kls_rebuild_numeric_pinv(solver);
      return 0;
    }
    solver->numeric->Pnum[k1 + k] = solver->symbolic->P[k1 + local_row];
  }
  kls_record_fast_repaired_block_stats(solver, block, k1, nk, rejected_pivot,
                                       repaired_tail_begin, repaired_tail_end,
                                       repaired_tail_skipped_columns,
                                       repaired_tail_skipped_work,
                                       repaired_tail_exact_mask,
                                       repaired_tail_etree_mask,
                                       old_pblock, pblock);
  free(old_pblock);

  if (!kls_rebuild_numeric_pinv(solver)) {
    solver->common.status = old_status;
    solver->common.numerical_rank = old_numerical_rank;
    solver->common.singular_col = old_singular_col;
    solver->common.noffdiag = old_noffdiag;
    if (scaled_all_refresh) {
      (void)kls_parallel_refactor_permute_scale(solver);
    }
    return 0;
  }

  int offdiag_refreshed = 0;
  solver->fast_repaired_last_offdiag_suffix_refresh = 0;
  if (tail_restart_used) {
    offdiag_refreshed =
      kls_recompute_offdiag_suffix_from_pinv(solver, numeric_values,
                                             k1 + repaired_tail_begin);
    if (offdiag_refreshed) {
      solver->fast_repaired_last_offdiag_suffix_refresh = 1;
      solver->fast_repaired_offdiag_suffix_refresh_count++;
    }
  }
  if (!offdiag_refreshed) {
    offdiag_refreshed = kls_recompute_offdiag_from_pinv(solver, numeric_values);
    if (offdiag_refreshed) {
      solver->fast_repaired_offdiag_full_refresh_count++;
    }
  }
  if (!offdiag_refreshed) {
    solver->common.status = old_status;
    solver->common.numerical_rank = old_numerical_rank;
    solver->common.singular_col = old_singular_col;
    solver->common.noffdiag = old_noffdiag;
    if (scaled_all_refresh) {
      (void)kls_parallel_refactor_permute_scale(solver);
    }
    return 0;
  }
  if (scaled && !kls_parallel_refactor_permute_scale(solver)) {
    solver->common.status = old_status;
    solver->common.numerical_rank = old_numerical_rank;
    solver->common.singular_col = old_singular_col;
    solver->common.noffdiag = old_noffdiag;
    if (scaled_all_refresh) {
      (void)kls_parallel_refactor_permute_scale(solver);
    }
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
  kls_invalidate_i32_solve(solver);
  solver->fast_block_restarts++;
  if (kls_block_restart_used) {
    solver->fast_kls_block_restarts++;
  }
  if (tail_restart_used) {
    solver->fast_tail_restarts++;
  }
  solver->stats.fast_block_restarts = solver->fast_block_restarts;
  solver->stats.fast_kls_block_restarts = solver->fast_kls_block_restarts;
  solver->stats.fast_kls_rebuild_restarts =
    solver->fast_kls_rebuild_restarts;
  solver->stats.fast_tail_restarts = solver->fast_tail_restarts;
  solver->stats.fast_repaired_last_offdiag_suffix_refresh =
    solver->fast_repaired_last_offdiag_suffix_refresh;
  solver->stats.fast_repaired_offdiag_suffix_refresh_count =
    (int64_t)solver->fast_repaired_offdiag_suffix_refresh_count;
  solver->stats.fast_repaired_offdiag_full_refresh_count =
    (int64_t)solver->fast_repaired_offdiag_full_refresh_count;
  free_refactor_map(solver);
  free_refactor_schedule(solver);
  free_row_refactor_pattern_preserve_diagnostics(solver);
  free_refactor_lu_pointer_cache(solver);
  free_snode_panels(solver);
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

static int kls_refactor_map_index32_env_enabled(void) {
  const char *value = getenv("KLS_ENABLE_REFACTOR_MAP_INDEX32");
  return value == NULL || value[0] == '\0' ||
         !(value[0] == '0' && value[1] == '\0');
}

static int kls_build_refactor_map_index32(kls_solver *solver) {
  if (!kls_refactor_map_index32_env_enabled()) {
    return 1;
  }
  if (solver != NULL && !solver->compact_map32_probe_active &&
      getenv("KLS_ENABLE_SETTLED_COMPACT_MAP32_REFACTOR") == NULL &&
      (kls_moderate_work_single_block_lean_policy_enabled(solver) ||
       ((solver->lean_choice > 0 || solver->lean_probe_arm > 0) &&
        kls_moderate_work_fragmented_btf_lean_policy_enabled(solver)))) {
    /* A selected (or currently probed) lean engine consumes the row-ordered
       64-bit map directly; neither its fused updates nor the native column
       solve reads this copy.  The moderate single-block compact walk showed
       no recurring benefit over its native stream, so it also keeps avoiding
       the duplicate.  For fragmented BTF, however, structural eligibility
       alone is insufficient: once the timed lean trial rejects that engine,
       the ordinary mapped incumbent benefits from this compact stream on
       every later update. */
    return 1;
  }
  if (solver == NULL || solver->refactor_row_idx == NULL ||
      solver->refactor_input_pos == NULL || solver->n <= 0 ||
      solver->nnz == 0u) {
    return 0;
  }
  if (solver->n > (UF_long)INT32_MAX ||
      solver->nnz > (UF_long)INT32_MAX ||
      solver->nnz > (UF_long)(SIZE_MAX / sizeof(int32_t))) {
    return 1;
  }
  if (solver->refactor_row_idx32 != NULL &&
      solver->refactor_input_pos32 != NULL &&
      solver->refactor_input_oldrow32 != NULL) {
    return 1;
  }

  int32_t *row_idx32 =
    (int32_t *)malloc((size_t)solver->nnz * sizeof(*row_idx32));
  int32_t *input_pos32 =
    (int32_t *)malloc((size_t)solver->nnz * sizeof(*input_pos32));
  int32_t *input_oldrow32 =
    (int32_t *)malloc((size_t)solver->nnz * sizeof(*input_oldrow32));
  if (row_idx32 == NULL || input_pos32 == NULL || input_oldrow32 == NULL) {
    free(row_idx32);
    free(input_pos32);
    free(input_oldrow32);
    return 1;
  }
  for (UF_long p = 0; p < solver->nnz; ++p) {
    const UF_long row = solver->refactor_row_idx[p];
    const UF_long input_pos = solver->refactor_input_pos[p];
    if (row > (UF_long)INT32_MAX || input_pos > (UF_long)INT32_MAX ||
        input_pos >= solver->nnz || solver->row_idx == NULL ||
        solver->row_idx[input_pos] > (UF_long)INT32_MAX) {
      free(row_idx32);
      free(input_pos32);
      free(input_oldrow32);
      return 0;
    }
    row_idx32[p] = (int32_t)row;
    input_pos32[p] = (int32_t)input_pos;
    input_oldrow32[p] = (int32_t)solver->row_idx[input_pos];
  }

  solver->refactor_row_idx32 = row_idx32;
  solver->refactor_input_pos32 = input_pos32;
  solver->refactor_input_oldrow32 = input_oldrow32;
  solver->refactor_map_indices32_count = solver->nnz;
  solver->refactor_map_index32_enabled = 1;
  return 1;
}

/* Transpose only the input-position metadata, not the matrix values.  Rows
   are disjoint scale reductions and positions remain in original CSC order,
   preserving KLU's summation order exactly while exposing all rows to the
   already-resident numeric worker crew. */
static int kls_build_refactor_scale_rows(kls_solver *solver) {
  if (solver == NULL || solver->row_idx == NULL || solver->n == 0u ||
      solver->nnz == 0u || solver->nnz > (UF_long)INT32_MAX ||
      solver->n > (UF_long)(SIZE_MAX / sizeof(UF_long)) - 1u ||
      solver->nnz > (UF_long)(SIZE_MAX / sizeof(int32_t))) {
    return 0;
  }
  if (solver->refactor_scale_row_ptr != NULL &&
      solver->refactor_scale_input_pos32 != NULL) {
    return 1;
  }

  UF_long *row_ptr = (UF_long *)calloc(
    (size_t)solver->n + 1u, sizeof(*row_ptr));
  UF_long *cursor = (UF_long *)malloc(
    (size_t)solver->n * sizeof(*cursor));
  int32_t *input_pos = (int32_t *)malloc(
    (size_t)solver->nnz * sizeof(*input_pos));
  if (row_ptr == NULL || cursor == NULL || input_pos == NULL) {
    free(row_ptr);
    free(cursor);
    free(input_pos);
    return 0;
  }
  for (UF_long p = 0u; p < solver->nnz; ++p) {
    const UF_long row = solver->row_idx[p];
    if (row >= solver->n) {
      free(row_ptr);
      free(cursor);
      free(input_pos);
      return 0;
    }
    row_ptr[row + 1u]++;
  }
  for (UF_long row = 0u; row < solver->n; ++row) {
    row_ptr[row + 1u] += row_ptr[row];
    cursor[row] = row_ptr[row];
  }
  for (UF_long p = 0u; p < solver->nnz; ++p) {
    const UF_long row = solver->row_idx[p];
    input_pos[cursor[row]++] = (int32_t)p;
  }
  free(cursor);
  solver->refactor_scale_row_ptr = row_ptr;
  solver->refactor_scale_input_pos32 = input_pos;
  return 1;
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
    return kls_build_refactor_map_index32(solver);
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
    return kls_build_refactor_map_index32(solver);
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
  return kls_build_refactor_map_index32(solver);
}

static int kls_build_refactor_user_input_pos32(kls_solver *solver) {
  if (solver == NULL || solver->refactor_input_pos == NULL ||
      solver->nnz == 0u ||
      solver->nnz > (UF_long)INT32_MAX ||
      solver->nnz > (UF_long)(SIZE_MAX / sizeof(int32_t))) {
    return 0;
  }
  /* A NULL input_to_csc is the common fast case: the public values already
     use the retained CSC order, so the inverse map is the identity. */
  const int input_ok = solver->input_to_csc == NULL ||
    kls_build_prepared_value_input_pos(solver);
  const int scale_ok =
    (solver->row_scale == NULL && solver->col_scale == NULL) ||
    kls_build_prepared_value_scale(solver);
  if (!input_ok || !scale_ok) {
    return 0;
  }
  if (solver->refactor_input_user_pos32 != NULL) {
    return 1;
  }
  int32_t *user_pos = (int32_t *)malloc(
    (size_t)solver->nnz * sizeof(*user_pos));
  if (user_pos == NULL) {
    return 0;
  }
  for (UF_long p = 0u; p < solver->nnz; ++p) {
    const UF_long internal = solver->refactor_input_pos[p];
    if (internal >= solver->nnz) {
      free(user_pos);
      return 0;
    }
    const UF_long user = solver->prepared_value_input_pos != NULL
      ? solver->prepared_value_input_pos[internal] : internal;
    if (user >= solver->nnz || user > (UF_long)INT32_MAX) {
      free(user_pos);
      return 0;
    }
    user_pos[p] = (int32_t)user;
  }
  solver->refactor_input_user_pos32 = user_pos;
  return 1;
}

static int kls_pts_direct_user_values_enabled(const kls_solver *solver) {
  return kls_pts_refactor_ready(solver) ||
    getenv("KLS_ENABLE_PTS_DIRECT_USER_VALUES") != NULL;
}

KLS_DEFINE_ENV_FLAG(kls_refactor_l_pattern_stats_env_enabled,
                    "KLS_ENABLE_REFACTOR_L_PATTERN_STATS")

KLS_DEFINE_ENV_FLAG(kls_refactor_stream_stats_env_enabled,
                    "KLS_ENABLE_REFACTOR_STREAM_STATS")

static int kls_refactor_l_index32_env_enabled(void) {
  const char *value = getenv("KLS_ENABLE_REFACTOR_L_INDEX32");
  return value == NULL || value[0] == '\0' ||
         !(value[0] == '0' && value[1] == '\0');
}

static int kls_refactor_u_index32_env_enabled(void) {
  const char *value = getenv("KLS_ENABLE_REFACTOR_U_INDEX32");
  return value == NULL || value[0] == '\0' ||
         !(value[0] == '0' && value[1] == '\0');
}

static int kls_ensure_refactor_l_index32_cache(kls_solver *solver) {
  if (!kls_refactor_l_index32_env_enabled()) {
    return 1;
  }

  if (solver == NULL || solver->numeric == NULL ||
      solver->symbolic == NULL || solver->symbolic->R == NULL ||
      solver->numeric->Llen == NULL ||
      solver->refactor_lu_pointer_count != solver->n ||
      solver->refactor_l_indices == NULL ||
      solver->n <= 0) {
    return 0;
  }
  if (solver->n > (UF_long)INT32_MAX) {
    return 1;
  }
  if (solver->refactor_l_indices32 != NULL) {
    return 1;
  }

  UF_long l_index32_entries = 0;
  for (UF_long col = 0; col < solver->n; ++col) {
    if (solver->refactor_l_indices[col] == NULL) {
      continue;
    }
    const UF_long len = solver->numeric->Llen[col];
    if (len > UF_long_max - l_index32_entries) {
      return 0;
    }
    l_index32_entries += len;
  }

  if (solver->n > (UF_long)(SIZE_MAX / sizeof(int32_t *))) {
    return 0;
  }
  int32_t **l_indices32 =
    (int32_t **)calloc((size_t)solver->n, sizeof(*l_indices32));
  if (l_indices32 == NULL) {
    return 0;
  }
  int32_t *l_indices32_storage = NULL;
  if (l_index32_entries > 0u) {
    if (l_index32_entries >
        (UF_long)(SIZE_MAX / sizeof(*l_indices32_storage))) {
      free(l_indices32);
      return 0;
    }
    l_indices32_storage =
      (int32_t *)malloc((size_t)l_index32_entries *
                        sizeof(*l_indices32_storage));
    if (l_indices32_storage == NULL) {
      free(l_indices32);
      return 0;
    }
  }

  UF_long offset = 0;
  for (UF_long col = 0; col < solver->n; ++col) {
    const UF_long len = solver->numeric->Llen[col];
    const UF_long *rows = solver->refactor_l_indices[col];
    if (len == 0u || rows == NULL) {
      l_indices32[col] = NULL;
      continue;
    }
    if (offset > l_index32_entries - len) {
      free(l_indices32_storage);
      free(l_indices32);
      return 0;
    }
    l_indices32[col] = l_indices32_storage + offset;
    for (UF_long p = 0; p < len; ++p) {
      if (rows[p] > (UF_long)INT32_MAX) {
        free(l_indices32_storage);
        free(l_indices32);
        return 0;
      }
      l_indices32_storage[offset + p] = (int32_t)rows[p];
    }
    offset += len;
  }

  solver->refactor_l_indices32 = l_indices32;
  solver->refactor_l_indices32_storage = l_indices32_storage;
  solver->refactor_l_indices32_count = l_index32_entries;
  solver->refactor_l_index32_enabled = 1;
  return 1;
}

static int kls_ensure_refactor_u_index32_cache(kls_solver *solver) {
  if (!kls_refactor_u_index32_env_enabled()) {
    return 1;
  }

  if (solver == NULL || solver->numeric == NULL ||
      solver->numeric->Ulen == NULL ||
      solver->refactor_lu_pointer_count != solver->n ||
      solver->refactor_u_indices == NULL ||
      solver->n <= 0) {
    return 0;
  }
  if (solver->n > (UF_long)INT32_MAX) {
    return 1;
  }
  if (solver->refactor_u_indices32 != NULL) {
    return 1;
  }

  UF_long u_index32_entries = 0;
  for (UF_long col = 0; col < solver->n; ++col) {
    if (solver->refactor_u_indices[col] == NULL) {
      continue;
    }
    const UF_long len = solver->numeric->Ulen[col];
    if (len > UF_long_max - u_index32_entries) {
      return 0;
    }
    u_index32_entries += len;
  }

  if (solver->n > (UF_long)(SIZE_MAX / sizeof(int32_t *))) {
    return 0;
  }
  int32_t **u_indices32 =
    (int32_t **)calloc((size_t)solver->n, sizeof(*u_indices32));
  if (u_indices32 == NULL) {
    return 0;
  }
  int32_t *u_indices32_storage = NULL;
  if (u_index32_entries > 0u) {
    if (u_index32_entries >
        (UF_long)(SIZE_MAX / sizeof(*u_indices32_storage))) {
      free(u_indices32);
      return 0;
    }
    u_indices32_storage =
      (int32_t *)malloc((size_t)u_index32_entries *
                        sizeof(*u_indices32_storage));
    if (u_indices32_storage == NULL) {
      free(u_indices32);
      return 0;
    }
  }

  UF_long offset = 0;
  for (UF_long col = 0; col < solver->n; ++col) {
    const UF_long len = solver->numeric->Ulen[col];
    const UF_long *rows = solver->refactor_u_indices[col];
    if (len == 0u || rows == NULL) {
      u_indices32[col] = NULL;
      continue;
    }
    if (offset > u_index32_entries - len) {
      free(u_indices32_storage);
      free(u_indices32);
      return 0;
    }
    u_indices32[col] = u_indices32_storage + offset;
    for (UF_long p = 0; p < len; ++p) {
      if (rows[p] > (UF_long)INT32_MAX) {
        free(u_indices32_storage);
        free(u_indices32);
        return 0;
      }
      u_indices32_storage[offset + p] = (int32_t)rows[p];
    }
    offset += len;
  }

  solver->refactor_u_indices32 = u_indices32;
  solver->refactor_u_indices32_storage = u_indices32_storage;
  solver->refactor_u_indices32_count = u_index32_entries;
  solver->refactor_u_index32_enabled = 1;
  return 1;
}

typedef struct kls_l_sorted_entry32 {
  int32_t row;
  int32_t pos;
} kls_l_sorted_entry32;

static int kls_compare_l_sorted_entry32(const void *a, const void *b) {
  const kls_l_sorted_entry32 *ea = (const kls_l_sorted_entry32 *)a;
  const kls_l_sorted_entry32 *eb = (const kls_l_sorted_entry32 *)b;
  if (ea->row < eb->row) {
    return -1;
  }
  if (ea->row > eb->row) {
    return 1;
  }
  return (ea->pos > eb->pos) - (ea->pos < eb->pos);
}

static void kls_free_refactor_l_sorted_cache(kls_solver *solver) {
  if (solver == NULL) {
    return;
  }
  free(solver->refactor_l_sorted_indices32);
  free(solver->refactor_l_sorted_indices32_storage);
  free(solver->refactor_l_sorted_pos32);
  free(solver->refactor_l_sorted_pos32_storage);
  free(solver->refactor_l_sorted_values);
  free(solver->refactor_l_sorted_values_storage);
  solver->refactor_l_sorted_indices32 = NULL;
  solver->refactor_l_sorted_indices32_storage = NULL;
  solver->refactor_l_sorted_pos32 = NULL;
  solver->refactor_l_sorted_pos32_storage = NULL;
  solver->refactor_l_sorted_values = NULL;
  solver->refactor_l_sorted_values_storage = NULL;
  solver->refactor_l_sorted_enabled = 0;
  solver->refactor_l_sorted_entries = 0u;
}

static int kls_ensure_refactor_l_sorted_cache(kls_solver *solver) {
  if (!kls_refactor_btf_scalar_run_exec_enabled(solver)) {
    return 1;
  }
  if (solver == NULL || solver->numeric == NULL ||
      solver->numeric->Llen == NULL ||
      solver->refactor_lu_pointer_count != solver->n ||
      solver->refactor_l_indices == NULL ||
      solver->refactor_l_values == NULL || solver->n <= 0) {
    return 0;
  }
  if (solver->n > (UF_long)INT32_MAX) {
    return 1;
  }
  if (solver->refactor_l_sorted_indices32 != NULL &&
      solver->refactor_l_sorted_pos32 != NULL &&
      solver->refactor_l_sorted_values != NULL) {
    solver->refactor_l_sorted_enabled = 1;
    return 1;
  }

  kls_free_refactor_l_sorted_cache(solver);
  UF_long l_entries = 0u;
  UF_long max_len = 0u;
  for (UF_long col = 0; col < solver->n; ++col) {
    const UF_long len = solver->numeric->Llen[col];
    const UF_long *rows = solver->refactor_l_indices[col];
    if (len == 0u || rows == NULL) {
      continue;
    }
    if (len > (UF_long)INT32_MAX || len > UF_long_max - l_entries) {
      return 0;
    }
    l_entries += len;
    if (len > max_len) {
      max_len = len;
    }
  }

  if (solver->n > (UF_long)(KLS_MAX_ALLOCATION / sizeof(int32_t *)) ||
      solver->n > (UF_long)(KLS_MAX_ALLOCATION / sizeof(double *))) {
    return 0;
  }
  const size_t n_size = (size_t)solver->n;
  if (n_size > (size_t)INT32_MAX) {
    return 0;
  }
  int32_t **sorted_rows =
    (int32_t **)calloc(n_size, sizeof(*sorted_rows));
  int32_t **sorted_pos =
    (int32_t **)calloc(n_size, sizeof(*sorted_pos));
  double **sorted_values =
    (double **)calloc(n_size, sizeof(*sorted_values));
  if (sorted_rows == NULL || sorted_pos == NULL || sorted_values == NULL) {
    free(sorted_rows);
    free(sorted_pos);
    free(sorted_values);
    return 0;
  }

  int32_t *row_storage = NULL;
  int32_t *pos_storage = NULL;
  double *value_storage = NULL;
  if (l_entries > 0u) {
    if (l_entries > (UF_long)(KLS_MAX_ALLOCATION / sizeof(*row_storage)) ||
        l_entries > (UF_long)(KLS_MAX_ALLOCATION / sizeof(*pos_storage)) ||
        l_entries > (UF_long)(KLS_MAX_ALLOCATION / sizeof(*value_storage))) {
      free(sorted_rows);
      free(sorted_pos);
      free(sorted_values);
      return 0;
    }
    row_storage =
      (int32_t *)malloc((size_t)l_entries * sizeof(*row_storage));
    pos_storage =
      (int32_t *)malloc((size_t)l_entries * sizeof(*pos_storage));
    value_storage =
      (double *)calloc((size_t)l_entries, sizeof(*value_storage));
    if (row_storage == NULL || pos_storage == NULL || value_storage == NULL) {
      free(row_storage);
      free(pos_storage);
      free(value_storage);
      free(sorted_rows);
      free(sorted_pos);
      free(sorted_values);
      return 0;
    }
  }

  kls_l_sorted_entry32 *scratch = NULL;
  if (max_len > 0u) {
    if (max_len > (UF_long)(KLS_MAX_ALLOCATION / sizeof(*scratch))) {
      free(row_storage);
      free(pos_storage);
      free(value_storage);
      free(sorted_rows);
      free(sorted_pos);
      free(sorted_values);
      return 0;
    }
    scratch =
      (kls_l_sorted_entry32 *)malloc((size_t)max_len * sizeof(*scratch));
    if (scratch == NULL) {
      free(row_storage);
      free(pos_storage);
      free(value_storage);
      free(sorted_rows);
      free(sorted_pos);
      free(sorted_values);
      return 0;
    }
  }

  UF_long offset = 0u;
  for (UF_long col = 0; col < solver->n; ++col) {
    const UF_long len = solver->numeric->Llen[col];
    const UF_long *rows = solver->refactor_l_indices[col];
    if (len == 0u || rows == NULL) {
      continue;
    }
    if (offset > l_entries - len || len > (UF_long)INT32_MAX) {
      free(scratch);
      free(row_storage);
      free(pos_storage);
      free(value_storage);
      free(sorted_rows);
      free(sorted_pos);
      free(sorted_values);
      return 0;
    }
    UF_long block = KLS_KLU_EMPTY;
    if (solver->symbolic->nblocks == 1u) {
      block = 0u;
    } else if (solver->refactor_col_block != NULL) {
      block = solver->refactor_col_block[col];
    } else {
      block = kls_block_for_pivot(solver, col);
    }
    if (block >= solver->symbolic->nblocks ||
        col < solver->symbolic->R[block]) {
      free(scratch);
      free(row_storage);
      free(pos_storage);
      free(value_storage);
      free(sorted_rows);
      free(sorted_pos);
      free(sorted_values);
      return 0;
    }
    const UF_long local_col = col - solver->symbolic->R[block];
    for (UF_long p = 0u; p < len; ++p) {
      if (rows[p] > (UF_long)INT32_MAX || p > (UF_long)INT32_MAX) {
        free(scratch);
        free(row_storage);
        free(pos_storage);
        free(value_storage);
        free(sorted_rows);
        free(sorted_pos);
        free(sorted_values);
        return 0;
      }
      if (rows[p] <= local_col) {
        free(scratch);
        free(row_storage);
        free(pos_storage);
        free(value_storage);
        free(sorted_rows);
        free(sorted_pos);
        free(sorted_values);
        return 0;
      }
      scratch[p].row = (int32_t)rows[p];
      scratch[p].pos = (int32_t)p;
    }
    qsort(scratch, (size_t)len, sizeof(*scratch),
          kls_compare_l_sorted_entry32);
    sorted_rows[col] = row_storage + offset;
    sorted_pos[col] = pos_storage + offset;
    sorted_values[col] = value_storage + offset;
    for (UF_long p = 0u; p < len; ++p) {
      if (p > 0u && scratch[p].row <= scratch[p - 1u].row) {
        free(scratch);
        free(row_storage);
        free(pos_storage);
        free(value_storage);
        free(sorted_rows);
        free(sorted_pos);
        free(sorted_values);
        return 0;
      }
      row_storage[offset + p] = scratch[p].row;
      pos_storage[offset + p] = scratch[p].pos;
    }
    offset += len;
  }
  free(scratch);

  solver->refactor_l_sorted_indices32 = sorted_rows;
  solver->refactor_l_sorted_indices32_storage = row_storage;
  solver->refactor_l_sorted_pos32 = sorted_pos;
  solver->refactor_l_sorted_pos32_storage = pos_storage;
  solver->refactor_l_sorted_values = sorted_values;
  solver->refactor_l_sorted_values_storage = value_storage;
  solver->refactor_l_sorted_entries = l_entries;
  solver->refactor_l_sorted_enabled = 1;
  return 1;
}

typedef struct kls_btf_scalar_run_group_item {
  UF_long block;
  UF_long dep0_global;
  UF_long current_global;
  UF_long up;
  UF_long rows;
  UF_long entries;
} kls_btf_scalar_run_group_item;

typedef struct kls_btf_scalar_run_group_wake_item {
  UF_long run_end_global;
  UF_long group;
  UF_long member;
} kls_btf_scalar_run_group_wake_item;

typedef struct kls_btf_scalar_run_group_step_item {
  UF_long producer_global;
  UF_long group;
  UF_long member_begin;
  UF_long member_count;
} kls_btf_scalar_run_group_step_item;

static int kls_compare_btf_scalar_run_group_item(const void *a,
                                                 const void *b) {
  const kls_btf_scalar_run_group_item *ea =
    (const kls_btf_scalar_run_group_item *)a;
  const kls_btf_scalar_run_group_item *eb =
    (const kls_btf_scalar_run_group_item *)b;
  if (ea->block < eb->block) {
    return -1;
  }
  if (ea->block > eb->block) {
    return 1;
  }
  if (ea->dep0_global < eb->dep0_global) {
    return -1;
  }
  if (ea->dep0_global > eb->dep0_global) {
    return 1;
  }
  if (ea->current_global < eb->current_global) {
    return -1;
  }
  if (ea->current_global > eb->current_global) {
    return 1;
  }
  return (ea->up > eb->up) - (ea->up < eb->up);
}

static int kls_compare_btf_scalar_run_group_wake_item(const void *a,
                                                      const void *b) {
  const kls_btf_scalar_run_group_wake_item *ea =
    (const kls_btf_scalar_run_group_wake_item *)a;
  const kls_btf_scalar_run_group_wake_item *eb =
    (const kls_btf_scalar_run_group_wake_item *)b;
  if (ea->run_end_global < eb->run_end_global) {
    return -1;
  }
  if (ea->run_end_global > eb->run_end_global) {
    return 1;
  }
  if (ea->group < eb->group) {
    return -1;
  }
  if (ea->group > eb->group) {
    return 1;
  }
  return (ea->member > eb->member) - (ea->member < eb->member);
}

static int kls_compare_btf_scalar_run_group_step_item(const void *a,
                                                      const void *b) {
  const kls_btf_scalar_run_group_step_item *ea =
    (const kls_btf_scalar_run_group_step_item *)a;
  const kls_btf_scalar_run_group_step_item *eb =
    (const kls_btf_scalar_run_group_step_item *)b;
  if (ea->producer_global < eb->producer_global) {
    return -1;
  }
  if (ea->producer_global > eb->producer_global) {
    return 1;
  }
  return (ea->group > eb->group) - (ea->group < eb->group);
}

static void kls_clear_refactor_btf_scalar_run_group_descriptor(
  kls_solver *solver) {
  if (solver == NULL) {
    return;
  }
  free(solver->refactor_btf_scalar_run_group_ptr);
  free(solver->refactor_btf_scalar_run_group_dep0);
  free(solver->refactor_btf_scalar_run_group_max_rows_by_group);
  free(solver->refactor_btf_scalar_run_group_currents);
  free(solver->refactor_btf_scalar_run_group_current_up);
  free(solver->refactor_btf_scalar_run_group_current_rows);
  free(solver->refactor_btf_scalar_run_group_step_ptr);
  free(solver->refactor_btf_scalar_run_group_step_dep);
  free(solver->refactor_btf_scalar_run_group_step_groups);
  free(solver->refactor_btf_scalar_run_group_step_member_ptr);
  free(solver->refactor_btf_scalar_run_group_step_members);
  free(solver->refactor_btf_scalar_run_group_live_step_ptr);
  free(solver->refactor_btf_scalar_run_group_live_step_row_idx);
  free(solver->refactor_btf_scalar_run_group_wake_ptr);
  free(solver->refactor_btf_scalar_run_group_wake_dep);
  free(solver->refactor_btf_scalar_run_group_wake_groups);
  free(solver->refactor_btf_scalar_run_group_wake_members);
  free(solver->refactor_btf_scalar_run_group_state_ptr);
  free(solver->refactor_btf_scalar_run_group_state_rows);
  free(solver->refactor_btf_scalar_run_group_state_current_ptr);
  free(solver->refactor_btf_scalar_run_group_state_current_members);
  free(solver->refactor_btf_scalar_run_group_state_current_best_member);
  free(solver->refactor_btf_scalar_run_group_state_value_base);
  free(solver->refactor_btf_scalar_run_group_live);
  free(solver->refactor_btf_scalar_run_group_wake_live);
  free(solver->refactor_btf_scalar_run_group_state_next);
  solver->refactor_btf_scalar_run_group_ptr = NULL;
  solver->refactor_btf_scalar_run_group_dep0 = NULL;
  solver->refactor_btf_scalar_run_group_max_rows_by_group = NULL;
  solver->refactor_btf_scalar_run_group_currents = NULL;
  solver->refactor_btf_scalar_run_group_current_up = NULL;
  solver->refactor_btf_scalar_run_group_current_rows = NULL;
  solver->refactor_btf_scalar_run_group_step_ptr = NULL;
  solver->refactor_btf_scalar_run_group_step_dep = NULL;
  solver->refactor_btf_scalar_run_group_step_groups = NULL;
  solver->refactor_btf_scalar_run_group_step_member_ptr = NULL;
  solver->refactor_btf_scalar_run_group_step_members = NULL;
  solver->refactor_btf_scalar_run_group_step_dep_count = 0u;
  solver->refactor_btf_scalar_run_group_step_count = 0u;
  solver->refactor_btf_scalar_run_group_step_member_total = 0u;
  solver->refactor_btf_scalar_run_group_live_step_ptr = NULL;
  solver->refactor_btf_scalar_run_group_live_step_row_idx = NULL;
  solver->refactor_btf_scalar_run_group_live_step_count = 0u;
  solver->refactor_btf_scalar_run_group_live_step_current_total = 0u;
  solver->refactor_btf_scalar_run_group_live_step_rows = 0u;
  solver->refactor_btf_scalar_run_group_live_step_unique_rows = 0u;
  solver->refactor_btf_scalar_run_group_live_step_reused_rows = 0u;
  solver->refactor_btf_scalar_run_group_live_step_max_currents = 0u;
  solver->refactor_btf_scalar_run_group_live_step_max_rows = 0u;
  solver->refactor_btf_scalar_run_group_live_step_max_unique_rows = 0u;
  solver->refactor_btf_scalar_run_group_live_step_stored_rows = 0u;
  solver->refactor_btf_scalar_run_group_live_step_storage_limited = 0u;
  solver->refactor_btf_scalar_run_group_live_step_runtime_full_step_count = 0u;
  solver->refactor_btf_scalar_run_group_live_step_runtime_full_current_count =
    0u;
  solver->refactor_btf_scalar_run_group_live_step_runtime_full_rows = 0u;
  solver->refactor_btf_scalar_run_group_live_step_runtime_partial_step_count =
    0u;
  solver
    ->refactor_btf_scalar_run_group_live_step_runtime_partial_current_count =
      0u;
  solver
    ->refactor_btf_scalar_run_group_live_step_runtime_partial_member_count =
      0u;
  solver->refactor_btf_scalar_run_group_live_step_runtime_partial_rows = 0u;
  solver->refactor_btf_scalar_run_group_wake_ptr = NULL;
  solver->refactor_btf_scalar_run_group_wake_dep = NULL;
  solver->refactor_btf_scalar_run_group_wake_groups = NULL;
  solver->refactor_btf_scalar_run_group_wake_members = NULL;
  solver->refactor_btf_scalar_run_group_wake_count = 0u;
  solver->refactor_btf_scalar_run_group_wake_member_total = 0u;
  solver->refactor_btf_scalar_run_group_state_ptr = NULL;
  solver->refactor_btf_scalar_run_group_state_rows = NULL;
  solver->refactor_btf_scalar_run_group_state_current_ptr = NULL;
  solver->refactor_btf_scalar_run_group_state_current_members = NULL;
  solver->refactor_btf_scalar_run_group_state_current_best_member = NULL;
  solver->refactor_btf_scalar_run_group_state_value_base = NULL;
  solver->refactor_btf_scalar_run_group_state_value_base_count = 0u;
  solver->refactor_btf_scalar_run_group_state_rows_total = 0u;
  solver->refactor_btf_scalar_run_group_state_max_rows = 0u;
  solver->refactor_btf_scalar_run_group_state_current_count = 0u;
  solver->refactor_btf_scalar_run_group_state_best_skip_total = 0u;
  solver->refactor_btf_scalar_run_group_state_max_best_skip = 0u;
  solver->refactor_btf_scalar_run_group_state_guard_lower_bound_rows = 0u;
  solver->refactor_btf_scalar_run_group_state_guard_lower_bound_rejected = 0u;
  solver->refactor_btf_scalar_run_group_live = NULL;
  solver->refactor_btf_scalar_run_group_live_count = 0u;
  solver->refactor_btf_scalar_run_group_wake_live = NULL;
  solver->refactor_btf_scalar_run_group_wake_live_count = 0u;
  solver->refactor_btf_scalar_run_group_state_next = NULL;
  solver->refactor_btf_scalar_run_group_state_next_count = 0u;
  solver->refactor_btf_scalar_run_group_built = 0;
}

static void kls_reset_refactor_btf_scalar_run_group_stats(kls_solver *solver) {
  if (solver == NULL) {
    return;
  }
  solver->refactor_btf_scalar_run_group_count = 0;
  solver->refactor_btf_scalar_run_group_current_total = 0;
  solver->refactor_btf_scalar_run_group_multi_count = 0;
  solver->refactor_btf_scalar_run_group_multi_current_total = 0;
  solver->refactor_btf_scalar_run_group_rows = 0;
  solver->refactor_btf_scalar_run_group_reused_rows = 0;
  solver->refactor_btf_scalar_run_group_entries = 0;
  solver->refactor_btf_scalar_run_group_reused_entries = 0;
  solver->refactor_btf_scalar_run_group_max_currents = 0;
  solver->refactor_btf_scalar_run_group_max_rows = 0;
  solver->refactor_btf_scalar_run_group_producer_step_count = 0;
  solver->refactor_btf_scalar_run_group_producer_step_multi_count = 0;
  solver->refactor_btf_scalar_run_group_producer_step_active_currents = 0;
  solver->refactor_btf_scalar_run_group_producer_step_unique_entries = 0;
  solver->refactor_btf_scalar_run_group_producer_step_duplicate_entries = 0;
  solver->refactor_btf_scalar_run_group_producer_step_reused_entries = 0;
  solver->refactor_btf_scalar_run_group_producer_step_max_currents = 0;
  solver->refactor_btf_scalar_run_group_producer_index_count = 0;
  solver->refactor_btf_scalar_run_group_producer_index_max_steps = 0;
  solver->refactor_btf_scalar_run_group_live_state_groups = 0;
  solver->refactor_btf_scalar_run_group_live_state_currents = 0;
  solver->refactor_btf_scalar_run_group_live_state_rows = 0;
  solver->refactor_btf_scalar_run_group_live_state_unique_rows = 0;
  solver->refactor_btf_scalar_run_group_live_state_reused_rows = 0;
  solver->refactor_btf_scalar_run_group_live_state_max_currents = 0;
  solver->refactor_btf_scalar_run_group_live_state_max_rows = 0;
  solver->refactor_btf_scalar_run_group_live_state_max_unique_rows = 0;
  solver->refactor_btf_scalar_run_group_live_step_count = 0;
  solver->refactor_btf_scalar_run_group_live_step_current_total = 0;
  solver->refactor_btf_scalar_run_group_live_step_rows = 0;
  solver->refactor_btf_scalar_run_group_live_step_unique_rows = 0;
  solver->refactor_btf_scalar_run_group_live_step_reused_rows = 0;
  solver->refactor_btf_scalar_run_group_live_step_max_currents = 0;
  solver->refactor_btf_scalar_run_group_live_step_max_rows = 0;
  solver->refactor_btf_scalar_run_group_live_step_max_unique_rows = 0;
  solver->refactor_btf_scalar_run_group_live_step_stored_rows = 0;
  solver->refactor_btf_scalar_run_group_live_step_storage_limited = 0;
  solver->refactor_btf_scalar_run_group_live_step_runtime_full_step_count = 0;
  solver->refactor_btf_scalar_run_group_live_step_runtime_full_current_count =
    0;
  solver->refactor_btf_scalar_run_group_live_step_runtime_full_rows = 0;
  solver->refactor_btf_scalar_run_group_live_step_runtime_partial_step_count =
    0;
  solver
    ->refactor_btf_scalar_run_group_live_step_runtime_partial_current_count =
      0;
  solver
    ->refactor_btf_scalar_run_group_live_step_runtime_partial_member_count =
      0;
  solver->refactor_btf_scalar_run_group_live_step_runtime_partial_rows = 0;
  solver->refactor_btf_scalar_run_group_live_step_runtime_partial_active_rows =
    0;
  solver
    ->refactor_btf_scalar_run_group_live_step_runtime_partial_active_unique_rows =
      0;
  solver->refactor_btf_scalar_run_group_wake_count = 0;
  solver->refactor_btf_scalar_run_group_wake_member_total = 0;
  solver->refactor_btf_scalar_run_group_state_rows_total = 0;
  solver->refactor_btf_scalar_run_group_state_max_rows = 0;
  solver->refactor_btf_scalar_run_group_state_current_count = 0;
  solver->refactor_btf_scalar_run_group_state_best_skip_total = 0;
  solver->refactor_btf_scalar_run_group_state_max_best_skip = 0;
  solver->refactor_btf_scalar_run_group_state_guard_lower_bound_rows = 0;
  solver->refactor_btf_scalar_run_group_state_guard_lower_bound_rejected = 0;
  solver->refactor_btf_scalar_run_group_state_materialized_current_count = 0;
  solver->refactor_btf_scalar_run_group_state_materialized_rows = 0;
  solver->refactor_btf_scalar_run_group_state_materialized_prefix_deps = 0;
  solver->refactor_btf_scalar_run_group_state_advanced_current_count = 0;
  solver->refactor_btf_scalar_run_group_state_advanced_rows = 0;
  solver->refactor_btf_scalar_run_group_state_advanced_entries = 0;
  solver->refactor_btf_scalar_run_group_state_advance_batch_group_count = 0;
  solver->refactor_btf_scalar_run_group_state_advance_batch_current_count = 0;
  solver->refactor_btf_scalar_run_group_state_advance_batch_unique_entries = 0;
  solver
    ->refactor_btf_scalar_run_group_state_advance_batch_duplicate_entries = 0;
  solver->refactor_btf_scalar_run_group_state_advance_batch_max_currents = 0;
  solver->refactor_btf_scalar_run_group_state_step_advance_trigger_count = 0;
  solver->refactor_btf_scalar_run_group_state_step_advance_steps = 0;
  solver->refactor_btf_scalar_run_group_state_step_advance_current_count = 0;
  solver->refactor_btf_scalar_run_group_state_step_advance_rows = 0;
  solver->refactor_btf_scalar_run_group_state_step_advance_entries = 0;
  solver
    ->refactor_btf_scalar_run_group_state_step_advance_ready_current_count =
      0;
  solver->refactor_btf_scalar_run_group_state_step_advance_reject_count = 0;
  solver->refactor_btf_scalar_run_group_state_step_batch_step_count = 0;
  solver->refactor_btf_scalar_run_group_state_step_batch_current_count = 0;
  solver->refactor_btf_scalar_run_group_state_step_batch_entries = 0;
  solver->refactor_btf_scalar_run_group_state_step_batch_state_rows = 0;
  solver->refactor_btf_scalar_run_group_state_step_batch_unique_state_rows = 0;
  solver->refactor_btf_scalar_run_group_state_step_window_round_count = 0;
  solver->refactor_btf_scalar_run_group_state_step_window_current_count = 0;
  solver->refactor_btf_scalar_run_group_state_step_window_entries = 0;
  solver->refactor_btf_scalar_run_group_state_reject_count = 0;
  solver->refactor_btf_scalar_run_group_state_exec_current_count = 0;
  solver->refactor_btf_scalar_run_group_state_exec_skipped_deps = 0;
  solver->refactor_btf_scalar_run_group_state_exec_restored_rows = 0;
  solver->refactor_btf_scalar_run_group_state_exec_terminal_current_count = 0;
  solver
    ->refactor_btf_scalar_run_group_state_exec_wake_terminal_current_count = 0;
  solver
    ->refactor_btf_scalar_run_group_state_exec_wake_terminal_candidate_current_count =
      0;
  solver
    ->refactor_btf_scalar_run_group_state_exec_wake_terminal_nonpipeline_current_count =
      0;
  solver
    ->refactor_btf_scalar_run_group_state_exec_wake_terminal_dependency_miss_current_count =
      0;
  solver
    ->refactor_btf_scalar_run_group_state_exec_wake_terminal_claim_miss_current_count =
      0;
  solver->refactor_btf_scalar_run_group_state_exec_remaining_deps = 0;
  solver->refactor_btf_scalar_run_group_state_exec_remaining_entries = 0;
  solver->refactor_btf_scalar_run_group_state_exec_reject_count = 0;
  solver->refactor_btf_scalar_run_group_state_exec_dispatch_bypass_current_count =
    0;
  solver
    ->refactor_btf_scalar_run_group_state_exec_dispatch_bypass_ready_current_count =
      0;
  solver->refactor_btf_scalar_run_group_state_exec_not_ready_current_count =
    0;
  solver->refactor_btf_scalar_run_group_state_exec_owned_ready_current_count =
    0;
  solver->refactor_btf_scalar_run_group_state_exec_late_ready_current_count =
    0;
}

static UF_long kls_saturating_add_uf_long(UF_long a, UF_long b) {
  return b > UF_long_max - a ? UF_long_max : a + b;
}

static UF_long kls_saturating_mul_uf_long(UF_long a, UF_long b) {
  if (a == 0u || b == 0u) {
    return 0u;
  }
  return a > UF_long_max / b ? UF_long_max : a * b;
}

static UF_long kls_btf_scalar_run_structural_len(
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

static UF_long kls_btf_scalar_run_l_entries(const kls_solver *solver,
                                            UF_long dep0_global,
                                            UF_long rows) {
  if (solver == NULL || solver->numeric == NULL ||
      solver->numeric->Llen == NULL || dep0_global > solver->n ||
      rows > solver->n - dep0_global) {
    return 0u;
  }
  UF_long entries = 0u;
  for (UF_long p = 0u; p < rows; ++p) {
    const UF_long len = solver->numeric->Llen[dep0_global + p];
    if (len > UF_long_max - entries) {
      return UF_long_max;
    }
    entries += len;
  }
  return entries;
}

static atomic_uint *kls_prepare_refactor_btf_scalar_run_group_live(
  kls_solver *solver) {
  if (solver == NULL || !solver->refactor_btf_scalar_run_group_built ||
      solver->refactor_btf_scalar_run_group_multi_count == 0u) {
    return NULL;
  }
  const UF_long count = solver->refactor_btf_scalar_run_group_multi_count;
  if (count > (UF_long)(KLS_MAX_ALLOCATION / sizeof(atomic_uint))) {
    return NULL;
  }
  if (solver->refactor_btf_scalar_run_group_live != NULL &&
      solver->refactor_btf_scalar_run_group_live_count != count) {
    free(solver->refactor_btf_scalar_run_group_live);
    solver->refactor_btf_scalar_run_group_live = NULL;
    solver->refactor_btf_scalar_run_group_live_count = 0u;
  }
  int fresh = 0;
  if (solver->refactor_btf_scalar_run_group_live == NULL) {
    solver->refactor_btf_scalar_run_group_live =
      (atomic_uint *)malloc((size_t)count *
                            sizeof(*solver->refactor_btf_scalar_run_group_live));
    if (solver->refactor_btf_scalar_run_group_live == NULL) {
      return NULL;
    }
    solver->refactor_btf_scalar_run_group_live_count = count;
    fresh = 1;
  }
  for (UF_long group = 0u; group < count; ++group) {
    if (fresh) {
      atomic_init(&solver->refactor_btf_scalar_run_group_live[group], 0u);
    } else {
      atomic_store_explicit(
        &solver->refactor_btf_scalar_run_group_live[group], 0u,
        memory_order_relaxed);
    }
  }
  return solver->refactor_btf_scalar_run_group_live;
}

static atomic_uint *kls_prepare_refactor_btf_scalar_run_group_wake_live(
  kls_solver *solver) {
  if (solver == NULL || !solver->refactor_btf_scalar_run_group_built ||
      solver->refactor_btf_scalar_run_group_multi_current_total == 0u) {
    return NULL;
  }
  const UF_long count =
    solver->refactor_btf_scalar_run_group_multi_current_total;
  if (count > (UF_long)(KLS_MAX_ALLOCATION / sizeof(atomic_uint))) {
    return NULL;
  }
  if (solver->refactor_btf_scalar_run_group_wake_live != NULL &&
      solver->refactor_btf_scalar_run_group_wake_live_count != count) {
    free(solver->refactor_btf_scalar_run_group_wake_live);
    solver->refactor_btf_scalar_run_group_wake_live = NULL;
    solver->refactor_btf_scalar_run_group_wake_live_count = 0u;
  }
  int fresh = 0;
  if (solver->refactor_btf_scalar_run_group_wake_live == NULL) {
    solver->refactor_btf_scalar_run_group_wake_live =
      (atomic_uint *)malloc(
        (size_t)count *
        sizeof(*solver->refactor_btf_scalar_run_group_wake_live));
    if (solver->refactor_btf_scalar_run_group_wake_live == NULL) {
      return NULL;
    }
    solver->refactor_btf_scalar_run_group_wake_live_count = count;
    fresh = 1;
  }
  for (UF_long member = 0u; member < count; ++member) {
    if (fresh) {
      atomic_init(&solver->refactor_btf_scalar_run_group_wake_live[member],
                  0u);
    } else {
      atomic_store_explicit(
        &solver->refactor_btf_scalar_run_group_wake_live[member], 0u,
        memory_order_relaxed);
    }
  }
  return solver->refactor_btf_scalar_run_group_wake_live;
}

static UF_long *kls_prepare_refactor_btf_scalar_run_group_state_next(
  kls_solver *solver) {
  if (solver == NULL || !solver->refactor_btf_scalar_run_group_built ||
      solver->refactor_btf_scalar_run_group_multi_current_total == 0u) {
    return NULL;
  }
  const UF_long count =
    solver->refactor_btf_scalar_run_group_multi_current_total;
  if (count > (UF_long)(KLS_MAX_ALLOCATION / sizeof(UF_long))) {
    return NULL;
  }
  if (solver->refactor_btf_scalar_run_group_state_next != NULL &&
      solver->refactor_btf_scalar_run_group_state_next_count != count) {
    free(solver->refactor_btf_scalar_run_group_state_next);
    solver->refactor_btf_scalar_run_group_state_next = NULL;
    solver->refactor_btf_scalar_run_group_state_next_count = 0u;
  }
  if (solver->refactor_btf_scalar_run_group_state_next == NULL) {
    solver->refactor_btf_scalar_run_group_state_next =
      (UF_long *)malloc(
        (size_t)count *
        sizeof(*solver->refactor_btf_scalar_run_group_state_next));
    if (solver->refactor_btf_scalar_run_group_state_next == NULL) {
      return NULL;
    }
    solver->refactor_btf_scalar_run_group_state_next_count = count;
  }
  memset(solver->refactor_btf_scalar_run_group_state_next, 0,
         (size_t)count *
           sizeof(*solver->refactor_btf_scalar_run_group_state_next));
  return solver->refactor_btf_scalar_run_group_state_next;
}

static double *kls_prepare_refactor_btf_scalar_run_group_state_values(
  kls_solver *solver,
  int need_u_values,
  int compact_selected,
  const UF_long **value_base_out,
  UF_long *value_base_count_out,
  double **u_values_out) {
  KLS_SET_OPTIONAL_OUTPUT(value_base_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(value_base_count_out, 0u);
  KLS_SET_OPTIONAL_OUTPUT(u_values_out, NULL);
  if (solver == NULL || !solver->refactor_btf_scalar_run_group_built ||
      solver->refactor_btf_scalar_run_group_state_rows_total == 0u ||
      solver->refactor_btf_scalar_run_group_state_ptr == NULL ||
      solver->refactor_btf_scalar_run_group_state_rows == NULL) {
    return NULL;
  }
  UF_long count = solver->refactor_btf_scalar_run_group_state_rows_total;
  if (compact_selected) {
    if (solver->refactor_btf_scalar_run_group_state_current_best_member ==
          NULL ||
        solver->refactor_btf_scalar_run_group_multi_current_total == 0u ||
        solver->n <= 0) {
      return NULL;
    }
    const UF_long member_count =
      solver->refactor_btf_scalar_run_group_multi_current_total;
    if (member_count > (UF_long)(KLS_MAX_ALLOCATION / sizeof(UF_long))) {
      return NULL;
    }
    if (solver->refactor_btf_scalar_run_group_state_value_base != NULL &&
        solver->refactor_btf_scalar_run_group_state_value_base_count !=
          member_count) {
      free(solver->refactor_btf_scalar_run_group_state_value_base);
      solver->refactor_btf_scalar_run_group_state_value_base = NULL;
      solver->refactor_btf_scalar_run_group_state_value_base_count = 0u;
    }
    if (solver->refactor_btf_scalar_run_group_state_value_base == NULL) {
      solver->refactor_btf_scalar_run_group_state_value_base =
        (UF_long *)malloc(
          (size_t)member_count *
          sizeof(*solver->refactor_btf_scalar_run_group_state_value_base));
      if (solver->refactor_btf_scalar_run_group_state_value_base == NULL) {
        return NULL;
      }
      solver->refactor_btf_scalar_run_group_state_value_base_count =
        member_count;
    }
    UF_long *value_base =
      solver->refactor_btf_scalar_run_group_state_value_base;
    for (UF_long member = 0u; member < member_count; ++member) {
      value_base[member] = KLS_KLU_EMPTY;
    }
    count = 0u;
    for (UF_long member = 0u; member < member_count; ++member) {
      if (solver->refactor_btf_scalar_run_group_currents == NULL) {
        return NULL;
      }
      const UF_long current =
        solver->refactor_btf_scalar_run_group_currents[member];
      if (current >= solver->n ||
          solver->refactor_btf_scalar_run_group_state_current_best_member[
            current] == KLS_KLU_EMPTY) {
        continue;
      }
      if (value_base[member] != KLS_KLU_EMPTY ||
          solver->refactor_btf_scalar_run_group_state_ptr[member] >
            solver->refactor_btf_scalar_run_group_state_ptr[member + 1u] ||
          solver->refactor_btf_scalar_run_group_state_ptr[member + 1u] >
            solver->refactor_btf_scalar_run_group_state_rows_total) {
        return NULL;
      }
      const UF_long rows =
        solver->refactor_btf_scalar_run_group_state_ptr[member + 1u] -
        solver->refactor_btf_scalar_run_group_state_ptr[member];
      if (rows > UF_long_max - count) {
        return NULL;
      }
      value_base[member] = count;
      count += rows;
    }
    if (count == 0u) {
      return NULL;
    }
    if (value_base_out != NULL) {
      *value_base_out = value_base;
    }
    if (value_base_count_out != NULL) {
      *value_base_count_out = member_count;
    }
  } else {
    free(solver->refactor_btf_scalar_run_group_state_value_base);
    solver->refactor_btf_scalar_run_group_state_value_base = NULL;
    solver->refactor_btf_scalar_run_group_state_value_base_count = 0u;
  }
  if (count > (UF_long)(KLS_MAX_ALLOCATION / sizeof(double))) {
    return NULL;
  }
  if (solver->egraph_btf_scalar_run_group_state_values == NULL ||
      solver->egraph_btf_scalar_run_group_state_size < count) {
    double *values =
      (double *)realloc(solver->egraph_btf_scalar_run_group_state_values,
                        (size_t)count *
                          sizeof(*solver
                                     ->egraph_btf_scalar_run_group_state_values));
    if (values == NULL) {
      return NULL;
    }
    solver->egraph_btf_scalar_run_group_state_values = values;
    solver->egraph_btf_scalar_run_group_state_size = count;
  }
  if (need_u_values) {
    double *u_values =
      (double *)realloc(solver->egraph_btf_scalar_run_group_state_u_values,
                        (size_t)count *
                          sizeof(*solver
                                     ->egraph_btf_scalar_run_group_state_u_values));
    if (u_values == NULL) {
      return NULL;
    }
    solver->egraph_btf_scalar_run_group_state_u_values = u_values;
  } else {
    free(solver->egraph_btf_scalar_run_group_state_u_values);
    solver->egraph_btf_scalar_run_group_state_u_values = NULL;
  }
  memset(solver->egraph_btf_scalar_run_group_state_values, 0,
         (size_t)count *
           sizeof(*solver->egraph_btf_scalar_run_group_state_values));
  if (need_u_values) {
    memset(solver->egraph_btf_scalar_run_group_state_u_values, 0,
           (size_t)count *
             sizeof(*solver->egraph_btf_scalar_run_group_state_u_values));
    if (u_values_out != NULL) {
      *u_values_out = solver->egraph_btf_scalar_run_group_state_u_values;
    }
  }
  return solver->egraph_btf_scalar_run_group_state_values;
}

static int kls_mark_btf_scalar_run_group_state_l_rows(
  const kls_solver *solver,
  UF_long dep_global,
  UF_long block_size,
  UF_long *stamp_workspace,
  UF_long stamp,
  UF_long *rows,
  UF_long row_capacity,
  UF_long *row_count_io) {
  if (solver == NULL || solver->numeric == NULL ||
      solver->numeric->Llen == NULL || dep_global >= solver->n ||
      stamp_workspace == NULL || row_count_io == NULL ||
      block_size == 0u) {
    return 0;
  }
  const UF_long len = solver->numeric->Llen[dep_global];
  const int32_t *rows32 =
    solver->refactor_l_indices32 != NULL
      ? solver->refactor_l_indices32[dep_global] : NULL;
  const UF_long *lrows =
    rows32 == NULL && solver->refactor_l_indices != NULL
      ? solver->refactor_l_indices[dep_global] : NULL;
  if (len > 0u && rows32 == NULL && lrows == NULL) {
    return 0;
  }
  for (UF_long p = 0u; p < len; ++p) {
    const UF_long row = rows32 != NULL ? (UF_long)rows32[p] : lrows[p];
    if (row >= block_size ||
        !kls_mark_unique_state_row(
          row, stamp_workspace, block_size, stamp, rows, row_capacity,
          row_count_io)) {
      return 0;
    }
  }
  return 1;
}

static int kls_collect_btf_scalar_run_group_state_rows(
  const kls_solver *solver,
  UF_long member,
  UF_long *stamp_workspace,
  UF_long stamp,
  UF_long *rows,
  UF_long row_capacity,
  UF_long *row_count_out) {
  KLS_SET_OPTIONAL_OUTPUT(row_count_out, 0u);
  if (solver == NULL || solver->symbolic == NULL ||
      solver->symbolic->R == NULL || solver->numeric == NULL ||
      solver->numeric->Ulen == NULL || solver->numeric->Llen == NULL ||
      solver->refactor_col_ptr == NULL ||
      solver->refactor_block_start == NULL ||
      solver->refactor_btf_scalar_run_group_currents == NULL ||
      solver->refactor_btf_scalar_run_group_current_up == NULL ||
      solver->refactor_btf_scalar_run_group_current_rows == NULL ||
      stamp_workspace == NULL || stamp == 0u || row_count_out == NULL ||
      member >= solver->refactor_btf_scalar_run_group_multi_current_total) {
    return 0;
  }
  const UF_long current =
    solver->refactor_btf_scalar_run_group_currents[member];
  if (current >= solver->n) {
    return 0;
  }
  UF_long block = 0u;
  if (solver->symbolic->nblocks != 1u) {
    if (solver->refactor_col_block == NULL) {
      return 0;
    }
    block = solver->refactor_col_block[current];
  }
  if (block >= solver->symbolic->nblocks) {
    return 0;
  }
  const UF_long k1 = solver->symbolic->R[block];
  const UF_long k2 = solver->symbolic->R[block + 1u];
  if (k2 <= k1 || current < k1 || current >= k2) {
    return 0;
  }
  const UF_long block_size = k2 - k1;
  const UF_long current_local = current - k1;
  UF_long row_count = 0u;

  const UF_long input_begin = solver->refactor_block_start[current];
  const UF_long input_end = solver->refactor_col_ptr[current + 1u];
  if (input_begin > input_end || input_end > solver->nnz) {
    return 0;
  }
  if (solver->refactor_row_idx32 != NULL) {
    const int32_t *input_rows32 = solver->refactor_row_idx32;
    for (UF_long p = input_begin; p < input_end; ++p) {
      const UF_long row_global = (UF_long)input_rows32[p];
      if (row_global < k1 || row_global >= k2) {
        return 0;
      }
      if (!kls_mark_unique_state_row(
            row_global - k1, stamp_workspace, block_size, stamp, rows,
            row_capacity, &row_count)) {
        return 0;
      }
    }
  } else {
    if (solver->refactor_row_idx == NULL) {
      return 0;
    }
    const UF_long *input_rows = solver->refactor_row_idx;
    for (UF_long p = input_begin; p < input_end; ++p) {
      const UF_long row_global = input_rows[p];
      if (row_global < k1 || row_global >= k2) {
        return 0;
      }
      if (!kls_mark_unique_state_row(
            row_global - k1, stamp_workspace, block_size, stamp, rows,
            row_capacity, &row_count)) {
        return 0;
      }
    }
  }

  const UF_long ulen = solver->numeric->Ulen[current];
  const int32_t *ui32 =
    solver->refactor_u_indices32 != NULL
      ? solver->refactor_u_indices32[current] : NULL;
  const UF_long *ui =
    ui32 == NULL && solver->refactor_u_indices != NULL
      ? solver->refactor_u_indices[current] : NULL;
  if (ulen > 0u && ui32 == NULL && ui == NULL) {
    return 0;
  }
  for (UF_long up = 0u; up < ulen; ++up) {
    const UF_long dep_local = ui32 != NULL ? (UF_long)ui32[up] : ui[up];
    if (dep_local >= current_local || dep_local >= block_size ||
        dep_local > UF_long_max - k1 ||
        !kls_mark_unique_state_row(
          dep_local, stamp_workspace, block_size, stamp, rows, row_capacity,
          &row_count) ||
        !kls_mark_btf_scalar_run_group_state_l_rows(
          solver, k1 + dep_local, block_size, stamp_workspace, stamp, rows,
          row_capacity, &row_count)) {
      return 0;
    }
  }

  if (!kls_mark_unique_state_row(
        current_local, stamp_workspace, block_size, stamp, rows, row_capacity,
        &row_count) ||
      !kls_mark_btf_scalar_run_group_state_l_rows(
        solver, current, block_size, stamp_workspace, stamp, rows,
        row_capacity, &row_count)) {
    return 0;
  }

  *row_count_out = row_count;
  return 1;
}

static int kls_count_refactor_btf_scalar_run_group_live_state_stats(
  kls_solver *solver) {
  if (solver == NULL ||
      !kls_refactor_btf_scalar_run_group_live_state_stats_env_enabled() ||
      !solver->refactor_btf_scalar_run_group_built ||
      solver->refactor_btf_scalar_run_group_multi_count == 0u ||
      solver->refactor_btf_scalar_run_group_multi_current_total == 0u ||
      solver->refactor_btf_scalar_run_group_ptr == NULL ||
      solver->refactor_btf_scalar_run_group_currents == NULL) {
    return 1;
  }
  const UF_long stamp_size = solver->n > 0u ? solver->n : 1u;
  if (stamp_size > (UF_long)(KLS_MAX_ALLOCATION / sizeof(UF_long))) {
    return 1;
  }
  UF_long *group_stamp_workspace =
    (UF_long *)calloc((size_t)stamp_size, sizeof(*group_stamp_workspace));
  UF_long *member_stamp_workspace =
    (UF_long *)calloc((size_t)stamp_size, sizeof(*member_stamp_workspace));
  if (group_stamp_workspace == NULL || member_stamp_workspace == NULL) {
    free(group_stamp_workspace);
    free(member_stamp_workspace);
    return 1;
  }

  UF_long groups = 0u;
  UF_long currents = 0u;
  UF_long rows = 0u;
  UF_long unique_rows = 0u;
  UF_long reused_rows = 0u;
  UF_long max_currents = 0u;
  UF_long max_rows = 0u;
  UF_long max_unique_rows = 0u;
  UF_long group_stamp_value = 0u;
  UF_long member_stamp_value = 0u;

  for (UF_long group = 0u;
       group < solver->refactor_btf_scalar_run_group_multi_count; ++group) {
    const UF_long begin = solver->refactor_btf_scalar_run_group_ptr[group];
    const UF_long end = solver->refactor_btf_scalar_run_group_ptr[group + 1u];
    if (begin > end ||
        end > solver->refactor_btf_scalar_run_group_multi_current_total) {
      free(group_stamp_workspace);
      free(member_stamp_workspace);
      return 1;
    }
    const UF_long group_currents = end - begin;
    if (group_currents <= 1u) {
      continue;
    }
    UF_long group_stamp = 0u;
    if (!kls_refactor_supernode_algorithm5_next_target_stamp(
          group_stamp_workspace, stamp_size, &group_stamp_value,
          &group_stamp)) {
      free(group_stamp_workspace);
      free(member_stamp_workspace);
      return 1;
    }
    UF_long group_rows = 0u;
    UF_long group_unique_rows = 0u;
    for (UF_long member = begin; member < end; ++member) {
      UF_long member_stamp = 0u;
      UF_long member_rows = 0u;
      UF_long member_unique_rows = 0u;
      if (!kls_refactor_supernode_algorithm5_next_target_stamp(
            member_stamp_workspace, stamp_size, &member_stamp_value,
            &member_stamp) ||
          !kls_collect_btf_scalar_run_group_state_rows(
            solver, member, member_stamp_workspace, member_stamp, NULL, 0u,
            &member_rows) ||
          !kls_collect_btf_scalar_run_group_state_rows(
            solver, member, group_stamp_workspace, group_stamp, NULL, 0u,
            &member_unique_rows)) {
        free(group_stamp_workspace);
        free(member_stamp_workspace);
        return 1;
      }
      group_rows = kls_saturating_add_uf_long(group_rows, member_rows);
      group_unique_rows =
        kls_saturating_add_uf_long(group_unique_rows, member_unique_rows);
    }

    groups = kls_saturating_add_uf_long(groups, 1u);
    currents = kls_saturating_add_uf_long(currents, group_currents);
    rows = kls_saturating_add_uf_long(rows, group_rows);
    unique_rows = kls_saturating_add_uf_long(unique_rows, group_unique_rows);
    if (group_rows > group_unique_rows) {
      reused_rows =
        kls_saturating_add_uf_long(reused_rows,
                                   group_rows - group_unique_rows);
    }
    if (group_currents > max_currents) {
      max_currents = group_currents;
    }
    if (group_rows > max_rows) {
      max_rows = group_rows;
    }
    if (group_unique_rows > max_unique_rows) {
      max_unique_rows = group_unique_rows;
    }
  }

  free(group_stamp_workspace);
  free(member_stamp_workspace);
  solver->refactor_btf_scalar_run_group_live_state_groups = groups;
  solver->refactor_btf_scalar_run_group_live_state_currents = currents;
  solver->refactor_btf_scalar_run_group_live_state_rows = rows;
  solver->refactor_btf_scalar_run_group_live_state_unique_rows = unique_rows;
  solver->refactor_btf_scalar_run_group_live_state_reused_rows = reused_rows;
  solver->refactor_btf_scalar_run_group_live_state_max_currents =
    max_currents;
  solver->refactor_btf_scalar_run_group_live_state_max_rows = max_rows;
  solver->refactor_btf_scalar_run_group_live_state_max_unique_rows =
    max_unique_rows;
  return 1;
}

static int kls_count_refactor_btf_scalar_run_group_state_skip_summary(
  const kls_solver *solver,
  UF_long *state_current_count_out,
  UF_long *state_best_skip_total_out,
  UF_long *state_max_best_skip_out,
  UF_long *member_skip_lower_bound_out) {
  KLS_SET_OPTIONAL_OUTPUT(state_current_count_out, 0u);
  KLS_SET_OPTIONAL_OUTPUT(state_best_skip_total_out, 0u);
  KLS_SET_OPTIONAL_OUTPUT(state_max_best_skip_out, 0u);
  KLS_SET_OPTIONAL_OUTPUT(member_skip_lower_bound_out, 0u);
  if (solver == NULL || state_current_count_out == NULL ||
      state_best_skip_total_out == NULL ||
      state_max_best_skip_out == NULL ||
      member_skip_lower_bound_out == NULL ||
      solver->refactor_btf_scalar_run_group_multi_current_total == 0u ||
      solver->refactor_btf_scalar_run_group_currents == NULL ||
      solver->refactor_btf_scalar_run_group_current_up == NULL ||
      solver->refactor_btf_scalar_run_group_current_rows == NULL ||
      solver->n == 0u ||
      solver->n > (UF_long)(KLS_MAX_ALLOCATION / sizeof(UF_long))) {
    return 0;
  }

  UF_long *best_skip_by_current =
    (UF_long *)calloc((size_t)solver->n, sizeof(*best_skip_by_current));
  if (best_skip_by_current == NULL) {
    return 0;
  }

  UF_long member_skip_lower_bound = 0u;
  const UF_long member_total =
    solver->refactor_btf_scalar_run_group_multi_current_total;
  for (UF_long member = 0u; member < member_total; ++member) {
    const UF_long current =
      solver->refactor_btf_scalar_run_group_currents[member];
    const UF_long current_up =
      solver->refactor_btf_scalar_run_group_current_up[member];
    const UF_long rows =
      solver->refactor_btf_scalar_run_group_current_rows[member];
    if (current >= solver->n || rows > UF_long_max - current_up) {
      free(best_skip_by_current);
      return 0;
    }
    const UF_long skip = current_up + rows;
    member_skip_lower_bound =
      kls_saturating_add_uf_long(member_skip_lower_bound, skip);
    if (skip > best_skip_by_current[current]) {
      best_skip_by_current[current] = skip;
    }
  }

  UF_long state_current_count = 0u;
  UF_long state_best_skip_total = 0u;
  UF_long state_max_best_skip = 0u;
  for (UF_long current = 0u; current < solver->n; ++current) {
    const UF_long skip = best_skip_by_current[current];
    if (skip == 0u) {
      continue;
    }
    state_current_count =
      kls_saturating_add_uf_long(state_current_count, 1u);
    state_best_skip_total =
      kls_saturating_add_uf_long(state_best_skip_total, skip);
    if (skip > state_max_best_skip) {
      state_max_best_skip = skip;
    }
  }

  free(best_skip_by_current);
  *state_current_count_out = state_current_count;
  *state_best_skip_total_out = state_best_skip_total;
  *state_max_best_skip_out = state_max_best_skip;
  *member_skip_lower_bound_out = member_skip_lower_bound;
  return 1;
}

static int kls_build_refactor_btf_scalar_run_group_state_plan(
  kls_solver *solver) {
  if (solver == NULL) {
    return 1;
  }
  free(solver->refactor_btf_scalar_run_group_state_ptr);
  free(solver->refactor_btf_scalar_run_group_state_rows);
  free(solver->refactor_btf_scalar_run_group_state_current_ptr);
  free(solver->refactor_btf_scalar_run_group_state_current_members);
  free(solver->refactor_btf_scalar_run_group_state_current_best_member);
  free(solver->refactor_btf_scalar_run_group_state_value_base);
  solver->refactor_btf_scalar_run_group_state_ptr = NULL;
  solver->refactor_btf_scalar_run_group_state_rows = NULL;
  solver->refactor_btf_scalar_run_group_state_current_ptr = NULL;
  solver->refactor_btf_scalar_run_group_state_current_members = NULL;
  solver->refactor_btf_scalar_run_group_state_current_best_member = NULL;
  solver->refactor_btf_scalar_run_group_state_value_base = NULL;
  solver->refactor_btf_scalar_run_group_state_value_base_count = 0u;
  solver->refactor_btf_scalar_run_group_state_rows_total = 0u;
  solver->refactor_btf_scalar_run_group_state_max_rows = 0u;
  solver->refactor_btf_scalar_run_group_state_current_count = 0u;
  solver->refactor_btf_scalar_run_group_state_best_skip_total = 0u;
  solver->refactor_btf_scalar_run_group_state_max_best_skip = 0u;
  solver->refactor_btf_scalar_run_group_state_guard_lower_bound_rows = 0u;
  solver->refactor_btf_scalar_run_group_state_guard_lower_bound_rejected = 0u;
  solver->refactor_btf_scalar_run_group_live_state_groups = 0u;
  solver->refactor_btf_scalar_run_group_live_state_currents = 0u;
  solver->refactor_btf_scalar_run_group_live_state_rows = 0u;
  solver->refactor_btf_scalar_run_group_live_state_unique_rows = 0u;
  solver->refactor_btf_scalar_run_group_live_state_reused_rows = 0u;
  solver->refactor_btf_scalar_run_group_live_state_max_currents = 0u;
  solver->refactor_btf_scalar_run_group_live_state_max_rows = 0u;
  solver->refactor_btf_scalar_run_group_live_state_max_unique_rows = 0u;

  (void)kls_count_refactor_btf_scalar_run_group_live_state_stats(solver);

  if (!kls_refactor_btf_scalar_run_group_state_stats_env_enabled()) {
    return 1;
  }
  if (solver == NULL || !solver->refactor_btf_scalar_run_group_built ||
      solver->symbolic == NULL || solver->symbolic->R == NULL ||
      solver->refactor_btf_scalar_run_group_multi_current_total == 0u ||
      solver->refactor_btf_scalar_run_group_currents == NULL) {
    return 1;
  }
  const UF_long member_total =
    solver->refactor_btf_scalar_run_group_multi_current_total;
  UF_long cheap_state_current_count = 0u;
  UF_long cheap_state_best_skip_total = 0u;
  UF_long cheap_state_max_best_skip = 0u;
  UF_long member_skip_lower_bound = 0u;
  if (kls_count_refactor_btf_scalar_run_group_state_skip_summary(
        solver, &cheap_state_current_count, &cheap_state_best_skip_total,
        &cheap_state_max_best_skip, &member_skip_lower_bound)) {
    solver->refactor_btf_scalar_run_group_state_current_count =
      cheap_state_current_count;
    solver->refactor_btf_scalar_run_group_state_best_skip_total =
      cheap_state_best_skip_total;
    solver->refactor_btf_scalar_run_group_state_max_best_skip =
      cheap_state_max_best_skip;
    solver->refactor_btf_scalar_run_group_state_guard_lower_bound_rows =
      member_skip_lower_bound;
    if (!kls_refactor_btf_scalar_run_group_state_exec_unguarded_env_enabled() &&
        cheap_state_best_skip_total < member_skip_lower_bound) {
      solver
        ->refactor_btf_scalar_run_group_state_guard_lower_bound_rejected = 1u;
      return 1;
    }
  }
  if (member_total > (UF_long)(KLS_MAX_ALLOCATION / sizeof(UF_long)) - 1u ||
      solver->n > (UF_long)(KLS_MAX_ALLOCATION / sizeof(UF_long)) - 1u) {
    return 1;
  }
  UF_long *state_ptr =
    (UF_long *)calloc((size_t)member_total + 1u, sizeof(*state_ptr));
  UF_long *stamp_workspace =
    (UF_long *)calloc((size_t)(solver->n > 0u ? solver->n : 1u),
                      sizeof(*stamp_workspace));
  UF_long *current_ptr =
    (UF_long *)calloc((size_t)solver->n + 1u, sizeof(*current_ptr));
  UF_long *current_members =
    (UF_long *)malloc((size_t)member_total * sizeof(*current_members));
  UF_long *current_best_member =
    (UF_long *)malloc((size_t)solver->n * sizeof(*current_best_member));
  if (state_ptr == NULL || stamp_workspace == NULL ||
      current_ptr == NULL || current_members == NULL ||
      current_best_member == NULL) {
    free(state_ptr);
    free(stamp_workspace);
    free(current_ptr);
    free(current_members);
    free(current_best_member);
    return 1;
  }
  for (UF_long current = 0u; current < solver->n; ++current) {
    current_best_member[current] = KLS_KLU_EMPTY;
  }

  UF_long stamp_value = 0u;
  UF_long rows_total = 0u;
  UF_long max_rows = 0u;
  for (UF_long member = 0u; member < member_total; ++member) {
    UF_long stamp = 0u;
    UF_long rows = 0u;
    if (!kls_refactor_supernode_algorithm5_next_target_stamp(
          stamp_workspace, solver->n, &stamp_value, &stamp) ||
        !kls_collect_btf_scalar_run_group_state_rows(
          solver, member, stamp_workspace, stamp, NULL, 0u, &rows) ||
        rows > UF_long_max - rows_total) {
      free(state_ptr);
      free(stamp_workspace);
      free(current_ptr);
      free(current_members);
      free(current_best_member);
      return 1;
    }
    const UF_long current =
      solver->refactor_btf_scalar_run_group_currents[member];
    if (current >= solver->n ||
        current_ptr[current + 1u] == UF_long_max) {
      free(state_ptr);
      free(stamp_workspace);
      free(current_ptr);
      free(current_members);
      free(current_best_member);
      return 1;
    }
    current_ptr[current + 1u]++;
    state_ptr[member + 1u] = rows;
    rows_total += rows;
    if (rows > max_rows) {
      max_rows = rows;
    }
  }

  UF_long prefix = 0u;
  for (UF_long member = 0u; member < member_total; ++member) {
    const UF_long rows = state_ptr[member + 1u];
    state_ptr[member] = prefix;
    prefix += rows;
    state_ptr[member + 1u] = prefix;
  }
  if (prefix != rows_total ||
      rows_total > (UF_long)(KLS_MAX_ALLOCATION / sizeof(UF_long))) {
    free(state_ptr);
    free(stamp_workspace);
    free(current_ptr);
    free(current_members);
    free(current_best_member);
    return 1;
  }
  for (UF_long current = 0u; current < solver->n; ++current) {
    if (current_ptr[current] > UF_long_max - current_ptr[current + 1u]) {
      free(state_ptr);
      free(stamp_workspace);
      free(current_ptr);
      free(current_members);
      free(current_best_member);
      return 1;
    }
    current_ptr[current + 1u] += current_ptr[current];
  }
  if (current_ptr[solver->n] != member_total) {
    free(state_ptr);
    free(stamp_workspace);
    free(current_ptr);
    free(current_members);
    free(current_best_member);
    return 1;
  }
  memcpy(stamp_workspace, current_ptr,
         (size_t)solver->n * sizeof(*stamp_workspace));
  for (UF_long member = 0u; member < member_total; ++member) {
    const UF_long current =
      solver->refactor_btf_scalar_run_group_currents[member];
    if (current >= solver->n || stamp_workspace[current] >= member_total) {
      free(state_ptr);
      free(stamp_workspace);
      free(current_ptr);
      free(current_members);
      free(current_best_member);
      return 1;
    }
    current_members[stamp_workspace[current]++] = member;
  }
  UF_long state_current_count = 0u;
  UF_long state_best_skip_total = 0u;
  UF_long state_max_best_skip = 0u;
  for (UF_long current = 0u; current < solver->n; ++current) {
    const UF_long begin = current_ptr[current];
    const UF_long end = current_ptr[current + 1u];
    UF_long best_member = KLS_KLU_EMPTY;
    UF_long best_skip = 0u;
    for (UF_long p = begin; p < end; ++p) {
      const UF_long member = current_members[p];
      const UF_long current_up =
        solver->refactor_btf_scalar_run_group_current_up[member];
      const UF_long rows =
        solver->refactor_btf_scalar_run_group_current_rows[member];
      if (rows > UF_long_max - current_up) {
        free(state_ptr);
        free(stamp_workspace);
        free(current_ptr);
        free(current_members);
        free(current_best_member);
        return 1;
      }
      const UF_long skip = current_up + rows;
      if (skip > best_skip) {
        best_skip = skip;
        best_member = member;
      }
    }
    current_best_member[current] = best_member;
    if (best_member != KLS_KLU_EMPTY) {
      state_current_count =
        kls_saturating_add_uf_long(state_current_count, 1u);
      state_best_skip_total =
        kls_saturating_add_uf_long(state_best_skip_total, best_skip);
      if (best_skip > state_max_best_skip) {
        state_max_best_skip = best_skip;
      }
    }
  }
  UF_long *state_rows =
    rows_total > 0u
      ? (UF_long *)malloc((size_t)rows_total * sizeof(*state_rows))
      : NULL;
  if (rows_total > 0u && state_rows == NULL) {
    free(state_ptr);
    free(stamp_workspace);
    free(current_ptr);
    free(current_members);
    free(current_best_member);
    return 1;
  }

  for (UF_long member = 0u; member < member_total; ++member) {
    const UF_long begin = state_ptr[member];
    const UF_long end = state_ptr[member + 1u];
    UF_long stamp = 0u;
    UF_long filled = 0u;
    if (end < begin ||
        !kls_refactor_supernode_algorithm5_next_target_stamp(
          stamp_workspace, solver->n, &stamp_value, &stamp) ||
        !kls_collect_btf_scalar_run_group_state_rows(
          solver, member, stamp_workspace, stamp, state_rows + begin,
          end - begin, &filled) ||
        filled != end - begin) {
      free(state_ptr);
      free(state_rows);
      free(stamp_workspace);
      free(current_ptr);
      free(current_members);
      free(current_best_member);
      return 1;
    }
    if (filled > 1u) {
      qsort(state_rows + begin, (size_t)filled, sizeof(*state_rows),
            kls_compare_uf_long);
    }
  }

  free(stamp_workspace);
  solver->refactor_btf_scalar_run_group_state_ptr = state_ptr;
  solver->refactor_btf_scalar_run_group_state_rows = state_rows;
  solver->refactor_btf_scalar_run_group_state_current_ptr = current_ptr;
  solver->refactor_btf_scalar_run_group_state_current_members =
    current_members;
  solver->refactor_btf_scalar_run_group_state_current_best_member =
    current_best_member;
  solver->refactor_btf_scalar_run_group_state_rows_total = rows_total;
  solver->refactor_btf_scalar_run_group_state_max_rows = max_rows;
  solver->refactor_btf_scalar_run_group_state_current_count =
    state_current_count;
  solver->refactor_btf_scalar_run_group_state_best_skip_total =
    state_best_skip_total;
  solver->refactor_btf_scalar_run_group_state_max_best_skip =
    state_max_best_skip;
  return 1;
}

static int kls_refactor_btf_scalar_run_group_state_plan_ready(
  const kls_solver *solver) {
  return solver != NULL && solver->refactor_btf_scalar_run_group_built &&
         solver->refactor_btf_scalar_run_group_multi_current_total > 0u &&
         solver->refactor_btf_scalar_run_group_state_ptr != NULL &&
         solver->refactor_btf_scalar_run_group_state_rows != NULL &&
         solver->refactor_btf_scalar_run_group_state_current_ptr != NULL &&
         solver->refactor_btf_scalar_run_group_state_current_members != NULL &&
         solver->refactor_btf_scalar_run_group_state_current_best_member !=
           NULL &&
         solver->refactor_btf_scalar_run_group_state_rows_total > 0u &&
         solver->refactor_btf_scalar_run_group_state_current_count > 0u;
}

static UF_long kls_refactor_btf_scalar_run_group_for_current(
  const kls_solver *solver,
  UF_long dep0_global,
  UF_long current_global,
  UF_long up,
  UF_long *rows_out) {
  KLS_SET_OPTIONAL_OUTPUT(rows_out, 0u);
  if (solver == NULL || !solver->refactor_btf_scalar_run_group_built ||
      solver->refactor_btf_scalar_run_group_multi_count == 0u ||
      solver->refactor_btf_scalar_run_group_ptr == NULL ||
      solver->refactor_btf_scalar_run_group_dep0 == NULL ||
      solver->refactor_btf_scalar_run_group_currents == NULL ||
      solver->refactor_btf_scalar_run_group_current_up == NULL ||
      solver->refactor_btf_scalar_run_group_current_rows == NULL) {
    return KLS_KLU_EMPTY;
  }
  const UF_long group_count =
    solver->refactor_btf_scalar_run_group_multi_count;
  UF_long lo = 0u;
  UF_long hi = group_count;
  while (lo < hi) {
    const UF_long mid = lo + (hi - lo) / 2u;
    if (solver->refactor_btf_scalar_run_group_dep0[mid] < dep0_global) {
      lo = mid + 1u;
    } else {
      hi = mid;
    }
  }
  for (UF_long group = lo; group < group_count &&
       solver->refactor_btf_scalar_run_group_dep0[group] == dep0_global;
       ++group) {
    const UF_long begin = solver->refactor_btf_scalar_run_group_ptr[group];
    const UF_long end =
      solver->refactor_btf_scalar_run_group_ptr[group + 1u];
    for (UF_long pos = begin; pos < end; ++pos) {
      if (solver->refactor_btf_scalar_run_group_currents[pos] ==
            current_global &&
          solver->refactor_btf_scalar_run_group_current_up[pos] == up) {
        if (rows_out != NULL) {
          *rows_out =
            solver->refactor_btf_scalar_run_group_current_rows[pos];
        }
        return group;
      }
    }
  }
  return KLS_KLU_EMPTY;
}

static int kls_refactor_btf_scalar_run_group_range_for_dep0(
  const kls_solver *solver,
  UF_long dep0_global,
  UF_long *begin_out,
  UF_long *end_out) {
  KLS_SET_OPTIONAL_OUTPUT(begin_out, 0u);
  KLS_SET_OPTIONAL_OUTPUT(end_out, 0u);
  if (solver == NULL || !solver->refactor_btf_scalar_run_group_built ||
      solver->refactor_btf_scalar_run_group_multi_count == 0u ||
      solver->refactor_btf_scalar_run_group_dep0 == NULL) {
    return 0;
  }
  const UF_long group_count =
    solver->refactor_btf_scalar_run_group_multi_count;
  UF_long lo = 0u;
  UF_long hi = group_count;
  while (lo < hi) {
    const UF_long mid = lo + (hi - lo) / 2u;
    if (solver->refactor_btf_scalar_run_group_dep0[mid] < dep0_global) {
      lo = mid + 1u;
    } else {
      hi = mid;
    }
  }
  if (lo >= group_count ||
      solver->refactor_btf_scalar_run_group_dep0[lo] != dep0_global) {
    return 0;
  }
  UF_long end = lo + 1u;
  while (end < group_count &&
         solver->refactor_btf_scalar_run_group_dep0[end] == dep0_global) {
    end++;
  }
  if (begin_out != NULL) {
    *begin_out = lo;
  }
  if (end_out != NULL) {
    *end_out = end;
  }
  return 1;
}

static int kls_refactor_btf_scalar_run_group_wake_range_for_dep(
  const kls_solver *solver,
  UF_long producer_global,
  UF_long *begin_out,
  UF_long *end_out) {
  KLS_SET_OPTIONAL_OUTPUT(begin_out, 0u);
  KLS_SET_OPTIONAL_OUTPUT(end_out, 0u);
  if (solver == NULL || !solver->refactor_btf_scalar_run_group_built ||
      solver->refactor_btf_scalar_run_group_wake_count == 0u ||
      solver->refactor_btf_scalar_run_group_wake_dep == NULL) {
    return 0;
  }
  const UF_long wake_count = solver->refactor_btf_scalar_run_group_wake_count;
  UF_long lo = 0u;
  UF_long hi = wake_count;
  while (lo < hi) {
    const UF_long mid = lo + (hi - lo) / 2u;
    if (solver->refactor_btf_scalar_run_group_wake_dep[mid] <
        producer_global) {
      lo = mid + 1u;
    } else {
      hi = mid;
    }
  }
  if (lo >= wake_count ||
      solver->refactor_btf_scalar_run_group_wake_dep[lo] !=
        producer_global) {
    return 0;
  }
  UF_long end = lo + 1u;
  while (end < wake_count &&
         solver->refactor_btf_scalar_run_group_wake_dep[end] ==
           producer_global) {
    end++;
  }
  if (begin_out != NULL) {
    *begin_out = lo;
  }
  if (end_out != NULL) {
    *end_out = end;
  }
  return 1;
}

static int kls_refactor_btf_scalar_run_group_step_range_for_dep(
  const kls_solver *solver,
  UF_long producer_global,
  UF_long *begin_out,
  UF_long *end_out) {
  KLS_SET_OPTIONAL_OUTPUT(begin_out, 0u);
  KLS_SET_OPTIONAL_OUTPUT(end_out, 0u);
  if (solver == NULL || !solver->refactor_btf_scalar_run_group_built ||
      solver->refactor_btf_scalar_run_group_step_dep_count == 0u ||
      solver->refactor_btf_scalar_run_group_step_dep == NULL ||
      solver->refactor_btf_scalar_run_group_step_ptr == NULL) {
    return 0;
  }
  const UF_long dep_count =
    solver->refactor_btf_scalar_run_group_step_dep_count;
  UF_long lo = 0u;
  UF_long hi = dep_count;
  while (lo < hi) {
    const UF_long mid = lo + (hi - lo) / 2u;
    if (solver->refactor_btf_scalar_run_group_step_dep[mid] <
        producer_global) {
      lo = mid + 1u;
    } else {
      hi = mid;
    }
  }
  if (lo >= dep_count ||
      solver->refactor_btf_scalar_run_group_step_dep[lo] !=
        producer_global) {
    return 0;
  }
  if (begin_out != NULL) {
    *begin_out = solver->refactor_btf_scalar_run_group_step_ptr[lo];
  }
  if (end_out != NULL) {
    *end_out = solver->refactor_btf_scalar_run_group_step_ptr[lo + 1u];
  }
  return 1;
}

static void kls_clear_refactor_btf_scalar_run_group_live_step_plan(
  kls_solver *solver) {
  if (solver == NULL) {
    return;
  }
  free(solver->refactor_btf_scalar_run_group_live_step_ptr);
  free(solver->refactor_btf_scalar_run_group_live_step_row_idx);
  solver->refactor_btf_scalar_run_group_live_step_ptr = NULL;
  solver->refactor_btf_scalar_run_group_live_step_row_idx = NULL;
  solver->refactor_btf_scalar_run_group_live_step_count = 0u;
  solver->refactor_btf_scalar_run_group_live_step_current_total = 0u;
  solver->refactor_btf_scalar_run_group_live_step_rows = 0u;
  solver->refactor_btf_scalar_run_group_live_step_unique_rows = 0u;
  solver->refactor_btf_scalar_run_group_live_step_reused_rows = 0u;
  solver->refactor_btf_scalar_run_group_live_step_max_currents = 0u;
  solver->refactor_btf_scalar_run_group_live_step_max_rows = 0u;
  solver->refactor_btf_scalar_run_group_live_step_max_unique_rows = 0u;
  solver->refactor_btf_scalar_run_group_live_step_stored_rows = 0u;
  solver->refactor_btf_scalar_run_group_live_step_storage_limited = 0u;
}

static int kls_fill_refactor_btf_scalar_run_group_live_step_rows(
  kls_solver *solver,
  const UF_long *step_ptr,
  const UF_long *group_step_ptr,
  const UF_long *group_row_step,
  UF_long *step_row_idx,
  UF_long row_capacity) {
  if (solver == NULL || step_ptr == NULL || group_step_ptr == NULL ||
      group_row_step == NULL || step_row_idx == NULL ||
      row_capacity == 0u ||
      !solver->refactor_btf_scalar_run_group_built ||
      solver->refactor_btf_scalar_run_group_ptr == NULL ||
      solver->refactor_btf_scalar_run_group_max_rows_by_group == NULL ||
      solver->refactor_btf_scalar_run_group_current_rows == NULL) {
    return 0;
  }
  const UF_long group_count =
    solver->refactor_btf_scalar_run_group_multi_count;
  const UF_long member_total =
    solver->refactor_btf_scalar_run_group_multi_current_total;
  const UF_long step_count =
    solver->refactor_btf_scalar_run_group_step_count;
  const UF_long stamp_size = solver->n > 0u ? solver->n : 1u;
  if (group_count == 0u || member_total == 0u || step_count == 0u ||
      stamp_size > (UF_long)(KLS_MAX_ALLOCATION / sizeof(UF_long))) {
    return 0;
  }

  UF_long *member_stamp_workspace =
    (UF_long *)calloc((size_t)stamp_size, sizeof(*member_stamp_workspace));
  UF_long *step_stamp_workspace =
    (UF_long *)calloc((size_t)stamp_size, sizeof(*step_stamp_workspace));
  if (member_stamp_workspace == NULL || step_stamp_workspace == NULL) {
    free(member_stamp_workspace);
    free(step_stamp_workspace);
    return 0;
  }

  UF_long member_stamp_value = 0u;
  UF_long step_stamp_value = 0u;
  int ok = 1;

  for (UF_long group = 0u; group < group_count && ok; ++group) {
    const UF_long member_begin =
      solver->refactor_btf_scalar_run_group_ptr[group];
    const UF_long member_end =
      solver->refactor_btf_scalar_run_group_ptr[group + 1u];
    const UF_long group_rows =
      solver->refactor_btf_scalar_run_group_max_rows_by_group[group];
    if (member_begin > member_end || member_end > member_total ||
        group_rows == 0u) {
      ok = 0;
      break;
    }
    const UF_long member_count = member_end - member_begin;
    if (member_count <= 1u ||
        member_count >
          (UF_long)(KLS_MAX_ALLOCATION / sizeof(UF_long)) - 1u) {
      continue;
    }

    UF_long *member_ptr =
      (UF_long *)calloc((size_t)member_count + 1u, sizeof(*member_ptr));
    if (member_ptr == NULL) {
      ok = 0;
      break;
    }
    UF_long member_row_total = 0u;
    for (UF_long local = 0u; local < member_count; ++local) {
      const UF_long member = member_begin + local;
      UF_long member_stamp = 0u;
      UF_long member_rows = 0u;
      if (!kls_refactor_supernode_algorithm5_next_target_stamp(
            member_stamp_workspace, stamp_size, &member_stamp_value,
            &member_stamp) ||
          !kls_collect_btf_scalar_run_group_state_rows(
            solver, member, member_stamp_workspace, member_stamp, NULL, 0u,
            &member_rows) ||
          member_rows > UF_long_max - member_row_total) {
        ok = 0;
        break;
      }
      member_ptr[local] = member_row_total;
      member_row_total += member_rows;
      member_ptr[local + 1u] = member_row_total;
    }
    if (!ok ||
        member_row_total >
          (UF_long)(KLS_MAX_ALLOCATION / sizeof(UF_long))) {
      free(member_ptr);
      ok = 0;
      break;
    }

    UF_long *member_rows =
      member_row_total > 0u
        ? (UF_long *)malloc((size_t)member_row_total * sizeof(*member_rows))
        : NULL;
    if (member_row_total > 0u && member_rows == NULL) {
      free(member_ptr);
      ok = 0;
      break;
    }
    for (UF_long local = 0u; local < member_count; ++local) {
      const UF_long member = member_begin + local;
      const UF_long begin = member_ptr[local];
      const UF_long end = member_ptr[local + 1u];
      UF_long *out_rows = end > begin ? member_rows + begin : NULL;
      UF_long member_stamp = 0u;
      UF_long filled = 0u;
      if (!kls_refactor_supernode_algorithm5_next_target_stamp(
            member_stamp_workspace, stamp_size, &member_stamp_value,
            &member_stamp) ||
          !kls_collect_btf_scalar_run_group_state_rows(
            solver, member, member_stamp_workspace, member_stamp,
            out_rows, end - begin, &filled) ||
          filled != end - begin) {
        ok = 0;
        break;
      }
    }

    for (UF_long row = 0u; row < group_rows && ok; ++row) {
      UF_long active_currents = 0u;
      for (UF_long local = 0u; local < member_count; ++local) {
        const UF_long member = member_begin + local;
        if (solver->refactor_btf_scalar_run_group_current_rows[member] >
            row) {
          active_currents++;
        }
      }
      if (active_currents <= 1u) {
        continue;
      }
      const UF_long step = group_row_step[group_step_ptr[group] + row];
      if (step == KLS_KLU_EMPTY || step >= step_count ||
          step_ptr[step] > step_ptr[step + 1u] ||
          step_ptr[step + 1u] > row_capacity) {
        ok = 0;
        break;
      }
      const UF_long begin = step_ptr[step];
      const UF_long end = step_ptr[step + 1u];
      UF_long unique_rows = 0u;
      UF_long step_stamp = 0u;
      if (begin == end ||
          !kls_refactor_supernode_algorithm5_next_target_stamp(
            step_stamp_workspace, stamp_size, &step_stamp_value,
            &step_stamp)) {
        ok = 0;
        break;
      }
      for (UF_long local = 0u; local < member_count; ++local) {
        const UF_long member = member_begin + local;
        if (solver->refactor_btf_scalar_run_group_current_rows[member] <=
            row) {
          continue;
        }
        const UF_long member_row_begin = member_ptr[local];
        const UF_long member_row_end = member_ptr[local + 1u];
        for (UF_long p = member_row_begin; p < member_row_end; ++p) {
          if (!kls_mark_unique_state_row(
                member_rows[p], step_stamp_workspace, stamp_size, step_stamp,
                step_row_idx + begin, end - begin, &unique_rows)) {
            ok = 0;
            break;
          }
        }
        if (!ok) {
          break;
        }
      }
      if (unique_rows != end - begin) {
        ok = 0;
        break;
      }
      if (unique_rows > 1u) {
        qsort(step_row_idx + begin, (size_t)unique_rows,
              sizeof(*step_row_idx),
              kls_compare_uf_long);
      }
    }

    free(member_rows);
    free(member_ptr);
  }

  free(member_stamp_workspace);
  free(step_stamp_workspace);
  return ok;
}

static int kls_build_refactor_btf_scalar_run_group_step_index(
  kls_solver *solver) {
  if (solver == NULL) {
    return 1;
  }
  free(solver->refactor_btf_scalar_run_group_step_ptr);
  free(solver->refactor_btf_scalar_run_group_step_dep);
  free(solver->refactor_btf_scalar_run_group_step_groups);
  free(solver->refactor_btf_scalar_run_group_step_member_ptr);
  free(solver->refactor_btf_scalar_run_group_step_members);
  solver->refactor_btf_scalar_run_group_step_ptr = NULL;
  solver->refactor_btf_scalar_run_group_step_dep = NULL;
  solver->refactor_btf_scalar_run_group_step_groups = NULL;
  solver->refactor_btf_scalar_run_group_step_member_ptr = NULL;
  solver->refactor_btf_scalar_run_group_step_members = NULL;
  solver->refactor_btf_scalar_run_group_step_dep_count = 0u;
  solver->refactor_btf_scalar_run_group_step_count = 0u;
  solver->refactor_btf_scalar_run_group_step_member_total = 0u;
  solver->refactor_btf_scalar_run_group_producer_index_count = 0u;
  solver->refactor_btf_scalar_run_group_producer_index_max_steps = 0u;
  kls_clear_refactor_btf_scalar_run_group_live_step_plan(solver);

  if (!solver->refactor_btf_scalar_run_group_built ||
      solver->refactor_btf_scalar_run_group_multi_count == 0u ||
      solver->refactor_btf_scalar_run_group_producer_step_count == 0u ||
      solver->refactor_btf_scalar_run_group_producer_step_active_currents == 0u ||
      solver->refactor_btf_scalar_run_group_producer_step_count == UF_long_max ||
      solver->refactor_btf_scalar_run_group_producer_step_active_currents ==
        UF_long_max ||
      solver->refactor_btf_scalar_run_group_ptr == NULL ||
      solver->refactor_btf_scalar_run_group_dep0 == NULL ||
      solver->refactor_btf_scalar_run_group_max_rows_by_group == NULL ||
      solver->refactor_btf_scalar_run_group_current_rows == NULL) {
    return 1;
  }

  const UF_long step_count =
    solver->refactor_btf_scalar_run_group_producer_step_count;
  const UF_long member_total =
    solver->refactor_btf_scalar_run_group_producer_step_active_currents;
  if (step_count > (UF_long)(KLS_MAX_ALLOCATION /
                             sizeof(kls_btf_scalar_run_group_step_item)) ||
      member_total > (UF_long)(KLS_MAX_ALLOCATION / sizeof(UF_long))) {
    return 1;
  }

  kls_btf_scalar_run_group_step_item *items =
    (kls_btf_scalar_run_group_step_item *)malloc(
      (size_t)step_count * sizeof(*items));
  UF_long *member_scratch =
    (UF_long *)malloc((size_t)member_total * sizeof(*member_scratch));
  if (items == NULL || member_scratch == NULL) {
    free(items);
    free(member_scratch);
    return 0;
  }

  UF_long step = 0u;
  UF_long member_out = 0u;
  for (UF_long group = 0u;
       group < solver->refactor_btf_scalar_run_group_multi_count; ++group) {
    const UF_long dep0 = solver->refactor_btf_scalar_run_group_dep0[group];
    const UF_long rows =
      solver->refactor_btf_scalar_run_group_max_rows_by_group[group];
    const UF_long begin = solver->refactor_btf_scalar_run_group_ptr[group];
    const UF_long end = solver->refactor_btf_scalar_run_group_ptr[group + 1u];
    if (begin > end ||
        end > solver->refactor_btf_scalar_run_group_multi_current_total ||
        rows > solver->n || dep0 > solver->n ||
        rows > solver->n - dep0) {
      free(items);
      free(member_scratch);
      return 0;
    }
    for (UF_long row = 0u; row < rows; ++row) {
      const UF_long member_begin = member_out;
      for (UF_long pos = begin; pos < end; ++pos) {
        if (solver->refactor_btf_scalar_run_group_current_rows[pos] > row) {
          if (member_out >= member_total) {
            free(items);
            free(member_scratch);
            return 0;
          }
          member_scratch[member_out++] = pos;
        }
      }
      const UF_long member_count = member_out - member_begin;
      if (member_count == 0u) {
        continue;
      }
      if (step >= step_count) {
        free(items);
        free(member_scratch);
        return 0;
      }
      items[step].producer_global = dep0 + row;
      items[step].group = group;
      items[step].member_begin = member_begin;
      items[step].member_count = member_count;
      step++;
    }
  }
  if (step != step_count || member_out != member_total) {
    free(items);
    free(member_scratch);
    return 0;
  }

  qsort(items, (size_t)step_count, sizeof(*items),
        kls_compare_btf_scalar_run_group_step_item);

  UF_long dep_count = 0u;
  UF_long max_steps = 0u;
  for (UF_long p = 0u; p < step_count;) {
    UF_long q = p + 1u;
    while (q < step_count &&
           items[q].producer_global == items[p].producer_global) {
      q++;
    }
    const UF_long steps = q - p;
    if (steps > max_steps) {
      max_steps = steps;
    }
    dep_count++;
    p = q;
  }
  if (dep_count == 0u ||
      dep_count > (UF_long)(KLS_MAX_ALLOCATION / sizeof(UF_long)) - 1u ||
      step_count > (UF_long)(KLS_MAX_ALLOCATION / sizeof(UF_long)) - 1u) {
    free(items);
    free(member_scratch);
    return 0;
  }

  UF_long *step_ptr =
    (UF_long *)calloc((size_t)dep_count + 1u, sizeof(*step_ptr));
  UF_long *step_dep =
    (UF_long *)malloc((size_t)dep_count * sizeof(*step_dep));
  UF_long *step_groups =
    (UF_long *)malloc((size_t)step_count * sizeof(*step_groups));
  UF_long *step_member_ptr =
    (UF_long *)calloc((size_t)step_count + 1u, sizeof(*step_member_ptr));
  UF_long *step_members =
    (UF_long *)malloc((size_t)member_total * sizeof(*step_members));
  if (step_ptr == NULL || step_dep == NULL || step_groups == NULL ||
      step_member_ptr == NULL || step_members == NULL) {
    free(step_ptr);
    free(step_dep);
    free(step_groups);
    free(step_member_ptr);
    free(step_members);
    free(items);
    free(member_scratch);
    return 0;
  }

  UF_long dep = 0u;
  UF_long out_step = 0u;
  member_out = 0u;
  step_ptr[0] = 0u;
  step_member_ptr[0] = 0u;
  for (UF_long p = 0u; p < step_count;) {
    UF_long q = p + 1u;
    while (q < step_count &&
           items[q].producer_global == items[p].producer_global) {
      q++;
    }
    step_dep[dep] = items[p].producer_global;
    for (UF_long r = p; r < q; ++r) {
      const UF_long count = items[r].member_count;
      if (member_out > member_total || count > member_total - member_out) {
        free(step_ptr);
        free(step_dep);
        free(step_groups);
        free(step_member_ptr);
        free(step_members);
        free(items);
        free(member_scratch);
        return 0;
      }
      step_groups[out_step] = items[r].group;
      memcpy(step_members + member_out,
             member_scratch + items[r].member_begin,
             (size_t)count * sizeof(*step_members));
      member_out += count;
      out_step++;
      step_member_ptr[out_step] = member_out;
    }
    dep++;
    step_ptr[dep] = out_step;
    p = q;
  }
  free(items);
  free(member_scratch);
  if (dep != dep_count || out_step != step_count ||
      member_out != member_total) {
    free(step_ptr);
    free(step_dep);
    free(step_groups);
    free(step_member_ptr);
    free(step_members);
    return 0;
  }

  solver->refactor_btf_scalar_run_group_step_ptr = step_ptr;
  solver->refactor_btf_scalar_run_group_step_dep = step_dep;
  solver->refactor_btf_scalar_run_group_step_groups = step_groups;
  solver->refactor_btf_scalar_run_group_step_member_ptr = step_member_ptr;
  solver->refactor_btf_scalar_run_group_step_members = step_members;
  solver->refactor_btf_scalar_run_group_step_dep_count = dep_count;
  solver->refactor_btf_scalar_run_group_step_count = step_count;
  solver->refactor_btf_scalar_run_group_step_member_total = member_total;
  solver->refactor_btf_scalar_run_group_producer_index_count = dep_count;
  solver->refactor_btf_scalar_run_group_producer_index_max_steps = max_steps;
  return 1;
}

static int kls_build_refactor_btf_scalar_run_group_live_step_plan(
  kls_solver *solver) {
  if (solver == NULL) {
    return 1;
  }
  const int build_plan =
    kls_refactor_btf_scalar_run_group_live_step_plan_env_enabled();
  const int collect_stats =
    build_plan ||
    kls_refactor_btf_scalar_run_group_live_step_stats_env_enabled();
  kls_clear_refactor_btf_scalar_run_group_live_step_plan(solver);
  if (!collect_stats) {
    return 1;
  }
  if (!solver->refactor_btf_scalar_run_group_built ||
      solver->refactor_btf_scalar_run_group_multi_count == 0u ||
      solver->refactor_btf_scalar_run_group_multi_current_total == 0u ||
      solver->refactor_btf_scalar_run_group_step_count == 0u ||
      solver->refactor_btf_scalar_run_group_step_member_total == 0u ||
      solver->refactor_btf_scalar_run_group_ptr == NULL ||
      solver->refactor_btf_scalar_run_group_dep0 == NULL ||
      solver->refactor_btf_scalar_run_group_max_rows_by_group == NULL ||
      solver->refactor_btf_scalar_run_group_current_rows == NULL ||
      solver->refactor_btf_scalar_run_group_step_ptr == NULL ||
      solver->refactor_btf_scalar_run_group_step_dep == NULL ||
      solver->refactor_btf_scalar_run_group_step_groups == NULL ||
      solver->refactor_btf_scalar_run_group_step_member_ptr == NULL ||
      solver->refactor_btf_scalar_run_group_step_members == NULL ||
      solver->refactor_btf_scalar_run_group_currents == NULL) {
    return 1;
  }
  const UF_long group_count =
    solver->refactor_btf_scalar_run_group_multi_count;
  const UF_long step_count =
    solver->refactor_btf_scalar_run_group_step_count;
  const UF_long stamp_size = solver->n > 0u ? solver->n : 1u;
  if (stamp_size > (UF_long)(KLS_MAX_ALLOCATION / sizeof(UF_long)) ||
      step_count > (UF_long)(KLS_MAX_ALLOCATION / sizeof(UF_long)) - 1u ||
      group_count > (UF_long)(KLS_MAX_ALLOCATION / sizeof(UF_long)) - 1u) {
    return 1;
  }

  UF_long *member_stamp_workspace =
    (UF_long *)calloc((size_t)stamp_size, sizeof(*member_stamp_workspace));
  UF_long *step_stamp_workspace =
    (UF_long *)calloc((size_t)stamp_size, sizeof(*step_stamp_workspace));
  UF_long *step_ptr =
    build_plan
      ? (UF_long *)calloc((size_t)step_count + 1u, sizeof(*step_ptr))
      : NULL;
  UF_long *group_step_ptr =
    build_plan
      ? (UF_long *)calloc((size_t)group_count + 1u, sizeof(*group_step_ptr))
      : NULL;
  UF_long *group_row_step = NULL;
  if (member_stamp_workspace == NULL || step_stamp_workspace == NULL ||
      (build_plan && (step_ptr == NULL || group_step_ptr == NULL))) {
    free(member_stamp_workspace);
    free(step_stamp_workspace);
    free(step_ptr);
    free(group_step_ptr);
    return 1;
  }
  if (build_plan) {
    UF_long prefix = 0u;
    for (UF_long group = 0u; group < group_count; ++group) {
      const UF_long rows =
        solver->refactor_btf_scalar_run_group_max_rows_by_group[group];
      group_step_ptr[group] = prefix;
      if (rows > UF_long_max - prefix) {
        free(member_stamp_workspace);
        free(step_stamp_workspace);
        free(step_ptr);
        free(group_step_ptr);
        return 1;
      }
      prefix += rows;
    }
    group_step_ptr[group_count] = prefix;
    if (prefix != step_count ||
        prefix > (UF_long)(KLS_MAX_ALLOCATION / sizeof(UF_long))) {
      free(member_stamp_workspace);
      free(step_stamp_workspace);
      free(step_ptr);
      free(group_step_ptr);
      return 1;
    }
    group_row_step =
      prefix > 0u
        ? (UF_long *)malloc((size_t)prefix * sizeof(*group_row_step))
        : NULL;
    if (prefix > 0u && group_row_step == NULL) {
      free(member_stamp_workspace);
      free(step_stamp_workspace);
      free(step_ptr);
      free(group_step_ptr);
      return 1;
    }
    for (UF_long p = 0u; p < prefix; ++p) {
      group_row_step[p] = KLS_KLU_EMPTY;
    }
    for (UF_long dep = 0u;
         dep < solver->refactor_btf_scalar_run_group_step_dep_count; ++dep) {
      const UF_long producer =
        solver->refactor_btf_scalar_run_group_step_dep[dep];
      const UF_long begin =
        solver->refactor_btf_scalar_run_group_step_ptr[dep];
      const UF_long end =
        solver->refactor_btf_scalar_run_group_step_ptr[dep + 1u];
      if (begin > end || end > step_count) {
        free(member_stamp_workspace);
        free(step_stamp_workspace);
        free(step_ptr);
        free(group_step_ptr);
        free(group_row_step);
        return 1;
      }
      for (UF_long step = begin; step < end; ++step) {
        const UF_long group =
          solver->refactor_btf_scalar_run_group_step_groups[step];
        if (group >= group_count ||
            producer < solver->refactor_btf_scalar_run_group_dep0[group]) {
          free(member_stamp_workspace);
          free(step_stamp_workspace);
          free(step_ptr);
          free(group_step_ptr);
          free(group_row_step);
          return 1;
        }
        const UF_long row =
          producer - solver->refactor_btf_scalar_run_group_dep0[group];
        if (row >=
            solver->refactor_btf_scalar_run_group_max_rows_by_group[group]) {
          free(member_stamp_workspace);
          free(step_stamp_workspace);
          free(step_ptr);
          free(group_step_ptr);
          free(group_row_step);
          return 1;
        }
        group_row_step[group_step_ptr[group] + row] = step;
      }
    }
  }

  UF_long live_steps = 0u;
  UF_long live_currents = 0u;
  UF_long duplicate_rows_total = 0u;
  UF_long unique_rows_total = 0u;
  UF_long reused_rows_total = 0u;
  UF_long max_currents = 0u;
  UF_long max_rows = 0u;
  UF_long max_unique_rows = 0u;
  UF_long member_stamp_value = 0u;
  UF_long step_stamp_value = 0u;
  UF_long *step_row_idx = NULL;
  UF_long stored_rows = 0u;
  UF_long storage_limited = 0u;

  for (UF_long group = 0u; group < group_count; ++group) {
    const UF_long member_begin =
      solver->refactor_btf_scalar_run_group_ptr[group];
    const UF_long member_end =
      solver->refactor_btf_scalar_run_group_ptr[group + 1u];
    const UF_long group_rows =
      solver->refactor_btf_scalar_run_group_max_rows_by_group[group];
    if (member_begin > member_end ||
        member_end >
          solver->refactor_btf_scalar_run_group_multi_current_total ||
        group_rows == 0u) {
      free(member_stamp_workspace);
      free(step_stamp_workspace);
      free(step_ptr);
      free(group_step_ptr);
      free(group_row_step);
      return 1;
    }
    const UF_long member_count = member_end - member_begin;
    if (member_count <= 1u ||
        member_count >
          (UF_long)(KLS_MAX_ALLOCATION / sizeof(UF_long)) - 1u) {
      continue;
    }
    UF_long *member_ptr =
      (UF_long *)calloc((size_t)member_count + 1u, sizeof(*member_ptr));
    if (member_ptr == NULL) {
      free(member_stamp_workspace);
      free(step_stamp_workspace);
      free(step_ptr);
      free(group_step_ptr);
      free(group_row_step);
      return 1;
    }
    UF_long member_row_total = 0u;
    for (UF_long local = 0u; local < member_count; ++local) {
      const UF_long member = member_begin + local;
      UF_long member_stamp = 0u;
      UF_long member_rows = 0u;
      if (!kls_refactor_supernode_algorithm5_next_target_stamp(
            member_stamp_workspace, stamp_size, &member_stamp_value,
            &member_stamp) ||
          !kls_collect_btf_scalar_run_group_state_rows(
            solver, member, member_stamp_workspace, member_stamp, NULL, 0u,
            &member_rows) ||
          member_rows > UF_long_max - member_row_total) {
        free(member_ptr);
        free(member_stamp_workspace);
        free(step_stamp_workspace);
        free(step_ptr);
        free(group_step_ptr);
        free(group_row_step);
        return 1;
      }
      member_ptr[local] = member_row_total;
      member_row_total += member_rows;
      member_ptr[local + 1u] = member_row_total;
    }
    if (member_row_total >
        (UF_long)(KLS_MAX_ALLOCATION / sizeof(UF_long))) {
      free(member_ptr);
      free(member_stamp_workspace);
      free(step_stamp_workspace);
      free(step_ptr);
      free(group_step_ptr);
      free(group_row_step);
      return 1;
    }
    UF_long *member_rows =
      member_row_total > 0u
        ? (UF_long *)malloc((size_t)member_row_total * sizeof(*member_rows))
        : NULL;
    if (member_row_total > 0u && member_rows == NULL) {
      free(member_ptr);
      free(member_stamp_workspace);
      free(step_stamp_workspace);
      free(step_ptr);
      free(group_step_ptr);
      free(group_row_step);
      return 1;
    }
    for (UF_long local = 0u; local < member_count; ++local) {
      const UF_long member = member_begin + local;
      const UF_long begin = member_ptr[local];
      const UF_long end = member_ptr[local + 1u];
      UF_long *out_rows = end > begin ? member_rows + begin : NULL;
      UF_long member_stamp = 0u;
      UF_long filled = 0u;
      if (!kls_refactor_supernode_algorithm5_next_target_stamp(
            member_stamp_workspace, stamp_size, &member_stamp_value,
            &member_stamp) ||
          !kls_collect_btf_scalar_run_group_state_rows(
            solver, member, member_stamp_workspace, member_stamp,
            out_rows, end - begin, &filled) ||
          filled != end - begin) {
        free(member_rows);
        free(member_ptr);
        free(member_stamp_workspace);
        free(step_stamp_workspace);
        free(step_ptr);
        free(group_step_ptr);
        free(group_row_step);
        return 1;
      }
    }

    for (UF_long row = 0u; row < group_rows; ++row) {
      UF_long active_currents = 0u;
      UF_long duplicate_rows = 0u;
      UF_long unique_rows = 0u;
      UF_long step_stamp = 0u;
      if (!kls_refactor_supernode_algorithm5_next_target_stamp(
            step_stamp_workspace, stamp_size, &step_stamp_value,
            &step_stamp)) {
        free(member_rows);
        free(member_ptr);
        free(member_stamp_workspace);
        free(step_stamp_workspace);
        free(step_ptr);
        free(group_step_ptr);
        free(group_row_step);
        return 1;
      }
      for (UF_long local = 0u; local < member_count; ++local) {
        const UF_long member = member_begin + local;
        if (solver->refactor_btf_scalar_run_group_current_rows[member] <=
            row) {
          continue;
        }
        active_currents++;
        const UF_long begin = member_ptr[local];
        const UF_long end = member_ptr[local + 1u];
        duplicate_rows =
          kls_saturating_add_uf_long(duplicate_rows, end - begin);
        for (UF_long p = begin; p < end; ++p) {
          if (!kls_mark_unique_state_row(
                member_rows[p], step_stamp_workspace, stamp_size, step_stamp,
                NULL, 0u, &unique_rows)) {
            free(member_rows);
            free(member_ptr);
            free(member_stamp_workspace);
            free(step_stamp_workspace);
            free(step_ptr);
            free(group_step_ptr);
            free(group_row_step);
            return 1;
          }
        }
      }
      if (active_currents <= 1u) {
        continue;
      }
      live_steps = kls_saturating_add_uf_long(live_steps, 1u);
      live_currents =
        kls_saturating_add_uf_long(live_currents, active_currents);
      duplicate_rows_total =
        kls_saturating_add_uf_long(duplicate_rows_total, duplicate_rows);
      unique_rows_total =
        kls_saturating_add_uf_long(unique_rows_total, unique_rows);
      if (duplicate_rows > unique_rows) {
        reused_rows_total =
          kls_saturating_add_uf_long(reused_rows_total,
                                     duplicate_rows - unique_rows);
      }
      if (active_currents > max_currents) {
        max_currents = active_currents;
      }
      if (duplicate_rows > max_rows) {
        max_rows = duplicate_rows;
      }
      if (unique_rows > max_unique_rows) {
        max_unique_rows = unique_rows;
      }
      if (build_plan) {
        const UF_long step = group_row_step[group_step_ptr[group] + row];
        if (step == KLS_KLU_EMPTY || step >= step_count) {
          free(member_rows);
          free(member_ptr);
          free(member_stamp_workspace);
          free(step_stamp_workspace);
          free(step_ptr);
          free(group_step_ptr);
          free(group_row_step);
          return 1;
        }
        step_ptr[step + 1u] = unique_rows;
      }
    }

    free(member_rows);
    free(member_ptr);
  }

  if (build_plan) {
    UF_long prefix = 0u;
    for (UF_long step = 0u; step < step_count; ++step) {
      const UF_long rows = step_ptr[step + 1u];
      step_ptr[step] = prefix;
      if (rows > UF_long_max - prefix) {
        free(member_stamp_workspace);
        free(step_stamp_workspace);
        free(step_ptr);
        free(group_step_ptr);
        free(group_row_step);
        return 1;
      }
      prefix += rows;
      step_ptr[step + 1u] = prefix;
    }
    if (prefix != unique_rows_total) {
      free(member_stamp_workspace);
      free(step_stamp_workspace);
      free(step_ptr);
      free(group_step_ptr);
      free(group_row_step);
      return 1;
    }
    if (unique_rows_total > 0u) {
      if (unique_rows_total >
          (UF_long)(KLS_MAX_ALLOCATION / sizeof(*step_row_idx))) {
        storage_limited = 1u;
      } else {
        step_row_idx =
          (UF_long *)malloc((size_t)unique_rows_total *
                            sizeof(*step_row_idx));
        if (step_row_idx == NULL) {
          storage_limited = 1u;
        } else if (
          kls_fill_refactor_btf_scalar_run_group_live_step_rows(
            solver, step_ptr, group_step_ptr, group_row_step, step_row_idx,
            unique_rows_total)) {
          stored_rows = unique_rows_total;
        } else {
          free(step_row_idx);
          step_row_idx = NULL;
          storage_limited = 1u;
        }
      }
    }
  }

  free(member_stamp_workspace);
  free(step_stamp_workspace);
  free(group_step_ptr);
  free(group_row_step);
  solver->refactor_btf_scalar_run_group_live_step_ptr = step_ptr;
  solver->refactor_btf_scalar_run_group_live_step_row_idx = step_row_idx;
  solver->refactor_btf_scalar_run_group_live_step_count = live_steps;
  solver->refactor_btf_scalar_run_group_live_step_current_total =
    live_currents;
  solver->refactor_btf_scalar_run_group_live_step_rows =
    duplicate_rows_total;
  solver->refactor_btf_scalar_run_group_live_step_unique_rows =
    unique_rows_total;
  solver->refactor_btf_scalar_run_group_live_step_reused_rows =
    reused_rows_total;
  solver->refactor_btf_scalar_run_group_live_step_max_currents =
    max_currents;
  solver->refactor_btf_scalar_run_group_live_step_max_rows = max_rows;
  solver->refactor_btf_scalar_run_group_live_step_max_unique_rows =
    max_unique_rows;
  solver->refactor_btf_scalar_run_group_live_step_stored_rows =
    stored_rows;
  solver->refactor_btf_scalar_run_group_live_step_storage_limited =
    storage_limited;
  return 1;
}

static int kls_build_refactor_btf_scalar_run_group_wake_index(
  kls_solver *solver) {
  if (solver == NULL || !solver->refactor_btf_scalar_run_group_built ||
      solver->refactor_btf_scalar_run_group_multi_count == 0u ||
      solver->refactor_btf_scalar_run_group_multi_current_total == 0u ||
      solver->refactor_btf_scalar_run_group_ptr == NULL ||
      solver->refactor_btf_scalar_run_group_dep0 == NULL ||
      solver->refactor_btf_scalar_run_group_current_rows == NULL) {
    return 1;
  }
  const UF_long member_total =
    solver->refactor_btf_scalar_run_group_multi_current_total;
  if (member_total > (UF_long)(KLS_MAX_ALLOCATION /
                               sizeof(kls_btf_scalar_run_group_wake_item))) {
    return 0;
  }
  kls_btf_scalar_run_group_wake_item *items =
    (kls_btf_scalar_run_group_wake_item *)malloc(
      (size_t)member_total * sizeof(*items));
  if (items == NULL) {
    return 0;
  }
  UF_long item_count = 0u;
  for (UF_long group = 0u;
       group < solver->refactor_btf_scalar_run_group_multi_count; ++group) {
    const UF_long dep0 = solver->refactor_btf_scalar_run_group_dep0[group];
    const UF_long begin = solver->refactor_btf_scalar_run_group_ptr[group];
    const UF_long end =
      solver->refactor_btf_scalar_run_group_ptr[group + 1u];
    if (begin > end || end > member_total) {
      free(items);
      return 0;
    }
    for (UF_long pos = begin; pos < end; ++pos) {
      const UF_long rows =
        solver->refactor_btf_scalar_run_group_current_rows[pos];
      if (rows == 0u || rows - 1u > UF_long_max - dep0) {
        free(items);
        return 0;
      }
      items[item_count].run_end_global = dep0 + rows - 1u;
      items[item_count].group = group;
      items[item_count].member = pos;
      item_count++;
    }
  }
  if (item_count != member_total) {
    free(items);
    return 0;
  }
  qsort(items, (size_t)item_count, sizeof(*items),
        kls_compare_btf_scalar_run_group_wake_item);

  UF_long wake_count = 0u;
  for (UF_long p = 0u; p < item_count;) {
    UF_long q = p + 1u;
    while (q < item_count &&
           items[q].run_end_global == items[p].run_end_global) {
      q++;
    }
    wake_count++;
    p = q;
  }
  if (wake_count == 0u ||
      wake_count > (UF_long)(KLS_MAX_ALLOCATION / sizeof(UF_long)) - 1u ||
      item_count > (UF_long)(KLS_MAX_ALLOCATION / sizeof(UF_long))) {
    free(items);
    return 0;
  }
  UF_long *wake_ptr =
    (UF_long *)calloc((size_t)wake_count + 1u, sizeof(*wake_ptr));
  UF_long *wake_dep =
    (UF_long *)malloc((size_t)wake_count * sizeof(*wake_dep));
  UF_long *wake_groups =
    (UF_long *)malloc((size_t)item_count * sizeof(*wake_groups));
  UF_long *wake_members =
    (UF_long *)malloc((size_t)item_count * sizeof(*wake_members));
  if (wake_ptr == NULL || wake_dep == NULL || wake_groups == NULL ||
      wake_members == NULL) {
    free(wake_ptr);
    free(wake_dep);
    free(wake_groups);
    free(wake_members);
    free(items);
    return 0;
  }

  UF_long wake = 0u;
  UF_long out = 0u;
  wake_ptr[0] = 0u;
  for (UF_long p = 0u; p < item_count;) {
    UF_long q = p + 1u;
    while (q < item_count &&
           items[q].run_end_global == items[p].run_end_global) {
      q++;
    }
    wake_dep[wake] = items[p].run_end_global;
    for (UF_long r = p; r < q; ++r) {
      wake_groups[out] = items[r].group;
      wake_members[out] = items[r].member;
      out++;
    }
    wake++;
    wake_ptr[wake] = out;
    p = q;
  }
  free(items);
  if (wake != wake_count || out != item_count) {
    free(wake_ptr);
    free(wake_dep);
    free(wake_groups);
    free(wake_members);
    return 0;
  }
  solver->refactor_btf_scalar_run_group_wake_ptr = wake_ptr;
  solver->refactor_btf_scalar_run_group_wake_dep = wake_dep;
  solver->refactor_btf_scalar_run_group_wake_groups = wake_groups;
  solver->refactor_btf_scalar_run_group_wake_members = wake_members;
  solver->refactor_btf_scalar_run_group_wake_count = wake_count;
  solver->refactor_btf_scalar_run_group_wake_member_total = item_count;
  return 1;
}

static int kls_build_refactor_btf_scalar_run_groups(kls_solver *solver) {
  if (!kls_refactor_btf_scalar_run_groups_env_enabled()) {
    kls_clear_refactor_btf_scalar_run_group_descriptor(solver);
    kls_reset_refactor_btf_scalar_run_group_stats(solver);
    return 1;
  }
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      solver->symbolic->R == NULL || solver->numeric->Ulen == NULL ||
      solver->numeric->Llen == NULL ||
      solver->refactor_lu_pointer_count != solver->n ||
      solver->refactor_u_indices == NULL || solver->refactor_u_values == NULL ||
      solver->n <= 0) {
    return 1;
  }
  if (solver->refactor_btf_scalar_run_group_built) {
    (void)kls_build_refactor_btf_scalar_run_group_live_step_plan(solver);
    const int need_state_plan =
      kls_refactor_btf_scalar_run_group_state_stats_env_enabled() ||
      kls_refactor_btf_scalar_run_group_live_state_stats_env_enabled();
    if (need_state_plan) {
      if (!kls_refactor_btf_scalar_run_group_state_plan_ready(solver)) {
        (void)kls_build_refactor_btf_scalar_run_group_state_plan(solver);
      }
    } else {
      free(solver->refactor_btf_scalar_run_group_state_ptr);
      free(solver->refactor_btf_scalar_run_group_state_rows);
      free(solver->refactor_btf_scalar_run_group_state_current_ptr);
      free(solver->refactor_btf_scalar_run_group_state_current_members);
      free(solver->refactor_btf_scalar_run_group_state_current_best_member);
      free(solver->refactor_btf_scalar_run_group_state_value_base);
      solver->refactor_btf_scalar_run_group_state_ptr = NULL;
      solver->refactor_btf_scalar_run_group_state_rows = NULL;
      solver->refactor_btf_scalar_run_group_state_current_ptr = NULL;
      solver->refactor_btf_scalar_run_group_state_current_members = NULL;
      solver->refactor_btf_scalar_run_group_state_current_best_member = NULL;
      solver->refactor_btf_scalar_run_group_state_value_base = NULL;
      solver->refactor_btf_scalar_run_group_state_value_base_count = 0u;
      solver->refactor_btf_scalar_run_group_state_rows_total = 0u;
      solver->refactor_btf_scalar_run_group_state_max_rows = 0u;
      solver->refactor_btf_scalar_run_group_state_current_count = 0u;
      solver->refactor_btf_scalar_run_group_state_best_skip_total = 0u;
      solver->refactor_btf_scalar_run_group_state_max_best_skip = 0u;
      solver->refactor_btf_scalar_run_group_state_guard_lower_bound_rows = 0u;
      solver
        ->refactor_btf_scalar_run_group_state_guard_lower_bound_rejected = 0u;
      solver->refactor_btf_scalar_run_group_live_state_groups = 0u;
      solver->refactor_btf_scalar_run_group_live_state_currents = 0u;
      solver->refactor_btf_scalar_run_group_live_state_rows = 0u;
      solver->refactor_btf_scalar_run_group_live_state_unique_rows = 0u;
      solver->refactor_btf_scalar_run_group_live_state_reused_rows = 0u;
      solver->refactor_btf_scalar_run_group_live_state_max_currents = 0u;
      solver->refactor_btf_scalar_run_group_live_state_max_rows = 0u;
      solver->refactor_btf_scalar_run_group_live_state_max_unique_rows = 0u;
    }
    return 1;
  }
  kls_clear_refactor_btf_scalar_run_group_descriptor(solver);
  kls_reset_refactor_btf_scalar_run_group_stats(solver);

  UF_long item_count = 0u;
  for (UF_long k = 0u; k < solver->n; ++k) {
    UF_long block = 0u;
    if (solver->symbolic->nblocks != 1u) {
      if (solver->refactor_col_block == NULL) {
        return 1;
      }
      block = solver->refactor_col_block[k];
    }
    if (block >= solver->symbolic->nblocks) {
      continue;
    }
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    if (k < k1 || k >= k2 || k2 <= k1 + 1u) {
      continue;
    }
    const UF_long local_k = k - k1;
    const UF_long ucol_len = solver->numeric->Ulen[k];
    const UF_long *ui = solver->refactor_u_indices[k];
    const int32_t *ui32 =
      solver->refactor_u_indices32 != NULL
        ? solver->refactor_u_indices32[k] : NULL;
    for (UF_long up = 0u; up < ucol_len;) {
      const UF_long rows =
        kls_btf_scalar_run_structural_len(up, ucol_len, ui, ui32, local_k);
      if (rows > 1u) {
        if (item_count == UF_long_max) {
          return 1;
        }
        item_count++;
        up += rows;
      } else {
        up++;
      }
    }
  }
  if (item_count == 0u) {
    solver->refactor_btf_scalar_run_group_built = 1;
    return 1;
  }
  if (item_count > (UF_long)(KLS_MAX_ALLOCATION / sizeof(
                     kls_btf_scalar_run_group_item))) {
    return 1;
  }
  kls_btf_scalar_run_group_item *items =
    (kls_btf_scalar_run_group_item *)malloc(
      (size_t)item_count * sizeof(*items));
  if (items == NULL) {
    return 1;
  }

  UF_long item = 0u;
  for (UF_long k = 0u; k < solver->n; ++k) {
    UF_long block = 0u;
    if (solver->symbolic->nblocks != 1u) {
      block = solver->refactor_col_block[k];
    }
    if (block >= solver->symbolic->nblocks) {
      continue;
    }
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    if (k < k1 || k >= k2 || k2 <= k1 + 1u) {
      continue;
    }
    const UF_long local_k = k - k1;
    const UF_long ucol_len = solver->numeric->Ulen[k];
    const UF_long *ui = solver->refactor_u_indices[k];
    const int32_t *ui32 =
      solver->refactor_u_indices32 != NULL
        ? solver->refactor_u_indices32[k] : NULL;
    for (UF_long up = 0u; up < ucol_len;) {
      const UF_long rows =
        kls_btf_scalar_run_structural_len(up, ucol_len, ui, ui32, local_k);
      if (rows > 1u) {
        if (item >= item_count) {
          free(items);
          return 1;
        }
        const UF_long dep0 =
          ui32 != NULL ? (UF_long)ui32[up] : ui[up];
        const UF_long dep0_global = k1 + dep0;
        items[item].block = block;
        items[item].dep0_global = dep0_global;
        items[item].current_global = k;
        items[item].up = up;
        items[item].rows = rows;
        items[item].entries =
          kls_btf_scalar_run_l_entries(solver, dep0_global, rows);
        item++;
        up += rows;
      } else {
        up++;
      }
    }
  }
  if (item != item_count) {
    free(items);
    return 1;
  }

  qsort(items, (size_t)item_count, sizeof(*items),
        kls_compare_btf_scalar_run_group_item);

  UF_long group_count = 0u;
  UF_long multi_count = 0u;
  UF_long multi_current_total = 0u;
  UF_long rows_total = 0u;
  UF_long reused_rows = 0u;
  UF_long entries_total = 0u;
  UF_long reused_entries = 0u;
  UF_long max_currents = 0u;
  UF_long max_rows = 0u;
  UF_long producer_step_count = 0u;
  UF_long producer_step_multi_count = 0u;
  UF_long producer_step_active_currents = 0u;
  UF_long producer_step_unique_entries = 0u;
  UF_long producer_step_duplicate_entries = 0u;
  UF_long producer_step_reused_entries = 0u;
  UF_long producer_step_max_currents = 0u;
  for (UF_long begin = 0u; begin < item_count;) {
    UF_long end = begin + 1u;
    while (end < item_count &&
           items[end].block == items[begin].block &&
           items[end].dep0_global == items[begin].dep0_global) {
      end++;
    }
    const UF_long currents = end - begin;
    UF_long group_rows = 0u;
    UF_long group_entries = 0u;
    UF_long group_max_rows = 0u;
    for (UF_long p = begin; p < end; ++p) {
      if (items[p].rows > group_max_rows) {
        group_max_rows = items[p].rows;
      }
      if (items[p].rows > UF_long_max - group_rows) {
        group_rows = UF_long_max;
      } else {
        group_rows += items[p].rows;
      }
      if (items[p].entries > UF_long_max - group_entries) {
        group_entries = UF_long_max;
      } else {
        group_entries += items[p].entries;
      }
    }
    const UF_long unique_entries =
      kls_btf_scalar_run_l_entries(
        solver, items[begin].dep0_global, group_max_rows);
    group_count++;
    if (group_rows > UF_long_max - rows_total) {
      rows_total = UF_long_max;
    } else {
      rows_total += group_rows;
    }
    if (group_entries > UF_long_max - entries_total) {
      entries_total = UF_long_max;
    } else {
      entries_total += group_entries;
    }
    if (currents > max_currents) {
      max_currents = currents;
    }
    if (group_max_rows > max_rows) {
      max_rows = group_max_rows;
    }
    if (currents > 1u) {
      multi_count++;
      if (currents > UF_long_max - multi_current_total) {
        multi_current_total = UF_long_max;
      } else {
        multi_current_total += currents;
      }
      if (group_rows > group_max_rows) {
        const UF_long saved_rows = group_rows - group_max_rows;
        if (saved_rows > UF_long_max - reused_rows) {
          reused_rows = UF_long_max;
        } else {
          reused_rows += saved_rows;
        }
      }
      if (group_entries > unique_entries) {
        const UF_long saved_entries = group_entries - unique_entries;
        if (saved_entries > UF_long_max - reused_entries) {
          reused_entries = UF_long_max;
        } else {
          reused_entries += saved_entries;
        }
      }
      for (UF_long row = 0u; row < group_max_rows; ++row) {
        UF_long fanout = 0u;
        for (UF_long p = begin; p < end; ++p) {
          if (items[p].rows > row) {
            fanout++;
          }
        }
        if (fanout == 0u) {
          continue;
        }
        producer_step_count =
          kls_saturating_add_uf_long(producer_step_count, 1u);
        producer_step_active_currents =
          kls_saturating_add_uf_long(producer_step_active_currents, fanout);
        if (fanout > producer_step_max_currents) {
          producer_step_max_currents = fanout;
        }
        if (fanout > 1u) {
          producer_step_multi_count =
            kls_saturating_add_uf_long(producer_step_multi_count, 1u);
        }
        if (items[begin].dep0_global > UF_long_max - row ||
            items[begin].dep0_global + row >= solver->n) {
          producer_step_unique_entries = UF_long_max;
          producer_step_duplicate_entries = UF_long_max;
          producer_step_reused_entries = UF_long_max;
          continue;
        }
        const UF_long dep_entries =
          solver->numeric->Llen[items[begin].dep0_global + row];
        const UF_long duplicate_entries =
          kls_saturating_mul_uf_long(dep_entries, fanout);
        producer_step_unique_entries =
          kls_saturating_add_uf_long(producer_step_unique_entries,
                                     dep_entries);
        producer_step_duplicate_entries =
          kls_saturating_add_uf_long(producer_step_duplicate_entries,
                                     duplicate_entries);
        if (fanout > 1u) {
          const UF_long reused =
            kls_saturating_mul_uf_long(dep_entries, fanout - 1u);
          producer_step_reused_entries =
            kls_saturating_add_uf_long(producer_step_reused_entries, reused);
        }
      }
    }
    begin = end;
  }

  solver->refactor_btf_scalar_run_group_count = group_count;
  solver->refactor_btf_scalar_run_group_current_total = item_count;
  solver->refactor_btf_scalar_run_group_multi_count = multi_count;
  solver->refactor_btf_scalar_run_group_multi_current_total =
    multi_current_total;
  solver->refactor_btf_scalar_run_group_rows = rows_total;
  solver->refactor_btf_scalar_run_group_reused_rows = reused_rows;
  solver->refactor_btf_scalar_run_group_entries = entries_total;
  solver->refactor_btf_scalar_run_group_reused_entries = reused_entries;
  solver->refactor_btf_scalar_run_group_max_currents = max_currents;
  solver->refactor_btf_scalar_run_group_max_rows = max_rows;
  solver->refactor_btf_scalar_run_group_producer_step_count =
    producer_step_count;
  solver->refactor_btf_scalar_run_group_producer_step_multi_count =
    producer_step_multi_count;
  solver->refactor_btf_scalar_run_group_producer_step_active_currents =
    producer_step_active_currents;
  solver->refactor_btf_scalar_run_group_producer_step_unique_entries =
    producer_step_unique_entries;
  solver->refactor_btf_scalar_run_group_producer_step_duplicate_entries =
    producer_step_duplicate_entries;
  solver->refactor_btf_scalar_run_group_producer_step_reused_entries =
    producer_step_reused_entries;
  solver->refactor_btf_scalar_run_group_producer_step_max_currents =
    producer_step_max_currents;

  if (multi_count == 0u) {
    solver->refactor_btf_scalar_run_group_built = 1;
    free(items);
    return 1;
  }
  if (multi_count > (UF_long)(KLS_MAX_ALLOCATION / sizeof(UF_long)) - 1u ||
      multi_current_total > (UF_long)(KLS_MAX_ALLOCATION / sizeof(UF_long))) {
    free(items);
    return 1;
  }
  UF_long *group_ptr =
    (UF_long *)calloc((size_t)multi_count + 1u, sizeof(*group_ptr));
  UF_long *group_dep0 =
    (UF_long *)malloc((size_t)multi_count * sizeof(*group_dep0));
  UF_long *group_max_rows =
    (UF_long *)malloc((size_t)multi_count * sizeof(*group_max_rows));
  UF_long *group_currents =
    (UF_long *)malloc((size_t)multi_current_total * sizeof(*group_currents));
  UF_long *group_current_up =
    (UF_long *)malloc((size_t)multi_current_total * sizeof(*group_current_up));
  UF_long *group_current_rows =
    (UF_long *)malloc((size_t)multi_current_total *
                      sizeof(*group_current_rows));
  if (group_ptr == NULL || group_dep0 == NULL || group_max_rows == NULL ||
      group_currents == NULL || group_current_up == NULL ||
      group_current_rows == NULL) {
    free(group_ptr);
    free(group_dep0);
    free(group_max_rows);
    free(group_currents);
    free(group_current_up);
    free(group_current_rows);
    free(items);
    return 1;
  }

  UF_long group = 0u;
  UF_long current_offset = 0u;
  group_ptr[0] = 0u;
  for (UF_long begin = 0u; begin < item_count;) {
    UF_long end = begin + 1u;
    while (end < item_count &&
           items[end].block == items[begin].block &&
           items[end].dep0_global == items[begin].dep0_global) {
      end++;
    }
    const UF_long currents = end - begin;
    if (currents > 1u) {
      UF_long group_max = 0u;
      group_dep0[group] = items[begin].dep0_global;
      for (UF_long p = begin; p < end; ++p) {
        if (items[p].rows > group_max) {
          group_max = items[p].rows;
        }
        group_currents[current_offset] = items[p].current_global;
        group_current_up[current_offset] = items[p].up;
        group_current_rows[current_offset] = items[p].rows;
        current_offset++;
      }
      group_max_rows[group] = group_max;
      group++;
      group_ptr[group] = current_offset;
    }
    begin = end;
  }
  if (group != multi_count || current_offset != multi_current_total) {
    free(group_ptr);
    free(group_dep0);
    free(group_max_rows);
    free(group_currents);
    free(group_current_up);
    free(group_current_rows);
    free(items);
    return 1;
  }

  solver->refactor_btf_scalar_run_group_ptr = group_ptr;
  solver->refactor_btf_scalar_run_group_dep0 = group_dep0;
  solver->refactor_btf_scalar_run_group_max_rows_by_group = group_max_rows;
  solver->refactor_btf_scalar_run_group_currents = group_currents;
  solver->refactor_btf_scalar_run_group_current_up = group_current_up;
  solver->refactor_btf_scalar_run_group_current_rows = group_current_rows;
  solver->refactor_btf_scalar_run_group_built = 1;
  (void)kls_build_refactor_btf_scalar_run_group_step_index(solver);
  (void)kls_build_refactor_btf_scalar_run_group_live_step_plan(solver);
  if (kls_refactor_btf_scalar_run_group_state_stats_env_enabled() ||
      kls_refactor_btf_scalar_run_group_live_state_stats_env_enabled()) {
    (void)kls_build_refactor_btf_scalar_run_group_state_plan(solver);
  }
  (void)kls_build_refactor_btf_scalar_run_group_wake_index(solver);
  free(items);
  return 1;
}

static int kls_fp32_refactor_env_state(void) {
  const char *value = getenv("KLS_ENABLE_FP32_REFACTOR");
  if (value == NULL || value[0] == '\0') {
    return 0; /* auto */
  }
  return (value[0] == '0' && value[1] == '\0') ? -1 : 1;
}

static int kls_fp32_tail_float_accum_enabled(void) {
  static int cached = -1;
  if (cached < 0) {
    const char *value = getenv("KLS_FP32_TAIL_FLOAT_ACCUM");
    cached = value == NULL || value[0] == '\0'
      ? 1
      : !(value[0] == '0' && value[1] == '\0');
  }
  return cached;
}

/* Float-mirrored values remain available as a hardware experiment, but AUTO
   no longer predicts them from a cold refactor/solve ratio.  The current
   EGraph kernels make the conversion-heavy FP32 pass slower on the former
   positive classes as well as on generic controls, and refinement then adds
   another solve.  An explicit request preserves the implementation for
   platforms on which a complete paired-cycle audit proves otherwise. */
static int kls_fp32_refactor_wanted(kls_solver *solver) {
  /* A failed numerical validation is authoritative even when the feature was
     requested through the experiment switch.  Otherwise every later call
     repeats a reduced-precision pass followed by the full double redo (seen
     as a stable 2x regression on ASIC_320k/ASIC_100ks). */
  if (solver != NULL && solver->fp32_decision < 0) {
    return 0;
  }
  const int env = kls_fp32_refactor_env_state();
  if (env <= 0) {
    if (solver != NULL) {
      solver->fp32_decision = -1;
    }
    if (env == 0 && getenv("KLS_TRACE_REFINE") != NULL) {
      fprintf(stderr, "KLS fp32 AUTO: disabled; use"
                      " KLS_ENABLE_FP32_REFACTOR=1 to experiment\n");
    }
    return 0;
  }
  if (solver != NULL && solver->fp32_decision == 0) {
    /* Positive means explicitly requested, not automatically predicted. */
    solver->fp32_decision = 1;
  }
  return 1;
}

/* Single-precision mirror of the L column values.  The refactorization
   tail updates stream producer columns tens of times each, so halving the
   value bytes halves the dominant memory traffic; accumulation stays in
   double and iterative refinement recovers the lost factor accuracy. */
static int kls_ensure_refactor_l_values32_cache(kls_solver *solver) {
  if (solver == NULL || solver->numeric == NULL || solver->symbolic == NULL ||
      solver->numeric->Llen == NULL ||
      solver->refactor_lu_pointer_count != solver->n ||
      solver->refactor_l_indices == NULL || solver->n <= 0) {
    if (solver != NULL && getenv("KLS_TRACE_REFINE") != NULL) {
      fprintf(stderr,
              "KLS fp32 cache: decline llen=%d ptrcount=%ld/%ld lidx=%d\n",
              solver->numeric == NULL || solver->numeric->Llen == NULL,
              (long)solver->refactor_lu_pointer_count, (long)solver->n,
              solver->refactor_l_indices == NULL);
    }
    return 0;
  }
  if (solver->refactor_l_values32 != NULL) {
    return 1;
  }
  UF_long entries = 0;
  for (UF_long col = 0; col < solver->n; ++col) {
    if (solver->refactor_l_indices[col] == NULL) {
      continue;
    }
    entries += solver->numeric->Llen[col];
  }
  float **values32 =
    (float **)calloc((size_t)solver->n, sizeof(*values32));
  float *storage = entries > 0
    ? (float *)malloc((size_t)entries * sizeof(*storage))
    : NULL;
  if (values32 == NULL || (entries > 0 && storage == NULL)) {
    free(values32);
    free(storage);
    return 0;
  }
  UF_long cursor = 0;
  for (UF_long col = 0; col < solver->n; ++col) {
    if (solver->refactor_l_indices[col] == NULL) {
      continue;
    }
    values32[col] = storage + cursor;
    cursor += solver->numeric->Llen[col];
  }
  solver->refactor_l_values32 = values32;
  solver->refactor_l_values32_storage = storage;
  return 1;
}

/* The float mirror is maintained by the store path, but double passes
   reach some columns through paths that do not mirror (measured: any
   fp32 pass directly after a non-fp32 pass reads stale/uninitialized
   floats and fails the validation probe at rmax ~3e-2 on ASIC_320k,
   while fp32-after-fp32 validates at 2e-10).  Refill the whole mirror
   from the current double numeric when engaging the mode after a
   non-fp32 pass; the sweep costs one pass over nnz(L) floats. */
static void kls_refresh_refactor_l_values32_from_numeric(
  kls_solver *solver) {
  if (solver == NULL || solver->refactor_l_values32 == NULL ||
      solver->refactor_l_values == NULL || solver->numeric == NULL ||
      solver->numeric->Llen == NULL) {
    return;
  }
  for (UF_long col = 0; col < solver->n; ++col) {
    float *dst = solver->refactor_l_values32[col];
    const double *src = solver->refactor_l_values[col];
    if (dst == NULL || src == NULL) {
      continue;
    }
    const UF_long len = solver->numeric->Llen[col];
    for (UF_long p = 0; p < len; ++p) {
      dst[p] = (float)src[p];
    }
  }
}

/* A retained pointer cache must describe the numeric the solver holds NOW.
   Numeric rebuilds can replace LU storage wholesale without passing through
   every invalidation site, so the reuse gate re-derives each column pointer
   from the live Lip/Uip/LUbx and rejects the cache on any mismatch. */
static int kls_refactor_lu_pointer_cache_matches(const kls_solver *solver) {
  const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  const trilinos_klu_l_numeric *numeric = solver->numeric;
  for (UF_long block = 0; block < symbolic->nblocks; ++block) {
    const UF_long k1 = symbolic->R[block];
    const UF_long k2 = symbolic->R[block + 1u];
    if (k2 - k1 <= 1u) {
      continue;
    }
    const double *lu = (const double *)numeric->LUbx[block];
    if (lu == NULL) {
      return 0;
    }
    for (UF_long k = k1; k < k2; ++k) {
      const double *lp = lu + numeric->Lip[k];
      const double *up = lu + numeric->Uip[k];
      if (solver->refactor_l_indices[k] != (UF_long *)lp ||
          solver->refactor_u_indices[k] != (UF_long *)up ||
          solver->refactor_l_values[k] !=
            (double *)(lp + kls_klu_units_for_indices(numeric->Llen[k])) ||
          solver->refactor_u_values[k] !=
            (double *)(up + kls_klu_units_for_indices(numeric->Ulen[k]))) {
        return 0;
      }
    }
  }
  return 1;
}

static int kls_ensure_authoritative_packed_l_values(kls_solver *solver) {
  if (solver == NULL || solver->numeric == NULL ||
      solver->numeric->Llen == NULL || solver->refactor_l_values == NULL) {
    return 0;
  }
  if (solver->refactor_l_packed_values != NULL &&
      solver->refactor_l_packed_storage != NULL) {
    return 1;
  }
  UF_long entries = 0u;
  for (UF_long col = 0u; col < solver->n; ++col) {
    if (solver->refactor_l_values[col] != NULL) {
      const UF_long len = solver->numeric->Llen[col];
      if (len > UF_long_max - entries) {
        return 0;
      }
      entries += len;
    }
  }
  if (entries == 0u || entries > (UF_long)(SIZE_MAX / sizeof(double))) {
    return 0;
  }
  double **values =
    (double **)calloc((size_t)solver->n, sizeof(*values));
  double *storage =
    (double *)malloc((size_t)entries * sizeof(*storage));
  if (values == NULL || storage == NULL) {
    free(values);
    free(storage);
    return 0;
  }
  UF_long offset = 0u;
  for (UF_long col = 0u; col < solver->n; ++col) {
    const UF_long len = solver->numeric->Llen[col];
    if (solver->refactor_l_values[col] == NULL || len == 0u) {
      continue;
    }
    values[col] = storage + offset;
    memcpy(values[col], solver->refactor_l_values[col],
           (size_t)len * sizeof(*storage));
    offset += len;
  }
  solver->refactor_l_packed_values = values;
  solver->refactor_l_packed_storage = storage;
  solver->refactor_l_packed_valid = 0;
  return 1;
}

static void kls_sync_authoritative_packed_l_values(kls_solver *solver) {
  if (solver == NULL || !solver->refactor_l_packed_valid ||
      solver->refactor_l_packed_values == NULL ||
      solver->refactor_l_values == NULL || solver->numeric == NULL ||
      solver->numeric->Llen == NULL) {
    return;
  }
  for (UF_long col = 0u; col < solver->n; ++col) {
    const UF_long len = solver->numeric->Llen[col];
    if (solver->refactor_l_values[col] != NULL &&
        solver->refactor_l_packed_values[col] != NULL && len > 0u) {
      memcpy(solver->refactor_l_values[col],
             solver->refactor_l_packed_values[col],
             (size_t)len * sizeof(double));
    }
  }
}

static int kls_ensure_authoritative_packed_u_values(kls_solver *solver) {
  if (solver == NULL || solver->numeric == NULL ||
      solver->numeric->Ulen == NULL || solver->refactor_u_values == NULL) {
    return 0;
  }
  if (solver->refactor_u_packed_values != NULL &&
      solver->refactor_u_packed_storage != NULL) {
    return 1;
  }
  UF_long entries = 0u;
  for (UF_long col = 0u; col < solver->n; ++col) {
    if (solver->refactor_u_values[col] != NULL) {
      const UF_long len = solver->numeric->Ulen[col];
      if (len > UF_long_max - entries) {
        return 0;
      }
      entries += len;
    }
  }
  if (entries == 0u || entries > (UF_long)(SIZE_MAX / sizeof(double))) {
    return 0;
  }
  double **values =
    (double **)calloc((size_t)solver->n, sizeof(*values));
  double *storage =
    (double *)malloc((size_t)entries * sizeof(*storage));
  if (values == NULL || storage == NULL) {
    free(values);
    free(storage);
    return 0;
  }
  UF_long offset = 0u;
  for (UF_long col = 0u; col < solver->n; ++col) {
    const UF_long len = solver->numeric->Ulen[col];
    if (solver->refactor_u_values[col] == NULL || len == 0u) {
      continue;
    }
    values[col] = storage + offset;
    memcpy(values[col], solver->refactor_u_values[col],
           (size_t)len * sizeof(*storage));
    offset += len;
  }
  solver->refactor_u_packed_values = values;
  solver->refactor_u_packed_storage = storage;
  solver->refactor_u_packed_valid = 0;
  return 1;
}

static void kls_sync_authoritative_packed_u_values(kls_solver *solver) {
  if (solver == NULL || !solver->refactor_u_packed_valid ||
      solver->refactor_u_packed_values == NULL ||
      solver->refactor_u_values == NULL || solver->numeric == NULL ||
      solver->numeric->Ulen == NULL) {
    return;
  }
  for (UF_long col = 0u; col < solver->n; ++col) {
    const UF_long len = solver->numeric->Ulen[col];
    if (solver->refactor_u_values[col] != NULL &&
        solver->refactor_u_packed_values[col] != NULL && len > 0u) {
      memcpy(solver->refactor_u_values[col],
             solver->refactor_u_packed_values[col],
             (size_t)len * sizeof(double));
    }
  }
}

static uint64_t kls_execution_domain_llc_bytes(void);

static int kls_authoritative_packed_lu_enabled(const kls_solver *solver) {
  if (solver == NULL) {
    return 0;
  }
  if (!kls_repeated_update_workload(&solver->options) ||
      solver->i32solve_state <= 0 || solver->numeric == NULL) {
    return 0;
  }
  if (solver->pts != NULL && solver->pts_ref_decision >= 0) {
    /* Do not make an EGraph trial authoritative while the independently
       measured subtree refactor tournament is still live or has won. */
    return 0;
  }
  /* This representation removes KLU's interleaved value traffic by keeping
     compact indices and both value decks resident together.  It wins only
     when that active deck fits comfortably in the calling core's LLC; once
     it spills, the duplicate representation increases bandwidth (the
     measured ASIC/G2 regime).  Query the actual execution domain rather
     than recognizing a CPU model or matrix.  Large-cache domains retain the
     already-fast native layout, which also avoids a one-time mirror build. */
#if defined(_SC_LEVEL3_CACHE_SIZE)
  const uint64_t llc = kls_execution_domain_llc_bytes();
  if (llc == 0u || llc > UINT64_C(64) * 1024u * 1024u) {
    return 0;
  }
  const UF_long factor_entries =
    solver->numeric->lnz + solver->numeric->unz;
  if (factor_entries > (UF_long)(SIZE_MAX / sizeof(double)) ||
      factor_entries * (UF_long)sizeof(double) >
        (UF_long)llc * 5u / 8u) {
    return 0;
  }
#else
  return 0;
#endif
  return 1;
}

static int kls_build_refactor_lu_pointer_cache(kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      solver->symbolic->R == NULL || solver->numeric->LUbx == NULL ||
      solver->numeric->Lip == NULL || solver->numeric->Llen == NULL ||
      solver->numeric->Uip == NULL || solver->numeric->Ulen == NULL ||
      solver->n <= 0) {
    return 0;
  }
  const int cache_present =
    solver->refactor_lu_pointer_count == solver->n &&
    solver->refactor_l_indices != NULL &&
    solver->refactor_l_values != NULL &&
    solver->refactor_u_indices != NULL &&
    solver->refactor_u_values != NULL;
  if (cache_present &&
      (solver->refactor_lu_pointer_egraph_certified ||
       kls_refactor_lu_pointer_cache_matches(solver))) {
    /* A successful EGraph pass mutates factor values but never its packed
       layout.  Every representation-changing path frees this cache, which
       also clears the certificate.  The first EGraph call after a cache build
       still performs the complete pointer comparison; subsequent calls avoid
       walking the whole factor merely to rediscover the same addresses. */
    return kls_ensure_refactor_l_index32_cache(solver) &&
           kls_ensure_refactor_u_index32_cache(solver) &&
           kls_ensure_refactor_l_sorted_cache(solver) &&
           kls_build_refactor_btf_scalar_run_groups(solver);
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

  const int record_l_pattern = kls_refactor_l_pattern_stats_env_enabled();
  const int record_l_sorted =
    kls_refactor_btf_scalar_run_exec_enabled(solver) ||
    kls_refactor_btf_scalar_run_group_state_exec_env_enabled();
  const int build_l_index32 =
    kls_refactor_l_index32_env_enabled() &&
    solver->n <= (UF_long)INT32_MAX &&
    1;
  const int build_u_index32 =
    kls_refactor_u_index32_env_enabled() &&
    solver->n <= (UF_long)INT32_MAX &&
    1;
  UF_long l_index32_entries = 0;
  UF_long u_index32_entries = 0;
  UF_long pattern_columns = 0;
  UF_long pattern_entries = 0;
  UF_long adjacent_run_count = 0;
  UF_long adjacent_run_entries = 0;
  UF_long adjacent_run_max_len = 0;
  UF_long contiguous_suffix_columns = 0;
  UF_long contiguous_suffix_entries = 0;
  UF_long contiguous_suffix_max_len = 0;
  int l_indices_sorted = 1;

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
      const UF_long *rows = l_indices[k1 + k];
      if (record_l_sorted && rows != NULL && len > 1u &&
          l_indices_sorted) {
        for (UF_long p = 1u; p < len; ++p) {
          if (rows[p] <= rows[p - 1u]) {
            l_indices_sorted = 0;
            break;
          }
        }
      }
      if (build_l_index32) {
        if (rows != NULL && len > UF_long_max - l_index32_entries) {
          free(l_indices);
          free(l_values);
          free(u_indices);
          free(u_values);
          return 0;
        }
        if (rows != NULL) {
          l_index32_entries += len;
        }
      }
      if (record_l_pattern && rows != NULL && len > 0u) {
        pattern_columns++;
        pattern_entries += len;
        UF_long run_len = 1u;
        for (UF_long p = 1u; p < len; ++p) {
          if (rows[p] == rows[p - 1u] + 1u) {
            run_len++;
          } else {
            if (run_len >= 2u) {
              adjacent_run_count++;
              adjacent_run_entries += run_len;
              if (run_len > adjacent_run_max_len) {
                adjacent_run_max_len = run_len;
              }
            }
            run_len = 1u;
          }
        }
        if (run_len >= 2u) {
          adjacent_run_count++;
          adjacent_run_entries += run_len;
          if (run_len > adjacent_run_max_len) {
            adjacent_run_max_len = run_len;
          }
        }
        UF_long suffix_len = 1u;
        UF_long p = len;
        while (p > 1u && rows[p - 1u] == rows[p - 2u] + 1u) {
          suffix_len++;
          p--;
        }
        if (suffix_len >= 2u) {
          contiguous_suffix_columns++;
          contiguous_suffix_entries += suffix_len;
          if (suffix_len > contiguous_suffix_max_len) {
            contiguous_suffix_max_len = suffix_len;
          }
        }
      }
      kls_klu_get_pointer(lu, uip, ulen, k,
                          &u_indices[k1 + k], &u_values[k1 + k], &len);
      const UF_long *u_rows = u_indices[k1 + k];
      if (build_u_index32) {
        if (u_rows != NULL && len > UF_long_max - u_index32_entries) {
          free(l_indices);
          free(l_values);
          free(u_indices);
          free(u_values);
          return 0;
        }
        if (u_rows != NULL) {
          u_index32_entries += len;
        }
      }
    }
  }

  int32_t **l_indices32 = NULL;
  int32_t *l_indices32_storage = NULL;
  if (build_l_index32) {
    if (solver->n > (UF_long)(SIZE_MAX / sizeof(*l_indices32))) {
      free(l_indices);
      free(l_values);
      free(u_indices);
      free(u_values);
      return 0;
    }
    l_indices32 =
      (int32_t **)calloc((size_t)solver->n, sizeof(*l_indices32));
    if (l_indices32 == NULL) {
      free(l_indices);
      free(l_values);
      free(u_indices);
      free(u_values);
      return 0;
    }
    if (l_index32_entries > 0u) {
      if (l_index32_entries >
          (UF_long)(SIZE_MAX / sizeof(*l_indices32_storage))) {
        free(l_indices32);
        free(l_indices);
        free(l_values);
        free(u_indices);
        free(u_values);
        return 0;
      }
      l_indices32_storage =
        (int32_t *)malloc((size_t)l_index32_entries *
                          sizeof(*l_indices32_storage));
      if (l_indices32_storage == NULL) {
        free(l_indices32);
        free(l_indices);
        free(l_values);
        free(u_indices);
        free(u_values);
        return 0;
      }
    }
    UF_long offset = 0;
    for (UF_long col = 0; col < solver->n; ++col) {
      const UF_long len = solver->numeric->Llen[col];
      const UF_long *rows = l_indices[col];
      if (len == 0u || rows == NULL) {
        l_indices32[col] = NULL;
        continue;
      }
      if (offset > l_index32_entries - len) {
        free(l_indices32_storage);
        free(l_indices32);
        free(l_indices);
        free(l_values);
        free(u_indices);
        free(u_values);
        return 0;
      }
      l_indices32[col] = l_indices32_storage + offset;
      for (UF_long p = 0; p < len; ++p) {
        if (rows[p] > (UF_long)INT32_MAX) {
          free(l_indices32_storage);
          free(l_indices32);
          free(l_indices);
          free(l_values);
          free(u_indices);
          free(u_values);
          return 0;
        }
        l_indices32_storage[offset + p] = (int32_t)rows[p];
      }
      offset += len;
    }
  }

  int32_t **u_indices32 = NULL;
  int32_t *u_indices32_storage = NULL;
  if (build_u_index32) {
    if (solver->n > (UF_long)(SIZE_MAX / sizeof(*u_indices32))) {
      free(l_indices32_storage);
      free(l_indices32);
      free(l_indices);
      free(l_values);
      free(u_indices);
      free(u_values);
      return 0;
    }
    u_indices32 =
      (int32_t **)calloc((size_t)solver->n, sizeof(*u_indices32));
    if (u_indices32 == NULL) {
      free(l_indices32_storage);
      free(l_indices32);
      free(l_indices);
      free(l_values);
      free(u_indices);
      free(u_values);
      return 0;
    }
    if (u_index32_entries > 0u) {
      if (u_index32_entries >
          (UF_long)(SIZE_MAX / sizeof(*u_indices32_storage))) {
        free(u_indices32);
        free(l_indices32_storage);
        free(l_indices32);
        free(l_indices);
        free(l_values);
        free(u_indices);
        free(u_values);
        return 0;
      }
      u_indices32_storage =
        (int32_t *)malloc((size_t)u_index32_entries *
                          sizeof(*u_indices32_storage));
      if (u_indices32_storage == NULL) {
        free(u_indices32);
        free(l_indices32_storage);
        free(l_indices32);
        free(l_indices);
        free(l_values);
        free(u_indices);
        free(u_values);
        return 0;
      }
    }
    UF_long offset = 0;
    for (UF_long col = 0; col < solver->n; ++col) {
      const UF_long len = solver->numeric->Ulen[col];
      const UF_long *rows = u_indices[col];
      if (len == 0u || rows == NULL) {
        u_indices32[col] = NULL;
        continue;
      }
      if (offset > u_index32_entries - len) {
        free(u_indices32_storage);
        free(u_indices32);
        free(l_indices32_storage);
        free(l_indices32);
        free(l_indices);
        free(l_values);
        free(u_indices);
        free(u_values);
        return 0;
      }
      u_indices32[col] = u_indices32_storage + offset;
      for (UF_long p = 0; p < len; ++p) {
        if (rows[p] > (UF_long)INT32_MAX) {
          free(u_indices32_storage);
          free(u_indices32);
          free(l_indices32_storage);
          free(l_indices32);
          free(l_indices);
          free(l_values);
          free(u_indices);
          free(u_values);
          return 0;
        }
        u_indices32_storage[offset + p] = (int32_t)rows[p];
      }
      offset += len;
    }
  }

  solver->refactor_l_indices = l_indices;
  solver->refactor_l_indices32 = l_indices32;
  solver->refactor_l_indices32_storage = l_indices32_storage;
  solver->refactor_l_values = l_values;
  solver->refactor_u_indices = u_indices;
  solver->refactor_u_indices32 = u_indices32;
  solver->refactor_u_indices32_storage = u_indices32_storage;
  solver->refactor_u_values = u_values;
  solver->refactor_lu_pointer_count = solver->n;
  solver->refactor_l_indices_sorted =
    record_l_sorted ? l_indices_sorted : 0;
  solver->refactor_l_indices32_count = l_index32_entries;
  solver->refactor_u_indices32_count = u_index32_entries;
  solver->refactor_l_index32_enabled = l_indices32 != NULL;
  solver->refactor_u_index32_enabled = u_indices32 != NULL;
  solver->refactor_l_pattern_columns = pattern_columns;
  solver->refactor_l_pattern_entries = pattern_entries;
  solver->refactor_l_adjacent_run_count = adjacent_run_count;
  solver->refactor_l_adjacent_run_entries = adjacent_run_entries;
  solver->refactor_l_adjacent_run_max_len = adjacent_run_max_len;
  solver->refactor_l_contiguous_suffix_columns = contiguous_suffix_columns;
  solver->refactor_l_contiguous_suffix_entries = contiguous_suffix_entries;
  solver->refactor_l_contiguous_suffix_max_len = contiguous_suffix_max_len;
  return kls_ensure_refactor_l_sorted_cache(solver) &&
         kls_build_refactor_btf_scalar_run_groups(solver);
}

static int kls_refactor_supernode_panel_block_info(
  const kls_solver *solver,
  UF_long start,
  UF_long end,
  UF_long *k1_out,
  UF_long *block_size_out,
  UF_long *local_start_out) {
  KLS_SET_OPTIONAL_OUTPUT(k1_out, 0u);
  KLS_SET_OPTIONAL_OUTPUT(block_size_out, 0u);
  KLS_SET_OPTIONAL_OUTPUT(local_start_out, 0u);
  if (solver == NULL || solver->symbolic == NULL ||
      solver->symbolic->R == NULL || start >= end || end > solver->n) {
    return 0;
  }
  UF_long k1 = 0u;
  UF_long k2 = solver->n;
  if (solver->symbolic->nblocks > 1u) {
    if (solver->refactor_col_block == NULL) {
      return 0;
    }
    const UF_long block = solver->refactor_col_block[start];
    if (block >= solver->symbolic->nblocks) {
      return 0;
    }
    k1 = solver->symbolic->R[block];
    k2 = solver->symbolic->R[block + 1u];
    if (start < k1 || end > k2) {
      return 0;
    }
    for (UF_long col = start + 1u; col < end; ++col) {
      if (solver->refactor_col_block[col] != block) {
        return 0;
      }
    }
  }
  if (k2 <= k1 || start < k1) {
    return 0;
  }
  if (k1_out != NULL) {
    *k1_out = k1;
  }
  if (block_size_out != NULL) {
    *block_size_out = k2 - k1;
  }
  if (local_start_out != NULL) {
    *local_start_out = start - k1;
  }
  return 1;
}

static int kls_refactor_supernode_panel_candidate(
  const kls_solver *solver,
  UF_long start,
  UF_long end,
  UF_long local_start,
  UF_long block_size,
  UF_long *trailing_len_out,
  UF_long *trailing_rows_out) {
  KLS_SET_OPTIONAL_OUTPUT(trailing_len_out, 0u);
  if (solver == NULL || solver->numeric == NULL ||
      solver->numeric->Llen == NULL ||
      solver->refactor_l_indices == NULL ||
      solver->refactor_l_values == NULL ||
      start >= end || end > solver->n || local_start >= block_size) {
    return 0;
  }
  const UF_long width = end - start;
  if (width <= 1u || width > block_size - local_start) {
    return 0;
  }
  const UF_long skip_begin = local_start;
  const UF_long skip_end = local_start + width;
  const UF_long *first_rows = solver->refactor_l_indices[start];
  const UF_long first_len = solver->numeric->Llen[start];
  if (first_rows == NULL && first_len > 0u) {
    return 0;
  }

  UF_long trailing_len = 0u;
  for (UF_long p = 0; p < first_len; ++p) {
    const UF_long row = first_rows[p];
    if (row >= block_size) {
      return 0;
    }
    if (row < skip_begin || row >= skip_end) {
      if (trailing_rows_out != NULL) {
        trailing_rows_out[trailing_len] = row;
      }
      trailing_len++;
    }
  }

  for (UF_long local = 0; local < width; ++local) {
    const UF_long dep_col = start + local;
    const UF_long *rows = solver->refactor_l_indices[dep_col];
    const UF_long length = solver->numeric->Llen[dep_col];
    if (rows == NULL && length > 0u) {
      return 0;
    }
    const UF_long expected_internal = width - local - 1u;
    UF_long internal_seen = 0u;
    UF_long trailing_seen = 0u;
    UF_long first_pos = 0u;
    for (UF_long p = 0; p < length; ++p) {
      const UF_long row = rows[p];
      if (row >= block_size) {
        return 0;
      }
      if (row >= skip_begin && row < skip_end) {
        const UF_long target = row - skip_begin;
        if (target <= local || target >= width) {
          return 0;
        }
        internal_seen++;
        continue;
      }
      if (local > 0u) {
        while (first_pos < first_len &&
               first_rows[first_pos] >= skip_begin &&
               first_rows[first_pos] < skip_end) {
          first_pos++;
        }
        if (first_pos >= first_len || first_rows[first_pos] != row) {
          return 0;
        }
        first_pos++;
      }
      trailing_seen++;
    }
    if (internal_seen != expected_internal ||
        trailing_seen != trailing_len) {
      return 0;
    }
    if (local > 0u) {
      while (first_pos < first_len) {
        const UF_long row = first_rows[first_pos++];
        if (row < skip_begin || row >= skip_end) {
          return 0;
        }
      }
    }
  }

  if (trailing_len_out != NULL) {
    *trailing_len_out = trailing_len;
  }
  return 1;
}

static void kls_free_refactor_u_supernode_l_cache(kls_solver *solver) {
  if (solver == NULL) {
    return;
  }
  free(solver->refactor_u_supernode_l_valid);
  free(solver->refactor_u_supernode_l_row_begin);
  free(solver->refactor_u_supernode_l_row_ptr);
  free(solver->refactor_u_supernode_l_cols);
  free(solver->refactor_u_supernode_l_dense_values);
  free(solver->refactor_u_supernode_l_values);
  free(solver->refactor_u_supernode_l_used);
  solver->refactor_u_supernode_l_valid = NULL;
  solver->refactor_u_supernode_l_row_begin = NULL;
  solver->refactor_u_supernode_l_row_ptr = NULL;
  solver->refactor_u_supernode_l_cols = NULL;
  solver->refactor_u_supernode_l_dense_values = NULL;
  solver->refactor_u_supernode_l_values = NULL;
  solver->refactor_u_supernode_l_used = NULL;
  solver->refactor_u_supernode_l_built = 0;
  solver->refactor_u_supernode_l_pruned = 0;
  solver->refactor_u_supernode_l_panel_count = 0;
  solver->refactor_u_supernode_l_dense_entries = 0;
  solver->refactor_u_supernode_l_trailing_entries = 0;
  solver->refactor_u_supernode_l_prune_count = 0;
  solver->refactor_u_supernode_l_pruned_panels = 0;
  solver->refactor_u_supernode_l_pruned_dense_entries = 0;
  solver->refactor_u_supernode_l_pruned_trailing_entries = 0;
  solver->refactor_last_u_supernode_l_update_runs = 0;
  solver->refactor_last_u_supernode_l_update_rows = 0;
  solver->refactor_last_u_supernode_l_update_entries = 0;
  solver->refactor_last_u_supernode_l_probe_attempts = 0;
  solver->refactor_last_u_supernode_l_panel_misses = 0;
  solver->refactor_last_u_supernode_l_short_rejects = 0;
  solver->refactor_last_u_supernode_l_stream_rejects = 0;
  solver->refactor_last_u_supernode_l_work_rejects = 0;
  solver->refactor_u_supernode_l_update_run_count = 0;
  solver->refactor_u_supernode_l_update_rows = 0;
  solver->refactor_u_supernode_l_update_entries = 0;
}

static int kls_build_refactor_u_supernode_l_cache(kls_solver *solver) {
  if (solver == NULL || !solver->refactor_u_supernode_pattern_built ||
      solver->numeric == NULL || solver->numeric->Llen == NULL ||
      solver->refactor_l_indices == NULL ||
      solver->refactor_l_values == NULL) {
    return solver != NULL && solver->refactor_u_supernode_l_built;
  }
  if (solver->refactor_u_supernode_l_built) {
    return 1;
  }
  if (solver->refactor_u_supernode_pattern_count == 0u) {
    solver->refactor_u_supernode_l_built = 1;
    return 1;
  }
  if (solver->refactor_u_supernode_start == NULL ||
      solver->refactor_u_supernode_local_start == NULL ||
      solver->refactor_u_supernode_width == NULL ||
      solver->refactor_u_supernode_dense_begin == NULL) {
    return 0;
  }

  kls_free_refactor_u_supernode_l_cache(solver);

  const UF_long pattern_count = solver->refactor_u_supernode_pattern_count;
  const UF_long pattern_rows = solver->refactor_u_supernode_pattern_rows;
  const UF_long dense_entries =
    solver->refactor_u_supernode_value_dense_entries;
  const size_t max_object = (size_t)PTRDIFF_MAX;
  if (pattern_count > (UF_long)(max_object / sizeof(unsigned char)) ||
      pattern_count > (UF_long)(max_object / sizeof(UF_long)) ||
      pattern_rows > (UF_long)((max_object / sizeof(UF_long)) - 1u) ||
      dense_entries > (UF_long)(max_object / sizeof(double))) {
    return 0;
  }

  unsigned char *valid =
    (unsigned char *)calloc((size_t)pattern_count, sizeof(*valid));
  atomic_uchar *used =
    (atomic_uchar *)malloc((size_t)pattern_count * sizeof(*used));
  UF_long *row_begin =
    (UF_long *)malloc((size_t)pattern_count * sizeof(*row_begin));
  UF_long *row_counts =
    (UF_long *)calloc((size_t)(pattern_rows + 1u), sizeof(*row_counts));
  double *dense_values =
    dense_entries > 0u
      ? (double *)calloc((size_t)dense_entries, sizeof(*dense_values))
      : NULL;
  if (valid == NULL || used == NULL || row_begin == NULL ||
      row_counts == NULL ||
      (dense_entries > 0u && dense_values == NULL)) {
    free(valid);
    free(used);
    free(row_begin);
    free(row_counts);
    free(dense_values);
    return 0;
  }
  for (UF_long id = 0u; id < pattern_count; ++id) {
    atomic_init(&used[id], 0u);
  }

  UF_long row_pos = 0u;
  UF_long panel_count = 0u;
  UF_long dense_valid_entries = 0u;
  UF_long trailing_total = 0u;
  for (UF_long id = 0; id < pattern_count; ++id) {
    const UF_long start = solver->refactor_u_supernode_start[id];
    const UF_long width = solver->refactor_u_supernode_width[id];
    const UF_long local_start =
      solver->refactor_u_supernode_local_start[id];
    row_begin[id] = row_pos;
    if (width <= 1u || width > pattern_rows - row_pos ||
        width > UF_long_max / width) {
      free(valid);
      free(used);
      free(row_begin);
      free(row_counts);
      free(dense_values);
      return 0;
    }

    UF_long k1 = 0u;
    UF_long block_size = 0u;
    UF_long checked_local_start = 0u;
    int ok = kls_refactor_supernode_panel_block_info(
      solver, start, start + width, &k1, &block_size,
      &checked_local_start);
    (void)k1;
    if (!ok || checked_local_start != local_start ||
        local_start > block_size || width > block_size - local_start) {
      row_pos += width;
      continue;
    }

    const UF_long skip_begin = local_start;
    const UF_long skip_end = local_start + width;
    UF_long pattern_trailing = 0u;
    ok = 1;
    for (UF_long local = 0u; local < width && ok; ++local) {
      const UF_long col = start + local;
      const UF_long *rows = solver->refactor_l_indices[col];
      const UF_long length = solver->numeric->Llen[col];
      if (rows == NULL && length > 0u) {
        ok = 0;
        break;
      }
      const UF_long expected_internal = width - local - 1u;
      UF_long internal_seen = 0u;
      UF_long trailing_seen = 0u;
      for (UF_long p = 0u; p < length; ++p) {
        const UF_long row = rows[p];
        if (row >= block_size) {
          ok = 0;
          break;
        }
        if (row >= skip_begin && row < skip_end) {
          const UF_long target = row - skip_begin;
          if (target <= local || target >= width) {
            ok = 0;
            break;
          }
          internal_seen++;
        } else {
          trailing_seen++;
        }
      }
      if (internal_seen != expected_internal) {
        ok = 0;
        break;
      }
      if (pattern_trailing > UF_long_max - trailing_total ||
          trailing_seen > UF_long_max - trailing_total - pattern_trailing) {
        ok = 0;
        break;
      }
      row_counts[row_pos + local] = trailing_seen;
      pattern_trailing += trailing_seen;
    }
    if (ok) {
      valid[id] = 1u;
      panel_count++;
      dense_valid_entries += width * width;
      trailing_total += pattern_trailing;
    } else {
      for (UF_long local = 0u; local < width; ++local) {
        row_counts[row_pos + local] = 0u;
      }
    }
    row_pos += width;
  }

  if (row_pos != pattern_rows ||
      trailing_total > (UF_long)(max_object / sizeof(UF_long)) ||
      trailing_total > (UF_long)(max_object / sizeof(double))) {
    free(valid);
    free(used);
    free(row_begin);
    free(row_counts);
    free(dense_values);
    return 0;
  }

  UF_long *row_ptr =
    (UF_long *)malloc((size_t)(pattern_rows + 1u) * sizeof(*row_ptr));
  if (row_ptr == NULL) {
    free(valid);
    free(used);
    free(row_begin);
    free(row_counts);
    free(dense_values);
    return 0;
  }
  row_ptr[0] = 0u;
  for (UF_long row = 0u; row < pattern_rows; ++row) {
    row_ptr[row + 1u] = row_ptr[row] + row_counts[row];
  }
  if (row_ptr[pattern_rows] != trailing_total) {
    free(valid);
    free(used);
    free(row_begin);
    free(row_counts);
    free(row_ptr);
    free(dense_values);
    return 0;
  }

  UF_long *cols =
    trailing_total > 0u
      ? (UF_long *)malloc((size_t)trailing_total * sizeof(*cols))
      : NULL;
  double *values =
    trailing_total > 0u
      ? (double *)calloc((size_t)trailing_total, sizeof(*values))
      : NULL;
  UF_long *cursor =
    pattern_rows > 0u
      ? (UF_long *)malloc((size_t)pattern_rows * sizeof(*cursor))
      : NULL;
  if ((trailing_total > 0u && (cols == NULL || values == NULL)) ||
      (pattern_rows > 0u && cursor == NULL)) {
    free(valid);
    free(used);
    free(row_begin);
    free(row_counts);
    free(row_ptr);
    free(cols);
    free(values);
    free(cursor);
    free(dense_values);
    return 0;
  }
  for (UF_long row = 0u; row < pattern_rows; ++row) {
    cursor[row] = row_ptr[row];
  }

  for (UF_long id = 0; id < pattern_count; ++id) {
    if (!valid[id]) {
      continue;
    }
    const UF_long start = solver->refactor_u_supernode_start[id];
    const UF_long width = solver->refactor_u_supernode_width[id];
    const UF_long local_start =
      solver->refactor_u_supernode_local_start[id];
    const UF_long skip_begin = local_start;
    const UF_long skip_end = local_start + width;
    const UF_long base = row_begin[id];
    for (UF_long local = 0u; local < width; ++local) {
      const UF_long col = start + local;
      const UF_long *rows = solver->refactor_l_indices[col];
      const UF_long length = solver->numeric->Llen[col];
      UF_long dst = cursor[base + local];
      const UF_long dst_end = row_ptr[base + local + 1u];
      for (UF_long p = 0u; p < length; ++p) {
        const UF_long row = rows[p];
        if (row >= skip_begin && row < skip_end) {
          continue;
        }
        if (dst >= dst_end || cols == NULL) {
          free(valid);
          free(used);
          free(row_begin);
          free(row_counts);
          free(row_ptr);
          free(cols);
          free(values);
          free(cursor);
          free(dense_values);
          return 0;
        }
        cols[dst++] = row;
      }
      if (dst != dst_end) {
        free(valid);
        free(used);
        free(row_begin);
        free(row_counts);
        free(row_ptr);
        free(cols);
        free(values);
        free(cursor);
        free(dense_values);
        return 0;
      }
      cursor[base + local] = dst;
    }
  }

  free(row_counts);
  free(cursor);
  solver->refactor_u_supernode_l_valid = valid;
  solver->refactor_u_supernode_l_row_begin = row_begin;
  solver->refactor_u_supernode_l_row_ptr = row_ptr;
  solver->refactor_u_supernode_l_cols = cols;
  solver->refactor_u_supernode_l_dense_values = dense_values;
  solver->refactor_u_supernode_l_values = values;
  solver->refactor_u_supernode_l_used = used;
  solver->refactor_u_supernode_l_built = 1;
  solver->refactor_u_supernode_l_pruned = 0;
  solver->refactor_u_supernode_l_panel_count = panel_count;
  solver->refactor_u_supernode_l_dense_entries = dense_valid_entries;
  solver->refactor_u_supernode_l_trailing_entries = trailing_total;
  return 1;
}

static int kls_egraph_ragged_u_supernode_update_allows(
  UF_long run_rows,
  UF_long dense_scatter_cols,
  UF_long trailing_entries);

static int kls_build_refactor_supernode_consumer_plan_group_l_cache(
  kls_solver *solver,
  int keep_publish_values,
  int keep_batch_candidate_values,
  int keep_advance_candidate_values) {
  if (solver == NULL || !solver->refactor_supernode_consumer_plan_built ||
      solver->numeric == NULL || solver->numeric->Llen == NULL ||
      solver->refactor_l_indices == NULL ||
      solver->refactor_l_values == NULL) {
    return solver != NULL &&
           solver->refactor_supernode_consumer_plan_group_l_built;
  }
  if (solver->refactor_supernode_consumer_plan_group_l_built) {
    return 1;
  }

  kls_free_refactor_supernode_consumer_plan_group_l_cache(solver);

  const UF_long group_count =
    solver->refactor_supernode_consumer_plan_shape_batch_count;
  if (group_count == 0u) {
  solver->refactor_supernode_consumer_plan_group_l_built = 1;
    return 1;
  }
  if (solver->refactor_supernode_consumer_plan_shape_group_panel_start ==
        NULL ||
      solver->refactor_supernode_consumer_plan_shape_group_runs == NULL ||
      solver->refactor_supernode_consumer_plan_shape_group_ptr == NULL ||
      solver->refactor_supernode_consumer_plan_shape_group_panel_offset ==
        NULL ||
      solver->refactor_supernode_consumer_plan_shape_group_rows == NULL) {
    return 0;
  }

  UF_long pattern_rows = 0u;
  for (UF_long group = 0u; group < group_count; ++group) {
    const UF_long width =
      solver->refactor_supernode_consumer_plan_shape_group_rows[group];
    if (width > UF_long_max - pattern_rows) {
      return 0;
    }
    pattern_rows += width;
  }

  const size_t max_object = (size_t)PTRDIFF_MAX;
  if (group_count > (UF_long)(max_object / sizeof(unsigned char)) ||
      group_count > (UF_long)(max_object / sizeof(UF_long)) ||
      solver->refactor_supernode_consumer_plan_run_count >
        (UF_long)(max_object / sizeof(UF_long)) ||
      pattern_rows > (UF_long)((max_object / sizeof(UF_long)) - 1u)) {
    return 0;
  }

  unsigned char *valid =
    (unsigned char *)calloc((size_t)group_count, sizeof(*valid));
  UF_long *dense_begin =
    (UF_long *)malloc((size_t)group_count * sizeof(*dense_begin));
  UF_long *row_begin =
    (UF_long *)malloc((size_t)group_count * sizeof(*row_begin));
  UF_long *row_counts =
    (UF_long *)calloc((size_t)(pattern_rows + 1u), sizeof(*row_counts));
#define KLS_FREE_GROUP_L_BASE() do { \
    free(valid);                         \
    free(dense_begin);                   \
    free(row_begin);                     \
    free(row_counts);                    \
  } while (0)
  if (valid == NULL || dense_begin == NULL || row_begin == NULL ||
      row_counts == NULL) {
    KLS_FREE_GROUP_L_BASE();
    return 0;
  }
  for (UF_long group = 0u; group < group_count; ++group) {
    dense_begin[group] = KLS_KLU_EMPTY;
  }

  UF_long row_pos = 0u;
  UF_long panel_count = 0u;
  UF_long valid_run_count = 0u;
  UF_long valid_rows = 0u;
  UF_long valid_run_rows = 0u;
  UF_long valid_exec_run_count = 0u;
  UF_long valid_exec_run_rows = 0u;
  UF_long batch_candidate_count = 0u;
  UF_long batch_candidate_run_count = 0u;
  UF_long batch_candidate_run_rows = 0u;
  UF_long batch_candidate_entries = 0u;
  UF_long batch_candidate_max_runs = 0u;
  UF_long batch_candidate_advance_count = 0u;
  UF_long batch_candidate_advance_run_count = 0u;
  UF_long batch_candidate_advance_run_rows = 0u;
  double batch_candidate_advance_work = 0.0;
  double batch_candidate_advance_max_work = 0.0;
  UF_long batch_candidate_payoff_count = 0u;
  UF_long batch_candidate_payoff_run_count = 0u;
  UF_long batch_candidate_payoff_run_rows = 0u;
  UF_long batch_candidate_payoff_entries = 0u;
  double batch_candidate_payoff_advance_work = 0.0;
  double batch_candidate_max_payoff_ratio = 0.0;
  UF_long dense_total = 0u;
  UF_long trailing_total = 0u;
  for (UF_long group = 0u; group < group_count; ++group) {
    const UF_long panel_start =
      solver->refactor_supernode_consumer_plan_shape_group_panel_start[group];
    const UF_long panel_offset =
      solver->refactor_supernode_consumer_plan_shape_group_panel_offset[group];
    const UF_long width =
      solver->refactor_supernode_consumer_plan_shape_group_rows[group];
    row_begin[group] = row_pos;
    if (width <= 1u || width > pattern_rows - row_pos ||
        panel_start > solver->n || panel_offset > solver->n - panel_start) {
      row_pos += width;
      continue;
    }
    const UF_long start = panel_start + panel_offset;
    if (start > solver->n || width > solver->n - start ||
        width > UF_long_max / width) {
      row_pos += width;
      continue;
    }

    UF_long k1 = 0u;
    UF_long block_size = 0u;
    UF_long local_start = 0u;
    int ok = kls_refactor_supernode_panel_block_info(
      solver, start, start + width, &k1, &block_size, &local_start);
    (void)k1;
    if (!ok || local_start > block_size ||
        width > block_size - local_start) {
      row_pos += width;
      continue;
    }

    const UF_long skip_begin = local_start;
    const UF_long skip_end = local_start + width;
    UF_long group_trailing = 0u;
    ok = 1;
    for (UF_long local = 0u; local < width && ok; ++local) {
      const UF_long col = start + local;
      const UF_long *rows = solver->refactor_l_indices[col];
      const UF_long length = solver->numeric->Llen[col];
      if (rows == NULL && length > 0u) {
        ok = 0;
        break;
      }
      const UF_long expected_internal = width - local - 1u;
      UF_long internal_seen = 0u;
      UF_long trailing_seen = 0u;
      for (UF_long p = 0u; p < length; ++p) {
        const UF_long row = rows[p];
        if (row >= block_size) {
          ok = 0;
          break;
        }
        if (row >= skip_begin && row < skip_end) {
          const UF_long target = row - skip_begin;
          if (target <= local || target >= width) {
            ok = 0;
            break;
          }
          internal_seen++;
        } else {
          trailing_seen++;
        }
      }
      if (internal_seen != expected_internal ||
          trailing_seen > UF_long_max - group_trailing) {
        ok = 0;
        break;
      }
      row_counts[row_pos + local] = trailing_seen;
      group_trailing += trailing_seen;
    }
    if (ok) {
      const UF_long dense_entries = width * width;
      const UF_long runs_begin =
        solver->refactor_supernode_consumer_plan_shape_group_ptr[group];
      const UF_long runs_end =
        solver->refactor_supernode_consumer_plan_shape_group_ptr[group + 1u];
      if (runs_begin > runs_end ||
          runs_end >
            solver->refactor_supernode_consumer_plan_shape_batch_run_count) {
        KLS_FREE_GROUP_L_BASE();
        return 0;
      }
      const UF_long runs = runs_end - runs_begin;
      if (runs > 0u && width > UF_long_max / runs) {
        KLS_FREE_GROUP_L_BASE();
        return 0;
      }
      const UF_long run_rows = runs * width;
      if (dense_entries > UF_long_max - dense_total ||
          runs > UF_long_max - valid_run_count ||
          width > UF_long_max - valid_rows ||
          run_rows > UF_long_max - valid_run_rows ||
          group_trailing > UF_long_max - trailing_total) {
        KLS_FREE_GROUP_L_BASE();
        return 0;
      }
      valid[group] = KLS_REFACTOR_PLAN_GROUP_L_VALID;
      dense_begin[group] = dense_total;
      panel_count++;
      valid_run_count += runs;
      valid_rows += width;
      valid_run_rows += run_rows;
      if (kls_egraph_ragged_u_supernode_update_allows(
            width, 0u, group_trailing)) {
        if (runs > UF_long_max - valid_exec_run_count ||
            run_rows > UF_long_max - valid_exec_run_rows) {
          KLS_FREE_GROUP_L_BASE();
          return 0;
        }
        valid_exec_run_count += runs;
        valid_exec_run_rows += run_rows;
      }
      if (runs >= 2u) {
        const UF_long trsv_entries = (width * (width - 1u)) / 2u;
        const UF_long per_run_entries = trsv_entries + group_trailing;
        const UF_long group_entries =
          (per_run_entries > 0u && runs > UF_long_max / per_run_entries)
            ? UF_long_max
            : runs * per_run_entries;
        if (batch_candidate_count == UF_long_max ||
            runs > UF_long_max - batch_candidate_run_count ||
            run_rows > UF_long_max - batch_candidate_run_rows ||
            group_entries > UF_long_max - batch_candidate_entries) {
          KLS_FREE_GROUP_L_BASE();
          return 0;
        }
        batch_candidate_count++;
        batch_candidate_run_count += runs;
        batch_candidate_run_rows += run_rows;
        batch_candidate_entries += group_entries;
        if (runs > batch_candidate_max_runs) {
          batch_candidate_max_runs = runs;
        }
        int advance_ok = 1;
        double group_advance_work = 0.0;
        for (UF_long pos = runs_begin; pos < runs_end; ++pos) {
          const UF_long run =
            solver->refactor_supernode_consumer_plan_shape_group_runs[pos];
          UF_long run_deps = 0u;
          double run_work = 0.0;
          if (!kls_supernode_consumer_plan_run_advance_cost(
                solver, solver->refactor_supernode_consumer_plan_current,
                solver->refactor_supernode_consumer_plan_dep,
                solver->refactor_supernode_consumer_plan_run_count, run,
                &run_deps, &run_work)) {
            advance_ok = 0;
            break;
          }
          (void)run_deps;
          group_advance_work += run_work;
        }
        if (advance_ok) {
          valid[group] |= KLS_REFACTOR_PLAN_GROUP_L_ADVANCE;
          if (batch_candidate_advance_count == UF_long_max ||
              runs > UF_long_max - batch_candidate_advance_run_count ||
              run_rows > UF_long_max - batch_candidate_advance_run_rows) {
            KLS_FREE_GROUP_L_BASE();
            return 0;
          }
          batch_candidate_advance_count++;
          batch_candidate_advance_run_count += runs;
          batch_candidate_advance_run_rows += run_rows;
          batch_candidate_advance_work += group_advance_work;
          if (group_advance_work > batch_candidate_advance_max_work) {
            batch_candidate_advance_max_work = group_advance_work;
          }
          const double group_entries_work = (double)group_entries;
          double payoff_ratio = 0.0;
          if (group_advance_work > 0.0) {
            payoff_ratio = group_entries_work / group_advance_work;
          } else if (group_entries > 0u) {
            payoff_ratio = DBL_MAX;
          }
          if (payoff_ratio > batch_candidate_max_payoff_ratio) {
            batch_candidate_max_payoff_ratio = payoff_ratio;
          }
          if (group_entries_work >= group_advance_work) {
            if (batch_candidate_payoff_count == UF_long_max ||
                runs > UF_long_max - batch_candidate_payoff_run_count ||
                run_rows > UF_long_max - batch_candidate_payoff_run_rows ||
                group_entries > UF_long_max - batch_candidate_payoff_entries) {
              KLS_FREE_GROUP_L_BASE();
              return 0;
            }
            batch_candidate_payoff_count++;
            batch_candidate_payoff_run_count += runs;
            batch_candidate_payoff_run_rows += run_rows;
            batch_candidate_payoff_entries += group_entries;
            batch_candidate_payoff_advance_work += group_advance_work;
            valid[group] |= KLS_REFACTOR_PLAN_GROUP_L_PAYOFF;
          }
        }
      }
      dense_total += dense_entries;
      trailing_total += group_trailing;
    } else {
      for (UF_long local = 0u; local < width; ++local) {
        row_counts[row_pos + local] = 0u;
      }
    }
    row_pos += width;
  }

  if (row_pos != pattern_rows ||
      dense_total > (UF_long)(max_object / sizeof(double)) ||
      trailing_total > (UF_long)(max_object / sizeof(UF_long)) ||
      trailing_total > (UF_long)(max_object / sizeof(double))) {
    KLS_FREE_GROUP_L_BASE();
    return 0;
  }

  if (!keep_publish_values &&
      valid_exec_run_count < KLS_REFACTOR_PLAN_GROUP_L_EXEC_MIN_RUNS &&
      (!keep_batch_candidate_values ||
       batch_candidate_payoff_run_count <
        KLS_REFACTOR_PLAN_GROUP_L_BATCH_EXEC_MIN_RUNS) &&
      (!keep_advance_candidate_values ||
       batch_candidate_advance_run_count <
        KLS_REFACTOR_PLAN_GROUP_L_BATCH_EXEC_MIN_RUNS)) {
    KLS_FREE_GROUP_L_BASE();
    solver->refactor_supernode_consumer_plan_group_l_built = 1;
    solver->refactor_supernode_consumer_plan_group_l_storage_limited = 0;
    solver->refactor_supernode_consumer_plan_group_l_panel_count =
      panel_count;
    solver->refactor_supernode_consumer_plan_group_l_run_count =
      valid_run_count;
    solver->refactor_supernode_consumer_plan_group_l_rows = valid_rows;
    solver->refactor_supernode_consumer_plan_group_l_run_rows =
      valid_run_rows;
    solver->refactor_supernode_consumer_plan_group_l_exec_run_count =
      valid_exec_run_count;
    solver->refactor_supernode_consumer_plan_group_l_exec_run_rows =
      valid_exec_run_rows;
    solver->refactor_supernode_consumer_plan_group_l_batch_candidate_count =
      batch_candidate_count;
    solver
      ->refactor_supernode_consumer_plan_group_l_batch_candidate_run_count =
        batch_candidate_run_count;
    solver
      ->refactor_supernode_consumer_plan_group_l_batch_candidate_run_rows =
        batch_candidate_run_rows;
    solver->refactor_supernode_consumer_plan_group_l_batch_candidate_entries =
      batch_candidate_entries;
    solver
      ->refactor_supernode_consumer_plan_group_l_batch_candidate_max_runs =
        batch_candidate_max_runs;
    solver
      ->refactor_supernode_consumer_plan_group_l_batch_candidate_advance_count =
        batch_candidate_advance_count;
    solver
      ->refactor_supernode_consumer_plan_group_l_batch_candidate_advance_run_count =
        batch_candidate_advance_run_count;
    solver
      ->refactor_supernode_consumer_plan_group_l_batch_candidate_advance_run_rows =
        batch_candidate_advance_run_rows;
    solver
      ->refactor_supernode_consumer_plan_group_l_batch_candidate_advance_work =
        batch_candidate_advance_work;
    solver
      ->refactor_supernode_consumer_plan_group_l_batch_candidate_advance_max_work =
        batch_candidate_advance_max_work;
    solver
      ->refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_count =
        batch_candidate_payoff_count;
    solver
      ->refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_run_count =
        batch_candidate_payoff_run_count;
    solver
      ->refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_run_rows =
        batch_candidate_payoff_run_rows;
    solver
      ->refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_entries =
        batch_candidate_payoff_entries;
    solver
      ->refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_advance_work =
        batch_candidate_payoff_advance_work;
    solver
      ->refactor_supernode_consumer_plan_group_l_batch_candidate_max_payoff_ratio =
        batch_candidate_max_payoff_ratio;
    solver->refactor_supernode_consumer_plan_group_l_dense_entries =
      dense_total;
    solver->refactor_supernode_consumer_plan_group_l_trailing_entries =
      trailing_total;
    solver->refactor_supernode_consumer_plan_group_l_bytes = 0u;
    return 1;
  }

  size_t cache_bytes = 0u;
#define KLS_ADD_CACHE_BYTES(count_, type_) \
  do { \
    const size_t kls_count__ = (size_t)(count_); \
    if (kls_count__ > SIZE_MAX / sizeof(type_) || \
        cache_bytes > SIZE_MAX - kls_count__ * sizeof(type_)) { \
      free(valid); \
      free(dense_begin); \
      free(row_begin); \
      free(row_counts); \
      return 0; \
    } \
    cache_bytes += kls_count__ * sizeof(type_); \
  } while (0)
  KLS_ADD_CACHE_BYTES(group_count, unsigned char);
  KLS_ADD_CACHE_BYTES(group_count, UF_long);
  KLS_ADD_CACHE_BYTES(group_count, UF_long);
  KLS_ADD_CACHE_BYTES(pattern_rows + 1u, UF_long);
  KLS_ADD_CACHE_BYTES(trailing_total, UF_long);
  KLS_ADD_CACHE_BYTES(solver->n + 1u, UF_long);
  KLS_ADD_CACHE_BYTES(valid_rows, UF_long);
  KLS_ADD_CACHE_BYTES(valid_rows, UF_long);
  KLS_ADD_CACHE_BYTES(solver->refactor_supernode_consumer_plan_run_count,
                      UF_long);
  KLS_ADD_CACHE_BYTES(dense_total, double);
  KLS_ADD_CACHE_BYTES(trailing_total, double);
#undef KLS_ADD_CACHE_BYTES

  if (cache_bytes > KLS_REFACTOR_PLAN_GROUP_L_CACHE_MAX_BYTES) {
    KLS_FREE_GROUP_L_BASE();
    solver->refactor_supernode_consumer_plan_group_l_built = 1;
    solver->refactor_supernode_consumer_plan_group_l_storage_limited = 1;
    solver->refactor_supernode_consumer_plan_group_l_panel_count =
      panel_count;
    solver->refactor_supernode_consumer_plan_group_l_run_count =
      valid_run_count;
    solver->refactor_supernode_consumer_plan_group_l_rows = valid_rows;
    solver->refactor_supernode_consumer_plan_group_l_run_rows =
      valid_run_rows;
    solver->refactor_supernode_consumer_plan_group_l_exec_run_count =
      valid_exec_run_count;
    solver->refactor_supernode_consumer_plan_group_l_exec_run_rows =
      valid_exec_run_rows;
    solver->refactor_supernode_consumer_plan_group_l_batch_candidate_count =
      batch_candidate_count;
    solver
      ->refactor_supernode_consumer_plan_group_l_batch_candidate_run_count =
        batch_candidate_run_count;
    solver
      ->refactor_supernode_consumer_plan_group_l_batch_candidate_run_rows =
        batch_candidate_run_rows;
    solver->refactor_supernode_consumer_plan_group_l_batch_candidate_entries =
      batch_candidate_entries;
    solver
      ->refactor_supernode_consumer_plan_group_l_batch_candidate_max_runs =
        batch_candidate_max_runs;
    solver
      ->refactor_supernode_consumer_plan_group_l_batch_candidate_advance_count =
        batch_candidate_advance_count;
    solver
      ->refactor_supernode_consumer_plan_group_l_batch_candidate_advance_run_count =
        batch_candidate_advance_run_count;
    solver
      ->refactor_supernode_consumer_plan_group_l_batch_candidate_advance_run_rows =
        batch_candidate_advance_run_rows;
    solver
      ->refactor_supernode_consumer_plan_group_l_batch_candidate_advance_work =
        batch_candidate_advance_work;
    solver
      ->refactor_supernode_consumer_plan_group_l_batch_candidate_advance_max_work =
        batch_candidate_advance_max_work;
    solver
      ->refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_count =
        batch_candidate_payoff_count;
    solver
      ->refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_run_count =
        batch_candidate_payoff_run_count;
    solver
      ->refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_run_rows =
        batch_candidate_payoff_run_rows;
    solver
      ->refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_entries =
        batch_candidate_payoff_entries;
    solver
      ->refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_advance_work =
        batch_candidate_payoff_advance_work;
    solver
      ->refactor_supernode_consumer_plan_group_l_batch_candidate_max_payoff_ratio =
        batch_candidate_max_payoff_ratio;
    solver->refactor_supernode_consumer_plan_group_l_dense_entries =
      dense_total;
    solver->refactor_supernode_consumer_plan_group_l_trailing_entries =
      trailing_total;
    solver->refactor_supernode_consumer_plan_group_l_bytes = cache_bytes;
    return 1;
  }

  UF_long *row_ptr =
    (UF_long *)malloc((size_t)(pattern_rows + 1u) * sizeof(*row_ptr));
  if (row_ptr == NULL) {
    KLS_FREE_GROUP_L_BASE();
    return 0;
  }
  row_ptr[0] = 0u;
  for (UF_long row = 0u; row < pattern_rows; ++row) {
    row_ptr[row + 1u] = row_ptr[row] + row_counts[row];
  }
  if (row_ptr[pattern_rows] != trailing_total) {
    KLS_FREE_GROUP_L_BASE();
    free(row_ptr);
    return 0;
  }

  UF_long *col_ptr =
    (UF_long *)calloc((size_t)(solver->n + 1u), sizeof(*col_ptr));
  if (col_ptr == NULL) {
    KLS_FREE_GROUP_L_BASE();
    free(row_ptr);
    return 0;
  }
  for (UF_long group = 0u; group < group_count; ++group) {
    if ((valid[group] & KLS_REFACTOR_PLAN_GROUP_L_VALID) == 0u) {
      continue;
    }
    const UF_long panel_start =
      solver->refactor_supernode_consumer_plan_shape_group_panel_start[group];
    const UF_long panel_offset =
      solver->refactor_supernode_consumer_plan_shape_group_panel_offset[group];
    const UF_long width =
      solver->refactor_supernode_consumer_plan_shape_group_rows[group];
    const UF_long start = panel_start + panel_offset;
    if (start > solver->n || width > solver->n - start) {
      KLS_FREE_GROUP_L_BASE();
      free(row_ptr);
      free(col_ptr);
      return 0;
    }
    for (UF_long local = 0u; local < width; ++local) {
      const UF_long col = start + local;
      if (col + 1u > solver->n ||
          col_ptr[col + 1u] == UF_long_max) {
        KLS_FREE_GROUP_L_BASE();
        free(row_ptr);
        free(col_ptr);
        return 0;
      }
      col_ptr[col + 1u]++;
    }
  }
  for (UF_long col = 0u; col < solver->n; ++col) {
    if (col_ptr[col + 1u] > UF_long_max - col_ptr[col]) {
      KLS_FREE_GROUP_L_BASE();
      free(row_ptr);
      free(col_ptr);
      return 0;
    }
    col_ptr[col + 1u] += col_ptr[col];
  }
  if (col_ptr[solver->n] != valid_rows) {
    KLS_FREE_GROUP_L_BASE();
    free(row_ptr);
    free(col_ptr);
    return 0;
  }
  UF_long *col_groups =
    valid_rows > 0u
      ? (UF_long *)malloc((size_t)valid_rows * sizeof(*col_groups))
      : NULL;
  UF_long *col_local =
    valid_rows > 0u
      ? (UF_long *)malloc((size_t)valid_rows * sizeof(*col_local))
      : NULL;
  UF_long *col_cursor =
    solver->n > 0u
      ? (UF_long *)malloc((size_t)solver->n * sizeof(*col_cursor))
      : NULL;
  if ((valid_rows > 0u && (col_groups == NULL || col_local == NULL)) ||
      (solver->n > 0u && col_cursor == NULL)) {
    KLS_FREE_GROUP_L_BASE();
    free(row_ptr);
    free(col_ptr);
    free(col_groups);
    free(col_local);
    free(col_cursor);
    return 0;
  }
  if (solver->n > 0u) {
    memcpy(col_cursor, col_ptr, (size_t)solver->n * sizeof(*col_cursor));
  }
  for (UF_long group = 0u; group < group_count; ++group) {
    if ((valid[group] & KLS_REFACTOR_PLAN_GROUP_L_VALID) == 0u) {
      continue;
    }
    const UF_long panel_start =
      solver->refactor_supernode_consumer_plan_shape_group_panel_start[group];
    const UF_long panel_offset =
      solver->refactor_supernode_consumer_plan_shape_group_panel_offset[group];
    const UF_long width =
      solver->refactor_supernode_consumer_plan_shape_group_rows[group];
    const UF_long start = panel_start + panel_offset;
    for (UF_long local = 0u; local < width; ++local) {
      const UF_long col = start + local;
      const UF_long dst = col_cursor[col]++;
      if (dst >= col_ptr[col + 1u] || col_groups == NULL ||
          col_local == NULL) {
        KLS_FREE_GROUP_L_BASE();
        free(row_ptr);
        free(col_ptr);
        free(col_groups);
        free(col_local);
        free(col_cursor);
        return 0;
      }
      col_groups[dst] = group;
      col_local[dst] = local;
    }
  }
  for (UF_long col = 0u; col < solver->n; ++col) {
    if (col_cursor[col] != col_ptr[col + 1u]) {
      KLS_FREE_GROUP_L_BASE();
      free(row_ptr);
      free(col_ptr);
      free(col_groups);
      free(col_local);
      free(col_cursor);
      return 0;
    }
  }
  free(col_cursor);

  UF_long *cols =
    trailing_total > 0u
      ? (UF_long *)malloc((size_t)trailing_total * sizeof(*cols))
      : NULL;
  double *dense_values =
    dense_total > 0u
      ? (double *)calloc((size_t)dense_total, sizeof(*dense_values))
      : NULL;
  double *values =
    trailing_total > 0u
      ? (double *)calloc((size_t)trailing_total, sizeof(*values))
      : NULL;
  UF_long *cursor =
    pattern_rows > 0u
      ? (UF_long *)malloc((size_t)pattern_rows * sizeof(*cursor))
      : NULL;
#define KLS_FREE_GROUP_L_STORAGE() do { \
    free(row_ptr);                          \
    free(col_ptr);                          \
    free(col_groups);                       \
    free(col_local);                        \
    free(cols);                             \
    free(dense_values);                     \
    free(values);                           \
    free(cursor);                           \
    KLS_FREE_GROUP_L_BASE();                \
  } while (0)
  if ((trailing_total > 0u && (cols == NULL || values == NULL)) ||
      (dense_total > 0u && dense_values == NULL) ||
      (pattern_rows > 0u && cursor == NULL)) {
    KLS_FREE_GROUP_L_STORAGE();
    return 0;
  }
  for (UF_long row = 0u; row < pattern_rows; ++row) {
    cursor[row] = row_ptr[row];
  }

  for (UF_long group = 0u; group < group_count; ++group) {
    if ((valid[group] & KLS_REFACTOR_PLAN_GROUP_L_VALID) == 0u) {
      continue;
    }
    const UF_long panel_start =
      solver->refactor_supernode_consumer_plan_shape_group_panel_start[group];
    const UF_long panel_offset =
      solver->refactor_supernode_consumer_plan_shape_group_panel_offset[group];
    const UF_long width =
      solver->refactor_supernode_consumer_plan_shape_group_rows[group];
    const UF_long start = panel_start + panel_offset;
    UF_long k1 = 0u;
    UF_long block_size = 0u;
    UF_long local_start = 0u;
    if (!kls_refactor_supernode_panel_block_info(
          solver, start, start + width, &k1, &block_size, &local_start)) {
      KLS_FREE_GROUP_L_STORAGE();
      return 0;
    }
    (void)k1;
    (void)block_size;
    const UF_long skip_begin = local_start;
    const UF_long skip_end = local_start + width;
    const UF_long base = row_begin[group];
    for (UF_long local = 0u; local < width; ++local) {
      const UF_long col = start + local;
      const UF_long *rows = solver->refactor_l_indices[col];
      const UF_long length = solver->numeric->Llen[col];
      UF_long dst = cursor[base + local];
      const UF_long dst_end = row_ptr[base + local + 1u];
      for (UF_long p = 0u; p < length; ++p) {
        const UF_long row = rows[p];
        if (row >= skip_begin && row < skip_end) {
          continue;
        }
        if (dst >= dst_end || cols == NULL) {
          KLS_FREE_GROUP_L_STORAGE();
          return 0;
        }
        cols[dst++] = row;
      }
      if (dst != dst_end) {
        KLS_FREE_GROUP_L_STORAGE();
        return 0;
      }
      cursor[base + local] = dst;
    }
  }

  const UF_long plan_run_count =
    solver->refactor_supernode_consumer_plan_run_count;
  UF_long *run_group =
    plan_run_count > 0u
      ? (UF_long *)malloc((size_t)plan_run_count * sizeof(*run_group))
      : NULL;
  if (plan_run_count > 0u && run_group == NULL) {
    KLS_FREE_GROUP_L_STORAGE();
    return 0;
  }
  for (UF_long run = 0u; run < plan_run_count; ++run) {
    run_group[run] = KLS_KLU_EMPTY;
  }
  for (UF_long group = 0u; group < group_count; ++group) {
    if ((valid[group] & KLS_REFACTOR_PLAN_GROUP_L_VALID) == 0u) {
      continue;
    }
    const UF_long width =
      solver->refactor_supernode_consumer_plan_shape_group_rows[group];
    const UF_long row_base = row_begin[group];
    if (width <= 1u || row_base + width > pattern_rows ||
        row_ptr[row_base + width] < row_ptr[row_base]) {
      continue;
    }
    const UF_long runs_begin =
      solver->refactor_supernode_consumer_plan_shape_group_ptr[group];
    const UF_long runs_end =
      solver->refactor_supernode_consumer_plan_shape_group_ptr[group + 1u];
    if (runs_begin > runs_end ||
        runs_end >
          solver->refactor_supernode_consumer_plan_shape_batch_run_count) {
      KLS_FREE_GROUP_L_STORAGE();
      free(run_group);
      return 0;
    }
    const UF_long runs = runs_end - runs_begin;
    const UF_long trailing_entries =
      row_ptr[row_base + width] - row_ptr[row_base];
    const int standalone_allowed =
      kls_egraph_ragged_u_supernode_update_allows(
        width, 0u, trailing_entries);
    const int payoff_batch_allowed =
      keep_batch_candidate_values && runs >= 2u &&
      (valid[group] & KLS_REFACTOR_PLAN_GROUP_L_PAYOFF) != 0u;
    const int advance_state_allowed =
      keep_advance_candidate_values && runs >= 2u &&
      (valid[group] & KLS_REFACTOR_PLAN_GROUP_L_ADVANCE) != 0u;
    if (!standalone_allowed && !payoff_batch_allowed &&
        !advance_state_allowed) {
      continue;
    }
    for (UF_long pos = runs_begin; pos < runs_end; ++pos) {
      const UF_long run =
        solver->refactor_supernode_consumer_plan_shape_group_runs[pos];
      if (run >= plan_run_count || run_group == NULL ||
          run_group[run] != KLS_KLU_EMPTY) {
        KLS_FREE_GROUP_L_STORAGE();
        free(run_group);
        return 0;
      }
      run_group[run] = group;
    }
  }

  free(row_counts);
  free(cursor);
  solver->refactor_supernode_consumer_plan_group_l_valid = valid;
  solver->refactor_supernode_consumer_plan_group_l_dense_begin =
    dense_begin;
  solver->refactor_supernode_consumer_plan_group_l_row_begin = row_begin;
  solver->refactor_supernode_consumer_plan_group_l_row_ptr = row_ptr;
  solver->refactor_supernode_consumer_plan_group_l_cols = cols;
  solver->refactor_supernode_consumer_plan_group_l_col_ptr = col_ptr;
  solver->refactor_supernode_consumer_plan_group_l_col_groups = col_groups;
  solver->refactor_supernode_consumer_plan_group_l_col_local = col_local;
  solver->refactor_supernode_consumer_plan_group_l_run_group = run_group;
  solver->refactor_supernode_consumer_plan_group_l_dense_values =
    dense_values;
  solver->refactor_supernode_consumer_plan_group_l_values = values;
#undef KLS_FREE_GROUP_L_STORAGE
#undef KLS_FREE_GROUP_L_BASE
  solver->refactor_supernode_consumer_plan_group_l_built = 1;
  solver->refactor_supernode_consumer_plan_group_l_storage_limited = 0;
  solver->refactor_supernode_consumer_plan_group_l_panel_count =
    panel_count;
  solver->refactor_supernode_consumer_plan_group_l_run_count =
    valid_run_count;
  solver->refactor_supernode_consumer_plan_group_l_rows = valid_rows;
  solver->refactor_supernode_consumer_plan_group_l_run_rows =
    valid_run_rows;
  solver->refactor_supernode_consumer_plan_group_l_exec_run_count =
    valid_exec_run_count;
  solver->refactor_supernode_consumer_plan_group_l_exec_run_rows =
    valid_exec_run_rows;
  solver->refactor_supernode_consumer_plan_group_l_batch_candidate_count =
    batch_candidate_count;
  solver->refactor_supernode_consumer_plan_group_l_batch_candidate_run_count =
    batch_candidate_run_count;
  solver->refactor_supernode_consumer_plan_group_l_batch_candidate_run_rows =
    batch_candidate_run_rows;
  solver->refactor_supernode_consumer_plan_group_l_batch_candidate_entries =
    batch_candidate_entries;
  solver->refactor_supernode_consumer_plan_group_l_batch_candidate_max_runs =
    batch_candidate_max_runs;
  solver
    ->refactor_supernode_consumer_plan_group_l_batch_candidate_advance_count =
      batch_candidate_advance_count;
  solver
    ->refactor_supernode_consumer_plan_group_l_batch_candidate_advance_run_count =
      batch_candidate_advance_run_count;
  solver
    ->refactor_supernode_consumer_plan_group_l_batch_candidate_advance_run_rows =
      batch_candidate_advance_run_rows;
  solver
    ->refactor_supernode_consumer_plan_group_l_batch_candidate_advance_work =
      batch_candidate_advance_work;
  solver
    ->refactor_supernode_consumer_plan_group_l_batch_candidate_advance_max_work =
      batch_candidate_advance_max_work;
  solver
    ->refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_count =
      batch_candidate_payoff_count;
  solver
    ->refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_run_count =
      batch_candidate_payoff_run_count;
  solver
    ->refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_run_rows =
      batch_candidate_payoff_run_rows;
  solver
    ->refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_entries =
      batch_candidate_payoff_entries;
  solver
    ->refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_advance_work =
      batch_candidate_payoff_advance_work;
  solver
    ->refactor_supernode_consumer_plan_group_l_batch_candidate_max_payoff_ratio =
      batch_candidate_max_payoff_ratio;
  solver->refactor_supernode_consumer_plan_group_l_dense_entries =
    dense_total;
  solver->refactor_supernode_consumer_plan_group_l_trailing_entries =
    trailing_total;
  solver->refactor_supernode_consumer_plan_group_l_bytes = cache_bytes;
  return 1;
}

static int
kls_collect_refactor_supernode_consumer_plan_group_l_state_rows(
  const kls_solver *solver,
  UF_long run,
  UF_long *stamp_workspace,
  UF_long stamp_workspace_size,
  UF_long stamp,
  UF_long *rows,
  UF_long row_capacity,
  UF_long *row_count_out) {
  KLS_SET_OPTIONAL_OUTPUT(row_count_out, 0u);
  if (solver == NULL || solver->symbolic == NULL ||
      solver->numeric == NULL || solver->numeric->Llen == NULL ||
      solver->numeric->Ulen == NULL ||
      solver->refactor_supernode_consumer_plan_current == NULL ||
      solver->refactor_supernode_consumer_plan_dep == NULL ||
      solver->refactor_supernode_consumer_plan_rows == NULL ||
      solver->refactor_supernode_consumer_plan_dep_pos == NULL ||
      solver->refactor_supernode_consumer_plan_block_start == NULL ||
      solver->refactor_supernode_consumer_plan_shape_group_panel_start == NULL ||
      solver->refactor_supernode_consumer_plan_shape_group_panel_offset == NULL ||
      solver->refactor_supernode_consumer_plan_shape_group_rows == NULL ||
      solver->refactor_supernode_consumer_plan_group_l_valid == NULL ||
      solver->refactor_supernode_consumer_plan_group_l_run_group == NULL ||
      stamp_workspace == NULL || stamp_workspace_size == 0u ||
      stamp == 0u || row_count_out == NULL ||
      run >= solver->refactor_supernode_consumer_plan_run_count) {
    return 0;
  }

  const UF_long group =
    solver->refactor_supernode_consumer_plan_group_l_run_group[run];
  if (group == KLS_KLU_EMPTY) {
    return 1;
  }
  if (group >= solver->refactor_supernode_consumer_plan_shape_batch_count) {
    return 0;
  }
  if ((solver->refactor_supernode_consumer_plan_group_l_valid[group] &
       KLS_REFACTOR_PLAN_GROUP_L_ADVANCE) == 0u) {
    return 1;
  }

  const UF_long panel_start =
    solver->refactor_supernode_consumer_plan_shape_group_panel_start[group];
  const UF_long panel_offset =
    solver->refactor_supernode_consumer_plan_shape_group_panel_offset[group];
  const UF_long width =
    solver->refactor_supernode_consumer_plan_shape_group_rows[group];
  const UF_long current =
    solver->refactor_supernode_consumer_plan_current[run];
  const UF_long dep =
    solver->refactor_supernode_consumer_plan_dep[run];
  const UF_long run_rows =
    solver->refactor_supernode_consumer_plan_rows[run];
  const UF_long dep_pos =
    solver->refactor_supernode_consumer_plan_dep_pos[run];
  const UF_long plan_block_start =
    solver->refactor_supernode_consumer_plan_block_start[run];
  if (width <= 1u || run_rows != width ||
      panel_start > solver->n || panel_offset > solver->n - panel_start ||
      current >= solver->n || dep >= current || dep_pos == KLS_KLU_EMPTY ||
      plan_block_start == KLS_KLU_EMPTY) {
    return 0;
  }
  const UF_long start = panel_start + panel_offset;
  if (start != dep || start > solver->n || width > solver->n - start) {
    return 0;
  }

  UF_long k1 = 0u;
  UF_long block_size = 0u;
  UF_long local_start = 0u;
  if (!kls_refactor_supernode_panel_block_info(
        solver, start, start + width, &k1, &block_size, &local_start)) {
    return 0;
  }
  if (k1 != plan_block_start || block_size == 0u ||
      block_size > stamp_workspace_size || local_start >= block_size ||
      width > block_size - local_start || current < k1 ||
      current >= k1 + block_size) {
    return 0;
  }
  const UF_long current_local = current - k1;
  if (local_start >= current_local ||
      width > current_local - local_start) {
    return 0;
  }

  UF_long row_count = 0u;
  if (solver->refactor_col_ptr == NULL ||
      (solver->refactor_row_idx == NULL &&
       solver->refactor_row_idx32 == NULL)) {
    return 0;
  }
  const UF_long scatter_col_begin = solver->refactor_col_ptr[current];
  const UF_long scatter_end = solver->refactor_col_ptr[current + 1u];
  if (scatter_col_begin > scatter_end || scatter_end > solver->nnz) {
    return 0;
  }
  UF_long scatter_begin = scatter_col_begin;
  if (solver->symbolic->nblocks != 1u) {
    if (solver->refactor_block_start == NULL ||
        solver->refactor_block_start[current] < scatter_col_begin ||
        solver->refactor_block_start[current] > scatter_end) {
      return 0;
    }
    scatter_begin = solver->refactor_block_start[current];
  }
  if (solver->refactor_row_idx32 != NULL) {
    const int32_t *input_rows32 = solver->refactor_row_idx32;
    for (UF_long p = scatter_begin; p < scatter_end; ++p) {
      const UF_long row = (UF_long)input_rows32[p];
      if (row >= k1 && row < k1 + block_size &&
          !kls_mark_unique_state_row(
            row - k1, stamp_workspace, stamp_workspace_size, stamp, rows,
            row_capacity, &row_count)) {
        return 0;
      }
    }
  } else {
    const UF_long *input_rows = solver->refactor_row_idx;
    for (UF_long p = scatter_begin; p < scatter_end; ++p) {
      const UF_long row = input_rows[p];
      if (row >= k1 && row < k1 + block_size &&
          !kls_mark_unique_state_row(
            row - k1, stamp_workspace, stamp_workspace_size, stamp, rows,
            row_capacity, &row_count)) {
        return 0;
      }
    }
  }

  const UF_long ulen = solver->numeric->Ulen[current];
  const int32_t *ui32 =
    solver->refactor_u_indices32 != NULL
      ? solver->refactor_u_indices32[current] : NULL;
  const UF_long *ui =
    ui32 == NULL && solver->refactor_u_indices != NULL
      ? solver->refactor_u_indices[current] : NULL;
  if (dep_pos > ulen || width > ulen - dep_pos ||
      (dep_pos + width > 0u && ui32 == NULL && ui == NULL)) {
    return 0;
  }
  const UF_long prefix_end = dep_pos + width;
  for (UF_long p = 0u; p < prefix_end; ++p) {
    const UF_long j = ui32 != NULL ? (UF_long)ui32[p] : ui[p];
    if (j >= current_local || j >= block_size) {
      return 0;
    }
    if (p >= dep_pos && j != local_start + (p - dep_pos)) {
      return 0;
    }
    if (!kls_mark_unique_state_row(
          j, stamp_workspace, stamp_workspace_size, stamp, rows,
          row_capacity, &row_count)) {
      return 0;
    }
    const UF_long col = k1 + j;
    if (col >= solver->n) {
      return 0;
    }
    const UF_long llen = solver->numeric->Llen[col];
    const int32_t *li32 =
      solver->refactor_l_indices32 != NULL
        ? solver->refactor_l_indices32[col] : NULL;
    const UF_long *li =
      li32 == NULL && solver->refactor_l_indices != NULL
        ? solver->refactor_l_indices[col] : NULL;
    if (llen > 0u && li32 == NULL && li == NULL) {
      return 0;
    }
    for (UF_long lp = 0u; lp < llen; ++lp) {
      const UF_long row = li32 != NULL ? (UF_long)li32[lp] : li[lp];
      if (row >= block_size ||
          !kls_mark_unique_state_row(
            row, stamp_workspace, stamp_workspace_size, stamp, rows,
            row_capacity, &row_count)) {
        return 0;
      }
    }
  }

  if (!kls_mark_unique_state_row(
        current_local, stamp_workspace, stamp_workspace_size, stamp, rows,
        row_capacity, &row_count)) {
    return 0;
  }
  *row_count_out = row_count;
  return 1;
}

static int kls_build_refactor_supernode_consumer_plan_group_l_state_plan(
  kls_solver *solver,
  int focus_largest_groups) {
  if (solver == NULL) {
    return 0;
  }
  if (solver->refactor_supernode_consumer_plan_group_l_state_built) {
    return 1;
  }

  kls_free_refactor_supernode_consumer_plan_group_l_state_plan(solver);
  if (!solver->refactor_supernode_consumer_plan_built ||
      !solver->refactor_supernode_consumer_plan_group_l_built ||
      solver->symbolic == NULL) {
    return 0;
  }

  const UF_long run_count =
    solver->refactor_supernode_consumer_plan_run_count;
  const UF_long group_count =
    solver->refactor_supernode_consumer_plan_shape_batch_count;
  if (run_count == 0u || group_count == 0u ||
      solver->refactor_supernode_consumer_plan_group_l_valid == NULL ||
      solver->refactor_supernode_consumer_plan_group_l_run_group == NULL) {
    solver->refactor_supernode_consumer_plan_group_l_state_built = 1;
    solver->refactor_supernode_consumer_plan_group_l_state_focus_enabled =
      focus_largest_groups ? 1 : 0;
    return 1;
  }

  const size_t max_object = (size_t)PTRDIFF_MAX;
  if (run_count > (UF_long)((max_object / sizeof(UF_long)) - 1u) ||
      run_count > (UF_long)(max_object / sizeof(UF_long)) ||
      group_count > (UF_long)(max_object / sizeof(UF_long)) ||
      group_count > (UF_long)(max_object / sizeof(unsigned char))) {
    return 0;
  }
  UF_long stamp_size = solver->symbolic->maxblock;
  if (stamp_size == 0u) {
    stamp_size = solver->n;
  }
  if (stamp_size == 0u ||
      stamp_size > (UF_long)(max_object / sizeof(UF_long))) {
    return 0;
  }

  UF_long *state_ptr =
    (UF_long *)malloc((size_t)(run_count + 1u) * sizeof(*state_ptr));
  UF_long *stamp_workspace =
    (UF_long *)calloc((size_t)stamp_size, sizeof(*stamp_workspace));
  UF_long *run_state_rows =
    (UF_long *)calloc((size_t)run_count, sizeof(*run_state_rows));
  UF_long *group_state_rows =
    (UF_long *)calloc((size_t)group_count, sizeof(*group_state_rows));
  UF_long *group_state_runs =
    (UF_long *)calloc((size_t)group_count, sizeof(*group_state_runs));
  UF_long *group_state_max_rows =
    (UF_long *)calloc((size_t)group_count, sizeof(*group_state_max_rows));
  unsigned char *group_selected =
    focus_largest_groups
      ? (unsigned char *)calloc((size_t)group_count, sizeof(*group_selected))
      : NULL;
  UF_long *state_row_idx = NULL;
  if (state_ptr == NULL || stamp_workspace == NULL ||
      run_state_rows == NULL || group_state_rows == NULL ||
      group_state_runs == NULL || group_state_max_rows == NULL ||
      (focus_largest_groups && group_selected == NULL)) {
    goto fail;
  }
  for (UF_long group = 0u; group < group_count; ++group) {
    solver->refactor_supernode_consumer_plan_group_l_valid[group] &=
      (unsigned char)~KLS_REFACTOR_PLAN_GROUP_L_STATE_FOCUS;
  }

  UF_long stamp = 0u;
  UF_long prefix = 0u;
  UF_long state_run_count = 0u;
  UF_long state_max_rows = 0u;
  UF_long state_candidate_group_count = 0u;
  state_ptr[0] = 0u;
  for (UF_long run = 0u; run < run_count; ++run) {
    UF_long run_stamp = 0u;
    UF_long row_count = 0u;
    if (!kls_refactor_supernode_algorithm5_next_target_stamp(
          stamp_workspace, stamp_size, &stamp, &run_stamp) ||
        !kls_collect_refactor_supernode_consumer_plan_group_l_state_rows(
          solver, run, stamp_workspace, stamp_size, run_stamp, NULL, 0u,
          &row_count) ||
        row_count > UF_long_max - prefix) {
      goto fail;
    }
    run_state_rows[run] = row_count;
    if (row_count > 0u) {
      const UF_long group =
        solver->refactor_supernode_consumer_plan_group_l_run_group[run];
      if (group == KLS_KLU_EMPTY || group >= group_count ||
          row_count > UF_long_max - group_state_rows[group]) {
        goto fail;
      }
      if (group_state_runs[group] == 0u) {
        if (state_candidate_group_count == UF_long_max) {
          goto fail;
        }
        state_candidate_group_count++;
      }
      group_state_rows[group] += row_count;
      group_state_runs[group]++;
      if (row_count > group_state_max_rows[group]) {
        group_state_max_rows[group] = row_count;
      }
      state_run_count++;
      if (row_count > state_max_rows) {
        state_max_rows = row_count;
      }
    }
    prefix += row_count;
    state_ptr[run + 1u] = focus_largest_groups ? 0u : prefix;
  }

  size_t state_bytes = 0u;
#define KLS_ADD_STATE_BYTES(count_, type_) \
  do { \
    const size_t kls_count__ = (size_t)(count_); \
    if ((UF_long)kls_count__ != (count_) || \
        kls_count__ > SIZE_MAX / sizeof(type_) || \
        state_bytes > SIZE_MAX - kls_count__ * sizeof(type_)) { \
      goto fail; \
    } \
    state_bytes += kls_count__ * sizeof(type_); \
  } while (0)
  KLS_ADD_STATE_BYTES(run_count + 1u, UF_long);
  KLS_ADD_STATE_BYTES(prefix, UF_long);
#undef KLS_ADD_STATE_BYTES

  if (!focus_largest_groups &&
      state_bytes > KLS_REFACTOR_PLAN_GROUP_L_STATE_MAX_BYTES) {
    free(state_ptr);
    free(stamp_workspace);
    free(run_state_rows);
    free(group_state_rows);
    free(group_state_runs);
    free(group_state_max_rows);
    solver->refactor_supernode_consumer_plan_group_l_state_built = 1;
    solver->refactor_supernode_consumer_plan_group_l_state_storage_limited =
      1;
    solver->refactor_supernode_consumer_plan_group_l_state_focus_enabled = 0;
    solver->refactor_supernode_consumer_plan_group_l_state_run_count =
      state_run_count;
    solver->refactor_supernode_consumer_plan_group_l_state_rows = prefix;
    solver->refactor_supernode_consumer_plan_group_l_state_max_rows =
      state_max_rows;
    solver->refactor_supernode_consumer_plan_group_l_state_bytes =
      state_bytes;
    solver
      ->refactor_supernode_consumer_plan_group_l_state_candidate_group_count =
        state_candidate_group_count;
    solver
      ->refactor_supernode_consumer_plan_group_l_state_candidate_run_count =
        state_run_count;
    solver->refactor_supernode_consumer_plan_group_l_state_candidate_rows =
      prefix;
    return 1;
  }

  UF_long selected_rows = prefix;
  UF_long selected_run_count = state_run_count;
  UF_long selected_group_count = state_candidate_group_count;
  UF_long selected_max_rows = state_max_rows;
  size_t selected_bytes = state_bytes;
  if (focus_largest_groups) {
    size_t ptr_bytes = 0u;
#define KLS_ADD_FOCUS_PTR_BYTES(count_, type_) \
    do { \
      const size_t kls_count__ = (size_t)(count_); \
      if ((UF_long)kls_count__ != (count_) || \
          kls_count__ > SIZE_MAX / sizeof(type_) || \
          ptr_bytes > SIZE_MAX - kls_count__ * sizeof(type_)) { \
        goto fail; \
      } \
      ptr_bytes += kls_count__ * sizeof(type_); \
    } while (0)
    KLS_ADD_FOCUS_PTR_BYTES(run_count + 1u, UF_long);
#undef KLS_ADD_FOCUS_PTR_BYTES
    selected_bytes = ptr_bytes;
    selected_rows = 0u;
    selected_run_count = 0u;
    selected_group_count = 0u;
    selected_max_rows = 0u;
    for (UF_long pick = 0u;
         pick < KLS_REFACTOR_PLAN_GROUP_L_STATE_FOCUS_MAX_GROUPS; ++pick) {
      UF_long best_group = KLS_KLU_EMPTY;
      UF_long best_rows = 0u;
      for (UF_long group = 0u; group < group_count; ++group) {
        if (group_selected[group] || group_state_rows[group] == 0u ||
            group_state_rows[group] > (UF_long)(SIZE_MAX / sizeof(UF_long))) {
          continue;
        }
        const size_t group_bytes =
          (size_t)group_state_rows[group] * sizeof(UF_long);
        if (selected_bytes > KLS_REFACTOR_PLAN_GROUP_L_STATE_FOCUS_MAX_BYTES ||
            group_bytes >
              KLS_REFACTOR_PLAN_GROUP_L_STATE_FOCUS_MAX_BYTES -
                selected_bytes) {
          continue;
        }
        if (best_group == KLS_KLU_EMPTY || group_state_rows[group] > best_rows) {
          best_group = group;
          best_rows = group_state_rows[group];
        }
      }
      if (best_group == KLS_KLU_EMPTY) {
        break;
      }
      group_selected[best_group] = 1u;
      solver->refactor_supernode_consumer_plan_group_l_valid[best_group] |=
        KLS_REFACTOR_PLAN_GROUP_L_STATE_FOCUS;
      selected_group_count++;
      selected_run_count += group_state_runs[best_group];
      selected_rows += group_state_rows[best_group];
      if (group_state_max_rows[best_group] > selected_max_rows) {
        selected_max_rows = group_state_max_rows[best_group];
      }
      selected_bytes += (size_t)group_state_rows[best_group] * sizeof(UF_long);
    }

    UF_long selected_prefix = 0u;
    state_ptr[0] = 0u;
    for (UF_long run = 0u; run < run_count; ++run) {
      const UF_long row_count = run_state_rows[run];
      if (row_count > 0u) {
        const UF_long group =
          solver->refactor_supernode_consumer_plan_group_l_run_group[run];
        if (group == KLS_KLU_EMPTY || group >= group_count) {
          goto fail;
        }
        if (group_selected[group]) {
          if (row_count > UF_long_max - selected_prefix) {
            goto fail;
          }
          selected_prefix += row_count;
        }
      }
      state_ptr[run + 1u] = selected_prefix;
    }
    if (selected_prefix != selected_rows) {
      goto fail;
    }
    state_bytes = selected_bytes;
  }

  state_row_idx =
    selected_rows > 0u
      ? (UF_long *)malloc((size_t)selected_rows * sizeof(*state_row_idx))
      : NULL;
  if (selected_rows > 0u && state_row_idx == NULL) {
    goto fail;
  }

  stamp = 0u;
  memset(stamp_workspace, 0, (size_t)stamp_size * sizeof(*stamp_workspace));
  for (UF_long run = 0u; run < run_count; ++run) {
    UF_long run_stamp = 0u;
    UF_long row_count = 0u;
    const UF_long begin = state_ptr[run];
    const UF_long end = state_ptr[run + 1u];
    if (begin > end || end > selected_rows) {
      goto fail;
    }
    if (begin == end) {
      continue;
    }
    if (
        !kls_refactor_supernode_algorithm5_next_target_stamp(
          stamp_workspace, stamp_size, &stamp, &run_stamp) ||
        !kls_collect_refactor_supernode_consumer_plan_group_l_state_rows(
          solver, run, stamp_workspace, stamp_size, run_stamp,
          state_row_idx != NULL ? state_row_idx + begin : NULL,
          end - begin, &row_count) ||
        row_count != end - begin) {
      goto fail;
    }
  }

  free(stamp_workspace);
  free(run_state_rows);
  free(group_state_rows);
  free(group_state_runs);
  free(group_state_max_rows);
  free(group_selected);
  solver->refactor_supernode_consumer_plan_group_l_state_ptr = state_ptr;
  solver->refactor_supernode_consumer_plan_group_l_state_row_idx =
    state_row_idx;
  solver->refactor_supernode_consumer_plan_group_l_state_built = 1;
  solver->refactor_supernode_consumer_plan_group_l_state_storage_limited =
    focus_largest_groups && selected_rows < prefix ? 1 : 0;
  solver->refactor_supernode_consumer_plan_group_l_state_focus_enabled =
    focus_largest_groups ? 1 : 0;
  solver->refactor_supernode_consumer_plan_group_l_state_run_count =
    selected_run_count;
  solver->refactor_supernode_consumer_plan_group_l_state_rows = selected_rows;
  solver->refactor_supernode_consumer_plan_group_l_state_max_rows =
    selected_max_rows;
  solver->refactor_supernode_consumer_plan_group_l_state_bytes = state_bytes;
  solver->refactor_supernode_consumer_plan_group_l_state_candidate_group_count =
    state_candidate_group_count;
  solver->refactor_supernode_consumer_plan_group_l_state_candidate_run_count =
    state_run_count;
  solver->refactor_supernode_consumer_plan_group_l_state_candidate_rows = prefix;
  solver->refactor_supernode_consumer_plan_group_l_state_selected_group_count =
    selected_group_count;
  solver->refactor_supernode_consumer_plan_group_l_state_selected_run_count =
    selected_run_count;
  solver->refactor_supernode_consumer_plan_group_l_state_selected_rows =
    selected_rows;
  solver->refactor_supernode_consumer_plan_group_l_state_selected_bytes =
    state_bytes;
  return 1;

fail:
  free(state_ptr);
  free(state_row_idx);
  free(stamp_workspace);
  free(run_state_rows);
  free(group_state_rows);
  free(group_state_runs);
  free(group_state_max_rows);
  free(group_selected);
  return 0;
}

static int kls_egraph_cached_supernode_update_allows(
  UF_long run_rows,
  UF_long dense_cols,
  UF_long trailing_len);

static int kls_mark_refactor_supernode_used_candidates_from_plan(
  const kls_solver *solver,
  const UF_long *candidate_start_id,
  const UF_long *candidate_start,
  const UF_long *candidate_width,
  const UF_long *candidate_trailing_len,
  UF_long candidate_count,
  unsigned char *used) {
  if (solver == NULL || candidate_start_id == NULL ||
      candidate_start == NULL || candidate_width == NULL ||
      candidate_trailing_len == NULL || used == NULL ||
      solver->refactor_supernode_consumer_plan_run_count == 0u ||
      solver->refactor_supernode_consumer_plan_current == NULL ||
      solver->refactor_supernode_consumer_plan_dep == NULL ||
      solver->refactor_supernode_consumer_plan_rows == NULL ||
      solver->refactor_supernode_consumer_plan_panel_start == NULL ||
      solver->refactor_supernode_consumer_plan_panel_offset == NULL) {
    return 0;
  }

  for (UF_long run = 0u;
       run < solver->refactor_supernode_consumer_plan_run_count; ++run) {
    const UF_long dep_global =
      solver->refactor_supernode_consumer_plan_dep[run];
    const UF_long current =
      solver->refactor_supernode_consumer_plan_current[run];
    const UF_long run_rows =
      solver->refactor_supernode_consumer_plan_rows[run];
    const UF_long plan_start =
      solver->refactor_supernode_consumer_plan_panel_start[run];
    const UF_long panel_offset =
      solver->refactor_supernode_consumer_plan_panel_offset[run];
    if (dep_global >= solver->n || current >= solver->n ||
        dep_global >= current || run_rows <= 1u ||
        plan_start >= solver->n || dep_global < plan_start ||
        panel_offset != dep_global - plan_start) {
      continue;
    }
    const UF_long candidate = candidate_start_id[plan_start];
    if (candidate == KLS_KLU_EMPTY || candidate >= candidate_count) {
      continue;
    }
    const UF_long start = candidate_start[candidate];
    const UF_long width = candidate_width[candidate];
    const UF_long trailing_len = candidate_trailing_len[candidate];
    if (start != plan_start || width <= 1u || start > dep_global ||
        dep_global >= start + width || start > solver->n ||
        width > solver->n - start || run_rows > width - panel_offset) {
      continue;
    }
    const UF_long available_end = panel_offset + run_rows;
    if (available_end > width || start + available_end > current) {
      continue;
    }
    const UF_long dense_scatter_len = width - available_end;
    if (kls_egraph_cached_supernode_update_allows(
          run_rows, dense_scatter_len, trailing_len)) {
      used[candidate] = 1u;
    }
  }
  return 1;
}

static int kls_mark_refactor_supernode_used_candidates(
  const kls_solver *solver,
  const UF_long *candidate_start_id,
  const UF_long *candidate_col_id,
  const UF_long *candidate_start,
  const UF_long *candidate_local_start,
  const UF_long *candidate_width,
  const UF_long *candidate_trailing_len,
  UF_long candidate_count,
  unsigned char *used) {
  if (solver == NULL || candidate_start_id == NULL ||
      candidate_col_id == NULL || candidate_start == NULL ||
      candidate_local_start == NULL || candidate_width == NULL ||
      candidate_trailing_len == NULL || used == NULL) {
    return 0;
  }
  if (kls_mark_refactor_supernode_used_candidates_from_plan(
        solver, candidate_start_id, candidate_start, candidate_width,
        candidate_trailing_len, candidate_count, used)) {
    return 1;
  }
  if (solver->symbolic == NULL || solver->symbolic->R == NULL ||
      solver->numeric == NULL || solver->numeric->Ulen == NULL ||
      (solver->refactor_u_indices == NULL &&
       solver->refactor_u_indices32 == NULL)) {
    return 0;
  }

  const UF_long nblocks = solver->symbolic->nblocks;
  for (UF_long current = 0; current < solver->n; ++current) {
    UF_long block = 0u;
    if (nblocks > 1u) {
      if (solver->refactor_col_block == NULL) {
        return 0;
      }
      block = solver->refactor_col_block[current];
      if (block >= nblocks) {
        continue;
      }
    } else if (nblocks == 0u) {
      continue;
    }
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    if (current < k1 || current >= k2) {
      continue;
    }
    const UF_long current_local = current - k1;
    const UF_long ulen = solver->numeric->Ulen[current];
    const int32_t *ui32 =
      solver->refactor_u_indices32 != NULL
        ? solver->refactor_u_indices32[current] : NULL;
    const UF_long *ui =
      ui32 == NULL && solver->refactor_u_indices != NULL
        ? solver->refactor_u_indices[current] : NULL;
    if (ulen == 0u || (ui32 == NULL && ui == NULL)) {
      continue;
    }

    for (UF_long up = 0; up < ulen; ++up) {
      const UF_long dep_local = ui32 != NULL ? (UF_long)ui32[up] : ui[up];
      if (dep_local >= current_local) {
        continue;
      }
      const UF_long dep_global = k1 + dep_local;
      if (dep_global >= current || dep_global >= solver->n) {
        continue;
      }
      const UF_long candidate = candidate_col_id[dep_global];
      if (candidate == KLS_KLU_EMPTY || candidate >= candidate_count) {
        continue;
      }
      const UF_long start = candidate_start[candidate];
      const UF_long local_start = candidate_local_start[candidate];
      const UF_long width = candidate_width[candidate];
      const UF_long trailing_len = candidate_trailing_len[candidate];
      if (width <= 1u || start > dep_global ||
          dep_global >= start + width || start > solver->n ||
          width > solver->n - start) {
        continue;
      }
      const UF_long panel_offset = dep_global - start;
      if (local_start > dep_local ||
          dep_local - local_start != panel_offset) {
        continue;
      }
      UF_long available_end = width;
      if (start + width > current) {
        available_end = current - start;
      }
      if (available_end > width) {
        available_end = width;
      }
      if (available_end <= panel_offset + 1u) {
        continue;
      }
      const UF_long run_rows = available_end - panel_offset;
      if (run_rows > ulen - up || run_rows > current_local - dep_local) {
        continue;
      }
      int contiguous = 1;
      for (UF_long local = 0; local < run_rows; ++local) {
        const UF_long run_dep_local =
          ui32 != NULL ? (UF_long)ui32[up + local] : ui[up + local];
        if (run_dep_local != dep_local + local) {
          contiguous = 0;
          break;
        }
      }
      const UF_long dense_scatter_len = width - available_end;
      if (contiguous &&
          kls_egraph_cached_supernode_update_allows(
            run_rows, dense_scatter_len, trailing_len)) {
        used[candidate] = 1u;
      }
    }
  }
  return 1;
}

static int kls_refactor_supernode_panel_cache_ready(kls_solver *solver) {
  if (solver == NULL || solver->n == 0u ||
      solver->refactor_lu_pointer_count != solver->n ||
      solver->refactor_supernode_pipeline_end == NULL ||
      solver->refactor_supernode_panel_count == 0u ||
      solver->refactor_supernode_panel_start_id == NULL ||
      solver->refactor_supernode_panel_col_id == NULL ||
      solver->refactor_supernode_panel_start == NULL ||
      solver->refactor_supernode_panel_local_start == NULL ||
      solver->refactor_supernode_panel_width == NULL ||
      solver->refactor_supernode_panel_dense_begin == NULL ||
      solver->refactor_supernode_panel_trailing_begin == NULL ||
      solver->refactor_supernode_panel_trailing_value_begin == NULL ||
      solver->refactor_supernode_panel_trailing_len == NULL ||
      solver->refactor_supernode_panel_dense_values == NULL) {
    return 0;
  }

  UF_long used_count = 0u;
  const UF_long panel_count = solver->refactor_supernode_panel_count;
  for (UF_long panel = 0; panel < panel_count; ++panel) {
    const UF_long start = solver->refactor_supernode_panel_start[panel];
    const UF_long local_start =
      solver->refactor_supernode_panel_local_start[panel];
    const UF_long width = solver->refactor_supernode_panel_width[panel];
    const UF_long trailing_len =
      solver->refactor_supernode_panel_trailing_len[panel];
    if (width == 0u) {
      continue;
    }
    if (width <= 1u || start > solver->n || width > solver->n - start ||
        local_start > solver->n ||
        (trailing_len > 0u &&
         (solver->refactor_supernode_panel_trailing_rows == NULL ||
          solver->refactor_supernode_panel_trailing_values == NULL))) {
      return 0;
    }
    used_count++;
  }
  solver->refactor_supernode_panel_used_count = used_count;
  return 1;
}

static void kls_count_consumer_plan_strict_cached_panel_candidates(
  kls_solver *solver) {
  if (solver == NULL) {
    return;
  }
  solver->refactor_supernode_consumer_plan_cached_panel_count = 0u;
  solver->refactor_supernode_consumer_plan_cached_panel_rows = 0u;
  solver->refactor_supernode_consumer_plan_strict_cached_panel_count = 0u;
  solver->refactor_supernode_consumer_plan_strict_cached_panel_rows = 0u;
  solver->refactor_supernode_consumer_plan_strict_cached_panel_trailing_entries =
    0u;
  if (!solver->refactor_supernode_consumer_plan_built ||
      solver->refactor_supernode_consumer_plan_ptr == NULL ||
      solver->refactor_supernode_consumer_plan_rows == NULL ||
      solver->refactor_supernode_consumer_plan_panel_start == NULL ||
      solver->refactor_supernode_consumer_plan_panel_offset == NULL ||
      solver->refactor_l_indices == NULL ||
      solver->refactor_l_values == NULL ||
      solver->numeric == NULL || solver->numeric->Llen == NULL) {
    return;
  }

  for (UF_long start = 0u; start < solver->n; ++start) {
    const UF_long begin =
      solver->refactor_supernode_consumer_plan_ptr[start];
    const UF_long end =
      solver->refactor_supernode_consumer_plan_ptr[start + 1u];
    if (begin > end ||
        end > solver->refactor_supernode_consumer_plan_run_count) {
      return;
    }
    UF_long width = 0u;
    for (UF_long run = begin; run < end; ++run) {
      if (solver->refactor_supernode_consumer_plan_panel_start[run] !=
          start) {
        continue;
      }
      const UF_long rows =
        solver->refactor_supernode_consumer_plan_rows[run];
      const UF_long offset =
        solver->refactor_supernode_consumer_plan_panel_offset[run];
      if (rows <= 1u || offset > solver->n - start ||
          rows > solver->n - start - offset) {
        continue;
      }
      const UF_long run_end = offset + rows;
      if (run_end > width) {
        width = run_end;
      }
    }
    if (width <= 1u) {
      continue;
    }
    solver->refactor_supernode_consumer_plan_cached_panel_count++;
    solver->refactor_supernode_consumer_plan_cached_panel_rows += width;

    UF_long block_size = 0u;
    UF_long local_start = 0u;
    if (!kls_refactor_supernode_panel_block_info(
          solver, start, start + width, NULL, &block_size, &local_start)) {
      continue;
    }
    UF_long trailing_len = 0u;
    if (!kls_refactor_supernode_panel_candidate(
          solver, start, start + width, local_start, block_size,
          &trailing_len, NULL)) {
      continue;
    }
    solver->refactor_supernode_consumer_plan_strict_cached_panel_count++;
    solver->refactor_supernode_consumer_plan_strict_cached_panel_rows +=
      width;
    if (trailing_len > 0u &&
        width >
          (UF_long_max -
           solver
             ->refactor_supernode_consumer_plan_strict_cached_panel_trailing_entries) /
            trailing_len) {
      return;
    }
    solver
      ->refactor_supernode_consumer_plan_strict_cached_panel_trailing_entries +=
      width * trailing_len;
  }
}

static int kls_build_refactor_supernode_panel_cache(kls_solver *solver) {
  if (solver == NULL || solver->n == 0u ||
      solver->refactor_supernode_pipeline_end == NULL ||
      solver->numeric == NULL || solver->numeric->Llen == NULL ||
      solver->refactor_lu_pointer_count != solver->n ||
      solver->refactor_l_indices == NULL ||
      solver->refactor_l_values == NULL) {
    return 0;
  }

  if (kls_refactor_supernode_panel_cache_ready(solver)) {
    return 1;
  }

  free_refactor_supernode_panel_cache(solver);
  kls_count_consumer_plan_strict_cached_panel_candidates(solver);

  UF_long candidate_count = 0u;
  for (UF_long start = 0; start < solver->n; ++start) {
    const UF_long end = solver->refactor_supernode_pipeline_end[start];
    if (end <= start + 1u || end > solver->n) {
      continue;
    }
    UF_long block_size = 0u;
    UF_long local_start = 0u;
    if (!kls_refactor_supernode_panel_block_info(
          solver, start, end, NULL, &block_size, &local_start)) {
      continue;
    }
    UF_long trailing_len = 0u;
    if (!kls_refactor_supernode_panel_candidate(
          solver, start, end, local_start, block_size, &trailing_len,
          NULL)) {
      continue;
    }
    candidate_count++;
  }
  if (candidate_count == 0u) {
    return 1;
  }
  if (solver->n > (UF_long)(SIZE_MAX / sizeof(UF_long)) ||
      candidate_count > (UF_long)(SIZE_MAX / sizeof(UF_long)) ||
      candidate_count > (UF_long)(PTRDIFF_MAX / sizeof(unsigned char))) {
    return 0;
  }

  UF_long *candidate_start_id =
    (UF_long *)malloc((size_t)solver->n * sizeof(*candidate_start_id));
  UF_long *candidate_col_id =
    (UF_long *)malloc((size_t)solver->n * sizeof(*candidate_col_id));
  UF_long *candidate_start =
    (UF_long *)malloc((size_t)candidate_count * sizeof(*candidate_start));
  UF_long *candidate_local_start =
    (UF_long *)malloc((size_t)candidate_count *
                      sizeof(*candidate_local_start));
  UF_long *candidate_width =
    (UF_long *)malloc((size_t)candidate_count * sizeof(*candidate_width));
  UF_long *candidate_trailing_len =
    (UF_long *)malloc((size_t)candidate_count *
                      sizeof(*candidate_trailing_len));
  unsigned char *candidate_used =
    (unsigned char *)calloc((size_t)candidate_count,
                            sizeof(*candidate_used));
#define KLS_FREE_PANEL_CANDIDATES() do { \
    free(candidate_start_id);                 \
    free(candidate_col_id);                   \
    free(candidate_start);                    \
    free(candidate_local_start);              \
    free(candidate_width);                    \
    free(candidate_trailing_len);             \
    free(candidate_used);                     \
  } while (0)
  if (candidate_start_id == NULL || candidate_col_id == NULL ||
      candidate_start == NULL || candidate_local_start == NULL ||
      candidate_width == NULL || candidate_trailing_len == NULL ||
      candidate_used == NULL) {
    KLS_FREE_PANEL_CANDIDATES();
    return 0;
  }
  for (UF_long k = 0; k < solver->n; ++k) {
    candidate_start_id[k] = KLS_KLU_EMPTY;
    candidate_col_id[k] = KLS_KLU_EMPTY;
  }

  UF_long candidate = 0u;
  for (UF_long start = 0; start < solver->n; ++start) {
    const UF_long end = solver->refactor_supernode_pipeline_end[start];
    if (end <= start + 1u || end > solver->n) {
      continue;
    }
    UF_long block_size = 0u;
    UF_long local_start = 0u;
    if (!kls_refactor_supernode_panel_block_info(
          solver, start, end, NULL, &block_size, &local_start)) {
      continue;
    }
    UF_long trailing_len = 0u;
    if (!kls_refactor_supernode_panel_candidate(
          solver, start, end, local_start, block_size, &trailing_len,
          NULL)) {
      continue;
    }
    if (candidate >= candidate_count) {
      KLS_FREE_PANEL_CANDIDATES();
      return 0;
    }
    const UF_long width = end - start;
    candidate_start[candidate] = start;
    candidate_local_start[candidate] = local_start;
    candidate_width[candidate] = width;
    candidate_trailing_len[candidate] = trailing_len;
    candidate_start_id[start] = candidate;
    for (UF_long col = start; col < end; ++col) {
      candidate_col_id[col] = candidate;
    }
    candidate++;
  }
  if (candidate != candidate_count) {
    KLS_FREE_PANEL_CANDIDATES();
    return 0;
  }
  if (!kls_mark_refactor_supernode_used_candidates(
        solver, candidate_start_id, candidate_col_id, candidate_start,
        candidate_local_start, candidate_width, candidate_trailing_len,
        candidate_count, candidate_used)) {
    memset(candidate_used, 1, (size_t)candidate_count *
                              sizeof(*candidate_used));
  }

  UF_long panel_count = 0u;
  UF_long dense_total = 0u;
  UF_long trailing_row_total = 0u;
  UF_long trailing_value_total = 0u;
  for (UF_long c = 0u; c < candidate_count; ++c) {
    if (!candidate_used[c]) {
      continue;
    }
    const UF_long width = candidate_width[c];
    const UF_long trailing_len = candidate_trailing_len[c];
    if (width > UF_long_max / width) {
      KLS_FREE_PANEL_CANDIDATES();
      return 0;
    }
    const UF_long dense_entries = width * width;
    if (trailing_len > 0u &&
        width > (UF_long_max - trailing_value_total) / trailing_len) {
      KLS_FREE_PANEL_CANDIDATES();
      return 0;
    }
    const UF_long trailing_entries = width * trailing_len;
    if (dense_entries > UF_long_max - dense_total ||
        trailing_len > UF_long_max - trailing_row_total ||
        trailing_entries > UF_long_max - trailing_value_total) {
      KLS_FREE_PANEL_CANDIDATES();
      return 0;
    }
    dense_total += dense_entries;
    trailing_row_total += trailing_len;
    trailing_value_total += trailing_entries;
    panel_count++;
  }
  if (panel_count == 0u) {
    KLS_FREE_PANEL_CANDIDATES();
    return 1;
  }
  if (panel_count > (UF_long)(SIZE_MAX / sizeof(UF_long)) ||
      dense_total > (UF_long)(SIZE_MAX / sizeof(double)) ||
      trailing_row_total > (UF_long)(SIZE_MAX / sizeof(UF_long)) ||
      trailing_value_total > (UF_long)(SIZE_MAX / sizeof(double))) {
    KLS_FREE_PANEL_CANDIDATES();
    return 0;
  }

  UF_long *start_id =
    (UF_long *)malloc((size_t)solver->n * sizeof(*start_id));
  UF_long *col_id =
    (UF_long *)malloc((size_t)solver->n * sizeof(*col_id));
  UF_long *panel_start =
    (UF_long *)malloc((size_t)panel_count * sizeof(*panel_start));
  UF_long *panel_local_start =
    (UF_long *)malloc((size_t)panel_count * sizeof(*panel_local_start));
  UF_long *panel_width =
    (UF_long *)malloc((size_t)panel_count * sizeof(*panel_width));
  UF_long *panel_dense_begin =
    (UF_long *)malloc((size_t)panel_count * sizeof(*panel_dense_begin));
  UF_long *panel_trailing_begin =
    (UF_long *)malloc((size_t)panel_count * sizeof(*panel_trailing_begin));
  UF_long *panel_trailing_value_begin =
    (UF_long *)malloc((size_t)panel_count *
                      sizeof(*panel_trailing_value_begin));
  UF_long *panel_trailing_len =
    (UF_long *)malloc((size_t)panel_count * sizeof(*panel_trailing_len));
  UF_long *trailing_rows = trailing_row_total > 0u
    ? (UF_long *)malloc((size_t)trailing_row_total *
                        sizeof(*trailing_rows)) : NULL;
  double *dense_values =
    (double *)calloc((size_t)dense_total, sizeof(*dense_values));
  double *trailing_values = trailing_value_total > 0u
    ? (double *)calloc((size_t)trailing_value_total,
                       sizeof(*trailing_values)) : NULL;
#define KLS_FREE_PANEL_STORAGE() do { \
    free(start_id);                       \
    free(col_id);                         \
    free(panel_start);                    \
    free(panel_local_start);              \
    free(panel_width);                    \
    free(panel_dense_begin);              \
    free(panel_trailing_begin);           \
    free(panel_trailing_value_begin);     \
    free(panel_trailing_len);             \
    free(trailing_rows);                  \
    free(dense_values);                   \
    free(trailing_values);                \
    KLS_FREE_PANEL_CANDIDATES();          \
  } while (0)
  if (start_id == NULL || col_id == NULL || panel_start == NULL ||
      panel_local_start == NULL || panel_width == NULL ||
      panel_dense_begin == NULL || panel_trailing_begin == NULL ||
      panel_trailing_value_begin == NULL || panel_trailing_len == NULL ||
      (trailing_row_total > 0u && trailing_rows == NULL) ||
      dense_values == NULL ||
      (trailing_value_total > 0u && trailing_values == NULL)) {
    KLS_FREE_PANEL_STORAGE();
    return 0;
  }
  for (UF_long k = 0; k < solver->n; ++k) {
    start_id[k] = KLS_KLU_EMPTY;
    col_id[k] = KLS_KLU_EMPTY;
  }

  UF_long panel = 0u;
  UF_long dense_pos = 0u;
  UF_long trailing_row_pos = 0u;
  UF_long trailing_value_pos = 0u;
  for (UF_long c = 0u; c < candidate_count; ++c) {
    if (!candidate_used[c]) {
      continue;
    }
    const UF_long start = candidate_start[c];
    const UF_long local_start = candidate_local_start[c];
    const UF_long width = candidate_width[c];
    UF_long block_size = 0u;
    UF_long checked_local_start = 0u;
    UF_long trailing_len = candidate_trailing_len[c];
    UF_long *rows_out = trailing_len > 0u ? trailing_rows + trailing_row_pos
                                          : NULL;
    if (!kls_refactor_supernode_panel_block_info(
          solver, start, start + width, NULL, &block_size,
          &checked_local_start) ||
        checked_local_start != local_start ||
        !kls_refactor_supernode_panel_candidate(
          solver, start, start + width, local_start,
          block_size, &trailing_len, rows_out) ||
        trailing_len != candidate_trailing_len[c]) {
      KLS_FREE_PANEL_STORAGE();
      return 0;
    }
    if (panel >= panel_count) {
      KLS_FREE_PANEL_STORAGE();
      return 0;
    }
    panel_start[panel] = start;
    panel_local_start[panel] = local_start;
    panel_width[panel] = width;
    panel_dense_begin[panel] = dense_pos;
    panel_trailing_begin[panel] = trailing_row_pos;
    panel_trailing_value_begin[panel] = trailing_value_pos;
    panel_trailing_len[panel] = trailing_len;
    start_id[start] = panel;
    for (UF_long col = start; col < start + width; ++col) {
      col_id[col] = panel;
    }
    dense_pos += width * width;
    trailing_row_pos += trailing_len;
    trailing_value_pos += width * trailing_len;
    panel++;
  }
  if (panel != panel_count || dense_pos != dense_total ||
      trailing_row_pos != trailing_row_total ||
      trailing_value_pos != trailing_value_total) {
    KLS_FREE_PANEL_STORAGE();
    return 0;
  }
  KLS_FREE_PANEL_CANDIDATES();

  solver->refactor_supernode_panel_start_id = start_id;
  solver->refactor_supernode_panel_col_id = col_id;
  solver->refactor_supernode_panel_start = panel_start;
  solver->refactor_supernode_panel_local_start = panel_local_start;
  solver->refactor_supernode_panel_width = panel_width;
  solver->refactor_supernode_panel_dense_begin = panel_dense_begin;
  solver->refactor_supernode_panel_trailing_begin = panel_trailing_begin;
  solver->refactor_supernode_panel_trailing_value_begin =
    panel_trailing_value_begin;
  solver->refactor_supernode_panel_trailing_len = panel_trailing_len;
  solver->refactor_supernode_panel_trailing_rows = trailing_rows;
  solver->refactor_supernode_panel_dense_values = dense_values;
  solver->refactor_supernode_panel_trailing_values = trailing_values;
  solver->refactor_supernode_panel_count = panel_count;
  solver->refactor_supernode_panel_used_count = panel_count;
#undef KLS_FREE_PANEL_STORAGE
#undef KLS_FREE_PANEL_CANDIDATES
  return 1;
}

static UF_long kls_numeric_l_len_in_block(const kls_solver *solver,
                                          UF_long k1,
                                          UF_long k2,
                                          UF_long k) {
  if (solver == NULL || solver->numeric == NULL ||
      solver->numeric->Llen == NULL || k < k1 || k >= k2 ||
      k2 - k1 <= 1u) {
    return 0u;
  }
  return solver->numeric->Llen[k];
}

static UF_long kls_numeric_u_len_in_block(const kls_solver *solver,
                                          UF_long k1,
                                          UF_long k2,
                                          UF_long k) {
  if (solver == NULL || solver->numeric == NULL ||
      solver->numeric->Ulen == NULL || k < k1 || k >= k2 ||
      k2 - k1 <= 1u) {
    return 0u;
  }
  return solver->numeric->Ulen[k];
}

static int kls_row_refactor_pattern_is_eligible(const kls_solver *solver) {
  const int scaled = solver != NULL && solver->common.scale > 0;
  const int btf = solver != NULL && solver->symbolic != NULL &&
                  solver->symbolic->nblocks > 1u;
  return solver != NULL && solver->symbolic != NULL &&
         solver->numeric != NULL && solver->symbolic->nblocks > 0u &&
         (!btf || (solver->symbolic->R != NULL &&
                   solver->refactor_block_start != NULL &&
                   solver->refactor_col_block != NULL &&
                   solver->numeric->Offp != NULL &&
                   solver->numeric->Offx != NULL)) &&
         ((scaled && solver->numeric->Rs != NULL &&
           solver->numeric->Pnum != NULL && solver->numeric->Xwork != NULL) ||
          (!scaled && solver->numeric->Rs == NULL)) &&
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

static int kls_build_row_refactor_input_cleanup_flags(
  UF_long n,
  const UF_long *l_ptr,
  const UF_long *l_cols,
  const UF_long *u_ptr,
  const UF_long *u_cols,
  const UF_long *input_ptr,
  const UF_long *input_cols,
  unsigned char **needs_cleanup_out,
  UF_long *cleanup_rows_out,
  UF_long *cleanup_entries_out) {
  KLS_SET_OPTIONAL_OUTPUT(needs_cleanup_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(cleanup_rows_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(cleanup_entries_out, 0);
  if (n < 0 || (uint64_t)n > (uint64_t)SIZE_MAX ||
      needs_cleanup_out == NULL || cleanup_rows_out == NULL ||
      cleanup_entries_out == NULL ||
      (n > 0u && (l_ptr == NULL || u_ptr == NULL ||
                  input_ptr == NULL))) {
    return 0;
  }
  if (n == 0u) {
    return 1;
  }

  unsigned char *needs_cleanup =
    (unsigned char *)calloc((size_t)n, sizeof(*needs_cleanup));
  if (needs_cleanup == NULL) {
    return 0;
  }

  UF_long cleanup_rows = 0;
  UF_long cleanup_entries = 0;
  for (UF_long row = 0; row < n; ++row) {
    UF_long lp = l_ptr[row];
    const UF_long lend = l_ptr[row + 1u];
    UF_long up = u_ptr[row];
    const UF_long uend = u_ptr[row + 1u];
    const UF_long ibegin = input_ptr[row];
    const UF_long iend = input_ptr[row + 1u];
    if ((lend > lp && l_cols == NULL) ||
        (uend > up && u_cols == NULL) ||
        (iend > ibegin && input_cols == NULL)) {
      free(needs_cleanup);
      return 0;
    }
    for (UF_long pos = ibegin; pos < iend; ++pos) {
      const UF_long col = input_cols[pos];
      int covered = col == row;
      if (col < row) {
        while (lp < lend && l_cols[lp] < col) {
          lp++;
        }
        covered = lp < lend && l_cols[lp] == col;
      } else if (col > row) {
        while (up < uend && u_cols[up] < col) {
          up++;
        }
        covered = up < uend && u_cols[up] == col;
      }
      if (!covered) {
        needs_cleanup[row] = 1u;
      }
    }
    if (needs_cleanup[row]) {
      cleanup_rows++;
      cleanup_entries += iend - ibegin;
    }
  }

  *needs_cleanup_out = needs_cleanup;
  *cleanup_rows_out = cleanup_rows;
  *cleanup_entries_out = cleanup_entries;
  return 1;
}

typedef struct {
  UF_long n;
  const UF_long *l_ptr;
  const UF_long *l_cols;
  const UF_long *u_ptr;
  const UF_long *u_cols;
  const UF_long *input_ptr;
  const UF_long *input_cols;
  unsigned char *needs_cleanup;
  UF_long cleanup_rows;
  UF_long cleanup_entries;
  double seconds;
  int ok;
} kls_row_input_cleanup_prep_job;

static void *kls_row_input_cleanup_prep_main(void *arg) {
  kls_row_input_cleanup_prep_job *job =
    (kls_row_input_cleanup_prep_job *)arg;
  if (job == NULL) {
    return NULL;
  }
  const double begin = kls_now_seconds();
  job->ok = kls_build_row_refactor_input_cleanup_flags(
    job->n, job->l_ptr, job->l_cols, job->u_ptr, job->u_cols,
    job->input_ptr, job->input_cols, &job->needs_cleanup,
    &job->cleanup_rows, &job->cleanup_entries);
  job->seconds = kls_now_seconds() - begin;
  return NULL;
}

typedef struct {
  UF_long n;
  const UF_long *l_ptr;
  const UF_long *l_cols;
  UF_long *levels;
  UF_long *level_ptr;
  UF_long *level_rows;
  UF_long level_count;
  UF_long max_width;
  double seconds;
  int ok;
} kls_row_level_prep_job;

static void *kls_row_level_prep_main(void *arg) {
  kls_row_level_prep_job *job = (kls_row_level_prep_job *)arg;
  if (job == NULL || job->n < 0 ||
      (uint64_t)job->n > (uint64_t)(SIZE_MAX / sizeof(UF_long))) {
    return NULL;
  }
  const double begin = kls_now_seconds();
  const UF_long n = job->n;
  UF_long *levels = n > 0u
    ? (UF_long *)calloc((size_t)n, sizeof(*levels)) : NULL;
  if (n > 0u && levels == NULL) {
    return NULL;
  }
  UF_long max_level = 0;
  for (UF_long row = 0; row < n; ++row) {
    UF_long level = 0;
    for (UF_long p = job->l_ptr[row]; p < job->l_ptr[row + 1u]; ++p) {
      const UF_long dep = job->l_cols[p];
      if (dep >= row) {
        free(levels);
        return NULL;
      }
      const UF_long dep_level = levels[dep] + 1u;
      if (dep_level > level) {
        level = dep_level;
      }
    }
    levels[row] = level;
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
    return NULL;
  }
  for (UF_long row = 0; row < n; ++row) {
    level_ptr[levels[row] + 1u]++;
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
    return NULL;
  }
  if (level_count > 0u) {
    memcpy(write_pos, level_ptr, (size_t)level_count * sizeof(*write_pos));
  }
  for (UF_long row = 0; row < n; ++row) {
    level_rows[write_pos[levels[row]]++] = row;
  }
  free(write_pos);
  UF_long max_width = 0;
  for (UF_long level = 0; level < level_count; ++level) {
    const UF_long width = level_ptr[level + 1u] - level_ptr[level];
    if (width > max_width) {
      max_width = width;
    }
  }
  job->levels = levels;
  job->level_ptr = level_ptr;
  job->level_rows = level_rows;
  job->level_count = level_count;
  job->max_width = max_width;
  job->seconds = kls_now_seconds() - begin;
  job->ok = 1;
  return NULL;
}

static int kls_row_solve_pattern_is_ready(const kls_solver *solver) {
  return solver != NULL &&
         solver->row_refactor_pattern_n == solver->n &&
         solver->row_refactor_l_ptr != NULL &&
         (solver->row_refactor_l_ptr[solver->n] == 0u ||
          (solver->row_refactor_l_cols != NULL &&
           solver->row_refactor_l_values != NULL &&
           solver->row_refactor_l_row_values != NULL)) &&
         solver->row_refactor_u_ptr != NULL &&
         (solver->row_refactor_u_ptr[solver->n] == 0u ||
          (solver->row_refactor_u_cols != NULL &&
           solver->row_refactor_u_values != NULL &&
           solver->row_refactor_u_row_values != NULL));
}

static UF_long kls_row_solve_dense_tail_start(UF_long n,
                                              const UF_long *ptr,
                                              UF_long *entries_out) {
  KLS_SET_OPTIONAL_OUTPUT(entries_out, 0);
  if (n == 0u || ptr == NULL || ptr[n] == 0u ||
      ptr[n] < KLS_ROW_SOLVE_DENSE_TAIL_MIN_NNZ) {
    return n;
  }

  const double total = (double)ptr[n];
  UF_long suffix = 0;
  for (UF_long remaining = n; remaining > 0u; --remaining) {
    const UF_long row = remaining - 1u;
    if (ptr[row] > ptr[row + 1u]) {
      return n;
    }
    suffix += ptr[row + 1u] - ptr[row];
    if (suffix >= KLS_ROW_SOLVE_DENSE_TAIL_MIN_NNZ &&
        (double)suffix >= KLS_ROW_SOLVE_DENSE_TAIL_MIN_FRACTION * total) {
      if (entries_out != NULL) {
        *entries_out = suffix;
      }
      return row;
    }
  }
  return n;
}

static UF_long *kls_build_row_solve_slice_bounds(UF_long n,
                                                 const UF_long *ptr,
                                                 UF_long start,
                                                 UF_long entries,
                                                 UF_long *max_entries_out) {
  KLS_SET_OPTIONAL_OUTPUT(max_entries_out, 0);
  if (ptr == NULL || start >= n || entries == 0u) {
    return NULL;
  }
  const UF_long slices = KLS_ROW_SOLVE_TRAPEZOID_SLICES;
  UF_long *bounds =
    (UF_long *)malloc(((size_t)slices + 1u) * sizeof(*bounds));
  if (bounds == NULL) {
    return NULL;
  }

  bounds[0] = start;
  bounds[slices] = n;
  UF_long row = start;
  const UF_long base = ptr[start];
  for (UF_long slice = 1u; slice < slices; ++slice) {
    const UF_long target =
      (entries / slices) * slice +
      ((entries % slices) * slice + slices - 1u) / slices;
    while (row < n && ptr[row] - base < target) {
      row++;
    }
    bounds[slice] = row;
  }

  UF_long max_entries = 0;
  for (UF_long slice = 0u; slice < slices; ++slice) {
    const UF_long begin = bounds[slice];
    const UF_long end = bounds[slice + 1u];
    if (begin > end || end > n || ptr[begin] > ptr[end]) {
      free(bounds);
      return NULL;
    }
    const UF_long slice_entries = ptr[end] - ptr[begin];
    if (slice_entries > max_entries) {
      max_entries = slice_entries;
    }
  }
  if (max_entries_out != NULL) {
    *max_entries_out = max_entries;
  }
  return bounds;
}

static UF_long kls_row_solve_col_lower_bound(const UF_long *cols,
                                             UF_long begin,
                                             UF_long end,
                                             UF_long threshold) {
  UF_long lo = begin;
  UF_long hi = end;
  while (lo < hi) {
    const UF_long mid = lo + (hi - lo) / 2u;
    if (cols[mid] < threshold) {
      lo = mid + 1u;
    } else {
      hi = mid;
    }
  }
  return lo;
}

static int kls_build_row_solve_segment_splits(UF_long n,
                                              const UF_long *ptr,
                                              const UF_long *cols,
                                              const UF_long *bounds,
                                              int upper,
                                              UF_long **split_out,
                                              UF_long *rows_out,
                                              UF_long *rect_entries_out,
                                              UF_long *tri_entries_out) {
  KLS_SET_OPTIONAL_OUTPUT(split_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(rows_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(rect_entries_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(tri_entries_out, 0);
  if (split_out == NULL || ptr == NULL || cols == NULL || bounds == NULL ||
      n == 0u) {
    return 1;
  }

  UF_long *split = (UF_long *)malloc((size_t)n * sizeof(*split));
  if (split == NULL) {
    return 0;
  }
  for (UF_long row = 0; row < n; ++row) {
    split[row] = ptr[row];
  }

  UF_long rows = 0;
  UF_long rect_entries = 0;
  UF_long tri_entries = 0;
  for (UF_long slice = 0u; slice < KLS_ROW_SOLVE_TRAPEZOID_SLICES; ++slice) {
    const UF_long slice_begin = bounds[slice];
    const UF_long slice_end = bounds[slice + 1u];
    if (slice_begin > slice_end || slice_end > n) {
      free(split);
      return 0;
    }
    const UF_long threshold = upper ? slice_end : slice_begin;
    for (UF_long row = slice_begin; row < slice_end; ++row) {
      const UF_long row_begin = ptr[row];
      const UF_long row_end = ptr[row + 1u];
      if (row_begin > row_end) {
        free(split);
        return 0;
      }
      for (UF_long p = row_begin; p < row_end; ++p) {
        const UF_long col = cols[p];
        if ((p > row_begin && cols[p - 1u] > col) ||
            (upper ? col <= row : col >= row)) {
          free(split);
          return 0;
        }
      }
      const UF_long cut =
        kls_row_solve_col_lower_bound(cols, row_begin, row_end, threshold);
      split[row] = cut;
      rows++;
      if (upper) {
        tri_entries += cut - row_begin;
        rect_entries += row_end - cut;
      } else {
        rect_entries += cut - row_begin;
        tri_entries += row_end - cut;
      }
    }
  }

  *split_out = split;
  if (rows_out != NULL) {
    *rows_out = rows;
  }
  if (rect_entries_out != NULL) {
    *rect_entries_out = rect_entries;
  }
  if (tri_entries_out != NULL) {
    *tri_entries_out = tri_entries;
  }
  return 1;
}

static UF_long kls_row_solve_rect_len(const UF_long *ptr,
                                      const UF_long *split,
                                      UF_long row,
                                      int upper) {
  return upper ? ptr[row + 1u] - split[row] : split[row] - ptr[row];
}

static int kls_build_one_row_solve_thread_bounds(
  UF_long n,
  const UF_long *ptr,
  const UF_long *split,
  const UF_long *slice_bounds,
  int upper,
  int thread_count,
  UF_long **bounds_out,
  UF_long *max_rect_entries_out) {
  KLS_SET_OPTIONAL_OUTPUT(bounds_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(max_rect_entries_out, 0);
  if (bounds_out == NULL || ptr == NULL || split == NULL ||
      slice_bounds == NULL || thread_count < 2) {
    return 0;
  }

  const UF_long slices = KLS_ROW_SOLVE_TRAPEZOID_SLICES;
  const UF_long workers = (UF_long)thread_count;
  UF_long *bounds =
    (UF_long *)malloc((size_t)slices * ((size_t)workers + 1u) *
                      sizeof(*bounds));
  if (bounds == NULL) {
    return 0;
  }

  UF_long max_rect_entries = 0;
  for (UF_long slice = 0u; slice < slices; ++slice) {
    const UF_long begin = slice_bounds[slice];
    const UF_long end = slice_bounds[slice + 1u];
    if (begin > end || end > n) {
      free(bounds);
      return 0;
    }
    UF_long total = 0;
    for (UF_long row = begin; row < end; ++row) {
      if (ptr[row] > split[row] || split[row] > ptr[row + 1u]) {
        free(bounds);
        return 0;
      }
      total += kls_row_solve_rect_len(ptr, split, row, upper);
    }

    UF_long *parts = bounds + slice * (workers + 1u);
    parts[0] = begin;
    parts[workers] = end;
    UF_long row = begin;
    UF_long prefix = 0;
    for (UF_long worker = 1u; worker < workers; ++worker) {
      const UF_long target =
        (total / workers) * worker +
        ((total % workers) * worker + workers - 1u) / workers;
      while (row < end && prefix < target) {
        prefix += kls_row_solve_rect_len(ptr, split, row, upper);
        row++;
      }
      parts[worker] = row;
    }

    for (UF_long worker = 0u; worker < workers; ++worker) {
      UF_long worker_entries = 0;
      for (UF_long r = parts[worker]; r < parts[worker + 1u]; ++r) {
        worker_entries += kls_row_solve_rect_len(ptr, split, r, upper);
      }
      if (worker_entries > max_rect_entries) {
        max_rect_entries = worker_entries;
      }
    }
  }

  *bounds_out = bounds;
  if (max_rect_entries_out != NULL) {
    *max_rect_entries_out = max_rect_entries;
  }
  return 1;
}

static void kls_clear_row_solve_thread_bounds(kls_solver *solver) {
  if (solver == NULL) {
    return;
  }
  free(solver->row_solve_l_thread_bounds);
  free(solver->row_solve_u_thread_bounds);
  solver->row_solve_l_thread_bounds = NULL;
  solver->row_solve_u_thread_bounds = NULL;
  solver->row_solve_thread_count = 0;
  solver->row_solve_l_thread_max_rect_entries = 0;
  solver->row_solve_u_thread_max_rect_entries = 0;
  solver->row_solve_l_sparse_cluster_levels = 0;
  solver->row_solve_u_sparse_cluster_levels = 0;
  solver->stats.row_solve_thread_count = 0;
  solver->stats.row_solve_l_thread_max_rect_entries = 0;
  solver->stats.row_solve_u_thread_max_rect_entries = 0;
  solver->stats.row_solve_l_sparse_cluster_levels = 0;
  solver->stats.row_solve_u_sparse_cluster_levels = 0;
}

static UF_long kls_row_solve_sparse_cluster_levels(
  const UF_long *level_ptr,
  UF_long level_count,
  int thread_count);

static int kls_build_row_solve_thread_bounds(kls_solver *solver,
                                             int thread_count) {
  if (solver == NULL || thread_count < 2 ||
      solver->row_solve_partition_ready == 0) {
    return 0;
  }
  const int have_l = solver->row_solve_l_slice_bounds != NULL &&
                     solver->row_solve_l_segment_split != NULL;
  const int have_u = solver->row_solve_u_slice_bounds != NULL &&
                     solver->row_solve_u_segment_split != NULL;
  const int have_l_sparse = solver->row_solve_l_sparse_level_ptr != NULL &&
                            solver->row_solve_l_sparse_level_rows != NULL;
  const int have_u_sparse = solver->row_solve_u_sparse_level_ptr != NULL &&
                            solver->row_solve_u_sparse_level_rows != NULL;
  if (!have_l && !have_u && !have_l_sparse && !have_u_sparse) {
    return 0;
  }
  if (solver->row_solve_thread_count == thread_count &&
      (!have_l || solver->row_solve_l_thread_bounds != NULL) &&
      (!have_u || solver->row_solve_u_thread_bounds != NULL)) {
    return 1;
  }

  kls_clear_row_solve_thread_bounds(solver);
  UF_long *l_bounds = NULL;
  UF_long *u_bounds = NULL;
  UF_long l_max_rect = 0;
  UF_long u_max_rect = 0;
  if (have_l) {
    if (!kls_build_one_row_solve_thread_bounds(
          solver->n, solver->row_refactor_l_ptr,
          solver->row_solve_l_segment_split, solver->row_solve_l_slice_bounds,
          0, thread_count, &l_bounds, &l_max_rect)) {
      free(l_bounds);
      return 0;
    }
  }
  if (have_u) {
    if (!kls_build_one_row_solve_thread_bounds(
          solver->n, solver->row_refactor_u_ptr,
          solver->row_solve_u_segment_split, solver->row_solve_u_slice_bounds,
          1, thread_count, &u_bounds, &u_max_rect)) {
      free(l_bounds);
      free(u_bounds);
      return 0;
    }
  }

  solver->row_solve_l_thread_bounds = l_bounds;
  solver->row_solve_u_thread_bounds = u_bounds;
  solver->row_solve_thread_count = thread_count;
  solver->row_solve_l_thread_max_rect_entries = l_max_rect;
  solver->row_solve_u_thread_max_rect_entries = u_max_rect;
  solver->row_solve_l_sparse_cluster_levels =
    kls_row_solve_sparse_cluster_levels(
      solver->row_solve_l_sparse_level_ptr,
      solver->row_solve_l_sparse_level_count, thread_count);
  solver->row_solve_u_sparse_cluster_levels =
    kls_row_solve_sparse_cluster_levels(
      solver->row_solve_u_sparse_level_ptr,
      solver->row_solve_u_sparse_level_count, thread_count);
  solver->stats.row_solve_thread_count = (int64_t)thread_count;
  solver->stats.row_solve_l_thread_max_rect_entries = (int64_t)l_max_rect;
  solver->stats.row_solve_u_thread_max_rect_entries = (int64_t)u_max_rect;
  solver->stats.row_solve_l_sparse_cluster_levels =
    (int64_t)solver->row_solve_l_sparse_cluster_levels;
  solver->stats.row_solve_u_sparse_cluster_levels =
    (int64_t)solver->row_solve_u_sparse_cluster_levels;
  return 1;
}

static int kls_build_row_solve_sparse_levels(
  UF_long n,
  const UF_long *ptr,
  const UF_long *cols,
  UF_long end,
  int upper,
  UF_long **level_ptr_out,
  UF_long **level_rows_out,
  UF_long *level_count_out,
  UF_long *level_max_width_out) {
  KLS_SET_OPTIONAL_OUTPUT(level_ptr_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(level_rows_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(level_count_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(level_max_width_out, 0);
  if (level_ptr_out == NULL || level_rows_out == NULL ||
      ptr == NULL || cols == NULL || end > n) {
    return 0;
  }
  if (end == 0u) {
    return 1;
  }

  UF_long *row_level = (UF_long *)malloc((size_t)end * sizeof(*row_level));
  UF_long *level_counts =
    (UF_long *)calloc((size_t)end, sizeof(*level_counts));
  if (row_level == NULL || level_counts == NULL) {
    free(row_level);
    free(level_counts);
    return 0;
  }

  UF_long max_level = 0;
  if (upper) {
    for (UF_long remaining = end; remaining > 0u; --remaining) {
      const UF_long row = remaining - 1u;
      UF_long level = 0;
      for (UF_long p = ptr[row]; p < ptr[row + 1u]; ++p) {
        const UF_long col = cols[p];
        if (col > row && col < end) {
          const UF_long dep_level = row_level[col] + 1u;
          if (dep_level > level) {
            level = dep_level;
          }
        }
      }
      row_level[row] = level;
      level_counts[level]++;
      if (level > max_level) {
        max_level = level;
      }
    }
  } else {
    for (UF_long row = 0u; row < end; ++row) {
      UF_long level = 0;
      for (UF_long p = ptr[row]; p < ptr[row + 1u]; ++p) {
        const UF_long col = cols[p];
        if (col < row) {
          const UF_long dep_level = row_level[col] + 1u;
          if (dep_level > level) {
            level = dep_level;
          }
        }
      }
      row_level[row] = level;
      level_counts[level]++;
      if (level > max_level) {
        max_level = level;
      }
    }
  }

  const UF_long level_count = max_level + 1u;
  UF_long *level_ptr =
    (UF_long *)malloc((size_t)(level_count + 1u) * sizeof(*level_ptr));
  UF_long *level_rows =
    (UF_long *)malloc((size_t)end * sizeof(*level_rows));
  if (level_ptr == NULL || level_rows == NULL) {
    free(level_ptr);
    free(level_rows);
    free(row_level);
    free(level_counts);
    return 0;
  }

  UF_long max_width = 0;
  level_ptr[0] = 0u;
  for (UF_long level = 0u; level < level_count; ++level) {
    if (level_counts[level] > max_width) {
      max_width = level_counts[level];
    }
    level_ptr[level + 1u] = level_ptr[level] + level_counts[level];
    level_counts[level] = level_ptr[level];
  }
  for (UF_long row = 0u; row < end; ++row) {
    const UF_long level = row_level[row];
    level_rows[level_counts[level]++] = row;
  }

  free(row_level);
  free(level_counts);
  *level_ptr_out = level_ptr;
  *level_rows_out = level_rows;
  if (level_count_out != NULL) {
    *level_count_out = level_count;
  }
  if (level_max_width_out != NULL) {
    *level_max_width_out = max_width;
  }
  return 1;
}

static UF_long kls_row_solve_sparse_cluster_levels(
  const UF_long *level_ptr,
  UF_long level_count,
  int thread_count) {
  if (level_ptr == NULL || thread_count < 2) {
    return 0;
  }
  const UF_long threshold =
    (UF_long)thread_count * KLS_ROW_SOLVE_CLUSTER_ALPHA_NUMERATOR;
  UF_long level = 0;
  for (; level < level_count; ++level) {
    const UF_long width = level_ptr[level + 1u] - level_ptr[level];
    if (width < threshold) {
      break;
    }
  }
  return level;
}

static int kls_build_one_row_solve_transpose_plan_impl(
  kls_solver *solver,
  int source_upper,
  kls_row_solve_transpose_plan *plan,
  int build_parallel_metadata) {
  if (solver == NULL || plan == NULL ||
      solver->row_refactor_pattern_n != solver->n) {
    return 0;
  }
  kls_clear_row_solve_transpose_plan(plan);

  const UF_long n = solver->n;
  const UF_long *src_ptr =
    source_upper ? solver->row_refactor_u_ptr : solver->row_refactor_l_ptr;
  const UF_long *src_cols =
    source_upper ? solver->row_refactor_u_cols : solver->row_refactor_l_cols;
  if (src_ptr == NULL || (src_ptr[n] > 0u && src_cols == NULL)) {
    return 0;
  }
  if (n > (UF_long)(SIZE_MAX / sizeof(*plan->ptr)) - 1u) {
    return 0;
  }

  UF_long *ptr = (UF_long *)calloc((size_t)n + 1u, sizeof(*ptr));
  if (ptr == NULL) {
    return 0;
  }

  for (UF_long row = 0; row < n; ++row) {
    if (src_ptr[row] > src_ptr[row + 1u]) {
      free(ptr);
      return 0;
    }
    for (UF_long p = src_ptr[row]; p < src_ptr[row + 1u]; ++p) {
      const UF_long col = src_cols[p];
      if (col >= n ||
          (source_upper ? col <= row : col >= row)) {
        free(ptr);
        return 0;
      }
      ptr[col + 1u]++;
    }
  }
  for (UF_long row = 0; row < n; ++row) {
    ptr[row + 1u] += ptr[row];
  }

  const UF_long nnz = ptr[n];
  UF_long *cols = NULL;
  UF_long *source_pos = NULL;
  uint16_t *cols16 = NULL;
  uint16_t *source_pos16 = NULL;
  UF_long *next = NULL;
  if (nnz > 0u) {
    const int compact_serial = !build_parallel_metadata &&
      n <= (UF_long)UINT16_MAX + 1u && nnz <= (UF_long)UINT16_MAX;
    if (nnz > (UF_long)(SIZE_MAX / sizeof(*cols))) {
      free(ptr);
      return 0;
    }
    cols = compact_serial ? NULL
      : (UF_long *)malloc((size_t)nnz * sizeof(*cols));
    source_pos = compact_serial ? NULL
      : (UF_long *)malloc((size_t)nnz * sizeof(*source_pos));
    cols16 = compact_serial
      ? (uint16_t *)malloc((size_t)nnz * sizeof(*cols16)) : NULL;
    source_pos16 = compact_serial
      ? (uint16_t *)malloc((size_t)nnz * sizeof(*source_pos16)) : NULL;
    next = (UF_long *)malloc((size_t)n * sizeof(*next));
    if ((!compact_serial && (cols == NULL || source_pos == NULL)) ||
        (compact_serial && (cols16 == NULL || source_pos16 == NULL)) ||
        next == NULL) {
      free(ptr);
      free(cols);
      free(source_pos);
      free(cols16);
      free(source_pos16);
      free(next);
      return 0;
    }
    memcpy(next, ptr, (size_t)n * sizeof(*next));
    for (UF_long row = 0; row < n; ++row) {
      for (UF_long p = src_ptr[row]; p < src_ptr[row + 1u]; ++p) {
        const UF_long col = src_cols[p];
        const UF_long dst = next[col]++;
        if (compact_serial) {
          cols16[dst] = (uint16_t)row;
          source_pos16[dst] = (uint16_t)p;
        } else {
          cols[dst] = row;
          source_pos[dst] = p;
        }
      }
    }
    free(next);
  }

  plan->ptr = ptr;
  plan->cols = cols;
  plan->source_pos = source_pos;
  plan->cols16 = cols16;
  plan->source_pos16 = source_pos16;
  plan->source_upper = source_upper ? 1 : 0;
  plan->upper = source_upper ? 0 : 1;
  plan->diagonal = source_upper ? 1 : 0;
  plan->serial_only = build_parallel_metadata ? 0 : 1;

  if (nnz == 0u || !build_parallel_metadata) {
    return 1;
  }

  UF_long tail_entries = 0;
  const UF_long tail_start =
    kls_row_solve_dense_tail_start(n, plan->ptr, &tail_entries);
  plan->dense_tail_start = tail_start;
  plan->dense_tail_rows = tail_start < n ? n - tail_start : 0u;
  plan->dense_tail_entries = tail_entries;
  (void)kls_build_row_solve_sparse_levels(
    n, plan->ptr, plan->cols, tail_start, plan->upper,
    &plan->sparse_level_ptr, &plan->sparse_level_rows,
    &plan->sparse_level_count, &plan->sparse_level_max_width);
  plan->slice_bounds =
    kls_build_row_solve_slice_bounds(n, plan->ptr, tail_start,
                                     tail_entries,
                                     &plan->slice_max_entries);
  if (plan->slice_bounds != NULL) {
    (void)kls_build_row_solve_segment_splits(
      n, plan->ptr, plan->cols, plan->slice_bounds, plan->upper,
      &plan->segment_split, &plan->segmented_rows, &plan->rect_entries,
      &plan->tri_entries);
  }
  return 1;
}

static int kls_build_one_row_solve_transpose_plan(
  kls_solver *solver,
  int source_upper,
  kls_row_solve_transpose_plan *plan) {
  return kls_build_one_row_solve_transpose_plan_impl(
    solver, source_upper, plan, 1);
}

static void kls_build_row_solve_transpose_plans(kls_solver *solver) {
  if (solver == NULL) {
    return;
  }
  (void)kls_build_one_row_solve_transpose_plan(
    solver, 1, &solver->row_solve_ut_plan);
  (void)kls_build_one_row_solve_transpose_plan(
    solver, 0, &solver->row_solve_lt_plan);
}

static void kls_clear_row_solve_transpose_thread_bounds(
  kls_row_solve_transpose_plan *plan) {
  if (plan == NULL) {
    return;
  }
  free(plan->thread_bounds);
  plan->thread_bounds = NULL;
  plan->thread_count = 0;
  plan->thread_max_rect_entries = 0;
  plan->sparse_cluster_levels = 0;
}

static int kls_build_one_row_solve_transpose_thread_bounds(
  kls_row_solve_transpose_plan *plan,
  UF_long n,
  int thread_count) {
  if (plan == NULL || plan->ptr == NULL || thread_count < 2) {
    return 0;
  }
  const int have_dense =
    plan->slice_bounds != NULL && plan->segment_split != NULL;
  const int have_sparse =
    plan->sparse_level_ptr != NULL && plan->sparse_level_rows != NULL;
  if (!have_dense && !have_sparse) {
    return 0;
  }
  if (plan->thread_count == thread_count &&
      (!have_dense || plan->thread_bounds != NULL)) {
    return 1;
  }

  kls_clear_row_solve_transpose_thread_bounds(plan);
  UF_long *bounds = NULL;
  UF_long max_rect = 0;
  if (have_dense &&
      !kls_build_one_row_solve_thread_bounds(
        n, plan->ptr, plan->segment_split, plan->slice_bounds,
        plan->upper, thread_count, &bounds, &max_rect)) {
    free(bounds);
    return 0;
  }

  plan->thread_bounds = bounds;
  plan->thread_count = thread_count;
  plan->thread_max_rect_entries = max_rect;
  plan->sparse_cluster_levels =
    kls_row_solve_sparse_cluster_levels(plan->sparse_level_ptr,
                                        plan->sparse_level_count,
                                        thread_count);
  return 1;
}

static int kls_build_row_solve_transpose_thread_bounds(kls_solver *solver,
                                                       int thread_count) {
  if (solver == NULL || thread_count < 2 ||
      solver->row_solve_partition_ready == 0) {
    return 0;
  }
  const int have_ut =
    solver->row_solve_ut_plan.ptr != NULL &&
    (solver->row_solve_ut_plan.segment_split != NULL ||
     solver->row_solve_ut_plan.sparse_level_ptr != NULL);
  const int have_lt =
    solver->row_solve_lt_plan.ptr != NULL &&
    (solver->row_solve_lt_plan.segment_split != NULL ||
     solver->row_solve_lt_plan.sparse_level_ptr != NULL);
  if (!have_ut && !have_lt) {
    return 0;
  }
  int ok = 1;
  if (have_ut) {
    ok = ok &&
      kls_build_one_row_solve_transpose_thread_bounds(
        &solver->row_solve_ut_plan, solver->n, thread_count);
  }
  if (have_lt) {
    ok = ok &&
      kls_build_one_row_solve_transpose_thread_bounds(
        &solver->row_solve_lt_plan, solver->n, thread_count);
  }
  return ok;
}

static void kls_build_one_row_solve_partition(kls_solver *solver,
                                              int upper) {
  const UF_long n = solver->n;
  const UF_long *ptr =
    upper ? solver->row_refactor_u_ptr : solver->row_refactor_l_ptr;
  const UF_long *cols =
    upper ? solver->row_refactor_u_cols : solver->row_refactor_l_cols;
  UF_long entries = 0;
  const UF_long start = kls_row_solve_dense_tail_start(n, ptr, &entries);
  if (upper) {
    solver->row_solve_u_dense_tail_start = start;
    solver->row_solve_u_dense_tail_rows = start < n ? n - start : 0u;
    solver->row_solve_u_dense_tail_entries = entries;
    (void)kls_build_row_solve_sparse_levels(
      n, ptr, cols, start, 1, &solver->row_solve_u_sparse_level_ptr,
      &solver->row_solve_u_sparse_level_rows,
      &solver->row_solve_u_sparse_level_count,
      &solver->row_solve_u_sparse_level_max_width);
    solver->row_solve_u_slice_bounds =
      kls_build_row_solve_slice_bounds(
        n, ptr, start, entries, &solver->row_solve_u_slice_max_entries);
    if (solver->row_solve_u_slice_bounds != NULL) {
      (void)kls_build_row_solve_segment_splits(
        n, ptr, cols, solver->row_solve_u_slice_bounds, 1,
        &solver->row_solve_u_segment_split,
        &solver->row_solve_u_segmented_rows,
        &solver->row_solve_u_rect_entries,
        &solver->row_solve_u_tri_entries);
    }
  } else {
    solver->row_solve_l_dense_tail_start = start;
    solver->row_solve_l_dense_tail_rows = start < n ? n - start : 0u;
    solver->row_solve_l_dense_tail_entries = entries;
    (void)kls_build_row_solve_sparse_levels(
      n, ptr, cols, start, 0, &solver->row_solve_l_sparse_level_ptr,
      &solver->row_solve_l_sparse_level_rows,
      &solver->row_solve_l_sparse_level_count,
      &solver->row_solve_l_sparse_level_max_width);
    solver->row_solve_l_slice_bounds =
      kls_build_row_solve_slice_bounds(
        n, ptr, start, entries, &solver->row_solve_l_slice_max_entries);
    if (solver->row_solve_l_slice_bounds != NULL) {
      (void)kls_build_row_solve_segment_splits(
        n, ptr, cols, solver->row_solve_l_slice_bounds, 0,
        &solver->row_solve_l_segment_split,
        &solver->row_solve_l_segmented_rows,
        &solver->row_solve_l_rect_entries,
        &solver->row_solve_l_tri_entries);
    }
  }
}

typedef struct {
  kls_solver *solver;
  int kind;
  double seconds;
} kls_row_solve_plan_build_job;

static void *kls_row_solve_plan_build_main(void *arg) {
  kls_row_solve_plan_build_job *job = (kls_row_solve_plan_build_job *)arg;
  if (job == NULL || job->solver == NULL) {
    return NULL;
  }
  const double begin = kls_now_seconds();
  switch (job->kind) {
    case 0:
      kls_build_one_row_solve_partition(job->solver, 0);
      break;
    case 1:
      kls_build_one_row_solve_partition(job->solver, 1);
      break;
    case 2:
      (void)kls_build_one_row_solve_transpose_plan(
        job->solver, 1, &job->solver->row_solve_ut_plan);
      break;
    case 3:
      (void)kls_build_one_row_solve_transpose_plan(
        job->solver, 0, &job->solver->row_solve_lt_plan);
      break;
    default:
      break;
  }
  job->seconds = kls_now_seconds() - begin;
  return NULL;
}

static void kls_record_row_solve_partition(kls_solver *solver) {
  if (solver == NULL) {
    return;
  }
  kls_clear_row_solve_partition(solver);
  if (solver->row_refactor_pattern_n != solver->n ||
      solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_u_ptr == NULL) {
    return;
  }

  solver->row_solve_partition_ready = 1;
  solver->row_solve_partition_slices = KLS_ROW_SOLVE_TRAPEZOID_SLICES;
  const int defer_transpose =
    getenv("KLS_DEFER_ROW_TRANSPOSE_PLANS") != NULL;
  const int parallel =
    getenv("KLS_PARALLEL_ROW_SOLVE_PARTITION_BUILD") != NULL &&
    solver->options.threads >= 4 && solver->n >= 10000u;
  if (!parallel) {
    kls_build_one_row_solve_partition(solver, 0);
    kls_build_one_row_solve_partition(solver, 1);
    if (!defer_transpose) {
      kls_build_row_solve_transpose_plans(solver);
    }
    return;
  }

  kls_row_solve_plan_build_job jobs[4] = {
    { solver, 0, 0.0 }, { solver, 1, 0.0 },
    { solver, 2, 0.0 }, { solver, 3, 0.0 }
  };
  pthread_t threads[4];
  unsigned char launched[4] = { 0u, 0u, 0u, 0u };
  const int job_count = defer_transpose ? 2 : 4;
  for (int kind = 1; kind < job_count; ++kind) {
    if (pthread_create(&threads[kind], NULL,
                       kls_row_solve_plan_build_main, &jobs[kind]) == 0) {
      launched[kind] = 1u;
    } else {
      (void)kls_row_solve_plan_build_main(&jobs[kind]);
    }
  }
  (void)kls_row_solve_plan_build_main(&jobs[0]);
  for (int kind = 1; kind < job_count; ++kind) {
    if (launched[kind]) {
      (void)pthread_join(threads[kind], NULL);
    }
  }
  if (getenv("KLS_TRACE_ROW_SOLVE_PARTITION_PREP") != NULL) {
    fprintf(stderr,
            "KLS row solve partition plans L %.3fms U %.3fms "
            "UT %.3fms LT %.3fms\n",
            1e3 * jobs[0].seconds, 1e3 * jobs[1].seconds,
            1e3 * jobs[2].seconds, 1e3 * jobs[3].seconds);
  }
}

typedef struct {
  kls_solver *solver;
} kls_row_solve_partition_prep_job;

static void *kls_row_solve_partition_prep_main(void *arg) {
  kls_row_solve_partition_prep_job *job =
    (kls_row_solve_partition_prep_job *)arg;
  if (job != NULL) {
    kls_record_row_solve_partition(job->solver);
  }
  return NULL;
}

static int kls_build_row_solve_pattern(kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      solver->symbolic->nblocks == 0u || solver->symbolic->R == NULL) {
    return 0;
  }
  if (kls_row_solve_pattern_is_ready(solver)) {
    kls_record_row_solve_partition(solver);
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
  UF_long *next = n > 0u
    ? (UF_long *)malloc((size_t)n * sizeof(*next)) : NULL;
  UF_long *l_cols = NULL;
  double **l_values = NULL;
  double *l_row_values = NULL;
  UF_long *u_cols = NULL;
  double **u_values = NULL;
  double *u_row_values = NULL;
  if (l_ptr == NULL || u_ptr == NULL || (n > 0u && next == NULL)) {
    free(l_ptr);
    free(u_ptr);
    free(next);
    return 0;
  }

  for (UF_long block = 0; block < solver->symbolic->nblocks; ++block) {
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    if (k1 > k2 || k2 > n) {
      goto fail;
    }
    const UF_long nk = k2 - k1;
    for (UF_long local_j = 0; local_j < nk; ++local_j) {
      const UF_long j = k1 + local_j;
      UF_long *li = solver->refactor_l_indices[j];
      const UF_long len = kls_numeric_l_len_in_block(solver, k1, k2, j);
      if (len > 0u && li == NULL) {
        goto fail;
      }
      for (UF_long p = 0; p < len; ++p) {
        const UF_long local_row = li[p];
        const UF_long row = k1 + local_row;
        if (local_row <= local_j || local_row >= nk || row >= n) {
          goto fail;
        }
        l_ptr[row + 1u]++;
      }
    }
  }
  for (UF_long block = 0; block < solver->symbolic->nblocks; ++block) {
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    const UF_long nk = k2 - k1;
    for (UF_long local_k = 0; local_k < nk; ++local_k) {
      const UF_long k = k1 + local_k;
      UF_long *ui = solver->refactor_u_indices[k];
      const UF_long len = kls_numeric_u_len_in_block(solver, k1, k2, k);
      if (len > 0u && ui == NULL) {
        goto fail;
      }
      for (UF_long p = 0; p < len; ++p) {
        const UF_long local_row = ui[p];
        const UF_long row = k1 + local_row;
        if (local_row >= local_k || local_row >= nk || row >= n) {
          goto fail;
        }
        u_ptr[row + 1u]++;
      }
    }
  }

  for (UF_long i = 0; i < n; ++i) {
    l_ptr[i + 1u] += l_ptr[i];
    u_ptr[i + 1u] += u_ptr[i];
  }
  const UF_long l_nnz = l_ptr[n];
  const UF_long u_nnz = u_ptr[n];
  l_cols = l_nnz > 0u
    ? (UF_long *)malloc((size_t)l_nnz * sizeof(*l_cols)) : NULL;
  l_values = l_nnz > 0u
    ? (double **)malloc((size_t)l_nnz * sizeof(*l_values)) : NULL;
  l_row_values = l_nnz > 0u
    ? (double *)malloc((size_t)l_nnz * sizeof(*l_row_values)) : NULL;
  u_cols = u_nnz > 0u
    ? (UF_long *)malloc((size_t)u_nnz * sizeof(*u_cols)) : NULL;
  u_values = u_nnz > 0u
    ? (double **)malloc((size_t)u_nnz * sizeof(*u_values)) : NULL;
  u_row_values = kls_aligned_double_values(u_nnz);
  if ((l_nnz > 0u &&
       (l_cols == NULL || l_values == NULL || l_row_values == NULL)) ||
      (u_nnz > 0u &&
       (u_cols == NULL || u_values == NULL || u_row_values == NULL))) {
    goto fail;
  }

  if (n > 0u) {
    memcpy(next, l_ptr, (size_t)n * sizeof(*next));
  }
  for (UF_long block = 0; block < solver->symbolic->nblocks; ++block) {
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    const UF_long nk = k2 - k1;
    for (UF_long local_j = 0; local_j < nk; ++local_j) {
      const UF_long j = k1 + local_j;
      UF_long *li = solver->refactor_l_indices[j];
      double *lx = solver->refactor_l_values[j];
      const UF_long len = kls_numeric_l_len_in_block(solver, k1, k2, j);
      for (UF_long p = 0; p < len; ++p) {
        const UF_long row = k1 + li[p];
        const UF_long dst = next[row]++;
        l_cols[dst] = j;
        l_values[dst] = lx + p;
      }
    }
  }
  if (n > 0u) {
    memcpy(next, u_ptr, (size_t)n * sizeof(*next));
  }
  for (UF_long block = 0; block < solver->symbolic->nblocks; ++block) {
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    const UF_long nk = k2 - k1;
    for (UF_long local_k = 0; local_k < nk; ++local_k) {
      const UF_long k = k1 + local_k;
      UF_long *ui = solver->refactor_u_indices[k];
      double *ux = solver->refactor_u_values[k];
      const UF_long len = kls_numeric_u_len_in_block(solver, k1, k2, k);
      for (UF_long p = 0; p < len; ++p) {
        const UF_long row = k1 + ui[p];
        const UF_long dst = next[row]++;
        u_cols[dst] = k;
        u_values[dst] = ux + p;
      }
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
  solver->row_refactor_pattern_n = n;
  kls_record_row_solve_partition(solver);
  return 1;

fail:
  free(l_ptr);
  free(u_ptr);
  free(next);
  free(l_cols);
  free(l_values);
  free(l_row_values);
  free(u_cols);
  free(u_values);
  free(u_row_values);
  return 0;
}

static int kls_build_row_refactor_group_successors(
  UF_long group_count,
  const UF_long *row_group,
  const UF_long *group_dep_ptr,
  const UF_long *group_dep_rows,
  UF_long **successor_ptr_out,
  UF_long **successor_groups_out,
  UF_long **pred_count_out,
  UF_long **root_groups_out,
  UF_long *edge_count_out,
  UF_long *root_count_out,
  UF_long *leaf_count_out,
  UF_long *max_fanout_out) {
  KLS_SET_OPTIONAL_OUTPUT(successor_ptr_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(successor_groups_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(pred_count_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(root_groups_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(edge_count_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(root_count_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(leaf_count_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(max_fanout_out, 0);
  if (group_count < 0 ||
      (uint64_t)group_count >=
        (uint64_t)(SIZE_MAX / sizeof(UF_long)) ||
      successor_ptr_out == NULL || successor_groups_out == NULL ||
      pred_count_out == NULL || root_groups_out == NULL ||
      edge_count_out == NULL || root_count_out == NULL ||
      leaf_count_out == NULL || max_fanout_out == NULL ||
      (group_count > 0u &&
       (row_group == NULL || group_dep_ptr == NULL))) {
    return 0;
  }
  if (group_count == 0u) {
    return 1;
  }

  UF_long *successor_ptr =
    (UF_long *)calloc((size_t)group_count + 1u, sizeof(*successor_ptr));
  UF_long *indegree =
    (UF_long *)calloc((size_t)group_count, sizeof(*indegree));
  UF_long *marks =
    (UF_long *)malloc((size_t)group_count * sizeof(*marks));
  if (successor_ptr == NULL || indegree == NULL || marks == NULL) {
    free(successor_ptr);
    free(indegree);
    free(marks);
    return 0;
  }
  for (UF_long g = 0; g < group_count; ++g) {
    marks[g] = KLS_KLU_EMPTY;
  }

  for (UF_long g = 0; g < group_count; ++g) {
    const UF_long begin = group_dep_ptr[g];
    const UF_long end = group_dep_ptr[g + 1u];
    if (begin > end || (end > begin && group_dep_rows == NULL)) {
      free(successor_ptr);
      free(indegree);
      free(marks);
      return 0;
    }
    for (UF_long pos = begin; pos < end; ++pos) {
      const UF_long dep_group = row_group[group_dep_rows[pos]];
      if (dep_group >= group_count || dep_group == g ||
          marks[dep_group] == g) {
        continue;
      }
      marks[dep_group] = g;
      successor_ptr[dep_group + 1u]++;
      indegree[g]++;
    }
  }
  for (UF_long g = 0; g < group_count; ++g) {
    successor_ptr[g + 1u] += successor_ptr[g];
    marks[g] = KLS_KLU_EMPTY;
  }

  const UF_long edge_count = successor_ptr[group_count];
  UF_long *successor_groups = edge_count > 0u
    ? (UF_long *)malloc((size_t)edge_count * sizeof(*successor_groups))
    : NULL;
  UF_long *write_pos =
    (UF_long *)malloc((size_t)group_count * sizeof(*write_pos));
  if ((edge_count > 0u && successor_groups == NULL) || write_pos == NULL) {
    free(successor_ptr);
    free(successor_groups);
    free(write_pos);
    free(indegree);
    free(marks);
    return 0;
  }
  memcpy(write_pos, successor_ptr, (size_t)group_count * sizeof(*write_pos));
  for (UF_long g = 0; g < group_count; ++g) {
    const UF_long begin = group_dep_ptr[g];
    const UF_long end = group_dep_ptr[g + 1u];
    for (UF_long pos = begin; pos < end; ++pos) {
      const UF_long dep_group = row_group[group_dep_rows[pos]];
      if (dep_group >= group_count || dep_group == g ||
          marks[dep_group] == g) {
        continue;
      }
      marks[dep_group] = g;
      successor_groups[write_pos[dep_group]++] = g;
    }
    for (UF_long pos = begin; pos < end; ++pos) {
      const UF_long dep_group = row_group[group_dep_rows[pos]];
      if (dep_group < group_count) {
        marks[dep_group] = KLS_KLU_EMPTY;
      }
    }
  }
  free(write_pos);
  free(marks);

  UF_long root_count = 0;
  UF_long leaf_count = 0;
  UF_long max_fanout = 0;
  for (UF_long g = 0; g < group_count; ++g) {
    if (indegree[g] == 0u) {
      root_count++;
    }
    const UF_long fanout = successor_ptr[g + 1u] - successor_ptr[g];
    if (fanout == 0u) {
      leaf_count++;
    }
    if (fanout > max_fanout) {
      max_fanout = fanout;
    }
  }
  UF_long *root_groups = root_count > 0u
    ? (UF_long *)malloc((size_t)root_count * sizeof(*root_groups))
    : NULL;
  if (root_count > 0u && root_groups == NULL) {
    free(successor_ptr);
    free(successor_groups);
    free(indegree);
    return 0;
  }
  UF_long root_pos = 0;
  for (UF_long g = 0; g < group_count; ++g) {
    if (indegree[g] == 0u) {
      if (root_pos >= root_count) {
        free(successor_ptr);
        free(successor_groups);
        free(indegree);
        free(root_groups);
        return 0;
      }
      root_groups[root_pos++] = g;
    }
  }

  *successor_ptr_out = successor_ptr;
  *successor_groups_out = successor_groups;
  *pred_count_out = indegree;
  *root_groups_out = root_groups;
  *edge_count_out = edge_count;
  *root_count_out = root_count;
  *leaf_count_out = leaf_count;
  *max_fanout_out = max_fanout;
  return 1;
}

typedef struct {
  UF_long group_count;
  const UF_long *row_group;
  const UF_long *group_dep_ptr;
  const UF_long *group_dep_rows;
  UF_long *successor_ptr;
  UF_long *successor_groups;
  UF_long *pred_count;
  UF_long *root_groups;
  UF_long edge_count;
  UF_long root_count;
  UF_long leaf_count;
  UF_long max_fanout;
  double seconds;
  int ok;
} kls_group_successor_prep_job;

static void *kls_group_successor_prep_main(void *arg) {
  kls_group_successor_prep_job *job =
    (kls_group_successor_prep_job *)arg;
  if (job == NULL) {
    return NULL;
  }
  const double begin = kls_now_seconds();
  job->ok = kls_build_row_refactor_group_successors(
    job->group_count, job->row_group, job->group_dep_ptr,
    job->group_dep_rows, &job->successor_ptr, &job->successor_groups,
    &job->pred_count, &job->root_groups, &job->edge_count,
    &job->root_count, &job->leaf_count, &job->max_fanout);
  job->seconds = kls_now_seconds() - begin;
  return NULL;
}

static void kls_free_group_successor_prep_job(
  kls_group_successor_prep_job *job) {
  if (job == NULL) {
    return;
  }
  free(job->successor_ptr);
  free(job->successor_groups);
  free(job->pred_count);
  free(job->root_groups);
  job->successor_ptr = NULL;
  job->successor_groups = NULL;
  job->pred_count = NULL;
  job->root_groups = NULL;
}

static int kls_row_refactor_u_row_extends_next(const UF_long *u_ptr,
                                               const UF_long *u_cols,
                                               UF_long n,
                                               UF_long prev,
                                               UF_long curr) {
  if (u_ptr == NULL || u_cols == NULL || prev >= n || curr >= n ||
      prev + 1u != curr) {
    return 0;
  }
  const UF_long prev_begin = u_ptr[prev];
  const UF_long prev_end = u_ptr[prev + 1u];
  const UF_long curr_begin = u_ptr[curr];
  const UF_long curr_end = u_ptr[curr + 1u];
  if (prev_begin > prev_end || curr_begin > curr_end) {
    return 0;
  }
  const UF_long prev_len = prev_end - prev_begin;
  const UF_long curr_len = curr_end - curr_begin;
  if (prev_len != curr_len + 1u || prev_len == 0u ||
      u_cols[prev_begin] != curr) {
    return 0;
  }
  return curr_len == 0u ||
         memcmp(u_cols + prev_begin + 1u, u_cols + curr_begin,
                (size_t)curr_len * sizeof(*u_cols)) == 0;
}

static int kls_row_refactor_row_depends_in_range(const UF_long *l_ptr,
                                                 const UF_long *l_cols,
                                                 UF_long row,
                                                 UF_long begin,
                                                 UF_long end) {
  if (l_ptr == NULL || begin >= end) {
    return 0;
  }
  const UF_long l_begin = l_ptr[row];
  const UF_long l_end = l_ptr[row + 1u];
  if (l_begin > l_end || (l_begin < l_end && l_cols == NULL)) {
    return 1;
  }
  for (UF_long p = l_begin; p < l_end; ++p) {
    const UF_long dep = l_cols[p];
    if (dep >= begin && dep < end) {
      return 1;
    }
  }
  return 0;
}

static int kls_row_refactor_group_has_internal_deps(const UF_long *l_ptr,
                                                    const UF_long *l_cols,
                                                    UF_long begin,
                                                    UF_long end) {
  if (begin >= end) {
    return 0;
  }
  for (UF_long row = begin; row < end; ++row) {
    if (kls_row_refactor_row_depends_in_range(l_ptr, l_cols, row,
                                              begin, row)) {
      return 1;
    }
  }
  return 0;
}

static void kls_try_build_row_refactor_compact_panel_storage(
  UF_long group_count,
  const UF_long *group_ptr,
  const unsigned char *group_kind,
  const UF_long *group_trailing_len,
  UF_long **begin_out,
  double **values_out,
  UF_long *groups_out,
  UF_long *entries_out) {
  KLS_SET_OPTIONAL_OUTPUT(begin_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(values_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(groups_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(entries_out, 0);
  if (group_count == 0u || group_ptr == NULL || group_kind == NULL ||
      group_trailing_len == NULL || begin_out == NULL || values_out == NULL ||
      groups_out == NULL || entries_out == NULL ||
      group_count >= (UF_long)(SIZE_MAX / sizeof(UF_long))) {
    return;
  }

  UF_long *begin =
    (UF_long *)calloc((size_t)group_count + 1u, sizeof(*begin));
  if (begin == NULL) {
    return;
  }

  UF_long total = 0;
  UF_long groups = 0;
  for (UF_long group = 0; group < group_count; ++group) {
    if (group_kind[group] != KLS_ROW_REFACTOR_GROUP_DENSE) {
      begin[group + 1u] = total;
      continue;
    }
    const UF_long width = group_ptr[group + 1u] - group_ptr[group];
    const UF_long trailing_len = group_trailing_len[group];
    if (!kls_row_refactor_prefers_compact_dense_panel(width, trailing_len)) {
      begin[group + 1u] = total;
      continue;
    }
    if (width == 0u || width > UF_long_max / width) {
      free(begin);
      return;
    }
    const UF_long dense_entries = width * width;
    if (trailing_len > 0u &&
        width > (UF_long_max - dense_entries) / trailing_len) {
      free(begin);
      return;
    }
    const UF_long entries = dense_entries + width * trailing_len;
    if (entries > UF_long_max - total) {
      free(begin);
      return;
    }
    total += entries;
    begin[group + 1u] = total;
    groups++;
  }
  if (total == 0u || total > (UF_long)(SIZE_MAX / sizeof(double))) {
    free(begin);
    return;
  }

  double *values = (double *)malloc((size_t)total * sizeof(*values));
  if (values == NULL) {
    free(begin);
    return;
  }
  *begin_out = begin;
  *values_out = values;
  *groups_out = groups;
  *entries_out = total;
}

static int kls_row_refactor_group_has_compact_panel(
  UF_long group_count,
  const UF_long *group_compact_panel_begin,
  UF_long group) {
  return group_compact_panel_begin != NULL && group < group_count &&
         group_compact_panel_begin[group + 1u] >
           group_compact_panel_begin[group];
}

static int kls_row_refactor_dense_producer_run_start(
  UF_long n,
  const UF_long *l_cols,
  const UF_long *row_group,
  UF_long group_count,
  const UF_long *group_ptr,
  const unsigned char *group_kind,
  const UF_long *group_compact_panel_begin,
  UF_long row,
  UF_long row_group_id,
  UF_long p,
  UF_long l_end,
  UF_long *group_out,
  UF_long *len_out,
  UF_long *suffix_out) {
  if (group_out != NULL) {
    *group_out = KLS_KLU_EMPTY;
  }
  KLS_SET_OPTIONAL_OUTPUT(len_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(suffix_out, 0);
  if (l_cols == NULL || row_group == NULL || group_ptr == NULL ||
      group_kind == NULL || p >= l_end || row >= n) {
    return 0;
  }
  const UF_long dep0 = l_cols[p];
  if (dep0 >= row || dep0 >= n) {
    return 0;
  }
  const UF_long dep_group = row_group[dep0];
  if (dep_group >= group_count || dep_group == row_group_id ||
      group_kind[dep_group] != KLS_ROW_REFACTOR_GROUP_DENSE ||
      !kls_row_refactor_group_has_compact_panel(
        group_count, group_compact_panel_begin, dep_group)) {
    return 0;
  }
  const UF_long group_begin = group_ptr[dep_group];
  const UF_long group_end = group_ptr[dep_group + 1u];
  if (group_begin >= group_end || dep0 < group_begin || dep0 >= group_end ||
      group_end > row) {
    return 0;
  }

  UF_long run_len = 0;
  UF_long expected_dep = dep0;
  while (p + run_len < l_end && expected_dep < group_end) {
    const UF_long dep = l_cols[p + run_len];
    if (dep != expected_dep || dep >= row ||
        row_group[dep] != dep_group) {
      break;
    }
    run_len++;
    expected_dep++;
  }
  if (run_len < 2u) {
    return 0;
  }
  if (group_out != NULL) {
    *group_out = dep_group;
  }
  if (len_out != NULL) {
    *len_out = run_len;
  }
  if (suffix_out != NULL) {
    *suffix_out = dep0 - group_begin;
  }
  return 1;
}

typedef struct {
  UF_long n;
  const UF_long *l_ptr;
  const UF_long *l_cols;
  const UF_long *row_group;
  UF_long group_count;
  const UF_long *group_ptr;
  const unsigned char *group_kind;
  const UF_long *group_compact_panel_begin;
  UF_long row_begin;
  UF_long row_end;
  UF_long *run_ptr;
  UF_long *run_l_begin;
  UF_long *run_group;
  UF_long *run_len;
  UF_long run_count;
  UF_long run_rows;
  UF_long run_dep_rows;
  UF_long run_max_per_row;
  UF_long full_suffix_run_count;
  UF_long full_suffix_rows;
  UF_long multi_run_rows;
  UF_long fragmented_rows;
  int fill;
  int ok;
} kls_dense_producer_run_build_job;

static void *kls_dense_producer_run_build_main(void *arg) {
  kls_dense_producer_run_build_job *job =
    (kls_dense_producer_run_build_job *)arg;
  if (job == NULL || job->l_ptr == NULL || job->l_cols == NULL ||
      job->row_group == NULL || job->group_ptr == NULL ||
      job->group_kind == NULL || job->group_compact_panel_begin == NULL ||
      job->run_ptr == NULL) {
    return NULL;
  }
  job->ok = 1;
  for (UF_long row = job->row_begin; row < job->row_end; ++row) {
    const UF_long row_group_id = job->row_group[row];
    const UF_long l_begin = job->l_ptr[row];
    const UF_long l_end = job->l_ptr[row + 1u];
    const UF_long row_l_len = l_end - l_begin;
    UF_long row_runs = 0u;
    UF_long covered_entries = 0u;
    UF_long row_full_suffix_runs = 0u;
    UF_long write = job->fill ? job->run_ptr[row] : 0u;
    UF_long p = l_begin;
    while (p < l_end) {
      UF_long dep_group = KLS_KLU_EMPTY;
      UF_long run_len = 0u;
      UF_long suffix = 0u;
      if (kls_row_refactor_dense_producer_run_start(
            job->n, job->l_cols, job->row_group, job->group_count,
            job->group_ptr, job->group_kind,
            job->group_compact_panel_begin, row, row_group_id, p, l_end,
            &dep_group, &run_len, &suffix)) {
        if (job->fill) {
          if (write >= job->run_ptr[row + 1u] ||
              job->run_l_begin == NULL || job->run_group == NULL ||
              job->run_len == NULL) {
            job->ok = 0;
            return NULL;
          }
          job->run_l_begin[write] = p;
          job->run_group[write] = dep_group;
          job->run_len[write] = run_len;
          write++;
        } else {
          const UF_long width =
            job->group_ptr[dep_group + 1u] - job->group_ptr[dep_group];
          row_runs++;
          job->run_count++;
          covered_entries += run_len;
          job->run_dep_rows += run_len;
          if (suffix + run_len == width) {
            row_full_suffix_runs++;
            job->full_suffix_run_count++;
          }
        }
        p += run_len;
      } else {
        p++;
      }
    }
    if (job->fill) {
      if (write != job->run_ptr[row + 1u]) {
        job->ok = 0;
        return NULL;
      }
      continue;
    }
    job->run_ptr[row + 1u] = row_runs;
    if (row_runs > 0u) {
      job->run_rows++;
      if (row_runs > job->run_max_per_row) {
        job->run_max_per_row = row_runs;
      }
      if (row_full_suffix_runs > 0u) {
        job->full_suffix_rows++;
      }
      if (row_runs > 1u) {
        job->multi_run_rows++;
      }
      if (covered_entries != row_l_len ||
          row_full_suffix_runs != row_runs) {
        job->fragmented_rows++;
      }
    }
  }
  return NULL;
}

static int kls_run_dense_producer_run_build_jobs(
  kls_dense_producer_run_build_job *jobs,
  int job_count) {
  if (jobs == NULL || job_count < 1 || job_count > 8) {
    return 0;
  }
  pthread_t threads[8];
  unsigned char launched[8] = { 0u };
  for (int tid = 1; tid < job_count; ++tid) {
    if (pthread_create(&threads[tid], NULL,
                       kls_dense_producer_run_build_main,
                       &jobs[tid]) == 0) {
      launched[tid] = 1u;
    } else {
      (void)kls_dense_producer_run_build_main(&jobs[tid]);
    }
  }
  (void)kls_dense_producer_run_build_main(&jobs[0]);
  int ok = jobs[0].ok;
  for (int tid = 1; tid < job_count; ++tid) {
    if (launched[tid]) {
      (void)pthread_join(threads[tid], NULL);
    }
    ok = ok && jobs[tid].ok;
  }
  return ok;
}

static UF_long kls_dense_producer_run_row_boundary(
  const UF_long *l_ptr,
  UF_long n,
  UF_long target) {
  UF_long lo = 0u;
  UF_long hi = n;
  while (lo < hi) {
    const UF_long mid = lo + (hi - lo) / 2u;
    if (l_ptr[mid] < target) {
      lo = mid + 1u;
    } else {
      hi = mid;
    }
  }
  return lo;
}

static void kls_try_build_row_refactor_dense_producer_runs(
  UF_long n,
  const UF_long *l_ptr,
  const UF_long *l_cols,
  const UF_long *row_group,
  UF_long group_count,
  const UF_long *group_ptr,
  const unsigned char *group_kind,
  const UF_long *group_compact_panel_begin,
  UF_long **run_ptr_out,
  UF_long **run_l_begin_out,
  UF_long **run_group_out,
  UF_long **run_len_out,
  UF_long *run_count_out,
  UF_long *run_rows_out,
  UF_long *run_dep_rows_out,
  UF_long *run_max_per_row_out,
  UF_long *full_suffix_run_count_out,
  UF_long *full_suffix_rows_out,
  UF_long *multi_run_rows_out,
  UF_long *fragmented_rows_out,
  int thread_count) {
  KLS_SET_OPTIONAL_OUTPUT(run_ptr_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(run_l_begin_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(run_group_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(run_len_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(run_count_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(run_rows_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(run_dep_rows_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(run_max_per_row_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(full_suffix_run_count_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(full_suffix_rows_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(multi_run_rows_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(fragmented_rows_out, 0);
  if (n == 0u || l_ptr == NULL || l_cols == NULL || row_group == NULL ||
      group_count == 0u || group_ptr == NULL || group_kind == NULL ||
      group_compact_panel_begin == NULL || run_ptr_out == NULL ||
      run_l_begin_out == NULL || run_group_out == NULL ||
      run_len_out == NULL ||
      n > (UF_long)((SIZE_MAX / sizeof(UF_long)) - 1u)) {
    return;
  }

  UF_long *run_ptr = (UF_long *)calloc((size_t)n + 1u, sizeof(*run_ptr));
  if (run_ptr == NULL) {
    return;
  }

  UF_long run_count = 0;
  UF_long run_rows = 0;
  UF_long run_dep_rows = 0;
  UF_long run_max_per_row = 0;
  UF_long full_suffix_run_count = 0;
  UF_long full_suffix_rows = 0;
  UF_long multi_run_rows = 0;
  UF_long fragmented_rows = 0;
  int run_job_count = 1;
  if (thread_count >= 2 && n >= 10000u) {
    run_job_count = thread_count;
    if (run_job_count > 8) {
      run_job_count = 8;
    }
  }
  kls_dense_producer_run_build_job run_jobs[8];
  memset(run_jobs, 0, sizeof(run_jobs));
  const UF_long l_nnz = l_ptr[n];
  for (int tid = 0; tid < run_job_count; ++tid) {
    const UF_long begin_target =
      (l_nnz / (UF_long)run_job_count) * (UF_long)tid +
      ((l_nnz % (UF_long)run_job_count) * (UF_long)tid) /
        (UF_long)run_job_count;
    const UF_long end_target =
      (l_nnz / (UF_long)run_job_count) * (UF_long)(tid + 1) +
      ((l_nnz % (UF_long)run_job_count) * (UF_long)(tid + 1)) /
        (UF_long)run_job_count;
    const UF_long row_begin = tid == 0
      ? 0u : kls_dense_producer_run_row_boundary(l_ptr, n, begin_target);
    const UF_long row_end = tid + 1 == run_job_count
      ? n : kls_dense_producer_run_row_boundary(l_ptr, n, end_target);
    run_jobs[tid] = (kls_dense_producer_run_build_job){
      .n = n,
      .l_ptr = l_ptr,
      .l_cols = l_cols,
      .row_group = row_group,
      .group_count = group_count,
      .group_ptr = group_ptr,
      .group_kind = group_kind,
      .group_compact_panel_begin = group_compact_panel_begin,
      .row_begin = row_begin,
      .row_end = row_end,
      .run_ptr = run_ptr
    };
  }
  if (!kls_run_dense_producer_run_build_jobs(run_jobs, run_job_count)) {
    free(run_ptr);
    return;
  }
  for (int tid = 0; tid < run_job_count; ++tid) {
    if (run_count > UF_long_max - run_jobs[tid].run_count ||
        run_dep_rows > UF_long_max - run_jobs[tid].run_dep_rows) {
      free(run_ptr);
      return;
    }
    run_count += run_jobs[tid].run_count;
    run_rows += run_jobs[tid].run_rows;
    run_dep_rows += run_jobs[tid].run_dep_rows;
    full_suffix_run_count += run_jobs[tid].full_suffix_run_count;
    full_suffix_rows += run_jobs[tid].full_suffix_rows;
    multi_run_rows += run_jobs[tid].multi_run_rows;
    fragmented_rows += run_jobs[tid].fragmented_rows;
    if (run_jobs[tid].run_max_per_row > run_max_per_row) {
      run_max_per_row = run_jobs[tid].run_max_per_row;
    }
  }

  if (run_count == 0u) {
    free(run_ptr);
    return;
  }
  if (run_count > (UF_long)(SIZE_MAX / sizeof(UF_long))) {
    free(run_ptr);
    return;
  }
  UF_long prefix = 0;
  for (UF_long row = 0; row < n; ++row) {
    const UF_long count = run_ptr[row + 1u];
    run_ptr[row] = prefix;
    prefix += count;
  }
  run_ptr[n] = prefix;
  if (prefix != run_count) {
    free(run_ptr);
    return;
  }

  UF_long *run_l_begin =
    (UF_long *)malloc((size_t)run_count * sizeof(*run_l_begin));
  UF_long *run_group =
    (UF_long *)malloc((size_t)run_count * sizeof(*run_group));
  UF_long *run_len =
    (UF_long *)malloc((size_t)run_count * sizeof(*run_len));
  if (run_l_begin == NULL || run_group == NULL || run_len == NULL) {
    free(run_ptr);
    free(run_l_begin);
    free(run_group);
    free(run_len);
    return;
  }

  for (int tid = 0; tid < run_job_count; ++tid) {
    run_jobs[tid].run_l_begin = run_l_begin;
    run_jobs[tid].run_group = run_group;
    run_jobs[tid].run_len = run_len;
    run_jobs[tid].fill = 1;
    run_jobs[tid].ok = 0;
  }
  if (!kls_run_dense_producer_run_build_jobs(run_jobs, run_job_count)) {
    free(run_ptr);
    free(run_l_begin);
    free(run_group);
    free(run_len);
    return;
  }

  *run_ptr_out = run_ptr;
  *run_l_begin_out = run_l_begin;
  *run_group_out = run_group;
  *run_len_out = run_len;
  if (run_count_out != NULL) {
    *run_count_out = run_count;
  }
  if (run_rows_out != NULL) {
    *run_rows_out = run_rows;
  }
  if (run_dep_rows_out != NULL) {
    *run_dep_rows_out = run_dep_rows;
  }
  if (run_max_per_row_out != NULL) {
    *run_max_per_row_out = run_max_per_row;
  }
  if (full_suffix_run_count_out != NULL) {
    *full_suffix_run_count_out = full_suffix_run_count;
  }
  if (full_suffix_rows_out != NULL) {
    *full_suffix_rows_out = full_suffix_rows;
  }
  if (multi_run_rows_out != NULL) {
    *multi_run_rows_out = multi_run_rows;
  }
  if (fragmented_rows_out != NULL) {
    *fragmented_rows_out = fragmented_rows;
  }
}

typedef struct {
  UF_long n;
  const UF_long *l_ptr;
  const UF_long *l_cols;
  const UF_long *l_internal_ptr;
  const UF_long *u_cols;
  const UF_long *row_group;
  UF_long group_count;
  const UF_long *group_ptr;
  const unsigned char *group_kind;
  const UF_long *group_trailing_begin;
  const UF_long *group_trailing_len;
  const UF_long *run_ptr;
  const UF_long *run_l_begin;
  const UF_long *run_group;
  const UF_long *run_len;
  const UF_long *target_ptr;
  unsigned char *target_kind;
  UF_long *target_pos;
  UF_long row_begin;
  UF_long row_end;
  UF_long none_count;
  UF_long external_count;
  UF_long dense_count;
  UF_long pivot_count;
  UF_long trailing_count;
  int ok;
} kls_dense_producer_target_job;

static void *kls_dense_producer_target_main(void *arg) {
  kls_dense_producer_target_job *job =
    (kls_dense_producer_target_job *)arg;
  if (job == NULL) {
    return NULL;
  }
  job->ok = 1;
  UF_long none_count = 0;
  UF_long external_count = 0;
  UF_long dense_count = 0;
  UF_long pivot_count = 0;
  UF_long trailing_count = 0;
  for (UF_long row = job->row_begin; row < job->row_end; ++row) {
    const UF_long current_group = job->row_group[row];
    if (current_group >= job->group_count) {
      job->ok = 0;
      return NULL;
    }
    const UF_long row_begin = job->group_ptr[current_group];
    const UF_long row_end = job->group_ptr[current_group + 1u];
    const UF_long l_begin = job->l_ptr[row];
    const int current_is_dense =
      job->group_kind[current_group] == KLS_ROW_REFACTOR_GROUP_DENSE;
    const UF_long l_end = job->l_ptr[row + 1u];
    const UF_long l_external_end =
      current_is_dense ? job->l_internal_ptr[row] : l_begin;
    const UF_long external_len =
      l_begin <= l_external_end && l_external_end <= l_end
        ? l_external_end - l_begin : 0u;
    const UF_long current_trailing_begin =
      job->group_trailing_begin[current_group];
    const UF_long current_trailing_len =
      job->group_trailing_len[current_group];

    for (UF_long run = job->run_ptr[row];
         run < job->run_ptr[row + 1u]; ++run) {
      const UF_long dep_group = job->run_group[run];
      if (dep_group >= job->group_count) {
        job->ok = 0;
        return NULL;
      }
      const UF_long dep_trailing_begin =
        job->group_trailing_begin[dep_group];
      const UF_long dep_trailing_len = job->group_trailing_len[dep_group];
      if (job->target_ptr[run + 1u] - job->target_ptr[run] !=
            dep_trailing_len ||
          job->run_l_begin[run] < l_begin) {
        job->ok = 0;
        return NULL;
      }
      const UF_long run_offset = job->run_l_begin[run] - l_begin;
      if (job->run_len[run] > UF_long_max - run_offset) {
        job->ok = 0;
        return NULL;
      }
      const UF_long run_end = run_offset + job->run_len[run];
      for (UF_long dep_col_pos = 0; dep_col_pos < dep_trailing_len;
           ++dep_col_pos) {
        const UF_long target_index = job->target_ptr[run] + dep_col_pos;
        const UF_long col = job->u_cols[dep_trailing_begin + dep_col_pos];
        unsigned char kind = KLS_FRAGMENTED_UPDATE_TARGET_NONE;
        UF_long pos = 0;
        if (current_is_dense && col < row_begin) {
          UF_long external_pos = KLS_KLU_EMPTY;
          if (kls_find_uflong_sorted(job->l_cols + l_begin, external_len,
                                     col, &external_pos) &&
              external_pos >= run_end) {
            kind = KLS_FRAGMENTED_UPDATE_TARGET_EXTERNAL;
            pos = external_pos;
          }
        } else if (current_is_dense && col >= row_begin && col < row_end) {
          pos = col - row_begin;
          kind = col == row ? KLS_FRAGMENTED_UPDATE_TARGET_PIVOT
                            : KLS_FRAGMENTED_UPDATE_TARGET_DENSE;
        } else if (current_is_dense && current_trailing_len > 0u) {
          UF_long trailing_pos = KLS_KLU_EMPTY;
          if (kls_find_uflong_sorted(
                job->u_cols + current_trailing_begin, current_trailing_len,
                col, &trailing_pos)) {
            kind = KLS_FRAGMENTED_UPDATE_TARGET_TRAILING;
            pos = trailing_pos;
          }
        }
        job->target_kind[target_index] = kind;
        job->target_pos[target_index] = pos;
        switch ((kls_fragmented_update_target_kind)kind) {
          case KLS_FRAGMENTED_UPDATE_TARGET_EXTERNAL:
            external_count++;
            break;
          case KLS_FRAGMENTED_UPDATE_TARGET_DENSE:
            dense_count++;
            break;
          case KLS_FRAGMENTED_UPDATE_TARGET_PIVOT:
            pivot_count++;
            break;
          case KLS_FRAGMENTED_UPDATE_TARGET_TRAILING:
            trailing_count++;
            break;
          case KLS_FRAGMENTED_UPDATE_TARGET_NONE:
          default:
            none_count++;
            break;
        }
      }
    }
  }
  job->none_count = none_count;
  job->external_count = external_count;
  job->dense_count = dense_count;
  job->pivot_count = pivot_count;
  job->trailing_count = trailing_count;
  return NULL;
}

static void kls_try_build_row_refactor_dense_producer_target_maps(
  UF_long n,
  const UF_long *l_ptr,
  const UF_long *l_cols,
  const UF_long *l_internal_ptr,
  const UF_long *u_cols,
  const UF_long *row_group,
  UF_long group_count,
  const UF_long *group_ptr,
  const unsigned char *group_kind,
  const UF_long *group_trailing_begin,
  const UF_long *group_trailing_len,
  const UF_long *run_ptr,
  const UF_long *run_l_begin,
  const UF_long *run_group,
  const UF_long *run_len,
  UF_long run_count,
  UF_long **target_ptr_out,
  unsigned char **target_kind_out,
  UF_long **target_pos_out,
  UF_long *target_count_out,
  UF_long *target_none_count_out,
  UF_long *target_external_count_out,
  UF_long *target_dense_count_out,
  UF_long *target_pivot_count_out,
  UF_long *target_trailing_count_out,
  int thread_count) {
  KLS_SET_OPTIONAL_OUTPUT(target_ptr_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(target_kind_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(target_pos_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(target_count_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(target_none_count_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(target_external_count_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(target_dense_count_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(target_pivot_count_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(target_trailing_count_out, 0);
  if (n == 0u || run_count == 0u || l_ptr == NULL || l_cols == NULL ||
      l_internal_ptr == NULL || u_cols == NULL || row_group == NULL ||
      group_count == 0u || group_ptr == NULL || group_kind == NULL ||
      group_trailing_begin == NULL || group_trailing_len == NULL ||
      run_ptr == NULL || run_l_begin == NULL || run_group == NULL ||
      run_len == NULL || target_ptr_out == NULL ||
      target_kind_out == NULL || target_pos_out == NULL ||
      run_count > (UF_long)((SIZE_MAX / sizeof(UF_long)) - 1u)) {
    return;
  }

  UF_long *target_ptr =
    (UF_long *)calloc((size_t)run_count + 1u, sizeof(*target_ptr));
  if (target_ptr == NULL) {
    return;
  }

  UF_long target_count = 0;
  for (UF_long run = 0; run < run_count; ++run) {
    const UF_long dep_group = run_group[run];
    if (dep_group >= group_count) {
      free(target_ptr);
      return;
    }
    const UF_long trailing_len = group_trailing_len[dep_group];
    if (target_count > UF_long_max - trailing_len) {
      free(target_ptr);
      return;
    }
    target_count += trailing_len;
    target_ptr[run + 1u] = target_count;
  }
  if (target_count <= 0) {
    free(target_ptr);
    return;
  }
  if (target_count > (UF_long)PTRDIFF_MAX) {
    free(target_ptr);
    return;
  }
  const size_t target_alloc_count = (size_t)target_count;
  if (target_alloc_count == 0u ||
      (UF_long)target_alloc_count != target_count ||
      target_alloc_count > SIZE_MAX / sizeof(UF_long)) {
    free(target_ptr);
    return;
  }

  unsigned char *target_kind =
    (unsigned char *)malloc(target_alloc_count * sizeof(*target_kind));
  UF_long *target_pos =
    (UF_long *)malloc(target_alloc_count * sizeof(*target_pos));
  if (target_kind == NULL || target_pos == NULL) {
    free(target_ptr);
    free(target_kind);
    free(target_pos);
    return;
  }

  kls_dense_producer_target_job common = {
    n, l_ptr, l_cols, l_internal_ptr, u_cols, row_group, group_count,
    group_ptr, group_kind, group_trailing_begin, group_trailing_len, run_ptr,
    run_l_begin, run_group, run_len, target_ptr, target_kind, target_pos,
    0u, n, 0u, 0u, 0u, 0u, 0u, 0
  };
  int workers =
    getenv("KLS_PARALLEL_DENSE_PRODUCER_TARGET_MAPS") != NULL
      ? thread_count : 1;
  if (workers < 1) {
    workers = 1;
  }
  if (workers > 32) {
    workers = 32;
  }
  if ((UF_long)workers > n) {
    workers = (int)n;
  }
  if (target_count < 262144u) {
    workers = 1;
  }

  UF_long none_count = 0;
  UF_long external_count = 0;
  UF_long dense_count = 0;
  UF_long pivot_count = 0;
  UF_long trailing_count = 0;
  int ok = 1;
  kls_dense_producer_target_job *jobs = workers > 1
    ? (kls_dense_producer_target_job *)calloc((size_t)workers,
                                               sizeof(*jobs)) : NULL;
  pthread_t *threads = workers > 1
    ? (pthread_t *)malloc((size_t)workers * sizeof(*threads)) : NULL;
  unsigned char *launched = workers > 1
    ? (unsigned char *)calloc((size_t)workers, sizeof(*launched)) : NULL;
  UF_long *row_bounds = workers > 1
    ? (UF_long *)malloc(((size_t)workers + 1u) * sizeof(*row_bounds)) : NULL;
  if (workers > 1 && jobs != NULL && threads != NULL && launched != NULL &&
      row_bounds != NULL) {
    row_bounds[0] = 0u;
    row_bounds[workers] = n;
    UF_long row = 0u;
    const UF_long target_per_worker = target_count / (UF_long)workers;
    const UF_long target_remainder = target_count % (UF_long)workers;
    for (int tid = 1; tid < workers; ++tid) {
      const UF_long tid_u = (UF_long)tid;
      const UF_long goal = target_per_worker * tid_u +
        (target_remainder * tid_u) / (UF_long)workers;
      while (row < n && target_ptr[run_ptr[row]] < goal) {
        row++;
      }
      row_bounds[tid] = row;
    }
    for (int tid = 0; tid < workers; ++tid) {
      jobs[tid] = common;
      jobs[tid].row_begin = row_bounds[tid];
      jobs[tid].row_end = row_bounds[tid + 1];
    }
    for (int tid = 1; tid < workers; ++tid) {
      if (pthread_create(&threads[tid], NULL,
                         kls_dense_producer_target_main, &jobs[tid]) == 0) {
        launched[tid] = 1u;
      } else {
        (void)kls_dense_producer_target_main(&jobs[tid]);
      }
    }
    (void)kls_dense_producer_target_main(&jobs[0]);
    for (int tid = 1; tid < workers; ++tid) {
      if (launched[tid]) {
        (void)pthread_join(threads[tid], NULL);
      }
    }
    for (int tid = 0; tid < workers; ++tid) {
      ok = ok && jobs[tid].ok;
      none_count += jobs[tid].none_count;
      external_count += jobs[tid].external_count;
      dense_count += jobs[tid].dense_count;
      pivot_count += jobs[tid].pivot_count;
      trailing_count += jobs[tid].trailing_count;
    }
  } else {
    common.row_begin = 0u;
    common.row_end = n;
    (void)kls_dense_producer_target_main(&common);
    ok = common.ok;
    none_count = common.none_count;
    external_count = common.external_count;
    dense_count = common.dense_count;
    pivot_count = common.pivot_count;
    trailing_count = common.trailing_count;
  }
  free(jobs);
  free(threads);
  free(launched);
  free(row_bounds);
  if (!ok) {
    free(target_ptr);
    free(target_kind);
    free(target_pos);
    return;
  }

  *target_ptr_out = target_ptr;
  *target_kind_out = target_kind;
  *target_pos_out = target_pos;
  if (target_count_out != NULL) {
    *target_count_out = target_count;
  }
  if (target_none_count_out != NULL) {
    *target_none_count_out = none_count;
  }
  if (target_external_count_out != NULL) {
    *target_external_count_out = external_count;
  }
  if (target_dense_count_out != NULL) {
    *target_dense_count_out = dense_count;
  }
  if (target_pivot_count_out != NULL) {
    *target_pivot_count_out = pivot_count;
  }
  if (target_trailing_count_out != NULL) {
    *target_trailing_count_out = trailing_count;
  }
}

static void kls_try_build_row_refactor_segment_input_targets(
  UF_long n,
  const UF_long *input_ptr,
  const UF_long *input_cols,
  const unsigned char *input_needs_cleanup,
  const UF_long *l_ptr,
  const UF_long *l_cols,
  const UF_long *u_ptr,
  const UF_long *u_cols,
  const UF_long *row_group,
  UF_long group_count,
  const UF_long *group_ptr,
  const unsigned char *group_kind,
  const UF_long *group_trailing_begin,
  const UF_long *group_trailing_len,
  unsigned char **row_ready_out,
  unsigned char **target_kind_out,
  UF_long **target_pos_out,
  UF_long **cleanup_ptr_out,
  UF_long **cleanup_cols_out,
  UF_long *target_rows_out,
  UF_long *target_entries_out,
  UF_long *cleanup_rows_out,
  UF_long *cleanup_entries_out) {
  KLS_SET_OPTIONAL_OUTPUT(row_ready_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(target_kind_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(target_pos_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(cleanup_ptr_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(cleanup_cols_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(target_rows_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(target_entries_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(cleanup_rows_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(cleanup_entries_out, 0);
  if (n <= 0 || input_ptr == NULL || input_cols == NULL ||
      l_ptr == NULL || l_cols == NULL || u_ptr == NULL || u_cols == NULL ||
      row_group == NULL || group_count <= 0 || group_ptr == NULL ||
      group_kind == NULL || group_trailing_begin == NULL ||
      group_trailing_len == NULL || row_ready_out == NULL ||
      target_kind_out == NULL || target_pos_out == NULL ||
      cleanup_ptr_out == NULL || cleanup_cols_out == NULL ||
      n > (UF_long)PTRDIFF_MAX) {
    return;
  }
  const size_t n_size = (size_t)n;
  if ((UF_long)n_size != n) {
    return;
  }

  const UF_long input_nnz = input_ptr[n];
  if (input_nnz < 0 || input_nnz > (UF_long)PTRDIFF_MAX) {
    return;
  }
  const size_t input_nnz_size = (size_t)input_nnz;
  if ((UF_long)input_nnz_size != input_nnz ||
      input_nnz_size > SIZE_MAX / sizeof(UF_long)) {
    return;
  }

  unsigned char *row_ready =
    (unsigned char *)calloc(n_size, sizeof(*row_ready));
  unsigned char *target_kind = input_nnz > 0u
    ? (unsigned char *)calloc(input_nnz_size, sizeof(*target_kind)) : NULL;
  UF_long *target_pos = input_nnz > 0u
    ? (UF_long *)malloc(input_nnz_size * sizeof(*target_pos)) : NULL;
  UF_long *cleanup_ptr =
    (UF_long *)calloc(n_size + 1u, sizeof(*cleanup_ptr));
  if (row_ready == NULL ||
      (input_nnz > 0u && (target_kind == NULL || target_pos == NULL)) ||
      cleanup_ptr == NULL) {
    free(row_ready);
    free(target_kind);
    free(target_pos);
    free(cleanup_ptr);
    return;
  }

  UF_long target_rows = 0;
  UF_long target_entries = 0;
  UF_long cleanup_rows = 0;
  UF_long cleanup_entries = 0;
  for (UF_long group = 0; group < group_count; ++group) {
    const unsigned char kind = group_kind[group];
    if (kind != KLS_ROW_REFACTOR_GROUP_SINGLE &&
        kind != KLS_ROW_REFACTOR_GROUP_GENERIC &&
        kind != KLS_ROW_REFACTOR_GROUP_DENSE &&
        kind != KLS_ROW_REFACTOR_GROUP_BATCH) {
      continue;
    }
    const UF_long row_begin = group_ptr[group];
    const UF_long row_end = group_ptr[group + 1u];
    if (row_begin < 0 || row_end < 0 || row_begin >= row_end ||
        row_end > n) {
      continue;
    }
    const int single_group = kind == KLS_ROW_REFACTOR_GROUP_SINGLE;
    const int batch_group = kind == KLS_ROW_REFACTOR_GROUP_BATCH;
    if (single_group && row_end != row_begin + 1u) {
      continue;
    }
    const UF_long trailing_begin_for_group = group_trailing_begin[group];
    const UF_long trailing_len_for_group = group_trailing_len[group];
    const UF_long *trailing_cols_for_group = trailing_len_for_group > 0u
      ? u_cols + trailing_begin_for_group : NULL;
    if (!single_group && trailing_len_for_group > 0u &&
        trailing_cols_for_group == NULL) {
      continue;
    }

    for (UF_long row = row_begin; row < row_end; ++row) {
      const UF_long u_begin = u_ptr[row];
      const UF_long u_end = u_ptr[row + 1u];
      if (u_begin > u_end) {
        continue;
      }
      const UF_long dense_len =
        (single_group || batch_group) ? 0u : row_end - row - 1u;
      const UF_long trailing_len =
        (single_group || batch_group) ? u_end - u_begin
                                      : trailing_len_for_group;
      const UF_long *trailing_cols = (single_group || batch_group)
        ? (trailing_len > 0u ? u_cols + u_begin : NULL)
        : trailing_cols_for_group;
      if (trailing_len > 0u && trailing_cols == NULL) {
        continue;
      }
      if (u_begin > u_end || u_end - u_begin != dense_len + trailing_len) {
        continue;
      }
      int row_supported = 1;
      for (UF_long offset = 0; row_supported && offset < dense_len; ++offset) {
        if (u_cols[u_begin + offset] != row + 1u + offset) {
          row_supported = 0;
        }
      }
      for (UF_long offset = 0; row_supported && offset < trailing_len;
           ++offset) {
        if (u_cols[u_begin + dense_len + offset] != trailing_cols[offset]) {
          row_supported = 0;
        }
      }
      const UF_long input_begin = input_ptr[row];
      const UF_long input_end = input_ptr[row + 1u];
      const UF_long l_begin = l_ptr[row];
      const UF_long l_end = l_ptr[row + 1u];
      for (UF_long p = input_begin; row_supported && p < input_end; ++p) {
        const UF_long col = input_cols[p];
        unsigned char target = KLS_ROW_SEGMENT_INPUT_TARGET_NONE;
        UF_long pos = 0;
        if (col >= n) {
          row_supported = 0;
          break;
        }
        if (batch_group) {
          if (col < row) {
            UF_long l_pos = KLS_KLU_EMPTY;
            if (!kls_find_uflong_sorted(l_cols + l_begin, l_end - l_begin,
                                        col, &l_pos) ||
                l_pos >= l_end - l_begin) {
              row_supported = 0;
              break;
            }
            target = KLS_ROW_SEGMENT_INPUT_TARGET_L;
            pos = l_begin + l_pos;
          } else if (col == row) {
            target = KLS_ROW_SEGMENT_INPUT_TARGET_PIVOT;
          } else {
            UF_long u_pos = KLS_KLU_EMPTY;
            if (trailing_len == 0u ||
                !kls_find_uflong_sorted(trailing_cols, trailing_len, col,
                                        &u_pos) ||
                u_pos >= trailing_len) {
              row_supported = 0;
              break;
            }
            target = KLS_ROW_SEGMENT_INPUT_TARGET_U;
            pos = u_pos;
          }
        } else if (single_group && col < row) {
          target = KLS_ROW_SEGMENT_INPUT_TARGET_EXTERNAL;
          pos = col;
        } else if (col < row_begin) {
          target = KLS_ROW_SEGMENT_INPUT_TARGET_EXTERNAL;
          pos = col;
        } else if (col < row) {
          UF_long l_pos = KLS_KLU_EMPTY;
          if (!kls_find_uflong_sorted(l_cols + l_begin, l_end - l_begin,
                                      col, &l_pos) ||
              l_pos >= l_end - l_begin) {
            row_supported = 0;
            break;
          }
          target = KLS_ROW_SEGMENT_INPUT_TARGET_L;
          pos = l_begin + l_pos;
        } else if (col == row) {
          target = KLS_ROW_SEGMENT_INPUT_TARGET_PIVOT;
        } else if (col < row_end) {
          target = KLS_ROW_SEGMENT_INPUT_TARGET_U;
          pos = col - row - 1u;
        } else {
          UF_long trailing_pos = KLS_KLU_EMPTY;
          if (trailing_len == 0u ||
              !kls_find_uflong_sorted(trailing_cols, trailing_len, col,
                                      &trailing_pos) ||
              trailing_pos >= trailing_len) {
            row_supported = 0;
            break;
          }
          target = KLS_ROW_SEGMENT_INPUT_TARGET_U;
          pos = dense_len + trailing_pos;
        }
        target_kind[p] = target;
        target_pos[p] = pos;
      }
      if (!row_supported) {
        continue;
      }
      row_ready[row] = 1u;
      target_rows++;
      target_entries += input_end - input_begin;
      if (input_needs_cleanup != NULL && input_needs_cleanup[row]) {
        UF_long row_cleanup_entries = 0;
        for (UF_long p = input_begin; p < input_end; ++p) {
          if (target_kind[p] == KLS_ROW_SEGMENT_INPUT_TARGET_EXTERNAL) {
            row_cleanup_entries++;
          }
        }
        if (row_cleanup_entries > 0u) {
          cleanup_ptr[row + 1u] = row_cleanup_entries;
          cleanup_rows++;
          cleanup_entries += row_cleanup_entries;
        }
      }
    }
  }

  if (target_rows == 0u) {
    free(row_ready);
    free(target_kind);
    free(target_pos);
    free(cleanup_ptr);
    return;
  }
  UF_long *cleanup_cols = NULL;
  if (cleanup_entries > 0u) {
    if (cleanup_entries > (UF_long)(SIZE_MAX / sizeof(*cleanup_cols))) {
      free(row_ready);
      free(target_kind);
      free(target_pos);
      free(cleanup_ptr);
      return;
    }
    for (UF_long row = 0; row < n; ++row) {
      cleanup_ptr[row + 1u] += cleanup_ptr[row];
    }
    cleanup_cols =
      (UF_long *)malloc((size_t)cleanup_entries * sizeof(*cleanup_cols));
    UF_long *cleanup_write =
      (UF_long *)malloc(n_size * sizeof(*cleanup_write));
    if (cleanup_cols == NULL || cleanup_write == NULL) {
      free(row_ready);
      free(target_kind);
      free(target_pos);
      free(cleanup_ptr);
      free(cleanup_cols);
      free(cleanup_write);
      return;
    }
    memcpy(cleanup_write, cleanup_ptr, n_size * sizeof(*cleanup_write));
    for (UF_long row = 0; row < n; ++row) {
      if (!row_ready[row] || cleanup_ptr[row] == cleanup_ptr[row + 1u]) {
        continue;
      }
      for (UF_long p = input_ptr[row]; p < input_ptr[row + 1u]; ++p) {
        if (target_kind[p] != KLS_ROW_SEGMENT_INPUT_TARGET_EXTERNAL) {
          continue;
        }
        UF_long write = cleanup_write[row]++;
        if (write >= cleanup_ptr[row + 1u]) {
          free(row_ready);
          free(target_kind);
          free(target_pos);
          free(cleanup_ptr);
          free(cleanup_cols);
          free(cleanup_write);
          return;
        }
        cleanup_cols[write] = target_pos[p];
      }
      if (cleanup_write[row] != cleanup_ptr[row + 1u]) {
        free(row_ready);
        free(target_kind);
        free(target_pos);
        free(cleanup_ptr);
        free(cleanup_cols);
        free(cleanup_write);
        return;
      }
    }
    free(cleanup_write);
  } else {
    free(cleanup_ptr);
    cleanup_ptr = NULL;
  }
  *row_ready_out = row_ready;
  *target_kind_out = target_kind;
  *target_pos_out = target_pos;
  *cleanup_ptr_out = cleanup_ptr;
  *cleanup_cols_out = cleanup_cols;
  if (target_rows_out != NULL) {
    *target_rows_out = target_rows;
  }
  if (target_entries_out != NULL) {
    *target_entries_out = target_entries;
  }
  if (cleanup_rows_out != NULL) {
    *cleanup_rows_out = cleanup_rows;
  }
  if (cleanup_entries_out != NULL) {
    *cleanup_entries_out = cleanup_entries;
  }
}

typedef struct {
  UF_long n;
  const UF_long *input_ptr;
  const UF_long *input_cols;
  const unsigned char *input_needs_cleanup;
  const UF_long *l_ptr;
  const UF_long *l_cols;
  const UF_long *u_ptr;
  const UF_long *u_cols;
  const UF_long *row_group;
  UF_long group_count;
  const UF_long *group_ptr;
  const unsigned char *group_kind;
  const UF_long *group_trailing_begin;
  const UF_long *group_trailing_len;
  unsigned char *row_ready;
  unsigned char *target_kind;
  UF_long *target_pos;
  UF_long *cleanup_ptr;
  UF_long *cleanup_cols;
  UF_long target_rows;
  UF_long target_entries;
  UF_long cleanup_rows;
  UF_long cleanup_entries;
  double seconds;
} kls_segment_input_target_prep_job;

static void *kls_segment_input_target_prep_main(void *arg) {
  kls_segment_input_target_prep_job *job =
    (kls_segment_input_target_prep_job *)arg;
  if (job == NULL) {
    return NULL;
  }
  const double begin = kls_now_seconds();
  kls_try_build_row_refactor_segment_input_targets(
    job->n, job->input_ptr, job->input_cols, job->input_needs_cleanup,
    job->l_ptr, job->l_cols, job->u_ptr, job->u_cols, job->row_group,
    job->group_count, job->group_ptr, job->group_kind,
    job->group_trailing_begin, job->group_trailing_len, &job->row_ready,
    &job->target_kind, &job->target_pos, &job->cleanup_ptr,
    &job->cleanup_cols, &job->target_rows, &job->target_entries,
    &job->cleanup_rows, &job->cleanup_entries);
  job->seconds = kls_now_seconds() - begin;
  return NULL;
}

typedef struct {
  UF_long group_count;
  const UF_long *group_ptr;
  const UF_long *l_ptr;
  const UF_long *l_cols;
  const UF_long *u_ptr;
  const UF_long *u_cols;
  UF_long *l_internal_ptr;
  UF_long *group_trailing_len;
  UF_long *group_trailing_begin;
  unsigned char *group_dense;
  unsigned char *group_kind;
  UF_long segment_count;
  UF_long segment_rows;
  UF_long segment_max_width;
  double segment_dense_entries;
  double segment_trailing_entries;
  UF_long dense_segment_count;
  UF_long dense_segment_rows;
  UF_long dense_segment_max_width;
  double dense_segment_dense_entries;
  double dense_segment_trailing_entries;
  UF_long compact_dense_panel_eligible_count;
  UF_long compact_dense_panel_eligible_rows;
  double compact_dense_panel_update_work;
  double compact_dense_panel_entries;
  double seconds;
} kls_dense_group_classification_job;

static void *kls_dense_group_classification_main(void *arg) {
  kls_dense_group_classification_job *job =
    (kls_dense_group_classification_job *)arg;
  if (job == NULL) {
    return NULL;
  }
  const double begin = kls_now_seconds();
  for (UF_long group = 0; group < job->group_count; ++group) {
    const UF_long row_begin = job->group_ptr[group];
    const UF_long row_end = job->group_ptr[group + 1u];
    const UF_long width = row_end - row_begin;
    if (width <= 1u) {
      continue;
    }
    if (!kls_row_refactor_group_has_internal_deps(
          job->l_ptr, job->l_cols, row_begin, row_end)) {
      job->group_kind[group] = KLS_ROW_REFACTOR_GROUP_BATCH;
      continue;
    }
    job->group_kind[group] = KLS_ROW_REFACTOR_GROUP_GENERIC;
    const UF_long trailing_begin = job->u_ptr[row_end - 1u];
    const UF_long trailing_len = job->u_ptr[row_end] - trailing_begin;
    const double dense_entries =
      (double)width * (double)(width - 1u) * 0.5;
    const double trailing_entries = (double)width * (double)trailing_len;
    job->group_trailing_len[group] = trailing_len;
    job->segment_count++;
    job->segment_rows += width;
    if (width > job->segment_max_width) {
      job->segment_max_width = width;
    }
    job->segment_dense_entries += dense_entries;
    job->segment_trailing_entries += trailing_entries;
    job->group_trailing_begin[group] = trailing_begin;

    if (dense_entries + trailing_entries < KLS_ROW_REFACTOR_DENSE_MIN_WORK) {
      continue;
    }
    int dense = 1;
    const UF_long *trailing_cols = job->u_cols + trailing_begin;
    for (UF_long row = row_begin; dense && row < row_end; ++row) {
      UF_long p = job->l_ptr[row];
      const UF_long l_end = job->l_ptr[row + 1u];
      while (p < l_end && job->l_cols[p] < row_begin) {
        p++;
      }
      job->l_internal_ptr[row] = p;
      for (UF_long dep = row_begin; dep < row; ++dep, ++p) {
        if (p >= l_end || job->l_cols[p] != dep) {
          dense = 0;
          break;
        }
      }
      if (dense && p != l_end) {
        dense = 0;
      }

      const UF_long dense_len = row_end - row - 1u;
      const UF_long u_begin = job->u_ptr[row];
      const UF_long u_end = job->u_ptr[row + 1u];
      if (dense && u_end - u_begin != dense_len + trailing_len) {
        dense = 0;
      }
      for (UF_long offset = 0; dense && offset < dense_len; ++offset) {
        if (job->u_cols[u_begin + offset] != row + 1u + offset) {
          dense = 0;
        }
      }
      for (UF_long offset = 0; dense && offset < trailing_len; ++offset) {
        if (job->u_cols[u_begin + dense_len + offset] !=
            trailing_cols[offset]) {
          dense = 0;
        }
      }
    }
    if (!dense) {
      continue;
    }
    job->group_dense[group] = 1u;
    job->group_kind[group] = KLS_ROW_REFACTOR_GROUP_DENSE;
    job->dense_segment_count++;
    job->dense_segment_rows += width;
    if (width > job->dense_segment_max_width) {
      job->dense_segment_max_width = width;
    }
    job->dense_segment_dense_entries += dense_entries;
    job->dense_segment_trailing_entries += trailing_entries;
    if (kls_row_refactor_prefers_compact_dense_panel(width, trailing_len)) {
      job->compact_dense_panel_eligible_count++;
      job->compact_dense_panel_eligible_rows += width;
      job->compact_dense_panel_update_work +=
        kls_row_refactor_dense_group_update_work(width, trailing_len);
      job->compact_dense_panel_entries +=
        kls_row_refactor_dense_group_panel_entries(width, trailing_len);
    }
  }
  job->seconds = kls_now_seconds() - begin;
  return NULL;
}

typedef struct {
  kls_solver *solver;
  int thread_count;
  int prewarm_solve_separator;
  int cached;
  int ordered_private;
  int solve_separator_ready;
  double seconds;
} kls_row_schedule_prewarm_job;

static int kls_row_refactor_private_thread_order_valid(
  const kls_solver *solver,
  const UF_long *private_groups,
  const UF_long *private_thread_ptr,
  const unsigned char *private_group_mask,
  int thread_count);
static int kls_update_row_solve_separator_plan_from_row_schedule(
  kls_solver *solver,
  const UF_long *private_groups,
  const UF_long *private_thread_ptr,
  const unsigned char *private_group_mask,
  int thread_count);
static void *kls_row_schedule_prewarm_main(void *arg);

typedef struct {
  UF_long n;
  const UF_long *l_ptr;
  const UF_long *l_cols;
  const UF_long *row_group;
  const UF_long *group_ptr;
  UF_long group_count;
  UF_long group_begin;
  UF_long group_end;
  UF_long *group_dep_ptr;
  UF_long *group_dep_rows;
  uint32_t *marks;
  uint32_t token_offset;
  int fill;
  int ok;
} kls_group_dep_build_job;

static void *kls_group_dep_build_main(void *arg) {
  kls_group_dep_build_job *job = (kls_group_dep_build_job *)arg;
  if (job == NULL || job->l_ptr == NULL || job->l_cols == NULL ||
      job->row_group == NULL || job->group_ptr == NULL ||
      job->group_dep_ptr == NULL ||
      job->marks == NULL) {
    return NULL;
  }
  job->ok = 1;
  for (UF_long group = job->group_begin; group < job->group_end; ++group) {
    const uint32_t token =
      job->token_offset + (uint32_t)group + 1u;
    UF_long count = 0u;
    UF_long write = job->fill ? job->group_dep_ptr[group] : 0u;
    const UF_long row_begin = job->group_ptr[group];
    const UF_long row_end = job->group_ptr[group + 1u];
    if (row_begin > row_end || row_end > job->n) {
      job->ok = 0;
      return NULL;
    }
    for (UF_long row = row_begin; row < row_end; ++row) {
      for (UF_long p = job->l_ptr[row]; p < job->l_ptr[row + 1u]; ++p) {
        const UF_long dep = job->l_cols[p];
        if (dep >= row || dep >= job->n) {
          job->ok = 0;
          return NULL;
        }
        if (job->row_group[dep] == group || job->marks[dep] == token) {
          continue;
        }
        job->marks[dep] = token;
        if (job->fill) {
          if (write >= job->group_dep_ptr[group + 1u] ||
              job->group_dep_rows == NULL) {
            job->ok = 0;
            return NULL;
          }
          job->group_dep_rows[write++] = dep;
        } else {
          count++;
        }
      }
    }
    if (job->fill) {
      if (write != job->group_dep_ptr[group + 1u]) {
        job->ok = 0;
        return NULL;
      }
    } else {
      job->group_dep_ptr[group + 1u] = count;
    }
  }
  return NULL;
}

static int kls_run_group_dep_build_jobs(kls_group_dep_build_job *jobs,
                                        int job_count) {
  if (jobs == NULL || job_count < 1 || job_count > 8) {
    return 0;
  }
  pthread_t threads[8];
  unsigned char launched[8] = { 0u };
  for (int tid = 1; tid < job_count; ++tid) {
    if (pthread_create(&threads[tid], NULL,
                       kls_group_dep_build_main, &jobs[tid]) == 0) {
      launched[tid] = 1u;
    } else {
      (void)kls_group_dep_build_main(&jobs[tid]);
    }
  }
  (void)kls_group_dep_build_main(&jobs[0]);
  int ok = jobs[0].ok;
  for (int tid = 1; tid < job_count; ++tid) {
    if (launched[tid]) {
      (void)pthread_join(threads[tid], NULL);
    }
    ok = ok && jobs[tid].ok;
  }
  return ok;
}

typedef struct {
  kls_solver *solver;
  UF_long group_begin;
  UF_long group_end;
} kls_group_work_build_job;

static void *kls_group_work_build_main(void *arg) {
  kls_group_work_build_job *job = (kls_group_work_build_job *)arg;
  if (job == NULL || job->solver == NULL ||
      job->solver->row_refactor_group_work == NULL) {
    return NULL;
  }
  for (UF_long group = job->group_begin; group < job->group_end; ++group) {
    job->solver->row_refactor_group_work[group] =
      kls_row_refactor_compute_group_work(job->solver, group);
  }
  return NULL;
}

static void kls_build_group_work_parallel(kls_solver *solver,
                                          UF_long group_count) {
  if (solver == NULL || solver->row_refactor_group_work == NULL) {
    return;
  }
  int job_count = 1;
  if (getenv("KLS_ENABLE_PARALLEL_GROUP_WORK_PREP") != NULL &&
      solver->options.threads >= 4 && group_count >= 10000u) {
    job_count = solver->options.threads - 1;
    if (job_count > 7) {
      job_count = 7;
    }
  }
  kls_group_work_build_job jobs[7];
  pthread_t threads[7];
  unsigned char launched[7] = { 0u };
  const UF_long groups_per_job =
    (group_count + (UF_long)job_count - 1u) / (UF_long)job_count;
  for (int tid = 0; tid < job_count; ++tid) {
    UF_long begin = (UF_long)tid * groups_per_job;
    UF_long end = begin + groups_per_job;
    if (begin > group_count) {
      begin = group_count;
    }
    if (end > group_count) {
      end = group_count;
    }
    jobs[tid] = (kls_group_work_build_job){ solver, begin, end };
  }
  for (int tid = 1; tid < job_count; ++tid) {
    if (pthread_create(&threads[tid], NULL,
                       kls_group_work_build_main, &jobs[tid]) == 0) {
      launched[tid] = 1u;
    } else {
      (void)kls_group_work_build_main(&jobs[tid]);
    }
  }
  (void)kls_group_work_build_main(&jobs[0]);
  for (int tid = 1; tid < job_count; ++tid) {
    if (launched[tid]) {
      (void)pthread_join(threads[tid], NULL);
    }
  }
}

static int kls_finish_row_refactor_pattern_from_arrays(
  kls_solver *solver,
  UF_long *l_ptr,
  UF_long *l_cols,
  double **l_values,
  double *l_row_values,
  UF_long *u_ptr,
  UF_long *u_cols,
  double **u_values,
  double *u_row_values,
  UF_long *input_ptr,
  UF_long *input_cols,
  UF_long *input_pos,
  UF_long *next) {
  const UF_long n = solver != NULL ? solver->n : 0u;
  const UF_long l_nnz = l_ptr != NULL ? l_ptr[n] : 0u;
  const UF_long u_nnz = u_ptr != NULL ? u_ptr[n] : 0u;
  const UF_long input_nnz = input_ptr != NULL ? input_ptr[n] : 0u;
  const int trace_finish =
    getenv("KLS_TRACE_ROW_PATTERN_FINISH") != NULL;
  const double trace_finish_begin =
    trace_finish ? kls_now_seconds() : 0.0;

  double trace_finish_mark = trace_finish_begin;

#define KLS_ROW_FINISH_MARK(name)                                      \
  do {                                                                 \
    if (trace_finish) {                                                 \
      const double trace_finish_now = kls_now_seconds();                \
      fprintf(stderr, "KLS row finish %-18s %.3fms total %.3fms\n",   \
              (name), 1e3 * (trace_finish_now - trace_finish_mark),     \
              1e3 * (trace_finish_now - trace_finish_begin));           \
      trace_finish_mark = trace_finish_now;                             \
    }                                                                  \
  } while (0)
  if (solver == NULL || l_ptr == NULL || u_ptr == NULL ||
      input_ptr == NULL || (n > 0u && next == NULL) ||
      (l_nnz > 0u &&
       (l_cols == NULL || l_values == NULL || l_row_values == NULL)) ||
      (u_nnz > 0u &&
       (u_cols == NULL || u_values == NULL || u_row_values == NULL)) ||
      (input_nnz > 0u && (input_cols == NULL || input_pos == NULL))) {
    free(l_ptr);
    free(l_cols);
    free(l_values);
    free(l_row_values);
    free(u_ptr);
    free(u_cols);
    free(u_values);
    free(u_row_values);
    free(input_ptr);
    free(input_cols);
    free(input_pos);
    free(next);
    return 0;
  }

  UF_long *successor_ptr =
    (UF_long *)calloc((size_t)n + 1u, sizeof(*successor_ptr));
  const int defer_row_successors =
    getenv("KLS_DEFER_ROW_SUCCESSORS") != NULL;
  UF_long *successor_rows = l_nnz > 0u && !defer_row_successors
    ? (UF_long *)malloc((size_t)l_nnz * sizeof(*successor_rows)) : NULL;
  unsigned char *input_needs_cleanup = NULL;
  UF_long input_cleanup_rows = 0;
  UF_long input_cleanup_entries = 0;
  if (successor_ptr == NULL ||
      (l_nnz > 0u && !defer_row_successors && successor_rows == NULL)) {
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
  if (!defer_row_successors) {
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
          free(input_needs_cleanup);
          return 0;
        }
        successor_ptr[dep + 1u]++;
      }
    }
  }
  KLS_ROW_FINISH_MARK("successor counts");

  kls_row_input_cleanup_prep_job input_cleanup_job = {
    .n = n,
    .l_ptr = l_ptr,
    .l_cols = l_cols,
    .u_ptr = u_ptr,
    .u_cols = u_cols,
    .input_ptr = input_ptr,
    .input_cols = input_cols
  };
  const int defer_input_cleanup =
    getenv("KLS_DEFER_ROW_INPUT_CLEANUP") != NULL;
  if (defer_input_cleanup) {
    input_cleanup_job.needs_cleanup =
      (unsigned char *)calloc((size_t)(n > 0u ? n : 1u),
                              sizeof(*input_cleanup_job.needs_cleanup));
    input_cleanup_job.ok = input_cleanup_job.needs_cleanup != NULL;
  }
  pthread_t input_cleanup_thread;
  const int input_cleanup_thread_started =
    !defer_input_cleanup &&
    getenv("KLS_PARALLEL_ROW_INPUT_CLEANUP_PREP") != NULL &&
    pthread_create(&input_cleanup_thread, NULL,
                   kls_row_input_cleanup_prep_main,
                   &input_cleanup_job) == 0;
  kls_row_level_prep_job row_level_job = {
    .n = n,
    .l_ptr = l_ptr,
    .l_cols = l_cols
  };
  const int disable_scalar_row_batches =
    getenv("KLS_DISABLE_SCALAR_ROW_BATCH_GROUPS") != NULL;
  if (disable_scalar_row_batches) {
    row_level_job.ok = 1;
  }
  pthread_t row_level_thread;
  const int row_level_thread_started =
    !disable_scalar_row_batches &&
    getenv("KLS_PARALLEL_ROW_LEVEL_PREP") != NULL &&
    pthread_create(&row_level_thread, NULL, kls_row_level_prep_main,
                   &row_level_job) == 0;
  if (!defer_row_successors) {
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
  }
  free(next);
  if (input_cleanup_thread_started) {
    (void)pthread_join(input_cleanup_thread, NULL);
  } else if (!defer_input_cleanup) {
    (void)kls_row_input_cleanup_prep_main(&input_cleanup_job);
  }
  if (row_level_thread_started) {
    (void)pthread_join(row_level_thread, NULL);
  } else if (!disable_scalar_row_batches) {
    (void)kls_row_level_prep_main(&row_level_job);
  }
  if (trace_finish) {
    fprintf(stderr,
            "KLS row finish prep jobs cleanup %.3fms levels %.3fms\n",
            1e3 * input_cleanup_job.seconds,
            1e3 * row_level_job.seconds);
  }
  input_needs_cleanup = input_cleanup_job.needs_cleanup;
  input_cleanup_rows = input_cleanup_job.cleanup_rows;
  input_cleanup_entries = input_cleanup_job.cleanup_entries;
  if (!input_cleanup_job.ok || !row_level_job.ok) {
    free(l_ptr);
    free(u_ptr);
    free(input_ptr);
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
    free(input_needs_cleanup);
    free(row_level_job.levels);
    free(row_level_job.level_ptr);
    free(row_level_job.level_rows);
    return 0;
  }
  KLS_ROW_FINISH_MARK("successors+cleanup");

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
  solver->row_refactor_input_needs_cleanup = input_needs_cleanup;
  solver->row_refactor_input_cleanup_rows = input_cleanup_rows;
  solver->row_refactor_input_cleanup_entries = input_cleanup_entries;
  solver->row_refactor_successor_ptr = successor_ptr;
  solver->row_refactor_successor_rows = successor_rows;
  solver->row_refactor_pattern_n = n;
  UF_long *levels = row_level_job.levels;
  const UF_long level_count = row_level_job.level_count;
  solver->row_refactor_level_ptr = row_level_job.level_ptr;
  solver->row_refactor_level_rows = row_level_job.level_rows;
  solver->row_refactor_level_count = level_count;
  solver->row_refactor_level_max_width = row_level_job.max_width;
  KLS_ROW_FINISH_MARK("row levels");

  UF_long *row_group = n > 0u
    ? (UF_long *)malloc((size_t)n * sizeof(*row_group)) : NULL;
  UF_long *group_ptr = n > 0u
    ? (UF_long *)malloc(((size_t)n + 1u) * sizeof(*group_ptr)) : NULL;
  if ((n > 0u && row_group == NULL) || (n > 0u && group_ptr == NULL)) {
    free(levels);
    free(row_group);
    free(group_ptr);
    free_row_refactor_pattern(solver);
    return 0;
  }
  UF_long group_count = 0;
  UF_long group_scalar_candidate_count = 0;
  UF_long group_scalar_candidate_rows = 0;
  UF_long group_scalar_short_count = 0;
  UF_long group_scalar_short_rows = 0;
  UF_long group_scalar_stop_level_mismatch_count = 0;
  UF_long group_scalar_stop_internal_dep_count = 0;
  UF_long group_scalar_stop_next_segment_count = 0;
  UF_long group_scalar_stop_max_width_count = 0;
  UF_long group_scalar_stop_matrix_end_count = 0;
  UF_long row = 0;
  while (row < n) {
    group_ptr[group_count] = row;
    UF_long end = row + 1u;
    int segment = 0;
    while (end < n) {
      if (!kls_row_refactor_u_row_extends_next(u_ptr, u_cols, n,
                                               end - 1u, end)) {
        break;
      }
      segment = 1;
      end++;
    }
    /* Coarsen only long scalar runs that live on the same dependency level.
       That trims ready-queue/task overhead without hiding useful small-width
       parallelism or delaying low-level rows behind later dependencies. */
    if (!segment && !disable_scalar_row_batches) {
      int stop_recorded = 0;
      group_scalar_candidate_count++;
      while (end < n && end - row < KLS_ROW_REFACTOR_BATCH_MAX_ROWS) {
        if (levels[end] != levels[row]) {
          group_scalar_stop_level_mismatch_count++;
          stop_recorded = 1;
          break;
        }
        if (end + 1u < n &&
            kls_row_refactor_u_row_extends_next(u_ptr, u_cols, n,
                                                end, end + 1u)) {
          group_scalar_stop_next_segment_count++;
          stop_recorded = 1;
          break;
        }
        if (kls_row_refactor_row_depends_in_range(l_ptr, l_cols, end,
                                                  row, end)) {
          group_scalar_stop_internal_dep_count++;
          stop_recorded = 1;
          break;
        }
        end++;
      }
      if (!stop_recorded) {
        if (end >= n) {
          group_scalar_stop_matrix_end_count++;
        } else {
          group_scalar_stop_max_width_count++;
        }
      }
      const UF_long candidate_rows = end - row;
      group_scalar_candidate_rows += candidate_rows;
      if (candidate_rows < KLS_ROW_REFACTOR_BATCH_MIN_ROWS) {
        group_scalar_short_count++;
        group_scalar_short_rows += candidate_rows;
        end = row + 1u;
      }
    } else if (!segment) {
      group_scalar_candidate_count++;
      group_scalar_candidate_rows++;
      group_scalar_short_count++;
      group_scalar_short_rows++;
    }
    for (UF_long i = row; i < end; ++i) {
      row_group[i] = group_count;
    }
    group_count++;
    row = end;
  }
  group_ptr[group_count] = n;
  free(levels);
  KLS_ROW_FINISH_MARK("group formation");

  if (solver->refactor_level_ptr != NULL &&
      solver->refactor_dependency_work > 0.0 && group_count > 0u &&
      getenv("KLS_DISABLE_ROW_TASK_FLOOR") == NULL) {
    /* Before constructing the remaining group DAG, account for the minimum
       work of the row executor that is already proven by these transposed
       streams.  Each scheduled group necessarily performs a ready/claim,
       completion publication and successor consultation in addition to its
       arithmetic; charge 96 work units for the claim/completion cache line
       plus the successor metadata that every independently scheduled group
       must touch
       irreducible task traffic.  If even this optimistic floor lies outside
       the same five-percent uncertainty band used by the completed row-work
       model, the measured column/EGraph challenger cannot be displaced.
       Aborting here avoids building all downstream schedules and solve
       partitions while preserving row candidates with sufficiently coarse
       work. */
    double row_task_floor = (double)input_nnz + (double)u_nnz;
    for (UF_long row_i = 0u; row_i < n; ++row_i) {
      for (UF_long p = l_ptr[row_i]; p < l_ptr[row_i + 1u]; ++p) {
        const UF_long dep = l_cols[p];
        if (dep < n) {
          row_task_floor += 1.0 +
            (double)(u_ptr[dep + 1u] - u_ptr[dep]);
        }
      }
    }
    row_task_floor += 96.0 * (double)group_count;
    const int reject =
      row_task_floor > 1.05 * solver->refactor_dependency_work;
    if (getenv("KLS_TRACE_ROW_TASK_FLOOR") != NULL) {
      fprintf(stderr,
              "KLS row task floor: work=%.6g groups=%ld floor=%.6g "
              "column=%.6g -> %s\n",
              row_task_floor - 96.0 * (double)group_count,
              (long)group_count, row_task_floor,
              solver->refactor_dependency_work,
              reject ? "reject" : "continue");
    }
    if (reject) {
      free(row_group);
      free(group_ptr);
      free_row_refactor_pattern(solver);
      solver->row_refactor_auto_lower_bound_work = row_task_floor;
      solver->row_refactor_auto_lower_bound_rejected = 1;
      return 0;
    }
  }

  UF_long *l_internal_ptr = n > 0u
    ? (UF_long *)calloc((size_t)n, sizeof(*l_internal_ptr)) : NULL;
  UF_long *group_trailing_len = group_count > 0u
    ? (UF_long *)calloc((size_t)group_count, sizeof(*group_trailing_len))
    : NULL;
  UF_long *group_trailing_begin = group_count > 0u
    ? (UF_long *)calloc((size_t)group_count, sizeof(*group_trailing_begin))
    : NULL;
  unsigned char *group_dense = group_count > 0u
    ? (unsigned char *)calloc((size_t)group_count, sizeof(*group_dense))
    : NULL;
  unsigned char *group_kind = group_count > 0u
    ? (unsigned char *)calloc((size_t)group_count, sizeof(*group_kind))
    : NULL;
  if ((n > 0u && l_internal_ptr == NULL) ||
      (group_count > 0u &&
       (group_trailing_len == NULL || group_trailing_begin == NULL ||
        group_dense == NULL ||
        group_kind == NULL))) {
    free(row_group);
    free(group_ptr);
    free(l_internal_ptr);
    free(group_trailing_len);
    free(group_trailing_begin);
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
    free(group_trailing_begin);
    free(group_dense);
    free(group_kind);
    free(group_dep_ptr);
    free(dep_marks);
    free_row_refactor_pattern(solver);
    return 0;
  }
  kls_dense_group_classification_job dense_classification_job = {
    .group_count = group_count,
    .group_ptr = group_ptr,
    .l_ptr = l_ptr,
    .l_cols = l_cols,
    .u_ptr = u_ptr,
    .u_cols = u_cols,
    .l_internal_ptr = l_internal_ptr,
    .group_trailing_len = group_trailing_len,
    .group_trailing_begin = group_trailing_begin,
    .group_dense = group_dense,
    .group_kind = group_kind
  };
  pthread_t dense_classification_thread;
  const int dense_classification_thread_started =
    getenv("KLS_PARALLEL_DENSE_GROUP_CLASSIFICATION_PREP") != NULL &&
    pthread_create(&dense_classification_thread, NULL,
                   kls_dense_group_classification_main,
                   &dense_classification_job) == 0;
  kls_group_dep_build_job group_dep_jobs[8];
  memset(group_dep_jobs, 0, sizeof(group_dep_jobs));
  int group_dep_job_count = solver->options.threads - 1;
  if (group_dep_job_count > 7) {
    group_dep_job_count = 7;
  }
  int parallel_group_deps =
    getenv("KLS_ENABLE_PARALLEL_GROUP_DEP_PREP") != NULL &&
    group_dep_job_count >= 2 && group_count >= 10000u &&
    group_count <= ((UF_long)UINT32_MAX - 1u) / 2u;
  if (parallel_group_deps) {
    const UF_long groups_per_job =
      (group_count + (UF_long)group_dep_job_count - 1u) /
      (UF_long)group_dep_job_count;
    for (int tid = 0; tid < group_dep_job_count; ++tid) {
      UF_long begin = (UF_long)tid * groups_per_job;
      UF_long end = begin + groups_per_job;
      if (begin > group_count) {
        begin = group_count;
      }
      if (end > group_count) {
        end = group_count;
      }
      group_dep_jobs[tid] = (kls_group_dep_build_job){
        .n = n,
        .l_ptr = l_ptr,
        .l_cols = l_cols,
        .row_group = row_group,
        .group_ptr = group_ptr,
        .group_count = group_count,
        .group_begin = begin,
        .group_end = end,
        .group_dep_ptr = group_dep_ptr,
        .marks = (uint32_t *)calloc(
          (size_t)(n > 0u ? n : 1u), sizeof(uint32_t))
      };
      if (group_dep_jobs[tid].marks == NULL) {
        parallel_group_deps = 0;
      }
    }
  }
  if (!parallel_group_deps) {
    for (int tid = 0; tid < group_dep_job_count; ++tid) {
      free(group_dep_jobs[tid].marks);
      group_dep_jobs[tid].marks = NULL;
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
  } else if (!kls_run_group_dep_build_jobs(group_dep_jobs,
                                            group_dep_job_count)) {
    if (dense_classification_thread_started) {
      (void)pthread_join(dense_classification_thread, NULL);
    }
    for (int tid = 0; tid < group_dep_job_count; ++tid) {
      free(group_dep_jobs[tid].marks);
    }
    free(row_group);
    free(group_ptr);
    free(l_internal_ptr);
    free(group_trailing_len);
    free(group_trailing_begin);
    free(group_dense);
    free(group_kind);
    free(group_dep_ptr);
    free(dep_marks);
    free_row_refactor_pattern(solver);
    return 0;
  }
  for (UF_long g = 0; g < group_count; ++g) {
    group_dep_ptr[g + 1u] += group_dep_ptr[g];
  }
  const UF_long group_dep_count = group_dep_ptr[group_count];
  group_dep_rows = group_dep_count > 0u
    ? (UF_long *)malloc((size_t)group_dep_count * sizeof(*group_dep_rows))
    : NULL;
  UF_long *group_dep_write = !parallel_group_deps && group_count > 0u
    ? (UF_long *)malloc((size_t)group_count * sizeof(*group_dep_write))
    : NULL;
  if ((group_dep_count > 0u && group_dep_rows == NULL) ||
      (!parallel_group_deps && group_count > 0u &&
       group_dep_write == NULL)) {
    if (dense_classification_thread_started) {
      (void)pthread_join(dense_classification_thread, NULL);
    }
    free(row_group);
    free(group_ptr);
    free(l_internal_ptr);
    free(group_trailing_len);
    free(group_trailing_begin);
    free(group_dense);
    free(group_kind);
    free(group_dep_ptr);
    free(group_dep_rows);
    free(group_dep_write);
    free(dep_marks);
    free_row_refactor_pattern(solver);
    return 0;
  }
  if (!parallel_group_deps) {
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
  } else {
    for (int tid = 0; tid < group_dep_job_count; ++tid) {
      group_dep_jobs[tid].group_dep_rows = group_dep_rows;
      group_dep_jobs[tid].token_offset = (uint32_t)group_count;
      group_dep_jobs[tid].fill = 1;
      group_dep_jobs[tid].ok = 0;
    }
    if (!kls_run_group_dep_build_jobs(group_dep_jobs,
                                      group_dep_job_count)) {
      if (dense_classification_thread_started) {
        (void)pthread_join(dense_classification_thread, NULL);
      }
      for (int tid = 0; tid < group_dep_job_count; ++tid) {
        free(group_dep_jobs[tid].marks);
      }
      free(row_group);
      free(group_ptr);
      free(l_internal_ptr);
      free(group_trailing_len);
      free(group_trailing_begin);
      free(group_dense);
      free(group_kind);
      free(group_dep_ptr);
      free(group_dep_rows);
      free(group_dep_write);
      free(dep_marks);
      free_row_refactor_pattern(solver);
      return 0;
    }
    for (int tid = 0; tid < group_dep_job_count; ++tid) {
      free(group_dep_jobs[tid].marks);
      group_dep_jobs[tid].marks = NULL;
    }
  }
  free(group_dep_write);
  free(dep_marks);
  KLS_ROW_FINISH_MARK("group dependencies");
  if (dense_classification_thread_started) {
    (void)pthread_join(dense_classification_thread, NULL);
  } else {
    (void)kls_dense_group_classification_main(&dense_classification_job);
  }
  const UF_long segment_count = dense_classification_job.segment_count;
  const UF_long segment_rows = dense_classification_job.segment_rows;
  const UF_long segment_max_width =
    dense_classification_job.segment_max_width;
  const double segment_dense_entries =
    dense_classification_job.segment_dense_entries;
  const double segment_trailing_entries =
    dense_classification_job.segment_trailing_entries;
  const UF_long dense_segment_count =
    dense_classification_job.dense_segment_count;
  const UF_long dense_segment_rows =
    dense_classification_job.dense_segment_rows;
  const UF_long dense_segment_max_width =
    dense_classification_job.dense_segment_max_width;
  const double dense_segment_dense_entries =
    dense_classification_job.dense_segment_dense_entries;
  const double dense_segment_trailing_entries =
    dense_classification_job.dense_segment_trailing_entries;
  const UF_long compact_dense_panel_eligible_count =
    dense_classification_job.compact_dense_panel_eligible_count;
  const UF_long compact_dense_panel_eligible_rows =
    dense_classification_job.compact_dense_panel_eligible_rows;
  const double compact_dense_panel_update_work =
    dense_classification_job.compact_dense_panel_update_work;
  const double compact_dense_panel_entries =
    dense_classification_job.compact_dense_panel_entries;
  UF_long *group_successor_ptr = NULL;
  UF_long *group_successor_groups = NULL;
  UF_long *group_pred_count = NULL;
  UF_long *group_roots = NULL;
  UF_long group_dependency_edges = 0;
  UF_long group_root_count = 0;
  UF_long group_leaf_count = 0;
  UF_long group_max_fanout = 0;
  KLS_ROW_FINISH_MARK("dense classification");

  kls_group_successor_prep_job group_successor_job = {
    .group_count = group_count,
    .row_group = row_group,
    .group_dep_ptr = group_dep_ptr,
    .group_dep_rows = group_dep_rows
  };
  pthread_t group_successor_thread;
  const int group_successor_thread_started =
    getenv("KLS_PARALLEL_GROUP_SUCCESSOR_PREP") != NULL &&
    pthread_create(&group_successor_thread, NULL,
                   kls_group_successor_prep_main,
                   &group_successor_job) == 0;

  UF_long *group_levels = group_count > 0u
    ? (UF_long *)calloc((size_t)group_count, sizeof(*group_levels)) : NULL;
  if (group_count > 0u && group_levels == NULL) {
    if (group_successor_thread_started) {
      (void)pthread_join(group_successor_thread, NULL);
    }
    kls_free_group_successor_prep_job(&group_successor_job);
    free(row_group);
    free(group_ptr);
    free(l_internal_ptr);
    free(group_trailing_len);
    free(group_trailing_begin);
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
          if (group_successor_thread_started) {
            (void)pthread_join(group_successor_thread, NULL);
          }
          kls_free_group_successor_prep_job(&group_successor_job);
          free(row_group);
          free(group_ptr);
          free(group_levels);
          free(l_internal_ptr);
          free(group_trailing_len);
          free(group_trailing_begin);
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
    if (group_successor_thread_started) {
      (void)pthread_join(group_successor_thread, NULL);
    }
    kls_free_group_successor_prep_job(&group_successor_job);
    free(row_group);
    free(group_ptr);
    free(group_levels);
    free(l_internal_ptr);
    free(group_trailing_len);
    free(group_trailing_begin);
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
    if (group_successor_thread_started) {
      (void)pthread_join(group_successor_thread, NULL);
    }
    kls_free_group_successor_prep_job(&group_successor_job);
    free(row_group);
    free(group_ptr);
    free(group_levels);
    free(l_internal_ptr);
    free(group_trailing_len);
    free(group_trailing_begin);
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
  KLS_ROW_FINISH_MARK("group level plan");
  if (group_successor_thread_started) {
    (void)pthread_join(group_successor_thread, NULL);
  } else {
    (void)kls_group_successor_prep_main(&group_successor_job);
  }
  group_successor_ptr = group_successor_job.successor_ptr;
  group_successor_groups = group_successor_job.successor_groups;
  group_pred_count = group_successor_job.pred_count;
  group_roots = group_successor_job.root_groups;
  group_dependency_edges = group_successor_job.edge_count;
  group_root_count = group_successor_job.root_count;
  group_leaf_count = group_successor_job.leaf_count;
  group_max_fanout = group_successor_job.max_fanout;
  if (!group_successor_job.ok) {
    free(row_group);
    free(group_ptr);
    free(group_levels);
    free(l_internal_ptr);
    free(group_trailing_len);
    free(group_trailing_begin);
    free(group_dense);
    free(group_kind);
    free(group_dep_ptr);
    free(group_dep_rows);
    free(group_level_ptr);
    free(level_groups);
    free(group_successor_ptr);
    free(group_successor_groups);
    free(group_pred_count);
    free(group_roots);
    free_row_refactor_pattern(solver);
    return 0;
  }
  free(group_levels);
  KLS_ROW_FINISH_MARK("group successors");

  /* Row-solve pattern construction only reads L/U.  Start it before the
     group work/sort passes so its destructive initialization completes
     before the schedule prewarm derives the separator solve partition. */
  pthread_t solve_partition_thread;
  kls_row_solve_partition_prep_job solve_partition_job = { solver };
  const int solve_partition_thread_started =
    getenv("KLS_PARALLEL_ROW_SOLVE_PARTITION_PREP") != NULL &&
    pthread_create(&solve_partition_thread, NULL,
                   kls_row_solve_partition_prep_main,
                   &solve_partition_job) == 0;
  int solve_partition_thread_joined = 0;

  /* Install the completed group DAG before building optional execution
     metadata.  This lets the reusable separator schedule run concurrently
     with producer-target construction below. */
  solver->row_refactor_group_ptr = group_ptr;
  solver->row_refactor_group_dep_ptr = group_dep_ptr;
  solver->row_refactor_group_dep_rows = group_dep_rows;
  solver->row_refactor_group_successor_ptr = group_successor_ptr;
  solver->row_refactor_group_successor_groups = group_successor_groups;
  solver->row_refactor_group_pred_count = group_pred_count;
  solver->row_refactor_group_roots = group_roots;
  solver->row_refactor_group_count = group_count;
  solver->row_refactor_group_dependency_edges = group_dependency_edges;
  solver->row_refactor_group_root_count = group_root_count;
  solver->row_refactor_group_leaf_count = group_leaf_count;
  solver->row_refactor_group_max_fanout = group_max_fanout;
  solver->row_refactor_group_level_ptr = group_level_ptr;
  solver->row_refactor_level_groups = level_groups;
  solver->row_refactor_row_group = row_group;
  solver->row_refactor_l_internal_ptr = l_internal_ptr;
  solver->row_refactor_group_trailing_len = group_trailing_len;
  solver->row_refactor_group_trailing_begin = group_trailing_begin;
  solver->row_refactor_group_dense = group_dense;
  solver->row_refactor_group_kind = group_kind;
  /* Group classification already proves the immutable triangular/tail
     geometry.  Retain that proof as a byte-per-group plan so the steady
     owner-compute solve can dispatch without repeating the O(width) shape
     walk on every RHS.  Allocation failure only disables the fast plan; the
     existing validating kernels remain authoritative. */
  /* One byte per group cannot overflow because group_count is already a
     validated matrix dimension. */
  if (group_count > 0u) {
    solver->row_refactor_group_shape_valid =
      (unsigned char *)malloc((size_t)group_count *
                              sizeof(*solver->row_refactor_group_shape_valid));
    if (solver->row_refactor_group_shape_valid != NULL) {
      for (UF_long group = 0u; group < group_count; ++group) {
        solver->row_refactor_group_shape_valid[group] =
          group_kind[group] == KLS_ROW_REFACTOR_GROUP_DENSE ? 1u : 0u;
      }
    }
  }
  solver->row_refactor_level_count = group_level_count;
  solver->row_refactor_level_max_width = group_max_width;
  if (group_count > 0u &&
      group_count <= (UF_long)(SIZE_MAX /
                               sizeof(*solver->row_refactor_group_work))) {
    solver->row_refactor_group_work =
      (double *)malloc((size_t)group_count *
                       sizeof(*solver->row_refactor_group_work));
    if (solver->row_refactor_group_work != NULL) {
      kls_build_group_work_parallel(solver, group_count);
    }
  }
  KLS_ROW_FINISH_MARK("group work");
  if ((getenv("KLS_SKIP_SUCCESSOR_SORT") == NULL)) {
    kls_sort_row_refactor_successors_by_work(solver);
  }
  KLS_ROW_FINISH_MARK("successor sort");
  if (solve_partition_thread_started) {
    (void)pthread_join(solve_partition_thread, NULL);
    solve_partition_thread_joined = 1;
  }

  kls_row_schedule_prewarm_job schedule_prewarm_job = {
    solver, solver->options.threads, solve_partition_thread_joined,
    0, 0, 0, 0.0
  };
  pthread_t schedule_prewarm_thread;
  const int schedule_prewarm_thread_started =
    getenv("KLS_ENABLE_ROW_SCHEDULE_PREWARM") != NULL &&
    solver->options.threads >= 4 &&
    pthread_create(&schedule_prewarm_thread, NULL,
                   kls_row_schedule_prewarm_main,
                   &schedule_prewarm_job) == 0;

  kls_segment_input_target_prep_job segment_input_job = {
    .n = n,
    .input_ptr = input_ptr,
    .input_cols = input_cols,
    .input_needs_cleanup = solver->row_refactor_input_needs_cleanup,
    .l_ptr = l_ptr,
    .l_cols = l_cols,
    .u_ptr = u_ptr,
    .u_cols = u_cols,
    .row_group = row_group,
    .group_count = group_count,
    .group_ptr = group_ptr,
    .group_kind = group_kind,
    .group_trailing_begin = group_trailing_begin,
    .group_trailing_len = group_trailing_len
  };
  pthread_t segment_input_thread;
  const int segment_input_thread_started =
    getenv("KLS_PARALLEL_SEGMENT_INPUT_TARGET_PREP") != NULL &&
    pthread_create(&segment_input_thread, NULL,
                   kls_segment_input_target_prep_main,
                   &segment_input_job) == 0;

  UF_long *group_compact_panel_begin = NULL;
  double *compact_panel_values = NULL;
  atomic_uint *compact_panel_valid = NULL;
  UF_long compact_dense_panel_persistent_groups = 0;
  UF_long compact_dense_panel_persistent_entries = 0;
  kls_try_build_row_refactor_compact_panel_storage(
    group_count, group_ptr, group_kind, group_trailing_len,
    &group_compact_panel_begin, &compact_panel_values,
    &compact_dense_panel_persistent_groups,
    &compact_dense_panel_persistent_entries);
  if (group_count > 0u &&
      group_count <= (UF_long)(SIZE_MAX / sizeof(*compact_panel_valid))) {
    compact_panel_valid =
      (atomic_uint *)malloc((size_t)group_count * sizeof(*compact_panel_valid));
    if (compact_panel_valid != NULL) {
      for (UF_long g = 0; g < group_count; ++g) {
        atomic_init(&compact_panel_valid[g], 0u);
      }
    }
  }
  KLS_ROW_FINISH_MARK("compact storage");

  UF_long *dense_producer_run_ptr = NULL;
  UF_long *dense_producer_run_l_begin = NULL;
  UF_long *dense_producer_run_group = NULL;
  UF_long *dense_producer_run_len = NULL;
  UF_long *dense_producer_target_ptr = NULL;
  unsigned char *dense_producer_target_kind = NULL;
  UF_long *dense_producer_target_pos = NULL;
  UF_long dense_producer_run_count = 0;
  UF_long dense_producer_run_rows = 0;
  UF_long dense_producer_run_dep_rows = 0;
  UF_long dense_producer_run_max_per_row = 0;
  UF_long dense_producer_full_suffix_run_count = 0;
  UF_long dense_producer_full_suffix_rows = 0;
  UF_long dense_producer_multi_run_rows = 0;
  UF_long dense_producer_fragmented_rows = 0;
  UF_long dense_producer_target_count = 0;
  UF_long dense_producer_target_none_count = 0;
  UF_long dense_producer_target_external_count = 0;
  UF_long dense_producer_target_dense_count = 0;
  UF_long dense_producer_target_pivot_count = 0;
  UF_long dense_producer_target_trailing_count = 0;
  if (!segment_input_thread_started) {
    (void)kls_segment_input_target_prep_main(&segment_input_job);
    KLS_ROW_FINISH_MARK("segment targets");
  }
  const int dense_producer_run_threads =
    getenv("KLS_ENABLE_PARALLEL_DENSE_PRODUCER_RUNS") != NULL &&
    solver->options.threads > 2
      ? solver->options.threads - 2 : 1;
  kls_try_build_row_refactor_dense_producer_runs(
    n, l_ptr, l_cols, row_group, group_count, group_ptr, group_kind,
    group_compact_panel_begin, &dense_producer_run_ptr,
    &dense_producer_run_l_begin, &dense_producer_run_group,
    &dense_producer_run_len, &dense_producer_run_count,
    &dense_producer_run_rows, &dense_producer_run_dep_rows,
    &dense_producer_run_max_per_row,
    &dense_producer_full_suffix_run_count,
    &dense_producer_full_suffix_rows,
    &dense_producer_multi_run_rows,
    &dense_producer_fragmented_rows,
    dense_producer_run_threads);
  KLS_ROW_FINISH_MARK("producer runs");
  kls_try_build_row_refactor_dense_producer_target_maps(
    n, l_ptr, l_cols, l_internal_ptr, u_cols, row_group, group_count,
    group_ptr, group_kind, group_trailing_begin, group_trailing_len,
    dense_producer_run_ptr, dense_producer_run_l_begin,
    dense_producer_run_group, dense_producer_run_len,
    dense_producer_run_count, &dense_producer_target_ptr,
    &dense_producer_target_kind, &dense_producer_target_pos,
    &dense_producer_target_count, &dense_producer_target_none_count,
    &dense_producer_target_external_count,
    &dense_producer_target_dense_count, &dense_producer_target_pivot_count,
    &dense_producer_target_trailing_count,
    solver->options.threads - (solve_partition_thread_started ? 1 : 0));
  KLS_ROW_FINISH_MARK("producer targets");

  UF_long group_single_count = 0u;
  UF_long group_batch_count = 0u;
  UF_long group_batch_rows = 0u;
  UF_long group_batch_max_width = 0u;
  UF_long group_batch_width_le_4_count = 0u;
  UF_long group_batch_width_le_8_count = 0u;
  UF_long group_generic_count = 0u;
  UF_long group_generic_rows = 0u;
  UF_long group_generic_max_width = 0u;
  UF_long group_dense_count = 0u;
  UF_long group_dense_rows = 0u;
  UF_long group_dense_max_width = 0u;
  for (UF_long group = 0; group < group_count; ++group) {
    const UF_long group_rows = group_ptr[group + 1u] - group_ptr[group];
    switch ((kls_row_refactor_group_kind)group_kind[group]) {
      case KLS_ROW_REFACTOR_GROUP_BATCH:
        group_batch_count++;
        group_batch_rows += group_rows;
        if (group_rows > group_batch_max_width) {
          group_batch_max_width = group_rows;
        }
        if (group_rows <= 4u) {
          group_batch_width_le_4_count++;
        }
        if (group_rows <= 8u) {
          group_batch_width_le_8_count++;
        }
        break;
      case KLS_ROW_REFACTOR_GROUP_GENERIC:
        group_generic_count++;
        group_generic_rows += group_rows;
        if (group_rows > group_generic_max_width) {
          group_generic_max_width = group_rows;
        }
        break;
      case KLS_ROW_REFACTOR_GROUP_DENSE:
        group_dense_count++;
        group_dense_rows += group_rows;
        if (group_rows > group_dense_max_width) {
          group_dense_max_width = group_rows;
        }
        break;
      case KLS_ROW_REFACTOR_GROUP_SINGLE:
      default:
        group_single_count++;
        break;
      }
  }
  if (segment_input_thread_started) {
    (void)pthread_join(segment_input_thread, NULL);
    KLS_ROW_FINISH_MARK("segment target join");
  }
  if (schedule_prewarm_thread_started) {
    (void)pthread_join(schedule_prewarm_thread, NULL);
  }
  if (trace_finish && schedule_prewarm_thread_started) {
    fprintf(stderr,
            "KLS row schedule prewarm %.3fms cached %d ordered %d solve %d\n",
            1e3 * schedule_prewarm_job.seconds,
            schedule_prewarm_job.cached,
            schedule_prewarm_job.ordered_private,
            schedule_prewarm_job.solve_separator_ready);
  }
  KLS_ROW_FINISH_MARK("schedule prewarm");

  solver->row_refactor_group_ptr = group_ptr;
  solver->row_refactor_group_dep_ptr = group_dep_ptr;
  solver->row_refactor_group_dep_rows = group_dep_rows;
  solver->row_refactor_group_successor_ptr = group_successor_ptr;
  solver->row_refactor_group_successor_groups = group_successor_groups;
  solver->row_refactor_group_pred_count = group_pred_count;
  solver->row_refactor_group_roots = group_roots;
  solver->row_refactor_group_count = group_count;
  solver->row_refactor_group_single_count = group_single_count;
  solver->row_refactor_group_batch_count = group_batch_count;
  solver->row_refactor_group_batch_rows = group_batch_rows;
  solver->row_refactor_group_batch_max_width = group_batch_max_width;
  solver->row_refactor_group_batch_width_le_4_count =
    group_batch_width_le_4_count;
  solver->row_refactor_group_batch_width_le_8_count =
    group_batch_width_le_8_count;
  solver->row_refactor_group_scalar_candidate_count =
    group_scalar_candidate_count;
  solver->row_refactor_group_scalar_candidate_rows =
    group_scalar_candidate_rows;
  solver->row_refactor_group_scalar_short_count =
    group_scalar_short_count;
  solver->row_refactor_group_scalar_short_rows =
    group_scalar_short_rows;
  solver->row_refactor_group_scalar_stop_level_mismatch_count =
    group_scalar_stop_level_mismatch_count;
  solver->row_refactor_group_scalar_stop_internal_dep_count =
    group_scalar_stop_internal_dep_count;
  solver->row_refactor_group_scalar_stop_next_segment_count =
    group_scalar_stop_next_segment_count;
  solver->row_refactor_group_scalar_stop_max_width_count =
    group_scalar_stop_max_width_count;
  solver->row_refactor_group_scalar_stop_matrix_end_count =
    group_scalar_stop_matrix_end_count;
  solver->row_refactor_group_generic_count = group_generic_count;
  solver->row_refactor_group_generic_rows = group_generic_rows;
  solver->row_refactor_group_generic_max_width =
    group_generic_max_width;
  solver->row_refactor_group_dense_count = group_dense_count;
  solver->row_refactor_group_dense_rows = group_dense_rows;
  solver->row_refactor_group_dense_max_width = group_dense_max_width;
  solver->row_refactor_group_dependency_edges = group_dependency_edges;
  solver->row_refactor_group_root_count = group_root_count;
  solver->row_refactor_group_leaf_count = group_leaf_count;
  solver->row_refactor_group_max_fanout = group_max_fanout;
  solver->row_refactor_group_level_ptr = group_level_ptr;
  solver->row_refactor_level_groups = level_groups;
  solver->row_refactor_row_group = row_group;
  solver->row_refactor_l_internal_ptr = l_internal_ptr;
  solver->row_refactor_group_trailing_len = group_trailing_len;
  solver->row_refactor_group_trailing_begin = group_trailing_begin;
  solver->row_refactor_group_dense = group_dense;
  solver->row_refactor_group_kind = group_kind;
  solver->row_refactor_group_compact_panel_begin =
    group_compact_panel_begin;
  solver->row_refactor_compact_panel_values =
    compact_panel_values;
  solver->row_refactor_compact_panel_valid =
    compact_panel_valid;
  solver->row_refactor_dense_producer_run_ptr =
    dense_producer_run_ptr;
  solver->row_refactor_dense_producer_run_l_begin =
    dense_producer_run_l_begin;
  solver->row_refactor_dense_producer_run_group =
    dense_producer_run_group;
  solver->row_refactor_dense_producer_run_len =
    dense_producer_run_len;
  solver->row_refactor_dense_producer_target_ptr =
    dense_producer_target_ptr;
  solver->row_refactor_dense_producer_target_kind =
    dense_producer_target_kind;
  solver->row_refactor_dense_producer_target_pos =
    dense_producer_target_pos;
  solver->row_refactor_segment_input_row_ready =
    segment_input_job.row_ready;
  solver->row_refactor_segment_input_target_kind =
    segment_input_job.target_kind;
  solver->row_refactor_segment_input_target_pos =
    segment_input_job.target_pos;
  solver->row_refactor_segment_input_cleanup_ptr =
    segment_input_job.cleanup_ptr;
  solver->row_refactor_segment_input_cleanup_cols =
    segment_input_job.cleanup_cols;
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
  solver->row_refactor_compact_dense_panel_eligible_count =
    compact_dense_panel_eligible_count;
  solver->row_refactor_compact_dense_panel_eligible_rows =
    compact_dense_panel_eligible_rows;
  solver->row_refactor_compact_dense_panel_update_work =
    compact_dense_panel_update_work;
  solver->row_refactor_compact_dense_panel_entries =
    compact_dense_panel_entries;
  solver->row_refactor_dense_producer_run_count =
    dense_producer_run_count;
  solver->row_refactor_dense_producer_run_rows =
    dense_producer_run_rows;
  solver->row_refactor_dense_producer_run_dep_rows =
    dense_producer_run_dep_rows;
  solver->row_refactor_dense_producer_run_max_per_row =
    dense_producer_run_max_per_row;
  solver->row_refactor_dense_producer_full_suffix_run_count =
    dense_producer_full_suffix_run_count;
  solver->row_refactor_dense_producer_full_suffix_rows =
    dense_producer_full_suffix_rows;
  solver->row_refactor_dense_producer_multi_run_rows =
    dense_producer_multi_run_rows;
  solver->row_refactor_dense_producer_fragmented_rows =
    dense_producer_fragmented_rows;
  solver->row_refactor_dense_producer_target_count =
    dense_producer_target_count;
  solver->row_refactor_dense_producer_target_none_count =
    dense_producer_target_none_count;
  solver->row_refactor_dense_producer_target_external_count =
    dense_producer_target_external_count;
  solver->row_refactor_dense_producer_target_dense_count =
    dense_producer_target_dense_count;
  solver->row_refactor_dense_producer_target_pivot_count =
    dense_producer_target_pivot_count;
  solver->row_refactor_dense_producer_target_trailing_count =
    dense_producer_target_trailing_count;
  solver->row_refactor_segment_input_target_rows =
    segment_input_job.target_rows;
  solver->row_refactor_segment_input_target_entries =
    segment_input_job.target_entries;
  solver->row_refactor_segment_input_cleanup_rows =
    segment_input_job.cleanup_rows;
  solver->row_refactor_segment_input_cleanup_entries =
    segment_input_job.cleanup_entries;
  solver->row_refactor_compact_dense_panel_persistent_groups =
    compact_dense_panel_persistent_groups;
  solver->row_refactor_compact_dense_panel_persistent_entries =
    compact_dense_panel_persistent_entries;
  if (solver->row_refactor_group_work == NULL && group_count > 0u &&
      group_count <= (UF_long)(SIZE_MAX / sizeof(double))) {
    solver->row_refactor_group_work =
      (double *)malloc((size_t)group_count *
                       sizeof(*solver->row_refactor_group_work));
    if (solver->row_refactor_group_work != NULL) {
      for (UF_long group = 0; group < group_count; ++group) {
        solver->row_refactor_group_work[group] =
          kls_row_refactor_compute_group_work(solver, group);
      }
    }
  }
  solver->row_refactor_group_single_work = 0.0;
  solver->row_refactor_group_batch_work = 0.0;
  solver->row_refactor_group_generic_work = 0.0;
  solver->row_refactor_group_dense_work = 0.0;
  for (UF_long group = 0; group < group_count; ++group) {
    const double work = kls_row_refactor_group_work(solver, group);
    switch ((kls_row_refactor_group_kind)group_kind[group]) {
      case KLS_ROW_REFACTOR_GROUP_BATCH:
        solver->row_refactor_group_batch_work += work;
        break;
      case KLS_ROW_REFACTOR_GROUP_GENERIC:
        solver->row_refactor_group_generic_work += work;
        break;
      case KLS_ROW_REFACTOR_GROUP_DENSE:
        solver->row_refactor_group_dense_work += work;
        break;
      case KLS_ROW_REFACTOR_GROUP_SINGLE:
      default:
        solver->row_refactor_group_single_work += work;
        break;
    }
  }
  KLS_ROW_FINISH_MARK("group stats");
  if (solve_partition_thread_started) {
    if (!solve_partition_thread_joined) {
      (void)pthread_join(solve_partition_thread, NULL);
    }
  } else {
    kls_record_row_solve_partition(solver);
  }
  KLS_ROW_FINISH_MARK("solve partition");
#undef KLS_ROW_FINISH_MARK
  return 1;
}

/* Diagnostic census of SubtreeLU-style strict row-supernodes over the
   row engine's U row patterns (KLS_TRACE_ROW_SNODE=1): consecutive rows
   where row i's U pattern equals row i-1's minus its leading entry i
   (diagonal-chain shift). Also measures dependency-segment coverage —
   how many of each row's L-dep updates arrive as complete supernode
   segments, i.e. the dense-consume batching a supernodal kernel gets. */
static void kls_row_snode_census(const kls_solver *solver) {
  const UF_long n = solver->n;
  const UF_long *up = solver->row_refactor_u_ptr;
  const UF_long *uc = solver->row_refactor_u_cols;
  const UF_long *lp = solver->row_refactor_l_ptr;
  const UF_long *lc = solver->row_refactor_l_cols;
  if (n == 0u || up == NULL || uc == NULL || lp == NULL || lc == NULL) {
    return;
  }
  UF_long *sn_of = (UF_long *)malloc((size_t)n * sizeof(*sn_of));
  UF_long *sn_end = (UF_long *)malloc((size_t)n * sizeof(*sn_end));
  if (sn_of == NULL || sn_end == NULL) {
    free(sn_of);
    free(sn_end);
    return;
  }
  UF_long unsorted = 0, diag_first = 0;
  for (UF_long i = 0; i < n; ++i) {
    UF_long b = up[i], e = up[i + 1u];
    if (e > b && uc[b] == i) {
      diag_first++;
    }
    for (UF_long p = b + 1u; p < e; ++p) {
      if (uc[p] <= uc[p - 1u]) {
        unsorted++;
        break;
      }
    }
  }
  UF_long nsn = 0, rows_in = 0, max_h = 0, u_in = 0, u_total = 0;
  UF_long start = 0;
  for (UF_long i = 0; i <= n; ++i) {
    int joins = 0;
    if (i > 0u && i < n) {
      const UF_long pb = up[i - 1u], pe = up[i];
      const UF_long cb = up[i], ce = up[i + 1u];
      const UF_long plen = pe - pb, clen = ce - cb;
      joins = plen == clen + 1u && plen > 0u && uc[pb] == i &&
              (clen == 0u ||
               memcmp(uc + pb + 1u, uc + cb,
                      (size_t)clen * sizeof(*uc)) == 0);
    }
    if (!joins && i > 0u) {
      const UF_long h = i - start;
      if (h >= 2u) {
        nsn++;
        rows_in += h;
        if (h > max_h) {
          max_h = h;
        }
        for (UF_long r = start; r < i; ++r) {
          u_in += up[r + 1u] - up[r];
        }
      }
      for (UF_long r = start; r < i; ++r) {
        sn_of[r] = start;
        sn_end[r] = i - 1u;
      }
      start = i;
    }
  }
  u_total = up[n];
  /* dependency-segment coverage: consume calls now (one per dep row)
     vs with supernode segments (one per maximal in-supernode run) */
  double dep_entries = 0.0, dep_in_runs = 0.0;
  double calls_now = 0.0, calls_sn = 0.0;
  for (UF_long i = 0; i < n; ++i) {
    const UF_long b = lp[i], e = lp[i + 1u];
    UF_long p = b;
    while (p < e) {
      const UF_long k = lc[p];
      UF_long q = p + 1u;
      if (k < n && sn_end[k] > k) {
        while (q < e && lc[q] == lc[q - 1u] + 1u && lc[q] <= sn_end[k] &&
               sn_of[lc[q]] == sn_of[k]) {
          q++;
        }
      }
      const UF_long run = q - p;
      double seg_entries = 0.0;
      for (UF_long t = p; t < q; ++t) {
        if (lc[t] < n) {
          seg_entries += (double)(up[lc[t] + 1u] - up[lc[t]]);
        }
      }
      dep_entries += seg_entries;
      calls_now += (double)run;
      calls_sn += 1.0;
      if (run >= 2u) {
        dep_in_runs += seg_entries;
      }
      p = q;
    }
  }
  fprintf(stderr,
          "KLS row-snode census: n=%ld usn=%ld rows_in=%ld (%.1f%%)"
          " mean_h=%.2f max_h=%ld u_in=%.1f%% of %ld |"
          " deps: entries=%.3g in_runs=%.1f%% calls %.3g->%.3g (%.2fx)"
          " | unsorted=%ld diag_first=%ld\n",
          (long)n, (long)nsn, (long)rows_in,
          n > 0u ? 100.0 * (double)rows_in / (double)n : 0.0,
          nsn > 0u ? (double)rows_in / (double)nsn : 0.0, (long)max_h,
          u_total > 0u ? 100.0 * (double)u_in / (double)u_total : 0.0,
          (long)u_total, dep_entries,
          dep_entries > 0.0 ? 100.0 * dep_in_runs / dep_entries : 0.0,
          calls_now, calls_sn,
          calls_sn > 0.0 ? calls_now / calls_sn : 0.0,
          (long)unsorted, (long)diag_first);
  free(sn_of);
  free(sn_end);
}

/* Strict U-row supernode map for the light-run consume: rows i-1 and i
   share a supernode when row i's U pattern equals row i-1's minus its
   leading entry i (diagonal-chain shift, the SubtreeLU/NICSLU strict
   criterion). Rebuilt with the pattern; single-threaded callers only. */
static void kls_row_refactor_ensure_sn_end(kls_solver *solver) {
  if (solver == NULL) {
    return;
  }
  {
    const char *trace = getenv("KLS_TRACE_ROW_SNODE");
    if (trace != NULL && atoi(trace) >= 2) {
      fprintf(stderr,
              "KLS light-snode: runs=%llu entries=%llu scalar=%llu"
              " decl(norun=%llu short=%llu wait=%llu ok=%llu)\n",
              (unsigned long long)atomic_load_explicit(
                &solver->row_light_snode_runs, memory_order_relaxed),
              (unsigned long long)atomic_load_explicit(
                &solver->row_light_snode_entries, memory_order_relaxed),
              (unsigned long long)atomic_load_explicit(
                &solver->row_scalar_dep_entries, memory_order_relaxed),
              (unsigned long long)atomic_load_explicit(
                &solver->lsn_decl_norun, memory_order_relaxed),
              (unsigned long long)atomic_load_explicit(
                &solver->lsn_decl_short, memory_order_relaxed),
              (unsigned long long)atomic_load_explicit(
                &solver->lsn_decl_wait, memory_order_relaxed),
              (unsigned long long)atomic_load_explicit(
                &solver->lsn_consumed, memory_order_relaxed));
    }
  }
  if (solver->row_refactor_sn_end != NULL ||
      solver->row_refactor_u_ptr == NULL ||
      solver->row_refactor_u_cols == NULL || solver->n == 0u ||
      getenv("KLS_DISABLE_ROW_LIGHT_SNODE") != NULL) {
    return;
  }
  const UF_long n = solver->n;
  const UF_long *up = solver->row_refactor_u_ptr;
  const UF_long *uc = solver->row_refactor_u_cols;
  UF_long *sn_end = (UF_long *)malloc((size_t)n * sizeof(*sn_end));
  if (sn_end == NULL) {
    return;
  }
  UF_long start = 0;
  for (UF_long i = 1; i <= n; ++i) {
    int joins = 0;
    if (i < n) {
      const UF_long pb = up[i - 1u], pe = up[i];
      const UF_long cb = up[i], ce = up[i + 1u];
      const UF_long plen = pe - pb, clen = ce - cb;
      joins = plen == clen + 1u && plen > 0u && uc[pb] == i &&
              (clen == 0u ||
               memcmp(uc + pb + 1u, uc + cb,
                      (size_t)clen * sizeof(*uc)) == 0);
    }
    if (!joins) {
      for (UF_long r = start; r < i; ++r) {
        sn_end[r] = i - 1u;
      }
      start = i;
    }
  }
  solver->row_refactor_sn_end = sn_end;
}

/* Short U-row supernodes are an executor property, not an input shape.  The
   candidate gate uses only the already accepted numeric and retained-row
   lifecycle.  The builder below then measures how much of the real L stream
   can use the fused executor before publishing the capability. */
static int kls_generic_hoisted_snode_worker_candidate(
  const kls_solver *solver) {
  if (solver == NULL ||
      getenv("KLS_DISABLE_GENERIC_HOISTED_SNODE_WORKER") != NULL ||
      solver->symbolic == NULL || solver->numeric == NULL ||
      !kls_repeated_update_workload(&solver->options) ||
      solver->options.threads < 2 || solver->common.scale > 0 ||
      solver->numeric->Rs != NULL ||
      solver->pivot_nudge_count != 0u ||
      solver->common.kls_perturb_count != 0u ||
      solver->numeric->lnz > UF_long_max - solver->numeric->unz) {
    return 0;
  }
  const UF_long factor_entries =
    solver->numeric->lnz + solver->numeric->unz;
  return factor_entries > 0u && solver->common.flops > 0.0 &&
    solver->common.flops >= 4.0 * (double)factor_entries;
}

static int kls_generic_hoisted_snode_worker_capable(
  const kls_solver *solver) {
  return kls_generic_hoisted_snode_worker_candidate(solver) &&
    solver->row_refactor_l_ptr != NULL &&
    solver->row_refactor_l_cols != NULL &&
    solver->row_refactor_l_values != NULL &&
    solver->row_refactor_u_ptr != NULL &&
    solver->row_refactor_u_cols != NULL &&
    solver->row_refactor_u_values != NULL &&
    solver->row_refactor_u_row_values != NULL &&
    solver->row_refactor_input_ptr != NULL &&
    solver->row_refactor_input_cols != NULL &&
    solver->row_refactor_input_pos != NULL &&
    solver->row_refactor_level_rows != NULL &&
    solver->row_refactor_level_max_width >=
      (UF_long)(2 * solver->options.threads) &&
    solver->lean_snode_run != NULL &&
    solver->lean_snode_worker_eligible > 0;
}

static void kls_lean_row_refactor_ensure_snode_runs(kls_solver *solver) {
  if (solver == NULL || solver->lean_snode_run != NULL ||
      solver->lean_snode_worker_eligible != 0 ||
      solver->numeric == NULL ||
      solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_l_cols == NULL) {
    return;
  }
  solver->lean_snode_worker_eligible = -1;
  solver->lean_snode_avx512_candidate = 0;
  kls_row_refactor_ensure_sn_end(solver);
  if (solver->row_refactor_sn_end == NULL) {
    solver->lean_snode_worker_eligible = 0;
    return;
  }
  const UF_long l_nnz = solver->row_refactor_l_ptr[solver->n];
  if (l_nnz == 0u || (uintmax_t)l_nnz > (uintmax_t)SIZE_MAX) {
    return;
  }
  unsigned char *runs = (unsigned char *)calloc((size_t)l_nnz,
                                                 sizeof(*runs));
  if (runs == NULL) {
    solver->lean_snode_worker_eligible = 0;
    return;
  }
  const UF_long *l_ptr = solver->row_refactor_l_ptr;
  const UF_long *l_cols = solver->row_refactor_l_cols;
  const UF_long *sn_end = solver->row_refactor_sn_end;
  const int only_width2 = getenv("KLS_LEAN_SNODE2_ONLY") != NULL;
  UF_long run_count[17] = {0};
  UF_long run_targets[17] = {0};
  UF_long covered_entries = 0u;
  UF_long fused_runs = 0u;
#if KLS_HAVE_AVX512_KERNELS
  const int width16_enabled = !only_width2 && kls_avx512_scatter_enabled();
#else
  const int width16_enabled = 0;
#endif
  for (UF_long row = 0u; row < solver->n; ++row) {
    UF_long p = l_ptr[row];
    const UF_long p_end = l_ptr[row + 1u];
    while (p < p_end) {
      const UF_long dep0 = l_cols[p];
      UF_long span = 1u;
      const UF_long span_limit = width16_enabled ? 16u : 8u;
      while (span < span_limit && span < p_end - p &&
             dep0 + span <= sn_end[dep0] &&
             l_cols[p + span] == dep0 + span) {
        span++;
      }
      /* A full width-16 tile can keep the common trailing workspace values
         live across two width-eight arithmetic groups.  Partial tiles retain
         the established 2--8 kernels, avoiding a dynamic-width hot loop. */
      const UF_long run = span == 16u ? 16u : (span < 8u ? span : 8u);
      const int complete = run >= 2u && (!only_width2 || run == 2u);
      if (complete) {
        const UF_long dep_end = dep0 + run - 1u;
        runs[p] = (unsigned char)run;
        run_count[run]++;
        covered_entries += run;
        fused_runs++;
        const UF_long target_length =
          solver->row_refactor_u_ptr[dep_end + 1u] -
          solver->row_refactor_u_ptr[dep_end];
        run_targets[run] += target_length;
        if ((run == 8u || run == 16u) &&
            target_length >= KLS_LEAN_SNODE_AVX512_MIN_TARGETS) {
          solver->lean_snode_avx512_candidate = 1;
        }
        p += run;
      } else {
        p++;
      }
    }
  }
  if (run_count[16] != 0u) {
    solver->lean_snode_avx512_candidate = 2;
  }
  /* The accepted path has its own hoisted worker.  The established broad
     contract still requires two-thirds coverage; compute-dense factors can
     amortize the run checks from one-third coverage.  Sparser run sets
     discard the descriptors and stay on the ordinary hoisted worker. */
  const UF_long broad_coverage = l_nnz - l_nnz / 3u;
  const UF_long dense_coverage =
    l_nnz / 3u + (l_nnz % 3u != 0u);
  const double factor_entries =
    (double)solver->numeric->lnz + (double)solver->numeric->unz;
  const int dense_work = factor_entries > 0.0 &&
    solver->common.flops >= 16.0 * factor_entries;
  const int eligible =
    fused_runs >= (UF_long)solver->options.threads &&
    (covered_entries >= broad_coverage ||
     (dense_work && covered_entries >= dense_coverage)) ? 1 : -1;
  solver->lean_snode_worker_eligible = eligible;
  if (eligible < 0 && getenv("KLS_ENABLE_LEAN_SNODE") == NULL) {
    free(runs);
    runs = NULL;
  }
  solver->lean_snode_run = runs;
  if (getenv("KLS_TRACE_LEAN_SNODE") != NULL) {
    fprintf(stderr,
            "KLS lean snode runs: 2=%ld/%ld 3=%ld/%ld 4=%ld/%ld"
            " 5=%ld/%ld 6=%ld/%ld 7=%ld/%ld 8=%ld/%ld"
            " 16=%ld/%ld coverage=%ld/%ld"
            " work=%.0f/%.0f dense=%d eligible=%d\n",
            (long)run_count[2], (long)run_targets[2],
            (long)run_count[3], (long)run_targets[3],
            (long)run_count[4], (long)run_targets[4],
            (long)run_count[5], (long)run_targets[5],
            (long)run_count[6], (long)run_targets[6],
            (long)run_count[7], (long)run_targets[7],
            (long)run_count[8], (long)run_targets[8],
            (long)run_count[16], (long)run_targets[16],
            (long)covered_entries, (long)l_nnz,
            solver->common.flops, factor_entries, dense_work, eligible);
  }
}

/* Rows owned by one pool worker are executed in increasing schedule-slot
   order.  Acquiring one producer slot therefore covers every earlier row
   from that producer for the rest of the consumer stream.  Predecode the
   advancing producer frontier once per retained schedule: the numeric loop
   consumes only compact grouped slots and producers publish only frontiers
   that a remote worker will actually acquire. */
static int kls_prepare_lean_snode_wait_masks(kls_solver *solver,
                                             int thread_count,
                                             const UF_long *rows) {
  if (solver == NULL) {
    return 0;
  }
  if (solver->lean_snode_wait_mask != NULL &&
      solver->lean_snode_wait_slot != NULL &&
      solver->lean_snode_wait_rows == rows &&
      solver->lean_snode_wait_thread_count == thread_count) {
    return 1;
  }
  free(solver->lean_snode_wait_mask);
  free(solver->lean_snode_wait_slot);
  solver->lean_snode_wait_mask = NULL;
  solver->lean_snode_wait_slot = NULL;
  solver->lean_snode_wait_rows = NULL;
  solver->lean_snode_wait_thread_count = 0;
  if (rows == NULL || thread_count < 2 || thread_count > 8 ||
      solver->lean_snode_run == NULL ||
      solver->lean_parallel_grouped_token == NULL ||
      solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_l_cols == NULL) {
    return 0;
  }
  const UF_long l_count = solver->row_refactor_l_ptr[solver->n];
  if (l_count == 0u ||
      (uintmax_t)l_count > (uintmax_t)PTRDIFF_MAX ||
      (uintmax_t)l_count >
        (uintmax_t)(SIZE_MAX / sizeof(uint32_t))) {
    return 0;
  }
  unsigned char *masks =
    (unsigned char *)calloc((size_t)l_count, sizeof(*masks));
  uint32_t *wait_slots =
    (uint32_t *)calloc((size_t)l_count, sizeof(*wait_slots));
  if (masks == NULL || wait_slots == NULL) {
    free(masks);
    free(wait_slots);
    return 0;
  }
  uint32_t acquired_slot[8][8] = {{0u}};
  for (UF_long row = 0u; row < solver->n; ++row) {
    solver->lean_parallel_grouped_token[row] &=
      ~KLS_LEAN_GROUPED_PUBLISH_FLAG;
  }
  int valid = 1;
  UF_long wait_count = 0u;
  UF_long publish_count = 0u;
  for (UF_long pos = 0u; pos < solver->n && valid; ++pos) {
    const UF_long row = rows[pos];
    if (row >= solver->n) {
      valid = 0;
      break;
    }
    const uint32_t consumer_token =
      solver->lean_parallel_grouped_token[row];
    const unsigned int consumer =
      consumer_token >> KLS_LEAN_GROUPED_OWNER_SHIFT;
    if (consumer >= (unsigned int)thread_count) {
      valid = 0;
      break;
    }
    UF_long p = solver->row_refactor_l_ptr[row];
    const UF_long p_end = solver->row_refactor_l_ptr[row + 1u];
    while (p < p_end) {
      const UF_long run = (UF_long)solver->lean_snode_run[p];
      if (run < 2u || (run > 8u && run != 16u) || p + run > p_end) {
        const UF_long dep = solver->row_refactor_l_cols[p];
        if (dep >= solver->n) {
          valid = 0;
          break;
        }
        const uint32_t token = solver->lean_parallel_grouped_token[dep];
        const unsigned int owner =
          token >> KLS_LEAN_GROUPED_OWNER_SHIFT;
        if (owner >= (unsigned int)thread_count) {
          valid = 0;
          break;
        }
        if (owner != consumer) {
          const uint32_t encoded_slot =
            (token & KLS_LEAN_GROUPED_SLOT_MASK) + 1u;
          if (encoded_slot > acquired_slot[consumer][owner]) {
            wait_slots[p] = encoded_slot;
            acquired_slot[consumer][owner] = encoded_slot;
            if ((solver->lean_parallel_grouped_token[dep] &
                 KLS_LEAN_GROUPED_PUBLISH_FLAG) == 0u) {
              publish_count++;
            }
            solver->lean_parallel_grouped_token[dep] |=
              KLS_LEAN_GROUPED_PUBLISH_FLAG;
            wait_count++;
          }
        }
        p++;
        continue;
      }
      /* Preserve the established dependency timing for a width-16 tile:
         acquire its first eight producers, form the intermediate values,
         then acquire the second eight.  This avoids extending a worker's
         blocking frontier and therefore cannot introduce a new static-owner
         rendezvous cycle. */
      const UF_long chunks = run == 16u ? 2u : 1u;
      for (UF_long chunk = 0u; chunk < chunks && valid; ++chunk) {
        const UF_long chunk_p = p + 8u * chunk;
        const UF_long chunk_run = run == 16u ? 8u : run;
        uint32_t owner_slots[8] = {0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u};
        unsigned char owner_locals[8] =
          {0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u};
        for (UF_long local = 0u; local < chunk_run; ++local) {
          const UF_long dep =
            solver->row_refactor_l_cols[chunk_p + local];
          if (dep >= solver->n) {
            valid = 0;
            break;
          }
          const uint32_t token = solver->lean_parallel_grouped_token[dep];
          const unsigned int owner =
            token >> KLS_LEAN_GROUPED_OWNER_SHIFT;
          if (owner >= (unsigned int)thread_count) {
            valid = 0;
            break;
          }
          if (owner == consumer) {
            continue;
          }
          const uint32_t encoded_slot =
            (token & KLS_LEAN_GROUPED_SLOT_MASK) + 1u;
          if (encoded_slot > owner_slots[owner]) {
            owner_slots[owner] = encoded_slot;
            owner_locals[owner] = (unsigned char)local;
          }
        }
        unsigned char mask = 0u;
        for (unsigned int owner = 0u;
             owner < (unsigned int)thread_count; ++owner) {
          const uint32_t encoded_slot = owner_slots[owner];
          if (encoded_slot > acquired_slot[consumer][owner]) {
            const UF_long local = (UF_long)owner_locals[owner];
            const UF_long dep =
              solver->row_refactor_l_cols[chunk_p + local];
            mask |= (unsigned char)(1u << local);
            wait_slots[chunk_p + local] = encoded_slot;
            acquired_slot[consumer][owner] = encoded_slot;
            if ((solver->lean_parallel_grouped_token[dep] &
                 KLS_LEAN_GROUPED_PUBLISH_FLAG) == 0u) {
              publish_count++;
            }
            solver->lean_parallel_grouped_token[dep] |=
              KLS_LEAN_GROUPED_PUBLISH_FLAG;
            wait_count++;
          }
        }
        masks[chunk_p] = mask;
      }
      p += run;
    }
  }
  if (!valid) {
    free(masks);
    free(wait_slots);
    for (UF_long row = 0u; row < solver->n; ++row) {
      solver->lean_parallel_grouped_token[row] |=
        KLS_LEAN_GROUPED_PUBLISH_FLAG;
    }
    return 0;
  }
  solver->lean_snode_wait_mask = masks;
  solver->lean_snode_wait_slot = wait_slots;
  solver->lean_snode_wait_rows = rows;
  solver->lean_snode_wait_thread_count = thread_count;
  if (getenv("KLS_TRACE_LEAN_SNODE") != NULL) {
    fprintf(stderr,
            "KLS lean snode synchronization: waits=%ld publishes=%ld/%ld\n",
            (long)wait_count, (long)publish_count, (long)solver->n);
  }
  return 1;
}

static int kls_finish_lean_row_refactor_pattern_from_arrays(
  kls_solver *solver,
  UF_long *l_ptr,
  UF_long *l_cols,
  double **l_values,
  double *l_row_values,
  UF_long *u_ptr,
  UF_long *u_cols,
  double **u_values,
  double *u_row_values,
  UF_long *input_ptr,
  UF_long *input_cols,
  UF_long *input_pos,
  uint32_t *packed_input,
  UF_long *next) {
  const UF_long n = solver != NULL ? solver->n : 0u;
  const UF_long l_nnz = l_ptr != NULL ? l_ptr[n] : 0u;
  const UF_long u_nnz = u_ptr != NULL ? u_ptr[n] : 0u;
  const UF_long input_nnz = input_ptr != NULL ? input_ptr[n] : 0u;
  if (solver == NULL || l_ptr == NULL || u_ptr == NULL || input_ptr == NULL ||
      (l_nnz > 0u &&
       (l_cols == NULL || l_values == NULL || l_row_values == NULL)) ||
      (u_nnz > 0u &&
       (u_cols == NULL || u_values == NULL || u_row_values == NULL)) ||
      (input_nnz > 0u && (input_cols == NULL || input_pos == NULL))) {
    free(l_ptr);
    free(l_cols);
    free(l_values);
    free(l_row_values);
    free(u_ptr);
    free(u_cols);
    free(u_values);
    free(u_row_values);
    free(input_ptr);
    free(input_cols);
    free(input_pos);
    free(packed_input);
    free(next);
    return 0;
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
  solver->row_refactor_input_col_user32 = packed_input;
  solver->row_refactor_pattern_n = n;

  UF_long *levels = n > 0u
    ? (UF_long *)calloc((size_t)n, sizeof(*levels)) : NULL;
  if (n > 0u && levels == NULL) {
    free_row_refactor_pattern(solver);
    return 0;
  }
  UF_long max_level = 0u;
  for (UF_long row = 0u; row < n; ++row) {
    UF_long level = 0u;
    for (UF_long p = l_ptr[row]; p < l_ptr[row + 1u]; ++p) {
      const UF_long dep = l_cols[p];
      if (dep >= row) {
        free(levels);
        free_row_refactor_pattern(solver);
        return 0;
      }
      const UF_long candidate = levels[dep] + 1u;
      if (candidate > level) {
        level = candidate;
      }
    }
    levels[row] = level;
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
  for (UF_long row = 0u; row < n; ++row) {
    level_ptr[levels[row] + 1u]++;
  }
  for (UF_long level = 0u; level < level_count; ++level) {
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
  if (level_count > 0u) {
    memcpy(write_pos, level_ptr,
           (size_t)level_count * sizeof(*write_pos));
  }
  for (UF_long row = 0u; row < n; ++row) {
    level_rows[write_pos[levels[row]]++] = row;
  }
  UF_long max_width = 0u;
  for (UF_long level = 0u; level < level_count; ++level) {
    const UF_long width = level_ptr[level + 1u] - level_ptr[level];
    if (width > max_width) {
      max_width = width;
    }
  }
  free(write_pos);
  free(levels);
  solver->row_refactor_level_ptr = level_ptr;
  solver->row_refactor_level_rows = level_rows;
  solver->row_refactor_level_count = level_count;
  solver->row_refactor_level_max_width = max_width;
  return 1;
}

static int kls_lean_packed_input_is_profitable(const kls_solver *solver,
                                               UF_long input_nnz) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      input_nnz == 0u) {
    return 0;
  }

  /* Packing replaces two full-width loads in every retained-row input
     scatter with one 32-bit load.  Once the caller has proved that both
     fields fit, the decode is a mask and shift while the unpacked stream
     consumes four times as much index bandwidth.  The old arithmetic-ratio
     gate rejected packing precisely when input scatter was a material share
     of the recurring numeric (for example, shallow circuit factors).  Keep
     only a small stream-length floor so construction remains worthwhile. */
  return input_nnz >= 1024u;
}

static void kls_lean_pattern_worker_run(
  kls_egraph_refactor_worker *worker) {
  if (worker == NULL || worker->shared == NULL || worker->x == NULL) {
    return;
  }
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  kls_lean_pattern_job *job = shared->lean_pattern_job;
  if (solver == NULL || job == NULL || solver->symbolic == NULL ||
      solver->numeric == NULL || shared->thread_count < 2 ||
      worker->tid < 0 || worker->tid >= shared->thread_count) {
    kls_egraph_refactor_record_invalid(shared);
    return;
  }

  UF_long *const default_cursor = (UF_long *)(void *)worker->x;
  const UF_long n = solver->n;
  const int all_phase =
    job->phase == KLS_LEAN_PATTERN_COUNT_ALL ||
    job->phase == KLS_LEAN_PATTERN_FILL_ALL;
  UF_long *const l_cursor = default_cursor;
  UF_long *const u_cursor = all_phase
    ? job->u_cursors + (size_t)worker->tid * (size_t)n
    : default_cursor;
  UF_long *const input_cursor = all_phase
    ? job->input_cursors + (size_t)worker->tid * (size_t)n
    : default_cursor;
  const int count_phase =
    job->phase == KLS_LEAN_PATTERN_COUNT_L ||
    job->phase == KLS_LEAN_PATTERN_COUNT_U ||
    job->phase == KLS_LEAN_PATTERN_COUNT_INPUT ||
    job->phase == KLS_LEAN_PATTERN_COUNT_ALL;
  if (count_phase) {
    memset(l_cursor, 0, (size_t)n * sizeof(*l_cursor));
    if (all_phase) {
      memset(u_cursor, 0, (size_t)n * sizeof(*u_cursor));
      memset(input_cursor, 0, (size_t)n * sizeof(*input_cursor));
    }
  }

  /* Contiguous column ownership makes the per-thread row segments globally
     sorted when the main thread concatenates them in tid order.  The numeric
     row walk relies on increasing predecessor columns. */
  const UF_long k_begin =
    n * (UF_long)worker->tid / (UF_long)shared->thread_count;
  const UF_long k_end =
    n * (UF_long)(worker->tid + 1) / (UF_long)shared->thread_count;
  UF_long direct_block = 0u;
  if (job->direct_numeric) {
    while (direct_block + 1u < solver->symbolic->nblocks &&
           solver->symbolic->R[direct_block + 1u] <= k_begin) {
      direct_block++;
    }
  }
  for (UF_long k = k_begin; k < k_end; ++k) {
    UF_long block = 0u;
    if (job->direct_numeric) {
      while (direct_block + 1u < solver->symbolic->nblocks &&
             solver->symbolic->R[direct_block + 1u] <= k) {
        direct_block++;
      }
      block = direct_block;
    } else if (solver->symbolic->nblocks > 1u) {
      if (solver->refactor_col_block == NULL) {
        kls_egraph_refactor_record_invalid(shared);
        return;
      }
      block = solver->refactor_col_block[k];
    }
    if (block >= solver->symbolic->nblocks) {
      kls_egraph_refactor_record_invalid(shared);
      return;
    }
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    if (k < k1 || k >= k2) {
      kls_egraph_refactor_record_invalid(shared);
      return;
    }

    if (job->phase == KLS_LEAN_PATTERN_COUNT_L ||
        job->phase == KLS_LEAN_PATTERN_FILL_L || all_phase) {
      UF_long *cursor = l_cursor;
      UF_long *out_cols = all_phase ? job->all_l_cols : job->cols;
      double **out_values = all_phase ? job->all_l_values : job->values;
      const UF_long capacity = all_phase
        ? job->all_l_capacity : job->capacity;
      const UF_long len =
        kls_numeric_l_len_in_block(solver, k1, k2, k);
      UF_long *direct_rows = NULL;
      double *direct_values = NULL;
      UF_long direct_len = 0u;
      if (job->direct_numeric && len > 0u) {
        double *lu = (double *)solver->numeric->LUbx[block];
        if (lu == NULL) {
          kls_egraph_refactor_record_invalid(shared);
          return;
        }
        kls_klu_get_pointer(
          lu, solver->numeric->Lip + k1, solver->numeric->Llen + k1,
          k - k1, &direct_rows, &direct_values, &direct_len);
      }
      const UF_long *rows = job->direct_numeric
        ? direct_rows : solver->refactor_l_indices[k];
      double *values = job->direct_numeric
        ? direct_values : solver->refactor_l_values[k];
      if (job->direct_numeric && direct_len != len) {
        kls_egraph_refactor_record_invalid(shared);
        return;
      }
      if (len > 0u && (rows == NULL || values == NULL)) {
        kls_egraph_refactor_record_invalid(shared);
        return;
      }
      for (UF_long p = 0u; p < len; ++p) {
        const UF_long row = k1 + rows[p];
        if (row <= k || row >= k2) {
          kls_egraph_refactor_record_invalid(shared);
          return;
        }
        if (count_phase) {
          cursor[row]++;
        } else {
          const UF_long dst = cursor[row]++;
          if (dst >= capacity) {
            kls_egraph_refactor_record_invalid(shared);
            return;
          }
          out_cols[dst] = k;
          out_values[dst] = values + p;
        }
      }
    }
    if (job->phase == KLS_LEAN_PATTERN_COUNT_U ||
        job->phase == KLS_LEAN_PATTERN_FILL_U || all_phase) {
      UF_long *cursor = u_cursor;
      UF_long *out_cols = all_phase ? job->all_u_cols : job->cols;
      double **out_values = all_phase ? job->all_u_values : job->values;
      const UF_long capacity = all_phase
        ? job->all_u_capacity : job->capacity;
      const UF_long len =
        kls_numeric_u_len_in_block(solver, k1, k2, k);
      UF_long *direct_rows = NULL;
      double *direct_values = NULL;
      UF_long direct_len = 0u;
      if (job->direct_numeric && len > 0u) {
        double *lu = (double *)solver->numeric->LUbx[block];
        if (lu == NULL) {
          kls_egraph_refactor_record_invalid(shared);
          return;
        }
        kls_klu_get_pointer(
          lu, solver->numeric->Uip + k1, solver->numeric->Ulen + k1,
          k - k1, &direct_rows, &direct_values, &direct_len);
      }
      const UF_long *rows = job->direct_numeric
        ? direct_rows : solver->refactor_u_indices[k];
      double *values = job->direct_numeric
        ? direct_values : solver->refactor_u_values[k];
      if (job->direct_numeric && direct_len != len) {
        kls_egraph_refactor_record_invalid(shared);
        return;
      }
      if (len > 0u && (rows == NULL || values == NULL)) {
        kls_egraph_refactor_record_invalid(shared);
        return;
      }
      for (UF_long p = 0u; p < len; ++p) {
        const UF_long row = k1 + rows[p];
        if (row >= k || row < k1) {
          kls_egraph_refactor_record_invalid(shared);
          return;
        }
        if (count_phase) {
          cursor[row]++;
        } else {
          const UF_long dst = cursor[row]++;
          if (dst >= capacity) {
            kls_egraph_refactor_record_invalid(shared);
            return;
          }
          out_cols[dst] = k;
          out_values[dst] = values + p;
        }
      }
    }
    if (job->phase == KLS_LEAN_PATTERN_COUNT_INPUT ||
        job->phase == KLS_LEAN_PATTERN_FILL_INPUT || all_phase) {
      UF_long *cursor = input_cursor;
      UF_long *out_cols = all_phase ? job->all_input_cols : job->cols;
      UF_long *out_input_pos = all_phase
        ? job->all_input_pos : job->input_pos;
      const UF_long capacity = all_phase
        ? job->all_input_capacity : job->capacity;
      if (job->direct_numeric) {
        const UF_long oldcol = solver->symbolic->Q[k];
        if (oldcol >= n) {
          kls_egraph_refactor_record_invalid(shared);
          return;
        }
        const UF_long begin = solver->col_ptr[oldcol];
        const UF_long end = solver->col_ptr[oldcol + 1u];
        if (begin > end || end > solver->nnz) {
          kls_egraph_refactor_record_invalid(shared);
          return;
        }
        UF_long poff = solver->numeric->Offp[k];
        const UF_long poff_end = solver->numeric->Offp[k + 1u];
        for (UF_long p = begin; p < end; ++p) {
          const UF_long oldrow = solver->row_idx[p];
          if (oldrow >= n) {
            kls_egraph_refactor_record_invalid(shared);
            return;
          }
          const UF_long row = solver->numeric->Pinv[oldrow];
          if (row >= k2) {
            kls_egraph_refactor_record_invalid(shared);
            return;
          }
          if (row < k1) {
            if (poff >= poff_end ||
                (job->offdiag_input_pos32 != NULL &&
                 p > (UF_long)UINT32_MAX) ||
                (job->offdiag_input_pos32 == NULL &&
                 p > (UF_long)UINT16_MAX)) {
              kls_egraph_refactor_record_invalid(shared);
              return;
            }
            if (!count_phase) {
              if (job->offdiag_input_pos32 != NULL) {
                job->offdiag_input_pos32[poff] = (uint32_t)p;
              } else if (job->offdiag_input_pos16 != NULL) {
                job->offdiag_input_pos16[poff] = (uint16_t)p;
              }
            }
            poff++;
            continue;
          }
          if (count_phase) {
            cursor[row]++;
          } else {
            const UF_long dst = cursor[row]++;
            if (dst >= capacity) {
              kls_egraph_refactor_record_invalid(shared);
              return;
            }
            out_cols[dst] = k;
            out_input_pos[dst] = p;
            if (job->packed_input != NULL &&
                k <= (UF_long)UINT16_MAX &&
                p <= (UF_long)UINT16_MAX) {
              job->packed_input[dst] =
                ((uint32_t)p << 16u) | (uint32_t)k;
            }
          }
        }
        if (poff != poff_end) {
          kls_egraph_refactor_record_invalid(shared);
          return;
        }
      } else {
        UF_long begin = solver->refactor_col_ptr[k];
        if (solver->symbolic->nblocks > 1u) {
          begin = solver->refactor_block_start[k];
        }
        const UF_long end = solver->refactor_col_ptr[k + 1u];
        if (begin > end) {
          kls_egraph_refactor_record_invalid(shared);
          return;
        }
        for (UF_long p = begin; p < end; ++p) {
          const UF_long row = solver->refactor_row_idx[p];
          if (row >= n || solver->refactor_input_pos[p] >= solver->nnz) {
            kls_egraph_refactor_record_invalid(shared);
            return;
          }
          if (count_phase) {
            cursor[row]++;
          } else {
            const UF_long dst = cursor[row]++;
            if (dst >= capacity) {
              kls_egraph_refactor_record_invalid(shared);
              return;
            }
            out_cols[dst] = k;
            out_input_pos[dst] = solver->refactor_input_pos[p];
            if (job->packed_input != NULL &&
                k <= (UF_long)UINT16_MAX &&
                solver->refactor_input_pos[p] <= (UF_long)UINT16_MAX) {
              job->packed_input[dst] =
                ((uint32_t)solver->refactor_input_pos[p] << 16u) |
                (uint32_t)k;
            }
          }
        }
      }
    }
  }
}

static int kls_lean_pattern_prepare_offsets(double **scratch,
                                            int thread_count,
                                            UF_long n,
                                            UF_long *ptr) {
  if (scratch == NULL || thread_count < 2 || ptr == NULL) {
    return 0;
  }
  ptr[0] = 0u;
  for (UF_long row = 0u; row < n; ++row) {
    UF_long total = 0u;
    for (int tid = 0; tid < thread_count; ++tid) {
      const UF_long count = ((UF_long *)(void *)scratch[tid])[row];
      if (count > UF_long_max - total) {
        return 0;
      }
      total += count;
    }
    if (total > UF_long_max - ptr[row]) {
      return 0;
    }
    ptr[row + 1u] = ptr[row] + total;
  }
  for (UF_long row = 0u; row < n; ++row) {
    UF_long offset = ptr[row];
    for (int tid = 0; tid < thread_count; ++tid) {
      UF_long *cursor = (UF_long *)(void *)scratch[tid];
      const UF_long count = cursor[row];
      cursor[row] = offset;
      offset += count;
    }
    if (offset != ptr[row + 1u]) {
      return 0;
    }
  }
  return 1;
}

static int kls_lean_pattern_prepare_offsets_base(
  UF_long *scratch,
  int thread_count,
  UF_long n,
  UF_long *ptr) {
  if (scratch == NULL || thread_count < 2 || ptr == NULL ||
      n > (UF_long)(SIZE_MAX / sizeof(*scratch)) ||
      (size_t)n > SIZE_MAX / (size_t)thread_count) {
    return 0;
  }
  ptr[0] = 0u;
  for (UF_long row = 0u; row < n; ++row) {
    UF_long total = 0u;
    for (int tid = 0; tid < thread_count; ++tid) {
      const UF_long count = scratch[(size_t)tid * (size_t)n + row];
      if (count > UF_long_max - total) return 0;
      total += count;
    }
    if (total > UF_long_max - ptr[row]) return 0;
    ptr[row + 1u] = ptr[row] + total;
  }
  for (UF_long row = 0u; row < n; ++row) {
    UF_long offset = ptr[row];
    for (int tid = 0; tid < thread_count; ++tid) {
      UF_long *cursor = scratch + (size_t)tid * (size_t)n;
      const UF_long count = cursor[row];
      cursor[row] = offset;
      offset += count;
    }
    if (offset != ptr[row + 1u]) return 0;
  }
  return 1;
}

static int kls_direct_numeric_lean_pattern_capable(
  const kls_solver *solver);

/* A direct row pattern whose complete immutable streams fit the compact
   worker representation benefits from using the same bounded crew while it
   is built and consumed.  Besides avoiding a worker-pool rebuild, this keeps
   the first and recurring passes on the same scratch/cache topology.  The
   admission is entirely representation based: the direct builder still
   validates every index before publishing the packed pattern. */
static int kls_compact_direct_numeric_row_pattern_capable(
  const kls_solver *solver) {
  return kls_direct_numeric_lean_pattern_capable(solver) &&
    solver->common.scale <= 0 && solver->numeric->Rs == NULL &&
    solver->n <= (UF_long)0x1000u &&
    solver->nnz <= (UF_long)UINT16_MAX &&
    solver->numeric->lnz <= (UF_long)UINT16_MAX &&
    solver->numeric->unz <= (UF_long)UINT16_MAX;
}

static int kls_build_lean_row_refactor_pattern_parallel(
  kls_solver *solver,
  int direct_numeric,
  int lean_only) {
  const int trace_pattern =
    getenv("KLS_TRACE_DIRECT_ROW_PATTERN") != NULL;
  const double trace_begin = trace_pattern ? kls_now_seconds() : 0.0;
  if (solver == NULL || solver->options.threads < 2 || solver->n == 0u ||
      solver->symbolic == NULL || solver->numeric == NULL ||
      (direct_numeric &&
       (solver->col_ptr == NULL || solver->row_idx == NULL ||
        solver->symbolic->Q == NULL || solver->symbolic->R == NULL ||
        solver->numeric->Pinv == NULL || solver->numeric->LUbx == NULL ||
        solver->numeric->Offp == NULL || solver->numeric->Offx == NULL ||
        solver->numeric->Lip == NULL || solver->numeric->Llen == NULL ||
        solver->numeric->Uip == NULL || solver->numeric->Ulen == NULL)) ||
      (!direct_numeric &&
       (solver->refactor_col_ptr == NULL ||
        solver->refactor_row_idx == NULL ||
        solver->refactor_input_pos == NULL ||
        solver->refactor_l_indices == NULL ||
        solver->refactor_l_values == NULL ||
        solver->refactor_u_indices == NULL ||
        solver->refactor_u_values == NULL))) {
    return 0;
  }
  const int compact_bounded_crew =
    (direct_numeric && lean_only &&
     kls_compact_direct_numeric_row_pattern_capable(solver));
  const int thread_count =
    compact_bounded_crew && solver->options.threads > 5
      ? 5 : solver->options.threads;
  double **scratch =
    ensure_egraph_worker_scratch(solver, thread_count, solver->n);
  kls_egraph_refactor_pool *pool =
    ensure_egraph_refactor_pool(solver, thread_count);
  if (scratch == NULL || pool == NULL) {
    return 0;
  }
  solver->lean_parallel_scratch_clean = 0;

  free_row_refactor_pattern(solver);
  const UF_long n = solver->n;
  UF_long *l_ptr = (UF_long *)calloc((size_t)n + 1u, sizeof(*l_ptr));
  UF_long *u_ptr = (UF_long *)calloc((size_t)n + 1u, sizeof(*u_ptr));
  UF_long *input_ptr =
    (UF_long *)calloc((size_t)n + 1u, sizeof(*input_ptr));
  UF_long *l_cols = NULL;
  double **l_values = NULL;
  double *l_row_values = NULL;
  UF_long *u_cols = NULL;
  double **u_values = NULL;
  double *u_row_values = NULL;
  UF_long *input_cols = NULL;
  UF_long *input_pos = NULL;
  uint32_t *packed_input = NULL;
  uint16_t *offdiag_input_pos16 = NULL;
  uint32_t *offdiag_input_pos32 = NULL;
  if (direct_numeric) {
    const UF_long offcount = solver->numeric->Offp[n];
    const int single_block = solver->symbolic->nblocks == 1u;
    const int mapped_user_frame =
      kls_direct_user_value_maps_capable(solver);
    /* The compact off-diagonal map is only needed for a multi-block direct
       pattern.  A single-block factor has no Offx entries, while its row
       input positions remain full-width UF_long values.  The old global-nnz
       guard therefore rejected large single-block factors even though no
       16-bit value was stored. */
    if ((single_block && offcount != 0u) ||
        (!single_block &&
         (solver->nnz > (UF_long)UINT32_MAX ||
          offcount > (UF_long)UINT32_MAX))) {
      free(l_ptr);
      free(u_ptr);
      free(input_ptr);
      return 0;
    }
    if (!single_block &&
        (mapped_user_frame || solver->input_to_csc != NULL ||
         solver->common.scale > 0 ||
         solver->nnz > (UF_long)UINT16_MAX + 1u ||
         offcount > (UF_long)UINT16_MAX)) {
      offdiag_input_pos32 = (uint32_t *)malloc(
        (size_t)(offcount > 0u ? offcount : 1u) *
        sizeof(*offdiag_input_pos32));
    } else {
      offdiag_input_pos16 = (uint16_t *)malloc(
        (size_t)(offcount > 0u ? offcount : 1u) *
        sizeof(*offdiag_input_pos16));
    }
  }
  if (l_ptr == NULL || u_ptr == NULL || input_ptr == NULL) {
    free(l_ptr);
    free(u_ptr);
    free(input_ptr);
    free(offdiag_input_pos16);
    free(offdiag_input_pos32);
    return 0;
  }
  if (direct_numeric && offdiag_input_pos16 == NULL &&
      offdiag_input_pos32 == NULL) {
    free(l_ptr);
    free(u_ptr);
    free(input_ptr);
    return 0;
  }
  UF_long *u_cursor_storage = NULL;
  UF_long *input_cursor_storage = NULL;
  int combined_pattern = 0;
  if (direct_numeric && solver->symbolic->nblocks > 1u &&
      n <= (UF_long)(SIZE_MAX / sizeof(UF_long)) &&
      (size_t)n <= SIZE_MAX / (size_t)thread_count) {
    const size_t cursor_count = (size_t)n * (size_t)thread_count;
    if (cursor_count <= SIZE_MAX / sizeof(UF_long)) {
      u_cursor_storage = (UF_long *)malloc(
        cursor_count * sizeof(*u_cursor_storage));
      input_cursor_storage = (UF_long *)malloc(
        cursor_count * sizeof(*input_cursor_storage));
      combined_pattern =
        u_cursor_storage != NULL && input_cursor_storage != NULL;
    }
  }
  if (!combined_pattern) {
    free(u_cursor_storage);
    free(input_cursor_storage);
    u_cursor_storage = NULL;
    input_cursor_storage = NULL;
  }

  kls_egraph_refactor_shared *shared = &pool->shared;
  kls_lean_pattern_job job;
  memset(&job, 0, sizeof(job));
  job.direct_numeric = direct_numeric;
  job.offdiag_input_pos16 = offdiag_input_pos16;
  job.offdiag_input_pos32 = offdiag_input_pos32;
  job.u_cursors = u_cursor_storage;
  job.input_cursors = input_cursor_storage;
  int ok = 1;
  pthread_mutex_lock(&shared->lock);
  if (pool->active_workers != 0) {
    ok = 0;
  } else {
    shared->solver = solver;
    shared->thread_count = thread_count;
    shared->lean_pattern_job = &job;
    shared->lean_pattern_mode = 1;
    shared->lean_refactor_mode = 0;
    shared->row_refactor_mode = 0;
    shared->row_solve_mode = 0;
    shared->row_publish_mode = 0;
    shared->invalid = 0;
    atomic_store_explicit(&shared->stop, 0, memory_order_release);
    for (int tid = 0; tid < thread_count; ++tid) {
      pool->workers[tid].shared = shared;
      pool->workers[tid].x = scratch[tid];
    }
  }

#define KLS_LEAN_PATTERN_DISPATCH(phase_value)                            \
  do {                                                                    \
    if (ok) {                                                             \
      job.phase = (phase_value);                                          \
      kls_egraph_pool_dispatch_and_spin_wait(pool, shared, thread_count); \
      ok = !shared->invalid;                                              \
    }                                                                     \
  } while (0)

  const double trace_dispatch_begin =
    trace_pattern ? kls_now_seconds() : 0.0;
  double trace_l_count = trace_dispatch_begin;
  double trace_l_fill = trace_dispatch_begin;
  double trace_u_count = trace_dispatch_begin;
  double trace_u_fill = trace_dispatch_begin;
  double trace_input_count = trace_dispatch_begin;
  double trace_input_fill = trace_dispatch_begin;
  UF_long l_nnz = 0u;
  UF_long u_nnz = 0u;
  UF_long input_nnz = 0u;
  if (combined_pattern) {
    KLS_LEAN_PATTERN_DISPATCH(KLS_LEAN_PATTERN_COUNT_ALL);
    trace_l_count = trace_pattern ? kls_now_seconds() : 0.0;
    if (ok) {
      ok = kls_lean_pattern_prepare_offsets(
        scratch, thread_count, n, l_ptr) &&
        kls_lean_pattern_prepare_offsets_base(
          u_cursor_storage, thread_count, n, u_ptr) &&
        kls_lean_pattern_prepare_offsets_base(
          input_cursor_storage, thread_count, n, input_ptr);
    }
    l_nnz = ok ? l_ptr[n] : 0u;
    u_nnz = ok ? u_ptr[n] : 0u;
    input_nnz = ok ? input_ptr[n] : 0u;
    if (ok && l_nnz > 0u) {
      l_cols = (UF_long *)malloc((size_t)l_nnz * sizeof(*l_cols));
      l_values = (double **)malloc((size_t)l_nnz * sizeof(*l_values));
      l_row_values = (double *)malloc(
        (size_t)l_nnz * sizeof(*l_row_values));
      ok = l_cols != NULL && l_values != NULL && l_row_values != NULL;
    }
    if (ok && u_nnz > 0u) {
      u_cols = (UF_long *)malloc((size_t)u_nnz * sizeof(*u_cols));
      u_values = (double **)malloc((size_t)u_nnz * sizeof(*u_values));
      u_row_values = kls_aligned_double_values(u_nnz);
      ok = u_cols != NULL && u_values != NULL && u_row_values != NULL;
    }
    if (ok && input_nnz > 0u) {
      input_cols = (UF_long *)malloc(
        (size_t)input_nnz * sizeof(*input_cols));
      input_pos = (UF_long *)malloc(
        (size_t)input_nnz * sizeof(*input_pos));
      ok = input_cols != NULL && input_pos != NULL;
      if (ok && lean_only &&
          getenv("KLS_DISABLE_LEAN_PACKED_INPUT") == NULL &&
          kls_lean_packed_input_is_profitable(solver, input_nnz) &&
          solver->n <= (UF_long)UINT16_MAX + 1u &&
          solver->nnz <= (UF_long)UINT16_MAX + 1u) {
        packed_input = (uint32_t *)malloc(
          (size_t)input_nnz * sizeof(*packed_input));
      }
    }
    job.all_l_cols = l_cols;
    job.all_l_values = l_values;
    job.all_l_capacity = l_nnz;
    job.all_u_cols = u_cols;
    job.all_u_values = u_values;
    job.all_u_capacity = u_nnz;
    job.all_input_cols = input_cols;
    job.all_input_pos = input_pos;
    job.all_input_capacity = input_nnz;
    job.packed_input = packed_input;
    trace_l_fill = trace_u_count = trace_u_fill = trace_input_count =
      trace_pattern ? kls_now_seconds() : 0.0;
    KLS_LEAN_PATTERN_DISPATCH(KLS_LEAN_PATTERN_FILL_ALL);
    trace_input_fill = trace_pattern ? kls_now_seconds() : 0.0;
  } else {
  KLS_LEAN_PATTERN_DISPATCH(KLS_LEAN_PATTERN_COUNT_L);
  trace_l_count = trace_pattern ? kls_now_seconds() : 0.0;
  if (ok) {
    ok = kls_lean_pattern_prepare_offsets(scratch, thread_count, n, l_ptr);
  }
  l_nnz = ok ? l_ptr[n] : 0u;
  if (ok && l_nnz > 0u) {
    l_cols = (UF_long *)malloc((size_t)l_nnz * sizeof(*l_cols));
    l_values = (double **)malloc((size_t)l_nnz * sizeof(*l_values));
    l_row_values =
      (double *)malloc((size_t)l_nnz * sizeof(*l_row_values));
    ok = l_cols != NULL && l_values != NULL && l_row_values != NULL;
  }
  job.cols = l_cols;
  job.values = l_values;
  job.input_pos = NULL;
  job.capacity = l_nnz;
  KLS_LEAN_PATTERN_DISPATCH(KLS_LEAN_PATTERN_FILL_L);
  trace_l_fill = trace_pattern ? kls_now_seconds() : 0.0;

  KLS_LEAN_PATTERN_DISPATCH(KLS_LEAN_PATTERN_COUNT_U);
  trace_u_count = trace_pattern ? kls_now_seconds() : 0.0;
  if (ok) {
    ok = kls_lean_pattern_prepare_offsets(scratch, thread_count, n, u_ptr);
  }
  u_nnz = ok ? u_ptr[n] : 0u;
  if (ok && u_nnz > 0u) {
    u_cols = (UF_long *)malloc((size_t)u_nnz * sizeof(*u_cols));
    u_values = (double **)malloc((size_t)u_nnz * sizeof(*u_values));
    u_row_values = kls_aligned_double_values(u_nnz);
    ok = u_cols != NULL && u_values != NULL && u_row_values != NULL;
  }
  job.cols = u_cols;
  job.values = u_values;
  job.input_pos = NULL;
  job.capacity = u_nnz;
  KLS_LEAN_PATTERN_DISPATCH(KLS_LEAN_PATTERN_FILL_U);
  trace_u_fill = trace_pattern ? kls_now_seconds() : 0.0;

  KLS_LEAN_PATTERN_DISPATCH(KLS_LEAN_PATTERN_COUNT_INPUT);
  trace_input_count = trace_pattern ? kls_now_seconds() : 0.0;
  if (ok) {
    ok = kls_lean_pattern_prepare_offsets(scratch, thread_count, n,
                                          input_ptr);
  }
  input_nnz = ok ? input_ptr[n] : 0u;
  if (ok && input_nnz > 0u) {
    input_cols =
      (UF_long *)malloc((size_t)input_nnz * sizeof(*input_cols));
    input_pos =
      (UF_long *)malloc((size_t)input_nnz * sizeof(*input_pos));
    ok = input_cols != NULL && input_pos != NULL;
    if (ok && lean_only &&
        getenv("KLS_DISABLE_LEAN_PACKED_INPUT") == NULL &&
        kls_lean_packed_input_is_profitable(solver, input_nnz) &&
        solver->n <= (UF_long)UINT16_MAX + 1u &&
        solver->nnz <= (UF_long)UINT16_MAX + 1u) {
      packed_input = (uint32_t *)malloc(
        (size_t)input_nnz * sizeof(*packed_input));
    }
  }
  job.cols = input_cols;
  job.values = NULL;
  job.input_pos = input_pos;
  job.capacity = input_nnz;
  job.packed_input = packed_input;
  KLS_LEAN_PATTERN_DISPATCH(KLS_LEAN_PATTERN_FILL_INPUT);
  trace_input_fill = trace_pattern ? kls_now_seconds() : 0.0;
  }
#undef KLS_LEAN_PATTERN_DISPATCH

  shared->lean_pattern_mode = 0;
  shared->lean_pattern_job = NULL;
  pthread_mutex_unlock(&shared->lock);
  free(u_cursor_storage);
  free(input_cursor_storage);
  if (!ok) {
    free(l_ptr);
    free(l_cols);
    free(l_values);
    free(l_row_values);
    free(u_ptr);
    free(u_cols);
    free(u_values);
    free(u_row_values);
    free(input_ptr);
    free(input_cols);
    free(input_pos);
    free(packed_input);
    free(offdiag_input_pos16);
    free(offdiag_input_pos32);
    return 0;
  }
  UF_long *next = !lean_only && n > 0u
    ? (UF_long *)malloc((size_t)n * sizeof(*next)) : NULL;
  const int finished = lean_only
    ? kls_finish_lean_row_refactor_pattern_from_arrays(
        solver, l_ptr, l_cols, l_values, l_row_values,
        u_ptr, u_cols, u_values, u_row_values,
        input_ptr, input_cols, input_pos, packed_input, NULL)
    : kls_finish_row_refactor_pattern_from_arrays(
        solver, l_ptr, l_cols, l_values, l_row_values,
        u_ptr, u_cols, u_values, u_row_values,
        input_ptr, input_cols, input_pos, next);
  if (finished && direct_numeric && solver->symbolic->nblocks > 1u) {
    if (offdiag_input_pos32 != NULL) {
      free(solver->lean_btf_off_input_pos);
      free(solver->lean_btf_off_input_runs);
      solver->lean_btf_off_input_pos = offdiag_input_pos32;
      solver->lean_btf_off_input_runs = NULL;
      solver->lean_btf_off_input_run_count = 0u;
      offdiag_input_pos32 = NULL;
    } else {
      solver->compact_match_offdiag_user_pos = offdiag_input_pos16;
      offdiag_input_pos16 = NULL;
    }
  }
  free(offdiag_input_pos16);
  free(offdiag_input_pos32);
  if (trace_pattern) {
    const double trace_end = kls_now_seconds();
    fprintf(stderr,
            "KLS direct row pattern: setup=%.3f Lcount=%.3f Lfill=%.3f"
            " Ucount=%.3f Ufill=%.3f Icount=%.3f Ifill=%.3f"
            " finish=%.3f total=%.3f ms\n",
            1e3 * (trace_dispatch_begin - trace_begin),
            1e3 * (trace_l_count - trace_dispatch_begin),
            1e3 * (trace_l_fill - trace_l_count),
            1e3 * (trace_u_count - trace_l_fill),
            1e3 * (trace_u_fill - trace_u_count),
            1e3 * (trace_input_count - trace_u_fill),
            1e3 * (trace_input_fill - trace_input_count),
            1e3 * (trace_end - trace_input_fill),
            1e3 * (trace_end - trace_begin));
  }
  return finished;
}

/* The retained KLU numeric already contains the immutable L/U index streams,
   their value addresses, the numeric row permutation, and the off-diagonal
   map needed by the lean row walk.  The caller values may either share that
   CSC frame or have a retained, compact, bijective user-to-internal map;
   constructing a generic refactor map first is redundant in both cases.
   Keep this as a representation/lifecycle capability: the row engine is
   selected independently, and the direct builder still validates every
   numeric and input index before publishing the pattern. */
static int kls_direct_numeric_lean_pattern_capable(
  const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL ||
      solver->numeric == NULL || solver->options.threads <= 1 ||
      !kls_repeated_update_workload(&solver->options) ||
      getenv("KLS_DISABLE_DIRECT_NUMERIC_LEAN_PATTERN") != NULL ||
      solver->input_format != KLS_INPUT_CSC ||
      solver->user_col_perm != NULL || solver->row_scale != NULL ||
      solver->col_scale != NULL || solver->pivot_nudge_count != 0u ||
      solver->common.kls_perturb_count != 0u) {
    return 0;
  }
  const int plain_direct_input_frame =
    solver->orientation == KLS_ORIENTATION_NORMAL &&
    solver->common.scale <= 0 && solver->numeric->Rs == NULL &&
    solver->input_to_csc == NULL && solver->row_perm == NULL;
  const int mapped_user_frame =
    kls_direct_user_value_maps_capable(solver);
  /* Any selected candidate with a retained public-to-internal permutation
     owns a complete internal CSC.  The direct builder reads that CSC (not the
     caller frame), while normal value preparation gathers updates into the
     same positions before the numeric walk.  This covers a transposed input
     without detecting its matrix shape or reconstructing its public rows.
     KLU row scaling is likewise numeric metadata: input positions remain
     unchanged and the row worker applies Rs through Pnum. */
  const int mapped_internal_frame =
    solver->input_to_csc != NULL && solver->row_perm == NULL &&
    ((solver->common.scale > 0 && solver->numeric->Rs != NULL) ||
     (solver->common.scale <= 0 && solver->numeric->Rs == NULL));
  return plain_direct_input_frame || mapped_user_frame ||
    mapped_internal_frame;
}

static int kls_build_row_refactor_pattern(kls_solver *solver,
                                          int lean_only) {
  const int trace_prep = getenv("KLS_TRACE_LEAN_PREP") != NULL;
  const double trace_begin = trace_prep ? kls_now_seconds() : 0.0;
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      solver->symbolic->nblocks == 0u || solver->symbolic->R == NULL) {
    return 0;
  }
  if (lean_only && solver->row_refactor_pattern_n == solver->n &&
      solver->row_refactor_l_ptr != NULL &&
      (solver->row_refactor_l_ptr[solver->n] == 0u ||
       solver->row_refactor_l_row_values != NULL) &&
      solver->row_refactor_u_ptr != NULL &&
      (solver->row_refactor_u_ptr[solver->n] == 0u ||
       solver->row_refactor_u_row_values != NULL) &&
      solver->row_refactor_input_ptr != NULL &&
      solver->row_refactor_level_ptr != NULL &&
      solver->row_refactor_level_rows != NULL) {
    return 1;
  }
  if (solver->row_refactor_pattern_n == solver->n &&
      solver->row_refactor_l_ptr != NULL &&
      (solver->row_refactor_l_ptr[solver->n] == 0u ||
       solver->row_refactor_l_row_values != NULL) &&
      solver->row_refactor_u_ptr != NULL &&
      (solver->row_refactor_u_ptr[solver->n] == 0u ||
       solver->row_refactor_u_row_values != NULL) &&
      solver->row_refactor_input_ptr != NULL &&
      (solver->n == 0u ||
       solver->row_refactor_input_needs_cleanup != NULL) &&
      solver->row_refactor_successor_ptr != NULL &&
      (solver->row_refactor_successor_ptr[solver->n] == 0u ||
       solver->row_refactor_successor_rows != NULL) &&
      solver->row_refactor_group_dep_ptr != NULL &&
      (solver->row_refactor_group_dep_ptr[solver->row_refactor_group_count] == 0u ||
       solver->row_refactor_group_dep_rows != NULL) &&
      solver->row_refactor_group_successor_ptr != NULL &&
      (solver->row_refactor_group_dependency_edges == 0u ||
       solver->row_refactor_group_successor_groups != NULL) &&
      solver->row_refactor_group_pred_count != NULL &&
      (solver->row_refactor_group_root_count == 0u ||
       solver->row_refactor_group_roots != NULL) &&
      (solver->row_refactor_segment_input_target_rows == 0u ||
       solver->row_refactor_segment_input_row_ready != NULL) &&
      (solver->row_refactor_input_ptr[solver->n] == 0u ||
       solver->row_refactor_segment_input_target_rows == 0u ||
       (solver->row_refactor_segment_input_target_kind != NULL &&
        solver->row_refactor_segment_input_target_pos != NULL)) &&
      (solver->row_refactor_segment_input_cleanup_entries == 0u ||
       (solver->row_refactor_segment_input_cleanup_ptr != NULL &&
        solver->row_refactor_segment_input_cleanup_cols != NULL))) {
    if (!solver->row_snode_census_done &&
        getenv("KLS_TRACE_ROW_SNODE") != NULL) {
      solver->row_snode_census_done = 1;
      kls_row_snode_census(solver);
    }
    return 1;
  }
  if (lean_only &&
      ((kls_direct_numeric_lean_pattern_capable(solver))) &&
      kls_build_lean_row_refactor_pattern_parallel(solver, 1, 1)) {
    if (trace_prep) {
      const double trace_end = kls_now_seconds();
      fprintf(stderr,
              "KLS lean prep: direct numeric pattern=%.3f total=%.3f ms\n",
              1e3 * (trace_end - trace_begin),
              1e3 * (trace_end - trace_begin));
    }
    return 1;
  }
  if (!lean_only && solver->options.threads > 1 &&
      ((getenv("KLS_ENABLE_PARALLEL_ROW_PATTERN") != NULL &&
        getenv("KLS_ENABLE_DIRECT_NUMERIC_ROW_PATTERN") != NULL) ||
       (kls_direct_numeric_lean_pattern_capable(solver))) &&
      kls_build_lean_row_refactor_pattern_parallel(solver, 1, 0)) {
    /* For a single BTF block, the numeric already contains every object the
       row transpose needs: Pinv maps input rows and LUbx supplies stable L/U
       value addresses.  Building the generic permuted CSC map and a second
       LU pointer cache first is redundant one-time work on large factors. */
    if (trace_prep) {
      const double trace_end = kls_now_seconds();
      fprintf(stderr,
              "KLS full row prep: direct-numeric-pattern=%.3f"
              " total=%.3f ms\n",
              1e3 * (trace_end - trace_begin),
              1e3 * (trace_end - trace_begin));
    }
    return 1;
  }
  if (!lean_only && solver->row_refactor_auto_lower_bound_rejected) {
    /* The parallel builder can decline a completed minimal pattern because
       its proven task floor cannot challenge EGraph.  That is an economic
       verdict, not a construction failure: do not rebuild the same streams
       serially only to reach the identical verdict. */
    return 0;
  }
  if (!kls_build_refactor_map(solver)) {
    return 0;
  }
  const double trace_map = trace_prep ? kls_now_seconds() : 0.0;
  if (!kls_build_refactor_lu_pointer_cache(solver)) {
    return 0;
  }
  const double trace_lu = trace_prep ? kls_now_seconds() : 0.0;
  if (!kls_row_refactor_pattern_is_eligible(solver)) {
    return 0;
  }
  const double trace_eligible = trace_prep ? kls_now_seconds() : 0.0;

  if (!lean_only && solver->options.threads > 1 &&
      getenv("KLS_ENABLE_PARALLEL_ROW_PATTERN") != NULL &&
      kls_build_lean_row_refactor_pattern_parallel(solver, 0, 0)) {
    if (trace_prep) {
      const double trace_end = kls_now_seconds();
      fprintf(stderr,
              "KLS full row prep: map=%.3f lu=%.3f eligible=%.3f"
              " parallel-pattern=%.3f total=%.3f ms\n",
              1e3 * (trace_map - trace_begin),
              1e3 * (trace_lu - trace_map),
              1e3 * (trace_eligible - trace_lu),
              1e3 * (trace_end - trace_eligible),
              1e3 * (trace_end - trace_begin));
    }
    return 1;
  }
  if (!lean_only && solver->row_refactor_auto_lower_bound_rejected) {
    return 0;
  }

  const UF_long factor_entries = solver->numeric->lnz + solver->numeric->unz;
  const int generic_parallel_lean_pattern =
    solver->options.threads > 1 &&
    solver->n >= 4096u && factor_entries >= 200000u;
  if (lean_only &&
      (generic_parallel_lean_pattern ||
       (getenv("KLS_DISABLE_GENERIC_PARALLEL_LEAN_PATTERN") == NULL &&
        getenv("KLS_ENABLE_GENERIC_PARALLEL_LEAN_PATTERN") != NULL)) &&
      kls_build_lean_row_refactor_pattern_parallel(solver, 0, 1)) {
    if (trace_prep) {
      const double trace_end = kls_now_seconds();
      fprintf(stderr,
              "KLS lean prep: map=%.3f lu=%.3f eligible=%.3f"
              " pattern=%.3f total=%.3f ms\n",
              1e3 * (trace_map - trace_begin),
              1e3 * (trace_lu - trace_map),
              1e3 * (trace_eligible - trace_lu),
              1e3 * (trace_end - trace_eligible),
              1e3 * (trace_end - trace_begin));
    }
    return 1;
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

  for (UF_long block = 0; block < solver->symbolic->nblocks; ++block) {
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    if (k1 > k2 || k2 > n) {
      free(l_ptr);
      free(u_ptr);
      free(input_ptr);
      free(next);
      return 0;
    }
    const UF_long nk = k2 - k1;
    for (UF_long local_j = 0; local_j < nk; ++local_j) {
      const UF_long j = k1 + local_j;
      UF_long *li = solver->refactor_l_indices[j];
      const UF_long len = kls_numeric_l_len_in_block(solver, k1, k2, j);
      if (len > 0u && li == NULL) {
        free(l_ptr);
        free(u_ptr);
        free(input_ptr);
        free(next);
        return 0;
      }
      for (UF_long p = 0; p < len; ++p) {
        const UF_long local_row = li[p];
        const UF_long row = k1 + local_row;
        if (local_row <= local_j || local_row >= nk || row >= n) {
          free(l_ptr);
          free(u_ptr);
          free(input_ptr);
          free(next);
          return 0;
        }
        l_ptr[row + 1u]++;
      }
    }
  }

  for (UF_long block = 0; block < solver->symbolic->nblocks; ++block) {
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    const UF_long nk = k2 - k1;
    for (UF_long local_k = 0; local_k < nk; ++local_k) {
      const UF_long k = k1 + local_k;
      UF_long *ui = solver->refactor_u_indices[k];
      const UF_long len = kls_numeric_u_len_in_block(solver, k1, k2, k);
      if (len > 0u && ui == NULL) {
        free(l_ptr);
        free(u_ptr);
        free(input_ptr);
        free(next);
        return 0;
      }
      for (UF_long p = 0; p < len; ++p) {
        const UF_long local_row = ui[p];
        const UF_long row = k1 + local_row;
        if (local_row >= local_k || local_row >= nk || row >= n) {
          free(l_ptr);
          free(u_ptr);
          free(input_ptr);
          free(next);
          return 0;
        }
        u_ptr[row + 1u]++;
      }
    }
  }

  for (UF_long k = 0; k < n; ++k) {
    UF_long input_begin = solver->refactor_col_ptr[k];
    if (solver->symbolic->nblocks > 1u) {
      if (solver->refactor_block_start == NULL ||
          solver->refactor_block_start[k] < input_begin ||
          solver->refactor_block_start[k] > solver->refactor_col_ptr[k + 1u]) {
        free(l_ptr);
        free(u_ptr);
        free(input_ptr);
        free(next);
        return 0;
      }
      input_begin = solver->refactor_block_start[k];
    }
    for (UF_long p = input_begin; p < solver->refactor_col_ptr[k + 1u]; ++p) {
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
  double *u_row_values = kls_aligned_double_values(u_nnz);
  UF_long *input_cols = input_nnz > 0u
    ? (UF_long *)malloc((size_t)input_nnz * sizeof(*input_cols)) : NULL;
  UF_long *input_pos = input_nnz > 0u
    ? (UF_long *)malloc((size_t)input_nnz * sizeof(*input_pos)) : NULL;
  uint32_t *packed_input = NULL;
  if (lean_only && input_nnz > 0u &&
      kls_lean_packed_input_is_profitable(solver, input_nnz) &&
      n <= (UF_long)UINT16_MAX + 1u &&
      solver->nnz <= (UF_long)UINT16_MAX + 1u) {
    packed_input =
      (uint32_t *)malloc((size_t)input_nnz * sizeof(*packed_input));
  }
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
    free(packed_input);
    return 0;
  }

  memcpy(next, l_ptr, (size_t)n * sizeof(*next));
  for (UF_long block = 0; block < solver->symbolic->nblocks; ++block) {
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    const UF_long nk = k2 - k1;
    for (UF_long local_j = 0; local_j < nk; ++local_j) {
      const UF_long j = k1 + local_j;
      UF_long *li = solver->refactor_l_indices[j];
      double *lx = solver->refactor_l_values[j];
      const UF_long len = kls_numeric_l_len_in_block(solver, k1, k2, j);
      for (UF_long p = 0; p < len; ++p) {
        const UF_long row = k1 + li[p];
        const UF_long dst = next[row]++;
        l_cols[dst] = j;
        l_values[dst] = lx + p;
      }
    }
  }
  memcpy(next, u_ptr, (size_t)n * sizeof(*next));
  for (UF_long block = 0; block < solver->symbolic->nblocks; ++block) {
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    const UF_long nk = k2 - k1;
    for (UF_long local_k = 0; local_k < nk; ++local_k) {
      const UF_long k = k1 + local_k;
      UF_long *ui = solver->refactor_u_indices[k];
      double *ux = solver->refactor_u_values[k];
      const UF_long len = kls_numeric_u_len_in_block(solver, k1, k2, k);
      for (UF_long p = 0; p < len; ++p) {
        const UF_long row = k1 + ui[p];
        const UF_long dst = next[row]++;
        u_cols[dst] = k;
        u_values[dst] = ux + p;
      }
    }
  }
  memcpy(next, input_ptr, (size_t)n * sizeof(*next));
  for (UF_long k = 0; k < n; ++k) {
    UF_long input_begin = solver->refactor_col_ptr[k];
    if (solver->symbolic->nblocks > 1u) {
      input_begin = solver->refactor_block_start[k];
    }
    for (UF_long p = input_begin; p < solver->refactor_col_ptr[k + 1u]; ++p) {
      const UF_long row = solver->refactor_row_idx[p];
      const UF_long dst = next[row]++;
      input_cols[dst] = k;
      input_pos[dst] = solver->refactor_input_pos[p];
      if (packed_input != NULL) {
        packed_input[dst] =
          ((uint32_t)solver->refactor_input_pos[p] << 16u) | (uint32_t)k;
      }
    }
  }
  if (lean_only) {
    return kls_finish_lean_row_refactor_pattern_from_arrays(
      solver, l_ptr, l_cols, l_values, l_row_values,
      u_ptr, u_cols, u_values, u_row_values,
      input_ptr, input_cols, input_pos, packed_input, next);
  }
  return kls_finish_row_refactor_pattern_from_arrays(
    solver, l_ptr, l_cols, l_values, l_row_values,
    u_ptr, u_cols, u_values, u_row_values,
    input_ptr, input_cols, input_pos, next);
}

typedef struct kls_row_refactor_group_work_entry {
  UF_long group;
  double work;
} kls_row_refactor_group_work_entry;

static int kls_compare_row_refactor_group_work_desc(const void *a,
                                                    const void *b) {
  const kls_row_refactor_group_work_entry *ga =
    (const kls_row_refactor_group_work_entry *)a;
  const kls_row_refactor_group_work_entry *gb =
    (const kls_row_refactor_group_work_entry *)b;
  if (ga->work > gb->work) {
    return -1;
  }
  if (ga->work < gb->work) {
    return 1;
  }
  if (ga->group < gb->group) {
    return -1;
  }
  if (ga->group > gb->group) {
    return 1;
  }
  return 0;
}

static KLS_ALWAYS_INLINE int kls_row_refactor_group_work_before(
  UF_long a_group,
  double a_work,
  UF_long b_group,
  double b_work) {
  if (a_work > b_work) {
    return 1;
  }
  if (a_work < b_work) {
    return 0;
  }
  return a_group < b_group;
}

static double kls_row_refactor_compute_group_work(const kls_solver *solver,
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
  const int dense_group =
    solver->row_refactor_group_dense != NULL &&
    solver->row_refactor_group_dense[group];

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
        if (dense_group && dep >= row_begin && dep < row) {
          continue;
        }
        work += 1.0 + (double)(solver->row_refactor_u_ptr[dep + 1u] -
                               solver->row_refactor_u_ptr[dep]);
      }
    }
  }

  if (dense_group) {
    const double width = (double)(row_end - row_begin);
    const double trailing = solver->row_refactor_group_trailing_len != NULL
      ? (double)solver->row_refactor_group_trailing_len[group] : 0.0;
    work += width * width * width / 3.0;
    work += width * width * trailing * 0.5;
  }
  return work;
}

static double kls_row_refactor_group_work(const kls_solver *solver,
                                          UF_long group) {
  if (solver != NULL &&
      solver->row_refactor_group_work != NULL &&
      group < solver->row_refactor_group_count) {
    return solver->row_refactor_group_work[group];
  }
  return kls_row_refactor_compute_group_work(solver, group);
}

static KLS_ALWAYS_INLINE int
kls_row_refactor_group_separator_component_fast(
  const kls_solver *solver,
  UF_long group,
  unsigned int *component_out) {
  if (component_out != NULL) {
    *component_out = UINT_MAX;
  }
  if (solver == NULL || component_out == NULL ||
      !kls_separator_analysis_has_global_range(&solver->separator) ||
      solver->separator.order_component == NULL ||
      solver->row_refactor_group_ptr == NULL ||
      group >= solver->row_refactor_group_count) {
    return 0;
  }
  const UF_long row = solver->row_refactor_group_ptr[group];
  if (row < solver->separator.global_begin ||
      row >= solver->separator.global_end) {
    return 0;
  }
  const UF_long local = row - solver->separator.global_begin;
  if (local >= solver->separator.n) {
    return 0;
  }
  const unsigned int component = solver->separator.order_component[local];
  if ((UF_long)component >= solver->separator.component_count) {
    return 0;
  }
  *component_out = component;
  return 1;
}

/* The separator partition needs a runtime proxy, not a FLOP count.  Dense
   row groups execute their cubic work through compact blocked updates, so
   raw symbolic work substantially overprices the largest components.  Keep
   the established model as mode zero and expose compressed probes while the
   paper-union calibration is being measured. */
static double kls_row_refactor_separator_schedule_work(
  const kls_solver *solver,
  UF_long group,
  int model) {
  const double work = kls_row_refactor_group_work(solver, group);
  if (model == 0 && solver != NULL &&
      solver->row_refactor_separator_component_calibrated &&
      solver->row_refactor_separator_component_scale != NULL &&
      solver->row_refactor_separator_component_scale_count ==
        solver->separator.component_count) {
    unsigned int component = UINT_MAX;
    if (kls_row_refactor_group_separator_component_fast(
          solver, group, &component)) {
      return work *
        solver->row_refactor_separator_component_scale[component];
    }
  }
  if (model == 1 && solver != NULL &&
      solver->row_refactor_group_ptr != NULL &&
      group < solver->row_refactor_group_count) {
    return 1.0 +
      (double)(solver->row_refactor_group_ptr[group + 1u] -
               solver->row_refactor_group_ptr[group]);
  }
  if (model == 2) {
    return sqrt(work > 1.0 ? work : 1.0);
  }
  if (model == 3) {
    return cbrt(work > 1.0 ? work : 1.0);
  }
  if (model == 4) {
    static double group_overhead = -1.0;
    if (group_overhead < 0.0) {
      group_overhead = 2000.0;
      const char *env = getenv("KLS_ROW_REFACTOR_GROUP_OVERHEAD");
      if (env != NULL && env[0] != '\0') {
        char *end = NULL;
        const double candidate = strtod(env, &end);
        if (end != env && *end == '\0' && isfinite(candidate) &&
            candidate >= 0.0 && candidate <= 1.0e9) {
          group_overhead = candidate;
        }
      }
    }
    return work + group_overhead;
  }
  return work;
}

static int kls_sort_row_refactor_group_list_by_work(
  const kls_solver *solver,
  UF_long *groups,
  UF_long count) {
  if (solver == NULL || groups == NULL || count <= 1u) {
    return 1;
  }
  if (count <= KLS_ROW_REFACTOR_SMALL_SORT_MAX) {
    for (UF_long i = 1u; i < count; ++i) {
      const UF_long group = groups[i];
      const double work = kls_row_refactor_group_work(solver, group);
      UF_long j = i;
      while (j > 0u) {
        const UF_long prev_group = groups[j - 1u];
        const double prev_work =
          kls_row_refactor_group_work(solver, prev_group);
        if (!kls_row_refactor_group_work_before(group, work,
                                                prev_group, prev_work)) {
          break;
        }
        groups[j] = prev_group;
        j--;
      }
      groups[j] = group;
    }
    return 1;
  }
  if (count > (UF_long)(SIZE_MAX / sizeof(kls_row_refactor_group_work_entry))) {
    return 0;
  }
  kls_row_refactor_group_work_entry *entries =
    (kls_row_refactor_group_work_entry *)
      malloc((size_t)count * sizeof(*entries));
  if (entries == NULL) {
    return 0;
  }
  for (UF_long i = 0; i < count; ++i) {
    entries[i].group = groups[i];
    entries[i].work = kls_row_refactor_group_work(solver, groups[i]);
  }
  qsort(entries, (size_t)count, sizeof(*entries),
        kls_compare_row_refactor_group_work_desc);
  for (UF_long i = 0; i < count; ++i) {
    groups[i] = entries[i].group;
  }
  free(entries);
  return 1;
}

static int kls_sort_row_refactor_group_list_by_cached_work(
  UF_long *groups,
  UF_long count,
  const double *group_work) {
  if (groups == NULL || group_work == NULL || count <= 1u) {
    return 1;
  }
  if (count <= KLS_ROW_REFACTOR_SMALL_SORT_MAX) {
    for (UF_long i = 1u; i < count; ++i) {
      const UF_long group = groups[i];
      const double work = group_work[group];
      UF_long j = i;
      while (j > 0u) {
        const UF_long prev_group = groups[j - 1u];
        const double prev_work = group_work[prev_group];
        if (!kls_row_refactor_group_work_before(group, work,
                                                prev_group, prev_work)) {
          break;
        }
        groups[j] = prev_group;
        j--;
      }
      groups[j] = group;
    }
    return 1;
  }
  if (count > (UF_long)(SIZE_MAX / sizeof(kls_row_refactor_group_work_entry))) {
    return 0;
  }
  kls_row_refactor_group_work_entry *entries =
    (kls_row_refactor_group_work_entry *)
      malloc((size_t)count * sizeof(*entries));
  if (entries == NULL) {
    return 0;
  }
  for (UF_long i = 0; i < count; ++i) {
    entries[i].group = groups[i];
    entries[i].work = group_work[groups[i]];
  }
  qsort(entries, (size_t)count, sizeof(*entries),
        kls_compare_row_refactor_group_work_desc);
  for (UF_long i = 0; i < count; ++i) {
    groups[i] = entries[i].group;
  }
  free(entries);
  return 1;
}

static int kls_build_row_refactor_private_thread_queues(
  const kls_solver *solver,
  const UF_long *source_groups,
  UF_long group_count,
  int thread_count,
  UF_long **private_groups_out,
  UF_long **thread_ptr_out) {
  KLS_SET_OPTIONAL_OUTPUT(private_groups_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(thread_ptr_out, NULL);
  if (private_groups_out == NULL || thread_ptr_out == NULL ||
      solver == NULL || source_groups == NULL || group_count == 0u ||
      thread_count <= 1) {
    return 0;
  }
  const size_t group_count_size = (size_t)group_count;
  const size_t thread_count_size = (size_t)thread_count;
  if ((UF_long)group_count_size != group_count ||
      group_count_size > SIZE_MAX / sizeof(UF_long) ||
      thread_count_size > (SIZE_MAX / sizeof(UF_long)) - 1u ||
      thread_count_size > SIZE_MAX / sizeof(double)) {
    return 0;
  }

  UF_long *private_groups =
    (UF_long *)malloc(group_count_size * sizeof(*private_groups));
  UF_long *thread_ptr =
    (UF_long *)calloc(thread_count_size + 1u, sizeof(*thread_ptr));
  UF_long *thread_counts =
    (UF_long *)calloc(thread_count_size, sizeof(*thread_counts));
  double *thread_work =
    (double *)calloc(thread_count_size, sizeof(*thread_work));
  if (private_groups == NULL || thread_ptr == NULL ||
      thread_counts == NULL || thread_work == NULL) {
    free(private_groups);
    free(thread_ptr);
    free(thread_counts);
    free(thread_work);
    return 0;
  }

  for (UF_long pos = 0; pos < group_count; ++pos) {
    const UF_long group = source_groups[pos];
    if (group >= solver->row_refactor_group_count) {
      free(private_groups);
      free(thread_ptr);
      free(thread_counts);
      free(thread_work);
      return 0;
    }
    int target = 0;
    double target_work = thread_work[0];
    for (int tid = 1; tid < thread_count; ++tid) {
      if (thread_work[tid] < target_work) {
        target = tid;
        target_work = thread_work[tid];
      }
    }
    thread_counts[target]++;
    thread_work[target] += kls_row_refactor_group_work(solver, group);
  }

  for (int tid = 0; tid < thread_count; ++tid) {
    thread_ptr[tid + 1] = thread_ptr[tid] + thread_counts[tid];
    thread_counts[tid] = thread_ptr[tid];
    thread_work[tid] = 0.0;
  }
  if (thread_ptr[thread_count] != group_count) {
    free(private_groups);
    free(thread_ptr);
    free(thread_counts);
    free(thread_work);
    return 0;
  }

  for (UF_long pos = 0; pos < group_count; ++pos) {
    const UF_long group = source_groups[pos];
    int target = 0;
    double target_work = thread_work[0];
    for (int tid = 1; tid < thread_count; ++tid) {
      if (thread_work[tid] < target_work) {
        target = tid;
        target_work = thread_work[tid];
      }
    }
    const UF_long dst = thread_counts[target]++;
    if (dst >= group_count) {
      free(private_groups);
      free(thread_ptr);
      free(thread_counts);
      free(thread_work);
      return 0;
    }
    private_groups[dst] = group;
    thread_work[target] += kls_row_refactor_group_work(solver, group);
  }

  for (int tid = 0; tid < thread_count; ++tid) {
    if (thread_counts[tid] != thread_ptr[tid + 1]) {
      free(private_groups);
      free(thread_ptr);
      free(thread_counts);
      free(thread_work);
      return 0;
    }
  }

  free(thread_counts);
  free(thread_work);
  *private_groups_out = private_groups;
  *thread_ptr_out = thread_ptr;
  return 1;
}

typedef struct kls_separator_component_work_entry {
  unsigned int component;
  unsigned char kind;
  double work;
} kls_separator_component_work_entry;

static int kls_compare_separator_component_work_desc(const void *a,
                                                     const void *b) {
  const kls_separator_component_work_entry *ca =
    (const kls_separator_component_work_entry *)a;
  const kls_separator_component_work_entry *cb =
    (const kls_separator_component_work_entry *)b;
  if (ca->kind < cb->kind) {
    return -1;
  }
  if (ca->kind > cb->kind) {
    return 1;
  }
  if (ca->work > cb->work) {
    return -1;
  }
  if (ca->work < cb->work) {
    return 1;
  }
  if (ca->component < cb->component) {
    return -1;
  }
  if (ca->component > cb->component) {
    return 1;
  }
  return 0;
}

static int kls_row_refactor_group_separator_component(
  const kls_solver *solver,
  UF_long group,
  unsigned int *component_out) {
  if (component_out != NULL) {
    *component_out = UINT_MAX;
  }
  if (solver == NULL || component_out == NULL ||
      !kls_separator_analysis_has_global_range(&solver->separator) ||
      solver->separator.component_count == 0u ||
      solver->separator.component_kind == NULL ||
      solver->separator.order_component == NULL ||
      solver->row_refactor_group_ptr == NULL ||
      group >= solver->row_refactor_group_count) {
    return 0;
  }
  const UF_long row_begin = solver->row_refactor_group_ptr[group];
  const UF_long row_end = solver->row_refactor_group_ptr[group + 1u];
  if (row_begin >= row_end || row_end > solver->n ||
      row_begin < solver->separator.global_begin ||
      row_end > solver->separator.global_end) {
    return 0;
  }
  const UF_long local_begin = row_begin - solver->separator.global_begin;
  const UF_long local_end = row_end - solver->separator.global_begin;
  if (local_end > solver->separator.n) {
    return 0;
  }
  const unsigned int component =
    solver->separator.order_component[local_begin];
  if ((UF_long)component >= solver->separator.component_count) {
    return 0;
  }
  for (UF_long row = local_begin + 1u; row < local_end; ++row) {
    if (solver->separator.order_component[row] != component) {
      return 0;
    }
  }
  *component_out = component;
  return 1;
}

static int kls_build_row_refactor_separator_private_thread_queues(
  const kls_solver *solver,
  const UF_long *source_groups,
  UF_long group_count,
  int thread_count,
  UF_long **private_groups_out,
  UF_long **thread_ptr_out,
  UF_long *component_count_out) {
  KLS_SET_OPTIONAL_OUTPUT(private_groups_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(thread_ptr_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(component_count_out, 0);
  if (private_groups_out == NULL || thread_ptr_out == NULL ||
      component_count_out == NULL || solver == NULL ||
      source_groups == NULL || group_count == 0u || thread_count <= 1 ||
      !kls_separator_analysis_has_global_range(&solver->separator) ||
      solver->separator.component_count <= 1u ||
      solver->separator.component_kind == NULL ||
      solver->separator.order_component == NULL) {
    return 0;
  }

  const UF_long component_count = solver->separator.component_count;
  const size_t group_count_size = (size_t)group_count;
  const size_t component_count_size = (size_t)component_count;
  const size_t thread_count_size = (size_t)thread_count;
  if ((UF_long)group_count_size != group_count ||
      (UF_long)component_count_size != component_count ||
      group_count_size > SIZE_MAX / sizeof(UF_long) ||
      component_count_size > SIZE_MAX / sizeof(UF_long) ||
      component_count_size > SIZE_MAX / sizeof(double) ||
      component_count_size > SIZE_MAX / sizeof(int) ||
      component_count_size >
        SIZE_MAX / sizeof(kls_separator_component_work_entry) ||
      thread_count_size > (SIZE_MAX / sizeof(UF_long)) - 1u ||
      thread_count_size > SIZE_MAX / sizeof(double)) {
    return 0;
  }

  UF_long *component_group_count =
    (UF_long *)calloc(component_count_size, sizeof(*component_group_count));
  double *component_work =
    (double *)calloc(component_count_size, sizeof(*component_work));
  int *component_thread =
    (int *)malloc(component_count_size * sizeof(*component_thread));
  kls_separator_component_work_entry *active_components =
    (kls_separator_component_work_entry *)
      malloc(component_count_size * sizeof(*active_components));
  UF_long *private_groups =
    (UF_long *)malloc(group_count_size * sizeof(*private_groups));
  UF_long *thread_ptr =
    (UF_long *)calloc(thread_count_size + 1u, sizeof(*thread_ptr));
  UF_long *thread_counts =
    (UF_long *)calloc(thread_count_size, sizeof(*thread_counts));
  double *thread_work =
    (double *)calloc(thread_count_size, sizeof(*thread_work));
  if (component_group_count == NULL || component_work == NULL ||
      component_thread == NULL || active_components == NULL ||
      private_groups == NULL || thread_ptr == NULL ||
      thread_counts == NULL || thread_work == NULL) {
    free(component_group_count);
    free(component_work);
    free(component_thread);
    free(active_components);
    free(private_groups);
    free(thread_ptr);
    free(thread_counts);
    free(thread_work);
    return 0;
  }
  for (UF_long component = 0; component < component_count; ++component) {
    component_thread[component] = -1;
  }

  UF_long active_count = 0;
  for (UF_long pos = 0; pos < group_count; ++pos) {
    const UF_long group = source_groups[pos];
    unsigned int component = UINT_MAX;
    if (!kls_row_refactor_group_separator_component(solver, group,
                                                    &component)) {
      free(component_group_count);
      free(component_work);
      free(component_thread);
      free(active_components);
      free(private_groups);
      free(thread_ptr);
      free(thread_counts);
      free(thread_work);
      return 0;
    }
    if (component_group_count[component] == 0u) {
      active_components[active_count].component = component;
      active_components[active_count].kind =
        solver->separator.component_kind[component];
      active_components[active_count].work = 0.0;
      active_count++;
    }
    component_group_count[component]++;
    component_work[component] += kls_row_refactor_group_work(solver, group);
  }
  if (active_count <= 1u) {
    free(component_group_count);
    free(component_work);
    free(component_thread);
    free(active_components);
    free(private_groups);
    free(thread_ptr);
    free(thread_counts);
    free(thread_work);
    return 0;
  }

  for (UF_long pos = 0; pos < active_count; ++pos) {
    active_components[pos].work =
      component_work[active_components[pos].component];
  }
  qsort(active_components, (size_t)active_count,
        sizeof(*active_components),
        kls_compare_separator_component_work_desc);

  for (UF_long pos = 0; pos < active_count; ++pos) {
    const unsigned int component = active_components[pos].component;
    int target = 0;
    double target_work = thread_work[0];
    for (int tid = 1; tid < thread_count; ++tid) {
      if (thread_work[tid] < target_work) {
        target = tid;
        target_work = thread_work[tid];
      }
    }
    component_thread[component] = target;
    thread_counts[target] += component_group_count[component];
    thread_work[target] += component_work[component];
  }

  for (int tid = 0; tid < thread_count; ++tid) {
    thread_ptr[tid + 1] = thread_ptr[tid] + thread_counts[tid];
    thread_counts[tid] = thread_ptr[tid];
  }
  if (thread_ptr[thread_count] != group_count) {
    free(component_group_count);
    free(component_work);
    free(component_thread);
    free(active_components);
    free(private_groups);
    free(thread_ptr);
    free(thread_counts);
    free(thread_work);
    return 0;
  }

  for (UF_long pos = 0; pos < group_count; ++pos) {
    const UF_long group = source_groups[pos];
    unsigned int component = UINT_MAX;
    if (!kls_row_refactor_group_separator_component(solver, group,
                                                    &component) ||
        component_thread[component] < 0) {
      free(component_group_count);
      free(component_work);
      free(component_thread);
      free(active_components);
      free(private_groups);
      free(thread_ptr);
      free(thread_counts);
      free(thread_work);
      return 0;
    }
    const int target = component_thread[component];
    const UF_long dst = thread_counts[target]++;
    if (dst >= group_count) {
      free(component_group_count);
      free(component_work);
      free(component_thread);
      free(active_components);
      free(private_groups);
      free(thread_ptr);
      free(thread_counts);
      free(thread_work);
      return 0;
    }
    private_groups[dst] = group;
  }

  for (int tid = 0; tid < thread_count; ++tid) {
    if (thread_counts[tid] != thread_ptr[tid + 1]) {
      free(component_group_count);
      free(component_work);
      free(component_thread);
      free(active_components);
      free(private_groups);
      free(thread_ptr);
      free(thread_counts);
      free(thread_work);
      return 0;
    }
  }

  free(component_group_count);
  free(component_work);
  free(component_thread);
  free(active_components);
  free(thread_counts);
  free(thread_work);
  *private_groups_out = private_groups;
  *thread_ptr_out = thread_ptr;
  *component_count_out = active_count;
  return 1;
}

static int kls_ensure_row_refactor_ready_queue_workspace(
  kls_solver *solver,
  UF_long group_count);

static int kls_prepare_row_refactor_separator_flop_ready_queue(
  kls_solver *solver,
  int thread_count,
  UF_long **ready_groups_out,
  atomic_uint **ready_slots_out,
  atomic_ulong **remaining_preds_out,
  unsigned char **tail_groups_out,
  UF_long **private_groups_out,
  UF_long **private_thread_ptr_out,
  unsigned char **private_group_mask_out,
  UF_long *tail_count_out,
  UF_long *initial_ready_out,
  UF_long *private_group_count_out,
  UF_long *pipeline_group_count_out,
  UF_long *closure_group_count_out,
  UF_long *component_count_out,
  UF_long *private_threads_out,
  UF_long *private_min_groups_out,
  UF_long *private_max_groups_out,
  double *private_min_work_out,
  double *private_max_work_out) {
  KLS_SET_OPTIONAL_OUTPUT(ready_groups_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(ready_slots_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(remaining_preds_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(tail_groups_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(private_groups_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(private_thread_ptr_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(private_group_mask_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(tail_count_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(initial_ready_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(private_group_count_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(pipeline_group_count_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(closure_group_count_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(component_count_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(private_threads_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(private_min_groups_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(private_max_groups_out, 0);
  if (private_min_work_out != NULL) {
    *private_min_work_out = 0.0;
  }
  if (private_max_work_out != NULL) {
    *private_max_work_out = 0.0;
  }
  if (solver == NULL || thread_count <= 1 ||
      ready_groups_out == NULL || ready_slots_out == NULL ||
      remaining_preds_out == NULL || tail_groups_out == NULL ||
      private_groups_out == NULL || private_thread_ptr_out == NULL ||
      private_group_mask_out == NULL || tail_count_out == NULL ||
      initial_ready_out == NULL || private_group_count_out == NULL ||
      pipeline_group_count_out == NULL || closure_group_count_out == NULL ||
      component_count_out == NULL || private_threads_out == NULL ||
      private_min_groups_out == NULL || private_max_groups_out == NULL ||
      private_min_work_out == NULL || private_max_work_out == NULL ||
      solver->row_refactor_group_count == 0u ||
      solver->row_refactor_group_ptr == NULL ||
      solver->row_refactor_group_pred_count == NULL ||
      solver->row_refactor_group_successor_ptr == NULL ||
      (solver->row_refactor_group_dependency_edges > 0u &&
       solver->row_refactor_group_successor_groups == NULL) ||
      !kls_separator_analysis_has_global_range(&solver->separator) ||
      solver->separator.component_count <= 1u ||
      solver->separator.component_kind == NULL ||
      solver->separator.order_component == NULL ||
      solver->separator.component_left_child == NULL ||
      solver->separator.component_right_child == NULL) {
    return 0;
  }

  const UF_long group_count = solver->row_refactor_group_count;
  const UF_long component_count = solver->separator.component_count;
  int separator_work_model = 0;
  {
    const char *model_env =
      getenv("KLS_ROW_REFACTOR_SEPARATOR_WORK_MODEL");
    if (model_env != NULL) {
      if (strcmp(model_env, "rows") == 0) {
        separator_work_model = 1;
      } else if (strcmp(model_env, "sqrt") == 0) {
        separator_work_model = 2;
      } else if (strcmp(model_env, "cbrt") == 0) {
        separator_work_model = 3;
      } else if (strcmp(model_env, "overhead") == 0) {
        separator_work_model = 4;
      }
    } else if (solver->prestatic_dense_spiked_match ||
               solver->dense_spiked_original_pivot_path) {
      /* Dense-spiked systems have a long separator chain whose
         compact-group setup and dependency traffic dominate the nominal
         flop estimate.  Charge each group its measured fixed scheduling
         surface as well as its rows; this was consistently better on the
         route's representative rajat30 shape.  Keep the environment knob
         authoritative for experiments and for callers that need another
         model. */
      separator_work_model = 4;
    }
  }
  double separator_balance_beta =
    KLS_ROW_REFACTOR_SEPARATOR_BALANCE_BETA;
  const char *separator_balance_beta_env =
    getenv("KLS_ROW_REFACTOR_SEPARATOR_BETA");
  if (separator_balance_beta_env != NULL &&
      separator_balance_beta_env[0] != '\0') {
    char *end = NULL;
    const double candidate =
      strtod(separator_balance_beta_env, &end);
    if (end != separator_balance_beta_env && *end == '\0' &&
        isfinite(candidate) && candidate >= 0.25 && candidate <= 4.0) {
      separator_balance_beta = candidate;
    }
  }
  const char *split_past_leaf_env =
    getenv("KLS_ROW_REFACTOR_SPLIT_PAST_LEAF");
  const int split_past_leaf =
    split_past_leaf_env != NULL && split_past_leaf_env[0] != '\0' &&
    !(split_past_leaf_env[0] == '0' && split_past_leaf_env[1] == '\0');
  const char *promote_overweight_leaf_env =
    getenv("KLS_ROW_REFACTOR_PROMOTE_OVERWEIGHT_LEAF");
  const int promote_overweight_leaf =
    promote_overweight_leaf_env != NULL &&
    promote_overweight_leaf_env[0] != '\0' &&
    !(promote_overweight_leaf_env[0] == '0' &&
      promote_overweight_leaf_env[1] == '\0');
  UF_long max_overweight_leaf_promotions = 1u;
  const char *max_leaf_promotions_env =
    getenv("KLS_ROW_REFACTOR_MAX_LEAF_PROMOTIONS");
  if (max_leaf_promotions_env != NULL &&
      max_leaf_promotions_env[0] != '\0') {
    char *end = NULL;
    const unsigned long candidate =
      strtoul(max_leaf_promotions_env, &end, 10);
    if (end != max_leaf_promotions_env && *end == '\0' &&
        candidate > 0ul && candidate <= (unsigned long)component_count) {
      max_overweight_leaf_promotions = (UF_long)candidate;
    }
  }
  const size_t group_count_size = (size_t)group_count;
  const size_t component_count_size = (size_t)component_count;
  const size_t thread_count_size = (size_t)thread_count;
  if ((UF_long)group_count_size != group_count ||
      (UF_long)component_count_size != component_count ||
      group_count_size > SIZE_MAX / sizeof(unsigned int) ||
      group_count_size > SIZE_MAX / sizeof(UF_long) ||
      group_count_size > SIZE_MAX / sizeof(unsigned char) ||
      component_count_size > SIZE_MAX / sizeof(UF_long) ||
      component_count_size > SIZE_MAX / sizeof(double) ||
      component_count_size > SIZE_MAX / sizeof(unsigned char) ||
      component_count_size > SIZE_MAX / sizeof(int) ||
      component_count_size >
        SIZE_MAX / sizeof(kls_separator_component_work_entry) ||
      thread_count_size > (SIZE_MAX / sizeof(UF_long)) - 1u ||
      thread_count_size > SIZE_MAX / sizeof(double)) {
    return 0;
  }

  unsigned int *group_component =
    (unsigned int *)malloc(group_count_size * sizeof(*group_component));
  unsigned char *group_forced_pipeline =
    (unsigned char *)calloc(group_count_size, sizeof(*group_forced_pipeline));
  UF_long *component_group_count =
    (UF_long *)calloc(component_count_size, sizeof(*component_group_count));
  double *component_work =
    (double *)calloc(component_count_size, sizeof(*component_work));
  UF_long *subtree_begin =
    (UF_long *)malloc(component_count_size * sizeof(*subtree_begin));
  UF_long *subtree_end =
    (UF_long *)malloc(component_count_size * sizeof(*subtree_end));
  double *subtree_work =
    (double *)calloc(component_count_size, sizeof(*subtree_work));
  UF_long *candidates =
    (UF_long *)malloc(component_count_size * sizeof(*candidates));
  unsigned char *component_pipeline =
    (unsigned char *)calloc(component_count_size, sizeof(*component_pipeline));
  int *component_thread =
    (int *)malloc(component_count_size * sizeof(*component_thread));
  kls_separator_component_work_entry *candidate_entries =
    (kls_separator_component_work_entry *)
      malloc(component_count_size * sizeof(*candidate_entries));
  UF_long *private_groups =
    (UF_long *)malloc(group_count_size * sizeof(*private_groups));
  UF_long *thread_ptr =
    (UF_long *)calloc(thread_count_size + 1u, sizeof(*thread_ptr));
  UF_long *thread_counts =
    (UF_long *)calloc(thread_count_size, sizeof(*thread_counts));
  double *thread_work =
    (double *)calloc(thread_count_size, sizeof(*thread_work));
  unsigned char *private_group_mask =
    (unsigned char *)calloc(group_count_size, sizeof(*private_group_mask));
#define KLS_FREE_SEPARATOR_QUEUE_TEMP() do { \
    free(group_component);                       \
    free(group_forced_pipeline);                 \
    free(component_group_count);                 \
    free(component_work);                        \
    free(subtree_begin);                         \
    free(subtree_end);                           \
    free(subtree_work);                          \
    free(candidates);                            \
    free(component_pipeline);                    \
    free(component_thread);                      \
    free(candidate_entries);                     \
    free(private_groups);                        \
    free(thread_ptr);                            \
    free(thread_counts);                         \
    free(thread_work);                           \
    free(private_group_mask);                    \
  } while (0)
  if (group_component == NULL || group_forced_pipeline == NULL ||
      component_group_count == NULL ||
      component_work == NULL || subtree_begin == NULL ||
      subtree_end == NULL || subtree_work == NULL ||
      candidates == NULL || component_pipeline == NULL ||
      component_thread == NULL || candidate_entries == NULL ||
      private_groups == NULL || thread_ptr == NULL ||
      thread_counts == NULL || thread_work == NULL ||
      private_group_mask == NULL) {
    KLS_FREE_SEPARATOR_QUEUE_TEMP();
    return 0;
  }
  for (UF_long component = 0; component < component_count; ++component) {
    component_thread[component] = -1;
  }

  for (UF_long group = 0; group < group_count; ++group) {
    unsigned int component = UINT_MAX;
    if (!kls_row_refactor_group_separator_component(solver, group,
                                                    &component)) {
      const UF_long row_begin = solver->row_refactor_group_ptr[group];
      const UF_long row_end = solver->row_refactor_group_ptr[group + 1u];
      if (row_begin >= row_end || row_end > solver->n ||
          row_begin < solver->separator.global_begin ||
          row_end > solver->separator.global_end) {
        KLS_FREE_SEPARATOR_QUEUE_TEMP();
        return 0;
      }
      const UF_long local_begin =
        row_begin - solver->separator.global_begin;
      if (local_begin >= solver->separator.n) {
        KLS_FREE_SEPARATOR_QUEUE_TEMP();
        return 0;
      }
      component = solver->separator.order_component[local_begin];
      if ((UF_long)component >= component_count) {
        KLS_FREE_SEPARATOR_QUEUE_TEMP();
        return 0;
      }
      group_forced_pipeline[group] = 1u;
    }
    group_component[group] = component;
    component_group_count[component]++;
    component_work[component] +=
      kls_row_refactor_separator_schedule_work(
        solver, group, separator_work_model);
  }

  for (UF_long component = 0; component < component_count; ++component) {
    UF_long begin = component;
    UF_long end = component + 1u;
    double work = component_work[component];
    const UF_long left = solver->separator.component_left_child[component];
    const UF_long right = solver->separator.component_right_child[component];
    if (left != KLS_KLU_EMPTY) {
      if (left >= component || left >= component_count) {
        KLS_FREE_SEPARATOR_QUEUE_TEMP();
        return 0;
      }
      begin = subtree_begin[left] < begin ? subtree_begin[left] : begin;
      end = subtree_end[left] > end ? subtree_end[left] : end;
      work += subtree_work[left];
    }
    if (right != KLS_KLU_EMPTY) {
      if (right >= component || right >= component_count) {
        KLS_FREE_SEPARATOR_QUEUE_TEMP();
        return 0;
      }
      begin = subtree_begin[right] < begin ? subtree_begin[right] : begin;
      end = subtree_end[right] > end ? subtree_end[right] : end;
      work += subtree_work[right];
    }
    subtree_begin[component] = begin;
    subtree_end[component] = end;
    subtree_work[component] = work;
  }

  const UF_long root_component = component_count - 1u;
  if (subtree_work[root_component] <= 0.0) {
    KLS_FREE_SEPARATOR_QUEUE_TEMP();
    return 0;
  }

  UF_long candidate_count = 1;
  candidates[0] = root_component;
  double candidate_work_sum = subtree_work[root_component];
  UF_long pipeline_component_count = 0;
  UF_long promoted_overweight_leaf_count = 0;
  while (candidate_count > 0u && candidate_work_sum > 0.0) {
    UF_long max_pos = KLS_KLU_EMPTY;
    double max_work = -1.0;
    for (UF_long pos = 0u; pos < candidate_count; ++pos) {
      if (split_past_leaf) {
        const UF_long component = candidates[pos];
        const UF_long left =
          solver->separator.component_left_child[component];
        const UF_long right =
          solver->separator.component_right_child[component];
        const int split_left =
          left != KLS_KLU_EMPTY && left < component_count &&
          subtree_work[left] > 0.0;
        const int split_right =
          right != KLS_KLU_EMPTY && right < component_count &&
          subtree_work[right] > 0.0;
        if (!split_left && !split_right) {
          continue;
        }
      }
      const double work = subtree_work[candidates[pos]];
      if (max_pos == KLS_KLU_EMPTY || work > max_work) {
        max_pos = pos;
        max_work = work;
      }
    }
    if (max_pos == KLS_KLU_EMPTY) {
      break;
    }
    if ((double)thread_count * max_work <
        separator_balance_beta * candidate_work_sum) {
      break;
    }

    const UF_long component = candidates[max_pos];
    const UF_long left = solver->separator.component_left_child[component];
    const UF_long right = solver->separator.component_right_child[component];
    const int split_left =
      left != KLS_KLU_EMPTY && left < component_count &&
      subtree_work[left] > 0.0;
    const int split_right =
      right != KLS_KLU_EMPTY && right < component_count &&
      subtree_work[right] > 0.0;
    if (!split_left && !split_right) {
      if (promote_overweight_leaf &&
          promoted_overweight_leaf_count < max_overweight_leaf_promotions) {
        candidates[max_pos] = candidates[candidate_count - 1u];
        candidate_count--;
        candidate_work_sum -= subtree_work[component];
        component_pipeline[component] = 1u;
        pipeline_component_count++;
        promoted_overweight_leaf_count++;
        continue;
      }
      break;
    }
    candidates[max_pos] = candidates[candidate_count - 1u];
    candidate_count--;
    candidate_work_sum -= subtree_work[component];
    component_pipeline[component] = 1u;
    pipeline_component_count++;

    if (split_left) {
      candidates[candidate_count++] = left;
      candidate_work_sum += subtree_work[left];
    }
    if (split_right) {
      candidates[candidate_count++] = right;
      candidate_work_sum += subtree_work[right];
    }
  }
  if (candidate_count == 0u || pipeline_component_count == 0u) {
    KLS_FREE_SEPARATOR_QUEUE_TEMP();
    return 0;
  }

  UF_long assigned_component_count = 0;
  for (UF_long component = 0; component < component_count; ++component) {
    if (component_group_count[component] == 0u ||
        component_pipeline[component]) {
      continue;
    }
    if (component > (UF_long)UINT_MAX) {
      KLS_FREE_SEPARATOR_QUEUE_TEMP();
      return 0;
    }
    candidate_entries[assigned_component_count].component =
      (unsigned int)component;
    candidate_entries[assigned_component_count].kind = 0u;
    candidate_entries[assigned_component_count].work =
      component_work[component];
    assigned_component_count++;
  }
  if (assigned_component_count == 0u) {
    KLS_FREE_SEPARATOR_QUEUE_TEMP();
    return 0;
  }
  qsort(candidate_entries, (size_t)assigned_component_count,
        sizeof(*candidate_entries),
        kls_compare_separator_component_work_desc);

  for (UF_long pos = 0; pos < assigned_component_count; ++pos) {
    const unsigned int component = candidate_entries[pos].component;
    int target = 0;
    double target_work = thread_work[0];
    for (int tid = 1; tid < thread_count; ++tid) {
      if (thread_work[tid] < target_work) {
        target = tid;
        target_work = thread_work[tid];
      }
    }
    component_thread[component] = target;
    thread_work[target] += component_work[component];
  }

  UF_long pipeline_group_count = 0;
  for (UF_long group = 0; group < group_count; ++group) {
    const unsigned int component = group_component[group];
    if ((UF_long)component >= component_count) {
      KLS_FREE_SEPARATOR_QUEUE_TEMP();
      return 0;
    }
    if (component_pipeline[component] || group_forced_pipeline[group]) {
      pipeline_group_count++;
      continue;
    }
    const int target = component_thread[component];
    if (target < 0 || target >= thread_count) {
      KLS_FREE_SEPARATOR_QUEUE_TEMP();
      return 0;
    }
    private_group_mask[group] = 1u;
  }

  int promoted_pipeline_group = 1;
  UF_long closure_group_count = 0;
  while (promoted_pipeline_group) {
    promoted_pipeline_group = 0;
    for (UF_long group = 0; group < group_count; ++group) {
      if (private_group_mask[group]) {
        continue;
      }
      const UF_long begin = solver->row_refactor_group_successor_ptr[group];
      const UF_long end =
        solver->row_refactor_group_successor_ptr[group + 1u];
      if (begin > end ||
          end > solver->row_refactor_group_dependency_edges) {
        KLS_FREE_SEPARATOR_QUEUE_TEMP();
        return 0;
      }
      for (UF_long pos = begin; pos < end; ++pos) {
        const UF_long successor =
          solver->row_refactor_group_successor_groups[pos];
        if (successor >= group_count) {
          KLS_FREE_SEPARATOR_QUEUE_TEMP();
          return 0;
        }
        if (private_group_mask[successor]) {
          private_group_mask[successor] = 0u;
          pipeline_group_count++;
          closure_group_count++;
          promoted_pipeline_group = 1;
        }
      }
    }
  }

  /* The initial component assignment uses the work of every group in each
     component.  Pipeline closure above can subsequently remove most (or all)
     of a component's groups, leaving the retained private queues badly
     imbalanced even though the pre-closure estimate was balanced.  Recompute
     component weights from the groups that actually remain private and run
     the same largest-first assignment again. */
  {
    const char *rebalance_env =
      getenv("KLS_ROW_REFACTOR_REBALANCE_PRIVATE");
    const int rebalance_private =
      rebalance_env != NULL
        ? (rebalance_env[0] != '\0' &&
           !(rebalance_env[0] == '0' && rebalance_env[1] == '\0'))
        : solver->dense_spiked_original_pivot_path;
    if (rebalance_private) {
      memset(component_group_count, 0,
             component_count_size * sizeof(*component_group_count));
      memset(component_work, 0,
             component_count_size * sizeof(*component_work));
      memset(thread_work, 0, thread_count_size * sizeof(*thread_work));
      for (UF_long component = 0; component < component_count; ++component) {
        component_thread[component] = -1;
      }

      UF_long retained_component_count = 0;
      for (UF_long group = 0; group < group_count; ++group) {
        if (!private_group_mask[group]) {
          continue;
        }
        const unsigned int component = group_component[group];
        if ((UF_long)component >= component_count) {
          KLS_FREE_SEPARATOR_QUEUE_TEMP();
          return 0;
        }
        if (component_group_count[component] == 0u) {
          candidate_entries[retained_component_count].component = component;
          candidate_entries[retained_component_count].kind = 0u;
          candidate_entries[retained_component_count].work = 0.0;
          retained_component_count++;
        }
        component_group_count[component]++;
        component_work[component] +=
          kls_row_refactor_separator_schedule_work(
            solver, group, separator_work_model);
      }
      for (UF_long pos = 0; pos < retained_component_count; ++pos) {
        candidate_entries[pos].work =
          component_work[candidate_entries[pos].component];
      }
      qsort(candidate_entries, (size_t)retained_component_count,
            sizeof(*candidate_entries),
            kls_compare_separator_component_work_desc);
      for (UF_long pos = 0; pos < retained_component_count; ++pos) {
        const unsigned int component = candidate_entries[pos].component;
        int target = 0;
        double target_work = thread_work[0];
        for (int tid = 1; tid < thread_count; ++tid) {
          if (thread_work[tid] < target_work) {
            target = tid;
            target_work = thread_work[tid];
          }
        }
        component_thread[component] = target;
        thread_work[target] += component_work[component];
      }
      assigned_component_count = retained_component_count;
    }
  }

  /* A separator component is a good ownership seed, but numeric L edges can
     still cross between two components assigned to different workers.  The
     dependency-ready implementation tolerates those edges by spinning.  For
     the static-private variant, promote each cross-thread target (and, in
     the same topological pass, all of its descendants) to the pipeline.  Row
     groups are numbered in factor order, so every successor is visited after
     its predecessor and one linear pass is sufficient. */
  const char *promote_cross_thread_env =
    getenv("KLS_ROW_REFACTOR_PROMOTE_CROSS_THREAD_PRIVATE");
  const int promote_cross_thread_private =
    promote_cross_thread_env != NULL &&
    promote_cross_thread_env[0] != '\0' &&
    !(promote_cross_thread_env[0] == '0' &&
      promote_cross_thread_env[1] == '\0');
  if (promote_cross_thread_private) {
    for (UF_long group = 0; group < group_count; ++group) {
      const int source_private = private_group_mask[group] ? 1 : 0;
      const int source_thread = source_private
        ? component_thread[group_component[group]] : -1;
      const UF_long begin = solver->row_refactor_group_successor_ptr[group];
      const UF_long end =
        solver->row_refactor_group_successor_ptr[group + 1u];
      if (begin > end ||
          end > solver->row_refactor_group_dependency_edges) {
        KLS_FREE_SEPARATOR_QUEUE_TEMP();
        return 0;
      }
      for (UF_long pos = begin; pos < end; ++pos) {
        const UF_long successor =
          solver->row_refactor_group_successor_groups[pos];
        if (successor <= group || successor >= group_count) {
          KLS_FREE_SEPARATOR_QUEUE_TEMP();
          return 0;
        }
        if (!private_group_mask[successor]) {
          continue;
        }
        const int successor_thread =
          component_thread[group_component[successor]];
        if (!source_private || source_thread != successor_thread) {
          private_group_mask[successor] = 0u;
          pipeline_group_count++;
          closure_group_count++;
        }
      }
    }
  }

  UF_long private_group_count = 0;
  for (UF_long group = 0; group < group_count; ++group) {
    if (!private_group_mask[group]) {
      continue;
    }
    const unsigned int component = group_component[group];
    if ((UF_long)component >= component_count) {
      KLS_FREE_SEPARATOR_QUEUE_TEMP();
      return 0;
    }
    const int target = component_thread[component];
    if (target < 0 || target >= thread_count) {
      KLS_FREE_SEPARATOR_QUEUE_TEMP();
      return 0;
    }
    thread_counts[target]++;
    private_group_count++;
  }
  if (private_group_count == 0u || private_group_count >= group_count) {
    KLS_FREE_SEPARATOR_QUEUE_TEMP();
    return 0;
  }
  for (int tid = 0; tid < thread_count; ++tid) {
    thread_ptr[tid + 1] = thread_ptr[tid] + thread_counts[tid];
    thread_counts[tid] = thread_ptr[tid];
  }
  if (thread_ptr[thread_count] != private_group_count) {
    KLS_FREE_SEPARATOR_QUEUE_TEMP();
    return 0;
  }
  for (UF_long group = 0; group < group_count; ++group) {
    if (!private_group_mask[group]) {
      continue;
    }
    const int target = component_thread[group_component[group]];
    const UF_long dst = thread_counts[target]++;
    if (dst >= private_group_count) {
      KLS_FREE_SEPARATOR_QUEUE_TEMP();
      return 0;
    }
    private_groups[dst] = group;
  }
  if (pipeline_group_count == 0u ||
      private_group_count + pipeline_group_count != group_count) {
    KLS_FREE_SEPARATOR_QUEUE_TEMP();
    return 0;
  }
  for (int tid = 0; tid < thread_count; ++tid) {
    if (thread_counts[tid] != thread_ptr[tid + 1]) {
      KLS_FREE_SEPARATOR_QUEUE_TEMP();
      return 0;
    }
  }

  UF_long private_threads = 0;
  UF_long private_min_groups = 0;
  UF_long private_max_groups = 0;
  double private_min_work = 0.0;
  double private_max_work = 0.0;
  for (int tid = 0; tid < thread_count; ++tid) {
    const UF_long begin = thread_ptr[tid];
    const UF_long end = thread_ptr[tid + 1];
    if (begin >= end) {
      continue;
    }
    double tid_work = 0.0;
    for (UF_long pos = begin; pos < end; ++pos) {
      tid_work += kls_row_refactor_separator_schedule_work(
        solver, private_groups[pos], separator_work_model);
    }
    const UF_long tid_groups = end - begin;
    if (private_threads == 0u || tid_groups < private_min_groups) {
      private_min_groups = tid_groups;
    }
    if (tid_groups > private_max_groups) {
      private_max_groups = tid_groups;
    }
    if (private_threads == 0u || tid_work < private_min_work) {
      private_min_work = tid_work;
    }
    if (tid_work > private_max_work) {
      private_max_work = tid_work;
    }
    private_threads++;
  }

  if (!kls_ensure_row_refactor_ready_queue_workspace(solver, group_count)) {
    KLS_FREE_SEPARATOR_QUEUE_TEMP();
    return 0;
  }
  UF_long *ready_groups = solver->row_refactor_ready_groups;
  atomic_uint *ready_slots = solver->row_refactor_ready_slots;
  atomic_ulong *remaining_preds = solver->row_refactor_remaining_preds;
  unsigned char *tail_groups = solver->row_refactor_ready_tail_groups;
  for (UF_long group = 0; group < group_count; ++group) {
    atomic_init(&ready_slots[group], 0u);
    atomic_init(&remaining_preds[group],
                (unsigned long)solver->row_refactor_group_pred_count[group]);
    tail_groups[group] = 1u;
  }
  UF_long ready_count = 0;
  for (UF_long group = 0; group < group_count; ++group) {
    if (private_group_mask[group]) {
      continue;
    }
    if (atomic_load_explicit(&remaining_preds[group],
                             memory_order_relaxed) == 0ul) {
      ready_groups[ready_count] = group;
      atomic_store_explicit(&ready_slots[ready_count], 1u,
                            memory_order_relaxed);
      ready_count++;
    }
  }

  free(group_component);
  free(group_forced_pipeline);
  free(component_group_count);
  free(component_work);
  free(subtree_begin);
  free(subtree_end);
  free(subtree_work);
  free(candidates);
  free(component_pipeline);
  free(component_thread);
  free(candidate_entries);
  free(thread_counts);
  free(thread_work);
  *ready_groups_out = ready_groups;
  *ready_slots_out = ready_slots;
  *remaining_preds_out = remaining_preds;
  *tail_groups_out = tail_groups;
  *private_groups_out = private_groups;
  *private_thread_ptr_out = thread_ptr;
  *private_group_mask_out = private_group_mask;
  *tail_count_out = group_count;
  *initial_ready_out = ready_count;
  *private_group_count_out = private_group_count;
  *pipeline_group_count_out = pipeline_group_count;
  *closure_group_count_out = closure_group_count;
  *component_count_out = assigned_component_count + pipeline_component_count;
  *private_threads_out = private_threads;
  *private_min_groups_out = private_min_groups;
  *private_max_groups_out = private_max_groups;
  *private_min_work_out = private_min_work;
  *private_max_work_out = private_max_work;
#undef KLS_FREE_SEPARATOR_QUEUE_TEMP
  return 1;
}

static void kls_invalidate_row_refactor_separator_cache(kls_solver *solver) {
  if (solver == NULL) {
    return;
  }
  free(solver->row_refactor_separator_cache_private_groups);
  free(solver->row_refactor_separator_cache_pipeline_groups);
  free(solver->row_refactor_separator_cache_thread_ptr);
  free(solver->row_refactor_separator_cache_private_mask);
  solver->row_refactor_separator_cache_private_groups = NULL;
  solver->row_refactor_separator_cache_pipeline_groups = NULL;
  solver->row_refactor_separator_cache_thread_ptr = NULL;
  solver->row_refactor_separator_cache_private_mask = NULL;
  solver->row_refactor_separator_cache_thread_count = 0;
  solver->row_refactor_separator_cache_order_checked = 0;
  solver->row_refactor_separator_cache_ordered_private = 0;
  solver->row_refactor_separator_cache_private_count = 0;
  solver->row_refactor_separator_cache_pipeline_count = 0;
  solver->row_refactor_separator_cache_closure_count = 0;
  solver->row_refactor_separator_cache_component_count = 0;
  solver->row_refactor_separator_cache_private_threads = 0;
  solver->row_refactor_separator_cache_min_groups = 0;
  solver->row_refactor_separator_cache_max_groups = 0;
  solver->row_refactor_separator_cache_min_work = 0.0;
  solver->row_refactor_separator_cache_max_work = 0.0;
}

static int kls_install_row_refactor_component_calibration(
  kls_solver *solver,
  const UF_long *private_groups,
  UF_long private_count,
  const double *component_seconds,
  UF_long component_count) {
  if (solver == NULL || private_groups == NULL || private_count == 0u ||
      component_seconds == NULL || component_count == 0u ||
      component_count != solver->separator.component_count ||
      component_count > (UF_long)(SIZE_MAX / sizeof(double))) {
    return 0;
  }
  double *component_work =
    (double *)calloc((size_t)component_count, sizeof(*component_work));
  double *component_scale =
    (double *)malloc((size_t)component_count * sizeof(*component_scale));
  if (component_work == NULL || component_scale == NULL) {
    free(component_work);
    free(component_scale);
    return 0;
  }
  for (UF_long pos = 0; pos < private_count; ++pos) {
    const UF_long group = private_groups[pos];
    unsigned int component = UINT_MAX;
    if (group >= solver->row_refactor_group_count ||
        !kls_row_refactor_group_separator_component_fast(
          solver, group, &component)) {
      free(component_work);
      free(component_scale);
      return 0;
    }
    component_work[component] += kls_row_refactor_group_work(solver, group);
  }

  double total_seconds = 0.0;
  double total_work = 0.0;
  UF_long measured = 0u;
  for (UF_long component = 0; component < component_count; ++component) {
    if (component_work[component] > 0.0 &&
        component_seconds[component] > 0.0 &&
        isfinite(component_seconds[component])) {
      total_seconds += component_seconds[component];
      total_work += component_work[component];
      measured++;
    }
  }
  if (measured < 2u || total_seconds <= 0.0 || total_work <= 0.0) {
    free(component_work);
    free(component_scale);
    return 0;
  }
  const double global_scale = total_seconds / total_work;
  const double min_scale = 0.20 * global_scale;
  const double max_scale = 5.00 * global_scale;
  for (UF_long component = 0; component < component_count; ++component) {
    double scale = global_scale;
    if (component_work[component] > 0.0 &&
        component_seconds[component] > 0.0 &&
        isfinite(component_seconds[component])) {
      scale = component_seconds[component] / component_work[component];
      if (scale < min_scale) {
        scale = min_scale;
      } else if (scale > max_scale) {
        scale = max_scale;
      }
    }
    component_scale[component] = scale;
  }

  free(component_work);
  free(solver->row_refactor_separator_component_scale);
  solver->row_refactor_separator_component_scale = component_scale;
  solver->row_refactor_separator_component_scale_count = component_count;
  solver->row_refactor_separator_component_calibrated = 1;
  if (getenv("KLS_TRACE_ROW_COMPONENT_CALIBRATION") != NULL) {
    fprintf(stderr,
            "KLS row component calibration: measured %ld/%ld, "
            "global %.9g seconds/work\n",
            (long)measured, (long)component_count, global_scale);
    for (UF_long component = 0; component < component_count; ++component) {
      if (component_seconds[component] > 0.0) {
        fprintf(stderr,
                "KLS row component calibration id %ld seconds %.6f "
                "scale %.9g\n",
                (long)component, component_seconds[component],
                component_scale[component]);
      }
    }
  }
  return 1;
}

static int kls_cache_row_refactor_separator_flop_plan(
  kls_solver *solver,
  int thread_count,
  const UF_long *private_groups,
  const UF_long *thread_ptr,
  const unsigned char *private_mask,
  UF_long private_count,
  UF_long pipeline_count,
  UF_long closure_count,
  UF_long component_count,
  UF_long private_threads,
  UF_long min_groups,
  UF_long max_groups,
  double min_work,
  double max_work) {
  if (solver == NULL || thread_count <= 1 || private_groups == NULL ||
      thread_ptr == NULL || private_mask == NULL || private_count == 0u ||
      private_count >= solver->row_refactor_group_count ||
      private_count + pipeline_count != solver->row_refactor_group_count ||
      solver->row_refactor_separator_cache_private_groups != NULL ||
      solver->row_refactor_separator_cache_pipeline_groups != NULL) {
    return 0;
  }
  const size_t private_size = (size_t)private_count;
  const size_t group_size = (size_t)solver->row_refactor_group_count;
  const size_t thread_size = (size_t)thread_count + 1u;
  if ((UF_long)private_size != private_count ||
      (UF_long)group_size != solver->row_refactor_group_count ||
      private_size > SIZE_MAX / sizeof(UF_long) ||
      group_size > SIZE_MAX / sizeof(unsigned char) ||
      thread_size > SIZE_MAX / sizeof(UF_long)) {
    return 0;
  }
  UF_long *groups =
    (UF_long *)malloc(private_size * sizeof(*groups));
  UF_long *pipeline_groups =
    (UF_long *)malloc((size_t)pipeline_count * sizeof(*pipeline_groups));
  UF_long *ptr = (UF_long *)malloc(thread_size * sizeof(*ptr));
  unsigned char *mask =
    (unsigned char *)malloc(group_size * sizeof(*mask));
  if (groups == NULL || pipeline_groups == NULL || ptr == NULL ||
      mask == NULL) {
    free(groups);
    free(pipeline_groups);
    free(ptr);
    free(mask);
    return 0;
  }
  memcpy(groups, private_groups, private_size * sizeof(*groups));
  memcpy(ptr, thread_ptr, thread_size * sizeof(*ptr));
  memcpy(mask, private_mask, group_size * sizeof(*mask));
  UF_long pipeline_pos = 0u;
  for (UF_long group = 0u; group < solver->row_refactor_group_count;
       ++group) {
    if (!private_mask[group]) {
      if (pipeline_pos >= pipeline_count) {
        free(groups);
        free(pipeline_groups);
        free(ptr);
        free(mask);
        return 0;
      }
      pipeline_groups[pipeline_pos++] = group;
    }
  }
  if (pipeline_pos != pipeline_count) {
    free(groups);
    free(pipeline_groups);
    free(ptr);
    free(mask);
    return 0;
  }
  solver->row_refactor_separator_cache_private_groups = groups;
  solver->row_refactor_separator_cache_pipeline_groups = pipeline_groups;
  solver->row_refactor_separator_cache_thread_ptr = ptr;
  solver->row_refactor_separator_cache_private_mask = mask;
  solver->row_refactor_separator_cache_thread_count = thread_count;
  solver->row_refactor_separator_cache_order_checked = 0;
  solver->row_refactor_separator_cache_ordered_private = 0;
  solver->row_refactor_separator_cache_private_count = private_count;
  solver->row_refactor_separator_cache_pipeline_count = pipeline_count;
  solver->row_refactor_separator_cache_closure_count = closure_count;
  solver->row_refactor_separator_cache_component_count = component_count;
  solver->row_refactor_separator_cache_private_threads = private_threads;
  solver->row_refactor_separator_cache_min_groups = min_groups;
  solver->row_refactor_separator_cache_max_groups = max_groups;
  solver->row_refactor_separator_cache_min_work = min_work;
  solver->row_refactor_separator_cache_max_work = max_work;
  return 1;
}

static void *kls_row_schedule_prewarm_main(void *arg) {
  kls_row_schedule_prewarm_job *job =
    (kls_row_schedule_prewarm_job *)arg;
  if (job == NULL || job->solver == NULL || job->thread_count <= 1) {
    return NULL;
  }
  const double begin = kls_now_seconds();
  kls_solver *solver = job->solver;
  const char *cache_env = getenv("KLS_ENABLE_ROW_SEPARATOR_PLAN_CACHE");
  if (cache_env == NULL || cache_env[0] == '\0' ||
      (cache_env[0] == '0' && cache_env[1] == '\0') ||
      solver->row_refactor_separator_cache_private_groups != NULL ||
      solver->row_refactor_separator_cache_pipeline_groups != NULL) {
    job->seconds = kls_now_seconds() - begin;
    return NULL;
  }

  UF_long *ready_groups = NULL;
  atomic_uint *ready_slots = NULL;
  atomic_ulong *remaining_preds = NULL;
  unsigned char *tail_groups = NULL;
  UF_long *private_groups = NULL;
  UF_long *private_thread_ptr = NULL;
  unsigned char *private_group_mask = NULL;
  UF_long tail_count = 0u;
  UF_long initial_ready = 0u;
  UF_long private_count = 0u;
  UF_long pipeline_count = 0u;
  UF_long closure_count = 0u;
  UF_long component_count = 0u;
  UF_long private_threads = 0u;
  UF_long min_groups = 0u;
  UF_long max_groups = 0u;
  double min_work = 0.0;
  double max_work = 0.0;
  const int prepared =
    kls_prepare_row_refactor_separator_flop_ready_queue(
      solver, job->thread_count, &ready_groups, &ready_slots,
      &remaining_preds, &tail_groups, &private_groups,
      &private_thread_ptr, &private_group_mask, &tail_count,
      &initial_ready, &private_count, &pipeline_count, &closure_count,
      &component_count, &private_threads, &min_groups, &max_groups,
      &min_work, &max_work);
  if (prepared) {
    job->ordered_private =
      kls_row_refactor_private_thread_order_valid(
        solver, private_groups, private_thread_ptr, private_group_mask,
        job->thread_count);
    job->cached = kls_cache_row_refactor_separator_flop_plan(
      solver, job->thread_count, private_groups, private_thread_ptr,
      private_group_mask, private_count, pipeline_count, closure_count,
      component_count, private_threads, min_groups, max_groups,
      min_work, max_work);
    if (job->cached) {
      solver->row_refactor_separator_cache_order_checked = 1;
      solver->row_refactor_separator_cache_ordered_private =
        job->ordered_private ? 1 : 0;
    }
    if (job->prewarm_solve_separator) {
      const char *solve_separator_env = getenv("KLS_ROW_SOLVE_SEPARATOR");
      const int solve_separator_enabled =
        solve_separator_env != NULL && solve_separator_env[0] != '\0' &&
        !(solve_separator_env[0] == '0' &&
          solve_separator_env[1] == '\0');
      if (solve_separator_enabled) {
        solver->row_solve_separator_schedule_attempted = 1;
        job->solve_separator_ready =
          kls_update_row_solve_separator_plan_from_row_schedule(
            solver, private_groups, private_thread_ptr,
            private_group_mask, job->thread_count);
      }
    }
  }
  free(private_groups);
  free(private_thread_ptr);
  free(private_group_mask);
  job->seconds = kls_now_seconds() - begin;
  return NULL;
}

static int kls_restore_row_refactor_separator_flop_plan(
  kls_solver *solver,
  int thread_count,
  UF_long **ready_groups_out,
  atomic_uint **ready_slots_out,
  atomic_ulong **remaining_preds_out,
  unsigned char **tail_groups_out,
  UF_long **private_groups_out,
  UF_long **thread_ptr_out,
  unsigned char **private_mask_out,
  UF_long *tail_count_out,
  UF_long *initial_ready_out,
  UF_long *private_count_out,
  UF_long *pipeline_count_out,
  UF_long *closure_count_out,
  UF_long *component_count_out,
  UF_long *private_threads_out,
  UF_long *min_groups_out,
  UF_long *max_groups_out,
  double *min_work_out,
  double *max_work_out) {
  if (solver == NULL || thread_count <= 1 ||
      solver->row_refactor_separator_cache_thread_count != thread_count ||
      solver->row_refactor_separator_cache_private_groups == NULL ||
      solver->row_refactor_separator_cache_pipeline_groups == NULL ||
      solver->row_refactor_separator_cache_thread_ptr == NULL ||
      solver->row_refactor_separator_cache_private_mask == NULL ||
      solver->row_refactor_group_pred_count == NULL ||
      solver->row_refactor_separator_cache_private_count == 0u ||
      solver->row_refactor_separator_cache_private_count +
          solver->row_refactor_separator_cache_pipeline_count !=
        solver->row_refactor_group_count ||
      !kls_ensure_row_refactor_ready_queue_workspace(
        solver, solver->row_refactor_group_count)) {
    return 0;
  }
  const UF_long private_count =
    solver->row_refactor_separator_cache_private_count;
  const UF_long group_count = solver->row_refactor_group_count;
  UF_long *private_groups =
    solver->row_refactor_separator_cache_private_groups;
  UF_long *thread_ptr = solver->row_refactor_separator_cache_thread_ptr;
  unsigned char *private_mask =
    solver->row_refactor_separator_cache_private_mask;
  /* Temporarily transfer ownership to the run.  The normal successful
     cleanup reattaches these arrays; any early error frees them through the
     existing local-plan cleanup and simply drops the cache. */
  solver->row_refactor_separator_cache_private_groups = NULL;
  solver->row_refactor_separator_cache_thread_ptr = NULL;
  solver->row_refactor_separator_cache_private_mask = NULL;

  UF_long *ready_groups = solver->row_refactor_ready_groups;
  atomic_uint *ready_slots = solver->row_refactor_ready_slots;
  atomic_ulong *remaining_preds = solver->row_refactor_remaining_preds;
  unsigned char *tail_groups = solver->row_refactor_ready_tail_groups;
  const char *pipeline_reset_env =
    getenv("KLS_ROW_REFACTOR_CACHE_RESET_PIPELINE_ONLY");
  const int pipeline_reset_only =
    pipeline_reset_env != NULL && pipeline_reset_env[0] != '\0' &&
    !(pipeline_reset_env[0] == '0' && pipeline_reset_env[1] == '\0');
  if (pipeline_reset_only) {
    const UF_long pipeline_count =
      solver->row_refactor_separator_cache_pipeline_count;
    for (UF_long pos = 0; pos < pipeline_count; ++pos) {
      atomic_store_explicit(&ready_slots[pos], 0u, memory_order_relaxed);
      const UF_long group =
        solver->row_refactor_separator_cache_pipeline_groups[pos];
      if (group >= group_count) {
        return 0;
      }
      atomic_store_explicit(
        &remaining_preds[group],
        (unsigned long)solver->row_refactor_group_pred_count[group],
        memory_order_relaxed);
      tail_groups[group] = 1u;
    }
  } else {
    for (UF_long group = 0; group < group_count; ++group) {
      atomic_store_explicit(&ready_slots[group], 0u, memory_order_relaxed);
      atomic_store_explicit(
        &remaining_preds[group],
        (unsigned long)solver->row_refactor_group_pred_count[group],
        memory_order_relaxed);
      tail_groups[group] = 1u;
    }
  }
  UF_long ready_count = 0u;
  if (pipeline_reset_only) {
    const UF_long pipeline_count =
      solver->row_refactor_separator_cache_pipeline_count;
    for (UF_long pos = 0; pos < pipeline_count; ++pos) {
      const UF_long group =
        solver->row_refactor_separator_cache_pipeline_groups[pos];
      if (solver->row_refactor_group_pred_count[group] != 0u) {
        continue;
      }
      ready_groups[ready_count] = group;
      atomic_store_explicit(&ready_slots[ready_count], 1u,
                            memory_order_relaxed);
      ready_count++;
    }
  } else {
    for (UF_long group = 0; group < group_count; ++group) {
      if (!private_mask[group] &&
          solver->row_refactor_group_pred_count[group] == 0u) {
        ready_groups[ready_count] = group;
        atomic_store_explicit(&ready_slots[ready_count], 1u,
                              memory_order_relaxed);
        ready_count++;
      }
    }
  }

  *ready_groups_out = ready_groups;
  *ready_slots_out = ready_slots;
  *remaining_preds_out = remaining_preds;
  *tail_groups_out = tail_groups;
  *private_groups_out = private_groups;
  *thread_ptr_out = thread_ptr;
  *private_mask_out = private_mask;
  *tail_count_out = group_count;
  *initial_ready_out = ready_count;
  *private_count_out = private_count;
  *pipeline_count_out =
    solver->row_refactor_separator_cache_pipeline_count;
  *closure_count_out = solver->row_refactor_separator_cache_closure_count;
  *component_count_out =
    solver->row_refactor_separator_cache_component_count;
  *private_threads_out =
    solver->row_refactor_separator_cache_private_threads;
  *min_groups_out = solver->row_refactor_separator_cache_min_groups;
  *max_groups_out = solver->row_refactor_separator_cache_max_groups;
  *min_work_out = solver->row_refactor_separator_cache_min_work;
  *max_work_out = solver->row_refactor_separator_cache_max_work;
  return 1;
}

static int kls_row_refactor_private_thread_order_valid(
  const kls_solver *solver,
  const UF_long *private_groups,
  const UF_long *private_thread_ptr,
  const unsigned char *private_group_mask,
  int thread_count) {
  if (solver == NULL || private_groups == NULL ||
      private_thread_ptr == NULL || private_group_mask == NULL ||
      thread_count <= 0 || solver->row_refactor_group_count == 0u ||
      solver->row_refactor_group_successor_ptr == NULL ||
      (solver->row_refactor_group_dependency_edges > 0u &&
       solver->row_refactor_group_successor_groups == NULL)) {
    return 0;
  }

  const UF_long group_count = solver->row_refactor_group_count;
  const size_t group_count_size = (size_t)group_count;
  if ((UF_long)group_count_size != group_count ||
      group_count_size > SIZE_MAX / sizeof(int) ||
      group_count_size > SIZE_MAX / sizeof(UF_long) ||
      private_thread_ptr[0] != 0u) {
    return 0;
  }

  int *group_thread =
    (int *)malloc(group_count_size * sizeof(*group_thread));
  UF_long *group_order =
    (UF_long *)malloc(group_count_size * sizeof(*group_order));
  if (group_thread == NULL || group_order == NULL) {
    free(group_thread);
    free(group_order);
    return 0;
  }
  for (UF_long group = 0; group < group_count; ++group) {
    group_thread[group] = -1;
    group_order[group] = KLS_KLU_EMPTY;
  }

  for (int tid = 0; tid < thread_count; ++tid) {
    const UF_long begin = private_thread_ptr[tid];
    const UF_long end = private_thread_ptr[tid + 1];
    if (begin > end || end > group_count) {
      free(group_thread);
      free(group_order);
      return 0;
    }
    for (UF_long pos = begin; pos < end; ++pos) {
      const UF_long group = private_groups[pos];
      if (group >= group_count || !private_group_mask[group] ||
          group_thread[group] >= 0) {
        free(group_thread);
        free(group_order);
        return 0;
      }
      group_thread[group] = tid;
      group_order[group] = pos;
    }
  }

  for (UF_long group = 0; group < group_count; ++group) {
    if (private_group_mask[group] && group_thread[group] < 0) {
      free(group_thread);
      free(group_order);
      return 0;
    }
  }

  for (UF_long group = 0; group < group_count; ++group) {
    const UF_long begin = solver->row_refactor_group_successor_ptr[group];
    const UF_long end =
      solver->row_refactor_group_successor_ptr[group + 1u];
    if (begin > end ||
        end > solver->row_refactor_group_dependency_edges) {
      free(group_thread);
      free(group_order);
      return 0;
    }
    for (UF_long pos = begin; pos < end; ++pos) {
      const UF_long successor =
        solver->row_refactor_group_successor_groups[pos];
      if (successor >= group_count) {
        free(group_thread);
        free(group_order);
        return 0;
      }
      if (!private_group_mask[successor]) {
        continue;
      }
      if (!private_group_mask[group] ||
          group_thread[group] != group_thread[successor] ||
          group_order[group] >= group_order[successor]) {
        free(group_thread);
        free(group_order);
        return 0;
      }
    }
  }

  free(group_thread);
  free(group_order);
  return 1;
}

/* Classify the small boundary where otherwise topological private queues
   depend on a group owned by another worker.  Same-thread edges must already
   follow queue order, and the private set must remain predecessor-closed.
   Callers can then wait directly only at marked targets while avoiding the
   per-edge ready counters used by the fully dynamic private path. */
static unsigned char *
kls_row_refactor_private_external_wait_mask(
  const kls_solver *solver,
  const UF_long *private_groups,
  const UF_long *private_thread_ptr,
  const unsigned char *private_group_mask,
  int thread_count,
  UF_long *external_target_count_out) {
  KLS_SET_OPTIONAL_OUTPUT(external_target_count_out, 0u);
  if (solver == NULL || private_groups == NULL ||
      private_thread_ptr == NULL || private_group_mask == NULL ||
      thread_count <= 0 || external_target_count_out == NULL ||
      solver->row_refactor_group_count == 0u ||
      solver->row_refactor_group_successor_ptr == NULL ||
      (solver->row_refactor_group_dependency_edges > 0u &&
       solver->row_refactor_group_successor_groups == NULL)) {
    return NULL;
  }

  const UF_long group_count = solver->row_refactor_group_count;
  const size_t count = (size_t)group_count;
  if ((UF_long)count != group_count ||
      count > SIZE_MAX / sizeof(int) ||
      count > SIZE_MAX / sizeof(UF_long) ||
      count > SIZE_MAX / sizeof(unsigned char) ||
      private_thread_ptr[0] != 0u ||
      private_thread_ptr[thread_count] > group_count) {
    return NULL;
  }
  int *group_thread = (int *)malloc(count * sizeof(*group_thread));
  UF_long *group_order = (UF_long *)malloc(count * sizeof(*group_order));
  unsigned char *external_wait =
    (unsigned char *)calloc(count, sizeof(*external_wait));
  if (group_thread == NULL || group_order == NULL || external_wait == NULL) {
    free(group_thread);
    free(group_order);
    free(external_wait);
    return NULL;
  }
  for (UF_long group = 0; group < group_count; ++group) {
    group_thread[group] = -1;
    group_order[group] = KLS_KLU_EMPTY;
  }
  for (int tid = 0; tid < thread_count; ++tid) {
    const UF_long begin = private_thread_ptr[tid];
    const UF_long end = private_thread_ptr[tid + 1];
    if (begin > end || end > private_thread_ptr[thread_count]) {
      free(group_thread);
      free(group_order);
      free(external_wait);
      return NULL;
    }
    for (UF_long pos = begin; pos < end; ++pos) {
      const UF_long group = private_groups[pos];
      if (group >= group_count || !private_group_mask[group] ||
          group_thread[group] >= 0) {
        free(group_thread);
        free(group_order);
        free(external_wait);
        return NULL;
      }
      group_thread[group] = tid;
      group_order[group] = pos;
    }
  }
  for (UF_long group = 0; group < group_count; ++group) {
    if (private_group_mask[group] && group_thread[group] < 0) {
      free(group_thread);
      free(group_order);
      free(external_wait);
      return NULL;
    }
  }

  UF_long external_targets = 0u;
  for (UF_long group = 0; group < group_count; ++group) {
    const UF_long begin = solver->row_refactor_group_successor_ptr[group];
    const UF_long end =
      solver->row_refactor_group_successor_ptr[group + 1u];
    if (begin > end ||
        end > solver->row_refactor_group_dependency_edges) {
      free(group_thread);
      free(group_order);
      free(external_wait);
      return NULL;
    }
    for (UF_long pos = begin; pos < end; ++pos) {
      const UF_long successor =
        solver->row_refactor_group_successor_groups[pos];
      if (successor >= group_count) {
        free(group_thread);
        free(group_order);
        free(external_wait);
        return NULL;
      }
      if (!private_group_mask[successor]) {
        continue;
      }
      if (!private_group_mask[group]) {
        free(group_thread);
        free(group_order);
        free(external_wait);
        return NULL;
      }
      if (group_thread[group] == group_thread[successor]) {
        if (group_order[group] >= group_order[successor]) {
          free(group_thread);
          free(group_order);
          free(external_wait);
          return NULL;
        }
      } else if (!external_wait[successor]) {
        external_wait[successor] = 1u;
        external_targets++;
      }
    }
  }
  free(group_thread);
  free(group_order);
  *external_target_count_out = external_targets;
  return external_wait;
}

/* Derive a triangular-solve partition from the row-refactor schedule's group
   ownership.  Refactor workers may use dependency-ready waits, so their queue
   order is not necessarily topological; solve queues are instead emitted in
   factor-row order.  Remove every private group whose L-or-U predecessor is
   in the pipeline or on another thread, then propagate that removal to its
   successors.  The retained private queues are closed under both factor
   graphs: L may run private-then-pipeline, while U may run pipeline-then-
   private in reverse order without cross-thread races. */
static int kls_update_row_solve_separator_plan_from_row_schedule(
  kls_solver *solver,
  const UF_long *private_groups,
  const UF_long *private_thread_ptr,
  const unsigned char *private_group_mask,
  int thread_count) {
  if (solver == NULL || private_groups == NULL ||
      private_thread_ptr == NULL || private_group_mask == NULL ||
      thread_count < 2 || solver->n == 0u ||
      solver->row_refactor_group_count == 0u ||
      solver->row_refactor_group_ptr == NULL ||
      solver->row_refactor_row_group == NULL ||
      solver->row_refactor_group_successor_ptr == NULL ||
      solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_l_cols == NULL ||
      solver->row_refactor_u_ptr == NULL ||
      solver->row_refactor_u_cols == NULL ||
      (solver->row_refactor_group_dependency_edges > 0u &&
       solver->row_refactor_group_successor_groups == NULL)) {
    return 0;
  }

  const UF_long n = solver->n;
  const UF_long group_count = solver->row_refactor_group_count;
  const size_t group_count_size = (size_t)group_count;
  if ((UF_long)group_count_size != group_count ||
      group_count_size > SIZE_MAX / sizeof(int) ||
      group_count_size > SIZE_MAX / sizeof(unsigned char) ||
      private_thread_ptr[0] != 0u ||
      private_thread_ptr[thread_count] > group_count ||
      solver->row_refactor_group_ptr[0] != 0u ||
      solver->row_refactor_group_ptr[group_count] != n) {
    return 0;
  }

  int *group_thread =
    (int *)malloc(group_count_size * sizeof(*group_thread));
  UF_long *remove_queue =
    (UF_long *)malloc(group_count_size * sizeof(*remove_queue));
  unsigned char *private_keep =
    (unsigned char *)malloc(group_count_size * sizeof(*private_keep));
  if (group_thread == NULL || remove_queue == NULL || private_keep == NULL) {
    free(group_thread);
    free(remove_queue);
    free(private_keep);
    return 0;
  }

  for (UF_long group = 0; group < group_count; ++group) {
    const UF_long row_begin = solver->row_refactor_group_ptr[group];
    const UF_long row_end = solver->row_refactor_group_ptr[group + 1u];
    if (row_begin > row_end || row_end > n) {
      free(group_thread);
      free(remove_queue);
      free(private_keep);
      return 0;
    }
    group_thread[group] = -1;
    private_keep[group] = private_group_mask[group] ? 1u : 0u;
  }
  for (int tid = 0; tid < thread_count; ++tid) {
    const UF_long begin = private_thread_ptr[tid];
    const UF_long end = private_thread_ptr[tid + 1];
    if (begin > end || end > private_thread_ptr[thread_count]) {
      free(group_thread);
      free(remove_queue);
      free(private_keep);
      return 0;
    }
    for (UF_long pos = begin; pos < end; ++pos) {
      const UF_long group = private_groups[pos];
      if (group >= group_count || !private_keep[group] ||
          group_thread[group] >= 0) {
        free(group_thread);
        free(remove_queue);
        free(private_keep);
        return 0;
      }
      group_thread[group] = tid;
    }
  }
  for (UF_long group = 0; group < group_count; ++group) {
    if (private_keep[group] && group_thread[group] < 0) {
      free(group_thread);
      free(remove_queue);
      free(private_keep);
      return 0;
    }
  }

  UF_long remove_head = 0u;
  UF_long remove_tail = 0u;
  for (UF_long group = 0; group < group_count; ++group) {
    const UF_long begin = solver->row_refactor_group_successor_ptr[group];
    const UF_long end =
      solver->row_refactor_group_successor_ptr[group + 1u];
    if (begin > end ||
        end > solver->row_refactor_group_dependency_edges) {
      remove_tail = group_count + 1u;
      break;
    }
    for (UF_long pos = begin; pos < end; ++pos) {
      const UF_long target_group =
        solver->row_refactor_group_successor_groups[pos];
      if (target_group >= group_count) {
        remove_tail = group_count + 1u;
        break;
      }
      if (private_keep[target_group] &&
          (!private_keep[group] ||
           group_thread[group] != group_thread[target_group] ||
           group >= target_group)) {
        private_keep[target_group] = 0u;
        remove_queue[remove_tail++] = target_group;
      }
    }
    if (remove_tail > group_count) {
      break;
    }
  }
  for (UF_long row = 0; row < n; ++row) {
    const UF_long source_group = solver->row_refactor_row_group[row];
    if (source_group >= group_count) {
      remove_tail = group_count + 1u;
      break;
    }
    const UF_long begin = solver->row_refactor_u_ptr[row];
    const UF_long end = solver->row_refactor_u_ptr[row + 1u];
    if (begin > end) {
      remove_tail = group_count + 1u;
      break;
    }
    for (UF_long pos = begin; pos < end; ++pos) {
      const UF_long col = solver->row_refactor_u_cols[pos];
      if (col <= row || col >= n) {
        remove_tail = group_count + 1u;
        break;
      }
      const UF_long target_group = solver->row_refactor_row_group[col];
      if (target_group >= group_count) {
        remove_tail = group_count + 1u;
        break;
      }
      if (target_group == source_group || !private_keep[target_group]) {
        continue;
      }
      if (!private_keep[source_group] ||
          group_thread[source_group] != group_thread[target_group] ||
          source_group >= target_group) {
        private_keep[target_group] = 0u;
        remove_queue[remove_tail++] = target_group;
      }
    }
    if (remove_tail > group_count) {
      break;
    }
  }
  if (remove_tail > group_count) {
    free(group_thread);
    free(remove_queue);
    free(private_keep);
    return 0;
  }

  while (remove_head < remove_tail) {
    const UF_long group = remove_queue[remove_head++];
    const UF_long l_begin = solver->row_refactor_group_successor_ptr[group];
    const UF_long l_end =
      solver->row_refactor_group_successor_ptr[group + 1u];
    if (l_begin > l_end ||
        l_end > solver->row_refactor_group_dependency_edges) {
      free(group_thread);
      free(remove_queue);
      free(private_keep);
      return 0;
    }
    for (UF_long pos = l_begin; pos < l_end; ++pos) {
      const UF_long successor =
        solver->row_refactor_group_successor_groups[pos];
      if (successor >= group_count) {
        free(group_thread);
        free(remove_queue);
        free(private_keep);
        return 0;
      }
      if (private_keep[successor]) {
        private_keep[successor] = 0u;
        remove_queue[remove_tail++] = successor;
      }
    }
    const UF_long row_begin = solver->row_refactor_group_ptr[group];
    const UF_long row_end = solver->row_refactor_group_ptr[group + 1u];
    for (UF_long row = row_begin; row < row_end; ++row) {
      for (UF_long pos = solver->row_refactor_u_ptr[row];
           pos < solver->row_refactor_u_ptr[row + 1u]; ++pos) {
        const UF_long col = solver->row_refactor_u_cols[pos];
        if (col >= n) {
          free(group_thread);
          free(remove_queue);
          free(private_keep);
          return 0;
        }
        const UF_long successor = solver->row_refactor_row_group[col];
        if (successor >= group_count) {
          free(group_thread);
          free(remove_queue);
          free(private_keep);
          return 0;
        }
        if (private_keep[successor]) {
          private_keep[successor] = 0u;
          remove_queue[remove_tail++] = successor;
        }
      }
    }
  }

  UF_long private_row_count = 0u;
  UF_long nonempty_threads = 0u;
  UF_long max_thread_rows = 0u;
  UF_long *thread_ptr =
    (UF_long *)calloc((size_t)thread_count + 1u, sizeof(*thread_ptr));
  if (thread_ptr == NULL) {
    free(group_thread);
    free(remove_queue);
    free(private_keep);
    return 0;
  }
  for (int tid = 0; tid < thread_count; ++tid) {
    UF_long thread_rows = 0u;
    for (UF_long group = 0; group < group_count; ++group) {
      if (private_keep[group] && group_thread[group] == tid) {
        thread_rows += solver->row_refactor_group_ptr[group + 1u] -
                       solver->row_refactor_group_ptr[group];
      }
    }
    if (thread_rows > 0u) {
      nonempty_threads++;
    }
    if (thread_rows > max_thread_rows) {
      max_thread_rows = thread_rows;
    }
    if (private_row_count > n - thread_rows) {
      free(thread_ptr);
      free(group_thread);
      free(remove_queue);
      free(private_keep);
      return 0;
    }
    private_row_count += thread_rows;
    thread_ptr[tid + 1] = private_row_count;
  }
  if (private_row_count == 0u || private_row_count > n ||
      nonempty_threads < 2u) {
    free(thread_ptr);
    free(group_thread);
    free(remove_queue);
    free(private_keep);
    return 0;
  }
  const UF_long pipeline_row_count = n - private_row_count;
  UF_long *private_rows =
    (UF_long *)malloc((size_t)private_row_count * sizeof(*private_rows));
  UF_long *pipeline_rows = pipeline_row_count > 0u
    ? (UF_long *)malloc((size_t)pipeline_row_count * sizeof(*pipeline_rows))
    : NULL;
  if (private_rows == NULL ||
      (pipeline_row_count > 0u && pipeline_rows == NULL)) {
    free(private_rows);
    free(pipeline_rows);
    free(thread_ptr);
    free(group_thread);
    free(remove_queue);
    free(private_keep);
    return 0;
  }

  UF_long private_pos = 0u;
  for (int tid = 0; tid < thread_count; ++tid) {
    for (UF_long group = 0; group < group_count; ++group) {
      if (!private_keep[group] || group_thread[group] != tid) {
        continue;
      }
      for (UF_long row = solver->row_refactor_group_ptr[group];
           row < solver->row_refactor_group_ptr[group + 1u]; ++row) {
        private_rows[private_pos++] = row;
      }
    }
  }
  UF_long pipeline_pos = 0u;
  for (UF_long group = 0; group < group_count; ++group) {
    if (private_keep[group]) {
      continue;
    }
    for (UF_long row = solver->row_refactor_group_ptr[group];
         row < solver->row_refactor_group_ptr[group + 1u]; ++row) {
      pipeline_rows[pipeline_pos++] = row;
    }
  }
  if (private_pos != private_row_count ||
      pipeline_pos != pipeline_row_count) {
    free(private_rows);
    free(pipeline_rows);
    free(thread_ptr);
    free(group_thread);
    free(remove_queue);
    free(private_keep);
    return 0;
  }

  free(solver->row_solve_separator_plan.private_rows);
  free(solver->row_solve_separator_plan.pipeline_rows);
  free(solver->row_solve_separator_plan.thread_ptr);
  memset(&solver->row_solve_separator_plan, 0,
         sizeof(solver->row_solve_separator_plan));
  solver->row_solve_separator_plan.private_rows = private_rows;
  solver->row_solve_separator_plan.pipeline_rows = pipeline_rows;
  solver->row_solve_separator_plan.thread_ptr = thread_ptr;
  solver->row_solve_separator_plan.private_rows_count = private_row_count;
  solver->row_solve_separator_plan.pipeline_rows_count = pipeline_row_count;
  solver->row_solve_separator_plan.nonempty_threads = nonempty_threads;
  solver->row_solve_separator_plan.max_thread_rows = max_thread_rows;
  solver->row_solve_separator_plan.partitioned = 1;
  solver->row_solve_separator_thread_count = thread_count;

  if (getenv("KLS_TRACE_ROW_SOLVE_SEPARATOR") != NULL) {
    fprintf(stderr,
            "KLS row-solve separator: private_rows=%" PRIu64
            " pipeline_rows=%" PRIu64 " threads=%" PRIu64
            " removed_groups=%" PRIu64 "\n",
            (uint64_t)private_row_count, (uint64_t)pipeline_row_count,
            (uint64_t)nonempty_threads, (uint64_t)remove_tail);
  }
  free(group_thread);
  free(remove_queue);
  free(private_keep);
  return 1;
}

typedef struct {
  kls_solver *solver;
  const double *group_work;
  UF_long group_begin;
  UF_long group_end;
  int ok;
} kls_row_successor_sort_job;

static void *kls_row_successor_sort_main(void *arg) {
  kls_row_successor_sort_job *job =
    (kls_row_successor_sort_job *)arg;
  if (job == NULL || job->solver == NULL || job->group_work == NULL) {
    return NULL;
  }
  kls_solver *solver = job->solver;
  job->ok = 1;
  for (UF_long group = job->group_begin; group < job->group_end; ++group) {
    const UF_long begin =
      solver->row_refactor_group_successor_ptr[group];
    const UF_long end =
      solver->row_refactor_group_successor_ptr[group + 1u];
    if (begin >= end) {
      continue;
    }
    if (!kls_sort_row_refactor_group_list_by_cached_work(
          solver->row_refactor_group_successor_groups + begin, end - begin,
          job->group_work)) {
      job->ok = 0;
      break;
    }
  }
  return NULL;
}

static void kls_sort_row_refactor_successors_by_work(kls_solver *solver) {
  if (solver == NULL ||
      solver->row_refactor_group_successor_ptr == NULL ||
      solver->row_refactor_group_successor_groups == NULL ||
      solver->row_refactor_group_count == 0u ||
      solver->row_refactor_group_max_fanout <= 1u) {
    return;
  }
  const double *group_work = solver->row_refactor_group_work;
  double *owned_group_work = NULL;
  if (group_work == NULL) {
    if (solver->row_refactor_group_count >
        (UF_long)(SIZE_MAX / sizeof(*owned_group_work))) {
      return;
    }
    owned_group_work =
      (double *)malloc((size_t)solver->row_refactor_group_count *
                       sizeof(*owned_group_work));
    if (owned_group_work == NULL) {
      return;
    }
    for (UF_long group = 0; group < solver->row_refactor_group_count;
         ++group) {
      owned_group_work[group] =
        kls_row_refactor_compute_group_work(solver, group);
    }
    group_work = owned_group_work;
  }
  int thread_count = 1;
  if (getenv("KLS_ENABLE_PARALLEL_SUCCESSOR_SORT") != NULL &&
      solver->options.threads >= 4 &&
      solver->row_refactor_group_count >= 10000u) {
    thread_count = solver->options.threads;
    if (thread_count > 8) {
      thread_count = 8;
    }
  }
  int sort_ok = 1;
  if (thread_count == 1) {
    kls_row_successor_sort_job job = {
      solver, group_work, 0u, solver->row_refactor_group_count, 0
    };
    (void)kls_row_successor_sort_main(&job);
    sort_ok = job.ok;
  } else {
    kls_row_successor_sort_job jobs[8];
    pthread_t threads[8];
    unsigned char launched[8] = { 0u };
    const UF_long groups_per_thread =
      (solver->row_refactor_group_count + (UF_long)thread_count - 1u) /
      (UF_long)thread_count;
    for (int tid = 0; tid < thread_count; ++tid) {
      UF_long begin = (UF_long)tid * groups_per_thread;
      UF_long end = begin + groups_per_thread;
      if (begin > solver->row_refactor_group_count) {
        begin = solver->row_refactor_group_count;
      }
      if (end > solver->row_refactor_group_count) {
        end = solver->row_refactor_group_count;
      }
      jobs[tid] = (kls_row_successor_sort_job){
        solver, group_work, begin, end, 0
      };
    }
    for (int tid = 1; tid < thread_count; ++tid) {
      if (pthread_create(&threads[tid], NULL,
                         kls_row_successor_sort_main, &jobs[tid]) == 0) {
        launched[tid] = 1u;
      } else {
        (void)kls_row_successor_sort_main(&jobs[tid]);
      }
    }
    (void)kls_row_successor_sort_main(&jobs[0]);
    for (int tid = 1; tid < thread_count; ++tid) {
      if (launched[tid]) {
        (void)pthread_join(threads[tid], NULL);
      }
    }
    for (int tid = 0; tid < thread_count; ++tid) {
      sort_ok = sort_ok && jobs[tid].ok;
    }
  }
  if (!sort_ok) {
    free(owned_group_work);
    return;
  }
  if (solver->row_refactor_group_roots != NULL &&
      solver->row_refactor_group_root_count > 1u) {
    (void)kls_sort_row_refactor_group_list_by_cached_work(
      solver->row_refactor_group_roots,
      solver->row_refactor_group_root_count,
      group_work);
  }
  free(owned_group_work);
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

static int kls_ensure_row_refactor_ready_queue_workspace(
  kls_solver *solver,
  UF_long group_count) {
  if (solver == NULL || group_count <= 0) {
    return 0;
  }
  if (solver->row_refactor_ready_queue_capacity >= group_count &&
      solver->row_refactor_ready_groups != NULL &&
      solver->row_refactor_ready_slots != NULL &&
      solver->row_refactor_remaining_preds != NULL &&
      solver->row_refactor_ready_tail_groups != NULL) {
    return 1;
  }
  const size_t count = (size_t)group_count;
  if ((UF_long)count != group_count ||
      count > SIZE_MAX / sizeof(*solver->row_refactor_ready_groups) ||
      count > SIZE_MAX / sizeof(*solver->row_refactor_ready_slots) ||
      count > SIZE_MAX / sizeof(*solver->row_refactor_remaining_preds) ||
      count > SIZE_MAX / sizeof(*solver->row_refactor_ready_tail_groups)) {
    return 0;
  }

  UF_long *ready_groups =
    (UF_long *)malloc(count * sizeof(*ready_groups));
  atomic_uint *ready_slots =
    (atomic_uint *)malloc(count * sizeof(*ready_slots));
  atomic_ulong *remaining_preds =
    (atomic_ulong *)malloc(count * sizeof(*remaining_preds));
  unsigned char *tail_groups =
    (unsigned char *)malloc(count * sizeof(*tail_groups));
  if (ready_groups == NULL || ready_slots == NULL ||
      remaining_preds == NULL || tail_groups == NULL) {
    free(ready_groups);
    free(ready_slots);
    free(remaining_preds);
    free(tail_groups);
    return 0;
  }

  free(solver->row_refactor_ready_groups);
  free(solver->row_refactor_ready_slots);
  free(solver->row_refactor_remaining_preds);
  free(solver->row_refactor_ready_tail_groups);
  solver->row_refactor_ready_groups = ready_groups;
  solver->row_refactor_ready_slots = ready_slots;
  solver->row_refactor_remaining_preds = remaining_preds;
  solver->row_refactor_ready_tail_groups = tail_groups;
  solver->row_refactor_ready_queue_capacity = group_count;
  return 1;
}

static int kls_prepare_row_refactor_ready_queue(
  kls_solver *solver,
  UF_long cluster_levels,
  UF_long **ready_groups_out,
  atomic_uint **ready_slots_out,
  atomic_ulong **remaining_preds_out,
  unsigned char **tail_groups_out,
  const UF_long **private_groups_out,
  UF_long *tail_count_out,
  UF_long *initial_ready_out) {
  KLS_SET_OPTIONAL_OUTPUT(ready_groups_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(ready_slots_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(remaining_preds_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(tail_groups_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(private_groups_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(tail_count_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(initial_ready_out, 0);
  if (solver == NULL || ready_groups_out == NULL ||
      ready_slots_out == NULL || remaining_preds_out == NULL ||
      tail_groups_out == NULL || private_groups_out == NULL ||
      tail_count_out == NULL ||
      initial_ready_out == NULL ||
      solver->row_refactor_group_count == 0u ||
      solver->row_refactor_group_level_ptr == NULL ||
      solver->row_refactor_level_groups == NULL ||
      solver->row_refactor_group_successor_ptr == NULL ||
      (solver->row_refactor_group_dependency_edges > 0u &&
       solver->row_refactor_group_successor_groups == NULL) ||
      cluster_levels >= solver->row_refactor_level_count) {
    return 0;
  }

  const UF_long group_count = solver->row_refactor_group_count;
  const UF_long tail_begin =
    solver->row_refactor_group_level_ptr[cluster_levels];
  if (tail_begin >= group_count) {
    return 0;
  }
  const UF_long tail_count = group_count - tail_begin;
  if (!kls_ensure_row_refactor_ready_queue_workspace(solver, group_count)) {
    return 0;
  }
  UF_long *ready_groups = solver->row_refactor_ready_groups;
  atomic_uint *ready_slots = solver->row_refactor_ready_slots;
  atomic_ulong *remaining_preds = solver->row_refactor_remaining_preds;
  unsigned char *tail_groups = solver->row_refactor_ready_tail_groups;

  for (UF_long pos = 0; pos < tail_count; ++pos) {
    atomic_init(&ready_slots[pos], 0u);
  }

  UF_long ready_count = 0;
  const int full_graph =
    tail_begin == 0u &&
    tail_count == group_count &&
    solver->row_refactor_group_pred_count != NULL &&
    solver->row_refactor_group_roots != NULL &&
    solver->row_refactor_group_root_count > 0u &&
    solver->row_refactor_group_root_count <= group_count;
  if (full_graph) {
    for (UF_long group = 0; group < group_count; ++group) {
      atomic_init(&remaining_preds[group],
                  (unsigned long)solver->row_refactor_group_pred_count[group]);
      tail_groups[group] = 1u;
    }
    ready_count = solver->row_refactor_group_root_count;
    *private_groups_out = solver->row_refactor_group_roots;
    for (UF_long pos = 0; pos < ready_count; ++pos) {
      const UF_long group = solver->row_refactor_group_roots[pos];
      if (group >= group_count) {
        return 0;
      }
    }
  } else {
    for (UF_long group = 0; group < group_count; ++group) {
      atomic_init(&remaining_preds[group], 0ul);
      tail_groups[group] = 0u;
    }
    for (UF_long pos = tail_begin; pos < group_count; ++pos) {
      const UF_long group = solver->row_refactor_level_groups[pos];
      if (group >= group_count) {
        return 0;
      }
      tail_groups[group] = 1u;
    }

    for (UF_long group = 0; group < group_count; ++group) {
      const UF_long begin = solver->row_refactor_group_successor_ptr[group];
      const UF_long end = solver->row_refactor_group_successor_ptr[group + 1u];
      if (begin > end ||
          end > solver->row_refactor_group_dependency_edges) {
        return 0;
      }
      if (!tail_groups[group]) {
        continue;
      }
      for (UF_long pos = begin; pos < end; ++pos) {
        const UF_long successor =
          solver->row_refactor_group_successor_groups[pos];
        if (successor >= group_count) {
          return 0;
        }
        if (tail_groups[successor]) {
          const unsigned long old =
            atomic_load_explicit(&remaining_preds[successor],
                                 memory_order_relaxed);
          atomic_store_explicit(&remaining_preds[successor], old + 1ul,
                                memory_order_relaxed);
        }
      }
    }

    for (UF_long pos = tail_begin; pos < group_count; ++pos) {
      const UF_long group = solver->row_refactor_level_groups[pos];
      if (atomic_load_explicit(&remaining_preds[group],
                               memory_order_relaxed) == 0ul) {
        ready_groups[ready_count] = group;
        atomic_store_explicit(&ready_slots[ready_count], 1u,
                              memory_order_relaxed);
        ready_count++;
      }
    }
  }
  if (ready_count == 0u) {
    return 0;
  }
  if (*private_groups_out == NULL &&
      !kls_sort_row_refactor_group_list_by_work(solver, ready_groups,
                                                ready_count)) {
    return 0;
  }

  *ready_groups_out = ready_groups;
  *ready_slots_out = ready_slots;
  *remaining_preds_out = remaining_preds;
  *tail_groups_out = tail_groups;
  *tail_count_out = tail_count;
  *initial_ready_out = ready_count;
  return 1;
}

static int kls_prepare_row_refactor_row_dep_ready_queue(
  kls_solver *solver,
  UF_long cluster_levels,
  UF_long **ready_groups_out,
  atomic_uint **ready_slots_out,
  atomic_ulong **remaining_preds_out,
  unsigned char **tail_groups_out,
  const UF_long **private_groups_out,
  UF_long *tail_count_out,
  UF_long *initial_ready_out) {
  KLS_SET_OPTIONAL_OUTPUT(ready_groups_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(ready_slots_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(remaining_preds_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(tail_groups_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(private_groups_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(tail_count_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(initial_ready_out, 0);
  if (solver == NULL || ready_groups_out == NULL ||
      ready_slots_out == NULL || remaining_preds_out == NULL ||
      tail_groups_out == NULL || private_groups_out == NULL ||
      tail_count_out == NULL || initial_ready_out == NULL ||
      solver->row_refactor_group_count == 0u ||
      solver->row_refactor_group_level_ptr == NULL ||
      solver->row_refactor_level_groups == NULL ||
      solver->row_refactor_group_ptr == NULL ||
      solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_l_cols == NULL ||
      solver->row_refactor_row_group == NULL ||
      cluster_levels >= solver->row_refactor_level_count) {
    return 0;
  }

  const UF_long group_count = solver->row_refactor_group_count;
  const UF_long tail_begin =
    solver->row_refactor_group_level_ptr[cluster_levels];
  if (tail_begin >= group_count) {
    return 0;
  }
  const UF_long tail_count = group_count - tail_begin;
  if (!kls_ensure_row_refactor_ready_queue_workspace(solver, group_count)) {
    return 0;
  }
  UF_long *ready_groups = solver->row_refactor_ready_groups;
  atomic_uint *ready_slots = solver->row_refactor_ready_slots;
  atomic_ulong *remaining_preds = solver->row_refactor_remaining_preds;
  unsigned char *tail_groups = solver->row_refactor_ready_tail_groups;

  for (UF_long pos = 0; pos < tail_count; ++pos) {
    atomic_init(&ready_slots[pos], 0u);
  }
  for (UF_long group = 0; group < group_count; ++group) {
    atomic_init(&remaining_preds[group], 0ul);
    tail_groups[group] = 0u;
  }
  for (UF_long pos = tail_begin; pos < group_count; ++pos) {
    const UF_long group = solver->row_refactor_level_groups[pos];
    if (group >= group_count) {
      return 0;
    }
    tail_groups[group] = 1u;
  }

  for (UF_long pos = tail_begin; pos < group_count; ++pos) {
    const UF_long group = solver->row_refactor_level_groups[pos];
    const UF_long row_begin = solver->row_refactor_group_ptr[group];
    const UF_long row_end = solver->row_refactor_group_ptr[group + 1u];
    if (row_begin >= row_end || row_end > solver->n) {
      return 0;
    }
    unsigned long pred_edges = 0ul;
    for (UF_long row = row_begin; row < row_end; ++row) {
      const UF_long l_begin = solver->row_refactor_l_ptr[row];
      const UF_long l_end = solver->row_refactor_l_ptr[row + 1u];
      if (l_begin > l_end) {
        return 0;
      }
      for (UF_long p = l_begin; p < l_end; ++p) {
        const UF_long dep = solver->row_refactor_l_cols[p];
        if (dep >= row || dep >= solver->n) {
          return 0;
        }
        const UF_long dep_group = solver->row_refactor_row_group[dep];
        if (dep_group >= group_count) {
          return 0;
        }
        if (dep_group != group && tail_groups[dep_group]) {
          pred_edges++;
        }
      }
    }
    atomic_store_explicit(&remaining_preds[group], pred_edges,
                          memory_order_relaxed);
  }

  UF_long ready_count = 0;
  for (UF_long pos = tail_begin; pos < group_count; ++pos) {
    const UF_long group = solver->row_refactor_level_groups[pos];
    if (atomic_load_explicit(&remaining_preds[group],
                             memory_order_relaxed) == 0ul) {
      ready_groups[ready_count] = group;
      atomic_store_explicit(&ready_slots[ready_count], 1u,
                            memory_order_relaxed);
      ready_count++;
    }
  }
  if (ready_count == 0u ||
      !kls_sort_row_refactor_group_list_by_work(solver, ready_groups,
                                                ready_count)) {
    return 0;
  }

  *ready_groups_out = ready_groups;
  *ready_slots_out = ready_slots;
  *remaining_preds_out = remaining_preds;
  *tail_groups_out = tail_groups;
  *tail_count_out = tail_count;
  *initial_ready_out = ready_count;
  return 1;
}

static int kls_prepare_row_refactor_row_dep_ready_queue_with_private(
  kls_solver *solver,
  const unsigned char *private_group_mask,
  UF_long **ready_groups_out,
  atomic_uint **ready_slots_out,
  atomic_ulong **remaining_preds_out,
  unsigned char **tail_groups_out,
  UF_long *tail_count_out,
  UF_long *initial_ready_out) {
  KLS_SET_OPTIONAL_OUTPUT(ready_groups_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(ready_slots_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(remaining_preds_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(tail_groups_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(tail_count_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(initial_ready_out, 0);
  if (solver == NULL || private_group_mask == NULL ||
      ready_groups_out == NULL || ready_slots_out == NULL ||
      remaining_preds_out == NULL || tail_groups_out == NULL ||
      tail_count_out == NULL || initial_ready_out == NULL ||
      solver->row_refactor_group_count == 0u ||
      solver->row_refactor_group_ptr == NULL ||
      solver->row_refactor_group_pred_count == NULL ||
      solver->row_refactor_group_successor_ptr == NULL ||
      (solver->row_refactor_group_dependency_edges > 0u &&
       solver->row_refactor_group_successor_groups == NULL) ||
      solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_l_cols == NULL ||
      solver->row_refactor_row_group == NULL) {
    return 0;
  }

  const UF_long group_count = solver->row_refactor_group_count;
  if (!kls_ensure_row_refactor_ready_queue_workspace(solver, group_count)) {
    return 0;
  }
  UF_long *ready_groups = solver->row_refactor_ready_groups;
  atomic_uint *ready_slots = solver->row_refactor_ready_slots;
  atomic_ulong *remaining_preds = solver->row_refactor_remaining_preds;
  unsigned char *tail_groups = solver->row_refactor_ready_tail_groups;

  for (UF_long group = 0; group < group_count; ++group) {
    atomic_init(&ready_slots[group], 0u);
    atomic_init(&remaining_preds[group], 0ul);
    tail_groups[group] = 1u;
  }

  for (UF_long group = 0; group < group_count; ++group) {
    const UF_long begin = solver->row_refactor_group_successor_ptr[group];
    const UF_long end = solver->row_refactor_group_successor_ptr[group + 1u];
    if (begin > end ||
        end > solver->row_refactor_group_dependency_edges) {
      return 0;
    }
    if (private_group_mask[group]) {
      continue;
    }
    for (UF_long pos = begin; pos < end; ++pos) {
      const UF_long successor =
        solver->row_refactor_group_successor_groups[pos];
      if (successor >= group_count) {
        return 0;
      }
      if (private_group_mask[successor]) {
        return 0;
      }
    }
  }

  for (UF_long group = 0; group < group_count; ++group) {
    if (private_group_mask[group]) {
      atomic_store_explicit(
        &remaining_preds[group],
        (unsigned long)solver->row_refactor_group_pred_count[group],
        memory_order_relaxed);
      continue;
    }
    const UF_long row_begin = solver->row_refactor_group_ptr[group];
    const UF_long row_end = solver->row_refactor_group_ptr[group + 1u];
    if (row_begin >= row_end || row_end > solver->n) {
      return 0;
    }
    unsigned long pred_edges = 0ul;
    for (UF_long row = row_begin; row < row_end; ++row) {
      const UF_long l_begin = solver->row_refactor_l_ptr[row];
      const UF_long l_end = solver->row_refactor_l_ptr[row + 1u];
      if (l_begin > l_end) {
        return 0;
      }
      for (UF_long p = l_begin; p < l_end; ++p) {
        const UF_long dep = solver->row_refactor_l_cols[p];
        if (dep >= row || dep >= solver->n) {
          return 0;
        }
        const UF_long dep_group = solver->row_refactor_row_group[dep];
        if (dep_group >= group_count) {
          return 0;
        }
        if (dep_group != group) {
          pred_edges++;
        }
      }
    }
    atomic_store_explicit(&remaining_preds[group], pred_edges,
                          memory_order_relaxed);
  }

  UF_long ready_count = 0;
  for (UF_long group = 0; group < group_count; ++group) {
    if (private_group_mask[group]) {
      continue;
    }
    if (atomic_load_explicit(&remaining_preds[group],
                             memory_order_relaxed) == 0ul) {
      ready_groups[ready_count] = group;
      atomic_store_explicit(&ready_slots[ready_count], 1u,
                            memory_order_relaxed);
      ready_count++;
    }
  }
  if (ready_count > 1u &&
      !kls_sort_row_refactor_group_list_by_work(solver, ready_groups,
                                                ready_count)) {
    return 0;
  }

  *ready_groups_out = ready_groups;
  *ready_slots_out = ready_slots;
  *remaining_preds_out = remaining_preds;
  *tail_groups_out = tail_groups;
  *tail_count_out = group_count;
  *initial_ready_out = ready_count;
  return 1;
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

static int kls_row_refactor_tail_has_large_partial_supernodes(
  const kls_solver *solver,
  UF_long cluster_levels,
  int thread_count,
  UF_long *group_count_out,
  UF_long *row_count_out) {
  KLS_SET_OPTIONAL_OUTPUT(group_count_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(row_count_out, 0);
  if (solver == NULL || group_count_out == NULL || row_count_out == NULL ||
      thread_count <= 1 ||
      solver->row_refactor_group_level_ptr == NULL ||
      solver->row_refactor_level_groups == NULL ||
      solver->row_refactor_group_ptr == NULL ||
      solver->row_refactor_group_kind == NULL ||
      solver->row_refactor_group_successor_ptr == NULL ||
      cluster_levels >= solver->row_refactor_level_count) {
    return 0;
  }

  const UF_long tail_begin =
    solver->row_refactor_group_level_ptr[cluster_levels];
  if (tail_begin >= solver->row_refactor_group_count) {
    return 0;
  }
  double tail_work = 0.0;
  for (UF_long pos = tail_begin; pos < solver->row_refactor_group_count;
       ++pos) {
    const UF_long group = solver->row_refactor_level_groups[pos];
    if (group >= solver->row_refactor_group_count) {
      return 0;
    }
    tail_work += kls_row_refactor_group_work(solver, group);
  }
  UF_long groups = 0;
  UF_long rows = 0;
  for (UF_long pos = tail_begin; pos < solver->row_refactor_group_count;
       ++pos) {
    const UF_long group = solver->row_refactor_level_groups[pos];
    if (group >= solver->row_refactor_group_count) {
      return 0;
    }
    const UF_long row_begin = solver->row_refactor_group_ptr[group];
    const UF_long row_end = solver->row_refactor_group_ptr[group + 1u];
    if (row_begin >= row_end || row_end > solver->n) {
      return 0;
    }
    const UF_long width = row_end - row_begin;
    const unsigned char kind = solver->row_refactor_group_kind[group];
    if (kind != KLS_ROW_REFACTOR_GROUP_DENSE &&
        kind != KLS_ROW_REFACTOR_GROUP_GENERIC) {
      continue;
    }
    if (solver->row_refactor_group_successor_ptr[group] ==
        solver->row_refactor_group_successor_ptr[group + 1u]) {
      continue;
    }
    /* Row-prefix publication replaces one group-level handoff with a row-edge
       graph.  Pay that cost only for a producer containing at least one
       worker's fair share of the tail work.  Unlike the former 2P-row rule,
       this remains meaningful when rows have very different update costs. */
    const double group_work = kls_row_refactor_group_work(solver, group);
    if (width < 2u ||
        group_work * (double)thread_count < tail_work) {
      continue;
    }
    groups++;
    rows += width;
  }
  *group_count_out = groups;
  *row_count_out = rows;
  return groups > 0u;
}
