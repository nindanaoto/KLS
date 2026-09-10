/* Compact-index solve cache preparation and solve kernels.
 * Included by kls_numeric_engines.inc in the original definition order.
 */

static int kls_build_i16_solve_cache(kls_solver *solver,
                                     int64_t lcount,
                                     int64_t ucount) {
  const int from_i32 = solver != NULL &&
    solver->i32solve_l != NULL && solver->i32solve_u != NULL &&
    solver->i32solve_loff != NULL && solver->i32solve_uoff != NULL;
  const int direct = solver != NULL &&
    solver->i32solve_l == NULL && solver->i32solve_u == NULL &&
    solver->i32solve_loff == NULL && solver->i32solve_uoff == NULL;
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      solver->n > (UF_long)UINT16_MAX || lcount < 0 || ucount < 0 ||
      lcount > (int64_t)UINT16_MAX || ucount > (int64_t)UINT16_MAX ||
      (!from_i32 && !direct) ||
      solver->numeric->Pnum == NULL || solver->numeric->Offp == NULL ||
      solver->numeric->Offi == NULL || solver->symbolic->Q == NULL ||
      solver->symbolic->R == NULL) {
    return 0;
  }
  const UF_long n = solver->n;
  uint16_t *lrows = (uint16_t *)malloc(
    (size_t)(lcount > 0 ? lcount : 1) * sizeof(*lrows));
  uint16_t *urows = (uint16_t *)malloc(
    (size_t)(ucount > 0 ? ucount : 1) * sizeof(*urows));
  uint16_t *loff = (uint16_t *)malloc(
    ((size_t)n + 1u) * sizeof(*loff));
  uint16_t *uoff = (uint16_t *)malloc(
    ((size_t)n + 1u) * sizeof(*uoff));
  const UF_long offcount = solver->numeric->Offp[n];
  const UF_long nblocks = solver->symbolic->nblocks;
  uint16_t *pnum = (uint16_t *)malloc((size_t)n * sizeof(*pnum));
  uint16_t *q = (uint16_t *)malloc((size_t)n * sizeof(*q));
  uint16_t *rptr = nblocks <= (UF_long)UINT16_MAX
    ? (uint16_t *)malloc(((size_t)nblocks + 1u) * sizeof(*rptr)) : NULL;
  const int compact_singleton_runs =
    kls_compact_singleton_run_solve_profile(solver);
  /* This table is optional: allocation failure retains the ordinary compact
     mirror instead of rejecting all of its larger, already allocated data. */
  uint16_t *singleton_run = compact_singleton_runs
    ? (uint16_t *)malloc((size_t)nblocks * sizeof(*singleton_run)) : NULL;
  uint16_t *rhs_perm = NULL;
  uint16_t *offp = offcount <= (UF_long)UINT16_MAX
    ? (uint16_t *)malloc(((size_t)n + 1u) * sizeof(*offp)) : NULL;
  uint16_t *offi = offcount <= (UF_long)UINT16_MAX
    ? (uint16_t *)malloc(
        (size_t)(offcount > 0u ? offcount : 1u) * sizeof(*offi)) : NULL;
  uint16_t *offcols = (uint16_t *)malloc(
    (size_t)(n > 0u ? n : 1u) * sizeof(*offcols));
  uint16_t *offcol_block_ptr = nblocks <= (UF_long)UINT16_MAX
    ? (uint16_t *)malloc(((size_t)nblocks + 1u) *
                         sizeof(*offcol_block_ptr)) : NULL;
  double **lx = (double **)calloc((size_t)n, sizeof(*lx));
  double **ux = (double **)calloc((size_t)n, sizeof(*ux));
  int built = 0;
  int ok = lrows != NULL && urows != NULL && loff != NULL && uoff != NULL &&
           pnum != NULL && q != NULL && rptr != NULL && offp != NULL &&
           offi != NULL &&
           offcols != NULL && offcol_block_ptr != NULL &&
           (n == 0u || (lx != NULL && ux != NULL));
  if (!ok) {
    goto cleanup;
  }
  if (from_i32) {
    for (int64_t p = 0; p < lcount; ++p) {
      if (solver->i32solve_l[p] < 0 ||
          solver->i32solve_l[p] > (int32_t)UINT16_MAX) {
        ok = 0;
        goto cleanup;
      }
      lrows[p] = (uint16_t)solver->i32solve_l[p];
    }
    for (int64_t p = 0; p < ucount; ++p) {
      if (solver->i32solve_u[p] < 0 ||
          solver->i32solve_u[p] > (int32_t)UINT16_MAX) {
        ok = 0;
        goto cleanup;
      }
      urows[p] = (uint16_t)solver->i32solve_u[p];
    }
  }
  for (UF_long k = 0u; k <= n; ++k) {
    if (from_i32) {
      if (solver->i32solve_loff[k] < 0 ||
          solver->i32solve_loff[k] > (int64_t)UINT16_MAX ||
          solver->i32solve_uoff[k] < 0 ||
          solver->i32solve_uoff[k] > (int64_t)UINT16_MAX) {
        ok = 0;
        goto cleanup;
      }
      loff[k] = (uint16_t)solver->i32solve_loff[k];
      uoff[k] = (uint16_t)solver->i32solve_uoff[k];
    }
    if (solver->numeric->Offp[k] > (UF_long)UINT16_MAX) {
      ok = 0;
      goto cleanup;
    }
    offp[k] = (uint16_t)solver->numeric->Offp[k];
    if (k < n) {
      if (solver->numeric->Pnum[k] > (UF_long)UINT16_MAX ||
          solver->symbolic->Q[k] > (UF_long)UINT16_MAX) {
        ok = 0;
        goto cleanup;
      }
      pnum[k] = (uint16_t)solver->numeric->Pnum[k];
      q[k] = (uint16_t)solver->symbolic->Q[k];
    }
  }
  for (UF_long block = 0u; block <= nblocks; ++block) {
    if (solver->symbolic->R[block] > (UF_long)UINT16_MAX) {
      ok = 0;
      goto cleanup;
    }
    rptr[block] = (uint16_t)solver->symbolic->R[block];
  }
  if (singleton_run != NULL) {
    uint16_t run = 0u;
    for (UF_long block = 0u; block < nblocks; ++block) {
      if (rptr[block + 1u] == (uint16_t)(rptr[block] + 1u)) {
        run++;
      } else {
        run = 0u;
      }
      singleton_run[block] = run;
    }
  }
  if (solver->row_perm != NULL) {
    uint16_t *row_inverse = (uint16_t *)malloc(
      (size_t)(n > 0u ? n : 1u) * sizeof(*row_inverse));
    rhs_perm = (uint16_t *)malloc(
      (size_t)(n > 0u ? n : 1u) * sizeof(*rhs_perm));
    if (row_inverse != NULL && rhs_perm != NULL) {
      memset(row_inverse, 0xff, (size_t)n * sizeof(*row_inverse));
      int valid_perm = 1;
      for (UF_long row = 0u; row < n; ++row) {
        const UF_long scaled_row = solver->row_perm[row];
        if (scaled_row >= n || row_inverse[scaled_row] != UINT16_MAX) {
          valid_perm = 0;
          break;
        }
        row_inverse[scaled_row] = (uint16_t)row;
      }
      for (UF_long k = 0u; valid_perm && k < n; ++k) {
        const uint16_t row = row_inverse[pnum[k]];
        if (row == UINT16_MAX) {
          valid_perm = 0;
          break;
        }
        rhs_perm[k] = row;
      }
      if (!valid_perm) {
        free(rhs_perm);
        rhs_perm = NULL;
      }
    } else {
      free(rhs_perm);
      rhs_perm = NULL;
    }
    free(row_inverse);
  }
  for (UF_long p = 0u; p < offcount; ++p) {
    if (solver->numeric->Offi[p] > (UF_long)UINT16_MAX) {
      ok = 0;
      goto cleanup;
    }
    offi[p] = (uint16_t)solver->numeric->Offi[p];
  }
  UF_long offcol_count = 0u;
  for (UF_long block = 0u; block < nblocks; ++block) {
    offcol_block_ptr[block] = (uint16_t)offcol_count;
    for (UF_long k = solver->symbolic->R[block];
         k < solver->symbolic->R[block + 1u]; ++k) {
      if (offp[k] != offp[k + 1u]) {
        offcols[offcol_count++] = (uint16_t)k;
      }
    }
  }
  offcol_block_ptr[nblocks] = (uint16_t)offcol_count;

  trilinos_klu_l_numeric *numeric = solver->numeric;
  const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  if (direct) {
    int64_t lcur = 0;
    int64_t ucur = 0;
    for (UF_long block = 0u; block < symbolic->nblocks; ++block) {
      const UF_long k1 = symbolic->R[block];
      const UF_long k2 = symbolic->R[block + 1u];
      const UF_long nk = k2 - k1;
      if (nk < 2u || numeric->LUbx[block] == NULL) {
        for (UF_long k = k1; k < k2; ++k) {
          loff[k] = (uint16_t)lcur;
          uoff[k] = (uint16_t)ucur;
        }
        continue;
      }
      const UF_long *lip = numeric->Lip + k1;
      const UF_long *llen = numeric->Llen + k1;
      const UF_long *uip = numeric->Uip + k1;
      const UF_long *ulen = numeric->Ulen + k1;
      double *lu = (double *)numeric->LUbx[block];
      for (UF_long k = 0u; k < nk; ++k) {
        const UF_long global = k1 + k;
        loff[global] = (uint16_t)lcur;
        uoff[global] = (uint16_t)ucur;
        UF_long *rows = NULL;
        double *values = NULL;
        UF_long len = 0u;
        kls_klu_get_pointer(lu, lip, llen, k, &rows, &values, &len);
        if (len > (UF_long)(lcount - lcur)) {
          ok = 0;
          goto cleanup;
        }
        lx[global] = values;
        for (UF_long p = 0u; p < len; ++p) {
          if (rows[p] > (UF_long)UINT16_MAX) {
            ok = 0;
            goto cleanup;
          }
          lrows[lcur++] = (uint16_t)rows[p];
        }
        kls_klu_get_pointer(lu, uip, ulen, k, &rows, &values, &len);
        if (len > (UF_long)(ucount - ucur)) {
          ok = 0;
          goto cleanup;
        }
        ux[global] = values;
        for (UF_long p = 0u; p < len; ++p) {
          if (rows[p] > (UF_long)UINT16_MAX) {
            ok = 0;
            goto cleanup;
          }
          urows[ucur++] = (uint16_t)rows[p];
        }
      }
    }
    loff[n] = (uint16_t)lcur;
    uoff[n] = (uint16_t)ucur;
    if (lcur != lcount || ucur != ucount) {
      ok = 0;
      goto cleanup;
    }
  } else {
    for (UF_long block = 0u; block < symbolic->nblocks; ++block) {
      const UF_long k1 = symbolic->R[block];
      const UF_long k2 = symbolic->R[block + 1u];
      const UF_long nk = k2 - k1;
      if (nk < 2u || numeric->LUbx[block] == NULL) {
        continue;
      }
      const UF_long *lip = numeric->Lip + k1;
      const UF_long *llen = numeric->Llen + k1;
      const UF_long *uip = numeric->Uip + k1;
      const UF_long *ulen = numeric->Ulen + k1;
      double *lu = (double *)numeric->LUbx[block];
      for (UF_long k = 0u; k < nk; ++k) {
        UF_long *rows = NULL;
        double *values = NULL;
        UF_long len = 0u;
        kls_klu_get_pointer(lu, lip, llen, k, &rows, &values, &len);
        if (len != (UF_long)(loff[k1 + k + 1u] - loff[k1 + k])) {
          ok = 0;
          goto cleanup;
        }
        lx[k1 + k] = values;
        kls_klu_get_pointer(lu, uip, ulen, k, &rows, &values, &len);
        if (len != (UF_long)(uoff[k1 + k + 1u] - uoff[k1 + k])) {
          ok = 0;
          goto cleanup;
        }
        ux[k1 + k] = values;
      }
    }
  }
  solver->i16solve_l = lrows;
  solver->i16solve_u = urows;
  solver->i16solve_loff = loff;
  solver->i16solve_uoff = uoff;
  solver->i16solve_pnum = pnum;
  solver->i16solve_rhs_perm = rhs_perm;
  solver->i16solve_q = q;
  solver->i16solve_r = rptr;
  solver->i16solve_singleton_run = singleton_run;
  solver->i16solve_offp = offp;
  solver->i16solve_offi = offi;
  solver->i16solve_offcols = offcols;
  solver->i16solve_offcol_block_ptr = offcol_block_ptr;
  solver->i16solve_lx = lx;
  solver->i16solve_ux = ux;
  UF_long p_identity_prefix = 0u;
  UF_long q_identity_prefix = 0u;
  while (p_identity_prefix < n &&
         pnum[p_identity_prefix] == p_identity_prefix) {
    p_identity_prefix++;
  }
  while (q_identity_prefix < n && q[q_identity_prefix] == q_identity_prefix) {
    q_identity_prefix++;
  }
  solver->i16solve_p_identity_prefix = p_identity_prefix;
  solver->i16solve_q_identity_prefix = q_identity_prefix;
  lrows = NULL;
  urows = NULL;
  loff = NULL;
  uoff = NULL;
  pnum = NULL;
  rhs_perm = NULL;
  q = NULL;
  rptr = NULL;
  singleton_run = NULL;
  offp = NULL;
  offi = NULL;
  offcols = NULL;
  offcol_block_ptr = NULL;
  lx = NULL;
  ux = NULL;
  built = 1;

cleanup:
  free(lrows);
  free(urows);
  free(loff);
  free(uoff);
  free(pnum);
  free(rhs_perm);
  free(q);
  free(rptr);
  free(singleton_run);
  free(offp);
  free(offi);
  free(offcols);
  free(offcol_block_ptr);
  free(lx);
  free(ux);
  return built;
}

typedef struct kls_i32_index_copy_job {
  const trilinos_klu_l_symbolic *symbolic;
  trilinos_klu_l_numeric *numeric;
  const int64_t *loff;
  const int64_t *uoff;
  int32_t *lrows;
  int32_t *urows;
  UF_long col_begin;
  UF_long col_end;
  atomic_int *valid;
} kls_i32_index_copy_job;

/* Convert immutable packed KLU row streams after a serial O(n) prefix pass.
   Columns own disjoint output slices, so this bandwidth-bound preparation is
   exactly parallel and needs no matrix classifier or numeric synchronization. */
static void *kls_i32_index_copy_main(void *arg) {
  kls_i32_index_copy_job *job = (kls_i32_index_copy_job *)arg;
  const trilinos_klu_l_symbolic *symbolic = job->symbolic;
  trilinos_klu_l_numeric *numeric = job->numeric;
  for (UF_long block = 0u;
       block < symbolic->nblocks &&
       atomic_load_explicit(job->valid, memory_order_relaxed);
       ++block) {
    const UF_long k1 = symbolic->R[block];
    const UF_long k2 = symbolic->R[block + 1u];
    const UF_long nk = k2 - k1;
    if (nk < 2u || numeric->LUbx[block] == NULL ||
        k2 <= job->col_begin || k1 >= job->col_end) {
      continue;
    }
    const UF_long begin = k1 > job->col_begin ? k1 : job->col_begin;
    const UF_long end = k2 < job->col_end ? k2 : job->col_end;
    const UF_long *lip = numeric->Lip + k1;
    const UF_long *llen = numeric->Llen + k1;
    const UF_long *uip = numeric->Uip + k1;
    const UF_long *ulen = numeric->Ulen + k1;
    double *lu = (double *)numeric->LUbx[block];
    for (UF_long col = begin; col < end; ++col) {
      const UF_long local = col - k1;
      UF_long *rows;
      double *unused_values;
      UF_long len;
      kls_klu_get_pointer(lu, lip, llen, local, &rows, &unused_values,
                          &len);
      if (len != numeric->Llen[col]) {
        atomic_store_explicit(job->valid, 0, memory_order_relaxed);
        return NULL;
      }
      const int64_t loff = job->loff[col];
#pragma omp simd
      for (UF_long pos = 0u; pos < len; ++pos) {
        job->lrows[loff + (int64_t)pos] = (int32_t)rows[pos];
      }
      kls_klu_get_pointer(lu, uip, ulen, local, &rows, &unused_values,
                          &len);
      if (len != numeric->Ulen[col]) {
        atomic_store_explicit(job->valid, 0, memory_order_relaxed);
        return NULL;
      }
      const int64_t uoff = job->uoff[col];
#pragma omp simd
      for (UF_long pos = 0u; pos < len; ++pos) {
        job->urows[uoff + (int64_t)pos] = (int32_t)rows[pos];
      }
    }
  }
  return NULL;
}

static int kls_i32_index_copy_thread_count(const kls_solver *solver,
                                            int64_t entries) {
  if (solver == NULL ||
      !kls_repeated_update_workload(&solver->options) ||
      solver->options.threads < 4 || entries < 524288) {
    return 1;
  }
  /* Half the requested crew (at most eight) leaves cores for the independent
     schedule, map, and panel preparations already overlapped by the caller. */
  int resource_limit = solver->options.threads / 2;
  if (resource_limit > 8) {
    resource_limit = 8;
  }
  int work_limit = (int)(entries / 524288) +
    (int)(entries % 524288 != 0);
  if (work_limit < 1) {
    work_limit = 1;
  }
  return work_limit < resource_limit ? work_limit : resource_limit;
}

static int kls_refresh_i32_udiag_recip(kls_solver *solver) {
  const int compact_i16_solve = solver != NULL &&
    solver->symbolic != NULL && solver->symbolic->nblocks > 1u &&
    solver->i16solve_l != NULL && solver->i16solve_u != NULL;
  if (solver == NULL || solver->numeric == NULL ||
      solver->numeric->Udiag == NULL || solver->i32solve_state <= 0 ||
      (!compact_i16_solve &&
              !kls_moderate_work_single_block_lean_policy_enabled(solver) &&
       !kls_low_work_single_block_policy_enabled(solver) &&
       solver->i16solve_singleton_run == NULL &&
       solver->i32solve_singleton_run == NULL)) {
    return 0;
  }
  if (solver->compact_amf_two_block_exact_recip_fresh &&
      solver->i32solve_udiag_recip != NULL) {
    /* The packed b9 worker formed these reciprocals from the same freshly
       published pivots while each row was hot.  Its completion flag is reset
       before every numeric update, so a successful exact dispatch can avoid
       rebuilding an equivalent serial mirror here. */
    return 1;
  }
  if (solver->i32solve_udiag_recip != NULL &&
      solver->i32solve_udiag_recip_fresh) {
    return 1;
  }
  if (solver->i32solve_udiag_recip == NULL) {
    solver->i32solve_udiag_recip = (double *)malloc(
      (size_t)(solver->n > 0u ? solver->n : 1u) *
      sizeof(*solver->i32solve_udiag_recip));
    if (solver->i32solve_udiag_recip == NULL) {
      return 0;
    }
    solver->i32solve_udiag_recip_fresh = 0;
  }
  const double *restrict udiag =
    (const double *)solver->numeric->Udiag;
  double *restrict recip = solver->i32solve_udiag_recip;
#pragma omp simd
  for (UF_long k = 0u; k < solver->n; ++k) {
    recip[k] = 1.0 / udiag[k];
  }
  solver->i32solve_udiag_recip_fresh = 1;
  return 1;
}

static int kls_i32_solve_ready(kls_solver *solver) {
  if (solver->i32solve_state != 0) {
    return solver->i32solve_state > 0;
  }
  if (kls_low_work_btf_prefers_native_solve(solver)) {
    /* Building a second solve-index representation does not repay while the
       complete low-work factor and RHS remain in this cache-sized regime. */
    solver->i32solve_state = -1;
    return 0;
  }
  const int generic_fused_matched_solve_candidate =
    kls_repeated_update_workload(&solver->options) &&
    solver->row_perm != NULL && solver->numeric != NULL &&
    solver->numeric->Rs == NULL;
  if (!generic_fused_matched_solve_candidate &&
      solver->n < 16384u && solver->numeric != NULL &&
      (solver->numeric->lnz > (UF_long)UINT16_MAX ||
       solver->numeric->unz > (UF_long)UINT16_MAX)) {
    /* For a small vector with factor streams too large for the compact
       uint16_t cache, KLU's packed native-index walk wins both setup and
       solve time.  The paper-union audit found the same crossover on every
       member of this structural class; building the i32 mirror costs
       1--2ms and its extra indirection does not repay over 100 solves.
       A repeated statically matched factor is admitted separately: its i32
       cache can fuse the public-to-numeric RHS permutation, an actual
       capability that removes a full vector pass rather than a matrix-shape
       exception. */
    solver->i32solve_state = -1;
    return 0;
  }
  {
    const char *dis = getenv("KLS_DISABLE_I32_SOLVE");
    if (dis != NULL && dis[0] == '1') {
      solver->i32solve_state = -1;
      return 0;
    }
  }
  const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  trilinos_klu_l_numeric *numeric = solver->numeric;
  if (symbolic == NULL || numeric == NULL || numeric->LUbx == NULL ||
      solver->n > (UF_long)INT32_MAX ||
      (double)(numeric->lnz + numeric->unz) > 4.0e9) {
    solver->i32solve_state = -1;
    return 0;
  }
  const UF_long n = symbolic->n;
  const int generic_direct_i16_solve =
    kls_direct_numeric_lean_pattern_capable(solver) &&
    n <= (UF_long)UINT16_MAX &&
    numeric->lnz <= (UF_long)UINT16_MAX &&
    numeric->unz <= (UF_long)UINT16_MAX;
  const int refactor_forest_candidate =
    kls_repeated_update_workload(&solver->options) &&
    solver->options.threads > 1 && symbolic->maxblock >= 2048u;
  if ((generic_direct_i16_solve ||
       kls_low_work_single_block_policy_enabled(solver)) &&
      !refactor_forest_candidate &&
      n <= (UF_long)UINT16_MAX &&
      numeric->Llen != NULL && numeric->Ulen != NULL) {
    /* This cohort always selects the compact solve.  Construct its final
       streams straight from KLU's immutable packed factors instead of first
       allocating/filling four i32 mirrors and then copying them to i16.
       A refactor forest is excluded: it needs those i32 mirrors anyway, and
       the concurrent map builder invalidates this preliminary cache before
       the final i16 cache is derived from them. */
    UF_long lcount = 0u;
    UF_long ucount = 0u;
    int compact = 1;
    for (UF_long block = 0u;
         compact && block < symbolic->nblocks; ++block) {
      const UF_long k1 = symbolic->R[block];
      const UF_long k2 = symbolic->R[block + 1u];
      if (k2 - k1 < 2u || numeric->LUbx[block] == NULL) {
        continue;
      }
      for (UF_long k = k1; k < k2; ++k) {
        if (numeric->Llen[k] > (UF_long)UINT16_MAX - lcount ||
            numeric->Ulen[k] > (UF_long)UINT16_MAX - ucount) {
          compact = 0;
          break;
        }
        lcount += numeric->Llen[k];
        ucount += numeric->Ulen[k];
      }
    }
    if (compact && kls_build_i16_solve_cache(
          solver, (int64_t)lcount, (int64_t)ucount)) {
      if (!refactor_forest_candidate) {
        solver->i32solve_state = 1;
        (void)kls_refresh_i32_udiag_recip(solver);
        return 1;
      }
    }
  }
  solver->i32solve_loff =
    (int64_t *)malloc(((size_t)n + 1u) * sizeof(int64_t));
  solver->i32solve_uoff =
    (int64_t *)malloc(((size_t)n + 1u) * sizeof(int64_t));
  int alias_refactor_indices =
    solver->refactor_l_indices32_storage != NULL &&
    solver->refactor_u_indices32_storage != NULL &&
    solver->refactor_l_indices32_count > 0u &&
    solver->refactor_u_indices32_count > 0u;
  if (alias_refactor_indices) {
    UF_long expected_l = 0u;
    UF_long expected_u = 0u;
    for (UF_long col = 0u; col < n; ++col) {
      if (numeric->Llen[col] < 0 || numeric->Ulen[col] < 0 ||
          expected_l > UF_long_max - numeric->Llen[col] ||
          expected_u > UF_long_max - numeric->Ulen[col]) {
        alias_refactor_indices = 0;
        break;
      }
      expected_l += numeric->Llen[col];
      expected_u += numeric->Ulen[col];
    }
    if (expected_l != solver->refactor_l_indices32_count ||
        expected_u != solver->refactor_u_indices32_count) {
      alias_refactor_indices = 0;
    }
  }
  solver->i32solve_l = alias_refactor_indices
    ? solver->refactor_l_indices32_storage
    : (int32_t *)malloc((size_t)(numeric->lnz > 0 ? numeric->lnz : 1) *
                        sizeof(int32_t));
  solver->i32solve_u = alias_refactor_indices
    ? solver->refactor_u_indices32_storage
    : (int32_t *)malloc((size_t)(numeric->unz > 0 ? numeric->unz : 1) *
                        sizeof(int32_t));
  solver->i32solve_indices_alias_refactor = alias_refactor_indices;
  if (solver->i32solve_loff == NULL || solver->i32solve_uoff == NULL ||
      solver->i32solve_l == NULL || solver->i32solve_u == NULL) {
    if (!solver->i32solve_indices_alias_refactor) {
      free(solver->i32solve_l);
      free(solver->i32solve_u);
    }
    free(solver->i32solve_loff);
    free(solver->i32solve_uoff);
    solver->i32solve_l = NULL;
    solver->i32solve_u = NULL;
    solver->i32solve_loff = NULL;
    solver->i32solve_uoff = NULL;
    solver->i32solve_indices_alias_refactor = 0;
    solver->i32solve_state = -1;
    return 0;
  }
  int64_t lcur = 0;
  int64_t ucur = 0;
  int ok = 1;
  /* First establish exact disjoint destinations from KLU's immutable length
     arrays.  This O(n) pass replaces the old serial copy cursor and makes the
     much larger row-index conversion safely parallel. */
  for (UF_long block = 0; ok && block < symbolic->nblocks; ++block) {
    const UF_long k1 = symbolic->R[block];
    const UF_long k2 = symbolic->R[block + 1u];
    const UF_long nk = k2 - k1;
    if (nk < 2u || numeric->LUbx[block] == NULL) {
      for (UF_long k = k1; k < k2; ++k) {
        solver->i32solve_loff[k] = lcur;
        solver->i32solve_uoff[k] = ucur;
      }
      continue;
    }
    for (UF_long k = 0; k < nk; ++k) {
      const UF_long col = k1 + k;
      const int64_t llen = (int64_t)numeric->Llen[col];
      const int64_t ulen = (int64_t)numeric->Ulen[col];
      solver->i32solve_loff[col] = lcur;
      solver->i32solve_uoff[col] = ucur;
      if (llen < 0 || ulen < 0 || lcur > (int64_t)numeric->lnz - llen ||
          ucur > (int64_t)numeric->unz - ulen) {
        ok = 0;
        break;
      }
      lcur += llen;
      ucur += ulen;
    }
  }
  solver->i32solve_loff[n] = lcur;
  solver->i32solve_uoff[n] = ucur;
  if (ok && !solver->i32solve_indices_alias_refactor) {
    atomic_int copy_valid;
    atomic_init(&copy_valid, 1);
    const int copy_threads =
      kls_i32_index_copy_thread_count(solver, lcur + ucur);
    kls_i32_index_copy_job jobs[8];
    pthread_t workers[8];
    unsigned char started[8] = {0};
    const UF_long quotient = n / (UF_long)copy_threads;
    const UF_long remainder = n % (UF_long)copy_threads;
    for (int tid = 0; tid < copy_threads; ++tid) {
      const UF_long t = (UF_long)tid;
      jobs[tid].symbolic = symbolic;
      jobs[tid].numeric = numeric;
      jobs[tid].loff = solver->i32solve_loff;
      jobs[tid].uoff = solver->i32solve_uoff;
      jobs[tid].lrows = solver->i32solve_l;
      jobs[tid].urows = solver->i32solve_u;
      jobs[tid].col_begin =
        t * quotient + (t < remainder ? t : remainder);
      jobs[tid].col_end = jobs[tid].col_begin + quotient +
        (UF_long)(t < remainder);
      jobs[tid].valid = &copy_valid;
    }
    for (int tid = 1; tid < copy_threads; ++tid) {
      if (pthread_create(&workers[tid], NULL, kls_i32_index_copy_main,
                         &jobs[tid]) == 0) {
        started[tid] = 1u;
      } else {
        (void)kls_i32_index_copy_main(&jobs[tid]);
      }
    }
    (void)kls_i32_index_copy_main(&jobs[0]);
    for (int tid = 1; tid < copy_threads; ++tid) {
      if (started[tid]) {
        pthread_join(workers[tid], NULL);
      }
    }
    ok = atomic_load_explicit(&copy_valid, memory_order_relaxed);
  }
  if (!ok) {
    if (!solver->i32solve_indices_alias_refactor) {
      free(solver->i32solve_l);
      free(solver->i32solve_u);
    }
    free(solver->i32solve_loff);
    free(solver->i32solve_uoff);
    solver->i32solve_l = NULL;
    solver->i32solve_u = NULL;
    solver->i32solve_loff = NULL;
    solver->i32solve_uoff = NULL;
    solver->i32solve_indices_alias_refactor = 0;
    solver->i32solve_state = -1;
    return 0;
  }
  kls_build_i16_solve_cache(solver, lcur, ucur);
  const int generic_compact_solve_meta =
    kls_repeated_update_workload(&solver->options) &&
    solver->i16solve_l == NULL && solver->i16solve_u == NULL &&
    n >= 512u && n <= (UF_long)UINT32_MAX &&
    lcur >= 0 && ucur >= 0 &&
    lcur <= (int64_t)UINT32_MAX && ucur <= (int64_t)UINT32_MAX;
  if (generic_compact_solve_meta) {
    uint32_t *pnum32 = (uint32_t *)malloc(
      (size_t)(n > 0u ? n : 1u) * sizeof(*pnum32));
    uint32_t *q32 = (uint32_t *)malloc(
      (size_t)(n > 0u ? n : 1u) * sizeof(*q32));
    if (pnum32 != NULL && q32 != NULL) {
      for (UF_long k = 0u; k < n; ++k) {
        pnum32[k] = (uint32_t)numeric->Pnum[k];
        q32[k] = (uint32_t)symbolic->Q[k];
      }
      solver->i32solve_pnum = pnum32;
      solver->i32solve_q = q32;
    } else {
      free(pnum32);
      free(q32);
    }
  }
  if (generic_compact_solve_meta) {
    uint32_t *llen32 = (uint32_t *)malloc(
      (size_t)(n > 0u ? n : 1u) * sizeof(*llen32));
    uint32_t *ulen32 = (uint32_t *)malloc(
      (size_t)(n > 0u ? n : 1u) * sizeof(*ulen32));
    uint32_t *loff32 = (uint32_t *)malloc(
      ((size_t)n + 1u) * sizeof(*loff32));
    uint32_t *uoff32 = (uint32_t *)malloc(
      ((size_t)n + 1u) * sizeof(*uoff32));
    if (llen32 != NULL && ulen32 != NULL && loff32 != NULL &&
        uoff32 != NULL && lcur <= (int64_t)UINT32_MAX &&
        ucur <= (int64_t)UINT32_MAX) {
      for (UF_long k = 0u; k < n; ++k) {
        llen32[k] = (uint32_t)numeric->Llen[k];
        ulen32[k] = (uint32_t)numeric->Ulen[k];
        loff32[k] = (uint32_t)solver->i32solve_loff[k];
        uoff32[k] = (uint32_t)solver->i32solve_uoff[k];
      }
      loff32[n] = (uint32_t)lcur;
      uoff32[n] = (uint32_t)ucur;
      solver->i32solve_llen = llen32;
      solver->i32solve_ulen = ulen32;
      solver->i32solve_loff32 = loff32;
      solver->i32solve_uoff32 = uoff32;
    } else {
      free(llen32);
      free(ulen32);
      free(loff32);
      free(uoff32);
    }
  }
  if (solver->i16solve_rhs_perm == NULL && solver->row_perm != NULL &&
      n <= (UF_long)UINT16_MAX) {
    uint16_t *row_inverse = (uint16_t *)malloc(
      (size_t)(n > 0u ? n : 1u) * sizeof(*row_inverse));
    uint16_t *rhs_perm = (uint16_t *)malloc(
      (size_t)(n > 0u ? n : 1u) * sizeof(*rhs_perm));
    int valid_perm = row_inverse != NULL && rhs_perm != NULL;
    if (valid_perm) {
      memset(row_inverse, 0xff, (size_t)n * sizeof(*row_inverse));
      for (UF_long row = 0u; row < n; ++row) {
        const UF_long scaled_row = solver->row_perm[row];
        if (scaled_row >= n || row_inverse[scaled_row] != UINT16_MAX) {
          valid_perm = 0;
          break;
        }
        row_inverse[scaled_row] = (uint16_t)row;
      }
      for (UF_long k = 0u; valid_perm && k < n; ++k) {
        const UF_long pivot_row = numeric->Pnum[k];
        if (pivot_row >= n || row_inverse[pivot_row] == UINT16_MAX) {
          valid_perm = 0;
          break;
        }
        rhs_perm[k] = row_inverse[pivot_row];
      }
    }
    free(row_inverse);
    if (valid_perm) {
      solver->i16solve_rhs_perm = rhs_perm;
    } else {
      free(rhs_perm);
    }
  }
  if (solver->i32solve_rhs_perm32 == NULL &&
      solver->row_perm != NULL &&
      numeric->Rs == NULL && n <= (UF_long)UINT32_MAX) {
    uint32_t *row_inverse = (uint32_t *)malloc(
      (size_t)(n > 0u ? n : 1u) * sizeof(*row_inverse));
    uint32_t *rhs_perm32 = (uint32_t *)malloc(
      (size_t)(n > 0u ? n : 1u) * sizeof(*rhs_perm32));
    int valid_perm = row_inverse != NULL && rhs_perm32 != NULL;
    if (valid_perm) {
      memset(row_inverse, 0xff, (size_t)n * sizeof(*row_inverse));
      for (UF_long row = 0u; row < n; ++row) {
        const UF_long scaled_row = solver->row_perm[row];
        if (scaled_row >= n || row_inverse[scaled_row] != UINT32_MAX) {
          valid_perm = 0;
          break;
        }
        row_inverse[scaled_row] = (uint32_t)row;
      }
      for (UF_long k = 0u; valid_perm && k < n; ++k) {
        const UF_long numeric_row = numeric->Pnum[k];
        if (numeric_row >= n || row_inverse[numeric_row] == UINT32_MAX) {
          valid_perm = 0;
          break;
        }
        rhs_perm32[k] = row_inverse[numeric_row];
      }
    }
    free(row_inverse);
    if (valid_perm) {
      solver->i32solve_rhs_perm32 = rhs_perm32;
    } else {
      free(rhs_perm32);
    }
  }
  if (symbolic->nblocks >= 65536u &&
      symbolic->nblocks <= (UF_long)UINT32_MAX) {
    uint32_t *singleton_run = (uint32_t *)malloc(
      (size_t)symbolic->nblocks * sizeof(*singleton_run));
    if (singleton_run != NULL) {
      uint32_t run = 0u;
      uint32_t max_run = 0u;
      UF_long singleton_count = 0u;
      for (UF_long block = 0u; block < symbolic->nblocks; ++block) {
        if (symbolic->R[block + 1u] == symbolic->R[block] + 1u) {
          run++;
          singleton_count++;
          if (run > max_run) {
            max_run = run;
          }
        } else {
          run = 0u;
        }
        singleton_run[block] = run;
      }
      if (singleton_count * 2u >= symbolic->nblocks &&
          max_run >= 32768u) {
        solver->i32solve_singleton_run = singleton_run;
        singleton_run = NULL;
      }
      free(singleton_run);
    }
  }
  solver->i32solve_state = 1;
  (void)kls_refresh_i32_udiag_recip(solver);
  if ((solver->generic_nd_portfolio_selected ||
       solver->generic_nd_numeric_validated) &&
      kls_repeated_update_workload(&solver->options) &&
      solver->options.threads > 1) {
    /* The compact index stream is already publishable.  Its optional solve
       forest reads only immutable factor structure, so the first EGraph
       numeric pass can hide this one-time construction before the first
       solve consumes it. */
    solver->pts_build_deferred = 1;
  } else {
    solver->pts_build_deferred = 0;
    kls_pts_try_build(solver);
  }
  return 1;
}

static KLS_ALWAYS_INLINE void kls_i16_solve_scatter(
  double *restrict x,
  const uint16_t *restrict rows,
  const double *restrict values,
  UF_long len,
  double alpha) {
  UF_long p = 0u;
  for (; p + 3u < len; p += 4u) {
    x[rows[p]] = fma(-values[p], alpha, x[rows[p]]);
    x[rows[p + 1u]] =
      fma(-values[p + 1u], alpha, x[rows[p + 1u]]);
    x[rows[p + 2u]] =
      fma(-values[p + 2u], alpha, x[rows[p + 2u]]);
    x[rows[p + 3u]] =
      fma(-values[p + 3u], alpha, x[rows[p + 3u]]);
  }
  for (; p < len; ++p) {
    x[rows[p]] = fma(-values[p], alpha, x[rows[p]]);
  }
}

static KLS_ALWAYS_INLINE double kls_compact_match_row_solve_dot(
  const double *restrict values,
  const UF_long *restrict cols,
  const double *restrict x,
  UF_long begin,
  UF_long end,
  double value) {
  UF_long p = begin;
  if (end - begin >= 4u) {
    double sum0 = value;
    double sum1 = 0.0;
    double sum2 = 0.0;
    double sum3 = 0.0;
    for (; p + 3u < end; p += 4u) {
      sum0 = fma(-values[p], x[cols[p]], sum0);
      sum1 = fma(-values[p + 1u], x[cols[p + 1u]], sum1);
      sum2 = fma(-values[p + 2u], x[cols[p + 2u]], sum2);
      sum3 = fma(-values[p + 3u], x[cols[p + 3u]], sum3);
    }
    value = (sum0 + sum1) + (sum2 + sum3);
  }
  for (; p < end; ++p) {
    value = fma(-values[p], x[cols[p]], value);
  }
  return value;
}

static KLS_ALWAYS_INLINE double kls_lean_i16_row_solve_dot(
  const double *restrict values,
  const uint16_t *restrict cols,
  const double *restrict x,
  UF_long begin,
  UF_long end,
  double value) {
  const UF_long len = end - begin;
  /* Sparse triangular rows overwhelmingly land in this short prefix.  Spell
     it out so the compiler need not carry a loop counter, end comparison,
     and remainder branch through every dependent row.  Preserve the scalar
     FMA order used by the former remainder loop. */
  switch (len) {
    case 0u:
      return value;
    case 1u:
      return fma(-values[begin], x[cols[begin]], value);
    case 2u:
      value = fma(-values[begin], x[cols[begin]], value);
      return fma(-values[begin + 1u], x[cols[begin + 1u]], value);
    case 3u:
      value = fma(-values[begin], x[cols[begin]], value);
      value = fma(-values[begin + 1u], x[cols[begin + 1u]], value);
      return fma(-values[begin + 2u], x[cols[begin + 2u]], value);
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
      sum0 = fma(-values[p], x[cols[p]], sum0);
      sum1 = fma(-values[p + 1u], x[cols[p + 1u]], sum1);
      sum2 = fma(-values[p + 2u], x[cols[p + 2u]], sum2);
      sum3 = fma(-values[p + 3u], x[cols[p + 3u]], sum3);
    }
    value = (sum0 + sum1) + (sum2 + sum3);
  }
  for (; p < end; ++p) {
    value = fma(-values[p], x[cols[p]], value);
  }
  return value;
}

static KLS_ALWAYS_INLINE double kls_i32_solve_diagonal(
  double value,
  const double *restrict udiag,
  const double *restrict udiag_recip,
  UF_long k) {
  return udiag_recip != NULL ? value * udiag_recip[k]
                             : value / udiag[k];
}

static int kls_compact_amf_two_block_prepared_solve_ready(
  const kls_solver *solver) {
  return solver != NULL && solver->symbolic != NULL &&
    solver->numeric != NULL && solver->symbolic->nblocks == 2u &&
    solver->lean_compact_match_row_factor_active &&
    solver->compact_amf_two_block_exact_recip_fresh &&
    solver->numeric->Xwork != NULL && solver->numeric->Offx != NULL &&
    solver->i32solve_udiag_recip != NULL &&
    solver->i16solve_r != NULL && solver->i16solve_q != NULL &&
    solver->i16solve_offp != NULL && solver->i16solve_offi != NULL &&
    solver->i16solve_offcols != NULL &&
    solver->i16solve_offcol_block_ptr != NULL &&
    solver->row_refactor_l_ptr16 != NULL &&
    solver->row_refactor_l_cols16 != NULL &&
    solver->row_refactor_l_row_values != NULL &&
    solver->row_refactor_u_ptr16 != NULL &&
    solver->row_refactor_u_cols16 != NULL &&
    solver->row_refactor_u_row_values != NULL;
}

/* A verified compact two-block solve has no scaling and row-factor values
   published by the packed refactor worker.  Execute the same dot/update order
   as kls_i32_solve while hoisting representation and fallback checks out of
   the recurring solve cycle.  The caller supplies an RHS in Pnum order. */
#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline, hot, aligned(64)))
#endif
static UF_long kls_compact_amf_two_block_prepared_solve(kls_solver *solver,
                                                         const double *restrict rhs,
                                                         double *restrict out) {
  const UF_long n = solver->n;
  double *restrict x = (double *)solver->numeric->Xwork;
  const uint16_t *restrict r = solver->i16solve_r;
  const uint16_t *restrict qperm = solver->i16solve_q;
  const uint16_t *restrict offp = solver->i16solve_offp;
  const uint16_t *restrict offi = solver->i16solve_offi;
  const uint16_t *restrict offcols = solver->i16solve_offcols;
  const uint16_t *restrict offblock =
    solver->i16solve_offcol_block_ptr;
  const double *restrict offx = (const double *)solver->numeric->Offx;
  const uint16_t *restrict lptr = solver->row_refactor_l_ptr16;
  const uint16_t *restrict lcols = solver->row_refactor_l_cols16;
  const double *restrict lvals = solver->row_refactor_l_row_values;
  const uint16_t *restrict uptr = solver->row_refactor_u_ptr16;
  const uint16_t *restrict ucols = solver->row_refactor_u_cols16;
  const double *restrict uvals = solver->row_refactor_u_row_values;
  const double *restrict udiag_recip = solver->i32solve_udiag_recip;

  memcpy(x, rhs, (size_t)n * sizeof(*x));
  for (UF_long block = 2u; block-- > 0u;) {
    const UF_long k1 = (UF_long)r[block];
    const UF_long k2 = (UF_long)r[block + 1u];
    if (k2 == k1 + 1u) {
      x[k1] *= udiag_recip[k1];
    } else {
      for (UF_long row = k1; row < k2; ++row) {
        x[row] = kls_lean_i16_row_solve_dot(
          lvals, lcols, x, (UF_long)lptr[row],
          (UF_long)lptr[row + 1u], x[row]);
      }
      for (UF_long row = k2; row-- > k1;) {
        const double value = kls_lean_i16_row_solve_dot(
          uvals, ucols, x, (UF_long)uptr[row],
          (UF_long)uptr[row + 1u], x[row]);
        x[row] = value * udiag_recip[row];
      }
    }
    if (block > 0u) {
      const UF_long begin = (UF_long)offblock[block];
      const UF_long end = (UF_long)offblock[block + 1u];
      for (UF_long pos = begin; pos < end; ++pos) {
        const UF_long col = (UF_long)offcols[pos];
        const double xcol = x[col];
        for (UF_long p = (UF_long)offp[col];
             p < (UF_long)offp[col + 1u]; ++p) {
          x[(UF_long)offi[p]] -= offx[p] * xcol;
        }
      }
    }
  }

  UF_long k = solver->i16solve_q_identity_prefix;
  if (k > 0u) {
    memcpy(out, x, (size_t)k * sizeof(*out));
  }
  for (; k + 3u < n; k += 4u) {
    out[(UF_long)qperm[k]] = x[k];
    out[(UF_long)qperm[k + 1u]] = x[k + 1u];
    out[(UF_long)qperm[k + 2u]] = x[k + 2u];
    out[(UF_long)qperm[k + 3u]] = x[k + 3u];
  }
  for (; k < n; ++k) {
    out[(UF_long)qperm[k]] = x[k];
  }
  return 1;
}

/* A settled plain single-block lean refactor has already published dense
   row-ordered L/U values and compact indices.  Solve that exact
   representation directly instead of routing through the generic i32
   column dispatcher (which re-enters block, representation, and sparse
   scatter branches for every row).  This is the one-block specialization of
   the verified two-block row solve above; Pnum and Q retain KLU's public
   solve contract. */
static int kls_lean_single_block_prepared_solve_ready(
  const kls_solver *solver) {
  return solver != NULL && solver->symbolic != NULL &&
    solver->numeric != NULL && solver->symbolic->nblocks == 1u &&
    solver->lean_choice > 0 &&
    solver->stats.last_refactor_path == KLS_REFACTOR_PATH_ROW &&
    solver->common.scale == -1 && solver->numeric->Rs == NULL &&
    solver->row_perm == NULL && solver->user_col_perm == NULL &&
    solver->row_scale == NULL && solver->col_scale == NULL &&
    solver->numeric->Xwork != NULL &&
    solver->i32solve_udiag_recip != NULL &&
    solver->i32solve_udiag_recip_fresh &&
    solver->i16solve_pnum != NULL && solver->i16solve_q != NULL &&
    solver->row_refactor_l_ptr16 != NULL &&
    solver->row_refactor_l_cols16 != NULL &&
    solver->row_refactor_l_row_values != NULL &&
    solver->row_refactor_u_ptr16 != NULL &&
    solver->row_refactor_u_cols16 != NULL &&
    solver->row_refactor_u_row_values != NULL;
}

#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline, hot, aligned(64)))
#endif
static UF_long kls_lean_single_block_prepared_solve(
  kls_solver *solver,
  const double *restrict rhs,
  double *restrict out) {
  const UF_long n = solver->n;
  double *restrict work = (double *)solver->numeric->Xwork;
  const uint16_t *restrict pnum = solver->i16solve_pnum;
  const uint16_t *restrict qperm = solver->i16solve_q;
  const uint16_t *restrict lptr = solver->row_refactor_l_ptr16;
  const uint16_t *restrict lcols = solver->row_refactor_l_cols16;
  const double *restrict lvals = solver->row_refactor_l_row_values;
  const uint16_t *restrict uptr = solver->row_refactor_u_ptr16;
  const uint16_t *restrict ucols = solver->row_refactor_u_cols16;
  const double *restrict uvals = solver->row_refactor_u_row_values;
  const double *restrict udiag_recip = solver->i32solve_udiag_recip;

  UF_long row = solver->i16solve_p_identity_prefix;
  if (row > 0u) {
    memcpy(work, rhs, (size_t)row * sizeof(*work));
  }
  for (; row + 3u < n; row += 4u) {
    work[row] = rhs[(UF_long)pnum[row]];
    work[row + 1u] = rhs[(UF_long)pnum[row + 1u]];
    work[row + 2u] = rhs[(UF_long)pnum[row + 2u]];
    work[row + 3u] = rhs[(UF_long)pnum[row + 3u]];
  }
  for (; row < n; ++row) {
    work[row] = rhs[(UF_long)pnum[row]];
  }
  for (row = 0u; row < n; ++row) {
    work[row] = kls_lean_i16_row_solve_dot(
      lvals, lcols, work, (UF_long)lptr[row],
      (UF_long)lptr[row + 1u], work[row]);
  }
  for (row = n; row-- > 0u;) {
    const double value = kls_lean_i16_row_solve_dot(
      uvals, ucols, work, (UF_long)uptr[row],
      (UF_long)uptr[row + 1u], work[row]);
    work[row] = value * udiag_recip[row];
  }
  row = solver->i16solve_q_identity_prefix;
  if (row > 0u) {
    memcpy(out, work, (size_t)row * sizeof(*out));
  }
  for (; row + 3u < n; row += 4u) {
    out[(UF_long)qperm[row]] = work[row];
    out[(UF_long)qperm[row + 1u]] = work[row + 1u];
    out[(UF_long)qperm[row + 2u]] = work[row + 2u];
    out[(UF_long)qperm[row + 3u]] = work[row + 3u];
  }
  for (; row < n; ++row) {
    out[(UF_long)qperm[row]] = work[row];
  }
  return 1u;
}

/* One-rhs solve with klu_l_solve semantics over the i32 streams. */
/* Keep this recurring sparse kernel's placement stable when unrelated
   diagnostics or setup code change earlier in the translation unit. */
#if defined(__GNUC__) || defined(__clang__)
__attribute__((aligned(64)))
#endif
static UF_long kls_i32_solve(kls_solver *solver,
                             const double *rhs,
                             double *out,
                             const uint16_t *rhs_perm,
                             const double *rhs_scale,
                             const double *out_scale,
                             int rhs_prepared) {
  const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  trilinos_klu_l_numeric *numeric = solver->numeric;
  const UF_long n = symbolic->n;
  const UF_long *Q = symbolic->Q;
  const UF_long *R = symbolic->R;
  const UF_long *pnum = numeric->Pnum;
  const uint32_t *restrict pnum32 = solver->i32solve_pnum;
  const uint32_t *restrict llen32 = solver->i32solve_llen;
  const uint32_t *restrict ulen32 = solver->i32solve_ulen;
  const uint32_t *restrict loff32 = solver->i32solve_loff32;
  const uint32_t *restrict uoff32 = solver->i32solve_uoff32;
  const UF_long *offp = numeric->Offp;
  const UF_long *offi = numeric->Offi;
  const double *offx = (const double *)numeric->Offx;
  const double *udiag = (const double *)numeric->Udiag;
  /* Descriptor construction may allocate this optional mirror before the
     final executor gates settle.  Pointer existence is therefore not a
     value-validity certificate: consume it only after either a serial
     refresh or the exact packed worker has published every reciprocal.
     Generic packed fallbacks retain exact division by Udiag. */
  const double *udiag_recip =
    (solver->i32solve_udiag_recip_fresh ||
     solver->compact_amf_two_block_exact_recip_fresh)
      ? solver->i32solve_udiag_recip : NULL;
  const uint32_t *restrict q32 = solver->i32solve_q;
  const double *rs = numeric->Rs;
  double *X = (double *)numeric->Xwork;
  if (X == NULL) {
    return 0;
  }
  const int i16_ready =
    solver->i16solve_l != NULL && solver->i16solve_u != NULL &&
    solver->i16solve_loff != NULL && solver->i16solve_uoff != NULL &&
    solver->i16solve_pnum != NULL && solver->i16solve_q != NULL &&
    solver->i16solve_r != NULL &&
    solver->i16solve_offp != NULL && solver->i16solve_offi != NULL &&
    solver->i16solve_offcols != NULL &&
    solver->i16solve_offcol_block_ptr != NULL &&
    solver->i16solve_lx != NULL && solver->i16solve_ux != NULL;
  const int compact_match_direct_row = i16_ready &&
    solver->lean_compact_match_row_factor_active &&
    solver->row_refactor_l_ptr != NULL &&
    solver->row_refactor_l_cols != NULL &&
    solver->row_refactor_l_row_values != NULL &&
    solver->row_refactor_u_ptr != NULL &&
    solver->row_refactor_u_cols != NULL &&
    solver->row_refactor_u_row_values != NULL;
  if (rhs_prepared == 1) {
    memcpy(X, rhs, (size_t)n * sizeof(*X));
  } else if (rhs_prepared == 2 &&
             solver->i32solve_rhs_perm32 != NULL && rs == NULL) {
    /* The 32-bit matched permutation selects a public-frame source, while
       an optional fixed row scale is still indexed by the numeric pivot
       row, exactly as in the compact 16-bit path below. */
    const uint32_t *restrict fused_perm = solver->i32solve_rhs_perm32;
    for (UF_long k = 0u; k < n; ++k) {
      const UF_long numeric_row = pnum32 != NULL
        ? (UF_long)pnum32[k] : pnum[k];
      X[k] = rhs[(UF_long)fused_perm[k]] *
        (rhs_scale != NULL ? rhs_scale[numeric_row] : 1.0);
    }
  } else if (i16_ready && rs == NULL) {
    /* A matched rhs_perm addresses the caller's row frame, while an explicit
       matching row scale remains indexed by the permuted internal row.  The
       numeric pivot row therefore selects the scale even though rhs_perm
       selects the source value. */
    const uint16_t *restrict p16 =
      rhs_perm != NULL ? rhs_perm : solver->i16solve_pnum;
    UF_long k = rhs_perm == NULL ? solver->i16solve_p_identity_prefix : 0u;
    if (k > 0u && rhs_scale == NULL) {
      memcpy(X, rhs, (size_t)k * sizeof(*X));
    } else if (k > 0u) {
      for (UF_long i = 0u; i < k; ++i) {
        X[i] = rhs[i] * rhs_scale[i];
      }
    }
    for (; k + 3u < n; k += 4u) {
      const UF_long p0 = (UF_long)p16[k];
      const UF_long p1 = (UF_long)p16[k + 1u];
      const UF_long p2 = (UF_long)p16[k + 2u];
      const UF_long p3 = (UF_long)p16[k + 3u];
      const UF_long s0 = rhs_perm != NULL
        ? (UF_long)solver->i16solve_pnum[k] : p0;
      const UF_long s1 = rhs_perm != NULL
        ? (UF_long)solver->i16solve_pnum[k + 1u] : p1;
      const UF_long s2 = rhs_perm != NULL
        ? (UF_long)solver->i16solve_pnum[k + 2u] : p2;
      const UF_long s3 = rhs_perm != NULL
        ? (UF_long)solver->i16solve_pnum[k + 3u] : p3;
      X[k] = rhs[p0] * (rhs_scale != NULL ? rhs_scale[s0] : 1.0);
      X[k + 1u] = rhs[p1] * (rhs_scale != NULL ? rhs_scale[s1] : 1.0);
      X[k + 2u] = rhs[p2] * (rhs_scale != NULL ? rhs_scale[s2] : 1.0);
      X[k + 3u] = rhs[p3] * (rhs_scale != NULL ? rhs_scale[s3] : 1.0);
    }
    for (; k < n; ++k) {
      const UF_long p = (UF_long)p16[k];
      const UF_long scale_row = rhs_perm != NULL
        ? (UF_long)solver->i16solve_pnum[k] : p;
      X[k] = rhs[p] *
        (rhs_scale != NULL ? rhs_scale[scale_row] : 1.0);
    }
  } else if (i16_ready) {
    const uint16_t *restrict p16 = solver->i16solve_pnum;
    for (UF_long k = 0; k < n; ++k) {
      const UF_long p = (UF_long)p16[k];
      X[k] = rhs[p] * (rhs_scale != NULL ? rhs_scale[p] : 1.0) / rs[k];
    }
  } else if (rs == NULL && pnum32 != NULL && rhs_perm == NULL) {
    if (rhs_scale != NULL ||
        !kls_parallel_i32_permute(solver, rhs, X, pnum32, 0)) {
      for (UF_long k = 0; k < n; ++k) {
        const UF_long p = (UF_long)pnum32[k];
        X[k] = rhs[p] * (rhs_scale != NULL ? rhs_scale[p] : 1.0);
      }
    }
  } else if (rs == NULL) {
    for (UF_long k = 0; k < n; ++k) {
      const UF_long p = rhs_perm != NULL ? (UF_long)rhs_perm[k] : pnum[k];
      X[k] = rhs[p] * (rhs_scale != NULL ? rhs_scale[p] : 1.0);
    }
  } else {
    for (UF_long k = 0; k < n; ++k) {
      const UF_long p = pnum[k];
      X[k] = rhs[p] * (rhs_scale != NULL ? rhs_scale[p] : 1.0) / rs[k];
    }
  }
  for (UF_long block = symbolic->nblocks; block-- > 0;) {
    const UF_long k1 = i16_ready
      ? (UF_long)solver->i16solve_r[block] : R[block];
    const UF_long k2 = i16_ready
      ? (UF_long)solver->i16solve_r[block + 1u] : R[block + 1u];
    const UF_long nk = k2 - k1;
    if (i16_ready && solver->i16solve_singleton_run != NULL &&
        solver->i16solve_singleton_run[block] != 0u) {
      const UF_long run =
        (UF_long)solver->i16solve_singleton_run[block];
      const UF_long first_block = block + 1u - run;
      const UF_long first_col =
        (UF_long)solver->i16solve_r[first_block];
      const uint16_t *restrict offp16 = solver->i16solve_offp;
      const uint16_t *restrict offi16 = solver->i16solve_offi;
      for (UF_long k = k2; k > first_col;) {
        --k;
        const double xk = kls_i32_solve_diagonal(
          X[k], udiag, udiag_recip, k);
        X[k] = xk;
        if (first_block != 0u || k != first_col) {
          for (UF_long p = (UF_long)offp16[k];
               p < (UF_long)offp16[k + 1u]; ++p) {
            X[offi16[p]] -= offx[p] * xk;
          }
        }
      }
      block = first_block;
      continue;
    }
    if (!i16_ready && solver->i32solve_singleton_run != NULL &&
        solver->i32solve_singleton_run[block] != 0u) {
      const UF_long run =
        (UF_long)solver->i32solve_singleton_run[block];
      const UF_long first_block = block + 1u - run;
      const UF_long first_col = R[first_block];
      for (UF_long k = k2; k > first_col;) {
        --k;
        const double xk = kls_i32_solve_diagonal(
          X[k], udiag, udiag_recip, k);
        X[k] = xk;
        if (first_block != 0u || k != first_col) {
          for (UF_long p = offp[k]; p < offp[k + 1u]; ++p) {
            X[offi[p]] -= offx[p] * xk;
          }
        }
      }
      block = first_block;
      continue;
    }
    if (i16_ready && nk == 1u) {
      const double xk = kls_i32_solve_diagonal(
        X[k1], udiag, udiag_recip, k1);
      X[k1] = xk;
      if (block > 0u) {
        const uint16_t *restrict offp16 = solver->i16solve_offp;
        const uint16_t *restrict offi16 = solver->i16solve_offi;
        for (UF_long p = (UF_long)offp16[k1];
             p < (UF_long)offp16[k1 + 1u]; ++p) {
          X[offi16[p]] -= offx[p] * xk;
        }
      }
      continue;
    }
    if (nk == 1u) {
      X[k1] = kls_i32_solve_diagonal(
        X[k1], udiag, udiag_recip, k1);
    } else {
      double *Xb = X + k1;
      if (compact_match_direct_row) {
        const UF_long *restrict lptr = solver->row_refactor_l_ptr;
        const uint16_t *restrict lptr16 =
          solver->row_refactor_l_ptr16;
        const UF_long *restrict lcols = solver->row_refactor_l_cols;
        const uint16_t *restrict lcols16 =
          solver->row_refactor_l_cols16;
        const double *restrict lvals = solver->row_refactor_l_row_values;
        const UF_long *restrict uptr = solver->row_refactor_u_ptr;
        const uint16_t *restrict uptr16 =
          solver->row_refactor_u_ptr16;
        const UF_long *restrict ucols = solver->row_refactor_u_cols;
        const uint16_t *restrict ucols16 =
          solver->row_refactor_u_cols16;
        const double *restrict uvals = solver->row_refactor_u_row_values;
        if (lcols16 != NULL && ucols16 != NULL &&
            lptr16 != NULL && uptr16 != NULL) {
          for (UF_long row = k1; row < k2; ++row) {
            X[row] = kls_lean_i16_row_solve_dot(
              lvals, lcols16, X,
              (UF_long)lptr16[row], (UF_long)lptr16[row + 1u], X[row]);
          }
          for (UF_long row = k2; row-- > k1;) {
            const double value = kls_lean_i16_row_solve_dot(
              uvals, ucols16, X,
              (UF_long)uptr16[row], (UF_long)uptr16[row + 1u], X[row]);
            X[row] = kls_i32_solve_diagonal(
              value, udiag, udiag_recip, row);
          }
        } else if (lcols16 != NULL && ucols16 != NULL) {
          for (UF_long row = k1; row < k2; ++row) {
            X[row] = kls_lean_i16_row_solve_dot(
              lvals, lcols16, X, lptr[row], lptr[row + 1u], X[row]);
          }
          for (UF_long row = k2; row-- > k1;) {
            const double value = kls_lean_i16_row_solve_dot(
              uvals, ucols16, X, uptr[row], uptr[row + 1u], X[row]);
            X[row] = kls_i32_solve_diagonal(
              value, udiag, udiag_recip, row);
          }
        } else {
          for (UF_long row = k1; row < k2; ++row) {
            X[row] = kls_compact_match_row_solve_dot(
              lvals, lcols, X, lptr[row], lptr[row + 1u], X[row]);
          }
          for (UF_long row = k2; row-- > k1;) {
            const double value = kls_compact_match_row_solve_dot(
              uvals, ucols, X, uptr[row], uptr[row + 1u], X[row]);
            X[row] = kls_i32_solve_diagonal(
              value, udiag, udiag_recip, row);
          }
        }
      } else if (i16_ready) {
        if (!solver->sparse_refinement_rhs_active) {
          for (UF_long k = 0u; k < nk; ++k) {
            const double xk = Xb[k];
            const UF_long global = k1 + k;
            const UF_long begin = (UF_long)solver->i16solve_loff[global];
            const UF_long end = (UF_long)solver->i16solve_loff[global + 1u];
            kls_i16_solve_scatter(Xb, solver->i16solve_l + begin,
                                  solver->i16solve_lx[global], end - begin,
                                  xk);
          }
          for (UF_long k = nk; k-- > 0u;) {
            const UF_long global = k1 + k;
            const double xk = kls_i32_solve_diagonal(
              Xb[k], udiag, udiag_recip, global);
            Xb[k] = xk;
            const UF_long begin = (UF_long)solver->i16solve_uoff[global];
            const UF_long end = (UF_long)solver->i16solve_uoff[global + 1u];
            kls_i16_solve_scatter(Xb, solver->i16solve_u + begin,
                                  solver->i16solve_ux[global], end - begin,
                                  xk);
          }
        } else {
          for (UF_long k = 0u; k < nk; ++k) {
            const double xk = Xb[k];
            const UF_long global = k1 + k;
            if (xk != 0.0) {
              const UF_long begin =
                (UF_long)solver->i16solve_loff[global];
              const UF_long end =
                (UF_long)solver->i16solve_loff[global + 1u];
              kls_i16_solve_scatter(Xb, solver->i16solve_l + begin,
                                    solver->i16solve_lx[global], end - begin,
                                    xk);
            }
          }
          for (UF_long k = nk; k-- > 0u;) {
            const UF_long global = k1 + k;
            const double xk = kls_i32_solve_diagonal(
              Xb[k], udiag, udiag_recip, global);
            Xb[k] = xk;
            if (xk != 0.0) {
              const UF_long begin =
                (UF_long)solver->i16solve_uoff[global];
              const UF_long end =
                (UF_long)solver->i16solve_uoff[global + 1u];
              kls_i16_solve_scatter(Xb, solver->i16solve_u + begin,
                                    solver->i16solve_ux[global], end - begin,
                                    xk);
            }
          }
        }
      } else {
        const UF_long *lip = numeric->Lip + k1;
        const UF_long *llen = numeric->Llen + k1;
        const UF_long *uip = numeric->Uip + k1;
        const UF_long *ulen = numeric->Ulen + k1;
        double *lu = (double *)numeric->LUbx[block];
        kls_pts *const pts = solver->pts;
        const int pts_candidate =
          pts != NULL && pts->block == block &&
          (pts->nk >= 16384u ||
           kls_repeated_update_workload(&solver->options)) &&
          pts->solve_ok;
        int use_pts = 0;
        int pts_used = 0;
        double pts_probe_start = 0.0;
        if (pts_candidate) {
          if (pts->solve_decision != 0) {
            use_pts = pts->solve_decision > 0;
          } else {
            /* Order the four samples serial,PTS,PTS,serial so neither arm
               owns both the cold-first nor warm-last position. */
            const int probe =
              pts->solve_probe_samples[0] + pts->solve_probe_samples[1];
            use_pts = probe == 1 || probe == 2;
          }
          pts_probe_start = kls_now_seconds();
          if (use_pts) {
            pts_used = kls_pts_solve_block(
              solver, Xb, lu, lip, llen, uip, ulen, udiag + k1,
              udiag_recip != NULL ? udiag_recip + k1 : NULL);
            if (!pts_used) {
              pts->solve_decision = -1;
              use_pts = 0;
            }
          }
        }
        if (!pts_used) {
          /* L-sweep, unit diagonal, block-local indices */
          for (UF_long k = 0; k < nk; ++k) {
            const double xk = Xb[k];
            if (xk != 0.0) {
              const UF_long len = llen32 != NULL
                ? (UF_long)llen32[k1 + k] : llen[k];
              const double *lx = solver->refactor_l_packed_valid &&
                  solver->refactor_l_packed_values != NULL &&
                  solver->refactor_l_packed_values[k1 + k] != NULL
                ? solver->refactor_l_packed_values[k1 + k]
                : lu + lip[k] + kls_klu_units_for_indices(len);
              const UF_long begin = loff32 != NULL
                ? (UF_long)loff32[k1 + k]
                : (UF_long)solver->i32solve_loff[k1 + k];
              const int32_t *li = solver->i32solve_l + begin;
              for (UF_long p = 0; p < len; ++p) {
                Xb[li[p]] -= lx[p] * xk;
              }
            }
          }
          /* U-sweep, descending, dividing by the diagonal */
          for (UF_long k = nk; k-- > 0;) {
            const double xk = kls_i32_solve_diagonal(
              Xb[k], udiag, udiag_recip, k1 + k);
            Xb[k] = xk;
            if (xk != 0.0) {
              const UF_long len = ulen32 != NULL
                ? (UF_long)ulen32[k1 + k] : ulen[k];
              const double *ux = solver->refactor_u_packed_valid &&
                  solver->refactor_u_packed_values != NULL &&
                  solver->refactor_u_packed_values[k1 + k] != NULL
                ? solver->refactor_u_packed_values[k1 + k]
                : lu + uip[k] + kls_klu_units_for_indices(len);
              const UF_long begin = uoff32 != NULL
                ? (UF_long)uoff32[k1 + k]
                : (UF_long)solver->i32solve_uoff[k1 + k];
              const int32_t *ui = solver->i32solve_u + begin;
              for (UF_long p = 0; p < len; ++p) {
                Xb[ui[p]] -= ux[p] * xk;
              }
            }
          }
        }
        if (pts_candidate && pts->solve_decision == 0) {
          const int arm = pts_used ? 1 : 0;
          const double elapsed = kls_now_seconds() - pts_probe_start;
          pts->solve_probe_samples[arm]++;
          if (pts->solve_probe_min[arm] <= 0.0 ||
              elapsed < pts->solve_probe_min[arm]) {
            pts->solve_probe_min[arm] = elapsed;
          }
          if (pts->solve_probe_samples[0] > 0 &&
              pts->solve_probe_samples[1] > 0 &&
              pts->solve_probe_min[1] > 1.25 * pts->solve_probe_min[0]) {
            pts->solve_decision = -1;
          } else if (pts->solve_probe_samples[0] >= 2 &&
                     pts->solve_probe_samples[1] >= 2) {
            pts->solve_decision =
              pts->solve_probe_min[1] < 0.98 * pts->solve_probe_min[0]
                ? 1 : -1;
          }
        }
      }
    }
    if (block > 0) {
      if (i16_ready) {
        const uint16_t *restrict offp16 = solver->i16solve_offp;
        const uint16_t *restrict offi16 = solver->i16solve_offi;
        const uint16_t *restrict offcols = solver->i16solve_offcols;
        const UF_long begin =
          (UF_long)solver->i16solve_offcol_block_ptr[block];
        const UF_long end =
          (UF_long)solver->i16solve_offcol_block_ptr[block + 1u];
        for (UF_long pos = begin; pos < end; ++pos) {
          const UF_long k = (UF_long)offcols[pos];
          const double xk = X[k];
          for (UF_long p = (UF_long)offp16[k];
               p < (UF_long)offp16[k + 1u]; ++p) {
            X[offi16[p]] -= offx[p] * xk;
          }
        }
      } else {
        for (UF_long k = k1; k < k2; ++k) {
          const double xk = X[k];
          if (xk == 0.0) continue;
          const UF_long pend = offp[k + 1u];
          for (UF_long p = offp[k]; p < pend; ++p) {
            X[offi[p]] -= offx[p] * xk;
          }
        }
      }
    }
  }
  if (i16_ready) {
    const uint16_t *restrict q16 = solver->i16solve_q;
    UF_long k = solver->i16solve_q_identity_prefix;
    if (k > 0u && out_scale == NULL) {
      memcpy(out, X, (size_t)k * sizeof(*out));
    } else if (k > 0u) {
      for (UF_long i = 0u; i < k; ++i) {
        out[i] = X[i] * out_scale[i];
      }
    }
    for (; k + 3u < n; k += 4u) {
      const UF_long q0 = (UF_long)q16[k];
      const UF_long q1 = (UF_long)q16[k + 1u];
      const UF_long q2 = (UF_long)q16[k + 2u];
      const UF_long q3 = (UF_long)q16[k + 3u];
      out[q0] = X[k] * (out_scale != NULL ? out_scale[q0] : 1.0);
      out[q1] = X[k + 1u] * (out_scale != NULL ? out_scale[q1] : 1.0);
      out[q2] = X[k + 2u] * (out_scale != NULL ? out_scale[q2] : 1.0);
      out[q3] = X[k + 3u] * (out_scale != NULL ? out_scale[q3] : 1.0);
    }
    for (; k < n; ++k) {
      const UF_long q = (UF_long)q16[k];
      out[q] = X[k] * (out_scale != NULL ? out_scale[q] : 1.0);
    }
  } else if (q32 != NULL) {
    if (out_scale != NULL ||
        !kls_parallel_i32_permute(solver, X, out, q32, 1)) {
      for (UF_long k = 0; k < n; ++k) {
        const UF_long q = (UF_long)q32[k];
        out[q] = X[k] * (out_scale != NULL ? out_scale[q] : 1.0);
      }
    }
  } else {
    for (UF_long k = 0; k < n; ++k) {
      const UF_long q = Q[k];
      out[q] = X[k] * (out_scale != NULL ? out_scale[q] : 1.0);
    }
  }
  return 1;
}

static void kls_run_deferred_factor_preps(kls_solver *solver,
                                           const double *numeric_values);

/* KLU's threshold pivoting prefers the symbolic diagonal when it is at least
   tol times the largest active-column entry.  Model max-row scaling in the
   actual symbolic BTF frame and count crossings of that same tolerance.  A
   positive net crossing is a numerical reason to scale before factorization,
   independent of matrix identity, dimension, or a fitted work threshold. */
static int kls_serial_scaled_pivots_preferred(const kls_solver *solver,
                                              const double *numeric_values) {
  if (solver == NULL || solver->symbolic == NULL || solver->col_ptr == NULL ||
      solver->row_idx == NULL || numeric_values == NULL ||
      solver->options.pivot_tolerance <= 0.0) {
    return 0;
  }
  const UF_long n = solver->n;
  const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  double *row_max = (double *)calloc((size_t)n, sizeof(*row_max));
  UF_long *pinv = (UF_long *)malloc((size_t)n * sizeof(*pinv));
  if (row_max == NULL || pinv == NULL) {
    free(row_max);
    free(pinv);
    return 0;
  }
  for (UF_long k = 0; k < n; ++k) {
    pinv[symbolic->P[k]] = k;
  }
  for (UF_long col = 0; col < n; ++col) {
    for (UF_long p = solver->col_ptr[col]; p < solver->col_ptr[col + 1u];
         ++p) {
      const UF_long row = solver->row_idx[p];
      const double a = fabs(numeric_values[p]);
      if (isfinite(a) && a > row_max[row]) {
        row_max[row] = a;
      }
    }
  }

  UF_long improved = 0;
  UF_long harmed = 0;
  for (UF_long block = 0; block < symbolic->nblocks; ++block) {
    const UF_long k1 = symbolic->R[block];
    const UF_long k2 = symbolic->R[block + 1u];
    if (k2 - k1 <= 1u) {
      continue;
    }
    for (UF_long k = k1; k < k2; ++k) {
      const UF_long oldcol = symbolic->Q[k];
      const UF_long olddiagrow = symbolic->P[k];
      double diag = 0.0;
      double raw_col_max = 0.0;
      double scaled_col_max = 0.0;
      for (UF_long p = solver->col_ptr[oldcol];
           p < solver->col_ptr[oldcol + 1u]; ++p) {
        const UF_long oldrow = solver->row_idx[p];
        const UF_long newrow = pinv[oldrow];
        if (newrow < k1 || newrow >= k2) {
          continue;
        }
        const double a = fabs(numeric_values[p]);
        if (a > raw_col_max) {
          raw_col_max = a;
        }
        if (row_max[oldrow] > 0.0) {
          const double scaled = a / row_max[oldrow];
          if (scaled > scaled_col_max) {
            scaled_col_max = scaled;
          }
        }
        if (oldrow == olddiagrow && a > diag) {
          diag = a;
        }
      }
      const double raw_ratio = raw_col_max > 0.0 ? diag / raw_col_max : 0.0;
      const double scaled_diag = row_max[olddiagrow] > 0.0
        ? diag / row_max[olddiagrow] : 0.0;
      const double scaled_ratio = scaled_col_max > 0.0
        ? scaled_diag / scaled_col_max : 0.0;
      const int raw_weak = raw_ratio < solver->options.pivot_tolerance;
      const int scaled_weak = scaled_ratio < solver->options.pivot_tolerance;
      improved += (UF_long)(raw_weak && !scaled_weak);
      harmed += (UF_long)(!raw_weak && scaled_weak);
    }
  }
  free(row_max);
  free(pinv);
  return improved > harmed;
}

static int kls_serial_factor(kls_solver *solver,
                             double *numeric_values) {
  if (solver == NULL || solver->symbolic == NULL || numeric_values == NULL) {
    return KLS_ERR_INVALID_ARGUMENT;
  }

  const double start = kls_now_seconds();
  free_numeric(solver);
  solver->factor_preps_deferred = 0;
  int auto_scale_policy = solver->options.scale;
  int auto_scale_policy_known =
    solver->options.scale != KLS_SCALE_AUTO;
  const int scaled_pivots_preferred =
    solver->options.scale == KLS_SCALE_AUTO &&
    kls_serial_scaled_pivots_preferred(solver, numeric_values);
  if (scaled_pivots_preferred) {
    auto_scale_policy = choose_auto_scale_from_values(solver, numeric_values);
    auto_scale_policy_known = 1;
  }
  const int preselect_max_scale =
    scaled_pivots_preferred && auto_scale_policy == 2;
  solver->common.scale = solver->options.scale == KLS_SCALE_AUTO
    ? (preselect_max_scale ? 2 : -1) : solver->options.scale;
  solver->common.tol = solver->options.pivot_tolerance;
  solver->common.status = TRILINOS_KLU_OK;
  solver->common.numerical_rank = KLS_KLU_EMPTY;
  solver->common.singular_col = KLS_KLU_EMPTY;
  solver->common.kls_dense_panels = 0;
  kls_klu_pipe_threads = 0;
  kls_klu_pipe_det = 0;
  kls_klu_dense_tail = 0;
  solver->dense_tail_cols = 0;
  solver->dense_tail_block = 0;
  solver->numeric = trilinos_klu_l_factor(solver->col_ptr,
                                          solver->row_idx,
                                          numeric_values,
                                          solver->symbolic,
                                          &solver->common);

  int diagnostics_have_flops = 0;
  int unscaled_work_inflated = 0;
  if (solver->options.scale == KLS_SCALE_AUTO && solver->common.scale != 2 &&
      solver->numeric != NULL &&
      solver->common.status >= 0 && solver->symbolic->est_flops > 0.0) {
    diagnostics_have_flops =
      trilinos_klu_l_flops(solver->symbolic, solver->numeric,
                           &solver->common) != 0;
    unscaled_work_inflated = diagnostics_have_flops &&
      solver->common.flops > 1.50 * solver->symbolic->est_flops;
  }

  /* A failed factor always gets KLU's robust max-row-scaled retry.  Do not
     throw away a successful unscaled numeric on the flop-estimate warning:
     inflated scalar work can still expose a much faster mapped Horizon
     cycle.  Require AUTO's independent value/structure policy to agree that
     max-row scaling is appropriate before paying for and adopting a
     replacement. */
  if (unscaled_work_inflated && !auto_scale_policy_known) {
    auto_scale_policy = choose_auto_scale_from_values(solver, numeric_values);
    auto_scale_policy_known = 1;
  }
  const int scale_policy_confirms_retry =
    unscaled_work_inflated && auto_scale_policy_known &&
    auto_scale_policy == 2;
  if (solver->options.scale == KLS_SCALE_AUTO &&
      solver->common.scale != 2 &&
      (solver->numeric == NULL || solver->common.status < 0 ||
       solver->common.status == TRILINOS_KLU_SINGULAR ||
       scale_policy_confirms_retry)) {
    free_numeric(solver);
    solver->common.scale = 2;
    solver->common.status = TRILINOS_KLU_OK;
    solver->common.numerical_rank = KLS_KLU_EMPTY;
    solver->common.singular_col = KLS_KLU_EMPTY;
    diagnostics_have_flops = 0;
    solver->numeric = trilinos_klu_l_factor(solver->col_ptr,
                                            solver->row_idx,
                                            numeric_values,
                                            solver->symbolic,
                                            &solver->common);
  }

  solver->stats.factor_seconds = kls_now_seconds() - start;
  kls_set_last_factor_path(solver, KLS_FACTOR_PATH_SERIAL);
  if (solver->numeric == NULL || solver->common.status < 0) {
    fill_numeric_stats(solver);
    return solver->common.status == TRILINOS_KLU_SINGULAR
      ? KLS_ERR_SINGULAR : KLS_ERR_FACTOR_FAILED;
  }
  if (diagnostics_have_flops) {
    /* The inflation guard already computed the retained factor's flop count;
       only reciprocal condition estimation remains. */
    kls_update_numeric_rcond(solver);
  } else {
    kls_update_numeric_diagnostics(solver, 1);
  }
  if (solver->n >= 512) {
    /* Preserve the lean one-shot, then let the first refactor pay once for
       the same mapped-pointer/index32 preparations used by native KLS. */
    solver->factor_preps_deferred = 1;
  }
  solver->stats.factor_seconds = kls_now_seconds() - start;
  fill_numeric_stats(solver);
  return solver->common.status == TRILINOS_KLU_SINGULAR
    ? KLS_ERR_SINGULAR : KLS_OK;
}

/* Classify the current numeric's solve-accuracy risk once, from pivot
   growth (one A+U pass over the internal-frame values, charged to the
   refactorization that ran it, never to the solves).  Plain-value factors,
   including pure matched row permutations, use the same permutation-
   invariant growth certificate.  Scaled, nudged, perturbed, predicted and
   reduced-precision numerics carry structural risk flags and probe their
   first solve instead unless a narrow raw-accuracy certificate covers the
   retained refactor path.  The Udiag-spread scan measured OUT as a signal:
   healthy rajat25 spans 3.1e-9 while the failing b2383 spans 7.2e-9 — no
   threshold separates them; reciprocal pivot growth does (b2383 reads
   1.3e-7). */
static void kls_solve_contract_classify(kls_solver *solver,
                                        const double *numeric_values) {
  if (solver == NULL || solver->solve_contract_probe != 0 ||
      solver->row_scale != NULL ||
      solver->col_scale != NULL || solver->pivot_nudge_count > 0 ||
      solver->common.kls_perturb_count > 0 ||
      solver->numeric_is_predicted ||
      solver->numeric == NULL || solver->symbolic == NULL ||
      numeric_values == NULL || solver->col_ptr == NULL ||
      solver->row_idx == NULL) {
    return;
  }
  kls_update_numeric_rcond_guard(solver);
  if (solver->options.backend != KLS_BACKEND_SERIAL &&
      solver->common.rcond > 0.0 &&
      solver->common.rcond < sqrt(DBL_EPSILON)) {
    /* Reciprocal pivot growth measures backward stability of the factor,
       but cannot certify the forward effect of an almost singular pivot
       sequence.  The inexpensive current-U diagonal ratio supplies that
       missing numeric condition and takes precedence over profile bypasses
       or a healthy rgrowth verdict. */
    solver->solve_contract_probe =
      solver->low_rcond_solve_contract_state == 1 ? 1 : 2;
    solver->solve_contract_verified =
      solver->low_rcond_solve_contract_state == 1;
    return;
  }

  if (solver->options.backend == KLS_BACKEND_SERIAL &&
      solver->stats.last_refactor_path == KLS_REFACTOR_PATH_KLU) {
    /* Explicit serial KLU-refactor mode keeps base KLU solve semantics.
       The factor and refactor are the same full-precision KLU operations;
       arming KLS's adaptive residual correction here turns every solve into
       an SpMV plus two triangular solves.  AUTO/KLS paths retain that extra
       contract protection, including every multi-thread execution. */
    solver->solve_contract_probe = 1;
    solver->solve_contract_verified = 1;
    return;
  }
  trilinos_klu_l_common growth_common = solver->common;
  if (!trilinos_klu_l_rgrowth(
        solver->col_ptr, solver->row_idx,
        (double *)numeric_values, solver->symbolic,
        solver->numeric, &growth_common)) {
    return;
  }
  const double thresh = 1.0e-6;
  solver->solve_contract_probe =
    growth_common.rgrowth < thresh && growth_common.rgrowth >= 0.0 ? 2 : 1;
  if (solver->solve_contract_probe == 2) {
    /* growth-armed corrections take the single-shot exit from the first
       solve on: every solve still measures the raw residual (the SpMV
       the correction is computed from), and a correction this class
       cannot repair shows up in the caller's residual check rather than
       costing every solve a verification sweep it always passes */
    solver->solve_contract_verified = 1;
  }
}

/* A public factor call starts a new numeric epoch just as a public refactor
   does, but it may reach that epoch through an in-place checked update or a
   complete replacement.  Classify every successful plain-frame factor and
   retain the values needed by an armed first solve.  Without this factor-side
   capture, a low-rcond replacement either refined against the preceding
   numeric's matrix or had no current matrix with which to verify the solve. */
static void kls_factor_solve_contract_classify(
  kls_solver *solver,
  const double *numeric_values) {
  if (solver == NULL || numeric_values == NULL ||
      solver->options.backend == KLS_BACKEND_SERIAL) {
    return;
  }
  kls_solve_contract_classify(solver, numeric_values);
  if (solver->solve_contract_probe != 2 || solver->nnz == 0u ||
      solver->nnz > (UF_long)(KLS_MAX_ALLOCATION / sizeof(double)) ||
      solver->row_perm != NULL || solver->row_scale != NULL ||
      solver->col_scale != NULL) {
    return;
  }
  if (solver->solve_refine_values == NULL) {
    solver->solve_refine_values = (double *)malloc(
      (size_t)solver->nnz * sizeof(*solver->solve_refine_values));
  }
  if (solver->solve_refine_values != NULL &&
      solver->solve_refine_values != numeric_values) {
    memcpy(solver->solve_refine_values, numeric_values,
           (size_t)solver->nnz * sizeof(*numeric_values));
  }
}

static void *kls_lean_prewarm_main(void *arg) {
  kls_lean_prewarm_job *job = (kls_lean_prewarm_job *)arg;
  if (job == NULL || job->solver == NULL || job->thread_count < 2) {
    return NULL;
  }
  if (job->build_pattern == 1) {
    (void)kls_build_row_refactor_pattern(job->solver, 1);
    /* These narrow descriptors are a second immutable view of the row
       pattern just published by this worker.  Their builder validates all
       representational bounds, so attempting it is safe for both the
       packed executor and the broader moderate-work row executor; a
       factor outside the narrow bounds simply retains native indices. */
    (void)kls_build_lean_row_i16_indices(job->solver);
    return NULL;
  }
  if (job->build_pattern == 2) {
    (void)kls_i32_solve_ready(job->solver);
    return NULL;
  }
  /* These structures depend only on n and the requested concurrency, so
     their allocation and worker startup can overlap the independent serial
     KLU first factor.  Numeric-dependent row mirrors are still built after
     the factor and remain fully charged. */
  (void)ensure_egraph_refactor_pool(job->solver, job->thread_count);
  (void)ensure_egraph_worker_scratch(job->solver, job->thread_count,
                                     job->solver->n);
  unsigned int generation = 0u;
  (void)ensure_lean_parallel_done(job->solver, &generation);
  return NULL;
}
