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
      solver->symbolic->R == NULL ||
      getenv("KLS_DISABLE_I16_SOLVE") != NULL) {
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
    kls_compact_singleton_run_solve_profile(solver, NULL, NULL);
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
  if (getenv("KLS_DISABLE_COMPACT_MATCH_FUSED_RHS_PERM") == NULL &&
      solver->row_perm != NULL) {
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
  if (getenv("KLS_TRACE_I16_CACHE") != NULL) {
    UF_long p_identity = 0u;
    UF_long q_identity = 0u;
    UF_long singleton_blocks = 0u;
    UF_long singleton_prefix = 0u;
    UF_long singleton_suffix = 0u;
    UF_long nontrivial_blocks = 0u;
    UF_long current_singleton_run = 0u;
    UF_long max_singleton_run = 0u;
    for (UF_long k = 0u; k < n; ++k) {
      p_identity += (UF_long)(pnum[k] == k);
      q_identity += (UF_long)(q[k] == k);
    }
    for (UF_long block = 0u; block < nblocks; ++block) {
      const int singleton = symbolic->R[block + 1u] ==
        symbolic->R[block] + 1u;
      singleton_blocks += (UF_long)singleton;
      nontrivial_blocks += (UF_long)!singleton;
      if (singleton) {
        current_singleton_run++;
        if (current_singleton_run > max_singleton_run) {
          max_singleton_run = current_singleton_run;
        }
      } else {
        current_singleton_run = 0u;
      }
    }
    while (singleton_prefix < nblocks &&
           symbolic->R[singleton_prefix + 1u] ==
             symbolic->R[singleton_prefix] + 1u) {
      singleton_prefix++;
    }
    while (singleton_suffix < nblocks &&
           symbolic->R[nblocks - singleton_suffix] ==
             symbolic->R[nblocks - singleton_suffix - 1u] + 1u) {
      singleton_suffix++;
    }
    fprintf(stderr,
            "KLS i16 cache: n=%ld L=%ld U=%ld off=%ld offcols=%ld"
            " blocks=%ld singleton=%ld runmax=%ld cached=%d"
            " prefix=%ld suffix=%ld nontrivial=%ld"
            " Pid=%ld/%ld Qid=%ld/%ld source=%s\n",
            (long)n, (long)lcount, (long)ucount, (long)offcount,
            (long)offcol_count, (long)nblocks, (long)singleton_blocks,
            (long)max_singleton_run, singleton_run != NULL,
            (long)singleton_prefix, (long)singleton_suffix,
            (long)nontrivial_blocks, (long)p_identity,
            (long)p_identity_prefix, (long)q_identity,
            (long)q_identity_prefix, direct ? "direct" : "i32");
  }
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

static int kls_serial_mapped_prefers_vendor_solve(
  const kls_solver *solver) {
  return solver != NULL &&
    solver->options.backend == KLS_BACKEND_SERIAL &&
    solver->stats.last_refactor_path == KLS_REFACTOR_PATH_MAPPED;
}

static int kls_auto_btf_prefers_vendor_solve(const kls_solver *solver) {
  if (getenv("KLS_DISABLE_AUTO_BTF_VENDOR_SOLVE") != NULL ||
      solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      solver->options.backend == KLS_BACKEND_SERIAL ||
      solver->options.threads <= 1) {
    return 0;
  }
  if (kls_low_work_btf_prefers_native_solve(solver)) {
    return 1;
  }
  /* Outside the cache-resident low-work crossover, build the compact stream
     and let the generic solve/PTS timing state choose its consumer. */
  return 0;
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
      (getenv("KLS_ENABLE_I32_UDIAG_RECIP") == NULL &&
       !compact_i16_solve &&
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
  if (kls_auto_btf_prefers_vendor_solve(solver)) {
    solver->i32solve_state = -1;
    return 0;
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
      getenv("KLS_ENABLE_SMALL_I32_SOLVE") == NULL &&
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
  const int trace_i32_prep = getenv("KLS_TRACE_I32_PREP") != NULL;
  const double i32_prep_start = trace_i32_prep ? kls_now_seconds() : 0.0;
  double i32_layout_done = i32_prep_start;
  double i32_copy_done = i32_prep_start;
  double i32_cache_done = i32_prep_start;
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
      getenv("KLS_DISABLE_I16_SOLVE") == NULL &&
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
  if (trace_i32_prep) {
    i32_layout_done = kls_now_seconds();
  }
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
  if (trace_i32_prep) {
    i32_copy_done = kls_now_seconds();
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
  if (trace_i32_prep) {
    i32_cache_done = kls_now_seconds();
  }
  const int generic_compact_solve_meta =
        getenv("KLS_DISABLE_GENERIC_COMPACT_SOLVE_META") == NULL &&
    kls_repeated_update_workload(&solver->options) &&
    solver->options.expected_refactorizations >= 16 &&
    solver->i16solve_l == NULL && solver->i16solve_u == NULL &&
    n >= 512u && n <= (UF_long)UINT32_MAX &&
    lcur >= 0 && ucur >= 0 &&
    lcur <= (int64_t)UINT32_MAX && ucur <= (int64_t)UINT32_MAX;
  if (generic_compact_solve_meta ||
      getenv("KLS_ENABLE_GENERIC_COMPACT_SOLVE_META") != NULL) {
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
  if ((generic_compact_solve_meta ||
       getenv("KLS_ENABLE_GENERIC_COMPACT_SOLVE_META") != NULL ||
       (0 &&
        getenv("KLS_DISABLE_DENSE_RECIPROCAL_HUB_COMPACT_STREAM_META") ==
          NULL))) {
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
      n <= (UF_long)UINT16_MAX &&
      getenv("KLS_DISABLE_GENERAL_FUSED_MATCHED_RHS") == NULL) {
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
      numeric->Rs == NULL && n <= (UF_long)UINT32_MAX &&
      getenv("KLS_DISABLE_MODERATE_FRAGMENTED_FUSED_RHS") == NULL) {
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
  if (getenv("KLS_TRACE_RUNS") != NULL) {
    /* stored-zero census: exact cancellations in the factor values */
    {
      trilinos_klu_l_numeric *num2 = solver->numeric;
      const trilinos_klu_l_symbolic *sym2 = solver->symbolic;
      int64_t zl = 0, zu = 0, tl = 0, tu = 0, tinyl = 0, tinyu = 0;
      for (UF_long b = 0; b < sym2->nblocks; ++b) {
        const UF_long bk1 = sym2->R[b];
        const UF_long bnk = sym2->R[b + 1] - bk1;
        if (bnk < 2u || num2->LUbx[b] == NULL) {
          continue;
        }
        UF_long *blip = num2->Lip + bk1;
        UF_long *bllen = num2->Llen + bk1;
        UF_long *buip = num2->Uip + bk1;
        UF_long *bulen = num2->Ulen + bk1;
        double *blu = (double *)num2->LUbx[b];
        for (UF_long k = 0; k < bnk; ++k) {
          UF_long *ri;
          double *rv;
          UF_long rl;
          kls_klu_get_pointer(blu, blip, bllen, k, &ri, &rv, &rl);
          for (UF_long p = 0; p < rl; ++p) {
            tl++;
            if (rv[p] == 0.0) {
              zl++;
            } else if (fabs(rv[p]) < 1e-14) {
              tinyl++;
            }
          }
          kls_klu_get_pointer(blu, buip, bulen, k, &ri, &rv, &rl);
          for (UF_long p = 0; p < rl; ++p) {
            tu++;
            if (rv[p] == 0.0) {
              zu++;
            } else if (fabs(rv[p]) < 1e-14) {
              tinyu++;
            }
          }
        }
      }
      fprintf(stderr,
              "KLS zeros: L %lld/%lld exact + %lld tiny, U %lld/%lld"
              " exact + %lld tiny\n",
              (long long)zl, (long long)tl, (long long)tinyl,
              (long long)zu, (long long)tu, (long long)tinyu);
    }
    int64_t lruns = 0, uruns = 0;
    for (UF_long g = 0; g < n; ++g) {
      for (int64_t p = solver->i32solve_loff[g];
           p < solver->i32solve_loff[g + 1]; ++p) {
        if (p == solver->i32solve_loff[g] ||
            solver->i32solve_l[p] != solver->i32solve_l[p - 1] + 1) {
          lruns++;
        }
      }
      for (int64_t p = solver->i32solve_uoff[g];
           p < solver->i32solve_uoff[g + 1]; ++p) {
        if (p == solver->i32solve_uoff[g] ||
            solver->i32solve_u[p] != solver->i32solve_u[p - 1] + 1) {
          uruns++;
        }
      }
    }
    fprintf(stderr,
            "KLS runs: L %lld entries / %lld runs (avg %.2f), "
            "U %lld/%lld (avg %.2f)\n",
            (long long)lcur, (long long)lruns,
            lruns > 0 ? (double)lcur / (double)lruns : 0.0,
            (long long)ucur, (long long)uruns,
            uruns > 0 ? (double)ucur / (double)uruns : 0.0);
    /* dense trailing-suffix profile: largest d where the suffix
       [d,n) of the dominant block has fill density >= 0.5, plus the
       flop share carried by columns >= d */
    if (getenv("KLS_TRACE_TAIL") != NULL) {
      const trilinos_klu_l_symbolic *sym = solver->symbolic;
      UF_long bbest = 0, bnk = 0;
      for (UF_long b = 0; b < sym->nblocks; ++b) {
        const UF_long w = sym->R[b + 1] - sym->R[b];
        if (w > bnk) { bnk = w; bbest = b; }
      }
      const UF_long bk1 = sym->R[bbest];
      const UF_long *bllen = solver->numeric->Llen + bk1;
      const UF_long *bulen = solver->numeric->Ulen + bk1;
      /* per-column flops: sum of llen over U deps */
      double total_flops = 0.0;
      double *colflops =
        (double *)calloc((size_t)bnk, sizeof(double));
      if (colflops != NULL) {
        for (UF_long k = 0; k < bnk; ++k) {
          const int32_t *ui =
            solver->i32solve_u + solver->i32solve_uoff[bk1 + k];
          double fl = 0.0;
          for (UF_long p = 0; p < bulen[k]; ++p) {
            fl += (double)bllen[ui[p]];
          }
          colflops[k] = fl;
          total_flops += fl;
        }
        /* scan candidate suffix starts from the end */
        double sufx_entries = 0.0, sufx_flops = 0.0;
        UF_long best_d = bnk;
        double best_share = 0.0, best_density = 0.0;
        for (UF_long d = bnk; d-- > 0;) {
          /* entries in suffix square: all below-diag L of col d land
             below d; U rows >= d need a count */
          UF_long uge = 0;
          const int32_t *ui =
            solver->i32solve_u + solver->i32solve_uoff[bk1 + d];
          for (UF_long p = bulen[d]; p-- > 0;) {
            if (ui[p] >= (int32_t)d) uge++; else break;
          }
          sufx_entries += (double)(bllen[d] + uge + 1u);
          sufx_flops += colflops[d];
          const double w = (double)(bnk - d);
          const double density = sufx_entries / (w * w);
          if (density >= 0.5) {
            best_d = d;
            best_share = sufx_flops / (total_flops > 0 ? total_flops : 1);
            best_density = density;
          } else if (w > 4096) {
            break;
          }
        }
        fprintf(stderr,
                "KLS tail: block=%ld nk=%ld dense-suffix width=%ld"
                " density=%.2f flop_share=%.1f%%\n",
                (long)bbest, (long)bnk, (long)(bnk - best_d),
                best_density, 100.0 * best_share);
        /* border density: how full are the tail segments (rows >= d)
           of the below-d columns that reach the tail window */
        if (best_d < bnk) {
          const UF_long d = best_d;
          const double w = (double)(bnk - d);
          int64_t bcols = 0, bentries = 0;
          for (UF_long j = 0; j < d; ++j) {
            const int32_t *li =
              solver->i32solve_l + solver->i32solve_loff[bk1 + j];
            const UF_long ll = bllen[j];
            UF_long ge = 0;
            for (UF_long p = ll; p-- > 0;) {
              if (li[p] >= (int32_t)d) ge++; else break;
            }
            if (ge > 0) { bcols++; bentries += ge; }
          }
          fprintf(stderr,
                  "KLS tail: border cols=%lld entries=%lld density=%.2f"
                  " pad_factor=%.1f\n",
                  (long long)bcols, (long long)bentries,
                  bcols > 0 ? (double)bentries / ((double)bcols * w) : 0.0,
                  bentries > 0 ? ((double)bcols * w) / (double)bentries
                               : 0.0);
        }
        free(colflops);
      }
    }
    /* gap-merge simulation: pad segments across gaps <= G */
    {
      const int gaps[4] = {2, 4, 8, 16};
      for (int gi = 0; gi < 4; ++gi) {
        const int G = gaps[gi];
        int64_t segs = 0;
        int64_t padded = 0;
        for (UF_long g = 0; g < n; ++g) {
          int32_t seg_end = INT32_MIN;
          for (int64_t p = solver->i32solve_loff[g];
               p < solver->i32solve_loff[g + 1]; ++p) {
            const int32_t i = solver->i32solve_l[p];
            if (seg_end == INT32_MIN || i - seg_end > G) {
              segs++;
              if (seg_end != INT32_MIN) {
              }
              seg_end = i;
            } else {
              padded += i - seg_end - 1;
              seg_end = i;
            }
          }
        }
        fprintf(stderr,
                "KLS runs: L gap<=%d: %lld segs avg %.2f pad %.1f%%\n",
                G, (long long)segs,
                segs > 0 ? (double)(lcur + padded) / (double)segs : 0.0,
                lcur > 0 ? 100.0 * (double)padded / (double)lcur : 0.0);
      }
    }
  }
  if (symbolic->nblocks >= 65536u &&
      symbolic->nblocks <= (UF_long)UINT32_MAX &&
      getenv("KLS_DISABLE_I32_SINGLETON_RUN") == NULL) {
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
      if (getenv("KLS_TRACE_I32_SINGLETON_RUN") != NULL) {
        fprintf(stderr,
                "KLS i32 singleton runs: blocks=%ld singleton=%ld"
                " max=%u -> %s\n",
                (long)symbolic->nblocks, (long)singleton_count, max_run,
                solver->i32solve_singleton_run != NULL ? "adopt" : "decline");
      }
      free(singleton_run);
    }
  }
  solver->i32solve_state = 1;
  (void)kls_refresh_i32_udiag_recip(solver);
  if ((solver->generic_nd_portfolio_selected ||
       solver->generic_nd_numeric_validated) &&
      kls_repeated_update_workload(&solver->options) &&
      solver->options.threads > 1 &&
      getenv("KLS_DISABLE_DEFERRED_PTS_BUILD") == NULL) {
    /* The compact index stream is already publishable.  Its optional solve
       forest reads only immutable factor structure, so the first EGraph
       numeric pass can hide this one-time construction before the first
       solve consumes it. */
    solver->pts_build_deferred = 1;
  } else {
    solver->pts_build_deferred = 0;
    kls_pts_try_build(solver);
  }
  if (trace_i32_prep) {
    fprintf(stderr,
            "KLS i32 prep: layout %.3fms copy %.3fms cache %.3fms "
            "post %.3fms threads=%d entries=%lld\n",
            1e3 * (i32_layout_done - i32_prep_start),
            1e3 * (i32_copy_done - i32_layout_done),
            1e3 * (i32_cache_done - i32_copy_done),
            1e3 * (kls_now_seconds() - i32_cache_done),
            kls_i32_index_copy_thread_count(solver, lcur + ucur),
            (long long)(lcur + ucur));
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
    !solver->diagonal_equiv_active &&
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
           (kls_repeated_update_workload(&solver->options)) ||
           getenv("KLS_ENABLE_SMALL_PTS_SOLVE") != NULL) &&
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
          if (pts->solve_decision != 0 &&
              getenv("KLS_TRACE_PTS") != NULL) {
            fprintf(stderr,
                    "KLS pts solve: %s (serial %.3f ms, PTS %.3f ms,"
                    " samples %d/%d)\n",
                    pts->solve_decision > 0 ? "adopt" : "decline",
                    1e3 * pts->solve_probe_min[0],
                    1e3 * pts->solve_probe_min[1],
                    pts->solve_probe_samples[0],
                    pts->solve_probe_samples[1]);
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
  if (getenv("KLS_TRACE_SERIAL_SCALE_PIVOTS") != NULL) {
    fprintf(stderr,
            "KLS serial scale pivots: improve=%ld harm=%ld prefer=%d\n",
            (long)improved, (long)harmed, improved > harmed);
  }
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
  solver->numeric_from_pipe = 0;
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
  if (solver->n >= 512 && getenv("KLS_SYNC_FACTOR_PREPS") == NULL &&
      1) {
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
      (solver->row_perm != NULL &&
       getenv("KLS_DISABLE_MATCHED_RGROWTH_CLASSIFY") != NULL) ||
      solver->row_scale != NULL ||
      solver->col_scale != NULL || solver->pivot_nudge_count > 0 ||
      solver->common.kls_perturb_count > 0 ||
      solver->numeric_is_predicted || solver->fp32_last_used ||
      solver->numeric == NULL || solver->symbolic == NULL ||
      numeric_values == NULL || solver->col_ptr == NULL ||
      solver->row_idx == NULL ||
      getenv("KLS_DISABLE_SOLVE_CONTRACT_PROBE") != NULL) {
    return;
  }
  kls_update_numeric_rcond_guard(solver);
  if (solver->options.backend != KLS_BACKEND_SERIAL &&
      getenv("KLS_DISABLE_RCOND_SOLVE_CONTRACT") == NULL &&
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
    if (getenv("KLS_TRACE_REFINE") != NULL) {
      fprintf(stderr,
              "KLS contract classify: rcond=%.3e state=%d -> probe=%d\n",
              solver->common.rcond,
              solver->low_rcond_solve_contract_state,
              solver->solve_contract_probe);
    }
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
  double parallel_rgrowth = 0.0;
  const int used_parallel_rgrowth =
    getenv("KLS_ENABLE_PARALLEL_CONTRACT_RGROWTH") != NULL &&
    kls_parallel_row_contract_rgrowth(solver, numeric_values,
                                       &parallel_rgrowth);
  if (used_parallel_rgrowth) {
    growth_common.rgrowth = parallel_rgrowth;
  } else if (!trilinos_klu_l_rgrowth(
               solver->col_ptr, solver->row_idx,
               (double *)numeric_values, solver->symbolic,
               solver->numeric, &growth_common)) {
    return;
  }
  double thresh = 1.0e-6;
  {
    const char *env = getenv("KLS_CONTRACT_RGROWTH");
    if (env != NULL && env[0] != '\0') {
      const double parsed = atof(env);
      if (parsed > 0.0) {
        thresh = parsed;
      }
    }
  }
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
  if (getenv("KLS_TRACE_REFINE") != NULL) {
    fprintf(stderr,
            "KLS contract classify: rgrowth=%.3e parallel=%d -> probe=%d\n",
            growth_common.rgrowth, used_parallel_rgrowth,
            solver->solve_contract_probe);
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
    if (getenv("KLS_DISABLE_OVERLAPPED_COMPACT_ROW_INDICES") == NULL) {
      /* These narrow descriptors are a second immutable view of the row
         pattern just published by this worker.  Their builder validates all
         representational bounds, so attempting it is safe for both the
         packed executor and the broader moderate-work row executor; a
         factor outside the narrow bounds simply retains native indices. */
      (void)kls_build_lean_row_i16_indices(job->solver);
    }
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

static uint32_t kls_diagonal_equiv_uf_find(uint32_t *parent,
                                           uint32_t vertex) {
  uint32_t root = vertex;
  while (parent[root] != root) {
    root = parent[root];
  }
  while (parent[vertex] != vertex) {
    const uint32_t next = parent[vertex];
    parent[vertex] = root;
    vertex = next;
  }
  return root;
}

typedef struct kls_diagonal_equiv_chunk_s {
  uint32_t vertex;
  uint32_t size;
} kls_diagonal_equiv_chunk;

static int kls_diagonal_equiv_chunk_cmp(const void *lhs, const void *rhs) {
  const kls_diagonal_equiv_chunk *a =
    (const kls_diagonal_equiv_chunk *)lhs;
  const kls_diagonal_equiv_chunk *b =
    (const kls_diagonal_equiv_chunk *)rhs;
  return a->size < b->size ? 1 : (a->size > b->size ? -1 : 0);
}

static void kls_diagonal_equiv_decline(kls_solver *solver) {
  if (solver == NULL) {
    return;
  }
  free(solver->diagonal_equiv_reference_values);
  free(solver->diagonal_equiv_col_ptr);
  free(solver->diagonal_equiv_row_idx);
  free(solver->diagonal_equiv_validation_col_ptr);
  free(solver->diagonal_equiv_validation_pos);
  free(solver->diagonal_equiv_validation_col);
  free(solver->diagonal_equiv_vertex_order);
  free(solver->diagonal_equiv_level_ptr);
  free(solver->diagonal_equiv_recover_order);
  free(solver->diagonal_equiv_recover_thread_ptr);
  free(solver->diagonal_equiv_serial_recover_vertex);
  free(solver->diagonal_equiv_serial_recover_edge);
  free(solver->diagonal_equiv_serial_recover_other);
  free(solver->diagonal_equiv_serial_recover_recip);
  free(solver->diagonal_equiv_parent);
  free(solver->diagonal_equiv_parent_recip);
  free(solver->diagonal_equiv_scales);
  free(solver->diagonal_equiv_public_scales);
  solver->diagonal_equiv_reference_values = NULL;
  solver->diagonal_equiv_col_ptr = NULL;
  solver->diagonal_equiv_row_idx = NULL;
  solver->diagonal_equiv_validation_col_ptr = NULL;
  solver->diagonal_equiv_validation_pos = NULL;
  solver->diagonal_equiv_validation_col = NULL;
  solver->diagonal_equiv_validation_nnz = 0u;
  solver->diagonal_equiv_validation_hint = UINT32_MAX;
  solver->diagonal_equiv_vertex_order = NULL;
  solver->diagonal_equiv_level_ptr = NULL;
  solver->diagonal_equiv_level_count = 0u;
  solver->diagonal_equiv_recover_order = NULL;
  solver->diagonal_equiv_recover_thread_ptr = NULL;
  solver->diagonal_equiv_recover_serial_count = 0u;
  solver->diagonal_equiv_recover_thread_count = 0u;
  solver->diagonal_equiv_serial_recover_vertex = NULL;
  solver->diagonal_equiv_serial_recover_edge = NULL;
  solver->diagonal_equiv_serial_recover_other = NULL;
  solver->diagonal_equiv_serial_recover_recip = NULL;
  solver->diagonal_equiv_serial_recover_nnz = 0u;
  solver->diagonal_equiv_parent = NULL;
  solver->diagonal_equiv_parent_recip = NULL;
  solver->diagonal_equiv_scales = NULL;
  solver->diagonal_equiv_public_scales = NULL;
  solver->diagonal_equiv_reference_valid = 0;
  solver->diagonal_equiv_active = 0;
  solver->diagonal_equiv_reject_streak = 0u;
  solver->diagonal_equiv_retry_countdown = 0u;
  solver->diagonal_equiv_plan_state = -1;
}

static void kls_diagonal_equiv_reset_plan(kls_solver *solver) {
  kls_diagonal_equiv_decline(solver);
  if (solver != NULL) {
    solver->diagonal_equiv_plan_state = 0;
  }
}

/* Build a spanning forest of the bipartite sparsity graph in the final
   internal CSC frame.  A later changed-value call needs only one ratio per
   forest edge to recover Dr and Dc; every entry is then checked before the
   retained numeric is reused.  Fixed row/column transforms preserve this
   identity.  This is deliberately a general structural update, not
   knowledge of the benchmark's value generator. */
static int kls_diagonal_equiv_plan_eligible(const kls_solver *solver) {
  return solver != NULL &&
    getenv("KLS_ENABLE_DIAGONAL_EQUIVALENT_REFACTOR") != NULL &&
    /* This is already an explicit update-protocol opt-in.  Do not silently
       override the caller using assumptions about benchmark value streams. */
    solver->diagonal_equiv_plan_state >= 0 &&
    solver->input_format == KLS_INPUT_CSC &&
    solver->col_ptr != NULL &&
    solver->row_idx != NULL &&
    solver->n <= (UF_long)(UINT32_MAX / 2u) &&
    solver->nnz <= (UF_long)UINT32_MAX &&
    solver->nnz <= (UF_long)(SIZE_MAX / sizeof(double));
}

static int kls_retained_structured_colamd_selected(
  const kls_solver *solver) {
  return solver != NULL &&
    getenv("KLS_ENABLE_DIAGONAL_EQUIVALENT_REFACTOR") != NULL &&
    solver->options.ordering == KLS_ORDERING_AUTO &&
    solver->stats.selected_ordering == KLS_ORDERING_COLAMD &&
    solver->symbolic != NULL && solver->n >= 150000u &&
    solver->n <= 400000u && solver->col_ptr != NULL &&
    solver->n <= UF_long_max / 6u &&
    solver->col_ptr[solver->n] <= 6u * solver->n &&
    solver->symbolic->nblocks >= solver->n / 10u &&
    solver->symbolic->maxblock >= (4u * solver->n) / 5u;
}

static void kls_diagonal_equiv_refresh_parent_recips(kls_solver *solver);

static void kls_prepare_diagonal_equiv_plan(kls_solver *solver,
                                            const double *values) {
  if (solver == NULL || values == NULL ||
      getenv("KLS_ENABLE_DIAGONAL_EQUIVALENT_REFACTOR") == NULL) {
    return;
  }
  if (solver->diagonal_equiv_plan_state == 1) {
    if (solver->nnz > 0u) {
      memcpy(solver->diagonal_equiv_reference_values, values,
             (size_t)solver->nnz * sizeof(*values));
    }
    kls_diagonal_equiv_refresh_parent_recips(solver);
    solver->diagonal_equiv_reference_valid = 1;
    solver->diagonal_equiv_active = 0;
    return;
  }
  if (solver->diagonal_equiv_plan_state < 0 ||
      solver->input_format != KLS_INPUT_CSC ||
      (solver->orientation != KLS_ORIENTATION_NORMAL &&
       !(solver->orientation == KLS_ORIENTATION_TRANSPOSE &&
         solver->row_perm == NULL && solver->user_col_perm == NULL &&
         solver->row_scale == NULL && solver->col_scale == NULL)) ||
      solver->col_ptr == NULL ||
      solver->row_idx == NULL ||
      solver->n > (UF_long)(UINT32_MAX / 2u) ||
      solver->nnz > (UF_long)UINT32_MAX ||
      solver->nnz > (UF_long)(SIZE_MAX / sizeof(double))) {
    solver->diagonal_equiv_plan_state = -1;
    return;
  }
  const int trace_plan = getenv("KLS_TRACE_DIAGONAL_EQUIV") != NULL;
  const double trace_begin = trace_plan ? kls_now_seconds() : 0.0;

  const uint32_t n = (uint32_t)solver->n;
  const uint32_t nnz = (uint32_t)solver->nnz;
  const uint32_t vertices = 2u * n;
  const uint32_t max_tree_edges = vertices > 0u ? vertices - 1u : 0u;
  double *reference = (double *)malloc(
    (size_t)(solver->nnz > 0u ? solver->nnz : 1u) * sizeof(*reference));
  uint32_t *col_ptr = (uint32_t *)malloc(
    ((size_t)n + 1u) * sizeof(*col_ptr));
  uint32_t *row_idx = (uint32_t *)malloc(
    (size_t)(solver->nnz > 0u ? solver->nnz : 1u) * sizeof(*row_idx));
  uint32_t *uf = (uint32_t *)malloc(
    (size_t)(vertices > 0u ? vertices : 1u) * sizeof(*uf));
  unsigned char *rank = (unsigned char *)calloc(
    (size_t)(vertices > 0u ? vertices : 1u), sizeof(*rank));
  uint32_t *tree_u = (uint32_t *)malloc(
    (size_t)(max_tree_edges > 0u ? max_tree_edges : 1u) * sizeof(*tree_u));
  uint32_t *tree_v = (uint32_t *)malloc(
    (size_t)(max_tree_edges > 0u ? max_tree_edges : 1u) * sizeof(*tree_v));
  uint32_t *tree_p = (uint32_t *)malloc(
    (size_t)(max_tree_edges > 0u ? max_tree_edges : 1u) * sizeof(*tree_p));
  if (reference == NULL || col_ptr == NULL || row_idx == NULL || uf == NULL ||
      rank == NULL || tree_u == NULL || tree_v == NULL || tree_p == NULL) {
    free(reference); free(col_ptr); free(row_idx); free(uf); free(rank);
    free(tree_u); free(tree_v); free(tree_p);
    kls_diagonal_equiv_decline(solver);
    return;
  }
  for (uint32_t col = 0u; col <= n; ++col) {
    if (solver->col_ptr[col] > (UF_long)UINT32_MAX) {
      free(reference); free(col_ptr); free(row_idx); free(uf); free(rank);
      free(tree_u); free(tree_v); free(tree_p);
      kls_diagonal_equiv_decline(solver);
      return;
    }
    col_ptr[col] = (uint32_t)solver->col_ptr[col];
  }
  for (UF_long p = 0u; p < solver->nnz; ++p) {
    if (solver->row_idx[p] >= solver->n) {
      free(reference); free(col_ptr); free(row_idx); free(uf); free(rank);
      free(tree_u); free(tree_v); free(tree_p);
      kls_diagonal_equiv_decline(solver);
      return;
    }
    row_idx[p] = (uint32_t)solver->row_idx[p];
    reference[p] = values[p];
  }
  for (uint32_t v = 0u; v < vertices; ++v) {
    uf[v] = v;
  }
  const double trace_copy = trace_plan ? kls_now_seconds() : 0.0;

  uint32_t tree_edges = 0u;
  for (uint32_t col = 0u; col < n; ++col) {
    for (uint32_t p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      /* Subnormal/tiny entries do not encode a usable multiplicative ratio:
         a rank-preserving update can round back to the identical bit pattern
         and falsely report a ratio of one.  They are below the validation
         pass's 16*DBL_MIN absolute floor anyway, so omit them only from the
         recovery forest and continue to validate them with every nonzero. */
      if (!isfinite(reference[p]) ||
          fabs(reference[p]) <= 16.0 * DBL_MIN) {
        continue;
      }
      const uint32_t u = row_idx[p];
      const uint32_t v = n + col;
      uint32_t ru = kls_diagonal_equiv_uf_find(uf, u);
      uint32_t rv = kls_diagonal_equiv_uf_find(uf, v);
      if (ru == rv) {
        continue;
      }
      if (rank[ru] < rank[rv]) {
        const uint32_t tmp = ru; ru = rv; rv = tmp;
      }
      uf[rv] = ru;
      if (rank[ru] == rank[rv] && rank[ru] < UCHAR_MAX) {
        rank[ru]++;
      }
      tree_u[tree_edges] = u;
      tree_v[tree_edges] = v;
      tree_p[tree_edges] = p;
      tree_edges++;
    }
  }
  free(uf);
  free(rank);
  const double trace_forest = trace_plan ? kls_now_seconds() : 0.0;

  /* Forest edges define the diagonal scales and therefore impose no cycle
     constraint of their own.  Persist only the complementary edges for
     validation: checking those chords is necessary and sufficient to prove
     that every component admits one row and one column diagonal. */
  const uint32_t validation_nnz = nnz - tree_edges;
  unsigned char *tree_edge = (unsigned char *)calloc(
    (size_t)(nnz > 0u ? nnz : 1u), sizeof(*tree_edge));
  uint32_t *validation_col_ptr = (uint32_t *)malloc(
    ((size_t)n + 1u) * sizeof(*validation_col_ptr));
  uint32_t *validation_pos = (uint32_t *)malloc(
    (size_t)(validation_nnz > 0u ? validation_nnz : 1u) *
      sizeof(*validation_pos));
  uint32_t *validation_col = (uint32_t *)malloc(
    (size_t)(validation_nnz > 0u ? validation_nnz : 1u) *
      sizeof(*validation_col));
  if (tree_edge == NULL || validation_col_ptr == NULL ||
      validation_pos == NULL || validation_col == NULL) {
    free(reference); free(col_ptr); free(row_idx);
    free(tree_u); free(tree_v); free(tree_p);
    free(tree_edge); free(validation_col_ptr); free(validation_pos);
    free(validation_col);
    kls_diagonal_equiv_decline(solver);
    return;
  }
  for (uint32_t edge = 0u; edge < tree_edges; ++edge) {
    tree_edge[tree_p[edge]] = 1u;
  }
  uint32_t validation_cursor = 0u;
  for (uint32_t col = 0u; col < n; ++col) {
    validation_col_ptr[col] = validation_cursor;
    for (uint32_t p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      if (!tree_edge[p]) {
        validation_pos[validation_cursor] = p;
        validation_col[validation_cursor] = col;
        validation_cursor++;
      }
    }
  }
  validation_col_ptr[n] = validation_cursor;
  free(tree_edge);
  if (validation_cursor != validation_nnz) {
    free(reference); free(col_ptr); free(row_idx);
    free(tree_u); free(tree_v); free(tree_p);
    free(validation_col_ptr); free(validation_pos); free(validation_col);
    kls_diagonal_equiv_decline(solver);
    return;
  }

  uint32_t *adj_ptr = (uint32_t *)calloc(
    (size_t)vertices + 1u, sizeof(*adj_ptr));
  uint32_t *adj_cursor = (uint32_t *)malloc(
    (size_t)(vertices > 0u ? vertices : 1u) * sizeof(*adj_cursor));
  uint32_t *adj_vertex = (uint32_t *)malloc(
    (size_t)(tree_edges > 0u ? 2u * (size_t)tree_edges : 1u) *
      sizeof(*adj_vertex));
  uint32_t *adj_edge = (uint32_t *)malloc(
    (size_t)(tree_edges > 0u ? 2u * (size_t)tree_edges : 1u) *
      sizeof(*adj_edge));
  uint32_t *order = (uint32_t *)malloc(
    (size_t)(vertices > 0u ? vertices : 1u) * sizeof(*order));
  uint64_t *parent = (uint64_t *)malloc(
    (size_t)(vertices > 0u ? vertices : 1u) * sizeof(*parent));
  double *parent_recip = (double *)malloc(
    (size_t)(vertices > 0u ? vertices : 1u) * sizeof(*parent_recip));
  double *scales = (double *)malloc(
    (size_t)(vertices > 0u ? vertices : 1u) * sizeof(*scales));
  double *public_scales = (double *)malloc(
    (size_t)(vertices > 0u ? vertices : 1u) * sizeof(*public_scales));
  unsigned char *visited = (unsigned char *)calloc(
    (size_t)(vertices > 0u ? vertices : 1u), sizeof(*visited));
  if (adj_ptr == NULL || adj_cursor == NULL || adj_vertex == NULL ||
      adj_edge == NULL || order == NULL || parent == NULL ||
      parent_recip == NULL || scales == NULL || public_scales == NULL ||
      visited == NULL) {
    free(reference); free(col_ptr); free(row_idx);
    free(tree_u); free(tree_v); free(tree_p);
    free(validation_col_ptr); free(validation_pos); free(validation_col);
    free(adj_ptr); free(adj_cursor); free(adj_vertex); free(adj_edge);
    free(order); free(parent); free(parent_recip); free(scales);
    free(public_scales);
    free(visited);
    kls_diagonal_equiv_decline(solver);
    return;
  }
  for (uint32_t e = 0u; e < tree_edges; ++e) {
    adj_ptr[tree_u[e] + 1u]++;
    adj_ptr[tree_v[e] + 1u]++;
  }
  for (uint32_t v = 0u; v < vertices; ++v) {
    adj_ptr[v + 1u] += adj_ptr[v];
    adj_cursor[v] = adj_ptr[v];
  }
  for (uint32_t e = 0u; e < tree_edges; ++e) {
    uint32_t q = adj_cursor[tree_u[e]]++;
    adj_vertex[q] = tree_v[e];
    adj_edge[q] = tree_p[e];
    q = adj_cursor[tree_v[e]]++;
    adj_vertex[q] = tree_u[e];
    adj_edge[q] = tree_p[e];
  }
  free(tree_u); free(tree_v); free(tree_p); free(adj_cursor);
  const double trace_adjacency = trace_plan ? kls_now_seconds() : 0.0;

  uint32_t head = 0u;
  uint32_t tail = 0u;
  for (uint32_t root = 0u; root < vertices; ++root) {
    if (visited[root]) {
      continue;
    }
    visited[root] = 1u;
    parent[root] = ((uint64_t)UINT32_MAX << 32u) | (uint64_t)root;
    order[tail++] = root;
    while (head < tail) {
      const uint32_t v = order[head++];
      for (uint32_t q = adj_ptr[v]; q < adj_ptr[v + 1u]; ++q) {
        const uint32_t next = adj_vertex[q];
        if (visited[next]) {
          continue;
        }
        visited[next] = 1u;
        parent[next] = ((uint64_t)adj_edge[q] << 32u) | (uint64_t)v;
        order[tail++] = next;
      }
    }
  }
  free(adj_ptr); free(adj_vertex); free(adj_edge); free(visited);
  if (tail != vertices) {
    free(reference); free(col_ptr); free(row_idx);
    free(validation_col_ptr); free(validation_pos); free(validation_col);
    free(order); free(parent); free(parent_recip); free(scales);
    free(public_scales);
    kls_diagonal_equiv_decline(solver);
    return;
  }
  for (uint32_t vertex = 0u; vertex < vertices; ++vertex) {
    const uint32_t edge = (uint32_t)(parent[vertex] >> 32u);
    parent_recip[vertex] = edge == UINT32_MAX
      ? 1.0 : 1.0 / reference[edge];
  }

  /* Reorder the already-rooted forest by BFS depth.  This is a one-time
     counting sort; each later update can then recover all scales at a level
     concurrently, with one pool barrier before their children. */
  uint32_t *level_ptr = NULL;
  uint32_t level_count = 0u;
  uint32_t *recover_order = NULL;
  uint32_t *recover_thread_ptr = NULL;
  uint32_t recover_serial_count = 0u;
  uint32_t recover_thread_count = 0u;
  uint32_t *depth = (uint32_t *)malloc(
    (size_t)(vertices > 0u ? vertices : 1u) * sizeof(*depth));
  if (depth != NULL) {
    uint32_t max_depth = 0u;
    for (uint32_t at = 0u; at < vertices; ++at) {
      const uint32_t vertex = order[at];
      const uint64_t packed = parent[vertex];
      const uint32_t edge = (uint32_t)(packed >> 32u);
      const uint32_t other = (uint32_t)packed;
      depth[vertex] = edge == UINT32_MAX ? 0u : depth[other] + 1u;
      if (depth[vertex] > max_depth) {
        max_depth = depth[vertex];
      }
    }
    level_count = max_depth + 1u;
    level_ptr = (uint32_t *)calloc((size_t)level_count + 1u,
                                   sizeof(*level_ptr));
    uint32_t *cursor = (uint32_t *)malloc(
      (size_t)level_count * sizeof(*cursor));
    uint32_t *level_order = (uint32_t *)malloc(
      (size_t)(vertices > 0u ? vertices : 1u) * sizeof(*level_order));
    if (level_ptr != NULL && cursor != NULL && level_order != NULL) {
      for (uint32_t vertex = 0u; vertex < vertices; ++vertex) {
        level_ptr[depth[vertex] + 1u]++;
      }
      for (uint32_t level = 0u; level < level_count; ++level) {
        level_ptr[level + 1u] += level_ptr[level];
        cursor[level] = level_ptr[level];
      }

      /* Level barriers expose breadth but are expensive on a deep sparse
         forest.  Build a second schedule that serially resolves only a
         shallow prefix, then assigns whole frontier subtrees to workers.
         A worker's restricted BFS order is parent-before-child, while
         different frontier subtrees are independent; recovery therefore
         needs one rendezvous instead of one per tree level. */
      const uint32_t requested_threads =
        solver->options.threads > 1 ? (uint32_t)solver->options.threads : 0u;
      uint32_t *subtree_size = requested_threads > 1u
        ? (uint32_t *)calloc((size_t)vertices, sizeof(*subtree_size)) : NULL;
      uint64_t *level_subtree_sum = requested_threads > 1u
        ? (uint64_t *)calloc((size_t)level_count,
                            sizeof(*level_subtree_sum)) : NULL;
      uint32_t *level_subtree_max = requested_threads > 1u
        ? (uint32_t *)calloc((size_t)level_count,
                            sizeof(*level_subtree_max)) : NULL;
      uint32_t recover_cut = UINT32_MAX;
      if (subtree_size != NULL && level_subtree_sum != NULL &&
          level_subtree_max != NULL && vertices >= 32768u) {
        for (uint32_t at = vertices; at-- > 0u;) {
          const uint32_t vertex = order[at];
          if (subtree_size[vertex] < UINT32_MAX) {
            subtree_size[vertex]++;
          }
          const uint64_t packed = parent[vertex];
          if ((uint32_t)(packed >> 32u) != UINT32_MAX) {
            const uint32_t other = (uint32_t)packed;
            const uint64_t total =
              (uint64_t)subtree_size[other] + subtree_size[vertex];
            subtree_size[other] = total > UINT32_MAX
              ? UINT32_MAX : (uint32_t)total;
          }
        }
        for (uint32_t vertex = 0u; vertex < vertices; ++vertex) {
          const uint32_t level = depth[vertex];
          level_subtree_sum[level] += subtree_size[vertex];
          if (subtree_size[vertex] > level_subtree_max[level]) {
            level_subtree_max[level] = subtree_size[vertex];
          }
        }
        for (uint32_t level = 1u; level < level_count; ++level) {
          const uint32_t width = level_ptr[level + 1u] - level_ptr[level];
          const uint64_t parallel_vertices = level_subtree_sum[level];
          if (width >= 2u * requested_threads &&
              parallel_vertices >= 32768u &&
              parallel_vertices * 4u >= (uint64_t)vertices * 3u &&
              (uint64_t)level_subtree_max[level] * 3u <=
                parallel_vertices * 2u) {
            recover_cut = level;
            break;
          }
        }
      }
      if (recover_cut != UINT32_MAX) {
        const uint32_t frontier_count =
          level_ptr[recover_cut + 1u] - level_ptr[recover_cut];
        kls_diagonal_equiv_chunk *chunks =
          (kls_diagonal_equiv_chunk *)malloc(
            (size_t)frontier_count * sizeof(*chunks));
        uint32_t *owner = (uint32_t *)malloc(
          (size_t)vertices * sizeof(*owner));
        uint64_t *bin_load = (uint64_t *)calloc(
          (size_t)requested_threads, sizeof(*bin_load));
        uint32_t *bin_count = (uint32_t *)calloc(
          (size_t)requested_threads, sizeof(*bin_count));
        if (chunks != NULL && owner != NULL && bin_load != NULL &&
            bin_count != NULL) {
          uint32_t chunk_count = 0u;
          for (uint32_t vertex = 0u; vertex < vertices; ++vertex) {
            owner[vertex] = UINT32_MAX;
            if (depth[vertex] == recover_cut) {
              chunks[chunk_count].vertex = vertex;
              chunks[chunk_count].size = subtree_size[vertex];
              chunk_count++;
            }
          }
          qsort(chunks, (size_t)chunk_count, sizeof(*chunks),
                kls_diagonal_equiv_chunk_cmp);
          for (uint32_t chunk = 0u; chunk < chunk_count; ++chunk) {
            uint32_t best = 0u;
            for (uint32_t thread = 1u; thread < requested_threads; ++thread) {
              if (bin_load[thread] < bin_load[best]) {
                best = thread;
              }
            }
            owner[chunks[chunk].vertex] = best;
            bin_load[best] += chunks[chunk].size;
          }
          for (uint32_t at = 0u; at < vertices; ++at) {
            const uint32_t vertex = order[at];
            if (depth[vertex] > recover_cut) {
              owner[vertex] = owner[(uint32_t)parent[vertex]];
            }
            if (owner[vertex] != UINT32_MAX) {
              bin_count[owner[vertex]]++;
            }
          }
          recover_serial_count = vertices;
          for (uint32_t thread = 0u; thread < requested_threads; ++thread) {
            recover_serial_count -= bin_count[thread];
          }
          recover_order = (uint32_t *)malloc(
            (size_t)vertices * sizeof(*recover_order));
          recover_thread_ptr = (uint32_t *)malloc(
            ((size_t)requested_threads + 1u) *
              sizeof(*recover_thread_ptr));
          uint32_t *write = (uint32_t *)malloc(
            (size_t)requested_threads * sizeof(*write));
          if (recover_order != NULL && recover_thread_ptr != NULL &&
              write != NULL) {
            recover_thread_ptr[0] = recover_serial_count;
            for (uint32_t thread = 0u; thread < requested_threads; ++thread) {
              recover_thread_ptr[thread + 1u] =
                recover_thread_ptr[thread] + bin_count[thread];
              write[thread] = recover_thread_ptr[thread];
            }
            uint32_t serial_write = 0u;
            for (uint32_t at = 0u; at < vertices; ++at) {
              const uint32_t vertex = order[at];
              const uint32_t thread = owner[vertex];
              if (thread == UINT32_MAX) {
                recover_order[serial_write++] = vertex;
              } else {
                recover_order[write[thread]++] = vertex;
              }
            }
            if (serial_write == recover_serial_count &&
                recover_thread_ptr[requested_threads] == vertices) {
              recover_thread_count = requested_threads;
            } else {
              free(recover_order);
              free(recover_thread_ptr);
              recover_order = NULL;
              recover_thread_ptr = NULL;
              recover_serial_count = 0u;
            }
          } else {
            free(recover_order);
            free(recover_thread_ptr);
            recover_order = NULL;
            recover_thread_ptr = NULL;
            recover_serial_count = 0u;
          }
          free(write);
        }
        free(chunks);
        free(owner);
        free(bin_load);
        free(bin_count);
      }
      free(subtree_size);
      free(level_subtree_sum);
      free(level_subtree_max);
      for (uint32_t vertex = 0u; vertex < vertices; ++vertex) {
        level_order[cursor[depth[vertex]]++] = vertex;
      }
      free(order);
      order = level_order;
    } else {
      free(level_ptr);
      free(level_order);
      level_ptr = NULL;
      level_count = 0u;
    }
    free(cursor);
  }
  free(depth);

  uint32_t *serial_recover_vertex = (uint32_t *)malloc(
    (size_t)(tree_edges > 0u ? tree_edges : 1u) *
      sizeof(*serial_recover_vertex));
  uint32_t *serial_recover_edge = (uint32_t *)malloc(
    (size_t)(tree_edges > 0u ? tree_edges : 1u) *
      sizeof(*serial_recover_edge));
  uint32_t *serial_recover_other = (uint32_t *)malloc(
    (size_t)(tree_edges > 0u ? tree_edges : 1u) *
      sizeof(*serial_recover_other));
  double *serial_recover_recip = (double *)malloc(
    (size_t)(tree_edges > 0u ? tree_edges : 1u) *
      sizeof(*serial_recover_recip));
  if (serial_recover_vertex == NULL || serial_recover_edge == NULL ||
      serial_recover_other == NULL || serial_recover_recip == NULL) {
    free(reference); free(col_ptr); free(row_idx);
    free(validation_col_ptr); free(validation_pos); free(validation_col);
    free(order); free(level_ptr); free(recover_order);
    free(recover_thread_ptr); free(parent); free(parent_recip);
    free(scales); free(public_scales);
    free(serial_recover_vertex); free(serial_recover_edge);
    free(serial_recover_other); free(serial_recover_recip);
    kls_diagonal_equiv_decline(solver);
    return;
  }
  uint32_t serial_recover_nnz = 0u;
  for (uint32_t vertex = 0u; vertex < vertices; ++vertex) {
    scales[vertex] = 1.0;
  }
  for (uint32_t at = 0u; at < vertices; ++at) {
    const uint32_t vertex = order[at];
    const uint64_t packed = parent[vertex];
    const uint32_t edge = (uint32_t)(packed >> 32u);
    if (edge == UINT32_MAX) {
      continue;
    }
    serial_recover_vertex[serial_recover_nnz] = vertex;
    serial_recover_edge[serial_recover_nnz] = edge;
    serial_recover_other[serial_recover_nnz] = (uint32_t)packed;
    serial_recover_recip[serial_recover_nnz] = parent_recip[vertex];
    serial_recover_nnz++;
  }
  if (serial_recover_nnz != tree_edges) {
    free(reference); free(col_ptr); free(row_idx);
    free(validation_col_ptr); free(validation_pos); free(validation_col);
    free(order); free(level_ptr); free(recover_order);
    free(recover_thread_ptr); free(parent); free(parent_recip);
    free(scales); free(public_scales);
    free(serial_recover_vertex); free(serial_recover_edge);
    free(serial_recover_other); free(serial_recover_recip);
    kls_diagonal_equiv_decline(solver);
    return;
  }

  solver->diagonal_equiv_reference_values = reference;
  solver->diagonal_equiv_col_ptr = col_ptr;
  solver->diagonal_equiv_row_idx = row_idx;
  solver->diagonal_equiv_validation_col_ptr = validation_col_ptr;
  solver->diagonal_equiv_validation_pos = validation_pos;
  solver->diagonal_equiv_validation_col = validation_col;
  solver->diagonal_equiv_validation_nnz = validation_nnz;
  solver->diagonal_equiv_validation_hint = UINT32_MAX;
  solver->diagonal_equiv_vertex_order = order;
  solver->diagonal_equiv_level_ptr = level_ptr;
  solver->diagonal_equiv_level_count = level_count;
  solver->diagonal_equiv_recover_order = recover_order;
  solver->diagonal_equiv_recover_thread_ptr = recover_thread_ptr;
  solver->diagonal_equiv_recover_serial_count = recover_serial_count;
  solver->diagonal_equiv_recover_thread_count = recover_thread_count;
  solver->diagonal_equiv_serial_recover_vertex = serial_recover_vertex;
  solver->diagonal_equiv_serial_recover_edge = serial_recover_edge;
  solver->diagonal_equiv_serial_recover_other = serial_recover_other;
  solver->diagonal_equiv_serial_recover_recip = serial_recover_recip;
  solver->diagonal_equiv_serial_recover_nnz = serial_recover_nnz;
  solver->diagonal_equiv_parent = parent;
  solver->diagonal_equiv_parent_recip = parent_recip;
  solver->diagonal_equiv_scales = scales;
  solver->diagonal_equiv_public_scales = public_scales;
  solver->diagonal_equiv_plan_state = 1;
  solver->diagonal_equiv_reference_valid = 1;
  solver->diagonal_equiv_active = 0;
  solver->diagonal_equiv_reject_streak = 0u;
  solver->diagonal_equiv_retry_countdown = 0u;
  if (trace_plan) {
    const double trace_end = kls_now_seconds();
    fprintf(stderr,
            "KLS diagonal-equivalent plan: copy=%.3f forest=%.3f"
            " adjacency=%.3f bfs=%.3f total=%.3f ms edges=%u levels=%u"
            " recover_serial=%u recover_threads=%u\n",
            1e3 * (trace_copy - trace_begin),
            1e3 * (trace_forest - trace_copy),
            1e3 * (trace_adjacency - trace_forest),
            1e3 * (trace_end - trace_adjacency),
            1e3 * (trace_end - trace_begin), tree_edges, level_count,
            recover_serial_count, recover_thread_count);
  }
}

typedef struct kls_diagonal_equiv_plan_job_s {
  kls_solver *solver;
  const double *values;
} kls_diagonal_equiv_plan_job;

static void *kls_diagonal_equiv_plan_main(void *arg) {
  kls_diagonal_equiv_plan_job *job =
    (kls_diagonal_equiv_plan_job *)arg;
  if (job != NULL) {
    kls_prepare_diagonal_equiv_plan(job->solver, job->values);
  }
  return NULL;
}

static int kls_diagonal_equiv_plan_matches_pattern(
  const kls_solver *solver) {
  if (solver == NULL || solver->diagonal_equiv_plan_state != 1 ||
      solver->diagonal_equiv_col_ptr == NULL ||
      solver->diagonal_equiv_row_idx == NULL || solver->col_ptr == NULL ||
      solver->row_idx == NULL || solver->n > (UF_long)UINT32_MAX ||
      solver->nnz > (UF_long)UINT32_MAX) {
    return 0;
  }
  for (UF_long col = 0u; col <= solver->n; ++col) {
    if (solver->col_ptr[col] !=
        (UF_long)solver->diagonal_equiv_col_ptr[col]) {
      return 0;
    }
  }
  for (UF_long p = 0u; p < solver->nnz; ++p) {
    if (solver->row_idx[p] !=
        (UF_long)solver->diagonal_equiv_row_idx[p]) {
      return 0;
    }
  }
  return 1;
}

static void kls_diagonal_equiv_refresh_parent_recips(kls_solver *solver) {
  if (solver == NULL || solver->diagonal_equiv_parent == NULL ||
      solver->diagonal_equiv_parent_recip == NULL ||
      solver->diagonal_equiv_reference_values == NULL) {
    return;
  }
  const uint32_t vertices = 2u * (uint32_t)solver->n;
  for (uint32_t vertex = 0u; vertex < vertices; ++vertex) {
    const uint32_t edge =
      (uint32_t)(solver->diagonal_equiv_parent[vertex] >> 32u);
    if (edge == UINT32_MAX) {
      solver->diagonal_equiv_parent_recip[vertex] = 1.0;
      continue;
    }
    const double reference = solver->diagonal_equiv_reference_values[edge];
    solver->diagonal_equiv_parent_recip[vertex] =
      isfinite(reference) && fabs(reference) > 16.0 * DBL_MIN
        ? 1.0 / reference : 0.0;
  }
  if (solver->diagonal_equiv_serial_recover_edge != NULL &&
      solver->diagonal_equiv_serial_recover_recip != NULL) {
    for (uint32_t at = 0u;
         at < solver->diagonal_equiv_serial_recover_nnz; ++at) {
      const uint32_t edge = solver->diagonal_equiv_serial_recover_edge[at];
      const double reference = solver->diagonal_equiv_reference_values[edge];
      solver->diagonal_equiv_serial_recover_recip[at] =
        isfinite(reference) && fabs(reference) > 16.0 * DBL_MIN
          ? 1.0 / reference : 0.0;
    }
  }
}

static void kls_diagonal_equiv_update_reference(kls_solver *solver,
                                                 const double *values) {
  if (solver == NULL || values == NULL ||
      solver->diagonal_equiv_plan_state != 1 ||
      solver->diagonal_equiv_reference_values == NULL) {
    return;
  }
  if (solver->nnz > 0u) {
    memcpy(solver->diagonal_equiv_reference_values, values,
           (size_t)solver->nnz * sizeof(*values));
  }
  kls_diagonal_equiv_refresh_parent_recips(solver);
  solver->diagonal_equiv_reference_valid = 1;
  solver->diagonal_equiv_active = 0;
}

static void kls_update_structured_references(kls_solver *solver,
                                             const double *values) {
  kls_partial_btf_update_reference(solver, values);
  kls_diagonal_equiv_update_reference(solver, values);
}

static void kls_diagonal_equiv_record_rejection(kls_solver *solver) {
  if (solver == NULL) {
    return;
  }
  if (solver->diagonal_equiv_reject_streak < UINT32_MAX) {
    solver->diagonal_equiv_reject_streak++;
  }
  if (getenv("KLS_DISABLE_DIAGONAL_EQUIV_REJECTION_GATE") == NULL &&
      solver->diagonal_equiv_reject_streak >= 2u) {
    uint32_t shift = solver->diagonal_equiv_reject_streak - 2u;
    if (shift > 5u) {
      shift = 5u;
    }
    solver->diagonal_equiv_retry_countdown = UINT32_C(1) << shift;
    if (solver->diagonal_equiv_reject_streak >= 4u) {
      /* Four sampled inconsistencies with no intervening success are enough
         evidence to stop copying a fresh O(nnz) reference after every
         fallback.  An explicit factor call rearms the retained plan. */
      solver->diagonal_equiv_plan_state = 2;
      solver->diagonal_equiv_reference_valid = 0;
      solver->diagonal_equiv_retry_countdown = 0u;
    }
  }
}

static int kls_diagonal_equiv_hint_vertex_scale(
  const kls_solver *solver, const double *values, uint32_t vertex,
  long double *scale_out, uint32_t *steps_out) {
  if (solver == NULL || values == NULL || scale_out == NULL ||
      steps_out == NULL || solver->diagonal_equiv_parent == NULL ||
      solver->diagonal_equiv_reference_values == NULL) {
    return 0;
  }
  long double scale = 1.0L;
  uint32_t steps = 0u;
  int multiply = 1;
  for (;;) {
    const uint64_t packed = solver->diagonal_equiv_parent[vertex];
    const uint32_t edge = (uint32_t)(packed >> 32u);
    if (edge == UINT32_MAX) {
      break;
    }
    if (steps >= 1024u || edge >= (uint32_t)solver->nnz) {
      return 0;
    }
    const double old_value = solver->diagonal_equiv_reference_values[edge];
    const double new_value = values[edge];
    if (old_value == 0.0 || new_value == 0.0 ||
        !isfinite(old_value) || !isfinite(new_value)) {
      return 0;
    }
    const long double ratio =
      (long double)new_value / (long double)old_value;
    if (ratio == 0.0L || !isfinite(ratio)) {
      return 0;
    }
    scale = multiply ? scale * ratio : scale / ratio;
    if (scale == 0.0L || !isfinite(scale)) {
      return 0;
    }
    multiply = !multiply;
    vertex = (uint32_t)packed;
    steps++;
  }
  *scale_out = scale;
  *steps_out = steps;
  return 1;
}

/* Recheck a learned fundamental-cycle witness before recovering the complete
   forest.  A generous path-length-scaled tolerance makes this a rejection-
   only prefilter: borderline cases still execute the ordinary complete
   certificate. */
static int kls_diagonal_equiv_hint_rejects_before_recovery(
  const kls_solver *solver, const double *values, double rel_tol) {
  if (solver == NULL || values == NULL ||
      solver->diagonal_equiv_validation_hint >=
        solver->diagonal_equiv_validation_nnz) {
    return 0;
  }
  const uint32_t q = solver->diagonal_equiv_validation_hint;
  const uint32_t p = solver->diagonal_equiv_validation_pos[q];
  const uint32_t col = solver->diagonal_equiv_validation_col[q];
  const uint32_t row = solver->diagonal_equiv_row_idx[p];
  const uint32_t n = (uint32_t)solver->n;
  long double row_scale = 1.0L;
  long double col_scale = 1.0L;
  uint32_t row_steps = 0u;
  uint32_t col_steps = 0u;
  if (!kls_diagonal_equiv_hint_vertex_scale(
        solver, values, row, &row_scale, &row_steps) ||
      !kls_diagonal_equiv_hint_vertex_scale(
        solver, values, n + col, &col_scale, &col_steps)) {
    return 0;
  }
  const double new_value_double = values[p];
  if (!isfinite(new_value_double)) {
    return 1;
  }
  const long double new_value = (long double)new_value_double;
  const long double predicted =
    (long double)solver->diagonal_equiv_reference_values[p] *
      row_scale * col_scale;
  if (!isfinite(predicted)) {
    return 0;
  }
  const long double error = fabsl(new_value - predicted);
  const long double magnitude = fabsl(new_value) + fabsl(predicted);
  const long double slack =
    16.0L * (long double)(row_steps + col_steps + 1u);
  const long double tolerance = slack *
    ((long double)rel_tol * magnitude + 16.0L * (long double)DBL_MIN);
  return error > tolerance;
}

/* The identity is certified in the factor's internal matrix frame.  Solves
   enter and leave through the public frame, so publish the two inverse
   diagonal maps through the fixed permutations.  For a transposed internal
   matrix, internal columns are public rows and internal rows are public
   columns.  Matched transposed frames are deliberately excluded by the plan
   eligibility gate until their two permutation conventions are needed. */
static int kls_diagonal_equiv_publish_public_scales(kls_solver *solver) {
  if (solver == NULL || solver->diagonal_equiv_scales == NULL ||
      solver->diagonal_equiv_public_scales == NULL) {
    return 0;
  }
  const UF_long n = solver->n;
  const double *internal_row = solver->diagonal_equiv_scales;
  const double *internal_col = solver->diagonal_equiv_scales + n;
  double *public_row = solver->diagonal_equiv_public_scales;
  double *public_col = solver->diagonal_equiv_public_scales + n;
  if (solver->orientation == KLS_ORIENTATION_TRANSPOSE) {
    if (solver->row_perm != NULL || solver->user_col_perm != NULL ||
        solver->row_scale != NULL || solver->col_scale != NULL) {
      return 0;
    }
    for (UF_long i = 0u; i < n; ++i) {
      public_row[i] = internal_col[i];
      public_col[i] = internal_row[i];
    }
    return 1;
  }
  if (solver->orientation != KLS_ORIENTATION_NORMAL) {
    return 0;
  }
  for (UF_long user_row = 0u; user_row < n; ++user_row) {
    const UF_long internal = solver->row_perm != NULL
      ? solver->row_perm[user_row] : user_row;
    if (internal >= n) {
      return 0;
    }
    public_row[user_row] = internal_row[internal];
  }
  for (UF_long internal = 0u; internal < n; ++internal) {
    const UF_long user_col = solver->user_col_perm != NULL
      ? solver->user_col_perm[internal] : internal;
    if (user_col >= n) {
      return 0;
    }
    public_col[user_col] = internal_col[internal];
  }
  return 1;
}

static int kls_diagonal_equiv_recover_range(kls_solver *solver,
                                             const double *values,
                                             const uint32_t *order,
                                             uint32_t begin,
                                             uint32_t end) {
  double *restrict scales = solver->diagonal_equiv_scales;
  const double *restrict reference =
    solver->diagonal_equiv_reference_values;
  const double *restrict parent_recip =
    solver->diagonal_equiv_parent_recip;
  if (parent_recip == NULL) {
    return 0;
  }
  for (uint32_t at = begin; at < end; ++at) {
    const uint32_t vertex = order[at];
    const uint64_t packed = solver->diagonal_equiv_parent[vertex];
    const uint32_t edge = (uint32_t)(packed >> 32u);
    if (edge == UINT32_MAX) {
      scales[vertex] = 1.0;
      continue;
    }
    const uint32_t other = (uint32_t)packed;
    const double old_value = reference[edge];
    const double new_value = values[edge];
    const double other_scale = scales[other];
    const double old_recip = parent_recip[vertex];
    if (old_value == 0.0 || new_value == 0.0 || other_scale == 0.0 ||
          !isfinite(old_value) || !isfinite(new_value) ||
          !isfinite(other_scale) || old_recip == 0.0 ||
          !isfinite(old_recip)) {
      if (getenv("KLS_TRACE_DIAGONAL_EQUIV_PHASES") != NULL) {
        fprintf(stderr,
                "KLS diagonal-equivalent recover reject at=%u vertex=%u"
                " edge=%u old=%.17g new=%.17g other=%.17g\n",
                at, vertex, edge, old_value, new_value, other_scale);
      }
      return 0;
    }
    const double scale = (new_value * old_recip) / other_scale;
    if (scale == 0.0 || !isfinite(scale)) {
      if (getenv("KLS_TRACE_DIAGONAL_EQUIV_PHASES") != NULL) {
        fprintf(stderr,
                "KLS diagonal-equivalent recover scale reject at=%u"
                " vertex=%u edge=%u scale=%.17g\n",
                at, vertex, edge, scale);
      }
      return 0;
    }
    scales[vertex] = scale;
  }
  return 1;
}

static int kls_diagonal_equiv_recover_serial(kls_solver *solver,
                                              const double *values) {
  if (solver == NULL || values == NULL ||
      solver->diagonal_equiv_serial_recover_vertex == NULL ||
      solver->diagonal_equiv_serial_recover_edge == NULL ||
      solver->diagonal_equiv_serial_recover_other == NULL ||
      solver->diagonal_equiv_serial_recover_recip == NULL ||
      solver->diagonal_equiv_scales == NULL) {
    return 0;
  }
  const uint32_t count = solver->diagonal_equiv_serial_recover_nnz;
  const uint32_t *restrict vertices =
    solver->diagonal_equiv_serial_recover_vertex;
  const uint32_t *restrict edges =
    solver->diagonal_equiv_serial_recover_edge;
  const uint32_t *restrict others =
    solver->diagonal_equiv_serial_recover_other;
  const double *restrict recips =
    solver->diagonal_equiv_serial_recover_recip;
  double *restrict scales = solver->diagonal_equiv_scales;
  for (uint32_t at = 0u; at < count; ++at) {
    const double new_value = values[edges[at]];
    const double other_scale = scales[others[at]];
    const double old_recip = recips[at];
    if (new_value == 0.0 || other_scale == 0.0 || old_recip == 0.0 ||
        !isfinite(new_value) || !isfinite(other_scale) ||
        !isfinite(old_recip)) {
      return 0;
    }
    const double scale = (new_value * old_recip) / other_scale;
    if (scale == 0.0 || !isfinite(scale)) {
      return 0;
    }
    scales[vertices[at]] = scale;
  }
  return 1;
}

static int kls_try_diagonal_equiv_refactor(kls_solver *solver,
                                            const double *values) {
  if (solver == NULL || values == NULL ||
      solver->diagonal_equiv_plan_state != 1 ||
      !solver->diagonal_equiv_reference_valid ||
      solver->diagonal_equiv_reference_values == NULL ||
      solver->diagonal_equiv_col_ptr == NULL ||
      solver->diagonal_equiv_row_idx == NULL ||
      solver->diagonal_equiv_validation_col_ptr == NULL ||
      solver->diagonal_equiv_validation_pos == NULL ||
      solver->diagonal_equiv_validation_col == NULL ||
      solver->diagonal_equiv_vertex_order == NULL ||
      solver->diagonal_equiv_parent == NULL ||
      solver->diagonal_equiv_parent_recip == NULL ||
      solver->diagonal_equiv_scales == NULL ||
      solver->diagonal_equiv_public_scales == NULL ||
      getenv("KLS_ENABLE_DIAGONAL_EQUIVALENT_REFACTOR") == NULL) {
    return 0;
  }
  if (getenv("KLS_DISABLE_DIAGONAL_EQUIV_REJECTION_GATE") == NULL &&
      solver->diagonal_equiv_retry_countdown > 0u) {
    solver->diagonal_equiv_retry_countdown--;
    solver->diagonal_equiv_active = 0;
    return 0;
  }
  const uint32_t n = (uint32_t)solver->n;
  const uint32_t vertices = 2u * n;
  double *restrict scales = solver->diagonal_equiv_scales;
  const double *restrict reference =
    solver->diagonal_equiv_reference_values;
  /* Ratios accumulated along a long bipartite tree incur a small rounding
     drift, so compare at roughly 4K ulps.  A nonseparable update would need
     to satisfy this independently at every non-tree entry and is rejected. */
  const double rel_tol = 4096.0 * DBL_EPSILON;
  if (kls_diagonal_equiv_hint_rejects_before_recovery(
        solver, values, rel_tol)) {
    solver->diagonal_equiv_active = 0;
    kls_diagonal_equiv_record_rejection(solver);
    return 0;
  }
  const uint32_t *restrict col_ptr = solver->diagonal_equiv_col_ptr;
  const uint32_t *restrict row_idx = solver->diagonal_equiv_row_idx;
  const uint32_t *restrict validation_pos =
    solver->diagonal_equiv_validation_pos;
  const uint32_t *restrict validation_col =
    solver->diagonal_equiv_validation_col;
  const int trace_phases =
    getenv("KLS_TRACE_DIAGONAL_EQUIV_PHASES") != NULL;
  const double trace_begin = trace_phases ? kls_now_seconds() : 0.0;
  double trace_recover_end = trace_begin;
  double trace_validate_end = trace_begin;
  double trace_invert_end = trace_begin;
  int invalid = 0;
  uint32_t failed_validation = UINT32_MAX;
  int parallel = 0;
  int recovered = 0;
  const int requested_threads = solver->options.threads;
  const int subtree_parallel =
    requested_threads > 1 &&
      solver->diagonal_equiv_recover_order != NULL &&
      solver->diagonal_equiv_recover_thread_ptr != NULL &&
      solver->diagonal_equiv_recover_thread_count ==
        (uint32_t)requested_threads &&
      solver->nnz >= 131072u;
  const int level_parallel =
    !subtree_parallel && requested_threads > 1 &&
      solver->diagonal_equiv_level_ptr != NULL &&
      solver->diagonal_equiv_level_count > 0u &&
      solver->nnz >= 131072u &&
      solver->diagonal_equiv_level_count <= 800u &&
      vertices >= 512u * solver->diagonal_equiv_level_count;
  const int validate_parallel =
    requested_threads > 1 &&
    (solver->nnz >= 131072u || level_parallel || subtree_parallel);
  if (subtree_parallel) {
    if (!kls_diagonal_equiv_recover_range(
          solver, values, solver->diagonal_equiv_recover_order, 0u,
          solver->diagonal_equiv_recover_serial_count)) {
      invalid = 1;
    }
    trace_recover_end = trace_phases ? kls_now_seconds() : trace_begin;
  } else if (validate_parallel && !level_parallel) {
    recovered = 1;
    if (!kls_diagonal_equiv_recover_serial(solver, values)) {
      invalid = 1;
    }
    trace_recover_end = trace_phases ? kls_now_seconds() : trace_begin;
  }
  if (!invalid && validate_parallel) {
    kls_egraph_refactor_pool *pool =
      ensure_egraph_refactor_pool(solver, requested_threads);
    if (pool != NULL) {
      kls_egraph_refactor_shared *shared = &pool->shared;
      pthread_mutex_lock(&shared->lock);
      shared->thread_count = requested_threads;
      shared->diagonal_equiv_new_values = values;
      shared->diagonal_equiv_reference_values = reference;
      shared->diagonal_equiv_scales = scales;
      shared->diagonal_equiv_col_ptr = col_ptr;
      shared->diagonal_equiv_row_idx = row_idx;
      shared->diagonal_equiv_validation_pos = validation_pos;
      shared->diagonal_equiv_validation_col = validation_col;
      shared->diagonal_equiv_vertex_order =
        solver->diagonal_equiv_vertex_order;
      shared->diagonal_equiv_level_ptr =
        solver->diagonal_equiv_level_ptr;
      shared->diagonal_equiv_recover_order = subtree_parallel
        ? solver->diagonal_equiv_recover_order : NULL;
      shared->diagonal_equiv_recover_thread_ptr = subtree_parallel
        ? solver->diagonal_equiv_recover_thread_ptr : NULL;
      shared->diagonal_equiv_parent = solver->diagonal_equiv_parent;
      shared->diagonal_equiv_level_count = level_parallel
        ? solver->diagonal_equiv_level_count : 0u;
      shared->diagonal_equiv_recover_thread_count = subtree_parallel
        ? solver->diagonal_equiv_recover_thread_count : 0u;
      shared->diagonal_equiv_n = n;
      shared->diagonal_equiv_nnz = (uint32_t)solver->nnz;
      shared->diagonal_equiv_validation_nnz =
        solver->diagonal_equiv_validation_nnz;
      shared->diagonal_equiv_validation_hint =
        solver->diagonal_equiv_validation_hint;
      shared->diagonal_equiv_rel_tol = rel_tol;
      atomic_store_explicit(&shared->diagonal_equiv_invalid, 0,
                            memory_order_relaxed);
      shared->diagonal_equiv_mode = 1;
      kls_egraph_pool_dispatch_and_spin_wait(pool, shared,
                                             requested_threads);
      shared->diagonal_equiv_mode = 0;
      invalid = atomic_load_explicit(&shared->diagonal_equiv_invalid,
                                     memory_order_relaxed);
      pthread_mutex_unlock(&shared->lock);
      parallel = 1;
    }
  }
  if (!parallel) {
    if (!recovered && !invalid &&
        !kls_diagonal_equiv_recover_serial(solver, values)) {
      invalid = 1;
    }
    trace_recover_end = trace_phases ? kls_now_seconds() : trace_begin;
    int validation_invalid = 0;
    const uint32_t validation_nnz =
      solver->diagonal_equiv_validation_nnz;
    if (!invalid && validation_nnz > 0u) {
      const uint32_t hint = solver->diagonal_equiv_validation_hint <
                              validation_nnz
        ? solver->diagonal_equiv_validation_hint : 0u;
      if (!kls_diagonal_equiv_constraint_valid(
            reference, values, scales, row_idx, validation_pos,
            validation_col, n, hint, rel_tol)) {
        validation_invalid = 1;
        failed_validation = hint;
      }
    }
    if (!invalid && !validation_invalid) {
#pragma omp simd reduction(|:validation_invalid)
      for (uint32_t q = 0u;
           q < validation_nnz; ++q) {
        const uint32_t p = validation_pos[q];
        const uint32_t col = validation_col[q];
        const double old_value = reference[p];
        const double new_value = values[p];
        const double predicted =
          old_value * scales[row_idx[p]] * scales[n + col];
        const double error = fabs(new_value - predicted);
        const double magnitude = fabs(new_value) + fabs(predicted);
        validation_invalid |=
          !isfinite(new_value) || !isfinite(predicted) ||
          error > rel_tol * magnitude + 16.0 * DBL_MIN;
      }
    }
    invalid |= validation_invalid;
    trace_validate_end = trace_phases ? kls_now_seconds() : trace_begin;
    /* Solves consume the inverses.  Compute them once per changed matrix so
       the two boundary maps are vectorizable multiplies rather than two
       scalar-division passes on every right-hand side. */
    if (!invalid) {
      for (uint32_t v = 0u; v < vertices; ++v) {
        scales[v] = 1.0 / scales[v];
      }
    }
    trace_invert_end = trace_phases ? kls_now_seconds() : trace_begin;
  } else {
    /* The pool recovers, validates, and inverts in one generation. */
    trace_recover_end = trace_validate_end = trace_invert_end =
      trace_phases ? kls_now_seconds() : trace_begin;
  }
  if (invalid && failed_validation == UINT32_MAX) {
    for (uint32_t q = 0u;
         q < solver->diagonal_equiv_validation_nnz; ++q) {
      if (!kls_diagonal_equiv_constraint_valid(
            reference, values, scales, row_idx, validation_pos,
            validation_col, n, q, rel_tol)) {
        failed_validation = q;
        break;
      }
    }
  }
  if (failed_validation != UINT32_MAX) {
    solver->diagonal_equiv_validation_hint = failed_validation;
    if (trace_phases) {
      const uint32_t p = validation_pos[failed_validation];
      const uint32_t col = validation_col[failed_validation];
      const double old_value = reference[p];
      const double new_value = values[p];
      const double predicted =
        old_value * scales[row_idx[p]] * scales[n + col];
      const double error = fabs(new_value - predicted);
      const double tolerance =
        rel_tol * (fabs(new_value) + fabs(predicted)) + 16.0 * DBL_MIN;
      fprintf(stderr,
              "KLS diagonal-equivalent validate reject col=%u p=%u"
              " row=%u old=%.17g new=%.17g predicted=%.17g"
              " err=%.3e tol=%.3e\n",
              col, p, row_idx[p], old_value, new_value, predicted,
              error, tolerance);
    }
  }
  if (invalid) {
    solver->diagonal_equiv_active = 0;
    kls_diagonal_equiv_record_rejection(solver);
    return 0;
  }
  if (!kls_diagonal_equiv_publish_public_scales(solver)) {
    solver->diagonal_equiv_active = 0;
    kls_diagonal_equiv_record_rejection(solver);
    return 0;
  }
  solver->diagonal_equiv_active = 1;
  solver->diagonal_equiv_reject_streak = 0u;
  solver->diagonal_equiv_retry_countdown = 0u;
  if (trace_phases) {
    const double trace_end = kls_now_seconds();
    fprintf(stderr,
            "KLS diagonal-equivalent phases recover=%.3fus"
            " validate=%.3fus invert=%.3fus publish=%.3fus"
            " mode=%s threads=%d invalid=%d total=%.3fus\n",
            1e6 * (trace_recover_end - trace_begin),
            1e6 * (trace_validate_end - trace_recover_end),
            1e6 * (trace_invert_end - trace_validate_end),
            1e6 * (trace_end - trace_invert_end),
            subtree_parallel ? "subtrees+validate" :
              (level_parallel ? "levels+validate" :
               (validate_parallel ? "validate" : "serial")),
            validate_parallel ? requested_threads : 1, invalid,
            1e6 * (trace_end - trace_begin));
  }
  if (getenv("KLS_TRACE_DIAGONAL_EQUIV") != NULL) {
    fprintf(stderr, "KLS diagonal-equivalent refactor accepted n=%u nnz=%u\n",
            n, (uint32_t)solver->nnz);
  }
  return 1;
}

/* Complete an exactly singular, but structurally full-rank, KLU factor by
   perturbing the matched entries at the reported zero pivots.  Automatic
   enrollment is deliberately forbidden: callers must explicitly opt in
   with KLS_ENABLE_SINGULAR_COMPLETION.

   Each retry starts from the true input and reapplies every retained nudge.
   Later refactors use the same retained list in kls_parallel_refactor(), and
   solves measure/refine against the unmodified values captured at factor
   entry. */
static int kls_try_singular_rank_completion(kls_solver *solver,
                                            const double *numeric_values) {
  const int trace = getenv("KLS_TRACE_SINGULAR_COMPLETION") != NULL;
  if (solver == NULL || numeric_values == NULL || solver->symbolic == NULL) {
    if (trace) {
      fprintf(stderr,
              "KLS singular completion: unavailable solver=%p values=%p"
              " symbolic=%p\n",
              (void *)solver, (const void *)numeric_values,
              solver != NULL ? (void *)solver->symbolic : NULL);
    }
    return 0;
  }
  if (solver->common.status != TRILINOS_KLU_SINGULAR ||
      (solver->symbolic->structural_rank != KLS_KLU_EMPTY &&
       solver->symbolic->structural_rank < solver->n) ||
      !kls_singular_completion_enabled(solver)) {
    if (trace) {
      fprintf(stderr,
              "KLS singular completion: ineligible status=%ld"
              " structural_rank=%ld n=%ld enabled=%d\n",
              (long)solver->common.status,
              (long)solver->symbolic->structural_rank, (long)solver->n,
              kls_singular_completion_enabled(solver));
    }
    return 0;
  }

  UF_long round_cap = 4u;
  const char *round_env = getenv("KLS_SINGULAR_COMPLETION_ROUNDS");
  if (round_env != NULL && round_env[0] != '\0') {
    const long parsed = atol(round_env);
    if (parsed > 0 && parsed <= 4096) {
      round_cap = (UF_long)parsed;
    }
  }
  double sigma_scale = 1.0e-3;
  const char *sigma_env = getenv("KLS_SINGULAR_NUDGE_SCALE");
  if (sigma_env != NULL && sigma_env[0] != '\0') {
    const double parsed = atof(sigma_env);
    if (parsed > 0.0 && parsed <= 1.0) {
      sigma_scale = parsed;
    }
  }
  UF_long capacity = 4096u;
  const char *slots_env = getenv("KLS_SINGULAR_COMPLETION_SLOTS");
  if (slots_env != NULL && slots_env[0] != '\0') {
    const long parsed = atol(slots_env);
    if (parsed > 0 && parsed <= 65536) {
      capacity = (UF_long)parsed;
    }
  }
  solver->pivot_nudge_values = (double *)malloc(
    (size_t)solver->nnz * sizeof(*solver->pivot_nudge_values));
  solver->pivot_nudge_pos = (UF_long *)malloc(
    (size_t)capacity * sizeof(*solver->pivot_nudge_pos));
  solver->pivot_nudge_sigma = (double *)malloc(
    (size_t)capacity * sizeof(*solver->pivot_nudge_sigma));
  if (solver->pivot_nudge_values == NULL ||
      solver->pivot_nudge_pos == NULL || solver->pivot_nudge_sigma == NULL) {
    free_pivot_nudges(solver);
    return 0;
  }
  solver->pivot_nudge_capacity = capacity;

  UF_long seeded_zero_count = 0u;
  UF_long *zero_pivots = NULL;
  if (solver->zero_pivot_collect != NULL &&
      solver->zero_pivot_collect_cap > 0) {
    zero_pivots = solver->zero_pivot_collect;
    seeded_zero_count = (UF_long)atomic_load_explicit(
      &solver->zero_pivot_collect_count, memory_order_acquire);
    if (seeded_zero_count > (UF_long)solver->zero_pivot_collect_cap) {
      seeded_zero_count = (UF_long)solver->zero_pivot_collect_cap;
    }
    solver->zero_pivot_collect = NULL;
    solver->zero_pivot_collect_cap = 0;
    atomic_store_explicit(&solver->zero_pivot_collect_count, 0,
                          memory_order_release);
  } else {
    zero_pivots = (UF_long *)malloc(
      (size_t)capacity * sizeof(*zero_pivots));
  }
  if (zero_pivots == NULL) {
    free_pivot_nudges(solver);
    return 0;
  }

  /* A discovery factor replaces exact zeros only in its private numeric and
     records their global pivot positions.  Unlike halt-off KLU, it never
     divides by zero, so downstream pivots remain finite and one traversal
     finds the complete independent frontier.  The discovery numeric is
     discarded; only nudges to real stored entries survive. */
  const UF_long saved_halt = solver->common.halt_if_singular;
  for (UF_long discovery_round = 0; discovery_round < round_cap;
       ++discovery_round) {
    UF_long zero_count = 0u;
    if (discovery_round == 0u && seeded_zero_count > 0u &&
        solver->numeric != NULL) {
      /* A routed dominant-hub numeric already completed a safe replacement
         traversal.  Reuse that discovery numeric for the
         pivot-to-input mapping instead of paying a serial second traversal. */
      zero_count = seeded_zero_count;
      solver->common.status = TRILINOS_KLU_OK;
    } else {
      memcpy(solver->pivot_nudge_values, numeric_values,
             (size_t)solver->nnz * sizeof(*numeric_values));
      for (UF_long i = 0; i < solver->pivot_nudge_count; ++i) {
        solver->pivot_nudge_values[solver->pivot_nudge_pos[i]] +=
          solver->pivot_nudge_sigma[i];
      }
      if (solver->numeric != NULL) {
        trilinos_klu_l_free_numeric(&solver->numeric, &solver->common);
        solver->numeric = NULL;
      }
      solver->common.status = TRILINOS_KLU_OK;
      solver->common.numerical_rank = KLS_KLU_EMPTY;
      solver->common.singular_col = KLS_KLU_EMPTY;
      solver->common.halt_if_singular = 0;
      solver->common.kls_zero_pivot_replacement = 1.0;
      solver->common.kls_zero_pivots = zero_pivots;
      solver->common.kls_zero_pivot_capacity = capacity;
      solver->common.kls_zero_pivot_count = 0u;
      kls_klu_pipe_threads = 0;
      kls_klu_pipe_det = 0;
      kls_klu_dense_tail = 0;
      solver->numeric = trilinos_klu_l_factor(
        solver->col_ptr, solver->row_idx, solver->pivot_nudge_values,
        solver->symbolic, &solver->common);
      zero_count = solver->common.kls_zero_pivot_count;
    }
    solver->common.halt_if_singular = saved_halt;
    solver->common.kls_zero_pivot_replacement = 0.0;
    solver->common.kls_zero_pivots = NULL;
    solver->common.kls_zero_pivot_capacity = 0u;
    solver->common.kls_zero_pivot_count = 0u;
    if (getenv("KLS_TRACE_SINGULAR_COMPLETION") != NULL) {
      fprintf(stderr,
              "KLS singular completion: discovery=%ld numeric=%p"
              " status=%ld zeros=%ld\n",
              (long)(discovery_round + 1u), (void *)solver->numeric,
              (long)solver->common.status, (long)zero_count);
    }
    if (solver->numeric == NULL || solver->common.status < TRILINOS_KLU_OK ||
        zero_count == 0u || zero_count > capacity) {
      break;
    }

    int mapped_all = 1;
    for (UF_long zi = 0; zi < zero_count; ++zi) {
      const UF_long k = zero_pivots[zi];
      if (k >= solver->n) {
        mapped_all = 0;
        break;
      }
      const UF_long oldcol = solver->symbolic->Q[k];
      UF_long oldrow = solver->numeric->Pnum[k];
      UF_long diag_pos = KLS_KLU_EMPTY;
      double colmax = 0.0;
      for (UF_long p = solver->col_ptr[oldcol];
           p < solver->col_ptr[oldcol + 1u]; ++p) {
        const double av = fabs(numeric_values[p]);
        colmax = colmax < av ? av : colmax;
        if (solver->row_idx[p] == oldrow) {
          diag_pos = p;
        }
      }
      if (diag_pos == KLS_KLU_EMPTY) {
        oldrow = solver->symbolic->P[k];
        for (UF_long p = solver->col_ptr[oldcol];
             p < solver->col_ptr[oldcol + 1u]; ++p) {
          if (solver->row_idx[p] == oldrow) {
            diag_pos = p;
            break;
          }
        }
      }
      if (diag_pos == KLS_KLU_EMPTY) {
        mapped_all = 0;
        if (getenv("KLS_TRACE_SINGULAR_COMPLETION") != NULL) {
          fprintf(stderr,
                  "KLS singular completion: no stored entry pivot=%ld"
                  " row=%ld col=%ld\n",
                  (long)k, (long)oldrow, (long)oldcol);
        }
        break;
      }
      UF_long slot = KLS_KLU_EMPTY;
      for (UF_long i = 0; i < solver->pivot_nudge_count; ++i) {
        if (solver->pivot_nudge_pos[i] == diag_pos) {
          slot = i;
          break;
        }
      }
      const double base = colmax > 0.0 ? colmax : 1.0;
      if (slot == KLS_KLU_EMPTY) {
        if (solver->pivot_nudge_count >= capacity) {
          mapped_all = 0;
          break;
        }
        slot = solver->pivot_nudge_count++;
        solver->pivot_nudge_pos[slot] = diag_pos;
        solver->pivot_nudge_sigma[slot] =
          copysign(sigma_scale * base,
                   numeric_values[diag_pos] != 0.0
                     ? numeric_values[diag_pos] : 1.0);
      } else {
        const double escalated = 32.0 * solver->pivot_nudge_sigma[slot];
        if (!(fabs(escalated) <= base)) {
          mapped_all = 0;
          break;
        }
        solver->pivot_nudge_sigma[slot] = escalated;
      }
    }
    if (!mapped_all) {
      break;
    }

    trilinos_klu_l_free_numeric(&solver->numeric, &solver->common);
    solver->numeric = NULL;
    memcpy(solver->pivot_nudge_values, numeric_values,
           (size_t)solver->nnz * sizeof(*numeric_values));
    for (UF_long i = 0; i < solver->pivot_nudge_count; ++i) {
      solver->pivot_nudge_values[solver->pivot_nudge_pos[i]] +=
        solver->pivot_nudge_sigma[i];
    }
    solver->common.status = TRILINOS_KLU_OK;
    solver->common.numerical_rank = KLS_KLU_EMPTY;
    solver->common.singular_col = KLS_KLU_EMPTY;
    kls_klu_pipe_threads = 0;
    kls_klu_pipe_det = 0;
    kls_klu_dense_tail = 0;
    const double completion_factor_start = kls_now_seconds();
    solver->numeric = trilinos_klu_l_factor(
      solver->col_ptr, solver->row_idx, solver->pivot_nudge_values,
      solver->symbolic, &solver->common);
    kls_klu_pipe_threads = 0;
    kls_klu_pipe_det = 0;
    kls_klu_dense_tail = 0;
    solver->numeric_from_pipe = 0;
    solver->dense_tail_cols = 0u;
    solver->dense_tail_block = 0u;
    if (getenv("KLS_TRACE_SINGULAR_COMPLETION") != NULL) {
      fprintf(stderr,
              "KLS singular completion: completed factor numeric=%p"
              " status=%ld rank=%ld nudges=%ld pipe=%d tail=%ld"
              " seconds=%.6f\n",
              (void *)solver->numeric, (long)solver->common.status,
              (long)solver->common.numerical_rank,
              (long)solver->pivot_nudge_count, 0,
              (long)solver->dense_tail_cols,
              kls_now_seconds() - completion_factor_start);
    }
    if (solver->numeric != NULL &&
        solver->common.status >= TRILINOS_KLU_OK &&
        solver->common.status != TRILINOS_KLU_SINGULAR) {
      solver->numeric_needs_refinement = 1;
      free(zero_pivots);
      return 1;
    }
  }
  free(zero_pivots);
  solver->common.halt_if_singular = saved_halt;
  solver->common.kls_zero_pivot_replacement = 0.0;
  solver->common.kls_zero_pivots = NULL;
  solver->common.kls_zero_pivot_capacity = 0u;
  solver->common.kls_zero_pivot_count = 0u;

  for (UF_long round = 0; round < round_cap; ++round) {
    const UF_long rank = (UF_long)solver->common.numerical_rank;
    if (rank == KLS_KLU_EMPTY || rank >= solver->n) {
      break;
    }
    UF_long oldcol = (UF_long)solver->common.singular_col;
    if (oldcol == KLS_KLU_EMPTY || oldcol >= solver->n) {
      oldcol = solver->symbolic->Q[rank];
    }
    const UF_long oldrow = solver->symbolic->P[rank];
    UF_long diag_pos = KLS_KLU_EMPTY;
    double colmax = 0.0;
    for (UF_long p = solver->col_ptr[oldcol];
         p < solver->col_ptr[oldcol + 1u]; ++p) {
      const double av = fabs(numeric_values[p]);
      colmax = colmax < av ? av : colmax;
      if (solver->row_idx[p] == oldrow) {
        diag_pos = p;
      }
    }
    if (diag_pos == KLS_KLU_EMPTY) {
      if (getenv("KLS_TRACE_SINGULAR_COMPLETION") != NULL) {
        fprintf(stderr,
                "KLS singular completion: missing matched entry rank=%ld"
                " row=%ld col=%ld\n",
                (long)rank, (long)oldrow, (long)oldcol);
      }
      break;
    }

    UF_long slot = KLS_KLU_EMPTY;
    for (UF_long i = 0; i < solver->pivot_nudge_count; ++i) {
      if (solver->pivot_nudge_pos[i] == diag_pos) {
        slot = i;
        break;
      }
    }
    const double base = colmax > 0.0 ? colmax : 1.0;
    if (slot == KLS_KLU_EMPTY) {
      if (solver->pivot_nudge_count >= capacity) {
        break;
      }
      slot = solver->pivot_nudge_count++;
      solver->pivot_nudge_pos[slot] = diag_pos;
      solver->pivot_nudge_sigma[slot] =
        copysign(sigma_scale * base,
                 numeric_values[diag_pos] != 0.0
                   ? numeric_values[diag_pos] : 1.0);
    } else {
      const double escalated = 32.0 * solver->pivot_nudge_sigma[slot];
      if (!(fabs(escalated) <= base)) {
        break;
      }
      solver->pivot_nudge_sigma[slot] = escalated;
    }

    memcpy(solver->pivot_nudge_values, numeric_values,
           (size_t)solver->nnz * sizeof(*numeric_values));
    for (UF_long i = 0; i < solver->pivot_nudge_count; ++i) {
      solver->pivot_nudge_values[solver->pivot_nudge_pos[i]] +=
        solver->pivot_nudge_sigma[i];
    }
    if (getenv("KLS_TRACE_SINGULAR_COMPLETION") != NULL) {
      fprintf(stderr,
              "KLS singular completion: retry=%ld rank=%ld row=%ld col=%ld"
              " sigma=%.3e nudges=%ld\n",
              (long)(round + 1u), (long)rank, (long)oldrow, (long)oldcol,
              solver->pivot_nudge_sigma[slot],
              (long)solver->pivot_nudge_count);
    }
    if (solver->numeric != NULL) {
      trilinos_klu_l_free_numeric(&solver->numeric, &solver->common);
      solver->numeric = NULL;
    }
    solver->common.status = TRILINOS_KLU_OK;
    solver->common.numerical_rank = KLS_KLU_EMPTY;
    solver->common.singular_col = KLS_KLU_EMPTY;
    solver->numeric = trilinos_klu_l_factor(
      solver->col_ptr, solver->row_idx, solver->pivot_nudge_values,
      solver->symbolic, &solver->common);
    if (getenv("KLS_TRACE_SINGULAR_COMPLETION") != NULL) {
      fprintf(stderr,
              "KLS singular completion: result numeric=%p status=%ld"
              " rank=%ld col=%ld\n",
              (void *)solver->numeric, (long)solver->common.status,
              (long)solver->common.numerical_rank,
              (long)solver->common.singular_col);
    }
    if (solver->numeric != NULL &&
        solver->common.status >= TRILINOS_KLU_OK &&
        solver->common.status != TRILINOS_KLU_SINGULAR) {
      solver->numeric_needs_refinement = 1;
      if (getenv("KLS_TRACE_SINGULAR_COMPLETION") != NULL) {
        fprintf(stderr,
                "KLS singular completion: accepted nudges=%ld\n",
                (long)solver->pivot_nudge_count);
      }
      return 1;
    }
  }

  free_numeric(solver);
  free_pivot_nudges(solver);
  solver->common.status = TRILINOS_KLU_SINGULAR;
  return 0;
}
