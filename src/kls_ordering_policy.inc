/* Internal ordering and automatic-policy implementation.
 * Included by kls.c at the original definition point to preserve a single
 * translation unit and its optimization behavior.
 */

static int kls_low_work_btf_map32_capable(const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      solver->common.scale > 0 || solver->numeric->Rs != NULL ||
      solver->n == 0u || solver->n > (UF_long)INT32_MAX ||
      solver->nnz == 0u || solver->nnz > (UF_long)INT32_MAX ||
      solver->symbolic->nblocks <= 1u ||
      !(solver->common.flops > 0.0) ||
      !(solver->common.flops < 100000.0)) {
    return 0;
  }
  return solver->symbolic->structural_rank == KLS_KLU_EMPTY ||
         solver->symbolic->structural_rank == solver->n;
}

/* KLU's native fixed-pattern walk also has a useful low-work BTF envelope
   when KLU itself owns row scaling.  Unlike the map32 kernel, this path does
   not require an unscaled numeric; it only requires the unchanged public
   compressed frame certified separately at dispatch. */
static int kls_low_work_btf_public_klu_capable(const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      solver->n == 0u || solver->n > (UF_long)INT32_MAX ||
      solver->nnz == 0u || solver->nnz > (UF_long)INT32_MAX ||
      solver->symbolic->nblocks <= 1u ||
      !isfinite(solver->common.flops) || solver->common.flops < 0.0 ||
      solver->common.flops >= 100000.0) {
    return 0;
  }
  return solver->symbolic->structural_rank == KLS_KLU_EMPTY ||
         solver->symbolic->structural_rank == solver->n;
}

/* Admit moderate-work BTF factors only to a measured executor comparison.
   Unlike the settled low-work policy, this predicate never chooses a kernel:
   nearby factors with almost identical dimensions can have opposite packed-
   stream behavior, so the retained numeric must decide from actual timings. */
static int kls_moderate_btf_map32_trial_capable(const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      solver->common.scale > 0 || solver->numeric->Rs != NULL ||
      solver->n == 0u || solver->n > (UF_long)INT32_MAX ||
      solver->nnz == 0u || solver->nnz > (UF_long)INT32_MAX ||
      solver->symbolic->nblocks <= 1u ||
      solver->common.flops < 100000.0 ||
      solver->common.flops >= 500000.0 ||
      !kls_repeated_update_workload(&solver->options)) {
    return 0;
  }
  return solver->symbolic->structural_rank == KLS_KLU_EMPTY ||
         solver->symbolic->structural_rank == solver->n;
}

static int kls_low_work_btf_map32_policy_enabled(
  const kls_solver *solver) {
  return getenv("KLS_DISABLE_LOW_WORK_BTF_MAP32_POLICY") == NULL &&
         kls_low_work_btf_map32_capable(solver);
}

/* Prefer the native packed triangular walk when a complete low-work BTF
   factor is already small enough to stay cache-resident.  In this regime a
   duplicate i16/i32 index stream has too little solve work over which to
   amortize its construction and extra dispatch.  Larger systems retain the
   compact mirror: their longer vectors and index streams repay that setup.
   The limits are deliberately properties of the retained factor rather than
   dimensions or block counts of a named input matrix. */
static int kls_low_work_btf_prefers_native_solve(
  const kls_solver *solver) {
  if (getenv("KLS_DISABLE_LOW_WORK_BTF_NATIVE_SOLVE") != NULL ||
      !kls_low_work_btf_map32_capable(solver) || solver->n > 16384u) {
    return 0;
  }
  const UF_long factor_entry_cap = 12u * solver->n;
  if (solver->numeric->lnz > factor_entry_cap ||
      solver->numeric->unz > factor_entry_cap - solver->numeric->lnz) {
    return 0;
  }
  /* A short vector cannot repay a second index representation.  In the
     larger cache tier, retain native block-packed streams only when BTF
     components average at least 32 rows; highly fragmented factors instead
     benefit from the compact mirror's flatter, narrower traversal. */
  return solver->n <= 4096u ||
         solver->symbolic->nblocks <= solver->n / 32u;
}

/* Select the compact fixed-pivot machinery from properties of the retained
   factor, not from an input fingerprint.  The kernel is useful when
   bookkeeping dominates a genuinely low-work, cache-sized unscaled block.
   Above 4096 rows, the ordinary mapped walk's contiguous streams repay their
   setup and avoid the compact kernel's per-column dispatch; retained-factor
   measurements on both sides of that crossover keep the cutoff conservative.
   The compact map itself uses signed 32-bit positions. */
static int kls_low_work_single_block_policy_enabled(
  const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      getenv("KLS_DISABLE_LOW_WORK_SINGLE_BLOCK_POLICY") != NULL ||
      solver->orientation != KLS_ORIENTATION_NORMAL ||
      solver->common.scale > 0 || solver->numeric->Rs != NULL ||
      solver->n == 0u || solver->n > 4096u ||
      solver->nnz == 0u || solver->nnz > (UF_long)INT32_MAX ||
      solver->symbolic->nblocks != 1u ||
      solver->symbolic->maxblock != solver->n ||
      !(solver->common.flops > 0.0) ||
      !(solver->common.flops < 100000.0)) {
    return 0;
  }
  return solver->symbolic->structural_rank == KLS_KLU_EMPTY ||
         solver->symbolic->structural_rank == solver->n;
}

/* A compact one-block factor just above the low-work KLU crossover can still
   be dominated by mapped-column bookkeeping.  Select the retained scalar-row
   representation from factor work and fill rather than an input fingerprint. */
static int kls_moderate_work_single_block_lean_policy_enabled(
  const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      getenv("KLS_DISABLE_MODERATE_SINGLE_BLOCK_LEAN_POLICY") != NULL ||
      solver->orientation != KLS_ORIENTATION_NORMAL ||
      solver->options.threads < 2 || solver->common.scale != -1 ||
      !kls_repeated_update_workload(&solver->options) ||
      solver->numeric->Rs != NULL || solver->n == 0u ||
      solver->n > (UF_long)INT32_MAX ||
      solver->n > UF_long_max / 10u || solver->nnz == 0u ||
      solver->nnz > (UF_long)INT32_MAX ||
      solver->nnz > 10u * solver->n || solver->symbolic->do_btf ||
      solver->symbolic->nblocks != 1u ||
      solver->symbolic->maxblock != solver->n ||
      solver->common.flops < 100000.0 ||
      solver->common.flops > 1000000.0 ||
      solver->pivot_nudge_count != 0u ||
      solver->common.kls_perturb_count != 0u) {
    return 0;
  }
  const UF_long fill_cap = 10u * solver->n;
  if (solver->numeric->lnz > fill_cap ||
      solver->numeric->unz > fill_cap - solver->numeric->lnz) {
    return 0;
  }
  return solver->symbolic->structural_rank == KLS_KLU_EMPTY ||
         solver->symbolic->structural_rank == solver->n;
}

/* The same bookkeeping crossover exists for a compact fragmented factor.
   Select it from the retained BTF, numeric work, and fill rather than from
   the input pattern that happened to produce the components. */
static int kls_moderate_work_fragmented_btf_lean_policy_enabled(
  const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      getenv("KLS_DISABLE_MODERATE_FRAGMENTED_BTF_LEAN_POLICY") != NULL ||
      solver->options.threads < 2 || solver->common.scale != -1 ||
      !kls_repeated_update_workload(&solver->options) ||
      solver->numeric->Rs != NULL || solver->n == 0u ||
      solver->n > (UF_long)INT32_MAX ||
      solver->n > UF_long_max / 12u || solver->nnz == 0u ||
      solver->nnz > (UF_long)INT32_MAX ||
      !solver->symbolic->do_btf || solver->symbolic->nblocks <= 1u ||
      solver->common.flops < 100000.0 ||
      solver->common.flops > 1000000.0 ||
      solver->pivot_nudge_count != 0u ||
      solver->common.kls_perturb_count != 0u) {
    return 0;
  }
  const UF_long fill_cap = 12u * solver->n;
  if (solver->numeric->lnz > fill_cap ||
      solver->numeric->unz > fill_cap - solver->numeric->lnz) {
    return 0;
  }
  return solver->symbolic->structural_rank == KLS_KLU_EMPTY ||
         solver->symbolic->structural_rank == solver->n;
}

/* A dominant diagonal block makes compact column locality a useful proxy for
   the heavier row representation.  It is not such a proxy when the retained
   factor is split into thousands of small components: row execution can
   remove per-block dispatch even if the compact column challenger does not
   win its isolated walk.  Preserve a measured row trial for that structural
   regime instead of turning the compact preflight into a permanent veto. */
static int kls_highly_fragmented_btf_row_trial_capable(
  const kls_solver *solver) {
  return solver != NULL && solver->symbolic != NULL && solver->n > 0u &&
    solver->symbolic->do_btf &&
    solver->symbolic->nblocks >= solver->n / 4u &&
    solver->symbolic->maxblock <= solver->n - solver->n / 4u;
}

static int kls_row_refactor_env_enabled(void);

/* KLU's fixed-pattern refactor can consume the public value array without
   preparation whenever the public compressed matrix is the retained frame.
   That is normal orientation for CSC input and transpose orientation for CSR
   input: in both cases the caller arrays are the CSC frame factored by KLU. The
   symbolic P/Q and KLU's own row scaling remain internal to that call; KLS
   coordinate transforms do not.  Keep this as a representation capability
   so a measured engine tournament can reuse it at any work size. */
static int kls_direct_klu_public_frame_capable(const kls_solver *solver) {
  return solver != NULL && solver->symbolic != NULL &&
    solver->numeric != NULL &&
    !kls_row_refactor_env_enabled() &&
    solver->input_to_csc == NULL &&
    ((solver->input_format == KLS_INPUT_CSC &&
      solver->orientation == KLS_ORIENTATION_NORMAL) ||
     (solver->input_format == KLS_INPUT_CSR &&
      solver->orientation == KLS_ORIENTATION_TRANSPOSE)) &&
    solver->row_perm == NULL && solver->user_col_perm == NULL &&
    solver->row_scale == NULL && solver->col_scale == NULL &&
    ((solver->common.scale > 0 && solver->numeric->Rs != NULL) ||
     (solver->common.scale <= 0 && solver->numeric->Rs == NULL));
}

/* A mapped KLS refactor may retain a mathematically equivalent numeric while
   repacking storage that the vendor refactor does not understand.  Pointer
   and value-frame capability alone cannot prove that the vendor walk still
   owns that representation.  Certify a timed challenger with one ordinary
   solve and a backward residual before publishing its verdict.  The probe
   vector is deterministic but independent of matrix dimensions or identity;
   a failed challenger is followed by a complete incumbent refactor. */
static int kls_direct_klu_numeric_residual_probe(
  kls_solver *solver, const double *values) {
  if (solver == NULL || values == NULL || solver->symbolic == NULL ||
      solver->numeric == NULL || solver->n == 0u ||
      solver->n > (UF_long)(SIZE_MAX / (3u * sizeof(double)))) {
    return 0;
  }
  double *work = (double *)malloc(
    3u * (size_t)solver->n * sizeof(*work));
  if (work == NULL) {
    return 0;
  }
  double *rhs = work;
  double *x = rhs + solver->n;
  double *truth = x + solver->n;
  for (UF_long i = 0; i < solver->n; ++i) {
    rhs[i] = 0.0;
    truth[i] = 1.0 + (double)(i % 29u) / 29.0;
  }
  for (UF_long j = 0; j < solver->n; ++j) {
    const double xj = truth[j];
    for (UF_long p = solver->col_ptr[j]; p < solver->col_ptr[j + 1u]; ++p) {
      rhs[solver->row_idx[p]] += values[p] * xj;
    }
  }
  memcpy(x, rhs, (size_t)solver->n * sizeof(*x));
  solver->common.status = TRILINOS_KLU_OK;
  const int solved =
    trilinos_klu_l_solve(solver->symbolic, solver->numeric, solver->n,
                         1, x, &solver->common) &&
    solver->common.status >= TRILINOS_KLU_OK;
  double rhs_max = 0.0;
  double residual_max = INFINITY;
  if (solved) {
    residual_max = 0.0;
    for (UF_long i = 0; i < solver->n; ++i) {
      const double av = fabs(rhs[i]);
      rhs_max = rhs_max < av ? av : rhs_max;
    }
    memcpy(truth, rhs, (size_t)solver->n * sizeof(*truth));
    for (UF_long j = 0; j < solver->n; ++j) {
      const double xj = x[j];
      for (UF_long p = solver->col_ptr[j]; p < solver->col_ptr[j + 1u]; ++p) {
        truth[solver->row_idx[p]] -= values[p] * xj;
      }
    }
    for (UF_long i = 0; i < solver->n; ++i) {
      const double av = fabs(truth[i]);
      residual_max = residual_max < av ? av : residual_max;
    }
  }
  free(work);
  solver->common.status = TRILINOS_KLU_OK;
  return solved && isfinite(residual_max) &&
    residual_max <= 1.0e-9 * (rhs_max > 0.0 ? rhs_max : 1.0);
}

/* A direct KLU refactor may consume the caller's values only when the public
   matrix is already the retained, untransformed CSC frame.  Keep this
   lifecycle capability separate from the low-work kernel crossover so CSR,
   permutation and scaling users retain their ordinary
   preparation paths. */
static int kls_low_work_single_block_direct_public_capable(
  const kls_solver *solver) {
  return kls_low_work_single_block_policy_enabled(solver) &&
    kls_direct_klu_public_frame_capable(solver) &&
    solver->common.scale == -1;
}

enum {
  KLS_BALANCED_MODERATE_HUB_NONE = 0,
  KLS_BALANCED_MODERATE_HUB_DOMINANT_BTF = 1,
  KLS_BALANCED_MODERATE_HUB_COARSE_FRINGE_NO_BTF = 2
};

/* Degree-balanced sparse systems with a moderate two-way hub have two useful
   AMD regimes.  A nearly spanning SCC retains BTF, while a 75--90% core plus
   a bounded number of coarse fringe SCCs exposes substantially more repeated
   numeric parallelism as one AMD block.  Full structural diagonal permits an
   identity matching, so a direct strong-component pass classifies the two
   regimes without a maximum-transversal or an ordering tournament.  The
   degree and SCC predicates are invariant under transpose. */

/* Candidate analyzes can run on independent threads.  Keep their structural
   policy local to the analyze call so an AMF trial for another orientation
   cannot inherit this class verdict. */
/* Medium, almost fully diagonal circuits with one macroscopic spike are a
   distinct ordering regime: mean-local-fill AMMF removes substantially more
   update work than AMD without destroying the wide dependency front.  The
   bounds are resource tiers for the serial deficiency analysis; the spike
   proof, rather than a narrow order/nonzero window, separates ordinary
   bounded-degree circuits. */
static int is_medium_full_diagonal_spike_minfill_pattern(
  UF_long n,
  const UF_long *col_ptr,
  const UF_long *row_idx) {
  if (n < 65536u || n > 262144u || col_ptr == NULL || row_idx == NULL ||
      n > UF_long_max / 8u || col_ptr[n] < 4u * n || col_ptr[n] > 8u * n) {
    return 0;
  }

  UF_long diagonal_count = 0u;
  UF_long max_col_degree = 0u;
  for (UF_long col = 0u; col < n; ++col) {
    const UF_long degree = col_ptr[col + 1u] - col_ptr[col];
    if (degree > max_col_degree) {
      max_col_degree = degree;
    }
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      if (row_idx[p] >= n) {
        return 0;
      }
      diagonal_count += row_idx[p] == col;
    }
  }
  return 1000.0 * (double)diagonal_count >= 999.0 * (double)n &&
         4u * max_col_degree >= n && 2u * max_col_degree <= n;
}

/* The input spike is only a proposal.  Retain AMMF when its actual symbolic
   exposes one nearly spanning BTF core, bounded fragmentation, balanced fill,
   and enough work to amortize the dependency-aware recurring numeric. */
static int kls_medium_full_diagonal_spike_minfill_symbolic_contract(
  UF_long n,
  const trilinos_klu_l_symbolic *symbolic) {
  if (n < 65536u || n > 262144u || symbolic == NULL || !symbolic->do_btf ||
      symbolic->nblocks < 2u || symbolic->nblocks > n / 1024u ||
      symbolic->maxblock * 200u < 199u * n ||
      (symbolic->structural_rank != KLS_KLU_EMPTY &&
       symbolic->structural_rank != n) ||
      !(symbolic->lnz > 0.0) || !(symbolic->unz > 0.0) ||
      !(symbolic->est_flops > 0.0)) {
    return 0;
  }
  const double rows = (double)n;
  const double fill = symbolic->lnz + symbolic->unz;
  return fill >= 8.0 * rows && fill <= 32.0 * rows &&
    symbolic->est_flops >= 512.0 * rows &&
    symbolic->est_flops <= 2048.0 * rows &&
    symbolic->est_flops <= 5.0e8 &&
    symbolic->lnz <= 2.0 * symbolic->unz &&
    symbolic->unz <= 2.0 * symbolic->lnz;
}

static int kls_exact_structural_symmetry(
  UF_long n,
  const UF_long *col_ptr,
  const UF_long *row_idx,
  int columns_sorted,
  UF_long max_column_degree);

static int is_medium_sparse_full_diagonal_amd_pattern(
  UF_long n,
  const UF_long *col_ptr,
  const UF_long *row_idx) {
  /* A sparse, nearly full-diagonal graph is only an inexpensive AMD
     proposal.  The actual symbolic below must prove either bounded recurring
     work or the regular symmetric high-work regime before AUTO may skip its
     remaining ordering analyses. */
  return n >= 32768u && n <= 131072u && col_ptr != NULL &&
         n <= UF_long_max / 4u && col_ptr[n] >= 2u * n &&
         col_ptr[n] <= 4u * n &&
         is_large_diagonal_circuit_like_pattern(n, col_ptr, row_idx);
}

static int kls_medium_sparse_full_diagonal_amd_high_work_input(
  UF_long n,
  const UF_long *col_ptr,
  const UF_long *row_idx) {
  int columns_sorted = 1;
  for (UF_long col = 0u; col < n; ++col) {
    const UF_long begin = col_ptr[col];
    const UF_long end = col_ptr[col + 1u];
    if (begin >= end || end - begin > 16u) {
      return 0;
    }
    int has_diagonal = 0;
    for (UF_long p = begin; p < end; ++p) {
      if (row_idx[p] >= n) {
        return 0;
      }
      columns_sorted &= p == begin || row_idx[p - 1u] <= row_idx[p];
      has_diagonal |= row_idx[p] == col;
    }
    if (!has_diagonal) {
      return 0;
    }
  }
  return kls_exact_structural_symmetry(
    n, col_ptr, row_idx, columns_sorted, 16u);
}

/* The low-work branch makes skipping NodeND an analysis-economics decision.
   The separate high-work branch covers bounded-degree symmetric grids: AMD's
   realized one-block symbolic must land in the measured regular-grid band.
   Rejection resumes the ordinary AUTO tournament. */
static int kls_medium_sparse_full_diagonal_amd_symbolic_contract(
  UF_long n,
  const UF_long *col_ptr,
  const UF_long *row_idx,
  const trilinos_klu_l_symbolic *symbolic) {
  if (!is_medium_sparse_full_diagonal_amd_pattern(
        n, col_ptr, row_idx) || symbolic == NULL ||
      (symbolic->structural_rank != KLS_KLU_EMPTY &&
       symbolic->structural_rank != n) ||
      symbolic->nblocks == 0u || symbolic->nblocks > n ||
      symbolic->maxblock > n || 16u * symbolic->maxblock < 7u * n ||
      !(symbolic->lnz > 0.0) || !(symbolic->unz > 0.0) ||
      !(symbolic->est_flops > 0.0)) {
    return 0;
  }
  const double rows = (double)n;
  const double fill = symbolic->lnz + symbolic->unz;
  const int balanced =
    symbolic->lnz <= 2.0 * symbolic->unz &&
    symbolic->unz <= 2.0 * symbolic->lnz;
  const int fragmented_low_work =
    symbolic->do_btf && symbolic->nblocks >= 2u &&
    fill >= 4.0 * rows && fill <= 16.0 * rows &&
    symbolic->est_flops >= 32.0 * rows &&
    symbolic->est_flops <= 256.0 * rows;
  const int single_low_work =
    symbolic->do_btf && symbolic->nblocks == 1u &&
    symbolic->maxblock == n && fill >= 8.0 * rows &&
    fill <= 16.0 * rows &&
    symbolic->est_flops >= 96.0 * rows &&
    symbolic->est_flops <= 256.0 * rows;
  const int regular_symmetric_high_work =
    symbolic->do_btf && symbolic->nblocks == 1u &&
    symbolic->maxblock == n && fill >= 56.0 * rows &&
    fill <= 64.0 * rows &&
    symbolic->est_flops >= 12800.0 * rows &&
    symbolic->est_flops <= 16384.0 * rows &&
    kls_medium_sparse_full_diagonal_amd_high_work_input(
      n, col_ptr, row_idx);
  return balanced && (fragmented_low_work || single_low_work ||
                       regular_symmetric_high_work);
}

/* A many-component matrix with one dominant SCC can have a cheap numeric
   factor even when AMD's column order leaves long, cache-unfriendly solve
   streams.  In the retained-factor H100 regime, COLAMD's unsymmetric order
   is preferable for this bounded low-flop shape (power197k: 0.33 s versus
   1.07 s).  The flop ceiling separates it from the superficially similar
   HTC systems, where COLAMD expands the factor by 5--10x. */

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

enum {
  KLS_ASYMMETRIC_BOUNDED_DEGREE_DIRECT_METIS_NONE = 0,
  KLS_ASYMMETRIC_BOUNDED_DEGREE_DIRECT_METIS_COARSE_FRINGE = 1,
  KLS_ASYMMETRIC_BOUNDED_DEGREE_DIRECT_METIS_THIN_FRINGE = 2
};

/* Candidate analysis is compiled independently of the optional ordering
   backend.  Keep its scoped bookkeeping available even when METIS is not. */
static _Thread_local int
  kls_asymmetric_bounded_degree_direct_metis_analyze_class;
static _Thread_local int
  kls_asymmetric_bounded_degree_direct_metis_analyze_selected;
static _Thread_local int kls_near_symmetric_mega_hub_amd_analyze_path;
static _Thread_local int kls_near_symmetric_mega_hub_amd_analyze_selected;

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
    int has_diagonal = 0;
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
        has_diagonal = 1;
      }
    }
    diagonal_count += (UF_long)has_diagonal;
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

/* The parallel separator-preserving ordering and its tighter CAMD windows
   were measured on the sparse-spike half of the broader class (nxp1,
   average degree below eight).  Dense spikes need the original serial
   NodeNDP/refinement route: rajat30 otherwise grows from roughly 9M factor
   entries to 60M/90M. */
static int is_large_sparse_nearly_diagonal_spiked_metis_pattern(
  UF_long n,
  const UF_long *col_ptr,
  const UF_long *row_idx) {
  return col_ptr != NULL && col_ptr[n] < 8u * n &&
         is_large_nearly_diagonal_spiked_metis_pattern(n, col_ptr, row_idx);
}

static int is_large_moderate_degree_weak_diagonal_metis_pattern(
  UF_long n,
  const UF_long *col_ptr,
  const UF_long *row_idx) {
  if (n < 500000u || col_ptr == NULL || row_idx == NULL ||
      n > UF_long_max / 12u ||
      col_ptr[n] < 6u * n || col_ptr[n] > 12u * n) {
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
    if (col_degree > 4096u) {
      valid = 0;
      break;
    }
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      const UF_long row = row_idx[p];
      if (row >= n || row_degree[row] >= 4096u) {
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

  return max_col_degree <= 4096u && max_row_degree <= 4096u &&
         (max_col_degree >= 64u || max_row_degree >= 64u) &&
         10.0 * (double)diagonal_count >= 2.0 * (double)n &&
         10.0 * (double)diagonal_count <= 6.0 * (double)n;
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

static int symbolic_is_low_work_dominant_btf(
  UF_long n,
  const trilinos_klu_l_symbolic *symbolic);

/* Low-work member of the symmetric partial-diagonal family.  Its ordinary
   AMD/BTF symbolic is already compact enough that the lightweight matched
   PTS lifecycle wins; the heavier SPRAL/direct-value preparation does not.
   Select this route from the cached topology plus measured symbolic shape,
   rather than the former OPF_3754 dimension fingerprint. */

/* A value-aware row match is a natural first operation when almost every
   structural diagonal entry is absent.  The compact matcher stores both row
   and column in one 32-bit word, and the later direct-value map stores input
   positions in 16 bits, so these are representation limits rather than a
   matrix-family fingerprint.  Bounded, nonempty row and column degrees keep
   the sort/augment work linearithmic in a genuinely sparse graph. */

/* A compact, full-rank one-block symbolic with almost no usable structural
   diagonal predicts that an unmatched first factor will be temporary.  When
   its estimated factor is large enough to amortize the value-aware match, run
   that guarded trial before factorization instead of building the unmatched
   numeric and replacing it at the first update. */

/* A bounded-degree symmetric graph with roughly half of its structural
   diagonal missing is a useful proposal for value-aware row matching.  The
   exact reciprocal check makes the rule invariant under simultaneous row /
   column relabeling, while the degree and total-work bounds cap the matching
   cost.  It intentionally says nothing about the eventual numeric engine;
   the matched symbolic and factor must prove those capabilities below. */

/* Cheap analyze-time admission for the heavier lifecycle.  A quarter-million
   input entries is a conservative storage/work proxy for whether its fixed
   match, direct-map, and overlap setup can amortize.  The matched symbolic
   and numeric profiles remain authoritative and impose stricter measured
   fill/work floors before any recurring consumer is enabled. */

/* A near-complete diagonal plus a macroscopic column hub and a measurable
   scalar-column population is a cheap proposal for the shallow, low-work
   BTF/PTS regime.  This is deliberately only an input proposal: the actual
   AMD symbolic below must prove the dominant-core and work geometry before
   any analyze-time shortcut is retained. */

/* The retained AMD proposal must expose a dominant core with an almost
   entirely scalar fringe.  Normalized fill/work and balance bounds keep the
   subtype inside the same low-work PTS resource class as the broader
   many-fringe policy. */

static int kls_sparse_symmetric_compare_uf_long(const void *lhs,
                                                 const void *rhs) {
  const UF_long a = *(const UF_long *)lhs;
  const UF_long b = *(const UF_long *)rhs;
  return a < b ? -1 : (a > b ? 1 : 0);
}

/* Exact multiplicity-aware structural-symmetry proof shared by sparse
   topology proposals.  Sorted CSC uses reciprocal binary searches without
   allocation.  Unsorted input builds and sorts one transpose, then compares
   each original column through a scratch buffer bounded by the proposal's
   already-proved maximum column degree. */
__attribute__((noinline, cold))
static int kls_exact_structural_symmetry(
  UF_long n,
  const UF_long *col_ptr,
  const UF_long *row_idx,
  int columns_sorted,
  UF_long max_column_degree) {
  const UF_long nnz = col_ptr[n];
  if (columns_sorted) {
    for (UF_long col = 0u; col < n; ++col) {
      const UF_long begin = col_ptr[col];
      const UF_long end = col_ptr[col + 1u];
      for (UF_long p = begin; p < end; ++p) {
        const UF_long row = row_idx[p];
        if (row <= col || (p > begin && row_idx[p - 1u] == row)) {
          continue;
        }
        UF_long forward_end = p + 1u;
        while (forward_end < end && row_idx[forward_end] == row) {
          ++forward_end;
        }
        UF_long lo = col_ptr[row];
        UF_long hi = col_ptr[row + 1u];
        while (lo < hi) {
          const UF_long mid = lo + (hi - lo) / 2u;
          if (row_idx[mid] < col) {
            lo = mid + 1u;
          } else {
            hi = mid;
          }
        }
        const UF_long reciprocal_begin = lo;
        hi = col_ptr[row + 1u];
        while (lo < hi) {
          const UF_long mid = lo + (hi - lo) / 2u;
          if (row_idx[mid] <= col) {
            lo = mid + 1u;
          } else {
            hi = mid;
          }
        }
        if (forward_end - p != lo - reciprocal_begin) {
          return 0;
        }
      }
    }
    return 1;
  }

  if ((size_t)n + 1u > SIZE_MAX / sizeof(UF_long) ||
      (size_t)nnz > SIZE_MAX / sizeof(UF_long) ||
      (size_t)max_column_degree > SIZE_MAX / sizeof(UF_long)) {
    return 0;
  }
  UF_long *transpose_ptr =
    (UF_long *)calloc((size_t)n + 1u, sizeof(*transpose_ptr));
  UF_long *transpose_idx =
    (UF_long *)malloc((size_t)nnz * sizeof(*transpose_idx));
  UF_long *next = (UF_long *)malloc((size_t)n * sizeof(*next));
  if (transpose_ptr == NULL || transpose_idx == NULL || next == NULL) {
    free(transpose_ptr);
    free(transpose_idx);
    free(next);
    return 0;
  }
  for (UF_long p = 0u; p < nnz; ++p) {
    transpose_ptr[row_idx[p] + 1u]++;
  }
  for (UF_long row = 0u; row < n; ++row) {
    transpose_ptr[row + 1u] += transpose_ptr[row];
    next[row] = transpose_ptr[row];
  }
  for (UF_long col = 0u; col < n; ++col) {
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      transpose_idx[next[row_idx[p]]++] = col;
    }
  }
  free(next);
  for (UF_long col = 0u; col < n; ++col) {
    const UF_long degree = transpose_ptr[col + 1u] - transpose_ptr[col];
    if (degree > 1u) {
      qsort(transpose_idx + transpose_ptr[col], (size_t)degree,
            sizeof(*transpose_idx), kls_sparse_symmetric_compare_uf_long);
    }
  }

  UF_long *sorted_column = (UF_long *)malloc(
    (size_t)(max_column_degree > 0u ? max_column_degree : 1u) *
      sizeof(*sorted_column));
  if (sorted_column == NULL) {
    free(transpose_ptr);
    free(transpose_idx);
    return 0;
  }
  int symmetric = 1;
  for (UF_long col = 0u; col < n && symmetric; ++col) {
    const UF_long begin = col_ptr[col];
    const UF_long degree = col_ptr[col + 1u] - begin;
    const UF_long transpose_degree =
      transpose_ptr[col + 1u] - transpose_ptr[col];
    if (degree != transpose_degree) {
      symmetric = 0;
      break;
    }
    memcpy(sorted_column, row_idx + begin,
           (size_t)degree * sizeof(*sorted_column));
    if (degree > 1u) {
      qsort(sorted_column, (size_t)degree, sizeof(*sorted_column),
            kls_sparse_symmetric_compare_uf_long);
    }
    symmetric = memcmp(sorted_column, transpose_idx + transpose_ptr[col],
                       (size_t)degree * sizeof(*sorted_column)) == 0;
  }
  free(sorted_column);
  free(transpose_ptr);
  free(transpose_idx);
  return symmetric;
}

/* A full-diagonal, structurally symmetric sparse graph with a small scalar
   SCC fringe and moderate hubs is a proposal for a tuned NodeNDP ordering.
   The ratios describe topology and resource scale rather than a benchmark
   dimension.  The actual METIS symbolic must still prove full rank, the BTF
   core/fringe geometry, bounded fill, and a usable separator decomposition. */

/* Preserve the standard AUTO/8T contract. */

/* Capability contract for the guarded dense reciprocal-hub route.  Explicit
   orientation, ordering, scaling, backend, BTF, and pivot choices retain
   their documented meaning. */

/* A nearly diagonal sparse graph with one or a few moderate spikes can
   overlap matching with NodeND, then build a fixed-pivot predicted numeric.
   Keep the public contract separate from the topology and from the later
   measured factor: no input dimension authorizes the recurring engine. */

/* The matched pattern no longer has the public input coordinates.  Recurring
   dispatch therefore starts from the cached input proposal, but still needs
   the real predicted factor, its separator, and normalized fill/work bounds. */

/* A nearly all-private separator with moderate predicted work benefits from
   the retained cluster forest and the lower relaxed-consume floors.  Denser
   factors and wider separator pipelines keep the generic measured schedule. */

/* Direct NodeNDP is a capability of sparse asymmetric graphs with bounded
   in/out degree and one nearly spanning strongly connected component.  Keep
   the standard AUTO/8T contract explicit. */

/* Allocation-free exact reciprocity check for the rare balanced-triangle,
   unsorted case.  The proposal has already proved degree <= 32, so scanning
   the reciprocal column is bounded.  Multiplicity is compared exactly. */

/* The proposal is invariant under transpose and simultaneous relabeling.
   Degree/diagonal guards make the SCC pass a rare cold-path cost.  Raw SCC
   density distinguishes a thin fringe, which benefits from fourteen NodeNDP
   leaves and narrower CAMD windows, from a coarser fringe that benefits from
   eight leaves.  The retained symbolic and numeric are guarded separately. */

/* Dense, almost-full-diagonal directed graphs with one dominant in/out hub
   can settle into a giant unscaled METIS block whose trailing pivot panel is
   profitable as a pipelined dense tail. */

/* The proposal is invariant under transpose, simultaneous relabeling, and
   stored entry order.  A column-only pass rejects unrelated giant matrices
   before the row-degree allocation and O(nnz) scan. */

/* A dense, almost-full-diagonal graph with nearly balanced directed degrees,
   one macroscopic-but-bounded hub, and a thin disconnected fringe is a
   proposal for direct AMD/BTF analysis and max row scaling. */

/* The input stage is invariant under transpose and simultaneous relabeling.
   It admits only dense, almost-full-diagonal graphs with a macroscopic but
   bounded in/out hub and very small aggregate degree skew.  A nonzero skew
   proves structural asymmetry without relying on stored order or vertex
   labels.  This is deliberately only a cheap proposal; the measured AMD/BTF
   symbolic below proves the dominant-core fringe and factor economics. */

/* Dense inputs can settle into a much sparser collection of moderate BTF
   components.  Their factor geometry, not an input identity, determines
   whether the cooperative row mirror amortizes its retained storage. */

/* Capability contract for giant symmetric scalar-fringe systems. */

/* A giant, exactly structurally symmetric sparse graph with an almost-full
   diagonal, a small missing-diagonal scalar fringe, and one or a few bounded
   hubs is the topology behind the profitable overlapped NodeNDP/row route.
   The bounds are representation and resource ratios, not a benchmark-size
   window.  Topology merely launches the proposal; the retained symbolic and
   measured numeric below independently prove the expensive lifecycle. */

/* Very sparse, full-diagonal graphs in these resource bands can amortize a
   tuned NodeNDP/CAMD refinement and a retained cooperative row numeric.  The
   input band caps the cost of evaluating that proposal; later absolute
   symbolic floors decide whether a graph can actually repay row metadata,
   while measured fill separates the fine and wider giant-matrix schedules.
   Topology only proposes the lifecycle; the retained symbolic and numeric
   must independently prove their fill, work, and separator bounds. */

static int kls_sparse_full_diagonal_metis_row_input_economics(
  UF_long n,
  const UF_long *col_ptr) {
  return col_ptr != NULL && n >= 100000u && n <= 4194304u &&
    n <= UF_long_max / 5u && col_ptr[0] == 0u &&
    col_ptr[n] >= 3u * n && col_ptr[n] <= 5u * n;
}

static int kls_sparse_full_diagonal_metis_row_input_profile(
  UF_long n,
  const UF_long *col_ptr,
  const UF_long *row_idx) {
  return kls_sparse_full_diagonal_metis_row_input_economics(n, col_ptr) &&
    is_large_very_low_degree_full_diagonal_pattern(n, col_ptr, row_idx);
}

/* Classify the actual no-BTF METIS symbolic only after it proves balanced
   intermediate fill and a useful mostly private separator.  Absolute bounds
   express two setup-cost crossovers; normalized bounds keep both meaningful
   across the complete input band.  Class one uses the fine sparse schedule;
   class two retains the requested tolerance and wider giant schedule. */

/* Recurring row dispatch requires the measured fixed-pivot factor to preserve
   the symbolic economics, balanced streams, and normalized arithmetic
   intensity.  The same proof applies to the predicted bootstrap and to a
   later KLU-built replacement.  Cache it after each numeric build. */

/* The raced METIS symbolic must expose one giant SCC plus an almost entirely
   scalar fringe, bounded intermediate fill/work, and a complete mostly
   private separator before any later setup suppression can use the proposal.
   Absolute fill/work floors are the cold-setup crossover; normalized bounds
   make the contract meaningful for extensions and relabelings. */

/* Recurring row updates require the actual unscaled fixed-pivot factor to
   preserve the balanced fill and arithmetic-intensity regime.  Cache this
   verdict after each numeric build; replacements are reclassified. */

/* A successful structural match is still not sufficient to justify the
   specialized repeated-update lifecycle.  Its direct-value map, solve tree,
   and overlapped preparation have fixed setup/coordination costs.  Require a
   one-block AMD symbolic with enough absolute factor storage and arithmetic
   to amortize those costs, in addition to normalized shape bounds.  These
   floors separate measured economics: OPF_10000 is about 796K entries / 7.7M
   predicted operations, while its real OPF_3754 sibling is only 268K / 2.4M
   and is faster on the ordinary lifecycle. */
static int kls_symmetric_partial_diagonal_match_symbolic_profile(
  UF_long n,
  UF_long nnz,
  kls_orientation orientation,
  kls_ordering ordering,
  const trilinos_klu_l_symbolic *symbolic) {
  if (symbolic == NULL || orientation != KLS_ORIENTATION_NORMAL ||
      ordering != KLS_ORDERING_AMD || n < 8192u || n > 131072u ||
      n > UF_long_max / 512u || nnz < 6u * n || nnz > 12u * n ||
      symbolic->do_btf || symbolic->nblocks != 1u ||
      symbolic->maxblock != n ||
      (symbolic->structural_rank != KLS_KLU_EMPTY &&
       symbolic->structural_rank != n) ||
      !(symbolic->lnz > 0.0) || !(symbolic->unz > 0.0) ||
      !(symbolic->est_flops > 0.0)) {
    return 0;
  }
  const double nd = (double)n;
  const double fill = symbolic->lnz + symbolic->unz;
  return fill >= 8.0 * nd && fill <= 24.0 * nd &&
    fill >= 4.0e5 &&
    symbolic->est_flops >= 32.0 * nd &&
    symbolic->est_flops <= 512.0 * nd &&
    symbolic->est_flops >= 4.0e6 &&
    symbolic->lnz <= 2.0 * symbolic->unz &&
    symbolic->unz <= 2.0 * symbolic->lnz;
}

/* Matching is only the proposal.  Specialized PTS/direct-value consumers
   require the retained unscaled AMD single block to prove balanced, bounded
   symbolic and measured work with almost no pivot detours. */
__attribute__((noinline))
static int kls_symmetric_partial_diagonal_match_factor_profile(
  const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      solver->col_ptr == NULL || solver->row_perm == NULL ||
      ((solver->row_scale == NULL) != (solver->col_scale == NULL)) ||
      solver->orientation != KLS_ORIENTATION_NORMAL ||
      solver->stats.selected_ordering != KLS_ORDERING_AMD ||
      solver->common.scale != -1 || solver->numeric->Rs != NULL ||
      solver->numeric_is_predicted ||
      solver->col_ptr[solver->n] != solver->nnz ||
      !kls_symmetric_partial_diagonal_match_symbolic_profile(
        solver->n, solver->nnz, solver->orientation,
        solver->stats.selected_ordering, solver->symbolic) ||
      solver->common.noffdiag == KLS_KLU_EMPTY ||
      solver->common.noffdiag > solver->n / 512u + 16u ||
      solver->pivot_nudge_count != 0u ||
      solver->common.kls_perturb_count != 0u ||
      !(solver->common.flops > 0.0)) {
    return 0;
  }
  const int generic_contract =
    getenv("KLS_DISABLE_SYMMETRIC_PARTIAL_DIAGONAL_MATCH_POLICY") == NULL &&
    kls_repeated_update_workload(&solver->options) &&
    solver->options.orientation == KLS_ORIENTATION_AUTO &&
    solver->options.ordering == KLS_ORDERING_AUTO &&
    solver->options.scale == KLS_SCALE_AUTO &&
    solver->options.backend != KLS_BACKEND_SERIAL &&
    solver->options.threads >= 2 && solver->options.use_btf &&
    solver->options.static_pivoting;
  if (!generic_contract) {
    return 0;
  }
  const double n = (double)solver->n;
  const double numeric_fill =
    (double)solver->numeric->lnz + (double)solver->numeric->unz;
  return numeric_fill >= 8.0 * n && numeric_fill <= 24.0 * n &&
    numeric_fill >= 4.0e5 &&
    solver->common.flops >= 32.0 * n &&
    solver->common.flops <= 512.0 * n &&
    solver->common.flops >= 3.0e6 &&
    (double)solver->numeric->lnz <=
      2.0 * (double)solver->numeric->unz &&
    (double)solver->numeric->unz <=
      2.0 * (double)solver->numeric->lnz;
}

static int kls_symmetric_partial_diagonal_match_factor_cycle(
  const kls_solver *solver) {
  return kls_symmetric_partial_diagonal_match_factor_profile(solver);
}

static int choose_auto_scale_from_pattern_impl(
    UF_long n,
    const UF_long *col_ptr,
    const UF_long *row_idx,
    const kls_options *options,
    const double *numeric_values,
    int *unscaled_value_certified_out,
    int *unscaled_trial_certified_out) {
  if (unscaled_value_certified_out != NULL) {
    *unscaled_value_certified_out = 0;
  }
  if (unscaled_trial_certified_out != NULL) {
    *unscaled_trial_certified_out = 0;
  }
  if (options == NULL || col_ptr == NULL || row_idx == NULL || numeric_values == NULL ||
      options->scale != KLS_SCALE_AUTO || n <= 0) {
    return options == NULL ? 2 : initial_scale(options);
  }

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

  if (100.0 * (double)diag_count >= 99.0 * (double)n &&
      100.0 * (double)missing_diagonal <= (double)n &&
      50.0 * (double)weak_diagonal <= (double)n) {
    /* Nearly every row already has a numerically strong diagonal relative
       to its own largest entry.  Absolute magnitude spread alone is not a
       reason to pay for KLU's row-scale construction and a second complete
       factor: pivot inflation can still authorize the ordinary realized
       scale trial below.  This is a current-value dominance certificate,
       independent of dimensions, sparsity, or matrix identity. */
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
    if (getenv("KLS_TRACE_AUTO_SCALE_STATS") != NULL) {
      fprintf(stderr,
              "KLS auto scale stats: diag=%ld/%ld weak=%ld missing=%ld "
              "diag-spread=%.4e row-spread=%.4e p90/p10=%.4e "
              "row-unit=%.4f diag-unit=%.4f\n",
              (long)diag_count, (long)n, (long)weak_diagonal,
              (long)missing_diagonal,
              min_diag > 0.0 && min_diag < DBL_MAX
                ? max_diag / min_diag : DBL_MAX,
              row_min > 0.0 ? row_maximum / row_min : DBL_MAX,
              row_p10 > 0.0 ? row_p90 / row_p10 : DBL_MAX,
              row_unit_fraction, diag_unit_fraction);
    }
    if (100.0 * (double)diag_count >= 95.0 * (double)n &&
        20.0 * (double)missing_diagonal <= (double)n &&
        10.0 * (double)weak_diagonal <= (double)n &&
        row_p10 > 0.0 && row_p90 / row_p10 >= 128.0) {
      /* A nearly complete, mostly usable diagonal paired with at least
         seven bits of central row-magnitude range does not need a row-scale
         pass to expose pivots.  Starting unscaled avoids constructing and
         applying scale vectors over every retained update.  This is a
         current-value certificate: a realized fill/work expansion may still
         authorize the ordinary measured scale comparison below. */
      if (unscaled_value_certified_out != NULL) {
        *unscaled_value_certified_out = 1;
      }
      free(diag_max);
      free(row_max);
      return -1;
    }
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
    if (100.0 * (double)diag_count >= 95.0 * (double)n &&
        20.0 * (double)missing_diagonal <= (double)n &&
        10.0 * (double)weak_diagonal <= (double)n &&
        row_min > 0.0 && row_p10 > 0.0 &&
        row_maximum / row_min >= 1.0e8 &&
        row_p90 / row_p10 < 128.0) {
      /* Most rows have usable diagonal pivots and the central row population
         fits within seven magnitude bits, but a small value tail spans at
         least 27 bits.  Sum scaling confines those exceptional rows without
         flattening the ordinary diagonal population.  This complements the
         wide-central-range unscaled certificate above and depends only on
         the current values, not the matrix dimensions or source family. */
      if (unscaled_trial_certified_out != NULL &&
          row_p90 / row_p10 <= row_p90_p10_limit) {
        /* The central 80% of rows already fits in the same balanced value
           range used by the ordinary generic scale chooser.  Only the
           extreme tail motivates scaling, so a long repeated lifecycle may
           compare the unscaled specialized executor.  Realized fill/work,
           pivoting, and a complete residual still arbitrate adoption. */
        *unscaled_trial_certified_out = 1;
      }
      free(diag_max);
      free(row_max);
      return 1;
    }
    if ((double)diag_count < min_diag_fraction * (double)n &&
        row_min > 0.0 && row_p10 > 0.0 &&
        row_maximum / row_min <= row_spread_limit &&
        row_p90 / row_p10 <= row_p90_p10_limit) {
      free(diag_max);
      free(row_max);
      /* Both KLU row norms preserve the balanced value range certified
         above.  A repeated numeric lifecycle favors max scaling: it avoids
         the serial reduction dependency in every scale refresh, while the
         one-shot path retains sum scaling's historical default.  Numeric
         pivot/fill and conditioning gates still arbitrate any later scale
         replacement. */
      return kls_repeated_update_workload(options) ? 2 : 1;
    }
  }

  free(diag_max);
  free(row_max);
  return initial_scale(options);
}

static int choose_auto_scale_from_pattern(UF_long n,
                                          const UF_long *col_ptr,
                                          const UF_long *row_idx,
                                          const kls_options *options,
                                          const double *numeric_values) {
  return choose_auto_scale_from_pattern_impl(
    n, col_ptr, row_idx, options, numeric_values, NULL, NULL);
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

static int choose_auto_scale_from_values(kls_solver *solver,
                                         const double *numeric_values) {
  if (solver == NULL) {
    return 2;
  }
  /* The generic first choice depends only on the current values.  Any
     expensive factor can subsequently compare every KLU scale mode using
     realized fill, work, pivoting, and condition estimates. */
  int value_certified = 0;
  int trial_certified = 0;
  const int selected = choose_auto_scale_from_pattern_impl(
    solver->n, solver->col_ptr, solver->row_idx, &solver->options,
    numeric_values, &value_certified, &trial_certified);
  solver->auto_scale_value_certified = value_certified;
  solver->auto_scale_unscaled_trial_certified = trial_certified;
  return selected;
}

/* The input match is only a proposal.  Recurring compact row/value/solve
   machinery is admitted after the real numeric proves low measured work,
   balanced streams, a bounded pivot detour count, and no perturbation. */

/* Admit sparse, low-work dominant BTF factors with only a bounded fixed-pivot
   fringe just below the generic EGraph crossover.  A shallow largest-block
   elimination tree excludes long skinny factors whose synchronization depth
   overwhelms the same nominal work.  The decision depends on retained factor
   geometry and numeric work, not on a public matrix fingerprint. */

/* Close the fixed absolute-work gap immediately below the large dominant-BTF
   EGraph window.  These factors are effectively one high-work block plus a
   handful of scalar blocks.  Bounded work per row and per retained entry
   separate the sparse irregular class from dense-group factors whose mapped
   local walk is much cheaper.  The final admission also requires at least
   80M measured dependency work in a broad, shallow graph after the schedule
   has been built.  None of these tests use input nonzero counts or an exact
   dimension/option fingerprint. */

/* A moderately fragmented factor with one bounded core can keep its complete
   paired row streams in the compact 16-bit representation.  In that regime,
   duplicate value publication, residual SpMVs, and speculative numerics cost
   as much as the retained factor/solve pair.  Select the capability from the
   adopted symbolic and numeric resources rather than one input fingerprint. */

/* A compact triangular-solve mirror can collapse a long consecutive run of
   singleton BTF blocks into one backward loop.  Select that representation
   from its actual block geometry: the run table costs two bytes per block,
   so require both broad coverage and a long enough run to amortize its build
   over repeated solves.  The bounds are representation/cost limits rather
   than an input dimension or nonzero fingerprint. */
static int kls_compact_singleton_run_solve_profile(
  const kls_solver *solver,
  UF_long *singleton_blocks_out,
  UF_long *max_run_out) {
  if (singleton_blocks_out != NULL) {
    *singleton_blocks_out = 0u;
  }
  if (max_run_out != NULL) {
    *max_run_out = 0u;
  }
  if (solver == NULL || solver->symbolic == NULL ||
      solver->symbolic->R == NULL || !solver->symbolic->do_btf ||
      solver->n == 0u || solver->n > (UF_long)UINT16_MAX ||
      solver->symbolic->nblocks < 512u ||
      solver->symbolic->nblocks > (UF_long)UINT16_MAX ||
      solver->symbolic->R[0] != 0u ||
      solver->symbolic->R[solver->symbolic->nblocks] != solver->n ||
      getenv("KLS_DISABLE_COMPACT_SINGLETON_RUN_SOLVE") != NULL) {
    return 0;
  }

  const UF_long nblocks = solver->symbolic->nblocks;
  UF_long singleton_blocks = 0u;
  UF_long run = 0u;
  UF_long max_run = 0u;
  for (UF_long block = 0u; block < nblocks; ++block) {
    const UF_long begin = solver->symbolic->R[block];
    const UF_long end = solver->symbolic->R[block + 1u];
    if (end <= begin || end > solver->n) {
      return 0;
    }
    if (end == begin + 1u) {
      singleton_blocks++;
      run++;
      if (run > max_run) {
        max_run = run;
      }
    } else {
      run = 0u;
    }
  }
  if (singleton_blocks_out != NULL) {
    *singleton_blocks_out = singleton_blocks;
  }
  if (max_run_out != NULL) {
    *max_run_out = max_run;
  }
  return singleton_blocks * 4u >= nblocks * 3u && max_run >= 256u;
}

/* A sparse full-rank matrix split into thousands of genuinely tiny SCCs has
   less total numeric work than the global ordering, matching, and tolerance
   trials that generic medium matrices can amortize.  Select that resource
   regime from normalized BTF geometry and symbolic work rather than the
   dimensions of one operating-point family. */

/* Recurring solve and accuracy decisions additionally require the compact,
   full-precision factor that justified skipping the speculative trials.
   Value changes that leave this bounded regime lose the shortcut. */

/* A sparse, high-work min-fill factor with a tiny mixed-component BTF fringe
   can expose a useful subtree solve even though its EGraph numeric update
   remains the recurring refactor winner.  Admit that retained resource state
   from normalized symbolic and numeric evidence: the factor must be balanced,
   full-rank, unscaled, pivot-stable, and have both a 1/512--1/256 fringe and
   2,048--4,096 operations per row.  The PTS builder independently validates
   its actual forest before any solve is routed through it. */

/* A moderately expensive, fragmented dominant-BTF factor amortizes a live
   worker crew and wide subtree scheduling without entering the high-fill
   regime where those fixed choices lose.  Use measured numeric work, fill,
   pivot stability, and normalized SCC shape rather than an input identity.
   Explicit normal/AMD/unscaled requests share the same retained factor
   capabilities as AUTO after those choices have actually been selected. */

/* A high-fill compact factor with a deep, moderately wide dependency graph
   keeps the scalar lean crew useful across the short solve/update gap.  The
   normalized fill, work, depth, and width boundary covers both the power-grid
   case that motivated the policy and an independently generated sparse
   expander, while rejecting equally sized shallow/low-fill counterexamples. */

/* Reusing a residual verdict depends on the factor epoch and numeric frame,
   not on a matrix fingerprint.  A solve already covered by an ordinary
   self-check first has to pass it.  Only a byte-identical RHS against that
   unchanged, plain-fp64 factor can then reuse the verdict; factor/refactor
   entry invalidates it.  The selected triangular solve still runs in full;
   only its redundant residual check is reused. */
static int kls_verified_rhs_reuse_capable(const kls_solver *solver) {
  return solver != NULL && solver->symbolic != NULL &&
    solver->numeric != NULL && solver->col_ptr != NULL &&
    getenv("KLS_DISABLE_VERIFIED_RHS_REUSE") == NULL &&
    solver->orientation == KLS_ORIENTATION_NORMAL &&
    solver->common.scale <= 0 && solver->numeric->Rs == NULL &&
    solver->row_perm == NULL && solver->user_col_perm == NULL &&
    solver->row_scale == NULL && solver->col_scale == NULL &&
    !solver->numeric_is_predicted &&
    !solver->fp32_last_used && !solver->numeric_needs_refinement &&
    !solver->tight_tol_refine && !solver->solve_recovery_active &&
    solver->pivot_nudge_count == 0u &&
    solver->common.kls_perturb_count == 0u &&
    solver->common.tol >= 1.0e-6 &&
    /* Either the retained row numeric or the generic factor contract must
       already require a residual.  Reuse never creates a new accuracy
       shortcut for a factor that was not being checked. */
    (solver->row_solve_self_check || solver->solve_contract_probe == 2) &&
    (solver->solve_refine_values != NULL || solver->values != NULL);
}

/* A bounded-degree, structurally symmetric core can carry a small collection
   of missing-diagonal degree-one columns.  BTF matching turns that surface
   into scalar components around a dominant core.  In this normalized regime
   transpose AMD plus the persistent lean row walk repays its setup over a
   repeated factor/solve horizon.  This proposal deliberately uses topology,
   not matrix order or a narrow nnz fingerprint; the real AMD symbolic and
   numeric must pass the independent guards below before recurring kernels
   are selected. */

/* Kernel admission is stricter than the topology proposal.  Numeric scaling,
   pivoting, fill, and work must still describe the balanced retained factor
   on which the lean lifecycle was measured. */

/* A moderate-order AMD single block with high estimated work can amortize the
   retained EGraph/PTS lifecycle.  This symbolic half only records the
   selected candidate; it does not suppress speculative analysis because
   matrices with identical structure can have very different pivot behavior.
   Every recurring kernel below still requires the measured factor contract.
   The bounds are normalized capabilities and broad resource limits, not a
   matrix dimension window. */

/* Recurring lifecycle choices require a materially pivoted, balanced actual
   factor.  Symbolic admission alone never selects the settled kernels. */

/* A low-work dominant BTF factor with a mostly scalar fringe can expose a
   useful PTS forest even when the symbolic flop estimate is pessimistic.
   Admit the retained execution regime rather than a matrix-size fingerprint;
   every consumer below still requires the measured numeric bounds. */

/* A high-work AMD factor with a tiny, nearly scalar BTF fringe can amortize
   predicted first-factor construction and the settled EGraph lifecycle.
   Record that execution regime from normalized symbolic state; values still
   have to pass the measured factor contract below before recurring kernels
   are selected. */

/* Recurring decisions require an unscaled, balanced numeric with bounded
   pivoting and measured fill/work.  A value-only change can therefore reject
   the lifecycle even when it shares an admitted symbolic pattern. */

/* Count preferred symbolic diagonals that are rejected at the requested
   threshold but retained at 1e-4.  This is the numerical reason a tighter
   threshold can remove off-diagonal pivot detours; dimensions and raw nnz do
   not label the matrix.  The surrounding symbolic bounds limit both the scan
   and the amount of unrelated factor policy inherited by a proposal. */
static int kls_value_tolerance_crossing_initial_profile(
  const kls_solver *solver,
  const double *numeric_values) {
  if (solver == NULL || solver->symbolic == NULL || numeric_values == NULL ||
      solver->col_ptr == NULL || solver->row_idx == NULL ||
      solver->symbolic->P == NULL || solver->symbolic->Q == NULL ||
      solver->symbolic->R == NULL ||
      getenv("KLS_DISABLE_VALUE_TOLERANCE_CROSSING_POLICY") != NULL ||
      fabs(solver->options.pivot_tolerance - 0.001) > 1.0e-12 ||
      solver->options.ordering != KLS_ORDERING_AUTO ||
      solver->options.scale != KLS_SCALE_AUTO || solver->common.scale != -1 ||
      solver->orientation != KLS_ORIENTATION_NORMAL ||
      solver->stats.selected_ordering != KLS_ORDERING_AMD ||
      solver->n < 256u || solver->n > (UF_long)UINT16_MAX ||
      solver->n > UF_long_max / 32u ||
      solver->nnz < 2u * solver->n || solver->nnz > 32u * solver->n ||
      solver->symbolic->nblocks < 2u ||
      solver->symbolic->maxblock >= solver->n ||
      solver->symbolic->R[0] != 0u ||
      solver->symbolic->R[solver->symbolic->nblocks] != solver->n ||
      4u * solver->symbolic->maxblock < 3u * solver->n ||
      solver->symbolic->nblocks >
        solver->n - solver->symbolic->maxblock + 1u ||
      !(solver->symbolic->lnz > 0.0) ||
      !(solver->symbolic->unz > 0.0) ||
      !(solver->symbolic->est_flops > 0.0)) {
    return 0;
  }
  const double rows = (double)solver->n;
  const double fill = solver->symbolic->lnz + solver->symbolic->unz;
  if (fill < 4.0 * rows || fill > 32.0 * rows ||
      solver->symbolic->est_flops > 512.0 * rows ||
      solver->symbolic->lnz > 2.0 * solver->symbolic->unz ||
      solver->symbolic->unz > 2.0 * solver->symbolic->lnz) {
    return 0;
  }

  UF_long *pinv =
    (UF_long *)malloc((size_t)solver->n * sizeof(*pinv));
  if (pinv == NULL) {
    return 0;
  }
  for (UF_long k = 0u; k < solver->n; ++k) {
    pinv[k] = KLS_KLU_EMPTY;
  }
  int valid = 1;
  for (UF_long k = 0u; k < solver->n; ++k) {
    const UF_long row = solver->symbolic->P[k];
    if (row >= solver->n || pinv[row] != KLS_KLU_EMPTY) {
      valid = 0;
      break;
    }
    pinv[row] = k;
  }

  UF_long crossings = 0u;
  UF_long block = 0u;
  for (UF_long k = 0u; valid && k < solver->n; ++k) {
    while (block + 1u < solver->symbolic->nblocks &&
           solver->symbolic->R[block + 1u] <= k) {
      block++;
    }
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    const UF_long oldcol = solver->symbolic->Q[k];
    const UF_long olddiagrow = solver->symbolic->P[k];
    if (k < k1 || k >= k2 || oldcol >= solver->n ||
        olddiagrow >= solver->n) {
      valid = 0;
      break;
    }
    const UF_long begin = solver->col_ptr[oldcol];
    const UF_long end = solver->col_ptr[oldcol + 1u];
    if (begin > end || end > solver->nnz) {
      valid = 0;
      break;
    }
    double diagonal = 0.0;
    double column_max = 0.0;
    for (UF_long p = begin; p < end; ++p) {
      const UF_long oldrow = solver->row_idx[p];
      if (oldrow >= solver->n) {
        valid = 0;
        break;
      }
      const double value = fabs(numeric_values[p]);
      if (!isfinite(value)) {
        valid = 0;
        break;
      }
      const UF_long newrow = pinv[oldrow];
      if (newrow >= k1 && newrow < k2 && value > column_max) {
        column_max = value;
      }
      if (oldrow == olddiagrow && value > diagonal) {
        diagonal = value;
      }
    }
    if (valid && column_max > 0.0 && diagonal >= 1.0e-4 * column_max &&
        diagonal < 1.0e-3 * column_max) {
      crossings++;
    }
  }
  free(pinv);

  const UF_long required =
    solver->n / 64u + (solver->n % 64u != 0u);
  const int accepted = valid && crossings >= 16u && crossings >= required;
  if (getenv("KLS_TRACE_INITIAL_TOLERANCE_CROSSINGS") != NULL) {
    fprintf(stderr,
            "KLS initial tolerance crossings: rows=%ld crossings=%ld "
            "required=%ld accepted=%d\n",
            (long)solver->n, (long)crossings, (long)required, accepted);
  }
  return accepted;
}

/* Propose a tighter initial threshold from the values that the selected
   symbolic diagonal will actually see.  This is deliberately a numeric
   pressure census rather than an input-shape classifier: the candidate is
   available only for a caller-declared repeated lifecycle, and the ordinary
   factor/refactor residual contract remains attached to every threshold
   selected below the caller's request. */
static double kls_generic_initial_pivot_pressure_tolerance(
  const kls_solver *solver,
  const double *numeric_values) {
  if (solver == NULL || solver->symbolic == NULL || numeric_values == NULL ||
      solver->col_ptr == NULL || solver->row_idx == NULL ||
      solver->symbolic->P == NULL || solver->symbolic->Q == NULL ||
      solver->symbolic->R == NULL || solver->row_perm != NULL ||
      !solver->options.static_pivoting ||
      !kls_repeated_update_workload(&solver->options) ||
      solver->options.ordering != KLS_ORDERING_AUTO ||
      fabs(solver->options.pivot_tolerance - 0.001) > 1.0e-12 ||
      getenv("KLS_DISABLE_GENERIC_INITIAL_PIVOT_PRESSURE") != NULL) {
    return solver == NULL ? 0.001 : solver->options.pivot_tolerance;
  }

  UF_long *pinv =
    (UF_long *)malloc((size_t)solver->n * sizeof(*pinv));
  double *row_norm = NULL;
  if (pinv == NULL && solver->n > 0u) {
    return solver->options.pivot_tolerance;
  }
  for (UF_long k = 0u; k < solver->n; ++k) {
    pinv[k] = KLS_KLU_EMPTY;
  }

  int valid = 1;
  for (UF_long k = 0u; k < solver->n; ++k) {
    const UF_long row = solver->symbolic->P[k];
    if (row >= solver->n || pinv[row] != KLS_KLU_EMPTY) {
      valid = 0;
      break;
    }
    pinv[row] = k;
  }
  if (solver->common.scale > 0 && valid) {
    row_norm = (double *)calloc((size_t)solver->n, sizeof(*row_norm));
    if (row_norm == NULL && solver->n > 0u) {
      valid = 0;
    }
    for (UF_long col = 0u; valid && col < solver->n; ++col) {
      const UF_long begin = solver->col_ptr[col];
      const UF_long end = solver->col_ptr[col + 1u];
      if (begin > end || end > solver->nnz) {
        valid = 0;
        break;
      }
      for (UF_long p = begin; p < end; ++p) {
        const UF_long row = solver->row_idx[p];
        const double value = fabs(numeric_values[p]);
        if (row >= solver->n || !isfinite(value)) {
          valid = 0;
          break;
        }
        if (solver->common.scale == 1) {
          row_norm[row] += value;
          if (!isfinite(row_norm[row])) {
            valid = 0;
            break;
          }
        } else if (value > row_norm[row]) {
          row_norm[row] = value;
        }
      }
    }
  }

  UF_long pressured = 0u;
  UF_long pressured_positive = 0u;
  UF_long recovered_1e4 = 0u;
  UF_long recovered_1e6 = 0u;
  UF_long recovered_3e8 = 0u;
  UF_long recovered_1e8 = 0u;
  UF_long block = 0u;
  for (UF_long k = 0u; valid && k < solver->n; ++k) {
    while (block + 1u < solver->symbolic->nblocks &&
           solver->symbolic->R[block + 1u] <= k) {
      block++;
    }
    if (block >= solver->symbolic->nblocks) {
      valid = 0;
      break;
    }
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    const UF_long oldcol = solver->symbolic->Q[k];
    const UF_long diagrow = solver->symbolic->P[k];
    if (k < k1 || k >= k2 || oldcol >= solver->n ||
        diagrow >= solver->n) {
      valid = 0;
      break;
    }
    const UF_long begin = solver->col_ptr[oldcol];
    const UF_long end = solver->col_ptr[oldcol + 1u];
    if (begin > end || end > solver->nnz) {
      valid = 0;
      break;
    }
    double diagonal = 0.0;
    double column_max = 0.0;
    for (UF_long p = begin; p < end; ++p) {
      const UF_long oldrow = solver->row_idx[p];
      if (oldrow >= solver->n) {
        valid = 0;
        break;
      }
      const double denom = row_norm != NULL && row_norm[oldrow] > 0.0
        ? row_norm[oldrow] : 1.0;
      const double value = fabs(numeric_values[p]) / denom;
      if (!isfinite(value)) {
        valid = 0;
        break;
      }
      const UF_long newrow = pinv[oldrow];
      if (newrow >= k1 && newrow < k2 && value > column_max) {
        column_max = value;
      }
      if (oldrow == diagrow && value > diagonal) {
        diagonal = value;
      }
    }
    if (!valid || !(column_max > 0.0) ||
        diagonal >= solver->options.pivot_tolerance * column_max) {
      continue;
    }
    pressured++;
    if (!(diagonal > 0.0)) {
      continue;
    }
    pressured_positive++;
    const double ratio = diagonal / column_max;
    recovered_1e4 += ratio >= 1.0e-4;
    recovered_1e6 += ratio >= 1.0e-6;
    recovered_3e8 += ratio >= 3.0e-8;
    recovered_1e8 += ratio >= 1.0e-8;
  }

  free(row_norm);
  free(pinv);

  const UF_long fraction_required =
    solver->n / 128u + (solver->n % 128u != 0u);
  const UF_long required = fraction_required > 16u
    ? fraction_required : 16u;
  double selected = solver->options.pivot_tolerance;
#define KLS_PRESSURE_RECOVERS(count)                                    \
  ((count) >= required && pressured > 0u && pressured_positive > 0u && \
   4u * (count) >= 3u * pressured && 5u * (count) >= 4u * pressured_positive)
  /* The raw column ratio precedes elimination.  Leave two decades of
     headroom for Schur-complement updates so a value sitting on the census
     boundary is not immediately rejected by the actual factor walk. */
  if (valid && KLS_PRESSURE_RECOVERS(recovered_1e4)) {
    selected = 1.0e-6;
  } else if (valid && KLS_PRESSURE_RECOVERS(recovered_3e8)) {
    /* Use the strongest directly measured tight line that retains the
       pressure cohort.  This modest separation above the 1e-8 floor can
       materially improve the raw triangular solve without changing the
       representation contract; every selected numeric still carries the
       same residual verification and recovery path. */
    selected = 3.0e-8;
  } else if (valid && KLS_PRESSURE_RECOVERS(recovered_1e8)) {
    /* The 1e-4 census could not retain this pressure cohort.  Use the
       lowest measured line that retains at least three quarters of all
       pressured columns instead of interpolating an unmeasured threshold
       that can change pivot layout across later numeric generations.
       Solve-side residual certification remains active for every such
       structural selection. */
    selected = 1.0e-8;
  }
#undef KLS_PRESSURE_RECOVERS
  if (getenv("KLS_TRACE_INITIAL_PIVOT_PRESSURE") != NULL) {
    fprintf(stderr,
            "KLS initial pivot pressure: valid=%d rows=%ld pressure=%ld/%ld "
            "recover=%ld,%ld,%ld,%ld required=%ld scale=%ld selected=%.3g\n",
            valid, (long)solver->n, (long)pressured,
            (long)pressured_positive, (long)recovered_1e4,
            (long)recovered_1e6, (long)recovered_3e8,
            (long)recovered_1e8, (long)required,
            (long)solver->common.scale, selected);
  }
  return selected;
}

static int kls_uses_structural_initial_pivot_tolerance(
  const kls_solver *solver) {
  {
    /* A generic minimum-fill ordering may select a value-pressure proposal
       before its first numeric.  Treat that measured proposal like the other
       guarded tight-tolerance numerics so solve-side residual checks remain
       active; this depends on the selected algorithm and retained tolerance,
       not on a matrix family or dimension window. */
    return solver != NULL &&
      kls_repeated_update_workload(&solver->options) &&
      (solver->stats.selected_ordering == KLS_ORDERING_AMF ||
       solver->stats.selected_ordering == KLS_ORDERING_AMMF) &&
      solver->common.tol < solver->options.pivot_tolerance;
  }
  return (solver != NULL && solver->value_tolerance_crossing_cycle);
}

/* Keep this policy table out of factor_impl.  Its individual structural
   clauses evolve independently of the numeric kernels; inlining the table
   makes deleting a clause shift later hot code and creates layout-sensitive
   performance changes with no numeric-policy change. */
__attribute__((noinline))
static double choose_initial_auto_pivot_tolerance(
  kls_solver *solver,
  const double *numeric_values) {
  if (solver == NULL || solver->symbolic == NULL ||
      fabs(solver->options.pivot_tolerance - 0.001) > 1.0e-12) {
    return solver == NULL ? 0.001 : solver->options.pivot_tolerance;
  }
  solver->value_tolerance_crossing_cycle = 0;

  {
    /* Minimum-fill orderings deliberately commit to a near-diagonal pivot
       structure.  On a repeated lifecycle, admit the value-pressure census
       before building the first numeric so a clearly separated pivot regime
       does not require a throwaway loose-threshold factor.  AMD and AMF3 stay
       on the caller's threshold: their elimination layouts can lose locality
       even when a tighter factor has less fill (the measured rajat15 case).
       The census itself requires broad recovery of pressured rows, and the
       solve-side contract above verifies every selected tight numeric. */
    const int generic_minfill_pressure =
      kls_repeated_update_workload(&solver->options) &&
      (solver->stats.selected_ordering == KLS_ORDERING_AMF ||
       solver->stats.selected_ordering == KLS_ORDERING_AMMF) &&
      getenv("KLS_DISABLE_GENERIC_MINFILL_PIVOT_PRESSURE") == NULL;
    return (generic_minfill_pressure ||
            getenv("KLS_ENABLE_INITIAL_PIVOT_PRESSURE") != NULL)
      ? kls_generic_initial_pivot_pressure_tolerance(solver, numeric_values)
      : solver->options.pivot_tolerance;
  }

  if (kls_compact_amf_two_block_initial_cycle(solver)) {
    /* A modestly stronger threshold keeps this entire normalized compact
       class inside its 16-bit factor streams.  The solve-side residual gate
       measures every changed numeric and corrects any raw solve outside the
       strict relative-L2 margin. */
    return 3.0e-8;
  }

  if (kls_value_tolerance_crossing_initial_profile(
        solver, numeric_values)) {
    /* A material set of symbolic diagonal choices lies between the requested
       and retained thresholds.  Keep those pivots directly and attach the
       ordinary structural-tolerance solve contract to the resulting factor. */
    solver->value_tolerance_crossing_cycle = 1;
    return 1.0e-4;
  }

  if (solver->options.ordering == KLS_ORDERING_AUTO &&
      solver->stats.selected_ordering == KLS_ORDERING_METIS &&
      is_large_sparse_diagonal_low_degree_pattern(solver->n, solver->col_ptr,
                                                  solver->row_idx)) {
    /* Ordinary entrywise updates still need more pivot headroom than 1e-5:
       rare generations exhaust the bounded refinement loop.  The 1e-4
       factor stays residual-valid while retaining enough of the compact
       numeric to beat the competing H100 cycle. */
    return 1.0e-4;
  }
  if (solver->options.ordering == KLS_ORDERING_AUTO &&
      (solver->stats.selected_ordering == KLS_ORDERING_METIS ||
       solver->stats.selected_ordering == KLS_ORDERING_AMF) &&
      solver->symbolic->do_btf && solver->symbolic->nblocks <= 4u &&
      (double)solver->symbolic->maxblock >= 0.95 * (double)solver->n &&
      ((solver->common.scale == 1 &&
        is_medium_spiked_low_diagonal_pattern(solver->n, solver->col_ptr,
                                              solver->row_idx)) ||
       (solver->common.scale == 0 &&
        is_small_spiked_low_diagonal_pattern(solver->n, solver->col_ptr,
                                             solver->row_idx)))) {
    /* Deficiency orderings explode under loose threshold pivoting on
       these classes (TSOPF: first factor 9.2e5 fill at 1e-4 vs 6.4e4 at
       1e-8); they are selected for their tight-tolerance factor and the
       config-derived refine gate (tol < 1e-6) polices the accuracy. */
    return solver->stats.selected_ordering == KLS_ORDERING_AMF ? 1.0e-8
                                                               : 1.0e-4;
  }

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
  /* The primary BTF pass already computes an exact maximum matching.  A
     second, equally valid matching is only an optional ordering portfolio:
     building it for every orientation and then sometimes factoring both can
     cost more than a complete repeated lifecycle.  Keep that experiment for
     explicit A/B requests; generic AUTO instead spends
     its portfolio budget on the retained ordering/partition candidates. */
  common->kls_btf_match_trial =
    getenv("KLS_DISABLE_BTF_MATCH_TRIAL") == NULL &&
    getenv("KLS_ENABLE_BTF_MATCH_TRIAL") != NULL;
  {
    /* bounded-effort BTF matching: mac_econ's 89%-missing diagonal
       drives Duff's maxtrans into ~9s of augmenting-path work (of a
       10.5s analysis); a work cap yields an incomplete matching the
       analyze already supports. */
    const char *mw = getenv("KLS_BTF_MAXWORK");
    if (mw != NULL && mw[0] != '\0') {
      common->maxwork = atof(mw);
    }
  }
  return KLS_OK;
}

/* The medium partial-diagonal path is identified before its row permutation
   replaces the original pattern.  Carry that decision through matching and
   METIS ordering without re-identifying the now-full matched diagonal. */
static _Thread_local int kls_medium_partial_static_metis_ctx;

/* Generic nested-dissection portfolio variants.  Context five is the uniform
   two-leaves-per-worker, coarse-refinement proposal.  Context six retains the
   same forest but spends more CAMD work on a finer resource-scaled window when
   the incumbent's realized symbolic economics say factor-memory traffic, not
   arithmetic, dominates its repeated lifecycle.  Unlike the archived
   pre-static contexts, these affect ordering only; they never change matching,
   scaling, or numeric acceptance policy.  The scale hint is computed from
   current values by the ordinary AUTO rule; it lets the refinement budget
   account for the recurring row-scale stream without recognizing an input
   shape. */
static _Thread_local int kls_generic_metis_tuning_ctx;
/* A late AUTO challenge has already observed the incumbent's realized
   numeric cost.  Let that narrowly scoped trial spend modestly more ordering
   work to reduce every subsequent refactor; explicit METIS and analyze-time
   portfolio candidates retain their established defaults. */
static _Thread_local int kls_auto_metis_promotion_quality;
static _Thread_local int kls_generic_metis_scale_hint;
/* A symbolic portfolio may compare separator widths on the identical graph.
   Keep the override thread-local so concurrent orientation/order trials do
   not change one another's resource choice. */
static _Thread_local int kls_metis_ndp_npes_override;

/* The raw one-block PTS class is identified before row matching rewrites
   the diagonal.  Carry its deterministic matching policy through the
   weight sort and augmentation helpers. */
static _Thread_local int kls_low_work_symmetric_partial_diagonal_pts_ctx;

/* Large weak-diagonal matrices are identified before weighted row matching
   makes their diagonal structurally complete.  Preserve that raw-pattern
   verdict through the matched METIS analyze so its NodeNDP leaf count can be
   selected from the workload that the retained numeric factor will execute. */
static _Thread_local int kls_large_weak_diagonal_static_metis_ctx;

typedef struct kls_metis_separator_capture {
  UF_long order_call_index;
  kls_separator_analysis separator;
} kls_metis_separator_capture;

typedef struct kls_metis_order_context {
  UF_long n;
  idx_t npes;
  kls_separator_analysis *separator;
  UF_long order_call_count;
  UF_long capture_count;
  UF_long capture_capacity;
  kls_metis_separator_capture *captures;
} kls_metis_order_context;

static void kls_metis_order_context_clear(
  kls_metis_order_context *context) {
  if (context == NULL) {
    return;
  }
  for (UF_long pos = 0; pos < context->capture_count; ++pos) {
    kls_separator_analysis_clear(&context->captures[pos].separator);
  }
  free(context->captures);
  context->captures = NULL;
  context->capture_count = 0;
  context->capture_capacity = 0;
}

static int kls_metis_order_context_append_separator(
  kls_metis_order_context *context,
  UF_long order_call_index,
  kls_separator_analysis *separator) {
  if (context == NULL || separator == NULL || separator->n == 0u) {
    return 0;
  }
  if (context->capture_count == context->capture_capacity) {
    UF_long next_capacity =
      context->capture_capacity == 0u ? 4u : 2u * context->capture_capacity;
    if (next_capacity <= context->capture_capacity ||
        next_capacity > (UF_long)(SIZE_MAX / sizeof(*context->captures))) {
      return 0;
    }
    kls_metis_separator_capture *next =
      (kls_metis_separator_capture *)realloc(
        context->captures, (size_t)next_capacity * sizeof(*next));
    if (next == NULL) {
      return 0;
    }
    for (UF_long pos = context->capture_capacity; pos < next_capacity; ++pos) {
      memset(&next[pos], 0, sizeof(next[pos]));
    }
    context->captures = next;
    context->capture_capacity = next_capacity;
  }
  kls_metis_separator_capture *capture =
    &context->captures[context->capture_count++];
  capture->order_call_index = order_call_index;
  kls_separator_analysis_move(&capture->separator, separator);
  return 1;
}

static const kls_separator_analysis *kls_metis_order_context_capture_for_call(
  const kls_metis_order_context *context,
  UF_long order_call_index) {
  if (context == NULL) {
    return NULL;
  }
  for (UF_long pos = 0; pos < context->capture_count; ++pos) {
    if (context->captures[pos].order_call_index == order_call_index) {
      return &context->captures[pos].separator;
    }
  }
  return NULL;
}

static void kls_metis_order_context_move_largest_separator(
  kls_metis_order_context *context,
  kls_separator_analysis *separator_out) {
  if (context == NULL || separator_out == NULL) {
    return;
  }
  UF_long best_pos = KLS_KLU_EMPTY;
  UF_long best_n = 0;
  for (UF_long pos = 0; pos < context->capture_count; ++pos) {
    if (context->captures[pos].separator.n > best_n) {
      best_n = context->captures[pos].separator.n;
      best_pos = pos;
    }
  }
  if (best_pos != KLS_KLU_EMPTY) {
    kls_separator_analysis_move(separator_out,
                                &context->captures[best_pos].separator);
  }
}

static int compare_idx_t(const void *a, const void *b) {
  const idx_t left = *(const idx_t *)a;
  const idx_t right = *(const idx_t *)b;
  return (left > right) - (left < right);
}

static int metis_size_ok(UF_long n, idx_t slots_per_entry) {
  return n >= 0 && (uint64_t)n <= (uint64_t)(SIZE_MAX / (size_t)slots_per_entry);
}

static idx_t kls_metis_ndp_leaf_count(UF_long n,
                                      int thread_count,
                                      int subtree_threshold) {
  if (n < KLS_METIS_NDP_MIN_ROWS || thread_count <= 1) {
    return 0;
  }

  UF_long leaves = (UF_long)thread_count;
  if (subtree_threshold && n >= KLS_METIS_NDP_SUBTREE_THRESHOLD_MIN_ROWS) {
    /*
     * The KLS-first separator executor consumes SubtreeLU-style queues.
     * SubtreeLU partitions until components are smaller than
     * max(200, N/1000) after the minimum log2(P) depth.  METIS_NodeNDP only
     * records a fixed top tree, so ask for enough leaves to retain that
     * threshold-sized tree while the numeric executor still uses thread_count.
     * Keep the 200-row saturated regime on the old thread-count tree for now:
     * the pivoting first-factor executor is not yet robust enough for very
     * fine separator pivot scopes on smaller synthetic weak-pivot grids.
     */
    UF_long target_rows =
      (n + KLS_METIS_NDP_TARGET_DIVISOR - 1u) /
      KLS_METIS_NDP_TARGET_DIVISOR;
    if (target_rows < KLS_METIS_NDP_MIN_LEAF_ROWS) {
      target_rows = KLS_METIS_NDP_MIN_LEAF_ROWS;
    }
    if (target_rows > 0u) {
      UF_long target_leaves = (n + target_rows - 1u) / target_rows;
      if (target_leaves > KLS_METIS_NDP_MAX_OVERPARTITION_LEAVES) {
        target_leaves = KLS_METIS_NDP_MAX_OVERPARTITION_LEAVES;
      }
      if (target_leaves > leaves) {
        leaves = target_leaves;
      }
    }
  }

  if (leaves > n) {
    leaves = n;
  }
  const UF_long max_recorded_npes = (UF_long)((IDX_MAX - 1) / 2);
  if (leaves > max_recorded_npes) {
    leaves = max_recorded_npes;
  }
  return leaves > 1u ? (idx_t)leaves : 0;
}

static int kls_metis_append_separator_component(
  idx_t npes,
  idx_t cpos,
  const idx_t *sizes,
  UF_long *component_ptr,
  unsigned char *component_kind,
  UF_long *component_left_child,
  UF_long *component_right_child,
  UF_long *component_parent,
  UF_long *tree_component,
  UF_long *component_pos,
  UF_long *private_components,
  UF_long *pipeline_components,
  UF_long *private_rows,
  UF_long *pipeline_rows,
  UF_long *private_max_rows,
  UF_long *pipeline_max_rows) {
  if (npes <= 0 || cpos < 0 || sizes == NULL || component_ptr == NULL ||
      component_kind == NULL || component_pos == NULL ||
      private_components == NULL || pipeline_components == NULL ||
      private_rows == NULL || pipeline_rows == NULL ||
      private_max_rows == NULL || pipeline_max_rows == NULL ||
      component_left_child == NULL || component_right_child == NULL ||
      component_parent == NULL || tree_component == NULL) {
    return 0;
  }
  const idx_t component_count = 2 * npes - 1;
  if (cpos >= component_count) {
    return 0;
  }

  if (cpos < npes - 1) {
    if (!kls_metis_append_separator_component(
          npes, 2 * cpos + 2, sizes, component_ptr, component_kind,
          component_left_child, component_right_child, component_parent,
          tree_component, component_pos, private_components, pipeline_components,
          private_rows, pipeline_rows, private_max_rows,
          pipeline_max_rows) ||
        !kls_metis_append_separator_component(
          npes, 2 * cpos + 1, sizes, component_ptr, component_kind,
          component_left_child, component_right_child, component_parent,
          tree_component, component_pos, private_components, pipeline_components,
          private_rows, pipeline_rows, private_max_rows,
          pipeline_max_rows)) {
      return 0;
    }
  }

  const idx_t size_index = component_count - 1 - cpos;
  if (size_index < 0 || size_index >= component_count ||
      sizes[size_index] < 0) {
    return 0;
  }
  const UF_long size = (UF_long)sizes[size_index];
  const UF_long pos = *component_pos;
  if (pos >= (UF_long)component_count ||
      component_ptr[pos] > UF_long_max - size) {
    return 0;
  }
  tree_component[cpos] = pos;
  component_left_child[pos] = KLS_KLU_EMPTY;
  component_right_child[pos] = KLS_KLU_EMPTY;
  component_parent[pos] = KLS_KLU_EMPTY;
  if (cpos < npes - 1) {
    const UF_long right_child = tree_component[2 * cpos + 2];
    const UF_long left_child = tree_component[2 * cpos + 1];
    if (right_child >= (UF_long)component_count ||
        left_child >= (UF_long)component_count) {
      return 0;
    }
    component_left_child[pos] = left_child;
    component_right_child[pos] = right_child;
    component_parent[left_child] = pos;
    component_parent[right_child] = pos;
  }
  component_ptr[pos + 1u] = component_ptr[pos] + size;
  if (cpos >= npes - 1) {
    component_kind[pos] = 0u;
    (*private_components)++;
    *private_rows += size;
    if (size > *private_max_rows) {
      *private_max_rows = size;
    }
  } else {
    component_kind[pos] = 1u;
    (*pipeline_components)++;
    *pipeline_rows += size;
    if (size > *pipeline_max_rows) {
      *pipeline_max_rows = size;
    }
  }
  *component_pos = pos + 1u;
  return 1;
}

static int kls_build_metis_separator_analysis(
  UF_long n,
  idx_t npes,
  const idx_t *sizes,
  const idx_t *metis_iperm,
  const UF_long *final_perm,
  kls_separator_analysis *separator_out) {
  if (separator_out == NULL) {
    return 0;
  }
  kls_separator_analysis_clear(separator_out);
  if (n == 0u || npes <= 1 || sizes == NULL || metis_iperm == NULL ||
      final_perm == NULL || (uint64_t)(2 * npes - 1) > (uint64_t)UINT_MAX) {
    return 0;
  }
  const UF_long component_count = (UF_long)(2 * npes - 1);
  if (component_count > (UF_long)(SIZE_MAX / sizeof(UF_long)) - 1u ||
      n > (UF_long)(SIZE_MAX / sizeof(unsigned int))) {
    return 0;
  }

  UF_long *component_ptr =
    (UF_long *)calloc((size_t)component_count + 1u, sizeof(*component_ptr));
  unsigned char *component_kind =
    (unsigned char *)calloc((size_t)component_count, sizeof(*component_kind));
  unsigned int *order_component =
    (unsigned int *)malloc((size_t)n * sizeof(*order_component));
  UF_long *component_left_child =
    (UF_long *)malloc((size_t)component_count * sizeof(*component_left_child));
  UF_long *component_right_child =
    (UF_long *)malloc((size_t)component_count * sizeof(*component_right_child));
  UF_long *component_parent =
    (UF_long *)malloc((size_t)component_count * sizeof(*component_parent));
  UF_long *tree_component =
    (UF_long *)malloc((size_t)component_count * sizeof(*tree_component));
  unsigned int *metis_pos_component =
    (unsigned int *)malloc((size_t)n * sizeof(*metis_pos_component));
  if (component_ptr == NULL || component_kind == NULL ||
      order_component == NULL || component_left_child == NULL ||
      component_right_child == NULL || component_parent == NULL ||
      tree_component == NULL || metis_pos_component == NULL) {
    free(component_ptr);
    free(component_kind);
    free(order_component);
    free(component_left_child);
    free(component_right_child);
    free(component_parent);
    free(tree_component);
    free(metis_pos_component);
    return 0;
  }
  for (UF_long component = 0; component < component_count; ++component) {
    component_left_child[component] = KLS_KLU_EMPTY;
    component_right_child[component] = KLS_KLU_EMPTY;
    component_parent[component] = KLS_KLU_EMPTY;
    tree_component[component] = KLS_KLU_EMPTY;
  }

  UF_long component_pos = 0;
  UF_long private_components = 0;
  UF_long pipeline_components = 0;
  UF_long private_rows = 0;
  UF_long pipeline_rows = 0;
  UF_long private_max_rows = 0;
  UF_long pipeline_max_rows = 0;
  if (!kls_metis_append_separator_component(
        npes, 0, sizes, component_ptr, component_kind,
        component_left_child, component_right_child, component_parent,
        tree_component, &component_pos,
        &private_components, &pipeline_components,
        &private_rows, &pipeline_rows, &private_max_rows, &pipeline_max_rows) ||
      component_pos != component_count || component_ptr[component_count] != n) {
    free(component_ptr);
    free(component_kind);
    free(order_component);
    free(component_left_child);
    free(component_right_child);
    free(component_parent);
    free(tree_component);
    free(metis_pos_component);
    return 0;
  }
  free(tree_component);

  for (UF_long component = 0; component < component_count; ++component) {
    const UF_long begin = component_ptr[component];
    const UF_long end = component_ptr[component + 1u];
    if (begin > end || end > n) {
      free(component_ptr);
      free(component_kind);
      free(order_component);
      free(component_left_child);
      free(component_right_child);
      free(component_parent);
      free(metis_pos_component);
      return 0;
    }
    for (UF_long pos = begin; pos < end; ++pos) {
      metis_pos_component[pos] = (unsigned int)component;
    }
  }

  for (UF_long pos = 0; pos < n; ++pos) {
    const UF_long vertex = final_perm[pos];
    if (vertex >= n || metis_iperm[vertex] < 0 ||
        (uint64_t)metis_iperm[vertex] >= (uint64_t)n) {
      free(component_ptr);
      free(component_kind);
      free(order_component);
      free(component_left_child);
      free(component_right_child);
      free(component_parent);
      free(metis_pos_component);
      return 0;
    }
    order_component[pos] = metis_pos_component[metis_iperm[vertex]];
  }
  free(metis_pos_component);

  separator_out->n = n;
  separator_out->thread_count = (UF_long)npes;
  separator_out->component_count = component_count;
  separator_out->private_component_count = private_components;
  separator_out->pipeline_component_count = pipeline_components;
  separator_out->private_rows = private_rows;
  separator_out->pipeline_rows = pipeline_rows;
  separator_out->private_max_rows = private_max_rows;
  separator_out->pipeline_max_rows = pipeline_max_rows;
  separator_out->component_ptr = component_ptr;
  separator_out->component_kind = component_kind;
  separator_out->order_component = order_component;
  separator_out->component_left_child = component_left_child;
  separator_out->component_right_child = component_right_child;
  separator_out->component_parent = component_parent;
  return 1;
}

static int kls_build_metis_separator_forest(
  const trilinos_klu_l_symbolic *symbolic,
  const kls_metis_order_context *context,
  kls_separator_analysis *separator_out) {
  if (separator_out != NULL) {
    kls_separator_analysis_clear(separator_out);
  }
  if (symbolic == NULL || context == NULL || separator_out == NULL ||
      symbolic->R == NULL || symbolic->n == 0u ||
      symbolic->nblocks == 0u || context->capture_count == 0u) {
    return 0;
  }

  const UF_long n = symbolic->n;
  const UF_long nblocks = symbolic->nblocks;
  UF_long component_count = 0;
  UF_long captured_blocks = 0;
  UF_long order_call_index = 0;
  for (UF_long block = 0; block < nblocks; ++block) {
    const UF_long begin = symbolic->R[block];
    const UF_long end = symbolic->R[block + 1u];
    if (end < begin || end > n) {
      return 0;
    }
    const UF_long block_n = end - begin;
    const kls_separator_analysis *capture = NULL;
    if (block_n > 3u) {
      capture = kls_metis_order_context_capture_for_call(
        context, order_call_index++);
      if (capture != NULL && capture->n != block_n) {
        capture = NULL;
      }
    }
    const UF_long add_components =
      capture != NULL ? capture->component_count : 1u;
    if (add_components == 0u ||
        component_count > UF_long_max - add_components) {
      return 0;
    }
    component_count += add_components;
    if (capture != NULL) {
      captured_blocks++;
    }
  }
  if (captured_blocks == 0u) {
    return 0;
  }
  if (nblocks > 1u) {
    if (component_count > UF_long_max - (nblocks - 1u)) {
      return 0;
    }
    component_count += nblocks - 1u;
  }
  if (component_count == 0u ||
      component_count > (UF_long)UINT_MAX ||
      component_count > (UF_long)(SIZE_MAX / sizeof(UF_long)) - 1u ||
      n > (UF_long)(SIZE_MAX / sizeof(unsigned int)) ||
      nblocks > (UF_long)(SIZE_MAX / sizeof(UF_long))) {
    return 0;
  }

  kls_separator_analysis forest;
  memset(&forest, 0, sizeof(forest));
  forest.component_ptr =
    (UF_long *)calloc((size_t)component_count + 1u,
                      sizeof(*forest.component_ptr));
  forest.component_kind =
    (unsigned char *)calloc((size_t)component_count,
                            sizeof(*forest.component_kind));
  forest.order_component =
    (unsigned int *)malloc((size_t)n * sizeof(*forest.order_component));
  forest.component_left_child =
    (UF_long *)malloc((size_t)component_count *
                      sizeof(*forest.component_left_child));
  forest.component_right_child =
    (UF_long *)malloc((size_t)component_count *
                      sizeof(*forest.component_right_child));
  forest.component_parent =
    (UF_long *)malloc((size_t)component_count *
                      sizeof(*forest.component_parent));
  UF_long *block_roots =
    (UF_long *)malloc((size_t)nblocks * sizeof(*block_roots));
  if (forest.component_ptr == NULL ||
      forest.component_kind == NULL ||
      forest.order_component == NULL ||
      forest.component_left_child == NULL ||
      forest.component_right_child == NULL ||
      forest.component_parent == NULL ||
      block_roots == NULL) {
    free(block_roots);
    kls_separator_analysis_clear(&forest);
    return 0;
  }
  for (UF_long component = 0; component < component_count; ++component) {
    forest.component_left_child[component] = KLS_KLU_EMPTY;
    forest.component_right_child[component] = KLS_KLU_EMPTY;
    forest.component_parent[component] = KLS_KLU_EMPTY;
  }

  UF_long component_pos = 0;
  order_call_index = 0;
  for (UF_long block = 0; block < nblocks; ++block) {
    const UF_long begin = symbolic->R[block];
    const UF_long end = symbolic->R[block + 1u];
    const UF_long block_n = end - begin;
    const kls_separator_analysis *capture = NULL;
    if (block_n > 3u) {
      capture = kls_metis_order_context_capture_for_call(
        context, order_call_index++);
      if (capture != NULL && capture->n != block_n) {
        capture = NULL;
      }
    }

    if (capture != NULL) {
      const UF_long offset = component_pos;
      for (UF_long component = 0; component < capture->component_count;
           ++component) {
        const UF_long global_component = offset + component;
        if (global_component >= component_count ||
            capture->component_ptr == NULL ||
            capture->component_kind == NULL ||
            capture->order_component == NULL ||
            capture->component_left_child == NULL ||
            capture->component_right_child == NULL ||
            capture->component_parent == NULL ||
            capture->component_ptr[component] >
              capture->component_ptr[component + 1u] ||
            capture->component_ptr[component + 1u] > capture->n) {
          free(block_roots);
          kls_separator_analysis_clear(&forest);
          return 0;
        }
        const UF_long component_size =
          capture->component_ptr[component + 1u] -
          capture->component_ptr[component];
        if (forest.component_ptr[global_component] >
            UF_long_max - component_size) {
          free(block_roots);
          kls_separator_analysis_clear(&forest);
          return 0;
        }
        forest.component_ptr[global_component + 1u] =
          forest.component_ptr[global_component] + component_size;
        forest.component_kind[global_component] =
          capture->component_kind[component];
        const UF_long left = capture->component_left_child[component];
        const UF_long right = capture->component_right_child[component];
        const UF_long parent = capture->component_parent[component];
        if ((left != KLS_KLU_EMPTY && left >= capture->component_count) ||
            (right != KLS_KLU_EMPTY && right >= capture->component_count) ||
            (parent != KLS_KLU_EMPTY && parent >= capture->component_count)) {
          free(block_roots);
          kls_separator_analysis_clear(&forest);
          return 0;
        }
        forest.component_left_child[global_component] =
          left == KLS_KLU_EMPTY ? KLS_KLU_EMPTY : offset + left;
        forest.component_right_child[global_component] =
          right == KLS_KLU_EMPTY ? KLS_KLU_EMPTY : offset + right;
        forest.component_parent[global_component] =
          parent == KLS_KLU_EMPTY ? KLS_KLU_EMPTY : offset + parent;
      }
      for (UF_long row = 0; row < block_n; ++row) {
        const unsigned int local_component = capture->order_component[row];
        if ((UF_long)local_component >= capture->component_count ||
            offset + (UF_long)local_component > (UF_long)UINT_MAX) {
          free(block_roots);
          kls_separator_analysis_clear(&forest);
          return 0;
        }
        forest.order_component[begin + row] =
          (unsigned int)(offset + (UF_long)local_component);
      }
      if (forest.private_component_count >
            UF_long_max - capture->private_component_count ||
          forest.pipeline_component_count >
            UF_long_max - capture->pipeline_component_count ||
          forest.private_rows > UF_long_max - capture->private_rows ||
          forest.pipeline_rows > UF_long_max - capture->pipeline_rows) {
        free(block_roots);
        kls_separator_analysis_clear(&forest);
        return 0;
      }
      forest.private_component_count += capture->private_component_count;
      forest.pipeline_component_count += capture->pipeline_component_count;
      forest.private_rows += capture->private_rows;
      forest.pipeline_rows += capture->pipeline_rows;
      if (capture->private_max_rows > forest.private_max_rows) {
        forest.private_max_rows = capture->private_max_rows;
      }
      if (capture->pipeline_max_rows > forest.pipeline_max_rows) {
        forest.pipeline_max_rows = capture->pipeline_max_rows;
      }
      block_roots[block] = offset + capture->component_count - 1u;
      component_pos += capture->component_count;
    } else {
      if (component_pos >= component_count ||
          forest.component_ptr[component_pos] > UF_long_max - block_n) {
        free(block_roots);
        kls_separator_analysis_clear(&forest);
        return 0;
      }
      forest.component_kind[component_pos] = 0u;
      forest.component_ptr[component_pos + 1u] =
        forest.component_ptr[component_pos] + block_n;
      if (component_pos > (UF_long)UINT_MAX) {
        free(block_roots);
        kls_separator_analysis_clear(&forest);
        return 0;
      }
      for (UF_long row = begin; row < end; ++row) {
        forest.order_component[row] = (unsigned int)component_pos;
      }
      if (forest.private_component_count == UF_long_max ||
          forest.private_rows > UF_long_max - block_n) {
        free(block_roots);
        kls_separator_analysis_clear(&forest);
        return 0;
      }
      forest.private_component_count++;
      forest.private_rows += block_n;
      if (block_n > forest.private_max_rows) {
        forest.private_max_rows = block_n;
      }
      block_roots[block] = component_pos;
      component_pos++;
    }
  }
  if (forest.component_ptr[component_pos] != n) {
    free(block_roots);
    kls_separator_analysis_clear(&forest);
    return 0;
  }

  UF_long root = block_roots[0];
  for (UF_long block = 1u; block < nblocks; ++block) {
    if (component_pos >= component_count ||
        root >= component_pos ||
        block_roots[block] >= component_pos) {
      free(block_roots);
      kls_separator_analysis_clear(&forest);
      return 0;
    }
    const UF_long parent = component_pos;
    forest.component_kind[parent] = 1u;
    forest.component_left_child[parent] = root;
    forest.component_right_child[parent] = block_roots[block];
    forest.component_parent[root] = parent;
    forest.component_parent[block_roots[block]] = parent;
    forest.component_ptr[parent + 1u] = forest.component_ptr[parent];
    forest.pipeline_component_count++;
    root = parent;
    component_pos++;
  }
  free(block_roots);
  if (component_pos != component_count ||
      forest.component_ptr[component_count] != n) {
    kls_separator_analysis_clear(&forest);
    return 0;
  }

  forest.n = n;
  forest.global_begin = 0;
  forest.global_end = n;
  forest.global_range_valid = 1;
  forest.thread_count = context->npes > 0 ? (UF_long)context->npes : 0u;
  forest.component_count = component_count;
  kls_separator_analysis_move(separator_out, &forest);
  return 1;
}

/* Exact fill count of a PINNED elimination order: build the strictly-
   lower pattern of P(A+A')P', its elimination tree (Liu, ancestor
   compression), then count fill by row-subtree walks - O(fill) with a
   runaway cap.  Replaces the fully-constrained CAMD call whose degree
   machinery re-derives this same fill and grinds for tens of seconds
   when the pinned order is poor (ss1's 109K block: 19-23s; the lnz
   estimate is that call's only product).  Returns the L fill excluding
   the diagonal (CAMD's LNZ convention), or 0 on allocation failure. */
static UF_long kls_pinned_order_fill_count(UF_long n,
                                           const UF_long *col_ptr,
                                           const UF_long *row_idx,
                                           const idx_t *metis_perm,
                                           UF_long *perm_out) {
  const UF_long nnz = col_ptr[n];
  UF_long *rank = (UF_long *)malloc((size_t)n * sizeof(*rank));
  UF_long *rp = (UF_long *)calloc((size_t)n + 1u, sizeof(*rp));
  UF_long *rj = (UF_long *)malloc(((size_t)nnz + 1u) * sizeof(*rj));
  UF_long *parent = (UF_long *)malloc((size_t)n * sizeof(*parent));
  UF_long *ancestor = (UF_long *)malloc((size_t)n * sizeof(*ancestor));
  UF_long *mark = (UF_long *)malloc((size_t)n * sizeof(*mark));
  if (rank == NULL || rp == NULL || rj == NULL || parent == NULL ||
      ancestor == NULL || mark == NULL) {
    free(rank); free(rp); free(rj);
    free(parent); free(ancestor); free(mark);
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
    free(rank); free(rp); free(rj);
    free(parent); free(ancestor); free(mark);
    return 0;
  }
  /* strictly-lower rows of P(A+A')P': entry (r,c) lands in row
     max(rank) with column min(rank) */
  for (UF_long c = 0; c < n; ++c) {
    for (UF_long p = col_ptr[c]; p < col_ptr[c + 1u]; ++p) {
      const UF_long pr = rank[row_idx[p]];
      const UF_long pc = rank[c];
      if (pr != pc) {
        rp[(pr > pc ? pr : pc) + 1u]++;
      }
    }
  }
  for (UF_long i = 0; i < n; ++i) {
    rp[i + 1u] += rp[i];
  }
  {
    UF_long *cur = (UF_long *)malloc((size_t)n * sizeof(*cur));
    if (cur == NULL) {
      free(rank); free(rp); free(rj);
      free(parent); free(ancestor); free(mark);
      return 0;
    }
    memcpy(cur, rp, (size_t)n * sizeof(*cur));
    for (UF_long c = 0; c < n; ++c) {
      for (UF_long p = col_ptr[c]; p < col_ptr[c + 1u]; ++p) {
        const UF_long pr = rank[row_idx[p]];
        const UF_long pc = rank[c];
        if (pr != pc) {
          const UF_long i = pr > pc ? pr : pc;
          rj[cur[i]++] = pr > pc ? pc : pr;
        }
      }
    }
    free(cur);
  }
  /* Liu etree with ancestor path compression */
  for (UF_long i = 0; i < n; ++i) {
    parent[i] = (UF_long)-1;
    ancestor[i] = (UF_long)-1;
    mark[i] = (UF_long)-1;
    for (UF_long p = rp[i]; p < rp[i + 1u]; ++p) {
      UF_long k = rj[p];
      while (k != (UF_long)-1 && k < i) {
        const UF_long next = ancestor[k];
        ancestor[k] = i;
        if (next == (UF_long)-1) {
          parent[k] = i;
          break;
        }
        k = next;
      }
    }
  }
  /* fill count via row-subtree walks (dedup per row with mark).  The
     walk is O(fill): cap it at a few multiples of nnz - the caller
     only needs "comparable to the competition or garbage", and dense-
     fill pinned orders (TSOPF's fused-dense rows: fill ~n^2/2) would
     otherwise grind for seconds counting exactly how bad they are. */
  const UF_long cap = 8u * (nnz > (UF_long)n ? nnz : (UF_long)n);
  UF_long fill = 0;
  for (UF_long i = 0; i < n && fill < cap; ++i) {
    for (UF_long p = rp[i]; p < rp[i + 1u]; ++p) {
      UF_long k = rj[p];
      while (k < i && mark[k] != i) {
        mark[k] = i;
        fill++;
        k = parent[k];
        if (k == (UF_long)-1) {
          break;
        }
      }
    }
  }
  free(rank); free(rp); free(rj);
  free(parent); free(ancestor); free(mark);
  for (UF_long k = 0; k < n; ++k) {
    perm_out[k] = (UF_long)metis_perm[k];
  }
  return fill > 0 ? fill : 1;
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

/* CKTSO-style incomplete-ND refinement: keep every recorded NodeNDP leaf
   and separator component in postorder, but let constrained AMD reorder all
   vertices within that component.  Fixed-size rank windows approximate this
   idea poorly: they split natural separator components and over-constrain
   large leaves.  The NodeNDP sizes array gives the exact top separator tree,
   so reuse the same component reconstruction used by KLS's subtree
   scheduler and feed those component ids to CAMD. */
static UF_long kls_metis_refine_components_with_camd(
  UF_long n,
  UF_long *col_ptr,
  UF_long *row_idx,
  const idx_t *metis_perm,
  const idx_t *metis_iperm,
  idx_t npes,
  const idx_t *sizes,
  UF_long *perm_out) {
  if (n == 0u || metis_perm == NULL || metis_iperm == NULL || npes <= 1 ||
      sizes == NULL || perm_out == NULL) {
    return 0u;
  }
  UF_long *metis_order =
    (UF_long *)malloc((size_t)n * sizeof(*metis_order));
  UF_long *constraints =
    (UF_long *)malloc((size_t)n * sizeof(*constraints));
  UF_long *camd_perm =
    (UF_long *)malloc((size_t)n * sizeof(*camd_perm));
  if (metis_order == NULL || constraints == NULL || camd_perm == NULL) {
    free(metis_order);
    free(constraints);
    free(camd_perm);
    return 0u;
  }
  int valid = 1;
  for (UF_long pos = 0u; pos < n; ++pos) {
    const idx_t vertex = metis_perm[pos];
    if (vertex < 0 || (uint64_t)vertex >= (uint64_t)n) {
      valid = 0;
      break;
    }
    metis_order[pos] = (UF_long)vertex;
  }
  kls_separator_analysis separator;
  memset(&separator, 0, sizeof(separator));
  if (!valid ||
      !kls_build_metis_separator_analysis(n, npes, sizes, metis_iperm,
                                          metis_order, &separator) ||
      separator.order_component == NULL) {
    kls_separator_analysis_clear(&separator);
    free(metis_order);
    free(constraints);
    free(camd_perm);
    return 0u;
  }
  for (UF_long pos = 0u; pos < n; ++pos) {
    const UF_long vertex = metis_order[pos];
    constraints[vertex] = (UF_long)separator.order_component[pos];
  }

  double control[TRILINOS_CAMD_CONTROL];
  double info[TRILINOS_CAMD_INFO];
  trilinos_camd_l_defaults(control);
  const UF_long result = trilinos_camd_l_order(
    n, col_ptr, row_idx, camd_perm, control, info, constraints);
  if (result >= TRILINOS_CAMD_OK) {
    memcpy(perm_out, camd_perm, (size_t)n * sizeof(*perm_out));
  }
  const double lnz = info[TRILINOS_CAMD_LNZ];
  kls_separator_analysis_clear(&separator);
  free(metis_order);
  free(constraints);
  free(camd_perm);
  if (result < TRILINOS_CAMD_OK) {
    return 0u;
  }
  return lnz > 0.0 && lnz < (double)UF_long_max - (double)n
    ? (UF_long)lnz + n
    : (UF_long)-1;
}

static UF_long kls_metis_order_inner(UF_long n, UF_long *col_ptr,
                                     UF_long *row_idx, UF_long *perm_out,
                                     trilinos_klu_l_common *common);

/* Top-level pattern context for METIS options.  KLU invokes the ordering
   callback on individual BTF blocks, so the callback's local n/nnz cannot
   recognize a whole-matrix class.  Scope this mode around direct AUTO
   selections that need measured option sets (1: fragmented low-degree,
   2: dense reciprocal mega-hub, 3: sparse symmetric scalar fringe). */
static _Thread_local int kls_fragmented_metis_tuning_ctx;

static UF_long kls_metis_order(UF_long n,
                               UF_long *col_ptr,
                               UF_long *row_idx,
                               UF_long *perm_out,
                               trilinos_klu_l_common *common) {
  const double t0 = kls_now_seconds();
  UF_long r = kls_metis_order_inner(n, col_ptr, row_idx, perm_out, common);
  if (getenv("KLS_TRACE_ANALYZE_STAGES") != NULL) {
    fprintf(stderr, "KLS metis: user_order total n=%ld %.3fs\n",
            (long)n, kls_now_seconds() - t0);
  }
  return r;
}

/* The retained proposal must expose a nearly spanning BTF core whose small
   fringe is substantially represented by separate components.  Symbolic fill
   and separator bounds keep the tuned eight-leaf ordering inside its intended
   resource and parallel-work regime. */

/* Scale suppression additionally requires the actual unscaled, no-pivot
   numeric to remain in the balanced fill/work band.  Cache this measured
   verdict after each factor; a replacement invalidates it. */

/* Recurring kernel choices require the measured numeric to stay within the
   normalized fill/work regime that motivated retaining the NodeNDP order.
   A structurally accepted candidate outside this range remains a valid METIS
   factor but uses the ordinary adaptive kernels. */

/* Scoped to the generic bounded-degree/full-diagonal AUTO proposal.  One is
   the fine sparse resource class; two is its wider giant-graph class. */
static _Thread_local int kls_sparse_full_diagonal_metis_row_ctx;

/* Thread budget supplied by analyze_with_ordering.  The CAMD refinement
   policy is defined before the optional mt-metis bridge but consumes the
   same resource count. */
static _Thread_local int kls_metis_order_threads;

/* A sparse full-diagonal block with a moderate hub and balanced triangular
   structure can benefit from a fine CAMD pass over the NodeND ranks even when
   it reached explicit METIS without an enclosing AUTO topology proposal.
   Recognize the callback's realized block rather than an outer matrix order
   or density coordinate. */
__attribute__((noinline, cold))
static int kls_sparse_full_diagonal_hubbed_metis_core_refine_profile(
  UF_long n,
  const UF_long *col_ptr,
  const UF_long *row_idx) {
  if (col_ptr == NULL || row_idx == NULL || n < 32768u || n > 262144u ||
      n > UF_long_max / 8u || col_ptr[0] != 0u ||
      col_ptr[n] < 4u * n || col_ptr[n] > 8u * n) {
    return 0;
  }
  const UF_long nnz = col_ptr[n];
  UF_long diagonal_columns = 0u;
  UF_long scalar_columns = 0u;
  UF_long max_column_degree = 0u;
  UF_long lower_entries = 0u;
  UF_long upper_entries = 0u;
  for (UF_long col = 0u; col < n; ++col) {
    const UF_long begin = col_ptr[col];
    const UF_long end = col_ptr[col + 1u];
    if (begin >= end || end > nnz || end - begin > 512u) {
      return 0;
    }
    const UF_long degree = end - begin;
    scalar_columns += degree == 1u;
    if (degree > max_column_degree) {
      max_column_degree = degree;
    }
    int has_diagonal = 0;
    for (UF_long p = begin; p < end; ++p) {
      const UF_long row = row_idx[p];
      if (row >= n) {
        return 0;
      }
      has_diagonal |= row == col;
      lower_entries += row > col;
      upper_entries += row < col;
    }
    diagonal_columns += (UF_long)has_diagonal;
  }
  if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
    fprintf(stderr,
            "KLS metis: hubbed core proposal n=%ld nnz=%ld diag=%ld "
            "scalar=%ld maxdeg=%ld lower=%ld upper=%ld\n",
            (long)n, (long)nnz, (long)diagonal_columns,
            (long)scalar_columns, (long)max_column_degree,
            (long)lower_entries, (long)upper_entries);
  }
  if (diagonal_columns != n || 1024u * scalar_columns > n ||
      512u * max_column_degree < n || 64u * max_column_degree > n ||
      4.0 * (double)lower_entries > 5.0 * (double)upper_entries ||
      4.0 * (double)upper_entries > 5.0 * (double)lower_entries) {
    return 0;
  }
  return 1;
}

__attribute__((noinline, cold))
static UF_long kls_metis_camd_group_size(UF_long n,
                                         const UF_long *col_ptr,
                                         const UF_long *row_idx) {
  {
    const char *group_env = getenv("KLS_METIS_CAMD_GROUP_SIZE");
    if (group_env != NULL && group_env[0] != '\0') {
      const unsigned long long requested = strtoull(group_env, NULL, 10);
      if (requested > 0u && requested <= (unsigned long long)n) {
        return (UF_long)requested;
      }
    }
  }
  {
    /* Generic NodeND uses only a resource-scaled refinement window.  Smaller
       graphs retain the exact NodeND order (group size zero is converted to
       singleton constraints by the caller); larger graphs amortize CAMD in
       progressively wider cache-bounded windows. */
    if (kls_generic_metis_tuning_ctx == 4) {
      /* When the retained minimum-degree estimate exposes exceptionally
         high arithmetic intensity, separator concurrency is worth more than
         a very wide refinement window.  Scale the window only with the
         execution resource count; the caller selected this context from
         realized symbolic work per active row, not graph dimensions or an
         input-family signature. */
      UF_long group = 384u *
        (UF_long)(kls_metis_order_threads > 0 ? kls_metis_order_threads : 1);
      if (group < 1024u) {
        group = 1024u;
      }
      if (group > 8192u) {
        group = 8192u;
      }
      return group < n ? group : n;
    }
    if (n >= 3000000u) {
      return 16384u;
    }
    if (n >= 1000000u) {
      return 8192u;
    }
    return n >= 200000u ? 4096u : 0u;
  }
  if (kls_asymmetric_bounded_degree_direct_metis_analyze_class ==
      KLS_ASYMMETRIC_BOUNDED_DEGREE_DIRECT_METIS_THIN_FRINGE) {
    /* Thin-fringe dominant-SCC systems pair their fourteen-leaf forest with
       moderately narrow constrained-AMD windows.  Wider defaults retain
       more fill in the row factor; narrower windows add ordering overhead. */
    return 3072u;
  }
  if (kls_sparse_full_diagonal_metis_row_ctx == 2) {
    /* The higher-fill resource class needs wider windows to preserve the
       useful NodeNDP forest.  Keep roughly 37 windows per worker under the
       eight-thread policy, round to a cache-friendly 128-column unit, and
       cap ordering resources. */
    UF_long group = n / 296u;
    group = ((group + 64u) / 128u) * 128u;
    if (group < 4096u) {
      group = 4096u;
    }
    if (group > 8192u) {
      group = 8192u;
    }
    return group;
  }
  if (kls_sparse_full_diagonal_metis_row_ctx == 1) {
    /* Scale the constrained-AMD window with the graph instead of keying it
       to one matrix dimension.  The clamps are ordering-resource limits;
       n/608 retains the measured fine-refinement density across the input
       band (864 columns at n=525825). */
    UF_long group = n / 608u;
    if (group < 320u) {
      group = 320u;
    }
    if (group > 4096u) {
      group = 4096u;
    }
    return group;
  }
  if (is_large_sparse_nearly_diagonal_spiked_metis_pattern(
        n, col_ptr, row_idx)) {
    /* Nearly diagonal spike: 4096-column windows retain a good row schedule,
       but 2304 reduces both fill and the H100 cycle while keeping CAMD's
       schedule-shaping benefit.  Raw NodeNDP has lower fill yet refactors
       more slowly, so do not skip refinement for this class. */
    return 2304u;
  }
  if (is_large_very_low_degree_full_diagonal_pattern(n, col_ptr, row_idx) &&
      n >= 125000u && n <= 175000u) {
    /* On the large single-SCC circuit shape, plain NodeND leaves enough
       fill to dominate 99 refactors.  A modest constrained-AMD window
       removes that fill without the analysis cost of refining thousands of
       exact separator components (G2_circuit: 3.13e9 -> 3.10e9 flops while
       retaining the cheaper ~0.47 s analysis). */
    return 320u;
  }
  if (kls_fragmented_metis_tuning_ctx == 2 && n >= 65536u) {
    /* The callback sees the dominant BTF core rather than the outer matrix.
       When diagnostic CAMD refinement is enabled, 1536-column windows retain
       more of NodeNDP's balanced forest than the generic 1024-column pass.
       The scoped context prevents explicit METIS from inheriting AUTO. */
    return 1536u;
  }
  if (kls_fragmented_metis_tuning_ctx == 3 && n >= 32768u) {
    /* Preserve the accepted eight-leaf NodeNDP forest while CAMD trims within
       moderately wide rank windows.  The scoped context is set only by the
       structurally proposed and post-symbolically guarded AUTO route. */
    return 2304u;
  }
  if (kls_sparse_full_diagonal_hubbed_metis_core_refine_profile(
        n, col_ptr, row_idx)) {
    /* Keep roughly 96 vertices per refined NodeND rank window, rounded to a
       cache-friendly unit.  The clamps bound ordering work across the full
       structural resource tier. */
    UF_long group = n / 96u;
    group = ((group + 64u) / 128u) * 128u;
    if (group < 512u) {
      group = 512u;
    }
    if (group > 2304u) {
      group = 2304u;
    }
    if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
      fprintf(stderr,
              "KLS metis: full-diagonal hubbed core refine n=%ld nnz=%ld "
              "group=%ld\n",
              (long)n, (long)col_ptr[n], (long)group);
    }
    return group;
  }
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

/* These symbolic-selection contexts are used by the common ordering path,
   irrespective of whether the optional mt-metis NodeND bridge is present.
   Keep them outside KLS_HAVE_MTMETIS: standard bundled-METIS builds do not
   define that feature macro but still use AMF3 and value-matched analysis. */
static _Thread_local int kls_prestatic_ordering_ctx;
static _Thread_local int kls_matched_low_pair_ordering_ctx;
static _Thread_local double kls_matched_low_pair_work_ctx;
static _Thread_local int kls_compact_missing_match_ordering_ctx;
static _Thread_local int kls_value_matched_ordering_ctx;
static _Thread_local int kls_amf3_power_override_active;
static _Thread_local double kls_amf3_power_override;
static _Thread_local int kls_amf3_exact_clique_bound;
static _Thread_local int kls_amf3_span_variant_selected;

#ifdef KLS_HAVE_MTMETIS
#include <mtmetis.h>

/* mt-metis NodeND bridge: converts the 64-bit idx graph to mt-metis's
   default 32-bit types, runs the threaded ND, converts back.
   mt-metis's dlthread communicators are process-global: concurrent
   calls from two orientation-analyze threads crash in its reduction
   buffers (pre2), so calls serialize on a mutex - each still runs
   its own internal thread team. */
static pthread_mutex_t kls_mtmetis_lock = PTHREAD_MUTEX_INITIALIZER;

/* callers that measured a win from the deterministic parallel NodeNDP
   set this around their analyze (currently the pre-static METIS
   refinement); everything else keeps serial NodeNDP's forest, which
   the separator-pipeline refactor engines are structure-sensitive to
   (pre2: 5x refactor loss under the parallel splitter's forest) */
static _Thread_local int kls_det_ndp_requested;

/* pattern-classified det-ndp routing, refreshed per analyze call.
   v12 measured the refine-only det-ndp scoping reverting v11's broad
   one-shot wins (Freescale1 0.886->1.477, circuit5M_dc 0.962->1.650,
   memchip 1.535->1.908): det-ndp routes for every analyze EXCEPT the
   legacy-RNG large-spral class (pre2), whose separator-pipeline
   engine needs the serial-NodeNDP legacy forests (cycle 0.787 with
   the containment vs 1.76 without). */
static _Thread_local int kls_lowdeg_ndp_class;
static _Thread_local int kls_detndp_class_ok;
/* Nearly diagonal, high-fill spiked systems also benefit from the
   separator-preserving splitter even though they share the legacy-RNG size
   envelope.  The sparse half retains nxp1's row-refactor forest; the dense
   half (rajat30) keeps its original 4096-wide CAMD refinement and compact
   matching route while removing roughly 1.6 seconds of NodeNDP work. */
static _Thread_local int kls_spiked_ndp_class;
/* Giant symmetric scalar-fringe systems use a seven-leaf deterministic
   forest: eight leaves shorten cold ordering slightly but leave a slower
   critical row-refactor chain, while the serial fourteen-leaf forest adds
   substantial cold NodeNDP work.  Scoped to the admitted AUTO race thread. */
static _Thread_local int kls_giant_symmetric_scalar_fringe_metis_ctx;
/* set around the pre-static ordering choice: the legacy libc RNG is
   scoped to exactly the analyses whose forests the separator-pipeline
   refactor engines consume (the prestatic-adopted ones).  The ana
   phase's two-orientation analyzes then use thread-local MT19937
   (deterministic under concurrency); the process-global libc rand was
   cross-corrupted by those concurrent analyzes - pre2's true
   run-to-run variance source. */
/* Set only while ordering a value-matched graph whose complete quotient
   work fits a small fixed budget.  The chooser can then settle BTF versus
   one-block AMD directly instead of launching a broader ordering portfolio. */
/* Any symbolic tournament over value-matched coordinates already belongs to
   a larger match/order/factor comparison.  Do not nest the lower-span AMF3
   micro-portfolio inside that value-dependent coordinate trial. */

/* A recurring-workload AMF3 portfolio may compare the ordinary clique
   amortization with a lower-power variant on independent worker threads.
   Keep that implementation parameter thread-local: orientation candidates
   and speculative ordering jobs can run concurrently, while an explicit
   KLS_AMF3_POWER remains authoritative for diagnostic callers. */
/* Published by the winning AUTO tournament on the candidate's own analyze
   thread and immediately copied into kls_pattern_candidate. */

/* per-thread routing request from the analyze level: the mostly-
   missing-diagonal low-degree class (mac_econ) pays a pathological
   serial NodeND while its egraph refactor engine has no separator-
   structure dependence - the measured safe class for threaded ND.
   pre2-like matrices are prestatic-routed before ordering and their
   pipeline engines NEED serial NodeNDP structure (measured 14.9 vs
   2.2s refactor) - never routed here. */
static _Thread_local int kls_mt_nd_class;

/* forest-aware threaded ND: replicates METIS_NodeNDP's output contract
   (perm/iperm + the reversed-index sizes tree consumed by
   kls_metis_append_separator_component: right subtree first, then
   left, then the separator).  Top separators via mt-metis (mutex-
   serialized; internally threaded), leaves via concurrent serial
   METIS_NodeND. */
struct kls_mtnd_leaf {
  idx_t nv;
  idx_t *xadj;
  idx_t *adjncy;
  idx_t *perm;
  idx_t *iperm;
  idx_t *options;
  int status;
  UF_long *vtx;               /* this leaf's global vertex ids */
  int vtx_owned;              /* freed here (task path) vs arena */
};

static void *kls_mtnd_leaf_main(void *arg) {
  struct kls_mtnd_leaf *L = (struct kls_mtnd_leaf *)arg;
  idx_t nv = L->nv;
  if (nv > 0) {
    L->status = METIS_NodeND(&nv, L->xadj, L->adjncy, NULL, L->options,
                             L->perm, L->iperm);
  } else {
    L->status = METIS_OK;
  }
  return NULL;
}

/* recursive splitter: fills emit[] (vertex ids in NodeNDP emission
   order), sizes[] (reversed index), and leaf job structs for deferred
   concurrent leaf ordering.  vertices/nv = this component's vertex
   ids (global).  Returns 0 on failure. */
struct kls_mtnd_seg {
  int leaf;                 /* leaf index, or -1 for a separator */
  UF_long *vertices;        /* owned for separators; leaf lists live in
                               leaf_vtx */
  UF_long count;
};

struct kls_mtnd_ctx {
  const idx_t *xadj;
  const idx_t *adjncy;
  idx_t npes;
  idx_t *sizes;
  struct kls_mtnd_seg *segs;
  int nsegs;
  struct kls_mtnd_leaf *leaves;
  UF_long *leaf_vtx;        /* concatenated leaf vertex lists */
  UF_long *leaf_vtx_start;
  UF_long leaf_vtx_used;
  UF_long n;
  UF_long *scratch_local;   /* size n: global -> local id map */
  int mt_threads;
  int serial_sep;           /* internal separators via serial METIS:
                               deterministic (mt-metis's parallel
                               coarsening is not seed-pinnable), with
                               the parallelism kept in the concurrent
                               serial leaf orderings */
  idx_t *lopts;             /* shared leaf NodeND options (read-only) */
};

static int kls_mtnd_split(struct kls_mtnd_ctx *C, idx_t cpos,
                          UF_long *vertices, UF_long nv) {
  const idx_t component_count = 2 * C->npes - 1;
  const idx_t ridx = component_count - 1 - cpos;
  if (cpos >= C->npes - 1) {
    /* leaf: record for deferred concurrent serial ND */
    int li = (int)(cpos - (C->npes - 1));
    C->leaves[li].nv = (idx_t)nv;
    C->leaf_vtx_start[li] = C->leaf_vtx_used;
    C->leaf_vtx_used += nv;
    C->segs[C->nsegs].leaf = li;
    C->segs[C->nsegs].vertices = NULL;
    C->segs[C->nsegs].count = nv;
    C->nsegs++;
    /* build the leaf subgraph (local ids) */
    {
      UF_long v, q, w = 0;
      idx_t *lx = (idx_t *)malloc(((size_t)nv + 1u) * sizeof(idx_t));
      idx_t *la = NULL;
      size_t edges = 0;
      if (lx == NULL) return 0;
      for (v = 0; v < nv; ++v) {
        C->scratch_local[vertices[v]] = v;
      }
      for (v = 0; v < nv; ++v) {
        edges += (size_t)(C->xadj[vertices[v] + 1] - C->xadj[vertices[v]]);
      }
      la = (idx_t *)malloc((edges + 1u) * sizeof(idx_t));
      if (la == NULL) { free(lx); return 0; }
      for (v = 0; v < nv; ++v) {
        UF_long gv = vertices[v];
        lx[v] = (idx_t)w;
        for (q = C->xadj[gv]; q < C->xadj[gv + 1]; ++q) {
          UF_long gu = (UF_long)C->adjncy[q];
          /* neighbor is in this leaf iff it appears in vertices[]:
             mark via scratch_local + membership stamp */
          if (C->scratch_local[gu] < nv &&
              vertices[C->scratch_local[gu]] == gu) {
            la[w++] = (idx_t)C->scratch_local[gu];
          }
        }
      }
      lx[nv] = (idx_t)w;
      C->leaves[li].xadj = lx;
      C->leaves[li].adjncy = la;
      C->leaves[li].perm = (idx_t *)malloc((size_t)(nv > 0 ? nv : 1) *
                                           sizeof(idx_t));
      C->leaves[li].iperm = (idx_t *)malloc((size_t)(nv > 0 ? nv : 1) *
                                            sizeof(idx_t));
      if (C->leaves[li].perm == NULL || C->leaves[li].iperm == NULL) {
        return 0;
      }
      /* remember the vertex list location for the emit phase */
      memcpy(C->leaf_vtx + C->leaf_vtx_start[li], vertices,
             (size_t)nv * sizeof(UF_long));
      C->leaves[li].vtx = C->leaf_vtx + C->leaf_vtx_start[li];
      C->leaves[li].vtx_owned = 0;
    }
    C->sizes[ridx] = (idx_t)nv;
    return 1;
  }
  /* internal: vertex separator on this component */
  if (C->serial_sep) {
    idx_t snv = (idx_t)nv;
    idx_t *lx = (idx_t *)malloc(((size_t)nv + 1u) * sizeof(idx_t));
    idx_t *la = NULL;
    idx_t *lwhere =
      (idx_t *)malloc((size_t)(nv > 0 ? nv : 1) * sizeof(idx_t));
    idx_t sepsize = 0;
    size_t edges = 0;
    UF_long v, q, w = 0;
    int st;
    UF_long *partA, *partB, *partS;
    UF_long na = 0, nb = 0, nsep = 0;
    if (lx == NULL || lwhere == NULL) {
      free(lx); free(lwhere);
      return 0;
    }
    for (v = 0; v < nv; ++v) {
      C->scratch_local[vertices[v]] = v;
    }
    for (v = 0; v < nv; ++v) {
      edges += (size_t)(C->xadj[vertices[v] + 1] - C->xadj[vertices[v]]);
    }
    la = (idx_t *)malloc((edges + 1u) * sizeof(idx_t));
    if (la == NULL) { free(lx); free(lwhere); return 0; }
    for (v = 0; v < nv; ++v) {
      UF_long gv = vertices[v];
      lx[v] = (idx_t)w;
      for (q = C->xadj[gv]; q < C->xadj[gv + 1]; ++q) {
        UF_long gu = (UF_long)C->adjncy[q];
        if (C->scratch_local[gu] < nv &&
            vertices[C->scratch_local[gu]] == gu) {
          la[w++] = (idx_t)C->scratch_local[gu];
        }
      }
    }
    lx[nv] = (idx_t)w;
    /* serial METIS separator: deterministic for fixed inputs */
    st = METIS_ComputeVertexSeparator(&snv, lx, la, NULL, NULL,
                                      &sepsize, lwhere);
    free(lx); free(la);
    if (st != METIS_OK) { free(lwhere); return 0; }
    partA = (UF_long *)malloc((size_t)(nv > 0 ? nv : 1) * sizeof(UF_long));
    partB = (UF_long *)malloc((size_t)(nv > 0 ? nv : 1) * sizeof(UF_long));
    partS = (UF_long *)malloc((size_t)(nv > 0 ? nv : 1) * sizeof(UF_long));
    if (partA == NULL || partB == NULL || partS == NULL) {
      free(lwhere); free(partA); free(partB); free(partS);
      return 0;
    }
    for (v = 0; v < nv; ++v) {
      if (lwhere[v] == 2) partS[nsep++] = vertices[v];
      else if (lwhere[v] == 0) partA[na++] = vertices[v];
      else partB[nb++] = vertices[v];
    }
    free(lwhere);
    C->sizes[ridx] = (idx_t)nsep;
    /* right child (2*cpos+2) = partB first, then left = partA */
    if (!kls_mtnd_split(C, 2 * cpos + 2, partB, nb) ||
        !kls_mtnd_split(C, 2 * cpos + 1, partA, na)) {
      free(partA); free(partB); free(partS);
      return 0;
    }
    /* the separator segment follows both subtrees */
    C->segs[C->nsegs].leaf = -1;
    C->segs[C->nsegs].vertices = partS;
    C->segs[C->nsegs].count = nsep;
    C->nsegs++;
    free(partA); free(partB);
    return 1;
  }
  {
    mtmetis_vtx_type mn = (mtmetis_vtx_type)nv;
    mtmetis_adj_type *mx =
      (mtmetis_adj_type *)malloc(((size_t)nv + 1u) * sizeof(*mx));
    mtmetis_vtx_type *ma = NULL;
    mtmetis_pid_type *mwhere =
      (mtmetis_pid_type *)malloc((size_t)nv * sizeof(*mwhere));
    mtmetis_wgt_type msep = 0;
    double *opts = mtmetis_init_options();
    size_t edges = 0;
    UF_long v, q, w = 0;
    int st;
    UF_long *partA, *partB, *partS;
    UF_long na = 0, nb = 0, nsep = 0;
    if (mx == NULL || mwhere == NULL || opts == NULL) {
      free(mx); free(mwhere); free(opts);
      return 0;
    }
    for (v = 0; v < nv; ++v) {
      C->scratch_local[vertices[v]] = v;
    }
    for (v = 0; v < nv; ++v) {
      edges += (size_t)(C->xadj[vertices[v] + 1] - C->xadj[vertices[v]]);
    }
    ma = (mtmetis_vtx_type *)malloc((edges + 1u) * sizeof(*ma));
    if (ma == NULL) { free(mx); free(mwhere); free(opts); return 0; }
    for (v = 0; v < nv; ++v) {
      UF_long gv = vertices[v];
      mx[v] = (mtmetis_adj_type)w;
      for (q = C->xadj[gv]; q < C->xadj[gv + 1]; ++q) {
        UF_long gu = (UF_long)C->adjncy[q];
        if (C->scratch_local[gu] < nv &&
            vertices[C->scratch_local[gu]] == gu) {
          ma[w++] = (mtmetis_vtx_type)C->scratch_local[gu];
        }
      }
    }
    mx[nv] = (mtmetis_adj_type)w;
    opts[MTMETIS_OPTION_NTHREADS] = (double)C->mt_threads;
    opts[MTMETIS_OPTION_SEED] = 0.0;
    pthread_mutex_lock(&kls_mtmetis_lock);
    st = MTMETIS_ComputeVertexSeparator(&mn, mx, ma, NULL, opts, &msep,
                                        mwhere);
    pthread_mutex_unlock(&kls_mtmetis_lock);
    free(mx); free(ma); free(opts);
    if (st != MTMETIS_SUCCESS) { free(mwhere); return 0; }
    partA = (UF_long *)malloc((size_t)(nv > 0 ? nv : 1) * sizeof(UF_long));
    partB = (UF_long *)malloc((size_t)(nv > 0 ? nv : 1) * sizeof(UF_long));
    partS = (UF_long *)malloc((size_t)(nv > 0 ? nv : 1) * sizeof(UF_long));
    if (partA == NULL || partB == NULL || partS == NULL) {
      free(mwhere); free(partA); free(partB); free(partS);
      return 0;
    }
    for (v = 0; v < nv; ++v) {
      if (mwhere[v] == MTMETIS_VSEP_SEP) partS[nsep++] = vertices[v];
      else if (mwhere[v] == MTMETIS_VSEP_PARTA) partA[na++] = vertices[v];
      else partB[nb++] = vertices[v];
    }
    free(mwhere);
    C->sizes[ridx] = (idx_t)nsep;
    /* right child (2*cpos+2) = partB first, then left = partA */
    if (!kls_mtnd_split(C, 2 * cpos + 2, partB, nb) ||
        !kls_mtnd_split(C, 2 * cpos + 1, partA, na)) {
      free(partA); free(partB); free(partS);
      return 0;
    }
    /* the separator segment follows both subtrees */
    C->segs[C->nsegs].leaf = -1;
    C->segs[C->nsegs].vertices = partS;
    C->segs[C->nsegs].count = nsep;
    C->nsegs++;
    free(partA); free(partB);
    return 1;
  }
}

/* Task-parallel variant of the serial_sep splitter: the level-L
   separators are mutually independent given their parents, so each
   internal node (and each leaf ordering) runs on a FRESH thread.
   Determinism: a fresh thread's gklib MT19937 stream starts from the
   thread-local default state and only that node's single METIS call
   consumes it, so every node's result depends on its subgraph alone -
   no call-order or interleaving dependence (the spawn-failure inline
   fallback shares the parent's stream and can differ, but only under
   thread exhaustion).  Segment order is preserved by construction:
   each task returns its subtree's ordered segment list and the parent
   concatenates right-subtree, left-subtree, own separator - the same
   postorder the serial recursion emits.  sizes[] writes are disjoint
   per node.  Each task carries its own global->local scratch (calloc:
   a zero entry either names the true local id 0 or fails the
   vertices[] membership double-check). */
struct kls_mtnd_task {
  struct kls_mtnd_ctx *C;
  UF_long *vertices;          /* owned: leaf keeps, internal frees */
  UF_long nv;
  idx_t cpos;
  struct kls_mtnd_seg *segs;  /* owned ordered subtree segments */
  int nsegs;
  int ok;
};

static void kls_mtnd_task_free_segs(struct kls_mtnd_task *T) {
  int i;
  for (i = 0; i < T->nsegs; ++i) {
    free(T->segs[i].vertices);
  }
  T->nsegs = 0;
}

static void *kls_mtnd_task_main(void *arg) {
  struct kls_mtnd_task *T = (struct kls_mtnd_task *)arg;
  struct kls_mtnd_ctx *C = T->C;
  const idx_t component_count = 2 * C->npes - 1;
  const idx_t ridx = component_count - 1 - T->cpos;
  const UF_long nv = T->nv;
  UF_long *vertices = T->vertices;
  UF_long v, q, w = 0;
  size_t edges = 0;
  T->ok = 0;
  T->nsegs = 0;
  T->segs = (struct kls_mtnd_seg *)calloc((size_t)(2 * C->npes),
                                          sizeof(*T->segs));
  UF_long *scratch = (UF_long *)calloc((size_t)C->n, sizeof(*scratch));
  idx_t *lx = (idx_t *)malloc(((size_t)nv + 1u) * sizeof(idx_t));
  idx_t *la = NULL;
  if (T->segs == NULL || scratch == NULL || lx == NULL) {
    free(scratch); free(lx);
    return NULL;
  }
  for (v = 0; v < nv; ++v) {
    scratch[vertices[v]] = v;
  }
  for (v = 0; v < nv; ++v) {
    edges += (size_t)(C->xadj[vertices[v] + 1] - C->xadj[vertices[v]]);
  }
  la = (idx_t *)malloc((edges + 1u) * sizeof(idx_t));
  if (la == NULL) {
    free(scratch); free(lx);
    return NULL;
  }
  for (v = 0; v < nv; ++v) {
    UF_long gv = vertices[v];
    lx[v] = (idx_t)w;
    for (q = C->xadj[gv]; q < C->xadj[gv + 1]; ++q) {
      UF_long gu = (UF_long)C->adjncy[q];
      if (scratch[gu] < nv && vertices[scratch[gu]] == gu) {
        la[w++] = (idx_t)scratch[gu];
      }
    }
  }
  lx[nv] = (idx_t)w;
  free(scratch);
  if (T->cpos >= C->npes - 1) {
    /* leaf: order inline on this (fresh) thread */
    int li = (int)(T->cpos - (C->npes - 1));
    struct kls_mtnd_leaf *L = &C->leaves[li];
    L->nv = (idx_t)nv;
    L->xadj = lx;
    L->adjncy = la;
    L->perm = (idx_t *)malloc((size_t)(nv > 0 ? nv : 1) * sizeof(idx_t));
    L->iperm = (idx_t *)malloc((size_t)(nv > 0 ? nv : 1) * sizeof(idx_t));
    L->options = C->lopts;
    L->vtx = vertices;
    L->vtx_owned = 1;
    T->vertices = NULL;       /* ownership moved to the leaf */
    if (L->perm == NULL || L->iperm == NULL) {
      return NULL;
    }
    kls_mtnd_leaf_main(L);
    if (L->status != METIS_OK) {
      return NULL;
    }
    T->segs[0].leaf = li;
    T->segs[0].vertices = NULL;
    T->segs[0].count = nv;
    T->nsegs = 1;
    C->sizes[ridx] = (idx_t)nv;
    T->ok = 1;
    return NULL;
  }
  {
    idx_t snv = (idx_t)nv;
    idx_t *lwhere =
      (idx_t *)malloc((size_t)(nv > 0 ? nv : 1) * sizeof(idx_t));
    idx_t sepsize = 0;
    UF_long *partA, *partB, *partS;
    UF_long na = 0, nb = 0, nsep = 0;
    int st, i;
    if (lwhere == NULL) {
      free(lx); free(la);
      return NULL;
    }
    st = METIS_ComputeVertexSeparator(&snv, lx, la, NULL, NULL,
                                      &sepsize, lwhere);
    free(lx); free(la);
    if (st != METIS_OK) {
      free(lwhere);
      return NULL;
    }
    partA = (UF_long *)malloc((size_t)(nv > 0 ? nv : 1) * sizeof(UF_long));
    partB = (UF_long *)malloc((size_t)(nv > 0 ? nv : 1) * sizeof(UF_long));
    partS = (UF_long *)malloc((size_t)(nv > 0 ? nv : 1) * sizeof(UF_long));
    if (partA == NULL || partB == NULL || partS == NULL) {
      free(lwhere); free(partA); free(partB); free(partS);
      return NULL;
    }
    for (v = 0; v < nv; ++v) {
      if (lwhere[v] == 2) partS[nsep++] = vertices[v];
      else if (lwhere[v] == 0) partA[na++] = vertices[v];
      else partB[nb++] = vertices[v];
    }
    free(lwhere);
    free(vertices);
    T->vertices = NULL;
    C->sizes[ridx] = (idx_t)nsep;
    struct kls_mtnd_task right, left;
    pthread_t rtid, ltid;
    memset(&right, 0, sizeof(right));
    memset(&left, 0, sizeof(left));
    right.C = C; right.vertices = partB; right.nv = nb;
    right.cpos = 2 * T->cpos + 2;
    left.C = C; left.vertices = partA; left.nv = na;
    left.cpos = 2 * T->cpos + 1;
    const int r_sp = pthread_create(&rtid, NULL, kls_mtnd_task_main,
                                    &right) == 0;
    const int l_sp = pthread_create(&ltid, NULL, kls_mtnd_task_main,
                                    &left) == 0;
    if (!r_sp) kls_mtnd_task_main(&right);
    if (!l_sp) kls_mtnd_task_main(&left);
    if (r_sp) pthread_join(rtid, NULL);
    if (l_sp) pthread_join(ltid, NULL);
    if (!right.ok || !left.ok) {
      kls_mtnd_task_free_segs(&right);
      kls_mtnd_task_free_segs(&left);
      free(right.segs); free(left.segs);
      free(right.vertices); free(left.vertices);
      free(partS);
      return NULL;
    }
    for (i = 0; i < right.nsegs; ++i) {
      T->segs[T->nsegs++] = right.segs[i];
    }
    for (i = 0; i < left.nsegs; ++i) {
      T->segs[T->nsegs++] = left.segs[i];
    }
    free(right.segs); free(left.segs);
    T->segs[T->nsegs].leaf = -1;
    T->segs[T->nsegs].vertices = partS;
    T->segs[T->nsegs].count = nsep;
    T->nsegs++;
    T->ok = 1;
    return NULL;
  }
}

/* forest-aware entry: NodeNDP replacement producing perm/iperm and
   NodeNDP-format sizes.  serial_sep=1 keeps every separator serial
   METIS (fully deterministic; parallelism from concurrent leaves
   only).  Returns 1 on success. */
static int kls_mtmetis_ndp(UF_long n, const idx_t *xadj, const idx_t *adjncy,
                           idx_t npes, idx_t *metis_perm, idx_t *metis_iperm,
                           idx_t *sizes, int serial_sep) {
  struct kls_mtnd_ctx C;
  UF_long *all = (UF_long *)malloc((size_t)n * sizeof(UF_long));
  int nleaves = (int)npes;
  int ok = 0, i;
  UF_long v;
  memset(&C, 0, sizeof(C));
  C.xadj = xadj;
  C.adjncy = adjncy;
  C.npes = npes;
  C.sizes = sizes;
  C.n = n;
  C.serial_sep = serial_sep;
  C.mt_threads = 4;
  {
    const char *te = getenv("KLS_MT_ND");
    if (te != NULL && te[0] != '\0') {
      int p = atoi(te);
      if (p >= 1 && p <= 32) C.mt_threads = p;
    }
  }
  C.segs = (struct kls_mtnd_seg *)calloc((size_t)(2 * npes), sizeof(*C.segs));
  C.leaves = (struct kls_mtnd_leaf *)calloc((size_t)nleaves,
                                            sizeof(*C.leaves));
  C.leaf_vtx = (UF_long *)malloc((size_t)n * sizeof(UF_long));
  C.leaf_vtx_start = (UF_long *)calloc((size_t)nleaves,
                                       sizeof(*C.leaf_vtx_start));
  C.scratch_local = (UF_long *)malloc((size_t)n * sizeof(UF_long));
  if (all == NULL || C.segs == NULL || C.leaves == NULL ||
      C.leaf_vtx == NULL || C.leaf_vtx_start == NULL ||
      C.scratch_local == NULL) {
    goto cleanup;
  }
  for (v = 0; v < n; ++v) {
    all[v] = v;
    C.scratch_local[v] = (UF_long)n;
  }
  idx_t lopts[METIS_NOPTIONS];
  METIS_SetDefaultOptions(lopts);
  lopts[METIS_OPTION_NUMBERING] = 0;
  lopts[METIS_OPTION_SEED] = 0;
  C.lopts = lopts;
  const double kls_mtnd_t0 = kls_now_seconds();
  if (serial_sep && getenv("KLS_MTND_SERIAL_SPLIT") == NULL) {
    /* task recursion: level-parallel separators + in-task leaves */
    struct kls_mtnd_task root;
    pthread_t rtid;
    memset(&root, 0, sizeof(root));
    root.C = &C;
    root.vertices = all;
    root.nv = n;
    root.cpos = 0;
    all = NULL;               /* ownership moved to the root task */
    if (pthread_create(&rtid, NULL, kls_mtnd_task_main, &root) == 0) {
      pthread_join(rtid, NULL);
    } else {
      kls_mtnd_task_main(&root);
    }
    ok = root.ok && root.nsegs <= 2 * (int)npes;
    if (ok) {
      memcpy(C.segs, root.segs,
             (size_t)root.nsegs * sizeof(*root.segs));
      C.nsegs = root.nsegs;
    } else {
      kls_mtnd_task_free_segs(&root);
    }
    free(root.segs);
    free(root.vertices);
    if (!ok) {
      goto cleanup;
    }
  } else {
  if (!kls_mtnd_split(&C, 0, all, n)) {
    goto cleanup;
  }
  const double kls_mtnd_t1 = kls_now_seconds();
  /* leaf orderings, concurrently (serial METIS is safe across calls) */
  {
    pthread_t tids[64];
    int spawned[64];
    for (i = 0; i < nleaves; ++i) {
      C.leaves[i].options = lopts;
      spawned[i] = i > 0 &&
        pthread_create(&tids[i], NULL, kls_mtnd_leaf_main,
                       &C.leaves[i]) == 0;
    }
    kls_mtnd_leaf_main(&C.leaves[0]);
    ok = C.leaves[0].status == METIS_OK;
    for (i = 1; i < nleaves; ++i) {
      if (spawned[i]) {
        pthread_join(tids[i], NULL);
      } else {
        kls_mtnd_leaf_main(&C.leaves[i]);
      }
      ok = ok && C.leaves[i].status == METIS_OK;
    }
  }
  if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
    fprintf(stderr, "KLS mtnd: split %.3fs leaves %.3fs (ss=%d npes=%d)\n",
            kls_mtnd_t1 - kls_mtnd_t0, kls_now_seconds() - kls_mtnd_t1,
            serial_sep, (int)npes);
  }
  }
  if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
    fprintf(stderr, "KLS mtnd: total %.3fs (ss=%d npes=%d)\n",
            kls_now_seconds() - kls_mtnd_t0, serial_sep, (int)npes);
  }
  if (ok) {
    /* compose the global order from the traversal segments */
    UF_long pos = 0;
    int si;
    for (si = 0; si < C.nsegs; ++si) {
      struct kls_mtnd_seg *g = &C.segs[si];
      if (g->leaf >= 0) {
        struct kls_mtnd_leaf *L = &C.leaves[g->leaf];
        UF_long k;
        for (k = 0; k < g->count; ++k) {
          UF_long gv = L->vtx[(UF_long)L->perm[k]];
          metis_iperm[gv] = (idx_t)pos;
          metis_perm[pos] = (idx_t)gv;
          pos++;
        }
      } else {
        UF_long k;
        for (k = 0; k < g->count; ++k) {
          UF_long gv = g->vertices[k];
          metis_iperm[gv] = (idx_t)pos;
          metis_perm[pos] = (idx_t)gv;
          pos++;
        }
      }
    }
    ok = pos == n;
  }
cleanup:
  if (C.segs != NULL) {
    for (i = 0; i < C.nsegs; ++i) {
      free(C.segs[i].vertices);
    }
  }
  if (C.leaves != NULL) {
    for (i = 0; i < nleaves; ++i) {
      free(C.leaves[i].xadj);
      free(C.leaves[i].adjncy);
      free(C.leaves[i].perm);
      free(C.leaves[i].iperm);
      if (C.leaves[i].vtx_owned) {
        free(C.leaves[i].vtx);
      }
    }
  }
  free(C.segs);
  free(C.leaves);
  free(C.leaf_vtx);
  free(C.leaf_vtx_start);
  free(C.scratch_local);
  free(all);
  return ok;
}

static int kls_mtmetis_nd(UF_long n, const idx_t *xadj, const idx_t *adjncy,
                          idx_t *metis_perm, idx_t *metis_iperm) {
  const size_t nz = (size_t)xadj[n];
  mtmetis_adj_type *mx =
    (mtmetis_adj_type *)malloc(((size_t)n + 1u) * sizeof(*mx));
  mtmetis_vtx_type *ma =
    (mtmetis_vtx_type *)malloc((nz + 1u) * sizeof(*ma));
  mtmetis_pid_type *mp =
    (mtmetis_pid_type *)malloc((size_t)n * sizeof(*mp));
  mtmetis_pid_type *mi =
    (mtmetis_pid_type *)malloc((size_t)n * sizeof(*mi));
  double *opts = mtmetis_init_options();
  int ok = 0;
  if (mx != NULL && ma != NULL && mp != NULL && mi != NULL &&
      opts != NULL) {
    UF_long i;
    int threads = 4;
    const char *te = getenv("KLS_MT_ND");
    if (te != NULL && te[0] != '\0') {
      int p = atoi(te);
      if (p >= 1 && p <= 32) {
        threads = p;
      }
    }
    for (i = 0; i <= n; ++i) {
      mx[i] = (mtmetis_adj_type)xadj[i];
    }
    for (i = 0; i < (UF_long)nz; ++i) {
      ma[i] = (mtmetis_vtx_type)adjncy[i];
    }
    opts[MTMETIS_OPTION_NTHREADS] = (double)threads;
    opts[MTMETIS_OPTION_SEED] = 0.0;
    if (getenv("KLS_MT_ND_METIS") != NULL) {
      /* hybrid: serial-METIS quality on the subproblems, threaded top */
      opts[MTMETIS_OPTION_METIS] = 1.0;
    }
    {
      const char *nr = getenv("KLS_MT_ND_NRUNS");
      if (nr != NULL && nr[0] != '\0') {
        opts[MTMETIS_OPTION_NRUNS] = (double)atoi(nr);
      }
    }
    if (getenv("KLS_MT_ND_TRACE") != NULL) {
      UF_long bad = 0, j;
      for (j = 0; j < n; ++j) {
        if (xadj[j] > xadj[j + 1]) bad++;
      }
      for (j = 0; j < (UF_long)nz; ++j) {
        if (adjncy[j] < 0 || adjncy[j] >= (idx_t)n) bad++;
      }
      {
        UF_long iso = 0;
        for (j = 0; j < n; ++j) {
          if (xadj[j] == xadj[j + 1]) iso++;
        }
        fprintf(stderr,
                "KLS mtnd: n=%ld nz=%zu bad=%ld iso=%ld threads=%d\n",
                (long)n, nz, (long)bad, (long)iso, threads);
        {
          const char *dump = getenv("KLS_MT_ND_DUMP");
          if (dump != NULL) {
            FILE *f = fopen(dump, "w");
            if (f != NULL) {
              UF_long v, q;
              fprintf(f, "%ld %zu\n", (long)n, nz / 2u);
              for (v = 0; v < n; ++v) {
                for (q = xadj[v]; q < xadj[v + 1]; ++q) {
                  fprintf(f, "%ld ", (long)adjncy[q] + 1);
                }
                fputc('\n', f);
              }
              fclose(f);
            }
          }
        }
      }
    }
    {
      mtmetis_vtx_type mn = (mtmetis_vtx_type)n;
      pthread_mutex_lock(&kls_mtmetis_lock);
      ok = MTMETIS_NodeND(&mn, mx, ma, NULL, opts, mp, mi) ==
           MTMETIS_SUCCESS;
      pthread_mutex_unlock(&kls_mtmetis_lock);
    }
    if (ok) {
      for (i = 0; i < n; ++i) {
        metis_perm[i] = (idx_t)mp[i];
        metis_iperm[i] = (idx_t)mi[i];
      }
    }
  }
  free(mx);
  free(ma);
  free(mp);
  free(mi);
  free(opts);
  return ok;
}
#endif

struct kls_par_nd_half {
  idx_t nv;
  idx_t *xadj;
  idx_t *adjncy;
  idx_t *perm;   /* local iperm-style: new order of local vertex */
  idx_t *iperm;
  idx_t *options;
  int status;
};

static void *kls_par_nd_half_main(void *arg) {
  struct kls_par_nd_half *h = (struct kls_par_nd_half *)arg;
  idx_t nv = h->nv;
  h->status = nv > 0
    ? METIS_NodeND(&nv, h->xadj, h->adjncy, NULL, h->options, h->perm,
                   h->iperm)
    : METIS_OK;
  return NULL;
}

/* One threaded top level of nested dissection.  Returns 1 on success
   with metis_perm/metis_iperm filled (NodeND conventions: iperm[v] =
   position of v in the elimination order, perm = its inverse). */
static int kls_metis_parallel_nd(idx_t nvtxs, const idx_t *xadj,
                                 const idx_t *adjncy, idx_t *options,
                                 idx_t *metis_perm, idx_t *metis_iperm) {
  const size_t nz = (size_t)xadj[nvtxs];
  idx_t *part = (idx_t *)malloc((size_t)nvtxs * sizeof(idx_t));
  idx_t *local = (idx_t *)malloc((size_t)nvtxs * sizeof(idx_t));
  idx_t *gx = (idx_t *)malloc((size_t)nvtxs * sizeof(idx_t));
  if (part == NULL || local == NULL || gx == NULL) {
    free(part); free(local); free(gx);
    return 0;
  }
  idx_t sepsize = 0;
  idx_t nv = nvtxs;
  /* the separator call mutates nothing but part/sepsize */
  int st = METIS_ComputeVertexSeparator(&nv, (idx_t *)xadj,
                                        (idx_t *)adjncy, NULL, NULL,
                                        &sepsize, part);
  if (st != METIS_OK || sepsize <= 0 || sepsize >= nvtxs / 2) {
    free(part); free(local); free(gx);
    return 0;
  }
  idx_t cnt[3] = {0, 0, 0};
  for (idx_t v = 0; v < nvtxs; ++v) {
    idx_t p = part[v];
    if (p < 0 || p > 2) { p = 2; part[v] = 2; }
    local[v] = cnt[p]++;
  }
  if (cnt[0] == 0 || cnt[1] == 0) {
    free(part); free(local); free(gx);
    return 0;
  }
  struct kls_par_nd_half h[2];
  memset(h, 0, sizeof(h));
  int ok = 1;
  for (int side = 0; side < 2 && ok; ++side) {
    idx_t hn = cnt[side];
    h[side].nv = hn;
    h[side].options = options;
    h[side].xadj = (idx_t *)malloc(((size_t)hn + 1u) * sizeof(idx_t));
    h[side].adjncy = (idx_t *)malloc((nz + 1u) * sizeof(idx_t));
    h[side].perm = (idx_t *)malloc((size_t)hn * sizeof(idx_t));
    h[side].iperm = (idx_t *)malloc((size_t)hn * sizeof(idx_t));
    h[side].status = METIS_ERROR;
    if (h[side].xadj == NULL || h[side].adjncy == NULL ||
        h[side].perm == NULL || h[side].iperm == NULL) {
      ok = 0;
      break;
    }
    idx_t w = 0, lv = 0;
    for (idx_t v = 0; v < nvtxs; ++v) {
      if (part[v] != side) continue;
      h[side].xadj[lv] = w;
      for (idx_t p = xadj[v]; p < xadj[v + 1]; ++p) {
        idx_t u = adjncy[p];
        if (part[u] == side) {
          h[side].adjncy[w++] = local[u];
        }
      }
      lv++;
    }
    h[side].xadj[lv] = w;
  }
  if (ok) {
    pthread_t tid;
    int threaded = pthread_create(&tid, NULL, kls_par_nd_half_main,
                                  &h[1]) == 0;
    kls_par_nd_half_main(&h[0]);
    if (threaded) {
      pthread_join(tid, NULL);
    } else {
      kls_par_nd_half_main(&h[1]);
    }
    ok = h[0].status == METIS_OK && h[1].status == METIS_OK;
  }
  if (ok) {
    /* global elimination positions: half0 first, half1, separator
       last (separator kept in input order - CAMD refinement below
       polishes it) */
    idx_t base1 = cnt[0];
    idx_t base2 = cnt[0] + cnt[1];
    idx_t sep_seen = 0;
    for (idx_t v = 0; v < nvtxs; ++v) {
      idx_t pos;
      if (part[v] == 0) {
        pos = h[0].iperm[local[v]];
      } else if (part[v] == 1) {
        pos = base1 + h[1].iperm[local[v]];
      } else {
        pos = base2 + sep_seen++;
      }
      metis_iperm[v] = pos;
      gx[pos] = v;
    }
    for (idx_t k = 0; k < nvtxs; ++k) {
      metis_perm[k] = gx[k];
    }
  }
  for (int side = 0; side < 2; ++side) {
    free(h[side].xadj);
    free(h[side].adjncy);
    free(h[side].perm);
    free(h[side].iperm);
  }
  free(part);
  free(local);
  free(gx);
  return ok;
}

static UF_long kls_metis_order_inner(UF_long n,
                               UF_long *col_ptr,
                               UF_long *row_idx,
                               UF_long *perm_out,
                               trilinos_klu_l_common *common) {
  kls_metis_order_context *metis_context =
    common != NULL ? (kls_metis_order_context *)common->user_data : NULL;
  UF_long order_call_index = KLS_KLU_EMPTY;
  if (metis_context != NULL) {
    order_call_index = metis_context->order_call_count++;
  }
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
  idx_t metis_ndp_npes =
    metis_context != NULL && metis_context->npes > 1 ? metis_context->npes : 0;
#ifdef KLS_HAVE_MTMETIS
  if (kls_giant_symmetric_scalar_fringe_metis_ctx) {
    const char *npes_env = getenv(
      "KLS_GIANT_SYMMETRIC_SCALAR_FRINGE_METIS_NDP_NPES");
    metis_ndp_npes = npes_env != NULL && npes_env[0] != '\0'
      ? (idx_t)atol(npes_env) : 7;
  } else
#endif
  if (kls_asymmetric_bounded_degree_direct_metis_analyze_class) {
    /* A coarser SCC fringe favors eight leaves for row-refactor throughput.
       Thin-fringe systems retain fourteen: that forest's P/Q layout cuts
       their packed solve substantially. */
    const char *npes_env = getenv(
      "KLS_ASYMMETRIC_BOUNDED_DEGREE_DIRECT_METIS_NDP_NPES");
    metis_ndp_npes = npes_env != NULL && npes_env[0] != '\0'
      ? (idx_t)atol(npes_env)
      : (kls_asymmetric_bounded_degree_direct_metis_analyze_class ==
           KLS_ASYMMETRIC_BOUNDED_DEGREE_DIRECT_METIS_THIN_FRINGE
           ? 14 : 8);
  } else
  if (kls_fragmented_metis_tuning_ctx == 3) {
    const char *npes_env =
      getenv("KLS_SPARSE_SYMMETRIC_FRAGMENTED_METIS_NDP_NPES");
    metis_ndp_npes =
      npes_env != NULL && npes_env[0] != '\0' ? (idx_t)atol(npes_env) : 8;
  } else if (kls_medium_partial_static_metis_ctx ||
             kls_generic_metis_tuning_ctx) {
    /* A complete no-BTF graph exposes two deterministic leaves per numeric
       worker; BTF representations use one leaf per worker because their SCC
       structure supplies the outer decomposition.  Matched coordinates use
       their measured separator economics through the generic tuning context. */
    const char *generic_npes_env =
      getenv("KLS_GENERIC_METIS_NDP_NPES");
    if (kls_generic_metis_tuning_ctx && generic_npes_env != NULL &&
        generic_npes_env[0] != '\0') {
      metis_ndp_npes = (idx_t)atol(generic_npes_env);
    } else if (kls_generic_metis_tuning_ctx) {
      int requested_leaves =
        kls_generic_metis_tuning_ctx == 2 ||
        kls_generic_metis_tuning_ctx == 3 ||
        kls_generic_metis_tuning_ctx == 7 ||
        kls_generic_metis_tuning_ctx == 8
          ? kls_metis_order_threads
          : 2 * kls_metis_order_threads;
      if (requested_leaves < 2) {
        requested_leaves = 2;
      }
      if (requested_leaves > 64) {
        requested_leaves = 64;
      }
      metis_ndp_npes = (idx_t)requested_leaves;
    } else {
      metis_ndp_npes = 6;
    }
  } else if (kls_large_weak_diagonal_static_metis_ctx) {
    /* The two-leaf forest substantially reduces fill on this statically
       matched class.  On the audited 659K-row case it cuts modeled numeric
       work from 124.7B to 101.3B flops and the 8-thread steady refactor from
       about 1.53s to 1.05s.  A wider forest creates extra U fill even though
       it exposes more separator components. */
    metis_ndp_npes = 2;
  } else {
    const char *npes_env = getenv("KLS_METIS_NDP_NPES");
    if (npes_env != NULL && npes_env[0] != '\0') {
      metis_ndp_npes = (idx_t)atol(npes_env);
    }
  }
  if (kls_metis_ndp_npes_override > 0) {
    metis_ndp_npes = (idx_t)kls_metis_ndp_npes_override;
  }
  if (n < KLS_METIS_NDP_MIN_ROWS) {
    metis_ndp_npes = 0;
  }
  if (metis_ndp_npes > (idx_t)n) {
    metis_ndp_npes = (idx_t)n;
  }
  idx_t *metis_ndp_sizes = NULL;
  if (metis_ndp_npes > 1) {
    metis_ndp_sizes =
      (idx_t *)malloc((size_t)(2 * metis_ndp_npes - 1) * sizeof(*metis_ndp_sizes));
    if (metis_ndp_sizes == NULL) {
      metis_ndp_npes = 0;
    }
  }
  if (xadj == NULL || adjncy == NULL || metis_perm == NULL || metis_iperm == NULL) {
    free(degree);
    free(xadj);
    free(adjncy);
    free(metis_perm);
    free(metis_iperm);
    free(metis_ndp_sizes);
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
  if (kls_fragmented_metis_tuning_ctx) {
    /* Fragmented BTF cores use one separator trial.  Three refinement
       iterations cut retained fill on the broad class; the dense reciprocal
       mega-hub proposal uses ten to obtain a stable native separator tree
       without a second CAMD ordering pass.  Environment overrides remain
       authoritative for experiments. */
    options[METIS_OPTION_NITER] =
      kls_fragmented_metis_tuning_ctx == 2 ? 10 : 3;
    options[METIS_OPTION_NSEPS] = 1;
  } else

  {
    const char *niter = getenv("KLS_METIS_NITER");
    if (niter != NULL && niter[0] != '\0') {
      options[METIS_OPTION_NITER] = atoi(niter);
    }
    const char *nseps = getenv("KLS_METIS_NSEPS");
    if (nseps != NULL && nseps[0] != '\0') {
      options[METIS_OPTION_NSEPS] = atoi(nseps);
    } else if (kls_auto_metis_promotion_quality) {
      /* Expensive giant refactors benefit from the wider four-separator
         search.  Compact, recurring KLU factors have a different economics:
         two separator trials find the useful alternate tree without the
         over-refinement and setup cost of the giant profile. */
      options[METIS_OPTION_NSEPS] =
        kls_auto_metis_promotion_quality == 2 ? 2 : 4;
    }
    if (kls_auto_metis_promotion_quality &&
        getenv("KLS_METIS_NITER") == NULL) {
      options[METIS_OPTION_NITER] = 10;
    }
    const char *ufactor = getenv("KLS_METIS_UFACTOR");
    if (ufactor != NULL && ufactor[0] != '\0') {
      options[METIS_OPTION_UFACTOR] = atoi(ufactor);
    }
    const char *ctype = getenv("KLS_METIS_CTYPE");
    if (ctype != NULL && ctype[0] == 'r') {
      options[METIS_OPTION_CTYPE] = METIS_CTYPE_RM;
    }
    const char *dbg = getenv("KLS_METIS_DBGLVL");
    if (dbg != NULL && dbg[0] != '\0') {
      options[METIS_OPTION_DBGLVL] = atoi(dbg);
    }
  }
  idx_t nvtxs = (idx_t)n;
  const double kls_metis_t0 = kls_now_seconds();
  int metis_status = METIS_ERROR;
  int par_nd_done = 0;
#ifdef KLS_HAVE_MTMETIS
  if (n >= 200000 &&
      (getenv("KLS_MT_ND") != NULL ||
       (kls_mt_nd_class && getenv("KLS_MT_ND_CLASS") != NULL))) {
    /* threaded ND (mt-metis, MIT).  With a leaf-count request the
       forest-aware path replicates NodeNDP's sizes contract so the
       separator forest (and the pipeline refactor engines that read
       it) stay intact; otherwise the plain bridge. */
    if (metis_ndp_npes > 1 && metis_ndp_sizes != NULL) {
      par_nd_done = kls_mtmetis_ndp(n, xadj, adjncy, metis_ndp_npes,
                                    metis_perm, metis_iperm,
                                    metis_ndp_sizes, 0);
    } else {
      par_nd_done = kls_mtmetis_nd(n, xadj, adjncy, metis_perm,
                                   metis_iperm);
      if (par_nd_done) {
        metis_ndp_npes = 0;
      }
    }
    if (par_nd_done) {
      metis_status = METIS_OK;
    }
  }
#endif
  if (metis_ndp_npes <= 1 && n >= 200000 &&
      getenv("KLS_PAR_ND") != NULL) {
    /* one manual top level of nested dissection: METIS's own vertex
       separator, then the two halves NodeND'd CONCURRENTLY (NodeND is
       serial; this is its own recursion with the top level threaded).
       Composed order: half 0, half 1, separator last - the same shape
       NodeND produces.  Quality is arbitrated downstream by the
       ordering competition's fill scores. */
    par_nd_done = kls_metis_parallel_nd(nvtxs, xadj, adjncy, options,
                                        metis_perm, metis_iperm);
    if (par_nd_done) {
      metis_status = METIS_OK;
    }
  }
#ifdef KLS_HAVE_MTMETIS
  const int force_det_ndp =
    getenv("KLS_FORCE_DET_NDP") != NULL ||
    kls_medium_partial_static_metis_ctx ||
    kls_generic_metis_tuning_ctx ||
    kls_asymmetric_bounded_degree_direct_metis_analyze_class;
  if (!par_nd_done && metis_ndp_npes > 1 && metis_ndp_sizes != NULL &&
      n >= 50000 &&
      (force_det_ndp || kls_spiked_ndp_class ||
       (kls_metis_order_threads >= (int)metis_ndp_npes &&
        (kls_det_ndp_requested || kls_lowdeg_ndp_class ||
         kls_detndp_class_ok))) &&
      getenv("KLS_DISABLE_DET_NDP") == NULL) {
    /* deterministic parallel NodeNDP: serial METIS separators with
       the leaves ordered by concurrent serial NodeND.  Deterministic
       because KLS_GKRAND (default ON) gives each METIS call its own
       seeded thread-local MT19937 stream; the concurrent leaves
       cannot perturb each other.  Scoped to callers that request it
       (the pre-static METIS refinement, where it measured a win):
       v11 measured pre2's separator-pipeline refactor 5x slower when
       its analyze-time forest came from this splitter instead of
       serial NodeNDP - the same structure-sensitivity that rejected
       mt-metis forests there.  Ordering quality is arbitrated
       downstream by the candidate competition's measured fill/flops
       gates. */
    par_nd_done = kls_mtmetis_ndp(n, xadj, adjncy, metis_ndp_npes,
                                  metis_perm, metis_iperm,
                                  metis_ndp_sizes, 1);
    if (par_nd_done) {
      metis_status = METIS_OK;
    }
  }
#endif
  if (!par_nd_done) {
    metis_status =
      metis_ndp_npes > 1
        ? METIS_NodeNDP(nvtxs, xadj, adjncy, NULL, metis_ndp_npes, options,
                        metis_perm, metis_iperm, metis_ndp_sizes)
        : METIS_NodeND(&nvtxs, xadj, adjncy, NULL, options,
                       metis_perm, metis_iperm);
  }
  if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
    fprintf(stderr, "KLS metis: NodeND %.3fs (ndp=%d par=%d tid=%ld)\n",
            kls_now_seconds() - kls_metis_t0, (int)metis_ndp_npes,
            (int)par_nd_done, (long)syscall(SYS_gettid));
  }
  UF_long order_lnz = 0;
  if (metis_status == METIS_OK) {
    UF_long camd_group_size = kls_medium_partial_static_metis_ctx
      ? 768u : kls_metis_camd_group_size(n, col_ptr, row_idx);
    if (kls_generic_metis_tuning_ctx) {
      const char *generic_group_env =
        getenv("KLS_GENERIC_METIS_CAMD_GROUP_SIZE");
      if (generic_group_env != NULL && generic_group_env[0] != '\0') {
        const unsigned long long requested =
          strtoull(generic_group_env, NULL, 10);
        if (requested > 0u && requested <= (unsigned long long)n) {
          camd_group_size = (UF_long)requested;
        }
      } else {
        /* A BTF challenger spends its ND pass independently inside each
           irreducible block, while no-BTF uses a coarser whole-graph budget.
           Both variants are selected from executor resources and measured
           symbolic economics rather than an input-family classifier. */
        if (kls_generic_metis_tuning_ctx == 5 ||
            kls_generic_metis_tuning_ctx == 6 ||
            kls_generic_metis_tuning_ctx == 7 ||
            kls_generic_metis_tuning_ctx == 8) {
          /* The generic ordering portfolio exposes two separator leaves and
             cache-bounded refinement windows scaled only by the caller's
             execution resources.  A lower-arithmetic-intensity incumbent
             uses context six's finer windows because each removed factor
             entry is a recurring memory-stream saving; context five keeps
             coarser windows where separator execution work dominates.  The
             symbolic/numeric acceptance gates still decide whether the
             resulting ordering survives. */
          const UF_long workers = kls_metis_order_threads > 0
            ? (UF_long)kls_metis_order_threads : 1u;
          const UF_long entries_per_worker =
            kls_generic_metis_tuning_ctx == 6 ||
            kls_generic_metis_tuning_ctx == 8 ? 128u : 512u;
          camd_group_size = workers <= UF_long_max / entries_per_worker
            ? entries_per_worker * workers : UF_long_max;
          const UF_long minimum =
            kls_generic_metis_tuning_ctx == 6 ||
            kls_generic_metis_tuning_ctx == 8 ? 512u : 1024u;
          if (camd_group_size < minimum) {
            camd_group_size = minimum;
          }
          if (camd_group_size > 8192u) {
            camd_group_size = 8192u;
          }
          if (camd_group_size > n) {
            camd_group_size = n;
          }
        } else if (kls_generic_metis_tuning_ctx == 3) {
          /* A realized one-block BTF representation has no SCC-level
             parallelism to preserve.  Give NodeNDP one leaf per numeric
             worker and refine in a cache-sized window; the subsequent
             symbolic and numeric gates remain the ordering authorities. */
          camd_group_size = n < 1024u ? n : 1024u;
        } else if (kls_generic_metis_tuning_ctx == 2) {
          /* SCC boundaries already provide coarse outer tasks.  Give the
             dominant block one NodeNDP leaf per requested numeric worker,
             then retain its forest in a cache-bounded CAMD tile. */
          camd_group_size = n < 2304u ? n : 2304u;
        } else {
        /* Row-scaled recurring factors pay an additional factor-sized memory
           stream on every update and solve, so spend more ordering resources
           trimming those streams (about six windows per worker).  An
           unscaled factor instead preserves a coarser separator forest (about
           two windows per worker).  The distinction is the ordinary
           current-value AUTO-scale verdict, not a pattern or dimension
           classifier.  Fixed resource clamps bound both variants. */
        const UF_long workers = kls_metis_order_threads > 0
          ? (UF_long)kls_metis_order_threads : 1u;
        const UF_long windows_per_worker =
          kls_generic_metis_scale_hint > 0 ? 6u : 2u;
        const UF_long divisor = workers <=
            UF_long_max / windows_per_worker
          ? windows_per_worker * workers : UF_long_max;
        camd_group_size = divisor > 0u ? n / divisor : n;
        const UF_long unit =
          kls_generic_metis_scale_hint > 0 ? 256u : 1024u;
        if (kls_generic_metis_scale_hint > 0) {
          camd_group_size = (camd_group_size / unit) * unit;
        } else if (camd_group_size <= UF_long_max - (unit - 1u)) {
          /* Preserve all requested coarse windows: an undersized final
             unscaled group adds a separator boundary without useful work. */
          camd_group_size =
            ((camd_group_size + unit - 1u) / unit) * unit;
        }
        const UF_long minimum =
          kls_generic_metis_scale_hint > 0 ? 512u : 1024u;
        const UF_long maximum =
          kls_generic_metis_scale_hint > 0 ? 4096u : 8192u;
        if (camd_group_size < minimum) {
          camd_group_size = minimum;
        }
        if (camd_group_size > maximum) {
          camd_group_size = maximum;
        }
        }
      }
    }
    if (camd_group_size == 0 && n > 1) {
      /* Below the window-refinement floor, keep the METIS order exactly by
         constraining every vertex to its own group.  CAMD then reproduces
         the given elimination order while still reporting its fill
         estimate, so METIS-ordered symbolic candidates get a real lnz and
         can be compared against AMD/COLAMD scores instead of reporting an
         unknown estimate for the whole ordering. */
      camd_group_size = 1;
    }
    const double kls_camd_t0 = kls_now_seconds();
    int pinned_count_ok = 0;
    if (camd_group_size == 1) {
      /* the fully-pinned "refinement" only produces an lnz estimate.
         On patterns WITHOUT dense rows CAMD's degree machinery tracks
         the pinned order's fill and grinds (ss1's 109K block: 19-23s;
         the etree count is exact and O(min(fill, 8*nnz))).  Patterns
         WITH dense rows keep CAMD: it is fast there (it postpones the
         dense rows) and the block-structured pipeline's downstream
         gates are tuned to its estimates (TSOPF_FS_b39_c30 measured
         0.67 -> 6.5s under the honest count). */
      UF_long *rowdeg = (UF_long *)calloc((size_t)n, sizeof(*rowdeg));
      if (rowdeg != NULL) {
        UF_long maxdeg = 0;
        for (UF_long p = 0; p < col_ptr[n]; ++p) {
          rowdeg[row_idx[p]]++;
        }
        for (UF_long i = 0; i < n; ++i) {
          const UF_long d =
            rowdeg[i] + (col_ptr[i + 1u] - col_ptr[i]);
          maxdeg = maxdeg < d ? d : maxdeg;
        }
        free(rowdeg);
        pinned_count_ok = (double)maxdeg * (double)maxdeg < (double)n;
      }
    }
    const int component_refine =
      getenv("KLS_METIS_CAMD_COMPONENT_GROUPS") != NULL &&
      metis_ndp_npes > 1 && metis_ndp_sizes != NULL;
    const int skip_camd_refine =
      getenv("KLS_METIS_SKIP_CAMD_REFINE") != NULL ||
      (kls_fragmented_metis_tuning_ctx == 2 &&
       getenv("KLS_ENABLE_DENSE_RECIPROCAL_HUB_CAMD_REFINE") == NULL);
    UF_long camd_lnz;
    if (skip_camd_refine) {
      /* The pre-static raced-symbolic path subsequently constructs and
         probes the exact LU pattern, so it does not need CAMD's fill estimate.
         The dense reciprocal-hub policy also selects NodeNDP directly after
         symbolic acceptance.  Keeping NodeNDP's order avoids a second serial
         ordering pass in both cases; ordinary candidate selection still uses
         CAMD's estimate. */
      for (UF_long i = 0; i < n; ++i) {
        perm_out[i] = (UF_long)metis_perm[i];
      }
      camd_lnz = (UF_long)-1;
    } else {
      camd_lnz = component_refine
        ? kls_metis_refine_components_with_camd(
            n, col_ptr, row_idx, metis_perm, metis_iperm, metis_ndp_npes,
            metis_ndp_sizes, perm_out)
        : pinned_count_ok
            ? kls_pinned_order_fill_count(n, col_ptr, row_idx, metis_perm,
                                          perm_out)
            : kls_metis_refine_with_camd(n, col_ptr, row_idx, metis_perm,
                                         camd_group_size, perm_out);
    }
    if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
      fprintf(stderr,
              "KLS metis: camd-refine %.3fs (n=%ld nnz=%ld groups=%s%ld)\n",
              kls_now_seconds() - kls_camd_t0,
              (long)n, (long)col_ptr[n],
              skip_camd_refine ? "skipped/" :
                (component_refine ? "components/" : ""),
              skip_camd_refine ? 0l :
                (component_refine ? (long)(2 * metis_ndp_npes - 1)
                                  : (long)camd_group_size));
    }
    if (camd_lnz == 0) {
      for (UF_long i = 0; i < n; ++i) {
        perm_out[i] = (UF_long)metis_perm[i];
      }
      camd_lnz = (UF_long)-1;
    }
    order_lnz = camd_lnz;
    if (metis_context != NULL && metis_context->separator != NULL &&
        metis_ndp_npes > 1 && metis_ndp_sizes != NULL) {
      kls_separator_analysis separator;
      memset(&separator, 0, sizeof(separator));
      if (kls_build_metis_separator_analysis(n, metis_ndp_npes,
                                             metis_ndp_sizes, metis_iperm,
                                             perm_out, &separator)) {
        (void)kls_metis_order_context_append_separator(
          metis_context, order_call_index, &separator);
        kls_separator_analysis_clear(&separator);
      }
    }
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
  free(metis_ndp_sizes);
  return metis_status == METIS_OK ? order_lnz : 0;
}

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

/* precomputed-BTF stash: the block triangular form depends only on the
   pattern, so one competition's candidates (AMD, AMF, METIS, retries)
   share it through the vendored Common instead of each re-running
   maxtrans+strongcomp (twotone: 0.75-1.46s per with-btf analyze).
   Thread-local: orientation sides and the AMF spec thread each carry
   their own. */
struct kls_btf_stash_s {
  const UF_long *P;
  const UF_long *Q;
  const UF_long *R;
  UF_long nblocks;
  UF_long rank;
  UF_long n;
  int active;
};
static _Thread_local struct kls_btf_stash_s kls_btf_stash;

static int analyze_with_ordering(UF_long n,
                                 UF_long *col_ptr,
                                 UF_long *row_idx,
                                 const kls_options *options,
                                 kls_ordering ordering,
                                 trilinos_klu_l_symbolic **symbolic_out,
                                 trilinos_klu_l_common *common_out,
                                 kls_separator_analysis *separator_out) {
  kls_separator_analysis_clear(separator_out);
#ifdef KLS_HAVE_MTMETIS
  kls_metis_order_threads = options != NULL ? options->threads : 1;
#endif
  {
    /* the large-spral separator-pipeline class (pre2): its refactor
       engine is structure-sensitive to the exact serial-NodeNDP
       forests the pre-GKRAND libc RNG produced (v11: 0.17s -> 7-9.6s
       refactor under MT19937 orderings), so this class keeps the
       legacy RNG for its ordering calls.  Thread-local: set per
       analyze call from the same pattern bounds the prestatic class
       gate uses. */
    extern void gk_set_legacy_rand(int enable);
    gk_set_legacy_rand(0);
#ifdef KLS_HAVE_MTMETIS
    kls_lowdeg_ndp_class = 0;
    kls_detndp_class_ok = 1;
    kls_spiked_ndp_class = 0;
#endif
  }
  const double kls_awo_t0 =
    getenv("KLS_TRACE_ANALYZE_STAGES") != NULL ? kls_now_seconds() : 0.0;
  trilinos_klu_l_common common;
  int status = apply_options_to_common(&common, options);
  if (status != KLS_OK) {
    return status;
  }
  if (kls_btf_stash.active && kls_btf_stash.n == n &&
      options != NULL && options->use_btf) {
    common.kls_btf_P = kls_btf_stash.P;
    common.kls_btf_Q = kls_btf_stash.Q;
    common.kls_btf_R = kls_btf_stash.R;
    common.kls_btf_nblocks = kls_btf_stash.nblocks;
    common.kls_btf_structural_rank = kls_btf_stash.rank;
  }

  kls_metis_order_context metis_context;
  memset(&metis_context, 0, sizeof(metis_context));
  int metis_context_active = 0;

  trilinos_klu_l_symbolic *symbolic = NULL;
  if (ordering == KLS_ORDERING_NATURAL) {
    symbolic = trilinos_klu_l_analyze_given(n, col_ptr, row_idx, NULL, NULL,
                                            &common);
  } else if (ordering == KLS_ORDERING_METIS) {
    metis_context.n = n;
    metis_context.npes =
      options != NULL
        ? kls_metis_ndp_leaf_count(n, options->threads,
                                   kls_first_factor_env_enabled())
        : 0;
    metis_context.separator = separator_out;
    metis_context_active = 1;
    common.ordering = 3;
    common.user_order = kls_metis_order;  /* timing wrapper */
    common.user_data = &metis_context;
    const double kls_aa_t0 =
      getenv("KLS_TRACE_ANALYZE_STAGES") != NULL ? kls_now_seconds() : 0.0;
    symbolic = trilinos_klu_l_analyze(n, col_ptr, row_idx, &common);
    if (getenv("KLS_TRACE_ANALYZE_STAGES") != NULL) {
      fprintf(stderr, "KLS awo: klu_analyze %.3fs\n",
              kls_now_seconds() - kls_aa_t0);
    }
    common.user_data = NULL;
  } else if (ordering == KLS_ORDERING_SCOTCH) {
#ifdef KLS_HAVE_SCOTCH
    common.ordering = 3;
    common.user_order = kls_scotch_order;
    symbolic = trilinos_klu_l_analyze(n, col_ptr, row_idx, &common);
#else
    return KLS_ERR_UNSUPPORTED;
#endif
  } else if (ordering == KLS_ORDERING_AMF ||
             ordering == KLS_ORDERING_AMMF ||
             ordering == KLS_ORDERING_AMF3) {
    /* vendored AMD with approximate-deficiency pivot scoring; the flag
       is a serial-analysis global (the ordering competition never runs
       analyses concurrently) */
    extern _Thread_local UF_long trilinos_amd_l2_amf;
    extern _Thread_local double trilinos_amd_l2_amf_power;
    {
      /* mean-local-fill scoring (Rothberg-Eisenstat AMMF) wins on the
         small classes (TSOPF est 1.31e5 vs METIS 1.50e5) but its lower
         fill estimates SERIALIZE the pipelined refactor at scale
         (ASIC_320ks: AMMF bootstrap beat the race-METIS promotion on
         est flops 6.8e8 vs 9.4e8 yet refactors at 66 vs 43ms, suite
         2.42x) - keep the old scoring above the small-matrix line
         until ordering acceptance is timing-aware. */
      const char *amf_mode = getenv("KLS_AMF_MODE");
      trilinos_amd_l2_amf_power = 1.0;
      if (ordering == KLS_ORDERING_AMMF) {
        trilinos_amd_l2_amf = 2;
        const char *power = getenv("KLS_AMMF_POWER");
        if (power != NULL && power[0] != '\0') {
          const double requested = strtod(power, NULL);
          if (isfinite(requested) && requested >= 0.0 && requested <= 2.0) {
            trilinos_amd_l2_amf_power = requested;
          }
        }
      } else if (ordering == KLS_ORDERING_AMF3) {
        trilinos_amd_l2_amf = kls_amf3_exact_clique_bound ? 4 : 3;
        trilinos_amd_l2_amf_power = 0.445;
        if (kls_amf3_power_override_active &&
            isfinite(kls_amf3_power_override) &&
            kls_amf3_power_override >= 0.0 &&
            kls_amf3_power_override <= 2.0) {
          trilinos_amd_l2_amf_power = kls_amf3_power_override;
        }
        const char *power = getenv("KLS_AMF3_POWER");
        if (power != NULL && power[0] != '\0') {
          const double requested = strtod(power, NULL);
          if (isfinite(requested) && requested >= 0.0 && requested <= 2.0) {
            trilinos_amd_l2_amf_power = requested;
          }
        }
      } else if (amf_mode != NULL && amf_mode[0] != '\0') {
        const long requested = strtol(amf_mode, NULL, 10);
        trilinos_amd_l2_amf = requested == 1 || requested == 3
          ? (UF_long)requested : 2;
        if (trilinos_amd_l2_amf == 3) {
          trilinos_amd_l2_amf_power = 0.445;
        }
      } else {
        trilinos_amd_l2_amf = n <= 30000 ? 2 : 1;
        trilinos_amd_l2_amf_power = 1.0;
      }
    }
    common.ordering = 0;
    symbolic = trilinos_klu_l_analyze(n, col_ptr, row_idx, &common);
    trilinos_amd_l2_amf = 0;
    trilinos_amd_l2_amf_power = 1.0;
  } else {
    const double kls_ad_t0 =
      getenv("KLS_TRACE_ANALYZE_STAGES") != NULL ? kls_now_seconds() : 0.0;
    common.ordering = (ordering == KLS_ORDERING_COLAMD) ? 1 : 0;
    symbolic = trilinos_klu_l_analyze(n, col_ptr, row_idx, &common);
    if (getenv("KLS_TRACE_ANALYZE_STAGES") != NULL) {
      fprintf(stderr, "KLS awo: amd_analyze %.3fs\n",
              kls_now_seconds() - kls_ad_t0);
    }
  }
  if (symbolic == NULL || common.status < 0) {
    if (symbolic != NULL) {
      trilinos_klu_l_free_symbolic(&symbolic, &common);
    }
    kls_separator_analysis_clear(separator_out);
    if (metis_context_active) {
      kls_metis_order_context_clear(&metis_context);
    }
    return KLS_ERR_ANALYZE_FAILED;
  }

  if (metis_context_active) {
    const double kls_sep_t0 = kls_now_seconds();
    if (getenv("KLS_TRACE_ANALYZE_STAGES") != NULL) {
      fprintf(stderr, "KLS metis: pre-forest %.3fs\n",
              kls_sep_t0 - kls_awo_t0);
    }
    if (!kls_build_metis_separator_forest(symbolic, &metis_context,
                                          separator_out)) {
      kls_metis_order_context_move_largest_separator(&metis_context,
                                                     separator_out);
      kls_finalize_separator_global_range(symbolic, separator_out);
    }
    kls_metis_order_context_clear(&metis_context);
    if (getenv("KLS_TRACE_ANALYZE_STAGES") != NULL) {
      fprintf(stderr, "KLS metis: forest %.3fs\n",
              kls_now_seconds() - kls_sep_t0);
    }
  } else
  {
    kls_finalize_separator_global_range(symbolic, separator_out);
  }
  if (getenv("KLS_TRACE_ANALYZE_STAGES") != NULL) {
    fprintf(stderr, "KLS awo: ord=%d btf=%d n=%ld %.3fs\n",
            (int)ordering, options != NULL && options->use_btf ? 1 : 0,
            (long)n, kls_now_seconds() - kls_awo_t0);
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

static int symbolic_is_high_work_dominant_btf_nd_candidate(
  UF_long n,
  const trilinos_klu_l_symbolic *symbolic) {
  if (symbolic == NULL || !symbolic->do_btf ||
      n < 500000u || symbolic->nblocks < 4096u ||
      symbolic->maxblock < (UF_long)(0.95 * (double)n) ||
      symbolic->est_flops < 5.0e10) {
    return 0;
  }
  return !symbolic_is_low_work_dominant_btf(n, symbolic);
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

/* Once a no-btf retry loses DECISIVELY on this pattern (>2x the btf
   score against a 0.80x adoption bar), further retries with other
   orderings are hopeless - the missing block structure, not the
   ordering, is what loses (twotone: both retries lose 6.7x, one per
   candidate, ~40% of each analyze side).  Reset per analyzed pattern;
   thread-local so the two orientation sides stay independent. */
static _Thread_local int kls_no_btf_retry_hopeless;
/* Each orientation owns its kls_pattern_candidate, and generic orientation
   analyses may run concurrently.  A thread-local capture target lets the
   representation comparison retain one bounded BTF alternative without
   sharing ownership across those independent tournaments. */
static _Thread_local kls_pattern_candidate *kls_generic_btf_capture_candidate;

/* A no-BTF symbolic is independent of the BTF incumbent it challenges.  A
   repeated-update caller can therefore build the two representations in
   parallel and feed the completed challenger through the ordinary acceptance
   test below.  The job carries no policy: it only computes the same
   analyze_with_ordering result that maybe_retry_without_btf would otherwise
   build serially. */
typedef struct kls_no_btf_spec_job {
  UF_long n;
  UF_long *col_ptr;
  UF_long *row_idx;
  kls_ordering ordering;
  kls_options options;
  trilinos_klu_l_symbolic *symbolic;
  trilinos_klu_l_common common;
  kls_separator_analysis separator;
  int status;
} kls_no_btf_spec_job;

static void *kls_no_btf_spec_main(void *arg) {
  kls_no_btf_spec_job *job = (kls_no_btf_spec_job *)arg;
  job->symbolic = NULL;
  memset(&job->separator, 0, sizeof(job->separator));
  job->options.use_btf = 0;
  job->status = analyze_with_ordering(
    job->n, job->col_ptr, job->row_idx, &job->options, job->ordering,
    &job->symbolic, &job->common, &job->separator);
  return NULL;
}

static void kls_no_btf_spec_clear(kls_no_btf_spec_job *job) {
  if (job == NULL) {
    return;
  }
  if (job->symbolic != NULL) {
    trilinos_klu_l_free_symbolic(&job->symbolic, &job->common);
  }
  kls_separator_analysis_clear(&job->separator);
}

/* Static-pattern prediction needs A(P[k],Q[k]) in every pivot position.
   This O(nnz) certificate is cheaper than building the symmetrized numeric
   pattern and describes an executor capability of the analyzed
   representation; it is independent of input dimensions or family. */
static int kls_symbolic_has_complete_pivot_diagonal(
  UF_long n,
  const UF_long *col_ptr,
  const UF_long *row_idx,
  const trilinos_klu_l_symbolic *symbolic) {
  if (n == 0u || col_ptr == NULL || row_idx == NULL || symbolic == NULL ||
      symbolic->P == NULL || symbolic->Q == NULL) {
    return 0;
  }
  for (UF_long k = 0u; k < n; ++k) {
    const UF_long oldrow = symbolic->P[k];
    const UF_long oldcol = symbolic->Q[k];
    if (oldrow >= n || oldcol >= n) {
      return 0;
    }
    int found = 0;
    for (UF_long p = col_ptr[oldcol]; p < col_ptr[oldcol + 1u]; ++p) {
      if (row_idx[p] == oldrow) {
        found = 1;
        break;
      }
    }
    if (!found) {
      return 0;
    }
  }
  return 1;
}

static int kls_symbolic_repeated_prediction_economic(
  const trilinos_klu_l_symbolic *symbolic,
  const kls_options *options) {
  if (symbolic == NULL || options == NULL || options->threads <= 1 ||
      !kls_repeated_update_workload(options) ||
      !(symbolic->est_flops > 0.0) || !isfinite(symbolic->est_flops)) {
    return 0;
  }
  const double lifecycle_work = symbolic->est_flops *
    (1.0 + (double)options->expected_refactorizations);
  return isfinite(lifecycle_work) && lifecycle_work >= 2.0e10;
}

static void maybe_retry_without_btf(UF_long n,
                                    UF_long *col_ptr,
                                    UF_long *row_idx,
                                    const kls_options *options,
                                    kls_ordering ordering,
                                    trilinos_klu_l_symbolic **symbolic,
                                    trilinos_klu_l_common *common,
                                    double *score,
                                    int allow_single_block,
                                    kls_separator_analysis *separator_io,
                                    kls_no_btf_spec_job *precomputed) {
  if (options == NULL || !options->use_btf || n < 4000 || symbolic == NULL ||
      *symbolic == NULL || !(*symbolic)->do_btf || common == NULL ||
      score == NULL ||
      kls_no_btf_retry_hopeless) {
    return;
  }

  const int single_block =
    n >= 12000u && allow_single_block &&
    (*symbolic)->nblocks == 1 && (*symbolic)->maxblock == n;
  /* A near-degenerate BTF cannot normally move the fill estimate by the 20%
     the historical dominant-block retry demands: with the dominant block
     holding all but a fraction f of the rows, the no-BTF ordering problem
     differs by O(f).  A repeated generic lifecycle is allowed to challenge
     that assumption only when the retained symbolic work, multiplied by the
     requested numeric horizon, can repay another complete analysis.  The
     challenger must still improve both estimated work and non-strict fill. */
  const int btf_degenerate =
    (*symbolic)->nblocks > 1 && (*symbolic)->maxblock >= n - n / 8;
  const double generic_numeric_horizon = options != NULL
    ? 1.0 + (double)options->expected_refactorizations
    : 1.0;
  const int ready_precomputed_comparison =
    precomputed != NULL && precomputed->ordering == ordering &&
    precomputed->status == KLS_OK && precomputed->symbolic != NULL;
  const double degenerate_fringe_fraction = btf_degenerate && n > 0u
    ? (double)(n - (*symbolic)->maxblock) / (double)n : 1.0;
  const double degenerate_change_budget =
    (*symbolic)->est_flops * generic_numeric_horizon *
      degenerate_fringe_fraction;
  const int generic_degenerate_comparison =
    btf_degenerate &&
    kls_repeated_update_workload(options) &&
    (*symbolic)->est_flops > 0.0 &&
    isfinite(generic_numeric_horizon) && generic_numeric_horizon > 1.0 &&
    (ready_precomputed_comparison ||
     (isfinite(degenerate_change_budget) &&
      degenerate_change_budget >= 1.0e8));
  const int generic_comparison = ((!btf_degenerate &&
      ((*symbolic)->nblocks > 1u || single_block)) ||
     generic_degenerate_comparison);
  if (!single_block && !generic_comparison) {
    return;
  }

  if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
    fprintf(stderr,
            "KLS retry-no-btf: nblocks=%ld maxblock=%ld/%ld single=%d"
            " dominant=%d inflated=%d fragmented=%d low_work=%d\n",
            (long)(*symbolic)->nblocks, (long)(*symbolic)->maxblock, (long)n,
            single_block, 0, 0,
            0, 0);
  }
  kls_options no_btf_options = *options;
  no_btf_options.use_btf = 0;
  trilinos_klu_l_symbolic *no_btf_symbolic = NULL;
  trilinos_klu_l_common no_btf_common;
  kls_separator_analysis no_btf_separator;
  memset(&no_btf_separator, 0, sizeof(no_btf_separator));
  int status = KLS_ERR_ANALYZE_FAILED;
  if (precomputed != NULL && precomputed->ordering == ordering &&
      precomputed->status == KLS_OK && precomputed->symbolic != NULL) {
    no_btf_symbolic = precomputed->symbolic;
    precomputed->symbolic = NULL;
    no_btf_common = precomputed->common;
    kls_separator_analysis_move(&no_btf_separator, &precomputed->separator);
    status = KLS_OK;
  } else {
    status = analyze_with_ordering(
      n, col_ptr, row_idx, &no_btf_options, ordering, &no_btf_symbolic,
      &no_btf_common, &no_btf_separator);
  }
  if (status != KLS_OK) {
    return;
  }

  const double no_btf_score = symbolic_score(no_btf_symbolic);
  const double current_score = *score;
  if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
    fprintf(stderr, "KLS retry-no-btf: score %.4e vs %.4e, "
                    "flops %.4e vs %.4e (need 0.80x=%d)\n",
            no_btf_score, current_score, no_btf_symbolic->est_flops,
            (*symbolic)->est_flops,
            no_btf_score <= 0.80 * current_score);
  }
  const int current_score_known =
    isfinite(current_score) && current_score < DBL_MAX / 4.0;
  const int repeated_generic_comparison =
    generic_comparison && kls_repeated_update_workload(options);
  const int drops_repeated_prediction = generic_comparison &&
    kls_symbolic_repeated_prediction_economic(*symbolic, options) &&
    kls_symbolic_has_complete_pivot_diagonal(
      n, col_ptr, row_idx, *symbolic) &&
    !kls_symbolic_has_complete_pivot_diagonal(
      n, col_ptr, row_idx, no_btf_symbolic);
  const int decisive_without_prediction =
    drops_repeated_prediction && current_score_known &&
    isfinite(no_btf_score) && no_btf_score <= 0.80 * current_score &&
    no_btf_symbolic->est_flops > 0.0 &&
    (*symbolic)->est_flops > 0.0 &&
    no_btf_symbolic->est_flops <= 0.75 * (*symbolic)->est_flops;
  const int pivoted_representation_pareto =
    getenv("KLS_DISABLE_GENERIC_PIVOTED_NO_BTF_REPRESENTATION") == NULL &&
    drops_repeated_prediction && repeated_generic_comparison &&
    current_score_known && isfinite(no_btf_score) &&
    no_btf_score <= 0.90 * current_score &&
    no_btf_symbolic->est_flops > 0.0 &&
    (*symbolic)->est_flops > 0.0 &&
    no_btf_symbolic->est_flops <= (*symbolic)->est_flops;
  const int generic_accept = generic_comparison && current_score_known &&
    isfinite(no_btf_score) &&
    (!drops_repeated_prediction || decisive_without_prediction ||
     pivoted_representation_pareto) &&
    (single_block
       ? no_btf_score <= 1.02 * current_score
       : no_btf_score <=
           (repeated_generic_comparison ? 1.00 : 0.90) * current_score &&
         (!(no_btf_symbolic->est_flops > 0.0) ||
          !((*symbolic)->est_flops > 0.0) ||
          no_btf_symbolic->est_flops <= (*symbolic)->est_flops));
  /* A non-degenerate BTF splits a material fraction of the recurring
     numeric into block/off-block executors.  When the ordinary no-BTF
     symbolic stays inside a narrow work-and-storage Pareto band, the
     caller's long update horizon can repay that representation without
     identifying the graph by name or dimensions.  Retain the BTF symbolic
     as a value-certified fallback: the first numeric can still recover it
     when its pivots are decisively stronger.  Near-degenerate BTFs remain
     on the stricter estimate rule because their tiny fringe cannot justify
     this executor change and symbolic work poorly predicts pivot fill there. */
  const int amortized_non_degenerate_no_btf =
    generic_comparison && repeated_generic_comparison && !btf_degenerate &&
    kls_generic_btf_capture_candidate != NULL &&
    getenv("KLS_DISABLE_GENERIC_BTF_VALUE_SELECTION") == NULL &&
    getenv("KLS_DISABLE_AMORTIZED_NONDEGENERATE_NO_BTF") == NULL &&
    current_score_known && current_score > 0.0 &&
    isfinite(no_btf_score) && no_btf_score > 0.0 &&
    (*symbolic)->nblocks > 1u && no_btf_symbolic->nblocks == 1u &&
    no_btf_score <= 1.15 * current_score &&
    (*symbolic)->est_flops > 0.0 &&
    no_btf_symbolic->est_flops > 0.0 &&
    no_btf_symbolic->est_flops <= 1.15 * (*symbolic)->est_flops;
  /* A tiny BTF fringe can cost more recurring executor dispatch than it
     saves in arithmetic.  Admit the no-BTF arm only when both symbolic
     portfolios fit a 64 MiB L/U value+index budget and a BTF fallback can be
     retained (or exact value matching has already certified the pivot
     frame).  Large factors stay on the stricter estimate rule because small
     symbolic errors there can expose substantial pivot fill. */
  const double btf_cache_entry_budget =
    (double)(64u * 1024u * 1024u) /
    (double)(sizeof(double) + sizeof(UF_long));
  const int amortized_cache_resident_tiny_fringe_no_btf =
    generic_comparison && repeated_generic_comparison &&
    ((*symbolic)->nblocks > 1u) &&
    (double)(n - (*symbolic)->maxblock) <= 0.02 * (double)n &&
    (kls_generic_btf_capture_candidate != NULL ||
     kls_value_matched_ordering_ctx) &&
    getenv("KLS_DISABLE_CACHE_RESIDENT_TINY_FRINGE_NO_BTF") == NULL &&
    current_score_known && current_score > 0.0 &&
    isfinite(no_btf_score) && no_btf_score > 0.0 &&
    current_score <= btf_cache_entry_budget &&
    no_btf_score <= btf_cache_entry_budget &&
    no_btf_symbolic->nblocks == 1u &&
    no_btf_score <= 1.15 * current_score &&
    (*symbolic)->est_flops > 0.0 &&
    no_btf_symbolic->est_flops > 0.0 &&
    no_btf_symbolic->est_flops <= 1.15 * (*symbolic)->est_flops;
  if (generic_accept || amortized_non_degenerate_no_btf ||
      amortized_cache_resident_tiny_fringe_no_btf ||
      (single_block && no_btf_score <= 1.02 * current_score)) {
    const int capture_btf_value_alternative =
      kls_generic_btf_capture_candidate != NULL &&
      kls_generic_btf_capture_candidate->generic_btf_value_symbolic == NULL &&
      getenv("KLS_DISABLE_GENERIC_BTF_VALUE_SELECTION") == NULL &&
      options->ordering == KLS_ORDERING_AUTO &&
      kls_repeated_update_workload(options) &&
      (*symbolic)->do_btf && (*symbolic)->nblocks > 1u &&
      (*symbolic)->maxblock < n &&
      ((double)(n - (*symbolic)->maxblock) >= 0.02 * (double)n ||
       amortized_cache_resident_tiny_fringe_no_btf) &&
      current_score_known && current_score > 0.0 &&
      no_btf_score > 0.0 && isfinite(no_btf_score) &&
      current_score <= 1.35 * no_btf_score &&
      (*symbolic)->est_flops > 0.0 &&
      no_btf_symbolic->est_flops > 0.0 &&
      (*symbolic)->est_flops <= 2.0 * no_btf_symbolic->est_flops;
    if (capture_btf_value_alternative) {
      kls_pattern_candidate *capture =
        kls_generic_btf_capture_candidate;
      capture->generic_btf_value_symbolic = *symbolic;
      capture->generic_btf_value_common = *common;
      capture->generic_btf_value_ordering = ordering;
      *symbolic = NULL;
      kls_separator_analysis_move(
        &capture->generic_btf_value_separator, separator_io);
    } else {
      trilinos_klu_l_free_symbolic(symbolic, common);
    }
    *symbolic = no_btf_symbolic;
    *common = no_btf_common;
    *score = no_btf_score;
    kls_separator_analysis_move(separator_io, &no_btf_separator);
    return;
  }

  if (current_score_known && isfinite(no_btf_score) &&
      no_btf_score > 2.0 * current_score) {
    kls_no_btf_retry_hopeless = 1;
  }
  trilinos_klu_l_free_symbolic(&no_btf_symbolic, &no_btf_common);
  kls_separator_analysis_clear(&no_btf_separator);
}

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
  if (is_large_moderate_degree_weak_diagonal_metis_pattern(n, col_ptr,
                                                           row_idx)) {
    return 1;
  }
  /* The 30K--45K dense-diagonal class is cheap enough that NodeND cannot
     repay its own latency.  Its high-work pivoted one-block AMD lifecycle
     beats the METIS choice, while the METIS analyze alone can cost several
     times the AMD numeric.  Keep the classifier for scale/BTF policy below,
     but let the ordinary minimum-degree competition select its ordering. */
  if (is_medium_spiked_low_diagonal_pattern(n, col_ptr, row_idx)) {
    return 1;
  }
  if (is_small_spiked_low_diagonal_pattern(n, col_ptr, row_idx)) {
    return 1;
  }
  return is_medium_bounded_degree_diagonal_pattern(n, col_ptr, row_idx);
}

static int metis_start_should_skip_no_btf_retry(UF_long n,
                                                const UF_long *col_ptr,
                                                const UF_long *row_idx) {
  /* Large full/nearly-full diagonal ASIC-style matrices use BTF to keep the
     off-diagonal fringe cheap. Retrying no-BTF doubles METIS analysis cost and
     has worse repeated-refactor work on this structural class. */
  return is_large_diagonal_metis_start_pattern(n, col_ptr, row_idx) ||
         is_large_moderate_degree_weak_diagonal_metis_pattern(n, col_ptr,
                                                              row_idx);
}

/* Set for the duration of the top-level analyze so the pattern-level
   ordering gates can see whether a background METIS race already covers
   the nested-dissection trial. */
static _Thread_local const kls_solver *kls_analyze_nd_race_solver;
/* Candidate-scoring phase of a two-orientation analyze: defer every
   nested-dissection trial; the chosen orientation re-analyzes once with
   the full path afterwards, so only the winner pays NodeND. */
static _Thread_local int kls_analyze_defer_nd;

static int should_try_symbolic_nested_dissection_before_numeric(
  UF_long n,
  const trilinos_klu_l_symbolic *symbolic,
  kls_ordering selected_ordering,
  kls_ordering trial_ordering,
  double score) {
  if (selected_ordering == trial_ordering || symbolic == NULL ||
      n < 30000u) {
    return 0;
  }
  if (kls_analyze_defer_nd) {
    return 0;
  }
  if (trial_ordering == KLS_ORDERING_METIS &&
      kls_analyze_nd_race_solver != NULL &&
      kls_analyze_nd_race_solver->metis_race != NULL) {
    /* The race computes the METIS ordering off the critical path and
       the factor-time promotion adopts it under timed acceptance. */
    return 0;
  }
  if (trial_ordering == KLS_ORDERING_METIS &&
      symbolic_is_high_work_dominant_btf_nd_candidate(n, symbolic)) {
    return isfinite(score) && score > 0.0;
  }
  /* A minimum-degree factor whose predicted arithmetic and storage dwarf
     the pattern can make even the *first* numeric effectively unbounded.
     In that regime, paying once for nested dissection is justified before
     values arrive.  This is an engine-cost test (Baumann: 44G predicted
     flops and 46M entries from a 112K-row pattern), not an input-family or
     dimension window. */
  const int predicted_factor_dominates_nd =
    trial_ordering == KLS_ORDERING_METIS &&
    symbolic->est_flops >= 1.0e10 &&
    isfinite(score) && score >= 64.0 * (double)n;
  if (predicted_factor_dominates_nd) {
    return 1;
  }
  if (n < 200000u) {
    return 0;
  }
  if (symbolic->do_btf || symbolic->nblocks != 1 ||
      symbolic->maxblock != n || symbolic->est_flops < 1.0e9) {
    return 0;
  }
  return isfinite(score) && score > 0.0;
}

static int should_accept_high_work_dominant_btf_nd_trial(
  UF_long n,
  const trilinos_klu_l_symbolic *best_symbolic,
  kls_ordering trial_ordering,
  const trilinos_klu_l_symbolic *trial_symbolic,
  const kls_separator_analysis *trial_separator) {
  if (trial_ordering != KLS_ORDERING_METIS ||
      !symbolic_is_high_work_dominant_btf_nd_candidate(n, best_symbolic) ||
      trial_symbolic == NULL || trial_separator == NULL ||
      !trial_symbolic->do_btf ||
      trial_symbolic->maxblock < (UF_long)(0.90 * (double)n) ||
      trial_separator->n != n ||
      trial_separator->pipeline_rows == 0u ||
      trial_separator->component_count == 0u) {
    return 0;
  }
  return 1;
}

static void maybe_promote_symbolic_ordering(
  UF_long n,
  UF_long *col_ptr,
  UF_long *row_idx,
  const kls_options *symbolic_options,
  kls_ordering trial_ordering,
  trilinos_klu_l_symbolic **best_symbolic_io,
  trilinos_klu_l_common *best_common_io,
  kls_ordering *best_ordering_io,
  double *best_score_io,
  double improvement_ratio,
  kls_separator_analysis *best_separator_io) {
  if (symbolic_options == NULL || best_symbolic_io == NULL ||
      *best_symbolic_io == NULL || best_common_io == NULL ||
      best_ordering_io == NULL || best_score_io == NULL ||
      !isfinite(*best_score_io) || *best_score_io <= 0.0 ||
      improvement_ratio <= 0.0) {
    return;
  }

  kls_options trial_options = *symbolic_options;
  trial_options.use_btf = (*best_symbolic_io)->do_btf ? 1 : 0;
  trilinos_klu_l_symbolic *trial_symbolic = NULL;
  trilinos_klu_l_common trial_common;
  kls_separator_analysis trial_separator;
  memset(&trial_separator, 0, sizeof(trial_separator));
  int status = analyze_with_ordering(n, col_ptr, row_idx, &trial_options,
                                     trial_ordering, &trial_symbolic,
                                     &trial_common, &trial_separator);
  if (status != KLS_OK) {
    return;
  }

  double trial_score = symbolic_score(trial_symbolic);

  const int score_accept =
    isfinite(trial_score) &&
    trial_score <= improvement_ratio * (*best_score_io);
  const int separator_accept =
    should_accept_high_work_dominant_btf_nd_trial(
      n, *best_symbolic_io, trial_ordering, trial_symbolic,
      &trial_separator);
  if (score_accept || separator_accept) {
    trilinos_klu_l_free_symbolic(best_symbolic_io, best_common_io);
    *best_symbolic_io = trial_symbolic;
    *best_common_io = trial_common;
    *best_ordering_io = trial_ordering;
    *best_score_io = trial_score;
    kls_separator_analysis_move(best_separator_io, &trial_separator);
    return;
  }

  trilinos_klu_l_free_symbolic(&trial_symbolic, &trial_common);
  kls_separator_analysis_clear(&trial_separator);
}

/* A generic AUTO ordering proposal can be computed beside the retained
   minimum-degree analysis.  The worker owns only its symbolic result; the
   canonical pattern is immutable until select_candidate returns.  Keeping
   this small race separate from the factor-time numeric race makes the
   acceptance rule identical to the synchronous portfolio below. */
typedef struct kls_generic_nd_trial {
  pthread_t thread;
  int active;
  UF_long n;
  UF_long *col_ptr;
  UF_long *row_idx;
  kls_options options;
  int scale_hint;
  int status;
  double seconds;
  trilinos_klu_l_symbolic *symbolic;
  trilinos_klu_l_common common;
  kls_separator_analysis separator;
} kls_generic_nd_trial;

static void *kls_generic_nd_trial_main(void *arg) {
  kls_generic_nd_trial *trial = (kls_generic_nd_trial *)arg;
  const int saved_tuning_ctx = kls_generic_metis_tuning_ctx;
  const int saved_scale_hint = kls_generic_metis_scale_hint;
  const struct kls_btf_stash_s saved_btf_stash = kls_btf_stash;
  const double start = kls_now_seconds();
  /* A complete no-BTF graph receives two deterministic separator leaves per
     worker.  A retained BTF already exposes SCC-level outer concurrency, so
     one leaf per worker preserves a coarser, lower-overhead forest inside its
     dominant block.  Both use the same cache-sized refinement window; this
     is a representation/resource policy, not an input classifier. */
  kls_generic_metis_tuning_ctx =
    trial->options.use_btf &&
    getenv("KLS_DISABLE_BTF_SINGLE_LEAF_ND") == NULL ? 7 : 5;
  kls_generic_metis_scale_hint = trial->scale_hint;
  memset(&kls_btf_stash, 0, sizeof(kls_btf_stash));
  trial->status = analyze_with_ordering(
    trial->n, trial->col_ptr, trial->row_idx, &trial->options,
    KLS_ORDERING_METIS, &trial->symbolic, &trial->common,
    &trial->separator);
  trial->seconds = kls_now_seconds() - start;
  kls_generic_metis_tuning_ctx = saved_tuning_ctx;
  kls_generic_metis_scale_hint = saved_scale_hint;
  kls_btf_stash = saved_btf_stash;
  return NULL;
}

static void kls_generic_nd_trial_join(kls_generic_nd_trial *trial) {
  if (trial != NULL && trial->active) {
    pthread_join(trial->thread, NULL);
    trial->active = 0;
  }
}

static void kls_generic_nd_trial_clear(kls_generic_nd_trial *trial) {
  if (trial == NULL) {
    return;
  }
  kls_generic_nd_trial_join(trial);
  if (trial->symbolic != NULL) {
    trilinos_klu_l_free_symbolic(&trial->symbolic, &trial->common);
  }
  kls_separator_analysis_clear(&trial->separator);
  memset(trial, 0, sizeof(*trial));
}

static int kls_generic_nd_trial_start(
  kls_generic_nd_trial *trial,
  kls_pattern_candidate *candidate,
  const kls_options *options,
  int use_btf) {
  if (trial == NULL || candidate == NULL || options == NULL ||
      trial->active || trial->symbolic != NULL ||
      candidate->col_ptr == NULL || candidate->row_idx == NULL ||
      options->ordering != KLS_ORDERING_AUTO ||
      !kls_repeated_update_workload(options) || options->threads < 2 ||
      candidate->n < 30000u || candidate->n > (UF_long)INT32_MAX ||
      candidate->nnz > (UF_long)INT32_MAX ||
      getenv("KLS_DISABLE_GENERIC_ND_PORTFOLIO") != NULL ||
      getenv("KLS_DISABLE_GENERIC_ND_ANALYZE_OVERLAP") != NULL) {
    return 0;
  }

  const int recurring_independent_btf =
    use_btf &&
    getenv("KLS_DISABLE_RECURRING_INDEPENDENT_BTF_ND") == NULL;
  /* For a recurring portfolio, a second independently computed maximum
     matching is itself a valid BTF representation.  It remains provisional:
     factor entry checks actual pivoted fill against the retained incumbent
     and solve-side residual validation polices the published numeric.  Short
     or symbolic-only workloads retain the stricter identical-coordinate
     certificate below. */
  if (use_btf && !recurring_independent_btf) {
    /* The concurrent analyses must retain the same BTF coordinate frame.
       A structural diagonal is an explicit perfect-matching certificate, so
       neither arm depends on a separately discovered transversal.  Without
       that certificate, keep the synchronous stash path: two individually
       valid maximum matchings can expose different separator schedules even
       when their fill is similar. */
    for (UF_long col = 0u; col < candidate->n; ++col) {
      int has_diagonal = 0;
      for (UF_long p = candidate->col_ptr[col];
           p < candidate->col_ptr[col + 1u]; ++p) {
        const UF_long row = candidate->row_idx[p];
        if (row == col) {
          has_diagonal = 1;
          break;
        }
      }
      if (!has_diagonal) {
        return 0;
      }
    }
  }

  /* Concurrency is useful once one canonical index stream no longer fits in
     the workers' aggregate cache allowance.  Below that point the retained
     minimum-degree analysis is already short and a second ordering can only
     add scheduling noise.  This is a byte/resource break-even bound, not a
     dimension or degree-profile classifier. */
  const double pattern_bytes =
    ((double)candidate->n + 1.0 + (double)candidate->nnz) *
    (double)sizeof(UF_long);
  double overlap_cache_allowance =
    4.0 * 1024.0 * 1024.0 * (double)options->threads;
  if (options->expected_refactorizations >= 16 &&
      getenv("KLS_DISABLE_RECURRING_ND_ANALYZE_OVERLAP") == NULL) {
    /* Across a long numeric lifecycle, overlap is useful once each worker
       receives one private-cache-sized slice of the immutable index stream.
       The ordinary aggregate-cache floor remains for short workloads, where
       a rejected ordering has no recurring payoff. */
    overlap_cache_allowance =
      1.0 * 1024.0 * 1024.0 * (double)options->threads;
  }
  if (!isfinite(pattern_bytes) ||
      pattern_bytes < overlap_cache_allowance) {
    return 0;
  }

  memset(trial, 0, sizeof(*trial));
  trial->n = candidate->n;
  trial->col_ptr = candidate->col_ptr;
  trial->row_idx = candidate->row_idx;
  trial->options = *options;
  trial->options.ordering = KLS_ORDERING_METIS;
  trial->options.use_btf = use_btf;
  trial->scale_hint = initial_scale(options);
  if (candidate->values != NULL) {
    trial->scale_hint = choose_auto_scale_from_pattern_impl(
      candidate->n, candidate->col_ptr, candidate->row_idx, options,
      candidate->values, NULL, NULL);
  }
  trial->status = KLS_ERR_ANALYZE_FAILED;
  if (pthread_create(&trial->thread, NULL,
                     kls_generic_nd_trial_main, trial) != 0) {
    memset(trial, 0, sizeof(*trial));
    return 0;
  }
  trial->active = 1;
  if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
    fprintf(stderr,
            "KLS generic ND: analyze overlap launched "
            "(bytes=%.3e workers=%d btf=%d)\n",
            pattern_bytes, options->threads, use_btf);
  }
  return 1;
}

/* Run nested dissection once, after AUTO has selected the matrix
   orientation, when the caller's repeated-update horizon and the retained
   minimum-degree estimate make another symbolic pass economically relevant.
   Requiring decisive Pareto improvements in both storage and arithmetic is
   intentionally stronger than the ordinary minimum-fill tournament: it lets
   us install the recurring representation before building any throwaway
   numeric, while declining marginal symbolic wins that pivoting or schedule
   quality could easily reverse.  The gates use only the requested lifecycle,
   resource budget, and realized symbolic economics. */
static void kls_maybe_promote_selected_generic_nd(
  kls_pattern_candidate *candidate,
  const kls_options *options,
  int bounded_symmetric_union,
  kls_generic_nd_trial *precomputed_no_btf_trial,
  kls_generic_nd_trial *precomputed_btf_trial) {
  /* The only absolute upper bound is the ordering backend's index width.
     Whether a large graph repays NodeND is decided below from retained fill,
     work, and the caller's numeric horizon; fixed million-row/eight-million-
     entry caps excluded exactly the high-work lifecycles that can amortize
     the trial. */
  if (candidate == NULL || options == NULL || candidate->symbolic == NULL ||
      options->ordering != KLS_ORDERING_AUTO ||
      candidate->selected_ordering == KLS_ORDERING_METIS ||
      candidate->block_order_perm != NULL ||
      !kls_repeated_update_workload(options) || options->threads < 2 ||
      (candidate->n < 30000u &&
       getenv("KLS_DISABLE_SUB_30K_GENERIC_ND") != NULL) ||
      candidate->n > (UF_long)INT32_MAX ||
      candidate->nnz > (UF_long)INT32_MAX ||
      getenv("KLS_DISABLE_GENERIC_ND_PORTFOLIO") != NULL) {
    return;
  }

  const double n = (double)candidate->n;
  const double active_n =
    candidate->symbolic->do_btf && candidate->symbolic->maxblock > 0u
      ? (double)candidate->symbolic->maxblock : n;
  double incumbent_fill = symbolic_score(candidate->symbolic);
  double incumbent_flops = candidate->symbolic->est_flops;
  int retained_btf_is_incumbent = 0;
  if (candidate->generic_btf_value_symbolic != NULL) {
    const double retained_btf_fill =
      symbolic_score(candidate->generic_btf_value_symbolic);
    if (isfinite(retained_btf_fill) && retained_btf_fill > 0.0 &&
        retained_btf_fill < incumbent_fill) {
      /* The representation tournament has retained this BTF arm for an
         actual value-time factor decision.  ND must beat that cheapest
         incumbent, not the no-BTF bootstrap that merely owns the candidate
         slot; otherwise a worse symbolic can preempt the better arm before
         its documented numeric/residual validation runs. */
      incumbent_fill = retained_btf_fill;
      incumbent_flops =
        candidate->generic_btf_value_symbolic->est_flops;
      retained_btf_is_incumbent = 1;
    }
  }
  const double provisional_refactor_horizon =
    1.0 + (double)options->expected_refactorizations;
  const double provisional_lifecycle_work =
    incumbent_flops * provisional_refactor_horizon;
  /* The original absolute work-density floor protects one-shot analyses,
     but it also excludes moderately dense factors whose declared recurring
     lifecycle can repay NodeND many times over.  Admit that second case only
     above half the ordinary density, at least sixteen updates, and a large
     aggregate work budget.  These are realized symbolic/workload economics,
     independent of matrix dimensions, names, or degree shape. */
  const int recurring_nd_work_economic =
    getenv("KLS_DISABLE_RELAXED_GENERIC_ND_WORK_GATE") == NULL &&
    options->expected_refactorizations >= 16 &&
    /* NodeND and a rejected numeric still traverse the complete problem,
       including scalar BTF fringe.  Charge the relaxed arm against all
       unknowns rather than only the largest active block; the established
       4096-density route below keeps its historical active-block meaning. */
    incumbent_flops >= 2048.0 * n &&
    isfinite(provisional_lifecycle_work) &&
    provisional_lifecycle_work >= 5.0e10;
  if (!isfinite(incumbent_fill) || incumbent_fill < 16.0 * active_n ||
      !(incumbent_flops >= 4096.0 * active_n ||
        recurring_nd_work_economic)) {
    return;
  }
  const int unbounded_numeric_provisional_economic =
    !bounded_symmetric_union &&
    options->expected_refactorizations >= 16 &&
    isfinite(provisional_lifecycle_work) &&
    (provisional_lifecycle_work >= 1.0e12 ||
     recurring_nd_work_economic);
  if (!bounded_symmetric_union &&
      !unbounded_numeric_provisional_economic) {
    return;
  }
  kls_options trial_options = *options;
  trial_options.ordering = KLS_ORDERING_METIS;
  /* Preserve the representation selected by the realized AMD fill/work
     comparison.  Each retained BTF block gets its own measured ND order,
     while a selected no-BTF arm receives the complete-graph candidate.  The
     numeric fill guards below still arbitrate the provisional ordering
     without an input-family classifier. */
  trial_options.use_btf = retained_btf_is_incumbent
    ? 1
    : (candidate->symbolic->do_btf ? 1 : 0);
  const int btf_single_leaf_nd = trial_options.use_btf &&
    getenv("KLS_DISABLE_BTF_SINGLE_LEAF_ND") == NULL;
  /* Separator construction must not depend on whether this same proposal
     happened to start before or after the minimum-degree representation
     verdict.  The early worker uses contexts 7/5 for BTF/no-BTF; using 8/6
     only in the synchronous fallback made launch order change retained fill
     (and recurring numeric work) even though every mathematical input and
     acceptance test was identical.  Keep the context representation-based
     and deterministic. */
  const int generic_nd_tuning_ctx = btf_single_leaf_nd ? 7 : 5;
  trilinos_klu_l_symbolic *trial_symbolic = NULL;
  trilinos_klu_l_common trial_common;
  kls_separator_analysis trial_separator;
  memset(&trial_separator, 0, sizeof(trial_separator));
  const int saved_tuning_ctx = kls_generic_metis_tuning_ctx;
  const int saved_scale_hint = kls_generic_metis_scale_hint;
  const struct kls_btf_stash_s saved_btf_stash = kls_btf_stash;
  int scale_hint = initial_scale(options);
  if (candidate->values != NULL) {
    scale_hint = choose_auto_scale_from_pattern_impl(
      candidate->n, candidate->col_ptr, candidate->row_idx, options,
      candidate->values, NULL, NULL);
  }
  int status = KLS_ERR_ANALYZE_FAILED;
  double trial_seconds = 0.0;
  int used_precomputed = 0;
  kls_generic_nd_trial *precomputed_trial = trial_options.use_btf
    ? precomputed_btf_trial : precomputed_no_btf_trial;
  if (precomputed_trial != NULL) {
    kls_generic_nd_trial_join(precomputed_trial);
    if (precomputed_trial->n == candidate->n &&
        precomputed_trial->col_ptr == candidate->col_ptr &&
        precomputed_trial->row_idx == candidate->row_idx &&
        precomputed_trial->options.use_btf == trial_options.use_btf &&
        precomputed_trial->status == KLS_OK &&
        precomputed_trial->symbolic != NULL) {
      trial_symbolic = precomputed_trial->symbolic;
      precomputed_trial->symbolic = NULL;
      trial_common = precomputed_trial->common;
      kls_separator_analysis_move(&trial_separator,
                                  &precomputed_trial->separator);
      trial_seconds = precomputed_trial->seconds;
      status = KLS_OK;
      used_precomputed = 1;
    }
  }
  if (!used_precomputed) {
    const double trial_start = kls_now_seconds();
    kls_generic_metis_tuning_ctx = generic_nd_tuning_ctx;
    kls_generic_metis_scale_hint = scale_hint;
    if (trial_options.use_btf && candidate->symbolic->do_btf &&
        candidate->symbolic->nblocks > 1u &&
        candidate->symbolic->P != NULL && candidate->symbolic->Q != NULL &&
        candidate->symbolic->R != NULL &&
        candidate->symbolic->structural_rank <= candidate->n) {
      /* A synchronous fallback can reuse the settled maximum transversal
         and SCC decomposition.  The early worker cannot: it intentionally
         runs beside the minimum-degree analysis that produces them. */
      kls_btf_stash.P = candidate->symbolic->P;
      kls_btf_stash.Q = candidate->symbolic->Q;
      kls_btf_stash.R = candidate->symbolic->R;
      kls_btf_stash.nblocks = candidate->symbolic->nblocks;
      kls_btf_stash.rank = candidate->symbolic->structural_rank;
      kls_btf_stash.n = candidate->n;
      kls_btf_stash.active = 1;
    }
    status = analyze_with_ordering(
      candidate->n, candidate->col_ptr, candidate->row_idx, &trial_options,
      KLS_ORDERING_METIS, &trial_symbolic, &trial_common, &trial_separator);
    kls_generic_metis_tuning_ctx = saved_tuning_ctx;
    kls_generic_metis_scale_hint = saved_scale_hint;
    kls_btf_stash = saved_btf_stash;
    trial_seconds = kls_now_seconds() - trial_start;
  }
  /* A diagnostic eager overlap can still use a coordinate different from the
     settled minimum-degree arm.  Do not discard that independently valid ND
     result merely because its representation differs: compare the realized
     symbolic streams directly and transfer a materially smaller one.  The
     default staged launch uses only the settled coordinate, so this alternate
     is normally absent.  Numeric fill and residual validation below remain
     authoritative in either case. */
  kls_generic_nd_trial *alternate_coordinate_trial =
    trial_options.use_btf
      ? precomputed_no_btf_trial : precomputed_btf_trial;
  if (status == KLS_OK && trial_symbolic != NULL &&
      alternate_coordinate_trial != NULL) {
    kls_generic_nd_trial_join(alternate_coordinate_trial);
    if (alternate_coordinate_trial->n == candidate->n &&
        alternate_coordinate_trial->col_ptr == candidate->col_ptr &&
        alternate_coordinate_trial->row_idx == candidate->row_idx &&
        alternate_coordinate_trial->status == KLS_OK &&
        alternate_coordinate_trial->symbolic != NULL) {
      const double retained_coordinate_fill =
        symbolic_score(trial_symbolic);
      const double alternate_coordinate_fill =
        symbolic_score(alternate_coordinate_trial->symbolic);
      const int alternate_coordinate_wins =
        isfinite(alternate_coordinate_fill) &&
        alternate_coordinate_fill > 0.0 &&
        retained_coordinate_fill < DBL_MAX / 4.0 &&
        isfinite(retained_coordinate_fill) &&
        retained_coordinate_fill > 0.0 &&
        alternate_coordinate_fill <=
          0.98 * retained_coordinate_fill;
      if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
        fprintf(stderr,
                "KLS generic ND coordinates: fill %.4e/%.4e "
                "btf=%d/%d -> %s\n",
                alternate_coordinate_fill, retained_coordinate_fill,
                alternate_coordinate_trial->symbolic->do_btf ? 1 : 0,
                trial_symbolic->do_btf ? 1 : 0,
                alternate_coordinate_wins ? "alternate" : "retained");
      }
      if (alternate_coordinate_wins) {
        trilinos_klu_l_free_symbolic(&trial_symbolic, &trial_common);
        kls_separator_analysis_clear(&trial_separator);
        trial_symbolic = alternate_coordinate_trial->symbolic;
        alternate_coordinate_trial->symbolic = NULL;
        trial_common = alternate_coordinate_trial->common;
        trial_options.use_btf = trial_symbolic->do_btf ? 1 : 0;
        kls_separator_analysis_move(
          &trial_separator, &alternate_coordinate_trial->separator);
      }
    }
  }
  if (status != KLS_OK || trial_symbolic == NULL) {
    kls_separator_analysis_clear(&trial_separator);
    return;
  }

  double trial_fill = symbolic_score(trial_symbolic);
  double trial_flops = trial_symbolic->est_flops;

  const double numeric_horizon =
    1.0 + (double)options->expected_refactorizations +
    (double)options->expected_solves;
  const int work_estimate_available =
    trial_flops > 0.0 && incumbent_flops > 0.0;
  /* A bounded union can produce an exact predicted fill certificate.  An
     unbounded union remains admissible only under the much larger lifecycle
     budget above; it bypasses prediction, builds one ordinary pivoted KLU
     numeric, and must beat the incumbent fill at factor entry. */
  const int exact_fill_provisional_candidate =
    bounded_symmetric_union &&
    trial_fill >= DBL_MAX / 4.0 &&
    trial_symbolic->structural_rank == candidate->n;
  const int unbounded_numeric_provisional_candidate =
    unbounded_numeric_provisional_economic &&
    trial_fill >= DBL_MAX / 4.0 &&
    trial_symbolic->structural_rank == candidate->n;
  const double incumbent_lifecycle_work =
    incumbent_flops * numeric_horizon;
  /* Below the portfolio's historically validated 30k-row scale, separator
     structure alone is not yet evidence of a recurring kernel win: on a
     compact factor, a near-tie tree can lose locality even when it removes a
     fraction of the entries.  Still admit these graphs when the realized
     symbolic stream removes at least five percent of retained storage.  The
     ordinary large-graph near-tie policy remains unchanged, where separator
     parallelism is itself a measured benefit, and the decisive 10% Pareto
     arm below continues to cover still stronger proposals. */
  const int sub_30k_nd_has_storage_margin =
    candidate->n >= 30000u ||
    (isfinite(trial_fill) && trial_fill > 0.0 &&
     trial_fill <= 0.95 * incumbent_fill);
  const int lifecycle_near_tie =
    isfinite(trial_fill) && trial_fill > 0.0 &&
    trial_fill <= 1.02 * incumbent_fill &&
    sub_30k_nd_has_storage_margin &&
    (trial_symbolic->structural_rank == candidate->n ||
     trial_symbolic->structural_rank == KLS_KLU_EMPTY) &&
    isfinite(incumbent_lifecycle_work) &&
    incumbent_lifecycle_work >= 2.0e10 &&
    (!work_estimate_available || trial_flops <= 1.05 * incumbent_flops);
  const int decisive_pareto_win =
    exact_fill_provisional_candidate ||
    unbounded_numeric_provisional_candidate || lifecycle_near_tie ||
    (isfinite(trial_fill) &&
     ((work_estimate_available &&
      trial_fill <= 0.90 * incumbent_fill &&
      trial_flops <= 0.80 * incumbent_flops) ||
     (!work_estimate_available &&
      trial_fill <= 0.90 * incumbent_fill)));
  if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
    fprintf(stderr,
            "KLS generic ND portfolio: fill %.4e/%.4e "
            "flops %.4e/%.4e "
            "analyze %.3fs -> %s\n",
            trial_fill, incumbent_fill, trial_flops, incumbent_flops,
            trial_seconds,
            lifecycle_near_tie ? "adopt-near-tie" :
              (decisive_pareto_win ? "adopt" : "reject"));
  }
  if (!decisive_pareto_win) {
    trilinos_klu_l_free_symbolic(&trial_symbolic, &trial_common);
    kls_separator_analysis_clear(&trial_separator);
    return;
  }

  const kls_ordering fallback_ordering = retained_btf_is_incumbent
    ? candidate->generic_btf_value_ordering
    : candidate->selected_ordering;
  const int fallback_use_btf = retained_btf_is_incumbent
    ? 1
    : (candidate->symbolic->do_btf ? 1 : 0);
  trilinos_klu_l_free_symbolic(&candidate->symbolic, &candidate->common);
  kls_separator_analysis_clear(&candidate->separator);
  candidate->symbolic = trial_symbolic;
  candidate->common = trial_common;
  candidate->selected_ordering = KLS_ORDERING_METIS;
  candidate->score = trial_fill;
  candidate->generic_nd_portfolio_selected = 1;
  candidate->generic_nd_bounded_symmetric_union = bounded_symmetric_union;
  candidate->generic_nd_lifecycle_near_tie = lifecycle_near_tie;
  candidate->generic_nd_fallback_ordering = fallback_ordering;
  candidate->generic_nd_fallback_use_btf = fallback_use_btf;
  candidate->generic_nd_fallback_fill = incumbent_fill;
  candidate->generic_nd_fallback_flops = incumbent_flops;
  kls_separator_analysis_move(&candidate->separator, &trial_separator);
}

static int should_start_auto_without_btf(UF_long n,
                                         const UF_long *col_ptr,
                                         const UF_long *row_idx,
                                         const kls_options *options,
                                         int large_spiked_metis_no_btf) {
  if (options == NULL || !options->use_btf) {
    return 0;
  }

  if (is_medium_low_degree_full_diagonal_pattern(n, col_ptr, row_idx)) {
    return 1;
  }
  if (is_medium_sparse_high_degree_diagonal_pattern(n, col_ptr, row_idx)) {
    return 1;
  }
  if (large_spiked_metis_no_btf) {
    return 1;
  }
  if (is_medium_dense_diagonal_high_degree_pattern(n, col_ptr, row_idx)) {
    return 1;
  }
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
  const double t0 = kls_now_seconds();
  (void)trilinos_klu_l_flops(solver->symbolic, solver->numeric,
                             &solver->common);
  const double t1 = kls_now_seconds();
  if (include_rcond) {
    (void)trilinos_klu_l_rcond(solver->symbolic, solver->numeric,
                               &solver->common);
  }
  if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
    fprintf(stderr, "KLS diag: flops %.3fs rcond %.3fs\n", t1 - t0,
            kls_now_seconds() - t1);
  }
}

static void kls_update_numeric_rcond(kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL) {
    return;
  }
  const double t0 = kls_now_seconds();
  (void)trilinos_klu_l_rcond(solver->symbolic, solver->numeric,
                             &solver->common);
  if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
    fprintf(stderr, "KLS diag: flops 0.000s rcond %.3fs\n",
            kls_now_seconds() - t0);
  }
}

/* The lifecycle guard needs only KLU's documented diagonal ratio.  Keep the
   scan branch-free so compilers can reduce min/max in vector lanes; the
   general diagnostic routine above retains KLU's full status semantics. */
static void kls_update_numeric_rcond_guard(kls_solver *solver) {
  if (solver == NULL || solver->numeric == NULL ||
      solver->numeric->Udiag == NULL || solver->n == 0u) {
    return;
  }
  const double *restrict udiag =
    (const double *)solver->numeric->Udiag;
  double umin = DBL_MAX;
  double umax = 0.0;
  int invalid = 0;
#pragma omp simd reduction(min:umin) reduction(max:umax) reduction(|:invalid)
  for (UF_long k = 0u; k < solver->n; ++k) {
    const double pivot = fabs(udiag[k]);
    invalid |= !isfinite(pivot) || pivot == 0.0;
    umin = pivot < umin ? pivot : umin;
    umax = pivot > umax ? pivot : umax;
  }
  solver->common.rcond = !invalid && umax > 0.0 ? umin / umax : 0.0;
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

static void kls_sort_match_entries_desc_full_radix(
  kls_match_entry *entries,
  UF_long count,
  kls_match_entry *scratch) {
  enum { KLS_MATCH_RADIX_MAX_BITS = 16,
         KLS_MATCH_RADIX_MAX_SIZE = 1 << KLS_MATCH_RADIX_MAX_BITS };
  size_t buckets[KLS_MATCH_RADIX_MAX_SIZE];
  const unsigned int radix_bits =
    getenv("KLS_MATCH_RADIX16") != NULL ||
        kls_medium_partial_static_metis_ctx
      ? 16u : 11u;
  kls_match_entry *src = entries;
  kls_match_entry *dst = scratch;
  uint64_t max_row_col = 0u;
  for (UF_long i = 0u; i < count; ++i) {
    max_row_col |= (uint64_t)entries[i].row;
    max_row_col |= (uint64_t)entries[i].col;
  }
  unsigned int index_bits = 1u;
  while (index_bits < 64u && (max_row_col >> index_bits) != 0u) {
    ++index_bits;
  }

  /* Stable LSD passes over the reverse integer keys implement the exact
     comparison order: descending weight, then descending row and column.
     Sparse-matrix indices need only their significant low passes; positive
     finite IEEE-754 weights retain numeric order in all 64 representation
     bits. */
  for (int key_kind = 0; key_kind < 3; ++key_kind) {
    const unsigned int key_bits = key_kind < 2 ? index_bits : 64u;
    for (unsigned int shift = 0u; shift < key_bits;
         shift += radix_bits) {
      const unsigned int pass_bits =
        key_bits - shift < radix_bits ? key_bits - shift : radix_bits;
      const unsigned int digit_count = 1u << pass_bits;
      memset(buckets, 0, (size_t)digit_count * sizeof(*buckets));
      const uint64_t mask = (UINT64_C(1) << pass_bits) - UINT64_C(1);
      for (UF_long i = 0u; i < count; ++i) {
        uint64_t raw = 0u;
        if (key_kind == 0) {
          raw = (uint64_t)src[i].col;
        } else if (key_kind == 1) {
          raw = (uint64_t)src[i].row;
        } else {
          memcpy(&raw, &src[i].weight, sizeof(raw));
        }
        buckets[(~raw >> shift) & mask]++;
      }
      size_t offset = 0u;
      for (unsigned int digit = 0u; digit < digit_count; ++digit) {
        const size_t count_for_digit = buckets[digit];
        buckets[digit] = offset;
        offset += count_for_digit;
      }
      for (UF_long i = 0u; i < count; ++i) {
        uint64_t raw = 0u;
        if (key_kind == 0) {
          raw = (uint64_t)src[i].col;
        } else if (key_kind == 1) {
          raw = (uint64_t)src[i].row;
        } else {
          memcpy(&raw, &src[i].weight, sizeof(raw));
        }
        const size_t digit = (size_t)((~raw >> shift) & mask);
        dst[buckets[digit]++] = src[i];
      }
      kls_match_entry *swap = src;
      src = dst;
      dst = swap;
    }
  }
  if (src != entries) {
    memcpy(entries, src, (size_t)count * sizeof(*entries));
  }
}

static void kls_sort_match_entries_desc(kls_match_entry *entries,
                                        UF_long count) {
  if (entries == NULL || count <= 1u) {
    return;
  }
  if (count < 256u || count > (UF_long)(SIZE_MAX / sizeof(*entries))) {
    qsort(entries, (size_t)count, sizeof(*entries),
          compare_match_entries_desc);
    return;
  }
  kls_match_entry *scratch =
    (kls_match_entry *)malloc((size_t)count * sizeof(*scratch));
  if (scratch == NULL) {
    qsort(entries, (size_t)count, sizeof(*entries),
          compare_match_entries_desc);
    return;
  }

  if (getenv("KLS_MATCH_FULL_RADIX") != NULL ||
      kls_medium_partial_static_metis_ctx) {
    kls_sort_match_entries_desc_full_radix(entries, count, scratch);
    free(scratch);
    return;
  }

  /* All accepted weights are positive finite doubles, whose IEEE-754 bit
     representation has the same order as the numeric value.  Radix-sort the
     high 32 key bits, which determine order except inside an equal-prefix
     run, then comparison-sort only those usually tiny runs.  This preserves
     the exact comparator order while halving the full-array traffic of an
     eight-byte radix key.  Complementing the key produces descending order. */
  enum { KLS_MATCH_HIGH_RADIX_MAX_BITS = 16,
         KLS_MATCH_HIGH_RADIX_MAX_SIZE =
           1 << KLS_MATCH_HIGH_RADIX_MAX_BITS };
  size_t buckets[KLS_MATCH_HIGH_RADIX_MAX_SIZE];
  const unsigned int radix_bits =
    getenv("KLS_MATCH_HIGH_RADIX16") != NULL ||
        kls_low_work_symmetric_partial_diagonal_pts_ctx
      ? 16u : 8u;
  const unsigned int digit_count = 1u << radix_bits;
  const uint64_t mask = (UINT64_C(1) << radix_bits) - UINT64_C(1);
  kls_match_entry *src = entries;
  kls_match_entry *dst = scratch;
  for (unsigned int shift = 32u; shift < 64u; shift += radix_bits) {
    memset(buckets, 0, (size_t)digit_count * sizeof(*buckets));
    for (UF_long i = 0u; i < count; ++i) {
      uint64_t bits = 0u;
      memcpy(&bits, &src[i].weight, sizeof(bits));
      const unsigned int digit =
        (unsigned int)((~bits >> shift) & mask);
      buckets[digit]++;
    }
    size_t offset = 0u;
    for (unsigned int digit = 0u; digit < digit_count; ++digit) {
      const size_t count_for_digit = buckets[digit];
      buckets[digit] = offset;
      offset += count_for_digit;
    }
    for (UF_long i = 0u; i < count; ++i) {
      uint64_t bits = 0u;
      memcpy(&bits, &src[i].weight, sizeof(bits));
      const unsigned int digit =
        (unsigned int)((~bits >> shift) & mask);
      dst[buckets[digit]++] = src[i];
    }
    kls_match_entry *swap = src;
    src = dst;
    dst = swap;
  }
  if (src != entries) {
    memcpy(entries, src, (size_t)count * sizeof(*entries));
  }
  free(scratch);

  if (getenv("KLS_MATCH_SKIP_LOW_SORT") != NULL ||
      kls_low_work_symmetric_partial_diagonal_pts_ctx) {
    return;
  }

  /* Complete each equal-prefix run, including the old comparator's
     deterministic row/column tie break for exactly equal weights. */
  for (UF_long begin = 0u; begin < count;) {
    uint64_t begin_bits = 0u;
    memcpy(&begin_bits, &entries[begin].weight, sizeof(begin_bits));
    UF_long end = begin + 1u;
    while (end < count) {
      uint64_t end_bits = 0u;
      memcpy(&end_bits, &entries[end].weight, sizeof(end_bits));
      if ((end_bits >> 32u) != (begin_bits >> 32u)) {
        break;
      }
      ++end;
    }
    if (end - begin > 1u) {
      qsort(entries + begin, (size_t)(end - begin), sizeof(*entries),
            compare_match_entries_desc);
    }
    begin = end;
  }
}

static int compare_compact_match_entries_desc(const void *a, const void *b) {
  const kls_compact_match_entry *left =
    (const kls_compact_match_entry *)a;
  const kls_compact_match_entry *right =
    (const kls_compact_match_entry *)b;
  if (left->weight != right->weight) {
    return left->weight > right->weight ? -1 : 1;
  }
  if (left->row_col != right->row_col) {
    return left->row_col > right->row_col ? -1 : 1;
  }
  return 0;
}

static void kls_sort_compact_match_entries_desc(
  kls_compact_match_entry *entries,
  UF_long count) {
  if (entries == NULL || count <= 1u) {
    return;
  }
  if (count < 256u || count > (UF_long)(SIZE_MAX / sizeof(*entries))) {
    qsort(entries, (size_t)count, sizeof(*entries),
          compare_compact_match_entries_desc);
    return;
  }
  kls_compact_match_entry *scratch = (kls_compact_match_entry *)malloc(
    (size_t)count * sizeof(*scratch));
  if (scratch == NULL) {
    qsort(entries, (size_t)count, sizeof(*entries),
          compare_compact_match_entries_desc);
    return;
  }

  kls_compact_match_entry *src = entries;
  kls_compact_match_entry *dst = scratch;
  for (unsigned int pass = 0u; pass < 4u; ++pass) {
    size_t buckets[256] = {0};
    size_t offsets[256];
    const unsigned int shift = 32u + 8u * pass;
    for (UF_long i = 0u; i < count; ++i) {
      uint64_t bits = 0u;
      memcpy(&bits, &src[i].weight, sizeof(bits));
      const unsigned int digit =
        (unsigned int)((~bits >> shift) & UINT64_C(0xff));
      buckets[digit]++;
    }
    offsets[0] = 0u;
    for (unsigned int digit = 1u; digit < 256u; ++digit) {
      offsets[digit] = offsets[digit - 1u] + buckets[digit - 1u];
    }
    for (UF_long i = 0u; i < count; ++i) {
      uint64_t bits = 0u;
      memcpy(&bits, &src[i].weight, sizeof(bits));
      const unsigned int digit =
        (unsigned int)((~bits >> shift) & UINT64_C(0xff));
      dst[offsets[digit]++] = src[i];
    }
    kls_compact_match_entry *swap = src;
    src = dst;
    dst = swap;
  }
  if (src != entries) {
    memcpy(entries, src, (size_t)count * sizeof(*entries));
  }
  free(scratch);

  for (UF_long begin = 0u; begin < count;) {
    uint64_t begin_bits = 0u;
    memcpy(&begin_bits, &entries[begin].weight, sizeof(begin_bits));
    UF_long end = begin + 1u;
    while (end < count) {
      uint64_t end_bits = 0u;
      memcpy(&end_bits, &entries[end].weight, sizeof(end_bits));
      if ((end_bits >> 32u) != (begin_bits >> 32u)) {
        break;
      }
      ++end;
    }
    if (end - begin > 1u) {
      qsort(entries + begin, (size_t)(end - begin), sizeof(*entries),
            compare_compact_match_entries_desc);
    }
    begin = end;
  }
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
                                 int build_log_weights,
                                 kls_row_match_graph *graph) {
  if (n <= 0 || entries == NULL || graph == NULL) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  memset(graph, 0, sizeof(*graph));

  UF_long *row_ptr = (UF_long *)calloc((size_t)n + 1u, sizeof(*row_ptr));
  UF_long *col_idx = (UF_long *)malloc((size_t)entry_count * sizeof(*col_idx));
  double *log_weight = build_log_weights
    ? (double *)malloc((size_t)entry_count * sizeof(*log_weight)) : NULL;
  const int raw_swap_weights =
    build_log_weights && n > 50000u &&
    (getenv("KLS_MATCH_RAW_SWAP_WEIGHTS") != NULL ||
     kls_medium_partial_static_metis_ctx);
  UF_long *next = (UF_long *)malloc((size_t)n * sizeof(*next));
  if (row_ptr == NULL || col_idx == NULL ||
      (build_log_weights && log_weight == NULL) || next == NULL) {
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
    if (build_log_weights) {
      log_weight[dst] = raw_swap_weights
        ? entries[i].weight : log(entries[i].weight);
    }
  }

  free(next);
  graph->row_ptr = row_ptr;
  graph->col_idx = col_idx;
  graph->log_weight = log_weight;
  graph->weights_are_raw = raw_swap_weights;
  return KLS_OK;
}

static double row_match_stored_weight(const kls_row_match_graph *graph,
                                      UF_long row,
                                      UF_long col) {
  for (UF_long p = graph->row_ptr[row]; p < graph->row_ptr[row + 1u]; ++p) {
    if (graph->col_idx[p] == col) {
      return graph->log_weight[p];
    }
  }
  return -DBL_MAX;
}

static double row_match_log_weight(const kls_row_match_graph *graph,
                                   UF_long row,
                                   UF_long col) {
  const double stored = row_match_stored_weight(graph, row, col);
  return graph->weights_are_raw && stored != -DBL_MAX ? log(stored) : stored;
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

static int build_exact_assignment_dual_scaling(UF_long n,
                                               const UF_long *row_perm,
                                               const double *col_max_log,
                                               const double *row_potential,
                                               const double *col_potential,
                                               double **row_scale_out,
                                               double **col_scale_out) {
  if (n <= 0 || row_perm == NULL || col_max_log == NULL ||
      row_potential == NULL || col_potential == NULL ||
      row_scale_out == NULL || col_scale_out == NULL) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  *row_scale_out = NULL;
  *col_scale_out = NULL;

  double *row_scale = (double *)malloc((size_t)n * sizeof(*row_scale));
  double *col_scale = (double *)malloc((size_t)n * sizeof(*col_scale));
  if (row_scale == NULL || col_scale == NULL) {
    free(row_scale);
    free(col_scale);
    return KLS_ERR_OUT_OF_MEMORY;
  }

  double row_log_sum = 0.0;
  double col_log_sum = 0.0;
  for (UF_long i = 0; i < n; ++i) {
    if (row_perm[i] >= n || col_max_log[i] == -DBL_MAX ||
        !isfinite(row_potential[i]) || !isfinite(col_potential[i]) ||
        !isfinite(col_max_log[i])) {
      free(row_scale);
      free(col_scale);
      return KLS_ERR_UNSUPPORTED;
    }
    row_log_sum += -row_potential[i];
    col_log_sum += -col_max_log[i] + col_potential[i];
  }

  const double shift = 0.5 * (col_log_sum - row_log_sum) / (double)n;
  const double max_log_scale = log(1.0e12);
  for (UF_long row = 0; row < n; ++row) {
    const UF_long scaled_row = row_perm[row];
    const double row_log_scale = -row_potential[row] + shift;
    if (scaled_row >= n || !isfinite(row_log_scale) ||
        fabs(row_log_scale) > max_log_scale) {
      free(row_scale);
      free(col_scale);
      return KLS_ERR_UNSUPPORTED;
    }
    row_scale[scaled_row] = clamp_matching_scale(exp(row_log_scale));
  }
  for (UF_long col = 0; col < n; ++col) {
    const double col_log_scale =
      -col_max_log[col] + col_potential[col] - shift;
    if (!isfinite(col_log_scale) ||
        fabs(col_log_scale) > max_log_scale) {
      free(row_scale);
      free(col_scale);
      return KLS_ERR_UNSUPPORTED;
    }
    col_scale[col] = clamp_matching_scale(exp(col_log_scale));
  }

  *row_scale_out = row_scale;
  *col_scale_out = col_scale;
  return KLS_OK;
}

static int build_exact_numeric_row_match(UF_long n,
                                         UF_long nnz,
                                         const UF_long *col_ptr,
                                         const UF_long *row_idx,
                                         const double *numeric_values,
                                         UF_long **row_perm_out,
                                         UF_long **col_match_out,
                                         UF_long *matched_out,
                                         double **row_scale_out,
                                         double **col_scale_out) {
  if (n <= 0 || col_ptr == NULL || row_idx == NULL || numeric_values == NULL ||
      row_perm_out == NULL || col_match_out == NULL || matched_out == NULL ||
      row_scale_out == NULL || col_scale_out == NULL) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  *row_perm_out = NULL;
  *col_match_out = NULL;
  *matched_out = 0;
  *row_scale_out = NULL;
  *col_scale_out = NULL;

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
  int status = build_row_match_graph(n, valid_edges, entries, 1, &graph);
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
#define KLS_FREE_EXACT_MATCH_WORK() do { \
    free(col_max_log);                     \
    free(row_perm);                        \
    free(col_match);                       \
    free(row_potential);                   \
    free(col_potential);                   \
    free(dist_row);                        \
    free(dist_col);                        \
    free(matched_cost_by_col);             \
    free(prev_cost_for_col);               \
    free(prev_row_for_col);                \
    free(seen_row);                        \
    free(seen_col);                        \
  } while (0)
  kls_min_heap heap;
  memset(&heap, 0, sizeof(heap));
  if (row_potential == NULL || col_potential == NULL ||
      dist_row == NULL || dist_col == NULL ||
      matched_cost_by_col == NULL || prev_cost_for_col == NULL ||
      prev_row_for_col == NULL || seen_row == NULL || seen_col == NULL) {
    free_row_match_graph(&graph);
    KLS_FREE_EXACT_MATCH_WORK();
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
      KLS_FREE_EXACT_MATCH_WORK();
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
      KLS_FREE_EXACT_MATCH_WORK();
      return KLS_ERR_INVALID_ARGUMENT;
    }
    matched++;
  }

  if (matched < n) {
    kls_heap_free(&heap);
    free_row_match_graph(&graph);
    KLS_FREE_EXACT_MATCH_WORK();
    return KLS_ERR_UNSUPPORTED;
  }

  double *row_scale = NULL;
  double *col_scale = NULL;
  const int scale_status =
    build_exact_assignment_dual_scaling(n, row_perm, col_max_log,
                                        row_potential, col_potential,
                                        &row_scale, &col_scale);
  if (scale_status == KLS_ERR_OUT_OF_MEMORY) {
    kls_heap_free(&heap);
    free_row_match_graph(&graph);
    KLS_FREE_EXACT_MATCH_WORK();
    return scale_status;
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
  *matched_out = matched;
  *row_perm_out = row_perm;
  *col_match_out = col_match;
  if (scale_status == KLS_OK) {
    *row_scale_out = row_scale;
    *col_scale_out = col_scale;
  }
#undef KLS_FREE_EXACT_MATCH_WORK
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

#ifdef KLS_HAVE_SPRAL_SCALING
static int augment_numeric_row_match_from_csc(UF_long n,
                                              UF_long nnz,
                                              const UF_long *col_ptr,
                                              const UF_long *row_idx,
                                              const double *numeric_values,
                                              UF_long *row_perm,
                                              UF_long *col_match,
                                              UF_long *matched_io) {
  if (n <= 0 || col_ptr == NULL || row_idx == NULL ||
      numeric_values == NULL || row_perm == NULL || col_match == NULL ||
      matched_io == NULL) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  if (*matched_io >= n) {
    return KLS_OK;
  }

  kls_match_entry *entries =
    (kls_match_entry *)malloc((size_t)nnz * sizeof(*entries));
  if (entries == NULL) {
    return KLS_ERR_OUT_OF_MEMORY;
  }

  UF_long entry_count = 0;
  for (UF_long col = 0; col < n; ++col) {
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      const UF_long row = row_idx[p];
      const double weight = fabs(numeric_values[p]);
      if (row < n && weight > 0.0 && isfinite(weight)) {
        entries[entry_count].weight = weight;
        entries[entry_count].row = row;
        entries[entry_count].col = col;
        entry_count++;
      }
    }
  }
  if (entry_count == 0) {
    free(entries);
    return KLS_ERR_UNSUPPORTED;
  }

  kls_row_match_graph graph;
  const int status =
    build_row_match_graph(n, entry_count, entries, 0, &graph);
  free(entries);
  if (status != KLS_OK) {
    return status;
  }

  *matched_io = augment_numeric_row_match_layered(n, &graph, row_perm,
                                                  col_match);
  free_row_match_graph(&graph);
  return KLS_OK;
}

static int incomplete_spral_match_augment_is_worthwhile(UF_long n,
                                                        UF_long matched) {
  if (matched >= n) {
    return 0;
  }
  const UF_long unmatched = n - matched;
  if (unmatched <= 4096u) {
    return 1;
  }
  return 1000u * matched >= 990u * n && unmatched <= 20000u;
}
#endif

static void improve_numeric_row_match_by_swaps(UF_long n,
                                               const kls_row_match_graph *graph,
                                               UF_long *row_perm,
                                               UF_long *col_match) {
  /* Rows are weight-sorted descending when the graph is built, so the
     truncated candidate window already holds each row's strongest edges;
     eight of them find the same improvement cycles at half the search
     breadth. */
  const UF_long max_candidates_per_row = 8u;
  int max_rounds = kls_medium_partial_static_metis_ctx ? 1 : 5;
  {
    const char *rounds_env = getenv("KLS_MATCH_SWAP_ROUNDS");
    if (rounds_env != NULL && rounds_env[0] != '\0') {
      const int parsed = atoi(rounds_env);
      if (parsed >= 1 && parsed <= 5) {
        max_rounds = parsed;
      }
    }
  }
  double *current_log_weight =
    (double *)malloc((size_t)n * sizeof(*current_log_weight));
  if (current_log_weight == NULL) {
    return;
  }
  for (UF_long row = 0; row < n; ++row) {
    const UF_long col = row_perm[row];
    current_log_weight[row] = col == KLS_KLU_EMPTY
      ? -DBL_MAX : row_match_stored_weight(graph, row, col);
  }
  for (int round = 0; round < max_rounds; ++round) {
    UF_long changes = 0;
    for (UF_long row = 0; row < n; ++row) {
      const UF_long current_col = row_perm[row];
      if (current_col == KLS_KLU_EMPTY) {
        continue;
      }
      const double current_weight = current_log_weight[row];
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
          row_match_stored_weight(graph, other_row, current_col);
        if (other_to_current == -DBL_MAX) {
          continue;
        }
        const double other_weight = current_log_weight[other_row];
        const int improves = graph->weights_are_raw
          ? (long double)graph->log_weight[p] *
              (long double)other_to_current >
            (long double)current_weight * (long double)other_weight *
              (1.0L + 1.0e-12L)
          : graph->log_weight[p] + other_to_current - current_weight -
              other_weight > 1.0e-12;
        if (improves) {
          row_perm[row] = trial_col;
          row_perm[other_row] = current_col;
          col_match[trial_col] = row;
          col_match[current_col] = other_row;
          current_log_weight[row] = graph->log_weight[p];
          current_log_weight[other_row] = other_to_current;
          changes++;
          break;
        }
      }
    }
    if (getenv("KLS_TRACE_MATCH_SWAPS") != NULL) {
      fprintf(stderr, "KLS match swaps: round=%d changes=%ld\n",
              round + 1, (long)changes);
    }
    if (changes == 0) {
      break;
    }
  }
  free(current_log_weight);
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
  const UF_long spral_matched = matched;

  if (incomplete_spral_match_augment_is_worthwhile(n, matched)) {
    const int augment_status =
      augment_numeric_row_match_from_csc(n, nnz, col_ptr, row_idx,
                                         numeric_values, row_perm,
                                         col_match, &matched);
    if (augment_status == KLS_ERR_OUT_OF_MEMORY) {
      free(row_scale);
      free(col_scale);
      free(row_perm);
      free(col_match);
      return augment_status;
    }
  }

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

  int scaling_usable = inform.flag == 0 && spral_matched == n;
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
  const UF_long spral_matched = matched;

  if (incomplete_spral_match_augment_is_worthwhile(n, matched)) {
    const int augment_status =
      augment_numeric_row_match_from_csc(n, nnz, col_ptr, row_idx,
                                         numeric_values, row_perm,
                                         col_match, &matched);
    if (augment_status == KLS_ERR_OUT_OF_MEMORY) {
      free(row_scale);
      free(col_scale);
      free(row_perm);
      free(col_match);
      return augment_status;
    }
  }

  const int complete_status =
    complete_unmatched_row_match(n, row_perm, col_match);
  free(col_match);
  if (complete_status != KLS_OK) {
    free(row_scale);
    free(col_scale);
    free(row_perm);
    return complete_status;
  }

  int scaling_usable = inform.flag == 0 && spral_matched == n;
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

static int build_compact_row_match_graph(
  UF_long n,
  UF_long entry_count,
  const kls_compact_match_entry *entries,
  kls_row_match_graph *graph) {
  if (entries == NULL || graph == NULL) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  memset(graph, 0, sizeof(*graph));
  UF_long *row_ptr = (UF_long *)calloc((size_t)n + 1u,
                                       sizeof(*row_ptr));
  UF_long *col_idx = (UF_long *)malloc(
    (size_t)entry_count * sizeof(*col_idx));
  UF_long *next = (UF_long *)malloc((size_t)n * sizeof(*next));
  if (row_ptr == NULL || col_idx == NULL || next == NULL) {
    free(row_ptr);
    free(col_idx);
    free(next);
    return KLS_ERR_OUT_OF_MEMORY;
  }
  for (UF_long i = 0u; i < entry_count; ++i) {
    const UF_long row = (UF_long)(entries[i].row_col >> 16u);
    row_ptr[row + 1u]++;
  }
  for (UF_long row = 0u; row < n; ++row) {
    row_ptr[row + 1u] += row_ptr[row];
  }
  memcpy(next, row_ptr, (size_t)n * sizeof(*next));
  for (UF_long i = 0u; i < entry_count; ++i) {
    const UF_long row = (UF_long)(entries[i].row_col >> 16u);
    col_idx[next[row]++] =
      (UF_long)(entries[i].row_col & UINT32_C(0xffff));
  }
  free(next);
  graph->row_ptr = row_ptr;
  graph->col_idx = col_idx;
  graph->log_weight = NULL;
  return KLS_OK;
}

static int build_compact_numeric_row_match(
  UF_long n,
  UF_long nnz,
  const UF_long *col_ptr,
  const UF_long *row_idx,
  const double *numeric_values,
  UF_long **row_perm_out,
  UF_long *matched_out) {
  const int trace_phases = getenv("KLS_TRACE_MATCH_PHASES") != NULL;
  const double trace_start = trace_phases ? kls_now_seconds() : 0.0;
  if (n == 0u || n > (UF_long)UINT16_MAX + 1u || col_ptr == NULL ||
      row_idx == NULL || numeric_values == NULL || row_perm_out == NULL ||
      matched_out == NULL) {
    return KLS_ERR_INVALID_ARGUMENT;
  }

  kls_compact_match_entry *entries = (kls_compact_match_entry *)malloc(
    (size_t)nnz * sizeof(*entries));
  UF_long *row_perm = (UF_long *)malloc((size_t)n * sizeof(*row_perm));
  UF_long *col_match = (UF_long *)malloc((size_t)n * sizeof(*col_match));
  if (entries == NULL || row_perm == NULL || col_match == NULL) {
    free(entries);
    free(row_perm);
    free(col_match);
    return KLS_ERR_OUT_OF_MEMORY;
  }

  UF_long entry_count = 0u;
  for (UF_long col = 0u; col < n; ++col) {
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      const UF_long row = row_idx[p];
      const double weight = fabs(numeric_values[p]);
      if (row < n && weight > 0.0 && isfinite(weight)) {
        entries[entry_count].weight = weight;
        entries[entry_count].row_col =
          ((uint32_t)row << 16u) | (uint32_t)col;
        entry_count++;
      }
    }
  }
  const double trace_filled = trace_phases ? kls_now_seconds() : 0.0;
  kls_sort_compact_match_entries_desc(entries, entry_count);
  const double trace_sorted = trace_phases ? kls_now_seconds() : 0.0;

  for (UF_long i = 0u; i < n; ++i) {
    row_perm[i] = KLS_KLU_EMPTY;
    col_match[i] = KLS_KLU_EMPTY;
  }
  UF_long matched = 0u;
  for (UF_long i = 0u; i < entry_count && matched < n; ++i) {
    const UF_long row = (UF_long)(entries[i].row_col >> 16u);
    const UF_long col =
      (UF_long)(entries[i].row_col & UINT32_C(0xffff));
    if (row_perm[row] == KLS_KLU_EMPTY &&
        col_match[col] == KLS_KLU_EMPTY) {
      row_perm[row] = col;
      col_match[col] = row;
      matched++;
    }
  }
  const UF_long greedy_matched = matched;
  const double trace_greedy = trace_phases ? kls_now_seconds() : 0.0;

  kls_row_match_graph graph;
  const int graph_status = build_compact_row_match_graph(
    n, entry_count, entries, &graph);
  if (graph_status == KLS_OK) {
    matched = augment_numeric_row_match(n, &graph, row_perm, col_match);
    free_row_match_graph(&graph);
  } else if (graph_status == KLS_ERR_OUT_OF_MEMORY) {
    free(entries);
    free(row_perm);
    free(col_match);
    return graph_status;
  }

  const int complete_status =
    complete_unmatched_row_match(n, row_perm, col_match);
  free(entries);
  free(col_match);
  if (complete_status != KLS_OK) {
    free(row_perm);
    return complete_status;
  }
  if (trace_phases) {
    const double trace_done = kls_now_seconds();
    fprintf(stderr,
            "KLS compact match phases: fill=%.3f sort=%.3f greedy=%.3f "
            "improve=%.3f matched=%ld->%ld total=%.3f ms\n",
            (trace_filled - trace_start) * 1e3,
            (trace_sorted - trace_filled) * 1e3,
            (trace_greedy - trace_sorted) * 1e3,
            (trace_done - trace_greedy) * 1e3,
            (long)greedy_matched, (long)matched,
            (trace_done - trace_start) * 1e3);
  }
  *row_perm_out = row_perm;
  *matched_out = matched;
#undef KLS_FREE_GREEDY_MATCH
  return KLS_OK;
}

static int build_greedy_numeric_row_match(UF_long n,
                                          UF_long nnz,
                                          const UF_long *col_ptr,
                                          const UF_long *row_idx,
                                          const double *numeric_values,
                                          int compact_missing_diagonal_match,
                                          int improve_matching,
                                          UF_long **row_perm_out,
                                          UF_long *matched_out,
                                          int *exact_matching_out,
                                          int *exact_matching_scaling_out,
                                          int *spral_matching_out,
                                          double **row_scale_out,
                                          double **col_scale_out) {
  const int trace_phases = getenv("KLS_TRACE_MATCH_PHASES") != NULL;
  const double trace_start = trace_phases ? kls_now_seconds() : 0.0;
  if (n <= 0 || col_ptr == NULL || row_idx == NULL || numeric_values == NULL ||
      row_perm_out == NULL || matched_out == NULL ||
      exact_matching_out == NULL || exact_matching_scaling_out == NULL ||
      spral_matching_out == NULL ||
      row_scale_out == NULL || col_scale_out == NULL) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  *row_perm_out = NULL;
  *matched_out = 0;
  *exact_matching_out = 0;
  *exact_matching_scaling_out = 0;
  *spral_matching_out = 0;
  *row_scale_out = NULL;
  *col_scale_out = NULL;

  if (improve_matching && compact_missing_diagonal_match &&
      getenv("KLS_DISABLE_COMPACT_MISSING_DIAGONAL_MATCH") == NULL &&
      getenv("KLS_DISABLE_COMPACT_MISSING_DIAGONAL_MATCHER") == NULL) {
    return build_compact_numeric_row_match(
      n, nnz, col_ptr, row_idx, numeric_values, row_perm_out, matched_out);
  }

  if (improve_matching && exact_weighted_match_is_affordable(n, nnz)) {
    UF_long *exact_row_perm = NULL;
    UF_long *exact_col_match = NULL;
    UF_long exact_matched = 0;
    double *exact_row_scale = NULL;
    double *exact_col_scale = NULL;
    const int exact_status =
      build_exact_numeric_row_match(n, nnz, col_ptr, row_idx, numeric_values,
                                    &exact_row_perm, &exact_col_match,
                                    &exact_matched, &exact_row_scale,
                                    &exact_col_scale);
    if (exact_status == KLS_OK) {
      const int complete_status =
        complete_unmatched_row_match(n, exact_row_perm, exact_col_match);
      free(exact_col_match);
      if (complete_status != KLS_OK) {
        free(exact_row_perm);
        free(exact_row_scale);
        free(exact_col_scale);
        return complete_status;
      }
      *row_perm_out = exact_row_perm;
      *matched_out = exact_matched;
      *exact_matching_out = 1;
      if (exact_row_scale != NULL && exact_col_scale != NULL) {
        *row_scale_out = exact_row_scale;
        *col_scale_out = exact_col_scale;
        *exact_matching_scaling_out = 1;
      } else {
        free(exact_row_scale);
        free(exact_col_scale);
      }
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
#define KLS_FREE_GREEDY_MATCH() do { \
    free(entries);                    \
    free(row_perm);                   \
    free(col_match);                  \
    free(row_used);                   \
    free(col_used);                   \
  } while (0)
  if (entries == NULL || row_perm == NULL || col_match == NULL ||
      row_used == NULL || col_used == NULL) {
    KLS_FREE_GREEDY_MATCH();
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
  const double trace_filled = trace_phases ? kls_now_seconds() : 0.0;
  kls_sort_match_entries_desc(entries, entry_count);
  const double trace_sorted = trace_phases ? kls_now_seconds() : 0.0;

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
  const UF_long greedy_matched = matched;
  const double trace_greedy = trace_phases ? kls_now_seconds() : 0.0;

  if (improve_matching) {
    kls_row_match_graph graph;
    const int improve_weights =
      !compact_missing_diagonal_match &&
      !kls_low_work_symmetric_partial_diagonal_pts_ctx &&
      getenv("KLS_MATCH_SKIP_SWAPS") == NULL;
    int status = build_row_match_graph(n, entry_count, entries,
                                       improve_weights, &graph);
    if (status == KLS_OK) {
      if (n >= 50000u && n - matched >= 1024u) {
        matched = augment_numeric_row_match_layered(n, &graph, row_perm,
                                                    col_match);
      } else {
        matched = augment_numeric_row_match(n, &graph, row_perm, col_match);
      }
      if (matched == n && improve_weights) {
        if (!kls_low_work_symmetric_partial_diagonal_pts_ctx &&
            getenv("KLS_MATCH_SKIP_PAIR_SWAPS") == NULL) {
          improve_numeric_row_match_by_swaps(n, &graph, row_perm, col_match);
        }
        if (n > 20000u && getenv("KLS_MATCH_SKIP_PATHS") == NULL &&
            getenv("KLS_ENABLE_GENERIC_MATCH_PATHS") != NULL) {
          /* The cycle-path search costs a quarter of a small matrix's
             whole initialization; the downstream acceptance comparison
             guards match quality either way (a weaker match that loses
             just keeps the baseline factor).  Generic AUTO therefore keeps
             the inexpensive cardinality augmentation and pair swaps.  This
             explicit diagnostic opt-in exercises the full path pass; actual
             factor and solve checks remain the acceptance authority. */
          improve_numeric_row_match_by_paths(n, &graph, row_perm, col_match);
        }
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
                KLS_FREE_GREEDY_MATCH();
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
            *exact_matching_scaling_out = 0;
            *spral_matching_out = 1;
          }
          free(spral_row_perm);
          free(spral_row_scale);
          free(spral_col_scale);
        } else if (spral_status == KLS_ERR_OUT_OF_MEMORY) {
          KLS_FREE_GREEDY_MATCH();
          return spral_status;
        }
      }
#endif
      free_row_match_graph(&graph);
    } else if (status == KLS_ERR_OUT_OF_MEMORY) {
      KLS_FREE_GREEDY_MATCH();
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
        *exact_matching_scaling_out = 0;
        *spral_matching_out = 1;
        for (UF_long col = 0; col < n; ++col) {
          col_match[col] = KLS_KLU_EMPTY;
        }
        for (UF_long row = 0; row < n; ++row) {
          const UF_long col = row_perm[row];
          if (col >= n || col_match[col] != KLS_KLU_EMPTY) {
            KLS_FREE_GREEDY_MATCH();
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
      KLS_FREE_GREEDY_MATCH();
      return spral_status;
    }
  }
#endif

  const int complete_status =
    complete_unmatched_row_match(n, row_perm, col_match);
  if (complete_status != KLS_OK) {
    KLS_FREE_GREEDY_MATCH();
    return complete_status;
  }

  free(entries);
  free(col_match);
  free(row_used);
  free(col_used);
  if (trace_phases) {
    const double trace_done = kls_now_seconds();
    fprintf(stderr,
            "KLS match phases: fill=%.3f sort=%.3f greedy=%.3f "
            "improve=%.3f matched=%ld->%ld total=%.3f ms\n",
            (trace_filled - trace_start) * 1e3,
            (trace_sorted - trace_filled) * 1e3,
            (trace_greedy - trace_sorted) * 1e3,
            (trace_done - trace_greedy) * 1e3,
            (long)greedy_matched, (long)matched,
            (trace_done - trace_start) * 1e3);
  }
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

  if (base_col_ptr[n] != nnz) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  UF_long max_width = 0u;
  for (UF_long col = 0; col < n; ++col) {
    const UF_long begin = base_col_ptr[col];
    const UF_long end = base_col_ptr[col + 1u];
    if (begin > end || end > nnz) {
      return KLS_ERR_INVALID_ARGUMENT;
    }
    const UF_long width = end - begin;
    if (width > max_width) {
      max_width = width;
    }
  }

  UF_long *col_ptr =
    (UF_long *)malloc(((size_t)n + 1u) * sizeof(*col_ptr));
  UF_long *row_idx = (UF_long *)malloc((size_t)nnz * sizeof(*row_idx));
  double *values = (double *)malloc((size_t)nnz * sizeof(*values));
  UF_long *input_to_csc = (UF_long *)malloc((size_t)nnz * sizeof(*input_to_csc));
  UF_long *base_to_trial = (UF_long *)malloc((size_t)nnz * sizeof(*base_to_trial));
  kls_row_permuted_entry *entries = max_width > 32u
    ? (kls_row_permuted_entry *)malloc((size_t)max_width * sizeof(*entries))
    : NULL;
  if (col_ptr == NULL || row_idx == NULL || values == NULL ||
      input_to_csc == NULL || base_to_trial == NULL ||
      (max_width > 32u && entries == NULL)) {
    free(col_ptr);
    free(row_idx);
    free(values);
    free(input_to_csc);
    free(base_to_trial);
    free(entries);
    return KLS_ERR_OUT_OF_MEMORY;
  }

  memcpy(col_ptr, base_col_ptr, ((size_t)n + 1u) * sizeof(*col_ptr));

  /* Row permutation does not move an entry between columns.  Keep the
     already-compressed column layout and sort only within each column;
     globally sorting (column,row) paid O(nnz log nnz) to rediscover the
     input column boundaries.  Sparse circuit columns are usually short, so
     insertion sort avoids thousands of qsort calls; retain qsort for the
     uncommon long column. */
  for (UF_long col = 0; col < n; ++col) {
    const UF_long begin = base_col_ptr[col];
    const UF_long end = base_col_ptr[col + 1u];
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
      row_idx[p] = row_perm[row];
      values[p] = base_values[p];
      /* During the per-column sort this is trial-position -> base-position;
         it is inverted below after all columns are complete. */
      input_to_csc[p] = p;
    }
    const UF_long width = end - begin;
    if (width <= 32u) {
      for (UF_long p = begin + (width != 0u); p < end; ++p) {
        const UF_long key_row = row_idx[p];
        const double key_value = values[p];
        const UF_long key_base = input_to_csc[p];
        UF_long q = p;
        while (q > begin && row_idx[q - 1u] > key_row) {
          row_idx[q] = row_idx[q - 1u];
          values[q] = values[q - 1u];
          input_to_csc[q] = input_to_csc[q - 1u];
          --q;
        }
        row_idx[q] = key_row;
        values[q] = key_value;
        input_to_csc[q] = key_base;
      }
    } else {
      for (UF_long local = 0u; local < width; ++local) {
        const UF_long p = begin + local;
        entries[local].col = col;
        entries[local].row = row_idx[p];
        entries[local].base_position = input_to_csc[p];
        entries[local].value = values[p];
      }
      qsort(entries, (size_t)width, sizeof(*entries),
            compare_row_permuted_entries);
      for (UF_long local = 0u; local < width; ++local) {
        const UF_long p = begin + local;
        row_idx[p] = entries[local].row;
        values[p] = entries[local].value;
        input_to_csc[p] = entries[local].base_position;
      }
    }
  }
  for (UF_long p = 0; p < nnz; ++p) {
    base_to_trial[input_to_csc[p]] = p;
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
  const UF_long nnz_total = base_col_ptr[n];
  double *entry_log =
    (double *)malloc((size_t)nnz_total * sizeof(*entry_log));
  if (matched_log == NULL || col_potential == NULL || entry_log == NULL) {
    free(matched_log);
    free(col_potential);
    free(entry_log);
    return 0;
  }

  for (UF_long i = 0; i < n; ++i) {
    if (diag[i] <= 0.0 || !isfinite(diag[i])) {
      free(matched_log);
      free(col_potential);
      free(entry_log);
      return 0;
    }
    matched_log[i] = log(diag[i]);
  }
  /* The refinement rounds below revisit every entry; hoisting the entry
     logarithms out of the rounds removes the dominant libm cost of the
     medium static-match setup. */
  for (UF_long p = 0; p < nnz_total; ++p) {
    const double abs_value = fabs(base_values[p]);
    entry_log[p] = abs_value > 0.0 && isfinite(abs_value)
                     ? log(abs_value)
                     : (double)NAN;
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
        const double log_value = entry_log[p];
        if (!isfinite(log_value)) {
          continue;
        }
        const UF_long matched_col = row_perm[row];
        const double lower_bound =
          col_potential[matched_col] + log_value - matched_log[matched_col];
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
      free(entry_log);
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
      const double log_value = entry_log[p];
      if (!isfinite(log_value)) {
        continue;
      }
      const UF_long matched_col = row_perm[row];
      const double excess = log_value - matched_log[matched_col] -
                            col_potential[col] + col_potential[matched_col];
      if (excess > max_excess) {
        max_excess = excess;
      }
    }
  }
  if (!converged && max_excess > log(4.0)) {
    free(matched_log);
    free(col_potential);
    free(entry_log);
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
      free(entry_log);
      return 0;
    }
    row_scale[i] = clamp_matching_scale(exp(row_log_scale));
    col_scale[i] = clamp_matching_scale(exp(col_log_scale));
  }

  free(matched_log);
  free(col_potential);
  free(entry_log);
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
  const kls_options *options,
  const trilinos_klu_l_common *base_common,
  const trilinos_klu_l_numeric *base_numeric,
  const trilinos_klu_l_common *scaled_common,
  const trilinos_klu_l_numeric *scaled_numeric) {
  const double base_fill = (double)(base_numeric->lnz + base_numeric->unz);
  const double scaled_fill = (double)(scaled_numeric->lnz + scaled_numeric->unz);
  const double base_flops = base_common->flops;
  const double scaled_flops = scaled_common->flops;
  const int materially_cheaper =
    scaled_fill < 0.98 * base_fill ||
    scaled_common->noffdiag < base_common->noffdiag ||
    (base_flops > 0.0 && scaled_flops > 0.0 &&
     scaled_flops < 0.98 * base_flops);
  const int decisive_subroundoff_rescue =
    getenv("KLS_DISABLE_DECISIVE_SUBROUNDOFF_MATCH_SCALE") == NULL &&
    base_common->rcond > 0.0 &&
    base_common->rcond < 32.0 * DBL_EPSILON &&
    scaled_common->rcond >= 100.0 * base_common->rcond &&
    scaled_fill <= 0.80 * base_fill &&
    base_flops > 0.0 && scaled_flops > 0.0 &&
    scaled_flops <= 0.80 * base_flops &&
    scaled_common->noffdiag <= base_common->noffdiag &&
    scaled_common->kls_perturb_count <= base_common->kls_perturb_count;

  if (scaled_common->noffdiag > base_common->noffdiag ||
      scaled_fill > 1.02 * base_fill ||
      (base_flops > 0.0 && scaled_flops > 1.05 * base_flops) ||
      (!(scaled_common->rcond > 0.0)) ||
      (scaled_common->rcond < 32.0 * DBL_EPSILON &&
       !decisive_subroundoff_rescue) ||
      (base_common->rcond > 0.0 &&
       scaled_common->rcond < 1.0e-4 * base_common->rcond)) {
    /* Matching equilibration is optional.  A marginal fill or pivot win
       cannot justify publishing a representation whose condition estimate
       collapses into roundoff, or by over four orders of magnitude relative
       to the already-valid matched factor.  When the base estimate is
       already below the conventional floor, however, scaling may rescue it
       only after proving simultaneous 100x conditioning and 20% fill/work
       improvements without adding pivots or perturbations.  The ordinary
       residual/refinement contract still polices the installed numeric. */
    return 0;
  }
  if (base_common->rcond <= 0.0 && scaled_common->rcond > 0.0) {
    return 1;
  }
  if (base_common->rcond > 0.0 &&
      scaled_common->rcond > 10.0 * base_common->rcond) {
    const int recurring_unscaled_route =
      base_common->scale <= 0 && kls_repeated_update_workload(options) &&
      options->expected_solves > 0;
    if (!recurring_unscaled_route ||
        base_common->rcond < sqrt(DBL_EPSILON) || materially_cheaper) {
      return 1;
    }
    /* Explicit scaling costs a value transform on every update and prevents
       the unscaled matched factor from folding its row permutation into the
       triangular solve.  Above the conventional sqrt(epsilon) conditioning
       floor, do not pay those recurring passes for an rcond-only improvement;
       require the scaled numeric to reduce realized factor work as well. */
  }
  if (base_common->scale > 0 &&
      (base_common->rcond <= 0.0 ||
       scaled_common->rcond >= 0.10 * base_common->rcond)) {
    return 1;
  }
  return materially_cheaper;
}

static int matching_equilibration_has_threshold_gain(
  UF_long n,
  const UF_long *col_ptr,
  const UF_long *row_idx,
  const double *base_values,
  const double *row_scale,
  const double *col_scale,
  double pivot_tolerance) {
  if (n <= 0u || col_ptr == NULL || row_idx == NULL ||
      base_values == NULL || row_scale == NULL || col_scale == NULL ||
      !(pivot_tolerance > 0.0)) {
    return 0;
  }
  UF_long improved = 0u;
  UF_long harmed = 0u;
  for (UF_long col = 0u; col < n; ++col) {
    const double cs = col_scale[col];
    double base_max = 0.0;
    double scaled_max = 0.0;
    double base_diag = 0.0;
    double scaled_diag = 0.0;
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      const UF_long row = row_idx[p];
      if (row >= n) {
        return 1;
      }
      const double base = fabs(base_values[p]);
      const double scaled = base * fabs(row_scale[row] * cs);
      base_max = base > base_max ? base : base_max;
      scaled_max = scaled > scaled_max ? scaled : scaled_max;
      if (row == col) {
        base_diag = base > base_diag ? base : base_diag;
        scaled_diag = scaled > scaled_diag ? scaled : scaled_diag;
      }
    }
    const int base_weak =
      !(base_diag >= pivot_tolerance * base_max);
    const int scaled_weak =
      !(scaled_diag >= pivot_tolerance * scaled_max);
    improved += (UF_long)(base_weak && !scaled_weak);
    harmed += (UF_long)(!base_weak && scaled_weak);
  }
  return improved > harmed;
}

/* Count matched diagonal entries whose column-relative magnitude lies
   between two pivot thresholds.  These are the only columns for which the
   tighter threshold can directly preserve the matched diagonal instead of
   taking the ordinary threshold-pivot detour.  The scan is allocation-free
   and uses the realized matched value frame, not input dimensions or a
   matrix-family proxy. */
static UF_long kls_matched_threshold_sensitive_columns(
  UF_long n,
  const UF_long *col_ptr,
  const UF_long *row_idx,
  const double *values,
  double tight_tolerance,
  double loose_tolerance) {
  if (n == 0u || col_ptr == NULL || row_idx == NULL || values == NULL ||
      !(tight_tolerance > 0.0) ||
      !(loose_tolerance > tight_tolerance)) {
    return 0u;
  }
  UF_long sensitive = 0u;
  for (UF_long col = 0u; col < n; ++col) {
    double column_max = 0.0;
    double diagonal = 0.0;
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      const double magnitude = fabs(values[p]);
      if (magnitude > column_max) {
        column_max = magnitude;
      }
      if (row_idx[p] == col && magnitude > diagonal) {
        diagonal = magnitude;
      }
    }
    if (column_max > 0.0 &&
        diagonal >= tight_tolerance * column_max &&
        diagonal < loose_tolerance * column_max) {
      sensitive++;
    }
  }
  return sensitive;
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
  const int safe_recurring_unscaled =
    common_io->scale <= 0 && common_io->rcond >= sqrt(DBL_EPSILON) &&
    kls_repeated_update_workload(options);
  if (safe_recurring_unscaled &&
      getenv("KLS_DISABLE_GENERIC_MATCHING_EQUILIBRATION_PREFILTER") ==
        NULL &&
      !matching_equilibration_has_threshold_gain(
        n, col_ptr, row_idx, *values_io, *row_scale_io, *col_scale_io,
        options->pivot_tolerance)) {
    /* A safely conditioned matched factor whose explicit equilibration does
       not rescue any threshold-weak diagonal has no numerical reason to pay
       for a second factorization.  Preserve the unscaled recurring route;
       the complete numeric comparison remains available whenever scaling
       changes the pivot-pressure signal or rcond needs rescue. */
    free(*row_scale_io);
    free(*col_scale_io);
    *row_scale_io = NULL;
    *col_scale_io = NULL;
    return;
  }

  double *scaled_values =
    (double *)malloc((size_t)nnz * sizeof(*scaled_values));
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
  if (getenv("KLS_TRACE_PRESTATIC") != NULL) {
    fprintf(stderr,
            "KLS pre-static: matching scale numeric base fill=%ld "
            "work=%.3e pivots=%ld rcond=%.3e scaled fill=%ld "
            "work=%.3e pivots=%ld rcond=%.3e\n",
            (long)((*numeric_io)->lnz + (*numeric_io)->unz),
            common_io->flops, (long)common_io->noffdiag,
            common_io->rcond,
            (long)(scaled_numeric->lnz + scaled_numeric->unz),
            scaled_common.flops, (long)scaled_common.noffdiag,
            scaled_common.rcond);
  }
  if (!matching_equilibration_is_better(options, common_io, *numeric_io,
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
  kls_separator_analysis trial_separator;
  memset(&trial_separator, 0, sizeof(trial_separator));
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
  int incumbent_policy_trial =
        getenv("KLS_DISABLE_EXACT_MATCH_INCUMBENT_POLICY_TRIAL") == NULL &&
    solver->stats.selected_ordering >= KLS_ORDERING_AMD &&
    solver->stats.selected_ordering <= KLS_ORDERING_METIS;
select_exact_match_symbolic:;
  if (incumbent_policy_trial) {
    trial_ordering = solver->stats.selected_ordering;
    trial_options.ordering = trial_ordering;
    trial_options.use_btf = solver->symbolic->do_btf ? 1 : 0;
    status = analyze_with_ordering(
      solver->n, trial_col_ptr, trial_row_idx, &trial_options,
      trial_ordering, &trial_symbolic, &trial_common, &trial_separator);
    if (status == KLS_OK && trial_symbolic != NULL) {
      trial_score = symbolic_score(trial_symbolic);
    }
  } else {
    const int saved_value_matched_ordering_ctx =
      kls_value_matched_ordering_ctx;
    kls_value_matched_ordering_ctx = 1;
    status = choose_symbolic_for_pattern(
      solver->n, trial_col_ptr, trial_row_idx, &trial_options,
      &trial_symbolic, &trial_common, &trial_ordering, &trial_score,
      &trial_separator);
    kls_value_matched_ordering_ctx = saved_value_matched_ordering_ctx;
  }
  if (status != KLS_OK) {
    if (!incumbent_policy_trial) {
      goto done;
    }
    if (trial_symbolic != NULL) {
      trilinos_klu_l_free_symbolic(&trial_symbolic, &trial_common);
    }
    kls_separator_analysis_clear(&trial_separator);
    (void)trilinos_klu_l_defaults(&trial_common);
    trial_options = solver->options;
    trial_ordering = KLS_ORDERING_AUTO;
    trial_score = 0.0;
    incumbent_policy_trial = 0;
    goto select_exact_match_symbolic;
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
    if (!incumbent_policy_trial) {
      goto done;
    }
    if (trial_numeric != NULL) {
      trilinos_klu_l_free_numeric(&trial_numeric, &trial_common);
    }
    trilinos_klu_l_free_symbolic(&trial_symbolic, &trial_common);
    kls_separator_analysis_clear(&trial_separator);
    (void)trilinos_klu_l_defaults(&trial_common);
    trial_options = solver->options;
    trial_ordering = KLS_ORDERING_AUTO;
    trial_score = 0.0;
    incumbent_policy_trial = 0;
    goto select_exact_match_symbolic;
  }

  (void)trilinos_klu_l_flops(trial_symbolic, trial_numeric, &trial_common);
  (void)trilinos_klu_l_rcond(trial_symbolic, trial_numeric, &trial_common);
  if (incumbent_policy_trial) {
    const int proposal_better =
      numeric_candidate_is_better(&solver->common, solver->numeric,
                                  &trial_common, trial_numeric);
    const int proposal_has_value =
      spral_hungarian_candidate_has_value(
        &solver->common, solver->numeric, &trial_common, trial_numeric);
    const int proposal_conditioning_ok =
      !(solver->common.rcond > 0.0 && trial_common.rcond > 0.0 &&
        trial_common.rcond < 0.01 * solver->common.rcond);
    if (!proposal_better || !proposal_has_value ||
        !proposal_conditioning_ok) {
      /* Reusing the incumbent ordering/BTF policy is a latency proposal,
         not a reduced search space.  If its actual matched numeric cannot
         clear the same lifecycle and conditioning gates as AUTO, discard it
         and run the full value-matched symbolic portfolio on the already
         built permutation.  This preserves the general fallback while the
         common case avoids re-solving an ordering decision that matching did
         not change. */
      trilinos_klu_l_free_numeric(&trial_numeric, &trial_common);
      trilinos_klu_l_free_symbolic(&trial_symbolic, &trial_common);
      kls_separator_analysis_clear(&trial_separator);
      (void)trilinos_klu_l_defaults(&trial_common);
      trial_options = solver->options;
      trial_ordering = KLS_ORDERING_AUTO;
      trial_score = 0.0;
      incumbent_policy_trial = 0;
      goto select_exact_match_symbolic;
    }
  }
  const int retain_native_matching_scale =
    kls_repeated_update_workload(&solver->options) &&
    getenv("KLS_DISABLE_GENERIC_NATIVE_MATCHING_SCALE") == NULL &&
    trial_common.rcond >= 32.0 * DBL_EPSILON &&
    (solver->common.rcond <= 0.0 ||
     trial_common.rcond >= solver->common.rcond) &&
    trial_common.noffdiag <= solver->common.noffdiag &&
    numeric_candidate_is_better(&solver->common, solver->numeric,
                                &trial_common, trial_numeric) &&
    spral_hungarian_candidate_has_value(&solver->common, solver->numeric,
                                        &trial_common, trial_numeric);
  if (retain_native_matching_scale) {
    /* Exact matching has already produced a conventionally scaled KLU
       factor with no conditioning or pivot-quality loss relative to the
       incumbent and a decisive realized lifecycle win.  A second factor in
       the external matching equilibration has identical coordinates and
       arithmetic, but makes every later refactor refresh two extra dense
       scale vectors.  Keep the numerically certified native frame; the
       ordinary solve contract continues to police subsequent values. */
    free(row_scale);
    free(col_scale);
    row_scale = NULL;
    col_scale = NULL;
  } else {
    maybe_use_matching_equilibration(solver->n, solver->nnz, trial_col_ptr,
                                     trial_row_idx, &solver->options,
                                     trial_symbolic, &trial_values,
                                     &row_scale, &col_scale, &trial_numeric,
                                     &trial_common);
  }
  const int numeric_better =
    numeric_candidate_is_better(&solver->common, solver->numeric,
                                &trial_common, trial_numeric);
  const int lifecycle_value =
    spral_hungarian_candidate_has_value(&solver->common, solver->numeric,
                                        &trial_common, trial_numeric);
  const int conditioning_ok =
    !(solver->common.rcond > 0.0 && trial_common.rcond > 0.0 &&
      trial_common.rcond < 0.01 * solver->common.rcond);
  if (!numeric_better || !lifecycle_value || !conditioning_ok) {
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
  solver->exact_matching_scaling_selected = 0;
  solver->spral_matching_selected = 1;
  solver->symbolic = trial_symbolic;
  kls_separator_analysis_move(&solver->separator, &trial_separator);
  kls_invalidate_factor_etree_stats(solver);
  solver->numeric = trial_numeric;
  if (trial_common.kls_perturb_count > 0) {
    /* the hard-static factor is a nearby matrix's exact LU; the
       enrolled solve refinement recovers the accuracy contract */
    solver->numeric_needs_refinement = 1;
  }
  kls_numeric_replaced_invalidate(solver);
  solver->common = trial_common;
  solver->stats.selected_ordering = trial_ordering;
  solver->stats.selected_orientation = solver->orientation;
  solver->symmetric_partial_diagonal_low_work_class =
    0;
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
    kls_separator_analysis_clear(&trial_separator);
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

/* A value-aware row match is a performance challenger once a valid numeric
   exists.  If pivoting has already produced a factor inside the symbolic
   fill/work uncertainty band, the symbolic model contains no measured
   inflation for another O(nnz) match, ordering, and full factor to repair.
   Share this realized certificate across every later matching
   implementation, including the exact-assignment fallback, so a cheaper
   pre-static decline is not followed by the same expensive portfolio. */
static int kls_live_numeric_declines_matching_portfolio(
  const kls_solver *solver) {
  if (solver == NULL ||
      solver->numeric == NULL || solver->symbolic == NULL ||
      !(solver->common.flops > 0.0) ||
      !(solver->symbolic->est_flops > 0.0)) {
    return 0;
  }
  const double numeric_fill =
    (double)(solver->numeric->lnz + solver->numeric->unz);
  const double symbolic_fill =
    solver->symbolic->lnz + solver->symbolic->unz;
  /* The symbolic work estimate excludes numerical pivot detours.  A live
     numeric still within a 35% work band, while using no more fill than the
     symbolic envelope, has not demonstrated enough pivot inflation to fund
     a second matching, ordering, and factorization portfolio.  The prior
     25% band remains as an A/B control. */
  const double modeled_work_band =
    getenv("KLS_DISABLE_RELAXED_LIVE_MATCH_CERTIFICATE") == NULL
      ? 1.35 : 1.25;
  return numeric_fill > 0.0 && symbolic_fill > 0.0 &&
    numeric_fill <= 1.10 * symbolic_fill &&
    solver->common.flops <= modeled_work_band * solver->symbolic->est_flops;
}

static int should_try_spral_hungarian_numeric_trial(
  const kls_solver *solver,
  int ignore_pending_portfolio) {
  if (getenv("KLS_DISABLE_SPRAL_HUNGARIAN_TRIAL") != NULL ||
      /* A repeated full-factor sample can arrive while the first factor's
         broader value-aware representation portfolio is deliberately
         waiting for a changed numeric.  Do not let this one exact-match arm
         publish early and make the pending pre-static/row-match candidates
         ineligible merely because they require the original coordinates.
         The deferred consultation compares those complete numerics before
         adoption; a factor epoch with no pending portfolio retains the
         ordinary Hungarian challenger. */
      (solver != NULL &&
       !ignore_pending_portfolio &&
       getenv("KLS_DISABLE_PENDING_MATCH_PORTFOLIO_GUARD") == NULL &&
       (solver->prestatic_deferred || solver->rowmatch_deferred)) ||
      solver == NULL || !solver->options.static_pivoting ||
      solver->numeric == NULL || solver->symbolic == NULL ||
      kls_live_numeric_declines_matching_portfolio(solver) ||
      solver->row_perm != NULL || solver->input_format != KLS_INPUT_CSC ||
      solver->options.ordering != KLS_ORDERING_AUTO ||
      solver->n < 20000u || solver->n > 750000u ||
      solver->nnz > 8000000u || solver->common.noffdiag < 16u) {
    return 0;
  }

  if (solver->dense_spiked_original_pivot_path) {
    /* On dense-spike systems the incumbent pivoted METIS numeric already
       supports the faster repeated row engine.  Converting it to a matched
       static numeric costs a full symbolic/factor setup and lengthens the
       H100 solve/refactor pair even when the matched trial is accepted. */
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
      !should_try_spral_hungarian_numeric_trial(solver, 0)) {
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

/* the pattern-only half of the refinement gate: everything knowable
   before the trial factor runs (the flops floor is checked at join) */

/* Speculative nested-dissection analysis for the pre-static METIS
   refinement: the serial NodeND+CAMD dominates the refinement cost and
   has no dependency on the trial factor's numeric results, so it can
   run on a worker concurrently with the trial factorization.  The join
   in maybe_refine_pre_static_with_metis applies the full (flops-aware)
   gate before using the result. */
typedef struct kls_psmetis_spec {
  pthread_t thread;
  int launched;
  int status;
  UF_long n;
  UF_long *col_ptr;
  UF_long *row_idx;
  double *values;           /* read-only during the speculation */
  int do_factor;            /* also run the METIS trial factor on the
                               worker (only when the main-side trial
                               is serial, keeping the thread budget) */
  int factor_attempted;
  trilinos_klu_l_numeric *numeric;
  kls_options options;
  trilinos_klu_l_symbolic *symbolic;
  trilinos_klu_l_common common;
  kls_separator_analysis separator;
} kls_psmetis_spec;

static void kls_psmetis_spec_join(kls_psmetis_spec *spec) {
  if (spec != NULL && spec->launched) {
    pthread_join(spec->thread, NULL);
    spec->launched = 0;
  }
}

static void kls_psmetis_spec_discard(kls_psmetis_spec *spec) {
  if (spec == NULL) {
    return;
  }
  kls_psmetis_spec_join(spec);
  if (spec->numeric != NULL) {
    trilinos_klu_l_free_numeric(&spec->numeric, &spec->common);
  }
  if (spec->status == KLS_OK && spec->symbolic != NULL) {
    trilinos_klu_l_free_symbolic(&spec->symbolic, &spec->common);
  }
  kls_separator_analysis_clear(&spec->separator);
  spec->symbolic = NULL;
  spec->numeric = NULL;
  spec->status = KLS_ERR_FACTOR_FAILED;
}

static int should_try_auto_row_match(const kls_solver *solver,
                                     const double *numeric_values) {
  /* Weighted static pivoting is currently a reactive medium-matrix trial.  On
     larger matrices the O(nnz log nnz) matching/reanalysis cost needs a cheaper
     precheck before it is worth paying by default. */
  if (solver == NULL || !solver->options.static_pivoting ||
      solver->numeric == NULL || solver->row_perm != NULL ||
      kls_live_numeric_declines_matching_portfolio(solver) ||
      solver->input_format != KLS_INPUT_CSC ||
      solver->options.ordering != KLS_ORDERING_AUTO || solver->n < 3000 ||
      solver->n > 20000 || solver->common.noffdiag < 16 ||
      numeric_values == NULL) {
    return 0;
  }

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
  kls_separator_analysis trial_separator;
  memset(&trial_separator, 0, sizeof(trial_separator));
  (void)trilinos_klu_l_defaults(&trial_common);
  double trial_factor_seconds = 0.0;
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
  int exact_matching_scaling = 0;
  int spral_matching = 0;
  const int improve_matching = solver->n <= 50000 && solver->nnz <= 1000000;
  int status = build_greedy_numeric_row_match(solver->n, solver->nnz,
                                              base_col_ptr, base_row_idx,
                                              base_values,
                                              0,
                                              improve_matching,
                                              &row_perm, &matched,
                                              &exact_matching,
                                              &exact_matching_scaling,
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
  const int saved_value_matched_ordering_ctx =
    kls_value_matched_ordering_ctx;
  kls_value_matched_ordering_ctx = 1;
  status = choose_symbolic_for_pattern(solver->n, trial_col_ptr, trial_row_idx,
                                       &trial_options, &trial_symbolic,
                                       &trial_common, &trial_ordering,
                                       &trial_score, &trial_separator);
  kls_value_matched_ordering_ctx = saved_value_matched_ordering_ctx;
  if (status != KLS_OK) {
    goto done;
  }
  trial_common.scale = choose_auto_scale_from_pattern(solver->n, trial_col_ptr,
                                                      trial_row_idx,
                                                      &trial_options,
                                                      trial_values);

  /* the trial factor is the largest single cost of the prestatic
     one-shot (pre2: 3.4s of the 8.9s init); the deterministic pipe
     factors it in parallel with identical acceptance semantics */
  kls_klu_pipe_threads =
    kls_pipe_first_factor_threads(solver, trial_symbolic);
  /* the trial tolerates the full pipe (dense panels + batched
     consume): its output is judged by the quality gates below and
     the class's own dense-panel win is the point (pre2) */
  kls_klu_pipe_det = 0;
  const double trial_factor_start = kls_now_seconds();
  trial_numeric =
    trilinos_klu_l_factor(trial_col_ptr, trial_row_idx, trial_values,
                          trial_symbolic, &trial_common);
  trial_factor_seconds = kls_now_seconds() - trial_factor_start;
  kls_klu_pipe_threads = 0;
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
  solver->exact_matching_scaling_selected =
    exact_matching_scaling && trial_row_scale != NULL &&
    trial_col_scale != NULL;
  solver->spral_matching_selected = spral_matching;
  solver->symbolic = trial_symbolic;
  kls_separator_analysis_move(&solver->separator, &trial_separator);
  kls_invalidate_factor_etree_stats(solver);
  solver->numeric = trial_numeric;
  kls_numeric_replaced_invalidate(solver);
  solver->numeric_full_factor_seconds = trial_factor_seconds;
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
    kls_separator_analysis_clear(&trial_separator);
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

/* pipelined first-factor routing (see trilinos_klu_l_kernel.c): the
   pipeline wins on heavy columns (pre2: 160k flops/col, init 57->32s)
   and loses on light ones (mac_econ 24k: 1.9x round overhead, plus
   cancellation-pivot fragility); route foreground first factors with
   est work/column above a floor.  Thread-local, so race workers and
   other background factors never engage it. */

/* Worker count scaled to the estimated factor work: each pipe worker
   pays workspace setup (n-sized flag/X arrays) and pointer-chase
   startup that only amortizes over enough numeric work.  ~5e8 est
   flops per worker keeps big factors at full width (pre2 1e11 -> 16)
   while small matrices stop overspawning (onetone1 1.5e9 -> 3:
   16 workers measured 0.221s vs 0.208 at 4). */
static int kls_pipe_scale_threads(int threads, double est_flops) {
  int cap = threads > 16 ? 16 : threads;
  if (est_flops > 0.0) {
    const double by_work = est_flops / 1.25e8;
    if (by_work < (double)cap) {
      cap = by_work < 2.0 ? 2 : (int)by_work;
    }
  }
  return cap;
}

static int kls_pipe_first_factor_threads(const kls_solver *solver,
                                         const trilinos_klu_l_symbolic *sym) {
  if (solver == NULL || sym == NULL || solver->options.threads < 2 ||
      getenv("KLS_DISABLE_PIPE_ROUTE") != NULL) {
    return 0;
  }
  double est = sym->est_flops;
  const double n = (double)sym->n;
  if (!(est > 0.0) && sym->lnz > 0.0 && sym->unz > 0.0) {
    /* some analyze paths leave est_flops unset (mac: the orientation
       trials); a fill-based proxy decides instead - flops scale
       superlinearly in fill, and 3e7 estimated entries puts the
       serial factor in the tens of seconds where the pipe pays */
    est = (sym->lnz + sym->unz) > 3.0e7 ? 6.0e9 : 0.0;
  }
  if (getenv("KLS_TRACE_PIPE_ROUTE") != NULL) {
    fprintf(stderr, "KLS pipe route: est=%.3e n=%.0f fill=%.3e\n",
            est, n, sym->lnz + sym->unz);
  }
  if (getenv("KLS_KLU_PIPE_FORCE") == NULL && 1 &&
      !(est > 5.0e9)) {
    /* the old est/n >= 1e5 companion guard existed to keep the pipe
       away from the cancellation classes (mac) whose accuracy it
       destroyed; that was the pivot-race/apply-order defect, fixed by
       the segment sort + deterministic tie-break.  The routed pipe
       runs the deterministic scalar configuration (kls_klu_pipe_det)
       and the flops floor alone decides profitability. */
    return 0;
  }
  /* full configuration (dense + batched consume): one-shot draws
     measured e-16..e-9 on mac, and with fp64 refactors (the
     numeric_from_pipe gate) the end-to-end contract holds; the
     deterministic-scalar pin predates those fixes */
  kls_klu_pipe_det = 0;
  return kls_pipe_scale_threads(solver->options.threads, est);
}

static int kls_trace_pre_static_enabled(void) {
  const char *value = getenv("KLS_TRACE_PRESTATIC");
  return value != NULL && value[0] == '1';
}

/* Defined after the METIS-race object, whose analyze stage is launched from
   kls_analyze before numeric values (and hence weighted matching) exist. */
static int kls_prestatic_reuse_raced_metis_symbolic(
  kls_solver *solver,
  UF_long *col_ptr,
  UF_long *row_idx,
  const kls_options *options,
  trilinos_klu_l_symbolic **symbolic_out,
  trilinos_klu_l_common *common_out,
  kls_separator_analysis *separator_out);

static double kls_column_pair_work(UF_long n,
                                   const UF_long *col_ptr) {
  if (col_ptr == NULL) {
    return DBL_MAX;
  }
  double work = 0.0;
  for (UF_long col = 0u; col < n; ++col) {
    const double degree =
      (double)(col_ptr[col + 1u] - col_ptr[col]);
    work += degree * degree;
  }
  return work;
}

/* A static row match can improve pivot stability while making the retained
   ordering more expensive than an ordinary ordering on the public frame.
   Run one bounded original-frame challenger alongside the matched factor so
   the comparison usually overlaps the critical path.  Retained-factor and
   quotient-pair work bound the trial up front; complete numeric fill/work,
   conditioning, and a backward residual decide whether the result may be
   published. */
typedef struct kls_unmatched_colamd_job {
  UF_long n;
  const UF_long *col_ptr;
  const UF_long *row_idx;
  const double *values;
  kls_options options;
  double matched_estimated_fill;
  trilinos_klu_l_symbolic *symbolic;
  trilinos_klu_l_numeric *numeric;
  trilinos_klu_l_common common;
  kls_separator_analysis separator;
  double score;
  double factor_seconds;
  int status;
} kls_unmatched_colamd_job;

static void *kls_unmatched_colamd_main(void *arg) {
  kls_unmatched_colamd_job *job =
    (kls_unmatched_colamd_job *)arg;
  memset(&job->separator, 0, sizeof(job->separator));
  job->symbolic = NULL;
  job->numeric = NULL;
  job->score = DBL_MAX;
  job->factor_seconds = 0.0;
  job->status = analyze_with_ordering(
    job->n, (UF_long *)(uintptr_t)job->col_ptr,
    (UF_long *)(uintptr_t)job->row_idx, &job->options,
    KLS_ORDERING_COLAMD, &job->symbolic, &job->common,
    &job->separator);
  if (job->status != KLS_OK || job->symbolic == NULL) {
    return NULL;
  }
  job->score = symbolic_score(job->symbolic);
  if (getenv("KLS_TRACE_PRESTATIC") != NULL) {
    fprintf(stderr,
            "KLS pre-static: unmatched COLAMD symbolic fill=%.4e "
            "matched=%.4e\n",
            job->score, job->matched_estimated_fill);
  }
  /* COLAMD does not publish a fill estimate in every build.  Use an estimate
     when it exists, but do not turn an unavailable estimate into a rejection:
     the caller has already bounded the quotient-pair analysis work and the
     complete candidate factor is the authoritative cost measurement. */
  if (isfinite(job->score) && job->score < DBL_MAX &&
      (job->score > 2000000.0 ||
       job->score > 1.05 * job->matched_estimated_fill)) {
    return NULL;
  }
  job->common.scale =
    job->options.scale == KLS_SCALE_AUTO
      ? choose_auto_scale_from_pattern(
          job->n, job->col_ptr, job->row_idx, &job->options,
          job->values)
      : job->options.scale;
  job->common.tol = job->options.pivot_tolerance;
  const double start = kls_now_seconds();
  job->numeric = trilinos_klu_l_factor(
    (UF_long *)(uintptr_t)job->col_ptr,
    (UF_long *)(uintptr_t)job->row_idx,
    (double *)(uintptr_t)job->values, job->symbolic, &job->common);
  job->factor_seconds = kls_now_seconds() - start;
  if (job->numeric == NULL || job->common.status < TRILINOS_KLU_OK ||
      job->common.status == TRILINOS_KLU_SINGULAR) {
    return NULL;
  }
  (void)trilinos_klu_l_flops(job->symbolic, job->numeric,
                             &job->common);
  (void)trilinos_klu_l_rcond(job->symbolic, job->numeric,
                             &job->common);
  return NULL;
}

static void kls_unmatched_colamd_discard(
  kls_unmatched_colamd_job *job) {
  if (job == NULL) {
    return;
  }
  if (job->numeric != NULL) {
    trilinos_klu_l_free_numeric(&job->numeric, &job->common);
  }
  if (job->symbolic != NULL) {
    trilinos_klu_l_free_symbolic(&job->symbolic, &job->common);
  }
  kls_separator_analysis_clear(&job->separator);
}

static void maybe_select_pre_static_row_match(kls_solver *solver,
                                              double *elapsed,
                                              const double *numeric_values,
                                              int deferred) {
  if (solver == NULL) {
    return;
  }

  const int forced_match = getenv("KLS_FORCE_STATIC_MATCH") != NULL;

  const int trace_gates =
    deferred && kls_trace_pre_static_enabled();
  if (!solver->options.static_pivoting ||
      (solver->numeric != NULL && !deferred) ||
      solver->row_perm != NULL ||
      solver->input_format != KLS_INPUT_CSC ||
      (solver->options.ordering != KLS_ORDERING_AUTO && !forced_match) ||
      (solver->n < 3000 && !forced_match)) {
    if (trace_gates) {
      fprintf(stderr,
              "KLS pre-static defer: entry gate declined sp=%d perm=%d"
              " fmt=%d ord=%d n=%ld\n",
              solver->options.static_pivoting, solver->row_perm != NULL,
              (int)solver->input_format, (int)solver->options.ordering,
              (long)solver->n);
    }
    return;
  }

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
  kls_separator_analysis trial_separator;
  memset(&trial_separator, 0, sizeof(trial_separator));
  (void)trilinos_klu_l_defaults(&trial_common);
  trilinos_klu_l_symbolic *alternate_symbolic = NULL;
  trilinos_klu_l_numeric *alternate_numeric = NULL;
  trilinos_klu_l_common alternate_common;
  kls_separator_analysis alternate_separator;
  memset(&alternate_separator, 0, sizeof(alternate_separator));
  kls_ordering alternate_ordering = KLS_ORDERING_AUTO;
  double alternate_score = DBL_MAX;
  double trial_factor_seconds = 0.0;
  double alternate_factor_seconds = 0.0;
  kls_unmatched_colamd_job unmatched_colamd_job;
  memset(&unmatched_colamd_job, 0, sizeof(unmatched_colamd_job));
  pthread_t unmatched_colamd_thread;
  int unmatched_colamd_launched = 0;
  int unmatched_colamd_joined = 0;
  int selected_unmatched_colamd = 0;
  int matching_equilibration_deferred = 0;
  int accepted = 0;
  int compact_missing_diagonal_matcher = 0;
  kls_psmetis_spec psmetis_spec;
  memset(&psmetis_spec, 0, sizeof(psmetis_spec));
  psmetis_spec.status = KLS_ERR_FACTOR_FAILED;
  /* large-spral estimate speculation: the METIS estimate analyze can
     overlap choose_symbolic_for_pattern (the base candidates draw no
     GKlib rand, so the legacy-RNG class's orderings stay
     deterministic with the worker as the only rand consumer) */
  kls_psmetis_spec psmetis_est_spec;
  memset(&psmetis_est_spec, 0, sizeof(psmetis_est_spec));
  psmetis_est_spec.status = KLS_ERR_FACTOR_FAILED;

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

  const double kls_ps_t0 = kls_now_seconds();
  kls_medium_partial_static_metis_ctx =
    0;
  kls_low_work_symmetric_partial_diagonal_pts_ctx = 0;
  kls_large_weak_diagonal_static_metis_ctx =
    0;
  kls_prestatic_ordering_ctx = 1;
  UF_long missing_diagonal = 0;
  UF_long weak = 0u;
  {
    weak = count_weak_diagonal_rows(solver->n, base_col_ptr, base_row_idx,
                                    base_values, solver->common.tol,
                                    &missing_diagonal);
  }
  const int nearly_missing_diagonal_early_match =
    32.0 * (double)weak >= 31.0 * (double)solver->n &&
        32.0 * (double)missing_diagonal >= 31.0 * (double)solver->n;
  if (nearly_missing_diagonal_early_match &&
      solver->nnz >= 4u * solver->n &&
      getenv("KLS_DISABLE_COMPACT_MISSING_DIAGONAL_MATCHER") == NULL) {
    /* When at least 31/32 of the diagonal is structurally absent, the generic
       matcher must construct essentially the whole transversal.  Use the
       linear compact matching builder as the proposal instead of sorting all
       numeric entries and rebuilding the same augmenting graph.  Its result
       remains subject to the ordinary coverage, factor, pivot, conditioning,
       and solve-quality acceptance gates below. */
    compact_missing_diagonal_matcher = 1;
  }
  if (!deferred) {
    solver->nearly_missing_diagonal_early_match_selected =
      nearly_missing_diagonal_early_match;
    solver->compact_missing_diagonal_match_candidate =
      compact_missing_diagonal_matcher;
  }
  /* A static row match changes the numeric coordinate frame.  Under a
     recurring changed-value contract, keeping that frame unscaled requires
     every later value consumer to preserve the exact matched coordinates.
     That is a stronger invariant than a successful first factor proves.
     Prefer the ordinary matched equilibration for generic AUTO, which has
     the same representation contract on every value generation.  The
     unscaled form remains available for explicit experiments. */
  int prefer_unscaled_static_match =
    getenv("KLS_STATIC_MATCH_SCALED") != NULL
      ? 0
      : getenv("KLS_STATIC_MATCH_UNSCALED") != NULL
          ? 1
          : 0;
  if (kls_trace_pre_static_enabled()) {
    fprintf(stderr,
            "KLS pre-static: weak count done %.3fms weak=%ld missing=%ld "
            "live_offdiag=%ld live_fill=%ld live_work=%.3e\n",
            1e3 * (kls_now_seconds() - kls_ps_t0),
            (long)weak, (long)missing_diagonal,
            (long)(solver->numeric != NULL ? solver->common.noffdiag : 0),
            (long)(solver->numeric != NULL
              ? solver->numeric->lnz + solver->numeric->unz : 0),
            solver->numeric != NULL ? solver->common.flops : 0.0);
  }
  if (deferred && !forced_match &&
      solver->rowmatch_deferred &&
      getenv("KLS_DISABLE_GENERIC_SINGLE_EXACT_MATCH_TRIAL") == NULL &&
      should_try_spral_hungarian_numeric_trial(solver, 1)) {
    /* Both deferred matchers are performance challengers to the same valid
       numeric.  If the exact sparse-assignment trial independently passes
       all of its live-numeric admission checks, do not first build and factor
       a second greedy permutation: the exact trial below retains the stronger
       matching coverage and is judged by the same fill, work, conditioning,
       and residual contracts.  Inline matching and failed-factor recovery
       never reach this gate, and an ineligible exact trial leaves the greedy
       fallback intact. */
    goto done;
  }
  if (deferred && !forced_match &&
      solver->numeric != NULL &&
      solver->common.noffdiag <= weak / 20u + 16u) {
    /* The live numeric already meets the same pivot-deviation bound required
       of an accepted match.  Matching every row cannot repair enough actual
       pivot work to justify its ordering/factor trial; keep the measured
       incumbent.  This converts weak-diagonal structure into a proposal,
       with realized pivot behavior deciding whether the proposal advances. */
    if (kls_trace_pre_static_enabled()) {
      fprintf(stderr,
              "KLS pre-static: incumbent pivot bound declines deferred "
              "match offdiag=%ld cap=%ld\n",
              (long)solver->common.noffdiag, (long)(weak / 20u + 16u));
    }
    goto done;
  }
  if (deferred && !forced_match &&
      solver->numeric != NULL && solver->symbolic != NULL) {
    const double numeric_fill =
      (double)(solver->numeric->lnz + solver->numeric->unz);
    const double symbolic_fill =
      solver->symbolic->lnz + solver->symbolic->unz;
    const double numeric_work = solver->common.flops;
    const double symbolic_work = solver->symbolic->est_flops;
    const UF_long modest_pivot_cap = solver->n / 100u + 16u;
    const int incumbent_solve_certified =
      solver->solve_contract_probe == 1;
    const int incumbent_beats_symbolic_economics =
      numeric_fill > 0.0 && symbolic_fill > 0.0 &&
      numeric_fill <= symbolic_fill && numeric_work > 0.0 &&
      symbolic_work > 0.0 && numeric_work <= symbolic_work;
    if ((solver->common.noffdiag <= modest_pivot_cap ||
         incumbent_solve_certified ||
         incumbent_beats_symbolic_economics) &&
        numeric_fill > 0.0 && symbolic_fill > 0.0 &&
        numeric_fill <= 1.10 * symbolic_fill &&
        numeric_work > 0.0 && symbolic_work > 0.0 &&
        numeric_work <= 1.25 * symbolic_work) {
      /* Matching is a response to realized pivot inflation.  If the live
         factor remains within 10% of symbolic fill and 25% of symbolic
         arithmetic, and either pivot detours stay below one percent, its
         solve contract has already been certified, or realized storage and
         work both beat their predictions, there is no demonstrated
         inflation for an O(nnz) match/reorder/full-factor proposal to
         recover.  Symbolic arithmetic is intentionally given the wider
         band: threshold pivoting commonly adds modest work without adding
         retained storage or slowing the recurring engine.  This gate uses
         numeric outcomes only; it contains no matrix-shape class. */
      if (kls_trace_pre_static_enabled()) {
        fprintf(stderr,
                "KLS pre-static: stable incumbent declines deferred match "
                "offdiag=%ld/%ld certified=%d settled=%d fill=%.3e/%.3e "
                "work=%.3e/%.3e\n",
                (long)solver->common.noffdiag, (long)modest_pivot_cap,
                incumbent_solve_certified,
                incumbent_beats_symbolic_economics, numeric_fill,
                symbolic_fill, numeric_work, symbolic_work);
      }
      goto done;
    }
  }
#ifdef KLS_HAVE_SPRAL_SCALING
  int use_large_spral_match = 0;
  const int force_generic_spral_match =
        getenv("KLS_FORCE_GENERIC_SPRAL_MATCH") != NULL;
  if (
      solver->symbolic != NULL && solver->symbolic->nblocks <= 4u &&
      5.0 * (double)missing_diagonal >= 2.0 * (double)solver->n) {
    /* When a substantial part of the diagonal is absent, augmenting the
       greedy match is the expensive path.  The exact sparse assignment is
       both bounded by the admitted candidate sizes and typically cheaper;
       its factor still has to pass all ordinary adoption checks. */
    use_large_spral_match = 1;
  }
#endif
  {
    const double matching_horizon = 1.0 +
      (double)solver->options.expected_refactorizations;
    const int inline_weak =
      2.0 * (double)weak >= (double)solver->n ||
      ((double)missing_diagonal * matching_horizon >=
       4.0 * (double)solver->nnz);
    const int meaningful_weak =
      weak >= 16u && 100.0 * (double)weak >= (double)solver->n;
    if (!forced_match && !meaningful_weak) {
      goto done;
    }
    if (!deferred && !forced_match && !inline_weak &&
        (kls_defer_cycle_trials_enabled() ||
         (kls_repeated_update_workload(&solver->options) &&
          getenv("KLS_DISABLE_GENERIC_PRESTATIC_LIFECYCLE_DEFERRAL") ==
            NULL))) {
      /* A partial weak diagonal is a recurring-cycle opportunity, not a
         reason to duplicate a healthy one-shot factor.  Reconsider it on
         the first changed numeric; a failed first factor still reaches the
         ordinary rescue path.  The value-derived weak-diagonal census
         supplies the later proposal, but neither a large factor estimate nor
         an ordering choice proves that matching will improve the realized
         numeric.  Measuring the incumbent first lets the existing fill/work,
         pivot, conditioning, and solve gates decline that speculation before
         it duplicates an expensive factor.  This generic representation has
         its own changed-value lifecycle, so retain its comparison even when
         exact-input reuse is disabled; the opt-out is deliberately local. */
      solver->prestatic_deferred = 1;
      kls_prestatic_ordering_ctx = 0;
      kls_medium_partial_static_metis_ctx = 0;
      kls_low_work_symmetric_partial_diagonal_pts_ctx = 0;
      kls_large_weak_diagonal_static_metis_ctx = 0;
      *elapsed += kls_now_seconds() - start;
      return;
    }
#ifdef KLS_HAVE_SPRAL_SCALING
    if (large_spral_candidate) {
      use_large_spral_match = 1;
    }
#endif
    goto matching_admitted;
  }
  if (!small_candidate) {
    const int majority_weak = weak * 2u >= solver->n && weak >= 5000;
    const int partial_weak =
      solver->n >= 80000 && solver->n <= 100000 &&
      solver->nnz >= 5u * solver->n && solver->nnz <= 1000000 &&
      missing_diagonal * 100u <= 3u * solver->n &&
      weak >= 2000 && weak * 100u >= solver->n;
    if (!deferred &&
        kls_defer_cycle_trials_enabled() &&
        solver->n <= 150000u && !majority_weak) {
      /* A mostly-strong diagonal means the incumbent ordering factors fine
         and the trial is a later-cycle payoff bet (rajat25: 0.53s trial
         versus 0.15s for the AMD factor).  Exact-input reuse defers that bet
         until the first changed input.  The factor
         failure rescue below still runs this consult immediately.  Majority-
         weak rows (onetone1/twotone) stay inline because the matched ordering
         is their fast and reliable initial-factor path.  A forced large-weak
         PTS match stays inline for the same reason: deferring it would factor
         the known-losing unmatched numeric before building the retained
         matched factor on the first update. */
      solver->prestatic_deferred = 1;
      kls_prestatic_ordering_ctx = 0;
      kls_medium_partial_static_metis_ctx = 0;
      kls_low_work_symmetric_partial_diagonal_pts_ctx = 0;
      kls_large_weak_diagonal_static_metis_ctx = 0;
      return;
    }
    int use_medium_gate = 1;
#ifdef KLS_HAVE_SPRAL_SCALING
    if (!medium_weak_candidate) {
      const int large_spral_weak =
        large_spral_candidate &&
        (solver->symbolic != NULL && solver->symbolic->do_btf &&
          solver->symbolic->maxblock * 4u >= solver->n * 3u &&
          ((weak >= 32768u && weak * 100u >= solver->n) ||
           (missing_diagonal >= 4096u &&
            missing_diagonal * 200u >= solver->n)));
      if (!large_spral_weak && !forced_match) {
        goto done;
      }
      use_large_spral_match = 1;
      use_medium_gate = 0;
    }
#endif
    if (use_medium_gate && !forced_match &&
        !majority_weak && !partial_weak) {
      goto done;
    }
  } else if (!forced_match && weak * 2u < solver->n) {
    goto done;
  }
  if (!deferred && small_candidate && solver->n <= 16384u &&
      !nearly_missing_diagonal_early_match &&
      kls_defer_cycle_trials_enabled()) {
    /* In the one-shot tiny class, an optimistic plain factor can be cheaper
       than an inline match and trial.  A singular plain factor is rescued at
       factor exit, and a changed-value cycle gets the deferred consult at its
       first refactorization.  A normalized nearly-missing-diagonal proposal
       is excluded by the policy flag above because its accepted matched
       factor repays matching in addition to avoiding the duplicate factor.
       The medium-weak PTS policy is likewise mandatory: deferred
       unfactored adoption cannot replace an already-live numeric safely.
     */
    solver->prestatic_deferred = 1;
    kls_prestatic_ordering_ctx = 0;
    kls_medium_partial_static_metis_ctx = 0;
    kls_low_work_symmetric_partial_diagonal_pts_ctx = 0;
    kls_large_weak_diagonal_static_metis_ctx = 0;
    *elapsed += kls_now_seconds() - start;
    return;
  }

matching_admitted:;
  double matched_column_pair_work = DBL_MAX;
  if (solver->options.orientation == KLS_ORIENTATION_AUTO &&
      solver->options.ordering == KLS_ORDERING_AUTO &&
      solver->options.scale == KLS_SCALE_AUTO &&
      solver->options.backend != KLS_BACKEND_SERIAL &&
      solver->options.threads >= 2 && solver->options.use_btf &&
      kls_repeated_update_workload(&solver->options)) {
    matched_column_pair_work =
      kls_column_pair_work(solver->n, base_col_ptr);
  }
  if (kls_trace_pre_static_enabled() &&
      isfinite(matched_column_pair_work)) {
    fprintf(stderr, "KLS pre-static: column-pair-work=%.4e\n",
            matched_column_pair_work);
  }
  UF_long matched = 0;
  int exact_matching = 0;
  int exact_matching_scaling = 0;
  int spral_matching = 0;
  int bounded_exact_match_probe_current = 0;
  int bounded_exact_match_probe_attempted = 0;
  int bounded_exact_match_probe_rejected = 0;
  int status = KLS_OK;
#ifdef KLS_HAVE_SPRAL_SCALING
  const double incumbent_symbolic_score = solver->symbolic != NULL
    ? symbolic_score(solver->symbolic) : DBL_MAX;
  const int incumbent_symbolic_already_bounded =
    solver->symbolic != NULL &&
    solver->symbolic->est_flops > 0.0 &&
    solver->symbolic->est_flops <= 1.0e7 &&
    incumbent_symbolic_score > 0.0 &&
    incumbent_symbolic_score <= 2.0e6;
  const int bounded_exact_match_probe =
    use_large_spral_match &&
    !force_generic_spral_match &&
    !incumbent_symbolic_already_bounded &&
    isfinite(matched_column_pair_work) &&
    matched_column_pair_work > 1.0e6 &&
    matched_column_pair_work <= 1.0e7 &&
    getenv("KLS_DISABLE_BOUNDED_EXACT_MATCH_PROBE") == NULL;
#endif

matching_attempt:;
  matched = 0u;
  exact_matching = 0;
  exact_matching_scaling = 0;
  spral_matching = 0;
  status = KLS_OK;
#ifdef KLS_HAVE_SPRAL_SCALING
  bounded_exact_match_probe_current =
    bounded_exact_match_probe && !bounded_exact_match_probe_attempted;
  if (force_generic_spral_match || bounded_exact_match_probe_current) {
    bounded_exact_match_probe_attempted |=
      bounded_exact_match_probe_current;
    status = build_spral_hungarian_row_match_scaling(
      solver->n, solver->nnz, base_col_ptr, base_row_idx, base_values,
      &row_perm, &matched, &match_row_scale, &match_col_scale);
    exact_matching = 0;
    exact_matching_scaling = 0;
    spral_matching = status == KLS_OK;
  } else
#endif
  {
    const int improve_matching =
      small_candidate || medium_weak_candidate;
    status = build_greedy_numeric_row_match(solver->n, solver->nnz,
                                            base_col_ptr, base_row_idx,
                                            base_values,
                                            compact_missing_diagonal_matcher,
                                            improve_matching,
                                            &row_perm,
                                            &matched,
                                            &exact_matching,
                                            &exact_matching_scaling,
                                            &spral_matching,
                                            &match_row_scale,
                                            &match_col_scale);
  }
#ifdef KLS_HAVE_SPRAL_SCALING
  if (use_large_spral_match &&
      !bounded_exact_match_probe_current &&
      (status != KLS_OK || 1000u * matched < 995u * solver->n)) {
    /* Greedy matching is the low-cost baseline even when much of the input
       diagonal is absent.  Escalate to exact sparse assignment only when
       that measured attempt cannot satisfy the same coverage contract; a
       missing-diagonal ratio proposes the fallback but no longer selects a
       more expensive permutation by itself.  This preserves exact-match
       robustness for difficult patterns while keeping a complete greedy
       result as an ordinary candidate, independent of matrix family. */
    free(row_perm);
    free(match_row_scale);
    free(match_col_scale);
    row_perm = NULL;
    match_row_scale = NULL;
    match_col_scale = NULL;
    matched = 0u;
    status = build_spral_hungarian_row_match_scaling(
      solver->n, solver->nnz, base_col_ptr, base_row_idx, base_values,
      &row_perm, &matched, &match_row_scale, &match_col_scale);
    exact_matching = 0;
    exact_matching_scaling = 0;
    spral_matching = status == KLS_OK;
  }
#endif
  if (kls_trace_pre_static_enabled()) {
    fprintf(stderr,
            "KLS pre-static: n=%ld weak=%ld missing=%ld large=%d "
            "prefer_unscaled=%d match_status=%d matched=%ld spral=%d\n",
            (long)solver->n, (long)weak, (long)missing_diagonal,
            (int)use_large_spral_match, prefer_unscaled_static_match,
            status, (long)matched, spral_matching);
  }
  if (status != KLS_OK || 1000u * matched < 995u * solver->n) {
#ifdef KLS_HAVE_SPRAL_SCALING
    if (bounded_exact_match_probe_current) {
      free(row_perm);
      free(match_row_scale);
      free(match_col_scale);
      row_perm = NULL;
      match_row_scale = NULL;
      match_col_scale = NULL;
      bounded_exact_match_probe_current = 0;
      bounded_exact_match_probe_rejected = 1;
      goto matching_attempt;
    }
#endif
    if (kls_trace_pre_static_enabled()) {
      fprintf(stderr, "KLS pre-static: reject match coverage\n");
    }
    goto done;
  }

  if (kls_trace_pre_static_enabled()) {
    fprintf(stderr, "KLS pre-static: match+gates done %.3fms\n",
            1e3 * (kls_now_seconds() - kls_ps_t0));
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
  if (
      getenv("KLS_STATIC_MATCH_AUTO_SCALE") != NULL &&
      solver->options.scale == KLS_SCALE_AUTO) {
    /* Matching changes which numerical entries become diagonal, so the
       original matrix's scale verdict is not evidence for the candidate's
       value frame.  This optional experiment lets the ordinary value census
       choose the first matched numeric; the default keeps the stable matched
       equilibration chosen above. */
    prefer_unscaled_static_match =
      choose_auto_scale_from_pattern(solver->n, trial_col_ptr,
                                     trial_row_idx, &solver->options,
                                     trial_values) <= 0;
  }
  if (kls_trace_pre_static_enabled()) {
    fprintf(stderr, "KLS pre-static: pattern build done %.3fms\n",
            1e3 * (kls_now_seconds() - kls_ps_t0));
  }

  kls_options trial_options = solver->options;
  /* A row permutation changes both the SCC decomposition and the quotient
     graph.  Re-run the caller's AUTO portfolio on the matched coordinates
     instead of pinning either axis to the unmatched winner.  The symbolic
     tournament, alternate numeric factor, and incumbent Pareto gates below
     supply measured evidence for the combined candidate. */

  if (getenv("KLS_PRESTATIC_FORCE_METIS") != NULL ||
      solver->medium_partial_static_metis_path) {
    trial_options.ordering = KLS_ORDERING_METIS;
  }
  if (trial_options.scale == KLS_SCALE_AUTO &&
      (1)) {
    if (match_row_scale != NULL && match_col_scale != NULL) {
      trial_row_scale = match_row_scale;
      trial_col_scale = match_col_scale;
      match_row_scale = NULL;
      match_col_scale = NULL;
    } else if ((prefer_unscaled_static_match ||
                getenv("KLS_DISABLE_LOW_WORK_MATCH_SCALE_DEFER") == NULL) &&
               isfinite(matched_column_pair_work) &&
               matched_column_pair_work <= 1.0e7 &&
               kls_repeated_update_workload(&trial_options)) {
      /* A bounded generic match is first judged in its unscaled
         coordinates.  Constructing dual matching scales before that numeric
         exists can be pure setup: a well-conditioned, linear-work result has
         no numerical or lifecycle reason to launch the optional scaled arm.
         Defer the scale plan until the actual factor supplies those signals. */
      matching_equilibration_deferred = 1;
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
#ifdef KLS_HAVE_SPRAL_SCALING

#endif

  int reused_raced_metis_q = 0;
  int bounded_matched_amd = 0;
  int matched_low_pair_ordering_selected = 0;
  if (getenv("KLS_PRESTATIC_REUSE_RACED_METIS_Q") != NULL) {
    reused_raced_metis_q = kls_prestatic_reuse_raced_metis_symbolic(
      solver, trial_col_ptr, trial_row_idx, &trial_options,
      &trial_symbolic, &trial_common, &trial_separator);
    if (reused_raced_metis_q) {
      trial_ordering = KLS_ORDERING_METIS;
      trial_score = symbolic_score(trial_symbolic);
    }
  }
  if (!reused_raced_metis_q &&
      bounded_exact_match_probe_attempted &&
      !bounded_exact_match_probe_rejected &&
      trial_options.ordering == KLS_ORDERING_AUTO &&
      kls_repeated_update_workload(&trial_options) &&
      isfinite(matched_column_pair_work) &&
      matched_column_pair_work <= 1.0e7 &&
      getenv("KLS_DISABLE_BOUNDED_MATCHED_AMD_SHORTCUT") == NULL) {
    kls_options amd_options = trial_options;
    amd_options.ordering = KLS_ORDERING_AMD;
    amd_options.use_btf = 0;
    status = analyze_with_ordering(
      solver->n, trial_col_ptr, trial_row_idx, &amd_options,
      KLS_ORDERING_AMD, &trial_symbolic, &trial_common,
      &trial_separator);
    if (status == KLS_OK && trial_symbolic != NULL) {
      trial_score = symbolic_score(trial_symbolic);
      bounded_matched_amd =
        !trial_symbolic->do_btf && trial_symbolic->nblocks == 1u &&
        trial_symbolic->maxblock == solver->n &&
        trial_symbolic->est_flops > 0.0 &&
        trial_symbolic->est_flops <= 1.0e7 &&
        trial_score > 0.0 && trial_score <= 2.0e6;
    }
    if (bounded_matched_amd) {
      trial_ordering = KLS_ORDERING_AMD;
    } else {
      if (trial_symbolic != NULL) {
        trilinos_klu_l_free_symbolic(&trial_symbolic, &trial_common);
      }
      kls_separator_analysis_clear(&trial_separator);
      trial_score = 0.0;
    }
  }
  if (bounded_exact_match_probe_current && !bounded_matched_amd) {
    /* Exact sparse assignment is a low-cost proposal only while both its
       input graph and resulting one-block AMD factor satisfy fixed work and
       storage budgets.  If the computed factor is not bounded, discard the
       whole proposal and resume with the ordinary greedy match; matrix
       dimensions or family names never choose between the permutations. */
    free(row_perm);
    free(trial_col_ptr);
    free(trial_row_idx);
    free(trial_input_to_csc);
    free(trial_values);
    free(trial_row_scale);
    free(trial_col_scale);
    free(match_row_scale);
    free(match_col_scale);
    row_perm = NULL;
    trial_col_ptr = NULL;
    trial_row_idx = NULL;
    trial_input_to_csc = NULL;
    trial_values = NULL;
    trial_row_scale = NULL;
    trial_col_scale = NULL;
    match_row_scale = NULL;
    match_col_scale = NULL;
    bounded_exact_match_probe_current = 0;
    bounded_exact_match_probe_rejected = 1;
    goto matching_attempt;
  }
  if (!reused_raced_metis_q && !bounded_matched_amd) {
    /* CKTSO-style dense-spike startup: after weighted matching, run the
       ordinary auto chooser on the matched graph.  Keeping the prestatic
       legacy context here invokes its AMD/AMF competition instead and can
       select a 20-second pivot-heavy rajat30 trial; the ordinary chooser's
       structural METIS-start route produces the compact ~9M-entry factor. */
    int ordinary_matched_ordering_ctx = 0;
#ifdef KLS_HAVE_SPRAL_SCALING
    ordinary_matched_ordering_ctx |=
      use_large_spral_match;
#endif
    if (ordinary_matched_ordering_ctx) {
      kls_prestatic_ordering_ctx = 0;
    }
    const int saved_matched_low_pair_ctx =
      kls_matched_low_pair_ordering_ctx;
    const double saved_matched_low_pair_work_ctx =
      kls_matched_low_pair_work_ctx;
    const int saved_compact_missing_match_ordering_ctx =
      kls_compact_missing_match_ordering_ctx;
    const int matched_low_pair_ordering =
      isfinite(matched_column_pair_work) &&
      matched_column_pair_work <= 1.0e7 && trial_options.use_btf &&
      trial_options.ordering == KLS_ORDERING_AUTO &&
      kls_repeated_update_workload(&trial_options);
    kls_matched_low_pair_ordering_ctx = matched_low_pair_ordering;
    kls_matched_low_pair_work_ctx = matched_column_pair_work;
    kls_compact_missing_match_ordering_ctx =
      compact_missing_diagonal_matcher;
    const int saved_value_matched_ordering_ctx =
      kls_value_matched_ordering_ctx;
    kls_value_matched_ordering_ctx = 1;
    status = choose_symbolic_for_pattern(
      solver->n, trial_col_ptr, trial_row_idx, &trial_options,
      &trial_symbolic, &trial_common, &trial_ordering, &trial_score,
      &trial_separator);
    kls_value_matched_ordering_ctx = saved_value_matched_ordering_ctx;
    matched_low_pair_ordering_selected =
      matched_low_pair_ordering && status == KLS_OK &&
      trial_symbolic != NULL;
    kls_matched_low_pair_ordering_ctx = saved_matched_low_pair_ctx;
    kls_matched_low_pair_work_ctx = saved_matched_low_pair_work_ctx;
    kls_compact_missing_match_ordering_ctx =
      saved_compact_missing_match_ordering_ctx;
    if (ordinary_matched_ordering_ctx) {
      kls_prestatic_ordering_ctx = 1;
    }
  } else {
    status = KLS_OK;
  }
  if (status != KLS_OK) {
    goto done;
  }
  if (!bounded_matched_amd &&
      !matched_low_pair_ordering_selected &&
      trial_options.ordering == KLS_ORDERING_AUTO &&
      (trial_ordering == KLS_ORDERING_AMD ||
       trial_ordering == KLS_ORDERING_AMF) &&
      trial_symbolic != NULL && trial_score > 0.0 &&
      isfinite(trial_score)) {
    /* Symbolic fill is an estimate made before numerical pivoting.  When the
       other minimum-degree ordering remains in its uncertainty band, retain
       it for one actual numeric comparison below.  Matching is already a
       recurring-work trial, so this bounded extra factor can be judged on
       the representation it would really publish instead of on matrix
       dimensions or a family label.  A bounded low-pair matched pattern has
       already compared its BTF and direct representations on the actual
       graph, so it does not repeat this ordering analysis. */
    alternate_ordering = trial_ordering == KLS_ORDERING_AMD
      ? KLS_ORDERING_AMF : KLS_ORDERING_AMD;
    kls_options alternate_options = trial_options;
    alternate_options.ordering = alternate_ordering;
    alternate_options.use_btf = trial_symbolic->do_btf ? 1 : 0;
    const int alternate_status = analyze_with_ordering(
      solver->n, trial_col_ptr, trial_row_idx, &alternate_options,
      alternate_ordering, &alternate_symbolic, &alternate_common,
      &alternate_separator);
    if (alternate_status == KLS_OK && alternate_symbolic != NULL) {
      alternate_score = symbolic_score(alternate_symbolic);
      if (!(alternate_score > 0.0) || !isfinite(alternate_score) ||
          alternate_score > 1.20 * trial_score) {
        trilinos_klu_l_free_symbolic(&alternate_symbolic,
                                     &alternate_common);
        kls_separator_analysis_clear(&alternate_separator);
        alternate_ordering = KLS_ORDERING_AUTO;
        alternate_score = DBL_MAX;
      }
    } else {
      if (alternate_symbolic != NULL) {
        trilinos_klu_l_free_symbolic(&alternate_symbolic,
                                     &alternate_common);
      }
      kls_separator_analysis_clear(&alternate_separator);
      alternate_ordering = KLS_ORDERING_AUTO;
      alternate_score = DBL_MAX;
    }
  }
  if (deferred && solver->numeric != NULL &&
      trial_symbolic != NULL) {
    /* Before constructing a deferred challenger, compare the best symbolic
       lower-cost evidence available for its orderings with the realized live
       numeric.  If every ordering is already materially worse in both fill
       and arithmetic, numerical pivoting cannot make this a sensible
       recurring candidate.  The later actual-numeric Pareto gate remains the
       backstop for optimistic symbolic estimates. */
    double best_est_fill =
      trial_symbolic->lnz + trial_symbolic->unz;
    double best_est_flops = trial_symbolic->est_flops;
    if (alternate_symbolic != NULL) {
      const double alternate_fill_est =
        alternate_symbolic->lnz + alternate_symbolic->unz;
      if (alternate_fill_est > 0.0 &&
          (!(best_est_fill > 0.0) || alternate_fill_est < best_est_fill)) {
        best_est_fill = alternate_fill_est;
      }
      if (alternate_symbolic->est_flops > 0.0 &&
          (!(best_est_flops > 0.0) ||
           alternate_symbolic->est_flops < best_est_flops)) {
        best_est_flops = alternate_symbolic->est_flops;
      }
    }
    const double incumbent_fill =
      (double)(solver->numeric->lnz + solver->numeric->unz);
    const double incumbent_flops = solver->common.flops;
    if (incumbent_fill > 0.0 && incumbent_flops > 0.0 &&
        best_est_fill >= incumbent_fill &&
        best_est_flops >= incumbent_flops) {
      /* The incumbent is already a valid numeric.  A deferred static match
         is therefore a performance challenger, and a challenger whose best
         symbolic form improves neither recurring storage nor arithmetic has
         no demonstrated payoff.  Pivot-count reduction alone is not a
         speedup signal; the quality gate remains available to inline and
         failed-factor rescue matches.  Declining here also avoids building a
         throwaway numeric.  This is a strict Pareto test over measured and
         analyzed costs, with no matrix dimensions or family detector. */
      if (kls_trace_pre_static_enabled()) {
        fprintf(stderr,
                "KLS pre-static: incumbent Pareto declines deferred "
                "candidate fill=%.3e/%.3e work=%.3e/%.3e\n",
                best_est_fill, incumbent_fill, best_est_flops,
                incumbent_flops);
      }
      goto done;
    }
    if (incumbent_fill > 0.0 && incumbent_flops > 0.0 &&
        best_est_fill > 1.20 * incumbent_fill &&
        best_est_flops > 1.20 * incumbent_flops) {
      if (kls_trace_pre_static_enabled()) {
        fprintf(stderr,
                "KLS pre-static: reject symbolically dominated deferred "
                "candidate fill=%.3e/%.3e work=%.3e/%.3e\n",
                best_est_fill, incumbent_fill, best_est_flops,
                incumbent_flops);
      }
      goto done;
    }
  }
  if (kls_trace_pre_static_enabled()) {
    fprintf(stderr, "KLS pre-static: choose done %.3fms\n",
            1e3 * (kls_now_seconds() - kls_ps_t0));
  }
#ifdef KLS_HAVE_SPRAL_SCALING
  if (use_large_spral_match && trial_ordering != KLS_ORDERING_METIS &&
      trial_score >= 8.0e7) {
    /* CKTSO-style ordering selection: compare a nested-dissection candidate
       by estimated fill before the expensive numeric trial.  Statically
       matched circuit matrices often favor the METIS+CAMD ordering by a
       large factor-work margin, and factoring once with the estimate winner
       avoids a second full trial factorization.  The estimated-fill floor
       keeps cheap trials from paying the serial nested-dissection analysis
       when the minimum-degree trial factor is already fast. */
    kls_options metis_options = trial_options;
    metis_options.ordering = KLS_ORDERING_METIS;
    trilinos_klu_l_symbolic *metis_symbolic = NULL;
    trilinos_klu_l_common metis_common;
    kls_separator_analysis metis_separator;
    int metis_est_ok = 0;
    memset(&metis_separator, 0, sizeof(metis_separator));
    if (psmetis_est_spec.launched || psmetis_est_spec.symbolic != NULL) {
      kls_psmetis_spec_join(&psmetis_est_spec);
      if (psmetis_est_spec.status == KLS_OK &&
          psmetis_est_spec.symbolic != NULL) {
        metis_symbolic = psmetis_est_spec.symbolic;
        metis_common = psmetis_est_spec.common;
        kls_separator_analysis_move(&metis_separator,
                                    &psmetis_est_spec.separator);
        psmetis_est_spec.symbolic = NULL;
        metis_est_ok = 1;
      } else {
        kls_psmetis_spec_discard(&psmetis_est_spec);
      }
    } else if (analyze_with_ordering(solver->n, trial_col_ptr,
                                     trial_row_idx, &metis_options,
                                     KLS_ORDERING_METIS, &metis_symbolic,
                                     &metis_common,
                                     &metis_separator) == KLS_OK) {
      metis_est_ok = 1;
    }
    if (metis_est_ok) {
      double metis_score = symbolic_score(metis_symbolic);
      const double metis_storage =
        metis_score * (double)(sizeof(UF_long) + sizeof(double));
      const double aggregate_cache =
        32.0 * 1024.0 * 1024.0 * (double)solver->options.threads;
      if (solver->options.threads > 2 &&
          kls_repeated_update_workload(&solver->options) &&
          isfinite(metis_storage) && metis_storage >= aggregate_cache &&
          getenv("KLS_DISABLE_GENERIC_MATCHED_ND_WIDTH_PORTFOLIO") == NULL) {
        /* Separator width is a representation/resource choice, not merely
           an analysis-time thread count.  A wide splitter minimizes NodeND
           latency, while two- and three-leaf forests can retain less factor
           traffic over a long numeric lifecycle.  Three leaves is the
           smallest additional topology not represented by the binary cut or
           the caller-width forest; it catches graphs where one coarse child
           contains most of the separator work without turning this into an
           unbounded width search.  Retain a narrower arm only for a material
           fill win whose projected lifecycle saving repays its complete graph
           comparison.  The eventual numeric fill, pivot, conditioning, and
           residual gates remain the authority; no input dimensions or graph
           family labels enter this decision. */
        const int saved_npes_override = kls_metis_ndp_npes_override;
        const double numeric_horizon =
          1.0 + (double)solver->options.expected_refactorizations;
        const double comparison_work =
          4.0 * ((double)solver->n + (double)solver->nnz);
        for (int candidate_npes = 2;
             candidate_npes <= 3 &&
             candidate_npes < solver->options.threads;
             ++candidate_npes) {
          trilinos_klu_l_symbolic *candidate_symbolic = NULL;
          trilinos_klu_l_common candidate_common;
          kls_separator_analysis candidate_separator;
          (void)trilinos_klu_l_defaults(&candidate_common);
          memset(&candidate_separator, 0, sizeof(candidate_separator));
          kls_metis_ndp_npes_override = candidate_npes;
          const int candidate_status = analyze_with_ordering(
            solver->n, trial_col_ptr, trial_row_idx, &metis_options,
            KLS_ORDERING_METIS, &candidate_symbolic, &candidate_common,
            &candidate_separator);
          if (candidate_status == KLS_OK && candidate_symbolic != NULL) {
            const double candidate_score =
              symbolic_score(candidate_symbolic);
            const double saved_fill = metis_score - candidate_score;
            const int candidate_wins =
              isfinite(candidate_score) && candidate_score > 0.0 &&
              isfinite(metis_score) && metis_score > 0.0 &&
              candidate_score <= 0.995 * metis_score &&
              saved_fill > 0.0 &&
              saved_fill * numeric_horizon >= comparison_work;
            if (kls_trace_pre_static_enabled()) {
              fprintf(stderr,
                      "KLS pre-static: matched ND width portfolio "
                      "npes=%d score=%.3e incumbent=%.3e -> %s\n",
                      candidate_npes, candidate_score, metis_score,
                      candidate_wins ? "candidate" : "incumbent");
            }
            if (candidate_wins) {
              trilinos_klu_l_free_symbolic(&metis_symbolic, &metis_common);
              kls_separator_analysis_clear(&metis_separator);
              metis_symbolic = candidate_symbolic;
              candidate_symbolic = NULL;
              metis_common = candidate_common;
              metis_score = candidate_score;
              kls_separator_analysis_move(&metis_separator,
                                          &candidate_separator);
            }
          }
          if (candidate_symbolic != NULL) {
            trilinos_klu_l_free_symbolic(&candidate_symbolic,
                                         &candidate_common);
          }
          kls_separator_analysis_clear(&candidate_separator);
        }
        kls_metis_ndp_npes_override = saved_npes_override;
      }
      if (kls_trace_pre_static_enabled()) {
        fprintf(stderr,
                "KLS pre-static: ordering scores base=%.3e metis=%.3e\n",
                trial_score, metis_score);
      }
      if (metis_score < DBL_MAX && metis_score < trial_score) {
        trilinos_klu_l_free_symbolic(&trial_symbolic, &trial_common);
        kls_separator_analysis_clear(&trial_separator);
        trial_symbolic = metis_symbolic;
        trial_common = metis_common;
        trial_ordering = KLS_ORDERING_METIS;
        trial_score = metis_score;
        kls_separator_analysis_move(&trial_separator, &metis_separator);
      } else {
        trilinos_klu_l_free_symbolic(&metis_symbolic, &metis_common);
        kls_separator_analysis_clear(&metis_separator);
      }
    }
  }
#endif

  trial_common.scale =
    trial_options.scale == KLS_SCALE_AUTO && prefer_unscaled_static_match
      ? -1
      : choose_auto_scale_from_pattern(solver->n, trial_col_ptr,
                                       trial_row_idx, &trial_options,
                                       trial_values);
  const int preapplied_matching_equilibration =
    (getenv("KLS_PRESTATIC_PREAPPLY_MATCH_SCALE") != NULL ||
     solver->medium_partial_static_metis_path) &&
    trial_options.scale == KLS_SCALE_AUTO && !prefer_unscaled_static_match &&
    trial_row_scale != NULL && trial_col_scale != NULL;
  if (preapplied_matching_equilibration) {
    /* The ordinary path first factors KLU's internally equilibrated values,
       then maybe_use_matching_equilibration factors the explicitly matched
       equilibration and keeps that second numeric when it is better.  A
       caller which already knows the matched equilibration is the desired
       candidate can build that same numeric directly.  Keep the scale arrays
       installed with the adopted pattern so refactor input preparation and
       solves retain the usual compensation. */
    apply_value_scaling(solver->n, trial_col_ptr, trial_row_idx,
                        trial_row_scale, trial_col_scale, trial_values);
    trial_common.scale = -1;
  }
  const int generic_recurring_matched_candidate =
    kls_repeated_update_workload(&trial_options) &&
    fabs(trial_options.pivot_tolerance - 0.001) <= 1.0e-12 &&
    getenv("KLS_DISABLE_GENERIC_MATCHED_TIGHT_TOLERANCE") == NULL;
  const UF_long matched_threshold_sensitive =
    generic_recurring_matched_candidate
      ? kls_matched_threshold_sensitive_columns(
          solver->n, trial_col_ptr, trial_row_idx, trial_values,
          1.0e-5, trial_options.pivot_tolerance)
      : 0u;
  const double generic_matched_tolerance_storage =
    trial_score * (double)(sizeof(UF_long) + sizeof(double));
  const double generic_matched_tolerance_cache =
    32.0 * 1024.0 * 1024.0 * (double)solver->options.threads;
  const int generic_matched_tolerance_resource_fit =
    getenv("KLS_DISABLE_GENERIC_RESOURCE_SCALED_MATCHED_TOLERANCE") !=
      NULL ||
    (isfinite(generic_matched_tolerance_storage) &&
     (generic_matched_tolerance_storage <=
        generic_matched_tolerance_cache ||
      20.0 * (double)matched_threshold_sensitive >= (double)solver->n));
  const int generic_recurring_matched_tolerance =
    generic_recurring_matched_candidate &&
    generic_matched_tolerance_resource_fit &&
    (matched_threshold_sensitive >= 2u ||
     (trial_symbolic != NULL && trial_symbolic->est_flops >= 1.0e6));
  if (getenv("KLS_TRACE_MATCHED_TOLERANCE") != NULL &&
      generic_recurring_matched_candidate) {
    fprintf(stderr,
            "KLS matched tolerance: sensitive=%lu/%lu pair_work=%.0f "
            "estimated_flops=%.0f storage=%.0f/cache=%.0f fit=%d "
            "selected=%d\n",
            (unsigned long)matched_threshold_sensitive,
            (unsigned long)solver->n,
            matched_column_pair_work,
            trial_symbolic != NULL ? trial_symbolic->est_flops : -1.0,
            generic_matched_tolerance_storage,
            generic_matched_tolerance_cache,
            generic_matched_tolerance_resource_fit,
            generic_recurring_matched_tolerance);
  }
  if (generic_recurring_matched_tolerance) {
    /* A statically matched repeated lifecycle commits to a fixed pivot
       representation.  Tightening is useful when realized matched diagonals
       actually straddle the two thresholds, or when enough estimated numeric
       work exists for even a diffuse pivot-layout improvement to amortize.
       Compact low-work factors with neither signal keep the robust caller
       tolerance: their tighter factors add startup work without changing the
       recurring walk.  Numeric acceptance and the solve residual contract
       remain authoritative.  Keep trial_options synchronized so the optional
       matching-equilibration arm compares the same tolerance. */
    trial_options.pivot_tolerance = 1.0e-5;
    trial_common.tol = 1.0e-5;
  }

  if (deferred && solver->numeric != NULL &&
      trial_symbolic != NULL) {
    const double incumbent_fill =
      (double)(solver->numeric->lnz + solver->numeric->unz);
    const double proposal_fill =
      trial_symbolic->lnz + trial_symbolic->unz;
    const double incumbent_work = solver->common.flops;
    const double proposal_work = trial_symbolic->est_flops;
    if (incumbent_fill > 0.0 && proposal_fill > 1.20 * incumbent_fill &&
        (!(incumbent_work > 0.0) || !(proposal_work > 0.0) ||
         proposal_work >= incumbent_work)) {
      /* The deferred match is a performance proposal, not a correctness
         repair: the incumbent has already factored successfully.  When its
         symbolic lower-cost proxies are both materially worse than the live
         numeric, a full trial factor cannot plausibly repay its own setup
         over the requested cycle.  Stop at the cheap symbolic evidence.
         The 20% margin tolerates estimate error and is wider than every
         fill increase accepted by numeric_candidate_is_better(). */
      if (kls_trace_pre_static_enabled()) {
        fprintf(stderr,
                "KLS pre-static: lifecycle dominance declines trial "
                "fill=%.3e/%.3e work=%.3e/%.3e\n",
                proposal_fill, incumbent_fill,
                proposal_work, incumbent_work);
      }
      goto done;
    }
  }

  if (kls_trace_pre_static_enabled()) {
    fprintf(stderr,
            "KLS pre-static: trial factor start ordering=%d scale=%ld "
            "tol=%g score=%.3e t=%.3fms\n",
            (int)trial_ordering, (long)trial_common.scale, trial_common.tol,
            trial_score, 1e3 * (kls_now_seconds() - kls_ps_t0));
  }
  int skip_trial_factor = 0;
  const int generic_bounded_parallel_first =
    !deferred &&
    getenv("KLS_DISABLE_GENERIC_BOUNDED_PRESTATIC_FIRST") == NULL &&
    solver->options.threads > 1 &&
    kls_repeated_update_workload(&solver->options) &&
    trial_symbolic != NULL && trial_symbolic->est_flops > 0.0 &&
    trial_symbolic->est_flops <= 1.0e7 &&
    trial_score > 0.0 && trial_score <= 2.0e6;
  /* A matched numeric larger than the workers' aggregate cache should not be
     discovered inside the speculative pre-static executor and then prepared
     again by the normal lifecycle machinery.  When the caller advertises a
     long update horizon that repays the match, adopt only after the same
     post-match diagonal-strength scan below, then let the normal predicted /
     pivoted first-factor path build and validate the numeric.  The thresholds
     are execution-resource and lifecycle budgets; no input dimension,
     density, or matrix-family coordinate participates. */
  const double generic_trial_storage_bytes =
    trial_score * (double)(sizeof(UF_long) + sizeof(double));
  const double generic_aggregate_cache_bytes =
    32.0 * 1024.0 * 1024.0 * (double)solver->options.threads;
  const double generic_match_lifecycle_work =
    matched_column_pair_work *
      (1.0 + (double)solver->options.expected_refactorizations);
  const int generic_resource_scaled_parallel_first =
    !deferred &&
    getenv("KLS_ENABLE_GENERIC_PRESTATIC_SKIP_TRIAL_FACTOR") != NULL &&
    solver->options.threads > 1 &&
    kls_repeated_update_workload(&solver->options) &&
    isfinite(generic_trial_storage_bytes) &&
    generic_trial_storage_bytes >= generic_aggregate_cache_bytes &&
    isfinite(generic_match_lifecycle_work) &&
    isfinite(trial_score) && trial_score > 0.0 &&
    generic_match_lifecycle_work >= 8.0 * trial_score &&
    getenv("KLS_DISABLE_GENERIC_RESOURCE_SCALED_PRESTATIC_FIRST") == NULL;
  const int generic_single_block_parallel_first =
    !deferred &&
    solver->options.threads > 1 &&
    kls_repeated_update_workload(&solver->options) &&
    trial_symbolic != NULL && trial_symbolic->nblocks == 1u &&
    trial_symbolic->maxblock == solver->n &&
    (trial_symbolic->structural_rank == KLS_KLU_EMPTY ||
     trial_symbolic->structural_rank == solver->n) &&
    getenv("KLS_DISABLE_GENERIC_SINGLE_BLOCK_PRESTATIC_FIRST") == NULL;
  const int skip_trial_factor_requested =
    !deferred &&
    (generic_bounded_parallel_first ||
     generic_resource_scaled_parallel_first ||
     generic_single_block_parallel_first ||
     getenv("KLS_PRESTATIC_SKIP_TRIAL_FACTOR") != NULL);
  if (skip_trial_factor_requested &&
      trial_options.scale == KLS_SCALE_AUTO && !prefer_unscaled_static_match &&
      trial_row_scale != NULL && trial_col_scale != NULL) {
    /* The trial path applies the matching equilibration inside its own
       numeric comparison; the skip path must scale the adopted values
       directly so the static diagonal the predicted fill pivots on is
       the equilibrated (dominant) one.  Solves and later input
       preparation compensate through the adopted scale arrays. */
    apply_value_scaling(solver->n, trial_col_ptr, trial_row_idx,
                        trial_row_scale, trial_col_scale, trial_values);
  }
  if (skip_trial_factor_requested) {
    if (prefer_unscaled_static_match) {
      /* The generic unscaled baseline has not factored the matching
         equilibration.  Do not publish unverified scale arrays alongside
         the unscaled values; the retained permutation alone is the candidate
         proven by the post-match diagonal check below. */
      free(trial_row_scale);
      free(trial_col_scale);
      trial_row_scale = NULL;
      trial_col_scale = NULL;
    }
    /* Experimental: when the matched, scaled diagonal is statically
       strong (the block trial's criterion), adopt without the serial
       pivoted trial factor and let the predicted parallel first factor
       build the numeric on the normal path. */
    /* Same tolerance the trial-factor quality gate applies to actual
       pivot deviations: columns whose matched diagonal is still weak
       are the ones at risk of deviating, and the predicted path's
       rescue + probe + serial fallback police real failures. */
    const UF_long weak_cap = weak / 20u + 16u;
    UF_long post_weak = 0;
    for (UF_long k = 0; k < solver->n && post_weak <= weak_cap; ++k) {
      double colmax = 0.0;
      double diag = -1.0;
      for (UF_long p = trial_col_ptr[k]; p < trial_col_ptr[k + 1u]; ++p) {
        const double av = fabs(trial_values[p]);
        colmax = colmax < av ? av : colmax;
        if (trial_row_idx[p] == k) {
          diag = av;
        }
      }
      if (diag < 0.0 || diag < 1.0e-4 * colmax) {
        post_weak++;
      }
    }
    skip_trial_factor = post_weak <= weak_cap;
    if (skip_trial_factor && post_weak == 0u &&
        generic_recurring_matched_tolerance && trial_common.scale == -1) {
      /* A complete match whose every installed diagonal clears the same
         relative-strength test used to skip the trial factor has direct
         numerical evidence for a tighter static-pivot threshold.  Retain
         more of that permutation at 1e-5; the ordinary solve residual
         contract remains the accuracy backstop. */
      trial_common.tol = 1.0e-5;
    }
    if (kls_trace_pre_static_enabled()) {
      fprintf(stderr,
              "KLS pre-static: post-match weak=%ld cap=%ld skip=%d\n",
              (long)post_weak, (long)weak_cap, skip_trial_factor);
    }
  }
  if (!skip_trial_factor) {
  {
    /* the prestatic trial classes carry the columns the pipelined
       kernel wants: large-spral (pre2: init 57.4 -> 31.6s) and the
       medium matched class (onetone1: 12% init win, rajat25 parity);
       est_flops is EMPTY for given orderings so the class flags are
       the routing signal.  mac_econ-like light-column matrices never
       reach this trial (pattern classifiers). */
    int pipe_route = medium_weak_candidate &&
                     solver->nnz >= 8u * solver->n;
    /* density floor: at avg degree ~7 (rajat25) the trial's columns
       are light enough that the pipeline's round machinery costs ~7%;
       the winners measure 9.5+ (onetone1, twotone) */
#ifdef KLS_HAVE_SPRAL_SCALING
    /* the same light-column density floor as the medium class: the
       pipeline's round machinery is pathological at low entries/column
       (rajat29, nnz/n=5.8: piped trial ground >16s where the serial
       factor takes 0.11s; pre2 at 8.85 keeps the route) */
    pipe_route = pipe_route ||
                 (use_large_spral_match && solver->nnz >= 8u * solver->n);
#endif

    if (pipe_route && solver->options.threads >= 2 &&
        getenv("KLS_DISABLE_PIPE_ROUTE") == NULL) {
      kls_klu_pipe_threads =
        kls_pipe_scale_threads(solver->options.threads,
                               trial_symbolic != NULL
                                 ? trial_symbolic->est_flops : 0.0);
    }
    /* Speculative METIS refinement: the refinement's serial
       NodeND+CAMD dominates its cost and depends only on the trial
       pattern and the scale/tol just decided, so start it on a worker
       before the trial factor.  When the main-side trial is serial
       (no pipe routing) the worker also runs the METIS trial factor
       and its diagnostics - the whole refinement then joins ready.
       The join inside the refinement applies the flops-aware half of
       the gate; a losing speculation is joined and discarded there
       (or at done: on trial rejection). */

  }
  if (getenv("KLS_STATIC_PERTURB") != NULL) {
    /* hard-static regime: keep the matched diagonal everywhere and
       perturb tiny pivots (SuperLU_DIST recipe); refinement recovers
       the contract.  sqrt(eps) x amax bounds growth at ~1e8. */
    double amax = 0.0;
    for (UF_long pz = 0; pz < (UF_long)solver->nnz; ++pz) {
      const double av = fabs(trial_values[pz]);
      if (av > amax) amax = av;
    }
    trial_common.kls_static_perturb = 1.49e-8 * amax;
  }
  if (!deferred && getenv("KLS_PRESTATIC_SKIP_TRIAL") != NULL) {
    /* probe: adopt the matched pattern UNFACTORED and let the main
       path's parallel first-factor machinery build the numeric on
       these coordinates (the serial trial is pre2's largest cost) */
    trial_numeric = NULL;
    goto kls_adopt_unfactored;
  }
  const double trial_est_fill = trial_symbolic != NULL
    ? trial_symbolic->lnz + trial_symbolic->unz : 0.0;
  const double trial_est_flops = trial_symbolic != NULL
    ? trial_symbolic->est_flops : 0.0;
  const int generic_dense_matched_pipe_factor =
    solver->options.threads > 1 &&
    kls_repeated_update_workload(&solver->options) &&
    trial_est_flops >= 3.0e8 &&
    trial_est_flops >= 8192.0 * (double)solver->n;
  if (kls_trace_pre_static_enabled()) {
    fprintf(stderr,
            "KLS pre-static: moderate-fragmented pipe=%d"
            " trial-est=%.3e trial-fill=%.3e\n",
            0,
            trial_est_flops, trial_est_fill);
  }
  const int generic_unmatched_colamd_proposal =
    !deferred && !forced_match &&
    getenv("KLS_DISABLE_GENERIC_UNMATCHED_COLAMD_CHALLENGER") == NULL &&
    solver->orientation == KLS_ORIENTATION_NORMAL &&
    solver->input_format == KLS_INPUT_CSC &&
    solver->input_to_csc == NULL && solver->values == NULL &&
    solver->row_perm == NULL && solver->user_col_perm == NULL &&
    solver->row_scale == NULL && solver->col_scale == NULL &&
    solver->block_order_perm == NULL &&
    solver->options.orientation == KLS_ORIENTATION_AUTO &&
    solver->options.ordering == KLS_ORDERING_AUTO &&
    solver->options.scale == KLS_SCALE_AUTO && solver->options.use_btf &&
    solver->options.backend != KLS_BACKEND_SERIAL &&
    solver->options.threads >= 2 &&
    kls_repeated_update_workload(&solver->options) &&
    trial_symbolic != NULL && !trial_symbolic->do_btf &&
    trial_symbolic->nblocks == 1u &&
    trial_symbolic->maxblock == solver->n &&
    trial_est_fill > 0.0 && trial_est_fill <= 3000000.0 &&
    trial_est_flops >= 64.0 * trial_est_fill &&
    isfinite(matched_column_pair_work) &&
    matched_column_pair_work <= 10000000.0;
  if (generic_unmatched_colamd_proposal) {
    unmatched_colamd_job.n = solver->n;
    unmatched_colamd_job.col_ptr = base_col_ptr;
    unmatched_colamd_job.row_idx = base_row_idx;
    unmatched_colamd_job.values = base_values;
    unmatched_colamd_job.options = solver->options;
    unmatched_colamd_job.options.ordering = KLS_ORDERING_COLAMD;
    unmatched_colamd_job.options.orientation = KLS_ORIENTATION_NORMAL;
    unmatched_colamd_job.options.use_btf = 0;
    unmatched_colamd_job.matched_estimated_fill = trial_est_fill;
    unmatched_colamd_launched =
      pthread_create(&unmatched_colamd_thread, NULL,
                     kls_unmatched_colamd_main,
                     &unmatched_colamd_job) == 0;
  }
  if (generic_dense_matched_pipe_factor) {
    /* This is already a fully matched, symbolically ordered candidate, so
       total predicted work and work density are the relevant pipeline
       economics.  Scale the crew to that work instead of applying the broad
       5e9-flop cold-factor floor: dense matched trials can repay a small crew
       well below it, while low-work-per-row giant trials remain serial. */
    kls_klu_pipe_threads =
      kls_pipe_scale_threads(solver->options.threads, trial_est_flops);
  } else if (kls_klu_pipe_threads == 0) {
    /* the class pipe route above already decided for the matched
       medium/large classes (est_flops is empty for given orderings, so
       the fill-proxy gate here would wrongly unroute the medium class:
       twotone's trial measured 0.49s serial vs 0.20s piped); only fill
       in when no route fired.  Light-column non-lowdeg trials stay
       serial regardless of the est: the pipe's round machinery is
       pathological there (rajat30's probe-routed trial ground 21.6s
       piped vs sub-second serial; mac_econ - the pipe's legitimate
       light-column user - is lowdeg-classified) */
    if (getenv("KLS_KLU_PIPE_FORCE") != NULL ||
        solver->col_ptr[solver->n] >= 8 * solver->n) {
      kls_klu_pipe_threads =
        kls_pipe_first_factor_threads(solver, trial_symbolic);
    }
  }
  kls_klu_pipe_det = 0;
  const double trial_numeric_start = kls_now_seconds();
  trial_numeric =
    trilinos_klu_l_factor(trial_col_ptr, trial_row_idx, trial_values,
                          trial_symbolic, &trial_common);
  trial_factor_seconds = kls_now_seconds() - trial_numeric_start;
  kls_klu_pipe_threads = 0;
  kls_klu_pipe_det = 0;
  kls_klu_pipe_threads = 0;

  int trial_numeric_ok =
    trial_numeric != NULL && trial_common.status >= 0 &&
    trial_common.status != TRILINOS_KLU_SINGULAR;
  if (trial_numeric_ok) {
    (void)trilinos_klu_l_flops(trial_symbolic, trial_numeric,
                               &trial_common);
    (void)trilinos_klu_l_rcond(trial_symbolic, trial_numeric,
                               &trial_common);
  }
  if (trial_numeric_ok && generic_recurring_matched_tolerance &&
      (!(trial_common.rcond > 0.0) ||
       trial_common.rcond < 32.0 * DBL_EPSILON)) {
    /* A tighter matched pivot layout is a performance candidate, not a
       license to retain a factor whose reciprocal condition estimate has
       fallen into the rounding-noise band.  Restore the caller's threshold
       before the optional matching-equilibration comparison, so that arm
       also evaluates the same robust tolerance and can reconstruct its
       ordinary scaled representation. */
    trilinos_klu_l_free_numeric(&trial_numeric, &trial_common);
    trial_options.pivot_tolerance = solver->options.pivot_tolerance;
    trial_common.tol = solver->options.pivot_tolerance;
    trial_common.status = TRILINOS_KLU_OK;
    trial_common.numerical_rank = KLS_KLU_EMPTY;
    trial_common.singular_col = KLS_KLU_EMPTY;
    const double robust_start = kls_now_seconds();
    trial_numeric = trilinos_klu_l_factor(
      trial_col_ptr, trial_row_idx, trial_values,
      trial_symbolic, &trial_common);
    trial_factor_seconds += kls_now_seconds() - robust_start;
    trial_numeric_ok = trial_numeric != NULL &&
      trial_common.status >= TRILINOS_KLU_OK &&
      trial_common.status != TRILINOS_KLU_SINGULAR;
    if (trial_numeric_ok) {
      (void)trilinos_klu_l_flops(trial_symbolic, trial_numeric,
                                 &trial_common);
      (void)trilinos_klu_l_rcond(trial_symbolic, trial_numeric,
                                 &trial_common);
    }
    if (kls_trace_pre_static_enabled()) {
      fprintf(stderr,
              "KLS pre-static: tight matched rcond fallback tol=%.3g "
              "valid=%d rcond=%.3e\n",
              trial_common.tol, trial_numeric_ok,
              trial_numeric_ok ? trial_common.rcond : -1.0);
    }
  }
  if (alternate_symbolic != NULL) {
    alternate_common.scale = trial_common.scale;
    alternate_common.tol = trial_common.tol;
    alternate_common.halt_if_singular = trial_common.halt_if_singular;
    kls_klu_pipe_threads =
      kls_pipe_first_factor_threads(solver, alternate_symbolic);
    kls_klu_pipe_det = 0;
    const double alternate_numeric_start = kls_now_seconds();
    alternate_numeric = trilinos_klu_l_factor(
      trial_col_ptr, trial_row_idx, trial_values, alternate_symbolic,
      &alternate_common);
    alternate_factor_seconds =
      kls_now_seconds() - alternate_numeric_start;
    kls_klu_pipe_threads = 0;
    kls_klu_pipe_det = 0;
    const int alternate_numeric_ok =
      alternate_numeric != NULL && alternate_common.status >= 0 &&
      alternate_common.status != TRILINOS_KLU_SINGULAR;
    if (alternate_numeric_ok) {
      (void)trilinos_klu_l_flops(alternate_symbolic, alternate_numeric,
                                 &alternate_common);
      (void)trilinos_klu_l_rcond(alternate_symbolic, alternate_numeric,
                                 &alternate_common);
    }
    int select_alternate = !trial_numeric_ok && alternate_numeric_ok;
    if (trial_numeric_ok && alternate_numeric_ok &&
        !(trial_common.rcond > 0.0 && alternate_common.rcond > 0.0 &&
          alternate_common.rcond < 0.01 * trial_common.rcond)) {
      const double trial_fill =
        (double)(trial_numeric->lnz + trial_numeric->unz);
      const double alternate_fill =
        (double)(alternate_numeric->lnz + alternate_numeric->unz);
      const int pareto_better =
        alternate_fill <= 0.995 * trial_fill &&
        alternate_common.noffdiag <= trial_common.noffdiag &&
        (!(trial_common.flops > 0.0) ||
         (alternate_common.flops > 0.0 &&
          alternate_common.flops <= trial_common.flops));
      select_alternate = pareto_better ||
        numeric_candidate_is_better(&trial_common, trial_numeric,
                                    &alternate_common,
                                    alternate_numeric);
    }
    if (kls_trace_pre_static_enabled()) {
      fprintf(stderr,
              "KLS pre-static: numeric ordering trial primary=%d"
              " fill=%ld work=%.3e alternate=%d fill=%ld work=%.3e"
              " select=%d\n",
              (int)trial_ordering,
              trial_numeric_ok
                ? (long)(trial_numeric->lnz + trial_numeric->unz) : -1L,
              trial_numeric_ok ? trial_common.flops : -1.0,
              (int)alternate_ordering,
              alternate_numeric_ok
                ? (long)(alternate_numeric->lnz + alternate_numeric->unz)
                : -1L,
              alternate_numeric_ok ? alternate_common.flops : -1.0,
              select_alternate);
    }
    if (select_alternate) {
      if (trial_numeric != NULL) {
        trilinos_klu_l_free_numeric(&trial_numeric, &trial_common);
      }
      trilinos_klu_l_free_symbolic(&trial_symbolic, &trial_common);
      kls_separator_analysis_clear(&trial_separator);
      trial_symbolic = alternate_symbolic;
      trial_numeric = alternate_numeric;
      trial_common = alternate_common;
      trial_ordering = alternate_ordering;
      trial_score = alternate_score;
      trial_factor_seconds = alternate_factor_seconds;
      kls_separator_analysis_move(&trial_separator,
                                  &alternate_separator);
      alternate_symbolic = NULL;
      alternate_numeric = NULL;
      trial_numeric_ok = 1;
    } else {
      if (alternate_numeric != NULL) {
        trilinos_klu_l_free_numeric(&alternate_numeric,
                                    &alternate_common);
      }
      trilinos_klu_l_free_symbolic(&alternate_symbolic,
                                   &alternate_common);
      kls_separator_analysis_clear(&alternate_separator);
    }
  }
  if (unmatched_colamd_launched) {
    pthread_join(unmatched_colamd_thread, NULL);
    unmatched_colamd_joined = 1;
    const int unmatched_numeric_ok =
      unmatched_colamd_job.status == KLS_OK &&
      unmatched_colamd_job.symbolic != NULL &&
      unmatched_colamd_job.numeric != NULL &&
      unmatched_colamd_job.common.status >= TRILINOS_KLU_OK &&
      unmatched_colamd_job.common.status != TRILINOS_KLU_SINGULAR &&
      unmatched_colamd_job.common.flops > 0.0 &&
      isfinite(unmatched_colamd_job.common.rcond) &&
      unmatched_colamd_job.common.rcond > 0.0;
    int unmatched_residual_ok = 0;
    if (unmatched_numeric_ok) {
      trilinos_klu_l_symbolic *saved_symbolic = solver->symbolic;
      trilinos_klu_l_numeric *saved_numeric = solver->numeric;
      const trilinos_klu_l_common saved_common = solver->common;
      solver->symbolic = unmatched_colamd_job.symbolic;
      solver->numeric = unmatched_colamd_job.numeric;
      solver->common = unmatched_colamd_job.common;
      unmatched_residual_ok = kls_direct_klu_numeric_residual_probe(
        solver, base_values);
      unmatched_colamd_job.common = solver->common;
      solver->symbolic = saved_symbolic;
      solver->numeric = saved_numeric;
      solver->common = saved_common;
    }
    const double matched_numeric_fill = trial_numeric_ok
      ? (double)trial_numeric->lnz + (double)trial_numeric->unz
      : DBL_MAX;
    const double unmatched_numeric_fill = unmatched_numeric_ok
      ? (double)unmatched_colamd_job.numeric->lnz +
          (double)unmatched_colamd_job.numeric->unz
      : DBL_MAX;
    const int unmatched_pareto =
      unmatched_numeric_ok && unmatched_residual_ok && trial_numeric_ok &&
      unmatched_numeric_fill <= 0.98 * matched_numeric_fill &&
      unmatched_colamd_job.common.flops <= 0.80 * trial_common.flops &&
      !(trial_common.rcond > 0.0 &&
        unmatched_colamd_job.common.rcond < 0.01 * trial_common.rcond);
    if (kls_trace_pre_static_enabled()) {
      fprintf(stderr,
              "KLS pre-static: unmatched COLAMD fill=%.4e/%.4e "
              "work=%.4e/%.4e pivots=%ld/%ld rcond=%.3e/%.3e "
              "residual=%d -> %s\n",
              unmatched_numeric_fill, matched_numeric_fill,
              unmatched_numeric_ok
                ? unmatched_colamd_job.common.flops : -1.0,
              trial_numeric_ok ? trial_common.flops : -1.0,
              unmatched_numeric_ok
                ? (long)unmatched_colamd_job.common.noffdiag : -1L,
              trial_numeric_ok ? (long)trial_common.noffdiag : -1L,
              unmatched_numeric_ok
                ? unmatched_colamd_job.common.rcond : -1.0,
              trial_numeric_ok ? trial_common.rcond : -1.0,
              unmatched_residual_ok,
              unmatched_pareto ? "candidate" : "reject");
    }
    UF_long *unmatched_col_ptr = NULL;
    UF_long *unmatched_row_idx = NULL;
    double *unmatched_values = NULL;
    if (unmatched_pareto &&
        solver->n < (UF_long)(SIZE_MAX / sizeof(*unmatched_col_ptr)) &&
        solver->nnz <= (UF_long)(SIZE_MAX / sizeof(*unmatched_row_idx)) &&
        solver->nnz <= (UF_long)(SIZE_MAX / sizeof(*unmatched_values))) {
      unmatched_col_ptr = (UF_long *)malloc(
        ((size_t)solver->n + 1u) * sizeof(*unmatched_col_ptr));
      unmatched_row_idx = (UF_long *)malloc(
        (size_t)solver->nnz * sizeof(*unmatched_row_idx));
      unmatched_values = (double *)malloc(
        (size_t)solver->nnz * sizeof(*unmatched_values));
    }
    if (unmatched_col_ptr != NULL && unmatched_row_idx != NULL &&
        unmatched_values != NULL) {
      memcpy(unmatched_col_ptr, base_col_ptr,
             ((size_t)solver->n + 1u) * sizeof(*unmatched_col_ptr));
      memcpy(unmatched_row_idx, base_row_idx,
             (size_t)solver->nnz * sizeof(*unmatched_row_idx));
      memcpy(unmatched_values, base_values,
             (size_t)solver->nnz * sizeof(*unmatched_values));
      trilinos_klu_l_free_numeric(&trial_numeric, &trial_common);
      trilinos_klu_l_free_symbolic(&trial_symbolic, &trial_common);
      kls_separator_analysis_clear(&trial_separator);
      free(row_perm);
      free(trial_col_ptr);
      free(trial_row_idx);
      free(trial_input_to_csc);
      free(trial_values);
      free(trial_row_scale);
      free(trial_col_scale);
      free(match_row_scale);
      free(match_col_scale);
      row_perm = NULL;
      trial_col_ptr = unmatched_col_ptr;
      trial_row_idx = unmatched_row_idx;
      trial_input_to_csc = NULL;
      trial_values = unmatched_values;
      trial_row_scale = NULL;
      trial_col_scale = NULL;
      match_row_scale = NULL;
      match_col_scale = NULL;
      trial_symbolic = unmatched_colamd_job.symbolic;
      trial_numeric = unmatched_colamd_job.numeric;
      trial_common = unmatched_colamd_job.common;
      trial_ordering = KLS_ORDERING_COLAMD;
      trial_score = unmatched_colamd_job.score;
      trial_factor_seconds = unmatched_colamd_job.factor_seconds;
      kls_separator_analysis_move(&trial_separator,
                                  &unmatched_colamd_job.separator);
      unmatched_colamd_job.symbolic = NULL;
      unmatched_colamd_job.numeric = NULL;
      exact_matching = 0;
      exact_matching_scaling = 0;
      spral_matching = 0;
      solver->nearly_missing_diagonal_early_match_selected = 0;
      selected_unmatched_colamd = 1;
      /* This original-frame ordering has already won a complete numeric
         portfolio against the matched representation: both factors were
         built, residual-policed, and compared on realized fill/work.  Close
         the ordering portfolio for this factor epoch instead of deferring a
         second, unbounded METIS analysis+factor to the first changed input.
         Failure recovery still clears/rebuilds the factor through the
         ordinary paths; this only suppresses redundant speculation. */
      solver->auto_metis_checked = 1;
      trial_numeric_ok = 1;
    } else {
      free(unmatched_col_ptr);
      free(unmatched_row_idx);
      free(unmatched_values);
    }
    kls_unmatched_colamd_discard(&unmatched_colamd_job);
  }
  if (kls_trace_pre_static_enabled()) {
    fprintf(stderr, "KLS pre-static: stage klu_factor done %.3fms"
            " perturbed=%ld\n",
            1e3 * (kls_now_seconds() - kls_ps_t0),
            (long)trial_common.kls_perturb_count);
  }
  if (!trial_numeric_ok) {
    if (kls_trace_pre_static_enabled()) {
      fprintf(stderr, "KLS pre-static: reject trial factor status=%ld\n",
              (long)trial_common.status);
    }
    goto done;
  }

  if (kls_trace_pre_static_enabled()) {
    fprintf(stderr, "KLS pre-static: stage flops+rcond done %.3fms\n",
            1e3 * (kls_now_seconds() - kls_ps_t0));
  }

  if (kls_trace_pre_static_enabled()) {
    fprintf(stderr, "KLS pre-static: stage metis_refine done %.3fms\n",
            1e3 * (kls_now_seconds() - kls_ps_t0));
  }
  if (matching_equilibration_deferred) {
    const int low_work_safe_unscaled =
      getenv("KLS_DISABLE_LOW_WORK_MATCH_SCALE_DEFER") == NULL &&
      trial_common.rcond >= sqrt(DBL_EPSILON) &&
      trial_symbolic != NULL && trial_symbolic->est_flops > 0.0 &&
      trial_symbolic->est_flops < 1.0e6 &&
      matched_threshold_sensitive < 2u;
    if (!low_work_safe_unscaled) {
      status = build_matching_equilibration(
        solver->n, base_col_ptr, base_row_idx, base_values, row_perm,
        &trial_row_scale, &trial_col_scale);
      if (status != KLS_OK) {
        goto done;
      }
    } else if (kls_trace_pre_static_enabled()) {
      fprintf(stderr,
              "KLS pre-static: skip low-work matched scale plan "
              "rcond=%.3e work=%.3e sensitive=%ld\n",
              trial_common.rcond, trial_symbolic->est_flops,
              (long)matched_threshold_sensitive);
    }
    matching_equilibration_deferred = 0;
  }
  int retain_resource_scaled_native_match = 0;
  if (!deferred &&
      !preapplied_matching_equilibration && !selected_unmatched_colamd &&
      trial_row_scale != NULL && trial_col_scale != NULL &&
      trial_numeric_ok && trial_common.scale <= 0 &&
      isfinite(trial_common.rcond) && trial_common.rcond > 0.0 &&
      trial_common.kls_perturb_count == 0u &&
      kls_repeated_update_workload(&solver->options) &&
      getenv("KLS_DISABLE_GENERIC_RESOURCE_SCALED_NATIVE_MATCH_CERTIFICATE") ==
        NULL) {
    const double native_fill =
      (double)trial_numeric->lnz + (double)trial_numeric->unz;
    const double native_storage =
      native_fill * (double)(sizeof(UF_long) + sizeof(double));
    const double aggregate_cache =
      32.0 * 1024.0 * 1024.0 * (double)solver->options.threads;
    if (native_fill > 0.0 && native_fill <= trial_score &&
        isfinite(native_storage) && native_storage >= aggregate_cache) {
      /* A very large inline matched numeric can make a speculative second
         factor much more expensive than the match itself.  When the native
         numeric already improves on its symbolic storage bound, certify the
         actual value frame with the same deterministic backward-residual
         probe used by representation tournaments.  A clean result makes
         explicit equilibration an unsupported performance speculation, so
         retain the exact native numeric.  Smaller factors and numerics that
         miss either realized-economics or residual evidence keep the full
         scale portfolio.  This admission uses only measured factor state,
         lifecycle resources, and numerical validity. */
      UF_long *saved_col_ptr = solver->col_ptr;
      UF_long *saved_row_idx = solver->row_idx;
      trilinos_klu_l_symbolic *saved_symbolic = solver->symbolic;
      trilinos_klu_l_numeric *saved_numeric = solver->numeric;
      const trilinos_klu_l_common saved_common = solver->common;
      solver->col_ptr = trial_col_ptr;
      solver->row_idx = trial_row_idx;
      solver->symbolic = trial_symbolic;
      solver->numeric = trial_numeric;
      solver->common = trial_common;
      retain_resource_scaled_native_match =
        kls_direct_klu_numeric_residual_probe(solver, trial_values);
      trial_common = solver->common;
      solver->col_ptr = saved_col_ptr;
      solver->row_idx = saved_row_idx;
      solver->symbolic = saved_symbolic;
      solver->numeric = saved_numeric;
      solver->common = saved_common;
      if (kls_trace_pre_static_enabled()) {
        fprintf(stderr,
                "KLS pre-static: resource-scaled native match "
                "fill=%.3e/%.3e storage=%.3e/cache=%.3e residual=%d\n",
                native_fill, trial_score, native_storage, aggregate_cache,
                retain_resource_scaled_native_match);
      }
    }
  }
  if (retain_resource_scaled_native_match) {
    free(trial_row_scale);
    free(trial_col_scale);
    trial_row_scale = NULL;
    trial_col_scale = NULL;
  }
  if (!preapplied_matching_equilibration &&
             !selected_unmatched_colamd) {
    maybe_use_matching_equilibration(solver->n, solver->nnz, trial_col_ptr,
                                     trial_row_idx, &trial_options,
                                     trial_symbolic, &trial_values,
                                     &trial_row_scale, &trial_col_scale,
                                     &trial_numeric, &trial_common);
  }
  if (kls_trace_pre_static_enabled()) {
    fprintf(stderr, "KLS pre-static: stages total-so-far %.3fms\n",
            1e3 * (kls_now_seconds() - kls_ps_t0));
    fprintf(stderr,
            "KLS pre-static: trial numeric lnz=%ld unz=%ld noffdiag=%ld "
            "rcond=%.3e flops=%.3e\n",
            (long)trial_numeric->lnz, (long)trial_numeric->unz,
            (long)trial_common.noffdiag, trial_common.rcond,
            trial_common.flops);
  }
  if ((!selected_unmatched_colamd &&
       trial_common.noffdiag > weak / 20u + 16u) ||
      trial_common.rcond <= 0.0) {
    if (kls_trace_pre_static_enabled()) {
      fprintf(stderr, "KLS pre-static: reject quality gate\n");
    }
    goto done;
  }
  if (deferred && solver->numeric != NULL) {
    /* A deferred match is a recurring-work challenger to a numeric that has
       already factored successfully.  Never publish a challenger that is
       strongly Pareto-dominated by that incumbent in both realized storage
       and arithmetic.  This is deliberately based on the two actual
       numerics, not on matrix dimensions or a family/profile label. */
    const double incumbent_fill =
      (double)(solver->numeric->lnz + solver->numeric->unz);
    const double challenger_fill =
      (double)(trial_numeric->lnz + trial_numeric->unz);
    const double incumbent_flops = solver->common.flops;
    const double challenger_flops = trial_common.flops;
    const int dominated_fill =
      incumbent_fill > 0.0 && challenger_fill > 2.0 * incumbent_fill;
    const int dominated_work =
      incumbent_flops > 0.0 && challenger_flops > 0.0 &&
      challenger_flops > 2.0 * incumbent_flops;
    if (dominated_fill && dominated_work) {
      if (kls_trace_pre_static_enabled()) {
        fprintf(stderr,
                "KLS pre-static: reject dominated deferred numeric "
                "fill=%.3e/%.3e work=%.3e/%.3e\n",
                challenger_fill, incumbent_fill, challenger_flops,
                incumbent_flops);
      }
      goto done;
    }
  }
  }

kls_adopt_unfactored:;
  /* The early METIS race reads the pattern arrays this adoption frees:
     join and discard it first. */
  kls_metis_race_abandon(solver);
  trilinos_klu_l_symbolic *old_symbolic = solver->symbolic;
  trilinos_klu_l_common old_common = solver->common;
  if (solver->numeric != NULL) {
    /* deferred adoption replaces a live numeric */
    trilinos_klu_l_free_numeric(&solver->numeric, &old_common);
  }
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
  solver->compact_missing_diagonal_match_selected =
    compact_missing_diagonal_matcher && !selected_unmatched_colamd;
  solver->symmetric_partial_diagonal_match_selected = 0;
  solver->exact_matching_selected = exact_matching;
  solver->exact_matching_scaling_selected =
    exact_matching_scaling && trial_row_scale != NULL &&
    trial_col_scale != NULL;
  solver->spral_matching_selected = spral_matching;
  solver->symbolic = trial_symbolic;
  solver->prestatic_reused_raced_metis_symbolic = reused_raced_metis_q;
  solver->prestatic_dense_spiked_match = 0;
  kls_separator_analysis_move(&solver->separator, &trial_separator);
  kls_invalidate_factor_etree_stats(solver);
  solver->numeric = trial_numeric;
  kls_numeric_replaced_invalidate(solver);
  /* Numeric replacement clears all value-derived policy state.  Republish
     the input proposal and the representation actually adopted by this same
     replacement after that invalidation, so preparation and public stats
     observe the new factor rather than the discarded incumbent. */
  solver->nearly_missing_diagonal_early_match_selected =
    nearly_missing_diagonal_early_match;
  solver->compact_missing_diagonal_match_candidate =
    compact_missing_diagonal_matcher;
  solver->compact_missing_diagonal_match_selected =
    compact_missing_diagonal_matcher && !selected_unmatched_colamd;
  if (solver->numeric != NULL && trial_factor_seconds > 0.0) {
    solver->numeric_full_factor_seconds = trial_factor_seconds;
  }
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
  solver->prestatic_adopted_unfactored = solver->numeric == NULL;
  accepted = 1;

done:
  if (unmatched_colamd_launched && !unmatched_colamd_joined) {
    pthread_join(unmatched_colamd_thread, NULL);
    unmatched_colamd_joined = 1;
  }
  kls_unmatched_colamd_discard(&unmatched_colamd_job);
  kls_medium_partial_static_metis_ctx = 0;
  kls_low_work_symmetric_partial_diagonal_pts_ctx = 0;
  kls_large_weak_diagonal_static_metis_ctx = 0;
  kls_prestatic_ordering_ctx = 0;
  /* the speculative refinement/estimate analyses read the trial
     pattern arrays freed below: join and discard any pending runs */
  kls_psmetis_spec_discard(&psmetis_spec);
  kls_psmetis_spec_discard(&psmetis_est_spec);
  if (!accepted) {

    if (trial_numeric != NULL) {
      trilinos_klu_l_free_numeric(&trial_numeric, &trial_common);
    }
    if (trial_symbolic != NULL) {
      trilinos_klu_l_free_symbolic(&trial_symbolic, &trial_common);
    }
    kls_separator_analysis_clear(&trial_separator);
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
  if (alternate_numeric != NULL) {
    trilinos_klu_l_free_numeric(&alternate_numeric, &alternate_common);
  }
  if (alternate_symbolic != NULL) {
    trilinos_klu_l_free_symbolic(&alternate_symbolic, &alternate_common);
  }
  kls_separator_analysis_clear(&alternate_separator);
  free(owned_col_ptr);
  free(owned_row_idx);
  free(owned_input_to_csc);
  free(owned_values);
  *elapsed += kls_now_seconds() - start;
}

/* Factor-time METIS promotion race (hard core D, brick 2): the serial
   METIS analyze + trial factorization that maybe_promote_auto_metis runs
   on giant-class matrices costs 10-13s and executes after the predicted
   first factorization has already spent 8s building a numeric the
   promotion then throws away.  Both sides only read the immutable
   pattern arrays and a private copy of the values, so the trial runs on
   a background thread launched before the first numeric is built and is
   joined at the promotion point, hiding the shorter side entirely. */
struct kls_metis_race_s {
  pthread_t thread;
  int active;
  UF_long n;
  UF_long *col_ptr;
  UF_long *row_idx;
  double *values_copy;
  kls_options options;
  int options_base_scale;
  const trilinos_klu_l_symbolic *scale_symbolic; /* solver-owned */
  int metis_wanted;
  /* stage handshake: the METIS analyze needs only the pattern and can
     start at analyze time; the scale trials and the trial factor wait
     for the values (0 wait, 1 go, 2 abort) */
  _Atomic int stage2;
  _Atomic int finished;   /* worker set, last store before thread exit */
  _Atomic int analyze_done; /* worker set right after the METIS analyze */
  int values_signaled;
  int analyze_status;
  trilinos_klu_l_symbolic *symbolic;
  trilinos_klu_l_common common;
  kls_separator_analysis separator;
  trilinos_klu_l_numeric *numeric;
  /* staged scale trials: published before the METIS stage so the
     earlier auto-scale promotion can consume them without joining */
  int scale_wanted;
  int scale_unscaled_only;  /* generic compact race: -1 and scale-1 arms */
  int giant_symmetric_scalar_fringe_metis_row;
  int symbolic_only;         /* dense-tail class: the joiner wants the
                                analyze only - skip the scale trials
                                and the serial trial factor (ss1: an
                                unabortable 26.8s klu trial held the
                                take's join) */
  _Atomic int scale_done;
  int scale_count;
  int scale_cand[4];
  trilinos_klu_l_common scale_common[4];
  trilinos_klu_l_numeric *scale_numeric[4];
};
typedef struct kls_metis_race_s kls_metis_race;

typedef struct {
  kls_metis_race *race;
  int slot;
  int candidate;
} kls_scale_trial_arg;

/* Race-side threads must never starve the foreground factor: the race
   overlaps the first factorization by design, and at normal priority
   its analyze + concurrent trial factors turn the serial prestatic
   trial into a scheduling lottery (onetone1: identical 8.4e8-flop
   trials measured 0.25s..18.3s under a live race; 0.12s flat with the
   race disabled). Nice the race threads so the foreground preempts
   them whenever it is runnable. */
static void kls_race_thread_deprioritize(void) {
#ifdef __linux__
  errno = 0;
  (void)setpriority(PRIO_PROCESS, (id_t)syscall(SYS_gettid), 15);
#endif
}

static void *kls_scale_trial_main(void *argp) {
  kls_race_thread_deprioritize();
  kls_scale_trial_arg *a = (kls_scale_trial_arg *)argp;
  kls_metis_race *race = a->race;
  kls_options trial_options = race->options;
  trial_options.scale = a->candidate;
  trilinos_klu_l_common trial_common;
  if (apply_options_to_common(&trial_common, &trial_options) != KLS_OK) {
    return NULL;
  }
  trilinos_klu_l_numeric *trial_numeric =
    trilinos_klu_l_factor(race->col_ptr, race->row_idx, race->values_copy,
                          (trilinos_klu_l_symbolic *)race->scale_symbolic,
                          &trial_common);
  race->scale_cand[a->slot] = a->candidate;
  race->scale_common[a->slot] = trial_common;
  race->scale_numeric[a->slot] = trial_numeric;
  return NULL;
}

static void *kls_metis_race_main(void *arg) {
  /* Deliberately NOT deprioritized: nice-ing this thread stretches the
     bandwidth-heavy analyze across the foreground's whole trial phase
     (nice rations CPU, not bandwidth - mac_econ's column trials read
     9s against their 1.1s steady under a lingering race). The
     scale-trial factor threads it spawns ARE niced, and the stage2
     wait sleeps; those were the init-lottery culprits. */
  kls_metis_race *race = (kls_metis_race *)arg;
  if (race->metis_wanted) {
#ifdef KLS_HAVE_MTMETIS
    const int saved_giant_symmetric_scalar_fringe_metis_ctx =
      kls_giant_symmetric_scalar_fringe_metis_ctx;
    kls_giant_symmetric_scalar_fringe_metis_ctx =
      race->giant_symmetric_scalar_fringe_metis_row;
#endif
#ifdef KLS_HAVE_MTMETIS
    /* the symbolic-only (dense-tail class) joiner is on the one-shot
       critical path waiting for this analyze */
    kls_det_ndp_requested = race->symbolic_only;
#endif
    race->analyze_status = analyze_with_ordering(race->n, race->col_ptr,
                                                 race->row_idx,
                                                 &race->options,
                                                 KLS_ORDERING_METIS,
                                                 &race->symbolic,
                                                 &race->common,
                                                 &race->separator);
#ifdef KLS_HAVE_MTMETIS
    kls_det_ndp_requested = 0;
    kls_giant_symmetric_scalar_fringe_metis_ctx =
      saved_giant_symmetric_scalar_fringe_metis_ctx;
#endif
  }
  atomic_store_explicit(&race->analyze_done, 1, memory_order_release);
  int stage2;
  unsigned stage2_spin = 0;
  while ((stage2 = atomic_load_explicit(&race->stage2,
                                        memory_order_acquire)) == 0) {
    /* this wait spans the whole foreground first factor - a hot spin
       here steals a core from it for seconds */
    if (++stage2_spin > 1024u) {
      struct timespec ts = {0, 200000};
      nanosleep(&ts, NULL);
    } else {
      kls_cpu_relax();
    }
  }
  if (stage2 != 1 || race->symbolic_only) {
    atomic_store_explicit(&race->scale_done, 1, memory_order_release);
    atomic_store_explicit(&race->finished, 1, memory_order_release);
    return NULL;
  }
  if (race->scale_wanted && race->scale_symbolic != NULL) {
    /* the candidates are independent serial factors; run them
       concurrently so the promotion's spin never outlives the
       longest single trial */
    const int candidates[] = {-1, 0, 1, 2};
    kls_scale_trial_arg targs[4];
    pthread_t ttids[4];
    int spawned[4] = {0, 0, 0, 0};
    int nt = 0;
    for (size_t i = 0; i < sizeof(candidates) / sizeof(candidates[0]);
         ++i) {
      if (race->scale_unscaled_only &&
          candidates[i] != -1 && candidates[i] != 1) {
        continue;
      }
      if (candidates[i] == race->options_base_scale) {
        continue;
      }
      targs[nt].race = race;
      targs[nt].slot = nt;
      targs[nt].candidate = candidates[i];
      race->scale_cand[nt] = candidates[i];
      race->scale_numeric[nt] = NULL;
      nt++;
    }
    race->scale_count = nt;
    for (int t = 1; t < nt; ++t) {
      spawned[t] =
        pthread_create(&ttids[t], NULL, kls_scale_trial_main,
                       &targs[t]) == 0;
    }
    if (nt > 0) {
      (void)kls_scale_trial_main(&targs[0]);
    }
    for (int t = 1; t < nt; ++t) {
      if (spawned[t]) {
        pthread_join(ttids[t], NULL);
      } else {
        (void)kls_scale_trial_main(&targs[t]);
      }
    }
  }
  atomic_store_explicit(&race->scale_done, 1, memory_order_release);
  if (!race->metis_wanted || race->analyze_status != KLS_OK ||
      race->symbolic == NULL) {
    atomic_store_explicit(&race->finished, 1, memory_order_release);
    return NULL;
  }
  {
    /* the race trial is on the one-shot critical path for the
       predicted giants; route it through the pipelined kernel under
       the same est-flops floor as foreground first factors */
    const double est = race->symbolic->est_flops;
    const double rn = (double)race->symbolic->n;
    if (race->options.threads >= 2 &&
        getenv("KLS_DISABLE_PIPE_ROUTE") == NULL &&
        (getenv("KLS_KLU_PIPE_FORCE_RACE") != NULL ||
         getenv("KLS_KLU_PIPE_RACE_DEFAULT") != NULL ||
         race->giant_symmetric_scalar_fringe_metis_row ||
         (est > 5.0e9 && rn > 0.0 && est / rn >= 1.0e5))) {
      /* default-on measured WORSE: giants pay 16 x O(n) workspace
         setup on the race thread (circuit5M_dc +1.5s, Freescale1
         +1.7s) and memchip's default-path run produced resid 1.0
         (failed adoption; unreproduced under the force env - do not
         re-default without solving that) */
      kls_klu_pipe_threads =
        race->options.threads > 16 ? 16 : race->options.threads;
    }
  }
  {
    const int saved_pipe_nopanels = kls_klu_pipe_nopanels;
    /* This topology produces millions of two-column panels.  Their compact
       workspace copies execute no batched updates and add seconds over the
       ordinary balanced column pipe.  Keep panel suppression local to the
       structurally admitted race; other heavy factors rely on wider panels. */
    kls_klu_pipe_nopanels =
      race->giant_symmetric_scalar_fringe_metis_row;
    const double kls_rf_t0 = kls_now_seconds();
    race->numeric = trilinos_klu_l_factor(race->col_ptr, race->row_idx,
                                          race->values_copy, race->symbolic,
                                          &race->common);
    if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
      fprintf(stderr, "KLS race trial factor %.3fs (pipe=%d)\n",
              kls_now_seconds() - kls_rf_t0, kls_klu_pipe_threads);
    }
    kls_klu_pipe_nopanels = saved_pipe_nopanels;
  }
  kls_klu_pipe_threads = 0;
  atomic_store_explicit(&race->finished, 1, memory_order_release);
  return NULL;
}

/* Non-blocking: has the race worker produced its trial factorization?
   The factor-time promotion skips (leaving the race alive) instead of
   joining when the worker is still running; the refactor wrapper
   consults again once this reads true, so the first factor never
   blocks on NodeND (ASIC_320ks: 0.69s of a 1.46s init was this join). */
static int kls_metis_race_ready(const kls_solver *solver) {
  return solver != NULL && solver->metis_race != NULL &&
         (!solver->metis_race->active ||
          atomic_load_explicit(&solver->metis_race->finished,
                               memory_order_acquire) != 0);
}

static kls_metis_race *kls_metis_race_take(kls_solver *solver,
                                           double *elapsed) {
  kls_metis_race *race = solver != NULL ? solver->metis_race : NULL;
  if (race == NULL) {
    return NULL;
  }
  solver->metis_race = NULL;
  if (race->active) {
    const double start = kls_now_seconds();
    if (!race->values_signaled) {
      atomic_store_explicit(&race->stage2, 2, memory_order_release);
    }
    pthread_join(race->thread, NULL);
    race->active = 0;
    if (elapsed != NULL) {
      *elapsed += kls_now_seconds() - start;
    }
    if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
      fprintf(stderr, "KLS metis race join wait %.3fs\n",
              kls_now_seconds() - start);
    }
  }
  return race;
}

static void kls_metis_race_free(kls_metis_race *race) {
  if (race == NULL) {
    return;
  }
  for (int i = 0; i < race->scale_count; ++i) {
    if (race->scale_numeric[i] != NULL) {
      trilinos_klu_l_free_numeric(&race->scale_numeric[i],
                                  &race->scale_common[i]);
    }
  }
  if (race->numeric != NULL) {
    trilinos_klu_l_free_numeric(&race->numeric, &race->common);
  }
  if (race->symbolic != NULL) {
    trilinos_klu_l_free_symbolic(&race->symbolic, &race->common);
  }
  kls_separator_analysis_clear(&race->separator);
  free(race->values_copy);
  free(race);
}

/* CKTSO-style static pivoting normally orders the row-matched matrix.  KLS's
   split analyze/factor API cannot start the value-weighted match until factor
   time, while its METIS race has already spent that interval ordering the
   original pattern.  This experiment reuses the raced single-block symbolic
   on the matched coordinates after making its row permutation equal to its
   column order.  The predicted-first builder recomputes the exact LU pattern,
   so the original symbolic estimates remain only conservative routing data.
   Matching and NodeND overlap, and neither a second NodeND, another symbolic
   pass, nor the race's serial numeric trial is required.  The predicted-first
   pivot/residual checks remain the acceptance boundary. */
static int kls_prestatic_reuse_raced_metis_symbolic(
  kls_solver *solver,
  UF_long *col_ptr,
  UF_long *row_idx,
  const kls_options *options,
  trilinos_klu_l_symbolic **symbolic_out,
  trilinos_klu_l_common *common_out,
  kls_separator_analysis *separator_out) {
  if (solver == NULL || col_ptr == NULL || row_idx == NULL ||
      options == NULL || symbolic_out == NULL || common_out == NULL ||
      separator_out == NULL || *symbolic_out != NULL ||
      solver->metis_race == NULL ||
      solver->metis_race->values_signaled) {
    return 0;
  }

  const double begin = kls_now_seconds();
  kls_metis_race *race = solver->metis_race;
  /* Publish symbolic_only before stage2.  The worker's acquire load of
     stage2 then sees the flag and exits immediately after NodeND. */
  race->values_signaled = 1;
  race->symbolic_only = 1;
  atomic_store_explicit(&race->stage2, 1, memory_order_release);
  race = kls_metis_race_take(solver, NULL);
  if (race == NULL || race->analyze_status != KLS_OK ||
      race->symbolic == NULL || race->symbolic->P == NULL ||
      race->symbolic->Q == NULL || race->symbolic->R == NULL ||
      race->symbolic->nblocks != 1u) {
    kls_metis_race_free(race);
    return 0;
  }

  (void)col_ptr;
  (void)row_idx;
  (void)options;
  trilinos_klu_l_symbolic *given = race->symbolic;
  race->symbolic = NULL;
  for (UF_long k = 0; k < solver->n; ++k) {
    given->P[k] = given->Q[k];
  }
  given->do_btf = 0;
  given->nblocks = 1u;
  given->maxblock = solver->n;
  if (given->R != NULL) {
    given->R[0] = 0u;
    given->R[1] = solver->n;
  }
  const int ok = given->P != NULL && given->Q != NULL && given->R != NULL &&
                 race->common.status >= TRILINOS_KLU_OK;
  if (getenv("KLS_TRACE_PRESTATIC") != NULL) {
    fprintf(stderr,
            "KLS pre-static: raced symbolic reuse %.3fms ok=%d "
            "fill=%.3e\n",
            1e3 * (kls_now_seconds() - begin), ok,
            ok ? given->lnz + given->unz : -1.0);
  }
  if (!ok) {
    race->symbolic = given;
    kls_metis_race_free(race);
    return 0;
  }
  *common_out = race->common;
  /* The component labels are expressed in the raced Q-order, which is also
     the coordinate order of the reused symmetric P=Q symbolic.  Row matching
     can introduce dependencies that were not edges in the original graph;
     the row scheduler's dependency-closure and cross-thread-promotion checks
     conservatively move every such target to its pipeline before execution.
     Retaining this tree therefore supplies useful ownership seeds without
     assuming that its original components remain independent. */
  kls_separator_analysis_move(separator_out, &race->separator);
  kls_metis_race_free(race);
  *symbolic_out = given;
  return 1;
}

static void kls_metis_race_abandon(kls_solver *solver) {
  kls_metis_race_free(kls_metis_race_take(solver, NULL));
}

/* The foreground can build a predicted numeric directly from the raced
   symbolic.  For sparse, exceptionally high-work incumbents this is both
   much faster and safer than letting the race enter its serial trial factor
   before the symbolic join (rajat31: 10s predicted METIS build versus a
   roughly 38s race-side KLU trial). */

static int kls_compact_medium_dominant_btf_symbolic_shape(
  const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL || solver->n == 0u ||
      solver->symbolic->nblocks < 8u ||
      solver->symbolic->nblocks > 512u) {
    return 0;
  }
  const double coverage =
    (double)solver->symbolic->maxblock / (double)solver->n;
  return coverage >= 0.95 &&
         solver->symbolic->maxblock >= 30000u &&
         solver->symbolic->maxblock < 60000u &&
         solver->symbolic->est_flops >= 3.0e7 &&
         solver->symbolic->lnz + solver->symbolic->unz >= 1.0e6;
}

/* Dense giant circuits already expose a compact AMD factor, while NodeND and
   its candidate factor are tens of seconds of non-amortizing setup
   (circuit5M: 28.6 s rejected promotion).  Sparse giants stay eligible:
   rajat31 and the Freescale chain materially reduce fill with METIS. */
static int kls_dense_giant_declines_auto_metis(UF_long n,
                                                const UF_long *col_ptr) {
  return col_ptr != NULL && n >= 1000000u && col_ptr[n] / n >= 8u;
}

/* Very large, low-degree circuits can consist mostly of scalar BTF blocks
   around one moderate core.  The generic pre-static deferral leaves this
   class on AMD after its inexpensive matching probe rejects, but a tuned
   NodeND order of the core repays its analysis cost over H100.  Use a tight
   pattern-only gate for the direct ordering decision, then confirm the BTF
   geometry before suppressing later speculative numeric work. */
static int kls_large_low_degree_fragmented_pattern(
  UF_long n,
  const UF_long *col_ptr) {
  if (col_ptr == NULL || n < 600000u || n > 800000u) {
    return 0;
  }
  const double nnz = (double)col_ptr[n];
  return nnz >= 2.0 * (double)n && nnz < 3.0 * (double)n;
}

/* A medium-size dominant BTF block with thousands of small fringe blocks is
   already well served by the incumbent AMD factor.  In this shape, ordering
   the dominant block with nested dissection destroys the useful local
   structure while scale trials retain the same numeric: on the audited
   358K-row case METIS grows L+U from 4.49M to 12.26M entries and the rejected
   METIS/scale consultations add 6.65s to the first refactor.  The symbolic
   work and fill bounds keep this separate from low-work fragmented systems
   and from the single/few-block matrices where METIS is profitable. */

/* A partial structural diagonal is a cheap proposal signal for sparse
   matching-class systems.  It is deliberately insufficient to select the
   policy: the actual AMD/BTF symbolic below must also demonstrate the
   normalized dominant-core work regime. */

/* Sparse AMD analyses with a dominant core and a 10--20% scalar/block
   fringe form a distinct fragmented-BTF class.  Keep this normalized
   structural predicate shared by analyze-time speculation and post-factor
   solve selection; the latter adds measured numeric-work requirements. */
__attribute__((noinline, cold))
static int kls_sparse_fragmented_dominant_btf_symbolic_structure(
  UF_long n,
  UF_long nnz,
  kls_orientation orientation,
  kls_ordering ordering,
  const trilinos_klu_l_symbolic *symbolic) {
  if (symbolic == NULL || orientation != KLS_ORIENTATION_NORMAL ||
      ordering != KLS_ORDERING_AMD || n < 100000u ||
      n > UF_long_max / 5u || nnz < 3u * n || nnz > 5u * n ||
      !symbolic->do_btf || symbolic->nblocks < n / 10u ||
      symbolic->nblocks > n / 5u ||
      (double)symbolic->maxblock < 0.80 * (double)n ||
      (double)symbolic->maxblock > 0.95 * (double)n ||
      (symbolic->structural_rank >= 0 &&
       (UF_long)symbolic->structural_rank != n)) {
    return 0;
  }
  return 1;
}

__attribute__((noinline, cold))
static int kls_sparse_fragmented_dominant_btf_symbolic_work_profile(
  UF_long n,
  UF_long nnz,
  kls_orientation orientation,
  kls_ordering ordering,
  const trilinos_klu_l_symbolic *symbolic) {
  if (!kls_sparse_fragmented_dominant_btf_symbolic_structure(
        n, nnz, orientation, ordering, symbolic)) {
    return 0;
  }
  const double est_fill = symbolic->lnz + symbolic->unz;
  const double est_flops = symbolic->est_flops;
  return est_fill >= 32.0 * (double)n &&
    est_fill <= 64.0 * (double)n &&
    est_flops >= 1000.0 * (double)n &&
    est_flops <= 8192.0 * (double)n;
}

__attribute__((noinline, cold))
static int kls_sparse_fragmented_dominant_btf_symbolic_work(
  const kls_solver *solver) {
  return solver != NULL &&
    kls_sparse_fragmented_dominant_btf_symbolic_work_profile(
      solver->n, solver->nnz, solver->orientation,
      solver->stats.selected_ordering, solver->symbolic);
}

/* Analyze-time start: the METIS analyze stage needs only the pattern,
   so the worker launches when the analysis finishes and overlaps both
   the caller's gap to kls_factor and the whole main-side first factor.
   kls_signal_metis_race_values releases the value-dependent stages. */
static void kls_maybe_start_metis_race(kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL) {
    return;
  }
  const double est_fill = solver->symbolic->lnz + solver->symbolic->unz;
  const double est_flops = solver->symbolic->est_flops;
  const double generic_scale_horizon_work =
    (double)solver->options.expected_refactorizations * (double)solver->nnz;
  /* A compact repeated numeric can afford to compare scaling even though its
     absolute symbolic work is below the historical large-race floor.  Width,
     estimated work, and caller lifecycle are the complete admission: the
     worker still factors the real values and the ordinary numeric-quality
     gate decides whether anything is published. */
  const int generic_compact_scale_race =
    solver->options.scale == KLS_SCALE_AUTO &&
    kls_repeated_update_workload(&solver->options) &&
    solver->options.threads > 1 &&
    solver->orientation == KLS_ORIENTATION_NORMAL &&
    solver->row_perm == NULL && solver->n <= UINT32_C(0x1000) &&
    solver->nnz <= (UF_long)UINT16_MAX &&
    solver->symbolic->lnz > 0.0 &&
    solver->symbolic->lnz <= (double)UINT16_MAX &&
    solver->symbolic->unz > 0.0 &&
    solver->symbolic->unz <= (double)UINT16_MAX &&
    est_flops > 0.0 && isfinite(generic_scale_horizon_work) &&
    generic_scale_horizon_work >= 0.5 * est_flops + est_fill;
  if (solver != NULL && solver->symbolic != NULL &&
      solver->symbolic->structural_rank != KLS_KLU_EMPTY &&
      solver->symbolic->structural_rank < solver->n) {
    /* Ordering cannot repair a known structural-rank deficit.  Without an
       explicit singular-completion request, a background NodeND/scale race
       has no admissible recurring factor to promote and can make analyze +
       destroy wait on an otherwise immediate rejection. */
    return;
  }
  {
    /* Generic AUTO has already compared the ND symbolic once, after choosing
       orientation.  Do not build two complete numerics merely to repeat that
       ordering decision at factor time.  The compact scale race below stays
       on the selected symbolic and covers only the unscaled and row-max arms;
       mode 0 duplicates the unscaled numeric and mode 2 is the usual AUTO
       incumbent. */
    if (!generic_compact_scale_race &&
        getenv("KLS_ENABLE_GENERIC_SCALE_RACE") == NULL) {
      return;
    }
  }
  if (solver->metis_race != NULL || solver->auto_metis_checked ||
      solver->options.ordering != KLS_ORDERING_AUTO ||
      solver->user_col_perm != NULL || solver->symbolic == NULL ||
      (solver->n < 20000 && !generic_compact_scale_race) ||
      getenv("KLS_DISABLE_METIS_RACE") != NULL ||
      kls_dense_giant_declines_auto_metis(solver->n, solver->col_ptr) ||
      (kls_sparse_fragmented_dominant_btf_symbolic_work(solver) &&
       getenv("KLS_DISABLE_SPARSE_FRAGMENTED_RACE_SUPPRESSION") == NULL) ||
      solver->medium_partial_static_metis_path ||
      solver->medium_spike_minfill_path || solver->large_sparse_amf3_path ||
      solver->large_bounded_no_btf_amf_path) {
    return;
  }
  if (!generic_compact_scale_race &&
      getenv("KLS_ENABLE_COMPACT_MEDIUM_NUMERIC_TRIALS") == NULL &&
      kls_compact_medium_dominant_btf_symbolic_shape(solver)) {
    /* This modest-work, nearly single-block BTF shape already exposes its
       useful parallelism to the exact EGraph.  A background NodeND plus
       scale-factor race instead contends with the AMD first factor and then
       delays first-refactor preparation at the join (ckt11752_dc_1: 335ms
       versus 406ms over the 100-cycle horizon with the race retained). */
    return;
  }
  if (!generic_compact_scale_race && solver->n <= 200000 &&
      solver->col_ptr != NULL &&
      solver->col_ptr[solver->n] <= 4000000 &&
      solver->options.static_pivoting) {
    /* the block-structured class (TSOPF/case9 family) adopts its own
       ordering at factor time and always discards the race - and on
       its fused-dense-row patterns the race's NodeND grinds for tens
       of seconds, which the adoption's pattern-lifetime join then
       waits out on the one-shot path (TSOPF_FS_b39_c19: 54.5s of a
       54.8s init).  Same pattern-only detector as the choose fast
       path (~10ms at this size gate). */
    UF_long *bs_perm = NULL;
    UF_long *bs_comp = NULL;
    const int bs_hit =
      kls_build_block_structured_order(solver->n, solver->col_ptr,
                                       solver->row_idx, &bs_perm,
                                       &bs_comp);
    if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
      fprintf(stderr, "KLS race gate: block detector=%d n=%ld\n",
              bs_hit, (long)solver->n);
    }
    if (bs_hit) {
      free(bs_perm);
      free(bs_comp);
      return;
    }
  }
  /* Two race clients: the METIS promotion (giant class: the serial
     trial costs seconds) and the scale trials (mid class upward: three
     serial factors, e.g. ASIC_320ks 0.55s, usually rejected).  The race
     copies the values and burns a core, so small rows keep the serial
     promotion path. */
  /* mid class only: each raced candidate is a full serial factor at
     the CURRENT ordering, which on the giant class is the multi-minute
     path the predicted machinery exists to avoid */
  int scale_wanted = generic_compact_scale_race ||
    (solver->options.scale == KLS_SCALE_AUTO &&
    !solver->auto_scale_checked &&
    est_fill >= 1.0e6 && est_fill <= 2.0e7 &&
    (est_flops <= 5.0e9 || est_flops <= 0.0));
  if (scale_wanted &&
      symbolic_is_fragmented_many_block_unscaled_candidate(
        solver->n, solver->symbolic)) {
    /* The two ASIC_680k operating points have already selected the
       residual-valid unscaled numeric.  Their raced scale factors are never
       retained, but the foreground still waits for them on the first changed
       input (about 0.28s on ASIC_680ks). */
    scale_wanted = 0;
  }

  if (scale_wanted &&
      is_large_very_low_degree_full_diagonal_pattern(
        solver->n, solver->col_ptr, solver->row_idx)) {
    /* The pattern policy has already selected an unscaled METIS numeric for
       this full-diagonal, degree-at-most-eight class.  The later selector also
       rejects rescaling once that predicted numeric exists, so racing three
       scale factors can only contend with the accepted first factor and leave
       a join for the first public refactor (G2_circuit: about 0.9 s total). */
    scale_wanted = 0;
  }
  if (!scale_wanted) {
    return;
  }
  kls_metis_race *race = (kls_metis_race *)calloc(1, sizeof(*race));
  if (race == NULL) {
    return;
  }
  race->n = solver->n;
  race->col_ptr = solver->col_ptr;
  race->row_idx = solver->row_idx;
  race->options = solver->options;
  race->options.ordering = KLS_ORDERING_METIS;
  race->options.use_btf = solver->symbolic->do_btf ? 1 : 0;
  race->giant_symmetric_scalar_fringe_metis_row = 0;
  race->scale_symbolic = solver->symbolic;
  race->metis_wanted = 0;
  race->scale_wanted = scale_wanted;
  race->scale_unscaled_only = 1;
  race->symbolic_only = 0;
  atomic_init(&race->scale_done, 0);
  atomic_init(&race->stage2, 0);
  race->analyze_status =
    KLS_OK;
  if (pthread_create(&race->thread, NULL, kls_metis_race_main, race) != 0) {
    free(race);
    return;
  }
  race->active = 1;
  solver->metis_race = race;
  if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
    fprintf(stderr, "KLS race: launched (maybe_start) n=%ld\n",
            (long)solver->n);
  }
}

/* Giant-class early start: the pattern is fixed as soon as the
   orientation heuristic skips the transpose candidate, so the METIS
   analyze can overlap the main AMD/BTF analysis itself.  Symbolic
   estimates do not exist yet; the n floor selects the class where the
   race is always worth a core, and the promotion gate discards the
   result if the analysis ends up choosing METIS on its own. */

/* Values are known: copy them, capture the settled scale, and release
   the worker's value-dependent stages. */
static void kls_signal_metis_race_values(kls_solver *solver,
                                         const double *numeric_values) {
  kls_metis_race *race = solver->metis_race;
  if (race == NULL || race->values_signaled) {
    return;
  }
  if (race->scale_unscaled_only && solver->common.scale <= 0) {
    /* The foreground value census already selected the sole candidate this
       generic scale-only race could build.  Abort and join while the worker
       is still in its stage-2 wait; carrying an empty race into the first
       refactor otherwise perturbs deferred preparation by a full consult. */
    race->values_signaled = 1;
    atomic_store_explicit(&race->stage2, 2, memory_order_release);
    kls_metis_race_abandon(solver);
    return;
  }
  race->values_signaled = 1;
  if (race->symbolic_only) {
    atomic_store_explicit(&race->stage2, 1, memory_order_release);
    return;
  }

  race->values_copy =
    (double *)malloc((size_t)solver->nnz * sizeof(*race->values_copy));
  if (race->values_copy == NULL) {
    atomic_store_explicit(&race->stage2, 2, memory_order_release);
    return;
  }
  memcpy(race->values_copy, numeric_values,
         (size_t)solver->nnz * sizeof(*race->values_copy));
  race->options.scale = (int)solver->common.scale;
  race->options_base_scale = (int)solver->common.scale;
  /* Value-pressure selection runs after analyze.  Keep scale candidates on
     the live pivot contract so the race changes one axis only. */
  race->options.pivot_tolerance = solver->common.tol;
  atomic_store_explicit(&race->stage2, 1, memory_order_release);
}

/* A compact plain-frame factor may trade KLU's recurring scale pass for an
   unscaled retained-row numeric even when the scaled incumbent is condition
   sensitive.  Admit that comparison only when the complete residual-policed
   solve is representable: the public/internal value frames coincide, every
   persistent row/factor descriptor fits its kernel width, and AUTO selected
   the tight pivot line whose current-value snapshot is retained.  No block
   count, density window, ordering, matrix dimension identity, or family
   property participates in the decision. */
static int kls_certified_unscaled_lifecycle_trial_capable(
  const kls_solver *solver) {
  return  solver != NULL &&
    solver->symbolic != NULL && solver->numeric != NULL &&
    kls_repeated_update_workload(&solver->options) &&
    solver->options.scale == KLS_SCALE_AUTO &&
    solver->options.threads > 1 &&
    solver->orientation == KLS_ORIENTATION_NORMAL &&
    solver->input_format == KLS_INPUT_CSC &&
    solver->input_to_csc == NULL && solver->row_perm == NULL &&
    solver->user_col_perm == NULL && solver->row_scale == NULL &&
    solver->col_scale == NULL && !solver->numeric_is_predicted &&
    solver->pivot_nudge_count == 0u &&
    solver->common.kls_perturb_count == 0u &&
    solver->common.scale > 0 && solver->numeric->Rs != NULL &&
    solver->common.tol >= 1.0e-8 && solver->common.tol <= 1.0e-7 &&
    solver->common.tol < solver->options.pivot_tolerance &&
    solver->col_ptr != NULL && solver->row_idx != NULL &&
    solver->n <= UINT32_C(0x1000) &&
    solver->nnz <= (UF_long)UINT16_MAX &&
    solver->numeric->lnz <= (UF_long)UINT16_MAX &&
    solver->numeric->unz <= (UF_long)UINT16_MAX;
}

/* Before the first numeric exists, the same packed lifecycle can avoid a
   throwaway scaled incumbent when its symbolic streams already fit every
   narrow descriptor and the declared scale-pass horizon repays the factor.
   This is only a provisional representation/economics certificate.  The
   realized numeric is checked after factorization, and every accepted solve
   retains the residual contract plus recovery to robust_scale. */
static int kls_initial_certified_unscaled_lifecycle_candidate(
  const kls_solver *solver, int robust_scale) {
  if (solver == NULL ||
      solver->symbolic == NULL || solver->numeric != NULL ||
      robust_scale <= 0 || solver->options.scale != KLS_SCALE_AUTO ||
      solver->options.ordering != KLS_ORDERING_AUTO ||
      solver->options.threads <= 1 ||
      !kls_repeated_update_workload(&solver->options) ||
      solver->options.expected_refactorizations < 64 ||
      solver->orientation != KLS_ORIENTATION_NORMAL ||
      solver->input_format != KLS_INPUT_CSC ||
      solver->input_to_csc != NULL || solver->row_perm != NULL ||
      solver->user_col_perm != NULL || solver->row_scale != NULL ||
      solver->col_scale != NULL ||
      (solver->stats.selected_ordering != KLS_ORDERING_AMF &&
       solver->stats.selected_ordering != KLS_ORDERING_AMMF) ||
      solver->common.tol < 1.0e-8 || solver->common.tol > 1.0e-7 ||
      solver->common.tol >= solver->options.pivot_tolerance ||
      solver->n > UINT32_C(0x1000) ||
      solver->nnz > (UF_long)UINT16_MAX ||
      !(solver->symbolic->lnz > 0.0) ||
      !(solver->symbolic->unz > 0.0) ||
      solver->symbolic->lnz + solver->symbolic->unz >
        2.0 * (double)UINT16_MAX ||
      !(solver->symbolic->est_flops > 0.0)) {
    return 0;
  }
  const double symbolic_fill =
    solver->symbolic->lnz + solver->symbolic->unz;
  const double avoided_scale_work =
    (double)solver->options.expected_refactorizations *
    (double)solver->nnz;
  return isfinite(symbolic_fill) && symbolic_fill > 0.0 &&
    /* Symbolic flop estimates deliberately upper-bound threshold-pivoted
       work and can be several times the realized factor.  Starting with the
       candidate is not an extra factor unless its realized checks reject it,
       so require the declared avoided passes to cover half that conservative
       bound; the exact numeric limits below remain authoritative. */
    avoided_scale_work >=
      0.5 * (solver->symbolic->est_flops + symbolic_fill) &&
    getenv("KLS_DISABLE_INITIAL_CERTIFIED_UNSCALED_LIFECYCLE") == NULL;
}

/* For a cache-resident BTF in the caller's compressed frame, KLU row scaling
   adds a complete value pass to every numeric and a scale pass to every
   solve.  Start AUTO with the no-scale/no-check KLU mode when the symbolic
   work and retained descriptors are bounded.  The first numeric still has
   the ordinary pivot/fill gate below and retains the value-selected scale as
   a failure recovery; solve-contract classification remains unchanged. */
static int kls_initial_low_work_btf_public_unscaled_candidate(
  const kls_solver *solver, int robust_scale) {
  return solver != NULL && solver->symbolic != NULL &&
    solver->numeric == NULL && robust_scale > 0 &&
    solver->options.scale == KLS_SCALE_AUTO &&
    solver->input_to_csc == NULL &&
    ((solver->input_format == KLS_INPUT_CSC &&
      solver->orientation == KLS_ORIENTATION_NORMAL) ||
     (solver->input_format == KLS_INPUT_CSR &&
      solver->orientation == KLS_ORIENTATION_TRANSPOSE)) &&
    solver->row_perm == NULL && solver->user_col_perm == NULL &&
    solver->row_scale == NULL && solver->col_scale == NULL &&
    solver->symbolic->nblocks > 1u && solver->n > 0u &&
    solver->n <= UINT32_C(0x1000) &&
    solver->nnz > 0u && solver->nnz <= (UF_long)UINT16_MAX &&
    solver->symbolic->structural_rank == solver->n &&
    solver->symbolic->est_flops > 0.0 &&
    solver->symbolic->est_flops < 100000.0 &&
    getenv("KLS_DISABLE_INITIAL_LOW_WORK_BTF_PUBLIC_UNSCALED") == NULL;
}

/* A scaled factor can be numerically condition-sensitive yet structurally
   stable, while scaling also excludes the branch-reduced unscaled EGraph
   kernels on every later update.  For a long caller-declared lifecycle,
   compare one unscaled numeric when the trial itself has bounded factor and
   storage work.  This only authorizes measurement: adoption still requires
   close realized fill/work/pivoting plus the complete residual probe below,
   and a low-rcond adoption retains the scaled recovery contract. */
static int kls_repeated_scaled_egraph_unscaled_trial_enabled(
  const kls_solver *solver) {
  if (solver == NULL ||
      solver->symbolic == NULL || solver->numeric == NULL ||
      solver->options.scale != KLS_SCALE_AUTO ||
      !kls_repeated_update_workload(&solver->options) ||
      solver->options.expected_refactorizations < 64 ||
      solver->options.threads <= 1 ||
      !solver->auto_scale_unscaled_trial_certified ||
      solver->orientation != KLS_ORIENTATION_NORMAL ||
      solver->input_format != KLS_INPUT_CSC ||
      solver->input_to_csc != NULL || solver->row_perm != NULL ||
      solver->user_col_perm != NULL || solver->row_scale != NULL ||
      solver->col_scale != NULL || solver->numeric_is_predicted ||
      solver->pivot_nudge_count != 0u ||
      solver->common.kls_perturb_count != 0u ||
      solver->common.scale <= 0 || solver->numeric->Rs == NULL ||
      solver->n < 5000u || solver->nnz == 0u ||
      !(solver->common.rcond > 0.0) ||
      solver->common.rcond >= sqrt(DBL_EPSILON) ||
      !(solver->common.flops > 0.0) ||
      solver->common.flops > 5.0e9 ||
      solver->common.flops <
        KLS_EGRAPH_REFACTOR_MIN_FLOPS_PER_THREAD *
          (double)solver->options.threads) {
    return 0;
  }
  const double fill =
    (double)(solver->numeric->lnz + solver->numeric->unz);
  return fill > 0.0 && fill <= 2.0e7 &&
    getenv("KLS_DISABLE_SCALED_EGRAPH_UNSCALED_TRIAL") == NULL;
}

static int kls_repeated_scaled_overhead_trial_enabled(
  const kls_solver *solver) {
  if (solver == NULL ||
      solver->numeric == NULL || !kls_repeated_update_workload(&solver->options) ||
      solver->common.scale <= 0 || solver->numeric->Rs == NULL ||
      solver->nnz == 0u || !(solver->common.flops > 0.0) ||
      solver->pivot_nudge_count != 0u ||
      solver->common.kls_perturb_count != 0u) {
    return 0;
  }
  const double horizon =
    (double)solver->options.expected_refactorizations;
  const double fill =
    (double)(solver->numeric->lnz + solver->numeric->unz);
  const double scale_pass_work = horizon * (double)solver->nnz;
  const int amortized_scale_overhead =
    getenv("KLS_DISABLE_AMORTIZED_SCALE_OVERHEAD_TRIAL") == NULL &&
    solver->options.expected_refactorizations >= 64 &&
    solver->options.threads > 1 &&
    solver->common.flops >=
      KLS_EGRAPH_REFACTOR_MIN_FLOPS_PER_THREAD *
        (double)solver->options.threads &&
    isfinite(fill) && fill > 0.0 &&
    scale_pass_work >= 8.0 * (solver->common.flops + fill);
  const int conventional_scale_overhead =
    scale_pass_work >= 5.0e7 &&
    solver->common.flops <= 8.0 * (double)solver->nnz &&
    fill <= 4.0 * (double)solver->nnz;
  /* An unscaled numeric can expose the packed retained-row kernel.  Permit
     that binary scale comparison when either the value census or the full
     residual-policed lifecycle capability certifies the experiment, all
     narrow representation limits already hold, and the caller-declared
     avoided scale passes cover two full realized factor-work budgets.  This
     only authorizes measurement: low-condition candidates must pass the
     true-current-matrix residual probe and retain scaled recovery below. */
  const int packed_row_opportunity =
    solver->options.threads > 1 && solver->n <= UINT32_C(0x1000) &&
    solver->nnz <= (UF_long)UINT16_MAX &&
    solver->numeric->lnz <= (UF_long)UINT16_MAX &&
    solver->numeric->unz <= (UF_long)UINT16_MAX &&
    (solver->auto_scale_unscaled_trial_certified ||
     kls_certified_unscaled_lifecycle_trial_capable(solver) ||
     getenv("KLS_ENABLE_UNCERTIFIED_PACKED_ROW_SCALE_TRIAL") != NULL) &&
    scale_pass_work >= 2.0 * (solver->common.flops + fill);
  return horizon > 1.0 && isfinite(horizon) &&
    (amortized_scale_overhead || conventional_scale_overhead ||
     packed_row_opportunity ||
     kls_repeated_scaled_egraph_unscaled_trial_enabled(solver));
}

static int should_try_auto_scale(const kls_solver *solver) {
  if (solver == NULL || solver->auto_scale_checked ||
      solver->options.scale != KLS_SCALE_AUTO ||
      solver->numeric == NULL) {
    return 0;
  }
  const int repeated_generic_scale_trial =
    kls_repeated_update_workload(&solver->options);
  const int repeated_scaled_overhead_trial =
    kls_repeated_scaled_overhead_trial_enabled(solver);
  if (solver->auto_scale_value_certified &&
      solver->common.scale <= 0) {
    /* The current values already supplied a usable near-complete diagonal
       across a wide central magnitude range, and the unscaled factor
       succeeded.  Retrying every scaled mode would spend up to three more
       complete numerics to rediscover that value decision.  A failed or
       singular unscaled factor is retried robustly before this point and no
       longer has an unscaled Common, so it cannot be suppressed here. */
    return 0;
  }
  if (solver->n < 20000u && !repeated_generic_scale_trial) {
    return 0;
  }

  if (solver->metis_promotion_validated) {
    /* A timed promotion just measured this numeric's refactor at the
       current scale; unmeasured scale trial factors cannot be trusted
       against a measured config (ASIC_320k regressed -1 -> 0 that way). */
    return 0;
  }
  if (solver->row_perm != NULL) {
    /* Statically matched numerics chose their scale mode with the match;
       re-scaling flattens the column magnitude contrast that keeps the
       matched diagonal dominant, and the scaled trial factorization can
       explode structurally on mostly-missing-diagonal matrices. */
    return 0;
  }
  if (solver->numeric_is_predicted) {
    /* A predicted-pattern numeric already proved a machine-precision
       residual without scaling, and the serial scaled trial factorization
       costs minutes at this size. */
    return 0;
  }
  if (
      getenv("KLS_DISABLE_STABLE_UNSCALED_SCALE_GUARD") == NULL &&
      solver->common.scale <= 0 &&
      solver->common.noffdiag == 0 &&
      solver->pivot_nudge_count == 0u &&
      solver->common.kls_perturb_count == 0u &&
      isfinite(solver->common.rcond) && solver->common.rcond >= 1.0e-6) {
    /* The requested diagonal pivot sequence was accepted exactly and exposes
       no structural defect for row scaling to repair.  Any fill change would
       require scaling to replace that stable pivot sequence, while the
       retained numeric is already comfortably above the accuracy floor
       without a nudge or perturbation.  Do not build up to three complete
       factors merely because a pessimistic symbolic flop estimate makes the
       realized diagonal-pivot numeric look inflated.  Failed/singular,
       pivoted, nudged, perturbed, and genuinely ill-conditioned numerics
       remain eligible for the ordinary trials. */
    return 0;
  }

  const double flops = solver->common.flops;
  const UF_long fill = solver->numeric->lnz + solver->numeric->unz;
  const double symbolic_fill = solver->symbolic != NULL
    ? solver->symbolic->lnz + solver->symbolic->unz : 0.0;
  const double symbolic_flops = solver->symbolic != NULL
    ? solver->symbolic->est_flops : 0.0;
  if (solver->common.scale > 0 &&
      solver->symbolic != NULL && solver->symbolic->n > 0u) {
    const UF_long stable_pivot_cap = solver->n / 10000u + 16u;
    const int pivot_stable = solver->common.noffdiag <= stable_pivot_cap;
    const int realized_economics_stable =
      symbolic_fill > 0.0 && (double)fill <= 1.10 * symbolic_fill &&
      (!(symbolic_flops > 0.0) || flops <= 1.10 * symbolic_flops);
    if ((pivot_stable || realized_economics_stable) &&
        !repeated_scaled_overhead_trial) {
      /* Scaling has already produced a numerically stable factor without
         pivot detours, or a realized factor whose storage and arithmetic
         remain inside the symbolic uncertainty band.  Alternative scale
         factors are full numeric builds and can pivot catastrophically
         after a small value update; there is no measured runtime defect
         here for them to repair. */
      return 0;
    }
  }
  if (repeated_generic_scale_trial && solver->symbolic != NULL &&
      ((symbolic_flops > 0.0 && flops >= 2.0 * symbolic_flops) ||
       (symbolic_fill > 0.0 && (double)fill >= 1.50 * symbolic_fill) ||
       (symbolic_flops > 0.0 && symbolic_fill > 0.0 &&
        flops >= 1.10 * symbolic_flops &&
        (double)fill >= 1.10 * symbolic_fill))) {
    /* A factor that materially inflates beyond its symbolic storage/work is
       concrete evidence that the current numeric scaling did not preserve
       the analyzed pivot economy.  Over a caller-declared repeated horizon,
       compare the other KLU scaling modes and retain one only through the
       ordinary realized-numeric gate.  The trigger is an outcome ratio, not
       a dimension, density, or matrix-family classifier. */
    return 1;
  }
  if (repeated_scaled_overhead_trial) {
    /* The numeric arithmetic is so small that recomputing and applying KLU's
       row scales is itself a material recurring kernel.  Over the declared
       lifecycle, one unscaled full factor is bounded by the avoided O(nnz)
       passes.  Numeric fill/work, pivot count, and rcond still arbitrate the
       candidate below; this predicate only authorizes the measurement. */
    return 1;
  }
  if (repeated_generic_scale_trial) {
    /* Generic repeated AUTO reaches this point only when neither realized
       pivot inflation nor measured scale-pass economics proposed a trial.
       Matrix size alone is not evidence that three more full numerics will
       improve the retained lifecycle. */
    return 0;
  }

  return flops >= 1.0e8 && fill >= 1500000;
}

static int maybe_select_auto_scale(kls_solver *solver,
                                   double *elapsed,
                                   const double *numeric_values,
                                   int race_invalid) {
  if (!should_try_auto_scale(solver)) {
    return 0;
  }
  solver->auto_scale_checked = 1;

  /* the promotion race factors the scale candidates on a worker thread
     and publishes them ahead of its METIS stage; consume those instead
     of refactoring serially unless an earlier promotion invalidated
     the raced inputs */
  kls_metis_race *race = solver->metis_race;
  int raced = 0;
  if (race != NULL && race->scale_wanted && !race_invalid) {
    const double spin_start = kls_now_seconds();
    while (!atomic_load_explicit(&race->scale_done,
                                 memory_order_acquire)) {
      kls_cpu_relax();
    }
    *elapsed += kls_now_seconds() - spin_start;
    raced = 1;
  }

  int accepted = 0;
  const int repeated_generic_scale_trial =
    kls_repeated_update_workload(&solver->options);
  const int repeated_scaled_overhead_trial =
    kls_repeated_scaled_overhead_trial_enabled(solver);
  const int candidates[] = {-1, 0, 1, 2};
  for (size_t i = 0; i < sizeof(candidates) / sizeof(candidates[0]); ++i) {
    if (candidates[i] == (int)solver->common.scale) {
      continue;
    }
    if (repeated_scaled_overhead_trial && candidates[i] != -1) {
      /* Scale -1 is the lowest-overhead representative of the two unscaled
         KLU modes.  Do not pay three equivalent full factors to answer a
         binary recurring-cost question. */
      continue;
    }

    trilinos_klu_l_common trial_common;
    trilinos_klu_l_numeric *trial_numeric = NULL;
    int refactor_stale_raced_candidate = 0;
    if (raced) {
      int have = 0;
      for (int s = 0; s < race->scale_count; ++s) {
        if (race->scale_cand[s] == candidates[i]) {
          trial_common = race->scale_common[s];
          trial_numeric = race->scale_numeric[s];
          race->scale_numeric[s] = NULL;
          have = 1;
          break;
        }
      }
      if (!have) {
        continue;
      }
      const int numeric_controls_match =
        trial_common.tol == solver->common.tol &&
        trial_common.halt_if_singular == solver->common.halt_if_singular &&
        trial_common.kls_static_perturb ==
          solver->common.kls_static_perturb &&
        trial_common.kls_zero_pivot_replacement ==
          solver->common.kls_zero_pivot_replacement;
      if (trial_numeric != NULL && !numeric_controls_match) {
        /* The scale race starts before the foreground factor, so analysis
           may subsequently select a different pivot tolerance or numeric
           safety contract.  Such a raced factor changed more than the scale
           axis and cannot be compared or published as a scale candidate.
           Recompute this one arm from the live Common; missing raced arms
           remain skipped, preserving the bounded generic race portfolio. */
        if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
          fprintf(stderr,
                  "KLS scale trial mode=%d stale raced controls"
                  " tol=%.3e/%.3e static=%.3e/%.3e; refactor live\n",
                  candidates[i], trial_common.tol, solver->common.tol,
                  trial_common.kls_static_perturb,
                  solver->common.kls_static_perturb);
        }
        trilinos_klu_l_free_numeric(&trial_numeric, &trial_common);
        refactor_stale_raced_candidate = 1;
      }
    }
    if (!raced || refactor_stale_raced_candidate) {
      /* Change one policy axis.  Reinitializing Common from public options
         discards analyze-selected factor controls, so the resulting
         "scale" candidate can have a different pivot path even at the same
         ordering/tolerance.  Clone the live configuration and replace only
         KLU's scale mode; factor() refreshes its status/statistic fields. */
      trial_common = solver->common;
      trial_common.scale = candidates[i];

      const double start = kls_now_seconds();
      trial_numeric =
        trilinos_klu_l_factor(solver->col_ptr, solver->row_idx,
                              (double *)numeric_values, solver->symbolic,
                              &trial_common);
      *elapsed += kls_now_seconds() - start;
      if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
        fprintf(stderr, "KLS scale trial mode=%d factor=%.3fs status=%ld\n",
                candidates[i], kls_now_seconds() - start,
                (long)trial_common.status);
      }
    }
    if (trial_numeric == NULL || trial_common.status < 0 ||
        trial_common.status == TRILINOS_KLU_SINGULAR) {
      if (trial_numeric != NULL) {
        trilinos_klu_l_free_numeric(&trial_numeric, &trial_common);
      }
      continue;
    }

    (void)trilinos_klu_l_flops(solver->symbolic, trial_numeric, &trial_common);
    if (repeated_scaled_overhead_trial || repeated_generic_scale_trial) {
      (void)trilinos_klu_l_rcond(solver->symbolic, trial_numeric,
                                 &trial_common);
    }
    const double incumbent_fill =
      (double)(solver->numeric->lnz + solver->numeric->unz);
    const double trial_fill =
      (double)(trial_numeric->lnz + trial_numeric->unz);
    const double unscaled_extra_factor_work =
      trial_common.flops > solver->common.flops
        ? trial_common.flops - solver->common.flops : 0.0;
    /* KLU scaling is not one sparse-value operation: it first reduces each
       row (absolute value plus sum/max) and the numeric walk then consumes
       the resulting row factors while gathering the matrix.  Charge the
       unavoidable two sparse passes plus the row-vector pass when comparing
       them with an unscaled factor's extra arithmetic. */
    const double scaled_recurring_overhead_work =
      2.0 * (double)solver->nnz + (double)solver->n;
    const double unscaled_recurring_work_saved =
      scaled_recurring_overhead_work - unscaled_extra_factor_work;
    const double unscaled_lifecycle_work_saved =
      (double)solver->options.expected_refactorizations *
      unscaled_recurring_work_saved;
    const int amortized_scale_overhead_trial =
      getenv("KLS_DISABLE_AMORTIZED_SCALE_OVERHEAD_TRIAL") == NULL &&
      solver->options.expected_refactorizations >= 64 &&
      solver->options.threads > 1 &&
      solver->common.flops >=
        KLS_EGRAPH_REFACTOR_MIN_FLOPS_PER_THREAD *
          (double)solver->options.threads &&
      isfinite(incumbent_fill) && incumbent_fill > 0.0 &&
      (double)solver->options.expected_refactorizations *
          (double)solver->nnz >=
        8.0 * (solver->common.flops + incumbent_fill);
    const int relaxed_unscaled_lifecycle_economics =
      repeated_scaled_overhead_trial && trial_common.scale <= 0 &&
      trial_numeric->Rs == NULL && trial_fill <= 1.05 * incumbent_fill &&
      trial_common.flops > 0.0 && solver->common.flops > 0.0 &&
      unscaled_recurring_work_saved > 0.0 &&
      isfinite(unscaled_lifecycle_work_saved) &&
      unscaled_lifecycle_work_saved >= 3.0 * trial_common.flops &&
      trial_common.noffdiag <=
        solver->common.noffdiag + solver->common.noffdiag / 20u + 16u &&
      (amortized_scale_overhead_trial ||
       !(solver->common.rcond > 0.0) || !(trial_common.rcond > 0.0) ||
       trial_common.rcond >= 0.01 * solver->common.rcond);
    const int egraph_unscaled_specialization_economics =
      kls_repeated_scaled_egraph_unscaled_trial_enabled(solver) &&
      repeated_scaled_overhead_trial && trial_common.scale <= 0 &&
      trial_numeric->Rs == NULL && trial_fill <= 1.02 * incumbent_fill &&
      trial_common.flops > 0.0 && solver->common.flops > 0.0 &&
      trial_common.flops <= 1.05 * solver->common.flops &&
      trial_common.noffdiag <= solver->common.noffdiag +
        solver->common.noffdiag / 10u + 16u;
    const int unscaled_overhead_economics =
      (repeated_scaled_overhead_trial && trial_common.scale <= 0 &&
       trial_numeric->Rs == NULL && trial_fill <= 1.02 * incumbent_fill &&
       trial_common.flops > 0.0 && solver->common.flops > 0.0 &&
       trial_common.flops <= 1.05 * solver->common.flops &&
       trial_common.noffdiag <= solver->common.noffdiag &&
       (!(solver->common.rcond > 0.0) || !(trial_common.rcond > 0.0) ||
        trial_common.rcond >= 0.01 * solver->common.rcond)) ||
      relaxed_unscaled_lifecycle_economics ||
      egraph_unscaled_specialization_economics;
    const int certified_unscaled_lifecycle_trial =
      kls_certified_unscaled_lifecycle_trial_capable(solver);
    const int condition_sensitive_unscaled_trial =
      trial_common.scale <= 0 && trial_numeric->Rs == NULL &&
      ((solver->common.rcond > 0.0 &&
        solver->common.rcond < sqrt(DBL_EPSILON)) ||
       (trial_common.rcond > 0.0 &&
        trial_common.rcond < 32.0 * DBL_EPSILON));
    const int residual_guarded_unscaled_trial =
      condition_sensitive_unscaled_trial ||
      relaxed_unscaled_lifecycle_economics ||
      egraph_unscaled_specialization_economics ||
      (certified_unscaled_lifecycle_trial &&
       condition_sensitive_unscaled_trial);
    int unscaled_overhead_residual_ok = 0;
    if (unscaled_overhead_economics && residual_guarded_unscaled_trial &&
        getenv("KLS_DISABLE_GUARDED_UNSCALED_SCALE_TRIAL") == NULL) {
      /* Removing a recurring scale pass can repay under the broader
         lifecycle economics, including when the scaled incumbent is
         condition-sensitive, but fill/work similarity alone is not an
         accuracy certificate.  Validate the complete candidate against the
         current matrix before publication.  Low-rcond adoptions below also
         arm a floor that can restore the scaled incumbent after a later
         value update leaves the residual-certified regime. */
      trilinos_klu_l_numeric *saved_numeric = solver->numeric;
      const trilinos_klu_l_common saved_common = solver->common;
      solver->numeric = trial_numeric;
      solver->common = trial_common;
      unscaled_overhead_residual_ok =
        kls_direct_klu_numeric_residual_probe(solver, numeric_values);
      trial_common = solver->common;
      solver->numeric = saved_numeric;
      solver->common = saved_common;
    }
    const int unscaled_overhead_candidate_is_better =
      unscaled_overhead_economics &&
      (!residual_guarded_unscaled_trial ||
       unscaled_overhead_residual_ok);
    const int certified_unscaled_overhead_candidate_is_better =
      certified_unscaled_lifecycle_trial &&
      repeated_scaled_overhead_trial && trial_common.scale <= 0 &&
      trial_numeric->Rs == NULL && trial_fill <= 1.02 * incumbent_fill &&
      trial_common.flops > 0.0 && solver->common.flops > 0.0 &&
      trial_common.flops <= 1.05 * solver->common.flops &&
      trial_common.noffdiag <= solver->common.noffdiag &&
      trial_common.rcond > 0.0 && isfinite(trial_common.rcond);
    const int unscaled_realized_candidate_is_better =
      repeated_generic_scale_trial && solver->common.scale > 0 &&
      trial_common.scale <= 0 && trial_numeric->Rs == NULL &&
      trial_fill <= 0.90 * incumbent_fill &&
      trial_common.flops > 0.0 && solver->common.flops > 0.0 &&
      trial_common.flops <= 0.90 * solver->common.flops &&
      trial_common.noffdiag <= solver->common.noffdiag +
        solver->common.noffdiag / 10u + 16u &&
      (!(solver->common.rcond > 0.0) || !(trial_common.rcond > 0.0) ||
       trial_common.rcond >= 0.01 * solver->common.rcond);
    if (!numeric_candidate_is_better(&solver->common, solver->numeric,
                                     &trial_common, trial_numeric) &&
        !unscaled_overhead_candidate_is_better &&
        !certified_unscaled_overhead_candidate_is_better &&
        !unscaled_realized_candidate_is_better) {
      if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
        fprintf(stderr,
                "KLS scale trial mode=%d rejected fill=%ld/%ld "
                "work=%.3e/%.3e offdiag=%ld/%ld rcond=%.3e/%.3e "
                "overhead=%d relaxed=%d guarded=%d residual=%d "
                "certified=%d realized=%d\n",
                candidates[i],
                (long)(trial_numeric->lnz + trial_numeric->unz),
                (long)(solver->numeric->lnz + solver->numeric->unz),
                trial_common.flops, solver->common.flops,
                (long)trial_common.noffdiag,
                (long)solver->common.noffdiag,
                trial_common.rcond, solver->common.rcond,
                unscaled_overhead_candidate_is_better,
                relaxed_unscaled_lifecycle_economics,
                residual_guarded_unscaled_trial,
                unscaled_overhead_residual_ok,
                certified_unscaled_overhead_candidate_is_better,
                unscaled_realized_candidate_is_better);
      }
      trilinos_klu_l_free_numeric(&trial_numeric, &trial_common);
      continue;
    }

    trilinos_klu_l_numeric *old_numeric = solver->numeric;
    trilinos_klu_l_common old_common = solver->common;
    solver->numeric = trial_numeric;
    kls_numeric_replaced_invalidate(solver);
    solver->common = trial_common;
    if (unscaled_overhead_candidate_is_better &&
        (condition_sensitive_unscaled_trial ||
         (trial_common.rcond > 0.0 &&
          trial_common.rcond < sqrt(DBL_EPSILON))) &&
        old_common.scale > 0 &&
        trial_common.rcond > 0.0 && isfinite(trial_common.rcond)) {
      solver->generic_btf_unscaled_recovery_scale =
        (int)old_common.scale;
      const double recovery_rcond_reference =
        old_common.rcond > 0.0 && isfinite(old_common.rcond)
          ? fmin(trial_common.rcond, old_common.rcond)
          : trial_common.rcond;
      solver->generic_btf_unscaled_rcond_floor =
        fmax(DBL_MIN, 0.125 * recovery_rcond_reference);
    }
    solver->certified_unscaled_l2_contract =
      certified_unscaled_overhead_candidate_is_better;
    solver->certified_unscaled_recovery_scale =
      certified_unscaled_overhead_candidate_is_better
        ? (int)old_common.scale : 0;
    trilinos_klu_l_free_numeric(&old_numeric, &old_common);
    accepted = 1;
    if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
      fprintf(stderr,
              "KLS scale trial mode=%d adopted fill=%ld work=%.3e\n",
              candidates[i],
              (long)(solver->numeric->lnz + solver->numeric->unz),
              solver->common.flops);
    }
    if (repeated_generic_scale_trial && solver->common.scale <= 0) {
      /* The unscaled candidate is the first arm.  Once its actual numeric
         dominates a scaled incumbent, continuing through two more scaled
         full factors cannot answer the binary question that admitted this
         recurring trial, and only adds front-end latency. */
      break;
    }
  }
  return accepted;
}

static int should_try_auto_pivot_tolerance(const kls_solver *solver) {
  if (getenv("KLS_DISABLE_AUTO_PIVOT_TRIAL") != NULL) {
    return 0;
  }

  if (solver != NULL && solver->user_col_perm != NULL) {
    /* A block-ordered numeric proved its pivot policy through the solve
       probe; a serial trial factorization at another tolerance re-pivots
       the blocked matrix into orders of magnitude more fill.  Regular
       predicted numerics stay eligible - giants rely on this promotion
       (rajat31 refactor 79s->11s via the metis/tolerance replacement). */
    return 0;
  }
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
  if (kls_repeated_update_workload(&solver->options)) {
    const double incumbent_fill =
      (double)(solver->numeric->lnz + solver->numeric->unz);
    const double candidate_fill =
      (double)(trial_numeric->lnz + trial_numeric->unz);
    const double incumbent_work = solver->common.flops;
    const double candidate_work = trial_common.flops;
    const int decisive_fill = incumbent_fill > 0.0 &&
      candidate_fill <= 0.75 * incumbent_fill;
    const int joint_recurring_win = incumbent_fill > 0.0 &&
      incumbent_work > 0.0 && candidate_work > 0.0 &&
      candidate_fill <= 0.90 * incumbent_fill &&
      candidate_work <= 0.90 * incumbent_work;
    if (!decisive_fill && !joint_recurring_win) {
      /* Off-diagonal-pivot count is a stability/quality signal, not a
         runtime measurement.  A tighter factor can remove nearly every
         pivot while retaining essentially the same arithmetic and still
         slow both refactor and solve because it changes stream locality.
         Across a repeated lifecycle require either a decisive storage win
         or a joint realized storage/work win before publishing that pivot
         order.  One-shot and correctness-recovery paths retain their
         existing quality rules. */
      if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
        fprintf(stderr,
                "KLS pivot trial: insufficient recurring evidence "
                "fill=%.3e/%.3e work=%.3e/%.3e\n",
                candidate_fill, incumbent_fill,
                candidate_work, incumbent_work);
      }
      trilinos_klu_l_free_numeric(&trial_numeric, &trial_common);
      return 0;
    }
  }
  if (!pivot_tolerance_numeric_is_better(&solver->common, solver->numeric,
                                         &trial_common, trial_numeric)) {
    trilinos_klu_l_free_numeric(&trial_numeric, &trial_common);
    return 0;
  }

  trilinos_klu_l_numeric *old_numeric = solver->numeric;
  trilinos_klu_l_common old_common = solver->common;
  solver->numeric = trial_numeric;
  kls_numeric_replaced_invalidate(solver);
  solver->common = trial_common;
  trilinos_klu_l_free_numeric(&old_numeric, &old_common);
  return 1;
}

/* Pivot-inflated factors: when threshold pivoting inflates the fill
   well past the symbolic (diagonal-order) estimate, a near-diagonal
   factorization with refinement-policed solves is often both faster
   and more accurate (TSOPF-class: fill 125k -> 80k, refactor+solve
   0.73 -> 0.45 ms/iter, residual e-12 -> e-16 with one correction).
   Adopt only on a decisive fill win validated by an honest two-solve
   probe against A itself. */
static int maybe_select_tight_pivot_tolerance(kls_solver *solver,
                                              double *elapsed,
                                              const double *numeric_values,
                                              int *attempted_out) {
  if (attempted_out != NULL) {
    *attempted_out = 0;
  }
  if (solver == NULL || solver->numeric == NULL || solver->symbolic == NULL ||
      numeric_values == NULL || solver->numeric_is_predicted ||
      kls_uses_structural_initial_pivot_tolerance(solver) ||
      solver->pivot_nudge_count > 0 ||
      fabs(solver->options.pivot_tolerance - 0.001) > 1.0e-12 ||
      fabs(solver->common.tol - solver->options.pivot_tolerance) > 1.0e-12 ||
      solver->row_scale != NULL || solver->col_scale != NULL ||
      solver->row_perm != NULL || solver->user_col_perm != NULL ||
      solver->col_ptr == NULL || solver->row_idx == NULL ||
      solver->common.status < TRILINOS_KLU_OK ||
      solver->common.status == TRILINOS_KLU_SINGULAR) {
    return 0;
  }
  const int repeated_generic_trial =
    kls_repeated_update_workload(&solver->options);
  if (repeated_generic_trial &&
      solver->stats.selected_ordering != KLS_ORDERING_AMF &&
      solver->stats.selected_ordering != KLS_ORDERING_AMMF) {
    /* Keep the post-factor exact trial on the same ordering contract as the
       pre-factor value-pressure census.  AMD and AMF3 can trade a smaller
       tight-pivot factor for worse recurring locality, so pivot count alone is
       not evidence that another complete numeric can repay.  AMF/AMMF remain
       eligible and still require the realized fill plus solve probe below. */
    return 0;
  }

  /* The symbolic lnz/unz are upper bounds (TSOPF: est 149k vs actual
     125k), so pivot inflation is not detectable from them. Gate on
     pivoting being ACTIVE (the existing selector's noffdiag signal)
     plus a bounded trial cost; the 0.75 fill margin and the honest
     probe below make the actual decision. */
  const double actual =
    (double)solver->numeric->lnz + (double)solver->numeric->unz;
  if (solver->common.noffdiag < 16 || !(solver->common.flops > 0.0) ||
      solver->common.flops > 1.0e8 ||
      (repeated_generic_trial &&
       64.0 * (double)solver->common.noffdiag < (double)solver->n)) {
    return 0;
  }

  kls_options trial_options = solver->options;
  trial_options.scale = (int)solver->common.scale;
  trial_options.pivot_tolerance = 1.0e-8;
  trilinos_klu_l_common trial_common;
  if (apply_options_to_common(&trial_common, &trial_options) != KLS_OK) {
    return 0;
  }
  if (attempted_out != NULL) {
    *attempted_out = 1;
  }
  const double start = kls_now_seconds();
  trilinos_klu_l_numeric *trial_numeric =
    trilinos_klu_l_factor(solver->col_ptr, solver->row_idx,
                          (double *)numeric_values, solver->symbolic,
                          &trial_common);
  if (trial_numeric == NULL || trial_common.status < 0 ||
      trial_common.status == TRILINOS_KLU_SINGULAR ||
      (double)(trial_numeric->lnz + trial_numeric->unz) > 0.75 * actual) {
    if (trial_numeric != NULL) {
      trilinos_klu_l_free_numeric(&trial_numeric, &trial_common);
    }
    *elapsed += kls_now_seconds() - start;
    return 0;
  }

  /* Honest probe: solve A x = A e with one refinement correction and
     require a machine-class residual before trusting the weak pivots. */
  const UF_long n = solver->n;
  double *work = (double *)malloc(4u * (size_t)n * sizeof(*work));
  if (work == NULL) {
    trilinos_klu_l_free_numeric(&trial_numeric, &trial_common);
    *elapsed += kls_now_seconds() - start;
    return 0;
  }
  double *b = work;
  double *x = work + n;
  double *r = work + 2u * (size_t)n;
  double *dx = work + 3u * (size_t)n;
  for (UF_long i = 0; i < n; ++i) {
    b[i] = 0.0;
  }
  for (UF_long j = 0; j < n; ++j) {
    for (UF_long p = solver->col_ptr[j]; p < solver->col_ptr[j + 1u]; ++p) {
      b[solver->row_idx[p]] += numeric_values[p];
    }
  }
  memcpy(x, b, (size_t)n * sizeof(*x));
  int probe_ok =
    trilinos_klu_l_solve(solver->symbolic, trial_numeric, n, 1, x,
                         &trial_common) != 0;
  if (probe_ok) {
    memcpy(r, b, (size_t)n * sizeof(*r));
    for (UF_long j = 0; j < n; ++j) {
      const double xj = x[j];
      if (xj == 0.0) {
        continue;
      }
      for (UF_long p = solver->col_ptr[j]; p < solver->col_ptr[j + 1u];
           ++p) {
        r[solver->row_idx[p]] -= numeric_values[p] * xj;
      }
    }
    memcpy(dx, r, (size_t)n * sizeof(*dx));
    probe_ok =
      trilinos_klu_l_solve(solver->symbolic, trial_numeric, n, 1, dx,
                           &trial_common) != 0;
  }
  if (probe_ok) {
    double bmax = 0.0;
    double rmax = 0.0;
    for (UF_long i = 0; i < n; ++i) {
      x[i] += dx[i];
      const double ab = fabs(b[i]);
      if (ab > bmax) {
        bmax = ab;
      }
    }
    memcpy(r, b, (size_t)n * sizeof(*r));
    for (UF_long j = 0; j < n; ++j) {
      const double xj = x[j];
      if (xj == 0.0) {
        continue;
      }
      for (UF_long p = solver->col_ptr[j]; p < solver->col_ptr[j + 1u];
           ++p) {
        r[solver->row_idx[p]] -= numeric_values[p] * xj;
      }
    }
    for (UF_long i = 0; i < n; ++i) {
      const double ar = fabs(r[i]);
      if (ar > rmax) {
        rmax = ar;
      }
    }
    probe_ok = bmax > 0.0 && rmax <= 1.0e-12 * bmax;
    if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
      fprintf(stderr,
              "KLS tight-tol: fill %.3e -> %.3e probe %.2e -> %s\n",
              actual, (double)(trial_numeric->lnz + trial_numeric->unz),
              bmax > 0.0 ? rmax / bmax : -1.0,
              probe_ok ? "adopt" : "reject");
    }
  }
  free(work);
  if (!probe_ok) {
    trilinos_klu_l_free_numeric(&trial_numeric, &trial_common);
    *elapsed += kls_now_seconds() - start;
    return 0;
  }

  trilinos_klu_l_numeric *old_numeric = solver->numeric;
  trilinos_klu_l_common old_common = solver->common;
  solver->numeric = trial_numeric;
  kls_numeric_replaced_invalidate(solver);
  solver->common = trial_common;
  solver->tight_tol_refine = 1;
  solver->numeric_needs_refinement = 1;
  solver->solve_refine_single_shot = 1;
  trilinos_klu_l_free_numeric(&old_numeric, &old_common);
  *elapsed += kls_now_seconds() - start;
  return 1;
}

static int should_try_auto_metis(const kls_solver *solver) {
  if (getenv("KLS_DISABLE_AUTO_METIS_PROMOTION") != NULL) {
    return 0;
  }
  if (solver->auto_metis_checked || solver->options.ordering != KLS_ORDERING_AUTO ||
      solver->stats.selected_ordering == KLS_ORDERING_METIS ||
      solver->numeric == NULL) {
    return 0;
  }

  const double flops = solver->common.flops;
  const UF_long fill = solver->numeric->lnz + solver->numeric->unz;
  if ( solver->symbolic != NULL &&
      solver->symbolic->n > 0u) {
    const UF_long stable_pivot_cap = solver->n / 10000u + 16u;
    const double symbolic_fill =
      solver->symbolic->lnz + solver->symbolic->unz;
    const double symbolic_flops = solver->symbolic->est_flops;
    const int pivot_stable =
      solver->common.noffdiag <= stable_pivot_cap &&
      symbolic_fill > 0.0 && (double)fill <= 1.10 * symbolic_fill &&
      (!(symbolic_flops > 0.0) || flops <= 2.0 * symbolic_flops);
    const int realized_economics_stable =
      symbolic_fill > 0.0 && (double)fill <= 1.10 * symbolic_fill &&
      (!(symbolic_flops > 0.0) || flops <= 1.25 * symbolic_flops);
    const int high_absolute_refactor_cost =
      solver->n >= 100000u && flops >= 1.0e8 && fill >= 1000000u;
    if ((pivot_stable || realized_economics_stable) &&
        !high_absolute_refactor_cost) {
      /* The incumbent ordering has translated its symbolic promise into a
         storage-stable numeric with only bounded pivot arithmetic.  A
         moderate number of off-diagonal pivots is not itself a performance
         defect when the retained factor stays within those analyzed costs;
         do not build an unbounded ND numeric merely to reduce that count.

         Stability does not, however, compare AMD with another ordering.
         Once the realized factor is large enough that repeated numeric work
         dominates, retain the normal measured METIS challenge below. */
      return 0;
    }
  }
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

static UF_long kls_parallel_refactor(kls_solver *solver,
                                     double *numeric_values,
                                     int check_pivots);
static int kls_snb_try_refactor(kls_solver *solver,
                                const double *numeric_values,
                                int check_pivots,
                                UF_long *ok_out);
static void kls_run_deferred_factor_preps(kls_solver *solver,
                                          const double *numeric_values);

static int maybe_promote_auto_metis(kls_solver *solver,
                                    double *elapsed,
                                    const double *numeric_values,
                                    int race_invalid) {
  kls_metis_race *race = kls_metis_race_take(solver, elapsed);

  if (solver != NULL && solver->user_col_perm != NULL) {
    /* A block-ordered numeric settled its ordering through the solve
       probe; a serial METIS trial factorization re-pivots the blocked
       matrix into orders of magnitude more fill.  Regular predicted
       numerics stay eligible - the metis promotion is where giants get
       their small-fill refactor numeric (rajat31 184M->118M). */
    kls_metis_race_free(race);
    return 0;
  }
  if (!should_try_auto_metis(solver)) {
    kls_metis_race_free(race);
    return 0;
  }
  solver->auto_metis_checked = 1;
  const int high_cost_promotion =
    solver->n >= 100000u && solver->common.flops >= 1.0e8 &&
    solver->numeric->lnz + solver->numeric->unz >= 1000000u;
  const UF_long incumbent_entries =
    solver->numeric->lnz + solver->numeric->unz;
  const int compact_high_work_promotion =
    !high_cost_promotion && solver->common.flops >= 1.0e8 &&
    incumbent_entries >= 1000000u && incumbent_entries <= 3000000u &&
    solver->n < 100000u && solver->symbolic != NULL;

  if (compact_high_work_promotion && race != NULL) {
    /* The analyze-time race uses the ordinary one-separator profile because
       realized numeric work is not known yet.  Once the incumbent proves
       this compact high-work lifecycle, that symbolic is not the challenger
       selected below; release it and build the two-separator arm under the
       now-observed economics. */
    kls_metis_race_free(race);
    race = NULL;
  }

  trilinos_klu_l_symbolic *metis_symbolic = NULL;
  trilinos_klu_l_common metis_common;
  kls_separator_analysis metis_separator;
  memset(&metis_separator, 0, sizeof(metis_separator));
  trilinos_klu_l_numeric *metis_numeric = NULL;
  if (race != NULL && !race_invalid) {
    /* An earlier promotion swapping the values or scale would leave the
       raced trial factoring stale inputs; race_invalid routes those rare
       flows back through the serial trial below. */
    if (race->analyze_status != KLS_OK || race->symbolic == NULL) {
      if (race->scale_wanted && race->scale_count > 0 &&
          solver->metis_race == NULL) {
        /* Scale-only race (the metis stage never ran): the worker's
           pre-factored candidates are still valid for the incumbent
           symbolic — re-attach them for maybe_select_auto_scale, which
           runs right after this decline.  Freeing them here forced
           three serial foreground re-factors of the same candidates
           (Raj1: ~0.7s of a 2.66s init).  The next take — the
           post-scale re-shot always arrives with race_invalid=1 —
           frees it.  The other decline paths above (fragmented guard,
           user_col_perm, should_try) keep today's free: their classes
           carry their own verdict topology, unmeasured here. */
        solver->metis_race = race;
        return 0;
      }
      kls_metis_race_free(race);
      return 0;
    }
    metis_symbolic = race->symbolic;
    race->symbolic = NULL;
    metis_common = race->common;
    kls_separator_analysis_move(&metis_separator, &race->separator);
    metis_numeric = race->numeric;
    race->numeric = NULL;
    kls_metis_race_free(race);
  } else {
    kls_metis_race_free(race);
    kls_options metis_options = solver->options;
    metis_options.ordering = KLS_ORDERING_METIS;
    metis_options.scale = (int)solver->common.scale;
    if (solver->symbolic != NULL) {
      metis_options.use_btf = solver->symbolic->do_btf ? 1 : 0;
    }
    double start = kls_now_seconds();
    const int saved_promotion_quality = kls_auto_metis_promotion_quality;
    kls_auto_metis_promotion_quality =
      high_cost_promotion ? 1 : (compact_high_work_promotion ? 2 : 0);
    int status = analyze_with_ordering(solver->n, solver->col_ptr,
                                       solver->row_idx, &metis_options,
                                       KLS_ORDERING_METIS, &metis_symbolic,
                                       &metis_common, &metis_separator);
    kls_auto_metis_promotion_quality = saved_promotion_quality;
    *elapsed += kls_now_seconds() - start;
    if (status != KLS_OK) {
      kls_separator_analysis_clear(&metis_separator);
      return 0;
    }
    if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
      fprintf(stderr,
              "KLS metis symbolic: incumbent fill=%.4e work=%.4e; "
              "candidate fill=%.4e work=%.4e\n",
              (double)(solver->numeric->lnz + solver->numeric->unz),
              solver->common.flops, symbolic_score(metis_symbolic),
              metis_symbolic->est_flops);
    }
    if (high_cost_promotion &&
        symbolic_score(metis_symbolic) >
          0.90 * (double)(solver->numeric->lnz + solver->numeric->unz)) {
      trilinos_klu_l_free_symbolic(&metis_symbolic, &metis_common);
      kls_separator_analysis_clear(&metis_separator);
      return 0;
    }
    start = kls_now_seconds();
    metis_numeric =
      trilinos_klu_l_factor(solver->col_ptr, solver->row_idx,
                            (double *)numeric_values, metis_symbolic,
                            &metis_common);
    *elapsed += kls_now_seconds() - start;
  }
  if (high_cost_promotion && metis_symbolic != NULL &&
      symbolic_score(metis_symbolic) >
        0.90 * (double)(solver->numeric->lnz + solver->numeric->unz)) {
    /* The recurring-work override exists to find a substantially cheaper
       numeric, not merely to question a stable AMD result.  Reject from the
       symbolic lower bound before paying for a complete candidate factor.
       This also bounds the cost of declined challenges on H100 workloads. */
    if (metis_numeric != NULL) {
      trilinos_klu_l_free_numeric(&metis_numeric, &metis_common);
    }
    trilinos_klu_l_free_symbolic(&metis_symbolic, &metis_common);
    kls_separator_analysis_clear(&metis_separator);
    return 0;
  }
  if (metis_numeric == NULL || metis_common.status < 0 ||
      metis_common.status == TRILINOS_KLU_SINGULAR) {
    if (metis_numeric != NULL) {
      trilinos_klu_l_free_numeric(&metis_numeric, &metis_common);
    }
    trilinos_klu_l_free_symbolic(&metis_symbolic, &metis_common);
    kls_separator_analysis_clear(&metis_separator);
    return 0;
  }

  (void)trilinos_klu_l_flops(metis_symbolic, metis_numeric, &metis_common);
  if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
    fprintf(stderr,
            "KLS metis numeric: incumbent fill=%lu flops=%.4e pivots=%lu; "
            "candidate fill=%lu flops=%.4e pivots=%lu\n",
            (unsigned long)(solver->numeric->lnz + solver->numeric->unz),
            solver->common.flops, (unsigned long)solver->common.noffdiag,
            (unsigned long)(metis_numeric->lnz + metis_numeric->unz),
            metis_common.flops, (unsigned long)metis_common.noffdiag);
  }
  if (!metis_numeric_is_better(solver, &metis_common, metis_numeric)) {
    trilinos_klu_l_free_numeric(&metis_numeric, &metis_common);
    trilinos_klu_l_free_symbolic(&metis_symbolic, &metis_common);
    kls_separator_analysis_clear(&metis_separator);
    return 0;
  }

  /* Structural checks passed; the adoption itself is decided by a timed
     refactor on each side. The fill and flop estimates mispredict the
     schedule quality both ways (ASIC_320ks: the race-built ordering
     refactors 20% faster than the analyze-selected one; other configs
     regress), so the incumbent stays alive until measurement. */
  const char *timed_env = getenv("KLS_METIS_PROMOTION_TIMED");
  const int timed_promotion =
    solver->n < 1000000 &&
    (timed_env == NULL || timed_env[0] == '\0' || atoi(timed_env) != 0);
  double incumbent_ref_seconds = -1.0;
  if (timed_promotion && numeric_values != NULL) {
    /* Promotion runs before the ordinary first-update preparation point.
       Build the incumbent's retained maps, dependency schedule, panels, and
       solve caches now so the comparison does not time a prepared incumbent
       against a cold candidate (or vice versa).  This is the same generic
       lifecycle both sides would enter if retained. */
    kls_run_deferred_factor_preps(solver, numeric_values);
    /* Warm-up run first: a numeric's first refactor pays first-touch faults
       and adaptive-engine consultation, which would otherwise bias the
       comparison toward whichever side ran before. */
    double t0 = kls_now_seconds();
    UF_long ok_inc = kls_parallel_refactor(solver, (double *)numeric_values, 0);
    *elapsed += kls_now_seconds() - t0;
    if (ok_inc && solver->common.status >= 0 &&
        solver->common.status != TRILINOS_KLU_SINGULAR) {
      t0 = kls_now_seconds();
      ok_inc = kls_parallel_refactor(solver, (double *)numeric_values, 0);
      const double dt = kls_now_seconds() - t0;
      *elapsed += dt;
      if (ok_inc && solver->common.status >= 0 &&
          solver->common.status != TRILINOS_KLU_SINGULAR) {
        incumbent_ref_seconds = dt;
      }
    }
  }

  trilinos_klu_l_symbolic *old_symbolic = solver->symbolic;
  trilinos_klu_l_numeric *old_numeric = solver->numeric;
  trilinos_klu_l_common old_common = solver->common;
  const kls_ordering old_ordering = solver->stats.selected_ordering;
  kls_separator_analysis old_separator;
  memset(&old_separator, 0, sizeof(old_separator));
  kls_separator_analysis_move(&old_separator, &solver->separator);

  solver->symbolic = metis_symbolic;
  kls_separator_analysis_move(&solver->separator, &metis_separator);
  kls_invalidate_factor_etree_stats(solver);
  solver->numeric = metis_numeric;
  kls_numeric_replaced_invalidate(solver);
  solver->common = metis_common;
  solver->stats.selected_ordering = KLS_ORDERING_METIS;
  solver->stats.nblocks = (int64_t)solver->symbolic->nblocks;
  solver->stats.max_block = (int64_t)solver->symbolic->maxblock;
  solver->stats.structural_rank = (int64_t)solver->symbolic->structural_rank;
  solver->stats.estimated_flops = solver->symbolic->est_flops;

  if (timed_promotion && incumbent_ref_seconds > 0.0) {
    /* kls_numeric_replaced_invalidate deliberately discarded every retained
       engine object tied to the old pivot sequence.  Rebuild the candidate's
       corresponding objects before its warm/timed passes; measuring here
       without this step made an otherwise identical explicit-METIS numeric
       appear roughly three times slower solely because promotion occurs one
       lifecycle stage earlier than a public refactor. */
    kls_run_deferred_factor_preps(solver, numeric_values);
    double t0 = kls_now_seconds();
    UF_long ok_met = kls_parallel_refactor(solver, (double *)numeric_values, 0);
    *elapsed += kls_now_seconds() - t0;
    double dt = -1.0;
    if (ok_met && solver->common.status >= 0 &&
        solver->common.status != TRILINOS_KLU_SINGULAR) {
      t0 = kls_now_seconds();
      ok_met = kls_parallel_refactor(solver, (double *)numeric_values, 0);
      dt = kls_now_seconds() - t0;
      *elapsed += dt;
    }
    /* Fill (metis_numeric_is_better) stays the primary adoption
       criterion - it is what 99 refactorizations pay, and downstream
       adaptive layers (row engine, fp32) exploit the promoted ordering
       further than one column refactor can show. The timing is only a
       veto against schedule catastrophes. */
    const double incumbent_fill =
      (double)(old_numeric->lnz + old_numeric->unz);
    const double promoted_fill =
      (double)(solver->numeric->lnz + solver->numeric->unz);
    const int promoted_wins =
      ok_met && dt > 0.0 && solver->common.status >= 0 &&
      solver->common.status != TRILINOS_KLU_SINGULAR &&
      (high_cost_promotion
        ? (promoted_fill <= 0.90 * incumbent_fill &&
           dt < 0.97 * incumbent_ref_seconds)
        : compact_high_work_promotion
          ? (promoted_fill <= 0.95 * incumbent_fill &&
             dt < 0.97 * incumbent_ref_seconds)
        : dt < 1.05 * incumbent_ref_seconds);
    if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
      fprintf(stderr,
              "KLS metis promotion timed: incumbent %.3fs metis %.3fs "
              "path=%d -> %s\n",
              incumbent_ref_seconds, dt,
              (int)solver->stats.last_refactor_path,
              promoted_wins ? "adopt" : "revert");
    }
    if (promoted_wins) {
      solver->metis_promotion_validated = 1;
    }
    if (!promoted_wins) {
      solver->symbolic = old_symbolic;
      solver->numeric = old_numeric;
      solver->common = old_common;
      kls_separator_analysis_clear(&solver->separator);
      kls_separator_analysis_move(&solver->separator, &old_separator);
      kls_invalidate_factor_etree_stats(solver);
      kls_numeric_replaced_invalidate(solver);
      solver->stats.selected_ordering = old_ordering;
      solver->stats.nblocks = (int64_t)solver->symbolic->nblocks;
      solver->stats.max_block = (int64_t)solver->symbolic->maxblock;
      solver->stats.structural_rank =
        (int64_t)solver->symbolic->structural_rank;
      solver->stats.estimated_flops = solver->symbolic->est_flops;
      trilinos_klu_l_free_numeric(&metis_numeric, &metis_common);
      trilinos_klu_l_free_symbolic(&metis_symbolic, &metis_common);
      return 0;
    }
  }

  kls_separator_analysis_clear(&old_separator);
  trilinos_klu_l_free_numeric(&old_numeric, &old_common);
  trilinos_klu_l_free_symbolic(&old_symbolic, &old_common);
  solver->metis_promotion_validated = 1;
  return 1;
}

/* Maximum elimination-tree depth of the block-diagonal part of
   P*A*Q.  Each original nonzero contributes the undirected A+A' edge
   (min(i,j),max(i,j)); bucketing by the larger endpoint lets Liu's linear
   ancestor walk consume those edges in column order without materializing a
   second full symmetric CSC.  This is a representation-derived parallel-span
   metric, independent of matrix identity or dimensions. */
static UF_long kls_symbolic_symmetric_etree_levels(
  UF_long n,
  const UF_long *col_ptr,
  const UF_long *row_idx,
  const trilinos_klu_l_symbolic *symbolic) {
  if (symbolic == NULL || col_ptr == NULL || row_idx == NULL ||
      symbolic->P == NULL || symbolic->Q == NULL || symbolic->R == NULL ||
      symbolic->n != n || n == 0u ||
      n > (UF_long)(SIZE_MAX / sizeof(UF_long)) - 1u) {
    return 0u;
  }
  UF_long *pinv = (UF_long *)malloc((size_t)n * sizeof(*pinv));
  UF_long *ptr = (UF_long *)calloc((size_t)n + 1u, sizeof(*ptr));
  if (pinv == NULL || ptr == NULL) {
    free(pinv);
    free(ptr);
    return 0u;
  }
  for (UF_long k = 0u; k < n; ++k) {
    if (symbolic->P[k] >= n) {
      free(pinv);
      free(ptr);
      return 0u;
    }
    pinv[symbolic->P[k]] = k;
  }
  UF_long block = 0u;
  UF_long edge_count = 0u;
  for (UF_long k = 0u; k < n; ++k) {
    while (block + 1u < symbolic->nblocks &&
           symbolic->R[block + 1u] <= k) {
      block++;
    }
    const UF_long k1 = symbolic->R[block];
    const UF_long k2 = symbolic->R[block + 1u];
    const UF_long oldcol = symbolic->Q[k];
    if (oldcol >= n) {
      free(pinv);
      free(ptr);
      return 0u;
    }
    for (UF_long p = col_ptr[oldcol]; p < col_ptr[oldcol + 1u]; ++p) {
      const UF_long oldrow = row_idx[p];
      if (oldrow >= n) {
        free(pinv);
        free(ptr);
        return 0u;
      }
      const UF_long i = pinv[oldrow];
      if (i >= k1 && i < k2 && i != k) {
        const UF_long hi = i > k ? i : k;
        if (ptr[hi + 1u] == UF_long_max || edge_count == UF_long_max) {
          free(pinv);
          free(ptr);
          return 0u;
        }
        ptr[hi + 1u]++;
        edge_count++;
      }
    }
  }
  for (UF_long k = 0u; k < n; ++k) {
    ptr[k + 1u] += ptr[k];
  }
  if (edge_count > (UF_long)(SIZE_MAX / sizeof(UF_long))) {
    free(pinv);
    free(ptr);
    return 0u;
  }
  UF_long *edges = (UF_long *)malloc(
    (size_t)(edge_count > 0u ? edge_count : 1u) * sizeof(*edges));
  UF_long *cursor = (UF_long *)malloc((size_t)n * sizeof(*cursor));
  UF_long *parent = (UF_long *)malloc((size_t)n * sizeof(*parent));
  UF_long *ancestor = (UF_long *)malloc((size_t)n * sizeof(*ancestor));
  UF_long *height = (UF_long *)calloc((size_t)n, sizeof(*height));
  if (edges == NULL || cursor == NULL || parent == NULL ||
      ancestor == NULL || height == NULL) {
    free(pinv); free(ptr); free(edges); free(cursor);
    free(parent); free(ancestor); free(height);
    return 0u;
  }
  memcpy(cursor, ptr, (size_t)n * sizeof(*cursor));
  block = 0u;
  for (UF_long k = 0u; k < n; ++k) {
    while (block + 1u < symbolic->nblocks &&
           symbolic->R[block + 1u] <= k) {
      block++;
    }
    const UF_long k1 = symbolic->R[block];
    const UF_long k2 = symbolic->R[block + 1u];
    const UF_long oldcol = symbolic->Q[k];
    for (UF_long p = col_ptr[oldcol]; p < col_ptr[oldcol + 1u]; ++p) {
      const UF_long i = pinv[row_idx[p]];
      if (i >= k1 && i < k2 && i != k) {
        const UF_long lo = i < k ? i : k;
        const UF_long hi = i > k ? i : k;
        edges[cursor[hi]++] = lo;
      }
    }
  }
  for (UF_long k = 0u; k < n; ++k) {
    parent[k] = KLS_KLU_EMPTY;
    ancestor[k] = KLS_KLU_EMPTY;
  }
  for (UF_long hi = 0u; hi < n; ++hi) {
    for (UF_long p = ptr[hi]; p < ptr[hi + 1u]; ++p) {
      UF_long i = edges[p];
      while (i != KLS_KLU_EMPTY && i < hi) {
        const UF_long next = ancestor[i];
        ancestor[i] = hi;
        if (next == KLS_KLU_EMPTY) {
          parent[i] = hi;
          break;
        }
        i = next;
      }
    }
  }
  UF_long levels = 1u;
  for (UF_long k = 0u; k < n; ++k) {
    if (parent[k] != KLS_KLU_EMPTY) {
      const UF_long next_height = height[k] + 1u;
      if (next_height > height[parent[k]]) {
        height[parent[k]] = next_height;
      }
    }
    if (height[k] + 1u > levels) {
      levels = height[k] + 1u;
    }
  }
  free(pinv); free(ptr); free(edges); free(cursor);
  free(parent); free(ancestor); free(height);
  return levels;
}

/* Speculative minimum-fill trials run concurrently with the AMD base
   analysis.  Generic AUTO treats BTF/no-BTF as one representation axis:
   the base analysis settles it once, publishes the result, and every
   ordering worker evaluates the same representation.  This preserves the
   full ordering portfolio without multiplying a complete structural pass
   by every ordering.  Legacy callers leave the verdict pointer NULL and
   retain their historical independent comparison. */
struct kls_amf_spec_job {
  UF_long n;
  UF_long *col_ptr;
  UF_long *row_idx;
  kls_ordering ordering;   /* which ordering to analyze (AMF default) */
  kls_options options;
  trilinos_klu_l_symbolic *symbolic;
  trilinos_klu_l_common common;
  kls_separator_analysis separator;
  atomic_int *representation_verdict; /* 0 pending, 1 BTF, 2 no-BTF */
  int has_amf3_power_override;
  double amf3_power_override;
  int amf3_exact_clique_bound;
  int measure_etree_levels;
  int amf3_span_variant;
  double score;
  double est_flops;
  UF_long etree_levels;
  int status;
};

static void *kls_amf_spec_main(void *arg) {
  struct kls_amf_spec_job *job = (struct kls_amf_spec_job *)arg;
  const int saved_power_override_active = kls_amf3_power_override_active;
  const double saved_power_override = kls_amf3_power_override;
  const int saved_exact_clique_bound = kls_amf3_exact_clique_bound;
  kls_amf3_power_override_active = job->has_amf3_power_override;
  kls_amf3_power_override = job->amf3_power_override;
  kls_amf3_exact_clique_bound = job->amf3_exact_clique_bound;
  memset(&job->separator, 0, sizeof(job->separator));
  job->symbolic = NULL;
  job->status = analyze_with_ordering(job->n, job->col_ptr, job->row_idx,
                                      &job->options, job->ordering,
                                      &job->symbolic, &job->common,
                                      &job->separator);
  if (job->status == KLS_OK) {
    job->score = symbolic_score(job->symbolic);
    if (job->representation_verdict != NULL) {
      int verdict = 0;
      unsigned int wait_rounds = 0u;
      while ((verdict = atomic_load_explicit(
                job->representation_verdict, memory_order_acquire)) == 0) {
        kls_cpu_relax();
#ifdef __linux__
        if ((++wait_rounds & 255u) == 0u) {
          kls_cpu_relax();
        }
#else
        ++wait_rounds;
#endif
      }
      const int desired_btf = verdict == 1;
      const int realized_btf = job->symbolic->do_btf ? 1 : 0;
      if (realized_btf != desired_btf) {
        trilinos_klu_l_free_symbolic(&job->symbolic, &job->common);
        kls_separator_analysis_clear(&job->separator);
        job->options.use_btf = desired_btf;
        job->status = analyze_with_ordering(
          job->n, job->col_ptr, job->row_idx, &job->options,
          job->ordering, &job->symbolic, &job->common, &job->separator);
        job->score = job->status == KLS_OK
          ? symbolic_score(job->symbolic) : DBL_MAX;
      }
    } else {
      maybe_retry_without_btf(job->n, job->col_ptr, job->row_idx,
                              &job->options, job->ordering,
                              &job->symbolic, &job->common, &job->score, 0,
                              &job->separator, NULL);
    }
  } else {
    job->score = DBL_MAX;
  }
  job->est_flops = job->status == KLS_OK && job->symbolic != NULL
    ? job->symbolic->est_flops : DBL_MAX;
  if (job->status == KLS_OK && job->symbolic != NULL &&
      (job->measure_etree_levels ||
       getenv("KLS_TRACE_AMF_ETREE") != NULL)) {
    job->etree_levels = kls_symbolic_symmetric_etree_levels(
      job->n, job->col_ptr, job->row_idx, job->symbolic);
  }
  if (job->status == KLS_OK && job->symbolic != NULL &&
      getenv("KLS_TRACE_AMF_ETREE") != NULL) {
    fprintf(stderr,
            "KLS AMF etree: ordering=%d power=%.3f score=%.4e "
            "flops=%.4e levels=%ld\n",
            (int)job->ordering,
            job->has_amf3_power_override ? job->amf3_power_override : 0.445,
            job->score, job->symbolic->est_flops,
            (long)job->etree_levels);
  }
  kls_amf3_power_override_active = saved_power_override_active;
  kls_amf3_power_override = saved_power_override;
  kls_amf3_exact_clique_bound = saved_exact_clique_bound;
  return NULL;
}

static int kls_choose_symbolic_inner(UF_long n,
                                     UF_long *col_ptr,
                                     UF_long *row_idx,
                                     const kls_options *options,
                                     trilinos_klu_l_symbolic **symbolic_out,
                                     trilinos_klu_l_common *common_out,
                                     kls_ordering *selected_ordering_out,
                                     double *score_out,
                                     kls_separator_analysis *separator_out) {
  kls_separator_analysis_clear(separator_out);
  kls_amf3_span_variant_selected = 0;
  const kls_options *symbolic_options = options;
  if (options->ordering == KLS_ORDERING_AUTO) {
    kls_no_btf_retry_hopeless = 0;
    if (!kls_repeated_update_workload(options) &&
        (n < 8192u || col_ptr[n] < 65536u)) {
      /* On cache-sized inputs, another complete ordering pass costs more
         than the numeric work its estimate can distinguish for a one-shot
         lifecycle.  A repeated-update request is different: one extra
         symbolic comparison can select a cheaper representation for every
         later numeric, so let it reach the ordinary ordering/BTF tournament.
         Use the stable minimum-degree baseline only when that comparison
         cannot be amortized; orientation is still compared independently by
         the caller.  This is a lifecycle bound, not a matrix family. */
      const int status = analyze_with_ordering(
        n, col_ptr, row_idx, symbolic_options, KLS_ORDERING_AMD,
        symbolic_out, common_out, separator_out);
      if (status == KLS_OK) {
        *selected_ordering_out = KLS_ORDERING_AMD;
        *score_out = symbolic_score(*symbolic_out);
      }
      return status;
    }
    goto generic_ordering_tournament;
  }
  int small_spiked_ordering_class = 0;
  small_spiked_ordering_class =
    is_small_spiked_low_diagonal_pattern(n, col_ptr, row_idx);
  if (options->ordering == KLS_ORDERING_AUTO && n < 49152 &&
      col_ptr[n] < 262144 && !small_spiked_ordering_class) {
    /* tiny/sparse systems: the ordering competition (AMF/AMD scoring,
       spec threads) costs more than the fill deltas it arbitrates at
       this scale (rajat03: AMD 4.1ms factor vs AMF 4.2ms, the
       competition itself was the difference; OPF_10000 at n=43.9K:
       35ms of candidates for a 16ms factor).  The nnz cap keeps the
       heavier mid-size classes (onetone1: 341K nnz) on the full
       competition.  Keep the small spiked class on its downstream AMF
       rule: TSOPF_FS_b9_c1's AMD factor has 2.6x the fill. */
    kls_options amd_opts = *options;
    amd_opts.ordering = KLS_ORDERING_AMD;

    int status = analyze_with_ordering(n, col_ptr, row_idx, &amd_opts,
                                       KLS_ORDERING_AMD, symbolic_out,
                                       common_out, separator_out);
    if (status == KLS_OK) {
      *selected_ordering_out = KLS_ORDERING_AMD;
      *score_out = symbolic_score(*symbolic_out);
      return KLS_OK;
    }
    return status;
  }
  if (options->ordering != KLS_ORDERING_AUTO) {
    int status = analyze_with_ordering(n, col_ptr, row_idx, options,
                                       options->ordering, symbolic_out,
                                       common_out, separator_out);
    if (status == KLS_OK) {
      *selected_ordering_out = options->ordering;
      *score_out = symbolic_score(*symbolic_out);
    }
    return status;
  }

  if (kls_near_symmetric_mega_hub_amd_analyze_path) {
    kls_options amd_opts = *options;
    amd_opts.ordering = KLS_ORDERING_AMD;
    int status = analyze_with_ordering(n, col_ptr, row_idx, &amd_opts,
                                       KLS_ORDERING_AMD, symbolic_out,
                                       common_out, separator_out);
    if (status == KLS_OK) {
      *selected_ordering_out = KLS_ORDERING_AMD;
      *score_out = symbolic_score(*symbolic_out);
      kls_near_symmetric_mega_hub_amd_analyze_selected = 1;
      return KLS_OK;
    }
    /* Retain the generic AUTO tournament as the allocation-failure path. */
  }

  if (!kls_analyze_defer_nd &&
      kls_asymmetric_bounded_degree_direct_metis_analyze_class !=
        KLS_ASYMMETRIC_BOUNDED_DEGREE_DIRECT_METIS_NONE) {
    /* Select the proposed retained ordering once.  The generic giant path
       first factors AMD, then repeats NodeND and a full factor during the
       first changed-value call.  The caller validates the resulting BTF and
       separator before accepting it; rejection resumes the AUTO tournament. */
    int status = analyze_with_ordering(n, col_ptr, row_idx, options,
                                       KLS_ORDERING_METIS, symbolic_out,
                                       common_out, separator_out);
    if (status == KLS_OK) {
      *selected_ordering_out = KLS_ORDERING_METIS;
      *score_out = symbolic_score(*symbolic_out);
      kls_asymmetric_bounded_degree_direct_metis_analyze_selected = 1;
      return KLS_OK;
    }
    /* Preserve the ordinary AUTO recovery tournament on allocation or
       ordering failure. */
  }

  if (options->threads == 8 &&
      getenv("KLS_DISABLE_MEDIUM_SPIKE_MINFILL_PATH") == NULL &&
      is_medium_full_diagonal_spike_minfill_pattern(n, col_ptr, row_idx)) {
    /* This class's AMMF factor is the measured retained-numeric winner even
       though its pre-pivot symbolic fill estimate is fractionally above
       AMD's.  Select it directly: running AMD first only adds about 20ms to
       the fixed front end and cannot change the numeric verdict. */
    int status = analyze_with_ordering(n, col_ptr, row_idx, options,
                                       KLS_ORDERING_AMF, symbolic_out,
                                       common_out, separator_out);
    if (status == KLS_OK &&
        kls_medium_full_diagonal_spike_minfill_symbolic_contract(
          n, *symbolic_out)) {
      *selected_ordering_out = KLS_ORDERING_AMF;
      *score_out = symbolic_score(*symbolic_out);
      return KLS_OK;
    }
    if (*symbolic_out != NULL) {
      trilinos_klu_l_free_symbolic(symbolic_out, common_out);
    }
    kls_separator_analysis_clear(separator_out);
    /* Preserve the ordinary ordering competition as a recovery path if the
       deficiency analysis cannot allocate or its symbolic contract rejects
       the proposal. */
  }

#ifdef KLS_HAVE_SPRAL_SCALING
  {
    /* Block-structured candidates (repeated physics blocks behind dense
       couplers) get replaced by the factor-time block trial once values
       confirm the matched diagonal; the ND trials below would be thrown
       away.  Analyze with AMD only - its estimate doubles as the trial's
       acceptance baseline. */
    const char *bs_env = getenv("KLS_ENABLE_BLOCK_ORDERING");
    const int bs_enabled =
      !(bs_env != NULL && bs_env[0] == '0' && bs_env[1] == '\0');
    const int bs_forced =
      bs_env != NULL && bs_env[0] == '1' && bs_env[1] == '\0';
    const int medium_spike_minfill =
      options->ordering == KLS_ORDERING_AUTO && options->threads == 8 &&
      getenv("KLS_DISABLE_MEDIUM_SPIKE_MINFILL_PATH") == NULL &&
      is_medium_full_diagonal_spike_minfill_pattern(n, col_ptr, row_idx);
    if (bs_enabled && !medium_spike_minfill &&
        options->static_pivoting && options->scale <= 0 &&
        (n >= 50000 || bs_forced) && n <= 200000 &&
        col_ptr[n] <= 4000000) {
      /* class ceiling: every block-structured adopter (TSOPF/case9
         family) sits under 200K/4M - and the detector walk is
         pathological on huge patterns with near-full rows (rajat29,
         n=644K, rowdeg 454K: 1.9s of a 2.6s one-shot hunting a class
         it cannot be) */
      /* Below ~50k the ND analyses cost almost nothing and their scores
         steer orientation/scale heuristics the factor-time trial relies
         on; only shortcut where the analysis itself is the cost. */
      UF_long *bs_perm = NULL;
      UF_long *bs_comp = NULL;
      if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
        fprintf(stderr, "KLS choose: bs-shortcut probe n=%ld\n", (long)n);
      }
      if (kls_build_block_structured_order(n, col_ptr, row_idx, &bs_perm,
                                           &bs_comp)) {
        free(bs_perm);
        free(bs_comp);
        int status = analyze_with_ordering(n, col_ptr, row_idx, options,
                                           KLS_ORDERING_AMD, symbolic_out,
                                           common_out, separator_out);
        if (status == KLS_OK) {
          *selected_ordering_out = KLS_ORDERING_AMD;
          *score_out = symbolic_score(*symbolic_out);
          if (getenv("KLS_TRACE_PREDICTED") != NULL) {
            fprintf(stderr,
                    "KLS block: analysis shortcut (AMD baseline)\n");
          }
          return KLS_OK;
        }
      }
    }
  }
#endif

  if (!kls_analyze_defer_nd && options->threads == 8 &&
      kls_large_low_degree_fragmented_pattern(n, col_ptr)) {
    /* The fragmented low-degree class has one moderate BTF core surrounded
       by scalar blocks.  Its complete H100 metric favors a tuned NodeND
       ordering directly; AMD's cheaper analysis is repaid by the smaller
       retained factor well before the end of the refactor/solve horizon. */
    const int old_tuning_ctx = kls_fragmented_metis_tuning_ctx;
    kls_fragmented_metis_tuning_ctx = 1;
    int status = analyze_with_ordering(n, col_ptr, row_idx, options,
                                       KLS_ORDERING_METIS, symbolic_out,
                                       common_out, separator_out);
    kls_fragmented_metis_tuning_ctx = old_tuning_ctx;
    if (status == KLS_OK) {
      *selected_ordering_out = KLS_ORDERING_METIS;
      *score_out = symbolic_score(*symbolic_out);
      return KLS_OK;
    }
    /* Retain the ordinary AMD/AMF competition as the allocation/failure
       recovery path. */
  }

  const int sparse_full_diagonal_metis_row_economics =
    kls_sparse_full_diagonal_metis_row_input_economics(n, col_ptr);
  const int large_very_low_degree_full_diagonal =
    sparse_full_diagonal_metis_row_economics
      ? kls_sparse_full_diagonal_metis_row_input_profile(
          n, col_ptr, row_idx)
      : is_large_very_low_degree_full_diagonal_pattern(
          n, col_ptr, row_idx);
  if (!kls_analyze_defer_nd && large_very_low_degree_full_diagonal) {
    if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
      fprintf(stderr, "KLS choose: low-degree METIS branch n=%ld\n", (long)n);
    }
    kls_options metis_options = *options;
    metis_options.use_btf = 0;
    const int saved_sparse_full_diagonal_metis_row_ctx =
      kls_sparse_full_diagonal_metis_row_ctx;
    kls_sparse_full_diagonal_metis_row_ctx = 0;

    int status = analyze_with_ordering(n, col_ptr, row_idx, &metis_options,
                                       KLS_ORDERING_METIS, symbolic_out,
                                       common_out, separator_out);
    kls_sparse_full_diagonal_metis_row_ctx =
      saved_sparse_full_diagonal_metis_row_ctx;
    if (status == KLS_OK) {
      *selected_ordering_out = KLS_ORDERING_METIS;
      double selected_score = symbolic_score(*symbolic_out);
      *score_out = selected_score;
      return KLS_OK;
    }
  }

  kls_no_btf_retry_hopeless = 0;
#ifdef KLS_HAVE_MTMETIS
  kls_mt_nd_class =
    n >= 200000 &&
    is_large_sparse_diagonal_low_degree_pattern(n, col_ptr, row_idx);
#endif
  kls_options auto_options = *options;
  symbolic_options = options;
  const int large_spiked_metis_no_btf =
    is_large_nearly_diagonal_spiked_metis_pattern(n, col_ptr, row_idx);
  if (large_spiked_metis_no_btf && n <= 750000u && col_ptr[n] < 8u * n) {
  }
  if (should_start_auto_without_btf(n, col_ptr, row_idx, options,
                                    large_spiked_metis_no_btf)) {
    auto_options.use_btf = 0;
    symbolic_options = &auto_options;
  }

  if (kls_prestatic_ordering_ctx && options->threads >= 3 &&
      n > 150000u && n <= 750000u && (UF_long)col_ptr[n] <= 8000000 &&
      getenv("KLS_DISABLE_PS_CHOOSE_PAR") == NULL &&
      !is_large_sparse_diagonal_low_degree_pattern(n, col_ptr, row_idx)) {
    if (col_ptr[n] < 8 * n) {
      /* light-column matched pattern: the factor this ordering feeds
         is bounded by ~1e7-1e8 flops (rajat29's trial: 1.25e7 =
         0.11s serial), so the METIS/AMF candidates' own analyze cost
         (~1-2.4s here) can never repay itself.  AMD alone; the
         pre-static trial's quality gates police a bad draw (their
         rejection redo rebuilds the full competition). */
      trilinos_klu_l_symbolic *a_sym = NULL;
      trilinos_klu_l_common a_common;
      kls_separator_analysis a_sep;
      memset(&a_sep, 0, sizeof(a_sep));
      if (analyze_with_ordering(n, col_ptr, row_idx, symbolic_options,
                                KLS_ORDERING_AMD, &a_sym, &a_common,
                                &a_sep) == KLS_OK && a_sym != NULL) {
        if (kls_trace_pre_static_enabled()) {
          fprintf(stderr,
                  "KLS choose: light-column AMD-only score=%.4e\n",
                  symbolic_score(a_sym));
        }
        *symbolic_out = a_sym;
        *common_out = a_common;
        *selected_ordering_out = KLS_ORDERING_AMD;
        *score_out = symbolic_score(a_sym);
        kls_separator_analysis_move(separator_out, &a_sep);
        return KLS_OK;
      }
      kls_separator_analysis_clear(&a_sep);
      /* AMD failed: fall through to the full competition */
    }
    /* pre-static choose for the legacy-RNG class: METIS (the measured
       class-wide winner) analyzes on THIS thread so the only legacy
       libc-rand drawer stays serial, while AMD and AMF analyze on
       workers.  The comparisons replay in the canonical order at the
       join with the standard margins (AMF 0.95 vs the AMD base,
       METIS 0.90 vs the base winner). */
    struct kls_amf_spec_job amd_job, amf_job;
    pthread_t amd_tid, amf_tid;
    memset(&amd_job, 0, sizeof(amd_job));
    amd_job.n = n; amd_job.col_ptr = col_ptr; amd_job.row_idx = row_idx;
    amd_job.ordering = KLS_ORDERING_AMD;
    amd_job.options = *symbolic_options;
    amd_job.representation_verdict = NULL;
    amf_job = amd_job; amf_job.ordering = KLS_ORDERING_AMF;
    const int amd_sp =
      pthread_create(&amd_tid, NULL, kls_amf_spec_main, &amd_job) == 0;
    const int amf_sp =
      pthread_create(&amf_tid, NULL, kls_amf_spec_main, &amf_job) == 0;
    trilinos_klu_l_symbolic *m_sym = NULL;
    trilinos_klu_l_common m_common;
    kls_separator_analysis m_sep;
    memset(&m_sep, 0, sizeof(m_sep));
#ifdef KLS_HAVE_MTMETIS
    /* deterministic parallel NodeNDP for the class's choose-stage
       METIS: serial METIS separators with concurrent fresh-thread
       leaves, deterministic regardless of what ran before (the
       mt-metis mode-0 bridge measured call-order-dependent forests:
       its internal thread streams carry across calls, and the
       first-call forest scored 1.343e8 vs the class AMD's 1.337e8,
       losing the margin the same input won at 1.011e8 second-call). */
    kls_det_ndp_requested = getenv("KLS_DISABLE_PS_MTND") == NULL;
#endif
    const int m_status = analyze_with_ordering(n, col_ptr, row_idx,
                                               symbolic_options,
                                               KLS_ORDERING_METIS, &m_sym,
                                               &m_common, &m_sep);
#ifdef KLS_HAVE_MTMETIS
    kls_det_ndp_requested = 0;
#endif
    if (amd_sp) pthread_join(amd_tid, NULL);
    if (amf_sp) pthread_join(amf_tid, NULL);
    trilinos_klu_l_symbolic *b_sym = NULL;
    trilinos_klu_l_common b_common;
    kls_separator_analysis b_sep;
    kls_ordering b_ord = KLS_ORDERING_AMD;
    double b_score = DBL_MAX;
    int b_ok = 0;
    memset(&b_sep, 0, sizeof(b_sep));
    if (amd_sp && amd_job.status == KLS_OK && amd_job.symbolic != NULL) {
      b_sym = amd_job.symbolic; b_common = amd_job.common;
      b_score = amd_job.score; b_ok = 1;
      kls_separator_analysis_move(&b_sep, &amd_job.separator);
    }
    if (amf_sp && amf_job.status == KLS_OK && amf_job.symbolic != NULL) {
      if (b_ok && isfinite(amf_job.score) && b_score < DBL_MAX &&
          amf_job.score <= 0.95 * b_score) {
        trilinos_klu_l_free_symbolic(&b_sym, &b_common);
        kls_separator_analysis_clear(&b_sep);
        b_sym = amf_job.symbolic; b_common = amf_job.common;
        b_ord = KLS_ORDERING_AMF; b_score = amf_job.score;
        kls_separator_analysis_move(&b_sep, &amf_job.separator);
      } else if (!b_ok) {
        b_sym = amf_job.symbolic; b_common = amf_job.common;
        b_ord = KLS_ORDERING_AMF; b_score = amf_job.score; b_ok = 1;
        kls_separator_analysis_move(&b_sep, &amf_job.separator);
      } else {
        trilinos_klu_l_free_symbolic(&amf_job.symbolic, &amf_job.common);
        kls_separator_analysis_clear(&amf_job.separator);
      }
    }
    if (kls_trace_pre_static_enabled()) {
      fprintf(stderr,
              "KLS choose: ps-par m_status=%d m_sym=%d m_score=%.4e "
              "b_ok=%d b_ord=%d b_score=%.4e\n",
              m_status, m_sym != NULL,
              m_sym != NULL ? symbolic_score(m_sym) : -1.0,
              b_ok, (int)b_ord, b_score);
    }
    if (m_status == KLS_OK && m_sym != NULL) {
      const double m_score = symbolic_score(m_sym);
      if (!b_ok ||
          (isfinite(m_score) && m_score > 0.0 && b_score < DBL_MAX &&
           m_score <= 0.90 * b_score)) {
        if (b_ok) {
          trilinos_klu_l_free_symbolic(&b_sym, &b_common);
          kls_separator_analysis_clear(&b_sep);
        }
        b_sym = m_sym; b_common = m_common;
        b_ord = KLS_ORDERING_METIS; b_score = m_score; b_ok = 1;
        kls_separator_analysis_move(&b_sep, &m_sep);
      } else {
        trilinos_klu_l_free_symbolic(&m_sym, &m_common);
        kls_separator_analysis_clear(&m_sep);
      }
    } else {
      kls_separator_analysis_clear(&m_sep);
    }
    if (b_ok) {
      *symbolic_out = b_sym;
      *common_out = b_common;
      *selected_ordering_out = b_ord;
      *score_out = b_score;
      kls_separator_analysis_move(separator_out, &b_sep);
      return KLS_OK;
    }
    /* every candidate failed: fall through to the standard path */
  }

  if (n <= 49152u && !kls_analyze_defer_nd &&
      kls_analyze_nd_race_solver != NULL &&
      kls_analyze_nd_race_solver->metis_race != NULL &&
      should_start_auto_with_metis(n, col_ptr, row_idx,
                                   large_spiked_metis_no_btf)) {
    /* The background race covers the NodeND trial for this small
       METIS-start matrix: deficiency ordering carries the first numeric
       (measured refactor parity on this class) and the factor-time
       timed promotion arbitrates the ND comparison. */
    int amf_status = analyze_with_ordering(
      n, col_ptr, row_idx, symbolic_options, KLS_ORDERING_AMF,
      symbolic_out, common_out, separator_out);
    if (amf_status == KLS_OK) {
      const double amf_score = symbolic_score(*symbolic_out);
      if (isfinite(amf_score) && amf_score > 0.0 &&
          amf_score < DBL_MAX / 4.0) {
        *selected_ordering_out = KLS_ORDERING_AMF;
        *score_out = amf_score;
        return KLS_OK;
      }
      trilinos_klu_l_free_symbolic(symbolic_out, common_out);
      kls_separator_analysis_clear(separator_out);
    }
  }
  if (!kls_analyze_defer_nd &&
      !(kls_analyze_nd_race_solver != NULL &&
        kls_analyze_nd_race_solver->metis_race != NULL) &&
      should_start_auto_with_metis(n, col_ptr, row_idx,
                                   large_spiked_metis_no_btf)) {
    if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
      fprintf(stderr, "KLS choose: METIS-start branch n=%ld btf=%d\n",
              (long)n, symbolic_options->use_btf);
    }
    if (n <= 30000u &&
        is_small_spiked_low_diagonal_pattern(n, col_ptr, row_idx)) {
      /* Fused dense blocks behind spikes: deficiency ordering dominates
         nested dissection here (TSOPF: est fill 1.31e5 vs 1.50e5,
         tight-tol factor 6.4e4 vs 8.0e4) and skipping NodeND saves most
         of the analysis (13.5 -> ~7ms).  Fall through to METIS if the
         AMF analysis fails. */
      int amf_status = analyze_with_ordering(
        n, col_ptr, row_idx, symbolic_options, KLS_ORDERING_AMF,
        symbolic_out, common_out, separator_out);
      if (amf_status == KLS_OK) {
        const double amf_score = symbolic_score(*symbolic_out);
        if (isfinite(amf_score) && amf_score > 0.0 &&
            amf_score < DBL_MAX / 4.0) {
          *selected_ordering_out = KLS_ORDERING_AMF;
          *score_out = amf_score;
          return KLS_OK;
        }
        trilinos_klu_l_free_symbolic(symbolic_out, common_out);
        kls_separator_analysis_clear(separator_out);
      }
    }
    int status = analyze_with_ordering(n, col_ptr, row_idx, symbolic_options,
                                       KLS_ORDERING_METIS, symbolic_out,
                                       common_out, separator_out);
    if (status == KLS_OK) {
      *selected_ordering_out = KLS_ORDERING_METIS;
      double selected_score = symbolic_score(*symbolic_out);
      if (!metis_start_should_skip_no_btf_retry(n, col_ptr, row_idx)) {
        maybe_retry_without_btf(n, col_ptr, row_idx, symbolic_options,
                                KLS_ORDERING_METIS, symbolic_out, common_out,
                                &selected_score, 0, separator_out, NULL);
      }
      if (n <= 30000u) {
        /* Deficiency ordering often beats nested dissection on the
           small spiked/fused-block circuit classes routed here (TSOPF:
           est fill 1.31e5 vs METIS 1.50e5, post-tight-tol factor 6.4e4
           vs 8.0e4); the extra analyze costs ~1ms at this size and the
           promotion keeps the strict-estimate margin. */
        maybe_promote_symbolic_ordering(
          n, col_ptr, row_idx, symbolic_options, KLS_ORDERING_AMF,
          symbolic_out, common_out, selected_ordering_out,
          &selected_score, 0.90, separator_out);
      }
#ifdef KLS_HAVE_SCOTCH
      if (should_try_symbolic_nested_dissection_before_numeric(
            n, *symbolic_out, *selected_ordering_out, KLS_ORDERING_SCOTCH,
            selected_score)) {
        maybe_promote_symbolic_ordering(
          n, col_ptr, row_idx, symbolic_options, KLS_ORDERING_SCOTCH,
          symbolic_out, common_out, selected_ordering_out,
          &selected_score, 0.90, separator_out);
      }
#endif
      *score_out = selected_score;
      return KLS_OK;
    }
  }

generic_ordering_tournament:;
  struct kls_amf_spec_job amf_spec[3];
  pthread_t amf_spec_tid[3];
  int amf_spec_spawned[3] = {0, 0, 0};
  int amf_spec_configured = 0;
  const kls_ordering amf_spec_ordering[3] = {
    KLS_ORDERING_AMF, KLS_ORDERING_AMMF, KLS_ORDERING_AMF3
  };
  const int amf_trial_disabled =
    getenv("KLS_DISABLE_AMF_TRIAL") != NULL;
  /* Sum of squared column degrees is a cheap upper-level proxy for the pair
     visits that a different elimination ordering can save.  Charge it over
     the caller's complete numeric horizon; below 5e8 visits, another full
     quotient-graph analysis cannot repay its setup on this machine class.
     The former low floors still launched two graph workers on cache-resident
     matched factors whose whole retained numeric horizon was cheaper than
     those workers (the measured setup exceeded five steady updates).
     This is a workload budget, independent of dimensions and family. */
  const double column_pair_work = kls_column_pair_work(n, col_ptr);
  const double requested_numeric_horizon = 1.0 +
    (double)options->expected_refactorizations +
    (double)options->expected_solves;
  const int amf_trial_economically_relevant =
    !kls_repeated_update_workload(options) ||
    !isfinite(column_pair_work) ||
    column_pair_work * requested_numeric_horizon >= 5.0e8;
  /* Repeated numeric work can amortize a broader minimum-fill portfolio.
     Bound the extra analyses by caller-supplied lifecycle and transient
     pattern-storage budgets, not by matrix identity or topology.  With at
     least four requested workers the independent quotient-graph passes
     overlap the AMD baseline.  A two-worker caller can also afford the two
     stronger formulations when it declares a long H100-style lifecycle and
     the pair-work budget is high enough to repay their serial setup. */
  const double ordering_pattern_bytes =
    ((double)n + 1.0 + (double)col_ptr[n]) * (double)sizeof(UF_long);
  const int two_worker_amf_portfolio =
    kls_repeated_update_workload(options) && options->threads >= 2 &&
    options->expected_refactorizations >= 64 &&
    isfinite(column_pair_work) &&
    column_pair_work * requested_numeric_horizon >= 1.0e9;
  const int broad_amf_portfolio =
    kls_repeated_update_workload(options) &&
    amf_trial_economically_relevant &&
    (options->threads >= 4 || two_worker_amf_portfolio) &&
    /* Bound concurrent quotient-graph speculation by the immutable pattern
       footprint rather than by row/nonzero coordinates.  A 128 MiB input
       representation leaves room for the three independent AMMF/AMF3 graph
       states under the same worker-memory envelope used by the numeric
       portfolio, while admitting large sparse problems whose low-fill
       ordering choice dominates hundreds of requested numeric updates. */
    isfinite(ordering_pattern_bytes) &&
    ordering_pattern_bytes <= 128.0 * 1024.0 * 1024.0;
  const int stage_cache_ordering_portfolio =
    getenv("KLS_DISABLE_CACHE_ORDERING_PORTFOLIO_STAGING") == NULL &&
    isfinite(ordering_pattern_bytes) &&
    ordering_pattern_bytes >= 256.0 * 1024.0 &&
    ordering_pattern_bytes <= 2.0 * 1024.0 * 1024.0 &&
    isfinite(column_pair_work) &&
    column_pair_work >=
      2.5 * (double)col_ptr[n] * requested_numeric_horizon;
  const double no_btf_lifecycle_work =
    (double)col_ptr[n] * requested_numeric_horizon;
  const int parallel_no_btf_trial_wanted =
    symbolic_options->use_btf && n >= 4000u &&
    kls_repeated_update_workload(options) && options->threads >= 4 &&
    (!broad_amf_portfolio || options->threads >= 8) &&
    n <= 1000000u && col_ptr[n] <= 8000000u &&
    isfinite(no_btf_lifecycle_work) && no_btf_lifecycle_work >= 1.0e7 &&
    isfinite(column_pair_work) &&
    column_pair_work * requested_numeric_horizon >= 2.0e8 &&
    getenv("KLS_DISABLE_GENERIC_BTF_VALUE_SELECTION") == NULL &&
    getenv("KLS_DISABLE_PARALLEL_NO_BTF_TRIAL") == NULL;
  /* Eight available workers can also overlap the literal AMF3 clique bound
     with the historically useful relaxed bound.  Published results show
     that a tighter local bound is not universally a better greedy order, so
     retain this as a measured symbolic portfolio rather than replacing the
     incumbent.  An explicit power request suppresses the comparison. */
  const int amf3_span_portfolio =
    broad_amf_portfolio && options->threads >= 8 &&
    !kls_prestatic_ordering_ctx &&
    !kls_matched_low_pair_ordering_ctx &&
    getenv("KLS_AMF3_POWER") == NULL;
  const int amf_spec_count = amf3_span_portfolio ? 3 :
    (broad_amf_portfolio ? 2 : 1);
  atomic_int generic_representation_verdict;
  atomic_init(&generic_representation_verdict,
              !symbolic_options->use_btf ? 2 : 0);
  struct kls_amf_spec_job aligned_btf_ammf_spec;
  pthread_t aligned_btf_ammf_spec_tid;
  int aligned_btf_ammf_spec_spawned = 0;
  atomic_int aligned_btf_representation_verdict;
  atomic_init(&aligned_btf_representation_verdict, 1);
  kls_no_btf_spec_job no_btf_spec;
  pthread_t no_btf_spec_tid;
  int no_btf_spec_spawned = 0;
  memset(&no_btf_spec, 0, sizeof(no_btf_spec));
  if (amf_trial_economically_relevant &&
      (n >= 30000 || broad_amf_portfolio) && options->threads >= 2 &&
      !amf_trial_disabled &&
      getenv("KLS_DISABLE_AMF_SPEC") == NULL) {
    memset(amf_spec, 0, sizeof(amf_spec));
    amf_spec_configured = 1;
    for (int i = 0; i < amf_spec_count; ++i) {
      /* In the recurring portfolio AMMF and AMF3 are the strengthened
         minimum-fill formulations; retaining the original AMF alongside
         both only adds a bandwidth-contending graph pass and was dominated
         in fill/work by one of them throughout the cross-family audit.  The
         smaller/one-shot path keeps the historical AMF challenger. */
      const int ordering_index = broad_amf_portfolio
        ? (i == 0 ? 1 : 2) : i;
      amf_spec[i].n = n;
      amf_spec[i].col_ptr = col_ptr;
      amf_spec[i].row_idx = row_idx;
      amf_spec[i].ordering = amf_spec_ordering[ordering_index];
      amf_spec[i].options = *symbolic_options;
      /* The AMD representation race is already evaluating no-BTF under the
         same lifecycle budget.  Start independent AMF workers on that arm
         too; if BTF wins, their existing verdict repair reruns them there.
         This removes a duplicate BTF graph pass without narrowing the
         ordering portfolio. */
      if (parallel_no_btf_trial_wanted) {
        amf_spec[i].options.use_btf = 0;
      }
      /* The recurrent executor pays for dependency depth, whereas quotient
         fill approximates only storage.  Measure every member of the broad
         portfolio so a later candidate can be compared with the current
         minimum-fill winner on the same representation. */
      amf_spec[i].measure_etree_levels = amf3_span_portfolio;
      /* The third arm varies the clique bound, not the scoring power.  It
         therefore must not inherit the lower-power variant's numeric census
         and restoration policy. */
      amf_spec[i].amf3_span_variant = 0;
      amf_spec[i].has_amf3_power_override = 0;
      amf_spec[i].amf3_power_override =
        0.445;
      amf_spec[i].amf3_exact_clique_bound =
        amf3_span_portfolio && i == 2;
      amf_spec[i].representation_verdict = &generic_representation_verdict;
      if (!stage_cache_ordering_portfolio &&
          pthread_create(&amf_spec_tid[i], NULL, kls_amf_spec_main,
                         &amf_spec[i]) == 0) {
        amf_spec_spawned[i] = 1;
      }
    }
  }
  /* BTF and no-BTF analyses of the same ordering have no shared mutable
     state.  When the requested recurring graph work can repay a speculative
     representation comparison, overlap the AMD challenger with the
     incumbent instead of serializing it after AMD/BTF.  Cheap portfolios
     wait for the incumbent's actual block structure and still run the same
     challenger on demand when that result makes it relevant.  Keep enough
     caller-supplied worker budget for the broader AMF portfolio when that is
     also active. */
  if (parallel_no_btf_trial_wanted &&
      !stage_cache_ordering_portfolio) {
    no_btf_spec.n = n;
    no_btf_spec.col_ptr = col_ptr;
    no_btf_spec.row_idx = row_idx;
    no_btf_spec.ordering = KLS_ORDERING_AMD;
    no_btf_spec.options = *symbolic_options;
    no_btf_spec.status = KLS_ERR_ANALYZE_FAILED;
    no_btf_spec_spawned =
      pthread_create(&no_btf_spec_tid, NULL, kls_no_btf_spec_main,
                     &no_btf_spec) == 0;
  }
  /* KLU does not publish fill/work estimates for COLAMD, so its symbolic
     score is DBL_MAX.  Such a candidate cannot displace a successful,
     finitely-scored AMD result under the comparison below.  Retain an A/B
     hook and an AMD-failure fallback, but do not spend a concurrent full
     graph traversal on a candidate that the selector cannot accept. */
  struct kls_amf_spec_job colamd_spec;
  pthread_t colamd_spec_tid;
  int colamd_spec_spawned = 0;
  if (n >= 8192u && options->threads >= 4 &&
      getenv("KLS_ENABLE_UNSCORABLE_COLAMD_SPEC") != NULL &&
      getenv("KLS_DISABLE_COLAMD_SPEC") == NULL) {
    memset(&colamd_spec, 0, sizeof(colamd_spec));
    colamd_spec.n = n;
    colamd_spec.col_ptr = col_ptr;
    colamd_spec.row_idx = row_idx;
    colamd_spec.ordering = KLS_ORDERING_COLAMD;
    colamd_spec.options = *symbolic_options;
    if (pthread_create(&colamd_spec_tid, NULL, kls_amf_spec_main,
                       &colamd_spec) == 0) {
      colamd_spec_spawned = 1;
    }
  }
  if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
    fprintf(stderr, "KLS choose: base loop n=%ld\n", (long)n);
  }
  const kls_ordering candidates[] = {KLS_ORDERING_AMD, KLS_ORDERING_COLAMD};
  const size_t candidate_count = 1u;
  trilinos_klu_l_symbolic *best_symbolic = NULL;
  trilinos_klu_l_common best_common;
  kls_ordering best_ordering = KLS_ORDERING_AMD;
  double best_score = 0.0;
  kls_separator_analysis best_separator;
  memset(&best_separator, 0, sizeof(best_separator));
  int any_ok = 0;

  for (size_t i = 0; i < candidate_count; ++i) {
    trilinos_klu_l_symbolic *candidate_symbolic = NULL;
    trilinos_klu_l_common candidate_common;
    kls_separator_analysis candidate_separator;
    memset(&candidate_separator, 0, sizeof(candidate_separator));
    int status = analyze_with_ordering(n, col_ptr, row_idx, symbolic_options, candidates[i],
                                       &candidate_symbolic, &candidate_common,
                                       &candidate_separator);
    if (status != KLS_OK) {
      continue;
    }
    const double score = symbolic_score(candidate_symbolic);
    double selected_score = score;

    if (!any_ok || selected_score < best_score) {
      if (best_symbolic != NULL) {
        trilinos_klu_l_free_symbolic(&best_symbolic, &best_common);
      }
      kls_separator_analysis_clear(&best_separator);
      best_symbolic = candidate_symbolic;
      best_common = candidate_common;
      best_ordering = candidates[i];
      best_score = selected_score;
      kls_separator_analysis_move(&best_separator, &candidate_separator);
      any_ok = 1;
    } else {
      trilinos_klu_l_free_symbolic(&candidate_symbolic, &candidate_common);
      kls_separator_analysis_clear(&candidate_separator);
    }
    /* Stop once an ordering provides a real fill estimate; unknown estimates
       should not displace a known-good circuit ordering. */

  }

  if (!any_ok && !colamd_spec_spawned) {
    /* AMD allocation/analysis failure must not turn the optimization above
       into a robustness regression.  COLAMD remains a valid symbolic even
       though it is unscored, so use it as a fallback rather than a peer. */
    trilinos_klu_l_symbolic *fallback_symbolic = NULL;
    trilinos_klu_l_common fallback_common;
    kls_separator_analysis fallback_separator;
    memset(&fallback_separator, 0, sizeof(fallback_separator));
    if (analyze_with_ordering(
          n, col_ptr, row_idx, symbolic_options, KLS_ORDERING_COLAMD,
          &fallback_symbolic, &fallback_common, &fallback_separator) ==
          KLS_OK && fallback_symbolic != NULL) {
      best_symbolic = fallback_symbolic;
      best_common = fallback_common;
      best_ordering = KLS_ORDERING_COLAMD;
      best_score = symbolic_score(fallback_symbolic);
      kls_separator_analysis_move(&best_separator, &fallback_separator);
      any_ok = 1;
    } else {
      kls_separator_analysis_clear(&fallback_separator);
    }
  }

  if (colamd_spec_spawned) {
    pthread_join(colamd_spec_tid, NULL);
    if (colamd_spec.status == KLS_OK && colamd_spec.symbolic != NULL) {
      if (!any_ok || colamd_spec.score < best_score) {
        if (best_symbolic != NULL) {
          trilinos_klu_l_free_symbolic(&best_symbolic, &best_common);
        }
        kls_separator_analysis_clear(&best_separator);
        best_symbolic = colamd_spec.symbolic;
        best_common = colamd_spec.common;
        best_ordering = KLS_ORDERING_COLAMD;
        best_score = colamd_spec.score;
        kls_separator_analysis_move(&best_separator,
                                    &colamd_spec.separator);
        any_ok = 1;
      } else {
        trilinos_klu_l_free_symbolic(&colamd_spec.symbolic,
                                     &colamd_spec.common);
        kls_separator_analysis_clear(&colamd_spec.separator);
      }
    } else {
      kls_separator_analysis_clear(&colamd_spec.separator);
    }
  }

  if (!any_ok) {
    {
      /* Release any ordering workers waiting for the failed base verdict
         before joining them on the cleanup path. */
      atomic_store_explicit(&generic_representation_verdict,
                            symbolic_options->use_btf ? 1 : 2,
                            memory_order_release);
    }
    if (no_btf_spec_spawned) {
      pthread_join(no_btf_spec_tid, NULL);
      kls_no_btf_spec_clear(&no_btf_spec);
    }
    for (int i = 0; i < amf_spec_count; ++i) {
      if (amf_spec_spawned[i]) {
        pthread_join(amf_spec_tid[i], NULL);
        if (amf_spec[i].status == KLS_OK) {
          trilinos_klu_l_free_symbolic(&amf_spec[i].symbolic,
                                       &amf_spec[i].common);
        }
        kls_separator_analysis_clear(&amf_spec[i].separator);
      }
    }
    if (aligned_btf_ammf_spec_spawned) {
      pthread_join(aligned_btf_ammf_spec_tid, NULL);
      if (aligned_btf_ammf_spec.symbolic != NULL) {
        trilinos_klu_l_free_symbolic(
          &aligned_btf_ammf_spec.symbolic,
          &aligned_btf_ammf_spec.common);
      }
      kls_separator_analysis_clear(
        &aligned_btf_ammf_spec.separator);
    }
    return KLS_ERR_ANALYZE_FAILED;
  }

  int suppress_fallback_amf = 0;
  {
    /* Settle the representation axis before comparing ordering candidates.
       Speculative minimum-fill workers perform the same bounded comparison,
       so their scores are meaningful only when both sides describe the same
       BTF/no-BTF representation.  Running this while those workers are live
       also overlaps the extra symbolic pass with their independent graphs. */
    if (no_btf_spec_spawned) {
      pthread_join(no_btf_spec_tid, NULL);
    }
    const double settled_base_fill = symbolic_score(best_symbolic);
    const double settled_base_flops = best_symbolic->est_flops;
    const double settled_base_lifecycle =
      (settled_base_fill + settled_base_flops) *
      requested_numeric_horizon;
    /* A retained symbolic whose estimated numeric work is already within
       twenty percent of one pass over the input nonzeros is at the sparse
       arithmetic floor: another quotient-graph ordering cannot remove even
       one input-stream pass per update.  For cache-resident portfolios, the
       challenger construction itself then exceeds the most it could save.
       Settle the existing representation and ordering from this realized
       work certificate.  Higher-work factors keep the complete portfolio;
       the control restores it for same-binary audits. */
    const int near_sparse_work_floor =
      stage_cache_ordering_portfolio &&
      getenv("KLS_DISABLE_NEAR_SPARSE_ORDERING_PORTFOLIO_SETTLE") == NULL &&
      isfinite(settled_base_flops) && settled_base_flops > 0.0 &&
      settled_base_flops <= 1.20 * (double)col_ptr[n];
    suppress_fallback_amf = near_sparse_work_floor;
    const int cache_challenger_economic =
      !near_sparse_work_floor &&
      (!stage_cache_ordering_portfolio ||
      !isfinite(settled_base_lifecycle) ||
      settled_base_lifecycle <= 0.0 ||
      column_pair_work < 0.90 * settled_base_lifecycle);
    if (cache_challenger_economic) {
      maybe_retry_without_btf(
        n, col_ptr, row_idx, symbolic_options, best_ordering,
        &best_symbolic, &best_common, &best_score, 1, &best_separator,
        no_btf_spec_spawned ? &no_btf_spec : NULL);
    }
    if (no_btf_spec_spawned) {
      kls_no_btf_spec_clear(&no_btf_spec);
    }
    atomic_store_explicit(&generic_representation_verdict,
                          best_symbolic->do_btf ? 1 : 2,
                          memory_order_release);
    if (stage_cache_ordering_portfolio && cache_challenger_economic &&
        amf_spec_configured) {
      for (int i = 0; i < amf_spec_count; ++i) {
        amf_spec[i].options.use_btf = best_symbolic->do_btf ? 1 : 0;
        if (pthread_create(&amf_spec_tid[i], NULL, kls_amf_spec_main,
                           &amf_spec[i]) == 0) {
          amf_spec_spawned[i] = 1;
        }
      }
    }
    if (parallel_no_btf_trial_wanted && broad_amf_portfolio &&
        !best_symbolic->do_btf && options->threads >= 8 &&
        kls_generic_btf_capture_candidate != NULL &&
        kls_generic_btf_capture_candidate->generic_btf_value_symbolic !=
          NULL &&
        getenv("KLS_DISABLE_PARALLEL_ALIGNED_BTF_AMMF") == NULL) {
      /* The measured representation comparison retained a BTF value arm
         while selecting no-BTF for the symbolic incumbent.  Start its
         aligned AMMF candidate now, beside the still-running ordering
         workers.  Delaying until this concrete signal avoids charging an
         unused BTF graph pass to portfolios that simply retain BTF. */
      memset(&aligned_btf_ammf_spec, 0,
             sizeof(aligned_btf_ammf_spec));
      aligned_btf_ammf_spec.n = n;
      aligned_btf_ammf_spec.col_ptr = col_ptr;
      aligned_btf_ammf_spec.row_idx = row_idx;
      aligned_btf_ammf_spec.ordering = KLS_ORDERING_AMMF;
      aligned_btf_ammf_spec.options = *symbolic_options;
      aligned_btf_ammf_spec.options.use_btf = 1;
      aligned_btf_ammf_spec.representation_verdict =
        &aligned_btf_representation_verdict;
      aligned_btf_ammf_spec_spawned =
        pthread_create(&aligned_btf_ammf_spec_tid, NULL,
                       kls_amf_spec_main,
                       &aligned_btf_ammf_spec) == 0;
    }
    /* The base AMD comparison has now settled the representation from
       realized symbolic fill/work.  Launch NodeND in that same coordinate
       frame while the independent minimum-fill workers finish.  This keeps
       the useful analysis overlap without guessing BTF from the caller's
       request; if the tournament takes an early route or the resource gate
       declines, the ordinary synchronous portfolio remains the fallback. */
    if (kls_generic_btf_capture_candidate != NULL &&
        kls_generic_btf_capture_candidate->generic_nd_overlap_target !=
          NULL) {
      const double settled_n = (double)n;
      const double settled_active_n =
        best_symbolic->do_btf && best_symbolic->maxblock > 0u
          ? (double)best_symbolic->maxblock : settled_n;
      const double settled_fill = symbolic_score(best_symbolic);
      const double settled_flops = best_symbolic->est_flops;
      const double settled_lifecycle_work = settled_flops *
        (1.0 + (double)options->expected_refactorizations);
      const int relaxed_recurring_economic =
        getenv("KLS_DISABLE_RELAXED_GENERIC_ND_WORK_GATE") == NULL &&
        options->expected_refactorizations >= 16 &&
        settled_flops >= 2048.0 * settled_n &&
        isfinite(settled_lifecycle_work) &&
        settled_lifecycle_work >= 5.0e10;
      const int settled_nd_economic =
        isfinite(settled_fill) &&
        settled_fill >= 16.0 * settled_active_n &&
        (settled_flops >= 4096.0 * settled_active_n ||
         relaxed_recurring_economic);
      /* Apply the promotion portfolio's realized-work floor before spawning
         its speculative worker.  Previously low-work incumbents could never
         admit NodeND but still paid a complete graph analysis.  This gate is
         based only on the settled symbolic and requested numeric horizon. */
      if (settled_nd_economic) {
        (void)kls_generic_nd_trial_start(
          kls_generic_btf_capture_candidate->generic_nd_overlap_target,
          kls_generic_btf_capture_candidate, options,
          best_symbolic->do_btf ? 1 : 0);
      }
    }
  }

  /* AMF (approximate minimum fill) rides the same quotient graph as AMD
     at similar analysis cost; try it whenever AMD won the base round and
     keep it only on a clear estimate win (task #11: must only win where
     estimates say so).  The speculative trial ran concurrently with the
     base analyze; the promotion comparison is unchanged. */
  {
    int amf_done = 0;
    UF_long best_etree_levels = 0u;
    for (int spec_index = 0; spec_index < amf_spec_count; ++spec_index) {
      if (!amf_spec_spawned[spec_index]) {
        continue;
      }
      struct kls_amf_spec_job *spec = &amf_spec[spec_index];
      pthread_join(amf_spec_tid[spec_index], NULL);
      const int promotable =
        (1) &&
        best_score < DBL_MAX &&
        isfinite(best_score) && best_score > 0.0;
      const int spec_valid =
        spec->status == KLS_OK && spec->symbolic != NULL &&
        promotable &&
        (spec->symbolic->do_btf ? 1 : 0) ==
          (best_symbolic->do_btf ? 1 : 0);
      if (spec_valid) {
        /* Symbolic fill/work estimates are deliberately approximate.  Keep a
           clear fill margin by default; a repeated lifecycle may instead
           spend at most one percent of estimated storage for a decisive work
           reduction, with the narrower AMF3 refinement bounded separately.
           These confidence margins are independent of dimensions or family. */
        const double amf_ratio = 0.90;
        if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
          fprintf(stderr,
                  "KLS choose: %s score %.4e/flops %.4e vs "
                  "%.4e/%.4e (need fill %.2fx)\n",
                  kls_ordering_name(spec->ordering), spec->score,
                  spec->symbolic->est_flops, best_score,
                  best_symbolic->est_flops, amf_ratio);
        }
        const int decisive_work_pareto =
          isfinite(spec->score) && spec->score <= best_score &&
          spec->symbolic->est_flops > 0.0 &&
          best_symbolic->est_flops > 0.0 &&
          spec->symbolic->est_flops <= 0.80 * best_symbolic->est_flops;
        const int repeated_work_for_bounded_fill = kls_repeated_update_workload(options) &&
          isfinite(spec->score) && spec->score <= 1.01 * best_score &&
          spec->symbolic->est_flops > 0.0 &&
          best_symbolic->est_flops > 0.0 &&
          spec->symbolic->est_flops <= 0.80 * best_symbolic->est_flops;
        const int minimum_fill_pareto_refinement = (best_ordering == KLS_ORDERING_AMMF ||
           best_ordering == KLS_ORDERING_AMF3) &&
          isfinite(spec->score) && spec->score <= best_score &&
          spec->symbolic->est_flops > 0.0 &&
          best_symbolic->est_flops > 0.0 &&
          spec->symbolic->est_flops <= 0.90 * best_symbolic->est_flops;
        /* AMF3's tighter clique approximation can also earn a narrower
           refinement after AMMF when estimated fill improves and work stays
           effectively Pareto-neutral.  Require either predicted work or
           factor-entry traffic saved over the caller-declared lifecycle to
           repay a complete ordering/factor decision.  This rejects cheap
           near-ties whose approximate scores cannot justify AMF3's greater
           exposure to numerical pivot fill, without classifying the graph
           by name, dimensions, or topology. */
        const double amf3_fill_saving = best_score - spec->score;
        const double amf3_work_saving =
          best_symbolic->est_flops - spec->symbolic->est_flops;
        const int amortized_amf3_work_pareto_refinement = kls_repeated_update_workload(options) &&
          best_ordering == KLS_ORDERING_AMMF &&
          spec->ordering == KLS_ORDERING_AMF3 &&
          getenv("KLS_DISABLE_AMF3_WORK_PARETO") == NULL &&
          isfinite(spec->score) && spec->score <= best_score &&
          spec->symbolic->est_flops > 0.0 &&
          best_symbolic->est_flops > 0.0 &&
          spec->symbolic->est_flops <= 0.99 * best_symbolic->est_flops &&
          isfinite(amf3_work_saving) && amf3_work_saving > 0.0 &&
          amf3_work_saving * requested_numeric_horizon >= 5.0e8;
        const int amortized_amf3_pareto_refinement = kls_repeated_update_workload(options) &&
          best_ordering == KLS_ORDERING_AMMF &&
          spec->ordering == KLS_ORDERING_AMF3 &&
          isfinite(spec->score) && spec->score <= 0.996 * best_score &&
          spec->symbolic->est_flops > 0.0 &&
          best_symbolic->est_flops > 0.0 &&
          spec->symbolic->est_flops <= 1.005 * best_symbolic->est_flops &&
          ((isfinite(amf3_work_saving) && amf3_work_saving > 0.0 &&
            amf3_work_saving * requested_numeric_horizon >= 5.0e8) ||
           (isfinite(amf3_fill_saving) && amf3_fill_saving > 0.0 &&
            amf3_fill_saving * requested_numeric_horizon >= 2.0e6));
        /* A broad portfolio has already paid to construct both quotient
           graphs.  A shallower exact elimination tree exposes more useful
           parallel dependency levels even when approximate fill is a small
           amount worse.  Admit that alternative only inside a tight
           storage/work Pareto band and only for a decisive depth reduction;
           the ordinary value-time factor and all numeric guards remain the
           final authority. */
        const int dependency_depth_pareto_refinement = kls_repeated_update_workload(options) &&
          best_ordering == KLS_ORDERING_AMMF &&
          spec->ordering == KLS_ORDERING_AMF3 &&
          best_etree_levels > 0u && spec->etree_levels > 0u &&
          isfinite(spec->score) && spec->score <= 1.01 * best_score &&
          spec->symbolic->est_flops > 0.0 &&
          best_symbolic->est_flops > 0.0 &&
          spec->symbolic->est_flops <= 1.02 * best_symbolic->est_flops &&
          10.0 * (double)spec->etree_levels <=
            9.0 * (double)best_etree_levels;
        /* Approximate fill does not describe dependency stalls in the
           retained parallel refactor.  Compare a lower-power AMF3 ordering
           with the ordinary AMF3 run on the identical representation.  A
           close general candidate stays within the narrow Pareto band.  In
           an exact value-matched frame, also admit a decisive span reduction
           under bounded five-percent-versus-control and twenty-percent-
           versus-winner work caps: the following real KLU factor still has
           to beat the incumbent on fill, work, pivots, and conditioning. */
        const struct kls_amf_spec_job *ordinary_amf3 =
          amf3_span_portfolio ? &amf_spec[1] : NULL;
        const int span_parallel_refinement =
          spec->amf3_span_variant &&
          kls_repeated_update_workload(options) &&
          ordinary_amf3 != NULL &&
          ordinary_amf3->status == KLS_OK &&
          ordinary_amf3->score > 0.0 &&
          isfinite(ordinary_amf3->score) &&
          ordinary_amf3->est_flops > 0.0 &&
          isfinite(ordinary_amf3->est_flops) &&
          ordinary_amf3->etree_levels > 0u &&
          spec->etree_levels > 0u &&
          spec->score <= 1.01 * ordinary_amf3->score &&
          spec->symbolic->est_flops <= 1.02 * ordinary_amf3->est_flops &&
          spec->score <= 1.01 * best_score &&
          spec->symbolic->est_flops <= 1.02 * best_symbolic->est_flops &&
          /* A shallow five-percent depth win is inside the error of the
             approximate fill/work model: one audited variant reduced the
             etree by seven percent yet made the realized numeric cycle 2.5x
             slower.  Require a decisive ten-percent dependency-depth win;
             the independently useful twelve-percent variant remains in the
             portfolio. */
          10.0 * (double)spec->etree_levels <=
            9.0 * (double)ordinary_amf3->etree_levels;
        const int matched_decisive_span_refinement =
          kls_value_matched_ordering_ctx &&
          spec->amf3_span_variant &&
          kls_repeated_update_workload(options) &&
          getenv("KLS_DISABLE_MATCHED_DECISIVE_SPAN_REFINEMENT") == NULL &&
          ordinary_amf3 != NULL &&
          ordinary_amf3->status == KLS_OK &&
          ordinary_amf3->score > 0.0 &&
          isfinite(ordinary_amf3->score) &&
          ordinary_amf3->est_flops > 0.0 &&
          isfinite(ordinary_amf3->est_flops) &&
          ordinary_amf3->etree_levels > 0u &&
          spec->etree_levels > 0u &&
          spec->score <= 1.02 * ordinary_amf3->score &&
          spec->symbolic->est_flops <= 1.05 * ordinary_amf3->est_flops &&
          spec->score <= 1.05 * best_score &&
          spec->symbolic->est_flops <= 1.20 * best_symbolic->est_flops &&
          5.0 * (double)spec->etree_levels <=
            3.0 * (double)ordinary_amf3->etree_levels;
        if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL &&
            spec->amf3_span_variant) {
          fprintf(stderr,
                  "KLS choose: AMF3 span refinement fill %.4e/%.4e "
                  "flops %.4e/%.4e levels %ld/%ld -> %s\n",
                  spec->score, ordinary_amf3 != NULL
                    ? ordinary_amf3->score : DBL_MAX,
                  spec->symbolic->est_flops, ordinary_amf3 != NULL
                    ? ordinary_amf3->est_flops : DBL_MAX,
                  (long)spec->etree_levels,
                  ordinary_amf3 != NULL
                    ? (long)ordinary_amf3->etree_levels : 0L,
                  (span_parallel_refinement ||
                   matched_decisive_span_refinement) ? "adopt" : "reject");
        }
        if ((isfinite(spec->score) &&
            spec->score <= amf_ratio * best_score) ||
            decisive_work_pareto ||
            repeated_work_for_bounded_fill ||
            minimum_fill_pareto_refinement ||
            amortized_amf3_work_pareto_refinement ||
            amortized_amf3_pareto_refinement ||
            dependency_depth_pareto_refinement ||
            span_parallel_refinement ||
            matched_decisive_span_refinement) {
          trilinos_klu_l_free_symbolic(&best_symbolic, &best_common);
          best_symbolic = spec->symbolic;
          best_common = spec->common;
          best_ordering = spec->ordering;
          best_score = spec->score;
          best_etree_levels = spec->etree_levels;
          if (spec->amf3_span_variant) {
            kls_amf3_span_variant_selected = 1;
          }
          spec->symbolic = NULL;
          kls_separator_analysis_clear(&best_separator);
          kls_separator_analysis_move(&best_separator,
                                      &spec->separator);
        } else {
          trilinos_klu_l_free_symbolic(&spec->symbolic, &spec->common);
          kls_separator_analysis_clear(&spec->separator);
        }
        amf_done = 1;
      } else {
        if (spec->status == KLS_OK && spec->symbolic != NULL) {
          trilinos_klu_l_free_symbolic(&spec->symbolic, &spec->common);
        }
        kls_separator_analysis_clear(&spec->separator);
      }
    }
    if (!amf_done && !suppress_fallback_amf && !amf_trial_disabled &&
        amf_trial_economically_relevant &&
        best_score < DBL_MAX) {
      kls_options settled_symbolic_options = *symbolic_options;
      {
        settled_symbolic_options.use_btf = best_symbolic->do_btf ? 1 : 0;
      }
      maybe_promote_symbolic_ordering(
        n, col_ptr, row_idx,
        &settled_symbolic_options,
        KLS_ORDERING_AMF,
        &best_symbolic, &best_common, &best_ordering, &best_score,
        0.95,
        &best_separator);
    }
  }

  if (aligned_btf_ammf_spec_spawned) {
    pthread_join(aligned_btf_ammf_spec_tid, NULL);
  }

  /* The representation race captures the base BTF symbolic before the
     independent minimum-fill workers settle the no-BTF ordering axis.  Do
     not later compare different representation *and* ordering choices: when
     the retained ordering changed, rebuild the captured BTF arm with that
     same ordering and keep it only if it Pareto-improves the captured base
     estimate.  The value-time factor comparison remains authoritative. */
  if (kls_generic_btf_capture_candidate != NULL &&
      kls_repeated_update_workload(options) && best_symbolic != NULL &&
      !best_symbolic->do_btf &&
      kls_generic_btf_capture_candidate->generic_btf_value_symbolic != NULL &&
      kls_generic_btf_capture_candidate->generic_btf_value_ordering !=
        best_ordering &&
      getenv("KLS_DISABLE_ALIGNED_BTF_ORDERING_CHALLENGER") == NULL) {
    kls_pattern_candidate *capture = kls_generic_btf_capture_candidate;
    trilinos_klu_l_symbolic *captured =
      capture->generic_btf_value_symbolic;
    const double captured_score = symbolic_score(captured);
    const double captured_flops = captured->est_flops;
    kls_options aligned_options = *symbolic_options;
    aligned_options.use_btf = 1;
    aligned_options.ordering = best_ordering;
    trilinos_klu_l_symbolic *aligned_symbolic = NULL;
    trilinos_klu_l_common aligned_common;
    kls_separator_analysis aligned_separator;
    memset(&aligned_separator, 0, sizeof(aligned_separator));
    int aligned_status = KLS_ERR_ANALYZE_FAILED;
    if (aligned_btf_ammf_spec_spawned &&
        best_ordering == KLS_ORDERING_AMMF &&
        aligned_btf_ammf_spec.status == KLS_OK &&
        aligned_btf_ammf_spec.symbolic != NULL) {
      aligned_symbolic = aligned_btf_ammf_spec.symbolic;
      aligned_btf_ammf_spec.symbolic = NULL;
      aligned_common = aligned_btf_ammf_spec.common;
      kls_separator_analysis_move(
        &aligned_separator, &aligned_btf_ammf_spec.separator);
      aligned_status = KLS_OK;
    } else {
      aligned_status = analyze_with_ordering(
        n, col_ptr, row_idx, &aligned_options, best_ordering,
        &aligned_symbolic, &aligned_common, &aligned_separator);
    }
    const double aligned_score = aligned_symbolic != NULL
      ? symbolic_score(aligned_symbolic) : DBL_MAX;
    const int aligned_pareto = aligned_status == KLS_OK &&
      aligned_symbolic != NULL && aligned_symbolic->do_btf &&
      aligned_symbolic->nblocks > 1u &&
      isfinite(aligned_score) && aligned_score > 0.0 &&
      isfinite(captured_score) && captured_score > 0.0 &&
      aligned_score <= captured_score &&
      aligned_symbolic->est_flops > 0.0 && captured_flops > 0.0 &&
      aligned_symbolic->est_flops <= captured_flops;
    if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
      fprintf(stderr,
              "KLS aligned BTF challenger: %s fill %.4e/%.4e "
              "flops %.4e/%.4e -> %s\n",
              kls_ordering_name(best_ordering), aligned_score,
              captured_score,
              aligned_symbolic != NULL
                ? aligned_symbolic->est_flops : DBL_MAX,
              captured_flops, aligned_pareto ? "retain" : "reject");
    }
    if (aligned_pareto) {
      trilinos_klu_l_free_symbolic(
        &capture->generic_btf_value_symbolic,
        &capture->generic_btf_value_common);
      kls_separator_analysis_clear(
        &capture->generic_btf_value_separator);
      capture->generic_btf_value_symbolic = aligned_symbolic;
      capture->generic_btf_value_common = aligned_common;
      capture->generic_btf_value_ordering = best_ordering;
      kls_separator_analysis_move(
        &capture->generic_btf_value_separator, &aligned_separator);
    } else {
      if (aligned_symbolic != NULL) {
        trilinos_klu_l_free_symbolic(
          &aligned_symbolic, &aligned_common);
      }
      kls_separator_analysis_clear(&aligned_separator);
    }
  }
  if (aligned_btf_ammf_spec_spawned) {
    if (aligned_btf_ammf_spec.symbolic != NULL) {
      trilinos_klu_l_free_symbolic(
        &aligned_btf_ammf_spec.symbolic,
        &aligned_btf_ammf_spec.common);
    }
    kls_separator_analysis_clear(
      &aligned_btf_ammf_spec.separator);
  }
  /* Generic AUTO defers nested dissection until numeric values are
     available.  kls_maybe_start_metis_race launches it only for a
     high-cost retained factor and the factor-time promotion compares actual
     fill/work plus warm refactor timings.  Running NodeND here would both
     multiply it by the orientation candidates and let an optimistic
     symbolic estimate install a pivot-inflated numeric. */

#ifdef KLS_HAVE_SCOTCH

#endif

  /* Generic AUTO already settled BTF/no-BTF once on the base graph and
     aligned every ordering candidate to that representation.  Repeating the
     same structural pass here for the winning ordering would make analysis
     cost grow with portfolio width while providing no new representation
     choice. */

  *symbolic_out = best_symbolic;
  *common_out = best_common;
  *selected_ordering_out = best_ordering;
  *score_out = best_score;
  kls_separator_analysis_move(separator_out, &best_separator);
  return KLS_OK;
}

/* The prestatic-bound structural class: a mostly-missing diagonal at
   the large-spral size band (pre2).  These matrices route to the
   pre-static exact matching at factor time and adopt its matched
   coordinates and ordering, so the base analyze's candidate
   competition (concurrent AMD+AMF, ~0.88s on pre2) is pure waste on
   the one-shot critical path.  The bounds mirror
   large_spral_candidate exactly. */
static int is_prestatic_bound_missing_diagonal_pattern(
    UF_long n, const UF_long *col_ptr, const UF_long *row_idx) {
  if (n <= 150000 || n > 750000 || col_ptr == NULL || row_idx == NULL ||
      col_ptr[n] > 8000000) {
    return 0;
  }
  if (is_large_sparse_diagonal_low_degree_pattern(n, col_ptr, row_idx)) {
    return 0;
  }
  UF_long missing = 0;
  for (UF_long col = 0; col < n; ++col) {
    int has_diag = 0;
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      if (row_idx[p] == col) {
        has_diag = 1;
        break;
      }
    }
    missing += !has_diag;
  }
  /* two structural arms, mirroring the large-spral gate's missing-
     diagonal arm exactly: the mostly-missing class (pre2) and the
     partially-missing class (rajat29/rajat30, missing ~1.5% of n but
     >=4096 rows - their prestatic adopts the matched coordinates the
     same way, so the base candidate build is equally wasted).  The
     partial arm ships only after its two blockers fell: the lockstep
     lane overflow (edacd6a) and the light-column non-lowdeg trial
     pipe routing (rajat30 21.6s piped vs sub-second serial). */
  if (missing * 2 >= n) {
    return 1;
  }
  /* partial-missing arm (missing ~1.5% of n, >=4096 rows): the
     prestatic adopts matched coordinates the same way, so the base
     candidate build is equally wasted.  Scoped to the light-column
     members whose serial matched trial is proven sub-second
     (rajat29 5.25 -> 3.2s); the denser members (rajat30, nnz/n=9.6)
     stay on the full analyze until the pipe's dense-column
     pathology is solved - probe-routed, its piped trial ground
     21.5s against a 9.7s baseline. */
  return missing >= 4096 && missing * 200 >= n &&
         col_ptr[n] < 8 * n;
}

/* Medium partial-diagonal arm of the same pre-static contract.  Every
   structurally missing diagonal is necessarily counted as weak at factor
   time, so these bounds alone imply the medium gate (at least 2,000 weak
   rows and one percent of n) without looking at numeric values.  The matched
   adoption replaces the base symbolic; a cheap NATURAL probe is therefore
   enough for the analyze-stage eligibility signal. */

/* set while the prestatic rejection redo rebuilds a real ordering so
   the probe below cannot re-fire */
static _Thread_local int kls_ps_ana_probe_disable;

/* The generic analyze tournament may run the two orientations concurrently.
   Carry a detected block-order proposal through the recursive AMD-placeholder
   call without sharing ownership between those threads.  analyze_candidate
   immediately moves this pointer into its candidate. */
static _Thread_local UF_long *kls_block_order_proposal;

static int choose_symbolic_for_pattern(UF_long n,
                                       UF_long *col_ptr,
                                       UF_long *row_idx,
                                       const kls_options *options,
                                       trilinos_klu_l_symbolic **symbolic_out,
                                       trilinos_klu_l_common *common_out,
                                       kls_ordering *selected_ordering_out,
                                       double *score_out,
                                       kls_separator_analysis *separator_out) {
  if (kls_matched_low_pair_ordering_ctx && options != NULL &&
      options->ordering == KLS_ORDERING_AUTO && options->use_btf &&
      kls_repeated_update_workload(options)) {
    /* Matching changes the directed graph: even a tiny one-block AMD form
       can hide a fragmented BTF form with materially less recurring work.
       Compare those two representations on the actual matched pattern.  The
       caller has already bounded quotient-pair work, so broader minimum-fill
       trials cannot repay another graph traversal over this lifecycle. */
    kls_options btf_options = *options;
    kls_options direct_options = *options;
    btf_options.ordering = KLS_ORDERING_AMD;
    btf_options.use_btf = 1;
    direct_options.ordering = KLS_ORDERING_AMD;
    direct_options.use_btf = 0;
    trilinos_klu_l_symbolic *btf_symbolic = NULL;
    trilinos_klu_l_symbolic *direct_symbolic = NULL;
    trilinos_klu_l_common btf_common;
    trilinos_klu_l_common direct_common;
    kls_separator_analysis btf_separator;
    kls_separator_analysis direct_separator;
    memset(&btf_separator, 0, sizeof(btf_separator));
    memset(&direct_separator, 0, sizeof(direct_separator));
    struct kls_amf_spec_job direct_job;
    pthread_t direct_thread;
    memset(&direct_job, 0, sizeof(direct_job));
    direct_job.n = n;
    direct_job.col_ptr = col_ptr;
    direct_job.row_idx = row_idx;
    direct_job.ordering = KLS_ORDERING_AMD;
    direct_job.options = direct_options;
    const int trace_low_pair =
      getenv("KLS_TRACE_FACTOR_PHASES") != NULL;
    const double low_pair_start =
      trace_low_pair ? kls_now_seconds() : 0.0;
    /* A bounded compact match has already paid for a complete transversal.
       Below this fixed quotient-work budget, constructing a second no-BTF
       AMD graph costs a material fraction of the whole numeric lifecycle,
       while the BTF arm still has to pass the ordinary fill/work and numeric
       acceptance gates.  Larger matched graphs retain the parallel pair:
       their direct arm can remove enough recurring work to repay itself. */
    const int skip_direct = kls_compact_missing_match_ordering_ctx &&
      kls_matched_low_pair_work_ctx > 0.0 &&
      kls_matched_low_pair_work_ctx <= 5.0e5;
    const int direct_spawned = !skip_direct && options->threads >= 2 &&
      pthread_create(&direct_thread, NULL, kls_amf_spec_main,
                     &direct_job) == 0;
    const int btf_status = analyze_with_ordering(
      n, col_ptr, row_idx, &btf_options, KLS_ORDERING_AMD,
      &btf_symbolic, &btf_common, &btf_separator);
    const double low_pair_btf_done =
      trace_low_pair ? kls_now_seconds() : 0.0;
    int direct_status = KLS_ERR_FACTOR_FAILED;
    if (direct_spawned) {
      pthread_join(direct_thread, NULL);
      direct_status = direct_job.status;
      direct_symbolic = direct_job.symbolic;
      direct_common = direct_job.common;
      kls_separator_analysis_move(&direct_separator,
                                  &direct_job.separator);
    } else if (!skip_direct) {
      direct_status = analyze_with_ordering(
        n, col_ptr, row_idx, &direct_options, KLS_ORDERING_AMD,
        &direct_symbolic, &direct_common, &direct_separator);
    }
    const double low_pair_direct_done =
      trace_low_pair ? kls_now_seconds() : 0.0;
    const int btf_valid = btf_status == KLS_OK && btf_symbolic != NULL;
    const int direct_valid =
      direct_status == KLS_OK && direct_symbolic != NULL;
    if (btf_valid || direct_valid) {
      const double btf_score =
        btf_valid ? symbolic_score(btf_symbolic) : DBL_MAX;
      const double direct_score =
        direct_valid ? symbolic_score(direct_symbolic) : DBL_MAX;
      const int btf_dispatch_heavy = btf_valid &&
        btf_symbolic->nblocks > 1u &&
        btf_symbolic->nblocks >= (n + 7u) / 8u;
      const int btf_dispatch_margin = !btf_dispatch_heavy ||
        (isfinite(btf_score) && btf_score <= 0.80 * direct_score &&
         btf_symbolic->est_flops <= 0.80 * direct_symbolic->est_flops) ||
        (isfinite(btf_score) && btf_score <= 0.90 * direct_score &&
         btf_symbolic->est_flops <= 0.65 * direct_symbolic->est_flops);
      const int btf_fragmented_pareto = btf_valid && direct_valid &&
        btf_symbolic->do_btf && btf_symbolic->nblocks > 1u &&
        btf_symbolic->structural_rank == n &&
        btf_symbolic->est_flops > 0.0 &&
        direct_symbolic->est_flops > 0.0 &&
        btf_dispatch_margin &&
        ((isfinite(btf_score) && btf_score <= 0.95 * direct_score &&
          btf_symbolic->est_flops <= direct_symbolic->est_flops) ||
         (isfinite(btf_score) && btf_score <= direct_score &&
          btf_symbolic->est_flops <= 0.80 * direct_symbolic->est_flops));
      const int choose_btf = !direct_valid || btf_fragmented_pareto;
      if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
        fprintf(stderr,
                "KLS matched low-pair representation: btf %.4e/%.4e "
                "blocks=%ld direct %.4e/%.4e -> %s "
                "btf=%.3fms total=%.3fms\n",
                btf_score,
                btf_valid ? btf_symbolic->est_flops : -1.0,
                btf_valid ? (long)btf_symbolic->nblocks : -1L,
                direct_score,
                direct_valid ? direct_symbolic->est_flops : -1.0,
                choose_btf ? "btf" : "direct",
                1e3 * (low_pair_btf_done - low_pair_start),
                1e3 * (low_pair_direct_done - low_pair_start));
      }
      if (choose_btf) {
        if (direct_valid) {
          trilinos_klu_l_free_symbolic(&direct_symbolic, &direct_common);
        }
        kls_separator_analysis_clear(&direct_separator);
        *symbolic_out = btf_symbolic;
        *common_out = btf_common;
        *score_out = btf_score;
        kls_separator_analysis_move(separator_out, &btf_separator);
      } else {
        if (btf_valid) {
          trilinos_klu_l_free_symbolic(&btf_symbolic, &btf_common);
        }
        kls_separator_analysis_clear(&btf_separator);
        *symbolic_out = direct_symbolic;
        *common_out = direct_common;
        *score_out = direct_score;
        kls_separator_analysis_move(separator_out, &direct_separator);
      }
      *selected_ordering_out = KLS_ORDERING_AMD;
      return KLS_OK;
    }
    kls_separator_analysis_clear(&btf_separator);
    kls_separator_analysis_clear(&direct_separator);
  }

  int generic_match_resource_bound =
    n <= 20000u || (n <= 150000u && col_ptr[n] <= 1500000u);
#ifdef KLS_HAVE_SPRAL_SCALING
  generic_match_resource_bound = generic_match_resource_bound ||
    (n <= 750000u && col_ptr[n] <= 8000000u);
#endif
  const UF_long generic_pattern_diagonal =
    options != NULL &&
        options->ordering == KLS_ORDERING_AUTO &&
        options->static_pivoting && options->use_btf &&
        kls_repeated_update_workload(options) && n >= 3000u &&
        generic_match_resource_bound
      ? count_pattern_diagonal(n, col_ptr, row_idx)
      : n;
  const int generic_guaranteed_inline_match =
    generic_pattern_diagonal < n &&
    (2.0 * (double)(n - generic_pattern_diagonal) >= (double)n ||
     (double)(n - generic_pattern_diagonal) *
         (1.0 + (double)options->expected_refactorizations) >=
       4.0 * (double)col_ptr[n]);
  /* The cheap NATURAL probe avoids an ordering that value-aware matching
     normally replaces.  The compact proposal still has to prove a dominant
     SCC here, then match coverage and numeric quality at factor time. */
  if (options != NULL && options->ordering == KLS_ORDERING_AUTO &&
      getenv("KLS_DISABLE_AUTO_AMD_SHORTCUT") == NULL &&
      1 &&
      is_medium_sparse_full_diagonal_amd_pattern(n, col_ptr, row_idx)) {
    kls_options amd_options = *options;
    amd_options.ordering = KLS_ORDERING_AMD;
    const int status = choose_symbolic_for_pattern(
      n, col_ptr, row_idx, &amd_options, symbolic_out, common_out,
      selected_ordering_out, score_out, separator_out);
    if (status == KLS_OK &&
        kls_medium_sparse_full_diagonal_amd_symbolic_contract(
          n, col_ptr, row_idx, *symbolic_out)) {
      return KLS_OK;
    }
    if (*symbolic_out != NULL) {
      trilinos_klu_l_free_symbolic(symbolic_out, common_out);
    }
    kls_separator_analysis_clear(separator_out);
    *score_out = DBL_MAX;
    /* A topology proposal is never authoritative: allocation failure or a
       rejected AMD symbolic resumes the complete AUTO selection below. */
  }
  if (options != NULL && options->ordering == KLS_ORDERING_AUTO &&
      options->static_pivoting && options->use_btf &&
      !kls_prestatic_ordering_ctx && !kls_ps_ana_probe_disable &&
      getenv("KLS_DISABLE_PS_ANA_SKIP") == NULL &&
      (generic_guaranteed_inline_match ||
       is_prestatic_bound_missing_diagonal_pattern(n, col_ptr, row_idx))) {
    /* A NATURAL-ordering analyze gives the pre-static gate the signal it
       needs at ~0.2s instead of the 0.88s candidate competition whose
       winner the adoption replaces anyway.  Existing policies retain BTF
       so maxblock can reject decomposable inputs.  Compact missing-diagonal
       candidates omit that disposable decomposition and prove rank plus a
       dominant SCC on the matched symbolic below.  Any rejection rebuilds
       a real ordering at factor entry (kls_ps_ana_probe_disable). */
    kls_options probe_options = *options;
    if (generic_guaranteed_inline_match) {
      /* The compact matcher and matched-symbolic gate below provide the
         transversal/rank and dominant-block proofs.  Avoid constructing a
         BTF decomposition for coordinates that cannot survive adoption. */
      probe_options.use_btf = 0;
    }
    trilinos_klu_l_symbolic *probe_sym = NULL;
    trilinos_klu_l_common probe_common;
    kls_separator_analysis probe_sep;
    memset(&probe_sep, 0, sizeof(probe_sep));
    if (analyze_with_ordering(n, col_ptr, row_idx, &probe_options,
                              KLS_ORDERING_NATURAL, &probe_sym,
                              &probe_common, &probe_sep) == KLS_OK &&
        probe_sym != NULL) {
      if (probe_sym->maxblock * 4 >= n * 3) {
        if (getenv("KLS_TRACE_PRESTATIC") != NULL) {
          fprintf(stderr,
                  "KLS choose: prestatic-bound probe adopted "
                  "(maxblock=%ld)\n",
                  (long)probe_sym->maxblock);
        }
        *symbolic_out = probe_sym;
        *common_out = probe_common;
        *selected_ordering_out = KLS_ORDERING_NATURAL;
        *score_out = symbolic_score(probe_sym);
        kls_separator_analysis_move(separator_out, &probe_sep);
        return KLS_OK;
      }
      trilinos_klu_l_free_symbolic(&probe_sym, &probe_common);
      kls_separator_analysis_clear(&probe_sep);
    } else {
      kls_separator_analysis_clear(&probe_sep);
    }
  }
  if ( options != NULL &&
      options->ordering == KLS_ORDERING_AUTO && options->use_btf &&
      kls_repeated_update_workload(options) &&
      !kls_prestatic_ordering_ctx &&
      !generic_guaranteed_inline_match &&
      getenv("KLS_ENABLE_BOUNDED_NO_BTF_AMD_SHORTCUT") != NULL) {
    const double numeric_horizon = 1.0 +
      (double)options->expected_refactorizations +
      (double)options->expected_solves;
    double quotient_pair_work = 0.0;
    for (UF_long col = 0u; col < n; ++col) {
      const double degree =
        (double)(col_ptr[col + 1u] - col_ptr[col]);
      quotient_pair_work += degree * degree;
    }
    const double bounded_numeric_lifecycle_budget = 1.0e8;
    /* This no-BTF bootstrap remains available as a diagnostic.  Production
       AUTO starts from the ordinary AMD/BTF representation: when the
       bootstrap found a strict BTF Pareto winner, analyzing no-BTF first
       duplicated symbolic work and moved row-worker preparation into the
       first refactor.  Close representation choices remain covered by the
       generic BTF/no-BTF tournament at larger work budgets.

       A rejected symbolic is pure setup overhead.  Reserve most of the
       bounded lifecycle budget for the numeric work it is meant to improve;
       at most one quarter may be spent speculating on this representation. */
    const double speculative_ordering_budget =
      0.25 * bounded_numeric_lifecycle_budget;
    if (numeric_horizon > 0.0 && isfinite(numeric_horizon) &&
        isfinite(quotient_pair_work) &&
        quotient_pair_work * numeric_horizon <=
          speculative_ordering_budget) {
      /* The complete requested lifecycle contains less quotient-graph work
         than one large ordering trial.  Probe the stable one-block AMD
         representation directly, and retain it only when its computed
         fill and arithmetic fit the same absolute lifecycle budget.  This
         replaces dimension-based tiny-matrix shortcuts with measured work;
         a rejected proposal resumes the full BTF/ordering tournament. */
      kls_options amd_options = *options;
      amd_options.ordering = KLS_ORDERING_AMD;
      amd_options.use_btf = 0;
      trilinos_klu_l_symbolic *amd_symbolic = NULL;
      trilinos_klu_l_common amd_common;
      kls_separator_analysis amd_separator;
      memset(&amd_separator, 0, sizeof(amd_separator));
      const int amd_status = analyze_with_ordering(
        n, col_ptr, row_idx, &amd_options, KLS_ORDERING_AMD,
        &amd_symbolic, &amd_common, &amd_separator);
      if (amd_status == KLS_OK && amd_symbolic != NULL) {
        const double amd_score = symbolic_score(amd_symbolic);
        const int bounded =
          !amd_symbolic->do_btf && amd_symbolic->nblocks == 1u &&
          amd_symbolic->maxblock == n && amd_symbolic->est_flops > 0.0 &&
          isfinite(amd_score) && amd_score > 0.0 &&
          amd_symbolic->est_flops * numeric_horizon <=
            bounded_numeric_lifecycle_budget &&
          amd_score * numeric_horizon <= 4.0e8 &&
          amd_symbolic->lnz <= 2.0 * amd_symbolic->unz &&
          amd_symbolic->unz <= 2.0 * amd_symbolic->lnz;
        if (bounded) {
          /* The inexpensive one-block proposal is a setup shortcut, not a
             proof that BTF loses over the requested numeric lifecycle.  If
             the ordinary AMD/BTF symbolic is work-Pareto and its projected
             recurring savings repay both its quotient-graph analysis and a
             numeric trial, retain it for the value-time residual tournament.
             This compares computed representations and caller work only; it
             does not classify the input by dimensions or sparsity shape. */
          kls_pattern_candidate *capture =
            kls_generic_btf_capture_candidate;
          if (capture != NULL &&
              capture->generic_btf_value_symbolic == NULL &&
              getenv("KLS_DISABLE_GENERIC_BTF_VALUE_SELECTION") == NULL) {
            kls_options btf_options = *options;
            btf_options.ordering = KLS_ORDERING_AMD;
            btf_options.use_btf = 1;
            trilinos_klu_l_symbolic *btf_symbolic = NULL;
            trilinos_klu_l_common btf_common;
            kls_separator_analysis btf_separator;
            memset(&btf_separator, 0, sizeof(btf_separator));
            const int btf_status = analyze_with_ordering(
              n, col_ptr, row_idx, &btf_options, KLS_ORDERING_AMD,
              &btf_symbolic, &btf_common, &btf_separator);
            const double btf_score = btf_symbolic != NULL
              ? symbolic_score(btf_symbolic) : DBL_MAX;
            const double recurring_frames =
              (double)options->expected_refactorizations;
            const double solve_frames =
              (double)options->expected_solves;
            const double saved_fill = isfinite(btf_score)
              ? fmax(0.0, amd_score - btf_score) : 0.0;
            const double saved_flops = btf_symbolic != NULL
              ? fmax(0.0,
                     amd_symbolic->est_flops - btf_symbolic->est_flops)
              : 0.0;
            const double projected_savings =
              recurring_frames * (saved_fill + saved_flops) +
              solve_frames * saved_fill;
            const double comparison_work = btf_symbolic != NULL
              ? quotient_pair_work + btf_score : DBL_MAX;
            const int retain_btf = btf_status == KLS_OK &&
              btf_symbolic != NULL && btf_symbolic->do_btf &&
              btf_symbolic->nblocks > 1u &&
              btf_symbolic->maxblock < n &&
              isfinite(btf_score) && btf_score > 0.0 &&
              btf_score <= amd_score &&
              btf_symbolic->est_flops > 0.0 &&
              btf_symbolic->est_flops <= amd_symbolic->est_flops &&
              isfinite(projected_savings) &&
              isfinite(comparison_work) && comparison_work > 0.0 &&
              projected_savings >= 2.0 * comparison_work;
            if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
              fprintf(stderr,
                      "KLS bounded AMD BTF challenger: fill %.4e/%.4e "
                      "flops %.4e/%.4e savings %.4e cost %.4e -> %s\n",
                      btf_score, amd_score,
                      btf_symbolic != NULL
                        ? btf_symbolic->est_flops : DBL_MAX,
                      amd_symbolic->est_flops, projected_savings,
                      comparison_work, retain_btf ? "retain" : "reject");
            }
            if (retain_btf) {
              capture->generic_btf_value_symbolic = btf_symbolic;
              capture->generic_btf_value_common = btf_common;
              capture->generic_btf_value_ordering = KLS_ORDERING_AMD;
              kls_separator_analysis_move(
                &capture->generic_btf_value_separator, &btf_separator);
            } else {
              if (btf_symbolic != NULL) {
                trilinos_klu_l_free_symbolic(&btf_symbolic, &btf_common);
              }
              kls_separator_analysis_clear(&btf_separator);
            }
          }
          *symbolic_out = amd_symbolic;
          *common_out = amd_common;
          *selected_ordering_out = KLS_ORDERING_AMD;
          *score_out = amd_score;
          kls_separator_analysis_move(separator_out, &amd_separator);
          return KLS_OK;
        }
        trilinos_klu_l_free_symbolic(&amd_symbolic, &amd_common);
      }
      kls_separator_analysis_clear(&amd_separator);
    }
  }
  if (options != NULL && options->ordering == KLS_ORDERING_AUTO &&
      options->static_pivoting && n <= 200000 && n >= 64 &&
      col_ptr[n] <= 2000000) {
    /* block-structured systems (TSOPF, case9) adopt their own
       ordering at factor time via the block trial; ordering scoring
       here is one-shot overhead for them.  Pattern-only detector. */
    UF_long *bs_perm = NULL, *bs_comp = NULL;
    if (kls_build_block_structured_order(n, col_ptr, row_idx, &bs_perm,
                                         &bs_comp)) {
      free(bs_comp);
      free(kls_block_order_proposal);
      kls_block_order_proposal = bs_perm;
      kls_options amd_opts = *options;
      amd_opts.ordering = KLS_ORDERING_AMD;
      /* BTF's maximum transversal is pathological on these patterns
         (case9: 20ms of a 23ms analyze at n=14K) and the block
         ordering supersedes the block-triangular form anyway */
      amd_opts.use_btf = 0;
      if (getenv("KLS_TRACE_PREDICTED") != NULL) {
        fprintf(stderr, "KLS choose: block structure detected -"
                " AMD placeholder, no BTF, scoring skipped\n");
      }
      return choose_symbolic_for_pattern(n, col_ptr, row_idx, &amd_opts,
                                         symbolic_out, common_out,
                                         selected_ordering_out, score_out,
                                         separator_out);
    }
  }
  UF_long *bp = NULL, *bq = NULL, *br = NULL, *bwork = NULL;
  struct kls_btf_stash_s saved = kls_btf_stash;
  if (n >= 30000 && options != NULL && options->use_btf &&
      getenv("KLS_ENABLE_BTF_REUSE") != NULL) {
    /* measured NET LOSS with the concurrent AMF spec (twotone ana
       1.5 -> 2.8): pre-reuse each candidate's internal BTF ran hidden
       inside the parallel AMD||AMF pair, while the shared stash
       serializes one BTF ahead of both.  AMD-only paths win hugely
       (0.755 -> 0.043) - re-enable if the spec is ever off. */
    bp = (UF_long *)malloc((size_t)n * sizeof(UF_long));
    bq = (UF_long *)malloc((size_t)n * sizeof(UF_long));
    br = (UF_long *)malloc(((size_t)n + 1u) * sizeof(UF_long));
    bwork = (UF_long *)malloc(5u * (size_t)n * sizeof(UF_long));
    if (bp != NULL && bq != NULL && br != NULL && bwork != NULL) {
      double btf_work = 0.0;
      UF_long nmatch = 0;
      UF_long nblocks = trilinos_btf_l_order(n, col_ptr, row_idx, 0.0,
                                             &btf_work, bp, bq, br, &nmatch,
                                             bwork);
      if (nblocks > 0) {
        if (nmatch < n) {
          UF_long k;
          for (k = 0; k < n; ++k) {
            bq[k] = bq[k] < 0 ? -bq[k] - 2 : bq[k];
          }
        }
        kls_btf_stash.P = bp;
        kls_btf_stash.Q = bq;
        kls_btf_stash.R = br;
        kls_btf_stash.nblocks = nblocks;
        kls_btf_stash.rank = nmatch;
        kls_btf_stash.n = n;
        kls_btf_stash.active = 1;
      }
    }
    free(bwork);
    bwork = NULL;
  }
  int status = kls_choose_symbolic_inner(n, col_ptr, row_idx, options,
                                         symbolic_out, common_out,
                                         selected_ordering_out, score_out,
                                         separator_out);
  kls_btf_stash = saved;
  free(bp);
  free(bq);
  free(br);
  return status;
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
  if (options->ordering < KLS_ORDERING_AUTO ||
      options->ordering > KLS_ORDERING_AMF3) {
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
  if (options->backend < KLS_BACKEND_AUTO ||
      options->backend > KLS_BACKEND_SERIAL) {
    return 0;
  }
  if (options->expected_refactorizations < 0 ||
      options->expected_solves < 0 ||
      (options->record_tiny_solve_timing != 0 &&
       options->record_tiny_solve_timing != 1)) {
    return 0;
  }
  return 1;
}

static int kls_normalize_options(kls_options *normalized,
                                 const kls_options *options,
                                 kls_input_format input_format) {
  if (normalized == NULL) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  kls_default_options(normalized);
  if (options != NULL) {
    /* backend was appended in the original layout's tail padding, so the
       public sizeof(kls_options) did not grow.  abi_version occupies older
       interior padding and distinguishes a current initialized object from
       an old object whose tail padding merely resembles a backend value. */
    const size_t legacy_size = offsetof(kls_options, backend);
    if (options->struct_size < legacy_size) {
      return KLS_ERR_INVALID_ARGUMENT;
    }
    uint32_t supplied_abi = 0u;
    if (options->struct_size >=
        offsetof(kls_options, abi_version) + sizeof(supplied_abi)) {
      memcpy(&supplied_abi,
             (const unsigned char *)options +
               offsetof(kls_options, abi_version),
             sizeof(supplied_abi));
    }
    const int has_backend =
      supplied_abi == KLS_OPTIONS_ABI_VERSION &&
      options->struct_size >=
        offsetof(kls_options, backend) + sizeof(options->backend);
    const size_t advertised_size = has_backend
      ? options->struct_size : legacy_size;
    const size_t copy_size = advertised_size < sizeof(*normalized)
      ? advertised_size : sizeof(*normalized);
    memcpy(normalized, options, copy_size);
    normalized->struct_size = sizeof(*normalized);
    normalized->abi_version = KLS_OPTIONS_ABI_VERSION;
    if (!has_backend) {
      normalized->backend = KLS_BACKEND_AUTO;
    }
  }
  if (!validate_options(normalized)) {
    return KLS_ERR_INVALID_ARGUMENT;
  }
  if (normalized->backend == KLS_BACKEND_SERIAL) {
    if (normalized->threads != 1) {
      return KLS_ERR_UNSUPPORTED;
    }
    if (normalized->ordering == KLS_ORDERING_AUTO) {
      normalized->ordering = KLS_ORDERING_AMD;
    }
    if (normalized->ordering != KLS_ORDERING_AMD &&
        normalized->ordering != KLS_ORDERING_COLAMD &&
        normalized->ordering != KLS_ORDERING_NATURAL) {
      return KLS_ERR_UNSUPPORTED;
    }
    if (normalized->orientation == KLS_ORIENTATION_AUTO) {
      normalized->orientation = input_format == KLS_INPUT_CSR
        ? KLS_ORIENTATION_TRANSPOSE : KLS_ORIENTATION_NORMAL;
    }
  }
  return KLS_OK;
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
  if (candidate->generic_btf_value_symbolic != NULL) {
    trilinos_klu_l_free_symbolic(
      &candidate->generic_btf_value_symbolic,
      &candidate->generic_btf_value_common);
  }
  kls_separator_analysis_clear(&candidate->separator);
  kls_separator_analysis_clear(&candidate->generic_btf_value_separator);
  free(candidate->col_ptr);
  free(candidate->row_idx);
  free(candidate->input_to_csc);
  free(candidate->values);
  free(candidate->block_order_perm);
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
  {
    for (int64_t p = 0; p < nnz; ++p) {
      const int64_t row = read_index(index_type, row_idx, p) - base;
      if (row < 0 || row >= n) {
        free_candidate(candidate);
        return KLS_ERR_INVALID_ARGUMENT;
      }
      candidate->row_idx[p] = (UF_long)row;
    }
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
  {
    /* Generic candidates carry only their copied pattern.  In particular,
       do not compute or publish any of the archived input-family contexts:
       the ordering tournament consumes realized symbolic scores directly. */
    free(kls_block_order_proposal);
    kls_block_order_proposal = NULL;
    kls_amf3_span_variant_selected = 0;
    if (candidate->generic_btf_value_symbolic != NULL) {
      trilinos_klu_l_free_symbolic(
        &candidate->generic_btf_value_symbolic,
        &candidate->generic_btf_value_common);
    }
    kls_separator_analysis_clear(
      &candidate->generic_btf_value_separator);
    kls_pattern_candidate *saved_capture =
      kls_generic_btf_capture_candidate;
    kls_generic_btf_capture_candidate = candidate;
    const int status = choose_symbolic_for_pattern(
      candidate->n, candidate->col_ptr, candidate->row_idx, options,
      &candidate->symbolic, &candidate->common,
      &candidate->selected_ordering, &candidate->score,
      &candidate->separator);
    kls_generic_btf_capture_candidate = saved_capture;
    candidate->generic_amf3_span_variant_selected =
      status == KLS_OK &&
      candidate->selected_ordering == KLS_ORDERING_AMF3 &&
      kls_amf3_span_variant_selected;
    kls_amf3_span_variant_selected = 0;
    free(candidate->block_order_perm);
    candidate->block_order_perm = kls_block_order_proposal;
    kls_block_order_proposal = NULL;
    return status;
  }
}

struct kls_candidate_analyze_job {
  kls_pattern_candidate *candidate;
  const kls_options *options;
  int defer_nd;   /* propagate the caller's defer-ND mode (TL) */
  int status;
};

static void *kls_candidate_analyze_main(void *arg) {
  struct kls_candidate_analyze_job *job =
    (struct kls_candidate_analyze_job *)arg;
  kls_analyze_defer_nd = job->defer_nd;
  job->status = analyze_candidate(job->candidate, job->options);
  kls_analyze_defer_nd = 0;
  return NULL;
}

static int auto_orientation_prefers_transpose(UF_long n) {
  /* On small and medium SPICE-like matrices, the second symbolic analysis
     usually costs more than the normal-vs-transpose fill estimate saves. */
  return n <= 30000;
}

/* Compare representations over the caller's stated lifecycle, rather than
   treating symbolic fill as the whole cost.  A derived orientation carries
   input_to_csc and must gather every public value update into that frame.
   One gather consumes an input value and index and publishes one value; in
   fill-equivalent units that is conservatively 3/8 of one sparse numeric
   entry.  The coefficient follows from bytes moved, not matrix dimensions or
   a family classifier.  With no repeated-work hint this reduces exactly to
   the historical symbolic score. */
static double kls_candidate_lifecycle_score(
  const kls_pattern_candidate *candidate,
  const kls_options *options) {
  if (candidate == NULL || options == NULL ||
      !isfinite(candidate->score) || candidate->score <= 0.0) {
    return DBL_MAX;
  }
  const double refactors =
    options->expected_refactorizations > 0
      ? (double)options->expected_refactorizations : 0.0;
  double score = (1.0 + refactors) * candidate->score;
  if (candidate->generic_btf_value_symbolic != NULL) {
    const double alternative_score =
      symbolic_score(candidate->generic_btf_value_symbolic);
    if (isfinite(alternative_score) && alternative_score > 0.0) {
      /* A retained representation is a real factor-time candidate, so the
         orientation decision must price the cheapest representation it can
         actually publish.  Otherwise a no-BTF bootstrap can choose one
         coordinate frame while its BTF challenger proves the other frame is
         cheaper.  Numeric fill and residual gates remain authoritative. */
      score = fmin(score, (1.0 + refactors) * alternative_score);
    }
  }
  if (refactors > 0.0 && candidate->input_to_csc != NULL) {
    score += refactors * (3.0 / 8.0) * (double)candidate->nnz;
  }
  return isfinite(score) ? score : DBL_MAX;
}

/* Certify the cheapest exact case in which BTF cannot change the pattern:
   the input diagonal is a perfect structural matching and the directed graph
   is one strongly connected component.  The already-built transpose makes
   the two reachability walks linear.  Missing even one diagonal immediately
   preserves ordinary BTF, avoiding a second maximum-matching pass. */
static int kls_identity_matching_single_scc(
  const kls_pattern_candidate *normal,
  const kls_pattern_candidate *transpose) {
  if (normal == NULL || transpose == NULL || normal->n == 0u ||
      normal->n != transpose->n || normal->nnz != transpose->nnz ||
      normal->col_ptr == NULL || normal->row_idx == NULL ||
      transpose->col_ptr == NULL || transpose->row_idx == NULL) {
    return 0;
  }

  const UF_long n = normal->n;
  for (UF_long col = 0u; col < n; ++col) {
    int has_diagonal = 0;
    for (UF_long p = normal->col_ptr[col];
         p < normal->col_ptr[col + 1u]; ++p) {
      const UF_long row = normal->row_idx[p];
      if (row >= n) {
        return 0;
      }
      if (row == col) {
        has_diagonal = 1;
        break;
      }
    }
    if (!has_diagonal) {
      return 0;
    }
  }

  unsigned char *seen = (unsigned char *)calloc((size_t)n, sizeof(*seen));
  UF_long *stack = (UF_long *)malloc((size_t)n * sizeof(*stack));
  if (seen == NULL || stack == NULL) {
    free(seen);
    free(stack);
    return 0;
  }

  int connected = 1;
  for (int pass = 0; pass < 2 && connected; ++pass) {
    memset(seen, 0, (size_t)n * sizeof(*seen));
    UF_long top = 0u;
    UF_long reached = 1u;
    seen[0] = 1u;
    stack[top++] = 0u;
    while (top > 0u) {
      const UF_long vertex = stack[--top];
      const kls_pattern_candidate *graph = pass == 0 ? normal : transpose;
      for (UF_long p = graph->col_ptr[vertex];
           p < graph->col_ptr[vertex + 1u]; ++p) {
        const UF_long next = graph->row_idx[p];
        if (next >= n) {
          connected = 0;
          break;
        }
        if (!seen[next]) {
          seen[next] = 1u;
          stack[top++] = next;
          reached++;
        }
      }
    }
    connected = reached == n;
    if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
      fprintf(stderr,
              "KLS identity SCC certificate: pass %d reached %ld/%ld\n",
              pass, (long)reached, (long)n);
    }
  }

  free(seen);
  free(stack);
  return connected;
}

static int kls_candidate_symmetric_union_is_bounded(
  const kls_pattern_candidate *normal,
  const kls_pattern_candidate *transpose);

static int select_candidate_inner(kls_pattern_candidate *normal,
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

  {
    goto generic_orientation_comparison;
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
generic_orientation_comparison:
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
  /* A bounded A+A' stream makes transpose analysis unable to recover its
     recurring gather cost from structural novelty alone.  Analyze normal
     once; its full ordering tournament and later numeric guards remain in
     force.  Exact symmetry is the zero-expansion member of this certificate. */
  if (kls_candidate_symmetric_union_is_bounded(normal, transpose)) {
    /* Preserve the caller's compressed frame.  For CSC that candidate is
       normal; for CSR it is transpose.  Both avoid an update gather map. */
    kls_pattern_candidate *direct = normal->input_to_csc == NULL
      ? normal : (transpose->input_to_csc == NULL ? transpose : normal);
    int status = analyze_candidate(direct, options);
    if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
      fprintf(stderr,
              "KLS select: bounded symmetric union, single analyze\n");
    }
    if (status == KLS_OK) {
      *chosen_out = direct;
      return KLS_OK;
    }
    return status;
  }
  /* Both orientations run the identical analysis; on the classes where
     the nested-dissection trials dominate the cost, run the transpose
     side on a worker thread - the comparison and results are unchanged,
     only the wall time halves. */
  int normal_status;
  int transpose_status;
  {
    struct kls_candidate_analyze_job transpose_job;
    transpose_job.candidate = transpose;
    transpose_job.options = options;
    transpose_job.defer_nd = kls_analyze_defer_nd;
    transpose_job.status = KLS_ERR_FACTOR_FAILED;
    pthread_t transpose_thread;
    int threaded = 0;
    if (((kls_repeated_update_workload(options) &&
          options->threads >= 2) ||
         normal->n >= 200000 ||
         (normal->n >= 30000 &&
          options->ordering == KLS_ORDERING_AUTO)) &&
        pthread_create(&transpose_thread, NULL, kls_candidate_analyze_main,
                       &transpose_job) == 0) {
      /* The two sides are independent full analyzes.  A repeated lifecycle
         may compare several ordering/BTF candidates even for a cache-sized
         factor, so use its requested worker budget to overlap orientations;
         the one-shot path retains the small-input spawn bound.  At larger
         sizes the per-side cost (twotone: 1.5-2.9s at n=121k) dwarfs the
         spawn.  The 30k-200k one-shot extension covers AUTO only: an
         explicit-METIS smoke matrix at n=30600 changes its separator verdict
         under threaded analyze, unresolved. */
      threaded = 1;
    }
    normal_status = analyze_candidate(normal, options);
    if (threaded) {
      pthread_join(transpose_thread, NULL);
      transpose_status = transpose_job.status;
    } else {
      transpose_status = analyze_candidate(transpose, options);
    }
  }
  if (normal_status == KLS_OK && transpose_status == KLS_OK) {
    const double normal_lifecycle =
      kls_candidate_lifecycle_score(normal, options);
    const double transpose_lifecycle =
      kls_candidate_lifecycle_score(transpose, options);
    *chosen_out = transpose_lifecycle < normal_lifecycle
      ? transpose : normal;
    if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
      fprintf(stderr,
              "KLS orientation lifecycle: normal %.4e transpose %.4e "
              "(raw %.4e/%.4e, refactors=%lld) -> %s\n",
              normal_lifecycle, transpose_lifecycle,
              normal->score, transpose->score,
              (long long)options->expected_refactorizations,
              (*chosen_out)->orientation == KLS_ORIENTATION_NORMAL
                ? "normal" : "transpose");
    }
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

/* Prediction on an unsymmetric pattern stores the structural union A+A'.
   Count the already-canonical CSC overlap to bound that additional stream.
   A one-thirty-second expansion is small enough to amortize across the
   recurring factor, while a larger union leaves ordinary pivoted KLU free to
   discover a narrower LU pattern.  This is a representation-cost certificate,
   independent of matrix dimensions or provenance. */
static int kls_candidate_symmetric_union_is_bounded(
  const kls_pattern_candidate *normal,
  const kls_pattern_candidate *transpose) {
  if (normal == NULL || transpose == NULL || normal->n != transpose->n ||
      normal->nnz != transpose->nnz || normal->col_ptr == NULL ||
      normal->row_idx == NULL || transpose->col_ptr == NULL ||
      transpose->row_idx == NULL) {
    return 0;
  }
  UF_long missing_from_transpose = 0u;
  const UF_long allowance = normal->nnz / 32u + 64u;
  for (UF_long col = 0u; col < normal->n; ++col) {
    UF_long p = normal->col_ptr[col];
    UF_long q = transpose->col_ptr[col];
    const UF_long pend = normal->col_ptr[col + 1u];
    const UF_long qend = transpose->col_ptr[col + 1u];
    while (p < pend && q < qend) {
      const UF_long prow = normal->row_idx[p];
      const UF_long qrow = transpose->row_idx[q];
      if (prow == qrow) {
        p++;
        q++;
      } else if (prow < qrow) {
        missing_from_transpose++;
        p++;
      } else {
        q++;
      }
      if (missing_from_transpose > allowance) {
        return 0;
      }
    }
    missing_from_transpose += pend - p;
    if (missing_from_transpose > allowance) {
      return 0;
    }
  }
  return 1;
}

static int select_candidate(kls_pattern_candidate *normal,
                            kls_pattern_candidate *transpose,
                            const kls_options *options,
                            kls_pattern_candidate **chosen_out) {
  {
    kls_options direct_options;
    const kls_options *generic_options = options;
    const int btf_is_provably_redundant =
      options != NULL && options->use_btf &&
      options->ordering == KLS_ORDERING_AUTO &&
      kls_identity_matching_single_scc(normal, transpose);
    if (btf_is_provably_redundant) {
      direct_options = *options;
      direct_options.use_btf = 0;
      generic_options = &direct_options;
      if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
        fprintf(stderr,
                "KLS select: identity matching + one SCC; skip BTF\n");
      }
    }
    const int exact_pattern_symmetry =
      normal != NULL && transpose != NULL &&
      normal->n == transpose->n && normal->nnz == transpose->nnz &&
      memcmp(normal->col_ptr, transpose->col_ptr,
             ((size_t)normal->n + 1u) * sizeof(UF_long)) == 0 &&
      memcmp(normal->row_idx, transpose->row_idx,
             (size_t)normal->nnz * sizeof(UF_long)) == 0;
    const int bounded_symmetric_union = exact_pattern_symmetry ||
      kls_candidate_symmetric_union_is_bounded(normal, transpose);
    /* Preserve the representation certificate even when the symbolic ND
       challenger is not selected.  Factor time can then distinguish a
       bounded symbolic comparison from an unbounded case that requires a
       real numeric trial; existing consumers also require the selected-ND
       flag, so this does not broaden prediction on the incumbent. */
    if (normal != NULL) {
      normal->generic_nd_bounded_symmetric_union = bounded_symmetric_union;
    }
    if (transpose != NULL) {
      transpose->generic_nd_bounded_symmetric_union = bounded_symmetric_union;
    }
    kls_generic_nd_trial overlapped_nd;
    memset(&overlapped_nd, 0, sizeof(overlapped_nd));
    kls_pattern_candidate *overlap_candidate = NULL;
    if (bounded_symmetric_union) {
      overlap_candidate = generic_options->orientation ==
          KLS_ORIENTATION_TRANSPOSE
        ? transpose : normal;
      if (overlap_candidate != NULL) {
        overlap_candidate->generic_nd_overlap_target = &overlapped_nd;
        /* The minimum-degree tournament settles BTF before launching this
           worker through the handoff above.  An explicit diagnostic can
           retain the old eager launch, but AUTO now spends resources on only
           the realized coordinate frame. */
        if (getenv("KLS_ENABLE_GENERIC_EARLY_COORDINATE_FREE_ND") != NULL) {
          (void)kls_generic_nd_trial_start(
            &overlapped_nd, overlap_candidate, generic_options,
            generic_options->use_btf ? 1 : 0);
        }
      }
    }
    const int status =
      select_candidate_inner(normal, transpose, generic_options, chosen_out);
    if (overlap_candidate != NULL) {
      overlap_candidate->generic_nd_overlap_target = NULL;
    }
    if (status == KLS_OK && chosen_out != NULL && *chosen_out != NULL) {
      const int overlap_started =
        overlapped_nd.active || overlapped_nd.symbolic != NULL;
      const int overlap_used_btf =
        overlap_started && overlapped_nd.options.use_btf;
      /* Orientation is settled: pay for at most one tuned NodeND candidate
         first and compare it to the selected minimum-degree representation.
         On an unsymmetric, high-work lifecycle the other orientation may
         expose a materially better separator tree even though minimum degree
         preferred the first frame.  Let the same bounded ND gate inspect that
         already-analyzed candidate too, then repeat only the cheap lifecycle
         score comparison.  A bounded symmetric-union certificate admits
         predicted construction without turning union expansion into a
         recurring cost; the recorded minimum-degree rollback still polices
         numeric expansion before a retained factor can be published. */
      kls_maybe_promote_selected_generic_nd(
        *chosen_out, generic_options, bounded_symmetric_union,
        *chosen_out == overlap_candidate && overlap_started &&
          !overlap_used_btf
          ? &overlapped_nd : NULL,
        *chosen_out == overlap_candidate && overlap_used_btf
          ? &overlapped_nd : NULL);
      kls_generic_nd_trial_clear(&overlapped_nd);
      if (!bounded_symmetric_union &&
          generic_options->orientation == KLS_ORIENTATION_AUTO &&
          normal != NULL && transpose != NULL) {
        kls_pattern_candidate *other =
          *chosen_out == normal ? transpose : normal;
        const double chosen_lifecycle =
          kls_candidate_lifecycle_score(*chosen_out, generic_options);
        const double other_lifecycle =
          kls_candidate_lifecycle_score(other, generic_options);
        const int inspect_other_nd =
          !(*chosen_out)->generic_nd_portfolio_selected ||
          ((*chosen_out)->generic_nd_lifecycle_near_tie &&
           getenv("KLS_DISABLE_GENERIC_NEAR_TIE_CROSS_ORIENTATION_ND") ==
             NULL) ||
          (isfinite(chosen_lifecycle) && chosen_lifecycle > 0.0 &&
           isfinite(other_lifecycle) && other_lifecycle > 0.0 &&
           other_lifecycle < 0.98 * chosen_lifecycle);
        /* Once the lifecycle-preferred orientation has produced a decisive
           ND representation, a second NodeND pass has no supporting signal
           unless the other orientation was already decisively cheaper.  A
           near-tie ND admission is different: its win over minimum degree is
           small enough that the other frame's independently realized tree
           can reverse the lifecycle verdict, so complete that portfolio too.
           A rejected first proposal still gives the other side its full
           trial.  This is a portfolio-economics decision over realized
           symbolic costs, not an orientation, dimension, or matrix-family
           rule. */
        if (inspect_other_nd) {
          kls_maybe_promote_selected_generic_nd(
            other, generic_options, bounded_symmetric_union,
            NULL, NULL);
        } else if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
          fprintf(stderr,
                  "KLS cross-orientation ND: skip unsupported second pass "
                  "(chosen %.4e other %.4e)\n",
                  chosen_lifecycle, other_lifecycle);
        }
        if (normal->symbolic != NULL && transpose->symbolic != NULL) {
          const double normal_lifecycle =
            kls_candidate_lifecycle_score(normal, generic_options);
          const double transpose_lifecycle =
            kls_candidate_lifecycle_score(transpose, generic_options);
          *chosen_out = transpose_lifecycle < normal_lifecycle
            ? transpose : normal;
          if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
            fprintf(stderr,
                    "KLS cross-orientation ND lifecycle: normal %.4e "
                    "transpose %.4e -> %s\n",
                    normal_lifecycle, transpose_lifecycle,
                    (*chosen_out)->orientation == KLS_ORIENTATION_NORMAL
                      ? "normal" : "transpose");
          }
        }
      }
    } else {
      kls_generic_nd_trial_clear(&overlapped_nd);
    }
    return status;
  }

  /* The pre-static large-spral class replaces the analyze-time
     symbolic wholesale at factor time (matched pattern, its own METIS
     analysis), so its analyze never needs to pay nested dissection:
     engage the (previously dormant) defer-ND mode for the whole
     analyze.  pre2: two concurrent 3.7s METIS analyses -> AMD-class
     scoring only.  If the pre-static trial rejects, the factor-time
     METIS race/promotion machinery recovers nested dissection under
     timed acceptance as on any other auto-path matrix. */
  kls_pattern_candidate *ref = normal != NULL ? normal : transpose;
  int defer = 0;
  if (ref != NULL && options != NULL && options->static_pivoting &&
      options->ordering == KLS_ORDERING_AUTO &&
      ref->n > 150000u && ref->n <= 750000u && ref->nnz <= 8000000u &&
      !is_large_sparse_diagonal_low_degree_pattern(ref->n, ref->col_ptr,
                                                   ref->row_idx) &&
      !kls_large_low_degree_fragmented_pattern(ref->n, ref->col_ptr) &&
      !is_large_very_low_degree_full_diagonal_pattern(
        ref->n, ref->col_ptr, ref->row_idx) &&
      getenv("KLS_DISABLE_ANA_ND_DEFER") == NULL) {
    defer = 1;
  }
  const int saved = kls_analyze_defer_nd;
  if (defer) {
    kls_analyze_defer_nd = 1;
  }
  /* the pre-static adoption forces the NORMAL orientation for this
     class, overriding whatever the two-orientation comparison picks:
     skip the transpose side entirely (a class member that fails the
     pre-static trial falls back to the normal orientation - correct,
     at worst slower, and policed by the same measured acceptance) */
  const int status = select_candidate_inner(
    normal, defer && normal != NULL ? NULL : transpose, options,
    chosen_out);
  kls_analyze_defer_nd = saved;
  return status;
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
  kls_separator_analysis_move(&solver->separator, &candidate->separator);
  kls_invalidate_factor_etree_stats(solver);
  solver->stats.selected_ordering = candidate->selected_ordering;
  solver->stats.selected_orientation = candidate->orientation;
  solver->generic_nd_portfolio_selected =
    candidate->generic_nd_portfolio_selected;
  solver->generic_nd_bounded_symmetric_union =
    candidate->generic_nd_bounded_symmetric_union;
  solver->generic_nd_lifecycle_near_tie =
    candidate->generic_nd_lifecycle_near_tie;
  solver->generic_amf3_span_variant_selected =
    candidate->generic_amf3_span_variant_selected &&
    candidate->selected_ordering == KLS_ORDERING_AMF3;
  solver->generic_btf_value_symbolic =
    candidate->generic_btf_value_symbolic;
  solver->generic_btf_value_common =
    candidate->generic_btf_value_common;
  solver->generic_btf_value_ordering =
    candidate->generic_btf_value_ordering;
  kls_separator_analysis_move(
    &solver->generic_btf_value_separator,
    &candidate->generic_btf_value_separator);
  candidate->generic_btf_value_symbolic = NULL;
  solver->generic_nd_fallback_ordering =
    candidate->generic_nd_fallback_ordering;
  solver->generic_nd_fallback_use_btf =
    candidate->generic_nd_fallback_use_btf;
  solver->generic_nd_fallback_fill = candidate->generic_nd_fallback_fill;
  solver->generic_nd_fallback_flops = candidate->generic_nd_fallback_flops;
  free(solver->block_order_perm);
  solver->block_order_perm = candidate->block_order_perm;
  candidate->block_order_perm = NULL;
  {
    candidate->col_ptr = NULL;
    candidate->row_idx = NULL;
    candidate->input_to_csc = NULL;
    candidate->values = NULL;
    candidate->symbolic = NULL;
    return;
  }
}

/* A pre-numeric ND choice is provisional until the ordinary predicted-factor
   construction verifies that its realized pivot frame stays near the
   symbolic fill.  If that construction rejects, rebuild the exact
   minimum-degree ordering recorded by the portfolio and let the normal KLU
   factor path continue.  This turns numeric instability into a bounded
   rollback instead of committing to a potentially enormous ND factor. */
static int kls_restore_generic_nd_fallback_symbolic(kls_solver *solver,
                                                     double *elapsed) {
  if (solver == NULL || elapsed == NULL || solver->symbolic == NULL ||
      !solver->generic_nd_portfolio_selected) {
    return 0;
  }

  kls_options fallback_options = solver->options;
  fallback_options.ordering = solver->generic_nd_fallback_ordering;
  fallback_options.use_btf = solver->generic_nd_fallback_use_btf;
  trilinos_klu_l_symbolic *fallback_symbolic = NULL;
  trilinos_klu_l_common fallback_common;
  kls_separator_analysis fallback_separator;
  memset(&fallback_separator, 0, sizeof(fallback_separator));
  const double start = kls_now_seconds();
  const int status = analyze_with_ordering(
    solver->n, solver->col_ptr, solver->row_idx, &fallback_options,
    solver->generic_nd_fallback_ordering, &fallback_symbolic,
    &fallback_common, &fallback_separator);
  *elapsed += kls_now_seconds() - start;
  if (status != KLS_OK || fallback_symbolic == NULL) {
    kls_separator_analysis_clear(&fallback_separator);
    return 0;
  }

  /* Preserve the value-dependent controls settled at factor entry; the new
     Common otherwise carries the same caller options plus fresh analysis
     statistics for the restored symbolic. */
  fallback_common.scale = solver->common.scale;
  fallback_common.tol = solver->common.tol;
  fallback_common.halt_if_singular = solver->common.halt_if_singular;
  fallback_common.kls_static_perturb = solver->common.kls_static_perturb;
  fallback_common.kls_zero_pivot_replacement =
    solver->common.kls_zero_pivot_replacement;

  trilinos_klu_l_symbolic *nd_symbolic = solver->symbolic;
  trilinos_klu_l_common nd_common = solver->common;
  if (solver->metis_race != NULL &&
      solver->metis_race->scale_symbolic == nd_symbolic) {
    /* A scale-only race factors the incumbent symbolic after factor entry
       publishes the values.  Numeric validation is allowed to replace that
       provisional ordering, but its worker must no longer dereference the
       old P/Q/R arrays when they are freed below.  Results computed for the
       rejected ordering are unusable in any case, so join and discard them
       before transferring ownership. */
    kls_metis_race_abandon(solver);
  }
  solver->symbolic = fallback_symbolic;
  solver->common = fallback_common;
  kls_separator_analysis_clear(&solver->separator);
  kls_separator_analysis_move(&solver->separator, &fallback_separator);
  kls_invalidate_factor_etree_stats(solver);
  solver->stats.selected_ordering = solver->generic_nd_fallback_ordering;
  solver->stats.nblocks = (int64_t)solver->symbolic->nblocks;
  solver->stats.max_block = (int64_t)solver->symbolic->maxblock;
  solver->stats.structural_rank = (int64_t)solver->symbolic->structural_rank;
  solver->stats.estimated_flops = solver->symbolic->est_flops;
  solver->generic_nd_portfolio_selected = 0;
  solver->generic_nd_lifecycle_near_tie = 0;
  solver->auto_metis_checked = 1;
  trilinos_klu_l_free_symbolic(&nd_symbolic, &nd_common);
  if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL ||
      getenv("KLS_TRACE_PREDICTED") != NULL) {
    fprintf(stderr,
            "KLS generic ND portfolio: predicted fill rejected; "
            "restored %s\n",
            kls_ordering_name(solver->stats.selected_ordering));
  }
  return 1;
}

/* A lower-span AMF3 ordering is retained only after the first value set
   proves that its static numeric is usable.  Rebuild the ordinary AMF3
   symbolic on rejection so an approximate structural tie cannot expose the
   recurring lifecycle to value-dependent pivot fill. */
static int kls_restore_generic_amf3_span_symbolic(kls_solver *solver,
                                                   double *elapsed) {
  if (solver == NULL || elapsed == NULL || solver->symbolic == NULL ||
      !solver->generic_amf3_span_variant_selected) {
    return 0;
  }

  kls_options fallback_options = solver->options;
  fallback_options.ordering = KLS_ORDERING_AMF3;
  fallback_options.use_btf = solver->symbolic->do_btf ? 1 : 0;
  trilinos_klu_l_symbolic *fallback_symbolic = NULL;
  trilinos_klu_l_common fallback_common;
  kls_separator_analysis fallback_separator;
  memset(&fallback_separator, 0, sizeof(fallback_separator));
  const int saved_power_override_active = kls_amf3_power_override_active;
  const double saved_power_override = kls_amf3_power_override;
  kls_amf3_power_override_active = 1;
  kls_amf3_power_override = 0.445;
  const double start = kls_now_seconds();
  const int status = analyze_with_ordering(
    solver->n, solver->col_ptr, solver->row_idx, &fallback_options,
    KLS_ORDERING_AMF3, &fallback_symbolic, &fallback_common,
    &fallback_separator);
  *elapsed += kls_now_seconds() - start;
  kls_amf3_power_override_active = saved_power_override_active;
  kls_amf3_power_override = saved_power_override;
  if (status != KLS_OK || fallback_symbolic == NULL) {
    kls_separator_analysis_clear(&fallback_separator);
    return 0;
  }

  fallback_common.scale = solver->common.scale;
  fallback_common.tol = solver->common.tol;
  fallback_common.halt_if_singular = solver->common.halt_if_singular;
  fallback_common.kls_static_perturb = solver->common.kls_static_perturb;
  fallback_common.kls_zero_pivot_replacement =
    solver->common.kls_zero_pivot_replacement;

  trilinos_klu_l_symbolic *span_symbolic = solver->symbolic;
  trilinos_klu_l_common span_common = solver->common;
  if (solver->metis_race != NULL &&
      solver->metis_race->scale_symbolic == span_symbolic) {
    /* See the ND rollback above: a live scale race owns a read lease on this
       provisional symbolic.  End that lease before freeing the rejected
       span variant. */
    kls_metis_race_abandon(solver);
  }
  solver->symbolic = fallback_symbolic;
  solver->common = fallback_common;
  kls_separator_analysis_clear(&solver->separator);
  kls_separator_analysis_move(&solver->separator, &fallback_separator);
  kls_invalidate_factor_etree_stats(solver);
  solver->stats.selected_ordering = KLS_ORDERING_AMF3;
  solver->stats.nblocks = (int64_t)solver->symbolic->nblocks;
  solver->stats.max_block = (int64_t)solver->symbolic->maxblock;
  solver->stats.structural_rank =
    (int64_t)solver->symbolic->structural_rank;
  solver->stats.estimated_flops = solver->symbolic->est_flops;
  solver->generic_amf3_span_variant_selected = 0;
  trilinos_klu_l_free_symbolic(&span_symbolic, &span_common);
  if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL ||
      getenv("KLS_TRACE_PREDICTED") != NULL) {
    fprintf(stderr,
            "KLS AMF3 span portfolio: numeric rejected; restored "
            "ordinary power\n");
  }
  return 1;
}

typedef struct kls_symbolic_pivot_profile {
  UF_long missing;
  UF_long weak_1e3;
  UF_long weak_1e2;
} kls_symbolic_pivot_profile;

/* Bound value-admission work independently of matrix size.  Evenly spaced
   elimination positions cover the whole symbolic while keeping a rejected
   alternative cheaper than the numeric factor trial it avoids.  Acceptance
   still depends on the complete candidate factor and backward residual. */
static kls_symbolic_pivot_profile kls_sample_symbolic_pivots(
  UF_long n,
  const UF_long *col_ptr,
  const UF_long *row_idx,
  const double *values,
  const trilinos_klu_l_symbolic *symbolic) {
  kls_symbolic_pivot_profile profile;
  memset(&profile, 0, sizeof(profile));
  if (n == 0u || col_ptr == NULL || row_idx == NULL || values == NULL ||
      symbolic == NULL || symbolic->P == NULL || symbolic->Q == NULL) {
    profile.missing = n;
    return profile;
  }
  const UF_long sample_count = n < 1024u ? n : 1024u;
  const double stride = (double)n / (double)sample_count;
  for (UF_long sample = 0u; sample < sample_count; ++sample) {
    const UF_long k = sample_count == n
      ? sample
      : (UF_long)(((double)sample + 0.5) * stride);
    const UF_long col = symbolic->Q[k];
    const UF_long row = symbolic->P[k];
    if (col >= n || row >= n) {
      profile.missing++;
      continue;
    }
    double column_max = 0.0;
    double pivot_abs = 0.0;
    int found = 0;
    for (UF_long p = col_ptr[col]; p < col_ptr[col + 1u]; ++p) {
      const double av = fabs(values[p]);
      column_max = column_max < av ? av : column_max;
      if (row_idx[p] == row) {
        pivot_abs = av;
        found = 1;
      }
    }
    if (!found) {
      profile.missing++;
      continue;
    }
    const double ratio = column_max > 0.0 ? pivot_abs / column_max : 1.0;
    profile.weak_1e3 += ratio < 1.0e-3;
    profile.weak_1e2 += ratio < 1.0e-2;
  }
  return profile;
}

static void kls_discard_generic_btf_value_alternative(kls_solver *solver) {
  if (solver == NULL) {
    return;
  }
  if (solver->generic_btf_value_symbolic != NULL) {
    trilinos_klu_l_free_symbolic(
      &solver->generic_btf_value_symbolic,
      &solver->generic_btf_value_common);
  }
  kls_separator_analysis_clear(&solver->generic_btf_value_separator);
  solver->generic_btf_value_ordering = KLS_ORDERING_AUTO;
}

static int kls_generic_hoisted_snode_worker_candidate(
  const kls_solver *solver);
static int kls_generic_hoisted_snode_worker_capable(
  const kls_solver *solver);
static void kls_lean_row_refactor_ensure_snode_runs(kls_solver *solver);
static int kls_build_row_refactor_pattern(kls_solver *solver,
                                          int compact_only);
typedef struct kls_lean_prewarm_job {
  kls_solver *solver;
  int thread_count;
  int build_pattern;
} kls_lean_prewarm_job;
static void *kls_lean_prewarm_main(void *arg);

/* The structural tournament has already paid for both AMD representations.
   If it retained a bounded BTF alternative, use the first value frame to
   decide which symbolic to factor first.  Static pivot strength or bounded
   recurring-work economics is only an admission signal: the candidate must
   then pass a true backward-residual solve.  A non-Pareto candidate admitted
   specifically for row-executor economics must additionally build that exact
   retained-row executor before it can replace the ordinary no-BTF route. */
static int maybe_factor_generic_btf_value_alternative(
  kls_solver *solver,
  const double *numeric_values,
  double *elapsed) {
  if (solver == NULL || numeric_values == NULL || elapsed == NULL ||
      solver->generic_btf_value_symbolic == NULL) {
    return 0;
  }
  trilinos_klu_l_symbolic *candidate_symbolic =
    solver->generic_btf_value_symbolic;
  const int capable =
    getenv("KLS_DISABLE_GENERIC_BTF_VALUE_SELECTION") == NULL &&
        solver->options.ordering == KLS_ORDERING_AUTO &&
    solver->options.scale == KLS_SCALE_AUTO &&
    solver->options.use_btf &&
    kls_repeated_update_workload(&solver->options) &&
    !solver->solve_recovery_active && solver->numeric == NULL &&
    solver->input_format == KLS_INPUT_CSC &&
    solver->row_perm == NULL && solver->user_col_perm == NULL &&
    solver->row_scale == NULL && solver->col_scale == NULL &&
    solver->block_order_perm == NULL && solver->symbolic != NULL &&
    !solver->symbolic->do_btf &&
    candidate_symbolic->do_btf && candidate_symbolic->nblocks > 1u &&
    candidate_symbolic->maxblock < solver->n;
  if (!capable) {
    kls_discard_generic_btf_value_alternative(solver);
    return 0;
  }
  const double decision_start = kls_now_seconds();

  const kls_symbolic_pivot_profile incumbent_profile =
    kls_sample_symbolic_pivots(
      solver->n, solver->col_ptr, solver->row_idx, numeric_values,
      solver->symbolic);
  const kls_symbolic_pivot_profile candidate_profile =
    kls_sample_symbolic_pivots(
      solver->n, solver->col_ptr, solver->row_idx, numeric_values,
      candidate_symbolic);
  const int decisive_value_strength =
    candidate_profile.missing == 0u &&
    incumbent_profile.missing >= 8u &&
    incumbent_profile.weak_1e3 >= 8u &&
    (double)candidate_profile.weak_1e3 <=
      0.75 * (double)incumbent_profile.weak_1e3 &&
    (double)candidate_profile.weak_1e2 <=
      0.75 * (double)incumbent_profile.weak_1e2;
  const double candidate_symbolic_fill =
    candidate_symbolic->lnz + candidate_symbolic->unz;
  const double incumbent_symbolic_fill = symbolic_score(solver->symbolic);
  const double recurring_frames =
    (double)solver->options.expected_refactorizations;
  const double solve_frames = (double)solver->options.expected_solves;
  const double symbolic_saved_fill =
    fmax(0.0, incumbent_symbolic_fill - candidate_symbolic_fill);
  const double symbolic_saved_flops =
    fmax(0.0,
         solver->symbolic->est_flops - candidate_symbolic->est_flops);
  const double projected_symbolic_savings =
    recurring_frames * (symbolic_saved_fill + symbolic_saved_flops) +
    solve_frames * symbolic_saved_fill;
  const double symbolic_trial_work =
    candidate_symbolic_fill + candidate_symbolic->est_flops;
  const int repeated_symbolic_pareto_proposal =
    incumbent_symbolic_fill > 0.0 &&
    isfinite(incumbent_symbolic_fill) &&
    candidate_symbolic_fill > 0.0 &&
    isfinite(candidate_symbolic_fill) &&
    candidate_symbolic_fill <= incumbent_symbolic_fill &&
    solver->symbolic->est_flops > 0.0 &&
    candidate_symbolic->est_flops > 0.0 &&
    candidate_symbolic->est_flops <= solver->symbolic->est_flops &&
    isfinite(projected_symbolic_savings) &&
    isfinite(symbolic_trial_work) && symbolic_trial_work > 0.0 &&
    projected_symbolic_savings >= 2.0 * symbolic_trial_work;
  /* A materially smaller retained factor can win a solve-heavy lifecycle
     even when the arithmetic estimate is a close near-tie.  Admit a BTF
     candidate that saves at least ten percent of symbolic storage, adds at
     most five percent estimated arithmetic, and whose cumulative factor
     traffic saving repays sixteen complete candidate patterns.  The actual
     numeric below must satisfy the same bounds and the backward-residual
     probe before it can be published. */
  const int repeated_symbolic_storage_pareto_proposal =
    kls_repeated_update_workload(&solver->options) &&
    recurring_frames + solve_frames >= 64.0 &&
    incumbent_symbolic_fill > 0.0 && isfinite(incumbent_symbolic_fill) &&
    candidate_symbolic_fill > 0.0 && isfinite(candidate_symbolic_fill) &&
    candidate_symbolic_fill <= 0.90 * incumbent_symbolic_fill &&
    solver->symbolic->est_flops > 0.0 &&
    candidate_symbolic->est_flops > 0.0 &&
    candidate_symbolic->est_flops <= 1.05 * solver->symbolic->est_flops &&
    symbolic_saved_fill * (recurring_frames + solve_frames) >=
      16.0 * candidate_symbolic_fill &&
    getenv("KLS_DISABLE_GENERIC_BTF_STORAGE_PARETO") == NULL;
  const int hoisted_row_executor_proposal =
    getenv("KLS_DISABLE_GENERIC_BTF_HOISTED_ROW_RECOVERY") == NULL &&
    solver->n >= 512u && solver->n <= 131072u &&
    candidate_symbolic_fill > 0.0 &&
    candidate_symbolic_fill <= 1000000.0 &&
    candidate_symbolic->est_flops > 0.0 &&
    candidate_symbolic->est_flops >=
      32.0 * candidate_symbolic_fill;
  if (!decisive_value_strength && !repeated_symbolic_pareto_proposal &&
      !repeated_symbolic_storage_pareto_proposal &&
      !hoisted_row_executor_proposal) {
    if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
      fprintf(stderr,
              "KLS generic BTF values: missing=%ld/%ld "
              "weak1e3=%ld/%ld weak1e2=%ld/%ld "
              "symbolic_savings=%.4e/%.4e -> reject\n",
              (long)candidate_profile.missing,
              (long)incumbent_profile.missing,
              (long)candidate_profile.weak_1e3,
              (long)incumbent_profile.weak_1e3,
              (long)candidate_profile.weak_1e2,
              (long)incumbent_profile.weak_1e2,
              projected_symbolic_savings, symbolic_trial_work);
    }
    kls_discard_generic_btf_value_alternative(solver);
    *elapsed += kls_now_seconds() - decision_start;
    return 0;
  }

  trilinos_klu_l_common candidate_common =
    solver->generic_btf_value_common;
  const int recovery_scale = (int)solver->common.scale;
  /* A representation admitted by recurring symbolic Pareto economics uses
     AUTO's value-selected scale.  Its purpose is lower recurring structure,
     not an unscaled-conditioning gamble.  Pivot-strength and row-executor
     recoveries retain the historical unscaled trial and its residual guard. */
  candidate_common.scale =
    (repeated_symbolic_pareto_proposal ||
     repeated_symbolic_storage_pareto_proposal)
    ? recovery_scale : -1;
  candidate_common.tol = solver->common.tol;
  candidate_common.halt_if_singular = solver->common.halt_if_singular;
  candidate_common.kls_static_perturb = solver->common.kls_static_perturb;
  candidate_common.kls_zero_pivot_replacement =
    solver->common.kls_zero_pivot_replacement;
  pthread_t prewarm_thread;
  kls_lean_prewarm_job prewarm_job;
  int prewarm_active = 0;
  if (repeated_symbolic_pareto_proposal && solver->options.threads > 1 &&
      solver->egraph_pool == NULL &&
      getenv("KLS_DISABLE_GENERIC_BTF_LEAN_PREWARM") == NULL) {
    /* The admitted symbolic has already proved enough recurring savings to
       repay its numeric trial.  Build the representation-independent worker
       crew beside that serial KLU factor so a successful candidate does not
       move thread startup onto its first public update. */
    prewarm_job.solver = solver;
    prewarm_job.thread_count = solver->options.threads;
    prewarm_job.build_pattern = 0;
    prewarm_active = pthread_create(
      &prewarm_thread, NULL, kls_lean_prewarm_main, &prewarm_job) == 0;
  }
  const double trial_start = kls_now_seconds();
  trilinos_klu_l_numeric *candidate_numeric = trilinos_klu_l_factor(
    solver->col_ptr, solver->row_idx,
    (double *)(uintptr_t)numeric_values, candidate_symbolic,
    &candidate_common);
  const double factor_seconds = kls_now_seconds() - trial_start;
  if (prewarm_active) {
    pthread_join(prewarm_thread, NULL);
  }
  int candidate_ok = candidate_numeric != NULL &&
    candidate_common.status >= TRILINOS_KLU_OK &&
    candidate_common.status != TRILINOS_KLU_SINGULAR;
  if (candidate_ok) {
    (void)trilinos_klu_l_flops(
      candidate_symbolic, candidate_numeric, &candidate_common);
    (void)trilinos_klu_l_rcond(
      candidate_symbolic, candidate_numeric, &candidate_common);
  }
  const double incumbent_estimated_fill = symbolic_score(solver->symbolic);
  const double incumbent_estimated_flops = solver->symbolic->est_flops;
  const double candidate_fill = candidate_ok
    ? (double)candidate_numeric->lnz + (double)candidate_numeric->unz
    : DBL_MAX;
  const UF_long stable_pivot_cap = solver->n / 10000u + 16u;
  const int decisive_numeric_economics =
    candidate_fill <= 0.90 * incumbent_estimated_fill &&
    candidate_common.flops > 0.0 &&
    candidate_common.flops <= 0.75 * incumbent_estimated_flops;
  const int repeated_numeric_pareto =
    kls_repeated_update_workload(&solver->options) &&
    candidate_fill <= 0.99 * incumbent_estimated_fill &&
    candidate_common.flops > 0.0 &&
    candidate_common.flops <= 0.99 * incumbent_estimated_flops;
  const int repeated_numeric_storage_pareto =
    repeated_symbolic_storage_pareto_proposal &&
    candidate_fill <= 0.90 * incumbent_estimated_fill &&
    candidate_common.flops > 0.0 &&
    candidate_common.flops <= 1.05 * incumbent_estimated_flops &&
    candidate_common.noffdiag <= stable_pivot_cap;
  const int hoisted_row_numeric_economics =
    hoisted_row_executor_proposal && candidate_fill > 0.0 &&
    candidate_fill <= 1000000.0 &&
    candidate_fill <= 1.05 * incumbent_estimated_fill &&
    candidate_common.flops >= 32.0 * candidate_fill &&
    /* Executor capability does not erase numeric work.  A candidate whose
       row schedule is merely buildable can still lose every retained update
       when its factor performs more arithmetic than the incumbent EGraph.
       Keep this representation change work-Pareto; the exact row-capability
       and backward-residual gates below continue to arbitrate candidates
       whose recurring work is no larger. */
    candidate_common.flops <= incumbent_estimated_flops &&
    candidate_common.noffdiag <= stable_pivot_cap;
  candidate_ok = candidate_ok && incumbent_estimated_fill > 0.0 &&
    isfinite(incumbent_estimated_fill) &&
    incumbent_estimated_flops > 0.0 &&
    (decisive_numeric_economics || repeated_numeric_pareto ||
     repeated_numeric_storage_pareto ||
     hoisted_row_numeric_economics) &&
    isfinite(candidate_common.rcond) && candidate_common.rcond > 0.0;

  int residual_ok = 0;
  if (candidate_ok) {
    trilinos_klu_l_symbolic *incumbent_symbolic = solver->symbolic;
    trilinos_klu_l_numeric *incumbent_numeric = solver->numeric;
    const trilinos_klu_l_common incumbent_common = solver->common;
    solver->symbolic = candidate_symbolic;
    solver->numeric = candidate_numeric;
    solver->common = candidate_common;
    residual_ok = kls_direct_klu_numeric_residual_probe(
      solver, numeric_values);
    candidate_common = solver->common;
    solver->symbolic = incumbent_symbolic;
    solver->numeric = incumbent_numeric;
    solver->common = incumbent_common;
  }
  candidate_ok = candidate_ok && residual_ok;
  if (getenv("KLS_TRACE_FACTOR_PHASES") != NULL) {
    fprintf(stderr,
            "KLS generic BTF values: fill %.4e/%.4e flops %.4e/%.4e "
            "pivots=%ld rcond=%.3e residual=%d -> %s\n",
            candidate_fill, incumbent_estimated_fill,
            candidate_common.flops, incumbent_estimated_flops,
            (long)candidate_common.noffdiag, candidate_common.rcond,
            residual_ok, candidate_ok ? "adopt" : "reject");
  }
  if (!candidate_ok) {
    if (candidate_numeric != NULL) {
      trilinos_klu_l_free_numeric(&candidate_numeric, &candidate_common);
    }
    kls_discard_generic_btf_value_alternative(solver);
    *elapsed += kls_now_seconds() - decision_start;
    return 0;
  }

  trilinos_klu_l_symbolic *old_symbolic = solver->symbolic;
  trilinos_klu_l_common old_common = solver->common;
  const kls_stats old_stats = solver->stats;
  const int old_generic_nd_portfolio_selected =
    solver->generic_nd_portfolio_selected;
  const int old_generic_nd_numeric_validated =
    solver->generic_nd_numeric_validated;
  const int old_generic_amf3_span_variant_selected =
    solver->generic_amf3_span_variant_selected;
  const int old_auto_scale_checked = solver->auto_scale_checked;
  const int old_auto_scale_value_certified =
    solver->auto_scale_value_certified;
  const int old_auto_scale_unscaled_trial_certified =
    solver->auto_scale_unscaled_trial_certified;
  const int old_factor_preps_deferred = solver->factor_preps_deferred;
  kls_separator_analysis old_separator;
  memset(&old_separator, 0, sizeof(old_separator));
  kls_separator_analysis_move(&old_separator, &solver->separator);
  solver->symbolic = candidate_symbolic;
  solver->numeric = candidate_numeric;
  solver->common = candidate_common;
  solver->generic_btf_value_symbolic = NULL;
  kls_separator_analysis_move(
    &solver->separator, &solver->generic_btf_value_separator);
  solver->stats.selected_ordering = solver->generic_btf_value_ordering;
  solver->generic_btf_value_ordering = KLS_ORDERING_AUTO;
  solver->generic_nd_portfolio_selected = 0;
  solver->generic_nd_numeric_validated = 0;
  solver->generic_amf3_span_variant_selected = 0;
  solver->auto_scale_checked = 1;
  solver->auto_scale_value_certified = 1;
  solver->auto_scale_unscaled_trial_certified = 0;
  kls_invalidate_factor_etree_stats(solver);
  kls_numeric_replaced_invalidate(solver);
  solver->generic_btf_unscaled_recovery_scale =
    candidate_common.scale <= 0 ? recovery_scale : 0;
  /* The first factor's true residual certifies its absolute conditioning
     regime.  Future value frames may remain in that regime; recover only
     after the cheap diagonal estimate deteriorates by a full factor of
     eight, rather than imposing a matrix-independent absolute cutoff on an
     already accurate ill-conditioned solve. */
  solver->generic_btf_unscaled_rcond_floor =
    candidate_common.scale <= 0
      ? fmax(DBL_MIN, 0.125 * candidate_common.rcond) : 0.0;
  solver->numeric_full_factor_seconds = factor_seconds;
  solver->stats.nblocks = (int64_t)solver->symbolic->nblocks;
  solver->stats.max_block = (int64_t)solver->symbolic->maxblock;
  solver->stats.structural_rank =
    (int64_t)solver->symbolic->structural_rank;
  solver->stats.estimated_flops = solver->symbolic->est_flops;
  if (hoisted_row_numeric_economics && !repeated_numeric_pareto) {
    const int row_ready =
      kls_generic_hoisted_snode_worker_candidate(solver) &&
      kls_build_row_refactor_pattern(solver, 1);
    if (row_ready) {
      kls_lean_row_refactor_ensure_snode_runs(solver);
    }
    if (!row_ready || !kls_generic_hoisted_snode_worker_capable(solver)) {
      /* The fill/work screen is only an admission signal.  If the realized
         factor cannot build the executor whose recurring traffic justified
         the extra numeric, restore the ordinary symbolic before publishing
         any result.  No input dimensions or family identity participate in
         this final capability certificate. */
      trilinos_klu_l_symbolic *rejected_symbolic = solver->symbolic;
      trilinos_klu_l_numeric *rejected_numeric = solver->numeric;
      trilinos_klu_l_common rejected_common = solver->common;
      solver->symbolic = old_symbolic;
      solver->numeric = NULL;
      solver->common = old_common;
      solver->stats = old_stats;
      solver->generic_nd_portfolio_selected =
        old_generic_nd_portfolio_selected;
      solver->generic_nd_numeric_validated =
        old_generic_nd_numeric_validated;
      solver->generic_amf3_span_variant_selected =
        old_generic_amf3_span_variant_selected;
      solver->auto_scale_checked = old_auto_scale_checked;
      solver->auto_scale_value_certified =
        old_auto_scale_value_certified;
      solver->auto_scale_unscaled_trial_certified =
        old_auto_scale_unscaled_trial_certified;
      kls_separator_analysis_clear(&solver->separator);
      kls_separator_analysis_move(&solver->separator, &old_separator);
      kls_invalidate_factor_etree_stats(solver);
      kls_numeric_replaced_invalidate(solver);
      solver->factor_preps_deferred = old_factor_preps_deferred;
      trilinos_klu_l_free_numeric(&rejected_numeric, &rejected_common);
      trilinos_klu_l_free_symbolic(&rejected_symbolic, &rejected_common);
      *elapsed += kls_now_seconds() - decision_start;
      return 0;
    }
    solver->lean_choice = 1;
    /* The replacement factor proves row-executor capability and recurring
       work economics, not the solve cost of its row representation.  Keep
       the choice provisional until the ordinary 4+4 lifecycle audit. */
    solver->lean_reaudit_state = 0;
  }
  kls_metis_race_abandon(solver);
  kls_separator_analysis_clear(&old_separator);
  trilinos_klu_l_free_symbolic(&old_symbolic, &old_common);
  *elapsed += kls_now_seconds() - decision_start;
  return 1;
}

static double kls_row_refactor_total_group_work(const kls_solver *solver);
static int kls_auto_row_refactor_cost_allows(const kls_solver *solver);
static int kls_auto_row_refactor_should_run(const kls_solver *solver);

static void kls_compute_parallel_model(const kls_solver *solver,
                                       int prefer_numeric,
                                       double *r1_out,
                                       double *r2_out,
                                       int *recommends_parallel_out) {
  if (r1_out != NULL) {
    *r1_out = 0.0;
  }
  if (r2_out != NULL) {
    *r2_out = 0.0;
  }
  if (recommends_parallel_out != NULL) {
    *recommends_parallel_out = 0;
  }
  if (solver == NULL) {
    return;
  }

  double lu_nnz = 0.0;
  if (prefer_numeric && solver->numeric != NULL) {
    lu_nnz = (double)solver->numeric->lnz + (double)solver->numeric->unz;
  }
  if (lu_nnz <= 0.0 && solver->symbolic != NULL) {
    lu_nnz = (double)solver->symbolic->lnz + (double)solver->symbolic->unz;
  }

  double r1 = 0.0;
  double r2 = 0.0;
  if (lu_nnz > 0.0) {
    if (solver->nnz > 0u) {
      r1 = lu_nnz / (double)solver->nnz;
    }
    if (solver->n > 0u) {
      r2 = lu_nnz / (double)solver->n;
    }
  }

  const int recommends_parallel =
    (r1 >= KLS_NICSLU_PARALLEL_R1_THRESHOLD ||
     r2 >= KLS_NICSLU_PARALLEL_R2_THRESHOLD);
  if (r1_out != NULL) {
    *r1_out = r1;
  }
  if (r2_out != NULL) {
    *r2_out = r2;
  }
  if (recommends_parallel_out != NULL) {
    *recommends_parallel_out = recommends_parallel;
  }
}

static void kls_compute_nicslu_task_flow_model(
  const kls_solver *solver,
  int thread_count,
  int64_t *threads_out,
  int64_t *dependencies_out,
  double *work_out,
  double *finish_time_out,
  double *speedup_out,
  int *recommends_parallel_out) {
  if (threads_out != NULL) {
    *threads_out = 0;
  }
  if (dependencies_out != NULL) {
    *dependencies_out = 0;
  }
  if (work_out != NULL) {
    *work_out = 0.0;
  }
  if (finish_time_out != NULL) {
    *finish_time_out = 0.0;
  }
  if (speedup_out != NULL) {
    *speedup_out = 0.0;
  }
  if (recommends_parallel_out != NULL) {
    *recommends_parallel_out = 0;
  }
  if (solver == NULL || solver->symbolic == NULL ||
      solver->numeric == NULL || solver->symbolic->R == NULL ||
      solver->numeric->LUbx == NULL || solver->numeric->Llen == NULL ||
      solver->numeric->Uip == NULL || solver->numeric->Ulen == NULL ||
      solver->n == 0u) {
    return;
  }

  int threads = thread_count > 0 ? thread_count : 1;
  if (threads < 1) {
    threads = 1;
  }
  if ((UF_long)threads > solver->n) {
    threads = (int)solver->n;
  }
  if (threads < 1 || solver->n > (UF_long)(SIZE_MAX / sizeof(double)) ||
      (size_t)threads > SIZE_MAX / sizeof(double)) {
    return;
  }

  double *finish =
    (double *)calloc((size_t)solver->n, sizeof(*finish));
  double *thread_end =
    (double *)calloc((size_t)threads, sizeof(*thread_end));
  if (finish == NULL || thread_end == NULL) {
    free(finish);
    free(thread_end);
    return;
  }

  double work = 0.0;
  double finish_time = 0.0;
  int64_t dependencies = 0;
  int valid = 1;
  for (UF_long block = 0;
       valid && block < solver->symbolic->nblocks; ++block) {
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    if (k1 > k2 || k2 > solver->n || block >= solver->numeric->nblocks) {
      valid = 0;
      break;
    }
    const UF_long nk = k2 - k1;
    double *lu = (double *)solver->numeric->LUbx[block];
    if (nk > 1u && lu == NULL) {
      valid = 0;
      break;
    }
    const UF_long *uip = solver->numeric->Uip + k1;
    const UF_long *ulen = solver->numeric->Ulen + k1;
    for (UF_long local_k = 0; valid && local_k < nk; ++local_k) {
      const UF_long k = k1 + local_k;
      int best_thread = 0;
      double best_end = thread_end[0];
      for (int p = 1; p < threads; ++p) {
        if (thread_end[p] < best_end) {
          best_end = thread_end[p];
          best_thread = p;
        }
      }

      double ft = best_end;
      UF_long *ui = NULL;
      double *ux = NULL;
      UF_long ucol_len = 0;
      if (lu != NULL) {
        kls_klu_get_pointer(lu, uip, ulen, local_k, &ui, &ux, &ucol_len);
      } else if (solver->numeric->Ulen[k] != 0u) {
        valid = 0;
        break;
      }
      (void)ux;
      for (UF_long p = 0; p < ucol_len; ++p) {
        const UF_long dep_local = ui[p];
        if (dep_local >= local_k) {
          valid = 0;
          break;
        }
        const UF_long dep = k1 + dep_local;
        const double update_work =
          2.0 * ((double)solver->numeric->Llen[dep] + 1.0);
        const double ready_time = finish[dep] > ft ? finish[dep] : ft;
        ft = ready_time + update_work + KLS_NICSLU_TASK_FLOW_SYNC_COST;
        work += update_work;
        dependencies++;
      }
      if (!valid) {
        break;
      }
      const double norm_work = (double)solver->numeric->Llen[k] + 1.0;
      ft += norm_work;
      work += norm_work;
      finish[k] = ft;
      thread_end[best_thread] = ft;
      if (ft > finish_time) {
        finish_time = ft;
      }
    }
  }

  if (valid && finish_time > 0.0 && work > 0.0) {
    const double speedup = work / finish_time;
    if (threads_out != NULL) {
      *threads_out = (int64_t)threads;
    }
    if (dependencies_out != NULL) {
      *dependencies_out = dependencies;
    }
    if (work_out != NULL) {
      *work_out = work;
    }
    if (finish_time_out != NULL) {
      *finish_time_out = finish_time;
    }
    if (speedup_out != NULL) {
      *speedup_out = speedup;
    }
    if (recommends_parallel_out != NULL) {
      *recommends_parallel_out = threads > 1 && speedup > 1.0;
    }
  }

  free(finish);
  free(thread_end);
}

static void kls_update_parallel_model_stats(kls_solver *solver,
                                            int prefer_numeric) {
  if (solver == NULL) {
    return;
  }
  kls_compute_parallel_model(
    solver, prefer_numeric, &solver->stats.parallel_model_r1,
    &solver->stats.parallel_model_r2,
    &solver->stats.parallel_model_recommends_parallel);
}

static void kls_update_nicslu_task_flow_model_stats(kls_solver *solver) {
  if (solver == NULL) {
    return;
  }
  kls_compute_nicslu_task_flow_model(
    solver, solver->options.threads,
    &solver->stats.parallel_task_flow_threads,
    &solver->stats.parallel_task_flow_dependencies,
    &solver->stats.parallel_task_flow_work,
    &solver->stats.parallel_task_flow_finish_time,
    &solver->stats.parallel_task_flow_speedup,
    &solver->stats.parallel_task_flow_recommends_parallel);
}

static void kls_fill_fast_kls_block_restart_pipeline_stats(kls_solver *solver) {
  if (solver == NULL) {
    return;
  }
  solver->stats.fast_kls_block_restart_last_row_pipeline =
    solver->fast_kls_block_restart_last_row_pipeline;
  solver->stats.fast_kls_block_restart_row_pipeline_count =
    (int64_t)solver->fast_kls_block_restart_row_pipeline_count;
  solver->stats.fast_kls_block_restart_last_row_pipeline_rows =
    (int64_t)solver->fast_kls_block_restart_last_row_pipeline_rows;
  solver->stats.fast_kls_block_restart_last_row_pipeline_threads =
    (int64_t)solver->fast_kls_block_restart_last_row_pipeline_threads;
  solver->stats.fast_kls_block_restart_last_row_pipeline_prefix_rows =
    (int64_t)solver->fast_kls_block_restart_last_row_pipeline_prefix_rows;
  solver->stats.fast_kls_block_restart_last_row_pipeline_suffix_rows =
    (int64_t)solver->fast_kls_block_restart_last_row_pipeline_suffix_rows;
  solver->stats.fast_kls_block_restart_last_row_pipeline_gap_rows =
    (int64_t)solver->fast_kls_block_restart_last_row_pipeline_gap_rows;
  solver->stats.fast_kls_block_restart_last_row_pipeline_etree_tail =
    solver->fast_kls_block_restart_last_row_pipeline_etree_tail;
  solver->stats.fast_kls_block_restart_row_pipeline_etree_tail_count =
    (int64_t)solver->fast_kls_block_restart_row_pipeline_etree_tail_count;
  solver->stats.fast_kls_block_restart_last_row_pipeline_etree_tail_rows =
    (int64_t)
      solver->fast_kls_block_restart_last_row_pipeline_etree_tail_rows;
  solver->stats.fast_kls_block_restart_last_row_pipeline_etree_tail_gap_rows =
    (int64_t)
      solver->fast_kls_block_restart_last_row_pipeline_etree_tail_gap_rows;
  solver->stats
    .fast_kls_block_restart_last_row_pipeline_etree_tail_exact_mask =
      solver
        ->fast_kls_block_restart_last_row_pipeline_etree_tail_exact_mask;
  solver->stats
    .fast_kls_block_restart_last_row_pipeline_etree_ready =
      solver->fast_kls_block_restart_last_row_pipeline_etree_ready;
  solver->stats
    .fast_kls_block_restart_row_pipeline_etree_ready_count =
      (int64_t)
        solver->fast_kls_block_restart_row_pipeline_etree_ready_count;
  solver->stats
    .fast_kls_block_restart_last_row_pipeline_etree_ready_rows =
      (int64_t)
        solver->fast_kls_block_restart_last_row_pipeline_etree_ready_rows;
  solver->stats
    .fast_kls_block_restart_last_row_pipeline_etree_ready_threads =
      (int64_t)
        solver->fast_kls_block_restart_last_row_pipeline_etree_ready_threads;
  solver->stats
    .fast_kls_block_restart_last_row_pipeline_etree_prefactor =
      solver->fast_kls_block_restart_last_row_pipeline_etree_prefactor;
  solver->stats
    .fast_kls_block_restart_row_pipeline_etree_prefactor_count =
      (int64_t)
        solver->fast_kls_block_restart_row_pipeline_etree_prefactor_count;
  solver->stats
    .fast_kls_block_restart_last_row_pipeline_etree_prefactor_rows =
      (int64_t)
        solver->fast_kls_block_restart_last_row_pipeline_etree_prefactor_rows;
  solver->stats
    .fast_kls_block_restart_last_row_pipeline_etree_prefactor_threads =
      (int64_t)
        solver->fast_kls_block_restart_last_row_pipeline_etree_prefactor_threads;
  solver->stats
    .fast_kls_block_restart_last_row_pipeline_etree_prefactor_wait_rows =
      (int64_t)
        solver
          ->fast_kls_block_restart_last_row_pipeline_etree_prefactor_wait_rows;
  solver->stats
    .fast_kls_block_restart_last_row_pipeline_etree_prefactor_wait_deps =
      (int64_t)
        solver
          ->fast_kls_block_restart_last_row_pipeline_etree_prefactor_wait_deps;
  solver->stats
    .fast_kls_block_restart_last_row_pipeline_separator_tail_scope =
      solver
        ->fast_kls_block_restart_last_row_pipeline_separator_tail_scope;
  solver->stats
    .fast_kls_block_restart_row_pipeline_separator_tail_scope_count =
      (int64_t)
        solver
          ->fast_kls_block_restart_row_pipeline_separator_tail_scope_count;
  solver->stats
    .fast_kls_block_restart_last_row_pipeline_separator_tail_scope_rows =
      (int64_t)
        solver
          ->fast_kls_block_restart_last_row_pipeline_separator_tail_scope_rows;
  solver->stats
    .fast_kls_block_restart_last_row_pipeline_separator_queue =
      solver->fast_kls_block_restart_last_row_pipeline_separator_queue;
  solver->stats
    .fast_kls_block_restart_row_pipeline_separator_queue_count =
      (int64_t)
        solver->fast_kls_block_restart_row_pipeline_separator_queue_count;
  solver->stats
    .fast_kls_block_restart_last_row_pipeline_separator_private_rows =
      (int64_t)
        solver->fast_kls_block_restart_last_row_pipeline_separator_private_rows;
  solver->stats
    .fast_kls_block_restart_last_row_pipeline_separator_pipeline_rows =
      (int64_t)
        solver->fast_kls_block_restart_last_row_pipeline_separator_pipeline_rows;
  solver->stats
    .fast_kls_block_restart_last_row_pipeline_separator_private_threads =
      (int64_t)
        solver
          ->fast_kls_block_restart_last_row_pipeline_separator_private_threads;
  solver->stats
    .fast_kls_block_restart_last_row_pipeline_separator_partitioned =
      solver->fast_kls_block_restart_last_row_pipeline_separator_partitioned;
  solver->stats
    .fast_kls_block_restart_last_row_pipeline_separator_split_components =
      (int64_t)
        solver
          ->fast_kls_block_restart_last_row_pipeline_separator_split_components;
  solver->stats
    .fast_kls_block_restart_last_row_pipeline_pivot_tail_rows =
      (int64_t)
        solver->fast_kls_block_restart_last_row_pipeline_pivot_tail_rows;
  solver->stats
    .fast_kls_block_restart_last_row_pipeline_pivot_restarts =
      (int64_t)
        solver->fast_kls_block_restart_last_row_pipeline_pivot_restarts;
  solver->stats
    .fast_kls_block_restart_last_row_pipeline_supernode_update_groups =
      (int64_t)
        solver
          ->fast_kls_block_restart_last_row_pipeline_supernode_update_groups;
  solver->stats
    .fast_kls_block_restart_last_row_pipeline_supernode_update_rows =
      (int64_t)
        solver->fast_kls_block_restart_last_row_pipeline_supernode_update_rows;
  solver->stats
    .fast_kls_block_restart_last_row_pipeline_supernode_panel_update_groups =
      (int64_t)
        solver
          ->fast_kls_block_restart_last_row_pipeline_supernode_panel_update_groups;
  solver->stats
    .fast_kls_block_restart_last_row_pipeline_supernode_panel_update_rows =
      (int64_t)
        solver
          ->fast_kls_block_restart_last_row_pipeline_supernode_panel_update_rows;
}

static void fill_symbolic_stats(kls_solver *solver, double elapsed) {
  solver->stats.struct_size = sizeof(solver->stats);
  fill_build_stats(&solver->stats);
  solver->stats.n = (int64_t)solver->n;
  solver->stats.nnz = (int64_t)solver->nnz;
  solver->stats.analysis_seconds = elapsed;
  solver->stats.selected_orientation = solver->orientation;
  solver->stats.selected_scale = (int)solver->common.scale;
  solver->stats.selected_pivot_tolerance = solver->common.tol;
  solver->stats.selected_static_pivoting = solver->row_perm != NULL;
  solver->stats.selected_exact_matching = solver->exact_matching_selected;
  solver->stats.selected_exact_matching_scaling =
    solver->exact_matching_scaling_selected;
  solver->stats.selected_spral_matching = solver->spral_matching_selected;
  solver->stats.fast_block_restarts = solver->fast_block_restarts;
  solver->stats.fast_kls_block_restarts = solver->fast_kls_block_restarts;
  solver->stats.fast_kls_rebuild_restarts =
    solver->fast_kls_rebuild_restarts;
  kls_fill_fast_kls_block_restart_pipeline_stats(solver);
  solver->stats.fast_tail_restarts = solver->fast_tail_restarts;
  solver->stats.fast_repaired_last_offdiag_suffix_refresh =
    solver->fast_repaired_last_offdiag_suffix_refresh;
  solver->stats.fast_repaired_offdiag_suffix_refresh_count =
    (int64_t)solver->fast_repaired_offdiag_suffix_refresh_count;
  solver->stats.fast_repaired_offdiag_full_refresh_count =
    (int64_t)solver->fast_repaired_offdiag_full_refresh_count;
  solver->stats.fast_repaired_parallel_tail_blocks =
    (int64_t)solver->fast_repaired_parallel_tail_blocks;
  solver->stats.fast_rejected_refresh_state =
    solver->fast_reject_refresh_state;
  solver->stats.fast_rejected_prefix_refresh_columns =
    (int64_t)solver->fast_rejected_prefix_refresh_columns;
  solver->stats.fast_rejected_prefix_refresh_count =
    (int64_t)solver->fast_rejected_prefix_refresh_count;
  solver->stats.row_refactor_total_group_work =
    kls_row_refactor_total_group_work(solver);
  solver->stats.row_refactor_auto_enabled =
    solver->row_refactor_auto_enabled;
  solver->stats.row_refactor_auto_values_ready =
    solver->row_refactor_values_ready;
  solver->stats.row_refactor_auto_work_allowed =
    kls_auto_row_refactor_cost_allows(solver);
  solver->stats.row_refactor_auto_should_run =
    kls_auto_row_refactor_should_run(solver);
  solver->stats.row_refactor_auto_lower_bound_work =
    solver->row_refactor_auto_lower_bound_work;
  solver->stats.row_refactor_auto_lower_bound_rejected =
    solver->row_refactor_auto_lower_bound_rejected;
  solver->stats.row_refactor_auto_pattern_build_failed =
    solver->row_refactor_auto_pattern_build_failed;
  solver->stats.row_refactor_auto_value_copy_failed =
    solver->row_refactor_auto_value_copy_failed;
  kls_fill_separator_stats(&solver->stats, &solver->separator);
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
  if (!solver->parallel_model_stats_valid) {
    /* Both models read only the numeric's pattern and the thread count;
       recomputing them on every refactorization costs ~10% of a small
       matrix's refactor. */
    kls_update_parallel_model_stats(solver, 0);
    kls_update_nicslu_task_flow_model_stats(solver);
    solver->parallel_model_stats_valid = 1;
  }
  kls_update_factor_etree_stats(solver);
  solver->stats.memory_bytes = solver->common.memusage;
  solver->stats.memory_peak_bytes = solver->common.mempeak;
}

static void fill_numeric_stats(kls_solver *solver) {
  {
    /* A numeric-stat refresh is an observability operation.  It must not
       silently classify the retained factor and arm a benchmark-family
       kernel for the next refactor/solve.  Use -1 (known ineligible), not
       zero (unknown), because deprecated policy fields must remain inert. */
    solver->symmetric_scalar_fringe_numeric_eligible = -1;
    solver->pivoted_high_work_single_block_numeric_eligible = -1;
    solver->low_work_many_fringe_btf_pts_numeric_eligible = -1;
    solver->sparse_symmetric_fragmented_metis_numeric_eligible = -1;
    solver->sparse_spiked_predicted_numeric_eligible = -1;
    solver->sparse_full_diagonal_metis_row_numeric_eligible = -1;
    solver->giant_symmetric_scalar_fringe_metis_row_numeric_eligible = -1;
    solver->asymmetric_bounded_degree_direct_metis_numeric_eligible = -1;
    solver->near_symmetric_mega_hub_amd_numeric_eligible = -1;
    solver->giant_dominant_hub_metis_dense_tail_numeric_eligible = -1;
    solver->hybrid_huge_single_egraph_numeric_eligible = -1;
    solver->high_work_tiny_scalar_fringe_amd_numeric_eligible = -1;
  }
  kls_dbg_snapshot_numeric(solver);
  fill_build_stats(&solver->stats);
  solver->stats.refactor_lean_choice = solver->lean_choice;
  solver->stats.egraph_worker_spin_iters = solver->egraph_pool != NULL
    ? (int64_t)solver->egraph_pool->worker_spin_iters : 0;
  solver->stats.last_kernel_status = (int)solver->common.status;
  solver->stats.selected_scale = (int)solver->common.scale;
  solver->stats.selected_pivot_tolerance = solver->common.tol;
  solver->stats.selected_static_pivoting = solver->row_perm != NULL;
  solver->stats.selected_exact_matching = solver->exact_matching_selected;
  solver->stats.selected_exact_matching_scaling =
    solver->exact_matching_scaling_selected;
  solver->stats.selected_spral_matching = solver->spral_matching_selected;
  solver->stats.fast_block_restarts = solver->fast_block_restarts;
  solver->stats.fast_kls_block_restarts = solver->fast_kls_block_restarts;
  solver->stats.fast_kls_rebuild_restarts =
    solver->fast_kls_rebuild_restarts;
  kls_fill_fast_kls_block_restart_pipeline_stats(solver);
  solver->stats.fast_tail_restarts = solver->fast_tail_restarts;
  solver->stats.fast_repaired_last_offdiag_suffix_refresh =
    solver->fast_repaired_last_offdiag_suffix_refresh;
  solver->stats.fast_repaired_offdiag_suffix_refresh_count =
    (int64_t)solver->fast_repaired_offdiag_suffix_refresh_count;
  solver->stats.fast_repaired_offdiag_full_refresh_count =
    (int64_t)solver->fast_repaired_offdiag_full_refresh_count;
  solver->stats.fast_repaired_parallel_tail_blocks =
    (int64_t)solver->fast_repaired_parallel_tail_blocks;
  solver->stats.fast_rejected_refresh_state =
    solver->fast_reject_refresh_state;
  solver->stats.fast_rejected_prefix_refresh_columns =
    (int64_t)solver->fast_rejected_prefix_refresh_columns;
  solver->stats.fast_rejected_prefix_refresh_count =
    (int64_t)solver->fast_rejected_prefix_refresh_count;
  solver->stats.row_refactor_total_group_work =
    kls_row_refactor_total_group_work(solver);
  solver->stats.row_refactor_auto_enabled =
    solver->row_refactor_auto_enabled;
  solver->stats.row_refactor_auto_values_ready =
    solver->row_refactor_values_ready;
  solver->stats.row_refactor_auto_work_allowed =
    kls_auto_row_refactor_cost_allows(solver);
  solver->stats.row_refactor_auto_should_run =
    kls_auto_row_refactor_should_run(solver);
  solver->stats.row_refactor_auto_lower_bound_work =
    solver->row_refactor_auto_lower_bound_work;
  solver->stats.row_refactor_auto_lower_bound_rejected =
    solver->row_refactor_auto_lower_bound_rejected;
  solver->stats.row_refactor_auto_pattern_build_failed =
    solver->row_refactor_auto_pattern_build_failed;
  solver->stats.row_refactor_auto_value_copy_failed =
    solver->row_refactor_auto_value_copy_failed;
  kls_fill_separator_stats(&solver->stats, &solver->separator);
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
  if (!solver->parallel_model_stats_valid) {
    kls_update_parallel_model_stats(solver, 1);
    kls_update_nicslu_task_flow_model_stats(solver);
    solver->parallel_model_stats_valid = 1;
  }
  kls_update_factor_etree_stats(solver);
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

  solver->stats.refactor_last_egraph_algorithm5_prefactor_columns =
    (int64_t)solver->refactor_last_egraph_algorithm5_prefactor_columns;
  solver->stats.refactor_last_egraph_algorithm5_prefactor_deps =
    (int64_t)solver->refactor_last_egraph_algorithm5_prefactor_deps;
  solver->stats.refactor_egraph_algorithm5_prefactor_column_count =
    (int64_t)solver->refactor_egraph_algorithm5_prefactor_column_count;
  solver->stats.refactor_egraph_algorithm5_prefactor_dep_count =
    (int64_t)solver->refactor_egraph_algorithm5_prefactor_dep_count;

  solver->stats.refactor_last_supernode_algorithm5_payoff_slot_accum_runs =
    0;
  solver->stats.refactor_last_supernode_algorithm5_payoff_slot_accum_rows =
    0;
  solver->stats
    .refactor_last_supernode_algorithm5_payoff_slot_accum_target_entries =
      0;
  solver->stats
    .refactor_last_supernode_algorithm5_payoff_slot_accum_target_slots =
      0;
  solver->stats.refactor_last_supernode_algorithm5_payoff_prefix_prep_runs =
    0;
  solver->stats.refactor_last_supernode_algorithm5_payoff_prefix_prep_rows =
    0;
  solver->stats
    .refactor_last_supernode_algorithm5_payoff_prefix_prep_target_entries =
      0;
  solver->stats
    .refactor_last_supernode_algorithm5_payoff_prefix_prep_target_slots =
      0;
  solver->stats.refactor_last_supernode_algorithm5_payoff_advance_seed_runs =
    0;
  solver->stats.refactor_last_supernode_algorithm5_payoff_advance_seed_deps =
    0;
  solver->stats.refactor_last_supernode_algorithm5_payoff_advance_seed_slots =
    0;
  solver->stats
    .refactor_last_supernode_algorithm5_payoff_current_state_seed_runs =
      0;
  solver->stats
    .refactor_last_supernode_algorithm5_payoff_current_state_seed_deps =
      0;
  solver->stats
    .refactor_last_supernode_algorithm5_payoff_current_state_seed_rows =
      0;
  solver->stats
    .refactor_last_supernode_algorithm5_payoff_final_trigger_batches =
      0;
  solver->stats
    .refactor_last_supernode_algorithm5_payoff_final_trigger_multi_batches =
      0;
  solver->stats
    .refactor_last_supernode_algorithm5_payoff_final_trigger_claims =
      0;
  solver->stats
    .refactor_last_supernode_algorithm5_payoff_final_trigger_suffix_deps =
      0;
  solver->stats
    .refactor_last_supernode_algorithm5_payoff_suffix_advance_slots =
      0;
  solver->stats
    .refactor_last_supernode_algorithm5_payoff_suffix_advance_deps =
      0;
  solver->stats
    .refactor_last_supernode_algorithm5_payoff_suffix_advance_updates =
      0;
  solver->stats
    .refactor_last_supernode_algorithm5_payoff_suffix_advance_finished =
      0;

  solver->stats.refactor_supernode_panel_count =
    (int64_t)solver->refactor_supernode_panel_count;
  solver->stats.refactor_supernode_panel_used_count =
    (int64_t)solver->refactor_supernode_panel_used_count;

  solver->stats.refactor_map_index32_enabled =
    solver->refactor_map_index32_enabled;
  solver->stats.refactor_map_index32_entries =
    (int64_t)solver->refactor_map_indices32_count;
  solver->stats.refactor_l_index32_enabled =
    solver->refactor_l_index32_enabled;
  solver->stats.refactor_l_index32_entries =
    (int64_t)solver->refactor_l_indices32_count;
  solver->stats.refactor_u_index32_enabled =
    solver->refactor_u_index32_enabled;
  solver->stats.refactor_u_index32_entries =
    (int64_t)solver->refactor_u_indices32_count;
  solver->stats.refactor_last_supernode_pipeline_tasks = 0;
  solver->stats.refactor_last_supernode_pipeline_columns = 0;
  solver->stats.refactor_supernode_pipeline_task_count = 0;
  solver->stats.refactor_supernode_pipeline_column_count = 0;
  solver->stats.refactor_last_supernode_update_runs =
    (int64_t)solver->refactor_last_supernode_update_runs;
  solver->stats.refactor_last_supernode_update_rows =
    (int64_t)solver->refactor_last_supernode_update_rows;
  solver->stats.refactor_last_supernode_update_entries =
    (int64_t)solver->refactor_last_supernode_update_entries;
  solver->stats.refactor_supernode_update_run_count =
    (int64_t)solver->refactor_supernode_update_run_count;
  solver->stats.refactor_supernode_update_rows =
    (int64_t)solver->refactor_supernode_update_rows;
  solver->stats.refactor_supernode_update_entries =
    (int64_t)solver->refactor_supernode_update_entries;
  solver->stats.refactor_last_supernode_cblas_update_runs =
    (int64_t)solver->refactor_last_supernode_cblas_update_runs;
  solver->stats.refactor_last_supernode_cblas_update_rows =
    (int64_t)solver->refactor_last_supernode_cblas_update_rows;
  solver->stats.refactor_last_supernode_cblas_update_entries =
    (int64_t)solver->refactor_last_supernode_cblas_update_entries;
  solver->stats.refactor_supernode_cblas_update_run_count =
    (int64_t)solver->refactor_supernode_cblas_update_run_count;
  solver->stats.refactor_supernode_cblas_update_rows =
    (int64_t)solver->refactor_supernode_cblas_update_rows;
  solver->stats.refactor_supernode_cblas_update_entries =
    (int64_t)solver->refactor_supernode_cblas_update_entries;
  solver->stats.refactor_last_supernode_blocked_update_runs =
    (int64_t)solver->refactor_last_supernode_blocked_update_runs;
  solver->stats.refactor_last_supernode_blocked_update_rows =
    (int64_t)solver->refactor_last_supernode_blocked_update_rows;
  solver->stats.refactor_last_supernode_blocked_update_entries =
    (int64_t)solver->refactor_last_supernode_blocked_update_entries;
  solver->stats.refactor_supernode_blocked_update_run_count =
    (int64_t)solver->refactor_supernode_blocked_update_run_count;
  solver->stats.refactor_supernode_blocked_update_rows =
    (int64_t)solver->refactor_supernode_blocked_update_rows;
  solver->stats.refactor_supernode_blocked_update_entries =
    (int64_t)solver->refactor_supernode_blocked_update_entries;
  solver->stats.refactor_last_supernode_cached_probe_attempts =
    (int64_t)solver->refactor_last_supernode_cached_probe_attempts;
  solver->stats.refactor_last_supernode_cached_probe_panel_hits =
    (int64_t)solver->refactor_last_supernode_cached_probe_panel_hits;
  solver->stats.refactor_last_supernode_cached_probe_contiguous =
    (int64_t)solver->refactor_last_supernode_cached_probe_contiguous;
  solver->stats.refactor_last_supernode_cached_probe_allowed =
    (int64_t)solver->refactor_last_supernode_cached_probe_allowed;
  solver->stats.refactor_last_supernode_cached_probe_allowed_rows =
    (int64_t)solver->refactor_last_supernode_cached_probe_allowed_rows;
  solver->stats.refactor_last_supernode_cached_probe_applied =
    (int64_t)solver->refactor_last_supernode_cached_probe_applied;
  solver->stats.refactor_last_supernode_cached_probe_applied_rows =
    (int64_t)solver->refactor_last_supernode_cached_probe_applied_rows;
  solver->stats.refactor_last_supernode_cached_probe_shape_rejects =
    (int64_t)solver->refactor_last_supernode_cached_probe_shape_rejects;
  solver->stats.refactor_last_supernode_cached_probe_shape_reject_rows =
    (int64_t)solver->refactor_last_supernode_cached_probe_shape_reject_rows;
  solver->stats.refactor_last_supernode_cached_probe_stream_rejects =
    (int64_t)solver->refactor_last_supernode_cached_probe_stream_rejects;
  solver->stats.refactor_last_supernode_cached_probe_stream_reject_rows =
    (int64_t)solver->refactor_last_supernode_cached_probe_stream_reject_rows;
  solver->stats.refactor_last_supernode_cached_probe_work_rejects =
    (int64_t)solver->refactor_last_supernode_cached_probe_work_rejects;
  solver->stats.refactor_last_supernode_cached_probe_work_reject_rows =
    (int64_t)solver->refactor_last_supernode_cached_probe_work_reject_rows;
  solver->stats.refactor_last_supernode_cached_probe_workspace_rejects =
    (int64_t)solver->refactor_last_supernode_cached_probe_workspace_rejects;
  solver->stats.refactor_last_supernode_cached_probe_workspace_reject_rows =
    (int64_t)solver->refactor_last_supernode_cached_probe_workspace_reject_rows;
  solver->stats.refactor_supernode_cached_probe_attempt_count =
    (int64_t)solver->refactor_supernode_cached_probe_attempt_count;
  solver->stats.refactor_supernode_cached_probe_panel_hits =
    (int64_t)solver->refactor_supernode_cached_probe_panel_hits;
  solver->stats.refactor_supernode_cached_probe_contiguous =
    (int64_t)solver->refactor_supernode_cached_probe_contiguous;
  solver->stats.refactor_supernode_cached_probe_allowed =
    (int64_t)solver->refactor_supernode_cached_probe_allowed;
  solver->stats.refactor_supernode_cached_probe_allowed_rows =
    (int64_t)solver->refactor_supernode_cached_probe_allowed_rows;
  solver->stats.refactor_supernode_cached_probe_applied =
    (int64_t)solver->refactor_supernode_cached_probe_applied;
  solver->stats.refactor_supernode_cached_probe_applied_rows =
    (int64_t)solver->refactor_supernode_cached_probe_applied_rows;
  solver->stats.refactor_supernode_cached_probe_shape_rejects =
    (int64_t)solver->refactor_supernode_cached_probe_shape_rejects;
  solver->stats.refactor_supernode_cached_probe_shape_reject_rows =
    (int64_t)solver->refactor_supernode_cached_probe_shape_reject_rows;
  solver->stats.refactor_supernode_cached_probe_stream_rejects =
    (int64_t)solver->refactor_supernode_cached_probe_stream_rejects;
  solver->stats.refactor_supernode_cached_probe_stream_reject_rows =
    (int64_t)solver->refactor_supernode_cached_probe_stream_reject_rows;
  solver->stats.refactor_supernode_cached_probe_work_rejects =
    (int64_t)solver->refactor_supernode_cached_probe_work_rejects;
  solver->stats.refactor_supernode_cached_probe_work_reject_rows =
    (int64_t)solver->refactor_supernode_cached_probe_work_reject_rows;
  solver->stats.refactor_supernode_cached_probe_workspace_rejects =
    (int64_t)solver->refactor_supernode_cached_probe_workspace_rejects;
  solver->stats.refactor_supernode_cached_probe_workspace_reject_rows =
    (int64_t)solver->refactor_supernode_cached_probe_workspace_reject_rows;

  solver->stats.refactor_last_btf_scalar_run_exec_runs =
    (int64_t)solver->refactor_last_btf_scalar_run_exec_runs;
  solver->stats.refactor_last_btf_scalar_run_exec_rows =
    (int64_t)solver->refactor_last_btf_scalar_run_exec_rows;
  solver->stats.refactor_last_btf_scalar_run_exec_entries =
    (int64_t)solver->refactor_last_btf_scalar_run_exec_entries;
  solver->stats.refactor_last_btf_scalar_run_exec_max_rows =
    (int64_t)solver->refactor_last_btf_scalar_run_exec_max_rows;
  solver->stats.refactor_btf_scalar_run_exec_count =
    (int64_t)solver->refactor_btf_scalar_run_exec_count;
  solver->stats.refactor_btf_scalar_run_exec_rows =
    (int64_t)solver->refactor_btf_scalar_run_exec_rows;
  solver->stats.refactor_btf_scalar_run_exec_entries =
    (int64_t)solver->refactor_btf_scalar_run_exec_entries;
  solver->stats.refactor_btf_scalar_run_exec_max_rows =
    (int64_t)solver->refactor_btf_scalar_run_exec_max_rows;

  solver->stats.refactor_supernode_cached_probe_disabled =
    solver->refactor_supernode_cached_probe_disabled;
  solver->stats.refactor_supernode_cached_probe_disable_count =
    (int64_t)solver->refactor_supernode_cached_probe_disable_count;
  solver->stats.refactor_supernode_update_disabled =
    solver->refactor_supernode_update_disabled;
  solver->stats.refactor_supernode_update_disable_count =
    (int64_t)solver->refactor_supernode_update_disable_count;
  solver->stats.refactor_last_ready_queue_columns =
    (int64_t)solver->refactor_last_ready_queue_columns;
  solver->stats.refactor_ready_queue_run_count =
    (int64_t)solver->refactor_ready_queue_run_count;
  solver->stats.row_refactor_group_count =
    (int64_t)solver->row_refactor_group_count;
  solver->stats.row_refactor_group_single_count =
    (int64_t)solver->row_refactor_group_single_count;
  solver->stats.row_refactor_group_batch_count =
    (int64_t)solver->row_refactor_group_batch_count;
  solver->stats.row_refactor_group_batch_rows =
    (int64_t)solver->row_refactor_group_batch_rows;
  solver->stats.row_refactor_group_batch_max_width =
    (int64_t)solver->row_refactor_group_batch_max_width;
  solver->stats.row_refactor_group_batch_width_le_4_count =
    (int64_t)solver->row_refactor_group_batch_width_le_4_count;
  solver->stats.row_refactor_group_batch_width_le_8_count =
    (int64_t)solver->row_refactor_group_batch_width_le_8_count;
  solver->stats.row_refactor_group_scalar_candidate_count =
    (int64_t)solver->row_refactor_group_scalar_candidate_count;
  solver->stats.row_refactor_group_scalar_candidate_rows =
    (int64_t)solver->row_refactor_group_scalar_candidate_rows;
  solver->stats.row_refactor_group_scalar_short_count =
    (int64_t)solver->row_refactor_group_scalar_short_count;
  solver->stats.row_refactor_group_scalar_short_rows =
    (int64_t)solver->row_refactor_group_scalar_short_rows;
  solver->stats.row_refactor_group_scalar_stop_level_mismatch_count =
    (int64_t)solver->row_refactor_group_scalar_stop_level_mismatch_count;
  solver->stats.row_refactor_group_scalar_stop_internal_dep_count =
    (int64_t)solver->row_refactor_group_scalar_stop_internal_dep_count;
  solver->stats.row_refactor_group_scalar_stop_next_segment_count =
    (int64_t)solver->row_refactor_group_scalar_stop_next_segment_count;
  solver->stats.row_refactor_group_scalar_stop_max_width_count =
    (int64_t)solver->row_refactor_group_scalar_stop_max_width_count;
  solver->stats.row_refactor_group_scalar_stop_matrix_end_count =
    (int64_t)solver->row_refactor_group_scalar_stop_matrix_end_count;
  solver->stats.row_refactor_group_generic_count =
    (int64_t)solver->row_refactor_group_generic_count;
  solver->stats.row_refactor_group_generic_rows =
    (int64_t)solver->row_refactor_group_generic_rows;
  solver->stats.row_refactor_group_generic_max_width =
    (int64_t)solver->row_refactor_group_generic_max_width;
  solver->stats.row_refactor_group_dense_count =
    (int64_t)solver->row_refactor_group_dense_count;
  solver->stats.row_refactor_group_dense_rows =
    (int64_t)solver->row_refactor_group_dense_rows;
  solver->stats.row_refactor_group_dense_max_width =
    (int64_t)solver->row_refactor_group_dense_max_width;
  solver->stats.row_refactor_group_single_work =
    solver->row_refactor_group_single_work;
  solver->stats.row_refactor_group_batch_work =
    solver->row_refactor_group_batch_work;
  solver->stats.row_refactor_group_generic_work =
    solver->row_refactor_group_generic_work;
  solver->stats.row_refactor_group_dense_work =
    solver->row_refactor_group_dense_work;
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
  solver->stats.row_refactor_group_dependency_edges =
    (int64_t)solver->row_refactor_group_dependency_edges;
  solver->stats.row_refactor_group_root_count =
    (int64_t)solver->row_refactor_group_root_count;
  solver->stats.row_refactor_group_leaf_count =
    (int64_t)solver->row_refactor_group_leaf_count;
  solver->stats.row_refactor_group_max_fanout =
    (int64_t)solver->row_refactor_group_max_fanout;
  solver->stats.row_refactor_last_ready_queue =
    solver->row_refactor_last_ready_queue;
  solver->stats.row_refactor_ready_queue_run_count =
    (int64_t)solver->row_refactor_ready_queue_run_count;
  solver->stats.row_refactor_ready_queue_group_count =
    (int64_t)solver->row_refactor_ready_queue_group_count;
  solver->stats.row_refactor_last_work_ready_queue =
    solver->row_refactor_last_work_ready_queue;
  solver->stats.row_refactor_work_ready_queue_run_count =
    (int64_t)solver->row_refactor_work_ready_queue_run_count;
  solver->stats.row_refactor_last_partial_supernode_pipeline =
    solver->row_refactor_last_partial_supernode_pipeline;
  solver->stats.row_refactor_last_partial_supernode_pipeline_groups =
    (int64_t)solver->row_refactor_last_partial_supernode_pipeline_groups;
  solver->stats.row_refactor_last_partial_supernode_pipeline_rows =
    (int64_t)solver->row_refactor_last_partial_supernode_pipeline_rows;
  solver->stats.row_refactor_partial_supernode_pipeline_run_count =
    (int64_t)solver->row_refactor_partial_supernode_pipeline_run_count;
  solver->stats.row_refactor_last_compact_dense_panel =
    solver->row_refactor_last_compact_dense_panel;
  solver->stats.row_refactor_compact_dense_panel_count =
    (int64_t)solver->row_refactor_compact_dense_panel_count;
  solver->stats.row_refactor_ready_queue_workspace_groups =
    (int64_t)solver->row_refactor_ready_queue_capacity;
  solver->stats.row_refactor_last_local_ready_groups =
    (int64_t)solver->row_refactor_last_local_ready_groups;
  solver->stats.row_refactor_local_ready_group_count =
    (int64_t)solver->row_refactor_local_ready_group_count;
  solver->stats.row_refactor_last_private_ready_groups =
    (int64_t)solver->row_refactor_last_private_ready_groups;
  solver->stats.row_refactor_private_ready_group_count =
    (int64_t)solver->row_refactor_private_ready_group_count;
  solver->stats.row_refactor_last_separator_private_queue =
    solver->row_refactor_last_separator_private_queue;
  solver->stats.row_refactor_separator_private_queue_run_count =
    (int64_t)solver->row_refactor_separator_private_queue_run_count;
  solver->stats.row_refactor_last_separator_private_components =
    (int64_t)solver->row_refactor_last_separator_private_components;
  solver->stats.row_refactor_separator_private_component_count =
    (int64_t)solver->row_refactor_separator_private_component_count;
  solver->stats.row_refactor_last_separator_flop_queue =
    solver->row_refactor_last_separator_flop_queue;
  solver->stats.row_refactor_separator_flop_queue_run_count =
    (int64_t)solver->row_refactor_separator_flop_queue_run_count;
  solver->stats.row_refactor_last_separator_flop_ordered_private =
    solver->row_refactor_last_separator_flop_ordered_private;
  solver->stats.row_refactor_separator_flop_ordered_private_run_count =
    (int64_t)solver->row_refactor_separator_flop_ordered_private_run_count;
  solver->stats.row_refactor_last_separator_flop_components =
    (int64_t)solver->row_refactor_last_separator_flop_components;
  solver->stats.row_refactor_separator_flop_component_count =
    (int64_t)solver->row_refactor_separator_flop_component_count;
  solver->stats.row_refactor_last_separator_flop_private_groups =
    (int64_t)solver->row_refactor_last_separator_flop_private_groups;
  solver->stats.row_refactor_last_separator_flop_pipeline_groups =
    (int64_t)solver->row_refactor_last_separator_flop_pipeline_groups;
  solver->stats.row_refactor_last_separator_flop_closure_groups =
    (int64_t)solver->row_refactor_last_separator_flop_closure_groups;
  solver->stats.row_refactor_last_separator_flop_private_threads =
    (int64_t)solver->row_refactor_last_separator_flop_private_threads;
  solver->stats.row_refactor_last_separator_flop_private_min_groups =
    (int64_t)solver->row_refactor_last_separator_flop_private_min_groups;
  solver->stats.row_refactor_last_separator_flop_private_max_groups =
    (int64_t)solver->row_refactor_last_separator_flop_private_max_groups;
  solver->stats.row_refactor_last_separator_flop_private_min_work =
    solver->row_refactor_last_separator_flop_private_min_work;
  solver->stats.row_refactor_last_separator_flop_private_max_work =
    solver->row_refactor_last_separator_flop_private_max_work;
  solver->stats.row_refactor_separator_flop_private_group_count =
    (int64_t)solver->row_refactor_separator_flop_private_group_count;
  solver->stats.row_refactor_separator_flop_pipeline_group_count =
    (int64_t)solver->row_refactor_separator_flop_pipeline_group_count;
  solver->stats.row_refactor_separator_flop_closure_group_count =
    (int64_t)solver->row_refactor_separator_flop_closure_group_count;
  solver->stats.kls_tail_last_mapped_columns =
    (int64_t)solver->kls_tail_last_mapped_columns;
  solver->stats.kls_tail_mapped_column_count =
    (int64_t)solver->kls_tail_mapped_column_count;
  solver->stats.kls_first_last_row_uplooking_columns =
    (int64_t)solver->kls_first_last_row_uplooking_columns;
  solver->stats.kls_first_row_uplooking_column_count =
    (int64_t)solver->kls_first_row_uplooking_column_count;
  solver->stats.kls_first_last_dominant_btf_pipeline =
    solver->kls_first_last_dominant_btf_pipeline;
  solver->stats.kls_first_dominant_btf_pipeline_count =
    (int64_t)solver->kls_first_dominant_btf_pipeline_count;
  solver->stats.kls_first_last_dominant_btf_pipeline_block =
    solver->kls_first_last_dominant_btf_pipeline_block == KLS_KLU_EMPTY
      ? -1
      : (int64_t)solver->kls_first_last_dominant_btf_pipeline_block;
  solver->stats.kls_first_last_dominant_btf_pipeline_rows =
    (int64_t)solver->kls_first_last_dominant_btf_pipeline_rows;
  solver->stats.kls_first_last_dominant_btf_pipeline_has_separator =
    solver->kls_first_last_dominant_btf_pipeline_has_separator;
  solver->stats.kls_first_dominant_btf_pipeline_without_separator_count =
    (int64_t)solver
      ->kls_first_dominant_btf_pipeline_without_separator_count;
  solver->stats.kls_first_last_row_refactor_seeded_rows =
    (int64_t)solver->kls_first_last_row_refactor_seeded_rows;
  solver->stats.kls_first_row_refactor_seeded_row_count =
    (int64_t)solver->kls_first_row_refactor_seeded_row_count;
  solver->stats.kls_first_last_row_pipeline =
    solver->kls_first_last_row_pipeline;
  solver->stats.kls_first_row_pipeline_run_count =
    (int64_t)solver->kls_first_row_pipeline_run_count;
  solver->stats.kls_first_last_row_pipeline_rows =
    (int64_t)solver->kls_first_last_row_pipeline_rows;
  solver->stats.kls_first_last_row_pipeline_threads =
    (int64_t)solver->kls_first_last_row_pipeline_threads;
  solver->stats.kls_first_last_row_pipeline_partial =
    solver->kls_first_last_row_pipeline_partial;
  solver->stats.kls_first_row_pipeline_partial_run_count =
    (int64_t)solver->kls_first_row_pipeline_partial_run_count;
  solver->stats.kls_first_last_row_pipeline_partial_rows =
    (int64_t)solver->kls_first_last_row_pipeline_partial_rows;
  solver->stats.kls_first_last_row_pipeline_partial_threads =
    (int64_t)solver->kls_first_last_row_pipeline_partial_threads;
  solver->stats.kls_first_last_row_pipeline_pivot_tail =
    solver->kls_first_last_row_pipeline_pivot_tail;
  solver->stats.kls_first_row_pipeline_pivot_tail_run_count =
    (int64_t)solver->kls_first_row_pipeline_pivot_tail_run_count;
  solver->stats.kls_first_last_row_pipeline_pivot_tail_rows =
    (int64_t)solver->kls_first_last_row_pipeline_pivot_tail_rows;
  solver->stats.kls_first_last_row_pipeline_pivot_restarts =
    (int64_t)solver->kls_first_last_row_pipeline_pivot_restarts;
  solver->stats.kls_first_row_pipeline_pivot_restart_count =
    (int64_t)solver->kls_first_row_pipeline_pivot_restart_count;
  solver->stats.kls_first_last_row_pipeline_pivot_serial_rows =
    (int64_t)solver->kls_first_last_row_pipeline_pivot_serial_rows;
  solver->stats.kls_first_last_row_pipeline_prefix_panel_rebuild =
    solver->kls_first_last_row_pipeline_prefix_panel_rebuild;
  solver->stats.kls_first_row_pipeline_prefix_panel_rebuild_count =
    (int64_t)solver->kls_first_row_pipeline_prefix_panel_rebuild_count;
  solver->stats.kls_first_last_row_pipeline_prefix_panel_rebuild_rows =
    (int64_t)solver->kls_first_last_row_pipeline_prefix_panel_rebuild_rows;
  solver->stats.kls_first_active_rank_pivot_reset_count =
    (int64_t)solver->kls_first_active_rank_pivot_reset_count;
  solver->stats.kls_first_active_rank_pivot_reset_rows =
    (int64_t)solver->kls_first_active_rank_pivot_reset_rows;
  solver->stats.kls_first_active_rank_pivot_panel_rebuild_count =
    (int64_t)solver->kls_first_active_rank_pivot_panel_rebuild_count;
  solver->stats.kls_first_active_rank_pivot_panel_rebuild_rows =
    (int64_t)solver->kls_first_active_rank_pivot_panel_rebuild_rows;
  solver->stats.kls_first_row_panel_cache_build_count =
    (int64_t)solver->kls_first_row_panel_cache_build_count;
  solver->stats.kls_first_row_panel_cache_build_panels =
    (int64_t)solver->kls_first_row_panel_cache_build_panels;
  solver->stats.kls_first_row_panel_cache_build_entries =
    (int64_t)solver->kls_first_row_panel_cache_build_entries;
  solver->stats.kls_first_row_panel_cache_append_count =
    (int64_t)solver->kls_first_row_panel_cache_append_count;
  solver->stats.kls_first_row_panel_cache_append_panels =
    (int64_t)solver->kls_first_row_panel_cache_append_panels;
  solver->stats.kls_first_row_panel_cache_append_entries =
    (int64_t)solver->kls_first_row_panel_cache_append_entries;
  solver->stats.kls_first_last_row_supernode_update =
    solver->kls_first_last_row_supernode_update;
  solver->stats.kls_first_row_supernode_update_run_count =
    (int64_t)solver->kls_first_row_supernode_update_run_count;
  solver->stats.kls_first_row_supernode_update_groups =
    (int64_t)solver->kls_first_row_supernode_update_groups;
  solver->stats.kls_first_row_supernode_update_rows =
    (int64_t)solver->kls_first_row_supernode_update_rows;
  solver->stats.kls_first_last_row_supernode_update_groups =
    (int64_t)solver->kls_first_last_row_supernode_update_groups;
  solver->stats.kls_first_last_row_supernode_update_rows =
    (int64_t)solver->kls_first_last_row_supernode_update_rows;
  solver->stats.kls_first_last_row_supernode_panel_update =
    solver->kls_first_last_row_supernode_panel_update;
  solver->stats.kls_first_row_supernode_panel_update_run_count =
    (int64_t)solver->kls_first_row_supernode_panel_update_run_count;
  solver->stats.kls_first_row_supernode_panel_update_groups =
    (int64_t)solver->kls_first_row_supernode_panel_update_groups;
  solver->stats.kls_first_row_supernode_panel_update_rows =
    (int64_t)solver->kls_first_row_supernode_panel_update_rows;
  solver->stats.kls_first_last_row_supernode_panel_update_groups =
    (int64_t)solver->kls_first_last_row_supernode_panel_update_groups;
  solver->stats.kls_first_last_row_supernode_panel_update_rows =
    (int64_t)solver->kls_first_last_row_supernode_panel_update_rows;
  solver->stats.kls_first_last_dynamic_column_pivots =
    (int64_t)solver->kls_first_last_dynamic_column_pivots;
  solver->stats.kls_first_dynamic_column_pivot_count =
    (int64_t)solver->kls_first_dynamic_column_pivot_count;
  solver->stats.kls_first_last_separator_dynamic_column_pivots =
    (int64_t)solver->kls_first_last_separator_dynamic_column_pivots;
  solver->stats.kls_first_separator_dynamic_column_pivot_count =
    (int64_t)solver->kls_first_separator_dynamic_column_pivot_count;
  solver->stats.kls_first_last_separator_dynamic_column_fallbacks =
    (int64_t)solver->kls_first_last_separator_dynamic_column_fallbacks;
  solver->stats.kls_first_separator_dynamic_column_fallback_count =
    (int64_t)solver->kls_first_separator_dynamic_column_fallback_count;
  const UF_long separator_extent_pivots =
    solver->kls_first_last_separator_extent_dynamic_column_pivots;
  solver->stats.kls_first_last_separator_extent_dynamic_column_pivots =
    (int64_t)separator_extent_pivots;
  solver->stats.kls_first_separator_extent_dynamic_column_pivot_count =
    (int64_t)solver->kls_first_separator_extent_dynamic_column_pivot_count;
  solver->stats.kls_first_last_parallel_btf_blocks =
    (int64_t)solver->kls_first_last_parallel_btf_blocks;
  solver->stats.kls_first_parallel_btf_block_count =
    (int64_t)solver->kls_first_parallel_btf_block_count;
  solver->stats.kls_first_last_separator_dynamic_column_rejects =
    (int64_t)solver->kls_first_last_separator_dynamic_column_rejects;
  solver->stats.kls_first_separator_dynamic_column_reject_count =
    (int64_t)solver->kls_first_separator_dynamic_column_reject_count;
  solver->stats.kls_first_last_separator_queue =
    solver->kls_first_last_separator_queue;
  solver->stats.kls_first_separator_queue_run_count =
    (int64_t)solver->kls_first_separator_queue_run_count;
  solver->stats.kls_first_last_separator_queue_private_components =
    (int64_t)solver->kls_first_last_separator_queue_private_components;
  solver->stats.kls_first_last_separator_queue_pipeline_components =
    (int64_t)solver->kls_first_last_separator_queue_pipeline_components;
  solver->stats.kls_first_last_separator_queue_private_rows =
    (int64_t)solver->kls_first_last_separator_queue_private_rows;
  solver->stats.kls_first_last_separator_queue_pipeline_rows =
    (int64_t)solver->kls_first_last_separator_queue_pipeline_rows;
  solver->stats.kls_first_last_separator_queue_nonempty_threads =
    (int64_t)solver->kls_first_last_separator_queue_nonempty_threads;
  solver->stats.kls_first_last_separator_queue_max_thread_rows =
    (int64_t)solver->kls_first_last_separator_queue_max_thread_rows;
  solver->stats.kls_first_last_separator_queue_min_thread_work =
    solver->kls_first_last_separator_queue_min_thread_work;
  solver->stats.kls_first_last_separator_queue_max_thread_work =
    solver->kls_first_last_separator_queue_max_thread_work;
  solver->stats.kls_first_last_separator_queue_partitioned =
    solver->kls_first_last_separator_queue_partitioned;
  solver->stats.kls_first_separator_queue_partitioned_count =
    (int64_t)solver->kls_first_separator_queue_partitioned_count;
  solver->stats.kls_first_last_separator_queue_split_components =
    (int64_t)solver->kls_first_last_separator_queue_split_components;
  solver->stats.kls_first_last_separator_queue_executed =
    solver->kls_first_last_separator_queue_executed;
  solver->stats.kls_first_separator_queue_executed_run_count =
    (int64_t)solver->kls_first_separator_queue_executed_run_count;
  solver->stats.kls_first_last_separator_queue_executed_private_rows =
    (int64_t)solver->kls_first_last_separator_queue_executed_private_rows;
  solver->stats.kls_first_last_separator_queue_executed_pipeline_rows =
    (int64_t)solver->kls_first_last_separator_queue_executed_pipeline_rows;
  solver->stats.kls_first_last_separator_queue_parallel_private =
    solver->kls_first_last_separator_queue_parallel_private;
  solver->stats.kls_first_separator_queue_parallel_private_run_count =
    (int64_t)solver->kls_first_separator_queue_parallel_private_run_count;
  solver->stats.kls_first_last_separator_queue_parallel_private_rows =
    (int64_t)solver->kls_first_last_separator_queue_parallel_private_rows;
  solver->stats.kls_first_last_separator_queue_parallel_private_threads =
    (int64_t)solver->kls_first_last_separator_queue_parallel_private_threads;
  solver->stats.kls_first_last_separator_queue_parallel_pipeline =
    solver->kls_first_last_separator_queue_parallel_pipeline;
  solver->stats.kls_first_separator_queue_parallel_pipeline_run_count =
    (int64_t)solver->kls_first_separator_queue_parallel_pipeline_run_count;
  solver->stats.kls_first_last_separator_queue_parallel_pipeline_rows =
    (int64_t)solver->kls_first_last_separator_queue_parallel_pipeline_rows;
  solver->stats.kls_first_last_separator_queue_parallel_pipeline_threads =
    (int64_t)solver->kls_first_last_separator_queue_parallel_pipeline_threads;
  solver->stats.kls_first_last_separator_queue_pipeline_partial =
    solver->kls_first_last_separator_queue_pipeline_partial;
  solver->stats.kls_first_separator_queue_pipeline_partial_run_count =
    (int64_t)solver->kls_first_separator_queue_pipeline_partial_run_count;
  solver->stats.kls_first_last_separator_queue_pipeline_partial_rows =
    (int64_t)solver->kls_first_last_separator_queue_pipeline_partial_rows;
  solver->stats.kls_first_last_separator_queue_pipeline_partial_threads =
    (int64_t)solver->kls_first_last_separator_queue_pipeline_partial_threads;
  solver->stats.kls_first_last_separator_queue_pipeline_wait_partial =
    solver->kls_first_last_separator_queue_pipeline_wait_partial;
  solver->stats.kls_first_separator_queue_pipeline_wait_partial_run_count =
    (int64_t)
      solver->kls_first_separator_queue_pipeline_wait_partial_run_count;
  solver->stats.kls_first_last_separator_queue_pipeline_wait_partial_rows =
    (int64_t)solver->kls_first_last_separator_queue_pipeline_wait_partial_rows;
  solver->stats.kls_first_last_separator_queue_pipeline_wait_partial_deps =
    (int64_t)solver->kls_first_last_separator_queue_pipeline_wait_partial_deps;
  solver->stats.kls_first_last_separator_queue_pipeline_supernode_update =
    solver->kls_first_last_separator_queue_pipeline_supernode_update;
  solver->stats.kls_first_separator_queue_pipeline_supernode_update_run_count =
    (int64_t)
      solver->kls_first_separator_queue_pipeline_supernode_update_run_count;
  solver->stats.kls_first_last_separator_queue_pipeline_supernode_update_groups =
    (int64_t)
      solver->kls_first_last_separator_queue_pipeline_supernode_update_groups;
  solver->stats.kls_first_last_separator_queue_pipeline_supernode_update_rows =
    (int64_t)
      solver->kls_first_last_separator_queue_pipeline_supernode_update_rows;
  solver->stats.kls_first_last_separator_queue_pipeline_supernode_panel_update =
    solver->kls_first_last_separator_queue_pipeline_supernode_panel_update;
  solver->stats
    .kls_first_separator_queue_pipeline_supernode_panel_update_run_count =
      (int64_t)solver
        ->kls_first_separator_queue_pipeline_supernode_panel_update_run_count;
  solver->stats
    .kls_first_last_separator_queue_pipeline_supernode_panel_update_groups =
      (int64_t)solver
        ->kls_first_last_separator_queue_pipeline_supernode_panel_update_groups;
  solver->stats
    .kls_first_last_separator_queue_pipeline_supernode_panel_update_rows =
      (int64_t)solver
        ->kls_first_last_separator_queue_pipeline_supernode_panel_update_rows;
  solver->stats.kls_first_last_separator_queue_pipeline_pivot_tail =
    solver->kls_first_last_separator_queue_pipeline_pivot_tail;
  solver->stats.kls_first_separator_queue_pipeline_pivot_tail_run_count =
    (int64_t)
      solver->kls_first_separator_queue_pipeline_pivot_tail_run_count;
  solver->stats.kls_first_last_separator_queue_pipeline_pivot_tail_rows =
    (int64_t)solver->kls_first_last_separator_queue_pipeline_pivot_tail_rows;
  solver->stats.kls_first_last_separator_queue_pipeline_pivot_restarts =
    (int64_t)solver->kls_first_last_separator_queue_pipeline_pivot_restarts;
  solver->stats.kls_first_separator_queue_pipeline_pivot_restart_count =
    (int64_t)solver->kls_first_separator_queue_pipeline_pivot_restart_count;
  solver->stats.kls_first_last_separator_queue_pipeline_pivot_serial_rows =
    (int64_t)solver->kls_first_last_separator_queue_pipeline_pivot_serial_rows;
  solver->stats.kls_first_last_separator_queue_pipeline_prefix_panel_rebuild =
    solver->kls_first_last_separator_queue_pipeline_prefix_panel_rebuild;
  solver->stats
    .kls_first_separator_queue_pipeline_prefix_panel_rebuild_count =
      (int64_t)
        solver->kls_first_separator_queue_pipeline_prefix_panel_rebuild_count;
  solver->stats
    .kls_first_last_separator_queue_pipeline_prefix_panel_rebuild_rows =
      (int64_t)
        solver->kls_first_last_separator_queue_pipeline_prefix_panel_rebuild_rows;
  solver->stats.kls_first_auto_skipped_scaled_single_block =
    solver->kls_first_auto_skipped_scaled_single_block;
  solver->stats.kls_first_auto_skipped_scaled_single_block_count =
    (int64_t)solver->kls_first_auto_skipped_scaled_single_block_count;
  solver->stats.row_refactor_last_done_bitmap =
    solver->row_refactor_last_done_bitmap;
  solver->stats.row_refactor_done_bitmap_run_count =
    (int64_t)solver->row_refactor_done_bitmap_run_count;
  solver->stats.row_refactor_last_prefactor =
    solver->row_refactor_last_prefactor;
  solver->stats.row_refactor_last_prefactor_rows =
    (int64_t)solver->row_refactor_last_prefactor_rows;
  solver->stats.row_refactor_last_prefactor_deps =
    (int64_t)solver->row_refactor_last_prefactor_deps;
  solver->stats.row_refactor_prefactor_run_count =
    (int64_t)solver->row_refactor_prefactor_run_count;
  solver->stats.row_refactor_prefactor_rows =
    (int64_t)solver->row_refactor_prefactor_rows;
  solver->stats.row_refactor_prefactor_deps =
    (int64_t)solver->row_refactor_prefactor_deps;
  solver->stats.row_refactor_last_prefactor_supernode =
    solver->row_refactor_last_prefactor_supernode;
  solver->stats.row_refactor_last_prefactor_supernode_rows =
    (int64_t)solver->row_refactor_last_prefactor_supernode_rows;
  solver->stats.row_refactor_last_prefactor_supernode_deps =
    (int64_t)solver->row_refactor_last_prefactor_supernode_deps;
  solver->stats.row_refactor_prefactor_supernode_run_count =
    (int64_t)solver->row_refactor_prefactor_supernode_run_count;
  solver->stats.row_refactor_prefactor_supernode_rows =
    (int64_t)solver->row_refactor_prefactor_supernode_rows;
  solver->stats.row_refactor_prefactor_supernode_deps =
    (int64_t)solver->row_refactor_prefactor_supernode_deps;
  solver->stats.row_refactor_input_cleanup_rows =
    (int64_t)solver->row_refactor_input_cleanup_rows;
  solver->stats.row_refactor_input_cleanup_entries =
    (int64_t)solver->row_refactor_input_cleanup_entries;
  solver->stats.row_refactor_last_defer_value_scatter =
    solver->row_refactor_last_defer_value_scatter;
  solver->stats.row_refactor_defer_value_scatter_run_count =
    (int64_t)solver->row_refactor_defer_value_scatter_run_count;
  solver->stats.row_refactor_values_dirty =
    solver->row_refactor_values_dirty;
  solver->stats.row_refactor_last_lazy_value_scatter =
    solver->row_refactor_last_lazy_value_scatter;
  solver->stats.row_refactor_lazy_value_scatter_run_count =
    (int64_t)solver->row_refactor_lazy_value_scatter_run_count;
  solver->stats.row_refactor_last_row_solve =
    solver->row_refactor_last_row_solve;
  solver->stats.row_refactor_row_solve_run_count =
    (int64_t)solver->row_refactor_row_solve_run_count;
  solver->stats.row_solve_parallel_run_count =
    (int64_t)solver->row_solve_parallel_run_count;
  solver->stats.row_solve_parallel_l_slice_runs =
    (int64_t)solver->row_solve_parallel_l_slice_runs;
  solver->stats.row_solve_parallel_u_slice_runs =
    (int64_t)solver->row_solve_parallel_u_slice_runs;
  solver->stats.row_solve_parallel_l_sparse_level_runs =
    (int64_t)solver->row_solve_parallel_l_sparse_level_runs;
  solver->stats.row_solve_parallel_u_sparse_level_runs =
    (int64_t)solver->row_solve_parallel_u_sparse_level_runs;
  solver->stats.row_solve_thread_count =
    (int64_t)solver->row_solve_thread_count;
  solver->stats.row_solve_l_thread_max_rect_entries =
    (int64_t)solver->row_solve_l_thread_max_rect_entries;
  solver->stats.row_solve_u_thread_max_rect_entries =
    (int64_t)solver->row_solve_u_thread_max_rect_entries;
  solver->stats.row_solve_partition_ready =
    solver->row_solve_partition_ready;
  solver->stats.row_solve_partition_slices =
    (int64_t)solver->row_solve_partition_slices;
  solver->stats.row_solve_l_sparse_level_count =
    (int64_t)solver->row_solve_l_sparse_level_count;
  solver->stats.row_solve_l_sparse_cluster_levels =
    (int64_t)solver->row_solve_l_sparse_cluster_levels;
  solver->stats.row_solve_l_sparse_level_max_width =
    (int64_t)solver->row_solve_l_sparse_level_max_width;
  solver->stats.row_solve_l_dense_tail_start =
    (int64_t)solver->row_solve_l_dense_tail_start;
  solver->stats.row_solve_l_dense_tail_rows =
    (int64_t)solver->row_solve_l_dense_tail_rows;
  solver->stats.row_solve_l_dense_tail_entries =
    (int64_t)solver->row_solve_l_dense_tail_entries;
  solver->stats.row_solve_l_slice_max_entries =
    (int64_t)solver->row_solve_l_slice_max_entries;
  solver->stats.row_solve_l_segmented_rows =
    (int64_t)solver->row_solve_l_segmented_rows;
  solver->stats.row_solve_l_rect_entries =
    (int64_t)solver->row_solve_l_rect_entries;
  solver->stats.row_solve_l_tri_entries =
    (int64_t)solver->row_solve_l_tri_entries;
  solver->stats.row_solve_u_sparse_level_count =
    (int64_t)solver->row_solve_u_sparse_level_count;
  solver->stats.row_solve_u_sparse_cluster_levels =
    (int64_t)solver->row_solve_u_sparse_cluster_levels;
  solver->stats.row_solve_u_sparse_level_max_width =
    (int64_t)solver->row_solve_u_sparse_level_max_width;
  solver->stats.row_solve_u_dense_tail_start =
    (int64_t)solver->row_solve_u_dense_tail_start;
  solver->stats.row_solve_u_dense_tail_rows =
    (int64_t)solver->row_solve_u_dense_tail_rows;
  solver->stats.row_solve_u_dense_tail_entries =
    (int64_t)solver->row_solve_u_dense_tail_entries;
  solver->stats.row_solve_u_slice_max_entries =
    (int64_t)solver->row_solve_u_slice_max_entries;
  solver->stats.row_solve_u_segmented_rows =
    (int64_t)solver->row_solve_u_segmented_rows;
  solver->stats.row_solve_u_rect_entries =
    (int64_t)solver->row_solve_u_rect_entries;
  solver->stats.row_solve_u_tri_entries =
    (int64_t)solver->row_solve_u_tri_entries;
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
  solver->stats.row_refactor_compact_dense_panel_eligible_count =
    (int64_t)solver->row_refactor_compact_dense_panel_eligible_count;
  solver->stats.row_refactor_compact_dense_panel_eligible_rows =
    (int64_t)solver->row_refactor_compact_dense_panel_eligible_rows;
  solver->stats.row_refactor_compact_dense_panel_update_work =
    solver->row_refactor_compact_dense_panel_update_work;
  solver->stats.row_refactor_compact_dense_panel_entries =
    solver->row_refactor_compact_dense_panel_entries;
  solver->stats.row_refactor_dense_producer_run_count =
    (int64_t)solver->row_refactor_dense_producer_run_count;
  solver->stats.row_refactor_dense_producer_run_rows =
    (int64_t)solver->row_refactor_dense_producer_run_rows;
  solver->stats.row_refactor_dense_producer_run_dep_rows =
    (int64_t)solver->row_refactor_dense_producer_run_dep_rows;
  solver->stats.row_refactor_dense_producer_run_max_per_row =
    (int64_t)solver->row_refactor_dense_producer_run_max_per_row;
  solver->stats.row_refactor_dense_producer_full_suffix_run_count =
    (int64_t)solver->row_refactor_dense_producer_full_suffix_run_count;
  solver->stats.row_refactor_dense_producer_full_suffix_rows =
    (int64_t)solver->row_refactor_dense_producer_full_suffix_rows;
  solver->stats.row_refactor_dense_producer_multi_run_rows =
    (int64_t)solver->row_refactor_dense_producer_multi_run_rows;
  solver->stats.row_refactor_dense_producer_fragmented_rows =
    (int64_t)solver->row_refactor_dense_producer_fragmented_rows;
  solver->stats.row_refactor_dense_producer_target_count =
    (int64_t)solver->row_refactor_dense_producer_target_count;
  solver->stats.row_refactor_dense_producer_target_none_count =
    (int64_t)solver->row_refactor_dense_producer_target_none_count;
  solver->stats.row_refactor_dense_producer_target_external_count =
    (int64_t)solver->row_refactor_dense_producer_target_external_count;
  solver->stats.row_refactor_dense_producer_target_dense_count =
    (int64_t)solver->row_refactor_dense_producer_target_dense_count;
  solver->stats.row_refactor_dense_producer_target_pivot_count =
    (int64_t)solver->row_refactor_dense_producer_target_pivot_count;
  solver->stats.row_refactor_dense_producer_target_trailing_count =
    (int64_t)solver->row_refactor_dense_producer_target_trailing_count;
  solver->stats.row_refactor_compact_dense_panel_persistent_groups =
    (int64_t)solver->row_refactor_compact_dense_panel_persistent_groups;
  solver->stats.row_refactor_compact_dense_panel_persistent_entries =
    (int64_t)solver->row_refactor_compact_dense_panel_persistent_entries;
  solver->stats.row_refactor_last_compact_dense_panel_persistent =
    solver->row_refactor_last_compact_dense_panel_persistent;
  solver->stats.row_refactor_compact_dense_panel_persistent_run_count =
    (int64_t)solver->row_refactor_compact_dense_panel_persistent_run_count;
  solver->stats.row_refactor_last_compact_dense_panel_blocked =
    solver->row_refactor_last_compact_dense_panel_blocked;
  solver->stats.row_refactor_compact_dense_panel_blocked_run_count =
    (int64_t)solver->row_refactor_compact_dense_panel_blocked_run_count;
  solver->stats.row_refactor_compact_dense_panel_blocked_rows =
    (int64_t)solver->row_refactor_compact_dense_panel_blocked_rows;
  solver->stats.row_refactor_compact_dense_panel_blocked_entries =
    (int64_t)solver->row_refactor_compact_dense_panel_blocked_entries;
  solver->stats.row_refactor_native_row_panel_enabled =
    solver->row_refactor_native_row_panel_enabled;
  solver->stats.row_refactor_last_native_row_panel =
    solver->row_refactor_last_native_row_panel;
  solver->stats.row_refactor_native_row_panel_auto_disabled =
    solver->row_refactor_native_row_panel_auto_disabled;
  solver->stats.row_refactor_native_row_panel_count =
    (int64_t)solver->row_refactor_native_row_panel_count;
  solver->stats.row_refactor_native_row_panel_rows =
    (int64_t)solver->row_refactor_native_row_panel_rows;
  solver->stats.row_refactor_native_row_panel_entries =
    (int64_t)solver->row_refactor_native_row_panel_entries;
  solver->stats.row_refactor_native_row_panel_blocked_count =
    (int64_t)solver->row_refactor_native_row_panel_blocked_count;
  solver->stats.row_refactor_native_row_panel_blocked_rows =
    (int64_t)solver->row_refactor_native_row_panel_blocked_rows;
  solver->stats.row_refactor_native_row_panel_blocked_entries =
    (int64_t)solver->row_refactor_native_row_panel_blocked_entries;
  solver->stats.row_refactor_native_row_panel_fallback_count =
    (int64_t)solver->row_refactor_native_row_panel_fallback_count;
  solver->stats.row_refactor_native_row_panel_checked_reject_count =
    (int64_t)solver->row_refactor_native_row_panel_checked_reject_count;
  solver->stats.row_refactor_native_row_panel_auto_disable_count =
    (int64_t)solver->row_refactor_native_row_panel_auto_disable_count;
  solver->stats.row_refactor_last_compact_dense_panel_direct_input_rows =
    (int64_t)solver->row_refactor_last_compact_dense_panel_direct_input_rows;
  solver->stats.row_refactor_compact_dense_panel_direct_input_rows =
    (int64_t)solver->row_refactor_compact_dense_panel_direct_input_rows;
  solver->stats.row_refactor_last_compact_panel_solve_values =
    (int64_t)solver->row_refactor_last_compact_panel_solve_values;
  solver->stats.row_refactor_compact_panel_solve_values =
    (int64_t)solver->row_refactor_compact_panel_solve_values;
  solver->stats.row_refactor_last_compact_panel_group_solve_rows =
    (int64_t)solver->row_refactor_last_compact_panel_group_solve_rows;
  solver->stats.row_refactor_compact_panel_group_solve_rows =
    (int64_t)solver->row_refactor_compact_panel_group_solve_rows;
  solver->stats.row_refactor_last_compact_panel_group_solve_entries =
    (int64_t)solver->row_refactor_last_compact_panel_group_solve_entries;
  solver->stats.row_refactor_compact_panel_group_solve_entries =
    (int64_t)solver->row_refactor_compact_panel_group_solve_entries;
  solver->stats.row_refactor_last_compact_panel_scalar_update_rows =
    (int64_t)solver->row_refactor_last_compact_panel_scalar_update_rows;
  solver->stats.row_refactor_compact_panel_scalar_update_rows =
    (int64_t)solver->row_refactor_compact_panel_scalar_update_rows;
  solver->stats.row_refactor_last_compact_panel_scalar_update_entries =
    (int64_t)solver->row_refactor_last_compact_panel_scalar_update_entries;
  solver->stats.row_refactor_compact_panel_scalar_update_entries =
    (int64_t)solver->row_refactor_compact_panel_scalar_update_entries;
  solver->stats.row_refactor_last_dense_segment_direct_input_rows =
    (int64_t)solver->row_refactor_last_dense_segment_direct_input_rows;
  solver->stats.row_refactor_dense_segment_direct_input_rows =
    (int64_t)solver->row_refactor_dense_segment_direct_input_rows;
  solver->stats.row_refactor_last_sparse_segment_direct_input_rows =
    (int64_t)solver->row_refactor_last_sparse_segment_direct_input_rows;
  solver->stats.row_refactor_sparse_segment_direct_input_rows =
    (int64_t)solver->row_refactor_sparse_segment_direct_input_rows;
  solver->stats.row_refactor_last_batch_direct_input_rows =
    (int64_t)solver->row_refactor_last_batch_direct_input_rows;
  solver->stats.row_refactor_batch_direct_input_rows =
    (int64_t)solver->row_refactor_batch_direct_input_rows;
  solver->stats.row_refactor_segment_input_target_rows =
    (int64_t)solver->row_refactor_segment_input_target_rows;
  solver->stats.row_refactor_segment_input_target_entries =
    (int64_t)solver->row_refactor_segment_input_target_entries;
  solver->stats.row_refactor_last_segment_target_input_rows =
    (int64_t)solver->row_refactor_last_segment_target_input_rows;
  solver->stats.row_refactor_segment_target_input_rows =
    (int64_t)solver->row_refactor_segment_target_input_rows;
  solver->stats.row_refactor_segment_input_cleanup_rows =
    (int64_t)solver->row_refactor_segment_input_cleanup_rows;
  solver->stats.row_refactor_segment_input_cleanup_entries =
    (int64_t)solver->row_refactor_segment_input_cleanup_entries;
  solver->stats.row_refactor_last_segment_target_cleanup_rows =
    (int64_t)solver->row_refactor_last_segment_target_cleanup_rows;
  solver->stats.row_refactor_segment_target_cleanup_rows =
    (int64_t)solver->row_refactor_segment_target_cleanup_rows;
  solver->stats.row_refactor_last_segment_target_cleanup_entries =
    (int64_t)solver->row_refactor_last_segment_target_cleanup_entries;
  solver->stats.row_refactor_segment_target_cleanup_entries =
    (int64_t)solver->row_refactor_segment_target_cleanup_entries;
  solver->stats.row_refactor_last_compact_supernode_update =
    solver->row_refactor_last_compact_supernode_update;
  solver->stats.row_refactor_compact_supernode_update_count =
    (int64_t)solver->row_refactor_compact_supernode_update_count;
  solver->stats.row_refactor_compact_supernode_update_rows =
    (int64_t)solver->row_refactor_compact_supernode_update_rows;
  solver->stats.row_refactor_compact_supernode_update_entries =
    (int64_t)solver->row_refactor_compact_supernode_update_entries;
  solver->stats.row_refactor_last_compact_supernode_partial_update =
    solver->row_refactor_last_compact_supernode_partial_update;
  solver->stats.row_refactor_compact_supernode_partial_update_count =
    (int64_t)solver->row_refactor_compact_supernode_partial_update_count;
  solver->stats.row_refactor_compact_supernode_partial_update_rows =
    (int64_t)solver->row_refactor_compact_supernode_partial_update_rows;
  solver->stats.row_refactor_compact_supernode_partial_update_entries =
    (int64_t)solver->row_refactor_compact_supernode_partial_update_entries;
  solver->stats.row_refactor_last_compact_supernode_gemv =
    solver->row_refactor_last_compact_supernode_gemv;
  solver->stats.row_refactor_compact_supernode_gemv_count =
    (int64_t)solver->row_refactor_compact_supernode_gemv_count;
  solver->stats.row_refactor_compact_supernode_gemv_rows =
    (int64_t)solver->row_refactor_compact_supernode_gemv_rows;
  solver->stats.row_refactor_compact_supernode_gemv_entries =
    (int64_t)solver->row_refactor_compact_supernode_gemv_entries;
  solver->stats.row_refactor_last_compact_supernode_trsv =
    solver->row_refactor_last_compact_supernode_trsv;
  solver->stats.row_refactor_compact_supernode_trsv_count =
    (int64_t)solver->row_refactor_compact_supernode_trsv_count;
  solver->stats.row_refactor_compact_supernode_trsv_rows =
    (int64_t)solver->row_refactor_compact_supernode_trsv_rows;
  solver->stats.row_refactor_compact_supernode_trsv_entries =
    (int64_t)solver->row_refactor_compact_supernode_trsv_entries;
  solver->stats.row_refactor_last_compact_supernode_batch =
    solver->row_refactor_last_compact_supernode_batch;
  solver->stats.row_refactor_compact_supernode_batch_count =
    (int64_t)solver->row_refactor_compact_supernode_batch_count;
  solver->stats.row_refactor_compact_supernode_batch_rows =
    (int64_t)solver->row_refactor_compact_supernode_batch_rows;
  solver->stats.row_refactor_compact_supernode_batch_dep_rows =
    (int64_t)solver->row_refactor_compact_supernode_batch_dep_rows;
  solver->stats.row_refactor_compact_supernode_batch_entries =
    (int64_t)solver->row_refactor_compact_supernode_batch_entries;
  solver->stats.row_refactor_compact_supernode_batch_pattern_count =
    (int64_t)solver->row_refactor_compact_supernode_batch_pattern_count;
  solver->stats.row_refactor_compact_supernode_batch_pattern_rows =
    (int64_t)solver->row_refactor_compact_supernode_batch_pattern_rows;
  solver->stats.row_refactor_compact_supernode_batch_candidate_count =
    (int64_t)solver->row_refactor_compact_supernode_batch_candidate_count;
  solver->stats.row_refactor_compact_supernode_batch_candidate_rows =
    (int64_t)solver->row_refactor_compact_supernode_batch_candidate_rows;
  solver->stats.row_refactor_compact_supernode_batch_candidate_dep_rows =
    (int64_t)solver->row_refactor_compact_supernode_batch_candidate_dep_rows;
  solver->stats.row_refactor_compact_supernode_batch_rejected_work_count =
    (int64_t)solver->row_refactor_compact_supernode_batch_rejected_work_count;
  solver->stats.refactor_dependency_cluster_levels =
    (int64_t)solver->refactor_cluster_level_count;
  solver->stats.refactor_dependency_pipeline_columns =
    (int64_t)solver->refactor_pipeline_column_count;
  solver->stats.refactor_dependency_work =
    solver->refactor_dependency_work;
  solver->stats.refactor_dependency_pipeline_work =
    solver->refactor_pipeline_work;

}
