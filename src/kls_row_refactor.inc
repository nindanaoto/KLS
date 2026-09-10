/* Row-oriented refactor kernels, compact row solves, and row execution policy.
 * Included by kls_numeric_engines.inc in the original definition order.
 */

KLS_DEFINE_ENV_FLAG(kls_row_refactor_env_enabled,
                    "KLS_ENABLE_ROW_REFACTOR")

static int kls_row_refactor_env_disabled(void) {
  const char *value = getenv("KLS_ENABLE_ROW_REFACTOR");
  return value != NULL && value[0] == '0' && value[1] == '\0';
}

/* Predicted-pattern numerics retain their column refactor engines; the
   ordinary row-policy gates below exclude predicted numeric factors. */

/* The row executor's setup is not worthwhile for ordinary sparse factors,
   but the column EGraph's irregular scatter stream becomes dominant for the
   very large dependency graphs.  These cases have enough repeated numeric
   work to amortize construction; retain the usual measured acceptance trial
   rather than forcing the row engine outright. */
static int kls_auto_row_refactor_high_work_candidate(
  const kls_solver *solver) {
  return solver != NULL && solver->options.threads > 1 &&
         /* The row schedule retains O(n) metadata.  On the measured
            cache-domain workload it pays back on the mid-sized, high-fill
            ASIC factors, but its preparation pushes multi-million-row
            factors beyond the per-case wall-time bound before the numeric
            executor can amortize it. */
         solver->n >= 300000u && solver->n <= 800000u &&
         solver->refactor_dependency_work >= 1.0e8;
}

KLS_DEFINE_ENV_FLAG(kls_checked_row_refactor_env_enabled,
                    "KLS_ENABLE_CHECKED_ROW_REFACTOR")

enum {
  KLS_NATIVE_ROW_PANEL_ENV_OFF = -1,
  KLS_NATIVE_ROW_PANEL_ENV_AUTO = 0
};

static int kls_native_row_panel_structural_auto_allows(
  const kls_solver *solver) {
  return solver != NULL &&
         solver->row_refactor_compact_dense_panel_persistent_groups > 0u &&
         solver->row_refactor_compact_dense_panel_eligible_rows > 1u &&
         solver->row_refactor_compact_dense_panel_update_work >=
           KLS_ROW_REFACTOR_NATIVE_PANEL_AUTO_MIN_WORK &&
         solver->row_refactor_compact_dense_panel_update_work >=
           KLS_ROW_REFACTOR_COMPACT_PANEL_MIN_WORK_PER_ENTRY *
             solver->row_refactor_compact_dense_panel_entries;
}

static int kls_native_row_panel_effective_state(const kls_solver *solver) {
  return solver != NULL &&
         (solver->row_refactor_auto_native_row_panel ||
          kls_native_row_panel_structural_auto_allows(solver))
           ? KLS_NATIVE_ROW_PANEL_ENV_AUTO
           : KLS_NATIVE_ROW_PANEL_ENV_OFF;
}

static int kls_native_row_panel_enabled_for_run(const kls_solver *solver,
                                                int env_state,
                                                int check_pivots);
static void kls_maybe_disable_native_row_panel_auto(
  kls_solver *solver,
  int check_pivots,
  int env_state,
  UF_long native_count_before,
  UF_long compact_update_count_before,
  UF_long partial_pipeline_count_before);

static int kls_egraph_cached_supernode_update_allows(
  UF_long run_rows,
  UF_long trailing_len) {
  if (run_rows < KLS_EGRAPH_CACHED_SUPERNODE_MIN_ROWS) {
    return 0;
  }
  const double work =
    0.5 * (double)run_rows * (double)(run_rows - 1u) +
    (double)run_rows * (double)trailing_len;
  const double copied_entries =
    (double)run_rows + (double)trailing_len;
  return work >= KLS_EGRAPH_CACHED_SUPERNODE_MIN_WORK &&
         work >= KLS_EGRAPH_CACHED_SUPERNODE_MIN_WORK_PER_ENTRY *
                   copied_entries;
}


static int kls_pipeline_refactor_fast_factor_repair_is_risky(
  const kls_solver *solver) {
  if (solver == NULL ||
      solver->refactor_dependency_work <
        KLS_FAST_FACTOR_PIPELINE_REFACTOR_MIN_WORK ||
      solver->refactor_pipeline_work <= 0.0 ||
      solver->refactor_level_ptr == NULL ||
      solver->options.threads <= 1) {
    return 0;
  }
  if (solver->common.noffdiag == 0u) {
    return 0;
  }
  return solver->refactor_pipeline_work >=
         KLS_FAST_FACTOR_PIPELINE_REFACTOR_MIN_SHARE *
           solver->refactor_dependency_work;
}


static double kls_row_refactor_total_group_work(const kls_solver *solver) {
  if (solver == NULL || solver->row_refactor_group_count == 0u) {
    return 0.0;
  }
  double work = 0.0;
  for (UF_long group = 0; group < solver->row_refactor_group_count; ++group) {
    work += kls_row_refactor_group_work(solver, group);
  }
  return work;
}

static int kls_auto_row_refactor_cost_allows(const kls_solver *solver) {
  if (solver == NULL || !solver->row_refactor_auto_enabled) {
    return 0;
  }

  if (solver->refactor_dependency_work <= 0.0) {
    return 1;
  }
  if (solver->row_refactor_pattern_n != solver->n ||
      solver->row_refactor_group_count == 0u) {
    return 1;
  }

  const double row_work = kls_row_refactor_total_group_work(solver);
  if (row_work <= 0.0) {
    return 1;
  }
  if (row_work <= 1.05 * solver->refactor_dependency_work) {
    /* The estimates count arithmetic but not the column engine's scatter and
       synchronization streams.  Treat a five-percent work band as genuinely
       undecided and let the paired runtime audit arbitrate it.  This also
       permits the first public update to be the row sample, avoiding a cold
       column pass whose only purpose would be consultation setup. */
    return 1;
  }
  return 0;
}

static int kls_auto_row_refactor_should_run(const kls_solver *solver) {
  if (!kls_row_refactor_env_enabled() &&
      !kls_auto_row_refactor_high_work_candidate(solver) &&
      (solver == NULL || solver->n > 64u)) {
    return 0;
  }
  if (solver != NULL && solver->numeric_is_predicted) {
    return 0;
  }
  return solver != NULL &&
         !kls_row_refactor_env_disabled() &&
         solver->row_refactor_auto_enabled &&
         solver->row_refactor_values_ready &&
         kls_auto_row_refactor_cost_allows(solver);
}

/* Measured per-numeric adoption of the row-refactor engine. Analytic
   models mispredict the cooperative engine (memchip: modeled hopeless,
   measured faster with half the solve time), so alternate real
   refactors between the column and row engines, pair each with its
   following solve, and decide on the measured refactor+solve pair —
   the quantity a SPICE cycle actually pays. */
static int kls_row_refactor_acceptance_structurally_ready(
  const kls_solver *solver) {
  return (kls_row_refactor_env_enabled() ||
          kls_auto_row_refactor_high_work_candidate(solver) ||
          (solver != NULL &&
           (solver->n <= 64u))) &&
         solver != NULL &&
         (!solver->numeric_is_predicted) &&
         !kls_row_refactor_env_disabled() &&
         solver->row_refactor_auto_enabled &&
         solver->row_refactor_pattern_n == solver->n &&
         solver->options.threads > 1;
}

static int kls_row_refactor_acceptance_wants_row(kls_solver *solver) {
  solver->row_trial_deadline = 0.0;

  if (kls_row_refactor_env_enabled()) {
    return 1;
  }
  if (!kls_row_refactor_acceptance_structurally_ready(solver)) {
    return 0;
  }
  if (solver->row_accept_decision != 0) {
    if (solver->row_accept_decision > 0 &&
        solver->row_accept_first_consult) {
      const double base = solver->adaptive_refactor_seconds > 0.0
        ? solver->adaptive_refactor_seconds
        : 0.5 * solver->stats.factor_seconds;
      solver->row_trial_deadline = base > 0.0 ? 4.0 * base + 0.010 : 0.0;
    }
    if (solver->row_accept_decision > 0 &&
        (solver->row_reaudit_state == 1 ||
         solver->row_reaudit_state == 4)) {
      /* Two forced column cycles refresh the loser's refactor and solve
         floors.  One cold column pass is not representative after eight
         row cycles have displaced its cache and scheduler state. */
      solver->row_reaudit_state =
        solver->row_reaudit_state == 1 ? 2 : 5;
      return 0;
    }
    if (solver->row_accept_decision > 0 &&
        (solver->row_reaudit_state == 3 ||
         solver->row_reaudit_state == 6)) {
      /* A refactor-only burst may postpone the solve.  Keep the published
         column representation until its paired solve arrives instead of
         silently switching the sample back to the row engine. */
      return 0;
    }
    return solver->row_accept_decision > 0;
  }
  const int row_n = solver->row_accept_ref_samples[1];
  const int col_n = solver->row_accept_ref_samples[0];
  int want;
  if (row_n == col_n) {
    want = kls_auto_row_refactor_cost_allows(solver) ? 1 : 0;
  } else {
    want = row_n < col_n;
  }
  if (want) {
    /* Bound the trial: a pathological row schedule (mac_econ-class runs
       100x the column time at equal flops) must cost at most a few
       column refactors before the verdict settles on the column side. */
    const double base = solver->adaptive_refactor_seconds > 0.0
      ? solver->adaptive_refactor_seconds
      : 0.5 * solver->stats.factor_seconds;
    solver->row_trial_deadline = base > 0.0 ? 4.0 * base + 0.010 : 0.0;
  }
  return want;
}

static void kls_row_refactor_acceptance_try_decide(kls_solver *solver) {
  if (solver->row_accept_decision != 0) {
    return;
  }
  const int col_refs = solver->row_accept_ref_samples[0];
  const int row_refs = solver->row_accept_ref_samples[1];
  const int col_solves = solver->row_accept_solve_samples[0];
  const int row_solves = solver->row_accept_solve_samples[1];
  if (col_refs >= 1 && row_refs >= 1 &&
      col_solves >= 1 && row_solves >= 1 &&
      solver->row_accept_ref_min[1] + solver->row_accept_solve_min[1] >
        1.25 * (solver->row_accept_ref_min[0] +
                solver->row_accept_solve_min[0])) {
    /* A complete pair that is already 25% slower cannot plausibly be timing
       noise.  Reject that arm after one solve per side; positive adoption
       still requires two samples per side so cold-column optimism cannot
       promote row. */
    solver->row_accept_decision = -1;
  } else if (col_refs >= 2 && row_refs >= 2 &&
             col_solves >= 2 && row_solves >= 2) {
    /* Decide on per-side minima: the sides' engines warm at different
       rates (a cold egraph sample can run 4x its steady rate), and a
       mean over two samples freezes that noise into the verdict. The
       minimum is each side's demonstrated steady rate. */
    const double col_pair =
      solver->row_accept_ref_min[0] + solver->row_accept_solve_min[0];
    const double row_pair =
      solver->row_accept_ref_min[1] + solver->row_accept_solve_min[1];
    solver->row_accept_decision = row_pair < col_pair * 0.98 ? 1 : -1;
    if (solver->row_accept_decision > 0) {
      /* The verdict pairs also measured both solve routes: when the
         published column solve decisively beat the row-value solve,
         publish even where the structural gate would keep row values
         (pre2: row solve 168ms vs column 64ms at +~20ms publish). */
      if (solver->row_accept_solve_min[0] <
          0.8 * solver->row_accept_solve_min[1]) {
        /* A faster column solve is only a candidate: reaching it requires
           publishing the complete row factor back into packed-column
           storage.  Measure that transfer on the next refactor and decide
           on publication-plus-solve instead of solve time alone. */
        solver->row_publish_experiment = 1;
        solver->row_publish_probe_seconds = 0.0;
      }
    }
  } else if (col_refs + row_refs >= 8 && col_refs >= 2 && row_refs >= 2) {
    /* Solves never arrived (pure refactor burst): adopt only on a
       decisive refactor margin, mirroring the pts_ref precedent. */
    solver->row_accept_decision =
      solver->row_accept_ref_min[1] < solver->row_accept_ref_min[0] * 0.75
        ? 1 : -1;
  }
}

static void kls_row_refactor_acceptance_record_refactor(kls_solver *solver,
                                                        double seconds) {
  if (solver == NULL ||
      kls_row_refactor_env_enabled() ||
      !kls_row_refactor_acceptance_structurally_ready(solver) ||
      solver->common.status < TRILINOS_KLU_OK ||
      solver->common.status == TRILINOS_KLU_SINGULAR) {
    return;
  }
  if (solver->row_accept_decision > 0) {
    /* Post-adoption steady audit, refactor side. The trial samples
       both engines during their warmup (mac_econ: cold column trials
       read 9.1s vs its true 1.1s steady while the row engine trials
       at 3.6s and steadies at 5s - the pair verdict is honest on
       corrupted data). After the adopted row engine demonstrates its
       steady floor, spend ONE forced column refactor to refresh the
       loser's floor and revert on a decisive margin. */
    if (solver->stats.last_refactor_path == KLS_REFACTOR_PATH_ROW) {
      if (solver->row_steady_ref_over == 0 ||
          seconds < solver->row_steady_ref_min) {
        solver->row_steady_ref_min = seconds;
      }
      if (++solver->row_steady_ref_over == 4 &&
          solver->row_reaudit_state == 0) {
        /* By the fourth real row cycle the solve-side audit has obtained its
           three-sample floor and, when useful, its publication probe.  The
           re-audit now takes two warmed column cycles, so retaining four more
           row calls before starting it adds cost without new evidence. */
        solver->row_reaudit_state = 1;
      }
    } else if (solver->row_reaudit_state == 2 ||
               solver->row_reaudit_state == 3 ||
               solver->row_reaudit_state == 5 ||
               solver->row_reaudit_state == 6) {
      /* A forced column cycle.  Retain its demonstrated minimum if a
         refactor-only burst performs more than one refresh before solving. */
      if (solver->row_accept_ref_min[0] <= 0.0 ||
          seconds < solver->row_accept_ref_min[0]) {
        solver->row_accept_ref_min[0] = seconds;
      }
      if (solver->row_accept_ref_min[0] * 1.5 <
          solver->row_steady_ref_min) {
        solver->row_accept_decision = -1;
        solver->row_reaudit_state = 7;
      } else if (solver->options.expected_solves <= 0) {
        solver->row_reaudit_state = 7;
      } else if (solver->row_reaudit_state == 2) {
        solver->row_reaudit_state = 3;
      } else if (solver->row_reaudit_state == 5) {
        solver->row_reaudit_state = 6;
      }
    }
    return;
  }
  if (solver->row_accept_decision != 0) {
    return;
  }
  const int side =
    solver->stats.last_refactor_path == KLS_REFACTOR_PATH_ROW ? 1 : 0;
  solver->row_accept_ref_seconds[side] += seconds;
  solver->row_accept_ref_samples[side]++;
  if (solver->row_accept_ref_samples[side] == 1 ||
      seconds < solver->row_accept_ref_min[side]) {
    solver->row_accept_ref_min[side] = seconds;
  }
  solver->row_accept_pending_side = side + 1;
  if (solver->row_accept_ref_samples[0] > 0 &&
      solver->row_accept_ref_samples[1] > 0) {
    /* Start sampling with the first real refactor and end a decisively
       losing row experiment after one sample per side.  The previous three
       unrecorded warm-up calls pushed the first row/column comparison into
       calls four and five; a short steady-state sample then extrapolated the
       losing trial across the whole SPICE horizon.  Positive row adoption
       still requires the existing two samples (and paired solves), so a cold
       first column sample cannot promote a pathological row engine. */
    const double col_mean = solver->row_accept_ref_seconds[0] /
                            (double)solver->row_accept_ref_samples[0];
    const double row_mean = solver->row_accept_ref_seconds[1] /
                            (double)solver->row_accept_ref_samples[1];
    if (row_mean > 2.0 * col_mean) {
      solver->row_accept_decision = -1;
      solver->row_accept_pending_side = 0;
      return;
    }
  }
  kls_row_refactor_acceptance_try_decide(solver);
}

static void kls_row_refactor_acceptance_record_solve(kls_solver *solver,
                                                     double seconds) {
  if (solver == NULL || solver->row_accept_decision != 0 ||
      solver->row_accept_pending_side == 0) {
    return;
  }
  const int side = solver->row_accept_pending_side - 1;
  solver->row_accept_solve_seconds[side] += seconds;
  solver->row_accept_solve_samples[side]++;
  if (solver->row_accept_solve_samples[side] == 1 ||
      seconds < solver->row_accept_solve_min[side]) {
    solver->row_accept_solve_min[side] = seconds;
  }
  solver->row_accept_pending_side = 0;
  if (side == 0 && solver->row_accept_ref_samples[1] > 0 &&
      solver->row_accept_ref_min[1] >
        solver->row_accept_ref_min[0] +
          solver->row_accept_solve_min[0]) {
    /* The sampled row refactor alone exceeds the complete measured column
       cycle.  Its unknown solve has a nonnegative cost, so no row lifecycle
       can win this comparison.  The omitted row solve itself supplies the
       noise cushion, so the strict inequality needs no tuned margin. */
    solver->row_accept_decision = -1;
    return;
  }
  kls_row_refactor_acceptance_try_decide(solver);
}

static void kls_selected_lean_reaudit_finish(kls_solver *solver) {
  const int solve_pair_ready =
    solver->lean_reaudit_column_solve_samples >=
      solver->lean_reaudit_column_samples &&
    solver->lean_reaudit_row_solve_samples >=
      solver->lean_reaudit_row_samples &&
    solver->lean_reaudit_column_solve_min > 0.0 &&
    solver->lean_reaudit_row_solve_min > 0.0;
  const double row_ref = solver->lean_reaudit_row_min > 0.0
    ? solver->lean_reaudit_row_min
    : (solver->lean_reaudit_samples > 0
         ? solver->lean_reaudit_seconds /
             (double)solver->lean_reaudit_samples
         : 0.0);
  const double row_cycle = solve_pair_ready &&
      solver->lean_reaudit_row_cycle_min > 0.0
    ? solver->lean_reaudit_row_cycle_min
    : row_ref +
        (solve_pair_ready ? solver->lean_reaudit_row_solve_min : 0.0);
  const double column_cycle = solve_pair_ready &&
      solver->lean_reaudit_column_cycle_min > 0.0
    ? solver->lean_reaudit_column_cycle_min
    : solver->lean_reaudit_column_min +
        (solve_pair_ready ? solver->lean_reaudit_column_solve_min : 0.0);
  const double margin = solve_pair_ready ? 0.98 : 0.95;
  if (row_cycle > 0.0 && column_cycle > 0.0 &&
      column_cycle < margin * row_cycle) {
    /* The challenger cycles replaced ordinary recurring work; they are not
       an extra setup bill that the remaining lifecycle must repay.  Once
       four complete cycles establish a guarded win, retaining the slower
       provisional row arm would only compound the mistake. */
    solver->lean_choice = -1;
    solver->lean_probe_arm = 0;
    solver->lean_pair_active = 0;
  }
  solver->lean_reaudit_pending_side = 0;
  solver->lean_reaudit_state = 5;
}

static void kls_declined_lean_reaudit_finish(kls_solver *solver) {
  const int solve_pair_ready =
    solver->lean_reaudit_column_solve_samples >=
      solver->lean_reaudit_column_samples &&
    solver->lean_reaudit_row_solve_samples >=
      solver->lean_reaudit_row_samples &&
    solver->lean_reaudit_column_solve_min > 0.0 &&
    solver->lean_reaudit_row_solve_min > 0.0;
  const double column_cycle = solve_pair_ready &&
      solver->lean_reaudit_column_cycle_min > 0.0
    ? solver->lean_reaudit_column_cycle_min
    : solver->lean_reaudit_column_min +
        (solve_pair_ready ? solver->lean_reaudit_column_solve_min : 0.0);
  const double row_cycle = solve_pair_ready &&
      solver->lean_reaudit_row_cycle_min > 0.0
    ? solver->lean_reaudit_row_cycle_min
    : solver->lean_reaudit_row_min +
        (solve_pair_ready ? solver->lean_reaudit_row_solve_min : 0.0);
  const double adoption_margin = solve_pair_ready ? 0.98 : 0.95;
  if (solver->lean_reaudit_column_samples >=
        KLS_LEAN_LIFECYCLE_REAUDIT_SAMPLES &&
      solver->lean_reaudit_row_samples >=
        KLS_LEAN_LIFECYCLE_REAUDIT_SAMPLES &&
      solver->lean_reaudit_column_min > 0.0 &&
      solver->lean_reaudit_row_min > 0.0 &&
      row_cycle < adoption_margin * column_cycle) {
    solver->lean_choice = solver->lean_reaudit_row_arm;
    solver->lean_pair_active = solver->lean_choice == 2;
  }
  solver->lean_reaudit_pending_side = 0;
  solver->lean_reaudit_state = 5;
}

/* Settle only a very wide complete-cycle separation before the ordinary
   four-by-four audit finishes.  The incumbent arm still receives its full
   four-sample warm window and the challenger receives two complete cycles;
   compare their established timing floors only across a 30-percent gap.
   Close outcomes retain the established audit unchanged. */
static int kls_lean_reaudit_try_decisive(kls_solver *solver) {
  if (solver == NULL ||
      solver->lean_reaudit_column_cycle_samples < 2 ||
      solver->lean_reaudit_row_cycle_samples < 2 ||
      solver->lean_reaudit_column_cycle_min <= 0.0 ||
      solver->lean_reaudit_row_cycle_min <= 0.0) {
    return 0;
  }
  const int selected_row = solver->lean_choice > 0;
  if ((selected_row && solver->lean_reaudit_row_cycle_samples <
                         KLS_LEAN_LIFECYCLE_REAUDIT_SAMPLES) ||
      (!selected_row && solver->lean_reaudit_column_cycle_samples <
                          KLS_LEAN_LIFECYCLE_REAUDIT_SAMPLES)) {
    return 0;
  }
  int choice = 0;
  if (solver->lean_reaudit_row_cycle_min <=
      0.70 * solver->lean_reaudit_column_cycle_min) {
    choice = solver->lean_choice > 0
      ? solver->lean_choice : solver->lean_reaudit_row_arm;
  } else if (solver->lean_reaudit_column_cycle_min <=
             0.70 * solver->lean_reaudit_row_cycle_min) {
    choice = -1;
  }
  if (choice == 0) {
    return 0;
  }
  solver->lean_choice = choice;
  solver->lean_pair_active = choice == 2;
  solver->lean_probe_arm = 0;
  solver->lean_reaudit_pending_side = 0;
  solver->lean_reaudit_pending_ref_seconds = 0.0;
  solver->lean_reaudit_state = 5;
  return 1;
}

/* Pair a lean lifecycle re-audit with the solve that immediately follows
   each sampled numeric.  A row representation can be neutral in
   refactor time yet materially cheaper to solve, while another row engine
   with identical factor work can be much worse.  Recording the realized
   public lifecycle lets the later verdict distinguish those cases without
   inferring anything from matrix dimensions or sparsity shape. */
static void kls_lean_reaudit_record_solve(kls_solver *solver,
                                          double seconds) {
  if (solver == NULL || solver->lean_reaudit_pending_side == 0 ||
      !(seconds > 0.0) || !isfinite(seconds)) {
    return;
  }
  const int side = solver->lean_reaudit_pending_side;
  const double cycle = solver->lean_reaudit_pending_ref_seconds + seconds;
  if (side == 1) {
    if (solver->lean_reaudit_column_solve_samples == 0 ||
        seconds < solver->lean_reaudit_column_solve_min) {
      solver->lean_reaudit_column_solve_min = seconds;
    }
    solver->lean_reaudit_column_solve_samples++;
    if (solver->lean_reaudit_column_cycle_samples == 0 ||
        cycle < solver->lean_reaudit_column_cycle_min) {
      solver->lean_reaudit_column_cycle_min = cycle;
    }
    solver->lean_reaudit_column_cycle_samples++;
  } else if (side == 2) {
    if (solver->lean_reaudit_row_solve_samples == 0 ||
        seconds < solver->lean_reaudit_row_solve_min) {
      solver->lean_reaudit_row_solve_min = seconds;
    }
    solver->lean_reaudit_row_solve_samples++;
    if (solver->lean_reaudit_row_cycle_samples == 0 ||
        cycle < solver->lean_reaudit_row_cycle_min) {
      solver->lean_reaudit_row_cycle_min = cycle;
    }
    solver->lean_reaudit_row_cycle_samples++;
  }
  solver->lean_reaudit_pending_side = 0;
  solver->lean_reaudit_pending_ref_seconds = 0.0;
  if (kls_lean_reaudit_try_decisive(solver)) {
    return;
  }
  if (solver->lean_reaudit_state == 6 && side == 1) {
    kls_selected_lean_reaudit_finish(solver);
  } else if (solver->lean_reaudit_state == 15 && side == 2) {
    kls_declined_lean_reaudit_finish(solver);
  }
}

/* Post-adoption audit of the row-value solve. The pair verdict is two
   samples per side; the steady row solve can later cost multiples of
   the column route (solve-partition mode shifts, pool dispatch at high
   thread counts) with nothing left to correct it — and the solve is
   the term a SPICE cycle charges 99x. Watch the adopted row solve's
   floor: flip to publish-for-solve when the measured column route
   stays decisively cheaper, and when the verdict never saw column
   solves, pay for one published refactor to measure the column route
   before letting the row solve stand. */
static void kls_row_solve_steady_audit(kls_solver *solver, double seconds) {
  if (solver == NULL || solver->row_accept_decision <= 0 ||
      kls_row_refactor_env_enabled()) {
    return;
  }
  if (solver->row_reaudit_state == 3 ||
      solver->row_reaudit_state == 6) {
    /* This solve follows one of the forced column refreshes above.  Compare
       complete recurring cycles: a row refactor can win on its own yet lose
       after its representation-specific solve (and conversely).  Repeating
       the column arm once avoids freezing its cold post-row cache state into
       the lifetime verdict. */
    if (!solver->row_refactor_last_row_solve) {
      if (solver->row_accept_solve_min[0] <= 0.0 ||
          seconds < solver->row_accept_solve_min[0]) {
        solver->row_accept_solve_min[0] = seconds;
      }
      if (solver->row_reaudit_state == 3) {
        solver->row_reaudit_state = 4;
      } else {
        const double column_pair =
          solver->row_accept_ref_min[0] +
          solver->row_accept_solve_min[0];
        const double row_publication =
          solver->row_accept_publish_preferred
            ? solver->row_publish_probe_seconds : 0.0;
        const double row_pair =
          solver->row_steady_ref_min + row_publication +
          solver->row_steady_cycle_solve_min;
        if (column_pair > 0.0 && row_pair > 0.0 &&
            column_pair < 0.98 * row_pair) {
          solver->row_accept_decision = -1;
        }
        solver->row_reaudit_state = 7;
      }
      return;
    }
  }
  if (solver->stats.last_refactor_path == KLS_REFACTOR_PATH_ROW &&
      solver->row_publish_experiment != 2) {
    /* Do not let the one-off published-column probe masquerade as the
       adopted row representation's solve floor.  Ordinary structural
       publication and a settled publish preference are recurring row-cycle
       costs and remain eligible here. */
    solver->row_steady_cycle_solve_samples++;
    if (solver->row_steady_cycle_solve_samples == 1 ||
        seconds < solver->row_steady_cycle_solve_min) {
      solver->row_steady_cycle_solve_min = seconds;
    }
  }
  if (solver->row_accept_publish_preferred ||
      solver->row_publish_experiment == 3) {
    return;
  }
  if (solver->row_refactor_last_row_solve) {
    solver->row_steady_solve_samples++;
    if (solver->row_steady_solve_samples == 1 ||
        seconds < solver->row_steady_solve_min) {
      solver->row_steady_solve_min = seconds;
    }
    if (solver->row_steady_solve_samples < 3 ||
        solver->row_publish_experiment != 0) {
      return;
    }
    if (solver->row_accept_solve_samples[0] > 0) {
      if (solver->row_accept_solve_min[0] <
          0.8 * solver->row_steady_solve_min) {
        solver->row_accept_publish_preferred = 1;
      }
      solver->row_publish_experiment = 3;
    } else {
      solver->row_publish_experiment = 1;
    }
    } else if (solver->row_publish_experiment == 2) {
    /* this solve ran the column route off the published values */
    solver->row_accept_publish_preferred =
      solver->row_publish_probe_seconds > 0.0 &&
      solver->row_publish_probe_seconds + seconds <
        0.98 * solver->row_steady_solve_min;
      solver->row_publish_experiment = 3;
      if (solver->row_accept_publish_preferred &&
          (solver->row_steady_cycle_solve_samples == 0 ||
           seconds < solver->row_steady_cycle_solve_min)) {
        /* Once publication wins, this formerly one-off column solve is the
           recurring solve representation of the adopted row factor. */
        solver->row_steady_cycle_solve_min = seconds;
        solver->row_steady_cycle_solve_samples++;
      }
  }
}

static int kls_row_refactor_should_defer_value_scatter(
  const kls_solver *solver,
  int check_pivots) {
  if (!check_pivots) {
    return 1;
  }
  return solver != NULL && solver->row_refactor_dense_segment_count > 0u;
}

/* Consume a complete suffix of a strict row supernode in one target-major
   pass.  The scalar row walk writes every workspace target once per
   dependency; this is the compact-row equivalent of SubtreeLU's trsv+gemv
   update and keeps the short (2-8 row) circuit supernodes in registers. */
static UF_long kls_lean_row_refactor_snode_step_chain(
  UF_long p,
  UF_long run,
  UF_long dep0,
  double *restrict x,
  double *restrict l_values,
  double **restrict l_lu,
  const UF_long *restrict u_ptr,
  const UF_long *restrict u_cols,
  const double *restrict u_values,
  const double *restrict udiag) {
  if (run < 2u || run > 8u) {
    return 0u;
  }
  const UF_long dep_end = dep0 + run - 1u;

  double multipliers[8];
  const double *trailing_values[8];
  for (UF_long local = 0u; local < run; ++local) {
    const UF_long dep = dep0 + local;
    double candidate = x[dep];
    for (UF_long previous = 0u; previous < local; ++previous) {
      candidate -= multipliers[previous] *
        u_values[u_ptr[dep0 + previous] + local - previous - 1u];
    }
    const double multiplier = candidate / udiag[dep];
    multipliers[local] = multiplier;
    if (l_values != NULL) {
      l_values[p + local] = multiplier;
    } else {
      *l_lu[p + local] = multiplier;
    }
    x[dep] = 0.0;
    trailing_values[local] =
      u_values + u_ptr[dep] + dep_end - dep;
  }

  const UF_long trailing_begin = u_ptr[dep_end];
  const UF_long trailing_len = u_ptr[dep_end + 1u] - trailing_begin;
#define KLS_LEAN_SNODE_CHAIN_TARGET_LOOP(...)                             \
  do {                                                                    \
    for (UF_long offset = 0u; offset < trailing_len; ++offset) {          \
      const UF_long target = u_cols[trailing_begin + offset];             \
      double value = x[target];                                           \
      __VA_ARGS__                                                         \
      x[target] = value;                                                  \
    }                                                                     \
  } while (0)
  switch (run) {
    case 2u:
      KLS_LEAN_SNODE_CHAIN_TARGET_LOOP(
        value -= multipliers[0] * trailing_values[0][offset];
        value -= multipliers[1] * trailing_values[1][offset];);
      break;
    case 3u:
      KLS_LEAN_SNODE_CHAIN_TARGET_LOOP(
        value -= multipliers[0] * trailing_values[0][offset];
        value -= multipliers[1] * trailing_values[1][offset];
        value -= multipliers[2] * trailing_values[2][offset];);
      break;
    case 4u:
      KLS_LEAN_SNODE_CHAIN_TARGET_LOOP(
        value -= multipliers[0] * trailing_values[0][offset];
        value -= multipliers[1] * trailing_values[1][offset];
        value -= multipliers[2] * trailing_values[2][offset];
        value -= multipliers[3] * trailing_values[3][offset];);
      break;
    case 5u:
      KLS_LEAN_SNODE_CHAIN_TARGET_LOOP(
        value -= multipliers[0] * trailing_values[0][offset];
        value -= multipliers[1] * trailing_values[1][offset];
        value -= multipliers[2] * trailing_values[2][offset];
        value -= multipliers[3] * trailing_values[3][offset];
        value -= multipliers[4] * trailing_values[4][offset];);
      break;
    case 6u:
      KLS_LEAN_SNODE_CHAIN_TARGET_LOOP(
        value -= multipliers[0] * trailing_values[0][offset];
        value -= multipliers[1] * trailing_values[1][offset];
        value -= multipliers[2] * trailing_values[2][offset];
        value -= multipliers[3] * trailing_values[3][offset];
        value -= multipliers[4] * trailing_values[4][offset];
        value -= multipliers[5] * trailing_values[5][offset];);
      break;
    case 7u:
      KLS_LEAN_SNODE_CHAIN_TARGET_LOOP(
        value -= multipliers[0] * trailing_values[0][offset];
        value -= multipliers[1] * trailing_values[1][offset];
        value -= multipliers[2] * trailing_values[2][offset];
        value -= multipliers[3] * trailing_values[3][offset];
        value -= multipliers[4] * trailing_values[4][offset];
        value -= multipliers[5] * trailing_values[5][offset];
        value -= multipliers[6] * trailing_values[6][offset];);
      break;
    case 8u:
      KLS_LEAN_SNODE_CHAIN_TARGET_LOOP(
        value -= multipliers[0] * trailing_values[0][offset];
        value -= multipliers[1] * trailing_values[1][offset];
        value -= multipliers[2] * trailing_values[2][offset];
        value -= multipliers[3] * trailing_values[3][offset];
        value -= multipliers[4] * trailing_values[4][offset];
        value -= multipliers[5] * trailing_values[5][offset];
        value -= multipliers[6] * trailing_values[6][offset];
        value -= multipliers[7] * trailing_values[7][offset];);
      break;
    default:
      return 0u;
  }
#undef KLS_LEAN_SNODE_CHAIN_TARGET_LOOP
  return run;
}

#if KLS_HAVE_AVX512_KERNELS
/* A fused width-eight run shares one sorted target stream across all producer
   rows.  Unlike a conventional indexed sparse update, this kernel therefore
   pays one gather and one scatter for several contiguous coefficient loads.
   It evaluates the same balanced tree as the scalar target loop below, lane
   by lane, and leaves the short remainder to that loop. */
__attribute__((target("avx512f"), noinline,
               section(".text.kls_avx512"))) static UF_long
kls_lean_snode_target_width8_avx512(
  UF_long run,
  double *restrict x,
  const UF_long *restrict targets,
  const double *const *restrict trailing_values,
  const double *restrict multipliers,
  UF_long length) {
  if (run != 8u) {
    return 0u;
  }
  const UF_long vector_length = length & ~(UF_long)7u;
  const __m512d m0 = _mm512_set1_pd(multipliers[0]);
  const __m512d mn1 = _mm512_set1_pd(-multipliers[1]);
  const __m512d mn2 = _mm512_set1_pd(-multipliers[2]);
  const __m512d mn3 = _mm512_set1_pd(-multipliers[3]);
  const __m512d m4 = _mm512_set1_pd(multipliers[4]);
  const __m512d m5 = _mm512_set1_pd(multipliers[5]);
  const __m512d m6 = _mm512_set1_pd(multipliers[6]);
  const __m512d m7 = _mm512_set1_pd(multipliers[7]);
  for (UF_long offset = 0u; offset < vector_length; offset += 8u) {
    const __m512i indices = _mm512_loadu_si512(
      (const void *)(targets + offset));
    __m512d value = _mm512_i64gather_pd(indices, x, 8);
    __m512d sum0 = _mm512_fnmadd_pd(
      _mm512_loadu_pd(trailing_values[0] + offset), m0, value);
    __m512d sum1 = _mm512_mul_pd(
      _mm512_loadu_pd(trailing_values[1] + offset), mn1);
    __m512d sum2 = _mm512_mul_pd(
      _mm512_loadu_pd(trailing_values[2] + offset), mn2);
    __m512d sum3 = _mm512_mul_pd(
      _mm512_loadu_pd(trailing_values[3] + offset), mn3);
    sum0 = _mm512_fnmadd_pd(
      _mm512_loadu_pd(trailing_values[4] + offset), m4, sum0);
    sum1 = _mm512_fnmadd_pd(
      _mm512_loadu_pd(trailing_values[5] + offset), m5, sum1);
    sum2 = _mm512_fnmadd_pd(
      _mm512_loadu_pd(trailing_values[6] + offset), m6, sum2);
    sum3 = _mm512_fnmadd_pd(
      _mm512_loadu_pd(trailing_values[7] + offset), m7, sum3);
    value = _mm512_add_pd(_mm512_add_pd(sum0, sum1),
                          _mm512_add_pd(sum2, sum3));
    _mm512_i64scatter_pd(x, indices, value, 8);
  }
  return vector_length;
}

/* The compact row mirror stores the same checked target stream in 16 bits.
   Widening eight targets from one 16-byte load cuts the index traffic of the
   irregular update without changing the gather/scatter arithmetic. */
__attribute__((target("avx512f,avx2"), noinline,
               section(".text.kls_i16_snode"))) static UF_long
kls_lean_snode_target_width8_i16_avx512(
  UF_long run,
  double *restrict x,
  const uint16_t *restrict targets,
  const double *const *restrict trailing_values,
  const double *restrict multipliers,
  UF_long length) {
  if (run != 8u) {
    return 0u;
  }
  const UF_long vector_length = length & ~(UF_long)7u;
  const __m512d m0 = _mm512_set1_pd(multipliers[0]);
  const __m512d mn1 = _mm512_set1_pd(-multipliers[1]);
  const __m512d mn2 = _mm512_set1_pd(-multipliers[2]);
  const __m512d mn3 = _mm512_set1_pd(-multipliers[3]);
  const __m512d m4 = _mm512_set1_pd(multipliers[4]);
  const __m512d m5 = _mm512_set1_pd(multipliers[5]);
  const __m512d m6 = _mm512_set1_pd(multipliers[6]);
  const __m512d m7 = _mm512_set1_pd(multipliers[7]);
  for (UF_long offset = 0u; offset < vector_length; offset += 8u) {
    const __m128i packed = _mm_loadu_si128(
      (const __m128i *)(const void *)(targets + offset));
    const __m256i indices = _mm256_cvtepu16_epi32(packed);
    __m512d value = _mm512_i32gather_pd(indices, x, 8);
    __m512d sum0 = _mm512_fnmadd_pd(
      _mm512_loadu_pd(trailing_values[0] + offset), m0, value);
    __m512d sum1 = _mm512_mul_pd(
      _mm512_loadu_pd(trailing_values[1] + offset), mn1);
    __m512d sum2 = _mm512_mul_pd(
      _mm512_loadu_pd(trailing_values[2] + offset), mn2);
    __m512d sum3 = _mm512_mul_pd(
      _mm512_loadu_pd(trailing_values[3] + offset), mn3);
    sum0 = _mm512_fnmadd_pd(
      _mm512_loadu_pd(trailing_values[4] + offset), m4, sum0);
    sum1 = _mm512_fnmadd_pd(
      _mm512_loadu_pd(trailing_values[5] + offset), m5, sum1);
    sum2 = _mm512_fnmadd_pd(
      _mm512_loadu_pd(trailing_values[6] + offset), m6, sum2);
    sum3 = _mm512_fnmadd_pd(
      _mm512_loadu_pd(trailing_values[7] + offset), m7, sum3);
    value = _mm512_add_pd(_mm512_add_pd(sum0, sum1),
                          _mm512_add_pd(sum2, sum3));
    _mm512_i32scatter_pd(x, indices, value, 8);
  }
  return vector_length;
}

/* Keep a common target vector live while applying two width-eight producer
   groups.  The arithmetic inside each group is identical to the established
   width-eight kernel, and the second group consumes the first group's result
   exactly as it did when the natural supernode was split into two calls. */
__attribute__((target("avx512f"), noinline,
               section(".kls_width16"))) static UF_long
kls_lean_snode_target_width16_avx512(
  double *restrict x,
  const UF_long *restrict targets,
  const double *const *restrict trailing_values,
  const double *restrict multipliers,
  UF_long length) {
  const UF_long vector_length = length & ~(UF_long)7u;
  __m512d m[16];
  for (UF_long local = 0u; local < 16u; ++local) {
    m[local] = _mm512_set1_pd(multipliers[local]);
  }
  for (UF_long offset = 0u; offset < vector_length; offset += 8u) {
    const __m512i indices = _mm512_loadu_si512(
      (const void *)(targets + offset));
    __m512d value = _mm512_i64gather_pd(indices, x, 8);
    for (UF_long group = 0u; group < 2u; ++group) {
      const UF_long base = 8u * group;
      __m512d sum0 = _mm512_fnmadd_pd(
        _mm512_loadu_pd(trailing_values[base] + offset), m[base], value);
      __m512d sum1 = _mm512_mul_pd(
        _mm512_loadu_pd(trailing_values[base + 1u] + offset),
        _mm512_sub_pd(_mm512_setzero_pd(), m[base + 1u]));
      __m512d sum2 = _mm512_mul_pd(
        _mm512_loadu_pd(trailing_values[base + 2u] + offset),
        _mm512_sub_pd(_mm512_setzero_pd(), m[base + 2u]));
      __m512d sum3 = _mm512_mul_pd(
        _mm512_loadu_pd(trailing_values[base + 3u] + offset),
        _mm512_sub_pd(_mm512_setzero_pd(), m[base + 3u]));
      sum0 = _mm512_fnmadd_pd(
        _mm512_loadu_pd(trailing_values[base + 4u] + offset),
        m[base + 4u], sum0);
      sum1 = _mm512_fnmadd_pd(
        _mm512_loadu_pd(trailing_values[base + 5u] + offset),
        m[base + 5u], sum1);
      sum2 = _mm512_fnmadd_pd(
        _mm512_loadu_pd(trailing_values[base + 6u] + offset),
        m[base + 6u], sum2);
      sum3 = _mm512_fnmadd_pd(
        _mm512_loadu_pd(trailing_values[base + 7u] + offset),
        m[base + 7u], sum3);
      value = _mm512_add_pd(_mm512_add_pd(sum0, sum1),
                            _mm512_add_pd(sum2, sum3));
    }
    _mm512_i64scatter_pd(x, indices, value, 8);
  }
  return vector_length;
}

__attribute__((target("avx512f,avx2"), noinline,
               section(".kls_width16"))) static UF_long
kls_lean_snode_target_width16_i16_avx512(
  double *restrict x,
  const uint16_t *restrict targets,
  const double *const *restrict trailing_values,
  const double *restrict multipliers,
  UF_long length) {
  const UF_long vector_length = length & ~(UF_long)7u;
  __m512d m[16];
  for (UF_long local = 0u; local < 16u; ++local) {
    m[local] = _mm512_set1_pd(multipliers[local]);
  }
  for (UF_long offset = 0u; offset < vector_length; offset += 8u) {
    const __m128i packed = _mm_loadu_si128(
      (const __m128i *)(const void *)(targets + offset));
    const __m256i indices = _mm256_cvtepu16_epi32(packed);
    __m512d value = _mm512_i32gather_pd(indices, x, 8);
    for (UF_long group = 0u; group < 2u; ++group) {
      const UF_long base = 8u * group;
      __m512d sum0 = _mm512_fnmadd_pd(
        _mm512_loadu_pd(trailing_values[base] + offset), m[base], value);
      __m512d sum1 = _mm512_mul_pd(
        _mm512_loadu_pd(trailing_values[base + 1u] + offset),
        _mm512_sub_pd(_mm512_setzero_pd(), m[base + 1u]));
      __m512d sum2 = _mm512_mul_pd(
        _mm512_loadu_pd(trailing_values[base + 2u] + offset),
        _mm512_sub_pd(_mm512_setzero_pd(), m[base + 2u]));
      __m512d sum3 = _mm512_mul_pd(
        _mm512_loadu_pd(trailing_values[base + 3u] + offset),
        _mm512_sub_pd(_mm512_setzero_pd(), m[base + 3u]));
      sum0 = _mm512_fnmadd_pd(
        _mm512_loadu_pd(trailing_values[base + 4u] + offset),
        m[base + 4u], sum0);
      sum1 = _mm512_fnmadd_pd(
        _mm512_loadu_pd(trailing_values[base + 5u] + offset),
        m[base + 5u], sum1);
      sum2 = _mm512_fnmadd_pd(
        _mm512_loadu_pd(trailing_values[base + 6u] + offset),
        m[base + 6u], sum2);
      sum3 = _mm512_fnmadd_pd(
        _mm512_loadu_pd(trailing_values[base + 7u] + offset),
        m[base + 7u], sum3);
      value = _mm512_add_pd(_mm512_add_pd(sum0, sum1),
                            _mm512_add_pd(sum2, sum3));
    }
    _mm512_i32scatter_pd(x, indices, value, 8);
  }
  return vector_length;
}
#endif

static KLS_ALWAYS_INLINE double kls_lean_snode_apply_width8_scalar(
  double value,
  const double *const *restrict trailing_values,
  const double *restrict multipliers,
  UF_long offset) {
  double sum0 = fma(-multipliers[0], trailing_values[0][offset], value);
  double sum1 = -multipliers[1] * trailing_values[1][offset];
  double sum2 = -multipliers[2] * trailing_values[2][offset];
  double sum3 = -multipliers[3] * trailing_values[3][offset];
  sum0 = fma(-multipliers[4], trailing_values[4][offset], sum0);
  sum1 = fma(-multipliers[5], trailing_values[5][offset], sum1);
  sum2 = fma(-multipliers[6], trailing_values[6][offset], sum2);
  sum3 = fma(-multipliers[7], trailing_values[7][offset], sum3);
  return (sum0 + sum1) + (sum2 + sum3);
}

static KLS_ALWAYS_INLINE int kls_lean_snode_wait_mask_ready(
  kls_egraph_refactor_shared *shared,
  atomic_uint *restrict grouped_done,
  unsigned int generation,
  unsigned int mask,
  const uint32_t *restrict wait_slots) {
  while (mask != 0u) {
#if defined(__GNUC__) || defined(__clang__)
    const unsigned int local = (unsigned int)__builtin_ctz(mask);
#else
    unsigned int local = 0u;
    while ((mask & (1u << local)) == 0u) {
      local++;
    }
#endif
    mask &= mask - 1u;
    const uint32_t slot = wait_slots[local] - 1u;
    unsigned spin = 0u;
    while (atomic_load_explicit(&grouped_done[slot],
                                memory_order_acquire) != generation) {
      if ((++spin & 255u) == 0u &&
          atomic_load_explicit(&shared->stop,
                               memory_order_acquire) != 0) {
        return 0;
      }
      kls_cpu_relax();
    }
  }
  return 1;
}

/* Consume a true width-16 producer supernode as two numerically identical
   width-eight groups while traversing the shared trailing workspace once.
   The second dependency frontier is acquired only after the first group has
   produced the same eight intermediate values as the former first call. */
#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline, hot, section(".kls_width16")))
#endif
static UF_long kls_lean_row_refactor_snode_step_width16(
  UF_long p,
  UF_long dep0,
  double *restrict x,
  double *restrict l_values,
  double **restrict l_lu,
  const UF_long *restrict u_ptr,
  const UF_long *restrict u_cols,
  const uint16_t *restrict u_cols16,
  const double *restrict u_values,
  const double *restrict udiag,
  int compact_targets,
  int avx512_target,
  kls_egraph_refactor_shared *shared,
  atomic_uint *restrict grouped_done,
  unsigned int generation,
  unsigned int second_wait_mask,
  const uint32_t *restrict second_wait_slots) {
  double multipliers[16];
  const double *trailing_values[16];

#if defined(__GNUC__) || defined(__clang__)
#pragma GCC unroll 8
#endif
  for (UF_long local = 0u; local < 8u; ++local) {
    const UF_long dep = dep0 + local;
    double candidate = x[dep];
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC unroll 8
#endif
    for (UF_long previous = 0u; previous < local; ++previous) {
      candidate -= multipliers[previous] *
        u_values[u_ptr[dep0 + previous] + local - previous - 1u];
    }
    const double multiplier = candidate / udiag[dep];
    multipliers[local] = multiplier;
    if (l_values != NULL) {
      l_values[p + local] = multiplier;
    } else {
      *l_lu[p + local] = multiplier;
    }
    x[dep] = 0.0;
    trailing_values[local] =
      u_values + u_ptr[dep] + 7u - local;
  }

  /* These are the dense triangular entries between the two tiles. */
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC unroll 8
#endif
  for (UF_long local = 0u; local < 8u; ++local) {
    const UF_long dep = dep0 + 8u + local;
    x[dep] = kls_lean_snode_apply_width8_scalar(
      x[dep], trailing_values, multipliers, local);
  }
  if (!kls_lean_snode_wait_mask_ready(
        shared, grouped_done, generation, second_wait_mask,
        second_wait_slots)) {
    return 0u;
  }

#if defined(__GNUC__) || defined(__clang__)
#pragma GCC unroll 8
#endif
  for (UF_long local = 0u; local < 8u; ++local) {
    const UF_long dep = dep0 + 8u + local;
    double candidate = x[dep];
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC unroll 8
#endif
    for (UF_long previous = 0u; previous < local; ++previous) {
      candidate -= multipliers[8u + previous] *
        u_values[u_ptr[dep0 + 8u + previous] +
                 local - previous - 1u];
    }
    const double multiplier = candidate / udiag[dep];
    multipliers[8u + local] = multiplier;
    if (l_values != NULL) {
      l_values[p + 8u + local] = multiplier;
    } else {
      *l_lu[p + 8u + local] = multiplier;
    }
    x[dep] = 0.0;
    trailing_values[8u + local] =
      u_values + u_ptr[dep] + 7u - local;
  }

  const UF_long trailing_begin = u_ptr[dep0 + 15u];
  const UF_long trailing_len = u_ptr[dep0 + 16u] - trailing_begin;
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC unroll 8
#endif
  for (UF_long local = 0u; local < 8u; ++local) {
    trailing_values[local] += 8u;
  }
  UF_long target_offset_begin = 0u;
#if KLS_HAVE_AVX512_KERNELS
  if (avx512_target &&
      trailing_len >= KLS_LEAN_SNODE_AVX512_MIN_TARGETS) {
    target_offset_begin = compact_targets
      ? kls_lean_snode_target_width16_i16_avx512(
          x, u_cols16 + trailing_begin, trailing_values,
          multipliers, trailing_len)
      : kls_lean_snode_target_width16_avx512(
          x, u_cols + trailing_begin, trailing_values,
          multipliers, trailing_len);
  }
#else
  (void)avx512_target;
#endif
  for (UF_long offset = target_offset_begin;
       offset < trailing_len; ++offset) {
    const UF_long target = compact_targets
      ? (UF_long)u_cols16[trailing_begin + offset]
      : u_cols[trailing_begin + offset];
    double value = kls_lean_snode_apply_width8_scalar(
      x[target], trailing_values, multipliers, offset);
    value = kls_lean_snode_apply_width8_scalar(
      value, trailing_values + 8u, multipliers + 8u, offset);
    x[target] = value;
  }
  return 16u;
}

static KLS_ALWAYS_INLINE UF_long kls_lean_row_refactor_snode_step_impl(
  UF_long p,
  UF_long run,
  UF_long dep0,
  double *restrict x,
  double *restrict l_values,
  double **restrict l_lu,
  const UF_long *restrict u_ptr,
  const UF_long *restrict u_cols,
  const uint16_t *restrict u_cols16,
  const double *restrict u_values,
  const double *restrict udiag,
  int avx512_target,
  int compact_targets) {
  if (run < 2u || run > 8u) {
    return 0u;
  }
  const UF_long dep_end = dep0 + run - 1u;

  double multipliers[8];
  const double *trailing_values[8];
  /* The dedicated callers make run constant, but GCC otherwise merges this
     bounded prefix back into a dynamic loop.  Unroll only the 2--8-wide
     triangular solve; the potentially long trailing update remains looped. */
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC unroll 8
#endif
  for (UF_long local = 0u; local < run; ++local) {
    const UF_long dep = dep0 + local;
    double candidate = x[dep];
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC unroll 8
#endif
    for (UF_long previous = 0u; previous < local; ++previous) {
      candidate -= multipliers[previous] *
        u_values[u_ptr[dep0 + previous] + local - previous - 1u];
    }
    const double multiplier = candidate / udiag[dep];
    multipliers[local] = multiplier;
    if (l_values != NULL) {
      l_values[p + local] = multiplier;
    } else {
      *l_lu[p + local] = multiplier;
    }
    x[dep] = 0.0;
    trailing_values[local] =
      u_values + u_ptr[dep] + dep_end - dep;
  }

  const UF_long trailing_begin = u_ptr[dep_end];
  const UF_long trailing_len = u_ptr[dep_end + 1u] - trailing_begin;
  UF_long target_offset_begin = 0u;
#if KLS_HAVE_AVX512_KERNELS
  if (avx512_target && run == 8u) {
    if (compact_targets) {
      target_offset_begin = kls_lean_snode_target_width8_i16_avx512(
        run, x, u_cols16 + trailing_begin, trailing_values,
        multipliers, trailing_len);
    } else {
      target_offset_begin = kls_lean_snode_target_width8_avx512(
        run, x, u_cols + trailing_begin, trailing_values,
        multipliers, trailing_len);
    }
  }
#else
  (void)avx512_target;
#endif
#define KLS_LEAN_SNODE_TARGET_LOOP(...)                                   \
  do {                                                                    \
    for (UF_long offset = target_offset_begin;                            \
         offset < trailing_len; ++offset) {                               \
      const UF_long target = compact_targets                              \
        ? (UF_long)u_cols16[trailing_begin + offset]                      \
        : u_cols[trailing_begin + offset];                                \
      double value = x[target];                                           \
      __VA_ARGS__                                                         \
      x[target] = value;                                                  \
    }                                                                     \
  } while (0)
  /* Dispatch once per producer supernode.  Keeping a distinct target loop
     for each short width lets the compiler retain all multipliers and row
     pointers in registers, matching CKTSO's width-specialized kernels. */
  switch (run) {
    case 2u:
      KLS_LEAN_SNODE_TARGET_LOOP(
        const double sum0 = fma(-multipliers[0],
                                trailing_values[0][offset], value);
        const double sum1 = -multipliers[1] * trailing_values[1][offset];
        value = sum0 + sum1;);
      break;
    case 3u:
      KLS_LEAN_SNODE_TARGET_LOOP(
        const double sum0 = fma(-multipliers[0],
                                trailing_values[0][offset], value);
        const double sum1 = -multipliers[1] * trailing_values[1][offset];
        const double sum2 = -multipliers[2] * trailing_values[2][offset];
        value = (sum0 + sum1) + sum2;);
      break;
    case 4u:
      KLS_LEAN_SNODE_TARGET_LOOP(
        const double sum0 = fma(-multipliers[0],
                                trailing_values[0][offset], value);
        const double sum1 = -multipliers[1] * trailing_values[1][offset];
        const double sum2 = -multipliers[2] * trailing_values[2][offset];
        const double sum3 = -multipliers[3] * trailing_values[3][offset];
        value = (sum0 + sum1) + (sum2 + sum3););
      break;
    case 5u:
      KLS_LEAN_SNODE_TARGET_LOOP(
        double sum0 = fma(-multipliers[0],
                          trailing_values[0][offset], value);
        const double sum1 = -multipliers[1] * trailing_values[1][offset];
        const double sum2 = -multipliers[2] * trailing_values[2][offset];
        const double sum3 = -multipliers[3] * trailing_values[3][offset];
        sum0 = fma(-multipliers[4], trailing_values[4][offset], sum0);
        value = (sum0 + sum1) + (sum2 + sum3););
      break;
    case 6u:
      KLS_LEAN_SNODE_TARGET_LOOP(
        double sum0 = fma(-multipliers[0],
                          trailing_values[0][offset], value);
        double sum1 = -multipliers[1] * trailing_values[1][offset];
        const double sum2 = -multipliers[2] * trailing_values[2][offset];
        const double sum3 = -multipliers[3] * trailing_values[3][offset];
        sum0 = fma(-multipliers[4], trailing_values[4][offset], sum0);
        sum1 = fma(-multipliers[5], trailing_values[5][offset], sum1);
        value = (sum0 + sum1) + (sum2 + sum3););
      break;
    case 7u:
      KLS_LEAN_SNODE_TARGET_LOOP(
        double sum0 = fma(-multipliers[0],
                          trailing_values[0][offset], value);
        double sum1 = -multipliers[1] * trailing_values[1][offset];
        double sum2 = -multipliers[2] * trailing_values[2][offset];
        const double sum3 = -multipliers[3] * trailing_values[3][offset];
        sum0 = fma(-multipliers[4], trailing_values[4][offset], sum0);
        sum1 = fma(-multipliers[5], trailing_values[5][offset], sum1);
        sum2 = fma(-multipliers[6], trailing_values[6][offset], sum2);
        value = (sum0 + sum1) + (sum2 + sum3););
      break;
    case 8u:
      KLS_LEAN_SNODE_TARGET_LOOP(
        double sum0 = fma(-multipliers[0],
                          trailing_values[0][offset], value);
        double sum1 = -multipliers[1] * trailing_values[1][offset];
        double sum2 = -multipliers[2] * trailing_values[2][offset];
        double sum3 = -multipliers[3] * trailing_values[3][offset];
        sum0 = fma(-multipliers[4], trailing_values[4][offset], sum0);
        sum1 = fma(-multipliers[5], trailing_values[5][offset], sum1);
        sum2 = fma(-multipliers[6], trailing_values[6][offset], sum2);
        sum3 = fma(-multipliers[7], trailing_values[7][offset], sum3);
        value = (sum0 + sum1) + (sum2 + sum3););
      break;
    default:
      return 0u;
  }
#undef KLS_LEAN_SNODE_TARGET_LOOP
  return run;
}

#if defined(__GNUC__) || defined(__clang__)
#define KLS_LEAN_SNODE_WIDTH_ATTR __attribute__((noinline, hot))
#else
#define KLS_LEAN_SNODE_WIDTH_ATTR
#endif
#define KLS_DEFINE_LEAN_SNODE_WIDTH(width)                                \
  KLS_LEAN_SNODE_WIDTH_ATTR                                               \
  static UF_long kls_lean_row_refactor_snode_step_width_##width(          \
    UF_long p,                                                            \
    UF_long dep0,                                                         \
    double *restrict x,                                                   \
    double *restrict l_values,                                            \
    double **restrict l_lu,                                               \
    const UF_long *restrict u_ptr,                                        \
    const UF_long *restrict u_cols,                                       \
    const double *restrict u_values,                                      \
    const double *restrict udiag) {                                       \
    return kls_lean_row_refactor_snode_step_impl(                         \
      p, (width), dep0, x, l_values, l_lu, u_ptr, u_cols, NULL, u_values,  \
      udiag, 0, 0);                                                    \
  }

KLS_DEFINE_LEAN_SNODE_WIDTH(2u)
KLS_DEFINE_LEAN_SNODE_WIDTH(3u)
KLS_DEFINE_LEAN_SNODE_WIDTH(4u)
KLS_DEFINE_LEAN_SNODE_WIDTH(5u)
KLS_DEFINE_LEAN_SNODE_WIDTH(6u)
KLS_DEFINE_LEAN_SNODE_WIDTH(7u)
KLS_DEFINE_LEAN_SNODE_WIDTH(8u)

#if KLS_HAVE_AVX512_KERNELS
__attribute__((noinline, section(".text.kls_avx512"))) static UF_long
kls_lean_row_refactor_snode_step_width_8u_avx512(
  UF_long p,
  UF_long dep0,
  double *restrict x,
  double *restrict l_values,
  double **restrict l_lu,
  const UF_long *restrict u_ptr,
  const UF_long *restrict u_cols,
  const double *restrict u_values,
  const double *restrict udiag) {
  return kls_lean_row_refactor_snode_step_impl(
    p, 8u, dep0, x, l_values, l_lu, u_ptr, u_cols, NULL, u_values, udiag, 1, 0);
}
#endif

#if defined(__GNUC__) || defined(__clang__)
#define KLS_LEAN_SNODE_I16_WIDTH_ATTR                                    \
  __attribute__((noinline, hot, section(".text.kls_i16_snode")))
#else
#define KLS_LEAN_SNODE_I16_WIDTH_ATTR
#endif
#define KLS_DEFINE_LEAN_SNODE_I16_WIDTH(width)                            \
  KLS_LEAN_SNODE_I16_WIDTH_ATTR                                           \
  static UF_long kls_lean_row_refactor_snode_step_i16_width_##width(      \
    UF_long p,                                                            \
    UF_long dep0,                                                         \
    double *restrict x,                                                   \
    double *restrict l_values,                                            \
    double **restrict l_lu,                                               \
    const UF_long *restrict u_ptr,                                        \
    const uint16_t *restrict u_cols16,                                    \
    const double *restrict u_values,                                      \
    const double *restrict udiag) {                                       \
    return kls_lean_row_refactor_snode_step_impl(                         \
      p, (width), dep0, x, l_values, l_lu, u_ptr, NULL, u_cols16, u_values,\
      udiag, 0, 1);                                                    \
  }

KLS_DEFINE_LEAN_SNODE_I16_WIDTH(2u)
KLS_DEFINE_LEAN_SNODE_I16_WIDTH(3u)
KLS_DEFINE_LEAN_SNODE_I16_WIDTH(4u)
KLS_DEFINE_LEAN_SNODE_I16_WIDTH(5u)
KLS_DEFINE_LEAN_SNODE_I16_WIDTH(6u)
KLS_DEFINE_LEAN_SNODE_I16_WIDTH(7u)
KLS_DEFINE_LEAN_SNODE_I16_WIDTH(8u)

#if KLS_HAVE_AVX512_KERNELS
__attribute__((noinline, section(".text.kls_i16_snode"))) static UF_long
kls_lean_row_refactor_snode_step_i16_width_8u_avx512(
  UF_long p,
  UF_long dep0,
  double *restrict x,
  double *restrict l_values,
  double **restrict l_lu,
  const UF_long *restrict u_ptr,
  const uint16_t *restrict u_cols16,
  const double *restrict u_values,
  const double *restrict udiag) {
  return kls_lean_row_refactor_snode_step_impl(
    p, 8u, dep0, x, l_values, l_lu, u_ptr, NULL, u_cols16, u_values, udiag, 1, 1);
}
#endif

#undef KLS_DEFINE_LEAN_SNODE_I16_WIDTH
#undef KLS_LEAN_SNODE_I16_WIDTH_ATTR

#undef KLS_DEFINE_LEAN_SNODE_WIDTH
#undef KLS_LEAN_SNODE_WIDTH_ATTR

static KLS_ALWAYS_INLINE UF_long
kls_lean_row_refactor_snode_step_width_specialized(
  UF_long p,
  UF_long run,
  UF_long dep0,
  double *restrict x,
  double *restrict l_values,
  double **restrict l_lu,
  const UF_long *restrict u_ptr,
  const UF_long *restrict u_cols,
  const double *restrict u_values,
  const double *restrict udiag,
  int avx512_target) {
#if !KLS_HAVE_AVX512_KERNELS
  (void)avx512_target;
#endif
#define KLS_LEAN_SNODE_WIDTH_CASE(width)                                  \
    case (width):                                                         \
      return kls_lean_row_refactor_snode_step_width_##width(              \
        p, dep0, x, l_values, l_lu, u_ptr, u_cols, u_values, udiag)
  switch (run) {
    KLS_LEAN_SNODE_WIDTH_CASE(2u);
    KLS_LEAN_SNODE_WIDTH_CASE(3u);
    KLS_LEAN_SNODE_WIDTH_CASE(4u);
    KLS_LEAN_SNODE_WIDTH_CASE(5u);
    KLS_LEAN_SNODE_WIDTH_CASE(6u);
    KLS_LEAN_SNODE_WIDTH_CASE(7u);
    case 8u:
#if KLS_HAVE_AVX512_KERNELS
      /* Eight producer streams amortize one irregular gather/scatter pair
         after four vectors.  Decide before entering the dedicated function,
         leaving short runs and every other width on byte-identical scalar
         kernels. */
      if (avx512_target &&
          u_ptr[dep0 + 8u] - u_ptr[dep0 + 7u] >=
            KLS_LEAN_SNODE_AVX512_MIN_TARGETS) {
        return kls_lean_row_refactor_snode_step_width_8u_avx512(
          p, dep0, x, l_values, l_lu, u_ptr, u_cols, u_values, udiag);
      }
#endif
      return kls_lean_row_refactor_snode_step_width_8u(
        p, dep0, x, l_values, l_lu, u_ptr, u_cols, u_values, udiag);
    default:
      return 0u;
  }
#undef KLS_LEAN_SNODE_WIDTH_CASE
}

static KLS_ALWAYS_INLINE UF_long
kls_lean_row_refactor_snode_step_i16_width_specialized(
  UF_long p,
  UF_long run,
  UF_long dep0,
  double *restrict x,
  double *restrict l_values,
  double **restrict l_lu,
  const UF_long *restrict u_ptr,
  const uint16_t *restrict u_cols16,
  const double *restrict u_values,
  const double *restrict udiag,
  int avx512_target) {
#if !KLS_HAVE_AVX512_KERNELS
  (void)avx512_target;
#endif
#define KLS_LEAN_SNODE_I16_WIDTH_CASE(width)                              \
    case (width):                                                         \
      return kls_lean_row_refactor_snode_step_i16_width_##width(          \
        p, dep0, x, l_values, l_lu, u_ptr, u_cols16, u_values, udiag)
  switch (run) {
    KLS_LEAN_SNODE_I16_WIDTH_CASE(2u);
    KLS_LEAN_SNODE_I16_WIDTH_CASE(3u);
    KLS_LEAN_SNODE_I16_WIDTH_CASE(4u);
    KLS_LEAN_SNODE_I16_WIDTH_CASE(5u);
    KLS_LEAN_SNODE_I16_WIDTH_CASE(6u);
    KLS_LEAN_SNODE_I16_WIDTH_CASE(7u);
    case 8u:
#if KLS_HAVE_AVX512_KERNELS
      if (avx512_target &&
          u_ptr[dep0 + 8u] - u_ptr[dep0 + 7u] >=
            KLS_LEAN_SNODE_AVX512_MIN_TARGETS) {
        return kls_lean_row_refactor_snode_step_i16_width_8u_avx512(
          p, dep0, x, l_values, l_lu, u_ptr, u_cols16, u_values, udiag);
      }
#endif
      return kls_lean_row_refactor_snode_step_i16_width_8u(
        p, dep0, x, l_values, l_lu, u_ptr, u_cols16, u_values, udiag);
    default:
      return 0u;
  }
#undef KLS_LEAN_SNODE_I16_WIDTH_CASE
}

static KLS_ALWAYS_INLINE UF_long kls_row_refactor_input_pos_at(
  const kls_solver *solver,
  UF_long pos) {
  return solver->row_refactor_input_pos32 != NULL
    ? (UF_long)solver->row_refactor_input_pos32[pos]
    : solver->row_refactor_input_pos[pos];
}

/* Compact input coordinates need not have the same width.  Circuit factors
   commonly keep row-local columns below 2^16 after their value positions
   have crossed that boundary, so retain two dense streams instead of
   inflating both coordinates to UF_long. */
static int kls_build_lean_mixed_input_indices(kls_solver *solver) {
  if (solver == NULL || solver->n == 0u ||
      solver->n > (UF_long)UINT16_MAX + 1u ||
      solver->row_refactor_input_ptr == NULL ||
      solver->row_refactor_input_cols == NULL ||
      solver->row_refactor_input_pos == NULL) {
    return 0;
  }
  if (solver->row_refactor_input_cols16 != NULL &&
      solver->row_refactor_input_pos32 != NULL) {
    return 1;
  }
  const UF_long count = solver->row_refactor_input_ptr[solver->n];
  if (count > (UF_long)UINT32_MAX ||
      count > (UF_long)(SIZE_MAX / sizeof(uint16_t)) ||
      count > (UF_long)(SIZE_MAX / sizeof(uint32_t))) {
    return 0;
  }
  uint16_t *cols = solver->row_refactor_input_cols16 == NULL
    ? (uint16_t *)malloc((size_t)(count > 0u ? count : 1u) * sizeof(*cols))
    : NULL;
  uint32_t *pos = solver->row_refactor_input_pos32 == NULL
    ? (uint32_t *)malloc((size_t)(count > 0u ? count : 1u) * sizeof(*pos))
    : NULL;
  if ((solver->row_refactor_input_cols16 == NULL && cols == NULL) ||
      (solver->row_refactor_input_pos32 == NULL && pos == NULL)) {
    free(cols);
    free(pos);
    return 0;
  }
  for (UF_long p = 0u; p < count; ++p) {
    if (solver->row_refactor_input_cols[p] > (UF_long)UINT16_MAX ||
        solver->row_refactor_input_pos[p] > (UF_long)UINT32_MAX) {
      free(cols);
      free(pos);
      return 0;
    }
    if (cols != NULL) cols[p] = (uint16_t)solver->row_refactor_input_cols[p];
    if (pos != NULL) pos[p] = (uint32_t)solver->row_refactor_input_pos[p];
  }
  if (cols != NULL) solver->row_refactor_input_cols16 = cols;
  if (pos != NULL) solver->row_refactor_input_pos32 = pos;
  return 1;
}

static int kls_build_lean_row_i16_indices(kls_solver *solver) {
  if (solver == NULL || solver->n == 0u ||
      solver->n > (UF_long)UINT16_MAX + 1u ||
      solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_u_ptr == NULL ||
      solver->row_refactor_input_ptr == NULL ||
      solver->row_refactor_l_cols == NULL ||
      solver->row_refactor_u_cols == NULL ||
      solver->row_refactor_level_rows == NULL) {
    return 0;
  }
  if (solver->row_refactor_l_cols16 != NULL &&
      solver->row_refactor_u_cols16 != NULL) {
    return kls_build_lean_mixed_input_indices(solver);
  }
  const UF_long n = solver->n;
  const UF_long lcount = solver->row_refactor_l_ptr[n];
  const UF_long ucount = solver->row_refactor_u_ptr[n];
  const UF_long input_count = solver->row_refactor_input_ptr[n];
  const int compact_ptrs =
    lcount <= (UF_long)UINT16_MAX &&
    ucount <= (UF_long)UINT16_MAX &&
    input_count <= (UF_long)UINT16_MAX;
  const int compact_ptrs32 =
    !compact_ptrs &&
    kls_repeated_update_workload(&solver->options) &&
    lcount <= (UF_long)UINT32_MAX &&
    ucount <= (UF_long)UINT32_MAX &&
    input_count <= (UF_long)UINT32_MAX;
  uint16_t *lcols = (uint16_t *)malloc(
    (size_t)(lcount > 0u ? lcount : 1u) * sizeof(*lcols));
  uint16_t *ucols = (uint16_t *)malloc(
    (size_t)(ucount > 0u ? ucount : 1u) * sizeof(*ucols));
  uint16_t *lptr = compact_ptrs ? (uint16_t *)malloc(
    ((size_t)n + 1u) * sizeof(*lptr)) : NULL;
  uint16_t *uptr = compact_ptrs ? (uint16_t *)malloc(
    ((size_t)n + 1u) * sizeof(*uptr)) : NULL;
  uint16_t *input_ptr = compact_ptrs ? (uint16_t *)malloc(
    ((size_t)n + 1u) * sizeof(*input_ptr)) : NULL;
  uint32_t *lptr32 = compact_ptrs32 ? (uint32_t *)malloc(
    ((size_t)n + 1u) * sizeof(*lptr32)) : NULL;
  uint32_t *uptr32 = compact_ptrs32 ? (uint32_t *)malloc(
    ((size_t)n + 1u) * sizeof(*uptr32)) : NULL;
  uint32_t *input_ptr32 = compact_ptrs32 ? (uint32_t *)malloc(
    ((size_t)n + 1u) * sizeof(*input_ptr32)) : NULL;
  uint16_t *level_rows =
    compact_ptrs
      ? (uint16_t *)malloc((size_t)n * sizeof(*level_rows)) : NULL;
  if (lcols == NULL || ucols == NULL ||
      (compact_ptrs && (lptr == NULL || uptr == NULL ||
                        input_ptr == NULL)) ||
      (compact_ptrs32 && (lptr32 == NULL || uptr32 == NULL ||
                          input_ptr32 == NULL))) {
    free(lcols);
    free(ucols);
    free(lptr);
    free(uptr);
    free(input_ptr);
    free(lptr32);
    free(uptr32);
    free(input_ptr32);
    free(level_rows);
    return 0;
  }
  if (compact_ptrs) {
    for (UF_long row = 0u; row <= n; ++row) {
      if (solver->row_refactor_l_ptr[row] > (UF_long)UINT16_MAX ||
          solver->row_refactor_u_ptr[row] > (UF_long)UINT16_MAX ||
          solver->row_refactor_input_ptr[row] > (UF_long)UINT16_MAX) {
        free(lcols);
        free(ucols);
        free(lptr);
        free(uptr);
        free(input_ptr);
        free(lptr32);
        free(uptr32);
        free(input_ptr32);
        free(level_rows);
        return 0;
      }
      lptr[row] = (uint16_t)solver->row_refactor_l_ptr[row];
      uptr[row] = (uint16_t)solver->row_refactor_u_ptr[row];
      input_ptr[row] = (uint16_t)solver->row_refactor_input_ptr[row];
    }
    if (level_rows != NULL) {
      for (UF_long pos = 0u; pos < n; ++pos) {
        const UF_long row = solver->row_refactor_level_rows[pos];
        if (row > (UF_long)UINT16_MAX) {
          free(level_rows);
          level_rows = NULL;
          break;
        }
        level_rows[pos] = (uint16_t)row;
      }
    }
  } else if (compact_ptrs32) {
    for (UF_long row = 0u; row <= n; ++row) {
      lptr32[row] = (uint32_t)solver->row_refactor_l_ptr[row];
      uptr32[row] = (uint32_t)solver->row_refactor_u_ptr[row];
      input_ptr32[row] = (uint32_t)solver->row_refactor_input_ptr[row];
    }
  }
  for (UF_long p = 0u; p < lcount; ++p) {
    const UF_long col = solver->row_refactor_l_cols[p];
    if (col > (UF_long)UINT16_MAX) {
      free(lcols);
      free(ucols);
      free(lptr);
      free(uptr);
      free(input_ptr);
      free(lptr32);
      free(uptr32);
      free(input_ptr32);
      free(level_rows);
      return 0;
    }
    lcols[p] = (uint16_t)col;
  }
  for (UF_long p = 0u; p < ucount; ++p) {
    const UF_long col = solver->row_refactor_u_cols[p];
    if (col > (UF_long)UINT16_MAX) {
      free(lcols);
      free(ucols);
      free(lptr);
      free(uptr);
      free(input_ptr);
      free(lptr32);
      free(uptr32);
      free(input_ptr32);
      free(level_rows);
      return 0;
    }
    ucols[p] = (uint16_t)col;
  }
  free(solver->row_refactor_l_cols16);
  free(solver->row_refactor_u_cols16);
  free(solver->row_refactor_l_ptr16);
  free(solver->row_refactor_u_ptr16);
  free(solver->row_refactor_input_ptr16);
  free(solver->row_refactor_l_ptr32);
  free(solver->row_refactor_u_ptr32);
  free(solver->row_refactor_input_ptr32);
  free(solver->row_refactor_level_rows16);
  solver->row_refactor_l_cols16 = lcols;
  solver->row_refactor_u_cols16 = ucols;
  solver->row_refactor_l_ptr16 = lptr;
  solver->row_refactor_u_ptr16 = uptr;
  solver->row_refactor_input_ptr16 = input_ptr;
  solver->row_refactor_l_ptr32 = lptr32;
  solver->row_refactor_u_ptr32 = uptr32;
  solver->row_refactor_input_ptr32 = input_ptr32;
  solver->row_refactor_level_rows16 = level_rows;
  return kls_build_lean_mixed_input_indices(solver);
}

/* Minimal persistent-pool variant of the lean row walk.  Rows are assigned
   from the true topological level order in a static round-robin stream.  A
   release marker after each row is the only scheduling metadata consumed at
   run time: unlike the general row engine this has no groups, ready queue, or
   per-level barriers. */
static KLS_ALWAYS_INLINE int kls_lean_parallel_process_row(
  kls_egraph_refactor_worker *worker,
  UF_long row,
  unsigned int generation) {
  if (worker == NULL || worker->shared == NULL || worker->x == NULL) {
    return 0;
  }
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  if (solver == NULL || row >= solver->n || shared->values == NULL ||
      solver->row_refactor_input_ptr == NULL ||
      solver->row_refactor_input_cols == NULL ||
      solver->row_refactor_input_pos == NULL ||
      solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_l_cols == NULL ||
      solver->row_refactor_l_values == NULL ||
      solver->row_refactor_l_row_values == NULL ||
      solver->row_refactor_u_ptr == NULL ||
      solver->row_refactor_u_cols == NULL ||
      solver->row_refactor_u_values == NULL ||
      solver->row_refactor_u_row_values == NULL ||
      solver->numeric == NULL || solver->numeric->Udiag == NULL ||
      solver->numeric->Pnum == NULL) {
    return 0;
  }
  double *restrict x = worker->x;
  const UF_long *restrict in_ptr = solver->row_refactor_input_ptr;
  const UF_long *restrict in_cols = solver->row_refactor_input_cols;
  const UF_long *restrict in_pos = solver->row_refactor_input_pos;
  const UF_long *restrict l_ptr = solver->row_refactor_l_ptr;
  const UF_long *restrict l_cols = solver->row_refactor_l_cols;
  double *restrict l_val = solver->row_refactor_l_row_values;
  double **restrict l_lu = solver->row_refactor_l_values;
  const UF_long *restrict u_ptr = solver->row_refactor_u_ptr;
  const UF_long *restrict u_cols = solver->row_refactor_u_cols;
  double *restrict u_val = solver->row_refactor_u_row_values;
  double **restrict u_lu = solver->row_refactor_u_values;
  double *restrict udiag = (double *)solver->numeric->Udiag;
  const int retain_l_row_values =
    shared->row_refactor_defer_value_scatter ||
    (!shared->lean_symmetric_scalar_fringe_mode &&
     !shared->lean_compact_match_mode);
  const double row_rs_inv = shared->rs != NULL
    ? 1.0 / shared->rs[solver->numeric->Pnum[row]] : 1.0;

  for (UF_long p = in_ptr[row]; p < in_ptr[row + 1u]; ++p) {
    x[in_cols[p]] = shared->values[in_pos[p]] * row_rs_inv;
  }
  const unsigned char *restrict snode_run = solver->lean_snode_run;
  for (UF_long p = l_ptr[row]; p < l_ptr[row + 1u];) {
    const UF_long run = snode_run != NULL ? (UF_long)snode_run[p] : 0u;
    if (run >= 2u && run <= 8u && p + run <= l_ptr[row + 1u]) {
      for (UF_long local = 0u; local < run; ++local) {
        const UF_long dep = l_cols[p + local];
        unsigned spin = 0u;
        if (dep >= row) {
          return 0;
        }
        const unsigned int owner = shared->lean_done[dep].owner;
        if (owner != (unsigned int)worker->tid) {
          while (atomic_load_explicit(&shared->lean_done[dep].generation,
                                      memory_order_acquire) != generation) {
            if ((++spin & 255u) == 0u &&
                atomic_load_explicit(&shared->stop,
                                     memory_order_acquire) != 0) {
              return 0;
            }
            kls_cpu_relax();
          }
        }
      }
      if (kls_lean_row_refactor_snode_step_chain(
            p, run, l_cols[p], x,
            retain_l_row_values ? l_val : NULL, l_lu,
            u_ptr, u_cols, u_val, udiag) != run) {
        return 0;
      }
      p += run;
      continue;
    }
    const UF_long dep = l_cols[p];
    unsigned spin = 0u;
    if (dep >= row) {
      return 0;
    }
    const unsigned int owner = shared->lean_done[dep].owner;
    if (owner != (unsigned int)worker->tid) {
      while (atomic_load_explicit(&shared->lean_done[dep].generation,
                                  memory_order_acquire) != generation) {
        if ((++spin & 255u) == 0u &&
            atomic_load_explicit(&shared->stop, memory_order_acquire) != 0) {
          return 0;
        }
        kls_cpu_relax();
      }
    }
    const double lik = x[dep] / udiag[dep];
    x[dep] = 0.0;
    if (retain_l_row_values) {
      l_val[p] = lik;
    }
    if (!shared->row_refactor_defer_value_scatter) {
      *l_lu[p] = lik;
    }
    if (lik != 0.0) {
      for (UF_long q = u_ptr[dep]; q < u_ptr[dep + 1u]; ++q) {
        x[u_cols[q]] -= lik * u_val[q];
      }
    }
    p++;
  }

  const double pivot = x[row];
  x[row] = 0.0;
  udiag[row] = pivot;
  for (UF_long q = u_ptr[row]; q < u_ptr[row + 1u]; ++q) {
    const UF_long col = u_cols[q];
    const double value = x[col];
    u_val[q] = value;
    if (!shared->row_refactor_defer_value_scatter) {
      *u_lu[q] = value;
    }
    x[col] = 0.0;
  }
  if (pivot == 0.0) {
    kls_egraph_refactor_record_singular(shared, row, row);
  }
  return 1;
}

/* A retained row permutation can consume caller values without rebuilding
   the permuted CSC frame when every user/internal position fits the compact
   inverse map.  The builders below validate the permutation bijection and
   every row/off-diagonal position before publishing the maps; no input-shape
   classification is required. */
static int kls_direct_user_value_maps_capable(const kls_solver *solver) {
  return solver != NULL && solver->symbolic != NULL &&
    solver->numeric != NULL &&
    kls_repeated_update_workload(&solver->options) &&
    solver->orientation == KLS_ORIENTATION_NORMAL &&
    solver->input_format == KLS_INPUT_CSC && solver->input_to_csc != NULL &&
    solver->row_perm != NULL && solver->row_scale == NULL &&
    solver->col_scale == NULL && solver->common.scale <= 0 &&
    solver->numeric->Rs == NULL && solver->nnz <= (UF_long)UINT16_MAX;
}

static int kls_deferred_lean_value_prep_capable(
  const kls_solver *solver) {
  return solver != NULL && solver->values != NULL &&
    solver->input_to_csc != NULL &&
    solver->prepared_value_input_pos != NULL &&
    solver->row_scale == NULL && solver->col_scale == NULL &&
    solver->common.scale <= 0 && solver->numeric != NULL &&
    solver->numeric->Rs == NULL &&
    kls_repeated_update_workload(&solver->options);
}

static int kls_build_compact_match_direct_value_maps(kls_solver *solver) {
  const int direct_capable =
    kls_direct_user_value_maps_capable(solver);
  /* The direct-numeric row builder records internal CSC positions without
     constructing the generic refactor-column map.  Its compact Offx array is
     therefore an internal-position staging map until the inverse user map
     below validates and translates it. */
  const int direct_pattern_positions =
    direct_capable && solver->refactor_input_pos == NULL &&
    (solver->symbolic->nblocks == 1u ||
     solver->lean_btf_off_input_pos != NULL ||
     solver->compact_match_offdiag_user_pos != NULL);
  const int generic_offdiag_positions =
    solver != NULL && solver->refactor_input_pos != NULL &&
    solver->refactor_col_ptr != NULL &&
    solver->refactor_block_start != NULL;
  if (solver == NULL || solver->input_to_csc == NULL ||
      solver->row_scale != NULL || solver->col_scale != NULL ||
      solver->row_refactor_input_ptr == NULL ||
      solver->row_refactor_input_pos == NULL ||
      solver->numeric == NULL || solver->numeric->Offp == NULL ||
      (!direct_pattern_positions && !generic_offdiag_positions) ||
      solver->nnz > (UF_long)UINT16_MAX ||
      (!direct_capable)) {
    return 0;
  }
  if (solver->row_refactor_input_user_pos != NULL &&
      solver->compact_match_offdiag_user_pos != NULL) {
    return 1;
  }
  const UF_long nnz = solver->nnz;
  const UF_long row_count = solver->row_refactor_input_ptr[solver->n];
  const UF_long offcount = solver->numeric->Offp[solver->n];
  uint16_t *inverse = (uint16_t *)malloc((size_t)nnz * sizeof(*inverse));
  uint16_t *row_user = (uint16_t *)malloc(
    (size_t)(row_count > 0u ? row_count : 1u) * sizeof(*row_user));
  uint16_t *offdiag_user = (uint16_t *)malloc(
    (size_t)(offcount > 0u ? offcount : 1u) * sizeof(*offdiag_user));
  if (inverse == NULL || row_user == NULL || offdiag_user == NULL) {
    free(inverse); free(row_user); free(offdiag_user);
    return 0;
  }
  for (UF_long p = 0u; p < nnz; ++p) {
    inverse[p] = UINT16_MAX;
  }
  for (UF_long user = 0u; user < nnz; ++user) {
    const UF_long internal = solver->input_to_csc[user];
    if (internal >= nnz || inverse[internal] != UINT16_MAX) {
      free(inverse); free(row_user); free(offdiag_user);
      return 0;
    }
    inverse[internal] = (uint16_t)user;
  }
  for (UF_long p = 0u; p < row_count; ++p) {
    const UF_long internal = solver->row_refactor_input_pos[p];
    if (internal >= nnz || inverse[internal] == UINT16_MAX) {
      free(inverse); free(row_user); free(offdiag_user);
      return 0;
    }
    row_user[p] = inverse[internal];
  }
  if (direct_pattern_positions) {
    for (UF_long poff = 0u; poff < offcount; ++poff) {
      const UF_long internal = solver->lean_btf_off_input_pos != NULL
        ? (UF_long)solver->lean_btf_off_input_pos[poff]
        : (UF_long)solver->compact_match_offdiag_user_pos[poff];
      if (internal >= nnz || inverse[internal] == UINT16_MAX) {
        free(inverse); free(row_user); free(offdiag_user);
        return 0;
      }
      offdiag_user[poff] = inverse[internal];
    }
  } else {
    for (UF_long k = 0u; k < solver->n; ++k) {
      UF_long poff = solver->numeric->Offp[k];
      const UF_long poff_end = solver->numeric->Offp[k + 1u];
      for (UF_long p = solver->refactor_col_ptr[k];
           p < solver->refactor_block_start[k]; ++p) {
        if (poff >= poff_end || p >= nnz) {
          free(inverse); free(row_user); free(offdiag_user);
          return 0;
        }
        const UF_long internal = solver->refactor_input_pos[p];
        if (internal >= nnz || inverse[internal] == UINT16_MAX) {
          free(inverse); free(row_user); free(offdiag_user);
          return 0;
        }
        offdiag_user[poff++] = inverse[internal];
      }
      if (poff != poff_end) {
        free(inverse); free(row_user); free(offdiag_user);
        return 0;
      }
    }
  }
  free(inverse);
  {
    uint32_t *packed = (uint32_t *)malloc(
      (size_t)(row_count > 0u ? row_count : 1u) * sizeof(*packed));
    if (packed != NULL) {
      UF_long p = 0u;
      for (; p < row_count; ++p) {
        const UF_long col = solver->row_refactor_input_cols[p];
        if (col > (UF_long)UINT16_MAX) {
          break;
        }
        packed[p] = ((uint32_t)row_user[p] << 16u) | (uint32_t)col;
      }
      if (p == row_count) {
        free(solver->row_refactor_input_col_direct32);
        solver->row_refactor_input_col_direct32 = packed;
      } else {
        free(packed);
      }
    }
  }
  solver->row_refactor_input_user_pos = row_user;
  if (direct_pattern_positions && solver->lean_btf_off_input_pos == NULL) {
    free(solver->compact_match_offdiag_user_pos);
  }
  solver->compact_match_offdiag_user_pos = offdiag_user;
  return 1;
}

static int kls_lean_parallel_refresh_offdiag_worker(
  kls_egraph_refactor_worker *worker,
  UF_long stride) {
  if (worker == NULL || worker->shared == NULL ||
      !worker->shared->lean_parallel_offdiag_mode) {
    return 1;
  }
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  if (solver == NULL || shared->values == NULL || stride == 0u ||
      solver->symbolic == NULL || solver->symbolic->nblocks <= 1u ||
      solver->numeric == NULL || solver->numeric->Offp == NULL ||
      solver->numeric->Offx == NULL) {
    return 0;
  }
  if (solver->lean_btf_off_input_pos != NULL &&
      ((solver->common.scale > 0 && shared->rs != NULL) ||
       (solver->common.scale <= 0 && shared->rs == NULL))) {
    const UF_long offcount = solver->numeric->Offp[solver->n];
    double *restrict offx = (double *)solver->numeric->Offx;
    const uint32_t *restrict input = solver->lean_btf_off_input_pos;
    for (UF_long p = (UF_long)worker->tid; p < offcount; p += stride) {
      const UF_long input_pos = (UF_long)input[p];
      if (input_pos >= solver->nnz) {
        return 0;
      }
      double value = shared->values[input_pos];
      if (shared->rs != NULL) {
        const UF_long row = solver->row_idx[input_pos];
        if (row >= solver->n || shared->rs[row] == 0.0) {
          return 0;
        }
        value /= shared->rs[row];
      }
      offx[p] = value;
    }
    return 1;
  }
  if (solver->input_to_csc == NULL && solver->row_perm == NULL &&
      solver->compact_match_offdiag_user_pos != NULL &&
      solver->common.scale <= 0 && solver->numeric->Rs == NULL) {
    const UF_long offcount = solver->numeric->Offp[solver->n];
    double *restrict offx = (double *)solver->numeric->Offx;
    const uint16_t *restrict input =
      solver->compact_match_offdiag_user_pos;
    for (UF_long p = (UF_long)worker->tid; p < offcount; p += stride) {
      if ((UF_long)input[p] >= solver->nnz) {
        return 0;
      }
      const UF_long input_pos = (UF_long)input[p];
      const double value = shared->values[input_pos];
      offx[p] = value;
    }
    return 1;
  }
  if (solver->refactor_col_ptr == NULL ||
      solver->refactor_block_start == NULL ||
      solver->refactor_input_pos == NULL || solver->row_idx == NULL) {
    return 0;
  }
  double *restrict offx = (double *)solver->numeric->Offx;
  for (UF_long k = (UF_long)worker->tid; k < solver->n; k += stride) {
    const UF_long begin = solver->refactor_col_ptr[k];
    const UF_long end = solver->refactor_block_start[k];
    if (end < begin || end > solver->refactor_col_ptr[k + 1u]) {
      return 0;
    }
    UF_long poff = solver->numeric->Offp[k];
    const UF_long poff_end = solver->numeric->Offp[k + 1u];
    for (UF_long p = begin; p < end; ++p) {
      const UF_long input = solver->refactor_input_pos[p];
      if (input >= solver->nnz || poff >= poff_end) {
        return 0;
      }
      double value = shared->values[input];
      if (shared->rs != NULL) {
        const UF_long row = solver->row_idx[input];
        if (row >= solver->n || shared->rs[row] == 0.0) {
          return 0;
        }
        value /= shared->rs[row];
      }
      offx[poff++] = value;
    }
    if (poff != poff_end) {
      return 0;
    }
  }
  return 1;
}

static int kls_parallel_offdiag_is_profitable(kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL ||
      solver->numeric == NULL || solver->numeric->Offp == NULL ||
      solver->symbolic->nblocks <= 1u || solver->options.threads < 2) {
    return 0;
  }
  if (solver->lean_parallel_offdiag_decision != 0) {
    return solver->lean_parallel_offdiag_decision > 0;
  }

  const UF_long *offp = solver->numeric->Offp;
  const UF_long total = offp[solver->n];
  const int threads = solver->options.threads;
  UF_long *worker_work = (UF_long *)calloc(
    (size_t)threads, sizeof(*worker_work));
  if (worker_work == NULL) {
    solver->lean_parallel_offdiag_decision = -1;
    return 0;
  }
  UF_long shared_boundaries = 0u;
  UF_long previous_end = 0u;
  int previous_owner = -1;
  int valid = 1;
  for (UF_long col = 0u; col < solver->n; ++col) {
    const UF_long begin = offp[col];
    const UF_long end = offp[col + 1u];
    if (begin > end || end > total) {
      valid = 0;
      break;
    }
    const UF_long length = end - begin;
    if (length == 0u) {
      continue;
    }
    const int owner = (int)(col % (UF_long)threads);
    if (length > UF_long_max - worker_work[owner]) {
      valid = 0;
      break;
    }
    worker_work[owner] += length;
    /* Adjacent streams owned by different workers write the same 64-byte
       line exactly when their boundary is not eight-double aligned. */
    if (previous_owner >= 0 && previous_end == begin &&
        previous_owner != owner && (begin & 7u) != 0u) {
      shared_boundaries++;
    }
    previous_end = end;
    previous_owner = owner;
  }
  UF_long max_work = 0u;
  for (int tid = 0; tid < threads; ++tid) {
    if (worker_work[tid] > max_work) {
      max_work = worker_work[tid];
    }
  }
  /* Parallel dispatch needs enough values per worker, a bounded critical
     stream, and at least 32 useful writes per contested cache-line boundary.
     These are direct costs of this representation, independent of matrix
     dimensions, density, or the relative size of any BTF block. */
  const int profitable = valid &&
    (double)total >= 512.0 * (double)threads &&
    (double)max_work * (double)threads <= 1.5 * (double)total &&
    (shared_boundaries == 0u || shared_boundaries <= total / 32u);
  free(worker_work);
  solver->lean_parallel_offdiag_decision = profitable ? 1 : -1;
  return profitable;
}

static KLS_ALWAYS_INLINE int kls_compact_match_lean_parallel_process_row(
  kls_egraph_refactor_worker *worker,
  UF_long row,
  unsigned int generation) {
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  double *restrict x = worker->x;
  const UF_long *restrict in_ptr = solver->row_refactor_input_ptr;
  const UF_long *restrict in_cols = solver->row_refactor_input_cols;
  const uint16_t *restrict in_cols16 =
    solver->row_refactor_input_cols16;
  const uint32_t *restrict in_packed =
    solver->row_refactor_input_col_direct32;
  const UF_long *restrict l_ptr = solver->row_refactor_l_ptr;
  const UF_long *restrict l_cols = solver->row_refactor_l_cols;
  double **restrict l_lu = solver->row_refactor_l_values;
  double *restrict l_val = solver->row_refactor_l_row_values;
  const UF_long *restrict u_ptr = solver->row_refactor_u_ptr;
  const UF_long *restrict u_cols = solver->row_refactor_u_cols;
  double *restrict u_val = solver->row_refactor_u_row_values;
  double **restrict u_lu = solver->row_refactor_u_values;
  double *restrict udiag = (double *)solver->numeric->Udiag;

  if (solver->lean_user_values_active) {
    const uint16_t *restrict in_user =
      solver->row_refactor_input_user_pos;
    if (in_packed != NULL) {
      for (UF_long p = in_ptr[row]; p < in_ptr[row + 1u]; ++p) {
        const uint32_t packed = in_packed[p];
        x[packed & UINT32_C(0xffff)] = shared->values[packed >> 16u];
      }
    } else if (in_cols16 != NULL) {
      for (UF_long p = in_ptr[row]; p < in_ptr[row + 1u]; ++p) {
        x[in_cols16[p]] = shared->values[in_user[p]];
      }
    } else {
      for (UF_long p = in_ptr[row]; p < in_ptr[row + 1u]; ++p) {
        x[in_cols[p]] = shared->values[in_user[p]];
      }
    }
  } else {
    const UF_long *restrict in_pos = solver->row_refactor_input_pos;
    for (UF_long p = in_ptr[row]; p < in_ptr[row + 1u]; ++p) {
      x[in_cols[p]] = shared->values[in_pos[p]];
    }
  }
  for (UF_long p = l_ptr[row]; p < l_ptr[row + 1u]; ++p) {
    const UF_long dep = l_cols[p];
    if (shared->lean_grouped_done_mode) {
      const uint32_t token = solver->lean_parallel_grouped_token[dep];
      if ((unsigned int)(token >> KLS_LEAN_GROUPED_OWNER_SHIFT) !=
          (unsigned int)worker->tid) {
        const uint32_t slot = token & KLS_LEAN_GROUPED_SLOT_MASK;
        while (atomic_load_explicit(
                 &solver->lean_parallel_grouped_done[slot],
                 memory_order_acquire) != generation) {
          kls_cpu_relax();
        }
      }
    } else if (shared->lean_done[dep].owner !=
               (unsigned int)worker->tid) {
      while (atomic_load_explicit(&shared->lean_done[dep].generation,
                                  memory_order_acquire) != generation) {
        kls_cpu_relax();
      }
    }
    const double lik = x[dep] / udiag[dep];
    x[dep] = 0.0;
    if (shared->lean_row_values_mode) {
      l_val[p] = lik;
    } else {
      *l_lu[p] = lik;
    }
    for (UF_long q = u_ptr[dep]; q < u_ptr[dep + 1u]; ++q) {
      const UF_long col = u_cols[q];
      x[col] = fma(-lik, u_val[q], x[col]);
    }
  }
  const double pivot = x[row];
  x[row] = 0.0;
  udiag[row] = pivot;
  for (UF_long q = u_ptr[row]; q < u_ptr[row + 1u]; ++q) {
    const UF_long col = u_cols[q];
    const double value = x[col];
    u_val[q] = value;
    if (!shared->lean_row_values_mode) {
      *u_lu[q] = value;
    }
    x[col] = 0.0;
  }
  if (pivot == 0.0) {
    kls_egraph_refactor_record_singular(shared, row, row);
  }
  return 1;
}

#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline, hot))
#endif
static void kls_compact_match_hoisted_worker_run(
  kls_egraph_refactor_worker *worker,
  unsigned int generation,
  const UF_long *restrict rows,
  UF_long stride) {
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  const UF_long n = solver->n;
  double *restrict x = worker->x;
  const double *restrict values = shared->values;
  const UF_long *restrict in_ptr = solver->row_refactor_input_ptr;
  const UF_long *restrict in_cols = solver->row_refactor_input_cols;
  const uint16_t *restrict in_cols16 = solver->row_refactor_input_cols16;
  const uint32_t *restrict in_packed =
    solver->row_refactor_input_col_direct32;
  const uint16_t *restrict in_user = solver->row_refactor_input_user_pos;
  const UF_long *restrict l_ptr = solver->row_refactor_l_ptr;
  const UF_long *restrict l_cols = solver->row_refactor_l_cols;
  double *restrict l_val = solver->row_refactor_l_row_values;
  const UF_long *restrict u_ptr = solver->row_refactor_u_ptr;
  const UF_long *restrict u_cols = solver->row_refactor_u_cols;
  double *restrict u_val = solver->row_refactor_u_row_values;
  double *restrict udiag = (double *)solver->numeric->Udiag;
  const uint32_t *restrict done_token =
    solver->lean_parallel_grouped_token;
  atomic_uint *restrict grouped_done = solver->lean_parallel_grouped_done;
  const unsigned int tid = (unsigned int)worker->tid;

  for (UF_long pos = (UF_long)worker->tid; pos < n; pos += stride) {
    const UF_long row = rows[pos];
    if (in_packed != NULL) {
      for (UF_long p = in_ptr[row]; p < in_ptr[row + 1u]; ++p) {
        const uint32_t packed = in_packed[p];
        x[packed & UINT32_C(0xffff)] = values[packed >> 16u];
      }
    } else if (in_cols16 != NULL) {
      for (UF_long p = in_ptr[row]; p < in_ptr[row + 1u]; ++p) {
        x[in_cols16[p]] = values[in_user[p]];
      }
    } else {
      for (UF_long p = in_ptr[row]; p < in_ptr[row + 1u]; ++p) {
        x[in_cols[p]] = values[in_user[p]];
      }
    }
    for (UF_long p = l_ptr[row]; p < l_ptr[row + 1u]; ++p) {
      const UF_long dep = l_cols[p];
      const uint32_t token = done_token[dep];
      if ((unsigned int)(token >> KLS_LEAN_GROUPED_OWNER_SHIFT) != tid) {
        const uint32_t slot = token & KLS_LEAN_GROUPED_SLOT_MASK;
        while (atomic_load_explicit(&grouped_done[slot],
                                    memory_order_acquire) != generation) {
          kls_cpu_relax();
        }
      }
      const double lik = x[dep] / udiag[dep];
      x[dep] = 0.0;
      l_val[p] = lik;
      for (UF_long q = u_ptr[dep]; q < u_ptr[dep + 1u]; ++q) {
        const UF_long col = u_cols[q];
        x[col] = fma(-lik, u_val[q], x[col]);
      }
    }
    const double pivot = x[row];
    x[row] = 0.0;
    udiag[row] = pivot;
    for (UF_long q = u_ptr[row]; q < u_ptr[row + 1u]; ++q) {
      const UF_long col = u_cols[q];
      const double value = x[col];
      u_val[q] = value;
      x[col] = 0.0;
    }
    if (pivot == 0.0) {
      kls_egraph_refactor_record_singular(shared, row, row);
    }
    const uint32_t slot = done_token[row] & KLS_LEAN_GROUPED_SLOT_MASK;
    atomic_store_explicit(&grouped_done[slot], generation,
                          memory_order_release);
  }

  double *restrict offx = (double *)solver->numeric->Offx;
  const uint16_t *restrict off_user =
    solver->compact_match_offdiag_user_pos;
  const UF_long offcount = solver->numeric->Offp[n];
  for (UF_long p = (UF_long)worker->tid; p < offcount; p += stride) {
    offx[p] = values[off_user[p]];
  }
}

static KLS_ALWAYS_INLINE int kls_lean_forward_initial_value(
  kls_egraph_refactor_shared *shared,
  kls_solver *solver,
  UF_long row,
  double *value) {
  if (shared == NULL || solver == NULL || value == NULL ||
      shared->lean_forward_solve_work == NULL) {
    if (value != NULL) *value = 0.0;
    return 1;
  }
  if (shared->lean_forward_solve_rhs != NULL) {
    if (solver->numeric == NULL || solver->numeric->Pnum == NULL) return 0;
    const UF_long source = solver->numeric->Pnum[row];
    if (source >= solver->n) return 0;
    *value = shared->lean_forward_solve_rhs[source];
  } else {
    *value = shared->lean_forward_solve_work[row];
  }
  return 1;
}

static KLS_ALWAYS_INLINE void
kls_symmetric_scalar_fringe_hoisted_worker_core(
  kls_egraph_refactor_worker *worker,
  unsigned int generation,
  const UF_long *restrict rows,
  UF_long stride,
  int avx512_target,
  int compact_cols,
  int wide_snodes) {
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  const UF_long n = solver->n;
  double *restrict x = worker->x;
  const double *restrict values = shared->values;
  const UF_long *restrict in_ptr = solver->row_refactor_input_ptr;
  const UF_long *restrict in_cols = solver->row_refactor_input_cols;
  const UF_long *restrict in_pos = solver->row_refactor_input_pos;
  const uint32_t *restrict in_packed =
    solver->row_refactor_input_col_user32;
  const UF_long *restrict l_ptr = solver->row_refactor_l_ptr;
  const UF_long *restrict l_cols = solver->row_refactor_l_cols;
  const uint16_t *restrict l_cols16 = solver->row_refactor_l_cols16;
  double *restrict l_val = solver->row_refactor_l_row_values;
  double **restrict l_lu = solver->row_refactor_l_values;
  const UF_long *restrict u_ptr = solver->row_refactor_u_ptr;
  const UF_long *restrict u_cols = solver->row_refactor_u_cols;
  const uint16_t *restrict u_cols16 = solver->row_refactor_u_cols16;
  double *restrict u_val = solver->row_refactor_u_row_values;
  double **restrict u_lu = solver->row_refactor_u_values;
  double *restrict udiag = (double *)solver->numeric->Udiag;
  const unsigned char *restrict snode_run = solver->lean_snode_run;
  const unsigned char *restrict snode_wait_mask =
    solver->lean_snode_wait_mask;
  const uint32_t *restrict snode_wait_slot =
    solver->lean_snode_wait_slot;
  const uint32_t *restrict done_token =
    solver->lean_parallel_grouped_token;
  atomic_uint *restrict grouped_done =
    solver->lean_parallel_grouped_done;
  const int row_values_mode = shared->lean_row_values_mode;
  double *restrict forward_work = shared->lean_forward_solve_work;
  const UF_long forward_begin = shared->lean_forward_solve_begin;

  for (UF_long pos = (UF_long)worker->tid; pos < n; pos += stride) {
    const UF_long row = rows[pos];
    const int forward_row = forward_work != NULL && row >= forward_begin;
    double forward_value = 0.0;
    if (forward_row &&
        !kls_lean_forward_initial_value(
          shared, solver, row, &forward_value)) {
      kls_egraph_refactor_record_invalid(shared);
      return;
    }
    if (in_packed != NULL) {
      for (UF_long p = in_ptr[row]; p < in_ptr[row + 1u]; ++p) {
        const uint32_t packed = in_packed[p];
        x[packed & UINT32_C(0xffff)] = values[packed >> 16u];
      }
    } else {
      for (UF_long p = in_ptr[row]; p < in_ptr[row + 1u]; ++p) {
        x[in_cols[p]] = values[in_pos[p]];
      }
    }
    for (UF_long p = l_ptr[row]; p < l_ptr[row + 1u];) {
      const UF_long run = snode_run != NULL ? (UF_long)snode_run[p] : 0u;
#if defined(__GNUC__) || defined(__clang__)
      if (wide_snodes && __builtin_expect(run == 16u, 0) &&
          p + 16u <= l_ptr[row + 1u]) {
#else
      if (wide_snodes && run == 16u && p + 16u <= l_ptr[row + 1u]) {
#endif
        if (!kls_lean_snode_wait_mask_ready(
              shared, grouped_done, generation,
              (unsigned int)snode_wait_mask[p], snode_wait_slot + p)) {
          return;
        }
        const UF_long dep0 = compact_cols
          ? (UF_long)l_cols16[p] : l_cols[p];
        const UF_long completed = kls_lean_row_refactor_snode_step_width16(
          p, dep0, x, row_values_mode ? l_val : NULL, l_lu,
          u_ptr, u_cols, u_cols16, u_val, udiag, compact_cols,
          avx512_target, shared, grouped_done, generation,
          (unsigned int)snode_wait_mask[p + 8u],
          snode_wait_slot + p + 8u);
        if (completed != 16u) {
          kls_egraph_refactor_record_invalid(shared);
          return;
        }
        if (forward_row) {
          for (UF_long local = 0u; local < 16u; ++local) {
            const double lij = row_values_mode
              ? l_val[p + local] : *l_lu[p + local];
            forward_value = fma(
              -lij, forward_work[dep0 + local], forward_value);
          }
        }
        p += 16u;
        continue;
      }
      if (run >= 2u && run <= 8u && p + run <= l_ptr[row + 1u]) {
        unsigned int mask = (unsigned int)snode_wait_mask[p];
        while (mask != 0u) {
#if defined(__GNUC__) || defined(__clang__)
          const unsigned int local = (unsigned int)__builtin_ctz(mask);
#else
          unsigned int local = 0u;
          while ((mask & (1u << local)) == 0u) {
            local++;
          }
#endif
          mask &= mask - 1u;
          const uint32_t slot = snode_wait_slot[p + local] - 1u;
          unsigned spin = 0u;
          while (atomic_load_explicit(&grouped_done[slot],
                                      memory_order_acquire) != generation) {
            if ((++spin & 255u) == 0u &&
                atomic_load_explicit(&shared->stop,
                                     memory_order_acquire) != 0) {
              return;
            }
            kls_cpu_relax();
          }
        }
        const UF_long dep0 = compact_cols
          ? (UF_long)l_cols16[p] : l_cols[p];
        const UF_long completed = compact_cols
          ? kls_lean_row_refactor_snode_step_i16_width_specialized(
              p, run, dep0, x, row_values_mode ? l_val : NULL, l_lu,
              u_ptr, u_cols16, u_val, udiag, avx512_target)
          : kls_lean_row_refactor_snode_step_width_specialized(
              p, run, dep0, x, row_values_mode ? l_val : NULL, l_lu,
              u_ptr, u_cols, u_val, udiag, avx512_target);
        if (completed != run) {
          kls_egraph_refactor_record_invalid(shared);
          return;
        }
        if (forward_row) {
          for (UF_long local = 0u; local < run; ++local) {
            const UF_long dep = dep0 + local;
            const double lij = row_values_mode
              ? l_val[p + local] : *l_lu[p + local];
            forward_value = fma(-lij, forward_work[dep], forward_value);
          }
        }
        p += run;
        continue;
      }
      const UF_long dep = compact_cols
        ? (UF_long)l_cols16[p] : l_cols[p];
      const uint32_t encoded_slot = snode_wait_slot[p];
      if (encoded_slot != 0u) {
        const uint32_t slot = encoded_slot - 1u;
        unsigned spin = 0u;
        while (atomic_load_explicit(&grouped_done[slot],
                                    memory_order_acquire) != generation) {
          if ((++spin & 255u) == 0u &&
              atomic_load_explicit(&shared->stop,
                                   memory_order_acquire) != 0) {
            return;
          }
          kls_cpu_relax();
        }
      }
      const double lik = x[dep] / udiag[dep];
      x[dep] = 0.0;
      if (row_values_mode) {
        l_val[p] = lik;
      } else {
        *l_lu[p] = lik;
      }
      if (forward_row) {
        forward_value = fma(-lik, forward_work[dep], forward_value);
      }
      if (lik != 0.0) {
        for (UF_long q = u_ptr[dep]; q < u_ptr[dep + 1u]; ++q) {
          const UF_long col = compact_cols
            ? (UF_long)u_cols16[q] : u_cols[q];
          x[col] = fma(-lik, u_val[q], x[col]);
        }
      }
      p++;
    }
    if (forward_row) {
      forward_work[row] = forward_value;
    }
    const double pivot = x[row];
    x[row] = 0.0;
    udiag[row] = pivot;
    const UF_long u_begin = u_ptr[row];
    const UF_long u_end = u_ptr[row + 1u];
    const uint32_t row_token = done_token[row];
    const int publishes =
      (row_token & KLS_LEAN_GROUPED_PUBLISH_FLAG) != 0u;
    if (row_values_mode || publishes) {
      /* Remote consumers and row-native solves read the contiguous mirror,
         never the packed-column destination.  Publish that dependency
         frontier as soon as its required values are ready; when the row
         representation remains authoritative, omit the packed scatter. */
      for (UF_long q = u_begin; q < u_end; ++q) {
        const UF_long col = compact_cols
          ? (UF_long)u_cols16[q] : u_cols[q];
        u_val[q] = x[col];
        x[col] = 0.0;
      }
    } else {
      for (UF_long q = u_begin; q < u_end; ++q) {
        const UF_long col = compact_cols
          ? (UF_long)u_cols16[q] : u_cols[q];
        const double value = x[col];
        u_val[q] = value;
        *u_lu[q] = value;
        x[col] = 0.0;
      }
    }
    if (pivot == 0.0) {
      kls_egraph_refactor_record_singular(shared, row, row);
    }
    if (publishes) {
      const uint32_t row_slot = row_token & KLS_LEAN_GROUPED_SLOT_MASK;
      atomic_store_explicit(&grouped_done[row_slot], generation,
                            memory_order_release);
    }
    if (publishes && !row_values_mode) {
      for (UF_long q = u_begin; q < u_end; ++q) {
        *u_lu[q] = u_val[q];
      }
    }
    if (atomic_load_explicit(&shared->stop, memory_order_acquire) != 0) {
      return;
    }
  }
}

#if defined(__GNUC__) || defined(__clang__)
#define KLS_SYMMETRIC_SNODE_WORKER_ATTR __attribute__((noinline, hot))
#else
#define KLS_SYMMETRIC_SNODE_WORKER_ATTR
#endif

KLS_SYMMETRIC_SNODE_WORKER_ATTR static void
kls_symmetric_scalar_fringe_hoisted_worker_run(
  kls_egraph_refactor_worker *worker,
  unsigned int generation,
  const UF_long *restrict rows,
  UF_long stride) {
  kls_symmetric_scalar_fringe_hoisted_worker_core(
    worker, generation, rows, stride, 0, 0, 0);
}

#if KLS_HAVE_AVX512_KERNELS
__attribute__((noinline, section(".text.kls_avx512"))) static void
kls_symmetric_scalar_fringe_hoisted_worker_run_avx512(
  kls_egraph_refactor_worker *worker,
  unsigned int generation,
  const UF_long *restrict rows,
  UF_long stride) {
  kls_symmetric_scalar_fringe_hoisted_worker_core(
    worker, generation, rows, stride, 1, 0, 0);
}
#endif

#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline, section(".text.kls_i16_snode")))
#endif
static void kls_symmetric_scalar_fringe_hoisted_worker_run_i16(
  kls_egraph_refactor_worker *worker,
  unsigned int generation,
  const UF_long *restrict rows,
  UF_long stride) {
  kls_symmetric_scalar_fringe_hoisted_worker_core(
    worker, generation, rows, stride, 0, 1, 0);
}

#if KLS_HAVE_AVX512_KERNELS
__attribute__((noinline, section(".text.kls_i16_snode"))) static void
kls_symmetric_scalar_fringe_hoisted_worker_run_i16_avx512(
  kls_egraph_refactor_worker *worker,
  unsigned int generation,
  const UF_long *restrict rows,
  UF_long stride) {
  kls_symmetric_scalar_fringe_hoisted_worker_core(
    worker, generation, rows, stride, 1, 1, 0);
}

/* Wide descriptors are selected structurally at setup and use dedicated
   workers.  Duplicating these two AVX-512 variants keeps the established
   workers' instruction layout independent of whether wide support is built. */
__attribute__((noinline, section(".kls_width16"))) static void
kls_symmetric_scalar_fringe_hoisted_worker_run_width16_avx512(
  kls_egraph_refactor_worker *worker,
  unsigned int generation,
  const UF_long *restrict rows,
  UF_long stride) {
  kls_symmetric_scalar_fringe_hoisted_worker_core(
    worker, generation, rows, stride, 1, 0, 1);
}

__attribute__((noinline, section(".kls_width16"))) static void
kls_symmetric_scalar_fringe_hoisted_worker_run_i16_width16_avx512(
  kls_egraph_refactor_worker *worker,
  unsigned int generation,
  const UF_long *restrict rows,
  UF_long stride) {
  kls_symmetric_scalar_fringe_hoisted_worker_core(
    worker, generation, rows, stride, 1, 1, 1);
}
#endif

#undef KLS_SYMMETRIC_SNODE_WORKER_ATTR

#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline, hot))
#endif
static void kls_generic_i16ptr_hoisted_worker_run(
  kls_egraph_refactor_worker *worker,
  unsigned int generation,
  const UF_long *restrict rows,
  const uint16_t *restrict rows16,
  UF_long stride) {
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  const UF_long n = solver->n;
  const unsigned int tid = (unsigned int)worker->tid;
  double *restrict x = worker->x;
  const double *restrict values = shared->values;
  const uint16_t *restrict in_ptr = solver->row_refactor_input_ptr16;
  const UF_long *restrict in_cols = solver->row_refactor_input_cols;
  const UF_long *restrict in_pos = solver->row_refactor_input_pos;
  const uint32_t *restrict in_packed =
    solver->row_refactor_input_col_user32;
  const uint16_t *restrict l_ptr = solver->row_refactor_l_ptr16;
  const uint16_t *restrict l_cols = solver->row_refactor_l_cols16;
  double **restrict l_lu = solver->row_refactor_l_values;
  double *restrict l_val = solver->row_refactor_l_row_values;
  const uint16_t *restrict u_ptr = solver->row_refactor_u_ptr16;
  const uint16_t *restrict u_cols = solver->row_refactor_u_cols16;
  double *restrict u_val = solver->row_refactor_u_row_values;
  double **restrict u_lu = solver->row_refactor_u_values;
  double *restrict udiag = (double *)solver->numeric->Udiag;
  kls_lean_done_slot *restrict done = shared->lean_done;
  const int row_values_mode = shared->lean_row_values_mode;
  const int grouped_done_mode = shared->lean_grouped_done_mode;
  const uint32_t *restrict done_token =
    solver->lean_parallel_grouped_token;
  atomic_uint *restrict grouped_done = solver->lean_parallel_grouped_done;
  double *restrict forward_work = shared->lean_forward_solve_work;
  for (UF_long pos = (UF_long)worker->tid; pos < n; pos += stride) {
    const UF_long row = rows16 != NULL
      ? (UF_long)rows16[pos] : rows[pos];
    if (row < solver->lean_scalar_btf_prefix) {
      continue;
    }
    double forward_value = 0.0;
    if (forward_work != NULL &&
        !kls_lean_forward_initial_value(
          shared, solver, row, &forward_value)) {
      kls_egraph_refactor_record_invalid(shared);
      return;
    }
    if (in_packed != NULL) {
      for (UF_long p = (UF_long)in_ptr[row];
           p < (UF_long)in_ptr[row + 1u]; ++p) {
        const uint32_t packed = in_packed[p];
        x[packed & UINT32_C(0xffff)] = values[packed >> 16u];
      }
    } else {
      for (UF_long p = (UF_long)in_ptr[row];
           p < (UF_long)in_ptr[row + 1u]; ++p) {
        x[in_cols[p]] = values[in_pos[p]];
      }
    }
    for (UF_long p = (UF_long)l_ptr[row];
         p < (UF_long)l_ptr[row + 1u]; ++p) {
      const UF_long dep = (UF_long)l_cols[p];
      if (grouped_done_mode) {
        const uint32_t token = done_token[dep];
        if ((unsigned int)(token >> KLS_LEAN_GROUPED_OWNER_SHIFT) != tid) {
          const uint32_t slot = token & KLS_LEAN_GROUPED_SLOT_MASK;
          while (atomic_load_explicit(&grouped_done[slot],
                                      memory_order_acquire) != generation) {
            kls_cpu_relax();
          }
        }
      } else if (done[dep].owner != tid) {
        unsigned spin = 0u;
        while (atomic_load_explicit(&done[dep].generation,
                                    memory_order_acquire) != generation) {
          if ((++spin & 255u) == 0u &&
              atomic_load_explicit(&shared->stop,
                                   memory_order_acquire) != 0) {
            return;
          }
          kls_cpu_relax();
        }
      }
      const double lik = x[dep] / udiag[dep];
      if (forward_work != NULL) {
        forward_value = fma(-lik, forward_work[dep], forward_value);
      }
      x[dep] = 0.0;
      if (row_values_mode) {
        l_val[p] = lik;
      } else {
        *l_lu[p] = lik;
      }
      if (lik != 0.0) {
        for (UF_long q = (UF_long)u_ptr[dep];
             q < (UF_long)u_ptr[dep + 1u]; ++q) {
          const UF_long col = (UF_long)u_cols[q];
          x[col] = fma(-lik, u_val[q], x[col]);
        }
      }
    }
    const double pivot = x[row];
    x[row] = 0.0;
    udiag[row] = pivot;
    for (UF_long q = (UF_long)u_ptr[row];
         q < (UF_long)u_ptr[row + 1u]; ++q) {
      const UF_long col = (UF_long)u_cols[q];
      const double value = x[col];
      u_val[q] = value;
      if (!row_values_mode) {
        *u_lu[q] = value;
      }
      x[col] = 0.0;
    }
    if (pivot == 0.0) {
      kls_egraph_refactor_record_singular(shared, row, row);
    }
    if (forward_work != NULL) forward_work[row] = forward_value;
    if (grouped_done_mode) {
      const uint32_t slot =
        done_token[row] & KLS_LEAN_GROUPED_SLOT_MASK;
      atomic_store_explicit(&grouped_done[slot], generation,
                            memory_order_release);
    } else {
      atomic_store_explicit(&done[row].generation, generation,
                            memory_order_release);
    }
    if (atomic_load_explicit(&shared->stop, memory_order_acquire) != 0) {
      return;
    }
  }
}

/* Generic unscaled packed-row worker.  Its execution stream removes random
   owner-token and U-pointer probes while retaining dependency-by-dependency
   overlap and the original arithmetic order. */
#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline, hot, aligned(64)))
#endif
static void kls_generic_packed_row_worker_run(
  kls_egraph_refactor_worker *worker,
  UF_long stride) {
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  const UF_long n = solver->n;
  double *restrict x = worker->x;
  const double *restrict values = shared->values;
  const uint32_t *restrict in_packed =
    solver->row_refactor_input_col_user32;
  const uint16_t *restrict in_cols16 = solver->row_refactor_input_cols16;
  const uint32_t *restrict in_pos32 = solver->row_refactor_input_pos32;
  double *restrict l_val = solver->row_refactor_l_row_values;
  double **restrict l_lu = solver->row_refactor_l_values;
  const uint16_t *restrict u_cols = solver->row_refactor_u_cols16;
  double *restrict u_val = solver->row_refactor_u_row_values;
  double **restrict u_lu = solver->row_refactor_u_values;
  double *restrict udiag = (double *)solver->numeric->Udiag;
  const kls_generic_packed_row_work32 *restrict work_rows =
    solver->lean_generic_packed_rows32;
  const kls_generic_packed_dependency32 *restrict dependencies =
    solver->lean_generic_packed_dependencies32;
  kls_lean_done_slot *restrict owner_frontier =
    solver->lean_parallel_owner_frontier;
  atomic_uint *restrict grouped_done = solver->lean_parallel_grouped_done;
  const uint32_t *restrict grouped_token = solver->lean_parallel_grouped_token;
  const unsigned char *restrict publish_mailbox =
    solver->lean_generic_packed_publish_mailbox;
  const uint32_t *restrict publish_begin =
    solver->lean_generic_packed_publish_begin;
  const unsigned int consumer = (unsigned int)worker->tid;
  const unsigned int sequence_bits =
    solver->lean_parallel_owner_frontier_sequence_bits;
  const unsigned int sequence_mask =
    solver->lean_parallel_owner_frontier_sequence_mask;
  const int row_values_mode = shared->lean_row_values_mode;
  const int mixed_input_mode = in_packed == NULL;
  double *restrict forward_work = shared->lean_forward_solve_work;
  const unsigned int generation = shared->pipeline_generation;
  unsigned int producer_sequence = 1u;
  UF_long publish_pos = (UF_long)publish_begin[consumer];
  for (UF_long pos = (UF_long)worker->tid; pos < n;
       pos += stride, ++producer_sequence) {
    const kls_generic_packed_row_work32 *restrict work = &work_rows[pos];
    const UF_long row = (UF_long)work->row;
    double forward_value = 0.0;
    if (forward_work != NULL &&
        !kls_lean_forward_initial_value(
          shared, solver, row, &forward_value)) {
      kls_egraph_refactor_record_invalid(shared);
      return;
    }
    UF_long p = (UF_long)work->input_begin;
    const UF_long input_end = (UF_long)work->input_end;
    if (in_packed != NULL) {
      for (; p + 3u < input_end; p += 4u) {
        const uint32_t packed0 = in_packed[p];
        const uint32_t packed1 = in_packed[p + 1u];
        const uint32_t packed2 = in_packed[p + 2u];
        const uint32_t packed3 = in_packed[p + 3u];
        x[packed0 & UINT32_C(0xffff)] = values[packed0 >> 16u];
        x[packed1 & UINT32_C(0xffff)] = values[packed1 >> 16u];
        x[packed2 & UINT32_C(0xffff)] = values[packed2 >> 16u];
        x[packed3 & UINT32_C(0xffff)] = values[packed3 >> 16u];
      }
      for (; p < input_end; ++p) {
        const uint32_t packed = in_packed[p];
        x[packed & UINT32_C(0xffff)] = values[packed >> 16u];
      }
    } else {
      for (; p + 3u < input_end; p += 4u) {
        x[in_cols16[p]] = values[in_pos32[p]];
        x[in_cols16[p + 1u]] = values[in_pos32[p + 1u]];
        x[in_cols16[p + 2u]] = values[in_pos32[p + 2u]];
        x[in_cols16[p + 3u]] = values[in_pos32[p + 3u]];
      }
      for (; p < input_end; ++p) x[in_cols16[p]] = values[in_pos32[p]];
    }
    p = (UF_long)work->dependency_begin;
    const UF_long dependency_end = (UF_long)work->dependency_end;
    for (; p < dependency_end; ++p) {
      const kls_generic_packed_dependency32 *restrict dependency =
        &dependencies[p];
      const uint32_t wait_frontier = dependency->wait_frontier;
      if (mixed_input_mode) {
        const uint32_t token = grouped_token[dependency->dep];
        if ((token >> KLS_LEAN_GROUPED_OWNER_SHIFT) != consumer) {
          const uint32_t slot = token & KLS_LEAN_GROUPED_SLOT_MASK;
          while (atomic_load_explicit(&grouped_done[slot],
                                      memory_order_acquire) != generation) {
            kls_cpu_relax();
          }
        }
      } else if (wait_frontier != 0u) {
        const unsigned int mailbox = wait_frontier >> sequence_bits;
        const unsigned int sequence = wait_frontier & sequence_mask;
        unsigned int spin = 0u;
        while (atomic_load_explicit(&owner_frontier[mailbox].generation,
                                    memory_order_acquire) < sequence) {
          if ((++spin & 255u) == 0u &&
              atomic_load_explicit(&shared->stop,
                                   memory_order_acquire) != 0) {
            return;
          }
          kls_cpu_relax();
        }
      }
      const UF_long dep = (UF_long)dependency->dep;
      const double lik = x[dep] / udiag[dep];
      if (forward_work != NULL) {
        forward_value = fma(-lik, forward_work[dep], forward_value);
      }
      x[dep] = 0.0;
      l_val[p] = lik;
      if (!row_values_mode) *l_lu[p] = lik;
      if (lik != 0.0) {
        for (UF_long q = (UF_long)dependency->update_begin;
             q < (UF_long)dependency->update_end; ++q) {
          const UF_long col = (UF_long)u_cols[q];
          x[col] = fma(-lik, u_val[q], x[col]);
        }
      }
    }
    const double pivot = x[row];
    x[row] = 0.0;
    udiag[row] = pivot;
    for (UF_long q = (UF_long)work->output_begin;
         q < (UF_long)work->output_end; ++q) {
      const UF_long col = (UF_long)u_cols[q];
      const double value = x[col];
      u_val[q] = value;
      if (!row_values_mode) *u_lu[q] = value;
      x[col] = 0.0;
    }
    if (pivot == 0.0) {
      kls_egraph_refactor_record_singular(shared, row, row);
    }
    if (forward_work != NULL) forward_work[row] = forward_value;
    if (mixed_input_mode) {
      const uint32_t slot = grouped_token[row] & KLS_LEAN_GROUPED_SLOT_MASK;
      atomic_store_explicit(&grouped_done[slot], generation,
                            memory_order_release);
    }
    const UF_long publish_end = (UF_long)work->publish_end;
    for (; publish_pos < publish_end; ++publish_pos) {
      const UF_long mailbox = (UF_long)publish_mailbox[publish_pos];
      atomic_store_explicit(&owner_frontier[mailbox].generation,
                            producer_sequence, memory_order_release);
    }
    if (atomic_load_explicit(&shared->stop, memory_order_acquire) != 0) {
      return;
    }
  }
}

#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline, hot, aligned(64)))
#endif
static void kls_scaled_generic_packed_row_worker_run(
  kls_egraph_refactor_worker *worker,
  UF_long stride) {
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  const UF_long n = solver->n;
  double *restrict x = worker->x;
  const double *restrict values = shared->values;
  double *restrict rs = (double *)(uintptr_t)shared->rs;
  const UF_long *restrict pnum = solver->numeric->Pnum;
  const uint32_t *restrict in_packed =
    solver->row_refactor_input_col_user32;
  const uint16_t *restrict in_cols16 = solver->row_refactor_input_cols16;
  const uint32_t *restrict in_pos32 = solver->row_refactor_input_pos32;
  double *restrict l_val = solver->row_refactor_l_row_values;
  double **restrict l_lu = solver->row_refactor_l_values;
  const uint16_t *restrict u_cols = solver->row_refactor_u_cols16;
  double *restrict u_val = solver->row_refactor_u_row_values;
  double **restrict u_lu = solver->row_refactor_u_values;
  double *restrict udiag = (double *)solver->numeric->Udiag;
  const kls_generic_packed_row_work32 *restrict work_rows =
    solver->lean_generic_packed_rows32;
  const kls_generic_packed_dependency32 *restrict dependencies =
    solver->lean_generic_packed_dependencies32;
  kls_lean_done_slot *restrict owner_frontier =
    solver->lean_parallel_owner_frontier;
  atomic_uint *restrict grouped_done = solver->lean_parallel_grouped_done;
  const uint32_t *restrict grouped_token = solver->lean_parallel_grouped_token;
  const unsigned char *restrict publish_mailbox =
    solver->lean_generic_packed_publish_mailbox;
  const uint32_t *restrict publish_begin =
    solver->lean_generic_packed_publish_begin;
  const unsigned int consumer = (unsigned int)worker->tid;
  const unsigned int sequence_bits =
    solver->lean_parallel_owner_frontier_sequence_bits;
  const unsigned int sequence_mask =
    solver->lean_parallel_owner_frontier_sequence_mask;
  const int fused_scale_mode = shared->lean_fused_scale_mode;
  const int sum_scale_mode = solver->common.scale == 1;
  const int row_values_mode = shared->lean_row_values_mode;
  const int mixed_input_mode = in_packed == NULL;
  const unsigned int generation = shared->pipeline_generation;
  unsigned int producer_sequence = 1u;
  UF_long publish_pos = (UF_long)publish_begin[consumer];
  for (UF_long pos = (UF_long)worker->tid; pos < n;
       pos += stride, ++producer_sequence) {
    const kls_generic_packed_row_work32 *restrict work = &work_rows[pos];
    const UF_long row = (UF_long)work->row;
    const UF_long input_begin = (UF_long)work->input_begin;
    const UF_long input_end = (UF_long)work->input_end;
    if (fused_scale_mode) {
      double row_scale = 0.0;
      for (UF_long p = input_begin; p < input_end; ++p) {
        const uint32_t col = in_packed != NULL
          ? in_packed[p] & UINT32_C(0xffff) : (uint32_t)in_cols16[p];
        const uint32_t source = in_packed != NULL
          ? in_packed[p] >> 16u : in_pos32[p];
        const double value = values[source];
        x[col] = value;
        const double magnitude = fabs(value);
        if (sum_scale_mode) {
          row_scale += magnitude;
        } else if (magnitude > row_scale) {
          row_scale = magnitude;
        }
      }
      if (row_scale == 0.0) {
        row_scale = 1.0;
      }
      rs[row] = row_scale;
      const double row_rs_inv = 1.0 / row_scale;
      for (UF_long p = input_begin; p < input_end; ++p) {
        const UF_long col = in_packed != NULL
          ? (UF_long)(in_packed[p] & UINT32_C(0xffff))
          : (UF_long)in_cols16[p];
        x[col] *= row_rs_inv;
      }
    } else {
      const double row_rs_inv = 1.0 / rs[pnum[row]];
      for (UF_long p = input_begin; p < input_end; ++p) {
        const uint32_t col = in_packed != NULL
          ? in_packed[p] & UINT32_C(0xffff) : (uint32_t)in_cols16[p];
        const uint32_t source = in_packed != NULL
          ? in_packed[p] >> 16u : in_pos32[p];
        x[col] = values[source] * row_rs_inv;
      }
    }
    for (UF_long p = (UF_long)work->dependency_begin;
         p < (UF_long)work->dependency_end; ++p) {
      const kls_generic_packed_dependency32 *restrict dependency =
        &dependencies[p];
      const uint32_t wait_frontier = dependency->wait_frontier;
      if (mixed_input_mode) {
        const uint32_t token = grouped_token[dependency->dep];
        if ((token >> KLS_LEAN_GROUPED_OWNER_SHIFT) != consumer) {
          const uint32_t slot = token & KLS_LEAN_GROUPED_SLOT_MASK;
          while (atomic_load_explicit(&grouped_done[slot],
                                      memory_order_acquire) != generation) {
            kls_cpu_relax();
          }
        }
      } else if (wait_frontier != 0u) {
        const unsigned int mailbox = wait_frontier >> sequence_bits;
        const unsigned int sequence = wait_frontier & sequence_mask;
        unsigned int spin = 0u;
        while (atomic_load_explicit(&owner_frontier[mailbox].generation,
                                    memory_order_acquire) < sequence) {
          if ((++spin & 255u) == 0u &&
              atomic_load_explicit(&shared->stop,
                                   memory_order_acquire) != 0) {
            return;
          }
          kls_cpu_relax();
        }
      }
      const UF_long dep = (UF_long)dependency->dep;
      const double lik = x[dep] / udiag[dep];
      x[dep] = 0.0;
      l_val[p] = lik;
      if (!row_values_mode) *l_lu[p] = lik;
      for (UF_long q = (UF_long)dependency->update_begin;
           q < (UF_long)dependency->update_end; ++q) {
        const UF_long col = (UF_long)u_cols[q];
        x[col] = fma(-lik, u_val[q], x[col]);
      }
    }
    const double pivot = x[row];
    x[row] = 0.0;
    udiag[row] = pivot;
    for (UF_long q = (UF_long)work->output_begin;
         q < (UF_long)work->output_end; ++q) {
      const UF_long col = (UF_long)u_cols[q];
      const double value = x[col];
      u_val[q] = value;
      if (!row_values_mode) *u_lu[q] = value;
      x[col] = 0.0;
    }
    if (pivot == 0.0) {
      kls_egraph_refactor_record_singular(shared, row, row);
    }
    if (mixed_input_mode) {
      const uint32_t slot = grouped_token[row] & KLS_LEAN_GROUPED_SLOT_MASK;
      atomic_store_explicit(&grouped_done[slot], generation,
                            memory_order_release);
    }
    const UF_long publish_end = (UF_long)work->publish_end;
    for (; publish_pos < publish_end; ++publish_pos) {
      const UF_long mailbox = (UF_long)publish_mailbox[publish_pos];
      atomic_store_explicit(&owner_frontier[mailbox].generation,
                            producer_sequence, memory_order_release);
    }
    if (atomic_load_explicit(&shared->stop, memory_order_acquire) != 0) {
      return;
    }
  }
}

/* Exact unscaled compact-row worker for the common grouped/row-value
   representation.  The preceding numeric's diagonal reciprocals are valid
   until their owning rows publish replacements, so consume that dependency
   stream with multiplies and form each new reciprocal while the pivot is
   hot.  This also removes the serial O(n) reciprocal refresh before the
   following compact solve.  Eligibility below proves every representation
   invariant; the general i16 worker retains all fallback combinations. */
#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline, hot, aligned(64)))
#endif
static void kls_generic_i16ptr_grouped_recip_worker_run(
  kls_egraph_refactor_worker *worker,
  unsigned int generation,
  const UF_long *restrict rows,
  const uint16_t *restrict rows16,
  UF_long stride) {
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  const UF_long n = solver->n;
  const unsigned int tid = (unsigned int)worker->tid;
  double *restrict x = worker->x;
  const double *restrict values = shared->values;
  const uint16_t *restrict in_ptr = solver->row_refactor_input_ptr16;
  const UF_long *restrict in_cols = solver->row_refactor_input_cols;
  const UF_long *restrict in_pos = solver->row_refactor_input_pos;
  const uint32_t *restrict in_packed =
    solver->row_refactor_input_col_user32;
  const uint16_t *restrict l_ptr = solver->row_refactor_l_ptr16;
  const uint16_t *restrict l_cols = solver->row_refactor_l_cols16;
  double *restrict l_val = solver->row_refactor_l_row_values;
  const uint16_t *restrict u_ptr = solver->row_refactor_u_ptr16;
  const uint16_t *restrict u_cols = solver->row_refactor_u_cols16;
  double *restrict u_val = solver->row_refactor_u_row_values;
  double *restrict udiag = (double *)solver->numeric->Udiag;
  double *restrict udiag_recip = shared->udiag_recip;
  const uint32_t *restrict done_token =
    solver->lean_parallel_grouped_token;
  atomic_uint *restrict grouped_done = solver->lean_parallel_grouped_done;
  double *restrict forward_work = shared->lean_forward_solve_work;

  for (UF_long pos = (UF_long)worker->tid; pos < n; pos += stride) {
    const UF_long row = rows16 != NULL
      ? (UF_long)rows16[pos] : rows[pos];
    if (row < solver->lean_scalar_btf_prefix) {
      continue;
    }
    double forward_value = 0.0;
    if (forward_work != NULL &&
        !kls_lean_forward_initial_value(
          shared, solver, row, &forward_value)) {
      kls_egraph_refactor_record_invalid(shared);
      return;
    }
    if (in_packed != NULL) {
      UF_long p = (UF_long)in_ptr[row];
      const UF_long input_end = (UF_long)in_ptr[row + 1u];
      for (; p + 3u < input_end; p += 4u) {
        const uint32_t packed0 = in_packed[p];
        const uint32_t packed1 = in_packed[p + 1u];
        const uint32_t packed2 = in_packed[p + 2u];
        const uint32_t packed3 = in_packed[p + 3u];
        x[packed0 & UINT32_C(0xffff)] = values[packed0 >> 16u];
        x[packed1 & UINT32_C(0xffff)] = values[packed1 >> 16u];
        x[packed2 & UINT32_C(0xffff)] = values[packed2 >> 16u];
        x[packed3 & UINT32_C(0xffff)] = values[packed3 >> 16u];
      }
      for (; p < input_end; ++p) {
        const uint32_t packed = in_packed[p];
        x[packed & UINT32_C(0xffff)] = values[packed >> 16u];
      }
    } else {
      UF_long p = (UF_long)in_ptr[row];
      const UF_long input_end = (UF_long)in_ptr[row + 1u];
      for (; p + 3u < input_end; p += 4u) {
        x[in_cols[p]] = values[in_pos[p]];
        x[in_cols[p + 1u]] = values[in_pos[p + 1u]];
        x[in_cols[p + 2u]] = values[in_pos[p + 2u]];
        x[in_cols[p + 3u]] = values[in_pos[p + 3u]];
      }
      for (; p < input_end; ++p) {
        x[in_cols[p]] = values[in_pos[p]];
      }
    }
    for (UF_long p = (UF_long)l_ptr[row];
         p < (UF_long)l_ptr[row + 1u]; ++p) {
      const UF_long dep = (UF_long)l_cols[p];
      const uint32_t token = done_token[dep];
      if ((unsigned int)(token >> KLS_LEAN_GROUPED_OWNER_SHIFT) != tid) {
        const uint32_t slot = token & KLS_LEAN_GROUPED_SLOT_MASK;
        while (atomic_load_explicit(&grouped_done[slot],
                                    memory_order_acquire) != generation) {
          kls_cpu_relax();
        }
      }
      const double lik = x[dep] * udiag_recip[dep];
      if (forward_work != NULL) {
        forward_value = fma(-lik, forward_work[dep], forward_value);
      }
      x[dep] = 0.0;
      l_val[p] = lik;
      if (lik != 0.0) {
        for (UF_long q = (UF_long)u_ptr[dep];
             q < (UF_long)u_ptr[dep + 1u]; ++q) {
          const UF_long col = (UF_long)u_cols[q];
          x[col] = fma(-lik, u_val[q], x[col]);
        }
      }
    }
    const double pivot = x[row];
    x[row] = 0.0;
    udiag[row] = pivot;
    UF_long q = (UF_long)u_ptr[row];
    const UF_long output_end = (UF_long)u_ptr[row + 1u];
    for (; q + 3u < output_end; q += 4u) {
      const UF_long col0 = (UF_long)u_cols[q];
      const UF_long col1 = (UF_long)u_cols[q + 1u];
      const UF_long col2 = (UF_long)u_cols[q + 2u];
      const UF_long col3 = (UF_long)u_cols[q + 3u];
      const double value0 = x[col0];
      const double value1 = x[col1];
      const double value2 = x[col2];
      const double value3 = x[col3];
      u_val[q] = value0;
      u_val[q + 1u] = value1;
      u_val[q + 2u] = value2;
      u_val[q + 3u] = value3;
      x[col0] = 0.0;
      x[col1] = 0.0;
      x[col2] = 0.0;
      x[col3] = 0.0;
    }
    for (; q < output_end; ++q) {
      const UF_long col = (UF_long)u_cols[q];
      const double value = x[col];
      u_val[q] = value;
      x[col] = 0.0;
    }
    if (pivot == 0.0) {
      udiag_recip[row] = 0.0;
      kls_egraph_refactor_record_singular(shared, row, row);
    } else {
      udiag_recip[row] = 1.0 / pivot;
    }
    if (forward_work != NULL) forward_work[row] = forward_value;
    const uint32_t slot =
      done_token[row] & KLS_LEAN_GROUPED_SLOT_MASK;
    atomic_store_explicit(&grouped_done[slot], generation,
                          memory_order_release);
  }
}

/* A retained compact AMF two-block factor can enter the generic worker with
   compact pointers/columns, row-value mirrors, grouped completion, an
   unscaled input, and a dependency-aware static schedule.  Hoist those
   verified representation invariants out of the numeric walk.  The general
   worker retains every fallback; this loop preserves its scatter,
   dependency, arithmetic, publication, and stop-check order exactly. */
#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline, hot))
#endif
static void kls_compact_amf_two_block_i16ptr_worker_run(
  kls_egraph_refactor_worker *worker,
  unsigned int generation,
  const UF_long *restrict rows,
  UF_long stride) {
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  const UF_long n = solver->n;
  const unsigned int tid = (unsigned int)worker->tid;
  double *restrict x = worker->x;
  const double *restrict values = shared->values;
  const uint16_t *restrict in_ptr = solver->row_refactor_input_ptr16;
  const uint32_t *restrict in_packed =
    solver->row_refactor_input_col_user32;
  const UF_long *restrict in_cols = solver->row_refactor_input_cols;
  const UF_long *restrict in_pos = solver->row_refactor_input_pos;
  const uint16_t *restrict l_ptr = solver->row_refactor_l_ptr16;
  const uint16_t *restrict l_cols = solver->row_refactor_l_cols16;
  double *restrict l_val = solver->row_refactor_l_row_values;
  const uint16_t *restrict u_ptr = solver->row_refactor_u_ptr16;
  const uint16_t *restrict u_cols = solver->row_refactor_u_cols16;
  double *restrict u_val = solver->row_refactor_u_row_values;
  double *restrict udiag = (double *)solver->numeric->Udiag;
  const uint32_t *restrict done_token =
    solver->lean_parallel_grouped_token;
  atomic_uint *restrict grouped_done = solver->lean_parallel_grouped_done;
  for (UF_long pos = (UF_long)worker->tid; pos < n; pos += stride) {
    const UF_long row = rows[pos];
    if (in_packed != NULL) {
      for (UF_long p = (UF_long)in_ptr[row];
           p < (UF_long)in_ptr[row + 1u]; ++p) {
        const uint32_t packed = in_packed[p];
        x[packed & UINT32_C(0xffff)] = values[packed >> 16u];
      }
    } else {
      for (UF_long p = (UF_long)in_ptr[row];
           p < (UF_long)in_ptr[row + 1u]; ++p) {
        x[in_cols[p]] = values[in_pos[p]];
      }
    }
    for (UF_long p = (UF_long)l_ptr[row];
         p < (UF_long)l_ptr[row + 1u]; ++p) {
      const UF_long dep = (UF_long)l_cols[p];
      const uint32_t token = done_token[dep];
      if ((unsigned int)(token >> KLS_LEAN_GROUPED_OWNER_SHIFT) != tid) {
        const uint32_t slot = token & KLS_LEAN_GROUPED_SLOT_MASK;
        while (atomic_load_explicit(&grouped_done[slot],
                                    memory_order_acquire) != generation) {
          kls_cpu_relax();
        }
      }
      const double lik = x[dep] / udiag[dep];
      x[dep] = 0.0;
      l_val[p] = lik;
      if (lik != 0.0) {
        for (UF_long q = (UF_long)u_ptr[dep];
             q < (UF_long)u_ptr[dep + 1u]; ++q) {
          const UF_long col = (UF_long)u_cols[q];
          x[col] = fma(-lik, u_val[q], x[col]);
        }
      }
    }
    const double pivot = x[row];
    x[row] = 0.0;
    udiag[row] = pivot;
    for (UF_long q = (UF_long)u_ptr[row];
         q < (UF_long)u_ptr[row + 1u]; ++q) {
      const UF_long col = (UF_long)u_cols[q];
      const double value = x[col];
      u_val[q] = value;
      x[col] = 0.0;
    }
    if (pivot == 0.0) {
      kls_egraph_refactor_record_singular(shared, row, row);
    }
    const uint32_t slot =
      done_token[row] & KLS_LEAN_GROUPED_SLOT_MASK;
    atomic_store_explicit(&grouped_done[slot], generation,
                          memory_order_release);
    if (atomic_load_explicit(&shared->stop, memory_order_acquire) != 0) {
      return;
    }
  }
}

/* Predecoded waits preserve the original dependency-by-dependency overlap,
   while replacing random token-owner loads with one sequential packed word.
   A zero upper half means an earlier acquire on this consumer stream already
   covers the dependency. */
static KLS_ALWAYS_INLINE int kls_compact_amf_two_block_wait_for_dependency(
  kls_egraph_refactor_shared *shared,
  const kls_lean_done_slot *owner_frontier,
  unsigned int frontier_sequence_bits,
  unsigned int frontier_sequence_mask,
  uint64_t descriptor) {
  const uint32_t encoded_frontier =
    (uint32_t)((descriptor >> 16u) & UINT64_C(0xffff));
#if defined(__GNUC__) || defined(__clang__)
  const int wait_required = __builtin_expect(encoded_frontier != 0u, 0);
#else
  const int wait_required = encoded_frontier != 0u;
#endif
  if (!wait_required) {
    return 1;
  }
  const unsigned int mailbox =
    encoded_frontier >> frontier_sequence_bits;
  const unsigned int sequence =
    encoded_frontier & frontier_sequence_mask;
  unsigned int spin = 0u;
  while (atomic_load_explicit(&owner_frontier[mailbox].generation,
                              memory_order_acquire) < sequence) {
    if ((++spin & 255u) == 0u &&
        atomic_load_explicit(&shared->stop, memory_order_acquire) != 0) {
      return 0;
    }
    kls_cpu_relax();
  }
  return 1;
}

static KLS_ALWAYS_INLINE void
kls_compact_amf_two_block_update_sparse_row_pair(
  double *restrict x,
  const uint16_t *restrict u_cols,
  const double *restrict u_val,
  UF_long begin,
  UF_long end,
  double multiplier) {
  UF_long q = begin;
  for (; q + 1u < end; q += 2u) {
    const UF_long col0 = (UF_long)u_cols[q];
    const UF_long col1 = (UF_long)u_cols[q + 1u];
    const double value0 = x[col0];
    const double value1 = x[col1];
    x[col0] = fma(-multiplier, u_val[q], value0);
    x[col1] = fma(-multiplier, u_val[q + 1u], value1);
  }
  if (q < end) {
    const UF_long col = (UF_long)u_cols[q];
    x[col] = fma(-multiplier, u_val[q], x[col]);
  }
}

#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline, hot, aligned(64)))
#endif
static void kls_compact_amf_two_block_packed_exact_worker_run(
  kls_egraph_refactor_worker *worker,
  UF_long stride) {
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  const UF_long n = solver->n;
  double *restrict x = worker->x;
  const double *restrict values = shared->values;
  const uint32_t *restrict in_packed =
    solver->row_refactor_input_col_user32;
  const uint64_t *restrict l_dep_work =
    solver->lean_parallel_l_dep_work64;
  double *restrict l_val = solver->row_refactor_l_row_values;
  const uint16_t *restrict u_cols = solver->row_refactor_u_cols16;
  double *restrict u_val = solver->row_refactor_u_row_values;
  double *restrict udiag = (double *)solver->numeric->Udiag;
  double *restrict udiag_inv = solver->i32solve_udiag_recip != NULL
    ? solver->i32solve_udiag_recip : solver->lean_parallel_udiag_inv;
  const kls_compact_amf_two_block_work_row16 *restrict work_rows =
    solver->lean_parallel_work_rows16;
  kls_lean_done_slot *restrict owner_frontier =
    solver->lean_parallel_owner_frontier;
  const unsigned char *restrict publish_mailbox =
    solver->lean_parallel_publish_mailbox;
  const uint16_t *restrict publish_begin =
    solver->lean_parallel_publish_begin;
  const unsigned int consumer = (unsigned int)worker->tid;
  const unsigned int frontier_sequence_bits =
    solver->lean_parallel_owner_frontier_sequence_bits;
  const unsigned int frontier_sequence_mask =
    solver->lean_parallel_owner_frontier_sequence_mask;
  unsigned int producer_sequence = 1u;
  UF_long publish_pos = (UF_long)publish_begin[consumer];
  for (UF_long pos = (UF_long)worker->tid; pos < n;
       pos += stride, ++producer_sequence) {
    const kls_compact_amf_two_block_work_row16 *restrict work = &work_rows[pos];
    const UF_long row = (UF_long)work->row;
    UF_long p = (UF_long)work->input_begin;
    const UF_long input_end = (UF_long)work->input_end;
    for (; p + 3u < input_end; p += 4u) {
      const uint32_t packed0 = in_packed[p];
      const uint32_t packed1 = in_packed[p + 1u];
      const uint32_t packed2 = in_packed[p + 2u];
      const uint32_t packed3 = in_packed[p + 3u];
      x[packed0 & UINT32_C(0xffff)] = values[packed0 >> 16u];
      x[packed1 & UINT32_C(0xffff)] = values[packed1 >> 16u];
      x[packed2 & UINT32_C(0xffff)] = values[packed2 >> 16u];
      x[packed3 & UINT32_C(0xffff)] = values[packed3 >> 16u];
    }
    for (; p < input_end; ++p) {
      const uint32_t packed = in_packed[p];
      x[packed & UINT32_C(0xffff)] = values[packed >> 16u];
    }
    p = (UF_long)work->dependency_begin;
    const UF_long dependency_end = (UF_long)work->dependency_end;
    while (p < dependency_end) {
      const uint64_t packed = l_dep_work[p];
      const UF_long dep = (UF_long)(packed & UINT64_C(0x0fff));
      const uint32_t run =
        (uint32_t)((packed >> 12u) & UINT64_C(0x0f));
      if (run < 2u) {
        if (!kls_compact_amf_two_block_wait_for_dependency(
              shared, owner_frontier, frontier_sequence_bits,
              frontier_sequence_mask, packed)) {
          return;
        }
        const double lik = x[dep] * udiag_inv[dep];
        x[dep] = 0.0;
        l_val[p] = lik;
        kls_compact_amf_two_block_update_sparse_row_pair(
          x, u_cols, u_val,
          (UF_long)((packed >> 32u) & UINT64_C(0xffff)),
          (UF_long)(packed >> 48u), lik);
        p++;
        continue;
      }
      if (run == 2u && p + 1u < dependency_end) {
        const uint64_t packed1 = l_dep_work[p + 1u];
        if (!kls_compact_amf_two_block_wait_for_dependency(
              shared, owner_frontier, frontier_sequence_bits,
              frontier_sequence_mask, packed) ||
            !kls_compact_amf_two_block_wait_for_dependency(
              shared, owner_frontier, frontier_sequence_bits,
              frontier_sequence_mask, packed1)) {
          return;
        }
        const UF_long dep1 =
          (UF_long)(packed1 & UINT64_C(0x0fff));
        const UF_long update0 =
          (UF_long)((packed >> 32u) & UINT64_C(0xffff));
        const UF_long update1 =
          (UF_long)((packed1 >> 32u) & UINT64_C(0xffff));
        const UF_long update1_end = (UF_long)(packed1 >> 48u);
        const double multiplier0 = x[dep] * udiag_inv[dep];
        const double candidate1 =
          fma(-multiplier0, u_val[update0], x[dep1]);
        const double multiplier1 = candidate1 * udiag_inv[dep1];
        x[dep] = 0.0;
        x[dep1] = 0.0;
        l_val[p] = multiplier0;
        l_val[p + 1u] = multiplier1;
        for (UF_long q = update1; q < update1_end; ++q) {
          const UF_long col = (UF_long)u_cols[q];
          const UF_long offset = q - update1;
          const double after0 =
            fma(-multiplier0, u_val[update0 + 1u + offset], x[col]);
          x[col] = fma(-multiplier1, u_val[q], after0);
        }
        p += 2u;
        continue;
      }
      if (run == 3u && p + 2u < dependency_end) {
        const uint64_t packed1 = l_dep_work[p + 1u];
        const uint64_t packed2 = l_dep_work[p + 2u];
        if (!kls_compact_amf_two_block_wait_for_dependency(
              shared, owner_frontier, frontier_sequence_bits,
              frontier_sequence_mask, packed) ||
            !kls_compact_amf_two_block_wait_for_dependency(
              shared, owner_frontier, frontier_sequence_bits,
              frontier_sequence_mask, packed1) ||
            !kls_compact_amf_two_block_wait_for_dependency(
              shared, owner_frontier, frontier_sequence_bits,
              frontier_sequence_mask, packed2)) {
          return;
        }
        const UF_long dep1 =
          (UF_long)(packed1 & UINT64_C(0x0fff));
        const UF_long dep2 =
          (UF_long)(packed2 & UINT64_C(0x0fff));
        const UF_long update0 =
          (UF_long)((packed >> 32u) & UINT64_C(0xffff));
        const UF_long update1 =
          (UF_long)((packed1 >> 32u) & UINT64_C(0xffff));
        const UF_long update2 =
          (UF_long)((packed2 >> 32u) & UINT64_C(0xffff));
        const UF_long update2_end = (UF_long)(packed2 >> 48u);
        const double multiplier0 = x[dep] * udiag_inv[dep];
        const double candidate1 =
          fma(-multiplier0, u_val[update0], x[dep1]);
        const double multiplier1 = candidate1 * udiag_inv[dep1];
        double candidate2 =
          fma(-multiplier0, u_val[update0 + 1u], x[dep2]);
        candidate2 = fma(-multiplier1, u_val[update1], candidate2);
        const double multiplier2 = candidate2 * udiag_inv[dep2];
        x[dep] = 0.0;
        x[dep1] = 0.0;
        x[dep2] = 0.0;
        l_val[p] = multiplier0;
        l_val[p + 1u] = multiplier1;
        l_val[p + 2u] = multiplier2;
        for (UF_long q = update2; q < update2_end; ++q) {
          const UF_long col = (UF_long)u_cols[q];
          const UF_long offset = q - update2;
          double value = fma(-multiplier0,
                             u_val[update0 + 2u + offset], x[col]);
          value = fma(-multiplier1,
                      u_val[update1 + 1u + offset], value);
          x[col] = fma(-multiplier2, u_val[q], value);
        }
        p += 3u;
        continue;
      }
      if (run == 4u && p + 3u < dependency_end) {
        const uint64_t packed1 = l_dep_work[p + 1u];
        const uint64_t packed2 = l_dep_work[p + 2u];
        const uint64_t packed3 = l_dep_work[p + 3u];
        if (!kls_compact_amf_two_block_wait_for_dependency(
              shared, owner_frontier, frontier_sequence_bits,
              frontier_sequence_mask, packed) ||
            !kls_compact_amf_two_block_wait_for_dependency(
              shared, owner_frontier, frontier_sequence_bits,
              frontier_sequence_mask, packed1) ||
            !kls_compact_amf_two_block_wait_for_dependency(
              shared, owner_frontier, frontier_sequence_bits,
              frontier_sequence_mask, packed2) ||
            !kls_compact_amf_two_block_wait_for_dependency(
              shared, owner_frontier, frontier_sequence_bits,
              frontier_sequence_mask, packed3)) {
          return;
        }
        const UF_long dep1 =
          (UF_long)(packed1 & UINT64_C(0x0fff));
        const UF_long dep2 =
          (UF_long)(packed2 & UINT64_C(0x0fff));
        const UF_long dep3 =
          (UF_long)(packed3 & UINT64_C(0x0fff));
        const UF_long update0 =
          (UF_long)((packed >> 32u) & UINT64_C(0xffff));
        const UF_long update1 =
          (UF_long)((packed1 >> 32u) & UINT64_C(0xffff));
        const UF_long update2 =
          (UF_long)((packed2 >> 32u) & UINT64_C(0xffff));
        const UF_long update3 =
          (UF_long)((packed3 >> 32u) & UINT64_C(0xffff));
        const UF_long update3_end = (UF_long)(packed3 >> 48u);
        const double multiplier0 = x[dep] * udiag_inv[dep];
        const double candidate1 =
          fma(-multiplier0, u_val[update0], x[dep1]);
        const double multiplier1 = candidate1 * udiag_inv[dep1];
        double candidate2 =
          fma(-multiplier0, u_val[update0 + 1u], x[dep2]);
        candidate2 = fma(-multiplier1, u_val[update1], candidate2);
        const double multiplier2 = candidate2 * udiag_inv[dep2];
        double candidate3 =
          fma(-multiplier0, u_val[update0 + 2u], x[dep3]);
        candidate3 =
          fma(-multiplier1, u_val[update1 + 1u], candidate3);
        candidate3 = fma(-multiplier2, u_val[update2], candidate3);
        const double multiplier3 = candidate3 * udiag_inv[dep3];
        x[dep] = 0.0;
        x[dep1] = 0.0;
        x[dep2] = 0.0;
        x[dep3] = 0.0;
        l_val[p] = multiplier0;
        l_val[p + 1u] = multiplier1;
        l_val[p + 2u] = multiplier2;
        l_val[p + 3u] = multiplier3;
        for (UF_long q = update3; q < update3_end; ++q) {
          const UF_long col = (UF_long)u_cols[q];
          const UF_long offset = q - update3;
          double value = fma(-multiplier0,
                             u_val[update0 + 3u + offset], x[col]);
          value = fma(-multiplier1,
                      u_val[update1 + 2u + offset], value);
          value = fma(-multiplier2,
                      u_val[update2 + 1u + offset], value);
          x[col] = fma(-multiplier3, u_val[q], value);
        }
        p += 4u;
        continue;
      }
      if (!kls_compact_amf_two_block_wait_for_dependency(
            shared, owner_frontier, frontier_sequence_bits,
            frontier_sequence_mask, packed)) {
        return;
      }
      const double lik = x[dep] * udiag_inv[dep];
      x[dep] = 0.0;
      l_val[p] = lik;
      kls_compact_amf_two_block_update_sparse_row_pair(
        x, u_cols, u_val,
        (UF_long)((packed >> 32u) & UINT64_C(0xffff)),
        (UF_long)(packed >> 48u), lik);
      p++;
    }
    const double pivot = x[row];
    x[row] = 0.0;
    udiag[row] = pivot;
    UF_long q = (UF_long)work->output_begin;
    const UF_long output_end = (UF_long)work->output_end;
    for (; q + 3u < output_end; q += 4u) {
      const UF_long col0 = (UF_long)u_cols[q];
      const UF_long col1 = (UF_long)u_cols[q + 1u];
      const UF_long col2 = (UF_long)u_cols[q + 2u];
      const UF_long col3 = (UF_long)u_cols[q + 3u];
      const double value0 = x[col0];
      const double value1 = x[col1];
      const double value2 = x[col2];
      const double value3 = x[col3];
      u_val[q] = value0;
      u_val[q + 1u] = value1;
      u_val[q + 2u] = value2;
      u_val[q + 3u] = value3;
      x[col0] = 0.0;
      x[col1] = 0.0;
      x[col2] = 0.0;
      x[col3] = 0.0;
    }
    for (; q < output_end; ++q) {
      const UF_long col = (UF_long)u_cols[q];
      const double value = x[col];
      u_val[q] = value;
      x[col] = 0.0;
    }
    if (pivot == 0.0) {
      udiag_inv[row] = 0.0;
      kls_egraph_refactor_record_singular(shared, row, row);
    } else {
      udiag_inv[row] = 1.0 / pivot;
    }
    const UF_long publish_end = (UF_long)work->publish_end;
    for (; publish_pos < publish_end; ++publish_pos) {
      const UF_long mailbox = (UF_long)publish_mailbox[publish_pos];
      atomic_store_explicit(&owner_frontier[mailbox].generation,
                            producer_sequence, memory_order_release);
    }
  }
}

#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline, hot))
#endif
static void kls_generic_i16_hoisted_worker_run(
  kls_egraph_refactor_worker *worker,
  unsigned int generation,
  const UF_long *restrict rows,
  UF_long stride) {
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  const UF_long n = solver->n;
  const unsigned int tid = (unsigned int)worker->tid;
  double *restrict x = worker->x;
  const double *restrict values = shared->values;
  const UF_long *restrict in_ptr = solver->row_refactor_input_ptr;
  const UF_long *restrict in_cols = solver->row_refactor_input_cols;
  const UF_long *restrict in_pos = solver->row_refactor_input_pos;
  const uint32_t *restrict in_packed =
    solver->row_refactor_input_col_user32;
  const UF_long *restrict l_ptr = solver->row_refactor_l_ptr;
  const uint16_t *restrict l_cols = solver->row_refactor_l_cols16;
  double **restrict l_lu = solver->row_refactor_l_values;
  double *restrict l_val = solver->row_refactor_l_row_values;
  const UF_long *restrict u_ptr = solver->row_refactor_u_ptr;
  const uint16_t *restrict u_cols = solver->row_refactor_u_cols16;
  double *restrict u_val = solver->row_refactor_u_row_values;
  double **restrict u_lu = solver->row_refactor_u_values;
  double *restrict udiag = (double *)solver->numeric->Udiag;
  kls_lean_done_slot *restrict done = shared->lean_done;
  const int row_values_mode = shared->lean_row_values_mode;
  const int grouped_done_mode = shared->lean_grouped_done_mode;
  const uint32_t *restrict done_token =
    solver->lean_parallel_grouped_token;
  atomic_uint *restrict grouped_done = solver->lean_parallel_grouped_done;
  double *restrict forward_work = shared->lean_forward_solve_work;

  for (UF_long pos = (UF_long)worker->tid; pos < n; pos += stride) {
    const UF_long row = rows[pos];
    if (row < solver->lean_scalar_btf_prefix) {
      continue;
    }
    double forward_value = 0.0;
    if (forward_work != NULL &&
        !kls_lean_forward_initial_value(
          shared, solver, row, &forward_value)) {
      kls_egraph_refactor_record_invalid(shared);
      return;
    }
    if (in_packed != NULL) {
      for (UF_long p = in_ptr[row]; p < in_ptr[row + 1u]; ++p) {
        const uint32_t packed = in_packed[p];
        x[packed & UINT32_C(0xffff)] = values[packed >> 16u];
      }
    } else {
      for (UF_long p = in_ptr[row]; p < in_ptr[row + 1u]; ++p) {
        x[in_cols[p]] = values[in_pos[p]];
      }
    }
    for (UF_long p = l_ptr[row]; p < l_ptr[row + 1u]; ++p) {
      const UF_long dep = (UF_long)l_cols[p];
      if (grouped_done_mode) {
        const uint32_t token = done_token[dep];
        if ((unsigned int)(token >> KLS_LEAN_GROUPED_OWNER_SHIFT) != tid) {
          const uint32_t slot = token & KLS_LEAN_GROUPED_SLOT_MASK;
          while (atomic_load_explicit(&grouped_done[slot],
                                      memory_order_acquire) != generation) {
            kls_cpu_relax();
          }
        }
      } else if (done[dep].owner != tid) {
        unsigned spin = 0u;
        while (atomic_load_explicit(&done[dep].generation,
                                    memory_order_acquire) != generation) {
          if ((++spin & 255u) == 0u &&
              atomic_load_explicit(&shared->stop,
                                   memory_order_acquire) != 0) {
            return;
          }
          kls_cpu_relax();
        }
      }
      const double lik = x[dep] / udiag[dep];
      if (forward_work != NULL) {
        forward_value = fma(-lik, forward_work[dep], forward_value);
      }
      x[dep] = 0.0;
      if (row_values_mode) {
        l_val[p] = lik;
      } else {
        *l_lu[p] = lik;
      }
      if (lik != 0.0) {
        for (UF_long q = u_ptr[dep]; q < u_ptr[dep + 1u]; ++q) {
          const UF_long col = (UF_long)u_cols[q];
          x[col] = fma(-lik, u_val[q], x[col]);
        }
      }
    }
    const double pivot = x[row];
    x[row] = 0.0;
    udiag[row] = pivot;
    for (UF_long q = u_ptr[row]; q < u_ptr[row + 1u]; ++q) {
      const UF_long col = (UF_long)u_cols[q];
      const double value = x[col];
      u_val[q] = value;
      if (!row_values_mode) {
        *u_lu[q] = value;
      }
      x[col] = 0.0;
    }
    if (pivot == 0.0) {
      kls_egraph_refactor_record_singular(shared, row, row);
    }
    if (forward_work != NULL) forward_work[row] = forward_value;
    if (grouped_done_mode) {
      const uint32_t slot =
        done_token[row] & KLS_LEAN_GROUPED_SLOT_MASK;
      atomic_store_explicit(&grouped_done[slot], generation,
                            memory_order_release);
    } else {
      atomic_store_explicit(&done[row].generation, generation,
                            memory_order_release);
    }
    if (atomic_load_explicit(&shared->stop, memory_order_acquire) != 0) {
      return;
    }
  }
}

#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline, hot))
#endif
static void kls_generic_hoisted_worker_run(
  kls_egraph_refactor_worker *worker,
  unsigned int generation,
  const UF_long *restrict rows,
  UF_long stride) {
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  const UF_long n = solver->n;
  const unsigned int tid = (unsigned int)worker->tid;
  double *restrict x = worker->x;
  const double *restrict values = shared->values;
  const UF_long *restrict in_ptr = solver->row_refactor_input_ptr;
  const UF_long *restrict in_cols = solver->row_refactor_input_cols;
  const UF_long *restrict in_pos = solver->row_refactor_input_pos;
  const uint32_t *restrict in_packed =
    solver->row_refactor_input_col_user32;
  const UF_long *restrict l_ptr = solver->row_refactor_l_ptr;
  const UF_long *restrict l_cols = solver->row_refactor_l_cols;
  double **restrict l_lu = solver->row_refactor_l_values;
  double *restrict l_val = solver->row_refactor_l_row_values;
  const UF_long *restrict u_ptr = solver->row_refactor_u_ptr;
  const UF_long *restrict u_cols = solver->row_refactor_u_cols;
  double *restrict u_val = solver->row_refactor_u_row_values;
  double **restrict u_lu = solver->row_refactor_u_values;
  double *restrict udiag = (double *)solver->numeric->Udiag;
  kls_lean_done_slot *restrict done = shared->lean_done;
  const int row_values_mode = shared->lean_row_values_mode;
  const int grouped_done_mode = shared->lean_grouped_done_mode;
  const uint32_t *restrict done_token =
    solver->lean_parallel_grouped_token;
  atomic_uint *restrict grouped_done = solver->lean_parallel_grouped_done;
  double *restrict forward_work = shared->lean_forward_solve_work;

  for (UF_long pos = (UF_long)worker->tid; pos < n; pos += stride) {
    const UF_long row = rows[pos];
    if (row < solver->lean_scalar_btf_prefix) {
      continue;
    }
    double forward_value = 0.0;
    if (forward_work != NULL &&
        !kls_lean_forward_initial_value(
          shared, solver, row, &forward_value)) {
      kls_egraph_refactor_record_invalid(shared);
      return;
    }
    if (in_packed != NULL) {
      for (UF_long p = in_ptr[row]; p < in_ptr[row + 1u]; ++p) {
        const uint32_t packed = in_packed[p];
        x[packed & UINT32_C(0xffff)] = values[packed >> 16u];
      }
    } else {
      for (UF_long p = in_ptr[row]; p < in_ptr[row + 1u]; ++p) {
        x[in_cols[p]] = values[in_pos[p]];
      }
    }
    for (UF_long p = l_ptr[row]; p < l_ptr[row + 1u]; ++p) {
      const UF_long dep = l_cols[p];
      if (grouped_done_mode) {
        const uint32_t token = done_token[dep];
        if ((unsigned int)(token >> KLS_LEAN_GROUPED_OWNER_SHIFT) != tid) {
          const uint32_t slot = token & KLS_LEAN_GROUPED_SLOT_MASK;
          while (atomic_load_explicit(&grouped_done[slot],
                                      memory_order_acquire) != generation) {
            kls_cpu_relax();
          }
        }
      } else if (done[dep].owner != tid) {
        unsigned spin = 0u;
        while (atomic_load_explicit(&done[dep].generation,
                                    memory_order_acquire) != generation) {
          if ((++spin & 255u) == 0u &&
              atomic_load_explicit(&shared->stop,
                                   memory_order_acquire) != 0) {
            return;
          }
          kls_cpu_relax();
        }
      }
      const double lik = x[dep] / udiag[dep];
      if (forward_work != NULL) {
        forward_value = fma(-lik, forward_work[dep], forward_value);
      }
      x[dep] = 0.0;
      if (row_values_mode) {
        l_val[p] = lik;
      } else {
        *l_lu[p] = lik;
      }
      if (lik != 0.0) {
        for (UF_long q = u_ptr[dep]; q < u_ptr[dep + 1u]; ++q) {
          const UF_long col = u_cols[q];
          x[col] -= lik * u_val[q];
        }
      }
    }
    const double pivot = x[row];
    x[row] = 0.0;
    udiag[row] = pivot;
    for (UF_long q = u_ptr[row]; q < u_ptr[row + 1u]; ++q) {
      const UF_long col = u_cols[q];
      const double value = x[col];
      u_val[q] = value;
      if (!row_values_mode) {
        *u_lu[q] = value;
      }
      x[col] = 0.0;
    }
    if (pivot == 0.0) {
      kls_egraph_refactor_record_singular(shared, row, row);
    }
    if (forward_work != NULL) forward_work[row] = forward_value;
    if (grouped_done_mode) {
      const uint32_t slot =
        done_token[row] & KLS_LEAN_GROUPED_SLOT_MASK;
      atomic_store_explicit(&grouped_done[slot], generation,
                            memory_order_release);
    } else {
      atomic_store_explicit(&done[row].generation, generation,
                            memory_order_release);
    }
    if (atomic_load_explicit(&shared->stop, memory_order_acquire) != 0) {
      return;
    }
  }
}

/* Scaling does not change the retained row graph: every input entry in one
   pivot row uses the same original-row scale selected by Pnum.  Keep the
   unscaled hot loop above unchanged, and specialize the scaled loop here so
   the scale-mode test, row-value scatter choice, and representation checks
   remain outside the dependency walk. */
static KLS_ALWAYS_INLINE void kls_scaled_generic_hoisted_worker_core(
  kls_egraph_refactor_worker *worker,
  unsigned int generation,
  const UF_long *restrict rows,
  UF_long stride,
  int compact_cols,
  int compact_ptrs,
  int row_values_mode,
  int grouped_done_mode,
  int fused_scale_mode,
  int sum_scale_mode,
  int branchless_update_mode) {
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  const UF_long n = solver->n;
  const unsigned int tid = (unsigned int)worker->tid;
  double *restrict x = worker->x;
  const double *restrict values = shared->values;
  double *restrict rs = (double *)(uintptr_t)shared->rs;
  const UF_long *restrict pnum = solver->numeric->Pnum;
  const UF_long *restrict in_ptr = solver->row_refactor_input_ptr;
  const uint32_t *restrict in_ptr32 = solver->row_refactor_input_ptr32;
  const uint16_t *restrict in_ptr16 = solver->row_refactor_input_ptr16;
  const UF_long *restrict in_cols = solver->row_refactor_input_cols;
  const UF_long *restrict in_pos = solver->row_refactor_input_pos;
  const uint32_t *restrict in_packed =
    solver->row_refactor_input_col_user32;
  const UF_long *restrict l_ptr = solver->row_refactor_l_ptr;
  const uint32_t *restrict l_ptr32 = solver->row_refactor_l_ptr32;
  const uint16_t *restrict l_ptr16 = solver->row_refactor_l_ptr16;
  const UF_long *restrict l_cols = solver->row_refactor_l_cols;
  const uint16_t *restrict l_cols16 = solver->row_refactor_l_cols16;
  double **restrict l_lu = solver->row_refactor_l_values;
  double *restrict l_val = solver->row_refactor_l_row_values;
  const UF_long *restrict u_ptr = solver->row_refactor_u_ptr;
  const uint32_t *restrict u_ptr32 = solver->row_refactor_u_ptr32;
  const uint16_t *restrict u_ptr16 = solver->row_refactor_u_ptr16;
  const UF_long *restrict u_cols = solver->row_refactor_u_cols;
  const uint16_t *restrict u_cols16 = solver->row_refactor_u_cols16;
  double *restrict u_val = solver->row_refactor_u_row_values;
  double **restrict u_lu = solver->row_refactor_u_values;
  double *restrict udiag = (double *)solver->numeric->Udiag;
  kls_lean_done_slot *restrict done = shared->lean_done;
  const uint32_t *restrict done_token =
    solver->lean_parallel_grouped_token;
  atomic_uint *restrict grouped_done = solver->lean_parallel_grouped_done;

  for (UF_long pos = (UF_long)worker->tid; pos < n; pos += stride) {
    const UF_long row = rows[pos];
    const UF_long input_begin = compact_ptrs == 1
      ? (UF_long)in_ptr32[row]
      : compact_ptrs == 2 ? (UF_long)in_ptr16[row] : in_ptr[row];
    const UF_long input_end = compact_ptrs == 1
      ? (UF_long)in_ptr32[row + 1u]
      : compact_ptrs == 2 ? (UF_long)in_ptr16[row + 1u] : in_ptr[row + 1u];
    const UF_long dependency_begin = compact_ptrs == 1
      ? (UF_long)l_ptr32[row]
      : compact_ptrs == 2 ? (UF_long)l_ptr16[row] : l_ptr[row];
    const UF_long dependency_end = compact_ptrs == 1
      ? (UF_long)l_ptr32[row + 1u]
      : compact_ptrs == 2 ? (UF_long)l_ptr16[row + 1u] : l_ptr[row + 1u];
    double row_scale = 0.0;
    if (fused_scale_mode) {
      if (in_packed != NULL) {
        for (UF_long p = input_begin; p < input_end; ++p) {
          const uint32_t packed = in_packed[p];
          const double value = values[packed >> 16u];
          x[packed & UINT32_C(0xffff)] = value;
          const double magnitude = fabs(value);
          if (sum_scale_mode) {
            row_scale += magnitude;
          } else if (magnitude > row_scale) {
            row_scale = magnitude;
          }
        }
      } else {
        for (UF_long p = input_begin; p < input_end; ++p) {
          const double value = values[in_pos[p]];
          x[in_cols[p]] = value;
          const double magnitude = fabs(value);
          if (sum_scale_mode) {
            row_scale += magnitude;
          } else if (magnitude > row_scale) {
            row_scale = magnitude;
          }
        }
      }
      if (row_scale == 0.0) {
        row_scale = 1.0;
      }
      rs[row] = row_scale;
      const double row_rs_inv = 1.0 / row_scale;
      if (in_packed != NULL) {
        for (UF_long p = input_begin; p < input_end; ++p) {
          x[in_packed[p] & UINT32_C(0xffff)] *= row_rs_inv;
        }
      } else {
        for (UF_long p = input_begin; p < input_end; ++p) {
          x[in_cols[p]] *= row_rs_inv;
        }
      }
    } else {
      row_scale = rs[pnum[row]];
      const double row_rs_inv = 1.0 / row_scale;
      if (in_packed != NULL) {
        for (UF_long p = input_begin; p < input_end; ++p) {
          const uint32_t packed = in_packed[p];
          x[packed & UINT32_C(0xffff)] =
            values[packed >> 16u] * row_rs_inv;
        }
      } else {
        for (UF_long p = input_begin; p < input_end; ++p) {
          x[in_cols[p]] = values[in_pos[p]] * row_rs_inv;
        }
      }
    }
    for (UF_long p = dependency_begin; p < dependency_end; ++p) {
      const UF_long dep = compact_cols
        ? (UF_long)l_cols16[p] : l_cols[p];
      if (grouped_done_mode) {
        const uint32_t token = done_token[dep];
        if ((unsigned int)(token >> KLS_LEAN_GROUPED_OWNER_SHIFT) != tid) {
          const uint32_t slot = token & KLS_LEAN_GROUPED_SLOT_MASK;
          while (atomic_load_explicit(&grouped_done[slot],
                                      memory_order_acquire) != generation) {
            kls_cpu_relax();
          }
        }
      } else if (done[dep].owner != tid) {
        unsigned spin = 0u;
        while (atomic_load_explicit(&done[dep].generation,
                                    memory_order_acquire) != generation) {
          if ((++spin & 255u) == 0u &&
              atomic_load_explicit(&shared->stop,
                                   memory_order_acquire) != 0) {
            return;
          }
          kls_cpu_relax();
        }
      }
      const double lik = x[dep] / udiag[dep];
      x[dep] = 0.0;
      if (row_values_mode) {
        l_val[p] = lik;
      } else {
        *l_lu[p] = lik;
      }
      if (branchless_update_mode || lik != 0.0) {
        const UF_long update_begin = compact_ptrs == 1
          ? (UF_long)u_ptr32[dep]
          : compact_ptrs == 2 ? (UF_long)u_ptr16[dep] : u_ptr[dep];
        const UF_long update_end = compact_ptrs == 1
          ? (UF_long)u_ptr32[dep + 1u]
          : compact_ptrs == 2 ? (UF_long)u_ptr16[dep + 1u] : u_ptr[dep + 1u];
        for (UF_long q = update_begin; q < update_end; ++q) {
          const UF_long col = compact_cols
            ? (UF_long)u_cols16[q] : u_cols[q];
          x[col] = fma(-lik, u_val[q], x[col]);
        }
      }
    }
    const double pivot = x[row];
    x[row] = 0.0;
    udiag[row] = pivot;
    const UF_long output_begin = compact_ptrs == 1
      ? (UF_long)u_ptr32[row]
      : compact_ptrs == 2 ? (UF_long)u_ptr16[row] : u_ptr[row];
    const UF_long output_end = compact_ptrs == 1
      ? (UF_long)u_ptr32[row + 1u]
      : compact_ptrs == 2 ? (UF_long)u_ptr16[row + 1u] : u_ptr[row + 1u];
    for (UF_long q = output_begin; q < output_end; ++q) {
      const UF_long col = compact_cols
        ? (UF_long)u_cols16[q] : u_cols[q];
      const double value = x[col];
      u_val[q] = value;
      if (!row_values_mode) {
        *u_lu[q] = value;
      }
      x[col] = 0.0;
    }
    if (pivot == 0.0) {
      kls_egraph_refactor_record_singular(shared, row, row);
    }
    if (grouped_done_mode) {
      const uint32_t slot =
        done_token[row] & KLS_LEAN_GROUPED_SLOT_MASK;
      atomic_store_explicit(&grouped_done[slot], generation,
                            memory_order_release);
    } else {
      atomic_store_explicit(&done[row].generation, generation,
                            memory_order_release);
    }
    if (atomic_load_explicit(&shared->stop, memory_order_acquire) != 0) {
      return;
    }
  }
}

#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline, hot))
#endif
static void kls_scaled_generic_hoisted_worker_run(
  kls_egraph_refactor_worker *worker,
  unsigned int generation,
  const UF_long *restrict rows,
  UF_long stride) {
  kls_scaled_generic_hoisted_worker_core(
    worker, generation, rows, stride, 0, 0,
    worker->shared->lean_row_values_mode,
    worker->shared->lean_grouped_done_mode,
    worker->shared->lean_fused_scale_mode,
    worker->shared->solver->common.scale == 1, 0);
}

#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline, hot))
#endif
static void kls_scaled_generic_i16_hoisted_worker_run(
  kls_egraph_refactor_worker *worker,
  unsigned int generation,
  const UF_long *restrict rows,
  UF_long stride) {
  kls_scaled_generic_hoisted_worker_core(
    worker, generation, rows, stride, 1, 0,
    worker->shared->lean_row_values_mode,
    worker->shared->lean_grouped_done_mode,
    worker->shared->lean_fused_scale_mode,
    worker->shared->solver->common.scale == 1, 0);
}

#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline, hot))
#endif
static void kls_scaled_generic_i16ptr_hoisted_worker_run(
  kls_egraph_refactor_worker *worker,
  unsigned int generation,
  const UF_long *restrict rows,
  UF_long stride) {
  kls_scaled_generic_hoisted_worker_core(
    worker, generation, rows, stride, 1, 2,
    worker->shared->lean_row_values_mode,
    worker->shared->lean_grouped_done_mode,
    worker->shared->lean_fused_scale_mode,
    worker->shared->solver->common.scale == 1, 0);
}

#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline, hot))
#endif
static void kls_scaled_generic_i16_u32ptr_hoisted_worker_run(
  kls_egraph_refactor_worker *worker,
  unsigned int generation,
  const UF_long *restrict rows,
  UF_long stride) {
  kls_scaled_generic_hoisted_worker_core(
    worker, generation, rows, stride, 1, 1,
    worker->shared->lean_row_values_mode,
    worker->shared->lean_grouped_done_mode,
    worker->shared->lean_fused_scale_mode,
    worker->shared->solver->common.scale == 1, 0);
}

#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline, hot))
#endif
static void kls_scaled_generic_i16_u32ptr_row_values_worker_run(
  kls_egraph_refactor_worker *worker,
  unsigned int generation,
  const UF_long *restrict rows,
  UF_long stride) {
  kls_scaled_generic_hoisted_worker_core(
    worker, generation, rows, stride, 1, 1, 1,
    worker->shared->lean_grouped_done_mode,
    worker->shared->lean_fused_scale_mode,
    worker->shared->solver->common.scale == 1, 0);
}

#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline, hot))
#endif
static void kls_scaled_generic_i16_u32ptr_row_values_grouped_branchless_worker_run(
  kls_egraph_refactor_worker *worker,
  unsigned int generation,
  const UF_long *restrict rows,
  UF_long stride) {
  kls_scaled_generic_hoisted_worker_core(
    worker, generation, rows, stride, 1, 1, 1, 1,
    worker->shared->lean_fused_scale_mode,
    worker->shared->solver->common.scale == 1, 1);
}

#if defined(__GNUC__) || defined(__clang__)
#define KLS_SCALED_EXACT_WORKER_ATTR __attribute__((noinline, hot))
#else
#define KLS_SCALED_EXACT_WORKER_ATTR
#endif

#define KLS_DEFINE_SCALED_EXACT_WORKER(name, sum_mode) \
  KLS_SCALED_EXACT_WORKER_ATTR static void name(                       \
    kls_egraph_refactor_worker *worker,                                \
    unsigned int generation,                                          \
    const UF_long *restrict rows,                                      \
    UF_long stride) {                                                  \
    kls_scaled_generic_hoisted_worker_core(                            \
      worker, generation, rows, stride, 1, 1, 1, 1,               \
      2, sum_mode, 1);                                                 \
  }

KLS_DEFINE_SCALED_EXACT_WORKER(
  kls_scaled_exact_direct_max_branchless_worker_run, 0)
KLS_DEFINE_SCALED_EXACT_WORKER(
  kls_scaled_exact_direct_sum_branchless_worker_run, 1)

#undef KLS_DEFINE_SCALED_EXACT_WORKER
#undef KLS_SCALED_EXACT_WORKER_ATTR

static void kls_lean_parallel_worker_run(
  kls_egraph_refactor_worker *worker) {
  if (worker == NULL || worker->shared == NULL || worker->x == NULL) {
    return;
  }
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  if (solver == NULL || shared->lean_done == NULL ||
      shared->pipeline_generation == 0u ||
      solver->row_refactor_level_ptr == NULL ||
      solver->row_refactor_level_rows == NULL ||
      shared->thread_count < 2 || worker->tid < 0 ||
      worker->tid >= shared->thread_count) {
    kls_egraph_refactor_record_invalid(shared);
    return;
  }
  const int symmetric_scalar_fringe_lean_class =
    shared->lean_symmetric_scalar_fringe_mode;
  if (!solver->lean_parallel_scratch_clean) {
    memset(worker->x, 0, (size_t)solver->n * sizeof(*worker->x));
  }
  if (shared->value_prep_input != NULL &&
      shared->value_prep_output != NULL &&
      shared->value_prep_input_to_csc != NULL) {
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
    pthread_barrier_wait(&shared->barrier);
    const UF_long prefix = solver->lean_scalar_btf_prefix;
    const UF_long *restrict in_ptr = solver->row_refactor_input_ptr;
    const UF_long *restrict in_cols = solver->row_refactor_input_cols;
    const UF_long *restrict in_pos = solver->row_refactor_input_pos;
    double *restrict udiag = (double *)solver->numeric->Udiag;
    for (UF_long row = (UF_long)worker->tid; row < prefix;
         row += (UF_long)shared->thread_count) {
      int found = 0;
      double pivot = 0.0;
      for (UF_long p = in_ptr[row]; p < in_ptr[row + 1u]; ++p) {
        if (in_cols[p] == row) {
          pivot = output[in_pos[p]];
          found = 1;
          break;
        }
      }
      if (!found) {
        kls_egraph_refactor_record_invalid(shared);
        break;
      }
      udiag[row] = pivot;
      if (pivot == 0.0) {
        kls_egraph_refactor_record_singular(shared, row, row);
      }
    }
    pthread_barrier_wait(&shared->barrier);
    if (shared->invalid ||
        (shared->singular && solver->common.halt_if_singular)) {
      return;
    }
  }
  const unsigned int generation = shared->pipeline_generation;
  const UF_long *restrict rows = shared->lean_rows != NULL
    ? shared->lean_rows : solver->row_refactor_level_rows;
  const UF_long stride = (UF_long)shared->thread_count;
  const int symmetric_scalar_fringe_hoisted_worker =
    symmetric_scalar_fringe_lean_class && shared->rs == NULL &&
    solver->row_refactor_input_ptr != NULL &&
    solver->row_refactor_input_cols != NULL &&
    solver->row_refactor_input_pos != NULL &&
    solver->row_refactor_l_ptr != NULL &&
    solver->row_refactor_l_cols != NULL &&
    solver->row_refactor_l_values != NULL &&
    solver->row_refactor_u_ptr != NULL &&
    solver->row_refactor_u_cols != NULL &&
    solver->row_refactor_u_values != NULL &&
    solver->row_refactor_u_row_values != NULL &&
    solver->numeric != NULL && solver->numeric->Udiag != NULL;
  if (symmetric_scalar_fringe_hoisted_worker) {
    const int compact_snode_cols =
      solver->row_refactor_l_cols16 != NULL &&
      solver->row_refactor_u_cols16 != NULL;
#if KLS_HAVE_AVX512_KERNELS
    if (compact_snode_cols && shared->lean_snode_avx512_enabled &&
        __builtin_cpu_supports("avx2")) {
      if (solver->lean_snode_avx512_candidate > 1) {
        kls_symmetric_scalar_fringe_hoisted_worker_run_i16_width16_avx512(
          worker, generation, rows, stride);
      } else {
        kls_symmetric_scalar_fringe_hoisted_worker_run_i16_avx512(
          worker, generation, rows, stride);
      }
    } else if (compact_snode_cols) {
      kls_symmetric_scalar_fringe_hoisted_worker_run_i16(
        worker, generation, rows, stride);
    } else if (shared->lean_snode_avx512_enabled) {
      if (solver->lean_snode_avx512_candidate > 1) {
        kls_symmetric_scalar_fringe_hoisted_worker_run_width16_avx512(
          worker, generation, rows, stride);
      } else {
        kls_symmetric_scalar_fringe_hoisted_worker_run_avx512(
          worker, generation, rows, stride);
      }
    } else {
      kls_symmetric_scalar_fringe_hoisted_worker_run(
        worker, generation, rows, stride);
    }
#else
    if (compact_snode_cols) {
      kls_symmetric_scalar_fringe_hoisted_worker_run_i16(
        worker, generation, rows, stride);
    } else {
      kls_symmetric_scalar_fringe_hoisted_worker_run(
        worker, generation, rows, stride);
    }
#endif
    if (!kls_lean_parallel_refresh_offdiag_worker(worker, stride)) {
      kls_egraph_refactor_record_invalid(shared);
    }
    return;
  }
  const int generic_hoisted_worker =
    !shared->lean_compact_match_mode &&
    !symmetric_scalar_fringe_lean_class &&
    solver->lean_snode_run == NULL &&
    solver->row_refactor_input_ptr != NULL &&
    solver->row_refactor_input_cols != NULL &&
    solver->row_refactor_input_pos != NULL &&
    solver->row_refactor_l_ptr != NULL &&
    solver->row_refactor_l_cols != NULL &&
    solver->row_refactor_l_values != NULL &&
    solver->row_refactor_u_ptr != NULL &&
    solver->row_refactor_u_cols != NULL &&
    solver->row_refactor_u_values != NULL &&
    solver->row_refactor_u_row_values != NULL &&
    solver->numeric != NULL && solver->numeric->Udiag != NULL &&
    (shared->rs == NULL || solver->numeric->Pnum != NULL);
  if (generic_hoisted_worker) {
    if (shared->rs != NULL) {
      if (shared->lean_generic_packed_mode &&
          solver->lean_scalar_btf_prefix == 0u &&
          solver->lean_generic_packed_rows32 != NULL &&
          solver->lean_generic_packed_dependencies32 != NULL &&
          solver->lean_generic_packed_publish_mailbox != NULL &&
          solver->lean_generic_packed_publish_begin != NULL &&
          solver->lean_parallel_owner_frontier != NULL &&
          solver->row_refactor_u_cols16 != NULL &&
          (solver->row_refactor_input_col_user32 != NULL ||
           (solver->row_refactor_input_cols16 != NULL &&
            solver->row_refactor_input_pos32 != NULL))) {
        kls_scaled_generic_packed_row_worker_run(worker, stride);
      } else if (solver->row_refactor_l_cols16 != NULL &&
          solver->row_refactor_u_cols16 != NULL &&
          solver->row_refactor_l_ptr16 != NULL &&
          solver->row_refactor_u_ptr16 != NULL &&
          solver->row_refactor_input_ptr16 != NULL) {
        kls_scaled_generic_i16ptr_hoisted_worker_run(
          worker, generation, rows, stride);
      } else if (solver->row_refactor_l_cols16 != NULL &&
          solver->row_refactor_u_cols16 != NULL &&
          solver->row_refactor_l_ptr32 != NULL &&
          solver->row_refactor_u_ptr32 != NULL &&
          solver->row_refactor_input_ptr32 != NULL) {
        if (shared->lean_row_values_mode &&
            shared->lean_grouped_done_mode) {
          if (shared->lean_fused_scale_mode > 1) {
            if (solver->common.scale == 1) {
              kls_scaled_exact_direct_sum_branchless_worker_run(
                worker, generation, rows, stride);
            } else {
              kls_scaled_exact_direct_max_branchless_worker_run(
                worker, generation, rows, stride);
            }
          } else {
            kls_scaled_generic_i16_u32ptr_row_values_grouped_branchless_worker_run(
              worker, generation, rows, stride);
          }
        } else if (shared->lean_row_values_mode) {
          kls_scaled_generic_i16_u32ptr_row_values_worker_run(
            worker, generation, rows, stride);
        } else {
          kls_scaled_generic_i16_u32ptr_hoisted_worker_run(
            worker, generation, rows, stride);
        }
      } else if (solver->row_refactor_l_cols16 != NULL &&
          solver->row_refactor_u_cols16 != NULL) {
        kls_scaled_generic_i16_hoisted_worker_run(
          worker, generation, rows, stride);
      } else {
        kls_scaled_generic_hoisted_worker_run(
          worker, generation, rows, stride);
      }
      if (!kls_lean_parallel_refresh_offdiag_worker(worker, stride)) {
        kls_egraph_refactor_record_invalid(shared);
      }
      return;
    }
    const int compact_amf_two_block_specialized_worker =
      kls_packed_row_worker_representation_capable(solver) &&
      shared->lean_grouped_done_mode && shared->lean_row_values_mode &&
      solver->lean_scalar_btf_prefix == 0u &&
      solver->row_refactor_input_ptr16 != NULL &&
      solver->row_refactor_input_col_user32 != NULL &&
      solver->row_refactor_l_ptr16 != NULL &&
      solver->row_refactor_l_cols16 != NULL &&
      solver->row_refactor_l_row_values != NULL &&
      solver->row_refactor_u_ptr16 != NULL &&
      solver->row_refactor_u_cols16 != NULL &&
      solver->lean_parallel_grouped_token != NULL &&
      solver->lean_parallel_grouped_done != NULL &&
      solver->lean_parallel_owner_frontier != NULL &&
      solver->lean_parallel_owner_frontier_count == shared->thread_count &&
      solver->lean_parallel_owner_frontier_sequence_bits > 0u &&
      solver->lean_parallel_owner_frontier_sequence_mask > 0u;
    if (compact_amf_two_block_specialized_worker) {
      if (solver->lean_parallel_l_dep_work64 != NULL &&
          solver->lean_parallel_work_rows16 != NULL &&
          solver->lean_parallel_publish_mailbox != NULL &&
          solver->lean_parallel_publish_begin != NULL &&
          solver->lean_parallel_udiag_inv != NULL) {
        kls_compact_amf_two_block_packed_exact_worker_run(
          worker, stride);
      } else {
        kls_compact_amf_two_block_i16ptr_worker_run(
          worker, generation, rows, stride);
      }
    } else if (shared->lean_generic_packed_mode &&
        solver->lean_scalar_btf_prefix == 0u &&
        solver->lean_generic_packed_rows32 != NULL &&
        solver->lean_generic_packed_dependencies32 != NULL &&
        solver->lean_generic_packed_publish_mailbox != NULL &&
        solver->lean_generic_packed_publish_begin != NULL &&
        solver->lean_parallel_owner_frontier != NULL) {
      kls_generic_packed_row_worker_run(worker, stride);
    } else if (solver->row_refactor_l_cols16 != NULL &&
        solver->row_refactor_u_cols16 != NULL) {
      if (solver->row_refactor_l_ptr16 != NULL &&
          solver->row_refactor_u_ptr16 != NULL &&
          solver->row_refactor_input_ptr16 != NULL) {
        const uint16_t *rows16 =
          rows == solver->row_refactor_level_rows
            ? solver->row_refactor_level_rows16 : NULL;
        if (shared->udiag_recip != NULL &&
            shared->lean_grouped_done_mode &&
            shared->lean_row_values_mode) {
          kls_generic_i16ptr_grouped_recip_worker_run(
            worker, generation, rows, rows16, stride);
        } else {
          kls_generic_i16ptr_hoisted_worker_run(
            worker, generation, rows, rows16, stride);
        }
      } else {
        kls_generic_i16_hoisted_worker_run(
          worker, generation, rows, stride);
      }
    } else {
      kls_generic_hoisted_worker_run(worker, generation, rows, stride);
    }
    if (!kls_lean_parallel_refresh_offdiag_worker(worker, stride)) {
      kls_egraph_refactor_record_invalid(shared);
    }
    return;
  }
  if (shared->lean_compact_match_mode) {
    const int hoisted_worker =
      solver->lean_user_values_active &&
      shared->lean_grouped_done_mode && shared->lean_row_values_mode &&
      solver->row_refactor_input_user_pos != NULL &&
      solver->lean_parallel_grouped_token != NULL &&
      solver->lean_parallel_grouped_done != NULL;
    if (hoisted_worker) {
      kls_compact_match_hoisted_worker_run(
        worker, generation, rows, stride);
      return;
    } else {
      for (UF_long pos = (UF_long)worker->tid; pos < solver->n;
           pos += stride) {
        const UF_long row = rows[pos];
        if (!kls_compact_match_lean_parallel_process_row(
              worker, row, generation)) {
          kls_egraph_refactor_record_invalid(shared);
          return;
        }
        if (shared->lean_grouped_done_mode) {
          const uint32_t slot =
            solver->lean_parallel_grouped_token[row] &
            KLS_LEAN_GROUPED_SLOT_MASK;
          atomic_store_explicit(&solver->lean_parallel_grouped_done[slot],
                                generation, memory_order_release);
        } else {
          atomic_store_explicit(&shared->lean_done[row].generation,
                                generation, memory_order_release);
        }
      }
    }
    if (solver->lean_user_values_active) {
      double *restrict offx = (double *)solver->numeric->Offx;
      const uint16_t *restrict off_user =
        solver->compact_match_offdiag_user_pos;
      const UF_long offcount = solver->numeric->Offp[solver->n];
      for (UF_long p = (UF_long)worker->tid; p < offcount; p += stride) {
        offx[p] = shared->values[off_user[p]];
      }
    } else if (!kls_lean_parallel_refresh_offdiag_worker(worker, stride)) {
      /* The internal-frame compact worker refreshes the diagonal blocks
         above, but BTF Offx entries live outside those row mirrors.  Leaving
         them at the preceding generation would make the raw solve an
         implicit old/new hybrid. */
      kls_egraph_refactor_record_invalid(shared);
    }
    return;
  }
  for (UF_long pos = (UF_long)worker->tid; pos < solver->n;
       pos += stride) {
    const UF_long row = rows[pos];
    if (row >= solver->n ||
        !kls_lean_parallel_process_row(worker, row, generation)) {
      kls_egraph_refactor_record_invalid(shared);
      return;
    }
    if (shared->lean_grouped_done_mode) {
      const uint32_t slot =
        solver->lean_parallel_grouped_token[row] &
        KLS_LEAN_GROUPED_SLOT_MASK;
      atomic_store_explicit(&solver->lean_parallel_grouped_done[slot],
                            generation, memory_order_release);
    } else {
      atomic_store_explicit(&shared->lean_done[row].generation, generation,
                            memory_order_release);
    }
    if (atomic_load_explicit(&shared->stop, memory_order_acquire) != 0) {
      return;
    }
  }
  if (!kls_lean_parallel_refresh_offdiag_worker(worker, stride)) {
    kls_egraph_refactor_record_invalid(shared);
  }
}

static int kls_lean_parallel_refactor_run(kls_solver *solver,
                                          double *numeric_values,
                                          const double *rs,
                                          int fused_scale_mode) {
  if (solver == NULL || numeric_values == NULL ||
      solver->row_refactor_level_rows == NULL ||
      solver->row_refactor_level_max_width < 2u ||
      solver->options.threads < 2) {
    return -1;
  }
  int thread_count = solver->options.threads;
  if (kls_packed_row_worker_representation_capable(solver)) {
    thread_count = kls_packed_row_worker_thread_count(solver, thread_count);
  } else if ((UF_long)thread_count > solver->row_refactor_level_max_width) {
    thread_count = (int)solver->row_refactor_level_max_width;
  }
  if (thread_count < 2) {
    return -1;
  }
  double **scratch =
    ensure_egraph_worker_scratch(solver, thread_count, solver->n);
  unsigned int generation = 0u;
  kls_lean_done_slot *done =
    ensure_lean_parallel_done(solver, &generation);
  kls_egraph_refactor_pool *pool =
    ensure_egraph_refactor_pool(solver, thread_count);
  if (scratch == NULL || done == NULL || generation == 0u || pool == NULL) {
    return -1;
  }
  if (kls_packed_row_worker_representation_capable(solver)) {
    if (pool->worker_spin_iters < 200000u) {
      /* Keep the already-created packed crew live across the bounded solve
         gap.  Its numeric phase is short enough that a futex wake is a
         material fraction of each update, independently of whether a
         separate Offx stream happens to be present. */
      pool->worker_spin_iters = 200000u;
    }
  }

  int grouped_done_ready =
    solver->lean_parallel_grouped_done != NULL &&
    solver->lean_parallel_grouped_token != NULL;
  const UF_long *rows =
    kls_prepare_lean_affinity_rows(solver, thread_count);
  if (rows == NULL) {
    return -1;
  }
  if (solver->lean_parallel_owner_thread_count != thread_count) {
    for (UF_long pos = 0u; pos < solver->n; ++pos) {
      const UF_long row = rows[pos];
      if (row >= solver->n) {
        return -1;
      }
      done[row].owner = (unsigned int)(pos % (UF_long)thread_count);
    }
    grouped_done_ready =
      kls_prepare_lean_grouped_done(solver, thread_count, rows);
    solver->lean_parallel_owner_thread_count = thread_count;
  }
  if (solver->lean_snode_run != NULL) {
    (void)kls_prepare_lean_snode_wait_masks(solver, thread_count, rows);
  }
  const int deferred_value_prep =
    solver->lean_deferred_value_prep_active && rs == NULL &&
    kls_deferred_lean_value_prep_capable(solver);
  kls_egraph_refactor_shared *shared = &pool->shared;
  pthread_mutex_lock(&shared->lock);
  if (pool->active_workers != 0) {
    pthread_mutex_unlock(&shared->lock);
    return -1;
  }
  shared->solver = solver;
  shared->values = deferred_value_prep ? solver->values : numeric_values;
  shared->rs = rs;
  shared->thread_count = thread_count;
  shared->lean_done = done;
  shared->lean_rows = rows;
  shared->pipeline_generation = generation;
  shared->lean_fused_scale_mode = fused_scale_mode;
  shared->row_refactor_defer_value_scatter = 0;
  shared->lean_compact_match_mode =
    (solver->lean_user_values_active &&
     solver->row_refactor_input_user_pos != NULL &&
     solver->compact_match_offdiag_user_pos != NULL);
  shared->lean_parallel_offdiag_mode =
    deferred_value_prep ||
    (!solver->lean_user_values_active &&
     kls_parallel_offdiag_is_profitable(solver));
  shared->value_prep_input =
    deferred_value_prep ? numeric_values : NULL;
  shared->value_prep_output =
    deferred_value_prep ? solver->values : NULL;
  shared->value_prep_input_to_csc =
    deferred_value_prep ? solver->prepared_value_input_pos : NULL;
  shared->value_prep_nnz = deferred_value_prep ? solver->nnz : 0u;
  const int symmetric_scalar_fringe_lean_mode =
    !shared->lean_compact_match_mode && grouped_done_ready &&
    solver->lean_snode_wait_mask != NULL &&
    solver->lean_snode_wait_slot != NULL &&
    (kls_generic_hoisted_snode_worker_capable(solver));
  shared->lean_symmetric_scalar_fringe_mode =
    symmetric_scalar_fringe_lean_mode;
  shared->lean_snode_avx512_enabled = 0;
#if KLS_HAVE_AVX512_KERNELS
  if (symmetric_scalar_fringe_lean_mode &&
      solver->lean_snode_avx512_candidate) {
    /* Resolve the ISA once on the dispatching thread.  Workers then read an
       immutable mode bit instead of consulting the runtime feature cache in
       every fused run. */
    shared->lean_snode_avx512_enabled = kls_avx512_scatter_enabled();
  }
#endif
  /* Grouped completion tokens are valid only for workers that consume that
     scoreboard.  The fused supernode worker consumes them independently of
     this generic-mode flag.  Compact-match workers have their own grouped
     consumer, while the generic hoisted worker requires the same
     no-supernode representation as its dispatch gate below. */
  const int generic_grouped_done_consumer =
    shared->lean_compact_match_mode ||
    (!symmetric_scalar_fringe_lean_mode &&
     solver->lean_snode_run == NULL);
  const int generic_grouped_done =
    grouped_done_ready &&
    generic_grouped_done_consumer &&
    kls_lean_grouped_done_is_profitable(solver, thread_count, rows);
  shared->lean_grouped_done_mode =
    grouped_done_ready &&
    (generic_grouped_done &&
          /* Deferred gathering changes which worker publishes the first
             ready generation.  The compact grouped-token proof assumes the
             already-prepared row stream and cannot represent that extra
             producer phase; combining them can publish a numerically
             incomplete dependency generation despite a successful return.
             Keep both optimizations independently, but never compose them. */
          !deferred_value_prep);
  const UF_long factor_entries =
    solver->numeric->lnz + solver->numeric->unz;
  const int generic_scaled_row_values =
    rs != NULL && factor_entries > 0u &&
    solver->common.flops >= 4.0 * (double)factor_entries;
  shared->lean_row_values_mode =
    shared->lean_compact_match_mode ||
     ((rs == NULL || generic_scaled_row_values) &&
      solver->i16solve_l != NULL && solver->i16solve_u != NULL &&
      solver->i16solve_loff != NULL && solver->i16solve_uoff != NULL);
  shared->lean_forward_solve_work =
    solver->fused_refactor_solve_requested &&
      solver->symbolic != NULL && solver->symbolic->nblocks == 1u &&
      solver->symbolic->R != NULL &&
      solver->common.scale <= 0 && solver->numeric->Rs == NULL &&
      solver->row_refactor_level_rows != NULL &&
      solver->row_refactor_l_ptr != NULL &&
      solver->row_refactor_l_cols != NULL &&
      solver->row_refactor_u_ptr != NULL &&
      solver->row_refactor_u_cols != NULL
    ? solver->fused_refactor_solve_work : NULL;
  shared->lean_forward_solve_rhs =
    shared->lean_forward_solve_work != NULL
      ? solver->fused_refactor_solve_rhs : NULL;
  shared->lean_forward_solve_begin =
    shared->lean_forward_solve_work != NULL
      ? 0u
      : solver->n;
  if (symmetric_scalar_fringe_lean_mode &&
      shared->lean_row_values_mode &&
      pool->worker_spin_iters < 400000u) {
    /* This executor removes the packed-value scatter and follows the
       numeric with a compact row-mirror solve.  Keep the persistent crew
       live through that bounded solve/verification gap; otherwise these
       shorter cycles intermittently cross the futex boundary and turn the
       next microsecond-scale numeric into a wakeup benchmark. */
    pool->worker_spin_iters = 400000u;
  }
  const int generic_i16ptr_grouped_recip_capable =
    rs == NULL &&
    !shared->lean_compact_match_mode &&
    !symmetric_scalar_fringe_lean_mode &&
    /* Multiplication by a rounded cached reciprocal is not equivalent to
       division when fixed-pivot refactor error is strongly amplified.  On
       such numerics the extra rounding can turn an otherwise backward-
       stable row update into a poor factor (the solve residual, not merely
       the factor entries, exposes the loss).  Keep the same grouped row
       representation and dependency schedule, but let the ordinary worker
       divide by Udiag directly whenever the retained factor's condition
       estimate is below the conventional 1e-6 accuracy-risk threshold.
       This is a numeric certificate independent of input dimensions,
       sparsity shape, ordering, or matrix identity. */
    solver->common.rcond >= 1.0e-6 &&
    solver->lean_snode_run == NULL &&
    shared->lean_grouped_done_mode && shared->lean_row_values_mode &&
    solver->row_refactor_input_ptr16 != NULL &&
    solver->row_refactor_l_ptr16 != NULL &&
    solver->row_refactor_l_cols16 != NULL &&
    solver->row_refactor_l_row_values != NULL &&
    solver->row_refactor_u_ptr16 != NULL &&
    solver->row_refactor_u_cols16 != NULL &&
    solver->row_refactor_u_row_values != NULL &&
    solver->lean_parallel_grouped_token != NULL &&
    solver->lean_parallel_grouped_done != NULL &&
    (solver->i32solve_rhs_perm32 != NULL ||
     solver->i16solve_rhs_perm != NULL);
  int generic_i16ptr_grouped_recip_worker = 0;

  shared->lean_refactor_mode = 1;
  shared->row_refactor_mode = 0;
  shared->row_solve_mode = 0;
  shared->row_publish_mode = 0;
  shared->invalid = 0;
  shared->pivot_rejected = 0;
  shared->singular = 0;
  shared->numerical_rank = UF_long_max;
  shared->singular_col = KLS_KLU_EMPTY;
  atomic_store_explicit(&shared->stop, 0, memory_order_release);
  solver->lean_scalar_btf_prefix = 0u;
  if (rs == NULL &&
      !shared->lean_compact_match_mode &&
      !symmetric_scalar_fringe_lean_mode &&
      solver->lean_snode_run == NULL && solver->symbolic != NULL &&
      solver->symbolic->R != NULL) {
    UF_long block = 0u;
    while (block < solver->symbolic->nblocks &&
           solver->symbolic->R[block + 1u] ==
             solver->symbolic->R[block] + 1u) {
      block++;
    }
    const UF_long prefix = solver->symbolic->R[block];
    if (prefix >= 64u) {
      if (deferred_value_prep) {
        solver->lean_scalar_btf_prefix = prefix;
      } else {
        const UF_long *restrict in_ptr = solver->row_refactor_input_ptr;
        const UF_long *restrict in_cols = solver->row_refactor_input_cols;
        const UF_long *restrict in_pos = solver->row_refactor_input_pos;
        double *restrict udiag = (double *)solver->numeric->Udiag;
        for (UF_long row = 0u; row < prefix; ++row) {
          int found = 0;
          double pivot = 0.0;
          for (UF_long p = in_ptr[row]; p < in_ptr[row + 1u]; ++p) {
            if (in_cols[p] == row) {
              pivot = numeric_values[in_pos[p]];
              found = 1;
              break;
            }
          }
          if (!found) {
            kls_egraph_refactor_record_invalid(shared);
            break;
          }
          udiag[row] = pivot;
          if (generic_i16ptr_grouped_recip_capable &&
              solver->i32solve_udiag_recip != NULL) {
            solver->i32solve_udiag_recip[row] = 1.0 / pivot;
          }
          if (pivot == 0.0) {
            kls_egraph_refactor_record_singular(shared, row, row);
          }
        }
        if (!shared->invalid) {
          solver->lean_scalar_btf_prefix = prefix;
        }
      }
    }
  }
  if (generic_i16ptr_grouped_recip_capable &&
      solver->lean_scalar_btf_prefix > 0u) {
    if (solver->i32solve_udiag_recip == NULL) {
      double *recip = (double *)malloc(
        (size_t)(solver->n > 0u ? solver->n : 1u) * sizeof(*recip));
      if (recip != NULL) {
        const double *restrict udiag =
          (const double *)solver->numeric->Udiag;
#pragma omp simd
        for (UF_long row = 0u; row < solver->n; ++row) {
          recip[row] = 1.0 / udiag[row];
        }
        solver->i32solve_udiag_recip = recip;
      }
    }
    generic_i16ptr_grouped_recip_worker =
      solver->i32solve_udiag_recip != NULL;
  }
  shared->udiag_recip = generic_i16ptr_grouped_recip_worker
    ? solver->i32solve_udiag_recip : NULL;
  if (generic_i16ptr_grouped_recip_worker &&
      pool->worker_spin_iters < 200000u) {
    /* This mixed scalar-prefix/grouped numeric is followed by a compact
       fused solve measured in only a few dozen microseconds.  Keep its
       already-created crew live across that bounded gap; a futex wake is
       otherwise larger than the complete recurring kernel and produces a
       bimodal update rate.  Eligibility is the exact executor/solve
       contract above, not an input-shape policy. */
    pool->worker_spin_iters = 200000u;
  }
  const int compact_amf_two_block_exact_worker =
    kls_packed_row_worker_representation_capable(solver) &&
    shared->lean_grouped_done_mode && shared->lean_row_values_mode &&
    solver->lean_scalar_btf_prefix == 0u &&
    solver->row_refactor_input_ptr16 != NULL &&
    solver->row_refactor_input_col_user32 != NULL &&
    solver->row_refactor_l_ptr16 != NULL &&
    solver->row_refactor_l_cols16 != NULL &&
    solver->row_refactor_l_row_values != NULL &&
    solver->row_refactor_u_ptr16 != NULL &&
    solver->row_refactor_u_cols16 != NULL &&
    solver->lean_parallel_grouped_token != NULL &&
    solver->lean_parallel_grouped_done != NULL &&
    solver->lean_parallel_owner_frontier != NULL &&
    solver->lean_parallel_owner_frontier_count == thread_count &&
    solver->lean_parallel_owner_frontier_sequence_bits > 0u &&
    solver->lean_parallel_owner_frontier_sequence_mask > 0u &&
    solver->lean_parallel_l_dep_work64 != NULL &&
    solver->lean_parallel_work_rows16 != NULL &&
    solver->lean_parallel_publish_mailbox != NULL &&
    solver->lean_parallel_publish_begin != NULL &&
    solver->lean_parallel_udiag_inv != NULL &&
    solver->i32solve_udiag_recip != NULL;
  const int generic_packed_row_worker =
    !compact_amf_two_block_exact_worker &&
    solver->lean_scalar_btf_prefix == 0u &&
    solver->lean_generic_packed_rows32 != NULL &&
    solver->lean_generic_packed_dependencies32 != NULL &&
    solver->lean_generic_packed_publish_mailbox != NULL &&
    solver->lean_generic_packed_publish_begin != NULL &&
    solver->lean_parallel_owner_frontier != NULL;
  const int scaled_generic_packed_row_worker =
    shared->rs != NULL &&
    solver->lean_scalar_btf_prefix == 0u &&
    solver->lean_generic_packed_rows32 != NULL &&
    solver->lean_generic_packed_dependencies32 != NULL &&
    solver->lean_generic_packed_publish_mailbox != NULL &&
    solver->lean_generic_packed_publish_begin != NULL &&
    solver->lean_parallel_owner_frontier != NULL;
  const int generic_packed_capable =
    generic_packed_row_worker || scaled_generic_packed_row_worker;
  int generic_packed_sample_candidate = 0;
  if (generic_packed_capable) {
    if (solver->lean_generic_packed_choice == 2) {
      shared->lean_generic_packed_mode = 1;
    } else if (solver->lean_generic_packed_choice == 0) {
      /* Interleave the alternatives to cancel warmup and slow drift. */
      static const unsigned char candidate_sequence[8] = {
        0u, 1u, 1u, 0u, 1u, 0u, 0u, 1u
      };
      const unsigned int sample = solver->lean_generic_packed_sample_count;
      generic_packed_sample_candidate = candidate_sequence[sample & 7u] != 0u;
      shared->lean_generic_packed_mode = generic_packed_sample_candidate;
    }
  }
  if (compact_amf_two_block_exact_worker ||
      (generic_packed_capable && shared->lean_generic_packed_mode)) {
    for (int mailbox = 0;
         mailbox < thread_count * thread_count; ++mailbox) {
      atomic_store_explicit(
        &solver->lean_parallel_owner_frontier[mailbox].generation, 0u,
        memory_order_relaxed);
    }
  }
  const int compact_parallel_refine_copy =
    compact_amf_two_block_exact_worker &&
    solver->nnz >= (UF_long)(1024 * thread_count);
  const int parallel_refine_copy =
    !solver->parallel_refine_values_copied &&
    solver->solve_refine_values != NULL &&
    solver->solve_refine_values != numeric_values &&
    (solver->nnz >= (UF_long)(16384 * thread_count) ||
     compact_parallel_refine_copy);
  shared->refine_copy_input =
    parallel_refine_copy ? numeric_values : NULL;
  shared->refine_copy_output =
    parallel_refine_copy ? solver->solve_refine_values : NULL;
  shared->refine_copy_nnz = parallel_refine_copy ? solver->nnz : 0u;
  for (int i = 0; i < thread_count; ++i) {
    pool->workers[i].shared = shared;
    pool->workers[i].x = scratch[i];
  }
  const int time_generic_packed_sample =
    generic_packed_capable && solver->lean_generic_packed_choice == 0;
  const double packed_sample_begin =
    time_generic_packed_sample ? kls_now_seconds() : 0.0;

  /* This path targets factors whose entire numeric phase is measured in
     microseconds.  The generic dispatcher puts the caller through a
     condition-variable wait; keep the caller spinning alongside the already
     spinning persistent workers instead. */
  kls_egraph_pool_dispatch_and_spin_wait(pool, shared, thread_count);
  const double packed_sample_end =
    time_generic_packed_sample ? kls_now_seconds() : 0.0;
  const int invalid = shared->invalid;
  const int singular = shared->singular;
  const UF_long numerical_rank = shared->numerical_rank;
  const UF_long singular_col = shared->singular_col;
  if (time_generic_packed_sample && !invalid && !singular) {
    const unsigned int sample = solver->lean_generic_packed_sample_count;
    const unsigned int slot = sample >> 1u;
    double *samples = generic_packed_sample_candidate
      ? solver->lean_generic_packed_candidate_seconds
      : solver->lean_generic_packed_legacy_seconds;
    if (slot < 8u) {
      samples[slot] = packed_sample_end - packed_sample_begin;
      solver->lean_generic_packed_sample_count = sample + 1u;
    }
    if (solver->lean_generic_packed_sample_count == 16u) {
      double paired_ratio[8];
      for (unsigned int i = 0u; i < 8u; ++i) {
        paired_ratio[i] = solver->lean_generic_packed_candidate_seconds[i] /
          solver->lean_generic_packed_legacy_seconds[i];
      }
      for (unsigned int i = 1u; i < 8u; ++i) {
        for (unsigned int j = i;
             j > 0u && paired_ratio[j] < paired_ratio[j - 1u]; --j) {
          const double tmp = paired_ratio[j];
          paired_ratio[j] = paired_ratio[j - 1u];
          paired_ratio[j - 1u] = tmp;
        }
      }
      /* Require a repeatable margin so noise defaults to established code.
         This is deliberately wider than timer noise on sub-microsecond rows. */
      solver->lean_generic_packed_choice =
        0.5 * (paired_ratio[3] + paired_ratio[4]) <= 0.985 ? 2 : 1;
    }
  }
  const int row_values_mode = shared->lean_row_values_mode;
  solver->fused_refactor_solve_computed =
    !invalid && !singular && shared->lean_forward_solve_work != NULL;
  solver->fused_refactor_solve_row_values =
    solver->fused_refactor_solve_computed && row_values_mode;
  solver->compact_amf_two_block_exact_recip_fresh =
    !invalid && !singular && compact_amf_two_block_exact_worker;
  if (!invalid && !singular && generic_i16ptr_grouped_recip_worker) {
    solver->i32solve_udiag_recip_fresh = 1;
  }
  if (!invalid && !singular && parallel_refine_copy) {
    solver->parallel_refine_values_copied = 1;
  }
  shared->refine_copy_input = NULL;
  shared->refine_copy_output = NULL;
  shared->refine_copy_nnz = 0u;
  solver->lean_parallel_scratch_clean =
    !invalid && !singular;
  solver->lean_compact_match_row_factor_active =
    row_values_mode && !invalid && !singular;
  shared->lean_refactor_mode = 0;
  shared->lean_compact_match_mode = 0;
  shared->lean_symmetric_scalar_fringe_mode = 0;
  shared->lean_snode_avx512_enabled = 0;
  shared->lean_parallel_offdiag_mode = 0;
  shared->lean_grouped_done_mode = 0;
  shared->lean_row_values_mode = 0;
  shared->lean_generic_packed_mode = 0;
  shared->lean_fused_scale_mode = 0;
  shared->lean_forward_solve_work = NULL;
  shared->lean_forward_solve_rhs = NULL;
  shared->lean_forward_solve_begin = 0u;
  shared->udiag_recip = NULL;
  shared->row_refactor_defer_value_scatter = 0;
  shared->value_prep_input = NULL;
  shared->value_prep_output = NULL;
  shared->value_prep_input_to_csc = NULL;
  shared->value_prep_nnz = 0u;
  shared->lean_done = NULL;
  shared->lean_rows = NULL;
  shared->pipeline_generation = 0u;
  pthread_mutex_unlock(&shared->lock);

  if (invalid) {
    solver->common.status = TRILINOS_KLU_INVALID;
    solver->egraph_worker_scratch_dirty = 1;
    return 0;
  }
  if (singular) {
    solver->common.status = TRILINOS_KLU_SINGULAR;
    solver->common.numerical_rank = numerical_rank;
    solver->common.singular_col = singular_col;
    return solver->common.halt_if_singular ? 0 : 1;
  }
  solver->common.status = TRILINOS_KLU_OK;
  return 1;
}

/* Lean serial up-looking refactorization over the row mirrors: one flat
   two-stream walk per row (ST Alg.1 refactor mode), no group machinery,
   no workspaces beyond Xwork, values written straight through the
   LUbx-pointer mirrors.  Built for the microsecond factors where the
   full row engine's per-row overheads lose 2.4-3.8x to the column
   engines and the column walk itself runs ~2x off CKTSO's rate.
   Returns 1 ok, 0 factor failure (status set), -1 decline. */
static int kls_lean_row_refactor_numeric(kls_solver *solver,
                                         double *numeric_values) {
  if (solver == NULL || solver->numeric_is_predicted) {
    return -1;
  }
  const int scaled = solver->common.scale > 0;
  if (numeric_values == NULL || solver->symbolic == NULL ||
      solver->numeric == NULL || solver->symbolic->nblocks == 0u ||
      (scaled && (solver->numeric->Rs == NULL ||
                  solver->numeric->Pnum == NULL)) ||
      (!scaled && solver->numeric->Rs != NULL) ||
      solver->numeric->Udiag == NULL || solver->numeric->Xwork == NULL ||
      !kls_build_row_refactor_pattern(solver, 1)) {
    return -1;
  }
  const int generic_supernode_candidate =
    kls_generic_hoisted_snode_worker_candidate(solver);
  if (generic_supernode_candidate) {
    kls_lean_row_refactor_ensure_snode_runs(solver);
  }
  const int generic_supernode_values =
    kls_generic_hoisted_snode_worker_capable(solver);
  const int generic_i16_capable =
    solver->n <= (UF_long)UINT16_MAX &&
    solver->nnz <= (UF_long)UINT16_MAX &&
    solver->numeric->lnz <= (UF_long)UINT16_MAX &&
    solver->numeric->unz <= (UF_long)UINT16_MAX &&
    !generic_supernode_values;
  const int generic_snode_i16_capable =
    generic_supernode_values &&
    /* Building both compact column streams is a one-time O(L+U) pass.  The
       measured crossover is above the generic recurring-workload floor, so
       require a lifecycle long enough to amortize that representation. */
    solver->options.expected_refactorizations >= 64 &&
    solver->options.threads > 1 &&
    solver->n <= (UF_long)UINT16_MAX + 1u;
  const int scaled_hoisted_i16_capable =
    scaled && !generic_supernode_values &&
    kls_repeated_update_workload(&solver->options) &&
    solver->n <= (UF_long)UINT16_MAX + 1u;
  if (generic_i16_capable || generic_snode_i16_capable ||
      scaled_hoisted_i16_capable) {
    (void)kls_build_lean_row_i16_indices(solver);
  }
  if (solver->input_to_csc != NULL && solver->values != NULL &&
      solver->row_scale == NULL && solver->col_scale == NULL &&
      !scaled && kls_repeated_update_workload(&solver->options)) {
    (void)kls_build_prepared_value_input_pos(solver);
  }
  if (solver->lean_user_values_active ||
      kls_direct_user_value_maps_capable(solver)) {
    (void)kls_build_compact_match_direct_value_maps(solver);
  }
  if (solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_l_cols == NULL ||
      solver->row_refactor_l_values == NULL ||
      solver->row_refactor_u_ptr == NULL ||
      solver->row_refactor_u_cols == NULL ||
      solver->row_refactor_u_values == NULL ||
      solver->row_refactor_input_ptr == NULL ||
      solver->row_refactor_input_cols == NULL ||
      solver->row_refactor_input_pos == NULL) {
    return -1;
  }
  trilinos_klu_l_common *common = &solver->common;
  trilinos_klu_l_numeric *numeric = solver->numeric;
  const UF_long n = solver->n;
  const int parallel_wanted =
    solver->options.threads > 1 && solver->n >= 800u &&
    solver->common.flops >= 100000.0 &&
    solver->row_refactor_level_max_width >=
      (UF_long)(2 * solver->options.threads);
  int fused_scale_mode = 0;
  if (scaled) {
    const int scale_cached =
      kls_lean_scale_input_is_unchanged(solver, numeric_values);
    fused_scale_mode =
      (!scale_cached && solver->lean_scale_input_state < 0 &&
      parallel_wanted && !generic_supernode_values &&
      solver->lean_snode_run == NULL &&
      solver->row_refactor_input_ptr[n] == solver->nnz &&
      solver->nnz >=
        (UF_long)4096u * (UF_long)solver->options.threads) ? 2 : 0;
    if (!scale_cached && !fused_scale_mode &&
        !kls_scale_from_row_refactor_stream(solver, numeric_values)) {
      if (!trilinos_klu_l_scale((UF_long)common->scale, n,
                                solver->col_ptr, solver->row_idx,
                                numeric_values, numeric->Rs, NULL, common)) {
        return 0;
      }
    }
    if (!scale_cached && !fused_scale_mode) {
      kls_arm_lean_scale_input_cache(solver, numeric_values);
    }
  }
  const int parallel_direct_offdiag =
    (solver->lean_deferred_value_prep_active ||
     (solver->lean_user_values_active &&
      solver->compact_match_offdiag_user_pos != NULL)) &&
    solver->options.threads > 1 &&
    solver->n >= 800u && solver->common.flops >= 100000.0 &&
    solver->row_refactor_level_max_width >=
      (UF_long)(2 * solver->options.threads);
  const int generic_parallel_offdiag =
    !solver->lean_user_values_active &&
    kls_parallel_offdiag_is_profitable(solver) &&
    solver->options.threads > 1 &&
    solver->n >= 800u && solver->common.flops >= 100000.0 &&
    solver->row_refactor_level_max_width >=
      (UF_long)(2 * solver->options.threads);
  if (!parallel_direct_offdiag && !generic_parallel_offdiag &&
      !kls_refresh_row_refactor_offdiag_values(solver, numeric_values)) {
    common->status = TRILINOS_KLU_INVALID;
    return 0;
  }
  common->status = TRILINOS_KLU_OK;
  common->numerical_rank = KLS_KLU_EMPTY;
  common->singular_col = KLS_KLU_EMPTY;
  double *x = (double *)numeric->Xwork;
  memset(x, 0, (size_t)n * sizeof(*x));
  const UF_long *in_ptr = solver->row_refactor_input_ptr;
  const uint16_t *in_ptr16 = solver->row_refactor_input_ptr16;
  const UF_long *in_cols = solver->row_refactor_input_cols;
  const UF_long *in_pos = solver->row_refactor_input_pos;
  const UF_long *l_ptr = solver->row_refactor_l_ptr;
  const uint16_t *l_ptr16 = solver->row_refactor_l_ptr16;
  const UF_long *l_cols = solver->row_refactor_l_cols;
  const uint16_t *l_cols16 = solver->row_refactor_l_cols16;
  double *l_val = solver->row_refactor_l_row_values;
  double **l_lu = solver->row_refactor_l_values;
  const UF_long *u_ptr = solver->row_refactor_u_ptr;
  const uint16_t *u_ptr16 = solver->row_refactor_u_ptr16;
  const UF_long *u_cols = solver->row_refactor_u_cols;
  const uint16_t *u_cols16 = solver->row_refactor_u_cols16;
  double *u_val = solver->row_refactor_u_row_values;
  double **u_lu = solver->row_refactor_u_values;
  double *udiag = (double *)numeric->Udiag;
  if (l_val == NULL || u_val == NULL || l_lu == NULL || u_lu == NULL) {
    return -1;
  }
  const unsigned char *snode_runs = generic_supernode_values
    ? solver->lean_snode_run : NULL;
  const double *rs = scaled ? numeric->Rs : NULL;
  const UF_long *pnum = numeric->Pnum;
  const int pair_mode = solver->lean_pair_active;
  solver->lean_parallel_last_used = 0;
  int parallel_status = -1;
  if (parallel_wanted) {
    parallel_status =
      kls_lean_parallel_refactor_run(
        solver, numeric_values, rs, fused_scale_mode);
    solver->lean_parallel_last_used = parallel_status > 0;
    if (parallel_status == 0 &&
        solver->lean_deferred_value_prep_active) {
      /* A synchronized gather changes the launch order seen by the retained
         dependency executor.  Rare valid factors reject that launch even
         though the gathered internal frame is exact.  Retry from the
         ordinary prepared frame so this optimization remains an optional
         executor fast path, not a new numerical contract or matrix policy. */
      double *prepared = NULL;
      if (prepare_numeric_values(solver, numeric_values, &prepared) != KLS_OK ||
          prepared == NULL) {
        common->status = TRILINOS_KLU_INVALID;
        return 0;
      }
      numeric_values = prepared;
      solver->lean_deferred_value_prep_active = 0;
      common->status = TRILINOS_KLU_OK;
      common->numerical_rank = KLS_KLU_EMPTY;
      common->singular_col = KLS_KLU_EMPTY;
      if (!kls_refresh_row_refactor_offdiag_values(
            solver, numeric_values)) {
        common->status = TRILINOS_KLU_INVALID;
        return 0;
      }
      parallel_status = kls_lean_parallel_refactor_run(
        solver, numeric_values, rs, fused_scale_mode);
      solver->lean_parallel_last_used = parallel_status > 0;
    }
    if (parallel_status == 0) {
      return 0;
    }
  }
  if (parallel_status < 0) {
    if (solver->lean_deferred_value_prep_active) {
      double *prepared = NULL;
      if (prepare_numeric_values(solver, numeric_values, &prepared) != KLS_OK ||
          prepared == NULL) {
        common->status = TRILINOS_KLU_INVALID;
        return 0;
      }
      numeric_values = prepared;
      solver->lean_deferred_value_prep_active = 0;
    }
    if (fused_scale_mode &&
        !trilinos_klu_l_scale((UF_long)common->scale, n,
                              solver->col_ptr, solver->row_idx,
                              numeric_values, numeric->Rs, NULL, common)) {
      return 0;
    }
    /* The coverage-selected descriptor is profitable in the hoisted pool
       executor, but its serial chain fallback adds more bookkeeping than it
       removes.  A resource failure or unselected pool therefore
       falls back to the ordinary scalar row walk. */
    if (generic_supernode_values) {
      snode_runs = NULL;
    }
    if ((parallel_direct_offdiag || generic_parallel_offdiag) &&
        !kls_refresh_row_refactor_offdiag_values(solver, numeric_values)) {
      common->status = TRILINOS_KLU_INVALID;
      return 0;
    }
#define KLS_LEAN_IN_PTR(row)                                               \
  (in_ptr16 != NULL ? (UF_long)in_ptr16[(row)] : in_ptr[(row)])
#define KLS_LEAN_L_PTR(row)                                                \
  (l_ptr16 != NULL ? (UF_long)l_ptr16[(row)] : l_ptr[(row)])
#define KLS_LEAN_U_PTR(row)                                                \
  (u_ptr16 != NULL ? (UF_long)u_ptr16[(row)] : u_ptr[(row)])
#define KLS_LEAN_L_COL(pos)                                                \
  (l_cols16 != NULL ? (UF_long)l_cols16[(pos)] : l_cols[(pos)])
#define KLS_LEAN_U_COL(pos)                                                \
  (u_cols16 != NULL ? (UF_long)u_cols16[(pos)] : u_cols[(pos)])
#define KLS_LEAN_ROW_SCATTER(row, xv)                                      \
  do {                                                                     \
    const double row_rs_inv =                                              \
      rs != NULL ? 1.0 / rs[pnum[(row)]] : 1.0;                            \
    for (UF_long p = KLS_LEAN_IN_PTR((row));                               \
         p < KLS_LEAN_IN_PTR((row) + 1u); ++p) {                           \
      (xv)[in_cols[p]] = numeric_values[in_pos[p]] * row_rs_inv;           \
    }                                                                      \
  } while (0)
#define KLS_LEAN_ROW_KSTEP(pvar, xv)                                       \
  do {                                                                     \
    const UF_long k = KLS_LEAN_L_COL((pvar));                              \
    const double lik = (xv)[k] / udiag[k];                                 \
    (xv)[k] = 0.0;                                                         \
    l_val[(pvar)] = lik;                                                   \
    *l_lu[(pvar)] = lik;                                                   \
    if (lik != 0.0) {                                                      \
      for (UF_long q = KLS_LEAN_U_PTR(k);                                  \
           q < KLS_LEAN_U_PTR(k + 1u); ++q) {                              \
        (xv)[KLS_LEAN_U_COL(q)] -= lik * u_val[q];                         \
      }                                                                    \
    }                                                                      \
    (pvar)++;                                                              \
  } while (0)
#define KLS_LEAN_ROW_SNODE2(pvar, xv)                                      \
  do {                                                                     \
    const UF_long dep0 = KLS_LEAN_L_COL((pvar));                           \
    const UF_long dep1 = dep0 + 1u;                                        \
    const double multiplier0 = (xv)[dep0] / udiag[dep0];                  \
    const double multiplier1 =                                            \
      ((xv)[dep1] - multiplier0 * u_val[KLS_LEAN_U_PTR(dep0)]) /          \
       udiag[dep1];                                                        \
    (xv)[dep0] = 0.0;                                                      \
    (xv)[dep1] = 0.0;                                                      \
    l_val[(pvar)] = multiplier0;                                           \
    l_val[(pvar) + 1u] = multiplier1;                                      \
    *l_lu[(pvar)] = multiplier0;                                           \
    *l_lu[(pvar) + 1u] = multiplier1;                                      \
    const UF_long trailing_begin = KLS_LEAN_U_PTR(dep1);                   \
    const UF_long trailing_end = KLS_LEAN_U_PTR(dep1 + 1u);                \
    const double *restrict u0 = u_val + KLS_LEAN_U_PTR(dep0) + 1u;         \
    const double *restrict u1 = u_val + trailing_begin;                    \
    for (UF_long q = trailing_begin; q < trailing_end; ++q) {             \
      const UF_long target = KLS_LEAN_U_COL(q);                            \
      const UF_long offset = q - trailing_begin;                           \
      double value = (xv)[target];                                         \
      value -= multiplier0 * u0[offset];                                   \
      value -= multiplier1 * u1[offset];                                   \
      (xv)[target] = value;                                                \
    }                                                                      \
    (pvar) += 2u;                                                          \
  } while (0)
#define KLS_LEAN_ROW_STEP(pvar, xv)                                        \
  do {                                                                     \
    const UF_long snode_run =                                              \
      snode_runs != NULL ? (UF_long)snode_runs[(pvar)] : 0u;              \
    if (snode_run == 2u) {                                                 \
      KLS_LEAN_ROW_SNODE2((pvar), (xv));                                   \
    } else if (snode_run > 0u) {                                           \
      (void)kls_lean_row_refactor_snode_step_chain(                        \
        (pvar), snode_run, KLS_LEAN_L_COL((pvar)), (xv), l_val, l_lu,     \
        u_ptr, u_cols, u_val, udiag);                                      \
      (pvar) += snode_run;                                                 \
    } else {                                                               \
      KLS_LEAN_ROW_KSTEP((pvar), (xv));                                    \
    }                                                                      \
  } while (0)
#define KLS_LEAN_ROW_FINISH(row, xv)                                       \
  do {                                                                     \
    const double piv = (xv)[(row)];                                        \
    (xv)[(row)] = 0.0;                                                     \
    udiag[(row)] = piv;                                                    \
    if (piv == 0.0) {                                                      \
      common->status = TRILINOS_KLU_SINGULAR;                              \
      if (common->numerical_rank == KLS_KLU_EMPTY) {                       \
        common->numerical_rank = (row);                                    \
        common->singular_col = (row);                                      \
      }                                                                    \
    }                                                                      \
    for (UF_long q = KLS_LEAN_U_PTR((row));                                \
         q < KLS_LEAN_U_PTR((row) + 1u); ++q) {                            \
      const UF_long j = KLS_LEAN_U_COL(q);                                 \
      const double v = (xv)[j];                                            \
      u_val[q] = v;                                                        \
      *u_lu[q] = v;                                                        \
      (xv)[j] = 0.0;                                                       \
    }                                                                      \
  } while (0)
  if (pair_mode) {
    /* Keep the numeric row order: the fixed L pattern is sufficient for
       the scalar walk but is not a safe topological schedule for every
       transformed input frame.  Adjacent rows can still be multiplexed
       when the second row has no direct dependency on the first. */
    if (solver->lean_row_x2 == NULL) {
      solver->lean_row_x2 =
        (double *)calloc((size_t)n, sizeof(*solver->lean_row_x2));
    }
    double *x2 = solver->lean_row_x2;
    if (x2 == NULL) {
      return -1;
    }
    memset(x2, 0, (size_t)n * sizeof(*x2));
    UF_long i1 = 0u;
    while (i1 < n) {
      const UF_long i2 = i1 + 1u;
      int independent = i2 < n;
      if (independent) {
        for (UF_long p = KLS_LEAN_L_PTR(i2);
             p < KLS_LEAN_L_PTR(i2 + 1u); ++p) {
          if (KLS_LEAN_L_COL(p) == i1) {
            independent = 0;
            break;
          }
        }
      }
      if (independent) {
        KLS_LEAN_ROW_SCATTER(i1, x);
        KLS_LEAN_ROW_SCATTER(i2, x2);
        UF_long p1 = KLS_LEAN_L_PTR(i1);
        const UF_long e1 = KLS_LEAN_L_PTR(i1 + 1u);
        UF_long p2 = KLS_LEAN_L_PTR(i2);
        const UF_long e2 = KLS_LEAN_L_PTR(i2 + 1u);
        while (p1 < e1 && p2 < e2) {
          KLS_LEAN_ROW_STEP(p1, x);
          KLS_LEAN_ROW_STEP(p2, x2);
        }
        while (p1 < e1) {
          KLS_LEAN_ROW_STEP(p1, x);
        }
        while (p2 < e2) {
          KLS_LEAN_ROW_STEP(p2, x2);
        }
        KLS_LEAN_ROW_FINISH(i1, x);
        KLS_LEAN_ROW_FINISH(i2, x2);
        i1 += 2u;
      } else {
        KLS_LEAN_ROW_SCATTER(i1, x);
        UF_long p1 = KLS_LEAN_L_PTR(i1);
        const UF_long e1 = KLS_LEAN_L_PTR(i1 + 1u);
        while (p1 < e1) {
          KLS_LEAN_ROW_STEP(p1, x);
        }
        KLS_LEAN_ROW_FINISH(i1, x);
        i1++;
      }
      if (common->status == TRILINOS_KLU_SINGULAR &&
          common->halt_if_singular) {
        return 0;
      }
    }
  } else {
    for (UF_long i = 0; i < n; ++i) {
      KLS_LEAN_ROW_SCATTER(i, x);
      UF_long p = KLS_LEAN_L_PTR(i);
      const UF_long pe = KLS_LEAN_L_PTR(i + 1u);
      while (p < pe) {
        KLS_LEAN_ROW_STEP(p, x);
      }
      KLS_LEAN_ROW_FINISH(i, x);
      if (common->status == TRILINOS_KLU_SINGULAR &&
          common->halt_if_singular) {
        return 0;
      }
    }
  }
#undef KLS_LEAN_ROW_SCATTER
#undef KLS_LEAN_ROW_KSTEP
#undef KLS_LEAN_ROW_SNODE2
#undef KLS_LEAN_ROW_STEP
#undef KLS_LEAN_ROW_FINISH
#undef KLS_LEAN_IN_PTR
#undef KLS_LEAN_L_PTR
#undef KLS_LEAN_U_PTR
#undef KLS_LEAN_L_COL
#undef KLS_LEAN_U_COL
  }
  if (scaled && fused_scale_mode < 2 &&
      !kls_parallel_refactor_permute_scale(solver)) {
    common->status = TRILINOS_KLU_INVALID;
    return 0;
  }
  /* Dual stores: the flat mirrors serve this walk's own consumes (each
     entry is re-read flops/fill times, so the reuse pays for the extra
     store) while LUbx stays current for the faster column solve paths. */
  solver->row_refactor_values_dirty =
    solver->lean_compact_match_row_factor_active;
  solver->stats.row_refactor_values_dirty =
    solver->lean_compact_match_row_factor_active;
  solver->row_refactor_values_ready = 0;
  solver->row_refactor_solve_direct_ready = 0;
  solver->row_refactor_solve_validated = 0;
  return 1;
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

static int kls_publish_row_refactor_values(kls_solver *solver) {
  if (solver == NULL || !solver->row_refactor_values_dirty) {
    return 1;
  }
  if (!kls_scatter_row_refactor_l_values(solver) ||
      !kls_scatter_row_refactor_u_values(solver)) {
    solver->common.status = TRILINOS_KLU_INVALID;
    return 0;
  }
  solver->row_refactor_values_dirty = 0;
  solver->lean_compact_match_row_factor_active = 0;
  solver->row_refactor_solve_direct_ready = 0;
  solver->stats.row_refactor_values_dirty = 0;
  return 1;
}

static int kls_copy_row_refactor_values_from_numeric(kls_solver *solver) {
  if (solver == NULL) {
    return 0;
  }
  const UF_long l_nnz =
    solver->row_refactor_l_ptr != NULL ? solver->row_refactor_l_ptr[solver->n]
                                       : 0u;
  if (l_nnz > 0u &&
      (solver->row_refactor_l_values == NULL ||
       solver->row_refactor_l_row_values == NULL)) {
    return 0;
  }
  for (UF_long p = 0; p < l_nnz; ++p) {
    if (solver->row_refactor_l_values[p] == NULL) {
      return 0;
    }
    solver->row_refactor_l_row_values[p] = *solver->row_refactor_l_values[p];
  }

  const UF_long u_nnz =
    solver->row_refactor_u_ptr != NULL ? solver->row_refactor_u_ptr[solver->n]
                                       : 0u;
  if (u_nnz > 0u &&
      (solver->row_refactor_u_values == NULL ||
       solver->row_refactor_u_row_values == NULL)) {
    return 0;
  }
  for (UF_long p = 0; p < u_nnz; ++p) {
    if (solver->row_refactor_u_values[p] == NULL) {
      return 0;
    }
    solver->row_refactor_u_row_values[p] = *solver->row_refactor_u_values[p];
  }

  solver->row_refactor_values_ready = 1;
  solver->row_refactor_solve_direct_ready = 0;
  solver->row_refactor_solve_validated = 0;
  return 1;
}

static int kls_seed_row_refactor_values_from_numeric(kls_solver *solver) {
  if (solver == NULL || !kls_build_row_refactor_pattern(solver, 0)) {
    return 0;
  }
  return kls_copy_row_refactor_values_from_numeric(solver);
}

static int kls_estimate_row_refactor_lower_bound_work(kls_solver *solver,
                                                      double *work_out) {
  if (work_out != NULL) {
    *work_out = 0.0;
  }
  if (solver == NULL || work_out == NULL || solver->symbolic == NULL ||
      solver->numeric == NULL || solver->symbolic->R == NULL ||
      solver->numeric->Llen == NULL || solver->numeric->Ulen == NULL ||
      !kls_build_refactor_map(solver) ||
      !kls_build_refactor_lu_pointer_cache(solver)) {
    return 0;
  }

  const UF_long n = solver->n;
  UF_long *u_row_counts =
    n > 0u ? (UF_long *)calloc((size_t)n, sizeof(*u_row_counts)) : NULL;
  if (n > 0u && u_row_counts == NULL) {
    return 0;
  }

  double work = 0.0;
  for (UF_long k = 0; k < n; ++k) {
    UF_long input_begin = solver->refactor_col_ptr[k];
    if (solver->symbolic->nblocks > 1u) {
      if (solver->refactor_block_start == NULL ||
          solver->refactor_block_start[k] < input_begin ||
          solver->refactor_block_start[k] > solver->refactor_col_ptr[k + 1u]) {
        free(u_row_counts);
        return 0;
      }
      input_begin = solver->refactor_block_start[k];
    }
    if (solver->refactor_col_ptr[k + 1u] < input_begin) {
      free(u_row_counts);
      return 0;
    }
    work += (double)(solver->refactor_col_ptr[k + 1u] - input_begin);
  }

  for (UF_long block = 0; block < solver->symbolic->nblocks; ++block) {
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    if (k1 > k2 || k2 > n) {
      free(u_row_counts);
      return 0;
    }
    for (UF_long local_k = 0; local_k < k2 - k1; ++local_k) {
      const UF_long k = k1 + local_k;
      UF_long *ui = solver->refactor_u_indices[k];
      const UF_long len = kls_numeric_u_len_in_block(solver, k1, k2, k);
      if (len > 0u && ui == NULL) {
        free(u_row_counts);
        return 0;
      }
      work += (double)len;
      for (UF_long p = 0; p < len; ++p) {
        const UF_long local_row = ui[p];
        const UF_long row = k1 + local_row;
        if (local_row >= local_k || row >= n) {
          free(u_row_counts);
          return 0;
        }
        u_row_counts[row]++;
      }
    }
  }

  for (UF_long block = 0; block < solver->symbolic->nblocks; ++block) {
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    for (UF_long local_j = 0; local_j < k2 - k1; ++local_j) {
      const UF_long j = k1 + local_j;
      UF_long *li = solver->refactor_l_indices[j];
      const UF_long len = kls_numeric_l_len_in_block(solver, k1, k2, j);
      if (len > 0u && li == NULL) {
        free(u_row_counts);
        return 0;
      }
      const double update_work = 1.0 + (double)u_row_counts[j];
      for (UF_long p = 0; p < len; ++p) {
        const UF_long local_row = li[p];
        const UF_long row = k1 + local_row;
        if (local_row <= local_j || local_row >= k2 - k1 || row >= n) {
          free(u_row_counts);
          return 0;
        }
        work += update_work;
      }
    }
  }

  free(u_row_counts);
  *work_out = work;
  return 1;
}

static void kls_reset_auto_row_refactor_prepare_stats(kls_solver *solver) {
  if (solver == NULL) {
    return;
  }
  solver->row_refactor_auto_native_row_panel = 0;
  solver->row_refactor_auto_lower_bound_work = 0.0;
  solver->row_refactor_auto_lower_bound_rejected = 0;
  solver->row_refactor_auto_pattern_build_failed = 0;
  solver->row_refactor_auto_value_copy_failed = 0;
}

static int kls_prepare_auto_row_refactor_from_numeric(kls_solver *solver) {
  if (solver == NULL) {
    return 0;
  }
  kls_reset_auto_row_refactor_prepare_stats(solver);
  if (solver->refactor_dependency_work > 0.0) {
    double lower_bound_work = 0.0;
    if (kls_estimate_row_refactor_lower_bound_work(solver,
                                                   &lower_bound_work)) {
      solver->row_refactor_auto_lower_bound_work = lower_bound_work;
      double lower_bound_allowance = solver->refactor_dependency_work;
      if (solver->options.threads > 1) {
        /* The bound is serial work; the cooperative row engine recovers
           up to thread_count of it, so keep candidates inside that
           envelope alive and let measured trials decide. */
        lower_bound_allowance *= (double)solver->options.threads;
      }
      if (lower_bound_work > lower_bound_allowance) {
        solver->row_refactor_auto_enabled = 0;
        solver->row_refactor_auto_native_row_panel = 0;
        solver->row_refactor_values_ready = 0;
        solver->row_refactor_solve_direct_ready = 0;
        solver->row_refactor_solve_validated = 0;
        solver->row_refactor_auto_lower_bound_rejected = 1;
        return 1;
      }
    }
  }
  if (!kls_build_refactor_map(solver) ||
      !kls_build_row_refactor_pattern(solver, 0)) {
    solver->row_refactor_auto_enabled = 0;
    solver->row_refactor_auto_native_row_panel = 0;
    solver->row_refactor_values_ready = 0;
    solver->row_refactor_solve_direct_ready = 0;
    solver->row_refactor_solve_validated = 0;
    solver->row_refactor_auto_pattern_build_failed = 1;
    return 0;
  }
  solver->row_refactor_auto_enabled = 1;
  solver->row_refactor_auto_native_row_panel =
    kls_native_row_panel_structural_auto_allows(solver) ? 1 : 0;
  if (!kls_auto_row_refactor_cost_allows(solver)) {
    solver->row_refactor_auto_native_row_panel = 0;
    solver->row_refactor_values_ready = 0;
    solver->row_refactor_solve_direct_ready = 0;
    solver->row_refactor_solve_validated = 0;
    return 1;
  }
  if (solver->row_refactor_values_ready &&
      !solver->row_refactor_values_dirty) {
    solver->row_refactor_solve_direct_ready = 0;
    solver->row_refactor_solve_validated = 0;
    return 1;
  }
  if (!kls_copy_row_refactor_values_from_numeric(solver)) {
    solver->row_refactor_auto_native_row_panel = 0;
    solver->row_refactor_auto_value_copy_failed = 1;
    return 0;
  }
  return 1;
}

static int kls_egraph_refactor_is_eligible(const kls_solver *solver);

static void kls_maybe_prepare_model_row_refactor_from_numeric(
  kls_solver *solver,
  double *elapsed) {
  if (solver == NULL || elapsed == NULL ||
      kls_row_refactor_env_disabled() ||
      (!kls_row_refactor_env_enabled() &&
       !kls_auto_row_refactor_high_work_candidate(solver)) ||
      solver->options.threads <= 1 ||
      solver->row_refactor_auto_enabled ||
      solver->row_refactor_values_dirty ||
      solver->numeric == NULL ||
      solver->common.status < TRILINOS_KLU_OK ||
      solver->common.status == TRILINOS_KLU_SINGULAR) {
    return;
  }

  /* A predicted numeric is not eligible for the automatic row-engine
     tournament unless the caller explicitly enables that experiment (see
     kls_row_refactor_acceptance_structurally_ready).  Building its complete
     row mirror here anyway used to strand the numeric in an undecided state:
     no row refactor could ever be sampled, while every solve consumed the
     unmeasured mirror and every column refactor paid to keep it coherent.
     Require the same eligibility at preparation time so an optional
     representation cannot become the default merely as a setup side effect. */
  if (solver->numeric_is_predicted) {
    return;
  }

  if (solver->refactor_level_ptr != NULL && solver->n > 0u &&
      solver->refactor_dependency_work > 0.0 &&
      solver->refactor_pipeline_work >=
        0.95 * solver->refactor_dependency_work &&
      solver->refactor_level_max_width >= solver->n / 4u &&
      solver->refactor_pipeline_column_count <= solver->n / 32u) {
    /* The retained column graph already certifies a broad executor whose
       bounded dependency tail contains nearly all arithmetic but fewer than
       one column in 32, while an ordinary level already exposes at least a
       quarter of all columns.  A row transpose cannot expose material
       additional concurrency here: it duplicates the complete factor
       streams merely to rediscover the same narrow tail.  This schedule
       certificate is intentionally independent of final EGraph admission,
       which is decided after row preparation.  Skip that O(fill)
       representation build from realized coverage, width, and work—not from
       an input family or dimension.  This also protects finite repeated-
       update lifecycles where a rejected row mirror costs several steady
       EGraph generations before its first timing sample could repay anything. */
    return;
  }

  /* Row preparation consumes the retained numeric's validated L/U index and
     value streams, not the algorithm that originally allocated them.  Once
     an origin is eligible, let the ordinary work model, lower-bound
     comparison, complete pattern builder, timed engine tournament, and
     numeric residual checks arbitrate it. */

  int recommends_parallel = 0;
  kls_compute_parallel_model(solver, 1, NULL, NULL, &recommends_parallel);
  int task_flow_recommends_parallel = 0;
  kls_compute_nicslu_task_flow_model(
    solver, solver->options.threads, NULL, NULL, NULL, NULL, NULL,
    &task_flow_recommends_parallel);
  if (!recommends_parallel && !task_flow_recommends_parallel) {
    return;
  }
  solver->stats.row_refactor_auto_model_recommended = 1;
  if (solver->refactor_level_ptr == NULL) {
    /* No dependency schedule: the egraph and scheduled column engines
       cannot run off it here, so the steady path is mapped/serial —
       with an explicit egraph-eligibility check as defense in depth
       against transiently-absent schedules.  Admit the prepare only
       when no engine already claimed the row (b9_c1's dense-help
       groups sit at ~100us; an unconditional re-enable displaced them
       8x in its tail cert) and the factor is large enough for the
       alternating acceptance trials to amortize.  The fill floor also
       structurally excludes every lean-consult row (lean caps at 1M
       fill).  b2383 class: serial steady 240-305ms, measured row
       acceptance 156ms — the deferred-consult wipe used to strand it
       on serial with no re-prepare. */
    if (solver->row_refactor_group_count > 0 ||
        kls_egraph_refactor_is_eligible(solver) ||
        solver->numeric->lnz + solver->numeric->unz < 2000000) {
      return;
    }
  } else if (solver->refactor_dependency_work < 1.0e7) {
    /* Below ~1e7 dependency work the plan build (~4ms on rajat03)
       exceeds anything the row engine could win back, and the timed
       acceptance has never adopted at that scale; the adopters all
       measure >=1e8. */
    return;
  }

  solver->stats.row_refactor_auto_model_attempted = 1;
  const double start = kls_now_seconds();
  (void)kls_prepare_auto_row_refactor_from_numeric(solver);
  *elapsed += kls_now_seconds() - start;
  solver->stats.row_refactor_auto_model_accepted =
    solver->row_refactor_auto_enabled &&
    solver->row_refactor_values_ready &&
    kls_auto_row_refactor_cost_allows(solver);
}

static void kls_maybe_reseed_auto_row_refactor_values(kls_solver *solver,
                                                      double *elapsed) {
  if (solver == NULL || elapsed == NULL ||
      !solver->row_refactor_auto_enabled ||
      solver->row_refactor_values_ready ||
      solver->row_accept_decision != 0 ||
      !kls_auto_row_refactor_cost_allows(solver) ||
      solver->numeric == NULL ||
      solver->common.status < TRILINOS_KLU_OK ||
      solver->common.status == TRILINOS_KLU_SINGULAR) {
    return;
  }
  if (kls_row_refactor_acceptance_structurally_ready(solver) &&
      solver->stats.last_refactor_path != KLS_REFACTOR_PATH_ROW) {
    /* During the paired engine audit, a column refactor must be followed by
       the column solve it actually enables.  Rebuilding the optional row
       value mirror here both charges an unrelated O(fill) copy and, for a
       column engine that refreshed a representation-specific numeric, can
       make the following solve consume a row replica that was not produced
       by the row engine.  Leave the column numeric published.  The next row
       trial overwrites every row-factor value directly and therefore needs
       the retained pattern, not a value seed from its competitor. */
    return;
  }
  const double start = kls_now_seconds();
  (void)kls_seed_row_refactor_values_from_numeric(solver);
  *elapsed += kls_now_seconds() - start;
}

static KLS_ALWAYS_INLINE int kls_row_refactor_solve_uses_direct_values(
  const kls_solver *solver) {
  return solver != NULL &&
         !solver->row_refactor_values_ready &&
         solver->row_refactor_solve_direct_ready &&
         !solver->row_refactor_values_dirty;
}

static int kls_row_refactor_solve_has_values(const kls_solver *solver,
                                             int upper) {
  if (solver == NULL) {
    return 0;
  }
  const UF_long *ptr =
    upper ? solver->row_refactor_u_ptr : solver->row_refactor_l_ptr;
  const UF_long nnz = ptr != NULL ? ptr[solver->n] : 0u;
  if (nnz == 0u) {
    return 1;
  }
  if (kls_row_refactor_solve_uses_direct_values(solver)) {
    double **values =
      upper ? solver->row_refactor_u_values : solver->row_refactor_l_values;
    return values != NULL;
  }
  if (!solver->row_refactor_values_ready) {
    return 0;
  }
  const double *row_values =
    upper ? solver->row_refactor_u_row_values
          : solver->row_refactor_l_row_values;
  return row_values != NULL;
}

static KLS_ALWAYS_INLINE double kls_row_refactor_l_solve_value(
  const kls_solver *solver,
  UF_long p) {
  if (kls_row_refactor_solve_uses_direct_values(solver)) {
    return *solver->row_refactor_l_values[p];
  }
  return solver->row_refactor_l_row_values[p];
}

static KLS_ALWAYS_INLINE double kls_row_refactor_u_solve_value(
  const kls_solver *solver,
  UF_long p) {
  if (kls_row_refactor_solve_uses_direct_values(solver)) {
    return *solver->row_refactor_u_values[p];
  }
  return solver->row_refactor_u_row_values[p];
}

static const double *kls_row_refactor_compact_panel_const(
  const kls_solver *solver,
  UF_long group,
  UF_long entry_count) {
  if (solver == NULL ||
      solver->row_refactor_group_compact_panel_begin == NULL ||
      solver->row_refactor_compact_panel_values == NULL ||
      group >= solver->row_refactor_group_count || entry_count == 0u) {
    return NULL;
  }
  const UF_long begin = solver->row_refactor_group_compact_panel_begin[group];
  const UF_long end =
    solver->row_refactor_group_compact_panel_begin[group + 1u];
  if (end < begin || end - begin != entry_count ||
      end > solver->row_refactor_compact_dense_panel_persistent_entries) {
    return NULL;
  }
  return solver->row_refactor_compact_panel_values + begin;
}

static int kls_row_refactor_try_compact_panel_scalar_update(
  kls_egraph_refactor_shared *shared,
  double *x,
  UF_long current_group,
  UF_long dep,
  double lij) {
  kls_solver *solver = shared != NULL ? shared->solver : NULL;
  if (solver == NULL || x == NULL ||
      dep >= solver->n ||
      solver->row_refactor_row_group == NULL ||
      solver->row_refactor_group_ptr == NULL ||
      solver->row_refactor_group_kind == NULL ||
      solver->row_refactor_group_trailing_len == NULL ||
      solver->row_refactor_group_trailing_begin == NULL ||
      solver->row_refactor_compact_panel_valid == NULL ||
      solver->row_refactor_u_ptr == NULL ||
      solver->row_refactor_u_cols == NULL) {
    return 0;
  }
  const UF_long group = solver->row_refactor_row_group[dep];
  if (group >= solver->row_refactor_group_count ||
      (current_group != KLS_KLU_EMPTY && group == current_group) ||
      solver->row_refactor_group_kind[group] !=
        KLS_ROW_REFACTOR_GROUP_DENSE) {
    return 0;
  }
  const UF_long row_begin = solver->row_refactor_group_ptr[group];
  const UF_long row_end = solver->row_refactor_group_ptr[group + 1u];
  if (row_begin >= row_end || dep < row_begin || dep >= row_end) {
    return 0;
  }
  const UF_long width = row_end - row_begin;
  const UF_long local_dep = dep - row_begin;
  if (kls_row_refactor_compact_panel_valid_prefix(solver, group) <=
        local_dep ||
      width > (UF_long)(SIZE_MAX / sizeof(double)) / width) {
    return 0;
  }
  const UF_long dense_entries = width * width;
  const UF_long trailing_len = solver->row_refactor_group_trailing_len[group];
  if (trailing_len > 0u &&
      width > ((UF_long)(SIZE_MAX / sizeof(double)) - dense_entries) /
                trailing_len) {
    return 0;
  }
  const UF_long panel_entries = dense_entries + width * trailing_len;
  const double *dense_panel =
    kls_row_refactor_compact_panel_const(solver, group, panel_entries);
  if (dense_panel == NULL) {
    return 0;
  }

  const UF_long u_begin = solver->row_refactor_u_ptr[dep];
  const UF_long u_end = solver->row_refactor_u_ptr[dep + 1u];
  const UF_long dense_len = row_end - dep - 1u;
  if (u_end < u_begin || u_end - u_begin != dense_len + trailing_len) {
    return 0;
  }
  if (dense_len + trailing_len == 0u) {
    return 0;
  }
  const UF_long *u_cols = solver->row_refactor_u_cols + u_begin;
  for (UF_long offset = 0; offset < dense_len; ++offset) {
    const UF_long col = u_cols[offset];
    if (col != dep + 1u + offset || col >= row_end) {
      return 0;
    }
  }
  const UF_long trailing_begin =
    solver->row_refactor_group_trailing_begin[group];
  const UF_long row_end_u = solver->row_refactor_u_ptr[row_end];
  if (trailing_len > 0u &&
      (trailing_begin > row_end_u ||
       trailing_len > row_end_u - trailing_begin)) {
    return 0;
  }
  for (UF_long offset = 0; offset < trailing_len; ++offset) {
    const UF_long col = u_cols[dense_len + offset];
    if (solver->row_refactor_u_cols[trailing_begin + offset] != col) {
      return 0;
    }
  }

  const double *row_dense_panel = dense_panel + local_dep * width;
  for (UF_long offset = 0; offset < dense_len; ++offset) {
    const UF_long col = u_cols[offset];
    x[col] -= lij * row_dense_panel[local_dep + 1u + offset];
  }
  if (trailing_len > 0u) {
    const double *trailing_panel =
      dense_panel + dense_entries + local_dep * trailing_len;
    for (UF_long offset = 0; offset < trailing_len; ++offset) {
      const UF_long col = u_cols[dense_len + offset];
      x[col] -= lij * trailing_panel[offset];
    }
  }
  kls_record_row_refactor_compact_panel_scalar_update(
    solver, dense_len + trailing_len);
  return 1;
}

static int kls_row_refactor_compact_panel_solve_value(
  const kls_solver *solver,
  int upper,
  UF_long row,
  UF_long p,
  double *value_out) {
  if (value_out != NULL) {
    *value_out = 0.0;
  }
  if (solver == NULL || value_out == NULL ||
      kls_row_refactor_solve_uses_direct_values(solver) ||
      row >= solver->n ||
      !solver->row_refactor_values_ready ||
      solver->row_refactor_row_group == NULL ||
      solver->row_refactor_group_ptr == NULL ||
      solver->row_refactor_group_kind == NULL ||
      solver->row_refactor_group_trailing_len == NULL ||
      solver->row_refactor_group_trailing_begin == NULL ||
      solver->row_refactor_compact_panel_valid == NULL ||
      solver->row_refactor_u_ptr == NULL ||
      solver->row_refactor_u_cols == NULL) {
    return 0;
  }
  const UF_long group = solver->row_refactor_row_group[row];
  if (group >= solver->row_refactor_group_count ||
      solver->row_refactor_group_kind[group] !=
        KLS_ROW_REFACTOR_GROUP_DENSE) {
    return 0;
  }
  const int shape_proven =
    solver->row_refactor_group_shape_valid != NULL &&
    solver->row_refactor_group_shape_valid[group] != 0u;
  const UF_long row_begin = solver->row_refactor_group_ptr[group];
  const UF_long row_end = solver->row_refactor_group_ptr[group + 1u];
  if (row_begin >= row_end || row < row_begin || row >= row_end) {
    return 0;
  }
  const UF_long width = row_end - row_begin;
  if (kls_row_refactor_compact_panel_valid_prefix(solver, group) < width ||
      width > (UF_long)(SIZE_MAX / sizeof(double)) / width) {
    return 0;
  }
  const UF_long dense_entries = width * width;
  const UF_long trailing_len = solver->row_refactor_group_trailing_len[group];
  if (trailing_len > 0u &&
      width > ((UF_long)(SIZE_MAX / sizeof(double)) - dense_entries) /
                trailing_len) {
    return 0;
  }
  const UF_long panel_entries = dense_entries + width * trailing_len;
  const double *dense_panel =
    kls_row_refactor_compact_panel_const(solver, group, panel_entries);
  if (dense_panel == NULL) {
    return 0;
  }
  const UF_long local_row = row - row_begin;
  const double *row_dense_panel = dense_panel + local_row * width;
  const double *trailing_panel =
    trailing_len > 0u ? dense_panel + dense_entries : NULL;

  if (!upper) {
    if (solver->row_refactor_l_ptr == NULL ||
        solver->row_refactor_l_cols == NULL ||
        p < solver->row_refactor_l_ptr[row] ||
        p >= solver->row_refactor_l_ptr[row + 1u]) {
      return 0;
    }
    const UF_long dep = solver->row_refactor_l_cols[p];
    if (shape_proven) {
      if (solver->row_refactor_l_internal_ptr == NULL) {
        return 0;
      }
      const UF_long internal = solver->row_refactor_l_internal_ptr[row];
      if (p < internal || dep < row_begin || dep >= row) {
        return 0;
      }
      *value_out = row_dense_panel[dep - row_begin];
      return 1;
    }
    if (dep < row_begin || dep >= row) {
      return 0;
    }
    *value_out = row_dense_panel[dep - row_begin];
    return 1;
  }

  if (p < solver->row_refactor_u_ptr[row] ||
      p >= solver->row_refactor_u_ptr[row + 1u]) {
    return 0;
  }
  const UF_long offset = p - solver->row_refactor_u_ptr[row];
  const UF_long dense_len = row_end - row - 1u;
  const UF_long col = solver->row_refactor_u_cols[p];
  if (offset < dense_len) {
    if (!shape_proven && (col != row + 1u + offset || col >= row_end)) {
      return 0;
    }
    *value_out = row_dense_panel[local_row + 1u + offset];
    return 1;
  }
  const UF_long trailing_offset = offset - dense_len;
  if (trailing_offset >= trailing_len || trailing_panel == NULL) {
    return 0;
  }
  const UF_long trailing_begin =
    solver->row_refactor_group_trailing_begin[group];
  const UF_long row_end_u = solver->row_refactor_u_ptr[row_end];
  if (trailing_begin > row_end_u ||
      solver->row_refactor_group_trailing_len[group] >
        row_end_u - trailing_begin) {
    return 0;
  }
  const UF_long trailing_index = trailing_begin + trailing_offset;
  if (trailing_index >= row_end_u ||
      (!shape_proven && solver->row_refactor_u_cols[trailing_index] != col)) {
    return 0;
  }
  *value_out = trailing_panel[local_row * trailing_len + trailing_offset];
  return 1;
}

static double kls_row_refactor_l_solve_row_value(
  const kls_solver *solver,
  UF_long row,
  UF_long p) {
  double value = 0.0;
  if (kls_row_refactor_compact_panel_solve_value(solver, 0, row, p,
                                                 &value)) {
    return value;
  }
  return kls_row_refactor_l_solve_value(solver, p);
}

static double kls_row_refactor_u_solve_row_value(
  const kls_solver *solver,
  UF_long row,
  UF_long p) {
  double value = 0.0;
  if (kls_row_refactor_compact_panel_solve_value(solver, 1, row, p,
                                                 &value)) {
    return value;
  }
  return kls_row_refactor_u_solve_value(solver, p);
}

static int kls_row_refactor_compact_panel_solve_group(
  const kls_solver *solver,
  UF_long group,
  const double **dense_panel_out,
  UF_long *row_begin_out,
  UF_long *row_end_out,
  UF_long *trailing_len_out) {
  KLS_SET_OPTIONAL_OUTPUT(dense_panel_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(row_begin_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(row_end_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(trailing_len_out, 0);
  if (solver == NULL ||
      kls_row_refactor_solve_uses_direct_values(solver) ||
      !solver->row_refactor_values_ready ||
      solver->row_refactor_group_ptr == NULL ||
      solver->row_refactor_group_kind == NULL ||
      solver->row_refactor_group_trailing_len == NULL ||
      solver->row_refactor_compact_panel_valid == NULL ||
      group >= solver->row_refactor_group_count ||
      solver->row_refactor_group_kind[group] !=
        KLS_ROW_REFACTOR_GROUP_DENSE) {
    return 0;
  }
  const UF_long row_begin = solver->row_refactor_group_ptr[group];
  const UF_long row_end = solver->row_refactor_group_ptr[group + 1u];
  if (row_begin >= row_end) {
    return 0;
  }
  const UF_long width = row_end - row_begin;
  if (kls_row_refactor_compact_panel_valid_prefix(solver, group) < width ||
      width > (UF_long)(SIZE_MAX / sizeof(double)) / width) {
    return 0;
  }
  const UF_long dense_entries = width * width;
  const UF_long trailing_len = solver->row_refactor_group_trailing_len[group];
  if (trailing_len > 0u &&
      width > ((UF_long)(SIZE_MAX / sizeof(double)) - dense_entries) /
                trailing_len) {
    return 0;
  }
  const UF_long panel_entries = dense_entries + width * trailing_len;
  const double *dense_panel =
    kls_row_refactor_compact_panel_const(solver, group, panel_entries);
  if (dense_panel == NULL) {
    return 0;
  }
  if (dense_panel_out != NULL) {
    *dense_panel_out = dense_panel;
  }
  if (row_begin_out != NULL) {
    *row_begin_out = row_begin;
  }
  if (row_end_out != NULL) {
    *row_end_out = row_end;
  }
  if (trailing_len_out != NULL) {
    *trailing_len_out = trailing_len;
  }
  return 1;
}

static int kls_try_row_refactor_compact_panel_lower_solve_group(
  kls_solver *solver,
  UF_long block_begin,
  UF_long block_end,
  UF_long row,
  double *work,
  UF_long *next_row_out) {
  if (next_row_out != NULL) {
    *next_row_out = row;
  }
  if (solver == NULL || work == NULL || row >= solver->n ||
      solver->row_refactor_row_group == NULL ||
      solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_l_cols == NULL ||
      solver->row_refactor_l_internal_ptr == NULL) {
    return 0;
  }
  const UF_long group = solver->row_refactor_row_group[row];
  const double *dense_panel = NULL;
  UF_long row_begin = 0;
  UF_long row_end = 0;
  if (!kls_row_refactor_compact_panel_solve_group(
        solver, group, &dense_panel, &row_begin, &row_end, NULL) ||
      row != row_begin || row_begin < block_begin || row_end > block_end) {
    return 0;
  }
  const UF_long width = row_end - row_begin;
  const int shape_proven =
    solver->row_refactor_group_shape_valid != NULL &&
    solver->row_refactor_group_shape_valid[group] != 0u;
  if (!shape_proven) {
    for (UF_long local_row = 0; local_row < width; ++local_row) {
      const UF_long r = row_begin + local_row;
      const UF_long l_begin = solver->row_refactor_l_ptr[r];
      const UF_long l_internal = solver->row_refactor_l_internal_ptr[r];
      const UF_long l_end = solver->row_refactor_l_ptr[r + 1u];
      if (l_begin > l_internal || l_internal > l_end ||
          l_end - l_internal != local_row) {
        return 0;
      }
      for (UF_long offset = 0; offset < local_row; ++offset) {
        if (solver->row_refactor_l_cols[l_internal + offset] !=
            row_begin + offset) {
          return 0;
        }
      }
    }
  }

  UF_long panel_entries = 0;
  for (UF_long local_row = 0; local_row < width; ++local_row) {
    const UF_long r = row_begin + local_row;
    double value = work[r];
    const UF_long l_begin = solver->row_refactor_l_ptr[r];
    const UF_long l_internal = solver->row_refactor_l_internal_ptr[r];
    for (UF_long p = l_begin; p < l_internal; ++p) {
      value -= kls_row_refactor_l_solve_row_value(solver, r, p) *
               work[solver->row_refactor_l_cols[p]];
    }
    const double *row_panel = dense_panel + local_row * width;
    for (UF_long dep = 0; dep < local_row; ++dep) {
      value -= row_panel[dep] * work[row_begin + dep];
    }
    panel_entries = panel_entries > UF_long_max - local_row
      ? UF_long_max : panel_entries + local_row;
    work[r] = value;
  }
  kls_record_row_refactor_compact_panel_group_solve(
    solver, width, panel_entries);
  if (next_row_out != NULL) {
    *next_row_out = row_end;
  }
  return 1;
}

static int kls_try_row_refactor_compact_panel_upper_solve_group(
  kls_solver *solver,
  UF_long block_begin,
  UF_long block_end,
  UF_long row_end,
  double *work,
  const double *udiag,
  UF_long *next_row_end_out) {
  if (next_row_end_out != NULL) {
    *next_row_end_out = row_end;
  }
  if (solver == NULL || work == NULL || udiag == NULL ||
      row_end == 0u || row_end > solver->n ||
      solver->row_refactor_row_group == NULL ||
      solver->row_refactor_u_ptr == NULL ||
      solver->row_refactor_u_cols == NULL ||
      solver->row_refactor_group_trailing_begin == NULL) {
    return 0;
  }
  const UF_long last_row = row_end - 1u;
  const UF_long group = solver->row_refactor_row_group[last_row];
  const double *dense_panel = NULL;
  UF_long row_begin = 0;
  UF_long group_end = 0;
  UF_long trailing_len = 0;
  if (!kls_row_refactor_compact_panel_solve_group(
        solver, group, &dense_panel, &row_begin, &group_end, &trailing_len) ||
      group_end != row_end || row_begin < block_begin || group_end > block_end) {
    return 0;
  }
  const UF_long width = group_end - row_begin;
  const int shape_proven =
    solver->row_refactor_group_shape_valid != NULL &&
    solver->row_refactor_group_shape_valid[group] != 0u;
  const UF_long dense_entries = width * width;
  const double *trailing_panel =
    trailing_len > 0u ? dense_panel + dense_entries : NULL;
  const UF_long trailing_begin =
    solver->row_refactor_group_trailing_begin[group];
  const UF_long row_end_u = solver->row_refactor_u_ptr[group_end];
  if (trailing_len > 0u &&
      (trailing_begin > row_end_u ||
       trailing_len > row_end_u - trailing_begin)) {
    return 0;
  }
  const UF_long *trailing_cols =
    trailing_len > 0u ? solver->row_refactor_u_cols + trailing_begin : NULL;
  if (!shape_proven) {
    for (UF_long local_row = 0; local_row < width; ++local_row) {
      const UF_long r = row_begin + local_row;
      const UF_long u_begin = solver->row_refactor_u_ptr[r];
      const UF_long u_end = solver->row_refactor_u_ptr[r + 1u];
      const UF_long dense_len = width - local_row - 1u;
      if (u_end < u_begin || u_end - u_begin != dense_len + trailing_len) {
        return 0;
      }
      const UF_long *u_cols = solver->row_refactor_u_cols + u_begin;
      for (UF_long offset = 0; offset < dense_len; ++offset) {
        if (u_cols[offset] != r + 1u + offset) {
          return 0;
        }
      }
      for (UF_long offset = 0; offset < trailing_len; ++offset) {
        if (u_cols[dense_len + offset] != trailing_cols[offset]) {
          return 0;
        }
      }
    }
  }

  UF_long panel_entries = 0;
  for (UF_long remaining = width; remaining > 0u; --remaining) {
    const UF_long local_row = remaining - 1u;
    const UF_long r = row_begin + local_row;
    double value = work[r];
    const double *row_panel = dense_panel + local_row * width;
    for (UF_long target = local_row + 1u; target < width; ++target) {
      value -= row_panel[target] * work[row_begin + target];
    }
    const UF_long dense_len = width - local_row - 1u;
    if (trailing_len > 0u) {
      const double *row_trailing = trailing_panel + local_row * trailing_len;
      for (UF_long offset = 0; offset < trailing_len; ++offset) {
        value -= row_trailing[offset] * work[trailing_cols[offset]];
      }
    }
    const UF_long row_entries =
      dense_len > UF_long_max - trailing_len
        ? UF_long_max : dense_len + trailing_len;
    panel_entries = panel_entries > UF_long_max - row_entries
      ? UF_long_max : panel_entries + row_entries;
    work[r] = value / udiag[r];
  }
  kls_record_row_refactor_compact_panel_group_solve(
    solver, width, panel_entries);
  if (next_row_end_out != NULL) {
    *next_row_end_out = row_begin;
  }
  return 1;
}

static int kls_try_row_refactor_compact_panel_lower_solve_group_multi(
  kls_solver *solver,
  UF_long block_begin,
  UF_long block_end,
  UF_long row,
  double *work,
  UF_long rhs_count,
  UF_long *next_row_out) {
  if (next_row_out != NULL) {
    *next_row_out = row;
  }
  if (solver == NULL || work == NULL || rhs_count == 0u ||
      row >= solver->n ||
      solver->row_refactor_row_group == NULL ||
      solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_l_cols == NULL ||
      solver->row_refactor_l_internal_ptr == NULL) {
    return 0;
  }
  const UF_long group = solver->row_refactor_row_group[row];
  const double *dense_panel = NULL;
  UF_long row_begin = 0;
  UF_long row_end = 0;
  if (!kls_row_refactor_compact_panel_solve_group(
        solver, group, &dense_panel, &row_begin, &row_end, NULL) ||
      row != row_begin || row_begin < block_begin || row_end > block_end) {
    return 0;
  }
  const UF_long width = row_end - row_begin;
  const int shape_proven =
    solver->row_refactor_group_shape_valid != NULL &&
    solver->row_refactor_group_shape_valid[group] != 0u;
  if (!shape_proven) {
    for (UF_long local_row = 0; local_row < width; ++local_row) {
      const UF_long r = row_begin + local_row;
      const UF_long l_begin = solver->row_refactor_l_ptr[r];
      const UF_long l_internal = solver->row_refactor_l_internal_ptr[r];
      const UF_long l_end = solver->row_refactor_l_ptr[r + 1u];
      if (l_begin > l_internal || l_internal > l_end ||
          l_end - l_internal != local_row) {
        return 0;
      }
      for (UF_long offset = 0; offset < local_row; ++offset) {
        if (solver->row_refactor_l_cols[l_internal + offset] !=
            row_begin + offset) {
          return 0;
        }
      }
    }
  }

  UF_long panel_entries = 0;
  for (UF_long local_row = 0; local_row < width; ++local_row) {
    const UF_long r = row_begin + local_row;
    double *wrow = work + r * rhs_count;
    const UF_long l_begin = solver->row_refactor_l_ptr[r];
    const UF_long l_internal = solver->row_refactor_l_internal_ptr[r];
    for (UF_long p = l_begin; p < l_internal; ++p) {
      const double lij = kls_row_refactor_l_solve_row_value(solver, r, p);
      const double *wcol =
        work + solver->row_refactor_l_cols[p] * rhs_count;
      for (UF_long rhs_offset = 0; rhs_offset < rhs_count; ++rhs_offset) {
        wrow[rhs_offset] -= lij * wcol[rhs_offset];
      }
    }
    const double *row_panel = dense_panel + local_row * width;
    for (UF_long dep = 0; dep < local_row; ++dep) {
      const double lij = row_panel[dep];
      const double *wcol = work + (row_begin + dep) * rhs_count;
      for (UF_long rhs_offset = 0; rhs_offset < rhs_count; ++rhs_offset) {
        wrow[rhs_offset] -= lij * wcol[rhs_offset];
      }
    }
    panel_entries = panel_entries > UF_long_max - local_row
      ? UF_long_max : panel_entries + local_row;
  }
  kls_record_row_refactor_compact_panel_group_solve(
    solver, width, panel_entries);
  if (next_row_out != NULL) {
    *next_row_out = row_end;
  }
  return 1;
}

static int kls_try_row_refactor_compact_panel_upper_solve_group_multi(
  kls_solver *solver,
  UF_long block_begin,
  UF_long block_end,
  UF_long row_end,
  double *work,
  UF_long rhs_count,
  const double *udiag,
  UF_long *next_row_end_out) {
  if (next_row_end_out != NULL) {
    *next_row_end_out = row_end;
  }
  if (solver == NULL || work == NULL || rhs_count == 0u || udiag == NULL ||
      row_end == 0u || row_end > solver->n ||
      solver->row_refactor_row_group == NULL ||
      solver->row_refactor_u_ptr == NULL ||
      solver->row_refactor_u_cols == NULL ||
      solver->row_refactor_group_trailing_begin == NULL) {
    return 0;
  }
  const UF_long last_row = row_end - 1u;
  const UF_long group = solver->row_refactor_row_group[last_row];
  const double *dense_panel = NULL;
  UF_long row_begin = 0;
  UF_long group_end = 0;
  UF_long trailing_len = 0;
  if (!kls_row_refactor_compact_panel_solve_group(
        solver, group, &dense_panel, &row_begin, &group_end, &trailing_len) ||
      group_end != row_end || row_begin < block_begin || group_end > block_end) {
    return 0;
  }
  const UF_long width = group_end - row_begin;
  const int shape_proven =
    solver->row_refactor_group_shape_valid != NULL &&
    solver->row_refactor_group_shape_valid[group] != 0u;
  const UF_long dense_entries = width * width;
  const double *trailing_panel =
    trailing_len > 0u ? dense_panel + dense_entries : NULL;
  const UF_long trailing_begin =
    solver->row_refactor_group_trailing_begin[group];
  const UF_long row_end_u = solver->row_refactor_u_ptr[group_end];
  if (trailing_len > 0u &&
      (trailing_begin > row_end_u ||
       trailing_len > row_end_u - trailing_begin)) {
    return 0;
  }
  const UF_long *trailing_cols =
    trailing_len > 0u ? solver->row_refactor_u_cols + trailing_begin : NULL;
  if (!shape_proven) {
    for (UF_long local_row = 0; local_row < width; ++local_row) {
      const UF_long r = row_begin + local_row;
      const UF_long u_begin = solver->row_refactor_u_ptr[r];
      const UF_long u_end = solver->row_refactor_u_ptr[r + 1u];
      const UF_long dense_len = width - local_row - 1u;
      if (u_end < u_begin || u_end - u_begin != dense_len + trailing_len) {
        return 0;
      }
      const UF_long *u_cols = solver->row_refactor_u_cols + u_begin;
      for (UF_long offset = 0; offset < dense_len; ++offset) {
        if (u_cols[offset] != r + 1u + offset) {
          return 0;
        }
      }
      for (UF_long offset = 0; offset < trailing_len; ++offset) {
        if (u_cols[dense_len + offset] != trailing_cols[offset]) {
          return 0;
        }
      }
    }
  }

  UF_long panel_entries = 0;
  for (UF_long remaining = width; remaining > 0u; --remaining) {
    const UF_long local_row = remaining - 1u;
    const UF_long r = row_begin + local_row;
    double *wrow = work + r * rhs_count;
    const double *row_panel = dense_panel + local_row * width;
    for (UF_long target = local_row + 1u; target < width; ++target) {
      const double uij = row_panel[target];
      const double *wcol = work + (row_begin + target) * rhs_count;
      for (UF_long rhs_offset = 0; rhs_offset < rhs_count; ++rhs_offset) {
        wrow[rhs_offset] -= uij * wcol[rhs_offset];
      }
    }
    const UF_long dense_len = width - local_row - 1u;
    if (trailing_len > 0u) {
      const double *row_trailing = trailing_panel + local_row * trailing_len;
      for (UF_long offset = 0; offset < trailing_len; ++offset) {
        const double uij = row_trailing[offset];
        const double *wcol = work + trailing_cols[offset] * rhs_count;
        for (UF_long rhs_offset = 0; rhs_offset < rhs_count; ++rhs_offset) {
          wrow[rhs_offset] -= uij * wcol[rhs_offset];
        }
      }
    }
    for (UF_long rhs_offset = 0; rhs_offset < rhs_count; ++rhs_offset) {
      wrow[rhs_offset] /= udiag[r];
    }
    const UF_long row_entries =
      dense_len > UF_long_max - trailing_len
        ? UF_long_max : dense_len + trailing_len;
    panel_entries = panel_entries > UF_long_max - row_entries
      ? UF_long_max : panel_entries + row_entries;
  }
  kls_record_row_refactor_compact_panel_group_solve(
    solver, width, panel_entries);
  if (next_row_end_out != NULL) {
    *next_row_end_out = row_begin;
  }
  return 1;
}

static int kls_try_row_refactor_compact_panel_u_transpose_solve_group(
  kls_solver *solver,
  UF_long block_begin,
  UF_long block_end,
  UF_long row,
  double *work,
  const double *udiag,
  UF_long *next_row_out) {
  if (next_row_out != NULL) {
    *next_row_out = row;
  }
  if (solver == NULL || work == NULL || udiag == NULL ||
      row >= solver->n ||
      solver->row_refactor_row_group == NULL ||
      solver->row_refactor_u_ptr == NULL ||
      solver->row_refactor_u_cols == NULL ||
      solver->row_refactor_group_trailing_begin == NULL) {
    return 0;
  }
  const UF_long group = solver->row_refactor_row_group[row];
  const double *dense_panel = NULL;
  UF_long row_begin = 0;
  UF_long row_end = 0;
  UF_long trailing_len = 0;
  if (!kls_row_refactor_compact_panel_solve_group(
        solver, group, &dense_panel, &row_begin, &row_end, &trailing_len) ||
      row != row_begin || row_begin < block_begin || row_end > block_end) {
    return 0;
  }
  const UF_long width = row_end - row_begin;
  const UF_long dense_entries = width * width;
  const double *trailing_panel =
    trailing_len > 0u ? dense_panel + dense_entries : NULL;
  const UF_long trailing_begin =
    solver->row_refactor_group_trailing_begin[group];
  const UF_long row_end_u = solver->row_refactor_u_ptr[row_end];
  if (trailing_len > 0u &&
      (trailing_begin > row_end_u ||
       trailing_len > row_end_u - trailing_begin)) {
    return 0;
  }
  const UF_long *trailing_cols =
    trailing_len > 0u ? solver->row_refactor_u_cols + trailing_begin : NULL;
  for (UF_long offset = 0; offset < trailing_len; ++offset) {
    if (trailing_cols[offset] < row_end ||
        trailing_cols[offset] >= block_end) {
      return 0;
    }
  }
  for (UF_long local_row = 0; local_row < width; ++local_row) {
    const UF_long r = row_begin + local_row;
    const UF_long u_begin = solver->row_refactor_u_ptr[r];
    const UF_long u_end = solver->row_refactor_u_ptr[r + 1u];
    const UF_long dense_len = width - local_row - 1u;
    if (u_end < u_begin || u_end - u_begin != dense_len + trailing_len) {
      return 0;
    }
    const UF_long *u_cols = solver->row_refactor_u_cols + u_begin;
    for (UF_long offset = 0; offset < dense_len; ++offset) {
      if (u_cols[offset] != r + 1u + offset) {
        return 0;
      }
    }
    for (UF_long offset = 0; offset < trailing_len; ++offset) {
      if (u_cols[dense_len + offset] != trailing_cols[offset]) {
        return 0;
      }
    }
  }

  UF_long panel_entries = 0;
  for (UF_long local_row = 0; local_row < width; ++local_row) {
    const UF_long r = row_begin + local_row;
    double value = work[r];
    for (UF_long dep = 0; dep < local_row; ++dep) {
      const double *dep_panel = dense_panel + dep * width;
      value -= dep_panel[local_row] * work[row_begin + dep];
    }
    value /= udiag[r];
    work[r] = value;
    if (trailing_len > 0u) {
      const double *row_trailing = trailing_panel + local_row * trailing_len;
      for (UF_long offset = 0; offset < trailing_len; ++offset) {
        work[trailing_cols[offset]] -= row_trailing[offset] * value;
      }
    }
    const UF_long row_entries =
      local_row > UF_long_max - trailing_len
        ? UF_long_max : local_row + trailing_len;
    panel_entries = panel_entries > UF_long_max - row_entries
      ? UF_long_max : panel_entries + row_entries;
  }
  kls_record_row_refactor_compact_panel_group_solve(
    solver, width, panel_entries);
  if (next_row_out != NULL) {
    *next_row_out = row_end;
  }
  return 1;
}

static int kls_try_row_refactor_compact_panel_l_transpose_solve_group(
  kls_solver *solver,
  UF_long block_begin,
  UF_long block_end,
  UF_long row_end,
  double *work,
  UF_long *next_row_end_out) {
  if (next_row_end_out != NULL) {
    *next_row_end_out = row_end;
  }
  if (solver == NULL || work == NULL ||
      row_end == 0u || row_end > solver->n ||
      solver->row_refactor_row_group == NULL ||
      solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_l_cols == NULL ||
      solver->row_refactor_l_internal_ptr == NULL) {
    return 0;
  }
  const UF_long last_row = row_end - 1u;
  const UF_long group = solver->row_refactor_row_group[last_row];
  const double *dense_panel = NULL;
  UF_long row_begin = 0;
  UF_long group_end = 0;
  if (!kls_row_refactor_compact_panel_solve_group(
        solver, group, &dense_panel, &row_begin, &group_end, NULL) ||
      group_end != row_end || row_begin < block_begin || group_end > block_end) {
    return 0;
  }
  const UF_long width = group_end - row_begin;
  for (UF_long local_row = 0; local_row < width; ++local_row) {
    const UF_long r = row_begin + local_row;
    const UF_long l_begin = solver->row_refactor_l_ptr[r];
    const UF_long l_internal = solver->row_refactor_l_internal_ptr[r];
    const UF_long l_end = solver->row_refactor_l_ptr[r + 1u];
    if (l_begin > l_internal || l_internal > l_end ||
        l_end - l_internal != local_row) {
      return 0;
    }
    for (UF_long p = l_begin; p < l_internal; ++p) {
      const UF_long col = solver->row_refactor_l_cols[p];
      if (col < block_begin || col >= row_begin) {
        return 0;
      }
    }
    for (UF_long offset = 0; offset < local_row; ++offset) {
      if (solver->row_refactor_l_cols[l_internal + offset] !=
          row_begin + offset) {
        return 0;
      }
    }
  }

  UF_long panel_entries = 0;
  for (UF_long remaining = width; remaining > 0u; --remaining) {
    const UF_long local_row = remaining - 1u;
    const UF_long r = row_begin + local_row;
    double value = work[r];
    for (UF_long target = local_row + 1u; target < width; ++target) {
      const double *target_panel = dense_panel + target * width;
      value -= target_panel[local_row] * work[row_begin + target];
    }
    work[r] = value;
    const UF_long l_begin = solver->row_refactor_l_ptr[r];
    const UF_long l_internal = solver->row_refactor_l_internal_ptr[r];
    for (UF_long p = l_begin; p < l_internal; ++p) {
      work[solver->row_refactor_l_cols[p]] -=
        kls_row_refactor_l_solve_row_value(solver, r, p) * value;
    }
    const UF_long row_entries = width - local_row - 1u;
    panel_entries = panel_entries > UF_long_max - row_entries
      ? UF_long_max : panel_entries + row_entries;
  }
  kls_record_row_refactor_compact_panel_group_solve(
    solver, width, panel_entries);
  if (next_row_end_out != NULL) {
    *next_row_end_out = row_begin;
  }
  return 1;
}

static int kls_try_row_refactor_compact_panel_u_transpose_solve_group_multi(
  kls_solver *solver,
  UF_long block_begin,
  UF_long block_end,
  UF_long row,
  double *work,
  UF_long rhs_count,
  const double *udiag,
  UF_long *next_row_out) {
  if (next_row_out != NULL) {
    *next_row_out = row;
  }
  if (solver == NULL || work == NULL || rhs_count == 0u ||
      udiag == NULL || row >= solver->n ||
      solver->row_refactor_row_group == NULL ||
      solver->row_refactor_u_ptr == NULL ||
      solver->row_refactor_u_cols == NULL ||
      solver->row_refactor_group_trailing_begin == NULL) {
    return 0;
  }
  const UF_long group = solver->row_refactor_row_group[row];
  const double *dense_panel = NULL;
  UF_long row_begin = 0;
  UF_long row_end = 0;
  UF_long trailing_len = 0;
  if (!kls_row_refactor_compact_panel_solve_group(
        solver, group, &dense_panel, &row_begin, &row_end, &trailing_len) ||
      row != row_begin || row_begin < block_begin || row_end > block_end) {
    return 0;
  }
  const UF_long width = row_end - row_begin;
  const UF_long dense_entries = width * width;
  const double *trailing_panel =
    trailing_len > 0u ? dense_panel + dense_entries : NULL;
  const UF_long trailing_begin =
    solver->row_refactor_group_trailing_begin[group];
  const UF_long row_end_u = solver->row_refactor_u_ptr[row_end];
  if (trailing_len > 0u &&
      (trailing_begin > row_end_u ||
       trailing_len > row_end_u - trailing_begin)) {
    return 0;
  }
  const UF_long *trailing_cols =
    trailing_len > 0u ? solver->row_refactor_u_cols + trailing_begin : NULL;
  for (UF_long offset = 0; offset < trailing_len; ++offset) {
    if (trailing_cols[offset] < row_end ||
        trailing_cols[offset] >= block_end) {
      return 0;
    }
  }
  for (UF_long local_row = 0; local_row < width; ++local_row) {
    const UF_long r = row_begin + local_row;
    const UF_long u_begin = solver->row_refactor_u_ptr[r];
    const UF_long u_end = solver->row_refactor_u_ptr[r + 1u];
    const UF_long dense_len = width - local_row - 1u;
    if (u_end < u_begin || u_end - u_begin != dense_len + trailing_len) {
      return 0;
    }
    const UF_long *u_cols = solver->row_refactor_u_cols + u_begin;
    for (UF_long offset = 0; offset < dense_len; ++offset) {
      if (u_cols[offset] != r + 1u + offset) {
        return 0;
      }
    }
    for (UF_long offset = 0; offset < trailing_len; ++offset) {
      if (u_cols[dense_len + offset] != trailing_cols[offset]) {
        return 0;
      }
    }
  }

  UF_long panel_entries = 0;
  for (UF_long local_row = 0; local_row < width; ++local_row) {
    const UF_long r = row_begin + local_row;
    double *wrow = work + r * rhs_count;
    for (UF_long dep = 0; dep < local_row; ++dep) {
      const double uij = dense_panel[dep * width + local_row];
      const double *wdep = work + (row_begin + dep) * rhs_count;
      for (UF_long rhs_offset = 0; rhs_offset < rhs_count; ++rhs_offset) {
        wrow[rhs_offset] -= uij * wdep[rhs_offset];
      }
    }
    for (UF_long rhs_offset = 0; rhs_offset < rhs_count; ++rhs_offset) {
      wrow[rhs_offset] /= udiag[r];
    }
    if (trailing_len > 0u) {
      const double *row_trailing = trailing_panel + local_row * trailing_len;
      for (UF_long offset = 0; offset < trailing_len; ++offset) {
        const double uij = row_trailing[offset];
        double *wcol = work + trailing_cols[offset] * rhs_count;
        for (UF_long rhs_offset = 0; rhs_offset < rhs_count; ++rhs_offset) {
          wcol[rhs_offset] -= uij * wrow[rhs_offset];
        }
      }
    }
    const UF_long row_entries =
      local_row > UF_long_max - trailing_len
        ? UF_long_max : local_row + trailing_len;
    panel_entries = panel_entries > UF_long_max - row_entries
      ? UF_long_max : panel_entries + row_entries;
  }
  kls_record_row_refactor_compact_panel_group_solve(
    solver, width, panel_entries);
  if (next_row_out != NULL) {
    *next_row_out = row_end;
  }
  return 1;
}

static int kls_try_row_refactor_compact_panel_l_transpose_solve_group_multi(
  kls_solver *solver,
  UF_long block_begin,
  UF_long block_end,
  UF_long row_end,
  double *work,
  UF_long rhs_count,
  UF_long *next_row_end_out) {
  if (next_row_end_out != NULL) {
    *next_row_end_out = row_end;
  }
  if (solver == NULL || work == NULL || rhs_count == 0u ||
      row_end == 0u || row_end > solver->n ||
      solver->row_refactor_row_group == NULL ||
      solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_l_cols == NULL ||
      solver->row_refactor_l_internal_ptr == NULL) {
    return 0;
  }
  const UF_long last_row = row_end - 1u;
  const UF_long group = solver->row_refactor_row_group[last_row];
  const double *dense_panel = NULL;
  UF_long row_begin = 0;
  UF_long group_end = 0;
  if (!kls_row_refactor_compact_panel_solve_group(
        solver, group, &dense_panel, &row_begin, &group_end, NULL) ||
      group_end != row_end || row_begin < block_begin || group_end > block_end) {
    return 0;
  }
  const UF_long width = group_end - row_begin;
  for (UF_long local_row = 0; local_row < width; ++local_row) {
    const UF_long r = row_begin + local_row;
    const UF_long l_begin = solver->row_refactor_l_ptr[r];
    const UF_long l_internal = solver->row_refactor_l_internal_ptr[r];
    const UF_long l_end = solver->row_refactor_l_ptr[r + 1u];
    if (l_begin > l_internal || l_internal > l_end ||
        l_end - l_internal != local_row) {
      return 0;
    }
    for (UF_long p = l_begin; p < l_internal; ++p) {
      const UF_long col = solver->row_refactor_l_cols[p];
      if (col < block_begin || col >= row_begin) {
        return 0;
      }
    }
    for (UF_long offset = 0; offset < local_row; ++offset) {
      if (solver->row_refactor_l_cols[l_internal + offset] !=
          row_begin + offset) {
        return 0;
      }
    }
  }

  UF_long panel_entries = 0;
  for (UF_long remaining = width; remaining > 0u; --remaining) {
    const UF_long local_row = remaining - 1u;
    const UF_long r = row_begin + local_row;
    double *wrow = work + r * rhs_count;
    for (UF_long target = local_row + 1u; target < width; ++target) {
      const double lij = dense_panel[target * width + local_row];
      const double *wtarget = work + (row_begin + target) * rhs_count;
      for (UF_long rhs_offset = 0; rhs_offset < rhs_count; ++rhs_offset) {
        wrow[rhs_offset] -= lij * wtarget[rhs_offset];
      }
    }
    const UF_long l_begin = solver->row_refactor_l_ptr[r];
    const UF_long l_internal = solver->row_refactor_l_internal_ptr[r];
    for (UF_long p = l_begin; p < l_internal; ++p) {
      const double lij = kls_row_refactor_l_solve_row_value(solver, r, p);
      double *wcol = work + solver->row_refactor_l_cols[p] * rhs_count;
      for (UF_long rhs_offset = 0; rhs_offset < rhs_count; ++rhs_offset) {
        wcol[rhs_offset] -= lij * wrow[rhs_offset];
      }
    }
    const UF_long row_entries = width - local_row - 1u;
    panel_entries = panel_entries > UF_long_max - row_entries
      ? UF_long_max : panel_entries + row_entries;
  }
  kls_record_row_refactor_compact_panel_group_solve(
    solver, width, panel_entries);
  if (next_row_end_out != NULL) {
    *next_row_end_out = row_begin;
  }
  return 1;
}

static UF_long kls_row_refactor_compact_panel_solve_value_count(
  const kls_solver *solver) {
  if (solver == NULL ||
      kls_row_refactor_solve_uses_direct_values(solver) ||
      !solver->row_refactor_values_ready ||
      solver->row_refactor_group_ptr == NULL ||
      solver->row_refactor_group_kind == NULL ||
      solver->row_refactor_group_trailing_len == NULL ||
      solver->row_refactor_compact_panel_valid == NULL) {
    return 0;
  }
  UF_long values = 0;
  for (UF_long group = 0; group < solver->row_refactor_group_count; ++group) {
    if (solver->row_refactor_group_kind[group] !=
          KLS_ROW_REFACTOR_GROUP_DENSE) {
      continue;
    }
    const UF_long begin = solver->row_refactor_group_ptr[group];
    const UF_long end = solver->row_refactor_group_ptr[group + 1u];
    if (begin >= end) {
      continue;
    }
    const UF_long width = end - begin;
    if (kls_row_refactor_compact_panel_valid_prefix(solver, group) < width) {
      continue;
    }
    const UF_long trailing_len = solver->row_refactor_group_trailing_len[group];
    if (width > UF_long_max / width) {
      values = UF_long_max;
      continue;
    }
    const UF_long dense_entries = width * width;
    const UF_long trailing_values =
      trailing_len > 0u && width > UF_long_max / trailing_len
        ? UF_long_max
        : width * trailing_len;
    if (dense_entries > UF_long_max - trailing_values) {
      values = UF_long_max;
      continue;
    }
    const UF_long panel_entries = dense_entries + trailing_values;
    if (kls_row_refactor_compact_panel_const(solver, group, panel_entries) ==
        NULL) {
      continue;
    }
    const UF_long dense_values = dense_entries - width;
    const UF_long panel_values =
      dense_values > UF_long_max - trailing_values
        ? UF_long_max : dense_values + trailing_values;
    values = values > UF_long_max - panel_values
      ? UF_long_max : values + panel_values;
  }
  return values;
}

static int kls_row_refactor_solve_is_eligible(const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      (!solver->row_refactor_values_ready &&
       !kls_row_refactor_solve_uses_direct_values(solver)) ||
      solver->symbolic->nblocks == 0u ||
      solver->symbolic->R == NULL ||
      solver->row_refactor_pattern_n != solver->n ||
      solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_u_ptr == NULL ||
      solver->numeric->Pnum == NULL || solver->symbolic->Q == NULL ||
      solver->numeric->Udiag == NULL || solver->numeric->Xwork == NULL) {
    return 0;
  }
  if (solver->symbolic->R[0] != 0 ||
      solver->symbolic->R[solver->symbolic->nblocks] != solver->n) {
    return 0;
  }
  if (solver->symbolic->nblocks > 1u &&
      (solver->numeric->Offp == NULL || solver->numeric->Offi == NULL ||
       solver->numeric->Offx == NULL)) {
    return 0;
  }
  if ((solver->common.scale > 0 && solver->numeric->Rs == NULL) ||
      (solver->common.scale <= 0 && solver->numeric->Rs != NULL)) {
    return 0;
  }
  if (solver->row_refactor_l_ptr[solver->n] > 0u &&
      (solver->row_refactor_l_cols == NULL ||
       !kls_row_refactor_solve_has_values(solver, 0))) {
    return 0;
  }
  if (solver->row_refactor_u_ptr[solver->n] > 0u &&
      (solver->row_refactor_u_cols == NULL ||
       !kls_row_refactor_solve_has_values(solver, 1))) {
    return 0;
  }
  return 1;
}

static int kls_validate_row_refactor_solve_metadata(kls_solver *solver) {
  if (solver == NULL) {
    return 0;
  }
  if (solver->row_refactor_solve_validated) {
    return 1;
  }
  const UF_long n = solver->n;
  const UF_long nblocks = solver->symbolic->nblocks;
  const UF_long *r = solver->symbolic->R;
  const UF_long *pnum = solver->numeric->Pnum;
  const UF_long *q = solver->symbolic->Q;
  const double *udiag = (const double *)solver->numeric->Udiag;
  const double *rs = solver->numeric->Rs;
  const UF_long *offp = solver->numeric->Offp;
  const int direct_values =
    kls_row_refactor_solve_uses_direct_values(solver);

  for (UF_long k = 0; k < n; ++k) {
    if (pnum[k] < 0 || pnum[k] >= n || q[k] < 0 || q[k] >= n) {
      return 0;
    }
    if (udiag[k] == 0.0 || !isfinite(udiag[k])) {
      return 0;
    }
    if (rs != NULL && (rs[k] == 0.0 || !isfinite(rs[k]))) {
      return 0;
    }
    if (nblocks > 1u && (offp[k] < 0 || offp[k] > offp[k + 1u])) {
      return 0;
    }
  }
  for (UF_long block = 0; block < nblocks; ++block) {
    if (r[block] < 0 || r[block] > r[block + 1u] ||
        r[block + 1u] > n) {
      return 0;
    }
  }
  if (direct_values) {
    for (UF_long p = 0; p < solver->row_refactor_l_ptr[n]; ++p) {
      if (solver->row_refactor_l_values == NULL ||
          solver->row_refactor_l_values[p] == NULL) {
        return 0;
      }
    }
    for (UF_long p = 0; p < solver->row_refactor_u_ptr[n]; ++p) {
      if (solver->row_refactor_u_values == NULL ||
          solver->row_refactor_u_values[p] == NULL) {
        return 0;
      }
    }
  }
  solver->row_refactor_solve_validated = 1;
  return 1;
}

/* Solve one normal BTF diagonal block and scatter its off-diagonal update
   into the earlier blocks.  Keeping this operation separate lets a matrix
   with one overwhelmingly dominant block run that block through the
   parallel row-solve engine while retaining the exact serial BTF ordering
   for its small fringe. */
static int kls_row_solve_serial_normal_block(kls_solver *solver,
                                             UF_long block,
                                             double *work) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      work == NULL || block >= solver->symbolic->nblocks) {
    return 0;
  }
  const UF_long *r = solver->symbolic->R;
  const UF_long k1 = r[block];
  const UF_long k2 = r[block + 1u];
  const double *udiag = (const double *)solver->numeric->Udiag;

  for (UF_long row = k1; row < k2;) {
    UF_long next_row = row;
    if (kls_try_row_refactor_compact_panel_lower_solve_group(
          solver, k1, k2, row, work, &next_row) &&
        next_row > row) {
      row = next_row;
      continue;
    }
    double value = work[row];
    const UF_long begin = solver->row_refactor_l_ptr[row];
    const UF_long end = solver->row_refactor_l_ptr[row + 1u];
    for (UF_long p = begin; p < end; ++p) {
      value -= kls_row_refactor_l_solve_row_value(solver, row, p) *
               work[solver->row_refactor_l_cols[p]];
    }
    work[row] = value;
    row++;
  }

  for (UF_long remaining = k2; remaining > k1;) {
    UF_long next_remaining = remaining;
    if (kls_try_row_refactor_compact_panel_upper_solve_group(
          solver, k1, k2, remaining, work, udiag, &next_remaining) &&
        next_remaining < remaining) {
      remaining = next_remaining;
      continue;
    }
    const UF_long row = remaining - 1u;
    double value = work[row];
    const UF_long begin = solver->row_refactor_u_ptr[row];
    const UF_long end = solver->row_refactor_u_ptr[row + 1u];
    for (UF_long p = begin; p < end; ++p) {
      value -= kls_row_refactor_u_solve_row_value(solver, row, p) *
               work[solver->row_refactor_u_cols[p]];
    }
    work[row] = value / udiag[row];
    remaining--;
  }

  if (block > 0u) {
    const UF_long *offp = solver->numeric->Offp;
    const UF_long *offi = solver->numeric->Offi;
    const double *offx = (const double *)solver->numeric->Offx;
    if (offp == NULL || offi == NULL || offx == NULL) {
      return 0;
    }
    for (UF_long k = k1; k < k2; ++k) {
      const double wk = work[k];
      for (UF_long p = offp[k]; p < offp[k + 1u]; ++p) {
        const UF_long row = offi[p];
        if (row >= k1) {
          return 0;
        }
        work[row] -= offx[p] * wk;
      }
    }
  }
  return 1;
}


static int kls_try_row_refactor_solve_one_rhs(kls_solver *solver,
                                              int kernel_transpose,
                                              double *x) {
  const UF_long n = solver->n;
  const UF_long nblocks = solver->symbolic->nblocks;
  const UF_long *r = solver->symbolic->R;
  const UF_long *pnum = solver->numeric->Pnum;
  const UF_long *q = solver->symbolic->Q;
  const double *udiag = (const double *)solver->numeric->Udiag;
  const double *rs = solver->numeric->Rs;
  const UF_long *offp = solver->numeric->Offp;
  const UF_long *offi = solver->numeric->Offi;
  const double *offx = (const double *)solver->numeric->Offx;
  double *work = (double *)solver->numeric->Xwork;

  if (!kernel_transpose) {
    if (kls_try_parallel_row_solve_one_rhs(solver, x)) {
      return 1;
    }

    for (UF_long k = 0; k < n; ++k) {
      const double rhs_value = x[pnum[k]];
      work[k] = rs != NULL ? rhs_value / rs[k] : rhs_value;
    }

    for (UF_long remaining_block = nblocks; remaining_block > 0u;
         --remaining_block) {
      if (!kls_row_solve_serial_normal_block(
            solver, remaining_block - 1u, work)) {
        return 0;
      }
    }

    for (UF_long k = 0; k < n; ++k) {
      x[q[k]] = work[k];
    }
  } else {
    if (kls_try_parallel_row_solve_transpose_one_rhs(solver, x)) {
      return 1;
    }

    for (UF_long k = 0; k < n; ++k) {
      work[k] = x[q[k]];
    }

    for (UF_long block = 0; block < nblocks; ++block) {
      const UF_long k1 = r[block];
      const UF_long k2 = r[block + 1u];
      if (block > 0u) {
        for (UF_long k = k1; k < k2; ++k) {
          double value = work[k];
          const UF_long begin = offp[k];
          const UF_long end = offp[k + 1u];
          for (UF_long p = begin; p < end; ++p) {
            const UF_long row = offi[p];
            if (row >= k1) {
              return 0;
            }
            value -= offx[p] * work[row];
          }
          work[k] = value;
        }
      }

      for (UF_long row = k1; row < k2;) {
        UF_long next_row = row;
        if (kls_try_row_refactor_compact_panel_u_transpose_solve_group(
              solver, k1, k2, row, work, udiag, &next_row) &&
            next_row > row) {
          row = next_row;
          continue;
        }
        const double value = work[row] / udiag[row];
        work[row] = value;
        const UF_long begin = solver->row_refactor_u_ptr[row];
        const UF_long end = solver->row_refactor_u_ptr[row + 1u];
        for (UF_long p = begin; p < end; ++p) {
          work[solver->row_refactor_u_cols[p]] -=
            kls_row_refactor_u_solve_row_value(solver, row, p) * value;
        }
        row++;
      }

      for (UF_long remaining = k2; remaining > k1;) {
        UF_long next_remaining = remaining;
        if (kls_try_row_refactor_compact_panel_l_transpose_solve_group(
              solver, k1, k2, remaining, work, &next_remaining) &&
            next_remaining < remaining) {
          remaining = next_remaining;
          continue;
        }
        const UF_long row = remaining - 1u;
        const double value = work[row];
        const UF_long begin = solver->row_refactor_l_ptr[row];
        const UF_long end = solver->row_refactor_l_ptr[row + 1u];
        for (UF_long p = begin; p < end; ++p) {
          work[solver->row_refactor_l_cols[p]] -=
            kls_row_refactor_l_solve_row_value(solver, row, p) * value;
        }
        remaining--;
      }
    }

    for (UF_long k = 0; k < n; ++k) {
      const double result = work[k];
      x[pnum[k]] = rs != NULL ? result / rs[k] : result;
    }
  }

  return 1;
}

/* A lean row refactor can deliberately leave KLU's column value streams
   stale while its flat row factor is current.  For a transposed kernel, the
   immutable transpose plans turn that row factor back into a dot-product
   triangular walk without first scattering every value through the column
   pointer mirrors.  BTF off-diagonal updates retain the ordinary block
   order; only the two diagonal solves consume the transposed row graph. */
static int kls_try_dirty_row_transpose_plan_solve_one_rhs(
  kls_solver *solver,
  double *x) {
  if (solver == NULL || x == NULL || solver->symbolic == NULL ||
      solver->numeric == NULL ||
      !solver->lean_compact_match_row_factor_active ||
      !solver->row_refactor_values_dirty ||
      solver->row_refactor_values_ready ||
      solver->row_refactor_pattern_n != solver->n ||
      solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_u_ptr == NULL ||
      solver->row_refactor_l_row_values == NULL ||
      solver->row_refactor_u_row_values == NULL ||
      solver->symbolic->R == NULL || solver->symbolic->Q == NULL ||
      solver->numeric->Pnum == NULL || solver->numeric->Udiag == NULL ||
      solver->numeric->Xwork == NULL ||
      solver->options.expected_solves < 16 ||
      solver->n > (UF_long)UINT16_MAX + 1u ||
      solver->row_refactor_l_ptr[solver->n] > (UF_long)UINT16_MAX ||
      solver->row_refactor_u_ptr[solver->n] > (UF_long)UINT16_MAX) {
    return 0;
  }

  if (solver->row_solve_ut_plan.ptr == NULL ||
      solver->row_solve_lt_plan.ptr == NULL) {
    (void)kls_build_one_row_solve_transpose_plan_impl(
      solver, 1, &solver->row_solve_ut_plan, 0);
    (void)kls_build_one_row_solve_transpose_plan_impl(
      solver, 0, &solver->row_solve_lt_plan, 0);
  }
  const kls_row_solve_transpose_plan *ut = &solver->row_solve_ut_plan;
  const kls_row_solve_transpose_plan *lt = &solver->row_solve_lt_plan;
  if (ut->ptr == NULL || lt->ptr == NULL ||
      (ut->ptr[solver->n] > 0u &&
       ((ut->cols == NULL || ut->source_pos == NULL) &&
        (ut->cols16 == NULL || ut->source_pos16 == NULL))) ||
      (lt->ptr[solver->n] > 0u &&
       ((lt->cols == NULL || lt->source_pos == NULL) &&
        (lt->cols16 == NULL || lt->source_pos16 == NULL)))) {
    return 0;
  }

  const UF_long n = solver->n;
  const UF_long nblocks = solver->symbolic->nblocks;
  const UF_long *restrict r = solver->symbolic->R;
  const UF_long *restrict q = solver->symbolic->Q;
  const UF_long *restrict pnum = solver->numeric->Pnum;
  const UF_long *restrict offp = solver->numeric->Offp;
  const UF_long *restrict offi = solver->numeric->Offi;
  const double *restrict offx = (const double *)solver->numeric->Offx;
  const double *restrict udiag =
    (const double *)solver->numeric->Udiag;
  const double *restrict rs = solver->numeric->Rs;
  const double *restrict u_values = solver->row_refactor_u_row_values;
  const double *restrict l_values = solver->row_refactor_l_row_values;
  double *restrict work = (double *)solver->numeric->Xwork;
  if (nblocks == 0u || r[0] != 0u || r[nblocks] != n ||
      (nblocks > 1u && (offp == NULL || offi == NULL || offx == NULL)) ||
      (solver->common.scale > 0 && rs == NULL) ||
      (solver->common.scale <= 0 && rs != NULL)) {
    return 0;
  }

  for (UF_long k = 0u; k < n; ++k) {
    work[k] = x[q[k]];
  }
  for (UF_long block = 0u; block < nblocks; ++block) {
    const UF_long k1 = r[block];
    const UF_long k2 = r[block + 1u];
    if (k1 > k2 || k2 > n) {
      return 0;
    }
    if (block > 0u) {
      for (UF_long k = k1; k < k2; ++k) {
        double value = work[k];
        const UF_long begin = offp[k];
        const UF_long end = offp[k + 1u];
        for (UF_long p = begin; p < end; ++p) {
          const UF_long row = offi[p];
          if (row >= k1) {
            return 0;
          }
          value -= offx[p] * work[row];
        }
        work[k] = value;
      }
    }

    if (ut->cols16 != NULL && ut->source_pos16 != NULL) {
      for (UF_long row = k1; row < k2; ++row) {
        double value = work[row];
        for (UF_long p = ut->ptr[row]; p < ut->ptr[row + 1u]; ++p) {
          value -= u_values[ut->source_pos16[p]] * work[ut->cols16[p]];
        }
        work[row] = value / udiag[row];
      }
    } else {
      for (UF_long row = k1; row < k2; ++row) {
        double value = work[row];
        for (UF_long p = ut->ptr[row]; p < ut->ptr[row + 1u]; ++p) {
          value -= u_values[ut->source_pos[p]] * work[ut->cols[p]];
        }
        work[row] = value / udiag[row];
      }
    }
    if (lt->cols16 != NULL && lt->source_pos16 != NULL) {
      for (UF_long remaining = k2; remaining > k1; --remaining) {
        const UF_long row = remaining - 1u;
        double value = work[row];
        for (UF_long p = lt->ptr[row]; p < lt->ptr[row + 1u]; ++p) {
          value -= l_values[lt->source_pos16[p]] * work[lt->cols16[p]];
        }
        work[row] = value;
      }
    } else {
      for (UF_long remaining = k2; remaining > k1; --remaining) {
        const UF_long row = remaining - 1u;
        double value = work[row];
        for (UF_long p = lt->ptr[row]; p < lt->ptr[row + 1u]; ++p) {
          value -= l_values[lt->source_pos[p]] * work[lt->cols[p]];
        }
        work[row] = value;
      }
    }
  }
  for (UF_long k = 0u; k < n; ++k) {
    const double result = work[k];
    x[pnum[k]] = rs != NULL ? result / rs[k] : result;
  }
  return 1;
}

static int kls_try_row_refactor_solve(kls_solver *solver,
                                      int kernel_transpose,
                                      int64_t nrhs,
                                      double *x,
                                      int64_t ldx) {
  if (!kls_row_refactor_solve_is_eligible(solver)) {
    return 0;
  }

  const UF_long n = solver->n;
  const UF_long nblocks = solver->symbolic->nblocks;
  const UF_long *r = solver->symbolic->R;
  const UF_long *pnum = solver->numeric->Pnum;
  const UF_long *q = solver->symbolic->Q;
  const double *udiag = (const double *)solver->numeric->Udiag;
  const double *rs = solver->numeric->Rs;
  const UF_long *offp = solver->numeric->Offp;
  const UF_long *offi = solver->numeric->Offi;
  const double *offx = (const double *)solver->numeric->Offx;
  double *work = (double *)solver->numeric->Xwork;

  if (!kls_validate_row_refactor_solve_metadata(solver)) {
    return 0;
  }

  solver->common.status = TRILINOS_KLU_OK;

  if (nrhs == 1) {
    if (!kls_try_row_refactor_solve_one_rhs(solver, kernel_transpose, x)) {
      return 0;
    }
    kls_record_row_refactor_row_solve(solver);
    return 1;
  }

  for (int64_t chunk = 0; chunk < nrhs; chunk += 4) {
    const UF_long rhs_count =
      (UF_long)(nrhs - chunk < 4 ? nrhs - chunk : 4);
    if (!kernel_transpose) {
      for (UF_long k = 0; k < n; ++k) {
        double *wk = work + k * rhs_count;
        for (UF_long rhs_offset = 0; rhs_offset < rhs_count; ++rhs_offset) {
          const double *rhs = x + (chunk + (int64_t)rhs_offset) * ldx;
          const double rhs_value = rhs[pnum[k]];
          wk[rhs_offset] = rs != NULL ? rhs_value / rs[k] : rhs_value;
        }
      }

      for (UF_long remaining_block = nblocks; remaining_block > 0u;
           --remaining_block) {
        const UF_long block = remaining_block - 1u;
        const UF_long k1 = r[block];
        const UF_long k2 = r[block + 1u];
        for (UF_long row = k1; row < k2;) {
          UF_long next_row = row;
          if (kls_try_row_refactor_compact_panel_lower_solve_group_multi(
                solver, k1, k2, row, work, rhs_count, &next_row) &&
              next_row > row) {
            row = next_row;
            continue;
          }
          double *wrow = work + row * rhs_count;
          const UF_long begin = solver->row_refactor_l_ptr[row];
          const UF_long end = solver->row_refactor_l_ptr[row + 1u];
          for (UF_long p = begin; p < end; ++p) {
            const UF_long col = solver->row_refactor_l_cols[p];
            if (col < k1 || col >= row) {
              return 0;
            }
            const double lij =
              kls_row_refactor_l_solve_row_value(solver, row, p);
            const double *wcol = work + col * rhs_count;
            for (UF_long rhs_offset = 0; rhs_offset < rhs_count;
                 ++rhs_offset) {
              wrow[rhs_offset] -= lij * wcol[rhs_offset];
            }
          }
          row++;
        }

        for (UF_long remaining = k2; remaining > k1;) {
          UF_long next_remaining = remaining;
          if (kls_try_row_refactor_compact_panel_upper_solve_group_multi(
                solver, k1, k2, remaining, work, rhs_count, udiag,
                &next_remaining) &&
              next_remaining < remaining) {
            remaining = next_remaining;
            continue;
          }
          const UF_long row = remaining - 1u;
          double *wrow = work + row * rhs_count;
          const UF_long begin = solver->row_refactor_u_ptr[row];
          const UF_long end = solver->row_refactor_u_ptr[row + 1u];
          for (UF_long p = begin; p < end; ++p) {
            const UF_long col = solver->row_refactor_u_cols[p];
            if (col <= row || col >= k2) {
              return 0;
            }
            const double uij =
              kls_row_refactor_u_solve_row_value(solver, row, p);
            const double *wcol = work + col * rhs_count;
            for (UF_long rhs_offset = 0; rhs_offset < rhs_count;
                 ++rhs_offset) {
              wrow[rhs_offset] -= uij * wcol[rhs_offset];
            }
          }
          for (UF_long rhs_offset = 0; rhs_offset < rhs_count; ++rhs_offset) {
            wrow[rhs_offset] /= udiag[row];
          }
          remaining--;
        }

        if (block > 0u) {
          for (UF_long k = k1; k < k2; ++k) {
            const UF_long begin = offp[k];
            const UF_long end = offp[k + 1u];
            const double *wk = work + k * rhs_count;
            for (UF_long p = begin; p < end; ++p) {
              const UF_long row = offi[p];
              if (row < 0 || row >= k1) {
                return 0;
              }
              double *wrow = work + row * rhs_count;
              const double value = offx[p];
              for (UF_long rhs_offset = 0; rhs_offset < rhs_count;
                   ++rhs_offset) {
                wrow[rhs_offset] -= value * wk[rhs_offset];
              }
            }
          }
        }
      }

      for (UF_long k = 0; k < n; ++k) {
        const double *wk = work + k * rhs_count;
        for (UF_long rhs_offset = 0; rhs_offset < rhs_count; ++rhs_offset) {
          double *rhs = x + (chunk + (int64_t)rhs_offset) * ldx;
          rhs[q[k]] = wk[rhs_offset];
        }
      }
    } else {
      for (UF_long k = 0; k < n; ++k) {
        double *wk = work + k * rhs_count;
        for (UF_long rhs_offset = 0; rhs_offset < rhs_count; ++rhs_offset) {
          const double *rhs = x + (chunk + (int64_t)rhs_offset) * ldx;
          wk[rhs_offset] = rhs[q[k]];
        }
      }

      for (UF_long block = 0; block < nblocks; ++block) {
        const UF_long k1 = r[block];
        const UF_long k2 = r[block + 1u];
        if (block > 0u) {
          for (UF_long k = k1; k < k2; ++k) {
            const UF_long begin = offp[k];
            const UF_long end = offp[k + 1u];
            double *wk = work + k * rhs_count;
            for (UF_long p = begin; p < end; ++p) {
              const UF_long row = offi[p];
              if (row < 0 || row >= k1) {
                return 0;
              }
              const double *wrow = work + row * rhs_count;
              const double value = offx[p];
              for (UF_long rhs_offset = 0; rhs_offset < rhs_count;
                   ++rhs_offset) {
                wk[rhs_offset] -= value * wrow[rhs_offset];
              }
            }
          }
        }

        for (UF_long row = k1; row < k2;) {
          UF_long next_row = row;
          if (kls_try_row_refactor_compact_panel_u_transpose_solve_group_multi(
                solver, k1, k2, row, work, rhs_count, udiag,
                &next_row) &&
              next_row > row) {
            row = next_row;
            continue;
          }
          double *wrow = work + row * rhs_count;
          for (UF_long rhs_offset = 0; rhs_offset < rhs_count; ++rhs_offset) {
            wrow[rhs_offset] /= udiag[row];
          }
          const UF_long begin = solver->row_refactor_u_ptr[row];
          const UF_long end = solver->row_refactor_u_ptr[row + 1u];
          for (UF_long p = begin; p < end; ++p) {
            const UF_long col = solver->row_refactor_u_cols[p];
            if (col <= row || col >= k2) {
              return 0;
            }
            const double uij =
              kls_row_refactor_u_solve_row_value(solver, row, p);
            double *wcol = work + col * rhs_count;
            for (UF_long rhs_offset = 0; rhs_offset < rhs_count;
                 ++rhs_offset) {
              wcol[rhs_offset] -= uij * wrow[rhs_offset];
            }
          }
          row++;
        }

        for (UF_long remaining = k2; remaining > k1;) {
          UF_long next_remaining = remaining;
          if (kls_try_row_refactor_compact_panel_l_transpose_solve_group_multi(
                solver, k1, k2, remaining, work, rhs_count,
                &next_remaining) &&
              next_remaining < remaining) {
            remaining = next_remaining;
            continue;
          }
          const UF_long row = remaining - 1u;
          const double *wrow = work + row * rhs_count;
          const UF_long begin = solver->row_refactor_l_ptr[row];
          const UF_long end = solver->row_refactor_l_ptr[row + 1u];
          for (UF_long p = begin; p < end; ++p) {
            const UF_long col = solver->row_refactor_l_cols[p];
            if (col < k1 || col >= row) {
              return 0;
            }
            const double lij =
              kls_row_refactor_l_solve_row_value(solver, row, p);
            double *wcol = work + col * rhs_count;
            for (UF_long rhs_offset = 0; rhs_offset < rhs_count;
                 ++rhs_offset) {
              wcol[rhs_offset] -= lij * wrow[rhs_offset];
            }
          }
          remaining--;
        }
      }

      for (UF_long k = 0; k < n; ++k) {
        const double *wk = work + k * rhs_count;
        for (UF_long rhs_offset = 0; rhs_offset < rhs_count; ++rhs_offset) {
          double *rhs = x + (chunk + (int64_t)rhs_offset) * ldx;
          const double result = wk[rhs_offset];
          rhs[pnum[k]] = rs != NULL ? result / rs[k] : result;
        }
      }
    }
  }

  kls_record_row_refactor_row_solve(solver);
  return 1;
}

static int kls_clear_row_refactor_segment_target_residuals(kls_solver *solver,
                                                           double *x,
                                                           UF_long row) {
  if (solver == NULL || x == NULL || row >= solver->n ||
      solver->row_refactor_segment_input_row_ready == NULL ||
      !solver->row_refactor_segment_input_row_ready[row] ||
      solver->row_refactor_segment_input_cleanup_ptr == NULL ||
      solver->row_refactor_segment_input_cleanup_cols == NULL) {
    return 0;
  }
  const UF_long begin = solver->row_refactor_segment_input_cleanup_ptr[row];
  const UF_long end =
    solver->row_refactor_segment_input_cleanup_ptr[row + 1u];
  if (begin >= end) {
    return 0;
  }
  if (end > solver->row_refactor_segment_input_cleanup_entries) {
    return 0;
  }
  for (UF_long p = begin; p < end; ++p) {
    const UF_long col = solver->row_refactor_segment_input_cleanup_cols[p];
    if (col >= solver->n) {
      return 0;
    }
    x[col] = 0.0;
  }
  kls_record_row_refactor_segment_target_cleanup(solver, end - begin);
  return 1;
}

static void kls_clear_row_refactor_input_residuals(kls_solver *solver,
                                                   double *x,
                                                   UF_long row) {
  if (solver == NULL || x == NULL ||
      solver->row_refactor_input_ptr == NULL ||
      solver->row_refactor_input_cols == NULL ||
      row >= solver->n) {
    return;
  }
  if (solver->row_refactor_input_needs_cleanup != NULL &&
      !solver->row_refactor_input_needs_cleanup[row]) {
    return;
  }
  if (kls_clear_row_refactor_segment_target_residuals(solver, x, row)) {
    return;
  }
  for (UF_long p = solver->row_refactor_input_ptr[row];
       p < solver->row_refactor_input_ptr[row + 1u]; ++p) {
    x[solver->row_refactor_input_cols[p]] = 0.0;
  }
}

static int kls_scatter_row_refactor_row_values(kls_solver *solver,
                                               UF_long row) {
  if (solver == NULL || row >= solver->n ||
      solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_l_values == NULL ||
      solver->row_refactor_l_row_values == NULL ||
      solver->row_refactor_u_ptr == NULL ||
      solver->row_refactor_u_values == NULL ||
      solver->row_refactor_u_row_values == NULL) {
    return 0;
  }
  for (UF_long p = solver->row_refactor_l_ptr[row];
       p < solver->row_refactor_l_ptr[row + 1u]; ++p) {
    *solver->row_refactor_l_values[p] =
      solver->row_refactor_l_row_values[p];
  }
  for (UF_long p = solver->row_refactor_u_ptr[row];
       p < solver->row_refactor_u_ptr[row + 1u]; ++p) {
    *solver->row_refactor_u_values[p] =
      solver->row_refactor_u_row_values[p];
  }
  return 1;
}

static KLS_ALWAYS_INLINE int
kls_row_refactor_multiplier_rejects(double lij, double tol) {
  if (tol <= DBL_MIN) {
    return 0;
  }
  const double lij_abs = fabs(lij);
  return !isfinite(lij_abs) || lij_abs * tol > 1.0 + 1.0e-12;
}

static int kls_row_refactor_pivot_rejects(double pivot_abs,
                                          double row_max_abs,
                                          double tol,
                                          double *ratio_out) {
  if (ratio_out != NULL) {
    *ratio_out = -1.0;
  }
  if (tol <= DBL_MIN) {
    return 0;
  }
  if (!isfinite(row_max_abs)) {
    if (ratio_out != NULL) {
      *ratio_out = DBL_MAX;
    }
    return 1;
  }
  if (row_max_abs <= 0.0) {
    return 0;
  }
  if (!isfinite(pivot_abs) || pivot_abs <= DBL_MIN) {
    if (ratio_out != NULL) {
      *ratio_out = DBL_MAX;
    }
    return 1;
  }
  const double ratio = row_max_abs / pivot_abs;
  if (ratio_out != NULL) {
    *ratio_out = isfinite(ratio) ? ratio : DBL_MAX;
  }
  return !isfinite(ratio) || ratio * tol > 1.0 + 1.0e-12;
}

static int kls_numeric_rowwise_u_pivots_pass_threshold_from_block(
  const kls_solver *solver,
  UF_long start_block,
  UF_long *rejected_pivot_out,
  UF_long *rejected_pivot_col_out,
  UF_long *rejected_row_out,
  double *rejected_ratio_abs_out,
  double *rejected_pivot_abs_out,
  double *rejected_row_max_abs_out) {
  if (rejected_pivot_out != NULL) {
    *rejected_pivot_out = KLS_KLU_EMPTY;
  }
  if (rejected_pivot_col_out != NULL) {
    *rejected_pivot_col_out = KLS_KLU_EMPTY;
  }
  if (rejected_row_out != NULL) {
    *rejected_row_out = KLS_KLU_EMPTY;
  }
  if (rejected_ratio_abs_out != NULL) {
    *rejected_ratio_abs_out = -1.0;
  }
  if (rejected_pivot_abs_out != NULL) {
    *rejected_pivot_abs_out = -1.0;
  }
  if (rejected_row_max_abs_out != NULL) {
    *rejected_row_max_abs_out = -1.0;
  }
  if (solver == NULL || solver->symbolic == NULL ||
      solver->numeric == NULL) {
    return -1;
  }
  const double tol = solver->common.tol;
  if (tol <= DBL_MIN) {
    return 1;
  }

  const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  const trilinos_klu_l_numeric *numeric = solver->numeric;
  if (symbolic->R == NULL || numeric->Udiag == NULL ||
      numeric->Uip == NULL || numeric->Ulen == NULL ||
      numeric->LUbx == NULL || start_block > symbolic->nblocks ||
      symbolic->nblocks > numeric->nblocks) {
    return -1;
  }

  const double *udiag = (const double *)numeric->Udiag;
  for (UF_long block = start_block; block < symbolic->nblocks; ++block) {
    const UF_long k1 = symbolic->R[block];
    const UF_long k2 = symbolic->R[block + 1u];
    if (k2 < k1 || k2 > solver->n) {
      return -1;
    }
    const UF_long nk = k2 - k1;
    if (nk == 0u) {
      continue;
    }
    if (nk > (UF_long)(SIZE_MAX / sizeof(double))) {
      return -1;
    }

    double *row_max =
      (double *)malloc((size_t)nk * sizeof(*row_max));
    if (row_max == NULL) {
      return -1;
    }
    for (UF_long k = 0; k < nk; ++k) {
      const double pivot_abs = fabs(udiag[k1 + k]);
      row_max[k] = isfinite(pivot_abs) ? pivot_abs : DBL_MAX;
    }

    if (nk > 1u) {
      double *lu = (double *)numeric->LUbx[block];
      if (lu == NULL) {
        free(row_max);
        return -1;
      }
      const UF_long *uip = numeric->Uip + k1;
      const UF_long *ulen = numeric->Ulen + k1;
      for (UF_long local_col = 0; local_col < nk; ++local_col) {
        UF_long *ui = NULL;
        double *ux = NULL;
        UF_long ucol_len = 0;
        kls_klu_get_pointer(lu, uip, ulen, local_col, &ui, &ux, &ucol_len);
        if ((ucol_len > 0u && (ui == NULL || ux == NULL)) ||
            ucol_len > nk) {
          free(row_max);
          return -1;
        }
        for (UF_long p = 0; p < ucol_len; ++p) {
          const UF_long local_row = ui[p];
          if (local_row >= nk || local_row >= local_col) {
            free(row_max);
            return -1;
          }
          const double value_abs = fabs(ux[p]);
          if (!isfinite(value_abs)) {
            row_max[local_row] = DBL_MAX;
          } else if (value_abs > row_max[local_row]) {
            row_max[local_row] = value_abs;
          }
        }
      }
    }

    for (UF_long local_row = 0; local_row < nk; ++local_row) {
      const UF_long global_row = k1 + local_row;
      const double pivot_abs = fabs(udiag[global_row]);
      double ratio_abs = -1.0;
      if (kls_row_refactor_pivot_rejects(
            pivot_abs, row_max[local_row], tol, &ratio_abs)) {
        if (rejected_pivot_out != NULL) {
          *rejected_pivot_out = global_row;
        }
        if (rejected_pivot_col_out != NULL) {
          *rejected_pivot_col_out =
            symbolic->Q != NULL ? symbolic->Q[global_row] : KLS_KLU_EMPTY;
        }
        if (rejected_row_out != NULL) {
          *rejected_row_out = global_row;
        }
        if (rejected_ratio_abs_out != NULL) {
          *rejected_ratio_abs_out = ratio_abs;
        }
        if (rejected_pivot_abs_out != NULL) {
          *rejected_pivot_abs_out =
            isfinite(pivot_abs) ? pivot_abs : DBL_MAX;
        }
        if (rejected_row_max_abs_out != NULL) {
          *rejected_row_max_abs_out = row_max[local_row];
        }
        free(row_max);
        return 0;
      }
    }
    free(row_max);
  }
  return 1;
}

static int kls_checked_refactor_accepts_rowwise_u(
  kls_solver *solver,
  UF_long start_block) {
  if (solver == NULL) {
    return 0;
  }
  if (solver->common.status < 0 ||
      solver->common.status == TRILINOS_KLU_SINGULAR) {
    return 1;
  }

  UF_long rejected_pivot = KLS_KLU_EMPTY;
  UF_long rejected_pivot_col = KLS_KLU_EMPTY;
  UF_long rejected_row = KLS_KLU_EMPTY;
  double rejected_ratio_abs = -1.0;
  double rejected_pivot_abs = -1.0;
  double rejected_row_max_abs = -1.0;
  const int status =
    kls_numeric_rowwise_u_pivots_pass_threshold_from_block(
      solver, start_block, &rejected_pivot, &rejected_pivot_col,
      &rejected_row, &rejected_ratio_abs, &rejected_pivot_abs,
      &rejected_row_max_abs);
  if (status > 0) {
    return 1;
  }
  if (status < 0 || rejected_pivot == KLS_KLU_EMPTY) {
    kls_record_fast_factor_failure(
      solver, KLS_FAST_FACTOR_FAIL_ROWWISE_U_INVALID,
      solver->common.status);
    solver->common.status = TRILINOS_KLU_INVALID;
    return 0;
  }

  solver->fast_reject_refresh_state = KLS_FAST_REJECT_REFRESH_ALL;
  kls_clear_fast_reject_tail_seed(solver);
  kls_record_fast_reject_detail(solver, rejected_pivot,
                                rejected_pivot_col, rejected_row,
                                rejected_ratio_abs, rejected_pivot_abs,
                                rejected_row_max_abs);
  solver->common.status = TRILINOS_KLU_OK;
  return 0;
}

static int kls_row_refactor_u_row_abs_max_from_direct(
  const kls_solver *solver,
  UF_long row,
  double pivot,
  const double *x,
  const double *row_u_values,
  int direct_input,
  double *row_max_abs_out) {
  if (row_max_abs_out != NULL) {
    *row_max_abs_out = -1.0;
  }
  if (solver == NULL || x == NULL || row_max_abs_out == NULL ||
      row >= solver->n || solver->row_refactor_u_ptr == NULL ||
      solver->row_refactor_u_cols == NULL ||
      (direct_input && row_u_values == NULL)) {
    return 0;
  }

  double row_max_abs = fabs(pivot);
  if (!isfinite(row_max_abs)) {
    *row_max_abs_out = DBL_MAX;
    return 1;
  }
  const UF_long begin = solver->row_refactor_u_ptr[row];
  const UF_long end = solver->row_refactor_u_ptr[row + 1u];
  if (begin > end) {
    return 0;
  }
  for (UF_long p = begin; p < end; ++p) {
    const UF_long offset = p - begin;
    const UF_long col = solver->row_refactor_u_cols[p];
    if (col >= solver->n) {
      return 0;
    }
    const double value =
      (direct_input ? row_u_values[offset] : 0.0) + x[col];
    const double value_abs = fabs(value);
    if (!isfinite(value_abs)) {
      *row_max_abs_out = DBL_MAX;
      return 1;
    }
    if (value_abs > row_max_abs) {
      row_max_abs = value_abs;
    }
  }
  *row_max_abs_out = row_max_abs;
  return 1;
}

static double kls_row_refactor_abs_max_from_values(double pivot,
                                                   const double *values,
                                                   UF_long count) {
  double row_max_abs = fabs(pivot);
  if (!isfinite(row_max_abs)) {
    return DBL_MAX;
  }
  for (UF_long p = 0; p < count; ++p) {
    const double value_abs = fabs(values[p]);
    if (!isfinite(value_abs)) {
      return DBL_MAX;
    }
    if (value_abs > row_max_abs) {
      row_max_abs = value_abs;
    }
  }
  return row_max_abs;
}

static int kls_find_uflong_sorted(const UF_long *cols,
                                  UF_long len,
                                  UF_long target,
                                  UF_long *pos_out) {
  if (pos_out != NULL) {
    *pos_out = KLS_KLU_EMPTY;
  }
  if (cols == NULL || len == 0u) {
    return 0;
  }
  UF_long lo = 0;
  UF_long hi = len;
  while (lo < hi) {
    const UF_long mid = lo + (hi - lo) / 2u;
    const UF_long col = cols[mid];
    if (col == target) {
      if (pos_out != NULL) {
        *pos_out = mid;
      }
      return 1;
    }
    if (col < target) {
      lo = mid + 1u;
    } else {
      hi = mid;
    }
  }
  return 0;
}

static int kls_parallel_row_refactor_process_group(
  kls_egraph_refactor_worker *worker,
  UF_long group,
  int wait_for_dependencies);
static double *kls_row_refactor_compact_group_panel(kls_solver *solver,
                                                    UF_long group,
                                                    UF_long entry_count);
static double *kls_egraph_worker_supernode_workspace(
  kls_egraph_refactor_worker *worker,
  UF_long entry_count);
static unsigned char *kls_egraph_worker_byte_workspace(
  kls_egraph_refactor_worker *worker,
  UF_long entry_count);
static int kls_row_refactor_try_compact_supernode_update(
  kls_egraph_refactor_worker *worker,
  UF_long row,
  UF_long row_group,
  UF_long *p_io,
  UF_long l_begin,
  UF_long l_end,
  int wait_for_dependencies,
  const unsigned char *applied);
static int kls_row_refactor_group_trailing_slice(
  const kls_solver *solver,
  UF_long group,
  UF_long row_end,
  UF_long *len_out,
  const UF_long **cols_out);
static int kls_row_refactor_segment_input_targets_available(
  const kls_solver *solver,
  UF_long row,
  UF_long input_begin,
  UF_long input_end);
static int kls_row_refactor_apply_segment_input_target(
  kls_egraph_refactor_shared *shared,
  double *x,
  UF_long row,
  UF_long input_pos,
  UF_long l_begin,
  UF_long l_end,
  UF_long row_u_len,
  double *row_u_values,
  double *udiag,
  double value);
static int kls_parallel_row_refactor_load_targeted_input_row(
  kls_egraph_refactor_shared *shared,
  double *x,
  UF_long row);

static int kls_serial_row_refactor_numeric(kls_solver *solver,
                                           double *numeric_values,
                                           int check_pivots) {
  if (solver != NULL && solver->numeric_is_predicted) {
    return -1;
  }
  const int scaled = solver != NULL && solver->common.scale > 0;
  if (numeric_values == NULL || solver == NULL || solver->symbolic == NULL ||
      solver->numeric == NULL || solver->symbolic->nblocks == 0u ||
      (scaled && (solver->numeric->Rs == NULL ||
                  solver->numeric->Pnum == NULL)) ||
      (!scaled && solver->numeric->Rs != NULL) ||
      solver->numeric->Udiag == NULL || solver->numeric->Xwork == NULL ||
      !kls_build_row_refactor_pattern(solver, 0)) {
    return -1;
  }
  kls_row_refactor_ensure_sn_end(solver);

  trilinos_klu_l_common *common = &solver->common;
  trilinos_klu_l_numeric *numeric = solver->numeric;
  double *x = (double *)numeric->Xwork;
  const UF_long n = solver->n;

  if (scaled &&
      !trilinos_klu_l_scale((UF_long)common->scale, n,
                            solver->col_ptr, solver->row_idx,
                            numeric_values, numeric->Rs, NULL, common)) {
    return 0;
  }
  if (!kls_refresh_row_refactor_offdiag_values(solver, numeric_values)) {
    common->status = TRILINOS_KLU_INVALID;
    return 0;
  }

  memset(x, 0, (size_t)n * sizeof(*x));
  common->status = TRILINOS_KLU_OK;
  common->numerical_rank = KLS_KLU_EMPTY;
  common->singular_col = KLS_KLU_EMPTY;
  common->nrealloc = 0;
  const int defer_value_scatter =
    kls_row_refactor_should_defer_value_scatter(solver, check_pivots);
  const int lazy_value_scatter = !check_pivots;
  const int native_row_panel_state =
    kls_native_row_panel_effective_state(solver);
  const UF_long native_row_panel_count_before =
    solver->row_refactor_native_row_panel_count;
  const UF_long compact_supernode_update_count_before =
    solver->row_refactor_compact_supernode_update_count;
  const UF_long partial_supernode_pipeline_count_before =
    solver->row_refactor_partial_supernode_pipeline_run_count;
  kls_row_refactor_record_pipeline_scope(solver,
                                         solver->row_refactor_level_count);
  kls_record_row_refactor_run(solver, check_pivots, 0);
  kls_reset_row_refactor_compact_panel_valid(solver);
  if (defer_value_scatter) {
    kls_record_row_refactor_defer_value_scatter_run(solver);
  }

  kls_egraph_refactor_shared shared;
  memset(&shared, 0, sizeof(shared));
  if (pthread_mutex_init(&shared.lock, NULL) != 0) {
    common->status = TRILINOS_KLU_INVALID;
    return 0;
  }
  shared.solver = solver;
  shared.values = numeric_values;
  shared.rs = scaled ? numeric->Rs : NULL;
  shared.check_pivots = check_pivots;
  shared.scale = (int)common->scale;
  shared.kernel = KLS_EGRAPH_REFACTOR_KERNEL_GENERIC;
  shared.thread_count = 1;
  shared.row_refactor_native_row_panel_state = native_row_panel_state;

  atomic_init(&shared.stop, 0);
  shared.rejected_pivot = KLS_KLU_EMPTY;
  shared.rejected_pivot_col = KLS_KLU_EMPTY;
  shared.rejected_row = KLS_KLU_EMPTY;
  shared.rejected_multiplier_abs = -1.0;
  shared.rejected_pivot_abs = -1.0;
  shared.rejected_candidate_abs = -1.0;
  shared.numerical_rank = UF_long_max;
  shared.singular_col = KLS_KLU_EMPTY;
  shared.row_refactor_mode = 1;
  shared.row_refactor_defer_value_scatter = defer_value_scatter;
  shared.row_refactor_lazy_value_scatter = lazy_value_scatter;
  shared.row_refactor_blocked_trailing_update = 1;
  shared.row_refactor_native_row_panel_active = 0;
  solver->row_refactor_native_row_panel_enabled =
    kls_native_row_panel_enabled_for_run(
      solver, native_row_panel_state, check_pivots);
  atomic_init(&shared.row_prefactor_rows, 0ul);
  atomic_init(&shared.row_prefactor_deps, 0ul);
  atomic_init(&shared.row_prefactor_supernode_rows, 0ul);
  atomic_init(&shared.row_prefactor_supernode_deps, 0ul);

  kls_egraph_refactor_worker worker;
  memset(&worker, 0, sizeof(worker));
  worker.shared = &shared;
  worker.tid = 0;
  worker.x = x;

  int groups_ok = 1;
  for (UF_long group = 0; group < solver->row_refactor_group_count; ++group) {
    if (!kls_parallel_row_refactor_process_group(&worker, group, 0)) {
      groups_ok = 0;
      break;
    }
  }
  free(worker.segment_panel);
  free(worker.supernode_workspace);
  free(worker.index_workspace);
  free(worker.byte_workspace);
  free(worker.object_workspace);

  if (!groups_ok || shared.invalid || shared.pivot_rejected ||
      (shared.singular && common->halt_if_singular)) {
    memset(x, 0, (size_t)n * sizeof(*x));
  }

  if (!groups_ok && !shared.invalid && !shared.pivot_rejected &&
      !(shared.singular && common->halt_if_singular)) {
    shared.invalid = 1;
  }
  if (shared.invalid) {
    pthread_mutex_destroy(&shared.lock);
    common->status = TRILINOS_KLU_INVALID;
    return 0;
  }
  if (shared.pivot_rejected) {
    const UF_long rejected_pivot = shared.rejected_pivot;
    if (defer_value_scatter &&
        !kls_scatter_row_refactor_prefix_values(solver, rejected_pivot)) {
      pthread_mutex_destroy(&shared.lock);
      common->status = TRILINOS_KLU_INVALID;
      return 0;
    }
    solver->fast_reject_refresh_state = KLS_FAST_REJECT_REFRESH_PREFIX;
    kls_record_fast_reject_block_suffix_tail_seed(solver, rejected_pivot);
    kls_record_fast_reject_detail(solver, rejected_pivot,
                                  shared.rejected_pivot_col,
                                  shared.rejected_row,
                                  shared.rejected_multiplier_abs,
                                  shared.rejected_pivot_abs,
                                  shared.rejected_candidate_abs);
    kls_fill_fast_reject_row_tail_candidate_stats(
      solver, numeric_values, rejected_pivot, x);
    pthread_mutex_destroy(&shared.lock);
    common->status = TRILINOS_KLU_OK;
    return 0;
  }
  if (shared.singular) {
    common->status = TRILINOS_KLU_SINGULAR;
    common->numerical_rank = shared.numerical_rank;
    common->singular_col = shared.singular_col;
    if (common->halt_if_singular) {
      pthread_mutex_destroy(&shared.lock);
      return 0;
    }
  } else {
    common->status = TRILINOS_KLU_OK;
  }
  pthread_mutex_destroy(&shared.lock);

  if (defer_value_scatter) {
    if (lazy_value_scatter) {
      kls_record_row_refactor_lazy_value_scatter_run(solver);
    } else if (!kls_scatter_row_refactor_l_values(solver) ||
               !kls_scatter_row_refactor_u_values(solver)) {
      common->status = TRILINOS_KLU_INVALID;
      return 0;
    }
  }
  if (scaled && !kls_parallel_refactor_permute_scale(solver)) {
    common->status = TRILINOS_KLU_INVALID;
    return 0;
  }
  solver->row_refactor_values_ready = 1;
  solver->row_refactor_solve_direct_ready = 0;
  solver->row_refactor_solve_validated = 0;
  solver->numeric_needs_refinement = 0;
  solver->solve_refine_single_shot = 0;
  kls_maybe_disable_native_row_panel_auto(
    solver, check_pivots, native_row_panel_state,
    native_row_panel_count_before, compact_supernode_update_count_before,
    partial_supernode_pipeline_count_before);
  return 1;
}

static KLS_ALWAYS_INLINE int kls_parallel_row_refactor_rejects_multiplier(
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

static KLS_ALWAYS_INLINE int kls_parallel_row_refactor_rejects_pivot(
  kls_egraph_refactor_worker *worker,
  UF_long row,
  double pivot,
  double row_max_abs) {
  if (worker == NULL || worker->shared == NULL ||
      !worker->shared->check_pivots) {
    return 0;
  }
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  if (solver == NULL || solver->symbolic == NULL ||
      solver->symbolic->Q == NULL || row >= solver->n) {
    return 0;
  }
  double ratio_abs = -1.0;
  const double pivot_abs = fabs(pivot);
  if (!kls_row_refactor_pivot_rejects(pivot_abs, row_max_abs,
                                      solver->common.tol, &ratio_abs)) {
    return 0;
  }
  kls_egraph_refactor_record_reject(
    shared, row, solver->symbolic->Q[row], row, ratio_abs,
    pivot_abs, row_max_abs);
  return 1;
}

/* A row mirror groups input entries by the numeric pivot row.  Consequently
   every entry in one row uses the same original-row scaling factor.  The
   generic input helper rediscovers that row through Ai and divides once per
   nonzero; the cooperative row engine can hoist both operations to once per
   row through Pnum, just as the lean row path does. */
static KLS_ALWAYS_INLINE int kls_row_refactor_input_multiplier(
  const kls_egraph_refactor_shared *shared,
  UF_long row,
  double *multiplier_out) {
  if (shared == NULL || multiplier_out == NULL || shared->solver == NULL ||
      row >= shared->solver->n) {
    return 0;
  }
  if (shared->scale <= 0) {
    *multiplier_out = 1.0;
    return 1;
  }
  const kls_solver *solver = shared->solver;
  if (shared->rs == NULL || solver->numeric == NULL ||
      solver->numeric->Pnum == NULL) {
    return 0;
  }
  const UF_long oldrow = solver->numeric->Pnum[row];
  if (oldrow >= solver->n || shared->rs[oldrow] == 0.0) {
    return 0;
  }
  *multiplier_out = 1.0 / shared->rs[oldrow];
  return 1;
}

static KLS_ALWAYS_INLINE int kls_row_refactor_scaled_input_value(
  const kls_egraph_refactor_shared *shared,
  UF_long input_pos,
  double multiplier,
  double *value_out) {
  if (shared == NULL || shared->solver == NULL || shared->values == NULL ||
      value_out == NULL || input_pos >= shared->solver->nnz) {
    return 0;
  }
  *value_out = shared->values[input_pos] * multiplier;
  return 1;
}

static int kls_parallel_row_refactor_load_input_row(
  kls_egraph_refactor_shared *shared,
  double *x,
  UF_long row) {
  kls_solver *solver = shared != NULL ? shared->solver : NULL;
  if (solver == NULL || x == NULL || row >= solver->n ||
      solver->row_refactor_input_ptr == NULL ||
      solver->row_refactor_input_cols == NULL ||
      solver->row_refactor_input_pos == NULL) {
    if (shared != NULL) {
      kls_egraph_refactor_record_invalid(shared);
    }
    return 0;
  }
  double input_multiplier = 1.0;
  if (!kls_row_refactor_input_multiplier(shared, row, &input_multiplier)) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }
  for (UF_long p = solver->row_refactor_input_ptr[row];
       p < solver->row_refactor_input_ptr[row + 1u]; ++p) {
    const UF_long col = solver->row_refactor_input_cols[p];
    if (col >= solver->n ||
        !kls_row_refactor_scaled_input_value(
          shared, kls_row_refactor_input_pos_at(solver, p), input_multiplier,
          x + col)) {
      kls_egraph_refactor_record_invalid(shared);
      return 0;
    }
  }
  return 1;
}

static int kls_parallel_row_refactor_load_batch_input_row_targets(
  kls_egraph_refactor_shared *shared,
  UF_long row,
  UF_long input_begin,
  UF_long input_end,
  UF_long dep_count,
  double *dep_values,
  double *pivot_out,
  UF_long u_count,
  double *u_values) {
  kls_solver *solver = shared != NULL ? shared->solver : NULL;
  if (solver == NULL || row >= solver->n || pivot_out == NULL ||
      (dep_count > 0u && dep_values == NULL) ||
      (u_count > 0u && u_values == NULL) ||
      solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_segment_input_target_kind == NULL ||
      solver->row_refactor_segment_input_target_pos == NULL) {
    return 0;
  }
  if (!kls_row_refactor_segment_input_targets_available(
        solver, row, input_begin, input_end)) {
    return 0;
  }
  const UF_long l_begin = solver->row_refactor_l_ptr[row];
  const UF_long l_end = solver->row_refactor_l_ptr[row + 1u];
  if (l_begin > l_end || dep_count != l_end - l_begin) {
    return 0;
  }

  for (UF_long p = input_begin; p < input_end; ++p) {
    const unsigned char target =
      solver->row_refactor_segment_input_target_kind[p];
    const UF_long pos = solver->row_refactor_segment_input_target_pos[p];
    switch ((kls_row_segment_input_target_kind)target) {
      case KLS_ROW_SEGMENT_INPUT_TARGET_L:
        if (pos < l_begin || pos >= l_end || pos - l_begin >= dep_count) {
          kls_egraph_refactor_record_invalid(shared);
          return -1;
        }
        break;
      case KLS_ROW_SEGMENT_INPUT_TARGET_PIVOT:
        break;
      case KLS_ROW_SEGMENT_INPUT_TARGET_U:
        if (pos >= u_count) {
          kls_egraph_refactor_record_invalid(shared);
          return -1;
        }
        break;
      case KLS_ROW_SEGMENT_INPUT_TARGET_EXTERNAL:
      case KLS_ROW_SEGMENT_INPUT_TARGET_NONE:
      default:
        return 0;
    }
  }

  if (dep_count > 0u) {
    memset(dep_values, 0, (size_t)dep_count * sizeof(*dep_values));
  }
  *pivot_out = 0.0;
  if (u_count > 0u) {
    memset(u_values, 0, (size_t)u_count * sizeof(*u_values));
  }

  double input_multiplier = 1.0;
  if (!kls_row_refactor_input_multiplier(shared, row, &input_multiplier)) {
    kls_egraph_refactor_record_invalid(shared);
    return -1;
  }
  for (UF_long p = input_begin; p < input_end; ++p) {
    double value = 0.0;
    if (!kls_row_refactor_scaled_input_value(
          shared, kls_row_refactor_input_pos_at(solver, p), input_multiplier,
          &value)) {
      kls_egraph_refactor_record_invalid(shared);
      return -1;
    }
    const unsigned char target =
      solver->row_refactor_segment_input_target_kind[p];
    const UF_long pos = solver->row_refactor_segment_input_target_pos[p];
    switch ((kls_row_segment_input_target_kind)target) {
      case KLS_ROW_SEGMENT_INPUT_TARGET_L:
        dep_values[pos - l_begin] = value;
        break;
      case KLS_ROW_SEGMENT_INPUT_TARGET_PIVOT:
        *pivot_out = value;
        break;
      case KLS_ROW_SEGMENT_INPUT_TARGET_U:
        u_values[pos] = value;
        break;
      case KLS_ROW_SEGMENT_INPUT_TARGET_EXTERNAL:
      case KLS_ROW_SEGMENT_INPUT_TARGET_NONE:
      default:
        kls_egraph_refactor_record_invalid(shared);
        return -1;
    }
  }

  kls_record_row_refactor_segment_target_input(solver);
  kls_record_row_refactor_batch_direct_input(solver);
  return 1;
}

static int kls_parallel_row_refactor_load_batch_input_row(
  kls_egraph_refactor_shared *shared,
  UF_long row,
  const UF_long *dep_cols,
  UF_long dep_count,
  double *dep_values,
  double *pivot_out,
  const UF_long *u_cols,
  UF_long u_count,
  double *u_values) {
  kls_solver *solver = shared != NULL ? shared->solver : NULL;
  if (solver == NULL || row >= solver->n || pivot_out == NULL ||
      (dep_count > 0u && (dep_cols == NULL || dep_values == NULL)) ||
      (u_count > 0u && (u_cols == NULL || u_values == NULL)) ||
      solver->row_refactor_input_ptr == NULL ||
      solver->row_refactor_input_cols == NULL ||
      solver->row_refactor_input_pos == NULL) {
    if (shared != NULL) {
      kls_egraph_refactor_record_invalid(shared);
    }
    return -1;
  }

  const UF_long input_begin = solver->row_refactor_input_ptr[row];
  const UF_long input_end = solver->row_refactor_input_ptr[row + 1u];
  {
    const int target_status =
      kls_parallel_row_refactor_load_batch_input_row_targets(
        shared, row, input_begin, input_end, dep_count, dep_values,
        pivot_out, u_count, u_values);
    if (target_status != 0) {
      return target_status;
    }
  }
  for (UF_long p = input_begin; p < input_end; ++p) {
    const UF_long col = solver->row_refactor_input_cols[p];
    if (col >= solver->n) {
      kls_egraph_refactor_record_invalid(shared);
      return -1;
    }
    if (col == row) {
      continue;
    }
    UF_long pos = KLS_KLU_EMPTY;
    if ((dep_count > 0u &&
         kls_find_uflong_sorted(dep_cols, dep_count, col, &pos)) ||
        (u_count > 0u &&
         kls_find_uflong_sorted(u_cols, u_count, col, &pos))) {
      continue;
    }
    return 0;
  }

  if (dep_count > 0u) {
    memset(dep_values, 0, (size_t)dep_count * sizeof(*dep_values));
  }
  *pivot_out = 0.0;
  if (u_count > 0u) {
    memset(u_values, 0, (size_t)u_count * sizeof(*u_values));
  }

  double input_multiplier = 1.0;
  if (!kls_row_refactor_input_multiplier(shared, row, &input_multiplier)) {
    kls_egraph_refactor_record_invalid(shared);
    return -1;
  }
  for (UF_long p = input_begin; p < input_end; ++p) {
    const UF_long col = solver->row_refactor_input_cols[p];
    double value = 0.0;
    if (!kls_row_refactor_scaled_input_value(
          shared, kls_row_refactor_input_pos_at(solver, p), input_multiplier,
          &value)) {
      kls_egraph_refactor_record_invalid(shared);
      return -1;
    }
    if (col == row) {
      *pivot_out = value;
      continue;
    }
    UF_long pos = KLS_KLU_EMPTY;
    if (dep_count > 0u &&
        kls_find_uflong_sorted(dep_cols, dep_count, col, &pos) &&
        pos < dep_count) {
      dep_values[pos] = value;
    } else if (u_count > 0u &&
               kls_find_uflong_sorted(u_cols, u_count, col, &pos) &&
               pos < u_count) {
      u_values[pos] = value;
    } else {
      kls_egraph_refactor_record_invalid(shared);
      return -1;
    }
  }

  kls_record_row_refactor_batch_direct_input(solver);
  return 1;
}

static int kls_row_refactor_apply_compact_dense_input_target(
  kls_egraph_refactor_shared *shared,
  double *x,
  UF_long row,
  UF_long input_pos,
  UF_long row_begin,
  UF_long row_end,
  UF_long trailing_len,
  double *row_dense_panel,
  double *row_panel,
  double value) {
  kls_solver *solver = shared != NULL ? shared->solver : NULL;
  if (solver == NULL || x == NULL || row < row_begin || row >= row_end ||
      row_dense_panel == NULL || row >= solver->n ||
      (trailing_len > 0u && row_panel == NULL) ||
      solver->row_refactor_segment_input_target_kind == NULL ||
      solver->row_refactor_segment_input_target_pos == NULL ||
      solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_l_cols == NULL) {
    if (shared != NULL) {
      kls_egraph_refactor_record_invalid(shared);
    }
    return 0;
  }
  const UF_long width = row_end - row_begin;
  const UF_long local_row = row - row_begin;
  const UF_long dense_len = row_end - row - 1u;
  const unsigned char target =
    solver->row_refactor_segment_input_target_kind[input_pos];
  const UF_long pos =
    solver->row_refactor_segment_input_target_pos[input_pos];
  switch ((kls_row_segment_input_target_kind)target) {
    case KLS_ROW_SEGMENT_INPUT_TARGET_EXTERNAL:
      if (pos >= solver->n) {
        kls_egraph_refactor_record_invalid(shared);
        return 0;
      }
      x[pos] = value;
      return 1;
    case KLS_ROW_SEGMENT_INPUT_TARGET_L:
      if (pos < solver->row_refactor_l_ptr[row] ||
          pos >= solver->row_refactor_l_ptr[row + 1u]) {
        kls_egraph_refactor_record_invalid(shared);
        return 0;
      }
      {
        const UF_long col = solver->row_refactor_l_cols[pos];
        if (col < row_begin || col >= row || col - row_begin >= width) {
          kls_egraph_refactor_record_invalid(shared);
          return 0;
        }
        row_dense_panel[col - row_begin] = value;
      }
      return 1;
    case KLS_ROW_SEGMENT_INPUT_TARGET_PIVOT:
      row_dense_panel[local_row] = value;
      return 1;
    case KLS_ROW_SEGMENT_INPUT_TARGET_U:
      if (pos < dense_len) {
        row_dense_panel[local_row + 1u + pos] = value;
        return 1;
      }
      if (pos >= dense_len + trailing_len || row_panel == NULL) {
        kls_egraph_refactor_record_invalid(shared);
        return 0;
      }
      row_panel[pos - dense_len] = value;
      return 1;
    case KLS_ROW_SEGMENT_INPUT_TARGET_NONE:
    default:
      kls_egraph_refactor_record_invalid(shared);
      return 0;
  }
}

static int kls_parallel_row_refactor_load_compact_dense_input_row(
  kls_egraph_refactor_shared *shared,
  double *x,
  UF_long row_begin,
  UF_long row_end,
  UF_long trailing_len,
  const UF_long *trailing_cols,
  double *dense_panel,
  double *trailing_panel,
  UF_long row) {
  kls_solver *solver = shared != NULL ? shared->solver : NULL;
  if (solver == NULL || x == NULL || row < row_begin || row >= row_end ||
      row_end > solver->n || dense_panel == NULL ||
      (trailing_len > 0u && (trailing_cols == NULL ||
                             trailing_panel == NULL)) ||
      solver->row_refactor_input_ptr == NULL ||
      solver->row_refactor_input_cols == NULL ||
      solver->row_refactor_input_pos == NULL) {
    if (shared != NULL) {
      kls_egraph_refactor_record_invalid(shared);
    }
    return -1;
  }

  const UF_long width = row_end - row_begin;
  const UF_long local_row = row - row_begin;
  double *row_dense_panel = dense_panel + local_row * width;
  double *row_panel = trailing_len > 0u
    ? trailing_panel + local_row * trailing_len : NULL;
  const UF_long input_begin = solver->row_refactor_input_ptr[row];
  const UF_long input_end = solver->row_refactor_input_ptr[row + 1u];
  const int use_retained_targets =
    kls_row_refactor_segment_input_targets_available(
      solver, row, input_begin, input_end);

  if (!use_retained_targets) {
    for (UF_long p = input_begin; p < input_end; ++p) {
      const UF_long col = solver->row_refactor_input_cols[p];
      if (col >= solver->n) {
        kls_egraph_refactor_record_invalid(shared);
        return -1;
      }
      if (col < row_begin || col < row_end) {
        continue;
      }
      UF_long trailing_pos = KLS_KLU_EMPTY;
      if (trailing_len == 0u ||
          !kls_find_uflong_sorted(trailing_cols, trailing_len, col,
                                  &trailing_pos)) {
        return 0;
      }
    }
  }

  memset(row_dense_panel, 0, (size_t)width * sizeof(*row_dense_panel));
  if (trailing_len > 0u) {
    memset(row_panel, 0, (size_t)trailing_len * sizeof(*row_panel));
  }

  double input_multiplier = 1.0;
  if (!kls_row_refactor_input_multiplier(shared, row, &input_multiplier)) {
    kls_egraph_refactor_record_invalid(shared);
    return -1;
  }
  for (UF_long p = input_begin; p < input_end; ++p) {
    const UF_long col = solver->row_refactor_input_cols[p];
    double value = 0.0;
    if (!kls_row_refactor_scaled_input_value(
          shared, kls_row_refactor_input_pos_at(solver, p), input_multiplier,
          &value)) {
      kls_egraph_refactor_record_invalid(shared);
      return -1;
    }
    if (use_retained_targets) {
      if (!kls_row_refactor_apply_compact_dense_input_target(
            shared, x, row, p, row_begin, row_end, trailing_len,
            row_dense_panel, row_panel, value)) {
        return -1;
      }
    } else {
      if (col < row_begin) {
        x[col] = value;
      } else if (col < row_end) {
        row_dense_panel[col - row_begin] = value;
      } else {
        UF_long trailing_pos = KLS_KLU_EMPTY;
        if (!kls_find_uflong_sorted(trailing_cols, trailing_len, col,
                                    &trailing_pos) ||
            trailing_pos >= trailing_len) {
          kls_egraph_refactor_record_invalid(shared);
          return -1;
        }
        row_panel[trailing_pos] = value;
      }
    }
  }
  if (use_retained_targets) {
    kls_record_row_refactor_segment_target_input(solver);
  }
  kls_record_row_refactor_compact_dense_panel_direct_input(solver);
  return 1;
}

static int kls_row_refactor_dependency_done_now(
  const kls_egraph_refactor_shared *shared,
  UF_long dep) {
  if (shared == NULL || shared->pipeline_done == NULL) {
    return 1;
  }
  return atomic_load_explicit(&shared->pipeline_done[dep],
                              memory_order_acquire) ==
         shared->pipeline_generation;
}

static int kls_row_refactor_u_row_has_col(const kls_solver *solver,
                                          UF_long dep,
                                          UF_long col) {
  if (solver == NULL || dep >= solver->n ||
      solver->row_refactor_u_ptr == NULL ||
      solver->row_refactor_u_cols == NULL) {
    return 1;
  }
  const UF_long begin = solver->row_refactor_u_ptr[dep];
  const UF_long end = solver->row_refactor_u_ptr[dep + 1u];
  if (begin > end) {
    return 1;
  }
  for (UF_long p = begin; p < end; ++p) {
    if (solver->row_refactor_u_cols[p] == col) {
      return 1;
    }
  }
  return 0;
}

static int kls_row_refactor_prefactor_dep_safe_with_run(
  const kls_solver *solver,
  const kls_egraph_refactor_shared *shared,
  UF_long l_begin,
  UF_long run_begin,
  UF_long pos,
  const unsigned char *applied,
  UF_long dep) {
  if (solver == NULL || shared == NULL || applied == NULL ||
      solver->row_refactor_l_cols == NULL ||
      !kls_row_refactor_dependency_done_now(shared, dep)) {
    return 0;
  }
  for (UF_long q = l_begin; q < pos; ++q) {
    if (applied[q - l_begin] || (q >= run_begin && q < pos)) {
      continue;
    }
    const UF_long prior = solver->row_refactor_l_cols[q];
    if (kls_row_refactor_u_row_has_col(solver, prior, dep)) {
      return 0;
    }
  }
  return 1;
}

static int kls_row_refactor_prefactor_dep_safe(
  const kls_solver *solver,
  const kls_egraph_refactor_shared *shared,
  UF_long l_begin,
  UF_long pos,
  const unsigned char *applied,
  UF_long dep) {
  return kls_row_refactor_prefactor_dep_safe_with_run(
    solver, shared, l_begin, pos, pos, applied, dep);
}

/* Light strict-supernode run consume: a maximal run of consecutive dep
   rows inside one strict U-row supernode updates x as h short
   triangular steps plus h contiguous AXPYs over the shared tail,
   gathered once into a dense accumulator. Replaces h scattered
   index-walked row updates (the strict-shift pattern guarantees row r
   holds columns r+1..d1 first, then row d1's tail, contiguously in
   u_row_values). Census: 93-99%% of dep entries on the paper-set gap
   rows arrive in such runs. Steady checkless refactors only. */
static int kls_row_refactor_try_light_supernode_run(
  kls_egraph_refactor_worker *worker,
  UF_long row,
  UF_long *p_io,
  UF_long l_begin,
  UF_long l_end,
  int wait_for_dependencies,
  const unsigned char *applied,
  double *x) {
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  const UF_long *sn_end = solver->row_refactor_sn_end;
  if (sn_end == NULL || shared->check_pivots ||
      solver->row_refactor_l_row_values == NULL ||
      solver->row_refactor_l_values == NULL) {
    return 0;
  }
  const UF_long p0 = *p_io;
  const UF_long *l_cols = solver->row_refactor_l_cols;
  const UF_long dep0 = l_cols[p0];
  if (dep0 >= row || dep0 >= solver->n) {
    return 0;
  }
  const UF_long e = sn_end[dep0];
  if (e <= dep0) {
    return 0;
  }
  UF_long q = p0 + 1u;
  while (q < l_end && l_cols[q] == l_cols[q - 1u] + 1u &&
         l_cols[q] <= e &&
         (applied == NULL || !applied[q - l_begin])) {
    q++;
  }
  const UF_long run = q - p0;
  if (run < 2u) {
    return 0;
  }
  const UF_long d1 = dep0 + run - 1u;
  if (wait_for_dependencies) {
    for (UF_long r = dep0; r <= d1; ++r) {
      if (!kls_row_refactor_dependency_done_now(shared, r)) {
        return 0;
      }
    }
  }
  const UF_long *up = solver->row_refactor_u_ptr;
  const UF_long *uc = solver->row_refactor_u_cols;
  const double *uvals = solver->row_refactor_u_row_values;
  const double *udiag = (const double *)solver->numeric->Udiag;
  const UF_long tail_begin = up[d1];
  const UF_long tail_len = up[d1 + 1u] - tail_begin;
  const UF_long *tail_cols = uc + tail_begin;
  double *acc = NULL;
  if (tail_len > 0u) {
    acc = kls_egraph_worker_supernode_workspace(worker, tail_len);
    if (acc == NULL) {
      return 0;
    }
    for (UF_long t = 0; t < tail_len; ++t) {
      acc[t] = x[tail_cols[t]];
    }
  }
  const int scatter_values = !shared->row_refactor_defer_value_scatter;
  for (UF_long r = dep0, p = p0; r <= d1; ++r, ++p) {
    const double lij = x[r] / udiag[r];
    solver->row_refactor_l_row_values[p] = lij;
    if (scatter_values) {
      *solver->row_refactor_l_values[p] = lij;
    }
    x[r] = 0.0;
    const double *uv = uvals + up[r];
    const UF_long tri = d1 - r;
    for (UF_long t = 0; t < tri; ++t) {
      x[uc[up[r] + t]] -= lij * uv[t];
    }
    const double *tv = uv + tri;
    for (UF_long t = 0; t < tail_len; ++t) {
      acc[t] -= lij * tv[t];
    }
  }
  for (UF_long t = 0; t < tail_len; ++t) {
    x[tail_cols[t]] = acc[t];
  }
  *p_io = q;
  return 1;
}

static int kls_parallel_row_refactor_apply_scalar_dep(
  kls_egraph_refactor_worker *worker,
  UF_long row,
  UF_long current_group,
  UF_long p,
  double *x,
  int record_reject) {
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  double *udiag = (double *)solver->numeric->Udiag;
  const UF_long dep = solver->row_refactor_l_cols[p];
  const double candidate = x[dep];
  const double lij = candidate / udiag[dep];
  if (record_reject) {
    if (kls_parallel_row_refactor_rejects_multiplier(
          worker, row, dep, candidate, lij)) {
      x[dep] = 0.0;
      return 0;
    }
  } else if (shared->check_pivots &&
             kls_row_refactor_multiplier_rejects(lij,
                                                 solver->common.tol)) {
    return 2;
  }
  solver->row_refactor_l_row_values[p] = lij;
  if (!shared->row_refactor_defer_value_scatter) {
    *solver->row_refactor_l_values[p] = lij;
  }
  x[dep] = 0.0;
  if (!kls_row_refactor_try_compact_panel_scalar_update(
        shared, x, current_group, dep, lij)) {
    const UF_long u_begin = solver->row_refactor_u_ptr[dep];
    const UF_long u_end = solver->row_refactor_u_ptr[dep + 1u];
    const UF_long *u_cols = solver->row_refactor_u_cols + u_begin;
    const double *u_values = solver->row_refactor_u_row_values + u_begin;
    kls_scatter_subtract(x, u_cols, u_values, u_end - u_begin, lij);
  }
  return 1;
}

static int kls_parallel_row_refactor_prefactor_compact_supernode_update(
  kls_egraph_refactor_worker *worker,
  UF_long row,
  UF_long current_group,
  UF_long external_dep_limit,
  UF_long l_begin,
  UF_long l_end,
  UF_long scan,
  unsigned char *applied,
  double *x,
  UF_long *run_rows_out) {
  KLS_SET_OPTIONAL_OUTPUT(run_rows_out, 0);
  if (worker == NULL || worker->shared == NULL || applied == NULL ||
      x == NULL || run_rows_out == NULL || scan >= l_end) {
    return 0;
  }
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  if (solver == NULL ||
      solver->row_refactor_l_cols == NULL ||
      solver->row_refactor_l_row_values == NULL ||
      solver->row_refactor_l_values == NULL ||
      solver->row_refactor_row_group == NULL ||
      solver->row_refactor_group_ptr == NULL ||
      solver->row_refactor_group_kind == NULL ||
      solver->row_refactor_group_trailing_len == NULL ||
      solver->row_refactor_group_compact_panel_begin == NULL ||
      solver->row_refactor_compact_panel_values == NULL ||
      solver->numeric == NULL || solver->numeric->Udiag == NULL) {
    return 0;
  }
  if (external_dep_limit > solver->n) {
    external_dep_limit = solver->n;
  }

  const UF_long dep0 = solver->row_refactor_l_cols[scan];
  if (dep0 >= external_dep_limit || dep0 >= row || dep0 >= solver->n) {
    return 0;
  }
  const UF_long dep_group = solver->row_refactor_row_group[dep0];
  UF_long valid_prefix =
    kls_row_refactor_compact_panel_valid_prefix(solver, dep_group);
  if (dep_group >= solver->row_refactor_group_count ||
      dep_group == current_group ||
      solver->row_refactor_group_kind[dep_group] !=
        KLS_ROW_REFACTOR_GROUP_DENSE ||
      valid_prefix == 0u) {
    return 0;
  }
  const UF_long group_begin = solver->row_refactor_group_ptr[dep_group];
  const UF_long group_end = solver->row_refactor_group_ptr[dep_group + 1u];
  if (group_begin >= group_end || dep0 < group_begin ||
      dep0 >= group_end || group_end > row) {
    return 0;
  }
  const UF_long width = group_end - group_begin;
  if (valid_prefix > width) {
    valid_prefix = width;
  }
  const UF_long valid_end = group_begin + valid_prefix;
  if (dep0 >= valid_end) {
    return 0;
  }

  UF_long run_end = scan;
  UF_long expected_dep = dep0;
  while (run_end < l_end && expected_dep < valid_end) {
    if (applied[run_end - l_begin] ||
        kls_egraph_refactor_should_stop(shared)) {
      break;
    }
    const UF_long dep = solver->row_refactor_l_cols[run_end];
    if (dep != expected_dep || dep >= external_dep_limit ||
        dep >= row || dep >= solver->n ||
        solver->row_refactor_row_group[dep] != dep_group) {
      break;
    }
    if (!kls_row_refactor_prefactor_dep_safe_with_run(
          solver, shared, l_begin, scan, run_end, applied, dep)) {
      break;
    }
    ++run_end;
    ++expected_dep;
  }

  const UF_long run_rows = run_end - scan;
  if (run_rows < 2u) {
    return 0;
  }
  const UF_long suffix_begin = dep0 - group_begin;
  const UF_long run_prefix_end = suffix_begin + run_rows;
  if (run_prefix_end > width || run_prefix_end > valid_prefix) {
    return 0;
  }

  UF_long dense_panel_entries = 0;
  if (width > (UF_long)(SIZE_MAX / sizeof(double)) / width) {
    kls_egraph_refactor_record_invalid(shared);
    return -1;
  }
  dense_panel_entries = width * width;
  const UF_long trailing_len =
    solver->row_refactor_group_trailing_len[dep_group];
  if (trailing_len > 0u &&
      width > ((UF_long)(SIZE_MAX / sizeof(double)) -
                 dense_panel_entries) / trailing_len) {
    kls_egraph_refactor_record_invalid(shared);
    return -1;
  }
  const UF_long trailing_panel_entries = width * trailing_len;
  const UF_long panel_entries = dense_panel_entries + trailing_panel_entries;
  double *dense_panel =
    kls_row_refactor_compact_group_panel(solver, dep_group, panel_entries);
  if (dense_panel == NULL) {
    return 0;
  }
  double *trailing_panel = trailing_len > 0u
    ? dense_panel + dense_panel_entries : NULL;
  const UF_long *trailing_cols = NULL;
  if (!kls_row_refactor_group_trailing_slice(
        solver, dep_group, group_end, NULL, &trailing_cols)) {
    kls_egraph_refactor_record_invalid(shared);
    return -1;
  }
  if (trailing_len > 0u && trailing_cols == NULL) {
    kls_egraph_refactor_record_invalid(shared);
    return -1;
  }

  const UF_long max_workspace_entries =
    (UF_long)(SIZE_MAX / sizeof(double));
  if (run_rows > max_workspace_entries ||
      trailing_len > max_workspace_entries - run_rows) {
    return 0;
  }
  double *workspace =
    kls_egraph_worker_supernode_workspace(worker, run_rows + trailing_len);
  if (workspace == NULL) {
    return 0;
  }
  double *multipliers = workspace;
  double *trailing_workspace =
    trailing_len > 0u ? workspace + run_rows : NULL;
  const double *udiag = (const double *)solver->numeric->Udiag;

  for (UF_long local = 0; local < run_rows; ++local) {
    const UF_long p = scan + local;
    const UF_long dep = solver->row_refactor_l_cols[p];
    multipliers[local] = x[dep];
  }

  UF_long trsv_entries = 0;
  for (UF_long local = 0; local < run_rows; ++local) {
    const UF_long p = scan + local;
    const UF_long dep = solver->row_refactor_l_cols[p];
    const UF_long local_dep = dep - group_begin;
    const double candidate = multipliers[local];
    const double lij = candidate / udiag[dep];
    if (shared->check_pivots &&
        kls_row_refactor_multiplier_rejects(lij, solver->common.tol)) {
      return 0;
    }
    multipliers[local] = lij;
    const double *dep_dense_panel = dense_panel + local_dep * width;
    for (UF_long target = local_dep + 1u; target < run_prefix_end;
         ++target) {
      multipliers[target - suffix_begin] -=
        lij * dep_dense_panel[target];
    }
    trsv_entries += run_prefix_end - local_dep - 1u;
  }

  if (trailing_len > 0u) {
    memset(trailing_workspace, 0,
           (size_t)trailing_len * sizeof(*trailing_workspace));
  }
  UF_long dense_tail_entries = 0;
  for (UF_long local = 0; local < run_rows; ++local) {
    const UF_long p = scan + local;
    const UF_long dep = solver->row_refactor_l_cols[p];
    const UF_long local_dep = dep - group_begin;
    const double lij = multipliers[local];
    solver->row_refactor_l_row_values[p] = lij;
    if (!shared->row_refactor_defer_value_scatter) {
      *solver->row_refactor_l_values[p] = lij;
    }
    x[dep] = 0.0;

    const double *dep_dense_panel = dense_panel + local_dep * width;
    for (UF_long target = run_prefix_end; target < width; ++target) {
      x[group_begin + target] -= lij * dep_dense_panel[target];
      dense_tail_entries++;
    }
    if (trailing_len > 0u) {
      const double *dep_panel = trailing_panel + local_dep * trailing_len;
      kls_accumulate_scaled_dense(trailing_workspace, dep_panel,
                                  trailing_len, lij);
    }
  }

  UF_long touched_entries = trsv_entries;
  touched_entries = touched_entries > UF_long_max - dense_tail_entries
    ? UF_long_max : touched_entries + dense_tail_entries;
  kls_record_row_refactor_compact_supernode_trsv(
    solver, run_rows, trsv_entries);
  if (trailing_len > 0u) {
    const UF_long gemv_entries =
      run_rows > 0u && trailing_len > UF_long_max / run_rows
        ? UF_long_max : run_rows * trailing_len;
    for (UF_long offset = 0; offset < trailing_len; ++offset) {
      x[trailing_cols[offset]] -= trailing_workspace[offset];
    }
    touched_entries = touched_entries > UF_long_max - gemv_entries
      ? UF_long_max : touched_entries + gemv_entries;
    kls_record_row_refactor_compact_supernode_gemv(
      solver, run_rows, gemv_entries);
  }
  kls_record_row_refactor_compact_supernode_update(
    solver, run_rows, touched_entries);
  if (run_prefix_end < width) {
    kls_record_row_refactor_compact_supernode_partial_update(
      solver, run_rows, touched_entries);
  }
  *run_rows_out = run_rows;
  return 1;
}

static int kls_parallel_row_refactor_prefactor_finished_deps(
  kls_egraph_refactor_worker *worker,
  UF_long row,
  UF_long current_group,
  UF_long external_dep_limit,
  UF_long l_begin,
  UF_long l_end,
  UF_long blocked_pos,
  unsigned char *applied,
  double *x,
  UF_long *prefactor_deps_out,
  UF_long *prefactor_supernode_deps_out) {
  if (worker == NULL || worker->shared == NULL || applied == NULL ||
      prefactor_deps_out == NULL) {
    return 1;
  }
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  if (external_dep_limit > solver->n) {
    external_dep_limit = solver->n;
  }
  for (UF_long scan = blocked_pos + 1u; scan < l_end; ++scan) {
    if (applied[scan - l_begin] ||
        kls_egraph_refactor_should_stop(shared)) {
      continue;
    }
    const UF_long dep = solver->row_refactor_l_cols[scan];
    if (dep >= external_dep_limit) {
      continue;
    }
    if (!kls_row_refactor_prefactor_dep_safe(
          solver, shared, l_begin, scan, applied, dep)) {
      continue;
    }
    UF_long compact_run_rows = 0;
    const int compact_status =
      kls_parallel_row_refactor_prefactor_compact_supernode_update(
        worker, row, current_group, external_dep_limit, l_begin, l_end,
        scan, applied, x, &compact_run_rows);
    if (compact_status < 0) {
      return 0;
    }
    if (compact_status > 0 && compact_run_rows > 0u &&
        compact_run_rows <= l_end - scan) {
      memset(applied + (scan - l_begin), 1, (size_t)compact_run_rows);
      *prefactor_deps_out += compact_run_rows;
      if (prefactor_supernode_deps_out != NULL) {
        *prefactor_supernode_deps_out += compact_run_rows;
      }
      scan += compact_run_rows - 1u;
      continue;
    }
    const int status =
      kls_parallel_row_refactor_apply_scalar_dep(
        worker, row, current_group, scan, x, 0);
    if (status == 0) {
      return 0;
    }
    if (status == 2) {
      continue;
    }
    applied[scan - l_begin] = 1u;
    (*prefactor_deps_out)++;
  }
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

  const int direct_input =
    kls_parallel_row_refactor_load_targeted_input_row(shared, x, row);
  if (direct_input < 0) {
    return 0;
  }
  if (!direct_input &&
      !kls_parallel_row_refactor_load_input_row(shared, x, row)) {
    return 0;
  }

  const UF_long current_group =
    solver->row_refactor_row_group != NULL && row < solver->n
      ? solver->row_refactor_row_group[row] : KLS_KLU_EMPTY;
  const UF_long l_begin = solver->row_refactor_l_ptr[row];
  UF_long p = l_begin;
  const UF_long l_end = solver->row_refactor_l_ptr[row + 1u];
  const UF_long dep_count = l_end - l_begin;
  unsigned char *applied = NULL;
  if (wait_for_dependencies && dep_count > 1u &&
      (uintmax_t)dep_count <= (uintmax_t)SIZE_MAX / sizeof(*applied)) {
    applied = kls_egraph_worker_byte_workspace(worker, dep_count);
    if (applied != NULL) {
      memset(applied, 0, (size_t)dep_count);
    }
  }
  UF_long prefactor_deps = 0;
  UF_long prefactor_supernode_deps = 0;
  while (p < l_end) {
    if (applied != NULL && applied[p - l_begin]) {
      ++p;
      continue;
    }
    const UF_long compact_begin = p;
    const int compact_status =
      kls_row_refactor_try_compact_supernode_update(
        worker, row, current_group, &p, l_begin, l_end,
        wait_for_dependencies, applied);
    if (compact_status < 0) {
      return 0;
    }
    if (compact_status > 0) {
      if (applied != NULL && p > compact_begin && p <= l_end) {
        memset(applied + (compact_begin - l_begin), 1,
               (size_t)(p - compact_begin));
      }
      continue;
    }
    if (kls_row_refactor_try_light_supernode_run(
          worker, row, &p, l_begin, l_end, wait_for_dependencies,
          applied, x)) {
      if (applied != NULL && p > compact_begin) {
        memset(applied + (compact_begin - l_begin), 1,
               (size_t)(p - compact_begin));
      }
      continue;
    }
    const UF_long dep = solver->row_refactor_l_cols[p];
    if (wait_for_dependencies &&
        !kls_row_refactor_dependency_done_now(shared, dep)) {
      if (applied != NULL &&
          !kls_parallel_row_refactor_prefactor_finished_deps(
            worker, row, current_group, row, l_begin, l_end, p, applied,
            x, &prefactor_deps, &prefactor_supernode_deps)) {
        return 0;
      }
      if (!kls_egraph_refactor_wait_done(shared, dep)) {
        x[dep] = 0.0;
        return 0;
      }
    }
    if (!kls_parallel_row_refactor_apply_scalar_dep(
          worker, row, current_group, p, x, 1)) {
      return 0;
    }
    if (applied != NULL) {
      applied[p - l_begin] = 1u;
    }
    p++;
  }

  const UF_long u_begin = solver->row_refactor_u_ptr[row];
  const UF_long u_end = solver->row_refactor_u_ptr[row + 1u];
  double *row_u_values = solver->row_refactor_u_row_values + u_begin;
  const double pivot = direct_input ? udiag[row] + x[row] : x[row];
  if (shared->check_pivots) {
    double row_max_abs = -1.0;
    if (!kls_row_refactor_u_row_abs_max_from_direct(
          solver, row, pivot, x, row_u_values, direct_input,
          &row_max_abs)) {
      kls_egraph_refactor_record_invalid(shared);
      return 0;
    }
    if (kls_parallel_row_refactor_rejects_pivot(worker, row, pivot,
                                               row_max_abs)) {
      memset(x, 0, (size_t)solver->n * sizeof(*x));
      return 0;
    }
  }
  x[row] = 0.0;
  if (pivot == 0.0) {
    kls_egraph_refactor_record_singular(shared, row, symbolic->Q[row]);
    if (solver->common.halt_if_singular) {
      memset(x, 0, (size_t)solver->n * sizeof(*x));
      return 0;
    }
  }
  udiag[row] = pivot;

  if (shared->row_refactor_defer_value_scatter) {
    for (UF_long p = u_begin; p < u_end; ++p) {
      const UF_long offset = p - u_begin;
      const UF_long col = solver->row_refactor_u_cols[p];
      row_u_values[offset] =
        direct_input ? row_u_values[offset] + x[col] : x[col];
      x[col] = 0.0;
    }
  } else {
    double **klu_u_values = solver->row_refactor_u_values + u_begin;
    for (UF_long p = u_begin; p < u_end; ++p) {
      const UF_long offset = p - u_begin;
      const UF_long col = solver->row_refactor_u_cols[p];
      const double value =
        direct_input ? row_u_values[offset] + x[col] : x[col];
      row_u_values[offset] = value;
      *klu_u_values[offset] = value;
      x[col] = 0.0;
    }
  }
  kls_clear_row_refactor_input_residuals(solver, x, row);
  if (prefactor_deps > 0u) {
    atomic_fetch_add_explicit(&shared->row_prefactor_rows, 1ul,
                              memory_order_relaxed);
    atomic_fetch_add_explicit(&shared->row_prefactor_deps,
                              (unsigned long)prefactor_deps,
                              memory_order_relaxed);
  }
  if (prefactor_supernode_deps > 0u) {
    atomic_fetch_add_explicit(&shared->row_prefactor_supernode_rows, 1ul,
                              memory_order_relaxed);
    atomic_fetch_add_explicit(&shared->row_prefactor_supernode_deps,
                              (unsigned long)prefactor_supernode_deps,
                              memory_order_relaxed);
  }
  kls_row_refactor_mark_row_done(worker, row);
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

static double *kls_egraph_worker_supernode_workspace(
  kls_egraph_refactor_worker *worker,
  UF_long entry_count) {
  if (entry_count == 0u) {
    return NULL;
  }
  if (worker == NULL ||
      entry_count >
        (UF_long)(SIZE_MAX / sizeof(*worker->supernode_workspace))) {
    return NULL;
  }
  if (worker->supernode_workspace != NULL &&
      worker->supernode_workspace_size >= entry_count) {
    return worker->supernode_workspace;
  }

  double *workspace =
    (double *)realloc(worker->supernode_workspace,
                      (size_t)entry_count * sizeof(*workspace));
  if (workspace == NULL) {
    return NULL;
  }
  worker->supernode_workspace = workspace;
  worker->supernode_workspace_size = entry_count;
  return workspace;
}

static UF_long *kls_egraph_worker_index_workspace(
  kls_egraph_refactor_worker *worker,
  UF_long entry_count) {
  if (entry_count == 0u) {
    return NULL;
  }
  if (worker == NULL ||
      entry_count > (UF_long)(SIZE_MAX / sizeof(*worker->index_workspace))) {
    return NULL;
  }
  if (worker->index_workspace != NULL &&
      worker->index_workspace_size >= entry_count) {
    return worker->index_workspace;
  }

  UF_long *workspace =
    (UF_long *)realloc(worker->index_workspace,
                       (size_t)entry_count * sizeof(*workspace));
  if (workspace == NULL) {
    return NULL;
  }
  worker->index_workspace = workspace;
  worker->index_workspace_size = entry_count;
  return workspace;
}

static unsigned char *kls_egraph_worker_byte_workspace(
  kls_egraph_refactor_worker *worker,
  UF_long entry_count) {
  if (entry_count == 0u) {
    return NULL;
  }
  if (worker == NULL ||
      (uintmax_t)entry_count >
        (uintmax_t)SIZE_MAX / sizeof(*worker->byte_workspace)) {
    return NULL;
  }
  if (worker->byte_workspace != NULL &&
      worker->byte_workspace_size >= entry_count) {
    return worker->byte_workspace;
  }

  unsigned char *workspace =
    (unsigned char *)realloc(worker->byte_workspace,
                             (size_t)entry_count * sizeof(*workspace));
  if (workspace == NULL) {
    return NULL;
  }
  worker->byte_workspace = workspace;
  worker->byte_workspace_size = entry_count;
  return workspace;
}

static void *kls_egraph_worker_object_workspace(
  kls_egraph_refactor_worker *worker,
  UF_long entry_count,
  size_t entry_size) {
  if (entry_count == 0u || entry_size == 0u) {
    return NULL;
  }
  if ((uintmax_t)entry_count >
      (uintmax_t)UF_long_max / (uintmax_t)entry_size) {
    return NULL;
  }
  const uintmax_t byte_count =
    (uintmax_t)entry_count * (uintmax_t)entry_size;
  if (byte_count > (uintmax_t)SIZE_MAX) {
    return NULL;
  }
  if (worker == NULL || byte_count > (uintmax_t)UF_long_max) {
    return NULL;
  }
  if (worker->object_workspace != NULL &&
      worker->object_workspace_size >= (UF_long)byte_count) {
    return worker->object_workspace;
  }

  void *workspace =
    realloc(worker->object_workspace, (size_t)byte_count);
  if (workspace == NULL) {
    return NULL;
  }
  worker->object_workspace = workspace;
  worker->object_workspace_size = (UF_long)byte_count;
  return workspace;
}

static double *kls_row_refactor_compact_group_panel(kls_solver *solver,
                                                    UF_long group,
                                                    UF_long entry_count) {
  if (solver == NULL ||
      solver->row_refactor_group_compact_panel_begin == NULL ||
      solver->row_refactor_compact_panel_values == NULL ||
      group >= solver->row_refactor_group_count || entry_count == 0u) {
    return NULL;
  }
  const UF_long begin = solver->row_refactor_group_compact_panel_begin[group];
  const UF_long end = solver->row_refactor_group_compact_panel_begin[group + 1u];
  if (end < begin || end - begin != entry_count ||
      end > solver->row_refactor_compact_dense_panel_persistent_entries) {
    return NULL;
  }
  return solver->row_refactor_compact_panel_values + begin;
}

typedef struct {
  UF_long group;
  UF_long len;
} kls_row_refactor_dense_producer_run;

static int kls_row_refactor_find_dense_producer_run_at(
  const kls_solver *solver,
  UF_long row,
  UF_long p0,
  kls_row_refactor_dense_producer_run *run_out) {
  if (run_out != NULL) {
    run_out->group = KLS_KLU_EMPTY;
    run_out->len = 0;
  }
  if (solver == NULL || row >= solver->n ||
      solver->row_refactor_dense_producer_run_ptr == NULL ||
      solver->row_refactor_dense_producer_run_l_begin == NULL ||
      solver->row_refactor_dense_producer_run_group == NULL ||
      solver->row_refactor_dense_producer_run_len == NULL) {
    return 0;
  }
  UF_long lo = solver->row_refactor_dense_producer_run_ptr[row];
  UF_long hi = solver->row_refactor_dense_producer_run_ptr[row + 1u];
  while (lo < hi) {
    const UF_long mid = lo + (hi - lo) / 2u;
    const UF_long run_begin =
      solver->row_refactor_dense_producer_run_l_begin[mid];
    if (run_begin == p0) {
      const UF_long group =
        solver->row_refactor_dense_producer_run_group[mid];
      const UF_long len = solver->row_refactor_dense_producer_run_len[mid];
      if (group >= solver->row_refactor_group_count || len == 0u) {
        return 0;
      }
      if (run_out != NULL) {
        run_out->group = group;
        run_out->len = len;
      }
      return 1;
    }
    if (run_begin < p0) {
      lo = mid + 1u;
    } else {
      hi = mid;
    }
  }
  return 0;
}

static int kls_row_refactor_try_compact_supernode_update(
  kls_egraph_refactor_worker *worker,
  UF_long row,
  UF_long row_group,
  UF_long *p_io,
  UF_long l_begin,
  UF_long l_end,
  int wait_for_dependencies,
  const unsigned char *applied) {
  if (worker == NULL || worker->shared == NULL || p_io == NULL) {
    return -1;
  }
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  if (solver == NULL || *p_io >= l_end || *p_io < l_begin ||
      solver->row_refactor_l_cols == NULL ||
      solver->row_refactor_l_row_values == NULL ||
      solver->row_refactor_l_values == NULL ||
      solver->row_refactor_row_group == NULL ||
      solver->row_refactor_group_ptr == NULL ||
      solver->row_refactor_group_kind == NULL ||
      solver->row_refactor_group_trailing_len == NULL ||
      solver->row_refactor_group_compact_panel_begin == NULL ||
      solver->row_refactor_compact_panel_values == NULL ||
      solver->numeric == NULL || solver->numeric->Udiag == NULL ||
      worker->x == NULL) {
    return 0;
  }

  const UF_long p0 = *p_io;
  const UF_long dep0 = solver->row_refactor_l_cols[p0];
  if (dep0 >= row || dep0 >= solver->n) {
    return 0;
  }
  UF_long dep_group = KLS_KLU_EMPTY;
  UF_long planned_run_len = 0;
  const int has_dense_producer_plan =
    solver->row_refactor_dense_producer_run_ptr != NULL;
  if (has_dense_producer_plan) {
    kls_row_refactor_dense_producer_run planned_run;
    if (!kls_row_refactor_find_dense_producer_run_at(
          solver, row, p0, &planned_run)) {
      return 0;
    }
    dep_group = planned_run.group;
    planned_run_len = planned_run.len;
  } else {
    dep_group = solver->row_refactor_row_group[dep0];
  }
  UF_long valid_prefix =
    kls_row_refactor_compact_panel_valid_prefix(solver, dep_group);
  if (dep_group >= solver->row_refactor_group_count ||
      dep_group == row_group ||
      solver->row_refactor_group_kind[dep_group] !=
        KLS_ROW_REFACTOR_GROUP_DENSE ||
      valid_prefix == 0u) {
    return 0;
  }
  const UF_long group_begin = solver->row_refactor_group_ptr[dep_group];
  const UF_long group_end = solver->row_refactor_group_ptr[dep_group + 1u];
  if (group_begin >= group_end || dep0 < group_begin || dep0 >= group_end ||
      group_end > row) {
    return 0;
  }
  const UF_long width = group_end - group_begin;
  if (valid_prefix > width) {
    valid_prefix = width;
  }
  const UF_long valid_end = group_begin + valid_prefix;
  if (dep0 >= valid_end) {
    return 0;
  }
  const UF_long trailing_len =
    solver->row_refactor_group_trailing_len[dep_group];

  UF_long run_end = p0;
  UF_long run_rows = 0;
  if (has_dense_producer_plan) {
    const UF_long suffix_begin = dep0 - group_begin;
    if (suffix_begin >= valid_prefix) {
      return 0;
    }
    run_rows = valid_prefix - suffix_begin;
    if (run_rows > planned_run_len || p0 > l_end ||
        run_rows > l_end - p0) {
      return 0;
    }
    run_end = p0 + run_rows;
    if (applied != NULL) {
      for (UF_long p = p0; p < run_end; ++p) {
        if (applied[p - l_begin]) {
          return 0;
        }
      }
    }
  } else {
    UF_long expected_dep = dep0;
    while (run_end < l_end && expected_dep < valid_end) {
      if (applied != NULL && applied[run_end - l_begin]) {
        break;
      }
      const UF_long dep = solver->row_refactor_l_cols[run_end];
      if (dep != expected_dep || dep >= row ||
          solver->row_refactor_row_group[dep] != dep_group) {
        break;
      }
      run_end++;
      expected_dep++;
    }
    if (expected_dep != valid_end) {
      return 0;
    }
    run_rows = run_end - p0;
  }
  if (run_rows < 2u) {
    return 0;
  }
  const int partial_update = valid_prefix < width;
  const UF_long suffix_begin = dep0 - group_begin;
  const UF_long run_prefix_end = suffix_begin + run_rows;
  if (run_prefix_end > width || run_prefix_end > valid_prefix) {
    return 0;
  }

  UF_long dense_panel_entries = 0;
  if (width > (UF_long)(SIZE_MAX / sizeof(double)) / width) {
    kls_egraph_refactor_record_invalid(shared);
    return -1;
  }
  dense_panel_entries = width * width;
  if (trailing_len > 0u &&
      width > ((UF_long)(SIZE_MAX / sizeof(double)) -
                 dense_panel_entries) / trailing_len) {
    kls_egraph_refactor_record_invalid(shared);
    return -1;
  }
  const UF_long trailing_panel_entries = width * trailing_len;
  const UF_long panel_entries = dense_panel_entries + trailing_panel_entries;
  double *dense_panel =
    kls_row_refactor_compact_group_panel(solver, dep_group, panel_entries);
  if (dense_panel == NULL) {
    return 0;
  }
  double *trailing_panel = trailing_len > 0u
    ? dense_panel + dense_panel_entries : NULL;
  const UF_long *trailing_cols = NULL;
  if (!kls_row_refactor_group_trailing_slice(
        solver, dep_group, group_end, NULL, &trailing_cols)) {
    kls_egraph_refactor_record_invalid(shared);
    return -1;
  }
  if (trailing_len > 0u && trailing_cols == NULL) {
    kls_egraph_refactor_record_invalid(shared);
    return -1;
  }

  double *x = worker->x;
  const double *udiag = (const double *)solver->numeric->Udiag;
  double *supernode_workspace = NULL;
  double *trailing_multipliers = NULL;
  double *trailing_workspace = NULL;
  int use_trailing_gemv = 0;
  const int want_blocked_trailing =
    shared->row_refactor_blocked_trailing_update &&
    run_rows >= KLS_ROW_REFACTOR_BLOCKED_TRAILING_MIN_ROWS &&
    trailing_len >= KLS_ROW_REFACTOR_BLOCKED_TRAILING_MIN_COLS;
  int use_blocked_trailing = 0;
  const UF_long max_workspace_entries =
    (UF_long)(SIZE_MAX / sizeof(*supernode_workspace));
  if (trailing_len > 0u &&
      run_rows <= max_workspace_entries &&
      trailing_len <= max_workspace_entries -
        (want_blocked_trailing ? run_rows : 0u)) {
    const UF_long multiplier_entries = want_blocked_trailing ? run_rows : 0u;
    const UF_long workspace_entries = multiplier_entries + trailing_len;
    supernode_workspace =
      kls_egraph_worker_supernode_workspace(worker, workspace_entries);
    if (supernode_workspace != NULL) {
      use_blocked_trailing = want_blocked_trailing;
      trailing_multipliers = use_blocked_trailing
        ? supernode_workspace : NULL;
      trailing_workspace = supernode_workspace + multiplier_entries;
      memset(trailing_workspace, 0,
             (size_t)trailing_len * sizeof(*trailing_workspace));
      use_trailing_gemv = 1;
    }
  }
  UF_long touched_entries = 0;
  for (UF_long p = p0; p < run_end; ++p) {
    const UF_long dep = solver->row_refactor_l_cols[p];
    if (wait_for_dependencies &&
        !kls_egraph_refactor_wait_done(shared, dep)) {
      x[dep] = 0.0;
      return -1;
    }
    const UF_long local_dep = dep - group_begin;
    const double candidate = x[dep];
    const double lij = candidate / udiag[dep];
    if (kls_parallel_row_refactor_rejects_multiplier(
          worker, row, dep, candidate, lij)) {
      x[dep] = 0.0;
      return -1;
    }
    solver->row_refactor_l_row_values[p] = lij;
    if (!shared->row_refactor_defer_value_scatter) {
      *solver->row_refactor_l_values[p] = lij;
    }
    x[dep] = 0.0;

    const double *dep_dense_panel = dense_panel + local_dep * width;
    kls_accumulate_scaled_dense(x + group_begin + local_dep + 1u,
                                dep_dense_panel + local_dep + 1u,
                                width - local_dep - 1u, -lij);
    touched_entries += width - local_dep - 1u;
    if (trailing_len > 0u && use_blocked_trailing) {
      trailing_multipliers[p - p0] = lij;
    } else if (trailing_len > 0u && use_trailing_gemv) {
      const double *dep_panel = trailing_panel + local_dep * trailing_len;
      kls_accumulate_scaled_dense(trailing_workspace, dep_panel,
                                  trailing_len, lij);
    } else if (trailing_len > 0u) {
      const double *dep_panel = trailing_panel + local_dep * trailing_len;
      for (UF_long offset = 0; offset < trailing_len; ++offset) {
        x[trailing_cols[offset]] -= lij * dep_panel[offset];
      }
      touched_entries += trailing_len;
    }
  }
  if (use_trailing_gemv) {
    const UF_long gemv_entries = run_rows * trailing_len;
    if (use_blocked_trailing) {
      kls_accumulate_scaled_dense_rows8(
        trailing_workspace,
        trailing_panel + suffix_begin * trailing_len,
        trailing_len, trailing_multipliers, run_rows, trailing_len);
    }
    for (UF_long offset = 0; offset < trailing_len; ++offset) {
      x[trailing_cols[offset]] -= trailing_workspace[offset];
    }
    touched_entries += gemv_entries;
    kls_record_row_refactor_compact_supernode_gemv(
      solver, run_rows, gemv_entries);
  }
  kls_record_row_refactor_compact_supernode_update(solver, run_rows,
                                                   touched_entries);
  if (partial_update) {
    kls_record_row_refactor_compact_supernode_partial_update(
      solver, run_rows, touched_entries);
  }
  *p_io = run_end;
  return 1;
}

static int kls_parallel_row_refactor_load_dense_native_input_row(
  kls_egraph_refactor_shared *shared,
  double *x,
  UF_long row_begin,
  UF_long row_end,
  UF_long trailing_len,
  const UF_long *trailing_cols,
  UF_long row) {
  kls_solver *solver = shared != NULL ? shared->solver : NULL;
  if (solver == NULL || x == NULL || row < row_begin || row >= row_end ||
      row_end > solver->n ||
      (trailing_len > 0u && trailing_cols == NULL) ||
      solver->numeric == NULL || solver->numeric->Udiag == NULL ||
      solver->row_refactor_input_ptr == NULL ||
      solver->row_refactor_input_cols == NULL ||
      solver->row_refactor_input_pos == NULL ||
      solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_l_internal_ptr == NULL ||
      solver->row_refactor_l_row_values == NULL ||
      solver->row_refactor_u_ptr == NULL ||
      solver->row_refactor_u_row_values == NULL) {
    if (shared != NULL) {
      kls_egraph_refactor_record_invalid(shared);
    }
    return -1;
  }

  const UF_long l_internal = solver->row_refactor_l_internal_ptr[row];
  const UF_long l_end = solver->row_refactor_l_ptr[row + 1u];
  if (l_internal > l_end || l_end - l_internal != row - row_begin) {
    kls_egraph_refactor_record_invalid(shared);
    return -1;
  }
  double *row_l_values = solver->row_refactor_l_row_values + l_internal;
  const UF_long input_begin = solver->row_refactor_input_ptr[row];
  const UF_long input_end = solver->row_refactor_input_ptr[row + 1u];
  const int use_retained_targets =
    kls_row_refactor_segment_input_targets_available(
      solver, row, input_begin, input_end);
  const UF_long dense_len = row_end - row - 1u;
  const UF_long u_begin = solver->row_refactor_u_ptr[row];
  const UF_long u_end = solver->row_refactor_u_ptr[row + 1u];
  if (u_end < u_begin || u_end - u_begin != dense_len + trailing_len) {
    kls_egraph_refactor_record_invalid(shared);
    return -1;
  }
  double *row_u_values = solver->row_refactor_u_row_values + u_begin;
  double *udiag = (double *)solver->numeric->Udiag;

  if (!use_retained_targets) {
    for (UF_long p = input_begin; p < input_end; ++p) {
      const UF_long col = solver->row_refactor_input_cols[p];
      if (col >= solver->n) {
        kls_egraph_refactor_record_invalid(shared);
        return -1;
      }
      if (col < row_begin || col < row_end) {
        continue;
      }
      UF_long trailing_pos = KLS_KLU_EMPTY;
      if (trailing_len == 0u ||
          !kls_find_uflong_sorted(trailing_cols, trailing_len, col,
                                  &trailing_pos)) {
        return 0;
      }
    }
  }

  if (row > row_begin) {
    memset(row_l_values, 0,
           (size_t)(row - row_begin) * sizeof(*row_l_values));
  }
  if (dense_len + trailing_len > 0u) {
    memset(row_u_values, 0,
           (size_t)(dense_len + trailing_len) * sizeof(*row_u_values));
  }
  udiag[row] = 0.0;

  double input_multiplier = 1.0;
  if (!kls_row_refactor_input_multiplier(shared, row, &input_multiplier)) {
    kls_egraph_refactor_record_invalid(shared);
    return -1;
  }
  for (UF_long p = input_begin; p < input_end; ++p) {
    const UF_long col = solver->row_refactor_input_cols[p];
    double value = 0.0;
    if (!kls_row_refactor_scaled_input_value(
          shared, kls_row_refactor_input_pos_at(solver, p), input_multiplier,
          &value)) {
      kls_egraph_refactor_record_invalid(shared);
      return -1;
    }
    if (use_retained_targets) {
      if (!kls_row_refactor_apply_segment_input_target(
            shared, x, row, p, l_internal, l_end,
            dense_len + trailing_len, row_u_values, udiag, value)) {
        return -1;
      }
    } else {
      if (col < row_begin) {
        x[col] = value;
      } else if (col < row) {
        row_l_values[col - row_begin] = value;
      } else if (col == row) {
        udiag[row] = value;
      } else if (col < row_end) {
        row_u_values[col - row - 1u] = value;
      } else {
        UF_long trailing_pos = KLS_KLU_EMPTY;
        if (!kls_find_uflong_sorted(trailing_cols, trailing_len, col,
                                    &trailing_pos) ||
            trailing_pos >= trailing_len) {
          kls_egraph_refactor_record_invalid(shared);
          return -1;
        }
        row_u_values[dense_len + trailing_pos] = value;
      }
    }
  }
  if (use_retained_targets) {
    kls_record_row_refactor_segment_target_input(solver);
  }
  kls_record_row_refactor_dense_segment_direct_input(solver);
  return 1;
}

static int kls_row_refactor_segment_input_targets_available(
  const kls_solver *solver,
  UF_long row,
  UF_long input_begin,
  UF_long input_end) {
  return solver != NULL && row < solver->n &&
         solver->row_refactor_segment_input_row_ready != NULL &&
         solver->row_refactor_segment_input_row_ready[row] &&
         (input_begin == input_end ||
          (solver->row_refactor_segment_input_target_kind != NULL &&
           solver->row_refactor_segment_input_target_pos != NULL));
}

static int kls_row_refactor_apply_segment_input_target(
  kls_egraph_refactor_shared *shared,
  double *x,
  UF_long row,
  UF_long input_pos,
  UF_long l_begin,
  UF_long l_end,
  UF_long row_u_len,
  double *row_u_values,
  double *udiag,
  double value) {
  kls_solver *solver = shared != NULL ? shared->solver : NULL;
  if (solver == NULL || x == NULL || row >= solver->n ||
      solver->row_refactor_segment_input_target_kind == NULL ||
      solver->row_refactor_segment_input_target_pos == NULL ||
      solver->row_refactor_l_row_values == NULL ||
      row_u_values == NULL || udiag == NULL) {
    if (shared != NULL) {
      kls_egraph_refactor_record_invalid(shared);
    }
    return 0;
  }
  const unsigned char target =
    solver->row_refactor_segment_input_target_kind[input_pos];
  const UF_long pos =
    solver->row_refactor_segment_input_target_pos[input_pos];
  switch ((kls_row_segment_input_target_kind)target) {
    case KLS_ROW_SEGMENT_INPUT_TARGET_EXTERNAL:
      if (pos >= solver->n) {
        kls_egraph_refactor_record_invalid(shared);
        return 0;
      }
      x[pos] = value;
      return 1;
    case KLS_ROW_SEGMENT_INPUT_TARGET_L:
      if (pos < l_begin || pos >= l_end) {
        kls_egraph_refactor_record_invalid(shared);
        return 0;
      }
      solver->row_refactor_l_row_values[pos] = value;
      return 1;
    case KLS_ROW_SEGMENT_INPUT_TARGET_PIVOT:
      udiag[row] = value;
      return 1;
    case KLS_ROW_SEGMENT_INPUT_TARGET_U:
      if (pos >= row_u_len) {
        kls_egraph_refactor_record_invalid(shared);
        return 0;
      }
      row_u_values[pos] = value;
      return 1;
    case KLS_ROW_SEGMENT_INPUT_TARGET_NONE:
    default:
      kls_egraph_refactor_record_invalid(shared);
      return 0;
  }
}

static int kls_parallel_row_refactor_load_targeted_input_row(
  kls_egraph_refactor_shared *shared,
  double *x,
  UF_long row) {
  kls_solver *solver = shared != NULL ? shared->solver : NULL;
  if (solver == NULL || x == NULL || row >= solver->n ||
      solver->numeric == NULL || solver->numeric->Udiag == NULL ||
      solver->row_refactor_input_ptr == NULL ||
      solver->row_refactor_input_cols == NULL ||
      solver->row_refactor_input_pos == NULL ||
      solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_u_ptr == NULL ||
      solver->row_refactor_u_row_values == NULL ||
      solver->row_refactor_segment_input_row_ready == NULL ||
      !solver->row_refactor_segment_input_row_ready[row]) {
    return 0;
  }
  const UF_long input_begin = solver->row_refactor_input_ptr[row];
  const UF_long input_end = solver->row_refactor_input_ptr[row + 1u];
  if (!kls_row_refactor_segment_input_targets_available(
        solver, row, input_begin, input_end)) {
    return 0;
  }
  if (solver->row_refactor_row_group != NULL &&
      solver->row_refactor_group_kind != NULL) {
    const UF_long group = solver->row_refactor_row_group[row];
    if (group < solver->row_refactor_group_count &&
        solver->row_refactor_group_kind[group] ==
          KLS_ROW_REFACTOR_GROUP_BATCH) {
      return 0;
    }
  }

  const UF_long l_begin = solver->row_refactor_l_ptr[row];
  const UF_long l_end = solver->row_refactor_l_ptr[row + 1u];
  const UF_long u_begin = solver->row_refactor_u_ptr[row];
  const UF_long u_end = solver->row_refactor_u_ptr[row + 1u];
  if (l_begin > l_end || u_begin > u_end) {
    kls_egraph_refactor_record_invalid(shared);
    return -1;
  }
  double *row_u_values = solver->row_refactor_u_row_values + u_begin;
  const UF_long row_u_len = u_end - u_begin;
  if (row_u_len > 0u) {
    memset(row_u_values, 0, (size_t)row_u_len * sizeof(*row_u_values));
  }
  double *udiag = (double *)solver->numeric->Udiag;
  udiag[row] = 0.0;

  double input_multiplier = 1.0;
  if (!kls_row_refactor_input_multiplier(shared, row, &input_multiplier)) {
    kls_egraph_refactor_record_invalid(shared);
    return -1;
  }
  for (UF_long p = input_begin; p < input_end; ++p) {
    double value = 0.0;
    if (!kls_row_refactor_scaled_input_value(
          shared, kls_row_refactor_input_pos_at(solver, p), input_multiplier,
          &value)) {
      kls_egraph_refactor_record_invalid(shared);
      return -1;
    }
    if (!kls_row_refactor_apply_segment_input_target(
          shared, x, row, p, l_begin, l_end, row_u_len,
          row_u_values, udiag, value)) {
      return -1;
    }
  }
  kls_record_row_refactor_segment_target_input(solver);
  return 1;
}

static int kls_parallel_row_refactor_load_sparse_segment_input_row(
  kls_egraph_refactor_shared *shared,
  double *x,
  UF_long row_begin,
  UF_long row_end,
  UF_long trailing_len,
  const UF_long *trailing_cols,
  UF_long row) {
  kls_solver *solver = shared != NULL ? shared->solver : NULL;
  if (solver == NULL || x == NULL || row < row_begin || row >= row_end ||
      row_end > solver->n ||
      (trailing_len > 0u && trailing_cols == NULL) ||
      solver->numeric == NULL || solver->numeric->Udiag == NULL ||
      solver->row_refactor_input_ptr == NULL ||
      solver->row_refactor_input_cols == NULL ||
      solver->row_refactor_input_pos == NULL ||
      solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_l_cols == NULL ||
      solver->row_refactor_l_row_values == NULL ||
      solver->row_refactor_u_ptr == NULL ||
      solver->row_refactor_u_cols == NULL ||
      solver->row_refactor_u_row_values == NULL) {
    if (shared != NULL) {
      kls_egraph_refactor_record_invalid(shared);
    }
    return -1;
  }

  const UF_long l_begin = solver->row_refactor_l_ptr[row];
  const UF_long l_end = solver->row_refactor_l_ptr[row + 1u];
  if (l_begin > l_end) {
    kls_egraph_refactor_record_invalid(shared);
    return -1;
  }
  const UF_long input_begin = solver->row_refactor_input_ptr[row];
  const UF_long input_end = solver->row_refactor_input_ptr[row + 1u];
  const int use_retained_targets =
    kls_row_refactor_segment_input_targets_available(
      solver, row, input_begin, input_end);
  const UF_long dense_len = row_end - row - 1u;
  const UF_long u_begin = solver->row_refactor_u_ptr[row];
  const UF_long u_end = solver->row_refactor_u_ptr[row + 1u];
  if (u_begin > u_end || u_end - u_begin != dense_len + trailing_len) {
    kls_egraph_refactor_record_invalid(shared);
    return -1;
  }
  for (UF_long offset = 0; offset < dense_len; ++offset) {
    if (solver->row_refactor_u_cols[u_begin + offset] != row + 1u + offset) {
      kls_egraph_refactor_record_invalid(shared);
      return -1;
    }
  }
  for (UF_long offset = 0; offset < trailing_len; ++offset) {
    if (solver->row_refactor_u_cols[u_begin + dense_len + offset] !=
        trailing_cols[offset]) {
      kls_egraph_refactor_record_invalid(shared);
      return -1;
    }
  }

  if (!use_retained_targets) {
    for (UF_long p = input_begin; p < input_end; ++p) {
      const UF_long col = solver->row_refactor_input_cols[p];
      if (col >= solver->n) {
        kls_egraph_refactor_record_invalid(shared);
        return -1;
      }
      if (col >= row_begin && col < row) {
        UF_long l_pos = KLS_KLU_EMPTY;
        if (!kls_find_uflong_sorted(solver->row_refactor_l_cols + l_begin,
                                    l_end - l_begin, col, &l_pos) ||
            l_pos >= l_end - l_begin) {
          return 0;
        }
      } else if (col >= row_end) {
        UF_long trailing_pos = KLS_KLU_EMPTY;
        if (trailing_len == 0u ||
            !kls_find_uflong_sorted(trailing_cols, trailing_len, col,
                                    &trailing_pos)) {
          return 0;
        }
      }
    }
  }

  for (UF_long p = l_begin; p < l_end; ++p) {
    const UF_long dep = solver->row_refactor_l_cols[p];
    if (dep >= row) {
      kls_egraph_refactor_record_invalid(shared);
      return -1;
    }
    if (dep >= row_begin) {
      solver->row_refactor_l_row_values[p] = 0.0;
    }
  }
  double *row_u_values = solver->row_refactor_u_row_values + u_begin;
  if (dense_len + trailing_len > 0u) {
    memset(row_u_values, 0,
           (size_t)(dense_len + trailing_len) * sizeof(*row_u_values));
  }
  double *udiag = (double *)solver->numeric->Udiag;
  udiag[row] = 0.0;

  double input_multiplier = 1.0;
  if (!kls_row_refactor_input_multiplier(shared, row, &input_multiplier)) {
    kls_egraph_refactor_record_invalid(shared);
    return -1;
  }
  for (UF_long p = input_begin; p < input_end; ++p) {
    const UF_long col = solver->row_refactor_input_cols[p];
    double value = 0.0;
    if (!kls_row_refactor_scaled_input_value(
          shared, kls_row_refactor_input_pos_at(solver, p), input_multiplier,
          &value)) {
      kls_egraph_refactor_record_invalid(shared);
      return -1;
    }
    if (use_retained_targets) {
      if (!kls_row_refactor_apply_segment_input_target(
            shared, x, row, p, l_begin, l_end,
            dense_len + trailing_len, row_u_values, udiag, value)) {
        return -1;
      }
    } else {
      if (col < row_begin) {
        x[col] = value;
      } else if (col < row) {
        UF_long l_pos = KLS_KLU_EMPTY;
        if (!kls_find_uflong_sorted(solver->row_refactor_l_cols + l_begin,
                                    l_end - l_begin, col, &l_pos) ||
            l_pos >= l_end - l_begin) {
          kls_egraph_refactor_record_invalid(shared);
          return -1;
        }
        solver->row_refactor_l_row_values[l_begin + l_pos] = value;
      } else if (col == row) {
        udiag[row] = value;
      } else if (col < row_end) {
        row_u_values[col - row - 1u] = value;
      } else {
        UF_long trailing_pos = KLS_KLU_EMPTY;
        if (!kls_find_uflong_sorted(trailing_cols, trailing_len, col,
                                    &trailing_pos) ||
            trailing_pos >= trailing_len) {
          kls_egraph_refactor_record_invalid(shared);
          return -1;
        }
        row_u_values[dense_len + trailing_pos] = value;
      }
    }
  }
  if (use_retained_targets) {
    kls_record_row_refactor_segment_target_input(solver);
  }
  kls_record_row_refactor_sparse_segment_direct_input(solver);
  return 1;
}

static int kls_row_refactor_segment_u_abs_max(
  UF_long row,
  UF_long row_end,
  UF_long trailing_len,
  const UF_long *trailing_cols,
  const double *x,
  const double *row_u_values,
  int direct_input,
  double pivot,
  double *row_max_abs_out) {
  if (row_max_abs_out != NULL) {
    *row_max_abs_out = -1.0;
  }
  if (row >= row_end || x == NULL || row_max_abs_out == NULL ||
      (direct_input && row_u_values == NULL) ||
      (trailing_len > 0u && trailing_cols == NULL)) {
    return 0;
  }

  double row_max_abs = fabs(pivot);
  if (!isfinite(row_max_abs)) {
    *row_max_abs_out = DBL_MAX;
    return 1;
  }
  const UF_long dense_len = row_end - row - 1u;
  for (UF_long offset = 0; offset < dense_len; ++offset) {
    const UF_long col = row + 1u + offset;
    const double value =
      (direct_input ? row_u_values[offset] : 0.0) + x[col];
    const double value_abs = fabs(value);
    if (!isfinite(value_abs)) {
      *row_max_abs_out = DBL_MAX;
      return 1;
    }
    if (value_abs > row_max_abs) {
      row_max_abs = value_abs;
    }
  }
  for (UF_long offset = 0; offset < trailing_len; ++offset) {
    const UF_long col = trailing_cols[offset];
    const double value =
      (direct_input ? row_u_values[dense_len + offset] : 0.0) + x[col];
    const double value_abs = fabs(value);
    if (!isfinite(value_abs)) {
      *row_max_abs_out = DBL_MAX;
      return 1;
    }
    if (value_abs > row_max_abs) {
      row_max_abs = value_abs;
    }
  }
  *row_max_abs_out = row_max_abs;
  return 1;
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
  if (solver->row_refactor_group_dense == NULL ||
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
  const UF_long width = row_end - row_begin;
  UF_long *row_prefactor_dep_counts = NULL;
  UF_long *row_prefactor_supernode_dep_counts = NULL;
  if (wait_for_dependencies && width > 0u) {
    row_prefactor_dep_counts =
      width <= (UF_long)(SIZE_MAX / (2u * sizeof(*row_prefactor_dep_counts)))
        ? kls_egraph_worker_index_workspace(worker, 2u * width) : NULL;
    if (row_prefactor_dep_counts != NULL) {
      row_prefactor_supernode_dep_counts =
        row_prefactor_dep_counts + width;
      memset(row_prefactor_dep_counts, 0,
             (size_t)(2u * width) * sizeof(*row_prefactor_dep_counts));
    }
  }

  for (UF_long row = row_begin; row < row_end; ++row) {
    const int direct_input =
      kls_parallel_row_refactor_load_dense_native_input_row(
        shared, x, row_begin, row_end, trailing_len, trailing_cols, row);
    if (direct_input < 0) {
      return 0;
    }
    if (!direct_input &&
        !kls_parallel_row_refactor_load_input_row(shared, x, row)) {
      return 0;
    }

    const UF_long l_begin = solver->row_refactor_l_ptr[row];
    UF_long lp = l_begin;
    const UF_long l_internal = solver->row_refactor_l_internal_ptr[row];
    unsigned char *applied = NULL;
    if (wait_for_dependencies && l_internal > l_begin + 1u &&
        (uintmax_t)(l_internal - l_begin) <=
          (uintmax_t)SIZE_MAX / sizeof(*applied)) {
      applied = kls_egraph_worker_byte_workspace(worker,
                                                 l_internal - l_begin);
      if (applied != NULL) {
        memset(applied, 0, (size_t)(l_internal - l_begin));
      }
    }
    UF_long prefactor_deps = 0;
    UF_long prefactor_supernode_deps = 0;
    while (lp < l_internal) {
      if (applied != NULL && applied[lp - l_begin]) {
        ++lp;
        continue;
      }
      const UF_long compact_begin = lp;
      const int compact_status =
        kls_row_refactor_try_compact_supernode_update(
          worker, row, group, &lp, l_begin, l_internal,
          wait_for_dependencies, applied);
      if (compact_status < 0) {
        x[row] = 0.0;
        return 0;
      }
      if (compact_status > 0) {
        if (applied != NULL) {
          memset(applied + (compact_begin - l_begin), 1,
                 (size_t)(lp - compact_begin));
        }
        continue;
      }
      const UF_long dep = solver->row_refactor_l_cols[lp];
      if (wait_for_dependencies &&
          !kls_row_refactor_dependency_done_now(shared, dep)) {
        if (applied != NULL &&
            !kls_parallel_row_refactor_prefactor_finished_deps(
              worker, row, group, row_begin, l_begin, l_internal, lp,
              applied, x, &prefactor_deps, &prefactor_supernode_deps)) {
          x[dep] = 0.0;
          return 0;
        }
        if (!kls_egraph_refactor_wait_done(shared, dep)) {
          x[dep] = 0.0;
          return 0;
        }
      }
      const int scalar_status =
        kls_parallel_row_refactor_apply_scalar_dep(
          worker, row, group, lp, x, 1);
      if (scalar_status != 1) {
        x[dep] = 0.0;
        return 0;
      }
      if (applied != NULL) {
        applied[lp - l_begin] = 1u;
      }
      lp++;
    }

    const UF_long l_end = solver->row_refactor_l_ptr[row + 1u];
    double *row_l_values = solver->row_refactor_l_row_values + l_internal;
    for (UF_long dep = row_begin; dep < row; ++dep, ++lp) {
      if (lp >= l_end || solver->row_refactor_l_cols[lp] != dep) {
        kls_egraph_refactor_record_invalid(shared);
        return 0;
      }
      if (direct_input) {
        row_l_values[dep - row_begin] += x[dep];
      } else {
        row_l_values[dep - row_begin] = x[dep];
      }
      x[dep] = 0.0;
    }
    if (lp != l_end) {
      kls_egraph_refactor_record_invalid(shared);
      return 0;
    }

    if (direct_input) {
      udiag[row] += x[row];
    } else {
      udiag[row] = x[row];
    }
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
      if (direct_input) {
        row_u_values[offset] += x[col];
      } else {
        row_u_values[offset] = x[col];
      }
      x[col] = 0.0;
    }
    for (UF_long offset = 0; offset < trailing_len; ++offset) {
      const UF_long col = trailing_cols[offset];
      if (direct_input) {
        row_u_values[dense_len + offset] += x[col];
      } else {
        row_u_values[dense_len + offset] = x[col];
      }
      x[col] = 0.0;
    }
    kls_clear_row_refactor_input_residuals(solver, x, row);
    if (row_prefactor_dep_counts != NULL) {
      row_prefactor_dep_counts[row - row_begin] = prefactor_deps;
      row_prefactor_supernode_dep_counts[row - row_begin] =
        prefactor_supernode_deps;
    }
  }

  for (UF_long row = row_begin; row < row_end; ++row) {
    const UF_long l_internal = solver->row_refactor_l_internal_ptr[row];
    double *row_l_values = solver->row_refactor_l_row_values + l_internal;
    const UF_long row_dense_len = row_end - row - 1u;
    const UF_long row_u_begin = solver->row_refactor_u_ptr[row];
    double *row_u_values =
      solver->row_refactor_u_row_values + row_u_begin;

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

      for (UF_long target = dep + 1u; target < row; ++target) {
        row_l_values[target - row_begin] -=
          lij * dep_u_values[target - dep - 1u];
      }
      udiag[row] -= lij * dep_u_values[row - dep - 1u];
      /* The current row stores its dense U suffix immediately before its
         trailing values.  The completed producer has the corresponding
         suffix in the same order, so both former loops are one contiguous
         AXPY and need only one target-cloned kernel dispatch. */
      const UF_long output_len = row_dense_len + trailing_len;
      const double *dep_output_values = dep_u_values + row - dep;
      kls_accumulate_scaled_dense(row_u_values, dep_output_values,
                                  output_len, -lij);
    }

    const double pivot = udiag[row];
    if (shared->check_pivots) {
      const double row_max_abs =
        kls_row_refactor_abs_max_from_values(
          pivot, row_u_values, row_dense_len + trailing_len);
      if (kls_parallel_row_refactor_rejects_pivot(worker, row, pivot,
                                                 row_max_abs)) {
        memset(x, 0, (size_t)solver->n * sizeof(*x));
        return 0;
      }
    }
    if (pivot == 0.0) {
      kls_egraph_refactor_record_singular(shared, row, symbolic->Q[row]);
      if (solver->common.halt_if_singular) {
        memset(x, 0, (size_t)solver->n * sizeof(*x));
        return 0;
      }
    }
    if (!shared->row_refactor_defer_value_scatter &&
        !kls_scatter_row_refactor_row_values(solver, row)) {
      kls_egraph_refactor_record_invalid(shared);
      return 0;
    }
    const UF_long row_prefactor_deps =
      row_prefactor_dep_counts != NULL
        ? row_prefactor_dep_counts[row - row_begin] : 0u;
    if (row_prefactor_deps > 0u) {
      atomic_fetch_add_explicit(&shared->row_prefactor_rows, 1ul,
                                memory_order_relaxed);
      atomic_fetch_add_explicit(&shared->row_prefactor_deps,
                                (unsigned long)row_prefactor_deps,
                                memory_order_relaxed);
    }
    const UF_long row_prefactor_supernode_deps =
      row_prefactor_supernode_dep_counts != NULL
        ? row_prefactor_supernode_dep_counts[row - row_begin] : 0u;
    if (row_prefactor_supernode_deps > 0u) {
      atomic_fetch_add_explicit(&shared->row_prefactor_supernode_rows, 1ul,
                                memory_order_relaxed);
      atomic_fetch_add_explicit(&shared->row_prefactor_supernode_deps,
                                (unsigned long)row_prefactor_supernode_deps,
                                memory_order_relaxed);
    }
    kls_row_refactor_mark_row_done(worker, row);
  }
  return 1;
}

static double kls_row_refactor_dense_group_update_work(UF_long width,
                                                       UF_long trailing_len) {
  const double w = (double)width;
  const double trailing = (double)trailing_len;
  return (w * w * w / 3.0) + (0.5 * w * w * trailing);
}

static double kls_row_refactor_dense_group_panel_entries(
  UF_long width,
  UF_long trailing_len) {
  const double w = (double)width;
  const double trailing = (double)trailing_len;
  return (w * w) + (w * trailing);
}

static int kls_row_refactor_prefers_compact_dense_panel(UF_long width,
                                                        UF_long trailing_len) {
  const double update_work =
    kls_row_refactor_dense_group_update_work(width, trailing_len);
  const double panel_entries =
    kls_row_refactor_dense_group_panel_entries(width, trailing_len);
  return update_work >= KLS_ROW_REFACTOR_COMPACT_PANEL_MIN_WORK &&
         update_work >=
           KLS_ROW_REFACTOR_COMPACT_PANEL_MIN_WORK_PER_ENTRY * panel_entries;
}

static int kls_native_row_panel_enabled_for_run(const kls_solver *solver,
                                                int env_state,
                                                int check_pivots) {
  if (env_state < 0) {
    return 0;
  }
  if (env_state == 0 && !check_pivots && solver != NULL &&
      solver->row_refactor_native_row_panel_auto_disabled) {
    return 0;
  }
  return 1;
}

static int kls_native_row_panel_should_try_for_solver(
  const kls_solver *solver,
  UF_long width,
  int env_state,
  int check_pivots) {
  if (!kls_native_row_panel_enabled_for_run(solver, env_state,
                                            check_pivots)) {
    return 0;
  }
  return width > 1u;
}

static void kls_maybe_disable_native_row_panel_auto(
  kls_solver *solver,
  int check_pivots,
  int env_state,
  UF_long native_count_before,
  UF_long compact_update_count_before,
  UF_long partial_pipeline_count_before) {
  if (solver == NULL || check_pivots || env_state != 0 ||
      solver->row_refactor_native_row_panel_auto_disabled ||
      solver->row_refactor_native_row_panel_count <= native_count_before ||
      solver->row_refactor_compact_supernode_update_count >
        compact_update_count_before ||
      solver->row_refactor_partial_supernode_pipeline_run_count >
        partial_pipeline_count_before) {
    return;
  }
  solver->row_refactor_native_row_panel_auto_disabled = 1;
  solver->row_refactor_native_row_panel_auto_disable_count++;
  solver->row_refactor_native_row_panel_enabled = 0;
}

static int kls_store_compact_dense_panel_row(
  kls_egraph_refactor_shared *shared,
  UF_long row_begin,
  UF_long row_end,
  UF_long trailing_len,
  double *dense_panel,
  double *trailing_panel,
  UF_long row,
  int scatter_values) {
  kls_solver *solver = shared != NULL ? shared->solver : NULL;
  if (solver == NULL || row < row_begin || row >= row_end ||
      solver->row_refactor_l_internal_ptr == NULL ||
      solver->row_refactor_l_row_values == NULL ||
      solver->row_refactor_l_values == NULL ||
      solver->row_refactor_u_ptr == NULL ||
      solver->row_refactor_u_row_values == NULL ||
      solver->row_refactor_u_values == NULL ||
      dense_panel == NULL ||
      (trailing_len > 0u && trailing_panel == NULL)) {
    return 0;
  }

  const UF_long width = row_end - row_begin;
  const UF_long local_row = row - row_begin;
  double *row_dense_panel = dense_panel + local_row * width;
  double *row_panel = trailing_len > 0u
    ? trailing_panel + local_row * trailing_len : NULL;
  const UF_long row_dense_len = row_end - row - 1u;

  const UF_long l_dense_begin =
    solver->row_refactor_l_internal_ptr[row];
  double *row_l_row_values =
    solver->row_refactor_l_row_values + l_dense_begin;
  double **row_l_values = solver->row_refactor_l_values + l_dense_begin;
  for (UF_long dep = row_begin; dep < row; ++dep) {
    const double value = row_dense_panel[dep - row_begin];
    row_l_row_values[dep - row_begin] = value;
    if (scatter_values && !shared->row_refactor_defer_value_scatter) {
      *row_l_values[dep - row_begin] = value;
    }
  }

  const UF_long row_u_begin = solver->row_refactor_u_ptr[row];
  double *row_u_row_values =
    solver->row_refactor_u_row_values + row_u_begin;
  if (!scatter_values || shared->row_refactor_defer_value_scatter) {
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
  return 1;
}

static UF_long kls_work_to_uflong(double work);

enum kls_dense_help_kind {
  KLS_DENSE_HELP_ASSEMBLE = 0,
  KLS_DENSE_HELP_PANEL_STEP = 1,
  KLS_DENSE_HELP_PANEL_BLOCK_TAIL = 2
};

struct kls_dense_help_ctx {
  long epoch;   /* written last before publish; claims bind to it */
  int kind;
  UF_long group;
  UF_long row_begin;
  UF_long row_end;
  UF_long trailing_len;
  const UF_long *trailing_cols;
  double *dense_panel;
  double *trailing_panel;
  int wait_for_dependencies;
  /* KLS_DENSE_HELP_PANEL_STEP: current block columns [k, kb) of the
     width*width dense panel; cursor claims local panel rows in
     [kb, width). */
  UF_long panel_k;
  UF_long panel_kb;
};

#define KLS_DENSE_HELP_TAIL_SLICE 128u

/* Diagonal-block row tails of a blocked dense-panel factorization form
   a forward substitution: row i's tail update reads the FINAL tails of
   rows k..i-1, so rows are sequential — but every tail column is
   independent. Each participant runs the whole row-sequential
   substitution on a disjoint column slice; the element-wise operation
   order matches the serial code exactly. Offsets [0, width-kb) map to
   the dense targets, the rest to the trailing panel. */
static void kls_compact_dense_panel_block_tail_slice(
  double *dense_panel,
  double *trailing_panel,
  UF_long width,
  UF_long trailing_len,
  UF_long k,
  UF_long kb,
  UF_long off_begin,
  UF_long off_end) {
  const UF_long dense_tail = width - kb;
  const UF_long dense_end = off_end < dense_tail ? off_end : dense_tail;
  const UF_long trail_begin =
    off_begin > dense_tail ? off_begin - dense_tail : 0u;
  const UF_long trail_end = off_end > dense_tail ? off_end - dense_tail : 0u;
  for (UF_long i = k + 1u; i < kb; ++i) {
    double *row_dense_panel = dense_panel + i * width;
    double *row_panel = trailing_len > 0u
      ? trailing_panel + i * trailing_len : NULL;
    for (UF_long dep = k; dep < i; ++dep) {
      const double lij = row_dense_panel[dep];
      const double *dep_dense_panel = dense_panel + dep * width;
      for (UF_long off = off_begin; off < dense_end; ++off) {
        row_dense_panel[kb + off] -= lij * dep_dense_panel[kb + off];
      }
      if (trail_end > trail_begin) {
        const double *dep_panel = trailing_panel + dep * trailing_len;
        for (UF_long t = trail_begin; t < trail_end; ++t) {
          row_panel[t] -= lij * dep_panel[t];
        }
      }
    }
  }
}

/* One row of a blocked dense-panel factorization step: solve the row's
   L values against block columns [k, kb) (rows above kb are final),
   then apply the block's outer-product update to the row's remaining
   dense and trailing targets. Rows are independent within a step. */
static void kls_compact_dense_panel_block_step_row(
  double *dense_panel,
  double *trailing_panel,
  const double *udiag_base,
  UF_long width,
  UF_long trailing_len,
  UF_long k,
  UF_long kb,
  UF_long row) {
  double *row_dense_panel = dense_panel + row * width;
  for (UF_long dep = k; dep < kb; ++dep) {
    double value = row_dense_panel[dep];
    for (UF_long prev = k; prev < dep; ++prev) {
      value -= row_dense_panel[prev] * dense_panel[prev * width + dep];
    }
    row_dense_panel[dep] = value / udiag_base[dep];
  }
  for (UF_long target = kb; target < width; ++target) {
    double update = 0.0;
    for (UF_long dep = k; dep < kb; ++dep) {
      update += row_dense_panel[dep] * dense_panel[dep * width + target];
    }
    row_dense_panel[target] -= update;
  }
  if (trailing_len > 0u) {
    double *row_panel = trailing_panel + row * trailing_len;
    for (UF_long offset = 0; offset < trailing_len; ++offset) {
      double update = 0.0;
      for (UF_long dep = k; dep < kb; ++dep) {
        update += row_dense_panel[dep] *
                  trailing_panel[dep * trailing_len + offset];
      }
      row_panel[offset] -= update;
    }
  }
}

static int kls_row_refactor_assemble_compact_dense_row(
  kls_egraph_refactor_worker *worker,
  UF_long group,
  UF_long row_begin,
  UF_long row_end,
  UF_long trailing_len,
  const UF_long *trailing_cols,
  double *dense_panel,
  double *trailing_panel,
  int wait_for_dependencies,
  UF_long row);


/* Claim the next work index of the current help session. The cursor
   packs (epoch << 32 | index); a claim succeeds only when the caller's
   session epoch matches, so a participant holding a stale view of the
   session (acquire loads may lag; only RMWs are guaranteed fresh) backs
   off WITHOUT consuming an index. A plain fetch-add here loses one
   index to exactly that straggler and deadlocks the owner drain on
   done < target (observed live: cursor exhausted, inflight 0, done one
   short). Returns the claimed index or -1. */
static long kls_dense_help_claim(kls_egraph_refactor_shared *shared,
                                 long epoch,
                                 long limit) {
  long v = atomic_load_explicit(&shared->dense_help_cursor,
                                memory_order_acquire);
  for (;;) {
    if ((v >> 32) != epoch) {
      return -1;
    }
    const long idx = v & 0xffffffffl;
    if (idx >= limit) {
      return -1;
    }
    if (atomic_compare_exchange_weak_explicit(
          &shared->dense_help_cursor, &v, v + 1,
          memory_order_acq_rel, memory_order_acquire)) {
      return idx;
    }
  }
}

static _Thread_local int kls_dense_help_busy;
/* Lowest row the current thread might still produce: the first row of the
   group it is processing. Rows below it are finished or owned by other
   threads, so groups ending at or before the floor are safe to help.
   Defaults to 0 (help nothing) until a worker enters a group. */
static _Thread_local long kls_dense_help_row_floor;

/* Assemble rows of another worker's dense group while blocked on it.
   row_bound is the smallest row the caller might still produce: helping
   is restricted to groups that end at or before it, so every dependency
   of a helped row lies strictly below any row the caller has in flight
   (wait chains strictly decrease in row order and cannot cycle). The
   helper runs on a private scratch vector so the caller's in-progress
   row accumulation in worker->x stays intact. Returns nonzero when at
   least one row was assembled. */

static int kls_compact_dense_panel_factor_blocked_unchecked(
  kls_egraph_refactor_worker *worker,
  UF_long group,
  UF_long row_begin,
  UF_long row_end,
  UF_long trailing_len,
  double *dense_panel,
  double *trailing_panel,
  int publish_prefix) {
  if (worker == NULL || worker->shared == NULL ||
      worker->shared->check_pivots ||
      dense_panel == NULL ||
      row_begin >= row_end ||
      (trailing_len > 0u && trailing_panel == NULL)) {
    return 0;
  }
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  if (solver == NULL || solver->numeric == NULL ||
      solver->numeric->Udiag == NULL ||
      solver->symbolic == NULL || solver->symbolic->Q == NULL ||
      row_end > solver->n) {
    kls_egraph_refactor_record_invalid(shared);
    return -1;
  }
  const UF_long width = row_end - row_begin;
  if (width <= 1u) {
    return 0;
  }

  double *udiag = (double *)solver->numeric->Udiag;
  for (UF_long i = 0; i < width; ++i) {
    dense_panel[i * width + i] = udiag[row_begin + i];
  }

  int panel_coop = 0;
  /* One ctx per session kind: helpers hold a published pointer, so a
     struct must never be written again once its session went active -
     the compiler is free to hoist plain stores for the NEXT session
     above the current drain loop otherwise (observed as a lost claim
     and a hung drain). */
  struct kls_dense_help_ctx panel_step_ctx;
  struct kls_dense_help_ctx panel_tail_ctx;
  if (shared->thread_count > 1 && width >= 16) {
    int expected = 0;
    if (atomic_compare_exchange_strong_explicit(
          &shared->dense_help_owner_lock, &expected, 1,
          memory_order_acq_rel, memory_order_relaxed)) {
      panel_coop = 1;
      panel_step_ctx.kind = KLS_DENSE_HELP_PANEL_STEP;
      panel_step_ctx.group = group;
      panel_step_ctx.row_begin = row_begin;
      panel_step_ctx.row_end = row_end;
      panel_step_ctx.trailing_len = trailing_len;
      panel_step_ctx.trailing_cols = NULL;
      panel_step_ctx.dense_panel = dense_panel;
      panel_step_ctx.trailing_panel = trailing_panel;
      panel_step_ctx.wait_for_dependencies = 0;
      panel_step_ctx.panel_k = 0;
      panel_step_ctx.panel_kb = 0;
      panel_tail_ctx = panel_step_ctx;
      panel_tail_ctx.kind = KLS_DENSE_HELP_PANEL_BLOCK_TAIL;
    }
  }

  double update_work = 0.0;
  for (UF_long k = 0; k < width; k += KLS_ROW_REFACTOR_DENSE_BLOCK_ROWS) {
    UF_long block = KLS_ROW_REFACTOR_DENSE_BLOCK_ROWS;
    if (block > width - k) {
      block = width - k;
    }
    const UF_long kb = k + block;

    const double tail_flops =
      0.5 * (double)block * (double)block *
      (double)(width - kb + trailing_len);
    const int coop_block_tail = panel_coop && kb < width &&
                                tail_flops >= 16384.0;
    for (UF_long i = k; i < kb; ++i) {
      double *row_dense_panel = dense_panel + i * width;
      double *row_panel = trailing_len > 0u
        ? trailing_panel + i * trailing_len : NULL;
      const UF_long triangle_end = coop_block_tail ? kb : width;
      for (UF_long dep = k; dep < i; ++dep) {
        const double pivot = udiag[row_begin + dep];
        const double lij = row_dense_panel[dep] / pivot;
        row_dense_panel[dep] = lij;
        const double *dep_dense_panel = dense_panel + dep * width;
        for (UF_long target = dep + 1u; target < triangle_end; ++target) {
          row_dense_panel[target] -= lij * dep_dense_panel[target];
        }
        update_work += (double)(width - dep - 1u);
        if (trailing_len > 0u && !coop_block_tail) {
          const double *dep_panel = trailing_panel + dep * trailing_len;
          for (UF_long offset = 0; offset < trailing_len; ++offset) {
            row_panel[offset] -= lij * dep_panel[offset];
          }
          update_work += (double)trailing_len;
        }
      }

      udiag[row_begin + i] = row_dense_panel[i];
      if (udiag[row_begin + i] == 0.0) {
        kls_egraph_refactor_record_singular(shared, row_begin + i,
                                            solver->symbolic->Q[row_begin + i]);
        if (solver->common.halt_if_singular) {
          if (panel_coop) {
            atomic_store_explicit(&shared->dense_help_owner_lock, 0,
                                  memory_order_release);
          }
          return -1;
        }
      }
    }
    if (coop_block_tail) {
      /* Block-row tails: a forward substitution across the block rows,
         parallel by disjoint column slices. */
      const UF_long tail_total = width - kb + trailing_len;
      const long slice_count =
        (long)((tail_total + KLS_DENSE_HELP_TAIL_SLICE - 1u) /
               KLS_DENSE_HELP_TAIL_SLICE);
      panel_tail_ctx.panel_k = k;
      panel_tail_ctx.panel_kb = kb;
      const long tail_epoch = ++shared->dense_help_epoch_seq;
      panel_tail_ctx.epoch = tail_epoch;
      atomic_store_explicit(&shared->dense_help_ctx, &panel_tail_ctx,
                            memory_order_relaxed);
      atomic_store_explicit(&shared->dense_help_cursor, tail_epoch << 32,
                            memory_order_release);
      atomic_store_explicit(&shared->dense_help_done, 0,
                            memory_order_relaxed);
      atomic_store_explicit(&shared->dense_help_active, 1,
                            memory_order_release);
      for (;;) {
        const long claim =
          kls_dense_help_claim(shared, tail_epoch, slice_count);
        if (claim < 0) {
          break;
        }
        const UF_long off_begin =
          (UF_long)claim * KLS_DENSE_HELP_TAIL_SLICE;
        const UF_long off_end =
          off_begin + KLS_DENSE_HELP_TAIL_SLICE < tail_total
            ? off_begin + KLS_DENSE_HELP_TAIL_SLICE : tail_total;
        kls_compact_dense_panel_block_tail_slice(
          dense_panel, trailing_panel, width, trailing_len, k, kb,
          off_begin, off_end);
        atomic_fetch_add_explicit(&shared->dense_help_done, 1,
                                  memory_order_release);
      }
      atomic_store_explicit(&shared->dense_help_active, 0,
                            memory_order_release);
      while (atomic_load_explicit(&shared->dense_help_inflight,
                                  memory_order_acquire) != 0 ||
             atomic_load_explicit(&shared->dense_help_done,
                                  memory_order_acquire) < slice_count) {
        kls_cpu_relax();
      }
      atomic_store_explicit(&shared->dense_help_ctx, NULL,
                            memory_order_release);
      update_work += tail_flops;
    }

    if (kb >= width) {
      continue;
    }

    const double step_flops =
      (double)(width - kb) *
      ((double)block * (double)(width - kb + trailing_len) +
       0.5 * (double)block * (double)block);
    if (panel_coop && step_flops >= 16384.0) {
      panel_step_ctx.panel_k = k;
      panel_step_ctx.panel_kb = kb;
      const long step_epoch = ++shared->dense_help_epoch_seq;
      panel_step_ctx.epoch = step_epoch;
      atomic_store_explicit(&shared->dense_help_ctx, &panel_step_ctx,
                            memory_order_relaxed);
      atomic_store_explicit(&shared->dense_help_cursor,
                            (step_epoch << 32) | (long)kb,
                            memory_order_release);
      atomic_store_explicit(&shared->dense_help_done, 0,
                            memory_order_relaxed);
      atomic_store_explicit(&shared->dense_help_fail, 0,
                            memory_order_relaxed);
      atomic_store_explicit(&shared->dense_help_active, 1,
                            memory_order_release);
      for (;;) {
        const long claim =
          kls_dense_help_claim(shared, step_epoch, (long)width);
        if (claim < 0) {
          break;
        }
        kls_compact_dense_panel_block_step_row(
          dense_panel, trailing_panel, udiag + row_begin, width,
          trailing_len, k, kb, (UF_long)claim);
        atomic_fetch_add_explicit(&shared->dense_help_done, 1,
                                  memory_order_release);
      }
      atomic_store_explicit(&shared->dense_help_active, 0,
                            memory_order_release);
      while (atomic_load_explicit(&shared->dense_help_inflight,
                                  memory_order_acquire) != 0 ||
             atomic_load_explicit(&shared->dense_help_done,
                                  memory_order_acquire) <
               (long)(width - kb)) {
        kls_cpu_relax();
      }
      atomic_store_explicit(&shared->dense_help_ctx, NULL,
                            memory_order_release);
      update_work += step_flops;
      continue;
    }

    for (UF_long row = kb; row < width; ++row) {
      double *row_dense_panel = dense_panel + row * width;
      for (UF_long dep = k; dep < kb; ++dep) {
        double value = row_dense_panel[dep];
        for (UF_long prev = k; prev < dep; ++prev) {
          value -= row_dense_panel[prev] * dense_panel[prev * width + dep];
          update_work += 1.0;
        }
        row_dense_panel[dep] = value / udiag[row_begin + dep];
      }
    }

    for (UF_long row = kb; row < width; ++row) {
      double *row_dense_panel = dense_panel + row * width;
      for (UF_long target = kb; target < width; ++target) {
        double update = 0.0;
        for (UF_long dep = k; dep < kb; ++dep) {
          update += row_dense_panel[dep] *
                    dense_panel[dep * width + target];
        }
        row_dense_panel[target] -= update;
      }
      update_work += (double)block * (double)(width - kb);
      if (trailing_len > 0u) {
        double *row_panel = trailing_panel + row * trailing_len;
        for (UF_long offset = 0; offset < trailing_len; ++offset) {
          double update = 0.0;
          for (UF_long dep = k; dep < kb; ++dep) {
            update += row_dense_panel[dep] *
                      trailing_panel[dep * trailing_len + offset];
          }
          row_panel[offset] -= update;
        }
        update_work += (double)block * (double)trailing_len;
      }
    }
  }

  if (panel_coop) {
    atomic_store_explicit(&shared->dense_help_owner_lock, 0,
                          memory_order_release);
  }

  for (UF_long row = row_begin; row < row_end; ++row) {
    if (udiag[row] == 0.0 && solver->common.halt_if_singular) {
      return -1;
    }
    if (!kls_store_compact_dense_panel_row(shared, row_begin, row_end,
                                          trailing_len, dense_panel,
                                          trailing_panel, row, 1)) {
      kls_egraph_refactor_record_invalid(shared);
      return -1;
    }
    if (publish_prefix) {
      kls_mark_row_refactor_compact_panel_prefix_valid(
        solver, group, row - row_begin + 1u);
    }
    kls_row_refactor_mark_row_done(worker, row);
  }
  const UF_long recorded_work = kls_work_to_uflong(update_work);
  kls_record_row_refactor_compact_dense_panel_blocked(
    solver, width, recorded_work);
  if (shared->row_refactor_native_row_panel_active) {
    kls_record_row_refactor_native_row_panel_blocked(
      solver, width, recorded_work);
  }
  return 1;
}

static int kls_compact_dense_panel_factor_blocked_checked(
  kls_egraph_refactor_worker *worker,
  UF_long group,
  UF_long row_begin,
  UF_long row_end,
  UF_long trailing_len,
  double *dense_panel,
  double *trailing_panel,
  int publish_prefix) {
  if (worker == NULL || worker->shared == NULL ||
      !worker->shared->check_pivots ||
      dense_panel == NULL ||
      row_begin >= row_end ||
      (trailing_len > 0u && trailing_panel == NULL)) {
    return 0;
  }
  const UF_long width = row_end - row_begin;
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  if (solver == NULL || solver->numeric == NULL ||
      solver->numeric->Udiag == NULL ||
      solver->symbolic == NULL || solver->symbolic->Q == NULL ||
      row_end > solver->n) {
    kls_egraph_refactor_record_invalid(shared);
    return -1;
  }
  if (width <= 1u) {
    return 0;
  }

  double *x = worker->x;
  double *udiag = (double *)solver->numeric->Udiag;
  for (UF_long i = 0; i < width; ++i) {
    dense_panel[i * width + i] = udiag[row_begin + i];
  }

  double update_work = 0.0;
  for (UF_long k = 0; k < width; k += KLS_ROW_REFACTOR_DENSE_BLOCK_ROWS) {
    UF_long block = KLS_ROW_REFACTOR_DENSE_BLOCK_ROWS;
    if (block > width - k) {
      block = width - k;
    }
    const UF_long kb = k + block;

    for (UF_long i = k; i < kb; ++i) {
      const UF_long row = row_begin + i;
      double *row_dense_panel = dense_panel + i * width;
      double *row_panel = trailing_len > 0u
        ? trailing_panel + i * trailing_len : NULL;
      for (UF_long dep = k; dep < i; ++dep) {
        const UF_long dep_row = row_begin + dep;
        const double candidate = row_dense_panel[dep];
        const double lij = candidate / udiag[dep_row];
        if (kls_parallel_row_refactor_rejects_multiplier(
              worker, row, dep_row, candidate, lij)) {
          return -1;
        }
        row_dense_panel[dep] = lij;
        const double *dep_dense_panel = dense_panel + dep * width;
        for (UF_long target = dep + 1u; target < width; ++target) {
          row_dense_panel[target] -= lij * dep_dense_panel[target];
        }
        update_work += (double)(width - dep - 1u);
        if (trailing_len > 0u) {
          const double *dep_panel = trailing_panel + dep * trailing_len;
          for (UF_long offset = 0; offset < trailing_len; ++offset) {
            row_panel[offset] -= lij * dep_panel[offset];
          }
          update_work += (double)trailing_len;
        }
      }

      udiag[row] = row_dense_panel[i];
      if (!kls_store_compact_dense_panel_row(shared, row_begin, row_end,
                                            trailing_len, dense_panel,
                                            trailing_panel, row, 0)) {
        kls_egraph_refactor_record_invalid(shared);
        return -1;
      }

      const UF_long row_dense_len = width - i - 1u;
      double row_max_abs =
        kls_row_refactor_abs_max_from_values(
          udiag[row], row_dense_panel + i + 1u, row_dense_len);
      if (trailing_len > 0u) {
        row_max_abs =
          kls_row_refactor_abs_max_from_values(row_max_abs, row_panel,
                                               trailing_len);
      }
      if (kls_parallel_row_refactor_rejects_pivot(worker, row, udiag[row],
                                                 row_max_abs)) {
        memset(x, 0, (size_t)solver->n * sizeof(*x));
        return -1;
      }
      if (udiag[row] == 0.0) {
        kls_egraph_refactor_record_singular(shared, row,
                                            solver->symbolic->Q[row]);
        if (solver->common.halt_if_singular) {
          memset(x, 0, (size_t)solver->n * sizeof(*x));
          return -1;
        }
      }
      if (!shared->row_refactor_defer_value_scatter &&
          !kls_scatter_row_refactor_row_values(solver, row)) {
        kls_egraph_refactor_record_invalid(shared);
        return -1;
      }
      if (publish_prefix) {
        kls_mark_row_refactor_compact_panel_prefix_valid(
          solver, group, i + 1u);
      }
      kls_row_refactor_mark_row_done(worker, row);
    }

    if (kb >= width) {
      continue;
    }

    for (UF_long row = kb; row < width; ++row) {
      const UF_long global_row = row_begin + row;
      double *row_dense_panel = dense_panel + row * width;
      for (UF_long dep = k; dep < kb; ++dep) {
        double candidate = row_dense_panel[dep];
        for (UF_long prev = k; prev < dep; ++prev) {
          candidate -= row_dense_panel[prev] *
                       dense_panel[prev * width + dep];
          update_work += 1.0;
        }
        const UF_long dep_row = row_begin + dep;
        const double lij = candidate / udiag[dep_row];
        if (kls_parallel_row_refactor_rejects_multiplier(
              worker, global_row, dep_row, candidate, lij)) {
          return -1;
        }
        row_dense_panel[dep] = lij;
      }
    }

    for (UF_long row = kb; row < width; ++row) {
      double *row_dense_panel = dense_panel + row * width;
      for (UF_long target = kb; target < width; ++target) {
        double update = 0.0;
        for (UF_long dep = k; dep < kb; ++dep) {
          update += row_dense_panel[dep] *
                    dense_panel[dep * width + target];
        }
        row_dense_panel[target] -= update;
      }
      update_work += (double)block * (double)(width - kb);
      if (trailing_len > 0u) {
        double *row_panel = trailing_panel + row * trailing_len;
        for (UF_long offset = 0; offset < trailing_len; ++offset) {
          double update = 0.0;
          for (UF_long dep = k; dep < kb; ++dep) {
            update += row_dense_panel[dep] *
                      trailing_panel[dep * trailing_len + offset];
          }
          row_panel[offset] -= update;
        }
        update_work += (double)block * (double)trailing_len;
      }
    }
  }

  const UF_long recorded_work = kls_work_to_uflong(update_work);
  kls_record_row_refactor_compact_dense_panel_blocked(
    solver, width, recorded_work);
  if (shared->row_refactor_native_row_panel_active) {
    kls_record_row_refactor_native_row_panel_blocked(
      solver, width, recorded_work);
  }
  return 1;
}

static UF_long kls_work_to_uflong(double work) {
  if (!(work > 0.0)) {
    return 0;
  }
  if (work >= (double)UF_long_max) {
    return UF_long_max;
  }
  return (UF_long)work;
}

typedef struct {
  UF_long group_begin;
  UF_long width;
  UF_long suffix_begin;
  UF_long dep_rows;
  UF_long external_offset;
  UF_long trailing_len;
  const UF_long *trailing_cols;
  double *dense_panel;
  double *trailing_panel;
} kls_batched_supernode_run;

static UF_long kls_compact_dense_group_external_pattern_end(
  const kls_solver *solver,
  UF_long row_begin,
  UF_long row_end,
  UF_long batch_begin) {
  if (solver == NULL ||
      solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_l_cols == NULL ||
      solver->row_refactor_l_internal_ptr == NULL ||
      batch_begin < row_begin ||
      batch_begin >= row_end) {
    return batch_begin;
  }
  const UF_long first_begin = solver->row_refactor_l_ptr[batch_begin];
  const UF_long first_internal =
    solver->row_refactor_l_internal_ptr[batch_begin];
  if (first_begin >= first_internal) {
    return batch_begin + 1u;
  }
  const UF_long external_len = first_internal - first_begin;
  UF_long batch_end = batch_begin + 1u;
  while (batch_end < row_end) {
    const UF_long begin = solver->row_refactor_l_ptr[batch_end];
    const UF_long internal =
      solver->row_refactor_l_internal_ptr[batch_end];
    if (begin > internal || internal - begin != external_len) {
      break;
    }
    int same = 1;
    for (UF_long offset = 0; offset < external_len; ++offset) {
      if (solver->row_refactor_l_cols[begin + offset] !=
          solver->row_refactor_l_cols[first_begin + offset]) {
        same = 0;
        break;
      }
    }
    if (!same) {
      break;
    }
    batch_end++;
  }
  return batch_end;
}

static int kls_compact_dense_group_single_producer_suffix(
  const kls_solver *solver,
  UF_long current_group,
  UF_long row_begin,
  UF_long row,
  UF_long *dep_group_out,
  UF_long *suffix_begin_out,
  UF_long *dep_rows_out) {
  if (dep_group_out != NULL) {
    *dep_group_out = KLS_KLU_EMPTY;
  }
  KLS_SET_OPTIONAL_OUTPUT(suffix_begin_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(dep_rows_out, 0);
  if (solver == NULL ||
      row >= solver->n ||
      solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_l_cols == NULL ||
      solver->row_refactor_l_internal_ptr == NULL ||
      solver->row_refactor_row_group == NULL ||
      solver->row_refactor_group_ptr == NULL ||
      solver->row_refactor_group_kind == NULL) {
    return 0;
  }

  const UF_long l_begin = solver->row_refactor_l_ptr[row];
  const UF_long l_internal = solver->row_refactor_l_internal_ptr[row];
  if (l_begin >= l_internal) {
    return 0;
  }
  const UF_long dep0 = solver->row_refactor_l_cols[l_begin];
  if (dep0 >= row_begin || dep0 >= row || dep0 >= solver->n) {
    return 0;
  }
  const UF_long dep_group = solver->row_refactor_row_group[dep0];
  if (dep_group >= solver->row_refactor_group_count ||
      dep_group == current_group ||
      solver->row_refactor_group_kind[dep_group] !=
        KLS_ROW_REFACTOR_GROUP_DENSE) {
    return 0;
  }
  const UF_long group_begin = solver->row_refactor_group_ptr[dep_group];
  const UF_long group_end = solver->row_refactor_group_ptr[dep_group + 1u];
  if (group_begin >= group_end ||
      group_end > row_begin ||
      dep0 < group_begin ||
      dep0 >= group_end) {
    return 0;
  }
  const UF_long width = group_end - group_begin;
  if (kls_row_refactor_compact_panel_valid_prefix(solver, dep_group) <
      width) {
    return 0;
  }

  UF_long p = l_begin;
  UF_long expected = dep0;
  while (p < l_internal && expected < group_end) {
    const UF_long dep = solver->row_refactor_l_cols[p];
    if (dep != expected ||
        dep >= row ||
        solver->row_refactor_row_group[dep] != dep_group) {
      return 0;
    }
    p++;
    expected++;
  }
  if (p != l_internal || expected != group_end) {
    return 0;
  }

  if (dep_group_out != NULL) {
    *dep_group_out = dep_group;
  }
  if (suffix_begin_out != NULL) {
    *suffix_begin_out = dep0 - group_begin;
  }
  if (dep_rows_out != NULL) {
    *dep_rows_out = l_internal - l_begin;
  }
  return 1;
}

static UF_long kls_compact_dense_group_ragged_supernode_end(
  const kls_solver *solver,
  UF_long group,
  UF_long row_begin,
  UF_long row_end,
  UF_long batch_begin,
  UF_long *dep_group_out,
  double *work_out,
  double *copied_entries_out,
  UF_long *dep_rows_out) {
  if (dep_group_out != NULL) {
    *dep_group_out = KLS_KLU_EMPTY;
  }
  if (work_out != NULL) {
    *work_out = 0.0;
  }
  if (copied_entries_out != NULL) {
    *copied_entries_out = 0.0;
  }
  KLS_SET_OPTIONAL_OUTPUT(dep_rows_out, 0);
  if (solver == NULL ||
      solver->row_refactor_group_trailing_len == NULL ||
      batch_begin < row_begin ||
      batch_begin >= row_end) {
    return batch_begin;
  }

  UF_long first_dep_group = KLS_KLU_EMPTY;
  UF_long dep_rows = 0;
  if (!kls_compact_dense_group_single_producer_suffix(
        solver, group, row_begin, batch_begin, &first_dep_group,
        NULL, &dep_rows)) {
    return batch_begin + 1u;
  }

  const UF_long dep_trailing_len =
    solver->row_refactor_group_trailing_len[first_dep_group];
  double work = 0.5 * (double)dep_rows * (double)(dep_rows - 1u) +
    (double)dep_rows * (double)dep_trailing_len;
  double copied = (double)dep_rows + (double)dep_trailing_len;
  double total_dep_rows = (double)dep_rows;
  UF_long batch_end = batch_begin + 1u;
  while (batch_end < row_end) {
    UF_long dep_group = KLS_KLU_EMPTY;
    UF_long row_dep_rows = 0;
    if (!kls_compact_dense_group_single_producer_suffix(
          solver, group, row_begin, batch_end, &dep_group,
          NULL, &row_dep_rows) ||
        dep_group != first_dep_group) {
      break;
    }
    work += 0.5 * (double)row_dep_rows * (double)(row_dep_rows - 1u) +
      (double)row_dep_rows * (double)dep_trailing_len;
    copied += (double)row_dep_rows + (double)dep_trailing_len;
    total_dep_rows += (double)row_dep_rows;
    batch_end++;
  }

  if (batch_end <= batch_begin + 1u) {
    return batch_begin + 1u;
  }
  if (dep_group_out != NULL) {
    *dep_group_out = first_dep_group;
  }
  if (work_out != NULL) {
    *work_out = work;
  }
  if (copied_entries_out != NULL) {
    *copied_entries_out = copied;
  }
  if (dep_rows_out != NULL) {
    *dep_rows_out = kls_work_to_uflong(total_dep_rows);
  }
  return batch_end;
}

static int kls_dense_help_try(kls_egraph_refactor_shared *shared,
                              long row_bound,
                              int idle) {
  if (!atomic_load_explicit(&shared->dense_help_active,
                            memory_order_acquire)) {
    return 0;
  }
  kls_egraph_refactor_worker *self = kls_dense_help_self;
  if (self == NULL || self->shared != shared || kls_dense_help_busy) {
    return 0;
  }
  {
    const struct kls_dense_help_ctx *peek =
      atomic_load_explicit(&shared->dense_help_ctx, memory_order_acquire);
    if (peek == NULL ||
        (peek->kind == KLS_DENSE_HELP_ASSEMBLE &&
         (long)peek->row_end > row_bound)) {
      return 0;
    }
    if (peek != NULL && peek->kind == KLS_DENSE_HELP_ASSEMBLE &&
        self->help_x == NULL) {
      kls_solver *solver = shared->solver;
      if (solver == NULL) {
        return 0;
      }
      self->help_x = (double *)calloc((size_t)solver->n, sizeof(double));
      if (self->help_x == NULL) {
        return 0;
      }
    }
  }
  int did = 0;
  atomic_fetch_add_explicit(&shared->dense_help_inflight, 1,
                            memory_order_acq_rel);
  if (atomic_load_explicit(&shared->dense_help_active,
                           memory_order_acquire) &&
      !atomic_load_explicit(&shared->dense_help_fail,
                            memory_order_acquire)) {
    const struct kls_dense_help_ctx *c =
      atomic_load_explicit(&shared->dense_help_ctx, memory_order_acquire);
    if (c != NULL && (c->kind == KLS_DENSE_HELP_PANEL_STEP ||
                      c->kind == KLS_DENSE_HELP_PANEL_BLOCK_TAIL)) {
      /* Panel rows read only finalized panel rows and write only the
         claimed row: no waits, so any blocked or idle thread may join
         regardless of its own position in the pipeline. */
      const int block_tail = c->kind == KLS_DENSE_HELP_PANEL_BLOCK_TAIL;
      const UF_long panel_width = c->row_end - c->row_begin;
      const UF_long tail_total = block_tail
        ? panel_width - c->panel_kb + c->trailing_len : 0u;
      const long claim_limit = block_tail
        ? (long)((tail_total + KLS_DENSE_HELP_TAIL_SLICE - 1u) /
                 KLS_DENSE_HELP_TAIL_SLICE)
        : (long)panel_width;
      const double *udiag_base =
        (const double *)shared->solver->numeric->Udiag + c->row_begin;
      kls_dense_help_busy = 1;
      const long panel_epoch = c->epoch;
      for (;;) {
        const long row =
          kls_dense_help_claim(shared, panel_epoch, claim_limit);
        if (row < 0) {
          break;
        }
        if (block_tail) {
          const UF_long off_begin =
            (UF_long)row * KLS_DENSE_HELP_TAIL_SLICE;
          const UF_long off_end =
            off_begin + KLS_DENSE_HELP_TAIL_SLICE < tail_total
              ? off_begin + KLS_DENSE_HELP_TAIL_SLICE : tail_total;
          kls_compact_dense_panel_block_tail_slice(
            c->dense_panel, c->trailing_panel, panel_width,
            c->trailing_len, c->panel_k, c->panel_kb, off_begin, off_end);
        } else {
          kls_compact_dense_panel_block_step_row(
            c->dense_panel, c->trailing_panel, udiag_base, panel_width,
            c->trailing_len, c->panel_k, c->panel_kb, (UF_long)row);
        }
        atomic_fetch_add_explicit(&shared->dense_help_done, 1,
                                  memory_order_release);
        did = 1;
        if (!atomic_load_explicit(&shared->dense_help_active,
                                  memory_order_acquire)) {
          break;
        }
      }
      kls_dense_help_busy = 0;
      atomic_fetch_sub_explicit(&shared->dense_help_inflight, 1,
                                memory_order_release);
      return did;
    }
    int allowed = c != NULL && c->kind == KLS_DENSE_HELP_ASSEMBLE &&
                  (long)c->row_end <= row_bound && self->help_x != NULL;
    int idle_counted = 0;
    if (allowed && idle && c->wait_for_dependencies) {
      /* Helped rows may block on unfinished dependencies. Keep at least
         one idle thread free to claim newly ready groups so a blocking
         group cannot sit unclaimed while every thread helps. */
      const int idlers = atomic_fetch_add_explicit(
        &shared->dense_help_idle_inflight, 1, memory_order_acq_rel) + 1;
      if (idlers > shared->thread_count - 2) {
        atomic_fetch_sub_explicit(&shared->dense_help_idle_inflight, 1,
                                  memory_order_release);
        allowed = 0;
      } else {
        idle_counted = 1;
      }
    }
    if (allowed) {
      double *saved_x = self->x;
      self->x = self->help_x;
      kls_dense_help_busy = 1;
      const long assemble_epoch = c->epoch;
      for (;;) {
        const long row =
          kls_dense_help_claim(shared, assemble_epoch, (long)c->row_end);
        if (row < 0) {
          break;
        }
        const int ok = kls_row_refactor_assemble_compact_dense_row(
          self, c->group, c->row_begin, c->row_end, c->trailing_len,
          c->trailing_cols, c->dense_panel, c->trailing_panel,
          c->wait_for_dependencies, (UF_long)row);
        atomic_fetch_add_explicit(&shared->dense_help_done, 1,
                                  memory_order_release);
        did = 1;
        if (!ok) {
          atomic_store_explicit(&shared->dense_help_fail, 1,
                                memory_order_release);
          memset(self->help_x, 0,
                 (size_t)shared->solver->n * sizeof(*self->help_x));
          break;
        }
        if (atomic_load_explicit(&shared->dense_help_fail,
                                 memory_order_acquire) ||
            !atomic_load_explicit(&shared->dense_help_active,
                                  memory_order_acquire)) {
          break;
        }
      }
      kls_dense_help_busy = 0;
      self->x = saved_x;
    }
    if (idle_counted) {
      atomic_fetch_sub_explicit(&shared->dense_help_idle_inflight, 1,
                                memory_order_release);
    }
  }
  atomic_fetch_sub_explicit(&shared->dense_help_inflight, 1,
                            memory_order_release);
  return did;
}

static int kls_compact_dense_group_try_batched_supernode_update(
  kls_egraph_refactor_worker *worker,
  UF_long group,
  UF_long row_begin,
  UF_long row_end,
  UF_long batch_begin,
  UF_long batch_end,
  UF_long trailing_len,
  const UF_long *trailing_cols,
  double *dense_panel,
  double *trailing_panel) {
  if (worker == NULL || worker->shared == NULL ||
      dense_panel == NULL ||
      row_begin >= row_end ||
      batch_begin < row_begin ||
      batch_begin >= batch_end ||
      batch_end > row_end) {
    return 0;
  }
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  if (solver == NULL ||
      solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_l_cols == NULL ||
      solver->row_refactor_l_values == NULL ||
      solver->row_refactor_l_row_values == NULL ||
      solver->row_refactor_l_internal_ptr == NULL ||
      solver->row_refactor_u_ptr == NULL ||
      solver->row_refactor_u_cols == NULL ||
      solver->row_refactor_row_group == NULL ||
      solver->row_refactor_group_ptr == NULL ||
      solver->row_refactor_group_kind == NULL ||
      solver->row_refactor_group_trailing_len == NULL ||
      solver->numeric == NULL ||
      solver->numeric->Udiag == NULL ||
      worker->x == NULL ||
      (trailing_len > 0u &&
       (trailing_cols == NULL || trailing_panel == NULL))) {
    return 0;
  }

  const UF_long batch_rows = batch_end - batch_begin;
  if (batch_rows == 0u || batch_rows > (UF_long)INT_MAX) {
    return 0;
  }
  const UF_long first_l_begin = solver->row_refactor_l_ptr[batch_begin];
  const UF_long first_l_internal =
    solver->row_refactor_l_internal_ptr[batch_begin];
  if (first_l_begin >= first_l_internal) {
    return 0;
  }
  const UF_long external_len = first_l_internal - first_l_begin;
  if (external_len == 0u || external_len > (UF_long)INT_MAX ||
      external_len > (UF_long)(SIZE_MAX / sizeof(UF_long)) ||
      external_len > (UF_long)(SIZE_MAX /
                               sizeof(kls_batched_supernode_run))) {
    return 0;
  }

  kls_batched_supernode_run *runs =
    (kls_batched_supernode_run *)kls_egraph_worker_object_workspace(
      worker, external_len, sizeof(*runs));
  UF_long *external_cols =
    kls_egraph_worker_index_workspace(worker, external_len);
  if (runs == NULL || external_cols == NULL) {
    return 0;
  }
  memset(runs, 0, (size_t)external_len * sizeof(*runs));

  int status = 0;
  UF_long run_count = 0;
  UF_long external_offset = 0;
  UF_long max_dep_trailing_len = 0;
  double update_work = 0.0;
  double copied_entries = 0.0;
  UF_long p = first_l_begin;
  while (p < first_l_internal) {
    const UF_long dep0 = solver->row_refactor_l_cols[p];
    if (dep0 >= row_begin || dep0 >= solver->n) {
      goto cleanup;
    }
    const UF_long dep_group = solver->row_refactor_row_group[dep0];
    if (dep_group >= solver->row_refactor_group_count ||
        dep_group == group ||
        solver->row_refactor_group_kind[dep_group] !=
          KLS_ROW_REFACTOR_GROUP_DENSE ||
        kls_row_refactor_compact_panel_valid_prefix(solver, dep_group) == 0u) {
      goto cleanup;
    }
    const UF_long dep_group_begin =
      solver->row_refactor_group_ptr[dep_group];
    const UF_long dep_group_end =
      solver->row_refactor_group_ptr[dep_group + 1u];
    if (dep_group_begin >= dep_group_end ||
        dep0 < dep_group_begin ||
        dep0 >= dep_group_end ||
        dep_group_end > row_begin) {
      goto cleanup;
    }

    const UF_long dep_width = dep_group_end - dep_group_begin;
    const UF_long dep_suffix_begin = dep0 - dep_group_begin;
    UF_long run_end = p;
    UF_long expected_dep = dep0;
    while (run_end < first_l_internal && expected_dep < dep_group_end) {
      const UF_long dep = solver->row_refactor_l_cols[run_end];
      if (dep != expected_dep ||
          solver->row_refactor_row_group[dep] != dep_group ||
          (external_offset > 0u &&
           dep <= external_cols[external_offset - 1u])) {
        break;
      }
      external_cols[external_offset++] = dep;
      run_end++;
      expected_dep++;
    }
    if (expected_dep != dep_group_end || run_end <= p) {
      goto cleanup;
    }
    const UF_long dep_rows = run_end - p;
    if (dep_rows < 2u ||
        dep_width > (UF_long)INT_MAX ||
        dep_rows > (UF_long)INT_MAX) {
      goto cleanup;
    }

    const UF_long dep_trailing_len =
      solver->row_refactor_group_trailing_len[dep_group];
    if (dep_trailing_len > (UF_long)INT_MAX) {
      goto cleanup;
    }
    UF_long dep_dense_entries = 0;
    if (dep_width > (UF_long)(SIZE_MAX / sizeof(double)) / dep_width) {
      kls_egraph_refactor_record_invalid(shared);
      status = -1;
      goto cleanup;
    }
    dep_dense_entries = dep_width * dep_width;
    if (dep_trailing_len > 0u &&
        dep_width > ((UF_long)(SIZE_MAX / sizeof(double)) -
                       dep_dense_entries) / dep_trailing_len) {
      kls_egraph_refactor_record_invalid(shared);
      status = -1;
      goto cleanup;
    }
    const UF_long dep_panel_entries =
      dep_dense_entries + dep_width * dep_trailing_len;
    double *dep_dense_panel =
      kls_row_refactor_compact_group_panel(solver, dep_group,
                                           dep_panel_entries);
    if (dep_dense_panel == NULL) {
      goto cleanup;
    }
    double *dep_trailing_panel = dep_trailing_len > 0u
      ? dep_dense_panel + dep_dense_entries : NULL;
    const UF_long *dep_trailing_cols = NULL;
    if (!kls_row_refactor_group_trailing_slice(
          solver, dep_group, dep_group_end, NULL,
          &dep_trailing_cols)) {
      kls_egraph_refactor_record_invalid(shared);
      status = -1;
      goto cleanup;
    }
    if (dep_trailing_len > 0u &&
        (dep_trailing_panel == NULL || dep_trailing_cols == NULL)) {
      kls_egraph_refactor_record_invalid(shared);
      status = -1;
      goto cleanup;
    }

    runs[run_count].group_begin = dep_group_begin;
    runs[run_count].width = dep_width;
    runs[run_count].suffix_begin = dep_suffix_begin;
    runs[run_count].dep_rows = dep_rows;
    runs[run_count].external_offset = p - first_l_begin;
    runs[run_count].trailing_len = dep_trailing_len;
    runs[run_count].trailing_cols = dep_trailing_cols;
    runs[run_count].dense_panel = dep_dense_panel;
    runs[run_count].trailing_panel = dep_trailing_panel;
    run_count++;

    update_work += (double)batch_rows *
      (0.5 * (double)dep_rows * (double)(dep_rows - 1u) +
       (double)dep_rows * (double)dep_trailing_len);
    copied_entries +=
      (double)batch_rows * (double)(dep_rows + dep_trailing_len);
    if (dep_trailing_len > max_dep_trailing_len) {
      max_dep_trailing_len = dep_trailing_len;
    }
    p = run_end;
  }
  if (run_count == 0u || external_offset != external_len) {
    goto cleanup;
  }

  for (UF_long row = batch_begin; row < batch_end; ++row) {
    const UF_long l_begin = solver->row_refactor_l_ptr[row];
    const UF_long l_internal = solver->row_refactor_l_internal_ptr[row];
    if (l_begin > l_internal || l_internal - l_begin != external_len) {
      goto cleanup;
    }
    for (UF_long local_dep = 0; local_dep < external_len; ++local_dep) {
      const UF_long dep = solver->row_refactor_l_cols[l_begin + local_dep];
      if (dep != external_cols[local_dep] || dep >= row) {
        goto cleanup;
      }
    }
  }

  for (UF_long run_pos = 0; run_pos < run_count; ++run_pos) {
    const kls_batched_supernode_run *run = runs + run_pos;
    const UF_long run_external_end =
      run->external_offset + run->dep_rows;
    for (UF_long dep_col_pos = 0; dep_col_pos < run->trailing_len;
         ++dep_col_pos) {
      const UF_long col = run->trailing_cols[dep_col_pos];
      UF_long external_pos = KLS_KLU_EMPTY;
      if (col < row_begin &&
          kls_find_uflong_sorted(external_cols, external_len, col,
                                 &external_pos) &&
          external_pos < run_external_end) {
        goto cleanup;
      }
    }
  }

  if (update_work < KLS_ROW_REFACTOR_BATCH_SUPERNODE_MIN_WORK ||
      update_work <
        KLS_ROW_REFACTOR_BATCH_SUPERNODE_MIN_WORK_PER_ENTRY *
          copied_entries) {
    const UF_long candidate_dep_rows =
      kls_work_to_uflong((double)batch_rows * (double)external_len);
    kls_record_row_refactor_compact_supernode_batch_candidate(
      solver, batch_rows, candidate_dep_rows, 1);
    goto cleanup;
  }
  const UF_long candidate_dep_rows =
    kls_work_to_uflong((double)batch_rows * (double)external_len);
  kls_record_row_refactor_compact_supernode_batch_candidate(
    solver, batch_rows, candidate_dep_rows, 0);

  const UF_long max_workspace_entries =
    (UF_long)(SIZE_MAX / sizeof(*worker->supernode_workspace));
  if (batch_rows > max_workspace_entries / external_len) {
    goto cleanup;
  }
  const UF_long multiplier_entries = batch_rows * external_len;
  UF_long update_entries = 0;
  if (max_dep_trailing_len > 0u) {
    if (batch_rows > (max_workspace_entries - multiplier_entries) /
                         max_dep_trailing_len) {
      goto cleanup;
    }
    update_entries = batch_rows * max_dep_trailing_len;
  }
  double *workspace =
    kls_egraph_worker_supernode_workspace(
      worker, multiplier_entries + update_entries);
  if (workspace == NULL) {
    goto cleanup;
  }
  double *multipliers = workspace;
  double *updates = max_dep_trailing_len > 0u
    ? multipliers + multiplier_entries : NULL;
  double *x = worker->x;
  double *udiag = (double *)solver->numeric->Udiag;
  const UF_long width = row_end - row_begin;

  for (UF_long batch_local = 0; batch_local < batch_rows; ++batch_local) {
    const UF_long row = batch_begin + batch_local;
    const UF_long group_local = row - row_begin;
    const int direct_input =
      kls_parallel_row_refactor_load_compact_dense_input_row(
        shared, x, row_begin, row_end, trailing_len, trailing_cols,
        dense_panel, trailing_panel, row);
    if (direct_input < 0) {
      status = -1;
      goto cleanup;
    }
    if (!direct_input &&
        !kls_parallel_row_refactor_load_input_row(shared, x, row)) {
      status = -1;
      goto cleanup;
    }

    const UF_long l_begin = solver->row_refactor_l_ptr[row];
    for (UF_long local_dep = 0; local_dep < external_len; ++local_dep) {
      const UF_long dep = external_cols[local_dep];
      multipliers[batch_local * external_len + local_dep] = x[dep];
      x[dep] = 0.0;
      solver->row_refactor_l_row_values[l_begin + local_dep] = 0.0;
    }

    UF_long lp = solver->row_refactor_l_internal_ptr[row];
    const UF_long l_end = solver->row_refactor_l_ptr[row + 1u];
    double *row_dense_panel = dense_panel + group_local * width;
    for (UF_long dep = row_begin; dep < row; ++dep, ++lp) {
      if (lp >= l_end || solver->row_refactor_l_cols[lp] != dep) {
        kls_egraph_refactor_record_invalid(shared);
        status = -1;
        goto cleanup;
      }
      if (direct_input) {
        row_dense_panel[dep - row_begin] += x[dep];
      } else {
        row_dense_panel[dep - row_begin] = x[dep];
      }
      x[dep] = 0.0;
    }
    if (lp != l_end) {
      kls_egraph_refactor_record_invalid(shared);
      status = -1;
      goto cleanup;
    }

    udiag[row] =
      direct_input ? row_dense_panel[group_local] + x[row] : x[row];
    x[row] = 0.0;

    const UF_long dense_len = row_end - row - 1u;
    const UF_long u_begin = solver->row_refactor_u_ptr[row];
    const UF_long u_end = solver->row_refactor_u_ptr[row + 1u];
    if (u_end - u_begin != dense_len + trailing_len) {
      kls_egraph_refactor_record_invalid(shared);
      status = -1;
      goto cleanup;
    }
    for (UF_long offset = 0; offset < dense_len; ++offset) {
      const UF_long col = row + 1u + offset;
      if (direct_input) {
        row_dense_panel[group_local + 1u + offset] += x[col];
      } else {
        row_dense_panel[group_local + 1u + offset] = x[col];
      }
      x[col] = 0.0;
    }
    if (trailing_len > 0u) {
      double *row_panel = trailing_panel + group_local * trailing_len;
      for (UF_long offset = 0; offset < trailing_len; ++offset) {
        const UF_long col = trailing_cols[offset];
        if (direct_input) {
          row_panel[offset] += x[col];
        } else {
          row_panel[offset] = x[col];
        }
        x[col] = 0.0;
      }
    }
    kls_clear_row_refactor_input_residuals(solver, x, row);
  }

  double trsm_work = 0.0;
  double gemm_work = 0.0;
  for (UF_long run_pos = 0; run_pos < run_count; ++run_pos) {
    kls_batched_supernode_run *run = runs + run_pos;
    for (UF_long local_dep = run->suffix_begin; local_dep < run->width;
         ++local_dep) {
      run->dense_panel[local_dep * run->width + local_dep] =
        udiag[run->group_begin + local_dep];
    }

    {
      for (UF_long batch_local = 0; batch_local < batch_rows; ++batch_local) {
        const UF_long row = batch_begin + batch_local;
        double *row_multipliers =
          multipliers + batch_local * external_len + run->external_offset;
        for (UF_long local_dep = 0; local_dep < run->dep_rows; ++local_dep) {
          const UF_long panel_row = run->suffix_begin + local_dep;
          const UF_long dep = run->group_begin + panel_row;
          const double pivot = udiag[run->group_begin + panel_row];
          const double lij = row_multipliers[local_dep] / pivot;
          const double *dep_dense_panel =
            run->dense_panel + panel_row * run->width;
          if (kls_parallel_row_refactor_rejects_multiplier(
                worker, row, dep, row_multipliers[local_dep], lij)) {
            status = -1;
            goto cleanup;
          }
          row_multipliers[local_dep] = lij;
          for (UF_long target = local_dep + 1u;
               target < run->dep_rows; ++target) {
            row_multipliers[target] -=
              lij * dep_dense_panel[run->suffix_begin + target];
          }
        }
      }
    }

    for (UF_long batch_local = 0; batch_local < batch_rows; ++batch_local) {
      const UF_long row = batch_begin + batch_local;
      const UF_long l_begin = solver->row_refactor_l_ptr[row];
      for (UF_long local_dep = 0; local_dep < run->dep_rows; ++local_dep) {
        const UF_long p =
          l_begin + run->external_offset + local_dep;
        const double lij =
          multipliers[batch_local * external_len +
                      run->external_offset + local_dep];
        solver->row_refactor_l_row_values[p] = lij;
        if (!shared->row_refactor_defer_value_scatter) {
          *solver->row_refactor_l_values[p] = lij;
        }
      }
    }
    trsm_work += (double)batch_rows * 0.5 *
      (double)run->dep_rows * (double)(run->dep_rows - 1u);

    if (run->trailing_len == 0u) {
      continue;
    }
    memset(updates, 0,
           (size_t)batch_rows * (size_t)run->trailing_len *
             sizeof(*updates));
    {
      for (long batch_local = 0; batch_local < (long)batch_rows;
           ++batch_local) {
        double *row_updates =
          updates + (size_t)batch_local * run->trailing_len;
        const double *row_multipliers =
          multipliers + (size_t)batch_local * external_len +
          run->external_offset;
        for (UF_long local_dep = 0; local_dep < run->dep_rows; ++local_dep) {
          const double lij = row_multipliers[local_dep];
          const double *dep_panel =
            run->trailing_panel +
            (size_t)(run->suffix_begin + local_dep) * run->trailing_len;
          for (UF_long offset = 0; offset < run->trailing_len; ++offset) {
            row_updates[offset] += lij * dep_panel[offset];
          }
        }
      }
    }
    gemm_work +=
      (double)batch_rows * (double)run->dep_rows *
      (double)run->trailing_len;

    for (UF_long dep_col_pos = 0; dep_col_pos < run->trailing_len;
         ++dep_col_pos) {
      const UF_long col = run->trailing_cols[dep_col_pos];
      UF_long external_pos = KLS_KLU_EMPTY;
      UF_long trailing_pos = KLS_KLU_EMPTY;
      const int in_external =
        col < row_begin &&
        kls_find_uflong_sorted(external_cols, external_len, col,
                               &external_pos);
      const int in_current_group = col >= row_begin && col < row_end;
      const int in_current_trailing =
        !in_external && !in_current_group &&
        kls_find_uflong_sorted(trailing_cols, trailing_len, col,
                               &trailing_pos);
      if (!in_external && !in_current_group && !in_current_trailing) {
        continue;
      }
      for (UF_long batch_local = 0; batch_local < batch_rows; ++batch_local) {
        const double update =
          updates[batch_local * run->trailing_len + dep_col_pos];
        if (update == 0.0) {
          continue;
        }
        if (in_external) {
          multipliers[batch_local * external_len + external_pos] -= update;
        } else {
          const UF_long row = batch_begin + batch_local;
          const UF_long group_local = row - row_begin;
          if (in_current_group) {
            const UF_long local_col = col - row_begin;
            if (col == row) {
              udiag[row] -= update;
            } else {
              dense_panel[group_local * width + local_col] -= update;
            }
          } else {
            trailing_panel[group_local * trailing_len + trailing_pos] -=
              update;
          }
        }
      }
    }
  }

  const UF_long applied_dep_rows =
    kls_work_to_uflong((double)batch_rows * (double)external_len);
  const UF_long trsm_entries = kls_work_to_uflong(trsm_work);
  const UF_long gemm_entries = kls_work_to_uflong(gemm_work);
  const UF_long total_entries =
    trsm_entries > UF_long_max - gemm_entries
      ? UF_long_max : trsm_entries + gemm_entries;
  kls_record_row_refactor_compact_supernode_batch(
    solver, batch_rows, applied_dep_rows, total_entries);
  kls_record_row_refactor_compact_supernode_update(
    solver, applied_dep_rows, total_entries);
  kls_record_row_refactor_compact_supernode_trsv(
    solver, applied_dep_rows, trsm_entries);
  if (gemm_entries > 0u) {
    kls_record_row_refactor_compact_supernode_gemv(
      solver, applied_dep_rows, gemm_entries);
  }
  status = 1;

cleanup:
  return status;
}

static int kls_compact_dense_group_try_ragged_supernode_update(
  kls_egraph_refactor_worker *worker,
  UF_long group,
  UF_long row_begin,
  UF_long row_end,
  UF_long batch_begin,
  UF_long batch_end,
  UF_long dep_group,
  double update_work,
  double copied_entries,
  UF_long candidate_dep_rows,
  UF_long trailing_len,
  const UF_long *trailing_cols,
  double *dense_panel,
  double *trailing_panel) {
  if (worker == NULL || worker->shared == NULL ||
      dense_panel == NULL ||
      row_begin >= row_end ||
      batch_begin < row_begin ||
      batch_begin >= batch_end ||
      batch_end > row_end ||
      candidate_dep_rows == 0u) {
    return 0;
  }
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  if (solver == NULL ||
      dep_group >= solver->row_refactor_group_count ||
      solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_l_cols == NULL ||
      solver->row_refactor_l_row_values == NULL ||
      solver->row_refactor_l_values == NULL ||
      solver->row_refactor_l_internal_ptr == NULL ||
      solver->row_refactor_u_ptr == NULL ||
      solver->row_refactor_u_cols == NULL ||
      solver->row_refactor_row_group == NULL ||
      solver->row_refactor_group_ptr == NULL ||
      solver->row_refactor_group_kind == NULL ||
      solver->row_refactor_group_trailing_len == NULL ||
      solver->numeric == NULL ||
      solver->numeric->Udiag == NULL ||
      worker->x == NULL ||
      (trailing_len > 0u &&
       (trailing_cols == NULL || trailing_panel == NULL))) {
    return 0;
  }
  if (update_work < KLS_ROW_REFACTOR_BATCH_SUPERNODE_MIN_WORK ||
      update_work <
        KLS_ROW_REFACTOR_BATCH_SUPERNODE_MIN_WORK_PER_ENTRY *
          copied_entries) {
    kls_record_row_refactor_compact_supernode_batch_candidate(
      solver, batch_end - batch_begin, candidate_dep_rows, 1);
    return 0;
  }
  kls_record_row_refactor_compact_supernode_batch_candidate(
    solver, batch_end - batch_begin, candidate_dep_rows, 0);

  const UF_long dep_group_begin = solver->row_refactor_group_ptr[dep_group];
  const UF_long dep_group_end = solver->row_refactor_group_ptr[dep_group + 1u];
  if (dep_group_begin >= dep_group_end || dep_group_end > row_begin) {
    return 0;
  }
  const UF_long dep_width = dep_group_end - dep_group_begin;
  if (dep_width == 0u ||
      kls_row_refactor_compact_panel_valid_prefix(solver, dep_group) <
        dep_width ||
      dep_width > (UF_long)INT_MAX) {
    return 0;
  }
  const UF_long dep_trailing_len =
    solver->row_refactor_group_trailing_len[dep_group];
  if (dep_trailing_len > (UF_long)INT_MAX) {
    return 0;
  }
  UF_long dep_dense_entries = 0;
  if (dep_width > (UF_long)(SIZE_MAX / sizeof(double)) / dep_width) {
    kls_egraph_refactor_record_invalid(shared);
    return -1;
  }
  dep_dense_entries = dep_width * dep_width;
  if (dep_trailing_len > 0u &&
      dep_width > ((UF_long)(SIZE_MAX / sizeof(double)) -
                     dep_dense_entries) / dep_trailing_len) {
    kls_egraph_refactor_record_invalid(shared);
    return -1;
  }
  const UF_long dep_panel_entries =
    dep_dense_entries + dep_width * dep_trailing_len;
  double *dep_dense_panel =
    kls_row_refactor_compact_group_panel(solver, dep_group,
                                         dep_panel_entries);
  if (dep_dense_panel == NULL) {
    return 0;
  }
  double *dep_trailing_panel = dep_trailing_len > 0u
    ? dep_dense_panel + dep_dense_entries : NULL;
  const UF_long *dep_trailing_cols = NULL;
  if (!kls_row_refactor_group_trailing_slice(
        solver, dep_group, dep_group_end, NULL,
        &dep_trailing_cols)) {
    kls_egraph_refactor_record_invalid(shared);
    return -1;
  }
  if (dep_trailing_len > 0u &&
      (dep_trailing_panel == NULL || dep_trailing_cols == NULL)) {
    kls_egraph_refactor_record_invalid(shared);
    return -1;
  }

  const UF_long batch_rows = batch_end - batch_begin;
  if (batch_rows == 0u ||
      batch_rows > (UF_long)(SIZE_MAX / sizeof(UF_long)) ||
      batch_rows > UF_long_max / 2u ||
      batch_rows > (UF_long)(SIZE_MAX / sizeof(*worker->supernode_workspace)) /
                     dep_width) {
    return 0;
  }
  UF_long *index_workspace =
    kls_egraph_worker_index_workspace(worker, 2u * batch_rows);
  if (index_workspace == NULL) {
    return 0;
  }
  UF_long *suffix_begins = index_workspace;
  UF_long *dep_rows_by_row = suffix_begins + batch_rows;

  int status = 0;
  double *workspace =
    kls_egraph_worker_supernode_workspace(worker, batch_rows * dep_width);
  if (workspace == NULL) {
    goto cleanup;
  }
  memset(workspace, 0,
         (size_t)batch_rows * (size_t)dep_width * sizeof(*workspace));

  double *x = worker->x;
  double *udiag = (double *)solver->numeric->Udiag;
  const UF_long width = row_end - row_begin;
  for (UF_long local_dep = 0; local_dep < dep_width; ++local_dep) {
    dep_dense_panel[local_dep * dep_width + local_dep] =
      udiag[dep_group_begin + local_dep];
  }

  for (UF_long batch_local = 0; batch_local < batch_rows; ++batch_local) {
    const UF_long row = batch_begin + batch_local;
    UF_long row_dep_group = KLS_KLU_EMPTY;
    UF_long suffix_begin = 0;
    UF_long row_dep_rows = 0;
    if (!kls_compact_dense_group_single_producer_suffix(
          solver, group, row_begin, row, &row_dep_group,
          &suffix_begin, &row_dep_rows) ||
        row_dep_group != dep_group ||
        row_dep_rows != dep_width - suffix_begin) {
      goto cleanup;
    }
    suffix_begins[batch_local] = suffix_begin;
    dep_rows_by_row[batch_local] = row_dep_rows;

    const UF_long local_row = row - row_begin;
    double *row_dense_panel = dense_panel + local_row * width;
    double *row_panel = trailing_panel != NULL
      ? trailing_panel + local_row * trailing_len : NULL;
    const int direct_input =
      kls_parallel_row_refactor_load_compact_dense_input_row(
        shared, x, row_begin, row_end, trailing_len, trailing_cols,
        dense_panel, trailing_panel, row);
    if (direct_input < 0) {
      status = -1;
      goto cleanup;
    }
    if (!direct_input &&
        !kls_parallel_row_refactor_load_input_row(shared, x, row)) {
      status = -1;
      goto cleanup;
    }

    const UF_long l_begin = solver->row_refactor_l_ptr[row];
    const UF_long l_internal = solver->row_refactor_l_internal_ptr[row];
    const UF_long l_end = solver->row_refactor_l_ptr[row + 1u];
    for (UF_long local = 0; local < row_dep_rows; ++local) {
      const UF_long p = l_begin + local;
      const UF_long dep = solver->row_refactor_l_cols[p];
      const UF_long dep_local = suffix_begin + local;
      if (dep != dep_group_begin + dep_local) {
        kls_egraph_refactor_record_invalid(shared);
        status = -1;
        goto cleanup;
      }
      workspace[batch_local * dep_width + dep_local] = x[dep];
      x[dep] = 0.0;
      solver->row_refactor_l_row_values[p] = 0.0;
    }

    UF_long lp = l_internal;
    for (UF_long dep = row_begin; dep < row; ++dep, ++lp) {
      if (lp >= l_end || solver->row_refactor_l_cols[lp] != dep) {
        kls_egraph_refactor_record_invalid(shared);
        status = -1;
        goto cleanup;
      }
      if (direct_input) {
        row_dense_panel[dep - row_begin] += x[dep];
      } else {
        row_dense_panel[dep - row_begin] = x[dep];
      }
      x[dep] = 0.0;
    }
    if (lp != l_end) {
      kls_egraph_refactor_record_invalid(shared);
      status = -1;
      goto cleanup;
    }

    if (direct_input) {
      udiag[row] = row_dense_panel[local_row] + x[row];
    } else {
      udiag[row] = x[row];
    }
    x[row] = 0.0;

    const UF_long dense_len = row_end - row - 1u;
    const UF_long u_begin = solver->row_refactor_u_ptr[row];
    const UF_long u_end = solver->row_refactor_u_ptr[row + 1u];
    if (u_end - u_begin != dense_len + trailing_len) {
      kls_egraph_refactor_record_invalid(shared);
      status = -1;
      goto cleanup;
    }
    for (UF_long offset = 0; offset < dense_len; ++offset) {
      const UF_long col = row + 1u + offset;
      if (direct_input) {
        row_dense_panel[local_row + 1u + offset] += x[col];
      } else {
        row_dense_panel[local_row + 1u + offset] = x[col];
      }
      x[col] = 0.0;
    }
    if (trailing_len > 0u) {
      for (UF_long offset = 0; offset < trailing_len; ++offset) {
        const UF_long col = trailing_cols[offset];
        if (direct_input) {
          row_panel[offset] += x[col];
        } else {
          row_panel[offset] = x[col];
        }
        x[col] = 0.0;
      }
    }
    kls_clear_row_refactor_input_residuals(solver, x, row);
  }

  double trsv_work = 0.0;
  for (UF_long batch_local = 0; batch_local < batch_rows; ++batch_local) {
    double *row_multipliers = workspace + batch_local * dep_width;
    const UF_long suffix_begin = suffix_begins[batch_local];
    const UF_long row = batch_begin + batch_local;
    const UF_long l_begin = solver->row_refactor_l_ptr[row];
    for (UF_long local_dep = suffix_begin; local_dep < dep_width;
         ++local_dep) {
      const UF_long dep = dep_group_begin + local_dep;
      const double candidate = row_multipliers[local_dep];
      const double lij = candidate / udiag[dep];
      if (kls_parallel_row_refactor_rejects_multiplier(
            worker, row, dep, candidate, lij)) {
        status = -1;
        goto cleanup;
      }
      row_multipliers[local_dep] = lij;
      for (UF_long target = local_dep + 1u; target < dep_width; ++target) {
        row_multipliers[target] -=
          lij * dep_dense_panel[local_dep * dep_width + target];
      }
      const UF_long p = l_begin + local_dep - suffix_begin;
      if (p >= solver->row_refactor_l_internal_ptr[row] ||
          solver->row_refactor_l_cols[p] != dep) {
        kls_egraph_refactor_record_invalid(shared);
        status = -1;
        goto cleanup;
      }
      solver->row_refactor_l_row_values[p] = lij;
      if (!shared->row_refactor_defer_value_scatter) {
        *solver->row_refactor_l_values[p] = lij;
      }
      trsv_work += (double)(dep_width - local_dep - 1u);
    }
  }

  double gemv_work = 0.0;
  if (dep_trailing_len > 0u) {
    for (UF_long dep_col_pos = 0; dep_col_pos < dep_trailing_len;
         ++dep_col_pos) {
      const UF_long col = dep_trailing_cols[dep_col_pos];
      UF_long trailing_pos = KLS_KLU_EMPTY;
      const int in_current_group = col >= row_begin && col < row_end;
      const int in_current_trailing =
        !in_current_group &&
        kls_find_uflong_sorted(trailing_cols, trailing_len, col,
                               &trailing_pos);
      if (!in_current_group && !in_current_trailing) {
        continue;
      }
      for (UF_long batch_local = 0; batch_local < batch_rows; ++batch_local) {
        const UF_long row = batch_begin + batch_local;
        const UF_long local_row = row - row_begin;
        const UF_long suffix_begin = suffix_begins[batch_local];
        const double *row_multipliers = workspace + batch_local * dep_width;
        double update = 0.0;
        for (UF_long local_dep = suffix_begin; local_dep < dep_width;
             ++local_dep) {
          update += row_multipliers[local_dep] *
            dep_trailing_panel[local_dep * dep_trailing_len + dep_col_pos];
        }
        if (update == 0.0) {
          continue;
        }
        if (in_current_group) {
          const UF_long local_col = col - row_begin;
          if (col == row) {
            udiag[row] -= update;
          } else {
            dense_panel[local_row * width + local_col] -= update;
          }
        } else {
          trailing_panel[local_row * trailing_len + trailing_pos] -= update;
        }
      }
    }
    for (UF_long batch_local = 0; batch_local < batch_rows; ++batch_local) {
      gemv_work += (double)dep_rows_by_row[batch_local] *
        (double)dep_trailing_len;
    }
  }

  const UF_long trsv_entries = kls_work_to_uflong(trsv_work);
  const UF_long gemv_entries = kls_work_to_uflong(gemv_work);
  const UF_long total_entries =
    trsv_entries > UF_long_max - gemv_entries
      ? UF_long_max : trsv_entries + gemv_entries;
  kls_record_row_refactor_compact_supernode_batch(
    solver, batch_rows, candidate_dep_rows, total_entries);
  kls_record_row_refactor_compact_supernode_update(
    solver, candidate_dep_rows, total_entries);
  kls_record_row_refactor_compact_supernode_trsv(
    solver, candidate_dep_rows, trsv_entries);
  if (gemv_entries > 0u) {
    kls_record_row_refactor_compact_supernode_gemv(
      solver, candidate_dep_rows, gemv_entries);
  }
  status = 1;

cleanup:
  return status;
}

static int kls_parallel_row_refactor_process_dense_group_compact_checked(
  kls_egraph_refactor_worker *worker,
  UF_long group,
  UF_long row_begin,
  UF_long row_end,
  UF_long trailing_len,
  double *dense_panel,
  double *trailing_panel,
  int publish_prefix) {
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  trilinos_klu_l_numeric *numeric = solver->numeric;
  const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  double *x = worker->x;
  double *udiag = (double *)numeric->Udiag;
  const UF_long width = row_end - row_begin;

  const int blocked_panel_status =
    kls_compact_dense_panel_factor_blocked_checked(
      worker, group, row_begin, row_end, trailing_len, dense_panel,
      trailing_panel, publish_prefix);
  if (blocked_panel_status > 0) {
    return 1;
  }
  if (blocked_panel_status < 0) {
    return 0;
  }

  for (UF_long row = row_begin; row < row_end; ++row) {
    const UF_long local_row = row - row_begin;
    double *row_dense_panel = dense_panel + local_row * width;
    double *row_panel = trailing_len > 0u
      ? trailing_panel + local_row * trailing_len : NULL;
    const UF_long row_dense_len = row_end - row - 1u;

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
        if (trailing_len > 0u) {
          double *dep_panel = trailing_panel + local_dep * trailing_len;
          for (UF_long offset = 0; offset < trailing_len; ++offset) {
            row_panel[offset] -= lij * dep_panel[offset];
          }
        }
      }

    if (!kls_store_compact_dense_panel_row(shared, row_begin, row_end,
                                          trailing_len, dense_panel,
                                          trailing_panel, row, 0)) {
      kls_egraph_refactor_record_invalid(shared);
      return 0;
    }

    double row_max_abs =
      kls_row_refactor_abs_max_from_values(
        udiag[row], row_dense_panel + local_row + 1u, row_dense_len);
    if (trailing_len > 0u) {
      row_max_abs =
        kls_row_refactor_abs_max_from_values(row_max_abs, row_panel,
                                             trailing_len);
    }
    if (kls_parallel_row_refactor_rejects_pivot(worker, row, udiag[row],
                                               row_max_abs)) {
      memset(x, 0, (size_t)solver->n * sizeof(*x));
      return 0;
    }
    if (udiag[row] == 0.0) {
      kls_egraph_refactor_record_singular(shared, row, symbolic->Q[row]);
      if (solver->common.halt_if_singular) {
        memset(x, 0, (size_t)solver->n * sizeof(*x));
        return 0;
      }
    }
    if (!shared->row_refactor_defer_value_scatter &&
        !kls_scatter_row_refactor_row_values(solver, row)) {
      kls_egraph_refactor_record_invalid(shared);
      return 0;
    }
    if (publish_prefix) {
      kls_mark_row_refactor_compact_panel_prefix_valid(
        solver, group, row - row_begin + 1u);
    }
    kls_row_refactor_mark_row_done(worker, row);
  }
  return 1;
}

static int kls_row_refactor_assemble_compact_dense_row(
  kls_egraph_refactor_worker *worker,
  UF_long group,
  UF_long row_begin,
  UF_long row_end,
  UF_long trailing_len,
  const UF_long *trailing_cols,
  double *dense_panel,
  double *trailing_panel,
  int wait_for_dependencies,
  UF_long row) {
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  double *x = worker->x;
  double *udiag = (double *)solver->numeric->Udiag;
  const UF_long width = row_end - row_begin;
  const UF_long local_row = row - row_begin;
  double *row_dense_panel = dense_panel + local_row * width;
  double *row_panel = trailing_panel != NULL
    ? trailing_panel + local_row * trailing_len : NULL;
  const int direct_input =
    kls_parallel_row_refactor_load_compact_dense_input_row(
      shared, x, row_begin, row_end, trailing_len, trailing_cols,
      dense_panel, trailing_panel, row);
  if (direct_input < 0) {
    return 0;
  }
  if (!direct_input &&
      !kls_parallel_row_refactor_load_input_row(shared, x, row)) {
    return 0;
  }

  const UF_long l_begin = solver->row_refactor_l_ptr[row];
  UF_long lp = l_begin;
  const UF_long l_internal = solver->row_refactor_l_internal_ptr[row];
  const UF_long l_end = solver->row_refactor_l_ptr[row + 1u];
  while (lp < l_internal) {
    const int compact_status =
      kls_row_refactor_try_compact_supernode_update(
        worker, row, group, &lp, l_begin, l_internal,
        wait_for_dependencies, NULL);
    if (compact_status < 0) {
      x[row] = 0.0;
      return 0;
    }
    if (compact_status > 0) {
      continue;
    }
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
    if (!shared->row_refactor_defer_value_scatter) {
      *solver->row_refactor_l_values[lp] = lij;
    }
    x[dep] = 0.0;
    if (!kls_row_refactor_try_compact_panel_scalar_update(
          shared, x, group, dep, lij)) {
      const UF_long u_begin = solver->row_refactor_u_ptr[dep];
      const UF_long u_end = solver->row_refactor_u_ptr[dep + 1u];
      const UF_long *u_cols = solver->row_refactor_u_cols + u_begin;
      const double *u_values =
        solver->row_refactor_u_row_values + u_begin;
      for (UF_long offset = 0; offset < u_end - u_begin; ++offset) {
        x[u_cols[offset]] -= lij * u_values[offset];

      }
    }
    lp++;
  }

  for (UF_long dep = row_begin; dep < row; ++dep, ++lp) {
    if (lp >= l_end || solver->row_refactor_l_cols[lp] != dep) {
      kls_egraph_refactor_record_invalid(shared);
      return 0;
    }
    if (direct_input) {
      row_dense_panel[dep - row_begin] += x[dep];
    } else {
      row_dense_panel[dep - row_begin] = x[dep];
    }
    x[dep] = 0.0;
  }
  if (lp != l_end) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }

  udiag[row] =
    direct_input ? row_dense_panel[local_row] + x[row] : x[row];
  x[row] = 0.0;

  const UF_long dense_len = row_end - row - 1u;
  const UF_long u_begin = solver->row_refactor_u_ptr[row];
  const UF_long u_end = solver->row_refactor_u_ptr[row + 1u];
  if (u_end - u_begin != dense_len + trailing_len) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }
  for (UF_long offset = 0; offset < dense_len; ++offset) {
    const UF_long col = row + 1u + offset;
    if (direct_input) {
      row_dense_panel[local_row + 1u + offset] += x[col];
    } else {
      row_dense_panel[local_row + 1u + offset] = x[col];
    }
    x[col] = 0.0;
  }
  for (UF_long offset = 0; offset < trailing_len; ++offset) {
    const UF_long col = trailing_cols[offset];
    if (direct_input) {
      row_panel[offset] += x[col];
    } else {
      row_panel[offset] = x[col];
    }
    x[col] = 0.0;
  }
  kls_clear_row_refactor_input_residuals(solver, x, row);
  return 1;
}

static int kls_parallel_row_refactor_process_dense_group_compact(
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
      solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_l_cols == NULL ||
      solver->row_refactor_l_row_values == NULL ||
      solver->row_refactor_l_values == NULL ||
      solver->row_refactor_u_ptr == NULL ||
      solver->row_refactor_u_row_values == NULL ||
      solver->row_refactor_u_values == NULL ||
      !solver->row_refactor_group_dense[group] ||
      solver->row_refactor_group_trailing_len[group] != trailing_len) {
    return -1;
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
  const UF_long panel_entries = dense_panel_entries + trailing_panel_entries;
  double *dense_panel =
    kls_row_refactor_compact_group_panel(solver, group, panel_entries);
  const int persistent_panel = dense_panel != NULL;
  if (dense_panel == NULL) {
    dense_panel = kls_egraph_worker_segment_panel(worker, panel_entries);
  }
  if (dense_panel == NULL) {
    return -1;
  }
  kls_record_row_refactor_compact_dense_panel(solver, persistent_panel);
  double *trailing_panel = NULL;
  if (trailing_len > 0u) {
    trailing_panel = dense_panel + dense_panel_entries;
  }

  {
    int coop = 0;
    if (shared->thread_count > 1 && width >= 8) {
      int expected = 0;
      if (atomic_compare_exchange_strong_explicit(
            &shared->dense_help_owner_lock, &expected, 1,
            memory_order_acq_rel, memory_order_relaxed)) {
        coop = 1;
      }
    }
    if (coop) {
      struct kls_dense_help_ctx help_ctx;
      help_ctx.kind = KLS_DENSE_HELP_ASSEMBLE;
      help_ctx.group = group;
      help_ctx.row_begin = row_begin;
      help_ctx.row_end = row_end;
      help_ctx.trailing_len = trailing_len;
      help_ctx.trailing_cols = trailing_cols;
      help_ctx.dense_panel = dense_panel;
      help_ctx.trailing_panel = trailing_panel;
      help_ctx.wait_for_dependencies = wait_for_dependencies;
      help_ctx.panel_k = 0;
      help_ctx.panel_kb = 0;
      const long assemble_epoch = ++shared->dense_help_epoch_seq;
      help_ctx.epoch = assemble_epoch;
      atomic_store_explicit(&shared->dense_help_ctx, &help_ctx,
                            memory_order_relaxed);
      atomic_store_explicit(&shared->dense_help_cursor,
                            (assemble_epoch << 32) | (long)row_begin,
                            memory_order_release);
      atomic_store_explicit(&shared->dense_help_done, 0,
                            memory_order_relaxed);
      atomic_store_explicit(&shared->dense_help_fail, 0,
                            memory_order_relaxed);
      atomic_store_explicit(&shared->dense_help_active, 1,
                            memory_order_release);
      for (;;) {
        if (atomic_load_explicit(&shared->dense_help_fail,
                                 memory_order_acquire)) {
          break;
        }
        const long claim =
          kls_dense_help_claim(shared, assemble_epoch, (long)row_end);
        if (claim < 0) {
          break;
        }
        const int ok = kls_row_refactor_assemble_compact_dense_row(
          worker, group, row_begin, row_end, trailing_len, trailing_cols,
          dense_panel, trailing_panel, wait_for_dependencies,
          (UF_long)claim);
        atomic_fetch_add_explicit(&shared->dense_help_done, 1,
                                  memory_order_release);
        if (!ok) {
          atomic_store_explicit(&shared->dense_help_fail, 1,
                                memory_order_release);
          break;
        }
      }
      atomic_store_explicit(&shared->dense_help_active, 0,
                            memory_order_release);
      for (;;) {
        const long done = atomic_load_explicit(&shared->dense_help_done,
                                               memory_order_acquire);
        const int failed = atomic_load_explicit(&shared->dense_help_fail,
                                                memory_order_acquire);
        if (atomic_load_explicit(&shared->dense_help_inflight,
                                 memory_order_acquire) == 0 &&
            (failed || done >= (long)width)) {
          break;
        }
        kls_cpu_relax();
      }
      atomic_store_explicit(&shared->dense_help_ctx, NULL,
                            memory_order_release);
      const int final_failed = atomic_load_explicit(
        &shared->dense_help_fail, memory_order_relaxed);
      atomic_store_explicit(&shared->dense_help_owner_lock, 0,
                            memory_order_release);
      if (final_failed) {
        return 0;
      }
    } else {
      UF_long row = row_begin;
      while (row < row_end) {
        const UF_long batch_end =
          kls_compact_dense_group_external_pattern_end(
            solver, row_begin, row_end, row);
        if (batch_end > row + 1u) {
          kls_record_row_refactor_compact_supernode_batch_pattern(
            solver, batch_end - row);
          const int batched_status =
            kls_compact_dense_group_try_batched_supernode_update(
              worker, group, row_begin, row_end, row, batch_end, trailing_len,
              trailing_cols, dense_panel, trailing_panel);
          if (batched_status < 0) {
            return 0;
          }
          if (batched_status > 0) {
            row = batch_end;
            continue;
          }
        }
        {
          UF_long ragged_dep_group = KLS_KLU_EMPTY;
          double ragged_work = 0.0;
          double ragged_copied_entries = 0.0;
          UF_long ragged_dep_rows = 0;
          const UF_long ragged_end =
            kls_compact_dense_group_ragged_supernode_end(
              solver, group, row_begin, row_end, row, &ragged_dep_group,
              &ragged_work, &ragged_copied_entries, &ragged_dep_rows);
          if (ragged_end > row + 1u) {
            const int ragged_status =
              kls_compact_dense_group_try_ragged_supernode_update(
                worker, group, row_begin, row_end, row, ragged_end,
                ragged_dep_group, ragged_work, ragged_copied_entries,
                ragged_dep_rows, trailing_len, trailing_cols, dense_panel,
                trailing_panel);
            if (ragged_status < 0) {
              return 0;
            }
            if (ragged_status > 0) {
              row = ragged_end;
              continue;
            }
          }
        }
        if (!kls_row_refactor_assemble_compact_dense_row(
              worker, group, row_begin, row_end, trailing_len, trailing_cols,
              dense_panel, trailing_panel, wait_for_dependencies, row)) {
          return 0;
        }
        row++;
      }
    }
  }

  if (shared->check_pivots) {
    const int ok = kls_parallel_row_refactor_process_dense_group_compact_checked(
      worker, group, row_begin, row_end, trailing_len, dense_panel,
      trailing_panel, persistent_panel);
    if (ok && persistent_panel) {
      kls_mark_row_refactor_compact_panel_valid(solver, group);
    }
    return ok;
  }

  const int blocked_panel_status =
    kls_compact_dense_panel_factor_blocked_unchecked(
      worker, group, row_begin, row_end, trailing_len, dense_panel,
      trailing_panel, persistent_panel);
  if (blocked_panel_status > 0) {
    if (persistent_panel) {
      kls_mark_row_refactor_compact_panel_valid(solver, group);
    }
    return 1;
  }
  if (blocked_panel_status < 0) {
    return 0;
  }

  for (UF_long row = row_begin; row < row_end; ++row) {
    const UF_long local_row = row - row_begin;
    double *row_dense_panel = dense_panel + local_row * width;
    double *row_panel = trailing_panel != NULL
      ? trailing_panel + local_row * trailing_len : NULL;

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
        if (trailing_len > 0u) {
          double *dep_panel = trailing_panel + local_dep * trailing_len;
          for (UF_long offset = 0; offset < trailing_len; ++offset) {
            row_panel[offset] -= lij * dep_panel[offset];
          }
        }
      }

    if (udiag[row] == 0.0) {
      kls_egraph_refactor_record_singular(shared, row, symbolic->Q[row]);
      if (solver->common.halt_if_singular) {
        memset(x, 0, (size_t)solver->n * sizeof(*x));
        return 0;
      }
    }

    if (!kls_store_compact_dense_panel_row(shared, row_begin, row_end,
                                          trailing_len, dense_panel,
                                          trailing_panel, row, 1)) {
      kls_egraph_refactor_record_invalid(shared);
      return 0;
    }
    if (persistent_panel) {
      kls_mark_row_refactor_compact_panel_prefix_valid(
        solver, group, row - row_begin + 1u);
    }
    kls_row_refactor_mark_row_done(worker, row);
  }
  if (persistent_panel) {
    kls_mark_row_refactor_compact_panel_valid(solver, group);
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
  const UF_long width = row_end - row_begin;
  kls_egraph_refactor_shared *shared =
    worker != NULL ? worker->shared : NULL;
  kls_solver *solver = shared != NULL ? shared->solver : NULL;
  const int native_env_state = shared != NULL
    ? shared->row_refactor_native_row_panel_state
    : KLS_NATIVE_ROW_PANEL_ENV_OFF;
  const int native_row_panel =
    shared != NULL &&
    kls_native_row_panel_should_try_for_solver(
      solver, width, native_env_state, shared->check_pivots);
  int compact_attempted = 0;

  if (native_row_panel) {
    const UF_long panel_entries = kls_work_to_uflong(
      kls_row_refactor_dense_group_panel_entries(width, trailing_len));
    kls_record_row_refactor_native_row_panel(solver, width, panel_entries);
    shared->row_refactor_native_row_panel_active = 1;
    compact_attempted = 1;
    const int compact_status =
      kls_parallel_row_refactor_process_dense_group_compact(
        worker, group, row_begin, row_end, trailing_len, trailing_cols,
        wait_for_dependencies);
    shared->row_refactor_native_row_panel_active = 0;
    if (compact_status >= 0) {
      if (compact_status == 0 && shared->pivot_rejected) {
        kls_record_row_refactor_native_row_panel_checked_reject(solver);
      }
      return compact_status;
    }
    kls_record_row_refactor_native_row_panel_fallback(solver);
  }

  const int native_status =
    kls_parallel_row_refactor_process_dense_group_native(
      worker, group, row_begin, row_end, trailing_len, trailing_cols,
      wait_for_dependencies);
  if (native_status >= 0) {
    return native_status;
  }

  if (!compact_attempted && native_env_state >= 0) {
    return kls_parallel_row_refactor_process_dense_group_compact(
      worker, group, row_begin, row_end, trailing_len, trailing_cols,
      wait_for_dependencies);
  }
  return -1;
}

static int kls_batch_group_single_producer_suffix(
  const kls_solver *solver,
  UF_long current_group,
  UF_long row,
  UF_long *dep_group_out,
  UF_long *suffix_begin_out,
  UF_long *dep_rows_out) {
  if (dep_group_out != NULL) {
    *dep_group_out = KLS_KLU_EMPTY;
  }
  KLS_SET_OPTIONAL_OUTPUT(suffix_begin_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(dep_rows_out, 0);
  if (solver == NULL ||
      row >= solver->n ||
      solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_l_cols == NULL ||
      solver->row_refactor_row_group == NULL ||
      solver->row_refactor_group_ptr == NULL ||
      solver->row_refactor_group_kind == NULL) {
    return 0;
  }
  const UF_long l_begin = solver->row_refactor_l_ptr[row];
  const UF_long l_end = solver->row_refactor_l_ptr[row + 1u];
  if (l_begin >= l_end) {
    return 0;
  }
  const UF_long dep0 = solver->row_refactor_l_cols[l_begin];
  if (dep0 >= row || dep0 >= solver->n) {
    return 0;
  }
  const UF_long dep_group = solver->row_refactor_row_group[dep0];
  if (dep_group >= solver->row_refactor_group_count ||
      dep_group == current_group ||
      solver->row_refactor_group_kind[dep_group] !=
        KLS_ROW_REFACTOR_GROUP_DENSE) {
    return 0;
  }
  const UF_long group_begin = solver->row_refactor_group_ptr[dep_group];
  const UF_long group_end = solver->row_refactor_group_ptr[dep_group + 1u];
  if (group_begin >= group_end ||
      dep0 < group_begin ||
      dep0 >= group_end ||
      group_end > row) {
    return 0;
  }
  const UF_long dep_width = group_end - group_begin;
  if (kls_row_refactor_compact_panel_valid_prefix(solver, dep_group) <
      dep_width) {
    return 0;
  }

  UF_long p = l_begin;
  UF_long expected = dep0;
  while (p < l_end && expected < group_end) {
    const UF_long dep = solver->row_refactor_l_cols[p];
    if (dep != expected ||
        dep >= row ||
        solver->row_refactor_row_group[dep] != dep_group) {
      return 0;
    }
    p++;
    expected++;
  }
  if (p != l_end || expected != group_end) {
    return 0;
  }
  if (dep_group_out != NULL) {
    *dep_group_out = dep_group;
  }
  if (suffix_begin_out != NULL) {
    *suffix_begin_out = dep0 - group_begin;
  }
  if (dep_rows_out != NULL) {
    *dep_rows_out = l_end - l_begin;
  }
  return 1;
}

static UF_long kls_batch_group_ragged_supernode_end(
  const kls_solver *solver,
  UF_long group,
  UF_long row_begin,
  UF_long row_end,
  UF_long batch_begin,
  UF_long *dep_group_out,
  double *work_out,
  double *copied_entries_out,
  UF_long *dep_rows_out) {
  if (dep_group_out != NULL) {
    *dep_group_out = KLS_KLU_EMPTY;
  }
  if (work_out != NULL) {
    *work_out = 0.0;
  }
  if (copied_entries_out != NULL) {
    *copied_entries_out = 0.0;
  }
  KLS_SET_OPTIONAL_OUTPUT(dep_rows_out, 0);
  if (solver == NULL ||
      solver->row_refactor_group_trailing_len == NULL ||
      batch_begin < row_begin ||
      batch_begin >= row_end) {
    return batch_begin;
  }

  UF_long first_dep_group = KLS_KLU_EMPTY;
  UF_long dep_rows = 0;
  if (!kls_batch_group_single_producer_suffix(
        solver, group, batch_begin, &first_dep_group, NULL,
        &dep_rows)) {
    return batch_begin + 1u;
  }

  const UF_long dep_trailing_len =
    solver->row_refactor_group_trailing_len[first_dep_group];
  double work = 0.5 * (double)dep_rows * (double)(dep_rows - 1u) +
    (double)dep_rows * (double)dep_trailing_len;
  double copied = (double)dep_rows + (double)dep_trailing_len;
  double total_dep_rows = (double)dep_rows;
  UF_long batch_end = batch_begin + 1u;
  while (batch_end < row_end) {
    UF_long dep_group = KLS_KLU_EMPTY;
    UF_long row_dep_rows = 0;
    if (!kls_batch_group_single_producer_suffix(
          solver, group, batch_end, &dep_group, NULL,
          &row_dep_rows) ||
        dep_group != first_dep_group) {
      break;
    }
    work += 0.5 * (double)row_dep_rows * (double)(row_dep_rows - 1u) +
      (double)row_dep_rows * (double)dep_trailing_len;
    copied += (double)row_dep_rows + (double)dep_trailing_len;
    total_dep_rows += (double)row_dep_rows;
    batch_end++;
  }
  if (batch_end <= batch_begin + 1u) {
    return batch_begin + 1u;
  }
  if (dep_group_out != NULL) {
    *dep_group_out = first_dep_group;
  }
  if (work_out != NULL) {
    *work_out = work;
  }
  if (copied_entries_out != NULL) {
    *copied_entries_out = copied;
  }
  if (dep_rows_out != NULL) {
    *dep_rows_out = kls_work_to_uflong(total_dep_rows);
  }
  return batch_end;
}

static int kls_batch_group_try_ragged_supernode_update(
  kls_egraph_refactor_worker *worker,
  UF_long group,
  UF_long batch_begin,
  UF_long batch_end,
  UF_long dep_group,
  double update_work,
  double copied_entries,
  UF_long candidate_dep_rows) {
  if (worker == NULL || worker->shared == NULL ||
      batch_begin >= batch_end ||
      candidate_dep_rows == 0u) {
    return 0;
  }
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  if (solver == NULL ||
      dep_group >= solver->row_refactor_group_count ||
      solver->symbolic == NULL ||
      solver->symbolic->Q == NULL ||
      solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_l_cols == NULL ||
      solver->row_refactor_l_row_values == NULL ||
      solver->row_refactor_l_values == NULL ||
      solver->row_refactor_u_ptr == NULL ||
      solver->row_refactor_u_cols == NULL ||
      solver->row_refactor_u_row_values == NULL ||
      solver->row_refactor_u_values == NULL ||
      solver->row_refactor_input_ptr == NULL ||
      solver->row_refactor_input_cols == NULL ||
      solver->row_refactor_row_group == NULL ||
      solver->row_refactor_group_ptr == NULL ||
      solver->row_refactor_group_kind == NULL ||
      solver->row_refactor_group_trailing_len == NULL ||
      solver->numeric == NULL ||
      solver->numeric->Udiag == NULL ||
      worker->x == NULL) {
    return 0;
  }
  const UF_long batch_rows = batch_end - batch_begin;
  if (batch_rows == 0u ||
      batch_rows > (UF_long)(SIZE_MAX / sizeof(UF_long))) {
    return 0;
  }
  if (update_work < KLS_ROW_REFACTOR_BATCH_SUPERNODE_MIN_WORK ||
      update_work <
        KLS_ROW_REFACTOR_BATCH_SUPERNODE_MIN_WORK_PER_ENTRY *
          copied_entries) {
    kls_record_row_refactor_compact_supernode_batch_candidate(
      solver, batch_rows, candidate_dep_rows, 1);
    return 0;
  }
  kls_record_row_refactor_compact_supernode_batch_candidate(
    solver, batch_rows, candidate_dep_rows, 0);

  const UF_long dep_group_begin = solver->row_refactor_group_ptr[dep_group];
  const UF_long dep_group_end = solver->row_refactor_group_ptr[dep_group + 1u];
  if (dep_group_begin >= dep_group_end) {
    return 0;
  }
  const UF_long dep_width = dep_group_end - dep_group_begin;
  if (dep_width == 0u ||
      kls_row_refactor_compact_panel_valid_prefix(solver, dep_group) <
        dep_width ||
      dep_width > (UF_long)INT_MAX ||
      batch_rows > (UF_long)(SIZE_MAX / sizeof(*worker->supernode_workspace)) /
                     dep_width) {
    return 0;
  }
  const UF_long dep_trailing_len =
    solver->row_refactor_group_trailing_len[dep_group];
  UF_long dep_dense_entries = 0;
  if (dep_width > (UF_long)(SIZE_MAX / sizeof(double)) / dep_width) {
    kls_egraph_refactor_record_invalid(shared);
    return -1;
  }
  dep_dense_entries = dep_width * dep_width;
  if (dep_trailing_len > 0u &&
      dep_width > ((UF_long)(SIZE_MAX / sizeof(double)) -
                     dep_dense_entries) / dep_trailing_len) {
    kls_egraph_refactor_record_invalid(shared);
    return -1;
  }
  const UF_long dep_panel_entries =
    dep_dense_entries + dep_width * dep_trailing_len;
  double *dep_dense_panel =
    kls_row_refactor_compact_group_panel(solver, dep_group,
                                         dep_panel_entries);
  if (dep_dense_panel == NULL) {
    return 0;
  }
  double *dep_trailing_panel = dep_trailing_len > 0u
    ? dep_dense_panel + dep_dense_entries : NULL;
  const UF_long *dep_trailing_cols = NULL;
  if (!kls_row_refactor_group_trailing_slice(
        solver, dep_group, dep_group_end, NULL,
        &dep_trailing_cols)) {
    kls_egraph_refactor_record_invalid(shared);
    return -1;
  }
  if (dep_trailing_len > 0u &&
      (dep_trailing_panel == NULL || dep_trailing_cols == NULL)) {
    kls_egraph_refactor_record_invalid(shared);
    return -1;
  }

  UF_long *suffix_begins =
    (UF_long *)malloc((size_t)batch_rows * sizeof(*suffix_begins));
  UF_long *dep_rows_by_row =
    (UF_long *)malloc((size_t)batch_rows * sizeof(*dep_rows_by_row));
  UF_long *u_offsets =
    (UF_long *)calloc((size_t)batch_rows + 1u, sizeof(*u_offsets));
  double *pivots =
    (double *)malloc((size_t)batch_rows * sizeof(*pivots));
  if (suffix_begins == NULL || dep_rows_by_row == NULL ||
      u_offsets == NULL || pivots == NULL) {
    free(suffix_begins);
    free(dep_rows_by_row);
    free(u_offsets);
    free(pivots);
    return 0;
  }

  int status = 0;
  double *u_workspace = NULL;
  for (UF_long batch_local = 0; batch_local < batch_rows; ++batch_local) {
    const UF_long row = batch_begin + batch_local;
    UF_long row_dep_group = KLS_KLU_EMPTY;
    UF_long suffix_begin = 0;
    UF_long row_dep_rows = 0;
    if (!kls_batch_group_single_producer_suffix(
          solver, group, row, &row_dep_group, &suffix_begin,
          &row_dep_rows) ||
        row_dep_group != dep_group ||
        row_dep_rows != dep_width - suffix_begin) {
      goto cleanup;
    }
    suffix_begins[batch_local] = suffix_begin;
    dep_rows_by_row[batch_local] = row_dep_rows;
    const UF_long u_len =
      solver->row_refactor_u_ptr[row + 1u] - solver->row_refactor_u_ptr[row];
    if (u_offsets[batch_local] > UF_long_max - u_len) {
      goto cleanup;
    }
    u_offsets[batch_local + 1u] = u_offsets[batch_local] + u_len;
  }
  const UF_long total_u_entries = u_offsets[batch_rows];
  if (total_u_entries > (UF_long)(SIZE_MAX / sizeof(double))) {
    goto cleanup;
  }
  u_workspace = total_u_entries > 0u
    ? (double *)malloc((size_t)total_u_entries * sizeof(*u_workspace))
    : NULL;
  if (total_u_entries > 0u && u_workspace == NULL) {
    goto cleanup;
  }

  double *workspace =
    kls_egraph_worker_supernode_workspace(worker, batch_rows * dep_width);
  if (workspace == NULL) {
    goto cleanup;
  }
  memset(workspace, 0,
         (size_t)batch_rows * (size_t)dep_width * sizeof(*workspace));

  double *x = worker->x;
  double *udiag = (double *)solver->numeric->Udiag;
  for (UF_long local_dep = 0; local_dep < dep_width; ++local_dep) {
    dep_dense_panel[local_dep * dep_width + local_dep] =
      udiag[dep_group_begin + local_dep];
  }

  for (UF_long batch_local = 0; batch_local < batch_rows; ++batch_local) {
    const UF_long row = batch_begin + batch_local;
    const UF_long l_begin = solver->row_refactor_l_ptr[row];
    const UF_long suffix_begin = suffix_begins[batch_local];
    const UF_long row_dep_rows = dep_rows_by_row[batch_local];
    for (UF_long local = 0; local < row_dep_rows; ++local) {
      const UF_long dep_local = suffix_begin + local;
      const UF_long dep = dep_group_begin + dep_local;
      const UF_long p = l_begin + local;
      if (solver->row_refactor_l_cols[p] != dep) {
        kls_egraph_refactor_record_invalid(shared);
        status = -1;
        goto cleanup;
      }
      solver->row_refactor_l_row_values[p] = 0.0;
    }
    const UF_long u_begin = solver->row_refactor_u_ptr[row];
    const UF_long u_end = solver->row_refactor_u_ptr[row + 1u];
    const UF_long u_len = u_end - u_begin;
    double *row_u_workspace = total_u_entries > 0u
      ? u_workspace + u_offsets[batch_local] : NULL;
    double *row_workspace = workspace + batch_local * dep_width;
    double *row_dep_values = row_workspace + suffix_begin;
    const UF_long *row_dep_cols = row_dep_rows > 0u
      ? solver->row_refactor_l_cols + l_begin : NULL;
    const UF_long *row_u_cols = u_len > 0u
      ? solver->row_refactor_u_cols + u_begin : NULL;
    const int direct_input =
      kls_parallel_row_refactor_load_batch_input_row(
        shared, row, row_dep_cols, row_dep_rows, row_dep_values,
        pivots + batch_local, row_u_cols, u_len, row_u_workspace);
    if (direct_input < 0) {
      status = -1;
      goto cleanup;
    }
    if (!direct_input) {
      if (!kls_parallel_row_refactor_load_input_row(shared, x, row)) {
        status = -1;
        goto cleanup;
      }
      for (UF_long local = 0; local < row_dep_rows; ++local) {
        const UF_long dep_local = suffix_begin + local;
        const UF_long dep = dep_group_begin + dep_local;
        row_workspace[dep_local] = x[dep];
        x[dep] = 0.0;
      }
      pivots[batch_local] = x[row];
      x[row] = 0.0;
      for (UF_long p = u_begin; p < u_end; ++p) {
        const UF_long offset = p - u_begin;
        const UF_long col = solver->row_refactor_u_cols[p];
        row_u_workspace[offset] = x[col];
        x[col] = 0.0;
      }
      kls_clear_row_refactor_input_residuals(solver, x, row);
    }
  }

  double trsv_work = 0.0;
  for (UF_long batch_local = 0; batch_local < batch_rows; ++batch_local) {
    const UF_long row = batch_begin + batch_local;
    const UF_long suffix_begin = suffix_begins[batch_local];
    const UF_long l_begin = solver->row_refactor_l_ptr[row];
    double *row_multipliers = workspace + batch_local * dep_width;
    for (UF_long local_dep = suffix_begin; local_dep < dep_width;
         ++local_dep) {
      const UF_long dep = dep_group_begin + local_dep;
      const double candidate = row_multipliers[local_dep];
      const double lij = candidate / udiag[dep];
      if (kls_parallel_row_refactor_rejects_multiplier(
            worker, row, dep, candidate, lij)) {
        status = -1;
        goto cleanup;
      }
      row_multipliers[local_dep] = lij;
      for (UF_long target = local_dep + 1u; target < dep_width; ++target) {
        row_multipliers[target] -=
          lij * dep_dense_panel[local_dep * dep_width + target];
      }
      const UF_long p = l_begin + local_dep - suffix_begin;
      solver->row_refactor_l_row_values[p] = lij;
      if (!shared->row_refactor_defer_value_scatter) {
        *solver->row_refactor_l_values[p] = lij;
      }
      trsv_work += (double)(dep_width - local_dep - 1u);
    }
  }

  double gemv_work = 0.0;
  if (dep_trailing_len > 0u) {
    for (UF_long dep_col_pos = 0; dep_col_pos < dep_trailing_len;
         ++dep_col_pos) {
      const UF_long col = dep_trailing_cols[dep_col_pos];
      for (UF_long batch_local = 0; batch_local < batch_rows; ++batch_local) {
        const UF_long row = batch_begin + batch_local;
        const UF_long suffix_begin = suffix_begins[batch_local];
        const double *row_multipliers = workspace + batch_local * dep_width;
        double update = 0.0;
        for (UF_long local_dep = suffix_begin; local_dep < dep_width;
             ++local_dep) {
          update += row_multipliers[local_dep] *
            dep_trailing_panel[local_dep * dep_trailing_len + dep_col_pos];
        }
        if (update == 0.0) {
          continue;
        }
        if (col == row) {
          pivots[batch_local] -= update;
          continue;
        }
        const UF_long u_begin = solver->row_refactor_u_ptr[row];
        const UF_long u_len =
          solver->row_refactor_u_ptr[row + 1u] - u_begin;
        UF_long u_pos = KLS_KLU_EMPTY;
        if (kls_find_uflong_sorted(solver->row_refactor_u_cols + u_begin,
                                   u_len, col, &u_pos) &&
            u_pos < u_len) {
          u_workspace[u_offsets[batch_local] + u_pos] -= update;
        }
      }
    }
    for (UF_long batch_local = 0; batch_local < batch_rows; ++batch_local) {
      gemv_work += (double)dep_rows_by_row[batch_local] *
        (double)dep_trailing_len;
    }
  }

  for (UF_long batch_local = 0; batch_local < batch_rows; ++batch_local) {
    const UF_long row = batch_begin + batch_local;
    const double pivot = pivots[batch_local];
    const UF_long u_begin = solver->row_refactor_u_ptr[row];
    const UF_long u_end = solver->row_refactor_u_ptr[row + 1u];
    const UF_long u_len = u_end - u_begin;
    double *row_u_workspace = total_u_entries > 0u
      ? u_workspace + u_offsets[batch_local] : NULL;
    const double row_max_abs =
      kls_row_refactor_abs_max_from_values(pivot, row_u_workspace, u_len);
    if (kls_parallel_row_refactor_rejects_pivot(worker, row, pivot,
                                               row_max_abs)) {
      status = -1;
      goto cleanup;
    }
    if (pivot == 0.0) {
      kls_egraph_refactor_record_singular(shared, row,
                                          solver->symbolic->Q[row]);
      if (solver->common.halt_if_singular) {
        status = -1;
        goto cleanup;
      }
    }
    udiag[row] = pivot;
    for (UF_long p = u_begin; p < u_end; ++p) {
      const UF_long offset = p - u_begin;
      solver->row_refactor_u_row_values[p] = row_u_workspace[offset];
      if (!shared->row_refactor_defer_value_scatter) {
        *solver->row_refactor_u_values[p] = row_u_workspace[offset];
      }
    }
    kls_row_refactor_mark_row_done(worker, row);
  }

  const UF_long trsv_entries = kls_work_to_uflong(trsv_work);
  const UF_long gemv_entries = kls_work_to_uflong(gemv_work);
  const UF_long total_entries =
    trsv_entries > UF_long_max - gemv_entries
      ? UF_long_max : trsv_entries + gemv_entries;
  kls_record_row_refactor_compact_supernode_batch(
    solver, batch_rows, candidate_dep_rows, total_entries);
  kls_record_row_refactor_compact_supernode_update(
    solver, candidate_dep_rows, total_entries);
  kls_record_row_refactor_compact_supernode_trsv(
    solver, candidate_dep_rows, trsv_entries);
  if (gemv_entries > 0u) {
    kls_record_row_refactor_compact_supernode_gemv(
      solver, candidate_dep_rows, gemv_entries);
  }
  status = 1;

cleanup:
  free(u_workspace);
  free(suffix_begins);
  free(dep_rows_by_row);
  free(u_offsets);
  free(pivots);
  return status;
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
    if (kind == KLS_ROW_REFACTOR_GROUP_BATCH) {
      return KLS_ROW_REFACTOR_GROUP_BATCH;
    }
  }
  if (solver != NULL && solver->row_refactor_group_dense != NULL &&
      group < solver->row_refactor_group_count &&
      solver->row_refactor_group_dense[group]) {
    return KLS_ROW_REFACTOR_GROUP_DENSE;
  }
  return KLS_ROW_REFACTOR_GROUP_GENERIC;
}

static int kls_row_refactor_group_trailing_slice(
  const kls_solver *solver,
  UF_long group,
  UF_long row_end,
  UF_long *len_out,
  const UF_long **cols_out) {
  KLS_SET_OPTIONAL_OUTPUT(len_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(cols_out, NULL);
  if (solver == NULL || group >= solver->row_refactor_group_count ||
      row_end == 0u || row_end > solver->n ||
      solver->row_refactor_u_ptr == NULL) {
    return 0;
  }

  const UF_long fallback_begin = solver->row_refactor_u_ptr[row_end - 1u];
  const UF_long row_end_ptr = solver->row_refactor_u_ptr[row_end];
  UF_long begin = fallback_begin;
  if (solver->row_refactor_group_trailing_begin != NULL) {
    begin = solver->row_refactor_group_trailing_begin[group];
    if (begin != fallback_begin) {
      return 0;
    }
  }
  if (begin > row_end_ptr) {
    return 0;
  }
  UF_long len = row_end_ptr - begin;
  if (solver->row_refactor_group_trailing_len != NULL) {
    len = solver->row_refactor_group_trailing_len[group];
    if (len != row_end_ptr - begin) {
      return 0;
    }
  }
  if (len > 0u && solver->row_refactor_u_cols == NULL) {
    return 0;
  }

  if (len_out != NULL) {
    *len_out = len;
  }
  if (cols_out != NULL) {
    *cols_out = len > 0u ? solver->row_refactor_u_cols + begin : NULL;
  }
  return 1;
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
  kls_dense_help_row_floor = (long)row_begin;
  const UF_long width = row_end - row_begin;
  const kls_row_refactor_group_kind group_kind =
    kls_row_refactor_group_kind_for(solver, group, width);
  if (group_kind == KLS_ROW_REFACTOR_GROUP_SINGLE) {
    return kls_parallel_row_refactor_process_row(worker, row_begin,
                                                wait_for_dependencies);
  }
  if (group_kind == KLS_ROW_REFACTOR_GROUP_BATCH) {
    UF_long row = row_begin;
    while (row < row_end) {
      UF_long ragged_dep_group = KLS_KLU_EMPTY;
      double ragged_work = 0.0;
      double ragged_copied_entries = 0.0;
      UF_long ragged_dep_rows = 0;
      const UF_long ragged_end =
        kls_batch_group_ragged_supernode_end(
          solver, group, row_begin, row_end, row, &ragged_dep_group,
          &ragged_work, &ragged_copied_entries, &ragged_dep_rows);
      if (ragged_end > row + 1u) {
        const int ragged_status =
          kls_batch_group_try_ragged_supernode_update(
            worker, group, row, ragged_end, ragged_dep_group,
            ragged_work, ragged_copied_entries, ragged_dep_rows);
        if (ragged_status < 0) {
          return 0;
        }
        if (ragged_status > 0) {
          row = ragged_end;
          continue;
        }
      }
      if (!kls_parallel_row_refactor_process_row(worker, row,
                                                wait_for_dependencies)) {
        return 0;
      }
      kls_egraph_refactor_mark_done(shared, row);
      row++;
    }
    return 1;
  }

  trilinos_klu_l_numeric *numeric = solver->numeric;
  const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  double *x = worker->x;
  double *udiag = (double *)numeric->Udiag;
  UF_long trailing_len = 0;
  const UF_long *trailing_cols = NULL;
  if (!kls_row_refactor_group_trailing_slice(
        solver, group, row_end, &trailing_len, &trailing_cols)) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }

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
    const int direct_input =
      kls_parallel_row_refactor_load_sparse_segment_input_row(
        shared, x, row_begin, row_end, trailing_len, trailing_cols, row);
    if (direct_input < 0) {
      return 0;
    }
    if (!direct_input &&
        !kls_parallel_row_refactor_load_input_row(shared, x, row)) {
      return 0;
    }

    const UF_long l_begin = solver->row_refactor_l_ptr[row];
    UF_long p = l_begin;
    const UF_long l_end = solver->row_refactor_l_ptr[row + 1u];
    const UF_long dep_count = l_end - l_begin;
    unsigned char *applied = NULL;
    if (wait_for_dependencies && dep_count > 1u &&
        (uintmax_t)dep_count <= (uintmax_t)SIZE_MAX / sizeof(*applied)) {
      applied = kls_egraph_worker_byte_workspace(worker, dep_count);
      if (applied != NULL) {
        memset(applied, 0, (size_t)dep_count);
      }
    }
    UF_long prefactor_deps = 0;
    UF_long prefactor_supernode_deps = 0;
    while (p < l_end) {
      if (applied != NULL && applied[p - l_begin]) {
        ++p;
        continue;
      }
      const UF_long compact_begin = p;
      const int compact_status =
        kls_row_refactor_try_compact_supernode_update(
          worker, row, group, &p, l_begin, l_end,
          wait_for_dependencies, applied);
      if (compact_status < 0) {
        x[row] = 0.0;
        return 0;
      }
      if (compact_status > 0) {
        if (applied != NULL && p > compact_begin && p <= l_end) {
          memset(applied + (compact_begin - l_begin), 1,
                 (size_t)(p - compact_begin));
        }
        continue;
      }
      const UF_long dep = solver->row_refactor_l_cols[p];
      if (wait_for_dependencies && dep < row_begin &&
          !kls_row_refactor_dependency_done_now(shared, dep)) {
        if (applied != NULL &&
            !kls_parallel_row_refactor_prefactor_finished_deps(
              worker, row, group, row_begin, l_begin, l_end, p, applied,
              x, &prefactor_deps, &prefactor_supernode_deps)) {
          x[dep] = 0.0;
          return 0;
        }
        if (!kls_egraph_refactor_wait_done(shared, dep)) {
          x[dep] = 0.0;
          return 0;
        }
      }
      if (dep < row_begin) {
        const int scalar_status =
          kls_parallel_row_refactor_apply_scalar_dep(
            worker, row, group, p, x, 1);
        if (scalar_status != 1) {
          x[dep] = 0.0;
          return 0;
        }
      } else {
        const double candidate =
          direct_input ? solver->row_refactor_l_row_values[p] + x[dep]
                       : x[dep];
        const double lij = candidate / udiag[dep];
        if (kls_parallel_row_refactor_rejects_multiplier(
              worker, row, dep, candidate, lij)) {
          x[dep] = 0.0;
          return 0;
        }
        solver->row_refactor_l_row_values[p] = lij;
        if (!shared->row_refactor_defer_value_scatter) {
          *solver->row_refactor_l_values[p] = lij;
        }
        x[dep] = 0.0;

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
      }
      if (applied != NULL) {
        applied[p - l_begin] = 1u;
      }
      p++;
    }

    const double pivot = direct_input ? udiag[row] + x[row] : x[row];
    const UF_long dense_len = row_end - row - 1u;
    const UF_long u_begin = solver->row_refactor_u_ptr[row];
    const UF_long u_end = solver->row_refactor_u_ptr[row + 1u];
    if (u_end - u_begin != dense_len + trailing_len) {
      kls_egraph_refactor_record_invalid(shared);
      return 0;
    }
    double *row_u_values = solver->row_refactor_u_row_values + u_begin;
    if (shared->check_pivots) {
      double row_max_abs = -1.0;
      if (!kls_row_refactor_segment_u_abs_max(
            row, row_end, trailing_len, trailing_cols, x, row_u_values,
            direct_input, pivot, &row_max_abs)) {
        kls_egraph_refactor_record_invalid(shared);
        return 0;
      }
      if (kls_parallel_row_refactor_rejects_pivot(worker, row, pivot,
                                                 row_max_abs)) {
        memset(x, 0, (size_t)solver->n * sizeof(*x));
        return 0;
      }
    }
    x[row] = 0.0;
    if (pivot == 0.0) {
      kls_egraph_refactor_record_singular(shared, row, symbolic->Q[row]);
      if (solver->common.halt_if_singular) {
        memset(x, 0, (size_t)solver->n * sizeof(*x));
        return 0;
      }
    }
    udiag[row] = pivot;

    if (shared->row_refactor_defer_value_scatter) {
      for (UF_long offset = 0; offset < dense_len; ++offset) {
        const UF_long col = row + 1u + offset;
        row_u_values[offset] =
          direct_input ? row_u_values[offset] + x[col] : x[col];
        x[col] = 0.0;
      }
      for (UF_long offset = 0; offset < trailing_len; ++offset) {
        const UF_long col = trailing_cols[offset];
        row_u_values[dense_len + offset] =
          direct_input ? row_u_values[dense_len + offset] + x[col] : x[col];
        x[col] = 0.0;
      }
    } else {
      double **klu_u_values = solver->row_refactor_u_values + u_begin;
      for (UF_long offset = 0; offset < dense_len; ++offset) {
        const UF_long col = row + 1u + offset;
        const double value =
          direct_input ? row_u_values[offset] + x[col] : x[col];
        row_u_values[offset] = value;
        *klu_u_values[offset] = value;
        x[col] = 0.0;
      }
      for (UF_long offset = 0; offset < trailing_len; ++offset) {
        const UF_long col = trailing_cols[offset];
        const double value =
          direct_input ? row_u_values[dense_len + offset] + x[col] : x[col];
        row_u_values[dense_len + offset] = value;
        *klu_u_values[dense_len + offset] = value;
        x[col] = 0.0;
      }
    }
    kls_clear_row_refactor_input_residuals(solver, x, row);
    if (prefactor_deps > 0u) {
      atomic_fetch_add_explicit(&shared->row_prefactor_rows, 1ul,
                                memory_order_relaxed);
      atomic_fetch_add_explicit(&shared->row_prefactor_deps,
                                (unsigned long)prefactor_deps,
                                memory_order_relaxed);
    }
    if (prefactor_supernode_deps > 0u) {
      atomic_fetch_add_explicit(&shared->row_prefactor_supernode_rows, 1ul,
                                memory_order_relaxed);
      atomic_fetch_add_explicit(&shared->row_prefactor_supernode_deps,
                                (unsigned long)prefactor_supernode_deps,
                                memory_order_relaxed);
    }
    kls_row_refactor_mark_row_done(worker, row);
  }
  return 1;
}

static void kls_row_refactor_refresh_missing_prefix(
  kls_egraph_refactor_shared *shared,
  double *x) {
  if (shared == NULL || shared->solver == NULL ||
      shared->pipeline_done == NULL || x == NULL ||
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

  UF_long rejected_pivot = shared->rejected_pivot;
  for (UF_long row = 0; row < rejected_pivot; ++row) {
    if (atomic_load_explicit(&shared->pipeline_done[row],
                             memory_order_acquire) ==
        shared->pipeline_generation) {
      continue;
    }
    memset(x, 0, (size_t)solver->n * sizeof(*x));
    if (!kls_parallel_row_refactor_process_row(&worker, row, 0)) {
      break;
    }
    kls_egraph_refactor_mark_done(shared, row);
    if (shared->rejected_pivot != KLS_KLU_EMPTY &&
        shared->rejected_pivot < rejected_pivot) {
      rejected_pivot = shared->rejected_pivot;
    }
  }
  memset(x, 0, (size_t)solver->n * sizeof(*x));
  free(worker.segment_panel);
  free(worker.supernode_workspace);
  free(worker.index_workspace);
  free(worker.byte_workspace);
  free(worker.object_workspace);
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

static int kls_row_refactor_enqueue_ready_group(
  kls_egraph_refactor_shared *shared,
  UF_long group) {
  if (shared == NULL || !shared->row_pipeline_ready_queue ||
      shared->solver == NULL ||
      shared->row_pipeline_ready_groups == NULL ||
      shared->row_pipeline_ready_slots == NULL ||
      shared->row_pipeline_tail_groups == NULL ||
      shared->row_pipeline_tail_count == 0u ||
      group >= shared->solver->row_refactor_group_count ||
      !shared->row_pipeline_tail_groups[group]) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }
  const unsigned long pos =
    atomic_fetch_add_explicit(&shared->row_pipeline_ready_tail, 1ul,
                              memory_order_acq_rel);
  if ((UF_long)pos >= shared->row_pipeline_tail_count) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }
  shared->row_pipeline_ready_groups[pos] = group;
  atomic_store_explicit(&shared->row_pipeline_ready_slots[pos], 1u,
                        memory_order_release);
  return 1;
}

static int kls_row_refactor_pop_ready_group(
  kls_egraph_refactor_worker *worker,
  UF_long *group_out) {
  if (group_out != NULL) {
    *group_out = KLS_KLU_EMPTY;
  }
  if (worker == NULL || worker->shared == NULL || group_out == NULL) {
    return 0;
  }
  kls_egraph_refactor_shared *shared = worker->shared;
  if (!shared->row_pipeline_ready_queue ||
      shared->row_pipeline_ready_groups == NULL ||
      shared->row_pipeline_ready_slots == NULL ||
      shared->row_pipeline_tail_count == 0u) {
    return 0;
  }

  unsigned spin = 0;
  for (;;) {
    if (kls_egraph_refactor_should_stop(shared)) {
      return 0;
    }
    if ((UF_long)atomic_load_explicit(
          &shared->row_pipeline_completed_groups, memory_order_acquire) >=
        shared->row_pipeline_tail_count) {
      return 0;
    }

    const unsigned long head =
      atomic_load_explicit(&shared->row_pipeline_ready_head,
                           memory_order_acquire);
    const unsigned long tail =
      atomic_load_explicit(&shared->row_pipeline_ready_tail,
                           memory_order_acquire);
    if (head >= tail) {
      if ((spin & 63u) == 1u && kls_dense_help_try(shared, LONG_MAX, 1)) {
        continue;
      }
      kls_cpu_relax();
      if ((spin++ & 1023u) == 0u &&
          kls_egraph_refactor_should_stop(shared)) {
        return 0;
      }
      continue;
    }

    unsigned long expected = head;
    if (!atomic_compare_exchange_weak_explicit(
          &shared->row_pipeline_ready_head, &expected, head + 1ul,
          memory_order_acq_rel, memory_order_acquire)) {
      continue;
    }
    if ((UF_long)head >= shared->row_pipeline_tail_count) {
      kls_egraph_refactor_record_invalid(shared);
      return 0;
    }
    while (atomic_load_explicit(&shared->row_pipeline_ready_slots[head],
                                memory_order_acquire) == 0u) {
                                  kls_cpu_relax();
      if (kls_egraph_refactor_should_stop(shared)) {
        return 0;
      }
    }
    *group_out = shared->row_pipeline_ready_groups[head];
    return 1;
  }
}

static int kls_row_refactor_finish_ready_group(
  kls_egraph_refactor_worker *worker,
  UF_long group,
  UF_long *local_ready_out) {
  if (local_ready_out != NULL) {
    *local_ready_out = KLS_KLU_EMPTY;
  }
  if (worker == NULL || worker->shared == NULL) {
    return 0;
  }
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  if (shared->row_pipeline_row_dep_ready_queue) {
    if (solver == NULL || group >= solver->row_refactor_group_count ||
        shared->row_pipeline_tail_groups == NULL ||
        !shared->row_pipeline_tail_groups[group]) {
      kls_egraph_refactor_record_invalid(shared);
      return 0;
    }
    if (shared->row_pipeline_private_group_mask != NULL) {
      if (solver->row_refactor_group_successor_ptr == NULL ||
          (solver->row_refactor_group_dependency_edges > 0u &&
           solver->row_refactor_group_successor_groups == NULL) ||
          shared->row_pipeline_remaining_preds == NULL) {
        kls_egraph_refactor_record_invalid(shared);
        return 0;
      }
      const UF_long begin = solver->row_refactor_group_successor_ptr[group];
      const UF_long end =
        solver->row_refactor_group_successor_ptr[group + 1u];
      if (begin > end ||
          end > solver->row_refactor_group_dependency_edges) {
        kls_egraph_refactor_record_invalid(shared);
        return 0;
      }
      for (UF_long pos = begin; pos < end; ++pos) {
        const UF_long successor =
          solver->row_refactor_group_successor_groups[pos];
        if (successor >= solver->row_refactor_group_count) {
          kls_egraph_refactor_record_invalid(shared);
          return 0;
        }
        if (!shared->row_pipeline_tail_groups[successor] ||
            !shared->row_pipeline_private_group_mask[successor]) {
          continue;
        }
        const unsigned long old =
          atomic_fetch_sub_explicit(
            &shared->row_pipeline_remaining_preds[successor], 1ul,
            memory_order_acq_rel);
        if (old == 0ul) {
          kls_egraph_refactor_record_invalid(shared);
          return 0;
        }
      }
    }
    atomic_fetch_add_explicit(&shared->row_pipeline_completed_groups, 1ul,
                              memory_order_release);
    return 1;
  }
  if (solver == NULL ||
      solver->row_refactor_group_successor_ptr == NULL ||
      (solver->row_refactor_group_dependency_edges > 0u &&
       solver->row_refactor_group_successor_groups == NULL) ||
      group >= solver->row_refactor_group_count ||
      shared->row_pipeline_remaining_preds == NULL ||
      shared->row_pipeline_tail_groups == NULL) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }

  const UF_long begin = solver->row_refactor_group_successor_ptr[group];
  const UF_long end = solver->row_refactor_group_successor_ptr[group + 1u];
  if (begin > end ||
      end > solver->row_refactor_group_dependency_edges) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }
  for (UF_long pos = begin; pos < end; ++pos) {
    const UF_long successor = solver->row_refactor_group_successor_groups[pos];
    if (successor >= solver->row_refactor_group_count) {
      kls_egraph_refactor_record_invalid(shared);
      return 0;
    }
    if (!shared->row_pipeline_tail_groups[successor]) {
      continue;
    }
    const unsigned long old =
      atomic_fetch_sub_explicit(
        &shared->row_pipeline_remaining_preds[successor], 1ul,
        memory_order_acq_rel);
    if (old == 0ul) {
      kls_egraph_refactor_record_invalid(shared);
      return 0;
    }
    if (old == 1ul) {
      if (shared->row_pipeline_private_group_mask != NULL &&
          shared->row_pipeline_private_group_mask[successor]) {
        continue;
      }
      if (local_ready_out != NULL &&
          *local_ready_out == KLS_KLU_EMPTY) {
        *local_ready_out = successor;
        atomic_fetch_add_explicit(&shared->row_pipeline_local_ready_groups,
                                  1ul, memory_order_relaxed);
      } else if (!kls_row_refactor_enqueue_ready_group(shared, successor)) {
        return 0;
      }
    }
  }

  atomic_fetch_add_explicit(&shared->row_pipeline_completed_groups, 1ul,
                            memory_order_release);
  return 1;
}

static int kls_row_refactor_finish_ordered_private_group(
  kls_egraph_refactor_worker *worker,
  UF_long group) {
  if (worker == NULL || worker->shared == NULL) {
    return 0;
  }
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  if (solver == NULL ||
      solver->row_refactor_group_successor_ptr == NULL ||
      (solver->row_refactor_group_dependency_edges > 0u &&
       solver->row_refactor_group_successor_groups == NULL) ||
      group >= solver->row_refactor_group_count ||
      shared->row_pipeline_remaining_preds == NULL ||
      shared->row_pipeline_tail_groups == NULL ||
      shared->row_pipeline_private_group_mask == NULL ||
      !shared->row_pipeline_private_group_mask[group] ||
      shared->row_pipeline_row_dep_ready_queue) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }

  const UF_long begin = solver->row_refactor_group_successor_ptr[group];
  const UF_long end = solver->row_refactor_group_successor_ptr[group + 1u];
  if (begin > end ||
      end > solver->row_refactor_group_dependency_edges) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }
  for (UF_long pos = begin; pos < end; ++pos) {
    const UF_long successor = solver->row_refactor_group_successor_groups[pos];
    if (successor >= solver->row_refactor_group_count) {
      kls_egraph_refactor_record_invalid(shared);
      return 0;
    }
    if (!shared->row_pipeline_tail_groups[successor] ||
        shared->row_pipeline_private_group_mask[successor]) {
      continue;
    }
    const unsigned long old =
      atomic_fetch_sub_explicit(
        &shared->row_pipeline_remaining_preds[successor], 1ul,
        memory_order_acq_rel);
    if (old == 0ul) {
      kls_egraph_refactor_record_invalid(shared);
      return 0;
    }
    if (old == 1ul &&
        !kls_row_refactor_enqueue_ready_group(shared, successor)) {
      return 0;
    }
  }

  atomic_fetch_add_explicit(&shared->row_pipeline_completed_groups, 1ul,
                            memory_order_release);
  return 1;
}

static int kls_row_refactor_wait_private_group_ready(
  kls_egraph_refactor_worker *worker,
  UF_long group) {
  if (worker == NULL || worker->shared == NULL) {
    return 0;
  }
  kls_egraph_refactor_shared *shared = worker->shared;
  if (shared->row_pipeline_private_group_mask == NULL) {
    return 1;
  }
  if (shared->solver == NULL ||
      group >= shared->solver->row_refactor_group_count ||
      shared->row_pipeline_remaining_preds == NULL ||
      !shared->row_pipeline_private_group_mask[group]) {
    kls_egraph_refactor_record_invalid(shared);
    return 0;
  }
  const long help_bound =
    shared->solver->row_refactor_group_ptr != NULL
      ? (long)shared->solver->row_refactor_group_ptr[group] : -1;
  unsigned spin = 0;
  while (atomic_load_explicit(&shared->row_pipeline_remaining_preds[group],
                              memory_order_acquire) != 0ul) {
    if ((spin & 63u) == 1u && kls_dense_help_try(shared, help_bound, 0)) {
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

static int kls_row_refactor_process_private_group(
  kls_egraph_refactor_worker *worker,
  UF_long group) {
  if (worker == NULL || worker->shared == NULL) {
    return 0;
  }
  if (!kls_row_refactor_wait_private_group_ready(worker, group)) {
    return 0;
  }
  if (!kls_parallel_row_refactor_process_group(worker, group, 0)) {
    return 0;
  }
  kls_parallel_row_refactor_mark_group_done(worker->shared,
                                           worker->shared->solver,
                                           group);
  return kls_row_refactor_finish_ready_group(worker, group, NULL);
}

static int kls_row_refactor_process_ordered_private_group(
  kls_egraph_refactor_worker *worker,
  UF_long group) {
  if (worker == NULL || worker->shared == NULL) {
    return 0;
  }
  if (!kls_parallel_row_refactor_process_group(worker, group, 0)) {
    return 0;
  }
  kls_parallel_row_refactor_mark_group_done(worker->shared,
                                           worker->shared->solver,
                                           group);
  return kls_row_refactor_finish_ordered_private_group(worker, group);
}

static int kls_row_refactor_process_ready_group_with_local_tail(
  kls_egraph_refactor_worker *worker,
  UF_long group) {
  if (worker == NULL || worker->shared == NULL) {
    return 0;
  }
  for (;;) {
    if (kls_egraph_refactor_should_stop(worker->shared)) {
      return 0;
    }
    const int wait_for_dependencies =
      worker->shared->row_pipeline_row_dep_ready_queue ? 1 : 0;
    if (!kls_parallel_row_refactor_process_group(
          worker, group, wait_for_dependencies)) {
      return 0;
    }
    kls_parallel_row_refactor_mark_group_done(worker->shared,
                                             worker->shared->solver,
                                             group);
    UF_long local_group = KLS_KLU_EMPTY;
    if (!kls_row_refactor_finish_ready_group(worker, group, &local_group)) {
      return 0;
    }
    if (local_group == KLS_KLU_EMPTY) {
      return 1;
    }
    group = local_group;
  }
}

static void kls_row_refactor_worker_run_ready_pipeline(
  kls_egraph_refactor_worker *worker) {
  if (worker != NULL && worker->shared != NULL &&
      worker->shared->row_pipeline_private_groups != NULL) {
    kls_egraph_refactor_shared *shared = worker->shared;
    if (shared->row_pipeline_private_thread_ptr != NULL &&
        worker->tid >= 0 && worker->tid < shared->thread_count) {
      const UF_long begin =
        shared->row_pipeline_private_thread_ptr[worker->tid];
      const UF_long end =
        shared->row_pipeline_private_thread_ptr[worker->tid + 1];
      if (begin <= end && end <= shared->row_pipeline_private_count) {
        for (UF_long pos = begin; pos < end; ++pos) {
          if (kls_egraph_refactor_should_stop(shared)) {
            return;
          }
          const UF_long group = shared->row_pipeline_private_groups[pos];
          int ok = 0;
          if (shared->row_pipeline_ordered_private) {
            ok = kls_row_refactor_process_ordered_private_group(worker,
                                                               group);
          } else if (shared->row_pipeline_private_group_mask != NULL) {
            ok = kls_row_refactor_process_private_group(worker, group);
          } else {
            ok = kls_row_refactor_process_ready_group_with_local_tail(
              worker, group);
          }
          if (!ok) {
            return;
          }
        }
      } else {
        kls_egraph_refactor_record_invalid(shared);
        return;
      }
    } else {
      for (;;) {
        if (kls_egraph_refactor_should_stop(shared)) {
          return;
        }
        const UF_long pos =
          (UF_long)atomic_fetch_add_explicit(
            &shared->row_pipeline_private_pos, 1ul, memory_order_relaxed);
        if (pos >= shared->row_pipeline_private_count) {
          break;
        }
        const UF_long group = shared->row_pipeline_private_groups[pos];
        int ok = 0;
        if (shared->row_pipeline_ordered_private) {
          ok = kls_row_refactor_process_ordered_private_group(worker,
                                                             group);
        } else if (shared->row_pipeline_private_group_mask != NULL) {
          ok = kls_row_refactor_process_private_group(worker, group);
        } else {
          ok = kls_row_refactor_process_ready_group_with_local_tail(worker,
                                                                   group);
        }
        if (!ok) {
          return;
        }
      }
    }
  }

  for (;;) {
    UF_long group = KLS_KLU_EMPTY;
    if (!kls_row_refactor_pop_ready_group(worker, &group)) {
      break;
    }
    if (!kls_row_refactor_process_ready_group_with_local_tail(worker,
                                                             group)) {
      break;
    }
  }
}

static void kls_row_refactor_worker_run(kls_egraph_refactor_worker *worker) {
  if (worker == NULL || worker->shared == NULL) {
    return;
  }
  kls_dense_help_self = worker;
  kls_dense_help_row_floor = LONG_MAX;
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
    }
    if ((shared->row_pipeline_ready_queue ||
         shared->pipeline_done != NULL) &&
        cluster_levels < solver->row_refactor_level_count &&
        !kls_egraph_refactor_should_stop(shared)) {
      if (shared->row_pipeline_ready_queue) {
        kls_row_refactor_worker_run_ready_pipeline(worker);
      } else {
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
            if (!kls_egraph_refactor_should_stop(shared)) {
              kls_egraph_refactor_record_invalid(shared);
            }
            break;
          }
          kls_parallel_row_refactor_mark_group_done(shared, solver, group);
        }
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
  }
}

/* Publish deferred row-major refactor values into the packed column
   factors when the follow-up solves are better served by the column
   engine (pts/i32) than by the row-value solve.  Multi-block factors can
   use the column route when its machinery is available. */
static int kls_row_refactor_should_publish_for_solve(
  const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL ||
      solver->options.threads <= 1) {
    return 0;
  }

  if (solver->row_accept_publish_preferred) {
    /* Trial solves measured the column route decisively faster; that
       measurement is the evidence, whichever column path produced it. */
    return 1;
  }
  if (solver->row_publish_experiment == 1 ||
      solver->row_publish_experiment == 2) {
    /* the steady audit needs one column-route solve to compare against
       an adoption that never measured its solves; keep publishing until
       that solve lands and concludes the probe */
    return 1;
  }
  if (solver->symbolic->nblocks == 1u && solver->row_solve_partition_ready) {
    return 0;
  }
  return kls_pts_solve_available(solver);
}

static void kls_row_publish_worker_run(kls_egraph_refactor_worker *worker) {
  if (worker == NULL || worker->shared == NULL) {
    return;
  }
  kls_egraph_refactor_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  if (solver == NULL || solver->row_refactor_l_ptr == NULL ||
      solver->row_refactor_u_ptr == NULL) {
    return;
  }
  const UF_long tid = (UF_long)worker->tid;
  const UF_long threads = (UF_long)shared->thread_count;
  if (threads == 0u || tid >= threads) {
    return;
  }
  const UF_long l_nnz = solver->row_refactor_l_ptr[solver->n];
  if (l_nnz > 0u && solver->row_refactor_l_values != NULL &&
      solver->row_refactor_l_row_values != NULL) {
    const UF_long begin = l_nnz * tid / threads;
    const UF_long end = l_nnz * (tid + 1u) / threads;
    double *const *values = solver->row_refactor_l_values;
    const double *row_values = solver->row_refactor_l_row_values;
    for (UF_long p = begin; p < end; ++p) {
      *values[p] = row_values[p];
    }
  }
  const UF_long u_nnz = solver->row_refactor_u_ptr[solver->n];
  if (u_nnz > 0u && solver->row_refactor_u_values != NULL &&
      solver->row_refactor_u_row_values != NULL) {
    const UF_long begin = u_nnz * tid / threads;
    const UF_long end = u_nnz * (tid + 1u) / threads;
    double *const *values = solver->row_refactor_u_values;
    const double *row_values = solver->row_refactor_u_row_values;
    for (UF_long p = begin; p < end; ++p) {
      *values[p] = row_values[p];
    }
  }
}

static int kls_threaded_row_refactor_numeric(kls_solver *solver,
                                             double *numeric_values,
                                             int check_pivots) {
  if (solver != NULL && solver->numeric_is_predicted) {
    return -1;
  }
  const int scaled = solver != NULL && solver->common.scale > 0;
  if (numeric_values == NULL || solver == NULL || solver->symbolic == NULL ||
      solver->numeric == NULL || solver->symbolic->nblocks == 0u ||
      (scaled && (solver->numeric->Rs == NULL ||
                  solver->numeric->Pnum == NULL ||
                  solver->numeric->Xwork == NULL)) ||
      (!scaled && solver->numeric->Rs != NULL) ||
      solver->numeric->Udiag == NULL || solver->options.threads <= 1 ||
      !kls_build_row_refactor_pattern(solver, 0) ||
      solver->row_refactor_level_max_width < 2u) {
    return -1;
  }
  kls_row_refactor_ensure_sn_end(solver);

  int thread_count = solver->options.threads;
  if ((UF_long)thread_count > solver->row_refactor_level_max_width) {
    thread_count = (int)solver->row_refactor_level_max_width;
  }
  if (thread_count < 2) {
    return -1;
  }
  int use_group_pipeline = 0;
  UF_long pipeline_group_begin = 0;
  UF_long cluster_levels = 0;
  const int have_group_dag =
    solver->row_refactor_group_level_ptr != NULL &&
    solver->row_refactor_level_groups != NULL &&
    solver->row_refactor_group_ptr != NULL &&
    solver->row_refactor_group_successor_ptr != NULL &&
    solver->row_refactor_group_count > 0u;
  if (have_group_dag) {
    cluster_levels =
      kls_row_refactor_choose_cluster_levels(solver, thread_count);
    if (cluster_levels > solver->row_refactor_level_count) {
      cluster_levels = solver->row_refactor_level_count;
    }
    use_group_pipeline = cluster_levels < solver->row_refactor_level_count;
  } else {
    cluster_levels = solver->row_refactor_level_count;
  }

  if (have_group_dag && cluster_levels > 0u &&
      !kls_build_row_refactor_group_thread_slices(solver, thread_count)) {
    return -1;
  }

  UF_long partial_supernode_groups = 0;
  UF_long partial_supernode_rows = 0;
  const int prefer_partial_supernode_pipeline =
    have_group_dag && use_group_pipeline &&
    kls_row_refactor_tail_has_large_partial_supernodes(
      solver, cluster_levels, thread_count, &partial_supernode_groups,
      &partial_supernode_rows);

  UF_long *row_ready_groups = NULL;
  atomic_uint *row_ready_slots = NULL;
  atomic_ulong *row_remaining_preds = NULL;
  unsigned char *row_tail_groups = NULL;
  const UF_long *row_private_groups = NULL;
  UF_long *row_owned_private_groups = NULL;
  UF_long *row_private_thread_ptr = NULL;
  unsigned char *row_private_group_mask = NULL;
  UF_long row_private_group_count = 0;
  UF_long row_separator_private_component_count = 0;
  int row_separator_private_queue = 0;
  int row_separator_flop_queue = 0;
  int row_ready_initial_queued = 0;
  int row_dep_ready_queue = 0;
  UF_long row_ready_tail_count = 0;
  UF_long row_initial_ready_count = 0;
  UF_long row_separator_flop_private_groups = 0;
  UF_long row_separator_flop_pipeline_groups = 0;
  UF_long row_separator_flop_closure_groups = 0;
  UF_long row_separator_flop_component_count = 0;
  UF_long row_separator_flop_private_threads = 0;
  UF_long row_separator_flop_private_min_groups = 0;
  UF_long row_separator_flop_private_max_groups = 0;
  double row_separator_flop_private_min_work = 0.0;
  double row_separator_flop_private_max_work = 0.0;
  int row_ordered_private_queue = 0;
  /* Preserve the retained separator-tree private/pipeline queue when it covers
     the row-group DAG.  When the separator schedule is pipeline-heavy and a
     work-significant dependent producer is present, wrap its pipeline side in
     SubtreeLU Algorithm 5's row-prefix release.  The work-share selector keeps
     the paper mechanism but avoids building a row-edge graph for producers too
     small to occupy one worker for their fair share of the tail.  This is
     prefix-safe for checked runs because each producer row is marked done only
     after its pivot check and row-value publication complete.  If preparation
     fails, keep the older cluster/tail schedule as a conservative fallback. */
  int use_row_ready_queue = 0;
  int prepared_separator_flop_queue = 0;
  int prepared_separator_row_dep_queue = 0;
  if (have_group_dag &&
      kls_prepare_row_refactor_separator_flop_ready_queue(
        solver, thread_count, &row_ready_groups, &row_ready_slots,
        &row_remaining_preds, &row_tail_groups, &row_owned_private_groups,
        &row_private_thread_ptr, &row_private_group_mask,
        &row_ready_tail_count, &row_initial_ready_count,
        &row_separator_flop_private_groups,
        &row_separator_flop_pipeline_groups,
        &row_separator_flop_closure_groups,
        &row_separator_flop_component_count,
        &row_separator_flop_private_threads,
        &row_separator_flop_private_min_groups,
        &row_separator_flop_private_max_groups,
        &row_separator_flop_private_min_work,
        &row_separator_flop_private_max_work)) {
    prepared_separator_flop_queue = 1;
  }
  if (prepared_separator_flop_queue &&
      prefer_partial_supernode_pipeline &&
      row_separator_flop_pipeline_groups >= row_separator_flop_private_groups) {
    prepared_separator_row_dep_queue =
      kls_prepare_row_refactor_row_dep_ready_queue_with_private(
          solver, row_private_group_mask, &row_ready_groups, &row_ready_slots,
          &row_remaining_preds, &row_tail_groups, &row_ready_tail_count,
          &row_initial_ready_count);
    if (!prepared_separator_row_dep_queue) {
      free(row_owned_private_groups);
      row_owned_private_groups = NULL;
      free(row_private_thread_ptr);
      row_private_thread_ptr = NULL;
      free(row_private_group_mask);
      row_private_group_mask = NULL;
      row_ready_groups = NULL;
      row_ready_slots = NULL;
      row_remaining_preds = NULL;
      row_tail_groups = NULL;
      row_ready_tail_count = 0;
      row_initial_ready_count = 0;
      row_separator_flop_private_groups = 0;
      row_separator_flop_pipeline_groups = 0;
      row_separator_flop_closure_groups = 0;
      row_separator_flop_component_count = 0;
      row_separator_flop_private_threads = 0;
      row_separator_flop_private_min_groups = 0;
      row_separator_flop_private_max_groups = 0;
      row_separator_flop_private_min_work = 0.0;
      row_separator_flop_private_max_work = 0.0;
      prepared_separator_flop_queue =
        kls_prepare_row_refactor_separator_flop_ready_queue(
          solver, thread_count, &row_ready_groups, &row_ready_slots,
          &row_remaining_preds, &row_tail_groups, &row_owned_private_groups,
          &row_private_thread_ptr, &row_private_group_mask,
          &row_ready_tail_count, &row_initial_ready_count,
          &row_separator_flop_private_groups,
          &row_separator_flop_pipeline_groups,
          &row_separator_flop_closure_groups,
          &row_separator_flop_component_count,
          &row_separator_flop_private_threads,
          &row_separator_flop_private_min_groups,
          &row_separator_flop_private_max_groups,
          &row_separator_flop_private_min_work,
          &row_separator_flop_private_max_work);
    }
  }
  if (prepared_separator_row_dep_queue) {
    use_row_ready_queue = 1;
    row_separator_flop_queue = 1;
    row_dep_ready_queue = 1;
    row_ready_initial_queued = 1;
    row_private_groups = row_owned_private_groups;
    row_private_group_count = row_separator_flop_private_groups;
    cluster_levels = 0;
    use_group_pipeline = solver->row_refactor_level_count > 0u;
  } else if (prepared_separator_flop_queue) {
    use_row_ready_queue = 1;
    row_separator_flop_queue = 1;
    row_ready_initial_queued = 1;
    row_private_groups = row_owned_private_groups;
    row_private_group_count = row_separator_flop_private_groups;
    cluster_levels = 0;
    use_group_pipeline = solver->row_refactor_level_count > 0u;
  }
  if (!use_row_ready_queue) {
    use_row_ready_queue =
      have_group_dag && use_group_pipeline &&
      (prefer_partial_supernode_pipeline
         ? (row_dep_ready_queue =
              kls_prepare_row_refactor_row_dep_ready_queue(
                solver, cluster_levels, &row_ready_groups, &row_ready_slots,
                &row_remaining_preds, &row_tail_groups, &row_private_groups,
                &row_ready_tail_count, &row_initial_ready_count))
         : kls_prepare_row_refactor_ready_queue(
             solver, cluster_levels, &row_ready_groups, &row_ready_slots,
             &row_remaining_preds, &row_tail_groups, &row_private_groups,
             &row_ready_tail_count, &row_initial_ready_count));
    row_ready_initial_queued = use_row_ready_queue &&
      row_private_groups == NULL;
  }

  if (have_group_dag && !use_row_ready_queue &&
      !prefer_partial_supernode_pipeline) {
    if (!kls_build_row_refactor_group_thread_slices(solver, thread_count)) {
      return -1;
    }
    cluster_levels =
      kls_row_refactor_choose_cluster_levels(solver, thread_count);
    use_group_pipeline = 0;
    pipeline_group_begin = 0;
    row_ready_groups = NULL;
    row_ready_slots = NULL;
    row_remaining_preds = NULL;
    row_tail_groups = NULL;
    row_private_groups = NULL;
    row_private_group_mask = NULL;
    row_private_group_count = 0;
    row_ready_tail_count = 0;
    row_initial_ready_count = 0;
    row_ready_initial_queued = 0;
  }

  if (solver->row_refactor_group_level_ptr != NULL &&
      solver->row_refactor_level_groups != NULL &&
      solver->row_refactor_group_ptr != NULL &&
      have_group_dag &&
      !use_row_ready_queue &&
      cluster_levels < solver->row_refactor_level_count) {
    pipeline_group_begin = solver->row_refactor_group_level_ptr[cluster_levels];
    use_group_pipeline = pipeline_group_begin < solver->row_refactor_group_count;
  }
  if (!use_group_pipeline) {
    cluster_levels = solver->row_refactor_level_count;
  }
  kls_row_refactor_record_pipeline_scope(solver, cluster_levels);
  if (use_group_pipeline && !use_row_ready_queue) {
    if (prefer_partial_supernode_pipeline) {
      use_row_ready_queue =
        kls_prepare_row_refactor_row_dep_ready_queue(
          solver, cluster_levels, &row_ready_groups, &row_ready_slots,
          &row_remaining_preds, &row_tail_groups, &row_private_groups,
          &row_ready_tail_count, &row_initial_ready_count);
      row_dep_ready_queue = use_row_ready_queue ? 1 : 0;
      row_ready_initial_queued = use_row_ready_queue &&
        row_private_groups == NULL;
    } else {
      use_row_ready_queue =
        kls_prepare_row_refactor_ready_queue(
          solver, cluster_levels, &row_ready_groups, &row_ready_slots,
          &row_remaining_preds, &row_tail_groups, &row_private_groups,
          &row_ready_tail_count, &row_initial_ready_count);
      row_ready_initial_queued = use_row_ready_queue &&
        row_private_groups == NULL;
    }
  }
  if (use_row_ready_queue && !row_separator_flop_queue &&
      row_private_groups != NULL &&
      row_initial_ready_count > 0u &&
      kls_build_row_refactor_separator_private_thread_queues(
        solver, row_private_groups, row_initial_ready_count, thread_count,
        &row_owned_private_groups, &row_private_thread_ptr,
        &row_separator_private_component_count)) {
    row_private_groups = row_owned_private_groups;
    row_separator_private_queue = 1;
    row_private_group_count = row_initial_ready_count;
    row_ready_initial_queued = 0;
  }
  if (use_row_ready_queue && !row_separator_flop_queue &&
      row_private_groups != NULL &&
      row_owned_private_groups == NULL &&
      row_initial_ready_count > 0u &&
      kls_build_row_refactor_private_thread_queues(
        solver, row_private_groups, row_initial_ready_count, thread_count,
        &row_owned_private_groups, &row_private_thread_ptr)) {
    row_private_groups = row_owned_private_groups;
    row_private_group_count = row_initial_ready_count;
    row_ready_initial_queued = 0;
  }
  if (use_row_ready_queue && !row_separator_flop_queue &&
      row_private_groups == NULL &&
      row_ready_groups != NULL && row_ready_slots != NULL &&
      row_initial_ready_count > 0u &&
      kls_build_row_refactor_separator_private_thread_queues(
        solver, row_ready_groups, row_initial_ready_count, thread_count,
        &row_owned_private_groups, &row_private_thread_ptr,
        &row_separator_private_component_count)) {
    row_private_groups = row_owned_private_groups;
    row_separator_private_queue = 1;
    row_private_group_count = row_initial_ready_count;
    row_ready_initial_queued = 0;
    for (UF_long pos = 0; pos < row_initial_ready_count; ++pos) {
      atomic_store_explicit(&row_ready_slots[pos], 0u,
                            memory_order_relaxed);
    }
  }
  if (use_row_ready_queue && !row_separator_flop_queue &&
      row_private_groups == NULL &&
      row_ready_groups != NULL && row_ready_slots != NULL &&
      row_initial_ready_count > 0u &&
      kls_build_row_refactor_private_thread_queues(
        solver, row_ready_groups, row_initial_ready_count, thread_count,
        &row_owned_private_groups, &row_private_thread_ptr)) {
    row_private_groups = row_owned_private_groups;
    row_private_group_count = row_initial_ready_count;
    row_ready_initial_queued = 0;
    for (UF_long pos = 0; pos < row_initial_ready_count; ++pos) {
      atomic_store_explicit(&row_ready_slots[pos], 0u,
                            memory_order_relaxed);
    }
  }
  if (use_row_ready_queue && row_private_groups != NULL &&
      row_private_group_count == 0u) {
    row_private_group_count = row_initial_ready_count;
  }
  if (use_row_ready_queue && row_separator_flop_queue &&
      !row_dep_ready_queue && row_private_groups != NULL &&
      row_private_thread_ptr != NULL && row_private_group_mask != NULL) {
    row_ordered_private_queue =
      kls_row_refactor_private_thread_order_valid(
        solver, row_private_groups, row_private_thread_ptr,
        row_private_group_mask, thread_count);
  }
  atomic_uint *pipeline_done = NULL;
  unsigned int pipeline_generation = 0;
  if (check_pivots || row_dep_ready_queue ||
      (use_group_pipeline && !use_row_ready_queue)) {
    pipeline_done =
      ensure_egraph_pipeline_done(solver, &pipeline_generation);
    if (pipeline_done == NULL) {
      free(row_owned_private_groups);
      free(row_private_thread_ptr);
      free(row_private_group_mask);
      return -1;
    }
  }

  double **scratch =
    ensure_egraph_worker_scratch(solver, thread_count, solver->n);
  if (scratch == NULL) {
    free(row_owned_private_groups);
    free(row_private_thread_ptr);
    free(row_private_group_mask);
    return -1;
  }
  kls_egraph_refactor_pool *pool =
    ensure_egraph_refactor_pool(solver, thread_count);
  if (pool == NULL) {
    free(row_owned_private_groups);
    free(row_private_thread_ptr);
    free(row_private_group_mask);
    return -1;
  }

  trilinos_klu_l_common *common = &solver->common;
  if (scaled &&
      !trilinos_klu_l_scale((UF_long)common->scale, solver->n,
                            solver->col_ptr, solver->row_idx,
                            numeric_values, solver->numeric->Rs, NULL,
                            common)) {
    free(row_owned_private_groups);
    free(row_private_thread_ptr);
    free(row_private_group_mask);
    return 0;
  }
  if (!kls_refresh_row_refactor_offdiag_values(solver, numeric_values)) {
    common->status = TRILINOS_KLU_INVALID;
    free(row_owned_private_groups);
    free(row_private_thread_ptr);
    free(row_private_group_mask);
    return 0;
  }
  common->status = TRILINOS_KLU_OK;
  common->numerical_rank = KLS_KLU_EMPTY;
  common->singular_col = KLS_KLU_EMPTY;
  common->nrealloc = 0;
  const int native_row_panel_state =
    kls_native_row_panel_effective_state(solver);
  const UF_long native_row_panel_count_before =
    solver->row_refactor_native_row_panel_count;
  const UF_long compact_supernode_update_count_before =
    solver->row_refactor_compact_supernode_update_count;
  const UF_long partial_supernode_pipeline_count_before =
    solver->row_refactor_partial_supernode_pipeline_run_count;

  kls_egraph_refactor_shared *shared = &pool->shared;
  pthread_mutex_lock(&shared->lock);
  if (pool->active_workers != 0) {
    pthread_mutex_unlock(&shared->lock);
    free(row_owned_private_groups);
    free(row_private_thread_ptr);
    free(row_private_group_mask);
    return -1;
  }

  shared->solver = solver;
  shared->values = numeric_values;
  shared->rs = scaled ? solver->numeric->Rs : NULL;
  shared->check_pivots = check_pivots;
  shared->scale = (int)common->scale;
  shared->kernel = KLS_EGRAPH_REFACTOR_KERNEL_GENERIC;
  shared->thread_count = thread_count;
  shared->row_refactor_mode = 1;
  shared->trial_deadline_seconds =
    solver->row_trial_deadline > 0.0
      ? kls_now_seconds() + solver->row_trial_deadline
      : 0.0;
  atomic_store_explicit(&shared->trial_deadline_hit, 0,
                        memory_order_relaxed);
  shared->row_solve_mode = 0;
  shared->row_solve_view = NULL;
  shared->row_solve_work = NULL;
  shared->row_refactor_defer_value_scatter =
    kls_row_refactor_should_defer_value_scatter(solver, check_pivots);
  shared->row_refactor_lazy_value_scatter = !check_pivots;
  shared->row_refactor_blocked_trailing_update = 1;
  shared->row_refactor_native_row_panel_state = native_row_panel_state;
  shared->row_refactor_native_row_panel_active = 0;

  solver->row_refactor_native_row_panel_enabled =
    kls_native_row_panel_enabled_for_run(
      solver, native_row_panel_state, check_pivots);
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
  shared->row_pipeline_ready_queue = use_row_ready_queue;
  shared->row_pipeline_row_dep_ready_queue =
    use_row_ready_queue && row_dep_ready_queue;
  shared->row_pipeline_ready_groups = row_ready_groups;
  shared->row_pipeline_ready_slots = row_ready_slots;
  shared->row_pipeline_remaining_preds = row_remaining_preds;
  shared->row_pipeline_tail_groups = row_tail_groups;
  shared->row_pipeline_tail_count =
    use_row_ready_queue ? row_ready_tail_count : 0u;
  shared->row_pipeline_private_groups =
    use_row_ready_queue ? row_private_groups : NULL;
  shared->row_pipeline_private_thread_ptr =
    use_row_ready_queue ? row_private_thread_ptr : NULL;
  shared->row_pipeline_private_group_mask =
    use_row_ready_queue ? row_private_group_mask : NULL;
  shared->row_pipeline_private_count =
    use_row_ready_queue && row_private_groups != NULL
      ? row_private_group_count : 0u;
  shared->row_pipeline_ordered_private =
    use_row_ready_queue ? row_ordered_private_queue : 0;
  atomic_store_explicit(&shared->row_pipeline_private_pos, 0ul,
                        memory_order_release);
  atomic_store_explicit(&shared->row_pipeline_local_ready_groups, 0ul,
                        memory_order_release);
  atomic_store_explicit(&shared->row_prefactor_rows, 0ul,
                        memory_order_release);
  atomic_store_explicit(&shared->row_prefactor_deps, 0ul,
                        memory_order_release);
  atomic_store_explicit(&shared->row_prefactor_supernode_rows, 0ul,
                        memory_order_release);
  atomic_store_explicit(&shared->row_prefactor_supernode_deps, 0ul,
                        memory_order_release);
  atomic_store_explicit(&shared->row_pipeline_ready_head, 0ul,
                        memory_order_release);
  atomic_store_explicit(&shared->row_pipeline_ready_tail,
                        (unsigned long)(use_row_ready_queue
                                          ? (row_ready_initial_queued
                                               ? row_initial_ready_count
                                               : 0u)
                                          : 0u),
                        memory_order_release);
  atomic_store_explicit(&shared->row_pipeline_completed_groups, 0ul,
                        memory_order_release);
  atomic_store_explicit(&shared->next_pipeline_pos,
                        (unsigned long)(use_group_pipeline
                                          ? pipeline_group_begin : 0u),
                        memory_order_release);
  shared->pipeline_pos_end =
    use_group_pipeline ? solver->row_refactor_group_count : 0u;
  shared->cluster_level_count = cluster_levels;
  kls_record_row_refactor_run(solver, check_pivots, 1);
  kls_reset_row_refactor_compact_panel_valid(solver);
  if (use_row_ready_queue) {
    kls_record_row_refactor_ready_queue_run(solver, row_ready_tail_count);
    kls_record_row_refactor_work_ready_queue_run(solver);
    if (row_dep_ready_queue) {
      kls_record_row_refactor_partial_supernode_pipeline_run(
        solver, partial_supernode_groups, partial_supernode_rows);
    }
    kls_record_row_refactor_private_ready_groups(
      solver, row_private_group_count);
    if (row_separator_private_queue) {
      kls_record_row_refactor_separator_private_queue(
        solver, row_separator_private_component_count);
    }
    if (row_separator_flop_queue) {
      kls_record_row_refactor_separator_flop_queue(
        solver, row_separator_flop_component_count,
        row_separator_flop_private_groups,
        row_separator_flop_pipeline_groups,
        row_separator_flop_closure_groups,
        row_ordered_private_queue,
        row_separator_flop_private_threads,
        row_separator_flop_private_min_groups,
        row_separator_flop_private_max_groups,
        row_separator_flop_private_min_work,
        row_separator_flop_private_max_work);
    }
  }
  if (pipeline_done != NULL) {
    kls_record_row_refactor_done_bitmap_run(solver);
  }
  if (shared->row_refactor_defer_value_scatter) {
    kls_record_row_refactor_defer_value_scatter_run(solver);
  }

  for (int i = 0; i < thread_count; ++i) {
    pool->workers[i].shared = shared;
    pool->workers[i].x = scratch[i];
    pool->workers[i].telemetry_compact_dense_input_rows = 0u;
    pool->workers[i].telemetry_dense_segment_input_rows = 0u;
    pool->workers[i].telemetry_sparse_segment_input_rows = 0u;
    pool->workers[i].telemetry_batch_input_rows = 0u;
    pool->workers[i].telemetry_segment_target_input_rows = 0u;
    pool->workers[i].telemetry_segment_target_cleanup_rows = 0u;
    pool->workers[i].telemetry_segment_target_cleanup_entries = 0u;
  }

  kls_egraph_pool_dispatch_and_wait(pool, shared, thread_count);
  for (int i = 0; i < thread_count; ++i) {
    const kls_egraph_refactor_worker *telemetry = &pool->workers[i];
    solver->row_refactor_last_compact_dense_panel_direct_input_rows +=
      telemetry->telemetry_compact_dense_input_rows;
    solver->row_refactor_compact_dense_panel_direct_input_rows +=
      telemetry->telemetry_compact_dense_input_rows;
    solver->row_refactor_last_dense_segment_direct_input_rows +=
      telemetry->telemetry_dense_segment_input_rows;
    solver->row_refactor_dense_segment_direct_input_rows +=
      telemetry->telemetry_dense_segment_input_rows;
    solver->row_refactor_last_sparse_segment_direct_input_rows +=
      telemetry->telemetry_sparse_segment_input_rows;
    solver->row_refactor_sparse_segment_direct_input_rows +=
      telemetry->telemetry_sparse_segment_input_rows;
    solver->row_refactor_last_batch_direct_input_rows +=
      telemetry->telemetry_batch_input_rows;
    solver->row_refactor_batch_direct_input_rows +=
      telemetry->telemetry_batch_input_rows;
    solver->row_refactor_last_segment_target_input_rows +=
      telemetry->telemetry_segment_target_input_rows;
    solver->row_refactor_segment_target_input_rows +=
      telemetry->telemetry_segment_target_input_rows;
    solver->row_refactor_last_segment_target_cleanup_rows +=
      telemetry->telemetry_segment_target_cleanup_rows;
    solver->row_refactor_segment_target_cleanup_rows +=
      telemetry->telemetry_segment_target_cleanup_rows;
    solver->row_refactor_last_segment_target_cleanup_entries +=
      telemetry->telemetry_segment_target_cleanup_entries;
    solver->row_refactor_segment_target_cleanup_entries +=
      telemetry->telemetry_segment_target_cleanup_entries;
  }
  const UF_long row_local_ready_groups =
    use_row_ready_queue
      ? (UF_long)atomic_load_explicit(
          &shared->row_pipeline_local_ready_groups, memory_order_acquire)
      : 0u;
  const UF_long row_prefactor_rows =
    (UF_long)atomic_load_explicit(&shared->row_prefactor_rows,
                                  memory_order_acquire);
  const UF_long row_prefactor_deps =
    (UF_long)atomic_load_explicit(&shared->row_prefactor_deps,
                                  memory_order_acquire);
  const UF_long row_prefactor_supernode_rows =
    (UF_long)atomic_load_explicit(&shared->row_prefactor_supernode_rows,
                                  memory_order_acquire);
  const UF_long row_prefactor_supernode_deps =
    (UF_long)atomic_load_explicit(&shared->row_prefactor_supernode_deps,
                                  memory_order_acquire);
  shared->row_refactor_mode = 0;
  shared->row_refactor_blocked_trailing_update = 0;
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
  pthread_mutex_unlock(&shared->lock);
  free(row_owned_private_groups);
  free(row_private_thread_ptr);
  free(row_private_group_mask);
  if (use_row_ready_queue) {
    kls_record_row_refactor_local_ready_groups(solver,
                                               row_local_ready_groups);
  }
  kls_record_row_refactor_prefactor_run(solver, row_prefactor_rows,
                                        row_prefactor_deps);
  kls_record_row_refactor_prefactor_supernode_run(
    solver, row_prefactor_supernode_rows, row_prefactor_supernode_deps);

  if (shared->invalid || shared->pivot_rejected ||
      (shared->singular && common->halt_if_singular)) {
    solver->egraph_worker_scratch_dirty = 1;
  }
  int unfinished_tail_seed_recorded = 0;
  if (shared->pivot_rejected && pipeline_done != NULL) {
    kls_record_fast_reject_unfinished_tail_seed(
      solver, shared, shared->rejected_pivot);
    unfinished_tail_seed_recorded =
      solver->fast_reject_tail_seed_valid ? 1 : 0;
  }
  if (shared->pivot_rejected && pipeline_done != NULL && scratch[0] != NULL) {
    kls_row_refactor_refresh_missing_prefix(shared, scratch[0]);
  }
  shared->trial_deadline_seconds = 0.0;
  if (atomic_load_explicit(&shared->trial_deadline_hit,
                           memory_order_acquire)) {
    /* Trial budget exhausted: hand the refactor to the column engine
       and settle the verdict. The column engine rewrites every factor
       value, so the partially updated numeric is safe to abandon. */
    solver->row_accept_decision = -1;
    solver->row_refactor_values_ready = 0;
    solver->row_refactor_solve_direct_ready = 0;
    solver->row_refactor_solve_validated = 0;
    solver->egraph_worker_scratch_dirty = 1;
    common->status = TRILINOS_KLU_OK;
    return -1;
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
    if (prefix_ready && shared->row_refactor_defer_value_scatter &&
        !kls_scatter_row_refactor_prefix_values(solver,
                                                shared->rejected_pivot)) {
      common->status = TRILINOS_KLU_INVALID;
      return 0;
    }
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
  } else {
    common->status = TRILINOS_KLU_OK;
  }
  int published_for_solve = 0;
  if (shared->row_refactor_defer_value_scatter) {
    if (shared->row_refactor_lazy_value_scatter &&
        kls_row_refactor_should_publish_for_solve(solver)) {
      const double publish_start = kls_now_seconds();
      pthread_mutex_lock(&shared->lock);
      shared->row_refactor_mode = 0;
      shared->row_publish_mode = 1;
      kls_egraph_pool_dispatch_and_wait(pool, shared, thread_count);
      shared->row_publish_mode = 0;
      shared->row_refactor_mode = 1;
      pthread_mutex_unlock(&shared->lock);
      solver->row_refactor_values_dirty = 0;
      solver->stats.row_refactor_values_dirty = 0;
      published_for_solve = 1;
      if (solver->row_publish_experiment == 1) {
        solver->row_publish_probe_seconds =
          kls_now_seconds() - publish_start;
        solver->row_publish_experiment = 2;
      }
    } else if (shared->row_refactor_lazy_value_scatter) {
      kls_record_row_refactor_lazy_value_scatter_run(solver);
    } else if (!kls_scatter_row_refactor_l_values(solver) ||
               !kls_scatter_row_refactor_u_values(solver)) {
      common->status = TRILINOS_KLU_INVALID;
      return 0;
    }
  }
  shared->row_refactor_lazy_value_scatter = 0;
  if (scaled && !kls_parallel_refactor_permute_scale(solver)) {
    common->status = TRILINOS_KLU_INVALID;
    return 0;
  }
  solver->row_refactor_values_ready = published_for_solve ? 0 : 1;
  solver->row_refactor_solve_direct_ready = 0;
  solver->row_refactor_solve_validated = 0;
  solver->numeric_needs_refinement = 0;
  solver->solve_refine_single_shot = 0;
  kls_maybe_disable_native_row_panel_auto(
    solver, check_pivots, native_row_panel_state,
    native_row_panel_count_before, compact_supernode_update_count_before,
    partial_supernode_pipeline_count_before);
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
  const int direct_user_values = solver->refactor_direct_user_values_active;
  const int32_t *direct_user_pos = solver->refactor_input_user_pos32;
  const double *prepared_scale = solver->prepared_value_scale;
  double *prepared_values = solver->values;
  if (direct_user_values &&
      (direct_user_pos == NULL || prepared_values == NULL)) {
    return -1;
  }

  common->status = TRILINOS_KLU_OK;
  common->numerical_rank = KLS_KLU_EMPTY;
  common->singular_col = KLS_KLU_EMPTY;
  common->nrealloc = 0;

  for (UF_long k = 0; k < solver->n; ++k) {
    for (UF_long p = solver->refactor_col_ptr[k];
         p < solver->refactor_col_ptr[k + 1u]; ++p) {
      const UF_long internal = solver->refactor_input_pos[p];
      double value = direct_user_values
        ? numeric_values[(UF_long)direct_user_pos[p]]
        : numeric_values[internal];
      if (direct_user_values && prepared_scale != NULL) {
        value *= prepared_scale[internal];
      }
      if (direct_user_values) {
        /* Keep the owned internal frame current for solve-contract and cache
           consumers.  This is the serial mapped counterpart of the fused PTS
           scatter; the public values are otherwise in a different order (and
           may require the explicit match scales). */
        prepared_values[internal] = value;
      }
      x[solver->refactor_row_idx[p]] = value;
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

      if (ujk != 0.0) {
        UF_long *li = NULL;
        double *lx = NULL;
        UF_long lcol_len = 0;
        kls_klu_get_pointer(lu, lip, llen, j, &li, &lx, &lcol_len);
        kls_scatter_subtract(x, li, lx, lcol_len, ujk);
      }
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
      kls_record_fast_reject_block_suffix_tail_seed(solver, k);
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

/* Guard-reduced fixed-pivot single-block walk.  Compact input-map positions
   halve the hot index traffic, and one reciprocal per pivot replaces every L
   divide.  If the compact solve is already live, publish those same exact
   reciprocals so the post-refactor solve refresh becomes a no-op. */
static int kls_lean_single_block_map32_refactor(kls_solver *solver,
                                                double *numeric_values) {
  if (numeric_values == NULL || solver == NULL || solver->symbolic == NULL ||
      solver->numeric == NULL || solver->symbolic->nblocks != 1u ||
      solver->numeric->Udiag == NULL || solver->numeric->Lip == NULL ||
      solver->numeric->Llen == NULL || solver->numeric->Uip == NULL ||
      solver->numeric->Ulen == NULL || solver->numeric->LUbx == NULL ||
      solver->numeric->LUbx[0] == NULL || solver->numeric->Xwork == NULL ||
      !kls_build_refactor_map(solver) ||
      solver->refactor_col_ptr == NULL ||
      solver->refactor_row_idx32 == NULL ||
      solver->refactor_input_pos32 == NULL ||
      (solver->common.scale > 0 &&
       (solver->numeric->Rs == NULL ||
        solver->numeric->Pnum == NULL ||
        solver->refactor_input_oldrow32 == NULL)) ||
      (solver->common.scale <= 0 && solver->numeric->Rs != NULL)) {
    return -1;
  }

  trilinos_klu_l_common *common = &solver->common;
  trilinos_klu_l_numeric *numeric = solver->numeric;
  const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  double *restrict x = (double *)numeric->Xwork;
  double *restrict udiag = (double *)numeric->Udiag;
  double *restrict diag_recip = solver->i32solve_udiag_recip;
  UF_long *lip = numeric->Lip;
  UF_long *llen = numeric->Llen;
  UF_long *uip = numeric->Uip;
  UF_long *ulen = numeric->Ulen;
  double *lu = (double *)numeric->LUbx[0];
  const UF_long *q = symbolic->Q;
  const int32_t *restrict map_row = solver->refactor_row_idx32;
  const int32_t *restrict map_input = solver->refactor_input_pos32;
  const int32_t *restrict map_oldrow = solver->refactor_input_oldrow32;
  const int scaled = common->scale > 0;
  const double *restrict rs = (const double *)numeric->Rs;

  common->status = TRILINOS_KLU_OK;
  common->numerical_rank = KLS_KLU_EMPTY;
  common->singular_col = KLS_KLU_EMPTY;
  common->nrealloc = 0;
  if (scaled && !trilinos_klu_l_scale(
        (UF_long)common->scale, solver->n, solver->col_ptr,
        solver->row_idx, numeric_values, (double *)numeric->Rs, NULL,
        common)) {
    return 0;
  }

  for (UF_long k = 0u; k < solver->n; ++k) {
    const UF_long input_begin = solver->refactor_col_ptr[k];
    const UF_long input_end = solver->refactor_col_ptr[k + 1u];
    for (UF_long p = input_begin; p < input_end; ++p) {
      double value = numeric_values[(UF_long)map_input[p]];
      if (scaled) {
        value /= rs[(UF_long)map_oldrow[p]];
      }
      x[(UF_long)map_row[p]] = value;
    }

    UF_long *ui = NULL;
    double *ux = NULL;
    UF_long ucol_len = 0u;
    kls_klu_get_pointer(lu, uip, ulen, k, &ui, &ux, &ucol_len);
    for (UF_long up = 0u; up < ucol_len; ++up) {
      const UF_long j = ui[up];
      const double ujk = x[j];
      x[j] = 0.0;
      ux[up] = ujk;
      if (ujk != 0.0) {
        UF_long *li = NULL;
        double *lx = NULL;
        UF_long lcol_len = 0u;
        kls_klu_get_pointer(lu, lip, llen, j, &li, &lx, &lcol_len);
        kls_scatter_subtract(x, li, lx, lcol_len, ujk);
      }
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
    const double reciprocal = 1.0 / ukk;
    if (diag_recip != NULL) {
      diag_recip[k] = reciprocal;
    }

    UF_long *li = NULL;
    double *lx = NULL;
    UF_long lcol_len = 0u;
    kls_klu_get_pointer(lu, lip, llen, k, &li, &lx, &lcol_len);
    for (UF_long p = 0u; p < lcol_len; ++p) {
      const UF_long i = li[p];
      lx[p] = x[i] * reciprocal;
      x[i] = 0.0;
    }
  }
  if (diag_recip != NULL) {
    solver->i32solve_udiag_recip_fresh = 1;
  }
  if (scaled) {
    double *restrict numeric_rs = (double *)numeric->Rs;
    for (UF_long k = 0u; k < solver->n; ++k) {
      x[k] = numeric_rs[numeric->Pnum[k]];
    }
    for (UF_long k = 0u; k < solver->n; ++k) {
      numeric_rs[k] = x[k];
      x[k] = 0.0;
    }
  }
  return 1;
}

/* Compare the vendor fixed-pattern walk with the compact predecoded stream
   only at a dispatch point that has already selected serial KLU.  Both arms
   recompute the same numeric on the same values; the compact arm is retained
   only after a deterministic solve/residual certificate and only when its
   projected recurring saving repays the one-time trial twice. */
static int kls_try_compact_map32_klu_tournament(
  kls_solver *solver,
  double *numeric_values,
  UF_long *ok_out,
  int *compact_used_out) {
  KLS_SET_OPTIONAL_OUTPUT(ok_out, 0u);
  KLS_SET_OPTIONAL_OUTPUT(compact_used_out, 0);
  if (solver == NULL || numeric_values == NULL || ok_out == NULL ||
      compact_used_out == NULL || solver->symbolic == NULL ||
      solver->numeric == NULL || solver->symbolic->nblocks != 1u ||
      solver->dense_tail_cols != 0u || solver->n > (UF_long)INT32_MAX ||
      solver->nnz > (UF_long)INT32_MAX ||
      !kls_direct_klu_public_frame_capable(solver)) {
    return 0;
  }

  const int settled_low_work_compact =
    kls_low_work_single_block_policy_enabled(solver) &&
    kls_repeated_update_workload(&solver->options);
  if (solver->compact_map32_choice == 0 &&
      settled_low_work_compact) {
    /* In this cache-resident regime the compact walk performs the same
       fixed-pivot arithmetic as KLU while replacing the hot 64-bit input
       maps with predecoded i32 streams and one divide per column with a
       reciprocal.  A six-update alternating tournament costs as much as
       the complete remaining difference on such short numerics, and its
       microsecond verdict is vulnerable to timer/cache noise.  Select the
       representation from the existing capability/work contract instead;
       the normal first-update solve contract still validates the resulting
       numeric before the later direct lifecycle can use it. */
    solver->compact_map32_choice = 1;
  }
  if (solver->compact_map32_choice < 0 ||
      !kls_repeated_update_workload(&solver->options)) {
    return 0;
  }

  if (solver->compact_map32_choice > 0) {
    solver->compact_map32_probe_active = 1;
    const int compact_ok =
      kls_lean_single_block_map32_refactor(solver, numeric_values);
    solver->compact_map32_probe_active = 0;
    if (compact_ok >= 0) {
      *ok_out = (UF_long)compact_ok;
      *compact_used_out = 1;
      return 1;
    }
    solver->compact_map32_choice = -1;
    return 0;
  }

  const int state = solver->compact_map32_trial_state;
  if (state == 0) {
    const double prep_start = kls_now_seconds();
    solver->compact_map32_probe_active = 1;
    const int map_ready = kls_build_refactor_map(solver);
    solver->compact_map32_probe_active = 0;
    solver->compact_map32_trial_overhead +=
      kls_now_seconds() - prep_start;
    if (!map_ready) {
      solver->compact_map32_choice = -1;
      return 0;
    }
  }

  /* Low-work cases settle above; remaining trials alternate compact/vendor. */
  const int compact_arm = (state & 1) == 0;
  solver->common.status = TRILINOS_KLU_OK;
  const double arm_start = kls_now_seconds();
  int arm_ok = 0;
  if (compact_arm) {
    arm_ok = kls_lean_single_block_map32_refactor(
      solver, numeric_values);
  } else {
    arm_ok = (int)trilinos_klu_l_refactor(
      solver->col_ptr, solver->row_idx, numeric_values, solver->symbolic,
      solver->numeric, &solver->common);
  }
  const double arm_seconds = kls_now_seconds() - arm_start;
  if (arm_ok <= 0 || solver->common.status < 0) {
    solver->compact_map32_choice = -1;
    if (compact_arm) {
      return 0;
    }
    *ok_out = (UF_long)arm_ok;
    return 1;
  }
  const int arm = compact_arm ? 0 : 1;
  solver->compact_map32_trial_samples[arm][state / 2] = arm_seconds;

  if (state == 2) {
    const double probe_start = kls_now_seconds();
    solver->compact_map32_trial_probe_ok =
      kls_direct_klu_numeric_residual_probe(solver, numeric_values);
    solver->compact_map32_trial_overhead +=
      kls_now_seconds() - probe_start;
  }
  solver->compact_map32_trial_state = state + 1;
  *ok_out = (UF_long)arm_ok;
  *compact_used_out = compact_arm;

  if (state == 5) {
    const double compact_a = solver->compact_map32_trial_samples[0][0];
    const double compact_b = solver->compact_map32_trial_samples[0][1];
    const double compact_c = solver->compact_map32_trial_samples[0][2];
    const double vendor_a = solver->compact_map32_trial_samples[1][0];
    const double vendor_b = solver->compact_map32_trial_samples[1][1];
    const double vendor_c = solver->compact_map32_trial_samples[1][2];
    /* Use the middle observation, not the optimistic minimum.  The latter
       can promote a challenger after a single lucky cache/scheduler sample,
       which is especially harmful on short fixed-pattern factorizations. */
    const double compact_seconds =
      compact_a + compact_b + compact_c -
      fmin(compact_a, fmin(compact_b, compact_c)) -
      fmax(compact_a, fmax(compact_b, compact_c));
    const double vendor_seconds =
      vendor_a + vendor_b + vendor_c -
      fmin(vendor_a, fmin(vendor_b, vendor_c)) -
      fmax(vendor_a, fmax(vendor_b, vendor_c));
    const double remaining =
      (double)(solver->options.expected_refactorizations - 6);
    const double projected_saving =
      remaining * (vendor_seconds - compact_seconds);
    /* Below this absolute delta, timer/scheduler noise can dominate the
       complete expected run and the maximum possible gain is immaterial. */
    const int adopt_compact =
      solver->compact_map32_trial_probe_ok && compact_seconds > 0.0 &&
      vendor_seconds > 0.0 && compact_seconds < 0.98 * vendor_seconds &&
      vendor_seconds - compact_seconds >= 1.0e-6 &&
      projected_saving >
        2.0 * solver->compact_map32_trial_overhead;
    solver->compact_map32_choice = adopt_compact ? 1 : -1;
  }
  return 1;
}
