# KLS Performance Plan

Goal (papers-union test set, medium 93 ∪ large 17 manifests, paired-measured,
SPICE-cycle metric `analysis + init_factor + solve + 99*(refactor+solve)` at 4
threads): **generally and meaningfully faster than both CKTSO and SubtreeLU,
worst case within 2x of both.**

NOTE: this container is not the machine that produced suite #9; absolute
ratios re-baselined here.  Suite #10 (post busy-wait pool + AMF, commit
3a0bf56) is running into `results/suite10_*.jsonl` — read it before starting
new work.  CKTSO needs `cktso.lic` next to `libcktso.so` or every row fails.

## Landed this round (2026-07-06, commits cc6b907..3a0bf56)

- **Supernodal prototype rebuilt in-tree** as `bench/supernodal_bench.c` +
  `bench/snb_refactor_impl.h` (rounds 2–7 of task #12).  Serial t1 vs the
  landed-kernel replica, same-moment: mc2depi 1.42–1.70x, ASIC_320ks
  1.01–1.18x (+9% more with fp32-L), rajat25 1.05–1.12x.  Design verdicts in
  memory note `kls-supernodal-findings` (exact U-row sublists, residual-gated
  verify, row-major GEMM temp, per-edge dst lists, postorder fill gate,
  fp32-L mirror; snode-pool parallelism is flat without pipelining).
- **Busy-wait BTF refactor pool** (queued #2): atomics/spin dispatch;
  coupled t4 refactor −5%, rajat03 neutral.
- **AMF ordering** (task #11): deficiency scoring inside vendored AMD,
  `--ordering amf` + auto-promotion when AMD wins the base round and AMF's
  estimate is ≥5% better.  Forced-AMF: ASIC_320ks −24% flops/−18% refactor;
  auto currently still prefers METIS there, so suite impact rides on
  AMD-selected rows.
- **Denoise pass** (queued #1): this-machine KLS/SubtreeLU min-cycle ratios:
  TSOPF_FS_b9_c1 3.39, mac_econ 2.81, coupled 2.12, rajat03 2.12,
  rajat28 1.87, transient 1.38.

## Task #12 — in-tree port of the supernodal refactor (NEXT, the big one)

Port bar (1.3x on 2/3 proxies) formally unmet — mc2depi 1.7x but ASIC ~1.2x,
rajat25 ~1.1x — however nothing regresses and the dominant gap rows
(rajat31, memchip, circuit5M, Freescale1) are mc2depi-class giants where the
win is 1.5–1.7x serial.  Port refactor path only, per-matrix gated like every
other KLS feature.  Design (validated in the prototype):

1. Build panels after first factor+sort where `snode_prepared`-style guards
   pass (mirror `kls_maybe_prepare_snode_panels`, kls.c:85765): greedy zeta
   amalgamation (zeta 0.2 default), exact per-edge U-row sublists, per-edge
   scatter-dst lists, panel = [U rows | diag | below rows] column-major.
2. Refactor execution replaces the per-column loop for prepared matrices:
   memset panel → A-scatter program → edges (scalar path for sub_count <
   WNARROW≈8, else gather/TRSM/writeback + GEMM into row-major temp + row
   scatter) → fused no-pivot getrf.  No prefetch in panel scatters (measured
   loss).  Keep the two kernel instantiations in separate TUs (codegen
   interference, ~16%).
3. Acceptance: first supernodal refactor is compared against the landed
   kernel (interleaved timing or reconstruction-residual check + time); fall
   back permanently on loss, like block-ordering guards.
4. Etree postorder behind a fill-acceptance gate (re-analyze via
   analyze_given with composed Pnum/Q; reject on lnz growth — rajat25 +767%).
5. Parallel: run panels through the EGraph pool with column-level pipelining
   (prototype showed snode-granularity task claiming alone is flat).
6. fp32-L mirror written in the getrf epilogue, behind the existing fp32
   validation machinery (full-fp32 panels are numerically unsafe on circuit
   value ranges).
7. Supernodal solve (round 6) after the refactor path lands, on the same
   panels — solve is charged 99x; 10–20% of cycle on mc2depi/memchip/G3/
   rajat31.

## Queued behind the port

1. **Predicted-METIS init race** (measured ~140s off rajat31's init on the
   old machine; matters once the port shrinks its refactor term).
2. **mac_econ init/analysis trim** (48s one-time vs SubtreeLU 5.3s on the
   old machine) — mac_econ's remaining gap is init + elimination structure;
   METIS stays its ordering (AMD-family, incl. AMF, times out there).
3. Re-check AMF's suite-wide effect in suite #10; consider widening the
   promotion beyond AMD-won rounds only if estimate quality proves out.

## Settled questions (do not revisit without new evidence)

- GEMM retrofits on KLU scattered storage lose to the landed batch kernel.
- TSOPF_FS_b9_c1 block detection cannot crack the fused 9-bus blocks; it is
  a task-#12-class row.
- mac_econ orderings: METIS best available (COLAMD 110M fill; AMD, Scotch,
  and AMF all time out / lose).
- Udiag/value comparisons on mc2depi/rajat25 are chaotic for ALL kernels —
  verify by (noise-floored componentwise reconstruction) residual.
- b9_c6/case9/c19/c30, gemat11/12: cleared by block ordering + plain-KLU
  trial + small-class trims; block ordering default-on with guards.
- Full-fp32 panels; fused GEMM-scatter; scatter prefetch in panel kernels;
  snode-granularity-only parallelism: all measured losses in the prototype.

## Measurement discipline

- `scripts/run_paired_suite.sh` PASSES≥2, medium TIMEOUT=200 / large 900;
  scorer = per-matrix min cycle across passes; never rebuild `build/` while
  a suite runs (use a second build dir); only same-run ratios are comparable
  across sessions.  Machine drifts ~20% thermally: A/B claims require
  interleaved same-moment runs (landed-vs-snb inside supernodal_bench is the
  drift-immune comparison).
