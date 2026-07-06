# KLS Performance Plan

Goal (papers-union test set, medium 93 ∪ large 17 manifests, paired-measured,
SPICE-cycle metric `analysis + init_factor + solve + 99*(refactor+solve)` at 4
threads): **generally and meaningfully faster than both CKTSO and SubtreeLU,
worst case within 2x of both.**

## Current standings (paired suite #9, 2026-07-06, commit 1246c73)

- vs CKTSO: geomean **0.838**, wins 53/106, over-2x = 6
  (rajat03 3.32†, ASIC_320ks 2.13, transient 2.09†, pre2 2.07, rajat28 2.01†,
  mac_econ 2.00)
- vs SubtreeLU: geomean **0.963** (first sub-1.0), wins 54/104, over-2x = 11
  (rajat31 2.94, TSOPF_FS_b9_c1 2.93, circuit5M_dc 2.83, mac_econ 2.69,
  memchip 2.59, Freescale1 2.54, rajat03 2.37†, ASIC_320ks 2.32, coupled
  2.14†, rajat28 2.13†, circuit5M 2.07)
- † = 0.05–3.5s cycles that churn ±30% between suites on identical binaries;
  re-verify with PASSES=3–4 before spending engineering on them.

Every non-noise row above traces to one mechanism: the refactor kernel is
2–3x behind SubtreeLU on KLU's scattered column storage.  Three retrofit
designs (dense mirror + GEMV, joint consumer-run, compact-row microkernel)
were prototyped and measured slower than the landed per-column batch kernel —
the fix requires supernodal storage (task #12 below).

## Task #12 — supernodal dense-panel storage rebuild (IN PROGRESS)

Prototype: `scratchpad/supernodal_bench.c` (round 2; round-1 backup
`supernodal_bench_r1.c.bak`).  Amalgamated supernodes (greedy zeta-padding,
wmax 32), COLUMN-major dense L panels as primary storage, per-edge UB gather
+ unit-lower TRSM + GEMM into row-major temp + column scatter, fused
no-pivot getrf/TRSM, hybrid scalar path for narrow supernodes (w < 4,
contiguous panel-column reads).  Verified vs klu_l_refactor to roundoff.

Serial t1 vs the landed production kernel (same-day interleaved):

| shape | round 1 | round 2 | note |
|---|---|---|---|
| mc2depi (long chains, avg_w 4) | 1.06x | **1.25x WIN** | giant-class proxy |
| ASIC_320ks (avg 8-entry cols) | 0.87x | 0.77x | structure resists amalgamation |
| rajat25 (spike rows, 640-entry cols) | 0.38x | 0.44x | wide-path edge overheads |

### Remaining rounds

- **Round 3 (started):** memset elimination landed (first producer column
  initializes the temp) — measured NOT to be rajat25's bottleneck.  Next
  diagnostic in flight: waste counters (GEMM rows whose scatter misses the
  consumer pattern → edge-level padding).  If hit-rate is low, precompute
  per-edge row intersection lists (the "refactor map" analogue) or accumulate
  directly into W columns per producer column with hoisted positions.  Also:
  UB gather/writeback are strided W accesses — hoistable.  Tune WNARROW
  (CLI arg 4) per shape.
- **Round 4: etree-postorder amalgamation.**  The greedy merge only takes
  *consecutive* columns; postordering the elimination tree (children adjacent
  to parents) is the standard way to widen supernodes from the same pattern.
  Target: ASIC_320ks avg width 2.9 → 6+, flipping its 0.77x.
- **Round 5: parallelize** at supernode granularity (the EGraph pool pattern;
  our t4 scaling is 2.6x vs SubtreeLU's ~3.2–3.9x; t8 shows 4.8x headroom).
- **Round 6: supernodal triangular solve** on the same panels (solve is
  charged 99x; 10–20% of cycle on mc2depi/memchip/G3/rajat31).
- **Round 7: FP32 panels** (dense layout doubles SIMD width in GEMM; the
  existing fp32 refinement/validation machinery applies).
- **Porting bar:** ≥1.3x over the landed kernel on 2 of the 3 proxy shapes
  before in-tree work begins.  In-tree port order: refactor path first
  (pattern fixed, no pivoting complications), factor and solve after.

## Task #11 — AMF ordering (AFTER #12 rounds, per direction 2026-07-06)

Approximate-minimum-fill ordering candidate: modify vendored AMD's
quotient-graph pivot scoring from approximate degree to approximate
deficiency; wire as `KLS_ORDERING_AMF` into the auto-ordering competition
(`choose_symbolic_for_pattern` candidates + promotion paths); validate on
both manifests (must only win where estimates say so).

Evidence for: mac_econ carries 1.46x SubtreeLU's flops at equal fill; the
rajat2x mid-family has similar ordering slack.  Caution: every AMD-family
variant tested on mac_econ (COLAMD, Scotch) was *worse* than METIS, and
SubtreeLU's advantage there is elimination *structure* at equal fill — a
minimum-fill score may not capture it.  Re-check the target list against a
denoised suite before starting.

## Queued behind those

1. **PASSES=3–4 denoising suite** for the † rows (cheap; run while
   prototyping — machine is otherwise idle).
2. **Busy-wait BTF pool** for small fragmented-BTF matrices (rajat03,
   coupled): convert `kls_refactor_pool` worker (~line 20250) + dispatch
   (~20544) to the atomics/spin pattern the EGraph pool uses; blocks are
   independent, dispatch latency dominates at 1ms scale.
3. **Predicted-METIS init race** (measured: ~140s off rajat31's init; only
   matters once #12 shrinks its refactor term — sequence after the port).
4. **mac_econ init/analysis trim** (48s of one-time cost vs SubtreeLU's 5.3s)
   — same dependency on #12 for the row to cross 2x.

## Settled questions (do not revisit without new evidence)

- GEMM retrofits on KLU scattered storage: three designs measured slower
  than the landed batch kernel (details in memory: dense-endgame-design).
- TSOPF_FS_b9_c1 block detection: its 9-bus blocks stay fused after dense
  vertex removal at any floor (interiors interconnect; one component of
  2447/2454).  It is a #12-class row, not a detection-class row.
- mac_econ orderings: METIS best available (COLAMD 110M fill; AMD/Scotch
  timeout).
- Udiag comparisons on mc2depi/rajat25 are chaotic for ALL kernels
  (summation order on ill-conditioned chains; landed kernel shows 1.58e+0
  too).  Verify those matrices by residual, not Udiag.
- b9_c6/case9/c19/c30, gemat11/12: cleared by block ordering + plain-KLU
  trial + small-class trims (suites #8/#9); block ordering is default-on
  with acceptance guards.

## Measurement discipline

- `scripts/run_paired_suite.sh` PASSES≥2, medium TIMEOUT=200 / large 900;
  scorer = per-matrix min cycle across passes; never rebuild `build/` while
  a suite runs; only same-run ratios are comparable across sessions.
- Machine drifts ~20% thermally: A/B claims require interleaved same-moment
  runs.
