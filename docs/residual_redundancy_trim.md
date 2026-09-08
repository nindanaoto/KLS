# Residual and lifecycle redundancy cleanup

Baseline: `da80097`, including the transformed-residual TSOPF recovery and
the earlier SNB SIMD repair. Both are retained.

## Changes

- Remove refinement-state resets in `clear_matrix` that duplicate its prior
  call to `free_solve_refine_workspace`. The owning helper still clears the
  pointers, row boundaries, and adaptive residual state.
- Establish workload/pool eligibility before choosing the CSR crew width.
  Remove the thread-count helper's redundant eligibility parameter and
  guards; it has one caller. All thread-selection thresholds and overrides
  retain their original meaning for eligible calls.
- Replace the panel-freeing `__LINE__` macro wrapper with a normal function.
  Only the diagnostic caller-line message is removed; numerical cleanup and
  other supernode diagnostics remain.
- Remove pool/crew checks immediately after successful CSR preparation.
  Preparation already proves the pool exists, its crew is complete and in
  range, and the cached CSR crew matches the selected width. The mutex and
  active-worker acquire check remain unchanged.
- Share the two transformed CSR row-loop bodies with a local, immediately
  undefined macro. Its two expansion sites retain the original 16-/32-bit
  column types, arithmetic order, and stats accumulation. The separate
  untransformed 32-bit loop remains; no runtime index-width branch is added
  inside the nonzero loop.

Net library-source reduction: 57 lines across `kls.c` and
`kls_egraph_refactor.inc`. No solver fields, accuracy contracts, numerical
thresholds, or residual engines are removed.

## Validation

Release and ASan/UBSan CTest each pass 6/6, including the transformed residual
test's two index widths, changed values, transpose and multiple-RHS paths.
`git diff --check` passes.

Artifacts are under `build/prep-trim-repair-dBcoha/residual-slim*`.
The focused comparison uses alternating-order pairs against the frozen
`da80097`-equivalent executable (`transformed-final-kls_bench`). The runner's
`original` label refers to this post-fix baseline, not pre-cleanup `254ea42`.
Build provenance is checked, and commands, binary/matrix hashes and source
diff are recorded. H100 uses entrywise changed refactors, automatic policy,
and 1e-8 residual verification. No timed run overlaps builds or tests.

The focused campaign completes with 196/196 valid launches across nine
matrices and eleven configurations, with 24 pairs per eight-thread TSOPF
CPU set. A further 80/80 launches recheck two short controls, for **276/276
valid launches** in total. This is not a full paper or Xyce rerun.

Median paired H100 lifecycle changes versus `da80097` (negative is faster):

| Matrix | Configuration | Change |
|---|---|---:|
| TSOPF_FS_b9_c6 | t8_32 | -0.88% |
| TSOPF_FS_b9_c6 | t8_96 | -0.75% |
| TSOPF_FS_b9_c6 | t1 | -0.09% |
| TSOPF_RS_b9_c6 | t8_32 | 2.23% |
| bcircuit | t8_32 | 0.42% |
| rajat03 | t8_96 | -1.42% |
| bips98_1142 | t8_96 | -0.86% |
| 1138_bus | t8_96 | 1.19% |
| circuit_4 | t8_32 | 0.34% |
| onetone1 | t1 | -0.66% |
| transient | t1 | -0.40% |

`t1` uses CPU 8; `t8_32` uses CPUs 8-15; `t8_96` uses CPUs 0-7.

The follow-up uses 24 pairs for TSOPF_RS_b9_c6 / t8_32 and 16 pairs for
1138_bus / t8_96. Their median changes are +1.03% and +1.83%, respectively.
These are small measured slowdowns, not a proof of identical performance.
The main TSOPF recovery is retained on both CPU sets; the large serial
onetone1 and transient controls remain slightly faster.

Candidate executable SHA-256:
`926f279cf46b4a14bc678df26b7de1d9815b8602bc14f548f6b543047d59f02b`.
