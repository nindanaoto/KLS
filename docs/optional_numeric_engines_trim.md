# Optional numeric engine removal

This batch, relative to `6b9d61c`, removes two opt-in experiments that the
paper configurations did not enable:

- Algorithm-5 speculative EGraph prefactor updates: the later-predecessor
  scan and safety helpers, applied-position bitmap branches in four kernels,
  bitmap parameters in batch/run helpers, dispatch exclusions, internal
  counters, and the environment switch.
- FP32 factor mirrors: allocation/refresh/free paths, float scalar and batch
  kernels, float stores and subtree-solve branches, adoption state,
  first-pass validation/retry/redo, and exclusive diagnostic switches.

`KLS_ENABLE_EGRAPH_ALGORITHM5_PREF_UPDATE`, `KLS_ENABLE_FP32_REFACTOR`, and
`KLS_FP32_TAIL_FLOAT_ACCUM` no longer affect execution. The README now says so.
Historical experiment reports are retained as historical evidence.

The four public EGraph prefactor statistics remain in their original ABI
positions and are explicitly zero. Smoke validation checks that property.
Ordinary double-precision refinement, residual validation, pivot repair,
GMRES/LSQR recovery, production row-prefactor workspace, narrow integer-index
mirrors, cached supernodes, plain/fused BTF dispatch, and worker affinity remain.
No matrix-name condition, production numeric threshold, or replacement engine
was introduced.

The production `src/` change removes **1,154 net lines**. In particular,
`kls_egraph_refactor.inc` drops from 12,168 to 11,398 lines (770 fewer).

## Validation protocol

Frozen binaries are compared against the pre-removal `6b9d61c` build with
matching linked-target build provenance. Before/after order alternates, with
one benchmark process at a time, on both eight-thread cache domains. Each
launch uses H100 entrywise updates of amplitude 0.001 and checks every refactor
at residual limit `1e-8`. Auxiliary library pools have one thread.

The isolated Algorithm-5 trial is in `build/algorithm5-trim-fxNvIg/`.
All 328 launches passed. TSOPF_FS_b9_c6 median paired lifecycle changes were
+0.18% at one thread, -0.74% on the 32MiB eight-thread domain, and +0.39% on
the 96MiB domain. Eight additional controls were tested on both domains.

The combined trial is in `build/optional-engines-trim-JIyt62/`, with source
diff, provenance, runner, frozen binaries, raw observations, and reductions.
It adds ASIC_100ks, ASIC_320k, ASIC_680ks, ckt11752_dc_1, circuit_4, transient,
and rajat20 to the previous controls, for 440 launches across 16 matrices.

Final Release and ASan/UBSan CTest each pass 5/5. A further Release smoke run
with all three retired switches set to 1 passes on CPUs 8-15. Removed-symbol
searches and `git diff --check` are clean. The rebuild that recompiled the
unchanged benchmark source emitted its existing pedantic function-pointer and
long-string warnings; the edited solver source emitted no warnings.

These are bounded matrix lifecycle checks, not a full paper-corpus or Xyce
SPICE rerun. They cannot establish full-corpus parity or exact performance
equality.

## Results

All **440/440 combined-trial launches passed**, exercising EGraph, row,
mapped, and KLU refactor paths. The maximum refactor residual was
`5.45371652e-9` in both versions, below the unchanged `1e-8` limit.

TSOPF_FS_b9_c6 received 20 one-thread pairs and 40 pairs per eight-thread
domain. Its median paired lifecycle changes were +0.08%, +0.56%, and +0.69%,
respectively. Neither version had an eight-thread launch above the 0.6ms
average steady-refactor cutoff used in the earlier affinity diagnosis.
After steady medians were 0.412ms and 0.416ms on the two domains.

Other matrices received four pairs per domain. Positive means slower:

| Matrix | 32MiB, 8 threads | 96MiB, 8 threads |
| --- | ---: | ---: |
| TSOPF_FS_b9_c1 | -1.20% | -1.53% |
| TSOPF_RS_b9_c6 | -1.77% | -0.23% |
| onetone1 | -0.06% | -1.79% |
| twotone | +0.48% | +0.09% |
| rajat25 | +0.12% | +0.26% |
| adder_dcop_01 | -1.50% | -8.11% |
| 1138_bus | -1.23% | -4.74% |
| bips98_1142 | -3.80% | -1.42% |
| ASIC_100ks | -0.91% | +0.75% |
| ASIC_320k | +0.13% | -0.21% |
| ASIC_680ks | -1.32% | +1.40% |
| ckt11752_dc_1 | -0.03% | +0.01% |
| circuit_4 | +0.64% | +1.17% |
| transient | +0.27% | -0.73% |
| rajat20 | -0.04% | -0.27% |

The three largest absolute control differences received independent 20-pair
rechecks in `build/optional-engines-recheck-bXqvQi/`. All 120 launches passed.
Changes were -1.97% for adder_dcop_01 on 96MiB, -0.62% for 1138_bus on 96MiB,
and +0.31% for bips98_1142 on 32MiB. The large apparent gains did not repeat
at their initial magnitudes.

In total, **888/888 benchmark launches passed**: 328 isolated Algorithm-5,
440 combined, and 120 combined rechecks. No large repeatable regression was
established. Both removals are retained. The final Release binary matches
the frozen combined and recheck binaries (SHA256
`b9b63a43c27fcf41db106480ca4632938d879a7f7a3ef92260e216b96c792b78`).
