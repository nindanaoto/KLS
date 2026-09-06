# Diagnostic and dead-policy cleanup (2026-09-06)

First batch of the expanded specialization audit, relative to `db9ff4c`.
Larger optional engines and threaded preparation variants are not removed
in this batch, so their effects can be measured separately.

## Removed

- The old opt-in scaled-lean rejection heuristic
  (`KLS_ENABLE_GENERIC_SCALED_LEAN_DECLINE`). The separate production
  measured-work policy and its disable control remain.
- The trusted-row residual-check bypass
  (`KLS_ENABLE_CONTRACT_CERTIFIED_TRUSTED_ROW_SOLVE`) and its exclusive
  predicates. Ordinary per-solve residual checking and refinement remain.
- Padded-panel diagnostic verification (`KLS_VERIFY_PADDED_PANELS`), not
  the automatically admitted padded-panel engine.
- Refactor microsecond tracing (`KLS_TRACE_REFACTOR_US`), exclusive private
  counters, and diagnostic exclusions from settled fast paths.
- U-diagonal, block-trial matrix, block-permutation, and separator-pipeline
  row dumps (`KLS_DUMP_UDIAG`, `KLS_DUMP_BLOCK_TRIAL`,
  `KLS_DUMP_BLOCK_PERMS`, `KLS_DUMP_SEP_PIPE_ROWS`).
- Unreachable AUTO policy following Algorithm-5's unconditional rejection,
  and the helper's now-unused solver argument. The opt-in Algorithm-5
  executor itself is deferred to a separate cleanup batch.

The source change removes 207 net lines across three implementation files.
There is no public C API or statistics-structure change. The private solver
structure loses its diagnostic counters. No padding or new specialization
was introduced. The paper's EGraph/panel ablation controls remain.

## Validation protocol

Release and ASan/UBSan CTest each pass 5/5. Both builds complete without new
warnings, removed-symbol searches are clean, and `git diff --check` passes.

Artifacts are in `build/diagnostic-trim-CK6sJs/`. The before binary matches
the previous cleanup's final binary byte-for-byte (SHA-256
`3c6255f9160a8f8496a6575b71b72dfe5a874b1c368fa2ed8fbd24c67eb80c43`).
Its archived build provenance is compared against the rebuilt after binary
for all linked targets, including optimized SPRAL. The benchmark freezes
both binaries and records source diffs, hashes, commands, and observations.

Each launch uses 100 systems with entrywise changes of amplitude 0.001 and
verification of every refactor at residual limit `1e-8`. Auxiliary library
pools use one thread. Before/after launches alternate, with one solver at
a time and fixed affinity. These are bounded regression checks, not a
restart of the full paper campaign.

## Paired performance checks

All 168 initial launches passed numerical verification. Positive lifecycle
changes below mean slower. TSOPF_FS_b9_c6 has 12 pairs per eight-thread
domain; the other matrices have four pairs per domain.

| Matrix | 32MiB LLC, 8 threads | 96MiB LLC, 8 threads |
| --- | ---: | ---: |
| TSOPF_FS_b9_c6 | +0.01% | +1.83% |
| TSOPF_FS_b9_c1 | -0.82% | -3.86% |
| TSOPF_RS_b9_c6 | +0.71% | -1.14% |
| onetone1 | +0.37% | +2.50% |
| twotone | +0.73% | +0.82% |
| rajat25 | +1.65% | +0.52% |

The initial 20-pair one-thread TSOPF_FS_b9_c6 check measured -3.78%:
median lifecycle 212.159ms before versus 204.247ms after, and first refactor
80.441ms versus 72.767ms. This approximately recovers the previous cleanup's
one-thread penalty; it does not identify that penalty's lower-level cause.

An independent 20-pair one-thread recheck also passed all 40 launches and
measured -3.98%. Median lifecycle was 212.953ms versus 204.633ms; first
refactor was 80.387ms versus 72.964ms. Steady refactor was 0.988ms versus
0.981ms. In total, all 208 benchmark launches passed. The recheck metadata,
raw observations, and phase summary are archived alongside the initial run.

These bounded measurements do not establish full-corpus parity or guarantee
that the small control-matrix differences are noise. Optional engines and
threaded preparation remain for subsequent, separately measured batches.
