# Remove unused batching, map and row-solve controls

Baseline: `9aac172`. Remove eight controls unused by tests, scripts and
benchmark drivers:

- `KLS_DISABLE_TINY_SINGLETON_SOLVE`
- `KLS_DISABLE_BATCH_CONSUME`
- `KLS_DISABLE_LEAN_BTF_CONTIGUOUS_OFF_RUNS`
- `KLS_DISABLE_LEAN_BTF_FLAT_OFF_MAP`
- `KLS_DISABLE_MAPPED_NATIVE_SHORT_L`
- `KLS_DISABLE_LARGE_WEAK_PTS_VALUE_MIRROR_ELISION`
- `KLS_ROW_SOLVE_FORCE_PARALLEL`
- `KLS_ROW_SOLVE_DISABLE_SPARSE_PARALLEL`

Delete the cached batching-disable helper and simplify its callers, including
one duplicate solver-null check in the fusion predicate. Keep default
batching and all numerical kernels, singleton eligibility, compact-map
allocation/range fallbacks, short-L capability checks and solve-certificate
requirements for omitting prepared-value mirrors. Keep parallel row-solve
coverage and synchronization-work admission and its sparse-level schedule.
Only the force/disable branches are removed.

The removed overrides are no longer supported; public C API unchanged.
Net reduction: 33 source lines across three files. Update the historical
batch-counter cleanup note after timing audits.

## Validation

Final Release and ASan builds pass without reported warnings. Release CTest
passes 6/6 in 4.10 s; ASan passes 6/6 in 17.13 s. `git diff --check` passes.
Executable text decreases by 3,664 bytes; data and BSS are unchanged.

All 579 focused H100 launches pass the provenance/result audit. Paired
lifecycle changes range from -2.224% to +1.060% versus accepted and from
-2.834% to +1.068% versus cumulative baseline `372f9d1`. No case reaches
the 2.5% review threshold. No alignment fix is needed.

Audits verify exact jobs, commands, source revision/diff, build provenance,
binary/matrix hashes, exits, residuals and recomputed medians. Source and
binaries remain fixed throughout timing. These focused results do not
establish full-paper or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/batch-map-control-trim-*`.
Runner: `validate-batch-map-control-trim.py`. Audit: `audit-focused-trim.py`.
