# Remove unused lean-pool and scale controls

Baseline: `de7af76`. Remove seven controls unused by tests, scripts,
benchmark drivers and existing documentation:

- `KLS_DISABLE_DECISIVE_LEAN_LIFECYCLE_AUDIT`
- `KLS_LEAN_PARALLEL_THREADS`
- `KLS_ENABLE_LEAN_PARALLEL`
- `KLS_DISABLE_SCALED_FUSED_SCALE`
- `KLS_DISABLE_SCALED_FUSED_DIRECT_RS`
- `KLS_DEFER_ROW_GROUP_THREAD_SLICES`
- `KLS_DISABLE_ROW_REFACTOR_SEPARATOR_FLOP_QUEUE`

Keep default lifecycle-audit admission, packed-row thread selection and
width capping, parallel-work admission, group-thread-slice preparation and
separator flop-queue construction. Keep all allocation/numerical fallbacks.

Fused scaling previously selected mode 1 only when the direct-Rs override
was disabled. Without that override the producer emits only 0 or 2; shared
state copies this value and the exact-worker wrappers explicitly pass 2.
Remove the two mode-1-only permuted scale stores, retaining direct row
stores in the fused arm and the original non-fused reads and permutation.
No arithmetic loop order or default structural threshold changes.

The removed environment overrides are no longer supported; the public C
API is unchanged. Net reduction: 33 source lines in one file.

## Validation

Release and ASan builds pass without reported warnings. Release CTest passes
6/6 in 4.12 s; ASan passes 6/6 in 17.20 s. `git diff --check` passes.
Executable text decreases by 1,280 bytes; data and BSS are unchanged.

All 579 focused H100 launches pass the provenance/result audit. Paired
lifecycle changes range from -1.137% to +1.533% versus accepted and from
-3.607% to +1.525% versus cumulative baseline `372f9d1`. No case reaches
the 2.5% review threshold. No alignment fix is needed.

Audits verify exact jobs, commands, source revision/diff, build provenance,
binary/matrix hashes, exits, residuals and recomputed medians. Source and
binaries remain fixed throughout timing. These focused results do not
establish full-paper or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/lean-pool-control-trim-*`.
Runner: `validate-lean-pool-control-trim.py`. Audit: `audit-focused-trim.py`.
