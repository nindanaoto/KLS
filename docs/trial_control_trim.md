# Remove unused trial controls

Baseline: `da5b0dd`. Remove nine controls unused by tests, scripts and
benchmark drivers:

- `KLS_LEAN_CHOICE`
- `KLS_DISABLE_LEAN_STEADY_REAUDIT`
- `KLS_DISABLE_DECLINED_LEAN_STEADY_REAUDIT`
- `KLS_DISABLE_DIRECT_KLU_TOURNAMENT`
- `KLS_LOW_WORK_BTF_DIRECT_PREWARM`
- `KLS_COMPACT_DENSE_SPIKE_DIRECT_PREWARM`
- `KLS_DISABLE_LEAN_PROBE`
- `KLS_DISABLE_LEAN_COLUMN_CONFIRMATION`
- `KLS_DISABLE_BATCH_FLOOR_CLOSE_REAUDIT`

Delete forced lean-selection parsing and share the identical default selected
and declined re-audit admission condition. Replace the configurable prewarm
loop with its default single call; retain status propagation, direct-value
activation and charged elapsed time. Keep automatic selection, all numerical
engines, recovery, structural guards and timed comparisons. These environment
overrides are no longer supported; the public C API is unchanged.
Net reduction: 42 source lines.

## Validation

Release and ASan builds pass without reported warnings. Release CTest passes
6/6 in 4.04 s; ASan passes 6/6 in 17.13 s. Executable text decreases by
1,152 bytes; data and BSS are unchanged. `git diff --check` passes.

All 579 focused H100 launches pass the provenance/result audit. Paired
lifecycle changes range from -1.718% to +1.737% versus accepted and from
-4.108% to +0.790% versus cumulative baseline `372f9d1`. No case crosses
the 2.5% regression review threshold. No alignment change is needed.

Audits verify exact jobs, commands, source revision/diff, build provenance,
binary/matrix hashes, exits, residuals and recomputed medians. Source and
binaries remain fixed throughout timing. These focused results do not
establish full-paper or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/trial-control-trim-*`.
Runner: `validate-trial-control-trim.py`. Audit: `audit-focused-trim.py`.
