# Remove unused affinity and SIMD controls

Baseline: `ed83cdc`. Remove nine controls unused by tests, scripts and
benchmark drivers:

- `KLS_SNODE_MIN_BATCH_OVERRIDE`
- `KLS_SNODE_MIN_BATCH_WORK_OVERRIDE`
- `KLS_AVX512_SCATTER`
- `KLS_AVX512_SCATTER_MIN_LENGTH`
- `KLS_DISABLE_GENERIC_LEAN_AFFINITY`
- `KLS_ENABLE_GENERIC_LEAN_AFFINITY`
- `KLS_ENABLE_COMPACT_LLC_AFFINITY`
- `KLS_DISABLE_COMPACT_LLC_AFFINITY`
- `KLS_DISABLE_PERFORMANCE_CORE_AFFINITY`

Keep automatic per-refactor batch floors and their defaults (3 and 192),
runtime AVX-512 feature detection and scalar fallback, and the default
1024-entry scatter cutoff. Delete three cached settings and the now-constant
cutoff helper. Keep topology discovery, caller cpuset constraints, automatic
worker placement, and lean scheduling profitability checks. The removed LLC
enable flag only caused extra topology-cache warming; it did not enable an
otherwise unreachable placement path. No numerical kernel is removed.
Net reduction: 48 source lines. Removed overrides are no longer supported;
the public C API is unchanged.

## Validation

Release and ASan builds pass without reported warnings. Release CTest passes
6/6 in 4.04 s; ASan passes 6/6 in 17.12 s. `git diff --check` passes.
Executable text decreases by 5,360 bytes; data and BSS are unchanged.

All 579 focused H100 launches pass the provenance/result audit. Paired
lifecycle changes range from -0.901% to +1.663% versus accepted and from
-2.898% to +1.101% versus cumulative baseline `372f9d1`. No case crosses
the 2.5% regression review threshold. No alignment change is needed.

Audits verify exact jobs, commands, source revision/diff, build provenance,
binary/matrix hashes, exits, residuals and recomputed medians. Source and
binaries remain fixed throughout timing. These focused results do not
establish full-paper or cross-machine equivalence, including hybrid-core
placement or execution on CPUs without AVX-512.

Artifacts: `build/prep-trim-repair-dBcoha/affinity-simd-control-trim-*`.
Runner: `validate-affinity-simd-control-trim.py`. Audit: `audit-focused-trim.py`.
