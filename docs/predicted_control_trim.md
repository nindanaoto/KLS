# Remove unused predicted-factor controls

Baseline: `113143d`. Remove eight controls unused by tests, scripts and
benchmark drivers:

- `KLS_DISABLE_PARALLEL_PREDICTED_BUILD`
- `KLS_PREDICTED_FORCE_SORT`
- `KLS_PREDICTED_NUDGE_SLOTS_8K`
- `KLS_PREDICTED_NUDGE_ROUNDS`
- `KLS_PREDICTED_NUDGE_SIGMA_SCALE`
- `KLS_PREDICTED_PROBE_ITERS`
- `KLS_STRICT_PREDICTED_REFINEMENT`
- `KLS_BLOCK_DENSE_FLOOR`

Keep default parallel construction and its allocation fallback, sorted-pattern
declaration, three fill rounds, 1024 nudge slots, unchanged sigma calculation,
probe limits (10 for block trials, three otherwise), 0.5 stall factor, raw
residual refinement decision and structural dense-row floor. Delete only
environment overrides and parsers; no default threshold or numeric kernel
changes. Net reduction: 47 source lines. Removed overrides are no longer
supported; the public C API is unchanged. Earlier diagnostic ablations using
these names describe historical executables, not current runtime options.

## Validation

Release and ASan builds pass without reported warnings. Release CTest passes
6/6 in 4.13 s; ASan passes 6/6 in 17.22 s. `git diff --check` passes.
Executable text decreases by 896 bytes; data and BSS are unchanged.

All 579 focused H100 launches pass the provenance/result audit. Paired
lifecycle changes range from -1.726% to +1.190% versus accepted and from
-4.092% to +0.214% versus cumulative baseline `372f9d1`. No case crosses
the 2.5% regression review threshold. No alignment change is needed.

Audits verify exact jobs, commands, source revision/diff, build provenance,
binary/matrix hashes, exits, residuals and recomputed medians. Source and
binaries remain fixed throughout timing. These focused results do not
establish full-paper or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/predicted-control-trim-*`.
Runner: `validate-predicted-control-trim.py`. Audit: `audit-focused-trim.py`.
