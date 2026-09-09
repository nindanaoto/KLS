# Remove unused factor-lifecycle controls

Baseline: `e200e0d`. Remove thirteen controls unused by tests, scripts,
benchmark drivers and existing documentation:

- `KLS_PRESTATIC_DEFER`
- `KLS_DISABLE_CYCLE_TRIAL_DEFERRAL`
- `KLS_DISABLE_GENERIC_ND_UNSCALED_FIRST`
- `KLS_DISABLE_PREDICTED_FULL_FACTOR_BUDGET_SEED`
- `KLS_DISABLE_UNBOUNDED_GENERIC_ND_NUMERIC_TRIAL`
- `KLS_DISABLE_DENSE_TAIL`
- `KLS_RACE_FULL_JOIN`
- `KLS_DISABLE_GENERIC_LEAN_PREWARM`
- `KLS_DISABLE_GENERIC_TIGHT_PIVOT_DEFERRAL`
- `KLS_DISABLE_OVERLAPPED_COMPACT_PATTERN`
- `KLS_DISABLE_OVERLAPPED_COMPACT_SOLVE`
- `KLS_DISABLE_LOW_WORK_DIRECT_DEFERRED_PREP_SKIP`
- `KLS_DISABLE_DIRECT_LEAN_PREP_OVERLAP`

Keep default ordering retries and portfolio admission, measured full-factor
repair-budget seeding, METIS join selection, prewarming, tight-pivot trial
deferral and preparation overlap. Keep all race exclusions, joins, allocation
failures and numerical fallback paths. Retain the test-used unchanged-input
control and synchronous-preparation control. Correct the obsolete deferred
preparation comment referring to the removed vendor-first preflight.

Removed environment overrides are no longer supported; public C API is
unchanged. Net reduction: 20 source lines in `src/kls.c`.

## Validation

Release and ASan builds pass without reported warnings. Release CTest passes
6/6 in 4.19 s; ASan passes 6/6 in 17.10 s. `git diff --check` passes.
Executable `size` text, data and BSS totals are unchanged.

All 579 focused H100 launches pass the provenance/result audit. Paired
lifecycle changes range from -2.054% to +0.861% versus accepted and from
-3.440% to +1.181% versus cumulative baseline `372f9d1`. No case reaches
the 2.5% review threshold. No alignment fix is needed.

Audits verify exact jobs, commands, source revision/diff, build provenance,
binary/matrix hashes, exits, residuals and recomputed medians. Source and
binaries remain fixed throughout timing. These focused results do not
establish full-paper or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/factor-lifecycle-control-trim-*`.
Runner: `validate-factor-lifecycle-control-trim.py`. Audit: `audit-focused-trim.py`.
