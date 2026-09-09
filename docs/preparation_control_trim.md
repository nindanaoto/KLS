# Remove unused pattern-preparation overrides

Baseline: `5ab8418`. Remove ten overrides with no references in tests,
benchmark drivers, scripts or existing documentation:

- `KLS_DEFER_ROW_TRANSPOSE_PLANS`
- `KLS_DEFER_ROW_SUCCESSORS`
- `KLS_DEFER_ROW_INPUT_CLEANUP`
- `KLS_DISABLE_SCALAR_ROW_BATCH_GROUPS`
- `KLS_SKIP_SUCCESSOR_SORT`
- `KLS_DISABLE_DIRECT_NUMERIC_LEAN_PATTERN`
- `KLS_ENABLE_PARALLEL_ROW_PATTERN`
- `KLS_ENABLE_DIRECT_NUMERIC_ROW_PATTERN`
- `KLS_DISABLE_GENERIC_PARALLEL_LEAN_PATTERN`
- `KLS_ENABLE_GENERIC_PARALLEL_LEAN_PATTERN`

Remove the corresponding deferral/skip branches, the override-only generic
full-pattern parallel call, and forced admission to the parallel builders.
Retain default successor construction and validation, cleanup analysis,
same-level batching, successor sorting, transpose-plan construction and
serial-only plan upgrading. Keep the default parallel builders, structural
admission thresholds, and allocation/numerical fallbacks. The removed
overrides are no longer supported; the public C API is unchanged.
Net reduction: 47 source lines across two files.

## Validation

Both builds succeed without reported warnings. Release CTest passes 6/6
in 4.11 s; ASan CTest passes 6/6 in 17.18 s. `git diff --check` passes.
Executable `size` text decreases by 5,328 bytes; data/BSS are unchanged.

All 579 focused H100 launches pass the provenance/result audit. Paired
lifecycle changes range from -2.496% to +0.498% versus accepted and from
-3.623% to +0.618% versus cumulative baseline `372f9d1`. No case reaches
the 2.5% review threshold. No alignment fix is needed.

Audits verify exact jobs, commands, source revision/diff, build provenance,
binary/matrix hashes, exits, residuals and recomputed medians. Source and
binaries remain fixed during timing. Focused results do not establish
full-paper or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/preparation-control-trim-*`.
Runner: `validate-preparation-control-trim.py`. Audit: `audit-focused-trim.py`.
