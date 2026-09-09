# Remove unused lean-worker controls

Baseline: `0c3a881`. Remove 24 environment controls unused by tests,
scripts and benchmark drivers:

- `KLS_DISABLE_COMPACT_MATCH_CLEAN_SCRATCH`
- `KLS_DISABLE_COMPACT_MATCH_HOISTED_WORKER`
- `KLS_DISABLE_COMPACT_MATCH_PACKED_INPUT`
- `KLS_DISABLE_COMPACT_MATCH_ROW_FACTOR`
- `KLS_DISABLE_GENERIC_GROUPED_DONE`
- `KLS_DISABLE_GENERIC_GROUPED_RECIP_WORKER`
- `KLS_DISABLE_GENERIC_HOISTED_WORKER`
- `KLS_DISABLE_GENERIC_LEAN_ROW_FACTOR`
- `KLS_DISABLE_GENERIC_PARALLEL_OFFDIAG`
- `KLS_DISABLE_LEAN_CLEAN_SCRATCH`
- `KLS_DISABLE_LEAN_I16_INDICES`
- `KLS_DISABLE_LEAN_I16_PTRS`
- `KLS_DISABLE_LEAN_I16_ROWS`
- `KLS_DISABLE_LEAN_SCALAR_PREFIX`
- `KLS_DISABLE_PARALLEL_REFINE_COPY`
- `KLS_DISABLE_SCALED_GENERIC_GROUPED_DONE`
- `KLS_DISABLE_SCALED_GENERIC_HOISTED_WORKER`
- `KLS_DISABLE_SCALED_LEAN_I16_INDICES`
- `KLS_DISABLE_SCALED_LEAN_U32_PTRS`
- `KLS_DISABLE_SYMMETRIC_SCALAR_FRINGE_CLEAN_SCRATCH`
- `KLS_DISABLE_SYMMETRIC_SCALAR_FRINGE_HOISTED_WORKER`
- `KLS_ENABLE_GENERIC_PARALLEL_OFFDIAG`
- `KLS_ENABLE_LEAN_I16_INDICES`
- `KLS_ENABLE_SCALED_LEAN_ROW_FACTOR`

Keep the default scratch-cleanliness guard, scaled and unscaled workers,
compact-index range checks and allocation fallbacks, row-value capability
checks, scalar-prefix handling, off-diagonal profitability checks and
parallel refinement-copy thresholds. Simplify tautological scaled/unscaled
admission checks and remove override-only force/disable branches. Retain
the deferred-value-preparation exclusion for grouped completion: those
producer phases still cannot be composed safely.

These environment overrides are no longer supported; the public C API is
unchanged. Update the historical compact-index experiment note to distinguish
its former diagnostic switch from the retained allocation/range fallback.
Net reduction: 41 source lines across three files.

## Validation

Release and ASan builds pass without reported warnings. Release CTest passes
6/6 in 4.13 s; ASan passes 6/6 in 17.25 s. `git diff --check` passes.
Executable text decreases by 1,672 bytes; data and BSS are unchanged.

All 579 focused H100 launches pass the provenance/result audit. Paired
lifecycle changes range from -1.178% to +1.038% versus accepted and from
-3.392% to +0.606% versus cumulative baseline `372f9d1`. No case reaches
the 2.5% review threshold. No alignment fix is needed.

Audits verify exact jobs, commands, source revision/diff, build provenance,
binary/matrix hashes, exits, residuals and recomputed medians. Source and
binaries remain fixed throughout timing; documentation is updated after
the audit. These focused results do not establish full-paper or
cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/lean-worker-control-trim-*`.
Runner: `validate-lean-worker-control-trim.py`. Audit: `audit-focused-trim.py`.
