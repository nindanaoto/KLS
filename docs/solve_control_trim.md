# Remove unused solve controls

Baseline: `6a3fb48`. Remove seven environment controls unused by tests,
scripts and benchmark drivers:

- `KLS_GENERIC_CONTRACT_THREADS`
- `KLS_DISABLE_GENERIC_PARALLEL_CONTRACT_RESIDUAL`
- `KLS_DISABLE_GENERIC_PARALLEL_CONTRACT_STATS`
- `KLS_DISABLE_COMPACT_MATCH_FUSED_RHS_PERM`
- `KLS_DISABLE_GENERAL_FUSED_RHS`
- `KLS_DISABLE_GENERAL_FUSED_MATCHED_RHS`
- `KLS_DISABLE_PLAIN_SOLVE_TOURNAMENT`

Keep default worker-count selection, residual preparation and dispatch,
RHS fusion, numerical guards, fallback paths and solve tournaments. Retain
`KLS_FORCE_GENERIC_PARALLEL_CONTRACT_RESIDUAL`, which exercises transformed
residuals in a correctness test. Update the historical disable-control
ablation note after validation. Removed overrides are no longer supported;
the public C API is unchanged. Net reduction: 16 source lines.

## Validation

Release and ASan builds pass without reported warnings. Release CTest passes
6/6 in 4.06 s; ASan passes 6/6 in 17.15 s. `git diff --check` passes.
Executable text increases by 2,304 bytes despite the source reduction;
data and BSS are unchanged.

All 579 focused H100 launches pass the provenance/result audit. Paired
lifecycle changes range from -1.501% to +1.536% versus accepted and from
-2.834% to +2.399% versus cumulative baseline `372f9d1`.

Repeat twotone / t8_32 because its cumulative result approaches the 2.5%
review threshold. All 27 additional launches pass the audit. The nine-triplet
confirmation measures +0.694% versus accepted and +1.047% versus cumulative;
the near-threshold slowdown does not repeat. No alignment change is needed.

Audits verify exact jobs, commands, source revision/diff, build provenance,
binary/matrix hashes, exits, residuals and recomputed medians. Source and
binaries remain fixed throughout timing. These focused results do not
establish full-paper or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/solve-control-trim-*` and
`solve-control-confirm-repeat-*`. Runner: `validate-solve-control-trim.py`.
Audits: `audit-focused-trim.py` and `audit-eligibility-repeat.py`, adapting
the latter's expected jobs/count to twotone / t8_32 / nine triplets /
27 launches.
