# Remove unused value-preparation and index-cache overrides

Baseline: `de3d0e7`. Remove eight controls with no references in tests,
benchmark drivers, scripts or existing documentation:

- `KLS_DISABLE_GENERIC_FUSED_VALUE_PREP`
- `KLS_ENABLE_FUSED_VALUE_PREP`
- `KLS_ENABLE_GATHER_FUSED_VALUE_PREP`
- `KLS_DISABLE_LEAN_SCALE_CACHE`
- `KLS_ENABLE_REFACTOR_MAP_INDEX32`
- `KLS_ENABLE_REFACTOR_L_INDEX32`
- `KLS_ENABLE_REFACTOR_U_INDEX32`
- `KLS_ENABLE_PTS_DIRECT_USER_VALUES`

Keep the default repeated-workload fused preparation and gather/scatter
fallbacks, lean scale snapshots, structural cache eligibility and 32-bit
width limits. Keep allocation-failure handling and 64-bit paths. Replace
the direct-user-values wrapper with its default PTS-readiness predicate.
The removed environment overrides are no longer supported; the public C
API is unchanged. Net reduction: 45 source lines across five files.

## Validation

Final Release and ASan builds pass. The first build exposed one additional
caller of the removed index-cache helpers in the parallel scheduler; it
was updated to the same default width check before final builds/tests.
Release CTest passes 6/6 in 4.19 s; ASan passes 6/6 in 17.34 s.
`git diff --check` passes. Executable text decreases by 1,280 bytes;
data and BSS are unchanged.

All 579 focused H100 launches pass the provenance/result audit. Paired
lifecycle changes range from -2.060% to +2.355% versus accepted and from
-2.857% to +0.674% versus cumulative baseline `372f9d1`.
LeGresley_2508 at eight threads on CPUs 8-15 was near the 2.5% review
threshold (+2.355%), so it received another 24 rotated triplets. Those
72 launches pass the audit and measure -0.959% versus accepted and
-2.994% versus cumulative baseline. The slowdown did not repeat.
No alignment fix is needed.

Audits verify exact jobs, commands, source revision/diff, build provenance,
binary/matrix hashes, exits, residuals and recomputed medians. Source and
binaries remain fixed throughout timing. These focused results do not
establish full-paper or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/value-index-control-trim-*` and
`value-index-control-confirm-repeat-*`. Runner:
`validate-value-index-control-trim.py`. Audits: `audit-focused-trim.py`
and `audit-eligibility-repeat.py`, with the latter's expected job and count
adapted to LeGresley_2508 / t8_32 / 24 triplets / 72 launches.
