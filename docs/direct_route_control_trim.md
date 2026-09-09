# Remove unused direct-route controls

Baseline: `ad7aa20`. Remove twelve controls unused by tests, scripts and
benchmark drivers:

- `KLS_DISABLE_MODERATE_SINGLE_BLOCK_LEAN_DIRECT_REFACTOR`
- `KLS_DISABLE_SETTLED_LEAN_DIRECT_REFACTOR`
- `KLS_DISABLE_LOW_WORK_BTF_PUBLIC_DIRECT_REFACTOR`
- `KLS_DISABLE_SETTLED_DIRECT_KLU_REFACTOR`
- `KLS_DISABLE_LOW_WORK_SINGLE_BLOCK_DIRECT_REFACTOR`
- `KLS_DISABLE_LOW_WORK_BTF_DIRECT_REFACTOR`
- `KLS_DISABLE_COMPACT_DENSE_SPIKE_FAST_REFACTOR`
- `KLS_DISABLE_COMPACT_MATCH_DIRECT_VALUES`
- `KLS_DISABLE_DEFERRED_LEAN_VALUE_PREP`
- `KLS_DISABLE_SYMMETRIC_PARTIAL_DIAGONAL_DIRECT_VALUES`
- `KLS_DISABLE_LOW_WORK_BTF_DIRECT_VALUES`
- `KLS_DISABLE_COMPACT_DENSE_SPIKE_DIRECT_VALUES`

Keep the default direct-value and direct-refactor predicates, including
map-readiness, scaling, lifecycle, numerical-certificate and fallback checks.
No kernel or structural threshold changes. The removed environment overrides
are no longer supported; public C API unchanged. Update the historical
direct-value isolation note after audits. Net reduction: 13 source lines.

## Validation

Release and ASan builds pass without reported warnings. Release CTest passes
6/6 in 4.13 s; ASan passes 6/6 in 17.15 s. `git diff --check` passes.
Executable text decreases by 896 bytes; data and BSS are unchanged.

All 579 focused H100 launches pass the provenance/result audit. Paired
lifecycle changes range from -2.384% to +2.024% versus accepted and from
-3.496% to +0.370% versus cumulative baseline `372f9d1`.

Repeat TSOPF_RS_b9_c6 / t8_32 because its +2.024% result was near the 2.5%
review threshold. All 72 additional launches pass the audit; the 24-triplet
confirmation measures +0.236% versus accepted and +0.370% versus cumulative.
The slowdown did not repeat. No alignment change is needed.

Audits verify exact jobs, commands, source revision/diff, build provenance,
binary/matrix hashes, exits, residuals and recomputed medians. Source and
binaries remain fixed throughout timing. These focused results do not
establish full-paper or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/direct-route-control-trim-*` and
`direct-route-control-confirm-repeat-*`. Runner:
`validate-direct-route-control-trim.py`. Audits: `audit-focused-trim.py`
and `audit-eligibility-repeat.py`, adapting the latter's expected job/count
to TSOPF_RS_b9_c6 / t8_32 / 24 triplets / 72 launches.
