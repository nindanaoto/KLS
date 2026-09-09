# Remove unused first-factor experiments

Baseline: `b230be1`. Remove two disabled-by-default experimental paths with
no test, script or benchmark-driver callers:

- `KLS_ENABLE_ROW_PIPELINE_SUPERNODE_PRODUCER_BATCH`
- `KLS_ENABLE_ROW_FIRST_STRICT_SEPARATOR_PIVOT_SCOPE`

Delete the supernode producer's candidate and execution helpers, its enabled
and active fields, call site and associated synchronization branches. Keep
the default scalar producer, including its live external-update guards and
condition waits, and keep consumer-side supernode arithmetic. Delete strict
separator pivot scoping and its private context field; retain default
retained-component-extent pivot selection and exact-component reporting.
Remove an empty block adjacent to initialization. Net reduction: 226 source
lines, including two helpers and three private fields. Removed overrides are
no longer supported; the public C API is unchanged.

## Validation

Release and ASan builds pass without reported warnings. Release CTest passes
6/6 in 4.12 s; ASan passes 6/6 in 17.30 s. `git diff --check` passes.
Executable text decreases by 3,424 bytes; data and BSS are unchanged.

All 579 focused H100 launches pass the provenance/result audit. Paired
lifecycle changes range from -1.885% to +1.441% versus accepted and from
-3.909% to +1.779% versus cumulative baseline `372f9d1`. No case crosses
the 2.5% regression review threshold. No alignment change is needed.

Audits verify exact jobs, commands, source revision/diff, build provenance,
binary/matrix hashes, exits, residuals and recomputed medians. Source and
binaries remain fixed throughout timing. These focused results do not
establish full-paper or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/producer-experiment-trim-*`.
Runner: `validate-producer-experiment-trim.py`. Audit: `audit-focused-trim.py`.
