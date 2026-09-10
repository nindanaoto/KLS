# Remove unused compact-policy controls

Baseline: `027d942`. Remove four controls unused by tests, scripts and
benchmark drivers:

- `KLS_DISABLE_COMPACT_AMF_TWO_BLOCK_POLICY`
- `KLS_DISABLE_COMPACT_AMF_TWO_BLOCK_SPECIALIZED_WORKER`
- `KLS_COMPACT_AMF_TWO_BLOCK_SCHEDULE_WEIGHTS`
- `KLS_DISABLE_COMPACT_AMF_TWO_BLOCK_PARALLEL_RESIDUAL`

Delete three environment helpers and the four-value scheduling-weight parser.
Keep default scheduling weights, symbolic and representation guards, packed
worker admission, and parallel residual computation. Make the unchanged
weights const. No default threshold or numerical kernel changes.
Net reduction: 35 source lines. Removed overrides are no longer supported;
the public C API is unchanged. Prior diagnostic ablations using these names
describe historical executables, not current runtime options.

## Validation

Release and ASan builds pass without reported warnings. Release CTest passes
6/6 in 4.13 s; ASan passes 6/6 in 17.24 s. `git diff --check` passes.
Executable text decreases by 408 bytes; data and BSS are unchanged.

All 579 focused H100 launches pass the provenance/result audit. Paired
lifecycle changes range from -2.571% to +2.261% versus accepted and from
-3.344% to +1.923% versus cumulative baseline `372f9d1`.

Repeat first-factor bcircuit / t8_32 because its +2.261% result approaches
the 2.5% review threshold. All 72 additional launches pass the audit. The
24-triplet confirmation measures -0.026% versus accepted and -0.535% versus
cumulative; the slowdown does not repeat. No alignment change is needed.

Audits verify exact jobs, commands, source revision/diff, build provenance,
binary/matrix hashes, exits, residuals and recomputed medians. Source and
binaries remain fixed throughout timing. These focused results do not
establish full-paper or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/compact-policy-control-trim-*` and
`compact-policy-control-confirm-first-*`. Runner:
`validate-compact-policy-control-trim.py`. Audits: `audit-focused-trim.py`
and `audit-eligibility-repeat.py`, adapting the latter's expected group,
jobs/count to first / bcircuit / t8_32 / 24 triplets / 72 launches.
