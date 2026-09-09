# Remove unused active-catchup experiment

Baseline: `d59b6f8`. Remove three controls unused by tests, scripts and
benchmark drivers:

- `KLS_ENABLE_ROW_PIPELINE_ACTIVE_CATCHUP_BATCH`
- `KLS_DISABLE_ROW_PIPELINE_PRODUCER_BATCH`
- `KLS_DISABLE_ROW_PIPELINE_PIVOT_SUPERNODE_REBASE`

Delete disabled-by-default active catch-up and its candidate retry. Remove
the now-unused candidate output flag and bounded stop-dependency parameter,
including three stopping/clipping branches. Fold the unbounded ready-row
wrapper into its implementation and delete the redundant declaration.
Delete three private enable fields. Keep default scalar producer batching,
its synchronization, and post-pivot supernode rebasing with all fallback
checks. Net reduction: 124 source lines. Removed overrides are no longer
supported; the public C API is unchanged.

## Validation

Release and ASan builds pass without reported warnings. Release CTest passes
6/6 in 4.05 s; ASan passes 6/6 in 17.21 s. `git diff --check` passes.
Executable text decreases by 1,552 bytes; data and BSS are unchanged.

All 579 focused H100 launches pass the provenance/result audit. Paired
lifecycle changes range from -2.715% to +1.113% versus accepted and from
-3.483% to +0.986% versus cumulative baseline `372f9d1`. No case crosses
the 2.5% regression review threshold. No alignment change is needed.

Audits verify exact jobs, commands, source revision/diff, build provenance,
binary/matrix hashes, exits, residuals and recomputed medians. Source and
binaries remain fixed throughout timing. These focused results do not
establish full-paper or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/catchup-experiment-trim-*`.
Runner: `validate-catchup-experiment-trim.py`. Audit: `audit-focused-trim.py`.
