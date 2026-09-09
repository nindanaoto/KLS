# Remove unused row-policy controls

Baseline: `f2c6bee`. Remove three controls unused by tests, scripts and
benchmark drivers:

- `KLS_ENABLE_AUTO_ROW_REFACTOR`
- `KLS_CBLAS_PANEL_BLOCK_ROWS`
- `KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_EXEC`

Delete the redundant auto-row environment helper and use the retained,
test-used row-refactor switch at its three call sites. Keep the independent
default high-work admission checks. Delete the CBLAS block-size parser and
cached settings, using the existing compile-time default at its sole caller.
Keep default BTF scalar-run selection and its structural/work guards, removing
only the forced enable/disable override. No numeric kernel or default
threshold changes. Net reduction: 34 source lines and two helpers.
Removed overrides are no longer supported; the public C API is unchanged.

## Validation

Release and ASan builds pass without reported warnings. Release CTest passes
6/6 in 4.11 s; ASan passes 6/6 in 17.22 s. `git diff --check` passes.
Executable text decreases by 264 bytes; data and BSS are unchanged.

All 579 focused H100 launches pass the provenance/result audit. Paired
lifecycle changes range from -3.043% to +0.877% versus accepted and from
-3.476% to +1.897% versus cumulative baseline `372f9d1`. No case crosses
the 2.5% regression review threshold. No alignment change is needed.

Audits verify exact jobs, commands, source revision/diff, build provenance,
binary/matrix hashes, exits, residuals and recomputed medians. Source and
binaries remain fixed throughout timing. These focused results do not
establish full-paper or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/row-policy-control-trim-*`.
Runner: `validate-row-policy-control-trim.py`. Audit: `audit-focused-trim.py`.
