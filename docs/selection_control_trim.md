# Remove unused refactor-selection controls

Baseline: `99d738d`. Remove eight controls unused by tests, scripts,
benchmark drivers and existing documentation:

- `KLS_DISABLE_GENERIC_SCALED_LEAN_PRESELECTION`
- `KLS_DISABLE_GENERIC_SCALED_LEAN_DECLINE`
- `KLS_DISABLE_LOW_INTENSITY_COLUMN_PRESELECTION`
- `KLS_DISABLE_GENERIC_PACKED_LEAN_PRESELECTION`
- `KLS_DISABLE_BATCH_FLOOR_PROBE`
- `KLS_ENABLE_PADDED_AFTER_LOW_FLOOR`
- `KLS_DISABLE_PADDED_PANEL_PROBE`
- `KLS_DISABLE_GENERIC_EARLY_BATCH_FLOOR_PROBE`

Keep all default selection cost models, dimension/index bounds, numerical
guards, batch-floor/padded-panel trials and warm-up counts. Merge the adjacent
single-block and lean-BTF low-work branches because both publish the same
floor/panel decisions; preserve predicate evaluation order and short-circuiting.
Remove only override behavior. Public C API unchanged; the removed environment
controls are no longer supported. Net reduction: 18 source lines in `src/kls.c`.

## Validation

Release and ASan builds pass without reported warnings. Release CTest passes
6/6 in 4.06 s; ASan passes 6/6 in 17.12 s. `git diff --check` passes.
Executable text decreases by 384 bytes; data and BSS are unchanged.

All 579 focused H100 launches pass the provenance/result audit. Paired
lifecycle changes range from -1.621% to +1.155% versus accepted and from
-3.667% to +0.722% versus cumulative baseline `372f9d1`. No case reaches
the 2.5% review threshold. No alignment fix is needed.

Audits verify exact jobs, commands, source revision/diff, build provenance,
binary/matrix hashes, exits, residuals and recomputed medians. Source and
binaries remain fixed throughout timing. These focused results do not
establish full-paper or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/selection-control-trim-*`.
Runner: `validate-selection-control-trim.py`. Audit: `audit-focused-trim.py`.
