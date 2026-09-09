# Remove compact-map and row-solve diagnostic overrides

Baseline: `5982907`. Remove four overrides with no references in tests,
benchmark drivers, scripts or existing documentation:

- `KLS_ENABLE_COMPACT_MATCH_INPUT_COLS16`
- `KLS_ENABLE_GENERIC_COMPACT_MATCH_DIRECT_VALUES`
- `KLS_ROW_SOLVE_CLAIM`
- `KLS_ROW_REFACTOR_PUBLISH`

Delete the opt-in duplicate compact-column builder and diagnostic map/solve
selection branches. Keep the standard mixed-width input builder, direct-map
capability validation, row-solve eligibility checks and automatic publication
decisions. Retain representation and allocation fallbacks. These overrides
are no longer supported; the public C API is unchanged. Net reduction:
38 source lines. Executable `size` text decreases by 680 bytes; data/BSS
are unchanged.

## Validation

Both builds succeed without reported warnings. Release CTest passes 6/6
in 4.11 s and ASan CTest passes 6/6 in 17.21 s. `git diff --check` passes.

All 579 focused H100 launches pass the provenance/result audit. Paired
lifecycle changes range from -1.101% to +1.739% versus accepted and from
-2.429% to +0.919% versus cumulative baseline `372f9d1`. The largest
incremental slowdown, mimo46x46_system on eight threads/96 MB placement,
receives another 24 rotated triplets (72 launches): +0.609% versus accepted
and -1.285% versus cumulative. All 651 launches are valid. No notable
regression is observed and no alignment fix is needed.

Audits verify exact jobs, commands, source revision/diff, build provenance,
binary/matrix hashes, exits, residuals and recomputed medians. Source and
binaries stay fixed during timing. These focused results do not establish
full-paper or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/map-solve-control-trim-*` and
`map-solve-control-confirm-repeat-*`.
Runner: `validate-map-solve-control-trim.py`. Audits: `audit-focused-trim.py`
and `audit-eligibility-repeat.py` adapted to the single mimo repeat job.
