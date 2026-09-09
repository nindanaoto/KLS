# Remove unused solve-recovery controls

Baseline: `6b4c5ed`. Remove thirteen controls unused by tests, scripts and
benchmark drivers:

- `KLS_DISABLE_GENERIC_UNSCALED_RCOND_GUARD`
- `KLS_DISABLE_GMRES_SOLVE_RECOVERY`
- `KLS_DISABLE_LSQR_SOLVE_RECOVERY`
- `KLS_DISABLE_RCOND_SOLVE_CONTRACT`
- `KLS_DISABLE_SOLVE_CONTRACT_PROBE`
- `KLS_DISABLE_ORDINARY_SELF_CHECK_L2_CONTRACT`
- `KLS_DISABLE_BOUNDED_PUBLIC_RAW_L2`
- `KLS_DISABLE_CERTIFIED_COMPACT_ROW_CONTRACT_SETTLE`
- `KLS_DISABLE_SETTLED_PTS_RAW_L2_CERTIFICATE`
- `KLS_DISABLE_SPARSE_REFINEMENT_RHS`
- `KLS_DISABLE_SPARSE_REFINEMENT_NORM_REUSE`
- `KLS_SPARSE_REFINEMENT_DROP_RELATIVE`
- `KLS_DISABLE_CONTRACT_SINGLE_SHOT`

Keep all default numerical algorithms, recovery checks, structural guards,
residual certificates and tolerances. Replace the drop-tolerance parser with
its default 1e-9 value at the existing multiplication. Retain the test-used
`KLS_DISABLE_PROMOTED_TOLERANCE_L2_RECOVERY` control. No recovery kernel is
removed. Net reduction: 22 source lines. Removed overrides are no longer
supported; the public C API is unchanged. Earlier diagnostic ablations using
these names describe historical executables, not current runtime options.

## Validation

Release and ASan builds pass without reported warnings. Release CTest passes
6/6 in 4.06 s; ASan passes 6/6 in 17.34 s. `git diff --check` passes.
Executable text decreases by 624 bytes; data and BSS are unchanged.

All 579 focused H100 launches pass the provenance/result audit. Paired
lifecycle changes range from -1.848% to +1.492% versus accepted and from
-4.343% to +1.227% versus cumulative baseline `372f9d1`. No case crosses
the 2.5% regression review threshold. No alignment change is needed.

Audits verify exact jobs, commands, source revision/diff, build provenance,
binary/matrix hashes, exits, residuals and recomputed medians. Source and
binaries remain fixed throughout timing. These focused results do not
establish full-paper or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/recovery-control-trim-*`.
Runner: `validate-recovery-control-trim.py`. Audit: `audit-focused-trim.py`.
