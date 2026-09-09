# Remove blocked-update diagnostic controls

Baseline: `b65d584`. Remove the unused environment overrides
`KLS_ENABLE_BLOCKED_TRAILING_UPDATE` and `KLS_BLOCKED_TRAILING_MIN_WORK`,
their parsing helpers, the private minimum-work field, and its work estimate.
The default target-major kernel, structural row/column gates, allocation
fallback, and dependency-major fallback remain. The setup/teardown enable
state remains because it is not constant throughout the shared object's life.
Net reduction: 38 source lines. These environment overrides are no longer
supported; no public C API changes.

## Validation

Release CTest passes 6/6 in 4.14 s and ASan CTest passes 6/6 in 17.27 s.
Builds succeed without reported warnings; `git diff --check` passes.
All 579 focused H100 launches pass the provenance/result audit. Paired
lifecycle changes range from -2.697% to +1.809% versus accepted and
-3.004% to +0.375% versus cumulative baseline `372f9d1`.

The largest incremental slowdown, TSOPF_RS_b9_c6 on eight threads/32 MB
placement, receives another 24 rotated triplets (72 launches): -0.282%
versus accepted and +0.478% versus cumulative. All 651 launches are valid;
no persistent notable regression is observed and no alignment fix is used.

Audits verify source revision/diff, build provenance, jobs, commands,
binary/matrix hashes, exits, residuals and recomputed medians. Source and
binaries stay fixed during timing. These focused tests do not establish
full-paper or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/blocked-control-trim-*` and
`blocked-control-confirm-repeat-*`. Runner: `validate-blocked-control-trim.py`.
Audits: `audit-focused-trim.py` and `audit-eligibility-repeat.py` adapted to
the single 24-triplet TSOPF_RS job.
