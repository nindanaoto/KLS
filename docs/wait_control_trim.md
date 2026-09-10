# Remove unused worker-wait overrides

Baseline: `67b1deb`. Remove three controls unused by tests, scripts and
benchmark drivers:

- `KLS_DISABLE_REFACTOR_POOL_BUSY_WAIT`
- `KLS_DISABLE_EGRAPH_BUSY_WAIT`
- `KLS_EGRAPH_WORKER_SPIN_OVERRIDE`

Keep default pool wait mode, spin durations, adaptive spin increases,
condition-variable fallback and shutdown synchronization. Delete the spin
parser and override checks; initialize the two wait-mode fields to their
existing default of one. The fixed fields and their branches remain for a
separate cleanup (subsequently completed in `fixed_wait_mode_trim.md`).
No default threshold or numerical kernel changes.
Net reduction: 19 source lines. Removed overrides are no longer supported;
the public C API is unchanged. Earlier diagnostic ablations using these
names describe historical executables, not current runtime options.

## Validation

Release and ASan builds pass without reported warnings. Release CTest passes
6/6 in 4.11 s; ASan passes 6/6 in 17.26 s. `git diff --check` passes.
Executable text increases by 128 bytes; data and BSS are unchanged.

All 579 focused H100 launches pass the provenance/result audit. Paired
lifecycle changes range from -1.511% to +0.810% versus accepted and from
-3.302% to +2.202% versus cumulative baseline `372f9d1`.

Repeat first-factor 1138_bus / t8_32 because its cumulative result approaches
the 2.5% review threshold. All 72 additional launches pass the audit. The
24-triplet confirmation measures +0.195% versus accepted and -1.861% versus
cumulative; the near-threshold slowdown does not repeat. No alignment change
is needed.

Audits verify exact jobs, commands, source revision/diff, build provenance,
binary/matrix hashes, exits, residuals and recomputed medians. Source and
binaries remain fixed throughout timing. These focused results do not
establish full-paper or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/wait-control-trim-*` and
`wait-control-confirm-first-*`. Runner: `validate-wait-control-trim.py`.
Audits: `audit-focused-trim.py` and `audit-eligibility-repeat.py`, adapting
the latter's expected group/jobs/count to first / 1138_bus / t8_32 /
24 triplets / 72 launches.
