# Remove optional dense-help waiting policy

Baseline: `6cefd37`. Remove five environment controls with no references in
tests, benchmark drivers, scripts or existing documentation:

- `KLS_DENSE_HELP`
- `KLS_DENSE_HELP_WAIT_LARGE`
- `KLS_DENSE_HELP_WAIT_LARGE_WIDTH`
- `KLS_DENSE_HELP_WAIT_LARGE_MIN_ENTRIES`
- `KLS_DENSE_HELP_WAIT_LARGE_SECONDS`

Delete their parsing/cache helpers, the opt-in timed owner-lock waiting loop,
and its private atomic waiter counter. Keep normal cooperative assembly and
panel work sharing, structural width/thread gates, the default one-shot lock
attempt, and the existing noncooperative fallback. The wait-large policy was
off by default; work sharing was on. Default lock memory orders, epoch claims,
dependency checks, helper draining and numerical kernels are unchanged.
The removed environment controls are no longer supported; the public C API
is unchanged. Net reduction: 110 source lines.

## Validation

Release and ASan builds succeed without reported warnings. Release CTest
passes 6/6 in 4.12 s; ASan CTest passes 6/6 in 17.24 s. `git diff --check`
passes. The release executable's `size` text column falls from 3,076,255 to
3,074,527 bytes; BSS falls by 32 bytes and data is unchanged.

All 579 focused H100 launches pass the provenance/result audit. Paired
lifecycle changes range from -3.666% to +1.450% versus accepted and from
-3.637% to +1.250% versus cumulative baseline `372f9d1`. No case reaches
the 2.5% review threshold. No alignment fix is needed. Timing gains are not
claimed as algorithmic improvements from this default-preserving cleanup.

The audit checks exact jobs, commands, source revision/diff, build provenance,
binary/matrix hashes, exits, residuals and recomputed medians. Source and
binaries remain fixed during timing. These focused results do not establish
full-paper or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/help-control-trim-*`.
Runner: `validate-help-control-trim.py`. Audit: `audit-focused-trim.py`.
