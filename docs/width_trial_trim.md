# Remove the opt-in EGraph width trial

Baseline: `882c770`. Remove five controls unused by tests, scripts and
benchmark drivers:

- `KLS_ENABLE_EGRAPH_THREAD_TRIAL`
- `KLS_DISABLE_EGRAPH_THREAD_TRIAL`
- `KLS_EGRAPH_WIDTH_MARGIN`
- `KLS_PREMARK_EGRAPH_CLUSTER`
- `KLS_DISABLE_EGRAPH_DOMINATED_QUAD_SKIP`

The thread-width consultation was opt-in only: absent its enable override,
the dispatcher always returned from the full-width branch. Remove the
unreachable-by-default narrow-arm dispatch, width verdict and margin parser,
including its static cache. Keep the caller's full team and the existing
separator-settlement and pipeline-fusion admission rules.

Retain the full-width first-touch warm-up and baseline timing consumed by
the pair/quad trial. Remove the unused width-count array, reduce its sample
and minimum arrays to scalars, and simplify the recorder to the sole live
baseline arm. Keep PTS and factor-reset behavior. Pair/quad selection,
dominated-quad skipping, scratch release and measured completion premarking
remain; only their unused overrides are removed.

Net reduction: 149 source lines across three files. Removed environment
overrides are no longer supported; the public C API is unchanged. The
historical width-margin note is updated after the timing audit.

## Validation

Release and ASan builds pass without reported warnings. Release CTest passes
6/6 in 4.06 s; ASan passes 6/6 in 17.10 s. `git diff --check` passes.
Executable text decreases by 1,088 bytes; data and BSS each decrease by
32 bytes.

All 579 focused H100 launches pass the provenance/result audit. Paired
lifecycle changes range from -5.200% to +0.770% versus accepted and from
-5.520% to +0.621% versus cumulative baseline `372f9d1`. No case reaches
the 2.5% review threshold. No alignment fix is needed. The largest gain,
forced-first-factor 1138_bus, is an observed timing result rather than a
separately established optimization claim.

Audits verify exact jobs, commands, source revision/diff, build provenance,
binary/matrix hashes, exits, residuals and recomputed medians. Source and
binaries remain fixed throughout timing. These focused results do not
establish full-paper or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/width-trial-trim-*`.
Runner: `validate-width-trial-trim.py`. Audit: `audit-focused-trim.py`.
