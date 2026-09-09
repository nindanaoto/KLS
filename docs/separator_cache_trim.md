# Remove opt-in separator-plan caching

Baseline: `ff6f217`. Remove the optional separator-plan cache as one unit:
its save, restore and invalidation helpers; sixteen private solver fields;
allocation/reset bookkeeping; cached ordering checks; ownership transfer;
and cache invalidation after component calibration. Keep component
calibration itself, the normal schedule builder, private-order validation,
dependency-ready fallback and local schedule cleanup.

The cache was disabled by default. Its two environment controls,
`KLS_ENABLE_ROW_SEPARATOR_PLAN_CACHE` and
`KLS_ROW_REFACTOR_CACHE_RESET_PIPELINE_ONLY`, had no references in tests,
benchmark drivers, scripts or existing documentation and are no longer
supported. There is no public C API change. Net reduction: 365 source lines.

## Validation

Both builds succeed without reported warnings. Release CTest passes 6/6
in 4.06 s; ASan CTest passes 6/6 in 17.09 s. `git diff --check` passes.
The release executable's `size` text column falls by 3,008 bytes; data/BSS
are unchanged. Searches confirm no cache fields or helper references remain
in the source.

All 579 focused H100 launches pass the provenance/result audit. Paired
lifecycle changes range from -1.685% to +0.819% versus accepted and from
-2.961% to +1.221% versus cumulative baseline `372f9d1`. No case reaches
the 2.5% review threshold; no alignment fix is needed.

The audit verifies exact jobs, commands, source revision/diff, build
provenance, binary/matrix hashes, exits, residuals and recomputed medians.
Source and binaries stay fixed during timing. These focused results do not
establish full-paper or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/separator-cache-trim-*`.
Runner: `validate-separator-cache-trim.py`. Audit: `audit-focused-trim.py`.
