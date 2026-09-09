# Remove output-only compact-cache diagnostics

Preceding baseline: `e9dc19b`; cumulative baseline: `372f9d1`.

Remove `KLS_TRACE_I16_CACHE` and its diagnostic-only local census loops,
plus `KLS_TRACE_I32_SINGLETON_RUN` and `KLS_TRACE_SERIAL_SCALE_PIVOTS`
output (62 source lines). Live identity-prefix metadata, cache construction,
singleton-run admission, and scale-pivot decisions remain unchanged.
No test or paper runner consumes these diagnostics. Existing alignment
hints remain unchanged; no new hint is added.

Release and ASan CTest pass 6/6 (3.65 and 15.32 seconds). All 579 focused
launches pass residual validation at 1e-8. Incremental paired-median
lifecycle changes range from -4.843% to +1.321%; cumulative changes range
from -4.168% to +1.914%. No case reaches the positive 2.5% review threshold.
TSOPF_FS incremental results are +0.643%, -0.479%, and +0.345% on CPUs
8–15, CPUs 0–7, and CPU 8. LeGresley_2508 measures -4.843%; forced-first
1138_bus measures +1.130%. These results do not establish zero loss or
cross-machine equivalence; no new full corpus screen was run.

Frozen binaries rotate through one hash-verified execution path with
pinned H100 entrywise-changed refactors and no concurrent builds or source
edits. The final audit verifies revision, exact pending diff, build
provenance, runner/binary/matrix hashes, expected jobs/counts, complete
commands, successful exits, residuals and recomputed medians.

Artifacts: `build/prep-trim-repair-dBcoha/compact-cache-trim-*`; runner:
`validate-compact-cache-trim.py`; preceding frozen executable:
`i32-prep-accepted-kls_bench`; audit: `audit-focused-trim.py`.

The broader low-risk cleanup goal remains active.
