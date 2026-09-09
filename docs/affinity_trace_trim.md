# Remove output-only lean-affinity diagnostics

Preceding baseline: `283b4bc`; cumulative baseline: `372f9d1`.

Remove the two `KLS_TRACE_LEAN_AFFINITY` print blocks (17 source lines).
Scheduling, work estimates, projected savings and acceptance decisions
remain unchanged. No test or paper runner consumes this output. Existing
alignment hints remain unchanged; none is added.

Release and ASan CTest pass 6/6 (3.66 and 15.40 seconds). All 579 focused
launches pass residual validation at 1e-8. Incremental paired-median
lifecycle changes range from -2.386% to +0.901%; cumulative changes range
from -2.876% to +1.007%. No case reaches the positive 2.5% review threshold.
TSOPF_FS incremental results are +0.901%, -0.456%, and -0.146% on CPUs
8–15, CPUs 0–7, and CPU 8. LeGresley_2508 measures -0.421%; forced-first
1138_bus measures -2.386%. These results do not establish zero loss or
cross-machine equivalence; no new full corpus screen was run.

Frozen binaries rotate through one hash-verified execution path with
pinned H100 entrywise-changed refactors and no concurrent builds or source
edits. The final audit verifies revision, exact pending diff, build
provenance, runner/binary/matrix hashes, expected jobs/counts, complete
commands, successful exits, residuals and recomputed medians.

Artifacts: `build/prep-trim-repair-dBcoha/affinity-trace-trim-*`; runner:
`validate-affinity-trace-trim.py`; preceding frozen executable:
`preselection-trace-accepted-kls_bench`.

The broader low-risk cleanup goal remains active.
