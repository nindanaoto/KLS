# Remove preparation-consult timing diagnostics

Preceding baseline: `35cc6c7`; cumulative baseline: `372f9d1`.

Remove `KLS_TRACE_PREP_CONSULT`, its local diagnostic timestamp, and the
timing macro/calls (16 source lines). The live preparation timer,
preparation order, thread creation/joins and numerical work remain intact.
No test or paper runner consumes this output. Existing alignment hints
remain unchanged; none is added.

Release and ASan CTest pass 6/6 (3.65 and 15.39 seconds). All 579 focused
launches pass residual validation at 1e-8. Incremental paired-median
lifecycle changes range from -3.097% to +1.643%; cumulative changes range
from -2.556% to +1.351%. No case reaches the positive 2.5% review threshold.
TSOPF_FS incremental results are -0.280%, +0.304%, and +0.237% on CPUs
8–15, CPUs 0–7, and CPU 8. LeGresley_2508 measures +0.753%; forced-first
1138_bus measures -2.668%. These results do not establish zero loss or
cross-machine equivalence; no new full corpus screen was run.

Frozen binaries rotate through one hash-verified execution path with
pinned H100 entrywise-changed refactors and no concurrent builds or source
edits. The final audit verifies revision, exact pending diff, build
provenance, runner/binary/matrix hashes, expected jobs/counts, complete
commands, successful exits, residuals and recomputed medians.

Artifacts: `build/prep-trim-repair-dBcoha/consult-trace-trim-*`; runner:
`validate-consult-trace-trim.py`; preceding frozen executable:
`packed-trace-accepted-kls_bench`.

The broader low-risk cleanup goal remains active.
