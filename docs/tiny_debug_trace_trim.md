# Remove small-system debug dumps

Preceding baseline: `b8e5460`; cumulative baseline: `372f9d1`.

Remove `KLS_TRACE_X` factor and solve dumps and the local trace flag
(67 source lines). Diagnostic loops only read matrix/factor/vector data;
no numerical operations, policies or recovery checks are removed.
No test or paper runner consumes this output. Existing alignment hints
remain unchanged; none is added.

Release and ASan CTest pass 6/6 (3.68 and 15.43 seconds). All 579 focused
launches pass residual validation at 1e-8. Incremental paired-median
lifecycle changes range from -0.802% to +1.352%; cumulative changes range
from -3.142% to +0.706%. No case reaches the positive 2.5% review threshold.
TSOPF_FS incremental results are +0.967%, -0.493%, and -0.209% on CPUs
8–15, CPUs 0–7, and CPU 8. LeGresley_2508 measures +0.779%; forced-first
1138_bus measures -0.519%. These results do not establish zero loss or
cross-machine equivalence; no new full corpus screen was run.

Frozen binaries rotate through one hash-verified execution path with
pinned H100 entrywise-changed refactors and no concurrent builds or source
edits. The final audit verifies revision, exact pending diff, build
provenance, runner/binary/matrix hashes, expected jobs/counts, complete
commands, successful exits, residuals and recomputed medians.

Artifacts: `build/prep-trim-repair-dBcoha/tiny-debug-trim-*`; runner:
`validate-tiny-debug-trim.py`; preceding frozen executable:
`solve-trace-accepted-kls_bench`.

The broader low-risk cleanup goal remains active.
