# Remove output-only solve dispatch diagnostics

Preceding baseline: `a52318d`; cumulative baseline: `372f9d1`.

Remove `KLS_TRACE_SOLVE_PATH`, `KLS_TRACE_TINY_SINGLETON_SOLVE`, and
`KLS_TRACE_REFINE_SUPPORT` print blocks and the local solve-path trace flag
(53 source lines). Dispatch, tiny-singleton preparation, live sparse RHS
support counting, refinement, acceptance timing and numerical work remain
unchanged. No test or paper runner consumes the removed output. Existing
alignment hints remain unchanged; none is added.

Release and ASan CTest pass 6/6 (3.60 and 15.32 seconds). All 579 focused
launches pass residual validation at 1e-8. Incremental paired-median
lifecycle changes range from -1.520% to +2.205%; cumulative changes range
from -3.475% to +1.122%. No case reaches the positive 2.5% review threshold.
The largest incremental increase is ASIC_680ks on CPUs 0–7 (+2.205%,
twelve pairs), not independently repeated. TSOPF_FS incremental results
are -0.229%, +0.217%, and -0.311% on CPUs 8–15, CPUs 0–7, and CPU 8.
LeGresley_2508 measures +0.425%; forced-first 1138_bus measures -0.798%.
These results do not establish zero loss or cross-machine equivalence;
no new full corpus screen was run.

Frozen binaries rotate through one hash-verified execution path with
pinned H100 entrywise-changed refactors and no concurrent builds or source
edits. The final audit verifies revision, exact pending diff, build
provenance, runner/binary/matrix hashes, expected jobs/counts, complete
commands, successful exits, residuals and recomputed medians.

Artifacts: `build/prep-trim-repair-dBcoha/solve-trace-trim-*`; runner:
`validate-solve-trace-trim.py`; preceding frozen executable:
`refine-trace-accepted-kls_bench`.

The broader low-risk cleanup goal remains active.
