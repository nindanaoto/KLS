# Remove output-only condition and residual diagnostics

Preceding baseline: `e1de96a`; cumulative baseline: `372f9d1`.

Remove `KLS_TRACE_GENERIC_UNSCALED_RCOND` and `KLS_TRACE_CONTRACT_RESIDUAL`
print blocks (27 source lines). Condition guards, residual probes,
recovery, timing samples and residual-executor selection remain unchanged.
No test or paper runner consumes these diagnostics. Existing alignment
hints remain unchanged; none is added.

Release and ASan CTest pass 6/6 (3.64 and 15.30 seconds). All 579 focused
launches pass residual validation at 1e-8. Incremental paired-median
lifecycle changes range from -1.843% to +1.101%; cumulative changes range
from -2.141% to +1.236%. No case reaches the positive 2.5% review threshold.
TSOPF_FS incremental results are -0.051%, +0.168%, and +0.098% on CPUs
8–15, CPUs 0–7, and CPU 8. LeGresley_2508 measures -0.286%; forced-first
1138_bus measures -0.114%. These results do not establish zero loss or
cross-machine equivalence; no new full corpus screen was run.

Frozen binaries rotate through one hash-verified execution path with
pinned H100 entrywise-changed refactors and no concurrent builds or source
edits. The final audit verifies revision, exact pending diff, build
provenance, runner/binary/matrix hashes, expected jobs/counts, complete
commands, successful exits, residuals and recomputed medians.

Artifacts: `build/prep-trim-repair-dBcoha/contract-trace-trim-*`; runner:
`validate-contract-trace-trim.py`; preceding frozen executable:
`tiny-debug-accepted-kls_bench`.

The broader low-risk cleanup goal remains active.
