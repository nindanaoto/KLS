# Remove remaining factor-phase logging

Preceding baseline: `0ddd870`; cumulative baseline: `372f9d1`.

Remove the remaining `KLS_TRACE_FACTOR_PHASES` output in factor entry,
analysis and ordering policy: 388 source lines, including diagnostic-only
timestamps and the ND trial duration field. Two logs shared with
`KLS_TRACE_PREDICTED` are also removed. Live numeric diagnostics, selection
policies, elapsed-time accounting and synchronization remain intact.
No test or paper runner consumes this output. Existing alignment hints
are unchanged; none is added.

Final builds are warning-free. Release and ASan CTest pass 6/6 (3.68 and
15.32 seconds). All 579 focused launches pass residual validation at 1e-8.
Incremental paired-median lifecycle changes range from -1.630% to +2.353%;
cumulative changes range from -1.666% to +1.445%. No case reaches the
positive 2.5% review threshold. The largest incremental increase is forced
first-factor 1138_bus (eight triplets), +0.179% against the cumulative
baseline. TSOPF_FS incremental results are +0.254%, -0.478%, and -0.183%
on CPUs 8–15, CPUs 0–7, and CPU 8; LeGresley_2508 measures +0.534%.
These results do not establish zero loss or cross-machine equivalence;
no new full corpus screen was run.

Frozen binaries rotate through one hash-verified execution path with
pinned H100 entrywise-changed refactors and no concurrent builds or source
edits. The final audit verifies revision, exact pending diff, build
provenance, runner/binary/matrix hashes, expected jobs/counts, complete
commands, successful exits, residuals and recomputed medians.

Artifacts: `build/prep-trim-repair-dBcoha/ordering-trace-trim-*`; runner:
`validate-ordering-trace-trim.py`; preceding frozen executable:
`factor-trace-accepted-kls_bench`.

The broader low-risk cleanup goal remains active.
