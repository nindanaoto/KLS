# Retire live row-pipeline tracing

Preceding baseline: `516c256`; cumulative baseline: `372f9d1`.

Remove 39 C-source lines: the trace flag, interval, printer and call sites.
Worker synchronization, progress state, pivot-tail recovery and numeric
accounting remain intact. Alignment hints are unchanged.

Remove the benchmark runner's trace-only failure reruns/options and its
unused imports. Analyze-only diagnostics remain the default; none disables
them. Update README accordingly. The standalone summarizer remains available
for archived logs. Runner --help and mocked none/analyze dispatch checks pass.

Warning-free Release and ASan builds pass CTest 6/6 each (3.61 and 15.34
seconds). All 579 focused launches pass residual validation at 1e-8.
Incremental paired-median lifecycle changes range from -1.028% to +2.814%;
cumulative changes range from -2.246% to +1.514%. TSOPF_FS incremental changes
are -0.385%, -0.292%, and +0.002% on CPUs 8–15, CPUs 0–7, and CPU 8.

1138_bus with forced first factor initially flags +2.814% incremental.
An independent 24-triplet repeat with identical frozen binaries passes all
72 launches and measures +1.839% incremental, +0.660% cumulative. The
initial magnitude does not reproduce above the 2.5% review threshold;
a smaller positive difference remains, not proof of exact parity. No
alignment intervention is retained. Total valid launches: 651.

Sequential pinned H100 comparisons rotate frozen binaries through one
hash-verified execution path, without overlapping builds or source edits.
Both audits verify revision, exact diff, build provenance, hashes, expected
jobs/counts, complete commands, exits, residuals and recomputed medians.
No new full corpus run was made; this is not cross-machine or zero-loss proof.

Artifacts: `build/prep-trim-repair-dBcoha/pipeline-trace-trim-*`, including
repeat artifacts; runner: `validate-pipeline-trace-trim.py`; repeat audit:
`audit-pipeline-repeat.py`; baseline: `pipeline-trace-accepted-kls_bench`.

The broader cleanup goal remains active.
