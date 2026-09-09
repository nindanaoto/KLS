# Remove matching and scaling diagnostics

Preceding baseline: `ba7df5c`; cumulative baseline: `372f9d1`.

Remove 100 net source lines of matching phase/swap logs, scaling and
pivot-pressure logs, pipe-route output, matched-tolerance output and AMF
etree diagnostics. Remove their diagnostic clocks, saved greedy counts and
one now-unused local dimension. Etree levels remain measured whenever the
selection policy requests them; only the diagnostic override is removed.
Matching algorithms, live policy statistics, acceptance checks and timing
remain intact. Existing alignment hints are unchanged.

Final builds are warning-free. Release and ASan CTest pass 6/6 (3.61 and
15.36 seconds). All 579 focused launches pass residual validation at 1e-8.
Incremental paired-median lifecycle changes range from -1.441% to +2.319%;
cumulative changes range from -2.153% to +1.418%. No case reaches the
positive 2.5% review threshold. Forced-first 1138_bus has the largest
incremental increase (eight triplets), while measuring -1.141% cumulatively.
TSOPF_FS incremental changes are -0.240%, -0.987%, and +0.286% on CPUs
8–15, CPUs 0–7, and CPU 8; LeGresley_2508 measures -1.330%. These results
do not establish zero loss or cross-machine equivalence; no new full corpus
screen was run.

Frozen binaries rotate through one hash-verified execution path with
pinned H100 entrywise-changed refactors and no concurrent builds or source
edits. The final audit verifies revision, exact pending diff, build
provenance, runner/binary/matrix hashes, expected jobs/counts, complete
commands, successful exits, residuals and recomputed medians.

Artifacts: `build/prep-trim-repair-dBcoha/match-trace-trim-*`; runner:
`validate-match-trace-trim.py`; preceding frozen executable:
`analyze-trace-accepted-kls_bench`.

The broader low-risk cleanup goal remains active.
