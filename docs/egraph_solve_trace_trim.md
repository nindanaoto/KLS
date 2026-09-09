# Remove remaining egraph solve diagnostics

Preceding baseline: `1f6d529`; cumulative baseline: `372f9d1`.

Remove 34 source lines: the row-solve phase trace flag, timestamps, macro
and calls, plus the plain-single-unscaled dispatch log. Solve ordering,
numerical work, live measurements and synchronization remain intact.
Existing alignment hints are unchanged. No test or paper runner consumes
the removed output.

Warning-free Release and ASan builds pass CTest 6/6 (3.65 and 15.38 seconds).
All 579 focused launches pass residual validation at 1e-8. Incremental
paired-median lifecycle changes range from -1.626% to +1.100%; cumulative
changes range from -5.741% to +1.849%. No case reaches the positive 2.5%
review threshold. TSOPF_FS incremental results are -0.674%, +0.532%, and
-0.034% on CPUs 8–15, CPUs 0–7, and CPU 8; LeGresley_2508 measures -0.782%.
These results do not establish zero loss or cross-machine equivalence;
no new full corpus screen was run.

Frozen binaries rotate through one hash-verified execution path with
pinned H100 entrywise-changed refactors and no concurrent builds or source
edits. The final audit verifies revision, exact pending diff, build
provenance, runner/binary/matrix hashes, expected jobs/counts, complete
commands, successful exits, residuals and recomputed medians.

Artifacts: `build/prep-trim-repair-dBcoha/egraph-trace-trim-*`; runner:
`validate-egraph-trace-trim.py`; preceding frozen executable:
`repair-trace-accepted-kls_bench`.

The broader low-risk cleanup goal remains active.
