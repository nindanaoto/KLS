# Remove worker-selection and affinity diagnostics

Preceding baseline: `148fcf9`; cumulative baseline: `372f9d1`.

Remove 105 source lines of output-only egraph thread-selection and affinity
logging, including the selected-CPU print loop. Thread policies, placement,
timing-based choices, memory release and synchronization remain intact.
No test or paper runner consumes the removed output. Existing alignment
hints are unchanged.

Warning-free Release and ASan builds pass CTest 6/6 (3.63 and 15.46 seconds).
All 579 focused launches pass residual validation at 1e-8. Incremental
paired-median lifecycle changes range from -2.012% to +1.692%; cumulative
changes range from -2.354% to +1.468%. No case reaches the positive 2.5%
review threshold. TSOPF_FS incremental results are +0.987%, -0.877%, and
+0.077% on CPUs 8–15, CPUs 0–7, and CPU 8; LeGresley_2508 measures -1.021%.
These results do not establish zero loss or cross-machine equivalence;
no new full corpus screen was run.

Frozen binaries rotate through one hash-verified execution path with
pinned H100 entrywise-changed refactors and no concurrent builds or source
edits. The final audit verifies revision, exact pending diff, build
provenance, runner/binary/matrix hashes, expected jobs/counts, complete
commands, successful exits, residuals and recomputed medians.

Artifacts: `build/prep-trim-repair-dBcoha/worker-trace-trim-*`; runner:
`validate-worker-trace-trim.py`; preceding frozen executable:
`match-trace-accepted-kls_bench`.

The broader low-risk cleanup goal remains active.
