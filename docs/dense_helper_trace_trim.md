# Remove dense-helper and compact-map32 diagnostics

Preceding baseline: `5ef83c4`; cumulative baseline: `372f9d1`.

Remove 75 source lines: four diagnostic atomic counters and their updates,
the dense-helper trace flag and owner-row counter, dense-helper output and
compact-map32 decision logs. Synchronization atomics, failure handling,
numerical work, timing and selection policies remain intact. No test or
paper runner consumes the removed output. Existing alignment hints are
unchanged.

Warning-free Release and ASan builds pass CTest 6/6 (3.60 and 15.38 seconds).
All 579 focused launches pass residual validation at 1e-8. Incremental
paired-median lifecycle changes range from -2.311% to +1.449%; cumulative
changes range from -2.442% to +0.935%. No case reaches the positive 2.5%
review threshold. TSOPF_FS incremental results are +0.174%, +1.449%, and
-0.020% on CPUs 8–15, CPUs 0–7, and CPU 8; LeGresley_2508 measures +0.499%.
These results do not establish zero loss or cross-machine equivalence;
no new full corpus screen was run.

Frozen binaries rotate through one hash-verified execution path with
pinned H100 entrywise-changed refactors and no concurrent builds or source
edits. The final audit verifies revision, exact pending diff, build
provenance, runner/binary/matrix hashes, expected jobs/counts, complete
commands, successful exits, residuals and recomputed medians.

Artifacts: `build/prep-trim-repair-dBcoha/dense-trace-trim-*`; runner:
`validate-dense-trace-trim.py`; preceding frozen executable:
`worker-trace-accepted-kls_bench`.

The broader low-risk cleanup goal remains active.
