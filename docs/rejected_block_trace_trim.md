# Remove rejected-block repair diagnostics

Preceding baseline: `4647f1d`; cumulative baseline: `372f9d1`.

Remove 52 net source lines of rejected-block stage, eligibility and failure
logging, plus two output-only state variables. Four misleadingly named
`trace_` variables participate in recovery decisions: retain them and rename
them to `repair_`. Recovery calls, failure handling, cleanup and live timing
remain intact. Existing alignment hints are unchanged.

Warning-free Release and ASan builds pass CTest 6/6 (3.67 and 15.35 seconds).
All 579 focused launches pass residual validation at 1e-8. Incremental
paired-median lifecycle changes range from -1.758% to +0.278%; cumulative
changes range from -2.316% to +0.614%. No case reaches the positive 2.5%
review threshold. TSOPF_FS incremental changes are +0.123%, +0.065%, and
-0.007% on CPUs 8–15, CPUs 0–7, and CPU 8; LeGresley_2508 measures -0.141%.
These results do not establish zero loss or cross-machine equivalence;
no new full corpus screen was run.

Frozen binaries rotate through one hash-verified execution path with
pinned H100 entrywise-changed refactors and no concurrent builds or source
edits. The final audit verifies revision, exact pending diff, build
provenance, runner/binary/matrix hashes, expected jobs/counts, complete
commands, successful exits, residuals and recomputed medians.

Artifacts: `build/prep-trim-repair-dBcoha/repair-trace-trim-*`; runner:
`validate-repair-trace-trim.py`; preceding frozen executable:
`preparation-trace-accepted-kls_bench`.

The broader low-risk cleanup goal remains active.
