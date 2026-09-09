# Remove first-factor decision logging

Preceding baseline: `70128ff`; cumulative baseline: `372f9d1`.

Remove 160 C-source lines: four first-factor logging blocks, the dominant
BTF logging helper and call, their cached environment flag, the now-unused
failure-name formatter, an output-only initial-plan snapshot, and failed
position/global-row fields and writes. Failed-row and failure-reason state
remain live for partial-commit safety. Ownership checks, numerical calls,
fallbacks and alignment hints are unchanged.

The benchmark failure-diagnostic runner no longer sets or reports the removed
KLS_TRACE_KLS_FIRST_FACTOR flag (two Python lines). Its separate row-pipeline
trace remains enabled. The runner's --help check passes.

Warning-free Release and ASan builds pass CTest 6/6 each (3.58 and 15.41
seconds). All 579 focused launches pass residual validation at 1e-8.
Incremental paired-median lifecycle changes range from -1.268% to +1.013%;
cumulative changes range from -2.870% to +1.983%. No result reaches the
positive 2.5% review threshold. TSOPF_FS changes are -0.310%, +0.316%,
and +0.279% on CPUs 8–15, CPUs 0–7, and CPU 8. LeGresley_2508 is -0.627%.

Sequential pinned H100 comparisons rotate frozen binaries through one
hash-verified execution path; no builds or source edits overlap timing.
The final audit verifies revision, exact diff, build provenance,
runner/binary/matrix hashes, expected jobs/counts, commands, exits,
residuals and recomputed medians. No new full corpus screen was run;
these focused results do not prove zero loss or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/firstlog-trace-trim-*`;
runner: `validate-firstlog-trace-trim.py`; preceding frozen executable:
`firstlog-trace-accepted-kls_bench`.

The broader cleanup goal remains active.
