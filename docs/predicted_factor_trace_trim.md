# Remove predicted-factor diagnostics

Preceding baseline: `30ef84d`; cumulative baseline: `372f9d1`.

Remove 327 C-source lines from kls_first_factor.inc: 54 flat diagnostic
blocks, three nested output-only dumps/censuses, their phase timers,
and output-only swap/veto counters. Numerical construction, fill gates,
pivot handling, residual acceptance, fallback decisions and live elapsed
accounting remain intact. Alignment hints are unchanged. No tests or
scripts consume KLS_TRACE_PREDICTED output.

Warning-free Release and ASan builds pass CTest 6/6 each (3.59 and 15.40
seconds). All 579 focused launches pass residual validation at 1e-8.
Incremental paired-median lifecycle changes range from -1.096% to +1.265%;
cumulative changes range from -3.564% to +1.706%. No result reaches the
positive 2.5% review threshold. TSOPF_FS changes are +0.109%, -0.840%,
and +0.187% on CPUs 8–15, CPUs 0–7, and CPU 8. LeGresley_2508 is +0.186%.

Sequential pinned H100 comparisons rotate frozen binaries through one
hash-verified execution path, without overlapping builds or source edits.
The final audit verifies revision, exact diff, build provenance,
runner/binary/matrix hashes, expected jobs/counts, commands, exits,
residuals and recomputed medians. This focused screen does not establish
zero loss or cross-machine equivalence; no new full corpus run was made.

Artifacts: `build/prep-trim-repair-dBcoha/predicted-trace-trim-*`;
runner: `validate-predicted-trace-trim.py`; preceding frozen executable:
`predicted-trace-accepted-kls_bench`.

The broader cleanup goal remains active.
