# Remove MTMETIS graph diagnostics

Preceding baseline: `f17c2c3`; cumulative baseline: `372f9d1`.

Remove 35 C-source lines: the opt-in graph anomaly/isolation census and
file dump under KLS_MT_ND_TRACE/KLS_MT_ND_DUMP. MTMETIS options, ordering
call, locking and result handling remain intact. No tests or scripts use
these diagnostics. KLS_HAVE_MTMETIS=1 is enabled in the Release build;
this is compiled code, not an excluded conditional section.

Warning-free Release and ASan builds pass CTest 6/6 each (3.68 and 15.46
seconds). All 579 focused launches pass residual validation at 1e-8.
Incremental paired-median lifecycle changes range from -1.194% to +1.840%;
cumulative changes range from -2.072% to +0.904%. No result reaches the
positive 2.5% review threshold. TSOPF_FS incremental changes are +1.024%,
+0.042%, and -0.262% on CPUs 8–15, CPUs 0–7, and CPU 8.
LeGresley_2508 measures -0.198%. Alignment hints are unchanged.

Sequential pinned H100 comparisons rotate frozen binaries through one
hash-verified execution path, with no overlapping builds or source edits.
The final audit verifies revision, exact diff, build provenance, hashes,
expected jobs/counts, complete commands, exits, residuals and recomputed
medians. No new full corpus run was made; this is not zero-loss or
cross-machine proof.

Artifacts: `build/prep-trim-repair-dBcoha/mttrace-trim-*`;
runner: `validate-mttrace-trim.py`; preceding frozen executable:
`mttrace-accepted-kls_bench`.

The broader cleanup goal remains active.
