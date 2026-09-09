# Remove standalone ownership and ordering diagnostics

Preceding baseline: `2155799`; cumulative baseline: `372f9d1`.

Remove 18 C-source lines: the separator-ownership conflict message and two
remaining predicted-ordering messages. Ownership rejection, ordering
selection, scoring and return behavior remain unchanged. No timers, numeric
thresholds or alignment hints change.

Warning-free Release and ASan builds pass CTest 6/6 each (3.66 and 15.31
seconds). All 579 focused launches pass residual validation at 1e-8.
Incremental paired-median lifecycle changes range from -1.768% to +0.846%;
cumulative changes range from -1.985% to +0.701%. No result reaches the
positive 2.5% review threshold. TSOPF_FS incremental changes are +0.337%,
-0.487%, and +0.399% on CPUs 8–15, CPUs 0–7, and CPU 8.
LeGresley_2508 measures -0.251%.

Sequential pinned H100 comparisons rotate frozen binaries through one
hash-verified execution path without overlapping builds or source edits.
The final audit verifies revision, exact diff, build provenance,
runner/binary/matrix hashes, expected jobs/counts, commands, exits,
residuals and recomputed medians. No new full corpus run was performed;
this focused screen does not prove zero loss or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/standalone-trace-trim-*`;
runner: `validate-standalone-trace-trim.py`; preceding frozen executable:
`standalone-trace-accepted-kls_bench`.

The broader cleanup goal remains active.
