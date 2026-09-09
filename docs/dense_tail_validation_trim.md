# Remove diagnostic dense-tail recomputation

Preceding baseline: `3b4073c`; cumulative baseline: `372f9d1`.

Remove the opt-in KLS_VERIFY_DENSE_TAIL_REFACTOR harness and its caller,
155 C-source lines. The harness read factor values, recomputed into private
temporary arrays and printed comparisons; it did not accept/reject numerical
results or write back factors. Production dense-tail refactoring, fallback,
residual checks and alignment hints remain unchanged. No tests or scripts
reference the removed harness.

Warning-free Release and ASan builds pass CTest 6/6 each (3.64 and 15.38
seconds). All 579 focused launches pass residual validation at 1e-8.
Incremental paired-median lifecycle changes range from -2.503% to +0.836%;
cumulative changes range from -3.903% to +0.771%. No result reaches the
positive 2.5% review threshold. TSOPF_FS incremental changes are +0.050%,
+0.224%, and +0.401% on CPUs 8–15, CPUs 0–7, and CPU 8.
LeGresley_2508 measures -0.358%.

Sequential pinned H100 comparisons rotate frozen binaries through one
hash-verified execution path without overlapping builds or source edits.
The final audit verifies revision, exact diff, build provenance, hashes,
expected jobs/counts, complete commands, exits, residuals and recomputed
medians. No new full corpus run was made; this is not zero-loss or
cross-machine proof.

Artifacts: `build/prep-trim-repair-dBcoha/tailverify-trace-trim-*`;
runner: `validate-tailverify-trace-trim.py`; preceding frozen executable:
`tailverify-trace-accepted-kls_bench`.

The broader cleanup goal remains active.
