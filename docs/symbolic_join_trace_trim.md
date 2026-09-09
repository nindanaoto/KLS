# Remove symbolic-join diagnostics

Preceding baseline: `6b3eab0`; cumulative baseline: `372f9d1`.

Remove the final three KLS_TRACE_* output blocks from kls.c (11 lines).
The symbolic wait timer remains live for elapsed-time accounting. Routing,
fallback, numerical work and alignment hints are unchanged.

Warning-free Release and ASan builds pass CTest 6/6 each (3.66 and 15.34
seconds). All 579 focused launches pass residual validation at 1e-8.
Incremental paired-median lifecycle changes range from -2.624% to +1.686%;
cumulative changes range from -2.086% to +3.209%. TSOPF_FS incremental
changes are -0.158%, -0.020%, and +0.021% on CPUs 8–15, CPUs 0–7, and CPU 8.

LeGresley_2508 initially flags +3.209% versus the cumulative baseline
(+1.686% incremental). An independent 24-triplet repeat using unchanged
frozen binaries passes all 72 launches and measures -1.045% cumulative,
+0.308% incremental. The initial slowdown does not reproduce; no alignment
change is justified by these results. Total valid launches: 651.

Sequential pinned H100 comparisons rotate three frozen binaries through
one hash-verified execution path, with no overlapping builds or source
edits. Both final audits verify revision, exact source diff, provenance,
runner/binary/matrix hashes, expected jobs/counts, commands, successful
exits, residuals and recomputed medians. This focused evidence does not
prove zero regression or cross-machine equivalence; no full corpus rerun
was performed for this chunk.

Artifacts: `build/prep-trim-repair-dBcoha/join-trace-trim-*`, including
`join-trace-trim-repeat*`. Runner: `validate-join-trace-trim.py`;
repeat audit: `audit-join-repeat.py`; preceding frozen executable:
`join-trace-accepted-kls_bench`.

The broader cleanup goal remains active.
