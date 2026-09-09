# Remove analysis-stage timing diagnostics

Preceding baseline: `e043bc2`; cumulative baseline: `372f9d1`.

Remove `KLS_TRACE_ANALYZE_STAGES` output, five diagnostic timestamps and
the now-redundant METIS timing wrapper. Rename the inner implementation
to the existing callback name. The result is a net reduction of 42 source
lines. Ordering calls, separator construction, error handling and live
timing remain intact. Existing alignment hints are unchanged.

Final builds are warning-free. Release and ASan CTest pass 6/6 (3.61 and
15.44 seconds). All 579 focused launches pass residual validation at 1e-8.
Incremental paired-median lifecycle changes range from -1.631% to +1.425%;
cumulative changes range from -3.474% to +1.693%. No case reaches the
positive 2.5% review threshold. TSOPF_FS incremental results are -0.251%,
-0.427%, and -0.108% on CPUs 8–15, CPUs 0–7, and CPU 8; LeGresley_2508
measures -1.202%. These results do not establish zero loss or cross-machine
equivalence; no new full corpus screen was run.

Frozen binaries rotate through one hash-verified execution path with
pinned H100 entrywise-changed refactors and no concurrent builds or source
edits. The final audit verifies revision, exact pending diff, build
provenance, runner/binary/matrix hashes, expected jobs/counts, complete
commands, successful exits, residuals and recomputed medians.

Artifacts: `build/prep-trim-repair-dBcoha/analyze-trace-trim-*`; runner:
`validate-analyze-trace-trim.py`; preceding frozen executable:
`prestatic-trace-accepted-kls_bench`.

The broader low-risk cleanup goal remains active.
