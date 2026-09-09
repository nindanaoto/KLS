# Remove unread supernode batch counters

Baseline: `d5464b6`. A source-wide reference audit found that every
`kls_snode_trace_*` counter had only its static definition and update sites:
none was printed, queried, or used to decide numerical work. Remove those
counters and the cached trace predicate used exclusively to update them.

The deletion covers ordinary, cached, paired, and multi-column batch paths
in `kls.c` and `kls_egraph_refactor.inc`. Actual batching, admission thresholds,
numeric updates and `KLS_DISABLE_BATCH_CONSUME` retained their original behavior
in this change. The batching override was later removed in the batch/map
control cleanup; default batching remains.
General `KLS_TRACE_SNODE` preparation messages remain. No solver fields or
public statistics are changed. Net source reduction: 80 lines.

Release and ASan/UBSan CTest each pass 6/6; no `kls_snode_trace_*` references
remain anywhere under `src`. The benchmark uses the same 282-launch focused
matrix/configuration set as the preceding compact-census chunk, rotating
three frozen executables: preceding commit (`original`), pushed `88a177e`
(`start`), and this candidate (`fix`). Automatic-policy H100 runs verify all
entrywise-changing refactors at a 1e-8 residual limit. Build provenance,
source diff, exact commands, hashes and results are recorded under
`build/prep-trim-repair-dBcoha/dead-counter*`. Timings do not overlap builds
or tests. The benchmark is not a full paper/Xyce rerun.

All **282/282 launches pass**. No notable regression was observed. TSOPF
eight-thread lifecycle changes are -0.40% / +1.02% versus the preceding
commit and +0.51% / -0.01% versus pushed `88a177e` (CPUs 8-15 / 0-7).
Serial onetone1 is -0.48% incrementally and -0.42% cumulatively; transient
is +0.06% and +0.12%. The paired samples do not establish zero regression
on every workload, but do not meet the stated stop threshold.

The broader low-risk audit is ongoing. Next, inspect first-factor pipeline
storage-accounting fields and their trace-only consumers before deciding
whether that accounting can be removed independently of worker telemetry
and allocation policy. Those candidates are not yet changed or validated.
