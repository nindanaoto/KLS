# Remove producer state-union diagnostics

Preceding baseline: `fac2821`; cumulative baseline: pushed `88a177e`.

Remove the trace-only scan that counted total and unique pattern rows across
producer-batch targets. This includes its generation-mark scratch allocation,
overflow reset, cleanup, two private trace counters, aggregation, and printed
fields. Neither count is used by numerical admission or execution. Producer
batching and general pipeline progress tracing remain unchanged. The summary
script still reads these counters in archived logs and documents that current
KLS no longer emits them. Net C-source reduction: 73 lines.

Release and ASan/UBSan CTest each pass 6/6, with no build warnings or errors.
The focused rotating three-binary H100
campaign uses the preceding and cumulative baselines and records provenance,
commands, hashes, source diff, results, and summaries under
`build/prep-trim-repair-dBcoha/state-union-*`. Every changed refactor is
numerically checked at a 1e-8 residual limit. Build/test work and timed runs
do not overlap. Private layouts change, so numerical tests alone are not
evidence of performance neutrality. This is not a full paper or Xyce rerun.

The trace-enabled release smoke test also passes, emitting 24 general
pipeline messages without the removed state-union fields. The summary helper
successfully reads that log. Source-wide checks find no remaining references
to the deleted fields or scratch state in library code, headers, or tests.

All 192 forced-KLS-first comparisons pass across 1138_bus (24 triplets on
each eight-core CPU set), circuit_4 and bcircuit (eight triplets each).
Median paired lifecycle changes span -0.49% to +0.99% incrementally and
-2.42% to +0.45% cumulatively. The forced fixtures exercise KLS-first but
do not establish timed parallel-pipeline coverage.

All 309 automatic-policy comparisons pass across seven matrices and ten
configurations, including 12 ASIC_680ks triplets. Incremental lifecycle
changes span -0.25% to +1.07%; cumulative changes span -1.17% to +2.04%.
TSOPF_FS_b9_c6 is +0.51% and +0.11% on the two eight-core CPU sets and
+0.31% on one thread versus the preceding commit. Serial onetone1 is +0.76%
and transient +0.23%. No measured case reaches the approximately 3% lifecycle
stop threshold. The source diff matches the timed campaign's recorded diff.
This chunk is retained; focused sampling does not prove zero regression on
all workloads.

Further audit found optional ETree repair phase timing suitable for the next
diagnostic-removal trial. In contrast, the producer rejection enum is not
entirely diagnostic: NOT_ROOT drives active catch-up scheduling. Do not
delete that behavior while removing trace accounting. Neither candidate is
changed in this chunk.
