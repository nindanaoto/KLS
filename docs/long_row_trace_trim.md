# Remove per-long-row tracing

Preceding baseline: `e1ef404`. Cumulative baseline: pushed `88a177e`.

Remove live and completed long-row diagnostic messages, their private state,
and `KLS_TRACE_ROW_PIPELINE_LONG_ROW_ENTRIES` parsing. None of the removed
values controls numerical work. The general pipeline progress trace remains.
The benchmark runner no longer accepts `--failure-trace-long-row-entries`
or injects that environment variable. The trace summary helper retains its
archived-log parser, now explicitly documented as historical support.
Historical experiment reports are not rewritten.

Reduction: 165 C-source lines and 17 benchmark-runner lines. Private structure
layout and generated code can change, so numerical equivalence alone is not
the performance gate.

Release and ASan/UBSan CTest each pass 6/6. Both Python scripts compile and
the runner's `--help` succeeds. A trace-enabled release smoke run passes,
emits 24 general pipeline trace messages, and emits no per-long-row messages
even with the removed environment variable set to 1. No removed field or
option references remain under `src`, `scripts`, or `tests`.

The focused H100 protocol rotates preceding (`original`), cumulative
(`start`), and candidate (`fix`) binaries, with pinned CPU sets, clean
environments, matching build provenance, and numerical checks on all
entrywise-changing refactors at a 1e-8 residual limit. Metadata, exact
commands, hashes, source diff, JSONL results and summaries live under
`build/prep-trim-repair-dBcoha/long-row-*`. Builds, tests, and timed campaigns
run sequentially. This is not a full paper or Xyce rerun.

All 192 forced-KLS-first launches pass across 1138_bus (24 triplets on
each eight-core CPU set), circuit_4 and bcircuit (eight triplets each).
Median paired lifecycle changes span -1.93% to +0.90% versus the preceding
commit and -2.10% to +0.55% versus the cumulative baseline. All report
`kls_first`; these fixtures do not establish timed parallel-pipeline coverage.

All 282 automatic-policy launches pass across seven matrices and ten
configurations. Incremental lifecycle changes span -0.23% to +2.18%, and
cumulative changes span -1.07% to +1.27%. TSOPF_FS_b9_c6 changes by +0.17%
and +0.50% on eight threads, and +0.43% on one thread versus the preceding
commit. Serial onetone1 is unchanged; transient is +0.44%.

ASIC_680ks had the largest incremental change (+2.18%) in its three-triplet
screen. A separate 12-triplet repeat passes 36/36 launches and measures
+1.38% incrementally, +0.37% cumulatively. No repeatable slowdown meets the
approximately 3% lifecycle stop threshold. This chunk is retained; these
focused results do not prove zero regression on all workloads.

The source diff recorded by the campaign matches the reviewed patch.
Next low-risk candidate: the trace-only producer state-union scan and its
scratch allocation. Reference inspection shows its two counts feed progress
output only; it has not been removed or performance-validated in this chunk.
