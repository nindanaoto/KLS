# Keep pipeline progress, remove scalar work counters

Preceding baseline: `d5fad7c`; resumed-pass baseline: `372f9d1`.
Remove the five remaining scalar dependency/run progress counters, their
private trace struct and aggregation helper, per-worker trace state/resets,
and trace-only parameters through the private update call chain. Remove
the diagnostic run-entry accumulation, but keep the operational run-row
count and output. This removes 72 net C-source lines.

Numerical update loops, validation guards, admission, solver statistics,
locking and failure diagnostics remain. Progress reporting keeps event,
completed/total and begin/end positions at the existing cadence. Archived
work-counter parsing remains and is labeled historical. The SNB entry
remains 64-byte aligned (0x64c80, size 0x1931).

## Validation

Release and ASan/UBSan builds are clean, and CTest passes 6/6 in each.
A release smoke enables tracing, active catch-up and supernode producer
batching. It passes and emits 24 progress-only records matching the expected
event/completed/begin/end format. The trace summarizer processes the log
successfully and produces valid JSON. This is not an execution-count claim
for producer paths, whose detailed counters are no longer emitted.

All 651 timed launches pass: 36 onetone controls, 495 automatic-policy
controls, 48 forced-first controls, and 72 repeat launches. Artifacts are
`build/prep-trim-repair-dBcoha/scalar-progress-*`; runner:
`validate-scalar-progress-trim.py`; preceding frozen binary:
`panel-trace-accepted-kls_bench`. All three versions launch through a shared
temporary executable path with copy/hash verification before each launch.
Rotating pinned H100 lifecycles verify every changed refactor at 1e-8.
Build provenance, commands, binary/matrix hashes and exact source diffs
are checked or recorded. No builds or tests overlap timings.

Selected median paired lifecycle changes versus preceding / resumed baseline:

| Case | Incremental | Cumulative |
|---|---:|---:|
| onetone1, CPU 8 | +0.14% | -0.22% |
| onetone1, CPU 0 | -0.10% | -0.06% |
| TSOPF_FS_b9_c6, CPUs 8-15 | -0.17% | +0.73% |
| TSOPF_FS_b9_c6, CPUs 0-7 | +0.21% | +0.21% |
| TSOPF_FS_b9_c6, CPU 8 | +0.01% | -0.14% |
| twotone, CPUs 8-15 | +0.28% | +1.04% |
| ASIC_680ks, CPUs 0-7 | +0.29% | -0.27% |
| mimo46x46_system, CPUs 0-7 | +1.48% | -0.64% |
| LeGresley_2508, CPUs 8-15 | +0.84% | -1.19% |
| circuit_4, forced first, CPUs 8-15, repeat | -0.19% | +0.24% |

The initial forced-first `circuit_4` result was -2.15% incremental but +2.73%
cumulative. Its independent 24-triplet repeat cleared that flag. No
repeatable approximately 3% stopping regression was established. Retain
this chunk; the broader cleanup goal remains active. These focused controls
do not constitute a new full-corpus or cross-machine campaign.
