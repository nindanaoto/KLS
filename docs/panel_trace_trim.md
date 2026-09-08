# Remove panel progress-trace diagnostics

Preceding baseline: `dae56c1` (numeric source `78a4c89`); resumed-pass
baseline: `372f9d1`. Remove five panel update/append progress-trace counters,
aggregation, printing and before/after work estimates. Their last consumer
also disappears for `kls_row_first_published_u_run_entries`, so remove that
private helper. This removes 72 net C-source lines. Separate solver panel
statistics, panel cache contents/publication, admission predicates,
numerical updates and locking remain unchanged. Archived progress-field
parsing is retained and labeled historical.

## Correctness

Final Release and ASan/UBSan builds are clean; CTest passes 6/6 in each.
The initial build identified the now-unused helper, removed before final
validation. A release smoke enables tracing, active catch-up and supernode
producer batching. It passes with 24 general trace records and no retired
panel fields in those records. It does not establish producer execution
counts, whose trace fields were removed in the preceding cleanup.

## Performance

All 687 timed launches pass: 36 onetone controls, 495 automatic-policy
controls, 48 forced-first controls, and 108 confirmation launches.
Automatic controls add `mimo46x46_system` on CPUs 0-7 and `LeGresley_2508`
on CPUs 8-15, based on the largest small increases in the corrected corpus
screen. These are additional controls, not matrix-specific solver changes.

Runner: `build/prep-trim-repair-dBcoha/validate-panel-trace-trim.py`.
Artifacts: `panel-trace-trim-*` and `panel-trace-repeat*` in that directory.
The preceding frozen binary is `producer-trace-accepted-kls_bench`.
All three versions execute through one shared temporary path, after a
sequential copy and hash check before each launch. Build provenance,
commands, matrix/binary hashes and exact source diffs are recorded or
verified. Rotating pinned H100 lifecycles check every changed refactor at
1e-8; builds/tests do not overlap timings.

Key median paired lifecycle changes versus preceding / resumed baseline:

| Case | Incremental | Cumulative |
|---|---:|---:|
| onetone1, CPU 8 | +0.07% | -0.14% |
| onetone1, CPU 0 | -0.18% | +0.15% |
| TSOPF_FS_b9_c6, CPUs 8-15 | -0.42% | +0.09% |
| TSOPF_FS_b9_c6, CPUs 0-7 | +0.28% | -0.07% |
| TSOPF_FS_b9_c6, CPU 8 | +0.06% | -0.21% |
| ASIC_680ks, CPUs 0-7 | -0.77% | -0.43% |
| mimo46x46_system, CPUs 0-7 | -1.38% | -1.40% |
| LeGresley_2508, CPUs 8-15, repeat | -1.61% | -2.01% |
| twotone, CPUs 8-15, repeat | -0.42% | -0.20% |

Two initial flags did not persist. `LeGresley_2508` initially measured
+4.42% / +2.74%; its 24-triplet repeat produced the result above.
`twotone` initially measured +2.86% / +3.18% with a large initial-factor
outlier; its 12-triplet repeat also cleared the flag. Retain this chunk:
no repeatable approximately 3% stopping regression was established.
This focused validation is not a new full-corpus or cross-machine claim.
