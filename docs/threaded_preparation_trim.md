# Optional threaded preparation cleanup (2026-09-06)

Relative to `84c7dae`, this batch removes thirteen explicitly opt-in
threading/prewarming variants from `src/kls_numeric_prepare_repair.inc`.
The paper campaign configurations do not enable these controls.

## Removed variants

| Preparation variant | Removed control |
| --- | --- |
| Internal row-solve partition threading | `KLS_PARALLEL_ROW_SOLVE_PARTITION_BUILD` |
| Dense-producer target-map threading | `KLS_PARALLEL_DENSE_PRODUCER_TARGET_MAPS` |
| Group-work threading | `KLS_ENABLE_PARALLEL_GROUP_WORK_PREP` |
| Row-input cleanup overlap | `KLS_PARALLEL_ROW_INPUT_CLEANUP_PREP` |
| Row-level preparation overlap | `KLS_PARALLEL_ROW_LEVEL_PREP` |
| Dense-group classification overlap | `KLS_PARALLEL_DENSE_GROUP_CLASSIFICATION_PREP` |
| Group-dependency threading | `KLS_ENABLE_PARALLEL_GROUP_DEP_PREP` |
| Group-successor preparation overlap | `KLS_PARALLEL_GROUP_SUCCESSOR_PREP` |
| Row-solve partition preparation overlap | `KLS_PARALLEL_ROW_SOLVE_PARTITION_PREP` |
| Row-schedule prewarming | `KLS_ENABLE_ROW_SCHEDULE_PREWARM` |
| Segment-input target preparation overlap | `KLS_PARALLEL_SEGMENT_INPUT_TARGET_PREP` |
| Dense-producer run construction threading | `KLS_ENABLE_PARALLEL_DENSE_PRODUCER_RUNS` |
| Successor sorting threading | `KLS_ENABLE_PARALLEL_SUCCESSOR_SORT` |

Their exclusive launchers, joins, per-worker allocations, partitioning,
reductions, and thread-error cleanup are removed. Retained preparation
helpers use typed direct calls rather than pthread callback adapters.
The exclusive `KLS_TRACE_ROW_SOLVE_PARTITION_PREP` timing report and prewarm
timing report are also removed.

The default serial construction order remains: dependency count/fill,
classification, level plan, successors, group work/sorting, segment targets,
dense-producer runs/targets, and solve partitions. Dense-producer construction
retains its count/prefix/fill passes and allocation/fill validation. Normal
separator schedule/cache construction, transpose plans, parallel numeric
execution, residual checking, and refinement remain. There are no public
API/statistics changes or new numerical thresholds.

The source diff removes **826 net lines**. The preparation file decreases
from 16,753 to 15,927 physical lines; total `src/` plus `include/` decreases
from 115,707 to 114,881 lines, including comments and blanks.

## Verification

Release and ASan/UBSan CTest each pass **5/5**, including existing row-pipeline,
group-statistics, and transpose parallel-row-solve coverage. Both builds
complete without new warnings. Searches find no remaining references to
the thirteen removed controls in source, tests, scripts, headers, or benchmark
configuration, and `git diff --check` passes.

Artifacts are in `build/prep-trim-f0RxjX/`. The before binary was frozen from
the committed Release build, with SHA-256
`c1eec839bcedbb482bdcd4b8da77a4ff43f7d1fb5fce26d332c4a93e5d45f003`.
The harness requires matching linked-target build provenance, including
optimized SPRAL, and archives frozen binaries, source diffs, hashes, commands,
raw observations, and phase summaries.

Each benchmark launch uses 100 systems, entrywise changes of amplitude
0.001, and verification of every refactor at residual limit `1e-8`. Auxiliary
library pools use one thread. Before/after launches alternate, one solver at
a time, with fixed CPU affinity. The full paper campaign remains stopped.

## Initial paired performance results

All 168 initial launches passed numerical verification. Positive lifecycle
changes mean slower. TSOPF_FS_b9_c6 has 12 pairs per eight-thread domain;
the other matrices have four pairs per domain.

| Matrix | 32MiB LLC, 8 threads | 96MiB LLC, 8 threads |
| --- | ---: | ---: |
| TSOPF_FS_b9_c6 | -8.98% | +7.13% |
| TSOPF_FS_b9_c1 | -1.32% | +0.21% |
| TSOPF_RS_b9_c6 | -3.66% | -3.40% |
| onetone1 | +0.78% | +0.48% |
| twotone | +0.68% | +0.99% |
| rajat25 | +0.85% | +0.74% |

The 20-pair one-thread TSOPF_FS_b9_c6 check measured -0.09%, with median
lifecycle 204.118ms before versus 204.055ms after. The earlier diagnostic
cleanup's one-thread improvement is retained in this check.

The eight-thread TSOPF difference is concentrated in steady refactor rather
than preparation: 0.490ms to 0.407ms on the 32MiB domain, and 0.419ms to
0.471ms on the 96MiB domain. Both retain the row-refactor path; the sets of
reported non-timing/non-accuracy scalar statistics are unchanged within each
domain. The initial observations alone do not establish a persistent
regression or its cause.

Independent 20-pair checks of each eight-thread domain are archived in
`build/prep-recheck-BWQDn8/`. All 80 launches passed. Paired lifecycle changes
were +0.25% (32MiB) and +0.28% (96MiB). Steady refactor medians were
0.415ms versus 0.414ms and 0.421ms versus 0.421ms, respectively. Neither
large initial difference reproduced. Timing variability, including possible
runtime scheduling or layout effects, has not been isolated to a cause.

A final independent confirmation is archived in `build/prep-confirm-WyjiTF/`.
All 80 confirmation launches passed, with paired lifecycle changes of
+0.29% (32MiB) and +0.77% (96MiB). The large initial differences again did
not reproduce. Across all three batches, **328/328 launches passed**.
No repeatable large regression was established in this sample; these bounded
checks do not establish full-corpus parity or prove exact performance equality.
