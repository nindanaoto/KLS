# Targeted trace and residual-statistics cleanup

Baseline: pushed commit `88a177e`. The TSOPF transformed-residual fix and
earlier numerical slimming remain intact.

## Removed redundancy

- Remove `KLS_TRACE_DENSE_HELP_GROUP`, its cached parser, selected-group
  timing, logging-only locals, and DHG output. Keep the dense-help ownership
  protocol, wait deadline, synchronization, and general `KLS_TRACE_DENSE_HELP`
  telemetry. Operational clocks are not diagnostic-only and remain.
- Remove `KLS_TRACE_ROW_TASK_ROOT`, its parser and RTG output. Keep general
  worker timing and maximum-group telemetry. Its existing timing predicate
  now suffices, so the nested duplicate telemetry-enabled check is removed.
- Remove the second pair of `auto_scale_deferred`/`tight_pivot_deferred`
  resets in `clear_matrix`; the earlier assignments remain.
- Remove the obsolete `KLS_DISABLE_ROW_SELF_CHECK_L2_CONTRACT` alias. The
  canonical `KLS_DISABLE_ORDINARY_SELF_CHECK_L2_CONTRACT` control remains.
  This is an intentional internal-environment compatibility change, not
  removal of the default L2 accuracy check.
- Replace three identical worker-statistics reductions with one always-inline
  helper. Thread traversal, floating-point addition order, maximum comparison
  semantics, and output-store order are preserved. Each caller still owns
  and releases its results buffer; no allocation policy changes.

Net library-source reduction: **109 lines**. No numerical engine, solver
field, public API, or policy threshold is removed.

## Validation

Release and ASan/UBSan CTest each pass 6/6. The transformed residual test
continues to exercise both index widths, changed values, transpose and
multiple RHS. No retired control references remain in source/tests/scripts.
`git diff --check` passes.

Artifacts: `build/prep-trim-repair-dBcoha/trace-trim*`. The alternating-order
paired H100 comparison uses the frozen `residual-slim-kls_bench` executable
equivalent to `88a177e`; the runner calls it `original`. It checks build
provenance and records commands, hashes and source diff. Automatic policies
and entrywise changes are used, with 1e-8 residual verification. Timings do
not overlap builds or tests.

The focused campaign completes with **232/232 valid launches** over eleven
matrices and thirteen configurations. TSOPF has 24 pairs per eight-thread
CPU set and eight at one thread. This is not a full paper or Xyce rerun.

Median paired H100 lifecycle changes versus `88a177e` (negative is faster):

| Matrix | Configuration | Change |
|---|---|---:|
| TSOPF_FS_b9_c6 | t8_32 | -0.51% |
| TSOPF_FS_b9_c6 | t8_96 | 0.57% |
| TSOPF_FS_b9_c6 | t1 | 0.41% |
| TSOPF_RS_b9_c6 | t8_32 | 0.15% |
| bcircuit | t8_32 | -0.52% |
| rajat03 | t8_96 | -1.38% |
| bips98_1142 | t8_96 | -0.78% |
| 1138_bus | t8_96 | -0.53% |
| circuit_4 | t8_32 | 1.56% |
| twotone | t8_32 | -0.31% |
| ASIC_680ks | t8_96 | 0.27% |
| onetone1 | t1 | -0.60% |
| transient | t1 | 0.29% |

`t1` uses CPU 8; `t8_32` uses CPUs 8-15; `t8_96` uses CPUs 0-7.
The TSOPF recovery remains intact. The measured control differences range
from -1.38% to +1.56%; these are focused samples, not a guarantee of zero
performance regression on every workload.

An additional smoke run passes with general dense-help and worker timing
enabled, while all three retired environment variables are also set.
It emits 118 general worker/dense-help trace lines and no DHG/RTG lines,
confirming the retired trace controls are ignored and general telemetry
still runs. This diagnostic run is excluded from performance comparisons.

Candidate executable SHA-256:
`0abe4f282ab02603cdecea93cd909e12534d29525ac9022f9548437b09db413b`.
