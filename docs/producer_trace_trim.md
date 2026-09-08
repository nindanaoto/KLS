# Remove remaining producer trace bookkeeping

Preceding baseline: `2a35dcd`; resumed-pass baseline: `372f9d1`.
Remove seven successful-batch and active-catch-up diagnostic counters,
aggregation and output, two unused producer trace pointers, and unused
work estimates. This removes 79 net C-source lines. Archived parsing is
retained and labeled historical. Actual producer batching and catch-up
remain: their admission predicates, retries, numerical calls, allocation
checks, ownership gates and locking are unchanged. The unrelated operational
`before_l_count` in pipeline catch-up remains. The SNB entry stays 64-byte
aligned (0x64d40, size 0x1931).

## Validation

Release and ASan/UBSan CTest each pass 6/6. An additional release smoke
enables tracing, active catch-up and supernode producer batching. It passes
with 24 general trace records and no retired producer fields. Since the
producer counters are removed, this smoke no longer establishes execution
counts for catch-up or producer batches.

All 507 timed launches pass: 36 onetone controls, 423 automatic-policy
controls and 48 forced-first controls. Artifacts:
`build/prep-trim-repair-dBcoha/producer-trim-*`; runner:
`validate-producer-trim.py`; preceding frozen binary:
`producer-candidate-accepted-kls_bench`. Build provenance, binary/matrix
hashes, exact source diffs and commands are checked or recorded. Three
binaries rotate through pinned H100 lifecycles; every entrywise-changed
refactor is verified at 1e-8. Timing does not overlap builds or tests.

Median paired lifecycle changes versus preceding / resumed-pass baseline:

| Case | Incremental | Cumulative |
|---|---:|---:|
| onetone1, CPU 8 | -0.47% | -0.27% |
| onetone1, CPU 0 | +0.16% | +0.26% |
| TSOPF_FS_b9_c6, CPUs 8-15 | -0.92% | -1.35% |
| TSOPF_FS_b9_c6, CPUs 0-7 | +0.36% | -0.06% |
| TSOPF_FS_b9_c6, CPU 8 | +0.44% | +0.29% |
| TSOPF_RS_b9_c6, CPUs 8-15 | -0.85% | -1.83% |
| ASIC_680ks, CPUs 0-7 | +0.47% | +1.69% |
| gemat11, CPUs 0-7 | -0.54% | +0.36% |
| rajat25, CPUs 8-15 | +0.15% | +0.18% |
| 1138_bus, forced first, CPUs 8-15 | +0.20% | +1.75% |

Across all cases, incremental changes span -1.06% to +0.47%, cumulative
changes -1.83% to +1.75%. No case reaches the 2.5% review trigger or the
approximately 3% notable-regression stop rule. Retain this chunk.

There are now 197 fewer C-source lines since the preceding medium-corpus
checkpoint (unreachable fallback, probe counters, candidate-work counters,
and this chunk). Repeat that broader cumulative check before further source
deletions. Focused validation alone does not establish corpus-wide neutrality.
