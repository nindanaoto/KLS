# Remove producer internal/output touch diagnostics

Preceding baseline: `e4039b9`; resumed-pass baseline: pushed `372f9d1`.
Retain the SNB alignment recovery and all previous accepted slimming.

Remove the producer internal/output touch counters, their private output
parameters, worker-trace fields, aggregation and printed values. The two
helper callers used these totals only to update trace accounting. The
inner-loop counts were computed even when tracing was disabled. Actual
updates, bounds checks, dependency discovery, locking, ownership flags,
detached target copies, failure cleanup, and catch-up admission are unchanged.
Net C-source reduction: 35 lines. The trace summary helper keeps archived
field support and documents that current code no longer emits them.

Release and ASan/UBSan CTest each pass 6/6. No removed counter/output names
remain in library code, headers, or tests. Source diff checking passes.
The SNB entry remains at 0x65000 with size 0x1931. A trace-enabled smoke
run and the trace-summary helper exercise current general progress output.

Focused H100 artifacts live under `build/prep-trim-repair-dBcoha/producer-touch-*`.
The rotating three-binary runner labels the preceding commit `original`,
the aligned resumed-pass baseline `start`, and the candidate `fix`. It
checks build provenance, records hashes/commands/source diff, pins CPU sets,
cleans the environment, and verifies each entrywise-changed refactor at a
1e-8 residual limit. Builds/tests do not overlap timed runs. This is not a
full paper/Xyce rerun or another-machine benchmark.

All **435/435 launches pass**: 36 onetone1 comparisons (six triplets on
each cache domain), 351 automatic-policy comparisons across six other
matrices/nine configurations, and 48 forced-KLS-first controls. Serial
TSOPF uses 24 triplets from the outset because its previous eight-triplet
screen was noisy. The campaign's recorded source diff matches this patch.

Median paired H100 lifecycle changes versus `e4039b9` / `372f9d1`:

| Matrix / configuration | Incremental | Resumed-pass cumulative |
|---|---:|---:|
| onetone1, CPU 8, one thread | -0.03% | +0.11% |
| onetone1, CPU 0, one thread | -0.17% | +0.09% |
| TSOPF_FS_b9_c6, CPUs 8-15, eight threads | +0.35% | -0.16% |
| TSOPF_FS_b9_c6, CPUs 0-7, eight threads | -0.41% | +0.79% |
| TSOPF_FS_b9_c6, one thread | -0.14% | +0.02% |
| transient, one thread | +0.26% | +0.07% |
| twotone, eight threads | -1.60% | +0.09% |
| ASIC_680ks, eight threads | +0.10% | +1.83% |

Forced-first controls span -0.79% to +0.61% incrementally. No case meets
the approximately 3% repeatable lifecycle stop condition. Retain this chunk;
focused samples do not establish zero regression on every workload or timed
parallel-pipeline coverage. The trace-enabled smoke passes with 24 general
pipeline messages and no retired producer fields; the parser accepts it.
Release/sanitizer builds have no warnings or errors.

The aligned SNB numerical body is unchanged; the only disassembly difference
within kls_snb_refactor is relocation of its call to the prelude, which moves
from 0x64be0 to 0x64bd0. The alignment recovery remains intact in timings.

Next audit candidate: producer rejection-reason diagnostic counters and
their enum plumbing. Only NOT_ROOT is consumed by an operational catch-up
decision; its semantics must be retained if the other diagnostic distinctions
are removed. This candidate is not changed or validated in this chunk.
