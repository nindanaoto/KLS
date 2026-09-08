# Remove producer candidate-work diagnostics

Preceding baseline: `cb2d7d9`; resumed-pass baseline: `372f9d1`.
Remove 58 C-source lines for four candidate-target/rejected-work counters,
their increments, aggregation and trace output. Remove the now-unused
scalar candidate-work estimate. Keep the supernode estimate because it
still feeds successful-batch work tracing. Archived-field parsing remains
available and is labeled historical in the summarizer.

Admission predicates, target counts, saved-stream calculations and batching
thresholds are unchanged, as are catch-up, allocation checks, ownership,
locking and numerical updates. The SNB entry remains 64-byte aligned
(0x64d40, size 0x1931).

## Validation

Release and ASan/UBSan CTest each pass 6/6. A release smoke with tracing,
active catch-up and supernode producer batching enabled passes with 24
general trace records and no retired fields. It observes 157 successful
catch-up targets from 157 attempts. Producer batch counts remain zero;
this smoke does not establish executed supernode batches.

All 507 timed launches pass: 36 onetone controls, 423 automatic-policy
controls and 48 forced-first controls. Artifacts and exact source diffs:
`build/prep-trim-repair-dBcoha/candidate-trim-*`; runner:
`validate-candidate-trim.py`. The preceding frozen binary is
`producer-probe-accepted-kls_bench`. Build provenance, commands and hashes
are recorded and checked. Three binaries rotate through pinned H100
lifecycles, with every entrywise-changed refactor verified at 1e-8. Builds
and correctness tests finish before timing begins.

Median paired lifecycle changes versus preceding / resumed-pass baseline:

| Case | Incremental | Cumulative |
|---|---:|---:|
| onetone1, CPU 8 | -0.28% | -0.26% |
| onetone1, CPU 0 | -0.21% | -0.02% |
| TSOPF_FS_b9_c6, CPUs 8-15 | +0.69% | +1.10% |
| TSOPF_FS_b9_c6, CPUs 0-7 | +0.81% | +0.80% |
| TSOPF_FS_b9_c6, CPU 8 | +0.31% | +0.15% |
| TSOPF_RS_b9_c6, CPUs 8-15 | +0.37% | -0.97% |
| ASIC_680ks, CPUs 0-7 | -2.08% | -0.77% |
| gemat11, CPUs 0-7 | -0.16% | +0.21% |
| rajat25, CPUs 8-15 | -0.24% | +0.51% |
| 1138_bus, forced first, CPUs 8-15 | -1.14% | -0.76% |

Across all tested cases, incremental changes span -2.08% to +0.94%, and
cumulative changes span -0.97% to +1.10%. No case reaches the conservative
2.5% repeat trigger or the approximately 3% notable-regression stop rule.
Retain this chunk. These focused measurements do not replace the wider
corpus screen or establish cross-machine neutrality. The goal remains active.
