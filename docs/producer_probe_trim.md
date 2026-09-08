# Remove producer-probe diagnostic counters

Preceding baseline: `fbb7a94`; resumed-pass baseline: `372f9d1`.
Remove 53 C-source lines: six producer-probe counters, their increments,
aggregation and trace output. The lookahead counter was never incremented.
Keep archived-field support in the trace summarizer and label it accordingly.
Admission predicates, target counts, saved-stream calculations, batching,
catch-up and numerical updates are unchanged. The SNB entry remains aligned
to 64 bytes (0x64d80, size 0x1931).

## Correctness

Release and ASan/UBSan CTest each pass 6/6. An additional release smoke
enables tracing, active catch-up and supernode producer batching. It passes,
emits 24 general trace records and none of the retired fields, and observes
160 catch-up attempts and 160 successful catch-up targets. Producer batch
counts are zero, so this smoke does not establish executed supernode batches.

## Performance

Artifacts: `build/prep-trim-repair-dBcoha/probe-trim-*`, with runner
`validate-probe-trim.py`. Frozen preceding binary:
`producer-fallback-accepted-kls_bench`. The runner checks build provenance,
records commands, hashes and exact source diff, and rotates three binaries
through pinned H100 lifecycles. Every entrywise-changed refactor is checked
at the 1e-8 residual limit. No builds or tests overlap timed runs.

All 651 launches pass: 36 onetone controls, 423 automatic-policy controls,
48 forced-first controls, and 144 repeat launches. Automatic controls include
the preceding focused set plus `gemat11` on CPUs 0-7 and `rajat25` on CPUs
8-15, selected from the completed medium-corpus screen.

Median paired lifecycle changes versus preceding / resumed-pass baseline:

| Case | Incremental | Cumulative |
|---|---:|---:|
| onetone1, CPU 8 | +0.32% | +0.04% |
| onetone1, CPU 0 | +0.08% | +0.20% |
| TSOPF_FS_b9_c6, CPUs 8-15 | +0.17% | -0.36% |
| TSOPF_FS_b9_c6, CPUs 0-7 | -0.63% | -0.57% |
| TSOPF_FS_b9_c6, CPU 8 | +0.12% | +0.06% |
| ASIC_680ks, CPUs 0-7 | +1.28% | +0.56% |
| gemat11, CPUs 0-7 | +0.08% | -0.64% |
| rajat25, CPUs 8-15 | +0.42% | +0.67% |

Two flags received separate 24-triplet repeats. `TSOPF_RS_b9_c6` initially
measured +2.61% incrementally with substantial sample variation; the repeat
was +0.28% incremental / +0.43% cumulative. Forced-first `1138_bus` initially
measured +3.09% cumulatively but -1.00% incrementally; its repeat was +0.78%
incremental / -0.60% cumulative. Neither flag persisted. Retain this chunk;
no repeatable approximately 3% lifecycle stopping regression was established.
This is focused validation, not a fresh full-corpus or cross-machine campaign.
