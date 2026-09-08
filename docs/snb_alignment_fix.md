# Stabilize serial SNB kernel placement after slimming

Baseline: `4b1d4d7` (numerical source `c4aa02e`). Cumulative performance
reference: `88a177e`. All accepted slimming remains intact.

## Cause and change

A rotating comparison of eight saved revisions found the onetone1 slowdown
in steady serial SNB refactorization, not initialization or solving. The
SNB numerical body retained an identical normalized 1,401-instruction
sequence, but its placement shifted as unrelated diagnostics disappeared.
Samples concentrated in the SIMD accumulation in `kls_snb_update_edge`,
inlined into `kls_snb_refactor`.

A controlled assembly-padding experiment changed only that routine's
address across four variants (all other function/data symbol addresses and
sizes fixed). A 12-round confirmation measured +2.53% H100 lifecycle for
the retained baseline, -0.05% for favorable placement, and +3.13% after
moving identical arithmetic another 16 bytes, versus `88a177e`. All 48
launches passed numerical checks. This establishes code-placement
sensitivity, without identifying a specific cache or branch-predictor
mechanism: hardware counters are unavailable in this environment.
Diagnostic artifacts: `build/onetone-cause-ilIWQW/`.

The production fix adds only guarded `__attribute__((aligned(64)))` to
`kls_snb_refactor`, using the GCC/Clang guard already used elsewhere in KLS.
It stabilizes the hot routine's entry alignment against changes in preceding
code size. There are no matrix names, size exceptions, CPU-model dispatch,
numerical-policy changes, or diagnostic restorations. Other compilers retain
their existing behavior. This is a placement hint, not a claim that every
compiler or machine will have identical performance. Internal code changes
can still change the hot loops' offsets and require ordinary regression tests.

## Validation

Release and ASan/UBSan CTest each pass 6/6. The release entry point moves
from `0x64fd0` to `0x65000`, with the same 0x1931-byte body and identical
normalized instruction sequence. Source diff checking passes.

The focused H100 campaign rotates preceding (`original`), cumulative
(`start`), and candidate (`fix`) executables with clean environments and
pinned CPU sets. Matching build provenance, source diff, commands, hashes,
JSONL results and summaries are under
`build/prep-trim-repair-dBcoha/snb-align-*`. Every entrywise-changing
refactor is checked at a 1e-8 residual limit. Timings do not overlap builds
or tests. The onetone1 check includes one thread on CPU 8 and on CPU 0,
covering this machine's two cache domains. This is not a full paper or Xyce
benchmark, nor validation on another machine.

All **375/375 launches pass**: 36 onetone1 comparisons (six rotating triplets
per CPU domain), 303 automatic-policy controls across six other matrices
and nine configurations, and 36 forced-KLS-first controls (1138_bus,
circuit_4, bcircuit). The recorded source diff matches the final source patch.

Median paired H100 lifecycle changes (negative means faster):

| Matrix / configuration | vs slimmed baseline | vs 88a177e |
|---|---:|---:|
| onetone1, 1 thread, CPU 8 / 32 MB L3 domain | -2.24% | -0.36% |
| onetone1, 1 thread, CPU 0 / 96 MB L3 domain | -3.16% | -0.52% |
| TSOPF_FS_b9_c6, 8 threads, CPUs 8-15 | +0.75% | +0.39% |
| TSOPF_FS_b9_c6, 8 threads, CPUs 0-7 | -0.14% | +1.27% |
| TSOPF_FS_b9_c6, 1 thread | -0.49% | +0.08% |
| transient, 1 thread | -1.49% | -0.02% |
| twotone, 8 threads | -1.94% | -2.84% |
| ASIC_680ks, 8 threads | -0.43% | -0.19% |

Onetone1's median steady refactor time falls from 56.744 to 55.378 ms
on CPU 8, and 57.762 to 55.721 ms on CPU 0. These are marginal medians;
the percentage table uses paired changes, not ratios of marginal medians.
The cumulative onetone1 slowdown is recovered on both cache domains.

No focused control reaches the approximately 3% lifecycle regression gate.
The largest incremental slowdown is +0.75%; the largest cumulative slowdown
is +1.27%. Forced-KLS-first changes span -1.39% to +0.25% incrementally.
These fixtures check first-factor correctness/layout effects but do not
establish timed parallel-row-pipeline coverage. No universal or cross-machine
performance guarantee is inferred from this focused validation.
