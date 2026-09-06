# TSOPF_FS_b9_c6 build-regression repair (2026-09-06)

## Cause and repair

The stopped slimming campaign compared different SPRAL builds. Both KLS C
compilations used the same Release flags, but the primary build cache had empty
`CMAKE_Fortran_FLAGS_RELEASE`, leaving bundled SPRAL unoptimized. The rebuilt
`4d32410` baseline used `-O3`. Matching runs during block-ordering preparation,
so this discrepancy appeared inside the initial-factor timer, not analysis.

Diagnostic library swaps transferred the initial-factor penalty between KLS
versions: at one thread, the baseline measured approximately 10.18ms with
optimized SPRAL and 12.46ms without optimization; the slimmed version measured
10.17ms and 12.57ms respectively. This establishes the build cause of the
repeatable approximately 2.4ms penalty. The earlier 30% lifecycle slowdown did
not reproduce consistently and must not be attributed wholly to this penalty.

The local cache was repaired and the Release artifacts rebuilt with:

```sh
cmake -S . -B build -U 'CMAKE_Fortran_FLAGS*' -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 8
ctest --test-dir build --output-on-failure
```

SPRAL now uses `-O3 -DNDEBUG -O3 -Jspral_modules -fPIC -march=native
-fno-tree-loop-vectorize`, identical to the baseline. The target-local
vectorization guard is preserved. No forced optimization was added to Debug
builds or to project defaults; this was a stale local cache override.

All slimming remains intact. The solver-source diff exactly matches the
previous campaign's preserved patch. No numerical kernel, threshold, worker
layout, specialization, or public C API changed.

## Targeted validation

Comparison: `4d32410` versus `9d8df18` plus the existing uncommitted slimming,
with matching build provenance across eight linked targets. Each launch runs
100 systems with entrywise updates of amplitude 0.001 and verifies every
refactor at a relative-residual limit of `1e-8`. Auxiliary library pools use
one thread; before/after launch order alternates. A launch has a 90-second
timeout and an 80-GiB address-space limit. No failures or outliers were removed.

All **312 launches passed**. Maximum observed refactor residual was
`5.45371652e-9` (TSOPF_FS_b9_c1); TSOPF_FS_b9_c6 stayed below `1.24e-14`.

TSOPF_FS_b9_c6 initial-factor times are marginal medians. Lifecycle changes
are medians of paired after/before ratios, not ratios of marginal medians.
Positive percentages mean slower.

| Configuration | Pairs | Initial factor before → after | Lifecycle change |
| --- | ---: | ---: | ---: |
| 1 thread, CPU 8 | 12 | 10.000 → 10.035ms | +1.42% |
| 8 threads, CPUs 8–15, 32MiB LLC | 40 | 12.680 → 12.702ms | −0.70% |
| 8 threads, CPUs 0–7, 96MiB LLC | 40 | 12.278 → 12.329ms | +0.46% |

Controls used eight pairs per cache domain:

| Matrix | 32MiB lifecycle change | 96MiB lifecycle change |
| --- | ---: | ---: |
| TSOPF_FS_b9_c1 | −0.93% | +0.02% |
| TSOPF_RS_b9_c6 | +1.58% | +0.59% |
| onetone1 | −0.43% | −0.25% |
| twotone | −0.51% | −0.97% |

Initial-factor performance is restored. Small lifecycle differences remain;
these results do not establish exact parity or justify speculative code tuning.
The full campaign remains stopped, and this is not a full-corpus no-regression
claim.

## Prevention and evidence

`scripts/run_slim_regression.py` now checks the benchmark and its linked
CMake targets, including Fortran flags, compiler identities, backend defines,
dependency revisions, external objects, and runtime-library hashes. It rejects
different or missing provenance and stale artifacts before any benchmark
launch. Only source/build-root paths are normalized; optimization flag order
is significant. This check requires Unix Makefiles build metadata, as did the
runner's original flags collection. The CLI and historical `--reduce-only`
behavior are unchanged.

Eight unit tests cover matching builds at different paths, Fortran-only
mismatches, flag ordering, missing metadata, stale artifacts, dependency
revision mismatches, and reduction without build metadata. Release CTest
passes all five registered tests, including this suite.

Local artifacts in `build/tsopf-build-fix-20260906/` contain frozen binaries,
their hashes, build provenance, the targeted harness, raw `observations.jsonl`,
`comparison.csv`, and `phases.json`. Baseline SHA-256:
`21909137572d86a54b0c1f01019b4397fa0b5f1532a1c4aa053355b6b7fce2da`.
Repaired SHA-256:
`9b013e1e6ee5dec51d6948e66628577f89745a70bd87dee459607ab3eda42f34`.

The stopped campaign's README now identifies its build confound. Its original
frozen binaries, observations, and summaries remain preserved, not corrected
retroactively or reused as isolated measurements of slimming.
