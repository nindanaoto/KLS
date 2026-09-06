# Three optional-branch removals (2026-09-06)

Removed relative to `2463f0d`:

- The opt-in parallel reciprocal-growth worker, its launcher, dispatch mode,
  and `KLS_ENABLE_PARALLEL_CONTRACT_RGROWTH` control. Ordinary KLU growth
  classification and refinement safeguards remain. The shared reduction
  scratch is still used by production solve/contract workers and remains.
- The extra opt-in serial NodeNDP portfolio, its thread-local override and
  exclusive argument, and both `KLS_ENABLE_GENERIC_SERIAL_NDP_PORTFOLIO` and
  `KLS_DISABLE_GENERIC_SERIAL_NDP_PORTFOLIO`. Production METIS ordering and
  deterministic parallel NodeNDP remain; the former permission gate is
  unconditionally true without the removed override.
- The diagnostic preselected-row/direct-KLU challenger and its
  `KLS_ENABLE_PRESELECTED_ROW_DIRECT_KLU_TOURNAMENT` control, including its
  exclusive cache-size predicate. The production column challenger remains.

The source change removes **293 net lines**. No public C API changed, and
the previous slimming is retained. These branches were disabled in the paper
configurations, but deleting their compiled code is not necessarily
performance-neutral.

## Verification

Final Release and ASan/UBSan CTest both pass **5/5**. The final builds have
no new compiler warnings, and `git diff --check` passes. Searches confirm no
remaining references to the removed controls, mode, functions, or TLS flag.

The targeted benchmark freezes the before/after binaries and compares all
linked-target build provenance, including optimized SPRAL. Each launch runs
100 systems, entrywise changes of amplitude 0.001, and verification of every
refactor at residual limit `1e-8`. Auxiliary library pools use one thread.
Before/after launches alternate and only one solver runs at a time.

All **144 initial launches** and **40 independent recheck launches** passed.
Initial paired lifecycle changes (positive means slower):

| Matrix | Pairs per domain | 32MiB LLC, 8 threads | 96MiB LLC, 8 threads |
| --- | ---: | ---: | ---: |
| TSOPF_FS_b9_c6 | 12 | +1.88% | +0.18% |
| TSOPF_FS_b9_c1 | 4 | −0.58% | +1.07% |
| TSOPF_RS_b9_c6 | 4 | −0.71% | −2.18% |
| onetone1 | 4 | +0.78% | +1.62% |
| twotone | 4 | +1.98% | +2.69% |
| rajat25 | 4 | +1.43% | +1.27% |

TSOPF_FS_b9_c6 at one thread showed **+3.75%** in eight initial pairs and
**+3.65%** in an independent twenty-pair recheck. Recheck phase medians:

| Phase | Before | After |
| --- | ---: | ---: |
| Initial factor | 10.090ms | 10.073ms |
| First refactor | 72.689ms | 80.314ms |
| Steady refactor | 0.984ms | 0.982ms |
| Full 100-system lifecycle | 204.740ms | 212.195ms |

The repeatable one-thread penalty is concentrated in the first refactor; its
lower-level cause has not been established. The removals are retained as
requested, but this is **not a no-performance-regression result**. No padding,
new specialization, or speculative kernel tuning was added to conceal it.
The small control sample sizes also do not establish full-corpus parity.

Artifacts: `build/optional-trim-20260906/` contains provenance, hashes,
frozen binaries, source diffs, the harnesses, raw observations, comparisons,
phase summaries, and the independent recheck. After the initial batch, a
stale comment was removed and the project rebuilt/retested. The measured
and final Release binaries have identical `.text` bytes; the recheck uses
the final binary. The full paper campaign remains stopped.
