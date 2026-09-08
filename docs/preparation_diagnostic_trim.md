# Preparation diagnostics and supernode override cleanup

Relative to `254ea42`, this batch removes 258 net source lines. No production
kernel, automatic eligibility threshold, scheduling dependency, scratch
alignment, or public statistics field is removed.

## Removed

- Unread `seconds` fields and their clock calls in group-successor,
  segment-input-target, and dense-group-classification preparation jobs.
- The diagnostic row-supernode census, its temporary allocations, and its
  once-only bookkeeping field. Production supernode discovery remains.
- Preparation phase timestamps and printing controlled by
  `KLS_TRACE_LEAN_PREP` and `KLS_TRACE_DIRECT_ROW_PATTERN`.
- `KLS_LEAN_SNODE2_ONLY` and `KLS_ENABLE_LEAN_SNODE`. Automatic supernode
  candidate/capability checks and width-2 through width-16 kernels remain.
  The former forced-path bypass now always follows the default policy.
- The trace-only separator-private work-imbalance scan and printing.
  Actual worker-work accounting and partitioning remain.
- The row/column lifecycle cycle maxima, their resets, updates, redundant
  guards, and trace arguments. The remaining lifecycle trace reports minima.

The six retired environment switches are the two manual overrides above and
`KLS_TRACE_ROW_SNODE`, `KLS_TRACE_LEAN_PREP`,
`KLS_TRACE_DIRECT_ROW_PATTERN`, and `KLS_TRACE_EGRAPH_SEPARATOR_PRIVATE`.
They are ignored. Historical reports mentioning them are not rewritten.

## Lifecycle guard equivalence

Each side initializes its minimum and maximum from the same first cycle.
Subsequent ordered samples can only decrease the minimum or increase the
maximum. Thus a nonpositive maximum implies a nonpositive minimum, and the
maximum's `<= 0` guard adds no rejection. Both sample-count guards remain.

If the initial cycle is NaN, both accumulators remain NaN and both comparisons
are false; later NaNs change neither accumulator. Infinities and signed zeros
also preserve guard equivalence. No new clock validity policy is introduced.
An exhaustive model check over five-sample sequences drawn from negative
infinity, -1, signed zeros, 1, positive infinity, and NaN passed all 67,228
prefix checks. Actual selection continues to use the same cycle minima,
sample windows, comparison margins, and state transitions.

## Validation

Artifacts: `build/prep-diagnostic-trim-MGfCA6/`.

- Release and ASan/UBSan CTest: 5/5 each; successful builds without warnings.
- Release smoke with all six retired switches set: passed, with no retired
  diagnostic messages.
- Initial focused sweep: 172/172 valid launches over nine matrices and
  fourteen matrix/CPU configurations.
- Rechecks: 60/60 valid launches; total 232/232 valid, maximum initial or
  changed-value relative residual 2.9406311e-12.

The baseline executable was frozen before editing, from `254ea42`. The
runner checks identical build/dependency/runtime provenance, records binary
and matrix hashes and exact commands, alternates before/after order, and
validates all 100 changed-value refactors with a residual limit of 1e-8.
CPU sets are 8 (one thread), 8-15 (eight threads, 32 MiB domain), and 0-7
(eight threads, 96 MiB domain). Each launch has a 90-second timeout and an
80 GiB address-space limit. Benchmarks run after builds and tests finish.

This is a focused regression check, not a full paper or Xyce campaign.

### Initial paired H100 lifecycle changes

Positive means slower. Percentages are medians of paired percentage changes,
not ratios of marginal medians. Baseline is `254ea42`, not an older full-paper
baseline; this does not establish recovery of any earlier regressions.

| Matrix | Configuration | Pairs | Change |
|---|---|---:|---:|
| TSOPF_FS_b9_c6 | t1 | 6 | +0.91% |
| TSOPF_FS_b9_c6 | t8, 32 MiB | 6 | +1.45% |
| TSOPF_FS_b9_c6 | t8, 96 MiB | 6 | +3.18% |
| 1138_bus | t8, 32 MiB | 16 | -0.76% |
| 1138_bus | t8, 96 MiB | 16 | +0.46% |
| transient | t1 | 6 | +2.68% |
| transient | t8, 32 MiB | 3 | -0.16% |
| transient | t8, 96 MiB | 3 | +0.66% |
| onetone1 | t1 | 6 | +7.04% |
| twotone | t1 | 6 | +5.98% |
| ASIC_320k | t1 | 3 | +2.65% |
| ASIC_680ks | t1 | 3 | +4.61% |
| bcircuit | t8, 32 MiB | 3 | +0.72% |
| rajat03 | t8, 32 MiB | 3 | +0.39% |

### Warning-case repeats

These use the identical binaries, not a modified follow-up implementation.

| Matrix | Configuration | Additional pairs | Change | Before H100 (s) | After H100 (s) |
|---|---|---:|---:|---:|---:|
| TSOPF_FS_b9_c6 | t8, 96 MiB | 16 | +2.49% | 0.077450 | 0.079136 |
| transient | t1 | 6 | +2.61% | 4.969013 | 5.096764 |
| onetone1 | t1 | 4 | +7.21% | 6.619444 | 7.089964 |
| twotone | t1 | 2 | +6.01% | 19.373092 | 20.537091 |
| ASIC_680ks | t1 | 2 | +4.58% | 12.197738 | 12.756072 |

The before/after time columns are marginal medians; the percentage column
uses paired changes as above. The smaller two-pair repeats reinforce the
initial samples but do not constitute a broad performance characterization.

SHA-256 of the frozen baseline (also matches the previous accepted binary):
`f177e76b8ce7e0a38fc45557b9f5afcbdfee15001ae7e2bc000aca2e97b963f7`.
SHA-256 of the cleanup binary (matches both benchmark snapshots and the
current Release binary):
`a2bfc9e8863385a343284864f23c00395c3595af82083ffe5d62ad541e3eb8cd`.

The requested removals are retained, but this version is **not verified
performance-neutral**. No unrelated policy tuning, matrix exceptions, or
padding were introduced to conceal the regressions. Their cause has not
been established by this cleanup validation.
