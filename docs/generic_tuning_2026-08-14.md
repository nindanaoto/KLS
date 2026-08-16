# Generic KLS tuning follow-up (2026-08-14)

This note records the first post-external-validation tuning revision.  The
external-validation comparison itself remains frozen; every result below is
development or post-reveal diagnostic evidence from the same machine.

## Attribution

The targeted six-arm route ablation covered `fd18`, `usps_norm_5NN`, `wang3`,
`ted_A_unscaled`, `k3plates`, `rma10`, and `airfoil_2d` at eight threads with
entrywise 0.001 value updates.  Disabling the row executor or fast first
factor changed the geometric mean by less than about one percent.  The larger
regret came from two generic policy failures:

1. AUTO excluded every graph below 30,000 rows from its already guarded
   repeated-lifecycle NodeND portfolio, regardless of realized symbolic work.
2. A stable unscaled numeric could run three complete scale trial factors on
   the first changed update even when it used the requested diagonal pivot
   sequence exactly and had no nudge, perturbation, or conditioning signal.

## Retained changes

The generic ND portfolio now admits sub-30k graphs to the existing lifecycle
economics test.  Because a compact separator tree can lose locality despite a
symbolic near-tie, the newly admitted path must remove at least five percent
of retained symbolic storage; the existing decisive ten-percent Pareto arm
and the established large-graph policy are unchanged.  A diagnostic
`KLS_DISABLE_SUB_30K_GENERIC_ND=1` switch restores the old cutoff.

AUTO scale trials are now skipped when the retained unscaled factor has zero
off-diagonal pivots, zero pivot nudges, zero KLS perturbations, and estimated
`rcond >= 1e-6`.  Difficult, pivoted, perturbed, or poorly conditioned
numerics remain eligible.  The diagnostic
`KLS_DISABLE_STABLE_UNSCALED_SCALE_GUARD=1` switch disables this guard.

## Targeted timing evidence

All times below are modeled H100 SPICE cycles in seconds.  The old column is
the median of the frozen post-reveal AUTO ablation; the new column is the
median of three runs after both changes.

| Matrix | Old AUTO | New AUTO | Ratio | Final route |
| --- | ---: | ---: | ---: | --- |
| `airfoil_2d` | 0.834 | 0.527 | 0.632 | METIS |
| `wang3` | 7.356 | 4.637 | 0.630 | METIS |
| `fd18` | 2.198 | 2.212 | 1.007 | AMD, unchanged |
| `usps_norm_5NN` | 6.081 | 5.809 | 0.955 | AMD, unchanged |
| `rma10` | 2.930 | 2.901 | 0.990 | METIS, unchanged |
| `ted_A_unscaled` | 0.264 | 0.267 | 1.010 | AMD, unchanged |
| `k3plates` | 0.477 | 0.519 | 1.088 | AMD, unchanged |

The `k3plates` timing difference is run noise rather than a retained route
change: an initial permissive version selected a marginal ND proposal and was
5.8 percent slower in seven counterbalanced pairs.  Its symbolic storage win
was only 0.7 percent, which motivated the frozen five-percent admission
margin.  The final policy rejects that proposal and restores exactly the AMD
fill and flop counts.

The scale guard was isolated with seven counterbalanced forced-METIS pairs on
`wang3`.  Enabling it reduced H100 from 8.851 to 5.230 seconds (median paired
ratio 0.590) and first changed-refactor time from 3.776 to 0.137 seconds.
Both arms retained scale -1, 7,699,744 factor entries, 4.594 billion reported
flops, and the same residual quality.

## Generalization evidence

The 24-family development manifest completed 16 matrices.  Three matrices
timed out at 120 seconds and five reported singular/factor failures; these
failures are retained in `build/generic-policy-dev-enabled.failures`.  Only
one successful sub-30k development matrix, `TEM27623`, reached the new
portfolio.  Its challenger was rejected and both old/new arms ended with the
same METIS numeric.  Seven counterbalanced pairs gave a 1.006 median ratio;
the new path's approximately 0.15-second rejected analysis is the measured
cost of broader exploration on that case.

After freezing the five-percent margin, the group-disjoint 24-family holdout
completed 24/24 without a crash or solver-status failure.  Two sub-30k
holdouts, `crystk02` and `poisson3Da`, selected METIS.  Their corresponding
historical H100-model rows were 7.92 and 4.78 seconds; the frozen-policy run
reported 4.88 and 2.92 seconds.  This historical comparison is directional,
not a counterbalanced A/B, so it is not used to tune the selector.

Every newly affected route was rerun for twenty entrywise updates with
`KLS_BENCH_VERIFY_EACH_REFACTOR=1`:

| Matrix | Max changed-update relative residual |
| --- | ---: |
| `airfoil_2d` | 5.16e-15 |
| `wang3` | 1.67e-15 |
| `crystk02` | 6.04e-14 |
| `poisson3Da` | 2.50e-15 |

The broad holdout harness does not enable that residual audit by default and
therefore must not be described as a 24-matrix accuracy certificate.  In
particular, its initial solve exposed poor quality on unrelated `shyy41`
(relative residual about 5.7e7) and `shermanACb` (about 2.4e-3).  Those routes
were not selected by either new policy, but they remain correctness work.

## Next kernel work

These changes reduce avoidable work; they do not close the CKTSO/SubtreeLU
gap.  The external profile still points to recurring numeric factorization as
the dominant cost.  The next high-value experiment should therefore measure
and improve the persistent refactor representation: cache-friendly row or
segment storage, supernode/panel update batching, and separator-aware task
granularity.  Work-reducing ordering remains important, but another matrix
shape classifier would be weaker evidence than improving those shared
kernels.

## Short-supernode executor selection (2026-08-16)

A later row-kernel profile found a representation-selection error rather than
a missing arithmetic specialization.  The short-supernode descriptor builder
published its run array even when its capability check rejected the fused
worker.  A non-null array then disabled the ordinary hoisted row worker, and
the fallback consumed the rejected descriptors.  Thus a capability rejection
could itself select a slower executor.

Rejected descriptors are now discarded, so those factors retain the ordinary
hoisted worker.  Accepted descriptors still use the dedicated fused worker.
The established broad contract accepts at least two-thirds coverage.  A
second, representation-based tier accepts at least one-third coverage only
when reported factor work is at least 16 times the numeric L+U entry count.
This is a reusable arithmetic-intensity condition rather than a matrix name,
dimension, or benchmark-family rule.  Allocation failures remain retryable,
whereas structural accept/reject decisions are cached.

Five counterbalanced entrywise-H100 pairs on each CCD put the new/old median
ratio at 0.869 and 0.855 for `add20`, 0.860 and 0.869 for `memplus`, and 0.970
and 0.966 for `cell2`.  Four- and two-thread checks retained the improvement
where the affected executor was selected.  `ACTIVSg70K` was the crossover
control: its 39.6-percent run coverage but only 14.9 work units per factor
entry fails the dense tier and preserves EGraph.  A 14-matrix remaining-gap
sweep completed 14/14 and had a 0.997 paired median geometric ratio, so the
targeted gains did not trade away aggregate performance.

Six boundary matrices completed 100 independently verified 0.1-percent
entrywise updates each.  The largest relative residual was 2.56e-13.  The
release and ASan/UBSan suites passed, as did explicit sanitized changed-value
loops for the rejected, accepted, and dense-tier row paths.

## Fused-run synchronization compression (2026-08-16)

Sampling the accepted short-supernode executor on `cell2`, `sts4098`, and
`thermal` put the padded per-row dependency scoreboard ahead of the fused
arithmetic in recurring samples.  Each worker executes its assigned schedule
slots in increasing order.  Consequently, when one fused run contains
several dependencies owned by the same remote worker, acquiring that owner's
latest required slot also makes every earlier slot from that owner visible.
The retained implementation precomputes one byte of wait bits per L-stream
position and polls only the last required dependency from each remote owner.
It preserves the padded, no-false-sharing scoreboard and changes neither the
factor arithmetic nor row publication order.  Allocation or representation
failure simply retains the complete wait sequence.

A complete scan of the 241 locally installed natural matrices found twelve
fused-run positives from eleven SuiteSparse groups.  Five counterbalanced
300-update pairs on CCD 0 put the compressed/full steady-refactor geometric
ratio at 0.985 and the complete H100 ratio at 0.991.  Eleven of twelve
per-matrix steady medians improved; `rajat03` was neutral at 1.001.  The
second CCD reproduced the result on eleven stable factors with ratios 0.983
and 0.991.  The pre-existing inaccurate `shyy41` route was noisy in both
arms, so a longer fifteen-pair run was used only for timing; its paired
medians were 0.979 for refactor and 0.986 for H100.  Two- and four-thread
tests on the original three factors improved steady refactor by 0.8--4.0
percent and H100 by 0.6--3.0 percent.

As a layout guard, the fourteen-matrix remaining-gap corpus also completed
three passes on both this revision and its parent with identical route
choices and no failures.  None accepted the fused-run executor (the sole
short-run candidate failed its generic coverage gate), so the optimization
does not add work to those unrelated paths.

The eleven numerically valid natural positives each completed 100 checked
entrywise updates at amplitudes 0.001, 0.01, and 0.1.  Another 600 checked
10-percent updates covered two- and four-thread schedules.  At amplitude 0.1
the compressed and complete-wait implementations produced exactly the same
maximum residual on every factor; the worst was 1.16e-10.  Release and
ASan/UBSan CTest passed, and sanitized 10-percent update loops covered all
twelve natural positives without a finding.  The rejected compact-scoreboard
prototype is not retained: despite identical residuals, it slowed `cell2`
from about 0.172 to 0.649 ms by increasing synchronization contention.

## Compile-time fused-width specialization (2026-08-16)

Profiles after wait compression still put recurring samples in the fused
short-supernode arithmetic.  Although the executor accepts only widths two
through eight, its single implementation received the width at run time.
That kept maximum-width local arrays and dynamic triangular loops in the hot
body.  The retained implementation dispatches once and instantiates the same
always-inlined arithmetic for each constant width from two through eight.
It changes neither the accepted factors nor their operation and reduction
order, and has no matrix policy, tuning threshold, environment selector, or
fallback implementation.  Every specialized width occurs in the natural
positive corpus.  The emitted hot body grows by about 1.7 KiB.

Five counterbalanced 300-update pairs over the twelve natural positives put
the specialized/dynamic steady-refactor geometric ratio at 0.986 on CCD 0.
The second CCD reproduced the gain at 0.992 over its eleven timing-stable
factors.  Two- and four-thread checks on `thermal`, `sts4098`, and `cell2`
gave steady-refactor ratios of 0.990 and 0.990.  A narrower width-two/eight
prototype was also a small win, but its 0.993 twelve-factor ratio left about
half of the full specialization's improvement unused even though widths
three through seven occur naturally.  Direct comparison of the final
executable with parent `e45bc5c` gave a 0.982 steady-refactor ratio over the
ten stable, non-tiny factors.  The fourteen-matrix remaining-gap control
suite completed 14/14 with identical routes and a neutral-to-positive 0.996
steady-refactor ratio; those factors do not select this executor.

The eleven numerically valid positives each completed 100 independently
verified entrywise updates at amplitudes 0.001, 0.01, and 0.1, for 3,300
checked updates.  Every run retained the row-refactor route.  The largest
per-update relative residual was 1.15e-10 and the largest terminal relative
residual was 3.79e-13.  Release and ASan/UBSan CTest passed.  Sanitized
10-percent entrywise loops added 220 independently checked updates over the
eleven valid positives and 20 memory-safety-only updates on the pre-existing
inaccurate `shyy41` route, without a finding.

Three adjacent redesigns were measured and removed.  An unrestricted
dependency/list schedule was two to five times slower and could change the
generic route selected by the cost model.  A dependency-level-constrained
schedule preserved the route but remained roughly 25--40 percent slower on
`sts4098` and `thermal`; reducing modeled waits did not compensate for lost
critical-path overlap.  Replacing the balanced triangular reduction with a
chain was also consistently slower, with representative steady-refactor
ratios of 1.004--1.021.  These results keep the optimization local to code
generation and preserve the established schedule and numerical ordering.
