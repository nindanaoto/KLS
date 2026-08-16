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
body.  The retained implementation dispatches once and presents each width
from two through eight as a constant to the same always-inlined arithmetic.
This specialized the trailing target loops, although a later assembly audit
found that GCC still merged the triangular prefix back into a dynamic loop.
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

## Dedicated fused-prefix unrolling (2026-08-16)

A fresh fine-grained profile and disassembly of `aad53a8` showed that GCC had
tail-merged the supposedly constant short triangular solves into one
runtime-bound loop.  The retained change emits one function for each natural
width two through eight and explicitly unrolls only the bounded triangular
prefix.  The potentially long trailing update stays looped and keeps its
balanced reduction tree.  The width switch is inlined into the worker, so
each run still pays one dispatch and one call.  There are no new selectors,
shape gates, matrix rules, or numerical operations.  The seven hot functions
total 3,533 bytes versus 3,620 bytes for the merged body; total executable
text grows by 560 bytes because of the inlined dispatch and alignment.

Five counterbalanced 300-update pairs over all twelve natural positives on
CCD 0 improved every steady-refactor median.  The dedicated/merged geometric
ratios were 0.948 for steady refactor and 0.973 for complete H100.  The other
CCD reproduced 0.946 and 0.971 over the eleven timing-stable factors.  On
`thermal`, `sts4098`, and `cell2`, two-thread ratios were 0.935 and 0.952,
and four-thread ratios were 0.929 and 0.956.  A simpler single-function
unrolled implementation was measured and removed: it was about 1.010 for
both steady refactor and H100 against the dedicated functions on all three
representatives, consistent with shared-prologue and register-pressure cost.

The fourteen-matrix remaining-gap control completed 14/14 with identical
route choices.  Its ratios were 1.001 for steady refactor and 1.002 for H100,
which is neutral short-run noise on executors that do not call these
functions.  The eleven numerically valid positives then completed 100
independently verified entrywise updates at amplitudes 0.001, 0.01, and 0.1:
3,300 checks, all on the row route, with a 1.15e-10 largest per-update and
3.79e-13 largest terminal relative residual.  These exactly match the prior
implementation's maxima.

Release and ASan/UBSan CTest passed.  Forcing the lean executor under the
sanitizers covered 220 checked 10-percent updates over the eleven valid
positives and 20 memory-safety-only updates on `shyy41`, without a finding.
A separate forced-cooperative-row run exposed a dense-help stack-context
lifetime race; the unmodified `aad53a8` parent reproduced it, so it is a
pre-existing defect outside this fused executor and is not attributed to the
unrolling change.

## Fused width-eight AVX-512 target update (2026-08-16)

The post-unrolling profile left the fused producer arithmetic and genuine
producer waits as the recurring costs.  Polling the dependency generation
less often was rejected: batches of eight polls measured 1.003 for steady
refactor and 0.999 for H100, while batches of four measured 1.005 and 1.002.
Those neutral-to-negative results indicate that the samples were observing
producer latency rather than excessive polling bookkeeping.

The retained change vectorizes the long trailing update of a width-eight
fused run.  This case is materially different from the previously rejected
per-dependency gather/scatter experiment: eight producer rows share one
sorted target stream, so one indexed gather and scatter are amortized over
eight contiguous coefficient streams.  Each vector lane uses the same FMA
placement and balanced reduction tree as the scalar width-eight kernel, and
the scalar kernel handles the remainder.

Activation depends only on the representation already built by the generic
fused executor: the producer width must be exactly eight and the common tail
must contain at least 32 targets, or four AVX-512 vectors.  The run census
records whether any such work exists, and the dispatching thread resolves the
existing AVX-512 feature gate once per refactor.  There is no matrix, size,
corpus, ordering, or policy rule and no new environment selector.
`KLS_AVX512_SCATTER=0` retains the existing diagnostic way to force the
scalar fallback.

The AVX-512 width-eight step and its worker copy live in a separate
`.text.kls_avx512` section.  This matters for generality as well as code
organization: the scalar width-two through width-eight functions retain the
same addresses and sizes as parent `10d167c`, and unrelated generic workers
retain their layout.  A baseline `-march=x86-64` GCC syntax build also passes,
so compilation does not depend on the build host exposing AVX-512.

Five counterbalanced 300-update pairs over all twelve natural positives gave
a 0.958 steady-refactor geometric ratio and 0.977 H100 ratio over the eleven
timing-stable factors.  Ten of eleven steady medians improved; the sole small
loss, `ACTIVSg10K` at 1.001, disappeared in a longer same-executable
AVX-on/off run (0.9997 steady and 0.992 H100 over fifteen pairs).  Two-thread
and four-thread checks on `thermal`, `sts4098`, and `cell2` gave geometric
ratios of 0.942/0.955 and 0.941/0.959 respectively.  The inaccurate `shyy41`
route remains excluded from aggregate claims.

The fourteen-matrix remaining-gap control completed 14/14 with identical
routes and ratios of 0.991 for both steady refactor and H100.  Disabling
AVX-512 over six representatives gave 0.994 and 0.997 against the parent,
showing that the scalar/non-feature path did not inherit a layout regression.
Broader all-width and 64-target-cutoff variants were measured and removed:
the former regressed `ACTIVSg10K`, while the latter was slower over the stable
corpus than the four-vector crossover.

The final release and ASan/UBSan CTest suites passed all four tests.  Forced
sanitized execution added 220 independently checked 10-percent entrywise
updates over the eleven numerically valid positives; every factor retained
the row route, the largest relative residual was 1.30e-11, and no sanitizer
finding occurred.  The emitted AVX section contains the intended
`vgatherqpd` and `vscatterqpd` instructions; non-AVX and non-x86 builds retain
the scalar implementation.

## Compact-index fused row worker (2026-08-16)

The next representation audit found that the fused worker still read every L
dependency and U target from 64-bit column streams even though the row engine
already had a checked 16-bit mirror builder.  The retained specialization
consumes those compact column streams throughout the fused worker.  Scalar
dependency waits, ordinary sparse updates, fused trailing updates, and U-row
publication all widen the stored `uint16_t` value only when it is used.  The
width-eight AVX-512 kernel widens eight targets from one 16-byte load with
`vpmovzxwd`, then uses `vgatherdpd` and `vscatterdpd`; the arithmetic and
balanced reduction tree are unchanged.

Admission is representation- and lifecycle-based.  The existing fused worker
must be structurally capable, the installed factor must have at most 65,536
rows, the caller must request more than one thread, and it must advertise at
least 64 expected refactorizations.  The last condition is deliberately more
conservative than the ordinary recurring-workload floor: building both compact
column streams is one O(L+U) pass.  Measured H16/H24/H32 geometric ratios were
1.005/1.002/0.999 when that setup was charged only to the candidate, while H64
and H100 were 0.992 and 0.989.  One-thread measurements were neutral apart from
the same setup cost, so automatic construction is confined to parallel
lifecycles.  Allocation or range-check failure retains the established 64-bit
worker, and `KLS_DISABLE_LEAN_I16_INDICES=1` remains the existing diagnostic
fallback.

The decisive comparison used one final executable on both sides, with compact
indices enabled normally versus disabled by that diagnostic control.  Ten
counterbalanced 100-update pairs over all twelve natural fused positives gave
a 0.975 steady-refactor and 0.990 complete-H100 geometric ratio over the eleven
numerically valid factors.  No valid matrix regressed by two percent in H100.
An independent parent/current comparison that charged compact construction
only to current gave 0.966 and 0.989.  Two- and four-thread parent/current
checks gave 0.953/0.976 and 0.962/0.982 for steady refactor/H100.  Disabling
AVX-512 still gave a 0.974 steady-refactor ratio and a neutral 0.998 H100,
showing that compact scalar stream traffic is independently useful.

The fourteen-matrix remaining-gap control retained identical executor routes
and completed 14/14; its same-binary compact-enabled/disabled ratios were 1.000
for steady refactor and 0.993 for H100.  Existing scalar fused functions and
the ordinary fused worker retain their parent addresses and exact code sizes;
the existing AVX-512 functions retain exact sizes.  New compact scalar and
vector functions live in `.text.kls_i16_snode`, so unrelated hot code is not
displaced.

Release and ASan/UBSan CTest passed all four tests.  Forced compact execution
under the sanitizers completed 440 independently checked 10-percent entrywise
updates over the eleven valid natural factors, all on the row-refactor route,
with a largest per-update relative residual of `1.29423314e-11` and no finding.
GCC `-march=x86-64` and Clang syntax builds also pass (apart from the existing
Clang uninitialized-thread warning), so the translation unit does not require
the build host's vector ISA.

## Frontier-predecoded fused synchronization (2026-08-16)

Sampling the compact fused worker after the column-stream change put its
largest remaining bucket in dependency polling.  The fused executor still
read one 64-byte padded completion record for each remote dependency even
though the retained static schedule already had a cache-line-partitioned
grouped scoreboard.  A direct substitution was measured and rejected: the
extra row-to-token lookup made H100 0.6 percent slower.  The retained design
instead predecodes the complete synchronization frontier once.

Each consumer worker executes its interleaved row subsequence in topological
order, and each producer executes its own subsequence in the same order.  An
acquire of producer slot `s` therefore covers all earlier slots from that
producer for every later row of the consumer.  Pattern setup records a grouped
slot only when a dependency advances that per-consumer/per-producer frontier.
Within a fused run it retains only the latest advancing slot per producer.
The same pass marks exactly the producer rows that some consumer will acquire;
the numeric worker omits every other completion publication.  Release/acquire
ordering and the numeric arithmetic are unchanged.  For example, `thermal`
needs 6,257 frontier acquires and 2,033 publications for 145,251 L entries and
3,456 rows.

This is a schedule-representation capability, not an input classifier.  It
uses the existing generic short-supernode eligibility and repeated-update
contract, supports two through eight workers, validates every row, owner, and
slot before publishing the stream, and falls back to the established generic
row worker if grouped-scoreboard or frontier allocation fails.  The new stream
costs four bytes per L entry.  It is stored at the cold tail of the private
solver object so existing hot-field offsets remain fixed, and there is no new
environment selector or matrix/size/corpus rule.

Fifteen counterbalanced 500-update parent/current pairs on the five higher-work
fused factors gave geometric ratios of 0.980 for steady refactor and 0.989 for
complete H100.  Four of five steady medians improved; `utm3060` was neutral at
1.001 H100.  A five-pass twelve-factor screen gave 0.991 H100 over the eleven
numerically valid factors with no two-percent regression.  With every tuning
environment variable removed, the ordinary public auto policy gave 0.993
H100 over the same eleven factors, with unchanged routes and no regression
above one percent.

The one-time frontier construction is already amortized at the generic
recurring-workload floor: fifteen H16 pairs on the five focus factors gave
0.968 steady-refactor and 0.992 complete-H16 ratios.  Two- and four-worker
checks on `cell2`, `sts4098`, and `thermal` gave H100 ratios of 0.969 and 0.987.
Disabling AVX-512 over the five focus factors still gave 0.989 H100.  The
fourteen-matrix remaining-gap control gave 0.997 H100 over its eleven valid
rows; a thirty-pair adjudication of its apparent small-factor outlier,
`circuit204`, was neutral at 1.004 H100 and 1.004 steady refactor.

The ordinary compact-match and generic row workers retain their parent
addresses and exact code sizes; only the four intended fused worker copies
change.  Release and ASan/UBSan CTest pass all four tests.  Forced sanitized
execution completed 440 independently checked 10-percent entrywise updates
over all eleven valid natural factors, all on `row_refactor`, with a largest
per-update relative residual of `1.15056277e-10` and no finding.  GCC
`-march=x86-64` and Clang syntax checks pass, apart from the existing Clang
uninitialized-thread warning.
