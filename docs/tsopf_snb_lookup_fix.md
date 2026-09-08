# TSOPF first-refactor SNB setup repair

The cleanup is committed separately as `ba52b9e`. This repair retains that
cleanup and the preceding removals in `27ef224`.

## Attribution

A same-session comparison of the frozen `d0ae5a6`, `27ef224`, and `ba52b9e`
binaries reproduced the one-thread first-refactor increase, approximately
72ms to 80ms. `KLS_TRACE_PREP_CONSULT` localized it to the SNB preparation and
acceptance call: approximately 63ms before, 71–72ms after. Other preparation
phases and the lean-engine consultation did not account for the difference.

The built-in sampling profiler and address-to-source resolution locate the
main preparation hotspot at the producer-edge scan in `kls_snb_prepare`'s
packed-U writeback construction. The same linear search was duplicated in
input scatter construction. For every entry it restarted at edge zero, even
though construction already sorts the producer edges. The numerical SNB trial
itself took about 1ms and was rejected in these traced runs. Thus this was
setup for an engine trial, not a slower recurring linear solve.

The old and slimmed binaries execute the same instruction sequence in that
hot scan, at different addresses, and produce the same 4,841 supernodes and
205 schedule levels. The original slowdown is not evidence of a removed
numeric algorithm. A precise hardware/cache/layout explanation for the
additional scan time is not established: hardware performance counters are
unavailable in this environment. The repair removes the expensive search
pattern rather than relying on padding or instruction placement.

## Change

One shared helper now binary-searches the producer and its subcolumn list.
The producer ordering follows from the sorted U-row union and the consecutive
column partition into supernodes. Missing entries still return failure, and
both callers retain their preparation-failure path. Numeric arithmetic, panel
contents, engine acceptance gates, threading, and public ABI are unchanged.
There are no new tuning thresholds, matrix-name checks, or special cases.

The existing deferred-sort smoke fixture now uses four-column panels, giving
16 producers per dense block, and checks a transpose solve as well as normal
changed-value solve accuracy. Release and ASan/UBSan CTest pass 5/5 each;
edited sources build without warnings.

## Benchmark protocol

Artifacts are in `build/tsopf-first-repair-RaQMmj/`. Frozen binaries are run
sequentially with rotating order, no tracing during timing, H100 entrywise
updates (amplitude 0.001), verification after every refactor at `1e-8`, an
80GiB address-space limit, and a 90-second launch timeout. Build provenance
matches exactly between cleanup and repair. Against the original binary,
compiler settings, linked targets, submodules, external objects, and runtime
libraries match; the intentional bundled-KLU changes in `27ef224` are recorded
as different dependency trees, not silently treated as identical sources.

Thirty one-thread triples and twenty triples on each eight-thread cache domain
compare the original, cleanup, and repaired versions. Fifteen other matrices
receive four cleanup/repair pairs on each of those three CPU configurations.
This is bounded matrix-lifecycle validation, not a full paper or Xyce rerun.

## TSOPF results

| One-thread median | Original | Cleanup | Repaired |
| --- | ---: | ---: | ---: |
| Initial factor | 10.06ms | 10.09ms | 10.07ms |
| First refactor | 72.62ms | 80.17ms | 24.30ms |
| Steady refactor | 0.987ms | 0.979ms | 1.015ms |
| H100 lifecycle | 205.02ms | 211.55ms | 159.46ms |

Median paired lifecycle changes are **−24.66% versus cleanup** and **−22.21%
versus original**. The small steady-refactor increase is included in the
lifecycle result, not hidden by reporting only the first-call improvement.
Eight-thread paired lifecycle changes versus cleanup are −0.20% and −0.03%;
versus original they are +0.55% and +1.18%. All 210 TSOPF launches pass.

## Controls and rechecks

All **570/570 main-trial launches pass**, with maximum relative residual
`5.45371652e-9`. No eight-thread TSOPF steady-refactor average exceeds 0.6ms
in any version (fixed maximum 0.493ms), preserving the earlier affinity repair.

The following are median paired lifecycle changes versus cleanup, four pairs
per cell; positive means slower. Small samples are not proof of exact parity.

| Matrix | 1 thread | 8 threads, 32MiB | 8 threads, 96MiB |
| --- | ---: | ---: | ---: |
| TSOPF_FS_b9_c1 | -2.15% | -1.24% | -2.91% |
| TSOPF_RS_b9_c6 | -0.21% | -2.79% | -0.79% |
| onetone1 | -3.96% | +1.60% | +1.77% |
| twotone | -1.38% | -1.30% | +1.27% |
| rajat25 | -8.78% | +1.01% | +0.93% |
| adder_dcop_01 | -0.69% | +2.54% | -2.54% |
| 1138_bus | +0.63% | +2.07% | -5.88% |
| bips98_1142 | +0.86% | +0.66% | -0.34% |
| ASIC_100ks | -5.03% | +1.11% | +2.27% |
| ASIC_320k | -3.92% | +0.11% | +0.86% |
| ASIC_680ks | -3.87% | -0.13% | -0.19% |
| ckt11752_dc_1 | -3.86% | +0.04% | +0.63% |
| circuit_4 | -30.60% | +0.05% | -0.23% |
| transient | -18.27% | -0.54% | +0.59% |
| rajat20 | -7.01% | +0.21% | +1.60% |

Every control initially slower by more than 2% received 16 fresh pairs:
adder_dcop_01 (32MiB) **+0.21%**, 1138_bus (32MiB) **-0.47%**, and
ASIC_100ks (96MiB) **+1.38%**. All **96/96 recheck launches pass**, for
**666/666 total**. The two short-case slowdowns did not repeat; a smaller
ASIC_100ks increase remains and is not claimed to be eliminated. The large
TSOPF setup improvement is therefore not a claim of universal speedup.

A separate traced smoke run confirms the modified fixture actually prepares
4,096 SNB supernodes over 16,384 columns, executes a successful SNB trial, and
passes its normal and transpose checks. The Release binary matches the frozen
repair binary, SHA256
`fa50ac1fdfa9b6d67c32ae7c3e35694d0d48029f6bd7cfa1587488ddf19b2268`.
