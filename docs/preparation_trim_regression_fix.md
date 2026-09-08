# Regression investigation after preparation diagnostic cleanup

The requested cleanup was committed as `da66a15`. All removals remain:
preparation timers, the row-supernode census, manual supernode overrides,
the separator-private diagnostic scan, and lifecycle cycle maxima.

## Retained repair

Four `omp simd` annotations describe independent output-column iterations
in the dense SNB scratch product: initialization and accumulation, in both
the shared serial/small-cooperative helper and the large cooperative path.
Scratch output and input products occupy distinct live allocations; each
iteration writes a distinct output column. Producer-column accumulation
order is unchanged. No horizontal reduction or alignment assumption is added.

The serial SNB disassembly loses 75 instructions, including repeated runtime
alias checks. No executor, matrix selector, numerical threshold, allocation
policy, worker policy, or public interface is added. A stale comment about
the retired manual supernode override is also removed.

## Evidence and limits of the diagnosis

The cleanup regressions include steady work, not just initialization.
One-thread `onetone1` steady refactors rose from approximately 59.73 ms to
64.56 ms in the earlier cleanup recheck. Profiling confirms identical SNB
geometry (20,180 supernodes, 32,219 columns, approximately 32.1 MB of panels)
and identifies the dense scratch product as a major sampled region.

Pre-cleanup and cleanup serial SNB disassemblies have the same 1,476
instructions after normalizing addresses, apart from shifted private field
offsets. Code/allocation layout changed, but the precise hardware cause of
the original timing difference has not been proved. The retained repair
removes unnecessary generated work rather than reconstructing old padding.

TSOPF's eight-thread path uses the generic row worker, not SNB. Its remaining
loss has both refactor and changed-system solve components. In particular,
`solve_seconds_avg` measures the initial solve, **not** the changed-system
solves: the latter's steady average rose from about 153 to 164 microseconds
in the earlier CPU 0-7 recheck. Earlier exploratory summaries that inferred
unchanged solve performance from the initial-solve metric were incomplete.
The validation runner now records both metrics separately.

Steady-phase profiling locates a solve hotspot in the CSC residual scatter,
including its per-nonzero inverse row mapping. This identifies expensive
work, but does not by itself prove why the cleanup slowed it. Worker wait
samples also dominate process-wide profiles, which include idle workers.

## Rejected experiments

Completion/pool alignment, sparse-update SIMD and unrolling, scatter-helper
reuse in the numeric worker, stop-poll changes, ownership-token checks,
local dependency frontiers, prefix-state hoisting, and a broader existing
spin budget did not consistently recover TSOPF without tradeoffs. None is
retained.

Computing the residual in internal row order and unpermuting once per row
passed 168/168 screen launches, but improved lifecycle time inconsistently
while requiring another n-double workspace slice. Reusing the existing
scatter helper in that experiment also passed 168/168 without a convincing
overall recovery. Both changes were removed; the original four-vector
workspace, residual arithmetic, and accuracy contracts remain intact.

## Validation protocol

Artifacts: `build/prep-trim-repair-dBcoha/`. Only `retained-*` denotes the
final source candidate. `final-*`, `accepted-*`, and other earlier tags are
exploratory or interrupted runs, not validation of the retained source.

Three frozen executables are compared sequentially in rotating order:
`original` (`254ea42`, pre-cleanup), `slim` (`da66a15`), and `fix` (retained
working-tree changes). Build/dependency provenance must match the original
campaign. Metadata records source diff, binary/matrix hashes, and commands.
The campaign uses H100, entrywise changed values, automatic policies, and
verifies changed-system residuals against 1e-8. CPU sets are 8, 8-15, and
0-7. Each launch has a 90-second timeout and an 80 GiB address-space limit.
Profiled/debugger runs are excluded from performance comparisons.

Lifecycle comparisons are medians of paired percentage differences;
positive means slower. The broader retained campaign contains 483 launches
across 20 matrices and 58 configurations, with extra repetitions for the
previous regression cases. This is not a full paper or Xyce rerun.

## Retained results

Release CTest and ASan/UBSan CTest each pass 5/5. All 483 broad
campaign launches pass numerical validation; none timed out or failed.

Selected H100 lifecycle results (negative is faster):

| Matrix | Threads / CPUs | vs pre-cleanup | vs cleanup |
|---|---|---:|---:|
| transient | 1 / 8 | -1.92% | -4.36% |
| TSOPF_FS_b9_c6 | 1 / 8 | 0.77% | 0.15% |
| TSOPF_FS_b9_c6 | 8 / 8-15 | 1.90% | -0.17% |
| TSOPF_FS_b9_c6 | 8 / 0-7 | 2.74% | 0.64% |
| onetone1 | 1 / 8 | -5.23% | -11.35% |
| twotone | 1 / 8 | -6.00% | -11.24% |
| ASIC_320k | 1 / 8 | -3.65% | -6.26% |
| ASIC_680ks | 1 / 8 | -4.84% | -9.25% |

The larger serial regressions are recovered. **The eight-thread TSOPF
regression is not fixed**: it remains 1.90–2.74% slower than pre-cleanup.
The TSOPF comparisons use 12 triplets at one thread and 16 triplets per
eight-thread CPU set, so they are not merely the two-triplet control screen.

Other controls generally have two triplets; Raj1 and mac_econ_fwd500 have
one per CPU set. These are screening results, not confidence bounds or proof
of no regressions. Flagged controls receive a separate repeated check.

Retained executable SHA-256:
`ce4afc3951a3b293455c3f0ded55e45580de405b3a8280b85fa06eefc7bf13cb`.

This is a partial repair, not completion of the whole regression request.
All slimming remains, and all speculative changes have been removed.

### Flagged-control rechecks

`retained-controls-*` freezes the identical executable. All 123 additional
launches pass, for **606/606 valid retained-candidate comparison launches**
including the two baselines. No timed performance run overlaps compilation,
sanitizer testing, or another comparison runner.

| Matrix / configuration | Triplets | vs pre-cleanup | vs cleanup |
|---|---:|---:|---:|
| rajat03 / 8 threads, CPUs 8-15 | 12 | -0.30% | -0.09% |
| adder_dcop_01 / 8 threads, CPUs 8-15 | 12 | -0.12% | +0.06% |
| bcircuit / 1 thread | 8 | -0.03% | -0.65% |
| onetone1 / 8 threads, CPUs 0-7 | 6 | -0.86% | -0.57% |
| Raj1 / 8 threads, CPUs 8-15 | 3 | +0.72% | +1.15% |

The earlier control warnings do not reproduce at their original magnitude.
Raj1's modest remaining difference has only three triplets; it should not
be called either a demonstrated recovery or a statistically established
regression. TSOPF's repeated eight-thread loss remains the clear unresolved
case. Its exact low-level cause needs stronger evidence before further
worker or residual-loop changes are justified.
