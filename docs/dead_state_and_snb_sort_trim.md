# Dead state, diagnostic paths, and SNB union-sort cleanup

Relative to `e01e3ca`, removed:

- Six private fields with no live nonzero state: `lean_pending`,
  `lean_probe_min`, `compact_partial_diagonal_column_fringe`,
  `egraph_tight_tol_state`, `prestatic_dense_spiked_match`, and
  `large_sparse_amf3_path`, including resets and constant condition terms.
- The historical opt-in `KLS_VETO_EGRAPH_TIGHT_TOL` veto. The numeric-replacement
  cache invalidation and current numerical recovery checks remain.
- `KLS_ENABLE_GENERIC_EARLY_COORDINATE_FREE_ND`, which launched ND before the
  coordinate frame was settled. The normal coordinated ND handoff remains.
- Seven private atomic light-supernode census counters, their update sites,
  their verbose reporting, and the now-unused shared trace flag. The numeric
  kernel, structural census, public statistics, and phase timers remain.

Corrected README documentation: the BTF scalar-run executor has automatic
eligibility as well as an explicit override; it is not universally off by
default. Its implementation and selection policy are unchanged.

## Regression discovered and repaired during cleanup

The first cleanup trial in `build/dead-state-trim-s1HjtE/` passed all 376
launches, including Raj1 and mac_econ_fwd500 checks on both cache domains.
However, one-thread TSOPF_FS_b9_c6 lifecycle time regressed by 4.78%:
first refactor rose from 24.33ms to 31.83ms, with steady refactors and solves
essentially unchanged. This intermediate cleanup was not accepted as final.

Preparation tracing attributed the extra time to SNB setup (approximately
15.6ms versus 23.2ms). Sampling and address-to-source resolution located the
hotspot at the custom insertion sort of U-row unions in `kls_snb_prepare`.
Its source algorithm had not changed; the precise machine/layout mechanism
behind the slowdown is not established.

Replaced that quadratic insertion loop with `qsort` and KLS's existing
integer comparator, preserving the already-sorted fast path. The union is
deduplicated before sorting; the resulting producer order is unchanged.
This is an algorithmic simplification, not a matrix-specific exception,
padding adjustment, or new numeric threshold. Together the changes remove
**118 net source lines**, excluding tests and documentation.

The deferred-sort smoke fixture retains its dense case and adds interleaved
even/odd subblocks. With explicit test-only amalgamation settings, successive
panel columns contribute different sorted U lists whose union is not sorted
(for example, [0,2] followed by [1,3]). Both cases check changed-value normal
and transpose solves, as well as rebuilding solve caches after sorting.

## Final validation protocol

Final artifacts are in `build/dead-state-final-xKm5yA/`. Frozen baseline and
after binaries have matching compiler options, linked targets, dependency
revisions, external objects, and runtime libraries. Runs are sequential with
alternating order, no tracing during timing, H100 entrywise updates of amplitude
0.001, and residual verification after every refactor at `1e-8`. Each launch
has an 80GiB address-space limit and a 90-second timeout.

TSOPF_FS_b9_c6 receives 20 pairs at one thread and on each eight-thread domain.
Fifteen controls receive four pairs on each eight-thread domain; rajat25,
circuit_4, and transient additionally receive two one-thread pairs. Raj1 and
mac_econ_fwd500 each receive one final pair per eight-thread domain, in addition
to their initial-cleanup checks. This is bounded matrix-lifecycle validation,
not a full paper-corpus or Xyce rerun.

## Final results and remaining caveat

Final Release and ASan/UBSan CTest each pass **5/5**, including both deferred-sort
fixture variants. A traced smoke run confirms that both variants build 4,096
supernodes and execute successful SNB trials; the interleaved fixture uses
test-only zeta 2 so adjacent parity classes actually merge into panels. Builds
are warning-free. Retired-switch smoke checks pass, and the removed runtime
census no longer appears in verbose output. The final fixture adjustment does
not change the benchmark binary; `final.diff` records the complete final patch.

All **380/380 final-pass launches pass**. TSOPF_FS_b9_c6 median paired lifecycle
changes versus `e01e3ca` are **-4.77%** at one thread, **-1.74%** on the 32MiB
eight-thread domain, and **-2.32%** on the 96MiB domain. One-thread first refactor
falls from **24.31ms to 16.65ms**, while steady refactor remains about 0.98ms.
The previous worker-affinity fix is retained.

The additional two-pair one-thread controls show rajat25 **-14.67%**, circuit_4
**-15.39%**, and transient **+6.80%**. These are small samples; the latter
regression was investigated further rather than dismissed. The only other
initial slowdown above 2% was TSOPF_RS_b9_c6 at 32MiB (**+2.41%**).

Sixteen fresh alternating pairs for each flagged configuration pass **64/64**.
TSOPF_RS_b9_c6 measures **+1.74%**, while transient measures **+6.60%**: the
transient regression repeats and remains unresolved. It is concentrated in
steady SNB execution, not sorting/setup. Separate profiled runs of the baseline,
initial cleanup, and qsort cleanup measure approximately 29.1ms, 33.2ms, and
33.2ms per steady refactor, respectively, all using SNB. Thus the slowdown
predates the qsort replacement. Sampling points to unchanged
`kls_snb_process_target` update arithmetic; a precise layout/hardware cause is
not established. The removals are retained as requested, not represented as
performance-neutral. No padding or matrix-specific compensating branch was
introduced to conceal the regression.

Total bounded validation across the intermediate cleanup, final pass, and
rechecks is **820/820 valid launches**. Raw observations, phase medians,
provenance, and profiles are retained in the two artifact directories. The
final Release benchmark matches its frozen after binary, SHA256
`300747e2a884cd89e128742c0f61811a86739aa007c40bcead6630d35abc31b8`.
