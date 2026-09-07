# TSOPF worker-placement repair

## Cause and change

The full slimming campaign measured a +23.02% median paired lifecycle change
for TSOPF_FS_b9_c6 on the 32MiB cache domain, but only +0.64% on 96MiB.
The difference was concentrated in steady refactoring, with individual slow
launches in both the old and slimmed binaries. It was not extra initial-factor
or preparation work.

The shared affinity planner treated a cpuset containing enough physical cores
in one homogeneous LLC as equivalent to explicit one-worker-per-core placement.
It returned without pinning the pool. A shared allowed set does not establish
that placement: workers can migrate or be delayed while dependency consumers
spin. The shortcut existed in both versions.

Read-only interventions on the frozen campaign binaries used the existing
explicit-placement policy on the same eight physical cores. At 100 systems,
the baseline/current unpinned cases had 8/40 and 4/40 launches above 0.6ms per
steady refactor; explicitly placed cases had none. A 1,000-system test confirmed
the pattern. All 240 diagnostic launches passed. Historical kernel scheduling
events were not recorded, so the exact co-location/migration sequence remains
unproven. Diagnostic evidence is in `build/tsopf-cause-K5E8nY/`.

The fix deletes the homogeneous-cpuset early return in
`kls_compact_llc_affinity_plan`. Both affected pools now use the existing
topology-aware placement policy even when the caller has constrained execution
to one cache domain. It retains the repeated-update admission rule, affinity
opt-out, topology/capacity checks, partial-failure fallback, and caller-affinity
restoration. No matrix-specific rule, numerical threshold, or padding is added;
all unrelated slimming remains. Production source decreases by 45 net lines.

## Correctness checks

Release and ASan/UBSan CTest each pass 5/5. The smoke harness now checks that
calls passing through `require_ok` leave the caller's CPU affinity unchanged
on Linux where affinity can be queried. The parallel transpose fixture runs
both with and without repeated-workload hints to exercise placement.
Builds have no new warnings and
`git diff --check` passes.

Additional constrained smoke runs on CPUs 8-15 and 0-7 pass and emit explicit
placement plans containing distinct CPUs wholly inside the allowed set. The
same test on CPUs 8-15 with `KLS_DISABLE_COMPACT_LLC_AFFINITY=1` passes without
any placement plan. Caller-affinity checks remain enabled in all three runs.

## Performance protocol

The before binary is committed `836b9d3`; after is that revision plus this fix.
The harness checks matching build provenance across all linked targets, freezes
the binaries, alternates paired launches, fixes CPU affinity, and verifies every
refactor at relative residual limit `1e-8`. Auxiliary library pools use one
thread. Entrywise updates have amplitude 0.001. No failed launches or timing
outliers are dropped.

Artifacts: `build/affinity-fix-hm5nbd/` (100-system targeted comparisons) and
`build/affinity-long-kxxcVI/` (1,000-system TSOPF comparisons). The completed
full campaign is preserved unchanged; these are separate repair checks.

In the 100-system TSOPF checks, 40 pairs per eight-thread domain measured
-0.74% lifecycle change on 32MiB and +0.45% on 96MiB. Steady-refactor outliers
above 0.6ms fell from 7/40 to 0/40 on 32MiB and from 3/40 to 0/40 on 96MiB.
The 0.6ms cutoff is descriptive only, not a solver policy. Fixed steady maxima
were 0.429ms and 0.441ms, versus 1.696ms and 1.688ms before.

All **328 targeted launches passed**. One-thread TSOPF changed -1.17% in
20 pairs. Eight control matrices (four pairs per cache domain) ranged from
-2.03% to +2.15%; they include TSOPF_FS_b9_c1, TSOPF_RS_b9_c6, onetone1,
twotone, rajat25, adder_dcop_01, 1138_bus, and bips98_1142. These small samples
do not establish exact performance equality, but no large short-lifecycle
penalty appeared.

All **80 long-lifecycle launches passed** (20 pairs per cache domain).
At 1,000 systems, paired lifecycle changes were -10.54% on 32MiB and -9.50%
on 96MiB. Median steady refactor fell from 0.486ms to 0.415ms and from 0.479ms
to 0.419ms. Fixed steady maxima were below 0.440ms on both domains, with no
launches above 0.6ms. The unchanged-before binary had four such launches on
32MiB and one on 96MiB.

In total **408/408 repair benchmark launches passed**. The repaired binary
has not yet undergone another full-corpus campaign; these results establish
the targeted placement repair, not a full-corpus no-regression guarantee.
