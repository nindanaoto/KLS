# Paper Ideas Audit

This note records which ideas from the reference papers are present in KLS and
which remain open. The intent is to keep KLS development focused on general
solver algorithms instead of tuning individual benchmark matrices.

## Current Conclusion

KLS has **not** implemented every paper idea that is still worth trying. It has
implemented the ideas that can be layered around the current KLU-derived
Gilbert-Peierls kernel: BTF, AMD/COLAMD/METIS ordering policy, explicit SCOTCH
ordering, CAMD refinement, auto scaling policy, pivot-checked reuse, static
row-pivoting trials with dual-potential matching-derived equilibration, and
BTF-block parallel refactorization with a solver-owned worker pool. KLS now
also keeps precomputed refactor scatter metadata for unscaled serial repeated
refactors, covering both single-block and serial BTF cases, records an exact
no-pivot EGraph level schedule from the numeric U pattern, and consumes that
schedule in a guarded cluster/pipeline intra-block refactor path for large
single-block matrices, including KLU row-scaled cases whose scale factors can
be recomputed and permuted safely, and inside large dominant BTF blocks whose
off-block entries can be refreshed from the retained map. The barriered EGraph
cluster levels now use FLOP-estimated per-thread slices instead of equal column
slices, which is a retained piece of the CKTSO/SubtreeLU load-balance idea.
Its static-pivot
preprocessing has a cheap exact sparse maximum-log-product assignment path for
small candidates and can improve medium row matchings with bounded alternating
cycles beyond the pair-swap pass. When optional SPRAL support is enabled, KLS
can also use BSD-licensed Hungarian matching/scaling as a pre-factor
MC64-adjacent candidate for large weak-diagonal dominant-block matrices. Fast
factorization can now repair an unsafe unscaled BTF diagonal block by
restarting that block with pivoting and then retrying the checked no-pivot
factorization.

The remaining CKTSO gap is large enough that it should be treated as a missing
major algorithm, not an ordering-package tuning problem. On the selected large
`pre2` case, CKTSO completed the factor-only comparison inside the 120s cap
with about 28.4s cycle time, while KLS timed out under AMD/BTF, METIS/no-BTF,
SCOTCH/BTF, and auto policy variants. The CKTSO ordering supplement reports
`pre2` operation counts for CKTSO nested dissection and METIS in the same
range, so the current evidence points more strongly at MC64-quality
matching/scaling and CKTSO's KLS-owned numeric/scheduling machinery than at
another separator package alone.

The remaining worthwhile ideas are not per-matrix tuning knobs. They require
new KLS-owned symbolic/numeric machinery:

- Production-scale MC64-equivalent maximum-weight matching with dual
  row/column scaling. KLS can now build BSD-licensed SPRAL Hungarian/auction
  matching from a pinned submodule or link to a system SPRAL install and uses
  it for bounded pre-factor large weak-diagonal trials, but this is still not
  a full production MC64-equivalent preprocessing stage. HSL MC64 and
  non-redistributable MC64 copies are out of scope for vendoring.
- A full intra-block parallel factor/refactor scheduler that consumes retained
  EGraph/ETree or separator-tree metadata with pivoting tail restart, beyond
  the current guarded no-pivot EGraph cluster/pipeline refactor.
- Full CKTSO-style fast factorization tail restart after a failed pivot check,
  beyond the current unscaled block-local repair.
- SubtreeLU-style private/pipeline scheduling from a retained separator tree.
- A structure-adaptive triangular solve built on LU storage that exposes cheap
  row/segment access.

Several smaller dispatch experiments were tried and rejected because they helped
some benchmark cases while regressing others. Those are documented below so the
project does not drift toward benchmark-name-specific heuristics.

## Paper-by-Paper Coverage

| Reference | Implemented in KLS | Partial or open coverage |
| --- | --- | --- |
| Algorithm 907 / KLU | BTF preprocessing, fill-reducing ordering, row scaling modes, Gilbert-Peierls factorization with partial pivoting, no-pivot refactorization, and block back substitution are present through the vendored KLU-derived kernel. KLS adds automatic policy selection, serial refactor scatter metadata, and exact EGraph level metadata around these pieces. | KLS still inherits KLU's fundamentally sequential intra-block numeric kernel. |
| NICSLU | AMD-style ordering, optional static-pivoting preprocessing, optional SPRAL Hungarian/scaling trials, and the idea that parallel kernels should be selected by general structural/numeric evidence are represented in KLS policies. KLS now records exact no-pivot EGraph levels from the numeric U pattern and uses them in a guarded large single-block or dominant-BTF-block cluster/pipeline refactor path, including KLU row-scaled cases and work-estimated cluster-level thread slices. | Full production MC64 matching/scaling is not implemented. Static symbolic R1/R2 performance modeling, ETree/EScheduler-guided intra-block factorization, and full pivoting-aware ETree scheduling are not implemented. Earlier broader EGraph prototypes were rejected because they were not general wins on the current kernel/storage. |
| CKTSO | METIS nested-dissection ordering, explicit SCOTCH ordering for experiments, constrained-minimum-degree-style CAMD refinement, combined ordering selection, pivot-checked fast factorization, unscaled block-local restart after a failed fast-factor pivot check, guarded EGraph cluster/pipeline no-pivot refactors for single and dominant BTF blocks, work-balanced cluster-level refactor slices, and dual-potential plus optional SPRAL matching-derived equilibration trials are implemented in KLS at the KLU-wrapper layer. | CKTSO's maximum-weight matching with dual scaling is still only partially approximated because KLS does not have an always-on production MC64-equivalent weighted assignment stage. Pipelined ETree-descendant tail restart with pivoting after a failed pivot check and structure-adaptive triangular solve are not implemented. Scaled or otherwise unsupported fast-factor failures still fall back to full pivoting factorization. |
| SubtreeLU | KLS vendors reproducible METIS/GKlib submodules and uses METIS plus CAMD refinement, which overlaps with SubtreeLU's nested-dissection and constrained-ordering motivation. | KLS does not retain a separator tree, collapse/partition it into private and pipeline task queues, constrain pivot search within separator-tree subdomains, perform FLOP-balanced refactor queue generation, or use SubtreeLU-style supernodes/BLAS updates. |

This means KLS has implemented or prototyped the ideas that can be layered
around the current KLU-derived data structures. It has **not** implemented all
paper ideas that are still worth trying. The remaining items are general solver
design work, not benchmark-specific tuning.

## Implemented

- KLU Algorithm 907 baseline: BTF preprocessing, AMD/COLAMD ordering,
  KLU-style row scaling modes, Gilbert-Peierls factorization with partial
  pivoting, no-pivot refactorization, block back substitution, diagnostics, and
  pivot tolerance controls are available through the vendored Trilinos
  SuiteSparse-derived sources.
- Nested-dissection ordering option: KLS vendors METIS/GKlib and SCOTCH as
  pinned submodules by default, and can optionally use compatible system METIS
  or SCOTCH installs.
- Combined ordering policy: KLS auto mode can choose AMD, COLAMD, or METIS,
  retry no-BTF symbolic analysis for structural cases where BTF is not useful,
  and promote expensive numeric factorizations to METIS when actual fill/flop
  evidence is better. The no-BTF retry covers both single-block BTF analyses
  and large dominant-block BTF analyses when the stripped fringe is small and
  the no-BTF symbolic estimate is materially better, matching KLU's warning
  that BTF can occasionally increase factor fill. Low-work dominant-BTF cases
  are excluded from this retry because keeping BTF is cheaper when the symbolic
  work estimate is already low. The METIS promotion gate also
  covers small BTF-dominant matrices when the first factorization shows many
  off-diagonal pivots and high actual fill/flop growth, and preserves the
  selected BTF/no-BTF mode when comparing a METIS promotion candidate. KLS also
  starts directly with METIS for narrow large-diagonal structural classes from
  the paper corpus: very low-degree full diagonals, sparse full diagonals with
  bounded but meaningful row/column degree, nearly full diagonals with a large
  dense degree spike, and large nearly diagonal sparse-spike or dense-spike
  patterns. The nearly diagonal spike class also starts without BTF so KLS does
  not first pay for an AMD symbolic analysis before settling on the same
  METIS/no-BTF numeric path. Medium spiked low-diagonal patterns with a large
  row/column degree spike also start directly with METIS and max scaling for
  TSOPF-style power-grid structures. For large high-work no-BTF single-block
  analyses outside those direct-start classes, auto can also try METIS
  symbolically before the first numeric factorization and keep it when the
  symbolic fill score is clearly lower, avoiding a delayed post-factor METIS
  promotion. This preserves the papers' nested-dissection motivation without
  naming individual matrices.
- Constrained nested-dissection refinement: KLS can refine METIS rank groups
  with CAMD constraints, preserving nested-dissection rank shape while reducing
  local fill and flops. Large METIS orderings now use coarse rank-group CAMD
  refinement by default instead of limiting the CKTSO-style refinement to one
  medium structural class.
- Scaling policy: KLS exposes KLU scale modes and has auto scale selection based
  on pattern and numeric evidence, including no-scale reuse where repeated
  SPICE refactorization benefits. Low-work dominant-BTF cases also start with
  no KLU row scaling, avoiding the scale recomputation cost while preserving
  the useful BTF decomposition. Large high-work METIS/no-BTF single-block paths
  that already selected max scaling skip redundant post-factor scale trial
  factorizations.
- Fast repeated factorization with pivot check: KLS reuses an existing numeric
  pattern, checks reused pivots against the selected threshold via L
  multipliers, and falls back to full pivoting factorization if the reused order
  is unsafe and cannot be repaired locally. For unscaled BTF patterns, KLS can
  restart the rejected diagonal block with pivoting, rebuild the global row
  permutation/off-diagonal entries around that block repair, and retry the
  checked fast factorization. KLS records the first rejected factor-order pivot,
  original matrix column, and number of block restarts so future tail-restart
  work can distinguish late-tail failures from early failures.
- Static pivoting trial: KLS has value-aware greedy row matching, layered
  augmenting-path search for larger weak-diagonal candidates, and swap
  improvement for weak or high-off-diagonal-pivot medium matrices. Medium
  static-pivot candidates also run a bounded alternating-cycle pass that can
  apply profitable three- and four-row exchanges missed by pair swaps, while
  larger candidates keep the lower-fill layered/swap pattern. KLS can also
  trial dual-potential matching-derived row/column equilibration and keeps the
  transformed candidate only when numeric quality and cost evidence justify it.
  With `KLS_ENABLE_SPRAL_SCALING=ON`, large weak-diagonal candidates whose BTF
  analysis leaves one dominant block can run BSD-licensed SPRAL Hungarian
  matching/scaling before the first factorization; many-block BTF cases are
  excluded because the existing BTF path is already efficient there.
  The pre-factor static-pivoting gate also covers moderately sized matrices
  whose input values show a majority of weak or missing diagonal entries, plus
  medium-large mostly diagonal matrices with thousands of weak diagonal rows.
  These are general numeric-structure rules used by the frequency-domain and
  Rajat-family paper cases without naming individual benchmarks.
- Initial KLS-owned parallelism: repeated refactorization can run across
  independent BTF diagonal blocks for large high-flop cases where the block work
  is wide enough to offset thread overhead. The threaded path keeps a persistent
  worker pool on the solver instance so repeated SPICE refactors reuse workers
  and scratch storage instead of relaunching threads each cycle.
- KLS-owned serial refactor metadata: for unscaled numeric patterns that are
  not handled by the threaded BTF worker pool, KLS precomputes the fixed
  scatter from factor-order columns to pivotal rows and input value positions.
  This removes repeated `Q`/`Pinv` structure lookups from SPICE refactor cycles
  while preserving the existing KLU-derived LU storage. It covers both
  single-block and serial BTF refactors. For serial BTF refactors, KLS also
  prepartitions each mapped column into off-block and diagonal-block entries, so
  repeated refactors do not reclassify the same BTF structure in the numeric
  scatter loop.
- KLS-owned EGraph metadata and guarded refactor consumption: for threaded
  numeric runs, KLS can levelize the exact no-pivot refactor dependency graph
  from the actual U pattern after factorization. It retains level pointers and
  column lists and reports the level count, maximum level width, and dependency
  edge count in `kls_stats` and benchmark JSON. It also reports the
  CKTSO-style cluster/pipeline split point, pipeline-column count, approximate
  no-pivot update work, and tail work implied by the current level widths.
  Large single-block matrices with enough dependency work and level width can
  now consume this schedule through a work-estimated per-thread cluster-mode
  refactor,
  then switch to a no-pivot pipeline tail where each worker waits only for
  actual U-pattern predecessors. The same block-aware EGraph kernel can run
  inside a large dominant BTF block and update Offx for entries above that
  block. Schedule construction is now limited to single-block or dominant-BTF
  shapes with enough numeric work to consume it, while non-dominant many-block
  BTF and low-work dominant-BTF cases skip the setup and stay on the BTF worker
  pool or mapped refactor paths.
  This is still narrower than CKTSO's production pivoting machinery, but it is
  the first retained intra-block EGraph cluster/pipeline refactor path.
- SPICE-cycle orientation policy: KLS can analyze normal and transposed storage
  orientations and select the faster internal form for repeated solve cycles.
- LGPL project licensing and third-party notices: KLS itself is
  LGPL-2.1-or-later, with vendored-source attribution separated from the KLS
  license. MC64-equivalent preprocessing must stay inside that licensing
  boundary: HSL MC64 and solver-tree copies that retain HSL redistribution
  restrictions are not vendorable, while KLS can build the BSD-licensed SPRAL
  scaling subset from `third_party/spral` or use independent KLS code. The
  Rust `rwl/mc64` package is also BSD-licensed and useful as a reference, but
  it is a partial SPRAL translation and is not a better fit than the pinned
  SPRAL submodule for KLS's C/Fortran build.
- Paper-derived benchmark manifests: the full public SuiteSparse union from
  the local KLU, NICSLU, SubtreeLU, CKTSO papers and CKTSO ordering supplement
  resolves to 110 matrices. The routine medium subset contains 93 matrices, and
  the large supplement records the 17 excluded public paper cases for deliberate
  overnight tuning.

## Partially Implemented

- CKTSO-style static pivoting is only partial. KLS has a practical weighted row
  permutation, dual-potential matching-derived row/column equilibration, and
  optional SPRAL Hungarian/scaling trials, but not an always-on production
  MC64-style maximum-product matching algorithm with assignment dual scaling.
- NICSLU/CKTSO parallel scheduling is present at BTF-block granularity and, for
  large single-block cases, as a guarded exact-EGraph cluster/pipeline
  no-pivot refactor. This path now covers both unscaled factors and KLU
  row-scaled factors whose scale vector is recomputed before the EGraph refactor
  and permuted afterward. KLS still does not have the ETree descendant scheduler
  used by CKTSO-style pivoting tail restart.
- CKTSO fast factorization is present as pivot-checked reuse plus an unscaled
  BTF-block repair path. KLS does not yet implement CKTSO's pipelined tail
  factorization that restarts from the ETree descendants after a failed pivot
  check, and scaled fast-factor failures still use full fallback.
- SubtreeLU-style nested-dissection metadata is only used indirectly through
  METIS orderings and CAMD refinement. KLS does not yet preserve a separator
  tree for private/pipeline task queues.

## Not Implemented Yet

- Full production MC64-equivalent weighted matching and assignment-dual
  row/column scaling across the broad matrix set.
- NICSLU static-symbolic R1/R2 performance model for choosing sequential versus
  parallel numeric kernels before factorization.
- Intra-block parallel factorization with pivoting scheduled by an ETree.
- ETree-based pivoting tail restart beyond the current guarded no-pivot
  EGraph cluster/pipeline path.
- CKTSO dual-mode cluster/pipeline fast factorization with pivot check.
- CKTSO pipelined ETree-descendant tail factorization with pivoting after a
  pivot-check failure.
- CKTSO structure-adaptive hybrid parallel triangular solve.
- SubtreeLU separator-tree collapse for pivoting factorization.
- SubtreeLU FLOP-balanced separator-tree partitioning for refactorization.
- SubtreeLU constrained pivot search within nested-dissection subdomains.
- SubtreeLU supernodal updates with BLAS-style kernels inside the sparse
  up-looking framework.

## Tried During This Audit

An intra-block levelized no-pivot refactor prototype was tested using the
existing U structure as the exact dependency graph, matching the NICSLU/CKTSO
EGraph refactorization idea. Broad versions were removed before commit because
they regressed dominant-block circuit cases where the current BTF-block
threaded refactor is faster. A later version replaced per-column mutex
scheduling with static per-thread level slices, closer to CKTSO cluster mode,
but still regressed the intended single-block case (`rajat15`) by about 16% on
the focused repeated-refactor sample. Keeping those broad dispatch rules would
have required case-specific tuning, which is not the desired direction for KLS.

A broader KLS-owned serial no-pivot refactor path was also prototyped by
reusing the threaded BTF-block refactor kernel when thread-level parallelism was
not eligible. It passed correctness tests and helped some no-scale cases, but it
regressed other representative circuit cases such as `bcircuit` and `rajat03`.
The broad prototype was removed because it did not represent a general
improvement. A later narrower version kept only the precomputed unscaled serial
scatter metadata; same-session suite testing showed that subset as a modest
general win.

The refactor scatter map was also tried in scaled and threaded-worker forms.
The scaled extension improved a focused hard case but regressed the full suite,
and passing the map through the worker pool regressed representative threaded
BTF cases. KLS therefore keeps the map on the unscaled serial path and leaves
the threaded worker-pool scatter path separate.

The unscaled serial BTF scatter map was extended with a per-column boundary
between off-block entries and diagonal-block entries. This is a retained
general metadata improvement: a 25-matrix, 5-pass same-session extended-suite
A/B run improved the geomean from 0.03665s to 0.03651s and the median ratio to
0.993, with the largest win on `ckt11752_tr_0`. The single-block map path was
left in its original one-pass layout because single-block refactors have no
off-block entries and the partitioning setup did not help them.

A second layer of persistent refactor metadata was prototyped by precomputing
per-column L/U value and index pointers for mapped serial refactors. Focused
tests helped some repeated-refactor losses, but full-suite same-session testing
regressed the geometric mean because low-arithmetic and memory-sensitive cases
lost more than the hard cases gained. A structural gate based on flop count,
BTF shape, orientation, and arithmetic density reduced the regression but still
did not beat the committed baseline.

A narrower scaled single-block scatter map was also retried with a structural
gate for high off-diagonal-pivot cases. It repeatedly improved the targeted
scaled single-block case, but the added scaled path still perturbed unscaled
single-block hot cases enough to fail the no-regression bar. That prototype was
removed as well; scaled refactors still use the KLU-derived path unless they are
handled by the existing threaded BTF worker path.

A work-balanced BTF refactor scheduler was tested by sorting independent BTF
blocks by an LU-length work estimate before launching worker threads. It was
removed because the extra scheduling work did not improve the dominant-block
cases where BTF-block threading is eligible, and it regressed representative
threaded refactor timings on `ckt11752_tr_0` and `circuit_4`.

The reactive static-pivoting gate was also widened from medium matrices to
larger sparse matrices with high off-diagonal pivot counts. The existing
acceptance checks rejected `rajat22`, but only after paying a large matching
and refactorization trial cost; `rajat27` accepted a better numeric pattern but
lost on the repeated-SPICE metric because the one-time trial cost dominated.
The gate was restored to avoid converting matching into a broad overhead.

An auction-style weighted assignment pass was prototyped to move the current
greedy row matching closer to MC64's maximum-weight matching. It was removed
because the only current extended-suite matrices that select static pivoting
(`gemat11` and `gemat12`) already have good enough matched patterns; the extra
auction and dual-scaling work increased first-factor cost without improving
fill, refactor time, or residuals. A cheap weighted-gap guard avoided the
largest regression, but the guarded implementation still did not improve the
suite enough to justify the additional code.

The BTF-block threaded refactor scheduler was also changed from a mutex-protected
block counter to a C11 atomic work counter. This reduced scheduler overhead on
some BTF-threaded samples such as `ckt11752_tr_0`, but same-session median
suite testing showed no aggregate improvement and small regressions on other
cases. The mutex scheduler was kept until a broader scheduling change has a
clearer win.

The BTF-block parallel eligibility gate was relaxed to try more medium and
dominant-block structures. This was a general structural dispatch experiment,
not a matrix-name rule, but the added thread scheduling overhead regressed the
focused repeated-refactor samples that motivated the test. The conservative
high-flop, many-block eligibility gate was kept.

The no-pivot refactor update loops were also split into pivot-checking and
non-checking variants to reduce branch work in the fast repeated-factorization
path. The change passed correctness tests, but focused A/B timings were mixed
and did not show a general win; the simpler shared loop remains in place until a
larger KLS-owned numeric kernel makes this separation worthwhile.

The no-BTF symbolic retry was initially widened from single-block BTF cases to
dominant-block BTF cases, including METIS-started auto orderings. This is
consistent with KLU's observation that BTF can occasionally hurt, but the broad
version paid extra symbolic-analysis cost on unrelated multi-block matrices and
regressed the full 25-matrix same-session extended suite, with the candidate
geomean at 0.03816s versus the baseline at 0.03799s. A later narrower version
was retained: it only retries no-BTF for large BTF analyses whose largest block
covers at least 95% of the matrix with a small stripped fringe, and it only
accepts the no-BTF symbolic when the score improves by at least 20%. That
rescued the KLU-paper `Raj1` and `rajat24` cases from timeout/missing status
under the paper-medium run while keeping ASIC-style dominant-block cases on BTF.
The METIS promotion trial was also fixed to preserve the selected BTF/no-BTF
mode, which lets `Raj1` promote from no-BTF AMD to no-BTF METIS.

The pre-factor static row-matching path was also tested without the augmenting
and swap-improvement pass, leaving only the initial greedy maximum-value
matching. This reduced setup cost on some static-pivot samples, but after
rebuilding a clean baseline the 25-matrix extended suite regressed from
0.03798s to 0.03869s. The full improvement pass remains enabled for pre-static
matching.

After adding the paper-derived benchmark corpus, the pre-factor static
row-matching gate was widened for moderately sized matrices with at least half
of rows having weak or missing diagonal entries. This retained the full
augmenting/swap improvement pass and is a general numeric-structure rule, not a
benchmark-name rule. On the AT&T `onetone1`/`onetone2` paper cases, the
focused SPICE-cycle geomean improved from about 99.5s to 9.3s and both cases
selected the static row permutation. The retained gate was later widened to
cover the AT&T `twotone` scale as well; `twotone` now completes the paper-medium
suite run and selects static pivoting, reducing off-diagonal pivots from about
9500 to about 1600. A still-broader attempt to cover `mac_econ_fwd500` scale
was not retained because it timed out with multi-GB memory use.

The same pre-factor static row-matching gate was then extended to
medium-large mostly diagonal matrices with thousands of weak diagonal rows and
very few missing diagonals. This retained a structural diagonal-completeness
guard, so earlier Rajat-family cases with fewer weak diagonals remain on the
normal dynamic-pivot path. On the paper-medium corpus, the retained gate
improved the KLS geomean from about 1.57s to 1.51s and reduced the affected
Rajat cases' off-diagonal pivots to 1-3. KLS is still much slower than CKTSO on
these cases, so this is a partial static-pivoting improvement rather than a
replacement for full MC64-style matching/scaling or a faster KLS-owned numeric
kernel.

The large-matrix static row-matching augment was then changed from one
independent breadth-first search per unmatched row to a Hopcroft-Karp-style
layered augment. This is a retained MC64-adjacent improvement to the existing
greedy matcher, not full weighted MC64 dual scaling. It improved the AT&T
`twotone` paper case in a same-session focused run from about 273s to about
124s on the 100-step SPICE-cycle estimate by finding a lower-fill static row
permutation, while leaving `onetone1`, `onetone2`, and the large Rajat static
cases essentially neutral. The 25-matrix extended suite still passed, with a
geomean of about 0.04409s versus about 0.04421s for the previous mainline run.
A prototype greedy-plus-layered perfect matching for `mac_econ_fwd500` reduced
matched-METIS refactor time from about 22s to about 16s, but that remained far
behind the CKTSO artifact at about 2.5s, so the `mac_econ_fwd500` gate was not
widened.

For medium-large static-row-matched matrices whose diagonal is both weak and
mostly missing, the matching-equilibration trial was narrowed to keep the row
permutation but prefer no numeric scaling. This is a structural rule for
80k-150k order, at most 1.5M nonzeros, at least half missing diagonal entries,
and at least half weak-or-missing diagonal rows. On AT&T `twotone`, it selected
`scale=-1`, reduced off-diagonal pivots from about 1514 to 219, skipped the
matching-equilibration setup, and improved the focused SPICE-cycle estimate
from about 124s to about 96-98s. It deliberately does not apply to the Rajat
static-pivot cases, whose diagonals are almost complete and were slower when
forced unscaled. This remains far behind the CKTSO artifact at about 26s on
the same SPICE-cycle formula.

After expanding the paper-medium corpus with the public CKTSO ordering
supplement cases, the post-factor METIS promotion gate was widened for small
BTF-dominant matrices whose initial AMD/COLAMD factorization produces both
many off-diagonal pivots and high fill/flops. This is a retained combined
ordering rule, not a matrix-name rule: the gate requires at most four BTF
blocks, a dominant block covering at least 90% of the matrix, at least 128
off-diagonal pivots, and enough actual factor work to amortize the METIS trial.
It selected METIS on the TSOPF/QY power-grid-style cases in the expanded paper
medium suite, improving the 87-common-row KLS geomean from about 0.333s to
0.322s and the common KLS/CKTSO ratio from about 1.48x to 1.43x. Guard cases
with no off-diagonal pivot pressure, BTF disabled, or many BTF blocks stayed on
the existing AMD path.

A medium spiked low-diagonal structural METIS start was retained for the large
TSOPF-style paper cases. The rule requires a 50k-125k order matrix, about
20-32 entries per column on average, a 20-35% diagonal fraction, no empty rows,
and one large row/column degree spike covering about 40-60% of the matrix. It
starts with METIS while keeping BTF enabled, and now leaves auto scaling on the
default max-scaling path rather than forcing no-scale. In same-session checks,
`TSOPF_FS_b39_c19` improved from the earlier about-568s SPICE-cycle estimate to
about 546s, and the widened bound rescued the large-supplement
`TSOPF_FS_b39_c30` from the 120s per-process timeout with a SPICE-cycle
estimate of about 920s. This is still much slower than CKTSO on `c30`, whose
same-session artifact is about 333s, but it turns a KLS timeout into a
completed structural paper case. Smaller TSOPF/QY guard cases retained their
previous policies.

Two additional low-degree structural METIS starts were retained to remove the
last current-auto 120s paper-medium timeouts. The very-low-degree full-diagonal
METIS start now allows row/column degree up to 8, covering `ss1` while still
leaving broader low-degree circuit families on their existing policies. A
separate sparse-diagonal low-degree class covers 150k-250k order matrices with
5-8 entries per column on average, 5-20% diagonal coverage, no empty rows or
columns, and row/column degree at most 64; it starts with METIS and sum scaling.
Both classes skip redundant post-factor auto scale trials. In same-session
120s-capped runs, `ss1` completed in about 104s wall time with a SPICE-cycle
estimate of about 3413s, and `mac_econ_fwd500` completed in about 107s wall
time with a SPICE-cycle estimate of about 2252s. This is a timeout/completeness
improvement, not a CKTSO win: the existing CKTSO artifacts are about 349s and
265s respectively on the same SPICE-cycle formula. The same scale-trial gate
also reduced `G2_circuit` initial factor time from about 2.8s to about 1.0s and
modestly improved `mc2depi`.

A selected large-supplement reconnaissance manifest was added for the smaller
large paper cases that are practical under a 120s per-process cap before an
overnight full-large run. In same-session KLS/CKTSO runs with one factor and
one refactor repeat, KLS beat CKTSO on `rajat29` by about 5%, but lost the
other common successful cases: about 3.1x on `ASIC_680k`, 4.0x on `rajat30`,
4.5x on `G3_circuit`, and 4.7x on `nxp1`. KLS also timed out on
`TSOPF_FS_b39_c30` and `pre2`, both of which CKTSO completed within the same
120s cap; both solvers timed out on `Hamrle3`. The later widened TSOPF
spiked-low-diagonal rule removes the `c30` timeout, but `pre2` remains an
unresolved large-case timeout. The result reinforces that the remaining
large-case gap is mostly very large single-block or near-single-block
numeric/refactor throughput, not the many-small-BTF-block class.

The remaining `pre2` selected-large timeout was then checked as a simple policy
question before attempting new numeric-kernel work. Eight 120s-capped variants
all timed out with one factor and one refactor repeat: AMD/max-scale,
AMD/no-scale, COLAMD/max-scale, METIS with max/sum/no-scale, AMD/no-BTF
max-scale, and METIS/no-BTF max-scale. Focused sweeps on the completed but
slow `nxp1` and `rajat30` cases also found no better simple dispatch: `nxp1`
needs the current auto no-BTF METIS path, while AMD/COLAMD time out and
BTF-enabled METIS is much worse; `rajat30` remains best under current auto
METIS/no-BTF max-scaling. These results make `pre2`, `nxp1`, and `rajat30`
poor candidates for another ordering/scale/BTF heuristic. They need the open
paper ideas around faster single-block numeric/refactor kernels, matching
quality, or EGraph/separator-tree scheduling.

SCOTCH was then added as a pinned, reproducible optional ordering package and
tested as an explicit `--ordering scotch` path. SCOTCH symbolic analysis
completed on `pre2`, but numeric factorization still timed out at 120s; a
`rajat30` SCOTCH numeric run was also not competitive before interruption near
the same cap. The result does not justify adding SCOTCH to auto ordering yet.
ParMETIS was deliberately not added in this pass because it is an
MPI/distributed-memory package, while the current KLS benchmark and solver path
is shared-memory and single-process.

The CKTSO paper's nested-dissection-plus-constrained-minimum-degree idea was
then strengthened for KLS METIS orderings by applying coarse rank-group CAMD
refinement to large METIS orderings, not only one medium structural class.
This helped the completed large cases: same-session `--ordering metis --no-btf`
runs finished `nxp1` with about 1.16s factor and 1.13s refactor averages, and
`rajat30` with about 0.70s factor and 0.70s refactor averages. Auto selected
the same METIS/no-BTF path for those cases. The change did not fix `pre2`:
both explicit METIS/no-BTF and auto still timed out at 120s. This reinforces
that `pre2` is a missing-major-algorithm case rather than a separator-package
case.

The matching-derived equilibration pass was then moved closer to the
MC64/NICSLU/CKTSO preprocessing contract by first solving dual-potential
scaling constraints for the current greedy row match. When the constraints are
consistent, the matched diagonal is normalized to one and nonmatched entries
are bounded by the matched diagonal; when the greedy match leaves large
positive-cycle evidence, KLS falls back to the older heuristic balancing pass.
Same-machine checks on the retained static-pivot representatives (`gemat11`,
`gemat12`, `twotone`, and `rajat25`) produced the same fill, flops,
off-diagonal pivot counts, conditioning, and residuals as the previous
mainline. A temporary large-gate experiment also let `pre2` try this
preprocessing path with layered matching, but the factor-only run still timed
out at 120s. This makes the dual-potential pass a bounded MC64-adjacent
preprocessing cleanup, not the missing CKTSO-scale algorithm.

A bounded alternating-cycle improvement was then added after the existing
layered cardinality augment and pair-swap pass. This is not a full MC64
shortest-augmenting-path implementation: it only accepts profitable local
cycles of up to four rows, and it is intentionally limited to `n <= 50000`.
Same-session A/B against the previous commit showed useful medium effects:
`gemat12` reduced off-diagonal pivots from 11 to 7 and improved reciprocal
condition evidence, and `onetone2` reduced actual fill/flops. The same ungated
pass was rejected for larger static-pivot cases because `twotone` and
`rajat25` increased fill/flops despite comparable or better pivot counts, so
large static-pivot candidates stay on the lower-fill layered/swap pattern. A
large `pre2` analyze-only check still shows a dominant 629628-row block and
about `2.08e11` estimated flops under the current AMD/BTF symbolic path, so the
remaining CKTSO gap is still a missing numeric-kernel/scheduler issue.

A cheap exact sparse assignment stage was then added ahead of the greedy
static row matcher for small candidates. It solves a min-cost augmenting-path
problem on `col_max - log(abs(a_ij))`, which is equivalent to maximum-product
matching when a full assignment is found. The path is intentionally gated to
`n <= 4000`, `nnz <= 75000`, and `n * nnz <= 5e7`; a same-session HB/gemat
experiment with a wider gate improved `gemat12` off-diagonal pivots from 7 to
3, but raised initial factor/preprocessing time from about 0.034s to about
5.2s. With the retained gate, `gemat11` and `gemat12` stay on the previous
fast matcher, while synthetic 3000-row static-pivot smoke tests exercise exact
matching. This closes a small piece of MC64 functionality but confirms that
KLS still needs an optimized production MC64-equivalent matcher, not the
straightforward min-cost implementation, before the CKTSO-scale gap can close.

The exact assignment implementation was then changed from a generic
source/sink residual graph to a KLS-owned sparse row/column
shortest-augmenting-path matcher. This removes the extra source/sink edges and
keeps the code LGPL-compatible, but the policy gate remains conservative. A
same-machine wider-gate retest with the corrected source/sink semantics
selected exact matching for `gemat11` and `gemat12`. `gemat12` improved from 7
to 3 off-diagonal pivots, but initial factor/preprocessing time was still about
4.24s; `gemat11` stayed at 0 off-diagonal pivots and also paid about 4.24s.
Restoring the small gate returned `gemat11`/`gemat12` to the prior fast path
(`selected_exact_matching=false`, 0 and 7 off-diagonal pivots), while
`onetone2`, `twotone`, and `rajat25` also stayed off the exact path. The
lesson is that maximum-product matching alone is too expensive to use as a
medium/large default; KLS still needs MC64-quality scaling and acceptance plus
a more optimized assignment implementation before this can close the CKTSO
gap.

An optional SPRAL hook was then added for BSD-licensed matching/scaling
support. It is deliberately not a default dependency and not a solver
replacement.
KLS first wired SPRAL auction matching as a large structural-deficit fallback,
then added SPRAL's MC64-like Hungarian unsymmetric matcher/scaler. The hook was
initially system-SPRAL only; KLS now also pins upstream SPRAL as a submodule and
builds just its scaling subset for reproducible LGPL-compatible MC64-adjacent
experiments. A wider same-cardinality Hungarian replacement was tested on the
AT&T `onetone2` and `twotone` static-pivot cases: it could reduce off-diagonal
pivots on `twotone`, but increased fill and refactor time, and it worsened
`onetone2`. The retained policy therefore uses optional SPRAL Hungarian when
it improves matching cardinality over the in-tree matcher, and can also factor
a same-cardinality SPRAL Hungarian candidate after the first numeric factor
only for expensive high-off-diagonal-pivot cases. The latter path is gated by
actual factor work/fill and still keeps the SPRAL candidate only when the
factored numeric evidence improves. This keeps a license-compatible
MC64-adjacent source available for hard structural-deficit and expensive
dynamic-pivot cases without letting a maximum-product match replace already
accepted KLS row matchings solely on weight.

The SPRAL Hungarian/scaling path was then promoted to a bounded pre-factor
large-matrix candidate for weak-diagonal matrices whose symbolic analysis
leaves one dominant BTF block. This is the closest retained KLS path to the
MC64 preprocessing described by NICSLU and CKTSO, while remaining
LGPL-compatible because it uses the BSD-licensed SPRAL scaling subset rather
than HSL MC64 or restricted solver-tree copies. Same-session checks show the
tradeoff clearly. With `KLS_ENABLE_SPRAL_SCALING=ON`, `mac_econ_fwd500`
completed a one-factor/one-refactor run inside the 120s cap, selecting static
pivoting and reducing repeated factor/refactor averages to about 5.22s; the
same command timed out at 120s in the no-SPRAL build. `rajat30` selected the
SPRAL static match, reduced off-diagonal pivots to one, and modestly improved
repeated factor/refactor times to about 0.74s/0.73s. A broad version also
matched `ASIC_680k`, but that regressed an already fast many-block BTF case
from about 0.15s repeated refactors to about 0.27s, so the retained gate now
requires no BTF or a dominant BTF block. `pre2` still timed out at 120s after
trying the SPRAL path; its analyze-only evidence remains a dominant
629628-row block with about `2.08e11` estimated flops. This confirms that
license-compatible MC64-adjacent preprocessing is useful and worth keeping, but
it does not replace the missing CKTSO/SubtreeLU-style intra-block
numeric/scheduling machinery.

The fast-factor pivot-check path was then made more diagnostic by recording the
first rejected factor-order pivot and original matrix column in `kls_stats` and
benchmark JSON. This does not implement CKTSO's pipelined tail factorization,
but it is a required prerequisite: KLS can now measure whether failed fast
factorizations reject near the tail, where an ETree-descendant restart could
avoid recomputing the whole matrix, or near the front, where full fallback is
still expected.

The same fast-factor path was then extended for scaled serial refactors. When a
scaled pattern is using fast factorization with pivot checks, KLS now runs its
own checked refactor loop and interrupts at the first unsafe multiplier instead
of running a complete KLU refactor and scanning the completed factors
afterward. This is still not CKTSO tail restart because fallback remains a full
pivoting factorization, but it narrows the wasted work before fallback and uses
the same rejected-pivot coordinate needed by a future ETree-descendant restart.

The unscaled fast-factor path was then extended with a block-local restart
primitive. When a checked no-pivot fast factorization rejects a reused pivot in
a BTF diagonal block, KLS can refactor only that BTF block with pivoting,
splice the block's new row order into the global numeric permutation, rebuild
the unscaled off-diagonal entries from the updated inverse permutation, and
retry the checked fast factorization. Smoke tests now cover both the unscaled
repair and the scaled fallback path, and representative static-pivot cases
(`gemat12`, `onetone2`, `twotone`, and `rajat25`) did not trigger unexpected
block restarts. This is useful CKTSO-aligned infrastructure, but it is still
not CKTSO's production tail restart: it does not retain an ETree/EGraph tail,
does not restart only descendant work inside a large single block, and does
not handle scaled repairs.

The selected-large KLU2 comparison was also run with the same 120s cap and one
factor/refactor repeat. KLU2 completed only `rajat29`, `rajat30`, and
`ASIC_680k`; it timed out on `G3_circuit`, `pre2`, `nxp1`, `Hamrle3`, and
`TSOPF_FS_b39_c30`. On successful common rows, KLU2 was faster than KLS on
`rajat29` but slower on `rajat30` and `ASIC_680k`. After the widened TSOPF
spiked rule, KLS also completes `TSOPF_FS_b39_c30` where KLU2 timed out. The
remaining large-case primary competitor is therefore CKTSO rather than KLU2.

Extending the existing mapped refactor metadata to scaled refactors was
retested against current mainline using a clean `HEAD` worktree. The scaled
map was KLU-semantics-compatible after recomputing row scale factors before
the mapped scatter and permuting them afterward, but it did not improve the
focused scaled large/medium cases. On the seven-case focused set the candidate
geomean regressed by about 1%, with no wins over 2%, so the experiment was
removed again. This keeps the retained refactor map limited to unscaled serial
patterns until a broader EGraph/separator-tree numeric kernel exists.

A narrower EGraph refactor consumer was then retained for the high-work
unscaled single-block class. It reuses the exact U-pattern level schedule,
assigns each level to static per-thread slices, and uses private dense scratch
per worker with a barrier between levels. The gate requires a single BTF block,
no KLU row scaling, at least 100k rows, at least `1e9` estimated no-pivot
dependency work, and level width at least four times the requested thread
count; scaled cases, many-block BTF cases, and smaller matrices continue using
the existing mapped or BTF-worker paths. In same-session checks on the default
build, `G3_circuit` improved from the earlier about 34s/34s factor/refactor
averages to about 27.9s/28.1s with valid residuals. Guard cases stayed on their
previous paths: `nxp1` and `rajat30` remained scaled, `ASIC_680k` remained
many-block BTF, and small `bcircuit` remained below the large-case threshold.
This is still not the full CKTSO cluster/pipeline scheduler, but it is a
retained general implementation of one paper idea where the current metadata
shows enough work to amortize synchronization.

The retained EGraph refactor was then extended with a CKTSO-style no-pivot
pipeline tail. Levels before the recorded split still use barriered cluster
mode. Tail columns are assigned to workers without per-level barriers, and a
column waits only for the actual U-pattern predecessors it consumes. Completion
is published with C11 atomics, so LU writes from predecessor columns are visible
before dependent columns update. On `G3_circuit`, the same one-factor and
one-refactor focused check improved from about 27.9s/28.1s after the first
EGraph path to about 19.5s/19.4s with valid residuals. `G2_circuit` also
completed correctly on the same path. Scaled single-block and many-block guard
cases (`nxp1`, `rajat30`, `ASIC_680k`, and small `bcircuit`) stayed on their
existing paths. This closes the no-pivot EGraph cluster/pipeline piece for the
current KLS-owned storage, but not CKTSO's ETree-descendant restart with
pivoting after a failed pivot check.

The same EGraph cluster/pipeline refactor was then extended to large
single-block KLU row-scaled factors. The scaled path recomputes KLU's row
scale vector before the EGraph numeric update, divides each mapped matrix entry
by the unpermuted scale for its original row, and then permutes `Rs` back into
pivot order after the refactor. The normal mapped refactor path still rejects
scaled factors, so scaled maps are built lazily only when the EGraph path runs.
On same-session focused checks, `nxp1` improved from about 1.02s/1.01s
factor/refactor averages to about 0.73s/0.71s, and `rajat30` improved from
about 0.65s/0.64s to about 0.38s/0.39s, all with valid residuals. Guard cases
remained on their existing paths: `ASIC_680k` stayed a many-block BTF case,
`rajat29` stayed below the high-work EGraph gate, and small `bcircuit` remained
below the size threshold.

The EGraph pipeline tail then replaced mutex-based stop polling with an atomic
stop flag while keeping error details protected by the existing mutex. This is
a small synchronization reduction in the common no-error path: workers no
longer take a mutex while polling predecessor completion or checking level
boundaries. Same-session focused checks improved `G3_circuit` to about
19.0s/18.9s factor/refactor, `nxp1` to about 0.72s/0.71s, and `rajat30` to
about 0.37s/0.36s, with valid residuals.

A follow-up attempt to cache the per-worker EGraph dense scratch vectors was
rejected. A dedicated solver-owned scratch cache improved repeated `rajat30`
factor/refactor averages to about 0.34s/0.34s, but `nxp1` regressed slightly
to about 0.72s/0.72s and the long `G3_circuit` guard regressed to about
19.0s/19.2s. Since the effect was not a general win, the experiment was
removed instead of adding a size or matrix-shape gate.

The MC64 compatibility boundary was rechecked after allowing existing code if
it remains LGPL-compatible. The retained vendored route is still SPRAL's
BSD-3-Clause scaling subset: it is redistribution-compatible with KLS's
LGPL-2.1-or-later license, whereas HSL MC64 itself and restricted MC64 copies
from solver trees remain out of scope for vendoring. The policy does not
require every MC64-style implementation to originate in KLS; it requires any
copied or vendored implementation to be redistributable inside an LGPL KLS
distribution. A small BSD Rust `mc64` crate exists as a partial SPRAL
translation, but it does not improve KLS's C integration story over the already
pinned SPRAL Fortran/C interface. Fresh SPRAL-enabled checks also confirm the
policy should stay guarded rather than become an unconditional default:
`rajat30` selected SPRAL matching, reduced off-diagonal pivots to one, and cut
initial factor time to about 5.7s, but its repeat-heavy factor/refactor
averages were about 0.67s/0.66s versus the faster current no-SPRAL EGraph path.
`nxp1` did not select SPRAL and stayed roughly neutral-to-slightly-worse. This
keeps license-compatible MC64-style code in KLS, but points the large remaining
CKTSO gap back to numeric scheduling and pivoting machinery rather than merely
importing another MC64 copy.

The EGraph refactor consumer was then generalized from hard-coded single-block
indices to block-local BTF indices. The same guarded cluster/pipeline schedule
can now run inside a large dominant BTF block, update the block-local LU
columns, and refresh Offx for entries above that BTF block. The dispatch remains
conservative: it requires a large block covering at least 75% of the matrix, so
ordinary many-block cases such as `ASIC_680k` stay on the existing BTF worker
pool. On same-session checks, the committed HEAD baseline timed out at 130s on
`mac_econ_fwd500` with one factor and one refactor, while the block-aware
EGraph path completed with valid residuals and about 14.9s/14.9s repeated
factor/refactor averages. Single-block guards remained in their previous
range: `nxp1` measured about 0.71s/0.71s and `rajat30` about 0.34s/0.33s in
focused repeated runs; `G3_circuit` stayed around 19.1s/19.1s. This is useful
dominant-block coverage, but not the missing major CKTSO algorithm: `pre2`
still timed out under a 125s cap.

A current large-recon run after the dominant-BTF change still shows the
remaining gap clearly. KLS completed six of the eight selected large cases and
timed out on `pre2` and `Hamrle3` under a 120s per-matrix cap. Against the
existing CKTSO four-thread run, KLS still lost all six common completed cases,
with the largest ratios on `nxp1`, `ASIC_680k`, `G3_circuit`, and `rajat30`.
The same run beat KLU2 on the larger common cases `ASIC_680k` and `rajat30`,
but still lost `rajat29`. A factor-only `pre2` probe also timed out under
125s, confirming that the unresolved `pre2` gap is first-factor numeric
machinery, not repeated-refactor scheduling.

A low-work dominant-BTF guard was then retained for Rajat-family large cases:
if BTF finds one block covering at least 95% of the matrix, the stripped fringe
is at most 5% but still nontrivial, and the symbolic flop estimate is below
`1e9`, auto keeps BTF and starts with no KLU row scaling. On `rajat29`, this
changed auto from AMD/no-BTF/max-scale to AMD/BTF/no-scale, reducing the
SPICE-cycle estimate from the old about 5.67s to about 3.74s. That now beats
the saved CKTSO artifact at about 5.28s and the saved KLU2 artifact at about
4.29s. Higher-work guards such as `rajat30`, `nxp1`, `Raj1`, and `rajat24`
still take the no-BTF retry when their symbolic evidence supports it.

A symbolic METIS retry was then retained before numeric factorization for
large high-work single-block no-BTF analyses. This is the same CKTSO/SubtreeLU
nested-dissection direction as the existing post-factor METIS promotion, but
it avoids first paying for an AMD numeric factorization when METIS already has
a clearly lower symbolic fill score. On same-session serial checks, `nxp1`
kept the same METIS/no-BTF numeric path but reduced the SPICE-cycle estimate
from the saved about 89.6s to about 78.9s, and `rajat30` moved from about
49.2s to about 42.4s. The low-work dominant-BTF guard stayed on `rajat29`, and
the many-block BTF guard stayed on `ASIC_680k`, so this is a general
symbolic-cost improvement rather than a benchmark-name policy.

The auto-scale gate was then tightened for the same large high-work
METIS/no-BTF single-block class. Earlier sweeps had already shown the completed
hard large cases should keep max row scaling, but auto mode still tried other
scale modes after the first factorization and rejected them. Skipping those
post-factor scale trials preserves the final factors and residuals while
cutting setup time: `rajat30` auto initial factor moved from about 4.06s to
about 1.89s, matching explicit METIS/no-BTF/max-scale, and `nxp1` moved to
about 2.77s while keeping the same METIS/no-BTF/max-scale EGraph path.
`rajat29` and `ASIC_680k` remained on their retained BTF policies.

The high-work METIS/no-BTF class was then moved one step earlier in analysis.
A structural direct-start gate now recognizes large nearly diagonal matrices
with a meaningful row/column spike in either a sparse-spike or dense-spike
density band. This lets `nxp1` and `rajat30` start directly with METIS/no-BTF
instead of first doing AMD symbolic analysis and then a METIS symbolic retry.
Same-session checks reduced `rajat30` analysis from about 3.9s to about 3.0s
while preserving the about 1.86s initial factor and about 0.31s refactor path,
and reduced `nxp1` analysis from about 2.6s to about 1.84s while preserving
the max-scale EGraph path. The low-work `rajat29` BTF guard and many-block
`ASIC_680k` guard stayed on their retained policies.

EGraph schedule construction was then tightened to the same structural class
as the retained EGraph consumer. KLS had been building dependency-level
metadata for many-block and low-work BTF cases that could not use the
cluster/pipeline path, including `ASIC_680k` and the retained `rajat29`
low-work dominant-BTF policy. The new gate keeps schedules for high-work
single-block cases such as `nxp1` and `rajat30`, but skips them for the BTF
worker-pool and low-work mapped paths. Same-session checks showed schedule
metrics dropping to zero for `ASIC_680k` and `rajat29`, with `ASIC_680k`
initial factor/setup moving from about 1.59s to about 1.54s and `rajat29`
from about 0.26s to about 0.24s, while the high-work EGraph guards retained
their schedule metadata.

A CKTSO-style dynamic atomic assignment prototype for the EGraph pipeline tail
was tested and rejected. It replaced the static per-thread tail stride with a
shared atomic cursor, but `nxp1` and `rajat30` were neutral-to-slightly worse
and `G3_circuit` regressed to about 19.45s/19.41s factor/refactor. The
existing static tail assignment therefore remains the better fit for KLS's
current column storage and scratch model.

The BTF worker pool was then adjusted to fetch small ranges of diagonal blocks
per mutex acquisition when a matrix has many thousands of non-dominant BTF
blocks. This targets scheduler overhead on ASIC-style matrices with hundreds
of thousands of tiny blocks, while dominant-block and smaller-BTF cases still
fetch one block at a time for load balance. On same-session checks against a
clean baseline, `ASIC_680k` moved from about 0.173s/0.157s repeated
factor/refactor averages to about 0.159s/0.156s, while `ASIC_680ks` moved from
about 0.163s/0.170s to about 0.173s/0.163s. The net SPICE-cycle effect is
small but positive on the many-block ASIC guards, so this is a retained
scheduler-overhead cleanup, not a CKTSO-scale algorithmic fix.

The EGraph cluster scheduler was then adjusted from equal per-thread column
slices to contiguous slices balanced by the existing no-pivot column-work
estimate. This is a small retained part of the CKTSO/SubtreeLU load-balance
idea, scoped to the barriered cluster levels and leaving the existing static
stride pipeline tail unchanged. Same-session clean-baseline checks showed
modest but consistent wins on the EGraph path: `rajat30` moved from about
0.335s/0.330s repeated factor/refactor averages to about 0.317s/0.314s,
`nxp1` moved from about 0.710s/0.706s to about 0.691s/0.669s, and
`G3_circuit` moved from about 19.40s/19.33s to about 18.41s/18.45s. This
improves the retained no-pivot refactor consumer, but it does not address the
remaining first-factor timeout on `pre2` or implement CKTSO's pivoting tail
restart.

The EGraph worker launch path was also retried with a solver-owned persistent
worker/scratch pool, analogous to the retained BTF refactor pool. This was
removed before commit because it was not a general win after the cluster
work-balancing change: `rajat30` improved slightly from about 0.325s/0.313s
to about 0.316s/0.310s factor/refactor averages, `nxp1` was neutral to
slightly worse at about 0.685s/0.686s versus 0.684s/0.685s, and `G3_circuit`
regressed from about 18.40s/18.49s to about 18.73s/18.69s. The evidence points
back to numeric update structure, pivoting-tail machinery, or separator-tree
scheduling rather than thread-launch overhead.

The EGraph cluster/pipeline split threshold was also checked after the
work-balanced cluster slices. The retained policy switches to pipeline mode at
the first level narrower than `2 * threads`. A later handoff at `1 * threads`
cut the `rajat30` pipeline tail to 1156 columns but worsened factor/refactor
averages to about 0.346s/0.330s. An earlier handoff at `4 * threads` expanded
the tail to 3715 columns and measured about 0.328s/0.317s. The current
middle split remains the better general setting in this quick check.

The generic EGraph column kernel was also tested with a single-block fast path
that bypassed the per-column BTF block lookup and `R` checks when
`nblocks == 1`. This was removed before commit because it did not help the
scaled single-block guard: `rajat30` was only slightly positive at about
0.317s/0.313s versus 0.322s/0.314s, while `nxp1` regressed from about
0.671s/0.675s to repeated samples around 0.693s/0.675s and 0.696s/0.693s.
The block-lookup branch is therefore not the next useful source of the CKTSO
gap.

The EGraph pipeline completion array was also tested with retained per-column
level labels so cluster-phase columns could skip atomic completion stores and
pipeline waits for predecessors known to be before the split level. This was
removed before commit because it helped `rajat30` and repeated-factor `nxp1`
but regressed the largest EGraph guard: `rajat30` moved from about
0.316s/0.314s to 0.310s/0.305s, `nxp1` moved from about 0.687s/0.673s to
0.665s/0.674s, but `G3_circuit` regressed from about 18.39s/18.38s to
18.73s/18.75s. The retained pipeline publication remains the simpler
per-completed-column atomic store until a fuller scheduler changes the tail
execution model.

The remaining `pre2` timeout was rechecked after the MC64 licensing boundary
was clarified. A SPRAL-enabled build, using the retained BSD-licensed
Hungarian/scaling path, still timed out under a 180s one-factor cap. Lowering
the initial pivot tolerance also did not provide a usable CKTSO-style
pivot-reuse substitute: `--pivot-tol 0` reached a singular factor quickly,
while `1e-8` and `1e-4` still timed out under 120s. On completed large guards,
lower tolerances were not a general win: `rajat30` at `1e-4` slowed initial
factor and repeated refactor, `rajat30` at `1e-8` improved only refactor while
worsening initial factor and conditioning, and both `nxp1` lower-tolerance
checks regressed. This keeps `pre2` in the missing pivoting-tail/numeric-kernel
bucket rather than the tuning bucket.

A narrow many-block BTF worker-map retry was also prototyped for ASIC-style
structures after the CKTSO comparison showed a large `ASIC_680k` gap. The
prototype made the worker map path scale-aware and built the precomputed
input-position map only for non-dominant BTF patterns with at least 100k blocks
and at least half as many blocks as rows. It was removed because it was not a
general win: `ASIC_680k` had only a small refactor improvement but slower
repeated factor and setup, while `ASIC_680ks` regressed in both factor and
refactor. This confirms the previous broad worker-map rejection and points
ASIC-style gaps toward a different small-block numeric/storage kernel rather
than passing the existing map through the worker pool.

## Recommended General Work

1. Evolve the retained EGraph metadata consumer into a fuller KLS-owned numeric
   engine: implement CKTSO's ETree-descendant pivoting tail restart or
   SubtreeLU's separator-tree private/pipeline scheduler, not isolated
   benchmark guards. The current unscaled block-local restart and large
   single-block no-pivot cluster/pipeline refactor are useful precursors, but
   the target is ETree-descendant tail restart inside large blocks.
2. Continue turning matching/scaling into a production MC64-equivalent stage,
   but keep it inside the LGPL-compatible boundary: use the BSD-licensed SPRAL
   scaling submodule or independent KLS code, not HSL MC64 or restricted MC64
   copies from other solver trees. The retained SPRAL path now helps large
   weak-diagonal dominant-block cases, but `pre2` still times out, so matching
   quality alone is not the remaining CKTSO-scale gap.
3. Add a structure-adaptive triangular solve only after the LU storage owned by
   KLS exposes row-oriented or segmented access cheaply.
4. Use static symbolic and numeric-cost models to decide whether a parallel
   kernel should run, so KLS avoids matrix-name-specific tuning.
