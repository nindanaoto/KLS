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
refactors and for a narrow scaled many-fringe dominant-BTF subset, records an
exact no-pivot EGraph level schedule from the numeric U pattern, and consumes
that schedule in a guarded cluster/pipeline intra-block refactor path for large
single-block matrices, including KLU row-scaled cases whose scale factors can
be recomputed and permuted safely, and inside large dominant BTF blocks whose
off-block entries can be refreshed from the retained map. A narrow many-fringe
dominant-BTF shape can also consume the same exact EGraph in an all-pipeline
mode, and very high-work unscaled single-block factors can do the same,
avoiding cluster barriers and waiting only on actual U-pattern predecessors.
KLS also now schedules extremely fragmented unscaled BTF refactors with one
substantial but non-dominant block through the exact EGraph path when measured
work is high enough, so the large block is not left serial behind hundreds of
thousands of singleton blocks. Auto scaling now starts those full-rank
fragmented BTF shapes in no-scale/no-recheck mode, because KLU row scaling
prevents the EGraph refactor from running and roughly doubled repeated refactor
time on the ASIC 680k/680ks structural class.
The row-permuted solve wrapper now keeps a solver-owned dense permutation
workspace instead of allocating it on every forward or transpose solve, trimming
repeated SPICE-cycle solve overhead for static-pivoting cases while leaving the
KLU triangular kernels unchanged.
Solve calls also now update only solve timing, status, and memory counters
instead of refreshing the full numeric stats snapshot after every triangular
solve.
The barriered EGraph
cluster levels now use FLOP-estimated per-thread slices instead of equal column
slices, and the no-pivot pipeline tail now uses an atomic dynamic work cursor
instead of static per-thread strides, which are retained pieces of the
CKTSO/SubtreeLU load-balance idea.
The same schedule pass now records adjacent row-major U-pattern supernode
candidate counts, covered rows, maximum width, dense-block entries, and shared
trailing entries. The gated row refactor also reports executable segment groups,
segment-covered rows, segment maximum width, and dense/trailing work. These are
diagnostic bridges toward SubtreeLU-style row/segment storage and BLAS-friendly
updates. Dense-eligible gated row-refactor groups now also compute their
internal dense `L` block, dense `U` block, and shared trailing `U` block
directly in KLS-owned row-major mirrors before scattering back to KLU on
success. Row-refactor groups now also retain a reverse group graph and report
group dependency edges, root groups, leaf groups, and maximum group fanout,
giving the row/segment layer the explicit task graph needed by future
SubtreeLU-style private/pipeline queues and CKTSO-style tail schedulers. The
experimental row pipeline can now consume that graph through a bounded
successor-ready queue for the pipeline tail, and reports whether the last row
run used that queue plus cumulative queued group counts. Unchecked queued row
tails also skip the per-row completion bitmap, while checked row runs still keep
it for prefix-reject validation and report that usage separately. Row-pattern
analysis also records when input columns are already structurally covered by
`L`, the pivot, or `U`, allowing row kernels to skip redundant residual cleanup
loops. Dense row segments now keep deferred all-value scatter only for checked
pivot probes, while unchecked row refactors scatter completed dense rows
directly. Dense and sparse row segments now retain per-input-entry target maps
for rows whose raw values can be placed exactly into external work-vector
slots, row-major `L` mirrors, the pivot, or in-segment/trailing `U` mirrors;
the direct loaders consume those maps on repeated refactors instead of
rediscovering the segment slots. The row ready queue now also orders initially
ready tail groups and newly released successor groups by the retained
FLOP-style group work estimate, which moves the row scheduler closer to
SubtreeLU's workload-balanced queue generation. Checked row fast-factor rejects
now also report a
conservative group-tail restart scope from the retained row-group successor
graph, giving future CKTSO-style pivoting tail work an explicit row/segment
task-tail measurement. Checked KLS-owned row fast-factor/refactor passes also
test the guessed diagonal against the maximum absolute value already present in
the current U row before publishing it, matching CKTSO's row-wise pivot check
rather than relying only on later L-multiplier growth. The row ready queue also
keeps solver-owned workspace across repeated row refactors, avoiding
queue/bitmap/predecessor allocation churn in the experimental row scheduler.
KLU-compatible checked fast-refactor paths now enforce the same row-wise U
acceptance rule after the finalized U rows are available; a violation is
reported as an all-current fast reject and enters the existing KLS block/tail
repair path. This closes the pivot-acceptance mismatch for the column kernels,
but it is still delayed validation rather than the row-major CKTSO numeric
executor that can avoid doing rejected work in the first place.
Full-graph queued row runs now consume cached
root groups through a private-root cursor before falling back to the shared
ready queue for successor-released groups, trimming the first wave of shared
queue traffic without changing the retained row DAG. When a completed group
releases successors, the worker now keeps one ready successor as a local
continuation and enqueues the rest, moving the queued DAG scheduler another
step toward SubtreeLU-style private/pipeline execution; stats now expose how
often that local continuation path is used. Checked queued rejects now refresh
any missing prefix rows before accepting the prefix-tail repair classification,
so work ordering cannot turn an already-repairable prefix into a scheduler-race
miss. Experimental row refactors now cover single-block
factors and BTF
diagonal blocks; BTF off-block values are refreshed into KLU `Offx` from the
retained input map. Unchecked row refactors can also hand dirty KLS-owned
row-major `L`/`U` mirrors directly to guarded forward and transpose
solves, leaving KLU's column values stale until later factorization or other
non-row fallback needs them. For BTF factors the dirty row solve follows KLU's
block order and uses refreshed `Offx` coupling directly. This avoids making
eligible row-refactor solves pay an immediate KLU publish. KLU row-scaled row
refactors can now also recompute `Rs`, use unpermuted row scales for
fixed-position input loads, and restore `Rs` to pivot order after an accepted
pass. Dirty row solves consume both unscaled and KLU row-scaled normal
row mirrors directly, using KLU's pivot-order `Rs` semantics for the
right-hand-side load and transpose output. The guarded row-major solve also
runs under KLS's external static-matching row permutation and matching-derived
row/column scaling because `solve_impl` wraps the kernel solve with the same
right-hand-side/result transforms used by KLU. These pieces still do not change
the current default KLU-column numeric kernel.
KLS now also records `initial_factor_path`, `last_factor_path`, and largest
ordered diagonal-block factor ETree counters in `kls_stats`/`kls_bench`. These
diagnostics expose whether a solve still entered the KLU-derived pivoting
kernel (`klu_first` or `klu_fallback`) and the ETree upper-bound shape CKTSO
uses for pivoting-tail scheduling; they are evidence for the remaining
row/up-looking first-factor work, not a substitute for that kernel.
An env-gated `KLS_ENABLE_KLS_FIRST_FACTOR=1` path can now allocate and assemble
KLU-compatible numeric storage itself for no-scale and KLU row-scaled first
factorizations, handle singleton BTF blocks directly, and use the KLS-owned
pivoted block kernel for multi-column BTF blocks. It reports `kls_first` when it
succeeds, seeds KLS-owned row-major `L`/`U` value mirrors for guarded solves and
unchecked and checked repeated refactors, and falls back to the KLU first-factor
path otherwise. This is a first KLS-owned factorization scaffold, not the
default production row-major CKTSO-style factorization.
When this row up-looking first-factor path performs a dynamic column pivot and
the retained METIS `NodeNDP` map covers the full factor order, it now uses
SubtreeLU Algorithm 4's scoped `N'` maximum inside the current collapsed
component extent instead of the whole updated row; if no safe scoped candidate
exists for a required pivot, it rejects the row-up attempt instead of crossing
the separator pivot domain. When a separator-pipeline row has been safely
pre-updated and then needs a scoped pivot, the pipeline now publishes that row
under the ordered lock, updates both the committed and phase-local prefix U
entries for the column exchange, advances a column-order epoch, and makes
speculative suffix rows that began under an older epoch discard and restart
inside the same worker phase instead of tearing down and relaunching the whole
pipeline suffix.
Benchmark stats separate separator-domain dynamic pivots, would-have-crossed
fallback candidates, and strict rejects.
Its static-pivot
preprocessing has a cheap exact sparse maximum-log-product assignment path for
small candidates, derives MC64-style row/column scales from that exact
assignment's dual potentials when the scaled candidate is accepted, and can
improve medium row matchings with bounded alternating cycles beyond the
pair-swap pass. KLS now enables the pinned BSD-licensed SPRAL
scaling subset by default, so it can use Hungarian matching/scaling as a
pre-factor MC64-adjacent candidate for small and medium weak-diagonal matrices,
SPRAL auction matching/scaling for very large weak-diagonal dominant-block
matrices, and a post-factor value-gated Hungarian trial for dense
high-off-diagonal-pivot cases while still allowing
`KLS_ENABLE_SPRAL_SCALING=OFF` builds. Fast factorization can now repair an
unsafe unscaled BTF diagonal block by
restarting that block with pivoting and then retrying the checked no-pivot
factorization when later columns may not have been refreshed. If the failed
pass had refreshed all columns, or if a serial BTF pass rejected in the final
block, KLS validates the pivoted repair block tail directly and skips the
redundant checked refactor pass.
The same block-local repair now also covers KLS-owned scaled checked-refactor
failures; when the rejected pass left a valid prefix-current state, KLS
recomputes row scales and continues with the threaded checked BTF pool over only
later BTF blocks, falling back to the serial checked continuation when the pool
is not applicable. Scaled prefix-current rejects can also use the conservative
serial suffix restart inside the rejected block when the repaired block
preserves a validated non-empty live prefix. The preserved-column refresh used
by masked ETree-tail repair now consumes the same input-row `Rs` state, so
scaled exact-tail attempts can make omitted independent columns current instead
of immediately widening to block rebuild. Root rejected pivots can now use that
same scaled refresh to derive the shorter reject-only ETree closure instead of
being pinned to the prefix-current full suffix. Scaled all-refresh KLU-refactor
rejects now recompute the row scale vector back to input-row order before
attempting the same block-local repair and validated non-root serial suffix
restart. The newer row-first rejected-block rebuild consumes that same
input-row scale state and permutes `Rs` back to pivot order after accepting the
KLS-owned dynamic-column-pivot block, so scaled checked rejects are no longer
limited to the serial KLU-compatible block repair before the KLU fallback. If
that local repair ladder still fails, scaled rejects can now attempt the same
quality-checked KLS-first whole-numeric rebuild as unscaled rejects.
The threaded BTF worker pool now keeps a per-run completed-block bitmap and
reports prefix-current only when every diagonal block before a checked pivot
reject has finished. That lets the existing block repair and later-block
checked continuation consume more safe parallel fast-factor rejects without
claiming the full CKTSO tail algorithm.
Checked EGraph refactors now also retain the completed-column bitmap for
barriered cluster-only runs, not only for all-pipeline or pipeline-tail
schedules. A cluster-mode pivot reject can therefore be classified as
prefix-current when every earlier factor-order column is proven finished,
while same-level out-of-order cases remain conservative.
Strict tail-restart readiness now also validates that the repaired block
preserved the old prefix pivot order and that KLS can reconstruct the live KLU
prefix state a pivoting tail kernel would need, including finalized-L row
unfinalization, live `P`/`Pinv`, and symmetric-pruning `Lpend` boundaries.
KLS now uses that proof to run a conservative serial suffix restart with
pivoting for unscaled and scaled prefix-current/all-current
rejected blocks with a non-empty reusable prefix before falling back to full
block repair. Accepted serial suffix restarts now refresh only the off-diagonal
column suffix whose inverse row permutation may have changed, including the
scaled all-current subset where KLS first recomputes the row-scale vector back
to input-row order before the final pivot-order permutation.
Root-of-block rejects can now use the same KLS-owned pivoted block kernel, but
they are deliberately counted as full KLS block restarts rather than
serial-tail-ready events because no contiguous prefix can be reused. This is
still not CKTSO's full pipelined ETree-descendant tail factorization.
SPRAL static-pivot matches that are already close to complete can now be
finished with KLS's existing nonzero structural augmenting-path graph before
the row permutation is accepted. The SPRAL scaling vectors are retained only
when SPRAL itself found a complete weighted match; structurally augmented
matches fall back to KLS's matching-equilibration path so the scales remain
consistent with the final diagonal choice.

The remaining CKTSO gap is large enough that it should be treated as a missing
major algorithm, not an ordering-package tuning problem. On the selected large
`pre2` case, CKTSO completed the factor-only comparison inside the 120s cap
with about 28.4s cycle time, while KLS timed out under AMD/BTF, METIS/no-BTF,
SCOTCH/BTF, and auto policy variants. The CKTSO ordering supplement reports
`pre2` operation counts for CKTSO nested dissection and METIS in the same
range. After switching the very large pre-static MC64-adjacent path from exact
SPRAL Hungarian matching to SPRAL auction matching/scaling, `pre2` still times
out under the same 120s cap. Forced zero pivot tolerance fails as singular, and
`1e-5`/`1e-4` pivot tolerances still time out. This makes full MC64-quality
matching/scaling worth keeping, but not sufficient as the next expected gap
closer by itself. The current evidence points most strongly at CKTSO's
KLS-owned row/up-looking numeric factorization, EGraph fast factor with pivot
checks, and ETree-descendant pipelined tail restart machinery.
After rereading the local papers in `refs/`, KLS now also uses METIS's
`METIS_NodeNDP` entry point for threaded METIS analyses of matrices with at
least 30,000 rows. This asks vendored METIS for at least `log2(threads)` nested
dissection levels, matching SubtreeLU's minimum-depth separator-tree setup more
closely than plain `METIS_NodeND`, while retaining the existing CAMD refinement
and auto-ordering policy. A focused 12-row CKTSO-gap run showed this is neutral
as a standalone ordering change: 3.1772s geomean versus the saved 3.1671s
baseline, with no failures. KLS now retains the accepted `NodeNDP` top-level
component sizes as private leaf domains plus pipeline separator components, so
this is no longer only an ordering precursor. The experimental row-refactor
ready queue and KLS-first row-up-looking pivot path can consume those retained
components even when the separator map covers the dominant BTF block rather
than the full matrix, provided the accepted symbolic range can be matched
uniquely. This is still not evidence that ordering alone closes the gap because
the retained separator queues are not yet consumed by a checked-tail numeric
kernel.
The CKTSO paper in `refs/` is explicit that CKTSO's core factorization is a
row-major sparse up-looking factorization, and that the fast path combines
guessed EGraph pivot-checked refactorization with an ETree-scheduled pipelined
tail factorization when repivoting is needed. KLS no longer always hands large
eligible first factorizations to KLU: the conservative automatic `kls_first`
path can run a row-up-looking first factor and now seeds row-refactor mirrors
directly from those row entries. The remaining gap is that accepted factors are
still packed into a KLU-compatible numeric object for fallback coherence, and
KLS does not yet have CKTSO's fully row-major parallel first-factor/EGraph
executor. A debug trace on `pre2` confirmed the timeout occurs in
`trilinos_klu_l_kernel` from the fallback first-factor call. The same trace
showed SPRAL auction static pivoting finds 633565 weighted matches for 659033
rows; structurally augmenting that to a full
permutation made the static-pivoted AMD symbolic estimate worse than the
baseline, so it did not expose a credible low-risk static-pivot-only fix for
this slow case. A broader nested-dissection retry on that static-pivoted
dominant-BTF shape spent the probe timeout inside ordering analysis, so it was
not retained as a default policy.
The KLS-owned row-up first-factor scaffold now retains per-column entry counts
while generating each block, so packing no longer rescans all generated `L` and
`U` entries just to compute KLU column lengths. This is useful staging work for
native row/segment storage: on sampled KLS-first probes it moved `nxp1` initial
factor time to about 3.50s and `G2_circuit` to about 1.97s. The same scaffold
now also reserves its generated row-entry buffers from KLU's symbolic block
fill estimate before starting numeric updates, avoiding repeated large
realloc/copy waves on high-fill blocks. A focused rerun moved `nxp1` initial
factor time further to about 2.93-3.01s; `G2_circuit` stayed in the same rough
range. It did not close the CKTSO paper gap. `pre2` still times out under a
120s factor-only cap with `KLS_ENABLE_KLS_FIRST_FACTOR=1`, and `ASIC_320k`
repeated refactor remains on the column EGraph path because the current
row-group work model is still higher than the exact EGraph work. Therefore the
next direct CKTSO-aligned step is still a production row/segment numeric
engine, not simply enabling the current row-up scaffold by default.
The KLS-first pivoted fallback can now factor independent BTF diagonal blocks
concurrently when the row-up scaffold rejects. Each worker uses private
KLU-kernel scatter/workspace, a private precomputed `Offp` seed to avoid
boundary races between adjacent BTF blocks, and commits only block-local
`LUbx`/pivot/statistics into the KLS numeric object. Final off-diagonal entries
are still rebuilt from the accepted global pivot order. This closes the safe
BTF-level "factorization with pivoting over independent work" gap for the
current KLU-compatible numeric object, but it is still not CKTSO's
single-large-block ETree-descendant pipelined tail factorization.

The retained broader SPRAL post-factor trial is deliberately value-gated. A
plain broad gate improved several MC64-sensitive cases but regressed the medium
corpus because rejected or small accepted trials added setup cost. The retained
gate requires dense off-diagonal pivot evidence and accepts only when the SPRAL
candidate removes substantial pivoting pressure or materially reduces
factorization work/fill. On the 93-matrix medium paper corpus this moved KLS
geomean from 0.32198s to 0.31410s with the same three known failures. The new
exact-match wins were `hvdc1`, `OPF_10000`, `LeGresley_87936`, `rajat22`,
`rajat23`, `rajat24`, `mult_dcop_03`, and `TSOPF_FS_b39_c19`.
Benchmark artifacts now report whether METIS, SCOTCH, and SPRAL scaling were
compiled into the tested binary, and the paper-gap benchmark runner can require
SPRAL scaling explicitly. This matters for interpreting CKTSO/NICSLU-style
static-pivoting evidence: a no-SPRAL build on `rajat24` falls back to the
slow `metis`/KLU-scaling path with thousands of off-diagonal pivots, while the
SPRAL-enabled build selects the low-fill AMD unscaled static-match path with
20 off-diagonal pivots. The remaining top-row losses after that correction are
still dominated by repeated refactor throughput, not by the missing matching
hook.
A refreshed 12-row CKTSO-gap run against the saved 93-row CKTSO medium-paper
artifact makes the same point on the current code path. With `--kls-first-factor
on` and row solves enabled, the common-row KLS/CKTSO SPICE-cycle geomean ratio
is 2.4855, with all 12 rows still losing. The largest ratios are
`ASIC_320ks` at 3.267, `ASIC_320k` at 3.225, `G2_circuit` at 2.810, and
`ASIC_100ks` at 2.800. These rows report `last_factor_path=kls_fast_refactor`;
most use `initial_factor_path=kls_first`, while several use the pre-static
first-factor path. Their cycle time is dominated by repeated numeric refactor
time, not triangular solve time or first-factor fallback. On the sampled hard
rows, ordering sweeps left `auto` already choosing the best available AMD,
METIS, or SCOTCH candidate, and the row-refactor auto gate correctly stayed off
because the current row-group work estimate exceeded the exact EGraph work.
That combination makes the missing piece a numeric/storage algorithm, not an
untried ordering package or a one-matrix policy rule.
A newer top-10 CKTSO-gap run on the default production path
(`kls_current_gap10_t4_r3_timeout120.jsonl`) sharpens the same diagnosis after
the compact row-panel work. The largest ratios are `ASIC_320k` at 3.013,
`ASIC_320ks` at 2.997, `gemat12` at 2.491, `G2_circuit` at 2.477,
`ASIC_100ks` at 2.452, `transient` at 2.337, `rajat28` at 2.233,
`onetone2` at 2.223, `onetone1` at 2.169, and `rajat24` at 2.155. The updated
gap decomposition labels nine of those rows as
`column_egraph_refactor_missing_row_engine`: they use
`last_factor_path=kls_fast_refactor`, spend about 96.5%-100% of measured
EGraph work in the pipeline tail, and report zero row-refactor group work and
zero compact-panel work. `gemat12` is labelled
`klu_first_factor_missing_row_engine` because its loss is dominated by the
prestatic KLU first-factor path. Forcing the current experimental row refactor
does not convert this into a flag-selection issue: `G2_circuit` can build the
row groups but its forced row refactor is slower than the default EGraph path,
and after the singleton-BTF fix `ASIC_100ks` and `onetone2` also build row
groups but remain slower than the default path. This matches
the local CKTSO paper's stated distinction: CKTSO's fast path is a row-major
sparse up-looking numeric engine scheduled from a guessed EGraph, with row-wise
pivot checks and an ETree-descendant pipelined pivoting tail if the guess fails.
KLS's default fast path still updates KLU-compatible column storage through an
exact no-pivot EGraph schedule. The clear paper-backed missing part for these
slow rows is therefore a production row/segment-oriented up-looking numeric
engine, followed by the CKTSO/SubtreeLU pivoting-tail and separator/private-
pipeline machinery; compact row-panel kernels layered onto the current
experimental row mirror are only a precursor.
The row mirror scaffold now also handles singleton-heavy BTF partitions when
building forced row-refactor and row-solve patterns. The previous failure mode
was a KLU-storage boundary issue: singleton BTF blocks have no per-column
`L`/`U` pointer slices in the retained refactor cache, so row-pattern builders
must treat their in-block lengths as zero instead of reading raw numeric
length arrays. With that fixed, forced row refactor builds and runs on
the sampled slow rows. The row task builder now also coalesces long runs of
same-level independent scalar rows into bounded row-batch tasks, preserving
existing U-chain row segments and small-width parallel ready-queue cases. That
reduces forced-row task count on `ASIC_100ks` from 83090 groups to 54297, and
on `onetone2` from 27994 groups to 24409. It is still slower than the current
default EGraph/column path on both sampled rows, and the auto work gate
correctly leaves it off because row-group work exceeds exact EGraph dependency
work. This confirms the row/up-looking scaffold now covers and coarsens the
BTF structures that blocked it before, but it does not change the main paper
diagnosis: the gap needs a production row/segment numeric engine, not a flag
flip.

The fragmented-BTF scale policy improved the large recon artifact geomean over
the preceding EGraph build from 30.58s to 27.37s on the five completed common
rows. A refreshed eight-row selected-large reconstruction at `e610226`, after
the scaled single-block EGraph specialization, scored KLS about 1.28x slower
geomean than CKTSO with 120s timeouts penalized as 1000s, while scoring about
2.88x faster geomean than the saved KLU2 artifact under the same failure
penalty. KLS now wins this selected-large set on `TSOPF_FS_b39_c30` and
`rajat29`, ties the shared `Hamrle3` timeout, and still loses materially on
`pre2`, `nxp1`, `G3_circuit`, `rajat30`, and `ASIC_680k`. KLS still times out
on `pre2` and `Hamrle3`; CKTSO completes `pre2` and times out on `Hamrle3`.
This remains far too large to explain as a missing ordering package alone.

The remaining worthwhile ideas are not per-matrix tuning knobs. They require
new KLS-owned symbolic/numeric machinery:

- Production-scale MC64-equivalent maximum-weight matching with dual
  row/column scaling. KLS now derives dual row/column scaling from its
  KLS-owned exact sparse assignment path for small accepted candidates, builds
  BSD-licensed SPRAL Hungarian/auction matching from a pinned submodule by
  default, or can link to a system SPRAL install, and uses SPRAL for bounded
  pre-factor large weak-diagonal trials plus a value-gated post-factor trial
  for dense high-off-diagonal-pivot cases. This is still not a full production
  MC64-equivalent preprocessing stage. Existing
  MC64-style code can be reused only when its license is LGPL-compatible,
  permits source and binary redistribution with KLS, and allows preservation of
  upstream notices in KLS's third-party notice file. HSL MC64 and
  non-redistributable MC64 copies are out of scope for vendoring.
- A full intra-block parallel factor/refactor scheduler that consumes retained
  EGraph/ETree or separator-tree metadata with pivoting tail restart, beyond
  the current guarded no-pivot EGraph cluster/pipeline refactor.
- Full CKTSO-style fast factorization tail restart after a failed pivot check,
  beyond the current KLS-owned block-local restart and serial
  prefix-current/all-current tail subset.
- A production KLS-owned first-factor engine with row/segment storage, rather
  than the current env-gated KLU-compatible scaffold.
- SubtreeLU-style private/pipeline pivoting factorization from a retained
  separator tree. The current threaded METIS `NodeNDP` call now preserves the
  accepted top-level separator component sequence and private/pipeline row
  split, the experimental row-refactor ready queue can use that map for
  separator-private initial thread queues, and KLS-first row-up factorization
  now first tries the same Algorithm 6-style separator-tree split used by the
  no-pivot row-refactor path: dominant subtrees are collapsed into pipeline
  roots, child subtrees become private-thread candidates, and the candidate
  subtrees are greedily assigned by structural row-input work. The KLS-first
  consumer validates the resulting private ownership against the original
  block row dependencies before remapping; if validation fails it falls back to
  the older retained-component private/pipeline queue. When accepted, it remaps
  the block-local row/column order and `Pnum` to that queue, runs validated
  private rows concurrently, consumes pipeline rows with an Algorithm 3 atomic
  counter, and pre-updates pipeline rows from already-published private
  predecessors and earlier published pipeline-prefix predecessors before the
  ordered publish step. The pipeline row updater now records conservative
  row-supernode
  membership from published row-major `U` patterns and consumes consecutive
  ready predecessor rows, or the already-finished prefix of such a run, through
  one guarded supernode-run executor. That
  executor preserves the triangular in-run discovery rule: applying row `k`
  may create row `k+1` as the next ready dependency before the run continues,
  and it now validates the dense prefix/common-trailing shape so the trailing
  contribution is accumulated once before scatter. This covers Algorithm 4's
  ready-supernode update branch at the scalar compact-kernel level without
  claiming production BLAS-backed panel storage. When a
  safely pre-updated pipeline row needs a scoped dynamic pivot, KLS now
  publishes that pivot row under the ordered pipeline lock, updates both
  committed and phase-local prefix U entries for the column exchange, advances
  a column-order epoch, and makes speculative suffix rows that began under an
  older epoch discard and retry inside the same guarded pipeline phase. KLS
  still lacks production coarse BLAS supernode storage and a checked-tail
  factor/refactor queue consumer.
- Broader supernodal row/segment updates in the sparse up-looking executor.
  KLS has exact-pattern, ragged single-producer, and default structural and
  work-gated fragmented multi-producer dense-panel updates, but these are
  still narrower than SubtreeLU's production supernodal coverage.

Several smaller dispatch experiments were tried and rejected because they helped
some benchmark cases while regressing others. Those are documented below so the
project does not drift toward benchmark-name-specific heuristics.

## Paper-by-Paper Coverage

| Reference | Implemented in KLS | Partial or open coverage |
| --- | --- | --- |
| Algorithm 907 / KLU | BTF preprocessing, fill-reducing ordering, row scaling modes, Gilbert-Peierls factorization with partial pivoting, no-pivot refactorization, and block back substitution are present through the vendored KLU-derived kernel. KLS adds automatic policy selection, serial refactor scatter metadata, and exact EGraph level metadata around these pieces. | KLS still inherits KLU's fundamentally sequential intra-block numeric kernel. |
| NICSLU | AMD-style ordering, optional static-pivoting preprocessing, optional SPRAL Hungarian/scaling trials, and the idea that parallel kernels should be selected by general structural/numeric evidence are represented in KLS policies. KLS now records exact no-pivot EGraph levels from the numeric U pattern and uses them in guarded large single-block, dominant-BTF-block, and fragmented non-dominant many-block refactor paths, including KLU row-scaled cases where scale handling is supported and work-estimated cluster-level thread slices. KLS also reports the NICSLU R1/R2 static parallel suitability model as `parallel_model_r1`, `parallel_model_r2`, and `parallel_model_recommends_parallel`, using the paper's 2.0 and 50.0 thresholds. After numeric factorization, KLS also evaluates NICSLU Algorithm 4's task-flow earliest-finish model on the actual U-dependency graph using `2*nnz(L(:,i))` update work, `nnz(L(:,k))` normalization work, and a unit dependency sync cost; benchmark artifacts expose the resulting work, finish time, speedup, dependency count, and recommendation. These models seed KLS-owned row/segment refactor preparation when an exact dependency schedule and the row-work gate agree. | Full production MC64 matching/scaling is not implemented. NICSLU's detailed ETree/EScheduler-guided intra-block factorization and full pivoting-aware ETree scheduling are not implemented. Earlier broader EGraph prototypes were rejected because they were not general wins on the current kernel/storage. |
| CKTSO | METIS nested-dissection ordering, guarded SCOTCH nested-dissection auto trials for large high-work symbolic candidates, constrained-minimum-degree-style CAMD refinement, combined ordering selection, pivot-checked fast factorization, default CKTSO-style row-wise guessed-diagonal checks in KLS-owned checked row fast/refactor passes, KLS-owned block-local restart after a failed fast-factor pivot check including root-of-block rejects, conservative serial prefix-current/all-current tail restart for validated non-root unscaled repaired blocks, range-aware row-first repair over retained topological pivoting-tail envelopes including non-contiguous active masks and compact retained ETree-tail row worklists, scaled block-local restart plus threaded checked continuation over later BTF blocks, scaled prefix-current/all-current block repair, scaled in-block serial tail restart after recomputing row scales to input-row order, threaded BTF worker-pool completed-block tracking for safe prefix-current rejects, guarded EGraph cluster/pipeline no-pivot refactors for single, dominant BTF, and selected fragmented many-block BTF shapes, work-balanced cluster-level refactor slices, cached row-permutation solve scratch, CKTSO Section V-style structure-adaptive triangular solve metadata/executor for normal and transpose single-RHS solves, and dual-potential plus optional SPRAL matching-derived equilibration trials are implemented in KLS at the KLU-wrapper layer. | CKTSO's maximum-weight matching with dual scaling is still only partially approximated because KLS does not have an always-on production MC64-equivalent weighted assignment stage. The retained ETree-tail executor exists inside KLS-owned block repair, but the full CKTSO clustered/pipelined fast-factor scheduler still is not implemented across the whole guessed-EGraph interruption path. Otherwise unsupported fast-factor failures still fall back to full pivoting factorization. |
| SubtreeLU | KLS vendors reproducible METIS/GKlib and SCOTCH submodules, uses METIS plus CAMD refinement, asks METIS `NodeNDP` for at least `log2(threads)` nested-dissection levels on larger threaded METIS analyses, and can keep SCOTCH from `auto` when its symbolic score is materially better on large high-work cases. KLS now retains accepted `NodeNDP` component sequences from METIS user-order callbacks, stitches them into a global BTF-aware separator forest with synthetic private components for blocks that did not run `NodeNDP`, reports the resulting global queue shape in stats/bench output, and uses the component map to build separator-private initial thread queues for the experimental row-refactor ready queue. No-pivot and checked row refactors can now consume the retained separator tree through a SubtreeLU Algorithm 6-style FLOP-balanced private/pipeline queue, with separator-crossing row groups forced into pipeline work. The KLS-first pivoting row-up-looking factorization now tries an Algorithm 6-style separator-tree split first: dominant subtrees become pipeline roots, child subtrees remain private candidates, candidates are assigned to private threads by block-local row-input work, and dependency validation falls back to the older retained-component queue if private ownership is unsafe. Accepted queues are remapped into the block-local row/column order and `Pnum`, execute private rows in worker-local entries, and consume pipeline rows through a guarded Algorithm 3 atomic counter. KLS also uses SubtreeLU Algorithm 4's scoped `N'` pivot maximum inside the current collapsed component extent, rejects unsafe cross-domain exchanges, applies ready predecessor row-supernode runs and ready prefixes during partial pipeline-row updates with in-run triangular discovery and a compact common-trailing accumulation, publishes compact-validated ready panels before the first consumer falls back to scalar/compact walking, and handles safely pre-updated weak-pivot pipeline rows with an ordered pivot publish plus column-order epoch retry for speculative suffix rows inside the same pipeline phase. This matches SubtreeLU's separator-domain pivot rule while preserving KLS's fallback to the pivoted block kernel for unsafe scoped rows. KLS also records row-major U-pattern supernode candidate diagnostics from the exact no-pivot refactor dependency pass, has scalar compact-panel producer and consumer updates, has KLS-owned scalar batched producer-to-consumer-row-subrange updates for exact multi-producer patterns, ragged single-producer suffix patterns in dense and independent row groups, and default structural/work-gated contiguous and fragmented multi-producer row-panel updates for dense producer suffixes. Optional CBLAS experiments remain separate for completed-supernode row updates and unchecked blocked producer-panel `dtrsm`/`dgemm`. | KLS still does not use production SubtreeLU-style coarse supernodes/BLAS updates broadly enough for the paper slow cases, and the native row/segment numeric engine remains a scaffold around KLU-compatible packing. CKTSO's full checked-tail/pivoting executor also remains open. |

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
  and large dominant-block BTF analyses when the stripped fringe is small; it
  also covers many-block analyses whose large block leaves an inflated symbolic
  estimate and whose no-BTF retry cuts the symbolic score at least in half,
  matching KLU's warning that BTF can occasionally increase factor fill.
  Dominant/inflated METIS-start retries now require a known BTF symbolic score;
  otherwise KLS preserves BTF, because an unknown BTF estimate can make a bad
  no-BTF single-block symbolic look artificially preferable.
  Low-work dominant-BTF cases are excluded from this retry because keeping BTF
  is cheaper when the symbolic work estimate is already low. The METIS promotion
  gate also
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
  TSOPF-style power-grid structures. The medium bounded-degree METIS start is
  limited to near-full diagonals, because Rommes/BIPS-style rows with a few
  percent missing diagonal entries measured faster on AMD despite fitting the
  same low-degree envelope. For large high-work no-BTF single-block analyses
  outside those direct-start classes, auto can also try METIS
  symbolically before the first numeric factorization and keep it when the
  symbolic fill score is clearly lower, avoiding a delayed post-factor METIS
  promotion. This preserves the papers' nested-dissection motivation without
  naming individual matrices.
- Constrained nested-dissection refinement: KLS can refine METIS rank groups
  with CAMD constraints, preserving nested-dissection rank shape while reducing
  local fill and flops. Large METIS orderings now use coarse rank-group CAMD
  refinement by default instead of limiting the CKTSO-style refinement to one
  medium structural class.
  Threaded METIS analyses for matrices with at least 30,000 rows now use
  `METIS_NodeNDP`, which enforces the top nested-dissection levels needed for a
  thread-count-sized separator-tree skeleton before the same CAMD refinement.
  Accepted METIS analyses also retain the `NodeNDP` top-level component sequence
  as private leaf domains plus pipeline separator components, including an
  ordering-position to component map. The row-refactor ready queue can consume
  that map to assign its first ready groups by separator component before the
  shared/local successor queue takes over.
- Scaling policy: KLS exposes KLU scale modes and has auto scale selection based
  on pattern and numeric evidence, including no-scale reuse where repeated
  SPICE refactorization benefits. Low-work dominant-BTF cases also start with
  no KLU row scaling, avoiding the scale recomputation cost while preserving
  the useful BTF decomposition. Large high-work METIS/no-BTF single-block paths
  that already selected max scaling skip redundant post-factor scale trial
  factorizations. Medium and large TSOPF-style spiked low-diagonal METIS starts
  begin with sum scaling and a lower `1e-4` pivot tolerance once their
  structural dominant-BTF shape is known, avoiding the earlier
  max-scale/default-tolerance discovery factorizations. Small spiked
  low-diagonal METIS starts use unscaled `0` mode with the same lower pivot
  tolerance and skip static-pivot trials that do not improve that structural
  class.
- Fast repeated factorization with pivot check: KLS reuses an existing numeric
  pattern, checks reused pivots against the selected threshold via L
  multipliers, and falls back to full pivoting factorization if the reused order
  is unsafe and cannot be repaired locally. For unscaled BTF patterns, KLS can
  restart the rejected diagonal block with pivoting, rebuild the global row
  permutation/off-diagonal entries around that block repair, and retry the
  checked fast factorization. KLS records the first rejected factor-order pivot,
  original matrix column, rejected block suffix, exact U-pattern descendant
  tail, ordered-block ETree successor path, sorted pivoting-tail worklist scope,
  and number of block restarts so future tail-restart work can distinguish
  late-tail failures from early failures and compare no-pivot versus
  pivoting-tail recomputation scopes.
- Static pivoting trial: KLS has value-aware greedy row matching, layered
  augmenting-path search for larger weak-diagonal candidates, and swap
  improvement for weak or high-off-diagonal-pivot medium matrices. Medium
  static-pivot candidates also run a bounded alternating-cycle pass that can
  apply profitable three- and four-row exchanges missed by pair swaps, while
  larger candidates keep the lower-fill layered/swap pattern. Small exact
  assignment candidates now also return assignment-dual row/column scaling, and
  KLS can still trial dual-potential matching-derived row/column equilibration
  for non-exact matches. It keeps the
  transformed candidate only when numeric quality and cost evidence justify it.
  With default SPRAL scaling enabled, large weak-diagonal candidates whose BTF
  analysis leaves one dominant block can run BSD-licensed SPRAL Hungarian
  matching/scaling before the first factorization; many-block BTF cases are
  excluded because the existing BTF path is already efficient there. After
  cheaper scale, ordering, and pivot-tolerance trials have run, KLS can also
  run a SPRAL Hungarian/scaling trial for dense high-off-diagonal-pivot cases.
  The accepted candidate must pass the normal numeric comparison plus a
  matching-specific value gate: it must remove substantial pivoting pressure or
  materially reduce factor work/fill, which rejects small cases where matching
  setup costs more than the repeated-refactor saving.
  The pre-factor static-pivoting gate also covers moderately sized matrices
  whose input values show a majority of weak or missing diagonal entries, plus
  medium-large mostly diagonal matrices with thousands of weak diagonal rows.
  These are general numeric-structure rules used by the frequency-domain and
  Rajat-family paper cases without naming individual benchmarks.
- Initial KLS-owned parallelism: repeated refactorization can run across
  independent BTF diagonal blocks for large high-flop cases where the block work
  is wide enough to offset thread overhead. The threaded path keeps a persistent
  worker pool on the solver instance so repeated SPICE refactors reuse workers
  and scratch storage instead of relaunching threads each cycle. Checked worker
  runs mark completed diagonal blocks and classify a pivot reject as
  prefix-current only when all earlier blocks completed.
- KLS-owned serial refactor metadata: for unscaled numeric patterns that are
  not handled by the threaded BTF worker pool, KLS precomputes the fixed
  scatter from factor-order columns to pivotal rows and input value positions.
  This removes repeated `Q`/`Pinv` structure lookups from SPICE refactor cycles
  while preserving the existing KLU-derived LU storage. It covers both
  single-block and serial BTF refactors. For serial BTF refactors, KLS also
  prepartitions each mapped column into off-block and diagonal-block entries, so
  repeated refactors do not reclassify the same BTF structure in the numeric
  scatter loop. A narrow scaled many-fringe dominant-BTF subset uses the same
  serial map after recomputing row scales and then permuting the scale vector
  back to pivot order.
- KLS-owned EGraph metadata and guarded refactor consumption: for threaded
  numeric runs, KLS can levelize the exact no-pivot refactor dependency graph
  from the actual U pattern after factorization. It retains level pointers and
  column lists and reports the level count, maximum level width, and dependency
  edge count in `kls_stats` and benchmark JSON. It also reports the
  CKTSO-style cluster/pipeline split point, pipeline-column count, approximate
  no-pivot update work, and tail work implied by the current level widths.
  Large single-block, high-flop dominant-BTF, and very fragmented unscaled
  many-block matrices with enough dependency work and level width can now
  consume this schedule through a work-estimated per-thread cluster-mode refactor,
  then switch to a no-pivot pipeline tail where each worker claims tail columns
  from an atomic cursor and waits only for actual U-pattern predecessors. The
  same block-aware EGraph kernel can run inside a large dominant BTF block and
  update Offx for entries above that block. The EGraph consumer now keeps a
  solver-owned worker pool and dense worker scratch across repeated refactors,
  plus generation-stamped solver-owned pipeline dependency markers, matching
  the CKTSO/NICSLU emphasis on retained scheduling state instead of relaunching
  threads or clearing a fresh done array on every SPICE step. Schedule
  construction is now
  limited to single-block, dominant-BTF, or high-work fragmented non-dominant
  many-block shapes with enough numeric or measured dependency work to consume
  it, while ordinary non-dominant many-block BTF and low-work dominant-BTF cases
  without enough dependency work skip the setup and stay on the BTF worker pool
  or mapped refactor paths.
  For the high-coverage many-fringe dominant-BTF class below the normal EGraph
  size floor, KLS can run the whole exact EGraph as an atomic topological
  pipeline with no cluster barriers. This is a retained SubtreeLU/CKTSO-aligned
  scheduler improvement for the current fixed-pivot LU storage.
  Very high-work single-block row-scaled factors can also use that full
  all-pipeline EGraph path once the existing scale recomputation/permutation
  support is available, so the scheduler does not leave a tiny barriered tail
  after hundreds of narrow cluster levels. Those scaled single-block EGraph
  refactors now use a dedicated column kernel that applies KLU row scales
  directly while loading the fixed input-position map, instead of paying the
  generic BTF-capable value-loader branch on every entry.
  A later fragmented many-block gate applies the same CKTSO-style lesson to
  unscaled BTF decompositions whose largest block covers less than half the
  matrix but still carries enough numeric work to justify intra-block
  scheduling. On the 93-matrix medium paper corpus this changed only
  `ASIC_680ks` schedule counters, cutting its cycle from 17.62s to 9.42s and
  moving KLS geomean from 0.28594s to 0.28392s with the same three known
  failures.
  A later retained low-work dominant-BTF gate lowers the EGraph schedule floor
  when measured dependency work is still material. On the 93-matrix medium
  paper corpus this moved KLS geomean from 0.32749s to 0.32198s with the same
  three known failures, mainly by cutting IBM `dc1/dc2/dc3/trans4/trans5`
  cycles by about 22-30% and `scircuit` by about 20%.
  This is still narrower than CKTSO's production pivoting machinery, but it is
  the first retained intra-block EGraph cluster/pipeline refactor path.
- SPICE-cycle orientation policy: KLS can analyze normal and transposed storage
  orientations and select the faster internal form for repeated solve cycles.
- LGPL project licensing and third-party notices: KLS itself is
  LGPL-2.1-or-later, with vendored-source attribution separated from the KLS
  license. MC64-equivalent preprocessing must stay inside that licensing
  boundary: HSL MC64 and solver-tree copies that retain HSL redistribution
  restrictions are not vendorable, while KLS can build the BSD-licensed SPRAL
  scaling subset from `third_party/spral`, use another LGPL-compatible
  redistributable MC64-style source with preserved notices, or use independent
  KLS code. The Rust `rwl/mc64` package is also BSD-licensed and useful as a
  reference, but it is a partial SPRAL translation and is not a better fit than
  the pinned SPRAL submodule for KLS's C/Fortran build.
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
  and permuted afterward. Failed checked passes now seed the retained pivoting
  tail from the interrupted guessed-EGraph unfinished set before prefix refresh,
  and KLS-owned block repair can consume that retained ETree tail as a compact
  topological row worklist. KLS still does not have CKTSO's full
  guessed-EGraph interruption scheduler around that executor.
- CKTSO fast factorization is present as pivot-checked reuse plus a KLS-owned
  BTF-block repair path. KLS validates the preserved-prefix live state and can
  enter a conservative serial pivoting-tail kernel for non-root rejects in that
  subset; threaded block repairs can consume retained ETree-tail rows through a
  compact topological row-pipeline worklist; root rejects can use the same
  KLS-owned pivoted kernel as full block restarts. KLS still does not implement
  the complete CKTSO fast-factor scheduler around that tail executor.
- SubtreeLU-style nested-dissection metadata is now retained from accepted
  METIS `NodeNDP` analyses as private/pipeline component queues, and the
  experimental row-refactor ready queue can consume the full-factor map for
  separator-private initial queues. KLS-first pivoting row-up factorization
  now consumes that map through a validated Algorithm 6-style
  private/pipeline splitter, but checked-tail factor/refactor work still does
  not use the retained queues.

## Not Implemented Yet

- Full production MC64-equivalent weighted matching and assignment-dual
  row/column scaling across the broad matrix set.
- NICSLU's ETree/EScheduler executor beyond the implemented R1/R2 and
  Algorithm 4 task-flow suitability counters.
- Intra-block parallel factorization with pivoting scheduled by an ETree.
- General ETree-based pivoting tail restart beyond the current guarded no-pivot
  EGraph cluster/pipeline path, KLS-owned retained-tail block repair, and serial
  prefix-current block-tail subset.
- CKTSO dual-mode cluster/pipeline fast factorization with pivot check.
- CKTSO's complete guessed-EGraph interruption scheduler around pipelined
  ETree-descendant tail factorization with pivoting after a pivot-check failure.
- SubtreeLU separator-tree collapse for checked-tail and refactor pivoting
  kernels outside the KLS-first row-up path.
- Broader/default SubtreeLU FLOP-balanced separator-tree partitioning for
  refactorization across BTF forests and checked-tail kernels.
- SubtreeLU constrained pivot search within nested-dissection subdomains.
- Production SubtreeLU supernodal updates with coarse BLAS kernels inside the
  sparse up-looking framework.

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
After the EGraph path was narrowed and retained, the pipeline tail scheduler
was changed from fixed per-thread strides to an atomic dynamic work cursor,
matching CKTSO's on-the-fly tail assignment more closely while preserving KLS's
existing fixed-pivot column kernel. On the nine-row hard focus set, same-session
4-thread runs moved from about 9.07s geomean for the saved KLS artifact to
about 8.91-9.04s. A 90-common-row paper-medium run kept the same failure set as
the recent KLS medium baseline and improved geomean on common rows versus that
baseline, while KLS still trailed the CKTSO artifact by about 1.29x geomean on
the 90 common medium rows.

An EGraph-specific persistent worker-pool prototype was then tested to avoid
recreating pthreads and scratch arrays on every no-pivot EGraph refactor. The
idea matched CKTSO's emphasis on retained scheduling machinery, but
same-session A/B checks against commit `f36bd4f` did not show a general win:
`ASIC_680k` moved only from about 0.06887s to 0.06784s average refactor time,
while `rajat30` moved from about 0.30470s to 0.30623s and `nxp1` from about
0.64476s to 0.65076s. The prototype was removed because the current EGraph
runtime is dominated by numeric scatter/update work rather than pthread launch
overhead on the large CKTSO-gap rows.

Two later EGraph scheduler tweaks were also rejected on the same basis. A
dynamic work-claim path inside wide cluster levels was neutral on `ASIC_100ks`
but slower on `rajat24` and `transient`. Forcing clustered cases into a full
all-pipeline EGraph schedule was also slower on `ASIC_100ks`, `rajat24`, and
`transient`. These tests indicate that KLS is not mainly missing another
cluster/pipeline dispatch tweak in the current column-storage EGraph kernel.
The CKTSO and SubtreeLU papers point instead to row-major sparse up-looking
storage and row/supernode updates as the remaining large lever.

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

The retry was later extended only for many-block BTF analyses with a large but
not overwhelming dominant block. This targets cases like the IPSO HTC matrices,
where BTF created about 29k blocks but still left an 87% dominant block and a
much larger symbolic estimate than the no-BTF structure. The retained gate
requires at least 1024 BTF blocks, at least 100k rows, a largest block covering
80-95% of the matrix, no low symbolic-work estimate, and a no-BTF symbolic score
at most half of the BTF score. The METIS-started auto path now runs this strict
retry before returning; it still skips the older single-block retry to avoid
extra symbolic work on already-good METIS/BTF single-block cases. In same-session
checks, `HTC_336_4438` moved from the saved METIS/BTF scale-1 path
(`initial=13.14s`, `refactor=0.082s`) to METIS/no-BTF no-scale
(`initial=3.52s`, `refactor=0.103s`), while `HTC_336_9129` moved from
`initial=8.37s` to `1.93s`. Guard cases kept their earlier BTF decisions:
`G2_circuit` remained METIS/BTF, `transient` and `power197k` remained AMD/BTF,
and `ASIC_680ks` remained METIS/BTF. A parallel SCOTCH sweep did not justify an
auto SCOTCH policy: it was much slower than METIS on `rajat30` and `nxp1`, and
mixed or worse on `G2_circuit`, `transient`, and `power197k`.

The direct METIS retry was then tightened after guard checks showed a harmful
interaction with the ASIC 320k family. Those matrices have a very useful BTF
decomposition, but their METIS/BTF symbolic score is unknown while the no-BTF
symbolic has a huge explicit fill estimate. Accepting that no-BTF candidate sent
`ASIC_320k` to a 27M+27M single-block symbolic and caused a timeout, while the
retained METIS/BTF path has about 2.0M+2.0M numeric fill and refactors near
0.10s. KLS now requires a known current BTF score before accepting dominant or
inflated many-block no-BTF retries, and the broad large-low-degree structural
shortcut no longer starts no-BTF blindly. Clean checks restored `ASIC_320k` and
`ASIC_320ks` to METIS/BTF, while `HTC_336_4438` still selects METIS/no-BTF and
`power197k` keeps AMD/BTF.

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
initially started with METIS while keeping BTF enabled, and left auto scaling
on the default max-scaling path rather than forcing no-scale. In same-session checks,
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
the same cap. The result did not justify broad SCOTCH auto ordering or a new
matrix-specific rule. Later auto support therefore kept SCOTCH restricted to a
guarded large single-block symbolic trial that is accepted only on a clear
fill/work score win.
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

A follow-up coarse-grouping probe reduced the large METIS CAMD group size from
4096 to 1024 as a general CKTSO-style constrained-ordering experiment. It was
rejected: `pre2` no-BTF METIS symbolic fill worsened from about 121.5M to
124.7M nonzeros, so the retained 4096 grouping remains the better large
default.

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

An LGPL-compatible medium-matrix SPRAL Hungarian-first static-pivot experiment
was also tested on the hard paper focus set. The structural gate targeted
50k-plus mostly complete diagonals with thousands of weak diagonal entries, so
majority-missing AT&T-style cases stayed on the retained greedy unscaled path.
It made `rajat28` select exact SPRAL matching, switch from KLU max scaling to
matching-derived no-scale values, and reduce off-diagonal pivots from 1 to 0,
but its refactor time worsened from about 0.162s to about 0.169s and the
9-row KLS/CKTSO focus ratio regressed from about 2.66x to about 2.72x. The
experiment was removed. This reinforces that BSD/LGPL-compatible MC64-style
matching is allowed and present through SPRAL, but widening it over medium
mostly diagonal cases is not enough to close the large CKTSO gap.

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

The exact sparse assignment path now keeps the assignment dual potentials when
the full match succeeds and converts them into row/column scaling for the
accepted candidate. For each original row matched to column `j`, KLS stores the
row scale at the permuted row `j` and the column scale at original column `j`;
the resulting scaled matched diagonal has unit magnitude, while all other
entries satisfy the dual reduced-cost bound. The smoke suite now requires the
3000-row scaled static-pivot fixture to report
`selected_exact_matching_scaling`, and benchmark JSON exports the same flag.
This fills the small exact-matching part of the NICSLU/CKTSO MC64 contract
without changing the conclusion above: larger paper cases still need an
optimized production MC64-equivalent implementation or a retained
LGPL-compatible SPRAL path, not a wider use of the straightforward exact
augmenting-path code.

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
accepted KLS row matchings solely on weight. KLS now reports accepted SPRAL
Hungarian or auction row permutations separately as
`selected_spral_matching`, because the auction path is LGPL-compatible but not
an exact assignment path and should not be conflated with
`selected_exact_matching`.

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
from about 0.15s repeated refactors to about 0.27s, so the retained gate was
narrowed toward dominant BTF structure. Later focused checks with the pinned
SPRAL submodule showed the same boundary on paper-medium cases: `power197k`
selected SPRAL matching, cut off-diagonal pivots from about 51k to about 1.8k,
and reduced repeated refactor from about 0.055s to about 0.008s, while
`HTC_336_4438` moved from the intended METIS/no-BTF path to AMD/BTF and
regressed from about 0.103s to about 0.190s. The large SPRAL pre-static gate
therefore now requires an existing BTF symbolic analysis with a dominant block,
preserving no-BTF ordering choices. `pre2` still timed out at 120s after trying
the SPRAL path; its analyze-only evidence remains a dominant 629628-row block
with about `2.08e11` estimated flops. This confirms that license-compatible
MC64-adjacent preprocessing is useful and worth keeping, but it does not
replace the missing CKTSO/SubtreeLU-style intra-block numeric/scheduling
machinery.

Because the retained SPRAL path is now guarded and materially improves a
paper-medium hard row, the pinned BSD scaling subset was promoted from an
opt-in component to the default build, while keeping
`KLS_ENABLE_SPRAL_SCALING=OFF` for C-only builds and
`KLS_USE_SYSTEM_SPRAL=ON` for system installations. Fresh default-vs-SPRAL
checks show why this is aligned with the solver goal: `power197k` moves from
the no-SPRAL default path with about 51k off-diagonal pivots and roughly 0.057s
refactors to the SPRAL static-match path with about 1.8k off-diagonal pivots
and about 0.008s refactors. The same guard keeps `HTC_336_4438` on its
METIS/no-BTF path, while `transient`, `onetone1`, `onetone2`, and `rajat28`
remain on their existing accepted paths.

The fast-factor pivot-check path was then made more diagnostic by recording the
first rejected factor-order pivot and original matrix column in `kls_stats` and
benchmark JSON. It now also records the rejected BTF block start/size, the
simple suffix length from the rejected pivot to the end of the block, and the
exact U-pattern descendant tail size/work inside that block. It now also
records an ordered-block ETree successor-path size/work estimate for the first
rejected pivot, and a sorted pivoting-tail worklist scope seeded from the
current refresh state, matching the CKTSO paper's pivoting-tail upper-bound
dependency idea more closely. This does not implement CKTSO's pipelined tail
factorization, but it is a required prerequisite: KLS can now measure whether
failed fast factorizations reject near the tail, where a pivoting-tail restart
could avoid recomputing the whole block, or near the front, where full fallback
is still expected, and can distinguish a true no-pivot dependency tail from a
broad suffix, a single ETree successor path, and the actual pivoting-tail plan.

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
retry the checked fast factorization when later columns may not have been
refreshed. When the rejected fast pass had already refreshed all columns, or
when a serial BTF pass rejected in the final block, KLS now validates the
repaired block tail directly and returns without a second full checked
refactor. Smoke tests now cover both the unscaled repair and the scaled
fallback path, and representative static-pivot cases (`gemat12`, `onetone2`,
`twotone`, and `rajat25`) did not trigger unexpected block restarts. This is
useful CKTSO-aligned infrastructure, but it is still not CKTSO's production
tail restart: it does not retain an ETree/EGraph tail or restart only
descendant work inside a large single block.

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
distribution. CMake now enforces that distinction for the system-SPRAL path by
requiring `KLS_SYSTEM_SPRAL_LGPL_COMPATIBLE=ON`; the bundled path checks the
pinned SPRAL `LICENCE` file before building its scaling subset. A small BSD
Rust `mc64` crate exists as a partial SPRAL translation, but it does not
improve KLS's C integration story over the already pinned SPRAL Fortran/C
interface. Fresh SPRAL-enabled checks also confirm the policy should stay
guarded rather than become an unconditional default:
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

The same "start with the already accepted policy" idea was applied to
TSOPF-style spiked low-diagonal matrices. Auto mode had already learned that
`TSOPF_FS_b39_c30` wanted METIS/BTF, sum scaling, and `1e-4` pivot tolerance,
but it reached that state by first paying for max-scale/default-tolerance
factorizations. The structural dominant-BTF gate now starts this class directly
with sum scaling and `1e-4` tolerance and skips the redundant scale and pivot
trials. A focused `TSOPF_FS_b39_c30` check reduced initial factor time from
about 86s to about 12s while preserving the about 2.6s refactor path and valid
residual, moving the SPICE-cycle estimate from about 370s to about 290s versus
the saved CKTSO artifact at about 333s.

The same structural idea was then extended downward to small spiked
low-diagonal TSOPF/QY cases. Across the medium manifest, the small rule matches
the `TSOPF_FS_b9_c1`, `TSOPF_FS_b9_c6`, and `case9` shapes: 2k-20k rows,
about 8-14 nonzeros per row, 20-35% diagonal coverage, and row/column degree
spikes around 40-60% of the matrix order. Focused sweeps showed these cases
prefer METIS/BTF, KLU's unscaled `0` mode, and `1e-4` pivot tolerance. The
static-pivot trial was also skipped for this class because it added setup cost
without being selected. Current auto checks moved `TSOPF_FS_b9_c1` to about
0.13s on the SPICE-cycle estimate, `TSOPF_FS_b9_c6` to about 1.9s, and
`case9` to about 1.9s with valid residuals. The saved CKTSO artifact is still
faster on `TSOPF_FS_b9_c1` at about 0.061s, but KLS is ahead on the larger
`b9_c6` and `case9` rows.

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

A CKTSO-style topological level-order all-pipeline cursor for huge single-block
EGraph refactors was also tested and rejected in favor of the current natural
column-order cursor. On same-session checks it regressed `nxp1` refactor time
to about 0.72s and `rajat30` to about 0.39s, so the missing gap is not simply
the absence of a level-ordered all-pipeline cursor. The retained improvement in
this area is narrower: scaled single-block EGraph refactors now dispatch to a
dedicated hot kernel that applies row scaling directly while loading the fixed
input-position map. Same-session probes improved `nxp1` repeated refactor from
about 0.301s to 0.287s and `rajat30` from about 0.264s to 0.218s with valid
residuals, while unscaled/BTF probes stayed valid.

Disabling the huge-single all-pipeline gate was then tested as a direct
cluster/pipeline split experiment. It regressed `nxp1` repeated refactor to
about 0.65s and `G3_circuit` to about 19.43s, so the current all-pipeline
natural cursor remains the better general choice for huge single-block factors.
A wider eight-way scalar scatter-subtract unroll was also tested and rejected:
it regressed `nxp1` repeated refactor to about 0.38s and did not provide a
general win on the quick large guards. The existing four-way generic scalar
scatter kernel is retained.

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

The dominant-BTF EGraph gate was then lowered for medium ASIC-style matrices
whose largest diagonal block is at least 90k rows, covers at least 95% of the
matrix, and has at least `5e8` actual factor flops. The EGraph consumer itself
still requires at least `2.5e8` no-pivot dependency work and enough level width,
so low-work dominant-BTF cases such as `rajat29` skip schedule construction.
This lets the retained cluster/pipeline refactor run inside the dominant block
of the `ASIC_100k` and `ASIC_320k` families while ordinary many-block
`ASIC_680k`/`ASIC_680ks` shapes remain on their existing paths. Focused checks
showed the largest wins on `ASIC_320k` and `ASIC_320ks`: repeated refactor time
moved from roughly 0.35s/0.29s to about 0.11s/0.08s with valid residuals.
`ASIC_100k` and `ASIC_100ks` improved more modestly, while the retained guard
kept `rajat29` schedule metrics at zero. This is a useful CKTSO-inspired
coverage extension, but CKTSO remains faster on these ASIC rows.

The same gate was then extended below the 90k dominant-block floor only for
high-work dominant-BTF shapes whose largest block covers at least 95% of the
matrix and whose measured factorization has at least `5e9` flops. This keeps
the rule tied to general work evidence rather than TSOPF names. On
`TSOPF_FS_b39_c19`, which has a 76215-row dominant block and about `1.54e10`
factor flops, the retained EGraph path reduced repeated factor/refactor
averages from about 4.9s/4.8s to about 1.5s/1.5s with valid residuals, cutting
the SPICE-cycle estimate to about 167s versus the saved CKTSO result near
603s. The previously rejected Rajat 80k-row class stays excluded because it
does not meet the high-work gate; focused `rajat28` checks still recorded zero
EGraph schedule metrics and stayed on the existing path.

The dominant-BTF EGraph gate was also probed down to 80k rows and 80k-row
largest blocks to see if the same policy should cover high-flop Rajat
dominant-BTF rows. It activated on `rajat20`, `rajat25`, and `rajat28`, but
slowed repeated refactors from the existing roughly 0.15s class to roughly
0.18-0.19s. The 90k retained floor remains the better general rule; the Rajat
rows need different numeric scheduling or pivoting work, not more EGraph
coverage with the current kernel.

The retained EGraph gate was later extended to a different dominant-BTF shape:
large-heavy blocks below the 95% coverage floor. The new branch requires the
largest block to cover at least 85% of the matrix, contain at least 100k rows,
have at most 20k total BTF blocks, and show at least `2e9` actual factor flops.
This activates on AT&T `twotone` but not on the previously rejected 80k-row
Rajat class. In clean checks, `twotone` built 2391 EGraph levels and moved
repeated refactor from the current worker-path class around 0.97-1.26s to about
0.34s, reducing the 100-cycle SPICE estimate from roughly 101s to roughly 38s.
Guards kept `rajat28` schedule metrics at zero, kept `ASIC_680ks` on the
many-small-block worker path, preserved `G2_circuit`, and restored
`ASIC_320k`/`ASIC_320ks` to their existing METIS/BTF EGraph path after the
no-BTF retry fix.

The same coverage-bounded idea was then extended to medium-heavy dominant BTF
blocks. The initial retained branch required 85-95% dominant-block coverage, at
least a 30k-row largest block, at most 5k BTF blocks, and at least `5e8` actual
factor flops; the upper coverage bound was deliberate so the previously
rejected 95%+ Rajat class remained excluded. In clean focused checks,
`onetone1` built 1705 EGraph levels and reduced repeated refactor from about
0.18s to about 0.06-0.065s, cutting the SPICE-cycle estimate from roughly 19s
to about 7s. `onetone2` stayed below that first work threshold with zero
schedule metrics, `rajat28` remained excluded, and `twotone` stayed on the
large-heavy EGraph path.

The high-coverage dominant-BTF schedule floor was then lowered separately from
the 85-95% medium-heavy branch. The retained rule still requires at least
`2e8` actual factor flops before building EGraph metadata and at least `1e8`
computed dependency work before the consumer runs. This activates a lower-work
95%+ dominant-block case without reopening the rejected Rajat class: `transient`
built 807 EGraph levels and reduced repeated refactor from about 0.057s to about
0.024s, cutting the focused SPICE-cycle estimate from about 7.0s to about
3.8s. At that point `onetone2` remained below the `2e8` high-coverage
factor-work floor with zero schedule metrics, `rajat28` remained excluded,
`onetone1` stayed on the medium-heavy EGraph branch, and `ASIC_320k` stayed on
the existing high-coverage path.

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

A separate worker-map policy was then retained for the opposite many-block
shape: dominant BTF matrices with a bounded block count and bounded input
size. KLS now builds the existing fixed-pivot input-position map for
pool-eligible dominant BTF cases when the largest block covers at least 75% of
the matrix, the block count is at most 20k, and the input has at most 3M
nonzeros. The worker-pool map path is scale-aware, so scaled Rajat cases can
reuse it without falling back to the original `Q`/`Pinv` scan. This is still
not the CKTSO numeric kernel, but it removes a general repeated-refactor
overhead from Rajat/AT&T-style dominant-block patterns while preserving the
previous ASIC tiny-block rejection. Focused one-pass checks improved
`twotone` from about 129s to about 101s on the SPICE-cycle estimate,
`onetone2` from about 5.3s to about 4.6s, and `rajat20`/`rajat25`/`rajat28`
from about 16.6s/18.2s/17.6s to about 14.7s/15.3s/15.5s. The large guards
remain on their intended paths: `rajat29` is excluded by the 3M-nnz gate and
stays near the 4.1s class, while `ASIC_680k` remains an extreme many-block
case. A full one-pass medium manifest with the retained policy completed 90 of
93 rows under the 120s cap, with the expected `ss1`/`mac_econ_fwd500`
timeouts and the known singular `bips07_1998`; common rows versus the previous
broad KLS artifact improved geomean by about 2%, though one-pass noise still
dominates many sub-millisecond rows.

The single-block EGraph schedule/consumer gate was then lowered from the
previous very-large-only floor to cover moderate single-block cases whose
actual factor work, LU fill, and measured dependency work are already high
enough to amortize schedule construction and worker scratch setup. This is a
paper-derived generalization of the CKTSO/NICSLU intra-block dependency-graph
idea, not a matrix-name rule: single-block scheduling now starts at about
`3e8` factor flops or 3M factor nonzeros, and the no-pivot EGraph consumer
requires about `1.5e8` dependency-work units. In same-session SPRAL-enabled
focused checks with four threads, `HTC_336_4438` changed from schedule-only
metadata to actual EGraph consumption and its repeated refactor average dropped
from about 0.105s to 0.041s; `rajat24` similarly dropped to about 0.096s
refactor average. A nine-row hard-focus JSONL improved geomean SPICE-cycle time
by about 1.17x versus the previous default-on artifact, with no loss over 2%.
The same comparison against CKTSO still leaves a large gap, about 2.91x
geomean on those focused rows, so this is retained as a useful KLS-owned
refactor threshold improvement rather than mistaken for CKTSO's full pivoting
tail scheduler.

A later moderate single-block refinement lowered the single-block EGraph floor
for unscaled one-block matrices between 30k and 100k rows when the factored
numeric object already has at least `5e7` measured flops, at least 1M LU
entries, enough level width, and at least `2e7` measured dependency-work units.
This is intentionally below the large single-block floor but still excludes
low-work single-block rows such as `bcircuit`, `ACTIVSg10K`, `ACTIVSg70K`, and
`OPF_10000`. In the current medium artifact the selector matches only
`rajat15`; it builds 633 levels and reduces repeated refactor from about
`0.0173s` to about `0.0101s`, moving the focused SPICE-cycle median from about
`2.16s` to about `1.47s`. A full one-pass medium run kept the same three
known failures and moved KLS geomean from the previous `0.3033s` artifact to
about `0.3023s`; the CKTSO ratio improved to about `1.145x` slower on the 90
common completed rows, and the KLU2 comparison improved to about a `2.02x`
geomean speedup on the 88 common completed rows. Large single-block guards
such as `G2_circuit` and `mc2depi` stay on the existing large EGraph path.

The medium-heavy dominant-BTF gate was then rechecked after the newer EGraph
cluster/pipeline scheduler and the single-block threshold work. The current
retained branch lowers the 85-95% coverage class to about `1.5e8` actual
factor flops and lets that class consume the EGraph path when measured
dependency work reaches about `8e7`. This brings `onetone2` into the same
structural policy as `onetone1`, without reopening the rejected 95%+ Rajat
class: `rajat28` still reports zero EGraph schedule metrics. In same-session
SPRAL-enabled focused checks with four threads, `onetone2` built 1010 EGraph
levels and repeated refactor dropped from about 0.040s to about 0.015s. The
nine-row hard-focus JSONL improved geomean SPICE-cycle time by about 1.09x
over the previous single-block-threshold artifact, and the KLS/CKTSO focused
geomean gap moved from about 2.91x to about 2.68x. CKTSO is still materially
faster, so the remaining gap still points to the larger pivoting-tail,
numeric-kernel, and solve-scheduler items rather than more ordering backends.

The single-block EGraph column fast path was then revisited with a narrower
scope than the earlier rejected scaled experiment. The retained version only
applies when the factor has one BTF block and no active KLU row scaling, so the
previous `nxp1`/`rajat30` scaled guards remain on the generic column kernel. In
that unscaled class, the hot loop bypasses repeated BTF-block lookup, BTF
off-block checks, and scaling branches. Same-session focused checks showed a
small net gain: the nine-row hard-focus JSONL improved by about 0.7% geomean
over the medium-dominant-block artifact, with `G2_circuit` moving from about
48.8s to about 47.8s and `HTC_336_4438` from about 8.64s to about 8.25s on the
SPICE-cycle estimate. A one-refactor `mc2depi` guard stayed in the same
refactor class, about 2.23s. This is retained as a minor EGraph kernel cleanup,
not as the missing CKTSO-scale scheduler.

Two follow-up EGraph hot-loop/scheduler probes were rejected after the retained
single-block fast path. First, the unscaled single-block fast path was changed
to trust the previously validated refactor map and U-pattern schedule, removing
per-entry bounds and dependency-order checks from the hot loop. This did not
help `G2_circuit`: refactor stayed in the same noisy 0.45s class, so the
validation branches are not the visible CKTSO gap. Second, the pipeline tail was
changed from the retained round-robin static assignment to contiguous
work-balanced ranges over the tail columns. This was a clear regression on
`G2_circuit`, moving refactor to about 0.94s because later contiguous ranges
wait behind earlier dependency ranges. The retained round-robin tail assignment
therefore remains the right fit for the current EGraph representation. A
separate LU pointer-cache prototype was also removed: caching L/U index and
value pointers with the schedule added memory and regressed the primary
unscaled EGraph guards (`G2_circuit` and `HTC_336_4438`) despite noise-driven
improvement on `rajat28`, which does not consume that cache. These results
make it unlikely that more small EGraph bookkeeping reductions will close the
remaining CKTSO gap.

Three later probes reached the same conclusion. Disabling the forced
all-pipeline path for huge unscaled single-block EGraph factors moved
`G2_circuit` and `mc2depi` backward, with only a small noisy `rajat30`
improvement, so the retained all-pipeline shape is still the better general
dispatch for that class. Retaining worker scratch buffers and stamped
pipeline-completion storage across refactors was numerically correct but
slower on the primary EGraph rows, which means per-refactor allocation is not
the main visible overhead. A focused `nxp1` unscaled-scale trial also did not
produce a safe general policy: no-scale candidates kept valid residuals and
sometimes lowered repeated refactor time, but the result was noisy, estimated
conditioning dropped by about five orders of magnitude versus max scaling, and
post-factor trial cost erased the possible cycle gain.

A follow-up `rajat28` policy sweep also confirmed that the remaining worst
focused-row gap is not a missing scale-mode or static-pivoting toggle. With the
retained static-pivoted AMD/BTF path, auto, no-scale, sum-scale, and max-scale
variants all stayed in the roughly 0.15-0.17s repeated-refactor class, while
CKTSO's saved four-thread artifact is about 0.011s. Disabling static pivoting
was worse: METIS/scale-auto paid about 10.37s initial factor and stayed around
0.21s repeated refactor, and METIS/no-scale paid about 13.93s initial factor
and about 1.54s repeated refactor. This keeps `rajat28` in the dominant-block
numeric-kernel bucket, not the MC64/preprocessing bucket. The MC64-compatible
boundary remains unchanged: use the pinned BSD-licensed SPRAL scaling subset,
a compatible system SPRAL, or independent KLS code; do not vendor HSL MC64 or
solver-tree copies that retain HSL redistribution restrictions.

The BTF worker-pool gate was then narrowed for low-work dominant decompositions
that have one 95%+ diagonal block, thousands of tiny fringe blocks, a largest
block below the current EGraph size floor, and less than `1e9` measured factor
flops. This is not matrix-name tuning: it is the structural case where the
dominant block still runs in one worker and the remaining fringe work is too
small to pay for threaded block scheduling. A broad first version also caught
`ckt11752_dc_1`, whose 172 BTF blocks still benefited from the pool, so the
retained gate requires at least 1024 BTF blocks. On the affected target set
(`ckt11752_dc_1`, `LeGresley_87936`, `rajat20`, `rajat25`, `rajat28`) the
narrow gate kept `ckt11752_dc_1` neutral, had no losses above 2%, and improved
geomean SPICE-cycle time by about 2.5% versus the prior dynamic-pipeline medium
artifact. The hard-focus run stayed essentially neutral against the second
dynamic-pipeline baseline (`9.07s` versus `9.04s` geomean) while improving
`rajat28` from about `17.36s` to about `15.84s` in that comparison. This is a
small policy cleanup. It does not change the main conclusion that the large
CKTSO gap on `rajat20/25/28`, `G2_circuit`, and similar rows requires a
CKTSO/SubtreeLU-style row-oriented numeric kernel, pivoting-tail restart, or
separator-tree/private-pipeline scheduler rather than another ordering or MC64
import.

To keep that diagnosis reproducible, `scripts/decompose_solver_gap.py` now
compares two benchmark JSONL files by phase contribution. On the current
nine-row hard-focus comparison against the saved CKTSO medium artifact, KLS is
still dominated by repeated refactorization: `rajat28` spends about 94% of its
SPICE-cycle estimate in repeated refactor work and that refactor component is
about 13.6x CKTSO's; `G2_circuit` spends about 93% there and is about 6.7x
CKTSO's; `onetone1`, `onetone2`, `rajat24`, and `transient` also have
refactor-component ratios around 2.7x to 3.4x. Solve ratios on those same rows
are only about 1.0x to 1.3x, and `twotone`/`power197k` solve is already faster
than CKTSO. This confirms that a CKTSO-style solve rewrite is secondary for the
current hard gap; the larger missing mechanism is the row-oriented
factor/refactor engine and its pivot-aware scheduler.

A narrow solve-workspace cache was tested and rejected after this diagnosis. The
prototype kept a solver-owned dense permutation buffer for the `row_perm` solve
path so static-pivot and exact-matching cases would not allocate/free an
`n`-entry buffer on every solve. A repeat-20 subset was tempting, with solve
geomean `0.945x` on `power197k`, `rajat20`, `onetone2`, `OPF_10000`,
`LeGresley_87936`, and a non-`row_perm` `G2_circuit` guard, but the broader
row-permutation paper set rejected it: `power197k`, `rajat20`, `rajat25`,
`rajat28`, `onetone1`, `onetone2`, `twotone`, `OPF_10000`,
`LeGresley_87936`, `rajat22`, `rajat23`, `rajat24`, `hvdc1`, and `hvdc2`
measured `1.007x` slower solve geomean against `a0583c5`, with the static
non-exact rows at `1.010x` slower. The current per-call buffer stays until the
larger KLS-owned LU storage makes a structure-adaptive solve worthwhile.

The scaled serial mapped BTF refactor was then enabled for the narrow
many-fringe dominant-BTF shape that the worker-pool narrowing intentionally
left serial: at least 1024 BTF blocks, a 95%+ largest block below the EGraph
floor, and at least `1e8` measured factor flops. The broad scale-aware mapped
prototype was rejected because it regressed unrelated low-work or small-fringe
scaled rows such as `dc1` and `ckt11752_dc_1`. The retained guard is structural
and, on the current medium paper artifact, selects only `rajat20` and
`rajat28`. On the targeted dominant-fringe set it improved geomean SPICE-cycle
time by about 2.4% with no losses over 2%; on the nine-row hard-focus set it was
about 1.2% faster than the previous many-fringe-serial artifact, mostly from
`rajat28`. The full 93-row medium manifest completed the same 90 rows as the
previous KLS artifact, with the known singular `bips07_1998` and the known
`ss1`/`mac_econ_fwd500` timeouts, and improved one-pass geomean from about
`0.3407s` to `0.3372s`. Against the saved CKTSO medium artifact, however, KLS
is still about `1.28x` slower geomean and `rajat28` remains about `9.3x` slower,
so this is another small fixed-pattern refactor cleanup rather than the missing
CKTSO-scale mechanism.

That same many-fringe dominant-BTF class was then revisited using the paper
idea that cluster barriers can dominate when exact dependencies are already
known. Instead of lowering the old barriered EGraph gate again, KLS now runs
this class through an all-pipeline exact-EGraph consumer: workers claim columns
in topological order and wait only for actual U-pattern predecessors. In the
current 93-row medium artifact this structural gate selects `rajat20`,
`rajat25`, and `rajat28`. On the repeat-heavy target set it improved geomean
SPICE-cycle time by about `1.71x` versus the previous guarded-scaled-map
artifact, with `rajat20/25/28` each moving to roughly 40% of their prior cycle
time. The full medium manifest kept the same 90 completed rows and the same
three failures, improving KLS one-pass geomean from about `0.3372s` to
`0.3275s`. Against CKTSO, the medium geomean gap moved from about `1.28x` to
about `1.24x`, and `rajat28` moved from about `9.3x` slower to about `3.8x`
slower. This is meaningful scheduler progress, but the remaining gap on
`G2_circuit`, ASIC rows, and the Rajat rows still points to the larger
row-oriented numeric kernel, pivot-aware scheduler, and separator/supernode
work described in the papers.

A small reactive static-pivoting payback gate was then retained. Earlier KLS
policy could launch the post-factor static row-matching trial on small,
low-work matrices where the accepted permutation improved pivoting and fill but
did not recover its setup cost over the 99-refactor SPICE-cycle estimate. The
new gate applies only after a first factorization has succeeded, only below
20k rows, and only when diagonal weakness is not severe; severe missing-diagonal
cases still use the pre-static path, and larger or higher-work cases keep the
existing matching policy. On the full medium paper manifest this kept the same
90 completed rows and the same known `bips07_1998`, `ss1`, and
`mac_econ_fwd500` failures, while improving KLS geomean from about `0.3141s`
to about `0.3096s`. Against the saved CKTSO artifact the geomean ratio moved
from about `1.190x` slower to about `1.173x` slower. The largest retained wins
were low-work `OPF_3754`, `bips98_*`, and `nopss_11k` cases; the phase
decomposition after this change still shows the remaining largest losses are
dominated by repeated refactor throughput (`G2_circuit`, ASIC, `mc2depi`, and
Rajat rows), not by another MC64-compatible matching import.

A follow-up low-work ordering/payback refinement tightened the medium
bounded-degree METIS-start class from 90% diagonal-present to 99%
diagonal-present. On the medium paper corpus this only affects the
`rajat03`/Rommes-BIPS/nopss structural class: `rajat03` remains a near-full
diagonal METIS start, while `bips98_606`, `bips98_1142`, `bips98_1450`, and
`nopss_11k` stay on AMD. The reactive static-match payback gate was also
extended to small many-block low-work cases with less than 5% weak and less
than 5% missing diagonal rows, preventing `bips98_1450` from replacing the
faster AMD factorization with a slower static-match candidate. A focused
five-row same-session check improved the geomean SPICE-cycle estimate from
about `0.0874s` under the previous auto policy to about `0.0564s`, while
keeping `rajat03` on METIS. The full 93-row medium paper run kept the same
three known failures (`bips07_1998`, `ss1`, and `mac_econ_fwd500`) and moved
KLS geomean from the previous saved `0.3096s` payback artifact to about
`0.3033s`; the CKTSO comparison ratio improved to about `1.149x` slower on the
90 common completed rows. Against the saved KLU2 artifact, KLS now wins 73 of
88 common completed rows with about a `2.01x` geomean speedup.

The EGraph floor was then lowered for a compact unscaled dominant-BTF shape:
95%+ largest-block coverage, 8-512 BTF blocks, a 10k-30k largest block, and at
least `2e7` measured factor flops. This is the smaller analogue of the retained
low-work dominant-BTF EGraph policy and intentionally excludes the TSOPF rows
with only two BTF blocks and the scaled `ckt11752_dc_1` shape. On the saved
medium artifact this selector matches only `coupled`. A focused three-pass run
with ten repeated refactors moved `coupled` from the prior `0.6345s` saved
SPICE-cycle estimate to a median `0.3783s`; the new row records 540 dependency
levels and about `1.19e7` dependency-work units. A one-pass full-medium run was
noise dominated overall (`0.3108s` geomean versus `0.3096s` in the previous
artifact), but the only row whose dependency schedule changed was `coupled`,
and that row improved from `0.6345s` to `0.3858s` in the full run. The retained
conclusion is narrow: compact dominant-BTF refactors can consume the exact
EGraph path, but this still does not address the much larger single-block and
ASIC refactor-kernel gap.

A scaled medium-dominant BTF EGraph gate was then retained for the adjacent
scaled shape that the unscaled compact gate deliberately skipped. The retained
selector requires BTF, 8-512 blocks, 95%+ largest-block coverage, a 30k-60k
largest block, active KLU scaling, at least `3e7` measured factor flops, at
least 1M LU entries, and at least `1.5e7` measured dependency-work units before
the consumer runs. In the current medium artifact this matches only
`ckt11752_dc_1`; a three-pass focused check built 1360 levels and reduced
repeated refactor from about `0.0100s` to about `0.0071s`, moving the focused
SPICE-cycle median from about `1.16s` to about `0.91s`. A full one-pass medium
run kept the same three known failures and moved geomean from about `0.3023s`
to about `0.3014s`; the CKTSO ratio improved to about `1.142x` slower and the
KLU2 comparison improved to about a `2.02x` geomean speedup. Nearby IBM
`dc*`/`trans*` scaled dominant-BTF guards keep their existing large-dominant
EGraph path.

The symbolic BTF retry policy was then extended to a large fragmented-BTF shape:
at least 100k rows, 8-512 BTF blocks, largest block below half the matrix, and
at least `5e7` estimated BTF flops. This is deliberately separate from the
dominant-BTF and inflated-many-block retries: it targets cases where BTF leaves
many off-block entries and the no-BTF symbolic score is at least 10% lower. In
the full medium artifact this changed only `hvdc2`, moving it from AMD/BTF to
AMD/no-BTF. The row's SPICE-cycle estimate dropped from about `4.88s` to about
`2.21s`, repeated refactor from about `0.0336s` to about `0.0154s`, and
off-diagonal pivots from 1019 to 8 with a valid residual. The full 93-row
medium run kept the same three known failures (`bips07_1998`, `ss1`, and
`mac_econ_fwd500`) and moved KLS geomean from about `0.3014s` to about
`0.2952s`; against saved CKTSO the ratio improved to about `1.118x` slower,
and against saved KLU2 the common-row geomean speedup improved to about
`2.06x`.

The BTF retry was then extended downward to a low-work medium many-block class:
4k-90k rows, at least 1024 BTF blocks, a 5-82% largest block, at most `2e7`
estimated BTF flops, and at least 95.5% structural diagonal coverage. The
no-BTF candidate is accepted only when symbolic fill stays within 1.75x and
estimated flops within 1.5x of the BTF symbolic. Focused checks showed forced
no-BTF wins for Bomhof `circuit_2/3/4`, Rommes BIPS/MIMO/NOPSS variants, and a
small `rajat22` improvement, while the diagonal-completeness gate kept the
noisy `rajat26/27` edge cases on the older BTF/static paths. The full 93-row
medium run again had the same three known failures and moved KLS geomean from
about `0.2952s` to about `0.2877s`. On the 90 rows common with saved CKTSO,
KLS is still about `1.09x` slower; on the 88 rows common with saved KLU2, KLS
is about `2.12x` faster. This is a retained general low-work BTF-overhead
reduction, not evidence that ordering alone closes the remaining CKTSO gap.

The largest-block ceiling was then widened from 80% to 82%, still under the
same low-work and diagonal-completeness guards. In the current medium corpus
this adds only `bips98_606`, another Rommes low-work matrix whose no-BTF
symbolic roughly halves the estimated flops without the tiny-block blow-up seen
on Sandia `mult_dcop_*` and `TSOPF_RS_b9_c6`. The full medium artifact again
kept the same three failures, moved KLS geomean to about `0.2853s`, narrowed
the CKTSO common-row ratio to about `1.081x`, and improved the KLU2 common-row
speedup to about `2.14x`.

The next retained policy pass addressed three low-work preprocessing costs
without changing the MC64 compatibility boundary. First, OPF-style bounded
degree matrices whose numeric diagonal is roughly half missing or weak now
start unscaled, and the score-gated single-block no-BTF retry is allowed down
to 12k rows for AMD/COLAMD auto candidates. This lets `OPF_3754` keep a
no-BTF/unscaled factor while still requiring symbolic evidence. Second, a
medium many-block, mostly diagonal, high-degree spike class whose largest BTF
block is 70-92% of the matrix now starts in KLU scale mode `0`; in the current
medium corpus the structural predicate matches only `rajat16`, `rajat17`,
`rajat18`, and `rajat26`, and it deliberately excludes the accepted-static
95%+ dominant Rajat rows. Third, the partial-weak pre-static matching gate now
requires at least five nonzeros per row on average, avoiding the expensive
rejected row-matching trial on sparse full-diagonal spike cases such as
`circuit_4` while preserving the denser Rajat static-match cases. The full
93-row medium run kept the same three failures, moved KLS geomean from about
`0.2853s` to about `0.2827s`, narrowed the CKTSO common-row ratio to about
`1.071x`, and improved the KLU2 common-row speedup to about `2.16x`. This is
still a preprocessing/payback refinement; the remaining large CKTSO gap is in
the KLS numeric refactor kernel and scheduling path.

The exact-EGraph all-pipeline mode was then extended from the retained
many-fringe dominant-BTF class to very high-work unscaled single-block
refactors. The selector is structural: one BTF block, no active KLU row
scaling, at least 100k rows, at least `1e9` measured factor flops, and at least
`1e9` measured no-pivot dependency-work before the EGraph consumer is allowed
to drop all cluster barriers. This targets the CKTSO/NICSLU observation that
pipeline mode can expose useful dependent-row overlap when the exact
dependency graph is large enough, while keeping lower-work unscaled HTC rows
and scaled `Raj1` on the existing barriered cluster/pipeline split. In focused
checks, `G2_circuit` switched from 501 cluster levels and 1617 pipeline-tail
columns to full-matrix all-pipeline execution, reducing repeated refactor from
about `0.447s` in the saved baseline artifact to about `0.392s`; `mc2depi`
also switched to full-matrix all-pipeline execution and moved from about
`2.16s` to about `2.10s` repeated refactor. A same-session full 93-row medium
comparison against a clean `HEAD` worktree kept the same known failures
(`bips07_1998`, `ss1`, and `mac_econ_fwd500`) and changed schedule metrics only
for `G2_circuit` and `mc2depi`; geomean moved from about `0.2869s` to about
`0.2859s`. The saved CKTSO comparison still shows KLS materially behind on
these rows, so this is a narrow scheduler improvement, not the missing
pivot-aware row-oriented numeric engine.

The same all-pipeline single-block selector was then widened to very high-work
row-scaled factors when KLU row scales and `Pnum` are already present, so KLS
can recompute `Rs`, execute the exact EGraph without cluster barriers, and
permute `Rs` back to pivot order after refactor. This keeps the same one-block,
100k-row, `1e9` factor-flop floor and does not affect lower-work scaled cases.
In the saved medium and large paper artifacts this structural gate matches only
`nxp1`. A same-session three-pass comparison against a clean baseline moved
`nxp1` median repeated refactor from about `0.677s` to about `0.649s`, median
fast-factor refactor from about `0.708s` to about `0.670s`, and the 100-step
SPICE-cycle estimate from about `75.7s` to about `73.4s`, while the solve time
stayed neutral.

Two follow-up scheduler probes were rejected after the scaled huge single-block
all-pipeline change. First, the all-pipeline work cursor was changed to claim
four-column dynamic chunks for modest average pipeline work. This reduced
atomic cursor traffic but delayed dependency publication inside each chunk and
was a major regression: on `rajat20`, `rajat25`, and `rajat28` repeated
refactor more than doubled, `rajat30` regressed by about 34%, and `nxp1` by
about 23%. Second, the fragmented non-dominant many-block EGraph shape used by
the ASIC 680k class was switched from the retained cluster/pipeline split to a
full all-pipeline run. That was neutral-to-worse: `ASIC_680k` was essentially
flat, `ASIC_680ks` regressed by about 4%, and the focused geomean regressed
slightly. The retained per-column dynamic cursor and barriered fragmented-BTF
split remain the better fit for the current fixed-pivot LU representation.

A third scheduler probe tried to skip `pipeline_done` publication for columns
with no successor in the retained U-pattern dependency graph. This was correct
but not useful on the current pipeline representation: a same-session focused
comparison against `23253eb` produced a `1.027x` repeated-refactor geomean
regression across `rajat20`, `rajat25`, `rajat28`, `rajat30`, `nxp1`,
`G2_circuit`, `mc2depi`, and `ASIC_680k`. The small wins on `rajat30`,
`G2_circuit`, and `mc2depi` did not offset `rajat25` at `1.088x`, `rajat28`
at `1.121x`, `nxp1` at `1.020x`, and `ASIC_680k` at `1.016x`. The extra
metadata load and changed publication pattern are therefore not a general
substitute for a larger CKTSO-style pivoting scheduler change.

A scaled refactor-map hot-loop probe was also rejected. The prototype stored the
original input row beside each retained refactor-map entry so scaled mapped
refactors could divide by `Rs[oldrow]` without chasing `row_idx[input_pos]`.
This looked like a cheap CKTSO-style row/segment metadata step, but the extra
map memory and load did not pay for itself on the current KLU storage. A
same-session focused comparison against `989360b` regressed repeated refactor by
`1.026x` geomean across `nxp1`, `rajat20`, `rajat28`, `Raj1`, `dc2`,
`G2_circuit`, `mc2depi`, `rajat25`, `rajat30`, and `ASIC_680k`; the scaled rows
alone regressed by `1.032x`, and the unscaled guards by `1.021x`. The prototype
was removed, leaving the retained refactor map limited to the map data that has
shown a general win.

The CKTSO paper was re-read after these scheduler probes because the remaining
gap is too large to explain by small EGraph bookkeeping. The missing mechanism
is larger and architectural: CKTSO's fast factorization is a row-oriented
up-looking kernel that first assumes the previous pivot order, schedules the
guessed EGraph in cluster/pipeline modes, checks the pivot against the current
row maximum, and, if the check fails, restarts only the ETree-descendant tail
with pivoting. KLS currently approximates only the no-pivot EGraph half on top
of KLU's column-oriented LU storage; the retained restart path invokes KLU's
serial pivoting kernel on the whole rejected BTF block. That is robust and
LGPL-compatible, but it cannot express CKTSO's partially recomputed
ETree-descendant tail or SubtreeLU's separator-tree private/pipeline queues.
The broad CKTSO gap on `G2_circuit`, `mc2depi`, `rajat20`, `rajat25`,
`rajat28`, `rajat30`, and `nxp1` should therefore be treated as evidence for a
new KLS-owned row/segment-oriented numeric layer, not another matrix-specific
ordering, scaling, MC64 import, or EGraph micro-optimization.

KLS now has a first KLS-owned scalar row-major no-pivot refactor scaffold behind
`KLS_ENABLE_ROW_REFACTOR=1` for unscaled single-block factors. It builds row
views of the fixed L pattern, U pattern, and factor-order input pattern, then
recomputes L row entries before writing each U row. Smoke checks on `add20` and
`G2_circuit` kept valid residuals, but the scalar row order is not a default
policy: on same-session one-thread checks, `add20` refactor was about
`0.00013s` versus the default `0.00006s`, and `G2_circuit` was about `1.19s`
versus `0.82s` over five repeated refactors. This narrows the missing paper
work: the row layer needs parallel scheduling and supernode/segment updates, not
just a scalar transposition of the KLU column kernel.

The row-major scaffold then gained an exact L-row dependency schedule and a
gated cluster-mode threaded execution path using the existing KLS worker pool.
This validates the CKTSO/SubtreeLU dependency direction without changing default
policy. Same-session four-thread checks kept valid residuals, but the
barriered scalar row tasks still lagged the retained column EGraph path:
`G2_circuit` measured about `0.86s` refactor with
`KLS_ENABLE_ROW_REFACTOR=1`, versus about `0.205s` for the default exact
U-pattern EGraph refactor. This points the next row-kernel work toward
supernode/segment updates and pipeline scheduling, not plain row-level cluster
barriers.

The row schedule was then upgraded to execute exact adjacent supernode-candidate
segments as single tasks, removing barriers between rows in the same dense-block
candidate. The grouped scheduler remained correct, but was still not a default
policy: on the same four-thread `G2_circuit` check, the grouped row path measured
about `0.88s` refactor, essentially no better than plain row-level barriers and
far behind the current column EGraph path. This confirms that merely coarsening
row tasks is not the missing paper mechanism; KLS needs a real segment kernel
that reuses shared trailing structure and eventually BLAS-style updates.

The grouped row path then gained a first true segment kernel: rows inside an
adjacent supernode-candidate segment use computed dense-block columns and one
shared trailing pattern instead of re-reading each row's full U pattern. This
kept valid residuals and moved `G2_circuit` from about `0.88s` to about `0.84s`
in the gated row path, but it remains far slower than the default column EGraph
path. The retained lesson is that the segment representation is now executable,
but the kernel still needs higher arithmetic intensity, such as dense triangular
mini-solves and batched trailing updates, before it can close the CKTSO gap.
The row-pattern builder now also retains the shared trailing `U` slice offset
for each executable multi-row group and the group dispatcher validates and
consumes that compact descriptor. This is deliberately a small SubtreeLU-style
row-supernode storage step, not a new tuning rule: the clear paper-level missing
piece for the slow cases remains a production row/up-looking factor/refactor
kernel with compact row-major supernode/segment updates and CKTSO's
ETree-descendant pivoting tail, rather than another ordering or scaling switch.

The row-segment path now exposes its executable workload separately from the
earlier adjacent-pattern candidate scan. `kls_stats`, `kls_bench` JSON/text, and
the gap decomposition script report row-refactor group count, group schedule
levels and maximum width, multirow segment count, covered segment rows, maximum
segment width, and dense/trailing segment work. These counters are only
populated when the gated row refactor pattern is built, so default benchmarks
continue to show zero. They let the next dense mini-solve or batched trailing
update experiment distinguish "no segment opportunity" from "segment
opportunity exists but the scalar kernel is still too weak." On `G2_circuit`
with four threads, the gated row path reported 113609 scheduled groups and 8035
multirow segments covering 44528 rows, max segment width 497, 606978 dense
entries, and 5033104 shared trailing entries. Its repeated refactor remained
about `0.906s` versus about `0.239s` for the default column EGraph path in the
same short run, so the counters expose real segment opportunity without
claiming the scalar row kernel is competitive yet.

The retained row-group pattern now also builds the reverse group dependency
graph, including unique group-to-group edges, root/leaf group counts, and
maximum group fanout. This is the row-segment analogue of the EGraph
root/leaf/fanout counters and is a direct bridge toward SubtreeLU-style
private/pipeline partitioning. It is intentionally metadata-only for now:
default KLS still uses the proven column/EGraph or KLU-backed refactor paths
unless the experimental row-refactor environment gates are enabled.

A first dense internal segment mini-solve was then added behind the same gated
row-refactor path. For exact multirow segments whose internal `L` pattern is a
dense lower block and whose dense/trailing work is large enough to amortize the
extra passes, KLS now gathers each row's external-update result into the
segment's own `L`/`U` value slots and applies the internal triangular updates
directly there. Smaller or non-dense-`L` segments keep the existing scalar
segment kernel. This is closer to the SubtreeLU supernode update shape than the
previous per-row scratch update, but it is still not default policy and not a
CKTSO-gap closer by itself: with four threads and five repeated refactors,
`G2_circuit` stayed valid but measured about `0.887s` for the gated row path
versus about `0.224s` for the default column EGraph path. A ten-refactor
`add20` smoke run also stayed valid at about `0.000236s`. The result indicates
that KLS now has an executable in-place dense segment solve scaffold, while the
larger missing piece is still compact row/segment numeric storage and batched
trailing updates with much higher arithmetic intensity.

The dense segment scaffold then stopped rediscovering segment shape on every
refactor. `kls_build_row_refactor_pattern` now precomputes dense-eligible group
flags, per-group trailing lengths, and the per-row internal-`L` offset used by
the gated mini-solve. `kls_stats`, `kls_bench`, and the gap decomposition script
also distinguish all executable row segments from the subset that can consume
the dense descriptor. This is a reusable row/segment descriptor step rather
than a CPU-specific hot-loop tweak. On a four-thread, five-refactor
`G2_circuit` check, KLS reported 8035 executable segments but only 662
dense-eligible segments covering 17070 rows; the gated row path stayed valid
and moved to about `0.846s`, while the default column EGraph path measured
about `0.201s`. On the `add20` ten-refactor smoke run, no segment crossed the
dense-work threshold and the path stayed valid at about `0.000223s`.

The dense segment mini-solve then split internal segment factorization from the
shared trailing-panel update. The dense path now first normalizes and applies
the internal lower/dense-upper block, checks pivots, and then runs the shared
trailing rows as a separate triangular update over the segment. This does not
change the flop count or make the gated row path competitive, but it matches
the paper direction more closely by isolating the panel operation that compact
row/segment storage or BLAS-style packing would later batch. Same-session
checks stayed valid: `G2_circuit` measured about `0.839s` repeated refactor for
the gated row path while the default column EGraph path measured about
`0.232s`, and `add20` stayed valid with no dense-eligible segments at about
`0.000274s`. The result is retained as row/segment scaffolding, not a
production dispatch candidate.

The separated trailing-panel step then gained worker-local compact storage.
Dense-eligible row segments now gather the shared trailing rows into a compact
row-major scratch panel, apply the triangular update there, and scatter the
finished panel back to the current KLU-owned value slots. This is still a
temporary bridge rather than persistent KLS-owned segment storage, but it
removes one layer of pointer chasing from the panel update and matches the
SubtreeLU supernode-storage direction more closely. Same-session checks stayed
valid: a clean `G2_circuit` gated row-refactor sample moved to about `0.798s`
repeated refactor with the same 662 dense-eligible segments covering 17070 rows,
while the default column EGraph path in a same-session check measured about
`0.223s`. The row path therefore remains gated, but this is the first retained
compact-panel step that measurably improves the experimental segment kernel.

The same compact scratch was then extended from the shared trailing panel to the
internal dense segment block. Dense-eligible row segments now gather the raw
internal lower block and dense upper block into one worker-local row-major
panel, perform the internal triangular update there, and scatter the normalized
`L`, dense `U`, and trailing `U` values back to the existing KLU slots only
after the segment is complete. This is still temporary scratch rather than a
persistent KLS numeric format, but it removes the main pointer-chasing loop from
the dense segment mini-solve. A clean four-thread `G2_circuit` sample stayed
valid and moved the gated row-refactor path to about `0.765s`; the default
column EGraph path in the same session measured about `0.205s`. This confirms
that compact segment storage is the right direction for the experimental row
kernel, while also confirming that the row kernel still needs much more work
before it can replace the current default.

The grouped row-refactor scheduler then gained work-balanced per-thread slices
inside each group level. The old grouped path assigned groups by static stride;
that can leave one thread with dense segment groups while others process mostly
singletons. KLS now builds persistent row-group level partitions from a
structural work estimate and reuses them across gated row refactors, mirroring
the work-balanced cluster slicing retained for the column EGraph path. This is
still a barriered level scheduler, not CKTSO's pivot-aware pipeline tail. It
stayed valid on the same focused checks: `G2_circuit` moved the gated row path
to about `0.739s` repeated refactor, while the same-session default column
EGraph path measured about `0.206s`. The small `add20` gated smoke remained
valid but measured about `0.000269s`, so this remains a large-row/segment
scaffolding step rather than a low-work dispatch policy.

The row refactor then gained a KLS-owned contiguous row-major `U` value mirror.
The retained KLU value slots are still updated so the existing triangular solve
path remains unchanged, but all gated row-refactor update loops now read
previous rows' `U` entries from KLS-owned row-major storage instead of following
one double pointer per value. This is a direct step toward the CKTSO/SubtreeLU
numeric format while preserving current semantics. Focused checks stayed valid:
the four-thread `G2_circuit` gated row path moved to about `0.532s` repeated
refactor, while the same-session default column EGraph path measured about
`0.234s`. The small `add20` gated smoke remained valid but moved from about
`0.000269s` to about `0.000291s`, confirming that this mirror is useful for
large row/update-heavy cases rather than as a low-work dispatch policy.

For patterns with dense-eligible row segments, the row-major `U` mirror then
became authoritative during the gated row refactor: KLS writes KLU's existing
`U` value slots only after the refactor succeeds. Sparse/no-dense row patterns
keep immediate KLU mirroring to avoid adding a final scatter pass where it is
not useful. This structurally gated policy preserves current solve semantics
while reducing pointer writes inside dense-segment hot loops. Focused checks
stayed valid: the four-thread `G2_circuit` gated row path measured about
`0.457s` repeated refactor, while longer no-dense checks remained stable
(`bcircuit` about `0.0050s`, `rajat22` about `0.00164s`, and `add20` about
`0.000209s` repeated refactor in same-session long-repeat samples).

The same dense-segment defer policy then gained a KLS-owned row-major `L` value
mirror. Sparse/no-dense patterns still write KLU `L` values immediately, but
dense-segment row refactors now write both `L` and `U` into KLS-owned mirrors
inside the hot loops and scatter both factor arrays back to KLU only after a
successful refactor. This keeps the current solve ABI while moving the
experimental row path another step toward native CKTSO/SubtreeLU-style numeric
storage. Focused checks stayed valid: the four-thread `G2_circuit` gated row
path measured about `0.441s` repeated refactor, and long no-dense samples
remained stable (`bcircuit` about `0.0046s`, `rajat22` about `0.00132s`, and
`add20` about `0.000175s` repeated refactor in same-session checks).

The dense-segment mini-solve then stopped using the worker-local dense panel as
the authoritative storage for dense-eligible row groups. For those groups, KLS
now writes raw internal `L` candidates, dense `U`, and trailing `U` values
directly into the persistent row-major mirrors, normalizes and updates those
mirrors in place, and leaves the older scratch-panel path as fallback for
non-deferred value storage. This removes one scratch-to-mirror copy layer while
preserving the final scatter to KLU for the current triangular solve ABI.
Focused checks stayed valid: four-thread `G2_circuit` measured about
`0.405s` repeated refactor in a two-repeat sample, down from the previous
retained `0.441s` row-mirror result. No-dense sanity checks still selected zero
dense segments and stayed numerically valid (`bcircuit`, `rajat22`, and
`add20`), while `mc2depi` exercised 1883 dense row segments with a valid
residual. This is a retained row/segment storage step, not the missing
CKTSO-style pivoting tail restart.

A selective dense-group scatter policy was then tested and rejected. The idea
was to let non-dense parallel row-refactor groups write KLU value slots
immediately while scattering only dense-group row mirrors after success. It was
structurally clean, but it moved the two intended dense-segment targets in the
wrong direction in same-session checks: `G2_circuit` regressed to about
`0.458s` repeated refactor from the retained native-mirror sample around
`0.405s`, and `mc2depi` regressed to about `2.39s` from about `2.19s`.
Keeping a uniform deferred scatter when dense row groups are present therefore
remains the better policy for the current hybrid KLS/KLU storage.

The native dense-group path was also tested with its repeated L/U shape checks
removed from the hot loop, trusting the precomputed dense-group descriptor in
the same way the older scratch-panel path does. This was also rejected:
same-session checks stayed numerically valid but moved `G2_circuit` to about
`0.445s` repeated refactor and `mc2depi` to about `2.28s`, both slower than
the retained native-mirror samples. The checks are therefore left in place until
KLS owns a fuller row/segment descriptor and numeric state that can make those
invariants cheaper to consume.

The experimental row/segment scheduler then gained a CKTSO-style dynamic
pipeline tail after barriered row-group cluster levels.
`kls_build_row_refactor_pattern` already stores levelized executable row
groups; the parallel row path now uses a `2 * threads` width split, processes
wide group levels with barriers, and lets workers claim the remaining
topological group sequence from an atomic cursor while waiting on predecessor
rows through the retained completion bitmap. Benchmark stats now report
row-refactor cluster levels, pipeline groups, pipeline rows, and pipeline work.
A constructed smoke matrix forces this path. A later instrumentation pass fixed
an important measurement ambiguity in this paragraph: enabling the checked row
fast-factor path does not by itself make the following unchecked `kls_refactor`
call use the row kernel. The row-pattern counters can therefore describe a
retained checked factor pass while the measured repeated refactor is still the
default column EGraph kernel. Benchmark JSON now reports
`row_refactor_last_run`, `row_refactor_last_checked`,
`row_refactor_last_parallel`, and row-refactor run counters, and `kls_bench`
accepts `--row-refactor env|off|refactor|checked|all` so future paper-suite
artifacts can separate checked-row factor probes from actual row-refactor
probes. The default KLS path is unchanged; this remains experimental
row-engine scaffolding until the row numeric factor/refactor path is validated
and faster across the broader paper matrix set.

`kls_bench` now has a deterministic diagonal-stress mode for measuring that
restart gap on real paper sparsity patterns without editing MatrixMarket files.
`--stress-diagonal-scale` and `--stress-diagonal-column` keep the first
factorization on the original values, then run repeated factor/refactor/solve
passes on values with selected diagonal entries scaled. With
`--repeat 1 --refactor-repeat 0`, the JSON `fast_rejected_*` fields can expose
the whole rejected BTF block, the suffix from the failed pivot, the exact
U-pattern descendant tail, the ordered-block ETree successor path, and the
refresh-state-seeded pivoting-tail scope that a CKTSO-style pivoting tail
restart would target. This is a diagnostic for architectural work, not a tuning
path for specific matrices.
The same diagnostic now also reports `fast_rejected_refresh_state`: unknown,
prefix-current, or all-current. This makes the benchmark artifact distinguish
KLS paths that can safely do block/tail continuation from parallel or KLU
all-refresh paths that must remain conservative until KLS owns the row/segment
numeric state needed for CKTSO-style ETree-descendant restart.

A follow-up removed redundant passes from that diagnostic path for unscaled
repairs whose surrounding numeric state is already current. After a failed
checked fast factorization, KLS tracks whether the failed pass refreshed all
columns, only a serial prefix, or stopped in an unknown parallel/EGraph state.
If all columns are current, or if the serial prefix reached the final BTF
block, KLS validates the repaired block and later block tail and returns
without immediately refactoring the same block again with the new fixed order;
unchanged earlier blocks are no longer rescanned on that shortcut path.
Same-session diagonal-stress checks kept valid residuals and reduced the
stressed `add20` factor path from the earlier `3.50s` to about `1.34s`, and
stressed `rajat03` from about `0.44s` to about `0.21s`. Other multi-block
repairs still retry the checked fast factorization because blocks after the
rejected block may not have been refreshed when the fast path stopped.

A block-level tail continuation was then retained for BTF cases. If the failed
checked pass stopped after a true prefix and the rejected block was not final,
KLS repairs the rejected block with pivoting and then refactors only later BTF
blocks with a checked continuation. The original unscaled serial continuation
was later opened to the threaded BTF worker pool, and the same input-row scale
invariant now lets KLU row-scaled repairs try that pool before falling back to
the serial checked scaled continuation. Same-session diagonal stress checks
against a clean `09d074a` baseline showed the intended non-final BTF benefit
with valid residuals: `coupled` moved from about `0.53-0.74s` to
`0.26-0.38s` factor time, and `circuit_1` moved from about
`0.012-0.014s` to `0.0096-0.0104s`. This is a useful CKTSO-aligned block-tail
step, but it still does not implement the larger missing single-block
ETree-descendant pivoting tail inside a KLS-owned row/segment numeric engine.

The experimental KLS-owned row refactor now runs the serial checked fast
factorization pass by default unless `KLS_ENABLE_CHECKED_ROW_REFACTOR=0`
explicitly disables it. It checks the
same no-pivot `L` multipliers used by the existing column fast factorization,
records the dependency pivot that fails, marks the refreshed state as a serial
prefix, and reuses the current block-restart machinery. This is intentionally
separate from `KLS_ENABLE_ROW_REFACTOR`, so the existing unchecked repeated
row-refactor path is unchanged and the checked CKTSO-style guessed-row pass is
the first checked fast-factor scheduler attempt rather than an opt-in
diagnostic. Focused checks stayed valid: normal `add20` used the serial
checked row path with no reject, while stressed `add20`
(`--stress-diagonal-scale 1e-9`) reported a prefix-current fast reject at
the intended dependency pivot.

That checked row pass can now use the experimental parallel row/segment
scheduler when multiple threads are available. The threaded row tasks share the
same earliest-rejected-pivot recorder as the EGraph refactor path; dense and
non-dense row-segment kernels all check the same multiplier predicate. A
parallel checked-row rejection now reports
`KLS_FAST_REJECT_REFRESH_PREFIX` only when the retained completion bitmap proves
all rows before the rejected pivot have finished. If dense row segments deferred
their writes in KLS-owned row-major mirrors, KLS scatters only that proven
prefix into the KLU numeric object before reporting prefix-current. Otherwise
it remains
`KLS_FAST_REJECT_REFRESH_UNKNOWN`, because rows in the active level may finish
out of factor-order prefix before the stop flag is observed. Smoke coverage now
includes both a known weak-pivot 3-row prefix plus independent diagonal work to
force a two-thread row schedule, and a generated dense-row-segment case whose
prefix values change before a later rejected multiplier. The dense case checks
that the pivoting-tail restart consumes the published prefix and leaves a
bounded residual. This is a real CKTSO-style fast-factor-with-pivot-check
execution slice, but it still falls back to full block repair when the prefix
proof is unavailable rather than consuming the full ETree tail worklist.

The row-refactor pattern can now also retain the ordered-block ETree parent
array lazily for single-block row-major checked rejects. Rejected-pivot
diagnostics consume that cached parent after the first checked row-refactor
failure, instead of rebuilding the same CKTSO tail scope on later failures, while
the normal unchecked row-refactor path does not pay the setup cost. The smoke
test forces `KLS_ENABLE_CHECKED_ROW_REFACTOR` through a deterministic
pivot-reject repair so the retained row metadata, block restart, and ETree-tail
stats stay covered. This is a scheduler-metadata step toward the missing
pivoting tail restart, not a replacement for the row/segment-owned pivoting
factor kernel.

The checked row-refactor metadata now also retains reverse row dependencies
derived from the row-major `L` mirror. When a checked row refactor rejects a
pivot, KLS marks the exact row-successor tail, compacts it into increasing
row order, and reports the `fast_rejected_row_tail_*` diagnostics from that
retained topological list. KLS-owned pivot checks also record the factor row,
multiplier magnitude, accepted pivot magnitude, and candidate entry magnitude
that tripped the reject. For checked row-major rejects, KLS also scans the
retained row-tail list using the current prefix state and records the strongest
tail candidate row/value plus the number of tail rows that challenge the
current pivot tolerance. The diagnostic also records the candidate's retained
row-tail position and a `fast_rejected_tail_repair_ready` flag when the
checked row-major prefix is current and the best row-tail candidate satisfies
the same pivot-tolerance predicate that rejected the reused pivot. This gives a
row/segment pivoting tail kernel the local row candidate, tail location, and
value comparison that the older pivot-only diagnostics lacked.
This keeps the CKTSO restart target tied to KLS row storage rather than only
the KLU-column U-pattern or an ETree upper bound; the remaining missing piece
is the full CKTSO-style pipelined row/segment tail kernel.

The retained row-tail list remains an executable fallback seed for the
conservative pivoting-tail envelope. For prefix-current checked row-major
rejects with no saved interrupted worker seed, KLS can mark the retained
block-local row tail, close it through the ordered-block ETree, and record the
seed size as `fast_rejected_pivoting_tail_row_seed_columns`; otherwise the
CKTSO-style unfinished guessed-EGraph seed is preferred before row-tail and
suffix fallback plans. This is still not CKTSO's full concurrent pivoting-tail
executor, but it moves the checked row/segment metadata from diagnostics into
the restart planner.

The fallback pivoting block repair now also records the repaired row selected
at the rejected pivot, whether it matches the retained row-tail candidate, the
first pivot whose row changed in the repaired block, and how many changed
pivots occur before versus at/after the rejected pivot. These repair-delta
fields separate cases where the robust full-block KLU repair already preserves
the checked prefix and chooses the row-tail candidate from cases that would
need broader repivoting than a CKTSO-style local tail restart can safely
provide. KLS summarizes the compatible serial subset as
`fast_repaired_tail_restart_ready`, which now requires a prefix-current or
all-current non-root reject, zero changed pivots before the rejected pivot, a
changed suffix pivot, a topological pivoting-tail scope containing the reject,
and a validated live prefix replay. The checked row-tail candidate and
candidate match fields remain diagnostics; they are no longer required before
attempting the serial tail restart. When that flag is true, KLS also records
the fallback full-block repair work, the actual serial suffix column/work
count, and the saved-work estimate the local restart targets.
`scripts/summarize_tail_restart_opportunities.py` now turns those JSONL fields
into suite-level counts, blocker reasons, and largest saved-work opportunities
so the tail-kernel work can be prioritized from broad benchmark evidence. It
also totals and prints row-tail scope work when `fast_rejected_row_tail_*`
fields are present, and compares executed serial suffix-tail work against the
retained CKTSO-style pivoting-tail work to expose where the current fallback
overcomputes the planned restart set.
`scripts/run_bench_suite.py` now forwards the deterministic
`--stress-diagonal-scale` and `--stress-diagonal-column` controls to
`kls_bench`, so these tail-restart opportunity scans can be generated across
paper manifests instead of only from one-off piped benchmark rows.
The first focused stress probe shows why that distinction matters: stressed
`add20` chose the retained row-tail candidate at the rejected pivot, but the
full KLU block repair also changed one pivot before the rejected pivot, so a
future local tail kernel cannot treat tail-candidate availability alone as a
safe replacement for the fallback repair; its strict tail-restart-ready flag is
therefore false even though the full repaired-block work is about `1.66e6`.
The ordinary non-reject `G2_circuit` row-refactor path leaves these fields at
their sentinels.

EGraph checked-refactor rejects now distinguish a provably current prefix from
an unknown partial refresh when the retained pipeline-done generation is
available. KLS marks the reject as prefix-current only if every factor-order
column before the rejected pivot completed in the same EGraph generation; the
existing block-repair path can then continue with the serial BTF block tail
instead of discarding the whole fast pass. Focused diagonal-stress checks on
`coupled` show that the first observed EGraph reject is schedule-dependent:
when the dominant-block-start reject is observed after all earlier columns are
done, KLS records a prefix-current block repair and avoids the later full
checked-pass retry; otherwise it remains unknown and keeps the conservative
fallback. Rows such as stressed `onetone2`, where the completed-prefix proof
fails, remain classified as unknown.

The row-tail diagnostic was then extended beyond the single-block
row-refactor mirror. When that mirror is unavailable, KLS now scans the
rejected BTF block's numeric `L` columns, follows exact local row-successor
edges, and reports the affected tail work with the same per-column work
estimator used by the U-descendant and ETree-tail counters. This does not
change the numeric repair yet, but it exposes the row-tail scope on real
multi-block EGraph rejects. Repeated focused stress checks kept valid
residuals and reported nonzero row-tail scopes for `coupled` and `onetone2`,
where the previous JSON fields were zero.
The same diagnostic now preserves a checked-pivot row/value as a tail-candidate
row when that row lies in the recorded row-tail scope. The candidate is marked
repair-ready only for prefix-current rejects; unknown partial-refresh rejects
keep the conservative not-ready classification. Per-reject candidate and repair
fields are reset before each new reject so multi-restart stress rows do not mix
diagnostics from different attempts. The KLU-storage checked refactor kernels
now record the strongest violating entry in the rejected column instead of the
first entry visited by the sparse `L` pattern. This does not change the
accepted success-path factors or the fallback decision, but it makes the
CKTSO-tail diagnostic sharper: if the fallback pivot row still differs from
this strongest current-column candidate, then the missing paper mechanism is
broader pivoting-tail state, not merely an early-exit artifact of KLS's
diagnostic loop. The tail-opportunity summarizer now separates missing or
non-tolerance-valid candidates from unknown refresh state, all-current reject
state, and other not-ready cases so the paper-derived restart blocker is not
misreported as a pivot-search failure. A focused stressed probe on `coupled`,
`onetone2`, and `hvdc1` kept valid residuals; the two large repaired blockers
reported `tail_row == repair_row` and were instead classified as
`unknown_refresh_state`, confirming that the next missing CKTSO mechanism is
prefix-safe pivoting-tail execution rather than merely finding a local row
candidate.

The EGraph refactor worker scratch allocation was then narrowed for BTF paths:
single-block refactors still allocate one dense `n`-entry vector per worker,
but BTF EGraph workers only need block-local indices and now allocate
`maxblock` entries. This targets the fragmented-BTF/ASIC class where `n` can be
far larger than the largest diagonal block, reducing transient allocation and
zeroing without changing the numeric kernel or matrix-specific dispatch.
Focused post-change checks kept valid residuals; `ASIC_680k` refactor stayed
around `0.07-0.075s` and `rajat28` around `0.051-0.056s`, so this should be
treated as a memory-footprint cleanup rather than the missing CKTSO-scale timing
fix.

A small fragmented-BTF scheduler fix was then retained. The non-dominant
many-block EGraph path had hundreds of thousands of singleton BTF blocks at
level 0, but those singleton columns had zero scheduling weight because their
numeric work is outside the multi-column block LU estimates. The clustered
level partitioner could therefore hand a huge count of trivial-but-not-free
singleton columns to one worker before the substantial block work started. KLS
now assigns a unit balancing weight to singleton BTF columns only for the
fragmented non-dominant many-block shape. A same-session comparison against a
clean `d9f4c29` baseline moved `ASIC_680k` repeated refactor from
`0.0689-0.0695s` to `0.0569-0.0578s` with valid residuals. The gate is not
enabled for dominant-BTF all-pipeline schedules; `rajat28` kept the old
dependency-work and pipeline-work counters and stayed in the same
`0.0547-0.0566s` timing band. This is a general scheduling balance fix for the
fragmented BTF class, not another ordering, scaling, or MC64 import.

The fragmented-BTF EGraph path then stopped binary-searching the BTF boundary
array for every factor-order column. KLS already retains a fixed scatter map
for BTF refactors, so that map now also stores each column's owning BTF block
and the block-aware EGraph column kernel uses it directly. This targets the
same general ASIC-style shape with hundreds of thousands of BTF blocks, without
touching single-block EGraph rows. A five-row ASIC focus comparison against the
prior non-dominant many-block EGraph artifact improved geomean cycle time by
about 1.09x (`9.05s -> 8.31s`), with the largest win on `ASIC_680ks`
(`9.42s -> 6.79s`) and no >2% losses in that focus. Spot checks on `nxp1`,
`rajat28`, `coupled`, and `onetone2` stayed in their expected timing bands.
This is still a scheduler hot-loop cleanup, not the missing CKTSO row-oriented
pivoting factorization engine.

The initial-factor wrapper path then stopped recomputing KLU diagnostics
unconditionally after every optional policy probe. Auto-scale, METIS promotion,
and pivot-tolerance probes now report whether they actually replaced the
numeric object, and KLS recomputes `flops`/`rcond` only when the current numeric
state needs fresh diagnostics for the next decision or final stats. This removes
redundant O(U-pattern) flop scans and O(n) reciprocal-condition scans from the
first-factor path without changing ordering, scaling, or numeric acceptance.
Same-session focused checks were correspondingly modest and noisy:
`ASIC_680k` initial factor moved from about `1.72s` to `1.70s`, `rajat28` from
about `0.676s` to `0.668s`, while `rajat30` and `coupled` were neutral within
noise. This is retained as wrapper overhead cleanup; it does not replace the
larger row/segment numeric engine still needed for the CKTSO-scale gap.

The row-permuted solve wrapper then stopped allocating its dense permutation
scratch on every solve call. Static-pivoting cases reuse one solver-owned
workspace across repeated forward and transpose solves, and the smoke test now
calls the scaled pre-static forward and transpose solves twice to cover the
cached path. Same-session checks against prior artifacts showed forward-solve
improvements on the row-permuted focused rows (`onetone2` about
`0.00116s -> 0.00107s`, `hvdc1` about `0.000315s -> 0.000265s`,
`LeGresley_87936` about `0.00208s -> 0.00157s`, and `rajat23` about
`0.00172s -> 0.00157s`). TSOPF-style rows remained dominated by triangular
work and were neutral within noise. This is a retained repeated-solve overhead
cleanup, not CKTSO's missing structure-adaptive triangular solve.

The solve wrapper then stopped calling the full numeric-stats refresh after
every triangular solve. Solving does not change fill, FLOP, dependency, pivot,
or condition-estimate metadata, so KLS now refreshes only the last kernel
status and memory counters on the solve path. Focused checks kept residuals
valid and showed the expected small wrapper gain: `onetone2` forward solve
averaged about `0.00086s` in a ten-repeat run versus the previous
`0.00101-0.00107s` band, `rajat30` stayed in the retained KLU solve band at
about `0.0336s`, `nxp1` stayed around `0.0445s`, and `ASIC_680k` stayed around
`0.0098s`. This is a small solve-path cleanup, not the missing CKTSO
structure-adaptive solve.

Scaled fast-factor pivot rejection then gained the same block-local repair
path already used for unscaled checked refactors when the rejection came from a
KLS-owned checked path where `Rs` is still in input row order. KLS rebuilds
off-diagonal entries with the same row scaling convention used by KLU, runs the
KLU block factor kernel with the unpermuted row scale vector, then re-permutes
the scale vector after the new pivot order is installed. Later follow-up work
also lets scaled prefix-current repairs resume checked refactorization over
later BTF blocks through the threaded pool before serial fallback, and lets
validated non-root scaled prefix-current repairs use the serial suffix tail
restart inside the rejected block. Scaled KLU-refactor all-refresh rejects now
recompute row scales back to input-row order before trying the same repair path.
Default focused checks (`onetone2`, `rajat30`) stayed in the previous timing
band. This is a CKTSO-aligned coverage improvement for pivot-check restarts,
not a solution to the large `pre2`/`nxp1` gap.

KLS then exposed the fast-reject refresh-state classification in `kls_stats`
and benchmark JSON, and corrected the serial single-block mapped refactor to
record its pivot-check reject as prefix-current instead of unknown. Smoke
coverage now checks the prefix-current state for single-block, scaled
checked-refactor, and BTF block-restart repairs. This does not make KLS faster
by itself, but it turns the CKTSO tail-restart boundary into explicit benchmark
evidence instead of hidden control-flow knowledge.

A direct single-block row-major solve prototype was then tested and rejected.
The prototype built reusable row views of KLU's L and U factors with offsets
back into the existing LU storage, then used row dot-products for single-RHS
non-transpose solves. This matched the CKTSO paper's row-major triangular-solve
direction at a small scope, but it regressed the large single-block losses:
`rajat30` solve time rose from the retained KLU scatter-solve band of about
`0.034-0.036s` to about `0.099s` with three repeats and still about `0.071s`
with ten repeats, while `nxp1` rose from about `0.045s` to about `0.125s`.
The extra row metadata also added large memory traffic for 16-23M factor
entries. The result says the next triangular-solve attempt should not simply
transpose KLU columns into row-offset lists; it needs CKTSO's full
structure-adaptive partitioning with sparse-block/rectangular-slice decisions
and a storage layout designed for that access pattern.

The EGraph refactor path then stopped treating every repeated refactor as a
fresh thread launch. KLS now keeps a solver-owned EGraph worker pool plus dense
worker scratch vectors and dispatches each accepted EGraph refactor as a new
generation through that pool. Scratch is cleared only after early-abort paths
that can leave dense entries live, while successful no-pivot columns retain the
existing invariant that touched slots are zeroed as they go. This is still not
the missing CKTSO row-oriented pivoting-tail engine, but it removes repeated
thread launch and scratch allocation from the current exact-EGraph consumer.
A four-row focused probe with three passes and three refactors per sample
(`onetone2`, `rajat28`, `transient`, `power197k`) improved geomean cycle time
from the scratch-only artifact's `3.5166s` to `3.3196s`; compared with the
saved pre-pool blockmap artifact, `rajat28` moved from `6.3837s` to `5.3433s`
and `onetone2` from `1.8011s` to `1.7070s`, while the mostly initial-factor
`power197k` row stayed neutral. This is retained as a general repeated-refactor
overhead cleanup for EGraph-active SPICE rows.

The same retained EGraph state was then extended to the pipeline dependency
markers. Instead of allocating and zero-initializing an `n`-entry
`pipeline_done` array for every pipeline-capable refactor, KLS now keeps a
solver-owned atomic generation array and marks completed columns with the
current generation. Waiting workers compare against that generation, and the
array is cleared only if the unsigned generation counter wraps. This preserves
the existing dependency protocol while removing a repeated O(n) setup pass from
all-pipeline and pipeline-tail EGraph refactors. On the same four-row focused
probe, geomean cycle time improved from the worker-pool artifact's `3.3196s`
to `3.2480s`; `onetone2` moved from `1.7070s` to `1.5491s`, `transient`
nudged from `3.5853s` to `3.5510s`, and `power197k` stayed neutral because it
does not use the EGraph pipeline path in this run.

All-pipeline EGraph schedule setup was then narrowed by skipping the
per-level thread-slice table. Full all-pipeline runs set the cluster split to
zero and never enter the barriered cluster loop that consumes those slices, so
building the table only adds setup work and memory traffic. This is a small
setup cleanup, not a new numeric kernel: `G2_circuit` moved from the saved
generation-marker sample's `43.2928s` to `43.0747s`, while a three-pass
`rajat28` median was neutral within noise at `5.5145s` versus the prior
`5.4418s` artifact.

The huge single-block all-pipeline schedule then stopped storing and reading
the level-ordered column list. In that shape, factor order is already
topological for the no-pivot U-pattern dependencies, and the existing
generation-stamped wait markers still enforce readiness, so workers can claim
columns directly by factor-order index. The optimization is not used for
dominant-BTF all-pipeline schedules: a first broader trial made `rajat28`
noisier in the suite metric, and direct samples showed the BTF level-order path
should remain the structural default. With the narrowed rule, five direct
`rajat28` samples stayed in the same kernel band (`0.0446s` median refactor
versus the saved `0.0457s` no-thread-slices sample), while three `G2_circuit`
direct samples moved repeated refactor to about `0.218-0.241s` versus the saved
`0.398s` no-thread-slices sample.

The huge single-block unscaled EGraph column kernel then stopped repeating
structural checks that are already proven by the retained refactor map, exact
U-pattern schedule, and EGraph eligibility gate before worker launch. Numeric
singularity and pivot-threshold checks remain in the kernel. Three direct
`G2_circuit` samples moved median repeated refactor from the saved
natural-order sample's `0.2383s` to `0.2188s`, while direct `rajat28` guard
samples stayed in the same scaled dominant-BTF band at about `0.048s`. A
broader attempt to split the BTF value/scaling helper was rejected because it
pushed the scaled `rajat28` guard outside its timing band; that path stays on
the conservative helper logic.

A fresh four-thread medium paper-suite run was recorded after those EGraph
cleanups so future work uses current evidence instead of stale artifacts. With
`--repeat 1 --refactor-repeat 3 --timeout 120`, current KLS completed 90 of 93
medium rows, with `bips07_1998` singular and `ss1`/`mac_econ_fwd500` timing out
as before. Against the saved CKTSO medium artifact on the 90 common completed
rows, KLS is now `1.022x` slower geomean (`0.2698s` versus `0.2639s`), with 40
wins and 49 losses over 2%. Against the saved KLU2 artifact on 88 common rows,
KLS is `2.26x` faster geomean (`0.2385s` versus `0.5399s`) and wins 82 rows
over 2%. The remaining CKTSO losses are still dominated by repeated refactor
throughput: `rajat28`, `rajat25`, `ASIC_320k/320ks`, `rajat20`, `G2_circuit`,
`rajat03`, `ASIC_100ks`, `onetone2`, and `transient` lead the current gap.
Those rows, plus the timeout rows, are now tracked in
`bench/suitesparse_cktso_gap_manifest.txt` for focused regression checks before
rerunning the full medium suite.

The CKTSO gap review then found a concrete missing combination rather than a
generic ordering problem: the medium Rajat dominant-BTF rows benefited from
static row matching plus METIS nested-dissection, but the default fast-factor
path was accepting the static match with AMD. KLS now runs a narrowly gated
numeric METIS refinement for medium weak-diagonal static-match candidates whose
AMD symbolic leaves a many-block BTF with one dominant block. The METIS
candidate is kept only when actual factor fill/flops improve materially and
the reciprocal-condition estimate does not collapse. This changes only
`rajat20`, `rajat25`, and `rajat28` on the full medium run: they switch from
AMD/static to METIS/static, cut numeric fill by about 12-16%, and move their
SPICE-cycle estimates from `5.23s`, `5.10s`, and `5.33s` to `3.80s`,
`3.52s`, and `4.08s`. A 20-row CKTSO-gap guard improved by about 6% geomean
against the previous KLS artifact. The full medium suite remains essentially
flat within one-pass noise (`1.001x` versus the previous KLS artifact) and is
still `1.023x` slower than the saved CKTSO artifact on the 90 common completed
rows, while staying about `2.26x` faster than KLU2 on common rows. The remaining
large CKTSO gap is therefore still repeated refactor throughput on
ASIC/G2/Onetone-style rows, not a missing LGPL-compatible MC64 import or a
single broad ordering switch.

The CKTSO cluster/pipeline split was then rechecked because the CKTSO paper's
`alpha * thread_count` handoff is one of the few remaining explicit scheduling
knobs. KLS keeps CKTSO's `alpha=2` default, but now allows a structural
`alpha=4` handoff for compact large dominant-BTF schedules: unscaled matrices
with 128-512 BTF blocks, a 95%+ dominant block, max block at least 90k, and at
least `5e8` factor flops. This moves the pipeline tail earlier only for that
shape. A direct same-session five-pass A/B on the four touched medium rows
(`kls_asic_cluster_alpha2_ab_t4_p5_r10_timeout180.jsonl` versus
`kls_asic_cluster_alpha4_ab_t4_p5_r10_timeout180.jsonl`) improved geomean by
about `1.2%`: `ASIC_100k` improved `3.0%`, `ASIC_100ks` improved `1.7%`,
`ASIC_320k` improved `0.3%`, and `ASIC_320ks` regressed `0.3%`. Wider
`alpha=6`/`alpha=8` probes were mixed, and a full all-pipeline probe regressed
the same family, so this is retained only as a small structural scheduler
refinement. It does not explain the multi-x CKTSO refactor gap; that still
points to a row/segment-oriented numeric kernel, pivot-checked tail restart,
and eventually a SubtreeLU-style separator-tree/private-pipeline scheduler.

The EGraph refactor kernel then cached KLU's packed L/U column index/value
pointers at numeric-object scope. The previous worker hot loop recomputed the
same `LUbx + Lip/Uip + aligned-index-length` split for every U dependency and
for the output L column. The cache is rebuilt from current `LUbx` storage before
an EGraph refactor first runs, reused across repeated SPICE refactors, and freed
with numeric state or after a pivoting block restart that reallocates `LUbx`.
This is still KLU-storage based, not the row-major numeric engine CKTSO uses,
but it removes repeated wrapper-level pointer arithmetic from the exact-EGraph
no-pivot refactor path. A 10-row focused three-pass repeated-refactor probe
(`kls_lu_pointer_cache_probe_t4_p3_r10_timeout180.jsonl`) improved geomean by
about `2.5%` versus the saved hot-gap baseline and about `2.2%` versus the
current alpha=4 medium artifact on the same rows. The full 93-row medium run
(`kls_lu_pointer_cache_medium93_t4_r3_timeout120.jsonl`) kept the same three
known failures, improved the 27 EGraph-active rows by about `1.1%` geomean, and
was neutral-to-slightly-positive on all 90 completed rows (`0.998x` versus the
previous KLS artifact). The largest retained downside is `onetone2`, which was
about `8.5%` slower in the one-pass full run and about `2%` slower in the
focused run; an attempted structural gate for only unscaled high-work rows was
rejected because the added fallback branch path regressed the focused set
overall and lost the useful low-work `coupled` win.

The compact dominant-BTF EGraph gate was then extended downward for genuinely
small but still nontrivial refactor rows: unscaled matrices with 8-512 BTF
blocks, a 90%+ dominant block, max block in the 3k-10k range, and at least
`2e6` measured factor flops. The retained EGraph consumer still requires
measured dependency work before it runs, so this does not turn every tiny BTF
case into a threaded refactor. On the current medium artifact this structural
predicate matches only `rajat03` and `ACTIVSg2000`, both of which have zero
off-diagonal pivots and therefore are not MC64/matching failures. A five-pass
focused probe with ten refactors per sample
(`kls_small_compact_egraph_probe_t4_p5_r10.jsonl`) moved `rajat03` repeated
refactor from the saved `0.00127s` class to median `0.000696s`, and
`ACTIVSg2000` from `0.000821s` to `0.000379s`. Nearby guard rows such as
`gemat12` and `TSOPF_FS_b9_c1` did not build an EGraph schedule under this
gate. A one-pass CKTSO-gap focus run remained noisy on unchanged large rows,
so this is recorded as a narrow scheduler coverage improvement rather than a
claim that EGraph micro-gating closes the broad CKTSO gap. The MC64 boundary is
unchanged: existing MC64-style code may be reused only when it is
redistributable with LGPL KLS, such as the pinned BSD SPRAL scaling subset or a
verified compatible system library, not HSL MC64 or restricted solver-tree
copies.

An unscaled BTF-specific EGraph column kernel was then tested and rejected. The
prototype bypassed the generic scaling/value helper for unscaled BTF rows and
used the already validated refactor map directly. That looked aligned with
CKTSO's row/segment hot-loop emphasis, but the three-pass focused probe
(`kls_btf_unscaled_egraph_probe_t4_p3_r10_timeout180.jsonl`) regressed the
unscaled BTF EGraph set by about `2.5%` geomean. It helped some low-work IBM
and Rajat rows, but slowed the ASIC family by roughly 3-6% and hurt `coupled`
by about 20%, so the source change was removed. The result reinforces that
duplicating the current KLU-storage column kernel is not the missing
row-oriented CKTSO engine.

The single-block EGraph gate was then extended in the opposite direction:
low-work, unscaled, nearly no-pivot single-block rows. The retained rule is
structural: one BTF block, no KLU row scaling, 15k-250k rows, no more than
eight off-diagonal pivots, `4e6`-`5e7` measured factor flops, and at least
100k numeric LU entries. The EGraph consumer still requires measured
dependency work before it runs. On the current medium artifact this matches
`ACTIVSg10K`, `ACTIVSg70K`, `bcircuit`, `hvdc2`, and `OPF_10000`. A five-pass
focused probe with ten refactors per sample
(`kls_low_work_single_egraph_probe_t4_p5_r10_timeout180.jsonl`) improved every
touched row's repeated-refactor median versus the saved pointer-cache artifact:
`ACTIVSg10K` `0.00136s` to `0.00122s`, `ACTIVSg70K` `0.00338s` to `0.00219s`,
`bcircuit` `0.00414s` to `0.00304s`, `hvdc2` `0.01528s` to `0.01320s`, and
`OPF_10000` `0.00213s` to `0.00197s`. Guard rows that do not satisfy the new
predicate kept zero new EGraph dependency work. This is still a threshold
coverage improvement on the existing no-pivot EGraph refactor, not a
replacement for CKTSO's pivot-aware row-oriented factorization.

A narrow dense-fringe dominant-BTF scale policy was then retained for the IBM
`dc*` shape that still selected KLU max scaling despite having stable pivots
and valid unscaled residuals. The structural predicate is 100k-150k rows, 8-64
BTF blocks, a 99%+ dominant block, and 5-8 nonzeros per row; it only overrides
the scale when the existing value-based pattern policy would otherwise keep
max scaling (`2`). This last condition keeps the same-structure `trans4` and
`trans5` rows on their prior unscaled `-1` path. A five-pass scale-policy probe
with ten refactors per sample moved `dc1` from max scaling to scale `0` and
cut repeated refactor from `0.00810s` to `0.00629s`; `dc2` moved from
`0.00772s` to `0.00645s`; and `dc3` moved from `0.00836s` to `0.00648s`.
Relative residuals stayed around `1e-11` and off-diagonal pivots stayed zero
except for `dc1`, where the faster scale-0 candidate introduced five
off-diagonal pivots but remained numerically valid. Nearby scaled guards
(`ckt11752_dc_1`, `Raj1`, `rajat20`, and `rajat28`) did not change selected
scale. This is a scale-cost cleanup for a repeated-refactor class, not an MC64
import or a substitute for the missing CKTSO row-oriented kernel.

The low-work dominant-BTF EGraph gate was then probed below the older 90k
largest-block floor. A broad first version for 80%+ dominant blocks also
activated on `ckt11752_tr_0` and regressed that guard, so it was narrowed to
the high-coverage no-pivot subcase: unscaled BTF, 30k-120k rows, at most 5000
BTF blocks, a 95%+ dominant block in the 60k-90k range, zero off-diagonal
pivots, and `1e7`-`3e7` measured factor flops. In the current medium evidence
this matches only `LeGresley_87936`. A seven-pass focused probe with ten
refactors per sample
(`kls_low_work_high_coverage_btf_probe_t4_p7_r10_timeout180.jsonl`) moved
`LeGresley_87936` repeated refactor from the current focused `0.00647s` median
to `0.00464s`, while `ckt11752_tr_0` kept zero EGraph dependency work under
the narrowed rule. This is another measured-work scheduler coverage step for
the existing no-pivot EGraph refactor; it does not change the conclusion that
the remaining ASIC/G2/Rajat gap needs a different row-oriented numeric kernel.

The MC64 source boundary was clarified again after accepting that existing
implementations can be reused when they are license-compatible with LGPL KLS.
The rule is license-based rather than authorship-based. The pinned SPRAL
submodule remains the retained vendored source because its BSD-3-Clause license
permits source and binary redistribution with preserved notices, while HSL's
current no-cost licence is personal-use only and does not allow redistribution
in source or binary form. A small BSD `mc64` translation of part of SPRAL is
compatible as a reference, but it does not improve integration over the already
pinned SPRAL C/Fortran interface. The current CMake guard therefore remains
correct: bundled SPRAL checks for redistribution-compatible BSD text, and a
system SPRAL/MC64-style library requires the builder to explicitly set
`KLS_SYSTEM_SPRAL_LGPL_COMPATIBLE=ON` after verifying the selected library.
This keeps MC64-quality preprocessing on the roadmap without admitting HSL
MC64 or restricted solver-tree copies into the LGPL distribution.

The latest discarded prototypes further support the larger-engine diagnosis.
A serial single-block mapped refactor variant reused cached `L`/`U` pointer
indices inside `kls_single_block_mapped_refactor`. It looked attractive as a
row/segment hot-loop cleanup, but same-session A/B checks regressed key
low-work rows (`OPF_3754` about `1.05x`, `xingo_afonso_itaipu` about `1.11x`,
and `ww_vref_6405` about `1.08x` versus the committed baseline), so the patch
was removed. An EGraph wait-loop CPU-relax prototype was also rejected. It
helped some large ASIC samples in one form but regressed `G2_circuit` and
`onetone2`, and narrowing the gate still left mixed results. These are
bookkeeping or spin-wait effects, not the missing CKTSO-scale mechanism.

A current four-thread medium-paper decomposition refresh now shows KLS only
about `1.03x` slower than the saved CKTSO artifact in geomean on the 90 common
completed rows, but the largest losses remain multi-x: `ASIC_320k`,
`ASIC_320ks`, `onetone2`, `ASIC_100ks`, and `G2_circuit` are all above
`2.6x` SPICE-cycle ratio and dominated by repeated refactor work. The refreshed
`bench/suitesparse_cktso_gap_manifest.txt` now tracks the top 40 current
losses plus the timeout rows `mac_econ_fwd500` and `ss1`, so future changes
are checked against the matrices that still expose the large gap instead of
only the older pre-EGraph focus rows.

The BTF EGraph numeric path was then given the same branch-light unscaled
column specialization that single-block EGraph already used. The retained
dispatch is structural: it only applies to unscaled BTF refactors whose largest
block is at least 30k columns, leaving the small compact BTF path on the older
generic kernel after an ungated probe regressed `coupled`. In a same-session
large-BTF A/B with three passes and ten refactors per sample, the retained
candidate improved the nine-row focused geomean by about `1.2%` versus the
committed baseline. Repeated refactor improved on the intended large-block
rows: `ASIC_320k` `0.1033s` to `0.0990s`, `ASIC_320ks` `0.0852s` to
`0.0813s`, `trans5` `0.00665s` to `0.00624s`, and `LeGresley_87936`
`0.00458s` to `0.00420s`. The important unchanged-path guard set
(`G2_circuit`, `onetone1`, `onetone2`, `rajat28`) stayed neutral-to-slightly
positive in a two-pass comparison. This is a real hot-loop cleanup for the
existing EGraph BTF refactor, but it is still an incremental KLU-storage
optimization rather than CKTSO's missing row-oriented pivoting-tail engine.

That BTF-specialized kernel was tightened once more by removing validation
branches already guaranteed by the EGraph dispatcher, refactor map builder, and
LU-pointer-cache builder. The hot path still checks Offx bounds while writing
off-block values, but trusts the already validated block metadata and
topological U pattern, matching the single-block EGraph kernel's trust model.
Against the previous retained specialization, the same nine-row BTF focus
improved by about `1.8%` geomean with no failed rows. The refactor medians
improved on the intended large-block cases: `ASIC_100k` `0.04441s` to
`0.04317s`, `ASIC_100ks` `0.04496s` to `0.04354s`, `ASIC_320k` `0.09904s`
to `0.09701s`, `ASIC_320ks` `0.08130s` to `0.07976s`, and `ASIC_680ks`
`0.04886s` to `0.04740s`. A broader attempt to hoist the
`shared->check_pivots` branch out of every L-update loop was rejected: it
helped some ASIC samples but regressed low-work/scaled guards such as
`LeGresley_87936`, `onetone2`, and `rajat28` under repeated same-session
probes.

Two follow-up BTF bookkeeping probes were also rejected. First, a cached
block-local row map avoided `global_row - block_start` inside the BTF EGraph
kernel, but the extra map memory did not pay back: the nine-row BTF focus
regressed by about `0.5%` geomean versus the retained validation-trim commit,
with clear losses on `LeGresley_87936` and `ckt11752_tr_0`. Second, removing
the remaining per-column LU-pointer-cache validation from the BTF specialized
kernel also regressed the same focus set by about `0.7%` geomean and slowed the
large ASIC rows. The retained BTF specialization therefore keeps the pointer
cache guard and computes block-local row indices directly. The broader signal
is unchanged: small KLU-storage bookkeeping trims are now yielding mixed or
single-digit effects, while the CKTSO gap still requires a row/segment numeric
engine with pivoting-tail restart semantics.

A corresponding single-block EGraph pointer-cache validation trim was also
rejected. Removing the per-column `refactor_lu_pointer_count`/pointer guard
helped `G2_circuit` slightly in a same-session two-pass focus run, but
regressed the Onetone rows by about `10-11%` in SPICE-cycle time and moved the
four-row guard geomean about `4.9%` slower. This confirms the remaining
single-block gap is not a simple pointer-cache guard branch; the retained
single-block kernel should keep its current validation shape until KLS owns a
different row/segment numeric representation.

The MC64 import policy was converted into an executable repository audit after
accepting that existing MC64-style code can be used when it is LGPL-compatible.
`scripts/audit_license_boundary.py` now checks that KLS remains LGPL, that the
only bundled MC64-adjacent implementation used by the build is SPRAL's
BSD-licensed scaling subset, and that first-party/vendor source outside the
vetted SPRAL path does not contain restricted HSL/MC64 markers. This does not
close a performance gap, but it prevents future matching/scaling work from
accidentally crossing the licensing boundary while keeping compatible existing
code available for KLS.

The refactor schedule diagnostics now also record dependency root columns,
leaf columns, and maximum successor fanout. These counters are computed while
KLS already walks the numeric U pattern to build exact EGraph levels. They are
intended to quantify the opportunity for a successor-driven ready scheduler or
row/segment engine on the CKTSO-gap rows, rather than to tune on matrix names
or infer the gap from total edge count alone.

A bounded successor-ready queue was then prototyped for modest EGraph pipeline
tails by retaining reverse U-pattern edges and running eligible pipeline
columns only after their pipeline predecessors completed. The design was
rejected. Same-session A/B against commit `ba8758c` with three passes and ten
refactors per pass showed broad regressions: `onetone1` cycle median
`6.06s -> 18.97s`, `onetone2` `1.60s -> 3.52s`, `ASIC_320k`
`13.62s -> 29.74s`, and `ASIC_320ks` `11.43s -> 25.02s`. Refactor medians
regressed similarly (`ASIC_320k` `0.097s -> 0.261s`, `ASIC_320ks`
`0.079s -> 0.216s`). This means the missing CKTSO mechanism is not a simple
reverse-edge ready queue layered on KLU column storage; KLS needs a different
row/segment numeric representation or separator/private-pipeline engine where
successor scheduling does not add another high-overhead synchronization layer.

The opposite scheduling probe was also rejected: disabling partial EGraph
pipeline tails and forcing the remaining exact-EGraph tail through level
barriers while preserving all-pipeline shapes. Against commit `eac177b`, the
intended ASIC rows regressed sharply in same-session tests with three passes
and ten refactors: `ASIC_320k` cycle median `13.58s -> 30.24s` and refactor
median `0.0976s -> 0.265s`; `ASIC_320ks` cycle median `11.43s -> 25.75s`
and refactor median `0.0797s -> 0.224s`. This confirms that KLS's current
pipeline tail is necessary on these rows, even though it remains much slower
than CKTSO. The remaining gap is therefore not fixed by choosing between
coarse level barriers and KLS's current fetch-and-wait tail; it points back to
the numeric representation and update granularity inside the heavy tail
columns.

All-pipeline single-block scheduling was also probed by replacing natural
column fetch order with exact EGraph level order for huge single-block
refactors. This was meant to expose more independent G2/mc2depi work without
adding a reverse-edge ready queue. It was rejected immediately on the intended
guard: with three passes and five refactors, `G2_circuit` regressed from a
`25.39s` cycle median and `0.220s` refactor median to `43.57s` and `0.404s`.
The natural all-pipeline order is therefore retained for huge single-block
cases. The result further narrows the scheduler diagnosis: changing the order
in which KLS feeds KLU-column kernels is not enough; the missing improvement
needs different per-column update granularity or storage.

The EGraph diagnostics were therefore extended again to report maximum
per-column work and maximum pipeline-tail column work. These fields quantify
when the outer schedule has exposed enough independent columns but one or a
few KLU-column updates still dominate the tail. They are intended to guide the
row/segment numeric-engine work rather than another column-ordering probe.

The immediate hot-loop follow-up was to centralize the repeated sparse
scatter update `x[row] -= L(row,j) * U(j,k)` and unroll it four ways in a
single inline helper used by the serial mapped, BTF, and EGraph refactor
paths. This does not change ordering, matching, pivot policy, or scheduler
semantics; it only reduces overhead in the per-column update loop exposed by
the new maximum-column-work diagnostics. Same-session A/B against commit
`d213eb8` retained the helper: on six heavy CKTSO-gap rows with three passes
and ten refactors per pass, SPICE-cycle geomean improved from `7.77s` to
`7.48s` with no row slower; refactor medians improved on `G2_circuit`,
`onetone1`, `onetone2`, `ASIC_100ks`, `ASIC_320k`, and `ASIC_320ks`. A
one-pass top-20 CKTSO-gap guard also completed all rows and improved geomean
from `3.64s` to `3.55s`; a follow-up three-pass loss check showed the
apparent `rajat28`, `onetone2`, and `ASIC_320ks` losses were noise or reversed
with longer repeats, while tiny `gemat12` remained cycle-noisy despite a
faster measured refactor. This is a small retained KLU-storage kernel cleanup,
not the missing CKTSO row/segment engine.

A wider scatter-unroll follow-up was rejected. Replacing the retained four-way
helper with an eight-way loop, and then with an eight-way-only-for-long-segments
variant, produced mixed same-session results: the six-row heavy focus was only
about `0.6%` positive, the top-20 CKTSO-gap guard was only about `0.3%`
positive, and repeated checks still showed row-level losses such as `rajat28`
and `ASIC_320ks` in refactor time. More importantly, this direction is
CPU-code-generation-specific rather than a solver algorithm improvement. KLS
therefore keeps the simpler four-way scatter helper and leaves further
progress to row/segment numeric storage, pivoting-tail restart, and general
cost-model work.

A branch-light scaled BTF EGraph column specialization was also rejected. The
prototype mirrored the retained unscaled BTF specialization but divided mapped
entries directly by `Rs[oldrow]`, bypassing the generic scaled value helper for
scaled multi-block EGraph refactors with large blocks. It was structurally
sound, built cleanly, and passed the smoke/license tests, but same-session
three-pass checks showed the wrong tradeoff: the scaled-focused eight-row
geomean was only about `0.9%` positive, while the hard scaled Rajat rows
regressed (`rajat20` cycle `3.55s -> 3.59s`, `rajat28` `3.54s -> 3.57s`).
The small `ckt11752_dc_1` win was not enough to justify a duplicate scaled BTF
kernel that worsens the remaining CKTSO-gap rows. KLS therefore keeps the
generic scaled EGraph value path until the larger row/segment numeric engine
can improve scaled and unscaled hard rows together.

A fallback row-segment compact-scratch probe was also rejected. Dense-eligible
segments already keep their compact dense/trailing panel path, but the broader
prototype packed every executable exact-U-pattern multirow segment into
worker-local row-major scratch so later rows could read internal and trailing
U values without KLU double-pointer indirection. It built cleanly and stayed
valid, but the intended four-thread `G2_circuit` guard moved the gated row
path back to about `0.797s` repeated refactor, versus about `0.739s` before
the probe and about `0.233s` for the default column EGraph path in the same
session. The result reinforces the current diagnosis: shallow packing around
KLU storage is not the CKTSO/SubtreeLU mechanism; KLS needs persistent
row/segment numeric storage and scheduler semantics that make segment updates
the native operation rather than a copied side path.

The existing block-restart fallback was then tightened without changing its
numerical semantics. `kls_pivot_restart_rejected_block` previously allocated a
temporary `n`-entry symbolic inverse row map and made multiple full passes over
it before calling the KLU block kernel. It now reuses `numeric->Pinv` as
restart scratch, fills it once from the symbolic row permutation, and rebuilds
the accepted numeric inverse after the repaired block permutation is published.
Focused stressed restart checks on `coupled`, `onetone2`, and `hvdc1` kept
valid residuals and the same restart diagnostics. This reduces current
full-block fallback overhead, but it is still a cleanup around KLU's block
repair path rather than the missing CKTSO pivoting-tail factorization.

The tail-restart diagnostic was then corrected to retain an executable
pivoting-tail plan, not only the old single ETree successor path. For
prefix-current rejects, KLS now seeds the `fast_rejected_pivoting_tail_*` plan
from the checked-refactor unfinished-node bitmap when it is available, and
falls back to the rejected block suffix otherwise; all-current or
unknown-refresh rejects seed the plan from the rejected pivot. The strict
tail-restart saved-work estimate uses this pivoting-tail scope. A focused
`coupled` diagonal-stress probe reported a prefix-current reject at pivot `24`,
a suffix/pivoting-tail scope of `11293` columns and about `1.193e7` work,
versus the older single ETree path of `8194` columns and about `1.191e7` work;
the full repaired-block work remained about `3.650e8`. This keeps the
opportunity estimate tied to the CKTSO restart set and gives the future
row/segment numeric kernel a concrete ordered worklist to consume.

That pivoting-tail plan is now checked as a real scheduler contract instead of
only a count/work estimate. KLS records the first and last global rows in the
retained plan, whether the plan includes the rejected pivot, and whether the
stored order is topologically safe with respect to the ordered-block ETree
parent links. The stressed `coupled` case is still classified as a
root-of-block reject with no reusable serial prefix, but its retained tail plan
is now verified as the ordered worklist a CKTSO-style row/segment tail kernel
would need to consume.

The strict tail-restart readiness gate was then strengthened to compare the
fallback repair against the old block pivot order and to replay KLU's live
prefix bookkeeping through the rejected pivot. The validator rebuilds local
final-to-live row maps, replays KLU `P`/`Pinv` pivot logging, unfinalizes stored
`L` row indices back to local row coordinates, and reconstructs prefix
`Lpend` pruning boundaries from the already-pruned KLU columns. A new smoke
case covers a strict-ready nonzero-prefix repaired block. The same focused
stressed `coupled` probe still reported one strict-ready reject, with about
`3.531e8` saved work over full-block repair. This removes another entry-state
ambiguity before broadening the serial pivoting tail path and before building
the full CKTSO-style pipelined row/segment tail kernel.

The executable serial pivoting-tail restart was then broadened to use the
validated prefix-current proof directly instead of requiring the checked
row-refactor tail-candidate diagnostic. That means default mapped checked
refactors can now run the local tail restart when the rejected pivot is not the
block root, the old prefix is replayable, and the retained pivoting-tail scope
is topological. A new smoke case keeps the default mapped-refactor distinction
explicit: with the checked row-refactor environment disabled, the block repair
uses one serial tail restart, preserves the prefix, and produces a bounded
solve. The benchmark summarizer now reports executed tail restarts separately
from remaining strict-ready opportunities and no longer classifies a missing
row-tail diagnostic as the first hard blocker.

The same serial suffix restart was then opened to unscaled all-current
postcheck rejects. Those failures have already refreshed the block with the
old pivot order, but the prefix before a non-root rejected pivot is still
replayable by the same live-prefix validator. KLS now tries the suffix restart
before full block repair in that state too. The `fast_repaired_tail_restart_*`
work counters now use the actual serial suffix work KLS executes, while the
`fast_rejected_pivoting_tail_*` fields remain the CKTSO-style ETree-tail
diagnostic target for a future row/segment scheduler.

The serial tail retry then stopped copying the full off-block row/value arrays
into scratch storage. The trial factor only needs mutable `Offp` offsets while
constructing columns; after a repair is accepted, KLS rebuilds `Offi`/`Offx`
from the final `Pinv` anyway. This removes a full `nzoff` copy from every
serial tail attempt without changing pivot choices or the published numeric
state.

Accepted serial tail retries then stopped doing a full global off-diagonal
rebuild when the already-refreshed prefix can be preserved. KLS now rebuilds the
off-diagonal suffix from the rejected column through the end of the matrix from
the rebuilt `Pinv`, preserving only the proven-good prefix `Offp` count and
falling back to the full rebuild if the suffix shape proof fails. The BTF
off-block tail smoke changes suffix off-block values between the base and
repaired matrices, so stale values in that suffix are observable in the solve.

The KLS-owned pivoted block restart was then opened to root-of-block rejects.
The serial suffix path still requires a reusable prefix and therefore does not
count root rejects as `fast_tail_restarts`, but `local_reject == 0` can now run
the same sparse pivoting kernel over the whole rejected block before considering
the KLU block-kernel fallback. `fast_kls_block_restarts` records these KLS-owned
block repairs separately from total block restarts, and the prior-pivot smoke
case now requires a root reject to use that path. This is a direct step toward a
KLS-owned factorization kernel, but not the CKTSO ETree-descendant pipelined
tail executor.

KLS now also gives the independent row-up-looking factorization kernel a chance
when a checked reject cannot be recovered by the exact in-place tail repair.
This moves that recovery branch toward the papers' "switch from invalid checked
refactorization to pivoting factorization" rule: exact block/tail repair remains
first, then KLS rebuilds the rejected BTF block with its row-first
dynamic-column-pivot executor, then it tries the quality-checked KLS-first
whole-numeric rebuild before the KLU block-kernel fallback. The block-local
executor now covers both unscaled and KLU row-scaled repair states by consuming
input-row `Rs` during row construction and permuting it back to pivot order
after an accepted scaled block. The scaled preserved-column refresh now uses
the same convention for masked ETree-tail attempts before this broader rebuild
is considered. The whole-numeric checked-reject recovery now covers unscaled
and KLU row-scaled states; it bypasses the normal automatic KLS-first cost gate
but still honors an explicit
`KLS_ENABLE_KLS_FIRST_FACTOR=0` disable. Successful row-first block repairs are
reported as `fast_kls_block_restarts`; successful whole-numeric recoveries are
reported as `fast_kls_rebuild_restarts`. This still stops short of CKTSO
Algorithm 5's ETree-descendant pipelined tail scheduler in place, but removes
another KLU-kernel step from the checked-reject fallback ladder.

The row-refactor benchmark controls were then made explicit after the
checked-row/refactor ambiguity above was found. `kls_bench` now accepts
`--row-refactor env|off|refactor|checked|all`, `run_bench_suite.py` forwards
the same option, and stats report both retained row-pattern shape and actual
last numeric row-kernel execution. A focused four-thread `G2_circuit` check
confirmed the distinction: `--row-refactor checked` ran one checked row
fast-factor pass but reported `row_refactor_last_run=0` after the subsequent
unchecked refactor, while `--row-refactor refactor` reported
`row_refactor_last_run=1` and `row_refactor_last_checked=0` for the final
refactor. In that same noisy sample the explicit row refactor was slower than
the default EGraph refactor, so this is retained as measurement cleanup and
benchmark reproducibility, not as a default policy change.

The EGraph refactor worker then stopped re-running the unscaled single-block
and large-BTF kernel dispatch inside every column. KLS now records the proven
column-kernel kind in the shared worker state before launching the solver-owned
pool, and the hot cluster/pipeline loops call the selected kernel directly. A
same-session three-pass A/B guard with five refactors per sample improved
median repeated refactor on `G2_circuit` by about `1.4%`, `ASIC_320k` by about
`1.2%`, and `onetone2` by about `6.0%`. This is retained as branch cleanup for
the existing exact-EGraph refactor consumer, not as evidence that dispatch
cleanup can close the CKTSO-scale row/segment gap.

The experimental parallel row-refactor pattern then started retaining a group
execution kind: single-row, generic multi-row, or dense multi-row. The row
scheduler consumes that persistent metadata directly and reuses each group's
precomputed trailing length, so it no longer probes the dense-group path for
groups that were already classified as generic. Same-session A/B checks of the
explicit `--row-refactor refactor` path showed `G2_circuit` essentially neutral
to improved at `0.981x` new/base median repeated-refactor time, `OPF_10000`
neutral at `1.000x`, and `xingo_afonso_itaipu` improved to `0.838x`. This is
retained as row-engine groundwork for the future CKTSO-style row/segment
numeric kernel; the row path remains explicit or environment-gated and is still
not a default solver policy.

The next row-scheduler metadata step retained each row group's external
dependency rows. During the dynamic pipeline tail, generic-only row patterns
can now wait once on those external rows before entering the group kernel,
instead of testing predecessor completion inside every row update. The broader
version that also pre-waited dense row-segment groups was rejected because it
regressed `G2_circuit`; dense segments keep their old in-kernel wait behavior.
With the generic-only gate, the same-session Release A/B on the explicit
`--row-refactor refactor` path improved `G2_circuit` to `0.969x` new/base
median repeated-refactor time, kept `OPF_10000` essentially neutral at
`1.003x`, and improved `xingo_afonso_itaipu` to `0.929x`. This is retained as
CKTSO-style scheduler metadata for the experimental row path, not as a claim
that the default KLU-storage solver has closed the CKTSO row-engine gap.

Two follow-up row-group scheduler metadata shortcuts were rejected. First,
retaining a precomputed per-group work estimate avoided repeated scans in
pipeline-scope reporting and level-slice construction, but same-session Release
A/B on the explicit row-refactor path regressed `OPF_10000` repeated refactor
to `1.157x` new/base and `xingo_afonso_itaipu` to `1.076x`, despite a small
`G2_circuit` win at `0.975x`. Second, retaining the generic-only prewait
eligibility byte avoided recomputing the prewait predicate in the pipeline
loop, but regressed `OPF_10000` to `1.110x` and `xingo_afonso_itaipu` to
`1.067x`, while `G2_circuit` was neutral at `1.001x`. Both probes are too
shallow: the remaining CKTSO gap still points to native row/segment numeric
storage and pivoting-tail restart semantics, not more cached scalar scheduler
decisions around the current KLU-backed row experiment.

The retained row-group reverse graph is now consumed by the experimental row
pipeline tail as a bounded successor-ready queue. Queue setup marks the groups
remaining after the cluster split, counts only tail-local predecessor groups,
starts from zero-predecessor tail roots, and enqueues successor groups when
their tail predecessor count reaches zero. A strengthened smoke fixture uses
three dense row groups with a lower coupling, proving a nonzero group edge and
queued execution while preserving a bounded residual. This is closer to the
SubtreeLU private/pipeline scheduler contract than the previous level-list
cursor, but it still runs the current KLU-backed row numeric updates and is not
the missing CKTSO row/segment kernel by itself.

The queued row tail then stopped allocating and publishing the per-row
completion bitmap for unchecked refactors. The queue already carries the
tail-local group predecessor counts, so the bitmap is only needed by checked
row fast-factor runs that may need prefix-current reject validation. Smoke
coverage now proves the split: the checked queued pass uses the bitmap, while
the later unchecked queued pass reports no last-run bitmap use. This removes a
KLU-wrapper bookkeeping artifact from the experimental row scheduler without
changing default KLS policy.

The row-pattern build then added per-row input-cleanup metadata. For each row,
KLS checks whether every input column is already cleared naturally by the row
numeric pass through an `L` dependency, the pivot, or a `U` entry. Rows with
full coverage skip the old final input-column cleanup loop. The smoke fixture
now proves zero cleanup rows on its dependent row-group pattern, and benchmark
JSON exposes cleanup rows and entries so broad runs can confirm whether the
optimization is structural on larger circuit matrices. This removes another
current-row-kernel bookkeeping pass without tuning on matrix names.

The dense row-segment path then narrowed deferred value scatter to checked
pivot probes. Previously, the presence of any dense row segment forced the
entire row refactor to store into row mirrors and copy all row-major `L`/`U`
values back into KLU storage at the end. The native dense group kernel now also
supports unchecked operation by scattering each completed dense row directly
after its in-group update is finalized. Smoke coverage proves the checked pass
still reports deferred scatter while the subsequent unchecked queued refactor
does not. This removes a broad post-pass copy from the experimental row
refactor without changing default KLS policy.

The row ready queue was then made work-aware using the existing group work
estimate. The retained successor graph is sorted once when the row pattern is
built, and each run sorts the initial tail-ready groups before publishing them
to worker threads. This follows SubtreeLU's queue-balancing motivation of
scheduling heavier refactor work earlier while preserving dependency readiness;
it is still running the experimental KLS row layer, not a full separator-tree
private/pipeline scheduler.

The row-pattern builder now retains the row-to-group map and uses it when a
checked row fast-factor pass rejects a pivot. KLS records the conservative
group-tail restart scope reachable from the rejected row's group: number of
groups, covered rows, and retained group-work estimate. This is a planning and
diagnostic bridge to CKTSO's "rows that need to be recomputed with pivoting"
step; it does not yet execute a parallel pivoting row-tail kernel.

The row ready queue workspace then moved from per-run heap allocation to
solver-owned storage tied to the row-pattern lifetime. The ready group array,
ready-slot atomics, tail-local predecessor counters, and tail-membership bitmap
are now resized only when the retained row task graph grows. This follows the
SPICE repeated-refactor requirement that scheduler metadata be retained rather
than rebuilt from scratch on every Newton step. The same chunk also made
checked queued rejects deterministic by publishing completed generic rows
inside multirow groups and by refreshing any missing rows before the rejected
pivot before the existing prefix proof runs.

Full-graph queued row runs then started consuming cached root groups through a
private-root cursor before falling back to the shared ready queue for
successor-released groups. The completing worker also keeps one ready successor
as a local continuation and spills the rest to the shared queue; benchmark JSON
now reports `row_refactor_last_local_ready_groups` and
`row_refactor_local_ready_group_count`, so broad runs can distinguish "ready
queue active" from actual worker-local continuation use. Focused probes showed
the path is heavily exercised: `G2_circuit` reported 64,962 local continuations
in the last row-refactor run and `rajat24` reported 212,824.

Cluster-prefix row tails now use the same private first-wave treatment when
the initially ready tail groups can be copied into per-thread queues, leaving
later successor releases on the shared/local pipeline queue. This is still a
row-group DAG scheduler, not SubtreeLU's full separator-tree partitioner, but
it removes avoidable shared-queue traffic from the private/pipeline boundary.
Benchmark JSON reports `row_refactor_last_private_ready_groups` and
`row_refactor_private_ready_group_count` for this path.

A follow-up attempt to replace the private-root cursor with static strided
root ownership was tested and rejected. Although this looked closer to a
SubtreeLU private-queue shape, it removed dynamic balancing from the first wave
of root tasks. In same-session focused probes it left `rajat24` essentially
flat (`0.27955s` versus the prior `0.27924s` average refactor) but regressed
`G2_circuit` from the prior `0.423s` range to `0.740s`. KLS therefore keeps the
dynamic private-root cursor until a true separator-tree or FLOP-balanced
private-queue partitioner exists.

The tail-restart summarizer then started reporting serial-suffix overcompute:
for executed local tail restarts it now totals the extra suffix columns and work
above the retained CKTSO-style pivoting-tail plan. This does not change solver
behavior, but it makes broad benchmark output rank the cases where replacing
the conservative suffix fallback with a real pipelined pivoting-tail executor
should remove the most work.

Those overcompute counters are now also published by the solver and benchmark
JSON as `fast_repaired_tail_restart_overcompute_columns` and
`fast_repaired_tail_restart_overcompute_work`, with the standalone summarizer
kept backward-compatible for older JSONL runs. This keeps the CKTSO-tail
executor target visible in every benchmark row instead of requiring a separate
postprocessing calculation.

The retained CKTSO-style pivoting-tail plan then gained an explicit executor
shape classification. Benchmark JSON now reports whether the ordered
ETree-descendant restart set is contiguous, whether it exactly matches the
contiguous suffix that KLS's current serial restart can execute, how many
columns are gaps inside the retained plan, and how much extra suffix work the
current serial executor would do above that plan. These fields do not execute
the non-contiguous tail; they make broad paper-suite runs identify the cases
where KLS needs the real CKTSO pipelined row-tail executor instead of another
safe KLU-style suffix cleanup.

KLS then started consuming a conservative executable subset of that retained
tail plan. When the ETree-descendant plan begins at the rejected pivot and its
last planned column is before the full block suffix end, the KLS-owned pivoted
block repair first tries to refactor only the contiguous envelope covering that
plan and preserve the later columns. Internal envelope gaps that are not in the
ETree-tail plan are locked to their old pivot rows, so a gap column cannot
introduce a new pivot change merely because KLS is using a serial envelope
rather than CKTSO's true non-contiguous pipeline. If pivoting selects a row
whose old pivot position lies outside the envelope, if a locked gap changes
pivot, if preserved `L` rows cannot be remapped to the new final pivot order,
or if the speculative attempt fails, KLS restores the original block pointer
metadata and falls back to the existing serial suffix restart. This is still
not the full CKTSO non-contiguous pipelined tail executor, but it is now an
executable tail-envelope approximation rather than only a diagnostic.

The pivoting-tail plan then stopped treating every prefix-current checked
reject as a full suffix when a narrower restart seed is available. KLS first
uses the retained row-refactor tail for checked row-major rejects; otherwise,
when the checked worker bitmap can identify unfinished nodes, it records those
unfinished local columns as `fast_rejected_pivoting_tail_seed_columns` and
closes only that seed set through the ordered-block ETree. The unfinished
bitmap seed is retained separately from the mutable planning queue, so it
remains available if the row-tail seed is rejected by the topological closure
check. This is still diagnostic/planning infrastructure, but it matches CKTSO's
restart-point determination more closely and exposes the true non-suffix
worklist that a pipelined pivoting-tail executor should consume.

The checked refactor-pool path then stopped discarding all unfinished-work
information when it can prove the prefix is current only at BTF-block
granularity. Those rejects now record the rejected block suffix as a
conservative unfinished seed before the ETree tail plan is built, and
benchmark output reports it separately as
`fast_rejected_pivoting_tail_block_seed_columns`. This does not make the pool
path a full CKTSO Algorithm 5 executor: row-tail seeds still take precedence
when they are valid, and the executable repair remains the existing serial
tail envelope. It does, however, preserve the paper's "unfinished guessed
EGraph nodes seed the ETree tail" information for the pool path instead of
collapsing directly to an untagged suffix fallback.

The same block-suffix seed now covers serial prefix-current checked rejects.
Serial row-major checked refactors, serial mapped checked refactors, and the
serial scaled checked continuation all have a precise interrupted set: the
rejected pivot through the end of its diagonal block. KLS now records that set
with `fast_rejected_pivoting_tail_block_seed_columns` before building the
ordered-ETree closure. A validated retained row-tail still wins when it gives a
narrower repair with preserved-gap refresh; otherwise these paths feed the same
conservative CKTSO Algorithm 5 planning fallback as the pool path instead of
clearing the seed and rediscovering an anonymous suffix later.

The retained CKTSO-style tail plan now also materializes the ETree forest
shape that a true pipelined pivoting-tail executor would consume. Once the
seeded worklist has been closed through ordered-block ETree parents and
validated as topological, KLS retains the tail-local parent, child-count, and
leaf-to-root level arrays and reports retained ETree edges, roots, leaves,
maximum fanout, level count, and maximum ready-level width as
`fast_rejected_pivoting_tail_etree_edges`,
`fast_rejected_pivoting_tail_etree_roots`,
`fast_rejected_pivoting_tail_etree_leaves`, and
`fast_rejected_pivoting_tail_etree_max_fanout`,
`fast_rejected_pivoting_tail_etree_levels`, and
`fast_rejected_pivoting_tail_etree_max_width`. Smoke tests require the retained
tail to satisfy the forest invariant `edges + roots == columns` and to expose a
nonempty ready-level structure. This is still scheduler state rather than the
full CKTSO tail executor, but it fills the paper-level representation gap
between a flat suffix/envelope retry and Algorithm 5's ETree-descendant ready
worklist.

The executable exact-mask repair paths now consume that forest representation
as a validation gate. Before the serial narrow-tail retry or the threaded
row-first active-mask repair accepts a retained tail mask, KLS verifies the mask
against the retained scheduler arrays and checks that its roots, leaves, edges,
maximum fanout, levels, and maximum width match the retained forest counters. If
the proof fails, KLS falls back to the existing wider repair path rather than
treating a generic range mask as a CKTSO-style ETree tail. Accepted repairs
report `fast_repaired_tail_restart_etree_mask` in addition to the older
exact-mask counter.

The serial suffix tail retry then stopped copying and mutating a private
`Offp` array. Tail-column construction now supports a discard-only off-block
mode, which is valid because accepted local repairs rebuild `Offi`/`Offx` from
the final `Pinv` before publishing numeric state. This removes one more
whole-matrix scratch allocation from the executable tail-restart path while
preserving the current conservative suffix semantics.

Unchecked row refactors then stopped publishing every completed row back into
KLU's column storage immediately. The KLS-owned row-major `L`/`U` mirrors now
remain authoritative across repeated unchecked row refactors, and a dirty flag
forces a single publish before a KLU solve, transpose solve, later
`kls_factor`, or non-row refactor fallback. This is a direct SPICE-cycle
optimization for the row/segment engine: repeated Newton refactors no longer
pay KLU scatter traffic on every step when the next step can consume the
row-major mirrors directly.
A focused four-thread `G2_circuit` row-refactor probe with three repeated
refactors reported `0.337s` average refactor and a `0.047s` solve that included
the delayed publish, with valid residuals. This is a row-engine throughput step,
not yet the missing exact pivoting-tail executor.

A guarded row-major solve then consumed those dirty row mirrors directly for
unscaled, normal-orientation, single-block, non-transpose solves with no
external KLS row/column scaling or permutation. It does not replace KLU's
general triangular solve; it only skips the lazy publish when the KLS row
mirrors are already authoritative. Focused checks kept residuals valid:
`G2_circuit` with `--row-refactor all` reported three row solves and reduced
the non-transpose solve average to about `0.020s` before the benchmark's
transpose solve forced the expected publish, while an unscaled `nxp1`
row-refactor probe used the row solve once with about `0.029s` solve time.

The pivoting-tail repair path then moved one step closer to CKTSO Algorithm 5.
When the retained ETree-descendant tail is non-contiguous, the preferred KLS
block repair now tries to copy already-finished gap columns from the current
numeric LU stream and numerically refactor only the marked tail columns. A gap
column is copied only if its pivot row is still unchanged, its L rows are still
unpivoted, and its U predecessors do not include a recomputed tail column. If
that copy is impossible because the gap depends on recomputed tail state, KLS
now promotes the gap into the active tail mask and recomputes it with normal
pivoting instead of forcing the old pivot inside a non-tail gap or immediately
falling back to the wider serial restart. Structural/mapping copy failures still
fall back. Benchmark JSON reports copied gap columns/work as
`fast_repaired_tail_restart_skipped_columns` and
`fast_repaired_tail_restart_skipped_work`, and
`fast_repaired_tail_restart_exact_mask` identifies successful repairs whose
executed columns exactly match the ETree-derived pivoting-tail mask; promoted
gap repairs deliberately clear that flag and show up as tail overcompute. This
is still a serial conservative subset of CKTSO's pipelined pivoting-tail
executor, not the full non-contiguous parallel row-tail algorithm, but successful
masked cases are no longer blended with serial gap refactors and dependency
blocked internal gaps no longer force an all-or-nothing suffix fallback.

The row-up first-factor dynamic column-pivot path then stopped relabeling
swapped columns by scanning every previously emitted U entry. U entries now keep
per-column linked lists, so a dynamic column exchange only touches entries in
the two swapped columns. This directly targets the CKTSO Algorithm 1 pivoting
step in the KLS-owned first-factor scaffold. Focused checks stayed valid: the
small dynamic-pivot smoke still exercises the row-up path, `nxp1` KLS-first
initial factor time dropped from about `6.37s` to about `3.72s` with 494
dynamic column pivots, and `rajat24` moved from about `47.2s` to about
`44.2s` with 2355 dynamic column pivots. The modest `rajat24` change confirms
that relabeling was not the dominant missing paper mechanism there; KLS still
needs a parallel row-up/EGraph first-factor executor rather than only cheaper
serial pivot bookkeeping.
The refactor gap remains; this is a storage-ownership bridge toward the
row/segment engine, not the missing CKTSO pivoting-tail executor.

The experimental row-refactor numeric path was then opened to KLU row-scaled
single-block factors. The row path now recomputes KLU row scales before the
serial or threaded row pass, divides input entries by the unpermuted row scale
using the same fixed input-position mapping as the EGraph refactor path, and
permutes `Rs` back to pivot order only after an accepted pass. Checked rejects
leave `Rs` in input-row order for the existing scaled repair/tail machinery, and
the row-tail candidate diagnostic uses the same scaled input loader. A smoke
case covers both one-thread and four-thread scaled row refactors, validates the
constructed solution, and now confirms that the scaled dirty row mirrors remain
authoritative through the normal solve by dividing the RHS through pivot-order
`Rs`, matching KLU's `P*(R\b)` solve setup. This is a general row/segment-engine
coverage step, not the full CKTSO ETree-descendant pivoting-tail executor.

The dirty row-major solve was then extended to transpose solves for the same
single-block eligibility class. It loads `Q' * b`, scatters along row-major `U`
for the `U'` solve, scatters backward along row-major `L` for the `L'` solve, and
writes through `Pnum`, dividing by pivot-order `Rs` when KLU row scaling is
active. The scaled row-refactor smoke now solves both forward and transpose
systems from dirty row mirrors and checks that no lazy publish occurs. This keeps
the benchmark-style transpose validation from forcing KLU column storage after a
row refactor, while external KLS row/column scaling and permutation still use the
publish-and-KLU fallback.

The row-refactor pattern and numeric pass were then opened from single-block
factors to BTF diagonal blocks. The retained refactor map now supplies only
diagonal-block input entries to the row-major update, translates KLU's local
block `L`/`U` row indices into global row order, and refreshes BTF off-block
`Offx` values separately from the original input positions. A smoke case covers
a reducible BTF matrix whose off-block entry changes across refactor, verifies
that the row refactor ran, and initially confirmed that BTF solves published the
dirty row mirrors before using KLU's general triangular solve. This is still a
prerequisite for CKTSO-style row/segment tails rather than the full
pivoting-tail executor.

The dirty row-major solve was then extended across BTF factors. Forward solves
now process BTF blocks in KLU order from last to first, solve each diagonal
block from KLS-owned row-major `L`/`U` mirrors, and subtract refreshed off-block
`Offx` columns from earlier block rows. Transpose solves process blocks from
first to last, apply `Offx'`, and then run the row-major `U'`/`L'` triangular
passes inside the block. The BTF row-refactor smoke now checks both forward and
transpose solves and verifies that dirty row mirrors remain authoritative
without publishing back to KLU storage. This is a storage-ownership step toward
the CKTSO row/segment engine, not the full pivoting-tail executor.

The dirty row-major solve then switched from one-RHS-at-a-time sparse traversal
to KLU-style chunks of up to four right-hand sides. The same BTF and
single-block row-mirror paths now load, update, and store a small RHS batch
while traversing each `L`, `U`, and `Offx` pattern once per chunk. The BTF smoke
now solves five forward and five transpose right-hand sides for both unscaled and
KLU row-scaled factors. This is a structure-adaptive solve-path cleanup enabled
by KLS-owned row storage; it does not replace the missing pivoting-tail
factorization engine.

The auto symbolic selector then started comparing SCOTCH as a guarded
nested-dissection candidate instead of keeping it explicit-only. KLS only pays
for this comparison on large single-block patterns with high estimated symbolic
work, and only keeps SCOTCH when its symbolic fill/work score is at least a
material win over the current AMD/COLAMD/METIS candidate. This follows the
CKTSO/SubtreeLU combined-ordering motivation without adding matrix-name
tuning, and it still leaves the row/segment pivoting-tail executor as the main
missing CKTSO-scale mechanism.

The checked row-tail diagnostic was then tightened for BTF rejects. When the
retained KLS row-successor graph can build a tail from the rejected global row
and every row in that tail stays inside the rejected BTF block, KLS now reports
that row-owned tail directly instead of first falling back to the older KLU
numeric `L` scan. If the retained graph is missing or crosses the block
boundary, the KLU scan remains the conservative fallback. This does not execute
CKTSO's pipelined pivoting tail, but it makes the available row/segment
worklist explicit for multi-block checked-row failures.

The KLS-owned pivoted block kernel was then wired into an experimental first
factorization scaffold behind `KLS_ENABLE_KLS_FIRST_FACTOR=1`. The scaffold
allocates the KLU-compatible numeric object itself, handles singleton BTF blocks
through `Udiag`/`Pnum`, runs the KLS-owned pivoted kernel for multi-column BTF
blocks, rebuilds `Pinv` and `Offp/Offi/Offx`, seeds KLS-owned row-major value
mirrors from the accepted numeric object, and reports
`initial_factor_path:"kls_first"` when it succeeds. The smoke suite covers a
matrix with a 2-column BTF block so this path cannot pass by singleton handling
alone, and it now requires subsequent forward and transpose solves to consume the
seeded row mirrors. The same KLS-owned pivoted tail now consumes the retained
factor-order input map when that map is valid for the current BTF block, so it
can iterate the already split in-block input slice instead of remapping every
original CSC row for each tail column. `kls_stats`, `kls_bench`, and the gap
decomposition script expose this through `kls_tail_last_mapped_columns` and
`kls_tail_mapped_column_count`; the KLS-first smoke case requires the mapped
tail to be exercised.

The scaffold then gained a more direct CKTSO Algorithm 1 bridge instead of only
wrapping KLU-style columns. Under the same `KLS_ENABLE_KLS_FIRST_FACTOR=1` gate,
eligible no-scale or KLU row-scaled matrices first try a KLS-owned sparse
row-major up-looking first factor with pivot checks over every BTF diagonal
block: each factor row scatters the permuted in-block input row, divides by
KLU's input-row `Rs` when row scaling is active, applies already computed
row-major U updates, checks the diagonal against the remaining U-row maximum,
records row-major L/U entries, and only then packs the accepted local block
factors into KLU-compatible numeric storage for the existing solve/refactor API.
After acceptance, scaled factors permute `Rs` through `Pnum` to preserve KLU
solve semantics. The bridge also implements the direct Algorithm 1 pivot
exchange for this eligibility class: when the diagonal fails the threshold
against the largest active U-tail entry, KLS swaps the current block-local
factor column with that entry, updates previously computed row-major U column
labels, publishes the accepted `Q` permutation, and continues rather than
falling straight back to KLU-compatible block tails. `kls_stats`, `kls_bench`,
and the gap
decomposition script report this through
`kls_first_last_row_uplooking_columns` and
`kls_first_row_uplooking_column_count`, and dynamic exchanges through
`kls_first_last_dynamic_column_pivots` and
`kls_first_dynamic_column_pivot_count`; smoke tests require both the no-exchange
row-major path, a weak-diagonal dynamic column-pivot case, a reducible two-block
BTF case with an off-diagonal coupling, and explicit scaled cases. This is the
first direct row-oriented first-factor kernel in KLS. It is still limited: it
falls back to the older KLS/KLU-compatible first-factor paths on unrecoverable
pivot rejection or unsupported scaling/static-pivoting state. The full CKTSO
production target still needs ETree cluster/pipeline scheduling for first
factorization and the ETree-descendant pivoting-tail executor.
When those mirrors were seeded by `kls_first`, unchecked `kls_refactor` now
automatically attempts the existing KLS-owned row-major refactor path even when
`KLS_ENABLE_ROW_REFACTOR=0`, so the scaffold drives the next SPICE-style numeric
update through KLS row storage instead of immediately returning to KLU's refactor
kernel. Checked fast-factor `kls_factor` calls now use the same automatic
row-major ownership when the mirrors remain current, even with
`KLS_ENABLE_CHECKED_ROW_REFACTOR=0`. The smoke suite forces both row-refactor env
gates off and verifies automatic unchecked and checked row updates plus
forward/transpose solves after each update.
The guarded row-major solve is now allowed when analysis selected internal
transpose orientation as well; solve dispatch already passes the required
internal `kernel_transpose` flag, so auto-oriented benchmark runs can consume
KLS-owned mirrors instead of publishing back to KLU solely because the stored
pattern is transposed. `kls_bench` and `run_bench_suite.py` expose
`--kls-first-factor env|off|on` so this experimental KLS-owned first-factor
path can be compared reproducibly across manifest chunks.
After a KLS-first factor has established automatic row-major ownership,
successful checked fast-factor pivot repairs now reseed the row-major mirrors
from the repaired numeric object when the repair path had to rebuild local LU
storage and drop stale mirror metadata. The guarded solve/refactor lifecycle
therefore stays on KLS-owned row storage after a local repivot, matching the
CKTSO paper's row-major fast-factor/recompute flow more closely even though the
full pipelined ETree pivoting-tail executor is still not implemented.
The automatic KLS-first row-refactor handoff then gained a static work gate:
explicit row-refactor controls still force the row engine for experiments, but
the automatic path now compares the retained row/group work estimate with the
exact EGraph refactor work estimate and skips row refactor when the row plan is
already more expensive. This follows the NICSLU/SubtreeLU recommendation to
select parallel kernels from structure and FLOP evidence rather than matrix
names, and prevents KLS-first scaffolding from replacing a cheaper existing
EGraph refactor with a slower row-major mirror update on hard rows.
The gate is now visible in `kls_stats` and benchmark JSON through
`row_refactor_total_group_work`, `row_refactor_auto_enabled`,
`row_refactor_auto_values_ready`, `row_refactor_auto_work_allowed`, and
`row_refactor_auto_should_run`, so future KLS-first comparisons can explain
whether the row-major path was skipped because the paper-style work model
rejected it.
A same-session 20-matrix CKTSO-gap probe compared this gated KLS-first mode
with `--kls-first-factor off` at 4 threads and 3 refactors. The geomean ratio
was about 1.01x, with wins on `ASIC_100k`, `rajat15`, `onetone1`, and
`transient`, but regressions on `G2_circuit`, `HTC_336_4438`, `ASIC_100ks`,
and `Raj1`. This confirms that the current scaffold should remain
experiment-gated; the papers point to the full row-major fast factorization and
pivoting-tail machinery as the missing general algorithm, not to blindly
enabling KLS-first on every pattern.
The row-major solve experiment was then decoupled from the KLS-first factor
experiment. `KLS_ENABLE_ROW_SOLVE_FROM_NUMERIC=1` now seeds KLS row-major
`L`/`U` solve mirrors from successful ordinary numeric factors/refactors when
the factor has no external KLS row/column permutation or scaling. `kls_bench`
and `run_bench_suite.py` expose this as `--row-solve env|off|on`, and a smoke
test forces ordinary `klu_first` factorization while verifying both forward and
transpose solves use the KLS row-solve path. This keeps CKTSO's row-oriented
solve idea measurable without implying that KLS is delegating to another solver
as a backend.
The row-solve path then gained a single-RHS scalar loop and cached solve
metadata validation. This matches the row-by-row triangular solve form in
CKTSO's Algorithm 2 for the common SPICE one-right-hand-side case, while
preserving the existing four-RHS batched path for wider solves. On the local
five-matrix smoke suite with `--row-solve on`, the geomean moved from about
`0.0161s` before the scalar path to about `0.0123s` after scalar solve plus
validation caching; the `--row-solve off` reference remained about `0.00946s`,
so row solve remains experiment-gated rather than default.
The solve-only seed then stopped building the full row-refactor group/segment
scheduler. It now builds just the row-major `L`/`U` pattern and value mirror
needed by triangular solve, while KLS-first and auto row-refactor paths still
request the full row-refactor pattern. On the same five-matrix smoke suite, the
`--row-solve on` geomean moved again to about `0.0111s`, and solve-only rows
report zero row-refactor groups/segments while still reporting ready row solve
mirrors.
The lean solve setup now also records CKTSO-style triangular partition
diagnostics without allocating the row-refactor scheduler: lower/upper dense
tail start, rows, and entries, plus the fixed eight trapezoid slices used by
the CKTSO paper. The dense tail criterion follows the paper's setup rule: at
least 70% of triangular entries and at least 300,000 entries in the suffix.
KLS now materializes those slice boundaries internally and reports the maximum
per-slice entry count for lower and upper triangular solves, so the parallel
rectangular-slice executor can gate on measured slice balance.
It also precomputes the CKTSO row segmentation step for dense-tail rows: lower
rows split at the slice start, and upper rows split at the slice end, yielding
rectangular versus within-slice triangular entry counts for both factors.
For one-RHS, single-block normal solves, KLS now runs the rectangular part of
each trapezoid slice in the persistent worker pool, with worker barriers around
the sequential triangular piece. It only enables that executor when the
rectangular entries reach the paper's 300,000-entry dense-tail work scale; the
`bcircuit` structure falls back to scalar row solve, while `G2_circuit` crosses
the gate and records parallel slice runs. Rectangular rows inside each slice are
partitioned by accumulated rectangular nonzeros, following CKTSO's thread
workload assignment rule, and diagnostics report the max per-thread rectangular
entries for lower and upper factors.
The next CKTSO Section V gap was the sparse triangular block before the dense
tail: a local `G2_circuit` probe showed that KLS could parallelize about 3.27M
rectangular entries but still left roughly 9.9M lower/upper prefix plus
triangular-piece entries serial. KLS now levelizes those sparse lower/upper
prefixes, runs wide levels in CKTSO-style cluster mode, and lets thread 0 solve
the remaining narrow levels sequentially. The cluster cutoff uses the same
`#threads * 2` width rule already used in the CKTSO-inspired row-refactor
cluster/pipeline split, and the JSON/text diagnostics report sparse level
counts, cluster levels, max widths, and sparse-level run counts.
Validation then showed that the straightforward pthread-barrier executor still
does not close the solve-time gap by itself: `bcircuit` and `G2_circuit` are
too small or synchronization-heavy, and `G3_circuit` ran the sparse-level path
but moved from about `0.199s` scalar row solve to about `0.213s` parallel row
solve. KLS therefore keeps the implementation but only activates it when the
parallelizable work is a substantial share of total solve work and the average
work per synchronization reaches the paper's 300,000-entry dense-tail scale.
These fields are visible in `kls_stats`, `kls_bench` JSON/text, and the gap
decomposition script, giving the parallel triangular-solve path a
structure-based gate instead of a matrix-name heuristic.
That cost model now also controls seeding row-solve structure from ordinary
numeric storage. Earlier `--row-solve on` probes copied `L`/`U` values into
row-major mirrors after each successful ordinary factor/refactor even when the
parallel row-solve executor never ran, adding an `O(nnz(L+U))` tax to the
slow repeated-refactor cases. The adaptive seed first builds the cheap
CKTSO-style partition diagnostics and only prepares row-solve values when the
predicted parallel row-solve work is large enough. On the local top-12 hard-case probe,
forced row-solve seeding had a `3.1772s` geomean, adaptive seeding had a
`3.0613s` geomean with no failures, and explicit `--row-solve off` was
`3.0426s`. This is a useful cleanup, but it is not the missing CKTSO numeric
engine.
The ordinary-numeric row-solve seed no longer performs that value copy even
when the gate accepts. It now keeps the row solve's structural metadata and
reads factor values through the same retained `double *` slots already used to
scatter KLS row-refactor mirrors back into KLU storage. A temporary
`4000 x 4000` lower-triangular probe with `--row-solve on` crossed both
CKTSO-style dense-tail gate terms and reported `row_refactor_last_row_solve=1`
with no row-refactor groups, proving the solve can run directly from ordinary
numeric storage. This removes an `O(nnz(L+U))` setup copy from the experiment,
but it still leaves the larger paper gap in the factor/refactor numeric engine.

The CKTSO Section V solve executor then gained the transpose-side equivalent
instead of using row-scatter loops for transpose solves. KLS now builds retained
transposed row views for `U^T` and `L^T`; each entry points back to the
authoritative KLU-compatible factor slot, so dirty row mirrors and ordinary
numeric value pointers keep the same ownership rule. Single-block, single-RHS
transpose solves use the same sparse-prefix levelization, dense-tail trapezoid
slices, rectangular/triangular split points, worker-pool dispatch, and
work/share gates as normal solves, but execute as dependency-gather triangular
solves to avoid concurrent scatter races. A `4000 x 4000` natural upper
triangular smoke fixture now solves `A^T x=b` through this path and requires a
recorded parallel row-solve run. This closes the obvious normal-vs-transpose
triangular-solve mismatch in the CKTSO Section V implementation; it does not
complete the separate compact-panel transpose group executor.

The same rule now applies to KLS-first row-refactor mirrors. KLS-first still
builds row-refactor metadata so `row_refactor_total_group_work` and
`row_refactor_auto_work_allowed` stay visible, but it copies values into
row-major mirrors only when the retained row-work model allows the row engine
to run. On the local top-six hard-case slice, this removed an unused setup copy
on the rows where `row_refactor_auto_work_allowed=0`; the same-session default
probe moved from `7.7455s` to `7.5856s` geomean. A top-12 guard with
`--row-solve on` completed all rows at `3.0847s`, within noise of the previous
adaptive `3.0613s` run and still faster than the older forced-seed `3.1772s`
baseline.
The row-major solve eligibility was then widened for static-matching factors.
The pre/post solve wrapper already maps normal solves through `R * P * b` and
returns `C * y`, and maps transpose solves through `C * b` and `P' * R * y`.
That means dirty KLS-owned row mirrors and adaptively seeded mirrors from
ordinary numeric storage can be used without publishing to KLU just because
the analysis selected the MC64-adjacent external row permutation or matching
equilibration. A smoke case now refactors a scaled pre-static permutation
matrix and requires both normal and transpose solves to consume the dirty
row-major mirrors.
The row-refactor scheduler was also corrected to preserve CKTSO's
cluster/pipeline split. It now runs the wide cluster prefix chosen by the
paper's `2 * threads` width rule before building the successor-ready queue for
the remaining narrow tail; the previous full-graph ready queue sent very wide
ASIC/G2 levels through hundreds of thousands of tiny atomic tasks. A forced
top-six row-refactor probe improved only from `11.4507s` to `11.3434s`
geomean, and remained much slower than the EGraph path. This confirms that the
current row engine still lacks the paper's production row/segment numeric
kernel and that automatic row-refactor activation should remain cost-gated.
The ready-queue row-refactor tail now also gives each worker a private initial
ready-root range balanced by cached row work before dependent successors enter
the shared queue. This is a small SubtreeLU-style private/pipeline scheduling
step that reduces root-queue contention without changing the row numeric
formulas or claiming the full separator-tree scheduler.

That scaffold was then extended to KLU row-scaled first factors. It computes
`Rs` in input-row order before constructing singleton and multi-column BTF
blocks, rebuilds off-block values before scale permutation, and finally permutes
`Rs` through `Pnum` for solve semantics. The smoke suite now runs the same
2-column BTF case with no scaling and with KLU max-row scaling.

The KLS-first pivoted-block path no longer allocates a fake LU payload before
root block factorization. The shared pivoted-block kernel now accepts an empty
block when the rejected prefix is zero and estimates its own initial LU memory;
non-root tail restarts still require reusable prefix LU. This removes a
wrapper-style placeholder from first factorization while preserving the existing
block-repair semantics.

Re-reading the local CKTSO, NICSLU, and SubtreeLU references leaves one clear
large missing part for the slow cases: KLS still does not own a complete
row/segment-oriented numeric factorization and refactorization engine. CKTSO's
documented advantage is not only METIS-style ordering or MC64-style matching; it
uses row-major sparse LU, EGraph/ETree cluster-pipeline scheduling, fast
factorization with pivot checks, and ETree-descendant tail factorization when a
pivot check fails. SubtreeLU pushes the same direction through separator-tree
private/pipeline queues and row/supernode updates. KLS has implemented the
ordering, LGPL-compatible matching/scaling boundary, diagnostics, and some
row-solve/refactor scaffolding, but the core repeated-iteration speedup in the
papers comes from doing numeric update in that owned row/segment schedule rather
than repeatedly adapting KLU-compatible numeric storage.

## Recommended General Work

1. Build a KLS-owned row/segment-oriented numeric engine instead of adding more
   wrapper-level gates around KLU storage. It should preserve enough row-major
   `L`/`U` access to run CKTSO-style no-pivot fast factorization with pivot
   checks, identify the ETree-descendant restart tail after a failed check, and
   later support SubtreeLU-style separator-tree private/pipeline queues. The
   current KLS-owned block-local pivot restart, suffix restart, and large
   single-block no-pivot cluster/pipeline refactor are useful precursors, but
   the target is a pivoting tail restart inside large blocks.
2. Continue turning matching/scaling into a production MC64-equivalent stage,
   but keep it inside the LGPL-compatible boundary: use the BSD-licensed SPRAL
   scaling submodule, system SPRAL, or independent KLS code, not HSL MC64 or
   restricted MC64 copies from other solver trees. The retained SPRAL path now
   helps large weak-diagonal dominant-BTF cases and avoids replacing no-BTF
   ordering wins, but `pre2` still times out, so matching quality alone is not
   the remaining CKTSO-scale gap.
3. Finish the structure-adaptive triangular solve only after the LU storage
   owned by KLS exposes row-oriented or segmented access cheaply enough that
   solve setup and value mirrors do not dominate repeated refactors.
4. Use static symbolic and numeric-cost models to decide whether a parallel
   kernel should run, so KLS avoids matrix-name-specific tuning. The NICSLU
   R1/R2 counters are now reported and can seed KLS-owned row/segment refactor
   preparation for scheduled numeric states. Benchmark artifacts distinguish
   the model recommendation, attempted preparation, and accepted preparation
   with `row_refactor_auto_model_recommended`,
   `row_refactor_auto_model_attempted`, and
   `row_refactor_auto_model_accepted`; the remaining work is to turn the same
   row/segment substrate into the full pivoting factor and tail-restart
   executor described by CKTSO/SubtreeLU.

The automatic KLS-first row-refactor handoff now applies that last rule before
building the full row-refactor pattern: it scans retained `L`/`U` structure and
input entries to form a conservative lower bound on row-update work, and skips
row metadata setup when even that bound exceeds the exact EGraph refactor work.
On the current 20-row CKTSO-gap focus guard with four threads and three
refactors, this pre-gate skipped row metadata on 11 rows, left nine ASIC/G2/DC
style rows to the exact post-build cost gate, and improved geomean cycle time
by about 2.6% versus the previous copy-gate artifact. This is useful static
cost discipline, but it also confirms the paper reading above: the worst
ASIC/G2/mc2 rows still need the larger CKTSO/SubtreeLU row/segment numeric
engine rather than another wrapper-level policy tweak.

The next forced-row measurements made that distinction sharper. On the top six
current CKTSO-gap rows, explicit `--row-refactor refactor` produced an
`11.31s` geomean cycle time versus `7.55s` for the auto-gated default, with
ASIC/G2 rows regressing by roughly `1.36x-2.20x`. Two same-session experiments
were rejected. First, specializing row input and off-block refresh for unscaled
values regressed the forced-row geomean to `11.75s`; it helped `ASIC_320k` by
about `2%` but lost `7%-11%` on `ASIC_320ks`, `ASIC_100ks`, and `G2_circuit`.
Second, running unchecked row refactors as a full retained group-DAG ready queue
instead of preserving the CKTSO-style cluster prefix regressed geomean to
`11.59s`. The full ready queue did improve `onetone2`, but the large rows
created tens of thousands of local ready continuations, so scheduler churn
outweighed barrier savings. Both tests were reverted. The retained lesson is
that the hard gap is inside the row/segment numeric representation and update
kernel, not a missing unscaled input branch or a blanket full-ready scheduler.

The saved large-paper reconnaissance pair tells the same story at larger scale.
With failures scored at `120s`, KLS was slightly ahead of the saved CKTSO large
artifact in geomean (`0.986x` candidate/reference), because it wins
`TSOPF_FS_b39_c30` and `rajat29` and both solvers struggle on `Hamrle3`.
However, the median ratio was still `1.18x` against KLS and the material losses
were refactor-heavy rows: `nxp1` (`1.67x`), `G3_circuit` (`1.56x`),
`rajat30` (`1.42x`), and `ASIC_680k` (`1.37x`). This keeps the next useful
implementation target aligned with the papers: persistent compact row/segment
numeric storage, batched trailing updates, and pivot-aware tail restart, rather
than another scheduling-only or input-copy-only tweak.
The phase split on those same artifacts rules out triangular solve and
rejected-pivot repair as the main explanation for these rows. Each of the four
losses reported `fast_block_restarts=0` and `fast_tail_restarts=0`. KLS solve
time was already faster than CKTSO on all four, while numeric factor/refactor
was slower: `nxp1` refactor `0.3065s` versus `0.1552s`, `G3_circuit`
`11.03s` versus `7.07s`, `rajat30` `0.2243s` versus `0.1375s`, and
`ASIC_680k` `0.0499s` versus `0.0309s`. Initial factorization was also much
slower on the same cases (`3.13s` versus `0.63s`, `44.16s` versus `5.28s`,
`2.21s` versus `0.58s`, and `1.84s` versus `0.16s`). This is why the local
CKTSO and SubtreeLU papers now point KLS toward the row-major up-looking
factor/refactor kernel with supernode/segment updates first; the full
ETree-descendant pivoting-tail restart remains necessary for stressed
repivoting cases, but it is not what explains these no-reject large losses.

As a small row-engine step after re-reading the CKTSO refactorization section,
the serial row-refactor path now executes the retained row-group processor with
a single local worker instead of keeping a separate scalar row loop. This makes
the dense/generic row-segment kernels available to one-thread refactors and to
single-thread checked fast-factor attempts, preserving the same prefix-scatter,
singular, and rejected-pivot bookkeeping. The smoke suite now asserts that the
dense checked-row prefix repair is a serial checked row-refactor over a dense
segment with deferred scatter. This is still only an incremental move toward
the paper target: KLS remains tied to KLU-compatible value storage for these
rows, so the larger missing item is still persistent compact row/segment
numeric storage and batched trailing updates.

The next CKTSO/SubtreeLU-aligned row-engine step made the compact scratch-panel
path active for unchecked dense row groups when the estimated internal dense and
shared trailing update work crosses a structural gate. The dense-group
dispatcher now tries the compact row-major panel before the native direct
row-mirror kernel for those unchecked groups, falls back to native if the
scratch panel cannot be allocated, and reports
`row_refactor_last_compact_dense_panel` plus
`row_refactor_compact_dense_panel_count`. A generated 48-by-48 dense
unchecked smoke case asserts that the compact-panel path is selected and keeps a
small residual.

A follow-up fixed the checked compact-panel semantics instead of keeping checked
dense groups gated to the native path. The failed trial showed that computing
the whole dense panel before testing any checked pivot can touch suffix state
that CKTSO-style prefix/tail repair expects to remain row-ordered. KLS now uses
a row-ordered compact checked path: after gathering the panel it updates one
row, publishes that row into the row-major mirrors, checks the pivot, and stops
immediately on reject. The dense checked-prefix smoke case now asserts both the
compact-panel marker and the tail repair residual. This is still not the full
pipelined CKTSO tail executor, but it closes a concrete semantic gap for
compact row-major segment updates with pivot-aware prefix repair.

The compact-panel gate was then tightened from a tiny absolute work floor into
a structural arithmetic-intensity check. A dense group now uses worker-local
compact scratch only when its estimated internal dense plus shared-trailing
update work is large enough overall and large enough per copied panel entry.
This avoids treating shallow panel packing as a SubtreeLU/CKTSO mechanism when
the copied scratch has too little update work to amortize it. Focused forced-row
checks stayed valid: `G2_circuit` still exercised 662 dense segments but compact
panel uses dropped from 1280 to about 500 across two refactors, while `mc2depi`
dropped from 3623 to about 1340 compact uses and improved in the same-session
sample. The retained lesson is still that KLS needs persistent compact
row/segment numeric storage, but the scratch bridge now follows a broader
work-per-byte rule instead of a low absolute threshold.

The compact-panel gate is now visible in benchmark artifacts. `kls_stats` and
`kls_bench` report compact-panel eligible dense-group count, eligible rows,
estimated update work, and copied panel entries in addition to the last-run
execution counter. This separates three paper-relevant cases in future suite
runs: no executable dense row segments, dense segments that are too shallow for
compact scratch, and dense segments whose arithmetic intensity is high enough to
exercise the SubtreeLU-style compact-panel bridge. The synthetic dense checked
and unchecked smoke cases now assert those eligibility counters, while small
`add20` still reports zero eligible compact panels.

The gap decomposition CSV now carries the same compact-panel eligibility and
execution counters. This keeps the medium/large CKTSO comparison workflow
aligned with the row-segment diagnostics, so future slow-row reviews can see
whether a refactor-heavy loss has no dense segment work, dense work rejected by
the compact arithmetic-intensity gate, or compact-panel execution that is still
too slow because KLS lacks persistent row/segment numeric storage.

Re-reading the local CKTSO Section IV algorithms also exposed one remaining
serial bridge in the existing block-repair path: after a prefix-current BTF
block repair, KLS recomputed all later diagonal blocks serially even though
those BTF blocks are independent. The solver-owned refactor pool can now start
at a later BTF block, marks earlier blocks as already current for prefix
classification, and the unscaled prefix-current fast-factor restart path tries
that pool before falling back to the old serial suffix. Benchmark stats report
`fast_repaired_parallel_tail_blocks`, and a reducible smoke fixture asserts that
an early repaired block continues over later BTF singleton blocks through the
pool. This follows CKTSO's "continue tail work in parallel" direction only at
BTF-block granularity; it does not implement Algorithm 5's ETree-descendant
row-tail factorization with pivoting, nor SubtreeLU's separator-tree
private/pipeline row queues.

The compact dense row-segment bridge then moved one step from transient scratch
toward persistent row/segment storage. When a dense row group passes the same
compact-panel arithmetic-intensity gate, the row-refactor pattern now allocates
a solver-owned compact panel slice for that group and reports the retained
groups/entries as
`row_refactor_compact_dense_panel_persistent_groups` and
`row_refactor_compact_dense_panel_persistent_entries`. The compact kernel uses
that retained slice before falling back to worker-local scratch, reports actual
retained-slice execution as
`row_refactor_last_compact_dense_panel_persistent` and
`row_refactor_compact_dense_panel_persistent_run_count`, and the dense
checked/unchecked smoke cases assert both execution and retained panel
consumption.
This still repacks current input values on each numeric pass and scatters back
to row mirrors, so it is not the full CKTSO/SubtreeLU row-major numeric storage
model; it does make the compact row-segment value lifetime solver-owned rather
than worker-scratch-owned.

The row solve path now consumes those retained compact panels directly when a
dense group has a complete valid prefix. Row-aware solve value accessors read
in-group `L`, in-group `U`, and trailing-panel entries from the solver-owned
row-major compact panel before falling back to the older sparse row-value
mirrors, and the public stats report last/cumulative compact-panel solve
values. This is still a storage-boundary step rather than the full paper row
numeric engine, but it removes another forced round trip through scattered
KLU-style value arrays from the rows that already have compact SubtreeLU-shaped
storage.

The one-RHS row solve then moved from per-entry compact-panel reads to
group-level compact-panel execution for complete dense groups. The normal
forward solve validates the dense group's lower suffix layout and walks the
row-major lower triangle directly; the backward solve validates the dense and
shared trailing U layout and walks the retained row-major upper/trailing panel
before dividing by the published pivots. Both paths fall back to the previous
sparse row loops if the panel or row layout is not exact, and stats now report
last/cumulative compact group-solve rows and panel entries. This is still a
serial one-RHS kernel, but it is a direct step from "compact panel as storage"
to "compact panel as triangular-solve executor."

The same compact group-solve executor now covers the normal multi-RHS row
solve chunks. KLS validates the retained dense group layout once per group and
updates the chunk's row-major work slab across all RHS columns before falling
back to sparse row loops for non-exact layouts. The compact-panel smoke fixture
now solves three RHS after the one-RHS solve and requires fresh last-run
compact group-solve rows, so the coverage is tied to actual multi-RHS executor
use rather than cumulative one-RHS statistics. Transpose solves still use the
row-aware per-entry accessors, so the solve-side panel executor is not yet
complete.

The refactor scalar fallback now follows the same storage rule for single-row
producer dependencies. If a completed dense producer row belongs to a retained
compact panel, scalar consumers validate that panel's dense/trailing U-row
layout and apply the update from the row-major panel before falling back to the
scattered row mirror. New stats report last/cumulative scalar rows and entries
served from compact panels, and a focused smoke fixture forces the one-row
producer case that the multi-row compact-supernode update intentionally skips.
This does not add the missing coarse BLAS row-panel engine, but it makes
retained compact panels a broader numeric source for the existing row executor.

KLS now retains METIS `NodeNDP` separator-tree queue metadata instead of
discarding it after ordering. For accepted METIS symbolic analyses, the
`NodeNDP` size tree is converted to a postorder private/pipeline component
sequence: leaf domains are private work components, internal separators are
pipeline components, and accepted ordering positions keep a component map for a
future numeric queue consumer. The stats and benchmark output report analyzed
rows, component counts, private/pipeline row totals, and max component sizes.
A focused `nxp1` analyze-only probe with four threads, METIS ordering, and BTF
disabled reported 414604 analyzed rows, seven components, four private leaf
components, three pipeline separator components, and a 414147/457
private/pipeline row split. This closes the direct state-retention gap against
SubtreeLU's separator-tree setup, but the gap that matters for CKTSO-scale slow
cases remains the numeric consumer: pivot-constrained separator work queues,
FLOP-balanced refactor queues, and row/supernode update kernels are still not
implemented.

The retained separator map is now consumed by one numeric path: the experimental
row-refactor ready queue. When the retained `NodeNDP` map covers the full
factor order and every initial ready row group stays inside one separator
component, KLS assigns those initial ready groups to per-thread private queues
by separator component before dependent successors enter the shared/local ready
queue. Benchmark stats report
`row_refactor_last_separator_private_queue`,
`row_refactor_separator_private_queue_run_count`,
`row_refactor_last_separator_private_components`, and
`row_refactor_separator_private_component_count`. A forced `G2_circuit`
METIS/KLS-first/row-refactor probe with four threads selected the
separator-private path, using four retained components for six initial private
groups, with a valid residual. This is a real SubtreeLU-aligned consumer of the
retained separator tree, but it is still narrower than the paper: it does not
constrain pivot search inside separator-tree subdomains, does not partition all
work by separator FLOP balance, and does not implement the pivoting
first-factor or checked-tail kernel that the slow CKTSO-gap rows still require.

KLS now consumes the same retained separator map more directly for no-pivot row
refactorization through a SubtreeLU Algorithm 6-style queue builder. When the
separator tree covers the row-refactor block, KLS computes retained row-group
work per separator component and subtree, repeatedly moves the dominant
subtree root component to a pipeline queue while returning its child subtrees
to the candidate set, then greedily assigns the remaining subtrees to private
thread queues by work. Row-refactor groups that cross separator-component
boundaries cannot be private-subtree work, so KLS classifies them as pipeline
groups instead of rejecting the whole separator queue. If any private subtree
group depends on a pipeline group, KLS promotes the dependent group into the
pipeline closure and keeps the separator schedule rather than falling back to
the generic ready queue. Benchmark stats report
`row_refactor_last_separator_flop_queue`,
`row_refactor_separator_flop_queue_run_count`,
`row_refactor_last_separator_flop_components`,
`row_refactor_last_separator_flop_private_groups`, and
`row_refactor_last_separator_flop_pipeline_groups`, plus closure promotions via
`row_refactor_last_separator_flop_closure_groups` and
`row_refactor_separator_flop_closure_group_count`. A forced `G2_circuit`
METIS/no-BTF/no-scale/no-pivot row-refactor probe with four threads selected
this queue on all three refactors, using seven separator components, 69,402
private groups and five pipeline groups per refactor, with residual
`3.32e-16` and a repeated-refactor average around `0.20s`. This fills the
paper's FLOP-balanced separator-queue mechanism for the eligible no-pivot row
engine; the remaining SubtreeLU gap is applying the same separator partition
inside checked-tail/pivoting kernels, BTF forests, and deeper supernodal BLAS
numeric storage.

KLS-first row-up factorization now builds and conservatively consumes the
corresponding first-factor queue shape from the retained separator map. For
each covered symbolic block, KLS now first applies the retained separator tree
as an Algorithm 6-style split/collapse: it computes component and subtree work
from permuted block-local row input counts, repeatedly moves the dominant
subtree root into the factor-order pipeline queue, returns child subtrees to
the private candidate set, and greedily assigns remaining subtrees to private
threads by work. That replaces the older direct leaf/private split when the
partition is dependency-safe. The numeric consumer validates the private phase
before threading it: if a private row would read another private thread's
mutable column domain, KLS discards the partitioned queue and falls back to the
legacy retained-component private/pipeline queue. When validation passes, each
thread factors its assigned private rows into local row-up `L`/`U` entries. KLS
now first remaps the block-local row/column order and published `Pnum` to that
private-then-pipeline queue before numeric assembly; this fixes the earlier
unsafe version that tried to execute queue order over the old KLU numeric row
indices. KLS then merges private entries and consumes the pipeline rows with an
Algorithm 3-style atomic task counter. Each pipeline worker now keeps local
scratch, applies all ready private-row predecessors against an immutable
private-U snapshot before waiting for its ordered publish slot, and can consume
already-published earlier pipeline-prefix rows while waiting for intervening
pipeline rows that it does not depend on. Published private and pipeline rows
are also grouped into conservative row-supernodes when adjacent row-major `U`
patterns have the full triangular extension and identical trailing pattern.
Partial pipeline-row updates consume such consecutive ready predecessors, or
the already-finished prefix of a retained predecessor supernode, through one
guarded supernode-run executor. The executor keeps the triangular in-run
discovery behavior, so an update from row `k` can create row `k+1` as the next
dependency before the run proceeds. It also validates the same
dense-prefix/common-trailing shape described for supernodes, solves the
in-supernode prefix in the current sparse row workspace, accumulates the shared
trailing contribution in worker-local scratch, and scatters each trailing column
once. This matches the Algorithm 4 update control flow and scalar compact
numeric shape while still using the current KLS row-entry storage rather than
production BLAS panel storage. The publish step can now accept a
scoped dynamic column exchange for a safely pre-updated pipeline row: while
holding the ordered pipeline lock, it swaps the committed U entries and the
phase-local prefix snapshot, publishes the pivot row, increments a
column-order epoch, and wakes workers. Any speculative suffix row that began
under the older epoch clears its local sparse workspace and recomputes the same
row under the new column order instead of forcing a whole suffix relaunch.
This is still not production BLAS-backed supernode storage or the CKTSO
checked-tail pivoting executor. Benchmark stats expose the last
planned queue through
`kls_first_last_separator_queue`,
`kls_first_last_separator_queue_private_components`,
`kls_first_last_separator_queue_pipeline_components`,
`kls_first_last_separator_queue_private_rows`,
`kls_first_last_separator_queue_pipeline_rows`,
`kls_first_last_separator_queue_nonempty_threads`, and
`kls_first_last_separator_queue_max_thread_rows`, with the structural work
range reported by `kls_first_last_separator_queue_min_thread_work` and
`kls_first_last_separator_queue_max_thread_work`. Algorithm 6-style partition
use is reported through `kls_first_last_separator_queue_partitioned`,
`kls_first_separator_queue_partitioned_count`, and
`kls_first_last_separator_queue_split_components`; scheduled consumption is
reported by `kls_first_last_separator_queue_executed`,
`kls_first_separator_queue_executed_run_count`,
`kls_first_last_separator_queue_executed_private_rows`, and
`kls_first_last_separator_queue_executed_pipeline_rows`; threaded private
execution is reported by `kls_first_last_separator_queue_parallel_private`,
`kls_first_separator_queue_parallel_private_run_count`,
`kls_first_last_separator_queue_parallel_private_rows`, and
`kls_first_last_separator_queue_parallel_private_threads`; threaded pipeline
queue consumption is reported by
`kls_first_last_separator_queue_parallel_pipeline`,
`kls_first_separator_queue_parallel_pipeline_run_count`,
`kls_first_last_separator_queue_parallel_pipeline_rows`, and
`kls_first_last_separator_queue_parallel_pipeline_threads`; private-predecessor
pipeline pre-updates are reported by
`kls_first_last_separator_queue_pipeline_partial`,
`kls_first_separator_queue_pipeline_partial_run_count`,
`kls_first_last_separator_queue_pipeline_partial_rows`, and
`kls_first_last_separator_queue_pipeline_partial_threads`; waiting-prefix
pipeline pre-updates are reported by
`kls_first_last_separator_queue_pipeline_wait_partial`,
`kls_first_separator_queue_pipeline_wait_partial_run_count`,
`kls_first_last_separator_queue_pipeline_wait_partial_rows`, and
`kls_first_last_separator_queue_pipeline_wait_partial_deps`; scalar
row-supernode partial updates are reported by
`kls_first_last_separator_queue_pipeline_supernode_update`,
`kls_first_separator_queue_pipeline_supernode_update_run_count`,
`kls_first_last_separator_queue_pipeline_supernode_update_groups`, and
`kls_first_last_separator_queue_pipeline_supernode_update_rows`; pipeline
pivot-tail restarts are reported by
`kls_first_last_separator_queue_pipeline_pivot_tail`,
`kls_first_separator_queue_pipeline_pivot_tail_run_count`, and
`kls_first_last_separator_queue_pipeline_pivot_tail_rows`, with restart and
actually serialized pivot-row counts reported by
`kls_first_last_separator_queue_pipeline_pivot_restarts`,
`kls_first_separator_queue_pipeline_pivot_restart_count`, and
`kls_first_last_separator_queue_pipeline_pivot_serial_rows`. Prefix
private-predecessor pre-updates, waiting-prefix pre-updates, and scalar
row-supernode updates completed before a pipeline pivot restart remain counted
in the same counters after the epoch retry, while
`kls_first_last_separator_queue_pipeline_pivot_serial_rows` remains reserved
for the older external serialized fallback. After a successful dynamic pivot
inside the separator pipeline, KLS now rebuilds the phase-local panel cache
from the post-exchange column order over the whole committed prefix instead of
disabling the cache for the rest of the phase; suffix rows that restart under
the new epoch can still consume validated dense/common-tail producer panels.
The rebuild is reported through
`kls_first_last_separator_queue_pipeline_prefix_panel_rebuild`,
`kls_first_separator_queue_pipeline_prefix_panel_rebuild_count`, and
`kls_first_last_separator_queue_pipeline_prefix_panel_rebuild_rows`. Smoke
tests cover the queue shape on a 30,000-row METIS-ordered tridiagonal
KLS-first factor and the epoch recovery path on a 30,600-row METIS-ordered grid
with separator-band weak diagonals. The queue-shape test requires more than
one private worker thread, verifies that all planned pipeline rows are consumed
by the guarded pipeline executor, verifies that the private-predecessor partial
pre-update path touches all pipeline rows, and now requires cached
panel-backed producer updates in the separator pipeline. The epoch test
requires a separator-pipeline pivot restart, a prefix-panel rebuild whose
covered rows exceed the private-prefix rows, cached separator-pipeline panel
updates after that restart, zero serialized pivot rows, and a clean solve
residual.
This closes Algorithm 3's queue shape and adds race-free scalar Algorithm 4
partial-update and scoped-pivot retry steps for dependency-safe retained
separator queues.

For KLS-first blocks without a retained separator order, the row-up first factor
now still uses the restartable Algorithm 5-style row pipeline when multiple
threads are available. The executor builds a natural factor-order row list,
lets workers perform partial dependency updates and ordered publication, and
uses the same epoch retry path for dynamic column exchanges. These runs are
reported separately from separator queues through
`kls_first_last_row_pipeline`, `kls_first_row_pipeline_run_count`,
`kls_first_last_row_pipeline_rows`, `kls_first_last_row_pipeline_threads`,
`kls_first_last_row_pipeline_partial`,
`kls_first_row_pipeline_partial_run_count`,
`kls_first_last_row_pipeline_partial_rows`, and
`kls_first_last_row_pipeline_partial_threads`. Dynamic pivot epochs in this
generic row pipeline are now reported separately through
`kls_first_last_row_pipeline_pivot_tail`,
`kls_first_row_pipeline_pivot_tail_run_count`,
`kls_first_last_row_pipeline_pivot_tail_rows`,
`kls_first_last_row_pipeline_pivot_restarts`,
`kls_first_row_pipeline_pivot_restart_count`,
`kls_first_last_row_pipeline_pivot_serial_rows`,
`kls_first_last_row_pipeline_prefix_panel_rebuild`,
`kls_first_row_pipeline_prefix_panel_rebuild_count`, and
`kls_first_last_row_pipeline_prefix_panel_rebuild_rows`, while separator
pipeline epochs stay on the separator-prefixed counters. Smoke coverage
verifies both a clean natural-order execution with no pivot epoch and a
two-thread weak-pivot natural-order execution that rebuilds the committed-prefix
panel cache without setting separator-pipeline counters. This does not add
CKTSO's ETree-descendant pivot-tail executor, but it removes the prior drop to
a purely serial KLS-first row factor whenever the separator queue was
unavailable and makes the generic Algorithm 5-style epoch behavior measurable.

The KLS-first row up-looking dynamic column pivot selector also now consumes
the retained separator map when it is available for the full factor order. On a
weak pivot, KLS now applies the SubtreeLU Algorithm 4 `N'` rule directly:
`N'` is the last retained factor row of the current collapsed component, the
diagonal is compared only with the strongest candidate in `i+1..N'`, and
outside-domain candidates no longer force a separator-scoped rejection through
the old global row maximum. If no safe scoped candidate exists for a required
pivot, KLS rejects the row-up attempt and lets the existing KLS-first fallback
path use the pivoted block kernel. This matches SubtreeLU's requirement that
private-mode pivots do not cross retained component domains. Stats report exact
separator-domain pivots through
`kls_first_last_separator_dynamic_column_pivots` and
`kls_first_separator_dynamic_column_pivot_count`, component-extent pivots
through `kls_first_last_separator_extent_dynamic_column_pivots` and
`kls_first_separator_extent_dynamic_column_pivot_count`, would-have-crossed
fallback candidates through `kls_first_last_separator_dynamic_column_fallbacks`
and `kls_first_separator_dynamic_column_fallback_count`, and strict rejections
through `kls_first_last_separator_dynamic_column_rejects` and
`kls_first_separator_dynamic_column_reject_count`. The existing smoke suite
still exercises KLS-first dynamic column pivoting. This should be treated as a
guarded Algorithm 4 pivot-domain step rather than evidence that the full
SubtreeLU private/pipeline pivoting executor is implemented.

The mapped EGraph fast-refactor path now mirrors the CKTSO checked-reject
prefix recovery that already existed in the row-refactor path. When a checked
EGraph worker rejects a pivot and the pipeline completion bitmap exists, KLS
serially recomputes any unfinished columns before the rejected pivot with the
same mapped EGraph column kernel, marks those columns done, and only then
classifies the reject as prefix-current or unknown. This fills a concrete
CKTSO Section IV gap: a rejected fast factorization no longer loses the valid
prefix just because some earlier EGraph tasks were not scheduled before the
stop flag. Benchmark JSON now reports
`fast_rejected_prefix_refresh_columns` and
`fast_rejected_prefix_refresh_count` so slow-case reruns can distinguish true
unknown-prefix failures from recovered prefix-current ETree-tail candidates.
The implementation is still a serial prefix recovery feeding the existing
serial pivoting-tail/block-repair path; it is not yet CKTSO's full parallel
ETree-scheduled pivoting tail executor.

The compact dense row-refactor kernel then removed a redundant panel sweep in
the unchecked refactor path. Previously the compact dense group first factored
the dense intra-segment panel and then made a second pass over the same
row/dependency pairs to update the trailing panel. The unchecked path now
updates the trailing panel while applying each dependency, matching the
checked compact kernel's dataflow and moving the implementation closer to the
supernode-style dense update direction in CKTSO/SubtreeLU. On the forced
row-refactor focus probe, `G2_circuit` repeated refactor improved from about
`0.507s` to `0.404-0.446s`, `onetone2` improved from about `0.0348s` to
`0.0329s`, and a repeated `ASIC_100ks` probe reported about `0.116s` versus
the previous one-pass `0.148s`, with residuals unchanged. Default automatic
selection remains cost-gated because the row engine is still not broadly
faster than the mapped EGraph refactor.

KLS now maps another explicit SubtreeLU Algorithm 5 detail into the
row-refactor scheduler. When large tail row groups with downstream successors
and width at least `2 * threads` dominate the tail by both group count and row
count, the scheduler prepares a row-dependency ready queue. Dense-group rows
are marked complete immediately after each row is numerically stored and, for
checked runs, after the pivot check passes. Consumers can therefore wait on the
exact finished row prefix instead of the whole supernode. This matches the
paper's large-unfinished-supernode split for no-pivot refactorization and the
checked row fast-factor/refactor path; `KLS_ENABLE_PARTIAL_SUPERNODE_PIPELINE=0`
remains a hard disable. Direct same-tree A/B probes on `coupled` and
`G2_circuit` showed earlier broad selectors were slower than the existing
ready queue on the current KLS row kernel, so the strict structural selector is
kept. It is still not the full CKTSO pivoting-tail factorization.

KLS now consumes retained compact dense row panels as external supernode
update sources during row refactorization. Once a compact dense group finishes
successfully in the current numeric pass, a per-group valid marker lets later
rows recognize a contiguous suffix of two or more dependencies from that same
group and update from the contiguous panel directly instead of expanding each
predecessor row through the sparse row mirror. This implements the
SubtreeLU/CKTSO supernode
update idea more directly than the earlier compact-panel work, which only used
the panel while factoring the producer group itself. The path preserves the
same multiplier checks and row-prefix publication rules as the scalar row
kernel, and falls back to scalar updates when a producer panel is absent,
stale, unfinished, or belongs to the current group. Focused forced-row probes
showed the path active with valid residuals: `coupled` used 1,770 compact
supernode updates over 62,162 predecessor rows, `G2_circuit` used 119,130
updates over 6,220,352 rows, and `ASIC_100ks` used 23,209 updates over
915,915 rows. This is still not full BLAS-backed supernodal factorization or
the CKTSO pivoting-tail executor, but it moves the current row engine from
mere compact storage toward actually consuming supernodes in later updates.

The partial supernode pipeline now closes a direct Algorithm 5-style gap in
that consumer path. Dense and generic producer groups publish the valid prefix
of their retained compact panel as each row is completed, and the ready-queue
scheduler can release successor groups at the SubtreeLU split point
(`width - threads`) instead of waiting for the whole producer group. A consumer
whose dependency suffix ends at the published prefix consumes the compact
panel immediately and uses the existing row-done waits for any unfinished tail
rows it later reaches. New diagnostics report
`row_refactor_last_compact_supernode_partial_update` and the corresponding
partial-update count, row, and entry totals. The smoke fixture forces a
240-row producer with a 2-thread split and requires a 238-row partial compact
supernode update with a residual-clean solve. This is still not CKTSO's full
pivoting-tail executor or a production BLAS supernodal numeric object, but it
fills the specific missing "start consumers from a finished producer prefix"
semantic instead of treating every compact supernode as all-or-nothing.

The scheduler now preserves SubtreeLU Algorithm 6's separator private/pipeline
partition when applying that Algorithm 5 prefix release. A new row-dependency
ready-queue preparation path can reuse the separator private-group mask, so
private separator groups still run through the retained FLOP-balanced queue
while pipeline groups may be released by row-prefix completion when the strict
large-unfinished-supernode test passes. Checked refactorization now uses the
same strict default selector: producer rows are not released until their pivot
checks complete, and the smoke fixture validates that default checked path. A
broader "any large supernode" selector was tested and rejected: on the six-row forced
row-refactor focus probe it became much slower than the separator-only run and
was interrupted after exceeding the normal short-run envelope. With the strict
selector, `build/kls_head_gap_focus6_forced_row_sep_alg5_strict_t4_r3_timeout120.jsonl`
completed all six rows without timeout, but the focus matrices stayed on the
separator queue (`row_refactor_last_partial_supernode_pipeline=0`) and reported
a `15.29s` SPICE-cycle geomean. This closes the scheduling semantic gap without
pretending it is a performance win on the current focus subset.

KLS then removed the remaining dominance selector from the ordinary SubtreeLU
Algorithm 5 row-tail path. The paper's condition is local to each unfinished
dependent supernode: if the supernode has at least `2P` rows, a consumer may
use the completed `k:(k'-P)` prefix instead of waiting for the whole producer.
KLS now selects the row-dependency queue whenever such a dependent
dense/generic producer exists in the pipeline tail; the large producer no
longer needs to dominate the number of tail groups or tail rows. Separator
Algorithm 6 private/pipeline queues still keep their stricter wrapper guard
until private-queue row-prefix release is made fully safe. Private queue groups
are protected from prefix-release enqueueing and receive normal row-edge
decrements. A new non-dominant-tail smoke fixture places one large dependent
producer behind many small dependent groups and requires the partial-supernode
pipeline by default, covering the paper condition that the large producer need
not dominate the tail.

The compact supernode consumer now applies the paper's matrix-vector update
shape for producer trailing panels. For a later row that
depends on a completed dense producer suffix, KLS still computes and checks the
`L` multipliers in row order, but it no longer scatters each producer row's
trailing contribution directly into the sparse work vector. Instead, it stores
the multipliers in worker scratch, accumulates the producer trailing panel into
a contiguous temporary vector, and scatters that vector once to the sparse
columns. This is an in-KLS BLAS-shaped `gemv` dataflow, not an external BLAS
dependency and not CPU-specific tuning. The smoke fixture now separates a dense
producer, lower consumer rows, and later trailing columns so it asserts this
path. Focused forced-row probes showed the accumulator active with valid
residuals: `coupled` accumulated about 1,770 updates over 62,000 rows and
6.6 million trailing entries with residual `8.97e-16`; `G2_circuit`
accumulated about 119,000 updates over 6.2 million rows and 2.0 billion
trailing entries with residual `3.32e-16`; and `ASIC_100ks` accumulated about
23,000 updates over 916,000 rows and 355 million trailing entries with residual
`1.63e-15`. This fills a
direct SubtreeLU supernode-update detail. The exact high-volume diagnostic
counts can vary slightly in threaded runs because these counters are telemetry,
not synchronization state. Full BLAS-backed supernodal
factorization, separator FLOP-balanced queues, and CKTSO's pipelined pivoting
tail executor remain open.

The default compact-supernode trailing accumulator was then tightened to avoid
an extra producer-suffix pass. The scalar suffix solve still uses sparse `x`
because that was faster on ASIC/G2-style probes, but KLS now accumulates the
contiguous trailing vector while each `L` multiplier is already live and only
delays the final scatter. Focused forced-row probes stayed residual-clean and
showed the default `gemv` path active with `trsv=0`: `coupled` refactor about
`0.0050s`, `ASIC_100ks` about `0.0866s`, and `G2_circuit` about `0.238s` in
same-session one-pass samples. This is a production-path cleanup of the
SubtreeLU-shaped update, not another env-only experiment.

The scalar compact-supernode `trsv` path is now the default algorithmic branch
for every ready producer run with at least two rows, instead of only an explicit
experiment or a work-threshold decision. With
`KLS_ENABLE_COMPACT_SUPERNODE_TRSV` unset, KLS copies the producer suffix from
sparse `x` into contiguous worker scratch, solves it with the producer dense
upper panel, and feeds it to the existing trailing-panel accumulator whenever a
ready compact producer run exists. Setting the variable to `0` disables the
automatic choice, and `1` still forces it for probes. Stats report
`row_refactor_last_compact_supernode_trsv`,
`row_refactor_compact_supernode_trsv_count`,
`row_refactor_compact_supernode_trsv_rows`, and
`row_refactor_compact_supernode_trsv_entries`, and the compact-panel smoke
fixture now verifies the unset automatic mode. The partial-prefix smoke now
also forces the same compact `trsv` executor while using SubtreeLU Algorithm 5's
large-supernode split: KLS solves only the completed producer prefix in
contiguous scratch, scatters dense updates that target unfinished producer rows
back to the sparse work row, and later consumes the tail through the normal
row-done waits. This closes the earlier gap where partial-prefix consumers
could start early but had to use scalar per-row producer updates. On the current
forced-row top slice, this
reduced the completed-row geomean from 20.5 s to 18.1 s; `G2_circuit`'s
explicit row-refactor average moved to about 0.236 s with a clean residual.
The default CKTSO-gap focus still keeps row refactor gated off on the large
ASIC/G2 rows, so this is a direct SubtreeLU supernode-kernel step rather than a
claim that the row engine is ready to replace the exact EGraph path.

KLS now also has an opt-in CBLAS supernode experiment. Configure with
`-DKLS_ENABLE_CBLAS_SUPERNODE=ON` and set
`KLS_ENABLE_CBLAS_SUPERNODE=1` at runtime to use standard CBLAS calls for the
row-major compact supernode updates. The consumer-side update now follows the
SubtreeLU text directly for a completed producer supernode: CBLAS `dtrsv`
solves the producer suffix multipliers and CBLAS `dgemv` applies the retained
trailing panel to the current sparse work row, with the same row-order
multiplier checks before publication. This consumer-side BLAS path is gated by
a structural work-per-copied-entry rule, because an ungated `onetone2` probe
issued thousands of tiny CBLAS calls and regressed badly despite producing a
valid residual. Larger gated probes did exercise the external CBLAS consumer
and stayed residual-clean, but they still did not beat the scalar row kernel in
same-session samples: `G2_circuit` was about `2.50s` versus `0.326s`, and
`ASIC_100ks` was about `0.181s` versus `0.156s`.

KLS then added the next, more paper-faithful batch shape: if an unchecked dense
consumer group, or a contiguous row subrange inside it, has the same ordered
list of completed dense producer suffixes as its external dependency pattern,
KLS gathers those consumer rows into the full consumer panel, solves each
producer's consumer-row multipliers as one right-side triangular solve, and
applies each producer trailing panel as one batched panel update before the
consumer panel is factored. The planner now lets an earlier producer update
later producer multiplier columns before those later suffixes are solved, so a
consumer group can batch across multiple completed producer supernodes instead
of requiring a single external producer, and no longer needs every row in the
dense consumer group to share the same external pattern. The original CBLAS
implementation remains an optional backend, but the same planner now has an
LGPL KLS-owned scalar executor for the triangular solves and trailing panel
updates. The smoke fixtures construct two dense producers feeding one dense
consumer and a separate one-producer case where only a 60-row suffix of an
80-row dense consumer group is batchable; both fixtures force the CBLAS runtime
flag off and require batch counters in normal and no-METIS builds. KLS also
reports batch-pattern, batch-candidate, and work-rejected counters so real
matrices can distinguish "no common row subrange" from "candidate too small for
the structural batch gate." A local generated `onetone2_mwmatch` probe after
adding those counters reported 62 compact panels and 17,995 scalar
compact-supernode updates over 417,868 dependency rows, but zero batch patterns
and zero batch candidates. This proves the direct SubtreeLU update shape is
executable in KLS for a broader producer/consumer pattern, while also showing
that the next real-matrix gap is heterogeneous-row producer batching or deeper
row/segment planning, not merely lowering a backend work threshold.

KLS then widened that consumer-side batch shape to the next direct paper case:
contiguous unchecked rows may share one completed dense producer supernode even
when each row starts at a different suffix of that producer. The dense
consumer-panel executor packs those ragged suffixes into a KLS-owned scalar
workspace, solves each row's producer suffix against the retained dense panel,
and applies the producer trailing panel into the current dense consumer panel
before the internal mini-solve. A second executor handles independent
`GROUP_BATCH` rows by keeping each row's pivot and U entries in a temporary row
panel, applying the same ragged producer update across the row batch, and then
publishing the completed row-major values. New smoke fixtures force both
shapes with `KLS_ENABLE_CBLAS_SUPERNODE=0`, so this is not a BLAS backend
dependency. A forced-row top-five CKTSO-gap probe still reported zero batch
candidates on `ASIC_320k`, `ASIC_320ks`, `onetone2`, `ASIC_100ks`, and
`G2_circuit` while retaining hundreds of thousands to millions of scalar
compact-supernode dependency rows. That narrows the remaining paper gap again:
the slow cases are not waiting for identical-pattern or single-producer ragged
suffix batching; they need broader multi-producer, non-contiguous row-panel
numeric execution over the row/segment storage described by CKTSO/SubtreeLU.

KLS also has a default scalar executor for the next contiguous independent-row
paper shape: unless `KLS_ENABLE_MULTI_PRODUCER_SUPERNODE=0` is set, a
`GROUP_BATCH` subrange whose full `L` row is an ordered sequence of completed
dense producer suffixes can use a temporary row panel. It solves each producer
suffix against the retained compact producer panel, applies that producer's
trailing panel into later producer multipliers, the pivot, or the row's `U`
workspace, and then publishes the independent rows as a batch. A smoke fixture
unsets `KLS_ENABLE_MULTI_PRODUCER_SUPERNODE` and validates two dense producer
groups feeding the same independent row batch without using CBLAS. The top-five
CKTSO-gap probe still reported zero batch candidates with this path both
disabled and enabled on `ASIC_320k`, `ASIC_320ks`, `onetone2`, `ASIC_100ks`,
and `G2_circuit`; the scalar compact-supernode update counters remained large.
This rules out contiguous independent multi-producer rows as the large missing
slow-case feature and leaves the direct paper gap at a more structural level:
KLS needs a broader CKTSO/SubtreeLU row/segment panel executor that can batch
heterogeneous, non-contiguous, and internal rows rather than only completed-
producer suffixes exposed by the current row-major storage.

KLS now records that row/segment structure explicitly during row-refactor
symbolic setup. After dense group and compact-panel detection, each row's `L`
pattern is scanned once for maximal contiguous references to completed compact
dense producer groups. The compact-supernode update kernel consumes this
persistent dense-producer run plan before falling back to the old runtime scan.
`kls_stats` and `kls_bench` expose the total planned producer runs, rows,
dependency rows, max runs per row, full-suffix producer runs, multi-run rows,
and fragmented rows. The two-producer smoke fixture requires these counters,
which means KLS now has a concrete symbolic substrate for the CKTSO/SubtreeLU
"for each contributing supernode, update the current row/range" algorithm
instead of rediscovering one suffix at a time inside the numeric kernel. The
remaining gap should be attacked by turning the fragmented and multi-run rows
seen on paper matrices into a broader row/segment panel executor. A forced-row
top-five CKTSO-gap probe after this change reported 4,057-59,571 planned dense
producer runs per matrix over 1,178-9,789 rows, with max 8-27 runs per row and
zero current batch candidates; fragmented rows equaled producer-run rows on all
five matrices. That is a direct paper-algorithm gap: the supernode producer
ranges are present, but KLS still lacks the row-panel executor that can consume
them together with the scalar/non-producer pieces of the same row.

KLS now enables that first fragmented dense-consumer executor by default under
its structural and work gates. Unless `KLS_ENABLE_MULTI_PRODUCER_SUPERNODE=0`
is set for an A/B run, unchecked compact dense consumer groups can batch
consecutive rows whose external prefixes mix scalar dependencies with multiple
completed dense producer suffixes. The executor uses the retained compact
dense input panel when possible, processes scalar gaps only up to the next
dense producer run, then solves and applies each planned producer suffix as a
coarse row-panel stage across the batch before moving to the next scalar gap.
Producer trailing updates can feed later external multipliers or the current
dense/trailing panel, and the completed consumer panel is still left for the
normal internal dense-group factor step. Smoke fixtures unset
`KLS_ENABLE_MULTI_PRODUCER_SUPERNODE`, construct one scalar external row, two
dense producers, and one dense consumer block, and vary producer suffix starts
by row so the older exact-pattern batch path cannot explain the result.

This closes only the structural-dispatch gap, not the performance gap. On the
same five forced-row CKTSO-gap cases, the path fired on all matrices
(`352-2,136` batches and `246,898-1,754,124` batched dependency rows), proving
that the hard rows are dense-consumer fragmented producer rows. The staged
row-panel and compact-direct-input version improved the five-case probe with
this executor enabled from about `33.2s` to about `32.1s`, and improved the
focused `ASIC_320k` refactor probe from about `0.347s` to about `0.319s`. But
the same five-case scalar baseline without this executor was about `13.5s`,
and the focused baseline was about `0.223s`. Promoting the executor to the
default path fills the paper-algorithm dispatch gap; the remaining paper-level
step is still native row/segment panel storage and updates instead of using
this scalar scaffold as if it were SubtreeLU's production supernodal kernel.

The fragmented dense-consumer executor now also prebuilds a bounded symbolic
target map for the active producer run's trailing updates. Each mapped target
is classified once as a later external multiplier, current dense-panel entry,
pivot, current trailing-panel entry, or no-op before the numeric panel update is
applied. This removes the numeric apply loop's destination rediscovery while
keeping memory bounded by `batch_rows * max_run_trailing_len`, not by all runs
in the fragmented row. A focused `ASIC_320k` probe was mixed (`~0.326s`
refactor versus `~0.319s` before the target map), but the same five forced-row
CKTSO-gap probe improved slightly from about `32.1s` to about `31.9s`
geomean. This is useful as a storage/scatter scaffold, but it confirms that the
main remaining paper gap is still native row/segment panel storage and a
production blocked update kernel, not symbolic target lookup alone.

The same executor now keeps its fragmented batch metadata in worker-owned
index/byte scratch, and keeps the per-batch pivot vector in the existing double
workspace, instead of allocating row offsets, run offsets, target positions,
target kinds, and pivots for each small batch. This is a production-storage
cleanup rather than a numeric algorithm change. The focused `ASIC_320k` probe
stayed about `0.326s`, but the five forced-row CKTSO-gap geomean improved
slightly again from about `31.9s` to about `31.7s`. A separate temporary
rectangular run-panel plus portable blocked multiply experiment was tested and
rejected in the same session because it regressed the five-case geomean to
about `32.7s`; copying tiny batches into a denser temporary panel is not enough.
The remaining path needs broader row-panel batches and native row/segment
storage that avoids the copy rather than just a local blocked multiply.

KLS now retains the dense-producer target maps symbolically with the producer
run metadata instead of rebuilding them inside each fragmented dense-consumer
batch. Each producer trailing column is classified once as no-op, later
external multiplier, current dense-panel entry, current pivot, or current
trailing-panel entry; the numeric executor stores the accepted global run IDs
for the batch and reads those persistent classifications directly. A synthetic
smoke fixture covers the retained-map path with two dense producers, a dense
consumer block, and producer trailing columns that update the consumer trailing
panel. This also fixed a signed-size guard that had been discarding otherwise
valid retained maps on 64-bit platforms. On `ASIC_320k`, the retained map count
is now nonzero (`3,205,236` entries: `1,742,094` no-op, `316,396` external,
`305,279` dense, `4,477` pivot, and `836,990` trailing), and the focused
four-thread refactor probe improved to about `0.307s`. The same five forced-row
CKTSO-gap probe improved from about `31.7s` to about `30.3s` geomean, still far
from the no-opt-in baseline near `13.5s`. This closes the symbolic target-map
gap from the paper audit but reinforces that the large remaining gap is the
native row/segment panel representation and production blocked update executor,
not another small matcher or CPU-specific kernel tweak.

The same fragmented dense-consumer executor now participates in checked row fast
factorization. The prior guard skipped this batch path whenever pivot checks were
active, even though the CKTSO/SubtreeLU row-supernode update idea still applies
to the already-completed external producer runs before the consumer group's
internal pivot-checked panel step. KLS now performs the same row-wise multiplier
rejection checks before publishing external `L` multipliers from scalar gaps and
producer suffix solves, then lets the existing checked compact-panel factor code
handle internal multipliers and pivots. The retained target-map smoke fixture now
runs the two-producer fragmented dense-consumer shape through checked
factorization with CBLAS disabled and requires the fragmented batch counters,
retained producer target counts, and residual-clean solve. This closes a direct
checked-mode paper coverage gap while leaving the larger native row/segment
numeric engine and CKTSO pipelined pivoting tail open.

The exact-pattern and ragged dense-consumer batch executors now have the same
checked-mode treatment. These executors cover the cleaner SubtreeLU-style cases
where a dense consumer row range depends on completed dense producer suffixes
with either identical external patterns or one ragged producer suffix. They were
previously disabled whenever pivot checks were active, so checked fast
factorization fell back to scalar external updates even though the later
consumer-panel factor step already performs checked internal multiplier and pivot
tests. KLS now runs the exact/ragged producer-suffix batch under checked row fast
factorization and applies row-wise multiplier rejection before each external
`L` value is published. The existing exact and ragged smoke fixtures each run a
fresh checked factorization with CBLAS disabled and require the compact
supernode batch counters plus a residual-clean solve. This fills another direct
paper-coverage hole without claiming to solve the remaining production
row/segment panel engine or CKTSO pivoting-tail executor.

The independent `GROUP_BATCH` row-supernode executors now have the same checked
coverage. The ragged single-producer suffix executor, the contiguous
multi-producer dense-suffix executor, and the fragmented scalar-gap plus
multi-producer executor no longer skip checked row fast factorization. Each
path rejects unsafe `L` multipliers before publication, tests the completed
row pivot against the row's pivot/U maximum before publishing `Udiag`, and only
then marks the independent row done. The smoke suite now runs checked
factorizations for ragged, contiguous multi-producer, and fragmented independent
row batches with CBLAS disabled and requires the compact-supernode batch and
dense-producer counters. This closes the checked-mode dispatch gap for the
paper row-supernode batch shapes while leaving the larger native row/segment
numeric engine and CKTSO pivoting-tail executor open.

An exact-match requirement was deliberately kept for fragmented dense-producer
batches. A common-prefix widening experiment was tried and rejected because it
increased the work on `ASIC_320k` (`~0.387s` focused refactor versus the prior
`~0.329s`) by admitting extra scalar suffix work into the batch. That result is
consistent with the paper lesson: broadened dispatch is not enough when the
underlying row-panel storage and update kernel are still the limiting pieces.

Three further local row-panel shortcuts were tested after retained target maps
and rejected. Streaming producer trailing updates through a one-row scratch
vector avoided the `batch_rows * trailing_len` temporary, but changed the access
pattern without reducing arithmetic: the five forced-row CKTSO-gap geomean
regressed from about `30.27s` to `30.55s`. Computing only retained-map targets
whose class was not no-op also regressed (`32.59s` geomean), and the aggregate
GEMV counters stayed essentially unchanged on the executed batches, showing
that the global no-op target count was not the active-batch bottleneck. A
portable batched TRSV over equal producer suffixes passed smoke but regressed
the focused `ASIC_320k` refactor from about `0.307s` to `0.314s`, because the
fragmented batches are too small and ragged for that loop interchange to
amortize its overhead. A temporary retained dense-input scatter-map experiment
for compact dense rows also regressed (`32.82s` five-case geomean), so repeated
input destination lookup alone was not the large missing paper mechanism. KLS
later kept the cleaned retained-target compact loader for storage-path
completeness and coverage, not as a claimed performance lever. These rejected
variants narrow the next useful implementation target: KLS needs a native
row/segment numeric representation that stores and updates producer/consumer
panels in the execution order directly, rather than more symbolic shortcuts
around the current KLU-shaped row mirrors.

The unchecked producer-panel refactor experiment also uses a blocked panel
algorithm: scalar code factors each diagonal block, `dtrsm` solves the
below-panel multiplier block, and `dgemm` updates both the dense right panel
and the shared trailing panel.
The scalar compact kernel remains the default because same-session `onetone2`
forced-row probes still favored it: the blocked CBLAS panel path was about
`0.072s` versus `0.029s` with the scalar fallback, with both runs
residual-clean. This improved on the earlier per-row CBLAS attempt (`0.219s`)
but confirms the paper gap more precisely: KLS needs batched
row-group/supernode consumer updates and a deeper blocked row-major numeric
layout, not BLAS calls wrapped around each current row in the existing compact
group shape.

KLS then filled a narrower but direct CKTSO tail-restart semantic gap for
prefix-current masked tails. A prefix-current root reject normally means every
later unstarted column belongs to the unfinished seed set, so the safe tail
remains a full suffix. When the reject-only ETree closure is shorter, however,
KLS can now
prove that preserved columns outside that closure have no U dependency on tail
columns, refresh those preserved columns with the existing mapped no-pivot
column kernel, and then execute the shorter pivoted tail envelope. The refresh
proof now works for single-block, BTF, unscaled, and KLU row-scaled diagonal
blocks by using the mapped global-column dispatcher, requiring retained
block/off-diagonal metadata before refreshing omitted block-local columns, and
feeding scaled attempts with input-row `Rs`. The smoke fixtures cover both a
weak root in a 2-by-2 dependent part plus an independent trailing singleton,
and a non-root weak pivot whose independent changed suffix column must be
refreshed before KLS preserves it outside the ETree-derived tail. Both assert
that the retained pivoting-tail plan is shorter than the block suffix and that
KLS counts the repair as `fast_tail_restarts=1`.
This is still a serial envelope rather than CKTSO's parallel Algorithm 5
executor, but it directly applies the paper's distinction between unfinished
EGraph nodes and ETree-descendant pivoting-tail work instead of treating all
safe prefix-current rejects as whole-block repairs.

The first-factor path was then moved closer to the paper instead of leaving it
as an opt-in experiment. Unset `KLS_ENABLE_KLS_FIRST_FACTOR` now means
conservative auto: large eligible symbolic blocks try the KLS-owned sparse
row-up-looking first factor before falling back to KLU, while `0` remains a hard
off switch and `1` still forces the attempt. More importantly for the
static-pivoting rows called out by the gap decomposition, an accepted pre-static
row-matching candidate can now replay its selected symbolic/value state through
the KLS row-up factor after the existing KLU trial has accepted the candidate.
If the replay fails, KLS restores the accepted trial numeric; if it succeeds,
`initial_factor_path` reports `kls_first` and row-major mirrors are prepared as
for an ordinary KLS-first factor. This does not remove the KLU trial used to
score static-pivot candidates, and it does not implement CKTSO's parallel
row-up/EGraph first-factor executor, but it removes a direct KLU-storage return
from accepted first-factor states and makes the row-up engine a production
candidate for the large cases the papers target.

The row-up first factor now also keeps its row-oriented product alive for the
next phase. While packing the accepted factors into the KLU-compatible numeric
object for fallback solves and pivot repairs, KLS records the same row entries,
factor-order input positions, and numeric value pointers into the shared
row-refactor CSR/group metadata finisher. A successful direct handoff is
reported by `kls_first_last_row_refactor_seeded_rows` and
`kls_first_row_refactor_seeded_row_count`; smoke coverage requires it for the
ordinary row-up, dynamic-column-pivot, BTF, large-auto, and pre-static replay
first-factor cases. This still is not CKTSO's fully row-major primary numeric
object, because BTF off-diagonal refresh and fallback coherence still require
the retained column-oriented refactor map, but it removes the previous
pack-then-reconstruct step for KLS-first row-refactor mirrors.

The pivot-repair path now attempts the retained non-contiguous CKTSO-style
ETree tail mask even when the pivoting tail reaches the end of the BTF block.
Columns outside the unfinished ETree closure are preserved only when the
existing structural dependency proof can lock their old pivots; otherwise KLS
falls back to the full suffix restart. This fills a direct semantic gap in the
current tail executor without pretending to implement CKTSO's full parallel
pipelined tail factorization.

The compact dense row-segment refactor path then removed another KLU-storage
adaptation step. For dense groups with retained compact panel slices, KLS now
loads current input values for in-panel dense columns and shared trailing
columns directly into that compact panel, leaving the sparse work vector only
for external dependency columns and their update deltas. Checked and unchecked
compact dense smoke cases require this direct-input path, and benchmark JSON
reports the last-run and cumulative direct compact-input row counts. This is a
small but direct move toward the SubtreeLU/CKTSO row/segment numeric-storage
model: dense panel values live in the row-segment panel for the numeric pass
instead of being staged through a KLU-shaped sparse accumulator first.

The same direct-input idea now covers native dense row segments that do not use
the retained compact-panel path. For dense groups below the compact arithmetic
intensity gate, KLS writes current in-group lower, diagonal, upper, and shared
trailing input values directly into the row-segment `L`/`U` mirrors and keeps
the sparse work vector for external dependencies and their update deltas. A
46-by-46 dense smoke fixture is deliberately large enough to form a dense
segment and deliberately below the compact-panel gate; it now requires native
dense direct-input rows and zero compact-panel executions. This removes another
KLU-shaped staging step from the default row-segment kernel while preserving
the conservative compact-panel gate.

The retained SubtreeLU separator-tree queue is now consumed by the checked
fast-factor row path as well as the unchecked no-pivot refactor path. Previously
the Algorithm 6-style FLOP-balanced private/pipeline queue was explicitly
disabled when pivot checks were active, leaving checked fast factorization on
the older generic ready queue even when a METIS `NodeNDP` separator map covered
the block. The checked row scheduler now keeps the same separator private
groups and pipeline groups while relying on the existing earliest-reject stop,
done bitmap, and prefix-refresh validation before any CKTSO tail repair is
accepted. A generated 30,600-row sparse-grid smoke fixture uses METIS,
no BTF, no scaling, and checked row fast factorization, and requires
`row_refactor_last_separator_flop_queue=1` with retained separator components,
private groups, pipeline groups, and a clean solve residual. This directly
fills the paper gap of applying SubtreeLU's separator private/pipeline
partition inside the pivot-aware fast path, although it is still not CKTSO's
full ETree-descendant pivoting-tail executor.

A current CKTSO-gap forced-row probe showed the closure rule enabling the
separator FLOP queue on `ASIC_320k`, `ASIC_320ks`, `ASIC_100ks`, and `rajat28`,
where earlier builds rejected the schedule because pipeline groups released
private successors. The effect was modestly positive for the ASIC cases and
neutral to slightly negative for `rajat28`, so this is best understood as
closing a paper-algorithm coverage gap rather than as the remaining large
performance lever by itself.

The separator map itself is now BTF-aware instead of depending on a unique
block-size match after KLU analysis. KLU invokes the METIS user-order callback
once per non-tiny BTF diagonal block, so KLS records each retained `NodeNDP`
component tree with that callback ordinal. After symbolic analysis returns the
accepted `R` block boundaries, KLS stitches the captured local separator trees
into one global postorder forest and synthesizes private components for blocks
that did not run `NodeNDP`. Artificial pipeline nodes connect the BTF roots so
the existing SubtreeLU Algorithm 6-style queue builder can consume one global
component map. A duplicate-size BTF smoke fixture builds two independent
30,600-row sparse-grid blocks; it requires the separator range to cover all
61,200 rows, reports at least 15 global separator components, and runs checked
row fast factorization through the separator FLOP queue with a clean solve.
This closes the earlier duplicate-block ambiguity and lets retained separator
metadata survive BTF forests, while still leaving production supernodal BLAS and
CKTSO's pivoting-tail executor as larger remaining paper gaps.

The KLS-first automatic selector now refuses one class where the implementation
was ahead of the paper coverage: very large scaled matrices whose accepted
symbolic state is one BTF block. The CKTSO paper treats this class with a
parallel row-up/ETree task factorization. KLS therefore keeps the automatic
skip for this class unless the accepted analysis retained a global separator
private/pipeline row queue that the KLS-first row-up executor can consume.
On the CKTSO-gap focus subset, `rajat24` was
the clear failure mode: before this guard it selected `kls_first`, spent about
40.1 s in initial factorization, and modeled at 66.1 s versus CKTSO's 4.68 s.
After the guard, default automatic mode reports
`kls_first_auto_skipped_scaled_single_block_count=1`, keeps the KLU/static first
factor, and models at 10.8 s. The 12-row focus geomean moved from 3.78 s
to 3.34 s, reducing the CKTSO ratio from 2.96x to 2.62x. This is not the final
paper algorithm; it prevents unpartitioned scaled single-block systems from
masking the actual missing piece, which is CKTSO's parallel row-up first factor
and the matching row-oriented refactor/solve engine.

The partial-supernode pipeline split is now a structural automatic row-engine
choice instead of an explicit experiment. With
`KLS_ENABLE_PARTIAL_SUPERNODE_PIPELINE` unset, row refactorization uses the
SubtreeLU/CKTSO split when the pipeline tail is dominated by large producer
groups; setting the variable to `0` remains a hard disable. The same
row-prefix publication is also allowed in checked row fast-factor/refactor
runs: a producer row is marked done only after its row values and pivot check
complete, and consumers still wait on each dependency row before applying the
partial compact-supernode update. This fills a direct Algorithm 5 coverage gap
for prefix-safe producer/consumer overlap, while the full CKTSO pivoting-tail
executor and production row-major first factor remain open.

The column EGraph refactor schedule now retains exact consecutive
supernode-candidate ranges instead of only counting them. Benchmark JSON and
gap decomposition output report the number of retained candidates and, when
enabled, the number of natural-order EGraph pipeline tasks and columns that
were coarsened. A guarded execution path exists behind
`KLS_ENABLE_EGRAPH_SUPERNODE_TASKS=1`: it only leases candidates whose width is
within the CKTSO-style `2 * threads` task-width bound, and still dispatches and
marks every column through the existing dependency and pivot checks. This is
off by default because the direct experiment showed that scheduler coarsening
without the matching supernodal numeric panel/update kernel is not the missing
CKTSO lever. On `G2_circuit` with four threads, default EGraph refactor stayed
near 0.201 s while the env-enabled coarsened path leased 6,991 tasks covering
23,328 columns and slowed to about 0.209 s. The 11 completed rows of the
CKTSO-gap focus subset were 1.0067x slower than the previous KLS run when the
task path was enabled; the `transient` timeout observed during that run was
also reproduced on the committed `d08ab35` baseline and is a pre-existing
EGraph pipeline flake. This narrows the remaining paper gap: the useful next
step is not more task scheduling, but production supernodal numeric storage
and BLAS-style panel/trailing updates for these retained ranges.

The intermittent 120 s `transient` timeout was traced to the clustered
EGraph/row-refactor barrier protocol, not to the supernode-range metadata. A
worker could observe `stop` after one level barrier and leave the clustered
phase while another worker had already entered the next level's barrier,
leaving the remaining workers asleep in `pthread_barrier_wait` and the main
thread waiting on the pool completion condition. Clustered column and row
refactor workers now drain all remaining level barriers collectively after a
stop; once stopped, they skip numeric work but still rendezvous with peers
until the level loop is complete. This preserves the existing fast
fetch-and-wait EGraph pipeline while removing the mismatched-barrier timeout
class. A 20-run `transient` stress loop with four threads completed without
timeout in default mode, with ready-queue columns reported as zero and refactor
time staying near the previous 0.026-0.030 s range.

A true column ready-queue scheduler is also available behind
`KLS_ENABLE_EGRAPH_READY_QUEUE=1`. It builds in-tail predecessor counts and
successor lists from the retained EGraph schedule and only dispatches columns
whose pipeline predecessors have completed. On `transient` it exercised all
1,037 pipeline columns and completed reliably, but repeated refactor time rose
to about 0.066-0.072 s, so it remains an off-by-default paper probe rather
than the production path. The ready-queue successor graph is now also built
only when that probe is enabled; a gated default 12-row CKTSO-gap focus run
completed with zero failures and a 3.216 s geomean, effectively unchanged from
the saved 3.206 s KLS reference while still about 2.52x slower than the saved
CKTSO reference on the same common rows. This confirms that the large CKTSO gap
is not closed by stricter task readiness alone; the next algorithmic gap is
still the numeric supernodal/panel update engine.

KLS now has a first column-EGraph numeric supernode update probe behind
`KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=1`. It reuses the retained consecutive
supernode candidates, detects a contiguous dependency run, validates the dense
internal L block, solves that dependency vector through a worker-local dense
triangular panel, and only commits to the grouped path when the producer L
columns also share one trailing row list so the update can be accumulated in a
single worker-local vector before one scatter. Benchmark JSON reports
`refactor_last_supernode_update_runs`, rows, entries, and cumulative
`refactor_supernode_update_*` totals. This fills the direct SubtreeLU/CKTSO
TRSV-plus-matrix-vector update semantics in the mapped EGraph refactor kernels,
but it is still off by default because it rebuilds the dense/trailing panel for
each consumer instead of publishing a persistent producer panel.

The same-session six-row CKTSO-gap focus check makes that missing storage piece
clear. Default KLS on the first six focus rows
(`build/kls_default_same_focus6_t4_r3_timeout120.jsonl`) had an 8.24 s
SPICE-cycle geomean. The common-trailing supernode update probe
(`build/kls_egraph_supernode_updates_common_focus6_t4_r3_timeout120.jsonl`)
completed the same rows but rose to 20.36 s. It did exercise real common
trailing panels: `G2_circuit` reported 424,395 grouped updates, 5,541,687 rows,
and 1.51e9 update entries; `ASIC_100ks` reported 160,646 updates and 3.06e8
entries. The large gap therefore is not a missing task queue or lack of
supernode detection anymore; it is the paper's production compact supernode
storage/publish step, so KLS can build a producer panel once and let many
consumers reuse it instead of reconstructing the panel at every dependency run.

KLS now has that persistent producer-panel step for the column EGraph probe.
When `KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=1` is set, the mapped EGraph refactor
builds a compact panel cache for retained supernode candidates whose L columns
have a dense internal block and one common trailing row list. Producer columns
publish their finalized L values into the cache before the existing
`pipeline_done` release, so downstream consumers reuse the panel after the
normal dependency wait. On the same six-row CKTSO-gap focus subset, the cached
path (`build/kls_egraph_supernode_panel_cache_focus6_t4_r3_timeout120.jsonl`)
improved the opt-in probe from 20.36 s to 10.43 s geomean, but remained slower
than the same-session default KLS geomean of 8.29 s
(`build/kls_default_panel_cache_base_focus6_t4_r3_timeout120.jsonl`). The
counters show the same high-volume updates are now covered without per-consumer
panel reconstruction: `G2_circuit` still applied 424,395 grouped updates and
1.51e9 entries in the last refactor. This fills the direct storage/publish
piece from the papers but does not close the slow-case gap by itself. The
remaining direct paper gap is compact or batched numerical kernels over these
panels, not more supernode detection.

The column EGraph cache now also consumes published supernode prefixes and
suffixes, matching the partial-publication idea already present in the
row-refactor compact panel path. The cached consumer uses `col_id` rather than
start-only lookup, solves the available dependency subrun inside the retained
dense panel, scatters any in-panel rows beyond the published prefix into the
current column work vector, and accumulates the shared trailing panel once. This
removes another direct paper gap: current columns inside a retained supernode
can reuse already published producer columns instead of rebuilding a partial
panel. On the same six-row focus subset, the corrected prefix/suffix path
(`build/kls_egraph_supernode_panel_prefix_correct_focus6_t4_r3_timeout120.jsonl`)
improved the persistent-panel probe from 10.43 s to 9.44 s geomean. The default
KLS column EGraph path remains faster at 8.29 s geomean, so this is still kept
behind `KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=1`. The next gap remains a
coarser compact/batched numeric update kernel; per-consumer scalar loops over
the retained panel are not enough.

The generic mapped EGraph kernel now consumes the same cached compact panels
before falling back to its scalar dependency loop. This closes a coverage gap
left by the earlier single-block and large unscaled-BTF specializations: scaled
BTF and smaller BTF refactors can use retained producer panels when
`KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=1` is set. A targeted scaled-BTF
`ckt11752_dc_1` probe reported 313 last-run compact supernode dependency
updates over 6,555 rows and 415,380 entries with the opt-in flag, while the
same command without the flag reported zero such updates. This is a direct
SubtreeLU/CKTSO compact-panel coverage step, but it still does not implement
the coarser batched numeric kernel needed to make the probe a default win.

The cached EGraph-panel consumer now has the same optional CBLAS shape as the
row-supernode experiment for sufficiently large retained panels. In CBLAS
builds with both `KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=1` and
`KLS_ENABLE_CBLAS_SUPERNODE=1`, an eligible cached dependency run solves the
unit-diagonal internal panel with `dtrsv`, applies any dense in-panel suffix
with `dgemv`, and applies the shared trailing rows with one more `dgemv`
before the existing scatter. Benchmark artifacts report that BLAS-taken subset
through `refactor_last_supernode_cblas_update_*` and cumulative
`refactor_supernode_cblas_update_*` counters. This is still deliberately behind
the opt-in supernode-update probe; it makes the compact panel mathematically
closer to SubtreeLU's BLAS update shape, but does not replace the remaining
need for broader batched producer/consumer kernels or CKTSO's pivoting tail
executor.

The KLS-owned row-up first factor now executes independent BTF diagonal blocks
in parallel when KLS-first factorization is requested with more than one
thread. Each worker uses private row-up scratch, performs dynamic column pivots
inside its assigned block, and then commits the block's KLU-compatible numeric
storage and row-refactor seed data under a short mutex. This closes the direct
paper gap where the KLS-first bridge was purely serial across BTF blocks even
though CKTSO/SubtreeLU treat first factorization as a row/task parallel
problem. The smoke suite now forces a two-block BTF row-up case with two
threads and requires `kls_first_last_parallel_btf_blocks=2`. This is still not
the full CKTSO ETree pipeline inside one large diagonal block; it is the safe
BTF-level parallel part of that first-factor algorithm.

The automatic scaled single-block guard now recognizes the separator row-up
executor added after the original guard. Very large scaled single-block states
are still skipped in automatic mode when they lack retained global separator
private/pipeline work, but a METIS `NodeNDP` analysis that preserves both
private and pipeline separator rows may now enter KLS-first automatically. This
fills a direct coverage gap against the papers' row-up task-factorization
precondition without enabling the old single-block bridge broadly. The smoke
suite keeps both sides covered: a natural-order 150,000-row scaled diagonal
remains on `klu_first` with `kls_first_auto_skipped_scaled_single_block=1`,
while a 150,000-row scaled METIS/no-BTF tridiagonal uses `kls_first` and
requires the separator private and pipeline phases to execute.

The row-refactor generic segment path now direct-loads raw input rows into
retained sparse row-major segment storage when the existing row pattern can
represent the row exactly: external dependencies remain in the work vector,
internal sparse `L` entries are placed in their row mirror slots, the pivot is
placed in `Udiag`, and in-segment/trailing `U` entries are placed in the row
mirror before numeric updates. Checked runs compute the CKTSO-style row pivot
test from the direct-loaded U entries plus residual updates, so this is a
storage/algorithm bridge rather than a pivot-policy shortcut. That path now
also keeps a retained per-input-entry target map for supported dense and sparse
segment rows. The symbolic row-pattern builder classifies each raw input entry
once as an external work-vector value, row-major `L` slot, pivot, or
in-segment/shared trailing `U` slot; repeated refactors then direct-load
through those retained destinations instead of searching the segment structures
inside the numeric loop. The dense direct-input smoke case now requires target
loading over the dense segment, and the generic sparse smoke case uses a full
upper shared U segment with only sparse subdiagonal L dependencies and requires
both `row_refactor_last_sparse_segment_direct_input_rows` and
`row_refactor_last_segment_target_input_rows` to cover the whole segment. This
narrows the SubtreeLU/CKTSO row-major storage gap for generic segments, but it
still leaves the larger paper work item open: a production compact/batched
row-major numeric engine and CKTSO's pipelined pivoting-tail executor.

The retained segment plan now also owns residual cleanup destinations for
supported target rows that still need post-row input cleanup. During symbolic
row-pattern construction, KLS records the external work-vector positions from
the retained per-input-entry target map in a compact row pointer/list. Numeric
row cleanup consumes that retained list first and only falls back to the old raw
input-column scan when no retained cleanup list exists. The diagnostics now
separate planned cleanup coverage
(`row_refactor_segment_input_cleanup_rows/entries`) from cleanup actually
performed through retained targets
(`row_refactor_last_segment_target_cleanup_rows/entries` and cumulative
companions). This removes another raw KLU/CSC-style rediscovery step from the
row-segment loop, but it remains scaffolding for the paper algorithm rather
than the full production row-major engine.

The independent row-batch executors now also direct-load refactor input rows
into their native batch panels. Ragged, multi-producer, and fragmented
multi-producer batches validate that each raw input entry maps exactly to the
batch's dependency vector, pivot, or row-major U workspace, then populate those
buffers without first staging through the sparse work vector. Rows that do not
fit the current batch shape still fall back to the old residual-safe path. The
smoke probes for all three independent batch shapes now require
`row_refactor_last_batch_direct_input_rows` to cover the executed batch rows.
This fills a direct storage-path gap versus the papers' row-major numeric
updates; it is not yet a replacement for CKTSO's full pipelined pivot-tail
factorization.

The retained SubtreeLU Algorithm 6 queue splitter then gained the missing
private-leaf guard. A focused checked METIS grid exposed that the earlier
FLOP-balance loop could keep promoting separator candidates after it reached
indivisible retained components, leaving no private subtree candidates and
falling back to the generic row-DAG queue. The splitter now only promotes a
candidate component into the pipeline when it has at least one positive-work
child subtree to return to the private candidate set. The checked separator
grid and duplicate-BTF separator forest smoke fixtures both select the
separator FLOP queue again, preserving the paper's private/pipeline structure
instead of turning unbalanced leaves into separator work.

KLS now fills the next direct CKTSO Algorithm 5 gap inside the checked
row-refactor pipeline. When a row blocks on an unfinished predecessor, the
worker scans later row dependencies and consumes any already-finished
predecessor whose value is provably final: every skipped earlier dependency must
lack a row-major `U` entry into that later dependency. Ambiguous pattern data,
missing row-major mirrors, and checked multipliers that would reject are left
for the normal in-order path, preserving the existing pivot-reject order. This
implements the paper's "use newly detected finished predecessors while waiting"
idea for the KLS-owned row-major executor and reports actual use through
`row_refactor_last_prefactor`, row/dependency counts, and cumulative
`row_refactor_prefactor_*` counters. It still does not claim CKTSO's complete
ETree-descendant pivoting-tail factorization: KLS has the guarded row-level
prefactor step, not the full tail scheduler that retopologizes all unfinished
descendants after a pivot exchange.

The guarded prefactor step now also consumes finished dense-producer runs as
runs, not only as scalar dependencies. If later dependencies after the blocked
predecessor are consecutive rows from a retained dense producer, the producer
panel prefix is already published, each dependency is finished now, and skipped
earlier dependencies have no `U` edge into the candidate run, the worker solves
the run's internal triangular prefix in workspace, checks all multipliers
without recording an out-of-order rejection, and then applies the dense suffix
and shared trailing update in one compact supernode pass. Any multiplier that
would fail the checked threshold cancels the compact prefactor attempt and
falls back to the existing in-order scalar path. The smoke fixture widens the
finished producer behind a blocked predecessor and now requires
`row_refactor_last_prefactor_supernode=1`, while benchmark output reports the
last and cumulative `row_refactor_prefactor_supernode_*` counters. This maps
the paper's Algorithm 5 "use supernode k:k'" waiting-phase operation more
directly, but it remains guarded row-level prefactoring rather than CKTSO's full
ETree-descendant pivoting-tail restart.

The KLS-first separator-pipeline factorization now has the analogous grouped
executor for SubtreeLU Algorithm 4's ready-supernode branch. Earlier code
identified published row-major `U` supernode runs and counted them, but consumed
the rows by repeatedly calling the scalar dependency updater. The new
`kls_row_first_partial_apply_supernode_run` helper owns the run: it applies the
first ready dependency, lets that triangular prefix update discover the next
in-run row if needed, pops each newly ready row in order, and reports a
supernode update only when more than one row is actually consumed. It then moved
one step closer to the paper's `trsv`/`gemv` shape: when the published rows
prove a dense upper-triangular prefix plus one common trailing row list, KLS
performs the in-supernode solve over the current row and accumulates the common
trailing update in reusable worker scratch before a single scatter over the
tail. If that compact shape is not validated, the helper falls back to the old
scalar row-by-row run. This fills a numeric-kernel gap in the KLS-owned
first-factor pipeline without pretending that the full paper supernodal storage
layer is finished; the values still flow through KLS row-entry storage rather
than a persistent BLAS panel.

That KLS-first ready-supernode selector now also releases the longest
already-finished prefix of a retained row-supernode instead of requiring the
whole supernode to be below the current ready frontier. The compact helper
treats later rows of the same producer supernode as part of the common tail,
so they remain pending dependencies in the current row while the ready prefix
still contributes immediately. This maps the paper's partial producer-release
idea into the KLS-first pipeline without adding a new scheduling threshold or a
matrix-specific rule.

The same KLS-first compact/prefix supernode executor is now used by ordinary
row-up dependency loops: natural/serial rows, separator private rows, serial
pivot-tail rows, and the BTF-level parallel first-factor worker. A new
conservative run-bound check extends a dependency only through adjacent
row-major `U` rows that are already owned by the same private worker, already
published in the shared ordered suffix, or already finished inside the current
BTF worker block, and that validate as the same supernode pattern. The row then
consumes the run through the shared triangular-prefix/trailing-update helper
instead of re-entering the scalar dependency loop for each producer row. This
closes the earlier mismatch where SubtreeLU-style supernode updates existed in
the separator pipeline but not in the private-mode or BTF-parallel row factor
paths. Stats now expose actual use through
`kls_first_last_row_supernode_update`,
`kls_first_last_row_supernode_update_groups`, and
`kls_first_last_row_supernode_update_rows`, separate from the older
separator-pipeline-only counter.

KLS now fills more of that persistent-panel storage gap. KLS-first row-up
producers publish a completed row-supernode as a dense upper-triangular panel
plus common trailing column list when the next row proves that the supernode has
ended. This covers parallel BTF workers plus the ordinary private and serial
row-up loops. Later rows try that cache before falling back to row-entry
validation, and benchmark output reports actual use through
`kls_first_last_row_supernode_panel_update`,
`kls_first_row_supernode_panel_update_run_count`,
`kls_first_last_row_supernode_panel_update_groups`, and
`kls_first_last_row_supernode_panel_update_rows`. At the start of each ordered
KLS-first separator pipeline phase, the
stable private-prefix `U` snapshot is scanned for validated row-supernodes and
published as phase-local dense upper-triangular panels plus a common trailing
column list. Pipeline workers try that cache before revalidating row entries,
so repeated consumers no longer rebuild the same private-prefix panel. A
dynamic column exchange rebuilds that phase cache from the post-exchange
private-prefix `U` snapshot; row-up producer panel caches are reset because
their open producer tail layout is tied to the old column order. Benchmark
output reports actual panel-backed use through
`kls_first_last_separator_queue_pipeline_supernode_panel_update`,
`kls_first_separator_queue_pipeline_supernode_panel_update_run_count`,
`kls_first_last_separator_queue_pipeline_supernode_panel_update_groups`, and
`kls_first_last_separator_queue_pipeline_supernode_panel_update_rows`. In
CBLAS builds, the existing `KLS_ENABLE_CBLAS_SUPERNODE=1` runtime gate can
consume eligible cached KLS-first panels with CBLAS `dtrsv` over the
non-unit upper-triangular panel and CBLAS `dgemv` for dense suffix/common-tail
updates; otherwise the cached panel uses the scalar in-panel solver. The cached
consumer no longer requires the whole requested dependency run to fit inside the
same panel: it consumes the published prefix and leaves the remaining suffix in
the dependency heap, matching the papers' private/pipeline rule that completed
producer prefixes can be used before the producer tail is available. The
ordinary row-up path now also publishes a compact-validated producer prefix
lazily when a consumer first uses it, so a later consumer can reuse the same
dense/common-tail panel even before the complete producer supernode is known.
The same lazy prefix publication now applies inside separator private/pipeline
row-up paths after the compact prefix has been validated under the scoped pivot
order; dynamic pivots reset or rebuild the affected caches from the committed
post-exchange prefix so stale column layouts are not reused. This is still not
the full paper storage layer: KLS does not proactively maintain mutable
open-supernode panels, and CKTSO's full ETree-descendant pivoting-tail
scheduler remains open. Separator pipeline
pivot-tail rows that
are serialized after a restart now use the same row-up producer panel cache for
completed-supernode publication, so the restarted suffix no longer loses those
completed panels just because a pipeline phase fell back to a serial pivot row.

The BTF-level parallel first-factor worker no longer carries its own older
row-up numeric loop. It now calls the shared `kls_row_first_factor_one_row`
executor used by the serial, private, separator, and pivot-tail row paths, then
performs only the worker-specific row pointer continuation, panel publication,
separator reject accounting, and block commit. This removes a remaining
algorithmic fork: dynamic column exchanges, scoped separator pivot decisions,
owned-supernode run detection, and cached panel consumption now have one
implementation across the KLS-first row-major executors. The paper gap that
remains is the larger one already noted above: KLS still needs the full CKTSO
ETree-descendant pivoting-tail scheduler and a production compact/batched
row-major numeric storage layer.

The retained row-refactor input-target map now covers the scalar/generic row
path instead of only dense and sparse segment kernels. Pattern construction
marks single-row groups whose raw input entries can be represented by the row's
external dependency work-vector slots, pivot, and retained row-major `U` tail.
The scalar row executor consumes that map by loading pivots and `U` entries
directly into row storage, leaving only true predecessor dependencies in the
work vector, and then combines those direct values with dependency updates
before the row-wise pivot check and publish. A new smoke fixture uses a
tridiagonal chain where all rows have retained targets but only a small subset
uses the sparse-segment direct-input kernel, proving the generic row path uses
the retained destinations too. This is another storage-path step toward the
CKTSO/SubtreeLU row-major numeric engine; it still does not replace the larger
compact/batched engine or CKTSO's ETree-descendant pivoting-tail scheduler.

The compact dense row-refactor panel loader now consumes the same retained
per-input-entry target map. When the symbolic row target data is available, raw
input values are placed directly into the compact row's external work-vector
slots, dense `L` panel, pivot, dense in-panel `U`, or trailing `U` panel instead
of rescanning trailing columns during numeric refactorization. The compact-panel
smoke coverage now requires both compact direct-input rows and retained target
consumption, so this closes another storage-path mismatch against the
CKTSO/SubtreeLU row-major algorithm description. The larger open paper gap is
unchanged: KLS still needs a production compact/batched row-major numeric engine
and CKTSO's pipelined ETree-descendant pivoting-tail executor.

The checked-reject KLS-owned block repair now attempts the restartable
Algorithm 5-style row pipeline first when multiple threads are requested, before
the older serial pivot-tail repair and KLU fallback. This path still works at
BTF-block repair scope rather than CKTSO's exact ETree-descendant tail mask, but
it removes a direct mismatch in the fast-reject recovery flow: a failed checked
factor can continue through the row-up pipeline and, if needed, serialize only
dynamic-pivot epochs inside that pipeline. Benchmark and
decomposition output expose the behavior through
`fast_kls_block_restart_last_row_pipeline`,
`fast_kls_block_restart_row_pipeline_count`,
`fast_kls_block_restart_last_row_pipeline_rows`,
`fast_kls_block_restart_last_row_pipeline_threads`,
`fast_kls_block_restart_last_row_pipeline_prefix_rows`,
`fast_kls_block_restart_last_row_pipeline_suffix_rows`,
`fast_kls_block_restart_last_row_pipeline_gap_rows`,
`fast_kls_block_restart_last_row_pipeline_pivot_tail_rows`,
`fast_kls_block_restart_last_row_pipeline_pivot_restarts`,
`fast_kls_block_restart_last_row_pipeline_supernode_update_groups`,
`fast_kls_block_restart_last_row_pipeline_supernode_update_rows`,
`fast_kls_block_restart_last_row_pipeline_supernode_panel_update_groups`, and
`fast_kls_block_restart_last_row_pipeline_supernode_panel_update_rows`; the
prior-pivot smoke fixture now requires the two-thread repair to use that
pipeline.

KLS now narrows that block-repair pipeline when the retained pivoting-tail plan
is an exact topological suffix ending at the end of the rejected BTF block. In
that case, the repair seeds the preserved prefix from the previous LU into
row-major entries, marks those prefix rows finished, and starts the restartable
row pipeline at the rejected suffix rather than at row zero. The activation is
deliberately conservative: if the prefix has non-identity row order or the
retained plan is a non-contiguous/shorter ETree mask, KLS leaves the repair on
the existing serial exact-mask path where preserved later rows are explicitly
validated. The BTF suffix smoke case now requires a nonzero preserved-prefix
count, showing that the multi-thread repair is no longer only a full-block
pipeline in exact-suffix cases.

The same row-pipeline repair now handles the next conservative CKTSO tail case:
a contiguous, topological pivoting-tail envelope that stops before the end of
the rejected BTF block. KLS builds the retained tail mask, refreshes preserved
non-tail block columns just like the serial exact-mask path, seeds both the
preserved prefix and suffix into row-major storage, and runs the restartable
pipeline only on the active envelope. During that envelope, dynamic pivot
selection is bounded to active tail columns; if a row would need a preserved
suffix column to satisfy the pivot check, the pipeline attempt fails and the
existing serial/KLU fallback ladder remains responsible. The root-independent
tail smoke fixture now runs with two threads and requires two active pipeline
rows plus one preserved suffix row for both unscaled and scaled repairs.

KLS now fills the next direct CKTSO tail gap by admitting non-contiguous
topological pivoting-tail masks into the same restartable row pipeline. The
pipeline now consumes the retained ETree-descendant tail worklist as its
compact topological row order, keeps preserved gap and suffix rows seeded from
the prior LU, and applies readiness through an active-rank map instead of raw
row-number thresholds. Suffix and contiguous retained tails also enter through
that worklist path instead of the earlier broad interval shortcut, so successful
KLS-owned block repairs now execute CKTSO Algorithm 5's retained restart-node
sequence inside the row-first pipeline before using the existing fallback ladder
for unsafe pivot exchanges.
Public and benchmark statistics now report
`fast_kls_block_restart_last_row_pipeline_gap_rows`,
`fast_kls_block_restart_last_row_pipeline_etree_tail`,
`fast_kls_block_restart_row_pipeline_etree_tail_count`,
`fast_kls_block_restart_last_row_pipeline_etree_tail_rows`,
`fast_kls_block_restart_last_row_pipeline_etree_tail_gap_rows`, and
`fast_kls_block_restart_last_row_pipeline_etree_tail_exact_mask`; the
non-contiguous gap smoke fixture requires a preserved gap row and no serial
tail restart when two threads are enabled.

The retained-tail row-pipeline repair now also validates SubtreeLU's component
pivot domain before accepting the active ETree mask. If separator analysis
covers the repaired BTF block, every retained active row must map to a valid
collapsed-separator component and dynamic pivot choice remains bounded by that
component extent intersected with the active tail mask; otherwise KLS falls back
to the existing serial/KLU repair ladder. Benchmark output records this direct
paper-algorithm overlap with
`fast_kls_block_restart_last_row_pipeline_separator_tail_scope`,
`fast_kls_block_restart_row_pipeline_separator_tail_scope_count`, and
`fast_kls_block_restart_last_row_pipeline_separator_tail_scope_rows`.

The serial fallback now consumes the same retained topological tail envelope
instead of forcing the failed pivot to be the restart boundary. When the
retained ETree-descendant plan is non-suffix, KLS builds the exact active mask
over that plan's first/last columns, reconstructs the current prefix only before
the retained tail begin, and records/refreshes accepted off-diagonal data from
that actual begin. This matters for unfinished parallel fast-factor states whose
retained descendants can begin before the pivot that finally failed. The
one-thread non-contiguous smoke variant requires the serial exact-mask restart,
while the two-thread variant still requires the row-pipeline path. This closes a
serial wrapper mismatch; the remaining paper gap is still CKTSO's full
multi-task ETree-descendant pivoting-tail scheduler and production row/segment
numeric storage.

The masked row-pipeline repair now makes the same ready-supernode branch
mask-aware. Packed active tail rows no longer force singleton producer rows:
the phase initializes supernode metadata from already-preserved rows, each
ordered active-row publish rebuilds the completed local supernode map across
preserved gaps and active rows, and readiness checks use the active-rank map to
stop a run at unfinished active rows while still admitting ready preserved gap
rows. The phase-local dense/common-tail panel cache now uses the same completed
row map in masked tails: preserved-row panels are built at phase start, active
row publication rebuilds the cache across preserved gaps and completed active
rows, and dynamic column exchanges refresh against the current row-major `U`
before speculative suffix rows retry. This closes the direct Algorithm 4/5
panel-use gap for completed dependencies in non-contiguous tail masks, while
still stopping at unfinished active rows through the active-rank guard. The fast
block-repair pipeline reports this through
`fast_kls_block_restart_last_row_pipeline_supernode_update_groups`,
`fast_kls_block_restart_last_row_pipeline_supernode_update_rows`,
`fast_kls_block_restart_last_row_pipeline_supernode_panel_update_groups`, and
`fast_kls_block_restart_last_row_pipeline_supernode_panel_update_rows`; the
non-contiguous gap smoke fixture now requires one grouped update and one
panel-backed update over two producer rows.

The same KLS-owned row-first tail repair now keeps the accepted off-diagonal
refresh inside the proven tail scope. After accepting a row-first repaired
block, KLS compares the repaired prefix `Pnum` and block-local `Q` order against
the saved prefix. When that prefix is unchanged, it refreshes `Offp`/`Offi`/
`Offx` only from the rejected pivot onward; if the proof or suffix refresh
fails, it falls back to the existing full offdiag rebuild. This applies the
same CKTSO-tail storage boundary that the serial tail restart already used to
the parallel row-first repair path. The non-contiguous mask and BTF suffix smoke
fixtures now require `fast_repaired_last_offdiag_suffix_refresh=1`, one suffix
refresh, and zero full offdiag refreshes.

The retained row-refactor input-target map now covers independent batch groups
as well as dense, sparse, compact-panel, and scalar/generic rows. Batch-group
symbolic setup classifies each raw input entry once as a batch dependency-vector
slot, pivot, or row-major `U` workspace slot. Ragged, contiguous multi-producer,
and fragmented multi-producer batch executors consume those destinations before
falling back to the older per-refactor pattern search. The existing batch smoke
fixtures now require `row_refactor_last_segment_target_input_rows` to cover the
executed batch rows in both checked and unchecked modes. This removes another
KLU-shaped rediscovery step from the KLS row/segment numeric path, while the
larger paper gaps remain the production compact/batched row-major numeric engine
and CKTSO's full ETree-descendant pivoting-tail scheduler.

The exact-pattern compact dense batch executor now consumes the same retained
compact-panel input targets before applying producer batches. Exact dense
consumer row batches no longer stage raw input through the sparse work vector
only to fill the compact panel; raw external dependencies still enter `x`,
while in-panel `L`, pivot, dense-`U`, and trailing-`U` values are direct-loaded
into the retained panel and combined with later residual updates. The exact and
subrange compact-supernode smoke probes require compact direct-input and
retained target counters to cover executed batch rows. This closes another
storage-path mismatch; the larger open items remain production batched
row-major kernels and the full CKTSO tail scheduler.

The same retained compact-panel input contract now covers ragged and fragmented
dense-consumer batches. The single-producer ragged path and the fragmented
multi-producer dense-group path both direct-load raw `L`, pivot, dense-`U`, and
trailing-`U` entries into the retained row-major panel, then merge any residual
work-vector value for that slot before applying producer suffix updates. Their
smoke fixtures now require compact direct-input rows and retained target-input
rows to cover the executed compact-supernode batch rows, including the checked
fragmented target-map case. This removes another KLU-shaped staging rule from
the SubtreeLU-style compact batch executors without claiming to close the
remaining production-kernel and CKTSO-tail work.

Compact dense row-refactor groups now have a portable blocked panel factor path
even when CBLAS is disabled. After the retained panel is loaded, KLS factors
diagonal panel blocks row-by-row, solves the below-block multipliers against the
just-factored upper block, and applies one dense/trailing right-looking update
before publishing the row-major values. The checked/pivoting variant preserves
the existing multiplier and pivot acceptance tests, publishes each accepted
diagonal row as soon as that row is final, and leaves fast-reject prefix
semantics intact. This fills a direct SubtreeLU/CKTSO algorithm gap: compact
supernodal panel arithmetic no longer exists only as a CBLAS probe or as per-row
scalar updates. Benchmark JSON and stats expose this through
`row_refactor_last_compact_dense_panel_blocked`,
`row_refactor_compact_dense_panel_blocked_run_count`,
`row_refactor_compact_dense_panel_blocked_rows`, and
`row_refactor_compact_dense_panel_blocked_entries`; both unchecked and checked
dense compact smoke coverage force CBLAS off and now require the blocked path.
The larger open gaps remain the production checked/pivoting tail executor and
broader row-major numeric engine.

KLS then removed a non-paper threshold from the scalar compact-supernode update
selector. SubtreeLU's row update branch treats a ready supernode as a triangular
solve plus trailing update; the previous default only used KLS's contiguous
worker-scratch `trsv` when a work-per-copied-entry gate said it would amortize.
The default now runs that compact `trsv` for every ready producer run with at
least two rows, while `KLS_ENABLE_COMPACT_SUPERNODE_TRSV=0` remains available as
an explicit A/B disable. Existing compact-panel smoke coverage already unsets
the variable and requires the `row_refactor_last_compact_supernode_trsv` counters
to fire, so the test now proves the paper-shaped default rather than a forced
probe.

KLS-first row-up and pivot-tail pipeline phases now also have a portable
cached-supernode panel executor. When a published U-row run has been retained as
a dense panel, the non-CBLAS path solves all ready predecessor multipliers as one
upper-triangular run and applies the common dense/trailing update in aggregate
before falling back to the older per-row panel walk. This moves the
SubtreeLU/CKTSO row-first path closer to the papers' supernode update branch
without claiming the full pipelined ETree-descendant tail factorization is
complete.

The row-first supernode consumer now publishes a compact-validated producer run
before its first consumer falls back to the scalar/compact walk. That makes the
same first consuming row use the cached dense/common-tail panel solve/update
path when the ready run already has the SubtreeLU-style panel shape; the new
single-consumer KLS-first smoke fixture requires exactly one panel-backed update
over the two ready producer rows. This fills the first-consumer half of the
paper ready-supernode branch. The remaining row-major gap is still larger:
production coarse supernode storage/kernels and CKTSO's full pivoting
ETree-tail executor are not complete.

The row-first failed-pivot repair then stopped requiring the retained
pivoting-tail plan to start exactly at the rejected pivot before the
restartable pipeline could use it. KLS now derives a local repair envelope from
the retained topological tail's first and last columns, validates that the
envelope contains the rejected pivot, and either runs the contiguous envelope or
uses the active tail mask for non-contiguous gaps. The reusable-prefix proof,
off-diagonal suffix refresh, and repaired-tail work accounting all start at the
actual retained-tail begin rather than blindly at the rejected pivot. The
non-contiguous gap smoke fixture now requires an exact repaired tail mask, zero
overcompute columns, and suffix-only off-diagonal refresh. This is a direct
CKTSO Algorithm 5 gap closure for the current KLS-owned row-first repair path;
it still is not the full parallel ETree-descendant pivoting-tail scheduler.

KLS now fills another CKTSO Algorithm 5 semantic gap in how that retained tail
is seeded. The CKTSO paper says that, after a pivot-check interruption, the
tail starts from the unfinished guessed-EGraph nodes and their ETree closure.
KLS previously refreshed missing prefix rows first and then built the seed from
the remaining suffix or a narrower row-tail heuristic, which could erase the
actual interrupted unfinished set. The checked row and mapped EGraph reject
paths now snapshot the whole rejected block's unfinished done-bitmap before
prefix refresh, guarantee that the rejected pivot is in the seed, allow seed
columns before the rejected pivot, and prefer that seed before row-tail/suffix
fallbacks. The parallel checked-row smoke fixture now requires the pivoting-tail
seed to come from that unfinished set, with row-tail seeding bypassed. This
still leaves the larger CKTSO executor gap open: the retained worklist is used
by KLS's guarded row-first repair path, not by a full production tail scheduler.

KLS now also fills the NICSLU Algorithm 4 performance-model gap. After a
numeric factorization, it walks the actual U-dependency graph in factor order,
assigns each column to the thread with the smallest current `END(p)`, charges
`2*nnz(L(:,i))` work for each dependency update, `nnz(L(:,k))` work for
normalization, and a unit sync cost per dependency, then reports the predicted
finish time and flop-only speedup. The model is exposed through
`parallel_task_flow_*` stats and benchmark JSON, and can seed the same
row/segment metadata preparation path as the earlier R1/R2 NICSLU suitability
counters. This is still a model and policy input, not NICSLU's full
ETree/EScheduler-guided factorization executor.

KLS-first pivoting row-up factorization now consumes the retained separator
tree through the same Algorithm 6-style split/collapse idea used by the
row-refactor separator queue. The first-factor planner computes block-local
row-input work per retained component, splits dominant subtrees into pipeline
roots plus private child-subtree candidates, assigns those candidate subtrees
to private thread queues by work, and validates private ownership against the
original row dependencies before remapping the block. Invalid partitions fall
back to the older retained-component queue instead of forcing the generic row
pipeline. Benchmark JSON and smoke coverage expose the path through
`kls_first_last_separator_queue_partitioned`,
`kls_first_separator_queue_partitioned_count`, and
`kls_first_last_separator_queue_split_components`. This closes the direct
KLS-first gap where the pivoting first factor had Algorithm 3 execution but not
the paper's Algorithm 6 separator split; checked-tail/refactor pivoting
consumers and production coarse supernode storage remain open.
