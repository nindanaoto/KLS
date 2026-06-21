# Paper Ideas Audit

This note records which ideas from the reference papers are present in KLS and
which remain open. The intent is to keep KLS development focused on general
solver algorithms instead of tuning individual benchmark matrices.

## Current Conclusion

KLS has **not** implemented every paper idea that is still worth trying. It has
implemented the ideas that can be layered around the current KLU-derived
Gilbert-Peierls kernel: BTF, AMD/COLAMD/METIS ordering policy, CAMD refinement,
auto scaling policy, pivot-checked reuse, static row-pivoting trials with
matching-derived equilibration, and BTF-block parallel refactorization.

The remaining worthwhile ideas are not per-matrix tuning knobs. They require
new KLS-owned symbolic/numeric machinery:

- Full MC64-equivalent maximum-weight matching with dual row/column scaling.
- Persistent EGraph/ETree or separator-tree metadata for intra-block parallel
  factor/refactor scheduling.
- CKTSO-style fast factorization tail restart after a failed pivot check, rather
  than full fallback from the beginning.
- SubtreeLU-style private/pipeline scheduling from a retained separator tree.
- A structure-adaptive triangular solve built on LU storage that exposes cheap
  row/segment access.

Several smaller dispatch experiments were tried and rejected because they helped
some benchmark cases while regressing others. Those are documented below so the
project does not drift toward benchmark-name-specific heuristics.

## Implemented

- KLU Algorithm 907 baseline: BTF preprocessing, AMD/COLAMD ordering,
  KLU-style row scaling modes, Gilbert-Peierls factorization with partial
  pivoting, no-pivot refactorization, block back substitution, diagnostics, and
  pivot tolerance controls are available through the vendored Trilinos
  SuiteSparse-derived sources.
- Nested-dissection ordering option: KLS vendors METIS/GKlib as pinned
  submodules by default and can optionally use a compatible system METIS.
- Combined ordering policy: KLS auto mode can choose AMD, COLAMD, or METIS,
  retry no-BTF symbolic analysis for structural cases where BTF is not useful,
  and promote expensive numeric factorizations to METIS when actual fill/flop
  evidence is better.
- Constrained nested-dissection refinement: KLS can refine METIS rank groups
  with CAMD constraints, preserving nested-dissection rank shape while reducing
  local fill and flops.
- Scaling policy: KLS exposes KLU scale modes and has auto scale selection based
  on pattern and numeric evidence, including no-scale reuse where repeated
  SPICE refactorization benefits.
- Fast repeated factorization with pivot check: KLS reuses an existing numeric
  pattern, checks reused pivots against the selected threshold via L
  multipliers, and falls back to full pivoting factorization if the reused order
  is unsafe.
- Static pivoting trial: KLS has value-aware greedy row matching and swap
  improvement for weak or high-off-diagonal-pivot medium matrices. It can also
  trial matching-derived row/column equilibration and keeps the transformed
  candidate only when numeric quality and cost evidence justify it.
- Initial KLS-owned parallelism: repeated refactorization can run across
  independent BTF diagonal blocks for large high-flop cases where the block work
  is wide enough to offset thread overhead.
- SPICE-cycle orientation policy: KLS can analyze normal and transposed storage
  orientations and select the faster internal form for repeated solve cycles.
- LGPL project licensing and third-party notices: KLS itself is
  LGPL-2.1-or-later, with vendored-source attribution separated from the KLS
  license.

## Partially Implemented

- CKTSO-style static pivoting is only partial. KLS has a practical weighted row
  permutation and matching-derived row/column equilibration, but not a full
  MC64-style maximum-product matching algorithm with assignment dual scaling.
- NICSLU/CKTSO parallel scheduling is only present at BTF-block granularity.
  KLS does not yet have an intra-block EGraph/ETree cluster/pipeline scheduler.
- CKTSO fast factorization is present as pivot-checked reuse plus full fallback.
  KLS does not yet implement CKTSO's pipelined tail factorization that restarts
  from the ETree descendants after a failed pivot check.
- SubtreeLU-style nested-dissection metadata is only used indirectly through
  METIS orderings and CAMD refinement. KLS does not yet preserve a separator
  tree for private/pipeline task queues.

## Not Implemented Yet

- Full MC64-equivalent weighted matching and assignment-dual row/column scaling.
- NICSLU static-symbolic R1/R2 performance model for choosing sequential versus
  parallel numeric kernels before factorization.
- Intra-block parallel factorization with pivoting scheduled by an ETree.
- Intra-block parallel no-pivot refactorization scheduled by an exact EGraph.
- CKTSO dual-mode cluster/pipeline fast factorization with pivot check.
- CKTSO pipelined tail factorization with pivoting after a pivot-check failure.
- CKTSO structure-adaptive hybrid parallel triangular solve.
- SubtreeLU separator-tree collapse for pivoting factorization.
- SubtreeLU FLOP-balanced separator-tree partitioning for refactorization.
- SubtreeLU constrained pivot search within nested-dissection subdomains.
- SubtreeLU supernodal updates with BLAS-style kernels inside the sparse
  up-looking framework.

## Tried During This Audit

An intra-block levelized no-pivot refactor prototype was tested using the
existing U structure as the exact dependency graph, matching the NICSLU/CKTSO
EGraph refactorization idea. It was removed before commit because the general
implementation regressed dominant-block circuit cases where the current
BTF-block threaded refactor is faster. Keeping it would have required
case-specific dispatch, which is not the desired direction for KLS.

A KLS-owned serial no-pivot refactor path was also prototyped by reusing the
threaded BTF-block refactor kernel when thread-level parallelism was not
eligible. It passed correctness tests and helped some no-scale cases, but it
regressed other representative circuit cases such as `bcircuit` and `rajat03`.
The prototype was removed because it did not represent a general improvement.

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

## Recommended General Work

1. Implement a real matching/scaling stage first, because MC64-style static
   pivoting appears in both NICSLU and CKTSO and can reduce dynamic pivoting
   before any parallel scheduler is added.
2. Build persistent dependency metadata for one KLS-owned numeric engine:
   either CKTSO's EGraph/ETree cluster/pipeline scheduler or SubtreeLU's
   separator-tree private/pipeline scheduler, not isolated benchmark guards.
3. Add a structure-adaptive triangular solve only after the LU storage owned by
   KLS exposes row-oriented or segmented access cheaply.
4. Use static symbolic and numeric-cost models to decide whether a parallel
   kernel should run, so KLS avoids matrix-name-specific tuning.
