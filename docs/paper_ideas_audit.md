# Paper Ideas Audit

This note records which ideas from the reference papers are present in KLS and
which remain open. The intent is to keep KLS development focused on general
solver algorithms instead of tuning individual benchmark matrices.

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
  improvement for weak or high-off-diagonal-pivot medium matrices, and keeps the
  permutation only when numeric quality and cost evidence justify it.
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
  permutation, but not a full MC64-style maximum product matching with both row
  and column scaling vectors.
- NICSLU/CKTSO parallel scheduling is only present at BTF-block granularity.
  KLS does not yet have an intra-block EGraph/ETree cluster/pipeline scheduler.
- CKTSO fast factorization is present as pivot-checked reuse plus full fallback.
  KLS does not yet implement CKTSO's pipelined tail factorization that restarts
  from the ETree descendants after a failed pivot check.
- SubtreeLU-style nested-dissection metadata is only used indirectly through
  METIS orderings and CAMD refinement. KLS does not yet preserve a separator
  tree for private/pipeline task queues.

## Not Implemented Yet

- Full MC64-equivalent weighted matching and diagonal row/column scaling.
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
