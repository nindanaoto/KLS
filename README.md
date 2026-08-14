# KLS

KLS is a standalone sparse direct solver project targeting SPICE-style
workloads: one symbolic analysis followed by many numeric factorizations,
refactorizations, and solves with a fixed sparse structure.

This repository currently contains the first working KLS implementation:

- A stable C API in `include/kls/kls.h`
- A vendored SuiteSparse-derived 64-bit symbolic/numeric engine
- CSC and CSR input paths with 32-bit or 64-bit index arrays
- Evidence-driven automatic symbolic ordering across AMD, COLAMD, AMMF, AMF3,
  and guarded nested dissection, with explicit AMD, COLAMD, natural, METIS,
  SCOTCH, AMF, AMMF, and AMF3 controls
- SPICE-cycle-oriented normal-vs-transpose internal orientation selection, with
  explicit orientation controls
- Factor, refactor, solve, transpose-solve, and statistics APIs
- Fast repeated factorization that reuses the existing numeric pattern, checks
  pivot quality, and can repair unscaled rejected BTF blocks before falling
  back to full pivoting factorization
- A pipelined parallel first factorization for large blocks (deterministic
  column pipeline with lockstep panels), with a BLAS3 dense-tail finish for
  extreme-fill patterns
- Value-aware static row pivoting and row/column equilibration trials for high
  off-diagonal-pivot cases
- A MatrixMarket benchmark tool
- A small correctness smoke test
- A SuiteSparse Matrix Collection downloader script for public benchmark cases

## AUTO policy

The production `AUTO` path is matrix-family agnostic. It compares affordable
orientation, ordering, BTF, matching, scaling, and numeric candidates, then
uses realized fill, factor work, pivot behavior, residual checks, dependency
parallelism, and timed recurring kernels to retain a choice. Small problems
use resource-based amortization floors so candidate setup cannot dominate the
numeric work. Fixed-width index limits and storage prerequisites remain hard
capability gates.

In particular, a retained row-matched factor can compare both minimum-degree
orderings and switch only on numeric Pareto evidence (fill, work, pivot detours,
and conditioning). A rejected checked-factor suffix compares its measured
repair lower bound with the measured full rebuild instead of entering a
family-specific repair engine. Verified subtree solve plans and equivalent
EGraph dispatch, fusion, indexed-scatter, and supernode-tail kernels are timed
on the realized factor; conservative margins keep the incumbent on close
draws. These consultations are cached for that numeric lifecycle.

For a caller-declared repeated lifecycle, AUTO also admits one post-orientation
nested-dissection symbolic candidate when the selected minimum-degree estimate
has enough storage and arithmetic to amortize it. The pre-numeric choice is
predicted only when the symmetric union is bounded, where the row/column
symbolic frame has a finite expansion certificate, and is skipped when a
pending generic block order will replace it. An unbounded frame can still enter
as a numeric-provisional candidate when at least 16 refactors and (10^{12})
estimated lifecycle flops can repay an ordinary pivoted trial. That candidate
must realize less than 95% of the recorded fallback fill; rejection restores
the minimum-degree symbolic and factors it normally. Within this high-work
provisional path, the incumbent's realized flops-per-fill ratio selects between
coarse and fine resource-scaled CAMD windows, while the exact factor remains
the acceptance authority. Bounded candidates retain their predicted-fill
validation. These are lifecycle, representation, and measured-factor contracts
rather than dimension or matrix-family selectors. Set
`KLS_DISABLE_GENERIC_ND_PORTFOLIO=1` for A/B runs.

When that repeated lifecycle has already selected nested dissection, a
value-aware static row-match candidate is allowed to settle before the initial
factor instead of waiting for the first changed numeric. It must still pass the
ordinary weak-diagonal admission and measured match/factor quality gates; the
ordering choice and caller lifecycle are the only extra evidence. Set
`KLS_DISABLE_REPEATED_ND_STATIC_SETTLEMENT=1` for a deferred A/B run.

The recurring numeric router is measurement-based as well. A timed lean-row
winner uses compact 16-bit streams whenever every retained index and pointer
fits that representation. After eight direct-worker samples, AUTO replays two
incumbent column numerics and reverses a cold-start row decision only when the
warmed incumbent is at least 2% faster. Symmetrically, a timed row decline
retains its already-built best arm; after eight real incumbent calls, AUTO
replays that arm twice and adopts it only for a 5% win. These audits use
realized runtimes and representation capability, not dimensions or a matrix
profile. The individual A/B controls are `KLS_DISABLE_LEAN_I16_INDICES=1`,
`KLS_DISABLE_LEAN_STEADY_REAUDIT=1`,
`KLS_DISABLE_DECLINED_LEAN_STEADY_REAUDIT=1`, and
`KLS_DISABLE_SETTLED_LEAN_DIRECT_REFACTOR=1`.

The cooperative row/column tournament is front-loaded into the first changed
numeric. A decisive cold result settles immediately; an inconclusive
model-selected row takes one identically positioned warm row/column pair.
Material warm losses seed the existing solve-aware tournament instead of
occupying several public refactors. Once the following column solve arrives,
row is rejected whenever its refactor alone exceeds the complete measured
column refactor-plus-solve cycle; otherwise both representations continue to
receive paired samples. One complete row pair that is already 25% slower may
also settle early, whereas positive row adoption still requires two solve
samples per side. Provisional winners are re-audited after four real row cycles
with two warmed column refactor/solve cycles. All of these decisions use
realized lifecycle time and representation validity, not matrix dimensions or
sparsity fingerprints.

The older profile selectors are retained only for reproducible A/B work. Set
`KLS_ENABLE_LEGACY_SHAPE_POLICIES=1` to enable them; they are off by default.
This compatibility switch covers the historical dimension/degree/BTF windows,
direct family routes, family-specific EGraph layouts, and residual-check
bypasses. Explicit API choices and explicit experiment environment variables
continue to take effect in either mode.

Callers that know their numeric lifecycle can set
`expected_refactorizations` and `expected_solves` in `kls_options`. These are
performance hints only; they do not weaken factor or solve correctness checks.
On heterogeneous-cache Linux systems, the implementation may keep a repeated
factor's cooperative worker pool on distinct physical cores in the largest LLC
domain when its retained working set exceeds a smaller eligible LLC. The
choice uses the caller's CPU set, hardware topology, cache capacity, realized
factor storage, and lifecycle hint—not matrix dimensions or a matrix-family
detector. It is a no-op on uniform-cache systems and when the work fits every
eligible LLC. Set `KLS_DISABLE_COMPACT_LLC_AFFINITY=1` for scheduler-controlled
placement, or `KLS_TRACE_AFFINITY=1` to report an adopted placement.

When a recurring solve contract needs a snapshot of the current coefficient
values, an already-active numeric worker pool also copies disjoint value
slices while completing the refactor. Set
`KLS_DISABLE_PARALLEL_REFINE_COPY=1` for a serial-copy A/B comparison.

Many later sections in this README document the development history of those
profile selectors. Unless a section explicitly describes a capability- or
measurement-based default, treat its named AUTO policy as compatibility-mode
documentation rather than the current production router.

## Build

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

The main build is reproducible from pinned submodules. KLS vendors the
SuiteSparse-derived KLU, AMD, COLAMD, and BTF C sources from Trilinos under
`third_party/suitesparse` as the current in-tree serial engine. METIS ordering
is enabled by default from pinned submodules under `third_party/metis` and
`third_party/gklib`, SCOTCH ordering is enabled by default from
`third_party/scotch`, and BSD-licensed SPRAL matching/scaling is enabled by
default from `third_party/spral`; initialize them with
`git submodule update --init --recursive` after cloning. To use a compatible
system METIS instead, configure with:

```sh
cmake -S . -B build -DKLS_USE_SYSTEM_METIS=ON
```

KLS builds METIS with 64-bit `idx_t` for compatibility with the 64-bit KLS/KLU
path. A system METIS install used this way must be ABI-compatible.

To use a compatible system SCOTCH instead of the pinned submodule, configure
with:

```sh
cmake -S . -B build -DKLS_USE_SYSTEM_SCOTCH=ON
```

A static system SCOTCH build may also require `scotcherr`; pass
`KLS_SYSTEM_SCOTCHERR_LIBRARY` if it is not discoverable.

KLS builds BSD-licensed SPRAL Hungarian/auction matching/scaling support from
the pinned `third_party/spral` submodule by default. To disable this MC64-adjacent
component and keep a C-only build, configure with:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DKLS_ENABLE_SPRAL_SCALING=OFF
cmake --build build
```

To use a system SPRAL library instead of the bundled submodule:

```sh
cmake -S . -B build -DKLS_ENABLE_SPRAL_SCALING=ON \
  -DCMAKE_BUILD_TYPE=Release \
  -DKLS_USE_SYSTEM_SPRAL=ON \
  -DKLS_SYSTEM_SPRAL_LGPL_COMPATIBLE=ON \
  -DKLS_SYSTEM_SPRAL_INCLUDE_DIR=/path/to/spral/include \
  -DKLS_SYSTEM_SPRAL_LIBRARY=/path/to/libspral.so
```

`KLS_SYSTEM_SPRAL_LGPL_COMPATIBLE=ON` is an explicit acknowledgement that the
selected system library is redistributable with LGPL-2.1-or-later KLS. Do not
use it for HSL MC64 itself or for solver-tree MC64 copies that retain HSL-style
redistribution restrictions.

This is not a solver replacement. KLS uses it as a
static-pivot matching fallback when the in-tree matcher is short of a full
cardinality match, as a guarded pre-factor Hungarian/scaling path for large
weak-diagonal dominant-BTF matrices, and as a guarded post-factor Hungarian
trial for dense high-off-diagonal-pivot cases after cheaper scale, ordering,
and pivot-tolerance fixes have run. Accepted factorization candidates still
have to pass KLS's normal numeric checks and a value gate that requires the
matching to remove substantial pivoting pressure or materially reduce factor
work/fill.

For a compact normal-AMD one-block symbolic with at least 31/32 of rows both
weak and structurally missing on the diagonal, KLS can run that guarded match
before the first numeric rather than deferring it to the first changed-value
refactor.  The timing policy is bounded to 4,096--16,384 rows, 3--8 input
entries per row, 16--64 estimated factor entries per row, balanced L/U, and
128--4,096 estimated flops per row inside a 1M--64M resource band.  Matching
coverage, candidate factorization, and the ordinary numeric acceptance gates
still decide whether the transformed candidate is retained.  Benchmark JSON
reports `nearly_missing_diagonal_early_match_selected`; set
`KLS_DISABLE_NEARLY_MISSING_DIAGONAL_EARLY_MATCH_POLICY=1` to restore deferred
evaluation.  The former Hamrle2-specific input check is removed.

AUTO inputs that can carry the complete compact matching lifecycle now use a
second staged policy instead of a `gemat11`/`gemat12` dimension and nonzero
fingerprint.  The input proposal requires 4,001--65,535 rows and at most
65,535 entries, 2--16 entries per row overall, no empty row or column, row and
column degrees at most 64, and at most 1/32 of columns represented on the
structural diagonal.  A lightweight NATURAL placeholder is replaced only
after a value-aware match covers at least 99.5% of rows and the matched
symbolic proves full rank with a block spanning at least 3/4 of the matrix.
Inputs below four entries per row retain the full weighted matcher, where
alternative-edge quality matters; denser inputs use the lower-overhead compact
matcher.  The accepted one-block holdout keeps the mapped update engine; the
specialized compact row/direct-value lifecycle additionally requires at least
64 BTF blocks, balanced L/U streams, at most 16 factor entries and 128 measured
flops per row, at most `n/64 + 16` off-diagonal pivots, and no scaling,
nudging, perturbation, or predicted numeric.  This separates matching
economics from recurring-engine eligibility.

`kls_stats` and benchmark JSON expose
`compact_missing_diagonal_match_candidate`,
`compact_missing_diagonal_match_selected`, and
`compact_missing_diagonal_factor_eligible`.  Set
`KLS_DISABLE_COMPACT_MISSING_DIAGONAL_MATCH=1` to disable the complete policy;
the former matrix-named `KLS_*GEMAT*` switches are no longer recognized.
`Hamrle2` is included in the extended SuiteSparse manifest as a
cross-family one-block holdout.  The exact Gemat selector and Gemat-named
internal representations are removed.

A separate symmetric partial-diagonal policy replaces the former OPF_3754 and
OPF_10000 size windows.  Its input proof covers 8,192--131,072 rows, 6--12
entries per row overall, 2--64 entries in every column, exact
multiplicity-aware structural symmetry, and structural diagonals in 40--60%
of columns.  The heavier match proposal starts only above 262,144 input
entries.  It is retained only when a value-aware match produces normal AMD,
one unblocked full-rank component, balanced symbolic L/U, 8--24 estimated
factor entries per row with at least 400,000 entries, and 32--512 estimated
flops per row with at least four million operations.  The final numeric must
independently meet the corresponding fill/work and pivot-repair bounds before
direct-value, PTS, preparation-overlap, or settled-probe consumers run.

Smaller members use an independently bounded AMD/BTF lifecycle: 100,000--
500,000 estimated factor entries and one--eight million estimated operations.
That immutable symbolic verdict is cached after adoption, so the generic proof
is not repeated in refactor or solve dispatch.  Direct updates consume public
values only after the first solve contract succeeds, and both PTS and the
serial one-block fallback apply the inverse input map and paired match scales.
Set `KLS_DISABLE_SYMMETRIC_PARTIAL_DIAGONAL_MATCH_POLICY=1` to disable the
whole policy or `KLS_DISABLE_SYMMETRIC_PARTIAL_DIAGONAL_DIRECT_VALUES=1` to
retain matching while disabling fused direct input.  `kls_stats` and benchmark
JSON expose `symmetric_partial_diagonal_match_candidate`,
`symmetric_partial_diagonal_match_selected`,
`symmetric_partial_diagonal_factor_eligible`, and
`symmetric_partial_diagonal_low_work_eligible`.

For an AUTO/8-thread normal-AMD symbolic whose full-rank BTF consists of many
genuinely tiny components, KLS avoids global block-ordering work that costs
more than the complete factor.  The structural contract covers
16,384--131,072 rows, 6--12 input entries per row, one BTF block per 2--8
rows, a largest block no wider than 256 or 1/128 of the matrix, 4--12
estimated factor entries per row, and 4--64 estimated flops per row below an
8M-work ceiling.  Recurring matching, tolerance, solve, and accuracy choices
use a stricter factor contract: no coordinate transform, no nudging or
perturbation, at most one off-diagonal pivot per 16 rows, 3--12 factor entries
and 3--64 flops per row.  A high-pivot value set therefore keeps the useful
matching trial even when its pattern shares the structural class.  Benchmark
JSON reports `low_work_tiny_block_btf_symbolic_eligible` and
`low_work_tiny_block_btf_policy_eligible`; set
`KLS_DISABLE_LOW_WORK_TINY_BLOCK_BTF_POLICY=1` to restore all generic trials.
The former Sandia operating-point input box is removed.

For an AUTO/8-thread normal-AMD, max-row-scaled factor with one moderate BTF
core and a fragmented fringe, KLS selects its compact paired-row lifecycle
only when the complete representation and measured work fit.  The contract
covers 8,192--65,535 rows, 3--8 input entries per row, `n/8`--`n/3` BTF
blocks, a largest block spanning one half to three quarters of the matrix,
8--20 estimated factor entries and 32--256 estimated flops per row, and
16-bit L/U index streams.  The retained factor must have 5--10 entries and
16--64 flops per row, no transforms, nudges, perturbations, or predicted
numeric, and at most one off-diagonal pivot per 32 rows.  Passing factors can
avoid redundant tolerance, value-publication, reciprocal, and residual work;
nearby scaled factors that exceed the representation or pivot bounds retain
the generic lifecycle.  Benchmark JSON and `kls_stats` report
`scaled_fragmented_compact_row_policy_eligible`; set
`KLS_DISABLE_SCALED_FRAGMENTED_COMPACT_ROW_POLICY=1` to restore the generic
trials, or
`KLS_DISABLE_SCALED_FRAGMENTED_COMPACT_ROW_DEFER_VALUE_SCATTER=1` to retain
the policy while disabling only deferred publication.  The former rajat27-
specific input box is removed.

Compact 16-bit triangular solves now cache consecutive singleton BTF runs
from the selected block structure instead of recognizing two Rommes input
boxes.  The two-byte-per-block table is built only when block pointers and
factor indices fit the 16-bit representation, there are at least 512 blocks,
at least 75% are singletons, and one consecutive run spans at least 256
blocks.  An adopted table also enables refreshed diagonal reciprocals from
the representation itself.  This applies to explicit as well as automatic
BTF choices and to scaled or unscaled numerics; factors that do not build the
compact mirror pay no table cost.  Set
`KLS_DISABLE_COMPACT_SINGLETON_RUN_SOLVE=1` to disable the capability.
`kls_stats` and benchmark JSON report structural eligibility, the number of
cached singleton blocks, and the longest cached run.  The former Itaipu and
MIMO8 dimension/nonzero/block fingerprints and their unrelated row-engine
overrides are removed.

For very large pre-static candidates, KLS uses SPRAL's auction
matching/scaling path instead of exact Hungarian matching, except on
mostly-missing-diagonal patterns where the auction predictably falls short of
the coverage floor and the exact Hungarian matching is both complete and
faster; those go straight to the exact matching. Both keep the MC64-adjacent
step on LGPL-compatible redistributed code.

KLS keeps this as the only vendored MC64-adjacent external implementation.
This is a license boundary, not an authorship boundary: existing MC64-style code
is acceptable when its license remains compatible with KLS's LGPL distribution
goal, allows redistribution in source and binary form with KLS, and allows KLS
to preserve the upstream license notices in `THIRD_PARTY_NOTICES.md`. HSL MC64
itself, and solver-tree copies that retain HSL redistribution restrictions, are
not compatible with that boundary. A permissively licensed translation of
SPRAL's scaling code can be used as a reference, but the pinned SPRAL submodule
is the preferred reproducible source for the C/Fortran build.

The boundary is checked by the normal CTest suite when Python is available and
can also be audited directly:

```sh
python3 scripts/audit_license_boundary.py
```

## Benchmark

```sh
./build/kls_bench ../cktso/demo/add20.mtx --repeat 20 --refactor-repeat 20 --orientation auto --json
```

The benchmark reports analysis, factorization, refactorization, solve,
transpose-solve, residual, selected orientation, BTF block/rank, fill, flop,
the initial and last factorization path, largest-BTF-block factor ETree
shape, refactor dependency-level metrics, dependency root/leaf/max-fanout
scheduler diagnostics, NICSLU-style `parallel_model_r1`,
`parallel_model_r2`, `parallel_model_recommends_parallel`, and numeric
task-flow model metrics, selected input-index width, and memory statistics. By
default `kls_bench` uses `--input-index auto`, which passes 32-bit CSC indices
when the MatrixMarket problem fits the public `KLS_INDEX_INT32` API and falls
back to 64-bit otherwise; use `--input-index 64` or `--input-index 32` for
forced A/B runs. Use `--analyze-only` to measure symbolic analysis and ordering
decisions without running numeric factorization. Use `--structure-only` to stop
after MatrixMarket cleanup/deduplication and report exact CSC diagonal coverage,
empty/scalar row and column counts, maximum row and column degrees, and total
row/column-degree mismatch. The two diagnostic modes are mutually exclusive.
JSON fields
`compact_solve_index_bytes` and `compact_solve_fused_rhs` report whether a
2- or 4-byte triangular-solve mirror was prepared and whether its public-RHS
permutation was precomposed; zero index bytes means the native factor storage
was retained.  `compact_solve_singleton_run_eligible`,
`compact_solve_singleton_run_blocks`, and
`compact_solve_singleton_run_max` expose the generic singleton-run candidate
and the representation actually adopted after a solve.

Use `--no-transpose-solve` when a comparison protocol scores only normal
solves.  The paired-suite runner supplies it because the CKTSO, SubtreeLU, and
KLU harnesses do not execute a second, unscored transpose-solve loop.  The
standalone benchmark keeps measuring both directions by default.

Exactly singular matrices report `KLS_ERR_SINGULAR` by default, regardless of
their dimensions or selected AUTO policy.  Callers that deliberately want a
checked rank-completion path can opt in with
`KLS_ENABLE_SINGULAR_COMPLETION=1`.  It records exact zero pivots, constrains
the corresponding degrees of freedom (through matched stored entries when
available), retains the true values for residual refinement, and carries those
constraints across refactorization.  `KLS_DISABLE_SINGULAR_COMPLETION=1`
overrides the opt-in. If analysis has already proved a structural-rank deficit
and completion is not enabled, AUTO also suppresses its background METIS race:
ordering cannot repair the deficit, and no resulting factor could be adopted.
The race remains available when singular completion is explicitly enabled.

### Opt-in lean serial backend

For one-thread SPICE workloads, `KLS_BACKEND_SERIAL` skips the adaptive
ordering/matching/promotion pipeline and uses a deliberately small policy:
direct AMD analysis by default, an unscaled pivoting factor first, and retained
mapped or structurally gated supernodal refactors. It is opt-in; `AUTO` and
`KLS` retain the existing adaptive behavior, including at multiple threads.

```c
kls_options options;
kls_default_options(&options);
options.backend = KLS_BACKEND_SERIAL;
options.threads = 1;
```

The same mode is available to the benchmark tools:

```sh
./build/kls_bench matrix.mtx --backend serial --threads 1 \
  --repeat 1 --refactor-repeat 5 --json
```

Serial mode rejects more than one thread and supports AMD, COLAMD, or natural
ordering. With automatic orientation it uses normal CSC and transpose CSR so
it does not pay for an orientation race. With automatic scaling it tries the
lean unscaled factor first and retries with row-sum scaling only after a factor
failure/singularity or measured numerical work inflation. See
[`docs/serial_backend.md`](docs/serial_backend.md) for the policy, validation
split, limitations, and reproducible corpus results.

When system SuiteSparse KLU headers and libraries are installed, the build also
provides `klu_width_compare` to compare system `klu_*` and `klu_l_*` on the same
MatrixMarket input. This is a diagnostic benchmark for deciding whether a
future 32-bit KLS backend is worth implementing:

```sh
./build/klu_width_compare matrix.mtx --repeat 3 --refactor-repeat 3 --json
```

An exact-source vendored 32/64-bit diagnostic is also available as an excluded
target, so it does not affect ordinary builds:

```sh
cmake --build build --target klu_width_compare_vendored
./build/klu_width_compare_vendored matrix.mtx --repeat 3 \
  --refactor-repeat 3 --json
```

JSON includes `initial_factor_path` and `last_factor_path`; values such as
`klu_first` or `klu_fallback` mean the factorization was handed to the
KLU-derived pivoting kernel, while `kls_fast_refactor` means KLS reused the
retained pattern through the checked fast path. Repeated refactor diagnostics
are reported separately as `last_refactor_path`, distinguishing row-refactor,
EGraph, mapped, pool, and serial KLU refactor branches.

For a repeated, unscaled, direct-CSC numeric whose retained factor is one
structurally full-rank block below the generic 100,000-flop work floor, AUTO
settles directly on the compact map32 fixed-pivot walk. This walk performs the
same arithmetic as the native KLU refactor while using predecoded 32-bit input
positions and one pivot reciprocal. Selecting it from this existing
capability/work contract avoids spending six microsecond-scale updates on the
otherwise generic alternating representation tournament; the ordinary
first-update solve contract still validates the numeric before direct reuse.
`KLS_DISABLE_LOW_WORK_SINGLE_BLOCK_COMPACT_MAP32_POLICY=1` restores the timed
representation tournament for A/B runs, while
`KLS_DISABLE_COMPACT_MAP32_TOURNAMENT=1` disables compact selection entirely.
The two current paper-suite members of this structural class both improved in
31-pass rank-preserving H100 measurements: `add32` from 7.181 ms to 6.850 ms
and `1138_bus` from 1.997 ms to 1.963 ms.

The parallel lean-row executor treats deferred public-to-internal value
gathering and grouped completion tokens as alternative optimizations. Deferred
gathering adds a producer phase before the row dependency stream, whereas the
compact grouped-token proof describes an already-prepared stream; composing
the two could publish an incomplete numeric generation while still returning
success. AUTO therefore retains deferred gathering but uses the ordinary
per-row completion records for that call. This is an executor-representation
contract, independent of matrix dimensions, ordering, or identity. A
24-matrix changed-value audit verified every refactor below the `1e-8` residual
limit and kept KLS ahead of CKTSO on all 24 mutually valid comparisons. In
particular, `gemat11` moved from a `4.44e-5` failing residual to a
`5.02e-16` worst-refactor residual and a 31-pass median H100 of 11.944 ms,
versus CKTSO's saved 17.883 ms.

Unset `KLS_ENABLE_KLS_FIRST_FACTOR` keeps the production cold first factor on the
KLU/static path for broad large cases. `KLS_ENABLE_KLS_FIRST_FACTOR=1` forces
the KLS-owned row-up-looking scaffold when possible, and
`KLS_ENABLE_KLS_FIRST_FACTOR=0` keeps the hard KLU-first behavior. Automatic
mode also keeps successful pre-static row-matching candidates on the accepted
KLU/static numeric instead of replaying them through the incomplete KLS-first
bridge; forcing `KLS_ENABLE_KLS_FIRST_FACTOR=1` keeps that bridge available for
experiments. The
scaffold remains useful for paper-algorithm experiments, but it is not the
default cold first-factor replacement for broad CKTSO-gap cases until it is
faster than the accepted KLU/static numeric there.
Benchmark JSON reports `kls_first` when the forced scaffold successfully
assembles KLU-compatible numeric
storage. For eligible no-scale or KLU row-scaled matrices, that path first
tries a KLS-owned sparse row-major up-looking first factor over each BTF
diagonal block, matching CKTSO Algorithm 1's row orientation before packing the
accepted factors for the existing solve/refactor interfaces. Accepted
pre-static row-matching candidates can also replay their selected symbolic and
value state through this row-up factor when the recovery gate allows it,
instead of returning with a KLU-built trial numeric that is slow or unsafe for
the repeated SPICE cycle. When a row diagonal fails the threshold against the
row's largest
U-tail entry, this bridge exchanges the active block-local column with that
largest entry, publishes the accepted `Q` order, and continues in KLS-owned
row-major storage. Benchmark stats report this direct bridge as
`kls_first_last_row_uplooking_columns` and
`kls_first_row_uplooking_column_count`. A successful row-up-looking first
factor also publishes its row-major `L`/`U` entries, input row positions, and
KLU numeric value pointers directly into the row-refactor metadata finisher;
that skip-over-numeric-scan handoff is reported as
`kls_first_last_row_refactor_seeded_rows` and
`kls_first_row_refactor_seeded_row_count`. BTF cases still build the retained
factor-order input map for off-diagonal refreshes and fallback coherence, but
the row-refactor pattern and value mirrors no longer have to be reconstructed
from packed KLU columns after a successful KLS-first row-up factor. Algorithm
4-style row-supernode dependency runs consumed by KLS-first row-up workers are
reported through `kls_first_last_row_supernode_update`,
`kls_first_row_supernode_update_run_count`,
`kls_first_last_row_supernode_update_groups`, and
`kls_first_last_row_supernode_update_rows`. KLS-first row-up producers publish
completed row-supernodes as cached dense/common-tail panels once the following
row proves the supernode ended; this covers parallel BTF workers plus the
ordinary private and serial row-up loops. In the ordinary row-up path, when a row
first consumes a compact validated producer run before the complete producer
supernode is known, KLS now publishes that ready prefix before falling back to
the older compact walk, so the same first consumer can use the cached panel
solve/update path. Separator private/pipeline row-up paths use the same lazy
publication once a compact prefix has been validated under the scoped pivot
order. Later rows also try those panels before falling back to row-entry
validation, and the cached consumer can use the published prefix of a longer
dependency run while leaving the remaining suffix in the dependency heap. That
matches the paper private/pipeline rule that finished producer prefixes should
be consumed before waiting for or continuing through the producer tail. The
broader row-up panel use is reported through
`kls_first_last_row_supernode_panel_update`,
`kls_first_row_supernode_panel_update_run_count`,
`kls_first_last_row_supernode_panel_update_groups`, and
`kls_first_last_row_supernode_panel_update_rows`. The KLS-first separator pipeline
also publishes phase-local dense/common-tail panels for the stable private
prefix and consumes those panels before falling back to row-entry validation;
actual panel-backed use is reported through
`kls_first_last_separator_queue_pipeline_supernode_panel_update`,
`kls_first_separator_queue_pipeline_supernode_panel_update_run_count`,
`kls_first_last_separator_queue_pipeline_supernode_panel_update_groups`, and
`kls_first_last_separator_queue_pipeline_supernode_panel_update_rows`.
When a retained METIS separator tree covers a KLS-first block, the planner now
first tries a SubtreeLU Algorithm 6-style split before using the older retained
component queue: dominant subtrees are collapsed into pipeline roots, child
subtrees become private-thread candidates, and the candidate subtrees are
assigned by block-local row-input work. The factor path validates private
ownership against the original row dependencies before remapping; unsafe
partitions fall back to the legacy retained-component queue. Benchmark output
reports this through `kls_first_last_separator_queue_partitioned`,
`kls_first_separator_queue_partitioned_count`, and
`kls_first_last_separator_queue_split_components`.
When no retained separator queue applies and multiple threads are available,
KLS-first runs the same restartable Algorithm 5-style row pipeline over the
block's natural row order instead of serializing every row. That non-separator
executor is reported separately through `kls_first_last_row_pipeline`,
`kls_first_row_pipeline_run_count`, `kls_first_last_row_pipeline_rows`,
`kls_first_last_row_pipeline_threads`,
`kls_first_last_row_pipeline_partial`,
`kls_first_row_pipeline_partial_run_count`,
`kls_first_last_row_pipeline_partial_rows`, and
`kls_first_last_row_pipeline_partial_threads`. Generic row-pipeline dynamic
pivot epochs are reported through `kls_first_last_row_pipeline_pivot_tail`,
`kls_first_row_pipeline_pivot_tail_run_count`,
`kls_first_last_row_pipeline_pivot_tail_rows`,
`kls_first_last_row_pipeline_pivot_restarts`,
`kls_first_row_pipeline_pivot_restart_count`,
`kls_first_last_row_pipeline_pivot_serial_rows`,
`kls_first_last_row_pipeline_prefix_panel_rebuild`,
`kls_first_row_pipeline_prefix_panel_rebuild_count`, and
`kls_first_last_row_pipeline_prefix_panel_rebuild_rows`; these rebuild counters
remain in the ABI, but the current pivot path avoids full prefix rebuilds while
holding the ordered pipeline lock. It advances a post-pivot supernode validity
floor, so stale pre-pivot supernodes are rejected while newly published
post-pivot rows can form fresh Algorithm 4-style supernode runs; cached panels
use the same floor after the pivot column exchange. Set
`KLS_DISABLE_ROW_PIPELINE_PIVOT_SUPERNODE_REBASE=1` to restore the older
post-pivot supernode disable behavior for A/B traces. Separator pipeline epochs
continue to use the separator-prefixed counters below. KLS-first
panel-cache staging is reported through
`kls_first_row_panel_cache_build_count`,
`kls_first_row_panel_cache_build_panels`,
`kls_first_row_panel_cache_build_entries`,
`kls_first_row_panel_cache_append_count`,
`kls_first_row_panel_cache_append_panels`, and
`kls_first_row_panel_cache_append_entries`, so large runs can distinguish
full prefix-cache rebuilds from panels appended as rows are published.
Builds configured with `-DKLS_ENABLE_CBLAS_SUPERNODE=ON` can use the same
runtime `KLS_ENABLE_CBLAS_SUPERNODE=1` gate to consume eligible KLS-first
cached panels with CBLAS `dtrsv` and `dgemv`; otherwise the cached panel uses
the scalar in-panel solver. KLS-first CBLAS consumption now also requires at
least 2048 producer rows, at least 512 dense/tail update columns, at least 50M
estimated update operations, and at least 16 estimated operations per copied
panel entry, so medium and fragmented Level-2 panel updates stay on the
portable compact kernel. Dynamic column
exchanges rebuild the phase-local
pipeline cache from the post-exchange column order over the whole committed
prefix and reset the row-up producer panel caches. Separator pipeline
pivot-tail rows that are serialized after a restart use the same row-up
producer cache, so a restarted suffix can still publish and consume completed
row-supernode panels. Algorithm
1-style dynamic column exchanges are reported as
`kls_first_last_dynamic_column_pivots` and
`kls_first_dynamic_column_pivot_count`. When a retained METIS `NodeNDP`
separator map covers the factor order, the dynamic selector follows
SubtreeLU Algorithm 4 more closely: it computes `N'` as the last retained
factor row of the current collapsed component, compares the diagonal only with
the strongest candidate in `i+1..N'`, and rejects unsafe cross-domain pivot
needs instead of using a global outside-component maximum. Stats report exact
separator pivots through
`kls_first_last_separator_dynamic_column_pivots` and
`kls_first_separator_dynamic_column_pivot_count`, component-extent pivots
through `kls_first_last_separator_extent_dynamic_column_pivots` and
`kls_first_separator_extent_dynamic_column_pivot_count`, and cross-domain
rejections through `kls_first_last_separator_dynamic_column_fallbacks` and
`kls_first_separator_dynamic_column_fallback_count`. Private separator workers
also restrict off-diagonal pivot candidates to columns owned by the same
private domain. If that private search still rejects a pivot, KLS does not
publish the completed private prefix: later dynamic column exchanges can
conflict with already-published private U rows unless a future deferred-swap
owner is available, so the block falls back to the ordinary full row pipeline.
If a guarded separator
pipeline row still needs a dynamic pivot, KLS preserves the completed prefix,
publishes the scoped pivot row while holding the ordered pipeline lock, updates
the phase-local prefix snapshot, advances a column-order epoch, and lets
speculative suffix rows that started under an older epoch discard and restart
inside the same pipeline phase. Benchmark stats count these epoch recoveries
with
`kls_first_last_separator_queue_pipeline_pivot_restarts` and
separate them from older external serialized rows with
`kls_first_last_separator_queue_pipeline_pivot_serial_rows`. Partial
pre-updates and scalar row-supernode updates completed before a pivot restart
remain counted in the same separator-pipeline counters; phase-local cached panel
updates after the restart remain visible through the separator-pipeline panel
counters. Separator pivot epochs follow the same no-rebuild policy: stale
prefix metadata is disabled, valid cached panels are exchanged or deactivated,
and future panels can still be appended as rows publish. The retained prefix
rebuild counters are therefore expected to stay zero on this path:
`kls_first_last_separator_queue_pipeline_prefix_panel_rebuild`,
`kls_first_separator_queue_pipeline_prefix_panel_rebuild_count`, and
`kls_first_last_separator_queue_pipeline_prefix_panel_rebuild_rows`.
Active-rank pivot resets that rebuild row-supernode metadata, if re-enabled for
diagnostics, are exposed separately through
`kls_first_active_rank_pivot_reset_count`,
`kls_first_active_rank_pivot_reset_rows`,
`kls_first_active_rank_pivot_panel_rebuild_count`, and
`kls_first_active_rank_pivot_panel_rebuild_rows`. If that
row-up-looking bridge is not eligible
or still rejects a pivot, the pivoted KLS block tail
still runs; it reuses KLS's retained factor-order input map when available,
reported as
`kls_tail_last_mapped_columns` and `kls_tail_mapped_column_count`, so the
experimental first factor and accepted tail repairs can skip repeated original
CSC row remapping for in-block entries. It also builds KLS-owned row-major
`L`/`U` metadata for guarded row refactor/solve experiments, but it
first runs a cheap lower-bound work scan and skips the heavier row metadata
setup when that scan already proves the row plan cannot beat the exact EGraph
refactor work estimate. When the cheaper scan is inconclusive, KLS builds the
row pattern and copies numeric values into those mirrors only when their
retained row-work estimate is no larger than the exact EGraph estimate. Repeated
unchecked `kls_refactor` calls and checked fast-factor `kls_factor` calls then
try those row-major update paths automatically while the mirrors remain current;
successful checked pivot repairs reseed those mirrors so subsequent solves do
not have to fall back to published KLU column storage solely because the repaired
block rebuilt its LU payload. Benchmark stats expose
`row_refactor_total_group_work`, `row_refactor_auto_enabled`,
`row_refactor_auto_values_ready`, `row_refactor_auto_work_allowed`, and
`row_refactor_auto_should_run` so this automatic handoff can be audited against
the exact EGraph work model. The same stats also report
`row_refactor_auto_lower_bound_work`,
`row_refactor_auto_lower_bound_rejected`,
`row_refactor_auto_pattern_build_failed`, and
`row_refactor_auto_value_copy_failed` to distinguish the cheap work lower-bound
gate from later row-pattern/value preparation failures. Ordinary successful
numeric factors can also enter the same owned row/segment refactor preparation
when the NICSLU R1/R2 model recommends parallel numeric work and an exact
dependency schedule is already available; the existing row-work gate still
rejects row setup when it is more expensive than the EGraph estimate.
`kls_bench` and
`run_bench_suite.py` also accept
`--kls-first-factor env|off|on` so this path can be compared reproducibly
without relying on an ambient environment variable. Benchmark JSON reports the
model bridge as `row_refactor_auto_model_recommended`,
`row_refactor_auto_model_attempted`, and `row_refactor_auto_model_accepted`.
`KLS_ENABLE_ROW_SOLVE_FROM_NUMERIC=1` separately seeds the same row-major
solve mirrors after successful ordinary numeric factors/refactors when the
factor has no external KLS row/column permutation or scaling. This is an
experiment switch for measuring KLS-owned row solves independent of the
KLS-first scaffold; `kls_bench` and `run_bench_suite.py` expose it as
`--row-solve env|off|on`.
`factor_etree_block_start`,
`factor_etree_block_size`, `factor_etree_levels`,
`factor_etree_max_width`, `factor_etree_edges`,
`factor_etree_root_columns`, `factor_etree_leaf_columns`, and
`factor_etree_max_fanout` summarize the largest ordered diagonal block's
pivoted-factor ETree upper bound, matching the dependency concept used by the
CKTSO paper for pivoting tail work.

To exercise fast-factor pivot rejection on an unchanged MatrixMarket sparsity
pattern, scale diagonal entries only in the repeated numeric phase:

```sh
./build/kls_bench matrix.mtx --repeat 1 --refactor-repeat 0 --orientation auto \
  --stress-diagonal-scale 1e-12 --stress-diagonal-column 0 --json
```

The initial factorization still uses the original values; the repeated
factor/refactor/solve/residual path uses the stressed values. JSON includes
`stress_diagonal_scale`, `stress_diagonal_column`, and
`stress_diagonal_entries`. Use `--repeat 1 --refactor-repeat 0` when inspecting
the first rejected fast-factor tail, because later repeated calls can overwrite
the first rejection diagnostics.

After a fresh full factorization on medium and larger matrices, KLS sorts the
retained KLU numeric columns into ascending row order and detects strict
L-supernode runs: consecutive pivot columns whose sorted patterns nest, so
each run column is stored as an in-batch prefix plus a shared extended tail.
Repeated refactors then consume consecutive run producers per U column as one
dense in-panel solve plus a chunked shared-tail panel update that reads a
single index stream, in the serial mapped/BTF-pool kernels and in the EGraph
column kernels (including pipeline mode, where a batch is taken only when
every producer in it is already published). The panel update is only used
when the batched work amortizes the staging, and the hot loop is
multi-versioned so an AVX2/FMA clone is selected at load time on capable
hosts. Set `KLS_DISABLE_SNODE_PANEL_REFACTOR=1` to keep the numeric unsorted
and stay on the scalar per-producer updates; `KLS_TRACE_SNODE=1` reports run
coverage and batch counters.

Use `--threads N` to enable KLS-owned parallel work where it is currently
available. The first threaded path is repeated numeric refactorization across
independent BTF diagonal blocks for large, high-flop cases, including existing
KLU row-scaling modes when reusable scale storage and enough independent BTF
block work are present. KLS keeps a solver-owned worker pool for this path, so
repeated SPICE refactors do not relaunch threads every cycle; for matrices with
hundreds of thousands of small BTF blocks, workers claim short block ranges to
reduce scheduler mutex traffic. Small matrices, factor kernels, solve kernels,
and BTF cases that are not wide enough for the threaded path still use serial
execution. For unscaled serial refactors, KLS
also precomputes a scatter map from the fixed pivot order so repeated
single-block and serial BTF refactors do not redo the same `Q` and `Pinv`
structure lookups every cycle. For serial BTF refactors, this map also keeps
each column's off-block entries separate from diagonal-block entries so the
numeric refactor loop avoids reclassifying the same structure on every SPICE
step. The retained BTF map also stores each factor-order column's owning BTF
block, so block-aware EGraph refactors do not binary-search the BTF boundaries
for every column in fragmented many-block matrices. Dominant many-block BTF
refactors with bounded block count and input size can reuse the same map inside
the worker-pool path, including row-scaled factors, while extreme tiny-block
ASIC-style cases stay on the unmapped worker or serial path unless they enter
the guarded EGraph path. For threaded runs, KLS also records the exact no-pivot
refactor EGraph level schedule from the numeric U pattern; this metadata is
exposed in
benchmark output and is used by a guarded large unscaled single-block or
high-work dominant-BTF level-sliced refactor path when there is enough
dependency work and level width to offset thread and scratch overhead. KLS only
builds this schedule for single-block, dominant-block, or extremely fragmented
many-block shapes with enough numeric work to consume it; ordinary many-block
BTF cases skip the setup and stay on their existing refactor paths. EGraph
consumers reuse a solver-owned worker pool and dense scratch vectors across
repeated SPICE refactors, so the retained schedule no longer relaunches threads
or reallocates per-worker scratch on every accepted refactor step. The pipeline
dependency markers are also solver-owned generation-stamped atomics, so
pipeline passes do not allocate and clear an `n`-entry done array every cycle.
All-pipeline EGraph schedules skip the per-level thread-slice table entirely,
because they never enter the barriered cluster-mode loop that consumes it.
For huge single-block all-pipeline schedules, KLS also skips the `n`-entry
level-column list and lets workers claim columns in natural factor order; the
U-pattern dependency waits still enforce correctness. Dominant-BTF
all-pipeline schedules keep the level-column list, because level order provides
better structural balance across the large block and its fringe.
The huge single-block unscaled EGraph column kernel also relies on that
prevalidated map and schedule instead of repeating structural pointer,
row-bound, and U-order checks inside every column. The EGraph worker records
the selected column-kernel kind once before dispatching the solver-owned worker
pool, so large unscaled single-block, scaled single-block, and large-BTF
cluster/pipeline loops do not re-run the same kernel selection branch for
every column. The scaled single-block kernel also applies KLU row scales
directly while loading the fixed input-position map, avoiding the generic
BTF-capable value loader on large scaled all-pipeline refactors.
Most low-work dominant-BTF cases also stay on the mapped serial path, but those with
enough measured dependency work can consume the same EGraph path. When a
low-work dominant BTF decomposition has one 95%+ block plus thousands of tiny
fringe blocks below the EGraph size floor, KLS also keeps repeated refactors on
the serial mapped path because the BTF worker pool has too little useful
off-dominant work to amortize synchronization. In the scaled subset of that
class, KLS recomputes KLU row scales, applies the same fixed input-position map,
and permutes the scale vector back to pivot order after each refactor. When that
same many-fringe class has enough measured dependency work, KLS can instead run
the whole exact EGraph as a no-barrier pipeline: workers claim columns in
topological order and wait only on actual U-pattern predecessors. Barriered
cluster levels are split
by KLS's per-column no-pivot work estimate instead of equal column counts, and
above the CKTSO-style split level KLS can switch to a no-pivot pipeline tail
that waits only for actual U-pattern predecessors. Pipeline-tail columns are
claimed through an atomic work cursor instead of a fixed per-thread stride, so
threads that finish a ready tail column can continue with later tail work
instead of idling behind another thread's blocked dependency. This avoids thousands of
narrow-level barriers on high-work tails while preserving the existing
fixed-pivot LU storage. The same EGraph consumer can also run inside a
dominant BTF block with a large enough diagonal block, or a smaller dominant
block whose measured factor work is high enough to amortize the schedule. A
large-heavy dominant-block class can also use this path when the largest block is
at least 100k rows, covers at least 85% of the matrix, and has high actual
factor work. A medium-heavy class covers smaller dominant blocks when coverage is
between 85% and 95%, the largest block has at least 30k rows, and measured
factor/dependency work is high enough; this excludes the previously rejected
95% Rajat class while covering lower-work AT&T-style dominant blocks. The
schedule floor is lower for 95%+ dominant-block cases whose measured factor work
is modest but whose dependency graph still has enough independent work to consume
the threaded EGraph path; this covers low-work dominant-BTF circuit matrices
such as the IBM `dc*`, `trans*`, and `scircuit` cases without naming them in the
dispatch policy. Compact unscaled dominant-BTF matrices with a 10k-30k largest
block, bounded fringe-block count, and at least 20M measured factor flops can
also use the same exact EGraph path once the dependency work is high enough,
covering smaller CKTSO-gap cases without enabling the overhead-prone TSOPF
shapes. A bounded scaled dominant-BTF class with a 30k-60k largest block can
also use the EGraph path when measured factor work, LU fill, and dependency
work are high enough to amortize KLU row-scale recomputation. Moderate unscaled
single-block cases can also consume the EGraph refactor once their measured
factor work, LU fill, and dependency work clear lower single-block floors,
which covers matrices such as `rajat15` without using matrix-name tuning.
Larger single-block cases keep the higher general floors, which cover matrices
such as `HTC_336_4438`. For very high-work single-block factors, including
row-scaled factors whose scale vector can be recomputed and restored safely,
KLS can also run the full exact EGraph through the no-barrier pipeline path
instead of stopping after a narrow pipeline tail. This is limited to measured
factor/refactor work large enough to amortize full topological waiting and
keeps lower-work single-block cases on the barriered cluster/pipeline split.
Extremely fragmented unscaled BTF decompositions with hundreds of thousands of
tiny blocks and one substantial but non-dominant block can also use the exact
EGraph refactor once measured work is high enough, which prevents the large
block from remaining serial behind an otherwise wide BTF fringe. Those
fragmented-BTF schedules assign a small synthetic weight to singleton diagonal
blocks when balancing clustered EGraph levels, while dominant-BTF all-pipeline
and single-block schedules keep the original work model. Unscaled
single-block EGraph refactors use a slimmer
column kernel that bypasses BTF and scaling checks in the hot loop. The path
updates that block's local LU and off-block
entries while leaving ordinary many-block BTF cases on the existing worker-pool
path. For large cases that use KLU row scaling, the same path recomputes the row
scale factors, divides the mapped entries by the unpermuted row scales, and then
permutes `Rs` back to pivot order after the refactor. Small cases and ordinary
non-dominant BTF cases still use the existing mapped or BTF-worker paths.

Use `--orientation auto|normal|transpose` to control KLS's internal storage
orientation. `auto` uses the transposed pattern directly for small and medium
matrices where avoiding a second symbolic analysis is usually faster in a
repeated SPICE solve cycle. Large sparse, strongly diagonal circuit-like
patterns use normal storage directly to avoid a redundant transposed symbolic
analysis. Other larger matrices still compare normal and transposed symbolic
fill estimates before choosing.

Use `--scale auto|-1|0|1|2` to compare KLS/KLU row scaling modes when studying
pivoting and refactorization behavior. KLS defaults to `auto`, which can start
unscaled using KLU's `-1` no-scale/no-recheck mode for patterns already
validated by KLS, or sum-scaled when the diagonal is sparse but row magnitudes
are already balanced. Otherwise it starts from KLU's max row scaling mode (`2`).
Large low-degree, nearly diagonal circuit patterns, METIS-started medium
near-full bounded-degree or dense-diagonal patterns, and low-work dominant-BTF
patterns also start without row scaling to avoid repeated scale work once KLS
has validated the structure. Medium and large TSOPF-style spiked low-diagonal
METIS starts use sum scaling, while small spiked low-diagonal METIS starts use
KLU's unscaled `0` mode and skip the static-pivot trial that only adds setup
cost for that shape. Medium many-block, mostly diagonal, high-degree spike
patterns whose largest BTF block is large but not dominant also start in KLU's
unscaled `0` mode when max scaling consistently increases repeated refactor
work. Very fragmented full-rank BTF decompositions with hundreds of thousands
of blocks and one substantial but non-dominant block start in no-scale/no-recheck
mode so the exact EGraph refactor path remains available; KLU row-scaled modes
hide that path behind repeated scale work for this structural class. OPF-style
bounded-degree matrices whose numeric diagonal is about half missing or weak
start without row scaling. Large sparse-diagonal low-degree patterns can also
start with sum scaling when that avoids more expensive max-scaling behavior.
For large
expensive cases, `auto` can still try other numeric scaling modes when actual
flop/fill evidence justifies the extra work, but structural METIS starts that
are already known to need a specific scale mode skip redundant scale trials.
Large high-work METIS/no-BTF single-block paths that already selected max
scaling also skip those trial factorizations, preserving the chosen numeric
factors without paying for rejected scale candidates. Explicit numeric scale
values remain fixed.

Use `--pivot-tol T` to benchmark the diagonal pivot tolerance exposed by the
KLS API. The default is `0.001`, matching the underlying KLU default. With that
default, KLS can trial a lower `1e-4` tolerance on large high-fill matrices with
noticeable off-diagonal pivoting, keeping it only when fill and pivoting improve
without a large reciprocal-condition drop. TSOPF-style spiked low-diagonal
METIS starts that already match accepted dominant-BTF policies start directly
at `1e-4` to avoid paying for a rejected default-tolerance factorization.
Benchmark JSON reports both requested and selected pivot tolerance.

When `fast_factor` reuses an existing factor pattern, KLS also checks the
resulting L multipliers against the selected pivot tolerance. On unscaled BTF
patterns, a rejected reused pivot can restart the rejected diagonal block with
pivoting, splice the repaired block permutation back into the numeric object,
and retry the checked fast factorization when some columns may not have been
refreshed. If the failed pass had already refreshed all columns, or if a serial
BTF refactor reached the final block before rejection, KLS validates the
repaired block tail and skips the redundant checked retry. Scaled fast-factor
calls still use the pivot-checking refactor kernel so they can stop at the
first unsafe multiplier; KLS can repair a KLS-owned scaled rejected block,
recompute row scales, and continue with the threaded checked BTF pool, falling
back to the serial checked refactor, over only later BTF blocks when the
rejected pass left a valid prefix-current state.
Scaled prefix-current rejects whose validated repair preserves a non-empty
live prefix state can also execute the conservative serial suffix restart of
the rejected block with pivoting; accepted scaled and unscaled serial suffix
restarts refresh only the off-diagonal column suffix whose inverse row
permutation can change. Scaled all-refresh KLU refactor rejects first
recompute the row scale vector back to input-row order, so the same block
repair and validated non-root serial suffix restart can run before the accepted
factor is restored to KLU's pivot-order row-scale convention.
When a BTF refactor rejects a non-final block, the repaired block can now be
followed by a checked continuation over only the later BTF blocks instead of
restarting from the first block; unscaled and KLU row-scaled states both try the
threaded BTF worker pool before serial continuation. The threaded BTF worker
pool also marks completed diagonal blocks and reports the same prefix-current
state only when every earlier block finished before a checked pivot reject. For
unscaled or scaled prefix-current/all-current rejects whose validated repair
preserves a non-empty live prefix state, KLS can execute the same conservative
serial suffix restart before falling back to full block repair. This executable tail path is
not limited to the checked row-refactor candidate diagnostic; the diagnostic is
retained only to explain row-major candidate quality when that metadata exists.
Root-of-block rejects can use the same KLS-owned pivoted block kernel. When a
prefix-current ETree closure leaves independent single-block columns outside the
tail, KLS first refreshes those preserved columns with the mapped no-pivot
column kernel, then runs the shorter pivoted tail envelope and counts it as a
tail restart. Internal envelope gaps that can be copied unchanged are skipped;
dependency-blocked internal gaps are promoted into the active tail and
recomputed with pivoting, while structural gap-copy failures still fall back.
Full-suffix root rejects still remain ordinary KLS block restarts. For any
retained ETree mask with more than one active tail row, including contiguous
unfinished tails, the KLS-owned repair now factors the boundary pivot row first
and then runs the remaining active descendant rows through the retained row
pipeline with an active ETree-rank map. When that retained mask covers the
whole local repair envelope, the preserved-column refresh is treated as a
successful no-op instead of rejecting the ETree-tail plan.
Those descendants use the same prefactorization/postfactorization loop as the
row-up-looking pipeline: finished dependencies are applied while earlier active
tail rows are still completing, then skipped dependencies are finished before
the row is published. If a descendant row still needs a dynamic column pivot,
the ETree-prefactor phase now serializes just that pivot row, refreshes the
row-up producer panel cache for the new column order, and resumes the remaining
active descendants in the same retained ETree-rank order. That guarded executor
reports
`fast_kls_block_restart_last_row_pipeline_etree_ready`,
`fast_kls_block_restart_row_pipeline_etree_ready_count`,
`fast_kls_block_restart_last_row_pipeline_etree_ready_rows`,
`fast_kls_block_restart_last_row_pipeline_etree_ready_threads`,
`fast_kls_block_restart_last_row_pipeline_etree_prefactor`,
`fast_kls_block_restart_row_pipeline_etree_prefactor_count`,
`fast_kls_block_restart_last_row_pipeline_etree_prefactor_rows`,
`fast_kls_block_restart_last_row_pipeline_etree_prefactor_threads`,
`fast_kls_block_restart_last_row_pipeline_etree_prefactor_wait_rows`, and
`fast_kls_block_restart_last_row_pipeline_etree_prefactor_wait_deps`, while
the existing pivot-epoch counters report the serialized descendant rows.
For rejected blocks fully covered by a retained separator tree, KLS now tries
the same SubtreeLU Algorithm 6-style private/pipeline split used by the
KLS-first path before installing the older ETree-tail active mask. The repair
planner remaps the block-local row/column order into private rows followed by
collapsed separator-pipeline rows, validates private ownership after that
remap, runs proven-private rows concurrently, and finishes the separator roots
with the restartable pivot-capable row pipeline. If private ownership or the
private phase cannot be proven safe, the block stays on the ordinary full
row-pipeline repair path in the separator order. Accepted separator-queue
repairs report
`fast_kls_block_restart_last_row_pipeline_separator_queue`,
`fast_kls_block_restart_row_pipeline_separator_queue_count`,
`fast_kls_block_restart_last_row_pipeline_separator_private_rows`,
`fast_kls_block_restart_last_row_pipeline_separator_pipeline_rows`,
`fast_kls_block_restart_last_row_pipeline_separator_private_threads`,
`fast_kls_block_restart_last_row_pipeline_separator_partitioned`, and
`fast_kls_block_restart_last_row_pipeline_separator_split_components`.
Suffix-exact retained ETree tails inside separator-covered blocks can also use
the ordered pivot-capable row pipeline after validating the separator pivot
scope and passing pre-commit numeric and rowwise-U pivot checks. Non-suffix
separator-covered tails stay on the full separator queue path because they can
leave later weak pivots outside the active ETree tail. This remains a guarded
CKTSO Algorithm 5 prefactor/postfactorization step plus a SubtreeLU separator
full-block repair, rather than the full production ETree-descendant scheduler.
The experimental KLS-owned row/segment checked fast-factor path is explicit
opt-in. Set `KLS_ENABLE_CHECKED_ROW_REFACTOR=1` to run it as the first checked
fast-factor attempt and to let it use the parallel row scheduler when multiple
threads are available. Leaving the variable unset keeps the production checked
fast-factor path on the column/EGraph repair ladder. A rejected
parallel row pass records the unsafe dependency pivot and falls back through
the same block-repair path. Checked KLS-owned row passes test the guessed
diagonal against the maximum absolute value in the current U row before
publishing that row, matching CKTSO's row-wise pivot criterion in addition to
the existing multiplier checks. The parallel pass now reuses the retained
completion bitmap to report a prefix-current refresh state only when all rows
before the rejected pivot have finished. When dense row segments defer writes
in KLS-owned row-major mirrors, KLS publishes only that proven prefix back into
the KLU numeric object before reporting prefix-current; otherwise it keeps the
conservative unknown refresh state because rows in the active level may have
completed out of factor-order prefix. KLU-compatible checked fast refactors
now also validate finalized U rows against the same CKTSO row-wise rule before
accepting the reused pivot order; those delayed rejects are reported as
all-current because the column refactor has already refreshed the factor.
The experimental row scheduler also
splits precomputed row-group levels into barriered cluster levels and a
dynamic topological pipeline tail. Row groups carry a retained execution kind
so single-row, generic multi-row, and dense multi-row groups dispatch without
rediscovering that shape during each numeric pass. Multi-row groups also retain
the shared trailing `U` slice offset used by the dense row-segment kernels, so
the executor consumes the compact row-segment descriptor instead of
re-deriving that slice from adjacent rows. Unchecked dense multi-row groups
whose internal dense and shared-trailing update work is large enough first use a
worker-local compact row-major panel before falling back to the direct
row-mirror kernel. The compact-panel gate requires both enough total structural
work and enough work per copied panel entry, so low-arithmetic-intensity dense
groups stay on the direct row-mirror kernel. Stats report
`row_refactor_last_compact_dense_panel` and
`row_refactor_compact_dense_panel_count` for executed compact panels, plus
`row_refactor_compact_dense_panel_eligible_count`,
`row_refactor_compact_dense_panel_eligible_rows`,
`row_refactor_compact_dense_panel_update_work`, and
`row_refactor_compact_dense_panel_entries` for the structural opportunity that
survives the compact-panel gate. Eligible dense groups also allocate
solver-owned compact panel slices reported as
`row_refactor_compact_dense_panel_persistent_groups` and
`row_refactor_compact_dense_panel_persistent_entries`; the compact kernel uses
those slices before falling back to worker-local scratch and reports actual
retained-slice execution as
`row_refactor_last_compact_dense_panel_persistent` and
`row_refactor_compact_dense_panel_persistent_run_count`. Stats also report
row-group shape counters such as
`row_refactor_group_batch_max_width`,
`row_refactor_group_batch_width_le_4_count`,
`row_refactor_group_batch_width_le_8_count`,
`row_refactor_group_scalar_candidate_count`,
`row_refactor_group_scalar_candidate_rows`,
`row_refactor_group_scalar_short_count`,
`row_refactor_group_scalar_short_rows`,
`row_refactor_group_scalar_stop_level_mismatch_count`,
`row_refactor_group_scalar_stop_internal_dep_count`,
`row_refactor_group_scalar_stop_next_segment_count`,
`row_refactor_group_scalar_stop_max_width_count`,
`row_refactor_group_scalar_stop_matrix_end_count`,
`row_refactor_group_generic_rows`,
`row_refactor_group_generic_max_width`,
`row_refactor_group_dense_rows`,
`row_refactor_group_dense_max_width`, and
`row_refactor_group_*_work` to distinguish small independent-batch overhead
from dense/generic panel work concentration. Later rows can consume
a current-pass retained compact panel as a supernode update source when they
have a contiguous suffix of dependencies on that completed dense group; stats report
`row_refactor_last_compact_supernode_update`,
`row_refactor_compact_supernode_update_count`,
`row_refactor_compact_supernode_update_rows`, and
`row_refactor_compact_supernode_update_entries`. The partial supernode pipeline
now also publishes compact producer-panel prefixes row by row, while the
ready-queue scheduler releases successor groups only after their exact
row-predecessor count reaches zero. A consumer with a contiguous dependency
suffix ending at a valid producer prefix can apply the compact supernode update
without rebuilding that panel; stats report
`row_refactor_last_compact_supernode_partial_update`,
`row_refactor_compact_supernode_partial_update_count`,
`row_refactor_compact_supernode_partial_update_rows`, and
`row_refactor_compact_supernode_partial_update_entries`. When the producer supernode
has a shared trailing panel, KLS now accumulates that contribution in contiguous
worker scratch and scatters it once, matching the matrix-vector update half of
the triangular-solve plus update shape described in the SubtreeLU paper without
adding an external BLAS dependency. The default path accumulates the trailing
vector while each producer multiplier is already live, so it avoids a second
suffix pass before the final scatter; stats report
`row_refactor_last_compact_supernode_gemv`,
`row_refactor_compact_supernode_gemv_count`,
`row_refactor_compact_supernode_gemv_rows`, and
`row_refactor_compact_supernode_gemv_entries`. Unset
`KLS_ENABLE_COMPACT_SUPERNODE_TRSV` keeps the triangular part on the direct
scalar dependency walk because current CKTSO-gap probes show the contiguous
worker-scratch triangular solve over-stages the row-refactor scaffold even after
large-run gates. Setting the variable to `1` forces that paper-shaped
triangular-solve probe for coverage and experiments. It reports
`row_refactor_last_compact_supernode_trsv`,
`row_refactor_compact_supernode_trsv_count`,
`row_refactor_compact_supernode_trsv_rows`, and
`row_refactor_compact_supernode_trsv_entries`. KLS-owned scalar
multi-producer row-panel updates are available as an opt-in
`KLS_ENABLE_MULTI_PRODUCER_SUPERNODE=1` experiment when the retained row/segment
structure and work gates accept them. The default stays on the scalar/compact
fallback because current CKTSO-gap probes show the fragmented batch scaffold can
over-stage these rows. These cover exact common producer suffixes, ragged
single-producer suffixes, contiguous independent-row producer suffixes, and
fragmented dense-consumer external prefixes without requiring CBLAS. The
independent-row multi-producer executor now keeps its metadata, pivots, U
scratch, and panel descriptors in per-worker scratch instead of allocating them
for each accepted batch. Dense-consumer exact, ragged, and fragmented producer
batches use the same worker scratch model for their run descriptors and index
metadata, while retaining a separate byte workspace for fallback target maps.
This is a storage cleanup for the paper-shaped executor, not a default policy
change. Current top-five
CKTSO-gap forced-row probes still accept zero independent compact-supernode
batches, so that part is covered by targeted smoke fixtures rather than by
those slow-case rows.
Builds
configured with `-DKLS_ENABLE_CBLAS_SUPERNODE=ON` also compile an opt-in
`KLS_ENABLE_CBLAS_SUPERNODE=1` supernode experiment that uses standard CBLAS
calls with the same scalar fallback and pivot checks. Completed producer
supernodes can update later rows with CBLAS `dtrsv` plus `dgemv`, matching the
paper's direct update shape over the retained row-major panel when the
structural update work and row/panel dimensions are large enough to amortize
BLAS calls. Smaller producer/consumer shapes stay on the KLS-owned scalar
compact kernels, and checked dense panels keep the native blocked kernel unless
at least one row update can pass the same large-work CBLAS gate. When a whole
unchecked dense consumer group, or a contiguous row subrange inside it, has the
same ordered list of completed dense producer suffixes as its external
dependency pattern, the CBLAS experiment can batch those producer updates
across those consumer rows with one `dtrsm` and one `dgemm` per producer.
Earlier producer updates can flow into later producer multiplier columns before
those later suffixes are solved; stats report
`row_refactor_last_compact_supernode_batch`,
`row_refactor_compact_supernode_batch_count`,
`row_refactor_compact_supernode_batch_rows`,
`row_refactor_compact_supernode_batch_dep_rows`, and
`row_refactor_compact_supernode_batch_entries`. The related
`row_refactor_compact_supernode_batch_pattern_*`,
`row_refactor_compact_supernode_batch_candidate_*`, and
`row_refactor_compact_supernode_batch_rejected_work_count` counters distinguish
missing same-pattern row subranges from candidates rejected by the structural
work gate. Dense row groups can pass through a named native row-panel selector
before falling back to the direct row-mirror kernel. Unset
`KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR` keeps the structural default selected by
row-refactor analysis, `0` forces the scalar row-major dense-group fallback,
`auto` uses the structural work gate, and `1` forces the retained compact panel
path for eligible dense groups. The selected path stores the
group in solver-owned row-major panel slices when available, applies the
portable blocked panel factor/update kernel, and preserves the checked
row-ordered update/check/publish loop so a rejected pivot leaves the same
prefix-visible row-major state as the direct row-mirror kernel. Stats report
`row_refactor_native_row_panel_enabled`,
`row_refactor_last_native_row_panel`,
`row_refactor_native_row_panel_count`,
`row_refactor_native_row_panel_rows`,
`row_refactor_native_row_panel_entries`,
`row_refactor_native_row_panel_blocked_count`,
`row_refactor_native_row_panel_blocked_rows`,
`row_refactor_native_row_panel_blocked_entries`,
`row_refactor_native_row_panel_fallback_count`, and
`row_refactor_native_row_panel_checked_reject_count` so CKTSO-gap runs can
distinguish native-panel use from unsupported groups and checked-pivot exits.
Row segments now also direct-load raw input entries
into the retained row-major segment mirrors when the `L` row, pivot,
in-segment `U`, and shared trailing `U` pattern can represent the row exactly.
The row-pattern builder retains a per-input-entry segment target map for dense
and sparse segment rows, so repeated refactors write raw values through
precomputed external, `L`, pivot, and `U` destinations instead of rediscovering
the slots in the numeric loop. The same retained target map now covers
single-row generic groups, compact dense panel rows, and independent batch
groups as well: scalar dependencies still enter the work vector for the
existing up-looking update, batch dependencies load into the batch dependency
vector, raw pivots and row-major `U` tails load directly into retained row
storage, and compact panel loads place `L`, pivot, dense-`U`, and trailing-`U`
entries through the retained destinations instead of searching the panel
columns again. Exact, ragged, and fragmented compact-supernode batch consumers
use the same retained compact-panel loader before applying producer updates, so
their `L`, pivot, dense-`U`, and trailing-`U` inputs no longer stage through
the sparse work vector just to populate the row-major panel; if a residual
work-vector value remains for the same row-major slot, it is merged into the
retained panel before publication. Stats report
`row_refactor_last_sparse_segment_direct_input_rows`,
`row_refactor_sparse_segment_direct_input_rows`,
`row_refactor_last_batch_direct_input_rows`,
`row_refactor_batch_direct_input_rows`,
`row_refactor_segment_input_target_rows`,
`row_refactor_segment_input_target_entries`,
`row_refactor_last_segment_target_input_rows`, and
`row_refactor_segment_target_input_rows`. For generic-only row patterns,
pipeline groups also
retain their external dependency rows so the scheduler can wait on them once
before running the group kernel; stats report
`row_refactor_group_cluster_levels`, `row_refactor_group_pipeline_groups`,
`row_refactor_group_pipeline_rows`, and
`row_refactor_group_pipeline_work`. The same row-group metadata now also retains
the reverse group graph and reports `row_refactor_group_dependency_edges`,
`row_refactor_group_root_count`, `row_refactor_group_leaf_count`, and
`row_refactor_group_max_fanout`, which are the row-segment task-graph counters
needed by future private/pipeline partitioning and tail-restart schedulers.
The column EGraph refactor path also has an opt-in paper probe behind
`KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=1`. It uses retained consecutive
supernode-candidate ranges to solve a dependency run through a dense internal
triangular panel and, only when the producer L columns share the same trailing
row list, accumulates that trailing contribution in worker scratch before one
scatter. Stats report `refactor_last_supernode_update_runs`,
`refactor_last_supernode_update_rows`, `refactor_last_supernode_update_entries`,
and cumulative `refactor_supernode_update_*` totals; CBLAS builds additionally
report the BLAS-taken subset as `refactor_last_supernode_cblas_update_runs`,
`refactor_last_supernode_cblas_update_rows`,
`refactor_last_supernode_cblas_update_entries`, and cumulative
`refactor_supernode_cblas_update_*` totals. Non-CBLAS or CBLAS-disabled runs
use a portable blocked cached-panel update and report that subset as
`refactor_last_supernode_blocked_update_runs`,
`refactor_last_supernode_blocked_update_rows`,
`refactor_last_supernode_blocked_update_entries`, and cumulative
`refactor_supernode_blocked_update_*` totals. This path is intentionally off
by default. It now builds persistent producer-side compact panels for
eligible retained EGraph supernodes and publishes panel rows as producer
columns finish, so later consumers reuse the dense/internal and shared trailing
values instead of reconstructing that structure for every dependency run. The
same cache can consume a published prefix or suffix of a retained panel, so a
column inside a supernode no longer has to rebuild a temporary partial panel
just to use already completed producer columns. The cached-panel consumer now
also runs from the generic mapped EGraph refactor kernel, so scaled BTF and
smaller BTF states can use the same retained producer panels instead of staying
on the scalar dependency loop. CBLAS builds can also apply an eligible cached
EGraph panel with a unit-diagonal `dtrsv` over the retained internal panel plus
`dgemv` updates for the dense suffix and shared trailing rows. Focused runs
show this removes the worst rebuild overhead, but the path is still slower than
default KLS because the broad default path still needs coarser batched panel
kernels.
When dependency waiting is active, the cached-panel consumer follows the
SubtreeLU large-supernode split rule: if a retained producer run has at least
`2 * threads` rows, it may consume the completed prefix ending `threads` rows
before the producer supernode tail and leave the tail for the next dependency
iteration. Set `KLS_ENABLE_EGRAPH_SUPERNODE_SPLIT=0` to disable this split for
A/B runs.
The retained consumer-plan diagnostics also report group-L batch candidate
payoff counters:
`refactor_supernode_consumer_plan_group_l_batch_candidate_advance_*` measures
the scalar prefix work needed to reach an exact retained producer panel, and
`refactor_supernode_consumer_plan_group_l_batch_candidate_payoff_*` measures
the subset where the grouped panel update work is at least that advance work.
These are diagnostics for future grouped producer-panel execution; the default
path does not build this cache unless the retained group-L cache experiment is
enabled. The opt-in group-L batch executor only maps payoff-positive groups, so
the batch path follows the same advance-work test instead of executing every
reusable shape candidate. If one clean numeric pass covers only a small
fraction of the retained group-L rows, later passes disable the batch executor
and keep the cache as diagnostic evidence.
The exact-shape bounded-owner diagnostic uses a default pre-panel advance
limit of 128 U dependencies and reports the active value through
`refactor_supernode_consumer_plan_shape_bounded_advance_dep_limit`. Set
`KLS_REFACTOR_PLAN_GROUP_L_BOUNDED_ADVANCE_MAX_DEPS=<n>` to rerun the same
retained-plan accounting with a different bound; this changes diagnostic
counters only and does not enable a numeric executor.
`KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_SHAPE_TARGETS=1` is an
off-by-default retained-target substrate for the grouped output accumulator
path. It builds the retained consumer plan, then stores a compact block-local
target-row map for each reusable producer shape while separately counting the
raw per-current publish entries. Benchmark JSON reports
`refactor_supernode_consumer_plan_shape_targets_built`,
`refactor_supernode_consumer_plan_shape_target_group_count`,
`refactor_supernode_consumer_plan_shape_target_run_count`,
`refactor_supernode_consumer_plan_shape_target_entries`,
`refactor_supernode_consumer_plan_shape_target_rows`, and
`refactor_supernode_consumer_plan_shape_target_max_rows`. This does not change
default numeric execution; it retains the sparse publish surface needed by a
future CKTSO/SubtreeLU-style multi-current producer executor.
Set `KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=cached` to isolate only the durable
cached-panel consumer and skip the per-consumer temporary panel reconstruction
fallback used by the full `=1` experiment.
Cached runs skip the heavier panel consumer unless the producer column belongs
to a retained panel, and JSON reports
`refactor_supernode_cached_probe_*` counters so profiling can distinguish panel
misses, contiguous producer runs, work-gate acceptance, and applied cached
updates.
In cached-only mode, if one completed numeric pass probes retained panels but
finds no work-gate-accepted cached updates, later passes with the same retained
panel cache skip the cached probe; JSON reports
`refactor_supernode_cached_probe_disabled` and
`refactor_supernode_cached_probe_disable_count`.
In the full `=1` experiment, if a clean numeric pass probes retained panels but
accepts less than one amortizable supernode-update window, later passes skip the
full supernode probe and report `refactor_supernode_update_disabled` plus
`refactor_supernode_update_disable_count`. Productive cached-panel work remains
eligible; this only trims the low-work full-probe case.
For paper-gap diagnosis without changing execution, set
`KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_STATS=1`. Schedule construction then
counts actual U-stream dependency runs that fall inside retained EGraph
supernodes and reports `refactor_supernode_consumer_run_count`,
`refactor_supernode_consumer_run_rows`,
`refactor_supernode_consumer_run_max_width`,
`refactor_supernode_consumer_suffix_count`,
`refactor_supernode_consumer_l_entries`, and
`refactor_supernode_consumer_internal_entries`. These counters estimate how
much scalar dependency work a future row-major/supernodal numeric object could
consume before the current common-trailing cached-panel gates are applied.
Set `KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN=1` to retain those contiguous
consumer runs as a producer-addressable plan. The plan stores each run's
current column, producer dependency, run length, producer-panel start, and
producer-panel offset, and reports the retained shape through
`refactor_supernode_consumer_plan_*` counters. Set
`KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_EXEC=1` to let the experimental
cached-panel executor consume that plan. This executor is intentionally
off-by-default: focused CKTSO-gap runs show that the current completed-panel
cache covers only a small fraction of retained rows, so the plan is primarily
staging for a future producer/consumer row-major numeric task. When this flag
is the only supernode update gate, KLS restricts cached-panel probing to
retained plan hits; explicit `KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES` modes keep
their broader opportunistic probes.
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_EXEC=1` is a stricter
Algorithm 5 probe: it builds the retained plan, selects payoff-positive
producer-prefix subsets, materializes those prefixes as the existing ragged
U-supernode L pattern, and runs the ragged-L executor only for selected runs.
For those selected runs, KLS trusts the retained Algorithm 5 payoff selector
instead of reapplying the generic ragged-U work gate; this makes the opt-in
path a direct scalar replay test for the paper selector, not a candidate
default.
That executor now addresses its scratch data through the retained
group/current workspace map, so the opt-in numeric path exercises the same
compact current-slot layout that a true multi-current batch kernel will need.
This is also intentionally experimental and off by default; it tests whether
the Algorithm 5 payoff surface is enough without a true multi-current batch
executor. Focused CKTSO-gap probes show that scalar selected-run replay is not
enough: making the selected ASIC runs execute regressed the top-five and
top-ten payoff-exec controls, so the remaining paper gap is still the grouped
multi-current numeric executor.
When `KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_PLAN=1` or the stricter
exec flag is set, KLS now also retains the selected producer groups as a
multi-current descriptor for that future executor. The retained internal layout
maps each selected run to its current-column slot and gives that slot a compact
workspace offset, so a batch kernel can gather/advance/publish current
workspaces without rediscovering the group shape. The descriptor also stores a
CSR run list for each retained current slot, so a group-triggered executor can
walk the selected runs owned by a current column without falling back to the
per-column consumer-plan lookup. Benchmark JSON reports
`refactor_supernode_algorithm5_payoff_group_current_total`,
`refactor_supernode_algorithm5_payoff_group_multi_current_count`,
`refactor_supernode_algorithm5_payoff_group_max_currents`,
`refactor_supernode_algorithm5_payoff_group_workspace_rows`,
`refactor_supernode_algorithm5_payoff_group_max_workspace_rows`,
`refactor_supernode_algorithm5_payoff_group_advance_deps`,
`refactor_supernode_algorithm5_payoff_group_max_advance_deps`, and
`refactor_supernode_algorithm5_payoff_group_zero_advance_runs`. It also reports
`refactor_supernode_algorithm5_payoff_group_target_entries`,
`refactor_supernode_algorithm5_payoff_group_max_target_entries`, and
`refactor_supernode_algorithm5_payoff_group_max_run_target_entries`, which
count the dense suffix and L-trailing update surface that a retained
multi-current accumulator would have to apply after the selected producer
prefix. The retained descriptor also stores each selected run's target-entry
offset, each run's exact target-slot offset, and each group's target-slot span;
the slot counters
`refactor_supernode_algorithm5_payoff_group_target_slots`,
`refactor_supernode_algorithm5_payoff_group_max_target_slots`, and
`refactor_supernode_algorithm5_payoff_group_max_run_target_slots` count the
deduplicated accumulator/publish positions retained in `group_target_cols`.
Benchmark output reports the pattern-width sum and maximum as
`refactor_supernode_algorithm5_payoff_group_pattern_width` and
`refactor_supernode_algorithm5_payoff_group_max_pattern_width`. It also reports
`refactor_supernode_algorithm5_payoff_current_state_rows` and
`refactor_supernode_algorithm5_payoff_current_state_max_rows` for the retained
compact local-row map of each current slot. That map is the union of input
scatter rows, current-column U/L rows, the pivot row, retained advance rows,
prefix rows, and retained target rows. These counters estimate the number of
distinct current-column workspaces, compact prefix workspace rows, addressable
target accumulator slots, retained update-entry work, prefix-advance
dependencies, and persistent current-state rows a real Algorithm 5 batch would
need. On the slow ASIC diagnostics, zero-advance selected runs are rare or
absent, so the next paper-aligned executor has to batch prefix advancement
before the shared producer update rather than relying on a first-dependency
shortcut or BLAS threshold tuning. These counters do not enable the old
per-current scalar replay by default.
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_PREFIX_PREP=1` is an
off-by-default retained-prefix replay probe. When a selected producer prefix is
published, KLS can prepare each eligible current slot by applying already
published advance dependencies into a private scratch row, retaining the prefix
U coefficients in the current workspace, retaining dense-tail/L-trailing
target deltas in the target slots, and later letting the ordinary scalar
ragged-L column path consume that ready prefix instead of recomputing it. This
is intentionally not a claimed grouped executor: the advance/scatter phase is
still per current column. The probe is useful because it validates the
producer-prefix/current-slot handoff that a CKTSO-style grouped executor needs,
but focused ASIC runs still show that batching several current workspaces
through the same producer panel is the missing performance step. Benchmark JSON
reports consumed retained prefixes through
`refactor_last_supernode_algorithm5_payoff_prefix_prep_runs`,
`refactor_last_supernode_algorithm5_payoff_prefix_prep_rows`,
`refactor_last_supernode_algorithm5_payoff_prefix_prep_target_entries`, and
`refactor_last_supernode_algorithm5_payoff_prefix_prep_target_slots`.
The scalar slot-accum and prefix-prep probes write retained dense suffix target
slots directly and only build row-stamp maps for irregular L-trailing targets.
This matches the retained Algorithm 5 layout more closely, but it does not
promote those probes: current top-five CKTSO-gap checks still reject per-current
scalar replay versus the default path.
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_PREFIX_PREP=1` is a
stronger grouped retained-prefix probe. It claims all currently eligible
current slots for a producer-prefix trigger, prepares their prefix workspaces,
and applies the retained dense and trailing target slots as one grouped pass
before the ordinary ragged-L path consumes the ready prefixes. The flag also
requests the retained Algorithm 5 payoff plan, runtime workspace, target slots,
and current-state flags; it remains opt-in and off by default. Current top-five
CKTSO-gap checks reject this retained-target executor as a default: it proves
the grouped current-slot handoff works, but the large ASIC rows materialize tens
of millions of retained target updates and are slower than the default EGraph
path. This points the next paper-aligned work at direct row-workspace
supernode-prefix/suffix application, matching CKTSO/SubtreeLU Algorithm 5
semantics, rather than BLAS threshold tuning.
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_PREP=1` is the
targetless grouped prefix probe. It uses the same selected Algorithm 5
producer-prefix groups and current-state flags, but the grouped preparation
pass retains only prefix workspaces and applies only the internal prefix
triangular work. When the ordinary ragged-L path later consumes a ready slot,
it streams the dense suffix and L-trailing contribution directly into the
current column workspace instead of replaying retained target slots. Benchmark
JSON still reports the skipped suffix/trailing work as
`refactor_last_supernode_algorithm5_payoff_prefix_prep_target_entries`, while
`refactor_last_supernode_algorithm5_payoff_prefix_prep_target_slots` is zero in
this direct mode. The flag is opt-in: focused CKTSO-gap checks show it removes
the retained target workspace and improves the retained-target grouped probe,
but it is still slower than the default EGraph path.
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_ADVANCE_SEED=1`
adds an additional opt-in probe on top of direct-prefix prep. The grouped
producer records the post-advance row workspace slots and skipped U
coefficients using the existing Algorithm 5 advance descriptors; the consuming
column can then restore that state and jump over any remaining advance
dependencies before consuming the prepared prefix. Benchmark JSON reports this
through `refactor_last_supernode_algorithm5_payoff_advance_seed_runs`,
`refactor_last_supernode_algorithm5_payoff_advance_seed_deps`, and
`refactor_last_supernode_algorithm5_payoff_advance_seed_slots`. This validates
the paper-level missing piece that direct-prefix prep was otherwise replaying,
but it remains rejected as a default: the focused ASIC cases skip many advance
dependencies, yet restoring hundreds of thousands of row slots is slower than
the targetless direct-prefix probe and the default EGraph path.
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_CURRENT_STATE=1`
is the fuller retained-current-row probe. It also implies direct-prefix prep,
stores values for the retained per-current state row map, and restores that
state plus skipped U coefficients before consuming the prepared prefix.
Benchmark JSON reports the restore through
`refactor_last_supernode_algorithm5_payoff_current_state_seed_runs`,
`refactor_last_supernode_algorithm5_payoff_current_state_seed_deps`, and
`refactor_last_supernode_algorithm5_payoff_current_state_seed_rows`. The
top-five CKTSO-gap check rejects this as a default too: it fills the
CKTSO/SubtreeLU current-row persistence gap, but the hard ASIC cases restore
hundreds of thousands to millions of state rows per run, making it slower than
the lighter advance-seed probe.
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_FINAL_STATE=1`
is the final retained-current-state probe. It also implies direct-prefix prep and
current-state storage. The grouped producer applies the prefix's dense suffix and
L-trailing effects into the retained current-state rows, so the consuming column
restores that final row state, writes the skipped advance and prefix U
coefficients, and jumps past the prepared prefix. Benchmark JSON reports this
through the same `refactor_last_supernode_algorithm5_payoff_current_state_*`
counters plus direct-prefix rows in
`refactor_last_supernode_algorithm5_payoff_prefix_prep_*`. Focused CKTSO-gap
checks reject this as a default too: it removes the later prefix replay, but the
large ASIC cases still copy millions of retained state rows.
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_SPARSE_RESTORE=1`
is a guarded variant of final retained-state mode. It also implies final-state
storage, but the consuming column restores only the retained rows that the
remaining scalar update can read or modify: remaining U dependencies, their L
rows, the pivot row, and the final L column rows. Input-scatter rows outside that
set are cleared before the scalar loop continues. A cheap L-pattern upper bound
skips sparse row-set construction unless the continuation can avoid restoring a
material share of the retained rows. The same
`refactor_last_supernode_algorithm5_payoff_current_state_seed_rows` counter then
reports restored rows rather than total retained rows, making the retained-copy
payoff directly measurable without changing the default executor.
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_COMPLETE=1` is
a safe direct-completion probe for the same final retained state. It also
implies direct-prefix prep/current-state/final-state mode. For an eligible
BTF-local column, KLS publishes terminal columns directly from the retained
final state when the prepared prefix covers all U dependencies; otherwise it
computes the remaining U stream, pivot, and L column from a private
retained-state copy and only commits factor storage after the pivot is nonzero.
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_SPARSE_DELTA=1`
is the same direct-complete setup with a sparse retained-state delta accumulator
instead of the private copy. The direct-complete suffix is now single-pass:
dependency readiness and retained-row mapping are checked while computing into
private/local state, so the suffix is not scanned once for validation and again
for arithmetic. Focused CKTSO-gap checks solve correctly but reject both as
defaults: single-pass terminal/private direct-complete measured `2.3314s`
SPICE-cycle geomean on the five-row slice, while sparse delta removed
retained-row copy but still measured `2.4101s`. The current default on the same
source measured `1.4629s`, so the next paper-level gap is a grouped live
Algorithm 5 workspace that advances the remaining dependency stream for
multiple current columns, not BLAS thresholding or retained-row copy alone.
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_COMPLETE=1` extends the
direct-complete probe with retained final-dependency triggers. The grouped
producer-prefix path still prepares the current-state workspaces in batch, but
KLS also builds a map from each current slot to the last U dependency in that
column. When that final dependency publishes, a worker can claim the prepared
current column and dispatch it through the existing direct-prefix completion
path instead of waiting for the ordinary pipeline position. This directly tests
the CKTSO/SubtreeLU scheduling gap that prefix-trigger-only completion missed.
For final-triggered claims whose dependency scan already proves every U
predecessor is complete, KLS now tries the prepared-state BTF completion helper
with dependency waits disabled before falling back to the ordinary dispatcher.
It remains guarded and off by default: the top-five CKTSO-gap run
`build/kls_alg5_group_complete_waitfree_final_gap5_t4_r1_ref3_timeout120.jsonl`
raised claimed prepared currents on the ASIC rows into the hundreds and trimmed
the prior final-trigger run, but the SPICE-cycle geomean was still `2.3694s`
versus the same-build default
`build/kls_default_after_waitfree_gap5_t4_r1_ref3_timeout120.jsonl` at
`1.4558s`. Claim waits and per-column prepared-state completion still dominate,
so the remaining paper gap is a true grouped numeric executor for several live
current workspaces, not delayed scalar dispatch alone.
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_LIVE_STATE=1`
is the rollback-safe version of the rejected in-place retained-state suffix
experiment. It implies the direct-prefix final-state setup, requires the
completion path to own the current slot, logs original retained-state values
before mutating them, restores the log on fallback/error, and clears the slot
only after a successful completion. Correctness passes, but the focused
top-five run
`build/kls_alg5_group_complete_live_owned_gap5_t4_r1_ref3_timeout120.jsonl`
measured `2.5316s` geomean versus `2.4617s` for the non-live owned-state
guarded path. It removes retained-row copy counters on the ASIC rows, so the
remaining gap is not state-copy overhead; it is still the absence of a grouped
multi-current suffix update.
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_STATE_RAGGED_SUFFIX=1`
is a guarded diagnostic for that remaining suffix: it tries to apply matching
ragged-L supernode runs directly against the retained state map before falling
back to scalar dependencies. It is intentionally not a default. `ASIC_320ks`
slowed from `16.6473s` to `22.3606s` because the broad probe issued `187,757`
ragged-L attempts with `175,399` panel misses; this confirms the missing paper
piece is retained multi-current workspace ownership, not one-dependency-at-a-time
suffix rediscovery.
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_CLAIMS=1` is a narrower
trigger-timing diagnostic. It still retains the payoff descriptor for all
selected groups, but it only enables the mapped numeric payoff pattern and
claim scheduling when the total current-column surface is small enough for the
current serial follow-on executor. Large retained groups are deliberately left
as plan-only diagnostics because the CKTSO/SubtreeLU paper shape requires a
shared multi-current work queue, not one producer thread eagerly replaying all
future current columns. When active, the existing
`refactor_supernode_consumer_plan_claimed_columns`,
`refactor_supernode_consumer_plan_claim_skip_count`, and
`refactor_supernode_consumer_plan_claim_wait_count` counters show how many
columns were claimed and later skipped by the ordinary pipeline.
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_QUEUE=1` is the shared-work
queue diagnostic for that same retained Algorithm 5 payoff descriptor. It builds
the payoff groups without requiring the older ragged-L payoff numeric executor,
publishes current-column candidates when their producer-prefix trigger column
finishes, preclaims only candidates whose U predecessors are already complete,
and lets bounded queue consumers execute those columns through the normal
dependency-checked EGraph dispatcher. Pipeline workers that encounter a
preclaimed column drain payoff-queue work while waiting, so the diagnostic does
not rely on a separate worker staying free. Candidates that are not dependency
ready are left to the ordinary pipeline. The flag remains off-by-default:
focused Sandia probes showed that the queue is correct and bounded but does not
close the CKTSO gap, while broader claim-on-publish variants can stall large
cases. The existing claimed/skipped/wait counters report how much queued work was
actually taken.
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_QUEUE_PREFETCH=1` implies the
same retained payoff queue and enables an experimental CKTSO Algorithm 5
prefactor/postfactor scheduling probe: after a producer-prefix trigger column
finishes, KLS may queue up to `8 * threads` unclaimed current-column hints whose
full U dependency list is not complete yet. Queue consumers claim those hints
only if the dependencies are complete by the time the hint is popped; otherwise
the ordinary pipeline remains responsible for the column. This path requires the
guarded Algorithm 5 EGraph prefactor state, leaves the default queue behavior
unchanged, and remains off by default because it still does not add the
still-missing grouped multi-current numeric executor.
Set `KLS_ENABLE_REFACTOR_U_SUPERNODE_PATTERN=1` to also retain the exact
row-major U-supernode structural object from the same schedule pass. Benchmark
JSON reports `refactor_u_supernode_pattern_count`,
`refactor_u_supernode_pattern_rows`,
`refactor_u_supernode_pattern_max_width`,
`refactor_u_supernode_pattern_right_entries`, and
`refactor_u_supernode_pattern_internal_entries`. This is still diagnostic
metadata: it preserves the producer supernode start/width and right-side column
pattern for the future CKTSO/SubtreeLU-style numeric executor, while the current
default refactor continues to use the scalar KLU-format dependency walk.
Set `KLS_ENABLE_REFACTOR_U_SUPERNODE_VALUES=1` to additionally allocate the
matching row-major dense and right-side value buffers and record U entries as
the EGraph refactor publishes them. Benchmark JSON reports
`refactor_u_supernode_value_dense_entries`,
`refactor_u_supernode_value_right_entries`,
`refactor_last_u_supernode_value_dense_writes`,
`refactor_last_u_supernode_value_right_writes`, and cumulative
`refactor_u_supernode_value_*_write_count` counters. This value cache is also a
diagnostic staging object, not a default speed path: it proves the numeric data
can be materialized in the retained paper-style layout, but a future executor
still has to consume those buffers without paying an extra scalar recording pass.
Set `KLS_ENABLE_REFACTOR_U_SUPERNODE_RAGGED_L=1` to enable the next
off-by-default executor prototype for the same retained U-supernode ranges. It
keeps a ragged L-panel structure for U-supernode producers, publishes fresh L
values as producer columns finish in each numeric refactor, and consumes
eligible contiguous U-dependency runs through the cached internal panel plus
per-producer trailing L rows. Benchmark JSON reports
`refactor_u_supernode_l_panel_count`,
`refactor_u_supernode_l_dense_entries`,
`refactor_u_supernode_l_trailing_entries`,
`refactor_last_u_supernode_l_update_*`, and cumulative
`refactor_u_supernode_l_update_*` counters. After one clean numeric pass, the
prototype disables ragged L panels that did not feed any dependency run and
reports that amortization guard through
`refactor_u_supernode_l_prune_count`,
`refactor_u_supernode_l_pruned_panels`,
`refactor_u_supernode_l_pruned_dense_entries`, and
`refactor_u_supernode_l_pruned_trailing_entries`. This path is intentionally
opt-in: it validates a broader paper-style producer/consumer executor than the
common-trailing cached panel, but current focused runs still show that KLS needs
coarser batching/reuse before this shape can beat the scalar EGraph walk.
Set `KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_STATS=1` to measure the default BTF
EGraph scalar tail without changing execution. The refactor loop then records
contiguous producer runs whose U-dependency positions and local row ids advance
together, reporting `refactor_last_btf_scalar_run_candidates`,
`refactor_last_btf_scalar_run_rows`,
`refactor_last_btf_scalar_run_entries`,
`refactor_last_btf_scalar_run_max_rows`, and cumulative
`refactor_btf_scalar_run_*` counters. This is an opt-in diagnostic for the
missing paper-level numeric owner: it answers whether the scalar tail has enough
producer-run surface for a future row-major/current-state executor to consume.
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUPS=1` builds the next structural object
for that executor without changing numeric execution. It scans the retained BTF
U patterns for the same contiguous producer runs, groups them by producer start,
and retains only producer starts used by multiple current columns. Benchmark JSON
reports `refactor_btf_scalar_run_group_built`,
`refactor_btf_scalar_run_group_count`,
`refactor_btf_scalar_run_group_current_total`,
`refactor_btf_scalar_run_group_multi_count`,
`refactor_btf_scalar_run_group_multi_current_total`,
`refactor_btf_scalar_run_group_rows`,
`refactor_btf_scalar_run_group_reused_rows`,
`refactor_btf_scalar_run_group_entries`,
`refactor_btf_scalar_run_group_reused_entries`,
`refactor_btf_scalar_run_group_max_currents`, and
`refactor_btf_scalar_run_group_max_rows`. It also reports
`refactor_btf_scalar_run_group_producer_step_count`,
`refactor_btf_scalar_run_group_producer_step_multi_count`,
`refactor_btf_scalar_run_group_producer_step_active_currents`,
`refactor_btf_scalar_run_group_producer_step_unique_entries`,
`refactor_btf_scalar_run_group_producer_step_duplicate_entries`,
`refactor_btf_scalar_run_group_producer_step_reused_entries`, and
`refactor_btf_scalar_run_group_producer_step_max_currents`, which measure the
per-producer-row fanout a live grouped current-state owner would consume.
`refactor_btf_scalar_run_group_producer_index_count` and
`refactor_btf_scalar_run_group_producer_index_max_steps` report the retained
completed-producer lookup built from that fanout: each producer index entry owns
one or more group steps, and each step owns the active grouped-current members
for that producer row. The reused-entry counters estimate the L-entry stream
that such an owner could read once while updating several current workspaces.
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_LIVE_STATE_STATS=1` adds a
descriptor-level sparse-state union diagnostic for the same grouped current
sets. It reports `refactor_btf_scalar_run_group_live_state_groups`,
`refactor_btf_scalar_run_group_live_state_currents`,
`refactor_btf_scalar_run_group_live_state_rows`,
`refactor_btf_scalar_run_group_live_state_unique_rows`,
`refactor_btf_scalar_run_group_live_state_reused_rows`,
`refactor_btf_scalar_run_group_live_state_max_currents`,
`refactor_btf_scalar_run_group_live_state_max_rows`, and
`refactor_btf_scalar_run_group_live_state_max_unique_rows`. These counters
estimate how much per-current retained sparse state a grouped live workspace
could collapse before enabling any retained-state numeric executor.
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_LIVE_STEP_STATS=1` measures the
same row-union surface at producer-step granularity, after the completed
producer index has selected the active current memberships for each producer
row. Benchmark JSON reports
`refactor_btf_scalar_run_group_live_step_count`,
`refactor_btf_scalar_run_group_live_step_current_total`,
`refactor_btf_scalar_run_group_live_step_rows`,
`refactor_btf_scalar_run_group_live_step_unique_rows`,
`refactor_btf_scalar_run_group_live_step_reused_rows`,
`refactor_btf_scalar_run_group_live_step_max_currents`,
`refactor_btf_scalar_run_group_live_step_max_rows`, and
`refactor_btf_scalar_run_group_live_step_max_unique_rows`.
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_LIVE_STEP_PLAN=1` additionally keeps
a per-step unique-row prefix plus the sorted row descriptor for each retained
step when memory allows. `refactor_btf_scalar_run_group_live_step_stored_rows`
reports the retained row-index entries, and
`refactor_btf_scalar_run_group_live_step_storage_limited` reports when KLS
counted the surface but could not retain the descriptor. When producer-step
retained-state advance is also enabled, KLS records runtime full-step coverage
through `refactor_last_btf_scalar_run_group_live_step_runtime_full_steps`,
`refactor_last_btf_scalar_run_group_live_step_runtime_full_currents`, and
`refactor_last_btf_scalar_run_group_live_step_runtime_full_rows`, plus matching
cumulative counters. It also records partial-step coverage through
`refactor_last_btf_scalar_run_group_live_step_runtime_partial_steps`,
`refactor_last_btf_scalar_run_group_live_step_runtime_partial_currents`,
`refactor_last_btf_scalar_run_group_live_step_runtime_partial_members`, and
`refactor_last_btf_scalar_run_group_live_step_runtime_partial_rows`, plus
matching cumulative counters. Full counters identify claimed producer batches
whose active members exactly match a retained live-step descriptor. Partial
counters identify descriptor surface reached by runtime producer batches whose
memberships are split across scheduler ownership. Numeric execution is still
unchanged.
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_LIVE_STEP_PARTIAL_OWNER_STATS=1`
also builds the live-step plan and, during retained-state step advance, measures
the duplicate and unique state rows for only the currently claimed members of
partial live steps through
`refactor_last_btf_scalar_run_group_live_step_runtime_partial_active_rows` and
`refactor_last_btf_scalar_run_group_live_step_runtime_partial_active_unique_rows`.
Those counters size the active subset that a future partial live-step owner
would gather, which can be much smaller than the full structural descriptor.
The retained group descriptor is tied to the LU pointer cache and is reused
across repeated numeric refactors until the numeric pattern changes.
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_WAIT_STATS=1` also builds that
descriptor and records live runtime overlap at the grouped producer-run wait
point. It does not change numeric execution. Benchmark JSON reports last and
cumulative grouped wait counts, overlapped wait counts, wait rows/entries, and
the maximum number of live waiters on a grouped producer run. These counters
answer whether the structural paper-style reuse surface is actually reached by
multiple current columns at the same time under the current BTF schedule.
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_CLAIMS=1` makes those grouped
producer runs schedule-visible without changing the numeric kernel. When a
producer column finishes, KLS claims grouped current columns whose full
dependency lists are already published, dispatches them on the same worker, and
marks them done through the normal pipeline path. Benchmark JSON reports last
and cumulative trigger/group/current surface counts plus the subset actually
claimed. This opt-in scheduler probe tests whether producer-side current
locality helps before adding a true multi-workspace grouped numeric owner.
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_PREFIX_STATS=1` keeps the same
producer-completion hook diagnostic-only and asks an earlier paper-level
question: for each grouped current membership, are all dependencies before the
contiguous producer run already published, and is the full producer run itself
ready? Benchmark JSON reports last and cumulative prefix-ready and run-ready
current/row/entry counts. A large prefix-ready count with a tiny run-ready count
means the next executor must retain live partial current state and wake it again
as later producers in the contiguous run complete; a full-column claim at the
producer-start hook is too late and too restrictive.
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_WAKE_STATS=1` adds the matching
run-end schedule probe. When the producer-start hook finds a prefix-ready
membership, KLS arms that membership in a per-refactor live array. The retained
descriptor also indexes memberships by the last producer in the contiguous run,
so the completion hook can count armed memberships whose full producer run has
become ready. Benchmark JSON reports last and cumulative wake-armed and
wake-ready current/row/entry counts. This remains diagnostic, but it is the
scheduler skeleton needed by a retained current-state owner that will store the
prefix state at arm time and apply shared producer-run updates at wake time.
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_STATE_STATS=1` extends that wake
probe with the first guarded retained-state materialization. The flag implies
the wake probe, allocates compact row-state storage for each grouped membership,
copies the prefix-ready current state when the membership is armed, and applies
the contiguous producer run into that compact state when the wake fires.
Benchmark JSON reports the planned state rows, the materialized current/row and
prefix-dependency counts, the advanced current/row/entry counts, and state
rejects. This still does not replace the normal numeric refactor; it validates
the paper-level live current-state shape before KLS commits it to the production
numeric path.
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_STATE_STEP_ADVANCE=1` makes that
retained state consume the retained producer-step schedule incrementally. The
flag implies state stats, tracks each materialized membership's next producer
row, and advances it from the producer-completion hook through the existing
per-member wake-state CAS. Benchmark JSON reports last and cumulative
step-advance trigger, step, current, row, entry, ready-current, and reject
counts, plus producer-batch step/current/entry counts when one completed
producer column is streamed across multiple claimed retained states. This is
still an opt-in retained-state experiment: it tests whether completed-producer
scheduling closes the paper gap before a broader grouped workspace owner
streams one producer L column across several current states.
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_STATE_STEP_WINDOW=1` extends that
producer-step path into a guarded grouped retained-state window. The flag
implies step advance, keeps successfully claimed retained states owned after
the first producer update, then advances additional ready producer rows when at
least two active retained states share the same completed producer. Benchmark
JSON reports last and cumulative grouped-window round, current, and entry
counts. The path remains opt-in while it is compared against CKTSO-gap cases.
`KLS_REFACTOR_BTF_SCALAR_RUN_GROUP_STATE_STEP_POSITION_MAX_ENTRIES=N` enables
an additional retained-state row-position matrix inside that opt-in step
window when the active producer update fits `N` stored positions per worker.
The default is `0`, so KLS keeps the direct grouped row stream unless this
experimental matrix is explicitly requested.
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_STATE_EXEC=1` turns that retained
state into a guarded numeric consumer. The flag implies state stats, stores the
U coefficients consumed while preparing the retained state, arms only the
farthest retained member per current column, restores that state when the
column is dispatched through either wait or no-wait plain-scalar BTF paths,
writes the skipped U entries, and lets the existing suffix, pivot check, and
L-store path finish the column. When the retained state is already terminal,
the executor publishes `U`, the diagonal, and `L` directly from the retained
state map, using the sorted retained rows as a merge stream when available, so
that terminal columns do not copy the sparse state back through the scalar
workspace. The same terminal publisher can also run at the producer wake point:
if all dependencies are already done, KLS claims the future current through the
normal claimed-column array and later level-synchronous workers skip it after
observing the done generation. Benchmark JSON reports last and cumulative
state-exec current, skipped-dependency, restored-row, terminal-current,
wake-terminal current, wake-terminal candidate, wake-terminal nonpipeline,
wake-terminal dependency-miss, wake-terminal claim-miss,
remaining-suffix-dependency, remaining-suffix-entry, and reject counts. It also
reports dispatch bypasses where the BTF column had a selected retained state but
the restore hook was gated off, and splits remaining misses into
not-ready-at-dispatch, ready-while-owned, and ready-after-done buckets. When
this executor is requested, KLS also records whether retained L row lists are
sorted and advances retained sparse states with a merge walk instead of a
per-entry binary search when that sorted shape holds.
This remains opt-in: the focused ASIC probes are residual-clean, but the state
owner is still too fine-grained to beat the default refactor path.
The executor is now additionally guarded by structural state/payoff counters:
`refactor_btf_scalar_run_group_wake_count`,
`refactor_btf_scalar_run_group_wake_member_total`,
`refactor_btf_scalar_run_group_state_current_count`,
`refactor_btf_scalar_run_group_state_best_skip_total`, and
`refactor_btf_scalar_run_group_state_max_best_skip`. Unless
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_STATE_EXEC_UNGUARDED=1` is set, the
state executor and step-window runtime are disabled when the retained best-skip
surface is smaller than the retained state-row surface. This keeps rejected
CKTSO/SubtreeLU state-owner probes from running the expensive wake/materialize
path when the paper-level state geometry is already upside down. The exact
state plan is retained across repeated refactors while the LU pointer cache is
valid, so repeated SPICE refactors do not rebuild the same state-row descriptor.
The guarded path also runs a cheaper lower-bound check before allocating exact
per-member state rows: if the best skipped U-dependency surface is already
smaller than the sum of member skip lengths, the exact retained-state plan
cannot pass the same payoff guard. Benchmark JSON reports this as
`refactor_btf_scalar_run_group_state_guard_lower_bound_rows` and
`refactor_btf_scalar_run_group_state_guard_lower_bound_rejected`.
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_STATE_COMPACT_VALUES=1` is an
additional opt-in for the retained-state executor. When combined with
`STATE_EXEC`, it stores retained state values only for memberships belonging to
columns that have a selected retained executor state, while keeping the full
structural row descriptors. This is a storage substrate for a future grouped
live-current owner, not a default execution mode.
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_STATE_ADVANCE_BATCH_STATS=1` is a
separate diagnostic for the next grouped-current owner. It implies state stats
and reports last/cumulative advance-batch group, current, unique-entry,
duplicate-entry, and max-current counts. When producer-step retained-state
advance is also enabled, it additionally reports the duplicate sparse state rows
covered by producer-step batches and the exact unique global state rows those
batches cover. The normal state executor leaves these counters at zero unless
this diagnostic is enabled, so timing runs do not pay for the extra wake-loop
measurement.
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_EXEC=1` enables a guarded executor for the
same BTF scalar producer-run shape. It first waits for every dependency in a
contiguous run, then applies the run through a local row workspace. Under the
same flag, KLS also retains an owned sorted L-row mirror: the in-run triangular
prefix is applied locally and the trailing suffix is sent through the existing
scatter primitive without rewriting the KLU-owned numeric factor.
Benchmark JSON reports the applied surface through
`refactor_last_btf_scalar_run_exec_runs`,
`refactor_last_btf_scalar_run_exec_rows`,
`refactor_last_btf_scalar_run_exec_entries`,
`refactor_last_btf_scalar_run_exec_max_rows`, and cumulative
`refactor_btf_scalar_run_exec_*` counters. This is intentionally off by default:
the direct unsorted-row prototype was residual-clean but slower on `ASIC_100ks`,
and the sorted-mirror executor activates the intended large producer-run surface
but still remains slower than the default scalar EGraph path on the focused ASIC
cases. The useful conclusion is that closing this refactor gap needs grouped
current-state arithmetic over this row-ordered representation, not another
wrapper around the scalar KLU scatter layout.
Eligible retained refactor-map row/input positions and L row-index arrays are
mirrored as 32-bit integers by default while leaving the KLU-owned numeric
factor and public index ABI unchanged. The EGraph value-scatter path and
refactor pool scatter kernels use the narrower mirrors when the matrix fits
32-bit local row/input indices and fall back to the original `UF_long` arrays
otherwise. The retained U dependency-index arrays use the same default 32-bit
mirror for eligible matrices, letting the EGraph scalar dependency walk read
narrower producer indices without changing the KLU-owned numeric factor. Set
`KLS_ENABLE_REFACTOR_MAP_INDEX32=0`, `KLS_ENABLE_REFACTOR_L_INDEX32=0`, or
`KLS_ENABLE_REFACTOR_U_INDEX32=0` to disable individual mirrors for A/B
comparisons. Benchmark JSON reports `refactor_map_index32_enabled` and
`refactor_map_index32_entries`, `refactor_l_index32_enabled` and
`refactor_l_index32_entries`, plus `refactor_u_index32_enabled` and
`refactor_u_index32_entries` so runs can verify whether each mirror was active.
The EGraph numeric kernels include a guarded CKTSO Algorithm 5-style prefactor
slice. The production generic policy leaves this scalar experiment disabled:
enabling it globally excludes the narrower fused BTF dispatch even when no
profitable speculative dependency is found. Set
`KLS_ENABLE_EGRAPH_ALGORITHM5_PREF_UPDATE=1` to force it or
`KLS_ENABLE_EGRAPH_ALGORITHM5_PREF_UPDATE=0` to disable it for A/B runs; legacy
compatibility mode retains its historical work-threshold selection. When a
pipeline column is
blocked on its current U predecessor, the kernel scans later U predecessors that
are already published, applies only those whose workspace entry cannot be
changed by any earlier unapplied predecessor, records them in an applied bitmap,
and skips them when the normal postfactor cursor reaches that position. This is
deliberately scalar. It only uses already-existing pipeline completion state, so
using this path does not create `pipeline_done` state or change the scheduler
shape by itself. The applied bitmap is allocated only after
a dependency actually blocks and is initialized with the already-consumed
prefix. The path now covers the single-block unscaled/scaled kernels, the
unscaled BTF kernel, and the generic scaled/fallback kernel. This fills more of
the paper's skip-unfinished prefactor semantics without pretending that the
grouped Algorithm 5 payoff descriptor already has a multi-current numeric
kernel. Benchmark JSON reports
`refactor_last_egraph_algorithm5_prefactor_columns`,
`refactor_last_egraph_algorithm5_prefactor_deps`,
`refactor_egraph_algorithm5_prefactor_column_count`, and
`refactor_egraph_algorithm5_prefactor_dep_count`.
The separate repeated-factor repair guard still blocks high-work pipeline
factors when pivot repair would be risky, but no longer blocks numerics whose
retained factor has zero off-diagonal pivots. Those no-pivot cases reuse the
checked EGraph fast refactor and still fall back to full KLU factorization if
the pivot check fails.
The experimental row pipeline preserves the CKTSO-style wide cluster prefix
selected by the `2 * threads` width rule, then consumes the remaining narrow
tail through a bounded successor-ready queue when explicit predecessor counts
are available.
If that queue cannot be prepared, it falls back to the older barriered cluster
levels plus queued tail. Stats report
`row_refactor_last_ready_queue`, `row_refactor_ready_queue_run_count`, and
`row_refactor_ready_queue_group_count` when the queued group scheduler is used.
Unchecked queued row refactors can run without the per-row completion bitmap;
stats report `row_refactor_last_done_bitmap` and
`row_refactor_done_bitmap_run_count` so checked pivot-prefix validation remains
visible. The partial row-dependency queue releases successor groups only when
their exact row-predecessor count reaches zero. An earlier speculative
Algorithm 5-style row-prefactor release could enqueue consumers before that
point and forced large dense groups into expensive wait/prefactor loops; it is
no longer used by the default row-refactor partial queue. The historical
prefactor diagnostics remain reported through `row_refactor_last_prefactor`,
`row_refactor_last_prefactor_rows`, `row_refactor_last_prefactor_deps`,
cumulative `row_refactor_prefactor_*`, and the corresponding
`row_refactor_prefactor_supernode_*` counters; exact row-dependency runs should
leave those counters at zero.

KLS-first row-up-looking factorization reports cumulative
`kls_first_row_supernode_update_groups` and
`kls_first_row_supernode_update_rows`, plus the panel-backed subset as
`kls_first_row_supernode_panel_update_groups` and
`kls_first_row_supernode_panel_update_rows`. These cumulative counters remain
useful when later fallback or refactor bookkeeping clears the volatile
`kls_first_last_row_supernode_*` fields.

Row-pattern analysis also records
`row_refactor_input_cleanup_rows`
and `row_refactor_input_cleanup_entries`; rows whose input columns are already
covered by `L`, the pivot, or `U` skip the redundant residual cleanup loop in
the row numeric kernels. Experimental row refactors handle single-block factors
and BTF diagonal blocks; BTF off-block values are refreshed into KLU `Offx`
from the retained input map. KLU row-scaled row refactors recompute `Rs`, scale
fixed-position input values by the unpermuted row scale, and permute `Rs` back
to pivot order after an accepted row pass. Row refactors can
keep KLS row-major `L`/`U` values authoritative across repeated unchecked
refactors and publish them to KLU storage lazily only when a KLU fallback needs
it. The guarded dirty row solve
handles unscaled and KLU row-scaled normal factors in forward and transpose
orientation by matching KLU's `P*(R\b)` and `Q'*b`/`P'*(R\x)` solve setup; for
BTF factors it follows KLU's block solve order and applies refreshed `Offx`
coupling without first publishing row mirrors to KLU `L`/`U` storage. The dirty
row solve traverses each sparse block once for up to four right-hand sides at a
time, matching KLU's small-RHS batching shape while reading KLS-owned row
mirrors. External KLS row/column scaling or permutation, later factorization,
and other non-row fallbacks still publish the dirty row mirrors before using KLU
storage.
Benchmark stats
report `row_refactor_last_defer_value_scatter`,
`row_refactor_defer_value_scatter_run_count`, `row_refactor_values_dirty`,
`row_refactor_last_lazy_value_scatter`, and
`row_refactor_lazy_value_scatter_run_count`, plus
`row_refactor_last_row_solve` and `row_refactor_row_solve_run_count` for this
row-storage solve handoff. When `KLS_ENABLE_ROW_SOLVE_FROM_NUMERIC=1`, KLS can
build row-solve structure from ordinary KLU-compatible numeric storage after a
successful factor/refactor, allowing the row solve to be benchmarked without
also enabling the experimental KLS-first factorization path. This automatic
seed is adaptive: KLS first builds the cheap row-solve partition diagnostics and,
when the CKTSO-style parallel row-solve executor has enough structural work,
reads values through retained KLU numeric value pointers instead of copying them
into row-major mirrors. Otherwise the ordinary numeric storage remains
authoritative and KLS avoids a repeated setup tax on refactor-heavy runs.
Single-RHS row solves use a scalar row-major loop with cached structural
validation, matching the common SPICE solve shape while keeping the four-RHS
batched path available for wider solves. Solve-only seeding builds only the
row-major `L`/`U` solve structure and KLU value-pointer arrays; the heavier
row-refactor group, segment, and scheduler metadata is left to the refactor
paths that actually need it.

For sparse, normal-orientation AMD factors with a dominant BTF core and a
10--20% block fringe, KLS selects the column-solve representation from retained
numeric state rather than matrix dimensions. Moderate-work unscaled factors
without static matching keep KLU's native packed indices; a retained static
matching permutation selects the 32-bit mirror and precomposes public RHS rows
with numeric rows. The class is bounded by normalized input density, core
coverage, factor fill/work, rank, and pivot-stability checks. Its analyze-time
counterpart also suppresses a speculative METIS/scale race only when symbolic
fill and work fall in the same moderate band. Set
`KLS_DISABLE_SPARSE_FRAGMENTED_RACE_SUPPRESSION=1`,
`KLS_DISABLE_MODERATE_FRAGMENTED_I32_SOLVE=1`, or
`KLS_DISABLE_MODERATE_FRAGMENTED_FUSED_RHS=1` for independent A/B diagnosis.

A separate retained-numeric profile supplies lifecycle defaults for stable,
moderate-work members of that fragmented dominant-BTF class. It requires
3--8 input entries per row, a BTF block count between 1/12 and 1/5 of the
order, an 80--95% dominant core, at most `n/128` off-diagonal pivots, no nudges
or perturbations, 12--40 retained factor entries per row, and 1,000--8,192
factor flops per row within a 100M--1B total-work band. The selected factor
must be normal AMD, full-rank, unscaled, and have no retained numeric `Rs`
array; explicit normal/AMD/unscaled requests are eligible after selecting the
same numeric state. Consumers that already admit EGraph can then reuse the
full worker crew, scalar scatter,
wide subtree/PTS setup, compact solve offsets where representable, and
overlapped preparations without an input-size fingerprint. The pre-static
factor-worker counterpart uses the trial symbolic's normalized SCC shape,
fill, and estimated work because no accepted numeric exists yet. Benchmark
JSON exposes `moderate_fragmented_policy_eligible`. Set
`KLS_DISABLE_MODERATE_FRAGMENTED_EGRAPH_POLICY=1` to disable the bundle.
The former `KLS_*ONETONE2*` environment spellings are no longer recognized.

A compact AMF two-block policy is selected from computed symbolic and numeric
state rather than an input fingerprint. The symbolic proposal requires
512--4096 rows, 8--16 input entries per row, full structural rank, exactly two
BTF blocks with a core covering at least 99% of the order, 40--64 estimated
factor entries per row, and 1,000--2,048 estimated flops per row. Each
estimated triangular stream must fit in 16 bits; the 4096-row ceiling is the
packed dependency descriptor's 12-bit row field. The recurring lifecycle is
enabled only after the actual unscaled normal-AMF numeric has no nudges or
perturbations, at most `n/64` off-diagonal pivots, 20--32 factor entries per
row, and 200--512 factor flops per row, with both triangular streams still
16-bit representable. Packed workers additionally validate every derived
pointer, index, schedule, and solve representation before dispatch. AUTO uses
the ordinary orientation/ordering selector, while explicit
normal/AMF/scale-0 requests are eligible after reaching the same state. Every
changed numeric retains the strict relative-L2 residual guard. Benchmark JSON
exposes `compact_amf_two_block_policy_eligible`. Set
`KLS_DISABLE_COMPACT_AMF_TWO_BLOCK_POLICY=1` to disable the whole policy or
`KLS_DISABLE_COMPACT_AMF_TWO_BLOCK_SPECIALIZED_WORKER=1` to retain the policy
without its specialized worker. The former `KLS_*TSOPF_B9*` controls are no
longer recognized.

A dense reciprocal-hub METIS policy replaces the former dense-ASIC size
window. Normal AUTO input with at least 65,536 rows only proposes the route
when it has 8--12 entries per row, at least a 99.5% structural diagonal, and a
vertex covering at least seven eighths of the graph in both row and column
directions. The actual NodeNDP symbolic must then be full-rank, have both its
BTF fringe and block count between `n/512` and `n/64`, and expose a complete
separator forest with useful private and pipeline components. Rejection
resumes ordinary AUTO orientation and ordering selection. Recurring kernel
defaults additionally require the measured unscaled factor to have no
off-diagonal pivots or numeric `Rs`, 24--48 factor entries per row, and
4,096--8,192 factor flops per row. Benchmark JSON exposes
`dense_reciprocal_hub_policy_eligible`; set
`KLS_DISABLE_DENSE_RECIPROCAL_HUB_METIS_POLICY=1` for an A/B fallback.

A symmetric scalar-fringe AMD policy replaces the former `rajat03` size and
nonzero window. Under the standard AUTO/8-thread BTF/static-pivoting contract,
it proposes transpose AMD only for bounded-degree, exactly structurally
symmetric inputs with 3--6 entries per row and a material but small surface of
missing-diagonal degree-one columns. The selected symbolic must be full-rank,
place a 1/256--1/28 fringe outside a dominant core, expose nearly all of that
surface as scalar BTF blocks, and stay within normalized fill/work and L/U
balance bounds. Recurring lean kernels additionally require an unscaled
numeric with no off-diagonal pivots, nudges, perturbations, or numeric `Rs`
array and the same measured factor bounds. Benchmark JSON exposes
`symmetric_scalar_fringe_policy_eligible`; set
`KLS_DISABLE_SYMMETRIC_SCALAR_FRINGE_AMD_LEAN_POLICY=1` for an A/B fallback.

A pivoted high-work single-block policy replaces the former `rajat15` order
and nonzero window. Under the standard AUTO/8-thread BTF/static-pivoting
contract, only a retained normal-AMD single block can qualify. The symbolic
state must have 6--16 input entries per row, balanced L/U estimates, 32--96
estimated factor entries per row, and 2,048--8,192 estimated flops per row
inside broad order and total-work resource bounds. Recurring kernels are
enabled only after the no-scale/no-recheck numeric has no `Rs`, nudges, or
perturbations; between `n/64` and `n/8` off-diagonal pivots; balanced L/U
storage; and normalized measured fill and work. Symbolic admission does not
suppress AUTO's ordinary ordering race because identical patterns can produce
different pivot behavior. Benchmark JSON exposes
`pivoted_high_work_single_block_policy_eligible`; set
`KLS_DISABLE_PIVOTED_HIGH_WORK_SINGLE_BLOCK_POLICY=1` for an A/B fallback.

A low-work many-fringe BTF/PTS policy replaces the former `rajat21` order,
nonzero, block-count, and core-size box. Under the standard AUTO/8-thread
BTF/static-pivoting contract, a retained normal-AMD factor may qualify from
131,072--1,048,576 rows and 3--8 input entries per row. It must be full-rank,
put between `n/64` and `n/16` vertices outside one dominant BTF core, expose
between half and all of that fringe as separate blocks, keep its symbolic L/U
estimates within 2:1, and meet normalized symbolic fill and work bounds.
Recurring PTS choices additionally require an unscaled numeric with no `Rs`,
nudges, or perturbations; at most `n/512` off-diagonal pivots; balanced L/U
storage; 4--16 retained factor entries per row; and 8--128 measured flops per
row inside a 2M--64M resource band. Benchmark JSON exposes
`low_work_many_fringe_btf_pts_policy_eligible`; set
`KLS_DISABLE_LOW_WORK_MANY_FRINGE_BTF_PTS_POLICY=1` for an A/B fallback.

A hubbed scalar-fringe subtype replaces the former `rajat29` order and
nonzero window.  A cheap input proposal requires the same broad order and
density range as the low-work many-fringe family, at least 31/32 of columns
to contain their diagonal, between `n/64` and `n/16` scalar columns, and one
column whose degree is between `n/8` and `7n/8`.  AUTO retains the proposed
normal-AMD symbolic only when it is full-rank, leaves an `n/64`--`n/16`
fringe around one dominant BTF core, represents at least seven eighths of
that fringe as separate blocks, and passes normalized fill, work, and L/U
balance bounds.  The measured numeric then inherits the low-work PTS guards;
this more specific topology admits a 192-flop-per-row symbolic floor and at
most `n/256` off-diagonal pivots.  Explicit orientation,
ordering, scale, backend, BTF, pivot, or thread choices remain authoritative.
Benchmark JSON exposes
`low_work_hubbed_scalar_fringe_pts_policy_eligible`.  Set
`KLS_DISABLE_LOW_WORK_HUBBED_SCALAR_FRINGE_PTS_POLICY=1` to disable only this
subtype, or use the many-fringe master switch above to disable the whole
family.

A high-work tiny scalar-fringe AMD policy replaces the former `Raj1` order,
nonzero, block-count, and core-size box. Under the standard AUTO/8-thread
BTF/static-pivoting contract, a retained normal-AMD symbolic from
131,072--1,048,576 rows and 3--12 input entries per row may qualify. It must
be full-rank, leave between `n/4096` and `n/256` vertices outside one dominant
core, represent at least three quarters of that fringe as separate BTF blocks,
keep estimated L/U within 2:1, and meet normalized symbolic fill and work
bounds. The compact-fill half may use probe-first predicted construction;
denser members retain the ordinary nudge sequence. Recurring EGraph choices
additionally require an unscaled numeric with no `Rs`, nudges, or
perturbations; either a verified predicted no-pivot factor or at most `n/512`
off-diagonal pivots; balanced L/U storage; and normalized measured fill and
work. Benchmark JSON exposes
`high_work_tiny_scalar_fringe_policy_eligible`; set
`KLS_DISABLE_HIGH_WORK_TINY_SCALAR_FRINGE_POLICY=1` for an A/B fallback.

The metamorphic holdout utility can produce deterministic simultaneous
row/column relabelings without checking generated matrices into the tree. For
example, the compact positive used by the policy audit is reproduced with:

```sh
python3 scripts/make_metamorphic_matrix.py \
  data/suitesparse/TSOPF/TSOPF_FS_b9_c1.mtx \
  /tmp/TSOPF_FS_b9_c1_local_holdout.mtx \
  --adjacent-swaps 100 --seed 211
```

Use `--shuffle --seed N` instead of `--adjacent-swaps` for a complete
deterministic random relabeling. `--append-diagonal-blocks COUNT` adds
independent unit diagonal components, while `--append-coupled-nodes COUNT`
adds weak reciprocal diagonal nodes, as separate transformation modes. These
transformations are useful for detecting classifiers or separator policies
that accidentally depend on the original vertex numbers, an exact matrix
order, or an exact component count.

For the corresponding AUTO analysis problem, a 40--60% structural diagonal
and three-to-five input entries per row only propose a direct AMD/BTF
candidate. KLS retains that candidate when the computed symbolic has full
structural rank, a 10--20% block fringe, an 80--95% dominant core, estimated
fill of 32--64 entries per row, and estimated work of 1,000--8,192 flops per
row. The accepted state stays unscaled and does not pay later row-matching,
scale, or METIS replacement trials. A rejected proposal resumes ordinary AUTO
ordering and orientation selection. Set
`KLS_DISABLE_SPARSE_PARTIAL_DIAGONAL_DIRECT_AMD=1` to disable the policy.

When row refactor has retained exact compact dense groups, normal and transpose
row solves consume complete groups as row-major triangular panels for one RHS
and four-RHS chunks, falling back to scalar row loops when a group does not
match the compact lower/upper layout exactly. The transpose path mirrors the
normal compact executor by solving `U^T` from the retained upper panel and
`L^T` from the retained lower panel instead of re-reading every in-panel entry
through scalar accessors.

It still records CKTSO-style triangular solve partition diagnostics using the
paper's dense-tail criteria, namely a suffix with at least 70% of
row-triangular entries and at least 300,000 entries, plus the eight trapezoid
slices CKTSO uses after a dense tail is found. KLS retains the slice boundaries
internally and reports the maximum per-slice triangular entries as a
load-balance diagnostic. It also
precomputes per-row rectangular/triangular split points for dense-tail rows and
reports the resulting lower/upper segment entry counts, matching CKTSO's
trapezoid-slice setup. Single-RHS, single-block normal and transpose solves now
use the persistent worker pool to parallelize the rectangular part of those
slices when the parallel rectangular work reaches the paper's 300,000-entry
dense-tail work scale. Transpose solves use retained transposed row views of
`U^T` and `L^T` with source positions back to the KLU-compatible factor values,
so they gather solved dependencies instead of racing through row-scatter
updates. Rectangular slice rows are assigned to workers by accumulated
rectangular-entry counts, and KLS reports the max per-thread rectangular entries
for lower/upper factors. It also levelizes the sparse triangular block before
the dense tail and solves wide prefix levels in the same persistent worker pool,
falling back to a sequential loop for the remaining narrow levels; the
within-slice triangular pieces remain sequential. The parallel row-solve
executor is conservatively gated by parallelizable work share and work per
synchronization, so cases where barrier overhead dominates keep the scalar
row-major solve.
The queued
row scheduler orders
ready groups by the retained FLOP-style group work estimate, including
successors released by completed groups, and reports
`row_refactor_last_work_ready_queue` plus
`row_refactor_work_ready_queue_run_count`. The ready queue keeps solver-owned
workspace across repeated row refactors and reports its capacity through
`row_refactor_ready_queue_workspace_groups`. When the queued scheduler keeps a
newly ready successor as a worker-local continuation instead of spilling it to
the shared queue, stats report `row_refactor_last_local_ready_groups` and
`row_refactor_local_ready_group_count`. Tail queued runs reuse the retained
group predecessor counts; when the cluster prefix is empty, full-graph queued
runs also reuse the retained root-group list and hand the root groups out
through a private-root cursor before using the shared queue for newly released
successors. Cluster-prefix tail queues now apply the same private first-wave
treatment to the initially ready tail groups when per-thread queues can be
built, reported as `row_refactor_last_private_ready_groups` and
`row_refactor_private_ready_group_count`. When retained METIS `NodeNDP`
separator trees cover the factor order, including BTF analyses where local
separator trees are stitched into a global forest, row refactors first try a
SubtreeLU Algorithm 6-style FLOP-balanced separator queue: the dominant
separator subtree is repeatedly split into a pipeline root plus child subtrees,
remaining active components are assigned to private thread queues by retained
component work, and separator-crossing row groups are forced into the pipeline
queue.
Indivisible retained components are left in the private candidate set, so an
unbalanced separator tree cannot collapse the Algorithm 6 queue into all
pipeline work with no private subtrees.
If a private subtree group depends on an already-pipeline group, KLS promotes
that dependent group into the pipeline closure instead of discarding the
separator schedule.
Stats report this path through
`row_refactor_last_separator_flop_queue`,
`row_refactor_separator_flop_queue_run_count`,
`row_refactor_last_separator_flop_ordered_private`,
`row_refactor_separator_flop_ordered_private_run_count`,
`row_refactor_last_separator_flop_components`,
`row_refactor_last_separator_flop_private_groups`, and
`row_refactor_last_separator_flop_pipeline_groups`, plus closure promotions via
`row_refactor_last_separator_flop_closure_groups` and
`row_refactor_separator_flop_closure_group_count`. When a completed group releases
multiple successors, the completing
worker keeps one local continuation and only spills the rest to the shared
queue. Checked queued
rejects refresh any missing prefix rows before accepting a prefix-tail repair
classification.
Focused CKTSO-gap probes show that this retained scheduler is not, by itself,
the missing row-refactor core: forced row refactors on the ASIC cases can place
almost every tiny group in private separator queues, while opt-in compact
supernode/native-panel runs over-stage more panel entries than the scalar row
work they replace. The default automatic row-refactor handoff therefore keeps
using the lower-bound gate, and the benchmark decomposer now labels these cases
as `row_refactor_private_scalar_scaffold`,
`row_refactor_panel_overstaged`, or `row_refactor_lower_bound_rejected`.
Those signals distinguish the remaining paper gap, a production row-major
supernode/current-state numeric executor, from the already-present scheduling
and panel probes.
Checked row fast-factor rejects also report the conservative row-group restart
tail through
`fast_rejected_group_tail_groups`, `fast_rejected_group_tail_rows`, and
`fast_rejected_group_tail_work`, giving the row/segment task graph a visible
scope for future CKTSO-style pivoting tail restart.
Because checked row fast-factor probes and
unchecked row refactors are selected independently, stats also report
`row_refactor_last_run`, `row_refactor_last_checked`,
`row_refactor_last_parallel`, and row-refactor run counters to identify the
last numeric kernel actually used.
Checked EGraph refactors also keep a completed-column bitmap for barriered
cluster-only runs, so a pivot reject is reported as prefix-current when every
earlier factor-order column is proven finished.
Benchmark stats report `fast_rejected_pivot`,
`fast_rejected_pivot_col`, `fast_rejected_row`,
`fast_rejected_multiplier_abs`, `fast_rejected_pivot_abs`,
`fast_rejected_candidate_abs`, `fast_rejected_tail_candidate_row`,
`fast_rejected_tail_candidate_abs`, `fast_rejected_tail_candidate_count`,
`fast_rejected_tail_candidate_position`,
`fast_rejected_tail_repair_ready`,
`fast_repaired_pivot_row`,
`fast_repaired_pivot_matches_tail_candidate`,
`fast_repaired_first_changed_pivot`,
`fast_repaired_prefix_changed_pivots`,
`fast_repaired_suffix_changed_pivots`,
`fast_repaired_tail_restart_ready`,
`fast_repaired_block_work`,
`fast_repaired_tail_restart_columns`,
`fast_repaired_tail_restart_work`,
`fast_repaired_tail_restart_saved_work`,
`fast_repaired_tail_restart_overcompute_columns`,
`fast_repaired_tail_restart_overcompute_work`,
`fast_repaired_tail_restart_exact_mask`,
`fast_repaired_tail_restart_etree_mask`,
`fast_rejected_block_start`, `fast_rejected_block_size`,
`fast_rejected_suffix_columns`,
`fast_rejected_descendant_columns`, `fast_rejected_descendant_work`,
`fast_rejected_row_tail_columns`, `fast_rejected_row_tail_work`,
`fast_rejected_etree_columns`, `fast_rejected_etree_work`,
`fast_rejected_pivoting_tail_columns`,
`fast_rejected_pivoting_tail_work`,
`fast_rejected_pivoting_tail_first`,
`fast_rejected_pivoting_tail_last`,
`fast_rejected_pivoting_tail_contains_reject`,
`fast_rejected_pivoting_tail_topological`,
`fast_rejected_pivoting_tail_seed_columns`,
`fast_rejected_pivoting_tail_row_seed_columns`,
`fast_rejected_pivoting_tail_block_seed_columns`,
`fast_rejected_pivoting_tail_contiguous`,
`fast_rejected_pivoting_tail_suffix_exact`,
`fast_rejected_pivoting_tail_gap_columns`,
`fast_rejected_pivoting_tail_suffix_overcompute_columns`,
`fast_rejected_pivoting_tail_suffix_overcompute_work`,
`fast_rejected_pivoting_tail_etree_edges`,
`fast_rejected_pivoting_tail_etree_roots`,
`fast_rejected_pivoting_tail_etree_leaves`,
`fast_rejected_pivoting_tail_etree_max_fanout`,
`fast_rejected_pivoting_tail_etree_levels`,
`fast_rejected_pivoting_tail_etree_max_width`,
`fast_rejected_refresh_state`, `fast_factor_fail_reason`,
`fast_factor_fail_status`, `fast_block_restarts`,
`fast_kls_block_restarts`, `fast_kls_rebuild_restarts`,
`fast_kls_block_restart_last_row_pipeline`,
`fast_kls_block_restart_row_pipeline_count`,
`fast_kls_block_restart_last_row_pipeline_rows`,
`fast_kls_block_restart_last_row_pipeline_threads`,
`fast_kls_block_restart_last_row_pipeline_prefix_rows`,
`fast_kls_block_restart_last_row_pipeline_suffix_rows`,
`fast_kls_block_restart_last_row_pipeline_gap_rows`,
`fast_kls_block_restart_last_row_pipeline_etree_tail`,
`fast_kls_block_restart_row_pipeline_etree_tail_count`,
`fast_kls_block_restart_last_row_pipeline_etree_tail_rows`,
`fast_kls_block_restart_last_row_pipeline_etree_tail_gap_rows`,
`fast_kls_block_restart_last_row_pipeline_etree_tail_exact_mask`,
`fast_kls_block_restart_last_row_pipeline_etree_ready`,
`fast_kls_block_restart_row_pipeline_etree_ready_count`,
`fast_kls_block_restart_last_row_pipeline_etree_ready_rows`,
`fast_kls_block_restart_last_row_pipeline_etree_ready_threads`,
`fast_kls_block_restart_last_row_pipeline_etree_prefactor`,
`fast_kls_block_restart_row_pipeline_etree_prefactor_count`,
`fast_kls_block_restart_last_row_pipeline_etree_prefactor_rows`,
`fast_kls_block_restart_last_row_pipeline_etree_prefactor_threads`,
`fast_kls_block_restart_last_row_pipeline_etree_prefactor_wait_rows`,
`fast_kls_block_restart_last_row_pipeline_etree_prefactor_wait_deps`,
`fast_kls_block_restart_last_row_pipeline_separator_tail_scope`,
`fast_kls_block_restart_row_pipeline_separator_tail_scope_count`,
`fast_kls_block_restart_last_row_pipeline_separator_tail_scope_rows`,
`fast_kls_block_restart_last_row_pipeline_separator_queue`,
`fast_kls_block_restart_row_pipeline_separator_queue_count`,
`fast_kls_block_restart_last_row_pipeline_separator_private_rows`,
`fast_kls_block_restart_last_row_pipeline_separator_pipeline_rows`,
`fast_kls_block_restart_last_row_pipeline_separator_private_threads`,
`fast_kls_block_restart_last_row_pipeline_separator_partitioned`,
`fast_kls_block_restart_last_row_pipeline_separator_split_components`,
`fast_kls_block_restart_last_row_pipeline_pivot_tail_rows`,
`fast_kls_block_restart_last_row_pipeline_pivot_restarts`,
`fast_kls_block_restart_last_row_pipeline_supernode_update_groups`,
`fast_kls_block_restart_last_row_pipeline_supernode_update_rows`,
`fast_kls_block_restart_last_row_pipeline_supernode_panel_update_groups`,
`fast_kls_block_restart_last_row_pipeline_supernode_panel_update_rows`,
`fast_tail_restarts`,
`fast_repaired_last_offdiag_suffix_refresh`,
`fast_repaired_offdiag_suffix_refresh_count`, and
`fast_repaired_offdiag_full_refresh_count`, and
`fast_repaired_parallel_tail_blocks`. `fast_factor_fail_reason` is `0` when no
checked fast-factor failure was recorded; nonzero values distinguish row-wise U
validation invalidation (`1`), EGraph invalidation (`2`), mapped-refactor
invalidation (`3`), pool invalid or allocation failure (`4`), KLU refactor
failure (`5`), exhausted repair without acceptance (`6`), invalid kernel status
(`7`), singular status (`8`), and the dominant-BTF fast-repair guard (`9`).
`fast_factor_fail_status` keeps the corresponding KLU status
value. Together these fields describe the first
rejected factor-order pivot, its original matrix column, the rejecting row, the
L-multiplier or row-maximum-to-pivot ratio that violated a KLS-owned pivot
check when available, the accepted pivot magnitude and candidate or row-maximum
magnitude at the reject, the checked-refactor unfinished seed used before
ordered-ETree tail closure, the best row-tail
candidate that can be computed from current prefix state in the checked
row-major path, its retained row-tail position, whether that prefix-current
candidate satisfies the same pivot-tolerance predicate that rejected the
original reused pivot, the row chosen by the pivoting block or tail repair at
that pivot, whether it matches the retained tail candidate when one was
available, the first pivot
whose row changed in the fallback repair, how many changed pivots were before
and at/after the rejected pivot, whether the robust repair outcome preserved
the old block prefix, whether KLS can reconstruct the live KLU prefix
`P`/`Pinv`/pruning state needed by a local serial tail restart, the estimated
full-block repair work and validated non-root serial suffix restart work/saved
work, the suffix work and columns currently overcomputed relative to the
retained CKTSO-style pivoting-tail worklist, the
rejected BTF block, the suffix from that pivot to the end
of the block, the exact U-pattern descendant tail inside that block, the
row-refactor successor tail when row-major metadata is available, including
block-local BTF tails retained by the checked row-major path, the
ordered-block ETree successor path that a pivoting tail-restart upper-bound
scheduler would at least have to revisit, the sorted pivoting-tail worklist
scope seeded from the interrupted guessed-EGraph unfinished set before any
prefix refresh, with row-tail, serial block-suffix, and pool block-suffix seeds
kept as fallbacks, the
first and last rows in that worklist, whether it includes the rejected pivot,
whether the retained order is topologically safe for a future tail kernel to
consume, whether the failed pass left an unknown, prefix-current, or all-current
numeric state, the
number of repaired BTF blocks, how many repairs used the KLS-owned pivoted block
kernel, whether the last multi-thread KLS-owned block repair first ran through
the restartable row pipeline before the older serial/KLU fallbacks, the number
of those block-repair pipeline runs, their processed rows, active worker count,
how many prefix and suffix rows were preserved around an exact tail-envelope
pipeline, grouped supernode and cached-panel update rows used inside the
pipeline, whether row-first block repair used the retained pivoting-tail
first/last row range with active masks for non-contiguous gap rows, whether
that repair consumed the retained ETree-descendant tail as its compact
topological row worklist, whether
accepted row-first or serial tail repairs refreshed only the off-diagonal suffix
after proving the prefix unchanged, and pivot-tail serial restart work, and
the number of serial tail restarts actually executed, plus whether
repaired serial-tail restarts refreshed only the off-diagonal suffix or rebuilt
all off-diagonal entries, and whether a non-contiguous serial repair exactly
matched the retained ETree tail mask or had to recompute promoted internal gap
columns as overcompute. These fields are
intended to guide fuller CKTSO-style tail-restart work without accepting an
unsafe reused pivot order.

Use `--no-static-pivoting` to disable KLS's value-aware static row-pivoting
trial. When enabled, `auto` can preemptively build a weighted row permutation
for medium matrices whose input values show a mostly weak or missing diagonal,
or react after a high off-diagonal pivot count. The permutation moves large
entries onto the diagonal. The reactive post-factor trial is skipped for small,
low-work cases where the first factorization already succeeded and the diagonal
weakness is not severe enough for the static-match setup cost to pay back over
the default repeated-refactor SPICE-cycle model. The pre-factor partial-weak
medium gate also requires enough nonzeros per row to pay for matching and a
trial factorization, so sparse full-diagonal spike cases avoid a rejected setup
trial. For larger weak-diagonal
candidates, the matching augment uses a layered bipartite search so KLS can
complete many independent augmenting paths per pass instead of restarting a
search from each unmatched row. For cheap small candidates, KLS can instead run
an exact sparse
augmenting-path assignment on transformed log magnitudes, equivalent to a
maximum-product diagonal match when a full matching is found. When that exact
assignment is accepted, KLS derives row/column scales from the assignment dual
potentials so matched diagonal entries normalize to unit magnitude and
nonmatched entries are bounded by the same transformed-cost inequalities. The
exact path is gated by matrix order and `n * nnz` work because even the KLS-owned
assignment-specific path is not a substitute for production MC64 acceptance
and scaling on larger SPICE matrices. Medium static-pivot matches then run a
bounded alternating-cycle improvement pass that can accept profitable three-
and four-row exchanges missed by the pair-swap pass. For matched diagonals with
large numeric spread, KLS can also trial matching-derived row/column
equilibration. The equilibration first tries a
dual-potential scaling pass that mirrors MC64's diagonal-normalization
conditions, rejecting cases whose greedy matching leaves large positive-cycle
evidence, then falls back to the older heuristic balancing pass. KLS keeps the
transformed candidate in no-scale mode only when the factorization succeeds
with acceptable pivoting and conditioning, or when numeric fill/flop, pivoting,
and reciprocal-condition evidence improve. Medium-large static-match
candidates with both majority missing and majority weak diagonals can keep the
row permutation but prefer unscaled values, avoiding matching-equilibration
setup when it would increase fill. For medium weak-diagonal static-match
candidates whose accepted AMD symbolic leaves a many-block BTF with one
dominant block, KLS can factor a METIS nested-dissection candidate and keep it
only when actual numeric fill/flops improve materially without unacceptable
conditioning loss. This keeps CKTSO-style METIS plus static-pivoting wins
available without forcing METIS on all weak-diagonal circuits. Benchmark JSON
reports both whether static pivoting was enabled, whether KLS selected it, and
whether the accepted static match used exact assignment. It also reports
`selected_exact_matching_scaling` when the accepted KLS-owned exact assignment
retained its dual-derived scaling, and `selected_spral_matching` when the
accepted row permutation came from the
LGPL-compatible SPRAL Hungarian or auction path rather than KLS's in-tree
matcher. KLS does not vendor HSL MC64 or the MC64 copies carried by some solver
projects. When
`KLS_ENABLE_SPRAL_SCALING=ON`, KLS builds or links BSD-licensed SPRAL
Hungarian/auction matching as an LGPL-compatible MC64-adjacent component.
SPRAL same-cardinality Hungarian
matches are not installed blindly; KLS uses them before factorization only for
large weak-diagonal dominant-BTF candidates, or after the first factorization
for dense high-off-diagonal-pivot cases, and keeps them only when the accepted
numeric path passes KLS's pivoting, conditioning, fill, flop, and setup-value
checks. The post-factor trial runs after cheaper KLS policy trials so rejected
SPRAL candidates do not mask a simpler scale or pivot-tolerance fix.

Use `--no-btf` to measure the same ordering/scaling policy without KLU's BTF
decomposition. With `--ordering auto` and BTF enabled, KLS can still bypass BTF
for large matrices where BTF finds a single block, one dominant block with a
small fringe, or a many-block decomposition whose large block produces an
inflated symbolic estimate, and a no-BTF symbolic retry gives enough evidence to
keep or improve the fill estimate. METIS-started auto paths skip the single-block
retry and only use the stricter dominant/inflated-block retries, avoiding extra
symbolic work on already-good single-block circuit cases. Dominant and inflated
many-block retries require a known BTF symbolic score; when BTF's score is
unknown, KLS keeps the decomposition instead of accepting a misleading no-BTF
single-block estimate. Score-gated AMD/COLAMD single-block retries also cover
medium OPF-style cases below the larger-circuit floor when the no-BTF symbolic
is essentially no worse. The retry is skipped for low-work dominant-BTF patterns
where keeping the decomposition is already the cheaper SPICE-cycle choice. For
large nearly diagonal spike patterns, medium low-degree full-diagonal patterns,
medium sparse high-degree mostly diagonal patterns, and for METIS-started medium
dense-diagonal high-degree patterns, `auto` can also start without BTF to avoid
analysis work it is likely to discard; benchmark JSON reports both
`requested_btf` and selected `btf`.

For large fragmented BTF analyses with a moderate block count and no dominant
diagonal block, `--ordering auto` can also run a no-BTF symbolic retry and keep
it when the symbolic fill score is substantially lower. This handles cases
where BTF exposes many off-block entries and worsens repeated refactor
throughput even though the fixed-pivot no-BTF factor remains stable.

For low-work medium many-block analyses, `auto` can also retry without BTF when
the largest BTF block is moderate, the diagonal is nearly complete, symbolic
work is bounded, and the no-BTF symbolic keeps fill growth and estimated flop
growth within guarded limits. This targets small circuit and Rommes-style
power-grid matrices where BTF's thousands of blocks add repeated-refactor
overhead without enough numeric work to amortize the decomposition. Dense
dominant-block and weak-diagonal cases stay on the existing BTF/static-matching
paths.

Use `--ordering metis` to force METIS nested-dissection ordering. The default
`--ordering auto` can start with METIS for medium, bounded-degree structures
only when the diagonal is nearly complete, preserving mesh-style
nested-dissection wins without forcing sparse-missing-diagonal Rommes-style
cases away from AMD. It also starts with METIS for medium dense-diagonal
high-degree patterns where delayed promotion would otherwise pay for an
avoidable first AMD factorization, and for small and medium spiked
low-diagonal patterns that resemble TSOPF-style paper cases where the symbolic
estimate understates nested-dissection benefit. Otherwise it compares AMD and
COLAMD by symbolic fill estimate. For large high-work no-BTF single-block
analyses, `auto` can also try METIS and SCOTCH symbolic candidates before
numeric factorization and keep one when its symbolic fill score is clearly
lower. This avoids paying for an AMD numeric factorization only to promote to a
nested-dissection ordering afterward. Dense-diagonal high-degree METIS
starts ask METIS for two separator attempts, then refine the nested-dissection
rank order with CAMD inside coarse rank constraints. This preserves the
separator-first shape while letting minimum degree reduce local fill/flops on
post-layout style circuits. Large METIS orderings also use coarse rank-group
CAMD refinement, matching the CKTSO paper's nested-dissection plus constrained
minimum-degree structure more closely than raw METIS. Medium near-full
bounded-degree METIS starts also use two separator attempts with random
matching coarsening to reduce nested-dissection factor work on mesh-like sparse
diagonals. Large full or nearly-full diagonal METIS starts keep the selected
BTF symbolic directly instead of paying for a no-BTF retry that is rejected or
slower on this ASIC-style shape. When METIS is enabled, `auto` can also promote
large, expensive first numeric factorizations to METIS if the trial
factorization materially reduces actual numeric flop/fill cost. This keeps
METIS available for hard nested-dissection cases without paying its analysis
cost on small circuit matrices. A narrower post-factor promotion also covers
small BTF-dominant matrices whose first AMD/COLAMD factorization shows both
many off-diagonal pivots and high actual fill/flop growth; this catches
power-grid-style cases where the symbolic estimate alone understates the
benefit of nested dissection.

Two analyze-time shortcuts use staged capability contracts rather than a
benchmark-sized input box. A 65,536--262,144-row, almost-full-diagonal input
with one macroscopic column spike may propose AMMF. KLS retains that proposal
only when the resulting symbolic has a nearly spanning BTF core, bounded
fragmentation, balanced fill, and 512--2,048 estimated operations per row
below the current setup-work ceiling; rejection resumes the ordinary ordering
tournament. A separate 32,768--131,072-row sparse-diagonal proposal may skip
the rest of the tournament after AMD proves one of three realized regimes:
fragmented bounded work, one-block bounded work, or a bounded-degree exactly
symmetric high-work grid. This second contract also requires a substantial
core, balanced fill, and full or unknown structural rank. Rejected AMD
symbolics are freed before normal AUTO selection continues. Set
`KLS_DISABLE_MEDIUM_SPIKE_MINFILL_PATH=1` or
`KLS_DISABLE_AUTO_AMD_SHORTCUT=1` to restore the respective generic fallback.
The public stats/benchmark JSON expose proposal and post-symbolic verdicts as
`medium_spike_minfill_candidate`/`medium_spike_minfill_symbolic_eligible`,
`sparse_broad_column_amf_no_btf_candidate`/
`sparse_broad_column_amf_no_btf_symbolic_eligible`, and
`bounded_degree_amf_no_btf_candidate`/
`bounded_degree_amf_no_btf_symbolic_eligible`. This lets corpus scans separate
a cheap input-stage near miss from a representation that can actually own the
recurring lifecycle.

For large paper-style diagonal patterns, `auto` can start directly with METIS
when the structure is a very-low-degree full diagonal, a sparse full diagonal
with bounded but nontrivial row/column degree, a sparse-diagonal low-degree
matrix with no empty rows or columns, a near-full diagonal with a large dense
row/column spike, or a large nearly diagonal matrix in a sparse-spike or
dense-spike density band. The nearly diagonal spike class also starts without
BTF so KLS does not first pay for an AMD symbolic pass before keeping the same
METIS/no-BTF numeric path. These predicates are structural, not
matrix-name-based, and are deliberately narrow so unrelated IBM `dc`/`trans`
cases and low-work dominant-BTF cases stay on the cheaper AMD/COLAMD path.

Use `--ordering scotch` to force SCOTCH nested-dissection ordering. In `auto`,
KLS can also compare SCOTCH as a guarded symbolic candidate for large
single-block, high-estimated-work patterns and keep it only when its symbolic
fill/work score materially beats the current AMD/COLAMD/METIS candidate. This
keeps SCOTCH available as a paper-backed nested-dissection alternative without
making it a broad default. ParMETIS is not wired into KLS yet; it is a
distributed-memory MPI package and should be treated as a separate future
component if KLS grows an MPI/distributed solver path.

To fetch public SuiteSparse Matrix Collection matrices listed in the manifest:

```sh
python3 scripts/fetch_suitesparse.py --manifest bench/suitesparse_circuit_manifest.txt --out data/suitesparse
```

For paper-driven tuning, `bench/suitesparse_paper_manifest.txt` is the
traceable union of public SuiteSparse matrices named in the local KLU, NICSLU,
SubtreeLU, CKTSO reference papers, and the CKTSO ordering supplement. It
currently resolves to 110 public SuiteSparse matrices. The CKTSO supplement also
names a few non-public or renamed labels that are not in the current
SuiteSparse index; those are documented as comments in the manifest. The full
manifest includes very large matrices. Use `scripts/audit_paper_manifest.py`
after editing the manifest to check the paper-derived public-name coverage. The
medium manifest is the default practical corpus:

```sh
python3 scripts/fetch_suitesparse.py --manifest bench/suitesparse_paper_medium_manifest.txt --out data/suitesparse-paper-medium
```

For deeper tuning against the remaining paper matrices, use the large
supplement. It contains the 17 public SuiteSparse paper cases excluded by the
medium caps, including the largest KLU/NICSLU/SubtreeLU labels:

```sh
python3 scripts/fetch_suitesparse.py --manifest bench/suitesparse_paper_large_manifest.txt --out data/suitesparse-paper-large
```

The paper corpus is a compatibility suite, not an independent generalization
test. `scripts/build_generalization_manifests.py` deterministically builds two
additional 24-matrix suites from SuiteSparse metadata. It excludes every
SuiteSparse group represented in the paper corpus, partitions the remaining
groups before selecting matrices, admits at most one matrix per group, and
round-robins across size, density, and structural-symmetry strata. Thus the
development and holdout manifests are SuiteSparse-group-disjoint from each
other and from the benchmark-derived corpus. The collection group is used as
a conservative, reproducible proxy for matrix family:

```sh
python3 scripts/build_generalization_manifests.py
python3 scripts/fetch_suitesparse.py \
  --manifest bench/suitesparse_generalization_dev_manifest.txt \
  --out data/suitesparse-generalization-dev
python3 scripts/fetch_suitesparse.py \
  --manifest bench/suitesparse_generalization_holdout_manifest.txt \
  --out data/suitesparse-generalization-holdout
```

Use `bench/suitesparse_generalization_dev_manifest.txt` while changing
selectors. Freeze the policy and thresholds before running
`bench/suitesparse_generalization_holdout_manifest.txt`, and report failures
as well as successful timing rows. The holdout is procedural rather than
secret: repeatedly tuning against it turns it into another development set.
Use `--index path/to/ssstats.csv` to regenerate from a pinned SuiteSparse
metadata snapshot. Generated manifests record the index timestamp and SHA-256;
`python3 scripts/build_generalization_manifests.py --check` verifies that the
committed files still match the selected metadata and seed.

For a one-use, source-frozen test outside every previously committed manifest,
use `bench/suitesparse_external_validation_v1_manifest.txt` and follow
`docs/external_validation.md`.  That protocol excludes the complete
SuiteSparse group of all earlier manifest entries, runs rank-preserving,
entrywise, and localized updates at 1/4/8 threads in counterbalanced solver
order, reports H10/H100/H1000 sensitivity without dropping failures, and keeps
post-reveal route ablations separate from the confirmatory score.

Two additional frozen manifests target the remaining large-policy evidence
gaps rather than the general suite score:
`bench/suitesparse_remaining_policy_natural_probe_manifest.txt` and the
disjoint `bench/suitesparse_remaining_policy_natural_holdout_manifest.txt`.
Both record their SuiteSparse-index snapshot and premeasurement selection
rules. Run them symbolically without paying for unrelated numeric failures:

```sh
python3 scripts/run_bench_suite.py --kls-bench build/kls_bench \
  --matrix-dir data/suitesparse \
  --manifest bench/suitesparse_remaining_policy_natural_holdout_manifest.txt \
  --analyze-only --threads 8 --timeout 300 \
  --jsonl build/remaining_policy_holdout_analyze.jsonl
```

In analyze-only suite mode, `selection_seconds` and the summary geometric mean
refer to `analysis_seconds`; ordinary suite runs retain the modeled SPICE-cycle
selection metric.

The default split caps matrices at 300,000 rows and one million entries so it
stays practical for routine development. A second deterministic 12+12 tier
extends coverage to 131,072--1,048,576 rows and up to 8,388,608 entries. Its
development and holdout groups are disjoint from each other and from every
paper-suite group, and CTest audits that boundary without network access:

```sh
python3 scripts/fetch_suitesparse.py \
  --manifest bench/suitesparse_generalization_large_dev_manifest.txt \
  --out data/suitesparse-generalization-large-dev
python3 scripts/fetch_suitesparse.py \
  --manifest bench/suitesparse_generalization_large_holdout_manifest.txt \
  --out data/suitesparse-generalization-large-holdout
```

Regenerate that tier from a pinned `ssstats.csv` with
`scripts/build_generalization_manifests.py`, counts of 12, row bounds
131,072--1,048,576, nonzero bounds 393,216--8,388,608, and explicit large
development/holdout output paths.

For SubtreeLU-specific tuning, `bench/suitesparse_subtreelu_manifest.txt`
contains the exact 46 public SuiteSparse circuit labels from SubtreeLU Fig. 7.
Fetch it into its own directory when running the full set, so the suite runners
see one matrix per basename:

```sh
python3 scripts/fetch_suitesparse.py --manifest bench/suitesparse_subtreelu_manifest.txt --out data/suitesparse-subtreelu
python3 scripts/audit_paper_manifest.py --manifest bench/suitesparse_subtreelu_manifest.txt --source SubtreeLU
```

Use the medium paper suite as the compatibility loop and the group-disjoint
development suite for solver-policy changes. Use the large paper supplement
as an overnight compatibility gate with explicit timeouts; a retained change
should improve a defensible structural class, not just one matrix name:

```sh
python3 scripts/run_bench_suite.py --kls-bench build/kls_bench --matrix-dir data/suitesparse-paper-large --orientation auto --threads 4 --timeout 900 --jsonl build/kls_paper_large.jsonl
```

For focused CKTSO-gap work, `bench/suitesparse_cktso_gap_manifest.txt` records
the current hard public medium rows from the latest KLS/CKTSO decomposition.
It is a tuning loop, not a replacement for the full medium manifest:

```sh
python3 scripts/run_bench_suite.py --kls-bench build/kls_bench \
  --matrix-dir data/suitesparse-paper-medium \
  --manifest bench/suitesparse_cktso_gap_manifest.txt \
  --orientation auto --threads 4 --repeat 1 --refactor-repeat 3 \
  --require-spral-scaling \
  --timeout 120 --jsonl build/kls_cktso_gap_focus.jsonl
```

For a shorter large-case reconnaissance before an overnight run, use the
selected large manifest. It covers the high-signal large paper cases that have
already shown KLS/CKTSO differences under a 120s per-process cap. Matrices
where CKTSO also times out under that cap are kept out of this tuning loop:

```sh
python3 scripts/run_bench_suite.py --kls-bench build/kls_bench --matrix-dir data/suitesparse-paper-large --manifest bench/suitesparse_paper_large_recon_manifest.txt --orientation auto --threads 4 --repeat 1 --refactor-repeat 1 --timeout 120 --jsonl build/kls_paper_large_recon.jsonl
python3 scripts/run_cktso_suite.py --cktso-compare build-cktso/cktso_compare --matrix-dir data/suitesparse-paper-large --manifest bench/suitesparse_paper_large_recon_manifest.txt --threads 4 --repeat 1 --refactor-repeat 1 --timeout 120 --jsonl build/cktso_paper_large_recon.jsonl
python3 scripts/run_klu2_suite.py --klu2-compare build-klu2/klu2_compare --matrix-dir data/suitesparse-paper-large --manifest bench/suitesparse_paper_large_recon_manifest.txt --repeat 1 --refactor-repeat 1 --timeout 120 --jsonl build/klu2_paper_large_recon.jsonl
```

Use `bench/suitesparse_stress_timeout_manifest.txt` for shared-hard timeout
cases such as `Hamrle3`. They are useful robustness and scalability checks, but
they should not drive CKTSO-relative tuning or geomean gap claims unless the
reference solver also produces a finite timing under the same cap.

To regenerate a metadata-bounded subset from the full paper corpus and record
the exact canonical SuiteSparse names used:

```sh
python3 scripts/fetch_suitesparse.py --manifest bench/suitesparse_paper_manifest.txt --out data/suitesparse-paper-medium --max-rows 700000 --max-cols 700000 --max-nnz 2500000 --write-resolved-manifest build/suitesparse_paper_medium_resolved.txt
```

To run every downloaded matrix and compute the SPICE-cycle geometric mean:

```sh
python3 scripts/run_bench_suite.py --kls-bench build/kls_bench --matrix-dir data/suitesparse --orientation auto --jsonl build/kls_suite.jsonl
```

Add `--passes N` to run each matrix multiple times and record the median
SPICE-cycle sample, which is useful when comparing small solver-policy changes.
Use an odd `N` when you need an exact median sample.
For broad paper suites, add `--timeout SECONDS` so one pathological matrix is
recorded as a failure instead of blocking the rest of the run.
When `--jsonl` is used, successful rows are written to that file and timeout or
process-failure records are written to the adjacent `.failures` sidecar. For
example, `build/kls_paper_large_recon.jsonl` records completed matrices and
`build/kls_paper_large_recon.failures` records timed-out or failed matrices.
`scripts/compare_bench_runs.py` reads those sidecars when present so large-case
comparisons show both timing ratios and missing rows caused by failures. By
default it still computes geomeans only over successful rows common to both
runs. Use `--include-manifest` and `--exclude-manifest` to keep the scored set
aligned with the intended tuning slice, for example excluding
`bench/suitesparse_stress_timeout_manifest.txt` from CKTSO-relative geomeans.
Add `--include-failures --failure-seconds SECONDS` when a hard-suite comparison
should score failed or missing rows with an explicit cycle-time penalty; using a
process timeout value as the penalty is only a lower bound for SPICE-cycle
comparisons.
For timeout diagnosis, `run_bench_suite.py --failure-diagnostics trace` reruns
each failed KLS matrix under a shorter row-pipeline trace cap, writes the trace
stderr under `<jsonl-stem>.failure-traces/`, and embeds the
`scripts/summarize_row_pipeline_trace.py` counters in the `.failures` record.
Use `--failure-diagnostics all` to keep the default analyze-only row alongside
that trace summary. For first-factor owner investigations, add
`--failure-trace-kls-first-factor on` so only the diagnostic run forces the
KLS-owned row pipeline. This is the reproducible way to recheck cases such as
`pre2`: analysis/order remains visible, while the failure sidecar also records
the scalar-output versus producer-target surface that drives the current
numeric-owner gap.

The suite metric is:

```text
analysis + initial_factor + final_state_solve
  + (first_refactor + its_paired_solve)
  + 98 * (steady_refactor + its_paired_solve)
```

The refactor-following solves are measured with their corresponding changed
numeric state rather than inferred from the final-state solve average.  JSON
records expose `refactor_solve_first_seconds` and
`refactor_solve_steady_seconds_avg`; this also charges any solve-side accuracy
recovery to the generation that required it.  Set
`BENCH_VERIFY_EACH_REFACTOR=1` to make every solver wrapper independently
check every changed-value generation and emit
`refactor_max_relative_residual`; the legacy
`KLS_BENCH_VERIFY_EACH_REFACTOR` spelling enables the same audit.

Compare two JSONL runs by matrix basename:

```sh
python3 scripts/compare_bench_runs.py --candidate build/kls_suite.jsonl --candidate-name kls-auto --reference build/klu_defaults.jsonl --reference-name klu-defaults
```

For timeout-heavy large reconnaissance runs, include sidecar failures with a
deliberate penalty:

```sh
python3 scripts/compare_bench_runs.py --candidate build/kls_paper_large_recon.jsonl --candidate-name kls-auto --reference build/cktso_paper_large_recon.jsonl --reference-name cktso --include-failures --failure-seconds 1000
```

To see which phase explains a solver gap, decompose the same JSONL pair into
analysis, initial factorization, repeated refactorization, and repeated solve
contributions. For KLS candidate rows, the report also includes EGraph
dependency levels, root/leaf/max-fanout counts, cluster levels, pipeline
columns, max per-column work, pipeline max per-column work, and
dependency-work estimates when those fields are present in the benchmark JSONL.
For row-engine experiments, `kls_bench` and `run_bench_suite.py` accept
`--row-refactor env|off|refactor|checked|all` and
`--kls-first-factor env|off|on`; they also accept
`--row-solve env|off|on` for the ordinary-factor row-solve seed gate and
`--input-index auto|32|64` for reproducing benchmark input-width choices. The
emitted `initial_factor_path`, `last_factor_path`, and `row_refactor_last_*`
fields show whether the first factorization, later numeric passes, and solves
really used KLS-owned paths. `--kls-first-factor env` preserves the library's
automatic first-factor decision, `off` disables it, and `on` forces the
KLS-owned attempt. Explicit row-refactor modes force the row engine, and
`--row-refactor off` also disables the automatic KLS-first row-refactor handoff
for reproducible column/EGraph baselines. With `env`, the automatic KLS-first
path uses the retained row/EGraph work estimates to skip row refactors whose
static work model is already worse than the existing exact EGraph schedule:

```sh
python3 scripts/decompose_solver_gap.py --candidate build/kls_suite.jsonl --candidate-name kls-auto --reference build/cktso_suite.jsonl --reference-name cktso
```

Add `--concise` for the usual solved-by-reference triage view: it keeps the
cycle ratio, dominant phase, paper-gap signal, path names, and the core EGraph
pipeline work counters while omitting the full diagnostic surface. When
producer-step advance-batch stats are present, it also keeps the duplicate
retained-state row count, exact unique row count, and row-collapse ratio so the
grouped-current owner opportunity remains visible in the short report. Use
`--include-manifest bench/suitesparse_cktso_gap_manifest.txt` with either
comparison script to rank only the current CKTSO-gap focus set, and
`--exclude-manifest bench/suitesparse_stress_timeout_manifest.txt` to keep
shared-hard stress rows such as `Hamrle3` out of CKTSO-relative tuning reports.
Manifest entries may be bare SuiteSparse names or `.mtx` basenames.

To summarize the CKTSO-style tail-restart opportunity fields across a KLS JSONL
run, use:

```sh
python3 scripts/summarize_tail_restart_opportunities.py --jsonl build/kls_suite.jsonl
```

For executed serial suffix restarts, the summary also reports how many columns
and how much work the suffix path does beyond the retained CKTSO-style
pivoting-tail plan. Large overcompute there marks cases where the full
pipelined pivoting-tail kernel should matter most.

The same script accepts `--jsonl -` for a single piped `kls_bench --json` row
when inspecting a focused fast-reject case.

To generate those diagnostics across a manifest, `run_bench_suite.py` forwards
the deterministic diagonal-stress controls accepted by `kls_bench`:

```sh
python3 scripts/run_bench_suite.py --kls-bench build/kls_bench \
  --matrix-dir data/suitesparse-paper-medium \
  --manifest bench/suitesparse_cktso_gap_manifest.txt \
  --orientation auto --threads 4 --repeat 1 --refactor-repeat 0 \
  --require-spral-scaling \
  --stress-diagonal-scale 1e-9 \
  --timeout 120 --jsonl build/kls_tail_stress_gap.jsonl

python3 scripts/summarize_tail_restart_opportunities.py \
  --jsonl build/kls_tail_stress_gap.jsonl
```

For long manifests, the KLS, CKTSO, and KLU2 suite runners preserve manifest
order and accept `--skip N --limit M`, so large or CKTSO-gap suites can be run
in reproducible chunks without treating unrun rows as solver failures in later
comparisons.

An optional CKTSO comparison tool can be built when you provide a local CKTSO
distribution:

```sh
cmake -S . -B build-cktso -DKLS_BUILD_CKTSO_COMPARE=ON -DCKTSO_ROOT=/path/to/cktso
cmake --build build-cktso -j --target cktso_compare
LD_LIBRARY_PATH=/path/to/cktso/rocky8_x64_gcc850 ./build-cktso/cktso_compare matrix.mtx 16 10 10
```

The CKTSO harness emits the same `spice_cycle_seconds` formula as `kls_bench`,
so JSON rows can be compared directly on the repeated-SPICE metric.

Run a CKTSO JSONL suite with:

```sh
python3 scripts/run_cktso_suite.py --cktso-compare build-cktso/cktso_compare --matrix-dir data/suitesparse --threads 1 --jsonl build/cktso_suite.jsonl
```

CKTSO must be licensed correctly according to its own distribution
requirements, usually by colocating the license file with the selected shared
library.

An optional Trilinos KLU2 comparison tool can be built against a local Trilinos
checkout. It calls the Amesos2 KLU2 headers directly and emits the same
SPICE-cycle JSON metric:

```sh
cmake -S . -B build-klu2 -DKLS_BUILD_KLU2_COMPARE=ON \
  -DKLU2_ROOT=/path/to/Trilinos/packages/amesos2/src/KLU2 \
  -DTEUCHOS_CORE_ROOT=/path/to/Trilinos/packages/rol/src/compatibility/teuchos-lite
cmake --build build-klu2 -j --target klu2_compare
./build-klu2/klu2_compare matrix.mtx 10 10
python3 scripts/run_klu2_suite.py --klu2-compare build-klu2/klu2_compare --matrix-dir data/suitesparse-paper-medium --jsonl build/klu2_suite.jsonl --timeout 120
```

## Status

This is a functional implementation with KLS-level analysis choices for
repeated SPICE-style solves and a KLS-owned threaded refactor path for BTF block
parallelism on a narrow class of large cases. KLS also has a precomputed
single-block and serial BTF refactor scatter path for unscaled repeated
refactors, experimental row-major row refactors for single-block factors and
BTF diagonal blocks, KLU row-scaled row refactors, and a narrow scaled
dominant-BTF subset, plus a KLS-owned block-local pivot restart for fast-factor
failures including root-of-block rejects. Benchmark stats also report
row-major U-pattern supernode candidates and detailed rejected-row/multiplier
coordinates from KLS-owned pivot checks, plus row-refactor cluster/pipeline
counters, compact dense-panel markers, and last-run markers, so the remaining
SubtreeLU/CKTSO row-segment work can be evaluated on the same slow-case
artifacts.
On the refreshed selected-large reconstruction, KLS is ahead of the saved KLU2
artifact but still trails the saved CKTSO artifact, with `pre2` still timing out.
The June 29, 2026 current-source rerun keeps the same shape:
`build/kls_current_large_recon_t4_r1_ref1_timeout120.jsonl` completes six of
eight selected large rows, wins `TSOPF_FS_b39_c30`, times out on `pre2` and
`Hamrle3`, and scores `1.06x` slower than the fresh CKTSO large-recon artifact
when both solver failures are charged as 1000s. Same-session `pre2` probes show
analysis completes quickly but
factor-only default AMD, forced KLS-first, transpose AMD, METIS, and
no-static-pivoting all exceed 120s, so the unresolved `pre2` gap is cold
first-factor numeric machinery rather than input width, ordering, orientation,
or repeated-refactor policy. A follow-up CBLAS-enabled build with
`KLS_ENABLE_CBLAS_SUPERNODE=1` did not change that conclusion: the matching
`pre2` forced KLS-first run still timed out at 120s with no JSON row, while
same-option `ASIC_680k` checks reported zero CBLAS update counters. The
existing BLAS gates are therefore not the current slow-case blocker.
The same current-source large run also rechecked the completed-but-slow
`nxp1`/`rajat30` pair. Forcing row refactor measured `151.3s`/`131.7s`, full
EGraph supernode updates measured `48.2s`/`40.6s`, and the consumer-plan
executor measured `49.3s`/`36.3s`, versus the default `36.2s`/`37.3s`. These
large-row probes reject bypassing the row lower-bound gate or promoting the
current scalar supernode executors; the needed paper-level work remains a
production row/supernode numeric representation, plus the separate `pre2`
first-factor bottleneck.
A current CBLAS-capable top-five CKTSO-gap check keeps that conclusion: the
same binary measured `1.44563s` geomean with `KLS_ENABLE_CBLAS_SUPERNODE=0`
and `1.49419s` with `=1`, while all focus rows reported zero CBLAS update
counters.
The latest `pre2` stack samples instead show the current KLS-first BTF-parallel
route leaving one worker to factor the dominant block while the other BTF
workers exit; the no-BTF intra-block route keeps workers alive but spends the
sampled active epoch rebuilding row-first panel-cache state under the pipeline.
KLS now grows the row-first panel cache with allocator-preserved `realloc`
rather than KLS-owned allocate/zero/copy/free loops, reducing avoidable storage
growth overhead on completed large paths, but `pre2` still exceeds the 120s
factor-only cap. A current CBLAS-capable rerun again kept the BLAS hypothesis
bounded: explicit `KLS_ENABLE_CBLAS_SUPERNODE=0` and `=1` `ASIC_680k` checks
both reported zero CBLAS update counters, and forced `pre2` with the runtime
gate on still timed out at 120s with no JSON row. KLS also removes a duplicate
row-first supernode reset before prefix panel-cache rebuilds; the rebuild
helper already performs the reset, so the call-site scan was redundant.
KLS now also reports row-first panel-cache build/append volume in stats and
benchmark JSON. A forced KLS-first `transient` probe reported `3,814`
appended panels and `1,075,649` stored panel-cache entries alongside `97,678`
cached panel update groups, while the forced KLS-first `ASIC_680k` stress
probe reported zero panel-cache activity. The matching forced `pre2` factor
probe still timed out at 120s, so that case remains unresolved, but future
slow-path runs can now separate panel-cache staging volume from cached-panel
consumption.
KLS also reports dominant-BTF row-up-looking first-factor coverage through
`kls_first_last_dominant_btf_pipeline`,
`kls_first_last_dominant_btf_pipeline_block`,
`kls_first_last_dominant_btf_pipeline_rows`, and
`kls_first_last_dominant_btf_pipeline_has_separator`; setting
`KLS_TRACE_KLS_FIRST_FACTOR=1` prints the same decision before a long factor
run can time out. Setting `KLS_TRACE_ROW_PIPELINE=1` additionally prints
row-first pipeline progress with committed rows, scalar dependency applications,
published-U entries scanned by the scalar update loop, scalar supernode fallback
rows, and local/shared row-entry reserve growth/copy volume. Setting
`KLS_TRACE_ROW_PIPELINE_LONG_ROW_ENTRIES=N` with row-pipeline tracing enabled
also prints the first live dependency and the commit summary for rows whose
scalar published-U scan crosses `N` entries. Setting
`KLS_TRACE_ROW_PIPELINE_COMPACT_WINDOW=1` adds a trace-only symbolic live-window
probe for current/future rows. `KLS_ROW_PIPELINE_COMPACT_WINDOW=<slots>`
controls the number of compact sparse states, and
`KLS_ROW_PIPELINE_COMPACT_WINDOW_MAX_ENTRIES=<entries>` caps one state.
`KLS_ENABLE_ROW_PIPELINE_COMPACT_EXEC=1` turns the same window into an opt-in
sparse numeric prototype for unreserved first-factor rows: producer publication
collects compact states whose root dependency is the just-published producer,
streams that producer row across the target batch, and lets a worker claim one
prepared state instead of replaying its prefix from the input row. Live compact
states are bucketed by root dependency, so producer publication no longer scans
the whole compact window before every grouped update. Each compact state also
keeps an internal row-to-position index, so wider compact-window probes can
append new fill entries without sorted-array insertion. When a persistent dense
compact group is active, later compact rows whose next root matches that group
can be merged into it instead of scattering immediately back to independent
sparse states; row-pipeline traces report
`compact_window_group_merge_{targets,cols,values}` for that surface. Dense
compact union materialization is guarded by a sparse-work envelope, and traces
report `compact_window_union_skip_{targets,cols,values,sparse_values}` when the
exact sparse batched update is used instead. Compact execution traces also
report
`compact_window_claim_{attempts,claims,misses,stale_clears,group_scatters,group_detaches}`
so reuse can be separated from eager-update cost.
`KLS_ENABLE_ROW_PIPELINE_COMPACT_SPARSE_GROUP=1` keeps an opt-in persistent
compact group in each member state's sparse row/index storage instead of the
dense target-by-union matrix. This is useful for paper-gap probes that need to
separate grouped-owner coverage from dense materialization cost; it remains off
by default and still requires `KLS_ENABLE_ROW_PIPELINE_COMPACT_EXEC=1`. It can
also be combined with delayed output, where grouped producer updates keep only
prefactor/internal columns and leave target-output columns for claim-time
	replay. Sparse grouped claims detach only the claimed row, so a current-row
	claim does not scatter the rest of the grouped owner.
	Sparse grouped-claim traces also report
	`compact_window_claim_group_run_{probes,states,consecutive,near64,near512,max_consecutive}`
	to separate strict commit-adjacent run length from the wider local grouped
	surface. On `pre2` with a 512-row compact window, the active grouped surface
	averaged `31.6` states per claim probe and `12.4` states within the next 64
	positions, but only `4.9` strictly consecutive states. This rejects a purely
	contiguous commit-owner as too narrow and points the next output-owner work at
	a wider local panel or a selective claim-run replay.
	`KLS_ENABLE_ROW_PIPELINE_COMPACT_DELAY_OUTPUT=1` adds an opt-in
	prefactor/postfactor split for that compact executor: producer publication
	updates only dependencies below the target row, while diagonal/output U updates
are replayed directly into the worker row workspace when the prepared compact
state is claimed. Traces report
`compact_window_delayed_output_{skips,replays,deps,entries,scan_entries,seek_skips}`.
With sparse compact groups enabled, traces also report
`compact_window_delayed_group_replay_{surfaces,states,deps,unique_deps,duplicate_deps,scan_entries,group_scan_entries}`
to measure how much claim-time delayed replay could share producer-row scans
across grouped states.
`KLS_ENABLE_ROW_PIPELINE_COMPACT_DEFERRED_OUTPUT_CACHE=1` adds an experimental
per-state postfactor cache for delayed output. When a compact state streams a
producer U row, output-column deltas can be accumulated into a separate sparse
map and applied at claim time before the normal delayed-output replay tail.
`KLS_ROW_PIPELINE_COMPACT_DEFERRED_OUTPUT_MAX_ENTRIES=<entries>` caps one
state's cached output map (`65536` by default), and
`KLS_ROW_PIPELINE_COMPACT_DEFERRED_OUTPUT_MIN_U_ENTRIES=<entries>` controls the
minimum producer U-row width that starts a cached tail (`1024` by default);
once a tail has started, later producer rows continue to cache so the covered L
range remains contiguous. This path remains off by default because current
probes show that producer-time hash aggregation can cost more than it saves:
on `G2_circuit`, no-trace factor time moved from `4.77s` to `5.49s` with a
128-entry start threshold, while preserving a `3.4e-16` relative residual. On a
60s traced `pre2` probe, the same threshold kept progress at `589824/629628`
rows and cut delayed-output replay entries from about `388M` to `304M`, but did
not advance past the timeout checkpoint. Traces report
`compact_window_deferred_output_cache_{stores,applies,l_entries,entries,disables,overflows}`
plus cache-entry and store-share ratios in
`scripts/summarize_row_pipeline_trace.py`.
`KLS_ENABLE_ROW_PIPELINE_COMPACT_GROUP_OUTPUT_CACHE=1` adds a second
experimental delayed-output cache for sparse compact groups. Instead of giving
each future row its own producer-time hash map, the sparse group owns one
row/column matrix of output deltas and materializes the claimed row into the
normal deferred-output replay path only when that row leaves the group.
`KLS_ROW_PIPELINE_COMPACT_GROUP_OUTPUT_MAX_SLOTS=<slots>` caps the grouped
matrix (`65536` by default), and
`KLS_ROW_PIPELINE_COMPACT_GROUP_OUTPUT_MIN_U_ENTRIES=<entries>` controls the
minimum producer U-row width that starts grouped output caching (`128` by
default). The existing claim-run/span-owner paths clear this cache and proceed
through their faster worker-owned replay paths, so this option remains a
conservative, off-by-default experiment. Current probes show correctness but
not a CKTSO-gap closer: `G2_circuit` no-trace factor time stayed in the same
range (`4.78s`, `3.4e-16` relative residual), while a 60s traced `pre2` probe
again reached `589824/629628` rows and reduced delayed-output replay entries
from the earlier roughly `387.8M` to `363.1M`. Traces report
`compact_window_group_output_cache_{stores,applies,l_entries,entries,cols,disables,overflows}`
plus grouped-cache entries/apply, columns/apply, and store-share ratios.
Claim-time sparse grouped-output probes additionally report
`compact_window_claim_group_output_{surfaces,states,deps,unique_deps,duplicate_deps,scan_entries,group_scan_entries}`.
These counters do not replay or materialize any future-row output; they estimate
the scan work a persistent grouped row/panel owner could share at each grouped
claim. On `pre2` with a 512-row compact window,
`build/kls_pre2_claim_group_output_surface_w512_trace45.stderr` timed out at
`510403/629628` rows under trace overhead, but measured `63,297` claim surfaces
with `31.6` states per surface. Only `5.95%` of delayed-output dependencies were
unique, and the grouped-scan estimate was `3.52%` of the per-state scan volume
(`662,930,248` versus `18,833,835,429` entries). This keeps the next CKTSO and
SubtreeLU-shaped implementation target on a streaming grouped owner that shares
producer-row scans without storing output into each future row.
The same grouped-claim probe now also reports
`compact_window_claim_group_run_{stealable,max_stealable,unreserved_near64}`.
These counters check whether strict commit-adjacent group rows have already been
reserved by other workers. In
`build/kls_pre2_claim_group_stealable_w512_trace45.stderr`, `pre2` reached
`589824/629628` rows and showed `3.86` stealable rows per claim probe, or
`79.4%` of the strict consecutive run, with a maximum stealable run of `121`.
Within the next 64 grouped positions, `91.9%` were still unreserved. That makes
a local claim-run owner feasible from a scheduler perspective; the remaining
constraint is to apply the shared output in a worker-owned scratch/run executor
instead of materializing it into future compact states.
`KLS_ENABLE_ROW_PIPELINE_COMPACT_CLAIM_RUN=1` adds an opt-in strict sparse-group
claim-run owner. It atomically reserves adjacent grouped positions when the
claim cursor reaches the first row, keeps the later rows in worker-owned compact
scratch, and consumes those reserved positions before taking new atomic work.
`KLS_ROW_PIPELINE_COMPACT_CLAIM_RUN_MAX_ROWS=<rows>` caps the run length
(`8` by default, clipped at `128`). Traces report
`compact_window_claim_run_{reservations,rows,recomputes,updates,update_targets}`
and the non-materializing owned-run delayed-output surface as
`compact_window_claim_run_output_{surfaces,states,deps,unique_deps,duplicate_deps,scan_entries,group_scan_entries}`.
On `pre2` with the same 512-row compact window, the first claim-run trace
reserved `19,031` runs and `69,727` rows (`3.66` rows per reservation) but
performed no owned-run producer updates and still timed out at `589824/629628`
rows. A follow-up ready-root catch-up pass performed `91,546` owned-run updates
over `255,424` targets, but it also stopped at `589824/629628` rows with
`407,796,098` delayed-output scan entries. The owned-run output-surface trace,
`build/kls_pre2_claim_run_output_surface_w512_trace45.stderr`, reached the same
checkpoint and measured `19,037` owned surfaces over `69,741` states. Strict
runs had real sharing (`965,697` unique dependencies out of `2,901,084`,
`33.3%`) and would cut scan volume to `86,861,893` grouped entries from
`324,714,822` per-state entries (`26.8%`), but this is much weaker than the
wider sparse-group surface (`3.7%` grouped-scan ratio in the same trace). This
keeps the claim-run path diagnostic and off by default: adjacent-row ownership
alone is not the missing CKTSO/SubtreeLU trailing-output workspace; the next
owner needs a wider local row/panel surface without storing output back into
future compact states.

`KLS_ENABLE_ROW_PIPELINE_COMPACT_CLAIM_STREAM_OWNER=1` adds a stricter
claim-run streaming experiment. While a strict claim run is reserved, newly
published producer rows can stream delayed-output contributions once into a
worker-owned panel and advance the claimed rows' delayed-output replay cutoff
when they are later consumed. The owner records each row's stream-start L-count
and replays any older delayed-output prefix before applying the streamed panel,
so the cutoff only advances across the streamed range.
`KLS_ROW_PIPELINE_COMPACT_CLAIM_STREAM_OWNER_MAX_SLOTS=<slots>` caps that panel
(`64` by default); lower values disable the owner instead of entering a
partial-stream fallback. This remains diagnostic: a focused `pre2` run with a
two-row claim run and 64 stream-owner slots,
`build/kls_pre2_claim_stream_owner_r2_s64_trace45.stderr`, did not reach the
first progress checkpoint in 45s, and the matching claim-run-only control also
did not reach it. The strict adjacent claim-run family is therefore not the
right `pre2` gap closer even when output replay is streamed incrementally.

The follow-up local-span probe records the unreserved grouped rows within the
next 64 factor positions as
`compact_window_claim_span_output_{surfaces,states,deps,unique_deps,duplicate_deps,scan_entries,group_scan_entries}`.
It also records the numeric owner shape as
`compact_window_claim_span_output_{value_entries,state_col_slots,unique_cols,panel_slots,dense_slots,max_states,max_unique_cols}`:
`value_entries` counts delayed-output additions, `state_col_slots` counts
distinct touched row/column cells, `unique_cols` counts the local panel column
surface, `panel_slots` counts the triangular row/panel cells with
`col >= row`, and `dense_slots` counts the full rectangular
`states * unique_cols` workspace.
`KLS_ROW_PIPELINE_COMPACT_CLAIM_SPAN_TRACE_ROWS=<rows>` changes that
trace-only local span width; `0` disables the span-output surface, and values
above the compact-window slot cap are clipped.
On a short `bcircuit` trace with a 16-row span, the owner-shape counters
recorded `15,974,574` delayed-output additions, `1,393,652` distinct
state/column cells, `1,839,073` triangular panel slots, and `2,183,001`
rectangular dense slots. That makes the next direct prototype a sparse or
masked local span/panel workspace, not just a scheduler reservation.
A 45s `pre2` shape trace with a 256-row span,
`build/kls_pre2_claim_span256_shape_trace45.stderr`, timed out earlier
(`262144/629628`) because the diagnostic is intentionally heavy, but its
workspace ratios are the key result: `8,251,520,214` delayed-output additions
collapsed to `119,231,230` distinct state/column cells, `137,503,076`
triangular panel slots, and `150,801,490` full dense slots. That is `69.2`
additions per touched cell and only `1.26x` dense slots per touched cell, so
the paper-aligned next step is a worker-owned masked/dense span workspace that
streams grouped delayed output into the local panel rather than materializing
it into each future compact state.
`KLS_ENABLE_ROW_PIPELINE_COMPACT_SPAN_OWNER=1` enables the first opt-in
implementation of that worker-owned span workspace. It reserves a local
contiguous factor-position span, owns the grouped compact rows inside that
span, streams currently-known delayed output into a worker-owned masked dense
panel, and applies each panel row only when that same worker consumes the
reserved row. `KLS_ROW_PIPELINE_COMPACT_SPAN_OWNER_ROWS=<rows>` controls the
span width (`4` by default, clipped at `256`). On `bcircuit`, the same
compact/delay setup improved from `0.625s` initial factor without the span owner
to `0.450s` at 16 rows and `0.384s` at 64 rows in quick 3-repeat probes; 256
rows was slightly worse at `0.392s`. The same prototype still timed out on
`pre2` after 45s, so it is retained as experimental paper-aligned substrate,
not a default policy.
`KLS_ROW_PIPELINE_COMPACT_SPAN_OWNER_MAX_SLOTS=<slots>` caps the worker-owned
masked panel (`65536` slots by default). The setup phase is capped separately:
`KLS_ROW_PIPELINE_COMPACT_SPAN_OWNER_MAX_LINKS=<links>` defaults to `1`
dependency links and
`KLS_ROW_PIPELINE_COMPACT_SPAN_OWNER_MAX_SCAN_ENTRIES=<entries>` defaults to
`524288` scanned U-row entries.
`KLS_ROW_PIPELINE_COMPACT_SPAN_OWNER_MIN_ROWS=<rows>` defaults to `4`, and
`KLS_ROW_PIPELINE_COMPACT_SPAN_OWNER_MIN_ENTRIES_PER_SCAN=<ratio>` defaults to
`2`, so tiny panels are rejected unless they own enough rows and are projected
to replay at least twice as many output entries as the grouped setup scan.
Candidates that exceed any cap or miss this payoff gate fall back to the
existing compact path before allocating or filling the owner panel. By default
the owner also requires every row in the reserved span to be an owned grouped
state, avoiding reserved holes that serialize recomputation inside the same
worker. `KLS_ROW_PIPELINE_COMPACT_SPAN_OWNER_ALLOW_HOLES=1` restores the earlier
experimental behavior. `KLS_ROW_PIPELINE_COMPACT_SPAN_OWNER_ALLOW_PREFIX=1`
accepts only the contiguous owned prefix of a holey candidate; it is also
experimental and remains off because it regressed `pre2`.
`KLS_ENABLE_ROW_PIPELINE_COMPACT_SPAN_SHARED_OWNER=1` enables a more
paper-shaped non-contiguous shared owner workspace: a grouped sparse claim can
build one shared delayed-output panel for nearby grouped states, tag those
states, and apply the panel only when each state is later claimed instead of
materializing output columns into future compact rows. This is opt-in because
the current implementation still does not close the `pre2` gap.
`KLS_ROW_PIPELINE_COMPACT_SPAN_SHARED_OWNER_SCAN_ROWS=<rows>` widens the
positional scan window used to choose non-contiguous shared-owner rows while
leaving `KLS_ROW_PIPELINE_COMPACT_SPAN_OWNER_ROWS` as the bounded owner row
capacity. The default keeps the old behavior by scanning only the owner-row
span.
`KLS_ROW_PIPELINE_COMPACT_SPAN_SHARED_OWNER_DEP_LIMIT=N` narrows that shared
owner to a capped delayed-output producer window. It picks up to `N` common
pending producer rows from nearby grouped states, streams only those producers
into the shared panel, and advances each tagged row through the contiguous
prefix covered by the selected producer set when claimed. The selected producer
list is now resizable and sorted, so `N` is no longer clipped at the old
64-producer probe cap; the usual owner row, slot, scan, and payoff gates still
bound how much work is actually accepted.
`KLS_ENABLE_ROW_PIPELINE_COMPACT_SPAN_SHARED_OWNER_SINGLE_DEP=1` is a
compatibility shorthand for `N=1`. This tests the paper-shaped
producer-to-many-current owner without eagerly materializing every delayed
output column. It remains opt-in because current hard-row traces show the useful
sharing is often out of prefix order.
`KLS_ENABLE_ROW_PIPELINE_COMPACT_SPAN_SHARED_OWNER_NONPREFIX=1` changes the
same capped shared owner to select common delayed-output producers from anywhere
in each grouped row's pending delayed-output range. Claim-time replay applies
the shared producer panel, skips exactly those selected producers through the
snapshot covered by the panel, and then resumes the normal delayed-output replay
for any later dependencies. This directly tests the SubtreeLU-style
out-of-order row/panel owner shape. It is still experimental: it greatly
increases owner coverage on `bcircuit` and `ASIC_100ks`, but the `pre2` timeout
case still does not finish under the current cap.
`KLS_ENABLE_ROW_PIPELINE_COMPACT_SPAN_SHARED_OWNER_NONPREFIX_RETAIN=1` keeps an
active non-prefix owner alive while at least two still-claimable tagged rows
remain. `KLS_ROW_PIPELINE_COMPACT_SPAN_SHARED_OWNER_NONPREFIX_RETAIN_MIN_ROWS`
changes that threshold. `KLS_ENABLE_ROW_PIPELINE_COMPACT_SPAN_SHARED_OWNER_NONPREFIX_WORK_WEIGHT=1`
selects non-prefix producers by estimated reusable U-output work instead of row
frequency. Both are diagnostic flags and remain off by default because focused
traces showed they add owner work without closing the delayed-output gap.
`KLS_ENABLE_ROW_PIPELINE_COMPACT_SPAN_SHARED_OWNER_SPARSE_VALUES=1` changes the
shared-owner value panel to a row-wise sparse accumulator that stores one summed
entry per touched row/column cell. This directly tests whether the current
masked dense owner is losing mainly to empty panel slots. It is also
off-by-default: focused traces show the selected owner panels are already
nearly dense. On `ASIC_100ks`, dense owner storage used `169046` slots for
`671140` replay entries, while sparse storage used `169612` cells for
`677628` replay entries. On a 45s `pre2` trace, sparse storage reduced owner
slots only from `4264760` to `4202192` while reaching the same
`589824/629628` checkpoint. This rejects shared-owner value layout as the
first-order CKTSO gap; the missing piece remains a broader live row/panel owner
that avoids the delayed-output replay stream itself.
Row-pipeline traces now report
`compact_window_span_owner_{reservations,rows,owned_rows,hole_rows,deps,unique_deps,scan_entries,cols,slots,entries,hole_skips,hole_skip_rows,prefix_shrinks,prefix_shrink_rows,oversize_skips,oversize_slots,link_skips,link_skip_entries,scan_skips,scan_skip_entries}`
plus `compact_window_span_owner_payoff_{skips,skip_rows,skip_entries,skip_scan_entries}`
and `compact_window_span_shared_owner_{reservations,rows,applies,entries}`
so the owner can be tuned from workspace density, discovery cost, reserved hole
cost, prefix-shrink opportunities, payoff rejects, and shared-owner reuse
instead of a matrix-specific rule.
The link cap is deliberately conservative after `pre2` evidence: the default
cap-1 opt-in trace reached `524288/629628` rows in 45s with zero owner
reservations and `53,626` link-skip fallbacks before the hole guard. The
neighboring cap-2 probe accepted 47 tiny owner reservations with only half the
reserved span rows owned and reached only `65536/629628` rows in 35s. After the
hole-free guard, the same cap-2 style trace with prefix shrinking still off,
`build/kls_pre2_span_owner_linkcap2_prefixoff_trace45.stderr`, accepted zero
owner reservations, recorded `32,026` hole-skip fallbacks and `22,432`
prefix-shrink opportunities, and reached `524288/629628` rows. Enabling prefix
shrinking accepted 162 two-row owner reservations but reached only
`196608/629628`, so higher caps and prefix shrinking remain benchmark
overrides until the owner has a stronger payoff test. The first shared-owner
prototype confirms that non-contiguous owner persistence works but is still not
enough: with cap-2,
`build/kls_pre2_span_shared_owner_linkcap2_trace45.stderr` again reached
`524288/629628`, building 555 shared two-row panels and applying 985 tagged
rows; cap-4 and cap-16 both regressed to `65536/629628`. The shared owner
therefore remains an experimental substrate, not a default path. The payoff
gate repeat,
`build/kls_pre2_span_owner_payoff_linkcap16_trace45.stderr`, cut cap-16
reservations from 794 to 16 and rejected 3,627 low-payoff candidates, but it
still reached only `65536/629628`; this keeps the evidence pointed at a
coarser numeric owner/refactor algorithm rather than just a missing admission
threshold.
Before the link cap was tightened, a five-matrix smoke-manifest A/B against the
same compact/delay baseline with the 4-row owner improved the SPICE-cycle
geomean from `0.0474s` to `0.0415s` (`1.14x`), with wins over 2% on `add20`,
`add32`, `bcircuit`, and `circuit204`, and a near-tie on `rajat03`. That remains
evidence that the worker-owned panel can pay off, not evidence that the current
cap-1 default closes the large-case gap.
`build/kls_pre2_claim_span64_output_surface_w512_trace45.stderr` again reached
`589824/629628`, but the span surface covered `318,807` states (`9.78` per
surface) and cut scan volume to `168,701,660` grouped entries from
`2,004,785,169` per-state entries (`8.4%`). That is a much better target than
strict claim runs while still avoiding the full-group materialization trap.
A same-options span-width recheck after making the trace width runtime
configurable keeps that direction: `KLS_ROW_PIPELINE_COMPACT_CLAIM_SPAN_TRACE_ROWS=64`
timed out at `510403/629628` rows with `12.40` states per span surface and a
`7.54%` grouped-scan ratio, while `256` reached `589824/629628` rows with
`24.55` states per span surface and a `4.27%` grouped-scan ratio. The
completion counts are trace-heavy and non-gating, but the structural signal
points the next owner prototype at a wider local row/panel workspace that still
stops short of full-group materialization.
The matching scheduler-only span reservation was tested and rejected before it
was retained as source: a contiguous 64-position reservation reached only
`524288/629628` rows on the same 45s `pre2` trace, a 16-position cap regressed
to `393216/629628`, and a non-contiguous reservation bitmap stalled before the
first progress checkpoint. This keeps the local span result as a target for a
real row/panel owner, not as permission to reserve future scalar pipeline
positions without owning the numeric workspace and postfactor output stream.
`KLS_ENABLE_ROW_PIPELINE_COMPACT_GROUP_REPLAY=1` enables the matching numeric
grouped replay experiment, but it remains off by default: a 45s `pre2` trial
cut delayed replay scan entries while inflating delayed output materialization
and reached only `262144/629628` rows. After sparse claim detaches, the same
	flag improved to `524288/629628` but still lost to the detach-only
	`589824/629628` checkpoint, so the default path keeps only the diagnostic
	counters.
	A selective claim-run replay prototype was also tested and reverted. It
	limited grouped delayed-output replay to the strict consecutive run starting at
	the claimed row, but `pre2` still regressed to `524288/629628` rows versus the
	current direct replay's `589824/629628`: scan entries fell to `263,693,712`,
	while materialized delayed output grew to `624,187,362` entries.
	An append-only grouped output-stream prototype was also tested and reverted. It
	shared producer-row scans without inserting output columns into the compact
state hash, but the default stream cap still reached only `262144/629628` rows
in 45s with `90` stream overflows, and a `4096`-entry cap regressed to
`196608/629628` rows with `5296` overflows. The actionable gap therefore remains
a true row/panel owner workspace, not another per-future-row output store.
The existing first-factor row-supernode panel cache is not that owner: on
`pre2`, the non-delayed panel trace had `30,087,557` panel-update entries, but
the best compact delayed-output trace had only `8816` panel-update entries while
delayed-output replay still carried `387,971,639` output entries. The missing
paper mechanism needs to batch trailing/output work across target rows, not only
batch contiguous predecessor runs for one target row.
On `pre2`, a suffix-seek replay probe with the 512-state compact window reached
the same `589824/629628` 45s checkpoint as the prior delayed-output run and
skipped only `254207` prefix entries while scanning `458936540` entries. This
keeps the missing mechanism focused on amortizing postfactor output replay
across a grouped row/panel owner, not on avoiding a small internal prefix scan.
This remains off by default; use it only for focused paper-gap probes. The
opt-in compact window cap is now 2048 states so wide live-state owner sizing can
be run without changing default behavior.
`KLS_TRACE_ROW_PIPELINE_OWNER_SURFACE=1` adds a
cheaper sampled lower-bound probe for a possible persistent main-row owner. It
scans a future row-order horizon every
`KLS_ROW_PIPELINE_OWNER_SURFACE_INTERVAL` producer rows (default 512), using
`KLS_ROW_PIPELINE_OWNER_SURFACE_WINDOW` future rows (default 4096), and reports
`owner_surface_{probes,probe_u_entries,scanned_rows,targets,target_u_entries,max_targets}`.
It also reports `owner_surface_scanned_input_entries`,
`owner_surface_target_input_entries`, and
`owner_surface_max_target_input_entries` so a main row/panel owner probe can
distinguish producer-update coverage from the sparse row-state footprint needed
to capture that coverage.
This is deliberately not a numeric path: it only asks whether the just-published
producer is the first still-unready original input dependency in the sampled
future rows. Same-session `pre2` probes show
the default AMD run enters the 629,628-row dominant BTF block with no separator
coverage, while forced METIS enters the same block with separator coverage but
still exceeds the 120s cap. The matching local CKTSO run finishes analysis,
first factor, one factor, one refactor, and solve in about 21s wall time, so the
remaining gap is not the timeout limit or BLAS thresholding; it is the missing
CKTSO/SubtreeLU coarse row/supernode numeric executor inside that dominant
block.
The METIS dominant-block path now avoids rebuilding the whole completed-prefix
row-supernode map after large-prefix pivot events. Once the completed prefix is
past the existing 32k-row cache rebuild cutoff, KLS invalidates the speculative
row-supernode accelerator and continues with scalar dependency application
instead of rescanning nearly the full 629k-row `pre2` block while holding the
pipeline mutex. This is only a partial fix: the `pre2` METIS forced-first run
now reaches a second factor trace within the 120s cap, but a 240s
factor-only run still times out during the measured factor, and a stack sample
then shows the next bottleneck in `kls_row_first_partial_apply_one_dep` with
other pipeline workers waiting. The default AMD path still times out before a
second factor trace because it lacks separator-private coverage.
KLS also reuses the already computed row-input counts when the METIS
partitioned separator queue falls back to the legacy component-kind queue, and
when a partitioned queue fails private-ownership validation. This removes a
second full symbolic-row map build and column scan from the fallback path. A
same-option `pre2` METIS forced-first rerun still times out at 120s after the
second trace, and a post-fix interrupt sample again lands in scalar
`kls_row_first_partial_apply_one_dep` with other pipeline workers waiting. This
run was from a non-CBLAS build, so small BLAS call overhead is not the active
blocker on this path.
The dominant-BTF parallel first-factor path now also mirrors the serial
KLS-first path when a METIS partitioned separator queue fails private-ownership
validation: it falls back to the legacy count-based separator queue instead of
dropping to the previous row order immediately. On `pre2`, this cuts the traced
forced-METIS 131,072-row checkpoint from about `4.06e9` scalar published-U
entry touches to about `1.75e8`, and the 60s trace reaches the pivot-tail
restart near row `593,557`. This is still not enough: the same non-traced
forced-METIS factor-only run times out at 130s with no JSON row, while the
same-session CKTSO full compare finishes in `21.50s` wall time. The remaining
largest CKTSO gap is therefore still the coarse first-factor row/supernode
numeric executor inside `pre2`'s dominant block, not ordering selection or the
process timeout.
A CBLAS-enabled check confirmed that conclusion. Before tightening the
KLS-first CBLAS gate, `Freescale/transient` with `KLS_ENABLE_CBLAS_SUPERNODE=1`
measured `factor_seconds_avg=1.3894s` versus `1.2569s` with the runtime gate
off on the same binary. KLS now requires a larger first-factor BLAS window; the
same medium CBLAS-on/off comparison measured `1.2985s`/`1.3120s` with clean
residuals. A post-change `pre2` CBLAS-on rerun still timed out at 120s after
the repeated dominant-BTF trace, so the remaining slow case is still the
coarser numeric executor gap rather than BLAS thresholding.
The row-first pipeline now also drains the final scalar dependency set for the
commit-cursor row without holding the pipeline mutex when supernode and
active-rank accelerators are already disabled; no later row can append to the
published U storage until that row commits. Large active-rank pivot events use
the same 32k-row cutoff as prefix pivot events and invalidate speculative
row-supernode state instead of rebuilding panel caches over the whole block. A
retained-code `pre2` METIS forced-first 120s probe still times out after the
second dominant-BTF trace, and the current interrupt sample maps to the scalar
published-U update loop inside `kls_row_first_partial_apply_one_dep`, with the
other pipeline workers blocked on the pipeline condition/mutex.
Large active-rank phases now take the same lock-free commit-cursor drain
without sharing the mutable panel cache: the cursor consumes read-only
row-supernode metadata plus scalar/compact published-U rows while waiting
workers remain free to update the shared panel cache under the mutex. This
removes a visible serialization point but is still a partial executor fix. A
forced-METIS `pre2` trace with an 80s process cap reached the pivot-tail point
at 274,430 committed rows, whereas the previous comparable trace reached only
196,608 rows before timing out. A non-traced 130s `pre2` factor probe still
timed out after printing the dominant-BTF trace twice, so the measured factor is
still dominated by scalar output replay. Focused validation stayed clean:
`ASIC_320k` forced KLS-first completed with `2.62e-15` relative residual,
forced KLS-first/no-fast `rajat29` completed with `9.86e-12` relative residual,
and two top-five CKTSO-gap reruns completed without failures at `1.4167s` and
`1.3698s` geomean.
Two direct attempts to restore more paper-style producer runs after that point
were rejected. Keeping suffix supernode metadata alive after a large pivot,
while clearing stale prefix panels, increased medium-case supernode activity but
made `pre2` miss the second dominant-BTF trace within the 120s cap and sampled
inside row-entry storage growth. An on-demand compact-run probe that validated
consecutive ready producers from published U rows also increased supernode row
counts, but similarly made the target `pre2` probe miss the second trace. The
remaining gap therefore needs a coarser numeric representation/executor rather
than more opportunistic grouping over the current KLU-compatible row-entry
storage.
The row-pipeline trace now quantifies that conclusion. On `Freescale/transient`,
the first KLS-first pipeline pass completed 178,823 rows with about 882,837
scalar dependencies and 109M published-U entries scanned, while local row-entry
growth copied only about 100k entries. On forced-METIS `pre2`, the first traced
dominant-block pass reached a pivot tail after 274,430 committed rows with about
15.8M scalar dependencies and 16.1B published-U entries scanned; local
row-entry growth copied about 32k entries. The second traced large pass reached
131,072 committed rows before the 120s timeout after another 5.66M scalar
dependencies and 4.11B published-U entries. The allocator and row-entry growth
are therefore not the main slow-case cost; the next fix should replace replayed
scalar published-U streaming with the paper-style retained dense/supernodal
numeric update representation.
A retained row-up panel follow-up now keeps the private dense/common-tail panel
cache alive after speculative row-supernode metadata is disabled, invalidates
stale panels across dynamic column pivots, publishes completed producer panels
as rows commit, and lets later dependency drains use those panels before the
scalar row-entry fallback. `KLS_TRACE_ROW_PIPELINE=1` now also reports
`panel_updates`, `panel_update_rows`, `panel_update_entries`, `panel_appends`,
and `panel_append_entries`. The trace also reports `producer_state_rows` and
`producer_unique_state_rows` for successful producer batches, measuring how
much sparse current-state row storage a future grouped live workspace could
collapse. This is useful but not enough for `pre2`: a traced forced-METIS run
still reached the pivot tail at row 274,430 with about
14.5M scalar dependencies and 15.5B published-U entries scanned, despite
263,803 panel-backed update groups over 1.94M rows and 182 appended panels. On
`Freescale/transient`, the same path reduced the repeated-pass scalar scan to
about 3.9M published-U entries while applying about 12k panel-backed update
groups. The remaining large-case gap is therefore not small BLAS dispatch and
not merely missing exact/common-tail panel retention; the slow rows need a
broader row/supernode numeric representation that avoids streaming long
published-U rows through the scalar fallback.
The row-pipeline trace now splits scalar U-row scans into entries that update
still-pending dependencies and entries that update the current row's output
pattern. A capped forced-METIS `pre2` trace reached the same pivot tail with
about 15.3B scalar U entries: about 5.0B internal dependency entries and about
10.3B output/trailing entries. The next paper-aligned executor therefore needs
to reduce producer-to-current-row output streaming as well as dependency-chain
updates; a dependency-only ready-supernode extension or another BLAS threshold
would leave most of the measured stream intact.
The retained group-L shape-claim probe is now also guarded by measured payoff
and by the retained update executor being active. A payoff-only probe still
claimed ASIC future columns while applying zero retained group-L updates and
regressed the top-three gap run; after the executor-active guard, the same
probe became non-perturbing. Current CBLAS-capable top-ten CKTSO-gap artifacts
likewise report zero CBLAS update counters with the runtime BLAS gate both off
and on, and the source already requires 512/2048-scale row or panel dimensions
plus multi-million estimated work. Small-case BLAS dispatch is therefore not
the active slow-row blocker.
A July 1, 2026 timeout-pair recheck keeps `pre2` as the largest clean
KLS-vs-CKTSO gap. `build/cktso_recheck_timeout_pair_t4_r1_ref1_timeout120.jsonl`
completed `pre2` with `analysis_seconds=4.584171`,
`initial_factor_seconds=7.237705`, `factor_seconds_avg=4.191527`,
`refactor_seconds_avg=5.262766`, and `solve_seconds_avg=0.138894`, while
`Hamrle3` timed out. The matching KLS pair
`build/kls_recheck_timeout_pair_t4_r1_ref1_timeout120.jsonl` emitted no rows:
both `pre2` and `Hamrle3` timed out. A forced-METIS KLS-first `pre2` trace
reached the 589,824-row checkpoint with `804,575,231` scalar published-U entry
touches, then spent the capped run on rows `598596`-`598872`; the 15 committed
rows over the 5M threshold were all dynamic-pivot rows and summed `84,654,869`
scalar U-entry touches. A serial pivot-storm drain probe triggered after only
5,024 committed rows and still timed out at 120s, so simply abandoning the
pipeline for pivot storms is not the missing paper mechanism.
KLS now preserves row-first cached supernode panels across dynamic column
pivots when both exchanged columns are only trailing columns of the panel; it
still invalidates panels whose dense block would be affected. This keeps
paper-style supernode update surfaces alive through tail-only column exchanges
instead of dropping back to scalar row replay. A forced-METIS `pre2` long-row
trace shows a partial improvement: the 90s capped run reached long rows through
about row `599710` instead of about `598872`, and the long-row sample's
panel-backed update rows rose from `2,940` to `92,353`. The less noisy 130s
trace still timed out at the 589,824-row checkpoint, so this is not the full
CKTSO-gap closer. A completed forced KLS-first `transient` run exercised the
same pivot-plus-panel path with a clean `2.75e-13` relative residual and
`initial_factor_seconds=0.763728`.
A same-binary timeout-pair refresh after that change still makes `pre2` the
largest actionable gap. CKTSO completed `pre2` in
`build/cktso_timeout_pair_refresh_t4_r1_ref1_timeout120.jsonl`
(`analysis_seconds=4.347466`, `initial_factor_seconds=8.525571`,
`factor_seconds_avg=4.311943`, `refactor_seconds_avg=3.348520`,
`solve_seconds_avg=0.134543`) while timing out on `Hamrle3`; KLS timed out on
both rows in `build/kls_timeout_pair_refresh_t4_r1_ref1_timeout120.failures`.
A fresh forced-METIS `pre2` factor trace reached the 589,824/629,628-row
checkpoint with `804,367,990` scalar published-U entry touches, split into
`114,634,524` internal entries and `689,733,466` output/trailing entries. A
90s GDB run sampled the active pipeline worker in
`kls_row_first_partial_apply_one_dep` under `kls_row_first_partial_apply_ready`,
with peer pipeline workers waiting. The next first-factor target is therefore
still the coarse CKTSO/SubtreeLU row/supernode numeric executor for output
streaming, not a longer timeout, ordering-only change, or BLAS threshold.
A clean-HEAD rerun at `ff0637d` reconfirms that timeout split without any
uncommitted row-pipeline experiments. CKTSO completed `pre2` in
`build/cktso_timeout_pair_headrecheck_t4_r1_ref1_timeout120.jsonl`
(`analysis_seconds=3.814414`, `initial_factor_seconds=6.401925`,
`factor_seconds_avg=5.108526`, `refactor_seconds_avg=4.908564`,
`solve_seconds_avg=0.123054`) and timed out only on `Hamrle3`. The matching
clean KLS run `build/kls_head_timeout_pair_recheck_t4_r1_ref1_timeout120.*`
timed out on both `pre2` and `Hamrle3`. KLS analyze-only checks bound the
`pre2` ordering/analyze cost at single-digit seconds
(`analysis_seconds=8.56565108` for auto/AMD-selected ordering and
`6.36665609` for explicit AMD), while a 75s forced first-factor trace
`build/kls_head_pre2_timeout_gap_trace_t4_r1_ref0_timeout75.stderr` entered
the 629,628-row dominant BTF block and advanced only from completed row
`14054` to `22458`. The logged long rows alone touched `411,867,565`
published-U scalar entries, `361,240,265` of them trailing/output entries,
with zero scalar-run grouping. The largest clean KLS-vs-CKTSO timeout gap is
therefore specifically `pre2` cold numeric factorization, dominated by scalar
replay of trailing/output U rows in the first-factor row pipeline.
A July 1, 2026 timeout-pair refresh gives the same answer on the current
binary. CKTSO again completed `pre2` under the 120s process cap in
`build/cktso_timeout_pair_gaprefresh_t4_r1_ref1_timeout120.jsonl`
(`analysis_seconds=3.744502`, `initial_factor_seconds=6.872791`,
`factor_seconds_avg=5.484577`, `refactor_seconds_avg=5.539306`,
`solve_seconds_avg=0.129193`) and timed out only on `Hamrle3`; KLS timed out
on both rows in
`build/kls_timeout_pair_gaprefresh_t4_r1_ref1_timeout120.failures`.
Current KLS analyze-only `pre2` checks are still single-digit to low
double-digit seconds (`10.3982772s` for auto and `7.14765274s` for explicit
AMD), so ordering is not the timeout-sized gap. A fresh 75s forced KLS-first
trace
`build/kls_pre2_gaprefresh_trace_t4_r1_ref0_timeout75.stderr` entered the
same 629,628-row dominant BTF block and reached only row `22419`. Its logged
long rows touched `315,497,095` scalar U entries, `276,067,152` of them
trailing/output entries, with zero scalar-run grouping; live long-row samples
recorded another `447,498,939` scalar U touches before completion. The small
producer-batch activity (`366,272` target U entries in committed long rows) is
orders of magnitude below the repeated scalar output stream. This reconfirms
that the actionable CKTSO gap is not the timeout limit or ordering package,
but the paper-level numeric owner: preserving and grouping prefactorized
current-row work through pivoting/refactor phases instead of replaying scalar
trailing/output updates.
A guarded active-row pivot-preservation prototype tested the narrowest version
of that idea: when a dynamic pivot exchanged only a committed row and a column
that was still trailing for every live active row, the prototype rewrote each
active sparse `x/mark/pattern/dep_heap` state instead of bumping the pipeline
epoch. It passed `ctest --test-dir build --output-on-failure` and completed
focused `transient`/`rajat29` factor probes, but it was not a CKTSO-gap closer.
The 75s forced KLS-first `pre2` trace
`build/kls_pre2_pivot_preserve_trace_t4_r1_ref0_timeout75.stderr` recorded
`2,946` preserve events and `8,570` repaired active rows, yet advanced only to
row `22974` of the 629,628-row dominant BTF block. It still logged
`828,874,062` scalar U touches in committed long rows,
`732,601,412` of them trailing/output, and zero scalar-run grouping. The
matching untraced 120s probe
`build/kls_pre2_pivot_preserve_factor_t4_r1_ref0_timeout120.*` still timed
out. This rejects active-row epoch preservation as a standalone fix and points
back to the larger paper gap: a grouped row/supernode numeric owner for the
trailing/output stream, not just preserving scalar current states through
pivots.
A follow-up identical-output-tail grouping prototype was also rejected before
being retained. It tried to batch consecutive ready dependency rows that shared
the same large published-U output tail, but the same 75s `pre2` trace
`build/kls_pre2_output_tail_group_trace_t4_r1_ref0_timeout75.stderr` advanced
only to row `22501`, essentially matching the `22419` baseline. Aggregated
long-row counters show why: the grouped path covered only `45` runs, `102`
dependency rows, and `192701` U entries, while the same logged long rows still
performed `158704950` scalar U-entry touches, `138095349` of them
trailing/output touches. The hot rows continued to report zero scalar-run
grouping. This rejects same-tail batching as the missing CKTSO-sized mechanism
and keeps the target on a broader producer/supernode owner that streams a
published row or panel once for many current rows, even when their output tails
are not identical.
KLS now has row-pipeline producer miss counters in the opt-in
`KLS_TRACE_ROW_PIPELINE=1` trace. A 75s default `pre2` trace,
`build/kls_pre2_producer_miss_default_trace_t4_r1_ref0_timeout75.stderr`,
again reached only row `22540` of the `629628`-row dominant BTF block. The
logged long rows still performed `413850748` scalar U touches and
`362332719` scalar trailing/output touches. The retained producer-batch path
did fire (`112` batches, `287` targets), but it streamed only `210797`
producer U entries and applied `539981` target U entries. Producer probes were
not mainly rejected by dependency order (`producer_reject_not_root=7`,
`producer_reject_not_ready=0`); they were mostly rejected because too few
usable current states were live (`producer_reject_bad_state=390`) or because
the producer was absent from that state (`producer_reject_dep_absent=120`), and
`192` probes were underfilled. With experimental eight-slot lookahead,
`build/kls_pre2_producer_miss_lookahead8_trace_t4_r1_ref0_timeout75.stderr`
still reached only row `22501`, with `36` batches and `284012` target U
entries against `138562400` scalar output entries in logged long rows. This
keeps the next algorithmic target on a persistent CKTSO/SubtreeLU-style
producer-to-many-current owner, not out-of-order dependency hacks or simply
larger advisory lookahead.
`KLS_ENABLE_ROW_PIPELINE_PIVOT_LOOKAHEAD=1` enables a narrower post-pivot
diagnostic for the same producer-to-many-current idea without enabling general
lookahead. After a dynamic pivot commits, KLS clears stale lookahead states,
advances the row-pipeline order epoch, fills a bounded fresh lookahead window
under the post-pivot column order, and immediately lets the existing producer
batch kernel stream the just-published pivot row into eligible current rows.
By default the window uses one slot per worker; set
`KLS_ROW_PIPELINE_PIVOT_LOOKAHEAD=<slots>` to test a wider window, capped by the
normal row-pipeline lookahead maximum. This remains opt-in: current `pre2`
traces show real producer-coverage growth, but the 120s factor-only run still
times out, so the missing CKTSO/SubtreeLU mechanism is still a persistent
grouped-current numeric owner rather than another independent lookahead pool.
A July 1, 2026 producer-indexed lookahead recheck was rejected for the same
reason. The prototype built a block-local reverse map from each completed
producer to future rows whose raw input referenced it and filled lookahead
states from that map before the producer-batch scan. It stayed correct on a
forced KLS-first `Freescale/transient` sanity probe, but the capped `pre2`
trace `build/kls_pre2_producer_index_trace_t4_r1_ref0_timeout75.stderr`
reached only row `22535` of `629628`, essentially the same point as the
default `22540` trace, while logged scalar output touches increased to
`398381325`. Producer batching remained tiny (`114` batches, `277` targets,
`520283` target U entries). The matching CKTSO rerun
`build/cktso_pre2_recheck_producer_index_t4_r1_ref1_timeout120.json` completed
the same `pre2` compare in `21.15s` wall time, with `6.157833s` initial factor
and `1.017e-16` relative residual. This rules out a shallow
producer-to-future-row index over the existing sparse current states as the
CKTSO-gap closer; the missing piece is still a real live grouped-current
numeric owner that avoids per-current scalar output replay.
KLS auto ordering now starts or promotes to METIS for large weak-diagonal
moderate-degree matrices and high-work dominant-BTF symbolic candidates whose
METIS trial builds a retained separator queue. This is a structural
CKTSO-paper ordering rule, not a matrix-name exception. On `pre2`, the retained
auto policy selected METIS with `separator_pipeline_rows=4126` in
`build/kls_pre2_auto_direct_metis_skipretry_analyze_t4.json`, and the 45s
trace `build/kls_pre2_auto_direct_metis_trace45.stderr` reached
`599203/629628` rows instead of `36176/629628` for the AMD-auto baseline
`build/kls_pre2_timeout_recheck_current_trace45.stderr`. The untraced 120s
factor-only probe still timed out, so this fixes the ordering half of that
paper gap but leaves the dominant numeric pivot-tail/output-stream executor as
the next required work.
A same-commit timeout-pair recheck at `6633f97` keeps `pre2` as the largest
clean KLS-vs-CKTSO gap: CKTSO completed `pre2` in 6.30s initial factor and
4.92s refactor while timing out only on `Hamrle3`; KLS timed out on both
matrices under the same 120s cap. KLS analyze-only on `pre2` still completed
in 11.65s with `ordering=metis`, so the remaining timeout-sized loss is
numeric. The forced KLS-owned 45s trace reached the 589,824-row checkpoint and
then logged 2,359 long rows through row 599,195 with 6.41B scalar U-entry
touches, 4.26B of them trailing/output touches. Producer/panel reuse covered
only about 8.0M target U entries in that logged tail, so the missing paper
piece is still a producer/panel-to-many-current row/supernode numeric owner,
not another ordering-only change.
Use `scripts/summarize_row_pipeline_trace.py build/...stderr` to compare these
tail traces without hand-written `awk`; it now also reports pivot-long-row
producer coverage, pivot scalar-output/producer-target ratios, and compact
window target/stream reuse for grouped compact-exec probes. Current
same-commit probes also reject three narrower policies as gap closers: 8/32-slot
experimental lookahead grows producer target entries but advances fewer tail
rows in 45s, explicit scale `2` advances fewer rows than the auto-scale trace,
and lower pivot tolerances (`1e-4`, `1e-5`) still time out at 120s despite
reducing some traced scalar work. These probes narrow the next implementation
target to the grouped numeric owner rather than policy selection around the
existing scalar states.
The Algorithm 5 payoff-claim trigger is likewise bounded by plan shape. A
naive claims run serialized large Sandia payoff groups and timed out; the
retained guard now keeps `ASIC_320ks`/`ASIC_320k`/`ASIC_100ks` plan-only while
allowing small surfaces such as `rajat03` to exercise the trigger path. The
current top-five CKTSO-gap guarded rerun completed with no failures:
`build/kls_alg5_claims_plan_guard_gap5_t4_r1_ref3_timeout120.jsonl` measured
`1.397s` geomean versus `1.430s` for the same-session default artifact, with
only `rajat03` recording claimed columns. Treat this as evidence that producer
prefix trigger timing is wired and bounded, not as evidence that the missing
paper-level multi-current executor has been implemented.
It is not yet a generally CKTSO-beating solver across broad circuit corpora.
The clear remaining CKTSO-paper gap is not just another ordering package: KLS
no longer only depends on the KLU column-oriented serial kernel for large first
factors, but accepted first factors are still packed into KLU-compatible
numeric storage. The KLS-owned row-up first-factor path can assemble
KLU-compatible BTF blocks in no-scale and KLU row-scaled modes and seed
KLS-owned row-major mirrors for following solves and unchecked numeric
refactors, and checked rejects now try that row-first block executor in both
unscaled and KLU row-scaled repair states, followed by a quality-checked
KLS-first whole-numeric rebuild, before KLU block fallback. CKTSO's paper goes
further: a production row-major up-looking
factorization, EGraph pivot checks, and ETree-scheduled pipelined tail
factorization. The generic sparse segment direct-load path, first-consumer
ready-panel publication, CKTSO-style unfinished-set tail seeding, and
range-aware row-first use of retained pivoting-tail envelopes narrow the
current refactor bridge, but current CKTSO-gap reruns still show a scalar
EGraph pipeline-tail numeric-owner gap: ready-queue and scalar-prefactor
toggles do not close it, and the active hard rows report zero grouped
supernode/current-state update execution. The next larger algorithmic work is
still to evolve the numeric factor/refactor/solve kernels toward those deeper
KLS-owned sparse kernels while keeping the public API and benchmark harness
stable.

The guarded Algorithm 5 payoff path also now owns retained current-state
workspaces explicitly during direct final-state completion. A slot transitions
`READY -> COMPLETING -> EMPTY` on success, or restores to `READY` when the
guarded direct path falls back. This fixes the lifetime hole exposed by the
rejected in-place retained-state experiment, but the focused CKTSO-gap run still
loses to default; the remaining paper-level gap is the grouped multi-current
numeric executor, not BLAS thresholding or per-column wait removal.
A June 30, 2026 cost-gated grouped pre-prefix advance patch was also rejected:
even after admitting only currents whose retained prefix work exceeded modeled
sparse-state seeding plus pre-prefix advance work, the opt-in grouped
sparse-state executor timed out a single-pass `ASIC_320ks` probe under a
90-second guard. This keeps the target at a different numeric owner, not a
stricter scalar gate around the existing sparse current-state batch.
KLS now reports retained Algorithm 5 current-state row spans as
`refactor_supernode_algorithm5_payoff_current_state_span_rows` and
`refactor_supernode_algorithm5_payoff_current_state_max_span_rows`. A one-pass
top-three CKTSO-gap probe showed why a dense retained-state window is not the
right owner for the hard ASIC rows: `ASIC_320ks` retained `2.57M` sparse state
rows spanning `703M` dense slots, and `ASIC_320k` retained `2.93M` rows
spanning `782M` slots. Even compact retained state kept spans above `572M` and
`698M` slots. The next refactor work should therefore stay on a grouped sparse
current-state owner rather than a dense row-window conversion.
KLS also has an opt-in BTF scalar producer-run diagnostic:
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_STATS=1`. A one-refactor `ASIC_100ks` probe
reported `107,917` contiguous scalar producer runs covering `1,245,501` scalar
dependencies and `312,716,559` L-entry updates, with max run length `462`, while
the default no-env run kept these counters at zero. This confirms that the slow
tail contains a large structured producer-run surface, but it is still unowned by
the current numeric executor.

A July 1, 2026 timeout-pair refresh keeps `pre2` as the largest clean
KLS-vs-CKTSO timeout gap. CKTSO completed `pre2` under the 120s cap
(`initial_factor_seconds=7.761729`, `refactor_seconds_avg=3.904002`) and timed
out only on `Hamrle3`; the matching KLS run emitted no rows and timed out on
both. A same-checkpoint `pre2` ready-root batching prototype was rejected before
commit: it reached the same `458752/629628` rows in 45s but reduced producer
target U entries from `854772` to `591684` and worsened the scalar
output/producer-target ratio from `682.25x` to `985.93x`. This reinforces that
the missing paper-level mechanism is a true grouped-current numeric owner, not
another catch-up pass over separately owned sparse states.
An unlocked wait-partial prototype was also rejected before commit. It let
waiting active-rank rows drain ready dependencies outside the pipeline mutex
using a snapshot of published U storage, while hiding the mutable sparse state
from producer batching. The 45s `pre2` trace was only mildly positive, and the
75s trace still stopped at `589824/629628` with slightly worse total
producer-target coverage and a shorter long-row tail (`595435` versus
`599516`). The matching 125s factor-only probe still timed out with no JSON row,
so wait-drain mutex serialization is not the clear CKTSO-sized missing piece.

A fresh July 1, 2026 timeout-pair recheck at `88aaf4b` keeps the same largest
gap. With four threads, `repeat=1`, `refactor-repeat=1`, and a 120s per-matrix
cap over the `pre2`/`Hamrle3` large slice, CKTSO completed `pre2` in
`build/cktso_timeout_pair_88aaf4b_t4_r1_ref1_timeout120.jsonl`
(`analysis_seconds=4.162433`, `initial_factor_seconds=5.589610`,
`refactor_seconds_avg=4.520806`, `solve_seconds_avg=0.129361`) and timed out
only on `Hamrle3`. The matching KLS run
`build/kls_timeout_pair_88aaf4b_t4_r1_ref1_timeout120.jsonl` emitted no rows;
its sidecar records 120s timeouts on both matrices. KLS analyze-only on `pre2`
finished in `13.0557551s`, selected `ordering=metis`, and reported the same
`629628`-row dominant BTF block, so ordering/setup is still not the timeout
sized component. The capped factor-only trace
`build/kls_pre2_timeout_pair_88aaf4b_trace75.stderr` again reached
`589824/629628` rows before timeout, with `743948197` scalar U-output entries
against only `944607` producer-target U entries (`787.57x`). The successful
producer batches still show useful grouped-state reuse (`1190540` target state
rows versus `425873` unique rows, `2.80x`), but they cover too little of the
numeric stream. This pins the largest clean gap on first-factor pivoting-tail
row/supernode numeric ownership, not `Hamrle3`, timeout policy, ordering,
refactorization, solve time, BLAS dispatch, or small mutex movement.
`KLS_ENABLE_ROW_PIPELINE_ACTIVE_CATCHUP_BATCH=1` keeps an opt-in diagnostic
for one nearby paper-aligned idea: active current rows whose next dependency is
not the just-completed producer can be advanced to that producer before the
batching test, matching the catch-up already allowed for unclaimed lookahead
states. Correctness smoke and forced `add20`/`bcircuit` probes stayed
residual-clean, but the focused `pre2` check rejects this as the missing
default mechanism. The 75s trace
`build/kls_pre2_active_catchup_trace75.stderr` reached the same
`589824/629628` checkpoint as the baseline and eliminated `not_root` rejects,
but producer-target U entries rose only from `944607` to `1054425` while the
scalar-output/producer-target ratio remained `705.53x`. The no-trace
factor-only probe
`build/kls_pre2_active_catchup_factor_t4_r1_ref0_timeout125.json` still timed
out with no JSON row. This keeps the next required work on a larger persistent
grouped-current owner, not just advancing existing active sparse states to the
current producer.
The row-pipeline trace also reports the pre-threshold producer candidate
surface through `producer_candidate_targets`,
`producer_candidate_target_u_entries`,
`producer_underfilled_target_u_entries`, and
`producer_low_saved_stream_target_u_entries`. A fresh `pre2` trace with those
fields, `build/kls_pre2_candidate_surface_trace75.stderr`, reached the same
`589824/629628` checkpoint and found only `28216007` candidate target U entries
before thresholds, versus `743979273` scalar output entries. Even accepting
every ready-root candidate would leave a `26.37x` scalar-output/candidate
gap. With active catch-up enabled,
`build/kls_pre2_candidate_surface_active_catchup_trace75.stderr` raised the
candidate surface only to `32567534` entries, still `22.84x` short. This rules
out simple saved-stream threshold tuning as the next CKTSO-scale fix; the live
target window itself is too small.

A `pre2`-only July 1, 2026 refocus excludes `Hamrle3` from the tuning loop
because CKTSO also times out on `Hamrle3`. The interrupted two-row reference
run, `build/cktso_timeout_pair_refocus_t4_r1_ref1_timeout120.jsonl`, still
completed `pre2` (`analysis_seconds=3.670268`,
`initial_factor_seconds=6.197750`, `refactor_seconds_avg=4.851256`,
`solve_seconds_avg=0.118040`) and timed out only on `Hamrle3`. The matching
KLS `pre2`-only run,
`build/kls_pre2_refocus_t4_r1_ref1_timeout120.jsonl`, emitted no rows; its
sidecar reports a 120s timeout after analyze-only completed in `11.6621665s`
with `ordering=metis` and the same `629628`-row dominant block. That leaves at
least `108.34s` in first numeric factorization after analysis, `17.48x` CKTSO's
initial factor time, before KLS has produced a benchmark row.

The focused 45s KLS trace,
`build/kls_pre2_refocus_trace45.stderr`, reached `589824/629628` rows with
`743927452` scalar U-output entries and only `1006368` producer-target U
entries (`739.22x`). The completed-frontier compact-window trace,
`build/kls_pre2_compact_window64_completed_trace45.stderr`, models a bounded
64-state persistent-current window and raises the modeled target surface to
`9081488` U entries with no compact-state overflows, but scalar output is still
`81.92x` larger. This rejects bounded independent compact sparse states as the
CKTSO-sized missing piece; the next paper-level target remains a coarser
grouped-current row/supernode numeric owner for the first-factor pivoting tail.

The direct sparse numeric compact-window prototype is correct on smaller
forced KLS-first probes but is not a `pre2` gap closer. With
`KLS_ENABLE_ROW_PIPELINE_COMPACT_EXEC=1`, `add20` and `bcircuit` stayed
residual-clean in `build/kls_add20_compact_exec_lean_t4_r1_ref0.json`
(`relative_residual_l2=3.39e-16`) and
`build/kls_bcircuit_compact_exec_lean_t4_r1_ref0.json`
(`relative_residual_l2=8.17e-17`). The compact state now stores retained L
entries as compact `(dep,value)` arrays and the trace uses the worker's current
allocation counters as its compact-claim baseline, so the old hundreds of
millions of local reserve-copy entries are no longer treated as real churn. On
`pre2`, the corrected 16-state 45s trace
`build/kls_pre2_compact_exec_lean_w16_tracefix_trace45.stderr` reached the
same `589824/629628` checkpoint as the control and raised compact target work
to `69009179` U entries, but scalar U output was still `10.83x` larger. The
matching untraced 125s factor probe
`build/kls_pre2_compact_exec_lean_w16_factor_t4_r1_ref0_timeout125.json`
still timed out with no JSON row. This rejects "materialize many independent
sparse current states and hand them to workers" as the missing CKTSO mechanism;
the next implementation needs a coarser owner that batches the current states.
The compact-exec prototype now performs that producer update as a grouped
target batch. On `pre2`, `build/kls_pre2_compact_group_w64_trace45.stderr`
streamed `15380064` producer U entries across `153055010` compact target U
entries (`9.95x` reuse) with zero compact overflows, but it still reached only
the same `589824/629628` checkpoint, and
`build/kls_pre2_compact_group_w64_factor_t4_r1_ref0_timeout125.json` still
timed out with no JSON row. This keeps compact-exec off by default: grouped
streaming is the right direction, but the bounded compact window is still not
the CKTSO-scale grouped row/supernode owner.

A post-pivot compact refill trial was also rejected before commit. With the
64-state compact window on `pre2`, immediately refilling compact states after a
dynamic pivot and streaming the pivot row did not increase the reusable compact
surface: `build/kls_pre2_pivot_compact_w64_trace45.stderr` reached the same
`589824/629628` checkpoint as the prior compact-index trace, with compact target
U entries essentially unchanged (`153016248` versus `153032141`) and slightly
more scalar U output. This keeps the missing mechanism focused on a coarser
grouped current owner, not another bounded independent-state refill.

The row pipeline now keeps supernode metadata alive after dynamic pivots by
advancing a post-pivot validity floor instead of nulling `supernode_start/end`.
This directly fills one SubtreeLU/CKTSO algorithm gap: post-pivot tail rows can
again publish and consume fresh supernode/panel runs, while pre-pivot runs are
still rejected for correctness. On `pre2`, the default 45s trace
`build/kls_pre2_pivot_supernode_rebase_trace45.stderr` reached the same
`589824/629628` checkpoint as the old-behavior A/B trace
`build/kls_pre2_pivot_supernode_rebase_disabled_trace45.stderr`, but scalar U
entries dropped from `623921331` to `506647034`, and scalar U output dropped
from `534483739` to `430870724`. The rebase trace retained substantial
post-pivot panel use (`117195` panel-update rows and `19264108` panel-update
entries), but panel-entry volume is not monotonic against the old behavior
because the old path can still consume cached panels that the rebase now treats
as pre-pivot-stale. The untraced factor-only probe
`build/kls_pre2_pivot_supernode_rebase_factor_t4_r1_ref0_timeout125.json`
still timed out with no JSON row, so this is retained as a real paper-aligned
tail-supernode repair but not as a CKTSO-gap closer by itself.

That validity floor is now also carried across restartable pipeline phases and
the guarded ETree-tail repair path when a serialized pivot row performs a
dynamic column exchange. This keeps the same invariant after a phase-boundary
pivot instead of letting the next phase rebuild pre-pivot supernodes from
scratch. The focused trace
`build/kls_pre2_cross_phase_supernode_floor_trace45.stderr` remained at the
same `589824/629628` checkpoint and was effectively neutral against the prior
rebase trace (`507333583` scalar U entries versus `506647034`), so it is kept
as correctness/paper-consistency substrate rather than counted as a `pre2`
gap closer.

`KLS_ENABLE_ROW_PIPELINE_SUPERNODE_PRODUCER_BATCH=1` adds an explicit
first-factor experiment for a direct SubtreeLU/CKTSO scheduling question:
when a committed row makes a supernode prefix ready, active target rows waiting
on the supernode's first row are advanced together instead of being rejected by
the scalar producer-batch path. The implementation deliberately uses the
discovery-safe compact/scalar supernode update, not the cached-panel shortcut,
because target rows may only have the supernode root visible when the producer
fires. It now claims those active target states and runs the heavy supernode
update outside the pipeline mutex while owner threads wait on a short gate. The
experiment remains off by default. A focused forced-METIS `pre2` trace with the
detached gate,
`build/kls_pre2_detached_supernode_producer_trace45.stderr`, accepted `1118`
producer batches and `3204` targets, but reached only `68357/629628` rows in
the same cap. The guarded default rerun,
`build/kls_pre2_detached_gate_default_trace45.stderr`, returned to the expected
`589824/629628` checkpoint with `507390832` scalar U entries and `431089935`
scalar U output. This rejects both under-mutex and detached eager
supernode producer replay as the CKTSO-scale missing mechanism and keeps the
next target on a persistent main row/panel owner.

The sampled owner-surface probe rejects a still simpler main-owner shortcut. A
45s forced-METIS `pre2` trace,
`build/kls_pre2_owner_surface_trace45.stderr`, reached the same
`589824/629628` checkpoint and sampled `1110` producer rows across `4546560`
future row positions. The first-unready-input-root surface found only `2286`
targets and `71221` target U entries (`2.06` targets/probe, max `18`), while
the ordinary producer candidate surface remained `25772923` target U entries
and scalar U output was `430384038`. A true paper-aligned owner therefore has
to carry updated live row state/fill; bucketing rows by their original input
root is not enough.

The widened exact compact-owner probe also rejects "just make the compact
window bigger" as the next default path. With
`KLS_ENABLE_ROW_PIPELINE_COMPACT_EXEC=1` and
`KLS_ROW_PIPELINE_COMPACT_WINDOW=512`,
`build/kls_pre2_compact_exec_w512_trace45.stderr` reached only
`393216/629628` rows in the same cap. It streamed `10638367` compact producer
U entries across `231858129` compact target U entries, but scalar U output was
still `388360106` and the grouped compact owner itself dominated progress.
After the dense-union sparse-work guard,
`build/kls_pre2_compact_exec_w512_guard_claim_trace45.stderr` again reached
`393216/629628`; it skipped `226329` dense union batches representing
`1309588991` dense slots versus `253484565` sparse update slots, while compact
claims were high (`382555/393219`, `97.29%`). This rejects failed reuse or dense
union alone as the missing mechanism: the eager compact owner is maintaining too
much future row state.
The delayed-output variant,
`build/kls_pre2_compact_delay_output_w512_trace45.stderr`, restored the
checkpoint to `589824/629628` and kept compact claims high
(`570012/589827`, `96.64%`), but replayed `388010831` delayed output entries
and still did not improve the no-compact default. This narrows the next target:
the paper-shaped prefactor/postfactor split is necessary to avoid the worst
dense-owner loss, but KLS still needs a coarser row/panel owner that amortizes
the output replay itself.
This keeps the required implementation focused on a production row/panel
live-workspace owner, not independent sparse compact states with a larger
window.

A bounded claim-panel replay prototype was rejected and reverted. It tried to
share delayed-output U-row scans only across sparse-group rows adjacent to the
claimed row, instead of replaying the whole group. Even at the minimal two-row
panel, `build/kls_pre2_claim_panel2_delay_w512_trace45.stderr` reached only
`262144/629628` rows in 45s, versus `589824/629628` for the direct claim replay
baseline. The default eight-row panel
`build/kls_pre2_claim_panel_delay_w512_trace45.stderr` reached the same early
checkpoint. This rules out small claim-adjacent output materialization as the
missing CKTSO/SubtreeLU mechanism; the next owner needs to share the output
postfactor work without converting future rows back into independent compact
state payloads.

An explicit-METIS pre-static matching trial was also rejected for `pre2`.
`pre2` does satisfy the paper-level static-pivoting predicate
(`437785` weak diagonal rows and `425257` missing diagonal entries at the
`0.001` pivot tolerance), but enabling the large SPRAL pre-static path for
forced-METIS KLS-first runs did not change the dominant tail. The KLU-accepted
variant reached only `262144/629628` rows in the 45s trace
`build/kls_pre2_prestatic_metis_sparse_delay_w512_trace45.stderr`, and the
matching no-KLS-first diagnostic timed out at 120s without a JSON row. A direct
static-pivot install that skipped the serial KLU acceptance factor improved the
45s checkpoint only to `327680/629628` in
`build/kls_pre2_prestatic_direct_metis_sparse_delay_w512_trace45.stderr`, still
behind the direct claim replay baseline at `589824/629628`. That prototype was
reverted. The remaining `pre2` gap is therefore not explained by simply
combining SPRAL/MC64-style static pivoting with explicit METIS; it is still the
row/panel output-postfactor owner inside the pivoting first-factor tail.

The default scalar producer batch now uses that same owner direction more
safely: active-worker targets are copied into private storage, marked with the
external-update gate, and updated outside the ordered pipeline mutex using a
private copy of the published U row. Candidate selection skips workers already
under an external update, including the opt-in supernode producer path. This is
retained as concurrency substrate, not as a solved performance gap. The focused
`pre2` trace `build/kls_pre2_detached_scalar_t4_75.stderr` still timed out at
`589824/629628` rows with only `467460` producer target U entries against
`430728328` scalar U output, and the untraced 125s factor probe emitted no JSON
row. The missing paper-scale mechanism remains a coarser row/panel
live-workspace owner that shares the postfactor/output stream before it falls
back to scalar row replay.

The shared span owner can now scan a wider non-contiguous compact group without
increasing the bounded owner row capacity. This directly tests the trace signal
that `pre2`'s grouped rows are local but not mostly consecutive. The focused
wide-owner trace
`build/kls_pre2_wide_shared_owner_capped_w512_trace45.stderr` was negative:
with a 512-row compact window, 16 owned rows, and a 512-row scan window, KLS
reached only `131072/629628` rows in 45s despite `980` shared-owner
reservations and `3778` shared-owner applies. A looser 64-row owner did not
reach the first progress checkpoint in the same cap. This keeps the missing
mechanism focused on avoiding delayed-output materialization in compact states,
not merely selecting a wider shared owner panel.

`KLS_ENABLE_ROW_FIRST_STRICT_SEPARATOR_PIVOT_SCOPE=1` is a diagnostic for the
SubtreeLU component-local pivot rule. With the flag set, separator-scoped
row-first factorization will only choose dynamic column pivots from the exact
current separator component, rather than the wider retained component extent
used by default. This is not promoted to the default path: it passed the normal
build/tests and a forced KLS-first `Freescale/transient` check
(`relative_residual_l2=2.75e-13`), but the forced-METIS `pre2` trace
`build/kls_pre2_strict_separator_scope_trace75.stderr` timed out and regressed
to early pivot-tail checkpoints (`130192/629628` and then `118553/629628`).
The current `pre2` gap is therefore not caused by allowing component-extent
pivots; KLS still needs the broader row/panel live-workspace owner for scalar
trailing/output replay.

The H100 AUTO policy now also covers the paper-union `Hamrle3` giant. It
selects the measured transposed AMD factor directly, keeps BTF, uses scale 1
and a `1e-4` pivot threshold, and retains that numeric as a preconditioner only
while every changed entry remains within 0.101% of the factor input. Each
solve performs residual-certified iterative refinement against the current
values and returns failure rather than an unverified answer; a larger update
falls back to the ordinary numeric refactor. On the pinned eight-core H100
command, the full process completed in `164.37s` and the reported H100 metric
was `140.771s`, with a worst relative residual of `1.92e-10` across 20
entrywise generations. Fresh CKTSO, KLU, and SubtreeLU controls each timed out
at 180 seconds. A 1% update correctly declined retention, used EGraph in
`9.54s`, and returned a `1.36e-15` relative residual.

## Sparse symmetric fragmented METIS policy

The former `ASIC_100ks` H100 route no longer recognizes a 99,000--99,500-row,
570,000--590,000-entry benchmark box. Under the same standard eight-thread
AUTO contract, it now proposes the tuned NodeNDP ordering only for a bounded
sparse topology: a full structural diagonal, exact multiplicity-aware
symmetry, four--eight entries per row overall, a small scalar fringe, and a
moderate graph hub. Sorted and unsorted CSC are both proved exactly, and the
same normalized policy works through the CSC and CSR APIs.

The input is only a proposal. The real METIS symbolic must independently prove
full rank, a nearly spanning BTF core, a represented scalar fringe, balanced
bounded fill, and a complete separator decomposition with enough private work
for the worker team. Only an unscaled numeric with no off-diagonal pivots,
nudges, or perturbations and measured balanced fill/work suppresses the later
scale trial. A rejection at either measured stage resumes ordinary AUTO.

`KLS_DISABLE_SPARSE_SYMMETRIC_FRAGMENTED_METIS_POLICY=1` provides a generic
same-binary control. The leaf-count override is
`KLS_SPARSE_SYMMETRIC_FRAGMENTED_METIS_NDP_NPES`.
`sparse_symmetric_fragmented_metis_symbolic_eligible` and
`sparse_symmetric_fragmented_metis_policy_eligible` expose the accepted
symbolic and numeric stages through `kls_stats` and benchmark JSON.

The independent smoke fixture has 65,696 rows, outside the old window, and
selects the policy through both public sparse formats; a one-edge reciprocity
break and the disable switch reject it. An appended-scalar `ASIC_100ks`
metamorph at 99,702 rows and simultaneous relabelings also select, while a
larger 101,238-row fringe crosses the normalized boundary and falls back.
Across the target and three metamorphs, 1,200 independently checked entrywise
generations through 10% amplitude had worst relative-L2 error `5.91e-15`.
Twenty alternating parent/current target pairs put the modeled-cycle ratio at
`0.9926` geometrically (`0.9932` by the means); twelve extension pairs measured
`0.8525`. Release and ASan/UBSan/LSan CTest pass all four tests, including
leak-enabled target, positive-extension, and boundary-rejection runs.

## Sparse-spike predicted lifecycle

The former `nxp1` H100 route no longer recognizes a 414,000--415,000-row,
2.64--2.67-million-entry benchmark box. Under the standard eight-thread AUTO
contract, it now starts from a bounded topology proposal: a 200,000--750,000-
row nearly diagonal graph with six--eight entries per row and a moderate row
or column spike. The proposal reuses the topology census already required by
AUTO rather than adding another full sparse-matrix pass.

Later capabilities have independent measured gates. The accepted factor must
be a matched, unscaled, one-block predicted METIS numeric with complete
separator coverage, no pivot repair, and normalized fill/work bounds before
KLS skips the scale trial and selects direct EGraph updates. Retaining the
cluster forest, full worker width, and lower relaxed-consume floors additionally
requires a nearly all-private, low-pipeline separator and tighter fill/work
bounds. A topology or factor that fails either measured stage follows the
ordinary AUTO schedule.

`KLS_DISABLE_SPARSE_SPIKED_PREDICTED_POLICY=1` is the generic same-binary
control. The candidate, factor, and clustered decisions are exposed through `kls_stats`
and benchmark JSON as `sparse_spiked_predicted_candidate`,
`sparse_spiked_predicted_factor_eligible`, and
`sparse_spiked_predicted_clustered_eligible`.

Validation includes an independent 200,000-row sparse-spike fixture and its
same-size, same-density unspiked control, a 1,024-row `nxp1` extension beyond
the old order window, and a simultaneous relabeling that proposes the topology
but is rejected by the measured factor gate. Nine hundred positive entrywise
generations through 10% amplitude had worst relative-L2 residual
`5.33e-11`; the measured rejection remained accurate through another 32
generations. Six alternating quiet-core target pairs put the complete modeled-
cycle ratio at `1.0074` versus the former exact policy, while the extension and
independent fixture measured `0.7016` and `0.6100`. Release and
ASan/UBSan/LSan CTest pass all four tests, including leak-enabled positive and
measured-rejection runs.

## Dense fragmented scaled-row lifecycle

The former `TSOPF_RS_b2383` row-update route no longer recognizes exactly
38,120 rows and 16,171,169 entries. It now makes no input-name, dimension,
orientation, or ordering decision. After the ordinary AUTO factor completes,
KLS admits direct cooperative-row updates only for a real scaled, full-rank BTF
numeric whose component count, largest component, pivot count, balanced fill,
work, and dense-input ratio all fall inside normalized bounds. A failed factor
gate retains the generic update engine.

`KLS_DISABLE_DENSE_FRAGMENTED_SCALED_ROW_POLICY=1` provides a generic
same-binary control. `dense_fragmented_scaled_row_factor_eligible` exposes the
measured decision through `kls_stats` and benchmark JSON.

The smoke test supplies an independent 4,096-row family with 16 dense
256-row components. It selects normal AMD rather than the target's transposed
METIS factor, works through both CSC and CSR, refactors changed values through
the row engine, and recovers an independently constructed solution. A
same-order, same-entry-count control merges the pattern into four 1,024-row
components and rejects. A 512-component TSOPF extension at 38,632 rows also
selects outside the old exact identity.

Across both real value sets, the extension, and the independent family, 1,200
independently checked entrywise generations through 10% amplitude had worst
relative-L2 residual `2.04e-11`. Six alternating target pairs put the modeled-
cycle ratio at `1.0021` versus the former exact policy. Four extension pairs
measured `0.4892`, and six independent-family pairs measured `0.1756`.
Release and ASan/UBSan/LSan CTest pass all four tests; leak-enabled target,
extension, and independent runs are clean.

## Sparse full-diagonal METIS row lifecycle

The former `mc2depi` row route and `G3_circuit` route no longer recognize,
respectively, a narrow benchmark box and one exact `(n, nnz)` pair. Under the
standard eight-thread AUTO contract, both now begin with one reusable topology
proposal: 100,000--4,194,304 rows, three--five stored entries per row, a full
structural diagonal, no empty row, and maximum row and column degree eight.
The existing AUTO topology scan supplies the verdict, so the larger class adds
no second full sparse-matrix pass.

Ordering resources and later capabilities are staged. Fine sparse candidates
retain an `n/608` constrained-AMD window. Candidates above 1,048,576 rows and
4.25 entries per row use about 37 windows per worker, rounded and bounded by
generic resource limits. The retained no-BTF METIS symbolic must then prove
full rank, balanced normalized and absolute fill, and a complete, almost
entirely private separator. The fine symbolic selects the measured `1e-6`
pivot tolerance; the high-resource symbolic retains the requested `1e-3`.

The actual fixed-pivot numeric must independently remain unscaled and
balanced, avoid pivot repair, and satisfy its class's normalized fill and work
bounds before KLS enables predicted-row preparation and recurring cooperative-
row updates. The final proof applies to a predicted bootstrap or a later KLU
replacement. A rejection at any stage resumes ordinary AUTO behavior.

`KLS_DISABLE_SPARSE_FULL_DIAGONAL_METIS_ROW_POLICY=1` provides a generic
same-binary control. The candidate, symbolic, and factor decisions are exposed through `kls_stats` and
benchmark JSON as
`sparse_full_diagonal_metis_row_candidate`,
`sparse_full_diagonal_metis_row_symbolic_eligible`, and
`sparse_full_diagonal_metis_row_factor_eligible`.

The smoke test now has independent 100,000-row cubic and quintic circulants.
Both propose the topology while deliberately failing the absolute symbolic
economics; CSC and CSR cover the cubic family, and the quintic family exercises
the newly admitted five-entry density. Removing one diagonal or setting the
generic or either legacy disable rejects the proposal. Screening all 110 local
SuiteSparse inputs finds only four complete topology proposals: `mc2depi` and
`G3_circuit` reach all measured stages, while `G2_circuit` is below the
absolute setup floor and `ss1` exceeds normalized fill.

The earlier `mc2depi` extensions remain positive, including a measured
`1/1/0` relabeling with a real pivot nudge. A new 4,096-row `G3_circuit`
extension outside the exact identity and an adjacent relabeling both report
`1/1/1`. Across the three high-resource matrices, 180 independently checked
entrywise generations through 10% amplitude had worst relative-L2 residual
`5.94e-13`; combined with the fine class, 1,080 positive generations remain
below `2.04e-12`.

Four alternating `G3_circuit` target pairs preserve predicted-first/row
execution and put the complete modeled-cycle ratio at `0.9816` versus the
former exact selector. Three extension pairs change the parent's EGraph
fallback to measured row updates and reduce the cycle to `0.5293`; final
residual improves from as high as `4.86e-13` to `4.21e-16`. Release and
ASan/UBSan/LSan CTest pass all four tests, and leak-enabled target, extension,
and relabeling runs are clean.

## Giant symmetric scalar-fringe METIS row lifecycle

The former `rajat31` AUTO/eight-thread exception no longer recognizes a narrow
4.68--4.70 million-row and 20.2--20.4 million-entry box. Its proposal is now
an exact topology proof: 1,048,576--8,388,608 rows, three--six entries per
row, exact multiplicity-aware structural symmetry, an almost-full diagonal,
a missing-diagonal scalar fringe between `n/8192` and `n/1024`, at most
`n/128` scalar columns overall, and one to sixteen bounded hubs. Sorted CSC
uses allocation-free reciprocal searches; unsorted input uses an exact
transpose comparison. CSC and CSR therefore receive the same verdict.

Topology may start an overlapped seven-leaf NodeNDP proposal and select scale
`-1`, pipelined first-factor routing, and narrow-panel suppression. It does not
authorize recurring row machinery. The retained METIS/BTF symbolic must
independently prove full structural rank, one giant block with a mostly scalar
normalized fringe, 32--80 estimated factor entries per row with an absolute
128--512 million-entry setup floor, balanced L/U, and a complete separator
with at least 99% private rows. Unknown predicted-symbolic flop estimates are
allowed because the numeric stage measures work directly.

The actual fixed-pivot numeric must then remain unscaled at the requested
`1e-3` tolerance, have no row permutation, scale vector, off-diagonal pivot,
nudge, or perturbation, retain balanced 40--64-entry-per-row fill, and perform
16,384--32,768 operations per row. Only this cached factor verdict enables
predicted-row preparation, direct cooperative row updates, large dense-group
priority, and packed-solve publication. Numeric replacement invalidates and
recomputes the verdict.

Set
`KLS_DISABLE_GIANT_SYMMETRIC_SCALAR_FRINGE_METIS_ROW_POLICY=1` for a generic
same-binary control. Candidate, symbolic, and numeric decisions are exposed in `kls_stats`
and benchmark JSON as
`giant_symmetric_scalar_fringe_metis_row_candidate`,
`giant_symmetric_scalar_fringe_metis_row_symbolic_eligible`, and
`giant_symmetric_scalar_fringe_metis_row_factor_eligible`.

The independent smoke family has a 1,048,576-row symmetric quintic core, one
517-entry hub, and 512 missing-diagonal scalar leaves. It reports `1/0/0` in
CSC, CSR, and deliberately unsorted CSC. Breaking one reciprocal edge while
preserving dimensions, density, degree counts, diagonal coverage, triangular
counts, and sorted storage rejects the proposal; the generic and legacy
controls reject it as well. Screening every local SuiteSparse input capable of
meeting the coarse resource bounds leaves only `rajat31` as a proposal.

Appending 10,240 independent diagonal blocks produces a 4,700,242-row holdout
outside the old dimension box. It reports `1/1/1`, grows the measured BTF
fringe from 2,502 to 12,742 rows, retains predicted-first/row execution, and
returns a worst `6.91e-12` relative-L2 residual across ten independently
checked 10%-amplitude updates. A 65,536-adjacent-swap relabeling also reports
`1/1/1`; its measured row update plus packed solve is faster than the parent's
factor-rejected EGraph update plus solve despite the relabeling's much slower
cold layout.

Three alternating target pairs preserve identical `1.81e-14` residuals and
have geometric current/parent ratios of `1.0244` for analysis, `0.9968` for
initial factor, `0.9921` for first refactor, `0.9152` for steady refactor,
`1.0005` for packed solve, and `0.9387` for the complete modeled cycle. Two
extension pairs reduce the modeled-cycle ratio to `0.5803`. Release and
ASan/UBSan/LSan CTest pass all four tests; leak-enabled target and extension
factor/refactor/solve runs are clean.

## Hybrid huge-single EGraph lifecycle

The former `G2_circuit` recurring scheduler no longer recognizes a
150,000--150,200-row, 726,000--727,500-entry identity followed by absolute
fill and work boxes. The replacement is a cached factor capability. It accepts
normal METIS factors with 32,768--524,288 rows, three--eight input entries per
row, one full-matrix component, no scaling transform or pivot repair, balanced
L/U streams totaling 40--128 entries per row, and 3,072--32,768 measured
operations per row. A no-BTF or predicted symbolic may retain the unknown
structural-rank sentinel; an observed rank must be full.

The separator is independently checked before the capability is admitted. It
must cover every row and component, expose at least one private component per
worker, place at least 31/32 of rows in private components, keep total pipeline
rows below 1/32, and bound the largest private and pipeline components by
`n/4` and `n/48`. Only this settled representation receives the hybrid
clustered-prefix/dependency-pipeline schedule, full-width EGraph dispatch, and
the associated settled floor/panel decisions. Numeric replacement invalidates
and recomputes the verdict.

`KLS_DISABLE_HYBRID_HUGE_SINGLE_EGRAPH_POLICY=1` provides a same-binary
control. `hybrid_huge_single_egraph_factor_eligible` exposes the accepted
factor through `kls_stats` and benchmark JSON. The old G2-only solve-contract
bypass is removed: predicted factors now run the ordinary first-solve residual
probe. On the target, that probe measured a `6.48e-11` maximum residual,
certified the numeric, and allowed later solves to use the normal clean-factor
fast path.

The smoke test constructs an independent 32,768-row five-point grid under
explicit METIS/no-BTF analysis, refactors changed values through EGraph, and
checks its solution entrywise. A same-order, same-density narrow strip has
lower measured fill/work and rejects; the master-disable control also rejects.
A 4,096-node reciprocally coupled extension of `G2_circuit` has 154,198 rows,
outside the former identity, and reports eligible with a 36-level prefix and
18,665-column dependency pipeline. Adjacent and full simultaneous relabelings
also retain the capability from their measured factors.

Six hundred independently verified entrywise generations across the target
and coupled extension at 0.1%, 1%, and 10% amplitude had worst relative-L2
residual `8.37e-13`. Eight alternating target pairs put the median paired
current/parent modeled-cycle ratio at `1.0035`; six extension pairs put the
median at `0.9462`, with first and steady-refactor medians of `0.9369` and
`0.9325`. Six same-binary grid pairs put the enabled/disabled modeled-cycle
ratio at `0.7951` geometrically (`0.7749` median) and the steady solve ratio at
`0.3252`. Release and ASan/UBSan/LSan tests pass, and leak-enabled target,
extension, and shuffled-relabeling runs are clean.

## Promoted-tolerance relative-L2 recovery

A solve-accuracy rule formerly recognized one matrix through exact dimension,
entry-count, BTF-block, ordering, orientation, worker-count, and input-format
coordinates. The rule is now attached to the numeric fact it needs: a plain
factor retained a positive pivot threshold below the caller's requested
threshold, and the solver retained current matrix values with which to measure
the true residual. For an out-of-place, single-RHS solve in the factor's
untransposed kernel frame, KLS accepts only a finite relative-L2 residual at or
below `1e-9`. It damps a correction only after measured stagnation and, if
stationary refinement cannot meet the contract, rebuilds at monotonically
stronger thresholds derived from the selected and requested thresholds. A
recovery solve retains a `5e-9` limit, still a 2x margin below the public
`1e-8` validity line.

The contract is independent of CSC versus equivalent CSR input, public solve
orientation, ordering, BTF geometry, and worker count. Recovery state persists
for the installed numeric even when its stronger threshold no longer compares
below the original request, and every changed-value refactor refreshes the
matrix values used by the residual. Set
`KLS_DISABLE_PROMOTED_TOLERANCE_L2_RECOVERY=1` for a same-binary control.
`promoted_tolerance_l2_recovery_eligible`,
`promoted_tolerance_l2_contract_run_count`, and
`promoted_tolerance_l2_recovery_count` expose the capability and its activity
through `kls_stats` and benchmark JSON.

The smoke regression uses an independently generated 150,000-row family with
150 strongly connected components, five entries per row and column, and only
10% diagonal coverage. It selects `1e-4` from a requested `1e-3`, verifies
changed-value CSC and equivalent-CSR/public-transpose solves with four workers,
and checks the live-factor disable. Appending 1,000 independent diagonal rows
to the SuiteSparse target moves it outside every former dimension, entry, and
BTF-block box while retaining eligibility. Sixty independently checked
entrywise generations per matrix across 0.1%, 1%, and 10% amplitudes had worst
relative-L2 residuals `7.43e-10` and `9.38e-10`; none required a recovery
rebuild. Leak-enabled ASan/UBSan/LSan factor/refactor/solve runs cover both
matrices.

## Bounded-degree retained-preconditioner lifecycle

The former `Hamrle3` lifecycle no longer recognizes a narrow order and entry-
count window or requires exactly eight workers. Under the ordinary parallel
AUTO/BTF/static-pivoting contract, a proposal now comes from sparse topology:
32,768--4,194,304 rows, two--six stored entries per row on average, no empty
row or column, maximum row and column degree 16, at most `n/32` scalar rows or
columns, and at most `n/128` columns with a structural diagonal. The proposal
uses transpose AMD with scale 1 and a `1e-4` pivot threshold, but cannot force
that representation by itself.

The measured BTF symbolic must prove full structural rank, at most `n/128+1`
blocks, a dominant block containing at least 31/32 of the matrix, 48--512
estimated factor entries per row, and at least 8,192 estimated operations per
row. A rejected proposal is discarded and ordinary two-sided AUTO analysis
resumes. The installed numeric must independently preserve the representation
and rank, avoid nudges and perturbations, keep L and U within 4x balance, retain
48--1,536 measured factor entries per row, and perform at least 8,192 measured
operations per row.

Only a numeric that passes all three stages captures its factor-time values.
A changed-value update may reuse that factor when every entry changes by no
more than `1.01e-3*abs(old)+64*DBL_MIN`. Every resulting solve checks the true
current-matrix relative-L2 residual and must reach `1e-9`; larger updates or a
reference allocation failure permanently return that numeric epoch to an
ordinary refactor. This makes retained-factor accuracy a measured contract,
not an inference from matrix identity.

Set
`KLS_DISABLE_BOUNDED_DEGREE_RETAINED_PRECONDITIONER_POLICY=1` for a generic
same-binary control. `kls_stats` and benchmark JSON expose
`bounded_degree_retained_preconditioner_candidate`,
`bounded_degree_retained_preconditioner_symbolic_eligible`,
`bounded_degree_retained_preconditioner_factor_eligible`, and
`bounded_degree_retained_preconditioner_reuse_count`.

The independent smoke fixture is a 32,768-row six-neighbor toroidal graph with
no diagonal and generated values. It qualifies at candidate, symbolic, and
factor stages under four workers, reuses the numeric for changed values, and
checks an independently formed solution. A same-order six-neighbor narrow-band
ring passes the input proposal, fails the symbolic economics, and completes
ordinary AUTO fallback. A scan of all 110 local SuiteSparse matrices leaves
only `Hamrle3` inside the coarse topology bounds; that is compatibility
evidence, while the independent fixture and extensions establish the
capability's generality.

Appending 1,000 independent diagonal blocks produces a 1,448,360-row holdout
outside the former order window. It reports `1/1/1`, retains a 1,447,360-row
dominant block, and completes 20 checked updates with worst relative-L2
residual `1.90e-10`. A 1% update rejects reuse, selects EGraph, and returns
`1.36e-15` with reuse count zero. Three pinned alternating target pairs keep
the same representation and residual while measuring geometric current/parent
ratios of `1.0093` for analysis, `1.0008` for steady update, `0.9997` for its
solve, and `1.0120` for the complete modeled cycle. The recurring path is
therefore effectively unchanged; the remaining cycle variation comes from the
layout-sensitive one-time factor. Release and leak-enabled
ASan/UBSan/LSan CTest pass all four tests; an explicit sanitized target
factor/update/solve also reports `1/1/1`, one retained reuse, and no findings.

## Asymmetric bounded-degree direct-METIS lifecycle

The direct NodeNDP selectors for Freescale1 and memchip no longer compare an
input with either matrix's exact order or entry count. Under the standard
AUTO/8-thread/BTF/static-pivoting contract, a proposal now comes from directed
topology: 131,072--8,388,608 rows, four--six stored entries per row on average,
no empty row or column, maximum in/out degree 32, an almost-full structural
diagonal, a bounded scalar population, either unequal strict-triangle counts as
a fast asymmetry proof or an exact bounded-degree reciprocity test, and a raw
strongly connected component covering at least 127/128 of the graph. There may
be at most `n/1024+1` raw SCCs.

Raw SCC fragmentation chooses tuning rather than matrix identity. At most
`n/8192+1` SCCs is thin-fringe class 2, which uses fourteen deterministic
NodeNDP leaves and 3,072-column constrained-AMD windows. A coarser fringe is
class 1, which uses eight leaves and the ordinary giant-window rule. The SCC
pass is invariant under transpose and simultaneous relabeling.

Topology only proposes direct METIS. The resulting symbolic must have at least
1,048,576 rows, full structural rank, BTF enabled, a 127/128 dominant block,
the same SCC-fragmentation class, a complete separator covering every row,
at least 127/128 private rows, bounded private/pipeline components, and—when
estimates exist—8--64 factor entries and 256--65,536 operations per row. A
rejected proposal is freed and the complete AUTO tournament resumes. The
installed numeric independently requires an unscaled, untransformed,
unperturbed fixed-pivot factor with balanced L/U streams, 8--64 measured
factor entries per row, and 256--16,384 operations per row.

Only the numeric-accepted capability receives the recurring giant-chain row
updates and setup suppressions. Thin-fringe numerics still use the deterministic
pipeline for first-pattern discovery, but their changed-value row factor is
published to the packed column solve through the existing factor-level
preference. This replaces the former hard-coded no-publish solve shortcut;
the coupled holdout showed that shortcut was not a generic SCC property.

Set
`KLS_DISABLE_ASYMMETRIC_BOUNDED_DEGREE_DIRECT_METIS_POLICY=1` for a generic
same-binary control. The NodeNDP leaf-count override is
`KLS_ASYMMETRIC_BOUNDED_DEGREE_DIRECT_METIS_NDP_NPES`. `kls_stats`
and benchmark JSON expose
`asymmetric_bounded_degree_direct_metis_candidate`,
`asymmetric_bounded_degree_direct_metis_tuning_class`,
`asymmetric_bounded_degree_direct_metis_symbolic_eligible`, and
`asymmetric_bounded_degree_direct_metis_factor_eligible`.

Independent smoke fixtures use directed circulant cores at 131,072, 133,120,
and 135,168 rows. They cover both tuning classes, a nearby-size coupled
extension, equivalent CSR/transposed topology, symmetric and degree controls,
symbolic rejection with ordinary AUTO fallback, and both disable switches.
Freescale2 rejects on fragmentation/degree and structurally symmetric
`circuit5M_dc` rejects at the exact asymmetry stage.

Real coupled holdouts also move both targets outside their former identities.
Freescale1 plus 1,024 reciprocal diagonal nodes reports `1/1/1`, preserves the
class-1 dominant core, and measures a 33.64-second modeled cycle versus 35.87
seconds for the frozen exact parent. Memchip plus 512 nodes reports class 2 and
`1/1/1`, measuring 39.43 versus 42.58 seconds. On the original matrices, a
pinned parent/current pair is 31.08 versus 30.94 seconds for Freescale1. Two
alternating memchip pairs give a geometric current/parent cycle ratio of
`0.8246`: publishing raises steady row-update time by about 20% but cuts the
steady solve to `0.2325x`, improving the complete workload. All verified
target and holdout residuals are below `4.13e-16`. Release and leak-enabled
ASan/UBSan/LSan CTest pass all four tests.

## Near-symmetric mega-hub AMD lifecycle

The former `circuit5M` selector no longer compares the input order and entry
count with a SuiteSparse coordinate. Under the standard AUTO, eight-worker,
BTF, static-pivoting contract, a topology proposal accepts 131,072--8,388,608
rows and eight--sixteen stored entries per row. Every row and column must be
nonempty, almost every column must contain its diagonal, scalar rows and
columns are bounded, and the maximum in- and out-degree must each lie between
`n/8` and `n/2`. A nonzero but bounded set of in/out-degree mismatches proves
directed asymmetry without depending on CSC order, transpose, or vertex
labels. A column-pointer pass rejects inputs with no plausible hub before the
classifier allocates row degrees or scans all nonzeros.

Topology only proposes normal AMD/BTF. The resulting symbolic must have at
least 1,048,576 rows, full structural rank, two through `n/1024+1` blocks, a
proper dominant core covering at least 1023/1024 of the matrix, and a fringe
between `n/4096` and `n/512+8` rows. Balanced symbolic L/U fill must total
32--64 entries per row and estimated work must be 512--4,096 operations per
row. If any measurement fails, KLS discards the trial and resumes the complete
AUTO tournament.

The installed numeric independently preserves the accepted symbolic identity,
normal AMD/BTF representation, caller pivot threshold, and max-row scaling.
It must remain nonsingular and have no row permutation, nudge, or perturbation,
few off-diagonal pivots, balanced L/U streams totaling 8--24 entries per row,
and measured work of 64--512 operations per row. Only that accepted numeric
suppresses redundant scale, METIS, predicted-pattern, and diagonal-equivalent
trials and retains the row-scale vector across changed-value EGraph refactors.

Set `KLS_DISABLE_NEAR_SYMMETRIC_MEGA_HUB_AMD_POLICY=1` for a generic
same-binary control. `kls_stats` and benchmark JSON expose
`near_symmetric_mega_hub_amd_candidate`,
`near_symmetric_mega_hub_amd_symbolic_eligible`, and
`near_symmetric_mega_hub_amd_factor_eligible`.

The independent smoke family uses generated 131,072- and 135,168-row
circulant cores, a quarter-order reciprocal hub, one unmatched directed edge,
and disconnected two-node fringe components. It covers equivalent CSR,
unsorted storage, a nearby size, exact-symmetric and oversized-hub negatives,
a connected proposal that fails the symbolic stage, and both disable
switches. These compact fixtures deliberately lie below the production
symbolic floor, so they prove invariant proposal and fallback behavior rather
than serving as a second large factor-positive family.

The development manifest adds `FullChip`, `Freescale2`, and `circuit5M_dc` as
near-family controls. The cross-family holdout manifest adds `kkt_power`,
`CurlCurl_3`, `StocF-1465`, `Transport`, and `wikipedia-20051105`; all reject
at the input stage for independent diagonal, hub, degree-balance, or density
reasons. `FullChip` also retains its existing METIS/no-BTF symbolic exactly;
a final frozen/current pair measured 7.107 versus 7.066 seconds of analysis,
showing that the early hub rejection removes material classifier overhead.

Appending 512 generated coupled nodes to the positive creates a 5,558,838-row
holdout outside the old identity. It reports `1/1/1`, retains a 5,556,264-row
dominant BTF core, and preserves identical 33,760,986/33,760,248 L/U counts.
With two factors and three independently checked refactors, its modeled cycle
falls from 53.59 seconds in the frozen exact parent to 42.01 seconds; the
worst relative-L2 residual is `1.26e-14`. On the original input, two pinned
alternating pairs keep identical factors and residuals, analysis within 0.3%,
and a geometric current/parent modeled-cycle ratio of about 1.03, within the
observed cold-factor/first-refactor variation. Thus the exact identity is gone
without materially changing the established workload, while the moved-size
positive demonstrates useful behavior outside the benchmark coordinate.

## Giant dominant-hub METIS dense-tail lifecycle

The former `FullChip` route no longer compares the input order and entry count
with a SuiteSparse coordinate. Under the standard AUTO, eight-worker, BTF,
static-pivoting contract, a topology proposal accepts 131,072--8,388,608 rows
and eight--twelve stored entries per row. Every row and column must be nonempty,
almost every column must contain its diagonal, and scalar rows and columns are
bounded. The maximum in- and out-degree must each lie between `2n/3` and
`7n/8`; a nonzero but bounded set of degree mismatches supplies a directed-
asymmetry proof. The predicates are invariant under transpose, simultaneous
relabeling, and stored entry order. A column-only pass rejects unrelated
matrices before allocating row degrees or scanning all entries.

Topology is only the first stage. The selected symbolic must contain at least
1,048,576 rows, use normal METIS without BTF, and consist of one full-size
block. Its balanced L/U estimate must total 32--96 entries per row. The
separator analysis must cover every row with at least two workers, at most four
components per worker, at least one pipeline component, at least 255/256 of
the rows in private components, and bounded private and pipeline maxima. An
unknown structural rank is accepted because no-BTF KLU does not compute one;
when rank or work estimates are available they must satisfy the contract.

The installed numeric supplies the final evidence. It must preserve that
symbolic identity, normal METIS/no-BTF representation, scale `-1`, and the
caller threshold; come from the pipelined KLU path with a 1,024--16,384-column
dense tail; have full numerical rank, bounded off-diagonal pivots and nudges,
no perturbation, balanced measured fill of 32--96 entries per row, and measured
work of 4,096--131,072 operations per row. Only the accepted numeric suppresses
the redundant automatic scale retry. Its verdict is cached across in-place
EGraph refactors and invalidated whenever the numeric or symbolic is replaced.

Rank completion is deliberately outside this automatic policy. `FullChip` is
exactly singular, and `KLS_ENABLE_SINGULAR_COMPLETION=1` still must be supplied
by the caller before KLS may constrain its zero-pivot degrees of freedom. The
general capability only routes that explicitly requested operation through
the pipelined dense-tail factor and reuses its checked zero-pivot discovery; it
does not silently change any ordinary singular problem.

Set
`KLS_DISABLE_GIANT_DOMINANT_HUB_METIS_DENSE_TAIL_POLICY=1` for a generic
same-binary control. Route controls are
`KLS_DISABLE_GIANT_DOMINANT_HUB_ROUTED_FACTOR`,
`KLS_DISABLE_GIANT_DOMINANT_HUB_ZERO_DISCOVERY`, and
`KLS_DISABLE_GIANT_DOMINANT_HUB_COMPLETION_PIPE`; the dense-tail override is
`KLS_GIANT_DOMINANT_HUB_DENSE_TAIL`. `kls_stats` and analyze/full benchmark JSON
expose `giant_dominant_hub_metis_dense_tail_candidate`,
`giant_dominant_hub_metis_dense_tail_symbolic_eligible`, and
`giant_dominant_hub_metis_dense_tail_factor_eligible`.

The independent smoke family uses 131,072- and 135,168-row reciprocal
nine-point circulants, a three-quarter-order reciprocal hub, and one unmatched
directed edge. Equivalent CSR, unsorted CSC, and nearby-size fixtures pass the
input stage; exact symmetry, half-order and nine-tenths-order hubs, and both
disable variables reject. These fixtures deliberately remain below the
million-row symbolic floor and therefore test invariant proposal and ordinary
fallback behavior rather than claiming an independent large numeric positive.

Separate development and cross-family holdout manifests add `circuit5M`,
`Freescale2`, `G3_circuit`, `rajat30`, `kkt_power`, `CurlCurl_3`, `StocF-1465`,
`Transport`, and `wikipedia-20051105` as controls. In live checks,
`G3_circuit`, `circuit5M`, and independently sourced `CurlCurl_3` all report
`0/0/0`, despite collectively covering large one-block METIS, similar density,
and dominant-hub shapes.

Appending 512 coupled nodes to the motivating matrix produces a
2,987,524-row, 26,623,519-entry holdout outside the old exact coordinate. It
reports `1/1/1`, retains normal METIS/no-BTF, factors in 14.74 seconds, and
solves with `1.97e-14` relative-L2 residual. The frozen exact parent did not
finish that factor within 180 seconds. On the original matrix, two alternating
pinned parent/current pairs retain identical 82,748,400/84,816,974 L/U counts
and residuals below `1.15e-14`; the generalized input proof adds roughly
0.1 second of analysis, while initial factors remain in the same 10--12-second
band. Refactor samples on the shared host remain more variable than the policy
change, so they are treated as dispersion rather than evidence for a new
kernel speedup. The limitation is explicit: the only full-lifecycle positive
outside the original coordinate is a coupled metamorphic extension; the
cross-family matrices are rejection controls, not natural positives.

## Fragmented-chain one-update selector retired

The residual Freescale H100 selector has been removed rather than widened. It
accepted any AUTO/8-thread input with 2.5--3.7 million rows and 12--18 million
stored entries, then forced unscaled fourteen-leaf METIS and a specialized
first-update lifecycle. Once Freescale1 and memchip had moved to the guarded
asymmetric bounded-degree capability, Freescale2 was the only SuiteSparse
matrix still inheriting this broad size/density box.

An attempted topology replacement identified the intended fragmented
dominant-BTF/pivot-boundary shape and also accepted a coupled moved-size
extension. The full horizon audit nevertheless falsified the optimization.
On 100 entrywise 0.1% updates, with every generation independently checked,
the old enabled path and the ordinary same-binary fallback measured:

| Metric | old specialized path | ordinary generic path |
| --- | ---: | ---: |
| analysis | 2.643 s | 2.662 s |
| initial factor | 7.877 s | 6.520 s |
| first refactor | 1.203 s | 4.411 s |
| steady refactor | 0.5730 s | 0.08190 s |
| first changed solve | 0.3789 s | 0.08366 s |
| steady changed solve | 0.2954 s | 0.04503 s |
| modeled 100-state SPICE cycle | 97.454 s | 26.160 s |
| worst relative-L2 update residual | `3.0466e-9` | `1.6377e-13` |

The specialized path won only the first refactor. Its 39,402,413/25,595,605
L/U entries, 2.528 billion measured operations, and 5,437 off-diagonal pivots
made the recurring pool update and solve much slower than the generic scale-2
factor's 14,782,212/14,655,388 entries, 376.9 million operations, and 1,696
off-diagonal pivots. The generic route costs only `0.268x` as much over the
actual horizon and is about four orders of magnitude more accurate. Promoting
the staged classifier would therefore have generalized a one-update benchmark
artifact; retiring the selector generalizes behavior instead.

Freescale2 now follows ordinary AUTO selection: normal METIS, scale 2, mapped
steady updates, and no family-specific race, row-refactor, scale, or deferred
setup inheritance. Its singular semantics are unchanged: without
`KLS_ENABLE_SINGULAR_COMPLETION=1` setup still returns the singular-matrix
error. Freescale1 and memchip retain scale `-1` only after their generic
direct-METIS symbolic capability is accepted. Frozen-parent/current checks
preserve their respective class 1/class 2 verdicts, exact L/U geometries
(24,801,399/24,801,399 and 29,216,967/29,216,949), row-refactor paths, and
residuals (`4.12e-16` and `7.37e-17`). No SuiteSparse coordinate replaces the
deleted box. Release and leak-enabled ASan/UBSan/LSan CTest pass all four
tests.

## High-work tiny-fringe BTF/PTS factor lifecycle

The former `ASIC_320k` post-factor selector is now a normalized resource
capability.  The old selector required scale 0, but the current AUTO factor of
the motivating matrix selects scale `-1`; instrumenting all 24 calls across
three changed-value updates showed that the old policy was dormant.  The new
gate therefore does more than rename an active benchmark exception.

Under the standard AUTO, eight-worker, BTF, static-pivoting contract, the
installed numeric must use normal AMD or AMF with no scale vector, cover
131,072--1,048,576 rows at five--seven stored entries per row, and have full
structural rank.  Its dominant BTF block must leave a fringe between `n/512`
and `n/256`; the remaining block count must be between one per sixteen fringe
rows and one per fringe row.  Both symbolic and measured L+U fill must be
balanced within 2x and total 10--16 entries per row.  Estimated and measured
work must each lie between 2,048 and 4,096 operations per row, with no pivot
nudge or perturbation and at most `n/512` known off-diagonal pivots.  A
predicted first factor may defer the last pivot count until it is installed.

That factor evidence authorizes only the retained decisions supported by the
component controls: skip redundant scale/METIS consultations, decline the
Algorithm-5 prefactor update, settle the EGraph width/floor probes, admit the
wider PTS forest and its 1.5 cut, select a validated PTS solve directly, and
store the solve permutations in 32-bit form.  The PTS builder still verifies
the actual elimination forest, shared top, and worker balance before the solve
can be selected.  ASIC-specific reciprocal and compact-stream metadata choices
were removed because isolated measurements showed no benefit.

Set `KLS_DISABLE_HIGH_WORK_TINY_FRINGE_BTF_PTS_POLICY=1` to disable this factor
lifecycle.  Because a valid large PTS solve plan is now an independent generic
capability, also set `KLS_DISABLE_VERIFIED_LARGE_PTS_SOLVE_POLICY=1` when the
comparison must restore its solve consultation.  Other component controls are
`KLS_DISABLE_HIGH_WORK_TINY_FRINGE_SETTLED_PROBES`,
`KLS_DISABLE_HIGH_WORK_TINY_FRINGE_PTS_CUT`, and
`KLS_DISABLE_HIGH_WORK_TINY_FRINGE_COMPACT_PERM`. `kls_stats` and both
benchmark JSON modes expose
`high_work_tiny_fringe_btf_pts_factor_eligible`.

The Sandia development family demonstrates independent rejections:
`ASIC_100k` and `ASIC_100ks` fall below the order floor and select METIS;
`ASIC_320ks` has one block; and the `ASIC_680k` variants have a much larger
fringe or select METIS.  Cross-family controls `scircuit`, `Raj1`, `transient`,
`rajat24`, `nxp1`, and `mac_econ_fwd500` all report zero for independent work,
fringe, block-count, fill, or scale reasons.  The corresponding development
and holdout manifests live under `bench/`.

Two generated positives remove the original coordinate.  A complete
simultaneous relabeling keeps 321,821 rows but changes AUTO from AMF to AMD and
changes measured fill/work to 2,120,493 L and U entries and 924.6 million
operations; it remains eligible and completes 100 checked updates with a
`3.71e-15` worst relative-L2 residual.  Appending 32,768 weakly coupled nodes
produces a 354,589-row, 2,030,132-entry matrix outside the old row ceiling.  It
selects AMD, has 105 blocks with a 353,691-row core, measures 2,290,549 L and U
entries and 1.152 billion operations, and also remains eligible.  In contrast,
appending the same number of independent diagonal blocks creates a 33,663-row
fringe and reports zero.  Appending 131,072 coupled nodes changes the selected
factor to a scale-2 single block and likewise rejects, showing that extension
alone is not sufficient.

Final pinned same-binary measurements use 100 deterministic entrywise 0.1%
updates and independently verify every generation.  Because refactor times on
the shared host have occasional multi-fold outliers in both arms, the table
reports medians rather than selected runs:

| matrix / policy | samples | steady refactor | steady changed solve | modeled 100-state cycle |
| --- | ---: | ---: | ---: | ---: |
| `ASIC_320k`, enabled | 7 | 13.88 ms | 3.65 ms | 2.99 s |
| `ASIC_320k`, disabled | 5 | 24.67 ms | 6.04 ms | 4.28 s |
| coupled +32,768, enabled | 6 | 16.29 ms | 4.20 ms | 2.90 s |
| coupled +32,768, disabled | 5 | 21.32 ms | 6.56 ms | 3.82 s |

Factor geometry is identical between enabled and disabled arms.  All six
enabled moved-size sweeps stay below `2.97e-15`; six of seven enabled target
sweeps stay below `2.80e-15`.  One target sweep reached `2.85e-9`, matching a
rare EGraph/refactor residual excursion also observed in frozen and disabled
controls, so it is recorded as shared variability rather than hidden from the
result.  The factor envelope remains deliberately narrow: the positives
outside the original coordinate are metamorphic, while the natural
SuiteSparse matrices reject that full factor lifecycle.  The solve-level
generalization below does not widen these factor thresholds; it starts from
the independently verified PTS plan that any installed numeric can build.

## Verified large-PTS solve adoption

Direct PTS solve selection no longer needs a matrix-family or factor-family
label once the generic PTS builder has validated the actual installed factor.
The builder checks ordered L/U streams, constructs and verifies an elimination
forest against every dependency, cuts independent subtrees, balances them over
the available workers, and bounds the shared ancestor top.  A solve-valid plan
whose dominant block has at least 16,384 columns is therefore selected
directly.  Valid smaller plans retain their existing measured family gates;
all other factors retain the two-serial/two-PTS runtime consultation.

This is intentionally a solve capability, not another approximation to the
`ASIC_320k` factor profile.  It covers normal AMD, AMF, and METIS factors,
single-block and fragmented BTF factors, scaled and unscaled numerics, and
different pivot regimes whenever their realized forest passes the same proof.
The refactor PTS gate remains separate and still uses its stricter
flop-weighted serial-top test.  Compact permutations were also tested as a
possible companion generalization, but a natural `epb3` A/B improved only
about 0.7%, so that behavior remains under its existing narrow policies.

`KLS_DISABLE_VERIFIED_LARGE_PTS_SOLVE_POLICY=1` restores the four-solve runtime
consultation for same-binary comparisons.  The appended
`kls_stats.verified_large_pts_solve_policy_eligible` field and both benchmark
JSON modes expose the verdict after solve metadata has been built (and report
zero in analyze-only mode).  Natural positives and structural controls are
recorded separately in
`bench/suitesparse_verified_large_pts_solve_positive_manifest.txt` and
`bench/suitesparse_verified_large_pts_solve_control_manifest.txt`.

The discovery pass covered 52 natural paper-union matrices and 75 matrices
from the generic development/holdout suites.  Every previously unclassified
solve-valid plan that completed the old independent timing consultation chose
PTS; none declined it.  That pass's positives span the Sandia family plus
`Raj1`, `G2_circuit`, `mc2depi`, `nxp1`, `ss1`, `scircuit`, and the natural
cross-family `Averous/epb3`.  Nearby controls such as `rajat21`, `rajat24`,
`transient`, `mac_econ_fwd500`, and `shallow_water1` stop at stream, forest,
shared-top, or plan-construction gates.

On `epb3`, the old consultation measured 1.724 ms for PTS versus 3.953 ms for
serial.  Five alternating 100-update same-binary samples put the median first
changed solve at 2.536 ms direct versus 5.360 ms with consultation, the steady
solve at 1.906 versus 1.933 ms, and the modeled cycle at 1.108 versus 1.123
seconds; every checked residual stayed near `5e-16`.  On `scircuit`, the old
probe selected PTS at 0.801 versus 1.927 ms, and direct adoption reduced the
median over 100 solves from 0.653 to 0.623 ms.  Three interleaved frozen-parent
and current `ASIC_320k` 100-update runs retained identical factor geometry and
`2.79e-15` residuals; modeled-cycle medians were 2.583 and 2.572 seconds,
respectively.  The relabeled and +32,768 coupled structural positives also
remain eligible and complete 100 checked updates below `3.72e-15`.

The remaining limitation is hardware portability rather than benchmark
identity: direct adoption relies on the conservative plan proof plus the
observed CPU corpus, while a radically different thread/runtime cost model may
still prefer consultation.  The master disable preserves that fallback.  The
change does not broaden factor selection, refactor routing, or small-plan
economics.

Moderately wide-top plans no longer need the former sparse-100K input box.
After the ordinary stream and forest proof, a previously unauthorized plan
may extend the solve top from 30% to 45% only when its block has at least
16,384 columns, its shared top spans at most 2% of those columns, it exposes at
least four independent chunks per worker, and its heaviest worker bin is at
most 1.20 times the average private work.  Existing family-authorized wide
plans are unchanged.  The resulting large plan still uses the same direct
verified-PTS policy; the stricter flop-weighted refactor gate remains separate.

The old 95,000--105,000-row/five--seven-entry sparse helper has been removed
from both PTS and METIS refinement.  An explicit METIS callback instead
recognizes a realized full-diagonal, bounded-hub core with balanced triangular
structure and scales its CAMD window with the core order.  A 105,390-row
coupled holdout, whose 105,003-column core is outside the old callback box,
independently satisfies both structural capabilities and completes 100 checked
updates below `3.06e-15`.  Five final AUTO H100 pairs reduced its median
modeled cycle from `1.4398s` to `1.3099s`.  On the original `ASIC_100ks`, 20
final AUTO samples preserved the exact factor and direct PTS route; the
generic/frozen medians were within 0.36%, while the 10%-trimmed mean favored
the generic build by 0.6%.

A predeclared twelve-matrix natural SuiteSparse probe set now supplies the
previously missing unrelated-family positive.  AUTO `Schenk_IBMNA/c-67` has a
57,293-column core, 3,595 independent chunks, 706 shared-top columns, 32.67%
top work, and an 8.42% heaviest bin.  It is outside the deleted input box and
selects the measured certificate directly.  Seven final alternating H100
pairs preserved the parent's exact 616,152/571,351-entry factor and reduced
the median changed solve from 0.732ms to 0.364ms; the paired modeled-cycle
ratio was 0.9436 and all 100-update residuals stayed below `1.61e-16`.  The
probe cohort is recorded in
`bench/suitesparse_balanced_wide_pts_natural_probe_manifest.txt`.

That same pass found a useful counterexample.  Explicit-METIS `rajat25`
initially passed the 1.25 private-bin bound with 39.09% top work but a 1.225x
heaviest-bin skew; direct PTS regressed nine paired modeled cycles by a 1.3093
median despite identical factors.  Tightening the generic balance proof to
1.20 rejects that direct-solve adoption while retaining `c-67`, `ASIC_100ks`,
and the coupled holdout.  Seven post-fix pairs restored `rajat25` to its
serial route and a 0.9936 paired cycle ratio.  Nine clean three-way target
rounds put the post-fix build at 0.9966 of the pre-tightening build and 0.9617
of the frozen parent, with identical `ASIC_100ks` factor geometry and PTS
eligibility.

The compact 1K circuit helper has also been removed from orientation,
tolerance, thread-count, and lean-refactor selection.  The ordinary small
AUTO policy now selects the faster normal frame for `orsirr_1`; eleven final
alternating H100 pairs measured a `0.8582` current/parent ratio.  The four
nearby natural controls (`circuit204`, `orsirr_2`, `rdb968`, and
`fpga_dcop_04`) remained within 0.6%, with unchanged factor decisions and
residuals.

## Resource-scaled portfolio completion and solve recovery

Three final AUTO gaps now use lifecycle evidence or a measured numerical
failure rather than matrix coordinates.  First, a matched pre-static proposal
whose estimated numeric storage exceeds a 32 MiB budget per requested worker
can skip its duplicate speculative serial factor when at least 16 advertised
refactors repay the match.  The ordinary post-match diagonal-strength scan is
still mandatory, and the normal predicted/pivoted first-factor path remains
the acceptance authority.  Set
`KLS_DISABLE_GENERIC_RESOURCE_SCALED_PRESTATIC_FIRST=1` to restore the
speculative factor for an A/B run.

Second, when nested dissection is only a near-tie lifecycle win in the first
AUTO orientation, KLS completes the independently realized ND portfolio in the
other orientation before comparing lifecycle scores.  This is limited by the
same ND admission and rollback contracts as the first candidate; it does not
select an orientation from dimensions, density, or matrix provenance.  Set
`KLS_DISABLE_GENERIC_NEAR_TIE_CROSS_ORIENTATION_ND=1` to retain the first-frame
shortcut.

Finally, an ordinary row-published factor or armed solve probe now verifies
the same relative-L2 quantity used by the benchmark validity rule, with a
`5e-9` acceptance line that leaves a factor-of-two margin below `1e-8`.
Stationary iterative refinement remains first.  If its measured residual
stalls, a cold restarted right-preconditioned GMRES recovery combines at most
eight applications of the installed LU and verifies the result against the
current stored matrix before returning it.  The recovery allocates its short
basis only after a real residual failure and is capability-gated to a plain
normal coordinate frame; transformed and transpose frames retain their
existing checked refinement paths.  The controls are
`KLS_DISABLE_ORDINARY_SELF_CHECK_L2_CONTRACT=1` and
`KLS_DISABLE_GMRES_SOLVE_RECOVERY=1`.

Plain AUTO factors whose current KLU diagonal-ratio estimate is below
`sqrt(DBL_EPSILON)` now receive one measured raw-solve certificate as well.
A clean first result settles that retained pivot family, so later changed
refactors do not inherit a residual-pass tax from a pessimistic `rcond` scalar.
If the raw solve and both LU-based correction methods fail, KLS uses a cold,
work-bounded LSQR recovery against the actual public CSC matrix. Twelve
two-sided norm-equilibration rounds precede the Krylov iteration, and every
accepted result is checked at the same `5e-9` relative-L2 margin. This is a
numeric outcome policy rather than a dimension, sparsity, ordering, or matrix
identity detector. `KLS_DISABLE_RCOND_SOLVE_CONTRACT=1` restores the old
uncertified plain solve, and `KLS_DISABLE_LSQR_SOLVE_RECOVERY=1` retains the
certificate but rejects after the cheaper recoveries fail. On
`fpga_dcop_01`, whose old full-precision LU answer had a `2.55e-2` residual,
all 20 rank-preserving generations now validate; the worst residual is
`4.99e-9`. The median five-pass H100 is 1.840 s. Both saved competitor runs
remain invalid on this numerically near-rank-deficient case (`3.48e-3` for CKTSO
and `2.43e-2` for SubtreeLU).

Focused eight-core entrywise-H100 measurements on the current machine give
the following development results.  They are one-pass engineering evidence,
not a replacement for the paper's counterbalanced 110-matrix campaign.

| case | prior generalized build | current build | accuracy / path result |
| --- | ---: | ---: | --- |
| `pre2` initial factor | 386.225 s | 9.005 s | relative residual `1.94e-16` |
| `mc2depi` H100 | 22.166 s | 17.861 s | worst generation `4.15e-13`; normal METIS + PTS |
| `FullChip` H100 | 237.224 s | 219.475 s | worst generation `3.23e-9`; stalled generation recovered in two GMRES directions |

The `mc2depi` result is below the saved CKTSO (`20.091 s`) and SubtreeLU
(`19.213 s`) measurements.  `FullChip` previously failed the paper's validity
rule at `7.72e-8`; the recovered run completes without a recovery refactor.
Five unrelated holdouts preserve their orientation and ordering decisions:
current/prior H100 ratios are `0.982` on `ASIC_100ks`, `1.007` on `Raj1`,
`0.967` on `rajat25`, `0.996` on `coupled`, and `0.980` on `G2_circuit`.
All remain below `2.35e-14` worst relative-L2 residual.  A seven-generation
ASan/UBSan `FullChip` replay exercised the two-direction GMRES path with no
sanitizer finding and the same `3.23e-9` externally audited maximum.

The large matched factor now applies the same resource argument to pivot and
scale speculation.  Tight matched pivots remain available when the estimated
numeric fits the workers' aggregate 32 MiB cache allowance, or when at least
five percent of the realized matched diagonals lie between the two pivot
thresholds.  Otherwise the ordinary `1e-3` factor is retained.  If that
native numeric is larger than aggregate cache but beats its symbolic storage
bound, an independent deterministic backward-residual probe can certify it
and avoid a second full matching-equilibration factor.  This uses realized
threshold pressure, factor storage, lifecycle length, and numerical evidence;
it does not use matrix dimensions or provenance.  The A/B controls are
`KLS_DISABLE_GENERIC_RESOURCE_SCALED_MATCHED_TOLERANCE=1` and
`KLS_DISABLE_GENERIC_RESOURCE_SCALED_NATIVE_MATCH_CERTIFICATE=1`.

Separator width is also part of the realized matched-ordering portfolio.  A
wide nested-dissection splitter remains the latency candidate.  When its
symbolic storage exceeds aggregate cache over a repeated lifecycle, KLS also
builds a two-leaf candidate and retains it only for a 0.5% symbolic-fill win
whose projected recurring saving repays a complete graph comparison.  The
normal numeric fill, pivot, conditioning, and residual gates still arbitrate
the selected ordering.  Set
`KLS_DISABLE_GENERIC_MATCHED_ND_WIDTH_PORTFOLIO=1` to retain only the wide
candidate.

On `pre2`, the default wide candidate has a `91.38M` symbolic score and
produces `75.45M` factor entries with `112.88B` factor flops.  The measured
two-leaf candidate scores `90.23M` and produces `72.81M` entries with
`101.29B` flops.  The final 20-update, verify-every-generation
rank-preserving run reports a `9.005s` initial factor, `0.991s` average
refactor, `0.0618s` solve, and `111.09s` modeled H100; its worst relative-L2
residual is `2.02e-16`.  The same-host saved CKTSO paper run is `296.52s`
H100, while SubtreeLU reports `251.32s` but fails the residual limit.  On the
more disruptive entrywise protocol, KLS is `111.87s` versus a fresh CKTSO
`122.54s`, with a `7.23e-13` worst checked residual.  Disabling the width
portfolio restores the `75.45M`-entry factor; disabling the native certificate
does not complete within 55 seconds because it rebuilds the very large scaled
candidate.  The retained rules therefore close the last valid paper-suite
performance gap while preserving the generic portfolio contract.

The broad recurring AMMF/AMF3 symbolic portfolio is now bounded by its actual
immutable ordering-pattern footprint rather than separate row and nonzero
ceilings. AUTO may overlap the three quotient-graph variants when that
footprint is at most 128 MiB, the caller supplies at least four workers, and
the existing column-pair-work/lifecycle test proves the analyses can repay.
This admits large sparse graphs without identifying their dimensions or
provenance and keeps denser graphs under the same memory budget.

`Hamrle3` is the new large positive: its 53.1 MiB pattern had previously been
excluded from the broad portfolio, leaving AUTO on a normal-AMF estimate of
141.1M factor entries and 500B operations and timing out beyond 300 seconds.
The resource-bounded portfolio selects the transpose side from independently
realized lifecycle scores, starts from a 126.1M-entry/358B-operation AMMF
symbolic, and ultimately installs a numerically validated METIS factor. A
strict 20-update rank-preserving run now completes within the same 300-second
process cap: 41.82 s analysis, 37.79 s initial factor, 5.996 s steady
refactor, 0.207 s changed solve, and a 775.10 s modeled H100. Every generation
was checked; the worst relative-L2 residual is `4.32e-10`. CKTSO and SubtreeLU
both still time out on this input. The four other paper matrices newly exposed
to the resource form (`TSOPF_RS_b2383`, its `_c1` variant, `G3_circuit`, and
`memchip`) preserve their prior METIS orientation and remain within 0.9% of
their preceding modeled cycles, with all residuals below `1.17e-12`.

## Generic matched-lifecycle portfolio completion

The remaining exact-update portfolio now settles from realized symbolic,
numeric, and recurring costs rather than matrix dimensions or provenance.  A
value-matched AMF3 factor may admit a lower-power ordering when its measured
elimination-tree span falls by at least 40%, while its estimated score and
work stay inside bounded control/incumbent caps.  The resulting real KLU
numeric must still pass the ordinary fill, work, pivot, conditioning, and
lifecycle gates.  An accepted exact match may retain native KLU scaling only
when it has no `rcond` or off-diagonal-pivot loss and independently clears the
same realized recurring-value test.  Pending representation trials now
settle before a full tight-pivot numeric is evaluated, so the lifecycle does
not optimize a factor that the next representation candidate discards.

The recurring executor portfolio follows the same rule.  A pair-fusion trial
more than 25% slower than the measured unfused EGraph executor cancels the
wider quad trial.  Generic repeated workloads can issue the batch-floor probe
after three warm samples; a close result receives one adjacent re-audit, and
a measured low-floor win settles the overlapping padded-panel choice.  The
PTS refactor challenger is adopted immediately only on its normal 5% win.  A
first result within 5% of parity receives one paired re-audit, which settles
from the two warm minima; a clearly slower challenger is rejected without a
second representation round-trip.  Once PTS has won, reproduced the numeric,
and passed a `1e-9` relative-L2 line in an unscaled full-precision frame, the
ordinary solve contract may retire its componentwise probe.  Unsettled
EGraph states retain the existing armed recovery contract.

The individual A/B controls are
`KLS_DISABLE_MATCHED_DECISIVE_SPAN_REFINEMENT`,
`KLS_DISABLE_GENERIC_NATIVE_MATCHING_SCALE`,
`KLS_DISABLE_GENERIC_TIGHT_PIVOT_DEFERRAL`,
`KLS_DISABLE_EGRAPH_DOMINATED_QUAD_SKIP`,
`KLS_DISABLE_GENERIC_EARLY_BATCH_FLOOR_PROBE`,
`KLS_DISABLE_BATCH_FLOOR_CLOSE_REAUDIT`,
`KLS_DISABLE_PTS_CLOSE_REAUDIT`, and
`KLS_DISABLE_SETTLED_PTS_RAW_L2_CERTIFICATE`.

On the final eight-core, 20-entrywise-update audit, two counterbalanced passes
over 18 development and unrelated holdout matrices completed all 72 checked
runs below the `1e-8` residual limit.  The combined current/control H100
geometric mean was `0.9142`: eight matrices improved by more than 2%, and none
regressed by more than 2%.  `ckt11752_dc_1` measured `0.485` of the control
cycle, while `ckt11752_tr_0` measured `0.997`.  A separate stability audit put
40/40 checked `ckt11752_dc_1` runs on the selected mapped path at a `0.3701 s`
median with a `3.39e-14` maximum generation residual; 20/20
`ckt11752_tr_0` runs retained EGraph at a `0.4091 s` median.

## Value-certified BTF lifecycle selection

AUTO repeated-update analysis can now retain a bounded BTF symbolic when the
ordinary structural tournament selects its no-BTF counterpart.  This does not
classify the input by name, dimensions, or a matrix-family signature.  A
fixed-cost sample of elimination positions must first show decisive static
pivot relief.  The complete BTF numeric must then beat the selected symbolic's
fill and work estimates, have a finite positive diagonal condition estimate,
and pass a deterministic backward-residual solve before it can be installed.
The independently analyzed AMF portfolio is preserved: when the existing AMD
race already speculates on no-BTF, AMF starts on that same representation and
reruns only if BTF wins.

An accepted unscaled factor records its initial condition regime.  Every
changed numeric receives a vectorized diagonal-ratio check; an eightfold
deterioration rebuilds that generation with the scaled AUTO mode that
authorized the trial.  If AUTO had already selected an unscaled incumbent,
the alternative introduces no new scaling risk.  The controls
`KLS_DISABLE_GENERIC_BTF_VALUE_SELECTION=1` and
`KLS_DISABLE_GENERIC_UNSCALED_RCOND_GUARD=1` restore the previous selector and
disable the changed-value guard, respectively.

On the final eight-core, 20-entrywise-update gap suite, the feature/control
geometric-mean ratio over 27 matrices was `0.9866`.  Fifteen-pass medians put
`rajat26` at `0.1420 s` versus `0.2773 s` with the selector disabled and
`rajat23` at `0.3595 s` versus `0.3620 s`; both are below the saved CKTSO
medians (`0.1433 s` and `0.3884 s`).  The same-session 27-matrix KLS geometric
mean was `0.1976 s`, versus `0.2096 s` for CKTSO and `0.5731 s` for SubtreeLU.
KLS won all 27 SubtreeLU comparisons and 15 of 27 CKTSO comparisons, so this is
an aggregate and targeted gap closure rather than an all-matrix dominance
claim.

The 93-matrix paper scan adopted only `rajat23` and `rajat26`; `bips98_606`
was sampled and rejected, with a 15-pass selector/control ratio of `1.0024`.
The 48 group-disjoint development/holdout scan found no adoption and one
value-profile rejection.  One hundred checked entrywise and localized updates
on both positives stayed below `2.13e-13` relative L2 residual at ordinary
amplitude.  At 90% entrywise perturbation, `rajat26` crossed its recorded
condition floor, rebuilt scaled, and stayed below `1.06e-11`; ASan/UBSan runs
exercised both the retained and recovery paths without findings.

## License

KLS is licensed under LGPL-2.1-or-later. The current in-tree solver engine
includes SuiteSparse-derived KLU, AMD, COLAMD, BTF, and UFconfig sources from
Trilinos, plus METIS/GKlib, SCOTCH, and SPRAL scaling support; see
`THIRD_PARTY_NOTICES.md` for attribution.

For a paper-by-paper implementation checklist, see
`docs/paper_ideas_audit.md`.
