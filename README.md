# KLS

KLS is a standalone sparse direct solver project targeting SPICE-style
workloads: one symbolic analysis followed by many numeric factorizations,
refactorizations, and solves with a fixed sparse structure.

This repository currently contains the first working KLS implementation:

- A stable C API in `include/kls/kls.h`
- A vendored SuiteSparse-derived 64-bit symbolic/numeric engine
- CSC and CSR input paths with 32-bit or 64-bit index arrays
- AMD-first automatic symbolic ordering with explicit AMD, COLAMD, natural,
  METIS, and SCOTCH nested-dissection ordering controls
- SPICE-cycle-oriented normal-vs-transpose internal orientation selection, with
  explicit orientation controls
- Factor, refactor, solve, transpose-solve, and statistics APIs
- Fast repeated factorization that reuses the existing numeric pattern, checks
  pivot quality, and can repair unscaled rejected BTF blocks before falling
  back to full pivoting factorization
- Value-aware static row pivoting and row/column equilibration trials for high
  off-diagonal-pivot cases
- A MatrixMarket benchmark tool
- A small correctness smoke test
- A SuiteSparse Matrix Collection downloader script for public benchmark cases

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

For very large pre-static candidates, KLS uses SPRAL's auction
matching/scaling path instead of exact Hungarian matching. This keeps the
MC64-adjacent step on LGPL-compatible redistributed code while avoiding an
unbounded exact-assignment setup cost on large circuit matrices.

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
decisions without running numeric factorization.

When system SuiteSparse KLU headers and libraries are installed, the build also
provides `klu_width_compare` to compare system `klu_*` and `klu_l_*` on the same
MatrixMarket input. This is a diagnostic benchmark for deciding whether a
future 32-bit KLS backend is worth implementing:

```sh
./build/klu_width_compare matrix.mtx --repeat 3 --refactor-repeat 3 --json
```

JSON includes `initial_factor_path` and `last_factor_path`; values such as
`klu_first` or `klu_fallback` mean the factorization was handed to the
KLU-derived pivoting kernel, while `kls_fast_refactor` means KLS reused the
retained pattern through the checked fast path. Repeated refactor diagnostics
are reported separately as `last_refactor_path`, distinguishing row-refactor,
EGraph, mapped, pool, and serial KLU refactor branches. Unset
`KLS_ENABLE_KLS_FIRST_FACTOR` keeps the production cold first factor on the
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
`kls_first_last_row_pipeline_prefix_panel_rebuild_rows`; separator pipeline
epochs continue to use the separator-prefixed counters below. KLS-first
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
`kls_first_separator_dynamic_column_fallback_count`. If a guarded separator
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
counters. Prefix panel rebuilds after a separator-pipeline pivot are reported
through `kls_first_last_separator_queue_pipeline_prefix_panel_rebuild`,
`kls_first_separator_queue_pipeline_prefix_panel_rebuild_count`, and
`kls_first_last_separator_queue_pipeline_prefix_panel_rebuild_rows`, matching
the retained prefix semantics of the epoch retry. Active-rank pivot resets that
still rebuild row-supernode metadata are exposed separately through
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
over-stage these rows. These cover contiguous independent-row producer suffixes
and fragmented dense-consumer external prefixes without requiring CBLAS. The
independent-row multi-producer executor now keeps its metadata, pivots, U
scratch, and panel descriptors in per-worker scratch instead of allocating them
for each accepted batch. Dense-consumer fragmented producer batches use the same
dedicated object workspace for their run-panel descriptors while retaining a
separate byte workspace for fallback target maps. This is a storage cleanup for
the paper-shaped executor, not a default policy change. Current top-five
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
That executor now addresses its scratch data through the retained
group/current workspace map, so the opt-in numeric path exercises the same
compact current-slot layout that a true multi-current batch kernel will need.
This is also intentionally experimental and off by default; it tests whether
the Algorithm 5 payoff surface is enough without a true multi-current batch
executor.
When `KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_PLAN=1` or the stricter
exec flag is set, KLS now also retains the selected producer groups as a
multi-current descriptor for that future executor. The retained internal layout
maps each selected run to its current-column slot and gives that slot a compact
workspace offset, so a batch kernel can gather/advance/publish current
workspaces without rediscovering the group shape. Benchmark JSON reports
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
`refactor_supernode_algorithm5_payoff_group_max_pattern_width`. These counters
estimate the number of distinct current-column workspaces, compact prefix
workspace rows, addressable target accumulator slots, retained update-entry
work, and prefix-advance dependencies a real Algorithm 5 batch would need. On
the slow ASIC diagnostics, zero-advance selected runs are rare or absent, so the
next paper-aligned executor has to batch prefix advancement before the shared
producer update rather than relying on a first-dependency shortcut or BLAS
threshold tuning. These counters do not enable the old per-current scalar replay
by default.
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
A separate CKTSO Algorithm 5-style probe tried letting the column EGraph
refactor consume later already-finished scalar dependencies while waiting for an
earlier unfinished dependency. The structural safety scan was correct but not a
usable default: KLU-first focused runs regressed `ASIC_320ks`, `ASIC_100ks`, and
`G2_circuit`, and `onetone2` timed out at the 120 s harness limit. The probe was
removed rather than kept behind another runtime flag; the same Algorithm 5 idea
remains implemented in the row-major row-refactor path where dependency metadata
is already row-oriented.
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
remaining subtrees are assigned to private thread queues by retained group
work, and separator-crossing row groups are forced into the pipeline queue.
Indivisible retained components are left in the private candidate set, so an
unbalanced separator tree cannot collapse the Algorithm 6 queue into all
pipeline work with no private subtrees.
If a private subtree group depends on an already-pipeline group, KLS promotes
that dependent group into the pipeline closure instead of discarding the
separator schedule.
Stats report this path through
`row_refactor_last_separator_flop_queue`,
`row_refactor_separator_flop_queue_run_count`,
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

For SubtreeLU-specific tuning, `bench/suitesparse_subtreelu_manifest.txt`
contains the exact 46 public SuiteSparse circuit labels from SubtreeLU Fig. 7.
Fetch it into its own directory when running the full set, so the suite runners
see one matrix per basename:

```sh
python3 scripts/fetch_suitesparse.py --manifest bench/suitesparse_subtreelu_manifest.txt --out data/suitesparse-subtreelu
python3 scripts/audit_paper_manifest.py --manifest bench/suitesparse_subtreelu_manifest.txt --source SubtreeLU
```

Use the medium paper suite as the normal inner loop for solver-policy changes.
Use the large supplement as an overnight or pre-merge generalization gate with
explicit timeouts; a retained change should improve the paper geomean or a
defensible structural class, not just one large matrix name:

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
already shown KLS/CKTSO differences under a 120s per-process cap:

```sh
python3 scripts/run_bench_suite.py --kls-bench build/kls_bench --matrix-dir data/suitesparse-paper-large --manifest bench/suitesparse_paper_large_recon_manifest.txt --orientation auto --threads 4 --repeat 1 --refactor-repeat 1 --timeout 120 --jsonl build/kls_paper_large_recon.jsonl
python3 scripts/run_cktso_suite.py --cktso-compare build-cktso/cktso_compare --matrix-dir data/suitesparse-paper-large --manifest bench/suitesparse_paper_large_recon_manifest.txt --threads 4 --repeat 1 --refactor-repeat 1 --timeout 120 --jsonl build/cktso_paper_large_recon.jsonl
python3 scripts/run_klu2_suite.py --klu2-compare build-klu2/klu2_compare --matrix-dir data/suitesparse-paper-large --manifest bench/suitesparse_paper_large_recon_manifest.txt --repeat 1 --refactor-repeat 1 --timeout 120 --jsonl build/klu2_paper_large_recon.jsonl
```

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
runs. Add `--include-failures --failure-seconds SECONDS` when a hard-suite
comparison should score failed or missing rows with an explicit cycle-time
penalty; using a process timeout value as the penalty is only a lower bound for
SPICE-cycle comparisons.

The suite metric is:

```text
analysis + initial_factor + solve + 99 * (refactor + solve)
```

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
rows, and local/shared row-entry reserve growth/copy volume. Same-session `pre2`
probes show the default AMD run enters the 629,628-row dominant BTF block with
no separator coverage, while forced METIS enters the same block with separator
coverage but still exceeds the 120s cap. The matching local CKTSO run finishes
analysis, first factor, one factor, one refactor, and solve in about 21s wall
time, so the remaining gap is not the timeout limit or BLAS thresholding; it is
the missing CKTSO/SubtreeLU coarse row/supernode numeric executor inside that
dominant block.
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
`panel_updates`, `panel_update_rows`, `panel_appends`, and
`panel_append_entries`. This is useful but not enough for `pre2`: a traced
forced-METIS run still reached the pivot tail at row 274,430 with about
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
current refactor bridge, but the next larger algorithmic work is still to evolve
the numeric factor/refactor/solve kernels toward those deeper KLS-owned sparse
kernels while keeping the public API and benchmark harness stable.

## License

KLS is licensed under LGPL-2.1-or-later. The current in-tree solver engine
includes SuiteSparse-derived KLU, AMD, COLAMD, BTF, and UFconfig sources from
Trilinos, plus METIS/GKlib, SCOTCH, and SPRAL scaling support; see
`THIRD_PARTY_NOTICES.md` for attribution.

For a paper-by-paper implementation checklist, see
`docs/paper_ideas_audit.md`.
