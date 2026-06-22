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

The main build is self-contained. KLS vendors the SuiteSparse-derived KLU,
AMD, COLAMD, and BTF C sources from Trilinos under `third_party/suitesparse`
as the current in-tree serial engine. METIS ordering is enabled by default from
pinned submodules under `third_party/metis` and `third_party/gklib`, and SCOTCH
ordering is enabled by default from `third_party/scotch`; initialize them with
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

KLS can optionally link to a compatible external SPRAL install for BSD-licensed
auction matching support:

```sh
cmake -S . -B build -DKLS_ENABLE_SPRAL_SCALING=ON \
  -DKLS_SYSTEM_SPRAL_INCLUDE_DIR=/path/to/spral/include \
  -DKLS_SYSTEM_SPRAL_LIBRARY=/path/to/libspral.so
```

This is not a solver backend and is not enabled by default. KLS uses it only
as a large-matrix fallback when the in-tree static matcher is far below the
normal cardinality acceptance threshold.

## Benchmark

```sh
./build/kls_bench ../cktso/demo/add20.mtx --repeat 20 --refactor-repeat 20 --orientation auto --json
```

The benchmark reports analysis, factorization, refactorization, solve,
transpose-solve, residual, selected orientation, BTF block/rank, fill, flop, and
memory statistics. Use `--analyze-only` to measure symbolic analysis and
ordering decisions without running numeric factorization.

Use `--threads N` to enable KLS-owned parallel work where it is currently
available. The first threaded path is repeated numeric refactorization across
independent BTF diagonal blocks for large, high-flop cases, including existing
KLU row-scaling modes when reusable scale storage and enough independent BTF
block work are present. KLS keeps a solver-owned worker pool for this path, so
repeated SPICE refactors do not relaunch threads every cycle. Small matrices,
factor kernels, solve kernels, and BTF cases that are not wide enough for the
threaded path still use serial execution. For unscaled serial refactors, KLS
also precomputes a scatter map from the fixed pivot order so repeated
single-block and serial BTF refactors do not redo the same `Q` and `Pinv`
structure lookups every cycle. For serial BTF refactors, this map also keeps
each column's off-block entries separate from diagonal-block entries so the
numeric refactor loop avoids reclassifying the same structure on every SPICE
step.

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
Large low-degree, nearly diagonal circuit patterns and METIS-started medium
bounded-degree or dense-diagonal patterns also start in no-scale mode to avoid
repeated row-scale work once KLS has validated the structure. TSOPF-style
spiked low-diagonal METIS starts use max scaling, and large sparse-diagonal
low-degree patterns can start with sum scaling when that avoids more expensive
max-scaling behavior. For large expensive
cases, `auto` can still try other numeric scaling modes when actual flop/fill
evidence justifies the extra work, but structural METIS starts that are already
known to need a specific scale mode skip redundant scale trials. Explicit
numeric scale values remain fixed.

Use `--pivot-tol T` to benchmark the diagonal pivot tolerance exposed by the
KLS API. The default is `0.001`, matching the underlying KLU default. With that
default, KLS can trial a lower `1e-4` tolerance on large high-fill matrices with
noticeable off-diagonal pivoting, keeping it only when fill and pivoting improve
without a large reciprocal-condition drop. Benchmark JSON reports both requested
and selected pivot tolerance.

When `fast_factor` reuses an existing factor pattern, KLS also checks the
resulting L multipliers against the selected pivot tolerance. On unscaled BTF
patterns, a rejected reused pivot can restart the rejected diagonal block with
pivoting, splice the repaired block permutation back into the numeric object,
and retry the checked fast factorization. Scaled fast-factor calls still use
the pivot-checking refactor kernel so they can stop at the first unsafe
multiplier, but they fall back to full pivoting factorization when a repair is
needed. Benchmark stats report `fast_rejected_pivot`,
`fast_rejected_pivot_col`, and `fast_block_restarts` for the first rejected
factor-order pivot, its original matrix column, and the number of repaired BTF
blocks. These fields are intended to guide fuller CKTSO-style tail-restart
work without accepting an unsafe reused pivot order.

Use `--no-static-pivoting` to disable KLS's value-aware static row-pivoting
trial. When enabled, `auto` can preemptively build a weighted row permutation
for medium matrices whose input values show a mostly weak or missing diagonal,
or react after a high off-diagonal pivot count. The permutation moves large
entries onto the diagonal. For larger weak-diagonal candidates, the matching
augment uses a layered bipartite search so KLS can complete many independent
augmenting paths per pass instead of restarting a search from each unmatched
row. For cheap small candidates, KLS can instead run an exact sparse
augmenting-path assignment on transformed log magnitudes, equivalent to a
maximum-product diagonal match when a full matching is found. The exact path is
gated by matrix order and `n * nnz` work because even the KLS-owned
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
setup when it would increase fill. Benchmark JSON reports both whether static
pivoting was enabled, whether KLS selected it, and whether the accepted static
match used exact assignment. KLS does not vendor HSL MC64 or the MC64 copies
carried by some solver projects. Optional SPRAL auction matching can be linked
from a compatible external install; otherwise MC64-equivalent code must be
LGPL-compatible in-tree KLS code.

Use `--no-btf` to measure the same ordering/scaling policy without KLU's BTF
decomposition. With `--ordering auto` and BTF enabled, KLS can still bypass BTF
for large matrices where BTF finds either a single block or one dominant block
with a small fringe, and a no-BTF symbolic retry gives enough evidence to keep
or improve the fill estimate. For large, low-degree, nearly
diagonal patterns, medium low-degree full-diagonal patterns, medium sparse
high-degree mostly diagonal patterns, and for METIS-started medium
dense-diagonal high-degree patterns, `auto` can also start without BTF to avoid
analysis work it is likely to discard; benchmark JSON reports both
`requested_btf` and selected `btf`.

Use `--ordering metis` to force METIS nested-dissection ordering. The default
`--ordering auto` can start with METIS for medium, bounded-degree structures
that resemble mesh-style nested-dissection cases and for medium dense-diagonal
high-degree patterns where delayed promotion would otherwise pay for an
avoidable first AMD factorization. It can also start with METIS for medium
spiked low-diagonal patterns that resemble TSOPF-style paper cases where the
symbolic estimate understates nested-dissection benefit. Otherwise it compares
AMD and COLAMD by symbolic fill estimate. Dense-diagonal high-degree METIS
starts ask METIS for two separator attempts, then refine the nested-dissection
rank order with CAMD inside coarse rank constraints. This preserves the
separator-first shape while letting minimum degree reduce local fill/flops on
post-layout style circuits. Large METIS orderings also use coarse rank-group
CAMD refinement, matching the CKTSO paper's nested-dissection plus constrained
minimum-degree structure more closely than raw METIS. Medium bounded-degree
METIS starts also use two separator attempts with random matching coarsening to
reduce nested-dissection factor work on mesh-like sparse diagonals. When METIS
is enabled, `auto` can also promote
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
matrix with no empty rows or columns, or a near-full diagonal with a large dense
row/column spike. These predicates are structural, not matrix-name-based, and
are deliberately narrow so unrelated IBM `dc`/`trans` cases stay on the cheaper
AMD/COLAMD path.

Use `--ordering scotch` to force SCOTCH nested-dissection ordering. SCOTCH is
kept as an explicit experimental backend rather than part of `auto` until it
shows a general win over the current AMD/METIS policy. ParMETIS is not wired
into KLS yet; it is a distributed-memory MPI package and should be treated as a
separate future backend if KLS grows an MPI/distributed solver path.

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
manifest includes very large matrices, so the medium manifest is the default
practical corpus:

```sh
python3 scripts/fetch_suitesparse.py --manifest bench/suitesparse_paper_medium_manifest.txt --out data/suitesparse-paper-medium
```

For deeper tuning against the remaining paper matrices, use the large
supplement. It contains the 17 public SuiteSparse paper cases excluded by the
medium caps, including the largest KLU/NICSLU/SubtreeLU labels:

```sh
python3 scripts/fetch_suitesparse.py --manifest bench/suitesparse_paper_large_manifest.txt --out data/suitesparse-paper-large
```

Use the medium paper suite as the normal inner loop for solver-policy changes.
Use the large supplement as an overnight or pre-merge generalization gate with
explicit timeouts; a retained change should improve the paper geomean or a
defensible structural class, not just one large matrix name:

```sh
python3 scripts/run_bench_suite.py --kls-bench build/kls_bench --matrix-dir data/suitesparse-paper-large --orientation auto --threads 4 --timeout 900 --jsonl build/kls_paper_large.jsonl
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

The suite metric is:

```text
analysis + initial_factor + solve + 99 * (refactor + solve)
```

Compare two JSONL runs by matrix basename:

```sh
python3 scripts/compare_bench_runs.py --candidate build/kls_suite.jsonl --candidate-name kls-auto --reference build/klu_defaults.jsonl --reference-name klu-defaults
```

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
refactors, plus an unscaled block-local pivot restart for fast-factor failures.
It is not yet a generally CKTSO-beating solver across broad circuit corpora.
The next algorithmic work is to evolve the numeric factor/refactor/solve
kernels toward deeper KLS-owned sparse kernels with better pivot reuse and
parallelism while keeping the public API and benchmark harness stable.

## License

KLS is licensed under LGPL-2.1-or-later. The current in-tree solver engine
includes SuiteSparse-derived KLU, AMD, COLAMD, BTF, and UFconfig sources from
Trilinos, plus optional METIS/GKlib and SCOTCH ordering support; see
`THIRD_PARTY_NOTICES.md` for attribution.

For a paper-by-paper implementation checklist, see
`docs/paper_ideas_audit.md`.
