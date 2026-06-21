# KLS

KLS is a standalone sparse direct solver project targeting SPICE-style
workloads: one symbolic analysis followed by many numeric factorizations,
refactorizations, and solves with a fixed sparse structure.

This repository currently contains the first working KLS implementation:

- A stable C API in `include/kls/kls.h`
- A vendored SuiteSparse-derived 64-bit symbolic/numeric engine
- CSC and CSR input paths with 32-bit or 64-bit index arrays
- AMD-first automatic symbolic ordering with explicit AMD, COLAMD, natural, and
  METIS nested-dissection ordering controls
- SPICE-cycle-oriented normal-vs-transpose internal orientation selection, with
  explicit orientation controls
- Factor, refactor, solve, transpose-solve, and statistics APIs
- Fast repeated factorization that reuses the existing numeric pattern before
  falling back to full pivoting factorization
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
pinned submodules under `third_party/metis` and `third_party/gklib`; initialize
them with `git submodule update --init --recursive` after cloning. To use a
compatible system METIS instead, configure with:

```sh
cmake -S . -B build -DKLS_USE_SYSTEM_METIS=ON
```

KLS builds METIS with 64-bit `idx_t` for compatibility with the 64-bit KLS/KLU
path. A system METIS install used this way must be ABI-compatible.

## Benchmark

```sh
./build/kls_bench ../cktso/demo/add20.mtx --repeat 20 --refactor-repeat 20 --orientation auto --json
```

The benchmark reports analysis, factorization, refactorization, solve,
transpose-solve, residual, selected orientation, BTF block/rank, fill, flop, and
memory statistics.

Use `--orientation auto|normal|transpose` to control KLS's internal storage
orientation. `auto` uses the transposed pattern directly for small and medium
matrices where avoiding a second symbolic analysis is usually faster in a
repeated SPICE solve cycle. Larger matrices still compare normal and transposed
symbolic fill estimates before choosing.

Use `--scale auto|-1|0|1|2` to compare KLS/KLU row scaling modes when studying
pivoting and refactorization behavior. KLS defaults to `auto`, which can start
unscaled for mostly diagonal circuit matrices with moderate diagonal magnitude
spread, or sum-scaled when the diagonal is sparse but row magnitudes are already
balanced. Otherwise it starts from KLU's max row scaling mode (`2`). For large
expensive cases, `auto` can still try other numeric scaling modes when actual
flop/fill evidence justifies the extra work. Explicit numeric scale values
remain fixed.

Use `--no-btf` to measure the same ordering/scaling policy without KLU's BTF
decomposition. With `--ordering auto` and BTF enabled, KLS can still bypass BTF
after analysis for large matrices where BTF finds a single block and a no-BTF
symbolic retry keeps the fill estimate within tolerance; benchmark JSON reports
both `requested_btf` and selected `btf`.

Use `--ordering metis` to force METIS nested-dissection ordering. The default
`--ordering auto` first compares AMD and COLAMD by symbolic fill estimate. When
METIS is enabled, `auto` can then promote large, expensive first numeric
factorizations to METIS if the trial factorization materially reduces actual
numeric flop/fill cost. This keeps METIS available for hard nested-dissection
cases without paying its analysis cost on small circuit matrices.

To fetch public SuiteSparse Matrix Collection matrices listed in the manifest:

```sh
python3 scripts/fetch_suitesparse.py --manifest bench/suitesparse_circuit_manifest.txt --out data/suitesparse
```

To run every downloaded matrix and compute the SPICE-cycle geometric mean:

```sh
python3 scripts/run_bench_suite.py --kls-bench build/kls_bench --matrix-dir data/suitesparse --orientation auto --jsonl build/kls_suite.jsonl
```

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

## Status

This is a functional serial implementation with KLS-level analysis choices for
repeated SPICE-style solves. It is not yet a generally CKTSO-beating solver
across broad circuit corpora. The next algorithmic work is to evolve the
numeric factor/refactor/solve kernels toward KLS-owned sparse kernels with
better pivot reuse and parallelism while keeping the public API and benchmark
harness stable.

The `threads` option is accepted for API stability but is not used by the
current serial KLS numeric engine.

## License

KLS is licensed under LGPL-2.1-or-later. The current in-tree solver engine
includes SuiteSparse-derived KLU, AMD, COLAMD, BTF, and UFconfig sources from
Trilinos, plus optional METIS/GKlib ordering support; see
`THIRD_PARTY_NOTICES.md` for attribution.
