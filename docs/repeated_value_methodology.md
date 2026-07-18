# Repeated-value methodology in the CKTSO and SubtreeLU papers

The papers distinguish the intended SPICE workload from their standalone
SuiteSparse timing loops.

## Intended SPICE workload

Both papers state that the coefficient matrix normally keeps its sparsity
pattern while its numerical values change between Newton--Raphson iterations.
CKTSO relies on those changes usually being smooth: it reuses the previous LU
structure and pivot order, checks pivots, and falls back to pivoting when the
prediction fails. SubtreeLU likewise describes refactorization as appropriate
when values vary little and factorization with pivoting as necessary when
stability is at risk.

CKTSO's real nonlinear operating-point simulations therefore exercise changing
values. Its separate linear transient experiment is a different workload: one
factorization followed by 1,000 triangular solves.

## Published SuiteSparse microbenchmarks

The standalone matrix experiments do not describe a numerical perturbation
between timed repetitions:

- CKTSO's benchmark section uses fixed SuiteSparse matrices. Its best-case
  factorization assumes no repivoting. Its synthetic general case repeats each
  matrix 100 times and injects a repivot event in a chosen percentage of runs;
  it does not say that the matrix values are changed.
- SubtreeLU says each solver is run 100 times per test case and reports the
  average factorization/refactorization time, without describing value changes.
- The released CKTSO demo calls `Factorize(ax, true)` and `Refactorize(ax)` 100
  times with the same `ax` array. The released SubtreeLU demo likewise calls
  `factorize(&ax[0])` or `refactorize(&ax[0])` repeatedly without modifying
  `ax`.

Thus byte-identical reuse matches the released standalone microbenchmark loops,
but it is not a sufficient proxy for nonlinear SPICE iterations. KLS should
report two repeated-workload checks:

1. an exact-value result for comparison with the published/released matrix
   microbenchmarks; and
2. a deterministic changed-value result to validate real numeric refactoring.

One-shot results must be collected in a separate process with zero preliminary
factor/refactor repetitions. Otherwise the measured solve may use caches or
layouts prepared by a preceding refactor and is not a cold H1 measurement.

## Cross-solver reporting protocol

The primary comparison with the previous works uses their released-demo shape:

1. analyze and perform one initial factorization;
2. call `Factorize` 100 times with the unchanged values array;
3. call `Refactorize` 100 times with that same array; and
4. report the average repeated-factor and repeated-refactor times separately.

There is no sleep, cache flush, or artificial delay between calls. A combined
100-iteration SPICE horizon is a useful KLS deployment metric, but it is not a
replacement for the two separate timing columns used by the previous works.
Use `scripts/score_repeated_factorization.py` with
`--min-factor-repeat 100 --min-refactor-repeat 100` to enforce this protocol;
the scorer rejects zero-repeat records so it cannot silently compare a cold
initial factor against a warm repeated-factor average.

The deterministic changed-value check remains secondary. It tests whether an
optimization that benefits the byte-identical released loops still performs
real numeric work correctly, rather than redefining the published comparison.

Run it with:

```sh
REFACTOR_VALUES=rank-preserving \
REFACTOR_VALUE_AMPLITUDE=0.001 \
  scripts/run_paired_suite.sh ...
```

At refactor generation `g`, every harness constructs the identical matrix

`A_g = D_row(g) A_0 D_col(g)`.

The diagonal multipliers are deterministic coordinate/generation hashes in
`[1-a, 1+a]`, where `a` is the amplitude. Because the command-line interface
requires `0 < a < 1`, both diagonal matrices are nonsingular and `A_g` has
exactly the same rank as `A_0`. The sequence is independent of matrix names,
values, dimensions, and the corpus. Each harness rebuilds the known-solution
right-hand side outside the solver timer before refactorization. This avoids
giving any solver credit for benchmark value generation while checking the
final residual against the actual final-generation matrix.

The default remains `REFACTOR_VALUES=unchanged`, preserving direct parity with
the released CKTSO and SubtreeLU demo loops. JSON output records both
`refactor_value_mode` and the effective amplitude so changed- and unchanged-
value records cannot be confused during review.

For a deterministic rejection workload, use
`REFACTOR_VALUES=entrywise`.  It applies an independent coordinate/generation
multiplier to every stored entry, so it is generally not expressible as
`D_row A_0 D_col`.  This mode exercises certificate rejection and ordinary
numeric fallback.  Unlike the rank-preserving mode, it does not guarantee
rank preservation; use a small amplitude and retain residual/status filters.

For a solver-neutral localized counter-workload, use
`REFACTOR_VALUES=localized-entrywise`. It uses the same independent
entrywise multiplier, but only within one fixed cyclic window containing
one in 1024 columns (at least one column). The window is defined solely by
the matrix dimension and a published constant, not by a solver ordering or
BTF decomposition, and is identical in every harness. This mode is intended
to measure localized numeric-update mechanisms; it has the same rank caveat
as the full entrywise mode.

## Opt-in partial-BTF retention

Set `KLS_ENABLE_PARTIAL_BTF_REFACTOR=1` before the reference factorization to
enable the localized numeric path. For an unscaled, multi-block BTF numeric,
KLS retains an exact reference value array and maps every stored entry either
to an independent diagonal block or to KLU's off-diagonal `Offx` storage. A
refactor compares every entry, copies changed couplings, and numerically
refreshes only changed diagonal blocks. It falls through to the ordinary full
refactor on any unsupported numeric state or failed refresh; scaled,
predicted, perturbed, reduced-precision, and refinement-dependent numerics are
deliberately ineligible.

The default work gate accepts a partial update only when changed blocks hold
at most half of the stored factor-entry work proxy. Override that fraction
with `KLS_PARTIAL_BTF_MAX_WORK_FRACTION` in `(0,1)`. Two consecutive
over-budget observations suppress later scans until an explicit
`kls_factor`; `KLS_DISABLE_PARTIAL_BTF_REJECTION_GATE=1` disables that
experimental backoff. `KLS_TRACE_PARTIAL_BTF=1` reports the changed-entry,
block, and work counts.

As a bounded paper-union check, `IBM_EDA/ckt11752_tr_0` has 199 BTF blocks and
the localized generator changes 196 entries in one block, representing
194,768 of 1,029,808 factor-work units. Seven solver-order-rotated processes
with eight cores, 20 refactors, 100 measured solves, and amplitude 0.001 give
median `H100` values of 0.2210 s for partial-BTF KLS, 0.2856 s for CKTSO, and
0.4161 s for SubtreeLU; full-refactor KLS takes 0.4092 s. All relative
residuals are at most 2.1e-15. This is a mechanism check on one localized
case, not a union-wide arbitrary-update ranking.

## Local primary material

- `../../refs/CKTSO_High-Performance_Parallel_Sparse_Linear_Solver_for_General_Circuit_Simulations.pdf`, Sections I, IV, and VI.
- `../../refs/SubtreeLU_High-Performance_Parallel_Sparse_LU_Factorization_for_Circuit_Simulation.pdf`, Sections I, III, and IV.
- `../../cktso/demo/benchmark.cpp`.
- `../../SubtreeLU/demo/demo.cpp`.
