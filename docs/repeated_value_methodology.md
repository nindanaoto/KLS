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

## Local primary material

- `../../refs/CKTSO_High-Performance_Parallel_Sparse_Linear_Solver_for_General_Circuit_Simulations.pdf`, Sections I, IV, and VI.
- `../../refs/SubtreeLU_High-Performance_Parallel_Sparse_LU_Factorization_for_Circuit_Simulation.pdf`, Sections I, III, and IV.
- `../../cktso/demo/benchmark.cpp`.
- `../../SubtreeLU/demo/demo.cpp`.
