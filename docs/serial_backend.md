# Lean serial backend

`KLS_BACKEND_SERIAL` is an explicit one-thread policy for fixed-pattern SPICE
workloads. Its goal is to remove KLS policy overhead from the cold serial path
while preserving KLS's faster repeated-factor machinery where measurements and
factor invariants justify it. It does not change `KLS_BACKEND_AUTO` or
`KLS_BACKEND_KLS`, so enabling it cannot silently alter an existing parallel
run.

## Policy

- Exactly one thread is required. A serial request with `threads != 1` returns
  `KLS_ERR_UNSUPPORTED`.
- Automatic ordering resolves directly to AMD. Explicit AMD, COLAMD, and
  natural ordering are supported; the serial policy does not launch METIS,
  SCOTCH, matching, or orientation competitions.
- Automatic orientation resolves to normal for CSC input and transpose for CSR
  input, avoiding an analysis race and an avoidable format conversion.
- The initial numeric factor uses the existing pivoting, 64-bit KLU-derived
  kernel. Automatic scaling starts unscaled, which removes a value pass and
  permits the mapped refactor path. It retries with KLU row-sum scaling after a
  failed/singular factor or when measured factor work exceeds the symbolic
  estimate by more than 50%. The latter is a numerical-fill signal, independent
  of matrix name or corpus membership.
- Refactor metadata is deferred until the first `kls_refactor` call so a
  factor-and-solve workload does not pay for machinery it never uses. Public
  `stats.refactor_seconds` includes that deferred work; serial engine A/B
  decisions use a separate kernel-only clock.
- Most successful unscaled factors use the existing position-mapped refactor.
  A supernodal A/B trial is considered only when the factor value arrays occupy
  at least 32 MiB, the pivoting factor has no off-diagonal pivots, and the
  largest BTF block covers at least 75% of the matrix. After that structural
  gate, SNB is retained only when its timed refactor is at least 10% faster.
  Parallel KLS keeps its existing, stricter 40% acceptance margin.
- If mapped invariants do not hold, the backend falls back to the existing KLU
  refactor. Correctness and singularity behavior remain those of the public KLS
  API.

Always initialize `kls_options` with `kls_default_options`. The initializer
sets the options ABI marker as well as defaults. The marker occupies padding in
the original layout, allowing a library with this backend field to recognize an
older binary's options object without interpreting old tail padding as a
backend request.

## Validation protocol

The optimization split was frozen before policy measurements from the 35-entry
`KLS-paper/evaluation/manifest.txt` corpus:

- Development: every fifth manifest entry, seven matrices.
- Holdout: the remaining 28 matrices, not consulted while selecting the first
  policy.
- Final: all successfully completed matrices from both sets.

The first holdout run exposed two unscaled numerical-fill explosions. The fix
was the factor-work/symbolic-work retry above; it contains no name, exact size,
or manifest-membership checks. The final run used a Release GCC 13.3 build with
CBLAS enabled on an AMD Ryzen 9 9950X3D, one pinned core,
`OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`, one initial factor, five measured
refactors, and the median of three passes. KLU, CKTSO, and SubtreeLU samples
were collected in rotating paired order. A horizon of `N` is one analysis, one
initial factor, `N` solves, and `N-1` refactors; the first refactor is accounted
separately from steady refactors.

On the 34-matrix common completed set, lower time is better and the table shows
`competitor time / serial-KLS time`:

| Horizon | KLU | CKTSO | SubtreeLU |
| ---: | ---: | ---: | ---: |
| 1 | 1.077x | 4.472x | 3.635x |
| 10 | 1.376x | 1.984x | 1.860x |
| 100 | 1.600x | 0.921x | 1.115x |

Thus the opt-in backend closes the reported serial KLU gap on this corpus and
is faster at all three horizons. It also beats SubtreeLU at horizon 100. CKTSO
remains about 8% faster at horizon 100, so this work does not claim to close
that repeated-refactor gap.

The horizon-100 KLU speedup was 2.82x on the seven-matrix development set and
1.38x on the 27 completed holdout matrices. The holdout KLU ratios at horizons
1, 10, and 100 were 0.983x, 1.205x, and 1.382x: the lean backend remained 1.7%
slower for a holdout one-shot, then won once refactor reuse mattered. This is
the more useful generalization check and avoids claiming a universal cold-start
win. Maximum relative residual was `3.62e-12`. Of the 34 final rows, 30 selected
mapped refactor, three selected SNB, and one used KLU refactor; 32 stayed
unscaled and two triggered the structural work-inflation retry to scale 2.

`ss1` is excluded from the common timing set because its initial factorization
exceeded the 60-second process cap in both serial KLS and KLU. CKTSO and
SubtreeLU completed it. This is a recorded coverage limitation, not a scored
win.

## Parallel safety check

The default backend was compared against the exact pre-change `ed53b6e` build
with matched compiler, BLAS, METIS, and SPRAL settings on the seven development
matrices. Three-pass paired whole-process measurements found no stable path or
cost regression: the final four-thread matrix-geomean ratio was 0.995
(`current / baseline`), while the earlier two- and eight-thread checks were
1.015 and 0.992 respectively. The two-thread median ratio was 1.002. A noisy
four-thread `onetone2` outlier was repeated for 21 alternating pairs; its
paired median ratio was 1.002. These measurements support leaving the default
parallel policy unchanged, but they should not be read as a claim that this
serial work improves parallel performance.

## Rejected width shortcut

The excluded `klu_width_compare_vendored` target builds the same vendored KLU
sources in 32-bit and 64-bit modes for an exact diagnostic. The 32-bit route
was not promoted: on equal-fill development cases its scaled refactor was
roughly 30--35% slower, and two cases changed AMD fill materially. Input CSC
indices may still be 32-bit through the public API, but the production numeric
kernel remains the established 64-bit implementation until a broader result
justifies another backend.

## Reproduction

Build and run the serial policy with:

```sh
cmake --build build -j
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 0 \
  python3 scripts/run_bench_suite.py \
    --kls-bench build/kls_bench \
    --matrix-dir data/suitesparse \
    --manifest ../KLS-paper/evaluation/manifest.txt \
    --backend serial --threads 1 --repeat 1 --refactor-repeat 5 \
    --passes 3 --timeout 60 --jsonl build/serial-kls.jsonl
```

The adjacent `.failures` file records capped rows. `kls_bench` JSON exposes the
selected `backend`, `initial_factor_path`, `last_refactor_path`, scale, residual,
and first/steady refactor timings used to audit every policy decision.
