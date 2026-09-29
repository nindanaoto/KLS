# Accumulated-error certificate qualification

Status: single-pass stage qualified on 2026-09-28; worker-budget-aware row
certification subsequently passed both full qualification gates. The existing
componentwise policy remains the default. Numerical kernels and correction
budgets are unchanged. Application-managed opt-out is not an accuracy-equivalent
replacement (MOS13 fails convergence in that mode).

## Qualified single-pass results

All six full accuracy audits passed: BJT had 1,146 certified successes per
thread configuration; MOS13 had 30,490 successes and one explicitly reported,
uniquely matched rejection per configuration; mem_plus had 1,201 successes per
configuration. Every successful answer passed the unchanged exact-rational
componentwise limit. Waveforms matched KLU within the declared tolerance; BJT
and mem_plus were identical. MOS13's maximum absolute difference was about
`8e-10`, with verified full-simulation recovery.

All 336 paired timing simulations (three screening and nine confirmation
blocks) completed with valid waveforms. Seven additional KLU reference runs
are not included in that count. Confirmed median paired improvements in total
simulation time versus the committed interval-only certificate follow:

| Circuit | One thread | Eight threads |
| --- | ---: | ---: |
| BJT | 16.62% | 22.43% |
| MOS13 | 28.51% | 28.41% |
| mem_plus | 26.79% | 26.88% |
| mux8 | 22.05% | 21.89% |
| DAC | 24.55% | 24.42% |
| Rallpack3 | 20.61% | 20.63% |
| ram2k | 24.81% | 24.73% |

Both benefit and regression gates passed. These are single-machine results,
not a claim of superiority to competing solvers. The frozen candidate replaced
only the certificate object; all other Xyce/KLS objects were held fixed.
The 37-test suite, GCC 14/Clang 18 generic/native differential and bound tests,
unsafe-flag rejection tests, and allocation-fault sanitizer checks passed.

Local provenance is under `/tmp/kls-refinement-first/`:
`fast-certificate-v2/combined-accuracy.json`,
`fast-certificate-v2-performance/confirmation-summary.json`, and the corresponding
build/observation manifests. The first mem_plus-8 capture was interrupted by
the storage-reserve guard; its failure record remains intact. A fresh scoped
retry supplied the missing audit, and the continuation verified and combined
the five prior audits with it before starting timings. Raw files reclaimed
during recovery are preserved in hash-verified archives with member receipts.

## Proof and fallback

For a unique-coordinate row with `k` entries, let
`S = |b| + sum |a*x|`, `u = 2^-64`, and `gamma_m = m*u/(1-m*u)`.
Each product and accumulation is rounded in extended precision. Both residual
error and magnitude error are at most `gamma_(2k)*S`. Therefore the computable
bound `gamma_(2k)/(1-gamma_(2k))*Shat` encloses both errors. The implementation
uses the larger exact coefficient `c = 8*(k+1)*u`, restricted to `c <= 1/4`.
This restriction implies that the integer conversion and coefficient evaluation
are exact. These constants establish the proof; they are not CPU tuning knobs.

An outward-rounded `c*Shat` bounds the residual from above and `S` from below.
Acceptance compares the upper residual against the lower denominator times a
lower enclosure of the exact rational tolerance `1/100000000`. A zero magnitude
requires exactly zero residual. The reported maximum ratio is rounded outward.

The fast proof requires IEEE binary64 inputs, binary radix, x87 extended
precision, and nearest rounding with full x87 precision. Products and sums stay
far inside the extended exponent range, including subnormal binary64 inputs.
Signed duplicates bypass this proof: their coefficients must be coalesced before
forming the magnitude denominator. Uncertain fast results use the old interval
checker. Correction residuals also retain the old interval computation. No
binary128 arithmetic is added.

The test-only probe checks fast acceptance and fallback independently. Its
hexadecimal long-double export lets Python compare the reported bound with exact
rational arithmetic without a lossy conversion through binary64.

## Reproduction

Run the full CTest suite and these standalone checks with GCC 14 and Clang 18:

```sh
CC=gcc-14 python3 tests/test_componentwise_certificate.py
CC=clang-18 python3 tests/test_componentwise_certificate.py
python3 tests/test_accuracy_build_flags.py gcc-14
python3 tests/test_accuracy_build_flags.py clang-18
python3 tests/test_certificate_campaign.py
```

Repeat differential tests with `CFLAGS='-O3 -march=native -ffp-contract=fast'`.
Run `tests/kls_componentwise_faults.c` both optimized and under address/undefined
behavior sanitizers. The production build must reject unsafe arithmetic flags.

`bench/qualify_certificate.py` uses explicit paths, not fixed temporary locations:

```sh
python3 bench/qualify_certificate.py run --phase accuracy \
  --support /absolute/path/KLS-paper/evaluation/spice \
  --inventory /absolute/path/inventory.json \
  --baseline /absolute/path/baseline/Xyce \
  --candidate /absolute/path/candidate/Xyce \
  --build-manifest /absolute/path/build.json \
  --output /absolute/path/new-evidence-root \
  --scratch /absolute/path/new-scratch-root
```

The build manifest must contain `baseline_sha256` and `candidate_sha256`; record
compiler flags, link command, source snapshot hashes, and compiler/library
provenance alongside them. The runner pins its support scripts, binaries,
manifest, inventory and recursive circuit inputs. The support directory supplies
the paper's staging, waveform comparison, resource controls and exact audit
helpers; their actual contents are hashed even if the checkout is dirty.

For timings, use `--phase performance` with separate new roots. The baseline and
candidate use the same explicit componentwise policy, compiler settings, and
unchanged numerical kernels. A certificate-object overlay is useful for isolating
the first stage, but its complete link provenance must be retained. Repeat a
stopped command with identical arguments to resume completed observations and
audited chunks. Failed/incomplete simulator runs are preserved and require a new
campaign root; resumption never silently turns them into successful timings.
An accuracy-only retry can use `--cases mos2_large_mem_plus --modes kls8` with
new roots. Its status explicitly reports incomplete full-suite coverage. Such a
retry must be combined with the unchanged, hash-verified prior audit evidence;
it cannot independently authorize performance promotion. The performance phase
rejects these subset options and always covers all seven cases at both thread
counts. Preserve the original runner version when extending a pinned campaign.

Accuracy captures are audited in groups of 20 attempts. Each chunk is compressed,
read back, hash-verified, and receipted before its explicitly listed raw files are
reclaimed. Archives use the evidence filesystem first and scratch second, with
the same 5 GiB reserve. Rejections, waveforms and audit records are retained.
All numerical jobs are serialized with a 24 GiB RSS ceiling and 8 GiB host-memory
reserve. Do not run performance measurements concurrently with builds or audits.

## Gates and remaining stages

The subsequent private-workspace experiment was rejected by its benefit gate.
All 336 timing runs completed; no regression gate failed, but confirmed target
improvements were only 0.52%/0.28% (BJT, one/eight threads), 0.57%/0.11%
(MOS13), and 0.31%/0.15% (mem_plus). The required incremental benefit was 2%.
Its full test suite and standalone fault/concurrency sanitizer checks passed,
but no fresh full capture audit was pursued for an unprofitable candidate.
The implementation was removed rather than retaining extra production state.
The frozen binary, source snapshot, rejected patch and observations remain under
`/tmp/kls-refinement-first/certificate-workspace*` for reproducibility.

The initial row-owned parallel prototype completed 336 comparisons but failed
the no-regression gate: eight-thread BJT slowed by 12.89% in the nine-block
confirmation, while the other six eight-thread cases improved. All waveforms
passed; MOS13 retained one recovered rejection per run and all other cases had
none. That prototype was not promoted; its timing checks were not an exact audit.
Evidence: `/tmp/kls-refinement-first/certificate-parallel-performance/`.

On 2026-09-29, diagnostic BJT runs localized the slowdown outside the short
certificate itself. A timed certificate overlay reported approximately 0.15 s
in parallel checking, while combined refactor/solve time increased much more.
Three-run diagnostic medians (seconds; profiling enabled, not qualification):

| Variant | Simulation | Refactor/solve |
| --- | ---: | ---: |
| Qualified baseline | 5.96735 | 3.12756 |
| Parallel prototype | 7.32322 | 4.47421 |
| Baseline, diagnostic `GOMP_SPINCOUNT=0` | 5.92034 | 3.05860 |
| Prototype, diagnostic `GOMP_SPINCOUNT=0` | 4.94010 | 2.07911 |
| Serial-row diagnostic overlay | 5.54246 | 2.68200 |

This intervention supports competing worker runtimes as the cause: disabling
OpenMP spinning removes the excess refactor/solve cost. The revised candidate
does not change process-wide OpenMP settings. Whenever KLS already owns created
refactor or elimination-graph workers, certification uses serial row traversal
and skips OpenMP timing trials entirely. Other multithreaded solves retain
bounded selection; one-thread solves retain the committed stateless checker.
The rule depends on actual worker ownership, not a circuit name, matrix-size
exception, or new numeric threshold. The accuracy proof is unchanged.

Diagnostic evidence is under `/tmp/kls-refinement-first/parallel-bjt-profile/`
and `parallel-bjt-interference/`. Revised-candidate evidence is under
`certificate-parallel-v2/` and `certificate-parallel-v2-performance/`. The
revised candidate passed both its full performance and exact accuracy gates.
Its nine fresh, alternating-order BJT pairs passed: median candidate/baseline
simulation-time ratio 0.92950 (7.05% faster), with all pairs faster, all waveforms
matching and no failed linear solves. All 37 CTests, GCC OpenMP exact-bound
tests, Clang serial-fallback tests and standalone sanitizer/fault checks passed.
The full seven-case performance campaign passed all 336 comparisons (plus seven
KLU reference runs). Confirmed eight-thread total-time improvements were 6.29%
BJT, 13.30% MOS13, 10.54% mem_plus, 19.83% mux8, 25.70% Rallpack3 and 10.64%
ram2k; DAC slowed 0.57%. Every one-thread change was within 1.29%. All waveforms
passed; MOS13 retained one recovered rejection per run and other cases had none.

Full exact accuracy qualification passed under
`/tmp/kls-refinement-first/certificate-parallel-v2-accuracy/`: all six case/thread
combinations had complete coverage. Each BJT run certified 1,146 successes;
each mem_plus run certified 1,201, without rejections. Each MOS13 run certified
30,490 successes and explicitly rejected one attempt with verified simulation
recovery. In total, 65,674 successes passed the unchanged exact `1e-8` gate,
with two permitted rejections. BJT and mem_plus waveforms matched KLU exactly;
MOS13 passed the waveform tolerance with maximum absolute error approximately
`8e-10`. The frozen production source hashes were checked against the working
tree before promotion. Its storage-only
driver (`qualify_parallel_accuracy.py` in that parent directory) runs the
unchanged exact auditor and may reuse an earlier archive only when every file
name and SHA-256 matches the fresh audited chunk. Before raw reclamation the
original archive helper verifies the archive hash, all archived members, and
all remaining raw files. Driver/runner/receipt hashes and storage provenance
are recorded in `certificate-parallel-v2-accuracy-storage/`; prior acceptance
results are not substituted for fresh audits. Approximately 3.5 GiB of old,
completed BJT raw captures were reclaimed using this verified, recoverable
archive mechanism before launching the new suite.

The standalone parallel-plan regression fixture covers allocation failures,
serial-only requests that must not launch trial workers, worker floating-point
environment rejection, value and pattern changes, thread-count changes, signed
duplicates, and independent concurrent solver contexts. Reproduce its sanitizer
check from the repository root (OpenMP development support required):

```sh
gcc-14 -O1 -g -fopenmp -fsanitize=address,undefined -fno-omit-frame-pointer \
  tests/kls_componentwise_parallel.c -pthread -lm -o /tmp/kls-certificate-test
taskset -c 0-7 /tmp/kls-certificate-test
CC=gcc-14 CFLAGS='-fopenmp -O3 -march=native' python3 tests/test_componentwise_certificate.py
```

- Accuracy: complete exact-rational audit of BJT, MOS13 and mem_plus at one and
  eight threads; every successful answer must satisfy the unchanged tolerance.
  Only MOS13 may have an explicitly matched rejection with verified recovery,
  and an increase above one rejection fails this campaign's promotion gate.
- Waveforms: KLU agreement with `atol=1e-10`, `rtol=1e-8` on all seven cases.
- Timings: three crossed screening blocks followed by nine fresh confirmation
  blocks on the three primary cases plus mux8, DAC, Rallpack3 and ram2k, at one
  and eight threads. Use median paired ratios, not ratios of medians.
- First-stage benefit: at least 10% total-time improvement on each primary
  case/thread combination and no slowdown above 2% anywhere.
- Later stages: independently test workspace reuse, then row-owned parallel
  checking with three paired serial/parallel selection trials. Each retained
  stage must improve at least one primary case by 2%, with no slowdown above 2%
  anywhere; use `--incremental` for this gate.

Performance gates are experiment criteria, not runtime accuracy thresholds.
Workspace reuse and parallel selection are not implemented by the first stage.
Do not begin their production integration until the first stage is qualified.
All retained stages require their own audit and independent commit. Profiling
and competitor comparison follow the last qualified stage; a speedup over the
baseline does not establish that KLS is fastest.
