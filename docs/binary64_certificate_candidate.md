# Guarded binary64 certificate: design and qualification history

Current status: the frozen candidate passed complete incremental SPICE
qualification (see the final sections). The experimental compile guard has
been removed, enabling the same row proof in normal supported builds. Its
GCC 14 object is byte-identical to the qualified object (SHA-256
`b04a2a0bda5f50356113277b48e8406523094a4ca3f0cf49462f81a5d571621b`).
GCC 14 and Clang 20 passed all eight exact-bound test groups without the flag;
normal-build full CTest passed all 40 tests. The earlier failed screens below remain
part of the history and are not replaced by the later successful confirmation.
The profiles show enough certificate cost that reducing arithmetic cost may be
useful in addition to reusing numerical workers, particularly at one thread.

The current proof accumulates in x87 extended precision to avoid exponent-range
problems for arbitrary binary64 inputs. A possible additional fast path could
use binary64 arithmetic only for rows where those problems are explicitly ruled
out, falling back to the existing extended/interval proof otherwise. This must
never bypass certification, change `1/100000000`, or change correction arithmetic.

## Proposed sufficient conditions

- Require binary64 arithmetic, nearest rounding, gradual underflow, and no
  denormals-are-zero mode. On the currently supported x86 fast path, check both
  x87 state and MXCSR; the existing x87-only environment check is insufficient
  for this new arithmetic. Unsupported environments retain the existing path.
- Validate all operands. For this simple first version, decline nonzero
  subnormal inputs, nonfinite inputs, nonfinite products/accumulations, and
  nonzero subnormal products or residual accumulations.
- Decline a zero product of two nonzero operands: it could be underflow, not
  an exact zero. With gradual underflow and representable intermediate operands,
  a zero subtraction is an exact cancellation; nonzero subnormal differences
  are declined. Zero products with a zero operand remain exact.
- That zero-subtraction argument requires a separately rounded product. Prevent
  multiply/subtract contraction in this version, or first supply a different
  proof covering an FMA result that underflows to zero. Merely checking the
  separately computed product does not constrain a contracted residual update.
- Keep the unique-coordinate requirement. Signed duplicates must still take
  the existing coalescing path, not use an inflated magnitude denominator.

## Bound to establish and independently test

Use a conservative relative error of `2u`, with `u=2^-53`, for each accepted
nonzero operation. This also covers a just-subnormal exact product that rounds
to the smallest normal value: its relative error is at most `u/(1-u) < 2u`,
but need not be at most `u`. Zero products from nonzero operands are rejected;
accepted zero subtractions are exact cancellations under gradual underflow.
There are at most `2k` operations along an accumulation path. Put `t=4ku`;
the accumulated bound is `gamma=t/(1-t)`, and converting from exact magnitude
to computed magnitude costs `gamma/(1-gamma)=t/(1-2t)`. The implemented
`c=8*(k+1)*u <= 1/4` ensures `t<=1/8`, hence `t/(1-2t)<=8ku<=c`.
Thus the existing coefficient covers this boundary without changing arithmetic
or a threshold. Its integer conversion is exact within the guard. These are
proof guards, not performance tuning thresholds. A targeted exact-fraction
test exercises `DBL_MIN * nextafter(1,0)` in both orientations.

Evaluate the error, numerator, denominator and reported ratio with outward
rounding. Bound evaluation itself may underflow, so it needs genuine outward
rounding under gradual underflow, not an assumed relative-error model. Use a
lower enclosure of the exact rational tolerance and explicitly require the
reported outward ratio not to exceed that enclosure. Any uncertain row falls
back to the existing proof. Do not replace an unsuccessful fast proof by an
uncertified success.

Before implementation, review this argument against zero/subnormal boundaries,
overflow, cancellation and compiler contraction. Qualification must include
exact-rational differential tests of every reported upper bound, adversarial
tolerance-boundary inputs, FTZ/DAZ and rounding-mode tests, GCC/Clang and FMA
variants, duplicates, both orientations, and the full existing SPICE accuracy
and performance gates. A faster binary64 calculation is not itself evidence
that it supplies a valid certificate.

## First prototype evidence and dispatch limitation

The exact-rational suite (eight test groups) passed with GCC 14 and Clang 20,
OpenMP, native optimization, and contraction enabled. Explicit tests exercise
fast acceptance, range rejection, FTZ/DAZ, and non-nearest rounding. This does
not replace the full captured-SPICE accuracy audit, which has not been run.

The first three-pair screen in
`/tmp/kls-refinement-first/binary64-performance/screen-summary.json` failed:
MOS13-8 median paired runtime increased 14.7%, RAM2k-8 increased 1.0%, and no
primary case met the benefit requirement. Confirmation was automatically
started by the harness; do not interpret the screen as final confirmation.

Important scope correction: `private_componentwise_finish` sends ordinary
one-thread solves directly to `kls_componentwise_certify`, which still uses
the unchanged extended scatter proof. This prototype modifies only the
planned row proof. Thus nearly unchanged one-thread timings do NOT measure
the performance of binary64 arithmetic. They are control measurements.

Before another performance candidate, measure row-proof hit/fallback counts
and exclusive costs in an isolated diagnostic build. The current per-term
classification and volatile operations may outweigh arithmetic savings; this
is a hypothesis, not yet measured. A separate one-thread integration would
require either a guarded scatter proof or a qualified row-plan dispatch change,
with its own allocation, mutation, and numerical tests. Do not silently conflate
that change with the first prototype or weaken any acceptance condition.

## Queued diagnostic workflow

The confirmation campaign is followed serially by these isolated helpers:

- `/tmp/kls-refinement-first/build_binary64_diagnostic.py` snapshots source,
  adds atomic row-proof hit/fallback counters only to that snapshot, and
  relinks a separate executable. It refuses to build while qualification is live.
- `/tmp/kls-refinement-first/run_binary64_diagnostic.py` checks MOS13, RAM2k,
  and BJT at eight threads against the saved KLU waveforms, plus MOS13 at one
  thread as an unchanged-dispatch control. Inputs and executable are hashed.
- `/tmp/kls-refinement-first/after_binary64_confirmation.py` waits for the
  specific confirmation process to exit and requires a completed campaign
  before invoking either helper. A failed performance gate is not itself an
  execution failure: diagnostics remain useful for a completed losing candidate.

The diagnostic results will be written to
`/tmp/kls-refinement-first/binary64-diagnostic-runs/results.json`.
Atomic counters perturb scheduling and execution cost; these runtimes must
never be used to qualify performance. Counts distinguish frequent fallback
from successful-proof overhead, but do not measure exclusive cost or prove
which guard dominates. Waveform agreement does not replace exact per-solve
accuracy qualification. These helpers have been syntax-checked; execution
results must be inspected before claiming diagnostic completion.

## Completed first confirmation and environment diagnosis

The original campaign completed all 343 runs and failed confirmation:
BJT-8 median paired ratio 1.023941, MOS13-8 1.018015, mem-plus-8
1.013939, and RAM2k-8 1.004312. No primary case met the benefit gate.
The earlier MOS13 screen's 14.7% slowdown was not its confirmation result.

The isolated counter campaign completed four waveform-valid runs. It found
zero binary64 hits, with 7,778,520 fallback rows in MOS13-8, 189,298,232 in
RAM2k-8, and 10,309,674 in BJT-8. The one-thread control correctly recorded
no row-proof attempts. A second diagnostic identified the MXCSR guard as
the cause for every attempted MOS13 and RAM2k row: both FTZ and DAZ were set.
Thus this first performance experiment measured an extra guard followed by
the unchanged proof, not the intended binary64 arithmetic.

A second default-disabled prototype now clears FTZ/DAZ once per certificate
worker or serial row traversal and restores the complete saved MXCSR afterward.
It retains the row environment/range guards and the extended fallback.
The eight exact-rational test groups pass with GCC 14 and Clang 20, including
new restoration checks for serial, OpenMP, and executor paths, success and
rejection, and both orientations. This is not full SPICE qualification.

`/tmp/kls-refinement-first/diagnose_binary64_scoped.py` runs the next isolated
counter campaign. Its first MOS13-8 run has 7,778,133 hits and 387 fallbacks,
with a valid reference waveform. Remaining cases and uninstrumented timing
must be inspected before drawing performance conclusions. The one-thread
production dispatch is still unchanged.

The scoped diagnostic subsequently completed all four waveform-valid runs:
RAM2k-8 had 189,298,232 hits and no fallbacks; BJT-8 had 10,318,536 hits
and 41 fallbacks. MOS13-1 again had no row-proof calls. These instrumented
runtimes are not performance evidence.

The uninstrumented follow-up is
`/tmp/kls-refinement-first/screen_binary64_scoped.py`, with artifacts under
`binary64-scoped-candidate` and `binary64-scoped-screen` in the same directory.
It retains the seven-case, two-mode, three-pair screen but explicitly omits
confirmation and labels the result as not full performance qualification.
This avoids spending nine more pairs confirming another losing prototype.
Successful screening still requires separate confirmation and accuracy audit.

The scoped screen completed all 91 runs and failed: no primary benefit,
mem-plus-8 median paired ratio 1.030559 (regression), MOS13-8 1.016232,
BJT-8 1.006675, RAM2k-8 1.008777. A high proof-hit rate therefore does not
establish a speedup. Do not promote this prototype or run confirmation solely
because it now takes the intended arithmetic path.

The next diagnostic, `/tmp/kls-refinement-first/profile_binary64_scoped.py`,
measures inclusive certificate clocks and main-thread stack samples on
mem-plus, RAM2k, and MOS13 at both thread counts. Results are kept separately
under `binary64-scoped-profiles`. Its timing instrumentation is diagnostic,
not a substitute for paired uninstrumented performance qualification.

All 12 scoped profile runs completed with valid waveforms. Inclusive
certificate seconds at eight threads were 0.363963 for mem-plus, 2.900727
for RAM2k, and 0.387330 for MOS13. Earlier qualified-baseline diagnostic
values were 2.661249 for RAM2k and 0.375865 for MOS13. These nonpaired
single diagnostic runs suggest certificate overhead, not a recovery increase:
RAM2k had zero certificate misses and MOS13 retained eight internal misses
and seven correction-RHS calls. One-thread values were essentially unchanged,
as expected from the unchanged dispatch.

Main-thread sampling does not resolve the expensive inner operations: MOS13-8
has only 19 samples and is not aggregate worker profiling. Before another
full SPICE screen, isolate row-kernel arithmetic/guard cost in a bounded
microbenchmark. Per-term classification and volatile intermediates remain
hypotheses, not proven hotspots. Any cheaper range guard must retain exact
zero/subnormal/nonfinite handling and separately rounded products.

## Synthetic row isolation and range-check candidate

`certificate_row_microbench.c` and `run_row_microbench.py` under the same
temporary evidence directory compare six row lengths (1–1024), both
orientations, five samples per invocation, in extended/binary64/binary64/extended
order on CPU 0. Inputs are synthetic, cache-resident, exact dyadic rows;
these timings neither represent the full circuit corpus nor qualify accuracy.
The scoped floating-classification prototype is 5–17% slower on short rows
and roughly tied on long rows (`row-microbench-results.json`).

The next experimental change replaces repeated floating classification with
a memcpy-based binary64 magnitude-bit check accepting only zero or finite
normal values. It retains every operand/intermediate check, zero-product
underflow rejection, and volatile separately rounded arithmetic. Tests compare
the classifier against the original predicate for both signs, all 2048
exponents, and four representative fraction patterns, alongside exact-rational
certificate tests and environment restoration tests.

`row-microbench-bitguard-results.json` shows candidate/extended ratios of
0.85–0.89 at degree four, 0.74–0.77 at sixteen, and 0.55–0.63 at 64–1024;
degree one remains about tied. This isolates classification cost as a useful
optimization target for these synthetic inputs. A new SPICE screen and full
accuracy qualification are still required; do not infer circuit gains.

The bit-guard SPICE screen completed 91 runs with mixed results: MOS13-8
median paired ratio 0.945419, mem-plus-8 0.967165, RAM2k-8 0.985476,
but DAC-8 1.024617 fails the no-regression gate. Its DAC pairs are
0.910617, 1.053411, and 1.024617; do not dismiss the regression as noise
without further evidence. The overall screen therefore failed despite
meeting the primary benefit requirement.

`confirm_binary64_bitguard.py` starts a separate nine-pair confirmation
across all seven cases and both modes, using exactly the frozen screened
binary. It does not rebuild or modify the candidate. Confirmation results
are separate from the failed screen, which remains preserved. No accuracy
qualification or production promotion is implied by launching confirmation.

Nine-pair confirmation completed all 259 runs and passed both benefit and
no-regression gates. Eight-thread median runtime changes versus the qualified
baseline were: BJT -1.240%, MOS13 -6.027%, mem-plus -2.170%, mux8 -0.069%,
DAC -3.558%, Rallpack3 -1.823%, RAM2k -1.155%. One-thread controls ranged
from -0.057% to +0.359%. The failed initial screen is retained; DAC's slowdown
did not persist in confirmation. These are improvements over KLS's baseline,
not evidence of beating every competitor or achieving the overall goal.

Captured-solve accuracy qualification has started with both BJT modes under
`binary64-bitguard-bjt-accuracy`. This is a scoped audit, not the full accuracy
suite. MOS13 and mem-plus still require their own captured-solve audits.
Storage is limited; preserve the capture reserve and archive verification.

The first BJT capture stopped at the 5 GiB storage reserve, not on a numerical
failure. Its partial exports were SHA-256 verified into
`/dev/shm/kls-bitguard-interrupted-evidence/bjt1.tar.gz` before raw removal;
the failed campaign remains failed. A fresh root is required for retry.
Recoverable compression of 272 competitor and 2,401 older performance waveforms
raised shared-memory free space to 8.3 GiB. Receipts are
`/tmp/kls-refinement-first/completed-waveform-compression.json` and
`historical-performance-waveform-compression.json`; `gzip -dk` restores each
original waveform. No accuracy acceptance condition was changed.

The GCC 14 exact-bound suite was rerun and all eight tests passed. The full
CMake test build is now being rebuilt with the experimental macro enabled;
its logs are `/tmp/kls-refinement-first/bitguard-ctest*.log`. Launching these
tests is not a pass or a substitute for captured-SPICE qualification.

Full macro-enabled CTest initially passed 39/40: the pool fixture required a
bit-identical upper bound from two different proof kernels. For its exact-zero
integer residual fixtures, it now checks the tolerance and compares pooled
execution with serial execution of the same planned proof. The original
stateless proof must still certify. All 40 tests then passed (see
`bitguard-ctest-retry.log`). The macro-enabled ASan/UBSan row-plan fault,
selection, lifecycle, and independent-concurrency test also passed.

After verified archival of both historical Rallpack captures (23,137 files
each) and hash-verified relocation of retained archives, shared-memory free
space reached approximately 11 GiB. BJT qualification was restarted in fresh
roots `binary64-bitguard-bjt-retry` under `/tmp/kls-refinement-first` and
`/dev/shm`, with the same frozen candidate and both thread modes. The original
storage-stopped campaign remains failed; the retry is not yet a qualification
pass. Old archive paths remain usable through symlinks to disk copies.

The retry's BJT one-thread result passed all 1,146 captured attempts with zero
explicit rejections and an exact KLU waveform match. Its eight-thread simulation
also completed and passed waveform comparison; the exact audit remains in
progress. See `binary64-bitguard-bjt-retry/PDE_2D_BJT_invertbjt_cir-kls1-audit.json`.
Do not infer the eight-thread proof result from the unchanged one-thread path.

BJT retry qualification subsequently completed and passed both modes: 1,146
certified successes per mode, zero explicit rejections, complete capture
coverage, and exact KLU waveform matches. This is 2,292 audited successes,
not full-suite qualification. The serialized MOS13 campaign has now started
under `binary64-bitguard-mos-accuracy`; the memory case remains pending.

MOS13 one-thread qualification passed with complete coverage: 30,490 certified
successes and one explicit rejection, with verified recovery and maximum
waveform difference approximately 8e-10. Its eight-thread capture completed
and is undergoing exact audit. A retained full memory-case archive contains
5,860,485,513 raw bytes (5.458 GiB), so its fresh capture needs over 10.458 GiB
free plus working headroom to preserve the unchanged 5 GiB capture reserve.

MOS13 qualification completed successfully in both modes: each reports 30,490
certified successes and one explicit rejection with verified simulation
recovery; both waveform comparisons passed with maximum difference about
8e-10. Together with BJT, this establishes 63,272 certified successes and two
permitted recovered rejections. Memory qualification still remains; this is
not yet the complete candidate accuracy gate.

Both memory modes subsequently passed: 1,201 certified successes each, zero
rejections, and exact KLU waveform matches. The combined verifier at
`/tmp/kls-refinement-first/verify_bitguard_qualification.py` checked all six
case/thread results, recomputed every successful residual fraction against
1/100000000, verified all 3,288 retained archive hashes, matched binary/input/
tool provenance, and recomputed the nine-pair performance medians. Its result
`binary64-bitguard-combined-qualification.json` reports 65,674 certified
successes and two permitted recovered rejections. Thus the frozen experimental
candidate passes incremental SPICE qualification. Production promotion and
fresh full-corpus competitor comparisons remain separate steps; this does not
establish that KLS is the fastest solver.
