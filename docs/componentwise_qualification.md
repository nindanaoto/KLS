# Componentwise acceptance qualification

## Contract and compatibility

The default policy is `KLS_ACCURACY_COMPONENTWISE_BACKWARD_ERROR`. Each
successful RHS satisfies, in the original caller coordinates,

`|b-Ax|[i] <= (|A||x|+|b|)[i] / 100000000` for every row.

A zero denominator requires an exactly zero residual. Signed duplicate
coordinates are combined before forming `|A|`. A conservative hardware
long-double interval certificate may decline to certify a mathematically
acceptable answer; it must not accept an answer outside the bound.

The policy tries the raw binary64 answer and at most three ordinary
double-factor corrections. Exhaustion returns failure without publishing a
partial answer, including for aliased and multiple-RHS calls. No heavy recovery
is enabled by this policy. Thresholds and correction budgets were not raised to
obtain the results below. The maxima optimization removes library calls without
changing the outward-rounded bounds or introducing circuit-specific conditions.

`KLS_ACCURACY_STRICT_RHS_L2` preserves the previous RHS-relative L2 contract
(absolute L2 for a zero RHS), including its cold recovery. Applications requiring
that contract must select it before the first successful analysis. Temporary
strict-recovery solvers explicitly inherit their parent's policy. The public
options and statistics layouts are unchanged.

Xyce accepts `KLS_ACCURACY_POLICY=DEFAULT`, `STRICT_RHS_L2`, or
`COMPONENTWISE_BACKWARD_ERROR` and reports the resolved policy. Omitting the
option is equivalent to `DEFAULT`. Updating adapter options destroys the old
solver and applies the new selection before analysis.

An explicit `KLS_ACCURACY_APPLICATION_MANAGED` opt-out is also available;
Xyce selects it with `KLS_ACCURACY_POLICY=APPLICATION_MANAGED`. It skips
residual certification and certificate-driven corrections/recovery, but retains
kernel failure reporting, finite input/output checks, and atomic publication.
It does not promise a residual bound. The default remains componentwise.
The 37-test suite covers the opt-out API and failure behavior. Initial opt-out
screens matched BJT and mux8 waveforms, but MOS13 failed DC operating-point
convergence at both thread counts (25 linear attempts, zero reported linear
failures). This failure is not a fast timing result and blocks default opt-in
promotion. Evidence: `/tmp/kls-refinement-first/application-managed/`.

## Numerical qualification, 2026-09-27

The optimized, explicit-componentwise candidate passed full exact-rational
audits at both one and eight threads:

| Circuit | Certified successes per run | Reported rejections | KLU waveform |
| --- | ---: | ---: | --- |
| BJT invertbjt | 1,146 | 0 | Identical |
| MOS13 invert50 | 30,490 | 1 | Within declared tolerance |
| mem_plus | 1,201 | 0 | Identical |

Each MOS13 rejection was uniquely matched between failed-system exports and
caller-system exports; simulator counts and warnings agreed. Every other
answer passed the unchanged componentwise bound. The simulation completed and
matched KLU: maximum absolute waveform difference approximately `8e-10`.
This is a user-approved **simulation recovery** allowance, not permission to
return an uncertified successful linear answer. BJT and mem_plus had no waiver.

Additional voter25, RC3000, ISSUE 575 and passivePatch waveform runs passed at
both thread counts. Those additional runs are not per-answer exact audits.
Waveform comparison uses `atol=1e-10`, `rtol=1e-8`; maximum absolute error alone
is not the waveform acceptance rule.

The certificate differential suite checks 1,500 random systems, 600 generated
near-solutions, and boundary/extreme/subnormal/duplicate-cancellation fixtures
against rational arithmetic. Certificate allocation and floating-point-mode
fault tests passed separately with GCC 14 and Clang 18. API tests cover
orientation, normal/transpose, CSC/CSR, scaling/BTF, numeric updates, aliasing,
strides, policy locking, and failed-output preservation. Caller fault tests
exercise first/second-RHS certificate OOM and correction exhaustion and trap
unexpected external lattice recovery. Forced exhaustion issues three correction
requests and seven certificate calls, without entering the strict certificate.

The fresh default-policy build passed all 36 CTest tests. Legacy fixtures that
assert strict-specific recovery or verified-RHS reuse explicitly select strict
mode; their assertions were retained. Xyce smoke tests passed for an omitted
policy, `DEFAULT`, explicit componentwise, and explicit strict. Default BJT,
MOS13 and mem_plus runs at both thread counts reproduced the previously
audited candidate's waveforms byte-for-byte, with identical solve and rejection
counts. These default smoke tests are not additional full linear-system audits.

## Performance scope

The optimized explicit policy passed 192 serialized control simulations:
three crossed screening blocks followed by nine fresh confirmation blocks,
four circuits, two thread counts, and two policies on the same frozen binary.
Every waveform passed and no linear rejection occurred. The acceptance gate
was no confirmed median paired slowdown above 2% on any control.

| Circuit | One-thread improvement vs strict | Eight-thread improvement vs strict |
| --- | ---: | ---: |
| mux8 | 4.35% | 4.28% |
| DAC | 0.52% | 0.47% |
| Rallpack3 | 2.93% | 3.12% |
| ram2k | 2.81% | 2.71% |

These are medians of paired elapsed-time ratios, not ratios of medians. They
isolate policy overhead; they do not establish cross-machine performance or
prove superiority to other solvers. Accuracy-capture timings are not performance
measurements. Prior paper timings must not be relabelled as results of this
implementation. Other solvers' waveform agreement does not establish that they
enforce this same linear accuracy contract.

The subsequent fresh default-build screen completed 60 simulations, all with
valid waveforms: three runs per available circuit/mode combination. Median
full-simulation seconds follow. SubtreeLU-8 MOS13 was stopped after 128 seconds
without reported operating-point progress and is excluded, not assigned a time.

| Circuit | KLU | KLS-1 | KLS-8 | CKTSO-1 | CKTSO-8 | SubtreeLU-1 | SubtreeLU-8 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| BJT | 51.032 | 10.555 | 7.634 | 7.450 | 4.349 | 10.057 | 4.243 |
| MOS13 | 0.973 | 1.959 | 1.952 | 0.783 | 1.684 | 0.823 | Excluded |
| mem_plus | 4.971 | 7.958 | 7.952 | 3.833 | 3.623 | 4.383 | 3.544 |

These are three-block discovery screens, not nine-block confirmations. KLS is
not the fastest solver on these cases. No claim of equivalent successful-answer
accuracy is made for the competitors; their answers were not independently
certified in this timing screen. The remaining performance gap is future work.

## Reproduction and evidence

Build with GCC 14 or later and run `ctest --test-dir BUILD --output-on-failure`.
Tests require strict floating-point semantics; fast-math/finite-only builds are
rejected. Clang 18 was also tested for the standalone certificate checks.

The local development evidence is under `/tmp/kls-refinement-first/`:

- `componentwise-maxima-requalification/`: full optimized-candidate audits.
- `componentwise-maxima-controls/`: 192-run policy comparison and paired ratios.
- `componentwise-maxima-faults/`: GCC/Clang certificate fault/differential tests.
- `componentwise-broader/`: additional waveform checks and 36-test suite.
- `componentwise-promotion-v2/`: fresh default build, 36-test suite, policy
  smoke tests and the start of the competitor screen. Its stopped status is
  from a manually terminated SubtreeLU-8 MOS13 run, not a KLS test failure.
- `componentwise-promotion-comparison-v2/`: continuation of the screen, with
  SubtreeLU-8 MOS13 explicitly excluded after 128 seconds without reported
  operating-point progress. The cause is unresolved; this is neither an OOM
  finding nor a completed performance measurement. Its other repetitions were
  not attempted. The continuation completed all 60 nonexcluded runs with valid
  waveforms; its terminal status records the declared exclusion.

These local paths are provenance, not a portable archived paper artifact.
Capture archives have member hashes and receipts; some are RAM-backed. The
last mem_plus archival hit a tmpfs reserve **after its audit passed** and was
recovered into a verified disk archive; `archive-recovery.json` records that
recovery without hiding the original stopped status.

## Interval-row numerical qualification, 2026-10-08

The current worktree stops a planned fast proof after a supported row rejects
it. The full interval pass still computes the correction RHS. A rejected
extended row proof also avoids repeating the same accumulated CSC proof;
unsupported workers or disagreeing timing trials retain the stateless path.

The full pass can now use the cached, duplicate-free row traversal in
`src/kls_componentwise_interval_rows.inc`. Stable scatter preserves each row's
CSC encounter order. Coefficient bounds, product bounds, endpoint and
denominator updates, tolerance comparisons, and correction midpoints retain
the original interval arithmetic. This removes coalescing scratch arrays and
residual scatter for eligible patterns. Duplicate entries, unsupported
environments/platforms, and a correction workspace overlapping any input array
retain the original CSC pass. Eligibility does not depend on matrix names or
performance thresholds. The row plan caches only traversal indices, never
values, RHSs, answers, or certificates.

After fixing a caller-MXCSR restoration failure in timing calibration, the
candidate passed all 42 CTest tests, both OpenMP-disabled certificate checks,
nine exact-rational audit tests with forced row/CSC bound and usable correction
comparisons, and AddressSanitizer/UndefinedBehaviorSanitizer checks on the raw
correction and parallel tests. Xyce rebuilt successfully. The initial failed
validation is retained under `build/early-rejection-20261008/failed-validation-mxcsr/`.
Added cases cover workspace overlap, subnormal corrections
with FTZ/DAZ, nonfinite inputs, and unrepresentable correction RHSs. The
acceptance tolerance and three-correction budget remain unchanged. The
completed frozen benchmark used the previous executable; its timings cannot
establish a speedup for these changes. All fourteen full-circuit phase profiles
passed the original waveform and step-count gates with no input drift.
Uninstrumented paired performance qualification remains pending.

## Binary64 register arithmetic qualification, 2026-10-08

A separate diagnostic sampled every 1024th certificate call when accepted,
replaying the validated binary64 row proof after stopping the certificate
timer and restoring MXCSR. All seven original circuits, with one and eight
threads, passed waveform and step-count gates with no input drift. All 598,912
sampled rows passed the binary64 proof. This is arithmetic-path attribution;
instrumented timings do not establish a solver speedup.

The binary64 row proof now uses explicit scalar SSE2 instructions instead of
volatile temporary-memory intermediates. This preserves separately rounded
multiplication, addition, subtraction, and division, including when the
including translation unit enables FMA. The operation order, range guards,
outward bounds, tolerance, and extended-precision fallback are unchanged.
The implementation remains inside the existing supported x86-64 proof gate;
other platforms retain their existing paths. No matrix-specific rule or new
performance threshold was added.

The candidate passed 42/42 CTest tests, both OpenMP-disabled certificate
checks, all nine exact-rational audit tests, and the raw correction/parallel
AddressSanitizer, leak, and UndefinedBehaviorSanitizer checks. A separate
audit compared the new and previous binary64 proof verdicts and hexadecimal
extended-precision bounds across all differential cases; all agreed. That
audit also passed with Intel assembler syntax. Xyce rebuilt successfully.
Validation artifacts are under `build/binary64-register-20261008/`.

The uninstrumented KLS-only experiment completed all 238 observations on the
seven original circuits, both cache domains, one/eight threads, and four
balanced passes with fresh KLU references. Every observation passed the
original correctness gates, and no input drift occurred. Candidate/baseline
geometric means of per-circuit median paired linear-solution-time ratios were
0.939621 (LLC32, one thread), 0.958537 (LLC32, eight threads), 0.935266 (LLC96,
one thread), and 0.956806 (LLC96, eight threads). These results compare KLS
register arithmetic with the frozen KLS interval-row executable. Current
comparison against CKTSO and SubtreeLU remains a separate requirement.

## Native-format row maximum qualification, 2026-10-08

An accepted binary64 row proof returns a finite nonnegative binary64 bound.
Converting such bounds to extended precision is exact, so taking their maximum
before conversion gives the same result as converting and reducing each row.
The row-range executor now retains a binary64 maximum, an extended maximum for
rows that need the extended proof, and converts the binary64 maximum once at
the end of a successful range. Worker partitioning and early rejection remain
unchanged. A rejected range still requests the original complete certificate
and correction RHS. No matrix-specific rule, performance threshold, tolerance,
or accuracy policy changed.

Independent comparisons against the previous compiled proof matched scalar
and planned verdicts and hexadecimal extended-precision bounds across all
audit cases, under both assembler dialects. New cases compare the range with
the previous scalar-row reduction for mixed binary64/extended rows, overflowing
binary64 magnitude sums, subnormal inputs, empty worker ranges, FTZ/DAZ, changed
SSE rounding controls, transpose, serial execution, OpenMP workers, and the
executor callback. The candidate passed 42/42 CTest tests, both OpenMP-disabled
checks, ten exact-rational audit tests, and the raw correction/parallel ASan,
leak, and UBSan checks. Xyce rebuilt successfully; source hashes stayed stable
during validation. Artifacts are under `build/binary64-row-maximum-20261008/`.
The uninstrumented KLS-only paired comparison completed all 238 observations
with every correctness gate passed and no input drift. Candidate/baseline
geometric means of per-circuit median paired linear-solution-time ratios were
0.905268 (LLC32, one thread), 0.951784 (LLC32, eight threads), 0.902512 (LLC96,
one thread), and 0.949832 (LLC96, eight threads). The baseline was the frozen
register-arithmetic executable. These measurements support an incremental KLS
improvement on the original circuit suite; they do not establish a current
win against CKTSO/SubtreeLU or broad-matrix performance.


## Caller-operator gather qualification, 2026-10-08

Cached caller-operator updates now copy each value's object representation
through a fixed-size byte-copy helper. This avoids the per-entry fortified
library call generated by the current compiler. The accuracy-policy branch
is outside the gather loop. Fresh allocation, exact structural/map comparisons,
source-map validation, and transactional publication retain their existing
contracts; an update may read the snapshot it replaces. There is no matrix
recognition, size-based dispatch, tolerance change, or ISA-specific kernel.

The expanded cache contract checks exact bytes for signed zeros, subnormals,
infinities, and quiet/signaling NaN payloads under all four rounding modes,
both orientations, and all permutation/map combinations. Native and generic
optimized builds and ASan/UBSan/leak checks passed. The integration build passed
42/42 CTest tests, and Xyce rebuilt successfully without source drift during
validation. The standalone gather experiment retained all 288 observations
across six sizes, three mappings, two policies, and four alternating passes;
all timed batches agreed byte-for-byte with the previous implementation.
For 1,024 through 1,048,576 entries, median candidate/baseline time ratios
ranged from approximately 0.273 to 0.652. These are isolated snapshot-update
measurements; SPICE and competitor qualification remains pending. Artifacts
are under `build/operator-gather-20261008/`.

The completed uninstrumented comparison retained all 462 observations on the
seven original circuits, two cache domains, one/eight threads, and four
balanced passes, including fresh KLU references. KLS and CKTSO observations
passed the unchanged correctness gates. Four SubtreeLU eight-thread MOS
observations reached the unchanged 1,800-second timeout; they remain in the
record. No input drift occurred. Geometric means of per-circuit median paired
linear-solution-time ratios (candidate divided by challenger) are:

| Cache domain | Threads | Previous KLS | CKTSO | SubtreeLU |
| --- | ---: | ---: | ---: | ---: |
| LLC32 | 1 | 0.907639 | 1.827398 | 1.049348 |
| LLC32 | 8 | 0.890062 | 0.699624 | unqualified |
| LLC96 | 1 | 0.897090 | 1.816522 | 1.103008 |
| LLC96 | 8 | 0.868810 | 0.859324 | unqualified |

These results establish a 9–13% aggregate improvement over the previous KLS
on this circuit suite. KLS remains slower than CKTSO and SubtreeLU in serial.
An all-circuit parallel SubtreeLU aggregate cannot be qualified with the
retained timeouts. No cases were removed to compute that aggregate. Broad
matrix performance and accuracy remain separate qualification requirements.
The full receipt is `build/xyce-operator-gather-campaign-20261008/summary.json`.

## Isolated binary64 budget candidate, 2026-10-08

The isolated candidate replaces per-entry normality checks with a conservative
row error budget that covers gradual underflow. Positive binary64 encoding
steps add the absolute underflow allowance without evaluating a subnormal
budget product. Unsupported environments and declined rows retain the existing
proof fallback. The tolerance, correction budget, duplicate handling, and
matrix-independent dispatch contract remain unchanged. This candidate is not
promoted to the production source tree.

The isolated integration passed 42/42 CTest tests, independent exact-rational
row audits in generic, native, and Intel assembler builds, OpenMP-disabled
checks, and correction/parallel sanitizer checks. The Xyce link used the tested
KLS archive and verified unchanged remaining link inputs. The complete SPICE
campaign retained all 462 observations with no input drift. All KLS and CKTSO
observations passed the unchanged waveform and step-count gates. One SubtreeLU
eight-thread MOS observation reached the unchanged 1,800-second timeout and
remains in the record. Candidate/challenger geometric means of per-circuit
median paired linear-solution-time ratios are:

| Cache domain | Threads | Caller-gather KLS baseline | CKTSO | SubtreeLU |
| --- | ---: | ---: | ---: | ---: |
| LLC32 | 1 | 0.774593 | 1.415520 | 0.814124 |
| LLC32 | 8 | 0.869958 | 0.609537 | unqualified |
| LLC96 | 1 | 0.763123 | 1.384074 | 0.836788 |
| LLC96 | 8 | 0.863529 | 0.786523 | 1.352972 |

The aggregate baseline improvements do not qualify promotion: serial CKTSO
and parallel SubtreeLU gaps remain. Eight-thread NEURON regressed against
baseline in both domains; the LLC96 median ratio was 1.124443, with all four
pairs slower. Its cause is unproven. The LLC32 parallel SubtreeLU aggregate
remains unqualified; no failed cases were removed. Full paired validation of
all 110 original broad-corpus matrices is running separately, retaining every
failure and timeout. No full-corpus conclusion is available yet.

The complete SPICE receipt is
`build/xyce-budget-bitsteps-campaign-20261008/summary.json`; the promotion
decision and independent integration evidence are under
`build/binary64-budget-bitsteps-integration-20261008/`. Broad validation is
under `build/binary64-budget-bitsteps-broad-20261008/`.

### Completed broad validation, 2026-10-09

The unchanged 110-matrix campaign completed all 1,760 observations across
both cache domains and four alternating paired passes. Input hashes remained
unchanged, with no candidate-only coverage regressions. All 80 invalid
observations remain: both arms failed every pass on hamrle3 (timeout), fullchip
(solve failure), bips07_1998 and fpga_dcop_01 (residual gate), and freescale2
(singular setup). Neither domain qualifies an all-matrix aggregate; no
successful-subset aggregate was substituted.

Five completed groups have candidate/baseline median ratios above one:
LLC32 mc2depi 1.004787; LLC96 g3_circuit 1.000388, mc2depi 1.012495,
htc_336_9129 1.007977, and mult_dcop_03 1.020361. Individual slower passes
also remain, including the LLC96 pre2 fourth-pair ratio of 1.282291 despite
its median of 0.997193. Its recorded refactor time increased; the cause is
unproven. These results do not qualify production promotion or establish
the overall competitor goal. The complete receipt is the broad
`summary.json`; `completed-broad-decision.json` records failures and median
regressions without dropping cases.

### Warm-selection comparison (2026-10-09): rejected

The isolated warm-selection candidate completed all 462 original SPICE
observations with unchanged input hashes. It retains one agreeing warm pair
before the existing three measured pairs and median rule. All KLS waveform,
step, end-marker and linear-failure checks passed. Four SubtreeLU MOS
eight-thread timeouts remain in the record: LLC32 pass 3 and LLC96 passes
0, 2 and 3, each under the original 1,800-second deadline.

Ratios below are candidate time divided by challenger time, using all seven
cases and four paired passes per case. Lower is faster.

| Domain | Threads | Integer baseline | CKTSO | SubtreeLU |
| --- | ---: | ---: | ---: | ---: |
| LLC32 | 1 | 1.009218 | 1.426657 | 0.818675 |
| LLC32 | 8 | 0.994948 | 0.604555 | unqualified |
| LLC96 | 1 | 1.007205 | 1.388020 | 0.842724 |
| LLC96 | 8 | 1.019206 | 0.782882 | unqualified |

The candidate regresses serially in both domains and at eight threads in
LLC96. It is rejected for promotion; its prepared broad campaign is not
launched. All slower groups remain, including NEURON at eight threads:
candidate/baseline medians 0.923147 in LLC32 and 1.098181 in LLC96. No
case-specific or domain-specific switch is introduced. The serial CKTSO gap
and overall competitor goal remain unresolved.

A separate revision restores the original 128-byte plan and 32-byte
direction stride by encoding warm state in the existing choice field. Its
11 public, cold-start, sanitizer and layout checks passed without input
drift. This establishes correctness for those checks, not performance; it
cannot inherit the rejected candidate's timings. Complete records are
`build/xyce-warm-selection-campaign-20261009/summary.json` and
`build/componentwise-warm-selection-integration-20261009/completed-spice-decision.json`.

### Cached supernodal destination runs (2026-10-09): rejected

The isolated cached-run update experiment completed all 103,680 observations
across 12,960 groups without input drift. The baseline uses the original
112-byte block layout; the candidate appends a run pointer, preserving all
original field offsets. Six arithmetic-audit configurations passed, but the
complete timing grid regressed:

| Build | Candidate/baseline geometric mean | Slower groups |
| --- | ---: | ---: |
| Generic | 1.116733 | 5,473 / 6,480 |
| Native | 1.331564 | 5,725 / 6,480 |

All twelve native destination-map class aggregates are slower. The
implementation is rejected; no case, build, shape, or winning-subset exception
is selected. Classification, preparation, and full production cache allocation
and lifetime costs are outside these isolated timings, so they cannot establish
a complete-solver gain. No cached-run source is promoted. This experiment is
independent of the caller-workspace reuse candidate and supplies no timing
evidence for it. The complete records are
`build/snb-cached-run-layout-corrected-performance-20261009/summary.json` and
`build/snb-cached-run-layout-corrected-decision-20261009/decision.json`.


### Caller workspace reuse (2026-10-10): rejected for promotion

The fresh recovery campaign completed all 462 observations with no input drift.
All 224 KLS observations qualified. Three SubtreeLU MOS eight-thread runs timed
out: LLC32 passes 0 and 2, and LLC96 pass 0. They remain in the record, and both
eight-thread SubtreeLU aggregates are unqualified. The interrupted original
campaign is retained separately and supplies no pooled observations.

Ratios below are candidate/challenger geometric means of the seven case
medians, each requiring all four qualified pairs; lower is faster.

| Domain | Threads | Integer baseline | CKTSO | SubtreeLU |
| --- | ---: | ---: | ---: | ---: |
| LLC32 | 1 | 1.004340 | 1.421838 | 0.815510 |
| LLC32 | 8 | 1.000653 | 0.647960 | unqualified |
| LLC96 | 1 | 0.999384 | 1.373590 | 0.837285 |
| LLC96 | 8 | 0.972854 | 0.729315 | unqualified |

Three baseline comparisons are near parity. The LLC96 eight-thread gain is
concentrated in NEURON (median 0.838482); its LLC32 median is 1.005283.
Serial CKTSO remains faster on every circuit in both domains. Workspace reuse
is rejected for promotion, with no case or domain switches. All paired ratios,
regressions and failed observations are preserved. Its prepared broad campaign
is not launched; no broad-generalization or overall competitor claim follows.
Native tests and sanitizer/lifecycle checks passed, but those checks do not
establish speed. Private caller buffer retention is also outside public KLU
memory accounting. The next independent prototype borrows caller RHS storage
and requires its own full correctness and performance qualification.

Complete records: `build/xyce-caller-workspace-cold-cache-recovery-campaign-20261009/summary.json`
and `build/caller-workspace-cold-cache-recovery-decision-20261010/decision.json`.

### Borrowed caller RHS (2026-10-10): rejected for promotion

The independent campaign completed all 462 observations. The actual controller
process was absent before summary review. All 224 KLS observations qualified;
460 observations qualified overall. Two SubtreeLU MOS13 eight-thread trials
timed out at the original 1800-second deadline: LLC32 pass 1 and LLC96 pass 0.
Both eight-thread SubtreeLU aggregates remain unqualified. No observations were
retried, excluded, or inherited from workspace reuse. Current input, controller,
binary, and archive hashes were checked, with no drift. Exact unique protocol
coverage and all twelve comparison blocks were independently recomputed.

Ratios are candidate/challenger geometric means of the seven case medians,
each requiring all four qualified pairs; lower is faster.

| Domain | Threads | Integer baseline | CKTSO | SubtreeLU |
| --- | ---: | ---: | ---: | ---: |
| LLC32 | 1 | 0.993796 | 1.405155 | 0.808350 |
| LLC32 | 8 | 0.998271 | 0.627671 | unqualified |
| LLC96 | 1 | 1.000170 | 1.373383 | 0.838533 |
| LLC96 | 8 | 1.012972 | 0.762829 | unqualified |

Borrowing does not establish a consistent improvement over the integer
baseline. The LLC96 eight-thread NEURON median is 1.091566, versus 0.983346
on LLC32. Serial CKTSO remains faster on every circuit in both domains.
Qualified SubtreeLU case regressions and every paired ratio remain recorded.
Promotion is rejected and the prepared broad campaign remains unlaunched.
There are no case or domain switches and no broad-generalization claim.

All four baseline/candidate native and O3 ASAN/UBSAN suites passed 44 tests,
315 guard/overlap/recovery cases, and 2430 shape cases. The shape fixtures
reported paths 0, 3, and 5; this is not an all-engine claim. The original O1
baseline EGraph bitwise fixture failure remains unresolved and retained.
Correctness qualification does not establish competitive performance. KLS
has not achieved the fastest-SPICE objective.

Complete records: `build/xyce-caller-rhs-borrowing-production-optimization-campaign-20261010/summary.json`
and `build/caller-rhs-borrowing-production-optimization-spice-decision-20261010/decision.json`.


## Direct indexed narrow-edge vector update (2026-10-10)

Rejected for promotion. The fixed structural grid completed all 103680
observations after the actual controller process exited. Independent review
verified all 12960 shape/build groups, four paired ratios per group, all 26
aggregates, balanced AB/BA order, byte and floating-point flag equality,
and frozen source, controller, executable, and Apptainer image hashes.
Both arms retained the original block layout and narrow-width threshold;
the candidate added eight-row masked gather/scatter without extra metadata.

| Build | Candidate / baseline geometric mean | Slower groups / 6480 |
| --- | ---: | ---: |
| Generic | 0.999793 | 3254 |
| Native | 1.602533 | 6417 |

Every native map class regressed. All observations and regressions are
retained; no fitted cutoff or map selector is introduced. Four extended
arithmetic builds each passed 616320 fixtures, including exceptional values,
rounding modes, legal injective maps, and guard pages. These are isolated
kernel checks, not full-solver or all-engine qualification. Production is
unchanged and full-solver testing of this prototype is rejected. The
fastest-SPICE objective remains unmet.

Complete records: `build/snb-direct-indexed-vector-performance-controller-20261010/summary.json`
and `build/snb-direct-indexed-vector-performance-decision-20261010/decision.json`.


### Scalar four-way narrow-edge update (2026-10-10)

The isolated scalar four-way unroll of the narrow SNB below-row loop is rejected.
It preserves row and producer order, zero skips, the original 112-byte block
layout, and the existing narrow-width boundary. It adds no map classification,
metadata, or benchmark-specific selector. Four correctness builds (generic,
native, separate arithmetic, and ASan/UBSan) each passed 616,320 fixtures with
full-panel byte equality, floating-point flag equality, and guard pages.

The actual timing writer (PID 3650009) was absent before review. All 103,680
observations cover the complete fixed grid, with four balanced paired trials
per shape in each build and 2,048 calls per observation. An independent review
verified the exact grid and ordering, recomputed all 12,960 groups and 26
aggregates, and checked frozen source, controller, executable, and Apptainer
image hashes. Every observation passed the byte and flag oracle.

| Build | Candidate / baseline geometric mean | Slower groups / 6,480 |
| --- | ---: | ---: |
| Generic | 1.0789286573170802 | 4,895 |
| Native | 1.1448015940758456 | 4,885 |

All twelve native map classes regress in aggregate. All observations and
regressions remain in the evidence; no fitted cutoff or selector is introduced.
The decision is recorded in
`build/snb-narrow-scalar-unroll4-performance-decision-20261010/decision.json`.
This is isolated kernel evidence, with clock and call overhead included; it
supports no full-solver, all-engine, or competitive SPICE claim. The prototype
is not promoted and no full-solver campaign is launched for it. KLS being
fastest against CKTSO and SubtreeLU remains unproven.


### Actual SNB invocation attribution (2026-10-10)

The isolated diagnostic in `build/snb-invocation-attribution-20261010` rebuilt
only the numeric owner object from the frozen integer-certificate baseline.
The other six archive members retain identical hashes. It adds owner-entry
scope timers to serial and parallel SNB refactor, including all early returns;
it adds no edge or row loop counters. The existing public API wrappers retain
final-path counts and certificate attribution. Parallel fallback can nest a
serial SNB call, so these timer totals must not be added together.

The actual controller (PID 3997741, session 92355) exited successfully before
review. All 28 configurations (seven original full circuits, two CPU domains,
one and eight threads) passed waveform, step, end-of-simulation, and diagnostic
checks. Independent review verified the unique complete grid, parsed both
profile records from every log, and checked all frozen preparation, build,
and campaign input hashes. Results and provenance are recorded in
`build/snb-invocation-attribution-20261010/attribution-review.json`.

In both domains, serial BJT invoked SNB 1,145 times and reported 1,144 final
public SNB refactors. Serial SNB entry wall time was 4.317 s in LLC32 and
4.350 s in LLC96, respectively 83.7% and 87.6% of diagnostic API wall time.
Serial mem_plus and RAM each invoked SNB once per domain (under 1.4 ms each),
with no final public SNB refactors. The other serial circuits and every
eight-thread configuration recorded zero SNB invocations; the parallel SNB
entry point recorded zero calls throughout.

These are actual diagnostic invocation records, rather than constructed
schedule census counts. Instrumentation can influence engine tournaments,
so the results do not guarantee uninstrumented dispatch and support no
competitive timing claim. They narrow SNB's measured role to repeated serial
BJT here, and motivate broader numeric work on mapped and EGraph paths.
No circuit-name selector, production change, or promotion follows from this
profile. The fastest-SPICE goal remains unmet.


### Actual mapped and EGraph invocation attribution (2026-10-10)

The isolated `build/mapped-egraph-invocation-attribution-20261010` package instruments six owner entry points against the frozen integer-certificate baseline. All 28 original full-circuit runs (seven circuits, two cache domains, one/eight threads) qualified with unchanged waveform, successful/failed step, simulation-end and error gates. The controller exited successfully; PID 4000902 was absent before review. Independent review verifies the exact unique grid, reparses both diagnostic log records, and checks all frozen preparation/build/run hashes. Only `kls.c.o` differs among seven archive members. Removing the six timer scopes and their helper reproduces every baseline source/include file byte for byte.

| Circuit | Threads | Timed entry | Calls per domain | Entry seconds / API wall seconds across domains |
|---|---:|---|---:|---:|
| invertbjt | 8 | egraph | 1145 | 72.08–73.08% |
| mem_plus | 1 | mapped | 1201 | 54.22–61.09% |
| mem_plus | 8 | mapped | 1202 | 60.22–67.16% |
| dac | 1 | mapped | 516 | 52.01–53.42% |
| dac | 8 | mapped | 516 | 56.93–57.86% |
| ram2k | 1 | mapped | 5797 | 49.12–55.64% |
| ram2k | 8 | mapped | 5798 | 56.11–64.49% |

These are attribution fractions, not competitive ratios. Outer `parallel` and inner engine timers overlap; rejected trials and early returns are included, and instrumentation can affect engine tournaments. The named `parallel` dispatcher also owns serial calls. The `singleblock` fallback is entered only twice for serial BJT; it is absent for mem_plus, dac and RAM. PTS entries in the eight-thread mapped cases return in microseconds and do not represent their recurring numeric work. MOS13, mux8 and NEURON spend negligible time in the instrumented refactor entries, so this evidence does not assign their remaining work to a particular kernel.

The recurring mapped region, including its BTF and other subordinate paths, is the next general optimization target for mem_plus, dac and RAM; eight-thread BJT remains an EGraph target. No circuit-specific dispatch, thresholds, numeric changes or diagnostic-build promotion follows from this review. The objective of beating both competitors remains unmet. Full entry data, receipts and limitations are retained in `attribution-review.json`, `source-review.json` and `build-review.json`.


### Generic mapped block constant-mode specialization, isolated candidate (2026-10-10)

The mapped/EGraph attribution review motivated an isolated compiler specialization of `kls_parallel_refactor_block`. `build/mapped-block-constant-mode-preparation-20261010/kls.c` has SHA256 `59373416ca5236cc1c1c0e52dbcb451d1081e8f4923b2124a0550deef69d7048`. An always-inline implementation receives constant unchecked/unscaled-map flags only when the existing generic branch would apply; all other states call the same body with runtime policy. No arithmetic, row order, input bounds, pivot rules, recovery behavior, caller eligibility, circuit selector or tuning threshold was deliberately changed. Independent inverse transformation reproduces the complete frozen integer baseline source, and all 384 Boolean admission cases preserve both input branches and the pivot-check condition. This is static equivalence evidence, not a bitwise floating-point or all-engine audit.

Two full native/sanitizer source preparation attempts were interrupted after confirmed live processes stalled in remote vendor file reads. The original `mapped-block-constant-mode-validation-20261010` package dereferenced vendor links and stalled at SCOTCH `tgt/h2.tgt`; the corrected package preserved links but stalled reading frozen KLU `trilinos_klu_memory.c`. Both packages, interruption records and partial copies are retained. Neither launched its validation controller; the four-suite native/O3-ASAN/UBSAN qualification remains pending. This is not a passing full-suite result.

The candidate instead compiled against already-qualified dependency objects in `build/mapped-block-constant-mode-api-qualification-20261010`. All six other archive members remain byte-identical; only `kls.c.o` differs. All 28 original full-circuit API-wrapper runs passed waveform rtol1e-8/atol1e-10, successful/failed step, completion and API/certificate/statistics error gates. The writer PID 4003664 terminated, its exec session exited zero, and review independently verifies the unique grid, diagnostic log records, and frozen preparation/build/run hashes. This establishes that simulation qualification only; wrapper timings are noncompetitive and may affect adaptive engine selection. The mapped function symbol grew from 6,989 to 12,169 bytes under the same native compile flags, leaving instruction-cache pressure as a measurement concern.

The separately linked uninstrumented executable `/tmp/kls-xyce-mapped-block-constant-mode-20261010` has SHA256 `ae14b9f9e2137c1c0a268b9c15486ff24e9098571cf143fa08696f47afa38156`. The frozen `mapped-block-constant-mode-spice-pair-controller-20261010` protocol has 238 observations: 14 fresh KLU references and 224 candidate/baseline runs across every original circuit, both cache domains, one/eight threads and four cyclic balanced passes. Each arm occupies each position twice. The writer PID 4005885 was verified live after metadata creation; its exec session is 10938. No complete summary or timing ratios are claimed here. Failures are retained without filtering or retries. No production promotion is authorized by this intermediate result, and beating both competitors remains unmet; broad/native/sanitizer and full competitive qualification are still required for a successful candidate.


### Generic mapped constant-mode specialization, complete timing decision (2026-10-10)

The frozen `xyce-mapped-block-constant-mode-pair-campaign-20261010` campaign completed all 238 observations, with all 112 candidate runs, 112 baseline runs and 14 fresh KLU references qualified. The writer PID 4005885 was absent, the exec session exited zero, and its terminal receipt matches the summary hash. `mapped-block-constant-mode-spice-pair-review-20261010/review.json` independently verifies the exact unique complete grid, frozen inputs, every reported timer/step log record, verification-netlist hashes and all four aggregate calculations. No retries, rejected trials or outliers were filtered.

| Cache domain | Threads | Candidate / baseline geometric mean |
|---|---:|---:|
| llc32 | 1 | 0.997124983 |
| llc32 | 8 | 1.000800954 |
| llc96 | 1 | 0.999678201 |
| llc96 | 8 | 1.025332497 |

The specialization is rejected. The eight-thread LLC96 aggregate is 2.53% slower, including a 1.187131048 candidate/baseline median on NEURON (18.7% slower, three of four ratios approximately 1.186–1.192). Serial geometric means improve only 0.29% and 0.03%, and LLC32 eight-thread performance is 0.08% slower. RAM medians improve approximately 1.0–2.4% across the four cells; they do not establish a general winning candidate. Function code size increased from 6,989 to 12,169 bytes, but the measured results do not by themselves assign causation to instruction-cache pressure or engine selection.

`mapped-block-constant-mode-performance-decision-20261010/decision.json` records no promotion and no further broad/competitor campaign for this candidate. The isolated source, original interrupted native/sanitizer preparation packages, all API qualification artifacts, timing observations and independent review are retained. No new circuit/shape selector follows from the isolated gains. The retained integer-certificate baseline stays unchanged. This comparison does not measure CKTSO/SubtreeLU superiority; the goal remains unmet. The next general mechanism must reduce storage/update cost beyond cloning the existing mapped block policies.


The postmortem also identifies a limit in the earlier mapped attribution: the outer `kls_mapped_refactor` timer includes lean-single-block, row, BTF and fallback executions, while the six-entry diagnostic did not time the lean-single-block or row entries. The outer dispatcher sets the public MAPPED label on any nonnegative return and can overwrite a subordinate label. The generic-block specialization therefore had a weaker invocation premise than initially stated; its full source-level timing comparison is valid, but the outer timing does not prove that the modified kernel ran repeatedly in the targeted cases. `mapped-block-constant-mode-performance-decision-20261010/attribution-limitation.json` retains the exact source hashes and this correction. Before another kernel change, owner-entry attribution must distinguish the omitted subordinate paths and the mapped BTF body.


### Actual mapped subordinate-path attribution (2026-10-10)

`build/mapped-subkernel-invocation-attribution-20261010` adds owner-entry timers and post-gate markers to the frozen integer-certificate baseline for outer mapped refactor, lean single-block, serial/threaded row entries and the mapped BTF body. Post-gate means that the documented structural/entry validation passed; it does not mean numeric acceptance or completion. There are no worker/inner-loop counters or timers. Removing the five scopes, five marks and reporting helper reproduces every baseline source/include file exactly. Six archive members are byte-identical, with only the diagnostic numeric owner changed.

All 28 original full-circuit observations qualified, preserving waveform rtol1e-8/atol1e-10, step counts, completion and error gates. Writer PID 4021307 was absent, its exec session exited zero, and `attribution-review.json` verifies the exact unique grid, both diagnostic records, post-gate/count/time relationships and all frozen preparation/build/run hashes. The full independent log/source/archive reviews are retained.

| Circuit | Threads | Numeric region | Calls across the two domains | Post-gate seconds / API wall seconds |
|---|---:|---|---:|---:|
| mem_plus | 1 | btf_body | 1201 | 53.92–60.83% |
| mem_plus | 8 | btf_body | 1202 | 60.42–68.57% |
| dac | 1 | btf_body | 515 | 52.84–53.52% |
| dac | 8 | btf_body | 516 | 56.87–57.92% |
| ram2k | 1 | btf_body | 5797 | 49.22–55.43% |
| ram2k | 8 | btf_body | 5798 | 56.68–64.04% |
| rallpack3 | 1 | lean_single | 7717 | 41.28–41.34% |
| rallpack3 | 8 | lean_single | 7719–7720 | 52.14–52.77% |

For mem_plus, dac and RAM, every mapped entry reaches the BTF body in both thread modes and cache domains. This directly resolves the earlier ambiguity and supports the generic block kernel as a recurring numeric target in these instrumented runs. Serial and threaded row-refactor entry counts are zero in all 28 observations. NEURON instead reaches the lean single-block routine thousands of times; its outer mapped entries decline before the post-map gate. BJT has only two outer mapped entries in the serial observations and none in the eight-thread observations, consistent with the separately reviewed serial SNB and eight-thread EGraph routes. MOS13 has three brief mapped/BTF entries per observation, while mux8 has none; this diagnostic does not assign their other work to a particular numeric kernel.

These fractions are attribution, not speedups. Parent/subordinate and API timers overlap, trials and declines remain counted, and instrumentation may alter engine tournaments. No uninstrumented dispatch-equivalence, all-engine correctness or competitor-superiority inference follows. The rejected constant-mode clone remains rejected despite this clearer invocation evidence. Future changes should address BTF block storage/update cost and lean single-block numeric work while retaining general validation, pivot, aliasing, ordering and recovery contracts. No diagnostic source was promoted, and the performance objective remains unmet.


### Producer-published L metadata descriptors: correctness progress (2026-10-10)

`build/producer-l-descriptor-preparation-formatted-20261010` is an isolated
prototype for the generic BTF block owner and lean single-block map32 owner.
Each producer publishes its packed L row pointer and length after completing
that column; a later dependency reads these descriptors only for `j < k`.
The table belongs to one numeric call, carries no values, and never survives
repivoting or a numeric rebuild. A fixed 64-entry stack buffer bounds local
stack use; larger blocks use optional heap storage. Allocation failure and
non-predecessor dependencies retain the original decoder. This resource
limit does not exclude larger inputs or choose an engine by circuit identity.

The static inverse review removes the descriptor mechanics and reproduces
both baseline files exactly. All thirteen post-allocation cleanup sites are
accounted for. No production numeric source has been replaced by this
prototype, and no persistent solver/worker fields or public ABI were added.

`build/producer-l-descriptor-runtime-validation-complete-20261010` reuses the
previously qualified native and optimized ASAN/UBSAN dependency libraries in
Apptainer. It recompiles the numeric owner, all 36 C test executables and the
two executables consumed by Python test gates. The benchmark-support and
benchmark-relink packages preserve those additional builds and ensure each
executable links the candidate numeric archive. Native validation passed all
44 existing CTest gates, including the guarded caller-RHS fixtures. Independent
native archive review confirms that six of seven members remain byte-identical
and only `kls.c.o` changed. The full optimized sanitizer suite is still running;
its eventual result is a separate required gate.

`build/producer-l-descriptor-differential-corrected-preparation-20261010`
completed a direct-owner audit independent of automatic dispatch. For each
build mode (native and optimized ASAN/UBSAN), baseline, descriptor and
forced-allocation-failure variants produced 6,804 snapshots each. Snapshots
include meaningful packed L/U rows and values, diagonal/off-diagonal/scale
storage, numeric workspace, floating-point exception flags and failure/exit
fields. Dense, banded and fragmented patterns cover dimensions
2/3/17/63/64/65/129, three scaling modes, four rounding modes, three numeric
epochs, singular exits and checked pivot rejection. All three variants match
byte-for-byte within each build mode; cross-build byte equality is not claimed.
Each descriptor variant attempted 1,512 optional allocations; the failure
variant forced every attempt to fail and still matched its baseline. Descriptor
read counters prove the new path was exercised. The audit controller exited
zero and is absent. Source inversion, executable/input hashes, snapshot hashes,
coverage counts and logs are retained in `completed-differential-review.json`.
This audit does not cover fused transpose, enabled exception traps, FTZ/DAZ or
concurrent scheduling; the full solver suite remains a separate gate.

The initial differential package failed in the baseline fixture before any
candidate run: it dereferenced the optional solver-owned value buffer. The
corrected fixture uses `prepare_numeric_values` and its returned pointer,
including the serial pass-through case. That failed package is retained.
Earlier validation preparation corrections (target count and notice-file
coverage) also remain retained, with no controller launch for those packages.

The full 28-observation API qualification, uninstrumented link and balanced
238-observation SPICE baseline pair controllers are prepared but have not run.
They retain every original full circuit, both cache domains, both thread counts
and all qualification gates. Competitive CKTSO/SubtreeLU comparison remains
required after a favorable uninstrumented result. There is no performance,
promotion or goal-completion claim from these correctness checks.


Producer-descriptor validation follow-up (2026-10-10): the optimized sanitizer
build and all 44 CTest gates completed successfully. Native and sanitizer suites
therefore passed 88/88 gates. PID 4023921 is absent and exec session 1495 ended
with exit zero. `completed-validation-review.json` verifies all 505 frozen
inputs, both seven-member archives (only `kls.c.o` changed), all 44 distinct test
names in each log, and the candidate numeric object used by both benchmark
support executables in each build. The direct differential snapshot framing
review independently recovered the complete unique 6,804-case grid in each
build mode, with 84 checked-pivot-rejected and 756 singular snapshots and zero
invalid snapshots. Full 28-circuit API qualification has now launched; its
results, uninstrumented competitive performance and promotion remain pending.


### Producer-published L metadata descriptors: rejected performance (2026-10-10)

Full API qualification passed 28/28 observations and its independent review.
The uninstrumented balanced baseline pair campaign then completed 238/238
qualified observations: 112 candidate, 112 baseline and 14 fresh KLU references.
It retained all seven original full circuits, both cache domains, one/eight
threads and four cyclic balanced paired passes. Actual writer PID 4030421 is
absent and exec session 26022 exited zero. The independent review confirms the
exact unique grid, every log timer/step count, frozen input hashes, verification
netlist hashes and all four complete aggregate calculations.

| Cache domain | Threads | Candidate / baseline geometric mean |
|---|---:|---:|
| LLC32 | 1 | 1.012734996225 |
| LLC32 | 8 | 1.024417297686 |
| LLC96 | 1 | 1.012300697274 |
| LLC96 | 8 | 1.015577728992 |

All four cells regress. The serial lean-single NEURON medians regress about
3.4–3.8%; mem_plus, dac and RAM also regress in all cells. All noisy parallel
NEURON observations remain retained. No case/domain subgroup rescues the
prototype and no circuit-specific eligibility is introduced. Additional
per-column stores and optional allocation are plausible costs, not proven
causes; this campaign did not profile those costs separately.

`build/producer-l-descriptor-performance-decision-20261010/decision.json`
pins the completed correctness, differential, API and paired timing reviews.
The prototype is rejected, no numeric source is promoted, and no broader
competitor campaign is justified for it. Passing 88 native/sanitizer tests and
all six differential variants establishes the stated correctness coverage,
not a speed benefit. The retained integer-certificate baseline and the
CKTSO/SubtreeLU performance gap remain unchanged. The full goal is unmet.


### Packed-pointer and workspace-alias code generation: no benefit (2026-10-10)

`build/packed-index-decoder-codegen-probe-20261010` compares the exact retained
packed-index helper with a direct equal-width address expression. GCC 15.2
inside the pinned Apptainer image emits identical instruction streams; the
compiler already removes the proposed alignment arithmetic. The alternative
is not promoted, and no weakened metadata bounds or alternate-width contract
is introduced. This is a code-generation diagnostic, not a numerical audit or
a measured speed claim.

`build/btf-workspace-restrict-codegen-preparation-20261010` builds two complete
numeric owners with the same qualified native recipe. The only proposed source
change declares the generic BTF worker workspace pointer `restrict`. Both
whole object files are byte-identical, including the 6,989-byte BTF function.
Session 1611 exited zero. Declaring ownership this way offers no code-generation
benefit in this build, so no correctness or performance campaign is justified
for the diagnostic and no numeric source is promoted.

The next diagnostic is an integer operation census in the generic BTF and
lean-single owners, with local counters merged once on every invocation exit.
It retains all original full circuits/configurations and waveform/step/error
gates; arithmetic and source loop order invert exactly to the retained
integer-certificate baseline. Instrumentation can affect engine selection and
code/cache layout, so eventual counts describe instrumented paths only. They
will not establish uninstrumented dispatch or competitive superiority, and
will not justify circuit identities, corpus thresholds or exclusion of losing
configurations. The initial preparation assertion miscounted the lean-body
loop substrings; it failed before any controller launch and is retained. The
corrected preparation and its direct numerical/FP/exit differential gate are
separate packages. The overall CKTSO/SubtreeLU goal remains unmet.


### Full-grid mapped-owner work census (2026-10-10)

Native and optimized ASAN/UBSAN direct-owner audits each reproduced all 6,804
baseline factor/FP-flag/exit snapshots byte-for-byte. Independent review checks
the exact unique domain, framing and failure coverage, and the integer counter
conservation identities. The full API/census then qualified all 28 original
full-circuit configurations. Independent review verifies source inversion,
all preparation/build/run hashes, the complete unique grid, log fields,
waveforms, step counts and terminal controllers. Writer PID 4046650 is absent
and exec session 25995 exited zero. Six of seven numeric archive members are
byte-identical to the retained baseline; only the diagnostic numeric owner
changed. The initial API build omitted its generated SCOTCH header and failed
before any circuit observation. That failed package is retained; the separate
ABI-corrected package supplies the exact qualified 64-bit SCOTCH header and
uses byte-identical differential-audited diagnostic numeric sources.

The following ranges include all four cache-domain/thread configurations,
with no filtering. Short means the existing zero-to-three-entry length class;
it defines an observation bin, not a new selector. Update-entry share counts
actual active scalar L updates, excluding batched consumers. The table omits
zero-count owners for readability; full 56 owner/configuration records remain
in `qualification-review.json` and the pinned evidence file. BJT and mux have
zero calls to these two owners in this instrumented census and remain in all
qualification and performance protocols. Their other kernels are not covered
by this census.

| Circuit | Owner | Mean active L length | Short dependency share | Short update-entry share | Batched dependency share |
|---|---|---:|---:|---:|---:|
| invert50_mos1 | btf | 1.89 | 100.00% | 100.00% | 0.00% |
| mem_plus | btf | 14.49 | 67.56% | 11.31% | 0.04% |
| dac | btf | 5.99 | 56.08% | 22.78% | 0.00% |
| rallpack3 | lean_single | 1.00 | 100.00% | 100.00% | 0.00% |
| ram2k | btf | 11.09 | 87.79% | 18.65% | 0.01% |

The recurring BTF cases have many short dependencies but most actual scalar
update entries come from longer columns: short columns account for about
11.3–22.8% of update entries. The exercised lean-single walk uses only
one-entry active L columns. These are integer operation counts, not phase
costs or proof of a wall-time bottleneck. Counts include initial/rejected
trials and instrumentation may alter engine selection or code/cache layout;
they do not establish uninstrumented dispatch or superiority.

`build/mapped-owner-work-census-evidence-20261010/evidence.json` pins all reviews
and all configuration metrics. The next question is whether the original-order
geometry and producer reuse of long L streams warrant a once-per-pattern
segmented representation. Existing wide-scatter/prefetch experiments remain
rejected. No circuit identity, corpus threshold, subgroup rescue or numeric
source promotion follows from these counts. General row order, floating-point
behavior, representation invalidation and scalar fallback remain requirements.
The overall CKTSO/SubtreeLU performance goal remains unmet.


### Original-order L runs: qualified structural opportunity (2026-10-10)

`build/factor-pattern-geometry-checked-preparation-20261010` captures L/U row
patterns after the first successful factor and first successful public
refactor/combined call. It adds only a read-only post-call serializer and API
wrappers; the original numerical functions invert exactly to the retained
baseline. The exact qualified generated SCOTCH header is retained, and six
of seven numeric archive members are unchanged. An initial preparation syntax
check caught a doubled quote before any build/controller launch; that package
is retained. All 28 full-circuit configurations qualify. Independent review
validates all 56 snapshot files, exact framing/EOF, block/index bounds, strict
triangular ranges and per-column index uniqueness, plus the full unique grid,
source/hash/waveform/step/error and terminal gates. Writer PID 4049018 is absent;
exec session 47349 exited zero.

The two snapshots differ in row-pattern geometry for BJT, mem_plus and RAM.
This makes row-order invalidation a required cache boundary even when factor
pointers and column lengths remain unchanged. Pointer-address comparison alone
cannot validate a future run plan. Canonical sorting, repivoting and any numeric
rebuild must retire plans before the stored row order can change. No run plan
is integrated or cached by the production solver at this stage.

`build/forward-l-runs-kernel-preparation-20261010` prototypes forward-run
updates without sorting or changing row order. Scalar gaps retain the exact
original helper; consecutive rows permit ordinary contiguous compiler
vectorization. Plans use only row geometry. Nearest/gradual/masked-SIMD mode
admits the vector path; other modes call the same original scalar entry point.
The isolated native, generic x86-64 and optimized ASAN/UBSAN audits pass
11,592 cases each (34,776 total), including all four rounding modes, preexisting
exception flags, NaN payloads, infinities, subnormals, repeated indices across
runs and guarded workspace boundaries. The ten FP-control guard variants
decline the vector path. `forward-l-runs-trap-audit-preparation-20261010` adds
32 first-fault/partial-workspace trace pairs in each build, across four enabled
trap classes and eight fault positions (96 pairs, 192 child invocations).
All match the original entry point. These are standalone kernel proofs; they
do not qualify a future solver cache, numeric-owner dispatch or thread-pool
exception environment. Sessions 95898 and 39956 exited zero and their
controllers are absent by exact process-argument scan.

`build/forward-l-runs-pattern-micro-preparation-20261010` completes 56
native/sanitized one-frame audits before timing every original 28 captured
post-refactor pattern in native and generic builds. Four balanced AB/BA/BA/AB
pairs per cell give 224 byte/flag-qualified pairs. Independent review verifies
the complete unique grid, every input/log/hash, replay frame count and terminal
gate. Writer PID 4051936 is absent; exec session 13897 exited zero. Whole
replays use fixed normal dyadic values and a nonzero scalar. Plans/allocation,
file IO, factor/pivot/solve work and numeric lifetime management are excluded.
Both modes run on CPU8; domain/thread labels identify captured input patterns,
not measured cache-domain or thread scaling. Original helpers are exported
noinline, so integrated code-generation/cache effects may differ.

The table retains all four captured configurations per circuit. Coverage is
whole four-entry contiguous groups weighted by every stored U dependency,
including possible runtime zeros and batching; it is not executed-work share.
Timing intervals are the minimum/maximum of the four paired-median cells.

| Capture circuit | Longest forward run | Modeled four-entry coverage | Native candidate / baseline | Generic candidate / baseline |
|---|---:|---:|---:|---:|
| invertbjt | 102 | 77.30% | 0.772–0.788 | 0.797–0.808 |
| invert50_mos1 | 2 | 0.00% | 1.060–1.073 | 0.992–1.015 |
| mem_plus | 104 | 76.65% | 0.730–0.753 | 0.776–0.780 |
| mux8 | 3 | 0.00% | 1.035–1.123 | 1.006–1.044 |
| dac | 7 | 7.65% | 1.086–1.099 | 1.045–1.050 |
| rallpack3 | 1 | 0.00% | 1.052–1.064 | 1.033–1.083 |
| ram2k | 126 | 79.28% | 0.721–0.726 | 0.752–0.756 |

Complete 28-cell geometric means of paired medians are 0.918771962088 native
and 0.912600470839 generic. Those component gains do not rescue regressions
on DAC or patterns with little/no run coverage, and do not establish whole
SPICE improvement. BJT's recurring selected kernel is outside the generic
BTF owner; its component result does not authorize a dispatch override.

`build/forward-l-runs-integration-evidence-20261010/evidence.json` preserves
all 28 geometry/work models and required integration gates. Its uncalibrated
integer traffic comparison counts row-index words eliminated versus three
segment-descriptor words per run and two offset words per long dependency.
There are no fitted coefficients or circuit/domain/thread identities. The
model is positive for the high-coverage patterns and negative/zero for the
others, but ignores branch, runtime-zero/batching, cache and setup costs; it
is grounds for a prototype, not proof of profitability. The added `trap-gate.json`
pins the isolated trap proof without weakening the integrated-owner requirement.

The next implementation must prepare plans once per immutable numeric pattern,
retain the original worker on missing/failed/unprofitable plans and unsupported
FP modes, audit every mutation/invalidation boundary, and charge setup and
dispatch in full-grid uninstrumented timings. Allocation/overflow fallbacks,
88 native/optimized-sanitizer existing tests, actual-owner byte/FP/first-fault
differentials, the full 28 API qualification, heldout generic controls and a
favorable complete baseline campaign remain required before another competitor
campaign. No production numeric source is promoted and the full goal is unmet.


### Immutable forward-L runs: integrated owner and lifetime gates (2026-10-10)

The original integration preparation remains frozen. The lifetime-corrected
preparation adds retirement before KLS block repair, predicted-pattern closure
and wholesale numeric replacement, and preserves `size_t` packed allocation
bounds without narrowing them to `UF_long`. The complete-dispatch preparation
connects all five remaining generic-worker callers in the serial/tail and PTS
fringe owners; preparation occurs before the corresponding worker loop or
thread dispatch. Existing engine selection and all other numerical owners
remain unchanged. `forward-l-runs-complete-source-review-20261010/review.json`
verifies that stripping marked additions and restoring dispatch calls recovers
all 30 existing frozen numeric source files exactly, and that the runtime
validation sources are byte-identical to this candidate.

The run threshold follows descriptor size: a descriptor consumes three index
words, so four entries is the smallest run eliminating more index words than
its own descriptor costs. Admission additionally charges two offset words per
stored long dependency and requires strict total word reduction. No coefficient,
circuit identity, dimension, domain or thread count was fitted. Runtime zeros,
batching, branch/cache costs and setup still require end-to-end measurement.

`forward-l-runs-owner-differential-preparation-20261010` completed native and
optimized ASAN/UBSAN baseline/candidate snapshots on the full 6,804-case grid.
Each candidate also ran all three separately forced cache-allocation failures.
All five snapshots per build are byte-identical (174,088,224 bytes each), including
meaningful factors, workspace, floating-point flags and pivot/failure state.
The fixture explicitly exercises the original native-index branch in both
arms and invokes the candidate's real builder and whole-worker dispatcher;
the unchanged lean owner is a control. Native and sanitized normal candidate
runs each record 270 optimized worker calls; all forced failures record zero.
Writer PID 4053230 is absent and session 78938 exited zero. This establishes
the integrated kernel differential, not public scheduling or profitability.

The complete-dispatch candidate passed all 44 native existing tests. Optimized
sanitizer tests, generic x86-64 control, and independent lifecycle/whole-worker
trap reviews are still underway. The first lifecycle/trap fixture's baseline
correctly did not overflow on `1 / DBL_MIN`; that reciprocal is finite. Its
failure and trace are retained. A separately prepared corrected fixture uses
`DBL_TRUE_MIN` and its native lifetime assertions and twenty whole-worker
first-fault traces pass and match the baseline. No failed preparation was
edited or silently rerun. Full 28-cell API qualification and uninstrumented
full-grid timings remain gated; no production numerical source is promoted.


The completed runtime review now confirms all 44 native and 44 optimized
ASAN/UBSAN tests, all support/test links against the candidate archive, every
frozen build input and exactly six unchanged archive members. Writer PID
4053463 is absent and session 3801 exited zero. The generic x86-64 owner audit
also passes the same 6,804-case grid and all three allocation failures; its
five 174,088,224-byte snapshots are identical. Writer PID 4055114 is absent
and session 87737 exited zero. The lifecycle/trap corrected review confirms
all assertions and forty baseline/candidate whole-worker first-fault pairs,
with exact saved MXCSR, signal class, partial LU/Udiag and workspace equality
per build. Every requested exception actually faults. Writer PID 4054818 is
absent and session 61890 exited zero. Public API qualification is the next
gate; these completed proofs establish no end-to-end speed claim.


All original28 API configurations now qualify with independently checked
waveforms, step/end markers, wrapper errors, source/build hashes and six
unchanged archive members. Writer PID4056673 is absent; session55077 exited
zero. These instrumented measurements are not performance evidence.

Before timing, a fresh alias-scope candidate confines the restricted dense
aliases to a block containing only the contiguous update. Scalar gaps and
other runs lie outside that block. The native numeric object compiled from
this physical source is byte-identical to the API-qualified numeric object
(SHA256 d069547571ff0e8d8eb2f4d29039ee6d40d0b58a551be94687f4557ff668b34c).
The read-only Apptainer mount was independently checked: the compiler's exact
virtual helper path contains the new scoped source, not the qualified original.
The existing native executable qualification therefore applies to the same
compiled code image. This bridge does not claim that the complete sanitizer
solver was rebuilt from the scope-only revision.

The separately scoped standalone kernel passes the retained34776 cases.
An expanded fixture adds a scalar gap overlapping both surrounding eight-entry
runs; native, generic and O3ASAN/UBSAN builds pass13248 cases each (39744 total),
including all prior values, NaN payloads, flags, rounding modes and guards.
Sessions62806,4649 and96511 exited zero and their controllers are absent.
No arithmetic, admission coefficients, dispatch or ownership changed. Original
preparations remain frozen. The alias-scoped uninstrumented link and same
all-original238 balanced baseline protocol are prepared; full timing review,
broad performance controls and competitor qualification remain outstanding.

The alias-scoped full238-run balanced campaign is now live (actual writer
PID4059820, exec session55018). All builds and diagnostic controllers finished
before launch. No heavy build/test/profile is authorized alongside its timing.
No partial ratios are reviewed; a terminal session, absent actual writer and
independent full-grid review remain required. The goal is still unmet.


The alias-scoped full238 campaign has completed. Actual writer PID4059820
is absent and exec session55018 exited zero. The independently prepared
review confirms the full unique grid, exact balanced execution order, frozen
inputs, all238 qualifications, timer/log agreement and aggregate reconstruction.
Candidate/baseline geometric means are 1.001095087 (llc32/1), 1.011683684
(llc32/8), 1.003049849 (llc96/1) and 1.003455819 (llc96/8). No configuration
shows a full-grid improvement; this candidate is rejected. Small per-case gains
do not justify subset selection or tuned admission. No production numerical
source is promoted, and broader/competitor timing for this rejected candidate
is not warranted. All original preparations, traces and timing controls remain
retained. Independent review and decision are under
`build/forward-l-runs-alias-scoped-spice-pair-review-preparation-20261010/review.json`
and `build/forward-l-runs-alias-scoped-pair-decision-20261010/decision.json`.
The fastest-versus-CKTSO/SubtreeLU goal remains unmet.


After that rejection and terminal review, a fresh baseline-derived SNB
destination-geometry diagnostic is launched under Apptainer. Actual controller
PID4075173, exec session4246, package
`build/snb-destination-geometry-preparation-20261010`. Existing preparations
remain frozen. A new KLSSNB02 schema adds native block boundaries, actual
row/edge/subcolumn/destination offsets and below-row lists, enabling the
independent reader to check constructor prefix sums and reconstruct every
destination from its global-row mapping. All numerical modules are exactly
the retained baseline, with an appended read-only integer-only serializer.
Public wrappers track the same handle at first successful factor, first
successful refactor/combined call and before destruction; absent SNB states
are valid captures. All original28 waveform/step/error gates remain required.
The conservative integer byte model uses actual 4-byte destinations, three
4-byte run fields, two8-byte offsets, narrow width*subcount and wide
ceil(width/4) decoding consumers. It includes zero-multiplier consumers and
is not an executed-operation or speed claim. No engine/width/partition policy
is changed. Instrumentation may affect adaptive choices. Full execution and
terminal independent review are outstanding; no optimization is promoted.


The SNB geometry controller completed all28 configurations and84 captures
with every API qualification passing. PID4075173 is absent and session4246
exited zero. The prepared reviewer also exited zero (session79182), confirming
all frozen inputs, unique full grid, snapshot framing/prefix sums/partition,
exact global-row destination mapping, source inverse equality and six unchanged
archive members. No nonnegative destination duplicates or monotonicity failures
were observed. Phase1 BJT serial has189441 positive destinations,123090 in
length-at-least4 runs,11417 runs and a10,159,764-byte positive conservative
index-traffic model. Memory serial has71838 positives,4884 run rows and477712
modeled bytes; RAM serial has126717 positives,1440 run rows and94428 modeled
bytes. All phase0 and eight-thread states are absent. Only BJT serial retains
an ordinary state at destruction; other states disappear. All absences remain
in the evidence. These observations are geometry, not executed work or speed.

The prior 20261009 cached-run prototype remains rejected (native GM1.331564).
It used a per-row streak table, including short runs. A potential new mechanism
uses flat three-int32 descriptors for runs longer than their descriptor, retains
original scalar gaps, and places contiguous aliases in a separate restricted
scope. Any admission derives only from checked integer index-byte counts, with
no measured threshold or circuit/owner/width policy changes. Full structural
controls and captured-pattern component evidence must precede integration;
no subset or prior losing-result rescue is authorized. Full reviewed evidence:
`build/snb-destination-geometry-preparation-20261010/geometry-review.json` and
`build/snb-destination-geometry-evidence-20261010/evidence.json`.


A fresh compact-descriptor isolated whole-edge prototype is under
`build/snb-compact-destination-runs-kernel-preparation-20261010`. It replaces
only narrow scatter and wide four-column/one-column scatter with scoped
contiguous updates, preserving all original scalar gaps, edge/row arithmetic
and whole-edge fallback. Admission is strict removed4-byte indices greater
than12-byte descriptors plus16-byte offsets; the common consumer multiplier
cancels per edge. Unsupported MXCSR contexts execute the original whole edge.
This differs from the rejected per-row streak table; no measured threshold,
case/owner/width override or new numeric-value cache is introduced. Production
allocation, ownership and lifetime are not implemented or qualified.

All2724736 fixtures per build pass in native, generic x86-64 and optimized
O3ASAN/UBSAN (8174208 total). This retains the full original narrow grid and
adds real wide-edge controls, repeated contiguous runs, scalar-gap overlap,
simultaneous producer/consumer NaN payloads and70 sticky-flag/trap-mask guard
checks. Whole-panel bytes and FP flags match, with guard pages. Every build
exercises151040 planned/default cases,2120576 unplanned cases and2043720
nondefault cases (categories overlap). Actual writerPID4078359 is absent and
exec session54086 exited zero. Independent source/hash/log/fixture review:
`build/snb-compact-destination-runs-audit-review-20261010/review.json`.
The audit does not prove production lifetime, enabled-trap first-fault equality
or speed. No numerical source is promoted.

The separate frozen performance package
`build/snb-compact-destination-runs-performance-preparation-20261010` launches
all12600 shapes per native/generic build with four balanced paired trials
(201600 observations). All original6480 shapes remain; wide edges and overlap
maps are added. Both compilers complete before any timing. A persistent
container benchmark server removes repeated container startup outside the
kernel timer; every request retains full byte/flag checks, resets, batch clocks
and call overhead. Baseline112-byte and candidate128-byte block layouts are
explicit. No classification or production construction/lifetime cost is timed.
Actual controllerPID4078862, session74885; no heavy jobs are allowed alongside
this isolated timing and no live/partial ratios are reviewed. Complete terminal
independent review, actual-capture evidence and full solver qualification remain
required; the original competitor objective remains unmet.


The compact-descriptor performance controller finished all201600 observations
(12600 groups per build), actualPID4078862 absent and exec session74885 exited
zero. Prepared independent review (session46386 exit0) confirms the full unique
grid, exact balanced order, all original6480 shapes plus every added control,
all per-request byte/flag/post-batch checks, unchanged inputs/binaries and exact
aggregate reconstruction. Native candidate/baseline GM1.090915601 with8401
slower groups; generic GM0.991037932 with7024 slower groups. The native
regression rejects this implementation. Winning map classes or generic-only
results are not selected into a fitted policy. No production integration,
actual-capture speed inference or superiority claim is warranted. Original
source and all prior failures/controls remain retained. Decision:
`build/snb-compact-destination-runs-decision-20261010/decision.json`. The
CKTSO/SubtreeLU fastest objective remains unmet.


After rejecting compact SNB runs, the full28 qualified API diagnostic records
attribute roughly19-38% of serial API time on most circuits to certification
(BJT about2.3-2.4%). This is instrumented attribution, not an uninstrumented
speed claim. A new mathematically guarded proof investigates exact integer
comparison of row-bound fractions, preserving each row's rounded-division
acceptance through the exact midpoint below an odd-significand tolerance,
and dividing only the exact maximum. The existing threshold multiply absorbs
any division inexact flag; power-of-two denominators have exact normal
quotients. Guard4*DBL_EPSILON<=tolerance<=1, odd significand, default masked
nearest/gradual context; unsupported cases must retain the original routine.
No fitted circuit/dimension/thread guard or numerical-value cache is proposed.

The first standalone proof run passed but its audit assembly lacked the
production volatile annotation. It remains retained and is not inherited as
production flag evidence. A fresh corrected exact-form package
`build/componentwise-deferred-division-volatile-proof-preparation-20261010`
passes1734062 tested and1676036 accepted primitive cases in each native,
generic x86-64 and O3ASAN/UBSAN build (5202186 tested,5028108 accepted total).
Every build retains190908 subnormal-denominator cases,265984 power-of-two
denominators,400000 near-cutoff attempts, all prior MXCSR status patterns and
exact rounded maximum agreement. An independent Python Fraction oracle passes
20000 arbitrary IEEE fraction comparisons and20000 midpoint comparisons.
ControllerPID4080251 is absent; session98429 exited zero. Completed hash/log/
scalar-instruction-form review is in that package's `completed-review.json`.
This establishes only the standalone guarded primitive proof. Actual row/range
order, extended fallback, all-zero and unsupported contexts, public certificate
execution, production build/lifetime and full performance remain outstanding.
No production source is changed or promoted. The goal remains unmet.


The deferred-division candidate now wraps the actual complete range routine
while retaining the unchanged scalar row API and original range fallback.
Its cloned budget prefix is arithmetic/order identical; exact midpoint
comparison makes each row's acceptance immediate, so early extended fallback
occurs on the same row. The exact fraction maximum gets one final division.
Tolerance admission is lazy so empty/all-zero ranges perform no extra cast.
Unsupported FP/compiler contexts enter the original routine directly; a
non-admitted tolerance replays the original unchanged range under masked traps
and sticky flags. No plan/worker/public ABI, partition/choice or factor numerical
source changes. Inverse source review is exact, all other modules unchanged.

The compiled actual-row/range differential passes native, generic x86-64,
Intel assembly syntax and separate-contraction O3 builds. Each retains2080
fixtures,29250 range comparisons,28866 direct status/ten-byte-bound/MXCSR/x87
flag comparisons,8656 empty-range checks,4720 unsupported-tolerance checks,
13509 exact-rational row-bound checks and112 real OpenMP comparisons. Total
117000 range and115464 direct comparisons; the scalar row helper remains byte
identical. Both O3ASAN/UBSAN arms produce identical102400-observation files,
including98304 direct context/tolerance/empty-range comparisons. Controller
PID4080796 is absent and session13729 exited zero; all input/output hashes,
source inverse, commands, logs and sanitizer byte files are reviewed in
`build/componentwise-deferred-division-range-preparation-20261010/completed-review.json`.
This is isolated actual component execution, not full production/API proof.

A separately frozen range micro campaign is live under
`build/componentwise-deferred-division-range-micro-preparation-20261010`,
controllerPID4089543, session45603. All source is exactly the range-qualified
candidate; both original and new range arms use the same plan layout and
shared call/clock overhead. Complete144 groups per native/generic build cover
three sizes, sparse/degree16/hub/diagonal/variable-degree/empty structures,
exact/zero/accepted-perturbed/rejected RHS configurations and both orientations
with four balanced paired passes (2304 observations). Every group verifies
status, bounds, MXCSR and x87 flags before timing, and status/bounds after
every batch. Construction excluded; compilers finish before timing. No heavy
job or partial ratio review is allowed during this campaign. Full terminal
review, API/full tests and original SPICE/competitor/broader performance remain
outstanding. No production source is promoted; the original goal is unmet.


The full2304-observation range micro campaign finished with controller
PID4089543 absent and session45603 exit0. Prepared independent review passes
all144 groups per build, unique full grid and exact balanced order, batch
qualification, frozen inputs/executables and aggregate reconstruction. Native
GM0.927302751 (71/144 groups slower), generic GM0.919591044 (77/144 slower).
All zero/rejected/empty/variable-degree and both-orientation groups remain.
These are component-only gains with setup excluded and shared call/clock
overhead; they do not prove public solver speed or generalization. The
mathematical guard and all selectors remain unchanged. Full native/O3san
solver tests, public API, complete original SPICE comparison, broader controls
and competitive timing are required before promotion. No winning subset or
fitted circuit/degree/size threshold is selected. Review and next gate:
`build/componentwise-deferred-division-range-micro-preparation-20261010/review.json`
and `build/componentwise-deferred-division-range-evidence-20261010/evidence.json`.
The overall competitor objective remains unmet.


Full solver qualification of the exact deferred-division range candidate is
now running in Apptainer: new package
`build/componentwise-deferred-division-runtime-validation-preparation-20261010`.
Only `kls_componentwise.c.o` is recompiled in each retained native/O3
ASAN-UBSAN archive; the six other archive members stay unchanged. All test
and support executables are linked after archive publication. The baseline
componentwise per-source compiler rule was checked and has no additional
overrides beyond its recorded C_DEFINES/C_INCLUDES/C_FLAGS. Prepared review
requires all44 native and all44 optimized sanitizer tests, exact range
candidate sources/headers, original unique CTest grids, archive-member
identity and input/build hashes after actual controller absence and terminal
exec. Controller PID4089861, exec session54959; no result inferred while live.

The next original28 public Xyce API package and independent review are
prepared, not launched:
`build/componentwise-deferred-division-api-qualification-preparation-20261010`
and `build/componentwise-deferred-division-api-review-preparation-20261010`.
It requires the completed full88 review before compilation, retains six
original Xyce archive members, and prepares separate instrumented and plain
executables with frozen original Xyce objects/dependencies. Every original
circuit/domain/thread case and waveform/step gate remains. Instrumented API
timings provide no performance claim. Full original SPICE, broader controls
and competitor checks remain required. No production source promotion; the
overall goal remains unmet.

Native full44 CTest execution passed against the deferred-division archive.
Controller PID4089861 and exec54959 remain live during the optimized
ASAN-UBSAN build; full88 terminal review remains pending. The next full238
original balanced baseline protocol and independent grid/order/log/hash/ratio
review are prepared in new packages
`build/componentwise-deferred-division-spice-pair-preparation-20261010` and
`build/componentwise-deferred-division-spice-pair-review-preparation-20261010`.
They require completed full88 and original28 API reviews, actual writer
absence/terminal exec, unchanged frozen inputs and absence of conflicting
heavy jobs before timing. All seven circuits, two domains, both thread counts,
four balanced passes and14 fresh KLU references remain; all failures retained.
Neither API nor performance has launched. Native passing alone is not the
full qualification gate. Production baseline retained; objective unmet.


The exact deferred-division candidate passes the full solver runtime gate:
44/44 native and44/44 optimized ASAN-UBSAN tests. Runtime controller
PID4089861 is absent and session54959 exited0. Prepared independent review
passes source/header equality with the qualified range candidate, all frozen
input/build hashes, both original unique44-test grids, candidate test/support
links and seven-member archive checks. Only `kls_componentwise.c.o` differs;
all six other members remain byte-identical in each mode. Evidence:
`build/componentwise-deferred-division-runtime-validation-preparation-20261010/completed-review.json`
and `build/componentwise-deferred-division-full-runtime-evidence-20261010/evidence.json`.

Original28 public Xyce API qualification is now launched with actual writer
PID4092854, exec43595. Its builder requires the completed full88 review,
retains frozen Xyce objects/dependencies and six original numerical archive
members, and builds profile/plain executables in Apptainer. These runs are
instrumented correctness qualification only. API terminal review and full
uninstrumented/broader/competitor gates remain pending; no promotion.

A read-only production integration audit confirms30 original source/header
files equal the retained integer baseline. The only difference is the
working tree's pre-budget `src/kls_componentwise_parallel.inc`. Both integer
budget and deferred division remain isolated, unpromoted changes. A later
promotion must carry the complete qualified delta while preserving dirty
worktree edits. Audit:
`build/componentwise-deferred-division-production-baseline-audit-20261010/audit.json`
and its `interpretation.json`. The overall goal remains unmet.


Original28 public Xyce API qualification finishes successfully for the exact
deferred-division candidate. Writer PID4092854 is absent; session43595
exited0. Prepared independent review verifies the complete unique original
7-circuit/two-domain/one-and-eight-thread grid, waveform/step/end/error/stats
gates, source and frozen input/output hashes, seven archive members and
exactly six unchanged baseline members. Only `kls_componentwise.c.o` differs.
API correctness evidence:
`build/componentwise-deferred-division-api-qualification-preparation-20261010/qualification-review.json`.
No performance inference is drawn from its instrumentation.

The new full238 original balanced uninstrumented baseline campaign is now
launched through the prepared wrapper. Actual writer PID4095044, exec90600.
All seven cases, both domains, one/eight threads, four balanced paired passes
and14 fresh KLU references remain; all failures retained. Actual API writer
absence/terminal review, full88 validation, frozen inputs, public-call timing
boundary and absence of conflicting heavy jobs are preflight requirements.
No heavy build/test/profile or partial ratio review while timing is live.
Terminal exec and writer absence precede the prepared independent full-grid/
order/log/hash/ratio review. Broad controls and competitive qualification
remain outstanding. No source promotion; the goal is still unmet.


The complete deferred-division full238 original balanced baseline campaign
finishes: controller PID4095044 absent, exec90600 exit0. Prepared independent
review passes full unique grid and balanced order, all238 correctness gates,
frozen hashes and independently reconstructed log timers/ratios. Nevertheless
all four full configurations regress; this candidate is rejected. No
favorable subset, circuit/size/degree/thread/build policy is selected. All
original evidence remains frozen. No broader/competitor timing or production
promotion is warranted for this rejected candidate. Review and decision:
`build/componentwise-deferred-division-spice-pair-review-preparation-20261010/review.json`
and `build/componentwise-deferred-division-spice-decision-20261010/decision.json`.
The retained integer baseline and working pre-budget row source stay intact.
The isolated range micro gain did not predict complete solver performance.
The overall fastest-solver objective remains unmet.
Full candidate/baseline GMs: llc32 t1 1.063463860, llc32 t8 1.044713741, llc96 t1 1.064235845, llc96 t8 1.041815504.


Fresh general cutoff-cache candidate follows source-level invariant study and
native assembly audit, independently of any favorable circuit subset. The
parent module executes admitted-tolerance comparisons and midpoint decoding
in the row loop (offsets14fa/1510/151a and1523/152d/1530). New per-range lazy
stack state reuses only that admitted cutoff. Every original per-row
tolerance cast, threshold arithmetic, immediate acceptance/fallback/replay,
scalar row interface, maximum, worker partition/choice and public/plan ABI
remain. Only three private files differ from the retained rejected parent;
inverse patch review recovers its exact source. No matrix values or state are
cached across calls; no circuit/size/degree/thread/build selector is added.

The fresh actual compiled range passes the complete original native/generic/
Intel-syntax/separate-contraction audit:117000 range comparisons,115464 direct
status/bounds/full-MXCSR/x87 exception comparisons, and102400 observations per
optimized leak-enabled ASAN-UBSAN arm, byte-identical. Unsupported tolerance
and FP contexts, zero/empty/rejected rows, independent rational bounds, scalar
identity and all partitions remain. Actual writer PID4110613 absent,
exec13473 exit0; new independent source/hash/log/statistics/sanitizer review:
`build/componentwise-deferred-division-cutoff-cache-range-preparation-20261010/completed-review.json`.

Both fresh complete2304-observation structural timing controls pass terminal
full-grid/order/input/binary reviews. Original in-module micro: native
GM0.905641593 (76/144 slower), generic GM0.888360947 (77/144 slower), writer
4119143 absent, exec38952 exit0. Additional opaque separately compiled non-LTO
shared-module control retains the exact full144-group input grid and balanced
order; frozen shared libraries and executables checked. Native GM0.904592207
(79/144 slower), generic GM0.889637884 (74/144 slower), writer4119389 absent,
exec74972 exit0. No subset rescues or full solver speed claim. Prepared
controllers compile both native/generic artifacts before timing; all known
jobs are terminal. Reviews:
`build/componentwise-deferred-division-cutoff-cache-range-micro-preparation-20261010/review.json`
and `build/componentwise-deferred-division-cutoff-cache-opaque-range-micro-preparation-20261010/review.json`.

The fresh cutoff-cache version still requires its own full native/O3san
solver tests, original28 API checks and full original238 uninstrumented
comparison. The previous parent's full solver/API evidence is not inherited
as proof for this version. Broad and competitor gates remain required before
promotion. Current evidence and exact next source:
`build/componentwise-deferred-division-cutoff-cache-evidence-20261010/evidence.json`.
Production source is untouched and the overall fastest-solver goal is unmet.


Fresh cutoff-cache full solver qualification now runs in Apptainer through
`build/componentwise-deferred-division-cutoff-cache-runtime-validation-preparation-20261010`: actual controller PID4119920, exec25201. Native44 and
optimized ASAN-UBSAN44 tests are freshly compiled and executed against the
exact range-qualified cutoff-cache source. Only the componentwise archive
member is recompiled; six other members retained. Every test and support tool
links after archive publication. New independent runtime review checks both
full44 grids, source/header equality, frozen hashes and archive-member
identity after actual writer absence and terminal exec. The parent's full
runtime/API evidence is not inherited.

New original28 API and original238 balanced comparison/launch/review packages
are prepared but not launched. Paths use the
`componentwise-deferred-division-cutoff-cache-` prefix. API build requires the
fresh completed full88 review; timing requires the fresh completed original28
API review, unchanged inputs, the retained public-call timing boundary and
absence of conflicting heavy jobs. All seven circuits, both CPU domains,
thread counts, four balanced passes, fresh KLU references and all failures
remain. No production promotion or speed claim; the goal remains unmet.


Fresh cutoff-cache full88 qualification completed successfully on 2026-10-10.
Actual controller PID4119920 is absent and exec25201 terminated with exit0.
Independent review verifies all44 native and all44 optimized ASAN-UBSAN tests,
the exact range-qualified source/headers, every frozen input, and only the
componentwise archive member changing. Native archive SHA256
4717c31a132a06119137e5c5855359d240cf79bbbb8cdc6598b685137b61a8e8;
sanitized archive SHA256
493cbc8f91934a2d59eff676876834358926ed6399178b5018270b6b3fac279c.
Review: `build/componentwise-deferred-division-cutoff-cache-runtime-validation-preparation-20261010/completed-review.json`.

The original28 API qualification is now launched through the frozen
`build/componentwise-deferred-division-cutoff-cache-api-launch-preparation-20261010/launch.py`.
Actual controller PID4122690, exec90030. The launch gate checks the freshfull88
review and terminal runtime receipts. API source and public wrappers remain
the prepared original grid; no timing inference from instrumentation.
The original238 uninstrumented timing is not launched and production source
remains unpromoted. The full CKTSO/SubtreeLU fastest-SPICE objective remains
unmet. An optional unimplemented exact rounding-bin study is retained under
`build/componentwise-rounded-maximum-study-20261010`; it is not this candidate
and does not replace any qualification gate.


Fresh cutoff-cache original28 API qualification completed successfully. Actual
controller PID4122690 is absent and exec90030 terminated with exit0. Independent
review checks all28 unique original configurations, waveform/step/end gates,
profile logs, exact candidate source, frozen build/run hashes and only the
componentwise archive member changing. Review:
`build/componentwise-deferred-division-cutoff-cache-api-qualification-preparation-20261010/qualification-review.json`.

The prepared original238 uninstrumented balanced campaign is now launched:
`build/xyce-componentwise-deferred-division-cutoff-cache-pair-campaign-20261010`.
Actual controller PID4124819, exec11382. It retains all seven original cases,
both domains, one/eight threads, four balanced passes per candidate/baseline
configuration and14 fresh KLU references. No partial ratios or filtering.
No heavy qualification runs coexist with timing. Completed independent review
requires actual controller absence plus both controller and exec terminal
receipts. No promotion or fastest-solver claim follows from the API gate;
broad and competitor verification remain required.


Fresh cutoff-cache full original238 campaign completed at 2026-10-10
19:35:42 UTC. Actual controller PID4124819 is absent, exec11382 terminated
with exit0, and independent review verifies complete unique original grid,
exact balanced order, frozen hashes, solver timers and all238 correctness
gates. Candidate/retained-baseline geometric means: llc32t1 1.04627243499,
llc32t8 1.01337339610, llc96t1 1.04984335514, llc96t8 1.02915865121.
All four regress: the cutoff-cache candidate is rejected, unpromoted, with
no selected-subset or circuit/size/degree/thread/build selector rescue.
The positive structural micro results did not establish full solver speed.
Review: `build/componentwise-deferred-division-cutoff-cache-spice-pair-review-preparation-20261010/review.json`.
Decision: `build/componentwise-deferred-division-cutoff-cache-spice-decision-20261010/decision.json`.
No broader/competitor candidate campaign is warranted for this version.
The general rounded-maximum prototype and fresh range audit are prepared
as separate packages; all correctness, FP-status, integration, API and full
performance gates remain fresh requirements. The full fastest-SPICE goal
against CKTSO/SubtreeLU remains unmet.


The isolated general rounded-maximum prototype completed the fresh existing
range qualification in Apptainer: native/generic/Intel-assembly/separate-
contraction O3 modes each passed29250 range and28866 direct status/bound/full
MXCSR/x87 flag comparisons. Both O3 ASAN/UBSAN arms passed102400 observations
with leak detection enabled and byte-identical complete snapshots. Actual
controller PID4140307 is absent, exec40938 terminalexit0. Independent review
verifies all frozen hashes, exact prototype source, unchanged public headers
and other sources, two changed private files and the exact inverse range
patch to the cutoff-cache parent. Review:
`build/componentwise-rounded-maximum-range-audit-preparation-20261010/completed-review.json`.
This covers the existing fixture grid only: fresh ascending/descending/same-
rounded-bin/boundary controls and compiled primitive oracle remain required.
No fullsolver/API evidence is inherited and no performance claim follows.
The rejected cutoff-cache campaign and its originals are retained. No
production promotion; the complete CKTSO/SubtreeLU objective remains unmet.


The general rounded-maximum prototype passed independent fresh boundary/order
qualification in four O3 build modes:514688 compiled exact rational checks,
170112 additional actual range status/bound/full MXCSR/x87 flags comparisons,
400000 optimized ASAN/UBSAN primitive comparisons. Corrected fixture package
retains the original origin1 cases and broadens the coefficient-origin grid
to ensure distinct same-bin rows; the first failed coverage harness remains
unaltered. Review:
`build/componentwise-rounded-maximum-boundary-order-corrected-preparation-20261010/completed-review.json`.

Original full144-group structural micro passed2304 observations, reviewed
only after actual controller PID4149546 absence and exec18379 terminalexit0.
Native candidate/baseline GM0.917173885074 (71/144 slower), generic
0.949622590465 (67/144 slower). Full additional42-group row-order grid passed
672 observations: all ascending, descending, equal, distinct-same-bin,
shuffled, zero-interleaved and final-rejected groups retained. Actual
controller PID4149800 absent, exec55188 terminalexit0; independent complete
grid/order/hash review. Native GM1.17105751863 (41/42 slower), generic
1.13817128164 (39/42 slower). The general order controls materially regress;
this version is rejected with no pooled-ratio, selected-order or circuit/
size/degree/thread/compiler selector rescue. Fresh fullsolver/API/SPICE
campaign is not warranted for this rejected version. No production
promotion; fastest-SPICE versus CKTSO/SubtreeLU remains unmet. Decision:
`build/componentwise-rounded-maximum-decision-20261010/decision.json`.


The isolated normal-product maximum screen preserved the original budget scalar
helper and range body and passed fresh qualification in Apptainer across native,
generic, Intel assembly syntax and separate contraction O3 modes. Independent
source/result review verified 117000 original range comparisons, 266880 additional
range status/bound/full MXCSR/x87 comparisons, 2560000 paired primitive FP-status
checks at product/exponent boundaries, and byte-identical optimized ASAN/UBSAN
observations (102400 per arm with leak detection). Mixed magnitude prefixes and
ascending/descending/equal/shuffled/zero/final-rejection/distinct-same-bin controls
were retained. Review:
`build/componentwise-normal-product-maximum-screen-qualification-source-review-preparation-20261010/review.json`.
Two failed reviewer scripts are retained; their failures were reconstruction
newline/syntax mistakes, corrected in a separate reviewer package.

All four complete micro grids were independently reviewed after actual controller
absence and terminal exec sessions. Original 144-group structural grid: native
GM 0.903930054013, generic 0.914042238020. Full 42-group order grid: native
0.897808507300, generic 0.882981281980. Opaque separately compiled shared-module
144-group grid: native 0.885898206106, generic 0.898238479001. Opaque full
42-group order control: native 1.084033007388 (31/42 slower), generic
1.099368281067 (32/42 slower). All 5952 observations are retained. The opaque
order control materially regresses in both builds, so this general version is
rejected. No selected fixture, circuit, dimension, degree, worker or compiler
rescue and no pooled-ratio promotion. Fresh fullsolver/API/SPICE work is not
warranted for this version. No production numerical source promotion; the full
fastest-SPICE objective against CKTSO/SubtreeLU remains unmet. Decision:
`build/componentwise-normal-product-maximum-screen-decision-20261010/decision.json`.


The general strict-product screen candidate caches the previously accepted
rounded quotient at maximum updates and uses n<RN(r*d), equivalent to comparing
against the product predecessor for binary64 n. The original normal-product
proof guards remain unchanged, and no selectors were added. Fresh Apptainer
qualification passed four O3 modes: 117000 original range comparisons, 266880
additional actual range/full FP-status comparisons, 2652160 paired primitive
FP-status checks, including exact products around DBL_MIN and all64 prior SSE
sticky flags, plus byte-identical optimized ASAN/UBSAN observations (102400 per
arm with leak detection). Source inverse review confirms only two private
files differ from the rejected parent, with the scalar helper and original
range retained. Review:
`build/componentwise-normal-product-strict-screen-qualification-review-preparation-20261010/review.json`.
The failed initial preparation controller performed no build because its receipt
was missing. Its package is retained. The first completed boundary run lacked a
recorded controller PID; the same complete proof was repeated in a separate
package with an automatic PID receipt, then independently reviewed after
actual PID absence and terminal exec sessions. No numerical-source correction
or timing selection occurred in these preparation/replay corrections.

All four complete timing grids (5952 observations) passed qualification and
independent terminal/grid/order/hash review. Original structural grid native
GM0.875865009870, generic0.894945269304; in-TU order grid native0.850868267456,
generic0.842186366062 (0/42 slower in each). Opaque structural grid native
0.887807789957, generic0.877137442892. Opaque full order grid native
1.033278266898, generic1.026928335026 (30/42 slower in each). The general
separate-module order control regresses, so the strict screen is rejected.
No selected circuit, size, degree, fixture, worker or compiler rescue; no
pooled-ratio promotion or new fullsolver/API/SPICE campaign for this version.
No production numerical-source promotion. The fastest-SPICE CKTSO/SubtreeLU
objective remains unmet. Decision:
`build/componentwise-normal-product-strict-screen-decision-20261010/decision.json`.
A separate unimplemented general lazy tolerance conversion study is retained at
`build/componentwise-lazy-tolerance-conversion-study-20261010/proposal.json`;
its correctness and performance remain unproven.


The initial lazy range-local tolerance conversion candidate passed the original
four O3 range/sanitizer gates and extended-tolerance/order controls. All four
complete timing grids improved: original structural nativeGM0.788021416048,
generic0.809850700595; in-TU order0.672477721177/0.685104841619; opaque
structural0.787810437867/0.811173687743; opaque order0.817332231300/
0.835679244473 (0/42 slower in each order build). All5952 observations remain
independently reviewed and retained. These gains do not override correctness.

A fresh lower-optimization zero-denominator proof exposed an O0 x87-status
mismatch: a row with magnitude3*eta and error2*eta has a zero downward budget
denominator. With raw extended tolerance2^-2000, result/bound/MXCSR match, but
original x87 sticky bits50 versus candidate34 omit underflow16. The initial
version is rejected for correctness; no compiler selector or timing rescue.
Failure and original prepared/executed packages are retained. Decision:
`build/componentwise-lazy-tolerance-conversion-decision-20261010/decision.json`.

A separate corrected candidate requests unchanged whole-range replay when its
budget denominator is nonpositive. This restores each compiled original cast
placement under the same masked FP environment without a compiler selector.
Fresh four O3 original range qualification passed117000 comparisons and both
optimized ASAN/UBSAN arms passed102400 byte-identical observations with leak
detection. Fresh full order/mixed/raw-tolerance/zero-denominator proof passed
all10 compiler modes: native/generic/Intel/separate O3, native/generic O0,
native O1/O2/Og/Os. Total3252720 additional actual status/bound/full MXCSR/x87
comparisons, including57 raw extended tolerance representations,14 prefix/
zero/magnitude fixtures, and all64 prior status patterns in each FP unit.
Independent frozen-source/output/terminal review:
`build/componentwise-lazy-tolerance-conversion-replay-qualification-review-preparation-20261010/review.json`.
Corrected source timing, fullsolver/API and full SPICE/competitor evidence
remain fresh requirements; none is inherited from the failed version. No
production numerical promotion. The fastest-SPICE objective remains unmet.


The corrected lazy tolerance conversion/replay source passed all four fresh
complete micro grids (5952 observations), independently reviewed after actual
controller absence and terminal exec sessions. Original structural144 groups
per build: nativeGM0.824914223884, generic0.843363698172. Full in-TU order42
per build:0.696457316205/0.706011662024, no slower groups in either build.
Opaque separate-module structural144:0.820015695831/0.844287072626; opaque
order42:0.843758001260/0.857745851542, no slower groups in either build.
All balanced pairs, original fixtures, rejection/zero controls and groups
remain included. Exact corrected source advances to fresh full runtime/API
and original SPICE qualification, with broader/competitor gates still required.
No production promotion and no fastest-SPICE claim. Evidence:
`build/componentwise-lazy-tolerance-conversion-replay-micro-evidence-20261010/evidence.json`.


Fresh full-runtime preparation initially omitted two SPRAL audit fixtures;
43 native tests passed and the required license-boundary gate failed. No
sanitizer pass or full88 review was inferred. Original process exited0 while
its explicit ctest status was8/all_passed false; originals and failure
receipts are retained. A new complete-source package adds only the exact
SPRAL LICENCE and scaling.f90 audit fixtures, preserving the qualified KLS
src/include bytes. Its controller now exits nonzero on a test failure. All44
native tests passed; native archiveSHA0175397995b2d81b892f8ce51528fbd1af43c7b49e5c57d5bf4bbc4ba8306a03
matches the initial native build. Optimized ASAN/UBSAN build remains active,
controllerPID1996/exec34784, verified from actual process and running handle.
No full88 review or API/performance launch until its actual terminal state
and independent source/archive/test review. Current package:
`build/componentwise-lazy-tolerance-conversion-replay-runtime-complete-source-preparation-20261010`.
The frozen complete28 API grid is prepared and gated on this full88 review;
no production promotion or SPICE performance claim. Goal remains unmet.


Corrected lazy tolerance replay, complete runtime (2026-10-10): controller 1996 is absent and exec session 34784 returned zero. All 44 native and 44 O3 ASAN/UBSAN tests passed. The original complete-source reviewer failed on three regenerated Python bytecode caches; it is retained unchanged. A new frozen derived-cache reviewer verified every other input hash and checked all three cache headers and code objects against freshly compiled unchanged pinned Python sources in the same Apptainer runtime. All numerical, test, archive, and terminal checks passed; completed-review.json records the three prepared/regenerated hashes. No numerical source change or promotion. The frozen complete-runtime API launcher subsequently started PID 5135 / session 16201 for the original full 28 configurations; qualification and full SPICE performance remain pending.


Corrected lazy tolerance replay, public API and SPICE launch (2026-10-10): actual API controller 5135 is absent; session 16201 returned zero. All original 28 configurations qualified, and the frozen complete-runtime API reviewer independently verified the complete grid, source/archive hashes, six unchanged KLS members, parsed public profile logs, step counts, finite waveforms, upstream rtol1e-8/atol1e-10 verification, completion, and no solve/certificate/stats errors. The new frozen full238 SPICE pair controller launched as actual PID7229 / exec session7578, with original seven complete circuits, both domains, threads1/8, four balanced passes, baseline/candidate and14 fresh KLU references. Source, launcher and reviewer packages are frozen; no subsets, retries, selectors, partial ratios or promotion. Performance and competitive fastest-SPICE goal remain pending.


SPICE preflight correction (2026-10-10): initial PID7229/session7578 exited1 before campaign directory or any timing, because controller expected a misnamed executable rather than the exact path in API build metadata. The failed frozen source/launcher/reviewer and terminal receipts are retained. A new complete-runtime SPICE pair package pins the actual qualified /tmp/kls-componentwise-lazy-tolerance-replay-complete-runtime-api-20261010-uninstrumented executable and unchanged numerical source and full238 protocol. Actual corrected controller PID7241/session58094 is live with frozen campaign metadata created, all88 and28 prerequisite reviews verified; no timing analysis or promotion.


Corrected lazy tolerance replay, complete SPICE pair (2026-10-10): actual controller7241 is absent; exec session58094 returned0 at21:40:49UTC. Frozen independent reviewer verified all238 unique observations in the original balanced order, all waveform/step/completion/error gates, input hashes, timer logs and four aggregates. All238 qualified. Candidate/budget-baseline geometric means are0.9853004337096772 (llc32/t1),0.9889797672891878 (llc32/t8),0.9858029257598009 (llc96/t1),0.9975940766338265 (llc96/t8). All individual ratios and slowdowns retained. Decision: advance exact general source to full broader matrix and fresh competitive SPICE qualification; no promotion or CKTSO/SubtreeLU fastest claim.


Corrected lazy tolerance replay, full broader matrix launch (2026-10-10): new isolated broad baseline/candidate links use one unchanged benchmark driver, identical dependency archives and six byte-identical retained budget-baseline KLS members. Only the exact full88-qualified native componentwise object differs. Frozen independent build reviewer passed. Original110 manifest matrices, two8-thread domains,100-system lifecycle,1e-8 per-refactor residual,900-second timeout and4 balanced pairs per matrix produce1760 observations; every failure/timeout retained, no filtered aggregate or rescue. Actual broad controller22906 / exec session85302 is live. Full broad and competitive qualification remain pending; no promotion.


Corrected lazy tolerance replay, competitive SPICE preparation (2026-10-10): frozen source/launcher/independent-review packages preserve the original full462 observations with baseline/candidate/CKTSO/SubtreeLU arms and14 fresh KLU references. Executable identity is exactly the terminal-qualified API metadata; CKTSO/SubtreeLU shared libraries are pinned. Broad full1760 terminal independent review without candidate coverage regressions is mandatory before launch. All competitor failures retained, no aggregate claim if any matrix/circuit coverage is incomplete. No competitive timing launched while broad controller22906/session85302 is live.

## Production certification default: paper review, 2026-10-11

The SubtreeLU paper separates known-solution accuracy validation from repeated
factorization/refactorization timing. The CKTSO paper times fast factorization
with pivot checks, which are distinct from solution residual certification.
CKTSO’s public demo explicitly treats its caller residual calculation as a
result check unnecessary in circuit simulation. SubtreeLU’s demo reports
solution error outside the repeated timed loop. The Xyce third-party adapter
calls both solve APIs directly; its residual reporting is debug-only. This
supports evaluating application-managed accuracy as a production default, but
does not prove that proprietary solver interiors contain no validation.

No production default changes in this review. Existing qualification notes
record MOS13 DC convergence failure with certification disabled at both thread
counts; their referenced temporary raw logs are absent and must be regenerated
before a default decision. Also, the current Xyce adapter checkout lacks the
documented accuracy-policy option, so a fresh candidate must prove the policy
selected by its exact linked adapter. Preserve independent benchmark accuracy
checks, retain every failure, use the full original SPICE grid, and avoid
circuit-specific policy selection. The live frozen campaign remains unchanged.
Source hashes and the bounded evidence review are retained in
`build/production-certification-paper-review-20261011/review.json`.

A separate full28 application-managed qualification package is now prepared at
`build/application-managed-complete-spice-qualification-preparation-20261011`.
Its generic creation wrapper sets and queries the public accuracy policy; the
plain candidate carries no diagnostic telemetry. The diagnostic candidate
requires zero componentwise certificate calls and successful policy selection
for every observed creation. It reuses the exact qualified archive without
changing numerical kernels, retains all seven circuits and both CPU domains
at one/eight threads, and preserves external waveform, step, end, and error
acceptance. The launch/build gates require actual original broad PID absence
and terminal exec exit zero; no heavy work has launched. Python syntax and all
pinned preparation hashes passed. The independent post-terminal reviewer
reparses the complete grid and reruns external waveform verification, retaining
all failures. C compilation, runtime qualification, and performance remain
unproven until those gates execute.


### Completed broad review and application-managed qualification (2026-10-11)

The original lazy-tolerance whole-range replay broad campaign completed all
1,760 observations; controller 22906 was absent and exec session 85302 exited 0.
The frozen independent reviewer reparsed the full grid and logs, checked inputs,
and found no candidate coverage regressions. All 80 invalid observations remain
included, so the broad overall performance aggregate is null.
Evidence: `build/componentwise-lazy-tolerance-conversion-replay-broad-review-preparation-20261010/review.json`.

The separate application-managed policy qualification completed all 28 original
SPICE cases (seven circuits, both CPU domains, threads 1/8). Controller 201786 was
absent and session 24574 exited 0. Independent review confirmed all cases qualified,
external waveform tolerances 1e-8/1e-10 and step/completion gates unchanged,
creation-time policy selection 2 verified for every instance, and zero certificate
calls. This fresh result includes MOS13 at both thread counts and CPU domains;
it does not erase the previously documented opt-out failure. It supplies SPICE
application acceptance evidence, not a certificate for every linear answer.
Evidence: `build/application-managed-complete-spice-qualification-preparation-20261011/qualification-review.json`.

A separate fixed full462 competitive campaign was prepared and launched as
controller 203776 / exec session 29644. It uses the uninstrumented generic
application-managed selector and unchanged qualified archive against certified
budget KLS, CKTSO, and SubtreeLU. All original circuits, four balanced passes,
two CPU domains, thread counts, fourteen fresh KLU references, and external
accuracy gates are retained. No retries, exclusions, overlapping heavy work,
or partial performance claims are allowed. The prepared reviewer independently
rechecks raw waveforms and qualification as well as the full unique order/grid,
logs, inputs, and aggregate calculations after actual process absence and terminal
session exit. Production defaults and numerical source promotion remain unchanged.
Evidence: `build/application-managed-competitive-spice-preparation-20261011/preparation-receipt.json`.

During the live competitive run, the 2026-10-11 09:13 UTC snapshot retained
55 observations, including 14 invalid CKTSO observations. All fourteen logs
reported `CKTSO linear solve failed (-99)`. The public CKTSO header maps -99
to a license error; this identifies the documented error category, not its
underlying cause or an expiry date. No license file was inspected. These failed
runs cannot support a CKTSO performance comparison. The controller continues
with the original order, timeout, and acceptance gates; failures remain in the
complete review. Evidence and source/log hashes:
`build/application-managed-competitive-ck-error-diagnostic-20261011/snapshot-20261011T091308.json`.

A separate general-policy build is prepared, but has not executed:
`build/application-managed-general-policy-build-preparation-20261011`.
Its required entry is
`build/application-managed-general-policy-build-launch-preparation-20261011/launch.py`,
which pins the public header and prepared sources. Both launch and build require
actual competitive controller 203776 absence and exact session 29644 terminal
exit zero; the build additionally requires independent full462 review and no
invalid KLS/KLU observations or overlapping heavy jobs. It preserves the exact
qualified broad archive and driver, copies the certified broad candidate as its
baseline, and adds only the generic creation-time application-managed selector
to the new plain and diagnostic binaries. Syntax and preparation hashes passed;
compilation, independent binary review, and general runtime validation remain
unproven. The existing broad driver independently verifies every refactor and
requires finite final/maximum residuals no greater than 1e-8. Any later general
validation must preserve all 110 matrices, both domains, 100 systems, four
balanced pairs, original owner selection, timeouts, and invalid observations.
This preparation does not change the production default or establish speed.

The separate post-build reviewer is prepared at
`build/application-managed-general-policy-build-review-preparation-20261011/review.py`.
It requires actual build PID absence and the exact terminal exec receipt, checks
every prepared source/output hash, compares all seven archive members and the
driver byte-for-byte, and requires the certified baseline executable to match
its original. It checks that both candidates contain the generic policy wrapper
and that only the diagnostic twin contains policy telemetry, with the expected
Apptainer compile/link commands. Syntax and preparation hashes passed; the
review itself has not executed. Binary provenance does not prove runtime policy
selection, numerical acceptance, or competitive speed.

### Completed application-managed SPICE comparison (2026-10-11)

Controller 203776 was absent and exec session 29644 exited zero after all 462
original observations completed. The frozen independent reviewer also exited
zero, rechecking the complete unique grid, balanced order, hashes, logs, raw
waveforms, upstream verification, qualification, and all twelve aggregates.
All 112 certified-budget KLS runs and 112 application-managed candidate runs
qualified. The 112 invalid CKTSO runs and two SubtreeLU MOS13 eight-thread
timeouts remain included; no CKTSO or complete eight-thread SubtreeLU aggregate
is available.

Candidate/certified-budget geometric-mean linear-solution-time ratios were
0.82014/0.85108 on llc32 at one/eight threads and 0.83124/0.87036 on llc96.
All 28 case/configuration medians improved. The candidate includes both the
qualified lazy-tolerance replay change and generic application-managed policy;
these ratios do not isolate the policy's contribution. Complete one-thread
candidate/SubtreeLU ratios were 0.66921 and 0.69720. Full overall superiority
remains unproven. Review:
`build/application-managed-competitive-spice-review-preparation-20261011/review.json`
(SHA256 `9b02b1e92a9ad034b2acd13660bca81d25ad972f95250340343544b9273519bc`).

The deferred general-policy build subsequently completed in Apptainer with
actual PID 234335 absent and exec session 59643 exiting zero. Its independent
provenance reviewer exited zero: all seven archive members and the driver match,
the certified baseline is byte-identical, both candidates contain the generic
selector, and the plain candidate has no policy telemetry symbols. This build
does not prove general runtime acceptance. Review:
`build/application-managed-general-policy-build-review-preparation-20261011/review.json`
(SHA256 `8a8b3833d9335e19458407bb0c392c910ec5af89806af6a1ec91920084cc2225`).
The decision advances complete general-policy validation while retaining every
failure and the production default. Numerical source promotion and the ultimate
fastest-in-SPICE objective remain unproven. Decision:
`build/application-managed-completed-spice-decision-20261011/decision.json`.

The separate complete general-policy campaign is now running as actual
controller 235660 / exec session 75152. Preparation:
`build/application-managed-general-policy-campaign-preparation-20261011`;
launcher and independent post-terminal reviewer use the corresponding
`application-managed-general-policy-campaign-{launch,review}-preparation-20261011`
directories. The original manifest, configurations, runner, and entire timing,
validation, and aggregation loop were checked byte-identical to the preceding
general campaign. All 110 matrices, two CPU domains, four balanced pairs,
100 changed-value systems, 900-second timeout, original owner/value/ordering
controls, and independent per-refactor residual acceptance at 1e-8 remain.
Both arms use the exact same qualified archive and driver; the candidate adds
only the generic application-managed creation selector. Launch gates checked
complete SPICE review, actual preceding writer absence and exact terminal
receipts, binary provenance, pinned sources, and no overlapping heavy work.
No retries, exclusions, partial speed review, or default change are allowed.
The prepared reviewer requires actual writer absence and terminal exec evidence,
then reparses every log/result, verifies the full unique grid/order/hashes,
checks finite residual acceptance, retains every failure, and independently
recomputes coverage regressions and aggregates. General runtime acceptance
remains unproven until this complete review finishes.

The general-policy campaign completed on 2026-10-11 with actual controller
235660 absent and exec session 75152 exiting zero. The independent reviewer
also exited zero (session 11851), verifying all 1760 observations, original
grid/order, inputs, logs and acceptance criteria. It retained 104 invalid
observations and identified 24 candidate-only coverage regressions: all four
trials in both CPU domains for `mac_econ_fwd500`, `power197k`, and
`tsopf_fs_b9_c1`. Overall general timing aggregates are null because coverage
is incomplete. Review:
`build/application-managed-general-policy-campaign-review-preparation-20261011/review.json`
(SHA256 `43b372c4afbfaf6e2d73b872ab0c0f68507631a4dcfd816d5efe30f648b2726e`).

The completed general decision rejects a global production default-off change
despite the separate successful SPICE qualification. The certified production
default remains unchanged; no matrix or circuit exceptions are introduced.
The next candidate caches exact coefficients derived only from sparse structure
and retains certification and original fallback behavior. Its isolated source,
coefficient/FP-state tests, optional-allocation-failure tests, launch guards and
reviewers are prepared but remain uncompiled and unqualified. Decision:
`build/application-managed-completed-general-policy-decision-20261011/decision.json`.
CKTSO runtime recovery and complete valid competitive evidence remain required
before any fastest-in-SPICE claim.

The structural coefficient-cache candidate subsequently completed isolated
qualification in Apptainer. The cache-present and forced optional-allocation-
failure campaigns each passed five builds (including O0), 146250 range
comparisons and 102400 leak-enabled ASan/UBSan calls per arm. The coefficient
audit passed 4485120 encoding checks and 30720 real-plan construction FP-state
checks. The complete retained ten-build tolerance-boundary sweep passed
3252720 status/bound/FP-state comparisons. Actual writers 412752, 422251 and
431238 are absent; exec sessions 85698, 49031 and 76673 exited zero, and all
three independent reviewers passed. Combined evidence:
`build/componentwise-structural-coefficient-cache-isolated-qualification-review-20261011/review.json`.
This qualifies the isolated numerical change only. Fixed-grid micro
measurement, full native/sanitized runtime and public lifecycle validation,
complete SPICE/general campaigns and valid competitive evidence remain
required. No production source, default or numerical promotion follows.

The structural coefficient cache was then rejected for performance after all
5952 fixed-grid micro observations completed and passed independent review.
Controls restore the qualified lazy-tolerance baseline row/range bodies by
inverse-checked identifier renames; fixtures, trial order and timers match the
preceding grids. All native/generic aggregates regress: structural ranges
1.021185/1.011900, row order 1.035553/1.026140, opaque structural ranges
1.017488/1.011326 and opaque row order 1.032644/1.024102 (candidate/baseline).
Every slower group is retained. The candidate is not promoted or advanced to
expensive full-solver performance testing. Decision:
`build/componentwise-structural-coefficient-cache-completed-micro-decision-20261011/decision.json`.
The next hypothesis simplifies exact coefficient arithmetic without cache
storage, preserving original arithmetic outside its proven exact domain.

The separate per-row guarded exact-coefficient arithmetic candidate passed
isolated qualification in ten compiler modes: 19169280 coefficient/full-FP
comparisons, 146250 range comparisons, 102400 leak-enabled sanitizer calls
per arm and the retained 3252720 tolerance-boundary comparisons. Reversing
the helper/one-expression edit recovers the qualified source exactly; plan,
allocation, original scalar and range-replay source remain unchanged.
Evidence:
`build/componentwise-exact-coefficient-arithmetic-isolated-qualification-review-20261011/review.json`.

All 5952 subsequent fixed-grid micro observations completed and passed
independent review. Results are mixed: structural ranges 1.012507/0.989065,
row order 0.988878/0.995857, opaque structural ranges 1.009498/0.977931 and
opaque row order 0.980904/0.994984 (native/generic candidate/baseline). Native
structural regressions and every slower group are retained; no promotion or
full-solver speed claim follows. Decision:
`build/componentwise-exact-coefficient-arithmetic-completed-micro-decision-20261011/decision.json`.
The next hypothesis moves the exact-domain proof to a range-level validated-
NNZ guard while retaining the qualified algorithm outside that domain.

The range-level coefficient-domain candidate completed isolated qualification
in ten compiler modes, including 14315520 coefficient/full-FP comparisons,
49920 actual routing diagnostics, 146250 range comparisons, 102400 sanitizer
calls per arm and 3252720 tolerance-boundary comparisons. A pre-compilation
boundary harness filename collision is retained in its original failed
package; a separately prepared corrected package passed. Evidence:
`build/componentwise-range-coefficient-domain-isolated-qualification-review-20261011/review.json`.

All 5952 fixed-grid micro observations passed independent review. Structural
ranges regress at 1.069438/1.052324 and opaque structural ranges at
1.044438/1.048537; row-order grids improve at 0.934196/0.919596 and opaque
row order at 0.935076/0.920469 (native/generic candidate/baseline). All slow
groups remain included. This candidate is not promoted or advanced to
full-solver performance testing. The order execution session expired across
context compaction; its successful terminal controller receipt and actual PID
absence were recovered, with that limitation recorded explicitly rather than
claiming a recovered PTY exit status. Decision:
`build/componentwise-range-coefficient-domain-completed-micro-decision-20261012/decision.json`.
The production certification default remains unchanged. Complete valid
competitive coverage is still required for the active SPICE performance goal.

An untimed Apptainer public-API stage probe uses the original CKTSO shared
library and the documented library directory. `CKTSO_CreateSolver` returns
-99 before analysis; stderr is empty. The completed competitive executable's
dynamic section directly requires the original `/work/cktso/.../libcktso.so`,
so a copied-library-path explanation is unsupported. This diagnostic does not
establish license validity, expiry or exact cause, and does not retroactively
identify the Xyce failure stage. No license-file contents or metadata were
inspected. No competitive retry or timing claim follows. Evidence:
`build/cktso-public-api-stage-probe-preparation-20261012/review.json`.

A separate candidate records the exact coefficient-degree domain once in the
existing private plan eligibility field, after full CSC validation and stable
row scatter. Plan layout and allocation sequence are unchanged; public
reuse retains dimension/NNZ reset checks, and the original valid eligibility
state takes the unchanged qualified fallback. There are no compiler, matrix,
circuit or history selectors. Source edits are inverse-checked against the
previous range-domain candidate.

The candidate passed isolated qualification in ten compiler modes: 14315520
coefficient/full-FP checks, 49920 actual routing diagnostics, 146250 range
comparisons, 102400 sanitizer calls per arm and 3252720 retained raw tolerance
boundary comparisons. Private large-NNZ controls test metadata dispatch;
they are not enormous valid public CSC fixtures. Evidence:
`build/componentwise-plan-coefficient-domain-isolated-qualification-review-20261012/review.json`.

All 5952 fixed-grid micro observations then passed independent review.
Native/generic candidate/baseline geometric means are structural ranges
0.940695/0.915627, row order 0.923865/0.911163, opaque structural ranges
0.936580/0.934867 and opaque row order 0.891233/0.918600. All slower groups
remain included. Plan construction is excluded from these isolated timers;
full solver setup and reuse still require measurement. Decision:
`build/componentwise-plan-coefficient-domain-completed-micro-decision-20261012/decision.json`.
This advances the candidate to runtime and public lifecycle qualification,
without production/default changes or a SPICE performance claim.

Runtime preparation's initial assertion-count mistake (five expected private
assertions instead of the actual four) stopped before compilation. The
separate corrected package changes only those four expected private positive
state values from 1 to 2, inverse-checked; all behavioral checks remain.
Its copied host then failed a package-name preflight because a separator was
missing, also before compilation. Both failures are retained. A separately
prepared corrected launcher uses the unchanged pinned source/plan/controller,
and a separate reviewer checks the actual successful controller's receipts.
Python bytecode generation is disabled for the new runtime tests.

The plan-construction degree-domain marker then completed all 88 runtime
checks: 44 native and 44 optimized ASan/UBSan tests, with leak detection and
halt-on-error enabled. The independent review checks all unique CTest gates,
prepared source hashes and archive identity; only `kls_componentwise.c.o`
differs from the pinned runtime baseline. Actual controller 477256 is absent
and exec session 27663 exited zero. Review:
`build/componentwise-plan-coefficient-domain-runtime-corrected-preparation-20261012/completed-review.json`
(SHA256 `202ebe936d1b5cc4729d1fb5ca271d299f121bb40c0a3f99308f7ff1e6d0698b`).
The successful control and reviewer are the separate runtime-launch-corrected
packages; the earlier failed prepared host remains preserved. The original
full28 public API campaign has been launched with a gated controller and an
independent reviewer prepared to rerun every waveform comparison. Public API,
full solver performance/general and valid competitive evidence remain pending;
production/default are unchanged and the SPICE goal remains active.

The plan-construction marker subsequently passed all original 28 public API
configurations: seven full circuits, two CPU domains and one/eight threads.
Every configuration retains positive certificate-call counts, zero certificate
and stats errors, matching steps/end/finite-output gates and unchanged
waveform tolerances. The independent reviewer reran all 28 upstream waveform
comparisons in separate folders. Actual controller 480119 is absent and exec
session 17638 exited zero. Evidence:
`build/componentwise-plan-coefficient-domain-api-qualification-preparation-20261012/qualification-review.json`
(SHA256 `676d5e9a590171e30e9988c307f8e678077ae42c5b67617e0127a42c7ca5f3f1`).

The first full-SPICE timing package stopped at executable-identity preflight:
a replacement accidentally made its BASELINE variable refer to the candidate.
Actual controller 482324/session 75515 exited one before campaign-output
creation or timing; that package is retained. A separately prepared corrected
package verifies actual parsed Python executable assignments against distinct
qualified API metadata, checks separate hashes, and restores the qualified
lazy-tolerance baseline. Both archives retain the same six other members and
both executables are uninstrumented. The original fixed fixture/order/timer
body is byte-identical. The corrected full238 campaign is now running;
full238 independent timer/log/ratio and waveform review is prepared and will
run only after actual controller absence and terminal execution. Canonical
controller package:
`build/componentwise-plan-coefficient-domain-spice-pair-corrected-preparation-20261012`.
There is no partial timing claim, production/default change or promotion.
Complete SPICE/general and valid competitive evidence remain required.

The corrected plan-construction marker full-SPICE campaign completed all 238
observations and passed independent review, including 238 waveform rechecks
against fresh KLU references and archived valid KLU outputs. Actual controller
482385/session 80608 exited zero and is absent. All observations qualified;
the full unique grid, balanced order, timer logs and input hashes were checked.
Candidate/qualified-lazy-baseline geometric means are 0.985692 (llc32,t1),
0.980100 (llc32,t8), 0.983335 (llc96,t1) and 0.997791 (llc96,t8). Two slower
medians in llc96,t8 remain included: BJT 1.000357 and neuron 1.060826. No
case/domain/compiler exceptions or qualified-subset claim are introduced.
Review:
`build/componentwise-plan-coefficient-domain-spice-pair-corrected-review-preparation-20261012/review.json`
(SHA256 `aca391fb989e36ac7de97a5789567c47c7bb5e7c2a60aed8a861745d3037ea48`).
Decision:
`build/componentwise-plan-coefficient-domain-completed-spice-decision-20261012/decision.json`.

The candidate advances to the unchanged complete general-matrix campaign.
Fresh isolated links use the independently qualified lazy-tolerance archive
as baseline, a common rebuilt driver, pinned dependencies and the full88
qualified native candidate component. Independent build review confirms only
`kls_componentwise.c.o` differs and six other archive members remain identical.
Preflight checks the parsed executable assignments, and the whole fixture,
qualification, trial-order and timer body is inverse-byte-identical to the
previous campaign. All original 110 matrices, two eight-thread CPU domains,
four balanced pairs and 100 changed-value systems per observation are kept:
1760 observations, per-system/max/final residual gates of 1e-8, and the original
900-second timeout. The canonical live campaign is:
`build/componentwise-plan-coefficient-domain-broad-campaign-preparation-20261012`.
The complete independent reviewer is prepared; partial aggregate claims are
not made and every failure/timeout will be retained. Production and the
certification default remain unchanged. Complete valid competitive coverage
is still required before claiming the active fastest-in-SPICE goal achieved.

The plan-construction marker general campaign completed all 1760 observations.
Actual controller 498929/session 45466 exited zero and is absent; independent
review session 5675 also exited zero. The reviewer checked the exact unique
grid and balanced order, source/executable hashes, retained stdout/stderr, and
independently reparsed every numerical result. There are zero candidate-only
qualification regressions. Eighty invalid observations are retained, sixteen
per matrix across both arms and domains: hamrle3 (timeouts), fullchip and
freescale2 (process failures), and bips07_1998 and fpga_dcop_01 (residual gates).
Both complete-grid timing aggregates remain null. The 37 slower qualified
matrix medians in llc32 and 42 in llc96 are retained without exceptions or a
qualified-subset overall claim. This does not establish successful solves on
the five common-invalid matrices.
Review: `build/componentwise-plan-coefficient-domain-broad-review-preparation-20261012/review.json`.
Decision: `build/componentwise-plan-coefficient-domain-completed-general-decision-20261012/decision.json`.
The candidate remains private and unpromoted; production certification remains
enabled. General value-buffer reuse can now proceed to isolated implementation
and qualification. Its static rejection-order requirements are recorded in
`build/original-operator-reuse-rejection-order-investigation-20261012/investigation.json`;
alias handling, rejection atomicity and full FP behavior still require runtime
proof. Complete valid competitive SPICE coverage remains unresolved.

The isolated original-operator reuse candidate now passes eight targeted cache
runs (native, generic, unoptimized generic, and ASAN/UBSAN, both arms) and
301,989,888 disjoint-update differential checks across three builds. The FP
checks compare status, full MXCSR and x87 state, published values and structure,
including nonfinite rejection and all tested rounding/precision/flush/flag
combinations. Exact layout checks remain mandatory. Disjoint inputs reuse the
existing values buffer; overlapping or unprovable spans retain allocation.
Application-managed updates validate every gathered value before publication.
These are isolated correctness results, not full solver or performance proof.
Reviews: `build/original-operator-disjoint-value-reuse-cache-review-20261012/review.json`
and `build/original-operator-disjoint-value-reuse-fp-review-20261012/review.json`.

The first full runtime run is retained as a failure: 43 native tests passed,
and `kls_accuracy_owner` required the removed allocation on a disjoint cache
hit; sanitizer testing had not started. The separate corrected package adds
an explicit allocation-free reuse assertion and exercises the existing
allocation-failure/epoch-invalidation checks with overlapping values. All six
mandatory structural allocation failures remain covered. The solver source
is unchanged between these runs. The corrected full native/sanitizer run is
in progress; public API qualification is prepared and gated on its independent
completed review. Production certification remains enabled and this candidate
remains unpromoted.
Failure review: `build/original-operator-disjoint-value-reuse-runtime-failure-review-20261012/review.json`.
Corrected preparation: `build/original-operator-disjoint-value-reuse-runtime-corrected-preparation-20261012/preparation-receipt.json`.

The corrected reuse runtime run completed with all 44 native and 44 optimized
ASAN/UBSAN tests passing. Actual controller 677938/session 63031 exited zero
and is absent; independent reviewer session 24347 exited zero. The review
checks source and executable hashes, every test log, and all archive members;
only `kls.c.o` differs from the qualified marker archive. The ownership test
retains overlapping-update OOM and all six structural OOM cases alongside
the new allocation-free reuse assertion. The failed first package remains
unchanged. Review: `build/original-operator-disjoint-value-reuse-runtime-corrected-preparation-20261012/completed-review.json`.
The prepared original 28 public API configurations have now been launched
with certification enabled. Their result, uninstrumented timing, general
qualification and complete peer coverage remain pending. No promotion or
production/default change follows from these correctness-only checks.

The first API preparation failed during compilation, before any circuit ran:
the inherited component-only include paths did not provide `scotch.h` for
`kls.c`. It is preserved with its logs and terminal status. A separate
corrected package uses the passing full-runtime dependency include paths and
pins the custom headers from that build's dependency file. Compiler defines
and options, solver sources, public wrappers, and the full28 controller are
retained. The generated Scotch header uses 64-bit numbers, matching the
original Xyce CMake cache. Actual corrected API controller 681115/session
75038 is live; its completed result remains pending. This is a preparation
correction, with no circuit or timing observations retried.
Failure: `build/original-operator-disjoint-value-reuse-api-build-failure-review-20261012/review.json`.
Corrected preflight: `build/original-operator-disjoint-value-reuse-api-corrected-preflight-review-20261012/review.json`.

The corrected full28 API run completed: controller 681115/session 75038
exited zero and is absent. Independent reviewer session 82461 exited zero,
reparsed every public profile/error/step log, verified all source and binary
hashes and six retained archive members, and reran all 28 waveform checks
in independent folders. Every configuration is qualified and certificate/API
error counts are zero. The original failed build remains unchanged.
Review: `build/original-operator-disjoint-value-reuse-api-qualification-corrected-preparation-20261012/qualification-review.json`.
Decision: `build/original-operator-disjoint-value-reuse-completed-correctness-decision-20261012/decision.json`.
This advances the general reuse candidate to uninstrumented fixed full-SPICE
and general matrix performance measurement, with failures and slower cases
retained. It establishes no speed advantage, makes no benchmark-specific
exception, and does not change production or the certification default.
The fastest-in-SPICE goal still requires complete valid competitive coverage.

The reuse candidate's uninstrumented fixed full238 SPICE comparison is now
running under actual controller 683356/session 2511. The baseline is the
qualified certified plan-marker solver; the candidate adds only value-buffer
reuse. Independent preflight verifies distinct executable assignments, exact
original fixture/trial-order/timer body, absence of profiling wrappers, six
unchanged archive members, and no overlapping heavy jobs. All seven original
circuits, two domains, t1/t8, four balanced pairs and fourteen fresh KLU
references are retained, with unchanged waveform and solver gates.
No partial ratios are reviewed. Actual writer absence and terminal execution
are required before independent full-grid log/ratio/waveform review. General
matrix qualification and complete valid peer coverage remain pending;
production/default/promotion remain unchanged.
Preparation: `build/original-operator-disjoint-value-reuse-spice-pair-preparation-20261012/preparation-receipt.json`.
Preflight: `build/original-operator-disjoint-value-reuse-spice-pair-preflight-review-20261012/review.json`.
Reviewer: `build/original-operator-disjoint-value-reuse-spice-pair-review-preparation-20261012/review.py`.

The reuse full238 campaign completed: actual controller 683356/session 2511
exited zero and is absent. Independent reviewer session 27116 exited zero,
verified the complete unique grid and exact balanced order, reparsed every
timer/error/step log, checked pinned inputs and reran all 238 waveform gates.
All observations qualify. Candidate/control full-case geometric means are
1.1360889435 (llc32 t1), 1.1495505077 (llc32 t8), 1.1586208133 (llc96 t1),
and 1.1580616282 (llc96 t8). The measured implementation is slower overall
in every configuration. All ratios and slower cases remain recorded,
including the small faster BJT median in llc32 t8. It is not promoted and
does not advance to general timing. Before attributing the cause or trying
another general implementation, a separate preparation recompiles the frozen
allocating baseline with the exact candidate compiler/options/dependency
settings for object and compiler-loop inspection. This is not a timing retry.
Review: `build/original-operator-disjoint-value-reuse-spice-pair-review-preparation-20261012/review.json`.
Hold decision: `build/original-operator-disjoint-value-reuse-completed-spice-hold-decision-20261012/decision.json`.
The cached gather's bounded input-span invariant is independently recorded
in `build/original-operator-disjoint-reuse-gather-span-review-20261012/review.json`;
that static proof changes no source and establishes no performance advantage.

After timing completed, the untimed original-library-directory CKTSO probe
also returned `create -99`; actual controller 698417 exited one and is absent.
Its working-directory change did not restore creation. The original binary,
library, CPU, container and 1x1 problem were retained. Neither diagnostic
establishes license validity, expiry, colocation or exact cause. License-file
contents and metadata remain uninspected; no version substitution or time
manipulation occurred. Competitive CKTSO coverage remains unavailable.
Review: `build/cktso-original-library-cwd-stage-probe-preparation-20261012/review.json`.
Production and certification default remain unchanged; the fastest-in-SPICE
goal remains unachieved.

The fresh allocating baseline object is byte-for-byte identical to the
qualified control's `kls.c.o`, including every code and relocation section.
The compiler/dependency refresh is therefore not a baseline mismatch in this
comparison. Native disassembly of the rejected reuse branch shows eight
separate byte loads and stores per copied double and snapshot/map pointer
reloads per entry, whereas the allocating baseline emits a word-sized copy.
This is a plausible contributor, not a measured attribution of the entire
regression. A further general candidate may hoist these pointers and express
the proven disjoint read/write regions with local `restrict` declarations;
owned map/header allocations are distinct from the destination. Exact
structure checks, finite-rejection order, object representations and alias
fallback must remain intact. No revised source, runtime or speed claim follows
from this static inspection.
Object review: `build/original-operator-reuse-fresh-baseline-object-review-20261012/review.json`.
Loop investigation: `build/original-operator-reuse-compiler-copy-loop-investigation-20261012/investigation.json`.

The first local-pointer preparation stopped at its source-replacement guard
before any edit, compilation or test: the prepared block's indentation did
not match. It is retained unchanged; a separate corrected source package
completed eight cache audits and 301,989,888 full-FP/published-byte checks,
all passing. Its production disassembly hoists snapshot pointers but still
emits eight byte loads/stores per value. It does not advance to full runtime
or timing and has no speed claim.
Review: `build/original-operator-disjoint-local-pointer-reuse-corrected-combined-qualification-review-20261012/review.json`.

A separate copy-contract candidate retains the validated local pointers and
adds `restrict` to the private byte-copy helper's parameters. Every helper
call copies to a fresh allocation or between already-validated disjoint
spans; the bounded gather and separately owned map/header establish that
contract without a stronger caller requirement. The byte-copy body, exact
structure keys, finite-rejection order, allocating alias fallback, layouts,
policies and public ABI remain unchanged. It adds no CPU/compiler/shape or
benchmark selector. This candidate passes eight targeted cache audits and
301,989,888 full-FP/published-byte checks across the same compiler grid.
Production native disassembly now uses word-sized copies without per-value
byte loads/stores. The independent review verifies terminal execution,
actual controller absence, binaries, logs and inverse-exact parameter delta.
Review: `build/original-operator-disjoint-copy-contract-reuse-combined-qualification-review-20261012/review.json`.
Its fresh full native/optimized-sanitizer runtime suite is now running. These
isolated results establish no full-solver speed advantage; public API, fixed
full-SPICE and general matrix qualification remain required. Production,
certification default and promotion remain unchanged.

The copy-contract full runtime run completed with all 44 native and 44
optimized ASAN/UBSAN tests passing. Actual controller 699818/session 16618
exited zero and is absent; independent reviewer session 27228 exited zero.
The reviewer verifies all tests and artifacts against preparation hashes,
checks exact candidate sources, and confirms only `kls.c.o` differs from the
qualified allocating-marker archive. Mandatory alias/structural OOM and
epoch checks remain covered alongside allocation-free disjoint reuse.
Review: `build/original-operator-disjoint-copy-contract-reuse-runtime-preparation-20261012/completed-review.json`.
The prepared full28 public API qualification has been launched with the
original circuits, domains, thread counts, public wrappers and waveform
gates, certification enabled, and known-good pinned dependency headers.
Completed API review and uninstrumented timing are pending. This candidate
remains private; production and certification default remain unchanged.

The copy-contract full28 API qualification completed: actual controller
702551/session 70554 exited zero and is absent. Independent reviewer session
19610 exited zero, checked exact candidate sources and six unchanged archive
members, reparsed every public profile/error/step log and reran all 28
waveform gates independently. All configurations qualify with certificate
and API error counts zero. Review:
`build/original-operator-disjoint-copy-contract-reuse-api-qualification-preparation-20261012/qualification-review.json`.
The completed correctness decision advances this general implementation to
fixed full-SPICE and general measurement, without a speed or promotion claim:
`build/original-operator-disjoint-copy-contract-reuse-completed-correctness-decision-20261012/decision.json`.

The copy-contract original full238 uninstrumented comparison is now running
under actual controller 704732/session 32941. Independent preflight checks
distinct executable identities, both unwrapped link commands, exact prior
fixture/trial-order/timer body, current isolated/full88/full28 correctness
gates, six unchanged archive members and absence of overlapping heavy jobs.
The qualified allocating marker remains the control. All seven original
circuits, both domains, t1/t8, four balanced pairs and fourteen fresh KLU
references retain unchanged waveform and solver gates. No partial ratios
are reviewed. Independent full-grid log/ratio/waveform review requires actual
writer absence and terminal execution. The earlier slower candidate and every
failed preparation remain unchanged. Production/default/promotion and the
fastest-in-SPICE goal remain unchanged.
Preparation: `build/original-operator-disjoint-copy-contract-reuse-spice-pair-preparation-20261012/preparation-receipt.json`.
Preflight: `build/original-operator-disjoint-copy-contract-reuse-spice-pair-preflight-review-20261012/review.json`.

After environment resume, the copy-contract full238 campaign is independently
reviewed complete: controller 704732/session 32941 exited zero and is absent;
reviewer session 86021 exited zero. Every observation and all 238 independent
waveform reruns qualify. Candidate/control full-case geometric means are
1.0021932861 (llc32 t1), 1.0063222241 (llc32 t8), 1.0113941058 (llc96 t1),
and 1.0088802017 (llc96 t8). The earlier large regression is removed, but
this implementation still has no measured aggregate advantage in any
configuration. All faster/slower cases and every pair remain retained; no
circuit/CPU/compiler exception or timing retry is introduced. It is held
unpromoted and does not advance to general timing.
Review: `build/original-operator-disjoint-copy-contract-reuse-spice-pair-review-preparation-20261012/review.json`.
Decision: `build/original-operator-disjoint-copy-contract-reuse-completed-spice-hold-decision-20261013/decision.json`.

A fresh bounded public CKTSO diagnostic after resume still returns
`create -99`. Controller 744051 exited one and is absent. Original binary,
library, working directory, image and CPU setup remain unchanged; no license
file contents/metadata or validity/expiry/cause inference is involved.
Review: `build/cktso-original-library-cwd-resume-stage-probe-preparation-20261013/review.json`.

Static investigation now targets the numeric-to-certificate dataflow:
an owned caller-order bulk copy could avoid gathering all values at every
numeric update, while composing cached row positions once. This is not yet
an implementation. Planned transpose, scalar/interval fallback, correction,
strict/cold recovery readers and mapping/format reset rules must be audited;
every logical read and arithmetic encounter order must remain equivalent.
Optional metadata failure must retain the complete canonical fallback without
an extra required NNZ allocation or caller-lifetime assumption.
Investigation: `build/original-operator-caller-order-certificate-dataflow-investigation-20261013/investigation.json`.
Production, default certification and promotion remain unchanged; valid peer
coverage and the fastest-in-SPICE goal remain unproven.

An isolated mapped-certificate prototype is now implemented privately;
snapshot storage and strict/cold recovery remain unchanged. Normal row plans
compose their value positions during stable scatter, and transpose plus
stateless/interval readers translate canonical indices. Existing canonical
entry-point signatures remain intact; new mapped entry points are private.
The first native probe failed on a reverse-mapped transpose planned check:
four tolerance-cache reads had not yet been translated. That source and
binary are retained. A separate corrected package fixes those four reads.
Failure: `build/componentwise-mapped-values-prototype-first-probe-failure-20261013/review.json`.
Corrected source: `build/componentwise-mapped-values-prototype-corrected-source-preparation-20261013/preparation-receipt.json`.

The corrected isolated audit completes 48,384 candidate canonical-versus-
mapped checks in each of O3 native, O3 generic, O0 generic and ASAN/UBSAN:
193,536 checks total. Actual audit controller exited zero and is absent;
exec session 57269 is terminal. Status, full SSE/x87 state, upper-bound
significand/exponent bytes and correction bytes match on the recorded grid
(three value permutations, seven finite/nonfinite fixtures, both directions,
scalar/interval/planned paths, four rounds, three x87 precisions, four
FTZ/DAZ settings and eight paired exception-flag seeds). This does not yet
prove equivalence to the original frozen certificate implementation or cover
all shapes, allocation faults, callbacks and mapping/format reset transitions.
Those checks are required before snapshot integration. No full solver,
performance, production/default or promotion claim follows.
Review: `build/componentwise-mapped-values-prototype-corrected-audit-review-20261013/review.json`.

The mapped prototype's canonical behavior now matches the original frozen
certificate implementation on the recorded fixture/FENV grid. Baseline and
candidate complete 70-byte status/FP/bound/correction records are byte-for-byte
equal in native, generic, unoptimized generic and sanitized builds; each trace
contains 48,384 records. Internal mapped-versus-canonical equality remains
checked. Exec session 14073 exited zero and its controller is absent.
Review: `build/componentwise-mapped-values-original-baseline-trace-review-20261013/review.json`.

Expanded shape/fault tests complete 2,184 checks per compiler configuration
(8,736 total), including duplicates, empty/zero-dimension and malformed
patterns, both directions, permutations, required/optional allocator failures,
successful/declined/after-work-failed executors, and mapped-to-canonical-to-
mapped plan transitions. Each run records 570 injected faults and 108 executor
dispatches. Two failed probe versions are retained: the first compared a
retained owner's allocation count against a fresh owner; the second read
unpromised correction scratch after OOM. Separate corrected probes compare
promised semantics for resets and only the first-n correction on uncertified
finite-bound returns, while paired fault allocation counts remain exact.
Certificate source is unchanged for these probe corrections. All four builds
pass; exec session 2116 exited zero and its controller is absent.
Review: `build/componentwise-mapped-shape-fault-audit-review-20261013/review.json`.
Decision: `build/componentwise-mapped-values-prototype-expanded-qualification-decision-20261013/decision.json`.
This advances to isolated owned snapshot integration, retaining canonical
metadata-failure and strict/cold semantics, mapping/format reset requirements,
and every publication/FENV/correction/OOM invariant. It establishes no full
solver or performance advantage; production/default/promotion remain unchanged.

Owned caller-order snapshot integration is now isolated privately. Successful
optional layout/map construction allows componentwise mode to own a fresh
bulk copy, with the validated source map describing canonical positions.
Strict-L2, application-managed and metadata-failure paths retain canonical
storage. Exact layout checks, fresh required allocation/OOM behavior and
atomic publication remain. Mapping/format changes reset row plans; an
identical map allocation replacement also resets before freeing its borrowed
pointer. That lifetime condition was caught in an unexecuted first source
preparation and corrected separately, without treating pointer identity as a
content proof. The corrected snapshot probe passes 1,096 logical matrix,
ownership, update, orientation/permutation, self-alias, OOM, metadata-fallback
and plan-lifetime checks in each of four compiler/sanitizer modes (4,384 total).
The integrated main solver also compiles with the pinned production options.
Isolated review: `build/caller-order-snapshot-integration-isolated-review-20261013/review.json`.
Lifetime review: `build/caller-order-snapshot-integration-map-lifetime-preflight-20261013/review.json`.

The first full runtime preparation rebuilt both `kls.c.o` and
`kls_componentwise.c.o` and retained the other five members. Forty-two native
tests passed; two test adapters failed: direct private differential calls
needed the nullable mapping argument, and the caller fault adapter intercepted
only the canonical certificate boundary. Its complete failure and terminal
status are preserved; sanitizer testing had not started. A separate corrected
package changes only those two adapters, retaining every caller assertion,
OOM/rejection counter and publication gate. Solver source and runtime
controller are byte-exact to the failed package. It is now running the fresh
full native/optimized-sanitizer suite; its result, full API/waveform
qualification and performance are pending. No production, default or
promotion change follows.
Failure: `build/caller-order-snapshot-integration-runtime-first-failure-review-20261013/review.json`.
Corrected preflight: `build/caller-order-snapshot-integration-runtime-adapter-preflight-review-20261013/review.json`.

## Qualification status, 2026-10-17

Main numerical source remains at `2754413`, including the EGraph startup and
atomic idle-budget fixes. Xyce's qualified cacheline import is committed as
`a15e69e`. Production componentwise verification remains enabled. The complete
application-managed comparison added 24 qualification failures, so its
default-off proposal was rejected.

The current-root certificate integration is committed separately on
`work/certificate-verdict-root-20261017` (`e3cc481`, `24dc56f`). Independent
review passes all 43 native and 43 optimized ASan/UBSan tests and all 28
circuit/API settings, covering 288,650 certificate calls. Its native
componentwise object is byte-identical to the earlier qualified verdict
object. A fresh complete 238-observation paired cost campaign is running;
current-root performance, broader-corpus and peer qualification remain pending.
The separate unused-L2-copy guards are committed on
`work/componentwise-l2-copy-guards-20261017` (`2f051ec`). All runtime/API checks
and both complete 238-observation campaigns qualify, with small per-case cost
regressions retained separately; the full 110-matrix stage remains pending.
Neither branch is promoted to production. Mapped-worker, caller-order snapshot
and packed-vector experiments remain separate; their gains are not combined
into a claimed speedup. The fastest-in-SPICE objective remains unmet.

The caller-order snapshot integration remains held after two added Hamrle3
timeouts in its complete generality campaign. The isolated verdict-witness
candidate completed all 1,760 observations and an independent raw-log review
with zero added failures. All 80 shared invalid observations remain retained;
its overall generality performance aggregate is null. Evidence is in
`build/componentwise-verdict-witness-completed-general-decision-20261016/decision.json`.

The fresh verdict-witness peer campaign completed all 462 observations.
Session 67958 exited zero, its controller and numerical children were absent,
and an independent reviewer rechecked every raw log, waveform, qualification,
comparison order and ratio. Every KLS run qualified. Two SubtreeLU MOS runs in
llc96 at eight threads timed out; those observations remain invalid and the
corresponding aggregate remains null. Of the 56 per-circuit/domain/thread/peer
comparisons, 23 favor KLS, 32 favor the peer, and one lacks qualified peer
coverage. The four suite ratios against the marker parent are 0.985845,
0.987395, 0.984019 and 0.978584 (llc32 t1/t8, then llc96 t1/t8).
These isolated gains do not establish the fastest-in-SPICE goal and must not
be multiplied by the separate mapped-worker gains. No incomplete or favorable
subset supplies an overall competitor comparison. The completed decision is
`build/componentwise-verdict-witness-completed-peer-decision-20261016/decision.json`.

CKTSO's installed license file was empty. Installing the packaged 719-byte
license restored successful API calls and full-circuit qualification without
changing the solver library, key contents, or clock. This was not evidence of
license expiration. The installation review is
`build/cktso-packaged-license-installation-20261015/review.json`.

### Finer refactor diagnostics and narrow-edge prototype

A diagnostic rebuilt only the main archive member and completed the full
28-circuit/domain/thread grid. After observed terminal exit and process absence,
an independent reviewer rechecked all counters, source/archive pins and fresh
waveforms. It separates mapped/lean stages and ordinary supernode scatter,
dependency waits, edge updates, division, diagonal updates and writeback.
BJT's serial supernode edge updates remain the largest measured stage.
Per-column clocks dominate several tiny mapped stages and alter low-work
engine selection, so those totals cannot establish removable production cost.
Worker elapsed sums also cannot be treated as wall time. This rejects a
reciprocal-normalization change as the first optimization target.
Review: `build/refactor-exclusive-stage-review-preparation-20261016/review.json`.

The next isolated prototype shares producer loads and destination decoding
across four independent target columns in narrow supernode edges. Each column
retains producer/destination update order; zero-upper handling, wide arithmetic,
scalar tails, allocations, barriers and publication remain unchanged. It adds
no matrix-identity or shape-policy selector. The first byte-exact proof found
a native-build NaN payload difference and is retained. The unchanged kernel
then passed the complete 2,724,736-fixture grid in native, generic and optimized
ASan/UBSan builds: finite results, infinities, signed zeros and exception flags
match exactly; same-class NaN arithmetic payload/sign differences are counted.
This permission does not apply to coefficient ownership or representation copies.
The independent review covers 8,174,208 fixtures, including duplicate destinations,
zeros, nonfinite/extreme values, rounding modes and guard pages. It supplies no
public API, production first-fault ordering, or performance guarantee.
Review: `build/snb-narrow-four-column-semantic-proof-corrected-review-preparation-20261016/review.json`.
The full original 201,600-observation isolated cost grid completed. Both arms
use the same 112-byte block layout. Candidate/baseline geometric mean cost is
1.098424 in the native build and 1.005495 in the generic build; native narrow
edges regress to 1.1686. The general prototype is rejected, without selecting
favorable shapes or adding a benchmark-specific admission rule. Decision:
`build/snb-narrow-four-column-completed-cost-rejection-20261016/decision.json`.

### Packed certificate result

A separate packed certificate prototype groups eight equal-degree rows while
retaining the current operator, RHS and solution and unchanged scalar fallback.
The final runtime-operand variant passed all seven arithmetic profiles with
3,096,576 comparisons and 6,646,304 exact-reference rows. Six native/generic/
sanitizer and literal/runtime-tolerance memory contexts also passed, including
allocation faults, fresh inputs, guards and concurrent immutable readers.
A compiler operand barrier was necessary to preserve scalar x87 exception flags
when a caller supplies a constant tolerance. Earlier failures remain retained.

The complete original 2,304-observation cost grid nevertheless rejects this
prototype. Native kernel cost is 1.040869 of baseline, with 68 of 144 groups
slower; generic kernel cost is 1.000529, with 63 slower. Including separately
measured construction in a modeled batch gives 1.128697 and 1.089416. These
modeled batches are not application timings. No favorable-subset selector or
production integration follows from the correctness proofs. Decision:
`build/componentwise-packed-verdict-runtime-operand-completed-cost-rejection-20261016/decision.json`.

### Adapter recovery and contiguous coefficient import

Xyce commit `7d578c3` checks the complete CSR graph before retaining analysis,
preserves a reanalysis request after a changed graph, and applies column-map
normalization in serial builds as well as MPI builds. Its permanent regression
now changes columns and storage at the same matrix address while retaining
row lengths. The original source fails the controls; corrected complete
adapters pass 2,304 public recovery solves in native and ASan/UBSan contexts.
A further permanent-test extension covers both optimized and unoptimized
Epetra storage: all eight expected contexts pass, including 2,688 corrected
solves and both retained original-source failure controls. Support/Epetra
archives remain release builds; MPI multirank is unqualified.

The separate coefficient-import candidates query fresh contiguous CSR pointers,
check the complete graph, and retain an owned byte-exact coefficient snapshot.
One retains per-row copying; the other compares and copies changed 64-byte
chunks aligned to the owned destination. Both retain the checked row-view
fallback. Neither selects circuits or caches coefficient pointers or answers.

All four frontend arms completed the full 112-run instrumented circuit grid,
with an independent waveform, step-count and raw-log review. Each arm makes
72,142 solver calls per domain/thread setting, all through the direct local
map path. Across the four settings, the cache-line candidate's value-update
phase takes 0.325–0.501 seconds versus 1.147–1.190 for the original and
2.003–2.434 for the graph-checked row-view baseline. These are instrumented
phase totals, not competitive application timings. All arms share the same
isolated qualified verdict-witness archive; the mapped-worker and ordinary-pool
prototypes are absent. Plain binaries have no profiling wrappers.
The balanced 462-observation uninstrumented frontend campaign completed.
Session 58968 exited zero; its controller and all 462 numerical children were
absent before independent review. All observations qualified. Cache-line /
original suite geometric-mean time ratios are 0.915480 and 0.912251 on llc32,
then 0.905409 and 0.896196 on llc96, at one and eight threads. Against the
correct graph-checked row-view baseline the ratios are 0.834658, 0.818478,
0.827663 and 0.818528. These aggregate improvements cover every original
circuit and four balanced passes. BJT's llc32 one-thread difference from the
original is only 0.02%, so it should be treated as effectively flat.

Xyce commit `a15e69e` integrates the exact qualified cache-line implementation
and the optimized/unoptimized storage regression. The decision is
`build/xyce-contiguous-cacheline-import-completed-decision-20261016/decision.json`.
The full timing review is
`build/xyce-normalized-csr-import-four-arm-spice-review-preparation-20261016/review.json`.
No circuit/size/degree/domain selector was introduced. The measured shared core
is still the isolated qualified verdict parent, so these results neither
integrate the separate core experiments nor establish superiority over CKTSO
and SubtreeLU. Fresh combined peer measurements remain required.
Reviews:
`build/xyce-normalized-csr-import-public-recovery-review-preparation-20261016/review.json`,
`build/xyce-csr-import-storage-regression-review-preparation-20261016/review.json`,
`build/xyce-normalized-csr-import-four-arm-api-review-preparation-20261016/review.json`.

### Ordinary-pool caller participation

An isolated core prototype gives ordinary BTF refactor worker zero to the
caller and creates N−1 background threads. It shares that dispatch protocol
with certificate callbacks, preserves block kernels and numerical dependency
barriers, and restores caller affinity. Engine admission and first-factor
scheduling remain unchanged. All 88 native/sanitizer runtime tests pass; an
independent review confirms that only the main archive member changes.
A separate lifecycle proof passes both builds: 48 total allocation failures,
24 partial-startup failures and 600 repeated dispatches, with caller identity,
background joins and affinity restoration checked. This is runtime evidence;
full circuits, complete lifecycle cost, generality and peer timing remain
required before integration.
Reviews:
`build/cooperative-ordinary-pool-runtime-review-preparation-20261016/review.json`,
`build/cooperative-ordinary-pool-lifecycle-review-preparation-20261016/review.json`.

The completed frontend diagnostic also constrains the next core experiment.
An independent reconstruction of all 112 reviewed raw logs finds zero
ordinary-pool public refactor events in every frontend arm. Each arm records
2,290 EGraph events, all on BJT at eight threads (1,145 per CPU domain).
Those counters cover public refactor/refactor-solve calls; first-factor/repair
pool construction and later certificate reuse have separate scope. Caller
cooperation alone therefore has no demonstrated benefit on the measured
recurring refactor paths. Prioritize sampled EGraph waits/worker imbalance
and measured general parallel admission for currently serial routes before
spending a full competitive campaign on that standalone prototype.
Evidence: `build/completed-frontend-engine-usage-review-preparation-20261016/review.json`.

The first full-28 sampled EGraph diagnostic failed its measurement-coverage
assertion: its guard incorrectly required optional packed factor storage.
Independent raw-log and fresh waveform review qualified all 28 numerical runs,
but absent samples supply no EGraph wait evidence. That failed package remains
preserved. The corrected guard admits both native and packed numerical work
while retaining auxiliary-mode exclusions, 128 warm jobs and one sample per
64 jobs. Native and ASan/UBSan positive controls each execute 640 real unpacked
EGraph refactors, cover crews of two and eight, assert positive sampling and
check a known solution. Their independent review also preserves the inherited
ordinary-pool, certificate, FENV and concurrent-owner assertions. An earlier
control with a missing explicit return after renaming main remains unqualified.

The corrected full-28 diagnostic completed with session 35577 exiting zero.
Its controller and executable were absent before independent review rechecked
all raw counters, input/archive pins, step counts and fresh waveforms. All 28
runs qualify, with 30 post-warmup numerical samples across the two domains.
The sampled EGraph pools each execute 1,145 numeric jobs on the eight-thread
BJT case. No sampled pipeline sleep occurs. Dependency-gate durations are
11.75% and 12.57% of worker elapsed sums on llc32 and llc96; cluster barriers
are 5.72% and 6.34%, and worker wake durations are 1.44% and 1.25%. These
categories overlap, gate time includes productive help, and sampled worker
sums are neither removable wall time nor competitive timing. They do not
establish that unsampled jobs never sleep or that barriers hide no imbalance.

Across the complete grid, ordinary-pool lifetime counters record only two
numeric and six certificate dispatches. Leave the existing sleep policy and
ordinary-pool production implementation unchanged. Prioritize numerical
kernel cost and general structural parallel admission for serial routes;
avoid repeating per-column clock instrumentation that changes low-work engine
selection. Retained optimizations still require complete combined peer timing
and independent confirmation before declaring the fastest-in-SPICE goal met.
Evidence:
`build/egraph-worker-sampling-failed-coverage-review-20261016/review.json`,
`build/egraph-native-worker-sampling-control-review-20261016/review.json`,
`build/egraph-native-worker-grid-review-20261016/review.json`,
`build/egraph-worker-completed-diagnostic-decision-20261016/decision.json`.

Host hardware counters are unavailable: `perf stat -e cycles,instructions --
true` fails with `perf_event_paranoid=4`. The optional diagnostic link object
`scripts/profile_cpu_samples.c` provides bounded process-CPU instruction
sampling on Linux x86-64 without inserting per-column clocks. Build it inside
the Apptainer image, link the object into the diagnostic executable with
`-ldl`, and set `KLS_CPU_SAMPLE_FILE` to a new capture path. After the process
exits, run `python3 scripts/summarize_cpu_samples.py EXACT_EXECUTABLE CAPTURE`.
Preserve that exact unstripped executable and the companion capture `.maps`
file. The reader checks capture completeness and saturation, preserves external
and unresolved samples, and records input hashes. Lock-free integer operations
are the signal handler's only shared-state operations; output occurs after
collection closes. An existing SIGPROF policy or timer disables the sampler.

Positive controls cover disabled sampling, live concurrent work, preserved
rounding modes, an existing signal policy, and PIE/non-PIE address resolution.
Both live controls resolve every captured instruction to the deliberately busy
function (178 and 202 samples). These checks qualify the tool, not solver
performance. ITIMER_PROF is process-directed: signals can coalesce, thread
selection can bias parallel profiles, and system CPU time is included. Do not
infer wall-time savings or precise parallel load balance from sample shares.
Control evidence:
`build/cpu-signal-sampling-positive-control-20261016/review.json`.

The first grid driver stopped on an unimported Path name in its capture-file
check and remains preserved as failed. A fresh corrected full-28 package
completed: session 88633 exited zero, its controller and executable were
absent, and independent review qualified every original waveform, step count,
API contract and complete unsaturated capture. It resolves 223,701 samples.
The qualified parent archive is byte-identical; only the diagnostic link
object is added. Against the completed plain-parent API profile, engine counts
differ only on NEURON's bounded consultations, by 1–17 of 7,740 refactors.
Those differences remain explicit in the review; sample shares do not prove
unchanged timing or selection.

The instruction addresses distinguish two different kinds of waiting.
Of 3,329 samples in `kls_egraph_refactor_pool_worker_main`, 3,325 fall in its
generation/shutdown polling and PAUSE loop, verified against the retained
executable's disassembly. This includes retained certificate crews on serial
numeric paths; it is not the numerical dependency-gate or pipeline-sleep cost
measured above. On llc96 ram2k at eight threads, all 1,076 worker-main samples
are in that idle loop. These counts establish executed idle work, not a
wall-time saving. The next isolated candidate should test job-specific idle
parking for certificate-only crew reuse while preserving bounded numerical
reuse and explicit tuning. Its shared spin budget also needs atomic access;
current workers read it concurrently with controller updates. No production
idle-policy change or speedup claim follows before complete qualification.
Evidence:
`build/cpu-signal-sampled-grid-review-20261016/review.json`,
`build/cpu-signal-sampled-grid-review-20261016/worker-main-hot-pc.txt`,
`build/cpu-signal-sampled-numeric-priority-decision-20261016/decision.json`.

### Cold EGraph startup and certificate-job parking

The parking experiment exposed a separate startup race. Before acquiring its
first generation, a background worker read `worker->shared`, which the
controller can rewrite while preparing the first numerical job. The wait loop
now uses the permanent `&pool->shared` descriptor. Per-job accesses retain their
existing generation publication and completion boundaries. The root fix is
separate from parking policy and changes no numerical kernel or accuracy gate.
The permanent `kls_egraph_numeric_startup` regression repeatedly creates crews
of two and eight, immediately performs a real numerical refactor, checks normal
and transpose solutions of a nonsymmetric matrix, and destroys each crew.
Its unfixed root control reproduces the race; the fixed root source passes
native, ASan/UBSan and memory-access TSan builds. The root componentwise unit is
rebuilt in those controls rather than borrowed from the isolated verdict core.
Review: `build/egraph-startup-root-fix-regression-preparation-20261016/review.json`.

GCC 15's generated IFUNC resolvers initially crashed before TSan initialization:
disassembly shows resolver calls to `__tsan_func_entry/exit` during loader
relocation. An independent minimal IFUNC program reproduces that failure.
`--param=tsan-instrument-func-entry-exit=0` avoids these startup hooks while
retaining memory-access instrumentation. A separate deliberate-race program
still produces the expected TSan report. This compatibility context supplies
memory-race evidence with reduced call-stack instrumentation; it does not
convert the original startup failures into passes. Controls:
`build/gcc-tsan-ifunc-entry-compatibility-controls-20261016/completed.json`.

The isolated parking candidate makes its shared numeric spin budget atomic
and captures a worker-local post-job hint before publishing completion.
Certificate jobs park unless an explicit tuning profile is active; numerical
jobs retain their bounded spin budget. Explicit profiles also bypass automatic
initial/PTS budget boosts. The initial version passed all 88 runtime tests but
failed the directed cold-start TSan check. It remains preserved. The corrected
candidate passes all 88 runtime tests and independent directed controls in
native, ASan/UBSan and memory-access TSan contexts. Together the three contexts
exercise 3,000 certificate jobs, 3,000 loop numeric jobs plus six cold numeric
jobs, concurrent atomic-budget updates, explicit profiles, 24 partial-startup
faults, cleanup and known solutions. These checks qualify concurrency behavior
in their tested scope. Full circuits, complete cost and generality remain
required before any parking promotion; added wakeup cost may outweigh idle
CPU savings for a crew that alternates numerical and certificate work.
Reviews:
`build/certificate-job-idle-parking-startup-safe-runtime-review-preparation-20261016/review.json`,
`build/certificate-job-idle-parking-startup-safe-policy-review-preparation-20261016/review.json`.

The corrected parking candidate subsequently passed all 28 original API/circuit
checks and the complete 238-observation paired cost campaign (14 fresh KLU
references and 224 timed observations). Independent review reconstructed all
476 child-process edges, verified immutable inputs and balanced order, rechecked
all 224 waveforms with the upstream verifier, and recomputed all four paired
passes per circuit/domain/thread. Every observation qualified. The controller
and reviewer both exited zero, with process absence recorded before review.

Parking is **rejected for production**. Candidate/baseline suite geometric mean
linear-solve-time ratios were 0.998652 and 1.012450 on llc32 at one and eight
threads, and 1.000751 and 1.001285 on llc96. Eight-thread BJT was 5.50% slower
on llc32; DAC was 3.13% and 3.22% slower across the two domains. The result
does not justify selecting favorable circuits, changing acceptance rules or
claiming idle CPU time as removable wall time. Preserve the separate startup
race fix and the entire unsuccessful parking experiment. The legacy control
retains documented C11 races: numerical qualification does not make it a
race-free reference. Neither full generality nor fresh peer leadership is
claimed for this rejected candidate.

Completed evidence:
`build/certificate-job-idle-parking-startup-safe-api-review-preparation-20261016/review.json`,
`build/certificate-job-idle-parking-startup-safe-spice-pair-review-preparation-20261016/review.json`,
`build/certificate-job-idle-parking-startup-safe-completed-decision-20261016/decision.json`.

### Eight-row register packet verdict experiment

A separate isolated candidate traverses eight original rows in vector registers
without coefficient packing or cached values. Masked gathers preserve each
row's encounter order; an opaque masked multiply prevents FMA contraction.
The original outward enclosure and lazy tolerance conversion remain. A failed
packet restores the entry floating-point environment and replays the complete
scalar range. Detailed-bound interfaces and current-operator ownership stay
unchanged. Subnormal finite-check lanes also replay: an initial directed test
found that bitwise vector finite checks omitted a scalar denormal flag. The
failed source/test packages remain preserved.

The corrected source passes 5,018 guarded comparisons in each of six contexts
(native, Intel assembler syntax, generic fallback, no contraction, x87, and
ASan/UBSan). Native exact-rational auditing covers 2,216 fixtures, 384,184
comparisons, 98,144 accepted rows and 1,880 observed vector acceptances. Initial
proof review checks the counter instrumentation inverse and unchanged scalar
reference. These controls do not establish enabled-trap or complete solver
qualification.

The complete 2,304-observation cost grid retains all 144 groups per native and
generic build, both orientations and all four numeric fixtures. Independent
review reconstructs balanced order, bounds, immutable inputs and aggregates.
Native candidate/baseline runtime geometric mean is **1.345856** (101 of 144
groups slower); generic fallback is 1.000866. Including shared plan construction
in each measured batch gives 1.180620 and 1.000934; this modeled metric is not
whole-application timing. The register packet candidate is rejected. Do not
rescue it by excluding failed verdicts, zero vectors or unfavorable shapes.
No performance implementation from this experiment enters production.

Evidence:
`build/componentwise-eight-row-register-packet-initial-proof-review-preparation-20261016/review.json`,
`build/componentwise-eight-row-register-packet-complete-cost-review-preparation-20261016/review.json`,
`build/componentwise-eight-row-register-packet-completed-decision-20261016/decision.json`.

### Atomic EGraph idle budget, independent of parking

The production EGraph pool now stores its adjustable worker spin budget in
`atomic_uint`. Construction uses `atomic_init`; idle polling and controller
updates use relaxed atomic accesses. The budget is only a waiting hint. The
existing generation release/acquire still publishes jobs; idle thresholds,
boost conditions, explicit-profile behavior and numerical bodies stay as before.
The immutable ordinary-pool budget is unchanged. This correctness fix does not
promote the rejected certificate parking policy or claim a speed improvement.

The permanent startup regression also updates the budget concurrently with
128 actual refactors across crews of two and eight, then checks known normal
and transpose solutions. A root control reverting only atomic budget access
reproduces the memory-TSan race between its writer and the idle worker. The
fixed source passes native, ASan/UBSan and memory-TSan controls. Root
componentwise sources are rebuilt in those controls; the five retained helper
archive members are independently compared. The previously documented GCC
IFUNC compatibility flag retains memory instrumentation and its deliberate
race positive control. This is tested concurrency evidence, not universal
race-freedom proof.

Fresh root CMake registration contains **43 tests**, including the enhanced
startup regression. All 43 pass in native and optimized ASan/UBSan contexts
(86 total). An initial controller incorrectly expected 45 after the native
43 passed, then stopped before sanitizer work. That failed wrapper is retained;
a separate corrected controller reuses the exact native build, reruns its full
43 and performs a fresh sanitizer build and full43. The historical isolated
44-test suite contains two unpromoted caller-RHS fixtures absent from root;
these counts must not be conflated. No root test is removed or weakened.
Independent review checks source/binary pins, all child edges, raw race reports,
full CTest logs and process absence after actual terminal exits.

Evidence:
`build/egraph-idle-budget-root-atomic-regression-preparation-20261016/completed.json`,
`build/egraph-idle-budget-root-full-runtime-count-corrected-preparation-20261016/completed.json`,
`build/egraph-idle-budget-root-atomic-review-preparation-20261016/review.json`.

### Unused inner L2 snapshots: full cost completed, candidate held

A separate candidate against committed root `2754413` adds six semantic guards
that skip duplicate snapshots used only by the disabled strict-RHS-L2 path
when the locked componentwise policy is active. It preserves allocations,
canonical owned-current-operator data, numeric bodies and public interfaces.
It does not stack the verdict-witness, mapped-worker or parking prototypes.
Native and optimized ASan/UBSan root suites pass all 43 registered tests each.
Poison/page-guard controls preserve the full 960-epoch grid per build/arm:
240 round-to-nearest successes and 720 existing unsupported-rounding
rejections with unchanged outputs. Exact rational review confirms every
successful residual and bit-identical baseline/candidate outputs. The original
fixture's incorrect assumption that all rounding modes must succeed, and its
diagnostics, remain preserved.

All 28 original circuit/domain/thread API settings qualify with 288,650
certificate calls. Fresh root-compatible diagnostic wrappers are confined to
the profiling executable; both plain arms share the committed cacheline Xyce
adapter and dependency archives. Independent review verifies that their link
commands differ only in core archive and output paths.

The complete 238-observation paired campaign then qualifies every observation:
14 fresh KLU references plus 224 timed observations, four balanced passes per
circuit/domain/thread. The actual controller exit is zero, all 238 child
processes are absent, and independent review reconstructs all 476 child edges,
raw logs, current input hashes, fresh waveform checks, step counts and ratios.
Candidate/baseline suite geometric mean linear-solve-time ratios are 0.989565
and 0.981619 on llc32 at one and eight threads, and 0.992559 and 0.997187 on
llc96. Twenty-seven of 28 cell medians improve, but llc96/eight-thread mem_plus
regresses: its four ratios are 0.993008, 1.032280, 1.055488 and 1.067632,
with median 1.043884. That case remains included in every aggregate.

The candidate is **held**, with no production promotion, broader-corpus or
peer-leadership claim. A matched full-28-per-arm API diagnostic is prepared to
compare engine counts under identical instrumentation. Such diagnostics can
change bounded timing consultations and cannot retrospectively prove the cause
of a plain timing regression. Any independent timing confirmation must use a
fresh complete balanced campaign; no favorable subset or case-specific selector
is authorized by these results.

Evidence:
`build/componentwise-inner-l2-copy-elision-root-runtime-review-preparation-20261017/review.json`,
`build/componentwise-inner-l2-copy-elision-root-fp-contract-review-preparation-20261016/review.json`,
`build/componentwise-inner-l2-copy-elision-root-api-review-preparation-20261017/review.json`,
`build/componentwise-inner-l2-copy-elision-root-spice-pair-review-preparation-20261017/review.json`,
`build/componentwise-inner-l2-copy-elision-root-completed-cost-decision-20261017/decision.json`.

The next isolated row-vector helper remains outside production. Raw controls
cover 195,842 comparisons in each of six native/generic/feature-subset/sanitizer
contexts (1,175,052 total), with explicit within-packet collision replay,
cross-packet repeats, signed indices, tails, IEEE values and protected pages.
Finite/signed-zero output bits and FP status match the production-contracted
scalar loop. Native and minimal F/CD/VL contexts each record 78,304 differing
written NaN payloads; no bit-identical NaN propagation is claimed. Source
instrumentation has an exact inverse, but real-edge construction, full numeric
integration, unmasked traps, performance and generality remain unqualified.
Raw review:
`build/snb-narrow-row-vector-first-proof-review-preparation-20261017/review.json`.

The matched full-28-per-arm diagnostic subsequently completed all 56 settings.
Both arms use one freshly compiled six-function wrapper object, the same
cacheline adapter and the same dependencies. Independent review reconstructs
all 118 compiler/linker/circuit child edges, source/binary pins, full-grid order,
raw counters and fresh waveforms. Every setting qualifies. Reported call counts
for mem_plus are identical between arms. Differences occur only in NEURON's
mapped/native-KLU path counts, by 1–14 of 7,740 refactors per setting. These
instrumented observations do not explain or retrospectively attribute the plain
4.39% regression. Copy removal remains held; any confirmation must retain the
whole original paired grid in a fresh campaign. Matched review:
`build/componentwise-inner-l2-copy-elision-root-matched-api-review-preparation-20261017/review.json`.

The separate row-vector prototype now also passes 2,728,320 whole-edge fixtures
in each of native, generic and optimized ASan/UBSan builds (8,184,960 total).
The source capsule extracts the current root edge/supernode/block declarations,
unchanged dense-product helper and complete edge arithmetic. Reversing the one
new narrow-row helper call restores the original edge body byte-for-byte. The
full original semantic grid is retained, with extra FTZ/DAZ cases at row heights
8, 9, 15 and 16 to exercise vectors and tails. Finite/infinity/signed-zero bits,
nonfinite quiet/signaling class, raw MXCSR/C exception flags, protected pages
and producer bytes agree. Native arithmetic differs in 103,860,144 NaN payload
words across 230,208 fixtures; those differences are counted, not treated as
bit identity. Generic and sanitizer comparisons report none. This establishes
raw bounded-edge arithmetic, not real-matrix construction, unmasked traps,
no-contract behavior or full solver correctness. The original unfiltered
201,600-observation cost grid is prepared/running separately; no kernel speed
claim follows before its actual terminal and independent full-grid review.
Whole-edge review:
`build/snb-narrow-row-vector-edge-review-preparation-20261017/review.json`.

The row-vector prototype's complete 201,600-observation cost run subsequently
finished with actual zero exit and all eight compiler/timing child edges closed
and absent. Independent review checks the exact original 12,600 shapes, four
balanced pairs, identical 112-byte root block types, raw flags/bytes, source and
binary pins, and every aggregate. Native candidate/baseline geometric mean cost
is 1.275525, with 10,293 of 12,600 shape medians slower; generic cost is 1.000627.
The entire preexisting narrow branch (7,560 shapes, without exclusions) costs
1.498855 natively; the unchanged wide control costs 1.001342. This general
helper is rejected. Correctness controls do not justify a production promotion,
a favorable-shape selector, a full solver claim or a peer-leadership claim.
Review and preserved decision:
`build/snb-narrow-row-vector-complete-cost-review-preparation-20261017/review.json`,
`build/snb-narrow-row-vector-completed-cost-rejection-20261017/decision.json`.

The six unused-inner-L2-copy guards remain held. A fresh independent complete
238-observation confirmation is prepared using the same qualified plain
binaries, reversed domain/circuit order and a shifted balanced arm phase.
It retains every circuit, both thread counts, fresh KLU references, default
componentwise verification and all failures. The first campaign's unexplained
mem_plus regression remains evidence regardless of the confirmation result;
no result is pooled or replaced. Broader-corpus and fresh peer qualification
remain outstanding before any production or fastest-in-SPICE claim.

The independent copy-guard confirmation is now running in its own package:
`build/componentwise-inner-l2-copy-elision-root-confirmation-preparation-20261017`.
The original full-110-matrix/1,760-observation broad controller, linker and
reviewer are prepared separately with exact root-qualified core archives and
one unchanged general driver. They are unexecuted; the closed confirmation
review and an advance-stage decision are prerequisites. No result follows
from preparing those sources.

The next independent core integration candidate rebases the earlier qualified
coefficient-domain/tolerance-cache/verdict work onto root `2754413`. The exact
certificate sources are preserved; the root EGraph declaration remains atomic
and its startup/publication and idle-budget implementation remain byte-exact.
This changes five solver source files (including two new private includes),
plus two existing test adapters: the caller fault hook follows the actual new
verdict boundary and internal eligibility assertions use its validated domain
state. Existing fault/correction/worker assertions and the root 43-test suite
remain. There are no unused-L2-copy guards, mapped specialization or parking
changes in this source. An initial source-only preflight incorrectly assumed
one idle-budget field; the preserved correction scopes the atomic declaration
to EGraph and leaves the ordinary pool unchanged. No build or numerical test
ran in that failed preparation.

Corrected source capsule:
`build/componentwise-verdict-witness-root-rebase-source-corrected-preparation-20261017/preparation.json`.
A separate full native/optimized-ASanUBSan 43+43 runtime controller is prepared,
not executed, under
`build/componentwise-verdict-witness-root-full-runtime-preparation-20261017`.
Historical isolated correctness, broad and peer results motivate this rebase;
they do not qualify the current-root combination or establish a speedup over
root. Fresh root runtime, API, complete timing, generality and peer checks
remain required. No concurrent build or proof may run during the live
uninstrumented confirmation.

The root certificate rebase is now also committed on the isolated branch
`work/certificate-verdict-root-20261017`, based on `2754413`, in two chunks:
`e3cc481` carries the validated binary64 degree/underflow machinery and
range-local tolerance reuse; `24dc56f` adds the private affirmative verdict
path when the caller does not consume an accepted bound. Every source,
public header and existing test file in its clean sparse worktree matches
the pinned corrected capsule. The branch retains both root EGraph fixes.
Main numerical code remains unchanged, and these commits are not a current-root
runtime, performance or generality qualification. Runtime, API, full-238 and
full-1,760 controllers/reviewers are prepared separately and remain unexecuted.
Branch/source receipt:
`build/componentwise-verdict-witness-root-rebase-source-corrected-preparation-20261017/integration-branch.json`.

A read-only SubtreeLU upstream/API check finds installed clean revision
`9930fd2cadff043f44bd0c9b48a63403f6ed218e` equal to a fresh official main-ref
lookup. Xyce supplies zero-based CSR and null permutation arguments for
analysis, then factorizes, refactorizes changed values and solves in the
sequence documented by the official repository. This limited source check
finds no newer main revision or sequence mismatch; it does not diagnose or
resolve the prior two MOS eight-thread timeouts. They and the incomplete peer
comparison remain unqualified. Evidence:
`build/subtreelu-upstream-and-api-source-check-20261017/source-check.json`.
Official guidance: https://github.com/THU-numbda/SubtreeLU#usage-of-the-library.

The independent copy-guard confirmation completed all 238 observations with
actual zero exit, controller absence and all 476 numerical child edges closed
and absent. Independent review reconstructs the reversed domain/circuit order,
shifted balanced arm phase, all fresh waveforms, step counts, raw logs and input
pins. Every observation qualifies. Candidate/baseline suite geometric means
are 0.991751 and 0.986136 on llc32 at one/eight threads, and 0.991316 and
0.989396 on llc96. Twenty-five of 28 cell medians improve. The original
llc96/eight-thread mem_plus regression does not recur (confirmation 0.995417),
but mux regresses by 0.53% on llc32/one thread, 1.12% on llc32/eight threads
and 0.66% on llc96/eight threads. The first campaign, these three regressions,
and every raw pair remain included separately; no campaign or case results
are pooled or replaced. No cell is slower in both campaigns. This does not
explain the original regression or establish absence of regressions.

Both complete campaigns qualify every observation and all eight separate suite
aggregates improve, so the whole unchanged candidate advances only to broader
validation. It remains held from production. The exact four-file/six-guard
source is committed independently as `2f051ec` on
`work/componentwise-l2-copy-guards-20261017`; it is not combined with the
certificate branch. Its full-110-matrix controller remains prepared, not run.
The independent root certificate native/ASanUBSan suite is executing first.
Evidence:
`build/componentwise-inner-l2-copy-elision-root-confirmation-review-preparation-20261017/review.json`,
`build/componentwise-inner-l2-copy-elision-root-two-campaign-decision-20261017/decision.json`,
`build/componentwise-inner-l2-copy-elision-root-two-campaign-decision-20261017/integration-branch.json`.

### General row-admission source finding, 2026-10-17

Read-only inspection identifies a remaining dimension-based automatic row gate:
`kls_auto_row_refactor_high_work_candidate` requires 300,000–800,000 rows and
at least 1e8 dependency work. A second preparation cutoff cites a measured
setup on rajat03 and the sizes of prior adopters. These cutoffs are not a
representation-validity argument and may prevent the existing paired
row/column tournament from examining other profitable structures. No speedup
from removing them is established.

The preparation path needs more than a changed predicate. Both general and
parallel pattern builders may free retained row metadata before constructing
the replacement, and neither accepts a preparation budget/deadline. Moving
construction later also requires checking packed/native factor authority
before seeding from retained value-pointer mirrors. A general replacement
must price complete preparation and publication over the caller-declared
lifecycle, preserve usable incumbent state when a bounded attempt declines,
and retain the existing measured refactor-plus-solve acceptance and re-audit.
Cancellation/OOM and factor freshness need directed controls before broader
admission. No numerical source or live timing input changes in this audit.
Evidence: `build/row-admission-general-resource-source-audit-20261017/audit.json`.

### Complete root certificate cost gate, 2026-10-17

The certificate branch (`e3cc481`, `24dc56f`) completes the fresh root-native
43-test suite, optimized ASan/UBSan 43-test suite, and all 28 API settings.
Its complete paired circuit campaign and independent review are now closed:
all 238 observations qualify, including 14 fresh KLU references and 224 timed
runs covering every circuit/domain/thread setting with four balanced pairs.
Candidate/baseline suite geometric means are 0.718491 and 0.835299 on llc32
at one/eight threads, and 0.702081 and 0.832481 on llc96. These are comparisons
with root2754413, not new gains over the earlier isolated verdict candidate
or comparisons with either peer solver.

Twenty-seven of 28 cell medians improve. NEURON/rallpack3 on llc96/eight
threads consistently regresses: median 1.034608, individual pair ratios
1.036760, 1.035606, 1.028744 and 1.033611. Its llc32/eight-thread result
is mixed (0.891502, 1.026211, 1.023063, 0.909794). Every pair remains in the
decision. Production promotion and broader qualification remain held while
a matched full-grid diagnostic examines certificate execution choices.
Cold calibration is a hypothesis, not an established cause. Verification
remains enabled; no circuit, dimension or measured-result selector is added.
Evidence:
`build/componentwise-verdict-witness-root-spice-pair-review-preparation-20261017/review.json`,
`build/componentwise-verdict-witness-root-completed-cost-decision-20261017/decision.json`.

### Identity-layout snapshot constructor controls, 2026-10-17

A separate caller-order snapshot prototype detects an identity value-source
map by checking every validated slot. Fresh owned caller bytes can then be
used as canonical bytes without indirect reads. Nonidentity maps retain
their existing interpretation. This prototype is separate from the verdict
and inner-L2-copy branches; its parent still needs fresh complete root
qualification.

The constructor fixture covers both orientations, identity/nonidentity slot
maps, all four rounding modes, three accuracy policies, fresh value epochs,
in-place map changes, allocation failures and canonical allocation fallback.
The first run failed its final expected coverage count: fallback checks had
incorrectly been counted as mapped. All preceding semantic assertions passed.
That failed package is preserved, including the unavailable original tool
terminal receipt; no successful controller exit is inferred. A corrected
preparation then stopped before compilation because its scratch path still
named the earlier run. It is also retained. A fresh corrected package passes
all six parent/candidate native, generic and optimized ASan/UBSan contexts;
independent review checks all 24 child edges, source pins and logs. Each
candidate context checks 144 identity and 144 nonidentity snapshots, with
32 caller-order and 256 canonical interpretations. This is constructor
qualification only; no complete-solver or speedup claim follows.
The exact source is preserved on `work/caller-order-identity-layout-root-20261017`
as separate parent snapshot commit `fea844d` and identity-access commit
`d4089e9`. Neither is merged into production.
Evidence:
`build/caller-order-snapshot-identity-layout-root-controls-preparation-20261017/failure-review.json`,
`build/caller-order-snapshot-identity-layout-root-coverage-scratch-corrected-controls-review-preparation-20261017/review.json`.

### Matched certificate plan diagnostic and bounded re-audit, 2026-10-17

All 56 matched baseline/candidate diagnostic runs complete and independently
qualify. The same seven-function wrapper object and shared dependencies are
used in both arms. Two exact, separately named private plan layouts are read
with `memcpy` after the certificate timer stops; the resolved private verdict
symbol determines the layout. The review reconstructs all 118 child edges,
all waveforms and step counts, raw counters and source/layout pins.

NEURON on llc96/eight threads records 7,739 parallel-choice states and two
undecided states in both arms. On llc32/eight threads, the baseline records
7,739 parallel-choice states whereas the candidate records 7,739 serial-choice
states, with two undecided states each. Those post-call states describe the
retained choice, not necessarily the execution of the call that just decided
it. The llc96/eight-thread mapped/KLU refactor counts also differ: 7,718/22
for baseline and 7,701/39 for candidate. Smaller mapped/KLU count differences
occur in the other NEURON settings. Other circuits retain matching numerical
refactor-path counts; llc32/eight-thread BJT has three additional candidate
certificate calls/misses. Instrumentation may affect adaptive decisions.
These observations do not establish the cause of the preceding plain timing
regression, and do not justify a NEURON-specific selector.

The next isolated candidate, `76b0aa3` on
`work/certificate-warm-selection-root-20261017`, adds general geometric
re-audits to the qualified certificate parent `24dc56f`. The initial choice
expires after three accepted uses, then after six, twelve, and successively
doubled use intervals. Each re-audit retains the existing three balanced
serial/parallel probe pairs and proof-disagreement checks. Rejected proofs
do not advance the use counter. Intervals saturate without wrapping, and
thread/executor/output-contract changes reset timing history. Serial-only
requests continue to avoid worker probes. No matrix name, dimension window,
case-specific threshold or verification relaxation is introduced.

The two-file source change and directed controls are committed separately.
Parent and candidate pass all six native, generic and optimized ASan/UBSan
plan contexts, independently reviewed. The new fixture checks a deterministic
cold/warm cost reversal using real proof arithmetic, changing current values,
fresh rejection, a logarithmic probe budget, counter saturation and declined
executor fallback. This is private plan qualification, not a speedup claim.
Fresh full solver native/ASanUBSan suites are the next gate. The first full
configure stopped because a Git archive leaves dependency gitlinks empty;
that package is retained. Its fresh replacement materializes the four exact
existing dependency commits and pins their archived contents. No dependency
version changes or solver results are inferred from the failed configure.
Evidence:
`build/componentwise-verdict-witness-root-matched-plan-review-20261017/review.json`,
`build/componentwise-verdict-witness-root-matched-plan-review-20261017/decision.json`,
`build/certificate-warm-selection-root-controls-review-preparation-20261017/review.json`.

The corrected full runtime subsequently completes with actual controller exit
zero and all recorded children absent. Independent review confirms all 43
native and 43 optimized ASan/UBSan tests in the unchanged root test grid,
including the new warm-selection fixture assertions. It also reconstructs
the committed candidate source and the archived contents of the four pinned
dependency submodules. The source is still held from production; the fresh
28-setting API gate and the separate full238 parent/child cost gate follow.
Evidence:
`build/certificate-warm-selection-root-full-runtime-review-preparation-20261017/review.json`.

The fresh API gate also completes and is independently reviewed: all 28
settings qualify, with 288,647 outer certificate calls through the private
verdict interface, zero legacy outer calls, and complete waveform/step checks.
Plain parent/child links differ only in the core archive and output paths;
both retain the same committed adapter and shared dependency archives, with
no diagnostic wrappers in either plain executable. The exact whole candidate
advances to a fresh full238 paired cost campaign against parent `24dc56f`.
This campaign is separate from the earlier root2754413 comparison. No cost
result, production promotion, generality or peer superiority is inferred.
Evidence:
`build/certificate-warm-selection-root-api-review-preparation-20261017/review.json`,
`build/certificate-warm-selection-root-spice-pair-preparation-20261017/stage-decision.json`.

### End-to-end row preparation and late equality source work, 2026-10-17

While the isolated re-audit candidate runs its complete paired timing grid,
read-only inspection refines the row-admission blocker. The lower-bound
estimator itself builds refactor maps and native LU pointer caches, allocates
row counts and scans the factors before row-pattern construction starts.
The row destructor also owns lean/compact input maps and scheduling state;
`compact_match_offdiag_user_pos`, for example, participates in direct caller
value admission and off-diagonal refresh. Parallel construction publishes a
stack job and the real solver into its retained worker pool. Its finisher
publishes base arrays before all auxiliary construction has completed. A
builder-only timeout or shallow solver clone therefore does not establish
bounded preparation with a usable incumbent. Dimension gates remain unchanged.
Source-pinned evidence:
`build/row-admission-end-to-end-ownership-audit-20261017/audit.json`.

A separate, smaller prototype uses the existing owned original-system snapshot
to re-examine exact equality after the raw unchanged-input cache declines its
first changed update. Failed scans back off geometrically. A positive result
compares every current caller value byte and the complete saved layout, and
requires the preceding successful public numeric epoch with no intervening
numeric/recovery entry. It dispatches the existing unchanged-factor executor;
current RHS/solution verification and fresh snapshot publication remain.
No additional per-update value snapshot is introduced. The prototype captures
x87 control and SSE control modes and requires a match; sticky exception flags
are excluded. Platforms without complete control capture retain ordinary
refactoring. It assumes canonical owned snapshot values and needs an explicit
adapter before combining with a caller-order snapshot variant.

The code and directed controls are committed as `356d27d` and `4394fdf` on
`work/late-unchanged-owned-snapshot-root-20261017`, based independently on
root2754413. Controls cover a first changed update followed by equality, fresh
RHSs, combined calls, CSR/CSC, both orientations, scaling, failed updates,
numeric/recovery epochs, control-mode changes and counter saturation. This is
source-only work: neither commit has been compiled or tested. The prepared
six-context comparator runner waits for the live paired timing and its
independent review to close. No compiler, test, profiler or numeric proof is
launched alongside that timing; its frozen binaries remain unchanged.
Evidence:
`build/late-unchanged-owned-snapshot-root-source-preparation-20261017/preparation.json`,
`build/late-unchanged-owned-snapshot-root-controls-preparation-20261017/run.py`.

### Complete re-audit cost and late-equality comparator controls, 2026-10-17

The geometric re-audit campaign closes with actual controller exit zero and
all 476 recorded child edges absent. Independent review qualifies every one
of the original 238 observations and reconstructs the complete balanced grid,
fresh references, waveforms, steps, source pins and raw ratios. Against parent
`24dc56f`, candidate `76b0aa3` has suite geometric means 0.992599/0.999433 on
llc96 at one/eight threads and 0.991717/0.985752 on llc32. Nineteen of 28 cell
medians improve. Every other cell and every raw pair remains included.

NEURON on llc96/eight threads is near parity (median 0.999354; pairs 1.008233,
0.995817, 1.002891, 0.995535). On llc32/eight threads its median is 0.974221,
with mixed pairs 0.882435, 0.940542, 1.123314 and 1.007900. This single
parent/child campaign does not establish a remedy for the earlier 3.46%
root2754413 regression. No ratios are multiplied or pooled with the earlier
root comparison. Production promotion remains held; independent confirmation,
direct root comparison and broader qualification are incomplete. The next
priority is qualification of the independent late-equality optimization.
Evidence:
`build/certificate-warm-selection-root-spice-pair-review-preparation-20261017/review.json`,
`build/certificate-warm-selection-root-completed-cost-decision-20261017/decision.json`.

After that timing and its review close, all six parent/late-equality standalone
comparator contexts pass (native, generic and optimized ASan/UBSan).
Independent review confirms their complete source archives, six-file delta,
compiler inputs, executable hashes, logs and all 24 child edges. They exercise
exact bit equality, slot/layout changes and rounding/exception preservation;
they do not establish numeric-epoch integration or a performance benefit.
The fresh complete solver controller subsequently passes its 43 native tests
and starts the optimized sanitizer build. Complete runtime qualification
remains pending that suite and independent review.
Evidence:
`build/late-unchanged-owned-snapshot-root-controls-review-preparation-20261017/review.json`.

The late-equality controller subsequently completes both complete solver
suites with actual exit zero and all recorded children absent. Independent
review confirms all 43 native and 43 optimized ASan/UBSan tests, the unchanged
root test grid, committed six-file source delta and exact dependency archives.
This includes the public-API reuse/RHS/layout/epoch controls. Separate native
archive inspection finds both `kls_componentwise.c.o` and
`kls_accuracy_certificate.c.o` byte-identical to their root2754413 objects.
That is evidence about those arithmetic implementations, not a substitute for
caller freshness or runtime checks. The fresh full28 API/waveform gate and
full238 paired cost gate remain the next stages; no performance, generality
or peer leadership follows from these runtime results.
Evidence:
`build/late-unchanged-owned-snapshot-root-full-runtime-review-preparation-20261017/review.json`,
`build/late-unchanged-owned-snapshot-root-certificate-object-identity-20261017/identity.json`.

### Late-equality application gate and next validation stage, 2026-10-17

The fresh full28 API gate qualifies every setting and all 288650 certificate
calls. Independent review closes the actual controller and children, waveform
and step checks, source pins and matched links. However, every setting records
zero unchanged-factor reuse. Every public refactor is a combined refactor/solve.
The committed Xyce adapter explains this: `updateValues_` compares every value
byte, and `refactorSolve_` already calls `kls_solve` directly when values are
unchanged. Changed values alone reach the combined refactor API.

Candidate `4394fdf` is therefore shelved for the current SPICE goal. Its source
and correctness evidence remain available; the prepared full238 timing is not
launched. Removing the existing adapter optimization would manufacture an
opportunity rather than improve the application. No performance benefit or
overhead has been measured for this prototype. Evidence:
`build/late-unchanged-owned-snapshot-root-api-review-preparation-20261017/review.json`,
`build/late-unchanged-owned-snapshot-root-no-application-opportunity-decision-20261017/decision.json`.

The next executable stage is broader qualification of the independent
inner-L2 copy-elision candidate `2f051ec`. Both complete SPICE campaigns have
already qualified all observations; all eight separate suite aggregates
improve, while their cell regressions remain recorded. Broad binaries now
link the exact qualified root/candidate archives with an identical driver and
identical dependencies inside Apptainer. This is preparation for the complete
110-matrix control, not broader qualification or production promotion.
Evidence:
`build/componentwise-inner-l2-copy-elision-root-two-campaign-decision-20261017/decision.json`,
`build/componentwise-inner-l2-copy-elision-root-broad-build-review-preparation-20261017/review.json`.

An independent full238 confirmation is also prepared for root certificate
candidate `24dc56f`. It uses the same qualified plain executables, all seven
circuits, both CPU domains, threads one/eight, four balanced pairs and fresh
KLU references. Domain/circuit order and balanced arm phase complement the
first root campaign. Its complete results must remain separate from that
campaign and from the `76b0aa3` parent/child experiment. In particular, the
first campaign's NEURON llc96/eight-thread regression remains unresolved.
Both freezing and execution require closure of the current copy-elision
broad run and its independent review. Only controller/reviewer source has
been prepared and parsed; no confirmation observations exist yet. Evidence:
`build/componentwise-verdict-witness-root-confirmation-preparation-20261017/source-preparation.json`.

### Complete independent root certificate confirmation, 2026-10-17

The fresh confirmation closes with actual controller exit zero, all 476 child
edges closed and the controller absent. Independent review subsequently closes
with exit zero and qualifies all 238 observations, including fresh waveform
and step checks and the complete unique balanced grid. The same qualified
plain executable contents, original netlists and verifier are used in both
root campaigns. All 28 confirmation cell medians improve. Suite ratios are
0.723104/0.833566 on llc32 at one/eight threads and 0.706363/0.819013 on llc96.

The original NEURON llc96/eight-thread median 1.034608, with all four pairs
slower, remains retained. Confirmation median is 0.907374, with pairs
0.898203, 0.900756, 1.032787 and 0.913991. That regression is not reproduced
at the median and remains unexplained. NEURON llc32/eight threads remains
mixed, with confirmation median 0.967237. No cell median is slower in both
campaigns; all eight separate suite aggregates improve. Neither raw pairs nor
campaign aggregates are pooled, replaced or multiplied with older experiments.

Candidate `24dc56f` has independently confirmed application benefit and is
eligible for broader validation. Production promotion, broader qualification
and current CKTSO/SubtreeLU leadership remain unproven. Short independent
arithmetic/FP controls for the dyadic child precede choosing the whole candidate
for expensive broader/peer validation; there is no source stacking or automatic
promotion. Evidence:
`build/componentwise-verdict-witness-root-confirmation-review-preparation-20261017/review.json`,
`build/componentwise-verdict-witness-root-two-campaign-decision-20261017/decision.json`.

### Dyadic private and post-loop controls completed, 2026-10-17

Candidate `d786a14` completes all six parent/candidate private control
contexts (native, generic and optimized ASan/UBSan). Actual controller exit is
zero and all 24 child edges close; independent review reconstructs the full
committed source archives and exact two-file delta. These are private fixture
controls, not a complete solver suite.

The separate parent post-loop probe subsequently completes native, generic
and sanitizer contexts at O0/O3, with actual exit zero and all 26 child edges
closed. Independent review confirms all 39813120 comparisons: each context
has 6635520 comparisons, 338688 newly admitted tuples and 6296832 declines.
Verdicts, bound bytes and x87/SSE sticky flags match the extracted parent
post-loop statements throughout. All six contexts produce the same digest.
Independent rational checks also verify the parent-acceptance inequalities
under the stated normal-rounding and ULP-step lemmas. This is not a formal
proof of the complete solver.

These controls assume an already-ready tolerance cache. Whole-row/range and
original-conversion controls, complete solver integration, application
opportunity and cost remain separate. The next stage is complete whole-range
microcost against parent `24dc56f`, preserving all fixture groups and both
compiler profiles. No production promotion or speed benefit follows from
the arithmetic controls. Evidence:
`build/certificate-dyadic-margin-root-controls-review-preparation-20261017/review.json`,
`build/certificate-dyadic-margin-root-post-proof-review-preparation-20261017/review.json`.

### Complete dyadic microcost, mixed result, 2026-10-17

Both native/generic binaries compile before either timing starts. Actual
controller exit is zero, all eight child edges close, and independent review
reconstructs all 2304 observations and the complete balanced grid. Baseline
uses the full committed parent verdict path, not its earlier detailed-bound
predecessor. All original 144 groups per build, iteration counts and timed
calls remain; symmetric entry-CSR recording and stronger untimed bound/FP
checks (precision clear/set) are explicit benchmark additions.

Candidate/parent geometric means are 0.897362 native and 0.869610 generic.
This is mixed: 79/144 native and 83/144 generic medians regress; 45 and 53
respectively exceed 1.02. The maximum medians are 1.284745 and 1.270780.
Every group and every raw pair remains included. All timed entries have SSE
precision sticky set naturally by the retained timing clock; cold-entry
performance is not established. No favorable fixture, shape, size,
orientation or compiler subset is used to rescue or select the implementation.

The unchanged whole candidate may advance to independent range/mixed-order,
raw-tolerance and zero-denominator controls, followed by complete solver/API
and paired application qualification. Aggregate microcost gains do not
override those gates or establish application benefit. No numerical source
is promoted. Evidence:
`build/certificate-dyadic-margin-root-full-grid-micro-review-preparation-20261017/review.json`,
`build/certificate-dyadic-margin-root-completed-micro-decision-20261017/decision.json`.

### Dyadic raw range and mixed-order controls completed, 2026-10-17

Fresh parent/child private-verdict probes complete all twelve compiler
contexts, including O0/O1/O2/O3/Og/Os, generic/native, x87 arithmetic,
Intel assembly, separate rounding and no-integer128 builds. Actual exit is
zero and all 96 child edges close. Independent review confirms 428752 range
comparisons in each context, totaling 5145024 matching status, bound-byte,
MXCSR and x87 sticky observations. Each context retains all 372 fixtures,
including the existing ordered, mixed, same-bin and zero-denominator families.

Raw coverage includes 57 actual 80-bit tolerances, all rounding directions,
x87 precision 24/53/64, and all 64 prior sticky combinations in each unit
separately plus both units all-set. Independent readers reconstruct complete
compact traces and actual argument bytes, and verify 443328 supported
diagonal error bounds using exact fractions. Unsupported environments and
raw tolerances establish parent/child equivalence, not new accuracy support.
This remains private-range qualification, separate from public correction,
complete solver, generality and application cost gates.

The exact committed child and four unchanged dependency archives are then
materialized for the fresh complete 43-native/43-sanitizer runtime gate.
That gate is running; no runtime result or production promotion is inferred.
Evidence:
`build/certificate-dyadic-margin-root-raw-range-review-preparation-20261017/review.json`,
`build/certificate-dyadic-margin-root-full-runtime-preparation-20261017/full-source-preparation.json`.

The complete dyadic solver gate subsequently closes with actual exit zero
and all twelve child edges absent. Independent review confirms every one
of the 43 native and 43 optimized ASan/UBSan tests, the unchanged root test
grid, exact committed source and four unchanged dependency archives. Both
builds retain the existing accuracy policy and root worker race fixes.
This is complete runtime qualification, not broader matrix or performance
qualification. The full28 Xyce API gate is now launched, using the current
child source headers, qualified parent core, shared dependency/adapter objects
and fresh KLU references from the completed root confirmation. Seven wrappers
are confined to the diagnostic executable; both plain arms have none.
Evidence:
`build/certificate-dyadic-margin-root-full-runtime-review-preparation-20261017/review.json`.

The full28 dyadic API gate subsequently closes with actual exit zero, its
build controller absent and all 64 recorded build/numeric child edges closed.
Independent review qualifies every setting, fresh-reference waveform/step
checks and all 288650 certificate calls (all verdict, zero legacy). Matched
plain links differ only in core archive and output paths; both have no
diagnostic wrappers. This establishes application numerical/API qualification,
not a speed benefit from the new helper.

The unchanged whole child now advances to the complete238 paired cost grid
against parent `24dc56f`: all seven circuits, both CPU domains, one/eight
threads, four balanced pairs per cell and fourteen fresh KLU references.
No microbenchmark group or case is used to select a different implementation.
No compiler, test or profiler runs alongside uninstrumented timing; assessment
waits for actual terminal, complete child closure and independent review.
Evidence:
`build/certificate-dyadic-margin-root-api-review-preparation-20261017/review.json`,
`build/certificate-dyadic-margin-root-spice-pair-preparation-20261017/stage-decision.json`.

### Dyadic whole-application cost rejected, 2026-10-17

The complete238 campaign closes with actual exit zero and all 476 recorded
child edges absent. Independent review qualifies every observation and
reconstructs the original complete balanced grid, fresh references, waveforms,
steps, pins and all raw pairs. Against parent `24dc56f`, child `d786a14` has
suite ratios 1.033170/1.032476 on llc96 at one/eight threads and
1.029167/1.031749 on llc32. Only 7/28 cell medians improve. The full results,
including every regression, remain preserved.

The whole child is rejected for production and does not advance to broad
validation. No favorable case, fixture, compiler or thread subset rescues it.
Numerical/FP qualification remains valid evidence for that scope, but its
warm-clock microcost gains did not transfer to the application. The cause is
not established. Parent `24dc56f`, whose application benefit was independently
confirmed, remains the best qualified certificate candidate. A fresh complete
CKTSO/SubtreeLU comparison is the next priority; neither older peer ranks nor
products of unrelated campaign ratios establish current leadership.
Evidence:
`build/certificate-dyadic-margin-root-spice-pair-review-preparation-20261017/review.json`,
`build/certificate-dyadic-margin-root-completed-cost-decision-20261017/decision.json`.

### Dyadic certificate margin prototype, source only, 2026-10-17

Isolated commit `d786a14`, based on `24dc56f`, adds an integer sufficient
acceptance test after the unchanged current-row residual/magnitude loop.
With dyadic `t <= tolerance64 <= 1`, it requires coefficient `c <= t/8`
and residual magnitude at most `t/4` times the lower power of two of the
computed magnitude. The proposed inherited enclosure then gives an error
ratio at most `2*t/3`. A declined test keeps the byte-identical ordinary
post-loop source. No circuit, thread, dimension or compiler selector is added.

The new test requires the tolerance cache already initialized at its original
conversion point. It also requires range-entry SSE inexact already sticky,
computed magnitude at least `2^-972`, and a normal or exactly zero residual.
These arithmetic restrictions are intended to preserve the skipped enclosure's
FP side effects. They remain proof obligations, not a claim of qualification.
Directed source controls cover exact threshold adjacency, degree/exponent
boundaries, subnormal/NaN declines and all 64 SSE sticky combinations.

Neither candidate nor controls have been compiled or executed. The prepared
native/generic/optimized-sanitizer runner waits for the current broad campaign
and its review, then independent root confirmation and its review. Exact
arithmetic, parent post-loop equivalence, all rounding/raw-tolerance controls
and complete cost evidence remain required. This prototype is not evidence of
a remedy for the root candidate's NEURON regression. No L2-copy guards,
geometric re-audit or parking changes are stacked. Source evidence:
`build/certificate-dyadic-margin-root-controls-preparation-20261017/source-review.json`,
`build/certificate-dyadic-margin-root-controls-preparation-20261017/source-only-policy.json`.

Additional post-loop equivalence controls are prepared, also unexecuted.
After the six private controls and independent review close, their source
preparer extracts the exact parent post-loop statements and checks they match
the candidate fallback. The probe covers all 64 SSE sticky combinations,
degree/exponent/tolerance and adjacent residual boundaries, including NaNs
and subnormals: 6635520 comparisons in each of native, generic and sanitizer
contexts at O0/O3. It compares verdict, bound bytes and x87/SSE sticky flags.
An independent rational check bounds the parent's rounded enclosure and
requires every newly admitted tuple to pass that parent enclosure. No C
compilation, algebra execution or probe result exists yet. These controls
concern an already-ready tolerance cache; complete rows, first conversion,
worker integration and full cost remain separate gates. Evidence:
`build/certificate-dyadic-margin-root-post-proof-preparation-20261017/source-only-policy.json`.

### Compute priority changed; broad attempt preserved incomplete, 2026-10-17

The copy-elision broad attempt is explicitly cancelled to prioritize the
substantive root certificate confirmation. The decision uses completed,
independently reviewed evidence: the first root certificate comparison has
suite ratios 0.702--0.835 and 27/28 improved medians, with its NEURON
regression unresolved; the two complete copy-elision comparisons have roughly
1--2% suite gains. No partial broad timing or qualification is assessed to
make this scheduling decision. The first certificate comparison is reviewed,
not independently repeat-confirmed; its complete confirmation is the next job.

Actual broad-controller exit is 130. Independent closure review verifies all
recorded children and observed descendants absent, no frozen native process
remaining, and all frozen campaign inputs unchanged. All 231 completed
observations, 463 ledger edges and the interrupted child remain preserved.
The last start has external cancellation closure rather than a fabricated
successful terminal edge. This incomplete attempt establishes neither broad
generality nor performance; copy-elision promotion remains held.

Before freezing, the unexecuted root confirmation preparation is revised to
require this cancellation review. Its original prepared versions are retained.
The complete238 grid, binaries, verification, references and balanced order
remain unchanged. Dyadic controls retain the subsequent root confirmation
gate; their unrelated broad-completion gate likewise becomes cancellation
closure. No test or timing job overlaps the stopped broad campaign. Evidence:
`build/componentwise-inner-l2-copy-elision-root-broad-priority-cancellation-20261017/review.json`,
`build/componentwise-verdict-witness-root-confirmation-preparation-20261017/source-preparation.json`.


### Fresh root24 full peer comparison launched, 2026-10-17

Qualified certificate candidate `24dc56f` now has a fresh complete four-arm
comparison running against production root `2754413`, CKTSO and SubtreeLU.
Its prerequisite reviews cover 86 runtime tests, API28 and two separate
complete238 root comparisons. The rejected dyadic child is not included.
The same matched plain root binaries and committed normalized/cacheline
adapter are used; production componentwise verification remains enabled.

The fixed grid retains all seven original circuits, both CPU domains,
threads one/eight, four balanced passes and fourteen fresh KLU references:
462 observations and 924 child ledger edges. Each arm occupies each position
once per cell. Original controls, strict waveform/step gates and 1800-second
timeouts remain. The restored 719-byte packaged CKTSO key and both original
peer libraries are pinned. No failure is omitted or retried for performance.

Actual exec session56426/controller1863947 owns this attempt. No compiler,
test or profiler overlaps its timing. Progress inspection is limited to
counts and process liveness. Result assessment requires the actual terminal,
controller/child absence and a separate independent review. Root24 broader
matrix qualification is still pending; neither this diagnostic nor previous
peer ranks establishes production readiness or all-cell leadership. Evidence:
`build/componentwise-verdict-witness-root-peer-spice-preparation-20261017/preparation-receipt.json`,
`build/componentwise-verdict-witness-root-peer-spice-review-preparation-20261017/review.py`.


### Complete current-root peer gap, 2026-10-17

The full462 root24 four-arm campaign has actual exit zero, all924 recorded
child edges closed and every recorded PID absent. Its prepared independent
review subsequently closes with actual exit zero and rechecks the complete
order, pins, raw logs, waveforms and ratios. All KLS observations qualify.
KLS wins22 of56 peer cells, loses32 and has2 unqualified peer cells. Three
SubtreeLU MOS13/eight-thread observations time out: one on llc96 and two on
llc32. All remain retained; both corresponding aggregates are null.

Against CKTSO, root24 suite ratios are1.197380/0.662748 on llc96 at one/eight
threads and1.236668/0.675070 on llc32. Against SubtreeLU, one-thread ratios
are0.723907/0.706286 on llc96/llc32; eight-thread complete aggregates are null.
Every one-thread CKTSO cell remains slower. Parallel RAM, memory, DAC and BJT
remain gaps; SubtreeLU also wins parallel NEURON. Neither isolated root gains
nor favorable aggregates establish the required56 independently confirmed wins.

The implementation priority is current-root coarse diagnostics, observed-cycle
funding for deferred engines and transactional row preparation before removing
size/work admission windows. Numeric/solve execution is a necessary target:
existing root24 diagnostics put about75--81% of parallel RAM/memory/DAC API time
outside solve and only8--12% in certification. These instrumented proportions
are diagnostic evidence, not removable wall-time guarantees. Default verification
remains enabled; root24 broader qualification and promotion remain pending.
Evidence:
`build/componentwise-verdict-witness-root-peer-spice-review-preparation-20261017/review.json`,
`build/componentwise-verdict-witness-root-peer-completed-decision-20261017/decision.json`.


### Current-root coarse execution diagnostics implemented, 2026-10-17

Commits `e53e8a7` and `274ccf8` add a standalone source preparer and diagnostic
helper, plus controls for nested scopes, early returns, three thread lifetimes,
allocation requests and all64 SSE sticky combinations under four rounding
modes. Native, generic and optimized ASan/UBSan controls pass. Integer clocks
preserve errno and tested x87/SSE state. These tools are not included in an
ordinary production build and change no public API or solver layout.

The first complete28 diagnostic passes its API and independent waveform checks.
It retains32 function/job scopes. Review then identifies a necessary distinction:
worker zero can execute on the public caller thread. That completed attempt and
its exact historical tool bytes remain retained. The extended helper tracks
same-thread public-API ancestry separately from background workers, adds API
root scopes and checks an exact partition: summed API-thread exclusive time
equals summed outermost API-scope time. Distinct thread counts describe thread
lifetimes entering scopes, not simultaneous occupancy or arithmetic utilization.

A fresh full28 extended diagnostic passes independent raw-log, waveform, step,
input/binary pin and process-ledger review. All1036 coarse records and all28
exact API-thread partitions reconstruct. Both full diagnostics close with
actual exit zero, all56 child edges closed and reviewers absent. The extended
source recovers all33 original src/include Git blobs after removing diagnostic
insertions; only kls.c.o changes in the seven-member root24 archive. Other core
objects, the committed adapter and dependencies remain pinned. Allocation
requests cover wrapped link boundaries, not retained or process memory.

In the extended instrumented run, mapped refactor occupies about57--68% of
parallel RAM/memory/DAC API-thread time. Their retained crew enters only a small
number of initial numeric trials while steady execution remains mapped. Serial
BJT spends about84--88% in SNB refactor; parallel BJT spends about20% in native
transpose solve. Operator snapshot preparation also requests roughly one value
allocation per update and occupies about7--13% in several remaining gaps.
These proportions identify execution targets, not removable-cost guarantees:
instrumentation can alter adaptive choices and independent API calls may overlap.
Worker totals must not be treated as caller wall time.

Next implementation: private observed-work accounting and a shared bounded
exploration ledger, preserving immutable caller hints and zero-hint initial
analysis policy. Existing declined-row complete-cycle re-audits already exist;
do not duplicate them or force a sampled loser. Budgeted row construction must
preserve row/lean/compact and off-diagonal map ownership, join real-owner pool
jobs, charge map/cache work before estimation, and establish native/packed value
authority before seeding. Broader admission remains held until rollback, OOM and
cancellation controls pass. No algorithmic speedup, generality, promotion or
fastest claim follows from these diagnostic implementations.
Evidence:
`build/root24-coarse-profile-api-lineage-controls-20261017/review.json`,
`build/root24-coarse-execution-api-lineage-api-review-preparation-20261017/review.json`,
`build/root24-coarse-execution-completed-priority-decision-20261017/decision.json`.

### Observed-work exploration funding implemented, 2026-10-17

The isolated branch `work/observed-exploration-budget-root-20261017` adds private
per-instance work accounting and a shared reservation ledger. Commits `46fa8cc`
and `a6df368` implement the ledger and its first SNB consumer; `5dbae38` corrects
the fixture constructor and exercises real zero-hint observation. The initial
integration attempt compiled the solver archive but failed while compiling the
new fixture. That attempt is retained and claims no runtime qualification.

The corrected exact candidate passes45 native and45 optimized ASan/UBSan suite
checks, retaining all43 parent tests in each context. The initial three focused
checks also pass in each context. Independent review checks source/dependency
archives, active assertions, unchanged public header/dedicated certificate
sources and every child edge. Five standalone budget contexts additionally
cover integer timing conversion, rounding/sticky flags, reservations, stale
leases, overflow and per-instance independence. The lifecycle fixture checks
successful external-call counting, internal-call exclusion, failures, both
orientations and formats, reset, immutable hints and actual zero-hint work.

Funding uses observed successful external work and optional immutable caller
hints. One active reservation fixes an experiment's allowance. Rejected work
and deadline overshoot remain charged; observation monitoring itself must be
funded. Ordinary bookkeeping reuses existing timings and adds no sampling
clock or allocation. SNB is the first consumer; other adaptive engines still
use their existing policy. This is not yet a global exploration-cost bound
across all engines, a memory cap or broader row admission.

All28 Xyce API settings pass independent waveform, step, certificate, binary
and source checks. The plain parent24/candidate5d binaries use the same
committed adapter and shared dependencies. Default current-operator verification
remains enabled. The first independent API reviewer completed but its exec
terminal handle was lost during context truncation; its package is preserved.
A fresh independent closure reviewer repeats the checks without repeating the
numerical campaign, reaches actual exit zero and is confirmed absent.

The complete238 paired application-cost comparison against parent24 and its
independent review finish with actual exit zero, all476 child edges closed,
and controllers/reviewer absent. All238 observations qualify. The candidate is
faster in23 of28 cells and slower in five, retaining every raw pair. Full-case
candidate/parent geometric means are0.990576 (llc96 t1),0.971570 (llc96 t8),
0.995024 (llc32 t1) and0.976490 (llc32 t8). All14 parallel medians improve.
Serial regressions remain: NEURON +0.673%/+0.134% and BJT +0.209%/+0.346%
(llc96/llc32), plus MOS llc32 +1.601%. No slower cell is removed or given a
special selector. The suite gains establish no causal decomposition or all-cell
advantage. This foundation remains unpromoted; generality, current-peer and
separate full confirmation gates remain outstanding. Transactional row ownership,
resource accounting and paired verified trials remain prerequisites to broader
row admission. Completed evidence:
`build/observed-exploration-budget-spice-pair-review-preparation-20261017/review.json`,
`build/observed-exploration-budget-completed-cost-decision-20261017/decision.json`.
 A separate parent24 snapshot-storage-reuse candidate `ed555c5`
was source-prepared, then shelved before compilation after finding the retained
Oct12/13 experiment of the same hypothesis. Even its improved copy-contract
implementation had no aggregate benefit in the complete238 comparison. The new
source branch is retained but has no new runtime qualification, and its prepared
controllers refuse execution. Do not repeat this rejected experiment without
new causal evidence. See `build/original-snapshot-reuse-source-only-shelving-20261017/decision.json`.
Evidence:
`build/observed-exploration-budget-controls-preparation-20261017/review.json`,
`build/observed-exploration-budget-fixture-corrected-full-runtime-review-preparation-20261017/review.json`,
`build/observed-exploration-budget-api-closure-review-preparation-20261017/review.json`.

### Detached row/lean ownership qualified, 2026-10-17

The isolated `work/row-metadata-owner-root-20261017` branch at `ca36e74`
implements a typed, move-only metadata owner. It exchanges the complete303-field
footprint of the row/solve-partition, input-map and LU-cache destructors while
leaving the solver and pool as their real owners. The footprint includes lean
scheduling aliases, both transpose plans and the off-diagonal maps replaced
outside the old row destructor. Disposal preserves the incumbent while freeing
candidate state and never dereferences borrowed factor values. This alone is
not a constructor/cancellation or factor-storage/value-epoch transaction.

The exact candidate passes46 native and46 optimized ASan/UBSan suite checks,
plus four initial focused checks in each context. Independent review confirms
all45 parent tests remain, every test name/log, source and dependency archives,
active new-fixture assertions, all303 direct destructor fields, unchanged public
header/dedicated certificate sources and all20 child edges. The native/SAN
fixture stages real lean/row constructors through intercepted allocation
failures across formats, orientations, threads and scaling, including coupled
two-block BTF maps. It checks the incumbent field by field and solves current A
after cleanup. Actual execution and review terminate zero with their recorded
PIDs absent. The kls.c change is only an inactive private include; there is no
automatic execution-path caller or wider admission. Evidence:
`build/row-metadata-owner-full-runtime-review-preparation-20261017/review.json`.

A separate construction-only child on `work/detached-row-constructor-root-20261017`
at `c84664f` uses generic serial metadata preparation to retain the incumbent
crew, while preserving the legacy parallel-capable builder body at its existing
callers. Its new fixture covers zero/repeated hints, native factor and crew
preservation, BTF off-diagonals and allocation faults. It also includes the
exploration header in tuning-profile fingerprints. All47 native and47 optimized
ASan/UBSan suite checks, plus five initial focused checks in each context, now
pass. Independent review confirms exact source/dependencies, every parent test,
active assertions, all20 child edges, the inverse-exact legacy builder body at
allow_parallel=1 and absence of any automatic constructor caller. Actual runtime
and reviewer execution terminate zero with recorded PIDs absent. Evidence:
`build/detached-row-constructor-full-runtime-review-preparation-20261017/review.json`.
This qualifies construction-only controls, not bounded resource/cancellation,
current-value/storage epochs, complete engine trials or a performance/admission
change. Its helper still has no automatic caller. A source-only inventory identifies21 directly
allocating helpers for the next explicit resource context; that inventory is
not a semantic proof of a complete cap. Bounds, cooperative cancellation,
epochs, verified full-cycle trials and broader admission remain outstanding.

### Caller-owned row preparation budgets implemented, 2026-10-17

The isolated `work/bounded-row-preparation-root-20261017` branch adds a
lease-bound preparation context (`d7e844c`) and real rollback controls
(`6fe3b74`). The serial path routes98 explicit allocation sites through checked
count/size accounting. Five helpers receive the context explicitly; solver-aware
helpers use a scoped cold-tail pointer and the row-level job carries its context.
No process-global/TLS state or new unlimited-path sampling clock is introduced.
Grants remain fixed, allocation attempts/temporary storage receive no refunds,
and measured cleanup and deadline overshoot remain debited. Unreadable/backward
clocks disable further optional funding; stale frames cannot mutate a new lease.
The scoped pointer remains installed through cleanup and blocks reentrant frames.

The real constructor fixture sweeps every recorded cumulative-payload boundary,
every clock cutoff including final accounting, and each intercepted allocator
failure across64 format/orientation/thread/scale/hint/BTF configurations. It
checks the unlimited pattern, all incumbent ownership fields, native LU/Offx/
permutations, crew/scratch identity and current caller solves after rejection.
The limited-constructor, unlimited-constructor and pure-budget fixtures all pass
native and optimized ASan/UBSan: six selected test executions. Independent review
checks committed source/dependency archives, assertions, all12 child edges and
unchanged ordinary public bodies/dedicated numeric/certificate sources. An
independent argument audit preserves the original allocation operands at all98
rewritten sites. The initial preparation stopped before archive/compilation on a
cwd-relative historical hash; it is retained. A fresh package resolves the two
known historical path prefixes explicitly without changing the numerical source.

A test-only follow-up `02373f2` passes all128 standard exception-mask/rounding
combinations in five compiler contexts (native, generic, O0, x87 and optimized
ASan/UBSan). Its resource headers match the focused constructor candidate. All
executions/reviewers exit zero and their recorded PIDs are absent. Evidence:
`build/bounded-row-preparation-focused-runtime-review-preparation-20261017/review.json`,
`build/bounded-row-preparation-allocation-arguments-review-20261017/review.json`,
`build/row-preparation-budget-fp-masks-controls-review-preparation-20261017/review.json`.

This is focused qualification, not a full49-test suite or Xyce/cost result.
Cancellation checks currently cover allocations and batched root/input-cleanup/
row-level scans; other helper phases/library calls can overshoot. Requested
payload excludes allocator bookkeeping, sorting scratch, stack use and RSS.
Complete cancellation/scratch coverage, storage/value epoch retirement and
verified complete-cycle trials remain before automatic admission. No dimension
gates, default verification, production promotion or peer-leadership claim change.


The isolated `work/bounded-row-preparation-scratch-root-20261017` follow-up
`5a6c3fa` bounds sorting scratch in the two reached preparation sort paths using
an allocation-free heap sort with cooperative checks. Unlimited callers retain
libc qsort. Cancellation also cleans up sorted-L cache temporaries and cached
work-sort entries. No public admission or default-verification policy changes.

All five standalone compiler/FP controls passed, including exact qsort comparison,
all observed clock cutoffs, record preservation and no sort allocations. Three
focused tests each passed native and optimized ASan/UBSan. The real-constructor
fixture now covers65 configurations: the added BTF case asserts sorted-L cache
construction with a nontrivial column using test-only work estimates restored
before public solves. Independent closure/source reviews passed:
`build/row-preparation-sort-controls-review-preparation-20261017/review.json` and
`build/bounded-row-preparation-sort-focused-runtime-review-preparation-20261017/review.json`.
All controller/child processes exited zero and are absent. This is limited
construction qualification, not a full49-test suite, API or performance result.
Production numerical code stays at2754413. Remaining integration requires
cancellation in other helper scans, distinct storage/value validity, verified
complete refactor-plus-solve trials, then the full peer grid and confirmation.


The isolated `work/row-preparation-scan-cancellation-root-20261017` follow-up
`8f353be` adds cooperative cancellation and local cleanup in seven map, L/U32-bit
cache, factor-pointer and transpose-plan helpers. Ordinary callers take no new
sampling clock. The limited constructor retains complete-candidate rejection;
this does not introduce an automatic caller or change verification policy.

All three focused tests passed native and optimized ASan/UBSan (six executions).
The real-constructor fixture now covers67 configurations, including two513-column
cases crossing batched scan boundaries. Failure sweeps retain the unlimited
pattern, incumbent metadata, native factor hash, crew and current public solves.
A fresh independent source/runtime review passes:
`build/bounded-row-preparation-scans-focused-runtime-review-preparation-20261017/review.json`.
All12 child edges and controllers closed zero and are absent. Separate source
reviews preserve all98 allocation operands and show inverse equivalence of the
seven changed helpers with byte-identical source outside those helpers:
`build/bounded-row-preparation-scans-source-review-20261017/review.json` and
`build/bounded-row-preparation-scans-inverse-review-20261017/review.json`.

This is focused preparation qualification, not complete helper-branch coverage,
a full49-test suite, API validation or a speed result. Remaining work includes
dependency/group, panel and partition scan review, storage/value validity,
verified complete refactor-plus-solve trials and general automatic admission,
then the full56 peer cells and a separate complete confirmation. Production
numerical code remains2754413; default verification stays enabled.


The isolated `work/row-preparation-group-cancellation-root-20261017` changes
`0680500` (implementation) and `d4a60ea` (direct controls) extend cancellation to
twelve group DAG, dependency, compact-panel, producer-run, segment-target,
classification and work/order helper bodies. Serial jobs and three leaf helpers
carry the caller-owned context explicitly. Unlimited callers sample no new
clock; production selection and default verification are unchanged.

Five direct metadata families now reach branches absent from sparse fixtures,
including a257-entry producer run, actual compact storage and target cleanup.
Successful outputs match unlimited helpers; input hashes stay unchanged. Each
native and optimized ASan/UBSan context passes1661 observed clock cutoffs,
18 payload boundaries and18 allocator cutoffs across these families, alongside
all67 real-constructor configurations. Three selected tests pass per context
(six executions). The independent review checks identical direct records,
source/four dependencies, all12 closed child edges, actual zero exits/absent
PIDs, active assertions, executable/log hashes and unchanged public/numeric/
certificate sources:
`build/bounded-row-preparation-groups-focused-runtime-review-preparation-20261017/review.json`.
Separate source audits preserve all98 allocation operands and inverse equivalence
of the twelve changed helper bodies, with only explicit context wiring outside:
`build/bounded-row-preparation-groups-allocation-review-20261017/review.json` and
`build/bounded-row-preparation-groups-inverse-review-20261017/review.json`.

Remaining work includes cancellation across finish-stage local/solver ownership
transitions and solve-partition helpers, storage/value validity, verified complete
refactor-plus-solve trials and general automatic admission. Full-suite/API/cost
qualification, all56 peer cells and a separate complete confirmation remain.
This helper qualification is not a performance result; production numerical code
remains2754413.


The isolated `work/row-preparation-finish-cancellation-root-20261017` implementation
`37f0b24` checks all28 for-loops and three while-loops in final row metadata
assembly. Cancellation cleanup distinguishes local arrays from arrays already
owned by the candidate. The bounded row-shape comparator checks in256-index
chunks; unlimited callers retain the prior memcmp path and allocation behavior.

The qualified2348b9f fixture covers69 constructor configurations, including two
actual33-row dense panel builds with scoped test-only thresholds restored before
public solves. Each dense case builds compact storage normally and accepts two
valid missing-panel fallbacks over86 intercepted allocation cuts. Every accepted
OOM candidate's base pattern matches the unlimited reference, while incumbent
metadata, native factor hash, crew and current public solves are preserved.
Six direct helper families include equal/late-mismatching257-entry shapes and
exercise1667 clock,18 payload and18 allocator cutoffs per context. All three
focused tests pass native and optimized ASan/UBSan (six executions), with equal
direct and fallback records:
`build/bounded-row-preparation-finish-optional-fallback-focused-runtime-review-preparation-20261017/review.json`.
Independent review also checks exact source/four dependencies, all12 closed
child edges, actual zero exits/absent PIDs, assertions, executable/log hashes and
unchanged public/numeric/certificate sources. Source audits preserve all98
allocation operands, finish operations and the unlimited comparison path:
`build/bounded-row-preparation-finish-allocation-review-20261017/review.json` and
`build/bounded-row-preparation-finish-inverse-review-20261017/review.json`.

The first executed fixture stopped on an overstrong assertion requiring optional
compact storage even during intercepted allocation failure. Its failed package
is retained; the correction keeps unfaulted coverage assertions and verifies the
accepted fallback patterns, with byte-identical numerical source. Details:
`build/bounded-row-preparation-finish-shapes-focused-runtime-preparation-20261017/failure-analysis.json`.
The corrected execution used fresh source/build/runtime packages. An earlier
archive-only package was never compiled.

Dense-tail/slice/segment/sparse-level partition cancellation, final preparation
closure review, storage/value validity, verified complete refactor-plus-solve
trials and general automatic admission remain unfinished. Full-suite/API/cost
gates, all56 peer cells and a separate complete confirmation are still required.
Production numerical code stays2754413 and default verification remains enabled;
this is focused construction qualification, with no speed claim.


The isolated `work/row-preparation-partition-cancellation-root-20261017` changes
`cc74ca8` (implementation) and `620e666` (controls) extend cancellation to the
remaining dense-tail, slice, segment and sparse-level preparation scans and stop
partition/transpose drivers between failed phases. Local arrays are released
before publication. Source inversion reproduces all seven changed helper bodies
apart from added checks/cleanup and pointer declaration moves, with byte-identical
source outside them; all98 allocation operands remain preserved.

Six new513-row direct controls exercise lower/upper levels and segment splits,
dense-tail selection and slice bounds, using valid test tuning and known expected
geometry plus unlimited comparison. All12 direct families pass1740 clock,29
payload and29 allocator cutoffs in each native and optimized ASan/UBSan context,
with identical records. The69 constructor configurations and two valid panel
fallbacks per dense case also pass. Three selected tests pass each context
(six executions). Independent review checks source/four dependencies, assertions,
executable/log hashes, all12 closed child edges, actual zero exits/absent PIDs
and unchanged public/numeric/certificate sources:
`build/bounded-row-preparation-partition-focused-runtime-review-preparation-20261017/review.json`.
Source evidence:
`build/bounded-row-preparation-partition-allocation-review-20261017/review.json` and
`build/bounded-row-preparation-partition-inverse-review-20261017/review.json`.

The known serial source slice has also been reviewed:57 helper bodies contain144
cooperatively checked loops and three halving binary-search leaves. Explicit
mode0/0 guards exclude parallel/lean builders, the legacy aligned allocator is
replaced in this slice, and bounded sorts avoid libc sorting scratch. Fixed
macros/atomic initialization and bulk memory/free operations remain explicit:
`build/bounded-row-preparation-serial-path-closure-review-20261017/review.json`.
This is a scoped source review, not a formal preprocessed callgraph proof,
exhaustive branch coverage, hard deadline or RSS bound; opaque library/callback
latency and cleanup overshoot remain charged. The retained first source-review
attempt stopped on a missing prior terminal receipt and was replaced after the
actual zero source-audit receipts were recorded; no numerical result was retried.

Storage/value validity, verified complete-cycle trials and general automatic
admission remain unfinished. Later worker/thread-bound setup is outside this
serial preparation slice. Full-suite/API/cost qualification, all56 peer cells and
a separate complete confirmation are still required. Production numerical code
remains2754413 and default verification remains enabled; no speedup is claimed.


The isolated `work/row-metadata-epochs-root-20261017` implementation3cf77a2 binds
detached metadata to separate storage/layout and value-attempt tokens. It uses
caller-thread hooks at live numeric release/adoption, sorting, block repair and
predicted repack/pivot fill; pure metadata cleanup does not advance them. Valid
numeric attempts conservatively invalidate copied values, including failed/no-op
attempts. Tokens persist through reanalysis and saturate fail-closed for optional
exchange, while ordinary numeric operations continue. Stale disposal remains
metadata-only and preserves the incumbent. No active ROW repair vehicle is
removed by the new hooks.

The4477456 hook/fixture source passes all49 runtime tests each native and optimized
ASan/UBSan (98 executions), plus four initial focused tests per context. Its new
owner controls exercise20 real lifecycle cases and128 FP mask/rounding states.
Native runs observe both value updates retaining storage and replacement reusing
the numeric object's address, with stale exchange rejected. Actual same-numeric
block LU replacement, repair OOM, release, reanalysis, metadata-only cleanup,
invalid calls and overflow are covered. The69 constructor configurations and12
direct preparation families remain qualified.

The finalc8cb8c3 refinement checks versions before comparing stored integer numeric
identities and uses integer address snapshots in tests. This exact two-file delta
changes only the private owner guard/test snapshots; public numeric hook files
are byte-identical to4477456. It separately passes four focused owner/constructor/
preparation/budget tests each native and optimized ASan/UBSan (eight executions).
The final refinement did not repeat the full49 suite. Independent reviews verify
source/four dependencies, assertions, executable/log hashes, all49 retained CTest
names, all20 full-run and12 final-focused child edges closed, actual zero exits
and absent PIDs:
`build/row-metadata-epochs-full-runtime-review-preparation-20261017/review.json` and
`build/row-metadata-epochs-integer-identity-focused-runtime-review-preparation-20261017/review.json`.
Source inversion preserves prior numerical operations, all303 payload fields,
swap/disposal bodies, public headers and dedicated certificate sources:
`build/row-metadata-epochs-source-review-20261017/review.json`.

Strict epoch rejection is not stale-value refresh or a trial coordinator. Reusing
geometry across updates, establishing native authority, explicit versioning for
direct private-kernel trials, incumbent-last recomputation and current-A/RHS
certification at unchanged tolerance remain required. Cleanup/recovery and later
worker setup must be charged. No automatic admission, API cost, Xyce, speed or
peer qualification follows. Production numerical code stays2754413, default
verification stays enabled, and all56 peer cells plus a separate complete
confirmation remain required.


Bounded detached value refresh is now implemented and qualified on experimental
branch `work/row-metadata-value-refresh-root-20261017`, code checkpoint13968f8
(worktree `build/observed-exploration-budget-root-worktree-20261017`). Native
authority and a current storage binding are explicit preconditions. The helper
allocates nothing and preserves live metadata/factors and crew/scratch state.
It copies row and sorted-L values, invalidates packed/panel authority and solve
readiness, accounts elapsed work, and publishes a current value binding only
after final budget/version checks. Interrupted work remains unbound and retryable.

The frozen076f28e source passes five focused checks per native and optimized
ASan/UBSan context. Its actual37 refresh configurations are recorded accurately
by review despite an incorrect69 summary in that test version. Final13968f8
adds hinted/unhinted coverage, portable-C fixes, storage boundaries, guards and
FP controls. Its refresh fixture passes separately in both contexts:69
configurations,128 FP mask/rounding combinations and646 cutoffs each for expiry,
clock failure and value-version changes, plus packed-address rejection and
stale retirement. Review confirms exact source/four dependencies, assertions,
all49 parent CTest names retained plus the new test, unchanged prior numerical
operations by source inversion, and all12 child edges closed per run with
actual zero exits and absent PIDs. This does not repeat a full50-test suite.
Evidence:
`build/row-metadata-value-refresh-fixture-focused-runtime-review-preparation-20261017/review.json`
and
`build/row-metadata-value-refresh-final-pinned-focused-runtime-review-preparation-20261017/review.json`.
A mistaken final archive pin was rejected before any compiler/test child; the
failed package is preserved and the corrected run uses a fresh package.

The next blocker is complete paired-trial coordination: explicit attempt
versioning for real private kernels, current-operator publication, saved caller
RHS ownership, funded worker setup and refresh, actual incumbent-last
recomputation, unchanged-tolerance current-A/RHS certification and recovery/debt
accounting. Metadata refresh alone neither restores numerical trial writes nor
proves row execution or a speedup. No automatic caller or broader selector is
added. Production numerical source remains2754413 and default verification stays
enabled. Complete candidate runtime/API/cost gates, all56 peer cells and a
separate complete confirmation remain required.


The private row paired-cycle path is now implemented and qualified on
`work/row-paired-cycles-root-20261017` (same experimental worktree; final code
checkpoint259c4d6). It borrows the enclosing API frame's saved RHS and separate
output, performs real incumbent/row/incumbent numerics on the current caller A,
and uses the unchanged caller certificate/recovery postlude. A private repeat
counter defeats generic/tiny unchanged shortcuts; forced input-cache lifecycle
now invalidates before writes and republishes after success, with A/B/A controls.
Empty L/U stream publication is fixed generally, including singleton and both
triangular directions. Metadata moves never substitute for the mandatory final
full recomputation. Expired/failed/stale restoration retires old metadata, and
cleanup plus recovery still run and are charged. Candidate answers stay private.

Shared public-path sourceee53d8e passes all51 CTests each native and optimized
ASan/UBSan (102 full-suite executions), plus one paired preflight per context.
Independent source inversion preserves numeric/caller postlude operations apart
from the explicit forced-input guards and empty-stream publication fix. The
final259c4d6 delta changes only the private pair module and fixture; public-path
source and CMake remain identical. It separately passes the paired fixture each
context, not the full51 again. Each final run covers214 paired cases:128 ordinary
format/orientation/thread/scale/hint/BTF/RHS configurations,32 triangular, two
singleton, two dense, five independent-component cases (four dispatching eight
effective row threads and a serial control),16 fault cases,28 observed expiry
edges and the reference cycle for that sweep. It additionally covers2068/1048
packed restoration expiry cutoffs on513-column/258-row dense fixtures and all
four authority combinations, without allocations or live-factor/metadata writes.
The collector retains151 stable-path pairs; post-certification storage/value
changes invalidate samples. These are injected test clocks, not speed evidence.

Reviews confirm exact sources/four dependencies, assertions, all50 parent CTest
names retained, actual CTest commands (including eight script-based tests),
executable/script/log hashes, all20 full-run and12 final-focused child edges
closed, actual zero exits and absent PIDs:
`build/row-paired-cycles-packed-empty-columns-full-runtime-command-review-preparation-20261017/review.json`,
`build/row-paired-cycles-final-workers-focused-runtime-review-preparation-20261017/review.json`,
and `build/row-paired-cycles-source-review-preparation-20261017/review.json`.
Failed compilation/preflight packages and mistaken review artifact assumptions
are preserved; corrected evidence comes from fresh packages. The flushed
preflight diagnostics identified singleton publication, not an expiry failure.

The next work is automatic integration: exclusive lease ownership, funded
persistent preparation, worker/thread-bound payload accounting, repeated
samples and payback using actual observed updates, and safe adoption/retirement.
Full elapsed setup/execution/recovery is charged by the prototype now, but later
worker allocations are not all inside the serial payload-byte cap; no hard
realtime/RSS bound is claimed. Actual SNB/EGraph/packed-incumbent campaign routes
need live coverage. No automatic selector, API cost, Xyce or peer qualification
is claimed. Production numerical source stays2754413, verification stays enabled,
and all56 peer cells plus separate complete confirmation remain required.


### Funded automatic row search, complete correctness qualification

The automatic caller is implemented on `work/row-funded-search-root-20261017`
in `build/observed-exploration-budget-root-worktree-20261017`. Final code is
4193d10; implementation/qualification notes are committed separately. This is
an experimental branch. Production numerical source remains2754413.

The outer combined API owns its saved RHS and output. An exclusive fixed lease
funds persistent serial preparation and four verified current-A row/incumbent
comparisons. Median savings must exceed the existing2%/1 microsecond margin
and repay twice the measured search debt using actual observed updates and
caller-declared reuse. Incumbent-engine changes discard waiting samples without
forgiving debt. Storage changes retire the decision; analysis disposes it.
Insufficient funding takes no search clocks or allocations. Adopted row failures
perform verified incumbent recovery even without optional credit. No circuit
identity, CPU-domain label or new dimension-range selector is added.

Forced comparison now mirrors the ordinary input-cache transition, preserving
the settled incumbent's normal dispatch. Qualification exposed and fixed stale
row copies after native KLU updates and analogous SNB/EGraph handoffs. EGraph
also retains existing packed authority independently of new-admission policy;
it cannot overwrite a fresh native update with an old valid packed copy.
Compact-panel copies are invalidated when their row values change. The1e-8
componentwise tolerance and1e-10 known-solution assertions are unchanged.

Final source4193d10 passes all52 CTests each native and optimized ASan/UBSan
(104 full-suite executions), plus one initial search preflight each. The fixture
covers128 adoption configurations and eight additional refusal/fault/cache/
incumbent-change modes, aliases/strides, real public adopted calls, mandatory
allocation-failure recovery, factor-version retirement and analyze reset. Six
actual SNB/packed-EGraph configurations execute three changed updates each.
Thirty-two independent components meet the unchanged eight-worker EGraph gate.
Injected clocks establish no performance result.

Independent review confirms exact source/four dependencies, assertions, real
CTest commands and executable/script hashes, all51 parent names retained,
all20 child edges closed, actual zero exits and absent PIDs:
`build/row-funded-search-verified-guard-full-runtime-review-preparation-20261017/review.json`.
The intermediate source7400eac also passed the complete native/sanitizer suite
and four initial pair/search checks:
`build/row-funded-search-packed-authority-full-runtime-review-preparation-20261017/review.json`.
Failed runs remain preserved; corrected fixtures use fresh packages. Failures
include an incorrect public-factor replacement assumption, missing test snapshot
frame, insufficient test graph width, and an unjustified first-trial sample
assumption, alongside the actual stale-row/packed-authority defects.

Measured preparation/setup/execution/recovery time is charged. State, serial
preparation and trial output use checked payload requests; worker/kernel
allocations retain their existing library limits. A complete search-memory cap
and hard realtime/RSS guarantee are not established. See the experimental
branch's docs/row_metadata_ownership.md for ownership and qualification details.

The exact final native archive is now running the complete28 circuit/API
diagnostic against parent24dc56f, with the same committed adapter/dependencies.
The complete238-observation cost comparison is prepared but unexecuted. Complete
API/cost qualification, all56 live peer cells and separate complete confirmation
remain required. Verification remains enabled. No fastest or production
promotion claim follows from the correctness controls.


The exact4193d10 native archive now passes all28 circuit/API settings, with
288808 certificate calls and unchanged tolerances. Independent waveform, step,
certificate, matched-link and process-closure review is complete:
`build/row-funded-search-api-closure-review-preparation-20261017/review.json`.
All56 API and six build child edges are closed; actual terminal exits are zero.

The reviewed profiles contain zero exposed accepted ROW dispatches across all28
settings. Automatic row selection has therefore demonstrated no workload
benefit here. The238-observation cost campaign is prepared but has not started;
it is held pending a full28 diagnostic of funding, construction, comparison and
adoption decisions. Diagnostic-only integer counter insertions have an exact
source inverse to4193d10 and use one end-of-process report; numerical bodies and
policy conditions remain unchanged. Adaptive timings may be affected, so those
diagnostics will not establish a speedup or authorize production promotion.


The full28 decision diagnostic now passes independent inverse-source, archive,
waveform and closure review:
`build/row-funded-search-decision-api-review-preparation-20261017/review.json`.
It records288540 entries,171836 funding refusals,85 constructor attempts and52
constructor time exhaustions. There are46 verified pairs and ten complete
sample sets rejected by the existing savings margin; no positive payback,
exchange or adopted cycle occurs. The counter does not distinguish a slower
median from failing the absolute/relative minimum. Instrumentation may affect
adaptive decisions, so these counts establish no uninstrumented speed result.
An initialized incumbent-result sentinel must not be labeled a public failure.

Preparation/comparison expense and incomplete or unprofitable samples now take
priority over wider admission or a cost campaign. General next work concerns
minimal/shared/resumable epoch-guarded metadata preparation, pricing mandatory
useful work once, and separating rejection by funding, stability and margin.
No tolerances, corpus coverage or production default change. The cost campaign
has not launched; production remains2754413. Completed decision:
`build/row-funded-search-decision-completed-priority-20261017/decision.json`.


### Native-authority serial row preparation

Code40f521c on `work/row-native-preparation-root-20261017` reads retained LU and
internal CSC directly, avoiding generic map/cache construction solely for the
row transpose. It retains the complete group/solve finisher, numerical kernels,
incumbent/crew ownership and original-constructor fallback. There is no new
circuit/domain/dimension/reuse selector. All new allocation and scan work uses
the existing checked preparation context.

All53 tests pass each native/optimizedASanUBSan (106 full-suite executions),
plus two initial search/native-preparation checks each. All52 parent names
remain. The71-configuration differential/fault fixture proves the same row
geometry, unchanged incumbent/native factor/crew before execution and three
actual changed-A two-RHS cycles, including eight effective workers. Requested
payload is18--30% lower than the old constructor. All19241 observed expiry,
5076 payload and5076 physical allocation cutoffs pass per context. This is a
payload-request result, not a peak-RSS or speed claim. Runtime review:
`build/row-native-preparation-full-runtime-review-preparation-20261017/review.json`.

The exact native archive also passes all28 original circuit/API settings,
with288823 certificate calls, unchanged tolerances, independently rechecked
waveforms/steps, same adapter/dependencies and matched root24 links. All56 API/
six build child edges, actual zero exits and absent PIDs are reviewed:
`build/row-native-preparation-api-closure-review-preparation-20261018/review.json`.
The profiles still expose zero accepted ROW dispatches. Whole-application cost
and peer superiority remain unproved; the cost campaign has not launched.
Completed decision:
`build/row-native-preparation-completed-priority-20261018/decision.json`.

The next code7ae90e3, on `work/row-phase-reservation-root-20261018` in the same
experimental worktree, separates one-time construction and recurring-comparison
reservation estimates. A completed candidate no longer requests its cold
construction high-watermark again. It uses existing budget readings and resets
phase estimates for a new factor-storage epoch, preserving all debt. Funding
fraction, deadlines, sample count, margins, payback and recovery remain as before.
Focused/full native checks pass; a fresh complete53-test sanitizer gate is in
progress. No performance conclusion follows. Production numerical source and
verification are unchanged; complete cost/all56 peers/confirmation remain.


Code7ae90e3 has completed qualification: all53 CTests pass in both native and
optimized ASan/UBSan builds (106 full-suite checks plus two focused checks).
All28 original circuit/domain/thread API settings also pass, with288892 current
operator certificate calls, unchanged waveform/step tolerances and zero exposed
accepted ROW dispatches. Independent closure reviews verify exact sources,
dependencies, actual controller exits and absent controller/child PIDs:
`build/row-phase-reservation-full-runtime-review-preparation-20261018/review.json`
and `build/row-phase-reservation-api-closure-review-preparation-20261018/review.json`.

A separate diagnostic derivative of that exact source now records the existing
four-sample medians and the original payback rejection branch, without adding
clocks or changing numerical/policy conditions. Its source edits reverse exactly
to7ae90e3; only the main member changes in the seven-member archive. All28 cells
pass independent waveform, step and certificate rechecks, with56 API and six
build child edges closed. Review and actual zero/absent-PID receipts are under
`build/row-phase-decision-state-api-review-preparation-20261018`.

All11 complete sample sets reject because the row path is slower than the
incumbent: sampled median ratios range1.9561--14.7307. None rejects solely on the
absolute1us or relative2% savings margin. This rules out lowering those margins
as a remedy for these completed sets; it does not characterize incomplete sets.
Across the diagnostic125 candidates start and44 pairs verify; remaining pair
reasons include funding, path changes and storage/restoration invalidation.
Reasons overlap, and the incumbent-failed counter includes an unexecuted
sentinel, not public solve failures. Diagnostic instrumentation can affect
adaptive choices; these are not uninstrumented application speed results.

The next implementation priority is numerical execution cost. Audit generic
serial worker initialization and scratch lifetime, and the row numeric work
against the optimized incumbent, before further reservation-only changes.
Preserve the complete kernel, general applicability, current-operator checks,
measured debt and mandatory recovery. The paper review supports amortizing
preparation and work-balanced subtree scheduling, but does not establish that
runtime verification can safely be disabled. Verification remains enabled.
No uninstrumented cost campaign has run for this source, and production
numerical source2754413 is unchanged. All56 peer cells and a separate complete
confirmation remain required for the fastest claim. Completed decision:
`build/row-phase-reservation-completed-priority-20261018/decision.json`.


Serial row execution now retains its five worker scratch capacities between
updates on `work/row-serial-scratch-root-20261018`: numerical codece2da22, final
NaN-poison fixture514a4bc. Each call still constructs fresh shared state and
telemetry and performs the original numerical operations. Only temporary buffer
pointers/capacities move into and back out of the local worker. free_numeric
releases and clears them; parallel and prefix-repair workers retain their
existing independent ownership. The new field sits at the cold tail: a source
layout proof verifies every original solver field declaration is an exact
prefix of the new structure, preserving existing offsets. No circuit, dimension,
CPU-domain or hint selector, new clock, funding grant or numerical shortcut is
introduced. Retained memory lasts until numeric cleanup; a process-RSS bound
is not claimed.

The corrected numerical code passes all54 CTests in native and optimized
ASan/UBSan builds (108 full-suite checks plus two initial focused checks).
The stronger test-only514a4bc source differs only in its fixture, with identical
numerical/dependency sources and flags. Both new executions pass24 configurations
covering CSR/CSC, normal/transpose, unscaled/scaled, BTF, dense and independent
component shapes. Each context checks192 certified RHS across96 changed-A
updates,120 failed buffer-growth attempts, NaN-poisoned floating scratch,
nonzero byte/index/object poison, zero warm numeric reallocations including
fresh realloc(NULL) requests, and free_numeric/re-analysis/destruction cleanup.
Independent reviews and actual zero/absent-PID receipts:
`build/row-serial-scratch-cold-tail-idle-full-runtime-review-preparation-20261018`
and `build/row-serial-scratch-nan-runtime-review-preparation-20261018`.

The exact corrected archive also passes all28 original circuit/API settings,
with288878 current-operator certificate calls and independent waveform/step
rechecks, matched links, unchanged adapter/dependencies/tolerances and all56 API
plus six build child edges closed. Review:
`build/row-serial-scratch-cold-tail-idle-api-closure-review-preparation-20261018/review.json`.
There are zero accepted ROW dispatches; an application speedup is not proved.
Some instrumented API qualification overlaps a test-only compiler, so those
instrumented times are not cost evidence. No uninstrumented cost/peer campaign
or production promotion follows. Completed decision:
`build/row-serial-scratch-completed-priority-20261018/decision.json`.

Performance investigation now returns to the default incumbent, with a separate
full28 SNB destination census on the best qualified root24dc56f. It describes
completed-build geometry and counts actual caller-side prelude/completion uses,
including unselected builds and internal cycles. The source reverses exactly to
root24, every other source file is identical, and it adds no timing reads or
worker hooks. The census asks whether wide edges have fully valid consecutive
producer-below destinations, allowing an experiment that computes the same
dense product and immediately subtracts it from the target panel. This could
remove T temporary store/load traffic while retaining producer accumulation and
edge order; it needs alias/exception/layout/crew proofs and complete correctness,
API and cost gates. The census is running in
`build/snb-destination-census-api-preparation-20261018`; results require actual
terminal/PID closure and independent review before interpretation. Geometry
counts indicate potential logical traffic, not measured CPU cost or speedup.
Production numerical source2754413 and default verification remain unchanged;
all56 peer cells and a separate complete confirmation remain required.


The full28 root24 SNB destination census has completed with independent source,
archive, process, waveform, step and certificate review (288650 certificate
calls). Sustained public SNB dispatch occurs in the two serial BJT cells:
1144 public updates and1145 completed SNB cycles per cell. Each retained layout
has1584926 wide-product output values per cycle, of which29456 (1.8585%) belong
to fully consecutive destination edges. Fully valid edges account for1489791
(93.9975%) of those values. A contiguity-only path would cover little of this
logical traffic. This is geometry evidence, not a measured CPU-cost bound.
Review and decision:
`build/snb-destination-census-api-review-preparation-20261018/review.json` and
`build/snb-destination-census-completed-priority-20261018/decision.json`.

A general sparse destination prototype is now code074066f on
`work/snb-fused-scatter-root-20261018` in
`build/certificate-verdict-root-worktree-20261017`. It computes the original
producer-ordered dense product and immediately updates arbitrary canonical
sparse destinations using masked gathers/subtractions/scatters. Holes still
participate in product arithmetic, preserving its exceptions; target arithmetic
and address accesses remain masked. A register-only compiler boundary preserves
the binary64 rounding formerly imposed by storing T, preventing contraction
of the final subtraction into a product. Distinct builder row-map destinations
and their separation from the solved B source justify immediate updates.
Ordinary and cooperative row-sliced calls use the helper with existing barriers.
The portable path retains the original product and scalar scatter. No clocks,
allocations, selectors, layouts, tolerances or verification policy change.

The corrected source passes all43 parent CTests in both native and optimized
ASan/UBSan builds (86 full checks plus four focused checks). The extended edge
fixture includes irregular destinations and holes; the guarded fixture executes
65280 cases per context with source/index/target page boundaries, four rounding
modes, five finite/nonfinite/range patterns and whole/split row slices, checking
stored bits and exception flags. Two additional private-kernel controls compile
for x86-64-v3 without AVX512F, pass the same guarded cases, and have no ZMM/opmask
registers in their product/scatter symbols. This is portable-kernel-path evidence,
not a complete portable-application qualification. The initial prototype package
had an erroneous fixture initialization bound; its native passes are explicitly
invalid. It was stopped and preserved, and only the fresh corrected074066f
package is used for qualification.

The corrected archive also passes all28 original circuit/API settings with
288647 certificate calls, independent waveform/step rechecks and matching
plain binary links. Reviews and actual zero/absent-PID receipts:
`build/snb-fused-scatter-corrected-full-runtime-review-preparation-20261018`,
`build/snb-fused-scatter-corrected-portable-review-preparation-20261018`, and
`build/snb-fused-scatter-corrected-api-closure-review-preparation-20261018`.
Completed correctness decision:
`build/snb-fused-scatter-corrected-completed-correctness-20261018/decision.json`.
The full238 paired cost package is prepared for four balanced passes against
root24 over every original case/domain/thread setting, plus14 KLU references.
No speedup, peer superiority or production promotion follows from these gates.
Production numerical source2754413 and default verification remain unchanged;
all56 peer cells and a separate complete confirmation remain mandatory.


The complete238 paired campaign for fused sparse destination code074066f is
finished and independently reviewed. All observations qualify, all476 child
edges close, actual controller/reviewer exits are zero and their PIDs are absent.
Every original28 setting has all four balanced pairs. The intended sustained
SNB target regresses: serial BJT candidate/baseline medians are1.0570807 on
llc96 and1.0649930 on llc32, with all four pairs slower in each domain. Across
the complete grid22 medians improve and six regress; those other medians do not
rescue a regression in the intended numerical work. Group geometric means are
1.0021325/0.9890547 (llc96 serial/eight threads) and1.0023697/1.0039417 (llc32).
Whole-application timing establishes the cost result, not an exclusive CPU
attribution to individual gather/scatter instructions. Code074066f is rejected
for promotion; its qualified artifacts and commits remain negative evidence.
Review and completed decision:
`build/snb-fused-scatter-corrected-spice-pair-review-preparation-20261018/review.json`
and `build/snb-fused-scatter-corrected-completed-cost-decision-20261018/decision.json`.

The next measurement returns to the best qualified root24dc56f. Count coverage
of eight-column product tiles and the corresponding logical L-value loads
versus existing four-column tiles, and characterize consecutive eight-row
subtiles inside irregular destination edges. This distinguishes useful load
reuse from the rejected masked destination update. Any next kernel keeps the
original scatter unless independent evidence supports replacing it, retains
producer order and exception behavior, and receives fresh full correctness,
API and uninstrumented cost gates. No circuit/domain/dimension rescue selector,
verification change, peer result, production promotion or fastest claim follows.


The completed caller-only product-tile census qualifies all28 original API
settings with independent waveform/step checks, exact instrumentation inverse,
and terminal controller/reviewer receipts. Serial BJT in each domain completes
1145 SNB cycles. Its weighted wide products number1814740270;1662347640
(91.6025%) fall in full eight-column tiles. The source-level L-load model falls
from10336108505 to5961044040 (42.3280%) for an eight/four/scalar column tiling.
Consecutive eight-row destination tiles cover56.0553% of product outputs, much
more than wholly contiguous edges, but no destination update is changed by the
next prototype. These are geometry and logical-load counts, not measured DRAM
traffic, exclusive CPU cost or demonstrated speedup. Census review and decision:
`build/snb-product-tiles-census-api-review-preparation-20261018/review.json` and
`build/snb-product-tiles-census-completed-priority-20261018/decision.json`.

The isolated branch `work/snb-eight-column-product-root-20261018` starts directly
at root24dc56f, excluding the rejected masked-scatter change. Commit1007710
reuses each L vector across eight independent columns, retaining producer order,
T scratch, original scatter and the existing four/scalar remainder. Commit593dcf0
adds a frozen root24 product oracle and protected-page differential tests for
rounding, FTZ/DAZ, padded B, gapped producers, row offsets and nonfinite values.
Native/sanitized qualification is pending. No source promotion, timing claim,
verification change or dimension/circuit/domain selector follows.


Candidate593dcf0 completes independent correctness qualification: all43 parent
test names pass native and optimized AddressSanitizer/UndefinedBehaviorSanitizer
(86 full tests), plus both initial SNB tests in each context. Each context
passes342720 protected frozen-parent comparisons across four rounding modes,
four FTZ/DAZ modes, row offsets, padded B, gapped producers and IEEE exceptional
inputs. Both report zero differing NaN payload words; finite/infinity/zero bits
and C/raw exception flags agree. The frozen oracle is the exact root24 function
with only its symbol renamed and one trailing separator blank line omitted.
Two private non-AVX512 controls pass with x86-64-v3, AVX512F absent and no ZMM or
opmask instructions in the checked product symbols. This is private kernel-path
evidence, not complete portable-application qualification.

All28 original API/circuit settings also qualify with288650 certificate calls,
independent waveform/step rechecks, the same seven API wrappers, and matched plain
links sharing the committed adapter and dependencies. All other src/include
bytes and solver/metadata layouts remain root24. Actual controller/reviewer
exits are zero and all child edges close with PIDs absent. Reviews:
`build/snb-eight-column-product-full-runtime-review-preparation-20261018`,
`build/snb-eight-column-product-portable-review-preparation-20261018`, and
`build/snb-eight-column-product-api-closure-review-preparation-20261018`.
Decision: `build/snb-eight-column-product-completed-correctness-20261018/decision.json`.
The candidate is qualified for full238 uninstrumented paired cost only. No
speedup or peer leadership is established, and production source2754413 and
default verification remain unchanged. All56 peer cells and a separate complete
confirmation remain necessary for the fastest goal.


The first complete238 paired cost of candidate593dcf0 is terminal and
independently reviewed: all238 observations qualify, all28 original settings
have four balanced pairs, all476 child edges close, and controller/reviewer
PIDs are absent after actual zero exits. Serial BJT improves in all four pairs
in each domain; candidate/parent medians are0.9811196712 (llc96) and0.9835459270
(llc32), reductions1.8880% and1.6454%. This is a modest whole-application gain,
not the42.328% modeled L-load reduction and not exclusive kernel CPU evidence.

The complete grid has16 faster cell medians and12 slower. Every MOS setting has
a slower median (1.00264 to1.01164), while serial RAM medians are1.00443/1.00589;
eight-thread NEURON pairs vary substantially. These regressions and variation
remain part of the decision. No aggregate or selected target can establish
unconditional promotion. Candidate593dcf0 remains isolated, with production
source2754413 and verification unchanged. A separate full238 confirmation is
prepared, followed by the complete56 peer grid and separate full confirmation
before any fastest claim. No circuit/domain/dimension rescue selector is added.
Review and completed decision:
`build/snb-eight-column-product-spice-pair-review-preparation-20261018/review.json`
and `build/snb-eight-column-product-completed-cost-decision-20261018/decision.json`.


The separate full238 confirmation of593dcf0 is terminal and independently
reviewed. All238 observations qualify, all476 child edges close, and the actual
controller/reviewer exits are zero with PIDs absent. Serial BJT medians are
0.9836230670 (llc96) and0.9865475913 (llc32): reductions1.6377% and1.3452%.
Together the two complete campaigns qualify476 application observations and
all16 serial target pairs improve. This supports a modest repeatable target
gain. The first campaign's MOS median regressions reverse in confirmation;
they are not established persistent regressions. Conversely llc32 RAM and
eight-thread mux medians are slower in both campaigns. Retain every setting
and pair; neither aggregate medians nor the target alone establishes a broad
speedup or unconditional production promotion. Eight-column products remain
a qualified building block, isolated from production. Review and decision:
`build/snb-eight-column-product-confirmation-review-preparation-20261018/review.json`
and `build/snb-eight-column-product-two-campaign-decision-20261018/decision.json`.

A1-2% serial SNB gain does not resolve the larger parallel mapped-refactor
blocker. The next diagnostic should measure exact caller-input bit changes
and their conservative closure through the retained column dependency graph.
The CKTSO paper's EGraph discussion (local text lines220-247) distinguishes
exact fixed-pivot dependencies from ETree bounds under changing pivots. This
motivates measuring retained-graph closure; neither paper establishes that
KLS may skip an unchanged column without further numeric proof.

The diagnostic must observe all original28 API/circuit settings, including
standalone refactor and combined refactor/solve entry points without counting
nested recovery twice. Compare each current input with the preceding successful
public numeric update. Validate the caller-to-internal value bijection, separate
BTF off-block value changes, seed diagonal-block columns with changed inputs,
and propagate through retained U-column predecessors in topological order.
For dynamic KLU row scaling, conservatively seed every diagonal-block input
using any row with a changed input; static external scaling does not justify
ignoring this coupling. Report both column counts and structural update-slot
weights, available/unavailable entry coverage, failures and overflow. These are
geometry/work models, not measured CPU savings or proof that reuse is safe.

Keep this census caller-only in an isolated root24-derived source package:
no worker hooks, solver/metadata layout changes, numerical writes, new clocks,
policy decisions, tolerance changes or verification changes. Prove its exact
source inverse and independently recheck all28 waveforms/steps/certificates.
Only useful complete-grid coverage justifies implementing a numerical reuse
path. Such a path must separately prove scaled/transformed coordinates, current
operator identity, pivot/repair invalidation, retained-factor consistency and
conservative rejection before full correctness and application cost gates.
The full56 peer cells and a separate complete confirmation remain mandatory
before the fastest goal; no new peer leadership or production promotion is
claimed by either the tile result or this diagnostic plan.


The caller-only input-closure census is implemented on isolated root24 branch
`work/mapped-input-closure-census-root-20261018`, commit7fbd73f. Native and
optimized ASan/UBSan private fixtures each pass8192 independent queue-reference
DAG cases, BTF off-block/scaling coupling, bad-map/DAG rejection, nested-call
suppression, failed-update invalidation, pointer reuse and unchanged FP flags.
Independent review proves the exact source inverse; all other src/include
bytes remain root24, with no solver-layout, worker, clock, numeric or policy
change. The fixture translation unit uses native support libraries; this is
not full sanitized-application qualification.

The full28 API/circuit campaign is terminal and independently reviewed. All
cells qualify with288647 certificate calls, exact seven-wrapper call-count
reconciliation, fresh waveform/step checks, only the main archive member
changed, all62 API/build edges closed, and actual zero controller/reviewer
exits with PIDs absent. Every MOS setting includes seven failed numeric
attempts followed by successful circuit recovery; they remain explicitly
recorded. No malformed map/DAG, tracking loss, overflow or OOM is reported.
Reviews: `build/mapped-input-closure-census-fixture-review-preparation-20261018`
and `build/mapped-input-closure-census-api-review-preparation-20261018`.

The modeled affected update-slot fraction in mapped mem_plus is99.1795%
(serial) and99.1893% (eight threads), DAC99.9459%, and RAM99.9049-99.9050%,
in both domains. Raw input changes already seed91-98% of their factor columns.
Do not implement selective skipping for these main mapped targets: little
potential remains even before scans, snapshots and bookkeeping. MOS, mux and
NEURON have no available retained entry graph in this census, so this is not
a proof about their reuse potential. Counts describe pre-call retained-graph
geometry; repair/selection may replace it, and the counts prove neither a
numerical reuse contract nor executed FLOPs, exclusive CPU cost or speedup.
Decision: `build/mapped-input-closure-census-completed-priority-20261018/decision.json`.

A separate source inspection finds that packed-LU admission rejects every LLC
above64MiB before checking whether its factor payload fits. Historical comment
8d5abf09^ records an empirical preference for native storage on large caches;
the current comment retains only the fit rationale. The next isolated candidate
will remove that absolute capacity ceiling and retain an overflow-safe payload
budget that is monotone in available cache. Existing packed authority must also
survive a later admission refusal, using the independently qualified ownership
fix from the row-funded research branch rather than silently writing native
values while stale packed values remain authoritative. No CPU/domain label,
circuit identity, new dimension band, tolerance or verification change is
introduced. This is an unqualified source proposal until fresh lifecycle/full
correctness, API and complete paired cost evidence supports it; all56 peer
cells and separate complete confirmation remain required for the fastest goal.


The packed-LU proposal is isolated on branch
`work/packed-lu-monotonic-cache-root-20261018`, based directly on root24.
Numerical source0db9dd2 removes the absolute LLC ceiling, preserves the existing
five-eighths payload budget with overflow-safe arithmetic, and ports the prior
qualified retained-owner guard. Wider admission executes existing cold deck
allocations and may retain more memory; all preparation must count in cost.
No new allocation call site, layout field, clock or verification change exists.

Fixture commitsfd9fe38 and937bd39 correct two direct-executor test assumptions:
retire/refresh solve metadata around changed numerics, and reestablish a
current packed copy before each tested refusal after a public transpose solve
legitimately synchronizes and retires its preceding owner. The unchanged fd
numeric source's ownership trace records that exact normal/transpose boundary.
Both failed initial packages and the diagnostic abort remain preserved; neither
supplies full correctness evidence. Trace:
`build/packed-lu-monotonic-cache-fixture-ownership-trace-20261018`.

The fresh lifecycle package for937bd39 passes both initial native checks:
621075 exact payload comparisons with overflow and monotonicity assertions,
16 ownership configurations (crews2/8, BTF0/1, scale0/2, normal/transposed
factor frame), four changed updates each, coefficient assertions and normal,
transpose, padded multiple-RHS known solutions, plus reanalysis cleanup;
the existing EGraph startup test also passes. Full44 native and optimized
ASan/UBSan qualification is running in
`build/packed-lu-monotonic-cache-lifecycle-full-runtime-preparation-20261018`.
No full-suite, sanitizer, API, application-cost, peer or promotion claim follows
from this native preflight. The production numerical source remains2754413.


Candidate937bd39 completes the fresh full44 native and44 optimized ASan/UBSan
suites (88 full tests), both initial focused checks each, and independent
artifact review. All43 parent names remain; all20 child edges close and actual
controller/reviewer exits are zero with PIDs absent. Each build executes all16
packed lifecycle contexts, four changed updates each, coefficient assertions,
normal/transpose/padded multiple-RHS known solutions and reanalysis, plus621075
exact payload comparisons with overflow and cache monotonicity checks.
Review: `build/packed-lu-monotonic-cache-lifecycle-full-runtime-review-preparation-20261018`.

A separate negative control restores only the old ownership guard in the same
corrected fixture/source. It aborts at the intended stale packed U-coefficient
assertion when admission is refused; the corrected candidate passes. The actual
control controller exits zero (expected child SIGABRT), all six child edges
close and PIDs are absent; an independent review proves the one-block source
inverse and unchanged fixture. This establishes that specific ownership refusal
counterexample, not performance or broader public reachability. Review:
`build/packed-lu-retained-authority-negative-control-review-preparation-20261018`.

The candidate also completes all28 original API/circuit settings with288650
certificate calls. Independent rechecks preserve waveforms/steps, seven API
wrappers, the committed adapter/shared dependencies and matched plain links.
All62 API/build edges close after actual zero exits and absent controller/child
PIDs. Review: `build/packed-lu-monotonic-cache-api-closure-review-preparation-20261018`.
Decision: `build/packed-lu-monotonic-cache-completed-correctness-20261018/decision.json`.
No uninstrumented application timing or promotion follows from these gates.

The next complete28 diagnostic resolves the shared mapped path label into
actual entry/return counts for lean BTF maps, generic mapped, EGraph, PTS,
compact-map/KLU tournament, serial KLU and block-pool dispatch. It also records
packed admission reasons, current ownership transitions and observed factor/
block size maxima. Source957d93d is isolated on
`work/packed-lu-backend-census-root-20261018`; a multi-module inverse restores
qualified937 exactly. No worker hooks, clocks, allocations, numeric writes,
layout or policy changes are added. Counters include internal candidates and
recovery and can perturb adaptive choices; counts/maxima cannot establish
exclusive CPU cost or speedup. All three public numeric entry totals must
reconcile with the original seven-wrapper profiler, and every original setting
must retain waveform/step/certificate checks. Production source2754413 and
default verification remain unchanged; all56 peer cells and a separate full
confirmation remain required for the fastest goal.


The complete28 backend census957d93d and independent review now close with
288645 certificate calls, exact reconciliation of all three public numeric
entry totals, unchanged waveform/step gates and all62 child edges closed.
Actual controller/reviewer exits are zero and their PIDs are absent. Review:
`build/packed-lu-backend-census-api-review-preparation-20261018`.
Every setting has zero accepted packed admissions and zero observed packed
owners. Eight-thread BJT enters EGraph1145 times per domain and every packed
admission declines as not_repeated: the default adapter supplies no lifecycle
forecast. Correctness-qualified937bd39 remains isolated; no cost campaign is
justified for an unused policy. Decision:
`build/packed-lu-backend-census-completed-priority-20261018/decision.json`.

The counters identify caller-level engines and factor/block maxima, not every
subordinate mapped body or exclusive CPU fractions. In particular the current
census omits lean-single-block, row and serial-tail bodies and a mapped-BTF
post-gate marker. Historical BTF/lean attribution is supporting context only.
RAM, mem_plus and DAC remain recurring mapped work with one dominant block;
PTS probes decline. Do not turn zero forecast hints into invented future work.

The next isolated exploratory control relaxes EGraph resource floors only
inside unchecked parallel refactors, restoring tuning on return while retaining
capability, pivot, repair and verification contracts. It tests engine potential;
a global floor bypass is not a production selection policy. Any eventual
automatic challenger must fund complete-cycle exploration from observed work
and demonstrate general selection, complete application cost and all56 peer
cells followed by separate full confirmation. Production numerical2754413,
default verification and tolerances remain unchanged.


Research controlc059cc3 on `work/egraph-resource-floor-control-root-20261018`
wraps the exact937 dispatcher body. Only unchecked parallel calls temporarily
set the existing EGraph size/work floors to zero; both tuning values restore
after every body return. Serial/checked dispatch and all structural capability,
pivot, repair and verification gates remain. No new clock, allocator, layout,
arithmetic, tolerance, forecast or case selector is added. This is an isolated
engine-potential experiment, not a production admission policy.

Full44 native and44 optimized ASan/UBSan tests plus four initial focused checks
pass. The startup fixture proves supported EGraph dispatch with deliberately
blocking floors, restoration, checked fallback and normal/transpose known
solutions for crews2/8. Independent review proves the exact dispatcher body,
all other source/include bytes against937, original43 parent test names plus
the packed lifecycle fixture, and all20 child edges closed after actual zero
controller/reviewer exits and absent PIDs. Review:
`build/egraph-resource-floor-control-full-runtime-review-preparation-20261018`.

All28 original API settings independently qualify with288650 certificates,
matching waveforms/steps, matched plain links, unchanged adapter/dependencies
and original seven wrappers; all62 child edges close. Review:
`build/egraph-resource-floor-control-api-review-preparation-20261018`.
Eight-thread mem_plus1200, DAC514, RAM5796 and NEURON7740 updates per domain
now select EGraph; BJT remains EGraph and small MOS/mux retain their earlier
paths. Counts establish coverage, not speed. Decision:
`build/egraph-resource-floor-control-completed-correctness-20261018/decision.json`.
The complete238 uninstrumented paired cost gate is next. Do not promote static
floor removal: a useful engine would still need a general observed-work funded
selection policy and complete cost validation. Production2754413 remains;
all56 peer cells and separate complete confirmation remain unmet.


The complete238 scoped-floor cost campaign now closes after actual controller
exit0, all476 child edges reconstructed and controller/child PIDs absent. Its
independent review also exits0 with PID absent, rechecks all waveforms, exact
frozen binaries, source pins and the complete unique balanced28-cell grid;
all238 observations qualify. Review:
`build/egraph-resource-floor-control-spice-pair-review-preparation-20261018`.
Decision: `build/egraph-resource-floor-control-completed-cost-decision-20261018/decision.json`.

Eight-thread candidate/root median ratios (four pairs each) are:

| Circuit | LLC96 | LLC32 |
| --- | ---: | ---: |
| RAM | 0.695093 | 0.771036 |
| mem_plus | 0.636371 | 0.755676 |
| DAC | 0.756580 | 0.795357 |
| NEURON | 2.805003 | 2.882453 |
| mux | 0.972556 | 0.994804 |
| MOS | 0.992941 | 0.982933 |
| BJT | 0.990775 | 0.996735 |

All four pairs improve on parallel RAM, mem_plus and DAC in each domain; all
four worsen on NEURON. Serial medians are near parity/mixed (0.988515--1.014617),
not a serial performance feature. Retain every cell rather than promoting an
aggregate. One complete campaign establishes this measured control outcome;
it neither confirms peer superiority nor justifies unconditional floor removal.

Implement the next candidate as a complete-cycle selector. Its allowance must
come from successful observed external work or explicit hints; trials cannot
finance themselves. Scope the floor override to the selected call and restore
tuning, charge cold preparation/comparisons/cleanup, compare verified cycles
using the same actual matrix/RHS, restore the incumbent after failure or budget
exhaustion, and discard evidence when the incumbent, factor geometry or RHS
frame changes. Preserve invalid-argument-before-mutation, all RHS shapes and
accuracy policies. Reuse qualified integer budget accounting where possible;
do not import unrelated unmeasured row-policy changes. The control remains
isolated and production2754413/default verification/tolerances stay unchanged.
All56 peer cells and a separate complete confirmation remain mandatory and
unmet for the fastest goal.


The selector is isolated as0a41f32 on
`work/egraph-observed-cycle-selection-root-20261018`, based on qualifiedc059.
The scoped override now requires a caller-owned cycle flag. The combined API
uses its validated, saved RHS and existing output scratch for complete cycles
against the same actual matrix/RHS. Every arm calls the unchanged caller
certificate body through a borrowed-RHS split; the trial returns the incumbent
last. A failed adopted numeric retires further optional exploration and follows
mandatory incumbent/full-factor recovery. Unsupported and unfunded calls retain
the ordinary path; all RHS layouts and accuracy policies retain their contract.

A cold tail state reuses the qualified integer exploration ledger exactly.
Successful existing factor/update/solve measurements fund at most2% optional
work; trials cannot fund themselves. The fixed lease is checked between whole
verified cycles, never inside numeric kernels: this is cooperative cancellation,
not a hard latency/RSS guarantee. Cold preparation, comparisons, certificates,
transition publication, cleanup and actual overshoot remain debt; only one
useful incumbent cycle can be credited. Four verified pairs and the existing
savings/payback margins govern adoption. Epoch, engine, crew, orientation,
accuracy and RHS-count changes discard comparison evidence without refunding
spent work. No case, domain, dimension band or invented forecast selects an arm.
The selector adds no allocation site or unfunded sampling clock; its appended
private state preserves the established hot solver prefix.

The native preflight reports all46 tests and four initial focused checks
passing, including real known solutions for crews2/8 under injected policy
clocks: adoption, measured loss, budget overshoot, clock failure, unsupported
challenger, incumbent restoration, invalid argument, padded RHS and epoch/frame
invalidation. The copied budget fixture covers an independent conversion oracle,
FP state, shared funding, stale leases, reset and overflow. Injected history and
clocks test branch correctness, not speed. The optimized ASan/UBSan run is still
active in `build/egraph-observed-cycle-selection-full-runtime-preparation-20261018`.
Actual complete runtime closure, independent artifact review, complete28 API,
application cost, full peer grid and separate confirmation are still required;
no promotion or fastest claim follows from this preflight. Production2754413,
default verification and tolerances remain unchanged.


Selector0a41f32 completes all46 native and46 optimized ASan/UBSan tests plus
eight focused checks and independent review. All44 parent names, hot-field
prefix bytes, exact dispatcher/certificate bodies, public headers and worker
modules are preserved; all20 child edges close after actual zero exits and
absent controller/reviewer/child PIDs. Review:
`build/egraph-observed-cycle-selection-full-runtime-review-preparation-20261018`.
All28 original API/circuit settings also independently qualify with288706
certificate calls, matching waveforms/steps and exact matched binaries;
all62 child edges close after actual zero exits/PID absence. Review:
`build/egraph-observed-cycle-selection-api-review-preparation-20261018`.

Observed eight-thread public EGraph counts are RAM5289/5301 and mem_plus694/730
on LLC32/96, respectively, out of5796 and1200 updates. NEURON remains on its
prior mapped/KLU paths; DAC stays mapped for all514 updates. BJT remains EGraph;
small MOS/mux retain their prior paths. These counts establish automatic
selection coverage and correctness, not uninstrumented performance.
Decision: `build/egraph-observed-cycle-selection-completed-correctness-20261018/decision.json`.

Before spending a cost campaign on0a41, source review finds a redundant leading
incumbent cycle in each optional comparison. The arms already refactor, solve
and certify the same actual A/RHS, and the final incumbent restores current
state and supplies the caller answer. The extra leading incumbent answer is
overwritten. Candidate062674d removes only that pass and adjusts the private
injected clock schedule, retaining publication before the challenger, mandatory
final incumbent/full-factor recovery, budget debt, evidence invalidation and
all policy thresholds. The previous successful external cycle identifies the
incumbent engine; a changed final engine discards the pair. No corpus selector,
forecast or threshold retuning is introduced. Fresh full92 plus8 focused gates
are running in `build/egraph-paired-cycle-selection-full-runtime-preparation-20261018`.
No cost campaign ran for0a41, and062674d has no qualification yet. Production
2754413/default verification/tolerances remain; all56 peer cells and separate
complete confirmation remain required and unmet.


Candidate062674d completes fresh full46 native and46 optimized ASan/UBSan tests
plus eight focused checks, with actual controller exit0, all20 child edges
closed and PIDs absent. Independent review exits0/PID absent and preserves all46
parent names, exact unchanged sources except the pair sequence, unchanged
thresholds, original dispatcher/certificate bodies, hot prefix and public/worker
interfaces. The fixture changes only its injected clock schedule. Review:
`build/egraph-paired-cycle-selection-full-runtime-review-preparation-20261018`.

All28 API settings independently qualify with288679 certificates and unchanged
waveforms, steps, seven wrappers, matched plain links and adapter/dependencies.
All62 child edges close after actual zero exits/PID absence. Review:
`build/egraph-paired-cycle-selection-api-review-preparation-20261018`.
Eight-thread public EGraph counts (LLC32/96) are mem_plus848/884 of1200, DAC0/180
of514, and RAM0/5467 of5796. NEURON retains mapped/KLU, BJT retains EGraph and
MOS/mux retain their prior paths. These are single instrumented-run selection
observations, not speed or reliable adoption guarantees; retain the varying RAM
and DAC outcomes explicitly. No static domain/circuit admission is justified.
Decision: `build/egraph-paired-cycle-selection-completed-correctness-20261018/decision.json`.
The next gate is the original full238 balanced uninstrumented cost comparison
against root24. No source changes, filtering or simultaneous heavy builds/tests
are permitted during it. Production2754413/default verification/tolerances stay;
all56 peer cells and separate full confirmation remain unmet.


Candidate062674d completes the full238 uninstrumented comparison; all observations
qualify and all28 balanced cells retain four pairs. Runtime and independent
review exit0 with controller/child PIDs absent. Review:
`build/egraph-paired-cycle-selection-spice-pair-review-preparation-20261018`.
Candidate/root24 median linear-solver time ratios (lower is faster):

| Circuit | LLC96, 1 thread | LLC96, 8 threads | LLC32, 1 thread | LLC32, 8 threads |
|---|---:|---:|---:|---:|
| CircuitSim90/MOS2_LARGE/ram2k.cir | 0.996809 | 0.709257 | 0.997849 | 0.779632 |
| NEURON/rallpack3.cir | 0.998197 | 1.003646 | 0.996073 | 1.004440 |
| CircuitSim90/MOS2_LARGE/dac.cir | 0.989919 | 0.914918 | 0.990565 | 0.927865 |
| CircuitSim90/MOS2/mux8.cir | 0.975897 | 1.002080 | 0.976928 | 0.986557 |
| CircuitSim90/MOS2_LARGE/mem_plus.cir | 0.990059 | 0.728418 | 0.986751 | 0.829487 |
| MOS13_IC/invert50_mos1.cir | 1.016802 | 1.000484 | 1.015647 | 0.999660 |
| PDE_2D_BJT/invertbjt.cir | 1.002099 | 1.004419 | 1.000373 | 1.001637 |

Hold promotion: serial MOS is slower in both domains. The selector preserves
the large parallel RAM/mem gains and avoids the rejected static-control NEURON
slowdown. Decision: `build/egraph-paired-cycle-selection-completed-cost-decision-20261018/decision.json`.
Candidatee4db221 only moves cold repeat checks after hot eligibility checks and
guards budget reads/debit accounting with outer parallel-call eligibility.
These are general inactive-path changes, with unchanged numerical bodies,
selector thresholds, tests and public interfaces. Avoidable source reads are
not established as the exclusive cause of MOS timing. Fresh full92 plus8
qualification precedes API and full-grid cost checks. Production2754413, default
verification and tolerances remain. All56 peer cells and separate full
confirmation remain required and unmet.


Candidatee4db221 passes all46 native and46 optimized ASan/UBSan tests plus eight
focused checks. Independent review confirms the exact two-file guard change,
unchanged tests, numerical bodies, policy thresholds, hot prefix and public/worker
interfaces. All20 runtime child edges close after actual exit0 and absent PIDs.
Review: `build/egraph-inactive-cycle-selection-full-runtime-review-preparation-20261018`.
All28 API settings also independently qualify with matching waveforms/steps,
original seven wrappers and matched plain links; all62 API/build child edges
close after actual exit0 and absent PIDs. Review:
`build/egraph-inactive-cycle-selection-api-review-preparation-20261018`.
Decision: `build/egraph-inactive-cycle-selection-completed-correctness-20261018/decision.json`.
The original full238 balanced uninstrumented cost grid is next, with source
frozen and no concurrent heavy work. No performance, peer or promotion claim
follows from these correctness gates. Production2754413 and verification remain.


Candidatee4db221 completes full238 cost and independent review with actual
exit0/controller/reviewer/child PIDs absent, all238 qualified and all28 balanced
cells complete. Review: `build/egraph-inactive-cycle-selection-spice-pair-review-preparation-20261018`.
Candidate/root24 median linear-solver time ratios (lower is faster):

| Circuit | LLC96, 1 thread | LLC96, 8 threads | LLC32, 1 thread | LLC32, 8 threads |
|---|---:|---:|---:|---:|
| CircuitSim90/MOS2_LARGE/ram2k.cir | 1.000673 | 0.720816 | 0.999118 | 0.761054 |
| NEURON/rallpack3.cir | 0.999681 | 1.001399 | 0.995198 | 1.003255 |
| CircuitSim90/MOS2_LARGE/dac.cir | 0.990300 | 0.925465 | 0.992233 | 0.935977 |
| CircuitSim90/MOS2/mux8.cir | 0.993011 | 0.986122 | 0.994921 | 0.993080 |
| CircuitSim90/MOS2_LARGE/mem_plus.cir | 0.991566 | 0.742242 | 1.002442 | 0.831727 |
| MOS13_IC/invert50_mos1.cir | 1.021179 | 1.006949 | 1.001189 | 1.000096 |
| PDE_2D_BJT/invertbjt.cir | 1.003531 | 1.000564 | 1.003527 | 0.991888 |

Hold promotion: serial MOS remains mixed across domains and is2.1% slower in
LLC96. The guard cleanup does not establish regression removal. Parallel RAM,
mem_plus and DAC gains persist; NEURON avoids the static-floor slowdown.
Decision: `build/egraph-inactive-cycle-selection-completed-cost-decision-20261018/decision.json`.
Native nm symbol sizes show kls_refactor_solve grows from490 bytes in root24
to5722 in e4db221, with no standalone kls_egraph_cycle_try symbol: the optional
selector is inlined into the common API. This establishes code growth, not
exclusive timing causality. The next general candidate extracts unchanged
parallel orchestration into a noinline helper, leaving the ordinary serial
body/finish sequence in the public wrapper. No circuit, domain, dimension or
future-use selector is added; numerical kernels and thresholds are unchanged.
Fresh full92 plus8 runtime,28 API and complete cost qualification are required.
Production2754413/default verification/tolerances remain; all56 peer cells and
separate full confirmation remain required and unmet.


Candidate3ba29c9 passes full46 native plus46 optimized ASan/UBSan tests and eight
focused checks. Independent source review proves exact extraction of the
qualified parallel orchestration, unchanged tests/numerical bodies/thresholds
and interfaces, plus native noinline helper retention. Native public wrapper
size is578 bytes versus5722 in e4db221 and490 in root24; this is code-size
evidence, not performance proof. All20 runtime child edges close after actual
exit0 and absent PIDs. Review:
`build/egraph-outlined-cycle-selection-full-runtime-review-preparation-20261018`.
All28 API settings independently qualify, retaining waveforms/steps, seven
original wrappers and matched plain links. All62 API/build child edges close
after actual exit0 and absent PIDs. Review:
`build/egraph-outlined-cycle-selection-api-review-preparation-20261018`.
Decision: `build/egraph-outlined-cycle-selection-completed-correctness-20261018/decision.json`.
Next is the original full238 balanced uninstrumented comparison, without source
changes or concurrent heavy work. Production2754413/default verification and
tolerances remain; all56 peer cells and separate full confirmation are unmet.


Candidate3ba29c9 completes full238 cost and independent review with actual
exit0 and absent controller/reviewer/child PIDs. All238 observations qualify;
all28 balanced cells retain four pairs. Review:
`build/egraph-outlined-cycle-selection-spice-pair-review-preparation-20261018`.
Candidate/root24 median linear-solver time ratios (lower is faster):

| Circuit | LLC96, 1 thread | LLC96, 8 threads | LLC32, 1 thread | LLC32, 8 threads |
|---|---:|---:|---:|---:|
| CircuitSim90/MOS2_LARGE/ram2k.cir | 0.997192 | 0.715738 | 1.000670 | 0.763937 |
| NEURON/rallpack3.cir | 1.003853 | 0.939570 | 1.001043 | 0.937859 |
| CircuitSim90/MOS2_LARGE/dac.cir | 0.986370 | 0.941770 | 0.983333 | 0.962339 |
| CircuitSim90/MOS2/mux8.cir | 0.990363 | 0.994252 | 1.006680 | 0.994212 |
| CircuitSim90/MOS2_LARGE/mem_plus.cir | 0.988053 | 0.723333 | 0.993109 | 0.833215 |
| MOS13_IC/invert50_mos1.cir | 1.009157 | 0.997805 | 0.990065 | 1.007148 |
| PDE_2D_BJT/invertbjt.cir | 1.001494 | 0.981441 | 1.003873 | 0.996218 |

Parallel RAM/mem/DAC gains persist. Serial MOS is0.9% slower LLC96 and1.0%
faster LLC32; small other regressions and variable NEURON pairs remain visible.
Native wrapper shrink is established, but universal improvement and exclusive
timing causality are not. Do not tune thresholds/case selectors to these cells.
Decision: `build/egraph-outlined-cycle-selection-completed-cost-decision-20261018/decision.json`.
Retain3ba29c9 as a qualified research candidate for the complete four-arm peer
diagnostic: root24/3ba29c9/CKTSO/SubtreeLU, all7 circuits, both domains, threads1/8,
four balanced passes plus14 fresh KLU references (462 observations). Preserve
all peer failures/timeouts and original1800-second limits; no inherited ratios.
This is a diagnostic, not promotion or confirmation. All56 qualified peer cells
and a separate full confirmation remain required and unmet. Production2754413,
default verification and tolerances remain unchanged.


Candidate3ba29c9 completes the original full462 four-arm peer diagnostic. Actual
runtime and independent reviewer exit0; controller/reviewer and all child groups
are absent, with all924 edges closed. All112 candidate,112 root24 baseline and
14 KLU observations qualify. Three retained SubtreeLU MOS eight-thread timeouts
leave two peer cells unqualified. Candidate wins25 peer cells and loses29; no
fastest or promotion claim follows. Review:
`build/egraph-outlined-cycle-selection-peer-spice-review-preparation-20261018`.
Median candidate/peer linear time ratios, below1 faster; an unqualified cell
is shown explicitly and never counted as a win:

llc96

| Circuit | CKTSO, 1 thread | SubtreeLU, 1 thread | CKTSO, 8 threads | SubtreeLU, 8 threads |
|---|---:|---:|---:|---:|
| CircuitSim90/MOS2_LARGE/ram2k.cir | 1.053689 | 0.634112 | 1.126695 | 0.889023 |
| NEURON/rallpack3.cir | 1.474429 | 0.796894 | 0.277816 | 1.578153 |
| CircuitSim90/MOS2_LARGE/dac.cir | 1.125498 | 0.600520 | 1.427934 | 1.233986 |
| CircuitSim90/MOS2/mux8.cir | 1.296948 | 0.829545 | 0.263578 | 0.519032 |
| CircuitSim90/MOS2_LARGE/mem_plus.cir | 1.165599 | 0.693703 | 1.142811 | 1.184682 |
| MOS13_IC/invert50_mos1.cir | 1.236322 | 0.870968 | 0.202653 | unqualified |
| PDE_2D_BJT/invertbjt.cir | 1.048992 | 0.675393 | 1.226773 | 1.309512 |

llc32

| Circuit | CKTSO, 1 thread | SubtreeLU, 1 thread | CKTSO, 8 threads | SubtreeLU, 8 threads |
|---|---:|---:|---:|---:|
| CircuitSim90/MOS2_LARGE/ram2k.cir | 1.083011 | 0.543024 | 0.980897 | 0.831441 |
| NEURON/rallpack3.cir | 1.431573 | 0.779637 | 0.307344 | 1.657988 |
| CircuitSim90/MOS2_LARGE/dac.cir | 1.161125 | 0.612947 | 1.060051 | 1.285125 |
| CircuitSim90/MOS2/mux8.cir | 1.294718 | 0.819214 | 0.240175 | 0.503975 |
| CircuitSim90/MOS2_LARGE/mem_plus.cir | 1.342029 | 0.674287 | 1.136632 | 1.432427 |
| MOS13_IC/invert50_mos1.cir | 1.244515 | 0.869051 | 0.192811 | unqualified |
| PDE_2D_BJT/invertbjt.cir | 1.157948 | 0.746933 | 1.192846 | 1.417201 |

All serial CKTSO cells remain losses. The largest remaining qualified gap is
NEURON against SubtreeLU at eight threads in LLC32: ratio1.657988, requiring
39.69% less time to reach parity. This rules out treating further wrapper
cleanup as a sufficient route to the goal. Decision:
`build/egraph-outlined-cycle-selection-peer-completed-decision-20261018/decision.json`.
Production2754413/default certification/tolerances remain unchanged.

### Explicit Xyce accuracy policy after peer-policy audit (2026-10-18)

The CKTSO user guide Version20260321 §5.6 documents forward/backward
substitution. Its demo computes a residual after Solve and describes that
check as unnecessary in circuit simulation. SubtreeLU's demo checks a known
solution externally; the paper's accuracy experiments are offline validation.
The Xyce third-party adapter adds no public residual certificate. Both peer
numerical libraries are closed binaries here, so undocumented internal checks
are not proven absent. These findings support a separately reported
application-managed experiment, not a claim of identical guarantees or a
production default change. Primary-source hashes and limits are retained in
`build/peer-accuracy-policy-audit-20261018/audit.json`.

Xyce research commit `67d9590` exposes `KLS_ACCURACY_POLICY` with the existing
COMPONENTWISE_BACKWARD_ERROR default and explicit STRICT_RHS_L2 and
APPLICATION_MANAGED modes. Invalid values fail before altering retained
options/factors. The policy is applied before analysis of each new solver;
changing options follows normal factor destruction/reanalysis. The adapter's
object layout and KLS numerical source are unchanged. Application-managed
mode retains finite/kernel error checks but provides no certified residual
bound or certificate-driven correction; integration tolerances are unchanged.

Fresh real Epetra adapter tests passed in native and optimized ASan/UBSan
builds: 192 fixtures and 1,344 solves per build, with all 768 analysis-time
policy readbacks checked. Coverage includes omission/explicit defaults,
quoted/lowercase names, rejected invalid requests preserving usable factors,
live policy changes/restoring omission, graph/value changes, refactor on/off,
factor reuse, normal/transpose, one/two RHS and both storage/map variants.
The metadata registration translation unit also compiled in both profiles.
The unchanged 3ba29c9 core's prior 46+46 tests are inherited, not fresh tests.
Xyce support dependencies remain release builds; MPI multi-rank, real-netlist
registration, full SPICE modes, cost comparisons and promotion remain pending.

Artifacts: `build/xyce-accuracy-policy-adapter-qualification-v2-20261018` and
`build/xyce-accuracy-policy-adapter-review-v2-20261018`. Actual controller and
reviewer exits and absent PIDs were checked before accepting results. The
original failed launcher and reviewer packages are preserved: one retained a
compiler in its run command; the other used an incorrect header-compare path.
Neither failure was a numerical test failure. No fastest claim is supported.

### Retained factor streams in the lean refactor candidate (2026-10-18)

Research commit `3230492` changes only the lean single-block map32 numerical
function and its private transpose fixture. When the existing compact solve
cache is ready, refactorization borrows its validated uint16 indices and
native factor-value pointers. It allocates no new representation, changes no
admission policy, and retains the native factor walk when that cache is
unavailable. Producer and destination order, pivot handling, current-value
scaling and the caller's accuracy policy are preserved. Numeric replacement
and panel teardown already invalidate these solve caches; block repair also
invalidates them. Public headers and solver object layout are unchanged.

Fresh focused native, optimized ASan/UBSan and generic-main-kernel checks
passed, followed by a separate review with actual exit/PID closure and hash
revalidation. The fixture compares initialized factor values, diagonals,
scales and fused partial results exactly against the native fallback. It also
checks known solutions across changed values, full factor rebuilds,
reanalysis, wide-index fallback and singular rejection. The generic context
uses unchanged native dependency/worker archives and is not a full portable
library qualification. Full core, Xyce API, cost and peer gates remain pending.
No speed benefit or promotion is established.

Artifacts: `build/lean-retained-factor-stream-focused-20261018` and
`build/lean-retained-factor-stream-focused-review-20261018`. These correctness
checks overlapped the live strict-policy diagnostic on logical CPU24. A later
topology check found CPU24 is CPU8's SMT sibling: disjoint masks did not give
physical isolation. All overlapping diagnostic timings are excluded from
speed comparisons. Earlier complete default/componentwise grids preceded
this overlap, but still require the campaign's independent review before
interpretation. The original executed packages remain unchanged; correction
and topology evidence are in
`build/lean-retained-factor-stream-overlap-topology-correction-20261018`.
Broader builds and uninstrumented timing are deferred until the live campaign
closes. Its source/binary/input hashes were rechecked and remain unchanged.

### Complete accuracy-mode checks and current phase profiles (2026-10-18)

The full112 campaign and independent review both exited zero, with controller
and reviewer absence, all224 runtime child edges and process-group closure
checked before interpretation. Qualification was28/28 for omitted/default,
28/28 for explicit componentwise,28/28 for application-managed, and20/28
for strict RHS-relative L2. All eight failures are retained: strict MOS at
both crews/domains logged five linear-solve failures despite matching step
counts and passing upstream waveform comparison; strict memory runs logged
a linear failure and changed step counts and waveforms. No failure was
converted into a win or repaired by weakening tolerances.

Application-managed remains explicitly selectable and separately labelled.
These full-circuit checks do not provide its missing certified residual
contract. The default remains COMPONENTWISE_BACKWARD_ERROR. Mode readbacks
were checked before every analysis; netlists differ only in the requested
accuracy option. Original complete circuits, waveform and step/end gates
remain in force.

The complete28-cell default phase grid below uses exclusive time within
caller API ancestry, including worker-zero execution and waits. Percentages
are diagnostic, rounded to integers; clock instrumentation can perturb
adaptive selection. Background worker intervals overlap and are not removable
wall time. All112 are instrumented diagnostics, never uninstrumented speed
evidence. Strict comparative timings are explicitly excluded because of
the recorded SMT overlap. The default and explicit componentwise grids
preceded that overlap. Each cell lists its two largest instrumented scopes.

| Case | LLC96,1 | LLC96,8 | LLC32,1 | LLC32,8 |
| --- | --- | --- | --- | --- |
| invertbjt | SNB 88%, transpose solve 8% | EGraph 72%, transpose solve 20% | SNB 84%, transpose solve 11% | EGraph 71%, transpose solve 21% |
| invert50_mos1 | refactor 33%, solve 24% | refactor 34%, certificate 23% | refactor 34%, certificate 23% | refactor 34%, certificate 23% |
| mem_plus | mapped 63%, certificate 16% | EGraph 30%, mapped 25% | mapped 56%, certificate 18% | EGraph 24%, mapped 23% |
| mux8 | refactor 36%, certificate 22% | refactor 36%, certificate 22% | refactor 37%, certificate 22% | refactor 37%, certificate 22% |
| dac | mapped 55%, certificate 18% | mapped 43%, solve 14% | mapped 54%, certificate 18% | mapped 43%, solve 14% |
| rallpack3 | lean 45%, certificate 31% | lean 45%, certificate 31% | lean 45%, certificate 32% | lean 46%, certificate 32% |
| ram2k | mapped 58%, certificate 19% | EGraph 41%, solve 24% | mapped 50%, certificate 19% | EGraph 37%, solve 21% |

This supports qualifying the existing lean factor-stream candidate3230492
and isolated SNB row-tile candidate. It also replaces the older assumption
that the eight-thread RAM/memory path remains wholly mapped: the current
selector exposes substantial EGraph work, alongside solve/snapshot costs.
Next representation/scheduling work must use a fresh actual-owner census;
the previously unused packed-policy admission stays parked meanwhile.

Artifacts: `build/outlined-accuracy-policy-full-modes-profile-20261018`,
`build/outlined-accuracy-policy-full-modes-review-20261018`, and
`build/outlined-accuracy-policy-completed-decision-20261018/decision.json`.
No numerical promotion, uninstrumented mode speedup or peer superiority is
established. All56 qualified peer wins plus separate confirmation remain
unproven.

### SNB row tiling focused qualification (2026-10-18)

Research commit `36d7ff0` prepares a per-edge row bound targeting a32KiB
streaming working set and shares a tiled product/scatter consumer between
serial and cooperative dense updates. It retains the original four-column
product body, producer/destination reduction order, scratch allocations,
allocator budgets and crew-barrier sites. Bounds depend on operation shape,
with no circuit or CPU-domain selector. The calibration fixture prepares the
same cached bound outside its timed kernel. Full calibrator build/validation
is still pending.

All12 focused contexts passed and were independently reviewed: product,
builder, masked-bound and triangular-update fixtures in native, optimized
ASan/UBSan and generic-main-kernel builds. The new oracle compares exact
outputs, standard floating-point exception flags and scratch canaries against
the original full-range product/scatter, across tile boundaries, partial
worker ranges, gapped subcolumns, absent destinations, four rounding modes
and available FTZ/DAZ states. Builder checks require actual edge coverage and
verify that early scatter cannot overwrite the current upper input block.
The product body was independently compared with the parent, and all original
barrier sites remain. The generic context inherits native worker/dependency
archives and is not full portable-library qualification.

The review classifies every compiler warning. A new GCC inliner array-bounds
warning occurs only in the triangular-only fixture's unused scratch path:
producer below=0, the caller returns before the dense consumer, and the
consumer independently returns for the empty range. All three runtime
contexts passed that guard-page/exception fixture. Existing feature-macro and
allocation-range warnings originate in unchanged core or allocator-wrapper
bodies. No warning was suppressed. Evidence and limits are in
`build/snb-cache-row-tile-focused-preparation-20261018` and
`build/snb-cache-row-tile-focused-review-20261018`.

Full core, portable library, Xyce API, complete cost/peer grids and confirmation
remain pending. No speed gain or promotion is established. The two numerical
research commits remain separate for causal cost comparisons.

### Portable batch FMA dispatch and corrected candidates (2026-10-18)

The first complete retained-factor-stream run at `3230492` passed all46
native and all46 optimized ASan/UBSan tests, but failed one of46 fully
portable tests: the unchanged exact batch-fusion comparison. The same
failure reproduces on parent `3ba29c9` with the portable flags and unchanged
worker/dependency archives. The137/138 result is retained as a failed full
qualification, not attributed to compact streams or counted as a pass.

On this GCC/x86 build, the independent cached consumer has a dispatched
x86-64-v4 clone with fused multiply/add instructions; the pair and multi
consumers previously had only generic separate multiply/add instructions.
Commit `5d1dfad` adds matching conditional target clones to those two fused
consumers. Arithmetic bodies, reduction order, the exact test comparison,
tolerances, admission policy and budgets are unchanged. Native FMA builds
retain direct calls. This addresses a general portable dispatch mismatch;
it has no circuit, matrix-size or benchmark selector.

The original exact batch test passes native, optimized ASan/UBSan and
portable focused contexts. Independent review checks the unchanged source
bodies, native direct symbols, portable fused clones/instructions, actual
process closure and hashes, and preserves the failing parent control.
Evidence is in `build/portable-batch-parent-reproduction-20261018`,
`build/lean-retained-factor-stream-portable-failure-review-v2-20261018`,
`build/portable-batch-fma-dispatch-focused-20261018` and
`build/portable-batch-fma-dispatch-focused-review-20261018`.

The fix is committed separately, then applied independently to compact
streams (`c2b2956`) and SNB row tiling (`1311dc7`). Fresh full138 packages
use frozen source archives and preserve the original failed/executed
packages. Full core qualification, API/circuit checks, held-out matrices
and uninstrumented comparisons remain required before either performance
change can be promoted. Certified production/default policy is unchanged.

### Retained-stream full correctness and API qualification (2026-10-18)

Candidate `c2b2956` passes all138 complete core tests:46 native,46 optimized
ASan/UBSan and46 fully portable, plus four initial guards in each profile.
The frozen source combines the separately committed retained compact stream
consumer with portable batch target consistency; the original exact batch
comparison is unchanged. All30 child edges/process groups and the controller
closed before independent artifact review. The reviewer then exited zero,
was absent, and all checked source/archive/binary hashes were revalidated.
Evidence: `build/lean-factor-stream-portable-fma-full-runtime-preparation-20261018`
and `build/lean-factor-stream-portable-fma-full-runtime-review-20261018`.

All28 original Xyce circuit/API settings also pass independent waveform,
step/end, linear-failure and seven-wrapper checks, with288678 current-operator
certificate calls. Matched plain links differ only in core/output paths from
parent3ba29c9. All62 runtime/build child edges and groups closed before the
reviewer exited zero and its hashes were rechecked. These instrumented runs
overlap explicitly recorded correctness/diagnostic work; their timings are
excluded from comparative performance evidence. Evidence:
`build/lean-factor-stream-portable-fma-api-preparation-20261018` and
`build/lean-factor-stream-portable-fma-api-review-20261018`.

The SNB candidate `1311dc7` remains separate. Its first launcher stopped
before controller startup/build because its required overlap-policy record
was missing. That attempt is preserved. The fresh v2 package records the
only allowed overlapping instrumented/correctness jobs and retains the same
frozen numerical source,46 tests, flags, budgets and qualification gates.
It makes no timing claims. Complete core/API review and uninstrumented cost
qualification remain pending for SNB. Neither candidate is promoted.

### SubtreeLU complete-call capture qualification (2026-10-18)

Commit `40a0ac9` adds optional diagnostic capture, strict replay and a small
fixture; it changes neither KLS nor the peer adapter/library. Eight fixture
configurations qualify original-library statuses, input/output bytes,
immediate FP state and errno, singular recovery and in-place solves. Strict
replay preserves original rejections/errors,26 malformed streams invoke no
solver API, and trailing entries remain explicitly unreturned. Supplemental
checks cover empty-nnz pointer flags. See `docs/subtreelu_callstream.md`.

The full original MOS diagnostic completes all16 predeclared runs: API team8,
OpenMP1/8, both domains, four trials. Independent review verifies all16
generations, raw stream/log boundaries, file extents and all11 completed
waveforms/step counts. Five capped runs retain a complete entry to refactor
call56 with no return. Each preserves refactor53 returning E_SINGULAR(-3),
factor54 success, and solve55 success. Stalls occur with both OpenMP settings;
changing that environment alone is not a fix. These are120-second diagnostic
caps with instrumentation/concurrent correctness work, not the original
1800-second competitive observations or a speed result.

The original review attempt omitted the independent verifier's input-netlist
copy and stopped. Its failure is retained; v2 only repairs verifier setup,
then rechecks the complete set. Runtime and reviewer zero exits, PID/group
absence and hashes are checked before interpretation. Evidence:
`build/subtreelu-mos13-callstream-diagnostic-20261018` and
`build/subtreelu-mos13-callstream-review-v2-20261018`.

A separate predeclared standalone diagnostic replays the lexicographically
first captured trailing stream with exact completed statuses/outputs, both
domains and OpenMP1/8. A fresh return from the captured unreturned call is
labelled as a new result, never a captured success. Constructor history,
cross-call addresses and internal scheduling remain outside replay scope.
No library/adapter recovery change, timeout win, peer promotion or fastest
claim follows from capture alone. All56 peer wins plus separate confirmation
remain required and unmet.

### Dense-help stack-context lifetime repair (2026-10-19)

The isolated SNB candidate `1311dc7` passed its46 native tests, then its
optimized ASan/UBSan run failed `kls_smoke` with a stack-use-after-return in
`kls_dense_help_try`. The helper inspected a published stack context before
announcing itself in `dense_help_inflight`; the publisher could withdraw,
observe no readers, and return first. That source is byte-identical to parent
`3ba29c9`, so this is an inherited protocol defect, without evidence that SNB
arithmetic caused it. Portable execution did not start. The failed package,
actual nonzero exit, child/process-group closure and independent failure
review are retained under
`build/snb-row-tile-portable-fma-full-runtime-v2-preparation-20261018` and
`build/dense-help-lifetime-failure-review-20261018`.

Separate commit `036aba3`, based on portable-dispatch fix `5d1dfad`, pins a
reader before any context inspection or help-workspace allocation. Reader
pins, active publication/withdrawal and all three publisher drains share a
sequentially consistent order. This prevents the reader and publisher from
both missing the other's announcement, which acquire/release operations on
separate atomics alone do not prevent. Arithmetic, row claims, dependency
rules, scratch budgets and public layouts are unchanged; the inactive fast
return remains.

A Linux regression fixture pauses the actual helper immediately after its
context-pointer load. It checks that the publisher cannot return until the
reader releases its pin. The parent negative control deterministically
reproduces ASan stack-use-after-return; fixed native, optimized ASan/UBSan
and fully portable contexts pass. The fixture adds no production hooks or
sleep-based scheduling assumptions. Focused execution and independent
source/artifact review are closed, including actual zero exits, process
absence and hash rechecks. Evidence:
`build/dense-help-context-lifetime-focused-v2-20261018` and
`build/dense-help-context-lifetime-focused-review-v2-20261019`. Earlier fixture
compile and reviewer warning-classification failures remain preserved in
separate packages.

The lifetime fix is applied separately to the retained-stream candidate
(`7addc94`) and SNB row-tile candidate (`4534485`). Fresh frozen full-core
qualification of corrected baseline `036aba3` and both candidates is running:
47 tests per native, optimized ASan/UBSan and fully portable profile, with
six initial guards per profile. Full141 results, independent review, fresh
28-setting Xyce API checks and separate uninstrumented238-observation paired
cost remain pending. The two candidate launchers initially rejected a
relative script path before starting any controller/build; those launch
failures are preserved, and fresh v2 candidate packages correct only that
check. No failed numerical run is retried in place. The old unfixed cost
package must not be used. Neither performance candidate is promoted and the
certified componentwise default remains enabled.

### SubtreeLU recovery candidate diagnostic (2026-10-19)

The standalone replay of an original trailing call stream reproduced one
stall after every completed status/output had matched. Six fresh executions
returned E_SINGULAR at the originally unreturned call; one stopped at an
earlier solution-byte mismatch. All outcomes remain explicit, without a
tolerance waiver. Replay excludes constructor history, cross-call addresses
and internal scheduling, so it does not uniquely identify the internal cause.
Evidence: `build/subtreelu-captured-stall-replay-20261018` and its independent
`build/subtreelu-captured-stall-replay-review-20261018`.

Separate Xyce candidate `343c527` rebuilds SubtreeLU analysis before full
factorization whenever a fixed-pattern refactor fails. It preserves the
current matrix/options and applies to every such failure, with no circuit,
size or thread selector. Matched builds differ only in the replacement
adapter object; the baseline executable is byte-identical to the original
peer driver. The same16-run MOS diagnostic now completes16 waveforms with
matching steps and no120-second cap, capturing80 generations including64
recoveries. Execution/build closure and independent reviews are complete:
`build/subtreelu-failed-refactor-reanalysis-build-review-v2-20261018` and
`build/subtreelu-reanalysis-recovery-mos13-review-20261018`.

These instrumented, overlapping diagnostic results establish neither cost
nor peer superiority. Broader native/sanitized adapter and full-circuit
qualification are still required before this recovery candidate can replace
the peer adapter. Original stalls and benchmark timeouts are retained; main
Xyce and the original linked archive remain unchanged. All56 qualified peer
wins and a separate complete confirmation remain unmet.

### Certificate work-balanced partition candidate (2026-10-19)

Inspection of corrected baseline `036aba3` confirms that certificate callbacks
already reuse existing numerical workers. Their contiguous ranges still
contain equal numbers of rows, even when row degrees differ substantially.
Independent candidate `01e2eea` balances an integer prefix of rows plus
nonzeros, using the existing immutable row offsets or CSC column pointers for
the requested direction. A binary search finds each boundary; no row is split,
no new metadata is allocated, and overflow retains the old row partition.
OpenMP uses its actual team size, including dynamic reductions. Per-row
arithmetic, acceptance, deterministic maximum, worker ownership, floating-point
environment handling, calibration and memory budgets are unchanged.

The expanded existing parallel certificate fixture checks exact row coverage,
contiguous ranges, a maximum-work bound, empty worker ranges, an irregular
pattern, both directions, exact serial/parallel certificate bounds and actual
dynamic teams, alongside its original allocator/environment faults, rejection,
correction, lifecycle and independent-context tests. Native, optimized
ASan/UBSan and portable profiles pass; the no-OpenMP portable profile checks
only the unchanged serial fallback. Independent review reconstructs the
original production source byte-for-byte after removing the helper and two
range substitutions, and verifies that the original fixture body is retained.
Runtime/reviewer zero exits, process absence and all hashes are closed.
Evidence: `build/certificate-work-partitions-focused-v2-20261019` and
`build/certificate-work-partitions-focused-review-20261019`.

The initial fixture stopped before partition execution because it assumed
ready eligibility exactly1, while the existing degree-domain proof may return2.
Its allocation also bypassed the fixture's accounting wrapper. The failed
package/worktree is preserved; v2 fixes fixture setup without changing the
production partition bytes. An explicit scope clarification records that the
executed command list contains one expanded parallel fixture in four profiles,
not separate standalone verdict/interval fixtures. Full47-profile integration,
Xyce API, held-out matrices and whole-application cost remain pending. The
candidate is separate from retained streams and SNB; it makes no speed claim
and is not promoted.

### Lifetime-fixed full core and broader recovery qualification (2026-10-19)

Corrected baseline `036aba3`, retained-stream candidate `7addc94` and SNB
row-tile candidate `4534485` each pass all141 full core tests:47 native,
47 optimized ASan/UBSan and47 fully portable, plus18 initial guards per
variant. All runtime controllers and90 child edges/process groups closed.
Each separate reviewer exited zero, was absent, and had every checked
source/archive/test/binary hash revalidated. Reviews also verify exact
composition of the separately focused changes and retention of every46
parent test name plus the new lifetime fixture. Existing warnings remain
visible; their six message classes already occur in the earlier qualified
full138 logs. Evidence is in `build/{baseline,lean-streams,snb-tile}-lifetime-full-runtime-review-20261019`
and `build/lifetime-fixed-core-warning-audit-20261019`.

SubtreeLU recovery `343c527` also passes the broader56 circuit observations:
all seven original circuits, both cache domains, one/eight requested threads,
with both native driver and optimized ASan/UBSan adapter builds. All steps,
completion markers and waveforms qualify; independent review reruns all56
waveform comparisons. The runtime's116 child edges/process groups and its
controller close, then the reviewer exits zero and all checked hashes are
revalidated. The closed peer library and existing application/dependency
objects remain uninstrumented; this is adapter-boundary and circuit
qualification, not a complete vendor memory-safety claim. Evidence:
`build/subtreelu-recovery-full-circuit-correctness-20261019` and
`build/subtreelu-recovery-full-circuit-correctness-review-20261019`.

These are correctness observations with explicitly recorded overlap, without
comparative timing evidence. Fresh full28 Xyce API checks for each corrected
KLS variant have started. The paired238 cost campaigns remain unexecuted;
these require closed API reviews and an otherwise quiet machine. The separate
work-partition candidate `01e2eea` has only its focused qualification; its
full-core source package is prepared and queued after these comparisons.
Production/default accuracy and numerical source are unchanged. No candidate
is promoted and the56-peer-win plus separate-confirmation goal remains unmet.
