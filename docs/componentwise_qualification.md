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

## Qualification status, 2026-10-16

The committed source contains the row/correction certificate and sleeping-worker
implementation described above. Later verdict-witness, mapped-worker,
caller-order snapshot, and packed-vector experiments remain isolated under
`build/`; committing this source does not integrate or promote those candidates.
Production componentwise verification remains enabled. The complete
application-managed comparison added 24 qualification failures, so its
default-off proposal was rejected.

The caller-order snapshot integration remains held after two added Hamrle3
timeouts in its complete generality campaign. The isolated verdict-witness
candidate completed all 1,760 observations and an independent raw-log review
with zero added failures. All 80 shared invalid observations remain retained;
its overall generality performance aggregate is null. Evidence is in
`build/componentwise-verdict-witness-completed-general-decision-20261016/decision.json`.

The fresh 462-observation verdict-witness peer campaign is running in session
67958. Qualification and performance conclusions must await its observed exit,
controller absence, and complete independent review. Its controller is
`build/componentwise-verdict-witness-restored-peer-spice-preparation-20261016/run_campaign.py`.
The latest completed peer evidence does not establish the fastest-in-SPICE goal.
No incomplete or favorable subset supplies an overall competitor comparison.

CKTSO's installed license file was empty. Installing the packaged 719-byte
license restored successful API calls and full-circuit qualification without
changing the solver library, key contents, or clock. This was not evidence of
license expiration. The installation review is
`build/cktso-packaged-license-installation-20261015/review.json`.
