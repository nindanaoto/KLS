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
