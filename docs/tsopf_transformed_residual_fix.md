# TSOPF recovery through the existing transformed residual worker

This repair follows `cc8b1e5` and retains all preparation cleanup and SNB
SIMD changes. No specialization, threshold, numeric executor, or public
API is introduced.

## Finding and repair

Changed-system solves, not initial-solve timings, exposed a significant
part of the remaining TSOPF loss. The existing `kls_csr32_residual_worker_run`
already maps internal rows through `solve_refine_rinv`, maps columns through
`user_col_perm`, and undoes row/column scaling. Its caller nevertheless
required all four transformation pointers to be null.

The existing dispatch has two policies: ordinary self-checks use an adaptive
CSC-versus-CSR trial; other admitted refinement workloads use the existing
work/horizon eligibility directly. Both were blocked for transformed inputs.
For ordinary self-checks this can mean building CSR metadata and executing
CSC on both trial arms. TSOPF uses the second policy: a debugger check at its
CSR worker entry shows tolerance 1e-4, row self-check off, and zero adaptive
trial samples. Its recovery does not depend on a timing-trial decision.

Remove those two restriction lines. Keep the ordinary eligibility checks,
the measured profitability decision where applicable, scalar
fallback on unavailable resources, transpose/multiple-RHS restrictions,
and all residual acceptance/correction rules. A comment documents why the
transformed case belongs in the trial. No numeric kernel changes are needed.

This restriction also existed before the cleanup. It is a demonstrated
dispatch mismatch and a general way to recover the lost performance, **not
proof of the low-level reason the cleanup originally slowed TSOPF**. That
earlier code/allocation-layout sensitivity remains a separate observation.

## Correctness coverage

The standalone `kls_transformed_residual` test uses row-permuted cyclic-band
matrices with independent row/column magnitude variation. It exercises
changed matrix values and RHS vectors, zero solution entries, scaling modes
0 and 2, ordinary solves, transpose solves, and multiple RHS. Sizes 4096 and
65536 cross the existing CSR column-index width boundary. Existing test-only
environment controls force matching, refinement, and parallel eligibility;
they are not used for reported benchmark timings. An additional SPRAL
matching case exercises the optional scaling implementation when available.

Debugger checks confirm actual entry into the CSR worker with non-null row
permutations for both index widths. These paused runs are coverage evidence,
not performance measurements. The optional SPRAL fixture can fall back when
its factor has no eligible retained crew; no claim of forced parallel scaled
execution is made for that fixture.

Release and ASan/UBSan CTest both pass 6/6, including the new test and the
existing smoke suite. The worker arithmetic and workspace layout are unchanged.

## Performance protocol

Artifacts: `build/prep-trim-repair-dBcoha/transformed-*`.

The initial screen (`transformed-csr-*`) has 168 valid launches against
pre-cleanup `254ea42` and cleanup `da66a15`. The final comparison uses the
same pre-cleanup binary, the committed `cc8b1e5` binary, and the new patch.
It rotates execution order and runs comparisons sequentially, without
overlapping builds or tests. H100 includes 99 entrywise-changing refactors
and 100 solves; every changed-system residual is verified at 1e-8.
Metadata retains binary/matrix hashes, source diff, build provenance, and
commands. Automatic policies are used without diagnostic overrides.

The final campaign contains 444 launches across 15 matrices and 35
configurations: 24 triplets per eight-thread TSOPF CPU set, eight at one
thread, three per eight-thread control configuration, and two per selected
large serial control. This is a focused regression campaign, not a full
paper or Xyce benchmark rerun. Percentages are medians of paired lifecycle
changes; negative means faster.

## Final results

All 444/444 final comparison launches are valid.

| TSOPF configuration | vs pre-cleanup | vs latest commit | Changed-system solve, latest → fix |
|---|---:|---:|---:|
| 8 threads / CPUs 8-15 | -6.35% | -8.11% | 166.1 → 97.2 μs |
| 8 threads / CPUs 0-7 | -7.16% | -9.06% | 163.7 → 97.2 μs |
| 1 thread / CPU 8 | 0.19% | 0.19% | 242.4 → 243.6 μs |

The eight-thread lifecycle regression is recovered on both tested CPU
sets. Refactor timings remain similar; faster changed-system residuals
account for the improvement. One-thread timing is essentially unchanged.

The retained serial SNB gains remain intact: transient, onetone1, twotone,
and ASIC_680ks are all within +0.40% of the latest committed version in the
two-triplet screen and remain faster than pre-cleanup. Those low-repetition
control results are not statistical guarantees for every workload.

Three control warnings were rechecked with 12 triplets each. All 108
additional launches pass; the original slowdowns do not reproduce:

| Control | CPUs | vs latest commit |
|---|---|---:|
| TSOPF_RS_b9_c6 | 8-15 | -1.02% |
| rajat03 | 0-7 | -0.10% |
| bips98_1142 | 0-7 | -0.97% |

### Same-executable ablation

Alternate the final executable's default dispatch against the existing
`KLS_DISABLE_GENERIC_PARALLEL_CONTRACT_RESIDUAL=1` control, using 12 pairs
per CPU set. All 48 launches pass. Default dispatch is 8.69% faster on CPUs
8-15 and 6.74% faster on CPUs 0-7. This comparison holds executable layout
constant and supports attributing the recovery to the newly reachable
parallel residual path, rather than incidental compilation layout changes.

Together the final campaign, control rechecks, and ablation contain
**600/600 valid launches** (including baseline executables), in addition to
the 168-launch initial screen. The trace/debugger checks are separate from
those timing samples. No full paper or Xyce rerun was performed.

Final executable SHA-256:
`a386eaa01c3b777c0ddcf713e4b740a913c2720bb2d4c9f8d68c02fc686da3b0`.
