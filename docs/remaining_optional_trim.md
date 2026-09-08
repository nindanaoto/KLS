# Remaining optional-path removal

Relative to `d0ae5a6`, this batch removes all seven candidate areas from the
preceding review:

- Opt-in EGraph supernode-update modes, cached-only execution, subtree splits,
  exclusive producer panel caches, and their internal atomic telemetry.
- Opt-in EGraph ready-queue construction, execution, state, and counters.
- Optional row L/U 32-bit index and target-position mirrors, their parallel
  builder, and exclusive numeric/solve consumers.
- Prestatic reuse of raced METIS symbolic analysis, including the now-unused
  symbolic-only race state.
- Static-match AUTO-scale override.
- Alternate BTF matching trial in bundled KLU: second analysis, matching
  fingerprints, second-factor arbitration, and alternate-symbolic cleanup.
- Constant predicates, an empty preprocessor block, the dead restored-BTF
  state path, the now-unreachable general BTF update branch, and a duplicated
  row scalar-update condition exposed by these removals.

Normal cached-supernode batching, first-factor panels and their shared
acceptance helper, normal BTF matching, the METIS race, compact row-input maps,
row/separator schedulers, double-precision refinement/recovery, and worker
affinity remain. In particular, the input-position mirror has a separate
production constructor; that constructor and its consumers are retained.

Public KLS statistics and bundled-KLU reserved fields retain their layouts.
Retired executor statistics are explicitly zero. The old supernode smoke
fixture now checks accurate EGraph execution with the retired switch set and
zero retired statistics. Historical reports remain untouched.

The batch removes **3,387 net lines from `src/`**, plus **314 net lines from
bundled KLU**. `kls_egraph_refactor.inc` shrinks from 11,398 to 9,651 lines.

## Validation

Final Release and ASan/UBSan CTest each pass **5/5**. Edited solver and KLU
sources compile without warnings. Recompiling the unchanged benchmark source
emits its existing pedantic function-pointer and long-string warnings.
Removed-switch searches in production source and `git diff --check` are clean.
A further Release smoke run on CPUs 8-15 with all seven retired switches
set to 1 also passes.

The paired matrix-lifecycle trial is recorded in
`build/remaining-optional-trim-ptCY1I/`. It uses frozen before/after binaries,
alternating order, one process at a time, H100 entrywise updates of amplitude
0.001, and verification of every refactor at residual limit `1e-8`.
TSOPF_FS_b9_c6 receives 20 one-thread pairs and 40 pairs on each eight-thread
cache domain. Fifteen controls receive four pairs per domain: 440 launches
across 16 matrices. Memory and per-launch time limits are 80GiB and 90 seconds.

The standard provenance checker correctly rejects dirty dependency sources.
Because bundled KLU is deliberately part of this change, a run-local copy
permits exactly the five reviewed KLU files and the exact dependency patch
SHA256 `acaa54cc5f3fcdb021e32f2309055513c36bcf513658307441e3ea7eca111f9c`.
It also requires rebuilding after those files changed. Compiler identities,
all linked-target flags (including SPRAL), submodule revisions, external
objects, and runtime libraries must still match. The global checker is
unchanged. The checker, runner, source patch, and binary hashes are recorded
with the raw results.

These bounded matrix tests are not a full paper-corpus or Xyce SPICE rerun;
they cannot establish full-corpus performance parity.

## Results and performance caveat

All **440 broad-trial launches passed**, with maximum refactor residual
`5.45371652e-9`. TSOPF_FS_b9_c6 median paired lifecycle changes were +3.80%
at one thread and +0.43% / +0.41% on the two eight-thread domains. Neither
version had an eight-thread TSOPF launch above the earlier 0.6ms average
steady-refactor outlier cutoff. The affinity repair remains effective in
this sample.

Other controls (four pairs per domain; positive means slower):

| Matrix | 32MiB, 8 threads | 96MiB, 8 threads |
| --- | ---: | ---: |
| TSOPF_FS_b9_c1 | +1.06% | +0.31% |
| TSOPF_RS_b9_c6 | -0.29% | +0.00% |
| onetone1 | -0.18% | +3.21% |
| twotone | +0.73% | -0.02% |
| rajat25 | +0.22% | +0.85% |
| adder_dcop_01 | -1.93% | -3.17% |
| 1138_bus | +0.97% | +2.49% |
| bips98_1142 | +0.53% | +1.03% |
| ASIC_100ks | -1.18% | +0.81% |
| ASIC_320k | +0.47% | +0.08% |
| ASIC_680ks | +1.90% | +1.99% |
| ckt11752_dc_1 | -0.35% | -0.76% |
| circuit_4 | +1.37% | -1.69% |
| transient | -0.44% | -0.25% |
| rajat20 | -0.22% | -0.99% |

Independent rechecks in `build/remaining-optional-recheck-HfGFgJ/` include
30 one-thread TSOPF pairs, 12 onetone1 pairs on 96MiB, and 20 adder_dcop_01
pairs on 96MiB. All 124 launches passed. Their median paired changes were
**+3.33%, -0.78%, and -0.39%**, respectively. The initial onetone1 slowdown
and large adder improvement did not repeat.

The one-thread TSOPF regression did repeat and is not claimed to be noise.
Its first-refactor median rose from 73.03ms to 80.53ms, while steady-refactor
medians were 0.993ms before and 0.988ms after, and solve medians were 0.463ms
and 0.460ms. Initial factor time stayed about 10.1ms. This localizes the
observed cost to the first-refactor phase; it does not identify the responsible
subroutine or prove a compiler-layout explanation. Further attribution remains
open. The removals are retained as requested, with this known performance
tradeoff, not represented as performance-neutral.

In total **564/564 benchmark launches passed**. The final Release binary
matches both frozen after binaries, SHA256
`e1199c8854e83b3c4694224f91b0dd22e15246c309389417cfa978a7c0f7475e`.
