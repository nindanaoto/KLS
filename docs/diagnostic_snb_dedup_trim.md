# Diagnostic cleanup and shared SNB edge arithmetic

Relative to `af89475`, this batch removes **211 net source lines**, excluding
tests and documentation. It retains the SNB scratch-alignment repair, union
sorting and producer lookup fixes, and worker-affinity behavior.

## Removed

- The unused `kls_snb_state::ub` and `gemm_temp` worker-zero aliases and their
  assignments. Callers already use the per-worker arrays directly.
- Trace-only synchronization wait/publication counts. The actual wait masks,
  acquired producer frontiers, publication flags, and memory ordering remain.
- Supernode width/target histograms and the lean-solve histogram. Width-16
  existence is retained as a boolean because it selects a production kernel;
  coverage, fused-run counts, and target lengths used by eligibility remain.
- `KLS_ROW_FINISH_MARK`, its phase-timing variables/calls, and the cleanup and
  row-level preparation jobs' diagnostic-only `seconds` fields and clocks.
- Unused-for-selection lifecycle projections (`remaining_saving` and
  `column_sample_work`) and their trace fields. The measured lifecycle
  comparison, adoption margins, and re-audit state machine remain.
- Five manual override hooks: `KLS_ENABLE_LEAN_ROW_CONSUME`,
  `KLS_LEAN_ROW_PAIR`, `KLS_ENABLE_PADDED_PANELS`, `KLS_PADDED_Z`, and
  `KLS_PADDED_MIN_WORK`. They are now ignored. The default padded union
  tolerance of 2 and minimum work of 1024 are unchanged.

The retired diagnostic switches are `KLS_TRACE_LEAN_SNODE`,
`KLS_TRACE_LEAN_SOLVE_HIST`, and `KLS_TRACE_ROW_PATTERN_FINISH`. Public
statistics and other retained trace facilities are unchanged. Paired-lean
and padded-panel engines remain automatically reachable through their
existing selection/trial paths; neither engine was removed.

## Shared arithmetic

`kls_snb_update_edge` is an always-inline helper shared by ordinary serial
edge execution and the cooperative executor's small-edge branch. It contains
the former serial narrow update and dense gather/TRSM/product/scatter body.
It takes the existing block, supernode, and edge metadata, the narrow-update
setting, and the two scratch buffers; immutable shape/address setup is local
to the helper rather than passed through fifteen separate arguments.
The accumulation order and row-major dense scratch layout are preserved.
The callers retain their dependency waits and crew barriers. Large cooperative
edges retain their existing row-sliced implementation.

For a zero-height dense product, the cooperative small-edge caller now skips
an unnecessary scratch transpose after the same TRSM and panel writeback.
That scratch has no consumer in this case. No matrix-specific route, tuning
threshold, padding field, or out-of-line hot-loop call is introduced.

## Tests

The existing deferred-sort fixture keeps all four serial variants and adds
two four-thread variants (narrow and dense updates). Their first dense block
is twice as long as the others, with a tall panel at columns 68–71. Once the
shorter blocks finish, that panel occupies a level alone and reaches the
cooperative executor. Its earlier producers have short tails, exercising the
shared small-edge body inside the crew barrier. Later edges also exercise
the existing large-edge implementation. Changed-value normal and transpose
solves remain checked.

Trace verification reports a **0.48 cooperative work share**, 4,096 supernodes,
32 levels, and successful four-thread SNB trials in both new variants. An
earlier fixture attempt reached only parallel, not cooperative, execution;
it was corrected before final validation.

Final Release and ASan/UBSan CTest each pass **5/5**, with warning-free builds.
Smoke tests also pass with all eight retired switches set, including extreme
values for the retired padded-panel parameters. None of the removed diagnostic
messages appears. The SNB edge helper is inlined in the Release binary.

## Performance validation

Artifacts are in `build/diagnostic-snb-trim-6Vq14C/`. The runner freezes the
baseline and candidate and checks matching compiler flags, dependency revisions,
linked objects, and runtime libraries. Runs are sequential and alternate order,
with no tracing during timing: H100, entrywise amplitude 0.001, verification
after every refactor at `1e-8`, 90-second timeout, and 80GiB address-space limit.
CPU sets are 8 for one thread, and 8–15 / 0–7 for the two eight-thread domains.

The initial fifteen-argument helper passed 18/18 probe and 284/284 broad-sweep
launches, but introduced steady SNB slowdowns: one-thread onetone1 +5.54%,
twotone +3.36%, ASIC_320k +2.72%, and ASIC_680ks +3.56% lifecycle time.
It was not accepted as the final implementation. Disassembly showed changed
register allocation despite preserved arithmetic; no precise hardware-level
cause is claimed.

Two non-overlap-annotation variants were tested and rejected. Annotating the
target and scratch pointers recovered onetone1 but slowed transient by 3.16%;
annotating only the target was worse on both. No such annotation remains.
The retained metadata-based interface is smaller and preserves the original
arithmetic without those tradeoffs in the initial checks: onetone1 +1.30%,
transient +0.08%, TSOPF_FS_b9_c6 +0.03%, and twotone +0.08%. All 28 launches
in that preliminary check pass.

The bounded sweep covers 20 matrices, including bcircuit and rajat03 to check
the retained padded/lean policies, plus the previous regression controls and
Raj1/mac_econ_fwd500. It is not a full paper-corpus or Xyce rerun.

## Retained-version results

The retained implementation passes **284/284** sweep launches. Twelve pairs
each cover one-thread transient and TSOPF_FS_b9_c6. TSOPF receives six pairs
per eight-thread domain, transient three, and the sixteen controls two pairs
at each of the three configurations. Raj1 and mac_econ_fwd500 receive one
pair per eight-thread domain, so those are accuracy and gross-regression
checks rather than precise performance estimates.

Median paired H100 lifecycle changes versus `af89475`:

| Case | Change |
| --- | ---: |
| transient, one thread | +0.28% |
| TSOPF_FS_b9_c6, one thread | -0.24% |
| TSOPF_FS_b9_c6, eight threads, 32MiB domain | -0.53% |
| TSOPF_FS_b9_c6, eight threads, 96MiB domain | +1.01% |
| onetone1, one thread | +1.38% |
| ASIC_320k, one thread | +1.39% |
| ASIC_680ks, one thread | +1.28% |

One-thread transient steady refactors measure 29.36ms before and 29.32ms
after; its former steady regression remains recovered. One-thread TSOPF
first refactor remains approximately 16.7ms. The sweep's one-thread twotone
sample is +2.24%, so it receives six fresh pairs. Four small eight-thread
comparisons above 2% each receive sixteen fresh pairs.

The small-case rechecks pass 128/128 launches. TSOPF_RS_b9_c6 measures -0.09%
and -0.07% on the two domains; 1138_bus is -0.03% on the 32MiB domain.
However, **1138_bus on the 96MiB domain repeats at +3.61%**. Its marginal
median H100 lifecycle rises from **1.982ms to 2.040ms** (about 0.058ms in
absolute time), using the mapped engine, not SNB. The percentage is the
median of paired changes, not the ratio of those two marginal medians.
This remains a measured caveat; the cleanup is not claimed to be strictly
performance-neutral. No case-specific workaround is added for it.

The six-pair twotone recheck passes 12/12 launches and measures **+0.90%**;
its earlier +2.24% sample does not repeat at that magnitude. The retained
implementation's sweep and rechecks total **424/424 valid launches**, with
maximum refactor residual **5.45371652e-9**.

The final Release binary matches the frozen retained/recheck binaries,
SHA256 `f177e76b8ce7e0a38fc45557b9f5afcbdfee15001ae7e2bc000aca2e97b963f7`.
The retained version's provenance, commands, source diff, observations, and
paired summaries are in the `accepted-*`/`accepted.jsonl` and
`recheck-*`/`recheck.jsonl` artifacts. The earlier fifteen-argument helper's
sweep uses the `final-*`/`final.jsonl` names; annotation and interface probes
are separate. Those earlier versions are not counted in the 424 retained-
version validation launches.
