# SNB scratch alignment after dead-state cleanup

The cleanup was committed as `7abe110` before this repair. All removals in
that commit remain, including the retired diagnostic paths. The SNB union
`qsort`, producer binary search, and worker-affinity fixes also remain.

## Investigation

The recorded regression is the one-thread `transient` H100 lifecycle slowdown
in [the cleanup report](dead_state_and_snb_sort_trim.md): sixteen paired runs
measured +6.60%, concentrated in steady SNB refactors (about 29.2ms to 33.3ms).
Both versions used SNB; the slowdown was already present before replacing the
union insertion sort with `qsort`.

Disassembly of the sampled dense-update region shows identical instructions
at different addresses in the pre-cleanup and cleanup binaries. Debugger
observations also show different scratch-buffer addresses after the private
state shrank, while the dominant panel's alignment was unchanged. Thus source
arithmetic did not acquire additional work, but executable and allocation
layout both changed. The exact hardware-level mechanism is not established.

Three rejected experiments are retained in `build/transient-repair-t9O9nk/`:

- Adding scratch-pointer `restrict` qualifiers did not recover performance.
- Changing the dense product to column-major scratch did not recover it either.
- Aligning both panels and scratch recovered `transient`, but slowed TSOPF
  steady refactors. Panel allocation was therefore restored.

Only the scratch-alignment change is retained. A debugger-time allocation
experiment was also attempted, but debugger pauses affect measured timings
and can affect the timed engine trial. It is not used as performance evidence
or proof of a particular hardware cause.

## Repair

Allocate the per-worker `ub_w` and `gemm_w` buffers through the existing
`kls_aligned_double_values` helper instead of `malloc`. This gives these
frequently reused dense-update buffers the same 64-byte alignment guarantee
as other KLS numeric workspaces. No buffer dimensions, arithmetic, scheduling,
engine eligibility, or measured engine-selection rules change.

The existing helper checks the byte count and allocation result. Partial
allocation failures still use the existing SNB cleanup/fallback path, and
the allocations remain compatible with `free`. No dead-field padding,
matrix-name check, size-band exception, or new tuning threshold is introduced.

The deferred-sort smoke fixture now runs both narrow and dense scratch-buffer
updates on its dense and interleaved even/odd matrices. Test-only narrow-update
settings are saved and restored. Each variant checks changed-value normal and
transpose solves, in addition to the solve before deferred sorting.

## Validation

Release and ASan/UBSan CTest each pass **5/5**, with warning-free builds.
A separate traced smoke run confirms successful SNB trials for all four
4,096-supernode fixture variants.

Performance artifacts are in `build/transient-repair-t9O9nk/`. The final
runner compares frozen `e01e3ca` (original), `7abe110` (slim), and repair
binaries in rotating order, sequentially, without tracing. Compiler options,
dependency revisions, linked objects, and runtime libraries match. Each
launch uses H100 entrywise updates of amplitude 0.001, verifies every refactor
at residual tolerance `1e-8`, and has a 90-second timeout and 80GiB address-space
limit. CPU sets are 8 (one thread), 8–15 and 0–7 (eight threads).

The bounded check covers `transient`, TSOPF_FS_b9_c6, fourteen controls, and
Raj1/mac_econ_fwd500 checks on both eight-thread cache domains. It is not a
full paper-corpus or Xyce rerun.

Twelve rotating-order triples for one-thread `transient` give these medians:

| Version | H100 lifecycle (s) | Steady refactor (ms) |
| --- | ---: | ---: |
| Before cleanup (`e01e3ca`) | 5.00003 | 29.147 |
| Committed cleanup (`7abe110`) | 5.32920 | 33.223 |
| Scratch-alignment repair | 4.95815 | 29.343 |

Median paired lifecycle changes are **-7.33% versus the cleanup** and
**-1.07% versus the original**. These are paired statistics, not ratios of
the marginal medians in the table. Steady refactors recover most of the
lost time and are within 0.20ms of the original; first-refactor time retains
the cleanup's approximately 70ms saving. All three versions use SNB.

TSOPF_FS_b9_c6 gets twelve triples at one thread and six on each eight-thread
domain. The repair changes lifecycle time versus the cleanup by **+0.11%**,
**+0.08%**, and **-0.47%**, respectively. Relative to the pre-cleanup version,
the retained improvements are **4.75%**, **1.72%**, and **1.76%**. One-thread
first refactor remains approximately 16.67ms.

All **390/390** final-sweep launches pass. The fourteen controls receive two
triples at one thread and on each eight-thread domain; eight-thread `transient`
receives three triples per domain. Raj1 and mac_econ_fwd500 each receive one
triple per eight-thread domain, so their checks are accuracy and gross-regression
checks rather than precise performance estimates.

The small initial samples flagged four comparisons above 2% against either
baseline. Sixteen fresh triples per flagged configuration pass **192/192**:

| Rechecked configuration | Lifecycle change vs cleanup | vs original |
| --- | ---: | ---: |
| TSOPF_RS_b9_c6, eight threads, 32MiB domain | -1.28% | -1.74% |
| TSOPF_RS_b9_c6, eight threads, 96MiB domain | +0.36% | +0.42% |
| adder_dcop_01, eight threads, 96MiB domain | -0.78% | +1.13% |
| 1138_bus, eight threads, 96MiB domain | +0.18% | +0.29% |

Thus none of the flagged slowdowns repeats above 2%. Other final-sweep
comparisons against the cleanup range from -20.38% to +1.60%. In particular,
one-thread onetone1 and twotone improve by approximately 14.5%, and ASIC_100ks,
ASIC_320k, and ASIC_680ks also improve. These control gains use smaller samples
than the primary regression case and should not be treated as universal gains.

The final sweep and rechecks total **582/582 valid launches**, with maximum
refactor residual **5.45371652e-9**. No material new regression is established
in this bounded check; this does not establish performance neutrality on
unmeasured machines or the entire paper corpus.

The final Release binary matches both frozen repair binaries, SHA256
`9dabab307f0647fb09d7213c695c3c9e62ba4284c65d8c313a5f9d5b3259ef15`.
`final-metadata.json`, `final.jsonl`, `final-summary.json` and the corresponding
`recheck-*`/`recheck.jsonl` files preserve provenance, commands, observations,
and paired summaries. Rejected exploratory trials are recorded separately
and are not counted in the 582 final validation launches.
