# Batched diagnostic cleanup

Baseline: `f23a157`. Remove 599 source lines across the vendored KLU kernel,
factor wrapper, condition diagnostics and 32-bit symbol wrapper. This batch
includes pipeline phase and wall-time profiling, row and level timing,
dense-tail timing, pivot/row logs, file dumps, fingerprints, their counters
and TLS state, the statistics-only supernode/fill probe, discarded row-kernel
min/max calculations, and two singular/NaN debug-print hooks.

Keep solver selection, thread and memory controls, compensated arithmetic,
deterministic ordering, pivot handling, singular status/rank reporting,
abort checks, allocation cleanup, and public statistics/ABI unchanged.
The two unconditional invalid-pivot error messages remain with their abort
checks. CPU-topology reads and live contract-statistics controls are not
diagnostic output and remain. No KLU optional PROF/TRACE/DUMP/STATS getenv
hooks remain in the inspected source tree.

## Correctness and provenance

Final release CTest passes 6/6 (4.07 s), and final ASan CTest passes 6/6
(17.25 s), including the existing singular-rejection tests. Both builds
are warning-free; `git diff --check` passes.

The 588-line version was tested before the final eleven lines of debug
prints were removed. Its reviewed dependency patch hash is
`ee83e25a4488071423e3b75b17d0b2da357f18dc71f6e3416dda24d8c2cb051c`.
The final 599-line version records
`09d3d50a88760e654850679b36f48249a084aaca7b8b59ddf9f3de92dc84ebbf`.
Each timing phase freezes source and binaries; no builds overlap timing.

## Performance

The 588-line version passes all 579 focused H100 launches against accepted
and cumulative `372f9d1`. Median paired lifecycle changes range from
-1.176% to +3.727% incremental. First-factor bcircuit initially crosses the
2.5% review threshold (+3.727% incremental / +2.918% cumulative).
An independent 171-launch confirmation gives bcircuit -0.756%, 1138_bus
-0.123%, and twotone +1.174% versus accepted; all validate.

After removing the final two debug prints, rerun release/ASan and 387
targeted launches against both accepted and the tested 588-line binary.
The final-build checks are:

| Case / mode | vs accepted | vs 588-line version |
| --- | ---: | ---: |
| bcircuit / first-factor | +0.378% | -0.679% |
| 1138_bus / first-factor | -0.095% | +0.704% |
| circuit_4 / first-factor | -0.553% | +0.107% |
| TSOPF_FS_b9_c6 / CPUs 8–15 | -0.291% | -0.011% |
| TSOPF_FS_b9_c6 / CPUs 0–7 | +0.247% | -0.100% |
| twotone / CPUs 8–15 | -0.575% | +0.163% |

All 1,137 launches validate. Audits recompute medians and verify exact jobs,
pinned commands, mode flags, source diffs, executable/matrix hashes, build
provenance, exits and residuals. The flagged slowdown does not reproduce;
no alignment fix is needed. This is a focused campaign plus final-build
targeted checks, not a full-paper or cross-machine performance guarantee.

Artifacts under `build/prep-trim-repair-dBcoha/`: `broad-diagnostic-trim-*`,
`broad-diagnostic-confirm-*`, and `broad-final-trim-*`. Runners:
`validate-broad-diagnostic-trim.py`, `validate-broad-final-trim.py`.
Main audit: `audit-broad-diagnostic-trim.py`; confirmation/final audits use
`audit-eligibility-repeat.py` with exact corresponding jobs, totals, mode
and reviewed dependency hash from each phase's metadata.
