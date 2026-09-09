# Remove factor-entry diagnostics

Preceding baseline: `c0f0337`; cumulative baseline: `372f9d1`.

Remove 64 lines of output-only factor-entry instrumentation: the local
trace flag, diagnostic clock and phase macro/calls, predicted-factor
zero-diagonal census, KLU-first metadata log and scale-adoption log.
Numerical diagnostics, solver decisions and live elapsed-time accounting
remain intact. Existing alignment hints are unchanged.

Release and ASan CTest pass 6/6 (3.65 and 15.33 seconds). All 579 focused
launches pass residual validation at 1e-8. Incremental paired-median
lifecycle changes range from -3.070% to +1.174%; cumulative changes range
from -4.756% to +1.817%. No case reaches the positive 2.5% review threshold.
TSOPF_FS incremental changes are +1.009%, -0.281%, and -0.013% on CPUs
8–15, CPUs 0–7, and CPU 8. LeGresley_2508 measures +1.174% incrementally
and -4.756% cumulatively. These results do not establish zero loss or
cross-machine equivalence; no new full corpus screen was run.

Frozen binaries rotate through one hash-verified execution path with
pinned H100 entrywise-changed refactors and no concurrent builds or source
edits. The final audit verifies revision, exact pending diff, build
provenance, runner/binary/matrix hashes, expected jobs/counts, complete
commands, successful exits, residuals and recomputed medians.

Artifacts: `build/prep-trim-repair-dBcoha/entry-trace-trim-*`; runner:
`validate-entry-trace-trim.py`; preceding frozen executable:
`consult-trace-accepted-kls_bench`.

The broader low-risk cleanup goal remains active.
