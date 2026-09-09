# Remove factor-preparation phase timing diagnostics

Preceding baseline: `f29a003`; cumulative baseline: `372f9d1`.

Remove 38 lines: the factor/forest phase trace macros and calls, their
local diagnostic timestamps and trace flag, and the factor-exit log.
The live forest preparation timer, numerical diagnostics, preparation
order, allocation and thread joins remain intact. Existing alignment
hints are unchanged.

Release and ASan CTest pass 6/6 (3.66 and 15.44 seconds). All 579 focused
launches pass residual validation at 1e-8. Incremental paired-median
lifecycle changes range from -1.849% to +2.059%; cumulative changes range
from -2.929% to +1.418%. No case reaches the positive 2.5% review threshold.
TSOPF_FS incremental results are +0.161%, -0.311%, and -0.077% on CPUs
8–15, CPUs 0–7, and CPU 8. TSOPF_RS measures +2.059% (12 triplets);
LeGresley_2508 measures -0.631%. These results do not establish zero loss
or cross-machine equivalence; no new full corpus screen was run.

Frozen binaries rotate through one hash-verified execution path with
pinned H100 entrywise-changed refactors and no concurrent builds or source
edits. The final audit verifies revision, exact pending diff, build
provenance, runner/binary/matrix hashes, expected jobs/counts, complete
commands, successful exits, residuals and recomputed medians.

Artifacts: `build/prep-trim-repair-dBcoha/factor-trace-trim-*`; runner:
`validate-factor-trace-trim.py`; preceding frozen executable:
`entry-trace-accepted-kls_bench`.

The broader low-risk cleanup goal remains active.
