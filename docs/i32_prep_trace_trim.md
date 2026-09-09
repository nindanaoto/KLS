# Remove diagnostic-only compact-index preparation timing

Preceding baseline: `dee85e2`; cumulative baseline: `372f9d1`.

Remove `KLS_TRACE_I32_PREP`, its four diagnostic timestamp locals,
conditional clock samples, and print block (25 source lines). Index
construction, copying, cache preparation, allocation/error handling and
solver selection remain unchanged. No test or paper runner consumes the
removed output. Existing alignment hints remain unchanged.

Release and ASan CTest pass 6/6 (3.59 and 15.30 seconds). All 579 focused
launches pass residual validation at 1e-8. Incremental paired-median
lifecycle changes range from -1.013% to +2.253%; cumulative changes range
from -1.266% to +1.277%. No case reaches the positive 2.5% review threshold.
The largest incremental increase is twotone on CPUs 8–15, with three
pairs; this was not independently repeated. TSOPF_FS incremental results
are +0.430%, +0.052%, and -0.003% on CPUs 8–15, CPUs 0–7, and CPU 8.
LeGresley_2508 measures -0.431%, and forced-first 1138_bus +1.186%.
These results do not establish zero loss or cross-machine equivalence;
no new full corpus screen was run.

Frozen binaries rotate through one hash-verified execution path with
pinned H100 entrywise-changed refactors and no concurrent builds or source
edits. The final audit verifies revision, exact pending diff, build
provenance, runner/binary/matrix hashes, expected jobs/counts, complete
commands, successful exits, residuals and recomputed medians.

Artifacts: `build/prep-trim-repair-dBcoha/i32-prep-trim-*`; runner:
`validate-i32-prep-trim.py`; preceding frozen executable:
`lean-probe-accepted-kls_bench`.

The broader low-risk cleanup goal remains active.
