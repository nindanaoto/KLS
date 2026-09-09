# Remove output-only lean-preselection diagnostics

Preceding baseline: `3849703`; cumulative baseline: `372f9d1`.

Remove scaled and packed lean-preselection print blocks (17 source lines).
Work estimates, projected savings, construction costs, and preselection
decisions remain unchanged. No test or paper runner consumes this output.
Existing alignment hints remain unchanged; none is added.

Release and ASan CTest pass 6/6 (3.65 and 15.37 seconds). All 579 focused
launches pass residual validation at 1e-8. Incremental paired-median
lifecycle changes range from -1.949% to +1.013%; cumulative changes range
from -1.781% to +1.186%. No case reaches the positive 2.5% review threshold.
TSOPF_FS incremental results are +1.013%, -0.535%, and +0.048% on CPUs
8–15, CPUs 0–7, and CPU 8. LeGresley_2508 measures -1.949%; forced-first
1138_bus measures -0.859%. These results do not establish zero loss or
cross-machine equivalence; no new full corpus screen was run.

Frozen binaries rotate through one hash-verified execution path with
pinned H100 entrywise-changed refactors and no concurrent builds or source
edits. The final audit verifies revision, exact pending diff, build
provenance, runner/binary/matrix hashes, expected jobs/counts, complete
commands, successful exits, residuals and recomputed medians.

Artifacts: `build/prep-trim-repair-dBcoha/preselection-trace-trim-*`;
runner: `validate-preselection-trace-trim.py`; preceding frozen executable:
`direct-trace-accepted-kls_bench`.

The broader low-risk cleanup goal remains active.
