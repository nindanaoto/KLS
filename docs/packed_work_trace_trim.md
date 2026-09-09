# Remove output-only packed-work diagnostics

Preceding baseline: `4b8e7f8`; cumulative baseline: `372f9d1`.

Remove the compact two-block trace helper, schedule output, and local
packed-work census (51 source lines). Live scheduling counts, descriptors,
reciprocal publication, allocation/error handling and numerical work remain
unchanged. No test or paper runner consumes this output. Existing alignment
hints remain unchanged; none is added.

Release and ASan CTest pass 6/6 (3.68 and 15.39 seconds). All 579 focused
launches pass residual validation at 1e-8. Incremental paired-median
lifecycle changes range from -1.181% to +1.611%; cumulative changes range
from -1.451% to +0.762%. No case reaches the positive 2.5% review threshold.
TSOPF_FS incremental results are -0.341%, -0.785%, and -0.048% on CPUs
8–15, CPUs 0–7, and CPU 8. LeGresley_2508 measures -0.842%; forced-first
1138_bus measures +1.028%. These results do not establish zero loss or
cross-machine equivalence; no new full corpus screen was run.

Frozen binaries rotate through one hash-verified execution path with
pinned H100 entrywise-changed refactors and no concurrent builds or source
edits. The final audit verifies revision, exact pending diff, build
provenance, runner/binary/matrix hashes, expected jobs/counts, complete
commands, successful exits, residuals and recomputed medians. An earlier
audit attempt correctly rejected the still-incomplete last control; the
final audit was run after normal process completion.

Artifacts: `build/prep-trim-repair-dBcoha/packed-trace-trim-*`; runner:
`validate-packed-trace-trim.py`; preceding frozen executable:
`affinity-trace-accepted-kls_bench`.

The broader low-risk cleanup goal remains active.
