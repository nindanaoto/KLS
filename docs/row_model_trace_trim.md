# Remove output-only row-model gate diagnostics

Preceding baseline: `ad11d75`; cumulative baseline: `372f9d1`.

Remove the 12-line `KLS_TRACE_ROW_MODEL_GATE` print block. The actual
admission conditions, eligibility checks, numerical work and thresholds
remain unchanged. No test or paper runner consumes this diagnostic.
The existing compact solve and SNB alignment hints remain intact.

Release and ASan CTest pass 6/6 (3.67 and 15.33 seconds). All 579 focused
launches pass numerical validation at 1e-8. Incremental paired-median
lifecycle changes range from -1.861% to +2.878%; cumulative changes range
from -2.041% to +0.427%. The sole 2.5% review flag is forced-first-factor
1138_bus on CPUs 8–15 (+2.878%, eight pairs). An independent 24-pair
confirmation measures -0.494% incremental and -1.904% cumulative, with
all 72 launches valid. Thus the flag is not repeatable; no new alignment
workaround is warranted. These results do not establish zero loss or
cross-machine equivalence, and no new full corpus screen was run.

All 651 launches use rotating frozen binaries through one hash-verified
execution path, pinned H100 entrywise-changed refactors, and no concurrent
builds or source edits. Audits check revision, pending diff, binary/runner/
matrix hashes, expected jobs, complete commands, successful exits,
residuals, and medians recomputed from raw records. The focused audit also
verifies build provenance against the established build configuration.

Artifacts: `build/prep-trim-repair-dBcoha/row-model-trim-*` and
`row-model-repeat*`; runner: `validate-row-model-trim.py`; preceding
frozen executable: `i32-aligned-accepted-kls_bench`.

The broader low-risk cleanup goal remains active.
