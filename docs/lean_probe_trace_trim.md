# Remove output-only lean-probe diagnostics

Preceding baseline: `4cb935a`; cumulative baseline: `372f9d1`.

Remove ten `KLS_TRACE_LEAN_PROBE` logging blocks (83 source lines) from
the main driver, row refactor and egraph refactor code. The reported
timing samples, lifecycle audits, choices, admission checks and numerical
work remain unchanged. No test or paper runner consumes this output.
Existing alignment hints are retained; none is added for this chunk.

Release and ASan CTest pass 6/6 (3.68 and 15.44 seconds). All 579 focused
launches pass residual validation at 1e-8. Incremental paired-median
lifecycle changes range from -1.270% to +4.863%; cumulative changes range
from -2.395% to +0.218%. The only 2.5% review flag is forced-first-factor
1138_bus on CPUs 8–15 (+4.863%, eight pairs).

Two independent 24-pair confirmations measure +1.649% and -0.373%
incremental, and -1.527% and -0.585% cumulative. All 144 confirmation
launches are valid. The pooled incremental median across those 48 pairs
is +0.660%. The initial flag therefore is not repeatable at the review
threshold. These results do not establish zero loss or cross-machine
equivalence; no new full corpus screen was run.

All 723 launches rotate frozen binaries through one hash-verified
execution path, with pinned H100 entrywise-changed refactors and no
concurrent builds or source edits. Audits verify revision, exact pending
diff, build provenance, runner/binary/matrix hashes, expected jobs/counts,
complete commands, successful exits, residuals and recomputed medians.

Artifacts: `build/prep-trim-repair-dBcoha/lean-probe-trim-*`,
`lean-probe-repeat*`, and `lean-probe-confirm*`; runner:
`validate-lean-probe-trim.py`; preceding frozen executable:
`row-model-accepted-kls_bench`.

The broader low-risk cleanup goal remains active.
