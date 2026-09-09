# Remove output-only direct-path diagnostics

Preceding baseline: `eb6e258`; cumulative baseline: `372f9d1`.

Remove seven direct-refactor/solve diagnostic blocks, including the
public-direct gate's diagnostic-only thread-local print counter
(55 source lines). Live gates, measured tournaments, timing samples,
state transitions and numerical work remain unchanged. No test or paper
runner consumes this output. Existing alignment hints remain unchanged;
none is added.

Release and ASan CTest pass 6/6 (3.59 and 15.30 seconds). All 579 focused
launches pass residual validation at 1e-8. Incremental paired-median
lifecycle changes range from -1.049% to +0.769%; cumulative changes range
from -0.658% to +4.599%. The only positive 2.5% review flag is TSOPF_RS
on CPUs 8–15 versus the cumulative baseline (+4.599%, twelve pairs),
while its incremental result is -0.343%.

An independent 24-pair confirmation measures -0.490% incremental and
+0.784% cumulative, with all 72 launches valid. The cumulative flag
does not repeat at the review threshold. These results do not establish
zero loss or cross-machine equivalence; no new full corpus screen was run.

All 651 launches rotate frozen binaries through one hash-verified
execution path with pinned H100 entrywise-changed refactors and no
concurrent builds or source edits. Audits verify revision, exact pending
diff, build provenance, runner/binary/matrix hashes, expected jobs/counts,
complete commands, successful exits, residuals and recomputed medians.

Artifacts: `build/prep-trim-repair-dBcoha/direct-trace-trim-*` and
`direct-trace-repeat*`; runner: `validate-direct-trace-trim.py`; preceding
frozen executable: `contract-trace-accepted-kls_bench`.

The broader low-risk cleanup goal remains active.
