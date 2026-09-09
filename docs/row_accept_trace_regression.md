# Slimming stop: row-acceptance diagnostic removal

Accepted source baseline: `12bc310`; resumed-pass baseline: `372f9d1`.

An experimental removal of `KLS_TRACE_ROW_ACCEPT` messages deleted 119
lines across `src/kls.c` and `src/kls_row_refactor.inc`, preserving runtime
decisions and samples. Release and ASan CTest passed 6/6 (3.81 and 15.62
seconds). The SNB kernel remained 64-byte aligned, moving from 0x641c0 to
0x64040 with unchanged size 0x1931. Alignment alone does not establish
performance equivalence.

All 579 initial focused launches were valid. LeGresley_2508 on CPUs 8–15
was the only positive 2.5% review flag: +2.791% versus the preceding
version and +4.016% versus the resumed-pass baseline over 12 repetitions.
An independent 24-repetition confirmation (72 additional launches) measured
+2.549% incrementally and +1.317% cumulatively. The candidate was slower
in 18/24 incremental pairs. Both runs exceed the established 2.5% review
threshold. No further rerun was used to seek a passing result.

In confirmation, median lifecycle times were 6.404 ms for the accepted
version and 6.586 ms for the candidate. Median steady refactor/solve time
rose from 10.908 us to 12.421 us; median analysis and initial factor times
did not increase. This points to recurring numerical execution, not just
initialization, as the affected portion. These are per-version medians;
the headline percentages are medians of paired ratios. The exact mechanism
(for example code placement or runtime selection) has not been established.

The full audit verified all 651 valid launches, expected jobs and counts,
revision, exact experimental source diff, build provenance, frozen runner/
binary/matrix hashes, complete pinned commands, successful exits, residual
validity at 1e-8, and every reported median recomputed from raw records.
H100 entrywise-changed refactors and rotating version order used the same
hash-verified execution path. No builds or source edits overlapped timing.

Artifacts: `build/prep-trim-repair-dBcoha/row-accept-trim-*` and
`row-accept-repeat*`; runner: `validate-row-accept-trim.py`. Each `.diff`
preserves the rejected candidate for reproduction. Frozen preceding binary:
`worker-time-accepted-kls_bench`.

The uncommitted 119-line candidate was reverted exactly to `12bc310`.
Previously validated slimming remains intact. This ends the requested pass
at a repeatable, roughly 3% incremental regression, rather than claiming
that every remaining cleanup candidate is exhausted. Further removal or
investigation of the affected code requires a new pass.
