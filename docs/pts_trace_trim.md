# Remove output-only PTS diagnostics

Preceding baseline: `1ff9ffc`; resumed-pass baseline: `372f9d1`.

Remove the private `KLS_TRACE_PTS` build, refactor, and solve messages;
diagnostic phase clocks and marks; rejection-stage labels; and the private
cut-scoring helper's trace argument. No trace references remain in source.
This removes 125 net C-source lines. Partition scores, validity certificates,
allocation, numerical work, and public state remain intact. In particular,
solve/refactor acceptance samples and two-stage decision timers are retained.
The variable no longer produces PTS debug output; no paper runner or test
consumes that output.

Release and ASan/UBSan builds succeed; CTest passes 6/6 in each (3.75 and
15.57 seconds). Stripped executable files differ; comparison copies are
under `build/pts-trace-fFLZ7l/`. The SNB kernel remains 64-byte aligned at
0x641c0, size 0x1931; no placement workaround was added.

All 651 launches are valid: 36 onetone, 495 AUTO, 48 forced-first controls,
and 72 confirmation launches. Three versions rotate through one hash-verified
execution path, with pinned H100 entrywise-changed refactors and residual
validation at 1e-8. No builds or edits overlap timing.

The forced-first 1138_bus control initially flags at +3.315% incrementally
and +1.214% cumulatively over eight repetitions. An independent 24-repetition
confirmation measures -0.078% / -1.480%, with no failures. The initial
slowdown does not persist. No other initial control reaches the positive
2.5% review trigger. TSOPF_FS_b9_c6 incremental changes are +0.613%, +0.559%,
and -0.725% on CPUs 8–15, CPUs 0–7, and CPU 8; cumulative changes are
-0.485%, -0.417%, and -0.398%. These measurements do not prove zero loss
or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/pts-trace-trim-*` and
`pts-trace-repeat*`; runner: `validate-pts-trace-trim.py`; preceding frozen
executable: `snode-trace-accepted-kls_bench`. The final audit verifies all
expected jobs/counts, exact pending source diff, build provenance,
runner/binary/matrix hashes, complete commands, successful exits, residual
validity, and all summary medians recomputed from raw records.

This is the third source chunk since the row-helper cumulative checkpoint.
The broader low-risk cleanup goal remains active.
