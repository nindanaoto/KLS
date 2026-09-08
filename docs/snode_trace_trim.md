# Remove remaining supernode preparation tracing

Preceding baseline: `e0bc72f`; resumed-pass baseline: `372f9d1`.

Remove all four output-only `KLS_TRACE_SNODE` preparation messages and its
README mention. No references remain in current source or README. The
coverage count still controls run acceptance, and the panel-offset count
still constructs retained storage; neither is removed. Sorting, eligibility,
allocation, fallback, arithmetic, and elapsed-time accounting remain intact.
This removes 27 C-source lines. The variable no longer prints preparation
diagnostics; no paper runner or test consumes those messages.

Release and ASan/UBSan builds succeed; CTest passes 6/6 in each (3.82 and
15.60 seconds). Stripped executables differ; comparison copies are under
`build/snode-trace-WZsS6c/`. The SNB kernel moves from 0x64540 to 0x64500,
remaining 64-byte aligned with size 0x1931. No placement padding or
matrix-specific workaround was added.

All 579 focused launches are valid: 36 onetone, 495 AUTO, and 48 forced-first.
Three versions rotate through one hash-verified execution path with pinned
H100 entrywise-changed refactors and residual validation at 1e-8. No builds
or source edits overlap timing.

Paired lifecycle median changes range from -5.317% to +0.981% incrementally
and from -2.277% to +1.213% against the resumed-pass baseline. No positive
change reaches the 2.5% review trigger. TSOPF_FS_b9_c6 incremental changes
are -0.417%, +0.867%, and +0.041% on CPUs 8–15, CPUs 0–7, and CPU 8;
cumulative changes are -0.266%, +0.718%, and -0.022%. These measurements
do not prove zero regression or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/snode-trace-trim-*`;
runner: `validate-snode-trace-trim.py`; preceding frozen executable:
`snb-trace-accepted-kls_bench`. The final audit verifies expected jobs/counts,
exact pending source/README diff, build provenance, runner/binary/matrix
hashes, complete commands, successful exits, residual validity, and all
summary medians against raw records.

This is the second source chunk since the completed row-helper cumulative
corpus checkpoint. The broader cleanup goal remains active.
