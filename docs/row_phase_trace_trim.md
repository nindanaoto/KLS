# Remove output-only row phase diagnostics

Source commit: `161a14d`; preceding baseline: `a201c01`;
resumed-pass baseline: `372f9d1`.

Remove `KLS_TRACE_LEAN_PARALLEL_SETUP`, `KLS_TRACE_REFACTOR_PHASES`,
and `KLS_TRACE_PARALLEL_I32_PERMUTE` output and diagnostic-only clocks.
This removes 78 net C-source lines. Packed-kernel comparison timestamps
remain, renamed to `packed_sample_begin` and `packed_sample_end`; their
sampling predicate, sample storage, and kernel-selection decisions remain.
No paper runner or test consumes the removed output.

Release and ASan builds succeed, with CTest passing 6/6 in each (3.69 and
15.38 seconds). The SNB kernel remains 64-byte aligned at 0x641c0, size
0x1931. No new alignment workaround or numerical policy change was added.

All 651 launches are valid: 36 onetone controls, 495 AUTO controls,
48 forced-first controls, and 72 confirmation launches. Three frozen
versions rotate through one hash-verified execution path with pinned CPUs,
H100 entrywise-changed refactors, and residual validation at 1e-8.
No source edits or builds overlap timing.

LeGresley_2508 on CPUs 8–15 initially flags at +3.182% incrementally and
+1.740% cumulatively over 12 repetitions. Independent confirmation over
24 repetitions measures -0.918% / -0.102%, with no failures. The initial
slowdown does not persist. No other focused case reaches the positive
2.5% review threshold. TSOPF_FS_b9_c6 incremental changes are -0.117%,
-0.037%, and -0.451% on CPUs 8–15, CPUs 0–7, and CPU 8; cumulative
changes are -0.663%, -1.012%, and -0.079%, respectively. These results
do not prove zero performance loss or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/row-phase-trim-*` and
`row-phase-repeat*`; runner: `validate-row-phase-trim.py`; preceding frozen
binary: `pts-trace-accepted-kls_bench`. The final audit verifies expected
jobs/counts, committed revision and empty source diff, build provenance,
runner/binary/matrix hashes, complete commands, successful exits, residual
validity, and all reported medians recomputed from raw records.

This is the fourth source cleanup since the row-helper cumulative corpus
checkpoint (292 net C-source lines across these four chunks). A new
cumulative corpus screen is required before further source deletion.
The broader low-risk cleanup goal remains active.
