# Remove output-only SNB tracing

Preceding source baseline: `eaec5d7`; resumed-pass baseline: `372f9d1`.

Remove the private `KLS_TRACE_SNB` predicate and all eight printing sites,
including the diagnostic-only local preparation census. No trace references
remain in the source. Solver-owned state, public statistics, trial timers,
adoption/retrial policies, allocation bounds, and numerical checks remain
unchanged. This removes 62 C-source lines. The environment variable no longer
produces SNB debug output; no paper runner or test consumes that output.

Release and ASan/UBSan builds succeed; CTest passes 6/6 in each (3.74 and
15.61 seconds). The smoke fixture exercises sparse/dense SNB updates and
cooperative trials. Stripped executable files differ; comparison copies
are under `build/snb-trace-QBnazZ/`. The SNB kernel remains 64-byte aligned
at 0x64540, size 0x1931.

All 579 focused launches are valid: 36 onetone, 495 AUTO, and 48 forced-first.
Three versions rotate through one hash-verified execution path, using pinned
H100 entrywise-changed refactors and residual validation at 1e-8. No builds
or edits overlap timing.

Paired lifecycle median changes range from -3.023% to +0.973% incrementally
and from -2.450% to +1.110% against the resumed-pass baseline. No positive
change reaches the 2.5% review trigger. TSOPF_FS_b9_c6 incremental changes
are +0.721%, -0.767%, and +0.021% on CPUs 8–15, CPUs 0–7, and CPU 8;
cumulative changes are +0.450%, +0.183%, and -0.284%. These measurements do
not prove zero regression or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/snb-trace-trim-*`;
runner: `validate-snb-trace-trim.py`; preceding frozen executable:
`row-accessor-accepted-kls_bench`. The final audit verifies expected jobs
and counts, exact pending source diff, build provenance, runner/binary/matrix
hashes, complete commands, successful exits, residual validity, and all
summary medians against raw records.

The broader cleanup goal remains active. This is the first source chunk
since the completed row-helper cumulative corpus checkpoint.
