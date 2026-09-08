# Remove redundant row access and pivot-check wrappers

Preceding baseline: `d8bdac8`; resumed-pass baseline: `372f9d1`.

Replace the private target-position accessor with direct reads of its sole
array at all four callers. Consolidate the rowwise-pivot acceptance check
under one private name with an explicit starting-block argument: seven
callers retain block zero, and both repaired-block callers retain their
nonzero-capable argument. Remove the unused PTS block-size local and its
void cast. Bounds checks, rejection recording, numerical predicates, and
the mixed-width input-position accessor remain unchanged. This removes
13 net C-source lines.

Release and ASan/UBSan builds succeed; CTest passes 6/6 in each (3.70 and
15.60 seconds). Stripped executable instructions differ, so this chunk
received a fresh focused timing comparison. Comparison copies are under
`build/row-accessor-88k0Bb/`. The SNB kernel retains 64-byte alignment at
0x64540, size 0x1931.

All 579 focused launches are valid: 36 onetone, 495 AUTO, and 48 forced-first.
Three versions rotate through one hash-verified execution path, using pinned
H100 changed-refactor workloads and residual validation at 1e-8. No builds
or source edits overlap timing.

Paired lifecycle median changes range from -1.914% to +0.943% incrementally
and from -2.980% to +1.428% against the resumed-pass baseline. No positive
change reaches the 2.5% review trigger. TSOPF_FS_b9_c6 incremental changes
are +0.761%, -0.605%, and +0.117% on CPUs 8–15, CPUs 0–7, and CPU 8;
cumulative changes are +0.707%, -1.253%, and +0.026%. These results do not
prove zero regression or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/row-accessor-trim-*`;
runner: `validate-row-accessor-trim.py`; preceding frozen executable:
`predicted-wrapper-accepted-kls_bench`. The final audit verifies expected
jobs/counts, pending source diff, build provenance, runner/binary/matrix
hashes, complete commands, successful exits, residual validity, and all
summary medians recomputed from raw records.

Four source chunks have removed 68 net C-source lines since the last
completed TRSV cumulative corpus checkpoint. Repeat that shared-path
corpus gate before further source deletions. The broader cleanup goal
remains active.
