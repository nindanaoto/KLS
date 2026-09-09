# Remove remaining inert cold policy fields

Baseline: `56016e7`. Remove eleven private cold-tail integer fields, their
zero resets, and obsolete comments (28 source lines). The removed state
covers sparse-spiked prediction, sparse-full-diagonal METIS, giant symmetric
scalar-fringe METIS, symmetric partial-diagonal matching, asymmetric
bounded-degree direct METIS, near-symmetric mega-hub AMD, and giant
dominant-hub METIS proposals. None of these private fields has a reader;
the only ordering assignments removed also stored literal zero.

Public statistics fields and their explicit zero-valued getter outputs are
unchanged. Live matching state and PTS accumulator alignment remain intact.

## Validation

Release CTest passes 6/6 (3.77 s), and ASan CTest passes 6/6 (15.77 s).
Both builds are warning-free. The focused rotating-order H100 experiment
compares this patch with accepted `56016e7` and cumulative `372f9d1`.
All 579 launches pass residual validation. Incremental lifecycle changes
range from -2.626% to +1.918%. TSOPF_FS_b9_c6 improves by 2.592% on
CPUs 8–15 and 1.309% on CPUs 0–7.

The initial mimo46x46_system result crosses the cumulative 2.5% review
threshold (+3.480%; +1.918% incremental). A separate 24-triplet confirmation
on CPUs 0–7 gives -0.136% incremental and +0.257% cumulative, with all
72 launches valid. The flag does not reproduce; no alignment fix is needed.

All 651 launches were audited for exact jobs, commands, source diff, build
provenance, binary and matrix hashes, exits, residuals, and recomputed
medians. No builds or source edits overlapped timing. This is focused
validation, not a full-paper or cross-machine performance guarantee.

Artifacts are under `build/prep-trim-repair-dBcoha/`, with prefixes
`remaining-cold-trim` and `remaining-cold-confirm-repeat`. The runner is
`validate-remaining-cold-trim.py`; the main audit is `audit-focused-trim.py`.
The confirmation uses `audit-eligibility-repeat.py` with the expected jobs
changed to `[["mimo46x46_system", "t8_96", 24]]` and total to 72.
