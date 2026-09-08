# Rejected ETree timing trial and cumulative stopping regression

Preceding baseline: `c4aa02e`; cumulative baseline: pushed `88a177e`.

The trial deleted `KLS_TRACE_FAST_REJECT_REPAIR_TIMING`, its local timer state, five
phase checkpoints, and final diagnostic print in
`kls_try_row_first_etree_ready_tail_repair`. None of these times is consumed
by repair admission, numerical work, failure handling, or performance
selection. Those behaviors remain intact. No public fields or private
structure layouts changed. Proposed C-source reduction: 39 lines. Historical
experiment notes mentioning the removed switch remain historical records.

Release and ASan/UBSan CTest each pass 6/6, with no build warnings or errors.
The focused H100 campaign rotates preceding (`original`), cumulative
(`start`), and candidate (`fix`) executables. It records matching build
provenance, commands, hashes, source diff, JSONL results, and summaries under
`build/prep-trim-repair-dBcoha/repair-timing-*`. Every entrywise-changed
refactor is numerically validated at a 1e-8 residual limit; build/test work
does not overlap timed campaigns. This is not a full paper or Xyce rerun.

All 192 forced-KLS-first and 309 automatic-policy launches pass numerical
validation. Forced-first median paired lifecycle changes range from -0.51%
to +0.19% incrementally, and -1.72% to +0.21% cumulatively. Automatic-policy
changes range from -1.31% to +1.06% incrementally, and -3.04% to +2.52%
cumulatively. TSOPF eight-thread changes are +1.06% and +0.21% incrementally.
The fixtures do not establish timed coverage of the ETree repair routine;
these are lifecycle/layout regression controls, not an ETree speed claim.

The onetone1 cumulative +2.52% screen prompted a 12-triplet serial repeat
(`repair-timing-onetone-recheck*`). All 36 launches pass numerically, but
the candidate is slower than pushed `88a177e` in **12/12 lifecycle pairs**:
median +2.77%, with paired changes ranging from +2.12% to +4.83%. A paired
percentile bootstrap of the median (10,000 resamples, Python random seed 42)
gives a descriptive 95% interval of +2.37% to +3.29%. Steady refactoring is
+3.23% in the median, positive in all 12 pairs; initial factorization and
refactor solve do not show a comparable slowdown. This is evidence of a
repeatable cumulative regression, not proof of its microarchitectural cause.

The latest candidate alone is +0.67% lifecycle versus `c4aa02e` (9/12 pairs
slower); the entire cumulative effect must not be attributed to 39 timer
lines. The original stop rule was approximately 3% repeatable lifecycle
regression, not an exact 3.000% cutoff. The +2.77% repeat, phase evidence,
and uniformly slower pairs meet that intended stopping point. Stop further
slimming instead of treating the near-threshold result as permission to
continue deleting code.

The uncommitted timer removal is reverted. All accepted chunks through
`c4aa02e` remain. This does **not** claim the retained baseline has no
cumulative slowdown: the preceding binary in this repeat is already +2.30%
slower than `88a177e` in median paired lifecycle time (12/12 pairs slower).
Locating and repairing that cumulative regression is a
separate next task; it is not silently folded into this stop-at-regression
cleanup goal. Remaining scalar/producer touch counters are not removed.

The candidate source diff matched the campaign's recorded diff before the
revert. Post-revert Release and ASan/UBSan CTest each pass 6/6. The rebuilt
release executable is byte-identical (SHA-256) to the frozen `c4aa02e`
baseline. Library and script sources have no uncommitted changes.

## Cleanup goal completion audit

The user asked for low-risk slimming in meaningful commits until a notable
regression is encountered. Seven accepted source commits from `55d1f6b`
through `c4aa02e` remove 868 net C-source lines versus `88a177e`; each has
its own committed scope and validation report. No numerical algorithm,
TSOPF residual recovery, or worker-selection policy was intentionally removed.
The repeatable approximately 3% cumulative onetone1 slowdown above supplies
the requested stopping condition. The uncommitted trial was rejected,
accepted work was retained and reverified, and no further deletion proceeded.
Thus this goal ends at the regression condition, not with a claim that every
remaining cleanup candidate has been exhausted. No commits were pushed.
