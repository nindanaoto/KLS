# Stop the low-risk slimming pass at the eligibility-field removal

Accepted source: `84b48b0`. The proposed deletion of twelve write-only
numeric-eligibility fields and their initialization writes removed 75 C
lines, but repeatedly slowed TSOPF_FS_b9_c6's post-refactor solves.
Withdraw that uncommitted candidate and all experimental alignment changes.
All previously committed slimming remains intact. This is a rollback of
one candidate, not a new optimization or an alignment fix.

## Evidence and limits

The original independent 144-launch confirmation measured lifecycle losses
of 2.473% on CPUs 8–15 and 1.738% on CPUs 0–7. Steady post-refactor solves
were 10.880% and 13.481% slower respectively. Unchanged numerical/statistical
fields in the inspected paired run and debugger observations do not show a
different numerical policy: steady solves select the existing PTS path.

The PTS worker retains the same address and size in accepted and rejected
binaries. Debugger observations show a 48-byte shift in numeric Xwork,
Udiag and compact permutation allocations. These are evidence of layout
sensitivity, not proof of a specific cache or false-sharing mechanism.
Hardware performance counters were unavailable, including through sudo.

Seven separate 64-byte alignment interventions did not recover the steady
solve loss: row worker entry, row solve dispatcher entry, refinement
workspace allocation, solver allocation, main solve entry, PTS solve entry,
and PTS accumulation allocation. Each pilot used 144 rotating-order launches
against accepted and rejected frozen binaries. The last four pilots were
audited during this investigation; the earlier three retain their prior
audits. No rejected alignment hint or allocation change remains.

This is the practical stopping point for the low-risk, alignment-only pass;
it does not establish that every possible alignment intervention must fail.
Further removal of these fields needs a separate, targeted investigation,
not acceptance of a persistent solve-phase loss hidden by lifecycle totals.

## Restored validation

The rebuilt Release executable is byte-for-byte identical to the frozen
accepted executable. Both SHA-256 hashes are:

`21c6c8bc0f689b1cb1fa83ecb80bcb95363f103f9eceed03b6b7df3e96946f86`

Release CTest passes 6/6 (3.74 s); sanitizer CTest passes 6/6 (15.65 s).
The final independent 144-launch H100 comparison uses entrywise-changing
values, automatic policies, one hash-verified execution path, rotating
order and pinned eight-thread CPU sets. No build or source edit overlaps
timing. All launches pass residual validation at 1e-8. The final audit
checks provenance, source diff, hashes, commands, complete jobs, exits,
residuals and recomputed medians.

| CPUs | Restored vs rejected lifecycle | Steady solve: rejected → restored |
| --- | ---: | ---: |
| 8–15 | -3.458% | 97.329 → 87.596 us |
| 0–7 | -1.034% | 96.831 → 86.564 us |

Against the byte-identical accepted binary, paired lifecycle changes are
+0.072% and -0.943%, illustrating measurement variation. No new full-corpus
or cross-machine benchmark was run.

Artifacts are under `build/prep-trim-repair-dBcoha/`: original
`eligibility-trim-repeat*`, alignment pilots `eligibility-align-pilot*`,
`eligibility-dispatch-align*`, `eligibility-workspace-align*`,
`eligibility-solver-align*`, `eligibility-solve-entry-align*`,
`eligibility-pts-align*`, `eligibility-pts-acc-align*`, and final
`eligibility-restored-repeat*`. Pilot/restoration runner:
`validate-eligibility-align-pilot.py`; audit: `audit-eligibility-repeat.py`.
In this runner, `original` is accepted, `start` is the rejected cleanup,
and `fix` is the tested intervention or restoration—not a cumulative baseline.
