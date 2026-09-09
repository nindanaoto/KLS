# Remove output-only refinement diagnostics

Preceding baseline: `95e37be`; cumulative baseline: `372f9d1`.

Remove twelve `KLS_TRACE_REFINE` print blocks (101 source lines).
Residual checks, iterative refinement, GMRES/LSQR recovery, tolerance
selection, contract classification and acceptance decisions remain intact.
No test or paper runner consumes this output. The separate refinement
support diagnostic remains outside this chunk. No alignment hint is added.

Release and ASan CTest pass 6/6 (3.69 and 15.35 seconds). All 579 focused
launches pass residual validation at 1e-8. Incremental paired-median
lifecycle changes range from -2.318% to +3.535%; cumulative changes range
from -3.274% to +0.833%. Two cases cross the positive 2.5% review threshold:
transient on CPU 8 (+3.535%, two pairs) and twotone on CPUs 8–15
(+3.251%, three pairs).

Independent confirmations triple the pair counts: transient measures
+0.239% incremental and +0.441% cumulative across six pairs; twotone
measures +0.710% and +0.315% across nine pairs. All 45 confirmation
launches pass. Neither initial flag repeats at the review threshold.
These results do not establish zero loss or cross-machine equivalence;
no new full corpus screen was run.

All 624 launches rotate frozen binaries through one hash-verified execution
path, using pinned H100 entrywise-changed refactors with no concurrent
builds or source edits. Audits verify revision, exact pending diff, build
provenance, runner/binary/matrix hashes, expected jobs/counts, complete
commands, successful exits, residuals and recomputed medians.

Artifacts: `build/prep-trim-repair-dBcoha/refine-trace-trim-*` and
`refine-trace-repeat*`; runner: `validate-refine-trace-trim.py`; preceding
frozen executable: `compact-cache-accepted-kls_bench`.

The broader low-risk cleanup goal remains active.
