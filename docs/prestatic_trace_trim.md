# Remove prestatic diagnostics

Preceding baseline: `d60f2cd`; cumulative baseline: `372f9d1`.

Remove 263 lines of output-only prestatic logging, its trace helper and
local diagnostic timestamp. Matching, quality gates, candidate selection,
live timing and synchronization remain intact. Existing alignment hints
are unchanged. No test or paper runner consumes the removed output.

Warning-free Release and ASan builds pass CTest 6/6 (3.61 and 15.50 seconds).
All 579 focused launches pass residual validation at 1e-8. Incremental
paired-median lifecycle changes range from -1.610% to +2.451%; cumulative
changes range from -2.306% to +1.133%. No case reaches the positive 2.5%
review threshold. TSOPF_FS incremental changes are +0.743%, -0.075%, and
-0.234% on CPUs 8–15, CPUs 0–7, and CPU 8. LeGresley_2508 measures -0.646%.

Because gemat11's initial +2.451% result was close to the threshold, repeat
it independently with 24 triplets on CPUs 0–7. The repeat measures -0.048%
incrementally and +0.138% cumulatively, with all 72 launches valid. The
initial slowdown does not reproduce. No alignment change is needed.
These results do not establish zero loss or cross-machine equivalence;
no new full corpus screen was run.

Frozen binaries rotate through one hash-verified execution path with
pinned H100 entrywise-changed refactors and no concurrent builds or source
edits. Audits of all 651 launches verify revision, exact pending diff,
build provenance, runner/binary/matrix hashes, expected jobs/counts,
complete commands, successful exits, residuals and recomputed medians.

Artifacts: `build/prep-trim-repair-dBcoha/prestatic-trace-trim-*` and
`prestatic-trace-repeat*`; runner: `validate-prestatic-trace-trim.py`;
preceding frozen executable: `ordering-trace-accepted-kls_bench`.

The broader low-risk cleanup goal remains active.
