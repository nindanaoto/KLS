# Remove preparation and repair-budget diagnostics

Preceding baseline: `2f7066a`; cumulative baseline: `372f9d1`.

Remove 54 net source lines of output-only repair dispatch, row task-floor,
component calibration, row-solve separator and fast-factor budget logging.
Discard only the diagnostic local copy of the separator-plan result; the
builder still runs. Calibration data, measured counts, numerical decisions,
live timing and synchronization remain intact. Existing alignment hints
are unchanged.

Warning-free Release and ASan builds pass CTest 6/6 (3.65 and 15.46 seconds).
All 579 focused launches pass residual validation at 1e-8. Incremental
paired-median lifecycle changes range from -2.811% to +0.468%; cumulative
changes range from -2.771% to +0.805%. No case reaches the positive 2.5%
review threshold. TSOPF_FS incremental changes are -0.037%, -0.177%, and
+0.037% on CPUs 8–15, CPUs 0–7, and CPU 8; LeGresley_2508 measures -1.135%.
These results do not establish zero loss or cross-machine equivalence;
no new full corpus screen was run.

Frozen binaries rotate through one hash-verified execution path with
pinned H100 entrywise-changed refactors and no concurrent builds or source
edits. The final audit verifies revision, exact pending diff, build
provenance, runner/binary/matrix hashes, expected jobs/counts, complete
commands, successful exits, residuals and recomputed medians.

Artifacts: `build/prep-trim-repair-dBcoha/preparation-trace-trim-*`; runner:
`validate-preparation-trace-trim.py`; preceding frozen executable:
`dense-trace-accepted-kls_bench`.

The broader low-risk cleanup goal remains active.
