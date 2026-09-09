# Remove unused fake clock

Preceding baseline: `0d6a4c0`; cumulative baseline: `372f9d1`.

Remove 13 C-source lines implementing KLS_FAKE_CLOCK and its cached flag
and atomic tick counter. CLOCK_MONOTONIC remains the clock source. No
tests or scripts use the removed switch. Real-time accounting, numerical
decisions and alignment hints are unchanged.

Warning-free Release and ASan builds pass CTest 6/6 each (3.60 and 15.47
seconds). All 579 focused launches pass residual validation at 1e-8.
Incremental paired-median lifecycle changes range from -1.823% to +0.907%;
cumulative changes range from -2.525% to +1.581%. No result reaches the
positive 2.5% review threshold. TSOPF_FS incremental changes are +0.427%,
-0.106%, and +0.010% on CPUs 8–15, CPUs 0–7, and CPU 8.
LeGresley_2508 measures +0.025%.

Sequential pinned H100 comparisons rotate frozen binaries through one
hash-verified execution path without overlapping builds or source edits.
The final audit verifies revision, exact diff, build provenance, hashes,
expected jobs/counts, complete commands, exits, residuals and medians.
The runner pins the new dependency tree and verifies the exact previously
reviewed tracker-removal patch against 93d2751. Incremental binaries share
that tree; the cumulative baseline predates the tracker removal, so that
comparison includes the documented dependency cleanup. All other build
provenance remains identical. No new full corpus run was made; this is
not zero-loss or cross-machine proof.

Artifacts: `build/prep-trim-repair-dBcoha/fakeclock-trim-*`;
runner: `validate-fakeclock-trim.py`; preceding frozen executable:
`fakeclock-accepted-kls_bench`.

The broader cleanup goal remains active.
