# Remove kernel construction and step profiling

Baseline: `73245ce`. Remove `KLS_CONSTRUCT_PROF` timing blocks and their
enable state, construction-time and step-time TLS counters, serial profile
logging, and four corresponding 32-bit wrapper aliases. Remove the deleted
construction-time value from the separate pipeline log while retaining its
live count fields. Net reduction: 43 source lines.

Symbolic solves, column construction, numeric solves, pivot selection,
pipeline scheduling, and numerical controls are unchanged. Pipeline phase
and wall-time profiling remain for separate review. Whole-source inspection
finds no remaining references to the removed profiler or counters.

## Validation

Release CTest passes 6/6 (4.13 s); ASan CTest passes 6/6 (17.27 s).
Both builds are warning-free; `git diff --check` passes.

All 579 focused H100 launches pass residual validation. Median paired
lifecycle changes range from -4.161% to +2.135% versus accepted and from
-2.515% to +1.807% versus cumulative `372f9d1`. TSOPF FS changes +0.223%
and -0.468% on the two eight-thread CPU sets; onetone1 changes -0.125%
and -0.153% on the two single-thread configurations.

The initial four-triplet circuit_4 first-factor result (+2.135%) receives
an independent 24-triplet first-factor confirmation on CPUs 8–15. All 72
launches pass. Repeated lifecycle change is -0.280% incremental and +0.292%
cumulative. The initial slowdown does not reproduce; no alignment fix is
needed.

All 651 launches have audited exact jobs, pinned commands (including
first-factor mode), source diff, build provenance, executable/matrix hashes,
exits, residuals, and recomputed medians. The reviewed dependency patch is
explicitly identified by SHA-256:
`d2b1ae14824c1ff0dd027810bb7311e08fce5b11d67086688fd5818b007df332`.
No builds or source edits overlapped timing. Focused validation does not
establish full-paper or cross-machine performance equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/kernel-construction-trim-*` and
`kernel-construction-confirm-first*`; runner `validate-kernel-construction-trim.py`.
Main audit: `audit-kernel-construction-trim.py`. Confirmation uses
`audit-eligibility-repeat.py` with group `first`, one 24-triplet circuit_4 job,
total 72, and the same explicit dependency-patch hash.
