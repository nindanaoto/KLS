# Remove redundant private policy wrappers

Preceding baseline: `e069934`; resumed-pass baseline: `372f9d1`.

Remove four private forwarding layers without changing the underlying
predicates, numerical thresholds, or public interface:

- Call the checked-row environment predicate directly from its two callers.
- Call the existing native-panel structural predicate directly when recording
  AUTO eligibility.
- Replace the lazy-scatter helper with `!check_pivots` at its two callers.
- Give the diagonal-match profile implementation its existing caller-facing
  name, eliminating the forwarding wrapper and retaining all its callers.

This removes 20 net C-source lines. Release and ASan/UBSan builds succeeded,
with CTest passing 6/6 in each (3.70 and 15.64 seconds). The stripped binaries
still differ; comparison copies are in `build/policy-wrapper-yvdVRM/`.
The SNB kernel retains its 64-byte alignment at 0x64540, size 0x1931.

All 579 focused launches are valid: 36 onetone, 495 AUTO, and 48 forced-first.
Three versions rotate through a single hash-verified execution path, using
pinned H100 changed-refactor workloads and residual validation at 1e-8.
No builds or source edits overlap timing.

Paired lifecycle median changes range from -4.276% to +2.198% incrementally
and from -3.529% to +0.825% against the resumed-pass baseline. No positive
change reaches the 2.5% review trigger. The largest incremental increase is
twotone (+2.198%, three repetitions), whose cumulative change is -0.403%.
TSOPF_FS_b9_c6 incremental changes are +1.092%, +0.293%, and +0.101% on
CPUs 8–15, CPUs 0–7, and CPU 8 respectively. Their cumulative changes are
-0.045%, +0.470%, and -0.008%. These results do not prove zero regression or
cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/policy-wrapper-trim-*`;
runner: `validate-policy-wrapper-trim.py`; preceding frozen executable:
`recovery-gate-accepted-kls_bench`. The final audit verifies expected jobs
and counts, exact pending source diff, build provenance, runner/binary/matrix
hashes, complete commands and shared execution path, successful exits,
residual validity, and all summary medians against the raw records.

The broader cleanup goal remains active. This is the second source cleanup
since the completed TRSV cumulative corpus checkpoint.
