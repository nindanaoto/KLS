# Simplify the opt-in first-factor recovery gate

Preceding source baseline: `6f16d4a`; resumed-pass baseline: `372f9d1`.

Remove the private `kls_should_try_first_factor_recovery` wrapper and replace
its two callers with `kls_first_factor_env_enabled()`. The wrapper returned
true only for that flag; its explicit-disable check and default branch both
returned false, and it never used its solver argument. Ordinary first-factor
AUTO selection is unchanged. Recovery remains opt-in. This removes 17 net
C-source lines without changing numerical thresholds or arithmetic.

Release and ASan/UBSan builds succeeded, with CTest passing 6/6 in each.
The focused comparison completed 579 valid launches: 36 onetone, 495 AUTO,
and 48 forced-first controls. Three versions rotate through one shared,
hash-verified execution path with pinned H100 changed-refactor workloads
and residual validation at 1e-8. No builds or edits overlap timing.

Paired lifecycle median changes range from -3.118% to +1.309% against the
preceding source and from -2.272% to +1.436% against the resumed-pass baseline.
No control reaches the positive 2.5% review trigger. These measurements do
not prove zero regression or cross-machine equivalence.

Artifacts are under `build/prep-trim-repair-dBcoha/recovery-gate-trim-*`;
runner: `validate-recovery-gate-trim.py`; preceding frozen executable:
`trsv-gate-accepted-kls_bench`. The final audit verifies launch counts, jobs,
runner/binary/matrix hashes, exact pending source diff, shared execution
path, forced-first arguments, successful exits, residual validity, and
paired lifecycle medians recomputed from raw records.
