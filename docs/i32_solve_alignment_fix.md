# Resolve the row-acceptance cleanup regression

Accepted source: `12bc310`; rejected experiment: `row-accept-repeat`;
cumulative baseline: `372f9d1`. This supersedes the stopping point recorded
in [row_accept_trace_regression.md](row_accept_trace_regression.md).

Restore the complete 119-line diagnostic removal and align `kls_i32_solve`
to 64 bytes on GCC/Clang. The alignment hint and explanatory comment add
five lines, for a net reduction of 114 C-source lines. No matrix, dimension,
thread-count, CPU-model, numerical threshold, or solver policy is added.
Other compilers retain their ordinary alignment behavior.

## Evidence for the cause

All non-timing/non-residual result fields match between the accepted and
rejected binaries in the original confirmation. Diagnostic launches also
select the same compact i16 column solve path, not the prepared single-block
row solve. The principal slowdown was the solve after each refactorization.
The benchmark's `refactor_solve_steady_seconds_avg` measures those solves
alone; it is not the combined refactor-plus-solve time.

Disassembly of `kls_i32_solve` contains the same 2,962 normalized instructions
and size 0x374c in all three binaries. Only relocation-dependent operands
are normalized (local branches retain their function-relative offsets):

| Version | Entry address | Offset within 64 bytes |
| --- | ---: | ---: |
| Accepted | 0x11acf0 | 48 |
| Rejected cleanup | 0x11aad0 | 16 |
| Cleanup plus alignment | 0x11ab00 | 0 |

The alignment-only intervention recovers the loss without changing numerical
execution. This supports code-placement sensitivity; the specific processor
front-end/cache mechanism has not been established. The hint stabilizes
entry alignment, not absolute address or all possible placement effects.

## Validation

The 144-launch pilot rotates accepted, rejected, and fixed binaries through
one hash-verified execution path, with 24 H100 triplets per cache domain.
All non-timing/non-residual fields match across the three versions.

| CPUs | Fixed vs rejected lifecycle | Fixed vs accepted lifecycle | Steady solve: accepted / rejected / fixed |
| --- | ---: | ---: | ---: |
| 8–15 | -4.630% | -1.521% | 10.888 / 12.429 / 10.856 us |
| 0–7 | -1.796% | +1.691% | 11.193 / 12.689 / 11.195 us |

The separate 579-launch focused suite passes without a positive 2.5% review
flag. LeGresley_2508 on CPUs 8–15 measures +0.117% versus accepted and
-1.580% versus the cumulative baseline. Across the suite, incremental
changes range from -1.397% to +2.331% (twotone, three repetitions);
cumulative changes range from -2.813% to +1.053%. This is not proof of
zero loss or cross-machine performance equivalence. A new full corpus
screen was not run for this fix.

Release and ASan CTest pass 6/6 each (3.73 and 15.56 seconds). All 723
timed launches pass residual validation at 1e-8. No builds or source edits
overlap timing. Final audits verify frozen hashes, build provenance, exact
pending source diff, expected jobs/counts, complete pinned commands,
successful exits, residual validity, and medians recomputed from raw data.
The existing SNB alignment remains intact at 0x64040, size 0x1931.

Artifacts: `build/prep-trim-repair-dBcoha/i32-align-*`, including
`i32-align-pilot*`; runners: `validate-i32-align.py` and
`validate-i32-align-pilot.py`. In pilot metadata, `original` is accepted,
`start` is rejected, and `fix` is aligned. In the focused suite, `start`
is the resumed-pass baseline instead.
