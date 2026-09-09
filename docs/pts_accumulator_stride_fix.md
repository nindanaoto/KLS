# Fix PTS accumulator false sharing and retain eligibility cleanup

This supersedes the rollback decision in
[eligibility_trim_regression.md](eligibility_trim_regression.md).
Baseline is `b12c76a`, whose production source and executable match `84b48b0`.

## Cause and repair

The twelve write-only eligibility fields were not needed for numerical
execution. Removing them reduced the solver allocation by 48 bytes and
perturbed subsequent heap placement. The PTS accumulation array was the
sensitive buffer: its per-thread slices were packed at `ntop` doubles.

TSOPF_FS_b9_c6 has 47 top entries, giving a 376-byte stride. In the rejected
layout, thread 0's slot 40 and thread 1's slot 0 share a 64-byte cache line.
The factor index streams map 1,998 and 1,999 potential updates to those
slots, respectively. Different threads update different values but contend
for ownership of the same line. The structural counts exclude no runtime
zero-xk skips. Other base offsets move the overlap to other thread pairs;
offset 48 happens to avoid overlap among written slots for this partition.

A 96-launch diagnostic confirmation changes only the accumulator offset
from 0 to 48 in the same rejected executable, with identical allocation
bases. Steady solves improve from 97.381 to 86.277 us on CPUs 8–15 and
99.328 to 88.968 us on CPUs 0–7. All results pass residual validation and
all non-timing/non-residual/non-memory result fields match. This pointer
intervention and the cross-thread write mapping identify false sharing;
hardware coherence counters were unavailable, so no cache-transfer count
is claimed. The temporary preload phase probe was unsuitable and is not
used as evidence.

The production repair does **not** select offset 48. It aligns the allocation
to 64 bytes and rounds every thread's stride to a cache-line multiple.
Both private accumulation and the existing ordered reduction use that
stride; their logical loops still use ntop. Thus no adjacent thread slices
share a 64-byte line, regardless of partition shape or allocator placement.
For this eight-thread case, the allocation grows from 3,008 to 3,072 bytes.

The stride helper checks rounding overflow, the allocation checks count and
byte-size overflow, and allocation failure retains the existing plan cleanup
and fallback. Zero-top plans retain a valid minimum slice. No numerical
arithmetic, ordering, acceptance threshold, thread policy, public option or
solver state field is added. Restore the entire 75-line eligibility cleanup;
including the repair, the net C-source reduction is 58 lines.

## Validation

Release CTest passes 6/6 (3.78 s); sanitizer CTest passes 6/6 (15.72 s).
Additional sanitizer H100 runs for TSOPF pass at 2, 4, 8 and 16 threads,
with every entrywise-changing system checked at relative residual 1e-8.
Both builds are warning-free.

The 144-launch native pilot compares accepted, rejected and fixed frozen
binaries with 24 rotating H100 triplets on each eight-thread CPU set:

| CPUs | Rejected steady solve | Fixed steady solve | Fixed lifecycle vs rejected | vs accepted |
| --- | ---: | ---: | ---: | ---: |
| 8–15 | 97.348 us | 87.483 us | -1.146% | +2.073% |
| 0–7 | 98.522 us | 86.665 us | -1.301% | -0.134% |

The remaining lifecycle variation in the first configuration is primarily
refactorization timing, not the recovered solve-phase loss.

The independent 579-launch focused suite compares accepted, cumulative
baseline `372f9d1`, and fixed. All launches are valid; none crosses the
positive 2.5% review threshold. Incremental lifecycle changes range from
-2.089% to +0.838%; cumulative changes range from -4.276% to +1.149%.
TSOPF eight-thread lifecycle changes are +0.838% and +0.665% versus accepted;
steady solve times are 87.051 and 86.117 us. One-thread TSOPF lifecycle is
-0.945%. This does not establish zero performance loss or cross-machine
equivalence, and it is not a new full-paper/full-corpus benchmark.

All 723 native timed launches pass audits for complete jobs, exact source
diff, frozen binary/matrix hashes, build provenance, pinned commands,
successful exits, residual validity and recomputed medians. No builds or
source edits overlap timing. Each launch uses one hash-verified execution
path and automatic policies without diagnostic overrides.

Artifacts: `build/prep-trim-repair-dBcoha/pts-padded-stride-*`.
Pilot runner: `validate-eligibility-align-pilot.py`; focused runner:
`validate-eligibility-trim.py`. Their `start` labels mean rejected cleanup
and cumulative baseline respectively. Audits: `audit-eligibility-repeat.py`
and `audit-focused-trim.py`. Causal diagnostic artifacts and the buffer
write-map analysis are recorded in `pts-accumulator-findings.md` under the
same build directory.
