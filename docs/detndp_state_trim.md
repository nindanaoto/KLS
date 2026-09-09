# Remove redundant deterministic NodeNDP routing state

Baseline: `3edfd2c`. Remove 39 net lines: `kls_detndp_class_ok`,
`kls_det_ndp_requested`, their stores and obsolete comments. The sole
installation of the private METIS ordering callback is in
`analyze_with_ordering`, which sets the former flag to one before invoking
the synchronous vendored KLU analysis. KLU calls the callback directly on
that thread. The flag has no other assignments or address escapes.
Consequently the request flag ORed with it cannot affect routing.

Keep thread-count/resource checks, generic context forcing,
`KLS_FORCE_DET_NDP`, and `KLS_DISABLE_DET_NDP`. Remove the ineffective
`KLS_DISABLE_PS_MTND` read associated solely with the redundant request flag.
No numerical or ordering decision changes.

## Validation

Release CTest passes 6/6 in 4.07 s; ASan passes 6/6 in 17.17 s.
Both builds succeed without reported warnings; `git diff --check` passes.
Code and TLS layout change, requiring timing checks.

All 579 initial focused H100 launches pass correctness and provenance audits.
Initial twotone changes +1.355% against accepted and +2.886% against cumulative
baseline `372f9d1`, triggering a repeat. ASIC_680ks is also repeated because
its initial accepted comparison is near the review threshold at +1.894%.

The longer 99-launch confirmation keeps source and binaries unchanged:

| Case | Rotating-order triplets | Change vs accepted | Change vs cumulative |
| --- | ---: | ---: | ---: |
| twotone | 9 | -0.665% | +0.309% |
| ASIC_680ks | 24 | -0.457% | +0.282% |

The initial slowdown does not persist. No alignment fix is needed. All 678
launches pass validation; audits verify jobs, commands, source revision/diff,
build provenance, binary/matrix hashes, successful exits, residuals, and
recomputed medians. Focused results do not establish full-paper or
cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/detndp-state-trim-*` and
`detndp-state-confirm-*`. Runner: `validate-detndp-state-trim.py`.
Audit: `audit-focused-trim.py`; confirmation uses `audit-eligibility-repeat.py`
with explicit twotone/ASIC jobs and total 99.
