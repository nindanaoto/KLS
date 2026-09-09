# Remove redundant retired-statistics overlays

Baseline: `dedcc68`. Remove 53 literal-zero overlays from `kls_get_stats`
(342 production source lines). Their stored public statistics members have
no writers: solver creation uses calloc, matrix clearing zeros the statistics,
and ordering rollback restores the solver's own saved statistics. The getter
already copies those zero bytes. Four identifiers also name private fields;
those are separate storage, not writers to the public statistics members.

All public fields, ABI sizes, bounded-copy behavior, and computed overlays
remain unchanged. In particular, compact solve, singleton-run, promoted
tolerance and verified-PTS statistics still use their existing calculations.

## Compatibility and correctness

The smoke test gains 38 lines checking all 53 retired fields with poisoned
output storage. It checks every requested byte length from the size header
through the full struct for writes beyond the requested prefix, and compares
partial copies in the retired tail with a full snapshot. These checks run
after creation, analysis, factorization, changing-value refactors, solves,
and reanalysis, across serial/KLS and normal/transposed orientations.

Release CTest passes 6/6 (4.07 s); ASan CTest passes 6/6 (17.28 s).
Both builds are warning-free and `git diff --check` passes.

## Performance validation

The focused 579-launch H100 experiment compares accepted `dedcc68` and
cumulative `372f9d1`. All runs pass residual validation. Median paired
lifecycle changes range from -1.346% to +1.951% versus accepted and from
-2.420% to +0.265% versus cumulative. TSOPF_FS_b9_c6 changes by +1.268%
on CPUs 8–15 and -0.300% on CPUs 0–7, remaining faster than cumulative.

Although no case crosses the 2.5% review threshold, the two near-2% results
receive additional confirmation: nine twotone triplets and 24 LeGresley_2508
triplets on CPUs 8–15. These 99 launches all pass. Twotone measures +0.266%
versus accepted / -1.024% cumulative; LeGresley measures -1.105% / -1.817%.
The initial near-2% slowdowns do not reproduce. No alignment fix is needed.

All 678 launches have audited commands, exact jobs, source diff, build
provenance, matrix/binary hashes, validity and recomputed medians. No source
edits or builds overlap timing. This focused validation is not a full-paper
or cross-machine guarantee.

Artifacts: `build/prep-trim-repair-dBcoha/stats-overlay-trim-*` and
`stats-overlay-confirm-repeat*`; runner `validate-stats-overlay-trim.py`.
Main audit: `audit-focused-trim.py`. Confirmation audit:
`audit-eligibility-repeat.py`, with expected jobs changed to
`[["twotone", "t8_32", 9], ["LeGresley_2508", "t8_32", 24]]` and total 99.
