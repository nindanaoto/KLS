# Remove unused private proposal state

Baseline: `ad451cd`. Remove eight write-only private fields and their
assignments, saving 33 net source lines. Six fields were only declared and
reset: `medium_spike_minfill_candidate`, `large_sparse_amf3_candidate`,
`large_bounded_no_btf_amf_candidate`,
`partial_diagonal_many_block_no_btf_cycle`,
`sparse_diagonal_row_hub_no_btf_cycle`, and
`large_reciprocal_hub_amd_btf_cycle`.

The other two, `nearly_missing_diagonal_early_match_selected` and
`compact_missing_diagonal_match_candidate`, only cached local proposal
flags. Their right-hand sides have no side effects. Remove the now-empty
conditional and update the adoption comment. Keep the live
`compact_missing_diagonal_match_selected` state and all matching decisions.
Public statistics fields and ABI are unchanged.

## Validation

Release CTest passes 6/6 (4.11 s); ASan CTest passes 6/6 (17.22 s).
Both builds are warning-free. All 579 focused H100 launches validate.
Median paired lifecycle changes range from -1.044% to +2.263% versus
accepted and -2.139% to +2.265% versus cumulative `372f9d1`.

An independent 243-launch confirmation repeats TSOPF FS on both CPU sets,
TSOPF RS on CPUs 8–15 (24 triplets each), and twotone on CPUs 8–15 (nine
triplets). All launches validate:

| Case / CPUs | vs accepted | vs cumulative |
| --- | ---: | ---: |
| TSOPF_FS_b9_c6 / 8–15 | +2.113% | +0.366% |
| TSOPF_FS_b9_c6 / 0–7 | +1.211% | -0.626% |
| TSOPF_RS_b9_c6 / 8–15 | +1.637% | +0.713% |
| twotone / 8–15 | -0.716% | +0.808% |

The TSOPF increases persist; they are not dismissed as noise. FS steady
refactor medians rise 2.958% / 1.167%, and steady refactor-solve solve
medians rise 3.192% / 3.482%. The overall lifecycle changes remain below
the established 2.5% review threshold. Accept this chunk with that measured
tradeoff recorded, not as a claim of identical performance. No alignment
adjustment is included.

All 822 launches have audited exact jobs, commands, hashes, source diff,
build provenance, exits, residuals, and recomputed medians. Source and
binaries stayed fixed throughout timing. These are focused comparisons,
not a full-paper or cross-machine performance guarantee.

Artifacts: `build/prep-trim-repair-dBcoha/private-proposal-trim-*` and
`private-proposal-confirm-repeat*`; runner `validate-private-proposal-trim.py`.
The main audit is `audit-focused-trim.py`. The confirmation uses
`audit-eligibility-repeat.py` with the four jobs listed above and total 243.
