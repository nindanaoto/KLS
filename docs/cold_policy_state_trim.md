# Remove write-only cold policy state

Baseline: `168ded1`. Remove 17 private solver fields and their resets from
the cold tail of `kls_solver`, and correct the obsolete symbolic-identity
comment. Each removed identifier occurred only in its declaration and a
zero/NULL reset in `clear_matrix`; none was read by the source. The net
C-source reduction is 40 lines. Public statistics fields and their ABI are
unchanged. The PTS per-thread accumulator padding remains intact.

Removed fields:

- `balanced_moderate_hub_amd_cycle`
- `sparse_partial_diagonal_amd_btf_cycle`
- `dense_reciprocal_hub_metis_cycle`
- `sparse_symmetric_fragmented_metis_cycle`
- `symmetric_scalar_fringe_amd_lean_cycle`
- `symmetric_scalar_fringe_symbolic_identity`
- `pivoted_high_work_single_block_symbolic_cycle`
- `pivoted_high_work_single_block_identity`
- `low_work_hubbed_scalar_fringe_input_class`
- `symmetric_partial_diagonal_match_input_class`
- `symmetric_partial_diagonal_match_symbolic_rejected`
- `high_work_tiny_scalar_fringe_amd_symbolic_cycle`
- `high_work_tiny_scalar_fringe_amd_symbolic_identity`
- `asymmetric_bounded_degree_direct_metis_class`
- `asymmetric_bounded_degree_direct_metis_symbolic_identity`
- `near_symmetric_mega_hub_amd_symbolic_identity`
- `giant_dominant_hub_metis_dense_tail_symbolic_identity`

## Validation

Release CTest passes 6/6 (3.74 s); sanitizer CTest passes 6/6 (15.70 s).
Both builds are warning-free.

The 579-launch focused H100 comparison uses the accepted baseline and the
cumulative baseline `372f9d1`. All launches pass relative residual validation
at 1e-8. Incremental lifecycle changes range from -2.016% to +2.062%.
The initial cumulative LeGresley result (+3.375%) crosses the 2.5% review
threshold; the +2.062% TSOPF and +1.920% three-repeat twotone results also
receive independent confirmation rather than being assumed harmless.

The 243-launch confirmation has 24 triplets per TSOPF eight-thread CPU set,
nine twotone triplets, and 24 LeGresley triplets. All are valid:

| Case / CPUs | vs accepted lifecycle | vs cumulative lifecycle |
| --- | ---: | ---: |
| TSOPF_FS_b9_c6 / 8–15 | -0.110% | -1.056% |
| TSOPF_FS_b9_c6 / 0–7 | +0.902% | -0.914% |
| twotone / 8–15 | -0.295% | -0.056% |
| LeGresley_2508 / 8–15 | -0.037% | -0.152% |

The original cumulative flag is not reproduced. TSOPF steady solves in the
confirmation are 85.068 and 86.012 us, versus 86.743 and 86.682 us for
accepted. There is no renewed accumulator false-sharing regression and no
alignment adjustment was needed for this chunk.

All 822 timed launches have audited jobs, pinned commands, binary/matrix
hashes, source diff, build provenance, exits, residual validity and medians
recomputed from raw records. Executables rotate through one hash-verified
path; no builds or source edits overlap timing. These focused comparisons
are not a full-corpus or cross-machine performance guarantee.

Artifacts: `build/prep-trim-repair-dBcoha/cold-policy-trim-*` and
`cold-policy-confirm-repeat*`. Runner: `validate-cold-policy-trim.py`;
audits: `audit-focused-trim.py` and `audit-cold-policy-confirm.py`.
