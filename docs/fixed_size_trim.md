# Remove fixed size arguments

Baseline: `d359c0f`. Remove two always-zero dense-column arguments from
cached-supernode/CBLAS vector eligibility helpers and eliminate their zero
terms. Remove the extra-count arguments from the LU/input seed reserve
helpers: both have one caller and always reserve one additional entry.
Retain null, count-overflow, allocation-size, and allocation-failure checks.
Net reduction: 12 source lines. No numerical policy or public API changes.

The improved scan includes cast-bearing calls. The symmetry helper's fixed
degree-bound argument is retained: it explicitly communicates the caller's
allocation bound, and removing it would only relocate that constant.

## Validation

Release CTest passes 6/6 in 4.13 s; ASan passes 6/6 in 17.25 s. Builds succeed
without reported warnings; `git diff --check` passes. Generated code differs,
so the final source receives the focused H100 timing campaign.

All 579 launches pass the provenance/result audit. Median paired lifecycle
changes range from -4.132% to +1.202% versus accepted and from -4.073% to
+1.625% versus cumulative baseline `372f9d1`. No case reaches the 2.5% review
threshold. The largest apparent gain is the two-triplet transient case;
do not interpret that small sample as a demonstrated optimization.
No alignment fix is needed.

The audit verifies exact jobs, commands, source revision/diff, provenance,
binary/matrix hashes, exits, residuals, and recomputed medians. Source and
binaries remain fixed during timing. Focused results do not establish
full-paper or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/fixed-size-trim-*`.
Runner: `validate-fixed-size-trim.py`. Audit: `audit-focused-trim.py`.
