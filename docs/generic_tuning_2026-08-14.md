# Generic KLS tuning follow-up (2026-08-14)

This note records the first post-external-validation tuning revision.  The
external-validation comparison itself remains frozen; every result below is
development or post-reveal diagnostic evidence from the same machine.

## Attribution

The targeted six-arm route ablation covered `fd18`, `usps_norm_5NN`, `wang3`,
`ted_A_unscaled`, `k3plates`, `rma10`, and `airfoil_2d` at eight threads with
entrywise 0.001 value updates.  Disabling the row executor or fast first
factor changed the geometric mean by less than about one percent.  The larger
regret came from two generic policy failures:

1. AUTO excluded every graph below 30,000 rows from its already guarded
   repeated-lifecycle NodeND portfolio, regardless of realized symbolic work.
2. A stable unscaled numeric could run three complete scale trial factors on
   the first changed update even when it used the requested diagonal pivot
   sequence exactly and had no nudge, perturbation, or conditioning signal.

## Retained changes

The generic ND portfolio now admits sub-30k graphs to the existing lifecycle
economics test.  Because a compact separator tree can lose locality despite a
symbolic near-tie, the newly admitted path must remove at least five percent
of retained symbolic storage; the existing decisive ten-percent Pareto arm
and the established large-graph policy are unchanged.  A diagnostic
`KLS_DISABLE_SUB_30K_GENERIC_ND=1` switch restores the old cutoff.

AUTO scale trials are now skipped when the retained unscaled factor has zero
off-diagonal pivots, zero pivot nudges, zero KLS perturbations, and estimated
`rcond >= 1e-6`.  Difficult, pivoted, perturbed, or poorly conditioned
numerics remain eligible.  The diagnostic
`KLS_DISABLE_STABLE_UNSCALED_SCALE_GUARD=1` switch disables this guard.

## Targeted timing evidence

All times below are modeled H100 SPICE cycles in seconds.  The old column is
the median of the frozen post-reveal AUTO ablation; the new column is the
median of three runs after both changes.

| Matrix | Old AUTO | New AUTO | Ratio | Final route |
| --- | ---: | ---: | ---: | --- |
| `airfoil_2d` | 0.834 | 0.527 | 0.632 | METIS |
| `wang3` | 7.356 | 4.637 | 0.630 | METIS |
| `fd18` | 2.198 | 2.212 | 1.007 | AMD, unchanged |
| `usps_norm_5NN` | 6.081 | 5.809 | 0.955 | AMD, unchanged |
| `rma10` | 2.930 | 2.901 | 0.990 | METIS, unchanged |
| `ted_A_unscaled` | 0.264 | 0.267 | 1.010 | AMD, unchanged |
| `k3plates` | 0.477 | 0.519 | 1.088 | AMD, unchanged |

The `k3plates` timing difference is run noise rather than a retained route
change: an initial permissive version selected a marginal ND proposal and was
5.8 percent slower in seven counterbalanced pairs.  Its symbolic storage win
was only 0.7 percent, which motivated the frozen five-percent admission
margin.  The final policy rejects that proposal and restores exactly the AMD
fill and flop counts.

The scale guard was isolated with seven counterbalanced forced-METIS pairs on
`wang3`.  Enabling it reduced H100 from 8.851 to 5.230 seconds (median paired
ratio 0.590) and first changed-refactor time from 3.776 to 0.137 seconds.
Both arms retained scale -1, 7,699,744 factor entries, 4.594 billion reported
flops, and the same residual quality.

## Generalization evidence

The 24-family development manifest completed 16 matrices.  Three matrices
timed out at 120 seconds and five reported singular/factor failures; these
failures are retained in `build/generic-policy-dev-enabled.failures`.  Only
one successful sub-30k development matrix, `TEM27623`, reached the new
portfolio.  Its challenger was rejected and both old/new arms ended with the
same METIS numeric.  Seven counterbalanced pairs gave a 1.006 median ratio;
the new path's approximately 0.15-second rejected analysis is the measured
cost of broader exploration on that case.

After freezing the five-percent margin, the group-disjoint 24-family holdout
completed 24/24 without a crash or solver-status failure.  Two sub-30k
holdouts, `crystk02` and `poisson3Da`, selected METIS.  Their corresponding
historical H100-model rows were 7.92 and 4.78 seconds; the frozen-policy run
reported 4.88 and 2.92 seconds.  This historical comparison is directional,
not a counterbalanced A/B, so it is not used to tune the selector.

Every newly affected route was rerun for twenty entrywise updates with
`KLS_BENCH_VERIFY_EACH_REFACTOR=1`:

| Matrix | Max changed-update relative residual |
| --- | ---: |
| `airfoil_2d` | 5.16e-15 |
| `wang3` | 1.67e-15 |
| `crystk02` | 6.04e-14 |
| `poisson3Da` | 2.50e-15 |

The broad holdout harness does not enable that residual audit by default and
therefore must not be described as a 24-matrix accuracy certificate.  In
particular, its initial solve exposed poor quality on unrelated `shyy41`
(relative residual about 5.7e7) and `shermanACb` (about 2.4e-3).  Those routes
were not selected by either new policy, but they remain correctness work.

## Next kernel work

These changes reduce avoidable work; they do not close the CKTSO/SubtreeLU
gap.  The external profile still points to recurring numeric factorization as
the dominant cost.  The next high-value experiment should therefore measure
and improve the persistent refactor representation: cache-friendly row or
segment storage, supernode/panel update batching, and separator-aware task
granularity.  Work-reducing ordering remains important, but another matrix
shape classifier would be weaker evidence than improving those shared
kernels.
