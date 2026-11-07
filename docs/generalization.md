# AUTO generalization qualification

The generalized AMD admission and optimized block detector are now the sole
implementations, enabled by default at the user's request. The build-time
opt-in switch and superseded implementations have been removed. Numerical
safeguards remain in place. This promotion is a user decision; it does not
change the recorded qualification results below. Remaining ordering rules and
stages 2 through 4 are still pending.

The control is the post-cleanup source frozen in
`build/generalization-20261106/control/source.json`. Its hash includes the
uncommitted cleanup; HEAD alone does not identify it. The adjacent snapshots,
raw observations, failed screens, and test logs are retained under
`build/generalization-20261106/`. `generalization_inventory.json` records the
implemented prototype, pending replacements, and capability checks.

The first prototype replaces the AMD input envelope, narrow symbolic fill/work
bands, and one-shot row/nonzero shortcut with a per-analysis cost reservation.
It retains the initial AMD allocation when another symbolic portfolio cannot
amortize setup. Both orientations share the ledger. Current BTF and block admission remain
baseline selection steps for this AMD change: when matching creates a different
one-block coordinate frame, the current comparison of actual BTF/no-BTF
symbolics settles the incumbent before optional ordering exploration. Each
coordinate frame is analyzed once; a congruent retained ordering does not need
a duplicate AMD pass. A concrete block proposal can avoid a disposable BTF,
and discovery decisions are carried into factor entry instead of repeated. Its
broader admission and rollback protocol remain separate pending work.

The block detector uses a transient byte limit, degree selection, and the
actual component decomposition. It rejects unsuitable degree distributions
before allocating symmetric adjacency and rejects impossible components before
RCM construction. Current production block admission is retained in both build modes until the
separate matching/block stage qualifies. The experimental detector primitive
supports a byte-bounded larger domain, but the AMD stage does not activate that
wider proposal domain. Proposals that discarded one-block matching were rejected, including a
later exact-symmetry variant: neither established numeric payback. The current
prototype retains the actual representation comparison. Explicit orderings
retain their requested behavior.

The portable cost rates are estimates, not a guarantee of 2% overhead. An
uninterruptible ordering can exceed its reservation; elapsed cost exhausts later
credit. Missing arithmetic estimates remain unknown, and
supplied-order KLU analysis now initializes that field consistently in both
integer widths. Initial factorization, first changed-value setup, warm
refactors, solves, input preparation, and certification have separate model
terms. Allocation and analysis failure still permit ordinary symbolic recovery.
The remaining speculative paths must be integrated in later stages.

This command fits an optional model from generated sparse
grids, irregular graphs, and fragmented blocks:

```sh
kls-autotune --policy-only --threads 8 --cpus 0,1,2,3,4,5,6,7 --output policy.conf
```

A fixed eight-pair held-out synthetic validation uses a 100-system changed-value
lifecycle and both execution orders. It requires aggregate nonregression and
no greater-than-2% per-family loss before publishing. Failed calibration leaves
existing output untouched. The `policy_cost_version=1` block is bound to the
build, CPU, thread count, and eligible cache capacity. It changes cost estimates,
leaving kernel thresholds at their defaults. Calibration is optional. The
ordinary autotuner modes and existing supported profile syntax remain available.
`--policy-only` is available in the default build.

`scripts/run_generalization_campaign.py` records paired direct-lifecycle runs
against the original and optional previous-stage controls. It verifies changed
solves with the benchmark harness and records failures and timeouts. Timing uses
the workspace lock shared with the Xyce campaign and rejects active compiler,
test, or solver processes. Matrix, image, executable, harness, and explicitly supplied tuning-profile hashes
identify the tested inputs. Calibration, diagnostics, tests, and instrumented runs are not
qualification timings. `scripts/qualify_spice_variant.py`
selects the prototype consistently across its native, portable, ASan/UBSan,
and scoped TSan builds.
Those build/test qualifications are separate from performance promotion.

`scripts/generalization_gate.py` rejects missing or duplicate observations and
candidate-only failures. The current acceptance policy, `generalization-v2`,
implements the November 7 requirement: individual cases may regress provided
KLS preserves each existing win against CKTSO and SubtreeLU. For each matrix,
configuration, and KLS control, the gate computes the paired median of
control/peer lifecycle ratios. If it is below 1, the candidate/peer paired
median must be at most 1. This protects each peer independently, including
cases where the control beat only one of them. Equality is allowed for the
candidate. The combined review protects the union of prior wins identified in
either campaign: the candidate must preserve each such win in both campaigns.
A missing comparison cannot confirm a protected win. This prevents a changing
control/peer ordering between runs from hiding a lost win. Cases already slower
than a peer in both campaigns need not acquire a new win.
Geometric means must still be at most 1.00 for every configuration/KLS control.
The old 1.02 case cap is reported for comparison and does not veto this policy.

Both peer executables and their external library/license inputs are frozen
alongside the KLS executables. Missing or invalid peer measurements on an
otherwise complete KLS comparison prevent acceptance. Shared invalid KLS rows
and invalid references remain disclosed exclusions; they prevent a
complete-corpus claim. A separate, opposite-order campaign with matching inputs
is still required. The Xyce runner and its independent auditor implement the
same peer-win rule across all 28 settings for both linear and wall time, while
requiring aggregate nonregression in both metrics. Use
`--acceptance-policy generalization-v2 --peers` for that runner. The original
`generalization-v1` gate remains available to reproduce historical decisions.

The matrix runner defaults to the new policy and requires `--ck`, `--st`, and
`--peer-library` inputs in addition to the KLS controls. The unchanged revision I
solver binaries can be evaluated under this new policy without rebuilding KLS.
The original screen lacks peer observations, so its results cannot establish
preservation of peer wins or be relabeled as a passing v2 campaign.

Development has already rejected several prototypes. The AMD-only revision C
passed the targeted H1/H10/H100 screens, but its broader 48-matrix H1 screen
regressed on discarded block work and BTF representation choice. Revision D
removed the large block-work regressions while retaining other losses. Revision
E improved the aggregate over 46 valid matrices by about 10%, but its unrestricted
one-block matching removal regressed badly. Revision F's exact-symmetry variant
improved the aggregate by about 13% and still had per-case losses above 2%, so it
also failed the gate. Both one-block shortcuts have been discarded. Aggregate
improvement does not qualify these revisions. Revision G also exposed a
one-sided pending block placeholder changing the orientation decision; its
expanded block domain has been disabled for the AMD stage. Two matrices in that screen were
invalid for both controls; those are not successful coverage.

The latest fixed revision, **I**, retained current block admission and removed
provisional budget credit. Its source digest is
`ae23665ac2132d33e3355faae6425b608d8e47a10ae47290c5896f8ac9f32bf1`.
The complete receipt is
`build/generalization-20261106/stage1i-checkpoint.json`; the matching-input
combined gate is
`build/generalization-20261106/stage1i-paper-small-initial/confirmed.json`.
The runtime source matched the frozen snapshot after both campaigns. Only this
qualification documentation and inventory were updated afterward.

Both native build modes passed all 74 tests. The prototype also passed five
scoped portable tests and five scoped ASan/UBSan tests. The shared-budget/model
helpers passed their scoped TSan check; that is not a full-solver TSan result.
All 12 profile-provenance and matrix-gate Python tests passed after the final
provenance update. The public header remains byte-identical to the control.

The fixed 48-matrix H1 development screen used eight threads on each cache
domain, eight balanced pairs, and a separate opposite-order confirmation.
Each campaign recorded 2,304 observations across the original control,
prototype, and current default build. The `previous` arm means the current
default build, not a previously accepted generalization stage. Lifecycle
geometric-mean ratios were:

| Campaign | Original, 96 MiB | Original, 32 MiB | Default, 96 MiB | Default, 32 MiB |
| --- | ---: | ---: | ---: | ---: |
| Initial | 0.8937 | 0.8972 | 0.8891 | 0.8937 |
| Confirmation | 0.8934 | 0.8937 | 0.8872 | 0.8912 |

**Both campaigns failed the original v1 qualification.** Eleven cell/reference comparisons
exceeded 1.02 initially and thirteen did so in confirmation. The largest
recorded ratios were 1.1304 initially and 1.0353 in confirmation. No
candidate-only correctness failures occurred. `bips07_1998` and `fpga_dcop_01`
were invalid for all three arms, leaving 46 valid matrices in each aggregate;
they are disclosed exclusions, not complete coverage. The two optional
synthetic calibration attempts also failed their performance gate and
published no profile. None of these results qualifies the prototype or proves
performance preservation. These are historical v1 decisions; the amended peer
criterion requires new measurements and does not change those recorded results.

The November 7 peer checkpoint reuses the unchanged revision I KLS binaries,
with freshly built CKTSO and SubtreeLU harnesses. Two independent campaigns
each retained 3,840 observations across five arms, using the same 48 matrices,
eight threads on both cache domains, and eight balanced pairs per setting.
The receipt is `build/generalization-20261107/stage1i-peer-checkpoint.json`.
The final combined review is
`build/generalization-20261107/stage1i-peer-initial/confirmed.json`.

| Campaign | Original, 96 MiB | Original, 32 MiB | Default, 96 MiB | Default, 32 MiB |
| --- | ---: | ---: | ---: | ---: |
| Initial | 0.8876 | 0.8863 | 0.8901 | 0.8886 |
| Confirmation | 0.8917 | 0.8950 | 0.8896 | 0.8936 |

Using each campaign's own prior-winning set, all 317 initial and 320
confirmation wins were preserved. Seven initial and twelve confirmation
cell/reference comparisons exceeded the old 1.02 cap; that cap no longer vetoes
acceptance. However, the combined review flags **`gemat12` versus SubtreeLU on
the 96 MiB domain**. In confirmation, the original and current default KLS
controls beat SubtreeLU (ratios 0.9649 and 0.9722), making both protected wins.
In the initial run, the candidate/SubtreeLU ratio was 1.0519, so the union of
prior wins is not preserved in both campaigns. The two control comparisons
refer to the same matrix/configuration/peer issue, not two independent cases.
The initial controls also trailed SubtreeLU (1.0184 and 1.1118); this is an
unstable peer ordering, not evidence that the candidate alone caused the flip.

Phase inspection of all sixteen `gemat12` observations on this cache domain
localizes most of the changing peer comparison to SubtreeLU analysis:

| Median phase/lifecycle | Initial | Confirmation |
| --- | ---: | ---: |
| Prototype KLS lifecycle | 12.522 ms | 12.636 ms |
| Prototype KLS analysis | 3.767 ms | 3.782 ms |
| Prototype KLS initial factor | 7.625 ms | 7.550 ms |
| SubtreeLU lifecycle | 12.166 ms | 12.858 ms |
| SubtreeLU analysis | 11.497 ms | 12.195 ms |

SubtreeLU analysis spans 9.044–12.204 ms in the initial campaign and
11.921–12.861 ms in confirmation. KLS consistently selects normal AMD with
27,215 L entries, 27,367 U entries, and 355,836 factor flops. SubtreeLU's factor
structure is also constant. This does not show an ordering-policy flip in KLS.
The peer analysis timing spread can change the winner at this small margin;
the underlying cause (for example dispatch, scheduling, or allocation cost)
has not been established. A focused fixed diagnostic campaign with phase CPU
and elapsed time is needed before attributing it to a solver defect or
changing production policy. The existing qualification flag remains unresolved;
no matrix-specific exception or selected retry has been applied.

SubtreeLU additionally failed initial factorization with status -3 on
`bips98_1142`, `bips98_1450`, and `bips98_606` in both campaigns. These leave
twelve matrix/configuration/control/peer comparisons incomplete per campaign.
They are not counted as peer wins or candidate failures. The two shared-invalid
matrices from the original screen remain excluded, leaving 46 valid KLS
matrices and 43 matrices valid for all five arms. There were no candidate-only
correctness failures. The overall v2 gate remains unpassed because of the
cross-campaign peer issue and missing peer evidence, rather than the old case
cap. Both measurement harnesses, raw observations, and the original per-run-only
review (`confirmed-per-run.json`) are retained. The final stricter reviewer is
separately hashed in the receipt. All 22 matrix-gate tests and 12 SPICE campaign
tests pass, including independent reconstruction of the SPICE confirmation rule.
Full Xyce performance qualification under v2 has not been run.

The final qualification still requires the complete 110-matrix paper corpus,
prior development/regression cases, sensitivity configurations, all 28 Xyce
settings, and a fresh independently selected unseen corpus. As checked on
November 6, 2026, [SuiteSparse reports that matrix downloads are offline and no
mirror is known](https://sparse.tamu.edu/). The missing development/holdout
matrices are not installed locally. Cached paper matrices and generated cases
permit further development, but cannot replace that evidence or justify
a claim of complete performance qualification. Later matching, scale/pivot, and compact
executor stages remain pending.
