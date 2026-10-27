# Complete SPICE qualification

`scripts/run_spice_campaign.py` compares complete Xyce executions on the seven
existing SPICE circuits, at 1 and 8 solver threads, on both Ryzen cache domains.
The primary measurement is the **last** `Total Linear Solution Time` in each
log. Process elapsed time is a separate regression guard. Instrumented runs
and smoke tests cannot pass the performance gates.

## Build and correctness

Use `scripts/qualify_spice_variant.py --help` to freeze each candidate's actual
source bytes and build native, portable, ASan/UBSan, and scoped TSan variants.
It verifies the four external dependency gitlinks against clean local sources,
records commands and exit codes, and retains failed logs. Scoped TSan covers
public solver ownership, the pthread certificate team, and the KLU allocator
boundary. Its clean/known-race controls must both behave as expected. GCC
function entry/exit hooks are disabled for the target-clone resolver startup
issue; memory instrumentation remains enabled. Direct OpenMP reductions are
outside this scope: the pinned uninstrumented libgomp reports a reduction race
on unchanged main as well as candidates. This is not full Xyce race coverage.
A failed sanitizer run
requires investigation; it is not a performance qualification.

`scripts/build_spice_links.py` rebuilds the KLS and third-party Xyce adapters
against the selected source. Build the control first; pass its `libxyce.a` as
`--common-archive` and its KLS build as `--dependency-build` for candidate links.
This keeps the adapters and dependencies identical across solver comparisons.
The public header must match. The link command and all explicit link inputs are
hashed. `bench/spice_link.json` specifies the pinned GCC 15 container ABI.
`--coarse-profile` is solely for the isolated source produced by
`prepare_coarse_profile.py`; its timings are diagnostic.

Before timing, run `--diagnostic --passes 1` on the complete grid. The runner
requires successful termination, finite waveforms, matching accepted/rejected
step counts, and the upstream Xyce verifier at `rtol=1e-8, atol=1e-10`.
Following the existing complete-circuit benchmark protocol, strict pointwise
and grid differences are retained separately as diagnostics: adaptive grids
can differ even for unchanged main. The upstream verifier performs its usual
transient interpolation. The solver's default componentwise
accuracy contract remains unchanged. `--smoke --passes 1 --peers` is a shorter
mux8 check of binaries, licenses, adapters, and verification plumbing.

## Timing

Supply absolute paths for `--workspace`, `--regression`, `--image`, `--license`,
`--baseline`, `--candidate`, `--output`, `--provenance`, and each
`--peer-library`. The regression tree must contain `Netlists` and
`TestScripts/xyce_verify.pl`. Use a fresh output directory for every campaign.

Run four passes with `--order initial`; repeat independently using
`--order confirmation`. Each pair-only campaign has 238 observations, including
14 fresh KLU references. Adding `--peers` compares CKTSO and SubtreeLU and
produces 462 observations. Launch positions are balanced. Benchmark sessions
are serialized, and ordinary campaigns reject active compiler/test/solver
processes before starting. Do not run builds or modify pinned inputs during
timing. Failed runs remain in the ledger; there is no automatic retry.

An incremental candidate must complete all 28 settings, improve at least one
linear-time setting by 2%, and regress neither linear time nor application
elapsed time by more than 2% anywhere. Final success additionally requires a
2% linear-time advantage in **all 56 peer comparisons**. Ratios use medians of
matched passes. Both gates require a separate confirmation campaign.

Recheck retained evidence with:

```sh
python3 scripts/run_spice_campaign.py --review /absolute/campaign
python3 scripts/run_spice_campaign.py --confirm /absolute/initial /absolute/confirmation
```

Review rechecks hashes, logs, waveforms, grid completeness, and launch order.
Confirmation requires matching binaries and pinned inputs and opposite orders.
`all56_peer_wins` in a single `summary.json` is provisional; only the combined
`confirmation.json` establishes the replicated gate. This mechanical review
does not substitute for human review of numerical changes.
