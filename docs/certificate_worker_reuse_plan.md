# Certificate execution on existing numerical workers

Unqualified prototype implemented; do not conflate it with the original-operator cache trial.
The fresh profiles in `/tmp/kls-refinement-first/qualified-profiles-resumed/`
measured 9.37 s of certification in ram2k-8, 1.16 s in mem_plus-8 and 0.71 s
in BJT-8. Existing serial row checking prevents the competing OpenMP-team
regression, but these costs still leave substantial potential benefit.

## Inspected integration boundaries

- The ordinary refactor pool owns `thread_count` pthread workers; the caller
  waits rather than acting as worker zero. Its generation release/acquire and
  `active_workers` completion counter already provide a synchronous dispatch.
- The elimination-graph pool owns `thread_count - 1` pthread workers plus the
  caller. `kls_egraph_pool_dispatch_and_spin_wait` publishes a generation,
  executes worker zero, and joins each worker's completed generation.
- An elimination-graph certificate job must intercept
  `kls_egraph_refactor_pool_run_worker`, not just
  `kls_egraph_refactor_worker_run`: the outer function dispatches other jobs and
  can execute a post-refactor permutation barrier. Accidentally entering that
  barrier during a certificate-only job would be incorrect.
- The ordinary pool's certificate branch must bypass numeric block traversal
  but reach the existing completion notification. It must not mark numeric
  workspaces dirty or change numeric worker failure/telemetry fields.

## Implementation and acceptance requirements

1. Add an internal synchronous row-executor hook to the componentwise checker.
   Keep the existing per-row term order, outward bound, worker floating-point
   environment check, and serial fallback. Each worker writes its own result;
   the caller combines validity and maximum bounds only after all workers join.
2. Publish a short-lived job descriptor using the existing pool lock and
   generation handshake. Never create a second team or change process-wide
   OpenMP policy. Preserve existing hot solver-field offsets. Do not overwrite
   any numerical job state that a later solve/refactor expects to retain.
3. Acquire only an idle, complete pool. Return to serial checking when no safe
   executor is available; do not guess thread ownership from configured counts.
   Avoid using an unrelated retained pool when multiple numerical teams exist.
4. Retain bounded serial/parallel comparison using the existing probe count.
   Invalidate timing selection when the executor changes. No circuit names,
   matrix-size exceptions, new acceptance limits, or cached solve answers.
5. Test both pool kinds, differing actual team sizes, repeated solves/refactors,
   normal/transpose, multiple RHS, worker FP rejection, dispatch failure,
   independent simultaneous instances, and destruction after a dispatched job.
6. Run this independently against a frozen qualified baseline. Require the
   full current performance and exact accuracy gates before promotion. The
   original-operator cache must first be evaluated independently so its effect
   is not confused with worker reuse.

The objective remains fastest verified performance across the current SPICE
suite, not merely passing these incremental experiment gates. Competitor
comparisons must be refreshed after retained optimizations.

The working prototype adds a synchronous executor hook with per-worker result
slots and uses the existing completion handshakes for both pool kinds. It
declines ambiguous ownership when both numerical pools are retained. Tests for
pool dispatch, numeric continuation, concurrent owners, executor failure and
unsupported worker rounding are implemented but have not run yet. Verification
is queued by `/tmp/kls-refinement-first/test_worker_reuse_after_cache.py`; it
watches the live frozen campaign processes and checks prototype source hashes
before building. Evidence will be under `worker-reuse-tests/`. Do not treat
test implementation or a queued run as a passing result.

Initial execution completed: after fixing disabled assertions in the separate
cache fixture, all 40 CTests passed, as did the executor-plan ASan/UBSan fixture
and GCC/OpenMP exact-bound checks. Evidence is preserved under
`/tmp/kls-refinement-first/worker-reuse-tests-v2/`. This does not qualify worker
reuse for SPICE performance or accuracy.

Review found that the pool fixture's public solve still used a one-thread
configuration. The expanded fixture now enables the owned team for public
normal and transpose solves, tests two RHS with padding, checks padding is
untouched, and verifies that the pool generation advances during these solves.
Both pool families and team sizes remain covered. This extra coverage is
queued under `worker-reuse-tests-v3/` behind the live cache-only exact audit;
it has not yet passed. No production arithmetic changed in this coverage step.

The expanded v3 run subsequently passed all 40 CTests, the sanitizer fixture,
and exact-bound checks after cache accuracy qualification finished. Evidence:
`/tmp/kls-refinement-first/worker-reuse-tests-v3/status.json`. Worker reuse still
requires its own frozen-binary performance comparison and SPICE exact audit.

The frozen worker-reuse candidate at `/dev/shm/kls-worker-reuse-_xxfo3er/Xyce`
passed the three-block performance screen against the qualified cache-only
baseline. Eight-thread median total-time reductions were RAM2k 15.31%,
mem_plus 14.29%, DAC 11.17%, and BJT 7.70%. Rallpack3 improved 0.92%; mux8
was essentially unchanged and MOS13 was 0.14% slower. One-thread medians were
within 0.60% of baseline. Evidence is
`/tmp/kls-refinement-first/worker-reuse-performance/screen-summary.json`.
Nine fresh confirmation blocks are now running. These screening results are
not final qualification or evidence of competitor leadership. The worker
candidate still needs its own SPICE exact audit after timing confirmation.

## Confirmed performance (2026-09-29)

All 343 runs completed: seven KLU references and 336 paired comparison runs
across the screen and nine confirmation blocks. Every comparison observation
passed its waveform checks. The confirmation passed the existing incremental
benefit and 2% median-regression gates; this is not a claim of zero variation.

| Circuit | One-thread runtime change | Eight-thread runtime change |
| --- | ---: | ---: |
| BJT | +0.09% | -8.28% |
| MOS13 | +0.00% | +0.09% |
| mem_plus | -0.20% | -15.19% |
| mux8 | +0.05% | +0.28% |
| DAC | -0.01% | -9.71% |
| Rallpack3 | -0.05% | -0.20% |
| RAM2k | +0.22% | -15.63% |

Changes are median paired candidate/baseline ratios minus one, not ratios of
unpaired medians. Baseline is the independently qualified original-operator
cache candidate. Frozen worker binary SHA-256:
`e6d4ccdf6f870b19efba9352ce8648e2a2fd33cf7cef3b40e4fcb6d30089c18e`.
Evidence: `/tmp/kls-refinement-first/worker-reuse-performance/` contains
`confirmation-summary.json`, `observations.json`, and terminal `status.json`;
`/tmp/kls-refinement-first/worker-reuse/build.json` pins the build.

Fresh six-case/mode exact qualification has started under
`/tmp/kls-refinement-first/worker-reuse-accuracy/`, using
`/tmp/kls-refinement-first/qualify_worker_accuracy.py`. It retains the unchanged
exact 1/100000000 componentwise acceptance limit and explicit-rejection recovery
policy. Archive reuse is byte-identical storage reuse only, never reuse of an
earlier acceptance decision. Accuracy qualification remains pending, as do
fresh competitor comparisons and broader-corpus evaluation.

## Accuracy outcome

All six required circuit/mode audits subsequently passed: BJT and mem_plus
had respectively 1,146 and 1,201 certified successes per mode, no rejections,
and exact waveform matches. MOS13 had 30,490 certified successes and one
explicit rejection per mode, with verified recovery and passing waveforms
(maximum absolute difference about 8e-10 under the unchanged comparison rule).
Total coverage is 65,674 certified successes and two permitted rejections.

The initial campaign stopped during mem_plus capture at the 5 GiB storage
reserve, not on a numerical failure. That failed observation remains preserved.
Its interrupted raw exports were hash-verified into an archive before reclamation;
two older archives were relocated with verified hashes and original-path links.
The separate `/tmp/kls-refinement-first/worker-reuse-mem-retry/` campaign then
passed both mem_plus modes with the same frozen binary and unchanged limits.
Its `full_accuracy_suite: false` correctly describes the retry alone; six-mode
coverage combines its two audit records with the four completed records under
`worker-reuse-accuracy/`. Neither the failed original campaign nor the scoped
retry should be represented as an independently completed full-suite run.

Worker reuse now passes the focused performance and accuracy qualification.
Fresh competitor comparisons and broader-corpus evaluation remain required;
these incremental results do not establish the overall fastest-KLS objective.
