# Cancellation checks and numeric-generation certificates

Two general correctness fixes in `src/kls.c` resolve the discovered SPICE
simulation failures without circuit-specific dispatch or relaxed tolerances:

- If the ordinary relative-L2 check fails for an internally transposed
  solve, recompute the residual with long-double accumulation before
  accepting/rejecting the answer or refining it. The 5e-9 limit is unchanged.
  This cold check prevents cancellation in the residual computation from
  falsely rejecting a valid solution.
- Refactorization refreshes the condition estimate and every retained matrix
  snapshot. Per-answer/RHS validity is not reused for changed values, and
  failed refinement discards the pivot-family policy as well. A settled
  policy that the diagonal-ratio heuristic was pessimistic is retained for
  the same fixed-pivot family; it is not an answer or RHS certificate.

`tests/kls_cancellation.c` checks a nearly floating three-node resistor
chain, two power-of-two rescalings, a zero-RHS certificate followed by a
changed matrix, a second RHS in each epoch, and the refreshed condition
estimate. The test fails against
the pre-fix implementation and passes after the fix.

## Validation of the original two-circuit fix (d3400b5)

GCC 14 native Release, same machine and default solver settings:

| Circuit | KLS-1 and KLS-8 | Maximum output difference from KLU |
| --- | --- | --- |
| MOS13_IC/invert50_mos1 | 3/3 complete each | 8.0e-10 V |
| NEURON/passivePatch-level6 | 3/3 complete each | 3.0e-10 V |

Every printed waveform value passes the existing comparison tolerance
`1e-10 + 1e-8 * max(abs(reference), abs(test))`. Xyce's own waveform verifier
also passes; the neuron additionally passes its analytical RC reference.

The inverter reports one rejected intermediate DC linear solve, **as does
KLU**, then converges and takes the same 9,875 successful steps and 251
rejected time steps. The fix does not hide genuinely inaccurate linear
answers or weaken their acceptance criterion. The neuron reports no failed
linear solves and completes 68 successful steps.

All eight standalone CTest tests and Xyce's adapter test pass. Full paper
performance benchmarking has not been rerun. Raw observations and binary/
source provenance are in the sibling paper repository under
`evaluation/spice/results/failure-fix-final-20260914/`.

An experimental opposite-orientation recovery was tested but not retained:
the two fixes above suffice for complete, validated simulations. There is
no new fallback solver, backend substitution or matrix-size exception.

## BJT regression follow-up

Commit d3400b5 cleared the retained pivot-family policy as well as the
per-answer validity state. Repeating the pessimistic low-rcond probe on every
changed numeric made `PDE_2D_BJT/invertbjt` fail DC convergence. The follow-up
separates those lifetimes: it preserves the established fixed-pivot policy,
refreshes the diagonal ratio, and copies every retained residual snapshot
from the current values even when checking is dormant. A failed check still
forces fresh classification; a full factorization resets the family policy.

No accuracy threshold is changed. Experimental public-limit certification
and transposed GMRES extensions were not retained. The existing cancellation
fix and failed-refinement invalidation remain in place.

The sibling paper repository now contains a reproducible four-circuit check:

```
python3 evaluation/spice/validate_spice_correctness.py --results /path/to/new/results
```

It runs fresh KLU references plus three KLS-1 and KLS-8 repetitions of the
BJT inverter, MOS1 inverter, passive neuron, and Rallpack3. Every printed
waveform value must match at the existing tolerance. The regression circuits
are test inputs only, not solver configuration or dispatch rules.

Follow-up results (`results/bjt-policy-validation-20260914/` in the paper
repository): all 24 waveform comparisons pass. Median BJT linear time is
6.995 s with one thread and 2.114 s with eight threads; both modes have zero
failed linear solves and exactly match the printed KLU waveform. Rallpack3
also passes in both modes. The MOS1 inverter matches within 1e-9 V and
converges in every run; it reports one rejected DC solve with one thread
and three with eight threads, versus one for KLU. The passive neuron has
zero rejected solves and matches within 3e-10 V. Intermediate failures are
recorded rather than hidden. No full paper performance rerun was performed.
