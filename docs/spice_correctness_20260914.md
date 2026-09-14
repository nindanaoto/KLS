# Cancellation checks and numeric-generation certificates

Two general correctness fixes in `src/kls.c` resolve the discovered SPICE
simulation failures without circuit-specific dispatch or relaxed tolerances:

- If the ordinary relative-L2 check fails for an internally transposed
  solve, recompute the residual with long-double accumulation before
  accepting/rejecting the answer or refining it. The 5e-9 limit is unchanged.
  This cold check prevents cancellation in the residual computation from
  falsely rejecting a valid solution.
- Refactorization invalidates a previous condition-based certificate or
  failed-refinement state. Classification then refreshes the condition
  estimate and, when required, the current matrix snapshot. A certificate
  for one set of values must not certify a later numeric generation.

`tests/kls_cancellation.c` checks a nearly floating three-node resistor
chain, two power-of-two rescalings, a zero-RHS certificate followed by a
changed matrix, and the refreshed condition estimate. The test fails against
the pre-fix implementation and passes after the fix.

## Validation

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
