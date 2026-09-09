# Remove unreachable predicted-factor thresholds

Baseline: `99cd8d4`.

Remove the unreachable size, sparsity, BTF and flop threshold checks after
the unconditional return in kls_predicted_probe_first_threshold_shape.
Flatten the redundant scope around that return while retaining the comment
explaining the live generic lifecycle candidate. Net reduction: 19 C-source
lines. The generic predicate and mandatory residual acceptance are unchanged.

Warning-free Release and ASan builds pass CTest 6/6 each (3.69 and 15.42
seconds). The rebuilt Release executable's .text is byte-identical to the
frozen preceding executable. More strongly, complete executables compare
byte-identical after objcopy --strip-debug --remove-section=.note.gnu.build-id
is applied to each. Thus this local Release build has identical runtime
code/data/layout; only debug/build-identity content differs. No timing rerun
is needed to attribute performance for this dead-code-only chunk. This is
not a claim about arbitrary other compiler builds.

Original SHA-256: 4ffbe37c2da4c20524cf8932009b00c87a2dfa0e1977d2430c0c376e782605ec.
Rebuilt SHA-256: 50c59d220e2179256a4f0473697a181d53a0b5832908ffdf449ff9a6acfe7173.
Artifacts: `build/prep-trim-repair-dBcoha/dead-threshold-*` (frozen baseline,
before/after .text sections and stripped executables).

The broader cleanup goal remains active.
