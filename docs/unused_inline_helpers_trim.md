# Remove unused inline helpers

Baseline: `14baacf`. Remove three unreferenced static inline functions,
119 source lines in total:

- `kls_scatter_subtract_skip_range`
- `kls_row_refactor_apply_segment_target_decoded`
- `kls_row_refactor_apply_compact_dense_target_decoded`

Whole-core token-reference inspection and macro-use inspection found no
callers. A repeated scan after removal finds no remaining single-reference
static `kls_` helper candidates; this is not proof that all other forms of
redundancy have been exhausted. Public interfaces and live paths are unchanged.

## Validation

After the release build completed, CTest passed 6/6 in 4.12 s. The rebuilt
ASan configuration passed 6/6 in 17.13 s. Builds completed successfully
without reported warnings; `git diff --check` passed.

Compared all 29 ELF allocated sections of the rebuilt release benchmark
against `build/prep-trim-repair-dBcoha/core-state-trim-auto-kls_bench`.
Every section has the same address and size. Contents match except for
`.note.gnu.build-id`; NOBITS sections match in address and size. Thus the
executable code and allocated runtime data are unchanged. No additional
timing campaign or alignment fix is warranted for this dead-inline removal.

Whole-file SHA-256 hashes (debug/build metadata may differ):

- Accepted: `794da8250271452834224fc23a97c40b69a24b143257d297d54728c045887ba2`
- Rebuilt: `f126d5cd4f2e84002cd0353a7f5cdff905977b586beb94af5a2c44a670ed3e03`

The preceding focused campaign remains the performance evidence for this
machine code; this does not extend its coverage to the full paper corpus
or other machines.
