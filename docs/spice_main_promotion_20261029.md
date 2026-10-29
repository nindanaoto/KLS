# Confirmed incremental SPICE improvement

Merged qualified numerical source `4c9c2873127159e5c8fbf2daa2df747c12858676` from main `56abf844fb07931370ef386ff724124dfa5d965c`.

Both SPICE orders passed: 924 correct observations, all 28 main comparisons improved in each order, and independent audits matched every cell. The fresh 1,760-observation broad comparison and independent review found no candidate-only correctness loss or greater-than-2% median lifecycle regression in any fully qualified matrix/domain group.

There are 40 shared failing pairs. These remain unresolved and excluded from performance aggregation; no complete-corpus speed aggregate is claimed. Standalone broad process-wall differences are retained as diagnostics in the evidence receipt. The SPICE process-wall guard passed both complete Xyce timing orders. This is an incremental main improvement, not a claim of beating both peers in every SPICE setting.

Evidence archive SHA-256: `2c60fe06c7f55297b11e51dda8ab6f8f335ed418067e76f2f3470884387037fb`. See `build/spice-evidence-20261027/main-improvement-broad-v2.json` for every qualified group and shared failure.
