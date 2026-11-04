"""Versioned worst-peer-v1 gate for complete paired SPICE campaigns.

Final success depends on peer margins, not improvement against an older KLS.
Confirmation and provenance remain the campaign runner's responsibility.
"""
import math

from auto_btf_merge_gate import CASES, DOMAINS


def evaluate(cells, *, all_qualified):
    expected = {(d, c, t, p) for d in DOMAINS for c in CASES
                for t in (1, 8) for p in ('baseline', 'ck', 'st')}
    invalid = dict(passed=False, complete=False, incremental_passed=False,
                   worst_ratio=None, peer_wins=0, regressions=[])
    indexed = {}
    if all_qualified is not True:
        return invalid
    for cell in cells:
        key = tuple(cell.get(k) for k in ('domain', 'case', 'threads', 'peer'))
        fields = ['median', 'wall_median']
        if key[-1] != 'baseline':
            fields.append('baseline_peer_median')
        if (key not in expected or key in indexed or cell.get('complete') is not True
                or any(type(cell.get(f)) not in (int, float)
                       or not math.isfinite(cell[f]) or cell[f] <= 0 for f in fields)):
            return invalid
        indexed[key] = cell
    if set(indexed) != expected:
        return invalid
    regressions = []
    for d, c, t, p in sorted(expected):
        if p != 'baseline':
            continue
        baseline = indexed[d, c, t, p]
        peers = [indexed[d, c, t, peer] for peer in ('ck', 'st')]
        if all(peer['median'] <= .98 for peer in peers):
            continue
        lost = [peer['peer'] for peer in peers
                if peer['baseline_peer_median'] <= 1 and peer['median'] > 1]
        if baseline['median'] > 1.02 or baseline['wall_median'] > 1.02 or lost:
            regressions.append(dict(domain=d, case=c, threads=t, lost_peer_wins=lost))
    peers = [v['median'] for k, v in indexed.items() if k[-1] != 'baseline']
    improved = any(v['median'] <= .98 for k, v in indexed.items() if k[-1] == 'baseline')
    return dict(passed=max(peers) <= .98, complete=True,
                incremental_passed=not regressions and improved,
                worst_ratio=max(peers), peer_wins=sum(v <= .98 for v in peers),
                regressions=regressions)
