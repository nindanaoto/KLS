"""User-approved gate: protect existing peer wins and bound existing losses."""
import math
from auto_btf_merge_gate import CASES, DOMAINS


def evaluate(cells, *, all_qualified):
    expected = {(d, c, t, p) for d in DOMAINS for c in CASES
                for t in (1, 8) for p in ('baseline', 'ck', 'st')}
    indexed = {}
    valid = all_qualified is True
    for cell in cells:
        key = (cell['domain'], cell['case'], cell['threads'], cell['peer'])
        valid &= key not in indexed and bool(cell.get('complete'))
        indexed[key] = cell
        fields = ('median', 'wall_median') + (() if cell['peer'] == 'baseline' else ('baseline_peer_median',))
        for field in fields:
            value = cell.get(field)
            valid &= type(value) in (int, float) and math.isfinite(value) and value > 0
    if not valid or set(indexed) != expected:
        return dict(passed=False, complete=False, targets=[], regressions=[], geometric_mean=None)
    targets = [indexed[d, 'MOS13_IC/invert50_mos1.cir', 1, 'baseline'] for d in DOMAINS]
    regressions, ratios = [], []
    for d, case, threads, peer in sorted(expected):
        if peer != 'baseline':
            continue
        main = indexed[d, case, threads, peer]
        peers = [indexed[d, case, threads, p] for p in ('ck', 'st')]
        # Exact paired main/peer medians, never a quotient of two medians.
        wins = [p for p in peers if p['baseline_peer_median'] <= 1.0]
        limit = 1.02 if wins else 1.05
        lost = [p['peer'] for p in wins if p['median'] > 1.0]
        ratios.append(main['median'])
        if main['median'] > limit or main['wall_median'] > limit or lost:
            regressions.append(dict(domain=d, case=case, threads=threads,
                limit=limit, lost_peer_wins=lost, median=main['median'], wall_median=main['wall_median']))
    geometric_mean = math.exp(math.fsum(math.log(x) for x in ratios) / len(ratios))
    return dict(passed=all(c['median'] <= .98 for c in targets) and not regressions and geometric_mean <= 1.0,
                complete=True, targets=[{k:c[k] for k in ('domain','case','threads','median','wall_median')}
                                        for c in targets], regressions=regressions, geometric_mean=geometric_mean)
