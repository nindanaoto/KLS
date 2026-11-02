"""Focused AUTO merge gate; consumes independently qualified campaign cells."""
CASES = frozenset(('PDE_2D_BJT/invertbjt.cir', 'MOS13_IC/invert50_mos1.cir',
    'CircuitSim90/MOS2_LARGE/mem_plus.cir', 'CircuitSim90/MOS2/mux8.cir',
    'CircuitSim90/MOS2_LARGE/dac.cir', 'NEURON/rallpack3.cir',
    'CircuitSim90/MOS2_LARGE/ram2k.cir'))
DOMAINS = ('llc96', 'llc32')


def evaluate(cells, *, all_qualified):
    """Require target benefit; permit material main regressions only with peer wins.

    Callers must establish the complete original grid, pinned inputs, four balanced
    passes, and qualified observations before supplying all_qualified=True.
    """
    import math
    expected = {(d, c, t, p) for d in DOMAINS for c in CASES
                for t in (1, 8) for p in ('baseline', 'ck', 'st')}
    indexed = {}
    valid = all_qualified is True
    for cell in cells:
        key = (cell['domain'], cell['case'], cell['threads'], cell['peer'])
        if key in indexed:
            valid = False
        indexed[key] = cell
        valid &= bool(cell.get('complete'))
        for name in ('median', 'wall_median'):
            value = cell.get(name)
            valid &= type(value) in (int, float) and math.isfinite(value) and value > 0
    valid &= set(indexed) == expected
    if not valid:
        return dict(passed=False, complete=False, targets=[], regressions=[])
    targets = [indexed[d, 'MOS13_IC/invert50_mos1.cir', 1, 'baseline'] for d in DOMAINS]
    regressions = []
    for d, c, t, p in sorted(expected):
        if p != 'baseline':
            continue
        main = indexed[d, c, t, p]
        if main['median'] > 1.02 or main['wall_median'] > 1.02:
            peer_wins = all(indexed[d, c, t, peer]['median'] <= .98 for peer in ('ck', 'st'))
            regressions.append(dict(domain=d, case=c, threads=t, peer_exception=peer_wins))
    return dict(passed=all(c['median'] <= .98 for c in targets) and
                all(c['peer_exception'] for c in regressions), complete=True,
                targets=[{k:c[k] for k in ('domain','case','threads','median','wall_median')} for c in targets], regressions=regressions)
