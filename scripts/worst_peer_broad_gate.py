"""Peer-aware standalone guard over independently audited paired evidence.

The caller must bind the original 110-matrix campaign and fresh peer samples
to identical candidate binaries and inputs. This reducer does not authenticate
receipts. Original regressions are never erased by a later baseline retest.
"""
import math
import statistics


def positive(value):
    return type(value) in (int, float) and math.isfinite(value) and value > 0


def evaluate(groups, peer_groups, *, expected_cells, coverage_qualified,
             peer_evidence_qualified):
    invalid = dict(passed=False, complete=False, failures=[], peer_exceptions=[])
    if coverage_qualified is not True:
        return invalid
    expected = set(expected_cells)
    if not expected or len(expected) != len(expected_cells):
        return invalid
    indexed = {}
    for group in groups:
        key = (group.get('config'), group.get('matrix'))
        values = group.get('ratios')
        count = group.get('valid_pairs')
        if (key not in expected or key in indexed or type(count) is not int
                or not 0 <= count <= 4 or not isinstance(values, list)
                or len(values) != count or not all(positive(v) for v in values)):
            return invalid
        median = statistics.median(values) if count == 4 else None
        if group.get('median') != median:
            return invalid
        indexed[key] = median
    if set(indexed) != expected:
        return invalid
    regressed = {key for key, median in indexed.items()
                 if median is not None and median > 1.02}
    peers = {}
    for group in peer_groups:
        key = (group.get('config'), group.get('matrix'))
        ratios = group.get('ratios')
        if (key not in regressed or key in peers or not isinstance(ratios, dict)
                or set(ratios) != {'ck', 'st'}):
            return invalid
        if any(not isinstance(values, list) or len(values) != 4
               or not all(positive(v) for v in values) for values in ratios.values()):
            return invalid
        peers[key] = {name: statistics.median(values) for name, values in ratios.items()}
    if set(peers) != regressed or (regressed and peer_evidence_qualified is not True):
        return invalid
    failures, exceptions = [], []
    for key in sorted(regressed):
        record = dict(config=key[0], matrix=key[1], baseline_median=indexed[key],
                      peer_medians=peers[key])
        (exceptions if all(v <= .98 for v in peers[key].values()) else failures).append(record)
    return dict(passed=not failures, complete=True, failures=failures,
                peer_exceptions=exceptions)
