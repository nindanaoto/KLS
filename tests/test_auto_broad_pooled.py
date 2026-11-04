import copy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from reassess_auto_broad_pooled import evaluate


class PooledBroadGateTests(unittest.TestCase):
    def setUp(self):
        self.groups = [dict(config='a', matrix='unchanged', median=.99),
                       dict(config='a', matrix='regressed', median=1.09),
                       dict(config='b', matrix='peer_win', median=1.03),
                       dict(config='a', matrix='shared_failure', median=None)]
        self.pooled = [dict(config='a', matrix='regressed', pairs=40, median=1.009,
                            peer_medians=dict(ck=1.04, st=.7)),
                       dict(config='b', matrix='peer_win', pairs=40, median=1.022,
                            peer_medians=dict(ck=.8, st=.9))]

    def test_approved_complete_evidence(self):
        self.assertTrue(evaluate(self.groups, self.pooled)['passed'])

    def test_unresolved_regression_rejected(self):
        self.pooled[0]['median'] = 1.020001
        self.assertFalse(evaluate(self.groups, self.pooled)['passed'])

    def test_tie_is_not_a_peer_win(self):
        self.pooled[1]['peer_medians']['ck'] = 1.
        self.assertFalse(evaluate(self.groups, self.pooled)['passed'])

    def test_all_original_failures_required(self):
        self.assertFalse(evaluate(self.groups, self.pooled[1:])['passed'])
        self.groups.append(dict(config='c', matrix='missed', median=1.04))
        self.assertFalse(evaluate(self.groups, self.pooled)['passed'])

    def test_duplicate_or_unrequested_cell_rejected(self):
        self.assertFalse(evaluate(self.groups, self.pooled + [copy.deepcopy(self.pooled[0])])['passed'])
        self.pooled.append(dict(config='x', matrix='extra', pairs=40, median=.9,
                                peer_medians=dict(ck=.8, st=.7)))
        self.assertFalse(evaluate(self.groups, self.pooled)['passed'])

    def test_missing_pairs_or_peers_rejected(self):
        original = copy.deepcopy(self.pooled)
        self.pooled[0]['pairs'] = 39
        self.assertFalse(evaluate(self.groups, self.pooled)['passed'])
        self.pooled = original
        del self.pooled[0]['peer_medians']['ck']
        self.assertFalse(evaluate(self.groups, self.pooled)['passed'])

    def test_nonfinite_evidence_rejected(self):
        self.pooled[0]['median'] = float('nan')
        self.assertFalse(evaluate(self.groups, self.pooled)['passed'])


if __name__ == '__main__':
    unittest.main()
