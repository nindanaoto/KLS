import copy
import unittest
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from auto_btf_merge_gate import CASES, DOMAINS, evaluate


class MergeGateTests(unittest.TestCase):
    def setUp(self):
        self.cells = [dict(domain=d, case=c, threads=t, peer=p, complete=True,
                           median=.97 if 'invert50' in c and t == 1 and p == 'baseline' else 1.,
                           wall_median=1.)
                      for d in DOMAINS for c in CASES for t in (1, 8)
                      for p in ('baseline', 'ck', 'st')]

    def run_gate(self):
        return evaluate(self.cells, all_qualified=True)

    def test_existing_peer_gaps_do_not_block(self):
        self.assertTrue(self.run_gate()['passed'])

    def test_regression_needs_both_peer_wins(self):
        cell = next(c for c in self.cells if c['case'].endswith('/dac.cir') and c['peer'] == 'baseline')
        cell['median'] = 1.1
        self.assertFalse(self.run_gate()['passed'])
        for peer in self.cells:
            if (peer['domain'], peer['case'], peer['threads']) == (cell['domain'], cell['case'], cell['threads']) and peer['peer'] != 'baseline':
                peer['median'] = .98
        self.assertTrue(self.run_gate()['passed'])

    def test_wall_regression_also_requires_exception(self):
        next(c for c in self.cells if c['peer'] == 'baseline')['wall_median'] = 1.03
        self.assertFalse(self.run_gate()['passed'])

    def test_noise_boundary(self):
        next(c for c in self.cells if c['case'].endswith('/dac.cir') and c['peer'] == 'baseline')['median'] = 1.02
        self.assertTrue(self.run_gate()['passed'])

    def test_target_cannot_use_peer_exception(self):
        next(c for c in self.cells if 'invert50' in c['case'] and c['threads'] == 1 and c['peer'] == 'baseline')['median'] = .981
        self.assertFalse(self.run_gate()['passed'])

    def test_incomplete_duplicate_invalid_and_unqualified_rejected(self):
        original = copy.deepcopy(self.cells)
        for mutate in (lambda: self.cells.pop(), lambda: self.cells.append(copy.deepcopy(self.cells[0])),
                       lambda: self.cells[0].update(median=float('nan')),
                       lambda: self.cells[0].update(complete=False)):
            self.cells = copy.deepcopy(original)
            mutate()
            self.assertFalse(self.run_gate()['passed'])
        self.assertFalse(evaluate(original, all_qualified=False)['passed'])


if __name__ == '__main__':
    unittest.main()
