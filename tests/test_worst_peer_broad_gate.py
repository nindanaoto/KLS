import copy
import math
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from worst_peer_broad_gate import evaluate


class BroadPeerGuardTests(unittest.TestCase):
    def setUp(self):
        self.groups = [dict(config='llc96', matrix='slower', valid_pairs=4,
                            ratios=[1.3] * 4, median=1.3),
                       dict(config='llc32', matrix='old_peer_loss', valid_pairs=4,
                            ratios=[1.02] * 4, median=1.02),
                       dict(config='llc96', matrix='shared_failure', valid_pairs=0,
                            ratios=[], median=None)]
        self.peers = [dict(config='llc96', matrix='slower', ratios=dict(ck=[.98]*4, st=[.9]*4))]
        self.expected = [(g['config'], g['matrix']) for g in self.groups]

    def decision(self, **kw):
        args = dict(expected_cells=self.expected, coverage_qualified=True,
                    peer_evidence_qualified=True)
        args.update(kw)
        return evaluate(self.groups, self.peers, **args)

    def test_large_main_regression_permitted_only_with_two_margin_wins(self):
        self.assertTrue(self.decision()['passed'])
        for peer in ('ck', 'st'):
            original = copy.deepcopy(self.peers)
            self.peers[0]['ratios'][peer] = [math.nextafter(.98, 1)]*4
            result = self.decision()
            self.assertTrue(result['complete'])
            self.assertFalse(result['passed'])
            self.peers = original

    def test_existing_peer_losses_need_not_be_fixed(self):
        self.groups = self.groups[1:]
        self.expected = self.expected[1:]
        self.peers = []
        self.assertTrue(self.decision(peer_evidence_qualified=False)['passed'])

    def test_coverage_and_fresh_peer_evidence_are_mandatory(self):
        self.assertFalse(self.decision(coverage_qualified=False)['passed'])
        self.assertFalse(self.decision(peer_evidence_qualified=False)['passed'])
        self.peers = []
        self.assertFalse(self.decision()['complete'])

    def test_missing_duplicate_or_extra_cells_fail_closed(self):
        original = copy.deepcopy(self.groups)
        for groups in (original[:-1], original + original[:1],
                       original + [dict(config='extra', matrix='extra')]):
            self.groups = groups
            self.assertFalse(self.decision()['complete'])
        self.groups = original
        self.peers += copy.deepcopy(self.peers)
        self.assertFalse(self.decision()['complete'])

    def test_invalid_measurements_and_wrong_pair_counts_are_rejected(self):
        for value in (0, -1, True, float('nan'), float('inf')):
            self.peers[0]['ratios']['ck'] = [value]*4
            self.assertFalse(self.decision()['complete'])
        self.peers[0]['ratios']['ck'] = [.9]*3
        self.assertFalse(self.decision()['complete'])

    def test_original_paired_median_cannot_be_replaced_by_later_results(self):
        self.groups[0]['median'] = 1.0
        self.assertFalse(self.decision()['complete'])
        self.groups[0]['median'] = 1.3
        self.peers[0]['ratios']['ck'] = [.8, .9, 1.05, 1.2]
        self.assertTrue(self.decision()['passed'])  # paired median .975
        self.peers[0]['ratios']['st'] = [.99]*4
        self.assertFalse(self.decision()['passed'])


if __name__ == '__main__':
    unittest.main()
