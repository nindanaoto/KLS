import copy
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from worst_peer_gate import evaluate
from audit_worst_peer_campaign import hybrid_acceptance
import run_worst_peer_campaign as campaign
import test_spice_campaign as fixture_module
from confirm_worst_peer import matching_campaigns


class WorstPeerTests(unittest.TestCase):
    def fixture(self):
        metadata, rows = fixture_module.CampaignTests().fixture()
        metadata['acceptance_policy'] = 'worst-peer-v1'
        return metadata, rows

    def check(self, metadata, rows, expected):
        result = campaign.reduce_results(metadata, rows)
        self.assertEqual(result['all56_peer_wins'], expected)
        self.assertEqual(result['acceptance_gate'], expected)
        official = result['all_qualified'] and not metadata.get('diagnostic_only', False)
        self.assertEqual(result['hybrid'], hybrid_acceptance(result['cells'], official))
        return result

    def test_peer_wins_override_main_and_wall_regression(self):
        m, rows = self.fixture()
        for row in rows:
            if row['side'] == 'baseline':
                row['linear_seconds'] = [.1]
            if row['side'] == 'candidate':
                row['wall_seconds'] = 3.
        result = self.check(m, rows, True)
        self.assertEqual(result['worst_ratio'], .9)
        self.assertEqual(result['peer_wins'], 56)
        self.assertFalse(result['incremental_gate'])

    def test_exact_margin_and_one_losing_cell(self):
        m, rows = self.fixture()
        for row in rows:
            if row['side'] == 'candidate':
                row['linear_seconds'] = [.98]
        self.check(m, rows, True)
        for row in rows:
            if row['side'] == 'candidate' and row['case'] == campaign.CASES[0]:
                row['linear_seconds'] = [.980001]
        self.check(m, rows, False)

    def test_invalid_ledger_cannot_win(self):
        m, rows = self.fixture()
        for bad in (rows[:-1], rows + [rows[-1]], rows[1:]):
            self.check(m, bad, False)
        for value in (0., -1., float('nan'), float('inf'), True):
            bad = copy.deepcopy(rows)
            next(r for r in bad if r['side'] == 'ck')['linear_seconds'] = [value]
            self.check(m, bad, False)
        m['diagnostic_only'] = True
        self.check(m, rows, False)

    def test_intermediate_regression_and_lost_win(self):
        m, rows = self.fixture()
        cells = campaign.reduce_results(m, rows)['cells']
        target = cells[0]
        peers = [c for c in cells if c['domain'] == target['domain'] and
                 c['case'] == target['case'] and c['threads'] == target['threads'] and c['peer'] != 'baseline']
        peers[0].update(median=1.2, baseline_peer_median=1.3)
        target.update(median=1.02, wall_median=1.02)
        self.assertTrue(evaluate(cells, all_qualified=True)['incremental_passed'])
        target['wall_median'] = 1.020001
        self.assertFalse(evaluate(cells, all_qualified=True)['incremental_passed'])
        target['wall_median'] = 1.
        peers[0]['baseline_peer_median'] = 1.
        self.assertFalse(evaluate(cells, all_qualified=True)['incremental_passed'])
        self.assertEqual(evaluate(cells, all_qualified=True), hybrid_acceptance(cells, True))

    def test_paired_ratios_and_malformed_cells(self):
        m, rows = self.fixture()
        for row in rows:
            if row['side'] == 'ck':
                row['linear_seconds'] = [.5 if row['trial'] % 2 else 2.]
        result = self.check(m, rows, False)
        self.assertTrue(all(c['baseline_peer_median'] == 1.25 for c in result['cells'] if c['peer'] == 'ck'))
        for change in ('duplicate', 'missing', 'incomplete'):
            cells = copy.deepcopy(result['cells'])
            if change == 'duplicate': cells.append(cells[0])
            elif change == 'missing': cells.pop()
            else: cells[0]['complete'] = False
            self.assertFalse(evaluate(cells, all_qualified=True)['complete'])
            self.assertEqual(evaluate(cells, all_qualified=True), hybrid_acceptance(cells, True))

    def test_confirmation_requires_same_inputs_and_opposite_orders(self):
        with tempfile.TemporaryDirectory() as temporary:
            roots = [Path(temporary) / name for name in ('initial', 'confirmation')]
            metadata = []
            for root in roots:
                (root / 'binaries').mkdir(parents=True)
                for side in ('baseline', 'candidate'):
                    (root / 'binaries' / side).write_bytes(side.encode())
                metadata.append(dict(order=root.name, sides=['baseline', 'candidate', 'ck', 'st'],
                    cases=list(campaign.CASES), domains=list(campaign.DOMAINS), threads=[1, 8],
                    passes=4, timeout=1800, primary_metric='linear', wall_metric='wall',
                    accuracy='componentwise', complete_grid=True, diagnostic_only=False,
                    input_sha256={'input': 'same'}))
            for field in ('cases', 'domains'):
                metadata[1][field].reverse()
            self.assertTrue(matching_campaigns(metadata, roots))
            original = copy.deepcopy(metadata)
            for field, value in (('passes', 8), ('accuracy', 'application-managed'),
                                 ('order', 'initial'), ('input_sha256', {'input': 'different'})):
                metadata = copy.deepcopy(original)
                metadata[1][field] = value
                self.assertFalse(matching_campaigns(metadata, roots))
            (roots[1] / 'binaries/candidate').write_bytes(b'changed')
            self.assertFalse(matching_campaigns(original, roots))


if __name__ == '__main__':
    unittest.main()
