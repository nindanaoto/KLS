"""Campaign integrity: missing/failed cells cannot become performance wins."""
import copy
import importlib.util
from pathlib import Path
import tempfile
import subprocess
import sys
import unittest

spec = importlib.util.spec_from_file_location('campaign', Path(__file__).resolve().parents[1] / 'scripts/run_spice_campaign.py')
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)


class CampaignTests(unittest.TestCase):
    def test_child_timeout_and_normal_exit(self):
        child = subprocess.Popen([sys.executable, '-c', 'pass'], start_new_session=True)
        self.assertEqual(c.wait_child(child, 10), (0, False))
        child = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(10)'], start_new_session=True)
        code, timeout = c.wait_child(child, .05)
        self.assertTrue(timeout)
        self.assertNotEqual(code, 0)
        self.assertIsNotNone(child.returncode)

    def fixture(self):
        metadata = dict(cases=list(c.CASES), domains=list(c.DOMAINS), threads=[1, 8],
                        passes=4, sides=['baseline', 'candidate', 'ck', 'st'], complete_grid=True)
        rows = []
        for domain, _ in metadata['domains']:
            for case in metadata['cases']:
                rows.append(dict(domain=domain, case=case, threads=1, trial=-1, side='klu', qualified=True))
                for t in metadata['threads']:
                    for trial in range(4):
                        for side in metadata['sides']:
                            rows.append(dict(domain=domain, case=case, threads=t, trial=trial, side=side,
                                             qualified=True, linear_seconds=[.9 if side == 'candidate' else 1.],
                                             wall_seconds=.9 if side == 'candidate' else 1.))
        return metadata, rows

    def test_complete_and_failed_grid(self):
        m, rows = self.fixture()
        self.assertEqual(len(rows), 462)
        self.assertTrue(c.reduce_results(m, rows)['all56_peer_wins'])
        for changed in (rows[:-1], rows + [rows[-1]], rows[1:]):
            self.assertFalse(c.reduce_results(m, changed)['all56_peer_wins'])
        for index in (0, 1, 2, len(rows)-1):
            bad = copy.deepcopy(rows)
            bad[index]['qualified'] = False
            self.assertFalse(c.reduce_results(m, bad)['all56_peer_wins'])
        bad = copy.deepcopy(rows)
        bad[0]['case'] = 'unexpected.cir'
        self.assertFalse(c.reduce_results(m, bad)['all56_peer_wins'])

    def test_margin_and_wall_guard(self):
        m, rows = self.fixture()
        for row in rows:
            if row['side'] == 'candidate':
                row['linear_seconds'] = [.99]
        self.assertFalse(c.reduce_results(m, rows)['all56_peer_wins'])
        for row in rows:
            if row['side'] == 'candidate':
                row['linear_seconds'] = [.9]
                row['wall_seconds'] = 1.03
        self.assertFalse(c.reduce_results(m, rows)['all56_peer_wins'])

    def test_balanced_positions(self):
        for sides in (['baseline', 'candidate'], ['baseline', 'candidate', 'ck', 'st']):
            for confirmation in (True, False):
                for case in range(7):
                    orders = [c.launch_order(sides, case, t, confirmation) for t in range(4)]
                    for side in sides:
                        for pos in range(len(sides)):
                            self.assertEqual(sum(o[pos] == side for o in orders), 4 // len(sides))

    def test_waveform_and_staging(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)
            a, b = p/'a.prn', p/'b.prn'
            a.write_text('Index TIME V(1)\n0 0 1\n1 1 2\n')
            b.write_text(a.read_text())
            self.assertEqual(c.compare_waveforms(a,b)['status'], 'pass')
            for content in ('Index TIME V(1)\n0 0 nan\n1 1 2\n', 'Index TIME V(1)\n0 0 1\n1 2 2\n',
                            'Index TIME V(1)\n0 0 1\n1 1 3\n'):
                b.write_text(content)
                self.assertNotEqual(c.compare_waveforms(a,b)['status'], 'pass')
            source, license_file = p/'source.cir', p/'license'
            original = '* test\nR1 1 0 10\n.options linsol type=KLU OTHER=4\n.print tran V(1)\n.end\n'
            source.write_text(original)
            license_file.write_text('fixture')
            c.stage(source,p/'run','KLS',8,license_file)
            result = (p/'run/input.cir').read_text()
            self.assertEqual(source.read_text(), original)
            self.assertIn('OTHER=4',result)
            self.assertIn('KLS_THREADS=8',result)
            self.assertEqual(result.lower().count('type='),1)
            self.assertIn('R1 1 0 10',result)


if __name__ == '__main__':
    unittest.main()
