import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from worst_peer_broad_evidence import bind_original, sha


class EvidenceBindingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.campaign = self.root / 'campaign'
        self.campaign.mkdir()
        self.audit = self.root / 'audit.json'
        self.binary = self.root / 'candidate'
        self.binary.write_bytes(b'candidate')
        matrices = {}
        for i in range(110):
            p = self.root / f'm{i}.mtx'
            p.write_text(str(i))
            matrices[f'm{i}'] = str(p)
        configs = [dict(id=f'headline_llc{size}_t8', cpus=cpus, threads=8,
                        lifecycle_systems=100, values='entrywise', amplitude=.001,
                        env={}) for size, cpus in ((32, '8-15'), (96, '0-7'))]
        groups = [dict(config=c['id'], matrix=m, valid_pairs=4,
                       ratios=[1.] * 4, median=1.) for c in configs for m in matrices]
        self.metadata = dict(expected=1760, matrix_count=110, passes=4,
            lifecycle_systems=100, timeout_seconds=900, residual_limit=1e-8,
            configs=configs, matrix_paths=matrices,
            binaries=dict(candidate=str(self.binary)),
            input_sha256={p: sha(p) for p in [str(self.binary), *matrices.values()]})
        self.summary = dict(observations=1760, expected=1760, changed_inputs=[],
                            coverage_regressions=[], invalid=[], groups=groups)
        self.save()

    def save(self):
        for name, data in (('metadata.json', self.metadata),
                           ('summary.json', self.summary)):
            (self.campaign / name).write_text(json.dumps(data))
        (self.campaign / 'observations.jsonl').write_text('retained ledger')
        auditor = Path(__file__).resolve().parents[1] / 'scripts/audit_auto_broad.py'
        pins = {str(self.campaign / n): sha(self.campaign / n) for n in
                ('metadata.json', 'summary.json', 'observations.jsonl')}
        pins[str(auditor)] = sha(auditor)
        self.audit.write_text(json.dumps(dict(observations=1760,
            original_manifest_complete=True, independent_exact_agreement=True,
            coverage_regressions=[], invalid=self.summary['invalid'],
            groups=self.summary['groups'], input_sha256=pins)))

    def test_original_regressions_remain_targets(self):
        self.summary['groups'][0].update(ratios=[1.3]*4, median=1.3)
        self.save()
        result = bind_original(self.campaign, self.audit)
        self.assertEqual(len(result['expected_cells']), 220)
        self.assertEqual(result['targets'][0]['median'], 1.3)

    def test_changed_candidate_or_matrix_is_rejected(self):
        for p in (self.binary, Path(self.metadata['matrix_paths']['m0'])):
            original = p.read_bytes()
            p.write_bytes(original + b'changed')
            with self.assertRaises(AssertionError):
                bind_original(self.campaign, self.audit)
            p.write_bytes(original)

    def test_missing_group_and_forged_median_are_rejected(self):
        group = self.summary['groups'].pop()
        self.save()
        with self.assertRaises(AssertionError):
            bind_original(self.campaign, self.audit)
        self.summary['groups'].append(group)
        group['median'] = .5
        self.save()
        with self.assertRaises(AssertionError):
            bind_original(self.campaign, self.audit)

    def test_changed_ledger_and_coverage_failure_are_rejected(self):
        (self.campaign / 'observations.jsonl').write_text('changed ledger')
        with self.assertRaises(AssertionError):
            bind_original(self.campaign, self.audit)
        self.summary['coverage_regressions'] = [dict(matrix='m0')]
        self.save()
        with self.assertRaises(AssertionError):
            bind_original(self.campaign, self.audit)


if __name__ == '__main__':
    unittest.main()
