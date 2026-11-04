import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from audit_worst_peer_broad import audit
from worst_peer_broad_evidence import bind_original, sha
import test_worst_peer_broad_evidence as fixtures


def save(path, data):
    path.write_text(json.dumps(data))


class RawAuditTests(unittest.TestCase):
    def setUp(self):
        f = fixtures.EvidenceBindingTests()
        f.setUp()
        self.addCleanup(f.doCleanups)
        self.f = f
        f.metadata.update(workspace=str(f.root), candidate_core=str(f.root / 'core'))
        f.summary['groups'][0].update(ratios=[1.2] * 4, median=1.2)
        f.save()
        root = self.root = f.root / 'fresh'
        root.mkdir()
        (root / 'original-audit').mkdir()
        (root / 'original-audit/audit.json').write_bytes(f.audit.read_bytes())
        bound = bind_original(f.campaign, root / 'original-audit/audit.json')
        build = f.root / 'peers'
        build.mkdir()
        save(build / 'commands.json', [])
        binaries = dict(candidate=str(f.binary))
        for side in ('ck', 'st'):
            binary = f.root / side
            binary.write_text(side)
            binaries[side] = str(binary)
        save(build / 'build.json', dict(passed=True, candidate_core=f.metadata['candidate_core'],
            binaries={s: binaries[s] for s in ('ck', 'st')}, input_sha256={},
            commands_sha256=sha(build / 'commands.json')))
        prefix = ['apptainer', 'exec', '--bind', f'{f.root}:/work', '--bind',
                  f'{f.root}:{f.root}', '--bind', '/tmp:/tmp', '--bind',
                  '/license:/peer-license:ro', '--pwd', '/peer-license',
                  str(f.root / 'KLS/build/kls-xyce-private.sif')]
        config = f.metadata['configs'][0]
        self.rows = []
        for trial in range(4):
            order = ['candidate', 'ck', 'st']
            if trial % 2:
                order.reverse()
            for position, side in enumerate(order):
                dest = root / config['id'] / 'm0' / str(trial) / side
                dest.mkdir(parents=True)
                binary, matrix = binaries[side], f.metadata['matrix_paths']['m0']
                tail = ([binary, matrix, '--lifecycle-systems', '100', '--threads', '8',
                    '--backend', 'auto', '--ordering', 'auto', '--orientation', 'auto',
                    '--scale', 'auto', '--input-index', '64', '--refactor-values',
                    'entrywise', '--refactor-value-amplitude', '0.001',
                    '--no-transpose-solve', '--json'] if side == 'candidate' else
                    [binary, matrix, '8', '1', '99', '0', 'entrywise', '0.001',
                     '--lifecycle-systems', '100'])
                command = ['taskset', '-c', config['cpus'], *prefix, *tail]
                result = dict(status=0, verify_each_refactor=True, lifecycle_mode='direct',
                    lifecycle_systems=100, relative_residual_l2=1e-12,
                    refactor_max_relative_residual=1e-12,
                    measured_lifecycle_seconds=.9 if side == 'candidate' else 1.)
                save(dest / 'stdout.log', result)
                (dest / 'stderr.log').write_text('')
                save(dest / 'started.json', dict(command=command))
                save(dest / 'numeric-exit.json', dict(returncode=0, timed_out=False))
                save(dest / 'closure.json', dict(all_session_processes_absent=True))
                self.rows.append(dict(config=config['id'], matrix='m0', trial=trial,
                    side=side, position=position, directory=str(dest), command=command,
                    returncode=0, timed_out=False, result=result, valid=True,
                    files_sha256={name: sha(dest / name) for name in (
                        'stdout.log', 'stderr.log', 'started.json', 'numeric-exit.json', 'closure.json')}))
        pins = dict(bound['input_sha256'])
        pins.update({path: sha(path) for path in binaries.values()})
        pins.update({str(build / n): sha(build / n) for n in ('build.json', 'commands.json')})
        save(root / 'metadata.json', dict(original_campaign=str(f.campaign),
            peer_build=str(build), binaries=binaries, targets=[[config['id'], 'm0']],
            expected=12, input_sha256=pins, configs=f.metadata['configs'],
            matrix_paths=f.metadata['matrix_paths'], prefix=prefix,
            environment_overrides=dict(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
                MKL_NUM_THREADS='1', KLS_BENCH_VERIFY_EACH_REFACTOR='1',
                BENCH_VERIFY_EACH_REFACTOR='1'), timeout_seconds=900,
            lifecycle_systems=100, residual_limit=1e-8, acceptance_decision=False))
        self.persist()

    def persist(self):
        (self.root / 'observations.jsonl').write_text(''.join(json.dumps(r) + '\n' for r in self.rows))
        save(self.root / 'completed.json', dict(observations=len(self.rows), acceptance_decision=False,
            metadata_sha256=sha(self.root / 'metadata.json'),
            ledger_sha256=sha(self.root / 'observations.jsonl')))

    def test_complete_peer_wins_admit_original_regression(self):
        result = audit(self.root)
        self.assertTrue(result['decision']['passed'])
        self.assertEqual(len(result['decision']['peer_exceptions']), 1)
        self.assertEqual(result['decision']['peer_exceptions'][0]['baseline_median'], 1.2)

    def test_forged_valid_flag_cannot_hide_bad_residual(self):
        row = self.rows[0]
        row['result']['refactor_max_relative_residual'] = .01
        dest = Path(row['directory'])
        save(dest / 'stdout.log', row['result'])
        row['files_sha256']['stdout.log'] = sha(dest / 'stdout.log')
        self.persist()
        with self.assertRaises(AssertionError):
            audit(self.root)

    def test_changed_peer_binary_is_rejected(self):
        meta = json.loads((self.root / 'metadata.json').read_text())
        Path(meta['binaries']['st']).write_text('changed peer')
        with self.assertRaises(AssertionError):
            audit(self.root)

    def test_missing_arm_and_reordered_ledger_are_rejected(self):
        self.rows[0], self.rows[1] = self.rows[1], self.rows[0]
        self.persist()
        with self.assertRaises(AssertionError):
            audit(self.root)
        self.rows.pop()
        self.persist()
        with self.assertRaises(AssertionError):
            audit(self.root)

    def test_correctly_reported_peer_failure_cannot_grant_exception(self):
        row = next(r for r in self.rows if r['side'] == 'st')
        row['result']['refactor_max_relative_residual'] = .01
        row['valid'] = False
        dest = Path(row['directory'])
        save(dest / 'stdout.log', row['result'])
        row['files_sha256']['stdout.log'] = sha(dest / 'stdout.log')
        self.persist()
        self.assertFalse(audit(self.root)['decision']['passed'])


if __name__ == '__main__':
    unittest.main()
