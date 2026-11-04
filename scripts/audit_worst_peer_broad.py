"""Reconstruct fresh standalone peer exceptions from pinned raw logs."""
import argparse
import json
import math
from pathlib import Path

from worst_peer_broad_evidence import bind_original, read, sha
from worst_peer_broad_gate import evaluate


def raw_result(directory, exit_code, timed_out):
    if timed_out or exit_code != 0:
        return None, False
    result = None
    for line in reversed((directory / 'stdout.log').read_text().splitlines()):
        try:
            value = json.loads(line)
        except json.JSONDecodeError:
            continue
        if isinstance(value, dict):
            result = value
            break
    if result is None:
        return None, False
    valid = (result.get('status') in (None, 0)
             and result.get('verify_each_refactor') is True
             and result.get('lifecycle_mode') == 'direct'
             and result.get('lifecycle_systems') == 100)
    for name in ('relative_residual_l2', 'refactor_max_relative_residual',
                 'measured_lifecycle_seconds'):
        try:
            number = float(result[name])
            valid &= math.isfinite(number) and (
                number > 0 if name == 'measured_lifecycle_seconds'
                else 0 <= number <= 1e-8)
        except (KeyError, TypeError, ValueError, OverflowError):
            valid = False
    return result, bool(valid)


def audit(campaign):
    campaign = Path(campaign).resolve()
    meta, done = read(campaign / 'metadata.json'), read(campaign / 'completed.json')
    assert done['metadata_sha256'] == sha(campaign / 'metadata.json')
    assert done['ledger_sha256'] == sha(campaign / 'observations.jsonl')
    assert meta['acceptance_decision'] is done['acceptance_decision'] is False
    assert meta['timeout_seconds'] == 900 and meta['lifecycle_systems'] == 100
    assert meta['residual_limit'] == 1e-8
    assert all(sha(p) == h for p, h in meta['input_sha256'].items())
    bound = bind_original(meta['original_campaign'], campaign / 'original-audit/audit.json')
    targets = [(g['config'], g['matrix']) for g in bound['targets']]
    assert meta['targets'] == [list(t) for t in targets]
    original = bound['metadata']
    assert meta['configs'] == original['configs']
    assert meta['matrix_paths'] == original['matrix_paths']
    assert meta['binaries']['candidate'] == bound['candidate']
    peers = read(Path(meta['peer_build']) / 'build.json')
    for name in ('build.json', 'commands.json'):
        path = str(Path(meta['peer_build']) / name)
        assert meta['input_sha256'][path] == sha(path)
    assert peers['passed'] is True and peers['candidate_core'] == original['candidate_core']
    assert set(meta['binaries']) == {'candidate', 'ck', 'st'}
    assert all(meta['input_sha256'].get(path) == sha(path)
               for path in meta['binaries'].values())
    assert all(meta['binaries'][s] == peers['binaries'][s] for s in ('ck', 'st'))
    assert peers['commands_sha256'] == sha(Path(meta['peer_build']) / 'commands.json')
    for path, digest in {**bound['input_sha256'], **peers['input_sha256']}.items():
        assert meta['input_sha256'][path] == digest == sha(path)
    assert meta['environment_overrides'] == dict(OPENBLAS_NUM_THREADS='1',
        OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', KLS_BENCH_VERIFY_EACH_REFACTOR='1',
        BENCH_VERIFY_EACH_REFACTOR='1')
    workspace = Path(original['workspace'])
    image = workspace / 'KLS/build/kls-xyce-private.sif'
    prefix = meta['prefix']
    assert len(prefix) == 13 and prefix[:3] == ['apptainer', 'exec', '--bind']
    assert prefix[3:9] == [f'{workspace}:/work', '--bind', f'{workspace}:{workspace}',
                          '--bind', '/tmp:/tmp', '--bind']
    assert prefix[9].endswith(':/peer-license:ro')
    assert prefix[10:] == ['--pwd', '/peer-license', str(image)], prefix
    rows = [json.loads(line) for line in (campaign / 'observations.jsonl').read_text().splitlines()]
    # Reconstruct schedule and commands without importing the launcher/parser.
    expected = []
    for index, (config, matrix) in enumerate(targets):
        for trial in range(4):
            order = ['candidate', 'ck', 'st']
            if (index + trial) % 2:
                order.reverse()
            expected.extend((config, matrix, trial, side, pos)
                            for pos, side in enumerate(order))
    assert len(rows) == meta['expected'] == done['observations'] == len(expected)
    configs = {c['id']: c for c in meta['configs']}
    indexed, invalid = {}, []
    for row, identity in zip(rows, expected):
        assert tuple(row[k] for k in ('config', 'matrix', 'trial', 'side', 'position')) == identity
        config, matrix, trial, side, _ = identity
        directory = campaign / config / matrix / str(trial) / side
        assert row['directory'] == str(directory)
        assert set(row['files_sha256']) == {'stdout.log', 'stderr.log', 'started.json',
                                           'numeric-exit.json', 'closure.json'}
        assert all(sha(directory / n) == h for n, h in row['files_sha256'].items())
        binary, path = meta['binaries'][side], meta['matrix_paths'][matrix]
        tail = ([binary, path, '--lifecycle-systems', '100', '--threads', '8',
                 '--backend', 'auto', '--ordering', 'auto', '--orientation', 'auto',
                 '--scale', 'auto', '--input-index', '64', '--refactor-values',
                 'entrywise', '--refactor-value-amplitude', '0.001',
                 '--no-transpose-solve', '--json'] if side == 'candidate' else
                [binary, path, '8', '1', '99', '0', 'entrywise', '0.001',
                 '--lifecycle-systems', '100'])
        command = ['taskset', '-c', configs[config]['cpus'], *prefix, *tail]
        assert row['command'] == read(directory / 'started.json')['command'] == command
        exited = read(directory / 'numeric-exit.json')
        assert row['returncode'] == exited['returncode']
        assert row['timed_out'] == exited['timed_out']
        assert read(directory / 'closure.json')['all_session_processes_absent'] is True
        result, valid = raw_result(directory, row['returncode'], row['timed_out'])
        assert result == row['result'] and valid == row['valid']
        if not valid:
            invalid.append(identity)
        indexed[config, matrix, trial, side] = (result, valid)
    groups = []
    for config, matrix in targets:
        ratios = {'ck': [], 'st': []}
        for trial in range(4):
            candidate, valid = indexed[config, matrix, trial, 'candidate']
            for side in ('ck', 'st'):
                peer, peer_valid = indexed[config, matrix, trial, side]
                if valid and peer_valid:
                    ratios[side].append(float(candidate['measured_lifecycle_seconds']) /
                                        float(peer['measured_lifecycle_seconds']))
        groups.append(dict(config=config, matrix=matrix, ratios=ratios))
    decision = evaluate(bound['groups'], groups, expected_cells=bound['expected_cells'],
                        coverage_qualified=True, peer_evidence_qualified=not invalid)
    pins = dict(meta['input_sha256'])
    for name in ('worst_peer_broad_evidence.py', 'worst_peer_broad_gate.py'):
        path = Path(__file__).with_name(name).resolve()
        pins[str(path)] = sha(path)
    pins.update({str(p): sha(p) for p in (
        Path(__file__).resolve(), campaign / 'metadata.json',
        campaign / 'observations.jsonl', campaign / 'completed.json')})
    return dict(decision=decision, invalid=invalid, peer_groups=groups,
                independent_raw_log_audit=True, original_observations=1760,
                fresh_observations=len(rows), input_sha256=pins)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.campaign)
    args.output.mkdir(exist_ok=False)
    (args.output / 'audit.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result['decision']))


if __name__ == '__main__':
    main()
