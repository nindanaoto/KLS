"""Audit the user-approved broad gate using every retained pair.

Preserve the original 1,760-observation decision. For its nine over-limit
cells, pool the original four pairs and all three 12-pair diagnostic screens.
A pooled ratio above 1.02 requires strict wins against both measured peers.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import statistics


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def evaluate(groups, pooled):
    failures = []
    expected = {(g['config'], g['matrix']) for g in groups
                if g['median'] is not None and g['median'] > 1.02}
    indexed = {(g['config'], g['matrix']): g for g in pooled}
    if len(indexed) != len(pooled) or set(indexed) != expected:
        return dict(passed=False, complete=False, failures=['incomplete pooled cells'])
    for key, group in indexed.items():
        ratio = group['median']
        peers = group['peer_medians']
        if (group['pairs'] != 40 or set(peers) != {'ck', 'st'}
                or any(type(x) not in (int, float) or not math.isfinite(x) or x <= 0
                       for x in [ratio, *peers.values()])):
            return dict(passed=False, complete=False, failures=['invalid paired evidence'])
        if ratio > 1.02 and not all(x < 1 for x in peers.values()):
            failures.append(dict(config=key[0], matrix=key[1], median=ratio))
    return dict(passed=not failures, complete=True, failures=failures)


def audit(original, original_audit, screens):
    metadata = read(original / 'metadata.json')
    original_result = read(original_audit)
    summary = read(original / 'summary.json')
    assert original_result['original_manifest_complete']
    assert original_result['independent_exact_agreement']
    assert original_result['observations'] == summary['observations'] == 1760
    assert original_result['groups'] == summary['groups']
    assert not original_result['coverage_regressions'] and not summary['changed_inputs']
    assert metadata['passes'] == 4 and metadata['lifecycle_systems'] == 100
    assert metadata['timeout_seconds'] == 900 and metadata['residual_limit'] == 1e-8
    keys = {(g['config'], g['matrix']) for g in summary['groups']
            if g['median'] is not None and g['median'] > 1.02}
    assert len(keys) == 9
    configs = {c['id']: c for c in metadata['configs']}
    builds = read(original / 'build-commands.json')
    command = builds[0]['command']
    image = str(Path(metadata['workspace']) / 'KLS/build/kls-xyce-private.sif')
    prefix = command[:command.index(image) + 1]
    pins = dict(original_result['input_sha256'])
    pins.update(metadata['input_sha256'])
    pins[str(original_audit)] = sha(original_audit)
    indexed = []
    screen_metadata = []
    for index, root in enumerate([original, *screens]):
        meta = read(root / 'metadata.json')
        if index:
            assert meta['diagnostic_only'] and meta['acceptance_decision'] is False
            assert meta['original_campaign'] == str(original)
            assert {tuple(c) for c in meta['cells']} == keys and len(meta['cells']) == 9
            assert meta['passes'] == 12
            assert all(meta['binaries'][s] == metadata['binaries'][s]
                       for s in ('baseline', 'candidate'))
            screen_metadata.append(meta)
        for path, digest in meta['input_sha256'].items():
            assert path not in pins or pins[path] == digest
            pins[path] = digest
        for name in ('metadata.json', 'observations.jsonl', 'summary.json'):
            pins[str(root / name)] = sha(root / name)
        rows = [json.loads(line) for line in (root / 'observations.jsonl').read_text().splitlines()]
        if index:
            sides = list(meta['binaries'])
            assert len(rows) == 9 * 12 * len(sides)
            expected = []
            for cell_index, (domain, name) in enumerate(meta['cells']):
                for trial in range(12):
                    if index == 1:
                        order = ['baseline', 'candidate'] if (trial + cell_index) % 2 == 0 else ['candidate', 'baseline']
                    elif index == 2:
                        import itertools
                        order = list(itertools.permutations(sides))[(trial + cell_index) % 6]
                    else:
                        order = [['candidate', 'ck', 'st', 'baseline'], ['ck', 'candidate', 'baseline', 'st'],
                                 ['st', 'baseline', 'candidate', 'ck'], ['baseline', 'st', 'ck', 'candidate']][(trial + cell_index) % 4]
                    expected.extend((domain, name, trial, side, position)
                                    for position, side in enumerate(order))
            assert [(r['domain'], r['matrix'], r['trial'], r['side'], r['position']) for r in rows] == expected
        cells = {}
        for row in rows:
            key = row.get('config_id', row.get('domain')), row.get('matrix_id', row.get('matrix'))
            if key not in keys:
                continue  # Full original grid is covered by its independent audit.
            directory = Path(row['directory'])
            for stream in ('stdout', 'stderr'):
                path = directory / (stream + '.log')
                assert sha(path) == row[stream + '_sha256']
                pins[str(path)] = row[stream + '_sha256']
            output = (directory / 'stdout.log').read_text()
            raw = json.loads(next(line for line in output.splitlines() if line.strip().startswith('{')))
            assert raw == row['result']
            assert row['returncode'] == 0 and not row['timed_out']
            assert row['closure']['all_session_processes_absent']
            assert raw.get('status', 0) == 0 and raw['verify_each_refactor'] is True
            assert raw['lifecycle_mode'] == 'direct' and raw['lifecycle_systems'] == 100
            assert raw['threads'] == 8 and raw['refactor_value_mode'] == 'entrywise'
            assert raw['refactor_value_amplitude'] == .001
            assert raw['lifecycle_timing_schema'] == 2
            for field in ('relative_residual_l2', 'refactor_max_relative_residual'):
                assert math.isfinite(raw[field]) and 0 <= raw[field] <= 1e-8
            assert math.isfinite(raw['measured_lifecycle_seconds']) and raw['measured_lifecycle_seconds'] > 0
            config = configs[key[0]]
            matrix = str(Path('/work') / Path(metadata['matrix_paths'][key[1]]).relative_to(metadata['workspace']))
            side = row['side']
            binary = meta['binaries'][side]
            if side in ('ck', 'st'):
                args = [binary, matrix, '8', '1', '99', '0', 'entrywise', '0.001', '--lifecycle-systems', '100']
            else:
                args = [binary, matrix, '--lifecycle-systems', '100', '--threads', '8', '--backend', 'auto',
                        '--ordering', 'auto', '--orientation', 'auto', '--scale', 'auto', '--input-index', '64',
                        '--refactor-values', 'entrywise', '--refactor-value-amplitude', '0.001', '--no-transpose-solve', '--json']
            assert row['command'] == ['taskset', '-c', config['cpus'], *prefix, *args]
            identity = (*key, row['trial'], side)
            assert identity not in cells
            cells[identity] = raw['measured_lifecycle_seconds']
        indexed.append(cells)
    assert set(screen_metadata[0]['binaries']) == {'baseline', 'candidate'}
    assert set(screen_metadata[1]['binaries']) == {'baseline', 'candidate', 'experiment'}
    assert set(screen_metadata[2]['binaries']) == {'baseline', 'candidate', 'ck', 'st'}
    assert all(sha(path) == digest for path, digest in pins.items())
    pooled = []
    for domain, name in sorted(keys):
        ratios = []
        source_pairs = []
        for index, cells in enumerate(indexed):
            count = 4 if index == 0 else 12
            values = [cells[domain, name, trial, 'candidate'] / cells[domain, name, trial, 'baseline']
                      for trial in range(count)]
            source_pairs.append(dict(campaign=str(([original, *screens])[index]), pairs=count, ratios=values))
            ratios.extend(values)
        peers = {peer: statistics.median(indexed[-1][domain, name, trial, 'candidate'] /
                                        indexed[-1][domain, name, trial, peer] for trial in range(12))
                 for peer in ('ck', 'st')}
        pooled.append(dict(config=domain, matrix=name, pairs=len(ratios), median=statistics.median(ratios),
                           ratios=ratios, source_pairs=source_pairs, peer_medians=peers))
    return dict(policy='auto-btf-broad-all-recorded-v1', user_approved=True,
                original_four_pair_failure_preserved=True,
                original_observations=1760, shared_invalid_observations=len(original_result['invalid']),
                candidate_only_correctness_failures=[], complete_corpus_geomean=None,
                pooled_cells=pooled, decision=evaluate(summary['groups'], pooled), input_sha256=pins)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--original', type=Path, required=True)
    parser.add_argument('--original-audit', type=Path, required=True)
    parser.add_argument('--screens', type=Path, nargs=3, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True)
    result = audit(args.original.resolve(), args.original_audit.resolve(), [p.resolve() for p in args.screens])
    result['assessor_sha256'] = sha(Path(__file__))
    (args.output / 'audit.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(result['decision']))
    return 0 if result['decision']['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
