"""Independently reconstruct the complete AUTO broad regression ledger.

Imports no campaign code and never accepts the launcher's valid flags as proof.
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--prefix',action='store_true',help='Audit an immutable complete-matrix prefix without a merge decision')
    args = parser.parse_args()
    metadata = json.loads((args.campaign / 'metadata.json').read_text())
    assert all(sha(p) == h for p, h in metadata['input_sha256'].items())
    assert metadata['expected'] == 1760 and metadata['passes'] == 4
    assert metadata['timeout_seconds'] == 900 and metadata['residual_limit'] == 1e-8
    assert metadata['lifecycle_systems'] == 100 and metadata['matched_benchmark_object']
    build_path = args.campaign / 'build-commands.json'
    assert sha(build_path) == metadata['build_commands_sha256']
    builds = json.loads(build_path.read_text())
    assert [r['stage'] for r in builds] == ['compile', 'baseline-link', 'candidate-link']
    assert all(r['exit_code'] == 0 and r['closure']['all_session_processes_absent'] for r in builds)
    obj = metadata['benchmark_object']
    assert builds[0]['command'].count(obj) == 1 and '-c' in builds[0]['command']
    driver = Path(metadata['candidate_core']) / 'source/bench/kls_bench.c'
    assert sha(driver) == metadata['driver_sha256'] and builds[0]['command'].count(str(driver)) == 1
    normalized = []
    for side, record in zip(('baseline', 'candidate'), builds[1:]):
        core = Path(metadata[side + '_core'])
        replacements = {metadata['binaries'][side]: '{binary}', str(core / 'native/libkls.a'): '{core}',
            metadata['ordering_archives'][side]: '{ordering}'}
        assert record['command'].count(obj) == 1
        assert all(record['command'].count(token) == 1 for token in replacements)
        normalized.append(['{depfile}' if token.startswith('-Wl,--dependency-file=') else
            replacements.get(token, token) for token in record['command']])
    assert normalized[0] == normalized[1], 'benchmark object and common dependencies must match'
    names = []
    for line in Path(metadata['manifest']).read_text().splitlines():
        text = line.split('#')[0].strip()
        if text:
            names.append(Path(text.split()[0]).stem)
    assert len(names) == len(set(names)) == 110
    configs = json.loads(Path(metadata['configs_path']).read_text())
    assert configs == metadata['configs']
    assert {(c['cpus'], c['threads'], c['lifecycle_systems'], c['values'], c['amplitude']) for c in configs} == {
        ('8-15', 8, 100, 'entrywise', .001), ('0-7', 8, 100, 'entrywise', .001)}
    expected = []
    for config in configs:
        for index, name in enumerate(names):
            for trial in range(4):
                sides = ['baseline', 'candidate'] if (index + trial) % 2 == 0 else ['candidate', 'baseline']
                for position, side in enumerate(sides):
                    expected.append((config['id'], name, trial, side, position))
    rows = [json.loads(line) for line in (args.campaign / 'observations.jsonl').read_text().splitlines()]
    if args.prefix:
        assert 0<len(rows)<1760 and len(rows)%8==0
        expected=expected[:len(rows)]
    assert [(r['config_id'], r['matrix_id'], r['trial'], r['side'], r['position']) for r in rows] == expected
    assert len(rows)==len(expected)
    inherited=0
    if 'continuation' in metadata:
        c=metadata['continuation'];origin=Path(c['source'])
        original=json.loads((origin/'metadata.json').read_text())
        original_bytes=(origin/'observations.jsonl').read_bytes()
        original_rows=[json.loads(line) for line in original_bytes.splitlines()]
        inherited=c['observations']
        assert inherited==len(original_rows) and inherited%8==0 and inherited>0
        assert hashlib.sha256(original_bytes).hexdigest()==c['prefix_sha256']
        assert (args.campaign/'observations.jsonl').read_bytes().startswith(original_bytes)
        assert rows[:inherited]==original_rows
        assert c['identical_binaries'] and c['complete_matrix_boundary'] and c['all_recorded_failures_retained']
        assert c['abandoned_unrecorded_attempt_included_in_performance'] is False
        assert all(original[k]==metadata[k] for k in ('binaries','benchmark_object','ordering_archives',
            'baseline_core','candidate_core','dependency_build','driver_sha256','build_commands_sha256',
            'configs','matrix_paths','lifecycle_systems','passes','timeout_seconds','residual_limit'))
        for path in (origin/'metadata.json',origin/'observations.jsonl',origin/'build-commands.json',Path(c['failure_receipt'])):
            assert metadata['input_sha256'][str(path)]==sha(path)
        assert all(metadata['input_sha256'][p]==h for p,h in original['input_sha256'].items())
        failure=json.loads(Path(c['failure_receipt']).read_text())
        assert failure['campaign']==str(origin) and failure['completed_observations']==inherited
        assert failure['terminal_exit_code']!=0 and failure['ledger_sha256']==sha(origin/'observations.jsonl')
        assert failure['later_strict_cleanup']['all_session_processes_absent']
        assert failure['failure']=='RuntimeError: unrecognized process in completed child session'
        st=failure['stage'];nxt=expected[inherited]
        assert (st['config'],st['matrix'],st['trial'],st['side'],0)==nxt
        assert st['side']=='baseline' and st['trial']==0
        assert st['numeric_outcome']=='not persisted before cleanup raised'

    by_config = {c['id']: c for c in configs}
    recomputed = {}
    invalid = []
    for ordinal,row in enumerate(rows):
        directory = Path(row['directory'])
        assert sha(directory / 'stdout.log') == row['stdout_sha256']
        assert sha(directory / 'stderr.log') == row['stderr_sha256']
        command = row['command']
        assert command[:3] == ['taskset', '-c', by_config[row['config_id']]['cpus']]
        binary = metadata['binaries'][row['side']]
        assert command.count(binary) == 1
        tail = command[command.index(binary) + 1:]
        matrix = Path('/work') / Path(metadata['matrix_paths'][row['matrix_id']]).relative_to(metadata['workspace'])
        assert tail[0] == str(matrix)
        assert tail[1:] == ['--lifecycle-systems', '100', '--threads', '8', '--backend', 'auto',
            '--ordering', 'auto', '--orientation', 'auto', '--scale', 'auto', '--input-index', '64',
            '--refactor-values', 'entrywise', '--refactor-value-amplitude', '0.001', '--no-transpose-solve', '--json']
        assert row['closure']['all_session_processes_absent']
        if ordinal>=inherited and not args.prefix:
            assert set(row['lifecycle_receipt_sha256'])=={'started.json','numeric-exit.json','closure.json'}
            assert all(sha(directory/name)==digest for name,digest in row['lifecycle_receipt_sha256'].items())
            started=json.loads((directory/'started.json').read_text())
            exited=json.loads((directory/'numeric-exit.json').read_text())
            assert started['command']==command and started['timeout_seconds']==900
            assert started['pid']==exited['pid']==row['closure']['session_id']
            assert exited['returncode']==row['returncode'] and exited['timed_out']==row['timed_out']
            assert json.loads((directory/'closure.json').read_text())==row['closure']
        result = None
        for line in (directory / 'stdout.log').read_text().splitlines():
            try:
                value = json.loads(line)
                if isinstance(value, dict):
                    result = value
            except json.JSONDecodeError:
                pass
        valid = row['returncode'] == 0 and not row['timed_out'] and result is not None
        if valid:
            valid = result.get('status') in (None, 0) and result.get('verify_each_refactor') is True
            valid = valid and result.get('lifecycle_mode') == 'direct' and result.get('lifecycle_systems') == 100
            for field in ('relative_residual_l2', 'refactor_max_relative_residual', 'measured_lifecycle_seconds'):
                try:
                    number = float(result[field])
                    valid = valid and math.isfinite(number) and (number > 0 if field == 'measured_lifecycle_seconds' else 0 <= number <= 1e-8)
                except (KeyError, ValueError, TypeError):
                    valid = False
        assert valid == row['valid']
        if row['returncode'] == 0:
            assert result == row['result']
        key = (row['config_id'], row['matrix_id'], row['trial'], row['side'])
        assert key not in recomputed
        recomputed[key] = dict(valid=valid, result=result)
        if not valid:
            invalid.append(key)
    groups, regressions = [], []
    for config in configs:
        for name in names:
            if args.prefix and (config['id'],name,3,'candidate') not in recomputed:continue
            ratios = []
            for trial in range(4):
                a, b = (recomputed[config['id'], name, trial, side] for side in ('candidate', 'baseline'))
                if b['valid'] and not a['valid']:
                    regressions.append((config['id'], name, trial))
                if a['valid'] and b['valid']:
                    ratios.append(a['result']['measured_lifecycle_seconds'] / b['result']['measured_lifecycle_seconds'])
            groups.append(dict(config=config['id'], matrix=name, valid_pairs=len(ratios), ratios=ratios,
                median=statistics.median(ratios) if len(ratios) == 4 else None))
    if args.prefix:
        args.output.mkdir(exist_ok=False)
        (args.output/'audit.json').write_text(json.dumps(dict(observations=len(rows),
            raw_prefix_revalidated=True,complete_matrix_boundary=True,merge_decision=False,
            coverage_regressions=regressions,invalid=invalid,groups=groups,
            input_sha256={str(p):sha(p) for p in (Path(__file__),args.campaign/'metadata.json',
                args.campaign/'observations.jsonl')}),indent=2)+'\n')
        print('Validated complete matrix prefix:',len(rows),'observations; no merge decision.')
        return
    summary = json.loads((args.campaign / 'summary.json').read_text())
    assert groups == summary['groups']
    assert regressions == [(r['config'], r['matrix'], r['trial']) for r in summary['coverage_regressions']]
    assert len(invalid) == len(summary['invalid'])
    assert summary['all_matrices_qualified'] == (not invalid)
    geomean = math.exp(statistics.mean(math.log(g['median']) for g in groups)) if not invalid else None
    assert summary['candidate_over_baseline_geomean'] == geomean
    assert summary['observations'] == summary['expected'] == 1760 and not summary['changed_inputs']
    args.output.mkdir(exist_ok=False)
    (args.output / 'audit.json').write_text(json.dumps(dict(observations=1760, original_manifest_complete=True,
        independent_exact_agreement=True, coverage_regressions=regressions, invalid=invalid,
        all_matrices_qualified=not invalid, groups=groups,
        input_sha256={str(p): sha(p) for p in [Path(__file__), args.campaign / 'metadata.json',
            args.campaign / 'observations.jsonl', args.campaign / 'summary.json']}), indent=2) + '\n')


if __name__ == '__main__':
    main()
