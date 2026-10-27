#!/usr/bin/env python3
"""Run/reduce complete, balanced Xyce SPICE comparisons without timing retries."""
from __future__ import annotations

import argparse
import fcntl
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import statistics
import subprocess
import time

CASES = (
    'PDE_2D_BJT/invertbjt.cir', 'MOS13_IC/invert50_mos1.cir',
    'CircuitSim90/MOS2_LARGE/mem_plus.cir', 'CircuitSim90/MOS2/mux8.cir',
    'CircuitSim90/MOS2_LARGE/dac.cir', 'NEURON/rallpack3.cir',
    'CircuitSim90/MOS2_LARGE/ram2k.cir',
)
DOMAINS = (('llc96', '0-7'), ('llc32', '8-15'))


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


def stage(source, destination, solver, threads, license_file):
    destination.mkdir(parents=True)
    text = source.read_text()
    if re.search(r'^\s*\.(?:include|inc|lib)\b', text, re.I | re.M):
        raise ValueError('campaign requires the original self-contained circuits')
    option = f'.options linsol type={solver}'
    if solver != 'KLU':
        option += f' {solver}_THREADS={threads}'
    if re.search(r'^\s*\.options\s+linsol\b', text, re.I | re.M):
        def replace_options(match):
            line = re.sub(r'\btype\s*=\s*\w+', '', match[0], flags=re.I)
            line = re.sub(r'\b(?:KLS|CKTSO|SUBTREELU)_THREADS\s*=\s*\d+', '', line, flags=re.I)
            return line + option[len('.options linsol'):]
        text = re.sub(r'^\s*\.options\s+linsol\b[^\n]*', replace_options, text, flags=re.I | re.M)
    else:
        text, count = re.subn(r'^(\s*\.end\s*)$', lambda m: option + '\n' + m[0], text, flags=re.I | re.M)
        if count != 1:
            raise ValueError('expected exactly one .end')
    def print_line(match):
        line = re.sub(r'\b(?:precision|width)\s*=\s*\d+', '', match[0], flags=re.I)
        return re.sub(r'^(\s*\.print\s+\w+)', r'\1 PRECISION=17 WIDTH=26', line, flags=re.I)
    text = re.sub(r'^\s*\.print\s+[^\n]*', print_line, text, flags=re.I | re.M)
    (destination / 'input.cir').write_text(text)
    shutil.copyfile(license_file, destination / 'cktso.lic')


def waveform(path):
    if not path.exists():
        return []
    result = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if fields and fields[0].isdigit():
            result.append([float(x) for x in fields[1:]])
    return result


def compare_waveforms(reference, output):
    a, b = waveform(reference), waveform(output)
    if not a or len(a) != len(b):
        return {'status': 'row_count_mismatch'}
    maximum = 0.0
    for ra, rb in zip(a, b):
        if len(ra) != len(rb) or len(ra) < 2:
            return {'status': 'column_count_mismatch'}
        if not all(math.isfinite(v) for v in (*ra, *rb)):
            return {'status': 'nonfinite'}
        # Adaptive stepping accumulates floating-point rounding even when the
        # accepted/rejected steps agree. Allow 128 binary64 rounding units;
        # never interpolate a different grid to make a comparison pass.
        if abs(ra[0] - rb[0]) > 128 * math.ulp(1.0) * max(abs(ra[0]), abs(rb[0]), 1e-12):
            return {'status': 'grid_mismatch'}
        for x, y in zip(ra[1:], rb[1:]):
            maximum = max(maximum, abs(x-y) / (1e-10 + 1e-8 * abs(x)))
    return {'status': 'pass' if maximum <= 1 else 'waveform_difference',
            'maximum_scaled_error': maximum, 'rows': len(a)}


def parse_log(text):
    fields = {}
    for key, pattern in {
        'linear_seconds': r'Total Linear Solution Time:\s*([\d.eE+-]+)',
        'successful_steps': r'Number Successful Steps Taken:\s*(\d+)',
        'failed_steps': r'Number Failed Steps Attempted:\s*(\d+)',
    }.items():
        values = re.findall(pattern, text)
        fields[key] = [float(v) if key == 'linear_seconds' else int(v) for v in values]
    fields['linear_failure_messages'] = len(re.findall('linear solve failed', text, re.I))
    fields['end_of_simulation'] = 'End of Xyce(TM) Simulation' in text
    return fields


def launch_order(sides, case_index, trial, confirmation):
    shift = (case_index + trial + int(not confirmation)) % len(sides)
    return sides[shift:] + sides[:shift]


def reduce_results(metadata, rows):
    indexed = {}
    duplicates = []
    for r in rows:
        key = (r['domain'], r['case'], r['threads'], r['trial'], r['side'])
        if key in indexed:
            duplicates.append(key)
        indexed[key] = r
    cells = []
    for domain, _ in metadata['domains']:
        for case in metadata['cases']:
            for threads in metadata['threads']:
                for peer in metadata['sides']:
                    if peer == 'candidate':
                        continue
                    ratios, wall = [], []
                    for trial in range(metadata['passes']):
                        prefix = (domain, case, threads, trial)
                        a, b = indexed.get((*prefix, 'candidate')), indexed.get((*prefix, peer))
                        if a and b and a['qualified'] and b['qualified']:
                            ratios.append(a['linear_seconds'][-1] / b['linear_seconds'][-1])
                            wall.append(a['wall_seconds'] / b['wall_seconds'])
                    complete = len(ratios) == metadata['passes']
                    cells.append(dict(domain=domain, case=case, threads=threads, peer=peer,
                                      complete=complete, ratios=ratios,
                                      median=statistics.median(ratios) if complete else None,
                                      wall_median=statistics.median(wall) if complete else None))
    expected = len(metadata['domains']) * len(metadata['cases']) * (
        1 + len(metadata['threads']) * metadata['passes'] * len(metadata['sides']))
    expected_keys = set()
    for domain, _ in metadata['domains']:
        for case in metadata['cases']:
            expected_keys.add((domain, case, 1, -1, 'klu'))
            for threads in metadata['threads']:
                for trial in range(metadata['passes']):
                    for side in metadata['sides']:
                        expected_keys.add((domain, case, threads, trial, side))
    all_qualified = set(indexed) == expected_keys and not duplicates and all(r['qualified'] for r in rows)
    baseline = [c for c in cells if c['peer'] == 'baseline']
    peers = [c for c in cells if c['peer'] in ('ck', 'st')]
    official = (metadata['complete_grid'] and all_qualified and set(metadata['cases']) == set(CASES)
                and set(map(tuple, metadata['domains'])) == set(DOMAINS) and metadata['threads'] == [1, 8]
                and metadata['passes'] >= 4 and metadata['passes'] % 4 == 0)
    nonregression = official and len(baseline) == 28 and all(c['complete'] and c['median'] <= 1.02 and c['wall_median'] <= 1.02 for c in baseline)
    incremental = nonregression and any(c['median'] <= .98 for c in baseline)
    return dict(expected=expected, observations=len(rows), all_qualified=all_qualified,
                duplicates=duplicates, cells=cells, incremental_gate=incremental,
                all56_peer_wins=nonregression and len(peers) == 56 and all(c['complete'] and c['median'] <= .98 for c in peers),
                independent_confirmation_required=True,
                unqualified=[r for r in rows if not r['qualified']])


def run(args):
    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location('cleanup', root / 'bench/process_session_cleanup.py')
    cleanup = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cleanup)
    if not args.smoke and not args.diagnostic:
        busy = []
        for process in Path('/proc').iterdir():
            if not process.name.isdigit() or int(process.name) == os.getpid():
                continue
            try:
                name = (process / 'comm').read_text().strip()
                if name in ('cc1', 'cc1plus', 'f951', 'ctest', 'cmake', 'Xyce', 'kls_bench'):
                    busy.append((int(process.name), name))
            except (FileNotFoundError, ProcessLookupError, PermissionError):
                continue
        if busy:
            raise RuntimeError(f'timing requires idle build/test/solver processes: {busy}')
    output = args.output.resolve()
    output.mkdir(parents=True)  # Never overwrite/retry an earlier campaign.
    lock = open(args.workspace / '.spice-timing.lock', 'w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    binaries = {'baseline': args.baseline.resolve(), 'candidate': args.candidate.resolve()}
    sides = ['baseline', 'candidate'] + (['ck', 'st'] if args.peers else [])
    cases = list(CASES) if not args.smoke else [CASES[3]]
    domains = list(DOMAINS)
    if args.order == 'confirmation':
        domains.reverse()
    else:
        cases.reverse()
    threads = [1, 8] if not args.smoke else [1]
    frozen = output / 'binaries'
    frozen.mkdir()
    hashes = {}
    for name, original in binaries.items():
        hashes[str(original)] = sha(original)
        shutil.copy2(original, frozen / name)
        binaries[name] = frozen / name
        hashes[str(binaries[name])] = sha(binaries[name])
    inputs = [Path(__file__), root / 'bench/process_session_cleanup.py', args.image,
              args.license, args.regression / 'TestScripts/xyce_verify.pl',
              *[args.regression / 'Netlists' / c for c in cases], *args.peer_library]
    if args.provenance:
        inputs.append(args.provenance)
    hashes.update({str(p.resolve()): sha(p) for p in inputs})
    metadata = dict(cases=cases, domains=domains, threads=threads, passes=args.passes,
                    sides=sides, order=args.order, complete_grid=not args.smoke and not args.diagnostic and args.passes >= 4 and args.passes % 4 == 0,
                    diagnostic_only=args.diagnostic or args.smoke,
                    timeout=1800, input_sha256=hashes, cpu_topology=subprocess.check_output(['lscpu'], text=True),
                    primary_metric='final Total Linear Solution Time', wall_metric='process elapsed before teardown',
                    accuracy='default componentwise; waveform rtol=1e-8, atol=1e-10; exact steps and grid',
                    provenance=json.loads(args.provenance.read_text()) if args.provenance else None)
    write_json(output / 'metadata.json', metadata)
    env = {k: v for k, v in os.environ.items() if not k.startswith(
        ('LD_', 'KLS_', 'OMP_', 'GOMP_', 'OPENBLAS_', 'APPTAINERENV_', 'ASAN_', 'UBSAN_', 'TSAN_'))}
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', XYCE_NO_TRACKING='1')
    rows = []

    def execute(domain, cpus, case, count, trial, side, reference=None):
        dest = output / domain / Path(case).stem / f't{count}-p{trial}-{side}'
        solver = {'baseline': 'KLS', 'candidate': 'KLS', 'ck': 'CKTSO', 'st': 'SUBTREELU', 'klu': 'KLU'}[side]
        stage(args.regression / 'Netlists' / case, dest, solver, count, args.license)
        executable = binaries['candidate' if side == 'candidate' else 'baseline']
        cmd = ['taskset', '-c', cpus, 'apptainer', 'exec', '--bind', f'{args.workspace}:/work',
               '--bind', f'{args.workspace}:{args.workspace}', '--bind', '/tmp:/tmp',
               '--pwd', str(dest), str(args.image), str(executable), 'input.cir']
        start = time.monotonic()
        with (dest / 'stdout.log').open('x') as stream:
            child = subprocess.Popen(cmd, stdout=stream, stderr=subprocess.STDOUT, env=env, start_new_session=True)
            timeout = False
            try:
                code = child.wait(timeout=1800)
            except subprocess.TimeoutExpired:
                timeout = True
                os.killpg(child.pid, signal.SIGKILL)
                code = child.wait()
            except BaseException:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
                cleanup.drain_exited_apptainer_session(child.pid, args.image)
                raise
        wall = time.monotonic() - start
        closure = cleanup.drain_exited_apptainer_session(child.pid, args.image)
        row = dict(domain=domain, cpus=cpus, case=case, threads=count, trial=trial, side=side,
                   directory=str(dest), command=cmd, pid=child.pid, exit_code=code, timed_out=timeout,
                   wall_seconds=wall, closure=closure, netlist_sha256=sha(dest / 'input.cir'))
        row.update(parse_log((dest / 'stdout.log').read_text()))
        values = waveform(dest / 'input.cir.prn')
        row['qualified'] = bool(code == 0 and not timeout and closure['all_session_processes_absent'] and row['end_of_simulation'] and
                                not row['linear_failure_messages'] and row['successful_steps'] and row['failed_steps'] and
                                row['linear_seconds'] and all(math.isfinite(v) and v > 0 for v in row['linear_seconds']) and
                                values and all(math.isfinite(v) for r in values for v in r))
        if reference is not None:
            row['waveform'] = compare_waveforms(Path(reference['directory']) / 'input.cir.prn', dest / 'input.cir.prn')
            ref_output = Path(reference['directory']) / 'input.cir.prn'
            test_output = dest / 'input.cir.prn'
            row['upstream_verifier'] = None
            if ref_output.exists() and test_output.exists():
                names = ref_output.read_text().splitlines()[0].split()[2:]
                netlist = re.sub(r'^\s*\*comp[^\n]*\n', '', (dest / 'input.cir').read_text(), flags=re.I | re.M)
                comments = ''.join(f'*comp {name} reltol=1e-8 abstol=1e-10 zerotol=0 absdifftol=0 numfail=0\n' for name in names)
                netlist = re.sub(r'^(\s*\.end\s*)$', lambda m: comments + m[0], netlist, flags=re.I | re.M)
                verification = dest / 'verification.cir'
                verification.write_text(netlist)
                try:
                    with (dest / 'verification.log').open('x') as log:
                        verified = subprocess.run(['perl', str(args.regression / 'TestScripts/xyce_verify.pl'),
                            '--verbose', str(verification), str(ref_output), str(test_output)],
                            stdout=log, stderr=subprocess.STDOUT, timeout=60)
                    row['upstream_verifier'] = verified.returncode
                except subprocess.TimeoutExpired:
                    row['upstream_verifier'] = 'timeout'
            row['qualified'] &= reference['qualified'] and row['waveform']['status'] == 'pass' and all(
                row[k] == reference[k] for k in ('successful_steps', 'failed_steps')) and row['upstream_verifier'] == 0
        row['output_sha256'] = {str(path): sha(path) for path in
                                (dest/'stdout.log', dest/'input.cir.prn', dest/'verification.log') if path.exists()}
        rows.append(row)
        with (output / 'observations.jsonl').open('a') as f:
            f.write(json.dumps(row, allow_nan=False) + '\n')
        print(domain, Path(case).stem, count, trial, side, 'PASS' if row['qualified'] else 'FAIL', f'{wall:.3f}s', flush=True)
        return row

    for domain, cpus in domains:
        for index, case in enumerate(cases):
            reference = execute(domain, cpus, case, 1, -1, 'klu')
            for count in threads:
                for trial in range(args.passes):
                    for side in launch_order(sides, index, trial, args.order == 'confirmation'):
                        execute(domain, cpus, case, count, trial, side, reference)
    changed = [p for p, h in hashes.items() if not Path(p).exists() or sha(p) != h]
    summary = reduce_results(metadata, rows)
    summary['changed_inputs'] = changed
    if changed:
        summary['all_qualified'] = summary['incremental_gate'] = summary['all56_peer_wins'] = False
    write_json(output / 'summary.json', summary)
    return 0 if summary['all_qualified'] else 1


def review(output):
    """Rebuild the complete decision from retained logs, waveforms and pins."""
    metadata = json.loads((output / 'metadata.json').read_text())
    rows = [json.loads(line) for line in (output / 'observations.jsonl').read_text().splitlines()]
    references = {(r['domain'], r['case']): r for r in rows if r['side'] == 'klu'}
    for row in rows:
        dest = Path(row['directory'])
        parsed = parse_log((dest / 'stdout.log').read_text())
        if any(row.get(k) != v for k, v in parsed.items()):
            row['qualified'] = False
        if sha(dest / 'input.cir') != row['netlist_sha256']:
            row['qualified'] = False
        if any(not Path(path).exists() or sha(path)!=h for path,h in row['output_sha256'].items()):
            row['qualified'] = False
        if row['exit_code'] != 0 or row['timed_out'] or not row['closure']['all_session_processes_absent']:
            row['qualified'] = False
        if row['side'] != 'klu':
            ref = references.get((row['domain'], row['case']))
            if ref is None:
                row['qualified'] = False
            else:
                comparison = compare_waveforms(Path(ref['directory'])/'input.cir.prn', dest/'input.cir.prn')
                row['qualified'] &= comparison['status'] == 'pass' and row.get('upstream_verifier') == 0 and all(
                    row[k] == ref[k] for k in ('successful_steps', 'failed_steps'))
    summary = reduce_results(metadata, rows)
    expected_order = []
    for domain, _ in metadata['domains']:
        for index, case in enumerate(metadata['cases']):
            expected_order.append((domain, case, 1, -1, 'klu'))
            for count in metadata['threads']:
                for trial in range(metadata['passes']):
                    for side in launch_order(metadata['sides'], index, trial, metadata['order']=='confirmation'):
                        expected_order.append((domain, case, count, trial, side))
    actual_order = [(r['domain'],r['case'],r['threads'],r['trial'],r['side']) for r in rows]
    summary['launch_order_verified'] = actual_order == expected_order
    changed = [p for p,h in metadata['input_sha256'].items() if not Path(p).exists() or sha(p)!=h]
    summary['changed_inputs'] = changed
    if changed or not summary['launch_order_verified']:
        summary['all_qualified'] = summary['incremental_gate'] = summary['all56_peer_wins'] = False
    write_json(output/'review.json', summary)
    return 0 if summary['all_qualified'] else 1


def main():
    import sys
    if len(sys.argv) == 3 and sys.argv[1] == '--review':
        return review(Path(sys.argv[2]).resolve())
    if len(sys.argv) == 4 and sys.argv[1] == '--confirm':
        outputs = [Path(x).resolve() for x in sys.argv[2:]]
        for output in outputs:
            review(output)
        metadata = [json.loads((x/'metadata.json').read_text()) for x in outputs]
        summaries = [json.loads((x/'review.json').read_text()) for x in outputs]
        same_binaries = all(sha(outputs[0]/'binaries'/side) == sha(outputs[1]/'binaries'/side)
                            for side in ('baseline', 'candidate'))
        def common_inputs(m, output):
            return {p:h for p,h in m['input_sha256'].items()
                    if not Path(p).is_relative_to(output/'binaries')}
        matching = (outputs[0] != outputs[1] and same_binaries and
                    [m['order'] for m in metadata] == ['initial', 'confirmation'] and
                    common_inputs(metadata[0], outputs[0]) == common_inputs(metadata[1], outputs[1]))
        result = dict(campaigns=list(map(str, outputs)), matching_inputs=matching,
                      incremental_gate=matching and all(s['incremental_gate'] for s in summaries),
                      all56_peer_wins=matching and all(s['all56_peer_wins'] for s in summaries))
        write_json(outputs[1]/'confirmation.json', result)
        print(json.dumps(result, indent=2))
        return 0 if result['incremental_gate'] else 1
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--workspace', type=Path, required=True)
    p.add_argument('--regression', type=Path, required=True)
    p.add_argument('--image', type=Path, required=True)
    p.add_argument('--license', type=Path, required=True)
    p.add_argument('--baseline', type=Path, required=True)
    p.add_argument('--candidate', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--peer-library', type=Path, action='append', default=[])
    p.add_argument('--provenance', type=Path)
    p.add_argument('--peers', action='store_true')
    p.add_argument('--passes', type=int, default=4)
    p.add_argument('--order', choices=['initial', 'confirmation'], default='initial')
    p.add_argument('--smoke', action='store_true', help='Incomplete mux8/one-thread diagnostic only')
    p.add_argument('--diagnostic', action='store_true', help='Full-grid correctness run; all timings excluded')
    args = p.parse_args()
    if args.passes < 1 or (not args.smoke and not args.diagnostic and args.passes % 4):
        p.error('complete campaigns require a positive multiple of four passes')
    return run(args)


if __name__ == '__main__':
    raise SystemExit(main())
