#!/usr/bin/env python3
"""Independently reconstruct a full SPICE campaign from its retained artifacts.

Does not import the launcher, its parser, its reducer, or its qualified flags.
Run between timing campaigns, never alongside one.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import statistics
import subprocess

CASES = {
    'PDE_2D_BJT/invertbjt.cir', 'MOS13_IC/invert50_mos1.cir',
    'CircuitSim90/MOS2_LARGE/mem_plus.cir', 'CircuitSim90/MOS2/mux8.cir',
    'CircuitSim90/MOS2_LARGE/dac.cir', 'NEURON/rallpack3.cir',
    'CircuitSim90/MOS2_LARGE/ram2k.cir',
}


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def output_shape(path):
    count = 0
    width = None
    with path.open() as f:
        header = next(f).split()
        for line in f:
            fields = line.split()
            if not fields or not fields[0].isdigit():
                continue
            if int(fields[0]) != count or len(fields) < 3:
                raise ValueError('nonsequential or empty output')
            if width is None:
                width = len(fields)
            if width != len(fields) or not all(math.isfinite(float(v)) for v in fields[1:]):
                raise ValueError('nonfinite or inconsistent output')
            count += 1
    if not count or len(header) != width:
        raise ValueError('empty output or inconsistent header')
    return count, header


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True)
    metadata_path = args.campaign / 'metadata.json'
    ledger_path = args.campaign / 'observations.jsonl'
    metadata = json.loads(metadata_path.read_text())
    ledger = [json.loads(line) for line in ledger_path.read_text().splitlines()]
    pins = dict(metadata['input_sha256'])
    pins.update({str(p.resolve()): digest(p) for p in (Path(__file__), metadata_path, ledger_path)})
    verifier = [Path(p) for p in pins if Path(p).name == 'xyce_verify.pl']
    assert len(verifier) == 1, 'ambiguous verifier'
    assert set(metadata['cases']) == CASES and len(metadata['cases']) == 7
    assert set(map(tuple, metadata['domains'])) == {('llc96', '0-7'), ('llc32', '8-15')}
    assert metadata['threads'] == [1, 8]
    sides = metadata['sides']
    assert sides in (['baseline', 'candidate'], ['baseline', 'candidate', 'ck', 'st'])
    expected_order = []
    for domain, _ in metadata['domains']:
        for index, case in enumerate(metadata['cases']):
            expected_order.append((domain, case, 1, -1, 'klu'))
            for threads in (1, 8):
                for trial in range(metadata['passes']):
                    offset = (index + trial + (metadata['order'] == 'initial')) % len(sides)
                    for side in sides[offset:] + sides[:offset]:
                        expected_order.append((domain, case, threads, trial, side))
    actual_order = [(r['domain'], r['case'], r['threads'], r['trial'], r['side']) for r in ledger]
    errors = []
    if actual_order != expected_order:
        errors.append('missing, duplicate, unexpected, or out-of-order observations')
    observations = {}
    references = {}
    for index, row in enumerate(ledger):
        key = actual_order[index]
        folder = Path(row['directory'])
        log_path = folder / 'stdout.log'
        text = log_path.read_text()
        times = [float(v) for v in re.findall(r'Total Linear Solution Time:\s*([\d.eE+-]+)', text)]
        steps = tuple(tuple(int(v) for v in re.findall(pattern, text)) for pattern in (
            r'Number Successful Steps Taken:\s*(\d+)', r'Number Failed Steps Attempted:\s*(\d+)'))
        row_errors = []
        expected_cpus = dict(metadata['domains'])[row['domain']]
        expected_binary = args.campaign.resolve() / 'binaries' / ('candidate' if row['side'] == 'candidate' else 'baseline')
        if (row['cpus'] != expected_cpus or row['command'][:3] != ['taskset', '-c', expected_cpus] or
                Path(row['command'][-2]).resolve() != expected_binary or row['command'][-1] != 'input.cir'):
            row_errors.append('wrong affinity, executable, or input command')
        if (row['exit_code'] != 0 or row['timed_out'] or
                not row['closure']['all_session_processes_absent'] or
                'End of Xyce(TM) Simulation' not in text or re.search('linear solve failed', text, re.I) or
                not times or not all(math.isfinite(v) and v > 0 for v in times) or not all(steps) or
                not math.isfinite(row['wall_seconds']) or row['wall_seconds'] <= 0):
            row_errors.append('execution or statistics failure')
        for path, expected_hash in row['output_sha256'].items():
            if not Path(path).exists() or digest(path) != expected_hash:
                row_errors.append('changed output: ' + path)
        if digest(folder / 'input.cir') != row['netlist_sha256']:
            row_errors.append('changed staged input')
        solver = {'klu':'KLU', 'baseline':'KLS', 'candidate':'KLS', 'ck':'CKTSO', 'st':'SUBTREELU'}[row['side']]
        options = '\n'.join(re.findall(r'^\s*\.options\s+linsol\b[^\n]*',
                                      (folder/'input.cir').read_text(), re.I | re.M))
        types = re.findall(r'\btype\s*=\s*(\w+)', options, re.I)
        if not types or any(t.upper() != solver for t in types):
            row_errors.append('wrong solver selection')
        if row['side'] != 'klu':
            counts = re.findall(r'\b' + solver + r'_THREADS\s*=\s*(\d+)', options, re.I)
            if not counts or any(int(t) != row['threads'] for t in counts):
                row_errors.append('wrong solver thread setting')
            if row.get('upstream_verifier') != 0:
                row_errors.append('original upstream verifier failed; a rerun cannot redeem it')
        try:
            shape = output_shape(folder / 'input.cir.prn')
            if row['side'] == 'klu':
                references[key[:2]] = (folder, steps, shape, not row_errors)
            else:
                reference, ref_steps, ref_shape, ref_valid = references[key[:2]]
                if not ref_valid or steps != ref_steps or shape != ref_shape:
                    row_errors.append('reference, step, row, or signal mismatch')
                source = re.sub(r'^\s*\*comp[^\n]*\n', '', (folder / 'input.cir').read_text(), flags=re.I | re.M)
                comparisons = ''.join(f'*comp {name} reltol=1e-8 abstol=1e-10 zerotol=0 absdifftol=0 numfail=0\n'
                                      for name in shape[1][2:])
                source, count = re.subn(r'^(\s*\.end\s*)$', lambda m: comparisons + m[0], source, flags=re.I | re.M)
                assert count == 1
                verification = args.output / f'{index:04d}.cir'
                verification.write_text(source)
                with (args.output / f'{index:04d}.log').open('x') as stream:
                    checked = subprocess.run(['perl', str(verifier[0]), '--verbose', str(verification),
                        str(reference / 'input.cir.prn'), str(folder / 'input.cir.prn')],
                        stdout=stream, stderr=subprocess.STDOUT, timeout=60)
                if checked.returncode != 0:
                    row_errors.append('fresh upstream verification failed')
        except (ValueError, KeyError, FileNotFoundError, StopIteration, subprocess.TimeoutExpired) as exc:
            row_errors.append(repr(exc))
        observations[key] = dict(valid=not row_errors, linear=times[-1] if times else None, wall=row['wall_seconds'])
        if row_errors:
            errors.append(dict(observation=index, key=key, errors=row_errors))
    changed = [p for p, h in pins.items() if not Path(p).exists() or digest(p) != h]
    if changed:
        errors.append(dict(changed_inputs=changed))
    cells = []
    if not errors:
        for domain, _ in metadata['domains']:
            for case in sorted(CASES):
                for threads in (1, 8):
                    for peer in sides:
                        if peer == 'candidate':
                            continue
                        pairs = [(observations[(domain, case, threads, t, 'candidate')],
                                  observations[(domain, case, threads, t, peer)]) for t in range(metadata['passes'])]
                        cells.append(dict(domain=domain, case=case, threads=threads, peer=peer,
                            median=statistics.median(a['linear']/b['linear'] for a,b in pairs),
                            wall_median=statistics.median(a['wall']/b['wall'] for a,b in pairs)))
    official = (not errors and not metadata['diagnostic_only'] and metadata['complete_grid'] and
                metadata['passes'] >= 4 and metadata['passes'] % 4 == 0)
    baseline = [c for c in cells if c['peer'] == 'baseline']
    peers = [c for c in cells if c['peer'] in ('ck', 'st')]
    stable = official and len(baseline) == 28 and all(c['median'] <= 1.02 and c['wall_median'] <= 1.02 for c in baseline)
    result = dict(all_qualified=not errors, observations=len(ledger), errors=errors, cells=cells,
                  incremental_gate=stable and any(c['median'] <= .98 for c in baseline),
                  all56_peer_wins=stable and len(peers) == 56 and all(c['median'] <= .98 for c in peers),
                  input_sha256=pins, separate_confirmation_required=True)
    (args.output / 'audit.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k:result[k] for k in ('all_qualified', 'observations', 'incremental_gate', 'all56_peer_wins')}))
    return 0 if not errors else 1


if __name__ == '__main__':
    raise SystemExit(main())
