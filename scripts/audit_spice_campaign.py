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


def hybrid_acceptance(cells, official):
    # Independent implementation: never import the campaign reducer or gate.
    required={(domain,case,threads,peer) for domain in ('llc96','llc32')
              for case in CASES for threads in (1,8) for peer in ('baseline','ck','st')}
    by_key={(c['domain'],c['case'],c['threads'],c['peer']):c for c in cells}
    valid=official is True and len(cells)==len(by_key)==84 and set(by_key)==required
    valid=valid and all(type(c.get(k)) in (int,float) and math.isfinite(c[k]) and c[k]>0
                        for c in cells for k in ('median','wall_median'))
    if not valid:return dict(passed=False,complete=False,targets=[],regressions=[])
    targets=[by_key[d,'MOS13_IC/invert50_mos1.cir',1,'baseline'] for d in ('llc96','llc32')]
    regressions=[]
    for d,case,threads,peer in sorted(required):
        if peer!='baseline':continue
        row=by_key[d,case,threads,peer]
        if row['median']>1.02 or row['wall_median']>1.02:
            exception=by_key[d,case,threads,'ck']['median']<=.98 and by_key[d,case,threads,'st']['median']<=.98
            regressions.append(dict(domain=d,case=case,threads=threads,peer_exception=exception))
    return dict(passed=all(c['median']<=.98 for c in targets) and all(r['peer_exception'] for r in regressions),
                complete=True,targets=[{k:c[k] for k in ('domain','case','threads','median','wall_median')} for c in targets],regressions=regressions)


def peer_generalization_acceptance(cells, official, passes):
    # Reconstruct the new rule independently of the launcher and its reducer.
    keys={(d,case,t,p) for d in ('llc96','llc32') for case in CASES
          for t in (1,8) for p in ('baseline','ck','st')}
    records={(c['domain'],c['case'],c['threads'],c['peer']):c for c in cells}
    valid=(official is True and passes>=8 and passes%4==0 and
           len(records)==len(cells)==84 and set(records)==keys and
           all(c.get('complete',True) is True for c in cells))
    for key,c in records.items():
        fields=['median','wall_median']
        if key[-1]!='baseline':fields+=['baseline_median','baseline_wall_median']
        valid=valid and all(type(c.get(f)) in (int,float) and math.isfinite(c[f]) and c[f]>0 for f in fields)
    if not valid:return dict(passed=False,complete=False,protected_wins=0,lost_wins=[])
    controls=[c for key,c in records.items() if key[-1]=='baseline']
    linear=math.exp(sum(math.log(c['median']) for c in controls)/28)
    wall=math.exp(sum(math.log(c['wall_median']) for c in controls)/28)
    count=0;losses=[]
    for key in sorted(keys):
        if key[-1]=='baseline':continue
        c=records[key]
        for metric,base,new in [('linear','baseline_median','median'),('wall','baseline_wall_median','wall_median')]:
            if c[base]>=1:continue
            count+=1
            if c[new]>1:
                losses.append(dict(domain=key[0],case=key[1],threads=key[2],peer=key[3],metric=metric,
                                   reference_ratio=c[base],candidate_ratio=c[new]))
    return dict(passed=not losses and linear<=1 and wall<=1,complete=True,
                protected_wins=count,lost_wins=losses,linear_geomean=linear,wall_geomean=wall)


def confirmation_peer_wins(summaries):
    # Independent cross-campaign reconstruction, including incomplete cells.
    keys=('domain','case','threads','peer')
    records=[];protected=set()
    for summary in summaries:
        run={}
        for c in summary['cells']:
            if c['peer'] not in ('ck','st'):continue
            key=tuple(c[k] for k in keys);run[key]=c
            for metric,prior in (('linear','baseline_median'),('wall','baseline_wall_median')):
                if type(c.get(prior)) in (float,int) and 0<c[prior]<1:protected.add((key,metric))
        records.append(run)
    failures=[];absent=[]
    for i,run in enumerate(records):
        for key,metric in sorted(protected):
            c=run.get(key,{})
            value=c.get({'linear':'median','wall':'wall_median'}[metric])
            detail=dict(zip(keys,key),campaign=('first','second')[i],metric=metric)
            if type(value) not in (int,float) or not math.isfinite(value) or value<=0:absent.append(detail)
            elif value>1:failures.append(dict(detail,candidate_ratio=value))
    return dict(passed=not failures and not absent,protected_wins=len(protected),
                lost_peer_wins=failures,missing_peer_comparisons=absent)


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
                        if peer in ('ck','st'):
                            prior=[(observations[(domain,case,threads,t,'baseline')],
                                    observations[(domain,case,threads,t,peer)]) for t in range(metadata['passes'])]
                            cells[-1].update(baseline_median=statistics.median(a['linear']/b['linear'] for a,b in prior),
                                             baseline_wall_median=statistics.median(a['wall']/b['wall'] for a,b in prior))
    official = (not errors and not metadata['diagnostic_only'] and metadata['complete_grid'] and
                metadata['passes'] >= 4 and metadata['passes'] % 4 == 0)
    baseline = [c for c in cells if c['peer'] == 'baseline']
    peers = [c for c in cells if c['peer'] in ('ck', 'st')]
    stable = official and len(baseline) == 28 and all(c['median'] <= 1.02 and c['wall_median'] <= 1.02 for c in baseline)
    hybrid=hybrid_acceptance(cells,official)
    policy=metadata.get('acceptance_policy','incremental')
    incremental=stable and any(c['median']<=.98 for c in baseline)
    generalization=(stable and metadata['passes']>=8 and
      math.exp(sum(math.log(c['median']) for c in baseline)/len(baseline))<=1 and
      math.exp(sum(math.log(c['wall_median']) for c in baseline)/len(baseline))<=1)
    peer_generalization=peer_generalization_acceptance(cells,official,metadata['passes'])
    result = dict(acceptance_policy=policy,
                  acceptance_gate=(peer_generalization['passed'] if policy=='generalization-v2' else
                                   generalization if policy=='generalization-v1' else hybrid['passed'] if policy=='auto-btf' else incremental),
                  generalization_gate=generalization,
                  peer_generalization_gate=peer_generalization['passed'],peer_generalization=peer_generalization,
                  hybrid_gate=hybrid['passed'],hybrid=hybrid,
                  all_qualified=not errors, observations=len(ledger), errors=errors, cells=cells,
                  incremental_gate=stable and any(c['median'] <= .98 for c in baseline),
                  all56_peer_wins=stable and len(peers) == 56 and all(c['median'] <= .98 for c in peers),
                  input_sha256=pins, separate_confirmation_required=True)
    (args.output / 'audit.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k:result[k] for k in ('all_qualified', 'observations', 'incremental_gate', 'all56_peer_wins')}))
    return 0 if not errors else 1


if __name__ == '__main__':
    raise SystemExit(main())
