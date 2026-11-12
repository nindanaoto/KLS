#!/usr/bin/env python3
"""Per-phase linear-solver time for KLS and both peers on SPICE decks.

Diagnostic only: the adapters' profile timers are on, so these numbers must not
be used for wins or gates.  The output answers where time goes per solver call
(import of values, analysis, factorization, solve) and how much fill each
solver's ordering produces, which is what decides the order of solver work.
"""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess

from run_spice_campaign import CASES, DOMAINS, parse_log, stage

SOLVERS = {'kls': ('KLS', 'KLS_PROFILE=1'), 'ck': ('CKTSO', 'CKTSO_PROFILE=1'),
           'st': ('SUBTREELU', 'SUBTREELU_PROFILE=1')}
PROFILE = re.compile(r'^(?:KLS|CKTSO|SubtreeLU) profile (.+?) = ([-\d.eE+]+)\s*$', re.M)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('workspace', 'regression', 'image', 'license', 'xyce', 'output'):
        p.add_argument('--' + name, type=Path, required=True)
    p.add_argument('--manifest', type=Path)
    p.add_argument('--split', default='dev')
    p.add_argument('--threads', type=lambda s: [int(x) for x in s.split(',')], default=[1, 8])
    p.add_argument('--cpus', default=dict(DOMAINS)['llc96'])
    p.add_argument('--timeout', type=float, default=1800)
    args = p.parse_args()
    cases = list(CASES)
    if args.manifest:
        cases = [c['case'] for c in json.loads(args.manifest.read_text())['cases'] if c['split'] == args.split]
    args.output.mkdir(parents=True)
    env = {k: v for k, v in os.environ.items() if not k.startswith(('LD_', 'KLS_', 'OMP_'))}
    env.update(OMP_NUM_THREADS='1', XYCE_NO_TRACKING='1')
    rows = []
    for case in cases:
        for threads in args.threads:
            for side, (solver, extra) in SOLVERS.items():
                dest = args.output / Path(case).stem / f't{threads}-{side}'
                stage(args.regression / 'Netlists' / case, dest, solver, threads, args.license, extra)
                cmd = ['taskset', '-c', args.cpus, 'apptainer', 'exec', '--bind', f'{args.workspace}:/work',
                       '--bind', f'{args.workspace}:{args.workspace}', '--pwd', str(dest),
                       str(args.image), str(args.xyce), 'input.cir']
                with (dest / 'stdout.log').open('w') as stream:
                    try:
                        code = subprocess.run(cmd, stdout=stream, stderr=subprocess.STDOUT, env=env,
                                              timeout=args.timeout).returncode
                    except subprocess.TimeoutExpired:
                        code = 'timeout'
                text = (dest / 'stdout.log').read_text(errors='replace')
                log = parse_log(text)
                row = dict(case=Path(case).stem, threads=threads, solver=side, exit_code=code,
                           linear_seconds=log['linear_seconds'][-1] if log['linear_seconds'] else None,
                           profile={k: float(v) for k, v in PROFILE.findall(text)})
                rows.append(row)
                calls = row['profile'].get('calls') or 0
                per_call = {k[:-len(' seconds')]: 1e6 * v / calls for k, v in row['profile'].items()
                            if k.endswith(' seconds') and calls}
                print(row['case'], f't{threads}', side, code, row['linear_seconds'], f'calls={calls:.0f}',
                      ' '.join(f'{k}={v:.1f}us' for k, v in per_call.items()), flush=True)
    (args.output / 'profile.json').write_text(json.dumps(rows, indent=1) + '\n')


if __name__ == '__main__':
    main()
