#!/usr/bin/env python3
"""Screen Xyce regression decks into a held-out SPICE manifest.

The seven campaign circuits are the tuning set: every policy decision in KLS has
seen them.  This script finds other self-contained transient decks, measures
them once with the KLU reference solver only (so no KLS or peer timing can
influence membership), keeps those with enough linear-solver work to time, and
splits them into `dev` and `holdout` by a hash of the deck path.  `dev` may be
inspected freely; `holdout` is reported only in aggregate at milestones.
"""
import argparse
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

from run_spice_campaign import CASES, parse_log, sha, stage

# Analyses that change what is solved; post-processing (.measure, .four, .fft) is fine.
ANALYSES = r'step|sampling|embeddedsampling|pce|hb|mpde|ac|noise|dc|sens|lin|result'


def candidate(path):
    """Plain single-run transient decks that stage() can retarget."""
    try:
        text = path.read_text(errors='strict')
    except (UnicodeDecodeError, OSError):
        return False
    flags = re.I | re.M
    return bool(re.search(r'^\s*\.tran\b', text, flags) and
                re.search(r'^\s*\.print\s+tran\b', text, flags) and
                len(re.findall(r'^\s*\.end\s*$', text, flags)) == 1 and
                not re.search(r'^\s*\.(?:include|inc|lib)\b', text, flags) and
                not re.search(rf'^\s*\.(?:{ANALYSES})\b', text, flags) and
                not re.search(r'\bfile\s*=', text, re.I))


def split_of(case):
    return 'holdout' if hashlib.sha256(case.encode()).digest()[0] & 1 else 'dev'


def screen(args, case):
    dest = args.work / hashlib.sha256(case.encode()).hexdigest()[:16]
    row = dict(case=case, sha256=sha(args.regression / 'Netlists' / case))
    try:
        stage(args.regression / 'Netlists' / case, dest, 'KLU', 1, args.license)
    except (ValueError, OSError) as error:
        return dict(row, status='unstageable', detail=str(error))
    cmd = ['apptainer', 'exec', '--bind', f'{args.workspace}:/work',
           '--bind', f'{args.workspace}:{args.workspace}', '--pwd', str(dest),
           str(args.image), str(args.xyce), 'input.cir']
    env = {k: v for k, v in os.environ.items() if not k.startswith(('LD_', 'KLS_', 'OMP_'))}
    env.update(OMP_NUM_THREADS='1', XYCE_NO_TRACKING='1')
    start = time.monotonic()
    try:
        with (dest / 'stdout.log').open('w') as stream:
            code = subprocess.run(cmd, stdout=stream, stderr=subprocess.STDOUT, env=env,
                                  timeout=args.timeout).returncode
    except subprocess.TimeoutExpired:
        shutil.rmtree(dest, ignore_errors=True)
        return dict(row, status='timeout')
    wall = time.monotonic() - start
    log = parse_log((dest / 'stdout.log').read_text(errors='replace'))
    ok = (code == 0 and log['end_of_simulation'] and log['linear_seconds'] and
          log['successful_steps'] and not log['linear_failure_messages'] and
          (dest / 'input.cir.prn').exists())
    shutil.rmtree(dest, ignore_errors=True)
    if not ok:
        return dict(row, status='failed', exit_code=code)
    return dict(row, status='ran', wall_seconds=wall, klu_linear_seconds=log['linear_seconds'][-1],
                successful_steps=log['successful_steps'][-1], failed_steps=log['failed_steps'][-1])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('workspace', 'regression', 'image', 'xyce', 'license', 'work', 'output'):
        p.add_argument('--' + name, type=Path, required=True)
    p.add_argument('--min-linear-seconds', type=float, default=0.2)
    p.add_argument('--timeout', type=float, default=300)
    p.add_argument('--per-directory', type=int, default=3,
                   help='Cap per regression directory so one deck family cannot dominate')
    p.add_argument('--jobs', type=int, default=4)
    p.add_argument('--rows', type=Path, help='Reuse or record the raw KLU screening rows')
    args = p.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    root = args.regression / 'Netlists'
    cases = sorted(str(path.relative_to(root)) for path in root.rglob('*.cir')
                   if str(path.relative_to(root)) not in CASES and candidate(path))
    if args.rows and args.rows.exists():
        rows = json.loads(args.rows.read_text())
    else:
        with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
            rows = list(pool.map(lambda case: screen(args, case), cases))
        if args.rows:
            args.rows.write_text(json.dumps(rows, indent=1) + '\n')
    eligible = [r for r in rows if r['status'] == 'ran' and
                r['klu_linear_seconds'] >= args.min_linear_seconds]
    # Deterministic, timing-independent choice inside a directory.
    eligible.sort(key=lambda r: hashlib.sha256(r['case'].encode()).hexdigest())
    kept, seen, contents = [], {}, set()
    for r in eligible:
        directory = str(Path(r['case']).parent)
        if r['sha256'] in contents or seen.get(directory, 0) >= args.per_directory:
            continue
        contents.add(r['sha256'])
        seen[directory] = seen.get(directory, 0) + 1
        kept.append(dict(r, split=split_of(r['case'])))
    kept.sort(key=lambda r: r['case'])
    statuses = {}
    for r in rows:
        statuses[r['status']] = statuses.get(r['status'], 0) + 1
    manifest = dict(
        schema=1,
        purpose='SPICE decks never used to tune KLS; membership decided from KLU reference runs only',
        regression_commit=subprocess.check_output(
            ['git', '-C', str(args.regression), 'rev-parse', 'HEAD'], text=True).strip(),
        tuning_cases=list(CASES), min_linear_seconds=args.min_linear_seconds,
        per_directory=args.per_directory, screened=len(rows), screen_statuses=statuses,
        eligible=len(eligible), cases=kept)
    args.output.write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(dict(screened=len(rows), statuses=statuses, eligible=len(eligible),
                          kept=len(kept), dev=sum(r['split'] == 'dev' for r in kept),
                          holdout=sum(r['split'] == 'holdout' for r in kept)), indent=2))


if __name__ == '__main__':
    main()
