#!/usr/bin/env python3
"""Draw a fresh sealed SPICE set from decks no KLS decision has seen.

The first held-out manifest (`bench/spice_holdout_manifest.json`) has been scored
at many milestones, so its holdout half can no longer be treated as unseen.  This
script draws a second set from regression decks the first screening could not use:

* decks that read other files (.include, .lib, file=), staged with the files
  beside them;
* decks that ran but were too short to time, with the transient's final time
  multiplied by a whole factor chosen from the KLU reference time alone.

Decks already in the first manifest, the tuning circuits, and every deck in a
directory that holds a tuning or dev deck are excluded.  Membership is decided from
KLU reference runs only, exactly as before; every case is `holdout` (sealed), and
this script prints counts, never names.
"""
import argparse
import concurrent.futures
import hashlib
import json
import math
from pathlib import Path
import subprocess

from build_spice_holdout_manifest import (STABILITY_FACTORS, STABILITY_STEP_RATIO, candidate, screen,
                                          stability)
from run_spice_campaign import CASES

MAX_TRAN_SCALE = 64


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for name in ('workspace', 'regression', 'image', 'xyce', 'license', 'work', 'output', 'seen_manifest',
                 'plain_rows', 'rows'):
        p.add_argument('--' + name.replace('_', '-'), type=Path, required=True)
    p.add_argument('--min-linear-seconds', type=float, default=0.05)
    p.add_argument('--extend-from', type=float, default=0.01,
                   help='Shortest KLU linear time a deck may have and still be lengthened')
    p.add_argument('--timeout', type=float, default=300)
    p.add_argument('--per-directory', type=int, default=2)
    p.add_argument('--jobs', type=int, default=4)
    p.add_argument('--cpus', help='taskset list for the reference runs')
    args = p.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    root = args.regression / 'Netlists'
    seen = json.loads(args.seen_manifest.read_text())
    used = set(CASES) | {c['case'] for c in seen['cases']}
    inspected = {str(Path(c).parent) for c in CASES} | {
        str(Path(c['case']).parent) for c in seen['cases'] if c['split'] == 'dev'}

    def fresh(case):
        return case not in used and str(Path(case).parent) not in inspected

    rows = json.loads(args.rows.read_text()) if args.rows.exists() else {}

    def save():
        args.rows.write_text(json.dumps(rows, indent=1) + '\n')

    def run_all(keys, function):
        todo = [k for k in keys if k not in rows]
        with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
            for key, result in zip(todo, pool.map(function, todo)):
                rows[key] = result
        save()

    # First look: plain decks were measured by the first screening; companion decks are new.
    first = {r['case']: dict(r, companions=False) for r in json.loads(args.plain_rows.read_text())
             if fresh(r['case'])}
    companion_cases = sorted(str(path.relative_to(root)) for path in root.rglob('*.cir')
                             if candidate(path, companions=True) and fresh(str(path.relative_to(root))))
    run_all([f'first|{c}' for c in companion_cases],
            lambda key: screen(args, key.split('|', 1)[1], companions=True))
    first.update({c: dict(rows[f'first|{c}'], companions=True) for c in companion_cases})

    # Decide the transient length from the reference time, then measure that deck.
    plan = {}
    for case, r in first.items():
        if r['status'] != 'ran' or r['ill_posed'] or r['klu_linear_seconds'] < args.extend_from:
            continue
        scale = 1 if r['klu_linear_seconds'] >= args.min_linear_seconds else math.ceil(
            2 * args.min_linear_seconds / r['klu_linear_seconds'])
        if scale <= MAX_TRAN_SCALE:
            plan[case] = dict(companions=r['companions'], tran_scale=scale)
    scaled = [c for c, s in plan.items() if s['tran_scale'] != 1]
    run_all([f'scaled|{c}' for c in scaled], lambda key: screen(args, key.split('|', 1)[1],
                                                                **plan[key.split('|', 1)[1]]))
    timed = {}
    for case, staging in plan.items():
        r = rows[f'scaled|{case}'] if staging['tran_scale'] != 1 else first[case]
        if r['status'] == 'ran' and not r['ill_posed'] and r['klu_linear_seconds'] >= args.min_linear_seconds:
            timed[case] = dict(r, **staging)
    run_all([f'stable|{c}' for c in timed], lambda key: stability(
        args, key.split('|', 1)[1], companions=timed[key.split('|', 1)[1]]['companions'],
        tran_scale=timed[key.split('|', 1)[1]]['tran_scale']))
    eligible = [dict(r, **rows[f'stable|{c}']) for c, r in timed.items()]
    eligible = [r for r in eligible if not r['unstable']]
    eligible.sort(key=lambda r: hashlib.sha256(r['case'].encode()).hexdigest())
    kept, per_directory, contents = [], {}, set()
    for r in eligible:
        directory = str(Path(r['case']).parent)
        if r['sha256'] in contents or per_directory.get(directory, 0) >= args.per_directory:
            continue
        contents.add(r['sha256'])
        per_directory[directory] = per_directory.get(directory, 0) + 1
        kept.append(dict(r, split='holdout'))
    kept.sort(key=lambda r: r['case'])
    manifest = dict(
        schema=2,
        purpose='Second sealed SPICE set; membership decided from KLU reference runs only',
        regression_commit=subprocess.check_output(
            ['git', '-C', str(args.regression), 'rev-parse', 'HEAD'], text=True).strip(),
        tuning_cases=list(CASES), seen_manifest=args.seen_manifest.name,
        min_linear_seconds=args.min_linear_seconds, extend_from=args.extend_from,
        max_tran_scale=MAX_TRAN_SCALE, per_directory=args.per_directory,
        stability=dict(reltol_factors=STABILITY_FACTORS, step_ratio=STABILITY_STEP_RATIO),
        considered=len(first), measured=len(plan), timed=len(timed), eligible=len(eligible), cases=kept)
    args.output.write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(dict(considered=len(first), measured=len(plan), timed=len(timed),
                          eligible=len(eligible), kept=len(kept),
                          with_companions=sum(r['companions'] for r in kept),
                          lengthened=sum(r['tran_scale'] != 1 for r in kept),
                          directories=len(per_directory)), indent=2))


if __name__ == '__main__':
    main()
