#!/usr/bin/env python3
"""Time KLS against the peers on generated circuit families and check invariants.

The tuning circuits and the held-out decks are single points.  A family is the same
circuit structure at doubling sizes, so it shows what single points cannot:

  * where along the size axis KLS loses to the best peer configuration,
  * size cliffs: between neighbouring sizes KLS's per-solve time may not grow more
    than CLIFF times faster than the KLU reference's does,
  * thread scaling: N threads may not be slower than one (runner's rule),
  * start-up cost: selected sizes also run with the transient lengthened, so a
    policy that only pays off in long runs, or only looks good in short ones, shows.

Decks come from `make_spice_family.py` (deterministic) and from the hand-written
tiny circuits in `bench/spice_tiny`.  Timing, peer-best scoring and correctness
checks are `run_spice_generalization.py`'s; this script only prepares its input and
reads its output.  Nothing here is sealed: every number may be inspected.
"""
import argparse
import collections
import json
import math
from pathlib import Path
import shutil
import statistics
import subprocess
import sys

from make_spice_family import FAMILIES
from run_spice_campaign import scale_tran

CLIFF = 1.3
GATED_MODE = 'kls-app'


def prepare(args):
    decks = args.output / 'decks'
    decks.mkdir(parents=True)
    cases = []
    for family in args.family or sorted(FAMILIES):
        for exponent in range(args.min_log2, args.max_log2 + 1):
            text = FAMILIES[family](2 ** exponent, None)
            name = f'{family}_{exponent:02d}.cir'
            (decks / name).write_text(text)
            cases.append(dict(case=name, split='dev', family=family, log2_size=exponent, length=1))
            if exponent in args.long_log2 and args.long_scale > 1:
                name = f'{family}_{exponent:02d}_x{args.long_scale}.cir'
                (decks / name).write_text(scale_tran(text, args.long_scale))
                cases.append(dict(case=name, split='dev', family=family, log2_size=exponent,
                                  length=args.long_scale))
    if args.tiny and not args.family:
        for path in sorted(args.tiny.glob('*.cir')):
            shutil.copyfile(path, decks / path.name)
            cases.append(dict(case=path.name, split='dev', family='tiny', log2_size=None, length=1))
    manifest = args.output / 'family_manifest.json'
    manifest.write_text(json.dumps(dict(schema=1, cases=cases), indent=1) + '\n')
    return decks, manifest, cases


def geomean(values):
    return math.exp(sum(math.log(v) for v in values) / len(values))


def report(args, cases):
    run = args.output / 'run'
    detail = json.loads((run / 'cells.json').read_text())
    info = {c['case']: c for c in cases}
    reference, kls = {}, collections.defaultdict(list)
    for line in (run / 'observations.jsonl').read_text().splitlines():
        row = json.loads(line)
        if not row.get('per_solve'):
            continue
        if row['side'] == 'klu':
            reference[row['case']] = row['per_solve']
        elif row['side'] == GATED_MODE and (row.get('qualified') or row.get('tolerant')):
            kls[(row['case'], row['threads'])].append(row['per_solve'])
    kls = {k: statistics.median(v) for k, v in kls.items()}
    cells = [c for c in detail['cells'] if c['mode'] == GATED_MODE and c['ratio'] is not None]
    groups = collections.defaultdict(list)
    for c in cells:
        groups[(info[c['case']]['family'], c['peer'], c['threads'])].append(c['ratio'])
    families = {f'{family}/{peer}-t{threads}': dict(count=len(v), geomean=geomean(v), worst=max(v),
                                                     losses=sum(x > 1 for x in v))
                for (family, peer, threads), v in sorted(groups.items())}
    losses = sorted((dict(case=c['case'], threads=c['threads'], peer_best=c['peer_best'], ratio=c['ratio'])
                     for c in cells if c['ratio'] > 1), key=lambda c: -c['ratio'])
    cliffs = []
    by_family = collections.defaultdict(list)
    for c in cases:
        if c['length'] == 1 and c['log2_size'] is not None:
            by_family[c['family']].append(c)
    threads = sorted({t for _, t in kls})
    for family, members in by_family.items():
        members.sort(key=lambda c: c['log2_size'])
        for small, large in zip(members, members[1:]):
            for count in threads:
                a, b = kls.get((small['case'], count)), kls.get((large['case'], count))
                ra, rb = reference.get(small['case']), reference.get(large['case'])
                if a and b and ra and rb:
                    growth = (b / a) / (rb / ra)
                    if growth > CLIFF:
                        cliffs.append(dict(family=family, threads=count, small=small['case'],
                                           large=large['case'], growth_over_reference=growth))
    slower = [s for s in detail['thread_scaling'] if s['mode'] == GATED_MODE and s['ratio'] > 1.03]
    result = dict(candidate=str(args.candidate), cases=len(cases), scored_cells=len(cells),
                  kls_failures=len(detail['kls_failures']),
                  beats_peers=not losses and not detail['kls_failures'], no_cliff=not cliffs,
                  thread_monotone=not slower, families=families, losses=losses, cliffs=cliffs,
                  slower_with_threads=sorted(slower, key=lambda s: -s['ratio']))
    (args.output / 'family_report.json').write_text(json.dumps(result, indent=1) + '\n')
    for name, g in families.items():
        print(f"{name:24s} n={g['count']:2d} geomean {g['geomean']:.3f} worst {g['worst']:.3f} losses {g['losses']}")
    for c in losses:
        print(f"loss  {c['case']:28s} t{c['threads']} {c['peer_best']:12s} {c['ratio']:.3f}")
    for c in cliffs:
        print(f"cliff {c['small']} -> {c['large']} t{c['threads']} x{c['growth_over_reference']:.2f} of reference growth")
    for s in result['slower_with_threads']:
        print(f"threads {s['case']:28s} t{s['threads']} {s['ratio']:.3f} of 1 thread")
    print(json.dumps({k: result[k] for k in ('beats_peers', 'no_cliff', 'thread_monotone', 'kls_failures')}))
    return 0 if result['beats_peers'] and result['no_cliff'] and result['thread_monotone'] else 1


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for name in ('workspace', 'regression', 'image', 'license', 'candidate', 'output'):
        p.add_argument('--' + name, type=Path, required=True)
    p.add_argument('--baseline', type=Path)
    p.add_argument('--family', action='append', choices=sorted(FAMILIES), default=[])
    p.add_argument('--min-log2', type=int, default=6)
    p.add_argument('--max-log2', type=int, default=13)
    p.add_argument('--long-log2', type=lambda s: [int(x) for x in s.split(',') if x], default=[8, 12],
                   help='Sizes that also run with a lengthened transient')
    p.add_argument('--long-scale', type=int, default=8)
    p.add_argument('--tiny', type=Path, default=Path(__file__).resolve().parents[1] / 'bench/spice_tiny')
    p.add_argument('--passes', type=int, default=2)
    p.add_argument('--threads', default='1,8')
    p.add_argument('--domain', default='llc96')
    p.add_argument('--order', choices=['initial', 'confirmation'], default='initial')
    p.add_argument('--timeout', type=float, default=300)
    p.add_argument('--report-only', action='store_true', help='Re-read a finished run')
    args = p.parse_args()
    if args.report_only:
        cases = json.loads((args.output / 'family_manifest.json').read_text())['cases']
        return report(args, cases)
    args.output.mkdir(parents=True)
    decks, manifest, cases = prepare(args)
    cmd = [sys.executable, str(Path(__file__).with_name('run_spice_generalization.py')),
           '--workspace', str(args.workspace), '--regression', str(args.regression),
           '--image', str(args.image), '--license', str(args.license), '--candidate', str(args.candidate),
           '--output', str(args.output / 'run'), '--manifest', str(manifest), '--case-root', str(decks),
           '--set', 'dev', '--passes', str(args.passes), '--threads', args.threads, '--domain', args.domain,
           '--order', args.order, '--timeout', str(args.timeout)]
    if args.baseline:
        cmd += ['--baseline', str(args.baseline)]
    with (args.output / 'run.log').open('w') as log:
        subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    if not (args.output / 'run' / 'cells.json').exists():
        raise SystemExit(f'family run did not finish; see {args.output / "run.log"}')
    return report(args, cases)


if __name__ == '__main__':
    sys.exit(main())
