#!/usr/bin/env python3
"""Score a paired paper-manifest run: per-row cycle + phase ratios vs CKTSO/SubtreeLU.

cycle = ana + init + solve + (rf + solve) + (H-2)*(rs + solve), H=100 solves.
Per-matrix minimum cycle across passes per side. Usage:
  score_paper_suite.py KLS.jsonl CK.jsonl [ST.jsonl]
"""
import json, math, sys
from collections import defaultdict

H = 100


def load(path):
    best = {}
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        d = json.loads(line)
        name = d.get('matrix', '')
        name = name.rsplit('/', 1)[-1].replace('.mtx', '').lower()
        if d.get('timeout_or_fail') or 'analysis_seconds' not in d:
            best.setdefault(name, None)
            continue
        ana = d['analysis_seconds']
        init = d.get('initial_factor_seconds', d.get('factor_seconds', 0.0))
        rf = d.get('refactor_first_seconds', d.get('refactor_seconds_avg', 0.0))
        rs = d.get('refactor_steady_seconds_avg', d.get('refactor_seconds_avg', 0.0))
        sv = d['solve_seconds_avg']
        cyc = ana + init + sv + (rf + sv) + (H - 2) * (rs + sv)
        rec = dict(cycle=cyc, ana=ana, init=init, rf=rf, rs=rs, solve=sv,
                   resid=d.get('relative_residual_l2', d.get('max_relative_residual', 0.0)),
                   path=d.get('last_refactor_path', ''), ordering=d.get('ordering', ''))
        cur = best.get(name)
        if cur is None or cyc < cur['cycle']:
            best[name] = rec
    return best


def gm(xs):
    xs = [x for x in xs if x and x > 0]
    return math.exp(sum(math.log(x) for x in xs) / len(xs)) if xs else float('nan')


def main():
    kls = load(sys.argv[1])
    comps = {'CK': load(sys.argv[2])}
    if len(sys.argv) > 3:
        comps['ST'] = load(sys.argv[3])
    for cname, comp in comps.items():
        rows = []
        for m, k in sorted(kls.items()):
            c = comp.get(m)
            if k is None or c is None:
                print(f"  {m}: SKIP (kls={'ok' if k else 'fail'} {cname.lower()}={'ok' if c else 'fail'})")
                continue
            rows.append((k['cycle'] / c['cycle'], m, k, c))
        rows.sort(reverse=True)
        print(f"\n== KLS vs {cname}: cycle gm={gm([r[0] for r in rows]):.3f} "
              f"wins {sum(1 for r in rows if r[0] < 1)}/{len(rows)} ==")
        print(f"{'matrix':16} {'cyc':>6} {'rs':>6} {'solve':>6} {'ana+init':>8}  path/ord      resid")
        for rat, m, k, c in rows:
            fe = (k['ana'] + k['init']) / max(c['ana'] + c['init'], 1e-12)
            print(f"{m:16} {rat:6.2f} {k['rs']/max(c['rs'],1e-12):6.2f} "
                  f"{k['solve']/max(c['solve'],1e-12):6.2f} {fe:8.2f}  "
                  f"{k['path']}/{k['ordering']:6} {k['resid']:.1e}")


if __name__ == '__main__':
    main()
