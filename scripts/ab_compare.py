#!/usr/bin/env python3
"""Compare two kls_bench jsonl sides (A=pre, B=post) from ab_kls_pair.sh.

Per-matrix minimum cycle across passes per side; prints per-row cycle,
steady-refactor and solve deltas plus path changes.
Usage: ab_compare.py PRE.jsonl POST.jsonl
"""
import math, sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from score_paper_suite import load, gm


def main():
    pre = load(sys.argv[1])
    post = load(sys.argv[2])
    rows = []
    for m in sorted(pre):
        a, b = pre.get(m), post.get(m)
        if a is None or b is None:
            print(f"  {m}: SKIP (pre={'ok' if a else 'fail'} post={'ok' if b else 'fail'})")
            continue
        rows.append((b['cycle'] / a['cycle'], m, a, b))
    rows.sort()
    print(f"\npost/pre cycle gm={gm([r[0] for r in rows]):.4f} "
          f"improved {sum(1 for r in rows if r[0] < 0.98)} "
          f"regressed {sum(1 for r in rows if r[0] > 1.02)} of {len(rows)}")
    print(f"{'matrix':16} {'cyc':>6} {'rs':>6} {'solve':>6}  pre-path -> post-path")
    for rat, m, a, b in rows:
        print(f"{m:16} {rat:6.3f} {b['rs']/max(a['rs'],1e-12):6.2f} "
              f"{b['solve']/max(a['solve'],1e-12):6.2f}  "
              f"{a['path']} -> {b['path']}"
              f"{'   resid ' + format(b['resid'], '.1e') if b['resid'] > 1e-8 else ''}")


if __name__ == '__main__':
    main()
