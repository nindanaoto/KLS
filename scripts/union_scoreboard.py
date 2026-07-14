#!/usr/bin/env python3
"""Union-goal scoreboard: beat BOTH CKTSO and SubtreeLU on every row.

Joins one or more (kls, ck, st) jsonl triples, computes each row's worst-side
T100 cycle ratio (CK everywhere, ST only where ST is residual-valid), and
tiers the rows:

  WON   worst < 1.0        (margin-hardened when worst <= 0.90)
  BAND  worst < 1.15
  NEAR  worst < 1.5
  MID   worst < 2.5
  DEEP  worst >= 2.5

Rows where KLS is invalid while a competitor is valid are BLOCKED
(automatic losses).  Rows where every solver fails are DISPOSITION rows.
Rows where KLS is valid and both competitors fail count as WON by default.

Usage:
  union_scoreboard.py --set med kls.jsonl ck.jsonl st.jsonl \
                      --set lrg kls_l.jsonl ck_l.jsonl st_l.jsonl \
                      [--selection median|min] [--residual-threshold 1e-8]
"""

import argparse
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from score_paper_suite import load, geometric_mean  # noqa: E402

TIERS = (
    ("WON", 1.0),
    ("BAND", 1.15),
    ("NEAR", 1.5),
    ("MID", 2.5),
    ("DEEP", float("inf")),
)

HORIZON = 100


def tier_of(worst):
    for name, ceiling in TIERS:
        if worst < ceiling:
            return name
    return "DEEP"


def dominant_term(kls, ref):
    """Which weighted term contributes most excess cycle time vs ref."""
    weights = (
        ("rs", HORIZON - 2),
        ("solve", HORIZON),
        ("front", 1.0),
    )
    kls_front = kls["ana"] + kls["init"] + kls["rf"]
    ref_front = ref["ana"] + ref["init"] + ref["rf"]
    excess = {
        "rs": (kls["rs"] - ref["rs"]) * (HORIZON - 2),
        "solve": (kls["solve"] - ref["solve"]) * HORIZON,
        "front": kls_front - ref_front,
    }
    return max(excess, key=excess.get)


def ratio(kls, ref, key):
    if key == "front":
        num = kls["ana"] + kls["init"]
        den = ref["ana"] + ref["init"]
    else:
        num, den = kls[key], ref[key]
    return num / max(den, 1e-12)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--set",
        dest="sets",
        nargs=4,
        action="append",
        metavar=("LABEL", "KLS", "CK", "ST"),
        required=True,
        help="one manifest group: label + kls/ck/st jsonl paths",
    )
    parser.add_argument("--residual-threshold", type=float, default=1e-8)
    parser.add_argument("--selection", choices=("median", "min"), default="median")
    parser.add_argument(
        "--won-margin",
        type=float,
        default=0.90,
        help="worst-side ratio at or under which a win counts as hardened",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    rows = {}          # name -> dict
    blocked = {}       # name -> reason
    dispositions = []  # all-solvers-fail
    default_wins = []  # kls valid, both competitors invalid

    for label, kls_path, ck_path, st_path in args.sets:
        kls, kls_fail = load(kls_path, args.residual_threshold, args.selection)
        ck, _ = load(ck_path, args.residual_threshold, args.selection)
        st, _ = load(st_path, args.residual_threshold, args.selection)
        for name in kls.keys() | ck.keys() | st.keys():
            k, c, s = kls.get(name), ck.get(name), st.get(name)
            if k is None:
                reason = ", ".join(kls_fail.get(name, ["missing"]))
                if c is None and s is None:
                    dispositions.append((name, label))
                else:
                    blocked[name] = (label, reason,
                                     "ck+st" if (c and s) else ("ck" if c else "st"))
                continue
            if c is None and s is None:
                default_wins.append((name, label))
                continue
            ck_ratio = k["cycle"] / c["cycle"] if c else None
            st_ratio = k["cycle"] / s["cycle"] if s else None
            candidates = [(r, ref, side) for r, ref, side in
                          ((ck_ratio, c, "CK"), (st_ratio, s, "ST")) if r is not None]
            worst, ref, side = max(candidates, key=lambda item: item[0])
            rows[name] = {
                "set": label,
                "ck": ck_ratio,
                "st": st_ratio,
                "worst": worst,
                "side": side,
                "rs": ratio(k, ref, "rs"),
                "solve": ratio(k, ref, "solve"),
                "front": ratio(k, ref, "front"),
                "dom": dominant_term(k, ref) if worst >= 1.0 else "",
                "path": f"{k['path']}/{k['ordering']}",
                "resid": k["resid"],
            }

    by_tier = {name: [] for name, _ in TIERS}
    for name, row in rows.items():
        by_tier[tier_of(row["worst"])].append((row["worst"], name, row))

    total = len(rows) + len(blocked) + len(dispositions) + len(default_wins)
    won = len(by_tier["WON"]) + len(default_wins)
    remaining = len(rows) - len(by_tier["WON"]) + len(blocked)
    print(f"UNION SCOREBOARD  rows={total}  won={won}  remaining={remaining}  "
          f"dispositions={len(dispositions)}  "
          f"(selection={args.selection}, resid<={args.residual_threshold:.0e})")
    ck_all = [r["ck"] for r in rows.values() if r["ck"]]
    st_all = [r["st"] for r in rows.values() if r["st"]]
    print(f"geomean vs CK {geometric_mean(ck_all):.3f} ({len(ck_all)} rows), "
          f"vs valid-ST {geometric_mean(st_all):.3f} ({len(st_all)} rows)\n")

    for tier_name, _ in TIERS:
        entries = sorted(by_tier[tier_name], reverse=True)
        if not entries:
            continue
        print(f"== {tier_name} ({len(entries)}) ==")
        for worst, name, row in entries:
            thin = ""
            if tier_name == "WON":
                thin = "  [thin]" if worst > args.won_margin else "  [hard]"
            ck_txt = f"{row['ck']:.2f}" if row["ck"] else "--"
            st_txt = f"{row['st']:.2f}" if row["st"] else "--"
            print(f"  {name:22s} {row['set']:3s} worst={worst:5.2f}({row['side']}) "
                  f"ck={ck_txt:>5} st={st_txt:>5} "
                  f"rs={row['rs']:5.2f} solve={row['solve']:5.2f} "
                  f"front={row['front']:5.2f} dom={row['dom']:5s} "
                  f"{row['path']} {row['resid']:.0e}{thin}")
        print()

    if default_wins:
        print(f"== WON by default (both competitors invalid) "
              f"({len(default_wins)}) ==")
        for name, label in sorted(default_wins):
            print(f"  {name} ({label})")
        print()
    if blocked:
        print(f"== BLOCKED: KLS invalid, competitor valid ({len(blocked)}) ==")
        for name, (label, reason, sides) in sorted(blocked.items()):
            print(f"  {name} ({label}): {reason} [loses to {sides}]")
        print()
    if dispositions:
        print(f"== DISPOSITIONS: all solvers fail ({len(dispositions)}) ==")
        for name, label in sorted(dispositions):
            print(f"  {name} ({label})")


if __name__ == "__main__":
    main()
