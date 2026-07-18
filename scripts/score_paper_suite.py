#!/usr/bin/env python3
"""Score paired paper-manifest runs at a requested solve horizon.

H1   = analysis + initial factor + solve
T100 = analysis + initial factor + first refactor + 98 * steady refactor
       + 100 * solve

By default, each matrix is represented by the median valid pass. A pass is
valid only when it completed and its reported relative residual is finite and
at most 1e-8.
"""

import argparse
import json
import math
from pathlib import Path


HORIZON = 100


def value_protocol(path):
    protocols = set()
    with open(path, encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip():
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f"{path}:{line_number}: invalid JSON: {exc}") from exc
            if record.get("timeout_or_fail") or "analysis_seconds" not in record:
                continue
            mode = record.get("refactor_value_mode", "unchanged")
            try:
                amplitude = float(record.get("refactor_value_amplitude", 0.0))
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"{path}:{line_number}: invalid refactor value amplitude"
                ) from exc
            if mode not in ("unchanged", "rank-preserving", "entrywise") or not math.isfinite(
                amplitude
            ):
                raise ValueError(
                    f"{path}:{line_number}: invalid refactor value protocol"
                )
            protocols.add((mode, amplitude))
    if len(protocols) > 1:
        rendered = ", ".join(f"{mode}@{amplitude:g}" for mode, amplitude in sorted(protocols))
        raise ValueError(f"{path}: mixed refactor value protocols: {rendered}")
    return next(iter(protocols), None)


def matrix_name(record):
    return Path(record.get("matrix", "")).stem.lower()


def relative_residual(record):
    value = record.get("relative_residual_l2")
    if value is None:
        value = record.get("max_relative_residual")
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def score_record(record, residual_threshold, horizon=HORIZON):
    required = ("analysis_seconds", "solve_seconds_avg")
    if record.get("timeout_or_fail") or any(key not in record for key in required):
        return None, "timeout/fail"

    residual = relative_residual(record)
    if residual is None:
        return None, "missing/non-finite residual"
    if residual > residual_threshold:
        return None, f"residual {residual:.3e} > {residual_threshold:.3e}"

    try:
        analysis = float(record["analysis_seconds"])
        initial = float(
            record.get("initial_factor_seconds", record.get("factor_seconds", 0.0))
        )
        first = float(
            record.get("refactor_first_seconds", record.get("refactor_seconds_avg", 0.0))
        )
        steady = float(
            record.get(
                "refactor_steady_seconds_avg", record.get("refactor_seconds_avg", 0.0)
            )
        )
        solve = float(record["solve_seconds_avg"])
    except (TypeError, ValueError):
        return None, "non-numeric timing"

    timings = (analysis, initial, first, steady, solve)
    if any(not math.isfinite(value) or value < 0.0 for value in timings):
        return None, "invalid timing"

    if horizon == 1:
        cycle = analysis + initial + solve
    else:
        cycle = analysis + initial + first + (horizon - 2) * steady + horizon * solve
    return {
        "cycle": cycle,
        "ana": analysis,
        "init": initial,
        "rf": first,
        "rs": steady,
        "solve": solve,
        "resid": residual,
        "path": record.get("last_refactor_path", ""),
        "ordering": record.get("ordering", ""),
    }, None


def load(path, residual_threshold, selection, horizon=HORIZON):
    samples = {}
    failures = {}
    with open(path, encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            line = line.strip()
            if not line:
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f"{path}:{line_number}: invalid JSON: {exc}") from exc
            name = matrix_name(record)
            if not name:
                continue
            scored, reason = score_record(record, residual_threshold, horizon)
            if scored is None:
                failures.setdefault(name, []).append(reason)
            else:
                samples.setdefault(name, []).append(scored)

    selected = {}
    for name in samples.keys() | failures.keys():
        valid = sorted(samples.get(name, []), key=lambda item: item["cycle"])
        if not valid:
            selected[name] = None
        elif selection == "min":
            selected[name] = valid[0]
        else:
            selected[name] = valid[len(valid) // 2]
    return selected, failures


def geometric_mean(values):
    values = [value for value in values if math.isfinite(value) and value > 0.0]
    return (
        math.exp(sum(math.log(value) for value in values) / len(values))
        if values
        else float("nan")
    )


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kls", help="KLS JSONL run")
    parser.add_argument("ck", help="CKTSO JSONL run")
    parser.add_argument("st", nargs="?", help="optional SubtreeLU JSONL run")
    parser.add_argument(
        "--residual-threshold",
        type=float,
        default=1e-8,
        help="maximum accepted relative residual (default: 1e-8)",
    )
    parser.add_argument(
        "--selection",
        choices=("median", "min"),
        default="median",
        help="pass selection per matrix (default: median)",
    )
    parser.add_argument(
        "--horizon",
        type=int,
        default=HORIZON,
        help="number of solve iterations to score; 1 is a cold one-shot (default: 100)",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if args.horizon < 1:
        raise SystemExit("--horizon must be at least 1")
    protocol_paths = [("KLS", args.kls), ("CK", args.ck)]
    if args.st:
        protocol_paths.append(("ST", args.st))
    protocols = [(name, value_protocol(path)) for name, path in protocol_paths]
    present_protocols = {protocol for _, protocol in protocols if protocol is not None}
    if len(present_protocols) > 1:
        rendered = ", ".join(
            f"{name}={protocol[0]}@{protocol[1]:g}"
            for name, protocol in protocols
            if protocol is not None
        )
        raise SystemExit(f"refactor value protocol mismatch: {rendered}")
    protocol = next(iter(present_protocols), ("unknown", 0.0))
    kls, kls_failures = load(
        args.kls, args.residual_threshold, args.selection, args.horizon
    )
    comparisons = {
        "CK": load(args.ck, args.residual_threshold, args.selection, args.horizon)
    }
    if args.st:
        comparisons["ST"] = load(
            args.st, args.residual_threshold, args.selection, args.horizon
        )

    kls_valid = sum(record is not None for record in kls.values())
    print(
        f"KLS valid: {kls_valid}/{len(kls)} "
        f"(selection={args.selection}, residual<={args.residual_threshold:.1e}, "
        f"values={protocol[0]}@{protocol[1]:g})"
    )
    for comparison_name, (comparison, comparison_failures) in comparisons.items():
        comparison_valid = sum(record is not None for record in comparison.values())
        print(f"{comparison_name} valid: {comparison_valid}/{len(comparison)}")
        rows = []
        for name in sorted(kls.keys() | comparison.keys()):
            candidate = kls.get(name)
            reference = comparison.get(name)
            if candidate is None or reference is None:
                candidate_reason = ", ".join(kls_failures.get(name, ["missing"]))
                reference_reason = ", ".join(
                    comparison_failures.get(name, ["missing"])
                )
                print(
                    f"  {name}: SKIP "
                    f"(kls={'ok' if candidate else candidate_reason}; "
                    f"{comparison_name.lower()}="
                    f"{'ok' if reference else reference_reason})"
                )
                continue
            rows.append((candidate["cycle"] / reference["cycle"], name, candidate, reference))

        rows.sort(reverse=True)
        print(
            f"\n== KLS vs {comparison_name}: "
            f"H{args.horizon} gm={geometric_mean([row[0] for row in rows]):.3f} "
            f"wins {sum(row[0] < 1.0 for row in rows)}/{len(rows)} =="
        )
        print(
            f"{'matrix':16} {f'H{args.horizon}':>6} {'rs':>6} {'solve':>6} "
            f"{'ana+init':>8}  path/ord      resid"
        )
        for ratio, name, candidate, reference in rows:
            front_end = (candidate["ana"] + candidate["init"]) / max(
                reference["ana"] + reference["init"], 1e-12
            )
            print(
                f"{name:16} {ratio:6.2f} "
                f"{candidate['rs'] / max(reference['rs'], 1e-12):6.2f} "
                f"{candidate['solve'] / max(reference['solve'], 1e-12):6.2f} "
                f"{front_end:8.2f}  "
                f"{candidate['path']}/{candidate['ordering']:6} "
                f"{candidate['resid']:.1e}"
            )


if __name__ == "__main__":
    main()
