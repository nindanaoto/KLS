#!/usr/bin/env python3
"""Strict paired reduction for run_paper_campaign.py output."""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import io
import json
import math
import pathlib
import random
import re
import statistics
import sys
from typing import Any


SOLVERS = ("kls", "ck", "st", "klu")
SCHEMA = "kls-paper-observation-v2"


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def manifest_names(path: pathlib.Path) -> list[str]:
    result = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        token = raw.split("#", 1)[0].strip()
        if token:
            result.append(pathlib.Path(token.split()[0]).stem)
    return result


def load_jsonl(path: pathlib.Path) -> list[dict[str, Any]]:
    rows = []
    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not raw.strip():
            continue
        try:
            row = json.loads(raw)
        except json.JSONDecodeError as exc:
            raise ValueError(f"{path}:{line_number}: {exc}") from exc
        if not isinstance(row, dict):
            raise ValueError(f"{path}:{line_number}: expected object")
        rows.append(row)
    return rows


def time_of(row: dict[str, Any]) -> float:
    value = float(row["result"]["measured_lifecycle_seconds"])
    if not (value > 0.0 and math.isfinite(value)):
        raise ValueError("invalid lifecycle time")
    return value


def independently_valid(row: dict[str, Any], limit: float,
                        systems: int) -> tuple[bool, str]:
    if row.get("valid") is not True:
        return False, f"runner:{row.get('validity_reason', 'invalid')}"
    if row.get("timed_out") is True or row.get("returncode") != 0:
        return False, "runner-status"
    result = row.get("result")
    if not isinstance(result, dict):
        return False, "missing-result"
    if result.get("status") not in (None, 0):
        return False, f"solver-status:{result.get('status')}"
    if result.get("verify_each_refactor") is not True:
        return False, "refactors-not-verified"
    if (result.get("lifecycle_mode") != "direct" or
            result.get("lifecycle_systems") != systems):
        return False, "wrong-lifecycle"
    for key in ("relative_residual_l2", "refactor_max_relative_residual"):
        try:
            value = float(result[key])
        except (KeyError, TypeError, ValueError):
            return False, f"missing:{key}"
        if not (math.isfinite(value) and 0.0 <= value <= limit):
            return False, f"residual:{key}={value:.6g}"
    try:
        elapsed = float(result["measured_lifecycle_seconds"])
    except (KeyError, TypeError, ValueError):
        return False, "missing:measured_lifecycle_seconds"
    if not (math.isfinite(elapsed) and elapsed > 0.0):
        return False, "invalid:measured_lifecycle_seconds"
    return True, "valid"


def geomean(values: list[float]) -> float | None:
    if not values:
        return None
    return math.exp(statistics.fmean(math.log(value) for value in values))


def bootstrap_ci(values: list[float], samples: int, seed: int) -> list[float] | None:
    if not values:
        return None
    rng = random.Random(seed)
    estimates = []
    for _ in range(samples):
        draw = [values[rng.randrange(len(values))] for _ in values]
        estimates.append(geomean(draw))
    estimates.sort()
    lo = estimates[int(0.025 * (samples - 1))]
    hi = estimates[int(0.975 * (samples - 1))]
    return [lo, hi]


def paired_summary(ratios: dict[str, float], bootstrap: int, seed: int) -> dict[str, Any]:
    values = list(ratios.values())
    return {
        "valid_pairs": len(values),
        "kls_wins": sum(value > 1.02 for value in values),
        "ties_within_2pct": sum(0.98 <= value <= 1.02 for value in values),
        "kls_losses": sum(value < 0.98 for value in values),
        "challenger_over_kls_geomean": geomean(values),
        "geomean_bootstrap_95ci": bootstrap_ci(values, bootstrap, seed),
        "challenger_over_kls_median": statistics.median(values) if values else None,
        "matrix_ratios": dict(sorted(ratios.items())),
    }


def configuration_summary(ratios: dict[str, float], bootstrap: int,
                          seed: int) -> dict[str, Any]:
    values = list(ratios.values())
    inverse = [1.0 / value for value in values]
    return {
        "valid_pairs": len(values),
        "arm_faster": sum(value < 0.98 for value in values),
        "ties_within_2pct": sum(0.98 <= value <= 1.02 for value in values),
        "baseline_faster": sum(value > 1.02 for value in values),
        "arm_over_baseline_geomean": geomean(values),
        "arm_over_baseline_bootstrap_95ci": bootstrap_ci(values, bootstrap, seed),
        "baseline_over_arm_geomean": geomean(inverse),
        "arm_over_baseline_median": statistics.median(values) if values else None,
        "matrix_arm_over_baseline": dict(sorted(ratios.items())),
    }


def markdown(summary: dict[str, Any]) -> str:
    lines = [
        "# KLS paper campaign", "",
        "All matrix results require every configured pass to be present and valid. ",
        "Ratios are formed within pass before matrix-level and corpus-level aggregation.", "",
    ]
    for config_id, config in summary["configurations"].items():
        lines += [f"## {config_id}", "", "| Solver | Complete valid matrices | Invalid matrices | Paired GM | 95% CI | W/T/L |",
                  "| --- | ---: | ---: | ---: | ---: | ---: |"]
        for solver in config["solvers"]:
            coverage = config["coverage"][solver]
            if solver == "kls":
                lines.append(f"| KLS | {coverage['valid_matrices']} | {coverage['invalid_matrices']} | — | — | — |")
                continue
            paired = config["pairwise"].get(solver)
            if paired is None:
                lines.append(f"| {solver} | {coverage['valid_matrices']} | {coverage['invalid_matrices']} | — | — | — |")
                continue
            gm = paired["challenger_over_kls_geomean"]
            ci = paired["geomean_bootstrap_95ci"]
            gm_text = "—" if gm is None else f"{gm:.3f}"
            ci_text = "—" if ci is None else f"[{ci[0]:.3f}, {ci[1]:.3f}]"
            wtl = f"{paired['kls_wins']}/{paired['ties_within_2pct']}/{paired['kls_losses']}"
            lines.append(f"| {solver} | {coverage['valid_matrices']} | {coverage['invalid_matrices']} | {gm_text} | {ci_text} | {wtl} |")
        lines.append("")
    if summary["incomplete_observations"]:
        lines += ["## Incomplete observations", "", *[f"- `{item}`" for item in summary["incomplete_observations"]], ""]
    if summary.get("configuration_comparisons"):
        lines += ["## Paired configuration comparisons", "",
                  "Ratios below one favor the arm; speedup is baseline/arm.", "",
                  "| Arm | Baseline | Solver | Pairs | Arm/base GM | 95% CI | Speedup | A/T/B |",
                  "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |"]
        for arm, comparison in summary["configuration_comparisons"].items():
            for solver, paired in comparison["per_solver"].items():
                ci = paired["arm_over_baseline_bootstrap_95ci"]
                ci_text = "—" if ci is None else f"[{ci[0]:.3f}, {ci[1]:.3f}]"
                ratio = paired["arm_over_baseline_geomean"]
                speedup = paired["baseline_over_arm_geomean"]
                ratio_text = "—" if ratio is None else f"{ratio:.3f}"
                speedup_text = "—" if speedup is None else f"{speedup:.3f}"
                atb = (f"{paired['arm_faster']}/{paired['ties_within_2pct']}/"
                       f"{paired['baseline_faster']}")
                lines.append(f"| {arm} | {comparison['baseline']} | {solver} | "
                             f"{paired['valid_pairs']} | {ratio_text} | {ci_text} | "
                             f"{speedup_text} | {atb} |")
        lines.append("")
    return "\n".join(lines)


def latex(summary: dict[str, Any]) -> str:
    """Return reproducible table rows for inclusion in a paper table."""
    lines = [
        "% Generated by scripts/reduce_paper_campaign.py; do not edit.",
        "% config & solver & valid & invalid & paired GM & 95% CI & W/T/L",
    ]
    for config_id, config in summary["configurations"].items():
        safe_config = re.sub(r"([%&#_$])", r"\\\1", config_id)
        for solver in config["solvers"]:
            coverage = config["coverage"][solver]
            paired = config["pairwise"].get(solver)
            if solver == "kls" or paired is None:
                gm_text = ci_text = wtl = "--"
            else:
                gm = paired["challenger_over_kls_geomean"]
                ci = paired["geomean_bootstrap_95ci"]
                gm_text = "--" if gm is None else f"{gm:.3f}"
                ci_text = "--" if ci is None else f"[{ci[0]:.3f}, {ci[1]:.3f}]"
                wtl = (f"{paired['kls_wins']}/{paired['ties_within_2pct']}/"
                       f"{paired['kls_losses']}")
            lines.append(
                f"\\texttt{{{safe_config}}} & {solver.upper()} & "
                f"{coverage['valid_matrices']} & {coverage['invalid_matrices']} & "
                f"{gm_text} & {ci_text} & {wtl} \\\\"
            )
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=pathlib.Path, required=True)
    parser.add_argument("--json", type=pathlib.Path)
    parser.add_argument("--markdown", type=pathlib.Path)
    parser.add_argument("--csv", type=pathlib.Path)
    parser.add_argument("--latex", type=pathlib.Path)
    parser.add_argument("--bootstrap", type=int, default=10000)
    parser.add_argument("--only", action="append", default=[])
    parser.add_argument("--allow-incomplete", action="store_true")
    args = parser.parse_args()
    if args.bootstrap < 100:
        parser.error("--bootstrap must be at least 100")
    results = args.results.resolve()
    try:
        campaign = json.loads((results / "campaign.json").read_text(encoding="utf-8"))
        if campaign.get("schema") != "kls-paper-campaign-v2":
            raise ValueError("unsupported campaign schema")
        config_path = pathlib.Path(campaign["configs"])
        manifest_path = pathlib.Path(campaign["manifest"])
        if config_path.exists() and sha256(config_path) != campaign["configs_sha256"]:
            raise ValueError("configuration file changed after campaign creation")
        if manifest_path.exists() and sha256(manifest_path) != campaign["manifest_sha256"]:
            raise ValueError("manifest changed after campaign creation")
        all_configs = campaign.get("configuration_spec")
        if all_configs is None:
            all_configs = json.loads(config_path.read_text(encoding="utf-8"))
        all_names = campaign.get("manifest_entries")
        if all_names is None:
            all_names = manifest_names(manifest_path)
        configs = list(all_configs)
        if args.only:
            selected_ids = {
                str(config["id"]) for config in configs
                if any(token in str(config["id"]) for token in args.only)
            }
            if not selected_ids:
                raise ValueError("--only did not select any configuration")
            by_id = {str(config["id"]): config for config in configs}
            pending = list(selected_ids)
            while pending:
                baseline = by_id[pending.pop()].get("baseline")
                if baseline is not None and str(baseline) not in selected_ids:
                    selected_ids.add(str(baseline))
                    pending.append(str(baseline))
            configs = [config for config in configs
                       if str(config["id"]) in selected_ids]
        passes = int(campaign["passes"])
        limit = float(campaign["residual_limit"])
        observations = load_jsonl(results / "observations.jsonl")
        all_by_id = {str(config["id"]): config for config in all_configs}
        indexed: dict[tuple[str, str, int, str], dict[str, Any]] = {}
        for row in observations:
            if row.get("schema") != SCHEMA:
                raise ValueError(f"unsupported observation schema: {row.get('schema')}")
            config_id = str(row["config_id"])
            config = all_by_id.get(config_id)
            if config is None:
                raise ValueError(f"unknown observation configuration: {config_id}")
            matrix_id = str(row["matrix_id"])
            if matrix_id not in config.get("matrices", all_names):
                raise ValueError(f"{config_id}: unexpected matrix {matrix_id}")
            pass_id = int(row["pass_id"])
            if not 0 <= pass_id < passes:
                raise ValueError(f"{config_id}/{matrix_id}: invalid pass {pass_id}")
            solver = str(row["solver"])
            if solver not in config.get("solvers", SOLVERS):
                raise ValueError(f"{config_id}/{matrix_id}: unexpected solver {solver}")
            key = (str(row["config_id"]), str(row["matrix_id"]),
                   int(row["pass_id"]), str(row["solver"]))
            if key in indexed:
                raise ValueError(f"duplicate observation: {key}")
            valid, reason = independently_valid(
                row, limit, int(config.get("lifecycle_systems", 100))
            )
            row["_reducer_valid"] = valid
            row["_reducer_reason"] = reason
            indexed[key] = row
        summary: dict[str, Any] = {
            "campaign": campaign, "configurations": {},
            "incomplete_observations": [],
        }
        csv_rows: list[dict[str, Any]] = []
        valid_samples: dict[tuple[str, str, str], dict[int, dict[str, Any]]] = {}
        for config_index, config in enumerate(configs):
            config_id = str(config["id"])
            names = list(config.get("matrices", all_names))
            solvers = list(config.get("solvers", SOLVERS))
            coverage: dict[str, Any] = {}
            for solver in solvers:
                valid_count = 0
                failures: dict[str, list[str]] = {}
                for matrix in names:
                    samples: dict[int, dict[str, Any]] = {}
                    reasons = []
                    for pass_id in range(passes):
                        row = indexed.get((config_id, matrix, pass_id, solver))
                        if row is None:
                            item = f"{config_id}/{matrix}/p{pass_id}/{solver}"
                            summary["incomplete_observations"].append(item)
                            reasons.append(f"p{pass_id}:missing")
                        elif row.get("_reducer_valid") is not True:
                            reasons.append(f"p{pass_id}:{row.get('_reducer_reason')}")
                        else:
                            samples[pass_id] = row
                    if reasons:
                        failures[matrix] = reasons
                    else:
                        valid_count += 1
                        valid_samples[(config_id, matrix, solver)] = samples
                coverage[solver] = {
                    "valid_matrices": valid_count,
                    "invalid_matrices": len(failures), "failures": failures,
                }
            comparisons: dict[str, Any] = {}
            if "kls" in solvers:
                for challenger in solvers:
                    if challenger == "kls":
                        continue
                    ratios: dict[str, float] = {}
                    for matrix in names:
                        left = valid_samples.get((config_id, matrix, "kls"))
                        right = valid_samples.get((config_id, matrix, challenger))
                        if left is None or right is None:
                            continue
                        pass_ratios = [time_of(right[p]) / time_of(left[p])
                                       for p in range(passes)]
                        ratios[matrix] = math.exp(statistics.median(
                            [math.log(value) for value in pass_ratios]
                        ))
                        csv_rows.append({"config_id": config_id, "matrix": matrix,
                                         "challenger": challenger,
                                         "challenger_over_kls": ratios[matrix]})
                    comparisons[challenger] = paired_summary(
                        ratios, args.bootstrap, 0x4B4C53 + config_index * 17 + len(challenger)
                    )
            summary["configurations"][config_id] = {
                "solvers": solvers, "matrix_count": len(names),
                "coverage": coverage, "pairwise": comparisons,
            }

        # KLS-only ablations name their AUTO configuration in `baseline`.
        ablations: dict[str, Any] = {}
        for config_index, config in enumerate(configs):
            baseline = config.get("baseline")
            if not baseline:
                continue
            config_id = str(config["id"])
            names = list(config.get("matrices", all_names))
            ratios = {}
            for matrix in names:
                base = valid_samples.get((str(baseline), matrix, "kls"))
                arm = valid_samples.get((config_id, matrix, "kls"))
                if base is None or arm is None:
                    continue
                pass_ratios = [time_of(arm[p]) / time_of(base[p]) for p in range(passes)]
                ratios[matrix] = math.exp(statistics.median(
                    [math.log(value) for value in pass_ratios]
                ))
            ablations[config_id] = {
                "baseline": baseline,
                "arm_over_auto": paired_summary(
                    ratios, args.bootstrap, 0x41424C + config_index
                ),
            }
        summary["ablations"] = ablations

        # Any arm naming a baseline gets a strict within-pass comparison.  This
        # covers scaling, workload/horizon sensitivity, and ablations without
        # conflating their ratios with cross-solver headline ratios.
        config_ids = {str(config["id"]) for config in configs}
        config_comparisons: dict[str, Any] = {}
        for config_index, config in enumerate(configs):
            baseline = config.get("baseline")
            if not baseline:
                continue
            if str(baseline) not in config_ids:
                raise ValueError(f"{config['id']}: unknown baseline {baseline}")
            config_id = str(config["id"])
            names = list(config.get("matrices", all_names))
            baseline_config = next(item for item in configs
                                   if str(item["id"]) == str(baseline))
            solvers = sorted(set(config.get("solvers", SOLVERS)).intersection(
                baseline_config.get("solvers", SOLVERS)))
            per_solver: dict[str, Any] = {}
            for solver in solvers:
                ratios: dict[str, float] = {}
                for matrix in names:
                    base = valid_samples.get((str(baseline), matrix, solver))
                    arm = valid_samples.get((config_id, matrix, solver))
                    if base is None or arm is None:
                        continue
                    pass_ratios = [time_of(arm[p]) / time_of(base[p])
                                   for p in range(passes)]
                    ratios[matrix] = math.exp(statistics.median(
                        [math.log(value) for value in pass_ratios]
                    ))
                per_solver[solver] = configuration_summary(
                    ratios, args.bootstrap, 0x434647 + config_index * 19 + len(solver)
                )
            config_comparisons[config_id] = {
                "baseline": str(baseline), "per_solver": per_solver,
            }
        summary["configuration_comparisons"] = config_comparisons
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError) as exc:
        print(exc, file=sys.stderr)
        return 1

    if args.json:
        args.json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    report = markdown(summary)
    if args.markdown:
        args.markdown.write_text(report, encoding="utf-8")
    else:
        print(report, end="")
    if args.csv:
        buffer = io.StringIO()
        writer = csv.DictWriter(buffer, fieldnames=("config_id", "matrix", "challenger", "challenger_over_kls"))
        writer.writeheader()
        writer.writerows(csv_rows)
        args.csv.write_text(buffer.getvalue(), encoding="utf-8")
    if args.latex:
        args.latex.write_text(latex(summary), encoding="utf-8")
    if summary["incomplete_observations"] and not args.allow_incomplete:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
