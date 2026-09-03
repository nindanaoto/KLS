#!/usr/bin/env python3
"""Select a one-use validation suite outside every committed KLS manifest.

The selector operates only on a pinned SuiteSparse metadata snapshot.  It
resolves all matrices already named by KLS manifests.  Protocol v1 excludes
their complete SuiteSparse groups; v2 excludes the exact exposed matrices,
then chooses at most one unseen matrix per group.  Both balance the result
across size, density, and structural-symmetry strata.
"""

from __future__ import annotations

import argparse
import collections
import csv
from dataclasses import dataclass
import hashlib
import json
import os
import pathlib
import subprocess
import sys


PROTOCOL = "kls-external-validation-v1"
CIRCUIT_TERMS = (
    "circuit",
    "semiconductor",
    "power network",
    "electromagnetics",
)


@dataclass(frozen=True)
class Matrix:
    group: str
    name: str
    rows: int
    cols: int
    nnz: int
    real: bool
    symmetry: float
    kind: str

    @property
    def canonical(self) -> str:
        return f"{self.group}/{self.name}"


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def stable_score(seed: str, purpose: str, value: str) -> bytes:
    return hashlib.sha256(
        f"{PROTOCOL}\0{seed}\0{purpose}\0{value.casefold()}".encode()
    ).digest()


def load_index(path: pathlib.Path) -> tuple[str, str, list[Matrix]]:
    payload = path.read_bytes()
    lines = payload.decode("utf-8", errors="replace").splitlines()
    if len(lines) < 3:
        raise ValueError(f"{path}: invalid SuiteSparse index")
    matrices: list[Matrix] = []
    for row in csv.reader(lines[2:]):
        if len(row) < 12:
            continue
        try:
            matrix = Matrix(
                group=row[0].strip(),
                name=row[1].strip(),
                rows=int(float(row[2])),
                cols=int(float(row[3])),
                nnz=int(float(row[4])),
                real=int(float(row[5])) == 1,
                symmetry=float(row[9]),
                kind=row[11].strip(),
            )
        except ValueError:
            continue
        if matrix.group and matrix.name:
            matrices.append(matrix)
    if not matrices:
        raise ValueError(f"{path}: no usable matrix records")
    return lines[1].strip(), sha256_bytes(payload), matrices


def manifest_names(path: pathlib.Path) -> list[str]:
    names: list[str] = []
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if raw.lstrip().startswith("{"):
            try:
                row = json.loads(raw)
            except json.JSONDecodeError:
                continue
            matrix = row.get("matrix") if isinstance(row, dict) else None
            if matrix:
                names.append(pathlib.Path(str(matrix)).stem)
            continue
        token = raw.split("#", 1)[0].strip()
        if token:
            names.append(pathlib.Path(token.split()[0]).stem)
    return names


def resolve_excluded_groups(
    matrices: list[Matrix], manifests: list[pathlib.Path], strict: bool = True
) -> tuple[set[str], set[str], int]:
    by_name: dict[str, list[Matrix]] = collections.defaultdict(list)
    for matrix in matrices:
        by_name[matrix.name.casefold()].append(matrix)
    groups: set[str] = set()
    canonical: set[str] = set()
    count = 0
    for manifest in manifests:
        for name in manifest_names(manifest):
            matches = by_name.get(name.casefold(), [])
            if strict and len(matches) != 1:
                raise ValueError(
                    f"{manifest}: {name!r} resolves to {len(matches)} index rows"
                )
            for match in matches:
                groups.add(match.group.casefold())
                canonical.add(match.canonical.casefold())
            count += int(bool(matches))
    return groups, canonical, count


def stratum(matrix: Matrix) -> tuple[str, str, str]:
    size = "medium" if matrix.rows < 100_000 else "large"
    degree = matrix.nnz / matrix.rows
    density = "sparse" if degree < 5.0 else "moderate" if degree < 20.0 else "dense"
    symmetry = "symmetric" if matrix.symmetry >= 0.9 else "unsymmetric"
    return size, density, symmetry


def is_circuit_like(matrix: Matrix) -> bool:
    kind = matrix.kind.casefold()
    return any(term in kind for term in CIRCUIT_TERMS)


def choose_representatives(
    candidates: list[Matrix], seed: str
) -> list[Matrix]:
    grouped: dict[str, list[Matrix]] = collections.defaultdict(list)
    for matrix in candidates:
        grouped[matrix.group.casefold()].append(matrix)
    return [
        min(
            group,
            key=lambda matrix: stable_score(seed, "representative", matrix.canonical),
        )
        for group in grouped.values()
    ]


def balanced_sample(pool: list[Matrix], count: int, seed: str) -> list[Matrix]:
    buckets: dict[tuple[str, str, str], list[Matrix]] = collections.defaultdict(list)
    for matrix in pool:
        buckets[stratum(matrix)].append(matrix)
    for key, bucket in buckets.items():
        bucket.sort(
            key=lambda matrix: stable_score(
                seed, "sample/" + "/".join(key), matrix.canonical
            ),
            reverse=True,
        )
    selected: list[Matrix] = []
    while len(selected) < count:
        progressed = False
        for key in sorted(buckets):
            if buckets[key] and len(selected) < count:
                selected.append(buckets[key].pop())
                progressed = True
        if not progressed:
            raise ValueError(
                f"only {len(selected)} unused groups satisfy the requested filters"
            )
    return selected


def git_output(root: pathlib.Path, *args: str) -> str:
    return subprocess.check_output(
        ["git", "-C", str(root), *args], text=True
    ).strip()


def source_freeze(root: pathlib.Path, scientific_stack: bool = False) -> tuple[str, str]:
    source_paths = ("CMakeLists.txt", "include", "src", "third_party")
    if scientific_stack:
        source_paths += ("bench", "scripts")
    dirty = subprocess.run(
        ["git", "-C", str(root), "diff", "--quiet", "HEAD", "--", *source_paths],
        check=False,
    )
    if dirty.returncode != 0:
        raise ValueError("solver source is dirty; freeze or revert it before reveal")
    if scientific_stack:
        untracked = git_output(
            root, "ls-files", "--others", "--exclude-standard", "--",
            "bench", "scripts",
        ).splitlines()
        unsafe_untracked = [
            path for path in untracked
            if ((path.startswith("bench/") and
                 pathlib.Path(path).suffix in {".c", ".cc", ".cpp", ".h", ".hpp"})
                or (path.startswith("scripts/") and
                    pathlib.Path(path).suffix in {".py", ".sh"})
                or path in {
                    "bench/paper_campaign_configs.json",
                    "bench/paper_scaling_subset_manifest.txt",
                })
        ]
        if unsafe_untracked:
            raise ValueError(
                "scientific stack has untracked executable/configuration files; "
                "commit it before reveal: " + ", ".join(unsafe_untracked)
            )
    # Benchmark scripts/docs may be committed after reveal.  Record the last
    # commit that actually changed the solver tree, not an unrelated future
    # repository HEAD.
    commit = git_output(root, "log", "-1", "--format=%H", "--", *source_paths)
    tree_listing = git_output(root, "ls-tree", "-r", "HEAD", "--", *source_paths)
    return commit, sha256_bytes((tree_listing + "\n").encode())


def render(
    selected: list[Matrix],
    *,
    timestamp: str,
    index_digest: str,
    exclusion_digest: str,
    excluded_groups: int,
    excluded_entries: int,
    source_commit: str,
    source_digest: str,
    args: argparse.Namespace,
) -> str:
    strata = collections.Counter(stratum(matrix) for matrix in selected)
    lines = [
        "# One-use KLS external validation manifest. Do not tune after reveal.",
        f"# Protocol: {args.protocol}",
        f"# Seed: {args.seed}",
        f"# SuiteSparse index timestamp: {timestamp}",
        f"# SuiteSparse index SHA256: {index_digest}",
        f"# Exclusion-set SHA256: {exclusion_digest}",
        f"# Excluded manifest entries: {excluded_entries}",
        (f"# Referenced SuiteSparse groups: {excluded_groups}"
         if args.protocol.endswith("v2")
         else f"# Excluded SuiteSparse groups: {excluded_groups}"),
        f"# Frozen KLS commit: {source_commit}",
        f"# Frozen solver-tree SHA256: {source_digest}",
        (
            "# Filters: real square, unique basename, one matrix per unused group; "
            f"rows={args.min_rows}..{args.max_rows}; nnz={args.min_nnz}..{args.max_nnz}."
        ),
        ("# Selection: one previously unseen matrix per group, circuit-like first, "
         "then balanced metadata strata."
         if args.protocol.endswith("v2")
         else "# Selection: all available circuit-like unused groups, then balanced metadata strata."),
        "# Strata: "
        + ", ".join(
            f"{'/'.join(key)}={value}" for key, value in sorted(strata.items())
        ),
        "",
    ]
    for matrix in selected:
        kind = " ".join(matrix.kind.replace("#", " ").split())
        lines.append(
            f"{matrix.name}  # group={matrix.group} rows={matrix.rows} "
            f"nnz={matrix.nnz} stratum={'/'.join(stratum(matrix))} "
            f"circuit_like={int(is_circuit_like(matrix))} kind={kind}"
        )
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--index", type=pathlib.Path, required=True)
    parser.add_argument("--root", type=pathlib.Path, default=pathlib.Path("."))
    parser.add_argument(
        "--protocol", choices=("kls-external-validation-v1", "kls-external-validation-v2"),
        default="kls-external-validation-v1",
    )
    parser.add_argument(
        "--exclude-glob", action="append", default=[],
        help="glob, relative to --root, containing every previously exposed suite",
    )
    parser.add_argument(
        "--output",
        type=pathlib.Path,
        default=None,
    )
    parser.add_argument("--count", type=int, default=24)
    parser.add_argument("--min-rows", type=int, default=10_000)
    parser.add_argument("--max-rows", type=int, default=1_000_000)
    parser.add_argument("--min-nnz", type=int, default=30_000)
    parser.add_argument("--max-nnz", type=int, default=12_000_000)
    parser.add_argument("--seed", default="kls-external-validation-2026-08-v1")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    global PROTOCOL
    PROTOCOL = args.protocol
    if args.count <= 0:
        parser.error("--count must be positive")
    if args.check and args.dry_run:
        parser.error("--check and --dry-run are mutually exclusive")
    root = args.root.resolve()
    index = args.index.resolve()
    if args.output is None:
        suffix = "v2" if args.protocol.endswith("v2") else "v1"
        args.output = pathlib.Path(f"bench/suitesparse_external_validation_{suffix}_manifest.txt")
    output = args.output if args.output.is_absolute() else root / args.output
    try:
        timestamp, index_digest, matrices = load_index(index)
        exclusion_globs = args.exclude_glob or ["bench/*manifest*.txt"]
        manifests = sorted({
            path.resolve()
            for pattern in exclusion_globs for path in root.glob(pattern)
            if path.is_file() and path.resolve() != output.resolve()
        })
        if not manifests:
            raise ValueError(f"no exclusion manifests match {exclusion_globs!r}")
        excluded_groups, excluded_matrices, excluded_entries = resolve_excluded_groups(
            matrices, manifests, strict=not args.protocol.endswith("v2")
        )
        exclusion_payload = b"".join(
            pathlib.PurePath(os.path.relpath(path, root)).as_posix().encode()
            + b"\0" + path.read_bytes() + b"\0"
            for path in manifests
        )
        # Freeze before deriving the selection, not merely before writing it.
        # Thus a dirty invocation never computes or emits the one-use corpus.
        source_commit, source_digest = source_freeze(
            root, scientific_stack=args.protocol.endswith("v2")
        )
        name_counts = collections.Counter(matrix.name.casefold() for matrix in matrices)
        candidates = [
            matrix for matrix in matrices
            if matrix.real
            and matrix.rows == matrix.cols
            and args.min_rows <= matrix.rows <= args.max_rows
            and args.min_nnz <= matrix.nnz <= args.max_nnz
            and (matrix.canonical.casefold() not in excluded_matrices
                 if args.protocol.endswith("v2")
                 else matrix.group.casefold() not in excluded_groups)
            and name_counts[matrix.name.casefold()] == 1
        ]
        representatives = choose_representatives(candidates, args.seed)
        circuit = sorted(
            (matrix for matrix in representatives if is_circuit_like(matrix)),
            key=lambda matrix: stable_score(args.seed, "circuit", matrix.canonical),
        )
        if len(circuit) > args.count:
            circuit = circuit[: args.count]
        remaining = [matrix for matrix in representatives if matrix not in circuit]
        selected = circuit + balanced_sample(remaining, args.count - len(circuit), args.seed)
        selected.sort(key=lambda matrix: matrix.canonical.casefold())
        text = render(
            selected,
            timestamp=timestamp,
            index_digest=index_digest,
            exclusion_digest=sha256_bytes(exclusion_payload),
            excluded_groups=len(excluded_groups),
            excluded_entries=excluded_entries,
            source_commit=source_commit,
            source_digest=source_digest,
            args=args,
        )
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        print(exc, file=sys.stderr)
        return 1

    if args.dry_run:
        print(text, end="")
    elif args.check:
        try:
            actual = output.read_text(encoding="utf-8")
        except OSError:
            actual = ""
        if actual != text:
            print(f"stale external validation manifest: {output}", file=sys.stderr)
            return 1
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(text, encoding="utf-8")
    print(
        f"external validation: selected={len(selected)} "
        f"unused_groups={len(representatives)} circuit_like={len(circuit)} "
        f"excluded_groups={len(excluded_groups)}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
