#!/usr/bin/env python3
"""Fetch MatrixMarket matrices from the SuiteSparse Matrix Collection."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import io
import pathlib
import tarfile
import urllib.request


SSTATS_URL = "https://sparse.tamu.edu/files/ssstats.csv"
MM_URL = "https://sparse.tamu.edu/MM/{group}/{name}.tar.gz"


@dataclass(frozen=True)
class MatrixInfo:
    group: str
    name: str
    nrows: int
    ncols: int
    nnz: int


def read_manifest(path: pathlib.Path) -> list[str]:
    names: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        names.append(line)
    return names


def load_index() -> dict[str, MatrixInfo]:
    with urllib.request.urlopen(SSTATS_URL, timeout=60) as response:
        text = response.read().decode("utf-8", errors="replace")
    result: dict[str, MatrixInfo] = {}
    reader = csv.reader(io.StringIO("\n".join(text.splitlines()[2:])))
    for row in reader:
        if len(row) < 5:
            continue
        group, name = row[0].strip(), row[1].strip()
        if not group or not name:
            continue
        try:
            nrows = int(float(row[2]))
            ncols = int(float(row[3]))
            nnz = int(float(row[4]))
        except ValueError:
            continue
        result[name.lower()] = MatrixInfo(group, name, nrows, ncols, nnz)
    return result


def fetch_matrix(group: str, name: str, out_dir: pathlib.Path) -> pathlib.Path:
    url = MM_URL.format(group=group, name=name)
    out_dir.mkdir(parents=True, exist_ok=True)
    target_dir = out_dir / group
    target_dir.mkdir(parents=True, exist_ok=True)
    target = target_dir / f"{name}.mtx"
    if target.exists():
        return target

    print(f"fetch {group}/{name}")
    with urllib.request.urlopen(url, timeout=300) as response:
        payload = response.read()
    with tarfile.open(fileobj=io.BytesIO(payload), mode="r:gz") as archive:
        member_name = f"{name}/{name}.mtx"
        member = archive.getmember(member_name)
        extracted = archive.extractfile(member)
        if extracted is None:
            raise RuntimeError(f"{member_name} missing in archive")
        target.write_bytes(extracted.read())
    return target


def skip_reason(info: MatrixInfo, args: argparse.Namespace) -> str | None:
    if args.max_rows is not None and info.nrows > args.max_rows:
        return f"rows {info.nrows} > {args.max_rows}"
    if args.max_cols is not None and info.ncols > args.max_cols:
        return f"cols {info.ncols} > {args.max_cols}"
    if args.max_nnz is not None and info.nnz > args.max_nnz:
        return f"nnz {info.nnz} > {args.max_nnz}"
    return None


def describe(info: MatrixInfo) -> str:
    return (
        f"{info.group}/{info.name} "
        f"rows={info.nrows} cols={info.ncols} nnz={info.nnz}"
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--max-rows", type=int)
    parser.add_argument("--max-cols", type=int)
    parser.add_argument("--max-nnz", type=int)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--write-resolved-manifest", type=pathlib.Path)
    args = parser.parse_args()

    wanted = read_manifest(args.manifest)
    index = load_index()
    missing: list[str] = []
    skipped: list[tuple[MatrixInfo, str]] = []
    selected: list[MatrixInfo] = []
    for requested in wanted:
        resolved = index.get(requested.lower())
        if resolved is None:
            missing.append(requested)
            continue
        reason = skip_reason(resolved, args)
        if reason is not None:
            skipped.append((resolved, reason))
            print(f"skip {describe(resolved)}: {reason}")
            continue
        if args.limit is not None and len(selected) >= args.limit:
            skipped.append((resolved, f"limit {args.limit} reached"))
            print(f"skip {describe(resolved)}: limit {args.limit} reached")
            continue
        selected.append(resolved)
        if args.dry_run:
            print(f"select {describe(resolved)}")
        else:
            target = fetch_matrix(resolved.group, resolved.name, args.out)
            print(target)

    if args.write_resolved_manifest is not None:
        args.write_resolved_manifest.parent.mkdir(parents=True, exist_ok=True)
        with args.write_resolved_manifest.open("w", encoding="utf-8") as out:
            for info in selected:
                out.write(f"{info.name}\n")

    if missing:
        print("missing from SuiteSparse index:")
        for name in missing:
            print(f"  {name}")
        return 1
    print(
        "summary: "
        f"requested={len(wanted)} selected={len(selected)} "
        f"skipped={len(skipped)} missing={len(missing)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
