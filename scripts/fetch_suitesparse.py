#!/usr/bin/env python3
"""Fetch MatrixMarket matrices from the SuiteSparse Matrix Collection."""

from __future__ import annotations

import argparse
import io
import pathlib
import tarfile
import urllib.request


SSTATS_URL = "https://sparse.tamu.edu/files/ssstats.csv"
MM_URL = "https://sparse.tamu.edu/MM/{group}/{name}.tar.gz"


def read_manifest(path: pathlib.Path) -> list[str]:
    names: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        names.append(line)
    return names


def load_index() -> dict[str, tuple[str, str]]:
    with urllib.request.urlopen(SSTATS_URL, timeout=60) as response:
        text = response.read().decode("utf-8", errors="replace")
    result: dict[str, tuple[str, str]] = {}
    for line in text.splitlines()[2:]:
        parts = line.split(",", 2)
        if len(parts) < 2:
            continue
        group, name = parts[0].strip(), parts[1].strip()
        if not group or not name:
            continue
        result[name.lower()] = (group, name)
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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    wanted = read_manifest(args.manifest)
    index = load_index()
    missing: list[str] = []
    for requested in wanted:
        resolved = index.get(requested.lower())
        if resolved is None:
            missing.append(requested)
            continue
        group, canonical_name = resolved
        target = fetch_matrix(group, canonical_name, args.out)
        print(target)

    if missing:
        print("missing from SuiteSparse index:")
        for name in missing:
            print(f"  {name}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
