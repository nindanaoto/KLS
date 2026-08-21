#!/usr/bin/env python3
"""List LLC domains or select one logical CPU per physical core."""

from __future__ import annotations

import argparse
import os
import pathlib
import re
import sys


SYS_CPU = pathlib.Path("/sys/devices/system/cpu")


def read(path: pathlib.Path) -> str:
    return path.read_text(encoding="ascii").strip()


def cache_bytes(value: str) -> int:
    match = re.fullmatch(r"(\d+)([KMG]?)", value, re.IGNORECASE)
    if match is None:
        raise ValueError(f"unrecognized cache size: {value}")
    scale = {"": 1, "K": 1024, "M": 1024**2, "G": 1024**3}
    return int(match.group(1)) * scale[match.group(2).upper()]


def cpu_llc(cpu: int) -> tuple[str, int]:
    cache_root = SYS_CPU / f"cpu{cpu}" / "cache"
    candidates: list[tuple[int, pathlib.Path]] = []
    for index in cache_root.glob("index*"):
        try:
            cache_type = read(index / "type")
            level = int(read(index / "level"))
        except (OSError, ValueError):
            continue
        if cache_type in ("Unified", "Data"):
            candidates.append((level, index))
    if not candidates:
        raise RuntimeError(f"no data cache topology for CPU {cpu}")
    _, index = max(candidates, key=lambda item: item[0])
    package = read(SYS_CPU / f"cpu{cpu}" / "topology" / "physical_package_id")
    shared = read(index / "shared_cpu_list")
    return f"package={package},cpus={shared}", cache_bytes(read(index / "size"))


def domains() -> list[dict[str, object]]:
    allowed = sorted(os.sched_getaffinity(0))
    grouped: dict[str, dict[str, object]] = {}
    for cpu in allowed:
        try:
            key, size = cpu_llc(cpu)
            core = int(read(SYS_CPU / f"cpu{cpu}" / "topology" / "core_id"))
            package = int(read(
                SYS_CPU / f"cpu{cpu}" / "topology" / "physical_package_id"
            ))
        except (OSError, RuntimeError, ValueError) as exc:
            print(f"warning: CPU {cpu}: {exc}", file=sys.stderr)
            continue
        domain = grouped.setdefault(
            key, {"key": key, "bytes": size, "cores": {}, "package": package}
        )
        cores = domain["cores"]
        assert isinstance(cores, dict)
        cores.setdefault(core, cpu)
    result = list(grouped.values())
    result.sort(key=lambda item: (int(item["bytes"]), str(item["key"])))
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--select", choices=("smallest", "largest"),
        help="print a comma-separated physical-core CPU set",
    )
    parser.add_argument("--threads", type=int, default=0)
    args = parser.parse_args()
    found = domains()
    if not found:
        parser.error("no LLC domains are visible in the current CPU affinity")
    if args.select is None:
        print("domain\tLLC_MiB\tphysical_cores\tselected_logical_CPUs")
        for index, domain in enumerate(found):
            cores = domain["cores"]
            assert isinstance(cores, dict)
            cpus = [cores[key] for key in sorted(cores)]
            print(f"{index}\t{int(domain['bytes']) / 1024**2:.1f}\t"
                  f"{len(cpus)}\t{','.join(map(str, cpus))}")
        return 0
    domain = found[0] if args.select == "smallest" else found[-1]
    cores = domain["cores"]
    assert isinstance(cores, dict)
    cpus = [cores[key] for key in sorted(cores)]
    count = args.threads or len(cpus)
    if count < 1 or count > len(cpus):
        parser.error(
            f"selected domain has {len(cpus)} physical cores, cannot use {count}"
        )
    print(",".join(map(str, cpus[:count])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
