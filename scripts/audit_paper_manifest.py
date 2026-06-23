#!/usr/bin/env python3
"""Check that the paper-derived SuiteSparse manifest covers known public cases."""

from __future__ import annotations

import argparse
import pathlib
import re
import sys


def normalize(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


EXPECTED_BY_SOURCE: dict[str, set[str]] = {
    "KLU Algorithm 907": {
        "asic_100ks",
        "asic_320k",
        "asic_320ks",
        "asic_680k",
        "asic_680ks",
        "ckt11752_dc_1",
        "freescale1",
        "raj1",
        "rajat20",
        "rajat24",
        "rajat25",
        "rajat28",
        "rajat30",
        "tsopf_rs_b2383",
        "tsopf_rs_b2383_c1",
    },
    "NICSLU": {
        "add20",
        "add32",
        "asic_100k",
        "circuit5m",
        "circuit5m_dc",
        "circuit_1",
        "circuit_2",
        "circuit_3",
        "circuit_4",
        "coupled",
        "dc1",
        "g2_circuit",
        "g3_circuit",
        "hamrle3",
        "hcircuit",
        "mac_econ_fwd500",
        "mc2depi",
        "memchip",
        "onetone1",
        "onetone2",
        "pre2",
        "raj1",
        "rajat03",
        "rajat30",
        "rajat31",
        "trans4",
        "transient",
        "twotone",
    },
    "SubtreeLU": {
        "asic_100k",
        "asic_100ks",
        "asic_320k",
        "asic_320ks",
        "asic_680k",
        "asic_680ks",
        "bcircuit",
        "circuit5m",
        "circuit5m_dc",
        "circuit_4",
        "ckt11752_dc_1",
        "ckt11752_tr_0",
        "dc1",
        "dc2",
        "dc3",
        "freescale1",
        "g2_circuit",
        "g3_circuit",
        "hcircuit",
        "memchip",
        "nxp1",
        "onetone1",
        "onetone2",
        "pre2",
        "raj1",
        "rajat15",
        "rajat16",
        "rajat17",
        "rajat18",
        "rajat20",
        "rajat21",
        "rajat22",
        "rajat23",
        "rajat24",
        "rajat25",
        "rajat26",
        "rajat28",
        "rajat29",
        "rajat30",
        "rajat31",
        "scircuit",
        "ss1",
        "trans4",
        "trans5",
        "transient",
        "twotone",
    },
    "CKTSO visible/supplement public names": {
        "1138_bus",
        "activsg10k",
        "activsg2000",
        "activsg70k",
        "adder_dcop_01",
        "adder_trans_01",
        "bips07_1693",
        "bips07_1998",
        "bips98_1142",
        "bips98_1450",
        "bips98_606",
        "case9",
        "circuit204",
        "fpga_dcop_01",
        "fpga_trans_01",
        "freescale2",
        "fullchip",
        "gemat11",
        "gemat12",
        "hamrle2",
        "htc_336_4438",
        "htc_336_9129",
        "hvdc1",
        "hvdc2",
        "legresley_2508",
        "legresley_4908",
        "legresley_87936",
        "meg1",
        "memplus",
        "mimo28x28_system",
        "mimo46x46_system",
        "mimo8x8_system",
        "mult_dcop_01",
        "mult_dcop_02",
        "mult_dcop_03",
        "nopss_11k",
        "opf_10000",
        "opf_3754",
        "power197k",
        "powersim",
        "qh1484",
        "rajat12",
        "rajat13",
        "rajat27",
        "tsopf_fs_b39_c19",
        "tsopf_fs_b39_c30",
        "tsopf_fs_b9_c1",
        "tsopf_fs_b9_c6",
        "tsopf_rs_b9_c6",
        "ww_vref_6405",
    },
}


KNOWN_UNRESOLVED_LABELS = {
    "dinesh",
    "hans",
    "leo",
    "steffen",
    "steve_mem",
    "sushil",
    "totyo_bsimcmg",
}


def load_manifest(path: pathlib.Path) -> set[str]:
    names: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            names.add(normalize(text))
    return names


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        type=pathlib.Path,
        default=pathlib.Path("bench/suitesparse_paper_manifest.txt"),
    )
    args = parser.parse_args()

    manifest_names = load_manifest(args.manifest)
    failed = False
    for source, expected in EXPECTED_BY_SOURCE.items():
        missing = sorted({normalize(name) for name in expected} - manifest_names)
        if missing:
            failed = True
            print(f"{source}: missing {', '.join(missing)}", file=sys.stderr)
        else:
            print(f"{source}: covered {len(expected)} names")

    present_unresolved = sorted(KNOWN_UNRESOLVED_LABELS & manifest_names)
    if present_unresolved:
        failed = True
        print(
            "unresolved CKTSO labels unexpectedly present: "
            + ", ".join(present_unresolved),
            file=sys.stderr,
        )

    if failed:
        return 1
    print(f"manifest {args.manifest} covers the audited public paper names")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
