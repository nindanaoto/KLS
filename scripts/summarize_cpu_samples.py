#!/usr/bin/env python3
"""Resolve profile_cpu_samples.c instruction samples against an unstripped ELF."""
import argparse
import bisect
import collections
import hashlib
import json
from pathlib import Path
import re
import subprocess


def summarize(executable, samples):
    records = [json.loads(line) for line in samples.read_text().splitlines()]
    header, *records = records
    if (header['version'] != 1 or header['saturated'] or
            len(records) != header['samples'] or not records):
        raise ValueError('Missing, incomplete, empty or saturated sample capture')
    base = int(header['executable_base'], 16)
    elf = subprocess.check_output(['readelf', '-h', str(executable)], text=True)
    kind = re.search(r'^\s*Type:\s+(\w+)', elf, re.M).group(1)
    if kind not in ('EXEC', 'DYN'):
        raise ValueError('Expected an executable ELF')
    symbols = []
    nm = subprocess.check_output(
        ['nm', '-n', '-S', '--defined-only', str(executable)], text=True)
    for line in nm.splitlines():
        match = re.fullmatch(r'([\da-fA-F]+) ([\da-fA-F]+) [tTwW] (.+)', line)
        if match:
            start, size, name = match.groups()
            if int(size, 16):
                symbols.append((int(start, 16), int(size, 16), name))
    symbols.sort()
    starts = [item[0] for item in symbols]
    maps_path = Path(str(samples) + '.maps')
    mappings = []
    for line in maps_path.read_text().splitlines():
        parts = line.split(maxsplit=5)
        start, end = (int(v, 16) for v in parts[0].split('-'))
        path = parts[5] if len(parts) == 6 else '[anonymous]'
        mappings.append((start, end, parts[1], path))
    # dladdr supplies the executable's mapping base, even for non-PIE ELFs.
    own_paths = {path for start, end, _, path in mappings if start <= base < end}
    if len(own_paths) != 1:
        raise ValueError('Cannot identify executable mapping from capture base')
    counts = collections.Counter()
    for record in records:
        pc = int(record['pc'], 16)
        mapping = next((m for m in mappings if m[0] <= pc < m[1]), None)
        if mapping is None:
            label = '[unmapped]'
        elif mapping[3] not in own_paths:
            label = '[external] ' + mapping[3]
        else:
            address = pc - base if kind == 'DYN' else pc
            at = bisect.bisect_right(starts, address) - 1
            if at >= 0 and address < symbols[at][0] + symbols[at][1]:
                label = symbols[at][2]
            else:
                label = '[unresolved executable]'
        counts[label] += 1
    return {
        'samples': len(records), 'period_us': header['period_us'],
        'elf_type': kind,
        'sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                   for p in (executable, samples, maps_path)},
        'functions': [{'name': name, 'samples': n,
                       'fraction': n / len(records)}
                      for name, n in counts.most_common()],
        'scope': 'Process-CPU instruction samples, including system time; '
                 'signals can coalesce. Counts locate hot code and do not '
                 'measure wall time or establish a performance improvement.'
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('samples', type=Path)
    args = parser.parse_args()
    print(json.dumps(summarize(args.executable, args.samples), indent=2))
