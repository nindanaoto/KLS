#!/usr/bin/env python3
"""Ratchet on matrix-shape gates in the solver.

A shape gate is a branch that compares a matrix dimension, nonzero count, block
count or thread count with a fixed number.  Such a window says nothing about
cost; it only records where a benchmark happened to sit.  New decisions must
come from a cost estimate (flops, factor entries, measured time) or a
family-calibrated tuning entry, so this audit:

  * fails when a gate appears that is not in the allowlist, and
  * fails when an allowlisted gate no longer exists (or occurs less often),
    until the allowlist is regenerated with --update.

The allowlist can therefore only shrink through deliberate edits.  Entries are
keyed by file and normalized expression, not by line number.
"""
import argparse
import json
from pathlib import Path
import re
import sys

SHAPE = r'(?:n|nnz|nz|nblocks|ncols|nrows|factor_entries|lnz|unz|anz|threads|nthreads)'
VAR = rf'(?<![\w.>])(?:[A-Za-z_]\w*(?:\[[^\]]*\])?(?:->|\.))*{SHAPE}\b'
CAST = r'(?:\(\s*[A-Za-z_][\w ]*\*?\s*\)\s*)?'
NUMBER = r'(\d+(?:\.\d*)?(?:[eE][+-]?\d+)?)[uUlLfF]*\b'
SIZED_NAME = (r'((?:KLS_[A-Z0-9_]*(?:MIN|MAX)_(?:SIZE|ROWS|COLS|NNZ|N|DIM|BLOCKS)\b)|'
              r'(?:(?:[A-Za-z_]\w*(?:->|\.))*tuning\.\w*(?:min|max)_(?:size|rows|cols|nnz|n|blocks)\b))')
OP = r'(<=|>=|<|>|==|!=)'
LITERAL = re.compile(rf'({VAR})\s*{OP}\s*{CAST}{NUMBER}')
REVERSED = re.compile(rf'(?<![\w.]){NUMBER}\s*{OP}\s*({VAR})')
NAMED = re.compile(rf'({VAR})\s*{OP}\s*{CAST}{SIZED_NAME}')
COMMENT = re.compile(r'/\*.*?\*/|//[^\n]*', re.S)


def structural(value, variable, operator):
    """Comparisons that separate empty or serial cases are not windows."""
    if float(value) <= 1.0:
        return True
    serial = variable.endswith('threads') and float(value) == 2.0
    return serial and operator in ('<', '>=')


def scan(root):
    found = {}
    for path in sorted([*root.glob('src/*.c'), *root.glob('src/*.inc'), *root.glob('src/*.h')]):
        text = COMMENT.sub(lambda m: re.sub(r'[^\n]', ' ', m[0]), path.read_text(errors='replace'))
        name = str(path.relative_to(root))
        for pattern, value_group, variable_group, flip in ((LITERAL, 3, 1, False), (REVERSED, 1, 3, True)):
            for match in pattern.finditer(text):
                operator = match[2]
                if flip:
                    operator = {'<': '>', '>': '<', '<=': '>=', '>=': '<='}.get(operator, operator)
                if structural(match[value_group], match[variable_group], operator):
                    continue
                key = name + ': ' + re.sub(r'\s+', '', match[0])
                found[key] = found.get(key, 0) + 1
        for match in NAMED.finditer(text):
            key = name + ': ' + re.sub(r'\s+', '', match[0])
            found[key] = found.get(key, 0) + 1
    return found


def main():
    root = Path(__file__).resolve().parents[1]
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=root)
    p.add_argument('--allowlist', type=Path, default=root / 'bench/shape_gate_allowlist.json')
    p.add_argument('--update', action='store_true',
                   help='Rewrite the allowlist; refuses to add gates or raise counts')
    p.add_argument('--initialize', action='store_true', help='Create the first allowlist')
    args = p.parse_args()
    found = scan(args.root)
    if args.initialize:
        if args.allowlist.exists():
            p.error('allowlist already exists; it may only shrink (--update)')
        args.allowlist.write_text(json.dumps(dict(schema=1, gates=dict(sorted(found.items()))), indent=1) + '\n')
        print(f'initialized {len(found)} gates, {sum(found.values())} occurrences')
        return 0
    allowed = json.loads(args.allowlist.read_text())['gates']
    added = {k: v for k, v in found.items() if v > allowed.get(k, 0)}
    stale = {k: v for k, v in allowed.items() if found.get(k, 0) < v}
    if added:
        print('New shape gates (replace with a cost model or a family-calibrated tuning entry):')
        for key, count in sorted(added.items()):
            print(f'  {key}  x{count} (allowed {allowed.get(key, 0)})')
        return 1
    if args.update:
        args.allowlist.write_text(json.dumps(dict(schema=1, gates=dict(sorted(found.items()))), indent=1) + '\n')
        print(f'allowlist now {len(found)} gates, {sum(found.values())} occurrences')
        return 0
    if stale:
        print('Allowlisted gates that were removed; shrink the allowlist with --update:')
        for key, count in sorted(stale.items()):
            print(f'  {key}  allowed x{count}, found x{found.get(key, 0)}')
        return 1
    print(f'{len(found)} shape gates, {sum(found.values())} occurrences; none new')
    return 0


if __name__ == '__main__':
    sys.exit(main())
