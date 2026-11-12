#!/usr/bin/env python3
"""Generate parametric SPICE deck families for size sweeps.

Each family produces the same circuit structure at geometrically increasing
sizes.  The decks are not benchmark targets; they exist so solver policy can
be checked for properties that hold on any workload:

  * linear-solver time grows smoothly with size (no cliff between neighbours),
  * renaming nodes does not change the chosen path (`--shuffle`),
  * more threads are never slower than one.

Families
  rcmesh     2-D resistor grid with a capacitor at every node (linear, mesh fill)
  rcladder   long RC line (linear, tridiagonal, no fill)
  invchain   CMOS inverter chain driven by a pulse (nonlinear, banded)
  ringarray  ring oscillators coupled through a resistive supply grid
  memarray   word/bit-line array of access transistors and storage capacitors
             (a few very dense rows and columns, like RAM decks)
"""
import argparse
import json
import math
from pathlib import Path
import random

MODELS = """\
.model nch nmos level=1 vto=0.7 kp=1.1e-4 gamma=0.4 phi=0.65 lambda=0.02 cgso=2e-10 cgdo=2e-10 cbd=2e-15 cbs=2e-15
.model pch pmos level=1 vto=-0.7 kp=4.0e-5 gamma=0.4 phi=0.65 lambda=0.02 cgso=2e-10 cgdo=2e-10 cbd=2e-15 cbs=2e-15
"""


class Deck:
    def __init__(self, title, shuffle_seed=None):
        self.lines = [title]
        self.names = {}
        self.random = random.Random(shuffle_seed) if shuffle_seed is not None else None
        self.count = 0

    def node(self, name):
        """Map a structural node name to a netlist name; shuffled when requested."""
        if name == '0':
            return '0'
        if name not in self.names:
            self.names[name] = f'n{self.random.getrandbits(40):010x}' if self.random else name
        return self.names[name]

    def add(self, kind, nodes, tail):
        self.count += 1
        self.lines.append(f'{kind}{self.count} ' + ' '.join(self.node(n) for n in nodes) + ' ' + tail)

    def raw(self, text):
        self.lines.append(text)

    def finish(self, stop, step, probes):
        self.raw(f'.tran {step:g} {stop:g}')
        self.raw('.print tran ' + ' '.join(f'v({self.node(p)})' for p in probes))
        self.raw('.end')
        if self.random:
            head, body, tail = self.lines[:1], self.lines[1:-3], self.lines[-3:]
            self.random.shuffle(body)
            self.lines = head + body + tail
        return '\n'.join(self.lines) + '\n'


def rcmesh(size, seed):
    side = max(2, round(math.sqrt(size)))
    d = Deck(f'rcmesh {side}x{side}', seed)
    d.add('V', ['in', '0'], 'pulse(0 1 0 1n 1n 20n 50n)')
    d.add('R', ['in', 'm0_0'], '10')
    for i in range(side):
        for j in range(side):
            d.add('C', [f'm{i}_{j}', '0'], '1p')
            if i + 1 < side:
                d.add('R', [f'm{i}_{j}', f'm{i+1}_{j}'], '100')
            if j + 1 < side:
                d.add('R', [f'm{i}_{j}', f'm{i}_{j+1}'], '100')
    d.add('R', [f'm{side-1}_{side-1}', '0'], '1k')
    return d.finish(200e-9, 0.2e-9, ['m0_0', f'm{side//2}_{side//2}', f'm{side-1}_{side-1}'])


def rcladder(size, seed):
    d = Deck(f'rcladder {size}', seed)
    d.add('V', ['in', '0'], 'pulse(0 1 0 1n 1n 20n 50n)')
    previous = 'in'
    for i in range(size):
        d.add('R', [previous, f'l{i}'], '10')
        d.add('C', [f'l{i}', '0'], '0.1p')
        previous = f'l{i}'
    d.add('R', [previous, '0'], '1k')
    return d.finish(200e-9, 0.2e-9, ['l0', f'l{size//2}', f'l{size-1}'])


def inverter(d, source, sink, supply):
    d.add('M', [sink, source, '0', '0'], 'nch l=1u w=4u')
    d.add('M', [sink, source, supply, supply], 'pch l=1u w=10u')
    d.add('C', [sink, '0'], '20f')


def invchain(size, seed):
    d = Deck(f'invchain {size}', seed)
    d.raw(MODELS.rstrip())
    d.add('V', ['vdd', '0'], '3.3')
    d.add('V', ['in', '0'], 'pulse(0 3.3 1n 0.2n 0.2n 10n 20n)')
    previous = 'in'
    for i in range(size):
        inverter(d, previous, f'c{i}', 'vdd')
        previous = f'c{i}'
    return d.finish(60e-9, 0.05e-9, ['c0', f'c{size//2}', f'c{size-1}'])


def ringarray(size, seed):
    stages = 5
    rings = max(2, size // (stages + 1))
    side = max(1, round(math.sqrt(rings)))
    d = Deck(f'ringarray {rings}x{stages}', seed)
    d.raw(MODELS.rstrip())
    d.add('V', ['vdd', '0'], 'pulse(0 3.3 0 1n 1n 1 2)')
    for r in range(rings):
        supply = f'p{r}'
        # Local supply nodes form a resistive grid fed from one corner.
        d.add('C', [supply, '0'], '50f')
        if r % side:
            d.add('R', [supply, f'p{r-1}'], '5')
        if r >= side:
            d.add('R', [supply, f'p{r-side}'], '5')
        for s in range(stages):
            inverter(d, f'r{r}_{(s - 1) % stages}', f'r{r}_{s}', supply)
        # Stagger the start so the rings are not identical.
        d.add('I', ['0', f'r{r}_0'], f'pulse(0 {20 + (r % 7) * 5}u {1 + (r % 5) * 0.2}n 0.1n 0.1n 1n 1)')
    d.add('R', ['vdd', 'p0'], '1')
    return d.finish(40e-9, 0.05e-9, ['r0_0', f'r{rings//2}_0', f'r{rings-1}_0'])


def memarray(size, seed):
    side = max(2, round(math.sqrt(size / 1.0)))
    d = Deck(f'memarray {side}x{side}', seed)
    d.raw(MODELS.rstrip())
    for i in range(side):
        d.add('V', [f'wd{i}', '0'], f'pulse(0 3.3 {2 + (i % 8) * 4}n 0.2n 0.2n 3n 40n)')
        d.add('R', [f'wd{i}', f'w{i}'], '50')
        d.add('C', [f'w{i}', '0'], '10f')
    for j in range(side):
        d.add('V', [f'bd{j}', '0'], f'pulse(0 3.3 {1 + (j % 4) * 8}n 0.2n 0.2n 6n 40n)')
        d.add('R', [f'bd{j}', f'b{j}'], '200')
        d.add('C', [f'b{j}', '0'], '30f')
    for i in range(side):
        for j in range(side):
            d.add('M', [f'b{j}', f'w{i}', f's{i}_{j}', '0'], 'nch l=1u w=2u')
            d.add('C', [f's{i}_{j}', '0'], '15f')
    return d.finish(80e-9, 0.05e-9, ['s0_0', f's{side//2}_{side//2}', f'b{side-1}'])


FAMILIES = dict(rcmesh=rcmesh, rcladder=rcladder, invchain=invchain, ringarray=ringarray, memarray=memarray)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--family', action='append', choices=sorted(FAMILIES), default=[])
    p.add_argument('--min-log2', type=int, default=8)
    p.add_argument('--max-log2', type=int, default=17)
    p.add_argument('--shuffle', type=int, metavar='SEED',
                   help='Also write a node-renamed, line-shuffled twin of every deck')
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    cases = []
    for family in args.family or sorted(FAMILIES):
        for exponent in range(args.min_log2, args.max_log2 + 1):
            for seed in (None,) + ((args.shuffle,) if args.shuffle is not None else ()):
                name = f'{family}_{exponent:02d}' + ('' if seed is None else '_shuffled') + '.cir'
                (args.output / name).write_text(FAMILIES[family](2 ** exponent, seed))
                cases.append(dict(case=name, family=family, log2_size=exponent, shuffled=seed is not None))
    (args.output / 'family_manifest.json').write_text(json.dumps(dict(schema=1, cases=cases), indent=1) + '\n')
    print(f'wrote {len(cases)} decks to {args.output}')


if __name__ == '__main__':
    main()
