#!/usr/bin/env python3
"""Generate reversible integer-only census copies of the real KLU kernels.

The caller must compile these copies in a diagnostic-only archive. Removing
all inserted spans must recover each exact supported source byte for byte.
Counters include trials and failed calls; they do not establish plain choices.
"""
import argparse
import hashlib
import json
from pathlib import Path


def instrument(source, name):
    if name not in ('refactor', 'tsolve'):
        raise ValueError('unsupported kernel')
    edits = []

    def insert(anchor, addition, count=1):
        nonlocal source
        if source.count(anchor) != count:
            raise ValueError('unexpected source revision: ' + anchor)
        source = source.replace(anchor, anchor + addition)
        edits.append(addition)

    prefix = '\n#include <stdatomic.h>\n#include <stdint.h>\n#include <stdio.h>\n'
    fields = (['calls', 'scaled_calls', 'len0', 'len1', 'len2', 'long', 'entries']
              if name == 'refactor' else
              ['calls', 'one_rhs', 'other_rhs', 'blocks', 'singleton', 'nonsingleton'])
    prefix += 'static _Atomic uint64_t kls_census_' + name + '[' + str(len(fields)) + '];\n'
    prefix += ('#define CENSUS_ADD(i,n) atomic_fetch_add_explicit('
               'kls_census_' + name + '+(i),(uint64_t)(n),memory_order_relaxed)\n')
    prefix += '__attribute__((destructor)) static void census_report(void) {\n'
    prefix += ' fprintf(stderr,"klu_' + name + '_census'
    prefix += ''.join(' ' + f + '=%llu' for f in fields) + '\\n"'
    prefix += ''.join(', (unsigned long long)atomic_load(kls_census_' + name + '+' + str(i) + ')' for i in range(len(fields))) + ');\n}\n'
    insert('#include "trilinos_klu_internal.h"', prefix)
    if name == 'refactor':
        insert('scale = Common->scale ;', '\n    CENSUS_ADD(0,1); CENSUS_ADD(1,scale>0);')
        insert('GET_POINTER (LU, Lip, Llen, Li, Lx, j, llen) ;',
               '\n            CENSUS_ADD(llen<=2?2+llen:5,1); CENSUS_ADD(6,llen);', 2)
    else:
        insert('nblocks = Symbolic->nblocks ;', '\n    CENSUS_ADD(0,1); CENSUS_ADD(nrhs==1?1:2,1);')
        insert('nk = k2 - k1 ;', '\n        CENSUS_ADD(3,1); CENSUS_ADD(nk==1?4:5,1);')
    restored = source
    for addition in reversed(edits):
        restored = restored.replace(addition, '')
    return source, restored


def generate(root, output):
    output.mkdir()
    manifest = {}
    for name in ('refactor', 'tsolve'):
        file = root / 'third_party/suitesparse/KLU/Source' / ('trilinos_klu_' + name + '.c')
        original = file.read_text()
        derived, restored = instrument(original, name)
        if restored != original:
            raise ValueError('restoration mismatch')
        target = output / file.name
        target.write_text(derived)
        manifest[name] = {'input': str(file), 'input_sha256': hashlib.sha256(file.read_bytes()).hexdigest(),
                          'derived': str(target), 'derived_sha256': hashlib.sha256(target.read_bytes()).hexdigest(),
                          'byte_restores_original': True}
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    generate(args.source.resolve(), args.output.resolve())
