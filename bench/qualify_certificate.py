"""Resumable, serialized SPICE qualification of a frozen certificate candidate.

Requires the paper's SPICE support directory, an input inventory, and two built
Xyce executables. Nothing is inferred from /tmp. Accuracy capture timings are
never included in performance results. Partial/failed observations are retained
and cannot be silently retried or accepted. Use a new output root for a retry.
"""
import argparse
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import statistics
import sys
import tarfile

PRIMARY = ('PDE_2D_BJT_invertbjt_cir', 'MOS13_IC_invert50_mos1_cir',
           'mos2_large_mem_plus')
CONTROLS = ('mos2_mux8', 'mos2_large_dac', 'neuron_rallpack3_1000level1s',
            'mos2_large_ram2k')
GIB = 1024**3
POLICY = 'KLS_ACCURACY_POLICY=COMPONENTWISE_BACKWARD_ERROR'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def summary(rows, phase, incremental):
    result = []
    for case in PRIMARY + CONTROLS:
        for mode in ('kls1', 'kls8'):
            selected = [r for r in rows if r['campaign_phase'] == phase
                        and r['circuit'] == case and r['mode'] == mode]
            ratios = []
            for block in sorted({r['block'] for r in selected}):
                pair = {r['variant']: r['xyce_elapsed_s'] for r in selected
                        if r['block'] == block}
                require(set(pair) == {'baseline', 'candidate'}, 'incomplete timing pair')
                require(all(math.isfinite(t) and t > 0 for t in pair.values()), 'invalid elapsed time')
                ratios.append(pair['candidate'] / pair['baseline'])
            require(len(ratios) == (3 if phase == 'screen' else 9), 'missing timing blocks')
            result.append(dict(case=case, mode=mode, paired_ratios=ratios,
                               median_ratio=statistics.median(ratios)))
    no_regression = all(r['median_ratio'] <= 1.02 for r in result)
    targets = [r for r in result if r['case'] in PRIMARY]
    benefit = (any(r['median_ratio'] <= .98 for r in targets) if incremental else
               all(r['median_ratio'] <= .90 for r in targets))
    return dict(passed=no_regression and benefit, no_regression=no_regression,
                benefit=benefit, results=result)


def audit_chunk(args):
    import numpy as np
    from scipy.io import mmread
    from audit_backward_error import exact_errors, sha256
    from qualify_recovered_spice import load_system, system_key
    from discovery_support import write_json
    directory = args.directory
    rejected = json.loads((directory / 'rejections.json').read_text())
    by_key = {r['system_key']: r['id'] for r in rejected}
    require(len(by_key) == len(rejected), 'ambiguous rejected systems')
    records = []
    for epoch in range(args.first, args.last + 1):
        paths = [directory / f'Transformed_{kind}{epoch}.mm'
                 for kind in ('Matrix', 'RHS', 'Soln')]
        a, b = load_system(*paths[:2])
        # The exported caller matrices are unique-coordinate matrices. Refuse
        # silent rounded duplicate coalescing in the independent audit reader.
        raw = mmread(paths[0])
        require(raw.nnz == a.nnz, 'duplicate coordinates in capture require exact coalescing')
        x = np.asarray(mmread(paths[2]), dtype=float).reshape(-1)
        require(len(x) == len(b) and np.all(np.isfinite(x)), 'invalid captured solution')
        key = system_key(a, b)
        rejection = by_key.get(key)
        certificate = None
        if rejection is not None:
            require(np.all(x == 0), 'nonzero rejection export')
            status = 'explicit_rejection'
        else:
            metrics = exact_errors(a, b, x)
            passed = Fraction(int(metrics['componentwise_numerator']),
                              int(metrics['componentwise_denominator'])) <= Fraction(1, 10**8)
            certificate = dict(passed=passed, metrics=metrics)
            status = 'certified_success' if passed else 'uncertified_success'
        records.append(dict(epoch=epoch, status=status, rejection_id=rejection,
                            system_key=key, certificate=certificate,
                            hashes={p.name: sha256(p) for p in paths}))
    write_json(args.output, records)
    require(all(r['status'] != 'uncertified_success' for r in records),
            'uncertified successful answer')


def archive_chunk(directory, records, archive, receipt):
    """Verify a small archive before reclaiming its explicitly listed raw files."""
    from discovery_support import sha256, write_json
    hashes = {name: digest for r in records for name, digest in r['hashes'].items()}
    def verify():
        seen = set()
        with tarfile.open(archive, 'r|gz') as stream:
            for member in stream:
                require(member.isfile() and member.name in hashes and member.name not in seen,
                        'unexpected archive member')
                h = hashlib.sha256()
                with stream.extractfile(member) as data:
                    for block in iter(lambda: data.read(1024*1024), b''):
                        h.update(block)
                require(h.hexdigest() == hashes[member.name], 'archive hash mismatch')
                seen.add(member.name)
        require(seen == set(hashes), 'incomplete archive')
    if receipt.exists():
        old = json.loads(receipt.read_text())
        require(old['files'] == hashes and old['archive'] == str(archive)
                and old['sha256'] == sha256(archive), 'changed archived chunk')
    else:
        require(not archive.exists(), 'unreceipted archive: preserve and investigate')
        archive.parent.mkdir(parents=True, exist_ok=True)
        with tarfile.open(archive, 'w:gz', compresslevel=1) as stream:
            for name, expected in hashes.items():
                path = directory / name
                require(path.is_file() and not path.is_symlink() and sha256(path) == expected,
                        'raw capture changed')
                stream.add(path, arcname=name, recursive=False)
        verify()
        write_json(receipt, dict(archive=str(archive), sha256=sha256(archive),
                                files=hashes, verified=True, raw_removed=False))
    # Recheck on resumption, including after a crash during raw reclamation.
    verify()
    for name, expected in hashes.items():
        path = directory / name
        if path.exists():
            require(not path.is_symlink() and sha256(path) == expected, 'raw capture changed')
            path.unlink()
    data = json.loads(receipt.read_text())
    data['raw_removed'] = True
    write_json(receipt, data)


def campaign(args):
    from discovery_support import bounded_run, compare_file, run_observation, sha256, write_json
    from qualify_recovered_spice import load_system, numbered, system_key, validate_coverage
    args.output.mkdir(parents=True, exist_ok=True)
    args.scratch.mkdir(parents=True, exist_ok=True)
    require(args.output != args.scratch, 'keep evidence and scratch roots distinct')
    require(args.phase=='accuracy' or (args.cases is None and args.modes is None),
            'performance qualification cannot restrict the case or thread coverage')
    ids = tuple(args.cases or (PRIMARY+CONTROLS if args.phase=='performance' else PRIMARY))
    modes = tuple(args.modes or ('kls1','kls8'))
    require(len(ids)==len(set(ids)) and len(modes)==len(set(modes)), 'duplicate case or mode')
    build = json.loads(args.build_manifest.read_text())
    require(build['candidate_sha256']==sha256(args.candidate) and
            build['baseline_sha256']==sha256(args.baseline), 'build manifest does not match binaries')
    support_files = [args.support / name for name in ('discovery_support.py', 'run_spice.py',
                      'qualify_recovered_spice.py', 'audit_backward_error.py')]
    pinned = [Path(__file__).resolve(), args.inventory, args.baseline, args.candidate,
              args.build_manifest, *support_files]
    identity = dict(schema=1, phase=args.phase, incremental=args.incremental,
                    scratch=str(args.scratch), support=str(args.support),
                    hashes={str(p): sha256(p) for p in pinned},
                    accuracy_contract='componentwise <= 1/100000000',
                    cases=list(ids), modes=list(modes),
                    screen_blocks=3, confirmation_blocks=9)
    manifest = args.output / 'manifest.json'
    if manifest.exists():
        require(json.loads(manifest.read_text()) == identity, 'campaign inputs changed')
    else:
        require(not any(args.output.iterdir()), 'new campaign needs an empty output root')
        require(not any(args.scratch.iterdir()), 'new campaign needs an empty scratch root')
        write_json(manifest, identity)
    cases = {c['id']: c for c in json.loads(args.inventory.read_text())}
    state = dict(complete=False, phase='starting', completed_runs=0)
    status = args.output / 'status.json'
    def save(**updates):
        state.update(updates)
        write_json(status, state)
    def check():
        require(all(sha256(Path(p)) == h for p,h in identity['hashes'].items()),
                'pinned binary or script changed')
        for cid in ids:
            inputs = cases[cid]['inputs']
            require(all(sha256(Path(inputs['root']) / p) == h for p,h in inputs['files'].items()),
                    'circuit input changed')
    def observation(cid, mode, block, phase, binary, extra=''):
        check()
        checkpoint = args.output / 'observations' / f'{phase}-{cid}-{mode}-{block}.json'
        if checkpoint.exists():
            row = json.loads(checkpoint.read_text())
        else:
            run_dir = args.scratch / 'runs' / phase / cid / f'rep-{block:02d}' / mode
            require(not run_dir.exists(), 'interrupted run preserved; use a new root to retry')
            save(phase=phase, case=cid, mode=mode, block=block)
            row = run_observation(cases[cid], mode, block, phase, binary, args.scratch,
                                  1800, 24*GIB, 8*GIB, extra=extra)
            write_json(checkpoint, row)
        require(row['completed'] and row.get('outputs'), f'unsuccessful run: {checkpoint}')
        require(all(sha256(args.scratch / item['path']) == item['sha256']
                    for item in row['outputs'].values()), 'changed waveform output')
        if mode.startswith('kls'):
            log = (args.scratch / row['log']).read_text()
            require('KLS accuracy policy = COMPONENTWISE_BACKWARD_ERROR' in log,
                    'accuracy policy not attested')
            limit = 1 if cid == PRIMARY[1] else 0
            require(isinstance(row['failed_linear_solves'], int) and
                    0 <= row['failed_linear_solves'] <= limit, 'new linear rejection')
        state['completed_runs'] += 1
        return row
    def waves(row, reference):
        require(set(row['outputs']) == set(reference['outputs']), 'waveform coverage differs')
        checks = {name: compare_file(args.scratch / item['path'],
                                    args.scratch / reference['outputs'][name]['path'])
                  for name,item in row['outputs'].items()}
        require(checks and all(c['valid'] for c in checks.values()), 'waveform mismatch')
        return checks
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    try:
        refs = {cid: observation(cid, 'klu', 0, 'reference', args.baseline) for cid in ids}
        if args.phase == 'performance':
            rows = []
            for phase,count in (('screen',3), ('confirmation',9)):
                for block in range(count):
                    order = ids[block % len(ids):] + ids[:block % len(ids)]
                    for ci,cid in enumerate(order):
                        for mi,mode in enumerate(modes):
                            variants = ('baseline','candidate') if (block+ci+mi)%2 == 0 else ('candidate','baseline')
                            for variant in variants:
                                row = observation(cid, mode, block, phase+'-'+variant,
                                                  getattr(args,variant), POLICY)
                                row.update(campaign_phase=phase, block=block, variant=variant,
                                           waveform_checks=waves(row,refs[cid]))
                                rows.append(row)
                                write_json(args.output / 'observations.json', rows)
                                print(phase, block, cid, mode, variant, row['xyce_elapsed_s'], flush=True)
                result = summary(rows,phase,args.incremental)
                write_json(args.output / (phase+'-summary.json'),result)
            save(complete=True, phase='complete', passed=result['passed'])
        else:
            results = []
            for cid in ids:
                for mode in modes:
                    label = cid+'-'+mode
                    row = observation(cid,mode,0,'accuracy',args.candidate,
                                      POLICY+' OUTPUT_LS=1 OUTPUT_FAILED_LS=1')
                    waveform_checks = waves(row,refs[cid])
                    directory = (args.scratch / row['log']).parent
                    expected,failed = row['linear_solves'],row['failed_linear_solves']
                    require(isinstance(expected,int) and expected>failed, 'invalid attempt count')
                    log = (directory / 'stdout.txt').read_text()
                    require(log.count('KLS linear solve failed:')==failed, 'warning count mismatch')
                    preparation = args.output / (label+'-capture.json')
                    if not preparation.exists():
                        for kind in ('Matrix','RHS','Soln'):
                            require(numbered(directory,'Transformed_'+kind)==list(range(1,expected+1)),
                                    'incomplete capture sequence')
                        rejected = []
                        for kind in ('Matrix','RHS'):
                            require(numbered(directory,'Failed_'+kind)==list(range(1,failed+1)),
                                    'incomplete failed sequence')
                        for ident in range(1,failed+1):
                            paths = [directory / f'Failed_{kind}{ident}.mm' for kind in ('Matrix','RHS')]
                            a,b = load_system(*paths)
                            rejected.append(dict(id=ident, system_key=system_key(a,b),
                                                 hashes={p.name:sha256(p) for p in paths}))
                        write_json(directory/'rejections.json', rejected)
                        write_json(preparation,dict(expected=expected, failed=failed,
                                                    rejections_sha256=sha256(directory/'rejections.json')))
                    prep = json.loads(preparation.read_text())
                    require(prep['expected']==expected and prep['failed']==failed and
                            prep['rejections_sha256']==sha256(directory/'rejections.json'),
                            'capture metadata changed')
                    for rejected in json.loads((directory/'rejections.json').read_text()):
                        require(all(sha256(directory/name)==digest for name,digest in rejected['hashes'].items()),
                                'failed-system export changed')
                    records = []
                    for first in range(1,expected+1,20):
                        check()
                        save(phase='exact-audit', case=cid, mode=mode,
                             audited=len(records), expected=expected)
                        chunk = args.output / 'chunks' / f'{label}-{first}.json'
                        chunk.parent.mkdir(exist_ok=True)
                        if not chunk.exists():
                            command = ['taskset','-c','0',sys.executable,str(Path(__file__).resolve()),
                                       'audit-chunk','--support',str(args.support),
                                       '--directory',str(directory),'--first',str(first),
                                       '--last',str(min(first+19,expected)),'--output',str(chunk)]
                            outcome = bounded_run(command,directory,env,120,24*GIB,8*GIB,
                                                  chunk.with_suffix('.stdout'),chunk.with_suffix('.stderr'))
                            require(outcome['returncode']==0 and not outcome['termination'],
                                    f'audit worker failed: {outcome}')
                        part = json.loads(chunk.read_text())
                        require([r['epoch'] for r in part]==list(range(first,min(first+19,expected)+1)),
                                'wrong audit chunk')
                        require(all(r['status'] in ('certified_success','explicit_rejection') and
                                    (r['status']!='certified_success' or r['certificate']['passed'])
                                    for r in part), 'failed audit chunk')
                        receipt = chunk.with_suffix('.archive.json')
                        if receipt.exists():
                            archive = Path(json.loads(receipt.read_text())['archive'])
                        else:
                            # Worst-case uncompressed allowance prevents writing below
                            # either filesystem's reserve. Small verified chunks let
                            # raw reclamation keep pace with compressed storage.
                            size = sum((directory/name).stat().st_size for r in part for name in r['hashes'])
                            archive_root = next((p for p in (args.output,args.scratch)
                                                 if shutil.disk_usage(p).free > 5*GIB+size+1024**2),None)
                            require(archive_root is not None, 'no safe archive space; resume after freeing space')
                            archive = archive_root / 'archives' / (chunk.stem+'.tar.gz')
                        archive_chunk(directory,part,archive,receipt)
                        records.extend(part)
                        print(label,len(records),'/',expected,flush=True)
                    coverage = validate_coverage(records,expected,failed)
                    result = dict(case=cid,mode=mode,passed=True,coverage=coverage,
                                  waveforms=waveform_checks)
                    write_json(args.output/(label+'-audit.json'),result)
                    results.append(result)
            save(complete=True,phase='complete',passed=True,results=results,
                 full_accuracy_suite=set(ids)==set(PRIMARY) and set(modes)=={'kls1','kls8'})
    except BaseException as error:
        save(complete=False,phase='stopped',passed=False,error=repr(error))
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command',required=True)
    run = sub.add_parser('run')
    run.add_argument('--phase',choices=('accuracy','performance'),required=True)
    for name in ('support','inventory','baseline','candidate','build-manifest','output','scratch'):
        run.add_argument('--'+name,type=Path,required=True)
    run.add_argument('--incremental',action='store_true')
    run.add_argument('--cases',nargs='+',choices=PRIMARY,
                     help='accuracy-only retry scope; not a full-suite qualification')
    run.add_argument('--modes',nargs='+',choices=('kls1','kls8'),
                     help='accuracy-only retry scope; performance always tests both modes')
    worker = sub.add_parser('audit-chunk')
    for name in ('support','directory','output'):
        worker.add_argument('--'+name,type=Path,required=True)
    for name in ('first','last'):
        worker.add_argument('--'+name,type=int,required=True)
    args = parser.parse_args()
    for key,value in vars(args).items():
        if isinstance(value,Path):setattr(args,key,value.resolve())
    sys.path.insert(0,str(args.support))
    audit_chunk(args) if args.command=='audit-chunk' else campaign(args)


if __name__=='__main__':main()
