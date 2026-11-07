#!/usr/bin/env python3
"""Counterbalanced direct-lifecycle comparisons of frozen KLS binaries."""
from __future__ import annotations
import argparse
import fcntl
import hashlib
import importlib.util
import json
import os
import signal
import subprocess
import time
from pathlib import Path
from run_paper_campaign import clean_environment, command_for, manifest_names, resolve_matrices, parse_json_output
from generalization_gate import evaluate, POLICY, PEER_POLICY

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def tuning_profiles(configs, arms):
    paths=set()
    for config in configs:
        for arm in arms:
            if arm in ('ck','st'):continue
            arguments=config.get(arm+'_args',config.get('kls_args',[]))
            if not isinstance(arguments,list) or any(not isinstance(v,str) for v in arguments):
                raise ValueError('solver arguments must be a list of strings')
            for i,value in enumerate(arguments):
                if value=='--tuning-profile':
                    if i+1==len(arguments):raise ValueError('missing tuning profile path')
                    paths.add(Path(arguments[i+1]).resolve(strict=True))
                elif value.startswith('--tuning-profile='):
                    paths.add(Path(value.split('=',1)[1]).resolve(strict=True))
    return paths

def launch_order(arms,pass_number,confirmation=False):
    # Adjacent reverse pairs balance every comparison, including five arms.
    shift=(pass_number//2)%len(arms)
    order=arms[shift:]+arms[:shift]
    return list(reversed(order)) if bool(pass_number%2)!=confirmation else order

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('baseline','candidate','previous','image','manifest','matrix-dir','configs','output'):
        p.add_argument('--'+name,type=Path,required=name!='previous')
    p.add_argument('--ck',type=Path)
    p.add_argument('--st',type=Path)
    p.add_argument('--peer-library',type=Path,action='append',default=[],help='Freeze external peer libraries and license inputs')
    p.add_argument('--acceptance-policy',choices=[POLICY,PEER_POLICY],default=PEER_POLICY)
    p.add_argument('--passes',type=int,default=8)
    p.add_argument('--timeout',type=float,default=300)
    p.add_argument('--order',choices=['initial','confirmation'],default='initial')
    a=p.parse_args()
    if a.passes<8 or a.passes%4 or a.timeout<=0:p.error('use at least eight passes in blocks of four and a positive timeout')
    if bool(a.ck)!=bool(a.st):p.error('supply both --ck and --st')
    if a.acceptance_policy==PEER_POLICY and not (a.ck and a.st and a.peer_library):
        p.error('generalization-v2 requires --ck, --st, and frozen --peer-library inputs')
    root=Path(__file__).resolve().parents[1];workspace=root.parent
    binaries={'baseline':a.baseline.resolve(),'candidate':a.candidate.resolve()}
    if a.previous:binaries['previous']=a.previous.resolve()
    if a.ck:binaries.update(ck=a.ck.resolve(),st=a.st.resolve())
    names=manifest_names(a.manifest);matrices=resolve_matrices(a.matrix_dir,names)
    configs=json.loads(a.configs.read_text())
    if not configs or len({c['id'] for c in configs})!=len(configs):p.error('nonempty, unique configurations required')
    if any(int(c['threads'])<=0 or int(c['lifecycle_systems'])<=0 for c in configs):p.error('invalid lifecycle/thread count')
    cleanup_path=root/'bench/process_session_cleanup.py'
    loader=importlib.util.spec_from_file_location('generalization_cleanup',cleanup_path)
    cleanup=importlib.util.module_from_spec(loader);loader.loader.exec_module(cleanup)
    scripts=[Path(__file__),root/'scripts/generalization_gate.py',root/'scripts/run_paper_campaign.py',cleanup_path]
    profiles=sorted(tuning_profiles(configs,binaries),key=str)
    external=sorted({p.resolve(strict=True) for p in a.peer_library},key=str)
    hashes={str(b):sha(b) for b in [*binaries.values(),*scripts,*profiles,*external]}
    spec=dict(policy=a.acceptance_policy,arms=list(binaries),configs=configs,matrices=names,
              passes=a.passes,order=a.order,timeout=a.timeout,
              binaries={arm:{'path':str(b),'sha256':sha(b)} for arm,b in binaries.items()},
              scripts={str(s):sha(s) for s in scripts},image_sha256=sha(a.image),
              profile_hashes={str(p):sha(p) for p in profiles},
              external_hashes={str(p):sha(p) for p in external},
              matrix_hashes={name:sha(path) for name,path in matrices.items()})
    lock=(workspace/'.spice-timing.lock').open('a+')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    # Timing must not overlap builds, tests or another solver campaign.
    active=[]
    for line in subprocess.check_output(['ps','-eo','comm='],text=True).splitlines():
        if line.strip() in {'cc1','cc1plus','f951','cmake','ctest','kls_bench','cktso_compare','subtreelu_compare','subtreelu_compa','Xyce'}:active.append(line.strip())
    if active:raise RuntimeError('timing rejected while build/test/solver processes are active: '+str(active))
    a.output.mkdir(parents=True,exist_ok=False)
    (a.output/'campaign.json').write_text(json.dumps(spec,indent=2)+'\n')
    rows=[]
    with (a.output/'observations.jsonl').open('x') as ledger:
        for name in (list(reversed(names)) if a.order=='confirmation' else names):
            if sha(matrices[name])!=spec['matrix_hashes'][name]:raise RuntimeError('frozen matrix changed: '+name)
            for config in (list(reversed(configs)) if a.order=='confirmation' else configs):
                for pass_number in range(a.passes):
                    order=launch_order(list(binaries),pass_number,a.order=='confirmation')
                    for arm in order:
                        if any(sha(path)!=digest for path,digest in hashes.items()):raise RuntimeError('frozen timing input changed')
                        side_config=dict(config)
                        side_config['kls_args']=config.get(arm+'_args',config.get('kls_args',[]))
                        command=command_for(arm if arm in ('ck','st') else 'kls',binaries[arm],matrices[name],side_config)
                        command=['apptainer','exec','--bind',f'{workspace}:{workspace}',str(a.image.resolve()),*command]
                        if config.get('cpus'):command=['taskset','-c',config['cpus'],*command]
                        env=clean_environment(config.get('env',{}));start=time.monotonic()
                        child=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,env=env,start_new_session=True)
                        timed_out=False
                        try:out,err=child.communicate(timeout=a.timeout)
                        except subprocess.TimeoutExpired:
                            timed_out=True;os.killpg(child.pid,signal.SIGKILL);out,err=child.communicate()
                        closure=cleanup.drain_exited_apptainer_session(child.pid,a.image)
                        try:
                            result=parse_json_output(out.decode())
                            # Keep nonfinite peer diagnostics in the ledger as
                            # invalid evidence rather than aborting the campaign.
                            result=json.loads(json.dumps(result),parse_constant=lambda value:value)
                        except (ValueError,UnicodeDecodeError):result=None
                        row=dict(config=config['id'],matrix=name,pass_number=pass_number,arm=arm,
                                 returncode=child.returncode,timed_out=timed_out,result=result,closure=closure,
                                 wall_seconds=time.monotonic()-start,command=command,stderr=err.decode(errors='replace')[-4000:])
                        row['pass']=row.pop('pass_number');rows.append(row)
                        ledger.write(json.dumps(row,allow_nan=False)+'\n');ledger.flush()
                print(config['id'],name,flush=True)
            if sha(matrices[name])!=spec['matrix_hashes'][name]:raise RuntimeError('frozen matrix changed during execution: '+name)
    result=evaluate(rows,spec);(a.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print('gate',result['passed'],flush=True)
    return 0 if result['passed'] else 1
if __name__=='__main__':raise SystemExit(main())
