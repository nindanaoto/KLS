"""Link frozen main/candidate cores and qualify unhinted API and owned execution.

All times collected here are diagnostic. Ordinary performance campaigns use
candidate-link/Xyce only, never activation/Xyce.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import signal
import subprocess
import sys


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def save(path,value):
    Path(path).write_text(json.dumps(value,indent=2)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('workspace','candidate-core','baseline-core','output'):
        p.add_argument('--'+name,required=True,type=Path)
    a=p.parse_args();workspace=a.workspace.resolve();out=a.output.resolve()
    cores={'candidate':a.candidate_core.resolve(),'baseline':a.baseline_core.resolve()}
    repo=Path(__file__).resolve().parents[1];image=workspace/'KLS/build/kls-xyce-private.sif'
    pins={str(path):sha(path) for path in (Path(__file__),image,
        repo/'scripts/build_spice_links.py',repo/'scripts/run_spice_campaign.py',
        repo/'scripts/audit_spice_campaign.py',repo/'scripts/auto_btf_merge_gate.py')}
    for side,core in cores.items():
        done=read(core/'completed.json');source=read(core/'source.json')
        assert done['passed'] and 'native' in done['profiles']
        if side=='candidate':assert done['profiles']==['native','portable','asan','tsan']
        assert done['source_sha256']==sha(core/'source.json')
        assert done['commands_sha256']==sha(core/'commands.json')
        assert done['libraries']['native']==sha(core/'native/libkls.a')
        pins.update(source['input_sha256'])
        for path in (core/'completed.json',core/'source.json',core/'commands.json'):
            pins[str(path)]=sha(path)
        for path in list((core/'native').glob('*.a'))+list((core/'native/lib').glob('*.a')):
            pins[str(path)]=sha(path)
    headers=lambda core:{str(p.relative_to(core/'source/include')):sha(p) for p in (core/'source/include').rglob('*.h')}
    assert headers(cores['candidate'])==headers(cores['baseline'])
    assert read(cores['candidate']/'source.json')['dependencies']==read(cores['baseline']/'source.json')['dependencies']
    out.mkdir(exist_ok=False);save(out/'inputs.json',dict(diagnostic_only=True,input_sha256=pins))
    spec=importlib.util.spec_from_file_location('cleanup',repo/'bench/process_session_cleanup.py')
    cleanup=importlib.util.module_from_spec(spec);spec.loader.exec_module(cleanup)
    env={k:v for k,v in os.environ.items() if not k.startswith(('LD_','KLS_','OMP_','GOMP_',
        'OPENBLAS_','APPTAINERENV_','ASAN_','UBSAN_','TSAN_'))}
    env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    prefix=['apptainer','exec','--bind',f'{workspace}:/work','--bind',f'{workspace}:{workspace}',
            '--bind','/tmp:/tmp',str(image)]
    commands=[]

    def verify():
        assert all(Path(path).exists() and sha(path)==digest for path,digest in pins.items())

    def run(stage,command,timeout=7200):
        verify();command=list(map(str,command));print('START',stage,flush=True)
        with (out/(stage+'.log')).open('x') as log:
            child=subprocess.Popen(command,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            try:code=child.wait(timeout=timeout)
            except BaseException:
                os.killpg(child.pid,signal.SIGKILL);child.wait()
                cleanup.drain_exited_apptainer_session(child.pid,image);raise
        closure=cleanup.drain_exited_apptainer_session(child.pid,image)
        commands.append(dict(stage=stage,command=command,exit_code=code,closure=closure))
        save(out/'commands.json',commands);print('END',stage,code,flush=True)
        assert code==0 and closure['all_session_processes_absent'],stage
        verify()

    for side,core in cores.items():
        dest=out/(side+'-link')
        run(side+'-link',[sys.executable,repo/'scripts/build_spice_links.py','--workspace',workspace,
            '--image',image,'--artifacts',workspace/'KLS/build/xyce-testing-bump-20261023/artifacts',
            '--provenance',workspace/'KLS/build/xyce-testing-bump-20261023/build-provenance',
            '--generated','/tmp/kls-xyce-testing-20261023','--trilinos','/tmp/kls-trilinos-install',
            '--xyce-source',workspace/'Xyce-spice-fast','--kls-source',core/'source',
            '--kls-build',core/'native','--dependency-build',core/'native',
            '--common-archive','/tmp/kls-spice-20261027/control-link-v3/libxyce.a','--output',dest])
        build=read(dest/'build.json');pins.update(build['input_sha256'])
        pins.update(build['output_sha256'])
        for name in ('build.json','commands.json'):pins[str(dest/name)]=sha(dest/name)
        save(out/'inputs.json',dict(diagnostic_only=True,input_sha256=pins))

    def campaign(label,binary,provenance):
        dest=out/label;audit=out/(label+'-independent')
        run(label,[sys.executable,repo/'scripts/run_spice_campaign.py','--workspace',workspace,
            '--regression','/tmp/kls-xyce-regression','--image',image,'--license','/tmp/kls-goal-cktso/cktso.lic',
            '--baseline',out/'baseline-link/Xyce','--candidate',binary,'--provenance',provenance,
            '--diagnostic','--passes','1','--output',dest])
        run(label+'-review',[sys.executable,repo/'scripts/run_spice_campaign.py','--review',dest])
        run(label+'-audit',[sys.executable,repo/'scripts/audit_spice_campaign.py','--campaign',dest,'--output',audit])
        r=read(dest/'review.json');q=read(audit/'audit.json')
        assert r['all_qualified'] and q['all_qualified'] and r['observations']==q['observations']==70
        key=lambda c:(c['domain'],c['case'],c['threads'],c['peer'])
        assert len(r['cells'])==len(q['cells'])==28
        assert {key(c):(c['median'],c['wall_median']) for c in r['cells']}=={key(c):(c['median'],c['wall_median']) for c in q['cells']}
        save(out/(label+'-agreement.json'),dict(all70_qualified=True,independent_exact_agreement=True))
        return dest

    campaign('api',out/'candidate-link/Xyce',out/'candidate-link/build.json')
    core=cores['candidate'];diagnostic=out/'activation-source'
    shutil.copytree(core/'source/src',diagnostic);insertions={}

    def insert(name,needle,added):
        path=diagnostic/name;before=path.read_text();assert before.count(needle)==1,(name,needle)
        path.write_text(before.replace(needle,needle+added,1));insertions.setdefault(name,[]).append(added)

    insert('kls.c','struct kls_solver {','\n  uint64_t diagnostic_btf_factor,diagnostic_btf_solve,diagnostic_btf_direct,diagnostic_btf_accepted;\n')
    insert('kls.c','#include <stdatomic.h>', '\n#define KLS_COMPACT_BTF_FACTOR_DIAGNOSTIC_EVENT(s) (++(s)->diagnostic_btf_factor)\n'
        '#define KLS_COMPACT_BTF_SOLVE_DIAGNOSTIC_EVENT(s) (++(s)->diagnostic_btf_solve)\n'
        '#define KLS_DIRECT_BTF_SOLVE_TEST_EVENT(s) (++(s)->diagnostic_btf_direct)\n')
    insert('kls.c','  ++s->private_caller_depth;','''
  const int diagnostic_owned=kls_small_btf_owner_matches(s);
  const uint64_t diagnostic_factor=s->diagnostic_btf_factor,diagnostic_solve=s->diagnostic_btf_solve,
    diagnostic_direct=s->diagnostic_btf_direct;
''')
    insert('kls.c','  status = private_caller_finish(s,0,nrhs,saved,candidate,0,status);','''
  if(status==KLS_OK && diagnostic_owned && kls_small_btf_owner_matches(s) &&
     s->options.ordering==KLS_ORDERING_AUTO && s->stats.selected_ordering==KLS_ORDERING_METIS &&
     s->numeric && !s->numeric->Rs && s->common.scale<=0 &&
     s->diagnostic_btf_factor>diagnostic_factor && s->diagnostic_btf_solve>diagnostic_solve &&
     s->diagnostic_btf_direct>diagnostic_direct)++s->diagnostic_btf_accepted;
''')
    insert('kls.c','void kls_destroy(kls_solver *solver) {','''
  if(solver)fprintf(stderr,"SMALL_BTF_ACTIVATION n=%llu threads=%d requested=%d selected=%d retained=%d factor=%llu solve=%llu direct=%llu accepted=%llu scale=%d rs=%d\\n",
    (unsigned long long)solver->n,solver->options.threads,(int)solver->options.ordering,
    (int)solver->stats.selected_ordering,kls_small_btf_owner_matches(solver),
    (unsigned long long)solver->diagnostic_btf_factor,(unsigned long long)solver->diagnostic_btf_solve,
    (unsigned long long)solver->diagnostic_btf_direct,(unsigned long long)solver->diagnostic_btf_accepted,
    (int)solver->common.scale,solver->numeric && solver->numeric->Rs!=NULL);
''')
    for name,adds in insertions.items():
        text=(diagnostic/name).read_text()
        for added in reversed(adds):assert text.count(added)==1;text=text.replace(added,'',1)
        assert text==(core/'source/src'/name).read_text()
    save(out/'activation-insertions.json',insertions)
    flags=[]
    for line in (core/'native/CMakeFiles/kls.dir/flags.make').read_text().splitlines():
        if line.startswith(('C_DEFINES =','C_INCLUDES =','C_FLAGS =')):flags+=shlex.split(line.split(' = ',1)[1])
    activation=out/'activation';activation.mkdir();obj=activation/'kls.c.o'
    run('activation-compile',prefix+['gcc',*flags,'-c',diagnostic/'kls.c','-o',obj])
    shutil.copyfile(core/'native/libkls.a',activation/'libkls.a')
    run('activation-archive',prefix+['ar','r',activation/'libkls.a',obj])
    run('activation-index',prefix+['ranlib',activation/'libkls.a'])
    ordinary=out/'candidate-link';build=read(ordinary/'build.json')
    command=build['commands'][-1]['command']
    assert command.count(str(core/'native/libkls.a'))==1 and command.count(str(ordinary/'Xyce'))==1
    command=[str(activation/'libkls.a') if x==str(core/'native/libkls.a') else str(activation/'Xyce') if x==str(ordinary/'Xyce') else x for x in command]
    run('activation-link',command)
    activation_pins={str(p):sha(p) for p in diagnostic.rglob('*') if p.is_file()}
    activation_pins.update({str(p):sha(p) for p in (activation/'Xyce',activation/'libkls.a',out/'activation-insertions.json')})
    save(activation/'build.json',dict(instrumented=True,input_sha256={**pins,**activation_pins}))
    pins.update(activation_pins);save(out/'inputs.json',dict(diagnostic_only=True,input_sha256=pins))
    dest=campaign('activation-api',activation/'Xyce',activation/'build.json')
    rows=[json.loads(s) for s in (dest/'observations.jsonl').read_text().splitlines()];events=[]
    for row in rows:
        if row['side']!='candidate':continue
        log=(Path(row['directory'])/'stdout.log').read_text()
        entries=[dict(zip(('n','threads','requested','selected','retained','factor','solve','direct','accepted','scale','rs'),map(int,m)))
                 for m in re.findall(r'SMALL_BTF_ACTIVATION n=(\d+) threads=(\d+) requested=(\d+) selected=(\d+) retained=(\d+) factor=(\d+) solve=(\d+) direct=(\d+) accepted=(\d+) scale=(-?\d+) rs=(\d+)',log)]
        if row['case']=='MOS13_IC/invert50_mos1.cir':
            assert any(e['requested']==0 and e['selected']==4 and e['retained']==1 and e['accepted']>0 and e['scale']<=0 and e['rs']==0 and
                       all(e[k]>=e['accepted'] for k in ('factor','solve','direct')) for e in entries)
        if row['case']=='CircuitSim90/MOS2/mux8.cir':assert not any(e['retained'] for e in entries)
        events.append(dict(domain=row['domain'],case=row['case'],threads=row['threads'],owners=entries))
    assert len(events)==28
    save(out/'activation-events.json',events);verify()
    save(out/'completed.json',dict(passed=True,api70_qualified=True,activation70_qualified=True,
        retained_execution_all_four_inverter_cells=True,independent_exact_agreement=True,
        performance_qualified=False,source_pins_unchanged=True,
        ordinary_candidate=str(out/'candidate-link/Xyce'),ordinary_baseline=str(out/'baseline-link/Xyce'),
        input_manifest_sha256=sha(out/'inputs.json'),commands_sha256=sha(out/'commands.json')))
    print('API and retained AUTO execution qualified; performance remains separate.',flush=True)


if __name__=='__main__':main()
