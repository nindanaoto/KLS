from common import *
import shlex,shutil,re,math
A=R/'build/subtree-execution-20261024/applications-v2';P=json.loads((R/'bench/xyce_testing.json').read_text());W=Path('/tmp/kls-main-recovery-xyce-20261024');W.mkdir();runner=Runner('xyce','0-7');S=source('control');flags=(R/'build/xyce-testing-bump-20261023/build-provenance/src/CMakeFiles/XyceLib.dir/flags.make').read_text();args=[]
for key in ['CXX_DEFINES','CXX_INCLUDES','CXX_FLAGS']:args+=shlex.split(next(x for x in flags.splitlines() if x.startswith(key+' = ')).split(' = ',1)[1])
args=[x.replace('/tmp/kls-three-control-v1-source-20261022/include',str(S/'include')) for x in args];assert '-I'+str(S/'include') in args
adapter=WORK/'Xyce/src/LinearAlgebraServicesPKG/N_LAS_KLSSolver.C';obj=W/'N_LAS_KLSSolver.C.o';runner.run('adapter-compile',['g++',*args,'-c',adapter,'-o',obj]);archive=W/'libxyce.a';shutil.copy2(R/P['xyce_archive'],archive);runner.run('archive',['ar','r',archive,obj])
template=next(x['command'] for x in json.loads((A/'build.json').read_text())['commands'] if x['side']=='candidate');old='/tmp/kls-subtree-candidate-v2-full-20261024-native'
for role in ['control','candidate']:
 cmd=[str(archive) if x=='/tmp/kls-subtree-applications-v2-20261024/libxyce.a' else x.replace(old,str(build(role))) for x in template if not x.startswith('-Wl,-Map=')];cmd[cmd.index('-o')+1]=str(W/role);runner.run(role+'-link',cmd)
spec=importlib.util.spec_from_file_location('whole',R/'build/klu-short-producer-v2-spice-pair-20261020/run_campaign.py');whole=importlib.util.module_from_spec(spec);spec.loader.exec_module(whole)
oldrows=json.loads((A/'auto.json').read_text())['rows'];cases=sorted(set(x['case'] for x in oldrows));results=[]
for case in cases:
 ref=None
 for role,threads in [('klu',1),('control',1),('candidate',1),('control',8),('candidate',8)]:
  dest=O/'applications'/Path(case).stem/(role+'-t'+str(threads));dest.mkdir(parents=True)
  prev=next(x for x in oldrows if x['case']==case and x['domain']=='llc96' and x['threads']==threads and x['side']==('klu' if role=='klu' else 'candidate'));shutil.copy2(Path(prev['directory'])/'input.cir',dest/'input.cir')
  terminal=runner.run(Path(case).stem+'-'+role+'-t'+str(threads),[W/('control' if role=='klu' else role),'input.cir'],cwd=dest,required=False,timeout=300)
  text=Path(terminal['log']).read_text();success=re.findall(r'Number Successful Steps Taken:\s*(\d+)',text);failed=re.findall(r'Number Failed Steps Attempted:\s*(\d+)',text)
  waveform={} if role=='klu' else whole.verified_waveform(Path(ref['directory'])/'input.cir.prn',dest/'input.cir.prn',dest)
  good=terminal['exit_code']==0 and 'End of Xyce(TM) Simulation' in text and not re.search('linear solve failed',text,re.I) and bool(success) and bool(failed) and (dest/'input.cir.prn').exists()
  if role!='klu':good=good and success==ref['successful_steps'] and failed==ref['failed_steps'] and waveform['status']=='pass'
  row=dict(case=case,side=role,threads=threads,directory=str(dest),qualified=bool(good),successful_steps=success,failed_steps=failed,waveform=waveform,terminal=terminal);results.append(row);(O/'xyce-observations.json').write_text(json.dumps(results,indent=2));print(case,role,threads,good,flush=True)
  if role=='klu':ref=row
assert len(results)==35
(O/'xyce-completed.json').write_text(json.dumps(dict(all_qualified=all(x['qualified'] for x in results),acceptance_timing=False,observations=results,processes=runner.rows,binaries={role:dict(path=str(W/role),sha256=sha(W/role)) for role in ['control','candidate']},xyce_revision=P['xyce_revision']),indent=2))
