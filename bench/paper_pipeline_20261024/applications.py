from common import *
import sys,shlex,shutil,re
version=sys.argv[1];resume='--resume' in sys.argv[2:];W=Path('/tmp/kls-paper-recovery-xyce-'+version+'-20261024');W.mkdir(exist_ok=resume);runner=Runner('recovery-'+version+'-xyce','0-7')
A=R/'build/subtree-execution-20261024/applications-v2';template=next(x['command'] for x in json.loads((A/'build.json').read_text())['commands'] if x['side']=='candidate');old='/tmp/kls-subtree-candidate-v2-full-20261024-native';archive=Path('/tmp/kls-main-recovery-xyce-20261024/libxyce.a')
cmd=[str(archive) if x=='/tmp/kls-subtree-applications-v2-20261024/libxyce.a' else x.replace(old,str(build('recovery',version))) for x in template if not x.startswith('-Wl,-Map=')];cmd[cmd.index('-o')+1]=str(W/'candidate');
if not resume:runner.run('candidate-link',cmd)
else:runner.rows=json.loads((O/(runner.name+'-processes.json')).read_text());assert (W/'candidate').exists()
spec=importlib.util.spec_from_file_location('whole',R/'build/klu-short-producer-v2-spice-pair-20261020/run_campaign.py');whole=importlib.util.module_from_spec(spec);spec.loader.exec_module(whole)
oldrows=json.loads((A/'auto.json').read_text())['rows'];cases=sorted(set(x['case'] for x in oldrows));results=json.loads((O/(runner.name+'-observations.json')).read_text()) if resume else []
if resume:assert len(results)==14 and all(x['domain']=='llc96' for x in results)
for domain,cpus in ([('llc32','8-15')] if resume else [('llc96','0-7'),('llc32','8-15')]):
 runner.cpus=cpus
 for case in cases:
  for threads in [1,8]:
   dest=O/('applications-'+version)/domain/Path(case).stem/('candidate-t'+str(threads));dest.mkdir(parents=True,exist_ok=resume)
   prev=next(x for x in oldrows if x['case']==case and x['domain']==domain and x['threads']==threads and x['side']=='candidate');shutil.copy2(Path(prev['directory'])/'input.cir',dest/'input.cir')
   ref=R/'build/scaled-main-20261024/applications'/Path(case).stem/'klu-t1';reftext=(R/'build/scaled-main-20261024'/('xyce-'+Path(case).stem+'-klu-t1.log')).read_text()
   terminal=runner.run(domain+'-'+Path(case).stem+'-t'+str(threads),[W/'candidate','input.cir'],cwd=dest,required=False,timeout=300)
   text=Path(terminal['log']).read_text();success=re.findall(r'Number Successful Steps Taken:\s*(\d+)',text);failed=re.findall(r'Number Failed Steps Attempted:\s*(\d+)',text)
   waveform=whole.verified_waveform(ref/'input.cir.prn',dest/'input.cir.prn',dest) if (dest/'input.cir.prn').exists() else {'status':'missing'}
   good=terminal['exit_code']==0 and 'End of Xyce(TM) Simulation' in text and not re.search('linear solve failed',text,re.I) and bool(success) and bool(failed) and success==re.findall(r'Number Successful Steps Taken:\s*(\d+)',reftext) and failed==re.findall(r'Number Failed Steps Attempted:\s*(\d+)',reftext) and waveform['status']=='pass'
   row=dict(case=case,domain=domain,threads=threads,qualified=bool(good),waveform=waveform,successful_steps=success,failed_steps=failed,terminal=terminal);results.append(row);(O/(runner.name+'-observations.json')).write_text(json.dumps(results,indent=2));print(case,domain,threads,good,flush=True)
assert len(results)==28
(O/(runner.name+'-completed.json')).write_text(json.dumps(dict(all_qualified=all(x['qualified'] for x in results),acceptance_timing=False,observations=results,processes=runner.rows,binary_sha256=sha(W/'candidate'),archive_sha256=sha(archive)),indent=2))
