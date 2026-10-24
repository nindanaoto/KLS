from common import *
import sys,shlex,statistics,math
version,campaign=sys.argv[1:];assert campaign in ['initial','confirmation'];order=['control','candidate','candidate','control'] if campaign=='initial' else ['candidate','control','control','candidate'];W=Path('/tmp/kls-paper-healthy-'+version+'-'+campaign+'-20261024');W.mkdir();runner=Runner('healthy-'+version+'-'+campaign,'16-23');rows=[]
for role in ['control','candidate']:
 B=Path('/tmp/kls-main-recovery-candidate-native-20261024') if role=='control' else build('recovery',version)
 flags=[];text=(B/'CMakeFiles/kls.dir/flags.make').read_text()
 for key in ['C_DEFINES','C_INCLUDES','C_FLAGS']:flags+=shlex.split(next(x for x in text.splitlines() if x.startswith(key+' = ')).split(' = ',1)[1])
 obj=W/(role+'.o');runner.run(role+'-compile',['gcc',*flags,'-c',WORK/'KLS-scaled-main/bench/scaled_main_20261024/cost.c','-o',obj])
 cmd=shlex.split((B/'CMakeFiles/kls_smoke.dir/link.txt').read_text());cmd=[str(obj) if x=='CMakeFiles/kls_smoke.dir/tests/kls_smoke.c.o' else x for x in cmd if not x.startswith('-Wl,--dependency-file=')];cmd[cmd.index('-o')+1]=str(W/role);runner.run(role+'-link',cmd,cwd=B)
(O/(runner.name+'-protocol.json')).write_text(json.dumps(dict(campaign=campaign,order=order,blocks=4,sizes=[128,4096,32768],threads=[1,8],gate='Every result valid; each setting median candidate/control <=1.02; both campaigns required'),indent=2))
for n in [128,4096,32768]:
 for threads in [1,8]:
  for block in range(4):
   for pos,role in enumerate(order):
    terminal=runner.run(str(n)+'-t'+str(threads)+'-'+str(block)+'-'+str(pos)+'-'+role,[W/role,'healthy',str(n),str(threads)],timeout=120,memory_gb=8)
    metrics=json.loads(next(x for x in Path(terminal['log']).read_text().splitlines() if x.startswith('{')));rows.append(dict(n=n,threads=threads,side=role,block=block,position=pos,metrics=metrics,terminal=terminal));(O/(runner.name+'-observations.json')).write_text(json.dumps(rows,indent=2))
summary=[]
for n in [128,4096,32768]:
 for threads in [1,8]:
  chosen=[x for x in rows if x['n']==n and x['threads']==threads];med={role:statistics.median(x['metrics']['mean_repeated_seconds'] for x in chosen if x['side']==role) for role in ['control','candidate']};valid=all(x['metrics']['calls']==x['metrics']['successes'] for x in chosen);summary.append(dict(n=n,threads=threads,valid=valid,median=med,ratio=med['candidate']/med['control']))
(O/(runner.name+'-completed.json')).write_text(json.dumps(dict(passed=all(x['valid'] and x['ratio']<=1.02 for x in summary),summary=summary,processes=runner.rows),indent=2));print('HEALTHY',json.dumps(summary),flush=True)
