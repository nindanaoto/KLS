from common import *
import sys,shlex,statistics,math
version=sys.argv[1];campaign=sys.argv[2];assert campaign in ['initial','confirmation'];order=['control','candidate','candidate','control'] if campaign=='initial' else ['candidate','control','control','candidate'];W=Path('/tmp/kls-paper-recovery-epochs-'+version+'-'+campaign+'-20261024');W.mkdir();runner=Runner('recovery-epochs-'+version+'-'+campaign,'8-15');rows=[]
for role in ['control','candidate']:
 B=Path('/tmp/kls-main-recovery-candidate-native-20261024') if role=='control' else build('recovery',version)
 S=WORK/'KLS-scaled-main' if role=='control' else source('recovery',version)
 f=(B/'CMakeFiles/kls.dir/flags.make').read_text();flags=[]
 for key in ['C_DEFINES','C_INCLUDES','C_FLAGS']:flags+=shlex.split(next(x for x in f.splitlines() if x.startswith(key+' = ')).split(' = ',1)[1])
 if role=='candidate':flags+=['-DPAPER_RECOVERY_REUSE']
 runner.run(role+'-compile',['gcc',*flags,'-DKLS_CORE_SOURCE="'+str(S/'src/kls.c')+'"','-c',Path(__file__).with_name('recovery_epochs.c'),'-o',W/(role+'.o')])
 cmd=shlex.split((B/'CMakeFiles/kls_scaled_recovery.dir/link.txt').read_text());cmd=[str(W/(role+'.o')) if x=='CMakeFiles/kls_scaled_recovery.dir/tests/kls_scaled_recovery.c.o' else x for x in cmd if not x.startswith('-Wl,--dependency-file=')];cmd[cmd.index('-o')+1]=str(W/role);runner.run(role+'-link',cmd,cwd=B)
T=R/'build/subtree-investigation-20261024';frames=[('invert50',T/'followup/invert50_mos1/failed-frame.bin'),('invertbjt',T/'nd/invertbjt/failed-frame.bin'),('invertbjt-later',R/'build/scaled-recovery-20261024/remaining/failed-frame.bin')]
(O/(runner.name+'-protocol.json')).write_text(json.dumps(dict(value_families=['row powers of two','column powers of two','entrywise 1e-4 deterministic, analytic RHS'],updates=8,order=order,blocks=2,campaign=campaign,primary_target='row-rescaled captured operators, both implementations certified; >=2% faster in each frame',maximum_healthy_ratio=1.02,frames={name:sha(p) for name,p in frames}),indent=2))
for name,p in frames:
 for family in [0,1,2]:
  for block in range(2):
   for pos,role in enumerate(order):
    terminal=runner.run(name+'-f'+str(family)+'-'+str(block)+'-'+str(pos)+'-'+role,[W/role,p,str(family)],required=False,timeout=180,memory_gb=8)
    text=Path(terminal['log']).read_text();metrics=next((json.loads(x) for x in text.splitlines() if x.startswith('{')),None)
    rows.append(dict(case=name,family=family,block=block,position=pos,side=role,metrics=metrics,terminal=terminal));(O/(runner.name+'-observations.json')).write_text(json.dumps(rows,indent=2))
summary=[]
for name,p in frames:
 for family in [0,1,2]:
  selected=[x for x in rows if x['case']==name and x['family']==family]
  valid={role:all(x['metrics'] and x['metrics']['valid']==x['metrics']['epochs'] for x in selected if x['side']==role) for role in ['control','candidate']}
  med={role:statistics.median(x['metrics']['mean_update_seconds'] for x in selected if x['side']==role and x['metrics']) for role in ['control','candidate']}
  summary.append(dict(case=name,family=family,valid=valid,median_update_seconds=med,ratio=med['candidate']/med['control'],minimum_reuse_accepts=min(x['metrics']['reuse_accepts'] for x in selected if x['side']=='candidate' and x['metrics'])))
(O/(runner.name+'-completed.json')).write_text(json.dumps(dict(summary=summary,observations=rows,processes=runner.rows,diagnostic_only=False),indent=2));print('SUMMARY',json.dumps(summary),flush=True)
