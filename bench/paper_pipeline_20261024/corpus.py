from common import *
import sys
role,version,thread=sys.argv[1:];threads=int(thread);assert role in ['subtree','recovery'] and threads in [1,8]
manifest=WORK/'KLS-scaled-main/bench/scaled_main_20261024/manifest.json';d=json.loads(manifest.read_text());runner=Runner(role+'-'+version+'-corpus-t'+thread,'0-7' if threads==1 else '8-15');results=[]
for i,m in enumerate(d['rows']):
 p=R/m['path'];assert sha(p)==m['sha256']
 for side in (['control','candidate'] if i%2==0 else ['candidate','control']):
  B=Path('/tmp/kls-main-recovery-candidate-native-20261024') if side=='control' else build(role,version)
  terminal=runner.run(str(i)+'-'+side,[B/'kls_bench',p,'--lifecycle-systems',str(d['systems']),'--threads',thread,'--json'],timeout=d['timeout_seconds'],required=False,memory_gb=d['address_space_gb'])
  metrics=next((json.loads(x) for x in Path(terminal['log']).read_text().splitlines() if x.startswith('{')),None)
  results.append(dict(matrix=m['path'],threads=threads,side=side,terminal=terminal,metrics=metrics,success=terminal['exit_code']==0 and metrics is not None and metrics.get('status')==0));(O/(runner.name+'-observations.json')).write_text(json.dumps(results,indent=2))
regressions=[]
for m in d['rows']:
 a=next(r for r in results if r['matrix']==m['path'] and r['side']=='control');b=next(r for r in results if r['matrix']==m['path'] and r['side']=='candidate')
 if a['success'] and not b['success']:regressions.append(m['path'])
(O/(runner.name+'-completed.json')).write_text(json.dumps(dict(all_preserved=not regressions,new_regressions=regressions,observations=results,processes=runner.rows,manifest_sha256=sha(manifest)),indent=2));print('CORPUS',role,threads,regressions,flush=True)
