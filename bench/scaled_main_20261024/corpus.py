from common import *
import sys
manifest=Path(__file__).with_name('manifest.json');d=json.loads(manifest.read_text());threads=int(sys.argv[1]);assert threads in [1,8];runner=Runner('corpus-t'+str(threads),'0-7' if threads==1 else '8-15');results=[]
for i,m in enumerate(d['rows']):
 p=R/m['path'];assert sha(p)==m['sha256']
 for role in (['control','candidate'] if i%2==0 else ['candidate','control']):
  label=str(i)+'-'+role;row=runner.run(label,[build(role)/'kls_bench',p,'--lifecycle-systems',str(d['systems']),'--threads',str(threads),'--json'],timeout=d['timeout_seconds'],required=False,memory_gb=d['address_space_gb'])
  metrics=None
  for line in Path(row['log']).read_text().splitlines():
   if line.startswith('{'):
    try:metrics=json.loads(line)
    except json.JSONDecodeError:pass
  results.append(dict(matrix=m['path'],threads=threads,side=role,terminal=row,metrics=metrics,success=row['exit_code']==0 and metrics is not None and metrics.get('status')==0))
  (O/(runner.name+'-observations.json')).write_text(json.dumps(results,indent=2));print('matrix',i+1,'of',len(d['rows']),role,results[-1]['success'],flush=True)
regressions=[];recoveries=[]
for m in d['rows']:
 a=next(r for r in results if r['matrix']==m['path'] and r['side']=='control');b=next(r for r in results if r['matrix']==m['path'] and r['side']=='candidate')
 if a['success'] and not b['success']:regressions.append(m['path'])
 if not a['success'] and b['success']:recoveries.append(m['path'])
(O/(runner.name+'-completed.json')).write_text(json.dumps(dict(manifest_sha256=sha(manifest),all_preserved=not regressions,new_regressions=regressions,new_successes=recoveries,observations=results,processes=runner.rows),indent=2));print('CORPUS',len(d['rows']),'regressions',regressions,'recoveries',recoveries,flush=True)
