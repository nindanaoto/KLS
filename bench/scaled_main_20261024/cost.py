from common import *
import statistics,math
# Frozen protocol: four ABBA blocks for healthy paths, two for failed-frame
# recovery. CPUs are reserved from this task's other campaigns. Every sample
# and failure is retained; only certified results count as recovered answers.
runner=Runner('cost','16-23');W=Path('/tmp/kls-main-recovery-cost-20261024');cases=[]
for n in [128,4096,32768]:
 for threads in [1,8]:cases.append(dict(name=f'healthy-{n}-t{threads}',kind='healthy',args=['healthy',str(n),str(threads)],blocks=4))
T=R/'build/subtree-investigation-20261024';frames=[T/'followup/invert50_mos1/failed-frame.bin',T/'nd/invertbjt/failed-frame.bin',R/'build/scaled-recovery-20261024/remaining/failed-frame.bin']
for name,p in zip(['invert50','invertbjt','invertbjt-later'],frames):cases.append(dict(name=name,kind='frame',args=['frame',str(p)],frame_sha256=sha(p),blocks=2))
protocol=dict(cases=cases,healthy_gate='Geometric mean candidate/control <=1.02 and each case <=1.05, using per-role median process means over four ABBA blocks. Every healthy call must succeed. Failed-frame control errors are not speed references.',isolation='No other task-owned campaigns active; same image, compiler and CPUs, alternating order.',recovery_gate='All three exact failed frames and their repeated/new-epoch calls certify. Record cold latency and peak RSS; never compare an error against a valid solution as a speedup.')
(O/'cost-protocol.json').write_text(json.dumps(protocol,indent=2));rows=[]
for case in cases:
 for block in range(case['blocks']):
  for pos,role in enumerate(['control','candidate','candidate','control']):
   terminal=runner.run(case['name']+'-'+str(block)+'-'+str(pos)+'-'+role,[W/role,*case['args']],timeout=120,memory_gb=8)
   metrics=json.loads(next(x for x in Path(terminal['log']).read_text().splitlines() if x.startswith('{')))
   rows.append(dict(case=case['name'],kind=case['kind'],block=block,position=pos,side=role,metrics=metrics,terminal=terminal));(O/'cost-observations.json').write_text(json.dumps(rows,indent=2))
healthy=[];recovery=[]
for case in cases:
 selected=[r for r in rows if r['case']==case['name']];by={role:[r['metrics'] for r in selected if r['side']==role] for role in ['control','candidate']}
 if case['kind']=='healthy':
  assert all(m['successes']==m['calls'] for metrics in by.values() for m in metrics)
  med={role:statistics.median(m['mean_repeated_seconds'] for m in metrics) for role,metrics in by.items()};healthy.append(dict(case=case['name'],median_mean_seconds=med,ratio=med['candidate']/med['control']))
 else:
  assert all(m['successes']==m['calls'] for m in by['candidate'])
  recovery.append(dict(case=case['name'],candidate_cold_seconds=statistics.median(m['cold_seconds'] for m in by['candidate']),candidate_repeated_seconds=statistics.median(m['median_repeated_seconds'] for m in by['candidate']),candidate_new_epoch_seconds=statistics.median(m['new_epoch_seconds'] for m in by['candidate']),peak_rss_kb={role:statistics.median(m['peak_final_kb'] for m in metrics) for role,metrics in by.items()},control_initial_statuses=[m['initial_status'] for m in by['control']]))
gmean=math.exp(statistics.mean(math.log(r['ratio']) for r in healthy));gate=gmean<=1.02 and all(r['ratio']<=1.05 for r in healthy)
(O/'cost-completed.json').write_text(json.dumps(dict(healthy_gate_passed=gate,geomean_ratio=gmean,healthy=healthy,recovery=recovery,observations=rows,processes=runner.rows),indent=2));print('COST COMPLETE',gmean,gate,json.dumps(recovery),flush=True)
