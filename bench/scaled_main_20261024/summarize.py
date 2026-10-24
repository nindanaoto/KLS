"""Seal the focused main qualification before fast-forwarding main."""
from common import *
import re
C=source('candidate');B=source('control');target=C/'docs/scaled_recovery_main_qualification_20261024.json';assert not target.exists()
def git(path,*args):return subprocess.check_output(['git','-C',str(path),*args],text=True).strip()
base=git(B,'rev-parse','HEAD');assert base=='e8077c22ad729cf4aaf80069cae1537cc8e4a707'
commit=git(C,'rev-parse','HEAD');assert commit.startswith('a33ec89')
assert git(R,'rev-parse','main')==base
assert not git(B,'status','--porcelain')
assert not git(C,'diff','--name-only','HEAD','--','src','include','tests','CMakeLists.txt')
expected={'CMakeLists.txt','src/kls.c','src/kls_componentwise_recovery.inc','tests/kls_componentwise_caller.c','tests/kls_scaled_recovery.c','tests/kls_scaled_recovery_fixtures.h'}
assert set(git(C,'diff','--name-only',base,commit).splitlines())==expected
assert sha(C/'include/kls/kls.h')==sha(B/'include/kls/kls.h')
assert git(WORK/'Xyce','rev-parse','HEAD')=='1be7670d4872c2334beb6d00e8ec59885176eca0'
checks={};processes={};tests={}
def retain(p):checks[str(Path(p).resolve())]=sha(p)
def audit(row):
 assert row['closure']['all_session_processes_absent'] and not Path('/proc',str(row['pid'])).exists(),row
 processes[row['pid']]=dict(exit_code=row['exit_code'],timed_out=row.get('timed_out',False))
for role,mode,count in [('control','native',43),('candidate','native',44),('candidate','portable',44),('candidate','asan',3)]:
 p=O/(role+'-'+mode+'-completed.json');d=json.loads(p.read_text());assert d['passed'] and sha(d['library'])==d['library_sha256'];retain(d['library'])
 assert f'100% tests passed, 0 tests failed out of {count}' in (O/(role+'-'+mode+'-tests.log')).read_text()
 for row in d['processes']:assert row['exit_code']==0;audit(row)
 tests[role+'-'+mode]=count
 for path in [build(role,mode)/'CMakeCache.txt',build(role,mode)/'CMakeFiles/kls.dir/flags.make']:retain(path)
corpus=[]
for t in [1,8]:
 d=json.loads((O/(f'corpus-t{t}-completed.json')).read_text());assert d['all_preserved'] and len(d['observations'])==138
 assert d['manifest_sha256']==sha(Path(__file__).with_name('manifest.json'))
 pairs=[]
 for item in json.loads(Path(__file__).with_name('manifest.json').read_text())['rows']:
  p=R/item['path'];assert sha(p)==item['sha256'];retain(p)
  a=next(x for x in d['observations'] if x['matrix']==item['path'] and x['side']=='control');b=next(x for x in d['observations'] if x['matrix']==item['path'] and x['side']=='candidate')
  pairs.append(dict(matrix=item['path'],control_success=a['success'],candidate_success=b['success'],control_exit=a['terminal']['exit_code'],candidate_exit=b['terminal']['exit_code']))
 for row in d['processes']:audit(row)
 corpus.append(dict(threads=t,matrices=69,control_successes=sum(x['control_success'] for x in pairs),candidate_successes=sum(x['candidate_success'] for x in pairs),new_regressions=d['new_regressions'],new_successes=d['new_successes'],pairs=pairs))
xyce=json.loads((O/'xyce-completed.json').read_text());assert xyce['all_qualified'] and len(xyce['observations'])==35
for row in xyce['processes']:assert row['exit_code']==0;audit(row)
for role,entry in xyce['binaries'].items():assert sha(entry['path'])==entry['sha256'];retain(entry['path'])
cost=json.loads((O/'cost-completed.json').read_text());assert cost['healthy_gate_passed'] and len(cost['observations'])==120
for row in cost['processes']:assert row['exit_code']==0;audit(row)
cb=json.loads((O/'cost-build-completed.json').read_text());assert cb['passed']
for row in cb['processes']:assert row['exit_code']==0;audit(row)
for entry in cb['binaries'].values():assert sha(entry['path'])==entry['sha256'];retain(entry['path'])
# Keep the failed affinity launches in the audit, never in the numerical gate.
for p in O.rglob('*processes.json'):
 for row in json.loads(p.read_text()):audit(row)
for root in [C/'src',C/'include',C/'tests',C/'bench/scaled_main_20261024',B/'src',B/'include',O]:
 for p in root.rglob('*'):
  if p.is_file() and '__pycache__' not in p.parts:retain(p)
for p in [C/'CMakeLists.txt',B/'CMakeLists.txt',IMAGE]:retain(p)
protocol=json.loads((O/'cost-protocol.json').read_text())
for case in protocol['cases']:
 if case['kind']=='frame':p=Path(case['args'][1]);assert sha(p)==case['frame_sha256'];retain(p)
# Recipes were executed from the research checkout and mirrored byte-for-byte
# here so the merged main tree owns the exact source of the qualification.
for p in (C/'bench/scaled_main_20261024').iterdir():
 if p.is_file():assert sha(p)==sha(R/'bench/scaled_main_20261024'/p.name)
report=dict(date='2026-10-24',status='qualified_for_main_as_robustness_fix',base_commit=base,candidate_source_commit=commit,scope=sorted(expected),public_header_unchanged=True,hot_field_offsets_preserved_by_tail_append=True,production_performance_claim=False,
 tests=tests,corpus=corpus,xyce=dict(all_qualified=True,runs=35,klu_references=7,control_settings=14,candidate_settings=14,threads=[1,8],fresh_main_header_adapter=True,fresh_main_dependency_links=True,waveforms_and_step_counts_match=True,acceptance_timings=False,revision=xyce['xyce_revision']),
 cost=dict(protocol=protocol,healthy_gate_passed=True,geometric_mean_ratio=cost['geomean_ratio'],healthy=cost['healthy'],recovery=cost['recovery'],timing_scope='Isolated task-owned campaigns; same Apptainer image, GCC flags, CPUs and ABBA order. Healthy samples mix normal, transpose and combined refactor/solve calls. Recovery controls fail, so their timings are not valid-solution speed references.'),
 behavior=dict(trigger='Observed DEFAULT componentwise failure after bounded correction, AUTO scaling requested and an unscaled incumbent. Kernel failures are not certified.',attempts='At most MAX then SUM scaling, with one isolated proposal live at a time and recursive recovery disabled.',acceptance='All RHSs independently pass the unchanged componentwise check against the current original operator before any publication or factor retention.',reuse='The incumbent is tried first on every RHS. A recovery factor is reused only after actual failure, within one numeric epoch. Public numeric updates and reanalysis discard it.',preserved='Explicit scaling requests, original pivot options, accuracy policy, atomic multi-RHS output and failure rollback. No ordering or parallel-executor admission changes.'),
 limits=['The unchanged bips07_1998 factorization failure remains at both thread counts; all failures were retained.','Both official SuiteSparse download routes returned HTTP 403 for all 24 predefined generalization-manifest entries. The numerical corpus is the frozen 69-matrix available set, not an unseen holdout.','Initial invalid CPU-affinity launches never reached the solver; archived under invalid-affinity folders and excluded from numerical results.','Recovery may retain an additional factor and pays cold factorization again after a numeric update. Reported RSS and latency are observations on three captured operators, not a global memory or wall-time bound.','This patch establishes robustness and measured healthy-path nonregression, not that KLS is fastest. Full forced-ND BJT previously remained unqualified; no forced-ND performance or production admission is claimed here.'],
 process_audit=dict(unique_children=len(processes),all_sessions_closed=True,nonzero_exits=sum(x['exit_code']!=0 for x in processes.values()),raw_receipts_location=str(O)),checked_sha256=checks)
for p,h in checks.items():assert sha(p)==h,p
target.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(report=str(target),source=commit,tests=tests,corpus=[dict(threads=x['threads'],successes=x['candidate_successes'],new_regressions=x['new_regressions']) for x in corpus],healthy_ratio=cost['geomean_ratio'],xyce_runs=35),indent=2),flush=True)
