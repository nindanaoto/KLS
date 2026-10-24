from common import *
import re,sys,statistics,struct

def read(name):return json.loads((O/(name+'.json')).read_text())
def receipt(name):
 d=read(name+'-completed');log=O/(name+'-tests.log');counts=re.findall(r'0 tests failed out of (\d+)',log.read_text()) if log.exists() else []
 return dict(passed=d.get('passed'),tests=int(counts[-1]) if counts else len(d.get('scope',[])),receipt=name+'-completed.json',sha256=sha(O/(name+'-completed.json')))
base='7d2c909c6ab77a7201bad8559b024f4e97625840'
frames=[]
for name,p in [('invert50',R/'build/subtree-investigation-20261024/followup/invert50_mos1/failed-frame.bin'),('invertbjt',R/'build/subtree-investigation-20261024/nd/invertbjt/failed-frame.bin'),('invertbjt-later',R/'build/scaled-recovery-20261024/remaining/failed-frame.bin')]:
 with p.open('rb') as f:header=f.read(128)
 fields=['struct_size','threads','ordering','orientation','use_btf','scale','abi_version','pivot_tolerance','memory_growth','halt_if_singular','fast_factor','static_pivoting','backend','expected_refactorizations','expected_solves','record_tiny_solve_timing','tuning_profile_path']
 options=dict(zip(fields,struct.unpack_from('<Q6i2d4i2qi4xQ',header,32)));assert options['struct_size']==96 and options['tuning_profile_path']==0
 magic,n,z,csr=struct.unpack_from('<4Q',header);assert magic==0x4b4c534449414731 and csr==1
 frames.append(dict(case=name,path=str(p),sha256=sha(p),n=n,nnz=z,options=options))
sources={}
for role,v in [('recovery','v6'),('subtree','compact-v11')]:
 d=read(role+'-'+v+'-source');w=WORK/('KLS-paper-'+role);bad=[name for name,h in d['sha256'].items() if sha(w/name)!=h];assert not bad,bad
 sources[role]=dict(tested_source_commit=d['head'],worktree=str(w),branch=subprocess.check_output(['git','-C',str(w),'branch','--show-current'],text=True).strip(),snapshot_files=len(d['sha256']),all_snapshot_inputs_match=True,snapshot_manifest_sha256=sha(O/(role+'-'+v+'-source.json')),commits=subprocess.check_output(['git','-C',str(w),'log','--format=%h %s',base+'..HEAD'],text=True).splitlines())
qualification={role:{mode:receipt(role+'-'+v+'-'+mode+('-full' if mode in ['native','portable'] else '')) for mode in ['native','portable','asan','tsan']} for role,v in [('recovery','v6'),('subtree','compact-v11')]}
for role,v in [('recovery','v6'),('subtree','compact-v11')]:
 name=role+'-'+v+'-asan-allocator-completed';d=read(name);qualification[role]['allocator_asan']=dict(passed=d['passed'],tests=1,receipt=name+'.json',sha256=sha(O/(name+'.json')))
corpus={}
for threads in [1,8]:
 d=read('recovery-v6-corpus-t'+str(threads)+'-completed');rows=d['observations'];corpus[str(threads)]=dict(count_per_side=69,passes={side:sum(x['success'] for x in rows if x['side']==side) for side in ['control','candidate']},new_regressions=d['new_regressions'],failures=[dict(matrix=x['matrix'],side=x['side'],status=x['metrics'].get('status') if x['metrics'] else None) for x in rows if not x['success']],receipt_sha256=sha(O/('recovery-v6-corpus-t'+str(threads)+'-completed.json')))
apps=read('recovery-v6-xyce-completed');xyce=dict(settings=len(apps['observations']),passes=sum(x['qualified'] for x in apps['observations']),all_qualified=apps['all_qualified'],binary_sha256=apps['binary_sha256'],archive_sha256=apps['archive_sha256'],domains=sorted(set(x['domain'] for x in apps['observations'])),acceptance_timing=False,harness_interruption='After 14 valid settings, a metadata lookup used nonexistent label llc97. Corrected to stored llc32; resumed the 14 missing settings without repeating completed settings.')
epochs={};healthy={};primary=True;paired_no_regression=True
for campaign in ['initial','confirmation']:
 d=read('recovery-epochs-v6-'+campaign+'-completed');epochs[campaign]=[]
 for raw in d['summary']:
  x=dict(raw);x['speed_reference_valid']=all(x['valid'].values());x['qualified_ratio']=x['ratio'] if x['speed_reference_valid'] else None
  if not x['speed_reference_valid']:x['unqualified_elapsed_ratio']=x.pop('ratio')
  x['median_peak_rss_kb']={side:statistics.median(r['metrics']['peak_kb'] for r in d['observations'] if r['case']==x['case'] and r['family']==x['family'] and r['side']==side) for side in ['control','candidate']}
  epochs[campaign].append(x)
 for x in d['summary']:
  if x['family']==0:primary &= all(x['valid'].values()) and x['ratio']<=.98 and x['minimum_reuse_accepts']>0
  if all(x['valid'].values()):paired_no_regression &= x['ratio']<=1.02
  if x['valid']['control'] and not x['valid']['candidate']:paired_no_regression=False
 d=read('healthy-v6-'+campaign+'-completed');healthy[campaign]=dict(passed=d['passed'],summary=d['summary'])
healthy_gate=all(x['passed'] for x in healthy.values());correctness=all(x['passed'] for modes in qualification.values() for x in modes.values()) and all(not x['new_regressions'] for x in corpus.values()) and apps['all_qualified']
probe_version=sys.argv[1] if len(sys.argv)>1 else 'compact-v9';probes=[]
for matched in [0,1]:
 for threads in [2,8]:
  label=('matched' if matched else 'original')+'-t'+str(threads);log=O/('subtree-probe-'+probe_version+'-'+label+'.log');line=next(x for x in log.read_text().splitlines() if x.startswith('probe '));fields=dict(re.findall(r'(\w+)=([^ ]+)',line));probes.append(dict(fields=fields,log=log.name,sha256=sha(log)))
processes=[]
for p in sorted(O.glob('*-processes.json')):
 rows=json.loads(p.read_text());assert all(x['closure']['all_session_processes_absent'] for x in rows),p
 processes.append(dict(file=p.name,sha256=sha(p),stages=len(rows),nonzero_stages=[dict(stage=x['stage'],exit_code=x['exit_code'],timed_out=x['timed_out']) for x in rows if x['exit_code'] or x['timed_out']]))
papers=[]
for p in sorted((WORK/'refs').glob('*.pdf')):papers.append(dict(path=str(p.relative_to(WORK)),sha256=sha(p)))
report=dict(baseline=base,main_unchanged=subprocess.check_output(['git','-C',str(WORK/'KLS-scaled-main'),'rev-parse','HEAD'],text=True).strip()==base,sources=sources,policy=read('policy') if (O/'policy.json').exists() else json.loads(Path(__file__).with_name('policy.json').read_text()),qualification=qualification,recovery=dict(timing_scope="Captured METIS operators at one thread with deterministic numeric changes; not an end-to-end AUTO Xyce speed claim.",frames=frames,corpus=corpus,xyce=xyce,epoch_timing=epochs,healthy_timing=healthy,gates=dict(correctness=bool(correctness),primary_2_percent_both_orders=bool(primary),every_valid_paired_epoch_within_2_percent=bool(paired_no_regression),healthy_within_2_percent_both_orders=healthy_gate),promotion_recommended=bool(correctness and primary and paired_no_regression and healthy_gate)),subtree=dict(probes=probes,unchanged_heap_limit_bytes=21850128,matching='Owned diagnostic SPRAL weighted-matching row permutation before ND, outside the numeric ledger; no claim that full matching+factor proposal meets the cap.',decision='Stop before BTF-preserving AUTO admission: matched pipeline workspace cannot fit the unchanged numeric cap. Scheduler work gated off.',auto_admission_implemented=False,scheduler_implemented=False,production_promotion_recommended=False),sanitizer_development_failures_preserved=['recovery-v4-tsan: cache preparation race; fixed by joining builders','subtree-compact-v9-tsan: independent KLU analysis writes global AMD allocator hooks; fixed by call-owned allocators'],papers=papers,paper_mapping=dict(SubtreeLU='Matching before ND, private scratch ownership, mandatory private/pipeline first-factor barrier; fixed-pattern refactor scheduling deferred until active target.',CKTSO='Checked retained-pivot refactors motivate retained scaled proposal; every new operator/RHS is certified independently.',EScheduler='Dependency-aware scheduling remains gated by native activation and measured avoidable waits.',HYLU_and_Advancing='Kernel/scheduling and randomized outer-solve ideas remain deferred; no relaxed correctness policy.'),process_audit=processes,raw_artifact_root=str(O),image=dict(path=str(IMAGE),sha256=sha(IMAGE)),no_main_merge_or_push=True)
artifacts={p.name:sha(p) for p in sorted(O.iterdir()) if p.is_file() and p.name not in ['summary.json','artifact-manifest.json'] and not p.name.endswith('-progress.log')}
(O/'artifact-manifest.json').write_text(json.dumps(artifacts,indent=2)+'\n');report['artifact_manifest_sha256']=sha(O/'artifact-manifest.json');(O/'summary.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(gates=report['recovery']['gates'],promotion_recommended=report['recovery']['promotion_recommended'],main_unchanged=report['main_unchanged']),indent=2))
