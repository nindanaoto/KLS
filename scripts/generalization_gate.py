#!/usr/bin/env python3
"""Independent nonregression gate; shared failures never become successful times."""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
import math
import statistics
from pathlib import Path

POLICY = 'generalization-v1'
PEER_POLICY = 'generalization-v2'
PEERS = ('ck', 'st')

def confirmation_peer_wins(first,second):
    """A win seen in either control run remains protected in both runs."""
    if not all(d.get('complete') for d in (first,second)):
        return dict(passed=False,protected_wins=0,lost_peer_wins=[],missing_peer_comparisons=[],reason='incomplete-campaign')
    fields=('config','matrix','reference','peer')
    keyed=[{tuple(c[k] for k in fields):c for c in d['peer_comparisons']} for d in (first,second)]
    wins={key for run in keyed for key,c in run.items() if c['reference_ratio']<1}
    lost=[];missing=[]
    for label,run in zip(('first','second'),keyed):
        for key in sorted(wins):
            if key not in run:missing.append(dict(zip(fields,key),campaign=label))
            elif run[key]['candidate_ratio']>1:lost.append(dict(run[key],campaign=label))
    return dict(passed=not lost and not missing,protected_wins=len(wins),
                lost_peer_wins=lost,missing_peer_comparisons=missing)

def matching_confirmation(first, second):
    fields=('policy','arms','configs','matrices','passes','binaries',
            'matrix_hashes','scripts','image_sha256','timeout')
    profiled=any(value=='--tuning-profile' or value.startswith('--tuning-profile=')
                 for config in first.get('configs',[]) for arm in first.get('arms',[])
                 if arm not in PEERS
                 for value in config.get(arm+'_args',config.get('kls_args',[])))
    if profiled and (not first.get('profile_hashes') or not second.get('profile_hashes')):
        return False
    return (first.get('policy') in (POLICY,PEER_POLICY) and second.get('policy')==first.get('policy') and
            all(k in first and k in second and first[k]==second[k] for k in fields) and
            first.get('profile_hashes',{})==second.get('profile_hashes',{}) and
            (first.get('policy')!=PEER_POLICY or
             ('external_hashes' in first and 'external_hashes' in second and
              first['external_hashes']==second['external_hashes'])) and
            {first.get('order'),second.get('order')}=={'initial','confirmation'})

def valid(row, systems):
    result = row.get('result')
    if (row.get('returncode') != 0 or row.get('timed_out') or not isinstance(result, dict) or
        row.get('closure',{}).get('all_session_processes_absent') is not True):
        return False
    # Peer harnesses report failures through the process status and residuals;
    # their successful JSON has no separate solver-status field.
    if result.get('status',0 if row.get('arm') in PEERS else None) != 0 or result.get('lifecycle_mode') != 'direct':
        return False
    if result.get('lifecycle_systems') != systems or result.get('verify_each_refactor') is not True:
        return False
    for field in ('relative_residual_l2', 'refactor_max_relative_residual'):
        value = result.get(field)
        if type(value) not in (float, int) or not math.isfinite(value) or not 0 <= value <= 1e-8:
            return False
    value = result.get('measured_lifecycle_seconds')
    return type(value) in (int, float) and math.isfinite(value) and value > 0

def evaluate(rows, specification):
    policy=specification.get('policy',POLICY)
    passes = specification['passes']
    arms = specification['arms']
    references=[a for a in arms if a not in ('candidate',*PEERS)]
    configs = {c['id']: c for c in specification['configs']}
    matrices = specification['matrices']
    if (type(passes) is not int or passes < 8 or passes % 4 or
        len(set(arms))!=len(arms) or not matrices or len(set(matrices))!=len(matrices) or
        not configs or len(configs)!=len(specification['configs']) or
        any(type(c.get('lifecycle_systems')) is not int or c['lifecycle_systems']<1
            for c in configs.values()) or
        'candidate' not in arms or 'baseline' not in arms or
        policy not in (POLICY,PEER_POLICY) or
        (policy==PEER_POLICY and (not set(PEERS)<=set(arms) or
          not set(arms)<=set(('baseline','candidate','previous',*PEERS))))):
        return dict(policy=policy, passed=False, reason='invalid-specification')
    expected = set(itertools.product(configs, matrices, range(passes), arms))
    indexed = {}
    for row in rows:
        key = (row.get('config'), row.get('matrix'), row.get('pass'), row.get('arm'))
        if key in indexed:
            return dict(policy=policy, passed=False, reason='duplicate-observation')
        indexed[key] = row
    if set(indexed) != expected:
        return dict(policy=policy, passed=False, reason='incomplete-or-unexpected-observations',
                    missing=len(expected-set(indexed)), unexpected=len(set(indexed)-expected))
    cells, failures, shared, reference_invalid = [], [], [], []
    for config, matrix, reference in itertools.product(configs, matrices, references):
        ratios, shared_passes = [], []
        for p in range(passes):
            a = indexed[config, matrix, p, reference]
            b = indexed[config, matrix, p, 'candidate']
            va, vb = valid(a, configs[config]['lifecycle_systems']), valid(b, configs[config]['lifecycle_systems'])
            if va and vb:
                ratios.append(b['result']['measured_lifecycle_seconds']/a['result']['measured_lifecycle_seconds'])
            elif va:
                failures.append(dict(config=config, matrix=matrix, reference=reference, pass_number=p))
            elif not vb:
                shared_passes.append(p)
            else:
                reference_invalid.append(dict(config=config, matrix=matrix,
                                               reference=reference, pass_number=p))
        if shared_passes:
            shared.append(dict(config=config, matrix=matrix, reference=reference, passes=shared_passes))
        if len(ratios) == passes:
            cells.append(dict(config=config, matrix=matrix, reference=reference,
                              ratio=statistics.median(ratios)))
    aggregates = []
    for config, reference in itertools.product(configs, references):
        values = [r['ratio'] for r in cells if r['config']==config and r['reference']==reference]
        aggregates.append(dict(config=config, reference=reference, covered=len(values), total=len(matrices),
                               geomean=math.exp(statistics.mean(map(math.log, values))) if values else None))
    regressions = [r for r in cells if r['ratio']>1.02]
    peer_cells, missing_peers, lost_wins = [], [], []
    if policy==PEER_POLICY:
        for cell,peer in itertools.product(cells,PEERS):
            config,matrix,reference=(cell[k] for k in ('config','matrix','reference'))
            peers=[indexed[config,matrix,p,peer] for p in range(passes)]
            invalid=[p for p,row in enumerate(peers) if not valid(row,configs[config]['lifecycle_systems'])]
            if invalid:
                missing_peers.append(dict(config=config,matrix=matrix,reference=reference,peer=peer,passes=invalid))
                continue
            times=[r['result']['measured_lifecycle_seconds'] for r in peers]
            prior=statistics.median(indexed[config,matrix,p,reference]['result']['measured_lifecycle_seconds']/t
                                    for p,t in enumerate(times))
            current=statistics.median(indexed[config,matrix,p,'candidate']['result']['measured_lifecycle_seconds']/t
                                      for p,t in enumerate(times))
            comparison=dict(config=config,matrix=matrix,reference=reference,peer=peer,
                            reference_ratio=prior,candidate_ratio=current,previously_winning=prior<1)
            peer_cells.append(comparison)
            if prior<1 and current>1:lost_wins.append(comparison)
    per_case_ok=(not lost_wins and not missing_peers) if policy==PEER_POLICY else not regressions
    passed = not failures and per_case_ok and all(a['geomean'] is not None and a['geomean']<=1 for a in aggregates)
    return dict(policy=policy, passed=passed, complete=True, cells=cells, aggregates=aggregates,
                candidate_only_failures=failures, shared_invalid=shared,
                reference_invalid=reference_invalid, regressions=regressions,
                peer_comparisons=peer_cells, missing_peer_evidence=missing_peers,lost_peer_wins=lost_wins,
                protected_peer_wins=sum(c['previously_winning'] for c in peer_cells),
                complete_valid_coverage=(len(cells)==len(configs)*len(matrices)*len(references) and not missing_peers))

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--campaign', type=Path, required=True)
    p.add_argument('--confirmation', type=Path)
    a=p.parse_args()
    specification=json.loads((a.campaign/'campaign.json').read_text())
    rows=[json.loads(line) for line in (a.campaign/'observations.jsonl').read_text().splitlines()]
    result=evaluate(rows,specification)
    if a.confirmation:
        other=json.loads((a.confirmation/'campaign.json').read_text())
        matching=matching_confirmation(specification,other)
        confirmation=evaluate([json.loads(line) for line in (a.confirmation/'observations.jsonl').read_text().splitlines()],other)
        result.update(confirmation=confirmation, matching_confirmation=matching)
        result['passed'] &= matching and confirmation['passed']
        if specification.get('policy')==PEER_POLICY:
            result['cross_campaign_peer_check']=confirmation_peer_wins(result,confirmation)
            result['passed'] &= result['cross_campaign_peer_check']['passed']
    else:
        result['confirmation_required']=True
    # The retained observations also support a stricter independent review
    # without relabeling the launcher's original per-run summaries.
    result['reviewer_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    path=a.campaign/('confirmed.json' if a.confirmation else 'summary.json')
    path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result.get(k) for k in ('policy','passed','complete','confirmation_required')}))
    return 0 if result['passed'] else 1
if __name__=='__main__':
    raise SystemExit(main())
