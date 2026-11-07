import copy
import sys
import unittest
import tempfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from generalization_gate import evaluate, matching_confirmation, confirmation_peer_wins, PEER_POLICY
from run_generalization_campaign import tuning_profiles, launch_order

class GeneralizationGateTests(unittest.TestCase):
    def setUp(self):
        self.spec=dict(passes=8,arms=['baseline','candidate','previous'],configs=[dict(id='t1',lifecycle_systems=10)],matrices=['a','b'])
        self.rows=[dict(config='t1',matrix=m,arm=arm,**{'pass':p},returncode=0,timed_out=False,closure={'all_session_processes_absent':True},
                       result=dict(status=0,lifecycle_mode='direct',lifecycle_systems=10,verify_each_refactor=True,
                                   relative_residual_l2=1e-14,refactor_max_relative_residual=1e-14,
                                   measured_lifecycle_seconds=.99 if arm=='candidate' else 1))
                   for m in ('a','b') for p in range(8) for arm in self.spec['arms']]
    def result(self):return evaluate(self.rows,self.spec)
    def test_complete_improvement(self):self.assertTrue(self.result()['passed'])
    def test_confirmation_requires_same_profile_and_opposite_order(self):
        first=dict(self.spec,policy='generalization-v1',binaries={},matrix_hashes={},scripts={},
                   image_sha256='image',timeout=60,order='initial',profile_hashes={'model':'abc'})
        second=dict(first,order='confirmation')
        self.assertTrue(matching_confirmation(first,second))
        second['profile_hashes']={'model':'changed'}
        self.assertFalse(matching_confirmation(first,second))
        second=dict(first)
        self.assertFalse(matching_confirmation(first,second))
        second=dict(first,order='unknown')
        self.assertFalse(matching_confirmation(first,second))
        first['configs']=[dict(self.spec['configs'][0],candidate_args=['--tuning-profile','model'])]
        second=dict(first,order='confirmation')
        self.assertTrue(matching_confirmation(first,second))
        first.pop('profile_hashes');second.pop('profile_hashes')
        self.assertFalse(matching_confirmation(first,second))
    def test_only_used_profile_arguments_are_frozen(self):
        with tempfile.TemporaryDirectory() as directory:
            common=Path(directory)/'common';common.touch()
            candidate=Path(directory)/'candidate';candidate.touch()
            config=dict(kls_args=['--tuning-profile',str(common)],
                        candidate_args=['--tuning-profile',str(candidate)],
                        previous_args=['--tuning-profile',str(Path(directory)/'unused')])
            self.assertEqual(tuning_profiles([config],['baseline','candidate']),{common,candidate})
            with self.assertRaises(ValueError):
                tuning_profiles([dict(kls_args=['--tuning-profile'])],['candidate'])
    def test_empty_or_invalid_configuration_cannot_pass(self):
        self.assertFalse(evaluate([],dict(self.spec,configs=[]))['passed'])
        self.spec['configs'][0]['lifecycle_systems']=0
        self.assertFalse(self.result()['passed'])
    def test_invalid_reference_is_disclosed(self):
        for r in self.rows:
            if r['matrix']=='a' and r['arm']=='baseline':r['returncode']=1
        result=self.result()
        self.assertTrue(result['passed'])
        self.assertFalse(result['complete_valid_coverage'])
        self.assertEqual(len(result['reference_invalid']),8)
    def test_incomplete_and_duplicate(self):
        self.rows.pop();self.assertFalse(self.result()['passed'])
        self.rows.append(copy.deepcopy(self.rows[0]));self.assertFalse(self.result()['passed'])
    def test_regression_not_hidden_by_geomean(self):
        for r in self.rows:
            if r['arm']=='candidate':r['result']['measured_lifecycle_seconds']=1.03 if r['matrix']=='a' else .8
        self.assertFalse(self.result()['passed'])
    def test_aggregate_slowdown_below_per_case_cap(self):
        for r in self.rows:
            if r['arm']=='candidate':r['result']['measured_lifecycle_seconds']=1.001
        self.assertFalse(self.result()['passed'])
    def test_candidate_only_failure_and_timeout(self):
        r=next(r for r in self.rows if r['arm']=='candidate');r['timed_out']=True
        self.assertFalse(self.result()['passed'])
    def test_shared_failures_disclosed(self):
        for r in self.rows:
            if r['matrix']=='a':r['returncode']=1
        result=self.result();self.assertTrue(result['passed']);self.assertFalse(result['complete_valid_coverage']);self.assertTrue(result['shared_invalid'])
    def test_previous_stage_protects_cumulative_loss(self):
        for r in self.rows:
            if r['arm']=='previous':r['result']['measured_lifecycle_seconds']=.9
        self.assertFalse(self.result()['passed'])
    def test_nonfinite_and_unverified_rejected(self):
        for value in (float('nan'),float('inf'),-1):
            original=copy.deepcopy(self.rows)
            next(r for r in self.rows if r['arm']=='candidate')['result']['relative_residual_l2']=value
            self.assertFalse(self.result()['passed']);self.rows=original
        next(r for r in self.rows if r['arm']=='candidate')['result']['verify_each_refactor']=False
        self.assertFalse(self.result()['passed'])

    def use_peer_policy(self):
        self.spec['policy']=PEER_POLICY
        self.spec['arms']+=['ck','st']
        for row in list(self.rows):
            if row['arm']!='baseline':continue
            for peer in ('ck','st'):
                new=copy.deepcopy(row);new['arm']=peer
                new['result'].pop('status')
                new['result']['measured_lifecycle_seconds']=1.5
                self.rows.append(new)

    def test_peer_policy_allows_slowdown_inside_previous_winning_margin(self):
        self.use_peer_policy()
        for row in self.rows:
            if row['arm']=='candidate':row['result']['measured_lifecycle_seconds']=1.2 if row['matrix']=='a' else .7
        result=self.result()
        self.assertTrue(result['passed']);self.assertTrue(result['regressions'])
        self.assertEqual(result['lost_peer_wins'],[])
        self.assertEqual(result['protected_peer_wins'],8)
        for row in self.rows:
            if row['matrix']=='a' and row['arm']=='candidate':row['result']['measured_lifecycle_seconds']=1.5
            if row['matrix']=='b' and row['arm']=='candidate':row['result']['measured_lifecycle_seconds']=.5
        self.assertTrue(self.result()['passed'])  # Equal to the peer is not slower.

    def test_peer_policy_does_not_require_new_wins_on_previously_losing_cases(self):
        self.use_peer_policy()
        for row in self.rows:
            if row['arm'] in ('ck','st'):row['result']['measured_lifecycle_seconds']=.8
            elif row['arm']=='candidate':row['result']['measured_lifecycle_seconds']=1.2 if row['matrix']=='a' else .7
        result=self.result()
        self.assertTrue(result['passed']);self.assertEqual(result['protected_peer_wins'],0)

    def test_peer_policy_protects_each_peer_individually(self):
        self.use_peer_policy()
        for row in self.rows:
            if row['matrix']=='a':
                row['result']['measured_lifecycle_seconds']={'candidate':1.3,'baseline':1.,'previous':1.,'ck':1.2,'st':.8}[row['arm']]
            elif row['arm']=='candidate':row['result']['measured_lifecycle_seconds']=.5
        result=self.result()
        self.assertFalse(result['passed'])
        self.assertEqual({c['peer'] for c in result['lost_peer_wins']},{'ck'})

    def test_peer_policy_preserves_previous_control_wins(self):
        self.use_peer_policy()
        for row in self.rows:
            if row['arm']=='previous':row['result']['measured_lifecycle_seconds']=.8
            elif row['arm'] in ('ck','st'):row['result']['measured_lifecycle_seconds']=.9
            elif row['arm']=='candidate':row['result']['measured_lifecycle_seconds']=.95 if row['matrix']=='a' else .5
        result=self.result()
        self.assertFalse(result['passed'])
        self.assertEqual({c['reference'] for c in result['lost_peer_wins']},{'previous'})

    def test_peer_policy_keeps_aggregate_requirement(self):
        self.use_peer_policy()
        for row in self.rows:
            if row['arm']=='candidate':row['result']['measured_lifecycle_seconds']=1.01
        self.assertFalse(self.result()['passed'])

    def test_peer_policy_does_not_invent_missing_peer_evidence(self):
        self.spec['policy']=PEER_POLICY
        self.assertFalse(self.result()['passed'])
        self.use_peer_policy()
        next(r for r in self.rows if r['arm']=='ck')['returncode']=1
        result=self.result()
        self.assertFalse(result['passed']);self.assertTrue(result['missing_peer_evidence'])
        self.assertFalse(result['complete_valid_coverage'])

    def test_peer_policy_needs_matched_peer_inputs_in_confirmation(self):
        self.use_peer_policy()
        first=dict(self.spec,binaries={},matrix_hashes={},scripts={},image_sha256='image',timeout=60,
                   order='initial',external_hashes={'peer-library':'abc'})
        second=dict(first,order='confirmation')
        self.assertTrue(matching_confirmation(first,second))
        second['external_hashes']={'peer-library':'changed'}
        self.assertFalse(matching_confirmation(first,second))
        second=dict(first,order='confirmation');second.pop('external_hashes')
        self.assertFalse(matching_confirmation(first,second))

    def test_all_arm_counts_have_balanced_pairwise_launch_order(self):
        import itertools
        for count in (2,3,4,5):
            arms=list('abcde')[:count]
            for confirmation in (False,True):
                orders=[launch_order(arms,p,confirmation) for p in range(8)]
                for a,b in itertools.combinations(arms,2):
                    self.assertEqual(sum(o.index(a)<o.index(b) for o in orders),4)
                self.assertEqual(arms,list('abcde')[:count])

    def test_confirmation_protects_a_win_found_in_either_run(self):
        self.use_peer_policy()
        first=copy.deepcopy(self.rows);second=copy.deepcopy(self.rows)
        for rows,prior,current in ((first,1.1,1.03),(second,.9,.95)):
            for row in rows:
                if row['matrix']=='a':
                    row['result']['measured_lifecycle_seconds']={'baseline':prior,'previous':prior,
                        'candidate':current,'ck':1.,'st':1.5}[row['arm']]
                elif row['arm']=='candidate':row['result']['measured_lifecycle_seconds']=.7
        a,b=evaluate(first,self.spec),evaluate(second,self.spec)
        self.assertTrue(a['passed']);self.assertTrue(b['passed'])
        for x,y in ((a,b),(b,a)):
            check=confirmation_peer_wins(x,y)
            self.assertFalse(check['passed'])
            self.assertEqual(len(check['lost_peer_wins']),2)
            self.assertEqual({c['peer'] for c in check['lost_peer_wins']},{'ck'})

    def test_confirmation_requires_an_observation_for_each_protected_win(self):
        self.use_peer_policy();a=self.result();b=copy.deepcopy(a)
        b['peer_comparisons'].pop()
        check=confirmation_peer_wins(a,b)
        self.assertFalse(check['passed']);self.assertTrue(check['missing_peer_comparisons'])
        self.assertFalse(confirmation_peer_wins(a,{'passed':False})['passed'])
if __name__=='__main__':unittest.main()
