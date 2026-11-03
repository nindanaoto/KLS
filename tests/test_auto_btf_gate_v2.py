import copy
import importlib.util
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from auto_btf_merge_gate_v2 import evaluate, CASES, DOMAINS
import run_auto_btf_campaign_v2 as campaign
from audit_auto_btf_campaign_v2 import hybrid_acceptance


class RevisedGateTests(unittest.TestCase):
    def setUp(self):
        self.cells = [dict(domain=d, case=c, threads=t, peer=p, complete=True,
            median=(.97 if 'invert50' in c and t==1 else .99) if p=='baseline' else 1.2,
            wall_median=1., baseline_peer_median=1.3)
            for d in DOMAINS for c in CASES for t in (1,8) for p in ('baseline','ck','st')]
        self.row=next(c for c in self.cells if c['case'].endswith('/dac.cir') and c['peer']=='baseline')

    def check(self, expected):
        primary=evaluate(self.cells,all_qualified=True)
        self.assertEqual(primary,hybrid_acceptance(self.cells,True))
        self.assertEqual(primary['passed'],expected)

    def peers(self):
        return [c for c in self.cells if (c['domain'],c['case'],c['threads'])==
                (self.row['domain'],self.row['case'],self.row['threads']) and c['peer']!='baseline']

    def test_small_existing_loss_regression_accepted(self):
        self.row['median']=1.0257
        self.check(True)
        self.row['median']=1.05
        self.check(True)
        self.row['median']=1.050001
        self.check(False)

    def test_existing_win_preserved_and_two_percent_limit(self):
        p=self.peers()[0];p.update(baseline_peer_median=.9,median=.95)
        self.row['median']=1.02
        self.check(True)
        self.row['median']=1.020001
        self.check(False)
        self.row['median']=1.
        p['median']=1.000001
        self.check(False)

    def test_wall_uses_same_class_limit(self):
        self.row['wall_median']=1.05
        self.check(True)
        self.row['wall_median']=1.050001
        self.check(False)
        self.row['wall_median']=1.03
        self.peers()[0].update(baseline_peer_median=.9,median=.95)
        self.check(False)

    def test_overall_slowdown_rejected(self):
        for c in self.cells:
            if c['peer']=='baseline' and 'invert50' not in c['case']:c['median']=1.04
        self.check(False)

    def test_target_still_requires_two_percent(self):
        next(c for c in self.cells if c['peer']=='baseline' and 'invert50' in c['case'] and c['threads']==1)['median']=.981
        self.check(False)

    def test_incomplete_or_invalid_evidence_rejected(self):
        original=copy.deepcopy(self.cells)
        for mutate in (lambda:self.cells.pop(),lambda:self.cells.append(copy.deepcopy(self.cells[0])),
                       lambda:self.peers()[0].pop('baseline_peer_median'),
                       lambda:self.peers()[0].update(baseline_peer_median=float('nan')),
                       lambda:self.row.update(complete=False)):
            self.cells=copy.deepcopy(original)
            self.row=next(c for c in self.cells if c['case'].endswith('/dac.cir') and c['peer']=='baseline')
            mutate()
            self.assertFalse(evaluate(self.cells,all_qualified=True)['passed'])
        self.assertFalse(evaluate(original,all_qualified=False)['passed'])

    def test_reducer_computes_paired_main_peer_median(self):
        # Import the unchanged full-grid fixture; only the new reducer is used.
        spec=importlib.util.spec_from_file_location('fixture',Path(__file__).with_name('test_spice_campaign.py'))
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        metadata,rows=module.CampaignTests().fixture();metadata['acceptance_policy']='auto-btf-v2'
        for r in rows:
            if r['side']=='ck' and r['trial']>=0:r['linear_seconds']=[.5 if r['trial']%2 else 2.]
        result=campaign.reduce_results(metadata,rows)
        self.assertTrue(all(c['baseline_peer_median']==1.25 for c in result['cells'] if c['peer']=='ck'))
        self.assertEqual(result['hybrid'],hybrid_acceptance(result['cells'],True))
        metadata['diagnostic_only']=True
        self.assertFalse(campaign.reduce_results(metadata,rows)['acceptance_gate'])


if __name__=='__main__':unittest.main()
