"""Qualification gates and recoverable archive behavior, without running Xyce."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

path=Path(__file__).resolve().parents[1]/'bench/qualify_certificate.py'
spec=importlib.util.spec_from_file_location('campaign',path)
campaign=importlib.util.module_from_spec(spec)
spec.loader.exec_module(campaign)


class CampaignTests(unittest.TestCase):
    def rows(self,ratio):
        return [dict(campaign_phase='confirmation',circuit=cid,mode=mode,block=block,
                     variant=variant,xyce_elapsed_s=10*(ratio if variant=='candidate' else 1))
                for cid in campaign.PRIMARY+campaign.CONTROLS for mode in ('kls1','kls8')
                for block in range(9) for variant in ('baseline','candidate')]

    def test_performance_gates(self):
        self.assertTrue(campaign.summary(self.rows(.85),'confirmation',False)['passed'])
        self.assertFalse(campaign.summary(self.rows(.95),'confirmation',False)['passed'])
        self.assertTrue(campaign.summary(self.rows(.95),'confirmation',True)['passed'])
        self.assertFalse(campaign.summary(self.rows(1.03),'confirmation',True)['passed'])
        rows=self.rows(.85)
        for row in rows:
            if row['circuit']==campaign.CONTROLS[0] and row['variant']=='candidate':
                row['xyce_elapsed_s']=10.3
        self.assertFalse(campaign.summary(rows,'confirmation',False)['passed'])
        with self.assertRaises(ValueError):campaign.summary(rows[:-1],'confirmation',False)
        for bad in (float('inf'),float('nan'),0.,-1.):
            rows=self.rows(.85);rows[0]['xyce_elapsed_s']=bad
            with self.assertRaises(ValueError):campaign.summary(rows,'confirmation',False)

    def test_archive_resume_and_tamper(self):
        support=types.ModuleType('discovery_support')
        support.sha256=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
        support.write_json=lambda p,d:p.write_text(json.dumps(d))
        with tempfile.TemporaryDirectory() as temp,patch.dict(sys.modules,{'discovery_support':support}):
            root=Path(temp);raw=root/'Transformed_Matrix1.mm';raw.write_text('test capture\n')
            records=[dict(hashes={raw.name:support.sha256(raw)})]
            archive=root/'chunk.tar.gz';receipt=root/'receipt.json'
            campaign.archive_chunk(root,records,archive,receipt)
            self.assertFalse(raw.exists())
            self.assertTrue(json.loads(receipt.read_text())['raw_removed'])
            campaign.archive_chunk(root,records,archive,receipt)
            raw.write_text('changed capture\n')
            with self.assertRaises(ValueError):campaign.archive_chunk(root,records,archive,receipt)
            self.assertTrue(raw.exists())


if __name__=='__main__':unittest.main()
