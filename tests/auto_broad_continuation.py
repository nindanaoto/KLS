"""Reject changed inputs and selective reuse of a recorded broad prefix."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('broad',ROOT/'scripts/run_auto_broad.py')
BROAD=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(BROAD)


class PrefixContract(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(prefix='kls-prefix-contract-')
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.required={'expected':1760,'passes':4,'timeout_seconds':900,'residual_limit':1e-8}
        self.expected=[];self.rows=[]
        for index,name in enumerate(('first','second','fullchip')):
            for trial in range(4):
                order=('baseline','candidate') if (index+trial)%2==0 else ('candidate','baseline')
                for position,side in enumerate(order):
                    identity=('headline_llc32_t8',name,trial,side,position)
                    self.expected.append(identity)
                    if index==2:continue
                    directory=self.root/f'{index}-{trial}-{side}';directory.mkdir()
                    for filename in ('stdout.log','stderr.log'):(directory/filename).write_text(filename)
                    self.rows.append(dict(zip(('config_id','matrix_id','trial','side','position'),identity),
                        valid=trial!=1,timed_out=trial==1,directory=str(directory),
                        closure={'all_session_processes_absent':True},command=[str(self.root/side)],
                        stdout_sha256=BROAD.sha(directory/'stdout.log'),stderr_sha256=BROAD.sha(directory/'stderr.log')))
        for side in ('baseline','candidate'):(self.root/side).write_text('frozen '+side)
        BROAD.save(self.root/'build-commands.json',[])
        BROAD.save(self.root/'metadata.json',dict(self.required,
            binaries={side:str(self.root/side) for side in ('baseline','candidate')},
            build_commands_sha256=BROAD.sha(self.root/'build-commands.json'),
            input_sha256={str(self.root/side):BROAD.sha(self.root/side) for side in ('baseline','candidate')}))
        (self.root/'observations.jsonl').write_text(''.join(json.dumps(row)+'\n' for row in self.rows))
        self.failure=self.root/'failure.json'
        BROAD.save(self.failure,dict(campaign=str(self.root),terminal_exit_code=1,
            completed_observations=len(self.rows),ledger_sha256=BROAD.sha(self.root/'observations.jsonl'),
            failure='RuntimeError: unrecognized process in completed child session',
            stage=dict(config='headline_llc32_t8',matrix='fullchip',trial=0,side='baseline',
                session_id=987654321,numeric_outcome='not persisted before cleanup raised'),
            later_strict_cleanup=dict(session_id=987654321,all_session_processes_absent=True)))

    def load(self):
        return BROAD.load_completed_prefix(self.root,self.expected,self.required,self.failure)

    def test_all_failed_rows_and_original_bytes_survive(self):
        _,rows,data=self.load()
        self.assertEqual(rows,self.rows)
        self.assertEqual(data,(self.root/'observations.jsonl').read_bytes())
        self.assertEqual(sum(not r['valid'] for r in rows),4)

    def test_changed_binary_is_rejected(self):
        (self.root/'candidate').write_text('different candidate')
        with self.assertRaises(AssertionError):self.load()

    def test_changed_output_is_rejected(self):
        (Path(self.rows[0]['directory'])/'stdout.log').write_text('changed')
        with self.assertRaises(AssertionError):self.load()

    def test_partial_matrix_is_rejected(self):
        p=self.root/'observations.jsonl';p.write_text(''.join(json.dumps(r)+'\n' for r in self.rows[:-1]))
        with self.assertRaises(AssertionError):self.load()

    def test_selective_or_reordered_rows_are_rejected(self):
        p=self.root/'observations.jsonl';rows=self.rows.copy();rows[0],rows[1]=rows[1],rows[0]
        p.write_text(''.join(json.dumps(r)+'\n' for r in rows))
        with self.assertRaises(AssertionError):self.load()


if __name__=='__main__':unittest.main()
