"""Exercise real post-exit sessions without signalling unknown descendants."""
import ctypes
import importlib.util
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import threading
import time
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('settle', ROOT / 'bench/process_session_settle.py')
SETTLE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SETTLE)

LEADER = r'''
import ctypes,json,os,signal,sys,time
from pathlib import Path
root=Path(sys.argv[1]); mode=sys.argv[2]
def mark(kind):
    with (root/'events.jsonl').open('a') as f:
        f.write(json.dumps(dict(kind=kind,time=time.monotonic_ns()))+'\n')
children=[]
pid=os.fork()
if pid==0:
    image=root/('wrong' if mode=='wrong-image' else 'image')
    fd=os.open(image,os.O_RDONLY)
    os.dup2(fd,3)
    if fd!=3:os.close(fd)
    def stop(*args):mark('helper-term');os._exit(0)
    signal.signal(signal.SIGTERM,stop)
    ctypes.CDLL(None).prctl(15,ctypes.c_char_p(b'squashfuse_ll'),0,0,0)
    while True:time.sleep(1)
children.append(pid)
if mode!='wrong-image':
    pid=os.fork()
    if pid==0:
        ctypes.CDLL(None).prctl(15,ctypes.c_char_p(b'exit-pending'),0,0,0)
        time.sleep(0.4 if mode=='transient' else 30)
        mark('unknown-natural-exit');os._exit(0)
    children.append(pid)
(root/'children.json').write_text(json.dumps(children))
'''


@unittest.skipUnless(sys.platform.startswith('linux') and hasattr(os, 'pidfd_open')
                     and hasattr(signal, 'pidfd_send_signal'), 'requires Linux pidfds')
class SessionSettle(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.libc = ctypes.CDLL(None, use_errno=True)
        cls.prior = ctypes.c_int()
        if cls.libc.prctl(37, ctypes.byref(cls.prior), 0, 0, 0):
            raise OSError(ctypes.get_errno(), 'get child subreaper')
        if cls.libc.prctl(36, 1, 0, 0, 0):
            raise OSError(ctypes.get_errno(), 'set child subreaper')

    @classmethod
    def tearDownClass(cls):
        cls.libc.prctl(36, cls.prior.value, 0, 0, 0)

    def exercise(self, mode, action):
        with tempfile.TemporaryDirectory(prefix='kls-session-settle-') as directory:
            root = Path(directory)
            (root/'image').write_bytes(b'owned image')
            (root/'wrong').write_bytes(b'other image')
            children, pidfds = [], []
            stop = threading.Event()
            reaper = None
            with (root/'leader.stderr').open('w') as errors:
                leader = subprocess.Popen([sys.executable, '-c', LEADER, str(root), mode],
                    start_new_session=True, stdout=subprocess.DEVNULL, stderr=errors)
                try:
                    leader.wait(timeout=5)
                    self.assertEqual(leader.returncode, 0, (root/'leader.stderr').read_text())
                    children = json.loads((root/'children.json').read_text())
                    for pid in children:
                        try:pidfds.append(os.pidfd_open(pid, 0))
                        except ProcessLookupError:pass
                    def reap():
                        while not stop.is_set():
                            for pid in children:
                                try:os.waitpid(pid, os.WNOHANG)
                                except ChildProcessError:pass
                            stop.wait(.005)
                    reaper = threading.Thread(target=reap)
                    reaper.start()
                    action(leader.pid, root)
                finally:
                    if leader.poll() is None:
                        leader.kill();leader.wait()
                    if not children and (root/'children.json').exists():
                        children = json.loads((root/'children.json').read_text())
                        for pid in children:
                            try:pidfds.append(os.pidfd_open(pid, 0))
                            except ProcessLookupError:pass
                    for pidfd in pidfds:
                        try:signal.pidfd_send_signal(pidfd, signal.SIGKILL, None, 0)
                        except ProcessLookupError:pass
                        os.close(pidfd)
                    deadline = time.monotonic()+5
                    while SETTLE._base.session_members(leader.pid) and time.monotonic()<deadline:
                        for pid in children:
                            try:os.waitpid(pid, os.WNOHANG)
                            except ChildProcessError:pass
                        time.sleep(.01)
                    stop.set()
                    if reaper:reaper.join()
                    self.assertEqual(SETTLE._base.session_members(leader.pid), [])

    def test_transient_unknown_exits_before_helper_signal(self):
        def action(sid, root):
            result = SETTLE.drain_exited_apptainer_session(sid, root/'image', settle_seconds=2)
            self.assertTrue(result['all_session_processes_absent'])
            self.assertTrue(result['unrecognized_processes_waited_for'])
            events = [json.loads(s) for s in (root/'events.jsonl').read_text().splitlines()]
            times = {e['kind']:e['time'] for e in events}
            self.assertLess(times['unknown-natural-exit'], times['helper-term'])
        self.exercise('transient', action)

    def test_persistent_unknown_is_not_signalled(self):
        def action(sid, root):
            with self.assertRaisesRegex(RuntimeError, 'did not exit during teardown'):
                SETTLE.drain_exited_apptainer_session(sid, root/'image', settle_seconds=.1)
            self.assertFalse((root/'events.jsonl').exists())
        self.exercise('persistent', action)

    def test_other_image_helper_is_not_signalled(self):
        def action(sid, root):
            with self.assertRaisesRegex(RuntimeError, 'holds another image'):
                SETTLE.drain_exited_apptainer_session(sid, root/'image', settle_seconds=2)
            self.assertFalse((root/'events.jsonl').exists())
        self.exercise('wrong-image', action)


if __name__ == '__main__':
    unittest.main()
