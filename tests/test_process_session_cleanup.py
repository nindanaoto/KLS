"""Real Linux sessions reproduce helpers escaping a launcher's process group.

The helpers only hold an image file; they never create a filesystem mount.
"""
import importlib.util
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('session_cleanup', ROOT / 'bench/process_session_cleanup.py')
cleanup = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cleanup)

HELPER = r'''
import ctypes,json,os,signal,sys,time
from pathlib import Path
pid=os.fork()
if pid:
    while not Path(sys.argv[2]).exists():time.sleep(0.005)
    os._exit(0)
os.setpgid(0,0)
fd=os.open(sys.argv[1],os.O_RDONLY)
os.dup2(fd,3)
if fd!=3:os.close(fd)
assert ctypes.CDLL(None).prctl(15,sys.argv[3].encode(),0,0,0)==0
signal.signal(signal.SIGTERM,lambda *_:os._exit(0))
out=Path(sys.argv[2]);pending=out.with_suffix('.tmp')
fields=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()
pending.write_text(json.dumps(dict(pid=os.getpid(),pgid=os.getpgrp(),sid=os.getsid(0),starttime=fields[19])))
os.replace(pending,out)
while True:signal.pause()
'''


class SessionCleanupTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='kls-session-cleanup-test-')
        self.folder = Path(self.temporary.name)
        self.image = ROOT / 'build/kls-xyce-private.sif'
        self.helpers = []

    def spawn(self, *, image=None, comm='squashfuse_ll'):
        output = self.folder / ('helper' + str(len(self.helpers)) + '.json')
        proc = subprocess.Popen([sys.executable, '-c', HELPER, str(image or self.image),
                                 str(output), comm], start_new_session=True,
                                stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                                stderr=subprocess.DEVNULL)
        self.assertEqual(proc.wait(timeout=5), 0)
        item = json.loads(output.read_text())
        self.helpers.append(item)
        self.assertEqual(item['sid'], proc.pid)
        self.assertNotEqual(item['pgid'], proc.pid)
        return item

    def tearDown(self):
        for item in self.helpers:
            process = Path('/proc', str(item['pid']))
            try:
                fields = (process / 'stat').read_text().rsplit(')', 1)[1].split()
                if fields[19] == item['starttime']:
                    os.kill(item['pid'], signal.SIGTERM)
            except (FileNotFoundError, ProcessLookupError):
                pass
        deadline = time.monotonic() + 5
        while any(Path('/proc', str(item['pid'])).exists() for item in self.helpers):
            self.assertLess(time.monotonic(), deadline, 'test helper not reaped')
            time.sleep(0.01)
        self.temporary.cleanup()

    def test_reaps_mount_helper_in_another_process_group(self):
        item = self.spawn()
        result = cleanup.drain_exited_apptainer_session(item['sid'], self.image)
        self.assertTrue(result['all_session_processes_absent'])
        self.assertEqual([p['pid'] for p in result['signalled_mount_helpers']], [item['pid']])
        self.assertFalse(Path('/proc', str(item['pid'])).exists())

    def test_rejects_another_image(self):
        other = self.folder / 'other.img'
        other.write_bytes(b'not the requested image')
        item = self.spawn(image=other)
        with self.assertRaisesRegex(RuntimeError, 'another image'):
            cleanup.drain_exited_apptainer_session(item['sid'], self.image)
        self.assertTrue(Path('/proc', str(item['pid'])).exists())

    def test_rejects_nonhelper_session_member(self):
        item = self.spawn(comm='not_a_helper')
        with self.assertRaisesRegex(RuntimeError, 'unrecognized process'):
            cleanup.drain_exited_apptainer_session(item['sid'], self.image)
        self.assertTrue(Path('/proc', str(item['pid'])).exists())

    def test_other_session_is_untouched(self):
        first, other = self.spawn(), self.spawn()
        cleanup.drain_exited_apptainer_session(first['sid'], self.image)
        self.assertTrue(Path('/proc', str(other['pid'])).exists())

    def test_empty_exited_session(self):
        proc = subprocess.Popen([sys.executable, '-c', 'pass'], start_new_session=True)
        self.assertEqual(proc.wait(timeout=5), 0)
        self.assertEqual(cleanup.drain_exited_apptainer_session(proc.pid, self.image),
                         dict(session_id=proc.pid, all_session_processes_absent=True,
                              signalled_mount_helpers=[]))


if __name__ == '__main__':
    unittest.main(verbosity=2)
