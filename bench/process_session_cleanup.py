"""Close an exited Apptainer child's session without leaving mount helpers.

The numerical child must already have been waited for. Unknown processes are
never signalled: only same-user squashfuse helpers holding the requested image
on fd 3 are eligible. A pidfd binds each signal to the validated process.
"""
import os
from pathlib import Path
import signal
import time


def session_members(session_id):
    members = []
    for process in Path('/proc').iterdir():
        if not process.name.isdigit():
            continue
        try:
            fields = (process / 'stat').read_text().rsplit(')', 1)[1].split()
            if int(fields[3]) == session_id:
                members.append(dict(pid=int(process.name), pgid=int(fields[2]),
                                    sid=int(fields[3]), state=fields[0],
                                    starttime=fields[19],
                                    comm=(process / 'comm').read_text().strip()))
        except (FileNotFoundError, ProcessLookupError):
            pass
    return sorted(members, key=lambda item: item['pid'])


def _validated_pidfd(member, image_identity):
    pidfd = os.pidfd_open(member['pid'], 0)
    try:
        process = Path('/proc', str(member['pid']))
        fields = (process / 'stat').read_text().rsplit(')', 1)[1].split()
        uid = int(next(line for line in (process / 'status').read_text().splitlines()
                       if line.startswith('Uid:')).split()[1])
        if (fields[19] != member['starttime'] or int(fields[3]) != member['sid'] or
                (process / 'comm').read_text().strip() != 'squashfuse_ll' or
                uid != os.getuid()):
            raise RuntimeError('unrecognized process in completed child session')
        image = (process / 'fd/3').stat()
        if (image.st_dev, image.st_ino) != image_identity:
            raise RuntimeError('mount helper holds another image')
        return pidfd
    except BaseException:
        os.close(pidfd)
        raise


def drain_exited_apptainer_session(session_id, image, *, grace_seconds=0.25,
                                    termination_seconds=5.0):
    """Return only when every PID/PGID in this session has disappeared.

    This is teardown after numerical completion/timeout, not extra numerical
    time. The caller retains the original exit status and timeout outcome.
    """
    if not isinstance(session_id, int) or session_id <= 0:
        raise ValueError('invalid child session ID')
    if session_id == os.getsid(0) or Path('/proc', str(session_id)).exists():
        raise RuntimeError('child session leader has not been reaped')
    deadline = time.monotonic() + grace_seconds
    members = session_members(session_id)
    while members and time.monotonic() < deadline:
        time.sleep(0.01)
        members = session_members(session_id)
    signalled = []
    if members:
        image_stat = Path(image).stat()
        identity = image_stat.st_dev, image_stat.st_ino
        pidfds = []
        try:
            # Validate the whole remaining set before sending any signal.
            for member in members:
                try:
                    pidfds.append((member, _validated_pidfd(member, identity)))
                except (FileNotFoundError, ProcessLookupError):
                    continue
            for member, pidfd in pidfds:
                try:
                    signal.pidfd_send_signal(pidfd, signal.SIGTERM, None, 0)
                    signalled.append(member)
                except ProcessLookupError:
                    pass
        finally:
            for _, pidfd in pidfds:
                os.close(pidfd)
        deadline = time.monotonic() + termination_seconds
        members = session_members(session_id)
        while members and time.monotonic() < deadline:
            time.sleep(0.01)
            members = session_members(session_id)
        if members:
            raise RuntimeError('completed child session still has processes: ' + repr(members))
    return dict(session_id=session_id, all_session_processes_absent=True,
                signalled_mount_helpers=signalled)
