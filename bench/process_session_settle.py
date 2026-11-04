"""Wait for post-exit descendants, then use the strict mount-helper cleanup.

Unknown processes are never signalled. The numerical deadline is unchanged:
this bounded wait begins only after the session leader has been reaped.
"""
import importlib.util
import os
from pathlib import Path
import time

_spec = importlib.util.spec_from_file_location(
    'strict_session_cleanup', Path(__file__).with_name('process_session_cleanup.py'))
_base = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_base)


def drain_exited_apptainer_session(session_id, image, *, settle_seconds=10.0,
                                    grace_seconds=0.25, termination_seconds=5.0):
    if not isinstance(session_id, int) or session_id <= 0:
        raise ValueError('invalid child session ID')
    if session_id == os.getsid(0) or Path('/proc', str(session_id)).exists():
        raise RuntimeError('child session leader has not been reaped')
    if settle_seconds < 0:
        raise ValueError('invalid post-exit settle interval')
    image_stat = Path(image).stat()
    identity = image_stat.st_dev, image_stat.st_ino
    start = time.monotonic()
    deadline = start + settle_seconds
    observed = {}
    while True:
        unknown = []
        for member in _base.session_members(session_id):
            try:
                pidfd = _base._validated_pidfd(member, identity)
            except (FileNotFoundError, ProcessLookupError):
                continue
            except RuntimeError as error:
                if str(error) != 'unrecognized process in completed child session':
                    raise
                unknown.append(member)
                observed[member['pid'], member['starttime']] = member
            else:
                os.close(pidfd)
        if not unknown:
            try:
                result = _base.drain_exited_apptainer_session(
                    session_id, image, grace_seconds=grace_seconds,
                    termination_seconds=termination_seconds)
            except RuntimeError as error:
                # The strict helper validates the entire set before signalling.
                # A newly observed unknown process therefore still receives no
                # signal. Any other validation/termination failure remains fatal.
                if str(error) != 'unrecognized process in completed child session':
                    raise
            else:
                result['unrecognized_processes_waited_for'] = list(observed.values())
                result['post_exit_settle_seconds'] = time.monotonic() - start
                return result
        if time.monotonic() >= deadline:
            raise RuntimeError('unrecognized processes did not exit during teardown: ' + repr(unknown))
        time.sleep(min(0.01, max(0, deadline - time.monotonic())))
