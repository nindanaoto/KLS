"""Unsafe flags must fail at the certificate's guard, not merely fail to link."""
from pathlib import Path
import subprocess
import sys
import tempfile


def main():
    if len(sys.argv) != 2:
        raise SystemExit('usage: test_accuracy_build_flags.py C_COMPILER')
    source = Path(__file__).resolve().parents[1] / 'src/kls_accuracy_certificate.c'
    marker = 'KLS accuracy certificate requires strict finite/nonfinite semantics'
    with tempfile.TemporaryDirectory(prefix='kls-accuracy-build-flags-') as tmp:
        target = Path(tmp) / 'certificate.o'
        for flags in ([], ['-ffast-math'], ['-ffinite-math-only']):
            command = [sys.argv[1], '-std=c11', '-O2', *flags,
                       '-c', str(source), '-o', str(target)]
            result = subprocess.run(command, capture_output=True, text=True, timeout=30)
            if not flags:
                if result.returncode:
                    raise RuntimeError('strict certificate did not compile:\n' + result.stderr)
            elif result.returncode == 0 or marker not in result.stderr:
                raise RuntimeError(f'{flags}: missing strict-arithmetic rejection:\n' + result.stderr)
    print('PASS certificate build flags: strict accepted; fast/finite-only rejected')


if __name__ == '__main__':
    main()
