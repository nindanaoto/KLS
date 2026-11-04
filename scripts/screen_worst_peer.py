"""Focused paired experiment; never supplies full-grid acceptance evidence.

Pass --case CASE followed by the ordinary campaign arguments. Run the screen
with --review OUTPUT to reparse its retained observations.
"""
from pathlib import Path
import sys

import run_worst_peer_campaign as campaign


def main():
    if '--case' not in sys.argv:
        raise SystemExit('--case is required, including for review')
    index = sys.argv.index('--case')
    case = sys.argv[index + 1]
    if case not in campaign.CASES:
        raise SystemExit('case must be an existing campaign circuit')
    del sys.argv[index:index + 2]
    campaign.CASES = (case,)
    original = campaign.write_json

    def write(path, value):
        if Path(path).name == 'metadata.json':
            value['complete_grid'] = False
            value['screen_only'] = True
            value['input_sha256'][str(Path(__file__).resolve())] = campaign.sha(__file__)
        return original(path, value)

    campaign.write_json = write
    return campaign.main()


if __name__ == '__main__':
    raise SystemExit(main())
