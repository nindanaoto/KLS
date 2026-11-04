#!/usr/bin/env python3
"""Reassess two retained opposite-order campaigns without modifying them.

Run after independent audit_auto_btf_campaign_v3.py audits, between timing runs.
"""
import argparse
import json
from pathlib import Path

from auto_btf_merge_gate_v3 import evaluate
from run_auto_btf_campaign_v2 import reduce_results, sha


def reassess(campaign, audit_path):
    metadata_path = campaign / 'metadata.json'
    ledger_path = campaign / 'observations.jsonl'
    metadata = json.loads(metadata_path.read_text())
    rows = [json.loads(line) for line in ledger_path.read_text().splitlines()]
    audit = json.loads(audit_path.read_text())
    assert audit['acceptance_policy'] == 'auto-btf-v3'
    assert audit['all_qualified'] and not audit['errors']
    for path, expected in audit['input_sha256'].items():
        assert sha(Path(path)) == expected, f'Changed audited input: {path}'
    for path in (metadata_path, ledger_path):
        assert audit['input_sha256'][str(path.resolve())] == sha(path)
    assert audit['assessor_sha256'] == sha(Path(__file__).with_name('audit_auto_btf_campaign_v3.py'))
    primary = reduce_results(metadata, rows)
    fields = ('median', 'wall_median', 'baseline_peer_median')
    def keyed(cells):
        return {(c['domain'], c['case'], c['threads'], c['peer']):
                tuple(c[field] for field in fields) for c in cells}
    assert keyed(primary['cells']) == keyed(audit['cells']), 'Independent cell disagreement'
    decision = evaluate(primary['cells'], all_qualified=primary['all_qualified'])
    assert decision == audit['hybrid'], 'Independent gate disagreement'
    return metadata, dict(campaign=str(campaign), acceptance_policy='auto-btf-v3',
                         recorded_acceptance_policy=metadata.get('acceptance_policy'),
                         policy_reassessment=True, decision=decision,
                         metadata_sha256=sha(metadata_path), ledger_sha256=sha(ledger_path),
                         independent_audit=str(audit_path), audit_sha256=sha(audit_path))


def matching_campaigns(metadata, campaigns):
    protocol_fields = ('threads', 'passes', 'timeout',
                       'primary_metric', 'wall_metric', 'accuracy',
                       'complete_grid', 'diagnostic_only')
    def common_inputs(m, root):
        return {p: h for p, h in m['input_sha256'].items()
                if not Path(p).is_relative_to(root / 'binaries')}
    return (campaigns[0] != campaigns[1]
            and [m['order'] for m in metadata] == ['initial', 'confirmation']
            and metadata[0]['sides'] == metadata[1]['sides'] == ['baseline', 'candidate', 'ck', 'st']
            and all(field in metadata[0] and field in metadata[1]
                    and metadata[0][field] == list(reversed(metadata[1][field]))
                    for field in ('cases', 'domains'))
            and all(field in metadata[0] and field in metadata[1]
                    and metadata[0][field] == metadata[1][field]
                    for field in protocol_fields)
            and common_inputs(metadata[0], campaigns[0]) == common_inputs(metadata[1], campaigns[1])
            and all(sha(campaigns[0] / 'binaries' / side) == sha(campaigns[1] / 'binaries' / side)
                    for side in ('baseline', 'candidate')))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaigns', nargs=2, type=Path, required=True)
    parser.add_argument('--audits', nargs=2, type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    campaigns = [p.resolve() for p in args.campaigns]
    # A fresh destination prevents replacing any earlier decision or raw campaign.
    args.output.mkdir(parents=True)
    results = [reassess(c, a.resolve()) for c, a in zip(campaigns, args.audits)]
    metadata, reports = zip(*results)
    matching = matching_campaigns(metadata, campaigns)
    result = dict(acceptance_policy='auto-btf-v3', policy_reassessment=True,
                  acceptance_gate=matching and all(r['decision']['passed'] for r in reports),
                  matching_inputs=matching, campaigns=reports,
                  tool_sha256={p.name: sha(p) for p in (
                      Path(__file__), Path(__file__).with_name('auto_btf_merge_gate_v3.py'),
                      Path(__file__).with_name('run_auto_btf_campaign_v2.py'))})
    (args.output / 'confirmation.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(dict(acceptance_gate=result['acceptance_gate'], matching_inputs=matching)))
    return 0 if result['acceptance_gate'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
