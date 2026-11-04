"""Bind an original standalone campaign to its independent audit.

This validates retained receipts and pins; it does not run the auditor or
collect fresh peer timings. A final acceptance caller must do both.
"""
import hashlib
import json
import statistics
from pathlib import Path


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def bind_original(campaign, audit):
    campaign, audit = Path(campaign).resolve(), Path(audit).resolve()
    metadata = read(campaign / 'metadata.json')
    summary = read(campaign / 'summary.json')
    checked = read(audit)
    assert checked['observations'] == summary['observations'] == 1760
    assert summary['expected'] == metadata['expected'] == 1760
    assert checked['original_manifest_complete'] is True
    assert checked['independent_exact_agreement'] is True
    assert not summary['changed_inputs']
    assert not checked['coverage_regressions'] and not summary['coverage_regressions']
    assert checked['groups'] == summary['groups']
    auditor = Path(__file__).with_name('audit_auto_broad.py').resolve()
    assert checked['input_sha256'].get(str(auditor)) == sha(auditor)
    assert len(checked['invalid']) == len(summary['invalid'])
    assert metadata['matrix_count'] == 110 and metadata['passes'] == 4
    assert metadata['lifecycle_systems'] == 100
    assert metadata['timeout_seconds'] == 900
    assert metadata['residual_limit'] == 1e-8
    configs = metadata['configs']
    assert len(configs) == 2
    assert {(c['id'], c['cpus']) for c in configs} == {
        ('headline_llc32_t8', '8-15'), ('headline_llc96_t8', '0-7')}
    assert all(c['threads'] == 8 and c['lifecycle_systems'] == 100
               and c['values'] == 'entrywise' and c['amplitude'] == .001
               and c['env'] == {} for c in configs)
    matrices = metadata['matrix_paths']
    assert len(matrices) == 110
    expected = {(c['id'], name) for c in configs for name in matrices}
    groups = {(g['config'], g['matrix']): g for g in checked['groups']}
    assert len(groups) == len(checked['groups']) == 220
    assert set(groups) == expected
    from worst_peer_broad_gate import positive
    for group in groups.values():
        count, ratios = group['valid_pairs'], group['ratios']
        assert type(count) is int and 0 <= count <= 4
        assert isinstance(ratios, list) and len(ratios) == count
        assert all(positive(ratio) for ratio in ratios)
        median = statistics.median(ratios) if count == 4 else None
        assert group['median'] == median
    pins = dict(metadata['input_sha256'])
    for path, digest in checked['input_sha256'].items():
        assert path not in pins or pins[path] == digest
        pins[path] = digest
    for name in ('metadata.json', 'observations.jsonl', 'summary.json'):
        path = str(campaign / name)
        assert checked['input_sha256'].get(path) == sha(path)
    for path, digest in pins.items():
        assert sha(path) == digest, ('changed evidence', path)
    candidate = Path(metadata['binaries']['candidate'])
    assert pins.get(str(candidate)) == sha(candidate)
    for path in matrices.values():
        assert str(Path(path)) in pins
    pins[str(audit)] = sha(audit)
    pins[str(Path(__file__).resolve())] = sha(__file__)
    # Preserve original ratios. Fresh peer samples can justify an exception,
    # but cannot replace these with a more favorable baseline retest.
    from worst_peer_broad_gate import evaluate
    validation = evaluate(checked['groups'], [], expected_cells=sorted(expected),
                          coverage_qualified=True, peer_evidence_qualified=False)
    targets = [g for g in checked['groups']
               if g['median'] is not None and g['median'] > 1.02]
    if not targets:
        assert validation['complete'] and validation['passed']
    return dict(metadata=metadata, groups=checked['groups'],
                expected_cells=sorted(expected), targets=targets,
                candidate=str(candidate), input_sha256=pins)
