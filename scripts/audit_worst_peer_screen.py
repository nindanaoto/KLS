#!/usr/bin/env python3
"""Recheck focused diagnostics with the established Xyce waveform protocol.

Preserves the original ledger and strict pointwise/grid diagnostics. This is
not full-grid acceptance evidence; final campaigns and their audits are required.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import subprocess


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--screen', type=Path, required=True)
    parser.add_argument('--ledger', default='observations.jsonl')
    parser.add_argument('--domain', choices=('llc96', 'llc32'))
    parser.add_argument('--reference-campaign', type=Path, action='append', required=True)
    parser.add_argument('--regression', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    runner = Path(__file__).with_name('run_worst_peer_campaign.py')
    spec = importlib.util.spec_from_file_location('campaign', runner)
    campaign = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(campaign)
    verifier = args.regression / 'TestScripts/xyce_verify.pl'
    reference = {}
    for root in args.reference_campaign:
        for line in (root / 'observations.jsonl').read_text().splitlines():
            row = json.loads(line)
            if row['side'] == 'klu':
                key = (row['case'], row['domain'])
                if key in reference:
                    raise ValueError('duplicate reference: ' + repr(key))
                if not row['qualified']:
                    raise ValueError('unqualified reference: ' + repr(key))
                reference[key] = row
    original = args.screen / args.ledger
    args.output.mkdir(parents=True)
    script_archive = args.output / 'audit_source.py'
    script_archive.write_text(Path(__file__).read_text())
    results = []
    for index, line in enumerate(original.read_text().splitlines()):
        row = json.loads(line)
        case = row.get('case')
        if case is None:
            cases = {key[0] for key in reference}
            if len(cases) != 1:
                raise ValueError('case missing from ambiguous diagnostic')
            case = next(iter(cases))
        domain = row.get('domain', args.domain)
        if domain is None:
            raise ValueError('domain missing; supply --domain')
        ref = reference[(case, domain)]
        destinations = {Path(x) for x in row['command']
                        if str(x).startswith(str(args.screen) + '/') and Path(x).is_dir()}
        if len(destinations) != 1:
            raise ValueError('ambiguous run directory: ' + repr(destinations))
        dest = destinations.pop()
        log = dest / 'stdout.log'
        observed = dest / 'input.cir.prn'
        expected = Path(ref['directory']) / 'input.cir.prn'
        fields = campaign.parse_log(log.read_text())
        pointwise = campaign.compare_waveforms(expected, observed)
        run = args.output / str(index)
        run.mkdir()
        names = expected.read_text().splitlines()[0].split()[2:]
        text = re.sub(r'^\s*\*comp[^\n]*\n', '', (dest / 'input.cir').read_text(), flags=re.I | re.M)
        comments = ''.join(f'*comp {name} reltol=1e-8 abstol=1e-10 zerotol=0 absdifftol=0 numfail=0\n' for name in names)
        text, count = re.subn(r'^(\s*\.end\s*)$', lambda m: comments + m[0], text, flags=re.I | re.M)
        if count != 1:
            raise ValueError('expected one .end')
        netlist = run / 'verification.cir'
        netlist.write_text(text)
        command = ['perl', str(verifier), '--verbose', str(netlist), str(expected), str(observed)]
        with (run / 'verification.log').open('x') as stream:
            try:
                code = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, timeout=60).returncode
            except subprocess.TimeoutExpired:
                code = 'timeout'
        values = fields['linear_seconds']
        # Missing exit/timeout fields in old diagnostic ledgers cannot prove
        # full qualification. Verify waveforms without inventing that evidence.
        execution_proven = row.get('exit_code') == 0 and row.get('timed_out') is False
        waveform_verified = bool(code == 0 and
            all(fields[k] == ref[k] for k in ('successful_steps', 'failed_steps')) and
            pointwise['status'] not in ('row_count_mismatch', 'column_count_mismatch', 'nonfinite'))
        qualified = bool(execution_proven and row['closure']['all_session_processes_absent']
                         and fields['end_of_simulation'] and not fields['linear_failure_messages']
                         and values and all(math.isfinite(v) and v > 0 for v in values)
                         and all(fields[k] == ref[k] for k in ('successful_steps', 'failed_steps'))
                         and pointwise['status'] not in ('row_count_mismatch', 'column_count_mismatch', 'nonfinite')
                         and code == 0)
        result = dict(index=index, case=case, domain=domain, threads=row.get('threads', 8),
                      trial=row['trial'], side=row.get('side', str(row['threads'])), directory=str(dest),
                      legacy_valid=row.get('valid'), qualified=qualified, pointwise=pointwise,
                      waveform_verified=waveform_verified, execution_proven=execution_proven,
                      upstream_verifier=code, command=command,
                      pins={str(p): sha(p) for p in (log, observed, expected, netlist, run / 'verification.log')})
        results.append(result)
        with (args.output / 'observations.jsonl').open('a') as stream:
            stream.write(json.dumps(result) + '\n')
    if not results:
        raise ValueError('empty diagnostic ledger')
    summary = dict(diagnostic_only=True, acceptance_evidence=False,
                   all_qualified=all(r['qualified'] for r in results),
                   all_waveforms_verified=all(r['waveform_verified'] for r in results), observations=len(results),
                   corrected_grid_labels=sum(r['legacy_valid'] is False and r['waveform_verified'] for r in results),
                   pins={str(p): sha(p) for p in (original, verifier, runner, script_archive,
                         *(root / 'observations.jsonl' for root in args.reference_campaign))},
                   unqualified=[r['index'] for r in results if not r['qualified']])
    (args.output / 'audit.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))
    return 0 if summary['all_waveforms_verified'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
