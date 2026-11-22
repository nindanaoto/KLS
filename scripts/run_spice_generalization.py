#!/usr/bin/env python3
"""Peer-best SPICE comparison on the tuning, development or held-out decks.

This complements run_spice_campaign.py, whose gates select on the seven tuning
circuits.  Differences:

  * The deck list comes from a set: `tune` (the seven campaign circuits), or
    `dev` / `holdout` from a manifest built by build_spice_holdout_manifest.py.
  * KLS runs library defaults and a thread count, in two accuracy modes:
    `kls` certifies every solve (the library default) and `kls-app` uses the
    application-managed mode, which is how the peers run under Xyce's Newton
    loop.  The win gate applies to `kls-app`; `kls` is always reported.
  * Each peer is scored at its best.  For a KLS run at T threads the peer
    score is the most favourable (for the peer) of all its documented
    configurations at T threads and at one thread, so a peer is never charged
    for a thread count or wait policy that makes it slower.
  * `holdout` results are sealed: per-deck numbers go to cells.sealed.json and
    only aggregates are printed or written to summary.json.

Ratios are medians over passes of paired KLS/peer linear-solution times; a
ratio below one means KLS is faster.  Correctness of every run is checked
against a fresh KLU reference (identical step counts and the upstream waveform
verifier), exactly as in run_spice_campaign.py.

On held-out sets a deck whose step sequence depends on rounding would fail that
check for every solver, peers included.  There a run also counts when its waveform
stays within TRAJECTORY_TOLERANCE integration-tolerance units of KLU's, and all
held-out ratios use time per linear solve, which equals the total-time ratio
whenever the step counts agree.  The tuning set keeps the strict rule.
"""
import argparse
import fcntl
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import signal
import statistics
import subprocess
import time

from run_spice_campaign import (CASES, DOMAINS, compare_waveforms, parse_log, sha, stage,
                                wait_child, waveform, write_json)

# side -> (Xyce solver type, extra `.options linsol` text, environment overrides)
KLS_SIDES = {'kls': ('KLS', '', {}), 'baseline': ('KLS', '', {}),
             'kls-app': ('KLS', 'KLS_ACCURACY=APPLICATION', {})}
GATED_MODE = 'kls-app'
PEER_CONFIGS = {
    'ck': {
        'ck': ('CKTSO', '', {}),
        'ck-auto': ('CKTSO', 'CKTSO_AUTO_THREADS=1', {}),
        'ck-busy': ('CKTSO', 'CKTSO_THREAD_SYNC=-1', {}),
    },
    'st': {
        'st': ('SUBTREELU', '', {}),
        'st-active': ('SUBTREELU', '', {'OMP_WAIT_POLICY': 'active'}),
    },
}
# Configurations that only differ from the base one when several threads run.
MULTITHREAD_ONLY = {'ck-auto', 'ck-busy', 'st-active'}
WIN = 0.98
HELD_OUT_WORST = 1.10
THREAD_MONOTONE = 1.03


def peer_sides(threads, peers):
    return [side for peer in peers for side in PEER_CONFIGS[peer]
            if threads > 1 or side not in MULTITHREAD_ONLY]


def sides_for(threads, peers, baseline):
    return ['kls', 'kls-app'] + (['baseline'] if baseline else []) + peer_sides(threads, peers)


def rotate(sides, case_index, trial, confirmation):
    shift = (case_index + trial + int(not confirmation)) % len(sides)
    return sides[shift:] + sides[:shift]


def geomean(values):
    return math.exp(statistics.mean(math.log(v) for v in values))


# Held-out decks whose step sequence depends on rounding are still scored: a run that
# stays within this many integration-tolerance units of the KLU waveform counts, and
# such decks are timed per linear solve because the solvers do different numbers.
TRAJECTORY_TOLERANCE = 10.0


def trajectory_error(reference, output):
    """Largest per-signal RMS distance between two runs that took different steps.

    The candidate is interpolated onto the reference time points and the distance is
    expressed in units of Xyce's default integration tolerance for that signal
    (reltol 1e-3 of its largest magnitude plus abstol 1e-6).  Returns None when the
    runs do not cover the same time span or the same signals.
    """
    a, b = waveform(reference), waveform(output)
    if not a or not b or len(a[0]) < 2 or any(len(r) != len(a[0]) for r in a + b):
        return None
    span = a[-1][0] - a[0][0]
    if not span > 0 or abs(a[-1][0] - b[-1][0]) > 1e-9 * span or abs(a[0][0] - b[0][0]) > 1e-9 * span:
        return None
    times = [r[0] for r in b]
    if any(y < x for x, y in zip(times, times[1:])):
        return None
    worst = 0.0
    for column in range(1, len(a[0])):
        total, j = 0.0, 0
        for row in a:
            t = min(max(row[0], times[0]), times[-1])
            while j + 2 < len(times) and times[j + 1] < t:
                j += 1
            t0, t1 = times[j], times[min(j + 1, len(times) - 1)]
            y0, y1 = b[j][column], b[min(j + 1, len(b) - 1)][column]
            y = y0 if t1 <= t0 else y0 + (y1 - y0) * min(max((t - t0) / (t1 - t0), 0.0), 1.0)
            total += (y - row[column]) ** 2
        scale = 1e-3 * max(abs(r[column]) for r in a) + 1e-6
        worst = max(worst, math.sqrt(total / len(a)) / scale)
    return worst if math.isfinite(worst) else None


def reduce_results(metadata, rows):
    index = {(r['domain'], r['case'], r['threads'], r['trial'], r['side']): r for r in rows}  # last attempt wins
    passes, peers = metadata['passes'], metadata['peers']
    held_out = metadata['set'] != 'tune'
    time_metric = 'per_solve' if held_out else 'linear'

    def usable(row):
        return bool(row.get('qualified') or (held_out and row.get('tolerant')))

    def paired(domain, case, threads, side, other_threads, other_side, metric):
        values = []
        for trial in range(passes):
            a = index.get((domain, case, threads, trial, side))
            b = index.get((domain, case, other_threads, trial, other_side))
            if a and b and usable(a) and usable(b) and a.get(metric) and b.get(metric):
                values.append(a[metric] / b[metric])
        return statistics.median(values) if len(values) == passes else None

    cells, kls_failures, scaling, baseline = [], [], [], []
    modes = [m for m in ('kls', 'kls-app') if any(r['side'] == m for r in rows)]
    # A circuit on which no solver under test (KLS, baseline or any peer configuration)
    # reproduced the reference in any run measures the circuit, not a solver: it is
    # reported and left out of the score for everyone instead of failing KLS alone.
    measured = {(r['domain'], r['case']) for r in rows if r['side'] != 'klu'}
    verified = {(r['domain'], r['case']) for r in rows if r['side'] != 'klu' and usable(r)}
    has_peer = {(r['domain'], r['case']) for r in rows if r['side'] not in KLS_SIDES and r['side'] != 'klu'}
    unverifiable = sorted((measured & has_peer) - verified)
    for domain, _ in metadata['domains']:
        for case in metadata['cases']:
            if (domain, case) in unverifiable:
                continue
            for threads in metadata['threads']:
                for mode in modes:
                    if any(not usable(index.get((domain, case, threads, t, mode), {}))
                           for t in range(passes)):
                        kls_failures.append(dict(domain=domain, case=case, threads=threads, mode=mode))
                        continue
                    if threads > 1 and 1 in metadata['threads']:
                        ratio = paired(domain, case, threads, mode, 1, mode, time_metric)
                        if ratio is not None:
                            scaling.append(dict(domain=domain, case=case, threads=threads, mode=mode, ratio=ratio))
                    if metadata['baseline'] and mode == 'kls':
                        ratio = paired(domain, case, threads, 'kls', threads, 'baseline', time_metric)
                        wall = paired(domain, case, threads, 'kls', threads, 'baseline', 'wall_seconds')
                        baseline.append(dict(domain=domain, case=case, threads=threads, ratio=ratio, wall=wall))
                    for peer in peers:
                        candidates = {}
                        for other in sorted({threads, 1} & set(metadata['threads'])):
                            for side in peer_sides(other, [peer]):
                                ratio = paired(domain, case, threads, mode, other, side, time_metric)
                                if ratio is not None:
                                    candidates[f'{side}@t{other}'] = ratio
                        cell = dict(domain=domain, case=case, threads=threads, peer=peer, mode=mode,
                                    configurations=candidates)
                        if candidates:
                            best = max(candidates, key=candidates.get)
                            cell.update(ratio=candidates[best], peer_best=best)
                        else:
                            cell.update(ratio=None, peer_best=None)
                        cells.append(cell)
    gated_mode = GATED_MODE if GATED_MODE in modes else 'kls'

    def aggregate(selection):
        values = sorted(c['ratio'] for c in selection if c['ratio'] is not None)
        if not values:
            return dict(count=0)
        return dict(count=len(values), geomean=geomean(values), worst=values[-1],
                    p90=values[min(len(values) - 1, math.ceil(0.9 * len(values)) - 1)],
                    wins=sum(v <= WIN for v in values), losses=sum(v > 1 for v in values))

    def groups_for(mode):
        return {f'{peer}-t{threads}': aggregate([c for c in cells if c['mode'] == mode and c['peer'] == peer
                                                 and c['threads'] == threads])
                for peer in peers for threads in metadata['threads']}

    groups = groups_for(gated_mode)
    scored = [c for c in cells if c['mode'] == gated_mode and c['ratio'] is not None]
    unscored = [dict(domain=c['domain'], case=c['case'], threads=c['threads'], peer=c['peer'], mode=c['mode'])
                for c in cells if c['ratio'] is None]
    monotone = all(s['ratio'] <= THREAD_MONOTONE for s in scaling)
    complete = not kls_failures and bool(scored)
    if metadata['set'] == 'tune':
        gate = complete and all(c['ratio'] <= WIN for c in scored)
    else:
        gate = complete and all(g.get('count', 0) and g['geomean'] < 1 and g['worst'] <= HELD_OUT_WORST
                                for g in groups.values())
    summary = dict(set=metadata['set'], sealed=metadata['sealed'], peer_policy='peer-best',
                   cases=len(metadata['cases']), passes=passes, complete=complete,
                   gated_mode=gated_mode, gate=gate, thread_monotone=monotone,
                   groups=groups, overall=aggregate(scored),
                   certified_groups=groups_for('kls') if gated_mode != 'kls' else None,
                   certified_overall=aggregate([c for c in cells if c['mode'] == 'kls']) if gated_mode != 'kls' else None,
                   kls_failures=len(kls_failures), time_metric=time_metric,
                   unverifiable_cases=len(unverifiable),
                   tolerant_runs=sum(bool(r.get('tolerant')) and not r.get('qualified') for r in index.values()),
                   peer_cells_without_valid_configuration=len(unscored),
                   thread_scaling=aggregate(scaling),
                   peer_hangs=sum(bool(r.get('timed_out')) and r['side'] not in KLS_SIDES for r in rows),
                   baseline=aggregate([b for b in baseline]) if baseline else None)
    detail = dict(cells=cells, kls_failures=kls_failures, unscored_peer_cells=unscored,
                  unverifiable_cases=[dict(domain=d, case=c) for d, c in unverifiable],
                  thread_scaling=scaling, baseline=baseline)
    return summary, detail


def run(args):
    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location('cleanup', root / 'bench/process_session_cleanup.py')
    cleanup = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cleanup)
    staging = {}
    case_root = args.case_root or args.regression / 'Netlists'
    if args.set == 'tune':
        cases = list(CASES)
    else:
        manifest = json.loads(args.manifest.read_text())
        cases = [c['case'] for c in manifest['cases'] if c['split'] == args.set]
        staging = {c['case']: dict(companions=bool(c.get('companions')), tran_scale=c.get('tran_scale', 1))
                   for c in manifest['cases']}
    if args.case:
        if args.set == 'holdout':
            raise SystemExit('held-out decks are never run individually')
        cases = [c for c in cases if any(pattern in c for pattern in args.case)]
    if not cases:
        raise SystemExit('no decks selected')
    if not args.diagnostic:
        busy = []
        for process in Path('/proc').iterdir():
            if process.name.isdigit() and int(process.name) != os.getpid():
                try:
                    name = (process / 'comm').read_text().strip()
                except (FileNotFoundError, ProcessLookupError, PermissionError):
                    continue
                if name in ('cc1', 'cc1plus', 'f951', 'ctest', 'cmake', 'Xyce', 'kls_bench'):
                    busy.append((int(process.name), name))
        if busy:
            raise RuntimeError(f'timing requires idle build/test/solver processes: {busy}')
    output = args.output.resolve()
    output.mkdir(parents=True)
    lock = open(args.workspace / '.spice-timing.lock', 'w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    binaries = {'kls': args.candidate.resolve()}
    if args.baseline:
        binaries['baseline'] = args.baseline.resolve()
    frozen = output / 'binaries'
    frozen.mkdir()
    hashes = {}
    for name, original in list(binaries.items()):
        hashes[str(original)] = sha(original)
        target = frozen / name
        target.write_bytes(original.read_bytes())
        target.chmod(0o755)
        binaries[name] = target
    inputs = [Path(__file__), Path(__file__).with_name('run_spice_campaign.py'), args.image, args.license,
              args.regression / 'TestScripts/xyce_verify.pl', *args.peer_library,
              *[case_root / c for c in cases]]
    if args.manifest:
        inputs.append(args.manifest)
    hashes.update({str(p.resolve()): sha(p) for p in inputs})
    domains = [d for d in DOMAINS if d[0] in args.domain]
    if args.order == 'confirmation':
        domains.reverse()
    else:
        cases.reverse()
    peers = [] if args.no_peers else list(PEER_CONFIGS)
    metadata = dict(set=args.set, sealed=args.set == 'holdout', cases=cases, domains=domains,
                    threads=args.threads, passes=args.passes, peers=peers, baseline=bool(args.baseline),
                    order=args.order, diagnostic=args.diagnostic, timeout=args.timeout,
                    peer_configurations={p: {s: list(v[:2]) + [v[2]] for s, v in c.items()}
                                         for p, c in PEER_CONFIGS.items()},
                    input_sha256=hashes, cpu_topology=subprocess.check_output(['lscpu'], text=True),
                    primary_metric='final Total Linear Solution Time')
    write_json(output / 'metadata.json', metadata)
    base_env = {k: v for k, v in os.environ.items() if not k.startswith(
        ('LD_', 'KLS_', 'OMP_', 'GOMP_', 'OPENBLAS_', 'APPTAINERENV_', 'ASAN_', 'UBSAN_', 'TSAN_'))}
    base_env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', XYCE_NO_TRACKING='1')
    rows = []
    quiet = metadata['sealed']

    def execute(domain, cpus, case, count, trial, side, reference=None, attempt=0):
        # A hung solver must not hold the machine for the full timeout: after the KLU
        # reference has run, every other run gets a bound scaled from its wall time.
        limit = args.timeout if reference is None else min(
            args.timeout, max(args.min_timeout, args.timeout_factor * reference['wall_seconds']))
        label = hashes[str((case_root / case).resolve())][:12] if quiet else Path(case).stem
        dest = output / domain / label / (f't{count}-p{trial}-{side}' + (f'-retry{attempt}' if attempt else ''))
        if side == 'klu':
            solver, extra, overrides = 'KLU', '', {}
        elif side in KLS_SIDES:
            solver, extra, overrides = KLS_SIDES[side]
        else:
            solver, extra, overrides = next(c[side] for c in PEER_CONFIGS.values() if side in c)
        stage(case_root / case, dest, solver, count, args.license, extra,
              **staging.get(case, {}))
        executable = binaries['baseline' if side == 'baseline' and 'baseline' in binaries else 'kls']
        cmd = ['taskset', '-c', cpus, 'apptainer', 'exec', '--bind', f'{args.workspace}:/work',
               '--bind', f'{args.workspace}:{args.workspace}', '--bind', '/tmp:/tmp',
               '--pwd', str(dest), str(args.image), str(executable), 'input.cir']
        env = dict(base_env, **overrides)
        start = time.monotonic()
        with (dest / 'stdout.log').open('x') as stream:
            child = subprocess.Popen(cmd, stdout=stream, stderr=subprocess.STDOUT, env=env, start_new_session=True)
            try:
                code, timeout = wait_child(child, limit)
            except BaseException:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
                cleanup.drain_exited_apptainer_session(child.pid, args.image)
                raise
        wall = time.monotonic() - start
        closure = cleanup.drain_exited_apptainer_session(child.pid, args.image)
        row = dict(domain=domain, case=case, threads=count, trial=trial, side=side, directory=str(dest),
                   exit_code=code, timed_out=timeout, wall_seconds=wall, attempt=attempt)
        log = parse_log((dest / 'stdout.log').read_text(errors='replace'))
        values = waveform(dest / 'input.cir.prn')
        row['qualified'] = bool(
            code == 0 and not timeout and closure['all_session_processes_absent'] and log['end_of_simulation'] and
            not log['linear_failure_messages'] and log['successful_steps'] and log['failed_steps'] and
            log['linear_seconds'] and all(math.isfinite(v) and v > 0 for v in log['linear_seconds']) and
            values and all(math.isfinite(v) for r in values for v in r))
        row['linear'] = log['linear_seconds'][-1] if log['linear_seconds'] else None
        row['steps'] = [log['successful_steps'][-1] if log['successful_steps'] else None,
                        log['failed_steps'][-1] if log['failed_steps'] else None]
        solves = re.findall(r'Number Linear Solves:\s*(\d+)', (dest / 'stdout.log').read_text(errors='replace'))
        row['solves'] = int(solves[-1]) if solves else None
        row['per_solve'] = row['linear'] / row['solves'] if row['linear'] and row['solves'] else None
        row['tolerant'] = False
        if reference is not None and row['qualified']:
            ref_output, test_output = Path(reference['directory']) / 'input.cir.prn', dest / 'input.cir.prn'
            row['waveform'] = compare_waveforms(ref_output, test_output)['status']
            names = ref_output.read_text().splitlines()[0].split()[2:]
            netlist = re.sub(r'^\s*\*comp[^\n]*\n', '', (dest / 'input.cir').read_text(), flags=re.I | re.M)
            comments = ''.join(f'*comp {name} reltol=1e-8 abstol=1e-10 zerotol=0 absdifftol=0 numfail=0\n'
                               for name in names)
            netlist = re.sub(r'^(\s*\.end\s*)$', lambda m: comments + m[0], netlist, flags=re.I | re.M)
            (dest / 'verification.cir').write_text(netlist)
            try:
                with (dest / 'verification.log').open('x') as stream:
                    row['verifier'] = subprocess.run(
                        ['perl', str(args.regression / 'TestScripts/xyce_verify.pl'), '--verbose',
                         str(dest / 'verification.cir'), str(ref_output), str(test_output)],
                        stdout=stream, stderr=subprocess.STDOUT, timeout=120).returncode
            except subprocess.TimeoutExpired:
                row['verifier'] = 'timeout'
            row['qualified'] = bool(reference['qualified'] and row['steps'] == reference['steps'] and
                                    row['verifier'] == 0)
            if not row['qualified'] and reference['qualified'] and metadata['set'] != 'tune':
                row['trajectory_error'] = trajectory_error(ref_output, test_output)
                row['tolerant'] = (row['trajectory_error'] is not None and
                                   row['trajectory_error'] <= TRAJECTORY_TOLERANCE)
        elif reference is not None:
            row['qualified'] = False
        rows.append(row)
        with (output / 'observations.jsonl').open('a') as stream:
            stream.write(json.dumps(row, allow_nan=False) + '\n')
        if quiet:
            print('.', end='', flush=True)
        else:
            print(domain, Path(case).stem, count, trial, side, 'PASS' if row['qualified'] else 'TOLERANT' if row['tolerant'] else 'FAIL',
                  f"{row['linear']}", flush=True)
        return row

    for domain, cpus in domains:
        for index, case in enumerate(cases):
            reference = execute(domain, cpus, case, 1, -1, 'klu')
            for count in args.threads:
                sides = sides_for(count, peers, bool(args.baseline))
                for trial in range(args.passes):
                    for side in rotate(sides, index, trial, args.order == 'confirmation'):
                        row = execute(domain, cpus, case, count, trial, side, reference)
                        # Peers get retries after a hang so their flakiness is recorded, not
                        # scored as a KLS win; KLS itself never gets a second attempt.
                        attempt = 0
                        while row['timed_out'] and side not in KLS_SIDES and attempt < args.peer_retries:
                            attempt += 1
                            row = execute(domain, cpus, case, count, trial, side, reference, attempt)
    changed = [p for p, h in hashes.items() if not Path(p).exists() or sha(p) != h]
    summary, detail = reduce_results(metadata, rows)
    summary['changed_inputs'] = changed
    if changed or args.diagnostic:
        summary['gate'] = False
    write_json(output / 'summary.json', summary)
    write_json(output / ('cells.sealed.json' if metadata['sealed'] else 'cells.json'), detail)
    print()
    print(json.dumps(summary, indent=2))
    return 0 if summary['complete'] else 1


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for name in ('workspace', 'regression', 'image', 'license', 'candidate', 'output'):
        p.add_argument('--' + name, type=Path, required=True)
    p.add_argument('--baseline', type=Path, help='Optional second KLS Xyce binary for an A/B ratio')
    p.add_argument('--manifest', type=Path, help='Held-out manifest (required for dev and holdout)')
    p.add_argument('--set', choices=['tune', 'dev', 'holdout'], required=True)
    p.add_argument('--case', action='append', default=[], help='Substring filter; not allowed for holdout')
    p.add_argument('--case-root', type=Path,
                   help='Directory the manifest cases are relative to (default: the regression Netlists)')
    p.add_argument('--peer-library', type=Path, action='append', default=[])
    p.add_argument('--no-peers', action='store_true')
    p.add_argument('--passes', type=int, default=4)
    p.add_argument('--threads', type=lambda s: [int(x) for x in s.split(',')], default=[1, 8])
    p.add_argument('--domain', action='append', choices=[d[0] for d in DOMAINS], default=[])
    p.add_argument('--order', choices=['initial', 'confirmation'], default='initial')
    p.add_argument('--timeout', type=float, default=1800)
    p.add_argument('--timeout-factor', type=float, default=30,
                   help='Per-run bound as a multiple of the KLU reference wall time')
    p.add_argument('--min-timeout', type=float, default=60)
    p.add_argument('--peer-retries', type=int, default=3, help='Retries of a peer run that hung')
    p.add_argument('--diagnostic', action='store_true', help='Allow a busy machine; the gate cannot pass')
    args = p.parse_args()
    if args.set != 'tune' and not args.manifest:
        p.error('--manifest is required for dev and holdout')
    if not args.domain:
        args.domain = [d[0] for d in DOMAINS]
    if args.passes < 1:
        p.error('--passes must be positive')
    return run(args)


if __name__ == '__main__':
    raise SystemExit(main())
