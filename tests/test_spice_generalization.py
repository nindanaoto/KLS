#!/usr/bin/env python3
"""Scoring rules of run_spice_generalization.py on synthetic observations."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import run_spice_generalization as g  # noqa: E402


def rows_for(times, passes=2, bad=()):
    rows = []
    for (case, threads, side), seconds in times.items():
        for trial in range(passes):
            rows.append(dict(domain='llc96', case=case, threads=threads, trial=trial, side=side,
                             linear=seconds, per_solve=seconds, wall_seconds=seconds + 1,
                             qualified=(case, threads, side) not in bad))
    return rows


def metadata(cases, name='dev', peers=('ck', 'st'), baseline=False):
    return dict(set=name, sealed=name == 'holdout', cases=cases, domains=[('llc96', '0-7')],
                threads=[1, 8], passes=2, peers=list(peers), baseline=baseline)


def cell(detail, case, threads, peer):
    return next(c for c in detail['cells'] if (c['case'], c['threads'], c['peer']) == (case, threads, peer))


def test_peer_is_scored_at_its_best_configuration():
    times = {('a', 1, 'kls'): 1.0, ('a', 8, 'kls'): 0.5,
             ('a', 1, 'ck'): 0.8, ('a', 1, 'st'): 2.0,
             # CKTSO collapses at eight threads unless auto-threads is on.
             ('a', 8, 'ck'): 4.0, ('a', 8, 'ck-auto'): 0.9, ('a', 8, 'ck-busy'): 3.0,
             ('a', 8, 'st'): 0.6, ('a', 8, 'st-active'): 0.4}
    summary, detail = g.reduce_results(metadata(['a']), rows_for(times))
    ck8 = cell(detail, 'a', 8, 'ck')
    # The peer's own one-thread run (0.8 s) beats all of its 8-thread settings.
    assert ck8['peer_best'] == 'ck@t1' and abs(ck8['ratio'] - 0.5 / 0.8) < 1e-12
    st8 = cell(detail, 'a', 8, 'st')
    assert st8['peer_best'] == 'st-active@t8' and abs(st8['ratio'] - 0.5 / 0.4) < 1e-12
    assert abs(cell(detail, 'a', 1, 'ck')['ratio'] - 1.25) < 1e-12
    assert summary['groups']['st-t8']['losses'] == 1 and not summary['gate']
    assert summary['thread_monotone']


def test_failed_peer_configuration_is_dropped_not_fatal():
    times = {('a', 1, 'kls'): 1.0, ('a', 8, 'kls'): 1.0, ('a', 1, 'ck'): 2.0, ('a', 1, 'st'): 2.0,
             ('a', 8, 'ck'): 0.1, ('a', 8, 'ck-auto'): 2.0, ('a', 8, 'ck-busy'): 2.0,
             ('a', 8, 'st'): 2.0, ('a', 8, 'st-active'): 2.0}
    summary, detail = g.reduce_results(metadata(['a']), rows_for(times, bad={('a', 8, 'ck')}))
    assert 'ck@t8' not in cell(detail, 'a', 8, 'ck')['configurations']
    assert cell(detail, 'a', 8, 'ck')['ratio'] == 0.5 and summary['complete'] and summary['gate']


def test_kls_failure_and_thread_regression_are_reported():
    times = {('a', 1, 'kls'): 1.0, ('a', 8, 'kls'): 1.2, ('a', 1, 'ck'): 3.0, ('a', 1, 'st'): 3.0,
             ('a', 8, 'ck'): 3.0, ('a', 8, 'ck-auto'): 3.0, ('a', 8, 'ck-busy'): 3.0,
             ('a', 8, 'st'): 3.0, ('a', 8, 'st-active'): 3.0,
             ('b', 1, 'kls'): 1.0, ('b', 8, 'kls'): 1.0}
    summary, _ = g.reduce_results(metadata(['a', 'b']), rows_for(times, bad={('b', 8, 'kls')}))
    assert not summary['thread_monotone']          # 8 threads 20% slower than one
    assert summary['kls_failures'] == 1 and not summary['complete'] and not summary['gate']


def test_tune_gate_needs_every_cell_and_holdout_needs_geomean_and_cap():
    fast = {('a', 1, 'kls'): 1.0, ('a', 8, 'kls'): 1.0}
    fast.update({('a', t, s): 1.05 for t in (1, 8) for s in g.peer_sides(t, ['ck', 'st'])})
    assert g.reduce_results(metadata(['a'], 'tune'), rows_for(fast))[0]['gate']
    marginal = dict(fast)
    marginal[('a', 1, 'ck')] = 1.01       # KLS only 1% ahead of CKTSO at one thread
    assert not g.reduce_results(metadata(['a'], 'tune'), rows_for(marginal))[0]['gate']
    summary, _ = g.reduce_results(metadata(['a'], 'holdout'), rows_for(marginal))
    assert summary['gate'] and summary['sealed'] and 'cells' not in summary
    outlier = dict(fast)
    outlier[('a', 8, 'st')] = 0.8          # KLS 25% slower than one peer on one deck
    assert not g.reduce_results(metadata(['a'], 'holdout'), rows_for(outlier))[0]['gate']


def test_gate_applies_to_application_managed_mode_and_reports_certified():
    times = {('a', 1, 'kls'): 1.2, ('a', 8, 'kls'): 1.2, ('a', 1, 'kls-app'): 1.0, ('a', 8, 'kls-app'): 1.0}
    times.update({('a', t, s): 1.1 for t in (1, 8) for s in g.peer_sides(t, ['ck', 'st'])})
    summary, detail = g.reduce_results(metadata(['a'], 'tune'), rows_for(times))
    assert summary['gated_mode'] == 'kls-app' and summary['gate']
    assert summary['certified_overall']['losses'] == 4 and summary['overall']['losses'] == 0
    assert {c['mode'] for c in detail['cells']} == {'kls', 'kls-app'}


def test_held_out_decks_with_differing_steps_are_scored_per_solve():
    times = {('a', t, s): 1.0 for t in (1, 8) for s in g.sides_for(t, ['ck', 'st'], False)}
    rows = rows_for(times, bad={('a', 1, 'kls-app'), ('a', 1, 'ck')})
    for r in rows:
        if not r['qualified']:
            r['tolerant'] = True
        if r['side'] == 'kls-app' and r['threads'] == 1:
            r.update(linear=3.0, per_solve=0.5)       # more solves, each one faster
    summary, detail = g.reduce_results(metadata(['a']), rows)
    assert summary['complete'] and summary['time_metric'] == 'per_solve' and summary['tolerant_runs'] == 4
    app = next(c for c in detail['cells'] if (c['mode'], c['threads'], c['peer']) == ('kls-app', 1, 'ck'))
    assert abs(app['ratio'] - 0.5) < 1e-12
    # The tuning set keeps the strict rule: a tolerant run there is a failure.
    summary, _ = g.reduce_results(metadata(['a'], 'tune'), rows)
    assert not summary['complete'] and summary['kls_failures'] == 1


def test_trajectory_error_is_zero_for_identical_runs_and_scaled_by_tolerance(tmp=None):
    import tempfile
    with tempfile.TemporaryDirectory() as directory:
        a, b, c = (Path(directory) / n for n in 'abc')
        a.write_text('Index TIME V(1)\n' + ''.join(f'{i} {i * 0.5} {float(i)}\n' for i in range(5)))
        b.write_text('Index TIME V(1)\n' + ''.join(f'{i} {i * 1.0} {i * 2.0}\n' for i in range(3)))
        c.write_text('Index TIME V(1)\n' + ''.join(f'{i} {i * 0.5} {i + 0.004}\n' for i in range(5)))
        assert g.trajectory_error(a, a) == 0
        assert g.trajectory_error(a, b) < 1e-9          # same line, coarser steps
        assert abs(g.trajectory_error(a, c) - 0.004 / (4e-3 + 1e-6)) < 1e-9
        b.write_text('Index TIME V(1)\n0 0 0\n1 1 1\n')   # stops early
        assert g.trajectory_error(a, b) is None


def test_sides_and_rotation():
    assert g.sides_for(1, ['ck', 'st'], False) == ['kls', 'kls-app', 'ck', 'st']
    assert g.sides_for(8, ['ck', 'st'], True) == ['kls', 'kls-app', 'baseline', 'ck', 'ck-auto', 'ck-busy', 'st', 'st-active']
    sides = g.sides_for(8, ['ck', 'st'], False)
    firsts = {g.rotate(sides, 0, trial, False)[0] for trial in range(len(sides))}
    assert firsts == set(sides)            # every side leads once per cycle


if __name__ == '__main__':
    for name, function in sorted(globals().items()):
        if name.startswith('test_'):
            function()
    print('ok')
