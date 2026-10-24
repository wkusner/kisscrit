"""Tests of scripts/measure/sample_full_size.py (the sampled size of the `full` run).

Status: untrusted (tests of measurement tooling; no part of the proof).

What is checked, and against what:
  * the constants of the script against the repository files they were taken from (gen files, the two pilot
    summaries, the checker ratio in docs/decisions.md), so a stale constant fails here;
  * the estimator against a toy population whose total is known: it is unbiased over all samples and its standard
    error formula equals the true variance of the estimator (exact enumeration of every sample);
  * the tiling counter against a direct count of `tilings` on a real triangulation, and the code hashes against
    `candidates_of_triangulation` itself;
  * the wall-clock cap with fake enumerators (one that yields, one that never yields), including that the alarm
    and its handler are restored;
  * the calibration against the gen files: synthetic mutations must be detected, and the first seconds of the real
    enumeration of indices 1, 5 and 7 must contain only codes the gen files list;
  * the whole driver with a fake measurement (fast): sample indices, the overlap, exit statuses 2 and 3, the closing
    sentence, and that no output carries a host name, user name or path.
"""
import importlib.util
import itertools
import json
import os
import pathlib
import re
import signal
import socket
import getpass

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
# The word reserved for claims backed by a certificate and an accepted checker run (docs/decisions.md D-034); measurement
# output must not use it.  This is the only place the word is written out in this file.
RESERVED_WORD = 'certified'
_spec = importlib.util.spec_from_file_location('sample_full_size', ROOT / 'scripts' / 'measure' / 'sample_full_size.py')
sfs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sfs)


@pytest.fixture(scope='module')
def gen():
    """({index: hashes of the codes the gen files attribute to it}, {index: done-marker count})."""
    return sfs.load_gen_cond()


@pytest.fixture(scope='module')
def triangulations():
    return sfs.load_triangulations()


# ------------------------------------------------------------------ constants against the repository

def test_constants_match_the_repository_files(gen):
    by_tri, markers = gen
    assert len(set().union(*by_tri.values())) == sfs.COND_CODES == 45260
    assert markers == {0: 9689, 1: 38523, 5: 6596, 7: 36028}
    assert {i: len(s) for i, s in by_tri.items()} == markers
    # How often a cond code occurs among the triangulations (STATUS.md, "computed"): 1x 9,140, 2x 26,664, 3x 9,456.
    multiplicity = {}
    for s in by_tri.values():
        for h in s:
            multiplicity[h] = multiplicity.get(h, 0) + 1
    histogram = {m: sum(1 for v in multiplicity.values() if v == m) for m in (1, 2, 3)}
    assert histogram == {1: 9140, 2: 26664, 3: 9456}
    desktop = json.loads((ROOT / 'results' / 'completeness' / 'summary-cond.json').read_text())
    server = json.loads((ROOT / 'results' / 'pilot' / 'cond-2026-10-22' / 'summary-cond.json').read_text())
    assert desktop['generated'] == server['generated'] == sfs.PILOT_CODES
    assert round(desktop['cpu_hours'], 2) == sfs.DESKTOP_CPU_H
    assert round(server['cpu_hours'], 2) == sfs.SERVER_CPU_H


def test_checker_ratio_is_the_one_in_decisions_md():
    text = (ROOT / 'docs' / 'decisions.md').read_text()
    line = [l for l in text.splitlines() if l.startswith('- %s' % sfs.CHECKER_RATIO_DATE) and 'F4 performance gate' in l]
    assert len(line) == 1
    assert 'measured ratio %s' % sfs.CHECKER_RATIO[0] in line[0] and 'about %s' % sfs.CHECKER_RATIO[1] in line[0]


# ------------------------------------------------------------------ sampling and estimators

def test_draw_sample_is_seeded_sorted_and_avoids_cond():
    population = [i for i in range(7595) if i not in sfs.COND]
    a, b = sfs.draw_sample(5, 150, population), sfs.draw_sample(5, 150, population)
    assert a == b == sorted(set(a)) and len(a) == 150 and not set(a) & set(sfs.COND)
    assert sfs.draw_sample(6, 150, population) != a
    assert sfs.draw_sample(1, len(population), population) == population
    for bad in (0, len(population) + 1):
        with pytest.raises(ValueError):
            sfs.draw_sample(1, bad, population)


def test_estimator_is_unbiased_and_its_standard_error_is_exact_on_a_toy_population():
    population = [1, 2, 3, 4, 5, 6, 20]
    M, n = len(population), 3
    estimates = [sfs.fpc_total(list(s), M) for s in itertools.combinations(population, n)]
    totals = [e['total'] for e in estimates]
    assert sum(totals) / len(totals) == pytest.approx(sum(population))
    true_variance = sum((t - sum(population)) ** 2 for t in totals) / len(totals)
    assert sum(e['se'] ** 2 for e in estimates) / len(estimates) == pytest.approx(true_variance)
    # A worked sample: values 1, 3, 5 of a population of 6: mean 3, total 18, s^2 = 4, se = 6 sqrt(0.5 x 4 / 3).
    e = sfs.fpc_total([1, 3, 5], 6)
    assert e['total'] == 18 and e['se'] == pytest.approx(6 * (0.5 * 4 / 3) ** 0.5)
    assert sfs.fpc_total([7], 10)['se'] is None
    assert sfs.fpc_total([1, 2, 3], 3)['se'] == 0


def test_bootstrap_interval_is_seeded_shifted_by_the_offset_and_degenerate_for_constants():
    values = [3, 5, 8, 13, 21, 34]
    lo, hi = sfs.bootstrap_interval(values, 100, reps=2000, seed=4)
    assert (lo, hi) == sfs.bootstrap_interval(values, 100, reps=2000, seed=4)
    assert lo < 100 * sum(values) / len(values) < hi
    lo2, hi2 = sfs.bootstrap_interval(values, 100, offset=1000.0, reps=2000, seed=4)
    assert (lo2, hi2) == pytest.approx((lo + 1000.0, hi + 1000.0))
    assert sfs.bootstrap_interval([4, 4, 4], 50, reps=500) == (200.0, 200.0)


def test_tail_summary():
    t = sfs.tail_summary([1, 2, 3, 4, 5, 100])
    assert (t['median'], t['max'], t['n']) == (3.5, 100, 6)
    assert t['mean'] == pytest.approx(115 / 6) and t['top5_share'] == pytest.approx(114 / 115)


def test_capped_records_are_counted_as_recorded_or_dropped_and_reweighted():
    records = [dict(status='done', codes=10), dict(status='done', codes=20), dict(status='done', codes=30),
               dict(status='capped', codes=100)]
    both = sfs.estimate_two_ways(records, 'codes', 100, 1000, seed=1)
    assert both['n_capped'] == 1
    assert both['as_recorded']['total'] == 1000 + 100 * 40       # mean of all four, the capped one as recorded
    assert both['capped_excluded']['total'] == 1000 + 100 * 20   # mean of the three done, reweighted to 100
    none = sfs.estimate_two_ways(records[:3], 'codes', 100, 1000, seed=1)
    assert none['n_capped'] == 0 and none['capped_excluded'] is None


def test_overlap_summary_counts_multiplicities():
    a, b, c = (sfs.code_hash(x) for x in 'abc')
    ov = sfs.overlap_summary([{a, b}, {b, c}, {c}])
    assert (ov['code_incidences'], ov['distinct_codes'], ov['codes_in_more_than_one']) == (5, 3, 2)
    assert ov['multiplicity_histogram'] == {'1': 1, '2': 2}
    assert ov['duplicate_fraction'] == pytest.approx(1 - 3 / 5)


def test_projection_arithmetic():
    p = sfs.projection(45260)
    assert p['search_cpu_h'] == pytest.approx([3.29, 9.28])
    assert p['checker_cpu_h'] == pytest.approx([3.29 * 0.24, 9.28 * 0.83])
    assert p['total_cpu_h'] == pytest.approx([3.29 * 1.24, 9.28 * 1.83])
    assert p['days'] == pytest.approx([3.29 * 1.24 / 64 / 24, 9.28 * 1.83 / 64 / 24])


# ------------------------------------------------------------------ the measurement of one triangulation

SMALL = 7594      # a triangulation whose enumeration takes under a second


def test_counts_of_a_real_triangulation_agree_with_the_enumerator(triangulations):
    from kisscrit.completeness import candidates as C
    F = triangulations[SMALL]
    record, hashes = sfs.enumerate_triangulation(F, None)
    assert record['status'] == 'done'
    assert record['tilings'] == sum(1 for _ in C.tilings(F))
    yielded = [code for code, _ in C.candidates_of_triangulation(F)]
    assert record['codes'] == len(yielded) == len(hashes) and hashes == {sfs.code_hash(c) for c in yielded}
    assert record['codes_seen'] >= record['codes'] > 0
    again, hashes_again = sfs.enumerate_triangulation(F, None)
    assert hashes_again == hashes and again['tilings'] == record['tilings']


def _never_yields(F, seen):
    while True:
        pass


def _yields_forever(F, seen):
    k = 0
    while True:
        assert k not in seen
        seen.add(k)
        k += 1
        yield str(k), []


def _finite(F, seen):
    for k in range(5):
        assert k not in seen
        seen.add(k)
        yield 'code%d' % k, []


@pytest.mark.parametrize('enumerator', [_never_yields, _yields_forever])
def test_cap_interrupts_the_enumeration_and_restores_the_alarm(enumerator):
    handler = signal.getsignal(signal.SIGALRM)
    record, hashes = sfs.enumerate_triangulation(None, 0.3, enumerator=enumerator)
    assert record['status'] == 'capped' and 0.25 <= record['wall_s'] < 5
    assert record['codes'] == len(hashes) <= record['tilings']
    assert signal.getitimer(signal.ITIMER_REAL) == (0.0, 0.0) and signal.getsignal(signal.SIGALRM) == handler


def test_a_finished_enumeration_under_a_cap_is_done_and_leaves_no_alarm():
    handler = signal.getsignal(signal.SIGALRM)
    record, hashes = sfs.enumerate_triangulation(None, 30.0, enumerator=_finite)
    assert (record['status'], record['tilings'], record['codes']) == ('done', 5, 5) and len(hashes) == 5
    assert signal.getitimer(signal.ITIMER_REAL) == (0.0, 0.0) and signal.getsignal(signal.SIGALRM) == handler


# ------------------------------------------------------------------ calibration

H = sfs.code_hash
SYNTH_GEN = ({1: {H('a'), H('b')}, 5: {H('c')}, 7: {H('d'), H('e')}}, {1: 2, 5: 1, 7: 2})


def _run(status, *codes):
    return dict(status=status), {H(c) for c in codes}


def _exact():
    # Index 5 contains a code of index 1 ('a'): the gen run attributes it to 1 only, an independent run finds it again.
    return {1: _run('done', 'a', 'b'), 5: _run('done', 'a', 'c'), 7: _run('done', 'd', 'e')}


def test_calibration_accepts_an_exact_reproduction():
    details, ok = sfs.calibrate(_exact(), *SYNTH_GEN)
    assert ok and {i: d['verdict'] for i, d in details.items()} == {1: 'match', 5: 'match', 7: 'match'}


def test_calibration_accepts_a_capped_prefix():
    results = _exact()
    results[7] = _run('capped', 'd')
    details, ok = sfs.calibrate(results, *SYNTH_GEN)
    assert ok and details[7]['verdict'] == 'prefix-consistent' and details[1]['verdict'] == 'match'


@pytest.mark.parametrize('mutation', ['extra code', 'missing code', 'foreign code at 5', 'capped with foreign', 'marker', 'error'])
def test_calibration_detects_every_mutation(mutation):
    results, (by_tri, markers) = _exact(), SYNTH_GEN
    if mutation == 'extra code':
        results[1] = _run('done', 'a', 'b', 'zzz')
    elif mutation == 'missing code':
        results[7] = _run('done', 'd')
    elif mutation == 'foreign code at 5':
        results[5] = _run('done', 'a', 'c', 'd')          # 'd' belongs to another shard
    elif mutation == 'capped with foreign':
        results[1] = _run('capped', 'a', 'zzz')
    elif mutation == 'marker':
        markers = {**markers, 7: 3}                       # the file's done marker disagrees with the codes it lists
    elif mutation == 'error':
        results[5] = (dict(status='error', error='ValueError'), set())
    details, ok = sfs.calibrate(results, by_tri, markers)
    assert not ok and {d['verdict'] for d in details.values()} & {'mismatch', 'error'}


@pytest.mark.parametrize('index', sfs.KNOWN)
def test_real_enumeration_prefix_contains_only_gen_file_codes(triangulations, gen, index):
    record, hashes = sfs.enumerate_triangulation(triangulations[index], 2.0)
    details, ok = sfs.calibrate({index: (record, hashes)}, *gen)
    assert len(hashes) > 100 and ok and details[index]['not_in_gen_file'] == 0
    assert details[index]['verdict'] in ('prefix-consistent', 'match')


# ------------------------------------------------------------------ the driver with a fake measurement

class FakeMeasurement:
    """Replaces `measure_index`: deterministic per-index code sets, one code shared by every triangulation,
    optional overrides for calibration runs, one optional failing index.  Records the indices it was asked for."""

    def __init__(self, overrides=None, failing=None):
        self.calls, self.overrides, self.failing = [], overrides or {}, failing

    def __call__(self, index, cap_seconds):
        self.calls.append(index)
        if index == self.failing:
            raise RuntimeError('boom at some/path')
        if index in self.overrides:
            hashes = set(self.overrides[index])
        else:
            hashes = {sfs.code_hash('%d-%d' % (index, k)) for k in range(5 + index % 9)} | {sfs.code_hash('shared')}
        status = 'capped' if index % 11 == 0 and index not in self.overrides else 'done'
        return dict(status=status, tilings=3 * len(hashes), codes_seen=len(hashes), codes=len(hashes), wall_s=0.01, cpu_s=0.02), hashes


def _drive(monkeypatch, tmp_path, fake, extra=(), seed='3', n='12'):
    monkeypatch.setattr(sfs, 'measure_index', fake)
    out = tmp_path / 'result.json'
    status = sfs.main(['--out', str(out), '--n', n, '--seed', seed, '--workers', '1', '--bootstrap-reps', '200'] + list(extra))
    return status, json.loads(out.read_text()), (tmp_path / 'result.summary.txt').read_text()


def test_driver_writes_estimates_overlap_and_the_closing_sentence(monkeypatch, tmp_path, capsys):
    fake = FakeMeasurement()
    status, result, summary = _drive(monkeypatch, tmp_path, fake)
    assert status == 0 and fake.calls == result['header']['sample_indices']
    indices = result['header']['sample_indices']
    assert len(indices) == 12 and indices == sorted(set(indices)) and not set(indices) & set(sfs.COND)
    assert result['header']['population_size'] == 7591 and result['header']['seed'] == 3
    assert [r['index'] for r in result['triangulations']] == indices
    # Every fake triangulation holds 'shared', so exactly that one code occurs in all twelve.
    ov = result['overlap']
    assert ov['multiplicity_histogram']['12'] == 1 and ov['codes_in_more_than_one'] == 1
    assert ov['code_incidences'] == sum(r['codes'] for r in result['triangulations'])
    assert result['estimates']['codes']['n_capped'] == sum(1 for i in indices if i % 11 == 0)
    assert summary.strip().splitlines()[-1].startswith('The full run generates about ')
    assert summary.strip().endswith('(projection from a sample of 12 triangulations).')
    assert result['summary'] == summary.strip() and summary.strip() in capsys.readouterr().out
    assert 'calibration' not in result


def test_driver_is_reproducible_from_the_seed(monkeypatch, tmp_path):
    first = _drive(monkeypatch, tmp_path, FakeMeasurement())[1]['header']['sample_indices']
    again = _drive(monkeypatch, tmp_path, FakeMeasurement())[1]['header']['sample_indices']
    other = _drive(monkeypatch, tmp_path, FakeMeasurement(), seed='4')[1]['header']['sample_indices']
    assert first == again != other


def test_driver_with_an_exact_calibration_goes_on_to_sample(monkeypatch, tmp_path, gen):
    by_tri, _ = gen
    overrides = {1: by_tri[1], 5: by_tri[5] | set(list(by_tri[1])[:100]), 7: by_tri[7]}
    fake = FakeMeasurement(overrides)
    status, result, summary = _drive(monkeypatch, tmp_path, fake, extra=['--include-known'])
    assert status == 0 and result['calibration']['ok'] and fake.calls[:3] == [1, 5, 7]       # one worker: calibration runs in order
    assert {d['verdict'] for d in result['calibration']['indices'].values()} == {'match'}
    example = result['calibration']['cond_overlap_example']
    assert example['union'] == 45260 and example['union_equals_repo_count'] and example['factor'] > 1
    assert re.search(r'index 1 match \(38,523 codes, [\d.]+ s wall\), index 5 match \([\d,]+ codes, [\d.]+ s wall\), '
                     r'index 7 match \(36,028 codes, [\d.]+ s wall\)', summary)
    assert result['calibration']['indices']['7']['tilings'] == 3 * 36028      # the fake's tiling count, carried into the JSON


def test_driver_stops_before_sampling_when_the_calibration_mismatches(monkeypatch, tmp_path, gen):
    by_tri, _ = gen
    broken = set(by_tri[7])
    broken.pop()                                           # one code of index 7 is missing from the run
    fake = FakeMeasurement({1: by_tri[1], 5: by_tri[5], 7: broken})
    status, result, summary = _drive(monkeypatch, tmp_path, fake, extra=['--include-known'])
    assert status == 2 and result['calibration']['ok'] is False
    assert result['calibration']['indices']['7']['verdict'] == 'mismatch' and result['calibration']['indices']['7']['missing_from_run'] == 1
    assert 'estimates' not in result and set(fake.calls) <= set(sfs.KNOWN)        # nothing was sampled
    assert 'FAILED' in summary


def test_driver_reports_an_exception_without_an_estimate(monkeypatch, tmp_path):
    indices = _drive(monkeypatch, tmp_path, FakeMeasurement())[1]['header']['sample_indices']
    status, result, summary = _drive(monkeypatch, tmp_path, FakeMeasurement(failing=indices[4]))
    assert status == 3 and 'estimates' not in result
    bad = [r for r in result['triangulations'] if r['status'] == 'error']
    assert bad == [dict(index=indices[4], status='error', error='RuntimeError')]       # the class name only
    assert 'some/path' not in json.dumps(result)


def test_a_bad_sample_size_is_refused(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(sfs, 'measure_index', FakeMeasurement())
    assert sfs.main(['--out', str(tmp_path / 'r.json'), '--n', '0']) == 1
    assert 'n must be between' in capsys.readouterr().out and not (tmp_path / 'r.json').exists()


# ------------------------------------------------------------------ no private strings, no reserved word

def _private_strings(tmp_path):
    values = {os.getcwd(), str(ROOT), str(tmp_path), os.path.expanduser('~'), os.environ.get('HOME'), socket.gethostname()}
    try:
        values.add(getpass.getuser())
    except Exception:
        pass
    return {v for v in values if v and len(v) >= 4}


def _leaks(text, tmp_path):
    found = [v for v in _private_strings(tmp_path) if re.search(r'(?<![\w.-])' + re.escape(v) + r'(?![\w.-])', text)]
    if re.search(r'(?<![\w./-])/(?:home|Users|root|tmp|var|mnt|opt|private)/', text):
        found.append('an absolute path')
    return found


def test_outputs_carry_no_host_user_or_path_and_never_use_the_reserved_word(monkeypatch, tmp_path, capsys, gen):
    by_tri, _ = gen
    fake = FakeMeasurement({1: by_tri[1], 5: by_tri[5], 7: by_tri[7]}, failing=None)
    _drive(monkeypatch, tmp_path, fake, extra=['--include-known'])
    printed = capsys.readouterr()
    outputs = [printed.out, printed.err, (tmp_path / 'result.json').read_text(), (tmp_path / 'result.summary.txt').read_text()]
    for text in outputs:
        assert not _leaks(text, tmp_path)
        assert RESERVED_WORD not in text.lower()
    # The same check on a failed run, whose messages are written by a different branch of the driver.
    fake = FakeMeasurement({1: by_tri[1], 5: by_tri[5], 7: set()})
    _drive(monkeypatch, tmp_path, fake, extra=['--include-known'])
    printed = capsys.readouterr()
    for text in (printed.out, printed.err, (tmp_path / 'result.json').read_text()):
        assert not _leaks(text, tmp_path)


def test_the_script_source_has_no_absolute_path_and_never_uses_the_reserved_word():
    text = (ROOT / 'scripts' / 'measure' / 'sample_full_size.py').read_text()
    assert not re.search(r'(?<![\w./-])/(?:home|Users|root|tmp|var|mnt|opt|private)/', text)
    assert RESERVED_WORD not in text.lower()
