#!/usr/bin/env python3
"""Sampled size of the `full` completeness run: candidates, CPU-hours and days, as a labelled projection.

Status: untrusted (measurement tooling; no part of the proof; every number it reports is a projection from a
random sample of triangulations, never the count of record).

What it measures.  The `full` run enumerates, for each of the 7,595 triangulations of the sphere with 12 vertices,
the tilings by triangles, quadrilaterals and pentagons that pass the vertex pruning, keeps the distinct canonical
codes that pass the linear relaxation (`candidates_of_triangulation` in kisscrit/completeness/candidates.py,
used here unchanged), and then searches and checks each distinct code.  The four `cond` triangulations (indices 0,
1, 5, 7) are already enumerated in results/completeness/gen.  This script draws a seeded simple random sample of the
other 7,591 triangulations, runs the unchanged enumerator on each under a wall-clock cap, and scales the sample up.

What the numbers mean (read before quoting any of them).
  * "codes" of a triangulation = the distinct canonical codes the enumerator yields for it alone, i.e. what the gen
    files record per triangulation.  "codes_seen" counts the distinct codes before the linear relaxation, "tilings"
    the tilings produced by the depth-first search.
  * The estimate  45,260 + 7,591 x (sample mean of codes)  adds the codes of every triangulation separately, so a code
    that arises from several triangulations is counted several times: the sum it estimates is an upper bound for the
    number of distinct candidates, and the sampling error comes on top.  (In the cond gen files the 45,260 distinct codes occur in 1, 2 or 3 of the four triangulations: 9,140,
    26,664 and 9,456 codes.)  The overlap measured inside the sample understates the overlap in the whole population (a
    code's other triangulations are each in the sample only with probability n/7,591), so it is reported, not extrapolated.
  * A triangulation that hits the cap records only the tilings enumerated so far.  Counted as recorded, the estimate
    is a lower bound for the same sum; dropped and reweighted it is biased low whenever the capped ones are the large
    ones.  With a capped triangulation present, the total is therefore bracketed from neither side.
  * The CPU-hour projection multiplies the candidate count by two measured pilot rates (a server and an older desktop
    run) and the checker by the measured checker/search time ratio of docs/decisions.md (2026-10-20); it is crude
    because the cost per candidate varies widely.  Days assume perfect use of 64 cores.
  * No tolerance is introduced here.  Counts are integers produced by the enumerator, whose floating-point steps
    (the cell masks with slack 1e-9 and the HiGHS discard) are the assumptions G2 and G1 of the table in STATUS.md.

Calibration (--include-known).  Before sampling, indices 1, 5 and 7 are enumerated and compared with the gen files
code by code (hashes of the canonical code strings).  A completed index must reproduce the file exactly; a capped one
must at least find only codes that the file lists.  Any other outcome stops the run with exit status 2.

  python scripts/measure/sample_full_size.py --out sample.json [--seed 1] [--n 150] [--cap-seconds 1200]
         [--workers N] [--include-known] [--calib-cap-seconds 3600]

Writes <out> (JSON: header, calibration, per-triangulation records, estimates, overlap, projection) and
<out minus .json>.summary.txt (the text summary, also printed).  Linux or macOS only (the cap uses SIGALRM).
Exit status: 0 done; 2 calibration mismatch; 3 an enumeration raised an exception.
"""
import os

# One thread per process: the enumerator is single-threaded Python plus one small LP per candidate, and library
# threads in 64 worker processes would oversubscribe the cores and distort the CPU-time measurement.
for _var in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_var, '1')

import argparse
import collections
import gzip
import hashlib
import json
import math
import multiprocessing
import pathlib
import platform
import random
import signal
import subprocess
import sys
import time

import numpy

# The repository root is found from this file, so the script runs from any working directory and from a clean clone.
ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

TOOL_VERSION = 'sample_full_size 1.0'
COND = (0, 1, 5, 7)                  # the triangulations of the `cond` route (run_shard.py ROUTES['cond'])
KNOWN = (1, 5, 7)                    # calibration indices; index 0 is not re-run, its codes are read from the gen file
COND_CODES = 45260                   # distinct codes of results/completeness/gen (the test re-derives it from the files)

# Pilot rates in search CPU-hours per 45,260 candidates, both from files in this repository:
SERVER_CPU_H = 3.29                  # 64-core server, results/pilot/cond-2026-10-22/summary-cond.json (3.287)
DESKTOP_CPU_H = 9.28                 # older desktop run, results/completeness/summary-cond.json (9.283)
PILOT_CODES = 45260
# Checker CPU time / search CPU time on the same candidates, docs/decisions.md 2026-10-20 (F4 performance gate):
# 0.24 on 242 fresh `cond` runs, about 0.83 without the two numba-compile candidates (2-core VM).
CHECKER_RATIO = (0.24, 0.83)
CHECKER_RATIO_DATE = '2026-10-20'
CORES_FOR_DAYS = 64


# --------------------------------------------------------------------------------------------- data

def load_triangulations():
    """The 7,595 triangulations (lists of vertex triples) used by every phase of the run."""
    with open(ROOT / 'data' / 'completeness' / 'triangulations12.json') as f:
        return json.load(f)


def code_hash(code):
    """16-byte digest of a canonical code string; a worker returns these instead of the strings to keep the
    messages small, and 128 bits make an accidental collision negligible at the sizes involved."""
    return hashlib.blake2b(code.encode(), digest_size=16).digest()


def load_gen_cond():
    """Read results/completeness/gen/gen-cond-*.jsonl.gz.  Returns ({index: set of hashes of the codes the file
    attributes to that triangulation}, {index: 'new' count of its done marker}).  Within a shard the gen run shares one
    `seen` set, so a code found by an earlier triangulation of the shard is attributed to it only (index 5 after 1)."""
    by_tri, markers = collections.defaultdict(set), {}
    for path in sorted((ROOT / 'results' / 'completeness' / 'gen').glob('gen-cond-*.jsonl.gz')):
        with gzip.open(path, 'rt') as f:
            for line in f:
                d = json.loads(line)
                if 'code' in d:
                    by_tri[d['tri']].add(code_hash(d['code']))
                else:
                    markers[d['done_triangulation']] = d['new']
    return dict(by_tri), markers


# --------------------------------------------------------------------------------------------- measuring one triangulation

class _CountingSet(set):
    """A set that counts membership tests.  `candidates_of_triangulation` makes exactly one `code in seen` test per
    tiling that the depth-first search produces, so the number of tests is the tiling count and the enumerator itself
    stays untouched (tests/test_sample_full_size.py compares it with a direct count of `tilings`)."""

    def __init__(self):
        super().__init__()
        self.tests = 0

    def __contains__(self, item):
        self.tests += 1
        return super().__contains__(item)


class _Capped(Exception):
    """Raised by the alarm handler inside the enumeration when the wall-clock cap is reached."""


def _raise_capped(signum, frame):
    raise _Capped()


def enumerate_triangulation(triangulation, cap_seconds, enumerator=None):
    """Run the enumerator on one triangulation.  Returns (record, hashes): the record holds status ('done' or
    'capped'), tilings, codes_seen (distinct codes before the linear relaxation), codes (distinct codes yielded),
    wall_s and cpu_s; hashes is the set of code hashes yielded.  A capped run keeps what was found before the cap and
    is never extrapolated.  The cap is a real-time alarm, so it also interrupts a stretch of the depth-first search
    that yields nothing.  A cap of None or 0 means no cap.  `enumerator` exists only so tests can substitute a fake."""
    if enumerator is None:
        from kisscrit.completeness.candidates import candidates_of_triangulation as enumerator
    seen = _CountingSet()
    hashes, yielded, status = set(), 0, 'done'
    use_cap = bool(cap_seconds) and cap_seconds > 0
    wall0, cpu0 = time.perf_counter(), time.process_time()
    old_handler = signal.signal(signal.SIGALRM, _raise_capped) if use_cap else None
    try:
        try:
            if use_cap:
                signal.setitimer(signal.ITIMER_REAL, cap_seconds)
            for code, _faces in enumerator(triangulation, seen):
                hashes.add(code_hash(code))
                yielded += 1
        finally:
            # Cancelling inside the inner `finally` keeps a late alarm from firing after the loop; if it fires during
            # the cancellation itself, the outer handler below records the run as capped (conservative).
            if use_cap:
                signal.setitimer(signal.ITIMER_REAL, 0)
    except _Capped:
        status = 'capped'
    finally:
        if use_cap:
            signal.signal(signal.SIGALRM, old_handler)
    if yielded != len(hashes):
        raise RuntimeError('the enumerator yielded a code twice or two codes share a hash')
    record = dict(status=status, tilings=seen.tests, codes_seen=len(seen), codes=len(hashes),
                  wall_s=round(time.perf_counter() - wall0, 3), cpu_s=round(time.process_time() - cpu0, 3))
    return record, hashes


_TRIANGULATIONS = None


def _init_worker():
    """Load the triangulation list once per worker process."""
    global _TRIANGULATIONS
    _TRIANGULATIONS = load_triangulations()


def measure_index(index, cap_seconds):
    """Measure triangulation number `index`; the unit of work of one worker."""
    global _TRIANGULATIONS
    if _TRIANGULATIONS is None:
        _init_worker()
    return enumerate_triangulation(_TRIANGULATIONS[index], cap_seconds)


def _task(args):
    """Pool task: never lets an exception cross the process boundary with a traceback (it could carry local paths);
    the record keeps only the exception class name and the run stops with status 3 after writing what it has."""
    index, cap_seconds = args
    try:
        record, hashes = measure_index(index, cap_seconds)
    except Exception as exc:
        return dict(index=index, status='error', error=type(exc).__name__), set()
    record['index'] = index
    return record, hashes


def run_indices(indices, cap_seconds, workers, label, log):
    """Measure the given triangulation indices on `workers` processes; returns {index: (record, hashes)}.
    One progress line per finished triangulation (count, elapsed seconds, outcome)."""
    tasks = [(i, cap_seconds) for i in indices]
    out, t0 = {}, time.time()
    if workers <= 1:
        results = map(_task, tasks)
        pool = None
    else:
        pool = multiprocessing.Pool(workers, initializer=_init_worker)
        results = pool.imap_unordered(_task, tasks, chunksize=1)
    try:
        for done, (record, hashes) in enumerate(results, 1):
            out[record['index']] = (record, hashes)
            log('%s [%d/%d] elapsed %ds  index %d  %s  tilings %s  codes %s  wall %ss' % (
                label, done, len(tasks), time.time() - t0, record['index'], record['status'],
                record.get('tilings', '-'), record.get('codes', '-'), record.get('wall_s', '-')))
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    return out


# --------------------------------------------------------------------------------------------- sampling and estimators

def draw_sample(seed, n, population):
    """A simple random sample of n indices without replacement from `population`, sorted.  The indices are written
    to the output header, so the sample can be reproduced from the file whatever the Python version."""
    if not 0 < n <= len(population):
        raise ValueError('n must be between 1 and the population size %d' % len(population))
    return sorted(random.Random(seed).sample(sorted(population), n))


def fpc_total(values, population_size):
    """Estimate of the population total from a simple random sample without replacement:
    total = M x mean, standard error = M x sqrt((1 - n/M) s^2 / n) with s^2 the sample variance (n - 1) and the
    finite-population factor (1 - n/M); the standard error is None for a single value."""
    n = len(values)
    if n == 0:
        raise ValueError('empty sample')
    mean = sum(values) / n
    se = None
    if n > 1:
        var = sum((v - mean) ** 2 for v in values) / (n - 1)
        se = population_size * math.sqrt((1 - n / population_size) * var / n)
    return dict(n=n, mean=mean, total=population_size * mean, se=se)


def bootstrap_interval(values, population_size, offset=0.0, reps=10000, seed=0, level=95.0):
    """Percentile bootstrap interval for offset + M x mean: resample the sample with replacement `reps` times.
    With heavy-tailed counts and a few hundred values the interval tends to be too short on the upper side; the
    summary says so."""
    rng = numpy.random.default_rng([seed, 1])
    v = numpy.asarray(values, dtype=float)
    draws = rng.integers(0, len(v), size=(reps, len(v)))
    totals = offset + population_size * v[draws].mean(axis=1)
    tail = (100.0 - level) / 2.0
    lo, hi = numpy.percentile(totals, [tail, 100.0 - tail])
    return float(lo), float(hi)


def tail_summary(values):
    """Median, mean, maximum and the share of the sample total carried by the five largest values."""
    s = sorted(values)
    n, total = len(s), sum(s)
    return dict(n=n, median=(s[(n - 1) // 2] + s[n // 2]) / 2, mean=total / n, max=s[-1],
                top5_share=(sum(s[-5:]) / total) if total else None)


def estimate_two_ways(records, key, population_size, offset, seed):
    """Estimate offset + M x mean(key) from the sample records in two ways.  `as_recorded` uses every record, a
    capped one at what was enumerated before the cap: it estimates the sum of the counts found within the cap, which
    is at most the sum of the full counts.  `capped_excluded` drops the capped records and reweights the rest to
    M: not a bound, biased low if the capped triangulations are the large ones.  Each entry has the estimate, its
    standard error (without `offset`) and a bootstrap 95% interval; `capped_excluded` is None if nothing is done."""
    def one(values):
        est = fpc_total(values, population_size)
        lo, hi = bootstrap_interval(values, population_size, offset=offset, seed=seed)
        return dict(n=est['n'], mean=est['mean'], total=offset + est['total'], se=est['se'], ci95=[lo, hi])
    done = [r[key] for r in records if r['status'] == 'done']
    n_capped = len(records) - len(done)
    return dict(n_capped=n_capped, as_recorded=one([r[key] for r in records]),
                capped_excluded=one(done) if n_capped and done else None)


# --------------------------------------------------------------------------------------------- duplicates and calibration

def overlap_summary(hash_sets):
    """Overlap of the code sets of the sampled triangulations: code incidences (the sum of the per-triangulation
    counts), distinct codes, how many distinct codes occur in more than one sampled triangulation, and the
    histogram of multiplicities.  It is measured on the sample and understates the overlap of the population."""
    multiplicity = collections.Counter()
    for hs in hash_sets:
        multiplicity.update(hs)
    histogram = collections.Counter(multiplicity.values())
    incidences = sum(multiplicity.values())
    return dict(triangulations=len(hash_sets), code_incidences=incidences, distinct_codes=len(multiplicity),
                codes_in_more_than_one=sum(c for m, c in histogram.items() if m > 1),
                duplicate_fraction=(1 - len(multiplicity) / incidences) if incidences else None,
                multiplicity_histogram={str(m): histogram[m] for m in sorted(histogram)})


def calibrate(results, by_tri, markers):
    """Compare runs of the cond indices 1, 5, 7 (`results` = {index: (record, hashes)}) with the gen files.
    The file lists, per index, the codes that index found and the earlier index of its shard had not (5 after 1).
    A completed run must contain every code the file attributes to it and nothing outside what its shard lists
    (for 1 and 7 that is exact equality; for 5 it is equality up to the codes shared with 1).  A capped run must
    contain nothing outside it.  Verdicts: 'match', 'prefix-consistent' (capped, nothing foreign), 'mismatch',
    'error'.  Returns ({index: details}, ok) with ok False if any verdict is 'mismatch' or 'error'."""
    allowed = {1: by_tri.get(1, set()), 5: by_tri.get(1, set()) | by_tri.get(5, set()), 7: by_tri.get(7, set())}
    details, ok = {}, True
    for index in sorted(results):
        record, hashes = results[index]
        if record['status'] == 'error':
            details[index] = dict(verdict='error', error=record.get('error'))
            ok = False
            continue
        own = by_tri.get(index, set())
        foreign = len(hashes - allowed[index])
        missing = len(own - hashes) if record['status'] == 'done' else None
        file_consistent = (len(own) == markers.get(index))
        if foreign or missing or not file_consistent:
            verdict = 'mismatch'
        else:
            verdict = 'match' if record['status'] == 'done' else 'prefix-consistent'
        ok = ok and verdict != 'mismatch'
        details[index] = dict(verdict=verdict, status=record['status'], codes=len(hashes), gen_file_codes=len(own),
                              not_in_gen_file=foreign, missing_from_run=missing, gen_marker_matches=file_consistent,
                              tilings=record.get('tilings'), codes_seen=record.get('codes_seen'),
                              wall_s=record.get('wall_s'), cpu_s=record.get('cpu_s'))
    return details, ok


def cond_overlap_example(results, by_tri, cond_hashes):
    """Overlap among the four cond triangulations as a measured example: the sum of their independent code counts
    against their union.  Index 0 is taken from the gen file (its shard starts with an empty `seen`, so the file holds
    exactly its codes).  Returns None unless 1, 5 and 7 were completed."""
    if not all(i in results and results[i][0]['status'] == 'done' for i in KNOWN):
        return None
    sets = [by_tri[0]] + [results[i][1] for i in KNOWN]
    union = set().union(*sets)
    total = sum(len(s) for s in sets)
    return dict(indices=list(COND), sum_of_counts=total, union=len(union), union_equals_repo_count=(len(union) == COND_CODES),
                factor=total / len(union), cond_codes_in_repo=len(cond_hashes))


# --------------------------------------------------------------------------------------------- projection and text

def projection(n_codes):
    """CPU-hours and days for a number of candidates: search at the server and desktop pilot rates, checker at the
    two measured ratios, days on CORES_FOR_DAYS cores at perfect efficiency.  `low` pairs the server rate with the lower
    ratio, `high` the desktop rate with the higher one; the true cost per candidate varies, so these are a range."""
    search_low = n_codes * SERVER_CPU_H / PILOT_CODES
    search_high = n_codes * DESKTOP_CPU_H / PILOT_CODES
    check_low, check_high = search_low * CHECKER_RATIO[0], search_high * CHECKER_RATIO[1]
    total_low, total_high = search_low + check_low, search_high + check_high
    return dict(candidates=n_codes, search_cpu_h=[search_low, search_high], checker_cpu_h=[check_low, check_high],
                total_cpu_h=[total_low, total_high], days=[total_low / CORES_FOR_DAYS / 24, total_high / CORES_FOR_DAYS / 24])


def _fmt(x, digits=3):
    """`x` to `digits` significant figures, with thousands separators from 1,000 up."""
    if x == 0:
        return '0'
    magnitude = math.floor(math.log10(abs(x)))
    r = round(x, digits - 1 - magnitude)
    if abs(r) >= 1000:
        return '{:,.0f}'.format(r)
    return ('{:.%df}' % max(digits - 1 - magnitude, 0)).format(r)


def _range(pair):
    return '%s to %s' % (_fmt(pair[0]), _fmt(pair[1]))


def build_summary(result):
    """The human-readable summary, ending with the one-sentence projection.  Pure function of the result dict."""
    h, est, ov = result['header'], result['estimates'], result['overlap']
    n, M = h['n'], h['population_size']
    cal = result.get('calibration')
    a, b = est['codes'], est['codes_not_in_cond']
    capped = a['n_capped']
    lines = ['Sampled size of the `full` run (%s)' % TOOL_VERSION,
             'Tree %s%s   seed %d   n %d of %s non-cond triangulations   cap %s s   workers %d' % (
                 h['commit'], ' (modified)' if h['tree_modified'] else '', h['seed'], n, '{:,}'.format(M),
                 '{:,.0f}'.format(h['cap_seconds']), h['workers']),
             'LABEL: projection from a sample of %d triangulations; not a count of record (the full enumeration replaces it).' % n, '']
    if cal is None:
        lines.append('Calibration against the gen files: not run (--include-known was not given).')
    else:
        lines.append('Calibration against the gen files: ' + ', '.join(
            'index %d %s (%s codes, %s s wall)' % (i, cal['indices'][str(i)]['verdict'], '{:,}'.format(cal['indices'][str(i)]['codes']),
                                                   _fmt(cal['indices'][str(i)]['wall_s'] or 0))
            for i in KNOWN if str(i) in cal['indices']) + '.')
        if cal.get('cond_overlap_example'):
            ex = cal['cond_overlap_example']
            lines.append('  Overlap among the four cond triangulations (measured example, not representative): per-triangulation counts add up to %s, '
                         'their union is %s, factor %.2f.' % ('{:,}'.format(ex['sum_of_counts']), '{:,}'.format(ex['union']), ex['factor']))
    wall = est['wall_s_tail']
    lines += ['Sample outcome: %d done, %d capped.  Wall seconds per triangulation: median %s, mean %s, maximum %s (the longest one bounds the '
              'sample phase from below at any number of cores).' % (n - capped, capped, _fmt(wall['median']), _fmt(wall['mean']), _fmt(wall['max'])), '']
    t = est['codes_tail']
    lines += ['Codes per triangulation (as recorded): median %s, mean %s, maximum %s; the 5 largest carry %.1f%% of the sample total.' % (
        _fmt(t['median']), _fmt(t['mean']), '{:,}'.format(t['max']), 100 * (t['top5_share'] or 0)), '',
        'Candidates, summed over triangulations (a code arising from several triangulations is counted each time):',
        '  A. %s (cond, from the repo) + %s x sample mean = %s  (standard error of the sample part %s; bootstrap 95%% interval %s)' % (
            '{:,}'.format(COND_CODES), '{:,}'.format(M), _fmt(a['as_recorded']['total']), _fmt(a['as_recorded']['se'] or 0),
            _range(a['as_recorded']['ci95'])),
        '  B. same, with the codes already in the cond set removed from the sample: %s  (bootstrap 95%% interval %s)' % (
            _fmt(b['as_recorded']['total']), _range(b['as_recorded']['ci95']))]
    if capped:
        lines.append('  %d sampled triangulation(s) hit the cap. A and B count them as recorded (a lower bound for this sum); '
                     'with them dropped and reweighted A is %s (biased low if the capped ones are the large ones).' % (
                         capped, _fmt(a['capped_excluded']['total']) if a['capped_excluded'] else 'not available'))
    lines += ['  The bootstrap interval is a guide only: with heavy-tailed counts and n = %d it tends to be too short on the upper side.' % n, '',
              'Overlap measured on the sample: %s code incidences, %s distinct codes, %s codes in more than one sampled triangulation '
              '(duplicate fraction %.3f); %s of the sampled incidences are codes already in the cond set.' % (
                  '{:,}'.format(ov['code_incidences']), '{:,}'.format(ov['distinct_codes']), '{:,}'.format(ov['codes_in_more_than_one']),
                  ov['duplicate_fraction'] or 0, '{:,}'.format(ov['incidences_in_cond'])),
              'This understates the overlap of all 7,591 triangulations (a code\'s other triangulations are in the sample only with probability %s), '
              'so A estimates a sum that is an upper bound for the distinct candidates%s (the sampling error comes on top).' % (_fmt(n / M), ', provided nothing is capped' if capped else ''), '']
    gen = est['enumeration_cpu_h']
    lines += ['Enumeration CPU-hours for the non-cond triangulations (as recorded%s): %s (95%% interval %s).' % (
        ', a lower bound' if capped else '', _fmt(gen['total']), _range(gen['ci95'])), '']
    p = result['projection']
    lines += ['Projection (crude: search at %s and %s CPU-h per %s candidates, the server and desktop pilots; checker at %s to %s of the search time, '
              'docs/decisions.md %s; days at %d cores with perfect parallel efficiency):' % (
                  _fmt(SERVER_CPU_H), _fmt(DESKTOP_CPU_H), '{:,}'.format(PILOT_CODES), CHECKER_RATIO[0], CHECKER_RATIO[1],
                  CHECKER_RATIO_DATE, CORES_FOR_DAYS),
              '  %-22s %14s %26s %26s %20s' % ('candidates (A)', 'N', 'search + check CPU-h', 'of which checker', 'days')]
    for name in ('interval_low', 'point', 'interval_high'):
        row = p[name]
        lines.append('  %-22s %14s %26s %26s %20s' % (name.replace('_', ' '), _fmt(row['candidates']), _range(row['total_cpu_h']),
                                                      _range(row['checker_cpu_h']), _range(row['days'])))
    point, low, high = p['point'], p['interval_low'], p['interval_high']
    caveat = ('N is summed over triangulations: duplicates across triangulations are not removed.' if not capped else
              'N is bracketed from neither side: capped triangulations undercount, duplicates across triangulations overcount.')
    lines += ['',
              'H = search plus check CPU-hours. The ranges pair the server rate with the lower checker ratio and the desktop rate with the higher ratio '
              '(at the point estimate of N); the figures in brackets carry the 95%% bootstrap interval of N through the same pairing. %s' % caveat,
              'The full run generates about %s candidates (%s), about %s CPU-hours (%s), about %s days (%s) on the available servers '
              '(projection from a sample of %d triangulations).' % (
                  _fmt(point['candidates']), _range([low['candidates'], high['candidates']]),
                  _range(point['total_cpu_h']), _range([low['total_cpu_h'][0], high['total_cpu_h'][1]]),
                  _range(point['days']), _range([low['days'][0], high['days'][1]]), n)]
    return '\n'.join(lines)


# --------------------------------------------------------------------------------------------- driver

def _git(*args):
    """Output of a git command in the repository, or None if git or the repository is not available."""
    try:
        done = subprocess.run(['git', '-C', str(ROOT)] + list(args), capture_output=True, text=True, check=True)
        return done.stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def build_header(args, indices, n_triangulations):
    """Everything needed to reproduce the run, with no host, user or path strings."""
    import scipy, networkx
    commit = _git('rev-parse', 'HEAD')
    return dict(tool=TOOL_VERSION, commit=commit or 'unknown', tree_modified=bool(_git('status', '--porcelain', '--untracked-files=no')),
                seed=args.seed, n=args.n, cap_seconds=args.cap_seconds, calib_cap_seconds=args.calib_cap_seconds, workers=args.workers,
                include_known=args.include_known, population_size=n_triangulations - len(COND), triangulations=n_triangulations,
                cond_indices=list(COND), sample_indices=indices, bootstrap_reps=args.bootstrap_reps,
                versions=dict(python=platform.python_version(), numpy=numpy.__version__, scipy=scipy.__version__, networkx=networkx.__version__),
                platform=dict(system=platform.system(), machine=platform.machine(), cpu_count=os.cpu_count()),
                tolerances='none introduced here; counts are integers from the unchanged enumerator (its floating-point steps are assumptions G1 and G2 in STATUS.md)')


def _write(out_path, result, summary):
    """Write the JSON and the summary next to it."""
    out = pathlib.Path(out_path)
    out.write_text(json.dumps(result, indent=1, sort_keys=True) + '\n')
    out.with_suffix('.summary.txt').write_text(summary + '\n')


def parse_args(argv):
    ap = argparse.ArgumentParser(description='Sampled size of the `full` completeness run (a labelled projection).')
    ap.add_argument('--out', required=True, help='JSON output file; the summary goes next to it as <name>.summary.txt')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--n', type=int, default=150, help='number of sampled triangulations (default 150)')
    ap.add_argument('--cap-seconds', type=float, default=1200.0, help='wall-clock cap per sampled triangulation')
    ap.add_argument('--calib-cap-seconds', type=float, default=3600.0, help='wall-clock cap per calibration triangulation')
    ap.add_argument('--workers', type=int, default=os.cpu_count() or 1)
    ap.add_argument('--include-known', action='store_true', help='first re-enumerate indices 1, 5, 7 and compare with the gen files')
    ap.add_argument('--bootstrap-reps', type=int, default=10000)
    return ap.parse_args(argv)


def _stop(args, result, summary, status, log):
    """Record why the run stopped before an estimate, write the files, and return the exit status."""
    result['summary'] = summary
    _write(args.out, result, summary)
    log(summary)
    return status


def estimate_and_project(records, hash_sets, seed, population_size):
    """The estimates, the overlap and the projection from the finished sample.  `records` are the per-triangulation
    records with `codes_in_cond` filled in; `hash_sets` are the code hashes of the same triangulations."""
    for r in records:
        r['codes_not_in_cond'] = r['codes'] - r['codes_in_cond']
    est = dict(codes=estimate_two_ways(records, 'codes', population_size, COND_CODES, seed),
               codes_not_in_cond=estimate_two_ways(records, 'codes_not_in_cond', population_size, COND_CODES, seed),
               codes_tail=tail_summary([r['codes'] for r in records]),
               wall_s_tail=tail_summary([r['wall_s'] for r in records]))
    cpu_hours = [r['cpu_s'] / 3600.0 for r in records]
    cpu = fpc_total(cpu_hours, population_size)
    est['enumeration_cpu_h'] = dict(total=cpu['total'], se=cpu['se'],
                                    ci95=list(bootstrap_interval(cpu_hours, population_size, seed=seed)))
    overlap = overlap_summary(hash_sets)
    overlap['incidences_in_cond'] = sum(r['codes_in_cond'] for r in records)
    lo, hi = est['codes']['as_recorded']['ci95']
    proj = dict(point=projection(est['codes']['as_recorded']['total']), interval_low=projection(lo), interval_high=projection(hi))
    return est, overlap, proj


def main(argv=None):
    """Calibrate (optional), sample, estimate, project, write.  Returns the exit status."""
    args = parse_args(argv)
    log = lambda s: print(s, flush=True)
    triangulations = load_triangulations()
    population = [i for i in range(len(triangulations)) if i not in COND]
    try:
        indices = draw_sample(args.seed, args.n, population)
    except ValueError as exc:
        log(str(exc))
        return 1
    result = dict(header=build_header(args, indices, len(triangulations)))
    by_tri, markers = load_gen_cond()
    cond_hashes = set().union(*by_tri.values())
    if len(cond_hashes) != COND_CODES:
        # The estimate adds COND_CODES to the sample part; a gen directory that disagrees means the constant or the data is stale.
        return _stop(args, result, 'The gen files hold %d distinct codes, not the expected %d; stopped before measuring.' % (
            len(cond_hashes), COND_CODES), 2, log)

    if args.include_known:
        runs = run_indices(KNOWN, args.calib_cap_seconds, min(args.workers, len(KNOWN)), 'calibration', log)
        details, ok = calibrate(runs, by_tri, markers)
        result['calibration'] = dict(ok=ok, indices={str(i): d for i, d in details.items()},
                                     cond_overlap_example=cond_overlap_example(runs, by_tri, cond_hashes) if ok else None)
        if not ok:
            return _stop(args, result, 'Calibration against the gen files FAILED (%s); stopped before sampling. See the JSON.' % ', '.join(
                'index %d %s' % (i, d['verdict']) for i, d in sorted(details.items())), 2, log)

    runs = run_indices(indices, args.cap_seconds, args.workers, 'sample', log)
    if any(runs[i][0]['status'] == 'error' for i in indices):
        result['triangulations'] = [runs[i][0] for i in indices]
        return _stop(args, result, 'At least one enumeration raised an exception (see the records in the JSON); no estimate was made.', 3, log)
    records = [dict(runs[i][0], codes_in_cond=len(runs[i][1] & cond_hashes)) for i in indices]
    result['estimates'], result['overlap'], result['projection'] = estimate_and_project(
        records, [runs[i][1] for i in indices], args.seed, result['header']['population_size'])
    result['triangulations'] = records
    return _stop(args, result, build_summary(result), 0, log)


if __name__ == '__main__':
    sys.exit(main())
