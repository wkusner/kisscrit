#!/usr/bin/env python3
"""Independent check of completeness certificates (kisscrit/core/checker.py), with the manifest step.

Status: trusted. Part of the trusted base (transfer prompt v1.2, code standards).
Implements: the checker front end: V5 (manifest) and V7 (summary counters) of the checker spec v0.2; paper Appendix D.

Checks every certificate in the given files (default: all of results/completeness/certificates/), in parallel,
skipping those already checked.  Output, one line per candidate, appended to
results/completeness/checks/<checker>-<name of the certificate file>:
    {"code", "verified", "structural": [...], "failed": [...], "terminal": [...], "open", "not_run",
     "counts", "seconds", "checker"}

  python scripts/completeness/check_certificates.py [--route cond] [--workers 4] [files ...]
  python scripts/completeness/check_certificates.py --summary [--route cond]
  python scripts/completeness/check_certificates.py --manifest --route cond

Works on partial results: run it again after `git pull` to check newly pushed shards.

--summary (V7) reports verified / not verified counts and the checker's own work: totals of bisect,
linear-hint, linear-bound, linear-parent and the codes that needed any bisection.
--manifest (V5) compares the code sets of gen, cert, certificates and checks for a route: counts of missing,
duplicate and unverified codes.  Completeness can be claimed only if the four sets are equal, no code is
duplicated, and every code is verified or is `terminal` with every terminal box identified (checklist C4).

The gen faces are read from results/completeness/gen/ so that V4 (faces match the gen record) can be checked;
if the gen files are absent the code identity is still checked from the certificate's own faces.
"""
import argparse, collections, glob, gzip, json, os, sys, time
from multiprocessing import Pool

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
BASE = os.path.join(ROOT, 'results', 'completeness')
CERTS = os.path.join(BASE, 'certificates')
CHECKS = os.path.join(BASE, 'checks')
CHECKER_VERSION = 'checker-2-reference'      # the Python reference with V1-V5; 'checker-1' was the old checker
_GEN = {}


def gen_faces(route):
    """code -> faces from the gen files of a route (empty if none)."""
    out = {}
    for f in sorted(glob.glob(os.path.join(BASE, 'gen', 'gen-%s-*.jsonl.gz' % route))):
        with gzip.open(f, 'rt') as fh:
            for line in fh:
                d = json.loads(line)
                if 'code' in d: out.setdefault(d['code'], d.get('faces'))
    return out


def _init(gen):
    """Pool initializer: the gen faces by code."""
    global _GEN
    _GEN = gen


def check_line(line):
    """Check one certificate line; a malformed certificate is reported as a failure, never accepted."""
    from kisscrit.core.checker import check_certificate
    d = json.loads(line); t0 = time.time()
    try:
        r = check_certificate(d['certificate'], code=d['code'], gen_faces=_GEN.get(d['code']))
    except Exception as e:                                     # a malformed certificate is a failure
        r = dict(verified=False, error=repr(e), structural=[], failed=[], terminal=[], open=0, not_run=0, counts={})
    r['code'] = d['code']; r['seconds'] = round(time.time() - t0, 2); r['checker'] = CHECKER_VERSION
    return r


def _codes_of(path, key='code', gz=False):
    """The list of `key` values of the JSON lines of a file (gzip if gz)."""
    opener = (lambda p: gzip.open(p, 'rt')) if gz else open
    out = []
    with opener(path) as fh:
        for line in fh:
            d = json.loads(line)
            if key in d: out.append(d[key])
    return out


def manifest(route, base=BASE):
    """V5: the code sets of gen, cert, certificates and checks for a route, with duplicates, missing and unverified counts; `complete` only if all sets agree, nothing is duplicated, everything is verified and no terminal box is unidentified."""
    BASE, CERTS, CHECKS = base, os.path.join(base, 'certificates'), os.path.join(base, 'checks')
    sets = {}
    sets['gen'] = sum((_codes_of(f, gz=True) for f in sorted(glob.glob(os.path.join(BASE, 'gen', 'gen-%s-*.jsonl.gz' % route)))), [])
    sets['cert'] = sum((_codes_of(f) for f in sorted(glob.glob(os.path.join(BASE, 'cert', 'cert-%s-*.jsonl' % route)))), [])
    sets['certificates'] = sum((_codes_of(f) for f in sorted(glob.glob(os.path.join(CERTS, 'certificates-%s-*.jsonl' % route)))), [])
    checks = []
    verified = set(); terminal = {}
    for f in sorted(glob.glob(os.path.join(CHECKS, '%s-certificates-%s-*.jsonl' % (CHECKER_VERSION, route)))):
        for line in open(f):
            r = json.loads(line); checks.append(r['code'])
            if r['verified']: verified.add(r['code'])
            if r.get('terminal'): terminal[r['code']] = len(r['terminal'])
    sets['checks'] = checks
    status = collections.Counter()
    for f in sorted(glob.glob(os.path.join(BASE, 'cert', 'cert-%s-*.jsonl' % route))):
        for line in open(f):
            d = json.loads(line)
            if 'status' in d: status[d['status']] += 1
    distinct = {k: set(v) for k, v in sets.items()}
    report = dict(route=route, checker=CHECKER_VERSION, status=dict(status))
    for k, v in sets.items():
        c = collections.Counter(v)
        report[k] = dict(lines=len(v), distinct=len(c), duplicates=sum(n - 1 for n in c.values() if n > 1))
    base = distinct['gen']          # the candidate list; without it nothing can be called complete
    for k in ('cert', 'certificates', 'checks'):
        report['missing_from_' + k] = len(base - distinct[k]); report['extra_in_' + k] = len(distinct[k] - base)
    report['verified'] = len(verified & base)
    report['not_verified'] = len(base - verified)
    report['terminal_unidentified'] = len([c for c in terminal if c in base])
    equal = all(distinct[k] == base for k in ('cert', 'certificates', 'checks')) and \
        all(report[k]['duplicates'] == 0 for k in ('cert', 'certificates', 'checks'))
    report['sets_equal_no_duplicates'] = equal
    report['complete'] = bool(base) and equal and report['not_verified'] == 0 and report['terminal_unidentified'] == 0
    report['note'] = ('gen lines count duplicates across shards; the distinct count is the candidate count. '
                      'Completeness requires sets_equal_no_duplicates, not_verified == 0 and every terminal box '
                      'identified (checklist C4); terminal identification is not yet implemented, so '
                      'terminal_unidentified counts candidates with terminal boxes.')
    return report


def summary(files):
    """V7: verified and failed counts plus the checker's own work (bisections, hint and fallback counts)."""
    tot = dict(certificates=0, checked=0, verified=0, structural_failures=0, failed_boxes=0, terminal_boxes=0,
               open_leaves=0, not_verified=[], bisected=[], counts=collections.Counter(), seconds=0.0)
    for f in files:
        tot['certificates'] += len(_codes_of(f))
        out = os.path.join(CHECKS, '%s-%s' % (CHECKER_VERSION, os.path.basename(f)))
        if not os.path.exists(out): continue
        for l in open(out):
            r = json.loads(l); tot['checked'] += 1; tot['verified'] += bool(r['verified'])
            tot['structural_failures'] += bool(r.get('structural'))
            tot['failed_boxes'] += len(r['failed']); tot['terminal_boxes'] += len(r['terminal'])
            tot['open_leaves'] += r.get('open', 0); tot['seconds'] += r.get('seconds', 0)
            for k, v in r.get('counts', {}).items(): tot['counts'][k] += v
            if r.get('counts', {}).get('bisect'): tot['bisected'].append(r['code'])
            if not r['verified']: tot['not_verified'].append(r['code'])
    tot['counts'] = dict(tot['counts'])
    tot['n_not_verified'] = len(tot.pop('not_verified')); tot['n_bisected'] = len(tot['bisected'])
    tot['bisected'] = tot['bisected'][:20]
    return tot


def main():
    """Command line entry point."""
    ap = argparse.ArgumentParser()
    ap.add_argument('files', nargs='*')
    ap.add_argument('--route', default='cond')
    ap.add_argument('--workers', type=int, default=os.cpu_count())
    ap.add_argument('--summary', action='store_true')
    ap.add_argument('--manifest', action='store_true')
    a = ap.parse_args()
    os.makedirs(CHECKS, exist_ok=True)
    files = a.files or sorted(glob.glob(os.path.join(CERTS, 'certificates-%s-*.jsonl' % a.route)))
    if a.manifest:
        print(json.dumps(manifest(a.route), indent=1)); return
    if a.summary:
        print(json.dumps(summary(files))); return
    gen = gen_faces(a.route)
    for f in files:
        out = os.path.join(CHECKS, '%s-%s' % (CHECKER_VERSION, os.path.basename(f)))
        done = set()
        if os.path.exists(out):
            done = {json.loads(l)['code'] for l in open(out)}
        todo = [l for l in open(f) if json.loads(l)['code'] not in done]
        if not todo: continue
        print('%s: %d to check' % (os.path.basename(f), len(todo)), flush=True)
        t0 = time.time(); n = 0; bad = 0
        with open(out, 'a') as fo, Pool(a.workers, initializer=_init, initargs=(gen,)) as pool:
            for r in pool.imap_unordered(check_line, todo, chunksize=1):
                fo.write(json.dumps(r) + '\n'); fo.flush(); n += 1; bad += (not r['verified'])
                if n % 20 == 0 or not r['verified']:
                    print('  %d/%d checked, %d not verified (%.0fs)' % (n, len(todo), bad, time.time() - t0), flush=True)
    print('done')


if __name__ == '__main__':
    main()
