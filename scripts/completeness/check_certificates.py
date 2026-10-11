#!/usr/bin/env python3
"""Independent check of completeness certificates (kisscrit/core/checker.py).

Checks every certificate in the given files (default: all of results/completeness/certificates/), in parallel,
skipping those already checked.  Output, one line per candidate, appended to
results/completeness/checks/checks-<name of the certificate file>:
    {"code", "verified", "failed": [...], "terminal": [...], "open", "not_run", "counts", "seconds", "checker"}

  python scripts/completeness/check_certificates.py [--workers 4] [files ...]
  python scripts/completeness/check_certificates.py --summary
Works on partial results: run it again after `git pull` to check newly pushed shards.
"""
import argparse, glob, json, os, sys, time
from multiprocessing import Pool

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
CERTS = os.path.join(ROOT, 'results', 'completeness', 'certificates')
CHECKS = os.path.join(ROOT, 'results', 'completeness', 'checks')
CHECKER_VERSION = 'checker-1'


def check_line(line):
    from kisscrit.core.checker import check_certificate
    d = json.loads(line); t0 = time.time()
    try:
        r = check_certificate(d['certificate'])
    except Exception as e:                                     # a malformed certificate is a failure
        r = dict(verified=False, error=repr(e), failed=[], terminal=[], open=0, not_run=0, counts={})
    r['code'] = d['code']; r['seconds'] = round(time.time() - t0, 2); r['checker'] = CHECKER_VERSION
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('files', nargs='*')
    ap.add_argument('--workers', type=int, default=os.cpu_count())
    ap.add_argument('--summary', action='store_true')
    a = ap.parse_args()
    os.makedirs(CHECKS, exist_ok=True)
    files = a.files or sorted(glob.glob(os.path.join(CERTS, 'certificates-*.jsonl')))
    if a.summary:
        tot = dict(certificates=0, checked=0, verified=0, failed_boxes=0, terminal_boxes=0, not_verified=[])
        for f in files:
            codes = [json.loads(l)['code'] for l in open(f)]
            tot['certificates'] += len(codes)
            out = os.path.join(CHECKS, 'checks-' + os.path.basename(f))
            if not os.path.exists(out): continue
            for l in open(out):
                r = json.loads(l); tot['checked'] += 1; tot['verified'] += bool(r['verified'])
                tot['failed_boxes'] += len(r['failed']); tot['terminal_boxes'] += len(r['terminal'])
                if not r['verified']: tot['not_verified'].append(r['code'])
        print(json.dumps({k: (v if k != 'not_verified' else len(v)) for k, v in tot.items()}))
        return
    for f in files:
        out = os.path.join(CHECKS, 'checks-' + os.path.basename(f))
        done = set()
        if os.path.exists(out):
            done = {json.loads(l)['code'] for l in open(out)}
        todo = [l for l in open(f) if json.loads(l)['code'] not in done]
        if not todo: continue
        print('%s: %d to check' % (os.path.basename(f), len(todo)), flush=True)
        t0 = time.time(); n = 0; bad = 0
        with open(out, 'a') as fo, Pool(a.workers) as pool:
            for r in pool.imap_unordered(check_line, todo, chunksize=1):
                fo.write(json.dumps(r) + '\n'); fo.flush(); n += 1; bad += (not r['verified'])
                if n % 20 == 0 or not r['verified']:
                    print('  %d/%d checked, %d not verified (%.0fs)' % (n, len(todo), bad, time.time() - t0), flush=True)
    print('done')


if __name__ == '__main__':
    main()
