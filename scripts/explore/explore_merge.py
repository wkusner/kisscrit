"""Merge the census files of the exploration workers into one census.

Status: untrusted (exploration or tooling; no part of the proof).

usage: python scripts/explore/explore_merge.py --out MERGED.json [--lower-dir DIR] [--no-euler] PATH [PATH ...]
PATH is a worker census file or a directory searched for w*.json.

Two records are the same type when kisscrit.search.type_key gives the same key (the key of Census.add and
Census.load; no key is defined here).  That key is a heuristic invariant, so two different types can in principle
be merged into one; the merged census then keeps the record found first.  For each key the merged record is the
copy with the earliest first-found time (ties: lower worker index, then file name); "hits" is the sum over all
copies; "surg" (a surgery call used the type) is the smallest value present.  A record without "sigma" gets the
max-min balanced stress computed from X and the contact pairs (explore_common.kkt_metrics), so every record has the
fields of data/n12_above1.json.

The merge starts from the worker files every time, so running it again on the same files gives the same output;
its result is written to a temporary name and renamed.  A file that cannot be read (for example a truncated copy),
and a file that is itself a merge output, is skipped and listed; the types of the other files are unaffected.
All counts are of numerical, unverified census records.
"""
import argparse, glob, os, sys
from collections import defaultdict
from explore_common import atomic_write_json, kkt_metrics, load_census_file
from kisscrit.search import Census, chi_M0N, euler_check

DEFAULT_LOWER = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'results')


def find_inputs(paths):
    """The worker census files named by paths: files as given, directories searched recursively for w*.json."""
    files = []
    for p in paths:
        files += sorted(glob.glob(os.path.join(p, '**', 'w*.json'), recursive=True)) if os.path.isdir(p) else [p]
    return sorted(set(files))


def _rank(rec, name):
    """Order that decides which copy of a type is the one kept: first-found time, worker index, file name."""
    prov = rec.get('prov') or {}
    return (prov.get('first_found_utc', '~'), prov.get('worker', 1 << 60), name)


def merge_files(files, N=None, d=None):
    """Merge the census files.  Returns (merged Census, info dict); info lists inputs, skipped files and counts."""
    copies = defaultdict(list); inputs = []; skipped = []
    for path in files:
        name = os.path.basename(path)
        try:
            C, meta = load_census_file(path)
        except ValueError as e:
            skipped.append(dict(file=name, reason=str(e)[:200])); continue
        if meta.get('merged'):
            skipped.append(dict(file=name, reason='a merge output, not a worker file (merging it would count hits twice)')); continue
        if N is not None and (C.N, C.d) != (N, d):
            raise ValueError('%s is a census for N=%d d=%d, not N=%d d=%d' % (name, C.N, C.d, N, d))
        N, d = C.N, C.d
        w = meta.get('worker', {})
        for key, rec in C.types.items():
            copies[key].append((_rank(rec, name), rec))
        own = sum(1 for t in C.types.values() if (t.get('prov') or {}).get('worker') == w.get('index') and t['hits'] > 0)
        inputs.append(dict(file=name, worker=w.get('index'), batch=w.get('batch'), types=len(C.types), found=own,
                           keys=list(C.types)))
    if N is None:
        raise ValueError('no readable census file among %d input(s)' % len(files))
    holders = {key: len(cs) for key, cs in copies.items()}   # a key appears at most once per file
    merged = Census(N, d)
    for key in sorted(copies):
        cs = sorted(copies[key], key=lambda c: c[0])
        rec = dict(cs[0][1]); rec['hits'] = sum(r['hits'] for _, r in cs)
        marks = [r['surg'] for _, r in cs if 'surg' in r]
        if marks:
            rec['surg'] = min(marks)
        if 'sigma' not in rec:
            rec['sigma'] = kkt_metrics(rec)['sigma']
        merged.types[key] = rec
    for row in inputs:
        row['only_here'] = sum(1 for k in row['keys'] if holders[k] == 1)
        del row['keys']
    info = dict(inputs=inputs, skipped=skipped, files_read=len(inputs), types_in_files=sum(r['types'] for r in inputs),
                types_merged=len(merged.types))
    return merged, info


def euler_report(merged, lower_dir):
    """Signed count of the merged census against chi(M_{0,N}) (kisscrit.search.euler_check), d = 3 only.

    Types on fewer points come from lower_dir/census_n{n}.json for n < N when those files exist.  The count is
    complete only if the census is; a difference is expected while it is not, and equality does not prove it."""
    if merged.d != 3 or merged.N < 4:
        return dict(status='not applicable (needs d = 3 and N >= 4)')
    cens = {merged.N: merged}; missing = []
    for n in range(3, merged.N):
        p = os.path.join(lower_dir or '', 'census_n%d.json' % n)
        if os.path.exists(p):
            cens[n] = Census.load(p)
        else:
            missing.append(n)
    total, skipped = euler_check(cens, merged.N)
    if missing:
        status = 'partial: lower censuses n=%s not found' % ','.join(map(str, missing))
    else:
        status = 'lower censuses n=3..%d read' % (merged.N - 1)
    return dict(status=status, signed_count=int(total), chi=int(chi_M0N(merged.N)), equal=bool(total == chi_M0N(merged.N)),
                skipped_types=len(skipped))


def euler_line(e):
    """One line for the Euler report e of euler_report."""
    if 'signed_count' not in e:
        return 'Euler check: %s' % e['status']
    return ('Euler check (numerical): signed count %d, chi(M_0,N) = %d, %s; %d types skipped (rattlers or degenerate); %s'
            % (e['signed_count'], e['chi'], 'equal' if e['equal'] else 'not equal', e['skipped_types'], e['status']))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('paths', nargs='+', help='worker census files, or directories searched for w*.json')
    ap.add_argument('--out', required=True, help='merged census to write')
    ap.add_argument('--lower-dir', default=DEFAULT_LOWER, help='directory with census_n{n}.json for the Euler check')
    ap.add_argument('--no-euler', action='store_true', help='skip the Euler check')
    a = ap.parse_args(argv)
    files = find_inputs(a.paths)
    try:
        merged, info = merge_files(files)
    except ValueError as e:
        print('explore_merge: %s' % e, file=sys.stderr); return 1
    info['euler'] = dict(status='skipped (--no-euler)') if a.no_euler else euler_report(merged, a.lower_dir)
    order = sorted(merged.types, key=lambda k: (-merged.types[k]['r'], k))   # decreasing r; the key breaks ties
    atomic_write_json(dict(N=merged.N, d=merged.d, types=[merged.types[k] for k in order], merged=True, merge_info=info), a.out)
    print('merge: %d files read, %d skipped; types in the files (repeats counted) %d -> merged %d'
          % (info['files_read'], len(info['skipped']), info['types_in_files'], info['types_merged']))
    for s in info['skipped']:
        print('  skipped %s: %s' % (s['file'], s['reason']))
    rows = info['inputs']
    if len(rows) > 24:                       # many files: one line per batch
        by = defaultdict(lambda: [0, 0, 0, 0])
        for r in rows:
            b = by[r['batch']]; b[0] += 1; b[1] += r['types']; b[2] += r['found']; b[3] += r['only_here']
        print('  per batch: batch files types found-by-own-worker found-only-there')
        for k in sorted(by, key=lambda x: (x is None, x)):
            print('    %5s %5d %6d %6d %6d' % ((k,) + tuple(by[k])))
    else:
        print('  per worker: file types found-by-own-worker found-only-there')
        for r in rows:
            print('    %s %d %d %d' % (r['file'], r['types'], r['found'], r['only_here']))
    print(euler_line(info['euler']))
    return 0


if __name__ == '__main__':
    sys.exit(main())
