"""Export a merged census for the refine stage: refine-input.json, refine-queue.tsv and summary.txt.

Status: untrusted (exploration or tooling; no part of the proof).

usage: python scripts/explore/explore_export.py MERGED.json --out-dir DIR [--known CENSUS.json ...]

* refine-input.json has the layout of data/n12_above1.json: N, d, note, and a list "types" whose records carry
  exactly the fields r, c, k, coindex, kernel, sdim, tmin, sym, sym_proper, pairs, sigma, X (sym and sym_proper
  are null for a type with a rattler, where the census does not compute them).  Every type of the merged census is
  in it, in the order of decreasing r; the queue says which ones to refine.
* refine-queue.tsv has one row per type: its position `index` in refine-input.json, `priority`, `flag`, `new`,
  hits, co-index, kernel, sdim (dimension of the space of self-stresses), rattlers, r, the number of contacts,
  the KKT residual, the minimum non-contact gap, the stress minimum, the symmetry order, and the provenance of
  the first find.  The diagnostics are recomputed from X
  and the contact pairs in double precision (explore_common.kkt_metrics).
  Flags: regular (kernel 0, a strictly positive max-min stress, no rattlers: the types the refine stage solves),
  degenerate and rattler (listed, never solved here), suspect (KKT residual or non-contact gap beyond the
  thresholds of explore_common; checked first).  `new` is yes when the type's key is in none of the --known census
  files (with none given, every type is new).  Priority: new types first, then fewer hits, then lower co-index,
  then larger r.
* summary.txt counts the flags and repeats the Euler line of the merge.
Every value is numerical (double precision); nothing is verified and no value is a proof input.
"""
import argparse, os, sys
from collections import Counter
from explore_common import (HEADER_LABEL, REFINE_FIELDS, KKT_RESID_MAX, GAP_MIN, STRESS_MIN, atomic_write_json,
                            atomic_write_text, classify, kkt_metrics, load_census_file)

TSV_COLUMNS = ('index', 'priority', 'flag', 'new', 'hits', 'coindex', 'kernel', 'sdim', 'rattlers', 'r', 'n_contacts',
               'kkt_resid', 'min_gap', 'stress_min', 'sym', 'first_found_utc', 'run_id', 'commit', 'worker')


def _fmt(v, spec):
    """v formatted with spec, or NA for a missing value."""
    return 'NA' if v is None else format(v, spec)


def build(merged_path, known_paths):
    """Rows of the queue and records of refine-input.json for the merged census at merged_path."""
    C, meta = load_census_file(merged_path)
    known = set()
    for p in known_paths:
        known |= set(load_census_file(p)[0].types)
    order = sorted(C.types, key=lambda k: (-C.types[k]['r'], k))
    rows, records = [], []
    for index, key in enumerate(order):
        rec = C.types[key]; m = kkt_metrics(rec); prov = rec.get('prov') or {}
        sigma = rec['sigma'] if rec.get('sigma') is not None else m['sigma']
        records.append({f: (sigma if f == 'sigma' else rec.get(f)) for f in REFINE_FIELDS})
        rows.append(dict(index=index, flag=classify(rec, m), new=key not in known, hits=rec.get('hits', 0), coindex=rec['coindex'],
                         kernel=rec['kernel'], sdim=rec['sdim'], rattlers=rec.get('rattlers', 0), r=rec['r'], n_contacts=m['n_contacts'],
                         kkt_resid=m['kkt_resid'], min_gap=m['min_gap'], stress_min=m['stress_min'], sym=rec.get('sym'),
                         first_found_utc=prov.get('first_found_utc'), run_id=prov.get('run_id'), commit=prov.get('commit'),
                         worker=prov.get('worker')))
    for rank, row in enumerate(sorted(rows, key=lambda r: (not r['new'], r['hits'], r['coindex'], -r['r'], r['index'])), 1):
        row['priority'] = rank
    return C, meta, rows, records


def _provenance(C):
    """Text naming the runs and kisscrit commits of the records, for headers."""
    prov = [t.get('prov') or {} for t in C.types.values()]
    runs = sorted({p['run_id'] for p in prov if p.get('run_id')}); commits = sorted({p['commit'] for p in prov if p.get('commit')})
    return 'run %s, kisscrit commit %s' % (', '.join(runs) or 'unknown', ', '.join(commits) or 'unknown')


def queue_text(C, rows, n_known, where):
    """Text of refine-queue.tsv: three comment lines, the column names, then the rows in priority order."""
    lines = ['# refine queue of %d points on S^%d (%s)' % (C.N, C.d - 1, where), '# ' + HEADER_LABEL,
             '# index = position in refine-input.json; new = key in none of the known censuses (%d given); '
             'priority: new first, fewer hits, lower coindex, larger r' % n_known, '\t'.join(TSV_COLUMNS)]
    for r in sorted(rows, key=lambda r: r['priority']):
        cells = [r['index'], r['priority'], r['flag'], 'yes' if r['new'] else 'no', r['hits'], r['coindex'], r['kernel'],
                 r['sdim'], r['rattlers'], repr(r['r']), r['n_contacts'], _fmt(r['kkt_resid'], '.3e'), _fmt(r['min_gap'], '.3e'),
                 _fmt(r['stress_min'], '.3e'), _fmt(r['sym'], 'd'), r['first_found_utc'] or 'NA', r['run_id'] or 'NA',
                 r['commit'] or 'NA', _fmt(r['worker'], 'd')]
        lines.append('\t'.join(map(str, cells)))
    return '\n'.join(lines) + '\n'


def summary_text(C, meta, rows, where):
    """Text of summary.txt: the counts per flag, the Euler line of the merge, the thresholds and the label."""
    flags = Counter(r['flag'] for r in rows)
    e = (meta.get('merge_info') or {}).get('euler') or {}
    if 'signed_count' in e:
        euler = ('Euler check (numerical): signed count %d, chi(M_0,N) = %d, %s; %d types skipped; %s'
                 % (e['signed_count'], e['chi'], 'equal' if e['equal'] else 'not equal', e['skipped_types'], e['status']))
    else:
        euler = 'Euler check (numerical): not available'
    return '\n'.join([
        'exploration summary: N=%d d=%d, %s' % (C.N, C.d, where),
        'types: %d (new: %d)' % (len(rows), sum(r['new'] for r in rows)),
        'flags: regular %d, degenerate %d, rattler %d, suspect %d' % tuple(flags[f] for f in ('regular', 'degenerate', 'rattler', 'suspect')),
        euler,
        'thresholds: KKT residual > %g or minimum non-contact gap < %g = suspect; stress <= %g = not strictly positive'
        % (KKT_RESID_MAX, GAP_MIN, STRESS_MIN),
        'label: every value is numerical (double precision); nothing is verified; no value is a proof input']) + '\n'


def write_outputs(C, meta, rows, records, known_paths, out_dir):
    """Write refine-input.json, refine-queue.tsv and summary.txt into out_dir (created if needed); returns the summary."""
    os.makedirs(out_dir, exist_ok=True)
    where = _provenance(C)
    note = ('Critical configurations of %d points on S^%d found by the exploration driver (%s): %d types. %s.'
            % (C.N, C.d - 1, where, len(records), HEADER_LABEL))
    atomic_write_json(dict(N=C.N, d=C.d, note=note, types=records), os.path.join(out_dir, 'refine-input.json'))
    atomic_write_text(queue_text(C, rows, len(known_paths), where), os.path.join(out_dir, 'refine-queue.tsv'))
    summary = summary_text(C, meta, rows, where)
    atomic_write_text(summary, os.path.join(out_dir, 'summary.txt'))
    return summary


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('merged', help='merged census written by explore_merge.py')
    ap.add_argument('--out-dir', required=True, help='directory for refine-input.json, refine-queue.tsv, summary.txt')
    ap.add_argument('--known', action='append', default=[], help='a census file of already known types (repeatable)')
    a = ap.parse_args(argv)
    try:
        C, meta, rows, records = build(a.merged, a.known)
    except (ValueError, OSError) as e:
        print('explore_export: %s' % e, file=sys.stderr); return 1
    print(write_outputs(C, meta, rows, records, a.known, a.out_dir), end='')
    return 0


if __name__ == '__main__':
    sys.exit(main())
