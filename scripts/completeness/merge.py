#!/usr/bin/env python3
"""Summarise the completeness computation from all result files in the repository.

Status: untrusted (exploration or tooling; no part of the proof).

  python scripts/completeness/merge.py --route cond [--of 40]
Prints: candidates generated, certified by status, missing shards, and the codes that need attention
(status 'terminal' = boxes left near a possible critical point, to be identified against the census;
status 'budget' = not resolved within the box budget).  Writes results/completeness/summary-<route>.json.
"""
import argparse, glob, gzip, json, os, collections

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--route', default='cond'); ap.add_argument('--of', type=int)
    a = ap.parse_args()
    cands = set()
    for p in glob.glob(os.path.join(ROOT, 'results/completeness/gen/gen-%s-*.jsonl.gz' % a.route)):
        with gzip.open(p, 'rt') as f:
            for line in f:
                d = json.loads(line)
                if 'code' in d: cands.add(d['code'])
    best = {}
    for p in glob.glob(os.path.join(ROOT, 'results/completeness/cert/cert-%s-*.jsonl' % a.route)):
        for line in open(p):
            d = json.loads(line); c = d['code']
            rank = {'excluded': 2, 'terminal': 2, 'budget': 1}[d['status']]
            if c not in best or (rank, d['budget']) > ({'excluded': 2, 'terminal': 2, 'budget': 1}[best[c]['status']], best[c]['budget']):
                best[c] = d
    st = collections.Counter(d['status'] for d in best.values())
    secs = sum(d['seconds'] for d in best.values())
    missing = len(cands - set(best))
    print('candidates generated: %d' % len(cands))
    print('certified: %d  %s   (CPU time %.1f h)' % (len(best), dict(st), secs / 3600))
    print('not yet run: %d' % missing)
    att = [dict(code=c, status=d['status'], E=d['E'], sizes=d['sizes'], n_terminal=d.get('n_terminal', 0), boxes=d['stats']['boxes'])
           for c, d in best.items() if d['status'] != 'excluded']
    json.dump(dict(route=a.route, generated=len(cands), certified=len(best), status=dict(st), missing=missing,
                   cpu_hours=secs / 3600, attention=att),
              open(os.path.join(ROOT, 'results/completeness/summary-%s.json' % a.route), 'w'), indent=1)
    print('needs attention: %d (listed in results/completeness/summary-%s.json)' % (len(att), a.route))


if __name__ == '__main__':
    main()
