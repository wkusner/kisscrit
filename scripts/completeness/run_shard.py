#!/usr/bin/env python3
"""Sharded, resumable driver for the completeness computation.

Phase A (gen):      enumerate candidates from triangulations, shard = triangulation indices i with i % K == k.
                    Output: results/completeness/gen/gen-<route>-<k>-of-<K>.jsonl.gz  (one candidate per line)
Phase B (certify):  certify deduplicated candidates, shard = candidates whose code hashes to k mod K.
                    Output: results/completeness/cert/cert-<route>-<k>-of-<K>.jsonl  (one result per line)

Both phases append and skip work already recorded, so a killed run is resumed by running the same command.
Routes: 'cond' = triangulations 0, 1, 5, 7 only (Conjecture 7.3(a), icosahedral-type hulls); 'full' = all 7595.

Examples:
  python scripts/completeness/run_shard.py gen --route cond --shard 0 --of 4
  python scripts/completeness/run_shard.py certify --route cond --shard 3 --of 40 --budget 300000 --workers 8
"""
import argparse, gzip, hashlib, json, os, sys, time, glob
from multiprocessing import Pool

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
GEN = os.path.join(ROOT, 'results', 'completeness', 'gen')
CERT = os.path.join(ROOT, 'results', 'completeness', 'cert')
ROUTES = {'cond': [0, 1, 5, 7], 'test': [3000]}
R_MIN_DEFAULT = 1.0001          # certify configurations with r >= R_MIN (theta >= theta(R_MIN))


def theta_of_r(r):
    import numpy as np
    return 2 * np.arcsin(r / (1 + r))


def tri_indices(route):
    if route in ROUTES: return ROUTES[route]
    T = json.load(open(os.path.join(ROOT, 'data', 'completeness', 'triangulations12.json')))
    return list(range(len(T)))


def shard_of(code, K):
    return int(hashlib.sha1(code.encode()).hexdigest(), 16) % K


# ------------------------------------------------------------------ phase A
def gen(args):
    from kisscrit.completeness.candidates import candidates_of_triangulation
    T = json.load(open(os.path.join(ROOT, 'data', 'completeness', 'triangulations12.json')))
    idx = [i for i in tri_indices(args.route) if i % args.of == args.shard]
    os.makedirs(GEN, exist_ok=True)
    out = os.path.join(GEN, 'gen-%s-%d-of-%d.jsonl.gz' % (args.route, args.shard, args.of))
    done = set(); seen = set()
    if os.path.exists(out):
        with gzip.open(out, 'rt') as f:
            for line in f:
                d = json.loads(line)
                if 'done_triangulation' in d: done.add(d['done_triangulation'])
                else: seen.add(d['code'])
    t0 = time.time()
    with gzip.open(out, 'at') as f:
        for i in idx:
            if i in done: continue
            n = 0
            for code, cyc in candidates_of_triangulation(T[i], seen):
                f.write(json.dumps(dict(code=code, faces=cyc, tri=i)) + '\n'); n += 1
            f.write(json.dumps(dict(done_triangulation=i, new=n)) + '\n'); f.flush()
            print('triangulation %d: %d new candidates (%d total, %.0fs)' % (i, n, len(seen), time.time() - t0), flush=True)
    print('gen shard done:', out)


# ------------------------------------------------------------------ phase B
def load_candidates(route):
    cands = {}
    for p in sorted(glob.glob(os.path.join(GEN, 'gen-%s-*.jsonl.gz' % route))):
        with gzip.open(p, 'rt') as f:
            for line in f:
                d = json.loads(line)
                if 'code' in d: cands.setdefault(d['code'], d['faces'])
    return cands


_ARGS = None


def _init(a):
    global _ARGS; _ARGS = a


def certify_one(item):
    import numpy as np
    from kisscrit.completeness.structure import Structure
    from kisscrit.completeness import runner
    code, faces = item; a = _ARGS
    S = Structure(faces)
    t0 = time.time()
    term, st = runner.run(S, theta_of_r(a.r_min), np.radians(63.4349489), tol=a.tol, max_boxes=a.budget)
    status = 'budget' if st['aborted'] else ('excluded' if not term else 'terminal')
    rec = dict(code=code, status=status, E=len(S.edges), sizes=sorted(S.sizes), stats=st,
               r_min=a.r_min, tol=a.tol, budget=a.budget, host=os.uname().nodename, seconds=round(time.time() - t0, 2))
    if term:
        rec['terminal'] = [dict(chart=ch, lo=list(map(float, lo)), hi=list(map(float, hi))) for ch, lo, hi in term[:200]]
        rec['n_terminal'] = len(term)
    rec['faces'] = faces
    return rec


def certify(args):
    cands = load_candidates(args.route)
    if not cands: sys.exit('no candidates for route %s: run the gen phase first (and git pull)' % args.route)
    mine = sorted(c for c in cands if shard_of(c, args.of) == args.shard)
    os.makedirs(CERT, exist_ok=True)
    out = os.path.join(CERT, 'cert-%s-%d-of-%d.jsonl' % (args.route, args.shard, args.of))
    done = set()
    if os.path.exists(out):
        for line in open(out):
            d = json.loads(line)
            if d.get('budget', 0) >= args.budget or d['status'] != 'budget': done.add(d['code'])
    todo = [(c, cands[c]) for c in mine if c not in done]
    print('shard %d/%d: %d candidates, %d done, %d to do' % (args.shard, args.of, len(mine), len(done), len(todo)), flush=True)
    t0 = time.time(); n = 0
    with open(out, 'a') as f, Pool(args.workers, initializer=_init, initargs=(args,)) as pool:
        for rec in pool.imap_unordered(certify_one, todo, chunksize=1):
            f.write(json.dumps(rec) + '\n'); f.flush(); n += 1
            if n % 20 == 0 or rec['status'] != 'excluded':
                print('%d/%d  %s  %s  %.1fs  (elapsed %.0fs)' % (n, len(todo), rec['status'], rec['code'][:24], rec['seconds'],
                                                              time.time() - t0), flush=True)
    print('certify shard done:', out)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('phase', choices=['gen', 'certify'])
    ap.add_argument('--route', default='cond', choices=['cond', 'full', 'test'])
    ap.add_argument('--shard', type=int, required=True)
    ap.add_argument('--of', type=int, required=True)
    ap.add_argument('--workers', type=int, default=os.cpu_count())
    ap.add_argument('--budget', type=int, default=300000, help='max boxes per candidate')
    ap.add_argument('--tol', type=float, default=1e-6, help='terminal box width (radians)')
    ap.add_argument('--r-min', type=float, default=R_MIN_DEFAULT, dest='r_min')
    a = ap.parse_args()
    gen(a) if a.phase == 'gen' else certify(a)
