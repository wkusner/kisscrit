#!/usr/bin/env python3
"""Differential run of the C checker kernel against the Python reference checker (checklist F4, spec v0.3
section 4.3 "Differential").

Status: untrusted (a test driver; it reads certificates and writes a report outside the repository's results).

Unit of comparison: one leaf or one contraction slab of one chart of one certificate (the nodes of
kisscrit/core/certificate.py:replay, as in tests/test_creplay.py).  For each unit it compares:
  1. the verdict of the kernel's chart driver (kc_check_chart through kisscrit/ckernel.py) with the verdict of
     the unchanged reference, kisscrit/core/checker.py:check_certificate run on that chart alone: a leaf or slab
     is rejected by the reference iff it appears in the reference's `failed` list (matched by kind and the IEEE
     bit patterns of its box), and accepted otherwise;
  2. the same leaf test on the same inputs, test by test, as IEEE bit patterns or rule names:
       untightened invalid/vsum   ckernel.CGraph.invalid_vsum(tighten=0)  vs  Checker.test_invalid_vsum on
                                  G.corners(lo, hi, tighten=False);
       tightened invalid/vsum     (charts with a pentagon) tighten=1      vs  Checker.evaluate + test_invalid_vsum;
       the hint's interval s      CGraph.context + CGraph.hint            vs  Checker.linear_context +
                                  checker.combination(sorted(mu.items()), Fm, JB, Y) with test_linear_mu's Y,
                                  on B = the contraction box (slabs) or the leaf's own box (linear leaves);
  3. per certificate, ckernel.check_certificate against the reference's results composed over the charts:
     verified, the failed boxes, the terminal boxes, open and leaf counts.
Every exception of the reference is recorded as 'error' and compared with the kernel's error codes.

Counts reported (spec v0.3 section 8, prompt "DIFFERENTIAL DESIGN"): both accept (the kernel closes, the
reference accepts); both reject (open leaves, terminal leaves wider than tol); kernel unclosed (fallback), split
by the reference's verdict; kernel accepts where the reference rejects (must be 0; any such unit is listed in
full).

Examples:
  python3 scripts/completeness/diff_kernel.py --sample 500 --seed 20261021 --workers 2 --out /tmp/x/diff-a
  python3 scripts/completeness/diff_kernel.py --files <scratch>/results/.../certificates-cond-174-of-200.jsonl --out ...
The output directory gets units.jsonl (one line per certificate: counts, fallback units with the reference's
verdict, mismatches) and summary.json.  Never point --out into the repository's results/.
"""
import argparse, collections, glob, json, os, random, struct, sys, time
from multiprocessing import Pool

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
CERTS = os.path.join(ROOT, 'results', 'completeness', 'certificates')
LIBPATH = None          # set in workers from --lib (a mutant build); None = build/libkisscheck.so


def _bits(xs):
    """IEEE bit patterns of a sequence of doubles (so 0.0 and -0.0 differ)."""
    return struct.pack('<%dd' % len(xs), *[float(x) for x in xs])


def _ibits(lo, hi):
    """Bit pattern of an interval."""
    return struct.pack('<dd', float(lo), float(hi)).hex()


RULE = {None: 0, 'invalid': 1, 'vsum': 2}


def _py_rule(K, G, lo, hi, tight):
    """The reference's untightened or tightened invalid/vsum on [lo, hi]: 0, 1, 2 or 'error'."""
    try:
        C, ok = (K.evaluate(lo, hi) if tight else G.corners(lo, hi, tighten=False))
        return RULE[K.test_invalid_vsum(C, ok)]
    except Exception as e:
        return 'error'


def _py_hint_s(K, G, B, lo, hi, mu, ctx_cache):
    """The reference's interval s of test_linear_mu (or 'error'), with the context of B cached per box."""
    from kisscrit.core import checker as ck
    from kisscrit.core import interval as I
    from kisscrit.core.interval import Iv
    try:
        key = (_bits(B[0]), _bits(B[1]))
        if key not in ctx_cache:
            CB, _ = G.corners(B[0], B[1], tighten=False)
            ctx_cache.clear(); ctx_cache[key] = K.linear_context(B[0], B[1], CB)
        _, _, m, Fm, JB = ctx_cache[key]
        Y = [I.hull(Iv(lo[j]) - Iv(m[j]), Iv(hi[j]) - Iv(m[j])) for j in range(G.npar)]
        s = ck.combination(sorted(mu.items()), Fm, JB, Y)
        return _ibits(s.lo, s.hi), (s.lo > 0 or s.hi < 0)
    except Exception:
        return 'error', None


def _c_hint_s(cg, B, lo, hi, mu, cctx):
    """The kernel's interval s (or 'error') for the same inputs."""
    key = (_bits(B[0]), _bits(B[1]))
    if cctx.get('key') != key:
        rc = cg.context(B[0], B[1]); cctx['key'] = key; cctx['rc'] = rc
    if cctx['rc'] != 0: return 'error', None
    rows = list(mu.keys()); vals = list(mu.values())
    rc, closed, s = cg.hint(lo, hi, (rows, vals))
    if rc != 0: return 'error', None
    return _ibits(*s), bool(closed)


def check_one(line):
    """Compare kernel and reference on every unit of one certificate line.  Returns a dict of counts and lists."""
    from kisscrit import ckernel as kk
    from kisscrit.core import checker as ck, model as am
    from kisscrit.core.certificate import replay
    d = json.loads(line); cert = d['certificate']
    out = dict(code=d['code'], charts=0, units=0, kinds=collections.Counter(), both_accept=0, both_reject=0,
               terminal_ok=0, unclosed=0, unclosed_py_accept=0, unclosed_py_reject=0, c_accept_py_reject=0,
               test_mismatch=collections.Counter(), tests=collections.Counter(), c_rules=collections.Counter(),
               fallback=[], violations=[], mismatches=[], chart_errors=[], whole=None,
               c_seconds=0.0, py_seconds=0.0)
    composed = dict(verified=True, failed=[], terminal=[], open=0, leaves=0)
    t_cert = float(cert['tol'])
    for ci, ch in enumerate(cert['charts']):
        if ch.get('not_run'):
            composed['verified'] = False; continue
        out['charts'] += 1
        G = am.Graph(cert['faces'], ch['modes']); K = ck.Checker(G)
        cg = kk.CGraph(cert['faces'], ch['modes'], LIBPATH)
        res = kk.check_chart(cg, ch, t_cert)
        out['c_seconds'] += res['seconds']
        # the reference on this chart alone
        t0 = time.perf_counter()
        try:
            ref = ck.check_certificate(dict(cert, charts=[ch]), structural=False); ref_err = None
        except Exception as e:
            ref = None; ref_err = repr(e)
        out['py_seconds'] += time.perf_counter() - t0
        if ref is None:
            composed['verified'] = False
        else:
            for k in ('failed', 'terminal'): composed[k] += ref[k]
            composed['open'] += ref['open']; composed['leaves'] += ref['leaves']
            if not ref['verified']: composed['verified'] = False
        failed = collections.Counter()
        if ref is not None:
            for f in ref['failed']: failed[(f['kind'], _bits(f['lo']), _bits(f['hi']))] += 1
        if res['rc'] != 0:
            out['chart_errors'].append(dict(chart=ci, c=res['error'], model_err=res['model_err'], py=ref_err))
        verdict = res['verdict']; u = 0
        ctx_py, ctx_c = {}, {}
        has_pent = any(k == 5 for k in G.sizes)
        try:
            nodes = list(replay(ch))
        except Exception as e:
            nodes = []; out['chart_errors'].append(dict(chart=ci, replay=repr(e)))
        for kind, lo, hi, extra in nodes:
            lo = [float(x) for x in lo]; hi = [float(x) for x in hi]
            if kind == 'split': continue
            if kind == 'contract':
                removed, _ = extra
                units = [('slab', [float(x) for x in s[0]], [float(x) for x in s[1]], s[6], (lo, hi)) for s in removed]
            else:
                units = [(kind, lo, hi, extra, (lo, hi))]
            for ukind, ulo, uhi, mu, B in units:
                if u >= len(verdict):           # the kernel stopped on an error before this unit
                    u += 1; continue
                v = int(verdict[u]); u += 1
                out['units'] += 1; out['kinds'][ukind] += 1; out['c_rules'][kk.VERDICT[v]] += 1
                key = (ukind, _bits(ulo), _bits(uhi))
                py_reject = ref is None or failed[key] > 0
                if ref is not None and failed[key] > 0: failed[key] -= 1
                if ukind == 'terminal':
                    if v == 8 and py_reject: out['both_reject'] += 1
                    elif v == 7 and not py_reject: out['terminal_ok'] += 1
                    else: out['mismatches'].append(dict(chart=ci, unit=u - 1, kind=ukind, c=v, py_reject=py_reject))
                    continue
                if ukind == 'open':
                    if v == 9: out['both_reject'] += 1
                    else: out['mismatches'].append(dict(chart=ci, unit=u - 1, kind=ukind, c=v))
                    continue
                if v in kk.CLOSED:
                    if py_reject:
                        out['c_accept_py_reject'] += 1
                        out['violations'].append(dict(code=d['code'], chart=ci, unit=u - 1, kind=ukind, lo=ulo,
                                                      hi=uhi, c=kk.VERDICT[v], py='error' if ref is None else 'failed',
                                                      ref_error=ref_err))
                    else:
                        out['both_accept'] += 1
                else:
                    out['unclosed'] += 1
                    out['unclosed_py_reject' if py_reject else 'unclosed_py_accept'] += 1
                    out['fallback'].append(dict(chart=ci, unit=u - 1, kind=ukind,
                                                py='reject' if py_reject else 'accept'))
                # the same tests on the same inputs
                for tight in ((0, 1) if has_pent else (0,)):
                    rc, rule = cg.invalid_vsum(ulo, uhi, tight)
                    c_r = rule if rc == 0 else 'error'
                    p_r = _py_rule(K, G, ulo, uhi, tight)
                    name = 'invalid_vsum_tight' if tight else 'invalid_vsum'
                    out['tests'][name] += 1
                    if c_r != p_r:
                        out['test_mismatch'][name] += 1
                        out['mismatches'].append(dict(chart=ci, unit=u - 1, test=name, c=c_r, py=p_r))
                if mu and (ukind == 'slab' or ukind == 'linear'):
                    pb, pc = _py_hint_s(K, G, B, ulo, uhi, mu, ctx_py)
                    cb, cc = _c_hint_s(cg, B, ulo, uhi, mu, ctx_c)
                    out['tests']['hint_s'] += 1
                    if pb != cb or pc != cc:
                        out['test_mismatch']['hint_s'] += 1
                        out['mismatches'].append(dict(chart=ci, unit=u - 1, test='hint_s', c=cb, py=pb))
        if res['rc'] == 0 and u != len(verdict):
            out['mismatches'].append(dict(chart=ci, units_c=len(verdict), units_py=u))
    # whole certificate: the front end (C + fallback) against the composed reference
    try:
        w = kk.check_certificate(cert, structural=False, path=LIBPATH)
        key = lambda fs: sorted((f['kind'], _bits(f['lo']), _bits(f['hi'])) for f in fs)
        out['whole'] = dict(same=(w['verified'] == composed['verified'] and key(w['failed']) == key(composed['failed'])
                                  and len(w['terminal']) == len(composed['terminal']) and w['open'] == composed['open']
                                  and w['leaves'] == composed['leaves']),
                            verified=w['verified'], fallback=w['kernel']['fallback'],
                            fallback_closed=w['kernel']['fallback_closed'])
    except Exception as e:
        out['whole'] = dict(same=False, error=repr(e))
    for k in ('kinds', 'test_mismatch', 'tests', 'c_rules'): out[k] = dict(out[k])
    return out


def _init(libpath):
    """Pool initializer: the library to test."""
    global LIBPATH
    LIBPATH = libpath


def sample_lines(n, seed):
    """n certificate lines drawn without replacement with random.Random(seed) from all lines of the `cond`
    certificate files (files in sorted order, then line order), as tests/test_creplay.py draws its sample."""
    index = []
    for p in sorted(glob.glob(os.path.join(CERTS, 'certificates-cond-*-of-40.jsonl'))):
        with open(p) as f:
            for i, _ in enumerate(f): index.append((p, i))
    pick = sorted(random.Random(seed).sample(range(len(index)), n))
    want = collections.defaultdict(set)
    for j in pick: want[index[j][0]].add(index[j][1])
    out = []
    for p in sorted(want):
        with open(p) as f:
            for i, line in enumerate(f):
                if i in want[p]: out.append(line)
    return out, len(index)


def summarize(rows):
    """Totals over certificates."""
    S = collections.Counter(); kinds = collections.Counter(); tests = collections.Counter()
    tm = collections.Counter(); rules = collections.Counter(); fb = collections.Counter()
    whole_same = whole_n = verified = 0; viol = []; chart_err = []; nmis = 0
    for r in rows:
        for k in ('charts', 'units', 'both_accept', 'both_reject', 'terminal_ok', 'unclosed', 'unclosed_py_accept',
                  'unclosed_py_reject', 'c_accept_py_reject'):
            S[k] += r[k]
        S['c_seconds'] += r['c_seconds']; S['py_seconds'] += r['py_seconds']
        kinds.update(r['kinds']); tests.update(r['tests']); tm.update(r['test_mismatch']); rules.update(r['c_rules'])
        for f in r['fallback']: fb[(f['kind'], f['py'])] += 1
        viol += r['violations']; chart_err += [dict(code=r['code'], **e) for e in r['chart_errors']]
        nmis += len(r['mismatches'])
        if r['whole'] is not None:
            whole_n += 1; whole_same += bool(r['whole'].get('same')); verified += bool(r['whole'].get('verified'))
    return dict(certificates=len(rows), **{k: S[k] for k in S}, units_by_kind=dict(kinds), c_verdicts=dict(rules),
                tests=dict(tests), test_mismatches=dict(tm), mismatch_records=nmis,
                fallback_by_kind_and_reference=[dict(kind=k, reference=p, n=n) for (k, p), n in sorted(fb.items())],
                whole_certificate_same=whole_same, whole_certificate_compared=whole_n,
                whole_certificate_verified_by_front_end=verified, violations=viol, chart_errors=chart_err)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sample', type=int, help='number of `cond` certificates to draw')
    ap.add_argument('--seed', type=int, default=20261021)
    ap.add_argument('--files', nargs='*', help='certificate files to compare in full instead of a sample')
    ap.add_argument('--workers', type=int, default=os.cpu_count())
    ap.add_argument('--lib', default=None, help='a different build of libkisscheck.so (mutation tests)')
    ap.add_argument('--out', required=True, help='output directory (outside the repository results/)')
    a = ap.parse_args()
    if os.path.abspath(a.out).startswith(os.path.join(ROOT, 'results')):
        sys.exit('refusing to write under the repository results/')
    os.makedirs(a.out, exist_ok=True)
    if a.files:
        lines = [l for p in a.files for l in open(p)]; total = None
    else:
        lines, total = sample_lines(a.sample, a.seed)
    t0 = time.time(); rows = []
    with open(os.path.join(a.out, 'units.jsonl'), 'w') as f, Pool(a.workers, initializer=_init, initargs=(a.lib,)) as pool:
        for i, r in enumerate(pool.imap_unordered(check_one, lines, chunksize=1)):
            f.write(json.dumps(r) + '\n'); f.flush(); rows.append(r)
            if (i + 1) % 50 == 0:
                print('%d/%d certificates, %.0fs' % (i + 1, len(lines), time.time() - t0), flush=True)
    s = summarize(rows)
    s.update(seed=None if a.files else a.seed, population=total, files=a.files, lib=a.lib, wall_seconds=round(time.time() - t0, 1))
    json.dump(s, open(os.path.join(a.out, 'summary.json'), 'w'), indent=1)
    print(json.dumps({k: v for k, v in s.items() if k not in ('violations',)}, indent=1))
    if s['violations']: print('VIOLATIONS: %d (kernel accepts where the reference rejects)' % len(s['violations']))


if __name__ == '__main__':
    main()
