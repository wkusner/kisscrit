"""The C kernel's leaf tests (csrc/kc_leaf.c, checklist F4) against the Python reference (kisscrit/core/checker.py).

Status: untrusted (a test).  Spec v0.3 section 4.3: 'The C kernel and the Python reference agree on ... the
cond certificates ... A disagreement of the form "C accepts, Python rejects" stops the work until it is
explained.'  The full differential run is scripts/completeness/diff_kernel.py; this test runs its comparison
(`check_one`) on a small sample so that it stays in the suite.

What is compared, per unit (leaf or contraction slab) and per certificate: see diff_kernel.py's docstring.  In
short: the kernel's verdict against the reference's (the kernel may leave a unit unclosed, never close one the
reference rejects); untightened and tightened invalid/vsum by rule; the hint's interval s by IEEE bit pattern;
and ckernel.check_certificate (C plus fallback) against the reference on the whole certificate.

Inputs:
  * test_differential_on_cond_sample: KISSCRIT_CKERNEL_CERTS (default 120) certificates drawn with seed
    KISSCRIT_CKERNEL_SEED (default 20261023) from all `cond` certificate lines, keeping those with at most
    KISSCRIT_CKERNEL_MAXEV (default 3000) events in all (to bound the suite's time; diff_kernel.py has no bound);
  * test_hint_edge_cases: hand-made hints (empty, a repeated row, zero, negative, NaN and infinite values, a row
    out of range) on real charts, comparing kc_leaf_hint's s with checker.combination by bit pattern;
  * test_errors_match_reference: malformed charts and a hint row out of range: the front end raises exactly
    when the reference raises;
  * test_terminal_and_open, test_limits: terminal width at the tolerance, open leaves, the depth and piece
    limits;
  * test_power: compiled mutants of kc_leaf.c (rounding direction, comparisons, summation order, the
    tightening flag, the validity flag), each of which the comparison must detect.
Skipped when build/libkisscheck.so is not built (run `make -C csrc`); test_power also needs a C compiler.
"""
import base64, collections, glob, importlib.util, json, math, os, pathlib, random, shutil, struct, subprocess, zlib
import numpy as np
import pytest
from kisscrit.core import carith, checker as ck, model as am
from kisscrit.core import interval as I
from kisscrit.core.interval import Iv
from kisscrit.core.certificate import replay, unpack

pytestmark = pytest.mark.skipif(not carith.available(), reason='C kernel not built (make -C csrc)')

ROOT = pathlib.Path(__file__).resolve().parents[1]
CERTS = ROOT / 'results' / 'completeness' / 'certificates'
NCERTS = int(os.environ.get('KISSCRIT_CKERNEL_CERTS', '120'))
SEED = int(os.environ.get('KISSCRIT_CKERNEL_SEED', '20261023'))
MAXEV = int(os.environ.get('KISSCRIT_CKERNEL_MAXEV', '3000'))
CFLAGS = ['-O2', '-std=c11', '-Wall', '-Wextra', '-fno-fast-math', '-ffp-contract=off', '-fPIC']   # csrc/Makefile


def _diff():
    """scripts/completeness/diff_kernel.py as a module."""
    spec = importlib.util.spec_from_file_location('diff_kernel', ROOT / 'scripts' / 'completeness' / 'diff_kernel.py')
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def _kk():
    from kisscrit import ckernel
    return ckernel


def _nev(line):
    """Total number of events over the run charts of a certificate line."""
    c = json.loads(line)['certificate']
    return sum(ch.get('n_events', 0) for ch in c['charts'] if not ch.get('not_run'))


def _sample(n, seed, maxev):
    lines = []
    for p in sorted(glob.glob(str(CERTS / 'certificates-cond-*-of-40.jsonl'))):
        with open(p) as f: lines += f.readlines()
    rng = random.Random(seed)
    out = []
    for j in rng.sample(range(len(lines)), min(len(lines), 20 * n)):
        if _nev(lines[j]) <= maxev: out.append(lines[j])
        if len(out) == n: break
    return out


def _problems(r):
    """Everything in a check_one result that counts as a disagreement."""
    p = []
    if r['c_accept_py_reject']: p.append(('C accepts, reference rejects', r['violations'][:2]))
    if r['test_mismatch']: p.append(('test mismatches', r['test_mismatch'], r['mismatches'][:2]))
    elif r['mismatches']: p.append(('mismatches', r['mismatches'][:2]))
    if not (r['whole'] or {}).get('same'): p.append(('whole certificate', r['whole']))
    if r['chart_errors']: p.append(('chart errors', r['chart_errors'][:2]))
    return p


@pytest.fixture(scope='module')
def sample():
    if not CERTS.exists(): pytest.skip('no certificates')
    return _sample(NCERTS, SEED, MAXEV)


def test_differential_on_cond_sample(sample):
    D = _diff(); D.LIBPATH = None
    tot = collections.Counter(); bad = []
    for line in sample:
        r = D.check_one(line)
        for k in ('units', 'both_accept', 'unclosed', 'c_accept_py_reject'): tot[k] += r[k]
        for k, v in r['tests'].items(): tot['test ' + k] += v
        if _problems(r): bad.append((r['code'][:30], _problems(r)))
    assert not bad, bad[:3]
    assert tot['units'] > 1000 and tot['test hint_s'] > 500, tot          # the sample exercises the tests
    assert tot['both_accept'] + tot['unclosed'] <= tot['units']


# ------------------------------------------------------------------------------------------- hand-made cases

def _first_chart_with_slabs(lines):
    """A certificate and chart whose replay has a contraction with slabs, with the first such slab and box."""
    for line in lines:
        cert = json.loads(line)['certificate']
        for ch in cert['charts']:
            if ch.get('not_run'): continue
            for kind, lo, hi, extra in replay(ch):
                if kind == 'contract' and extra[0]:
                    s = extra[0][0]
                    return cert, ch, [float(x) for x in lo], [float(x) for x in hi], \
                        [float(x) for x in s[0]], [float(x) for x in s[1]]
    pytest.skip('no slab found')


def _py_s(G, K, B, lo, hi, mu):
    CB, _ = G.corners(B[0], B[1], tighten=False)
    _, _, m, Fm, JB = K.linear_context(B[0], B[1], CB)
    Y = [I.hull(Iv(lo[j]) - Iv(m[j]), Iv(hi[j]) - Iv(m[j])) for j in range(G.npar)]
    return ck.combination(sorted(mu.items()), Fm, JB, Y)


def _f32(x):
    return float(np.float32(x))


def test_hint_edge_cases(sample):
    kk = _kk()
    cert, ch, Blo, Bhi, slo, shi = _first_chart_with_slabs(sample)
    G = am.Graph(cert['faces'], ch['modes']); K = ck.Checker(G); cg = kk.CGraph(cert['faces'], ch['modes'])
    assert cg.context(Blo, Bhi) == 0
    n = G.n
    cases = [
        ([0, 1, 0], [0.5, -0.25, 1.0]),                       # repeated row: the last value wins
        ([n - 1, 0], [_f32(0.1), _f32(-0.7)]),                 # out of order: sorted by row
        ([2], [0.0]),                                          # an exact zero multiplier
        ([1, 3], [math.inf, 1.0]), ([1, 3], [-math.inf, 1.0]), ([1], [math.nan]),
        (list(range(n)), [_f32((-1) ** v / (v + 1)) for v in range(n)]),
    ]
    for rows, vals in cases:
        mu = dict(zip(rows, [_f32(v) for v in vals]))
        s = _py_s(G, K, (Blo, Bhi), slo, shi, mu)
        rc, closed, cs = cg.hint(slo, shi, (rows, [_f32(v) for v in vals]))
        assert rc == 0
        assert struct.pack('<dd', s.lo, s.hi) == struct.pack('<dd', *cs), (rows, vals, s, cs)
        assert closed == int(s.lo > 0 or s.hi < 0)
    rc, closed, _ = cg.hint(slo, shi, ([], []))
    assert rc == 0 and closed == 0                              # `if not mu: return False`
    rc, closed, _ = cg.hint(slo, shi, ([n], [1.0]))
    assert rc == -21                                            # the reference raises IndexError here
    with pytest.raises(IndexError):
        _py_s(G, K, (Blo, Bhi), slo, shi, {n: 1.0})


def _encode(arr, dtype):
    return base64.b64encode(zlib.compress(np.asarray(arr, dtype=np.dtype(dtype).newbyteorder('<')).tobytes())).decode()


def _outcome(fun, cert):
    try:
        r = fun(cert)
        return ('ok', r['verified'], len(r['failed']), r['open'], len(r['terminal']), r['leaves'])
    except Exception as e:
        return ('raise', type(e).__name__)


def test_errors_match_reference(sample):
    kk = _kk()
    cert, ch, *_ = _first_chart_with_slabs(sample)
    ev = unpack(ch['events'], np.uint8); q = unpack(ch['q'], np.uint16); mu = unpack(ch['mu'], np.float32)
    variants = {
        'events truncated': dict(events=_encode(ev[:-1], np.uint8)),
        'events extended': dict(events=_encode(np.append(ev, 1), np.uint8)),
        'q truncated': dict(q=_encode(q[:-1], np.uint16)),
        'mu extended': dict(mu=_encode(np.append(mu, 1.0), np.float32)),
        'unknown event': dict(events=_encode(np.append(np.array([200], np.uint8), ev), np.uint8)),
    }
    # a hint row out of range: the first hint of the q stream is in the first contraction round; rewrite
    # every hint count c and its rows by replaying the reads is complex, so instead append a linear leaf with
    # a row >= n to a one-leaf chart
    n = len(ch['lo0'])
    G = am.Graph(cert['faces'], ch['modes'])
    leaf = dict(ch, events=_encode([3], np.uint8), q=_encode([1, G.n], np.uint16), mu=_encode([1.0], np.float32),
                n_events=1)
    for name, change in list(variants.items()) + [('hint row out of range', None)]:
        c2 = dict(ch, **change) if change is not None else leaf
        cert2 = dict(cert, charts=[c2])
        a = _outcome(lambda c: kk.check_certificate(c, structural=False), cert2)
        b = _outcome(lambda c: ck.check_certificate(c, structural=False), cert2)
        assert a == b, (name, a, b)


def test_terminal_and_open(sample):
    kk = _kk()
    cert, ch, *_ = _first_chart_with_slabs(sample)
    n = len(ch['lo0'])
    if n < 2: pytest.skip('chart with a single parameter')
    lo0 = [float(x) for x in ch['lo0']]
    # one terminal leaf whose widest coordinate is exactly tol, one just above, one open leaf
    for width_over in (0.0, 1e-3):
        tol = float(cert['tol'])
        hi0 = list(lo0); hi0[0] = lo0[0] + tol * (1 + width_over); hi0[-1] = lo0[-1] + tol / 2
        c2 = dict(ch, lo0=lo0, hi0=hi0, events=_encode([7], np.uint8), q=_encode([], np.uint16),
                  mu=_encode([], np.float32), n_events=1)
        cert2 = dict(cert, theta=[lo0[0], hi0[0]], charts=[c2])
        a = _outcome(lambda c: kk.check_certificate(c, structural=False), cert2)
        b = _outcome(lambda c: ck.check_certificate(c, structural=False), cert2)
        assert a == b, (width_over, a, b)
    c3 = dict(ch, events=_encode([8], np.uint8), q=_encode([], np.uint16), mu=_encode([], np.float32), n_events=1)
    cert3 = dict(cert, charts=[c3])
    assert _outcome(lambda c: kk.check_certificate(c, structural=False), cert3) == \
        _outcome(lambda c: ck.check_certificate(c, structural=False), cert3)


def test_limits(sample):
    kk = _kk()
    cert, ch, Blo, Bhi, *_ = _first_chart_with_slabs(sample)
    cg = kk.CGraph(cert['faces'], ch['modes'])
    assert cg.exclude(Blo, Bhi, 41, 400)[0] == -25 and cg.exclude(Blo, Bhi, 40, 401)[0] == -25
    assert cg.exclude(Blo, Bhi, -1, 400)[0] == -25
    rc, v, _ = cg.exclude(Blo, Bhi, 40, 0)                       # no piece allowed: nothing is tested
    assert rc == 0 and v == 0
    res = kk.check_chart(cg, ch, float(cert['tol']), 41, 400)
    assert res['rc'] == -25


# ------------------------------------------------------------------------------------------- test power

MUTANTS = {
    # rounding direction: the hint's sum steps its lower endpoint one ulp inward
    'hint lower endpoint one ulp inward':
        ('s = kc_add(s, kc_mul(mv, w->Fm[v]));                     /* s = s + mv * Fm[v] */',
         's = kc_add(s, kc_mul(mv, w->Fm[v])); s.lo = nextafter(s.lo, INFINITY);'),
    # rounding direction: Y_j from unrounded differences
    'Y without outward rounding':
        ('w->Y[j] = kc_hull(kc_sub(kc_point(lo[j]), kc_point(w->m[j])), kc_sub(kc_point(hi[j]), kc_point(w->m[j])));',
         'w->Y[j] = kc_hull(kc_point(lo[j] - w->m[j]), kc_point(hi[j] - w->m[j]));'),
    # summation order: hint rows in decreasing order
    'hint rows in decreasing order':
        ('for (int v = 0; v < w->n; v++) {\n        if (w->last[v] < 0) continue;',
         'for (int v = w->n - 1; v >= 0; v--) {\n        if (w->last[v] < 0) continue;'),
    # comparison: only the lower endpoint is tested
    'excludes_zero tests one side':
        ('return x.lo > 0.0 || x.hi < 0.0;', 'return x.lo > 0.0;'),
    # comparison: a slack that accepts intervals touching 0 (unsound)
    'excludes_zero with slack':
        ('return x.lo > 0.0 || x.hi < 0.0;', 'return x.lo > -1e-9 || x.hi < 1e-9;'),
    # the tightening flag inverted
    'tightening flag inverted':
        ('int rc = kc_corners(w->g, lo, hi, NULL, tighten, w->C, &ok);',
         'int rc = kc_corners(w->g, lo, hi, NULL, !tighten, w->C, &ok);'),
    # the validity flag ignored
    'validity flag ignored':
        ("if (!ok) { *rule = 1; return KC_OK; }", "(void)ok;"),
}


def _build_mutant(tmp, name, old, new):
    src = (ROOT / 'csrc' / 'kc_leaf.c').read_text()
    assert src.count(old) == 1, 'mutation site not found exactly once: %s' % name
    d = tmp / ('m%d' % abs(hash(name)))
    d.mkdir()
    for f in ('kc_arith.c', 'kc_arith.h', 'kc_model.c', 'kc_model.h', 'kc_replay.c', 'kc_replay.h', 'kc_leaf.h'):
        shutil.copy(ROOT / 'csrc' / f, d / f)
    (d / 'kc_leaf.c').write_text(src.replace(old, new))
    lib = d / 'libkisscheck.so'
    subprocess.run([shutil.which('cc') or 'gcc', *CFLAGS, '-shared', '-o', str(lib), 'kc_arith.c', 'kc_model.c',
                    'kc_replay.c', 'kc_leaf.c', '-lflint', '-lm'], cwd=d, check=True)
    return lib


@pytest.mark.skipif(shutil.which('cc') is None and shutil.which('gcc') is None, reason='no C compiler')
def test_power(sample, tmp_path):
    D = _diff()
    lines = sample                     # stops at the first detection per mutant
    detected = {}
    try:
        for name, (old, new) in MUTANTS.items():
            lib = _build_mutant(tmp_path, name, old, new)
            D.LIBPATH = str(lib)
            found = []
            for line in lines:
                found += _problems(D.check_one(line))
                if found: break
            detected[name] = bool(found)
    finally:
        D.LIBPATH = None
    assert all(detected.values()), detected
    assert len(detected) >= 5


# ------------------------------------------------------------------- V3: theta outside the model's domain

def test_v3_forged_theta_not_verified_by_kernel_front_end():
    """V3-FIX (audit B-1): the forgeries (a)-(d) of tests/test_certificates.py, a one-leaf certificate over the
    a-priori box at a theta range outside [THETA_FLOOR, THETA_CEIL] or not finite, are not verified by the
    front end (C kernel plus fallback): it either raises as the reference does, or reports verified False with
    the V3 domain problem.  Before the fix, (a) and (b) were verified True."""
    path = CERTS / 'certificates-cond-21-of-40.jsonl'
    if not path.exists(): pytest.skip('no certificates')
    with open(path) as f: d = json.loads(f.readlines()[1074])
    code, c0 = d['code'], d['certificate']
    r = ck.arb(ck.R_MIN)
    h = [ck.down(2 * (r / (1 + r)).asin()), ck.up(ck.arb(2).atan())]
    inf, nan = float('inf'), float('nan')
    pk = lambda a, dt: base64.b64encode(zlib.compress(np.asarray(a, dt).tobytes(), 9)).decode()
    K = _kk()
    for theta in ([-5.0, h[1]], [1.0471975511965976, 7.0], [-inf, h[1]], [nan, h[1]], [h[0], inf], [h[0], nan]):
        lo, hi, _ = ck.apriori_box(am.Graph(c0['faces'], {}), theta)
        lo = list(lo); hi = list(hi); lo[0], hi[0] = theta
        ch = dict(modes={}, lo0=lo, hi0=hi, events=pk([2], np.uint8), q=pk([], np.uint16), mu=pk([], np.float32),
                  n_events=1, complete=True)
        cert = dict(version=2, theta=list(theta), tol=c0['tol'], faces=c0['faces'], complete=True, charts=[ch])
        try:
            out = K.check_certificate(cert, code=code, gen_faces=cert['faces'])
        except ValueError:
            continue
        assert not out['verified'], theta
        assert any(p.startswith('V3: theta range outside the model domain') for p in out['structural']), theta
