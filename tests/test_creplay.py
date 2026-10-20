"""The C replay (csrc/kc_replay.c, checklist F3b) is bit-for-bit the Python reference (kisscrit/core/certificate.py:replay).

Status: untrusted (a test).  Spec v0.3 section 4.2: 'The replay reproduces certificate.replay bit for bit ... A
test compares the box sequences against the Python replay on 1,000 certificates.'

What is compared, node by node in preorder: the node kind; every endpoint of every box (node box; for a
contraction round also the new box and, per removed slab, the slab box and the box it was cut from) as IEEE bit
patterns, so 0.0 and -0.0 differ and the tolerance is 0; the split coordinate; each slab's coordinate k and side;
every multiplier hint (None, or the dict row -> value, with values compared as bit patterns); and how the replay
ends: complete, or the reference's exception (ValueError by message, or IndexError) against the C error code,
after the same nodes.

Inputs:
  * test_replay_matches_reference_on_cond_certificates: KISSCRIT_CREPLAY_CERTS (default 1000) certificates drawn
    without replacement, seed KISSCRIT_CREPLAY_SEED (default 20261020), from all lines of
    results/completeness/certificates/certificates-cond-*-of-40.jsonl (file order sorted, then line order), every
    run chart of each;
  * test_malformed_streams: hand-made charts for every point where the reference raises (each ValueError message
    and each IndexError read), the truncated and extended streams of tests/test_certificates.py
    test_v5_event_stream_truncated_or_extended, and edge cases (no mu stream, empty, short and repeated hints,
    a coordinate contracted twice in one round, point boxes, -0.0, infinities, NaN);
  * test_mutated_real_charts: random edits of events, q, mu and the starting box of small real charts;
  * test_power: copies of kc_replay.c with one deliberate change to the arithmetic or the order (an fma() grid
    value, a different grid formula, a different midpoint formula, upper slab before lower), compiled with the
    same flags, must disagree with the reference.
Skipped when build/libkisscheck.so is not built (run `make -C csrc`); test_power also needs a C compiler.
"""
import base64, ctypes, functools, glob, json, os, pathlib, random, shutil, struct, subprocess, zlib
import numpy as np
import pytest
from kisscrit.core import carith
from kisscrit.core.certificate import replay, unpack, LEAF

pytestmark = pytest.mark.skipif(not carith.available(), reason='C kernel not built (make -C csrc)')

ROOT = pathlib.Path(__file__).resolve().parents[1]
CERTS = ROOT / 'results' / 'completeness' / 'certificates'
NCERTS = int(os.environ.get('KISSCRIT_CREPLAY_CERTS', '1000'))
SEED = int(os.environ.get('KISSCRIT_CREPLAY_SEED', '20261020'))
CFLAGS = ['-O2', '-std=c11', '-Wall', '-Wextra', '-fno-fast-math', '-ffp-contract=off', '-fPIC']   # csrc/Makefile

# C return codes (kc_replay.h) and the reference's way of ending
ENDS = {0: 'done', -1: 'outside', -2: 'empty round', -3: 'ended inside', -4: 'unknown', -5: 'trailing',
        -6: 'index', -7: 'C limit', -8: 'C stack'}
MESSAGES = {'contraction outside the box': 'outside', 'empty contraction round': 'empty round',
            'event stream ended inside the tree': 'ended inside', 'trailing events': 'trailing'}


# ------------------------------------------------------------------------------------------- ctypes binding

class KcHint(ctypes.Structure):
    _fields_ = [('present', ctypes.c_int), ('count', ctypes.c_int64), ('rows', ctypes.c_void_p),
                ('vals', ctypes.c_void_p)]


class KcSlab(ctypes.Structure):
    _fields_ = [('lo', ctypes.c_void_p), ('hi', ctypes.c_void_p), ('cut_lo', ctypes.c_void_p),
                ('cut_hi', ctypes.c_void_p), ('k', ctypes.c_int), ('side', ctypes.c_int), ('hint', KcHint)]


class KcNode(ctypes.Structure):
    _fields_ = [('kind', ctypes.c_int), ('k', ctypes.c_int), ('lo', ctypes.c_void_p), ('hi', ctypes.c_void_p),
                ('new_lo', ctypes.c_void_p), ('new_hi', ctypes.c_void_p), ('nslabs', ctypes.c_int),
                ('slabs', ctypes.POINTER(KcSlab)), ('hint', KcHint)]


def _bind(lib):
    vp, i64 = ctypes.c_void_p, ctypes.c_int64
    lib.kc_replay_stack_need.restype = i64; lib.kc_replay_stack_need.argtypes = [vp, i64]
    lib.kc_replay_work_doubles.restype = i64; lib.kc_replay_work_doubles.argtypes = [ctypes.c_int, i64]
    lib.kc_replay_sizeof.restype = i64; lib.kc_replay_sizeof.argtypes = []
    lib.kc_replay_init.restype = ctypes.c_int
    lib.kc_replay_init.argtypes = [vp, ctypes.c_int, vp, vp, vp, i64, vp, i64, vp, i64, ctypes.c_int, vp, i64]
    lib.kc_replay_next.restype = ctypes.c_int; lib.kc_replay_next.argtypes = [vp, ctypes.POINTER(KcNode)]
    return lib


def _lib():
    lib = carith.load()
    if not getattr(lib, '_replay_ready', False):
        _bind(lib); lib._replay_ready = True
    return lib


# ------------------------------------------------------------------------------------------- canonical nodes

def _bits(a):
    return np.asarray(a, dtype='<f8').tobytes()


def _hint(d):
    """A hint as a comparable value: None, or the dict's items in order with values as double bit patterns."""
    if d is None: return None
    return tuple((int(r), struct.pack('<d', float(v))) for r, v in d.items())


def py_replay(chart):
    """The reference's nodes as comparable tuples, and how it ended."""
    out = []
    try:
        for kind, lo, hi, extra in replay(chart):
            if kind == 'split':
                out.append(('split', _bits(lo), _bits(hi), int(extra)))
            elif kind == 'contract':
                removed, (nlo, nhi) = extra
                slabs = tuple((_bits(a), _bits(b), _bits(c), _bits(d), int(k), int(side), _hint(mu))
                              for a, b, c, d, k, side, mu in removed)
                out.append(('contract', _bits(lo), _bits(hi), (slabs, _bits(nlo), _bits(nhi))))
            else:
                out.append((kind, _bits(lo), _bits(hi), _hint(extra)))
    except ValueError as e:
        msg = str(e)
        return out, MESSAGES.get(msg, 'unknown' if msg.startswith('unknown event') else 'ValueError: ' + msg)
    except IndexError:
        return out, 'index'
    return out, 'done'


def _streams(chart):
    ev = np.ascontiguousarray(unpack(chart['events'], np.uint8), dtype=np.uint8)
    q = np.ascontiguousarray(unpack(chart['q'], np.uint16), dtype=np.uint16)
    has_mu = 'mu' in chart
    mu = np.ascontiguousarray(unpack(chart['mu'], np.float32) if has_mu else np.zeros(0), dtype=np.float32)
    lo0 = np.ascontiguousarray(np.array(chart['lo0'], float)); hi0 = np.ascontiguousarray(np.array(chart['hi0'], float))
    return ev, q, mu, has_mu, lo0, hi0


def c_replay(chart, lib=None):
    """The C replay's nodes as the same comparable tuples, and how it ended."""
    lib = lib or _lib()
    ev, q, mu, has_mu, lo0, hi0 = _streams(chart)
    n = len(lo0); assert len(hi0) == n
    cap = lib.kc_replay_stack_need(ev.ctypes.data, len(ev))
    work = np.zeros(max(1, lib.kc_replay_work_doubles(n, cap)))
    state = np.zeros(lib.kc_replay_sizeof() // 8 + 1)              # 8-byte aligned storage for kc_replay
    sp = state.ctypes.data
    rc = lib.kc_replay_init(sp, n, lo0.ctypes.data, hi0.ctypes.data, ev.ctypes.data, len(ev), q.ctypes.data,
                            len(q), mu.ctypes.data, len(mu), int(has_mu), work.ctypes.data, cap)
    assert rc == 0, rc
    nb = 8 * n

    def b(p): return ctypes.string_at(p, nb) if nb else b''

    def h(x):
        if not x.present: return None
        if x.count == 0: return ()
        rows = np.frombuffer(ctypes.string_at(x.rows, 2 * x.count), dtype=np.uint16)
        vals = np.frombuffer(ctypes.string_at(x.vals, 4 * x.count), dtype=np.float32)
        return _hint(dict(zip([int(r) for r in rows], [float(v) for v in vals])))

    out = []; node = KcNode()
    while True:
        rc = lib.kc_replay_next(sp, ctypes.byref(node))
        if rc != 1: break
        if node.kind == 16:
            out.append(('split', b(node.lo), b(node.hi), node.k))
        elif node.kind == 15:
            slabs = tuple((b(s.lo), b(s.hi), b(s.cut_lo), b(s.cut_hi), s.k, s.side, h(s.hint))
                          for s in node.slabs[:node.nslabs])
            out.append(('contract', b(node.lo), b(node.hi), (slabs, b(node.new_lo), b(node.new_hi))))
        else:
            out.append((LEAF[node.kind], b(node.lo), b(node.hi), h(node.hint)))
    assert lib.kc_replay_next(sp, ctypes.byref(node)) == rc               # the end is sticky
    return out, ENDS.get(rc, 'code %d' % rc)


def compare(chart, lib=None):
    """(number of differing positions in the node sequences, reference result, C result)."""
    p_nodes, p_end = py_replay(chart); c_nodes, c_end = c_replay(chart, lib)
    diff = sum(1 for a, b in zip(p_nodes, c_nodes) if a != b) + abs(len(p_nodes) - len(c_nodes))
    return diff, (len(p_nodes), p_end), (len(c_nodes), c_end)


# ------------------------------------------------------------------------------------------- inputs

def _pack(arr):
    return base64.b64encode(zlib.compress(np.ascontiguousarray(arr).astype(arr.dtype.newbyteorder('<')).tobytes())).decode()


def chart(lo0, hi0, ev, q=(), mu=None):
    """A chart dict as certificates store it; mu=None leaves out the mu stream."""
    c = {'lo0': [float(x) for x in lo0], 'hi0': [float(x) for x in hi0],
         'events': _pack(np.array(ev, np.uint8)), 'q': _pack(np.array(q, np.uint16))}
    if mu is not None: c['mu'] = _pack(np.array(mu, np.float32))
    return c


@functools.lru_cache(maxsize=None)
def _all_lines():
    """(file, line index) of every cond certificate, files sorted, lines in order."""
    lines = []
    for f in sorted(glob.glob(str(CERTS / 'certificates-cond-*-of-40.jsonl'))):
        with open(f, 'rb') as fh:
            lines += [(f, i) for i, _ in enumerate(fh)]
    return tuple(lines)


def _sample(ncerts, seed):
    """ncerts certificate lines drawn without replacement (random.Random(seed).sample) from all cond lines."""
    files = sorted(glob.glob(str(CERTS / 'certificates-cond-*-of-40.jsonl')))
    lines = list(_all_lines())
    if not lines: pytest.skip('no cond certificates in results/')
    pick = sorted(random.Random(seed).sample(lines, min(ncerts, len(lines))))
    want = {}
    for f, i in pick: want.setdefault(f, set()).add(i)
    out = []
    for f in files:
        if f not in want: continue
        with open(f) as fh:
            for i, line in enumerate(fh):
                if i in want[f]: out.append(json.loads(line))
    return out


def _run_charts(certs):
    for r in certs:
        for ch in r['certificate']['charts']:
            if not ch.get('not_run'): yield ch


# ------------------------------------------------------------------------------------------- tests

def test_replay_matches_reference_on_cond_certificates():
    certs = _sample(NCERTS, SEED)
    ncharts = nnodes = 0; bad = []; ends = {}
    for ch in _run_charts(certs):
        d, p, c = compare(ch)
        ncharts += 1; nnodes += p[0]; ends[p[1]] = ends.get(p[1], 0) + 1
        if d or p != c: bad.append((d, p, c))
    print('\nC replay vs certificate.replay: %d certificates (seed %d), %d charts, %d nodes, ends %s, mismatching charts %d'
          % (len(certs), SEED, ncharts, nnodes, ends, len(bad)))
    assert len(certs) == min(NCERTS, len(certs)) and ncharts > 0
    assert not bad, bad[:5]


B2 = ([0.0, 0.0], [1.0, 1.0])

# (name, chart, expected end of the reference): each ValueError message and IndexError read, then edge cases
CASES = [
    ('outside: qlo above qhi', chart(*B2, [15, 64, 2], [40000, 35536, 0, 0]), 'outside'),
    ('outside: NaN box', chart([float('nan')], [1.0], [15, 64, 2], [1, 1, 0, 0]), 'outside'),
    ('outside: second coordinate', chart(*B2, [15, 64, 65, 2], [0, 0, 2, 65535, 0, 0], []), 'outside'),
    ('empty round before a leaf', chart(*B2, [15, 2]), 'empty round'),
    ('empty round at the end', chart(*B2, [15]), 'empty round'),
    ('second round empty', chart(*B2, [15, 64, 15, 2], [0, 1, 0], []), 'empty round'),
    ('ended inside: no events', chart(*B2, []), 'ended inside'),
    ('ended inside: split without children', chart(*B2, [16]), 'ended inside'),
    ('ended inside: one child', chart(*B2, [16, 2]), 'ended inside'),
    ('ended inside: after a round', chart(*B2, [15, 64], [0, 0]), 'ended inside'),
    ('unknown 0', chart(*B2, [0]), 'unknown'),
    ('unknown 9', chart(*B2, [9]), 'unknown'),
    ('unknown 14', chart(*B2, [14]), 'unknown'),
    ('unknown: contract without round', chart(*B2, [64], [0, 0]), 'unknown'),
    ('unknown 255', chart(*B2, [255]), 'unknown'),
    ('unknown in the second child', chart(*B2, [16, 2, 12]), 'unknown'),
    ('trailing event', chart(*B2, [2, 2]), 'trailing'),
    ('trailing q', chart(*B2, [2], [7]), 'trailing'),
    ('trailing mu', chart(*B2, [2], [], [1.0]), 'trailing'),
    ('short mu in a hint', chart(*B2, [3], [2, 0, 1], [0.5]), 'trailing'),
    ('hint without mu stream', chart(*B2, [15, 64, 2], [1, 1, 2, 7, 8, 0]), 'trailing'),
    ('index: split coordinate >= n', chart(*B2, [18]), 'index'),
    ('index: contract coordinate >= n', chart(*B2, [15, 69, 2], [0, 0]), 'index'),
    ('index: contract coordinate >= n after hints', chart(*B2, [15, 70, 2], [5, 5, 1, 3, 0], [1.0]), 'index'),
    ('index: q ends inside a pair', chart(*B2, [15, 64, 2], [0]), 'index'),
    ('index: q empty at a pair', chart(*B2, [15, 64, 2], []), 'index'),
    ('index: q ends before a hint count', chart(*B2, [15, 64, 2], [1, 0]), 'index'),
    ('index: q ends before the upper hint count', chart(*B2, [15, 64, 2], [1, 1, 0]), 'index'),
    ('index: linear leaf without q', chart(*B2, [3], [], []), 'index'),
    ('index: hint rows past the end, then a pair', chart(*B2, [15, 64, 65, 2], [1, 0, 9, 1], []), 'index'),
    # complete trees with edge cases
    ('linear leaf without mu stream', chart(*B2, [3]), 'done'),
    ('linear leaf with hint', chart(*B2, [3], [2, 5, 1], [0.5, -1.0]), 'done'),
    ('linear leaf, empty hint', chart(*B2, [3], [0], []), 'done'),
    ('repeated rows in a hint', chart(*B2, [3], [3, 4, 4, 2], [1.0, 2.0, 3.0]), 'done'),
    ('both slabs, two coordinates', chart(*B2, [15, 64, 65, 2], [100, 200, 1, 0, 1, 1, 7, 9, 2, 3, 4, 1, 5],
                                          [0.25, -1.0, 1.0, 0.5, -0.5]), 'done'),
    ('coordinate contracted twice in a round', chart(*B2, [15, 64, 64, 2], [100, 100, 1, 0, 1, 1, 30000, 0, 1, 2],
                                                    [1.0, 2.0, 3.0]), 'done'),
    ('lower hint kept, upper bound replaced', chart(*B2, [15, 64, 64, 2], [100, 0, 1, 3, 0, 100, 1, 4],
                                                   [1.0, 2.0]), 'done'),
    ('two rounds then a split', chart(*B2, [15, 65, 15, 64, 16, 1, 2], [0, 1000, 0, 7, 0, 0], []), 'done'),
    ('point box, q inside', chart([1.0, 2.0], [1.0, 2.0], [15, 64, 2], [5, 5, 0, 0]), 'done'),
    ('grid ends q = 0 and QMAX', chart(*B2, [15, 64, 2], [0, 0]), 'done'),
    ('negative zero', chart([-0.0, -1.0], [0.0, -0.0], [16, 17, 15, 64, 2, 2, 2], [1, 1, 0, 0]), 'done'),
    ('infinite box split', chart([-float('inf'), 0.0], [float('inf'), 1.0], [16, 2, 2]), 'done'),
    ('outside: hi - lo overflows', chart([-1.7e308, 0.0], [1.7e308, 1.0], [15, 64, 16, 2, 2], [12345, 23456, 0, 0]),
     'outside'),
    ('subnormal box', chart([5e-324, 0.0], [2e-323, 1.0], [15, 64, 16, 2, 2], [12345, 23456, 0, 0]), 'done'),
    ('wide box split (hi - lo inexact)', chart([-4.898619485211566, 0.0], [76091.34587177234, 1.0], [16, 2, 2]), 'done'),
    ('n = 0', chart([], [], [2]), 'done'),
    ('deep tree', chart(*B2, [16] * 40 + [2] * 41), 'done'),
]


@pytest.mark.parametrize('name,ch,expected', CASES, ids=[c[0] for c in CASES])
def test_malformed_streams(name, ch, expected):
    d, p, c = compare(ch)
    assert p[1] == expected, (name, p)                 # the case reaches the intended point of the reference
    assert d == 0 and p == c, (name, d, p, c)


def test_every_raise_site_is_covered():
    """Each way the reference can end (complete, its four ValueError messages plus 'unknown event', IndexError)
    occurs among CASES; together with the comments in kc_replay.c this covers every raise in certificate.py."""
    import inspect
    from kisscrit.core import certificate
    src = inspect.getsource(certificate.replay)
    raised = {m for m in MESSAGES if m in src} | ({'unknown'} if "'unknown event" in src else set())
    assert len(raised) == 5 and src.count('raise ValueError') == 5
    seen = {c[2] for c in CASES}
    assert {MESSAGES[m] for m in MESSAGES} | {'unknown', 'index', 'done'} <= seen


def _small_real_charts(k=12, max_events=400):
    out = []
    for f in sorted(glob.glob(str(CERTS / 'certificates-cond-*-of-40.jsonl')))[:3]:
        with open(f) as fh:
            for line in fh:
                for ch in json.loads(line)['certificate']['charts']:
                    if not ch.get('not_run') and 30 <= ch.get('n_events', 0) <= max_events and 'mu' in ch \
                            and len(unpack(ch['mu'], np.float32)) > 0:
                        out.append(ch)
                if len(out) >= 4 * k: break
    if not out: pytest.skip('no cond certificates in results/')
    return random.Random(SEED).sample(out, min(k, len(out)))


def test_truncated_and_extended_real_streams():
    """As tests/test_certificates.py test_v5_event_stream_truncated_or_extended, for events, q and mu."""
    ends = {}
    for ch in _small_real_charts():
        ev, q, mu, _, _, _ = _streams(ch)
        variants = [('events', ev[:-1]), ('events', np.concatenate([ev, np.array([2], np.uint8)])),
                    ('q', q[:-1]), ('q', np.concatenate([q, np.array([0], np.uint16)])),
                    ('mu', mu[:-1]), ('mu', np.concatenate([mu, np.array([1.0], np.float32)]))]
        for key, arr in variants:
            m = dict(ch); m[key] = _pack(arr)
            d, p, c = compare(m)
            ends[p[1]] = ends.get(p[1], 0) + 1
            assert p[1] != 'done' and d == 0 and p == c, (key, d, p, c)
    print('\ntruncated/extended real streams, reference ends:', ends)


def test_mutated_real_charts():
    rng = random.Random(SEED + 1); ends = {}; total = 0
    for ch in _small_real_charts():
        ev, q, mu, _, lo0, hi0 = _streams(ch)
        for _ in range(60):
            e, qq, mm, lo, hi = ev.copy(), q.copy(), mu.copy(), lo0.copy(), hi0.copy()
            what = rng.randrange(8)
            if what == 0: e[rng.randrange(len(e))] = rng.randrange(256)
            elif what == 1: e = np.delete(e, rng.randrange(len(e)))
            elif what == 2: e = np.insert(e, rng.randrange(len(e) + 1), rng.choice([2, 3, 15, 16, 17, 64, 65, 7]))
            elif what == 3 and len(qq): qq[rng.randrange(len(qq))] = rng.randrange(65536)
            elif what == 4 and len(qq): qq = np.delete(qq, rng.randrange(len(qq)))
            elif what == 5: qq = np.insert(qq, rng.randrange(len(qq) + 1), rng.randrange(4))
            elif what == 6 and len(mm): mm[rng.randrange(len(mm))] = np.float32(rng.uniform(-1, 1))
            else:
                j = rng.randrange(len(lo)); lo[j] = np.nextafter(lo[j], -np.inf) if rng.random() < 0.5 else -lo[j]
            m = dict(ch, events=_pack(e.astype(np.uint8)), q=_pack(qq.astype(np.uint16)), mu=_pack(mm.astype(np.float32)),
                     lo0=[float(x) for x in lo], hi0=[float(x) for x in hi])
            d, p, c = compare(m)
            total += 1; ends[p[1]] = ends.get(p[1], 0) + 1
            assert d == 0 and p == c, (what, d, p, c)
    print('\nmutated real charts: %d, reference ends %s' % (total, ends))


# ------------------------------------------------------------------------------------------- test power

MUTANTS = {
    'fma grid value': [('#include <string.h>\n', '#include <string.h>\n#include <math.h>\n'),
                       ('    double p = d * t;\n    double g = lo0 + p;\n', '    double g = fma(d, t, lo0);\n')],
    'other grid formula': [('    double p = d * t;\n    double g = lo0 + p;\n',
                            '    double g = (1.0 - t) * lo0 + t * hi0; (void)d;\n')],
    'other midpoint formula': [('    double s = lo + hi;\n    double m = 0.5 * s;\n',
                                '    double m = lo + 0.5 * (hi - lo);\n')],
    'midpoint one ulp low': [('#include <string.h>\n', '#include <string.h>\n#include <math.h>\n'),
                             ('    double m = 0.5 * s;\n', '    double m = nextafter(0.5 * s, -INFINITY);\n')],
    'upper slab first': [('for (int side = -1; side <= 1; side += 2)', 'for (int side = 1; side >= -1; side -= 2)')],
}


@pytest.mark.parametrize('mutant', list(MUTANTS))
def test_power(mutant, tmp_path):
    """The comparison detects a deliberate change: each mutant copy of kc_replay.c disagrees with the reference
    on the hand-made cases or on 60 sampled certificates.  (The formula lo + 0.5 (hi - lo) agrees with
    0.5 (lo + hi) whenever hi - lo is exact, e.g. by Sterbenz when hi/2 <= lo <= 2 hi; the case 'wide box split'
    is one where it is not.)"""
    cc = shutil.which(os.environ.get('CC', 'cc'))
    if not cc: pytest.skip('no C compiler')
    src = (ROOT / 'csrc' / 'kc_replay.c').read_text()
    for old, new in MUTANTS[mutant]:
        assert src.count(old) == 1, 'mutation target not found: %r' % old      # the test stays honest if the code changes
        src = src.replace(old, new)
    (tmp_path / 'kc_replay_mutant.c').write_text(src)
    so = tmp_path / 'mutant.so'
    subprocess.run([cc] + CFLAGS + ['-I', str(ROOT / 'csrc'), '-shared', '-o', str(so),
                                    str(tmp_path / 'kc_replay_mutant.c'), '-lm'], check=True)
    lib = _bind(ctypes.CDLL(str(so)))
    nodes = charts = 0
    for ch in [c[1] for c in CASES] + list(_run_charts(_sample(60, SEED))):
        d, p, c = compare(ch, lib)
        nodes += d; charts += bool(d or p != c)
    print('\nmutant %r: %d mismatching nodes in %d charts' % (mutant, nodes, charts))
    assert charts > 0
