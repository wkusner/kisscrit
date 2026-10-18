"""The C model (csrc/kc_model.c, checklist F3a) is bit-for-bit the Python reference (kisscrit/core/model.py).

Status: untrusted (a test).  Spec v0.3 section 7 step 3: 'the C model (corners, vertex sums, Jacobians,
tightening) ... diff them against Python value by value on random boxes'.

Inputs: candidates and charts from the `cond` certificates (faces, pentagon modes, a-priori boxes lo0/hi0).
For each chart the test evaluates, in both implementations:
  * the a-priori box, with and without monotonicity tightening;
  * random sub-boxes, each coordinate a random sub-interval of log-uniform width from 1e-12 up to the full width
    (so both the tightening branch, width > 0.01, and the narrow branch occur);
  * point boxes m = [m, m] with the clamp bound from an enclosing box (blo != lo), as the checker's linear test;
  * the first boxes of the chart's own replay (certificate.replay), which are the boxes the checker meets.
It compares every endpoint of every corner value and gradient, the validity flag, the vertex sums and their
gradients (the Jacobian rows), the edge directions, and errors (the reference's ArithmeticError against the C
error code).  "Bit-exact" means equal IEEE bit patterns (so 0.0 and -0.0 differ); the tolerance is 0.

Size: KISSCRIT_CMODEL_BOXES (default 1200) boxes in total; KISSCRIT_CMODEL_SEED (default 20261017).
Skipped when build/libkisscheck.so is not built (run `make -C csrc`).
"""
import collections, ctypes, json, math, os, pathlib, random, struct
import pytest
from kisscrit.core import model as M
from kisscrit.core import carith
from kisscrit.core.certificate import replay

pytestmark = pytest.mark.skipif(not carith.available(), reason='C kernel not built (make -C csrc)')

ROOT = pathlib.Path(__file__).resolve().parents[1]
CERTS = ROOT / 'results' / 'completeness' / 'certificates'
BOXES = int(os.environ.get('KISSCRIT_CMODEL_BOXES', '1200'))
SEED = int(os.environ.get('KISSCRIT_CMODEL_SEED', '20261017'))
MAXDEG = 32
ERR = {0: None, 1: 'orient', 2: 'rot', 3: 'sphere', 4: 'limit', 5: 'disjoint', 6: 'tighten', 7: 'mode'}


class KcAd(ctypes.Structure):
    _fields_ = [('v', carith.KcIv), ('g', carith.KcIv * 3)]


def _lib():
    lib = carith.load()
    if getattr(lib, '_model_ready', False): return lib
    P = ctypes.POINTER
    lib.kc_graph_new.restype = ctypes.c_void_p
    lib.kc_graph_new.argtypes = [ctypes.c_int, P(ctypes.c_int), P(ctypes.c_int), ctypes.c_char_p, P(ctypes.c_int)]
    lib.kc_graph_free.argtypes = [ctypes.c_void_p]
    for f in ('kc_graph_n', 'kc_graph_npar', 'kc_graph_nedges'):
        getattr(lib, f).restype = ctypes.c_int; getattr(lib, f).argtypes = [ctypes.c_void_p]
    lib.kc_graph_rot.restype = ctypes.c_int
    lib.kc_graph_rot.argtypes = [ctypes.c_void_p, ctypes.c_int, P(ctypes.c_int), P(ctypes.c_int), P(ctypes.c_int)]
    lib.kc_graph_edges.argtypes = [ctypes.c_void_p, P(ctypes.c_int)]
    lib.kc_corners.restype = ctypes.c_int
    lib.kc_corners.argtypes = [ctypes.c_void_p, P(ctypes.c_double), P(ctypes.c_double), P(ctypes.c_double),
                               ctypes.c_int, P(KcAd), P(ctypes.c_int)]
    lib.kc_vertex_sums.argtypes = [ctypes.c_void_p, P(KcAd), P(carith.KcIv), P(carith.KcIv)]
    lib.kc_directions.argtypes = [ctypes.c_void_p, P(KcAd), P(carith.KcIv)]
    lib._model_ready = True
    return lib


class CGraph:
    """The C graph for faces and modes; .err is the C error code (0 on success)."""
    def __init__(self, faces, modes):
        lib = _lib(); self.lib = lib
        self.nf = len(faces)
        sizes = (ctypes.c_int * self.nf)(*[len(f) for f in faces])
        flat = [v for f in faces for v in f]
        verts = (ctypes.c_int * len(flat))(*flat)
        md = bytes((ord(modes.get(str(i), modes.get(i, '\0'))) for i in range(self.nf)))
        err = ctypes.c_int(0)
        self.h = lib.kc_graph_new(self.nf, sizes, verts, md, ctypes.byref(err))
        self.err = err.value
        if self.h:
            self.n = lib.kc_graph_n(self.h); self.npar = lib.kc_graph_npar(self.h)

    def __del__(self):
        if getattr(self, 'h', None): self.lib.kc_graph_free(self.h)

    def rot(self):
        out = []
        for v in range(self.n):
            w = (ctypes.c_int * MAXDEG)(); fi = (ctypes.c_int * MAXDEG)(); pos = (ctypes.c_int * MAXDEG)()
            d = self.lib.kc_graph_rot(self.h, v, w, fi, pos)
            out.append([(w[i], fi[i], pos[i]) for i in range(d)])
        return out

    def edges(self):
        ne = self.lib.kc_graph_nedges(self.h); a = (ctypes.c_int * (2 * ne))()
        self.lib.kc_graph_edges(self.h, a)
        return [(a[2 * e], a[2 * e + 1]) for e in range(ne)]

    def evaluate(self, lo, hi, blo, tighten, sizes, rot):
        """(error name or None, corners as model.py's C, ok, vertex sums, directions)."""
        D = ctypes.c_double * self.npar
        C = (KcAd * (5 * self.nf))(); ok = ctypes.c_int(-1)
        rc = self.lib.kc_corners(self.h, D(*lo), D(*hi), D(*blo) if blo is not None else None, int(tighten),
                                 C, ctypes.byref(ok))
        if rc: return ERR[rc], None, None, None, None
        iv = lambda x: (x.lo, x.hi)
        corners = [[(iv(C[5 * f + q].v),) + tuple(iv(C[5 * f + q].g[i]) for i in range(3)) for q in range(sizes[f])]
                   for f in range(self.nf)]
        S = (carith.KcIv * self.n)(); G = (carith.KcIv * (self.n * self.npar))()
        self.lib.kc_vertex_sums(self.h, C, S, G)
        vs = [(iv(S[v]), [iv(G[v * self.npar + j]) for j in range(self.npar)]) for v in range(self.n)]
        Dd = (carith.KcIv * (self.n * MAXDEG))()
        self.lib.kc_directions(self.h, C, Dd)
        dirs = {(v, w): iv(Dd[v * MAXDEG + i]) for v in range(self.n) for i, (w, _, _) in enumerate(rot[v])}
        return None, corners, ok.value, vs, dirs


def py_evaluate(G, lo, hi, blo, tighten):
    """The reference on the same input, in the same shape as CGraph.evaluate."""
    iv = lambda x: (x.lo, x.hi)
    try:
        C, ok = G.corners(list(lo), list(hi), None if blo is None else list(blo), None, tighten=tighten)
    except ArithmeticError as e:
        return ('tighten' if 'tightening' in str(e) else 'disjoint'), None, None, None, None
    corners = [[tuple(iv(x) for x in c) for c in row] for row in C]
    vs = [(iv(s), [iv(x) for x in g]) for s, g in G.vertex_sums(C)]
    dirs = {k: iv(x) for k, x in G.directions(C).items()}
    return None, corners, int(ok), vs, dirs


def _flat(obj):
    """All doubles of a nested structure, in a fixed order."""
    if isinstance(obj, float): return [obj]
    if isinstance(obj, dict): return [x for k in sorted(obj) for x in _flat(obj[k])]
    return [x for o in obj for x in _flat(o)]


def _bits(x): return struct.pack('<d', x)


def _charts(limit):
    """(faces, modes, lo0, hi0, chart) from the first certificates of several cond shards, deterministic."""
    out = []
    files = sorted(CERTS.glob('certificates-cond-*-of-40.jsonl'))[:8]
    per_file = -(-limit // len(files))
    for p in files:
        got = 0
        with open(p) as fh:
            for line in fh:
                cert = json.loads(line)['certificate']
                for ch in cert['charts']:
                    out.append((cert['faces'], ch['modes'], ch['lo0'], ch['hi0'], ch)); got += 1
                if got >= per_file: break
    return out


class Tally:
    def __init__(self):
        self.boxes = 0; self.values = 0; self.mismatch = 0; self.maxdiff = 0.0; self.examples = []
        self.ok = {0: 0, 1: 0}; self.errors = 0

    def compare(self, tag, c, p):
        self.boxes += 1
        if p[0] is not None: self.errors += 1
        elif p[2] in self.ok: self.ok[p[2]] += 1
        if c[0] != p[0] or (c[0] is None and c[2] != p[2]):
            self.mismatch += 1; self.examples.append((tag, 'error/ok', c[0], p[0], c[2], p[2])); return
        if c[0] is not None: return
        a = _flat([c[1], c[3], c[4]]); b = _flat([p[1], p[3], p[4]])
        assert len(a) == len(b), tag
        self.values += len(a)
        for x, y in zip(a, b):
            if _bits(x) != _bits(y):
                self.mismatch += 1
                if math.isfinite(x) and math.isfinite(y): self.maxdiff = max(self.maxdiff, abs(x - y))
                else: self.maxdiff = math.inf
                if len(self.examples) < 5: self.examples.append((tag, x, y))


def _subbox(rng, lo0, hi0):
    lo, hi = [], []
    for a, b in zip(lo0, hi0):
        w = (b - a) * 10 ** rng.uniform(-12, 0)
        s = a + rng.random() * (b - a - w)
        lo.append(s); hi.append(min(s + w, b))
    return lo, hi


def test_lstraight_is_the_reference_double():
    assert ctypes.c_double.in_dll(_lib(), 'KC_LSTRAIGHT').value == M.LSTRAIGHT


def test_graph_structure_matches():
    for faces, modes, lo0, hi0, _ in _charts(40):
        G = M.Graph(faces, modes); CG = CGraph(faces, modes)
        assert CG.err == 0
        assert (CG.n, CG.npar) == (G.n, G.npar)
        assert CG.edges() == G.edges
        assert CG.rot() == [list(r) for r in G.rot]


def test_graph_errors_match():
    tet = [[0, 1, 2], [0, 2, 3], [0, 3, 1], [1, 3, 2]]
    assert CGraph(tet, {}).err == 0 and M.Graph(tet, {}).n == 4
    bad = [[0, 1, 2], [0, 1, 2], [0, 3, 1], [1, 3, 2]]          # an oriented edge in two faces
    with pytest.raises(ValueError, match='orientation'): M.Graph(bad, {})
    assert CGraph(bad, {}).err == 1
    two = [[0, 1, 2], [0, 2, 1]]                                 # a sphere with 3 vertices: V - E + F = 2
    assert CGraph(two, {}).err == 0 and M.Graph(two, {}).n == 3
    torus_like = [[0, 1, 2], [0, 2, 3], [0, 3, 1]]              # open disc: Euler characteristic 1
    with pytest.raises((ValueError, KeyError)): M.Graph(torus_like, {})
    assert CGraph(torus_like, {}).err in (2, 3)


def test_missing_pentagon_mode_is_an_error_in_both():
    faces, modes, lo0, hi0, _ = next(c for c in _charts(40) if any(len(f) == 5 for f in c[0]))
    G = M.Graph(faces, {}); CG = CGraph(faces, {})
    with pytest.raises(KeyError): G.corners(lo0, hi0)
    assert CG.evaluate(lo0, hi0, None, True, [len(f) for f in faces], [list(r) for r in G.rot])[0] == 'mode'


def test_corners_vertex_sums_directions_bit_exact():
    rng = random.Random(SEED)
    charts = _charts(max(10, BOXES // 30))
    t = Tally(); per = max(1, BOXES // len(charts)); narrowed = 0; kinds = collections.Counter()
    for ci, (faces, modes, lo0, hi0, ch) in enumerate(charts):
        G = M.Graph(faces, modes); CG = CGraph(faces, modes); assert CG.err == 0
        sizes = G.sizes; rot = [list(r) for r in G.rot]
        cases = [(lo0, hi0, None, True), (lo0, hi0, None, False)]
        for k, (kind, lo, hi, _x) in enumerate(replay(ch)):
            if k >= 4: break
            cases.append((list(lo), list(hi), None, True))
        while len(cases) < per:
            lo, hi = _subbox(rng, lo0, hi0)
            r = rng.random()
            if r < 0.2:
                m = [0.5 * (a + b) for a, b in zip(lo, hi)]
                cases.append((m, m, lo, rng.random() < 0.5))
            else:
                cases.append((lo, hi, None, r < 0.8))
        for j, (lo, hi, blo, tg) in enumerate(cases):
            ref = py_evaluate(G, lo, hi, blo, tg)
            t.compare((ci, j), CG.evaluate(lo, hi, blo, tg, sizes, rot), ref)
            if tg and ref[0] is None and _flat(ref[1]) != _flat(py_evaluate(G, lo, hi, blo, False)[1]):
                narrowed += 1
        kinds.update(len(f) for f in faces)
    print('\ncmodel diff: seed %d, %d charts, %d boxes (ok=1: %d, ok=0: %d, errors: %d, tightening narrowed: %d), '
          '%d doubles compared, %d mismatches, max abs diff %r'
          % (SEED, len(charts), t.boxes, t.ok[1], t.ok[0], t.errors, narrowed, t.values, t.mismatch, t.maxdiff))
    assert t.boxes >= min(BOXES, 1000)
    # not vacuous: both validity outcomes, the tightening branch, and quadrilaterals and pentagons all occur
    assert t.ok[0] > 0 and t.ok[1] > 0 and narrowed > 0, (t.ok, narrowed)
    assert kinds[4] > 0 and kinds[5] > 0, kinds
    assert t.mismatch == 0, t.examples
