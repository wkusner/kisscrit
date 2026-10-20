"""Front end of the C checker kernel (csrc/kc_leaf.c, checklist F4): a ctypes wrapper that runs the kernel
on every chart of a certificate and hands every box the kernel cannot close to the Python reference checker.

Status: untrusted, for now.  The trusted parts are the C kernel (csrc/kc_*.c) and the Python reference
(kisscrit/core/checker.py, unchanged).  This file is not only plumbing: it runs the structural checks and the
starting-box checks (by calling the reference), aggregates `verified`, matches the kernel's units to the
reference's units by position, and decides which units go to the fallback.  A bug here could accept a
certificate.  It is not used by check_certificates.py; before a later checklist item wires it in, it must enter
the trusted base (the TRUSTED list of tests/test_code_standards.py) or its aggregation must be re-derived by
trusted code.  Spec v0.3 section 4.1: the front end decodes the streams, the kernel returns per chart a verdict
per unit and counts per rule, and the Python reference handles the fallback boxes; a box that neither can
close is reported as failed, never accepted.
Implements: nothing mathematical.  The claim of `check_certificate` here is that it returns what
kisscrit/core/checker.py:check_certificate returns (verified, structural, failed, terminal, open, not_run,
leaves), with the work divided differently; `counts` are the exception (the kernel's own attribution for the
units it closes, the reference's for the fallback; informational only).  scripts/completeness/diff_kernel.py
and tests/test_ckernel.py compare the two (computed, not proved).

Units.  A unit is one leaf or one contraction slab of a chart's proof tree, in the order of
kisscrit/core/certificate.py:replay: for each node in preorder, a contraction's slabs in order, or the leaf.
The kernel's verdict codes are kc_leaf.h KC_V_*; their names are in VERDICT.

Fallback.  For a unit the kernel leaves unclosed, the reference's own calls are made, in the reference's order
(check_certificate's loop body): a linear leaf with a hint gets Checker.test_linear_mu on its own box and then
Checker.exclude; another leaf gets Checker.exclude; a slab gets Checker.exclude_slab with the context of its
contraction box.  If the kernel stops on an error (a malformed stream, a model error, a hint row out of
range), the whole chart is handed to checker.check_certificate unchanged, which raises or fails exactly as it
would alone.
"""
import ctypes, time
import numpy as np
from .core import carith, checker as ck, model as am
from .core.certificate import replay, unpack

VERDICT = {0: 'unclosed', 1: 'invalid', 2: 'vsum', 3: 'linear-hint', 4: 'invalid-tight', 5: 'vsum-tight',
           6: 'bisect', 7: 'terminal', 8: 'terminal-too-wide', 9: 'open'}
CLOSED = (1, 2, 3, 4, 5, 6)
COUNTS = ('invalid', 'vsum', 'linear-hint', 'invalid-tight', 'vsum-tight', 'bisect')   # kc_leaf.h KC_C_*
NCOUNT = len(COUNTS)
MAX_DEPTH, MAX_PIECES = 40, 400            # kc_leaf.h KC_MAX_DEPTH, KC_MAX_PIECES (= checker.py's)
MAXDEG = 32                                # kc_model.h KC_MAXDEG
ERRORS = {-1: 'outside', -2: 'empty round', -3: 'ended inside', -4: 'unknown event', -5: 'trailing',
          -6: 'index', -7: 'C limit', -8: 'C stack', -20: 'model', -21: 'hint row', -22: 'capacity',
          -23: 'size', -24: 'memory', -25: 'argument'}

P = ctypes.POINTER
vp, i64, c_int, c_double = ctypes.c_void_p, ctypes.c_int64, ctypes.c_int, ctypes.c_double


class KcAd(ctypes.Structure):
    """kc_ad: an interval value and its gradient in theta and the face's two local parameters."""
    _fields_ = [('v', carith.KcIv), ('g', carith.KcIv * 3)]


_LIBS = {}


def lib(path=None):
    """The library with the model and leaf functions bound, loaded once per path.  path=None is
    build/libkisscheck.so through carith.load(); another path (a test's mutant build) is loaded and initialised
    the same way (kc_init must return 0)."""
    key = str(path) if path is not None else None
    if key in _LIBS: return _LIBS[key]
    if path is None:
        L = carith.load()
    else:
        L = ctypes.CDLL(key)
        L.kc_init.restype = c_int
        rc = L.kc_init()
        if rc != 0: raise OSError('kc_init failed with code %d' % rc)
    L.kc_graph_new.restype = vp
    L.kc_graph_new.argtypes = [c_int, P(c_int), P(c_int), ctypes.c_char_p, P(c_int)]
    L.kc_graph_free.argtypes = [vp]
    for f in ('kc_graph_n', 'kc_graph_npar', 'kc_graph_nedges'):
        getattr(L, f).restype = c_int; getattr(L, f).argtypes = [vp]
    L.kc_leaf_ws_new.restype = vp; L.kc_leaf_ws_new.argtypes = [vp]
    L.kc_leaf_ws_free.argtypes = [vp]
    L.kc_leaf_invalid_vsum.restype = c_int; L.kc_leaf_invalid_vsum.argtypes = [vp, vp, vp, c_int, P(c_int)]
    L.kc_leaf_context.restype = c_int; L.kc_leaf_context.argtypes = [vp, vp, vp]
    L.kc_leaf_hint.restype = c_int
    L.kc_leaf_hint.argtypes = [vp, vp, vp, vp, vp, i64, P(c_int), P(carith.KcIv)]
    L.kc_leaf_exclude.restype = c_int; L.kc_leaf_exclude.argtypes = [vp, vp, vp, c_int, c_int, P(c_int), vp]
    L.kc_check_chart.restype = c_int
    L.kc_check_chart.argtypes = [vp, c_int, vp, vp, vp, i64, vp, i64, vp, i64, c_int, c_double, c_int, c_int,
                                 vp, vp, i64, P(i64), vp, vp, P(c_int)]
    _LIBS[key] = L
    return L


def _dbl(x):
    """A contiguous float64 array (the kernel reads doubles)."""
    return np.ascontiguousarray(np.asarray(x, dtype=np.float64))


class CGraph:
    """The C graph (kc_graph_new) of a candidate's faces and a chart's pentagon modes, with a workspace for
    the per-box functions.  .err is the model's error code (0 on success)."""

    def __init__(self, faces, modes, path=None):
        """Build the graph in the library at path (None: the default build)."""
        self.L = lib(path); self.h = None; self.ws = None
        nf = len(faces)
        sizes = (c_int * nf)(*[len(f) for f in faces])
        flat = [v for f in faces for v in f]
        verts = (c_int * max(1, len(flat)))(*flat)
        md = bytes(ord(modes.get(str(i), modes.get(i, '\0'))) for i in range(nf))
        err = c_int(0)
        self.h = self.L.kc_graph_new(nf, sizes, verts, md, ctypes.byref(err))
        self.err = err.value
        if self.h:
            self.n = self.L.kc_graph_n(self.h); self.npar = self.L.kc_graph_npar(self.h)
            self.ws = self.L.kc_leaf_ws_new(self.h)

    def __del__(self):
        """Free the workspace and the graph."""
        if getattr(self, 'ws', None): self.L.kc_leaf_ws_free(self.ws)
        if getattr(self, 'h', None): self.L.kc_graph_free(self.h)

    # ---- per-box functions (the differential test calls these one by one)
    def invalid_vsum(self, lo, hi, tighten):
        """(rc, rule): kc_leaf_invalid_vsum; rule 0, 1 'invalid', 2 'vsum'."""
        lo, hi = _dbl(lo), _dbl(hi); rule = c_int(-1)
        rc = self.L.kc_leaf_invalid_vsum(self.ws, lo.ctypes.data, hi.ctypes.data, int(tighten), ctypes.byref(rule))
        return rc, rule.value

    def context(self, Blo, Bhi):
        """rc of kc_leaf_context on B = [Blo, Bhi]."""
        Blo, Bhi = _dbl(Blo), _dbl(Bhi)
        return self.L.kc_leaf_context(self.ws, Blo.ctypes.data, Bhi.ctypes.data)

    def hint(self, lo, hi, mu):
        """(rc, closed, s): kc_leaf_hint with the hint given as parallel row and value sequences (rows,
        float32 values), in the order of the q and mu streams."""
        rows, vals = mu
        lo, hi = _dbl(lo), _dbl(hi)
        r = np.ascontiguousarray(np.asarray(rows, dtype=np.uint16))
        v = np.ascontiguousarray(np.asarray(vals, dtype=np.float32))
        closed = c_int(-1); s = carith.KcIv(0.0, 0.0)
        rc = self.L.kc_leaf_hint(self.ws, lo.ctypes.data, hi.ctypes.data, r.ctypes.data, v.ctypes.data, len(r),
                                 ctypes.byref(closed), ctypes.byref(s))
        return rc, closed.value, (s.lo, s.hi)

    def exclude(self, lo, hi, max_depth=MAX_DEPTH, max_pieces=MAX_PIECES):
        """(rc, verdict, counts): kc_leaf_exclude."""
        lo, hi = _dbl(lo), _dbl(hi); v = c_int(-1); cnt = np.zeros(NCOUNT, np.int64)
        rc = self.L.kc_leaf_exclude(self.ws, lo.ctypes.data, hi.ctypes.data, max_depth, max_pieces,
                                    ctypes.byref(v), cnt.ctypes.data)
        return rc, v.value, dict(zip(COUNTS, cnt.tolist()))


def streams(chart):
    """The decoded event, q and mu streams of a chart and whether it has a mu stream (certificate.unpack)."""
    ev = np.ascontiguousarray(unpack(chart['events'], np.uint8), dtype=np.uint8)
    q = np.ascontiguousarray(unpack(chart['q'], np.uint16), dtype=np.uint16)
    has_mu = 'mu' in chart
    mu = np.ascontiguousarray(unpack(chart['mu'], np.float32) if has_mu else np.zeros(0), dtype=np.float32)
    return ev, q, mu, has_mu


def check_chart(cg, chart, tol, max_depth=MAX_DEPTH, max_pieces=MAX_PIECES, timing=False, decoded=None):
    """Run kc_check_chart on one chart with the graph cg.  Returns a dict: rc (0 or an error code), error
    (its name), model_err, verdict and kind (int8 arrays, one entry per unit decided), counts (by the
    reference's names), unit_ns (int64 array, if timing), seconds (wall time of the C call)."""
    L = cg.L
    ev, q, mu, has_mu = decoded if decoded is not None else streams(chart)
    lo0, hi0 = _dbl(chart['lo0']), _dbl(chart['hi0'])
    n = len(lo0)
    cap = 2 * len(ev) + 1                    # a leaf event gives one unit, a contraction event at most two
    verdict = np.zeros(cap, np.int8); kind = np.zeros(cap, np.int8); counts = np.zeros(NCOUNT, np.int64)
    ns = np.zeros(cap, np.int64) if timing else None
    nunits = i64(0); merr = c_int(0)
    t0 = time.perf_counter()
    rc = L.kc_check_chart(cg.h, n, lo0.ctypes.data, hi0.ctypes.data, ev.ctypes.data, len(ev), q.ctypes.data,
                          len(q), mu.ctypes.data, len(mu), int(has_mu), float(tol), max_depth, max_pieces,
                          verdict.ctypes.data, kind.ctypes.data, cap, ctypes.byref(nunits), counts.ctypes.data,
                          ns.ctypes.data if timing else None, ctypes.byref(merr))
    dt = time.perf_counter() - t0
    k = nunits.value
    return dict(rc=rc, error=ERRORS.get(rc) if rc else None, model_err=merr.value, verdict=verdict[:k],
                kind=kind[:k], counts=dict(zip(COUNTS, counts.tolist())), unit_ns=ns[:k] if timing else None,
                seconds=dt)


def _merge_counts(dst, src):
    """Add the counters src into dst."""
    for k, v in src.items():
        if v: dst[k] = dst.get(k, 0) + v


def _fallback(cert, ch, G, res, out, kstats):
    """Walk the reference replay of the chart alongside the kernel's verdicts; make the reference's own calls
    for every unclosed unit, and record terminal and open leaves as check_certificate does."""
    K = ck.Checker(G); verdict = res['verdict']; u = 0
    for kind, lo, hi, extra in replay(ch):
        lo = [float(x) for x in lo]; hi = [float(x) for x in hi]
        if kind == 'split': continue
        if kind == 'contract':
            removed, _ = extra; lctx = None
            for sl in removed:
                v = int(verdict[u]); u += 1
                if v != 0: continue
                kstats['fallback'] += 1
                if lctx is None:                                   # check_certificate's two lines, as there
                    CB, _ = G.corners(lo, hi, tighten=False)
                    lctx = K.linear_context(lo, hi, CB)
                slo, shi, clo, chi, k, side, mu = sl
                sl = ([float(x) for x in slo], [float(x) for x in shi], [float(x) for x in clo],
                      [float(x) for x in chi], int(k), int(side), mu)
                if K.exclude_slab(sl, (lo, hi, CB, lctx)):
                    kstats['fallback_closed'] += 1
                else:
                    out['failed'].append(dict(kind='slab', lo=sl[0], hi=sl[1])); out['verified'] = False
            continue
        v = int(verdict[u]); u += 1
        out['leaves'] += 1
        if kind == 'terminal':
            if v == 8:
                out['failed'].append(dict(kind='terminal-too-wide', lo=lo, hi=hi)); out['verified'] = False
            out['terminal'].append(dict(modes=ch['modes'], lo=lo, hi=hi)); continue
        if kind == 'open':
            out['open'] += 1; out['verified'] = False; continue
        if v != 0: continue
        kstats['fallback'] += 1
        if kind == 'linear' and extra:
            C_, _ = G.corners(lo, hi, tighten=False)
            if K.test_linear_mu(lo, hi, lo, hi, C_, extra):
                K.bump('linear-hint'); kstats['fallback_closed'] += 1; continue
        if K.exclude(lo, hi):
            kstats['fallback_closed'] += 1
        else:
            out['failed'].append(dict(kind=kind, lo=lo, hi=hi)); out['verified'] = False
    if u != len(verdict):                       # never `assert`: it would vanish under python -O
        raise RuntimeError('kernel and reference disagree on the number of units (%d, %d)' % (len(verdict), u))
    _merge_counts(out['counts'], K.stats)


def check_certificate(cert, code=None, gen_faces=None, r_min=ck.R_MIN, structural=True,
                      max_depth=MAX_DEPTH, max_pieces=MAX_PIECES, path=None):
    """checker.check_certificate with the boxes split between the C kernel and the Python reference.  Returns
    the reference's dict (verified, structural, failed, terminal, open, not_run, counts, leaves) plus
    'kernel': units, closed by C, fallback units and how many of them the reference closed, charts the
    kernel stopped on (whole chart to the reference), C seconds and fallback seconds.  path: the library to
    use (None: build/libkisscheck.so; the mutation tests pass another build)."""
    out = dict(verified=True, structural=[], failed=[], terminal=[], open=0, not_run=0, counts={}, leaves=0)
    kstats = dict(units=0, closed_c=0, fallback=0, fallback_closed=0, chart_errors=[], c_seconds=0.0,
                  fallback_seconds=0.0)
    out['kernel'] = kstats
    if cert.get('version') not in (1, 2): raise ValueError('unknown certificate version')
    if structural:
        out['structural'] = ck.structural_check(cert, code, gen_faces, r_min)
        if out['structural']: out['verified'] = False
    for ch in cert['charts']:
        if ch.get('not_run'):
            out['not_run'] += 1; out['verified'] = False; continue
        G = am.Graph(cert['faces'], ch['modes'])                    # raises as the reference does
        if G.npar != len(ch['lo0']): raise ValueError('parameter count mismatch')
        if not (abs(ch['lo0'][0] - cert['theta'][0]) == 0 and abs(ch['hi0'][0] - cert['theta'][1]) == 0):
            raise ValueError('starting box does not span the stated theta range')
        cg = CGraph(cert['faces'], ch['modes'], path)
        res = check_chart(cg, ch, float(cert['tol']), max_depth, max_pieces) if cg.err == 0 else None
        if res is not None: kstats['c_seconds'] += res['seconds']
        t0 = time.perf_counter()
        if res is None or res['rc'] != 0:
            kstats['chart_errors'].append(dict(modes=ch['modes'], error=(res or {}).get('error', 'graph %d' % cg.err)))
            one = ck.check_certificate(dict(cert, charts=[ch]), structural=False)   # raises like the reference
            for key in ('failed', 'terminal'): out[key] += one[key]
            for key in ('open', 'not_run', 'leaves'): out[key] += one[key]
            if not one['verified']: out['verified'] = False
            _merge_counts(out['counts'], one['counts'])
        else:
            v = res['verdict']
            kstats['units'] += len(v); kstats['closed_c'] += int(np.isin(v, CLOSED).sum())
            _merge_counts(out['counts'], res['counts'])
            if np.isin(v, CLOSED).all():
                out['leaves'] += int((res['kind'] != 0).sum())
            else:
                _fallback(cert, ch, G, res, out, kstats)
        kstats['fallback_seconds'] += time.perf_counter() - t0
    return out
