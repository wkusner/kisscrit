"""Independent checker for completeness certificates (format: docs/certificate-format.md).

For each chart of a certificate the checker replays the proof tree (kisscrit.core.certificate.replay) and proves,
with its own interval evaluation of the model (kisscrit.core.model), that every excluded leaf and every region
removed by a contraction contains no critical configuration of the candidate graph.  It never trusts the
search: the search's reasons are only hints about which test to try first, multipliers and Gordan vectors
are recomputed here, and a box that no test excludes is bisected and the pieces are tested again, down to a
fixed depth.  A box that still fails is reported, never accepted.

The exclusion tests (each proves: no point of the box is a valid realization of the graph that carries a
strictly positive balanced stress):

  invalid   some face condition fails on the whole box (no convex equilateral face in this chart, a corner
            outside [alpha(theta), pi], a diagonal shorter than theta, a pentagon region condition);
  vsum      some vertex angle sum differs from 2 pi on the whole box;
  linear    for multipliers mu (any real vector), 0 is not in  mu.F(m) + sum_j (mu^T J(B))_j (B_j - m_j),
            where F = vertex sums - 2 pi, m is the centre of a box B containing the tested box, and J(B)
            encloses the (generalized) Jacobian of F on B.  By the mean value inclusion, F has no zero there;
  gordan    a vector y with  y_u . t_u(e) + y_w . t_w(e) > 0  for every edge e = uw on the whole box, where
            t_v(e) is the unit tangent at v towards the other end of e; then no nonzero nonnegative
            balanced stress exists (Gordan's alternative);
  motzkin   every edge term is exactly 0 (y_u = y_w = 0) or > 0 on the box, and some term is > 0; then no
            strictly positive balanced stress exists (Motzkin's alternative).

Floating point enters only through linear programs that choose mu and y; the inequalities themselves are
checked in ball arithmetic.
"""
import numpy as np
from scipy.optimize import linprog
from . import model as am
from . import interval as I
from .interval import Iv
from .certificate import replay

MAX_DEPTH = 40          # bisection depth allowed below a certificate leaf
MAX_PIECES = 400        # bisection pieces allowed per certificate leaf or slab before it is reported as failed


def _f(x):
    """Float midpoint of an interval (used only to set up linear programs)."""
    return 0.5 * (x.lo + x.hi)


def _r(x): return 0.5 * (x.hi - x.lo)


def combination(mu, Fm, JB, Y):
    """The interval  mu.F(m) + sum_j r_j Y_j  with  r_j = sum_v mu_v J_vj  (sum first, then multiply)."""
    n = len(Y)
    s = Iv(0.0); r = [Iv(0.0)] * n
    for v, val in mu:
        mv = Iv(float(val))
        s = s + mv * Fm[v]
        for j in range(n):
            r[j] = r[j] + mv * JB[v][j]
    for j in range(n):
        s = s + r[j] * Y[j]
    return s


class Checker:
    def __init__(self, graph):
        self.G = graph
        self.stats = {}

    def bump(self, k): self.stats[k] = self.stats.get(k, 0) + 1

    # ------------------------------------------------------------------ the tests
    def evaluate(self, lo, hi, blo=None, bhi=None):
        C, ok = self.G.corners(lo, hi, blo, bhi)
        return C, ok

    def test_invalid_vsum(self, C, ok):
        if not ok: return 'invalid'
        for s, g in self.G.vertex_sums(C):
            if s.lo > 0 or s.hi < 0: return 'vsum'
        return None

    def linear_context(self, Blo, Bhi, CB):
        """Centre m of B, F(m) (clamp bounds of B) and J(B), computed once per box B."""
        G = self.G
        m = [0.5 * (a + b) for a, b in zip(Blo, Bhi)]
        Cm, _ = G.corners(m, m, Blo, Bhi)
        Fm = [s for s, g in G.vertex_sums(Cm)]
        JB = [g for s, g in G.vertex_sums(CB)]
        return (Blo, Bhi, m, Fm, JB)

    def test_linear(self, lo, hi, Blo, Bhi, CB, lctx=None):
        """No zero of the vertex equations in [lo, hi], a sub-box of B = [Blo, Bhi]; CB = corners on B."""
        G = self.G; n = G.npar
        _, _, m, Fm, JB = lctx if lctx is not None else self.linear_context(Blo, Bhi, CB)
        rows = [v for v in range(G.n) if Fm[v].finite() and all(x.finite() for x in JB[v])]
        if not rows: return False
        Y = [Iv(lo[j]) - Iv(m[j]) if lo[j] == hi[j] else I.hull(Iv(lo[j]) - Iv(m[j]), Iv(hi[j]) - Iv(m[j]))
             for j in range(n)]
        ylo = np.array([y.lo for y in Y]); yhi = np.array([y.hi for y in Y])
        Jc = np.array([[_f(JB[v][j]) for j in range(n)] for v in rows])
        Jr = np.array([[_r(JB[v][j]) for j in range(n)] for v in rows])
        fc = np.array([_f(Fm[v]) for v in rows]); fr = np.array([_r(Fm[v]) for v in rows])
        R = fr + Jr @ np.maximum(np.abs(ylo), np.abs(yhi))
        V = len(rows)
        # min t  s.t.  Jc y - t <= R - fc,  -Jc y - t <= R + fc,  y in [ylo, yhi], t >= 0
        A = np.vstack([np.hstack([Jc, -np.ones((V, 1))]), np.hstack([-Jc, -np.ones((V, 1))])])
        b = np.concatenate([R - fc, R + fc])
        c = np.zeros(n + 1); c[-1] = 1
        res = linprog(c, A_ub=A, b_ub=b, bounds=list(zip(ylo, yhi)) + [(0, None)], method='highs')
        if res.status != 0 or res.fun <= 0: return False
        lam = -res.ineqlin.marginals
        mu = lam[:V] - lam[V:]
        if not np.any(mu): return False
        s = combination([(v, mu[i]) for i, v in enumerate(rows) if mu[i] != 0], Fm, JB, Y)
        return s.lo > 0 or s.hi < 0

    def test_linear_bound(self, lo, hi, cut_lo, cut_hi, k, side, Blo, Bhi, CB, lctx=None):
        """Slab [lo, hi] = part of the box [cut_lo, cut_hi] (inside B) beyond the new bound of coordinate k.
        The multipliers come from the same floating-point simplex code the search uses (min or max y_k over
        the linearization on B); they are only a choice of mu.  The proof is the interval inclusion: with
        y = p - m (m the centre of B), 0 is not in mu.F(m) + sum_j (mu^T J(B))_j Y_j over the slab."""
        from ..completeness import fast           # used only as a source of multipliers
        G = self.G; n = G.npar; V = G.n
        _, _, m, Fm, JB = lctx if lctx is not None else self.linear_context(Blo, Bhi, CB)
        mask = [Fm[v].finite() and all(x.finite() for x in JB[v]) for v in range(V)]
        if not any(mask): return False
        Jc = np.zeros((V, n)); R = np.zeros(V); fc = np.zeros(V)
        w = np.array(Bhi) - np.array(Blo)
        for v in range(V):
            if mask[v]:
                fc[v] = _f(Fm[v]); Jc[v] = [_f(x) for x in JB[v]]
                R[v] = (_r(Fm[v]) + sum(_r(JB[v][j]) * 0.5 * w[j] for j in range(n))) * (1 + 1e-9) + 1e-14
            else:
                R[v] = 1e30
        ylo = np.array(Blo) - np.array(m); yhi = np.array(Bhi) - np.array(m)
        T, basis = fast.lp_setup(Jc, fc, R, ylo, yhi - ylo)
        c = np.zeros(n + 1); c[n] = 1.0
        fast.lp_set_objective(T, basis, c)
        if fast._simplex(T, basis, 400) != 0: return False
        c = np.zeros(n + 1); c[k] = 1.0 if side < 0 else -1.0; c[n] = 1e3
        fast.lp_set_objective(T, basis, c)
        if fast._simplex(T, basis, 400) != 0: return False
        mu = fast.lp_multipliers(T, V, n)
        for v in range(V):
            if not mask[v]: mu[v] = 0.0
        if not np.any(mu): return False
        Y = [I.hull(Iv(lo[j]) - Iv(m[j]), Iv(hi[j]) - Iv(m[j])) for j in range(n)]
        s = combination([(v, mu[v]) for v in range(V) if mu[v] != 0], Fm, JB, Y)
        return s.lo > 0 or s.hi < 0

    def _edge_terms(self, C, y):
        G = self.G; d = G.directions(C); out = []
        for (u, w) in G.edges:
            t = Iv(0.0); zero = True
            for (a, b) in ((u, w), (w, u)):
                ya, yb = y[2 * a], y[2 * a + 1]
                if ya == 0 and yb == 0: continue
                zero = False
                ang = d[(a, b)]
                t = t + Iv(float(ya)) * I.cos(ang) + Iv(float(yb)) * I.sin(ang)
            out.append((zero, t))
        return out

    def _lp_vector(self, lo, hi, motzkin):
        G = self.G
        m = [0.5 * (a + b) for a, b in zip(lo, hi)]
        Cm, _ = G.corners(m, m, lo, hi)
        d = G.directions(Cm)
        E = len(G.edges); n2 = 2 * G.n
        A = np.zeros((E, n2))
        for k, (u, w) in enumerate(G.edges):
            for (a, b) in ((u, w), (w, u)):
                ang = _f(d[(a, b)]); A[k, 2 * a] = np.cos(ang); A[k, 2 * a + 1] = np.sin(ang)
        if not motzkin:     # max t: A y >= t
            Aub = np.hstack([-A, np.ones((E, 1))]); c = np.zeros(n2 + 1); c[-1] = -1
            res = linprog(c, A_ub=Aub, b_ub=np.zeros(E), bounds=[(-1, 1)] * n2 + [(None, 1)], method='highs')
            if res.status != 0 or res.x[-1] <= 1e-12: return None
            return res.x[:n2]
        Aub = np.hstack([-A, np.eye(E)]); c = np.concatenate([np.zeros(n2), -np.ones(E)])
        res = linprog(c, A_ub=Aub, b_ub=np.zeros(E), bounds=[(-1, 1)] * n2 + [(0, 1)] * E, method='highs')
        if res.status != 0 or -res.fun <= 1e-9: return None
        y = res.x[:n2].copy()
        for v in range(G.n):
            if abs(y[2 * v]) + abs(y[2 * v + 1]) < 1e-7: y[2 * v] = y[2 * v + 1] = 0.0
        return y

    def test_gordan(self, lo, hi, C):
        y = self._lp_vector(lo, hi, False)
        if y is None: return False
        return all((not z) and t.lo > 0 for z, t in self._edge_terms(C, y))

    def test_motzkin(self, lo, hi, C):
        y = self._lp_vector(lo, hi, True)
        if y is None: return False
        terms = self._edge_terms(C, y); npos = 0
        for z, t in terms:
            if z: continue
            if not t.lo > 0: return False
            npos += 1
        return npos > 0

    # ------------------------------------------------------------------ one box
    def test_linear_mu(self, lo, hi, Blo, Bhi, CB, mu, lctx=None):
        """The interval inclusion of test_linear with given multipliers mu ({vertex: value}, a hint)."""
        if not mu: return False
        G = self.G; n = G.npar
        _, _, m, Fm, JB = lctx if lctx is not None else self.linear_context(Blo, Bhi, CB)
        Y = [I.hull(Iv(lo[j]) - Iv(m[j]), Iv(hi[j]) - Iv(m[j])) for j in range(n)]
        s = combination(sorted(mu.items()), Fm, JB, Y)
        return s.lo > 0 or s.hi < 0

    def exclude_slab(self, slab, parent):
        """A region removed by a contraction: the search's multiplier hint first, then the other tests."""
        lo, hi, cut_lo, cut_hi, k, side, mu = slab
        Blo, Bhi, CB, lctx = parent
        if self.test_linear_mu(lo, hi, Blo, Bhi, CB, mu, lctx):
            self.bump('linear-hint'); return True
        if self.test_linear_bound(lo, hi, cut_lo, cut_hi, k, side, Blo, Bhi, CB, lctx):
            self.bump('linear-bound'); return True
        return self.exclude(lo, hi, parent=(Blo, Bhi, CB, lctx))

    def exclude(self, lo, hi, parent=None, depth=0):
        """Prove that [lo, hi] contains no critical configuration; bisect if needed.  Returns True/False."""
        if depth == 0: self.pieces = 0
        self.pieces += 1
        if self.pieces > MAX_PIECES: return False
        C, ok = self.evaluate(lo, hi)
        r = self.test_invalid_vsum(C, ok)
        if r: self.bump(r); return True
        if parent is not None:
            Blo, Bhi, CB, lctx = parent
            if self.test_linear(lo, hi, Blo, Bhi, CB, lctx): self.bump('linear-parent'); return True
        if self.test_linear(lo, hi, lo, hi, C): self.bump('linear'); return True
        if self.test_gordan(lo, hi, C): self.bump('gordan'); return True
        if self.test_motzkin(lo, hi, C): self.bump('motzkin'); return True
        if depth >= MAX_DEPTH: return False
        w = [b - a for a, b in zip(lo, hi)]; k = int(np.argmax(w))
        if w[k] <= 0: return False
        mid = 0.5 * (lo[k] + hi[k])
        lo1, hi1 = list(lo), list(hi); hi1[k] = mid
        lo2, hi2 = list(lo), list(hi); lo2[k] = mid
        self.bump('bisect')
        return self.exclude(lo1, hi1, None, depth + 1) and self.exclude(lo2, hi2, None, depth + 1)


def check_certificate(cert):
    """Check one certificate.  Returns a dict: verified (bool), per-reason counts, failed boxes, terminal boxes,
    number of open leaves and of charts not run."""
    out = dict(verified=True, failed=[], terminal=[], open=0, not_run=0, counts={}, leaves=0)
    if cert.get('version') not in (1, 2): raise ValueError('unknown certificate version')
    for ch in cert['charts']:
        if ch.get('not_run'):
            out['not_run'] += 1; out['verified'] = False; continue
        G = am.Graph(cert['faces'], ch['modes'])
        if G.npar != len(ch['lo0']): raise ValueError('parameter count mismatch')
        if not (abs(ch['lo0'][0] - cert['theta'][0]) == 0 and abs(ch['hi0'][0] - cert['theta'][1]) == 0):
            raise ValueError('starting box does not span the stated theta range')
        K = Checker(G)
        for kind, lo, hi, extra in replay(ch):
            lo = [float(x) for x in lo]; hi = [float(x) for x in hi]
            if kind == 'split': continue
            if kind == 'contract':
                removed, _ = extra
                CB, _ = G.corners(lo, hi)
                lctx = K.linear_context(lo, hi, CB)
                for sl in removed:
                    slo, shi, clo, chi, k, side, mu = sl
                    sl = ([float(x) for x in slo], [float(x) for x in shi], [float(x) for x in clo],
                          [float(x) for x in chi], int(k), int(side), mu)
                    if not K.exclude_slab(sl, (lo, hi, CB, lctx)):
                        out['failed'].append(dict(kind='slab', lo=sl[0], hi=sl[1])); out['verified'] = False
                continue
            out['leaves'] += 1
            if kind == 'terminal':
                out['terminal'].append(dict(modes=ch['modes'], lo=lo, hi=hi)); continue
            if kind == 'open':
                out['open'] += 1; out['verified'] = False; continue
            if kind == 'linear' and extra:
                C_, _ = G.corners(lo, hi)
                if K.test_linear_mu(lo, hi, lo, hi, C_, extra):
                    K.bump('linear-hint'); continue
            if not K.exclude(lo, hi):
                out['failed'].append(dict(kind=kind, lo=lo, hi=hi)); out['verified'] = False
        for k, v in K.stats.items(): out['counts'][k] = out['counts'].get(k, 0) + v
    return out
