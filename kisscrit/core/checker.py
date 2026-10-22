"""Independent checker for completeness certificates (format: docs/certificate-format.md).

Status: trusted. Part of the trusted base (transfer prompt v1.2, code standards).
Implements: the exclusion tests of docs/lemmas.md L15-L16 (necessary conditions, Gordan and Motzkin alternatives, mean-value inclusion) and the structural checks V1-V4 of the checker spec v0.2; paper Appendix A (candidate list) and Appendix D (certification).

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

Structural checks (spec v0.2, V1-V4), done before any box is examined:
  V1  the checker derives its own a-priori starting box per chart (Arb, outward rounded) from the lemmas and
      requires the certificate's box to contain it;
  V2  the checker recomputes which pentagon regions are required (Lemma P and L11) and requires each to be
      present and complete;
  V3  the theta range is tied to the claim: theta_lo <= 2 asin(r_min / (1 + r_min)), theta_hi >= arctan 2;
      and it lies in the model's domain: finite, THETA_FLOOR <= theta_lo < theta_hi <= THETA_CEIL (L17);
  V4  a checker-owned canonical code of the faces' rotation system equals the certificate's code, and the
      faces match the gen record up to cyclic rotation when the gen faces are supplied.
Dependencies: kisscrit.completeness is imported only inside test_linear_bound, as a source of multipliers
(any real mu is admissible in the interval inclusion); tests/test_checker_imports.py enforces this.
"""
import math
import numpy as np
from scipy.optimize import linprog
from flint import arb, fmpq
from . import model as am
from . import interval as I
from .interval import Iv, down, up
from .certificate import replay

R_MIN = fmpq(10001, 10000)      # the route constant r_min = 1.0001, exactly

MAX_DEPTH = 40          # bisection depth allowed below a certificate leaf
MAX_PIECES = 400        # bisection pieces allowed per certificate leaf or slab before it is reported as failed


def _f(x):
    """Float midpoint of an interval (used only to set up linear programs)."""
    return 0.5 * (x.lo + x.hi)


def _r(x):
    """Float radius of an interval (used only to set up linear programs)."""
    return 0.5 * (x.hi - x.lo)


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
        """A checker for one candidate graph (a kisscrit.core.model.Graph) with per-reason counters."""
        self.G = graph
        self.stats = {}

    def bump(self, k):
        """Count one exclusion of kind k."""
        self.stats[k] = self.stats.get(k, 0) + 1

    # ------------------------------------------------------------------ the tests
    def evaluate(self, lo, hi, blo=None, bhi=None):
        """Corner values with gradients on the box [lo, hi] (clamp bounds from [blo, bhi]), with monotonicity tightening."""
        C, ok = self.G.corners(lo, hi, blo, bhi)
        return C, ok

    def test_invalid_vsum(self, C, ok):
        """Claim: returns 'invalid' only if no point of the box describes valid faces in this chart (the model's ok flag, L3/L6-L8/L10), and 'vsum' only if some vertex angle sum differs from 2 pi on the whole box (the interval s = sum - 2 pi excludes 0), so no realization exists there (L15: the vertex sums of a realization are exactly 2 pi; L3, L6-L8 for the face conditions)."""
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
        if res.status == 0 and res.fun > 0:
            lam = -res.ineqlin.marginals
            mu = lam[:V] - lam[V:]
            if np.any(mu):
                s = combination([(v, mu[i]) for i, v in enumerate(rows) if mu[i] != 0], Fm, JB, Y)
                if s.lo > 0 or s.hi < 0: return True
        return self.test_linear_combo(rows, Fm, JB, Y, Jc, Jr, fc, fr, ylo, yhi)

    def test_linear_combo(self, rows, Fm, JB, Y, Jc, Jr, fc, fr, ylo, yhi):
        """A second choice of mu, by the 'sum first, then multiply' form the search uses: maximize over mu
        (with sum |mu_v| <= 1) the floating-point margin
            sign * sum_v mu_v fc_v - sum_j |sum_v mu_v Jc_vj| ymax_j - sum_v |mu_v| (fr_v + sum_j Jr_vj ymax_j),
        for both signs, and verify the result with the exact interval combination.  This exploits cancellation
        in sum_v mu_v J_vj across rows, which the row-wise LP of test_linear cannot."""
        V = len(rows); n = Jc.shape[1]
        ymax = np.maximum(np.abs(ylo), np.abs(yhi))
        # variables: mu+ (V), mu- (V), s (n) ; mu = mu+ - mu-, |mu| <= mu+ + mu-, s_j >= |sum_v mu_v Jc_vj|
        nv = 2 * V + n
        A_ub = []; b_ub = []
        for j in range(n):
            row = np.zeros(nv); row[:V] = Jc[:, j]; row[V:2 * V] = -Jc[:, j]; row[2 * V + j] = -1; A_ub.append(row); b_ub.append(0.0)
            row = np.zeros(nv); row[:V] = -Jc[:, j]; row[V:2 * V] = Jc[:, j]; row[2 * V + j] = -1; A_ub.append(row); b_ub.append(0.0)
        row = np.zeros(nv); row[:2 * V] = 1; A_ub.append(row); b_ub.append(1.0)
        pen = fr + Jr @ ymax                                   # per-row radius cost of |mu_v|
        for sign in (1.0, -1.0):
            c = np.zeros(nv)
            c[:V] = -(sign * fc) + pen; c[V:2 * V] = (sign * fc) + pen; c[2 * V:] = ymax
            res = linprog(c, A_ub=np.array(A_ub), b_ub=np.array(b_ub), bounds=[(0, None)] * nv, method='highs')
            if res.status != 0 or -res.fun <= 0: continue
            mu = res.x[:V] - res.x[V:2 * V]
            mu[np.abs(mu) < 1e-12 * np.abs(mu).max()] = 0.0
            if not np.any(mu): continue
            s = combination([(v, mu[i]) for i, v in enumerate(rows) if mu[i] != 0], Fm, JB, Y)
            if s.lo > 0 or s.hi < 0: return True
        return False

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
        """For each edge e = uw and a vector y (2 components per vertex), the interval enclosure of y_u . t_u(e) + y_w . t_w(e), with t_v(e) the unit tangent at v towards the other end, computed from the edge directions of the box in the per-vertex frames of model.Graph.directions (first neighbour at angle 0; the balance condition and both alternatives are invariant under the frame choice at each vertex, L16); `zero` marks edges whose two end vectors are exactly 0."""
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
        """A floating-point choice of y for the Gordan (motzkin=False) or Motzkin test at the centre of the box; None if the LP finds no candidate. Only a choice: the proof is the interval test that follows."""
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
        """Claim: returns True only if some y has y_u . t_u(e) + y_w . t_w(e) > 0 for every edge on the whole box, so by Gordan's alternative (L16) no nonzero nonnegative balanced stress exists at any point of the box."""
        y = self._lp_vector(lo, hi, False)
        if y is None: return False
        return all((not z) and t.lo > 0 for z, t in self._edge_terms(C, y))

    def test_motzkin(self, lo, hi, C):
        """Claim: returns True only if some y makes every edge term either exactly 0 (both end vectors zero) or > 0 on the whole box, with at least one > 0; by Motzkin's alternative (L16) no strictly positive balanced stress exists at any point of the box."""
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
        # cheap pass: corners without monotonicity tightening (gradients are unaffected, values are wider)
        C, ok = self.G.corners(lo, hi, tighten=False)
        r = self.test_invalid_vsum(C, ok)
        if r: self.bump(r); return True
        if parent is not None:
            Blo, Bhi, CB, lctx = parent
            if self.test_linear(lo, hi, Blo, Bhi, CB, lctx): self.bump('linear-parent'); return True
        if self.test_linear(lo, hi, lo, hi, C): self.bump('linear'); return True
        # full pass: tightened corner values for the tests that use them
        if any(k == 5 for k in self.G.sizes):
            C, ok = self.evaluate(lo, hi)
            r = self.test_invalid_vsum(C, ok)
            if r: self.bump(r + '-tight'); return True
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


# ---------------------------------------------------------------------- structural checks (V1-V4)
def _alpha(t):
    """alpha(t) = arccos(cos t / (1 + cos t)) as an Arb ball (the equilateral corner, docs/lemmas.md)."""
    c = arb(t).cos()
    return (c / (1 + c)).acos()


def theta_domain_ok(theta):
    """V3, domain part: theta_lo and theta_hi are finite and THETA_FLOOR <= theta_lo < theta_hi <= THETA_CEIL,
    decided in Arb on the certificate's doubles.  Every enclosure of the model and every bound of apriori_box
    is valid only for theta in [THETA_FLOOR, THETA_CEIL] (docs/lemmas.md L17, constants in model.py); outside
    it alpha is not monotone, alpha_range can be inverted, and a box can be excluded on that."""
    tlo, thi = (float(t) for t in theta)
    if not (math.isfinite(tlo) and math.isfinite(thi)): return False
    lo, hi = arb(tlo), arb(thi)
    return bool(arb(am.THETA_FLOOR) <= lo) and bool(lo < hi) and bool(hi <= arb(am.THETA_CEIL))


def theta_claim_ok(theta, r_min=R_MIN):
    """V3: [theta_lo, theta_hi] lies in the model's domain (theta_domain_ok, L17) and covers the claimed range,
    decided in Arb.  Two spheres of radius r touching a unit sphere have centres at angular distance
    2 asin(r / (1 + r)) when they touch each other, so r >= r_min is theta >= 2 asin(r_min / (1 + r_min)); the
    icosahedral angle arctan 2 is the largest minimum distance of 12 points (L9, uniqueness at theta_max)."""
    if not theta_domain_ok(theta): return False
    tlo, thi = theta
    r = arb(r_min)
    tmin = 2 * (r / (1 + r)).asin()
    tmax = arb(2).atan()
    return bool(arb(float(tlo)) <= tmin) and bool(arb(float(thi)) >= tmax)


def apriori_box(G, theta):
    """V1: the a-priori box of a chart, from the lemmas, rounded outward: theta in [theta_lo, theta_hi]; a
    rhombus corner in [alpha(theta_lo), 2 alpha(theta_hi)] (L7); a pentagon corner in [alpha(theta_lo), pi]
    (L3); in regions b and c the corner at boundary position 4 is >= LSTRAIGHT (L10'); at a degree-2 vertex
    that corner is exactly pi, enclosed by [down(pi), up(pi)] (L4).  Returns (lo, hi, fixed) with fixed the
    set of coordinates that are straight corners at degree-2 vertices.  The certificate's box must contain
    this box coordinatewise, so every valid configuration of the region lies in the chart's starting box."""
    tlo, thi = float(theta[0]), float(theta[1])
    a_lo = down(_alpha(tlo))
    lo = [tlo]; hi = [thi]; fixed = set()
    for fi, k in enumerate(G.sizes):
        if k == 3: continue
        if k == 4:
            lo.append(a_lo); hi.append(up(2 * _alpha(thi))); continue
        i1, i2 = G.pidx[fi]
        lo += [a_lo, a_lo]; hi += [up(I._PI), up(I._PI)]
        mode = G.modes[fi]; f = G.faces[fi]
        # the parameter at boundary position 4: formula position 2 in region b, formula position 1 in region c
        p4 = i2 if mode == 'b' else (i1 if mode == 'c' else None)
        if p4 is not None:
            lo[p4] = am.LSTRAIGHT                        # the region threshold, the same double as the clamp
            if len(G.rot[f[4]]) == 2:                    # L4: a degree-2 vertex has two straight corners
                lo[p4] = down(I._PI); fixed.add(p4)
    return lo, hi, fixed


def required_modes(G):
    """V2: the pentagon regions that must be present (docs/lemmas.md L10', L11').  With c_i the corner at
    boundary position i and L = LSTRAIGHT: region a is c4 <= L; region b is c4 >= L and c1 <= L; region c is
    c4 >= L and c1 >= L (then c2 < L by Lemma P).  Which regions can be nonempty follows from the degrees d4,
    d1 of the vertices at positions 4 and 1: a corner >= L sits at a vertex of degree <= 3 (L11), and a
    degree-2 vertex has a straight corner (L4).  The returned list is a superset of the nonempty regions (for
    d4 = 3, d1 = 2 it lists region b although c1 = pi makes it empty); an extra chart adds work, never removes
    a case.  Returns the list of mode dicts {face index: mode}."""
    import itertools
    pents = [fi for fi, k in enumerate(G.sizes) if k == 5]
    opts = []
    for fi in pents:
        f = G.faces[fi]; d4 = len(G.rot[f[4]]); d1 = len(G.rot[f[1]])
        if d4 == 2:
            opts.append('c' if d1 == 2 else ('b' if d1 >= 4 else 'bc'))
        else:
            opts.append('a' if d4 >= 4 else ('ab' if d1 >= 4 else 'abc'))
    return [dict(zip(pents, combo)) for combo in itertools.product(*opts)]


def canonical_code(faces):
    """V4: canonical code of the rotation system of the oriented faces, written from the definition in the
    canon module's docstring (and independently of it): for every directed start edge (u, w) and both
    orientations, number the vertices breadth-first, visiting the neighbours of each vertex in rotation order
    from the edge by which it was first reached; the code lists, per vertex, its degree and its neighbours'
    numbers; the canonical code is the lexicographic minimum.  Purpose: tie the certificate to the candidate
    list of L12-L13, whose entries are these codes.  Minimizing over both orientations identifies a plane graph
    with its mirror image on purpose: the reflection of a critical configuration is critical, so one
    certificate covers both."""
    nv = 1 + max(v for f in faces for v in f)
    nxt = [dict() for _ in range(nv)]            # nxt[v][b] = a for a corner (a, v, b): a follows b at v
    for f in faces:
        k = len(f)
        for pos in range(k):
            nxt[f[pos]][f[(pos + 1) % k]] = f[pos - 1]
    rot = []
    for v in range(nv):
        start = min(nxt[v]); seq = [start]
        while True:
            w = nxt[v][seq[-1]]
            if w == start: break
            seq.append(w)
        rot.append(seq)
    best = None
    for u in range(nv):
        for w in rot[u]:
            for rev in (False, True):
                num = {u: 0}; queue = [u]; via = {u: w}; code = []; head = 0
                while head < len(queue):
                    v = queue[head]; head += 1
                    r = rot[v][::-1] if rev else rot[v]
                    i = r.index(via[v]); r = r[i:] + r[:i]
                    code.append(len(r))
                    for x in r:
                        if x not in num:
                            num[x] = len(queue); queue.append(x); via[x] = v
                        code.append(num[x])
                if len(queue) == nv and (best is None or code < best): best = code
    return '.'.join(map(str, best))


def faces_match(cert_faces, gen_faces):
    """V4(b): every gen face is a cyclic rotation (same orientation) of exactly one certificate face."""
    def key(f):
        """Cyclic normal form of a face: rotated to start at its smallest vertex (orientation kept)."""
        k = len(f); i = f.index(min(f)); return tuple(f[i:] + f[:i])
    return sorted(key(list(f)) for f in cert_faces) == sorted(key(list(f)) for f in gen_faces)


def structural_check(cert, code=None, gen_faces=None, r_min=R_MIN):
    """V1, V1b, V2-V4 on a whole certificate; `code` is the candidate's canonical code from the record that
    carries the certificate, `gen_faces` the faces of its gen record.  Returns a list of problems (empty when
    all pass).
    V1b (docs/lemmas.md L31-L33), on every chart's start box [lo0, hi0] (every node box lies inside it): all
    endpoints finite; a rhombus coordinate has 0 < lo0 and hi0 <= pi, because the model's range of the corner
    O from two corners (model.py, `ex`) is valid only for A in (0, pi] (L32); a pentagon coordinate has
    0 <= lo0 and hi0 <= 2 pi, hygiene only, no pentagon bound is needed for soundness (L31).  The bounds on
    hi0 are decided in Arb as "certainly <=" against pi and 2 pi enclosed at 128 bits (radius about 2^-126,
    whatever the global precision): so hi0 = 3.141592653589793 (the double below pi) is accepted and the next
    double up is not; lo0 is compared as the double it is.  The theta coordinate is covered by V3, the
    straight-corner coordinate by V1 (its box contains pi)."""
    from flint import ctx
    with ctx.workprec(128):
        PI = arb.pi(); TWO_PI = 2 * PI
    problems = []
    if cert.get('version') not in (1, 2): return ['unknown certificate version']
    if not theta_domain_ok(cert['theta']):
        problems.append('V3: theta range outside the model domain [%r, %r] (docs/lemmas.md L17)'
                        % (am.THETA_FLOOR, am.THETA_CEIL))
    elif not theta_claim_ok(cert['theta'], r_min): problems.append('V3: theta range does not cover the claim')
    try:
        own = canonical_code(cert['faces'])
        if code is not None and own != code: problems.append('V4: canonical code mismatch')
    except Exception as e:
        problems.append('V4: faces are not a rotation system (%r)' % (e,))
    if gen_faces is not None and not faces_match(cert['faces'], gen_faces):
        problems.append('V4: faces differ from the gen record')
    G0 = am.Graph(cert['faces'], {})
    present = {}
    for ch in cert['charts']:
        present[tuple(sorted((int(k), v) for k, v in ch['modes'].items()))] = ch
    for req in required_modes(G0):
        ch = present.get(tuple(sorted(req.items())))
        if ch is None or ch.get('not_run') or not ch.get('complete', False):
            problems.append('V2: required chart %s missing or incomplete' % (req,))
    for ch in cert['charts']:
        if ch.get('not_run'): continue
        G = am.Graph(cert['faces'], ch['modes'])
        if G.npar != len(ch['lo0']) or G.npar != len(ch['hi0']):
            problems.append('V1: parameter count mismatch'); continue
        lo, hi, fixed = apriori_box(G, cert['theta'])
        for j in range(G.npar):
            if not (float(ch['lo0'][j]) <= lo[j] and float(ch['hi0'][j]) >= hi[j]):
                problems.append('V1: chart %s coordinate %d box [%r, %r] does not contain the a-priori box [%r, %r]'
                                % (ch['modes'], j, ch['lo0'][j], ch['hi0'][j], lo[j], hi[j]))
        rhombus = {G.pidx[fi][0] for fi, k in enumerate(G.sizes) if k == 4}
        pentagon = {i for fi, k in enumerate(G.sizes) if k == 5 for i in G.pidx[fi]}
        for j in range(G.npar):                                             # V1b (L31-L33)
            a, b = float(ch['lo0'][j]), float(ch['hi0'][j])
            ok = math.isfinite(a) and math.isfinite(b)
            if ok and j in rhombus: ok = a > 0 and bool(arb(b) <= PI)
            elif ok and j in pentagon: ok = a >= 0 and bool(arb(b) <= TWO_PI)
            if not ok:
                problems.append('V1b: chart %s coordinate %d box [%r, %r] outside the model domain '
                                '(docs/lemmas.md L32)' % (ch['modes'], j, ch['lo0'][j], ch['hi0'][j]))
    return problems


def check_certificate(cert, code=None, gen_faces=None, r_min=R_MIN, structural=True):
    """Check one certificate.  Returns a dict: verified (bool), structural problems, per-reason counts,
    failed boxes, terminal boxes, number of open leaves and of charts not run.  `code` and `gen_faces` are
    the candidate's code and gen faces (V4); without them the code identity is not checked.
    Claim: verified is True only if every structural check passes, every removed slab and every excluded leaf
    of every required chart is proved to contain no critical configuration, and every terminal leaf is no
    wider than the stated tolerance.  Terminal leaves are NOT excluded: a certificate with terminal boxes
    proves completeness only once each terminal box is identified with a census configuration (checklist C4,
    the manifest's terminal_unidentified count)."""
    out = dict(verified=True, structural=[], failed=[], terminal=[], open=0, not_run=0, counts={}, leaves=0)
    if cert.get('version') not in (1, 2): raise ValueError('unknown certificate version')
    if structural:
        out['structural'] = structural_check(cert, code, gen_faces, r_min)
        if out['structural']: out['verified'] = False
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
                CB, _ = G.corners(lo, hi, tighten=False)           # only its gradients are used
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
                # a terminal box is not an exclusion: it is recorded for identification (checklist C4) and its
                # width must not exceed the stated tolerance, so that the search cannot park a large region here
                if max(b - a for a, b in zip(lo, hi)) > float(cert['tol']):
                    out['failed'].append(dict(kind='terminal-too-wide', lo=lo, hi=hi)); out['verified'] = False
                out['terminal'].append(dict(modes=ch['modes'], lo=lo, hi=hi)); continue
            if kind == 'open':
                out['open'] += 1; out['verified'] = False; continue
            if kind == 'linear' and extra:
                C_, _ = G.corners(lo, hi, tighten=False)
                if K.test_linear_mu(lo, hi, lo, hi, C_, extra):
                    K.bump('linear-hint'); continue
            if not K.exclude(lo, hi):
                out['failed'].append(dict(kind=kind, lo=lo, hi=hi)); out['verified'] = False
        for k, v in K.stats.items(): out['counts'][k] = out['counts'].get(k, 0) + v
    return out
