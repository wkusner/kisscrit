"""Checks that the computer-assisted lemmas are not vacuous, and that Lemma P runs.

Lemma R uses the polygon vertices built by polygon_vertices and the inside tests det(P_m, P_m+1, c) >= 0.  If the
construction had the wrong orientation, the tests could pass for the wrong reason.  So: build sample polygons with
the same code, check that the sides have length theta, the polygons close up, the centroid passes every inside
test, and the deep-hole test fails just below the true covering radius.
"""
import math, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, 'scripts', 'completeness'))
import lemmas as L
from kisscrit.core.interval import Iv
from kisscrit.core import model as M

TH = math.radians(61.0)


def mid(v): return [0.5 * (c.lo + c.hi) for c in v]


def check_polygon(corners):
    P = L.polygon_vertices(Iv(TH), corners)
    Pm = [mid(p) for p in P]; k = len(Pm)
    for i in range(k - 1):                                   # the constructed sides
        d = math.acos(max(-1.0, min(1.0, sum(a * b for a, b in zip(Pm[i], Pm[i + 1])))))
        assert abs(d - TH) < 1e-9
    g = [sum(p[j] for p in Pm) for j in range(3)]; n = math.sqrt(sum(x * x for x in g)); g = [x / n for x in g]
    x = tuple(Iv(c) for c in g)
    assert all(L.det(P[i], P[(i + 1) % k], x).lo > 0 for i in range(k))
    phi = math.atan2(g[1], g[0]) % (2 * math.pi)
    assert 0 < phi < corners[0].lo
    return Pm


def test_triangle():
    al = M.alpha_range(TH, TH)
    Pm = check_polygon([al, al, al])
    assert abs(math.acos(sum(a * b for a, b in zip(Pm[2], Pm[0]))) - TH) < 1e-9     # closes up


def test_rhombus():
    C, ok = M.Graph([[0, 1, 2, 3], [3, 2, 1, 0]], {}).corners([TH, 1.6, 1.6], [TH, 1.6, 1.6])
    assert ok
    Pm = check_polygon([c[0] for c in C[0]])
    assert abs(math.acos(sum(a * b for a, b in zip(Pm[3], Pm[0]))) - TH) < 1e-9


def test_pentagons():
    for f1, f2 in [(2.2, 1.7), (2.0, 2.0), (1.7, 2.3)]:
        cs, ok = L.pentagon_values([(TH, TH), (f1, f1), (f2, f2)])
        assert ok
        assert sum(mid(cs)) > 3 * math.pi                        # positive spherical excess
        Pm = check_polygon(cs)
        assert abs(math.acos(sum(a * b for a, b in zip(Pm[4], Pm[0]))) - TH) < 1e-9


def test_lemma_P():
    L.lemma_P()                                              # raises SystemExit if not proved


def test_deep_hole_test_is_sharp():
    """no_deep_hole must fail when the threshold is below the true covering radius, and pass just above it."""
    for corners in ([M.alpha_range(TH, TH)] * 3,
                    L.pentagon_values([(TH, TH), (2.2, 2.2), (1.7, 1.7)])[0]):
        P = L.polygon_vertices(Iv(TH), corners)
        Pm = [mid(p) for p in P]
        # covering radius: max over a fine grid of the polygon of the distance to the nearest vertex
        best = 0.0; k = len(Pm)
        for i in range(1, k - 1):                     # fan triangles P0, P_i, P_{i+1}
            for u in range(41):
                for v in range(41 - u):
                    w = [(u * Pm[0][j] + v * Pm[i][j] + (40 - u - v) * Pm[i + 1][j]) / 40 for j in range(3)]
                    nw = math.sqrt(sum(t * t for t in w)); w = [t / nw for t in w]
                    best = max(best, min(math.acos(min(1.0, sum(a * b for a, b in zip(w, p)))) for p in Pm))
        assert best < TH
        assert not L.no_deep_hole(Iv(best - 1e-3), P)
        assert L.no_deep_hole(Iv(best + 2e-2), P)


def test_pentagon_vertex_chart():
    """pentagon_by_vertices builds an equilateral convex pentagon with the requested corners at P0 and P1, and
    the mean-value form (pentagon_quantities) agrees with it at a point."""
    from kisscrit.core import interval as I
    cosL = I.cos(Iv(L.LSTRAIGHT))
    for c0, f1 in [(2.3, 2.3), (2.0, 2.6), (2.5, 1.9), (2.2, 2.0)]:     # convex samples
        P, ok = L.pentagon_by_vertices(Iv(TH), Iv(c0), Iv(f1), cosL)
        assert ok
        Pm = [mid(p) for p in P]
        for i in range(5):                                   # sides theta, all turns left (convex)
            assert abs(math.acos(sum(a * b for a, b in zip(Pm[i], Pm[(i + 1) % 5]))) - TH) < 1e-9
            assert L.det(P[i - 1], P[i], P[(i + 1) % 5]).lo > 0
        # corners at P0 and P1
        def corner(i):
            A, B, C = Pm[i - 1], Pm[i], Pm[(i + 1) % 5]
            t1 = [a - sum(a * b for a, b in zip(A, B)) * b for a, b in zip(A, B)]
            t2 = [c - sum(c * b for c, b in zip(C, B)) * b for c, b in zip(C, B)]
            return math.acos(sum(x * y for x, y in zip(t1, t2)) / math.sqrt(sum(x * x for x in t1) * sum(x * x for x in t2)))
        assert abs(corner(0) - c0) < 1e-9 and abs(corner(1) - f1) < 1e-9
        q = L.pentagon_quantities(*[M.var(x, x, k) for k, x in enumerate((TH, c0, f1))], cosL)
        for i in range(5):
            t = L.det(P[i - 1], P[i], P[(i + 1) % 5])
            assert abs(q['turn'][i][0].lo - t.lo) < 1e-12
