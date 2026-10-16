"""Exact second-order data at the two kissing configurations FCC and HCP (r = 1).

Status: trusted. Part of the trusted base (transfer prompt v1.2, code standards).
Implements: paper Appendix B, Proposition B.5 (exact second-order data at FCC and HCP), in exact arithmetic.

For a balanced configuration C at r = 1 with contact graph G, let Q_w = (cD - W) (x) I_3 restricted to the tangent
space (the stress Laplacian of L17, the Hessian of -sum_ij w_ij <u_i, u_j> in the exponential chart).  The
second-order isolation theorem (paper, Appendix B; Project doc research-r1-second-order-isolation-2026-10-15)
needs a lower bound lambda for v^T Q_w v on the critical cone, uniformly over normalized stresses w.  When G
carries a full-support stress the critical cone is the space of nontrivial infinitesimal flexes of G on S^2.

This script computes, in exact arithmetic (sympy over Q(sqrt 2) and Q(sqrt 2, sqrt 3)):
  - the contact graph, the rank of the spherical rigidity matrix, and the nontrivial flex space;
  - for each edge e, the value f^T Q_e f along the unit flex f, so that f^T Q_w f = sum_e w_e (f^T Q_e f) for
    every normalized w, balanced or not;
  - the eigenvalue of Q_w on the flex space for the hexagon stresses and the uniform stress.

Results (exact): FCC: rank 20, one nontrivial flex (the jitterbug), per-edge value 1/6 on all 24 edges, so
lambda = 1/6.  HCP: rank 20, one nontrivial flex, per-edge values {1/12, 8/33, 3/22}, so lambda >= 1/12.

Ported into the repository on 2026-10-16 from the Project doc code-r1-exact-lambda-2026-10-15.py, because the
v4.1 commits that first carried it were never pushed.  Status: exact.
"""
import itertools
import sympy as sp

HALF = sp.Rational(1, 2)


def contact_edges(U):
    """Pairs at exactly 60 degrees (inner product exactly 1/2, decided symbolically)."""
    N = len(U)
    return [(i, j) for i in range(N) for j in range(i + 1, N) if sp.simplify(U[i].dot(U[j]) - HALF) == 0]


def tangent_bases(U):
    """An exact orthonormal basis of the tangent plane at each point."""
    T = []
    for u in U:
        a = u.cross(sp.Matrix([1, 0, 0])) if u[0] ** 2 != 1 else u.cross(sp.Matrix([0, 1, 0]))
        a = a / sp.sqrt(a.dot(a))
        T.append(sp.Matrix.hstack(a, u.cross(a)))
    return T


def rigidity(U, E, T):
    """Rows d/dt <u_i, u_j> in tangent coordinates (2 per point)."""
    N = len(U)
    R = sp.zeros(len(E), 2 * N)
    for k, (i, j) in enumerate(E):
        R[k, 2 * i:2 * i + 2] = (T[i].T * U[j]).T
        R[k, 2 * j:2 * j + 2] = (T[j].T * U[i]).T
    return R


def rotations(U, T):
    """The three infinitesimal rotations in tangent coordinates."""
    N = len(U)
    rot = sp.zeros(2 * N, 3)
    for a in range(3):
        e = sp.zeros(3, 1); e[a] = 1
        for i in range(N):
            rot[2 * i:2 * i + 2, a] = T[i].T * e.cross(U[i])
    return rot


def stress_laplacian(U, E, w):
    """(cD - W) (x) I_3 for the edge weights w (list indexed like E)."""
    N = len(U)
    W = sp.zeros(N, N)
    for k, (i, j) in enumerate(E):
        W[i, j] = W[j, i] = w[k]
    c = [sum(W[i, j] * U[i].dot(U[j]) for j in range(N)) for i in range(N)]
    return sp.kronecker_product(sp.diag(*c) - W, sp.eye(3))


def analyze(U, named_stresses):
    """Return a dict with rank, flex dimension, per-edge values along the flex and eigenvalues per stress."""
    N = len(U)
    E = contact_edges(U)
    T = tangent_bases(U)
    R = rigidity(U, E, T)
    rot = rotations(U, T)
    F = sp.Matrix.vstack(R, rot.T).nullspace()          # flexes orthogonal to the rotations
    Tb = sp.zeros(3 * N, 2 * N)
    for i in range(N):
        Tb[3 * i:3 * i + 3, 2 * i:2 * i + 2] = T[i]
    out = dict(contacts=len(E), rank=R.rank(), flex_dim=len(F), per_edge=set(), eigen={})
    if len(F) == 1:
        f = F[0] / sp.sqrt(F[0].dot(F[0]))
        for k in range(len(E)):
            w = [0] * len(E); w[k] = 1
            out['per_edge'].add(sp.nsimplify(sp.simplify((f.T * Tb.T * stress_laplacian(U, E, w) * Tb * f)[0])))
    Fm = sp.Matrix.hstack(*F)
    G = Fm.T * Fm
    for lab, w in named_stresses(E, U):
        Qf = Fm.T * Tb.T * stress_laplacian(U, E, w) * Tb * Fm
        out['eigen'][lab] = {sp.nsimplify(k): v for k, v in sp.simplify(G.inv() * Qf).eigenvals().items()}
    return out


# ---------------------------------------------------------------- the two configurations

def fcc():
    """The cuboctahedron (FCC kissing configuration): the 12 permutations of (+-1, +-1, 0)/sqrt 2."""
    s2 = sp.sqrt(2)
    return [sp.Matrix(p) / s2 for p in itertools.product([-1, 0, 1], repeat=3) if sum(abs(x) for x in p) == 2]


def fcc_stresses(E, U):
    """The four great-hexagon stresses of FCC (each a balanced 6-cycle) and the uniform stress, normalized."""
    out = []
    for n in [(1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1)]:
        n = sp.Matrix(n)
        w = [1 if (U[i].dot(n) == 0 and U[j].dot(n) == 0) else 0 for (i, j) in E]
        assert sum(w) == 6
        out.append(('hexagon%s' % str(tuple(n)), [sp.Rational(x, 6) for x in w]))
    out.append(('uniform', [sp.Rational(1, len(E))] * len(E)))
    return out


def hcp():
    """The triangular orthobicupola (HCP kissing configuration): a great hexagon and two aligned triangles at latitude arctan(sqrt 2)."""
    s3 = sp.sqrt(3)
    hexa = [sp.Matrix([sp.cos(k * sp.pi / 3), sp.sin(k * sp.pi / 3), 0]) for k in range(6)]
    cphi, sphi = 1 / s3, sp.sqrt(sp.Rational(2, 3))        # latitude arctan(sqrt 2): distance 60 deg to the hexagon
    tri = lambda z: [sp.Matrix([cphi * sp.cos(t), cphi * sp.sin(t), z]) for t in [sp.pi / 6, 5 * sp.pi / 6, 3 * sp.pi / 2]]
    return [sp.simplify(x) for x in hexa + tri(sphi) + tri(-sphi)]


def hcp_stresses(E, U):
    """The equatorial-hexagon stress of HCP and the uniform stress, normalized."""
    w = [1 if (U[i][2] == 0 and U[j][2] == 0) else 0 for (i, j) in E]
    assert sum(w) == 6
    return [('equator-hexagon', [sp.Rational(x, 6) for x in w]), ('uniform', [sp.Rational(1, len(E))] * len(E))]


if __name__ == '__main__':
    for name, U, S in [('FCC', fcc(), fcc_stresses), ('HCP', hcp(), hcp_stresses)]:
        r = analyze(U, S)
        print('%s: contacts %d, rank %d, nontrivial flexes %d' % (name, r['contacts'], r['rank'], r['flex_dim']))
        print('  per-edge values f^T Q_e f along the flex:', sorted(r['per_edge']))
        for lab, ev in r['eigen'].items():
            print('  %-18s eigenvalues of Q_w on Flex/SO(3): %s' % (lab, ev))
