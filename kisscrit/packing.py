"""Equal-ball packings in a container as a min-type problem: the same toolkit as on the sphere.

Status: untrusted (exploration or tooling; no part of the proof).

Configurations: n centers X in R^d.  The packing radius
    rho(X) = min over active functions g_k(X)
with
    pair   g = |x_i - x_j - lam| / 2        (lam = 0, or a lattice vector for the torus)
    wall   g = x_ia  or  1 - x_ia           (unit cube [0,1]^d)
    sphere g = 1 - |x_i|                    (unit ball)
rho is a min of smooth functions, so its critical points are the configurations whose active set carries a
nonnegative balanced stress, sum_k sigma_k grad g_k = 0, and with a fixed active set A they are the
solutions of the square system

    g_k(X) = t (k in A),   sum_k sigma_k grad g_k(X) = 0,   sum sigma = 1     (+ gauge)

with n d + 1 + |A| unknowns and |A| + n d + 1 equations.  The continuous symmetries are gauged:
translations for the torus (sum_i x_i = sum_i x0_i), rotations for the ball (sum_i x0_i ^ x_i = 0).

Containers: Cube(d), Torus(basis), Ball(d).
"""
import itertools, json, math
import numpy as np
import networkx as nx
from scipy.optimize import linprog


# ----------------------------------------------------------------------------------------- containers
class Container:
    d = 2
    name = ''
    def constraints(self, X):
        """All candidate constraints (tuples) whose value could be near the minimum."""
        raise NotImplementedError
    def gauge_rows(self, X0):
        """(rows, rhs): linear gauge equations  R z_X = rhs  on the flattened coordinates."""
        return np.zeros((0, X0.size)), np.zeros(0)
    def gauge_motions(self, X):
        return np.zeros((0, X.size))
    def random(self, n, rng):
        raise NotImplementedError
    def point_group(self):
        """Representatives of Iso(container) / continuous gauge group, as (matrix, center) acting x -> M (x-c) + c."""
        raise NotImplementedError
    def wrap(self, X):
        return X


class Cube(Container):
    def __init__(self, d=2):
        self.d = d; self.name = 'cube%d' % d
    def constraints(self, X):
        n = len(X); out = [('pair', i, j, None) for i in range(n) for j in range(i + 1, n)]
        out += [('wall', i, a, s) for i in range(n) for a in range(self.d) for s in (0, 1)]
        return out
    def random(self, n, rng):
        return rng.uniform(0.05, 0.95, size=(n, self.d))
    def point_group(self):
        d = self.d; G = []
        for perm in itertools.permutations(range(d)):
            for signs in itertools.product((1, -1), repeat=d):
                M = np.zeros((d, d))
                for a, p in enumerate(perm): M[a, p] = signs[a]
                G.append(M)
        return G, np.full(d, 0.5)


class Torus(Container):
    def __init__(self, basis, exact=None):
        """basis: rows are lattice vectors (floats).  exact: optional function returning the basis rows as
        mpmath numbers at the current precision (needed for high-precision refinement of irrational lattices)."""
        self.L = np.array(basis, float); self.d = len(self.L); self.exact = exact
        self.name = 'torus%d' % self.d
        self.images = [np.array(k) @ self.L for k in itertools.product((-1, 0, 1), repeat=self.d)]
        self.half_short = 0.5 * min(np.linalg.norm(v) for v in self.images if np.linalg.norm(v) > 0)
    def wrap(self, X):
        F = X @ np.linalg.inv(self.L); F -= np.floor(F); return F @ self.L
    def constraints(self, X):
        n = len(X); out = []
        for i in range(n):
            for j in range(i + 1, n):
                for lam in self.images:
                    out.append(('pair', i, j, tuple(lam)))
        return out
    def gauge_rows(self, X0):
        n, d = X0.shape; R = np.zeros((d, n * d))
        for a in range(d):
            R[a, a::d] = 1
        return R, X0.sum(0)
    def gauge_motions(self, X):
        n, d = X.shape; M = np.zeros((d, n * d))
        for a in range(d): M[a, a::d] = 1
        return M
    def random(self, n, rng):
        return rng.uniform(0, 1, size=(n, self.d)) @ self.L
    def point_group(self):
        """Linear maps preserving the lattice (integer matrices in the basis that are orthogonal)."""
        d = self.d; G = []
        B = self.L.T                        # columns = basis vectors
        for K in itertools.product((-1, 0, 1), repeat=d * d):
            K = np.array(K, float).reshape(d, d)
            if abs(abs(np.linalg.det(K)) - 1) > 1e-9: continue
            M = B @ K @ np.linalg.inv(B)
            if np.allclose(M @ M.T, np.eye(d), atol=1e-9):
                G.append(M)
        return G, np.zeros(d)


class Ball(Container):
    def __init__(self, d=2):
        self.d = d; self.name = 'ball%d' % d
    def constraints(self, X):
        n = len(X)
        return [('pair', i, j, None) for i in range(n) for j in range(i + 1, n)] + [('sphere', i) for i in range(n)]
    def gauge_rows(self, X0):
        n, d = X0.shape; rows = []
        for a in range(d):
            for b in range(a + 1, d):
                r = np.zeros(n * d)
                for i in range(n):
                    r[i * d + b] += X0[i, a]; r[i * d + a] -= X0[i, b]
                rows.append(r)
        return np.array(rows).reshape(-1, n * d), np.zeros(len(rows))
    def gauge_motions(self, X):
        n, d = X.shape; rows = []
        for a in range(d):
            for b in range(a + 1, d):
                W = np.zeros((n, d)); W[:, a] = -X[:, b]; W[:, b] = X[:, a]; rows.append(W.ravel())
        return np.array(rows).reshape(-1, n * d)
    def random(self, n, rng):
        Y = rng.normal(size=(n, self.d)); r = rng.uniform(0, 1, size=(n, 1)) ** (1 / self.d)
        return 0.9 * r * Y / np.linalg.norm(Y, axis=1, keepdims=True)
    def point_group(self):
        return None, None                   # handled separately (O(d) symmetry group of X)


# ----------------------------------------------------------------------------------------- functions
def g_val(X, con):
    if con[0] == 'pair':
        _, i, j, lam = con; v = X[i] - X[j] - (np.array(lam) if lam is not None else 0)
        return 0.5 * np.linalg.norm(v)
    if con[0] == 'wall':
        _, i, a, s = con; return X[i, a] if s == 0 else 1 - X[i, a]
    _, i = con; return 1 - np.linalg.norm(X[i])


def g_grad(X, con):
    n, d = X.shape; g = np.zeros(n * d)
    if con[0] == 'pair':
        _, i, j, lam = con; v = X[i] - X[j] - (np.array(lam) if lam is not None else 0)
        u = v / (2 * np.linalg.norm(v)); g[i * d:(i + 1) * d] += u; g[j * d:(j + 1) * d] -= u
    elif con[0] == 'wall':
        _, i, a, s = con; g[i * d + a] = 1 if s == 0 else -1
    else:
        _, i = con; g[i * d:(i + 1) * d] = -X[i] / np.linalg.norm(X[i])
    return g


def g_hess_add(H, X, con, w):
    n, d = X.shape
    if con[0] == 'pair':
        _, i, j, lam = con; v = X[i] - X[j] - (np.array(lam) if lam is not None else 0)
        r = np.linalg.norm(v); u = v / r; B = w * (np.eye(d) - np.outer(u, u)) / (2 * r)
        si, sj = slice(i * d, (i + 1) * d), slice(j * d, (j + 1) * d)
        H[si, si] += B; H[sj, sj] += B; H[si, sj] -= B; H[sj, si] -= B
    elif con[0] == 'sphere':
        _, i = con; r = np.linalg.norm(X[i]); u = X[i] / r
        si = slice(i * d, (i + 1) * d); H[si, si] -= w * (np.eye(d) - np.outer(u, u)) / r


def rho(cont, X):
    return min(g_val(X, k) for k in cont.constraints(X))


def active(cont, X, tol=1e-9):
    cons = cont.constraints(X); vals = np.array([g_val(X, k) for k in cons]); m = vals.min()
    return [k for k, v in zip(cons, vals) if v <= m + tol], float(m)


# ----------------------------------------------------------------------------------------- ascent
def ascend(cont, X, maxit=3000, delta0=1e-3, h0=0.05):
    X = cont.wrap(np.array(X, float)); n, d = X.shape; delta = delta0
    for it in range(maxit):
        cons = cont.constraints(X); vals = np.array([g_val(X, k) for k in cons]); m = vals.min()
        act = [k for k, v in zip(cons, vals) if v <= m + delta]
        G = np.array([g_grad(X, k) for k in act]); kk, N = G.shape
        # maximize tau s.t. G v >= tau, |v|_inf <= 1  (gauge directions are harmless)
        cost = np.zeros(N + 1); cost[-1] = -1
        res = linprog(cost, A_ub=np.hstack([-G, np.ones((kk, 1))]), b_ub=np.zeros(kk),
                      bounds=[(-1, 1)] * N + [(None, 1)], method='highs')
        if res.x is None:
            tau = 0.0; v = np.zeros((n, d))
        else:
            tau = res.x[-1]; v = res.x[:-1].reshape(n, d)
        if tau < 1e-9:
            if delta < 1e-12: break
            delta *= 0.1; continue
        h = min(h0, delta / tau); ok = False
        while h > 1e-14:
            Y = cont.wrap(X + h * v)
            if rho(cont, Y) > m + 0.3 * h * tau:
                ok = True; break
            h *= 0.5
        if not ok:
            delta *= 0.1
            if delta < 1e-13: break
            continue
        X = Y; delta = min(delta0, max(2 * delta, 1e-6))
    return X, rho(cont, X)


# ----------------------------------------------------------------------------------------- KKT
def kkt_solve(cont, X, A, s0=None, iters=150, tol=1e-13):
    X = np.array(X, float); n, d = X.shape; k = len(A); N = n * d
    R, rhs = cont.gauge_rows(X)
    G = np.array([g_grad(X, c) for c in A])
    if s0 is None:
        M = np.vstack([G.T, np.ones(k)]); b = np.zeros(N + 1); b[-1] = 1
        s0, *_ = np.linalg.lstsq(M, b, rcond=None)
    z = np.concatenate([X.ravel(), [np.mean([g_val(X, c) for c in A])], s0])

    def FJ(z):
        Y = z[:N].reshape(n, d); t = z[N]; s = z[N + 1:]
        G = np.array([g_grad(Y, c) for c in A])
        F1 = np.array([g_val(Y, c) for c in A]) - t
        H = np.zeros((N, N))
        for c, w in zip(A, s): g_hess_add(H, Y, c, w)
        F2 = G.T @ s; F3 = [s.sum() - 1]; F4 = R @ z[:N] - rhs
        F = np.concatenate([F1, F2, F3, F4])
        J = np.zeros((len(F), len(z)))
        J[:k, :N] = G; J[:k, N] = -1
        J[k:k + N, :N] = H; J[k:k + N, N + 1:] = G.T
        J[k + N, N + 1:] = 1
        J[k + N + 1:, :N] = R
        return F, J

    lam = 1e-8; F, J = FJ(z); nf = np.linalg.norm(F)
    for _ in range(iters):
        if nf < tol: break
        JTJ = J.T @ J; g = J.T @ F; D = np.diag(np.diag(JTJ)) + 1e-12 * np.eye(len(z)); improved = False
        for _ in range(14):
            try:
                dz = np.linalg.solve(JTJ + lam * D, -g)
            except np.linalg.LinAlgError:
                lam *= 10; continue
            F2, J2 = FJ(z + dz); n2 = np.linalg.norm(F2)
            if np.isfinite(n2) and n2 < nf:
                z, F, J, nf = z + dz, F2, J2, n2; lam = max(lam / 10, 1e-15); improved = True; break
            lam *= 10
        if not improved: break
    return z[:N].reshape(n, d), float(z[N]), z[N + 1:], float(nf)


def kkt_check(cont, X, t, s, A, nf, res_tol=1e-10, gap_tol=1e-8, s_tol=1e-9):
    if nf > res_tol or t <= 1e-6 or s.min() < s_tol:
        return False
    if isinstance(cont, Torus) and t >= cont.half_short - 1e-9:
        return False
    key = set(map(_ckey, A))
    others = [g_val(X, c) for c in cont.constraints(X) if _ckey(c) not in key]
    return min(others) > t + gap_tol if others else True


def _ckey(c):
    return tuple(tuple(np.round(x, 6)) if isinstance(x, tuple) else x for x in c)


# ----------------------------------------------------------------------------------------- analysis
def analyze(cont, X, tol=1e-9, tol_ev=1e-7):
    X = np.array(X, float); n, d = X.shape; N = n * d
    A, t = active(cont, X, tol)
    G = np.array([g_grad(X, c) for c in A]); k = len(A)
    cost = np.zeros(k + 1); cost[-1] = -1
    A_eq = np.zeros((N + 1, k + 1)); A_eq[:N, :k] = G.T; A_eq[N, :k] = 1; b_eq = np.zeros(N + 1); b_eq[N] = 1
    A_ub = np.zeros((k, k + 1)); A_ub[:, :k] = -np.eye(k); A_ub[:, -1] = 1
    res = linprog(cost, A_ub=A_ub, b_ub=np.zeros(k), A_eq=A_eq, b_eq=b_eq,
                  bounds=[(0, None)] * k + [(None, None)], method='highs')
    if res.status != 0:
        return dict(critical=False, rho=t, k=k)
    sig = res.x[:k]
    support = sig > 1e-10
    if res.x[-1] <= 1e-10:   # find a stress of maximal support
        support = np.zeros(k, bool)
        for e in range(k):
            c2 = np.zeros(k); c2[e] = -1
            r2 = linprog(c2, A_eq=A_eq[:, :k], b_eq=b_eq, bounds=[(0, None)] * k, method='highs')
            support[e] = r2.status == 0 and -r2.fun > 1e-9
        idx = np.where(support)[0]
        if len(idx):
            Gs = G[idx]; kk = len(idx)
            c3 = np.zeros(kk + 1); c3[-1] = -1
            Ae = np.zeros((N + 1, kk + 1)); Ae[:N, :kk] = Gs.T; Ae[N, :kk] = 1
            Au = np.zeros((kk, kk + 1)); Au[:, :kk] = -np.eye(kk); Au[:, -1] = 1
            r3 = linprog(c3, A_ub=Au, b_ub=np.zeros(kk), A_eq=Ae, b_eq=b_eq,
                         bounds=[(0, None)] * kk + [(None, None)], method='highs')
            sig = np.zeros(k); sig[idx] = r3.x[:kk]
    Asup = [c for c, s in zip(A, support) if s]; ssup = sig[support]
    touched = set()
    for c in Asup: touched.add(c[1]); (touched.add(c[2]) if c[0] == 'pair' else None)
    ratt = [i for i in range(n) if i not in touched]
    free = [i * d + a for i in range(n) if i not in ratt for a in range(d)]
    Gs = np.array([g_grad(X, c) for c in Asup])[:, free]
    M = np.hstack([Gs, -np.ones((len(Asup), 1))])
    U, S, Vt = np.linalg.svd(M); rk = int((S > 1e-9 * max(S[0], 1)).sum())
    K = Vt[rk:].T[:len(free)]            # motions of non-rattlers keeping supported values equal
    Gm = cont.gauge_motions(X)
    if Gm.size:
        Gm = Gm[:, free]
        if ratt and isinstance(cont, Ball):
            pass
        Q, _ = np.linalg.qr(Gm.T); rG = np.linalg.matrix_rank(Gm, 1e-9); Q = Q[:, :rG]
        K = K - Q @ (Q.T @ K)
    if K.shape[1]:
        U2, S2, _ = np.linalg.svd(K, full_matrices=False); K = U2[:, S2 > 1e-8]
    H = np.zeros((N, N))
    for c, w in zip(Asup, ssup): g_hess_add(H, X, c, w)
    H = H[np.ix_(free, free)]
    if K.shape[1]:
        ev, V = np.linalg.eigh(K.T @ H @ K)
        F0 = K @ V[:, np.abs(ev) <= tol_ev]
        flat = np.zeros((N, F0.shape[1])); flat[free] = F0
    else:
        ev = np.array([]); flat = np.zeros((N, 0))
    return dict(flat=flat, critical=True, rho=t, k=k, nsupport=int(support.sum()), rattlers=ratt,
                kernel=int(K.shape[1]), ev=ev.tolist(), coindex=int((ev > tol_ev).sum()),
                nzero=int((np.abs(ev) <= tol_ev).sum()), active=[list(map(_jsonable, c)) for c in Asup],
                sigma=ssup.tolist())


def _jsonable(x):
    if isinstance(x, tuple): return list(x)
    if isinstance(x, (np.integer,)): return int(x)
    return x


# ----------------------------------------------------------------------------------------- symmetry / counts
def stabilizer_order(cont, X, tol=1e-6):
    """|{(pi, g) : g in Iso/G0, g X = pi X mod G0}|."""
    X = np.array(X, float); n, d = X.shape
    if isinstance(cont, Ball):
        from .symmetry import symmetries
        rk = np.linalg.matrix_rank(X, 1e-8)
        if rk == d:
            S = symmetries(X, tol)
        elif rk == 0:
            return None                      # all at the center: infinite
        else:
            U, Sv, Vt = np.linalg.svd(X)
            Y = X @ Vt[:rk].T
            if rk == d - 1:
                if rk == 1:
                    # collinear in the plane: symmetries are reflections/rotation by pi preserving the set
                    vals = Y[:, 0]; cnt = 0
                    for sgn in (1, -1):
                        if sorted(np.round(sgn * vals, 6)) == sorted(np.round(vals, 6)): cnt += 1
                    return 2 * cnt
                S = symmetries(Y, tol); return 2 * len(S)
            return None
        return len(S)                        # |G_X| ; count = 2 n!/|G_X|
    G, ctr = cont.point_group(); cnt = 0
    for M in G:
        Y = (X - ctr) @ M.T + ctr
        if isinstance(cont, Torus):
            Linv = np.linalg.inv(cont.L)
            for j in range(n):
                tau = X[j] - Y[0]; Z = Y + tau
                D = (Z[:, None, :] - X[None, :, :]) @ Linv; D -= np.round(D)
                dist = np.linalg.norm(D @ cont.L, axis=2)
                p = dist.argmin(1)
                if np.all(dist[np.arange(n), p] < tol) and len(set(p.tolist())) == n:
                    cnt += 1
        else:
            dist = np.linalg.norm(Y[:, None, :] - X[None, :, :], axis=2); p = dist.argmin(1)
            if np.all(dist[np.arange(n), p] < tol) and len(set(p.tolist())) == n:
                cnt += 1
    return cnt


def family_count(cont, X, flat, tol=1e-6):
    """Number of labeled components of a Morse-Bott family of critical points through X whose tangent
    space is span(flat), assuming each component is the convex set X + span(flat) ∩ (critical set):
    n! |Iso/G0| / |setwise stabilizer|, the stabilizer being the (pi, g) with g X - pi X in span(flat)
    (mod the gauge group and the lattice).  Each component is counted with chi = 1 (contractible)."""
    X = np.array(X, float); n, d = X.shape
    G, ctr = cont.point_group()
    B = flat.copy()
    Gm = cont.gauge_motions(X)
    if Gm.size: B = np.hstack([B, Gm.T])
    Q = np.linalg.qr(B)[0][:, :np.linalg.matrix_rank(B, 1e-9)] if B.shape[1] else np.zeros((n * d, 0))
    cnt = 0
    for M in G:
        Y = (X - ctr) @ M.T + ctr
        for perm in itertools.permutations(range(n)):
            Z = Y[list(perm)]
            D = (Z - X)
            if isinstance(cont, Torus):
                Linv = np.linalg.inv(cont.L); F = D @ Linv; F -= np.round(F); D = F @ cont.L
            v = D.ravel(); res = v - Q @ (Q.T @ v)
            if np.linalg.norm(res) < tol:
                cnt += 1
    return math.factorial(n) * len(G) // cnt if cnt else None


def labeled_count(cont, X):
    n = len(X); st = stabilizer_order(cont, X)
    if st is None or st == 0: return None
    if isinstance(cont, Ball):
        return 2 * math.factorial(n) // st
    G, _ = cont.point_group()
    return math.factorial(n) * len(G) // st


def euler_target(cont, n):
    """chi of the configuration space modulo the gauge group (labeled points)."""
    if isinstance(cont, Cube):
        if cont.d == 2: return 1 if n == 1 else 0           # chi Conf_n(R^2)
        return math.factorial(n) if cont.d == 3 else None   # chi Conf_n(R^3) = n!
    if isinstance(cont, Torus):
        if cont.d == 2:                                     # Conf_n(T^2)/T^2 = Conf_{n-1}(T^2 - pt)
            return int((-1) ** (n - 1) * math.factorial(n - 1))
        if cont.d == 3:                                     # Conf_{n-1}(T^3 - pt), chi(T^3 - pt) = 1
            return math.factorial(n - 1)
        return None
    if isinstance(cont, Ball):
        if cont.d == 2:                                     # Conf_n(C)/S^1 ~ M_{0,n+1}
            p = 1
            for k in range(2, n): p *= (1 - k)
            return p
        return None
    return None


# ----------------------------------------------------------------------------------------- census
class PCensus:
    def __init__(self, cont, n):
        self.cont, self.n = cont, n; self.types = {}

    def key(self, rec, digits=8):
        G = nx.Graph(); G.add_nodes_from(range(self.n), kind='pt')
        for c in rec['active']:
            if c[0] == 'pair':
                G.add_edge(c[1], c[2]) if not G.has_edge(c[1], c[2]) else G.edges[c[1], c[2]].__setitem__('m', G.edges[c[1], c[2]].get('m', 1) + 1)
            else:
                w = ('w', c[1], len(G)); G.add_node(w, kind='wall'); G.add_edge(c[1], w)
        h = nx.weisfeiler_lehman_graph_hash(G, node_attr='kind', edge_attr=None, iterations=3)
        return json.dumps([round(rec['rho'], digits), rec['coindex'], len(rec['rattlers']), h])

    def add(self, X, src=''):
        a = analyze(self.cont, X)
        if not a.get('critical'): return False
        flat = a.pop('flat')
        a['X'] = np.asarray(X).tolist(); a['src'] = src
        k = self.key(a)
        if k in self.types:
            self.types[k]['hits'] += 1; return False
        a['hits'] = 1
        a['labeled'] = labeled_count(self.cont, X) if not a['rattlers'] else None
        if a['nzero'] and not a['rattlers'] and not isinstance(self.cont, Ball):
            a['family_count'] = family_count(self.cont, X, flat)
        self.types[k] = a; return True

    def sorted(self):
        return sorted(self.types.values(), key=lambda t: (-t['rho'], t['coindex']))

    def signed_count(self):
        tot = 0; skipped = []
        for t in self.sorted():
            if t['nzero'] and t.get('family_count') and not t['rattlers']:
                tot += (-1) ** t['coindex'] * t['family_count']; continue
            if t['rattlers'] or t['nzero'] or t['labeled'] is None:
                skipped.append((round(t['rho'], 8), t['coindex'], len(t['rattlers']), t['nzero'])); continue
            tot += (-1) ** t['coindex'] * t['labeled']
        return tot, skipped

    def dump(self, path):
        json.dump(dict(container=self.cont.name, n=self.n, types=self.sorted()), open(path, 'w'))


def explore(cont, n, rng, ascents=100, kkt_starts=2000, verbose=False):
    """Maxima by random ascent; other critical points by solving the square system for active sets
    sampled near random and near-critical configurations."""
    C = PCensus(cont, n)
    seeds = []
    for _ in range(ascents):
        X, r = ascend(cont, cont.random(n, rng))
        A, t = active(cont, X, 1e-7)
        Xs, t2, s, nf = kkt_solve(cont, X, A)
        if kkt_check(cont, Xs, t2, s, A, nf):
            C.add(Xs, 'ascent'); seeds.append(Xs)
    for it in range(kkt_starts):
        X = cont.random(n, rng)
        if seeds and rng.random() < 0.5:
            X = cont.wrap(seeds[rng.integers(len(seeds))] + rng.normal(scale=10 ** rng.uniform(-2.5, -0.7), size=X.shape))
        X, _ = ascend(cont, X, maxit=int(rng.integers(0, 40)))
        cons = cont.constraints(X); vals = np.array([g_val(X, k) for k in cons]); order = np.argsort(vals)
        m = int(rng.integers(2, min(len(cons), n * cont.d + 1) + 1))
        A = [cons[q] for q in order[:m]]
        Xs, t, s, nf = kkt_solve(cont, X, A)
        if kkt_check(cont, Xs, t, s, A, nf):
            if C.add(cont.wrap(Xs), 'kkt') and verbose:
                print('new', round(t, 8), flush=True)
    return C


def _to_con(c):
    c = list(c)
    if c[0] == 'pair':
        return ('pair', int(c[1]), int(c[2]), tuple(c[3]) if c[3] is not None else None)
    if c[0] == 'wall':
        return ('wall', int(c[1]), int(c[2]), int(c[3]))
    return ('sphere', int(c[1]))


def surgery(C, rounds=4, near=10, pairs_out=True, verbose=False):
    """Contact-set surgery for a PCensus: delete one or two active constraints, add one of the `near`
    nearest inactive ones, or swap one for one; solve the square system from the type's X."""
    cont = C.cont; queue = list(C.types.values()); seen = set()
    for _ in range(rounds):
        new = []
        for t in queue:
            X = np.array(t['X']); A = [_to_con(c) for c in t['active']]
            cons = cont.constraints(X); key = set(map(_ckey, A))
            inact = sorted([c for c in cons if _ckey(c) not in key], key=lambda c: g_val(X, c))[:near]
            cands = [[a for a in A if a is not b] for b in A]
            if pairs_out:
                cands += [[a for a in A if a is not b and a is not b2] for b, b2 in itertools.combinations(A, 2)]
            for f in inact:
                cands.append(A + [f])
                cands += [[a for a in A if a is not b] + [f] for b in A]
            for B in cands:
                if len(B) < 2: continue
                sig = (json.dumps(sorted(map(str, map(_ckey, B)))), round(t['rho'], 7))
                if sig in seen: continue
                seen.add(sig)
                Xs, tt, s, nf = kkt_solve(cont, X, B)
                if kkt_check(cont, Xs, tt, s, B, nf):
                    if C.add(cont.wrap(Xs), 'surgery'):
                        new.append(list(C.types.values())[-1])
                        if verbose: print('  surgery: new rho=%.8f' % tt, flush=True)
        if not new: break
        queue = new
    return C


# ----------------------------------------------------------------------------------------- rattler strata
def _free_mask(cont, X, r0, P):
    """Mask of grid points P where one more ball center fits strictly: all its constraint values > r0."""
    n, d = X.shape; ok = np.ones(len(P), bool)
    if isinstance(cont, Cube):
        ok &= (P > r0).all(1) & (P < 1 - r0).all(1)
    elif isinstance(cont, Ball):
        ok &= np.linalg.norm(P, axis=1) < 1 - r0
    for i in range(n):
        if isinstance(cont, Torus):
            Linv = np.linalg.inv(cont.L); D = (P - X[i]) @ Linv; D -= np.round(D); D = D @ cont.L
            dist = np.min([np.linalg.norm(D - lam, axis=1) for lam in cont.images], axis=0)   # true minimum image
        else:
            dist = np.linalg.norm(P - X[i], axis=1)
        ok &= dist > 2 * r0
    return ok


def _grid_axes(cont, h):
    """Regular grid in fractional coordinates; returns (points, shape, periodic)."""
    d = cont.d
    if isinstance(cont, Ball):
        m = int(round(2 / h)); ax = -1 + (np.arange(m) + 0.5) * (2 / m)
        P = np.array(list(itertools.product(ax, repeat=d))); return P, (m,) * d, False
    m = int(round(1 / h)); ax = (np.arange(m) + 0.5) / m
    F = np.array(list(itertools.product(ax, repeat=d)))
    if isinstance(cont, Torus):
        return F @ cont.L, (m,) * d, True
    return F, (m,) * d, False


def _family_samples(cont, X, flat, r0, ns=30, smax=0.6):
    """Points X + s v of a <=1-dim family (v = flat direction) that stay critical at the same level."""
    if flat.shape[1] == 0:
        return [X]
    v = flat[:, 0].reshape(X.shape); v = v / np.abs(v).max()
    out = []
    for s in np.linspace(-smax, smax, 2 * ns + 1):
        Y = cont.wrap(X + s * v)
        if abs(rho(cont, Y) - r0) < 1e-9 and len(active(cont, Y, 1e-9)[0]) >= 1:
            out.append(Y)
    return out


def _cubical_euler(O, periodic_axes):
    """Euler characteristic of the union of closed voxels {O == True} (any dimension), with the given
    axes periodic.  Counts every cell (vertex, edge, ...) touched by an occupied voxel."""
    D = O.ndim; chi = 0
    for kind in itertools.product((0, 1), repeat=D):   # 1 = the cell spans this axis (voxel extent), 0 = lower boundary
        # a cell of this kind indexed by i is touched by voxels i - e for e in {0,1}^{axes with kind 0}
        C = np.zeros_like(O)
        zero_axes = [a for a in range(D) if kind[a] == 0]
        for shift in itertools.product((0, 1), repeat=len(zero_axes)):
            S = O
            for a, sh in zip(zero_axes, shift):
                if sh == 0: continue
                if periodic_axes[a]:
                    S = np.roll(S, 1, axis=a)
                else:
                    S = np.concatenate([np.zeros_like(S.take([0], axis=a)), S.take(range(S.shape[a] - 1), axis=a)], axis=a)
            C = C | S
        ncell = int(C.sum())
        # non-periodic axes have one more boundary cell at the top end
        for a in zero_axes:
            if not periodic_axes[a]:
                top = O.take([O.shape[a] - 1], axis=a)
                others = [b for b in zero_axes if b != a]
                T = top
                for b in others:
                    if periodic_axes[b]:
                        T = T | np.roll(T, 1, axis=b)
                    else:
                        T = T | np.concatenate([np.zeros_like(T.take([0], axis=b)), T.take(range(T.shape[b] - 1), axis=b)], axis=b)
                ncell += int(T.sum())
        dim = sum(kind)
        chi += (-1) ** dim * ncell
    return chi


def stratum_euler(cont, X, flat, r0, h=None):
    """Euler characteristic of {(family point, rattler position)} for one extra ball, on a grid.
    Families of dimension <= 1.  Returns None if unsupported."""
    if flat.shape[1] > 1:
        return None
    if h is None:
        h = 0.01 if cont.d == 2 else 0.04
    P, shape, per = _grid_axes(cont, h)
    fam = _family_samples(cont, X, flat, r0)
    layers = [_free_mask(cont, Y, r0, P).reshape(shape) for Y in fam]
    if len(layers) == 1:
        O = layers[0]; periodic = [per] * cont.d
    else:
        O = np.stack(layers, axis=0); periodic = [False] + [per] * cont.d
    if not O.any():
        return 0
    return _cubical_euler(O, periodic)


def stratum_components(cont, X, flat, r0, m_extra=1, h=None):
    """Number of connected components of the one-rattler stratum (kept for diagnostics)."""
    from scipy import ndimage
    if flat.shape[1] > 1: return None
    if h is None: h = 0.01 if cont.d == 2 else 0.04
    P, shape, per = _grid_axes(cont, h)
    fam = _family_samples(cont, X, flat, r0)
    O = np.stack([_free_mask(cont, Y, r0, P).reshape(shape) for Y in fam], axis=0)
    lab, nc = ndimage.label(O)
    return nc


def euler_with_rattlers(censuses, n):
    """Signed count for n balls: nondegenerate types and Morse-Bott families (counted as contractible) on n
    balls, plus one-rattler strata over types and <=1-dim families on n-1 balls, each counted by the Euler
    characteristic of {(family point, rattler position)} computed on a grid.  censuses: {m: PCensus}."""
    cont = censuses[n].cont
    tot, skipped = 0, []
    for t in censuses[n].sorted():
        if t['rattlers']:
            continue                                        # counted from m = n-1
        if t['nzero']:
            if t.get('family_count') is not None:
                tot += (-1) ** t['coindex'] * t['family_count']
            else:
                skipped.append(('n', round(t['rho'], 8), 'nzero', t['nzero']))
            continue
        if t['labeled'] is None:
            skipped.append(('n', round(t['rho'], 8), 'labeled?')); continue
        tot += (-1) ** t['coindex'] * t['labeled']
    if n - 1 in censuses:
        for t in censuses[n - 1].sorted():
            if t['rattlers']:
                continue
            X = np.array(t['X']); a = analyze(cont, X)
            flat = a.get('flat', np.zeros((X.size, 0)))
            nc = stratum_euler(cont, X, flat, t['rho'])
            if nc is None:
                skipped.append(('n-1', round(t['rho'], 8), 'family dim > 1')); continue
            if nc == 0:
                continue
            base = t.get('family_count') if t['nzero'] else t['labeled']
            if base is None:
                skipped.append(('n-1', round(t['rho'], 8), 'base count?')); continue
            tot += (-1) ** t['coindex'] * base * n * nc
    return tot, skipped


# ----------------------------------------------------------------------------------------- high precision
def refine_mp(cont, X, A, s0, dps=100):
    """Newton on the square system of active set A at `dps` digits (mpmath). Returns (X, t, s, residual)."""
    import mpmath as mp
    mp.mp.dps = dps
    X = np.array(X, float); n, d = X.shape; N = n * d; k = len(A)
    R, rhs = cont.gauge_rows(X)
    z = [mp.mpf(v) for v in X.ravel()] + [mp.mpf(float(np.mean([g_val(X, c) for c in A])))] + [mp.mpf(float(v)) for v in s0]
    def parts(z):
        Y = [[z[i * d + a] for a in range(d)] for i in range(n)]
        vals, grads, hess = [], [], []
        for c in A:
            gr = [mp.mpf(0)] * N; H = {}
            if c[0] == 'pair':
                _, i, j, lam = c
                if lam is None:
                    lamx = [mp.mpf(0)] * d
                elif isinstance(cont, Torus) and cont.exact is not None:
                    kk = np.round(np.array(lam) @ np.linalg.inv(cont.L)).astype(int); Lx = cont.exact()
                    lamx = [mp.fsum(int(kk[b]) * Lx[b][a] for b in range(d)) for a in range(d)]
                else:
                    lamx = [mp.mpf(x) for x in lam]
                v = [Y[i][a] - Y[j][a] - lamx[a] for a in range(d)]; r = mp.sqrt(mp.fsum(x * x for x in v))
                vals.append(r / 2)
                for a in range(d):
                    gr[i * d + a] += v[a] / (2 * r); gr[j * d + a] -= v[a] / (2 * r)
                for a in range(d):
                    for b in range(d):
                        h = ((1 if a == b else 0) - v[a] * v[b] / r ** 2) / (2 * r)
                        for (p, q, sg) in ((i, i, 1), (j, j, 1), (i, j, -1), (j, i, -1)):
                            H[(p * d + a, q * d + b)] = H.get((p * d + a, q * d + b), 0) + sg * h
            elif c[0] == 'wall':
                _, i, a, sd = c; vals.append(Y[i][a] if sd == 0 else 1 - Y[i][a]); gr[i * d + a] = 1 if sd == 0 else -1
            else:
                _, i = c; r = mp.sqrt(mp.fsum(x * x for x in Y[i])); vals.append(1 - r)
                for a in range(d): gr[i * d + a] = -Y[i][a] / r
                for a in range(d):
                    for b in range(d):
                        H[(i * d + a, i * d + b)] = -((1 if a == b else 0) - Y[i][a] * Y[i][b] / r ** 2) / r
            grads.append(gr); hess.append(H)
        return vals, grads, hess
    for it in range(40):
        vals, grads, hess = parts(z); t = z[N]; s = z[N + 1:]
        m = k + N + 1 + R.shape[0]; F = mp.matrix(m, 1); J = mp.matrix(m, len(z))
        for e in range(k):
            F[e] = vals[e] - t
            for q in range(N): J[e, q] = grads[e][q]
            J[e, N] = -1
        for q in range(N):
            F[k + q] = mp.fsum(s[e] * grads[e][q] for e in range(k))
            for e in range(k): J[k + q, N + 1 + e] = grads[e][q]
        for e in range(k):
            for (p, q), h in hess[e].items(): J[k + p, q] += s[e] * h
        F[k + N] = mp.fsum(s) - 1
        for e in range(k): J[k + N, N + 1 + e] = 1
        for g in range(R.shape[0]):
            F[k + N + 1 + g] = mp.fsum(mp.mpf(R[g, q]) * z[q] for q in range(N)) - mp.mpf(rhs[g])
            for q in range(N): J[k + N + 1 + g, q] = mp.mpf(R[g, q])
        JTJ = J.T * J
        for q in range(len(z)): JTJ[q, q] += mp.mpf(10) ** (-dps)      # regularize stress non-uniqueness
        dz = mp.lu_solve(JTJ, -(J.T * F))
        z = [z[q] + dz[q] for q in range(len(z))]
        if mp.norm(F) < mp.mpf(10) ** (-dps + 10): break
    vals, grads, hess = parts(z)
    res = max(abs(v - z[N]) for v in vals)
    return [[z[i * d + a] for a in range(d)] for i in range(n)], z[N], z[N + 1:], res
