"""Searching for critical configurations: random ascent, KKT polishing, and contact-graph surgery.

A *type* is a critical configuration up to O(d) and relabeling.  Types are keyed by the radius
(rounded), the degree sequence and a Weisfeiler-Lehman hash of the stressed contact graph; the key
is a heuristic invariant, and collisions are resolved by comparing sorted Gram spectra.
"""
import itertools, json, math
import numpy as np
import networkx as nx
from .geometry import normalize, contacts, random_config, r_of_c
from .ascent import ascend
from .kkt import solve_kkt, check_solution
from .stress import analyze
from .symmetry import group_orders


def graph_key(N, pairs):
    G = nx.Graph(); G.add_nodes_from(range(N)); G.add_edges_from(map(tuple, pairs))
    return tuple(sorted(d for _, d in G.degree())), nx.weisfeiler_lehman_graph_hash(G, iterations=3)


def type_key(X, pairs, c, digits=9):
    return (round(float(r_of_c(c)), digits),) + graph_key(len(X), pairs)


def polish(X, tol_list=(1e-9, 1e-7, 1e-5, 1e-4)):
    """Turn an (approximately) critical X into an exact-to-double critical configuration.

    Tries contact graphs at several tolerances; returns (X, info) or (None, None)."""
    for tol in tol_list:
        A, _ = contacts(X, tol)
        if len(A) < 2:
            continue
        Xs, c, s, nf = solve_kkt(X, A)
        ch = check_solution(Xs, c, s, A, nf)
        if ch:
            return Xs, ch
    return None, None


def record(Xs, extra=None, with_symmetry=True):
    a = analyze(Xs)
    if not a.get('critical'):
        return None
    rec = dict(r=float(r_of_c(a['c'])), c=a['c'], k=a['k'], nsupport=a['nsupport'], sdim=a['sdim'],
               tmin=a['tmin'], rattlers=a['rattlers'], kernel=a['kernel'], coindex=a['coindex'],
               nzero=a['nzero'], ev=a['ev'], pairs=a['pairs'], X=np.asarray(Xs).tolist())
    if with_symmetry and a['rattlers'] == 0:
        try:
            g, gp = group_orders(Xs)
            # labeled copies in BConf(N) = Conf(N)/SO(d): N!/|G ∩ SO(d)| (mirror images included)
            rec.update(sym=g, sym_proper=gp, labeled=math.factorial(len(Xs)) // gp)
        except ValueError:
            pass
    if extra:
        rec.update(extra)
    return rec


class Census:
    """Collect critical types; JSON-serializable."""

    def __init__(self, N, d):
        self.N, self.d = N, d; self.types = {}

    def add(self, Xs, src=''):
        rec = record(Xs, dict(src=src))
        if rec is None:
            return None, False
        key = json.dumps(type_key(Xs, rec['pairs'], rec['c']))
        if key in self.types:
            self.types[key]['hits'] += 1
            return key, False
        rec['hits'] = 1
        self.types[key] = rec
        return key, True

    def sorted(self):
        return sorted(self.types.values(), key=lambda t: -t['r'])

    def dump(self, path):
        json.dump(dict(N=self.N, d=self.d, types=self.sorted()), open(path, 'w'))

    @classmethod
    def load(cls, path):
        D = json.load(open(path)); C = cls(D['N'], D['d'])
        for t in D['types']:
            C.types[json.dumps(type_key(np.array(t['X']), t['pairs'], t['c']))] = t
        return C


def random_ascent(census, runs, rng, maxit=5000):
    N, d = census.N, census.d
    for _ in range(runs):
        X, info = ascend(random_config(N, d, rng), maxit=maxit)
        Xs, ch = polish(X)
        if Xs is not None:
            census.add(Xs, 'ascent')


def surgery(census, max_rounds=3, near=6, verbose=False):
    """Contact-graph surgery: from each known type, delete one contact, add one of the `near` nearest
    non-contacts, or swap one for one; solve the KKT system for the new graph and keep solutions."""
    N = census.N; seen_graphs = set(); queue = [t for t in census.sorted()]
    for rnd in range(max_rounds):
        new = []
        for t in queue:
            X = np.array(t['X']); A = [tuple(p) for p in t['pairs']]
            C = X @ X.T
            non = sorted(((C[i, j], (i, j)) for i in range(N) for j in range(i + 1, N) if (i, j) not in set(A)),
                         reverse=True)[:near]
            cands = []
            for e in A:
                cands.append([f for f in A if f != e])
            for _, f in non:
                cands.append(A + [f])
                for e in A:
                    cands.append([g for g in A if g != e] + [f])
            for B in cands:
                if len(B) < 2:
                    continue
                key = (json.dumps(graph_key(N, B)), round(t['r'], 6))
                if key in seen_graphs:
                    continue
                seen_graphs.add(key)
                Xs, c, s, nf = solve_kkt(X, B)
                if check_solution(Xs, c, s, B, nf):
                    k, is_new = census.add(Xs, 'surgery')
                    if is_new:
                        new.append(census.types[k])
                        if verbose:
                            print('new type r=%.10f k=%d coindex=%d' % (census.types[k]['r'], census.types[k]['k'],
                                                                         census.types[k]['coindex']), flush=True)
        if not new:
            break
        queue = new
    return census


def graphs_on(n, min_degree=2):
    """All graphs on exactly n vertices (n <= 7, from the networkx atlas) with minimum degree >= min_degree,
    up to isomorphism."""
    from networkx.generators.atlas import graph_atlas_g
    return [G for G in graph_atlas_g() if G.number_of_nodes() == n and G.number_of_edges() > 0
            and min(dict(G.degree()).values()) >= min_degree]


def graph_start(X, A, rng, delta=0.05, eps=0.05):
    """Seed for the KKT solve of graph A: minimize sum_A (c_e - t)^2 + sum_{not A} relu(c_f - t + delta)^2
    + eps * t over (X, t) by L-BFGS (X normalized inside)."""
    from scipy.optimize import minimize
    n, d = X.shape; A = list(A); Aset = set(A)
    I = np.array([i for i, j in A]); J = np.array([j for i, j in A])
    NI, NJ = np.array([(i, j) for i in range(n) for j in range(i + 1, n) if (i, j) not in Aset]).reshape(-1, 2).T \
        if n * (n - 1) // 2 > len(A) else (np.array([], int), np.array([], int))

    def f(z):
        Y = z[:-1].reshape(n, d); t = z[-1]
        nr = np.linalg.norm(Y, axis=1, keepdims=True); U = Y / nr
        gU = np.zeros_like(U)
        ce = (U[I] * U[J]).sum(1) - t
        val = (ce ** 2).sum() + eps * t; gt = -2 * ce.sum() + eps
        np.add.at(gU, I, 2 * ce[:, None] * U[J]); np.add.at(gU, J, 2 * ce[:, None] * U[I])
        if len(NI):
            cf = np.maximum((U[NI] * U[NJ]).sum(1) - t + delta, 0)
            val += (cf ** 2).sum(); gt += -2 * cf.sum()
            np.add.at(gU, NI, 2 * cf[:, None] * U[NJ]); np.add.at(gU, NJ, 2 * cf[:, None] * U[NI])
        gY = (gU - (gU * U).sum(1, keepdims=True) * U) / nr
        return val, np.concatenate([gY.ravel(), [gt]])

    z0 = np.concatenate([X.ravel(), [0.0]])
    res = minimize(f, z0, jac=True, method='L-BFGS-B', options=dict(maxiter=2000))
    return normalize(res.x[:-1].reshape(n, d))


def solve_graph(G, d, rng, starts=30, sigma_tol=1e-8):
    """Find critical configurations of n = |V(G)| points in S^{d-1} whose stressed contact graph is
    exactly G (all edges positively stressed, all non-edges strictly longer). Returns list of X."""
    n = G.number_of_nodes(); A = sorted(tuple(sorted(e)) for e in G.edges())
    sols = []
    for k in range(starts):
        X = random_config(n, d, rng)
        if k % 2 == 1:
            X = graph_start(X, A, rng)
        Xs, c, s, nf = solve_kkt(X, A, iters=200)
        ch = check_solution(Xs, c, s, A, nf)
        if ch and s.min() > sigma_tol:
            sols.append(normalize(Xs))
    return sols


def rattler_room(X):
    """d = 3 only. Largest angular distance (degrees) from a point of S^2 to the configuration, i.e. the
    covering radius; a rattler fits strictly iff this exceeds the contact angle."""
    from scipy.spatial import SphericalVoronoi
    if np.linalg.matrix_rank(np.asarray(X, float), 1e-8) < 3:
        return 90.0   # configuration on a great circle: the poles are at distance 90 degrees
    sv = SphericalVoronoi(np.asarray(X, float), radius=1.0, threshold=1e-8)
    V = sv.vertices / np.linalg.norm(sv.vertices, axis=1, keepdims=True)
    return float(np.degrees(np.arccos(np.clip((V @ np.asarray(X).T).max(axis=1), -1, 1)).max()))


def rattler_components(X, c, n_grid=20000, return_diameters=False):
    """d = 3 only. Number of connected components of the open rattler region {y : <y, x_i> < c for all i},
    computed on a Fibonacci grid (numerical)."""
    from scipy.spatial import cKDTree
    import networkx as nx
    k = np.arange(n_grid) + 0.5
    phi = np.arccos(1 - 2 * k / n_grid); th = np.pi * (1 + 5 ** 0.5) * k
    P = np.column_stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)])
    free = P[(P @ np.asarray(X).T).max(axis=1) < c - 1e-9]
    if len(free) == 0:
        return 0
    h = 2.5 * np.sqrt(4 * np.pi / n_grid)
    pairs = cKDTree(free).query_pairs(h)
    G = nx.Graph(); G.add_nodes_from(range(len(free))); G.add_edges_from(pairs)
    comps = [cc for cc in nx.connected_components(G) if len(cc) >= 3]
    if return_diameters:
        diams = []
        for cc in comps:
            Q = free[list(cc)][:400]
            diams.append(float(np.degrees(np.arccos(np.clip((Q @ Q.T).min(), -1, 1)))))
        return len(comps), diams
    return len(comps)


def euler_check(census_by_n, N):
    """Signed count sum_types (-1)^coindex * labeled * C(N, n) * (N-n)! over nondegenerate types on n <= N
    vertices, with N - n = 0 or 1 rattlers, times the number of components of the rattler region (each assumed
    contractible). Returns the sum and the
    list of degenerate types it skipped. Compare chi(M_{0,N}) = prod_{k=2}^{N-2} (1-k)."""
    from math import comb, factorial
    tot = 0; skipped = []
    for n, C in census_by_n.items():
        m = N - n
        if m < 0:
            continue
        for t in C.sorted():
            ncomp = 1
            if m > 0:
                nc, diams = rattler_components(np.array(t['X']), t['c'], return_diameters=True)
                theta = np.degrees(np.arccos(t['c']))
                if any(dm >= theta for dm in diams):
                    skipped.append((n, t['r'], 'rattler component can hold two rattlers')); continue
                if nc < m:
                    continue
                ncomp = comb(nc, m)   # m rattlers in distinct components (labels: comb(N,n) * m!)
            if t['nzero'] > 0 or t.get('labeled') is None:
                skipped.append((n, t['r'])); continue
            tot += (-1) ** t['coindex'] * t['labeled'] * comb(N, n) * factorial(m) * ncomp
    return tot, skipped


def chi_M0N(N):
    p = 1
    for k in range(2, N - 1):
        p *= (1 - k)
    return p
