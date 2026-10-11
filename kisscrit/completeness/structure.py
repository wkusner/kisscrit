"""Combinatorial structure of a candidate stressed graph: oriented faces, corners, rotation systems,
and the face-parameter vector.

A Structure is built from oriented face cycles (each face listed counterclockwise as seen from outside the
sphere).  At a vertex v of a face (... a, v, b ...) the corner spans counterclockwise from the edge v->b to
the edge v->a, so in the counterclockwise rotation at v the neighbor a follows b.

Parameter vector p (radians):
    p[0] = theta;  one parameter per quadrilateral (corner at boundary positions 0, 2);
    two per pentagon (corners at boundary positions 1, 2).
"""
import itertools
import numpy as np
import networkx as nx
from . import faces as fc


class Structure:
    def __init__(self, face_cycles, n=12):
        self.n = n
        self.faces = [list(f) for f in face_cycles]
        self.sizes = [len(f) for f in self.faces]
        assert all(3 <= k <= 5 for k in self.sizes), 'faces must have 3..5 corners'
        E = set()
        for f in self.faces:
            for i in range(len(f)):
                E.add(tuple(sorted((f[i], f[(i + 1) % len(f)]))))
        self.edges = sorted(E)
        self.eidx = {e: k for k, e in enumerate(self.edges)}
        assert n - len(self.edges) + len(self.faces) == 2, 'not a sphere (Euler characteristic)'
        # parameters
        self.pidx = []; npar = 1
        for k in self.sizes:
            if k == 3: self.pidx.append(())
            elif k == 4: self.pidx.append((npar,)); npar += 1
            else: self.pidx.append((npar, npar + 1)); npar += 2
        self.npar = npar
        # corners and rotation systems
        succ = [dict() for _ in range(n)]      # succ[v][b] = (a, face, position)
        for fi, f in enumerate(self.faces):
            k = len(f)
            for pos in range(k):
                v, a, b = f[pos], f[pos - 1], f[(pos + 1) % k]
                assert b not in succ[v], 'inconsistent orientation'
                succ[v][b] = (a, fi, pos)
        self.rot = []                           # per vertex: list of (neighbor, face, pos) ccw, corner after neighbor
        for v in range(n):
            assert len(succ[v]) >= 2, 'vertex %d has degree < 2' % v
            start = min(succ[v]); seq = []; w = start
            for _ in range(len(succ[v])):
                a, fi, pos = succ[v][w]
                seq.append((w, fi, pos)); w = a
            assert w == start and len(seq) == len(succ[v]), 'rotation does not close'
            self.rot.append(seq)

    # ------------------------------------------------------------------ corners
    def corners(self, p):
        """List per face of its corners (tuple), and a validity flag (array or bool)."""
        th = p[0]; out = []; ok = True
        a = fc.alpha(th)
        for fi, k in enumerate(self.sizes):
            if k == 3:
                out.append((a, a, a))
            elif k == 4:
                A = p[self.pidx[fi][0]]; O, o = fc.quad_other(th, A); ok = ok & o
                out.append((A, O, A, O))
            else:
                i1, i2 = self.pidx[fi]
                cs, o = fc.pentagon(th, p[i1], p[i2]); ok = ok & o
                if cs is None: return out, False
                out.append(cs)
        return out, ok

    def vertex_sums(self, C):
        return [sum_list([C[fi][pos] for (_, fi, pos) in self.rot[v]]) for v in range(self.n)]

    def directions(self, C):
        """phi[v][w]: angle of the edge v->w at v, measured ccw from the first neighbor in rot[v]."""
        phis = []
        for v in range(self.n):
            d = {}; acc = None
            for (w, fi, pos) in self.rot[v]:
                d[w] = 0.0 * C[fi][pos] if acc is None else acc
                acc = C[fi][pos] if acc is None else acc + C[fi][pos]
            phis.append(d)
        return phis

    # ------------------------------------------------------------------ initial box
    def initial_box(self, th_lo, th_hi):
        a_lo = float(fc.alpha(th_lo)); a_hi = float(fc.alpha(th_hi))
        lo = np.zeros(self.npar); hi = np.zeros(self.npar)
        lo[0], hi[0] = th_lo, th_hi
        for fi, k in enumerate(self.sizes):
            for i in self.pidx[fi]:
                if k == 4: lo[i], hi[i] = a_lo, 2 * a_hi
                else: lo[i], hi[i] = a_lo, np.pi
        return lo, hi


def sum_list(xs):
    s = xs[0]
    for x in xs[1:]: s = s + x
    return s


# ---------------------------------------------------------------------- builders
def orient_triangulation(F):
    """Consistently oriented (outward ccw) triangles of a triangulation given as vertex triples."""
    G = nx.Graph()
    for a, b, c in F: G.add_edges_from([(a, b), (b, c), (a, c)])
    ok, emb = nx.check_planarity(G)
    tris = set(); out = []
    for u, v in emb.edges():
        f = emb.traverse_face(u, v)
        if len(f) == 3:
            key = frozenset(f)
            if key not in tris: tris.add(key); out.append(tuple(f))
    return out          # orientation is consistent; global sign irrelevant for our purposes


def merge_faces(oriented_tris, groups):
    """Merge groups of oriented triangles (tuples of indices into oriented_tris) into oriented polygons."""
    faces = []
    for g in groups:
        dedges = set()
        for t in g:
            a, b, c = oriented_tris[t]
            dedges |= {(a, b), (b, c), (c, a)}
        bd = {(u, v) for (u, v) in dedges if (v, u) not in dedges}
        nxt = {}
        for u, v in bd:
            assert u not in nxt, 'merged face is not a simple polygon'
            nxt[u] = v
        start = min(nxt); cyc = [start]; w = nxt[start]
        while w != start:
            cyc.append(w); w = nxt[w]
        assert len(cyc) == len(bd)
        faces.append(cyc)
    return faces


def from_configuration(X, edges):
    """Oriented faces of the geodesic plane graph (X, edges) on S^2 (edges must not cross)."""
    X = np.asarray(X, float); n = len(X)
    nb = [[] for _ in range(n)]
    for i, j in edges: nb[i].append(j); nb[j].append(i)
    rot = []
    for v in range(n):
        x = X[v]; e1 = X[nb[v][0]] - (X[nb[v][0]] @ x) * x; e1 /= np.linalg.norm(e1); e2 = np.cross(x, e1)
        ang = []
        for w in nb[v]:
            t = X[w] - (X[w] @ x) * x
            ang.append(np.arctan2(t @ e2, t @ e1))
        rot.append([w for _, w in sorted(zip(ang, nb[v]))])     # ccw around outward normal
    # trace faces: corner at v spans ccw from v->b to v->a, so a = ccw successor of b at v;
    # walking a -> v -> b: b = ccw predecessor of a at v
    seen = set(); faces = []
    for v in range(n):
        for a in nb[v]:
            if (a, v) in seen: continue
            cyc = []; u, w = a, v
            while (u, w) not in seen:
                seen.add((u, w)); cyc.append(w)
                r = rot[w]; i = r.index(u); b = r[(i - 1) % len(r)]
                u, w = w, b
            faces.append(cyc)
    return faces
