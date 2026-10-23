"""Candidate stressed graphs: plane graphs on 12 vertices whose faces are 3-, 4- or 5-gons, obtained from a
triangulation T of the sphere by merging adjacent triangles (a quadrilateral = 2 triangles, a pentagon = a dual
path of 3 triangles).  Every plane graph with faces of size 3..5 arises from some T (triangulate its faces by
diagonals), so running over all 7595 triangulations with 12 vertices enumerates every candidate.

Status: trusted. Part of the trusted base (listed in docs/code-tour.md section 0; code standards in docs/decisions.md, 2026-10-16).
Implements: candidate generation (tilings of the 7,595 triangulations, L12-L14 of docs/lemmas.md; paper Appendix A, proof of Corollary 6.3). Known gap G1: the LP discard is floating point until exact Farkas certificates are logged (checklist F6).

Pruning used during enumeration (each is a necessary condition for a critical configuration with
60 deg < theta <= 63.4349 deg; see docs/completeness-runbook.md, section "What is proved, what is assumed"):
  * every vertex has degree 2..5; a degree-2 vertex has two straight pentagon corners;
  * at each vertex the corners can sum to 360 deg for some common alpha = alpha(theta) in a cell of
    [70.5288, 72] deg, with triangle corners = alpha, quadrilateral corners in [alpha, 2 alpha],
    pentagon corners in [alpha, 180] (alpha is cut into 64 cells; a bit mask of feasible cells is kept);
  * a pentagon cannot be a fan of three triangles around a vertex of degree 3 (that vertex would be interior).
After enumeration each candidate passes a linear relaxation (lp_feasible) that uses only these same facts plus
"a spherical polygon has corner sum > (k - 2) 180 deg", and is identified by its canonical code.
"""
import itertools, collections
import numpy as np
from scipy.optimize import linprog
from .structure import Structure, orient_triangulation, merge_faces
from .canon import canonical_code, rotation_of, code_string

A60 = float(np.degrees(np.arccos(1 / 3))); A63 = 72.0
NC = 64
CELLS = [(A60 + (A63 - A60) * i / NC, A60 + (A63 - A60) * (i + 1) / NC) for i in range(NC)]
FULL = (1 << NC) - 1


def _vmask(t, q, p):
    """Bit mask over the alpha cells [A60, A63] of the vertex types (t triangles, q rhombi, p pentagons at a vertex) that can close to 360 degrees: a cell is kept if d*alpha_lo <= 360 and t*alpha_hi + 2q*alpha_hi + 180p >= 360 (L3, L7, L13: triangle corners equal alpha, rhombus corners at most 2 alpha, pentagon corners at most 180). Degree 2 only as two pentagon corners (L4). Floating point with a permissive slack of 1e-9; the exact version is checklist item G2 (F6)."""
    d = t + q + p
    if d < 2 or d > 5: return 0
    if d == 2: return FULL if (t == 0 and q == 0 and p == 2) else 0
    m = 0
    for i, (lo, hi) in enumerate(CELLS):
        if d * lo <= 360 + 1e-9 and t * hi + q * 2 * hi + p * 180 >= 360 - 1e-9:
            m |= 1 << i
    return m


VMASK = {(t, q, p): _vmask(t, q, p) for t in range(6) for q in range(6) for p in range(6)}


def tilings(F):
    """Yield the partitions of the triangles F (list of sorted vertex triples) into faces (tuples of triangle
    indices) that pass the vertex pruning."""
    F = [tuple(sorted(t)) for t in F]; nF = len(F)
    adj = collections.defaultdict(list); edge_faces = collections.defaultdict(list)
    for fi, t in enumerate(F):
        for e in itertools.combinations(t, 2): edge_faces[e].append(fi)
    for e, (a, b) in edge_faces.items():
        adj[a].append(b); adj[b].append(a)
    vtris = collections.defaultdict(list)
    for fi, t in enumerate(F):
        for v in t: vtris[v].append(fi)
    order = []; seen = set(); done_v = set()
    while len(order) < nF:
        best = None
        for v in range(12):
            if v in done_v: continue
            rem = [f for f in vtris[v] if f not in seen]
            if best is None or len(rem) < best[0] or (len(rem) == best[0] and len(rem) > 0):
                if rem or best is None: best = (len(rem), v)
        v = best[1]; done_v.add(v)
        for f in vtris[v]:
            if f not in seen: seen.add(f); order.append(f)
    pos = {f: i for i, f in enumerate(order)}
    vlast = {v: max(pos[f] for f in vtris[v]) for v in range(12)}
    face_of = [-1] * nF; faces = []

    def corners_at(v):
        """The (triangle, rhombus, pentagon) corner counts at vertex v of the current partial tiling."""
        t = q = p = 0
        for fid in set(face_of[f] for f in vtris[v]):
            k = len(faces[fid]) + 2
            if k == 3: t += 1
            elif k == 4: q += 1
            else: p += 1
        return t, q, p

    def check(k, mask):
        """True iff the vertex types at the vertices of triangle k are still feasible under the mask after assigning it."""
        for v in range(12):
            if vlast[v] == k:
                mask &= VMASK.get(corners_at(v), 0)
                if not mask: return 0
        return mask

    def rec(k, mask):
        """Depth-first assignment of triangle k to a face (alone, or merged with a neighbour into a rhombus or pentagon), pruned by the vertex masks."""
        while k < nF and face_of[order[k]] != -1:
            mask = check(k, mask)
            if not mask: return
            k += 1
        if k == nF:
            yield [tuple(f) for f in faces]; return
        f = order[k]
        options = [(f,)]
        for g in adj[f]:
            if face_of[g] == -1:
                options.append((f, g))
                for h in adj[g]:
                    if h != f and face_of[h] == -1 and h not in adj[f]:
                        options.append((f, g, h))
                for h in adj[f]:
                    if h != g and face_of[h] == -1 and h not in adj[g] and g < h:
                        options.append((g, f, h))
        for opt in options:
            if len(opt) == 3:
                common = set(F[opt[0]]) & set(F[opt[1]]) & set(F[opt[2]])
                if any(len(vtris[v]) == 3 for v in common):
                    continue
            fid = len(faces); faces.append(opt)
            for g in opt: face_of[g] = fid
            m2 = check(k, mask)
            if m2: yield from rec(k + 1, m2)
            for g in opt: face_of[g] = -1
            faces.pop()

    yield from rec(0, FULL)


def lp_feasible(S):
    """Linear relaxation on corner angles (degrees): vertex sums 360; triangle corner = alpha in [A60, 72];
    quadrilateral corners in [alpha, 2 alpha]; pentagon corners in [alpha, 180]; corner sum of a k-gon
    >= (k - 2) 180.  Only facts valid for every critical configuration are used."""
    nv = 1; idx = {}
    for fi, k in enumerate(S.sizes):
        if k > 3:
            for q in range(k): idx[(fi, q)] = nv; nv += 1
    A_ub = []; b_ub = []
    for (fi, q), c in idx.items():
        r = np.zeros(nv); r[0] = 1; r[c] = -1; A_ub.append(r); b_ub.append(0)          # corner >= alpha
        if S.sizes[fi] == 4:
            r = np.zeros(nv); r[c] = 1; r[0] = -2; A_ub.append(r); b_ub.append(0)      # corner <= 2 alpha
    for fi, k in enumerate(S.sizes):
        if k > 3:
            r = np.zeros(nv)
            for q in range(k): r[idx[(fi, q)]] = -1
            A_ub.append(r); b_ub.append(-(k - 2) * 180.0)
    A_eq = np.zeros((S.n, nv)); b_eq = np.full(S.n, 360.0)
    for v in range(S.n):
        for (w, fi, pos) in S.rot[v]:
            if S.sizes[fi] == 3: A_eq[v, 0] += 1
            else: A_eq[v, idx[(fi, pos)]] += 1
    bounds = [(A60, 72.0)] + [(0, 180.0)] * (nv - 1)
    res = linprog(np.zeros(nv), A_ub=np.array(A_ub) if A_ub else None, b_ub=b_ub if b_ub else None,
                  A_eq=A_eq, b_eq=b_eq, bounds=bounds, method='highs')
    return res.status == 0


def candidates_of_triangulation(F, seen=None):
    """Yield (code_string, face_cycles) for each new LP-feasible candidate from triangulation F."""
    F = [tuple(sorted(t)) for t in F]
    OT = orient_triangulation(F); key = {frozenset(t): i for i, t in enumerate(OT)}
    seen = set() if seen is None else seen
    for tiling in tilings(F):
        groups = [tuple(key[frozenset(F[t])] for t in tris) for tris in tiling]
        cyc = merge_faces(OT, groups)
        S = Structure(cyc)
        code = code_string(canonical_code(rotation_of(S)))
        if code in seen: continue
        seen.add(code)
        if lp_feasible(S):
            yield code, [list(map(int, c)) for c in cyc]
