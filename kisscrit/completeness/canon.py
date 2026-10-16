"""Canonical codes of embedded plane graphs (rotation systems), up to orientation-preserving and -reversing
isomorphism.

Status: trusted. Part of the trusted base (transfer prompt v1.2, code standards).
Implements: the canonical code of an embedded plane graph (L12: deduplication of candidates); the checker reimplements it independently in kisscrit/core/checker.py (V4).

For every directed edge (u, w) and both orientations, number the vertices in breadth-first order, visiting the
neighbours of each vertex in rotation order starting from the edge by which it was first reached (the start
vertex starts at w).  The code is the sequence of neighbour numbers; the canonical code is the lexicographic
minimum over all starts.  Two rotation systems have the same canonical code iff they are isomorphic as
embedded graphs, possibly after a reflection.  A reflection maps critical configurations to critical
configurations, so certifying one representative certifies both.
"""


def _code(rot, u, w, rev):
    """The code of the rotation system from the start edge (u, w) in the given orientation: per vertex in breadth-first order, its degree then its neighbours' numbers; None if the graph is disconnected."""
    n = len(rot)
    num = {u: 0}; order = [u]; first = {u: w}; out = []
    i = 0
    while i < len(order):
        v = order[i]; i += 1
        r = rot[v][::-1] if rev else rot[v]
        k = r.index(first[v]); seq = r[k:] + r[:k]
        out.append(len(seq))
        for x in seq:
            if x not in num:
                num[x] = len(order); order.append(x); first[x] = v
            out.append(num[x])
    return tuple(out) if len(order) == n else None


def canonical_code(rot):
    """rot[v] = list of neighbours of v in counterclockwise order."""
    best = None
    for u in range(len(rot)):
        for w in rot[u]:
            for rev in (False, True):
                c = _code(rot, u, w, rev)
                if c is not None and (best is None or c < best): best = c
    return best


def code_string(code):
    """Dot-joined text form of a code."""
    return '.'.join(map(str, code))


def rotation_of(S):
    """rot[v] = neighbours of v in counterclockwise order, from a Structure."""
    return [[w for (w, fi, pos) in S.rot[v]] for v in range(S.n)]
