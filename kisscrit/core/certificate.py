"""Reading certificates written by the completeness branch and bound (format: docs/certificate-format.md).

Status: trusted. Part of the trusted base (listed in docs/code-tour.md section 0; code standards in docs/decisions.md, 2026-10-16).
Implements: the certificate replay of docs/certificate-format.md section 3 (bisection and grid contraction, reproduced bit for bit).

This module only decodes and replays the proof tree; it uses no floating-point computation except the two
box operations the format defines (midpoint bisection and grid contraction), which are reproduced exactly.
"""
import base64, zlib
import numpy as np

QMAX = 65536
LEAF = {1: 'invalid', 2: 'vsum', 3: 'linear', 4: 'gordan-inherited', 5: 'gordan', 6: 'motzkin', 7: 'terminal',
        8: 'open'}
EV_ROUND, EV_SPLIT, EV_CONTRACT = 15, 16, 64


def unpack(s, dtype):
    """Decode a base64 + zlib stream into a little-endian numpy array of the given dtype."""
    return np.frombuffer(zlib.decompress(base64.b64decode(s)), dtype=np.dtype(dtype).newbyteorder('<'))


def grid_value(lo0, hi0, q):
    """Exactly as in the search: q == 0 -> lo0, q == QMAX -> hi0, else lo0 + (hi0 - lo0) * (q / 65536)."""
    if q <= 0: return lo0
    if q >= QMAX: return hi0
    return lo0 + (hi0 - lo0) * (q / 65536.0)


def replay(chart):
    """Yield (kind, lo, hi, extra) for every node of one chart's proof tree, in preorder.
    kind: 'split' (extra = k), 'contract' (extra = (slabs, new box); a slab is (lo, hi, cut_lo, cut_hi, k, side, mu):
    the part of the box [cut_lo, cut_hi] below (side -1) or above (side +1) the new bound of coordinate k, and
    mu the multiplier hint {vertex: value} that the search used for it; a 'linear' leaf carries its hint too),
    or a leaf name.
    Raises ValueError if the event stream does not describe a complete tree."""
    ev = unpack(chart['events'], np.uint8); q = unpack(chart['q'], np.uint16)
    mus = unpack(chart['mu'], np.float32) if 'mu' in chart else np.zeros(0, np.float32)
    pos = {'q': 0, 'm': 0}

    def read_mu():
        """A sparse multiplier hint: count and rows from the q stream, float32 values from the mu stream."""
        cnt = int(q[pos['q']]); pos['q'] += 1
        rows = [int(x) for x in q[pos['q']:pos['q'] + cnt]]; pos['q'] += cnt
        vals = [float(x) for x in mus[pos['m']:pos['m'] + cnt]]; pos['m'] += cnt
        return dict(zip(rows, vals))
    lo0 = np.array(chart['lo0'], float); hi0 = np.array(chart['hi0'], float)
    stack = [(lo0, hi0)]; i = 0; n = len(lo0)
    while stack:
        lo, hi = stack.pop()
        lo = lo.copy(); hi = hi.copy()
        # contraction rounds of this node, each relative to the box at the start of the round
        while i < len(ev) and ev[i] == EV_ROUND:
            i += 1
            base_lo, base_hi = lo.copy(), hi.copy()
            nchg = 0; hints = {}
            while i < len(ev) and ev[i] >= EV_CONTRACT:
                k = int(ev[i]) - EV_CONTRACT; i += 1; nchg += 1
                qlo = int(q[pos['q']]); qhi = QMAX - int(q[pos['q'] + 1]); pos['q'] += 2
                if qlo > 0: hints[(k, -1)] = read_mu()
                if qhi < QMAX: hints[(k, +1)] = read_mu()
                nlo = grid_value(base_lo[k], base_hi[k], qlo); nhi = grid_value(base_lo[k], base_hi[k], qhi)
                if not (base_lo[k] <= nlo <= nhi <= base_hi[k]):
                    raise ValueError('contraction outside the box')
                lo[k], hi[k] = nlo, nhi
            if nchg == 0: raise ValueError('empty contraction round')
            # removed region = round box minus new box, as slabs taken one coordinate at a time
            # each slab is listed with the box it was cut from (cur) and its coordinate and side
            removed = []; cur_lo, cur_hi = base_lo.copy(), base_hi.copy()
            for k in range(n):
                if cur_lo[k] < lo[k]:
                    s_lo, s_hi = cur_lo.copy(), cur_hi.copy(); s_hi[k] = lo[k]
                    removed.append((s_lo, s_hi, cur_lo.copy(), cur_hi.copy(), k, -1, hints.get((k, -1))))
                    cur_lo[k] = lo[k]
                if hi[k] < cur_hi[k]:
                    s_lo, s_hi = cur_lo.copy(), cur_hi.copy(); s_lo[k] = hi[k]
                    removed.append((s_lo, s_hi, cur_lo.copy(), cur_hi.copy(), k, +1, hints.get((k, +1))))
                    cur_hi[k] = hi[k]
            yield 'contract', base_lo, base_hi, (removed, (lo.copy(), hi.copy()))
        if i >= len(ev):
            raise ValueError('event stream ended inside the tree')
        e = int(ev[i]); i += 1
        if e in LEAF:
            yield LEAF[e], lo, hi, (read_mu() if e == 3 and 'mu' in chart else None)
        elif EV_SPLIT <= e < EV_CONTRACT:
            k = e - EV_SPLIT
            mid = 0.5 * (lo[k] + hi[k])
            lo_l, hi_l = lo.copy(), hi.copy(); hi_l[k] = mid
            lo_u, hi_u = lo.copy(), hi.copy(); lo_u[k] = mid
            yield 'split', lo, hi, k
            stack.append((lo_l, hi_l)); stack.append((lo_u, hi_u))      # upper child is next in the stream
        else:
            raise ValueError('unknown event %d' % e)
    if i != len(ev) or pos['q'] != len(q) or pos['m'] != len(mus):
        raise ValueError('trailing events')
