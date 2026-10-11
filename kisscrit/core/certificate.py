"""Reading certificates written by the completeness branch and bound (format: docs/certificate-format.md).

This module only decodes and replays the proof tree; it uses no floating-point computation except the two
box operations the format defines (midpoint bisection and grid contraction), which are reproduced exactly.
"""
import base64, zlib
import numpy as np

QMAX = 65536
LEAF = {1: 'invalid', 2: 'vsum', 3: 'linear', 4: 'gordan-inherited', 5: 'gordan', 6: 'motzkin', 7: 'terminal',
        8: 'open'}
EV_SPLIT, EV_CONTRACT = 16, 64


def unpack(s, dtype):
    return np.frombuffer(zlib.decompress(base64.b64decode(s)), dtype=np.dtype(dtype).newbyteorder('<'))


def grid_value(lo0, hi0, q):
    """Exactly as in the search: q == 0 -> lo0, q == QMAX -> hi0, else lo0 + (hi0 - lo0) * (q / 65536)."""
    if q <= 0: return lo0
    if q >= QMAX: return hi0
    return lo0 + (hi0 - lo0) * (q / 65536.0)


def replay(chart):
    """Yield (kind, lo, hi, extra) for every node of one chart's proof tree, in preorder.
    kind: 'split' (extra = k), 'contract' (extra = (removed slabs list, new box)), or a leaf name.
    Raises ValueError if the event stream does not describe a complete tree."""
    ev = unpack(chart['events'], np.uint8); q = unpack(chart['q'], np.uint16)
    lo0 = np.array(chart['lo0'], float); hi0 = np.array(chart['hi0'], float)
    stack = [(lo0, hi0)]; i = 0; iq = 0; n = len(lo0)
    while stack:
        lo, hi = stack.pop()
        lo = lo.copy(); hi = hi.copy()
        # contractions of this node
        base_lo, base_hi = lo.copy(), hi.copy()
        removed = []
        while i < len(ev) and ev[i] >= EV_CONTRACT:
            k = int(ev[i]) - EV_CONTRACT; i += 1
            qlo = int(q[iq]); qhi = QMAX - int(q[iq + 1]); iq += 2
            nlo = grid_value(base_lo[k], base_hi[k], qlo); nhi = grid_value(base_lo[k], base_hi[k], qhi)
            if not (base_lo[k] <= nlo <= nhi <= base_hi[k]):
                raise ValueError('contraction outside the box')
            lo[k], hi[k] = nlo, nhi
        if (lo != base_lo).any() or (hi != base_hi).any():
            # the removed region is the base box minus the new box: slabs, one coordinate at a time
            cur_lo, cur_hi = base_lo.copy(), base_hi.copy()
            for k in range(n):
                if cur_lo[k] < lo[k]:
                    s_lo, s_hi = cur_lo.copy(), cur_hi.copy(); s_hi[k] = lo[k]; removed.append((s_lo, s_hi))
                    cur_lo[k] = lo[k]
                if hi[k] < cur_hi[k]:
                    s_lo, s_hi = cur_lo.copy(), cur_hi.copy(); s_lo[k] = hi[k]; removed.append((s_lo, s_hi))
                    cur_hi[k] = hi[k]
            yield 'contract', base_lo, base_hi, (removed, (lo.copy(), hi.copy()))
        if i >= len(ev):
            raise ValueError('event stream ended inside the tree')
        e = int(ev[i]); i += 1
        if e in LEAF:
            yield LEAF[e], lo, hi, None
        elif EV_SPLIT <= e < EV_CONTRACT:
            k = e - EV_SPLIT
            mid = 0.5 * (lo[k] + hi[k])
            lo_l, hi_l = lo.copy(), hi.copy(); hi_l[k] = mid
            lo_u, hi_u = lo.copy(), hi.copy(); lo_u[k] = mid
            yield 'split', lo, hi, k
            stack.append((lo_l, hi_l)); stack.append((lo_u, hi_u))      # upper child is next in the stream
        else:
            raise ValueError('unknown event %d' % e)
    if i != len(ev) or iq != len(q):
        raise ValueError('trailing events')
