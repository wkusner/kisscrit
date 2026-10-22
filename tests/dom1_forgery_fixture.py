"""DOM-1 fixture: the DES-F2 rhombus forgery (docs/lemmas.md L32, V1b), rebuilt from a committed certificate line.

The forged certificate is the r1 saddle's candidate (E = 21, one chart, 11 faces, three pentagons and two rhombi) of
results/completeness/certificates/certificates-cond-11-of-40.jsonl, 0-based line 661 (DES-F2 writes "line 662",
counted from 1), with the start box of the coordinate of rhombus face 2 ([7, 9, 10, 11]) raised from the honest
[alpha(theta_lo), 2 alpha(theta_hi)] to [alpha(theta_lo), 2 pi + alpha(theta_lo)] and a proof tree of 1,791 events in
which every leaf is excluded by the checker's own direct tests (MAX_DEPTH = 0) on a box whose rhombus coordinate is
never split.  At 5da4bba check_certificate accepted it (verified True, no structural problem, 0 failed, 0 open, 0
terminal; counts vsum 744, invalid 126, vsum-tight 16, invalid-tight 10).

Stored, not regenerated: EVENTS_B64 is the zlib-compressed event stream (432 characters); forged_certificate() builds
the certificate from it in milliseconds.  regenerate_events() rebuilds the stream with the checker as oracle (about
35 s on 2 cores); tests/test_certificates.py runs it only with the environment variable DOM1_REGEN=1.
"""
import base64, json, math, os, zlib
import numpy as np

CERT_FILE = 'results/completeness/certificates/certificates-cond-11-of-40.jsonl'
CERT_LINE0 = 661                  # 0-based
RHOMBUS_FACE = 2                  # the rhombus [7, 9, 10, 11]; its coordinate is index 1 of the chart
N_EVENTS = 1791
EVENTS_B64 = (
    'eNqVVFuSxCAIrOojLJmo97/oIoqCYibzY5RnQ6AJ151ywR8+IL4jIaPUl4gggqa24mLtVb1YwwVV52lmBBMABhqb'
    'Dl5IKi6717QlCb0GbgYjipjYBiRXfn9crLoLo2ZzNujKNLSJ1YAciCsI8WsAdq8y09YZwrkLjlxx8KEFWxVrkFOx'
    'bTDoMlYfUeVqX9X1q4dN1158E5A+oeTCCTCpgCagXq17BOamFWNSGvg8/0kutWI5+NoHoF4wHS5fm/5+Lb6XPgaq'
    '9tb4Hjq0gPCZg9StRWEhNEU04fXGqEZHJHTx67PsSIrWxDYaodC6BClxEG9kEizzgTyOdAG3nl9269a+fU6s8Uwl'
    'hkn65d6bZOYVC7Su21nGPlfKecsXLyjCbRwe/uzreYzFCDiTfgZvgT1QhsO/FrAtw3mp1hrwDzSuUdM='
)


def pack(arr):
    return base64.b64encode(zlib.compress(np.ascontiguousarray(arr).tobytes(), 9)).decode()


def honest_theta():
    """[theta_min rounded down, arctan 2 rounded up], the range the checker's V3 asks for (decided in Arb)."""
    from kisscrit.core import checker as ck
    r = ck.arb(ck.R_MIN)
    return [ck.down(2 * (r / (1 + r)).asin()), ck.up(ck.arb(2).atan())]


def load_base():
    """(code, certificate) of the committed line, or None when the file is not in the checkout."""
    if not os.path.exists(CERT_FILE): return None
    with open(CERT_FILE) as f:
        for ln, line in enumerate(f):
            if ln == CERT_LINE0:
                d = json.loads(line); return d['code'], d['certificate']
    return None


def start_box(cert0, forged=True):
    """(lo, hi) of chart 0: the checker's own a-priori box at the honest theta (V1 passes by construction); with
    forged=True the rhombus coordinate's hi is 2 pi + alpha(theta_lo) rounded down, as DES-F2 did."""
    from kisscrit.core import checker as ck, model as am
    th = honest_theta(); G = am.Graph(cert0['faces'], cert0['charts'][0]['modes'])
    lo, hi, _ = ck.apriori_box(G, th); lo = list(lo); hi = list(hi); lo[0], hi[0] = th[0], th[1]
    j = G.pidx[RHOMBUS_FACE][0]
    if forged: hi[j] = 2 * math.pi + lo[j]
    return lo, hi, j


def forged_certificate(events=None):
    """(code, certificate, j): the forged certificate built from the stored events (or from `events`)."""
    code, cert0 = load_base()
    th = honest_theta(); lo, hi, j = start_box(cert0)
    ev = events if events is not None else np.frombuffer(zlib.decompress(base64.b64decode(EVENTS_B64)), np.uint8)
    ch = dict(modes=cert0['charts'][0]['modes'], lo0=lo, hi0=hi, events=pack(ev), q=pack(np.zeros(0, np.uint16)),
              mu=pack(np.zeros(0, np.float32)), n_events=len(ev), complete=True)
    cert = dict(version=2, theta=list(th), tol=cert0['tol'], faces=cert0['faces'], complete=True, charts=[ch])
    return code, cert, j


def regenerate_events():
    """The event stream, rebuilt as DES-F2's c4b did: bisect every coordinate except the rhombus one (widest first,
    upper child first, format section 3), closing each box by the checker's direct tests (MAX_DEPTH = 0, restored
    afterwards).  Returns (ok, events)."""
    from kisscrit.core import checker as ck, model as am
    from kisscrit.core.certificate import EV_SPLIT
    code, cert0 = load_base()
    lo, hi, j = start_box(cert0)
    G = am.Graph(cert0['faces'], cert0['charts'][0]['modes']); K = ck.Checker(G)
    events = []; failed = []
    saved = ck.MAX_DEPTH; ck.MAX_DEPTH = 0
    try:
        def rec(blo, bhi, depth):
            try: closed = K.exclude(list(blo), list(bhi))
            except ArithmeticError: closed = False
            if closed: events.append(2); return True
            w = [b - a for a, b in zip(blo, bhi)]; w[j] = 0.0
            k = int(np.argmax(w))
            if w[k] < 1e-9 or depth > 60 or failed:
                failed.append((list(blo), list(bhi))); events.append(7); return False
            mid = 0.5 * (blo[k] + bhi[k]); events.append(EV_SPLIT + k)
            u = list(blo); u[k] = mid; dn = list(bhi); dn[k] = mid
            a = rec(u, list(bhi), depth + 1); b = rec(list(blo), dn, depth + 1)
            return a and b
        ok = rec(lo, hi, 0)
    finally:
        ck.MAX_DEPTH = saved
    return ok, np.array(events, np.uint8)
