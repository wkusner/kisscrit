"""Certificates: the search's proof tree replays exactly; the independent checker accepts an honest certificate
for the local maximum M and rejects every mutation of it (spec v0.2 section 4, gaps V1-V5); the manifest step
detects a missing line.  Each test names the gap it covers."""
import base64, json, zlib, collections, os, gzip
import numpy as np
import pytest
import kisscrit as kc
from kisscrit.completeness.structure import Structure, from_configuration
from kisscrit.completeness import runner
from kisscrit.core.certificate import replay, unpack
from kisscrit.core import checker as ck
from kisscrit.core.checker import check_certificate
from kisscrit.completeness.canon import canonical_code as search_code, rotation_of, code_string

THETA = (np.radians(60.0033), np.radians(63.4349489))


def _pack(arr):
    return base64.b64encode(zlib.compress(np.ascontiguousarray(arr).tobytes(), 9)).decode()


def _rejected(cert, code):
    try:
        r = check_certificate(cert, code=code)
    except ValueError:
        return True
    return not r['verified']


@pytest.fixture(scope='module')
def m_certificate():
    D = json.load(open('data/n12_above1.json'))['types']
    t = D[1]                                            # M, r = 1.0010121
    X, c, s, nf = kc.solve_kkt(np.array(t['X']), [tuple(p) for p in t['pairs']])
    a = kc.analyze(X)
    sup = [tuple(p) for p, w in zip(a['pairs'], a['support']) if w]
    S = Structure(from_configuration(X, sup))
    term, st, cert = runner.run(S, *THETA, tol=1e-6, certificate=True)
    code = code_string(search_code(rotation_of(runner.oriented(S))))
    return term, st, json.loads(json.dumps(cert)), code


def test_replay_matches_search(m_certificate):
    term, st, cert, code = m_certificate
    cnt = collections.Counter(); tb = []
    for ch in cert['charts']:
        for kind, lo, hi, extra in replay(ch):
            cnt[kind] += 1
            if kind == 'terminal': tb.append((lo, hi))
    assert cnt['vsum'] == st['vsum'] and cnt['invalid'] == st['invalid'] and cnt['terminal'] == st['terminal']
    assert len(tb) == len(term) and all(np.array_equal(a[0], b[1]) and np.array_equal(a[1], b[2])
                                        for a, b in zip(tb, term))


def test_checker_accepts_honest_certificate(m_certificate):
    term, st, cert, code = m_certificate
    r = check_certificate(cert, code=code, gen_faces=cert['faces'])
    assert r['verified'] and not r['structural'] and not r['failed'] and len(r['terminal']) >= 1


# ---------------------------------------------------------------- mutations
def test_v1_shrunk_starting_box(m_certificate):
    cert = json.loads(json.dumps(m_certificate[2])); code = m_certificate[3]
    ch = cert['charts'][0]
    j = next(j for j in range(1, len(ch['lo0'])) if ch['lo0'][j] != ch['hi0'][j])
    from kisscrit.core import model as am
    own_lo, own_hi, fixed = ck.apriori_box(am.Graph(cert['faces'], ch['modes']), cert['theta'])
    assert ch['lo0'][j] <= own_lo[j] and ch['hi0'][j] >= own_hi[j]
    ch['lo0'][j] = float(np.nextafter(own_lo[j], 10))          # one ulp inside the a-priori box
    assert _rejected(cert, code)
    cert = json.loads(json.dumps(m_certificate[2])); ch = cert['charts'][0]
    ch['hi0'][j] = float(np.nextafter(own_hi[j], -10))
    assert _rejected(cert, code)


def test_v1_fixed_coordinate_must_be_a_straight_corner(m_certificate):
    cert = json.loads(json.dumps(m_certificate[2])); code = m_certificate[3]
    ch = cert['charts'][0]
    j = next(j for j in range(1, len(ch['lo0'])) if ch['lo0'][j] != ch['hi0'][j])
    ch['lo0'][j] = ch['hi0'][j] = 0.5 * (ch['lo0'][j] + ch['hi0'][j])
    assert _rejected(cert, code)


def test_v2_missing_required_chart(m_certificate):
    cert = json.loads(json.dumps(m_certificate[2])); code = m_certificate[3]
    assert len(cert['charts']) >= 1
    cert['charts'] = cert['charts'][1:]
    assert _rejected(cert, code)
    cert = json.loads(json.dumps(m_certificate[2]))
    cert['charts'][0]['complete'] = False
    assert _rejected(cert, code)


def test_v3_theta_range(m_certificate):
    cert = json.loads(json.dumps(m_certificate[2])); code = m_certificate[3]
    hi = float(np.nextafter(ck.up(ck.arb(2).atan()), 0))           # just below arctan 2
    cert['theta'][1] = hi
    for ch in cert['charts']: ch['hi0'][0] = hi
    assert _rejected(cert, code)
    cert = json.loads(json.dumps(m_certificate[2]))
    lo = float(ck.down(2 * (ck.arb(ck.R_MIN) / (1 + ck.arb(ck.R_MIN))).asin())) + 1e-12
    cert['theta'][0] = lo
    for ch in cert['charts']: ch['lo0'][0] = lo
    assert _rejected(cert, code)


def test_v4_code_and_faces(m_certificate):
    cert = json.loads(json.dumps(m_certificate[2])); code = m_certificate[3]
    assert _rejected(cert, code[:-1] + ('1' if code[-1] != '1' else '2'))
    cert['faces'][0] = cert['faces'][0][::-1]                 # reversed cycle: not the same rotation system
    assert _rejected(cert, code)
    cert = json.loads(json.dumps(m_certificate[2]))
    other = json.loads(json.dumps(cert['faces'])); other[0], other[1] = other[1], other[0]
    r = check_certificate(cert, code=code, gen_faces=other)   # a permutation of the face list is the same graph
    assert not r['structural']
    other[0] = other[0][::-1]
    r = check_certificate(cert, code=code, gen_faces=other, structural=True)
    assert any(p.startswith('V4') for p in r['structural'])


def test_v4_checker_code_matches_search_code_on_gen_records():
    """The checker's canonical code agrees with the search's on 50 stored gen records (and with their codes)."""
    f = sorted(__import__('glob').glob('results/completeness/gen/gen-cond-*.jsonl.gz'))
    if not f: pytest.skip('no gen records in the repository')
    n = 0
    with gzip.open(f[0], 'rt') as fh:
        for line in fh:
            d = json.loads(line)
            if 'code' not in d: continue
            S = Structure(d['faces'])
            assert ck.canonical_code(d['faces']) == code_string(search_code(rotation_of(S))) == d['code']
            n += 1
            if n >= 50: break
    assert n == 50


def test_v5_event_stream_truncated_or_extended(m_certificate):
    cert = json.loads(json.dumps(m_certificate[2])); code = m_certificate[3]
    ch = cert['charts'][0]; ev = unpack(ch['events'], np.uint8).copy()
    ch['events'] = _pack(ev[:-1]); assert _rejected(cert, code)
    ch['events'] = _pack(np.concatenate([ev, np.array([2], np.uint8)])); assert _rejected(cert, code)


@pytest.mark.parametrize('code', [1, 2, 5])
def test_forged_terminal_leaf(m_certificate, code):
    """Gap 6: a terminal leaf on the census point relabelled as excluded."""
    cert = json.loads(json.dumps(m_certificate[2])); ccode = m_certificate[3]
    for ch in cert['charts']:
        ev = unpack(ch['events'], np.uint8).copy()
        idx = np.nonzero(ev == 7)[0]
        if len(idx):
            ev[idx[0]] = code; ch['events'] = _pack(ev); break
    assert _rejected(cert, ccode)


def test_forged_open_leaf_and_linear_hint(m_certificate):
    """Gap 6: an open leaf (budget run) relabelled vsum; a terminal leaf relabelled linear with a forged hint."""
    D = json.load(open('data/n12_above1.json'))['types']
    t = D[1]
    X, c, s, nf = kc.solve_kkt(np.array(t['X']), [tuple(p) for p in t['pairs']])
    a = kc.analyze(X)
    sup = [tuple(p) for p, w in zip(a['pairs'], a['support']) if w]
    S = Structure(from_configuration(X, sup))
    term, st, cert = runner.run(S, *THETA, tol=1e-6, max_boxes=30, certificate=True)
    cert = json.loads(json.dumps(cert)); code = m_certificate[3]
    assert st['aborted']
    r = check_certificate(cert, code=code)
    assert not r['verified'] and (r['open'] > 0 or r['not_run'] > 0)
    changed = False
    for ch in cert['charts']:
        if ch.get('not_run'): continue
        ev = unpack(ch['events'], np.uint8).copy(); idx = np.nonzero(ev == 8)[0]
        if len(idx): ev[idx[0]] = 2; ch['events'] = _pack(ev); changed = True; break
    if changed:
        assert _rejected(cert, code)
    # a terminal leaf turned into 'linear' with a forged hint: the hint must be rejected and the fallback too
    cert = json.loads(json.dumps(m_certificate[2]))
    for ch in cert['charts']:
        ev = unpack(ch['events'], np.uint8).copy(); idx = np.nonzero(ev == 7)[0]
        if len(idx):
            ev[idx[0]] = 3; ch['events'] = _pack(ev)
            q = unpack(ch['q'], np.uint16).copy(); mu = unpack(ch['mu'], np.float32).copy()
            # the hint is read after the leaf: count 1, vertex 0, value 1.0 appended at the end of the streams
            # is only correct if this leaf is the last hint consumer; so instead give it count 0 (no rows)
            # by inserting a zero count at the right position is involved; use an empty-mu hint via count 0
            pos = _hint_position(ch, idx[0])
            q = np.concatenate([q[:pos], np.array([0], np.uint16), q[pos:]]); ch['q'] = _pack(q)
            break
    assert _rejected(cert, code)


def _hint_position(ch, leaf_index):
    """Position in the q stream where the hint of the leaf at event index `leaf_index` would be read."""
    ev = unpack(ch['events'], np.uint8); q = unpack(ch['q'], np.uint16)
    pos = 0; i = 0
    while i < leaf_index:
        e = int(ev[i]); i += 1
        if e == 15: continue
        if e >= 64:
            qlo = int(q[pos]); qhi = 65536 - int(q[pos + 1]); pos += 2
            for present in (qlo > 0, qhi < 65536):
                if present:
                    cnt = int(q[pos]); pos += 1 + cnt
        elif e == 3:
            cnt = int(q[pos]); pos += 1 + cnt
    return pos


def test_perturbed_hint_keeps_verdict(m_certificate):
    """Gap 8: scaling or perturbing a multiplier hint may move a leaf to the fallback test, never flip it."""
    cert = json.loads(json.dumps(m_certificate[2])); code = m_certificate[3]
    for ch in cert['charts']:
        mu = unpack(ch['mu'], np.float32).copy()
        if len(mu):
            mu = (mu * 0.37 + 0.01).astype(np.float32); ch['mu'] = _pack(mu)
    r = check_certificate(cert, code=code)
    assert r['verified']


def test_v5_manifest_detects_missing_line(tmp_path, m_certificate):
    from scripts.completeness import check_certificates as cc
    base = tmp_path; (base / 'gen').mkdir(); (base / 'cert').mkdir(); (base / 'certificates').mkdir(); (base / 'checks').mkdir()
    codes = ['a.1', 'b.2', 'c.3']
    with gzip.open(base / 'gen' / 'gen-t-0-of-1.jsonl.gz', 'wt') as f:
        for c in codes: f.write(json.dumps(dict(code=c, faces=[])) + '\n')
    with open(base / 'cert' / 'cert-t-0-of-1.jsonl', 'w') as f:
        for c in codes: f.write(json.dumps(dict(code=c, status='excluded')) + '\n')
    with open(base / 'certificates' / 'certificates-t-0-of-1.jsonl', 'w') as f:
        for c in codes[:-1]: f.write(json.dumps(dict(code=c, certificate={})) + '\n')      # one line removed
    with open(base / 'checks' / ('%s-certificates-t-0-of-1.jsonl' % cc.CHECKER_VERSION), 'w') as f:
        for c in codes: f.write(json.dumps(dict(code=c, verified=True, terminal=[])) + '\n')
    m = cc.manifest('t', base=str(base))
    assert m['missing_from_certificates'] == 1 and not m['sets_equal_no_duplicates'] and not m['complete']
    with open(base / 'certificates' / 'certificates-t-0-of-1.jsonl', 'a') as f:
        f.write(json.dumps(dict(code=codes[-1], certificate={})) + '\n')
    m = cc.manifest('t', base=str(base))
    assert m['sets_equal_no_duplicates'] and m['complete']


# ---------------------------------------------------------------- V3: the model's domain (docs/lemmas.md L17)
COND_RHOMBUS = ('results/completeness/certificates/certificates-cond-21-of-40.jsonl', 1074)   # 0-based line
INF, NAN = float('inf'), float('nan')


def _honest_theta():
    """[theta_min rounded down, arctan 2 rounded up], the range real certificates write (decided in Arb)."""
    r = ck.arb(ck.R_MIN)
    return [ck.down(2 * (r / (1 + r)).asin()), ck.up(ck.arb(2).atan())]


@pytest.fixture(scope='module')
def rhombus_candidate():
    """A real cond candidate: 18 triangles and one rhombus, one chart (theta, A).  Its terminal box contains the
    icosahedron with one unstressed edge (A = 72 deg at theta = arctan 2), a census critical configuration."""
    path, ln = COND_RHOMBUS
    if not os.path.exists(path): pytest.skip('no cond certificates in the repository')
    with open(path) as f: d = json.loads(f.readlines()[ln])
    assert [len(f) for f in d['certificate']['faces']].count(4) == 1 and len(d['certificate']['charts']) == 1
    return d['code'], d['certificate']


def _one_leaf(cert0, theta):
    """A forged certificate: cert0's faces, the stated theta, and one chart whose starting box is the checker's
    own a-priori box at that theta (so V1 passes by construction) with a single excluded leaf (event 2), the
    real tol and no terminal box.  It claims that the candidate has no critical configuration."""
    from kisscrit.core import model as am
    lo, hi, _ = ck.apriori_box(am.Graph(cert0['faces'], {}), theta)
    lo = list(lo); hi = list(hi); lo[0], hi[0] = theta[0], theta[1]
    ch = dict(modes={}, lo0=lo, hi0=hi, events=_pack(np.array([2], np.uint8)), q=_pack(np.zeros(0, np.uint16)),
              mu=_pack(np.zeros(0, np.float32)), n_events=1, complete=True)
    return dict(version=2, theta=list(theta), tol=cert0['tol'], faces=cert0['faces'], complete=True, charts=[ch])


FORGED_THETA = {                              # Step 1 of V3-FIX: (a), (b), (c), (d)
    'a-minus-5': lambda h: [-5.0, h[1]],
    'b-up-to-7': lambda h: [1.0471975511965976, 7.0],
    'c-minus-inf': lambda h: [-INF, h[1]],
    'c-nan': lambda h: [NAN, h[1]],
    'd-plus-inf': lambda h: [h[0], INF],
    'd-nan': lambda h: [h[0], NAN],
}


def test_v3_one_leaf_rejected_at_honest_theta(rhombus_candidate):
    """Control for the forgeries below: at the honest range the one-leaf certificate is rejected (its box
    contains the census configuration), and V3 raises no domain problem."""
    code, cert0 = rhombus_candidate
    cert = _one_leaf(cert0, _honest_theta())
    assert ck.theta_domain_ok(cert['theta'])
    assert not any(p.startswith('V3') for p in ck.structural_check(cert, code, cert['faces']))
    assert _rejected(cert, code)


@pytest.mark.parametrize('which', sorted(FORGED_THETA))
def test_v3_theta_outside_model_domain(rhombus_candidate, which):
    """Gap V3 (audit B-1): a theta range outside [THETA_FLOOR, THETA_CEIL], or not finite, is a structural
    problem.  Before the fix, (a) and (b) were accepted with verified True and no terminal box."""
    code, cert0 = rhombus_candidate
    cert = _one_leaf(cert0, FORGED_THETA[which](_honest_theta()))
    probs = ck.structural_check(cert, code, cert['faces'])
    assert any(p.startswith('V3: theta range outside the model domain') for p in probs), probs
    assert _rejected(cert, code)


def test_v3_domain_keeps_honest_certificate(m_certificate):
    """The honest certificate of the fixture is still inside the domain and still verifies."""
    term, st, cert, code = m_certificate
    assert ck.theta_domain_ok(cert['theta']) and ck.theta_claim_ok(cert['theta'])
    r = check_certificate(cert, code=code, gen_faces=cert['faces'])
    assert r['verified'] and not r['structural']


def test_v3_domain_boundary():
    """theta_lo = THETA_FLOOR exactly passes V3 and the next double below fails; the same at THETA_CEIL; an
    empty or reversed range fails."""
    from kisscrit.core import model as am
    h = _honest_theta()
    assert ck.theta_claim_ok([am.THETA_FLOOR, h[1]])
    assert not ck.theta_claim_ok([float(np.nextafter(am.THETA_FLOOR, -INF)), h[1]])
    assert ck.theta_claim_ok([h[0], am.THETA_CEIL])
    assert not ck.theta_claim_ok([h[0], float(np.nextafter(am.THETA_CEIL, INF))])
    assert not ck.theta_domain_ok([1.2, 1.2]) and not ck.theta_domain_ok([1.3, 1.2])
