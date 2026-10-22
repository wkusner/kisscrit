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


# ---------------------------------------------------------------- V1b: the start box's domain (docs/lemmas.md L31-L33)
# DOM-1 (2026-10-22).  The model's exact range of a rhombus corner from two corners is valid only for A in (0, pi]
# (L32); a start box outside it was accepted at 5da4bba (the DES-F2 forgery, reproduced there by the DOM-1 session).
# Numbers of the DES-F2 section 3 list: 4 (non-finite endpoints), 6 (the forged certificate), 7 (rhombus hi0 = 4.0),
# 8 (rhombus lo0 = -1.0), 9 (pentagon hi0 = 10.0).  Fixture: tests/dom1_forgery_fixture.py.
import glob, math
import dom1_forgery_fixture as F


@pytest.fixture(scope='module')
def r1_candidate():
    """The r1 saddle's candidate (E = 21, one chart, three pentagons, two rhombi) of the committed cond run."""
    base = F.load_base()
    if base is None: pytest.skip('no cond certificates in the repository')
    return base


def _kinds(cert0):
    """Coordinate indices of the chart 0 of cert0 by kind: theta, rhombus, pentagon (model.Graph's layout)."""
    from kisscrit.core import model as am
    G = am.Graph(cert0['faces'], cert0['charts'][0]['modes'])
    return dict(theta=[0], rhombus=[G.pidx[fi][0] for fi, k in enumerate(G.sizes) if k == 4],
                pentagon=[i for fi, k in enumerate(G.sizes) if k == 5 for i in G.pidx[fi]])


def _honest(cert0, events=None):
    """cert0 with honest theta and the checker's own a-priori start box (V1 and V3 pass by construction); with
    `events` (a uint8 array) the chart's proof tree is replaced, e.g. [2] = one excluded leaf."""
    cert = json.loads(json.dumps(cert0)); ch = cert['charts'][0]
    th = F.honest_theta(); lo, hi, _ = F.start_box(cert0, forged=False)
    cert['theta'] = list(th); ch['lo0'] = list(lo); ch['hi0'] = list(hi)
    if events is not None:
        ch.update(events=_pack(np.array(events, np.uint8)), q=_pack(np.zeros(0, np.uint16)),
                  mu=_pack(np.zeros(0, np.float32)), n_events=len(events))
    return cert


def _v1b(problems, j=None):
    return [p for p in problems if p.startswith('V1b') and (j is None or ' coordinate %d box' % j in p)]


def _rejected_any(cert, code):
    """Rejected: a ValueError or ArithmeticError (disjoint enclosures) is a rejection too, as in `_rejected`."""
    try:
        r = check_certificate(cert, code=code)
    except (ValueError, ArithmeticError):
        return True
    return not r['verified']


def test_v1b_honest_start_boxes_pass(r1_candidate):
    """Positive control: the checker's own a-priori start boxes, and the boxes of the committed line itself, raise
    no V1b problem; the rhombus and pentagon boundary doubles are decided as 'certainly <=' (pi enclosed at 128
    bits): hi0 = the double below pi passes, the next double up does not; the same for 2 pi on a pentagon."""
    code, cert0 = r1_candidate
    honest = _honest(cert0)
    assert ck.structural_check(honest, code, honest['faces']) == []
    assert not _v1b(ck.structural_check(cert0, code, cert0['faces']))          # the committed line: no V1b problem
    kinds = _kinds(cert0); jr, jp = kinds['rhombus'][0], kinds['pentagon'][0]
    for j, top in ((jr, math.pi), (jp, 2 * math.pi)):
        c = json.loads(json.dumps(honest)); c['charts'][0]['hi0'][j] = top
        assert not _v1b(ck.structural_check(c, code, c['faces']), j), (j, top)
        c['charts'][0]['hi0'][j] = float(np.nextafter(top, 10))
        assert _v1b(ck.structural_check(c, code, c['faces']), j), (j, top)      # the double above pi (2 pi)
    c = json.loads(json.dumps(honest)); c['charts'][0]['lo0'][jr] = 5e-324       # smallest positive double
    assert not _v1b(ck.structural_check(c, code, c['faces']), jr)
    c['charts'][0]['lo0'][jr] = 0.0                                              # 0 < lo0 is strict
    assert _v1b(ck.structural_check(c, code, c['faces']), jr)


@pytest.mark.parametrize('end', ['lo0', 'hi0'])
@pytest.mark.parametrize('bad', [NAN, INF, -INF], ids=['nan', 'plus-inf', 'minus-inf'])
@pytest.mark.parametrize('kind', ['theta', 'rhombus', 'pentagon'])
def test_v1b_non_finite_endpoints(r1_candidate, kind, bad, end):
    """DES-F2 item 4: NaN or +-inf in lo0 or hi0 of any coordinate is a V1b problem (before the fix, only some
    of them were a V1 problem, by accident of the comparisons)."""
    code, cert0 = r1_candidate
    j = _kinds(cert0)[kind][0]
    cert = _honest(cert0); cert['charts'][0][end][j] = bad
    probs = ck.structural_check(cert, code, cert['faces'])
    assert _v1b(probs, j), probs


def _one_leaf_r1(cert0, mutate):
    cert = _honest(cert0, events=[2]); mutate(cert['charts'][0]); return cert


@pytest.mark.parametrize('what', ['rhombus-hi0-4.0', 'rhombus-lo0-minus-1', 'pentagon-hi0-10'])
def test_v1b_domain_violations_on_one_leaf(r1_candidate, what):
    """DES-F2 items 7 (rhombus hi0 = 4.0), 8 (rhombus lo0 = -1.0; before the fix the replay raised an
    ArithmeticError, a rejection by accident), 9 (pentagon hi0 = 10.0, hygiene: before the fix accepted by V1 and
    harmless).  The structural check names the coordinate; the certificate is not verified.  Note: check_certificate
    still replays the tree after a structural problem, so it may raise ArithmeticError on the rhombus boxes."""
    code, cert0 = r1_candidate
    kinds = _kinds(cert0); jr, jp = kinds['rhombus'][0], kinds['pentagon'][0]
    j, end, val = {'rhombus-hi0-4.0': (jr, 'hi0', 4.0), 'rhombus-lo0-minus-1': (jr, 'lo0', -1.0),
                   'pentagon-hi0-10': (jp, 'hi0', 10.0)}[what]
    cert = _one_leaf_r1(cert0, lambda ch: ch[end].__setitem__(j, val))
    probs = ck.structural_check(cert, code, cert['faces'])
    assert len(_v1b(probs, j)) == 1 and 'docs/lemmas.md L32' in _v1b(probs, j)[0], probs
    assert not [p for p in probs if p.startswith('V1:')], probs                   # V1 itself passes: only V1b speaks
    assert _rejected_any(cert, code)


def test_v1b_forged_rhombus_certificate_rejected(r1_candidate):
    """DES-F2 item 6: the forged certificate (rhombus start box [alpha(theta_lo), 2 pi + alpha(theta_lo)], 1,791
    events, every leaf closed by the checker's own direct tests).  At 5da4bba it was accepted: verified True, no
    structural problem, 0 failed, 0 open, 0 terminal.  It is rejected now, by V1b, and only by V1b: the model is
    unchanged, so without the structural check the replay still accepts the tree (the first assertion pins that,
    and must be revisited if model.py's rhombus range is ever made sound on any box)."""
    code, cert, j = F.forged_certificate()
    assert j == 1 and cert['charts'][0]['hi0'][j] == 2 * math.pi + cert['charts'][0]['lo0'][j]
    assert cert['charts'][0]['n_events'] == F.N_EVENTS
    r0 = check_certificate(cert, code=code, gen_faces=cert['faces'], structural=False)
    assert r0['verified'] and not r0['failed'] and r0['open'] == 0 and not r0['terminal'], r0['counts']
    probs = ck.structural_check(cert, code, cert['faces'])
    assert len(probs) == 1 and probs[0].startswith('V1b: chart') and ' coordinate 1 box [' in probs[0], probs
    r = check_certificate(cert, code=code, gen_faces=cert['faces'])
    assert not r['verified'] and r['structural'] == probs and r['terminal'] == []


@pytest.mark.skipif(not __import__('kisscrit.core.carith', fromlist=['x']).available(), reason='C kernel not built')
def test_v1b_forged_certificate_rejected_by_kernel_front_end(r1_candidate):
    """The C kernel front end calls the same structural_check (kisscrit/ckernel.py), so it rejects the forgery too."""
    from kisscrit import ckernel
    code, cert, j = F.forged_certificate()
    r = ckernel.check_certificate(cert, code=code, gen_faces=cert['faces'])
    assert not r['verified'] and len(_v1b(r['structural'], 1)) == 1, r['structural']


@pytest.mark.skipif(os.environ.get('DOM1_REGEN') != '1', reason='set DOM1_REGEN=1 (about 35 s)')
def test_forgery_fixture_regenerates(r1_candidate):
    """The stored event stream is what the checker-as-oracle forger produces (deterministic)."""
    ok, ev = F.regenerate_events()
    assert ok and np.array_equal(ev, np.frombuffer(zlib.decompress(base64.b64decode(F.EVENTS_B64)), np.uint8))


def _scan_cond(files):
    """Read-only scan of cond certificate files: counts and extremes of every chart's start box by coordinate
    kind, and the number of coordinates violating the V1b rule (decided as the checker does, in Arb at 128 bits)."""
    from flint import arb, ctx
    from kisscrit.core import model as am
    with ctx.workprec(128): PI = arb.pi(); TPI = 2 * PI
    S = dict(lines=0, charts=0, coords=0, nonfinite=0, rh=0, pe=0, bad=0, rh_lo=math.inf, rh_hi=-math.inf)
    first = None
    for fn in files:
        with open(fn) as fh:
            for line in fh:
                d = json.loads(line); c = d['certificate']; S['lines'] += 1
                for ci, ch in enumerate(c['charts']):
                    if ch.get('not_run'): continue
                    G = am.Graph(c['faces'], ch['modes']); S['charts'] += 1
                    rh = {G.pidx[fi][0] for fi, k in enumerate(G.sizes) if k == 4}
                    pe = {i for fi, k in enumerate(G.sizes) if k == 5 for i in G.pidx[fi]}
                    assert G.npar == len(ch['lo0']) == len(ch['hi0'])
                    if first is None and rh: first = (d['code'], c, ci, sorted(rh)[0])
                    for j in range(G.npar):
                        a, b = float(ch['lo0'][j]), float(ch['hi0'][j]); S['coords'] += 1
                        if not (math.isfinite(a) and math.isfinite(b)): S['nonfinite'] += 1; S['bad'] += 1; continue
                        if j in rh:
                            S['rh'] += 1; S['rh_lo'] = min(S['rh_lo'], a); S['rh_hi'] = max(S['rh_hi'], b)
                            S['bad'] += not (a > 0 and bool(arb(b) <= PI))
                        elif j in pe:
                            S['pe'] += 1; S['bad'] += not (a >= 0 and bool(arb(b) <= TPI))
    return S, first


def test_v1b_cond_certificates_start_boxes_inside_domain():
    """Every rhombus start box of the scanned cond certificates lies in (0, pi], every endpoint is finite, and
    pentagon boxes lie in [0, 2 pi].  The suite scans every fifth file (8 of 40, about 8 s); DOM1_FULL_SCAN=1 scans
    all 40 (45,260 lines, about 35 s; result in docs/session-log.md, 2026-10-22).  The last assertions feed one
    scanned certificate to the checker with a rhombus box mutated to hi0 = 2 pi + lo0."""
    allf = sorted(glob.glob('results/completeness/certificates/certificates-cond-*-of-40.jsonl'),
                  key=lambda p: int(p.split('cond-')[1].split('-of')[0]))
    if not allf: pytest.skip('no cond certificates in the repository')
    files = allf if os.environ.get('DOM1_FULL_SCAN') == '1' else allf[::5]
    S, first = _scan_cond(files)
    assert S['lines'] > 5000 and S['rh'] > 10000 and S['pe'] > 10000, S
    assert S['nonfinite'] == 0 and S['bad'] == 0, S
    assert 0 < S['rh_lo'] and S['rh_hi'] <= math.pi, S
    if os.environ.get('DOM1_FULL_SCAN') == '1': assert S['lines'] == 45260, S
    code, cert, ci, j = first
    cert = json.loads(json.dumps(cert)); ch = cert['charts'][ci]
    assert not _v1b(ck.structural_check(cert, code, cert['faces']))
    ch['hi0'][j] = 2 * math.pi + ch['lo0'][j]
    assert _v1b(ck.structural_check(cert, code, cert['faces']), j)
