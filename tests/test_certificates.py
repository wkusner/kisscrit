"""Certificates: the search's proof tree replays exactly, the independent checker accepts an honest certificate
for the local maximum M, and rejects certificates in which the box around M is relabelled as excluded."""
import base64, json, zlib, collections
import numpy as np
import pytest
import kisscrit as kc
from kisscrit.completeness.structure import Structure, from_configuration
from kisscrit.completeness import runner
from kisscrit.core.certificate import replay, unpack
from kisscrit.core.checker import check_certificate

THETA = (np.radians(60.0033), np.radians(63.4349489))


@pytest.fixture(scope='module')
def m_certificate():
    D = json.load(open('data/n12_above1.json'))['types']
    t = D[1]                                            # M, r = 1.0010121
    X, c, s, nf = kc.solve_kkt(np.array(t['X']), [tuple(p) for p in t['pairs']])
    a = kc.analyze(X)
    sup = [tuple(p) for p, w in zip(a['pairs'], a['support']) if w]
    S = Structure(from_configuration(X, sup))
    term, st, cert = runner.run(S, *THETA, tol=1e-6, certificate=True)
    return term, st, json.loads(json.dumps(cert))


def test_replay_matches_search(m_certificate):
    term, st, cert = m_certificate
    cnt = collections.Counter(); tb = []
    for ch in cert['charts']:
        for kind, lo, hi, extra in replay(ch):
            cnt[kind] += 1
            if kind == 'terminal': tb.append((lo, hi))
    assert cnt['vsum'] == st['vsum'] and cnt['invalid'] == st['invalid'] and cnt['terminal'] == st['terminal']
    assert len(tb) == len(term) and all(np.array_equal(a[0], b[1]) and np.array_equal(a[1], b[2])
                                        for a, b in zip(tb, term))


def test_checker_accepts_honest_certificate(m_certificate):
    r = check_certificate(m_certificate[2])
    assert r['verified'] and not r['failed'] and len(r['terminal']) >= 1


@pytest.mark.parametrize('code', [1, 2, 5])
def test_checker_rejects_relabelled_terminal(m_certificate, code):
    cert = json.loads(json.dumps(m_certificate[2]))
    for ch in cert['charts']:
        ev = unpack(ch['events'], np.uint8).copy()
        idx = np.nonzero(ev == 7)[0]
        if len(idx):
            ev[idx[0]] = code
            ch['events'] = base64.b64encode(zlib.compress(ev.tobytes(), 9)).decode()
            break
    r = check_certificate(cert)
    assert not r['verified'] and len(r['failed']) >= 1
