"""Shared helpers of the exploration driver (worker, merge, export).

Status: untrusted (exploration or tooling; no part of the proof).

Everything here is numerical (double precision).  Nothing computed by the exploration driver is verified, and no
value it writes is an input to a proof.  This module holds what the three scripts share: atomic file writing,
reading a census file, the numerical quantities recomputed from X and the contact pairs, and the thresholds
that turn those quantities into the flags of the refine queue.
"""
import json, math, os, sys, tempfile
from datetime import datetime, timezone
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from kisscrit.search import Census
from kisscrit.stress import balanced_stress

# Thresholds of the refine queue.  They are stated in every output header.  Recomputed with kkt_metrics, the 12
# types of data/n12_above1.json have a KKT residual of at most 1.1e-11 and a smallest non-contact gap of 1.5e-4
# (double precision).  kisscrit/kkt.py accepts a solution at residual 1e-10 and gap 1e-8, and kisscrit/stress.py
# calls a stress strictly positive above 1e-10.
KKT_RESID_MAX = 1e-9    # above this the record is flagged suspect
GAP_MIN = 1e-6          # a non-contact cosine closer than this to the contact cosine makes it suspect
STRESS_MIN = 1e-10      # a max-min stress above this counts as strictly positive

# The fields of a record of data/n12_above1.json, in that file's order (the layout the refine scripts load).
REFINE_FIELDS = ('r', 'c', 'k', 'coindex', 'kernel', 'sdim', 'tmin', 'sym', 'sym_proper', 'pairs', 'sigma', 'X')

HEADER_LABEL = ('numerical (double precision); nothing in this file is verified and no value in it is a proof input; '
                'thresholds: KKT residual > %g or minimum non-contact gap < %g = suspect, max-min stress <= %g = '
                'not strictly positive' % (KKT_RESID_MAX, GAP_MIN, STRESS_MIN))


def utc_now():
    """Current UTC time as an ISO string with seconds; it sorts like the time it names."""
    return datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')


def atomic_write_text(text, path):
    """Write text to a temporary name in the same directory, flush it to disk, then rename it to path.

    A reader (or a crash) sees either the previous complete file or the new complete file, never half of one."""
    d = os.path.dirname(os.path.abspath(path))
    fd, tmp = tempfile.mkstemp(prefix=os.path.basename(path) + '.tmp.', dir=d)
    try:
        with os.fdopen(fd, 'w') as f:
            f.write(text)
            f.flush(); os.fsync(f.fileno())
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise


def atomic_write_json(obj, path):
    """atomic_write_text of the JSON text of obj."""
    atomic_write_text(json.dumps(obj), path)


def load_census_file(path):
    """Census.load of kisscrit.search, with every way a damaged file can fail turned into one ValueError.

    Returns (census, metadata dict of the file).  A truncated file fails here, so a caller can skip it."""
    try:
        with open(path) as f:
            D = json.load(f)
        C = Census.load(path)
        return C, {k: v for k, v in D.items() if k not in ('N', 'd', 'types')}
    except (OSError, ValueError, KeyError, TypeError, IndexError) as e:
        raise ValueError('%s: %s' % (os.path.basename(path), e))


def kkt_metrics(rec):
    """Numerical KKT quantities of one record, recomputed in double precision from X and the contact pairs only.

    Returns a dict with
      sigma       a nonnegative balanced stress summing to 1 that maximizes its smallest entry (kisscrit.stress.
                  balanced_stress, an LP), or None if no such stress exists;
      stress_min  that smallest entry (the max-min value; > STRESS_MIN means strictly positive);
      kkt_resid   the largest of three residuals: | |x_i| - 1 |, | <x_i, x_j> - c | over contacts, and the norm of the
                  tangential balance sum_j sigma_ij (x_j - <x_i, x_j> x_i) over vertices i;
      min_gap     c minus the largest non-contact cosine (None if every pair is a contact);
      n_contacts  the number of contact pairs."""
    X = np.asarray(rec['X'], float); pairs = [tuple(p) for p in rec['pairs']]; c = float(rec['c'])
    G = X @ X.T
    st = balanced_stress(X, pairs)
    sigma = st['sigma']
    resid = float(np.abs(np.linalg.norm(X, axis=1) - 1).max())
    resid = max(resid, max(abs(G[i, j] - c) for i, j in pairs))
    if sigma is None:
        resid = math.inf
    else:
        bal = np.zeros_like(X)
        for s, (i, j) in zip(sigma, pairs):
            bal[i] += s * (X[j] - G[i, j] * X[i]); bal[j] += s * (X[i] - G[i, j] * X[j])
        resid = max(resid, float(np.linalg.norm(bal, axis=1).max()))
    contact = set(pairs)
    others = [G[i, j] for i in range(len(X)) for j in range(i + 1, len(X)) if (i, j) not in contact]
    return dict(sigma=None if sigma is None else [float(s) for s in sigma], stress_min=float(st['tmin']),
                kkt_resid=resid, min_gap=float(c - max(others)) if others else None, n_contacts=len(pairs))


def classify(rec, m):
    """Flag of the refine queue from a record and its kkt_metrics m.  The first rule that applies wins:
    suspect (residual above KKT_RESID_MAX, or gap below GAP_MIN, or no balanced stress at all), rattler (the
    record has a vertex outside the stressed graph), regular (kernel 0 and max-min stress above STRESS_MIN),
    else degenerate.  'regular' is the only flag that asks for a solve in the refine stage."""
    if not m['kkt_resid'] <= KKT_RESID_MAX or (m['min_gap'] is not None and m['min_gap'] < GAP_MIN):
        return 'suspect'
    if rec.get('rattlers', 0) > 0:
        return 'rattler'
    if rec['kernel'] == 0 and m['stress_min'] > STRESS_MIN:
        return 'regular'
    return 'degenerate'
