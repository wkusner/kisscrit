"""Exact arithmetic in a real number field K = Q[t]/(P) with a chosen real embedding (EXACT-1).

Status: untrusted (exploration or tooling; no part of the completeness proof).

Every decision made here is exact:
  - P is an integer polynomial; irreducibility is checked by factoring over Z (FLINT).
  - The real root t* of P that defines the embedding is isolated in a rational interval [lo, hi] whose
    Sturm count is 1, with P(lo) and P(hi) nonzero; both facts are checked in exact rational arithmetic.
  - An element a of K is zero iff its reduction mod P is the zero polynomial (exact).
  - The sign of a nonzero element a at t* is decided by evaluating a in ball arithmetic (Arb) on the
    isolating interval, bisecting the interval exactly (by the sign of P at the midpoint) until the ball
    excludes zero.  Since a(t*) != 0 and the interval shrinks to t*, this terminates, and the answer is the
    sign of a(t*): a ball that excludes 0 has one sign on the whole interval, which contains t*.
  - Linear algebra (rank, null space) and the inertia of a symmetric matrix are done by Gaussian elimination
    over K with these exact zero and sign tests.

Floating point never decides anything here; mpmath values are used only as approximations to choose an
interval, and the interval is then checked exactly.
"""
from flint import fmpq_poly, fmpz_poly, fmpq, fmpz, arb, ctx


def _q(x):
    """A rational number as fmpq (from int, fmpq, fractions.Fraction or 'p/q' string)."""
    if isinstance(x, fmpq):
        return x
    if isinstance(x, str):
        if '/' in x:
            a, b = x.split('/')
            return fmpq(int(a), int(b))
        return fmpq(int(x))
    if hasattr(x, 'numerator'):
        return fmpq(int(x.numerator), int(x.denominator))
    return fmpq(x)


def qstr(x):
    """Canonical string of a rational: 'p' or 'p/q'."""
    x = _q(x)
    return str(x.p) if x.q == 1 else '%d/%d' % (int(x.p), int(x.q))


def sturm_count(P, lo, hi):
    """Number of distinct real roots of the squarefree rational polynomial P in (lo, hi], exactly."""
    seq = [P, P.derivative()]
    while seq[-1].degree() > 0:
        r = -(seq[-2] % seq[-1])
        if r.is_zero():
            break
        seq.append(r)

    def var(x):
        s = [p(x) for p in seq]
        s = [v for v in s if v != 0]
        return sum(1 for a, b in zip(s, s[1:]) if (a > 0) != (b > 0))
    return var(lo) - var(hi)


class Field:
    """K = Q(t*), t* the unique root of the irreducible integer polynomial P in the rational interval [lo, hi]."""

    def __init__(self, P, approx=None, lo=None, hi=None):
        """P: integer coefficients, constant term first.  Give either an mpmath approximation `approx` of the
        root (an interval is then chosen and checked exactly) or rational endpoints lo, hi (checked exactly)."""
        P = [int(a) for a in P]
        while P and P[-1] == 0:
            P.pop()
        self.Pz = fmpz_poly(P)
        self.P = fmpq_poly(P)
        self.d = self.Pz.degree()
        content, fac = self.Pz.factor()
        assert len(fac) == 1 and fac[0][1] == 1 and fac[0][0].degree() == self.d, 'P is not irreducible over Q'
        if lo is None:
            lo, hi = self._isolate(approx)
        self.lo, self.hi = _q(lo), _q(hi)
        assert self.lo < self.hi
        assert self.P(self.lo) != 0 and self.P(self.hi) != 0, 'endpoint is a root'
        assert sturm_count(self.P, self.lo, self.hi) == 1, 'interval does not isolate exactly one root'
        self.root_q = (-self.P[0] / self.P[1]) if self.d == 1 else None

    def _isolate(self, approx):
        """Choose rational endpoints around the approximation and shrink until the Sturm count is 1."""
        import mpmath as mp
        x = mp.mpf(approx)
        for k in range(8, 4000, 8):
            eps = mp.mpf(10) ** (-k)
            lo = self._rat(x - eps, k + 10)
            hi = self._rat(x + eps, k + 10)
            if self.P(lo) == 0 or self.P(hi) == 0:
                continue
            n = sturm_count(self.P, lo, hi)
            if n == 1:
                return lo, hi
            if n == 0:
                raise AssertionError('approximation is not near a root of P')
        raise AssertionError('could not isolate the root')

    @staticmethod
    def _rat(x, digits):
        """A rational approximation of the mpmath number x with `digits` decimal digits."""
        import mpmath as mp
        with mp.workdps(digits + 10):
            n = int(mp.nint(x * mp.mpf(10) ** digits))
        return fmpq(n, 10 ** digits)

    # ----- elements -----------------------------------------------------------------------------------
    def red(self, a):
        """Reduce a rational polynomial mod P."""
        if self.d > 1:
            return a % self.P
        return fmpq_poly([a(self.root_q)]) if not a.is_zero() else fmpq_poly([])

    def __call__(self, coeffs):
        """The element sum coeffs[k] t^k (coefficients rational, constant first)."""
        return KE(self, self.red(fmpq_poly([_q(c) for c in coeffs])))

    def const(self, x):
        """The rational constant x as an element of K."""
        return KE(self, fmpq_poly([_q(x)]))

    def gen(self):
        """The generator t (the contact cosine c)."""
        return KE(self, self.red(fmpq_poly([0, 1])))

    def refine(self):
        """Halve the isolating interval exactly."""
        m = (self.lo + self.hi) / 2
        pm = self.P(m)
        if pm == 0:  # impossible for irreducible P of degree > 1 unless d == 1
            self.lo = self.hi = m
            return
        if (self.P(self.lo) > 0) != (pm > 0):
            self.hi = m
        else:
            self.lo = m

    def ball(self):
        """An Arb ball containing the isolating interval (rigorous: arb of a rational is an enclosure)."""
        return arb(self.lo).union(arb(self.hi))

    def sign(self, a):
        """Exact sign (+1, 0, -1) of the element a at t*."""
        p = a.p if isinstance(a, KE) else a
        if p.is_zero():
            return 0
        if self.d == 1:
            v = p(self.root_q)
            return 1 if v > 0 else -1
        old = ctx.prec
        try:
            prec = 128
            while True:
                ctx.prec = prec
                x = self.ball()
                v = arb(0)
                for k in range(p.degree(), -1, -1):
                    v = v * x + arb(p[k])
                if v > 0:
                    return 1
                if v < 0:
                    return -1
                for _ in range(16):
                    self.refine()
                prec = min(prec * 2, 1 << 20) if prec < 4096 else prec + 1024
        finally:
            ctx.prec = old

    def eval_arb(self, a, prec=256):
        """An Arb enclosure of a(t*) (for reporting only)."""
        p = a.p if isinstance(a, KE) else a
        old = ctx.prec
        try:
            ctx.prec = prec
            while self.hi - self.lo > fmpq(1, 2 ** (prec // 2)):
                self.refine()
            x = self.ball()
            v = arb(0)
            for k in range(p.degree(), -1, -1):
                v = v * x + arb(p[k])
            return v
        finally:
            ctx.prec = old

    def coeffs(self, a):
        """Coefficient strings of a (constant term first, length d)."""
        p = a.p if isinstance(a, KE) else a
        c = [p[k] for k in range(self.d)]
        return [qstr(x) for x in c]


class KE:
    """An element of K, stored as a rational polynomial of degree < d reduced mod P."""
    __slots__ = ('F', 'p')

    def __init__(self, F, p):
        self.F, self.p = F, p

    def _c(self, o):
        if isinstance(o, KE):
            return o.p
        return fmpq_poly([_q(o)])

    def __add__(self, o):
        return KE(self.F, self.p + self._c(o))
    __radd__ = __add__

    def __sub__(self, o):
        return KE(self.F, self.p - self._c(o))

    def __rsub__(self, o):
        return KE(self.F, self._c(o) - self.p)

    def __neg__(self):
        return KE(self.F, -self.p)

    def __mul__(self, o):
        return KE(self.F, self.F.red(self.p * self._c(o)))
    __rmul__ = __mul__

    def inv(self):
        """Multiplicative inverse (exact, by the extended Euclidean algorithm with P)."""
        assert not self.p.is_zero(), 'division by zero in K'
        if self.F.d == 1:
            return KE(self.F, fmpq_poly([1 / self.F.red(self.p)[0]]))
        g, s, _ = self.p.xgcd(self.F.P)
        assert g.degree() == 0
        return KE(self.F, self.F.red(s / g[0]))

    def __truediv__(self, o):
        return self * (o.inv() if isinstance(o, KE) else KE(self.F, fmpq_poly([1 / _q(o)])))

    def is_zero(self):
        return self.p.is_zero()

    def __eq__(self, o):
        return (self - o).p.is_zero()

    def __hash__(self):
        return hash(str(self.p))

    def sign(self):
        return self.F.sign(self)

    def __repr__(self):
        return 'KE(%s)' % self.p


# ----- linear algebra over K ------------------------------------------------------------------------------

def rref(M):
    """Reduced row echelon form of a matrix over K (list of lists of KE); returns (R, pivot columns)."""
    R = [row[:] for row in M]
    m = len(R)
    n = len(R[0]) if m else 0
    piv = []
    r = 0
    for col in range(n):
        p = next((i for i in range(r, m) if not R[i][col].is_zero()), None)
        if p is None:
            continue
        R[r], R[p] = R[p], R[r]
        inv = R[r][col].inv()
        R[r] = [x * inv for x in R[r]]
        for i in range(m):
            if i != r and not R[i][col].is_zero():
                f = R[i][col]
                R[i] = [a - f * b for a, b in zip(R[i], R[r])]
        piv.append(col)
        r += 1
        if r == m:
            break
    return R, piv


def rank(M):
    """Rank of a matrix over K."""
    return len(rref(M)[1]) if M else 0


def nullspace(M, n=None):
    """A basis of {x : M x = 0} over K (columns of M are the unknowns)."""
    F = M[0][0].F if M else None
    n = n if n is not None else len(M[0])
    R, piv = rref(M)
    free = [j for j in range(n) if j not in piv]
    basis = []
    for f in free:
        v = [F.const(0) for _ in range(n)]
        v[f] = F.const(1)
        for i, pc in enumerate(piv):
            v[pc] = -R[i][f]
        basis.append(v)
    return basis


def matmul(A, B):
    """Product of matrices over K."""
    return [[sum((A[i][k] * B[k][j] for k in range(len(B))), A[0][0].F.const(0)) for j in range(len(B[0]))]
            for i in range(len(A))]


def det3(M):
    """Determinant of a 3x3 matrix over K."""
    return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
            + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))


def inv3(M):
    """Inverse of an invertible 3x3 matrix over K (adjugate over determinant)."""
    d = det3(M)
    assert not d.is_zero()
    di = d.inv()
    c = lambda i, j: (M[(i + 1) % 3][(j + 1) % 3] * M[(i + 2) % 3][(j + 2) % 3]
                      - M[(i + 1) % 3][(j + 2) % 3] * M[(i + 2) % 3][(j + 1) % 3])
    return [[c(j, i) * di for j in range(3)] for i in range(3)]


def inertia(S):
    """(n_plus, n_minus, n_zero) of a symmetric matrix over the ordered field K, exactly (Sylvester's law).

    Congruence diagonalization: a nonzero diagonal pivot is eliminated symmetrically; if every diagonal entry
    is zero but some off-diagonal a_kl is not, row/column l is added to k (a congruence), which makes the new
    diagonal entry 2 a_kl nonzero."""
    A = [row[:] for row in S]
    pos = neg = 0
    while A:
        n = len(A)
        k = next((i for i in range(n) if not A[i][i].is_zero()), None)
        if k is None:
            kl = next(((i, j) for i in range(n) for j in range(n) if i != j and not A[i][j].is_zero()), None)
            if kl is None:
                return pos, neg, len(A) + 0
            i, j = kl
            A[i] = [a + b for a, b in zip(A[i], A[j])]
            for r in range(n):
                A[r][i] = A[r][i] + A[r][j]
            k = i
        piv = A[k][k]
        s = piv.sign()
        pos += s > 0
        neg += s < 0
        inv = piv.inv()
        rest = [i for i in range(n) if i != k]
        A = [[A[i][j] - A[i][k] * A[k][j] * inv for j in rest] for i in rest]
    return pos, neg, 0
