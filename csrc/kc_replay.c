/* kc_replay.c: the replay of one chart's proof tree for the C checker kernel (contract in kc_replay.h).
 *
 * Status: trusted.
 * Implements: checklist F3b; docs/certificate-format.md section 3 (contraction rounds, grid values, multiplier
 * hints, removed slabs, leaf events, midpoint splits, order of the children).  A bit-for-bit reimplementation of
 * kisscrit/core/certificate.py:replay: the same reads of the three streams in the same order, the same double
 * operations in the same order, the same comparisons, and the same points of failure.  Comments of the form
 * "certificate.py: ..." quote the reference line that the code below them mirrors.
 *
 * Arithmetic: only the two box operations of the format, each a fixed sequence of IEEE double operations
 * (correctly rounded, round to nearest; kc_init asserts that mode and the absence of FMA contraction):
 *   grid value   t = q / 65536.0;  d = hi0 - lo0;  p = d * t;  g = lo0 + p   (four roundings, in this order);
 *   midpoint     s = lo + hi;  m = 0.5 * s                                    (two roundings, in this order).
 * Each is written as separate statements so that no reassociation is possible, and the build forbids FMA
 * contraction (-ffp-contract=off; with contraction the grid value would be fma(d, t, lo0), one rounding).
 * Comparisons follow IEEE semantics, as numpy's: every comparison with a NaN is false.
 * Build with -O2 -fno-fast-math -ffp-contract=off (csrc/Makefile).
 */
#include <string.h>
#include "kc_replay.h"

enum { PHASE_POP = 0, PHASE_NODE = 1, PHASE_END = 2 };

/* ===================================================================== the two box operations */

static double grid_value(double lo0, double hi0, int64_t q)
{   /* certificate.py grid_value: q <= 0 -> lo0; q >= QMAX -> hi0; else lo0 + (hi0 - lo0) * (q / 65536.0).
     * Claim: the same double as the reference for every q and every pair of doubles.  The reference
     * evaluates (q / 65536.0) first (an int divided by a float: q converted exactly, since q < 2^16, then one
     * correctly rounded division, which is exact as 65536 is a power of two), then hi0 - lo0 (numpy float64
     * subtraction), their product, and finally the sum; each step is one IEEE operation. */
    if (q <= 0) return lo0;
    if (q >= KC_QMAX) return hi0;
    double t = (double)q / 65536.0;
    double d = hi0 - lo0;
    double p = d * t;
    double g = lo0 + p;
    return g;
}

static double midpoint(double lo, double hi)
{   /* certificate.py: mid = 0.5 * (lo[k] + hi[k]).  Claim: the same double: one rounded sum, then one rounded
     * product by 0.5 (not lo + 0.5 * (hi - lo), which differs in some bits). */
    double s = lo + hi;
    double m = 0.5 * s;
    return m;
}

/* ===================================================================== stream reads */

static int read_hint(kc_replay *r, kc_hint *h)
{   /* certificate.py read_mu:
     *     cnt = int(q[pos['q']]); pos['q'] += 1
     *     rows = q[pos['q']:pos['q'] + cnt]; pos['q'] += cnt
     *     vals = mus[pos['m']:pos['m'] + cnt]; pos['m'] += cnt
     *     return dict(zip(rows, vals))
     * Claim: the same positions afterwards and the same (rows, vals) pairs.  q[pos] past the end raises
     * IndexError; the slices are cut at the end of their arrays without error, and zip stops at the shorter,
     * so count = min(cnt, rows available, values available). */
    if (r->qpos >= r->nq) return KC_RERR_INDEX;
    int64_t cnt = r->q[r->qpos]; r->qpos += 1;
    int64_t rows_avail = r->nq - r->qpos;                 /* >= 0 here, since qpos <= nq after the read */
    if (rows_avail > cnt) rows_avail = cnt;
    int64_t vals_avail = r->nmu - r->mpos;                /* may be negative: mpos can pass nmu */
    if (vals_avail < 0) vals_avail = 0;
    if (vals_avail > cnt) vals_avail = cnt;
    h->present = 1;
    h->count = rows_avail < vals_avail ? rows_avail : vals_avail;
    h->rows = h->count > 0 ? r->q + r->qpos : NULL;
    h->vals = h->count > 0 ? r->mu + r->mpos : NULL;
    r->qpos += cnt;
    r->mpos += cnt;
    return 0;
}

static kc_hint no_hint(void) { kc_hint h = { 0, 0, NULL, NULL }; return h; }

/* ===================================================================== sizes and set-up */

int64_t kc_replay_stack_need(const uint8_t *ev, int64_t nev)
{   /* Claim: the reference's stack never holds more entries.  It starts with one entry; each node pops one,
     * and a split pushes two, so its length is 1 + (splits so far) - (leaves so far) <= 1 + (split events). */
    int64_t s = 1;
    for (int64_t i = 0; i < nev; i++)
        if (ev[i] >= KC_EV_SPLIT && ev[i] < KC_EV_CONTRACT) s++;
    return s;
}

int64_t kc_replay_work_doubles(int n, int64_t stack_cap)
{   /* the layout of kc_replay_init: stack, then lo, hi, base_lo, base_hi, cur_lo, cur_hi, then the slabs */
    return 2 * (int64_t)n * stack_cap + 6 * (int64_t)n + 8 * (int64_t)n * n;
}

int64_t kc_replay_sizeof(void) { return (int64_t)sizeof(kc_replay); }

int kc_replay_init(kc_replay *r, int n, const double *lo0, const double *hi0,
                   const uint8_t *ev, int64_t nev, const uint16_t *q, int64_t nq,
                   const float *mu, int64_t nmu, int has_mu, double *work, int64_t stack_cap)
{   /* certificate.py: stack = [(lo0, hi0)]; i = 0; pos = {'q': 0, 'm': 0}; mus empty without a mu stream */
    memset(r, 0, sizeof *r);
    if (n < 0 || n > KC_REPLAY_MAXN || stack_cap < 1) { r->phase = PHASE_END; r->err = KC_RERR_LIMIT; return KC_RERR_LIMIT; }
    r->n = n; r->has_mu = has_mu ? 1 : 0;
    r->ev = ev; r->nev = nev; r->i = 0;
    r->q = q; r->nq = nq; r->qpos = 0;
    r->mu = r->has_mu ? mu : NULL; r->nmu = r->has_mu ? nmu : 0; r->mpos = 0;
    r->stack = work; r->stack_cap = stack_cap;
    double *w = work + 2 * (int64_t)n * stack_cap;
    r->lo = w; r->hi = w + n; r->base_lo = w + 2 * n; r->base_hi = w + 3 * n; r->cur_lo = w + 4 * n; r->cur_hi = w + 5 * n;
    r->slab_store = w + 6 * n;
    memcpy(r->stack, lo0, sizeof(double) * n);
    memcpy(r->stack + n, hi0, sizeof(double) * n);
    r->stack_len = 1;
    r->phase = PHASE_POP;
    return 0;
}

/* ===================================================================== one step */

static int finish(kc_replay *r, int code) { r->phase = PHASE_END; r->err = code; return code; }

static int push(kc_replay *r, const double *lo, const double *hi)
{   /* certificate.py: stack.append((lo, hi)) */
    if (r->stack_len >= r->stack_cap) return KC_RERR_STACK;
    double *e = r->stack + 2 * (int64_t)r->n * r->stack_len;
    memcpy(e, lo, sizeof(double) * r->n);
    memcpy(e + r->n, hi, sizeof(double) * r->n);
    r->stack_len++;
    return 0;
}

static int contraction_round(kc_replay *r, kc_node *out)
{   /* One round: the event 15 at r->i, the events 64 + k after it, then the removed slabs.  certificate.py:
     *     i += 1; base_lo, base_hi = lo.copy(), hi.copy(); nchg = 0; hints = {}
     *     while i < len(ev) and ev[i] >= EV_CONTRACT: ...
     *     if nchg == 0: raise ValueError('empty contraction round')
     * Claim: the same new box, the same slabs and hints, and the same errors at the same reads. */
    const int n = r->n;
    r->i += 1;
    memcpy(r->base_lo, r->lo, sizeof(double) * n);
    memcpy(r->base_hi, r->hi, sizeof(double) * n);
    for (int k = 0; k < n; k++) { r->hint_lower[k] = no_hint(); r->hint_upper[k] = no_hint(); }
    int64_t nchg = 0;
    while (r->i < r->nev && r->ev[r->i] >= KC_EV_CONTRACT) {
        int k = (int)r->ev[r->i] - KC_EV_CONTRACT; r->i += 1; nchg += 1;
        /* qlo = int(q[pos['q']]); qhi = QMAX - int(q[pos['q'] + 1]); pos['q'] += 2   (IndexError past the end) */
        if (r->qpos + 1 >= r->nq) return KC_RERR_INDEX;
        int64_t qlo = r->q[r->qpos];
        int64_t qhi = KC_QMAX - (int64_t)r->q[r->qpos + 1];
        r->qpos += 2;
        /* if qlo > 0: hints[(k, -1)] = read_mu();  if qhi < QMAX: hints[(k, +1)] = read_mu()
         * (a later event for the same k replaces the hint, as the dict assignment does) */
        kc_hint hl, hu; int rc;
        int have_lower = qlo > 0, have_upper = qhi < KC_QMAX;
        if (have_lower && (rc = read_hint(r, &hl)) != 0) return rc;
        if (have_upper && (rc = read_hint(r, &hu)) != 0) return rc;
        /* grid_value(base_lo[k], ...): IndexError for k >= n, raised after the hint reads above */
        if (k >= n) return KC_RERR_INDEX;
        if (have_lower) r->hint_lower[k] = hl;
        if (have_upper) r->hint_upper[k] = hu;
        double nlo = grid_value(r->base_lo[k], r->base_hi[k], qlo);
        double nhi = grid_value(r->base_lo[k], r->base_hi[k], qhi);
        /* if not (base_lo[k] <= nlo <= nhi <= base_hi[k]): raise ValueError('contraction outside the box') */
        if (!(r->base_lo[k] <= nlo && nlo <= nhi && nhi <= r->base_hi[k])) return KC_RERR_OUTSIDE;
        r->lo[k] = nlo; r->hi[k] = nhi;
    }
    if (nchg == 0) return KC_RERR_EMPTY_ROUND;

    /* The removed region, slabs one coordinate at a time in increasing k, lower before upper.  certificate.py:
     *     cur_lo, cur_hi = base_lo.copy(), base_hi.copy()
     *     for k in range(n):
     *         if cur_lo[k] < lo[k]:  slab = cur with s_hi[k] = lo[k]; cut from cur; side -1; cur_lo[k] = lo[k]
     *         if hi[k] < cur_hi[k]:  slab = cur with s_lo[k] = hi[k]; cut from cur; side +1; cur_hi[k] = hi[k]
     * Claim: the slabs are the same boxes in the same order (the second test sees the updated cur_lo[k]). */
    memcpy(r->cur_lo, r->base_lo, sizeof(double) * n);
    memcpy(r->cur_hi, r->base_hi, sizeof(double) * n);
    int ns = 0;
    for (int k = 0; k < n; k++) {
        for (int side = -1; side <= 1; side += 2) {
            int cut = side < 0 ? (r->cur_lo[k] < r->lo[k]) : (r->hi[k] < r->cur_hi[k]);
            if (!cut) continue;
            double *s_lo = r->slab_store + 4 * (int64_t)n * ns, *s_hi = s_lo + n, *c_lo = s_lo + 2 * n, *c_hi = s_lo + 3 * n;
            memcpy(s_lo, r->cur_lo, sizeof(double) * n); memcpy(s_hi, r->cur_hi, sizeof(double) * n);
            memcpy(c_lo, r->cur_lo, sizeof(double) * n); memcpy(c_hi, r->cur_hi, sizeof(double) * n);
            kc_slab *sl = &r->slabs[ns];
            sl->lo = s_lo; sl->hi = s_hi; sl->cut_lo = c_lo; sl->cut_hi = c_hi; sl->k = k; sl->side = side;
            if (side < 0) { s_hi[k] = r->lo[k]; sl->hint = r->hint_lower[k]; r->cur_lo[k] = r->lo[k]; }
            else          { s_lo[k] = r->hi[k]; sl->hint = r->hint_upper[k]; r->cur_hi[k] = r->hi[k]; }
            ns++;
        }
    }
    /* yield 'contract', base_lo, base_hi, (removed, (lo.copy(), hi.copy())) */
    out->kind = KC_NODE_CONTRACT; out->k = -1;
    out->lo = r->base_lo; out->hi = r->base_hi; out->new_lo = r->lo; out->new_hi = r->hi;
    out->nslabs = ns; out->slabs = r->slabs; out->hint = no_hint();
    return KC_REPLAY_NODE;
}

int kc_replay_next(kc_replay *r, kc_node *out)
{   /* One step of the reference's generator loop:
     *     while stack:
     *         lo, hi = stack.pop()
     *         while i < len(ev) and ev[i] == EV_ROUND: (round); yield 'contract', ...
     *         if i >= len(ev): raise ValueError('event stream ended inside the tree')
     *         e = int(ev[i]); i += 1
     *         leaf: yield; split 16 + k: yield 'split', push lower, push upper; else raise ValueError('unknown event')
     *     if i != len(ev) or pos['q'] != len(q) or pos['m'] != len(mus): raise ValueError('trailing events')
     * Claim: successive calls return the reference's nodes in order, then DONE or its error. */
    const int n = r->n;
    if (r->phase == PHASE_END) return r->err;
    memset(out, 0, sizeof *out);
    if (r->phase == PHASE_POP) {
        if (r->stack_len == 0) {
            if (r->i != r->nev || r->qpos != r->nq || r->mpos != r->nmu) return finish(r, KC_RERR_TRAILING);
            return finish(r, KC_REPLAY_DONE);
        }
        r->stack_len--;
        const double *e = r->stack + 2 * (int64_t)n * r->stack_len;
        memcpy(r->lo, e, sizeof(double) * n);
        memcpy(r->hi, e + n, sizeof(double) * n);
        r->phase = PHASE_NODE;
    }
    /* PHASE_NODE: the next contraction round of this node, or its own event */
    if (r->i < r->nev && r->ev[r->i] == KC_EV_ROUND) {
        int rc = contraction_round(r, out);
        if (rc != KC_REPLAY_NODE) return finish(r, rc);
        return KC_REPLAY_NODE;
    }
    if (r->i >= r->nev) return finish(r, KC_RERR_END_INSIDE);
    int e = r->ev[r->i]; r->i += 1;
    if (e >= 1 && e <= 8) {
        /* yield LEAF[e], lo, hi, (read_mu() if e == 3 and 'mu' in chart else None) */
        kc_hint h = no_hint();
        if (e == KC_LEAF_LINEAR && r->has_mu) { int rc = read_hint(r, &h); if (rc != 0) return finish(r, rc); }
        out->kind = e; out->k = -1; out->lo = r->lo; out->hi = r->hi; out->hint = h;
        r->phase = PHASE_POP;
        return KC_REPLAY_NODE;
    }
    if (e >= KC_EV_SPLIT && e < KC_EV_CONTRACT) {
        /* k = e - EV_SPLIT; mid = 0.5 * (lo[k] + hi[k]) (IndexError for k >= n, before the yield);
         * lower child [lo_k, mid], upper child [mid, hi_k]; stack.append(lower); stack.append(upper), so the
         * upper child is popped first: its subtree is next in the stream */
        int k = e - KC_EV_SPLIT;
        if (k >= n) return finish(r, KC_RERR_INDEX);
        double mid = midpoint(r->lo[k], r->hi[k]);
        double save_lo = r->lo[k], save_hi = r->hi[k];
        int rc;
        r->hi[k] = mid; rc = push(r, r->lo, r->hi); r->hi[k] = save_hi;           /* lower child */
        if (rc != 0) return finish(r, rc);
        r->lo[k] = mid; rc = push(r, r->lo, r->hi); r->lo[k] = save_lo;           /* upper child */
        if (rc != 0) return finish(r, rc);
        out->kind = KC_NODE_SPLIT; out->k = k; out->lo = r->lo; out->hi = r->hi; out->hint = no_hint();
        r->phase = PHASE_POP;
        return KC_REPLAY_NODE;
    }
    return finish(r, KC_RERR_UNKNOWN);
}
