/* kc_leaf.c: the leaf tests and the per-chart driver of the C checker kernel (contract in kc_leaf.h).
 *
 * Status: trusted.
 * Implements: checklist F4; docs/lemmas.md L15 (at a realization every vertex angle sum is exactly 2 pi), L16
 * (the mean-value inclusion with a given multiplier vector mu), the model's validity flag (L3, L6-L8, L10'),
 * the monotonicity tightening of L16 (through kc_corners), and bisection; paper Appendix D (certification); the
 * order of the tests in docs/code-tour.md section 0.1.  The non-LP tests of kisscrit/core/checker.py, written
 * with the same calls of the model (kc_model, checklist F3a), the same interval operations in the same order and
 * the same endpoint comparisons; that the results are equal bit for bit is computed (tests/test_ckernel.py,
 * scripts/completeness/diff_kernel.py), not proved.  Comments of the form "checker.py: ..." name the reference
 * line that the code below them mirrors.
 *
 * Arithmetic: only kc_arith (IEEE operations stepped one ulp outward; Arb for elementary functions, rounded
 * outward) and kc_model.  Outside those, plain double operations appear only where the reference has them, as
 * separate statements in the reference's order: the centre m = 0.5 (a + b) (one rounding, of the sum; the
 * halving is exact), the widths b - a and the bisection midpoint 0.5 (lo + hi) (which only choose boxes, never
 * decide a verdict), and the terminal width test.  A verdict is decided only by comparing interval endpoints
 * with 0, or by the model's validity flag.
 * Build with -O2 -fno-fast-math -ffp-contract=off (csrc/Makefile); kc_init asserts round-to-nearest and the
 * absence of FMA contraction.
 */
#define _POSIX_C_SOURCE 199309L        /* clock_gettime, for the optional per-unit timing only */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "kc_leaf.h"

struct kc_leaf_ws {
    const kc_graph *g;
    int n, npar, nfaces, has_pentagon;
    kc_ad *C;        /* corners of the box under test (5 per face) */
    kc_ad *Cm;       /* corners at the centre m of the context box B */
    kc_iv *s;        /* vertex sums (n) */
    kc_iv *grad;     /* vertex-sum gradients (n * npar), scratch */
    double *m;       /* linear context: the centre of B (npar) */
    kc_iv *Fm;       /* linear context: F(m) = vertex sums at m (n) */
    kc_iv *JB;       /* linear context: J(B) = vertex-sum gradients on B (n * npar) */
    int64_t *last;   /* hint: index of the last occurrence of each row (n), -1 if none */
    kc_iv *r;        /* hint: r_j = sum_v mu_v J_vj (npar) */
    kc_iv *Y;        /* hint: Y_j encloses p_j - m_j over the tested box (npar) */
    double *bis;     /* bisection: per depth level the two children, 4 * npar doubles */
    double *width;   /* bisection: the widths of the box (npar) */
    int pieces;      /* bisection: pieces visited below the current top-level box (checker.py self.pieces) */
};

void kc_leaf_ws_free(kc_leaf_ws *w)
{
    if (w == NULL) return;
    free(w->C); free(w->Cm); free(w->s); free(w->grad); free(w->m); free(w->Fm); free(w->JB);
    free(w->last); free(w->r); free(w->Y); free(w->bis); free(w->width);
    free(w);
}

kc_leaf_ws *kc_leaf_ws_new(const kc_graph *g)
{   /* allocation only; no mathematical content */
    kc_leaf_ws *w = calloc(1, sizeof *w);
    if (w == NULL) return NULL;
    w->g = g; w->n = g->n; w->npar = g->npar; w->nfaces = g->nfaces;
    for (int f = 0; f < g->nfaces; f++) if (g->size[f] == 5) w->has_pentagon = 1;
    size_t n = (size_t)g->n, np = (size_t)g->npar, nc = 5 * (size_t)g->nfaces;
    w->C = calloc(nc, sizeof(kc_ad)); w->Cm = calloc(nc, sizeof(kc_ad));
    w->s = calloc(n, sizeof(kc_iv)); w->grad = calloc(n * np, sizeof(kc_iv));
    w->m = calloc(np, sizeof(double)); w->Fm = calloc(n, sizeof(kc_iv)); w->JB = calloc(n * np, sizeof(kc_iv));
    w->last = calloc(n, sizeof(int64_t)); w->r = calloc(np, sizeof(kc_iv)); w->Y = calloc(np, sizeof(kc_iv));
    w->bis = calloc((size_t)(KC_MAX_DEPTH + 1) * 4 * np, sizeof(double));
    w->width = calloc(np, sizeof(double));
    if (!w->C || !w->Cm || !w->s || !w->grad || !w->m || !w->Fm || !w->JB || !w->last || !w->r || !w->Y ||
        !w->bis || !w->width) {
        kc_leaf_ws_free(w); return NULL;
    }
    return w;
}

/* ===================================================================== invalid and vertex sums */

static int excludes_zero(kc_iv x)
{   /* Claim: the interval [x.lo, x.hi] does not contain 0 when x.lo > 0 or x.hi < 0 (strict comparisons of
     * endpoints, as checker.py's `s.lo > 0 or s.hi < 0`).  A NaN endpoint makes both comparisons false, so a
     * NaN never excludes 0. */
    return x.lo > 0.0 || x.hi < 0.0;
}

int kc_leaf_invalid_vsum(kc_leaf_ws *w, const double *lo, const double *hi, int tighten, int *rule)
{   /* checker.py test_invalid_vsum(C, ok) with C, ok = G.corners(lo, hi, tighten=...).
     * Claim: *rule = 1 only if no point of the box is a valid configuration of the chart (the model's flag,
     * kc_model.h contract; L3, L6-L8, L10'); *rule = 2 only if for some vertex v the enclosure of
     * (sum of the corners at v) - 2 pi excludes 0 on the box, so no point of the box has all vertex sums equal
     * to 2 pi, which every realization has (L15).  Either way the box contains no critical configuration. */
    int ok = 0;
    int rc = kc_corners(w->g, lo, hi, NULL, tighten, w->C, &ok);
    if (rc != KC_OK) return rc;
    *rule = 0;
    if (!ok) { *rule = 1; return KC_OK; }                       /* checker.py: if not ok: return 'invalid' */
    kc_vertex_sums(w->g, w->C, w->s, w->grad);
    for (int v = 0; v < w->n; v++)                               /* checker.py: for s, g in G.vertex_sums(C) */
        if (excludes_zero(w->s[v])) { *rule = 2; return KC_OK; } /*   if s.lo > 0 or s.hi < 0: return 'vsum' */
    return KC_OK;
}

/* ===================================================================== the linear test with a hint */

int kc_leaf_context(kc_leaf_ws *w, const double *Blo, const double *Bhi)
{   /* checker.py linear_context(Blo, Bhi, CB) with CB = G.corners(Blo, Bhi, tighten=False), the caller's line
     * in check_certificate and in the linear-leaf branch.  The evaluation order is the reference's: CB first
     * (by the caller there), then the centre, then the corners at the centre.
     * Claim: Fm[v] encloses F_v(m) = (vertex sum at v) - 2 pi at the point m (corners at the point box [m, m],
     * clamp bound from Blo), and JB[v * npar + j] encloses dF_v/dp_j at every point of B where the clamped
     * corner functions are Lipschitz (kc_model.h contract; model.py docstring "Clamped functions"). */
    int ok = 0;
    int rc = kc_corners(w->g, Blo, Bhi, NULL, 0, w->C, &ok);     /* CB = G.corners(lo, hi, tighten=False) */
    if (rc != KC_OK) return rc;
    /* m lies in B, as the mean-value inclusion needs: for finite a <= b without overflow of a + b (the
     * coordinates are angles, |a|, |b| < 4), 2a <= a + b <= 2b exactly, rounding to nearest is monotone and
     * maps 2a and 2b to themselves, so 2a <= fl(a + b) <= 2b, and halving a double is exact. */
    for (int j = 0; j < w->npar; j++) {                          /* m = [0.5 * (a + b) for a, b in zip(Blo, Bhi)] */
        double t = Blo[j] + Bhi[j];
        w->m[j] = 0.5 * t;
    }
    rc = kc_corners(w->g, w->m, w->m, Blo, 1, w->Cm, &ok);       /* Cm, _ = G.corners(m, m, Blo, Bhi) */
    if (rc != KC_OK) return rc;
    kc_vertex_sums(w->g, w->Cm, w->Fm, w->grad);                 /* Fm = [s for s, g in G.vertex_sums(Cm)] */
    kc_vertex_sums(w->g, w->C, w->s, w->JB);                     /* JB = [g for s, g in G.vertex_sums(CB)] */
    return KC_OK;
}

int kc_leaf_hint(kc_leaf_ws *w, const double *lo, const double *hi, const uint16_t *rows, const float *vals,
                 int64_t count, int *closed, kc_iv *s_out)
{   /* checker.py test_linear_mu(lo, hi, Blo, Bhi, CB, mu, lctx) and combination(sorted(mu.items()), Fm, JB, Y).
     * Here F is the clamped vertex-sum map of the model (vertex sums minus 2 pi), which agrees with the true one
     * at every valid point of the chart's region (model.py docstring, "Clamped functions"; L10').
     * Claim: *closed = 1 only if no point p of [lo, hi] (a sub-box of the context's B) has F(p) = 0.
     * Argument (L16, the mean-value inclusion, in the "sum first" form that checker.combination evaluates):
     * fix real multipliers mu and let phi = mu.F, a scalar Lipschitz function on B (where the clamped model is
     * Lipschitz).  By Lebourg's mean value theorem, phi(p) - phi(m) = g.(p - m) for some g in the Clarke
     * generalized gradient of phi at a point of the segment [m, p], which lies in B (m is in B, see
     * kc_leaf_context, and B is convex).  By the chain rule for a linear functional of a Lipschitz map, that
     * gradient is contained in { mu^T A : A in the Clarke generalized Jacobian of F there }, and every such A
     * lies entrywise in J(B) (kc_model.h contract).  So g_j = sum_v mu_v A_vj lies in the interval
     * r_j = sum_v mu_v J_vj(B) computed below, phi(m) lies in sum_v mu_v Fm_v, p_j - m_j lies in Y_j, and
     * phi(p) lies in s.  If s excludes 0, phi(p) != 0, hence F(p) != 0.  (A row-by-row inclusion would give
     * the wider interval sum_v mu_v (Fm_v + sum_j J_vj Y_j); interval arithmetic is only subdistributive, so
     * the sum-first form needs the scalar argument above.  Precise references for the two theorems belong in
     * docs/lemmas.md L16; they are not restated here.)  Any mu is admissible, so the hint needs no trust.
     * Soundness does not need more; for information: where J has an unbounded entry (the collapsing diagonal,
     * where kc_model returns (-inf, inf)), the product with Y_j (never exactly [0, 0], since kc_sub steps
     * outward) is unbounded, so s cannot exclude 0. */
    *closed = 0;
    if (count <= 0) return KC_OK;                                /* checker.py: if not mu: return False */
    for (int64_t i = 0; i < count; i++)
        if ((int)rows[i] >= w->n) return KC_LERR_HINTROW;        /* the reference's Fm[v] raises IndexError */
    /* mu as a Python dict built from zip(rows, vals): a repeated row keeps its last value */
    for (int v = 0; v < w->n; v++) w->last[v] = -1;
    for (int64_t i = 0; i < count; i++) w->last[rows[i]] = i;
    int np = w->npar;
    /* Y = [I.hull(Iv(lo[j]) - Iv(m[j]), Iv(hi[j]) - Iv(m[j])) for j in range(n)]: each difference is an
     * outward-rounded enclosure of an endpoint minus m_j, and their hull encloses p_j - m_j for p_j in
     * [lo_j, hi_j] (p_j - m_j is increasing in p_j) */
    for (int j = 0; j < np; j++)
        w->Y[j] = kc_hull(kc_sub(kc_point(lo[j]), kc_point(w->m[j])), kc_sub(kc_point(hi[j]), kc_point(w->m[j])));
    /* combination: s = Iv(0.0); r = [Iv(0.0)] * n; then rows in increasing order (sorted(mu.items())) */
    kc_iv s = kc_point(0.0);
    for (int j = 0; j < np; j++) w->r[j] = kc_point(0.0);
    for (int v = 0; v < w->n; v++) {
        if (w->last[v] < 0) continue;
        kc_iv mv = kc_point((double)vals[w->last[v]]);           /* mv = Iv(float(val)): float32 -> double is exact */
        s = kc_add(s, kc_mul(mv, w->Fm[v]));                     /* s = s + mv * Fm[v] */
        const kc_iv *Jv = w->JB + (size_t)v * np;
        for (int j = 0; j < np; j++) w->r[j] = kc_add(w->r[j], kc_mul(mv, Jv[j]));   /* r[j] + mv * JB[v][j] */
    }
    for (int j = 0; j < np; j++) s = kc_add(s, kc_mul(w->r[j], w->Y[j]));            /* s = s + r[j] * Y[j] */
    if (s_out != NULL) *s_out = s;
    *closed = excludes_zero(s);                                  /* return s.lo > 0 or s.hi < 0 */
    return KC_OK;
}

/* ===================================================================== exclude without linear programs */

static int argmax_first(const double *x, int len)
{   /* int(np.argmax(w)) for a list of doubles: the first index of the maximum; a NaN counts as the maximum,
     * so the first NaN wins (computed with numpy 2.5.3: argmax([1, nan, 2, nan]) = 1, argmax([3, 1, 3]) = 0,
     * argmax([nan, 5]) = 0).  Only chooses the split coordinate. */
    int k = 0;
    for (int j = 0; j < len; j++) {
        if (isnan(x[j])) return j;
        if (x[j] > x[k]) k = j;
    }
    return k;
}

static int exclude_rec(kc_leaf_ws *w, const double *lo, const double *hi, int depth, int max_depth, int max_pieces,
                       int *verdict, int64_t *cnt)
{   /* checker.py exclude(lo, hi, parent, depth) without test_linear (both calls), test_gordan, test_motzkin.
     * Claim: *verdict != KC_V_UNCLOSED only if [lo, hi] contains no critical configuration: either one of the
     * tests of kc_leaf_invalid_vsum proves it on the whole box, or the box is covered by the two halves of
     * its widest coordinate (lower [lo_k, mid], upper [mid, hi_k]) and each half is proved the same way,
     * recursively.  Cover: for finite lo_k <= hi_k without overflow of lo_k + hi_k, lo_k <= mid <= hi_k (as for
     * m in kc_leaf_context), so the halves cover the box; if the sum overflowed to +inf, mid = +inf >= hi_k and
     * the lower half alone contains the box.  The argument assumes finite endpoints; with an infinite endpoint
     * the widths and midpoint can be inf or NaN and the cover is not claimed here (boxes of the format come
     * from a finite starting box by grid contraction and midpoint splits; a finiteness check on lo0, hi0 is
     * an open proposal for the front end and the reference together). */
    *verdict = KC_V_UNCLOSED;
    w->pieces += 1;                                              /* self.pieces += 1 */
    if (w->pieces > max_pieces) return KC_OK;                    /* if self.pieces > MAX_PIECES: return False */
    int rule = 0;
    int rc = kc_leaf_invalid_vsum(w, lo, hi, 0, &rule);          /* cheap pass: corners(lo, hi, tighten=False) */
    if (rc != KC_OK) return rc;
    if (rule) {
        *verdict = (rule == 1) ? KC_V_INVALID : KC_V_VSUM;
        cnt[(rule == 1) ? KC_C_INVALID : KC_C_VSUM] += 1;
        return KC_OK;
    }
    /* (the reference tries test_linear with the parent's context and with its own box here: an LP, left out) */
    if (w->has_pentagon) {                                       /* if any(k == 5 for k in self.G.sizes) */
        rc = kc_leaf_invalid_vsum(w, lo, hi, 1, &rule);          /* C, ok = self.evaluate(lo, hi) */
        if (rc != KC_OK) return rc;
        if (rule) {
            *verdict = (rule == 1) ? KC_V_INVALID_T : KC_V_VSUM_T;
            cnt[(rule == 1) ? KC_C_INVALID_T : KC_C_VSUM_T] += 1;
            return KC_OK;
        }
    }
    /* (the reference tries test_gordan and test_motzkin here: LPs, left out) */
    if (depth >= max_depth) return KC_OK;                        /* if depth >= MAX_DEPTH: return False */
    int np = w->npar;
    for (int j = 0; j < np; j++) w->width[j] = hi[j] - lo[j];    /* w = [b - a for a, b in zip(lo, hi)] */
    int k = argmax_first(w->width, np);                          /* k = int(np.argmax(w)) */
    if (w->width[k] <= 0.0) return KC_OK;                        /* if w[k] <= 0: return False */
    double t = lo[k] + hi[k];
    double mid = 0.5 * t;                                        /* mid = 0.5 * (lo[k] + hi[k]) */
    double *lo1 = w->bis + (size_t)depth * 4 * np, *hi1 = lo1 + np, *lo2 = hi1 + np, *hi2 = lo2 + np;
    memcpy(lo1, lo, sizeof(double) * np); memcpy(hi1, hi, sizeof(double) * np); hi1[k] = mid;
    memcpy(lo2, lo, sizeof(double) * np); memcpy(hi2, hi, sizeof(double) * np); lo2[k] = mid;
    cnt[KC_C_BISECT] += 1;                                       /* self.bump('bisect') */
    int v1 = KC_V_UNCLOSED, v2 = KC_V_UNCLOSED;
    rc = exclude_rec(w, lo1, hi1, depth + 1, max_depth, max_pieces, &v1, cnt);   /* lower half first */
    if (rc != KC_OK) return rc;
    if (v1 == KC_V_UNCLOSED) return KC_OK;                       /* `and` short-circuits */
    rc = exclude_rec(w, lo2, hi2, depth + 1, max_depth, max_pieces, &v2, cnt);
    if (rc != KC_OK) return rc;
    if (v2 == KC_V_UNCLOSED) return KC_OK;
    *verdict = KC_V_BISECT;
    return KC_OK;
}

int kc_leaf_exclude(kc_leaf_ws *w, const double *lo, const double *hi, int max_depth, int max_pieces,
                    int *verdict, int64_t *counts)
{   /* checker.py exclude(lo, hi) at depth 0 (pieces reset), without its LP tests; claim as exclude_rec.
     * The counters of a box that is not closed are discarded (kc_leaf.h). */
    *verdict = KC_V_UNCLOSED;
    if (max_depth < 0 || max_depth > KC_MAX_DEPTH || max_pieces < 0 || max_pieces > KC_MAX_PIECES)
        return KC_LERR_ARG;
    int64_t cnt[KC_NCOUNT] = { 0 };
    w->pieces = 0;                                               /* if depth == 0: self.pieces = 0 */
    int rc = exclude_rec(w, lo, hi, 0, max_depth, max_pieces, verdict, cnt);
    if (rc != KC_OK) { *verdict = KC_V_UNCLOSED; return rc; }
    if (*verdict != KC_V_UNCLOSED && counts != NULL)
        for (int i = 0; i < KC_NCOUNT; i++) counts[i] += cnt[i];
    return KC_OK;
}

/* ===================================================================== the chart */

static int64_t now_ns(void)
{   /* wall clock for the optional timing; never part of a verdict */
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (int64_t)ts.tv_sec * 1000000000 + ts.tv_nsec;
}

static int terminal_too_wide(const double *lo, const double *hi, int n, double tol)
{   /* checker.py: max(b - a for a, b in zip(lo, hi)) > float(cert['tol']).  Python's max keeps the first
     * item unless a later one compares strictly greater; then one strict comparison with tol.  This is a
     * bound on the search's terminal boxes (they are recorded for identification, never excluded). */
    if (n <= 0) return 0;                       /* the reference's max() of nothing raises; n >= 1 always (theta) */
    double mx = hi[0] - lo[0];
    for (int j = 1; j < n; j++) {
        double d = hi[j] - lo[j];
        if (d > mx) mx = d;
    }
    return mx > tol;
}

static int model_error(int rc, int *model_err)
{   /* map a kc_corners code to the driver's error, keeping the model's code */
    if (rc > 0) { if (model_err) *model_err = rc; return KC_LERR_MODEL; }
    return rc;
}

static int check_unit(kc_leaf_ws *w, const double *lo, const double *hi, const kc_hint *hint, int own_context,
                      int max_depth, int max_pieces, int *verdict, int64_t *counts)
{   /* One leaf (events 1-6) or one slab.  checker.py check_certificate (linear leaf): if kind == 'linear' and
     * extra: C_ = corners(lo, hi, tighten=False); if test_linear_mu(lo, hi, lo, hi, C_, extra): continue; then
     * exclude(lo, hi).  exclude_slab: if test_linear_mu(lo, hi, Blo, Bhi, CB, mu, lctx): return True; (the LP
     * test_linear_bound, left out); return exclude(lo, hi, parent).  A leaf's hint is present only for a linear
     * leaf (event 3) of a chart with a mu stream (kc_replay.h, kc_node), as the reference's
     * `kind == 'linear' and extra` requires.  own_context = 1 for the leaf (B is the
     * leaf's own box), 0 for a slab (the context of the contraction's box is already in the workspace).
     * Claim: *verdict != KC_V_UNCLOSED only if the box contains no critical configuration (the claims of
     * kc_leaf_hint and kc_leaf_exclude). */
    *verdict = KC_V_UNCLOSED;
    if (hint->present && hint->count > 0) {
        if (own_context) {
            int rc = kc_leaf_context(w, lo, hi);
            if (rc != KC_OK) return rc;
        }
        int closed = 0;
        int rc = kc_leaf_hint(w, lo, hi, hint->rows, hint->vals, hint->count, &closed, NULL);
        if (rc != KC_OK) return rc;
        if (closed) { *verdict = KC_V_HINT; counts[KC_C_HINT] += 1; return KC_OK; }
    }
    return kc_leaf_exclude(w, lo, hi, max_depth, max_pieces, verdict, counts);
}

int kc_check_chart(const kc_graph *g, int n, const double *lo0, const double *hi0,
                   const uint8_t *ev, int64_t nev, const uint16_t *q, int64_t nq,
                   const float *mu, int64_t nmu, int has_mu, double tol, int max_depth, int max_pieces,
                   int8_t *verdict, int8_t *unit_kind, int64_t cap, int64_t *nunits, int64_t *counts, int64_t *unit_ns,
                   int *model_err)
{   /* checker.py check_certificate, the body of `for ch in cert['charts']` after its two checks of the
     * starting box (done by the front end), without the LP tests.  Claim: every unit with a closed verdict
     * contains no critical configuration; terminal and open units are reported, never closed.  Together with
     * the replay's coverage (docs/certificate-format.md section 3: the leaves and slabs of a complete tree
     * cover the starting box), a chart all of whose units are closed or terminal within tol has no critical
     * configuration outside its terminal boxes. */
    *nunits = 0;
    if (model_err) *model_err = 0;
    for (int i = 0; i < KC_NCOUNT; i++) counts[i] = 0;
    if (n != g->npar) return KC_LERR_SIZE;                       /* 'parameter count mismatch' */
    if (max_depth < 0 || max_depth > KC_MAX_DEPTH || max_pieces < 0 || max_pieces > KC_MAX_PIECES)
        return KC_LERR_ARG;
    kc_leaf_ws *w = kc_leaf_ws_new(g);
    int64_t stack_cap = kc_replay_stack_need(ev, nev);
    int64_t nwork = kc_replay_work_doubles(n, stack_cap);
    double *work = malloc(sizeof(double) * (size_t)(nwork > 0 ? nwork : 1));
    kc_replay *r = malloc(sizeof *r);
    int result = 0;
    if (w == NULL || work == NULL || r == NULL) { result = KC_LERR_NOMEM; goto done; }
    result = kc_replay_init(r, n, lo0, hi0, ev, nev, q, nq, mu, nmu, has_mu, work, stack_cap);
    if (result != 0) goto done;
    kc_node node;
    for (;;) {
        int rc = kc_replay_next(r, &node);
        if (rc == KC_REPLAY_DONE) { result = 0; break; }
        if (rc < 0) { result = rc; break; }
        if (node.kind == KC_NODE_SPLIT) continue;
        if (node.kind == KC_NODE_CONTRACT) {
            /* CB, _ = G.corners(lo, hi, tighten=False); lctx = K.linear_context(lo, hi, CB) */
            rc = kc_leaf_context(w, node.lo, node.hi);
            if (rc != KC_OK) { result = model_error(rc, model_err); break; }
            for (int i = 0; i < node.nslabs; i++) {
                if (*nunits >= cap) { result = KC_LERR_CAP; goto done; }
                const kc_slab *sl = &node.slabs[i];
                int64_t t0 = unit_ns ? now_ns() : 0;
                int v = KC_V_UNCLOSED;
                rc = check_unit(w, sl->lo, sl->hi, &sl->hint, 0, max_depth, max_pieces, &v, counts);
                if (rc != KC_OK) { result = model_error(rc, model_err); goto done; }
                if (unit_ns) unit_ns[*nunits] = now_ns() - t0;
                if (unit_kind) unit_kind[*nunits] = 0;
                verdict[(*nunits)++] = (int8_t)v;
            }
            continue;
        }
        /* a leaf, kind 1..8 */
        if (*nunits >= cap) { result = KC_LERR_CAP; break; }
        int64_t t0 = unit_ns ? now_ns() : 0;
        int v;
        if (node.kind == 7)                                      /* terminal: recorded, width bounded by tol */
            v = terminal_too_wide(node.lo, node.hi, n, tol) ? KC_V_TERM_WIDE : KC_V_TERMINAL;
        else if (node.kind == 8)                                 /* open: never closed */
            v = KC_V_OPEN;
        else {
            rc = check_unit(w, node.lo, node.hi, &node.hint, 1, max_depth, max_pieces, &v, counts);
            if (rc != KC_OK) { result = model_error(rc, model_err); break; }
        }
        if (unit_ns) unit_ns[*nunits] = now_ns() - t0;
        if (unit_kind) unit_kind[*nunits] = (int8_t)node.kind;
        verdict[(*nunits)++] = (int8_t)v;
    }
done:
    kc_leaf_ws_free(w); free(work); free(r);
    return result;
}
