/* kc_model.c: the face-angle model for the C checker kernel (contract in kc_model.h).
 *
 * Status: trusted.
 * Implements: checklist F3a; lemmas L3, L6, L7, L8, L10'/L11', L15 and the monotonicity tightening of L16
 * (docs/lemmas.md); paper Appendix A, Lemmas A.2 and A.6-A.8.  A bit-for-bit reimplementation of
 * kisscrit/core/model.py: every function below has the reference's name (with a trailing underscore where C
 * needs one), performs the same interval operations in the same order, and so returns the same doubles.
 * tests/test_cmodel.py compares the two on candidates and charts from the `cond` certificates.
 *
 * Arithmetic: only kc_arith (IEEE operations stepped one ulp outward; Arb at KC_PREC = 64 bits for the
 * elementary functions, rounded outward).  Comparisons of endpoints decide flags and branches exactly as in
 * the reference; no floating-point value decides a verdict except through an enclosure.
 * Build with -O2 -fno-fast-math -ffp-contract=off (csrc/Makefile).
 */
#include <stdlib.h>
#include <string.h>
#include <flint/arb.h>
#include "kc_model.h"

/* LSTRAIGHT = pi - 25 degrees, the same double as model.py's math.pi - math.radians(25.0) (CPython computes
 * radians(x) as x * (pi / 180.0)).  It is a region threshold, not an enclosure: the pentagon charts are cut at
 * this double, and the reference uses the same one.  tests/test_cmodel.py checks equality. */
static const double PI_DOUBLE = 3.141592653589793;           /* math.pi */
const double KC_LSTRAIGHT = PI_DOUBLE - 25.0 * (PI_DOUBLE / 180.0);
#define LSTRAIGHT KC_LSTRAIGHT

static const int SHIFT_A = 0, SHIFT_B = 2, SHIFT_C = 3;      /* model.SHIFT: chart rotations of the pentagon */

/* ===================================================================== interval helpers (model.py top) */

static kc_iv Z(void) { return kc_point(0.0); }                /* the exact zero [0, 0] */
static kc_iv ZERO_ONE(void) { kc_iv r = { 0.0, 1.0 }; return r; }

static int intersect(kc_iv x, kc_iv y, kc_iv *out)
{   /* Intersection of two enclosures of the same quantity.  Two valid enclosures always meet, so disjoint
     * ones mean a bug in one of them: report it (the reference raises ArithmeticError). */
    double lo = (y.lo > x.lo) ? y.lo : x.lo;        /* Python max(x.lo, y.lo) */
    double hi = (y.hi < x.hi) ? y.hi : x.hi;        /* Python min(x.hi, y.hi) */
    if (lo > hi) return KC_ERR_DISJOINT;
    out->lo = lo; out->hi = hi;
    return KC_OK;
}

/* ===================================================================== forward AD, 3 local components */

static kc_ad var(double lo, double hi, int k)
{   /* the k-th local variable on [lo, hi]: value [lo, hi], gradient the k-th unit vector (exact) */
    kc_ad r; r.v.lo = lo; r.v.hi = hi;
    for (int i = 0; i < 3; i++) r.g[i] = Z();
    r.g[k] = kc_point(1.0);
    return r;
}

static kc_ad const_(kc_iv c)
{   /* a constant: the value c and a zero gradient */
    kc_ad r; r.v = c; r.g[0] = r.g[1] = r.g[2] = Z();
    return r;
}

static kc_ad add(kc_ad a, kc_ad b)
{   /* sum: value and gradient intervals add (each an enclosure of the sum of the enclosed quantities) */
    kc_ad r; r.v = kc_add(a.v, b.v);
    for (int i = 0; i < 3; i++) r.g[i] = kc_add(a.g[i], b.g[i]);
    return r;
}

static kc_ad sub(kc_ad a, kc_ad b)
{   /* difference */
    kc_ad r; r.v = kc_sub(a.v, b.v);
    for (int i = 0; i < 3; i++) r.g[i] = kc_sub(a.g[i], b.g[i]);
    return r;
}

static kc_ad scale(kc_ad a, double c)
{   /* a times the real constant c (c is a double, so [c, c] is exact) */
    kc_ad r; r.v = kc_mul(a.v, kc_point(c));
    for (int i = 0; i < 3; i++) r.g[i] = kc_mul(a.g[i], kc_point(c));
    return r;
}

static kc_ad shift(kc_ad a, double c)
{   /* a plus the real constant c; the gradient is unchanged */
    kc_ad r = a; r.v = kc_add(a.v, kc_point(c));
    return r;
}

static kc_iv sum0(kc_iv x, kc_iv y)
{   /* x + y keeping an exact zero exact, so structural zeros of gradients survive (model._sum0) */
    if (kc_is_zero(x)) return y;
    if (kc_is_zero(y)) return x;
    return kc_add(x, y);
}

static kc_ad mul(kc_ad a, kc_ad b)
{   /* product rule: value a b, gradient a' b + a b' (each term an interval enclosure) */
    kc_ad r; r.v = kc_mul(a.v, b.v);
    for (int i = 0; i < 3; i++) r.g[i] = sum0(kc_mul(a.g[i], b.v), kc_mul(b.g[i], a.v));
    return r;
}

static kc_ad div_(kc_ad a, kc_ad b)
{   /* quotient rule: value q = a / b, gradient (a' - b' q) / b */
    kc_ad r; kc_iv q = kc_div(a.v, b.v);
    r.v = q;
    for (int i = 0; i < 3; i++) r.g[i] = kc_div(sum0(a.g[i], kc_neg(kc_mul(b.g[i], q))), b.v);
    return r;
}

static kc_ad chain(kc_ad a, kc_iv v, kc_iv d)
{   /* chain rule: a function with value enclosure v and derivative enclosure d at a; gradient a' d */
    kc_ad r; r.v = v;
    for (int i = 0; i < 3; i++) r.g[i] = kc_mul(a.g[i], d);
    return r;
}

static kc_ad cos_(kc_ad a)
{   /* cos, derivative -sin */
    return chain(a, kc_cos(a.v), kc_neg(kc_sin(a.v)));
}

static kc_ad sin_(kc_ad a)
{   /* sin, derivative cos */
    return chain(a, kc_sin(a.v), kc_cos(a.v));
}

static kc_ad atan_(kc_ad a)
{   /* arctan, derivative 1 / (1 + a^2).  On a non-finite argument the derivative is enclosed by [0, 1], which
     * contains 1 / (1 + t^2) for every real t. */
    kc_iv x = a.v;
    if (!kc_finite(x)) return chain(a, kc_atan(x), ZERO_ONE());
    return chain(a, kc_atan(x), kc_div(kc_point(1.0), kc_add(kc_point(1.0), kc_square(x))));
}

static kc_ad clip(kc_ad a, double lo, double hi, int *ok)
{   /* clamp to [lo, hi]: the value is the image of the clamp; the gradient is multiplied by [0, 1] unless the
     * clamp certainly does not act on the box (then the function is the identity there).  *ok = 0 iff the
     * argument certainly lies outside [lo, hi]. */
    kc_iv x = a.v;
    *ok = !(x.hi < lo) && !(x.lo > hi);
    int inside = x.lo > lo && x.hi < hi;
    kc_ad r; r.v = kc_clamp(x, lo, hi);
    for (int i = 0; i < 3; i++) r.g[i] = inside ? a.g[i] : kc_mul(a.g[i], ZERO_ONE());
    return r;
}

static kc_ad asin_angle(kc_ad a, double lo, double hi, int *ok)
{   /* y = clamp(asin(clamp(x, -1, 1)), lo, hi) with [lo, hi] inside [-pi/2, pi/2); dy/dx = 1 / cos(y), times
     * [0, 1] where a clamp may act.  *ok = 0 iff x is certainly outside [-1, 1] or asin x certainly outside
     * [lo, hi]. */
    kc_iv x = a.v;
    int o = !(x.hi < -1.0) && !(x.lo > 1.0);
    kc_iv y0 = kc_asin(x);
    o = o && !(y0.hi < lo) && !(y0.lo > hi);
    int inside = x.lo > -1.0 && x.hi < 1.0 && y0.lo > lo && y0.hi < hi;
    kc_iv y = kc_clamp(y0, lo, hi);
    kc_iv d = kc_div(kc_point(1.0), kc_cos(y));
    if (!inside) d = kc_mul(d, ZERO_ONE());
    kc_ad r; r.v = y;
    for (int i = 0; i < 3; i++) r.g[i] = kc_mul(a.g[i], d);
    *ok = o;
    return r;
}

static kc_ad acos_(kc_ad a, int *ok)
{   /* arccos, derivative -1 / sin(arccos a); *ok = 0 iff the argument is certainly outside [-1, 1] */
    kc_iv x = a.v;
    *ok = !(x.hi < -1.0) && !(x.lo > 1.0);
    kc_iv y = kc_acos(x);
    return chain(a, y, kc_div(kc_point(-1.0), kc_sin(y)));
}

/* ===================================================================== faces */

static void alpha_arb(arb_t r, double t)
{   /* alpha(t) = arccos(cos t / (1 + cos t)) as an Arb ball at KC_PREC bits: the corner of the equilateral
     * spherical triangle with side t (docs/lemmas.md, setting; paper Appendix A).  Same operations as
     * model._alpha_arb: cos, add 1, divide, acos. */
    arb_t c, d; arb_init(c); arb_init(d);
    arb_set_d(c, t);
    arb_cos(c, c, KC_PREC);
    arb_add_si(d, c, 1, KC_PREC);
    arb_div(r, c, d, KC_PREC);
    arb_acos(r, r, KC_PREC);
    arb_clear(c); arb_clear(d);
}

static kc_iv alpha_range(double tlo, double thi)
{   /* alpha is increasing on (0, pi/2), so [alpha(tlo) rounded down, alpha(thi) rounded up] is its range */
    arb_t a; arb_init(a);
    kc_iv r;
    alpha_arb(a, tlo); r.lo = kc_down(a);
    alpha_arb(a, thi); r.hi = kc_up(a);
    arb_clear(a);
    return r;
}

static void beta_pt(arb_t r, double t, double A)
{   /* beta(t, A) = arctan(cos(A/2) / (cos t sin(A/2))) as an Arb ball at a point: the base angle of the
     * isosceles spherical triangle with legs t and apex A (L7; paper Appendix A, Lemma A.7).  A/2 is exact. */
    arb_t h, num, den, ct; arb_init(h); arb_init(num); arb_init(den); arb_init(ct);
    arb_set_d(h, A); arb_mul_2exp_si(h, h, -1);
    arb_cos(num, h, KC_PREC);
    arb_set_d(ct, t); arb_cos(ct, ct, KC_PREC);
    arb_sin(den, h, KC_PREC);
    arb_mul(den, ct, den, KC_PREC);
    arb_div(r, num, den, KC_PREC);
    arb_atan(r, r, KC_PREC);
    arb_clear(h); arb_clear(num); arb_clear(den); arb_clear(ct);
}

static kc_ad beta(kc_ad ct, kc_ad A)
{   /* base angle of the isosceles spherical triangle with legs theta and apex angle A:
     * beta = arctan(cos(A/2) / (cos theta sin(A/2))) (L7) */
    kc_ad h = scale(A, 0.5);
    return atan_(div_(cos_(h), mul(ct, sin_(h))));
}

static void pentagon(kc_ad ct, kc_ad st, double hmin, double hmax, kc_ad f1, kc_ad f2, kc_ad cs[5], int *ok)
{   /* Corners (formula positions 0..4) of the convex equilateral pentagon P0..P4 with corners f1, f2 at P1,
     * P2 (L8, paper Lemma A.8).  P0 P2 is the base of the isosceles triangle P0 P1 P2; P3 follows from the
     * angle g at P2 between P2 P0 and P2 P3; P4 is the apex of the isosceles triangle on the diagonal P0 P3,
     * on the convex side.  *ok = 0 only when no valid pentagon has parameters in the box. */
    int o, okk;
    /* d2/2 = asin(sin theta sin(f1/2)) (L8, step 1).  The clamp to [0, HALF_PI.lo] does not act at valid
     * points: d2 <= 2 theta < 127 deg, so d2/2 < pi/2 by a wide margin. */
    kc_ad hd2 = asin_angle(mul(st, sin_(scale(f1, 0.5))), 0.0, KC_HALF_PI.lo, &okk);
    kc_ad d2 = scale(hd2, 2.0);
    kc_ad b1 = beta(ct, f1);
    /* g = f2 - b1 lies in [0, pi] for a convex pentagon (L8, step 2).  The clamp to [0, PI.lo] does not act at
     * valid points: it would need two adjacent corners within an ulp of pi, excluded by L10. */
    kc_ad g = clip(sub(f2, b1), 0.0, KC_PI.lo, &o); okk = okk && o;
    kc_ad cg = cos_(g), sg = sin_(g);
    kc_ad cd2 = cos_(d2), sd2 = sin_(d2);
    /* spherical law of cosines in the triangle P0 P2 P3: cos d3 = cos d2 cos theta + sin d2 sin theta cos g */
    kc_ad cd3 = add(mul(cd2, ct), mul(mul(sd2, st), cg));
    okk = okk && !(cd3.v.lo > ct.v.hi);                 /* the diagonal P0 P3 is >= theta somewhere */
    cd3 = clip(cd3, -1.0, 1.0, &o);                      /* flag ignored, as in the reference */
    kc_ad d3 = acos_(cd3, &o); okk = okk && o;
    /* the angles at P0 and P3 of the triangle P0 P2 P3 (cotangent formulas, L8 step 3) */
    kc_ad a0 = sub(const_(KC_HALF_PI),
                   atan_(div_(sub(mul(sd2, ct), mul(mul(cd2, st), cg)), mul(sg, st))));
    kc_ad a3 = sub(const_(KC_HALF_PI),
                   atan_(div_(sub(mul(st, cd2), mul(mul(ct, sd2), cg)), mul(sg, sd2))));
    /* the apex: sin(f4/2) = sin(d3/2) / sin theta, with f4/2 clamped to the region [hmin, hmax] (L10') */
    kc_ad h = asin_angle(div_(sin_(scale(d3, 0.5)), st), hmin, hmax, &o); okk = okk && o;
    kc_ad f4 = scale(h, 2.0);
    kc_ad b3 = beta(ct, f4);
    cs[0] = add(add(b1, a0), b3);
    cs[1] = f1;
    cs[2] = f2;
    cs[3] = add(a3, b3);
    cs[4] = f4;
    *ok = okk;
}

static void pick(kc_iv g, double lo, double hi, int want_min, double out[2])
{   /* where a function with derivative enclosure g on [lo, hi] attains its minimum (want_min) or maximum: an
     * endpoint if g has a sign, else the whole interval */
    if (g.lo >= 0) { out[0] = out[1] = want_min ? lo : hi; return; }
    if (g.hi <= 0) { out[0] = out[1] = want_min ? hi : lo; return; }
    out[0] = lo; out[1] = hi;
}

static int tighten(kc_ad c, int j, double tlo, double thi, double alo, double ahi, double blo, double bhi,
                   double hmin, double hmax, kc_ad *out)
{   /* Monotonicity tightening of the computed pentagon corner j (L16): where a partial derivative enclosure
     * has a sign on the box, the clamped corner function (Lipschitz on the box) is monotone in that variable by
     * Lebourg's theorem, so its extreme values are attained on the corresponding face of the box; evaluate
     * there with the same clamp bounds.  The validity flags of the face evaluations are ignored on purpose:
     * the clamped formulas enclose the clamped function on every face.  Returns the corner with a narrower
     * value and the same gradient. */
    double lo = c.v.lo, hi = c.v.hi;
    for (int want_min = 1; want_min >= 0; want_min--) {
        double a[2], b[2], d[2];
        pick(c.g[0], tlo, thi, want_min, a);
        pick(c.g[1], alo, ahi, want_min, b);
        pick(c.g[2], blo, bhi, want_min, d);
        if (a[0] == tlo && a[1] == thi && b[0] == alo && b[1] == ahi && d[0] == blo && d[1] == bhi) continue;
        kc_ad th = var(a[0], a[1], 0);
        kc_ad ct = cos_(th), st = sin_(th);
        kc_ad cs[5]; int ignored;
        pentagon(ct, st, hmin, hmax, var(b[0], b[1], 1), var(d[0], d[1], 2), cs, &ignored);
        kc_iv v = cs[j].v;
        if (want_min) { if (v.lo > lo) lo = v.lo; }      /* Python max(lo, v.lo) */
        else          { if (v.hi < hi) hi = v.hi; }      /* Python min(hi, v.hi) */
    }
    if (lo > hi) return KC_ERR_TIGHTEN;
    out->v.lo = lo; out->v.hi = hi;
    out->g[0] = c.g[0]; out->g[1] = c.g[1]; out->g[2] = c.g[2];
    return KC_OK;
}

/* ===================================================================== the graph */

kc_graph *kc_graph_new(int nfaces, const int *sizes, const int *verts, const char *modes, int *err)
{   /* model.Graph.__init__: parameter indices; the rotation system at each vertex (successor map from the
     * oriented faces, walked from the smallest neighbour); the edge set; the Euler check n - E + F = 2. */
    *err = KC_ERR_LIMIT;
    if (nfaces <= 0 || nfaces > KC_MAXF) return NULL;
    kc_graph *g = calloc(1, sizeof(kc_graph));
    if (!g) return NULL;
    g->nfaces = nfaces;
    int off = 0, n = 1, nv = 0;
    for (int f = 0; f < nfaces; f++) {
        int k = sizes[f];
        if (k < 3 || k > 5) goto fail;
        g->size[f] = k;
        for (int q = 0; q < k; q++) {
            int v = verts[off + q];
            if (v < 0 || v >= KC_MAXV) goto fail;
            g->face[f][q] = v;
            if (v + 1 > nv) nv = v + 1;
        }
        off += k;
        g->mode[f] = modes ? modes[f] : 0;
        g->pidx[f][0] = g->pidx[f][1] = -1;
        if (k == 4) { g->pidx[f][0] = n; n += 1; }
        else if (k == 5) { g->pidx[f][0] = n; g->pidx[f][1] = n + 1; n += 2; }
    }
    g->npar = n;
    g->n = nv;
    /* successor maps: succ[v][b] = (a, fi, pos) for each face (..., a, v, b, ...) */
    static const int NONE = -1;
    int (*sa)[KC_MAXV] = malloc(sizeof(int[KC_MAXV][KC_MAXV]));
    int (*sf)[KC_MAXV] = malloc(sizeof(int[KC_MAXV][KC_MAXV]));
    int (*sp)[KC_MAXV] = malloc(sizeof(int[KC_MAXV][KC_MAXV]));
    if (!sa || !sf || !sp) { free(sa); free(sf); free(sp); goto fail; }
    for (int v = 0; v < KC_MAXV; v++) for (int w = 0; w < KC_MAXV; w++) sa[v][w] = NONE;
    for (int f = 0; f < nfaces; f++) {
        int k = g->size[f];
        for (int pos = 0; pos < k; pos++) {
            int v = g->face[f][pos], a = g->face[f][(pos + k - 1) % k], b = g->face[f][(pos + 1) % k];
            if (sa[v][b] != NONE) { *err = KC_ERR_ORIENT; goto fail_s; }
            sa[v][b] = a; sf[v][b] = f; sp[v][b] = pos;
        }
    }
    for (int v = 0; v < nv; v++) {
        int cnt = 0, start = -1;
        for (int w = 0; w < KC_MAXV; w++) if (sa[v][w] != NONE) { if (start < 0) start = w; cnt++; }
        if (cnt == 0) goto fail_s;                       /* min() of an empty map raises in the reference */
        if (cnt > KC_MAXDEG) goto fail_s;
        int w = start;
        for (int i = 0; i < cnt; i++) {
            if (sa[v][w] == NONE) { *err = KC_ERR_ROT; goto fail_s; }
            g->rot_w[v][i] = w; g->rot_f[v][i] = sf[v][w]; g->rot_p[v][i] = sp[v][w];
            w = sa[v][w];
        }
        if (w != start) { *err = KC_ERR_ROT; goto fail_s; }
        g->deg[v] = cnt;
    }
    /* edges: the set of sorted pairs, listed in lexicographic order */
    {
        unsigned char seen[KC_MAXV][KC_MAXV];
        memset(seen, 0, sizeof seen);
        for (int f = 0; f < nfaces; f++) {
            int k = g->size[f];
            for (int i = 0; i < k; i++) {
                int u = g->face[f][i], w = g->face[f][(i + 1) % k];
                if (u > w) { int t = u; u = w; w = t; }
                seen[u][w] = 1;
            }
        }
        int e = 0;
        for (int u = 0; u < KC_MAXV; u++)
            for (int w = 0; w < KC_MAXV; w++)
                if (seen[u][w]) { g->edge[e][0] = u; g->edge[e][1] = w; e++; }
        g->nedges = e;
    }
    free(sa); free(sf); free(sp);
    if (g->n - g->nedges + g->nfaces != 2) { *err = KC_ERR_SPHERE; free(g); return NULL; }
    *err = KC_OK;
    return g;
fail_s:
    free(sa); free(sf); free(sp);
fail:
    free(g);
    return NULL;
}

void kc_graph_free(kc_graph *g) { free(g); }

/* ===================================================================== evaluation */

int kc_corners(const kc_graph *g, const double *lo, const double *hi, const double *blo, int tighten_on,
               kc_ad *C, int *ok)
{   /* model.Graph.corners.  At every valid point of the chart's region in the box [lo, hi], each corner lies
     * in its value interval with its partial derivatives in the gradient intervals; *ok = 0 only when the box
     * certainly contains no valid point (L3, L6-L8, L10'). */
    if (blo == NULL) blo = lo;
    memset(C, 0, sizeof(kc_ad) * 5 * (size_t)g->nfaces);
    kc_ad th = var(lo[0], hi[0], 0);
    kc_ad ct = cos_(th), st = sin_(th);
    kc_iv al_v = alpha_range(lo[0], hi[0]);
    int o;
    /* alpha with its exact range as value and the AD gradient of arccos(cos t / (1 + cos t)) */
    kc_ad al = acos_(div_(ct, shift(ct, 1.0)), &o);
    al.v = al_v;
    arb_t ab; arb_init(ab);
    alpha_arb(ab, blo[0]);
    double hmin = 0.5 * kc_down(ab);                     /* exact halving of a double */
    arb_clear(ab);
    double hmax = LSTRAIGHT / 2;
    int okk = 1;
    for (int fi = 0; fi < g->nfaces; fi++) {
        kc_ad *row = C + 5 * fi;
        int k = g->size[fi];
        if (k == 3) { row[0] = row[1] = row[2] = al; continue; }
        if (k == 4) {
            int i = g->pidx[fi][0];
            kc_ad A = var(lo[i], hi[i], 1);
            kc_ad O = scale(beta(ct, A), 2.0);
            /* 2 beta is decreasing in A (L7) and increasing in theta (cot beta = cos theta tan(A/2), and cos
             * theta decreases on (0, pi/2)): the exact range comes from two corners of the box */
            arb_t b; arb_init(b);
            kc_iv ex;
            beta_pt(b, lo[0], hi[i]); arb_mul_2exp_si(b, b, 1); ex.lo = kc_down(b);
            beta_pt(b, hi[0], lo[i]); arb_mul_2exp_si(b, b, 1); ex.hi = kc_up(b);
            arb_clear(b);
            kc_iv Ov;
            if (intersect(ex, O.v, &Ov) != KC_OK) return KC_ERR_DISJOINT;
            O.v = Ov;
            row[0] = A; row[1] = O; row[2] = A; row[3] = O;
            continue;
        }
        int i1 = g->pidx[fi][0], i2 = g->pidx[fi][1];
        kc_ad f1 = var(lo[i1], hi[i1], 1), f2 = var(lo[i2], hi[i2], 2);
        kc_ad cs[5];
        pentagon(ct, st, hmin, hmax, f1, f2, cs, &o); okk = okk && o;
        if (tighten_on) {
            double w0 = hi[0] - lo[0], w1 = hi[i1] - lo[i1], w2 = hi[i2] - lo[i2];
            double wmax = w0;                            /* Python max(w0, w1, w2) */
            if (w1 > wmax) wmax = w1;
            if (w2 > wmax) wmax = w2;
            if (wmax > 0.01) {
                static const int TIGHTENED[3] = { 0, 3, 4 };
                for (int t = 0; t < 3; t++) {
                    int j = TIGHTENED[t];
                    int rc = tighten(cs[j], j, lo[0], hi[0], lo[i1], hi[i1], lo[i2], hi[i2], hmin, hmax, &cs[j]);
                    if (rc != KC_OK) return rc;
                }
            }
        }
        int s;
        switch (g->mode[fi]) {
            case 'a': s = SHIFT_A; break;
            case 'b': s = SHIFT_B; break;
            case 'c': s = SHIFT_C; break;
            default: return KC_ERR_MODE;
        }
        /* chart 'c': the apex corner (formula position 3 after rotation) reaches LSTRAIGHT (L11') */
        if (g->mode[fi] == 'c') okk = okk && !(cs[3].v.hi < LSTRAIGHT);
        for (int j = 0; j < 5; j++) row[(s + j) % 5] = cs[j];
    }
    /* every corner lies in [alpha(theta), pi] (L3) */
    for (int fi = 0; fi < g->nfaces; fi++)
        for (int q = 0; q < g->size[fi]; q++) {
            kc_iv v = C[5 * fi + q].v;
            if (v.lo > KC_PI.hi || v.hi < al_v.lo) okk = 0;
        }
    *ok = okk;
    return KC_OK;
}

void kc_vertex_sums(const kc_graph *g, const kc_ad *C, kc_iv *s, kc_iv *grad)
{   /* model.Graph.vertex_sums: the sum of the corners at v minus 2 pi (L15) and its gradient in the global
     * parameters, accumulated in rotation order with outward-rounded additions */
    for (int v = 0; v < g->n; v++) {
        kc_iv acc = kc_neg(KC_TWO_PI);
        kc_iv *gv = grad + (size_t)v * g->npar;
        for (int j = 0; j < g->npar; j++) gv[j] = Z();
        for (int i = 0; i < g->deg[v]; i++) {
            int fi = g->rot_f[v][i], pos = g->rot_p[v][i];
            const kc_ad *c = &C[5 * fi + pos];
            acc = kc_add(acc, c->v);
            gv[0] = kc_add(gv[0], c->g[0]);
            if (g->pidx[fi][0] >= 0) gv[g->pidx[fi][0]] = kc_add(gv[g->pidx[fi][0]], c->g[1]);
            if (g->pidx[fi][1] >= 0) gv[g->pidx[fi][1]] = kc_add(gv[g->pidx[fi][1]], c->g[2]);
        }
        s[v] = acc;
    }
}

void kc_directions(const kc_graph *g, const kc_ad *C, kc_iv *d)
{   /* model.Graph.directions: the direction of the i-th edge at v is the sum of the corners at v before it in
     * the rotation, measured from the first neighbour (outward-rounded additions from the exact zero) */
    for (int v = 0; v < g->n; v++) {
        kc_iv acc = Z();
        for (int i = 0; i < g->deg[v]; i++) {
            d[v * KC_MAXDEG + i] = acc;
            acc = kc_add(acc, C[5 * g->rot_f[v][i] + g->rot_p[v][i]].v);
        }
    }
}

/* ===================================================================== accessors */
int kc_graph_n(const kc_graph *g) { return g->n; }
int kc_graph_npar(const kc_graph *g) { return g->npar; }
int kc_graph_nedges(const kc_graph *g) { return g->nedges; }

int kc_graph_rot(const kc_graph *g, int v, int *w, int *fi, int *pos)
{
    for (int i = 0; i < g->deg[v]; i++) { w[i] = g->rot_w[v][i]; fi[i] = g->rot_f[v][i]; pos[i] = g->rot_p[v][i]; }
    return g->deg[v];
}

void kc_graph_edges(const kc_graph *g, int *uw)
{
    for (int e = 0; e < g->nedges; e++) { uw[2 * e] = g->edge[e][0]; uw[2 * e + 1] = g->edge[e][1]; }
}
