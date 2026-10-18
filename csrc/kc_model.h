/* kc_model.h: the face-angle model of a candidate stressed graph, for the C checker kernel.
 *
 * Status: trusted.  Part of the trusted base (transfer prompt v1.2, code standards; spec v0.3 sections 4.1-4.2).
 * Implements: checklist item F3a, a bit-for-bit reimplementation of kisscrit/core/model.py (the Python reference):
 * lemmas L3 (corners in [alpha, pi]), L6 (faces are convex equilateral 3-, 4-, 5-gons), L7 (rhombus corners),
 * L8 (pentagon formulas), L10'/L11' (pentagon regions), L15 (vertex sums), the monotonicity tightening of L16;
 * paper Appendix A, Lemmas A.2 and A.6-A.8.  docs/lemmas.md carries the proofs; model.py's module docstring
 * carries the conventions (orientation, corner positions, parameter vector, pentagon charts 'a', 'b', 'c',
 * clamped functions).  Every rule is restated in kc_model.c next to the code that applies it.
 *
 * Contract.  For a graph built by kc_graph_new and a box [lo, hi] of the parameter vector p (p[0] = theta,
 * one parameter per quadrilateral, two per pentagon), kc_corners returns for every corner an interval value
 * and an interval gradient with respect to theta and the face's local parameters, and a validity flag.  These
 * are exactly the doubles model.Graph.corners returns on the same input (tests/test_cmodel.py compares the
 * bit patterns), so every enclosure property proved for the reference holds here: at every valid point of the
 * chart's region in the box, each true corner lies in its value interval, and ok = 0 only when the box
 * certainly contains no valid point.
 */
#ifndef KC_MODEL_H
#define KC_MODEL_H
#include "kc_arith.h"

#define KC_MAXV   64        /* vertices of a candidate (N = 12 graphs have at most a few dozen) */
#define KC_MAXF   128       /* faces */
#define KC_MAXDEG 32        /* vertex degree */
#define KC_MAXPAR (1 + 2 * KC_MAXF)

/* the region threshold pi - 25 degrees, the same double as model.LSTRAIGHT */
extern const double KC_LSTRAIGHT;

/* error codes: each is a case where model.py raises an exception */
#define KC_OK            0
#define KC_ERR_ORIENT    1  /* 'inconsistent orientation' (an oriented edge in two faces) */
#define KC_ERR_ROT       2  /* 'rotation does not close', or a successor missing (KeyError in the reference) */
#define KC_ERR_SPHERE    3  /* 'not a sphere' (Euler characteristic != 2) */
#define KC_ERR_LIMIT     4  /* outside the fixed sizes above, a face size other than 3, 4, 5, or an isolated
                             * vertex label; the reference has no such limits (or fails later) */
#define KC_ERR_DISJOINT  5  /* intersect(): two enclosures of the same quantity are disjoint */
#define KC_ERR_TIGHTEN   6  /* _tighten produced an empty interval */
#define KC_ERR_MODE      7  /* a pentagon without a chart mode 'a', 'b' or 'c' (KeyError in the reference) */

/* an AD value: an interval value and its interval gradient with respect to theta and the face's two local
 * parameters (model.py's 4-tuple (v, g0, g1, g2)) */
typedef struct { kc_iv v; kc_iv g[3]; } kc_ad;

/* a candidate: faces, pentagon modes, parameter indices, rotation systems, edges (model.Graph) */
typedef struct {
    int nfaces, n, npar, nedges;
    int size[KC_MAXF];
    int face[KC_MAXF][5];
    char mode[KC_MAXF];                  /* 'a', 'b', 'c' for pentagons with a mode, else 0 */
    int pidx[KC_MAXF][2];                /* global parameter indices of the face's local parameters, -1 if none */
    int deg[KC_MAXV];
    int rot_w[KC_MAXV][KC_MAXDEG];       /* rotation at v: neighbour w, face fi, position of v in that face */
    int rot_f[KC_MAXV][KC_MAXDEG];
    int rot_p[KC_MAXV][KC_MAXDEG];
    int edge[KC_MAXV * KC_MAXDEG][2];    /* sorted pairs (u < w), in lexicographic order */
} kc_graph;

/* Build a graph from nfaces oriented faces: sizes[f] in {3, 4, 5}, the vertices of all faces concatenated in
 * verts, modes[f] in {'a', 'b', 'c', 0}.  Returns NULL and sets *err on failure (model.Graph.__init__ raising).
 * Free with kc_graph_free. */
kc_graph *kc_graph_new(int nfaces, const int *sizes, const int *verts, const char *modes, int *err);
void kc_graph_free(kc_graph *g);

/* Corner values with local gradients on the box [lo, hi] (arrays of length npar); hmin is taken from blo[0]
 * (blo = NULL means lo), as model.Graph.corners(lo, hi, blo, bhi, tighten).  C has room for nfaces * 5
 * entries: C[5 * f + q] is the corner at boundary position q of face f (positions >= size are zeroed).
 * *ok receives the validity flag.  Returns KC_OK or an error code. */
int kc_corners(const kc_graph *g, const double *lo, const double *hi, const double *blo, int tighten,
               kc_ad *C, int *ok);

/* For each vertex v: s[v] encloses (sum of its corners) - 2 pi, and grad[v * npar + j] its partial derivative
 * in the global parameter j (model.Graph.vertex_sums); these are the rows of the interval Jacobian. */
void kc_vertex_sums(const kc_graph *g, const kc_ad *C, kc_iv *s, kc_iv *grad);

/* The direction angle of the edge from v to rot_w[v][i] in the frame of v (first neighbour at angle 0):
 * d[v * KC_MAXDEG + i] (model.Graph.directions). */
void kc_directions(const kc_graph *g, const kc_ad *C, kc_iv *d);

/* small accessors for the ctypes front end and the tests */
int kc_graph_n(const kc_graph *g);
int kc_graph_npar(const kc_graph *g);
int kc_graph_nedges(const kc_graph *g);
int kc_graph_rot(const kc_graph *g, int v, int *w, int *fi, int *pos);   /* returns deg(v) */
void kc_graph_edges(const kc_graph *g, int *uw);                          /* 2 * nedges ints */

#endif
