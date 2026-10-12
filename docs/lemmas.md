# Lemma inventory for the completeness computation

This document lists every fact the completeness computation relies on, with its proof, reference or certificate script and its status. It is the place where "what the run proves" is decided.

Version 2, 2026-10-11 (v2: Lemma R restated for θ < θ_max with the icosahedral exception; circumcenter reduction; test E). Update this file whenever a rule is added to the search or the checker, or a status changes.

**Status keys**
- **proved:** a short proof is written here.
- **standard:** a textbook fact; the reference is given or still needed.
- **computer:** proved by an interval computation in the checker's arithmetic (`kisscrit/core/interval.py`); the script is named.
- **numerical:** only checked numerically, so not yet part of a proof.
- **assumption:** not proved; it is listed in `STATUS.md`, and the run's conclusion is conditional on it.

**Setting.** We consider 12 points u_1..u_12 on S² with minimum angular distance θ, where θ_min ≤ θ ≤ θ_max. Here θ_min = θ(r_min) for r_min = 1.0001 (θ_min ≈ 60.0033°) and θ_max = 63.4349° (the icosahedron). Write α(θ) = arccos(cos θ / (1 + cos θ)) for the corner of the equilateral spherical triangle with side θ. It is increasing in θ, with α(60°) = 70.53° and α(θ_max) = 72°.

A configuration is **critical** if it carries a balanced stress: weights w_ij ≥ 0, not all zero, on the pairs at distance θ, such that at each point the tangential component of Σ_j w_ij (u_i − u_j) vanishes. See the paper, Lemma 2.1 and the balance criterion of KKLS. The **stressed graph** G is the support of one such stress; its edges are geodesic arcs of length θ and every edge has positive weight.

## A. Structure of the stressed graph

**L1. G is a subgraph of the convex-hull graph** (proved). Let u_i, u_j be at the minimum distance θ. The open spherical cap centred at their midpoint with angular radius θ/2 contains no point of the configuration. A point in it would be at distance < θ/2 + θ/2 = θ from u_i or from u_j, by the triangle inequality, and strictly less by openness. A pair with an empty circumscribed cap is an edge of the spherical Delaunay triangulation, which is the radial projection of the convex hull of the points. So [u_i u_j] is an edge of the hull, or a diagonal of a non-triangular hull face; in the second case we triangulate that face so that it is an edge of the triangulation T.

**L2. G is planar** (proved). The edges are hull edges by L1, and hull edges project radially to non-crossing geodesic arcs.

**L3. All corners of G lie in [α(θ), π]** (proved).
- *Upper bound:* at a vertex v, the stressed edges at v carry positive weights and their tangential forces sum to zero. So the edge directions are not contained in an open half-plane of T_vS², and every angle between consecutive edges (a corner) is ≤ π.
- *Lower bound:* the two neighbours a, b bounding a corner are at distance ≥ θ from each other, and both are at distance θ from v. In the spherical triangle v a b, the angle at v is therefore ≥ the angle of the equilateral triangle with side θ, which is α(θ), by the spherical law of cosines (the angle at v increases with the opposite side).

**L4. Degrees are between 2 and 5, and a degree-2 vertex has two straight corners** (proved).
- *Degree ≤ 5:* corners are ≥ α(θ) > 70° and sum to 2π, so there are at most 5.
- *Degree ≥ 2:* a vertex of degree 1 cannot balance a single nonzero force.
- *Degree 2:* two corners ≤ π (L3) that sum to 2π are both π.
- *Vertices of degree 0* are rattlers; see L9.

**L5. G is connected and has no cut vertex** (proved, given the perimeter bound of L6).
- *No cut vertex:* let v be a cut vertex and S a component of G − v. For every edge, the torques of the two forces about the centre cancel, and the normal reactions k_i u_i have no torque. So the total torque on S is −(Σ_{i∈S} w_iv u_i) × u_v = 0. Hence the tangential forces from S at v balance by themselves, and the edges from S at v are not contained in an open half-plane. The same holds for the edges at v from the rest of G. But in a plane graph the edges at v to different components of G − v occupy disjoint angular sectors, and two disjoint sectors cannot both be wider than π. Contradiction.
- *Connected:* a component carrying a balanced positive stress is not contained in an open hemisphere. Look at its point that is extreme in the direction of the hemisphere's pole: all its stressed forces have a positive tangential component towards the pole, so they do not balance. But a second component would lie inside one face of the first, and faces are convex polygons of perimeter < 2π (L6), which lie in an open hemisphere.
- *Conclusion:* G is 2-connected, so every face is an open disk bounded by a simple cycle.

**L6. Faces are convex equilateral spherical polygons with 3, 4 or 5 corners** (proved, given a standard fact).
- *Convexity:* by L3, every corner of a face is ≤ π. A simple spherical polygon whose interior angles are all ≤ π, and which is not a lune or a hemisphere, is convex: it is the intersection of the hemispheres bounded by its edges and lies in an open hemisphere. **Reference needed:** the spherical analogue of the planar fact that a locally convex simple polygon is convex. A direct proof via the tangent-line argument on the hemisphere determined by any edge is planned.
- *Perimeter < 2π:* a convex spherical polygon in an open hemisphere has perimeter < 2π (standard; for example, it is contained in the great circle's hemisphere and its perimeter is less than that of the boundary circle). So kθ < 2π, and θ > π/3 gives k ≤ 5.
- *Triangulation:* a face contains no configuration point in its interior, by L9, so T triangulates it by diagonals. A k-gon face is therefore a union of k − 2 triangles of T, forming a path in the dual graph.

**L7. Rhombus faces have corners in [α(θ), 2α(θ)]** (proved). An equilateral spherical quadrilateral with side θ is a rhombus. Its corners are A, O, A, O with O = 2β(θ, A), where β(θ, A) = arctan(cos(A/2) / (cos θ sin(A/2))) is the base angle of the isosceles triangle with legs θ and apex A. 2β decreases in A. At A = 2α(θ), the rhombus is two equilateral triangles and O = α(θ). By L3, O ≥ α(θ), so A ≤ 2α(θ); A ≥ α(θ) is L3 again.

**L8. Pentagon faces: formulas and validity** (proved).
- *Formulas:* with corners f1, f2 at P1, P2, the corner functions of `kisscrit/core/model.py:pentagon` follow from the spherical law of cosines and the four-part (cotangent) formula, applied to the triangles P0P1P2, P0P2P3 and P0P3P4. Checked against direct construction on 707 random pentagons: `scripts/completeness/` pilot, error < 1e-12. A written derivation is in the docstring and is to be expanded in the appendix.
- *Validity:* every diagonal is ≥ θ, so in particular P0P3 ≥ θ, and every corner is ≥ α(θ) by L3.

**L9. No rattlers** (proved + computer: `scripts/completeness/lemmas.py R`). A rattler is a configuration point that is not a vertex of G. It lies in the interior of a face K of G, at distance ≥ θ from every vertex of K. The face is a convex equilateral triangle, rhombus or pentagon with side θ ∈ [60°, θ_max], corners in [α(θ), π] and diagonals ≥ θ (L3, L6–L8).

**Lemma R.** For θ < θ_max, no point of such a face is at distance ≥ θ from all its vertices. For θ = θ_max, the only exception is the centre of the regular pentagon, which is at distance exactly θ_max from its five vertices; that is the icosahedron.

**Consequence.** For θ < θ_max there are no rattlers. At θ = θ_max the configuration is the icosahedron, by the uniqueness of the optimal Tammes configuration for N = 12 (Fejes Tóth; **reference to be checked**), which is in the census.

**The bound is sharp at θ_max.** The five neighbours of an icosahedron vertex form a regular pentagon with side θ_max, and the vertex is at distance θ_max from all of them. So the margin of Lemma R tends to 0 as θ → θ_max, and a purely numerical test cannot reach θ_max. Test E below handles this.

*Reduction to circumcenters* (proved). Let f(x) = min over vertices v of d(x, v), and suppose max_K f ≥ θ, attained at x*.
1. x* is not on the boundary: a point of an edge of length θ is within θ/2 of an endpoint.
2. f(x*) < 90°: otherwise every vertex lies in the closed hemisphere {y : ⟨y, x*⟩ ≤ 0}, which is convex, so it contains K and hence x*, which is absurd.
3. So each d(·, v) is smooth near x* (x* ≠ ±v), and x* is an interior local maximum of their minimum over the active vertices A (those at distance f(x*)). By the first-order condition for a maximum of a minimum of smooth functions, 0 is in the convex hull of the gradients, which are the unit tangent vectors at x* pointing away from the vertices in A.
4. |A| = 1 is impossible (a single nonzero gradient). If |A| = 2, the two vertices u, w lie in opposite directions, so x* is the midpoint of a geodesic arc uw of length 2f(x*) < π. Moving from x* perpendicular to that arc increases both distances (cos d = cos f · cos t), so x* is not a local maximum.
5. Hence |A| ≥ 3: x* is equidistant, at distance R = f(x*) ∈ [θ, 90°), from three vertices a, b, c. So x* = c_abc, their circumcenter on the side where the circumradius is below 90°.
6. Conclusion: it suffices to show, for every triple of vertices and both normals ±n of the plane through them, that the corresponding circumcenter c satisfies one of:
   - A: its circumradius is ≥ 90°;
   - B: its circumradius is < θ;
   - C: it lies outside K;
   - D: another vertex is strictly closer.

   These are the tests A–D of `no_deep_hole` and `triple_excluded`. All of them are sign conditions on polynomials in the vertex coordinates.

*Test E (pentagons, near θ_max)* (proved). Suppose x* is a deep hole, as above, and all five vertices are within 70° of it.
- The triangles x*P_iP_{i+1} tile K, so their angles X_i at x* sum to 2π.
- The legs a = d(x*, P_i) and b = d(x*, P_{i+1}) lie in [θ, 70°], and the opposite side is θ. So cos X_i = F(a, b) = (cos θ − cos a cos b)/(sin a sin b).
- The sign of ∂F/∂a is the sign of cos b − cos a cos θ ≥ cos 70° − cos²60° > 0.34 − 0.25 > 0, and the same holds for b. So F(a, b) ≥ F(θ, θ) = cos θ/(1 + cos θ) = cos α(θ), with equality only at a = b = θ.
- Hence X_i ≤ α(θ) ≤ 72°, and 2π = ΣX_i ≤ 5α(θ) ≤ 2π. So θ = θ_max and every distance equals θ.

The script proves "all five vertices within 70°" as (n·P_m)² > cos²70°·|n|² with the correct sign. This test removes the neighbourhood of the regular θ_max pentagon, and it is the reason the statement of Lemma R distinguishes θ = θ_max. It is used only in the pentagon chart, where the threshold equals the side length (the sharpness test in `tests/test_lemmas.py` uses other thresholds and so calls `no_deep_hole` without E).

*Charts.*
- Triangles: one parameter θ.
- Rhombi: (θ, A), as in L7.
- Pentagons: the vertex chart (θ, c0, f1), with corners c0 at P0 and f1 at P1 in [α(θ0), π].
  - P0 is the north pole, P1 is at bearing 0, P4 is at bearing c0, and P2 is obtained by rotating P0 about P1 by −f1.
  - P3 is the apex over the diagonal P2P4, on the far side from P0: P3 = λ(P2 + P4) + μ(P4 × P2), with λ = cos θ/(1 + c), μ² = (1 − 2cos²θ/(1 + c))/(1 − c²) and c = ⟨P2, P4⟩.
  - Lemma R is invariant under cyclic relabelling. By Lemma P (L10), at most two corners are ≥ 155°, so we may take the corner at P3 to be ≤ 155°. This keeps the construction away from the singularity of μ (μ = 0 at a straight corner at P3).
  - Boxes are discarded when some turn det(P_{i−1}, P_i, P_{i+1}) is certainly negative (a corner > π) or some short diagonal is certainly < θ.
  - All quantities are enclosed in mean-value form: the value on the box, intersected with the value at the midpoint plus the interval gradient times the box radius.
- *Why not the (θ, f1, f2) chart:* the apex corner there has a square-root singularity, and plain interval evaluation of the chained rotations overestimates by a factor of about 50. A first version that searched over the point x as well (5 parameters) did not finish.

*Not vacuous* (`tests/test_lemmas.py`):
- the constructed polygons are equilateral, close up and are convex;
- the centroid passes every inside test;
- the test fails just below the true covering radius and passes just above it.

*Numerical margin.* At θ = 60° the largest covering radius of a valid pentagon is about θ − 3.9°. It tends to θ as θ → θ_max, at the regular pentagon (see above). A coarse grid at θ = 63.435° suggested a margin of 2.7°; that grid missed the regular pentagon. Triangles and rhombi have large margins (4α(θ) < 2π, so test E's argument alone excludes them).

## B. Pentagon charts

**L10. Lemma P** (computer: `scripts/completeness/lemmas.py P`, proved in 1319 boxes). No two adjacent corners of a face pentagon (θ ∈ [60°, 63.435°], corners ≥ α(θ), diagonal P0P3 ≥ θ) are both ≥ π − 25°. Together with the definitions of regions a, b and c (`docs/certificate-format.md`, `kisscrit/core/model.py`), this shows the three regions cover all pentagons.

**L11. A corner ≥ π − 25° sits at a vertex of degree ≤ 3** (proved). Otherwise the other three or more corners at that vertex are each ≥ α(θ) > 70.52°, and 155° + 3·70.52° > 360°. This is used to skip regions b and c (`kisscrit/completeness/runner.py:charts_of`).

## C. Candidate generation (`kisscrit/completeness/candidates.py`)

**L12. Enumeration of triangulations** (standard + computer). The 7595 combinatorial triangulations of the sphere with 12 vertices are generated by edge flips from the icosahedron (`data/completeness/triangulations12.json`).
- The flip graph is connected (Wagner 1936).
- The count agrees with OEIS A000109.
- The canonical-form deduplication is exact (`kisscrit/completeness/canon.py`).

**L13. Tilings** (proved, from L1, L4, L6). Every stressed graph arises from some T in L12 by merging triangles into faces of 1, 2 or 3 triangles. The vertex pruning in `tilings()` uses only:
- the corner ranges of L3 and L7;
- that triangle corners equal α(θ);
- that degrees are 2..5 with straight corners at degree 2 (L4), which can only occur in pentagons, since triangle corners are α < π and rhombus corners are ≤ 2α < π;
- that a pentagon is not a fan of three triangles around a vertex (which would be interior, contradicting L9).

**L14. Linear relaxation filter** (proved, from L3, L7, L13). The filter uses only:
- vertex sums of 360°;
- triangle corners equal to α;
- rhombus corners in [α, 2α] and pentagon corners in [α, 180°];
- that a k-gon has corner sum > (k − 2)·180°, by positivity of area (Gauss–Bonnet).

## D. The branch and bound

**L15. Necessary conditions** (proved). Every critical configuration with stressed graph G gives a parameter point in one of the pentagon regions of `charts_of(oriented(G))` at which:
- the vertex sums are 2π;
- all corners satisfy L3;
- the face conditions L7 and L8 hold;
- the edge directions admit a strictly positive balanced stress.

Only these necessary conditions are used to exclude boxes. Realizability, which needs the cone-metric argument, is not needed for exclusion.

**L16. Exclusion tests are sound** (proved; see `kisscrit/core/checker.py`).
- *Mean-value inclusion:* used for the clamped corner functions, which are Lipschitz with interval Clarke Jacobian enclosures.
- *Gordan alternative:* Ax > 0 solvable ⟺ there is no nonzero y ≥ 0 with Aᵀy = 0 (standard).
- *Motzkin alternative:* for strictly positive stresses (standard).
- *Monotonicity tightening:* a sign of a partial-derivative enclosure implies monotonicity on the box.

## E. Assumptions (listed in `STATUS.md`)

- **H (hull conjecture, part a).** For θ > 60°, the hull triangulation is one of T0, T1, T5, T7. This is the conditional route. It is checklist item C1; until it is proved, the run's conclusion is conditional on it.
- **R0 (range).** Only r ≥ r_min = 1.0001 is covered (checklist C5).

## Status of the computations

| lemma | script | status |
|---|---|---|
| P | `lemmas.py P` | proved (1319 boxes, about 1 s) |
| R, triangles | `lemmas.py R` | proved (1 box) |
| R, rhombi | `lemmas.py R` | proved (1061 boxes, about 1 s) |
| R, pentagons | `lemmas.py R` (or `R 0 2` and `R 1 2` in parallel) | proved (25409 + 32997 boxes, about 27 CPU-minutes) |

## Still to do

- **L6:** a reference or written proof for "a locally convex simple spherical polygon is convex".
- **L8:** a written derivation of the pentagon formulas, for the appendix.
- **L9:** a reference for the uniqueness of the icosahedron as the optimal Tammes configuration for N = 12.
