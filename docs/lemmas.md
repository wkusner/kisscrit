# Lemma inventory for the completeness computation

This document lists every fact the completeness computation relies on, with its proof, reference or certificate script and its status. It is the place where "what the run proves" is decided.

Version 3, 2026-10-12 (v3: L6 proved, L8 derived, uniqueness at θ_max via Lemma U; v2: Lemma R restated for θ < θ_max with the icosahedral exception; circumcenter reduction; test E). Update this file whenever a rule is added to the search or the checker, or a status changes.

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

**L6. Faces are convex equilateral spherical polygons with 3, 4 or 5 corners** (proved).

Let K be the closure of a face of G. By L5, K is a closed disk whose boundary γ is a simple geodesic polygon v_1 … v_k.
- Its edges e_i = [v_i v_{i+1}] have length θ < π.
- Its interior angles are β_i ∈ [α(θ), π] (L3).
- Orient γ with K on the left. Let n_i be the unit normal of the great circle C_i ⊃ e_i on the side of K, and H_i = {x : ⟨x, n_i⟩ ≥ 0}.

**Claim.** K = H_1 ∩ … ∩ H_k, and K lies in the open hemisphere {x : ⟨x, q⟩ > 0} with q = n_1 + … + n_k.

*Step 1 (the normals lie in an open hemisphere).*
- Gauss–Bonnet for the disk K gives Area(K) + Σ(π − β_i) = 2π. Since Area(K) > 0, the total turning τ = Σ(π − β_i) is < 2π.
- Consecutive normals n_{i−1}, n_i are both orthogonal to v_i, at angle π − β_i (the turning angle at v_i). So the closed polygon N = n_1 n_2 … n_k has length τ < 2π.
- A closed curve of length L < 2π on S² lies in an open hemisphere (the lemma in the proof of Fenchel's theorem):
  - Split it at two points q_1, q_2 into arcs of length L/2 each. Then d(q_1, q_2) ≤ L/2 < π, so the midpoint m is defined.
  - The rotation by π about m swaps q_1 and q_2 and sends x ⊥ m to −x. So for x on the great circle m^⊥, d(q_1, x) + d(x, q_2) = d(q_1, x) + d(−x, q_1) = π.
  - An arc through such an x would therefore have length ≥ π > L/2. So the curve misses m^⊥. It is connected and contains q_1 with ⟨q_1, m⟩ > 0, so it lies in {⟨·, m⟩ > 0}.
- Put p = m. Then ⟨p, n_i⟩ > 0 for every i, so neither p nor −p lies on any C_i, and γ avoids ±p.

*Step 2 (K is star-shaped from p).*
- Use polar coordinates (r, φ) about p, with r = d(p, x) ∈ (0, π) and φ the bearing.
- A great circle C not through ±p, traversed with p on its left, has strictly increasing bearing, and an arc of it of length < π changes φ by less than π.
- So φ increases strictly along γ, by 2πw in total with w ≥ 1. A simple closed curve in the annulus S² ∖ {±p} has winding number 0 or ±1 about it, so w = 1. Hence every open meridian from p to −p meets γ exactly once, and the two components of S² ∖ γ are the parts of the meridians before and after γ.
- At a point x of e_1, the geodesic from x towards p starts with ⟨·, n_1⟩ increasing, because ⟨p, n_1⟩ > 0. So it enters K. Hence K is the component containing p: K = {r ≤ r_γ(φ)}.

*Step 3 (support functions).*
- Put u = cot r. In the coordinates (φ, u):
  - a great circle C with pole n = cos δ · p + sin δ · e(φ_0), where 0 ≤ δ < π/2, is the graph u = L(φ) = −tan δ · cos(φ − φ_0), a solution of L″ + L = 0;
  - its p-side hemisphere {⟨x, n⟩ ≥ 0} is {u ≥ L(φ)}.
- So K = {u ≥ g(φ)}, where the 2π-periodic function g equals L_i on the φ-interval I_i of the edge e_i, and H_i = {u ≥ L_i}.
- At v_i the interior angle is ≤ π, so near v_i the set K is the wedge H_{i−1} ∩ H_i. Hence g = max(L_{i−1}, L_i) near the break point, and g′ jumps up by some κ_i ≥ 0 there.
- Fix j and put w = g − L_j. Then w″ + w = Σ_i κ_i δ_{φ_i} ≥ 0, and w ≡ 0 on I_j = [φ_s, φ_e]. Variation of constants gives:
  - w(φ) = Σ_{φ_e < φ_i ≤ φ} κ_i sin(φ − φ_i) ≥ 0 for φ ∈ [φ_e, φ_e + π];
  - w(φ) = Σ_{φ ≤ φ_i < φ_s} κ_i sin(φ_i − φ) ≥ 0 for φ ∈ [φ_s − π, φ_s].
- Since I_j has positive length, these two ranges cover the circle. So g ≥ L_j everywhere, i.e. K ⊆ H_j.
- As g equals some L_i at every φ, g = max_i L_i and K = ∩ H_i.

*Step 4 (open hemisphere).*
- If the n_i span R³, then for x ∈ K ∖ {0} we have ⟨x, q⟩ = Σ⟨x, n_i⟩ ≥ 0, with equality only if x ⊥ every n_i, i.e. x = 0. So K ⊂ {⟨·, q⟩ > 0}.
- Otherwise all the C_i pass through a pair ±c, and K = ∩ H_i is a lune or a hemisphere, with boundary length 2π = kθ. That is impossible here, because 2π/θ ∈ (5.67, 6) is not an integer for θ ∈ (60°, θ_max].

*Perimeter.*
- K is compact in the open hemisphere, so K ⊂ a cap B of radius ρ < π/2 about q/|q|.
- Cut B successively along the great circles C_1, …, C_k, keeping the side containing K. Each cut replaces a boundary arc by the minimizing geodesic between its endpoints, so the perimeter does not increase, and the result is K.
- Hence kθ = per(K) ≤ 2π sin ρ < 2π. With θ > π/3 this gives k ≤ 5.

*Triangulation.* For θ < θ_max, a face contains no configuration point in its interior (L9), so T triangulates it by diagonals. A k-gon face is a union of k − 2 triangles of T, forming a path in the dual graph. At θ = θ_max the configuration is the icosahedron (L9).

**L7. Rhombus faces have corners in [α(θ), 2α(θ)]** (proved). An equilateral spherical quadrilateral with side θ is a rhombus. Its corners are A, O, A, O with O = 2β(θ, A), where β(θ, A) = arctan(cos(A/2) / (cos θ sin(A/2))) is the base angle of the isosceles triangle with legs θ and apex A. 2β decreases in A. At A = 2α(θ), the rhombus is two equilateral triangles and O = α(θ). By L3, O ≥ α(θ), so A ≤ 2α(θ); A ≥ α(θ) is L3 again.

**L8. Pentagon faces: formulas and validity** (proved).

Let P0 … P4 be a convex equilateral pentagon with side θ and corners f1, f2 at P1, P2. By L6 it is convex, so the diagonals P0P2 and P0P3 lie inside it and cut it into the triangles P0P1P2, P0P2P3 and P0P3P4. Its corners are sums of the angles of these triangles. The code is `kisscrit/core/model.py:pentagon` (checker; the search side, `kisscrit/completeness/fast.py:pentagon`, is untrusted and plays no part in the proof). Write β(θ, A) = arctan(cos(A/2)/(cos θ sin(A/2))).

1. *Triangle P0P1P2* (legs θ, apex f1). Halve it into right triangles.
   - The base d2 = P0P2 satisfies sin(d2/2) = sin θ sin(f1/2) (Napier: the sine of a side is the sine of the hypotenuse times the sine of the opposite angle).
   - The base angles satisfy cos θ = cot(f1/2) cot b1, i.e. b1 = β(θ, f1).
2. *Angle at P2 in P0P2P3.* Convexity puts P2→P0 inside the corner f2 at P2, so the angle between P2→P0 and P2→P3 is g = f2 − b1.
3. *Triangle P0P2P3* (sides d2, θ with included angle g).
   - Law of cosines: cos d3 = cos d2 cos θ + sin d2 sin θ cos g, with d3 = P0P3.
   - The four-part formula, cot A sin C = cot a sin b − cos b cos C, applied with C = g:
     - at P0 (a = θ, b = d2): cot a0 = (cos θ sin d2 − sin θ cos d2 cos g)/(sin θ sin g);
     - at P3 (a = d2, b = θ): cot a3 = (sin θ cos d2 − cos θ sin d2 cos g)/(sin d2 sin g).

     The code writes these as a = π/2 − arctan(cot a), which is valid for a ∈ (0, π).
4. *Triangle P0P3P4* (legs θ, base d3).
   - The apex is f4 = 2h with sin h = sin(d3/2)/sin θ, h ∈ (0, π/2], since f4 ≤ π.
   - The base angles are b3 = β(θ, f4).
5. *Corners:* c0 = b1 + a0 + b3, c1 = f1, c2 = f2, c3 = a3 + b3, c4 = f4.

*Validity:*
- The triangle P0P2P3 is non-degenerate exactly when g ∈ (0, π); the code checks this (clip and ok flag).
- The pentagon exists exactly when d3 ≤ 2θ, i.e. when the asin argument is ≤ 1.
- Every diagonal is ≥ θ, so in particular d3 ≥ θ.
- Every corner is ≥ α(θ), by L3.
- Given θ, f1 and f2, the pentagon is determined up to isometry, so the chart (θ, f1, f2) covers all pentagons. In regions b and c it is used after a cyclic relabelling (L10).

*Check:* the formulas agree with direct construction (`polygon_vertices` in `scripts/completeness/lemmas.py`, and `tests/test_lemmas.py`) to 1e-12. The Lemma R pentagon chart builds the vertices directly instead (L9).

**L9. No rattlers** (proved + computer: `scripts/completeness/lemmas.py R`). A rattler is a configuration point that is not a vertex of G. It lies in the interior of a face K of G, at distance ≥ θ from every vertex of K. The face is a convex equilateral triangle, rhombus or pentagon with side θ ∈ [60°, θ_max], corners in [α(θ), π] and diagonals ≥ θ (L3, L6–L8).

**Lemma R.** For θ < θ_max, no point of such a face is at distance ≥ θ from all its vertices. For θ = θ_max, the only exception is the centre of the regular pentagon, which is at distance exactly θ_max from its five vertices; that is the icosahedron.

**Consequence.** For θ < θ_max there are no rattlers. At θ = θ_max the configuration is the icosahedron (below), which is in the census.

**Uniqueness at θ_max** (proved + computer: `scripts/completeness/lemmas.py U`). Twelve points with minimum distance θ_max = arctan 2 form a regular icosahedron. This is Fejes Tóth's argument [Fejes Tóth 1943], with the equality case made explicit.
1. Every hull face has an empty open circumcap of radius R. Its centre is at distance ≥ R from all twelve points. If R ≥ 60°, adding the centre would give 13 points with pairwise distances ≥ 60°, contradicting the kissing number 12 [Schütte–van der Waerden 1953]. So R < 60°. Triangulating any non-triangular hull face keeps the same circumcircle.
2. So the hull triangulation consists of 2·12 − 4 = 20 spherical triangles, with sides ≥ θ_max and circumradius < 60° (< 90°, so they tile the sphere). Their areas sum to 4π.
3. **Lemma U** (computer): such a triangle has area ≥ π/5 = 3·72° − 180°, with equality only for the equilateral triangle of side θ_max. Near that corner, the box is closed by monotonicity (all three partial derivatives of the area are positive).
4. Since 20 · π/5 = 4π, every triangle is equilateral with side θ_max and angles 72°. So every vertex has degree 5. The only 5-regular triangulation with 12 vertices is the icosahedron's, and developing from one triangle (each neighbour is its reflection across the shared edge) fixes the configuration up to isometry.

*The circumradius bound matters.* With circumradius θ_max, the isosceles triangle (θ_max, θ_max, 116.57°), half of the regular-pentagon fan, also has area exactly π/5. `tests/test_lemmas.py` checks that Lemma U fails on a box around it when the bound is raised to 64°.

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

**L10. Lemma P** (computer: `scripts/completeness/lemmas.py P`, proved in 1319 boxes). No two adjacent corners of a face pentagon (θ ∈ [60°, 63.435°], corners ≥ α(θ), diagonal P0P3 ≥ θ) are both ≥ L, where **L is the double `LSTRAIGHT` = `math.pi - math.radians(25.0)`** (about 155°). The computation tests the corners, as intervals, against this double, and every region inequality below is against the same double; the real number π − 25° never enters.

**L10′. The three pentagon regions and their cover** (proved, from L10 and L3). Label the corners of a face pentagon c0 … c4 in boundary order (the order of the face cycle in the certificate). The regions are

| region | defining inequalities | parameters (`p` entries) | clamped apex corner |
|---|---|---|---|
| a | c4 ≤ L | c1, c2 | c4 (formula position 4 = boundary 4) |
| b | c4 ≥ L and c1 ≤ L | c3, c4 | c1 (boundary 1) |
| c | c4 ≥ L and c1 ≥ L | c4, c0 | c2 (boundary 2) |

In the code, the "≥ L" inequalities are lower bounds of the a-priori box (`kisscrit/core/checker.py:apriori_box`: the parameter at boundary position 4 in regions b and c has lower end L), the "c1 ≥ L" of region c is the box condition `kisscrit/core/model.py:Graph.corners` (`mode == 'c'`: invalid if c1 is certainly < L), and the "≤ L" inequalities are the clamp of the apex corner at `LSTRAIGHT` in `model.py:pentagon` (through `hmax`).

*Cover.* Let P be a valid face pentagon. If c4 ≤ L it lies in region a. Otherwise c4 ≥ L; if c1 ≤ L it lies in region b, and if c1 ≥ L it lies in region c. So every valid pentagon lies in at least one region.

*The clamp is harmless.* In region a the clamped corner c4 is ≤ L by definition; in region b, c1 ≤ L by definition; in region c, c1 ≥ L and c2 is adjacent to c1, so c2 < L by Lemma P. So on each region the clamped pentagon formulas agree with the true corner functions, and `ok = False` (a box certainly containing no point of the region that is a valid pentagon) is a valid exclusion for that region's chart: a configuration whose parameters lie in the box but outside the region is covered by the chart of another region, and every required region's chart (L11′) covers that region's whole a-priori box.

*Degree-2 vertices.* By L4 a corner at a degree-2 vertex is exactly π. If the vertex at boundary position 4 has degree 2, region a is empty (c4 = π > L) and is not required; the coordinate for c4 is then the point π, which the certificate must store as an interval containing π (`[π⁻, π⁺]` in doubles), and the checker's a-priori box for it is `[down(π), up(π)]`.

**L11. A corner ≥ L sits at a vertex of degree ≤ 3** (proved). Otherwise the other three or more corners at that vertex are each ≥ α(θ) > 70.52°, and 155° + 3·70.52° > 360°.

**L11′. Required regions** (proved, from L10′, L11, L4). Let d4, d1 be the degrees of the vertices at boundary positions 4 and 1. If d4 ≥ 4 then c4 < L (L11), so region a alone covers; if d4 = 3 and d1 ≥ 4 then c1 < L, so regions a and b cover; if d4 = 3 and d1 ≤ 3 all three are needed. If d4 = 2 (c4 = π, L4) region a is empty: if d1 ≥ 4 region b alone covers, if d1 = 2 region c alone (c1 = π), otherwise b and c. The checker recomputes this list from the certificate's faces (`kisscrit/core/checker.py:required_modes`) and requires every listed chart to be present and complete (V2). The list may include a region that is in fact empty (for example d4 = 3, d1 = 2 lists region b although c1 = π ≥ L); an extra chart only adds work, never removes a case. The search's own list (`kisscrit/completeness/runner.py:charts_of`, untrusted) is the same rule.

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

**L15. Necessary conditions** (proved). Every critical configuration with stressed graph G gives a parameter point in one of the required pentagon regions of G (L11′, `kisscrit/core/checker.py:required_modes`) at which:
- the vertex sums are 2π;
- all corners satisfy L3;
- the face conditions L7 and L8 hold;
- the edge directions admit a strictly positive balanced stress.

Only these necessary conditions are used to exclude boxes. Realizability, which needs the cone-metric argument, is not needed for exclusion.

**L16. Exclusion tests are sound** (proved; see `kisscrit/core/checker.py`).
- *Mean-value inclusion:* for a function F that is Lipschitz on a box B with an interval enclosure J(B) of its Clarke generalized Jacobian, F(p) ∈ F(m) + J(B)(p − m) for all p, m ∈ B (Lebourg's mean value theorem). The clamped corner functions of `kisscrit/core/model.py` are Lipschitz on every box except at the degenerate boundary of the pentagon formulas where cos d3 → 1 (the diagonal P0P3 collapses); there the code's derivative enclosure is (−∞, ∞), the inclusion is vacuous, and the linear test proves nothing (it never accepts on an unbounded enclosure). Infinite endpoints in the arithmetic always mean "unbounded range of finite reals".
- *Gordan alternative:* Ax > 0 solvable ⟺ there is no nonzero y ≥ 0 with Aᵀy = 0 (standard).
- *Motzkin alternative:* for strictly positive stresses (standard). The checker demands every nonzero edge term to be > 0 on the box, which is a strengthening of the alternative (≥ 0 with one > 0 would do); sound, and simpler.
- *Edge directions and frames.* At each vertex v the checker places the first neighbour of v in the rotation system at direction 0 and the following neighbours at the cumulative corner sums (`model.py:Graph.directions`). By L2 and L5 the rotation system of the realization at v is the combinatorial one, so these are the true directions of the edges at v in some orthonormal frame of T_v S². The balance condition "Σ_e w_e t_v(e) = 0 at every v" is invariant under the choice of frame at each vertex, and so are the Gordan and Motzkin alternatives, since y_v ranges over all of R².
- *Monotonicity tightening:* if a partial derivative enclosure of a Lipschitz function has a sign on the box, the function is monotone in that variable there (Lebourg again), so its extreme values over the box are attained on the corresponding faces of the box. The tightening evaluates the clamped formulas on those faces; those evaluations enclose the clamped function whether or not the face contains valid pentagons, so the tightened interval encloses the corner on every valid point of the box.

**L17. Domain of the model's enclosures** (proved). Every enclosure of `kisscrit/core/model.py` and every bound of the V1 a-priori box is valid when θ ∈ [THETA_FLOOR, THETA_CEIL] = [1.0, 1.5] (radians; `model.py:THETA_FLOOR`, `THETA_CEIL`, both exact doubles, 0 < 1.0 < 1.5 < π/2). The claimed range [θ_min, θ_max] ≈ [1.04726, 1.10715] lies inside. The checker's V3 requires the certificate's θ_lo and θ_hi to be finite, θ_lo < θ_hi, and THETA_FLOOR ≤ θ_lo, θ_hi ≤ THETA_CEIL, decided in Arb on the certificate's doubles (`kisscrit/core/checker.py:theta_domain_ok`).

*Facts used, each true on (0, π/2) and hence on the domain.*
- (D1) *Positivity.* sin θ > 0, cos θ > 0 and 1 + cos θ > 0. These are the denominators of α′, of β (cos θ), and of the pentagon formulas (sin θ, in a0, a3 and the apex asin argument sin(d3/2)/sin θ).
- (D2) *α is increasing.* With c = cos θ ∈ (0, 1), c/(1 + c) = 1 − 1/(1 + c) increases with c, and c decreases with θ; arccos is decreasing. So α(θ) = arccos(c/(1 + c)) increases with θ, with values in (π/3, π/2).
- (D3) *β is increasing in θ.* For A ∈ (0, π), β(θ, A) = arctan(cot(A/2)/cos θ) with cot(A/2) > 0, and 1/cos θ increases with θ, so β increases with θ. (The rhombus range also uses that β decreases in A, which holds for A ∈ (0, 2π) because cot(A/2) decreases there and cos θ > 0 (D1). Both monotonicities together need A ∈ (0, π] on the box: a condition on the rhombus coordinate, not on θ; see below.)
- (D4) *The clamp bound.* `Graph.corners` clamps the half apex corner h = f4/2 of a pentagon below at hmin = α(blo[0])/2, with blo[0] the lower θ end of the enclosing box. At a valid point with θ ≥ blo[0], f4 ≥ α(θ) ≥ α(blo[0]) by L3 and D2, so the clamp does not act there. Also hmin < π/4 < hmax = LSTRAIGHT/2, so [hmin, hmax] lies inside [−π/2, π/2), as `asin_angle` requires.
- (D5) *The d2 clamp.* sin(d2/2) = sin θ sin(f1/2) ≤ sin θ, so d2/2 ≤ θ ≤ 1.5 < π/2 (rounded down) and the clamp of `pentagon` at π/2 does not act at valid points.

*Where the domain is used.*

| place | fact | covered by [1.0, 1.5] |
|---|---|---|
| `model.py:alpha_range` (exact range from the two ends) | D2 | yes |
| `model.py:Graph.corners`, gradient of α (`acos_(ct/(1+ct))`) | D1 | yes |
| `model.py:Graph.corners`, `hmin = α(blo[0])/2` | D2, D4 | yes |
| `model.py:Graph.corners`, rhombus range `ex` from `_beta_pt` at two corners | D3 (in θ) | yes for θ; also needs A ∈ (0, π] on the box |
| `model.py:beta`, `_beta_pt` (cos θ in the denominator) | D1 | yes |
| `model.py:pentagon` (sin θ in denominators; d2 clamp) | D1, D5 | yes |
| `model.py:Graph.corners`, "every corner ≥ α(θ)" against `alpha_range` | D2 | yes |
| `checker.py:apriori_box`, a_lo = α(θ_lo), rhombus hi = 2α(θ_hi) (L3, L7) | D2 | yes |

*Why this is needed.* Outside (0, π/2) the facts fail: α is even in θ, so for θ_lo < 0, α(θ_lo) can exceed α(θ_hi) and `alpha_range` is inverted; for θ_hi > π/2, α and β are not monotone. Before V3 was two-sided, a certificate with θ = [−5, θ_hi] or [1.0471975511965976, 7] and one excluded leaf was accepted for a candidate with a census critical configuration (computed, 2026-10-21, V3-FIX session; `tests/test_certificates.py::test_v3_theta_outside_model_domain`).

*Not covered by this entry.* The certificate's boxes for the corner coordinates are bounded only from inside by V1; D3 needs the rhombus coordinate of the box in (0, π]. This is L32 below: a rhombus start box outside (0, π] was accepted at 5da4bba (a forged certificate, reproduced by the DOM-1 session), and V1b (L33) now rejects it. Pentagon coordinates need no domain for soundness (L31).

**L31. Pentagon coordinates need no domain** (proved by DES-F2 by reading; re-derived by the AUD-2 audit (Fable 5.1, 2026-10-23, review document `admin-twelve-fable-audit2-review-2026-10-23.md`, item R3, available on request), which re-read `kisscrit/core/model.py` in full: the proof uses no fact about the position of a pentagon coordinate; AUD-2 ran no versioned script for it; not re-derived by the DOM-1 code session). Numbers L18–L30 are reserved for other DES-F2 documents (C3, C4) and may not exist yet; L31–L33 are used as proposed in `admin-twelve-fable-design-domain-2026-10-21.md` §2.3. For any box T (finite endpoints) in the parameter space of a chart, inside the θ-domain of L17 and with rhombus coordinates in (0, π], every corner enclosure of `kisscrit/core/model.py:Graph.corners` (with or without tightening) contains the true corners of every valid configuration of the chart's region whose parameters lie in T, and every validity flag is True at such a configuration. *Proof (DES-F2).* `pentagon` is a composition of interval operations, each of which encloses the range of its real counterpart on any input box, with `clip` and `asin_angle` returning the image of the clamp; at a valid point of the region the clamped formula equals the true corner (L8, L10′: the clamp is harmless), so the enclosure over any T containing the point contains it. The flags `ok` are each of the form "the argument's interval certainly misses the admissible range", which is False at a valid point, and are monotone under box inclusion. `_tighten` narrows a value only by evaluating the same clamped formula on a face of the box chosen by the sign of a derivative enclosure of the Clarke derivative on the box (L16, monotonicity tightening); its soundness does not depend on where the box lies. The corner-range check "every corner in [α(θ), π]" uses `alpha_range` (θ-domain) and L3. ∎ (A box far outside [α, π] in a pentagon coordinate is only wider; where cos d3 → 1 the enclosure is (−∞, ∞) and nothing is proved, as in L16.)

**L32. Rhombus coordinate domain (0, π]** (proved by DES-F2 by reading; re-derived by the AUD-2 audit (Fable 5.1, 2026-10-23, review document `admin-twelve-fable-audit2-review-2026-10-23.md`, item R3, available on request): agrees, no gap; also re-derived by the DOM-1 code session, 2026-10-22, in `session-log.md`: agrees, no gap). Let a rhombus coordinate's box be [A_lo, A_hi] with 0 < A_lo and A_hi ≤ π (reals; the checker decides this on the certificate's doubles, in Arb for A_hi). Then ex = [2β(θ_lo, A_hi), 2β(θ_hi, A_lo)] (each end rounded outward; `model.py`, the range `ex` of the corner O) encloses 2β(θ, A) for all (θ, A) in the box, with θ in the L17 domain. *Proof.* For A ∈ (0, π], A/2 ∈ (0, π/2], so sin(A/2) > 0, cos(A/2) ≥ 0, and cot(A/2) = cos(A/2)/sin(A/2) ≥ 0 is decreasing in A. By D1, cos θ > 0 and 1/cos θ is positive and increasing in θ on [1.0, 1.5] ⊂ (0, π/2). So cot(A/2)/cos θ is decreasing in A and nondecreasing in θ, and arctan is increasing: β(θ, A) = arctan(cos(A/2)/(cos θ sin(A/2))) (L7) has its minimum over the box at (θ_lo, A_hi) and its maximum at (θ_hi, A_lo). ∎ Outside the domain the claim fails: for A_hi ∈ (2π, 3π), β(θ_lo, A_hi) = β(θ_lo, A_hi − 2π) by periodicity, which is not a lower bound of 2β over A ∈ [A_lo, π]; the false range excludes the true value of O, and a vertex sum built from it can falsely exclude 2π. For A ∈ (π, 2π) the monotonicity in θ is reversed.

**L33. The check V1b** (proved by DES-F2 by reading; re-derived by the AUD-2 audit (Fable 5.1, 2026-10-23, review document `admin-twelve-fable-audit2-review-2026-10-23.md`, item R3, available on request): the code matches the statement; AUD-2 also replayed the DES-F2 forgery (accepted at cc375bc and 5da4bba, rejected at 111ca4b by the reference and the kernel front end) and re-scanned the 45,260 committed `cond` start boxes (0 violations), both in unversioned scratch scripts; the code is `kisscrit/core/checker.py:structural_check`). For every chart that is not `not_run` and every coordinate j: lo0[j] and hi0[j] are finite; for a rhombus coordinate 0 < lo0[j] and hi0[j] ≤ π; for a pentagon coordinate 0 ≤ lo0[j] and hi0[j] ≤ 2π (hygiene: L31 says no bound is needed for soundness, but a bounded box keeps the width check and the developing map meaningful). "≤ π" is decided as arb(hi0[j]) <= arb.pi() with π enclosed at 128 bits, i.e. "certainly ≤": the double below π (3.141592653589793) passes, the next double up does not; lo0 is compared as the double it is. The θ coordinate is covered by V3 (L17) and the straight-corner coordinate keeps V1's rule (its box contains π). Every node box of the tree lies inside the start box (certificate-format §3, exact replay), so the condition on the start box is the condition on every box. Together with V1 (lo0[j] ≤ a_lo[j], hi0[j] ≥ a_hi[j]) a rhombus start box lies in (0, α(θ_lo)↓] × [2α(θ_hi)↑, π] at its ends, where L32 applies. A message of the form `V1b: chart … coordinate j box [lo, hi] outside the model domain (docs/lemmas.md L32)` is a structural failure. The C kernel front end (`kisscrit/ckernel.py`) calls this `structural_check`; `csrc/` has no copy of the rule. Checked on the committed cond certificates: see the DOM-1 entry in `docs/session-log.md` (computed).

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
| U | `lemmas.py U` | proved (1805 boxes, about 1 s) |

## Still to do

Nothing for the lemmas themselves. For the appendix (D3), write L6, L8 and the uniqueness argument of L9 in LaTeX, and add the references Fejes Tóth 1943 and Schütte–van der Waerden 1953 (volume and page numbers to be checked against the originals before they enter a bibliography; not verified here).
