# Hull triangulations of twelve points with minimum angle above 60°

*2026-10-10. Numerical; nothing here is certified.*

## Question

Does every configuration of 12 points on S² with all pairwise angles ≥ 60° have an icosahedral convex hull? If so, every contact graph above r = 1 would be a subgraph of one fixed graph.

## Answer: no

The critical configurations above r = 1 are themselves counterexamples. Their hulls, computed from `data/n12_above1.json`, are:

| r | contacts | co-index | hull type | degrees | non-contact hull edges (max angle) |
|---|---|---|---|---|---|
| 1.10850854 (DOD) | 30 | 0 | T0 | 5¹² | none |
| 1.0010121 (M) | 23 | 0 | T5 | 4²5⁸6² | 7 (92.58°) |
| 1.00100739 | 21 | 1 | T1 | 4²5⁸6² | 9 (92.40°) |
| 1.0008452 | 21 | 1 | T5 | 4²5⁸6² | 9 (91.49°) |
| 1.00058342 | 20 | 2 | T5 | 4²5⁸6² | 10 (92.60°) |
| 1.00030629 | 21 | 1 | T7 | 4⁴5⁴6⁴ | 9 (92.57°) |
| 1.00030625 | 20 | 2 | T1 | 4²5⁸6² | 10 (92.57°) |
| 1.00030625 | 19 | 3 | T5 | 4²5⁸6² | 11 (92.57°) |
| 1.00023741 | 21 | 1 | T7 | 4⁴5⁴6⁴ | 9 (91.37°) |
| 1.00021856 | 20 | 2 | T7 | 4⁴5⁴6⁴ | 10 (91.60°) |
| 1.00021159 | 20 | 2 | T7 | 4⁴5⁴6⁴ | 10 (91.82°) |
| 1.00021116 | 19 | 3 | T7 | 4⁴5⁴6⁴ | 11 (91.81°) |

T_k is the k-th triangulation in `data/triangulations_n12.json`, which is in BFS order of flips from the icosahedron. T1 and T5 are the icosahedron with one edge flipped. T7 has two edges flipped. Every contact graph is a subgraph of the icosahedral graph (checked by VF2), so the flipped diagonal is never a contact.

## θ*(T)

θ*(T) is the supremum of the minimum angle over 12 points whose hull is combinatorially T, with faces allowed to become coplanar in the limit. `kisscrit.hull.theta_star` computes local maxima with SLSQP. It uses random Tutte-embedding starts, and also the critical configurations as starts. Every value it returns is a **lower bound** on θ*(T).

- T0: 63.43495° (the icosahedron).
- T1: ≥ 60.0633°. Random starts gave only 59.57°; seeding from a critical configuration was necessary.
- T5: ≥ 60.0398°. With every convexity slack forced ≥ 10⁻³, it is still 60.0392°, so the realization is strictly convex.
- T7: ≥ 60.0101°.
- Survey complete (all 7595 types). Of the other 7591 types:
  - 8 reach 60° but no more, to optimizer accuracy;
  - 2871 stay below 60°, the largest at 57.52°;
  - for 4712 the local method found no convex inscribed realization at all (a failure of the method, not a proof that none exists).

## Two more pieces of evidence for (a)

- **Flip walk** (`scripts/hull12/flipwalk.py`). Start from the optimal realizations of T0, T1, T5 and T7. Flip each hull edge in turn and re-optimize θ* for the flipped triangulation, seeded at the parent. Keep anything above 60° and repeat. The walk closes on {T0, T1, T5, T7}: no flip leads to a fifth type above 60°.
- **Sampling** (`scripts/hull12/sample_hulls.py`). Of 150 random starts ascended until the minimum angle passed 60.0005° (overshooting to r ≤ 1.002), 144 had icosahedral hulls, 3 had T1, and 3 did not reach the level. Generic configurations above 60° are icosahedral; the other three types occupy thin regions near the critical configurations.

## The corrected statement

**Conjecture.** If 12 points on S² have all pairwise angles > 60°, then:

(a) their hull triangulation is one of T0, T1, T5, T7;

(b) their contact graph is a subgraph of the icosahedral graph.

Given (a), part (b) is a finite check on three triangulations. Proving (a) needs certified upper bounds θ*(T) ≤ 60° for 7591 triangulations. A degree or area argument (Fejes Tóth) should dispose of most of them. The rest can be done by interval branch-and-bound on a small polynomial program.

## Candidate graphs (if the conjecture holds)

These are subgraphs of the icosahedral graph, counted up to its 120 symmetries by a Burnside count (`scripts/hull12/burn.c`):

- 30392 with minimum degree ≥ 3 (18 to 30 edges);
- 1229174 if degree 2 is allowed.

## Rigorous lemma

Hull faces have empty circumcaps of radius < π/3, by the kissing number 12. A grid computation puts the face angles in [35.8°, 140.6°], which excludes degree 11.
