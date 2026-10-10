# Square systems for other max–min problems

*Notes, 2026-10-10. A plan and a framework. Only the sphere case is implemented in `kisscrit.kkt`.*

## The general shape

Let V ⊂ ℝ^n be a smooth real algebraic variety, cut out by h_1 = … = h_m = 0 with independent gradients. Let g_1, …, g_M be polynomials, here the pair functions. Let G be a group acting on V and preserving each g_e up to permutation. The function

  F(x) = min_e g_e(x)   (or max, after a sign change)

is semialgebraic and locally Lipschitz. Its Clarke critical points with active set A are the solutions of

  g_e(x) = t (e ∈ A),  h_k(x) = 0,  Σ_{e∈A} σ_e ∇g_e(x) = Σ_k μ_k ∇h_k(x),  Σ σ_e = 1,

together with σ ≥ 0 and g_f(x) > t for f ∉ A.

The count is:

- unknowns: n (x), plus 1 (t), plus |A| (σ), plus m (μ);
- equations: |A| + m + n + 1.

The system is square for every active set A. A continuous symmetry group G of dimension g contributes g Noether identities: the balance equation paired with the generators of G vanishes identically. It also contributes a g-dimensional orbit of solutions. One gauge condition per generator restores a square system with isolated solutions in the generic case.

Two consequences follow, the same as for kissing spheres.

1. **Algebraicity.** Isolated critical values are real algebraic numbers. A degree bound is the Bézout number (or mixed volume) of the gauge-fixed system, maximized over A.
2. **Finiteness.** There are finitely many critical values, by Tarski–Seidenberg applied to the union over A; this is BDLS for this case. Only finitely many A occur, and for each one the critical-value set is a semialgebraic subset of ℝ with empty interior.

Inequality constraints, such as walls, are added the same way. Each active wall is one more equation and one more multiplier, with the sign condition λ ≥ 0. The system stays square.

## Instances

| problem | V | g_e | walls | symmetry (gauge) |
|---|---|---|---|---|
| kissing / Tammes on S^{d−1} (implemented) | (S^{d−1})^N | ⟨x_i, x_j⟩ (max) | none | O(d): d(d−1)/2 |
| spherical codes in RP^{d−1} (lines) | (S^{d−1})^N | ⟨x_i, x_j⟩² (max) | none | O(d) |
| circles in a square / points in a box [0,1]^d | ℝ^{Nd} | \|x_i − x_j\|² (min) | x_{ia} ≥ 0, x_{ia} ≤ 1 | finite (hyperoctahedral) only, so no gauge |
| equal balls in a ball | ℝ^{Nd} | \|x_i − x_j\|² (min) | \|x_i\|² ≤ (1 − ρ)² | O(d) |
| points on a flat torus ℝ^d/Λ | ℝ^{Nd} | \|x_i − x_j − λ\|² for the relevant λ ∈ Λ | none | translations: d |
| Grassmannian packings | products of Grassmannians as projector varieties | chordal distance tr(P_i P_j) | none | O(d) |
| packings on other surfaces (spheroid, torus of revolution) | products of the surface | ambient \|x_i − x_j\|² | none | the surface's isometries |

Points in a box show where the square count has to be done carefully. The objective min \|x_i − x_j\|² has two kinds of contact: pair contacts and wall contacts. Pair contacts carry the common value t. Wall contacts carry their own equation x_{ia} = 0 or 1 and a multiplier with no t. With that, the count is again |A_pairs| + |A_walls| + Nd + 1 equations in Nd + 1 + |A_pairs| + |A_walls| unknowns.

The torus is the case where the gauge is translations, and the stress balance is the classical periodic-framework equilibrium.

## Rattlers and degenerate strata

A point with no stressed contact is a rattler. It does not appear in the balance equations, so the solution set has its region as a factor. That makes it a Morse–Bott family (or worse), not an isolated point. The enumeration therefore runs over stressed graphs on n ≤ N vertices. Each solution is then given the rattler region for the remaining N − n points, together with the number of components of that region. `kisscrit.search.euler_check` does this for d = 3. For N = 6 it is what makes the signed count equal χ(M_{0,6}) = −6.

## What is needed to implement the other rows

`kkt.residual_and_jacobian` is written for the sphere. A general version needs three things: a callable list of (g_e, ∇g_e, ∇²g_e), a list of constraints h_k, and a list of gauge generators. The stress LP, the co-index computation (Lagrangian Hessian on the kernel of the active Jacobian, modulo the gauge), monotone ascent and graph enumeration use only those three. Packings in a square and in a ball are the best first tests, because best-known values are tabulated for them (E. Specht, packomania.com).
