# Exact descriptions of the N = 12 critical configurations found so far (EXACT-1)

Status: data and scripts written by the EXACT-1 code session, 2026-10-21. This page describes what has been
found so far. It is **not a completeness claim**: the list of critical configurations at r ≥ 1 is complete only
if the full-route branch and bound (all 7,595 triangulations) shows it, and that run is separate work.

## What "exactly verified" means here

An entry is *exactly verified* when an exact computation in the number field K = Q(c) accepted every listed
identity and inequality. Exact means: rational arithmetic and arithmetic modulo the minimal polynomial P of
c; zero tests are exact; the sign of a nonzero element of K is decided on a rational interval that contains
exactly one root of P (an exact Sturm count), by an Arb ball that excludes zero. Floating point never decides
anything. This is a statement about these files only; it is not the checker of the completeness proof
(`docs/certificate-format.md`, `scripts/completeness/check_certificates.py`).

## Conventions

- Contact angle θ, c = cos θ, s = sin(θ/2), touching radius r = s/(1 − s). Polynomials are stored with the
  constant term first, with integer coefficients, primitive, positive leading coefficient.
- Gram matrix G_ij = ⟨u_i, u_j⟩, entries in K, stored as coefficient lists in 1, t, …, t^(d−1) with t = c.
- Weight (stress) space K(C): weights w on the contacts with Σ_j w_ij (u_j − c u_i) = 0 at every point. In
  Gram form: Σ_j w_ij (G_bj − c G_bi) = 0 for every point i and each of three frame points b. Its
  nonnegative part is a polyhedral cone; its extreme rays are listed, normalized to sum 1.
- Rigidity rank: rank of the contact rows on the tangent space of (S²)¹². Flex dimension: motions keeping each
  contact length, modulo rotations. Kernel: motions keeping all contacts equal to a common, varying c, modulo
  rotations (the space of the bordered Hessian in `kisscrit/certify.py`).
- q(v) = 2 Σ_e w_e ⟨v_i, v_j⟩ − Σ_i λ_i |v_i|², λ_i = c Σ_j w_ij, Σ w = 1: the second derivative of c along
  the variety of equal contacts. **Co-index** = number of negative directions of q on the kernel (r increases);
  **descending index** = number of positive directions (r decreases). Nondegenerate means q has only the three
  rotation directions as null space. Where the stress is not unique and the kernel is not trivial (FCC, HCP),
  q depends on the stress and no single index is given.

## Files

| file | what |
|---|---|
| `scripts/exact/kfield.py` | K = Q[t]/(P): irreducibility, Sturm isolation, exact signs, rank, null space, inertia |
| `scripts/exact/hiprec.py` | guesses only: Gauss-Newton on the KKT system at high precision, LLL relations |
| `scripts/exact/verify.py` | exact checks: Gram, stress space, motions and index, symmetry, canonical graph, s and r |
| `scripts/exact/build_e2.py` | DOD, M, B, r1, A, the hexagonal antiprism band, FCC, HCP |
| `scripts/exact/build_e3.py` | the values without an exact polynomial before; quotient guesses for Gram entries |
| `scripts/exact/elim_system.py` | the symmetric-subspace KKT system for msolve (an elimination attempt; see below) |
| `scripts/exact/weights_e4.py` | weight spaces from the stored exact data: cones, extreme rays, orbits |
| `scripts/exact/table_md.py` | writes `data/exact/critical-configs-n12.md` |
| `data/exact/critical-configs-n12.json` | the exact data (strings), with per-entry script and commit |
| `data/exact/critical-configs-n12.md` | the human-readable table |
| `tests/test_exact_configs.py` | re-verifies every exact claim from the JSON strings; mutants must fail |

Reproduce: `python3 scripts/exact/build_e2.py` (about 45 s), `python3 scripts/exact/build_e3.py LABEL DPS`
for 1.0003062544 1500, 1.0010073939 3000, 1.0008452029 4000 (20 s, 12 min, 41 min), then
`python3 scripts/exact/weights_e4.py` (50 s) and `python3 scripts/exact/table_md.py`. Requires python-flint
(0.9.0 used), sympy, mpmath, networkx.

## How the data were obtained

1. A high-precision approximation of each configuration: Gauss-Newton on its KKT system from
   `data/n12_above1.json` (DOD, M and the saddles), or a closed form (antiprism band, FCC, HCP).
2. Guesses by LLL: the minimal polynomial of c, then each Gram entry in K. For P far from monic (degrees 46
   and 59, leading coefficients near 10^10 and 10^16) the entries are quotients −a(c)/b(c) of elements with
   small coefficients in the integral basis ω_k = a_d c^k + … + a_(d−k+1) c; their power-basis coefficients
   run to hundreds of digits, so they are guessed in quotient form and converted exactly.
3. Exact checks in K (`verify.py`), as listed in the tests. A wrong guess makes a check fail; the guessing
   step plays no role in correctness.

## Results (found so far)

Eleven configurations are exactly verified: DOD, M, 1.0010073939, 1.0008452029, B, r1, 1.0003062544, A,
the antiprism band, FCC and HCP. Their c polynomials equal those of the 2026-10-01 and 2026-10-08 notes
(up to sign); the high-precision relations of the 2026-10-08 note (degrees 28, 46, 59) are now exact.
All eleven have a strictly positive stress. The eleven contact graphs are pairwise non-isomorphic. See
`data/exact/critical-configs-n12.md` for the table.

Findings worth noting:

- **M** has a two-dimensional stress cone with two extreme rays, one vanishing on (0,1) and one on (5,9); the
  stored σ of `data/n12_above1.json` is the second ray. Its C₂-symmetric stress is strictly positive.
- **Antiprism band**: the uniform weight is not a stress; the D_6d-symmetric stress is positive with
  ring/cross ratio 2 + √3; 12 extreme rays.
- **DOD**: a nine-dimensional stress space whose nonnegative cone has 278 extreme rays, in I_h-orbits of
  sizes 60, 60, 60, 30, 30, 20, 12, 6 (each ray checked extreme; the set is closed under I_h; that the list is
  complete rests on the double description algorithm).
- **FCC, HCP**: four-dimensional stress spaces with 4 and 5 extreme rays; the uniform weight is a stress at FCC
  and not at HCP; the HCP witness of control item C5a (1/20, 11/150, 7/200) is a stress. This agrees with
  `scripts/r1/exact_lambda_stress_space.py`, by an independent code path.
- `data/n12_above1.json` gives **1.0003062547** a symmetry group of order 2, but its stored contact graph has
  only the trivial automorphism, so the stored symmetry and the stored contacts cannot both be right.

## Not found

- **1.0003062547, 1.0002374095, 1.0002115891, 1.0002111636**: no exact polynomial. For 1.0003062547 and
  1.0002374095, LLL at 8000 digits finds no integer relation of degree ≤ 150 with coefficients below about
  10^42 (high precision, not proof; details in `docs/session-log.md`). For 1.0002115891 and 1.0002111636 the
  known status is that of `docs/open-questions.md` section 3 (no relation of degree ≤ 90 at 1600 digits); the
  same degree-150 probe was still running when this page was written, and its result goes to the session log.
  The next tier (degree 200 at 12,000 digits) is estimated at 2.7–3.4 hours
  per LLL. Elimination with msolve 0.11.0 on the symmetric-subspace KKT system exceeded 7 GB of memory even on
  saddle B (degree 13). A dedicated later search is reserved (Wöden, 2026-10-21).
- **r = 1**: only FCC and HCP have data. The contact graphs with 15–23 contacts mentioned in the 2026-10-01
  note and the KKLS M₅ stratum have no stored data: not found, no data.

## Proposal (not run)

Once the full-route run has terminal boxes, compare each terminal box near r_min, and the `cond` terminal set,
with the canonical contact graphs and the exact c intervals stored here.
