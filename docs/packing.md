# Equal balls in containers (2026-10-10)

Purpose: test the min-type toolkit (`kisscrit.packing`) on disks and balls in containers. Everything here is numerical.

Containers: the unit square and unit cube (walls), the square and hexagonal flat 2-tori and the cubic 3-torus (modulo
translations), and the unit disk and ball (modulo rotations). The packing radius is
rho = min(pair half-distances, wall distances). It is a min-type function, and each active set gives a square system;
see `square-systems.md`.

Reproduce:

    python scripts/pack_maxima.py sq 2 10 60        # best radius, minimal polynomial, distinct maxima
    python scripts/pack_explore.py sq 4 60 6000     # census of critical points (ascent, sampled active sets, surgery)
    python scripts/pack_euler.py                    # Euler checks with Morse-Bott families and rattler strata
    python scripts/pack_tables.py                   # LaTeX tables used in the paper (Appendix B)

## Tables (LaTeX rows: maxima, then Euler checks)

```
square & 2 & $0.2928932188$ & $2\rho^{2}-4\rho+1$ & 5 & 1\\
 & 3 & $0.2543330950$ & $\rho^{4}-16\rho^{3}+20\rho^{2}-8\rho+1$ & 7 & 1\\
 & 4 & $0.2500000000$ & $4\rho-1$ & 12 & 1\\
 & 5 & $0.2071067812$ & $4\rho^{2}+4\rho-1$ & 12 & 2\\
 & 6 & $0.1876806011$ & $92\rho^{2}+52\rho-13$ & 13 & 7\\
 & 7 & $0.1744576302$ & $13\rho^{2}-8\rho+1$ & 14 & 7 (1 rattler)\\
 & 8 & $0.1705406887$ & $32\rho^{4}-32\rho^{3}-8\rho^{2}+8\rho-1$ & 20 & 8\\
 & 9 & $0.1666666667$ & $6\rho-1$ & 24 & 6\\
 & 10 & $0.1482043226$ & --- & 21 & 13\\
\hline
disk & 2 & $0.5000000000$ & $2\rho-1$ & 3 & 1\\
 & 3 & $0.4641016151$ & $\rho^{2}+6\rho-3$ & 6 & 1\\
 & 4 & $0.4142135624$ & $\rho^{2}+2\rho-1$ & 8 & 1\\
 & 5 & $0.3701919082$ & $\rho^{4}+20\rho^{3}+10\rho^{2}-20\rho+5$ & 10 & 1\\
 & 6 & $0.3333333333$ & $3\rho-1$ & 10 & 2 (degenerate)\\
 & 7 & $0.3333333333$ & $3\rho-1$ & 18 & 1\\
 & 8 & $0.3025933883$ & $\rho^{6}+42\rho^{5}+119\rho^{4}-84\rho^{3}-49\rho^{2}+42\rho-7$ & 14 & 1 (1 rattler)\\
 & 9 & $0.2767686539$ & $\rho^{4}+12\rho^{3}-2\rho^{2}-4\rho+1$ & 16 & 5 (1 rattler)\\
 & 10 & $0.2622589242$ & --- & 20 & 8\\
\hline
square torus & 2 & $0.3535533906$ & $8\rho^{2}-1$ & 4 & 1\\
 & 3 & $0.2588190451$ & $16\rho^{4}-16\rho^{2}+1$ & 5 & 1\\
 & 4 & $0.2588190451$ & $16\rho^{4}-16\rho^{2}+1$ & 10 & 1\\
 & 5 & $0.2236067977$ & $20\rho^{2}-1$ & 10 & 1\\
 & 6 & $0.2002027743$ & $36\rho^{4}-12\rho^{3}-14\rho^{2}-2\rho+1$ & 12 & 3\\
 & 7 & $0.1830127019$ & $8\rho^{2}+4\rho-1$ & 14 & 3\\
 & 8 & $0.1830127019$ & $8\rho^{2}+4\rho-1$ & 20 & 3\\
\hline
hexagonal torus & 2 & $0.2886751346$ & $12\rho^{2}-1$ & 3 & 1\\
 & 3 & $0.2886751346$ & $12\rho^{2}-1$ & 9 & 1\\
 & 4 & $0.2500000000$ & $4\rho-1$ & 12 & 1\\
 & 5 & $0.1980717478$ & $64\rho^{4}-28\rho^{2}+1$ & 9 & 1\\
 & 6 & $0.1889822365$ & $28\rho^{2}-1$ & 15 & 1\\
 & 7 & $0.1889822365$ & $28\rho^{2}-1$ & 21 & 1\\
 & 8 & $0.1666666667$ & $6\rho-1$ & 21 & 2\\
\hline
cube & 2 & $0.3169872981$ & $8\rho^{2}-12\rho+3$ & 7 & 1\\
 & 3 & $0.2928932188$ & $2\rho^{2}-4\rho+1$ & 12 & 1\\
 & 4 & $0.2928932188$ & $2\rho^{2}-4\rho+1$ & 18 & 3\\
 & 5 & $0.2639320225$ & $4\rho^{2}-20\rho+5$ & 18 & 6\\
 & 6 & $0.2573593129$ & $4\rho^{2}-36\rho+9$ & 24 & 10\\
 & 7 & $0.2501251005$ & --- & 23 & 10\\
 & 8 & $0.2500000000$ & $4\rho-1$ & 36 & 7\\
\hline
ball & 2 & $0.5000000000$ & $2\rho-1$ & 3 & 1\\
 & 3 & $0.4641016151$ & $\rho^{2}+6\rho-3$ & 6 & 1\\
 & 4 & $0.4494897428$ & $\rho^{2}+4\rho-2$ & 10 & 1\\
 & 5 & $0.4142135624$ & $\rho^{2}+2\rho-1$ & 11 & 1 (degenerate)\\
 & 6 & $0.4142135624$ & $\rho^{2}+2\rho-1$ & 18 & 1\\
 & 7 & $0.3859135574$ & $\rho^{6}+6\rho^{5}-3\rho^{4}-8\rho^{3}+12\rho^{2}-6\rho+1$ & 19 & 1\\
 & 8 & $0.3780248448$ & $\rho^{4}+8\rho^{3}+4\rho^{2}-8\rho+2$ & 24 & 1\\
\hline
cubic $3$-torus & 2 & $0.4330127019$ & $16\rho^{2}-3$ & 8 & 1\\
 & 3 & $0.3535533906$ & --- & 12 & 1\\
 & 4 & $0.3535533906$ & $8\rho^{2}-1$ & 24 & 1\\
 & 5 & $0.2844826943$ & --- & 13 & 11\\
\hline

square & 2 & 2 & 0 & 0 & yes\\
square & 3 & 3 & 0 & 0 & yes\\
square & 4 & 8 & 0 & 0 & yes\\
square & 5 & 7 & 1320 & 0 & no\\
disk & 2 & 1 & 1 & 1 & yes\\
disk & 3 & 2 & -1 & -1 & yes\\
disk & 4 & 3 & 18 & 2 & no\\
square torus & 2 & 2 & -1 & -1 & yes\\
square torus & 3 & 6 & 2 & 2 & yes\\
square torus & 4 & 9 & 37 & -6 & no\\
hexagonal torus & 2 & 2 & -1 & -1 & yes\\
hexagonal torus & 3 & 3 & 2 & 2 & yes\\
hexagonal torus & 4 & 8 & -6 & -6 & yes\\
cube & 2 & 3 & 2 & 2 & yes\\
cube & 3 & 6 & -18 & 6 & no\\
ball & 2 & 1 & 0 & --- & ---\\
ball & 3 & 2 & 1 & --- & ---\\
cubic $3$-torus & 2 & 3 & 1 & 1 & yes\\
cubic $3$-torus & 3 & 13 & 2 & 2 & yes\\
```

## Euler targets

- Square: chi(Conf_n(R^2)) = 0 for n >= 2. Cube: chi(Conf_n(R^3)) = n!.
- Flat 2-torus modulo translations: chi(Conf_{n-1}(T^2 - pt)) = (-1)^{n-1} (n-1)!.
- Cubic 3-torus modulo translations: (n-1)!.
- Disk modulo rotations: M_{0,n+1}, with chi = prod_{k=2}^{n-1} (1-k).
- Ball: no target, because SO(3) does not act freely.

## Known limitations

- Families on tori that wrap around are circles. `family_count` assumes contractible families, and its lattice representative is not always minimal (for example, it gave a non-integer count for the square torus with n = 4).
- Family counting is not implemented for the ball and disk (rotation gauge).
- Rattler strata are supported only over bases with families of dimension at most 1, and only with one rattler. Their Euler characteristic is computed on a grid of spacing 0.01 (2-D) or 0.04 (3-D).
