#!/usr/bin/env python3
"""
hunt_monomial_mismatch.py

Instead of just comparing convex hulls, compare the FULL monomial support
of the Burau and Lawrence charpolys. The Newton polytope (convex hull)
could agree even if the monomial supports differ at interior points.

Also: systematically look for braids where the monomial support is
"unusual" — few interior monomials relative to the hull.
"""

from sympy import symbols, eye, zeros, Matrix, cancel, Poly, Rational, fraction
import sys

q, x = symbols('q x')
Q_val = Rational(2)
a_val = Rational(3)

def build_Ni(i, gi, n, d):
    nb = n - 1; Id = eye(d); Zd = zeros(d)
    blocks = [[Id.copy() if r == c else Zd.copy() for c in range(nb)] for r in range(nb)]
    if i == 0:
        blocks[0][0] = -gi; blocks[0][1] = Zd; blocks[1][0] = Id; blocks[1][1] = Id
    elif i == n - 2:
        blocks[nb-2][nb-2] = Id; blocks[nb-2][nb-1] = gi; blocks[nb-1][nb-2] = Zd; blocks[nb-1][nb-1] = -gi
    else:
        p = i - 1
        blocks[p][p]=Id; blocks[p][p+1]=gi; blocks[p][p+2]=Zd
        blocks[p+1][p]=Zd; blocks[p+1][p+1]=-gi; blocks[p+1][p+2]=Zd
        blocks[p+2][p]=Zd; blocks[p+2][p+1]=Id; blocks[p+2][p+2]=Id
    return blocks
def left_multiply_blocks(blocks, M, nb):
    return [[M * blocks[r][c] for c in range(nb)] for r in range(nb)]
def multiply_block_matrices(A, B, nb, d):
    C = [[zeros(d) for _ in range(nb)] for _ in range(nb)]
    for r in range(nb):
        for c in range(nb):
            acc = zeros(d)
            for k in range(nb): acc = acc + A[r][k] * B[k][c]
            C[r][c] = acc
    return C
def assemble(blocks, nb, d):
    N = nb * d; M = zeros(N)
    for r in range(nb):
        for c in range(nb): M[r*d:(r+1)*d, c*d:(c+1)*d] = blocks[r][c]
    return M
def reduced_lm_colored(rho_g, rho_sigma, partition=None):
    n = len(rho_g); num_sigmas = len(rho_sigma); d = rho_g[0].shape[0]
    if partition is None: partition = [1] * n
    is_boundary = [partition[i] != partition[i+1] for i in range(num_sigmas)]
    nb = n - 1; rho_plus = []
    for i in range(num_sigmas):
        if not is_boundary[i]:
            blocks = build_Ni(i, rho_g[i], n, d)
            blocks = left_multiply_blocks(blocks, rho_sigma[i], nb)
            rho_plus.append(cancel(assemble(blocks, nb, d)))
        else:
            gi_conj = cancel(rho_g[i] * rho_g[i+1] * rho_g[i].inv())
            blocks_N = build_Ni(i, rho_g[i], n, d); blocks_Np = build_Ni(i, gi_conj, n, d)
            blocks_prod = multiply_block_matrices(blocks_Np, blocks_N, nb, d)
            si_sq = rho_sigma[i] * rho_sigma[i]
            blocks_prod = left_multiply_blocks(blocks_prod, si_sq, nb)
            rho_plus.append(cancel(assemble(blocks_prod, nb, d)))
    return rho_plus
def free_group_gens(rho_sigma, partition=None):
    n = len(rho_sigma)
    if partition is None: partition = [1] * (n + 1)
    last_is_boundary = (partition[n-1] != partition[n])
    gn = rho_sigma[n-1] if last_is_boundary else rho_sigma[n-1] * rho_sigma[n-1]
    rho_g = [None] * n; rho_g[n-1] = cancel(gn)
    for i in range(n-2, -1, -1):
        rho_g[i] = cancel(rho_sigma[i].inv() * rho_g[i+1] * rho_sigma[i])
    return rho_g
def build_burau_B4():
    return reduced_lm_colored([Matrix([[q]])]*4, [Matrix([[1]])]*3)
def build_lawrence_B4():
    L4Pre = reduced_lm_colored([Matrix([[Q_val]])]*4+[Matrix([[a_val]])], [Matrix([[1]])]*4, [1,1,1,1,2])
    L4PreFree = free_group_gens(L4Pre, [1,1,1,1,2])
    return reduced_lm_colored([g*q for g in L4PreFree], L4Pre[0:3])
def newton_polytope_xq(M):
    cp_expr = cancel(M.charpoly(x).as_expr()); num, den = fraction(cp_expr)
    p = Poly(num, x, q); monoms_num = p.monoms()
    q_shift = Poly(den, q).degree() if den.has(q) else 0
    return set((xd, qd - q_shift) for (xd, qd) in monoms_num)
def convex_hull_2d(points):
    points = sorted(set(points))
    if len(points) <= 2: return points
    def cross(O, A, B): return (A[0]-O[0])*(B[1]-O[1]) - (A[1]-O[1])*(B[0]-O[0])
    lower = []
    for p in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0: lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0: upper.pop()
        upper.append(p)
    return sorted(set(lower[:-1] + upper[:-1]))
def scale_polytope(vertices, s):
    return sorted(set((Rational(v[0])*s, Rational(v[1])*s) for v in vertices))
def eval_word(word, gens, inv_gens):
    M = eye(gens[0].shape[0])
    for w in word:
        if w > 0: M = M * gens[w-1]
        else: M = M * inv_gens[-w-1]
    return cancel(M)
def mat_power(M, n):
    if n == 0: return eye(M.shape[0])
    if n < 0: M = cancel(M.inv()); n = -n
    result = eye(M.shape[0]); base = M
    while n > 0:
        if n % 2 == 1: result = cancel(result * base)
        base = cancel(base * base); n //= 2
    return result

if __name__ == '__main__':
    print('Building representations...'); sys.stdout.flush()
    bg = build_burau_B4(); bgi = [cancel(g.inv()) for g in bg]
    lg = build_lawrence_B4(); lgi = [cancel(g.inv()) for g in lg]
    print('Done.\n'); sys.stdout.flush()

    scale = Rational(1, 4)
    counts = [0, 0]  # hull match, hull mismatch
    monomial_mismatches = 0

    braids = [
        # Sparse charpolys
        ('s1^2.s2.s3', [1, 1, 2, 3]),
        ('s1^3.s2.s1^-1.s3', [1, 1, 1, 2, -1, 3]),
        # Powers of sparse
        ('(s1^2.s2.s3)^2', [1, 1, 2, 3]*2),
        ('(s1^2.s2.s3)^3', [1, 1, 2, 3]*3),
        # Products of different sparse braids
        ('s1^2.s2.s3.s3^2.s2.s1', [1, 1, 2, 3, 3, 3, 2, 1]),
        # Zero writhe with lots of cancellation potential
        ('s1.s2.s3.s3^-1.s2^-1.s1^-1', [1,2,3,-3,-2,-1]),  # = identity
        ('s1.s2.s3.s2^-1.s1^-1.s3^-1', [1,2,3,-2,-1,-3]),
        # Braids where Alexander poly might factor
        ('(s1.s2)^3', [1,2]*3),  # half twist on 3 strands
        ('(s1.s2)^3.s3', [1,2]*3+[3]),
        ('(s1.s2)^3.s3^-1', [1,2]*3+[-3]),
        ('(s1.s2)^6', [1,2]*6),  # full twist on 3 strands
        ('(s1.s2)^6.s3', [1,2]*6+[3]),
        # Garside element and variants
        ('s1.s2.s1.s3.s2.s1', [1,2,1,3,2,1]),  # Delta
        ('s1.s2.s1.s3.s2.s1.s3', [1,2,1,3,2,1,3]),
        ('s1.s2.s1.s3.s2.s1.s3^-1', [1,2,1,3,2,1,-3]),
        # Birman-Menasco braids (known pseudo-Anosov)
        ('s1.s2^-1.s3', [1,-2,3]),  # minimal pseudo-Anosov for B4
        ('(s1.s2^-1.s3)^2', [1,-2,3]*2),
        ('(s1.s2^-1.s3)^3', [1,-2,3]*3),
        ('(s1.s2^-1.s3)^5', [1,-2,3]*5),
        ('(s1.s2^-1.s3)^10', [1,-2,3]*10),
        # Dehornoy braids (handle reduction-resistant)
        ('s1.s2.s1.s3^-1', [1,2,1,-3]),
        ('s1^2.s2.s3^-1.s2', [1,1,2,-3,2]),
        # Band generators: s_ij = s_i...s_{j-1} s_j^2 s_{j-1}^-1...s_i^-1
        # s_{13} = s1.s2.s3^2.s2^-1.s1^-1
        ('s_{13}.s2', [1,2,3,3,-2,-1,2]),
        ('s_{13}.s2^-1', [1,2,3,3,-2,-1,-2]),
        ('s_{13}^2', [1,2,3,3,-2,-1]*2),
        # Braids from fibered knots
        # Trefoil closure: (s1)^3 in B2, but in B4: s1^3.s3^k
        ('s1^3.s3', [1,1,1,3]),
        ('s1^3.s3^2', [1,1,1,3,3]),
        ('s1^3.s3^-1', [1,1,1,-3]),
        # Figure-8 knot: s1.s2^-1.s1.s2^-1 in B3, extend to B4
        ('s1.s2^-1.s1.s2^-1.s3', [1,-2,1,-2,3]),
        ('s1.s2^-1.s1.s2^-1.s3^2', [1,-2,1,-2,3,3]),
        ('s1.s2^-1.s1.s2^-1.s3^-1', [1,-2,1,-2,-3]),
        # Almost-reducible: products of Dehn twist-like elements
        # T_1 ~ s1^2, T_2 ~ s2^2, T_3 ~ s3^2
        ('s1^2.s2^2.s3^2', [1,1,2,2,3,3]),
        ('s1^2.s2^-2.s3^2', [1,1,-2,-2,3,3]),
        ('s1^4.s2^2.s3^2', [1,1,1,1,2,2,3,3]),
        ('s1^2.s2^4.s3^2', [1,1,2,2,2,2,3,3]),
        ('s1^2.s2^2.s3^4', [1,1,2,2,3,3,3,3]),
    ]

    for name, word in braids:
        Mb = eval_word(word, bg, bgi)
        Ml = eval_word(word, lg, lgi)

        pts_b = newton_polytope_xq(Mb)
        pts_l = newton_polytope_xq(Ml)

        hull_b = convex_hull_2d(list(pts_b))
        hull_l = convex_hull_2d(list(pts_l))
        hull_l_s = scale_polytope(hull_l, scale)

        # Scale Lawrence monomials by 1/4
        pts_l_s = set((Rational(p[0], 4), Rational(p[1], 4)) for p in pts_l)

        hull_match = (hull_b == hull_l_s)
        monom_match = (pts_b == pts_l_s)

        if hull_match and monom_match:
            print(f'{name:<45} hull=YES  monoms=YES  |B|={len(pts_b)} |L/4|={len(pts_l_s)}')
            counts[0] += 1
        elif hull_match and not monom_match:
            print(f'{name:<45} hull=YES  monoms=*** NO ***  |B|={len(pts_b)} |L/4|={len(pts_l_s)}')
            print(f'  Burau only:    {sorted(pts_b - pts_l_s)}')
            print(f'  Lawrence only: {sorted(pts_l_s - pts_b)}')
            counts[0] += 1
            monomial_mismatches += 1
        else:
            print(f'{name:<45} hull=*** NO ***  |B|={len(pts_b)} |L/4|={len(pts_l_s)}')
            print(f'  NP(Burau):      {hull_b}')
            print(f'  NP(Lawrence)/4: {hull_l_s}')
            counts[1] += 1
        sys.stdout.flush()

    print(f'\n=== Hull: {counts[0]} passed, {counts[1]} failed ===')
    print(f'=== Monomial mismatches (hull agreed but monomials differed): {monomial_mismatches} ===')
