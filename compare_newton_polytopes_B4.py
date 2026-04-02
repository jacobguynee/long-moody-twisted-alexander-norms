#!/usr/bin/env python3
"""
compare_newton_polytopes_B4.py

For alternating 4-braids (words in sigma_1, sigma_2^{-1}, sigma_3),
compare Newton polytopes of characteristic polynomials:
  - Reduced Burau (3x3, parameter q)
  - Lawrence k=2 (12x12, parameters Q, q, a)
    scaled by 1/4  (= dim_Burau / dim_Lawrence = 3/12)

To keep computation fast, substitute generic numerical values for Q, a
(Q=2, a=3) so charpolys are polynomials in (x, q) with rational coefficients.
Newton polytope is independent of generic parameter choices.
"""

from sympy import symbols, eye, zeros, Matrix, cancel, Poly, Rational, fraction
import sys

q, x = symbols('q x')

# Generic numerical values for Q and a — chosen to avoid accidental cancellation
Q_val = Rational(2)
a_val = Rational(3)


# ============================================================
# Reduced Long-Moody (exact port from .m files)
# ============================================================

def build_Ni(i, gi, n, d):
    nb = n - 1
    Id = eye(d); Zd = zeros(d)
    blocks = [[Id.copy() if r == c else Zd.copy() for c in range(nb)] for r in range(nb)]
    if i == 0:
        blocks[0][0] = -gi; blocks[0][1] = Zd
        blocks[1][0] = Id;  blocks[1][1] = Id
    elif i == n - 2:
        blocks[nb-2][nb-2] = Id;  blocks[nb-2][nb-1] = gi
        blocks[nb-1][nb-2] = Zd;  blocks[nb-1][nb-1] = -gi
    else:
        p = i - 1
        blocks[p][p] = Id;    blocks[p][p+1] = gi;    blocks[p][p+2] = Zd
        blocks[p+1][p] = Zd;  blocks[p+1][p+1] = -gi; blocks[p+1][p+2] = Zd
        blocks[p+2][p] = Zd;  blocks[p+2][p+1] = Id;  blocks[p+2][p+2] = Id
    return blocks

def left_multiply_blocks(blocks, M, nb):
    return [[M * blocks[r][c] for c in range(nb)] for r in range(nb)]

def multiply_block_matrices(A, B, nb, d):
    C = [[zeros(d) for _ in range(nb)] for _ in range(nb)]
    for r in range(nb):
        for c in range(nb):
            acc = zeros(d)
            for k in range(nb):
                acc = acc + A[r][k] * B[k][c]
            C[r][c] = acc
    return C

def assemble(blocks, nb, d):
    N = nb * d
    M = zeros(N)
    for r in range(nb):
        for c in range(nb):
            M[r*d:(r+1)*d, c*d:(c+1)*d] = blocks[r][c]
    return M

def reduced_lm_colored(rho_g, rho_sigma, partition=None):
    n = len(rho_g); num_sigmas = len(rho_sigma); d = rho_g[0].shape[0]
    if partition is None:
        partition = [1] * n
    is_boundary = [partition[i] != partition[i+1] for i in range(num_sigmas)]
    nb = n - 1; rho_plus = []
    for i in range(num_sigmas):
        if not is_boundary[i]:
            blocks = build_Ni(i, rho_g[i], n, d)
            blocks = left_multiply_blocks(blocks, rho_sigma[i], nb)
            rho_plus.append(cancel(assemble(blocks, nb, d)))
        else:
            gi_conj = cancel(rho_g[i] * rho_g[i+1] * rho_g[i].inv())
            blocks_N  = build_Ni(i, rho_g[i], n, d)
            blocks_Np = build_Ni(i, gi_conj, n, d)
            blocks_prod = multiply_block_matrices(blocks_Np, blocks_N, nb, d)
            si_sq = rho_sigma[i] * rho_sigma[i]
            blocks_prod = left_multiply_blocks(blocks_prod, si_sq, nb)
            rho_plus.append(cancel(assemble(blocks_prod, nb, d)))
    return rho_plus

def free_group_gens(rho_sigma, partition=None):
    n = len(rho_sigma)
    if partition is None:
        partition = [1] * (n + 1)
    last_is_boundary = (partition[n-1] != partition[n])
    gn = rho_sigma[n-1] if last_is_boundary else rho_sigma[n-1] * rho_sigma[n-1]
    rho_g = [None] * n
    rho_g[n-1] = cancel(gn)
    for i in range(n-2, -1, -1):
        rho_g[i] = cancel(rho_sigma[i].inv() * rho_g[i+1] * rho_sigma[i])
    return rho_g


# ============================================================
# Build representations for B_4
# ============================================================

def build_burau_B4():
    """Reduced Burau for B_4: rho(g_i)=q, rho(sigma_i)=1. 3x3 matrices."""
    rho_g = [Matrix([[q]])] * 4
    rho_sigma = [Matrix([[1]])] * 3
    gens = reduced_lm_colored(rho_g, rho_sigma)
    return gens  # [sigma_1, sigma_2, sigma_3]

def build_lawrence_B4():
    """Iterated reduced Long-Moody for B_4 (Lawrence k=2).
    Q and a are numerical; q is symbolic. Output: 12x12 matrices.
    Exact translation of verify_lawrence_conjugacy_B4.py."""
    # Step 1: F_5 x B_5 with partition [1 1 1 1 2]
    rho_g_seed = [Matrix([[Q_val]])] * 4 + [Matrix([[a_val]])]
    rho_s_seed = [Matrix([[1]])] * 4
    L4Pre = reduced_lm_colored(rho_g_seed, rho_s_seed, [1, 1, 1, 1, 2])

    # Step 2: recover free group generators
    L4PreFree = free_group_gens(L4Pre, [1, 1, 1, 1, 2])

    # Step 3: scale by q, second iteration
    rho_g2 = [g * q for g in L4PreFree]
    L4 = reduced_lm_colored(rho_g2, L4Pre[0:3])
    return L4  # [sigma_1, sigma_2, sigma_3]


# ============================================================
# Newton polytope in (x, q)
# ============================================================

def newton_polytope_xq(M):
    """Newton polytope of charpoly of M in (x, q) coordinates."""
    cp = M.charpoly(x)
    cp_expr = cancel(cp.as_expr())
    num, den = fraction(cp_expr)
    p = Poly(num, x, q)
    monoms_num = p.monoms()
    if den.has(q):
        q_shift = Poly(den, q).degree()
    else:
        q_shift = 0
    points = set()
    for (xd, qd) in monoms_num:
        points.add((xd, qd - q_shift))
    return points


def convex_hull_2d(points):
    points = sorted(set(points))
    if len(points) <= 2:
        return points
    def cross(O, A, B):
        return (A[0]-O[0])*(B[1]-O[1]) - (A[1]-O[1])*(B[0]-O[0])
    lower = []
    for p in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return sorted(set(lower[:-1] + upper[:-1]))


def scale_polytope(vertices, s):
    return sorted(set((Rational(v[0]) * s, Rational(v[1]) * s) for v in vertices))


# ============================================================
# Braid word evaluation
# ============================================================

def eval_braid_word(word, gens, inv_gens):
    n = gens[0].shape[0]
    M = eye(n)
    for w in word:
        if w > 0:
            M = M * gens[w - 1]
        else:
            M = M * inv_gens[-w - 1]
    return cancel(M)


# ============================================================
# Main
# ============================================================

if __name__ == '__main__':
    print(f'Generic numerical values: Q = {Q_val}, a = {a_val}')
    print(f'Symbolic variable: q\n')

    print('Building Burau generators for B_4 (3x3)...')
    sys.stdout.flush()
    b_gens_list = build_burau_B4()
    b_inv_list = [cancel(g.inv()) for g in b_gens_list]
    print(f'  Burau: {b_gens_list[0].shape[0]}x{b_gens_list[0].shape[1]}, {len(b_gens_list)} generators')

    print('Building Lawrence k=2 generators for B_4 (12x12)...')
    sys.stdout.flush()
    l_gens_list = build_lawrence_B4()
    l_inv_list = [cancel(g.inv()) for g in l_gens_list]
    print(f'  Lawrence: {l_gens_list[0].shape[0]}x{l_gens_list[0].shape[1]}, {len(l_gens_list)} generators')
    sys.stdout.flush()

    # Alternating 4-braids: words in sigma_1, sigma_2^{-1}, sigma_3
    # (alternating sign pattern by generator index)
    # All use all 3 generators.
    braids = [
        # Basic mixed words
        ('s1.s2^-1.s3',                    [1, -2, 3]),
        ('s3.s2^-1.s1',                    [3, -2, 1]),
        ('s1.s2^-1.s3.s2^-1',              [1, -2, 3, -2]),
        ('s1^2.s2^-1.s3',                  [1, 1, -2, 3]),
        ('s1.s2^-2.s3',                    [1, -2, -2, 3]),
        ('s1.s2^-1.s3^2',                  [1, -2, 3, 3]),
        # Higher complexity
        ('s1^2.s2^-2.s3^2',                [1, 1, -2, -2, 3, 3]),
        ('s1.s2^-1.s3.s1.s2^-1.s3',        [1, -2, 3, 1, -2, 3]),
        ('s1^2.s2^-1.s3.s2^-1.s1',          [1, 1, -2, 3, -2, 1]),
        ('s1.s2^-2.s3^2.s2^-1.s1',          [1, -2, -2, 3, 3, -2, 1]),
        ('s1^3.s2^-2.s3',                  [1, 1, 1, -2, -2, 3]),
        ('s1.s2^-3.s3^2',                  [1, -2, -2, -2, 3, 3]),
        # Longer words
        ('s1.s2^-1.s3.s1.s2^-1.s3.s1',      [1, -2, 3, 1, -2, 3, 1]),
        ('s1^2.s2^-1.s3^2.s2^-1.s1^2',      [1, 1, -2, 3, 3, -2, 1, 1]),
        ('s3.s2^-1.s1^2.s2^-1.s3.s2^-1',    [3, -2, 1, 1, -2, 3, -2]),
        ('s1.s2^-1.s3.s2^-1.s1.s2^-1.s3',  [1, -2, 3, -2, 1, -2, 3]),
    ]

    # Try both 1/3 and 1/4 scaling to see which (if either) works
    # n=4 so dim ratio = 12/3 = 4, expect 1/4
    scale = Rational(1, 4)

    print(f'\nNewton polytopes in (x, q). Scaling Lawrence by {scale}.\n')
    sys.stdout.flush()

    for name, word in braids:
        print(f'--- {name} ---', end='  ')
        sys.stdout.flush()

        # Burau
        M_b = eval_braid_word(word, b_gens_list, b_inv_list)
        pts_b = newton_polytope_xq(M_b)
        hull_b = convex_hull_2d(list(pts_b))

        # Lawrence
        M_l = eval_braid_word(word, l_gens_list, l_inv_list)
        pts_l = newton_polytope_xq(M_l)
        hull_l = convex_hull_2d(list(pts_l))

        hull_l_scaled = scale_polytope(hull_l, scale)
        match = (hull_b == hull_l_scaled)

        if match:
            print(f'Match: YES   NP(Burau) = {hull_b}')
        else:
            print(f'Match: NO')
            print(f'    NP(Burau):          {hull_b}')
            print(f'    NP(Lawrence)/{scale}: {hull_l_scaled}')
            print(f'    NP(Lawrence) raw:   {hull_l}')
        sys.stdout.flush()

    print(f'\n--- Done ---')
