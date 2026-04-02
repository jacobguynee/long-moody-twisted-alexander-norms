#!/usr/bin/env python3
"""
compare_newton_polytopes_B4_conjugated.py

For B_4 braids of the form f^a * g^b where:
  f = sigma_1
  g = sigma_2^2 * sigma_3 * sigma_2^{-2}
with a, b nonzero and |a|+|b| >= 10.

Compare NP(Lawrence k=2)/4 vs NP(Burau) in (x, q).
"""

from sympy import symbols, eye, zeros, Matrix, cancel, Poly, Rational, fraction
import sys

q, x = symbols('q x')
Q_val = Rational(2)
a_val = Rational(3)


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


def build_burau_B4():
    rho_g = [Matrix([[q]])] * 4
    rho_sigma = [Matrix([[1]])] * 3
    return reduced_lm_colored(rho_g, rho_sigma)

def build_lawrence_B4():
    rho_g_seed = [Matrix([[Q_val]])] * 4 + [Matrix([[a_val]])]
    rho_s_seed = [Matrix([[1]])] * 4
    L4Pre = reduced_lm_colored(rho_g_seed, rho_s_seed, [1, 1, 1, 1, 2])
    L4PreFree = free_group_gens(L4Pre, [1, 1, 1, 1, 2])
    rho_g2 = [g * q for g in L4PreFree]
    return reduced_lm_colored(rho_g2, L4Pre[0:3])


def newton_polytope_xq(M):
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


def mat_power(M, n):
    """Compute M^n for integer n (positive or negative)."""
    if n == 0:
        return eye(M.shape[0])
    if n < 0:
        M = cancel(M.inv())
        n = -n
    # Binary exponentiation
    result = eye(M.shape[0])
    base = M
    while n > 0:
        if n % 2 == 1:
            result = cancel(result * base)
        base = cancel(base * base)
        n //= 2
    return result


if __name__ == '__main__':
    print(f'Generic numerical values: Q = {Q_val}, a = {a_val}')
    print(f'f = sigma_1')
    print(f'g = sigma_2^2 * sigma_3 * sigma_2^{{-2}}\n')

    print('Building Burau generators for B_4 (3x3)...')
    sys.stdout.flush()
    b_gens = build_burau_B4()
    b_inv = [cancel(g.inv()) for g in b_gens]

    # Precompute f and g in Burau
    # f = sigma_1
    f_b = b_gens[0]
    f_b_inv = b_inv[0]
    # g = sigma_2^2 * sigma_3 * sigma_2^{-2}
    g_b = cancel(b_gens[1] * b_gens[1] * b_gens[2] * b_inv[1] * b_inv[1])
    g_b_inv = cancel(b_gens[1] * b_gens[1] * b_inv[2] * b_inv[1] * b_inv[1])
    print(f'  f_Burau, g_Burau precomputed')

    print('Building Lawrence k=2 generators for B_4 (12x12)...')
    sys.stdout.flush()
    l_gens = build_lawrence_B4()
    l_inv = [cancel(g.inv()) for g in l_gens]

    # Precompute f and g in Lawrence
    f_l = l_gens[0]
    f_l_inv = l_inv[0]
    g_l = cancel(l_gens[1] * l_gens[1] * l_gens[2] * l_inv[1] * l_inv[1])
    g_l_inv = cancel(l_gens[1] * l_gens[1] * l_inv[2] * l_inv[1] * l_inv[1])
    print(f'  f_Lawrence, g_Lawrence precomputed')
    sys.stdout.flush()

    # f^a g^b with |a|+|b| >= 10, a,b nonzero
    test_cases = [
        # Positive a, positive b
        (5, 5),
        (6, 4),
        (4, 6),
        (7, 3),
        (3, 7),
        (8, 2),
        (2, 8),
        (9, 1),
        (1, 9),
        # Mixed signs
        (5, -5),
        (-5, 5),
        (6, -4),
        (-4, 6),
        (7, -3),
        (-3, 7),
        (-6, -4),
        (-4, -6),
        # Higher complexity
        (8, -3),
        (-3, 8),
        (10, 1),
        (1, 10),
        (-5, -6),
    ]

    scale = Rational(1, 4)
    print(f'\nTesting f^a * g^b with |a|+|b| >= 10. Scaling Lawrence by {scale}.\n')
    sys.stdout.flush()

    pass_count = 0
    fail_count = 0

    for (aa, bb) in test_cases:
        label = f'f^{aa} g^{bb}'
        print(f'{label:<16}', end='')
        sys.stdout.flush()

        # Burau: f^a * g^b
        M_b = cancel(mat_power(f_b, aa) * mat_power(g_b, bb))
        pts_b = newton_polytope_xq(M_b)
        hull_b = convex_hull_2d(list(pts_b))

        # Lawrence: f^a * g^b
        M_l = cancel(mat_power(f_l, aa) * mat_power(g_l, bb))
        pts_l = newton_polytope_xq(M_l)
        hull_l = convex_hull_2d(list(pts_l))

        hull_l_scaled = scale_polytope(hull_l, scale)
        match = (hull_b == hull_l_scaled)

        if match:
            print(f'YES   NP(Burau) = {hull_b}')
            pass_count += 1
        else:
            print(f'NO')
            print(f'    NP(Burau):          {hull_b}')
            print(f'    NP(Lawrence)/4:     {hull_l_scaled}')
            fail_count += 1
        sys.stdout.flush()

    print(f'\n=== {pass_count} passed, {fail_count} failed out of {pass_count + fail_count} ===')
