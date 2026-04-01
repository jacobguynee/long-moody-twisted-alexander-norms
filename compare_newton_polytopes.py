#!/usr/bin/env python3
"""
compare_newton_polytopes.py

For alternating 3-braids (words in sigma_1 and sigma_2^{-1}), compare
the Newton polytope of the characteristic polynomial of:
  - Reduced Burau (2x2, parameter q)
  - Iterated reduced Long-Moody / Lawrence k=2 (6x6, parameters Q, q, a)
    scaled by 1/3

Newton polytope is with respect to (x, q). Coefficients in Q(Q, a)
are just checked nonzero.
"""

from sympy import symbols, eye, zeros, Matrix, cancel, Poly, Rational, fraction
import sys

Q, q, a, x = symbols('Q q a x')


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
# Build representations for B_3
# ============================================================

def build_burau_B3():
    """Reduced Burau for B_3: rho(g_i)=q, rho(sigma_i)=1."""
    rho_g = [Matrix([[q]]), Matrix([[q]]), Matrix([[q]])]
    rho_sigma = [Matrix([[1]]), Matrix([[1]])]
    gens = reduced_lm_colored(rho_g, rho_sigma)
    return gens[0], gens[1]

def build_lawrence_B3():
    """Iterated reduced Long-Moody for B_3 (Lawrence k=2). Parameters: Q, q, a."""
    rho_g_seed = [Matrix([[Q]]), Matrix([[Q]]), Matrix([[Q]]), Matrix([[a]])]
    rho_s_seed = [Matrix([[1]]), Matrix([[1]]), Matrix([[1]])]
    L3Pre = reduced_lm_colored(rho_g_seed, rho_s_seed, [1, 1, 1, 2])
    L3PreFree = free_group_gens(L3Pre, [1, 1, 1, 2])
    rho_g2 = [g * q for g in L3PreFree]
    L3 = reduced_lm_colored(rho_g2, L3Pre[0:2])
    return L3[0], L3[1]


# ============================================================
# Newton polytope in (x, q)
# ============================================================

def newton_polytope_xq(M):
    """Compute Newton polytope of charpoly of M in (x, q) coordinates.
    Uses Berkowitz algorithm via charpoly for speed.
    Returns set of (x_deg, q_deg) points with nonzero coefficients in Q(Q,a)."""
    n = M.shape[0]
    # charpoly returns polynomial in x with coefficients that are rational functions of q, Q, a
    cp = M.charpoly(x)
    # cp is a PurePoly; get it as an expression
    cp_expr = cp.as_expr()
    cp_expr = cancel(cp_expr)

    num, den = fraction(cp_expr)

    # Collect numerator as polynomial in x and q
    p = Poly(num, x, q)
    monoms_num = p.monoms()
    coeffs_num = p.coeffs()

    # q-degree of denominator
    if den.has(q):
        pd = Poly(den, q)
        q_shift = pd.degree()
    else:
        q_shift = 0

    points = set()
    for (xd, qd), coeff in zip(monoms_num, coeffs_num):
        if coeff != 0:
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
    print('Building Burau generators for B_3 (param: q)...')
    b1, b2 = build_burau_B3()
    b1_inv = cancel(b1.inv()); b2_inv = cancel(b2.inv())
    print(f'  Burau: {b1.shape[0]}x{b1.shape[1]}')

    print('Building Lawrence k=2 generators for B_3 (params: Q, q, a)...')
    sys.stdout.flush()
    l1, l2 = build_lawrence_B3()
    l1_inv = cancel(l1.inv()); l2_inv = cancel(l2.inv())
    print(f'  Lawrence: {l1.shape[0]}x{l1.shape[1]}')
    sys.stdout.flush()

    b_gens = [b1, b2]; b_inv = [b1_inv, b2_inv]
    l_gens = [l1, l2]; l_inv = [l1_inv, l2_inv]

    # Alternating 3-braids: words in sigma_1 (positive) and sigma_2^{-1}
    braids = [
        ('s1.s2^-1',              [1, -2]),
        ('s1^2.s2^-1',            [1, 1, -2]),
        ('s1.s2^-2',              [1, -2, -2]),
        ('s1^2.s2^-2',            [1, 1, -2, -2]),
        ('s1^3.s2^-1',            [1, 1, 1, -2]),
        ('s1.s2^-3',              [1, -2, -2, -2]),
        ('s1.s2^-1.s1.s2^-1',    [1, -2, 1, -2]),
    ]

    print(f'\nNewton polytopes in (x, q). Coefficients in Q(Q, a) checked nonzero.\n')
    sys.stdout.flush()

    for name, word in braids:
        print(f'--- {name} ---')
        sys.stdout.flush()

        # Burau
        M_b = eval_braid_word(word, b_gens, b_inv)
        pts_b = newton_polytope_xq(M_b)
        hull_b = convex_hull_2d(list(pts_b))

        # Lawrence
        M_l = eval_braid_word(word, l_gens, l_inv)
        print(f'  (computing 6x6 charpoly...)')
        sys.stdout.flush()
        pts_l = newton_polytope_xq(M_l)
        hull_l = convex_hull_2d(list(pts_l))

        hull_l_scaled = scale_polytope(hull_l, Rational(1, 3))

        match = (hull_b == hull_l_scaled)

        print(f'  NP(Burau):       {hull_b}')
        print(f'  NP(Lawrence)/3:  {hull_l_scaled}')
        print(f'  Match: {"YES" if match else "NO"}')
        if not match:
            print(f'  NP(Lawrence):    {hull_l}')
        print()
        sys.stdout.flush()
