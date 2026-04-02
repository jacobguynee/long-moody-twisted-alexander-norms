#!/usr/bin/env python3
"""
hunt_counterexample_B4.py

Hunt for braids in B_4 where NP(Lawrence k=2)/4 != NP(Burau).
Candidates: reducible braids, commutators, non-alternating mixed braids.

Uses symbolic computation with Q=2, a=3 generic substitution.
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
def mat_power(M, n):
    if n == 0: return eye(M.shape[0])
    if n < 0: M = cancel(M.inv()); n = -n
    result = eye(M.shape[0]); base = M
    while n > 0:
        if n % 2 == 1: result = cancel(result * base)
        base = cancel(base * base); n //= 2
    return result

def test_braid(name, M_b, M_l, scale=Rational(1,4)):
    pts_b = newton_polytope_xq(M_b); hull_b = convex_hull_2d(list(pts_b))
    pts_l = newton_polytope_xq(M_l); hull_l = convex_hull_2d(list(pts_l))
    hull_l_scaled = scale_polytope(hull_l, scale)
    match = (hull_b == hull_l_scaled)
    if match:
        print(f'{name:<45} YES  {hull_b}')
    else:
        print(f'{name:<45} *** NO ***')
        print(f'  NP(Burau):      {hull_b}')
        print(f'  NP(Lawrence)/4: {hull_l_scaled}')
        print(f'  NP(Lawrence):   {hull_l}')
    sys.stdout.flush()
    return match


if __name__ == '__main__':
    print('Building representations...'); sys.stdout.flush()
    b = build_burau_B4(); bi = [cancel(g.inv()) for g in b]
    l = build_lawrence_B4(); li = [cancel(g.inv()) for g in l]
    print('Done.\n'); sys.stdout.flush()

    # Precompute generators and useful elements
    s1b, s2b, s3b = b[0], b[1], b[2]
    s1i, s2i, s3i = bi[0], bi[1], bi[2]
    s1l, s2l, s3l = l[0], l[1], l[2]
    s1li, s2li, s3li = li[0], li[1], li[2]

    scale = Rational(1, 4)
    counts = [0, 0]  # [pass, fail]

    def run(name, Mb, Ml):
        if test_braid(name, Mb, Ml, scale):
            counts[0] += 1
        else:
            counts[1] += 1

    # === FAMILY 1: Reducible braids f=s1, g=s3 (commuting!) ===
    print('=== Reducible: f=s1, g=s3 (commuting) ===')
    sys.stdout.flush()
    for aa, bb in [(3,2), (2,3), (5,1), (1,5), (3,-2), (-2,3), (-3,-2),
                   (4,3), (5,5), (7,-3), (-5,5), (8,2), (2,8), (10,1)]:
        Mb = cancel(mat_power(s1b, aa) * mat_power(s3b, bb))
        Ml = cancel(mat_power(s1l, aa) * mat_power(s3l, bb))
        run(f's1^{aa} s3^{bb}', Mb, Ml)

    # === FAMILY 2: Commutators ===
    print('\n=== Commutators [s_i, s_j] and powers ===')
    sys.stdout.flush()
    # [s1, s2] = s1 s2 s1^-1 s2^-1
    comm12_b = cancel(s1b * s2b * s1i * s2i)
    comm12_l = cancel(s1l * s2l * s1li * s2li)
    run('[s1,s2]', comm12_b, comm12_l)

    # [s2, s3]
    comm23_b = cancel(s2b * s3b * s2i * s3i)
    comm23_l = cancel(s2l * s3l * s2li * s3li)
    run('[s2,s3]', comm23_b, comm23_l)

    # [s1, s2]^k
    for k in [2, 3, 4, 5]:
        run(f'[s1,s2]^{k}', mat_power(comm12_b, k), mat_power(comm12_l, k))

    # [s1,s2] * [s2,s3]
    Mb = cancel(comm12_b * comm23_b)
    Ml = cancel(comm12_l * comm23_l)
    run('[s1,s2][s2,s3]', Mb, Ml)

    # [s1,s2]^a * [s2,s3]^b
    for aa, bb in [(2,1), (1,2), (3,2), (2,3), (2,-1), (-1,2), (3,-1)]:
        Mb = cancel(mat_power(comm12_b, aa) * mat_power(comm23_b, bb))
        Ml = cancel(mat_power(comm12_l, aa) * mat_power(comm23_l, bb))
        run(f'[s1,s2]^{aa} [s2,s3]^{bb}', Mb, Ml)

    # === FAMILY 3: Non-alternating mixed braids ===
    print('\n=== Non-alternating mixed braids ===')
    sys.stdout.flush()

    def eval_word(word, gens, inv_gens):
        M = eye(gens[0].shape[0])
        for w in word:
            if w > 0: M = M * gens[w-1]
            else: M = M * inv_gens[-w-1]
        return cancel(M)

    mixed_braids = [
        # s1 s2 s1^-1 s3 (non-alternating: pos,pos,neg,pos)
        ('s1.s2.s1^-1.s3',            [1, 2, -1, 3]),
        # More complex non-alternating
        ('s1.s2^-1.s1.s3^-1',         [1, -2, 1, -3]),
        ('s1^-1.s2.s3^-1.s2',         [-1, 2, -3, 2]),
        ('s1.s2.s3^-1.s1^-1.s2.s3',   [1, 2, -3, -1, 2, 3]),
        ('s1^-1.s2^-1.s3.s1.s2.s3^-1',[-1, -2, 3, 1, 2, -3]),
        # Reducible-ish
        ('s1^2.s3^2.s2',              [1, 1, 3, 3, 2]),
        ('s1.s3.s1.s3.s2',            [1, 3, 1, 3, 2]),
        ('s2.s1^2.s2^-1.s3^2',        [2, 1, 1, -2, 3, 3]),
        # Long non-alternating
        ('s1.s2.s3.s1^-1.s2^-1.s3^-1.s1.s2',
                                       [1, 2, 3, -1, -2, -3, 1, 2]),
        ('s1^2.s2^-1.s3.s1^-1.s2.s3^-2',
                                       [1, 1, -2, 3, -1, 2, -3, -3]),
        ('s1.s3.s2^-1.s1^-1.s3^-1.s2.s1.s3',
                                       [1, 3, -2, -1, -3, 2, 1, 3]),
    ]

    for name, word in mixed_braids:
        Mb = eval_word(word, b, bi)
        Ml = eval_word(word, l, li)
        run(name, Mb, Ml)

    # === FAMILY 4: Products involving Garside ===
    print('\n=== Garside-related ===')
    sys.stdout.flush()
    # Half twist Delta = s1.s2.s1.s3.s2.s1 (Garside element of B4... actually
    # Delta = s1 s2 s3 s1 s2 s1 for standard Garside)
    # Standard Garside for B4: Delta = (s1)(s2 s1)(s3 s2 s1)
    delta_word = [1, 2, 1, 3, 2, 1]
    delta_b = eval_word(delta_word, b, bi)
    delta_l = eval_word(delta_word, l, li)

    run('Delta (Garside)', delta_b, delta_l)
    run('Delta^2 (full twist)', cancel(delta_b*delta_b), cancel(delta_l*delta_l))
    run('Delta^3', mat_power(delta_b, 3), mat_power(delta_l, 3))

    # Delta * s1
    run('Delta.s1', cancel(delta_b*s1b), cancel(delta_l*s1l))
    # Delta * s1^-1
    run('Delta.s1^-1', cancel(delta_b*s1i), cancel(delta_l*s1li))
    # s1 * Delta^-1
    delta_bi = cancel(delta_b.inv()); delta_li = cancel(delta_l.inv())
    run('s1.Delta^-1', cancel(s1b*delta_bi), cancel(s1l*delta_li))

    # Delta * s1^a * s3^b (reducible part conjugated by Garside)
    for aa, bb in [(2,1), (1,2), (3,-1), (-2,3)]:
        Mb = cancel(delta_b * mat_power(s1b, aa) * mat_power(s3b, bb))
        Ml = cancel(delta_l * mat_power(s1l, aa) * mat_power(s3l, bb))
        run(f'Delta.s1^{aa}.s3^{bb}', Mb, Ml)

    print(f'\n=== {counts[0]} passed, {counts[1]} failed out of {counts[0]+counts[1]} ===')
