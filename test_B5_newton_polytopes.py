#!/usr/bin/env python3
"""
test_B5_newton_polytopes.py

Test NP(Lawrence k=2)/5 = NP(Burau) for B_5 braids.
B5 has generators s1, s2, s3, s4.
Reduced Burau is 4x4, Lawrence k=2 is 4*5=20x20.
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

def build_burau_B5():
    return reduced_lm_colored([Matrix([[q]])]*5, [Matrix([[1]])]*4)

def build_lawrence_B5():
    L5Pre = reduced_lm_colored(
        [Matrix([[Q_val]])]*5 + [Matrix([[a_val]])],
        [Matrix([[1]])]*5,
        [1,1,1,1,1,2]
    )
    L5PreFree = free_group_gens(L5Pre, [1,1,1,1,1,2])
    return reduced_lm_colored([g*q for g in L5PreFree], L5Pre[0:4])

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

if __name__ == '__main__':
    print('Building B5 Burau (4x4)...'); sys.stdout.flush()
    bg = build_burau_B5(); bgi = [cancel(g.inv()) for g in bg]
    print(f'  Burau gen sizes: {[g.shape for g in bg]}')

    print('Building B5 Lawrence k=2...'); sys.stdout.flush()
    lg = build_lawrence_B5(); lgi = [cancel(g.inv()) for g in lg]
    print(f'  Lawrence gen sizes: {[g.shape for g in lg]}')
    print('Done.\n'); sys.stdout.flush()

    scale = Rational(1, 5)
    counts = [0, 0]

    def run(name, word):
        Mb = eval_word(word, bg, bgi)
        Ml = eval_word(word, lg, lgi)
        hull_b = convex_hull_2d(list(newton_polytope_xq(Mb)))
        hull_l = convex_hull_2d(list(newton_polytope_xq(Ml)))
        hull_l_s = scale_polytope(hull_l, scale)
        match = (hull_b == hull_l_s)
        if match:
            print(f'{name:<45} YES  {hull_b}')
            counts[0] += 1
        else:
            print(f'{name:<45} *** NO ***')
            print(f'  NP(Burau):      {hull_b}')
            print(f'  NP(Lawrence)/5: {hull_l_s}')
            counts[1] += 1
        sys.stdout.flush()

    # Simple alternating braids using all 4 generators
    print('=== B5 alternating braids ==='); sys.stdout.flush()
    braids = [
        ('s1.s2^-1.s3.s4^-1', [1,-2,3,-4]),
        ('s1.s2^-1.s3.s4^-1.s1', [1,-2,3,-4,1]),
        ('s1^2.s2^-1.s3.s4^-1', [1,1,-2,3,-4]),
        ('s1.s2^-1.s3^2.s4^-1', [1,-2,3,3,-4]),
        ('s1.s2^-1.s3.s4^-1.s3.s2^-1', [1,-2,3,-4,3,-2]),
        ('(s1.s2^-1)^2.s3.s4^-1', [1,-2,1,-2,3,-4]),
        ('s1.s2^-1.s3.s4^-1.s1.s2^-1', [1,-2,3,-4,1,-2]),
    ]
    for name, word in braids:
        run(name, word)

    # Mixed braids
    print('\n=== B5 mixed braids ==='); sys.stdout.flush()
    mixed = [
        ('s1.s2.s3.s4', [1,2,3,4]),
        ('s1.s2.s3.s4.s1', [1,2,3,4,1]),
        ('s1^2.s2.s3.s4', [1,1,2,3,4]),
        ('s1.s2^-1.s3^-1.s4', [1,-2,-3,4]),
        ('s1.s3.s2^-1.s4', [1,3,-2,4]),
        ('s1^-1.s2.s3^-1.s4.s2^-1', [-1,2,-3,4,-2]),
        ('s1.s2.s1^-1.s3.s4.s3^-1', [1,2,-1,3,4,-3]),
    ]
    for name, word in mixed:
        run(name, word)

    # Reducible and almost-reducible
    print('\n=== B5 reducible/almost-reducible ==='); sys.stdout.flush()
    reducible = [
        # s1,s2 commute with s4 (strands 1-3 decouple from strand 5)
        ('s1^2.s2.s4^3', [1,1,2,4,4,4]),
        ('s1.s2^-1.s4^2', [1,-2,4,4]),
        # Almost reducible: small s3 bridge
        ('s1^2.s2.s3.s4^2', [1,1,2,3,4,4]),
        ('s1^3.s3.s4^2', [1,1,1,3,4,4]),
        ('s1^2.s2^-1.s3.s4^-2', [1,1,-2,3,-4,-4]),
    ]
    for name, word in reducible:
        run(name, word)

    # Commutators
    print('\n=== B5 commutators ==='); sys.stdout.flush()
    comms = [
        ('[s1,s2]', [1,2,-1,-2]),
        ('[s1,s2].s3.s4', [1,2,-1,-2,3,4]),
        ('[s1,s2].[s3,s4]', [1,2,-1,-2,3,4,-3,-4]),
        ('[s2,s3]', [2,3,-2,-3]),
        ('[s1,s3].[s2,s4]', [1,3,-1,-3,2,4,-2,-4]),
    ]
    for name, word in comms:
        run(name, word)

    print(f'\n=== {counts[0]} passed, {counts[1]} failed out of {counts[0]+counts[1]} ===')
