#!/usr/bin/env python3
"""
hunt_counterexample_B4_v2.py

Second round: try braids that stress-test cancellations.
- Near-identity braids (conjugates of simple braids)
- Pseudo-Anosov with small dilatation
- Elements with zero writhe (equal pos/neg crossings)
- Random walks in B4
"""

from sympy import symbols, eye, zeros, Matrix, cancel, Poly, Rational, fraction
import sys
import random

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

def test_braid(name, M_b, M_l, scale=Rational(1,4)):
    pts_b = newton_polytope_xq(M_b); hull_b = convex_hull_2d(list(pts_b))
    pts_l = newton_polytope_xq(M_l); hull_l = convex_hull_2d(list(pts_l))
    hull_l_scaled = scale_polytope(hull_l, scale)
    match = (hull_b == hull_l_scaled)
    if match:
        print(f'{name:<55} YES  {hull_b}')
    else:
        print(f'{name:<55} *** NO ***')
        print(f'  NP(Burau):      {hull_b}')
        print(f'  NP(Lawrence)/4: {hull_l_scaled}')
    sys.stdout.flush()
    return match


if __name__ == '__main__':
    print('Building representations...'); sys.stdout.flush()
    b = build_burau_B4(); bi = [cancel(g.inv()) for g in b]
    l = build_lawrence_B4(); li = [cancel(g.inv()) for g in l]
    print('Done.\n'); sys.stdout.flush()

    scale = Rational(1, 4)
    counts = [0, 0]

    def run(name, word):
        Mb = eval_word(word, b, bi)
        Ml = eval_word(word, l, li)
        if test_braid(name, Mb, Ml, scale):
            counts[0] += 1
        else:
            counts[1] += 1

    # === Zero-writhe braids (equal positive and negative crossings) ===
    print('=== Zero-writhe braids ===')
    sys.stdout.flush()
    zero_writhe = [
        ('s1.s1^-1.s2.s3.s2^-1.s3^-1',      [1,-1,2,3,-2,-3]),
        ('s1.s2.s3.s1^-1.s2^-1.s3^-1',       [1,2,3,-1,-2,-3]),
        ('s1.s2.s1^-1.s3.s2^-1.s3^-1',       [1,2,-1,3,-2,-3]),
        ('s1.s3.s2^-1.s3^-1.s2.s1^-1',       [1,3,-2,-3,2,-1]),
        ('s2.s1.s3.s2^-1.s1^-1.s3^-1',       [2,1,3,-2,-1,-3]),
        ('s1.s2.s3^-1.s2^-1.s3.s1^-1',       [1,2,-3,-2,3,-1]),
        # longer zero-writhe
        ('s1^2.s2.s3.s1^-2.s2^-1.s3^-1',     [1,1,2,3,-1,-1,-2,-3]),
        ('s1.s2^2.s3.s1^-1.s2^-2.s3^-1',     [1,2,2,3,-1,-2,-2,-3]),
        ('s1.s2.s3^2.s1^-1.s2^-1.s3^-2',     [1,2,3,3,-1,-2,-3,-3]),
        ('s1.s2.s3.s2.s1^-1.s2^-1.s3^-1.s2^-1',
                                               [1,2,3,2,-1,-2,-3,-2]),
    ]
    for name, word in zero_writhe:
        run(name, word)

    # === Conjugates: w.s_i.w^-1 for complex w ===
    print('\n=== Conjugates of generators ===')
    sys.stdout.flush()
    conjugators = [
        ('s2.s3.s1.s2^-1',               [2,3,1,-2]),
        ('s1.s2.s3.s2^-1.s1^-1',         [1,2,3,-2,-1]),
        ('s3.s1.s2.s1^-1.s3^-1',         [3,1,2,-1,-3]),
        ('s1^-1.s2.s3.s2^-1.s1',         [-1,2,3,-2,1]),
    ]
    for conj_name, conj_word in conjugators:
        for gen in [1, 2, 3]:
            word = conj_word + [gen] + [-w for w in reversed(conj_word)]
            name = f'({conj_name}).s{gen}.({conj_name})^-1'
            run(name, word)

    # === Near-identity: products of conjugate pairs ===
    print('\n=== Products of conjugate pairs (near identity) ===')
    sys.stdout.flush()
    near_id = [
        # w.s1.w^-1 . v.s1^-1.v^-1 for different w,v
        ('w.s1.w^-1.v.s1^-1.v^-1 (w=s2,v=s3)',
         [2,1,-2, 3,-1,-3]),
        ('w.s1.w^-1.v.s1^-1.v^-1 (w=s2s3,v=s3s2)',
         [2,3,1,-3,-2, 3,2,-1,-2,-3]),
        ('w.s2.w^-1.v.s2^-1.v^-1 (w=s1s3,v=s3s1)',
         [1,3,2,-3,-1, 3,1,-2,-1,-3]),
        ('w.s1^2.w^-1.v.s1^-2.v^-1 (w=s2,v=s3)',
         [2,1,1,-2, 3,-1,-1,-3]),
        ('w.s1^3.w^-1.v.s1^-3.v^-1 (w=s2,v=s3)',
         [2,1,1,1,-2, 3,-1,-1,-1,-3]),
    ]
    for name, word in near_id:
        run(name, word)

    # === Random walks (seeded for reproducibility) ===
    print('\n=== Pseudo-random braids (seed=42) ===')
    sys.stdout.flush()
    random.seed(42)
    for trial in range(20):
        length = random.choice([6, 8, 10])
        word = []
        used_gens = set()
        while len(word) < length or len(used_gens) < 3:
            gen = random.choice([1, 2, 3])
            sign = random.choice([1, -1])
            word.append(sign * gen)
            used_gens.add(gen)
        name = '.'.join(f's{abs(w)}{"" if w>0 else "^-1"}' for w in word)
        if len(name) > 52:
            name = name[:49] + '...'
        run(name, word)

    print(f'\n=== {counts[0]} passed, {counts[1]} failed out of {counts[0]+counts[1]} ===')
