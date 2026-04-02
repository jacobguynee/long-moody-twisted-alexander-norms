#!/usr/bin/env python3
"""
hunt_exterior_cancellation.py

Extends v4 approach: binary exponentiation on precomputed f,g matrices.
Part A: 20 NEW diverse (f,g) families not tested in v4, with f^a g^b.
Part B: General words f^{a1}g^{b1}f^{a2}g^{b2}... using same families.
Part C: Random braid words built letter-by-letter from generators.
"""

from sympy import symbols, eye, zeros, Matrix, cancel, Poly, Rational, fraction
import sys
import random

random.seed(42)
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
def eval_word(word, gens, inv_gens):
    M = eye(gens[0].shape[0])
    for w in word:
        if w > 0: M = M * gens[w-1]
        else: M = M * inv_gens[-w-1]
    return cancel(M)

if __name__ == '__main__':
    print('Building representations...'); sys.stdout.flush()
    bg = build_burau_B4(); bgi = [cancel(g.inv()) for g in bg]
    lg = build_lawrence_B4(); lgi = [cancel(g.inv()) for g in lg]
    print('Done.\n'); sys.stdout.flush()

    scale = Rational(1, 4)
    counts = [0, 0]

    def run(label, Mb, Ml):
        pts_b = newton_polytope_xq(Mb); hull_b = convex_hull_2d(list(pts_b))
        pts_l = newton_polytope_xq(Ml); hull_l = convex_hull_2d(list(pts_l))
        hull_l_s = scale_polytope(hull_l, scale)
        match = (hull_b == hull_l_s)
        if match:
            print(f'  {label:<50} YES  {hull_b}')
            counts[0] += 1
        else:
            print(f'  {label:<50} *** NO ***')
            print(f'    NP(Burau):      {hull_b}')
            print(f'    NP(Lawrence)/4: {hull_l_s}')
            counts[1] += 1
        sys.stdout.flush()

    # ===== Part A: 20 NEW families not in v4, f^a g^b =====
    # v4 tested: s1.s2/s2^-1.s3, s1.s2/s3, s2.s1.s2^-1/s3,
    #   s2.s1.s2^-1/s2.s3.s2^-1, s1/s2.s3.s2^-1, s1^2/s2.s3.s2^-1,
    #   s1/s2.s3^2.s2^-1, s1.s3/s2, s1.s3^-1/s2, s1.s3/s2^-1,
    #   s1.s2.s1/s3, s1.s2.s1/s3^-1, [s1,s2]/s3, [s1,s2]/s3^-1,
    #   s1/[s2,s3], s1^-1/[s2,s3]
    # These are NEW:
    families_A = [
        ('f=s1, g=s3',                  [1],        [3]),
        ('f=s1, g=s3^-1',               [1],        [-3]),
        ('f=s1^-1, g=s3',               [-1],       [3]),
        ('f=s2, g=s1.s3',               [2],        [1, 3]),
        ('f=s2, g=s1.s3^-1',            [2],        [1, -3]),
        ('f=s2^-1, g=s1.s3',            [-2],       [1, 3]),
        ('f=s1.s2^-1, g=s3',            [1, -2],    [3]),
        ('f=s1.s2^-1, g=s2.s3',         [1, -2],    [2, 3]),
        ('f=s1.s2^-1, g=s3^-1',         [1, -2],    [-3]),
        ('f=s3.s2^-1, g=s1',            [3, -2],    [1]),
        ('f=s3.s2^-1, g=s1^-1',         [3, -2],    [-1]),
        ('f=s1.s3, g=s2.s1^-1',         [1, 3],     [2, -1]),
        ('f=s1.s3^-1, g=s2^-1',         [1, -3],    [-2]),
        ('f=s1.s2.s3, g=s2^-1',         [1, 2, 3],  [-2]),
        ('f=s1.s2.s3, g=s1^-1',         [1, 2, 3],  [-1]),
        ('f=s2.s3, g=s1^-1',            [2, 3],     [-1]),
        ('f=s2.s3, g=s1',               [2, 3],     [1]),
        ('f=s1.s2^-1.s3, g=s2',         [1, -2, 3], [2]),
        ('f=s1.s2^-1.s3, g=s2^-1',      [1, -2, 3], [-2]),
        ('f=s2.s1.s3, g=s2^-1.s1^-1.s3^-1', [2,1,3], [-2,-1,-3]),
    ]

    exps_A = [
        (5, -5), (5, 5), (10, -3), (3, -10),
        (8, -8), (10, 10), (15, -5), (20, -1),
    ]

    print('=== Part A: 20 new f^a g^b families ===\n')
    for fam_name, f_word, g_word in families_A:
        print(f'\n--- {fam_name} ---'); sys.stdout.flush()
        fb = eval_word(f_word, bg, bgi); gb = eval_word(g_word, bg, bgi)
        fl = eval_word(f_word, lg, lgi); gl = eval_word(g_word, lg, lgi)
        for aa, bb in exps_A:
            Mb = cancel(mat_power(fb, aa) * mat_power(gb, bb))
            Ml = cancel(mat_power(fl, aa) * mat_power(gl, bb))
            run(f'f^{aa} g^{bb}', Mb, Ml)

    print(f'\nPart A: {counts[0]}p {counts[1]}f\n')
    ca, cf = counts[0], counts[1]

    # ===== Part B: General words f^{a1}g^{b1}f^{a2}g^{b2}... =====
    print('=== Part B: General alternating words ===\n')

    # Use 8 families (mix from v4 and new), 4 general words each
    families_B = [
        ('s1/s3',         [1],     [3]),
        ('s1.s2/s3',      [1, 2],  [3]),
        ('s1/s2.s3.s2i',  [1],     [2, 3, -2]),
        ('s2/s1.s3',      [2],     [1, 3]),
        ('s1.s2i/s3',     [1, -2], [3]),
        ('s1.s3/s2',      [1, 3],  [2]),
        ('s1.s3i/s2',     [1, -3], [2]),
        ('s2.s3/s1i',     [2, 3],  [-1]),
    ]

    for fam_name, f_word, g_word in families_B:
        fb = eval_word(f_word, bg, bgi); gb = eval_word(g_word, bg, bgi)
        fl = eval_word(f_word, lg, lgi); gl = eval_word(g_word, lg, lgi)
        for trial in range(4):
            k = random.randint(2, 4)
            segs = [(random.choice([-1,1])*random.randint(2,5),
                     random.choice([-1,1])*random.randint(2,5)) for _ in range(k)]
            Mb = eye(3); Ml = eye(12)
            for ai, bi in segs:
                Mb = cancel(Mb * mat_power(fb, ai) * mat_power(gb, bi))
                Ml = cancel(Ml * mat_power(fl, ai) * mat_power(gl, bi))
            seg_str = '.'.join(f'f^{a}g^{b}' for a, b in segs)
            run(f'{fam_name}: {seg_str}', Mb, Ml)

    print(f'\nPart B: {counts[0]-ca}p {counts[1]-cf}f\n')
    cb, cbf = counts[0], counts[1]

    # ===== Part C: Random braid words =====
    print('=== Part C: Random braid words ===\n')

    for trial in range(15):
        length = random.randint(8, 18)
        word = [random.choice([-1,1])*random.randint(1,3) for _ in range(length)]
        used = set(abs(w) for w in word)
        if len(used) < 3: continue
        Mb = eval_word(word, bg, bgi)
        Ml = eval_word(word, lg, lgi)
        wstr = ''.join(f's{abs(w)}{"" if w>0 else "i"}' for w in word[:6]) + f'...(L={length})'
        run(wstr, Mb, Ml)

    print(f'\nPart C: {counts[0]-cb}p {counts[1]-cbf}f')
    print(f'\n=== TOTAL: {counts[0]} passed, {counts[1]} failed ===')
