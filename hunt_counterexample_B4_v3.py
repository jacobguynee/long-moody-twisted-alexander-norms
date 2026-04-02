#!/usr/bin/env python3
"""
hunt_counterexample_B4_v3.py

"Almost reducible" pseudo-Anosov braids in B_4.

Strategy: Take a reducible braid and perturb it by a small element
that connects the two reducing pieces. This produces pseudo-Anosov
braids with small dilatation that are "close to reducible" — exactly
where Burau and Lawrence might disagree.

Key construction:
  R = s1^a * s3^b  (reducible: strands 1-2 decouple from 3-4)
  P = s2^c          (connects the two pieces)
  beta = R * P      (almost reducible pseudo-Anosov)

Also try:
  beta = s1^a * s2^c * s3^b  (interleaved)
  beta = (s1^a * s2)(s3^b * s2^-1)  etc.
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

def test_braid(name, word, b, bi, l, li, scale=Rational(1,4)):
    M_b = eval_word(word, b, bi); M_l = eval_word(word, l, li)
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
        if test_braid(name, word, b, bi, l, li, scale):
            counts[0] += 1
        else:
            counts[1] += 1

    # === FAMILY 1: s1^a . s2 . s3^b — "almost reducible" ===
    # The s2 connects the otherwise decoupled s1 and s3 parts.
    # Large a,b with small connecting s2 = almost reducible.
    print('=== s1^a . s2 . s3^b (single s2 bridge) ===')
    sys.stdout.flush()
    for aa, bb in [(5,5), (5,-5), (-5,5), (10,1), (1,10),
                   (10,-1), (-1,10), (7,-3), (-3,7), (8,-8),
                   (10,10), (10,-10), (-10,10), (15,1), (1,15)]:
        word = [1]*abs(aa) if aa>0 else [-1]*abs(aa)
        word += [2]
        word += [3]*abs(bb) if bb>0 else [-3]*abs(bb)
        run(f's1^{aa}.s2.s3^{bb}', word)

    # === FAMILY 2: s1^a . s2^c . s3^b with various small c ===
    print('\n=== s1^a . s2^c . s3^b (small bridge) ===')
    sys.stdout.flush()
    for aa, cc, bb in [(5,2,5), (5,-1,5), (5,1,-5), (5,-2,5),
                       (8,1,2), (2,1,8), (8,-1,2), (2,-1,8),
                       (10,1,10), (10,-1,10), (10,2,10), (10,-2,10),
                       (5,1,5), (5,-1,-5), (-5,1,5), (-5,-1,-5)]:
        word = [1]*abs(aa) if aa>0 else [-1]*abs(aa)
        word += [2]*abs(cc) if cc>0 else [-2]*abs(cc)
        word += [3]*abs(bb) if bb>0 else [-3]*abs(bb)
        run(f's1^{aa}.s2^{cc}.s3^{bb}', word)

    # === FAMILY 3: (s1^a . s3^b) . s2 . (s1^c . s3^d) ===
    # Two reducible blocks connected by s2
    print('\n=== (s1^a.s3^b).s2.(s1^c.s3^d) ===')
    sys.stdout.flush()
    for (aa,bb,cc,dd) in [(3,2,2,3), (3,-2,2,-3), (5,1,-1,5),
                          (-3,2,2,-3), (4,4,-4,-4), (5,3,-3,-5),
                          (2,5,5,2), (1,1,-1,-1)]:
        word = ([1]*abs(aa) if aa>0 else [-1]*abs(aa))
        word += ([3]*abs(bb) if bb>0 else [-3]*abs(bb))
        word += [2]
        word += ([1]*abs(cc) if cc>0 else [-1]*abs(cc))
        word += ([3]*abs(dd) if dd>0 else [-3]*abs(dd))
        run(f's1^{aa}.s3^{bb}.s2.s1^{cc}.s3^{dd}', word)

    # === FAMILY 4: (s1^a . s2^{+-1})^k . s3^b ===
    # Pseudo-Anosov on strands 1-3, with s3 perturbation
    print('\n=== (s1^a . s2^e)^k . s3^b ===')
    sys.stdout.flush()
    for aa, e, k, bb in [(2,1,3,1), (2,-1,3,1), (3,1,2,2), (3,-1,2,2),
                         (2,1,3,-1), (2,-1,3,-1), (2,1,5,1), (2,-1,5,1),
                         (4,1,2,1), (4,-1,2,1), (3,1,3,3), (3,-1,3,3),
                         (2,1,4,-2), (2,-1,4,2)]:
        word = ([1]*abs(aa) + ([2] if e>0 else [-2])) * k
        word += [3]*abs(bb) if bb>0 else [-3]*abs(bb)
        run(f'(s1^{aa}.s2^{e})^{k}.s3^{bb}', word)

    # === FAMILY 5: s1^a . s2^{-1} . s3^b . s2 (almost reducible loop) ===
    print('\n=== s1^a . s2^-1 . s3^b . s2 ===')
    sys.stdout.flush()
    for aa, bb in [(3,3), (5,5), (5,-5), (-5,5), (8,2), (2,8),
                   (10,10), (7,-7), (3,-3), (-3,3)]:
        word = [1]*abs(aa) if aa>0 else [-1]*abs(aa)
        word += [-2]
        word += [3]*abs(bb) if bb>0 else [-3]*abs(bb)
        word += [2]
        run(f's1^{aa}.s2^-1.s3^{bb}.s2', word)

    print(f'\n=== {counts[0]} passed, {counts[1]} failed out of {counts[0]+counts[1]} ===')
