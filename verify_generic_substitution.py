#!/usr/bin/env python3
"""
verify_generic_substitution.py

Check that our Q=2, a=3 substitution doesn't accidentally cancel terms.
Test the same braids with multiple (Q, a) values and verify the Newton
polytopes match across all substitutions.
"""

from sympy import symbols, eye, zeros, Matrix, cancel, Poly, Rational, fraction
import sys

q, x = symbols('q x')

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

def build_lawrence_B4(Q_val, a_val):
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

if __name__ == '__main__':
    scale = Rational(1, 4)

    # Multiple (Q, a) pairs to test genericity
    qa_pairs = [
        (Rational(2), Rational(3)),
        (Rational(5), Rational(7)),
        (Rational(3), Rational(11)),
        (Rational(7), Rational(2)),
        (Rational(13), Rational(5)),
    ]

    # Test braids
    test_braids = [
        ('s1.s2.s1^-1.s3', [1, 2, -1, 3]),
        ('s1^3.s2^-1.s3^2', [1, 1, 1, -2, 3, 3]),
        ('[s1,s2].s3', [1, 2, -1, -2, 3]),
        ('s1^2.s2^-1.s3.s1^-1.s2.s3^-2', [1, 1, -2, 3, -1, 2, -3, -3]),
        ('s1.s2.s3.s1^-1.s2^-1.s3^-1', [1, 2, 3, -1, -2, -3]),
        ('s1^5.s3^-3', [1,1,1,1,1,-3,-3,-3]),
        ('s1.s2.s3.s2.s1^-1.s2^-1.s3^-1.s2^-1', [1,2,3,2,-1,-2,-3,-2]),
    ]

    print('Building Burau...'); sys.stdout.flush()
    bg = build_burau_B4(); bgi = [cancel(g.inv()) for g in bg]
    print('Done.\n')

    all_consistent = True

    for braid_name, word in test_braids:
        Mb = eval_word(word, bg, bgi)
        hull_b = convex_hull_2d(list(newton_polytope_xq(Mb)))

        hulls_l = []
        for Q_val, a_val in qa_pairs:
            print(f'  Building Lawrence (Q={Q_val}, a={a_val})...', end=''); sys.stdout.flush()
            lg = build_lawrence_B4(Q_val, a_val)
            lgi = [cancel(g.inv()) for g in lg]
            Ml = eval_word(word, lg, lgi)
            hull_l = convex_hull_2d(list(newton_polytope_xq(Ml)))
            hull_l_s = scale_polytope(hull_l, scale)
            hulls_l.append(hull_l_s)
            print(f' hull={hull_l_s}'); sys.stdout.flush()

        # Check all Lawrence hulls agree
        all_same = all(h == hulls_l[0] for h in hulls_l)
        burau_match = (hull_b == hulls_l[0])

        if all_same and burau_match:
            print(f'{braid_name}: ALL CONSISTENT, matches Burau  {hull_b}')
        elif all_same and not burau_match:
            print(f'{braid_name}: *** Lawrence hulls agree but DIFFER from Burau ***')
            print(f'  Burau:    {hull_b}')
            print(f'  Lawrence: {hulls_l[0]}')
            all_consistent = False
        else:
            print(f'{braid_name}: *** INCONSISTENT across (Q,a) values! ***')
            for i, (Q_val, a_val) in enumerate(qa_pairs):
                print(f'  (Q={Q_val}, a={a_val}): {hulls_l[i]}')
            all_consistent = False
        print(); sys.stdout.flush()

    if all_consistent:
        print('All braids consistent across all (Q,a) substitutions.')
    else:
        print('*** INCONSISTENCY FOUND ***')
