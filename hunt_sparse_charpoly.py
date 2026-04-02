#!/usr/bin/env python3
"""
hunt_sparse_charpoly.py

Look for B4 braids where the Burau characteristic polynomial has
few monomials (maximal cancellation). These are the most likely
candidates for Newton polytope disagreement, because cancellations
might not happen in the Lawrence representation.
"""

from sympy import symbols, eye, zeros, Matrix, cancel, Poly, Rational, fraction
import sys
import itertools

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

if __name__ == '__main__':
    print('Building representations...'); sys.stdout.flush()
    bg = build_burau_B4(); bgi = [cancel(g.inv()) for g in bg]
    lg = build_lawrence_B4(); lgi = [cancel(g.inv()) for g in lg]
    print('Done.\n'); sys.stdout.flush()

    scale = Rational(1, 4)
    counts = [0, 0]
    sparse_found = []

    # Enumerate all B4 braids of word length up to 7, using generators 1,2,3,-1,-2,-3
    # But that's 6^7 = 279936 braids. Too many.
    # Instead, enumerate up to length 5 (6^5 = 7776) and look for sparse charpolys.

    gens_list = [1, -1, 2, -2, 3, -3]

    print('Phase 1: Find braids with sparse Burau charpolys (word length 4-6)')
    sys.stdout.flush()

    for length in [4, 5, 6]:
        print(f'\n--- Word length {length} ---'); sys.stdout.flush()
        min_monoms = 999

        for word in itertools.product(gens_list, repeat=length):
            # Skip words that don't use all 3 generators
            used = set(abs(w) for w in word)
            if len(used) < 3: continue

            # Skip trivial cancellations (consecutive inverse pairs)
            skip = False
            for i in range(len(word)-1):
                if word[i] + word[i+1] == 0:
                    skip = True; break
            if skip: continue

            Mb = eval_word(list(word), bg, bgi)
            pts = newton_polytope_xq(Mb)
            n_monoms = len(pts)

            if n_monoms < min_monoms:
                min_monoms = n_monoms
                name = '.'.join(f's{abs(w)}{"" if w>0 else "^-1"}' for w in word)
                print(f'  New min: {n_monoms} monomials: {name}'); sys.stdout.flush()

            if n_monoms <= 5:  # Very sparse
                sparse_found.append((list(word), n_monoms, pts))

        print(f'  Min monomials at length {length}: {min_monoms}')
        sys.stdout.flush()

    # Now test the sparsest braids against Lawrence
    print(f'\n\nPhase 2: Test {len(sparse_found)} sparse braids against Lawrence')
    sys.stdout.flush()

    # Sort by monomial count, test sparsest first
    sparse_found.sort(key=lambda x: x[1])

    tested = set()
    for word, n_mon, pts_b in sparse_found[:50]:  # test up to 50
        word_key = tuple(word)
        if word_key in tested: continue
        tested.add(word_key)

        name = '.'.join(f's{abs(w)}{"" if w>0 else "^-1"}' for w in word)
        hull_b = convex_hull_2d(list(pts_b))

        Ml = eval_word(word, lg, lgi)
        pts_l = newton_polytope_xq(Ml)
        hull_l = convex_hull_2d(list(pts_l))
        hull_l_s = scale_polytope(hull_l, scale)

        match = (hull_b == hull_l_s)
        if match:
            print(f'  [{n_mon} mon] {name:<40} YES')
            counts[0] += 1
        else:
            print(f'  [{n_mon} mon] {name:<40} *** NO ***')
            print(f'    Burau monoms:   {sorted(pts_b)}')
            print(f'    NP(Burau):      {hull_b}')
            print(f'    NP(Lawrence)/4: {hull_l_s}')
            counts[1] += 1
        sys.stdout.flush()

    print(f'\n=== {counts[0]} passed, {counts[1]} failed out of {counts[0]+counts[1]} ===')
