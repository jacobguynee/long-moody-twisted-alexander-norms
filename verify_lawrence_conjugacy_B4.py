#!/usr/bin/env python3
"""
verify_lawrence_conjugacy_B4.py

Verifies that the iterated reduced Long-Moody representation of B_4
matches Lawrence's representation by checking:
1. Both satisfy B_4 braid relations
2. Corresponding generators have the same eigenvalues (characteristic polynomials)
"""

from sympy import symbols, eye, zeros, Matrix, cancel, factor, Poly, Symbol

Q, q, a = symbols('Q q a')
lam = Symbol('lambda')


def build_Ni(i, gi, n, d):
    nb = n - 1
    Id = eye(d)
    Zd = zeros(d)
    blocks = [[Id if r == c else Zd for c in range(nb)] for r in range(nb)]
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
    n = len(rho_g)
    num_sigmas = len(rho_sigma)
    d = rho_g[0].shape[0]
    if partition is None:
        partition = [1] * n
    is_boundary = [partition[i] != partition[i+1] for i in range(num_sigmas)]
    nb = n - 1
    rho_plus = []
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
    if last_is_boundary:
        gn = rho_sigma[n-1]
    else:
        gn = rho_sigma[n-1] * rho_sigma[n-1]
    rho_g = [None] * n
    rho_g[n-1] = cancel(gn)
    for i in range(n-2, -1, -1):
        rho_g[i] = cancel(rho_sigma[i].inv() * rho_g[i+1] * rho_sigma[i])
    return rho_g


def lawrence_generators(n, i, q1, q2, alpha):
    d = n * (n - 1)
    M = zeros(d)
    pairs = []
    idx = {}
    count = 0
    for la in range(1, n+1):
        for mu in range(1, n+1):
            if la != mu:
                pairs.append((la, mu))
                idx[(la, mu)] = count
                count += 1
    def p(a, b): return idx[(a, b)]
    si = i + 1; si1 = si + 1

    r = p(si, si1)
    M[r, p(si1, si)] = 1/(q1*q2)
    for j in range(1, si):
        M[r, p(si1, j)] = (1/q2 - 1)/q1
    for k in range(si1+1, n+1):
        M[r, p(si1, k)] = (1/q2 - 1)/q1

    r = p(si1, si)
    M[r, p(si, si1)] = 1/(q1*q2*alpha)
    for j in range(1, si):
        M[r, p(j, si1)] = (1/q1 - 1)/(q2*alpha)
    for k in range(si1+1, n+1):
        M[r, p(k, si1)] = (1/q1 - 1)/q2

    for j in range(1, si):
        M[p(si,j), p(si1,j)] = 1/q1
        M[p(j,si), p(j,si1)] = 1/q2
        M[p(si1,j), p(si,j)] = 1; M[p(si1,j), p(si1,j)] = 1-1/q1
        M[p(j,si1), p(j,si)] = 1; M[p(j,si1), p(j,si1)] = 1-1/q2

    for k in range(si1+1, n+1):
        M[p(si,k), p(si1,k)] = 1/q1
        M[p(k,si), p(k,si1)] = 1/q2
        M[p(si1,k), p(si,k)] = 1; M[p(si1,k), p(si1,k)] = 1-1/q1
        M[p(k,si1), p(k,si)] = 1; M[p(k,si1), p(k,si1)] = 1-1/q2

    for s_idx in range(d):
        la, mu = pairs[s_idx]
        if la != si and la != si1 and mu != si and mu != si1:
            M[s_idx, s_idx] = 1
    return M


def check_zero(M):
    return all(cancel(M[i, j]) == 0 for i in range(M.rows) for j in range(M.cols))


def charpoly(M):
    """Characteristic polynomial, factored."""
    return factor(M.charpoly(lam).as_expr())


if __name__ == '__main__':
    # === Step 1 ===
    print('=== Step 1: First reduced Long-Moody iteration ===')
    rho_g_seed = [Matrix([[Q]]), Matrix([[Q]]), Matrix([[Q]]), Matrix([[Q]]), Matrix([[a]])]
    rho_s_seed = [Matrix([[1]]), Matrix([[1]]), Matrix([[1]]), Matrix([[1]])]
    L4Pre = reduced_lm_colored(rho_g_seed, rho_s_seed, [1,1,1,1,2])
    print(f'  Output: {len(L4Pre)} matrices, each {L4Pre[0].shape[0]}x{L4Pre[0].shape[1]}')

    # === Step 2 ===
    print('\n=== Step 2: Recover free group generators ===')
    L4PreFree = free_group_gens(L4Pre, [1,1,1,1,2])
    print(f'  Recovered g_1,...,g_4 as {L4PreFree[0].shape[0]}x{L4PreFree[0].shape[1]} matrices')

    # === Step 3 ===
    print('\n=== Step 3: Second reduced Long-Moody iteration ===')
    rho_g2 = [g * q for g in L4PreFree]
    L4 = reduced_lm_colored(rho_g2, L4Pre[0:3])
    s1, s2, s3 = L4[0], L4[1], L4[2]
    print(f'  Output: {len(L4)} matrices, each {s1.shape[0]}x{s1.shape[1]}')

    # === Step 4 ===
    print('\n=== Step 4: Compute Lawrence representation for n=4 ===')
    L_s1 = lawrence_generators(4, 0, 1/Q, 1/q, 1/a)
    L_s2 = lawrence_generators(4, 1, 1/Q, 1/q, 1/a)
    L_s3 = lawrence_generators(4, 2, 1/Q, 1/q, 1/a)
    print(f'  Lawrence matrices: each {L_s1.shape[0]}x{L_s1.shape[1]}')

    # === Step 5: Braid relations ===
    print('\n=== Step 5: Braid relations ===')
    print('  Our representation:')
    print(f'    s1*s2*s1 = s2*s1*s2: {check_zero(cancel(s1*s2*s1 - s2*s1*s2))}')
    print(f'    s2*s3*s2 = s3*s2*s3: {check_zero(cancel(s2*s3*s2 - s3*s2*s3))}')
    print(f'    s1*s3 = s3*s1:       {check_zero(cancel(s1*s3 - s3*s1))}')

    print('  Lawrence representation:')
    print(f'    s1*s2*s1 = s2*s1*s2: {check_zero(cancel(L_s1*L_s2*L_s1 - L_s2*L_s1*L_s2))}')
    print(f'    s2*s3*s2 = s3*s2*s3: {check_zero(cancel(L_s2*L_s3*L_s2 - L_s3*L_s2*L_s3))}')
    print(f'    s1*s3 = s3*s1:       {check_zero(cancel(L_s1*L_s3 - L_s3*L_s1))}')

    # === Step 6: Compare characteristic polynomials ===
    print('\n=== Step 6: Compare characteristic polynomials ===')
    for name, ours, law in [('sigma_1', s1, L_s1), ('sigma_2', s2, L_s2), ('sigma_3', s3, L_s3)]:
        cp_ours = charpoly(ours)
        cp_law  = charpoly(law)
        diff = factor(cp_ours - cp_law)
        match = (diff == 0)
        print(f'  {name}: {"MATCH" if match else "MISMATCH"}')
        if not match:
            print(f'    Ours:     {cp_ours}')
            print(f'    Lawrence: {cp_law}')
            print(f'    Diff:     {diff}')

    print('\n=== Done. ===')
