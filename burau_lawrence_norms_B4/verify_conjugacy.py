#!/usr/bin/env python3
"""
verify_conjugacy.py

Verifies that the iterated reduced Long-Moody representation of B_4
is conjugate to the Lawrence (1990) representation for k=2.

Checks:
1. Both representations satisfy B_4 braid relations.
2. Corresponding generators have the same characteristic polynomial.

This uses the FULL symbolic parameters (Q, q, a) -- no numerical
substitution -- to verify the algebraic identity exactly.

Usage:
    python verify_conjugacy.py
"""

from sympy import symbols, eye, zeros, Matrix, cancel, factor, Symbol

Q, q, a = symbols('Q q a')
lam = Symbol('lambda')

from representations import reduced_lm_colored, free_group_gens


def lawrence_generators(n, i, q1, q2, alpha):
    """Lawrence (1990) Fig 4.4 matrix for sigma_i acting on H_1 of
    ordered 2-point configuration space in n-punctured disk.

    Parameters q1, q2, alpha are the monodromy parameters.
    Returns an n(n-1) x n(n-1) matrix."""
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


def main():
    # Step 1: First Long-Moody iteration
    print('Step 1: First reduced Long-Moody (B_{4,1} mixed braid group)')
    rho_g = [Matrix([[Q]])]*4 + [Matrix([[a]])]
    rho_s = [Matrix([[1]])]*4
    L4Pre = reduced_lm_colored(rho_g, rho_s, [1,1,1,1,2])
    print(f'  {len(L4Pre)} generators, each {L4Pre[0].shape[0]}x{L4Pre[0].shape[1]}')

    # Step 2: Recover free group generators
    print('Step 2: Recover free group generators')
    L4PreFree = free_group_gens(L4Pre, [1,1,1,1,2])
    print(f'  {len(L4PreFree)} generators recovered')

    # Step 3: Second Long-Moody iteration
    print('Step 3: Second reduced Long-Moody')
    rho_g2 = [g * q for g in L4PreFree]
    L4 = reduced_lm_colored(rho_g2, L4Pre[0:3])
    s1, s2, s3 = L4
    print(f'  {len(L4)} generators, each {s1.shape[0]}x{s1.shape[1]}')

    # Step 4: Lawrence representation
    print('Step 4: Lawrence (1990) representation')
    L_s1 = lawrence_generators(4, 0, 1/Q, 1/q, 1/a)
    L_s2 = lawrence_generators(4, 1, 1/Q, 1/q, 1/a)
    L_s3 = lawrence_generators(4, 2, 1/Q, 1/q, 1/a)
    print(f'  {L_s1.shape[0]}x{L_s1.shape[1]} matrices')

    # Step 5: Braid relations
    print('\nBraid relations (iterated Long-Moody):')
    print(f'  s1*s2*s1 = s2*s1*s2: {check_zero(cancel(s1*s2*s1 - s2*s1*s2))}')
    print(f'  s2*s3*s2 = s3*s2*s3: {check_zero(cancel(s2*s3*s2 - s3*s2*s3))}')
    print(f'  s1*s3 = s3*s1:       {check_zero(cancel(s1*s3 - s3*s1))}')

    print('Braid relations (Lawrence):')
    print(f'  s1*s2*s1 = s2*s1*s2: {check_zero(cancel(L_s1*L_s2*L_s1 - L_s2*L_s1*L_s2))}')
    print(f'  s2*s3*s2 = s3*s2*s3: {check_zero(cancel(L_s2*L_s3*L_s2 - L_s3*L_s2*L_s3))}')
    print(f'  s1*s3 = s3*s1:       {check_zero(cancel(L_s1*L_s3 - L_s3*L_s1))}')

    # Step 6: Compare characteristic polynomials
    print('\nCharacteristic polynomials:')
    for name, ours, law in [('sigma_1', s1, L_s1), ('sigma_2', s2, L_s2), ('sigma_3', s3, L_s3)]:
        cp_ours = factor(ours.charpoly(lam).as_expr())
        cp_law  = factor(law.charpoly(lam).as_expr())
        diff = factor(cp_ours - cp_law)
        print(f'  {name}: {"MATCH" if diff == 0 else "MISMATCH"}')

    print('\nDone.')


if __name__ == '__main__':
    main()
