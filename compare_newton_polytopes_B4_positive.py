#!/usr/bin/env python3
"""
compare_newton_polytopes_B4_positive.py

Positive 4-braids of word length 10-30, using all 3 generators.
Compare NP(Lawrence k=2)/4 vs NP(Burau) in (x, q).

FAST METHOD: Instead of symbolic charpoly over Q(q), substitute many
numerical values for q (rational primes) and reconstruct the Newton
polytope from the numerical charpolys. Each charpoly at fixed q is just
a numerical eigenvalue problem — instant.

Newton polytope detection: for each q-value, compute charpoly coefficients.
The coefficient of x^k is a Laurent polynomial in q. We evaluate at enough
q-values to determine which monomials q^j have nonzero coefficients, giving
the support. Then take convex hull.
"""

import numpy as np
from numpy.polynomial import polynomial as P
from fractions import Fraction
import sys


# ============================================================
# Build matrices numerically for a given q value
# ============================================================

def build_Ni_num(i, gi, n, d):
    nb = n - 1
    Id = np.eye(d); Zd = np.zeros((d, d))
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

def left_mult_num(blocks, M, nb):
    return [[M @ blocks[r][c] for c in range(nb)] for r in range(nb)]

def mult_blocks_num(A, B, nb, d):
    C = [[np.zeros((d,d)) for _ in range(nb)] for _ in range(nb)]
    for r in range(nb):
        for c in range(nb):
            acc = np.zeros((d,d))
            for k in range(nb): acc += A[r][k] @ B[k][c]
            C[r][c] = acc
    return C

def assemble_num(blocks, nb, d):
    N = nb * d; M = np.zeros((N, N))
    for r in range(nb):
        for c in range(nb): M[r*d:(r+1)*d, c*d:(c+1)*d] = blocks[r][c]
    return M

def reduced_lm_colored_num(rho_g, rho_sigma, partition=None):
    n = len(rho_g); ns = len(rho_sigma); d = rho_g[0].shape[0]
    if partition is None: partition = [1]*n
    is_bdy = [partition[i] != partition[i+1] for i in range(ns)]
    nb = n-1; out = []
    for i in range(ns):
        if not is_bdy[i]:
            blocks = build_Ni_num(i, rho_g[i], n, d)
            blocks = left_mult_num(blocks, rho_sigma[i], nb)
            out.append(assemble_num(blocks, nb, d))
        else:
            gi_conj = rho_g[i] @ rho_g[i+1] @ np.linalg.inv(rho_g[i])
            bN = build_Ni_num(i, rho_g[i], n, d)
            bNp = build_Ni_num(i, gi_conj, n, d)
            bp = mult_blocks_num(bNp, bN, nb, d)
            si2 = rho_sigma[i] @ rho_sigma[i]
            bp = left_mult_num(bp, si2, nb)
            out.append(assemble_num(bp, nb, d))
    return out

def free_group_gens_num(rho_sigma, partition=None):
    n = len(rho_sigma)
    if partition is None: partition = [1]*(n+1)
    last_bdy = partition[n-1] != partition[n]
    gn = rho_sigma[n-1] if last_bdy else rho_sigma[n-1] @ rho_sigma[n-1]
    rg = [None]*n; rg[n-1] = gn
    for i in range(n-2,-1,-1):
        si_inv = np.linalg.inv(rho_sigma[i])
        rg[i] = si_inv @ rg[i+1] @ rho_sigma[i]
    return rg


def build_burau_B4_num(qv):
    rg = [np.array([[qv]])]*4; rs = [np.array([[1.0]])]*3
    return reduced_lm_colored_num(rg, rs)

def build_lawrence_B4_num(qv, Qv=2.0, av=3.0):
    rg_s = [np.array([[Qv]])]*4 + [np.array([[av]])]
    rs_s = [np.array([[1.0]])]*4
    L4Pre = reduced_lm_colored_num(rg_s, rs_s, [1,1,1,1,2])
    L4Free = free_group_gens_num(L4Pre, [1,1,1,1,2])
    rg2 = [g*qv for g in L4Free]
    return reduced_lm_colored_num(rg2, L4Pre[0:3])


def eval_word_num(word, gens):
    M = np.eye(gens[0].shape[0])
    for w in word: M = M @ gens[w-1]
    return M


# ============================================================
# Newton polytope via numerical sampling
# ============================================================

def charpoly_coeffs(M):
    """Return coefficients of charpoly det(xI - M) as [c_0, c_1, ..., c_n]
    where charpoly = c_0 + c_1*x + ... + c_n*x^n."""
    return np.array(np.poly(M)[::-1])  # np.poly gives [c_n, ..., c_0]


def detect_newton_polytope(word, b_gens_func, l_gens_func, n_burau, n_law):
    """Detect Newton polytope by evaluating charpoly at multiple q values.

    For each x-degree k, the coefficient c_k(q) is a Laurent polynomial in q.
    We evaluate at many q values, fit via Vandermonde, and find the support
    (nonzero q-exponents).

    Returns (hull_burau, hull_lawrence) as sorted lists of (x_deg, q_deg) vertices.
    """
    # Use enough q-samples to determine Laurent polynomial coefficients
    # Max possible q-degree range: for word length L, q-exponents are in [-L*n, L*n] roughly
    # We'll use 200 sample points which is plenty

    # Sample at q = 2^(k/20) for varied magnitudes to avoid conditioning issues
    # Actually, use rational-ish points that are well-separated
    q_samples = [1.1 + 0.3*i for i in range(80)]

    # For each sample, compute charpoly coefficients
    burau_data = []  # list of (qv, coeffs_array)
    law_data = []

    for qv in q_samples:
        bg = b_gens_func(qv)
        lg = l_gens_func(qv)
        Mb = eval_word_num(word, bg)
        Ml = eval_word_num(word, lg)
        burau_data.append((qv, charpoly_coeffs(Mb)))
        law_data.append((qv, charpoly_coeffs(Ml)))

    hull_b = extract_newton_polytope(burau_data, n_burau)
    hull_l = extract_newton_polytope(law_data, n_law)

    return hull_b, hull_l


def extract_newton_polytope(data, mat_size):
    """Given [(qv, coeffs)] data, determine the Newton polytope.

    For each x-degree k (0..mat_size), coefficient c_k is c_k(q) = sum_j a_{k,j} q^j.
    We know c_{mat_size} = 1 (monic), c_0 = (-1)^n det(M).

    Strategy: for each k, evaluate c_k at many q values. The function c_k(q)
    is a Laurent polynomial. Find its q-support by fitting.

    We detect the min and max q-exponents by looking at log|c_k(q)| behavior.
    More precisely: divide c_k(q) by q^j for various j and check if the result
    is polynomial (bounded as q->inf). The Newton polytope is the convex hull
    of all (k, j) with nonzero a_{k,j}.
    """
    points = set()

    for k in range(mat_size + 1):
        # Extract c_k values at each q
        vals = [(qv, coeffs[k]) for qv, coeffs in data]

        # Find the q-exponent range by analyzing scaling
        # c_k(q) = sum_j a_j * q^j
        # Try to find the support by least-squares fitting

        # First, estimate the range of j by looking at log-log behavior
        # For large q, c_k(q) ~ a_{j_max} * q^{j_max}
        q_arr = np.array([v[0] for v in vals])
        c_arr = np.array([v[1] for v in vals])

        if np.all(np.abs(c_arr) < 1e-10):
            continue  # zero coefficient

        # Estimate max exponent from large-q behavior
        # log|c_k| ~ j_max * log(q) + const for large q
        mask = np.abs(c_arr) > 1e-15
        if np.sum(mask) < 10:
            # Very sparse — just mark the point with estimated exponent
            points.add((k, 0))
            continue

        log_q = np.log(q_arr[mask])
        log_c = np.log(np.abs(c_arr[mask]))

        # Fit line to estimate dominant exponent
        j_est = np.polyfit(log_q, log_c, 1)[0]
        j_max_est = int(np.round(j_est))

        # Now try fitting c_k(q) / q^{j_shift} as a polynomial
        # Search around the estimate
        search_range = max(50, abs(j_max_est) + 20)
        j_min_found = None
        j_max_found = None

        # Fit a Laurent polynomial: try range of j from j_max_est - search to j_max_est + 10
        # Use Vandermonde approach: c_k(q) = sum_{j=j_lo}^{j_hi} a_j q^j
        # => c_k(q) / q^{j_lo} = sum_{j=0}^{j_hi-j_lo} a_{j+j_lo} q^j
        # This is a polynomial fit

        j_lo = j_max_est - search_range
        j_hi = j_max_est + 10
        n_terms = j_hi - j_lo + 1

        # Normalize: y = c_k(q) / q^{j_lo}
        y = c_arr / (q_arr ** j_lo)

        # Vandermonde: fit polynomial of degree n_terms-1
        # V[i, j] = q_samples[i]^j
        if n_terms > len(q_arr):
            n_terms = len(q_arr)
            j_hi = j_lo + n_terms - 1

        V = np.vander(q_arr, N=n_terms, increasing=True)
        # Solve least squares
        coeffs_fit, _, _, _ = np.linalg.lstsq(V, y, rcond=None)

        # Find nonzero coefficients (above noise threshold)
        threshold = np.max(np.abs(coeffs_fit)) * 1e-8
        for idx_j, c in enumerate(coeffs_fit):
            if abs(c) > threshold:
                actual_j = j_lo + idx_j
                points.add((k, actual_j))

    return convex_hull_2d(list(points))


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
    return sorted(set((v[0]*s, v[1]*s) for v in vertices))


if __name__ == '__main__':
    print('Positive 4-braids, length 10-30. Fast numerical method.')
    print('Q=2, a=3. Scaling Lawrence by 1/4.\n')
    sys.stdout.flush()

    braids = [
        # Length 10
        ('(123)^3.1',            [1,2,3]*3+[1]),
        ('(321)^3.3',            [3,2,1]*3+[3]),
        ('1.3.2.1.3.2.1.3.2.1', [1,3,2,1,3,2,1,3,2,1]),
        ('1.2.1.2.3.2.3.1.2.3', [1,2,1,2,3,2,3,1,2,3]),
        # Length 12
        ('(123)^4',              [1,2,3]*4),
        ('(1213)^3',             [1,2,1,3]*3),
        ('(2132)^3',             [2,1,3,2]*3),
        ('(1323)^3',             [1,3,2,3]*3),
        # Length 15
        ('(12321)^3',            [1,2,3,2,1]*3),
        ('(13213)^3',            [1,3,2,1,3]*3),
        ('(21321)^3',            [2,1,3,2,1]*3),
        # Length 18
        ('(123)^6',              [1,2,3]*6),
        ('(123212)^3',           [1,2,3,2,1,2]*3),
        ('(132132)^3',           [1,3,2,1,3,2]*3),
        ('(213231)^3',           [2,1,3,2,3,1]*3),
        # Length 20
        ('(12132)^4',            [1,2,1,3,2]*4),
        ('(12321)^4',            [1,2,3,2,1]*4),
        ('(31213)^4',            [3,1,2,1,3]*4),
        ('(23123)^4',            [2,3,1,2,3]*4),
        # Length 24
        ('(123)^8',              [1,2,3]*8),
        ('(1213)^6',             [1,2,1,3]*6),
        ('(123213)^4',           [1,2,3,2,1,3]*4),
        ('(132312)^4',           [1,3,2,3,1,2]*4),
        ('(213123)^4',           [2,1,3,1,2,3]*4),
        ('(312132)^4',           [3,1,2,1,3,2]*4),
        # Length 25
        ('(12321)^5',            [1,2,3,2,1]*5),
        ('(31231)^5',            [3,1,2,3,1]*5),
        ('(21312)^5',            [2,1,3,1,2]*5),
        # Length 27
        ('(123)^9',              [1,2,3]*9),
        ('(132132132)^3',        [1,3,2,1,3,2,1,3,2]*3),
        # Length 28
        ('(1213231)^4',          [1,2,1,3,2,3,1]*4),
        ('(3212132)^4',          [3,2,1,2,1,3,2]*4),
        # Length 30
        ('(123)^10',             [1,2,3]*10),
        ('(12321)^6',            [1,2,3,2,1]*6),
        ('(123213)^5',           [1,2,3,2,1,3]*5),
        ('(213123)^5',           [2,1,3,1,2,3]*5),
        ('(132312)^5',           [1,3,2,3,1,2]*5),
        ('(1213)^7.12',          [1,2,1,3]*7+[1,2]),
    ]

    pass_count = 0; fail_count = 0
    scale = 0.25  # 1/4

    for name, word in braids:
        wlen = len(word)
        print(f'[{wlen:>2}] {name:<30}', end=''); sys.stdout.flush()

        hull_b, hull_l = detect_newton_polytope(
            word, build_burau_B4_num, build_lawrence_B4_num, 3, 12)

        hull_l_scaled = scale_polytope(hull_l, scale)
        # Convert to comparable format (round to handle floating point)
        hull_b_r = sorted(set((round(a), round(b)) for a,b in hull_b))
        hull_l_r = sorted(set((round(a*4)/4, round(b*4)/4) for a,b in hull_l_scaled))
        # Check: after scaling by 1/4, all coords should be k/4 integers
        # Burau coords are integers
        match = (hull_b_r == [(int(a), int(b)) for a,b in hull_l_r])

        if match:
            print(f'YES  {hull_b_r}'); pass_count += 1
        else:
            print(f'NO')
            print(f'  B:  {hull_b_r}')
            print(f'  L/4:{hull_l_r}')
            fail_count += 1
        sys.stdout.flush()

    print(f'\n=== {pass_count} passed, {fail_count} failed out of {pass_count+fail_count} ===')
