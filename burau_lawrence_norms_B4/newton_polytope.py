"""
newton_polytope.py

Utilities for computing Newton polytopes of characteristic polynomials
and comparing them between representations.
"""

from sympy import cancel, Poly, Rational, fraction, eye
from representations import q, x


def newton_polytope_xq(M):
    """Newton polytope of det(xI - M) in (x, q) coordinates.

    Returns the set of (x-degree, q-degree) pairs with nonzero coefficients.
    Handles the case where charpoly entries are Laurent polynomials in q
    (negative q-powers) by extracting numerator and shifting.
    """
    cp_expr = cancel(M.charpoly(x).as_expr())
    num, den = fraction(cp_expr)
    p = Poly(num, x, q)
    monoms_num = p.monoms()
    q_shift = Poly(den, q).degree() if den.has(q) else 0
    return set((xd, qd - q_shift) for (xd, qd) in monoms_num)


def convex_hull_2d(points):
    """Convex hull of 2D integer points via Andrew's monotone chain."""
    points = sorted(set(points))
    if len(points) <= 2:
        return points
    def cross(O, A, B):
        return (A[0]-O[0])*(B[1]-O[1]) - (A[1]-O[1])*(B[0]-O[0])
    lower = []
    for p in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return sorted(set(lower[:-1] + upper[:-1]))


def scale_polytope(vertices, s):
    """Scale all coordinates by rational number s."""
    return sorted(set((Rational(v[0]) * s, Rational(v[1]) * s) for v in vertices))


def compare_newton_polytopes(M_burau, M_lawrence, scale):
    """Compare NP(Burau charpoly) with NP(Lawrence charpoly) * scale.

    Returns (match, hull_burau, hull_lawrence_scaled).
    """
    pts_b = newton_polytope_xq(M_burau)
    hull_b = convex_hull_2d(list(pts_b))
    pts_l = newton_polytope_xq(M_lawrence)
    hull_l = convex_hull_2d(list(pts_l))
    hull_l_s = scale_polytope(hull_l, scale)
    return (hull_b == hull_l_s), hull_b, hull_l_s


def mat_power(M, n):
    """Binary exponentiation for symbolic matrices."""
    if n == 0:
        return eye(M.shape[0])
    if n < 0:
        M = cancel(M.inv()); n = -n
    result = eye(M.shape[0]); base = M
    while n > 0:
        if n % 2 == 1:
            result = cancel(result * base)
        base = cancel(base * base)
        n //= 2
    return result


def eval_word(word, gens, inv_gens):
    """Evaluate a braid word on generators.

    word: list of signed ints. +i means sigma_i, -i means sigma_i^{-1}.
    gens: [sigma_1_matrix, sigma_2_matrix, ...]
    inv_gens: [sigma_1_inv, sigma_2_inv, ...]
    """
    M = eye(gens[0].shape[0])
    for w in word:
        if w > 0:
            M = M * gens[w - 1]
        else:
            M = M * inv_gens[-w - 1]
    return cancel(M)
