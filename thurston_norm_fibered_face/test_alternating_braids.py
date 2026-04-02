#!/usr/bin/env python3
"""
test_alternating_braids.py

Verify the fibered face computation on alternating braids in B_3 and B_4.

An alternating braid uses all generators with alternating signs:
  B_3: words in {s1, s2^{-1}}
  B_4: words in {s1, s2^{-1}, s3}

When all generators appear, the braid is pseudo-Anosov, and the author's
theorem gives: Alexander norm = Thurston norm for these braids.

This script:
1. Computes the reduced Burau matrix for each braid.
2. Computes det(I - Burau(beta)) = Alexander polynomial (1-variable).
3. Computes the Alexander norm (breadth of Newton polytope).
4. Computes the 2D Newton polytope of the characteristic polynomial.
5. Verifies internal consistency and prints the fibered face data.

Usage:
    python test_alternating_braids.py
"""

from sympy import Rational, cancel, Poly, symbols, fraction

from alexander_thurston import (
    build_burau, eval_word, alexander_polynomial_1var,
    newton_polytope_1d, alexander_norm_1d,
    newton_polytope_2d, convex_hull_2d,
    fibered_face_from_alexander, print_fibered_face,
)

q = symbols('q')


def test_braid(name, word, n):
    """Test a single braid and return results."""
    data = fibered_face_from_alexander(word, n)
    alex_norm = data['alexander_norm']
    genus = data['genus']
    hull = data['charpoly_newton_2d']

    # Basic consistency checks
    assert alex_norm >= 0, f"Negative Alexander norm for {name}"

    # For pseudo-Anosov braids, alexander norm should be positive
    if len(set(abs(w) for w in word)) == n - 1:
        assert alex_norm > 0, f"Zero Alexander norm for pseudo-Anosov {name}"

    # Genus should be half the Alexander norm (knot case)
    assert genus == Rational(alex_norm, 2), f"Genus mismatch for {name}"

    # Newton polytope 1D endpoints should be symmetric about 0
    # (after normalization, which we don't enforce here)
    np1d = data['alexander_poly_newton_1d']
    breadth = np1d[1] - np1d[0]
    assert breadth == alex_norm, f"Breadth != alex_norm for {name}"

    return data


def main():
    pass_count = 0
    total = 0

    # =========================================================
    # B_3 alternating braids: words in {s1, s2^{-1}}
    # =========================================================
    print("=" * 65)
    print("B_3 alternating braids (words in s1, s2^-1)")
    print("=" * 65)

    b3_braids = [
        ("s1.s2^-1",                   [1, -2]),
        ("s2^-1.s1",                   [-2, 1]),
        ("(s1.s2^-1)^2",              [1, -2, 1, -2]),
        ("s1^2.s2^-1",                [1, 1, -2]),
        ("s1.s2^-2",                  [1, -2, -2]),
        ("s1^2.s2^-2",                [1, 1, -2, -2]),
        ("(s1.s2^-1)^3",              [1, -2, 1, -2, 1, -2]),
        ("s1^3.s2^-1",                [1, 1, 1, -2]),
        ("s1.s2^-3",                  [1, -2, -2, -2]),
        ("s1^2.s2^-1.s1.s2^-1",      [1, 1, -2, 1, -2]),
        ("s1.s2^-1.s1^2.s2^-2",      [1, -2, 1, 1, -2, -2]),
        ("(s1.s2^-1)^4",              [1, -2, 1, -2, 1, -2, 1, -2]),
        ("s1^3.s2^-3",                [1, 1, 1, -2, -2, -2]),
        ("s1^4.s2^-1",                [1, 1, 1, 1, -2]),
        ("s1.s2^-4",                  [1, -2, -2, -2, -2]),
    ]

    for name, word in b3_braids:
        total += 1
        try:
            data = test_braid(name, word, 3)
            an = data['alexander_norm']
            g = data['genus']
            np1d = data['alexander_poly_newton_1d']
            print(f"  {name:<40} norm={str(an):<4} genus={str(g):<6} NP=[{np1d[0]},{np1d[1]}]  OK")
            pass_count += 1
        except Exception as e:
            print(f"  {name:<40} FAIL: {e}")

    # =========================================================
    # B_4 alternating braids: words in {s1, s2^{-1}, s3}
    # =========================================================
    print()
    print("=" * 65)
    print("B_4 alternating braids (words in s1, s2^-1, s3)")
    print("=" * 65)

    b4_braids = [
        ("s1.s2^-1.s3",                [1, -2, 3]),
        ("s3.s2^-1.s1",                [3, -2, 1]),
        ("s1.s2^-1.s3.s2^-1",          [1, -2, 3, -2]),
        ("s1^2.s2^-1.s3",              [1, 1, -2, 3]),
        ("s1.s2^-2.s3",                [1, -2, -2, 3]),
        ("s1.s2^-1.s3^2",              [1, -2, 3, 3]),
        ("s1^2.s2^-2.s3^2",            [1, 1, -2, -2, 3, 3]),
        ("(s1.s2^-1.s3)^2",            [1, -2, 3, 1, -2, 3]),
        ("s1^2.s2^-1.s3.s2^-1.s1",     [1, 1, -2, 3, -2, 1]),
        ("s1.s2^-2.s3^2.s2^-1.s1",     [1, -2, -2, 3, 3, -2, 1]),
        ("s1^3.s2^-1.s3",              [1, 1, 1, -2, 3]),
        ("s1.s2^-1.s3^3",              [1, -2, 3, 3, 3]),
        ("s1.s2^-3.s3",                [1, -2, -2, -2, 3]),
        ("(s1.s2^-1.s3)^3",            [1, -2, 3, 1, -2, 3, 1, -2, 3]),
        ("s1^2.s2^-2.s3^2.s1.s2^-1",   [1, 1, -2, -2, 3, 3, 1, -2]),
        ("s1.s3.s2^-1.s1.s3",          [1, 3, -2, 1, 3]),
        ("s3^2.s2^-1.s1^2.s2^-1.s3",   [3, 3, -2, 1, 1, -2, 3]),
        ("s1.s2^-1.s3.s1.s2^-1.s3",    [1, -2, 3, 1, -2, 3]),
        ("s1^2.s3^2.s2^-2",            [1, 1, 3, 3, -2, -2]),
        ("s1.s2^-1.s3.s2^-1.s1.s3",    [1, -2, 3, -2, 1, 3]),
    ]

    for name, word in b4_braids:
        total += 1
        try:
            data = test_braid(name, word, 4)
            an = data['alexander_norm']
            g = data['genus']
            np1d = data['alexander_poly_newton_1d']
            hull = data['charpoly_newton_2d']
            print(f"  {name:<40} norm={str(an):<4} genus={str(g):<6} NP=[{np1d[0]},{np1d[1]}]  hull={len(hull)}v  OK")
            pass_count += 1
        except Exception as e:
            print(f"  {name:<40} FAIL: {e}")

    # =========================================================
    # Detailed output for a few examples
    # =========================================================
    print()
    print("=" * 65)
    print("Detailed examples")
    print("=" * 65)

    examples = [
        ("Figure-eight: s1.s2^-1.s1.s2^-1 in B_3", [1, -2, 1, -2], 3),
        ("Alternating B_4: s1.s2^-1.s3", [1, -2, 3], 4),
        ("Alternating B_4: (s1.s2^-1.s3)^2", [1, -2, 3, 1, -2, 3], 4),
    ]

    for name, word, n in examples:
        print(f"\n--- {name} ---")
        data = fibered_face_from_alexander(word, n)
        print_fibered_face(data)

    # =========================================================
    # Summary
    # =========================================================
    print()
    print("=" * 65)
    print(f"TOTAL: {pass_count}/{total} passed")
    if pass_count == total:
        print("All alternating braids: Alexander norm computed successfully.")
        print("For these braids, Alexander norm = Thurston norm (author's theorem).")
    print("=" * 65)


if __name__ == '__main__':
    main()
