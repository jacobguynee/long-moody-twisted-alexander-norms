#!/usr/bin/env python3
"""
test_alternating_braids.py

Compute Alexander norms for alternating braids in B_3 and B_4.

For alternating braids using all generators (pseudo-Anosov), the Alexander
norm equals the Thurston norm. This provides ground-truth values to validate
the veering Thurston norm computation (compute_fibered_face.py).

Validation strategy:
  1. Run this script to get Alexander norms (runs anywhere with SymPy).
  2. Run compute_fibered_face.py on the same braids (requires SageMath +
     snappy + regina + veering).
  3. Check that the Thurston norm from veering equals the Alexander norm
     from this script on every alternating braid.

Usage:
    python test_alternating_braids.py
"""

from sympy import Rational

from alexander_thurston import (
    compute_alexander_data,
    print_alexander_data,
    newton_polytope_1d,
    alexander_norm_1d,
)


def test_braid(name, word, n):
    """Compute Alexander norm for a single braid. Returns data dict."""
    data = compute_alexander_data(word, n)
    alex_norm = data['alexander_norm']
    genus = data['genus']

    assert alex_norm >= 0, f"Negative Alexander norm for {name}"

    # For pseudo-Anosov (uses all generators), norm should be positive
    if len(set(abs(w) for w in word)) == n - 1:
        assert alex_norm > 0, f"Zero Alexander norm for pseudo-Anosov {name}"

    assert genus == Rational(alex_norm, 2), f"Genus mismatch for {name}"

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
    print("Alexander norm = Thurston norm for these braids")
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
            print(f"  {name:<40} alex_norm={str(an):<4} genus={str(g):<6} NP=[{np1d[0]},{np1d[1]}]  OK")
            pass_count += 1
        except Exception as e:
            print(f"  {name:<40} FAIL: {e}")

    # =========================================================
    # B_4 alternating braids: words in {s1, s2^{-1}, s3}
    # =========================================================
    print()
    print("=" * 65)
    print("B_4 alternating braids (words in s1, s2^-1, s3)")
    print("Alexander norm = Thurston norm for these braids")
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
            print(f"  {name:<40} alex_norm={str(an):<4} genus={str(g):<6} NP=[{np1d[0]},{np1d[1]}]  hull={len(hull)}v  OK")
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
        data = compute_alexander_data(word, n)
        print_alexander_data(data)

    # =========================================================
    # Summary
    # =========================================================
    print()
    print("=" * 65)
    print(f"TOTAL: {pass_count}/{total} passed")
    if pass_count == total:
        print("All Alexander norms computed successfully.")
        print()
        print("To validate the veering Thurston norm computation:")
        print("  1. Install: sage -pip install snappy regina veering")
        print("  2. Run: sage -python compute_fibered_face.py")
        print("  3. Check: Thurston norm == Alexander norm for each braid above")
    print("=" * 65)


if __name__ == '__main__':
    main()
