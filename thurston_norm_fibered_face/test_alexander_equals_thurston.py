#!/usr/bin/env python3
"""
test_alexander_equals_thurston.py

For alternating braids in B_3 and B_4 using all generators (pseudo-Anosov),
the Alexander norm equals the Thurston norm. This script computes both
independently and checks they agree.

  - Thurston norm: via veering triangulations (compute_fibered_face.py)
  - Alexander norm: via reduced Burau representation (alexander_thurston.py)

Run with:  sage -python test_alexander_equals_thurston.py
"""

from compute_fibered_face import obvious_fibered_face
from alexander_thurston import compute_alexander_data


BRAIDS_B3 = [
    ("s1.s2^-1",              [1, -2],                 3),
    ("(s1.s2^-1)^2",          [1, -2, 1, -2],          3),
    ("s1^2.s2^-1",            [1, 1, -2],              3),
    ("s1.s2^-2",              [1, -2, -2],             3),
    ("s1^2.s2^-2",            [1, 1, -2, -2],          3),
    ("(s1.s2^-1)^3",          [1, -2, 1, -2, 1, -2],  3),
    ("s1^3.s2^-1",            [1, 1, 1, -2],           3),
    ("s1.s2^-3",              [1, -2, -2, -2],         3),
    ("s1^3.s2^-3",            [1, 1, 1, -2, -2, -2],  3),
]

BRAIDS_B4 = [
    ("s1.s2^-1.s3",           [1, -2, 3],              4),
    ("s1.s2^-1.s3.s2^-1",     [1, -2, 3, -2],          4),
    ("s1^2.s2^-1.s3",         [1, 1, -2, 3],           4),
    ("s1.s2^-2.s3",           [1, -2, -2, 3],          4),
    ("s1.s2^-1.s3^2",         [1, -2, 3, 3],           4),
    ("(s1.s2^-1.s3)^2",       [1, -2, 3, 1, -2, 3],   4),
    ("s1^2.s2^-2.s3^2",       [1, 1, -2, -2, 3, 3],   4),
    ("s1^3.s2^-1.s3",         [1, 1, 1, -2, 3],        4),
    ("s1.s2^-1.s3^3",         [1, -2, 3, 3, 3],        4),
    ("(s1.s2^-1.s3)^3",       [1, -2, 3, 1, -2, 3, 1, -2, 3], 4),
]


def thurston_norm_from_veering(word):
    """Extract the Thurston norm of the fiber from veering data.

    For a 1-cusped manifold (knot complement), the Thurston norm of the
    fiber equals the degree of the taut polynomial.
    """
    data = obvious_fibered_face(word)
    poly = data['taut_polynomial']
    # The taut polynomial degree gives the Thurston norm
    return poly.degree(), data


def main():
    passed = 0
    failed = 0

    for group_name, braids in [("B_3", BRAIDS_B3), ("B_4", BRAIDS_B4)]:
        print(f"\n{'='*60}")
        print(f"Alternating {group_name} braids")
        print(f"{'='*60}")

        for name, word, n in braids:
            alex_data = compute_alexander_data(word, n)
            alex_norm = alex_data['alexander_norm']

            try:
                thurston_deg, veering_data = thurston_norm_from_veering(word)
                if alex_norm == thurston_deg:
                    print(f"  {name:<30} alex={alex_norm}  thurston={thurston_deg}  PASS")
                    passed += 1
                else:
                    print(f"  {name:<30} alex={alex_norm}  thurston={thurston_deg}  FAIL (mismatch)")
                    failed += 1
            except Exception as e:
                print(f"  {name:<30} alex={alex_norm}  veering ERROR: {e}")
                failed += 1

    print(f"\n{'='*60}")
    print(f"Results: {passed} passed, {failed} failed, {passed + failed} total")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
