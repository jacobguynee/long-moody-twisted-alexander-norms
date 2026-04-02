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
from alexander_thurston import alexander_norm


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


def main():
    passed = 0
    failed = 0

    for group_name, braids in [("B_3", BRAIDS_B3), ("B_4", BRAIDS_B4)]:
        print(f"\n{'='*60}")
        print(f"Alternating {group_name} braids")
        print(f"{'='*60}")

        for name, word, n in braids:
            an = alexander_norm(word, n)

            try:
                data = obvious_fibered_face(word)
                theta = data['taut_polynomial']
                thurston = theta.degree()

                if an == thurston:
                    print(f"  {name:<30} alex={an}  thurston={thurston}  PASS")
                    passed += 1
                else:
                    print(f"  {name:<30} alex={an}  thurston={thurston}  FAIL")
                    print(f"    taut_isosig:     {data['taut_isosig']}")
                    print(f"    taut_polynomial: {theta}")
                    print(f"    face_rays:       {data['face_rays']}")
                    failed += 1
            except Exception as e:
                print(f"  {name:<30} alex={an}  veering ERROR: {e}")
                failed += 1

    print(f"\n{'='*60}")
    print(f"Results: {passed} passed, {failed} failed, {passed + failed} total")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
