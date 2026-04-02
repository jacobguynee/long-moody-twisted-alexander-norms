#!/usr/bin/env python3
"""
test_alexander_equals_thurston.py

Compute the obvious fibered face of the Thurston norm unit ball for
alternating braids in B_3 and B_4 using all generators (pseudo-Anosov).

Runs the veering computation on each braid and reports:
  - taut isosig of the layered veering triangulation
  - face rays (cone over the fibered face)
  - taut polynomial
  - Alexander norm (from Burau, for reference)

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
    succeeded = 0
    errored = 0

    for group_name, braids in [("B_3", BRAIDS_B3), ("B_4", BRAIDS_B4)]:
        print(f"\n{'='*60}")
        print(f"Alternating {group_name} braids")
        print(f"{'='*60}")

        for name, word, n in braids:
            an = alexander_norm(word, n)

            try:
                data = obvious_fibered_face(word)
                print(f"\n  {name}")
                print(f"    alexander_norm:  {an}")
                print(f"    taut_isosig:     {data['taut_isosig']}")
                print(f"    face_rays:       {data['face_rays']}")
                print(f"    taut_polynomial: {data['taut_polynomial']}")
                succeeded += 1
            except Exception as e:
                print(f"\n  {name}")
                print(f"    alexander_norm:  {an}")
                print(f"    ERROR: {e}")
                errored += 1

    print(f"\n{'='*60}")
    print(f"{succeeded} succeeded, {errored} errored, {succeeded + errored} total")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
