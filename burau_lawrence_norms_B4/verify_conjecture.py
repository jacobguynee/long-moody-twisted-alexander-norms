#!/usr/bin/env python3
"""
verify_conjecture.py

Main verification script for the conjecture:

    NP(charpoly(Lawrence k=2)) / 4  =  NP(charpoly(Burau))

for braids in B_4, where NP denotes Newton polytope in (eigenvalue x,
parameter q) coordinates.

Tests a diverse collection of braids:
  - Alternating braids (monoid in s1, s2^-1, s3)
  - f^a g^b families (almost-reducible pseudo-Anosovs)
  - General words f^{a1}g^{b1}...f^{ak}g^{bk}
  - Random braid words

Usage:
    python verify_conjecture.py
"""

import sys
import random
from sympy import cancel, Rational

from representations import build_burau_B4, build_lawrence_B4
from newton_polytope import (
    compare_newton_polytopes, mat_power, eval_word
)


def main():
    random.seed(42)

    print('Building B_4 representations...')
    sys.stdout.flush()
    bg = build_burau_B4()
    bgi = [cancel(g.inv()) for g in bg]
    lg = build_lawrence_B4()
    lgi = [cancel(g.inv()) for g in lg]
    print(f'  Burau:    {bg[0].shape[0]}x{bg[0].shape[1]}')
    print(f'  Lawrence: {lg[0].shape[0]}x{lg[0].shape[1]}')
    print('Done.\n')
    sys.stdout.flush()

    scale = Rational(1, 4)
    pass_count = 0
    fail_count = 0

    def test(name, word=None, Mb=None, Ml=None):
        nonlocal pass_count, fail_count
        if word is not None:
            Mb = eval_word(word, bg, bgi)
            Ml = eval_word(word, lg, lgi)
        match, hull_b, hull_l_s = compare_newton_polytopes(Mb, Ml, scale)
        if match:
            print(f'  {name:<55} YES  {hull_b}')
            pass_count += 1
        else:
            print(f'  {name:<55} *** NO ***')
            print(f'    NP(Burau):      {hull_b}')
            print(f'    NP(Lawrence)/4: {hull_l_s}')
            fail_count += 1
        sys.stdout.flush()

    # =========================================================
    # Part 1: Alternating braids (words in s1, s2^-1, s3)
    # =========================================================
    print('=== Part 1: Alternating braids ===\n')

    alternating = [
        ('s1.s2^-1.s3',              [1, -2, 3]),
        ('s3.s2^-1.s1',              [3, -2, 1]),
        ('s1.s2^-1.s3.s2^-1',        [1, -2, 3, -2]),
        ('s1^2.s2^-1.s3',            [1, 1, -2, 3]),
        ('s1.s2^-2.s3',              [1, -2, -2, 3]),
        ('s1.s2^-1.s3^2',            [1, -2, 3, 3]),
        ('s1^2.s2^-2.s3^2',          [1, 1, -2, -2, 3, 3]),
        ('(s1.s2^-1.s3)^2',          [1, -2, 3, 1, -2, 3]),
        ('s1^2.s2^-1.s3.s2^-1.s1',   [1, 1, -2, 3, -2, 1]),
        ('s1.s2^-2.s3^2.s2^-1.s1',   [1, -2, -2, 3, 3, -2, 1]),
    ]
    for name, word in alternating:
        test(name, word)

    # =========================================================
    # Part 2: f^a g^b families
    # =========================================================
    print('\n=== Part 2: f^a g^b families ===\n')

    families = [
        ('s1, s3',            [1],        [3]),
        ('s1, s2.s3.s2^-1',   [1],        [2, 3, -2]),
        ('s1.s2, s3',         [1, 2],     [3]),
        ('s1.s2, s2^-1.s3',   [1, 2],     [-2, 3]),
        ('s1.s3, s2',         [1, 3],     [2]),
        ('s2.s1.s2^-1, s3',   [2, 1, -2], [3]),
        ('[s1,s2], s3',       [1,2,-1,-2],[3]),
        ('s1.s2.s1, s3',      [1, 2, 1],  [3]),
    ]

    exp_pairs = [
        (5, -5), (5, 5), (10, -3), (3, -10), (10, 10), (15, -5), (20, -1),
    ]

    for fam_name, fw, gw in families:
        fb = eval_word(fw, bg, bgi); gb = eval_word(gw, bg, bgi)
        fl = eval_word(fw, lg, lgi); gl = eval_word(gw, lg, lgi)
        for a, b in exp_pairs:
            Mb = cancel(mat_power(fb, a) * mat_power(gb, b))
            Ml = cancel(mat_power(fl, a) * mat_power(gl, b))
            test(f'f={fam_name}: f^{a} g^{b}', Mb=Mb, Ml=Ml)

    # =========================================================
    # Part 3: General words f^{a1}g^{b1}f^{a2}g^{b2}...
    # =========================================================
    print('\n=== Part 3: General alternating words ===\n')

    from sympy import eye as sym_eye
    for fam_name, fw, gw in families[:4]:
        fb = eval_word(fw, bg, bgi); gb = eval_word(gw, bg, bgi)
        fl = eval_word(fw, lg, lgi); gl = eval_word(gw, lg, lgi)
        for trial in range(3):
            k = random.randint(2, 4)
            segs = [(random.choice([-1,1])*random.randint(2,5),
                     random.choice([-1,1])*random.randint(2,5)) for _ in range(k)]
            Mb = sym_eye(3); Ml = sym_eye(12)
            for ai, bi in segs:
                Mb = cancel(Mb * mat_power(fb, ai) * mat_power(gb, bi))
                Ml = cancel(Ml * mat_power(fl, ai) * mat_power(gl, bi))
            seg_str = '.'.join(f'f^{a}g^{b}' for a, b in segs)
            test(f'f={fam_name}: {seg_str}', Mb=Mb, Ml=Ml)

    # =========================================================
    # Part 4: Random braid words
    # =========================================================
    print('\n=== Part 4: Random braid words ===\n')

    for trial in range(10):
        length = random.randint(8, 16)
        word = [random.choice([-1, 1]) * random.randint(1, 3) for _ in range(length)]
        if len(set(abs(w) for w in word)) < 3:
            continue
        wstr = ''.join(f's{abs(w)}{"" if w>0 else "i"}' for w in word[:6]) + f'...(L={length})'
        test(wstr, word)

    # =========================================================
    # Summary
    # =========================================================
    print(f'\n{"="*60}')
    print(f'TOTAL: {pass_count} passed, {fail_count} failed')
    if fail_count == 0:
        print('All Newton polytopes match: NP(Lawrence k=2)/4 = NP(Burau)')
    print(f'{"="*60}')


if __name__ == '__main__':
    main()
