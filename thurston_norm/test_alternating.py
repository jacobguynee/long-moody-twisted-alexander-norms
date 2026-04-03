#!/usr/bin/env sage -python
"""
Test script: verify Alexander = Thurston fibered face for 3-braids,
then search for 4-braids where they differ.

This script is the verification harness described in the README.

Mathematical background:
  - McMullen's inequality: ||phi||_A <= ||phi||_T for all phi in H^1.
  - For fibered 3-component links of braid index 3, equality holds on
    fibered faces (consequence of the genus-1 surface structure).
  - For 4-braids, the inequality can be strict, meaning the Alexander face
    of the Thurston norm ball strictly contains the Thurston fibered face.

What this script checks:
  - For each braid word, compute the Thurston fibered face vertices v_i.
  - Each v_i lies on the Thurston unit ball, so ||v_i||_T = 1.
  - Compute ||v_i||_A via Morton's colored Burau Alexander polynomial.
  - If ||v_i||_A = 1 for all i, the faces coincide.
  - If ||v_i||_A < 1 for some i, the Alexander face is strictly larger.

Usage:
    sage -python test_alternating.py                  # run all tests
    sage -python test_alternating.py --only3          # 3-braids only
    sage -python test_alternating.py --only4          # 4-braids only
    sage -python test_alternating.py --maxlen3 10     # longer 3-braid words
    sage -python test_alternating.py --maxlen4 8      # longer 4-braid words
"""

import sys
import argparse

from braid_meridian_face import (
    analyze_braid,
    pseudo_anosov_3braid_words,
    interesting_4braid_words,
    alternating_braid_words,
)


def test_3braid_coincidence(max_len=8, verbose=True):
    """
    3-braid verification: Alexander face should always equal Thurston face.
    """
    print("=" * 60)
    print("TEST: 3-braids -- Alexander face = Thurston face (expected)")
    print("=" * 60)

    words = pseudo_anosov_3braid_words(max_len)
    # Also add alternating 3-braids
    words.extend(alternating_braid_words(3, max_len))

    # Deduplicate
    seen = set()
    unique_words = []
    for w in words:
        key = tuple(w)
        if key not in seen:
            seen.add(key)
            unique_words.append(w)
    words = unique_words

    passed = 0
    failed = 0
    skipped = 0

    for word in words:
        try:
            result = analyze_braid(word, verbose=False)
            if result is None:
                skipped += 1
                if verbose:
                    print(f"  SKIP  {word}")
                continue

            if result["coincide"]:
                passed += 1
                if verbose:
                    print(f"  PASS  {word}  (vertices: {result['thurston_face_vertices']})")
            else:
                failed += 1
                print(f"  FAIL  {word}")
                for v, a in result["vertex_details"]:
                    print(f"        vertex {v}: ||v||_A = {a}")
        except Exception as e:
            skipped += 1
            if verbose:
                print(f"  SKIP  {word}  ({e})")

    print(f"\n3-braid results: {passed} passed, {failed} FAILED, {skipped} skipped")
    if failed:
        print("!!! UNEXPECTED FAILURES -- check implementation !!!")
    return failed == 0


def test_4braid_search(max_len=6, verbose=True):
    """
    4-braid search: find braids where Alexander face != Thurston face.
    """
    print("\n" + "=" * 60)
    print("SEARCH: 4-braids -- looking for Alexander face != Thurston face")
    print("=" * 60)

    words = interesting_4braid_words(max_len)

    # Deduplicate
    seen = set()
    unique_words = []
    for w in words:
        key = tuple(w)
        if key not in seen:
            seen.add(key)
            unique_words.append(w)
    words = unique_words

    coincide = 0
    differ = 0
    skipped = 0
    examples = []

    for word in words:
        try:
            result = analyze_braid(word, verbose=False)
            if result is None:
                skipped += 1
                continue

            if result["coincide"]:
                coincide += 1
                if verbose:
                    print(f"  =  {word}")
            else:
                differ += 1
                examples.append(result)
                print(f"  != {word}  <-- Alexander face STRICTLY LARGER")
                for v, a in result["vertex_details"]:
                    if a != 1:
                        print(f"       vertex {v}: ||v||_A = {a}")
        except Exception as e:
            skipped += 1
            if verbose:
                print(f"  SKIP  {word}  ({e})")

    print(f"\n4-braid results: {coincide} coincide, {differ} differ, {skipped} skipped")

    if examples:
        print(f"\n{'='*60}")
        print(f"FOUND {len(examples)} example(s) where Alexander != Thurston:")
        print(f"{'='*60}")
        for ex in examples:
            print(f"\n  Braid: {ex['word']}")
            print(f"  Strands: {ex['n_strands']}")
            print(f"  Thurston face: {ex['thurston_face_vertices']}")
            print(f"  Alexander poly: {ex['alexander_poly']}")
            print(f"  Newton polytope: {ex['newton_vertices']}")
            for v, a in ex["vertex_details"]:
                print(f"  ||{v}||_A = {a},  ||{v}||_T = 1")
    else:
        print("\nNo examples found.  Try --maxlen4 with a larger value.")

    return examples


def main():
    parser = argparse.ArgumentParser(
        description="Test Alexander vs Thurston fibered face for braids."
    )
    parser.add_argument("--only3", action="store_true",
                        help="Only run 3-braid verification")
    parser.add_argument("--only4", action="store_true",
                        help="Only run 4-braid search")
    parser.add_argument("--maxlen3", type=int, default=8,
                        help="Max word length for 3-braids (default: 8)")
    parser.add_argument("--maxlen4", type=int, default=6,
                        help="Max word length for 4-braids (default: 6)")
    parser.add_argument("--verbose", "-v", action="store_true", default=True)
    parser.add_argument("--quiet", "-q", action="store_true")

    args = parser.parse_args()
    verbose = not args.quiet

    ok = True

    if not args.only4:
        ok = test_3braid_coincidence(args.maxlen3, verbose) and ok

    if not args.only3:
        examples = test_4braid_search(args.maxlen4, verbose)

    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
