#!/usr/bin/env sage -python
"""
Compute the obvious fibered face of the Thurston norm ball for the mapping
torus of a pseudo-Anosov braid, in the meridian homology basis.

Pipeline:
  1. Build mapping torus via SnapPy's Braid:[...] specification.
  2. Recover the closed-braid + axis link via exterior_to_link(check_answer=True).
  3. Use L.exterior() so peripheral structures (meridians) match Alexander variables.
  4. Use TNorm to compute the Thurston norm ball in the meridian-dual H_2 basis.
  5. Identify the obvious fibered face.
  6. Compare with the Alexander norm ball from Morton's colored Burau formula.

The meridian basis is the same one used for the multivariable Alexander
polynomial (colored Burau / Long-Moody), so both norms live in the same
coordinate system and can be compared directly.

Requirements (inside SageMath):
    sage -pip install tnorm snappy

Usage:
    sage -python braid_meridian_face.py --word 1 -2 1 -2
    sage -python braid_meridian_face.py --word 1 -2 1 -2 --verbose
    sage -python braid_meridian_face.py --test3 --maxlen 8
    sage -python braid_meridian_face.py --test4 --maxlen 6
    sage -python braid_meridian_face.py --search4 --maxlen 8
"""

import sys
import argparse
from itertools import product as iterproduct

from colored_burau import (
    braid_permutation,
    permutation_cycles,
    morton_alexander_poly,
    newton_polytope,
    alexander_norm,
    alexander_norm_ball,
)


# ---------------------------------------------------------------------------
# SnapPy / TNorm interface
# ---------------------------------------------------------------------------

def braid_spec(word):
    """Format a braid word for SnapPy's Braid:[...] input."""
    return "Braid:[" + ",".join(map(str, word)) + "]"


def build_link_exterior(word, verbose=False):
    """
    From a braid word, build the closed-braid + axis link exterior.

    SnapPy's Braid:[...] builds the mapping torus of the punctured-disk
    monodromy.  We recover the link diagram via exterior_to_link so that
    the peripheral structure (meridians) is canonical.

    Returns (M, L, n_strands, perm_cycles).
    """
    import snappy

    n = max(abs(g) for g in word) + 1

    if verbose:
        print(f"  {n}-strand braid: {word}")

    M0 = snappy.Manifold(braid_spec(word))

    if verbose:
        print(f"  Mapping torus cusps: {M0.num_cusps()}")

    # Recover link diagram preserving cusp/component order.
    L = M0.exterior_to_link(check_answer=True)

    # Use L.exterior() so peripheral basis matches the link components.
    M = L.exterior()

    perm = braid_permutation(word, n)
    cycles = permutation_cycles(perm)

    if verbose:
        print(f"  Link components: {L.num_components()}")
        print(f"  Permutation cycles: {cycles}")
        print(f"  Cusps: {M.num_cusps()}")

    return M, L, n, cycles


def compute_thurston_norm(M, verbose=False):
    """
    Compute the Thurston norm ball via TNorm.

    TNorm uses the peripheral (meridian) basis when given a link exterior,
    which is exactly the basis we need for comparison with the Alexander norm.

    Returns the TNorm object.
    """
    from tnorm import TNorm

    if verbose:
        print("  Computing Thurston norm ball...")

    T = TNorm(M)

    if verbose:
        print(f"  b_2 = {T.betti}")
        print(f"  Fibered faces: {len(T.fibered_faces)}")

    return T


def get_fibered_face_vertices(T, verbose=False):
    """
    Extract vertices of fibered faces from a TNorm result.

    Returns list of (face_index, vertices) for each fibered face.
    """
    results = []
    for i, face in enumerate(T.fibered_faces):
        verts = [tuple(v) for v in face.vertices]
        if verbose:
            print(f"  Fibered face {i}: {verts}")
        results.append((i, verts))
    return results


def get_obvious_fibered_face(T, verbose=False):
    """
    Get the 'obvious' fibered face -- for a mapping torus there is typically
    one natural fibered face corresponding to the fiber of the open book.

    If there are multiple fibered faces, return the top-dimensional one
    (most vertices).  For mapping tori of pseudo-Anosov braids, Thurston
    shows the fiber class lies in the interior of a unique top-dimensional
    fibered face.

    Returns vertices as list of tuples.
    """
    faces = get_fibered_face_vertices(T, verbose=verbose)

    if not faces:
        raise RuntimeError("No fibered faces found!")

    if len(faces) == 1:
        return faces[0][1]

    # Return the face with the most vertices (largest dimension).
    best = max(faces, key=lambda x: len(x[1]))
    return best[1]


# ---------------------------------------------------------------------------
# Face comparison
# ---------------------------------------------------------------------------

def faces_coincide(thurston_verts, alex_poly, delta, verbose=False):
    """
    Check whether the Thurston fibered face equals the Alexander face
    containing it.

    McMullen's inequality says ||phi||_A <= ||phi||_T for all phi in H^1.
    Equality on a face means the Alexander and Thurston faces coincide.

    Strategy: for each vertex v of the Thurston unit ball face, v has
    ||v||_T = 1.  Compute ||v||_A.  If ||v||_A = 1 for all vertices,
    the faces coincide.  If ||v||_A < 1 for some vertex, the Alexander
    face strictly contains the Thurston face.

    Returns (coincide, details) where details is a list of (vertex, alex_norm).
    """
    details = []
    coincide = True

    for v in thurston_verts:
        a = alexander_norm(delta, v)
        details.append((v, a))
        if a != 1:
            coincide = False
            if verbose:
                print(f"    vertex {v}: ||v||_A = {a} != 1")

    if verbose:
        if coincide:
            print("    All vertices have ||v||_A = 1: faces COINCIDE")
        else:
            print("    Some vertices have ||v||_A != 1: Alexander face is STRICTLY LARGER")

    return coincide, details


# ---------------------------------------------------------------------------
# Single braid analysis
# ---------------------------------------------------------------------------

def analyze_braid(word, verbose=False):
    """
    Full pipeline for a single braid word.

    Returns dict with all computed data, or None if the computation fails.
    """
    n = max(abs(g) for g in word) + 1

    if verbose:
        print(f"\n{'='*60}")
        print(f"Analyzing braid: {word}")
        print(f"{'='*60}")

    # Step 1: Alexander polynomial via Morton/colored Burau
    if verbose:
        print("\n--- Alexander polynomial (Morton) ---")

    delta, ring, info = morton_alexander_poly(word, n, include_axis=True)
    newton_verts, newton_poly = newton_polytope(delta)

    if verbose:
        print(f"  Delta = {delta}")
        print(f"  Ring = {ring}")
        print(f"  Closure components: {info['num_closure_components']}")
        print(f"  Cycles: {info['cycles']}")
        print(f"  Newton polytope vertices: {newton_verts}")

    # Step 2: Link exterior and Thurston norm
    if verbose:
        print("\n--- Thurston norm (TNorm) ---")

    try:
        M, L, n_strands, cycles = build_link_exterior(word, verbose=verbose)
        T = compute_thurston_norm(M, verbose=verbose)
        thurston_verts = get_obvious_fibered_face(T, verbose=verbose)
    except Exception as e:
        if verbose:
            print(f"  TNorm failed: {e}")
        return None

    # Step 3: Compare faces
    if verbose:
        print("\n--- Face comparison ---")
        print(f"  Thurston face vertices: {thurston_verts}")

    coincide, details = faces_coincide(thurston_verts, None, delta, verbose=verbose)

    return {
        "word": word,
        "n_strands": n,
        "cycles": cycles,
        "alexander_poly": delta,
        "newton_vertices": newton_verts,
        "thurston_face_vertices": thurston_verts,
        "coincide": coincide,
        "vertex_details": details,
        "manifold": M,
        "tnorm": T,
    }


# ---------------------------------------------------------------------------
# Braid generators for testing
# ---------------------------------------------------------------------------

def all_braid_words(n_strands, max_len, min_len=None):
    """
    Generate all braid words on n_strands strands up to max_len generators.

    Yields words as lists of nonzero integers in {-(n-1),...,-1, 1,...,n-1}.
    """
    if min_len is None:
        min_len = n_strands  # need at least n generators to be pseudo-Anosov

    gens = list(range(1, n_strands)) + list(range(-(n_strands - 1), 0))

    for length in range(min_len, max_len + 1):
        for w in iterproduct(gens, repeat=length):
            yield list(w)


def alternating_braid_words(n_strands, max_len):
    """
    Generate alternating braid words on n_strands strands.

    An alternating braid alternates the signs of consecutive generators.
    For 3-braids: sigma_1^{+/-1} sigma_2^{-/+1} sigma_1^{+/-1} ...
    """
    gens = list(range(1, n_strands))  # [1, ..., n-1]
    words = []

    for length in range(n_strands, max_len + 1):
        # Pattern: cycle through generators, alternating sign
        w = []
        for i in range(length):
            gen = gens[i % len(gens)]
            sign = 1 if i % 2 == 0 else -1
            w.append(sign * gen)
        words.append(w)

        # Opposite sign pattern
        w2 = [-g for g in w]
        words.append(w2)

    return words


def pseudo_anosov_3braid_words(max_len):
    """
    Generate pseudo-Anosov 3-braid words.

    A 3-braid sigma_1^{a_1} sigma_2^{b_1} ... is pseudo-Anosov if and only if
    all exponents have the same sign (positive or negative) and the word uses
    both generators.  This is a theorem of Birman-Menasco.

    For testing, we generate words of the form sigma_1^a sigma_2^b with
    a, b > 0, and also words with mixed positive exponent blocks.
    """
    words = []
    # All-positive exponent words: sigma_1^a1 sigma_2^b1 sigma_1^a2 sigma_2^b2 ...
    for total_len in range(3, max_len + 1):
        # Simple: alternating single generators
        w = []
        for i in range(total_len):
            w.append(1 if i % 2 == 0 else 2)
        words.append(w)

        # Negative version
        words.append([-g for g in w])

    return words


def interesting_4braid_words(max_len):
    """
    Generate 4-braid words to search for Alexander != Thurston.

    We try various patterns that are likely pseudo-Anosov.
    """
    words = []

    # Alternating braids
    words.extend(alternating_braid_words(4, max_len))

    # All-positive words mixing all three generators
    for length in range(4, max_len + 1):
        patterns = [
            [1, 2, 3],            # sigma_1 sigma_2 sigma_3 repeated
            [1, 3, 2],            # different ordering
            [1, 2, 1, 3],         # with repeats
            [1, 2, 3, 2],
            [1, 1, 2, 3],
            [1, 2, 2, 3],
            [1, 2, 3, 3],
            [1, 1, 2, 2, 3, 3],
        ]
        for pat in patterns:
            w = []
            for i in range(length):
                w.append(pat[i % len(pat)])
            if len(w) == length:
                words.append(w)
                words.append([-g for g in w])

    return words


# ---------------------------------------------------------------------------
# Test harnesses
# ---------------------------------------------------------------------------

def test_3braids(max_len=8, verbose=True):
    """
    Verify that for pseudo-Anosov 3-braids, the Alexander face equals the
    Thurston fibered face.

    For 3-braids this is known to be true (the Alexander and Thurston norms
    coincide on fibered faces for fibered links of braid index 3).
    """
    print("=" * 60)
    print("3-BRAID TEST: Alexander face = Thurston fibered face")
    print("=" * 60)

    words = pseudo_anosov_3braid_words(max_len)
    passed = 0
    failed = 0
    skipped = 0

    for word in words:
        label = f"3-braid {word}"
        try:
            result = analyze_braid(word, verbose=False)
            if result is None:
                skipped += 1
                if verbose:
                    print(f"  SKIP: {label}")
                continue

            if result["coincide"]:
                passed += 1
                if verbose:
                    print(f"  PASS: {label}")
                    print(f"        Thurston vertices: {result['thurston_face_vertices']}")
            else:
                failed += 1
                if verbose:
                    print(f"  FAIL: {label}")
                    for v, a in result["vertex_details"]:
                        print(f"        vertex {v}: ||v||_A = {a}")

        except Exception as e:
            skipped += 1
            if verbose:
                print(f"  SKIP: {label} -- {e}")

    print(f"\nResults: {passed} passed, {failed} failed, {skipped} skipped")
    if failed > 0:
        print("UNEXPECTED: Some 3-braids failed the coincidence test!")
        return False
    else:
        print("All 3-braid tests passed (as expected).")
        return True


def test_4braids(max_len=6, verbose=True):
    """
    Test 4-braids for Alexander vs Thurston face coincidence.

    Unlike 3-braids, some 4-braids may have Alexander face strictly larger
    than the Thurston face.  This function reports all cases.
    """
    print("=" * 60)
    print("4-BRAID TEST: Alexander face vs Thurston fibered face")
    print("=" * 60)

    words = interesting_4braid_words(max_len)
    coincide_count = 0
    differ_count = 0
    skipped = 0

    for word in words:
        label = f"4-braid {word}"
        try:
            result = analyze_braid(word, verbose=False)
            if result is None:
                skipped += 1
                continue

            if result["coincide"]:
                coincide_count += 1
                if verbose:
                    print(f"  COINCIDE: {label}")
            else:
                differ_count += 1
                print(f"  DIFFER:   {label}")
                for v, a in result["vertex_details"]:
                    if a != 1:
                        print(f"            vertex {v}: ||v||_A = {a} < 1")
                print(f"            Alexander face is STRICTLY LARGER than Thurston face")

        except Exception as e:
            skipped += 1
            if verbose:
                print(f"  SKIP: {label} -- {e}")

    print(f"\nResults: {coincide_count} coincide, {differ_count} differ, {skipped} skipped")
    if differ_count > 0:
        print(f"\nFound {differ_count} braid(s) where Alexander face != Thurston face!")
        return True
    else:
        print("\nNo examples found where faces differ.  Try increasing --maxlen.")
        return False


def search_4braids_exhaustive(max_len=6, verbose=True):
    """
    Exhaustive search over short 4-braid words to find one where the
    Alexander fibered face strictly contains the Thurston fibered face.
    """
    print("=" * 60)
    print("4-BRAID EXHAUSTIVE SEARCH")
    print("=" * 60)

    found = []
    tested = 0
    skipped = 0

    for word in all_braid_words(4, max_len, min_len=4):
        tested += 1
        if tested % 100 == 0 and verbose:
            print(f"  Tested {tested} braids, found {len(found)} examples...")

        try:
            result = analyze_braid(word, verbose=False)
            if result is None:
                skipped += 1
                continue

            if not result["coincide"]:
                found.append(result)
                print(f"\n  FOUND: {word}")
                for v, a in result["vertex_details"]:
                    print(f"    vertex {v}: ||v||_A = {a}")

        except Exception:
            skipped += 1

    print(f"\nTested {tested} braids, skipped {skipped}")
    print(f"Found {len(found)} examples where Alexander face != Thurston face")

    return found


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compare Thurston and Alexander fibered faces for braid mapping tori."
    )
    parser.add_argument("--word", nargs="+", type=int,
                        help="Braid word as space-separated integers")
    parser.add_argument("--verbose", "-v", action="store_true")
    parser.add_argument("--test3", action="store_true",
                        help="Run 3-braid verification (should all coincide)")
    parser.add_argument("--test4", action="store_true",
                        help="Run 4-braid comparison")
    parser.add_argument("--search4", action="store_true",
                        help="Exhaustive 4-braid search for non-coincidence")
    parser.add_argument("--maxlen", type=int, default=8,
                        help="Max braid word length (default: 8)")

    args = parser.parse_args()

    if args.test3:
        test_3braids(max_len=args.maxlen, verbose=args.verbose)
    elif args.test4:
        test_4braids(max_len=args.maxlen, verbose=args.verbose)
    elif args.search4:
        search_4braids_exhaustive(max_len=args.maxlen, verbose=args.verbose)
    elif args.word:
        result = analyze_braid(args.word, verbose=True)
        if result:
            print(f"\n{'='*60}")
            print(f"SUMMARY")
            print(f"{'='*60}")
            print(f"Braid: {result['word']}")
            print(f"Strands: {result['n_strands']}")
            print(f"Closure components: {len(result['cycles'])}")
            print(f"Thurston face vertices: {result['thurston_face_vertices']}")
            print(f"Alexander poly: {result['alexander_poly']}")
            print(f"Newton vertices: {result['newton_vertices']}")
            print(f"Faces coincide: {result['coincide']}")
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
