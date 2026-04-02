#!/usr/bin/env python3
"""
compute_fibered_face.py

Compute the obvious fibered face of the Thurston norm unit ball for a
pseudo-Anosov braid, using layered veering triangulations.

This is the actual Thurston norm computation, not an Alexander norm proxy.
It uses:
  - SnapPy: build the mapping torus of the braid
  - Regina: enumerate taut angle structures
  - veering: identify layered structures, compute the cone over the
    fibered face (Landry-Minsky-Taylor), compute the taut polynomial

REQUIREMENTS (run inside SageMath's Python):
    sage -pip install snappy
    sage -pip install regina
    sage -pip install veering
    sage -pip install flipper   # optional, for direct veering triangulations

See README.md for full installation instructions.

Usage:
    sage -python compute_fibered_face.py
    sage -python compute_fibered_face.py 1 -2 3

References:
    - Landry, Minsky, Taylor (2021): the carried cone of a layered
      veering triangulation is exactly the cone over the fibered face.
    - Agol (2011): ideal triangulations of pseudo-Anosov mapping tori.
"""

import sys

import snappy
import regina

from veering.taut import (
    taut_regina_angle_struct_to_taut_struct,
    isosig_from_tri_angle,
)
from veering.taut_polytope import is_layered, cone_in_homology
from veering.taut_polynomial import taut_polynomial_via_fox_calculus


def braid_spec(word):
    """Convert a braid word to SnapPy's braid notation.

    word: list of signed generator indices, e.g. [1, -2, 3].
    """
    return "Braid:[" + ",".join(map(str, word)) + "]"


def obvious_fibered_face(word):
    """Compute the obvious fibered face of the Thurston norm ball.

    For a pseudo-Anosov n-braid beta:
      1. Build the mapping torus M_beta using SnapPy.
      2. Convert to a Regina triangulation.
      3. Find layered taut angle structures (exist when beta is pA).
      4. The layered veering triangulation carries a cone that is exactly
         the cone over the fibered face (Landry-Minsky-Taylor).
      5. Compute the taut polynomial (= Teichmüller polynomial up to a
         unit in the layered case).

    Parameters
    ----------
    word : list of int
        Braid word as signed generator indices, e.g. [1, -2, 3].

    Returns
    -------
    dict with:
        manifold        : SnapPy Manifold object
        taut_isosig     : str, isomorphism signature of layered taut triangulation
        face_rays       : rays spanning the cone over the obvious fibered face
                          (in veering's homology basis — see caveat below)
        taut_polynomial : the taut polynomial via Fox calculus
        all_layered     : list of all layered taut isosigs found

    Notes
    -----
    The face_rays are in the veering homology basis, which is dual to the
    basis used for taut/veering polynomials. You may need a change-of-basis
    step to convert to geometric meridian/longitude coordinates.
    """
    # 1. Build mapping torus
    M = snappy.Manifold(braid_spec(word))

    # 2. Bridge SnapPy -> Regina
    T = regina.SnapPeaTriangulation(M._to_string())

    # 3. Find layered taut angle structures
    layered_sigs = []
    for regina_angle in regina.AngleStructures(T, True):  # tautOnly=True
        angle = taut_regina_angle_struct_to_taut_struct(regina_angle)
        if is_layered(T, angle):
            sig = isosig_from_tri_angle(T, angle)
            layered_sigs.append(sig)

    if not layered_sigs:
        raise RuntimeError(
            "No layered taut structure found on this triangulation. "
            "If the braid is pseudo-Anosov, try rebuilding with "
            "flipper.bundle() to get Agol's veering triangulation directly."
        )

    # 4. Cone over the fibered face
    sig = layered_sigs[0]
    rays = cone_in_homology(sig)

    # 5. Taut polynomial
    theta = taut_polynomial_via_fox_calculus(sig)

    return {
        "manifold": M,
        "taut_isosig": sig,
        "face_rays": rays,
        "taut_polynomial": theta,
        "all_layered": layered_sigs,
    }


def obvious_fibered_face_flipper(surface_word, surface='S_1_1'):
    """Alternative: use flipper to construct the veering triangulation directly.

    For a pseudo-Anosov mapping class, flipper.bundle() builds Agol's veering
    triangulation of the mapping torus without needing to search for layered
    structures.

    Parameters
    ----------
    surface_word : str
        Mapping class in flipper notation, e.g. 'abcABC'.
    surface : str
        Surface identifier for flipper.load(), e.g. 'S_1_1', 'S_0_4'.

    Returns
    -------
    dict with taut_isosig, face_rays, taut_polynomial.
    """
    import flipper

    S = flipper.load(surface)
    h = S.mapping_class(surface_word)
    T = h.bundle()

    sig = T.isosig()
    rays = cone_in_homology(sig)
    theta = taut_polynomial_via_fox_calculus(sig)

    return {
        "taut_isosig": sig,
        "face_rays": rays,
        "taut_polynomial": theta,
    }


# ---- CLI interface ----

def main():
    if len(sys.argv) > 1:
        word = [int(x) for x in sys.argv[1:]]
    else:
        # Default examples: alternating braids using all generators
        examples = [
            [1, -2],             # B_3, simplest alternating
            [1, -2, 3],          # B_4, simplest alternating using all gens
            [1, -2, 1, -2],      # figure-eight knot
            [1, -2, 3, 1, -2, 3],  # (s1.s2^-1.s3)^2
        ]
        for word in examples:
            print(f"\n{'='*60}")
            print(f"Braid word: {word}")
            print(f"{'='*60}")
            try:
                data = obvious_fibered_face(word)
                print(f"  Taut isosig:     {data['taut_isosig']}")
                print(f"  Face rays:       {data['face_rays']}")
                print(f"  Taut polynomial: {data['taut_polynomial']}")
                print(f"  # layered:       {len(data['all_layered'])}")
            except Exception as e:
                print(f"  Error: {e}")
        return

    print(f"Braid word: {word}")
    data = obvious_fibered_face(word)
    print(f"Taut isosig:     {data['taut_isosig']}")
    print(f"Face rays:       {data['face_rays']}")
    print(f"Taut polynomial: {data['taut_polynomial']}")
    print(f"# layered:       {len(data['all_layered'])}")


if __name__ == '__main__':
    main()
