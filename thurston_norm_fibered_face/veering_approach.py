"""
veering_approach.py

Compute the obvious fibered face using the veering triangulation machinery
(SnapPy + Regina + veering package).

REQUIRES: sage, snappy, regina, veering (and optionally flipper).
These are NOT available in all Python environments. This file serves as
a reference implementation. See README.md for installation instructions.

For pseudo-Anosov braids, the layered veering triangulation of the mapping
torus carries a cone that is exactly the cone over the fibered face of the
Thurston norm ball (Landry-Minsky-Taylor).
"""

try:
    import snappy
    import regina
    from veering.taut import (
        taut_regina_angle_struct_to_taut_struct,
        isosig_from_tri_angle,
    )
    from veering.taut_polytope import is_layered, cone_in_homology
    from veering.taut_polynomial import taut_polynomial_via_fox_calculus
    HAS_VEERING = True
except ImportError:
    HAS_VEERING = False


def braid_spec(word):
    """Convert braid word to SnapPy braid notation."""
    return "Braid:[" + ",".join(map(str, word)) + "]"


def obvious_fibered_face_veering(word):
    """Compute the obvious fibered face via veering triangulation.

    Parameters
    ----------
    word : list of int
        Braid word as signed generator indices, e.g. [1, -2, 3].

    Returns
    -------
    dict with keys:
        manifold : SnapPy Manifold
        taut_isosig : str
            Isomorphism signature of the layered taut triangulation.
        face_rays : list
            Rays spanning the cone over the obvious fibered face.
            These are in the veering homology basis.
        taut_polynomial : polynomial
            The taut polynomial (= Teichmuller polynomial up to unit
            in the layered case, by Landry-Minsky-Taylor).

    Raises
    ------
    ImportError
        If snappy/regina/veering are not installed.
    RuntimeError
        If no layered taut structure is found.
    """
    if not HAS_VEERING:
        raise ImportError(
            "This function requires snappy, regina, and veering. "
            "Install via: sage -pip install snappy regina veering"
        )

    # Build the mapping torus.
    M = snappy.Manifold(braid_spec(word))

    # Bridge SnapPy -> Regina.
    T = regina.SnapPeaTriangulation(M._to_string())

    # Find layered taut angle structures.
    layered = []
    for regina_angle in regina.AngleStructures(T, True):  # tautOnly=True
        angle = taut_regina_angle_struct_to_taut_struct(regina_angle)
        if is_layered(T, angle):
            sig = isosig_from_tri_angle(T, angle)
            layered.append(sig)

    if not layered:
        raise RuntimeError(
            "No layered taut structure found. If the braid is pseudo-Anosov, "
            "try rebuilding with flipper.bundle() for Agol's veering "
            "triangulation directly."
        )

    sig = layered[0]
    rays = cone_in_homology(sig)
    theta = taut_polynomial_via_fox_calculus(sig)

    return {
        "manifold": M,
        "taut_isosig": sig,
        "face_rays": rays,
        "taut_polynomial": theta,
    }


def obvious_fibered_face_flipper(surface_word):
    """Alternative: use flipper to build the veering triangulation directly.

    Parameters
    ----------
    surface_word : str
        Mapping class in flipper notation, e.g. 'abcABC' for Dehn twists.

    Returns
    -------
    dict with taut_isosig, face_rays, taut_polynomial.

    Notes
    -----
    flipper constructs Agol's veering triangulation directly for
    pseudo-Anosov mapping classes, bypassing the need to search for
    layered structures.
    """
    try:
        import flipper
    except ImportError:
        raise ImportError("flipper not installed. Install via: pip install flipper")

    if not HAS_VEERING:
        raise ImportError("veering package not installed.")

    S = flipper.load('S_1_1')  # adjust surface as needed
    h = S.mapping_class(surface_word)
    T = h.bundle()  # returns a veering triangulation

    sig = T.isosig()
    rays = cone_in_homology(sig)
    theta = taut_polynomial_via_fox_calculus(sig)

    return {
        "taut_isosig": sig,
        "face_rays": rays,
        "taut_polynomial": theta,
    }
