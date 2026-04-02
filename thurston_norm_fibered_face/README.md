# Thurston Norm Fibered Face Computation

Compute the obvious fibered face of the Thurston norm unit ball for
pseudo-Anosov braid closures via layered veering triangulations
(Landry-Minsky-Taylor).

## Requirements

Requires SageMath with topology packages:

```bash
sage -pip install snappy
sage -pip install regina
sage -pip install veering
sage -pip install flipper   # optional
```

## Files

- `compute_fibered_face.py` — Computes the fibered face via SnapPy + Regina + veering.
- `alexander_thurston.py` — Computes Alexander norm from the reduced Burau representation (pure SymPy).
- `test_alexander_equals_thurston.py` — Verifies Alexander norm = Thurston norm on alternating braids in B_3 and B_4.

## Usage

```bash
# Compute fibered face for a braid
sage -python compute_fibered_face.py 1 -2 3

# Run default examples
sage -python compute_fibered_face.py

# Run validation tests (Alexander norm vs Thurston norm)
sage -python test_alexander_equals_thurston.py
```

```python
from compute_fibered_face import obvious_fibered_face

data = obvious_fibered_face([1, -2, 3])
print(data['taut_isosig'])
print(data['face_rays'])
print(data['taut_polynomial'])
```

## How it works

1. SnapPy builds the mapping torus from the braid word.
2. Regina enumerates taut angle structures.
3. The veering package identifies layered structures and computes the
   cone over the fibered face via `cone_in_homology()` and the taut
   polynomial via `taut_polynomial_via_fox_calculus()`.

The face rays are in the veering homology basis (dual to the polynomial
basis). A change-of-basis may be needed for meridian/longitude coordinates.

## References

- Landry, Minsky, Taylor (2021). A polynomial invariant for veering triangulations.
- Agol (2011). Ideal triangulations of pseudo-Anosov mapping tori.
- McMullen (2002). The Alexander polynomial of a 3-manifold and the Thurston norm on cohomology.
