# Thurston Norm Fibered Face Computation

## Overview

Compute the obvious fibered face of the Thurston norm unit ball for braid
closures, using the fact that for alternating braids in B_3 and B_4 (using
all generators), Alexander norm = Thurston norm.

Two approaches are provided:

1. **`alexander_thurston.py`** — Pure SymPy. Computes the Alexander polynomial
   as det(I - Burau(beta)), then extracts the fibered face from its Newton
   polytope. Works in any Python environment with SymPy.

2. **`veering_approach.py`** — Reference implementation using SnapPy + Regina +
   veering package. Computes the fibered face via layered veering triangulations
   (Landry-Minsky-Taylor). Requires a SageMath environment with topology packages.

## Quick Start (alexander_thurston.py)

This works out of the box with just Python 3 + SymPy:

```bash
# Run the test suite (35 alternating braids in B_3 and B_4)
python test_alternating_braids.py

# Use interactively
python -c "
from alexander_thurston import fibered_face_from_alexander, print_fibered_face
data = fibered_face_from_alexander([1, -2, 3], 4)
print_fibered_face(data)
"
```

### What `fibered_face_from_alexander` returns

```python
{
    'braid_word': [1, -2, 3],
    'n_strands': 4,
    'writhe': 1,
    'burau_matrix': <SymPy Matrix>,        # Reduced Burau of braid word
    'alexander_poly': <SymPy expr>,         # det(I - Burau(beta)) in Z[q,q^-1]
    'alexander_poly_newton_1d': (-1, 2),    # (min_deg, max_deg) of Alexander poly
    'alexander_norm': 3,                    # breadth = max_deg - min_deg
    'genus': 3/2,                           # = alexander_norm / 2
    'charpoly_newton_2d': [(0,1), ...],     # convex hull of charpoly in (x, q)
    'thurston_norm_of_fiber': 3,            # = alexander norm (for alternating braids)
}
```

## Setting Up the Veering Approach (veering_approach.py)

The veering approach requires a SageMath environment. Here's how to set it up:

### Option 1: SageMath + pip (recommended)

```bash
# Install SageMath (Ubuntu/Debian)
sudo apt install sagemath

# Or via conda
conda create -n sage sage -c conda-forge
conda activate sage

# Install topology packages inside Sage's Python
sage -pip install snappy
sage -pip install regina
sage -pip install veering
sage -pip install flipper  # optional, for direct veering triangulations
```

### Option 2: Docker

```bash
# Use the computop/sage image which has snappy pre-installed
docker run -it computop/sage bash
sage -pip install veering flipper
```

### Using the veering approach

```python
# Inside a SageMath session or sage -python
from veering_approach import obvious_fibered_face_veering

# Braid word as signed generator indices
result = obvious_fibered_face_veering([1, -2, 3])
print(result['taut_isosig'])
print(result['face_rays'])       # Rays of the cone over the fibered face
print(result['taut_polynomial']) # = Teichmuller polynomial (layered case)

# Alternative: use flipper for pseudo-Anosov mapping classes
from veering_approach import obvious_fibered_face_flipper
result = obvious_fibered_face_flipper('abcABC')
```

### What the veering approach computes

1. Builds the mapping torus of the braid as a SnapPy manifold.
2. Converts to a Regina triangulation.
3. Searches for layered taut angle structures (these exist when the
   braid is pseudo-Anosov).
4. The layered veering triangulation carries a cone that is exactly
   the cone over the fibered face (Landry-Minsky-Taylor theorem).
5. Returns the face rays and the taut polynomial (= Teichmuller
   polynomial up to a unit in the layered case).

### Verifying Alexander norm = Thurston norm

Once the veering packages are installed, you can cross-check:

```python
from alexander_thurston import fibered_face_from_alexander
from veering_approach import obvious_fibered_face_veering

word = [1, -2, 3, 1, -2, 3]  # alternating B_4 braid

# Alexander approach
alex = fibered_face_from_alexander(word, 4)
print(f"Alexander norm: {alex['alexander_norm']}")

# Veering approach
veer = obvious_fibered_face_veering(word)
print(f"Taut polynomial: {veer['taut_polynomial']}")
print(f"Face rays: {veer['face_rays']}")
```

## Mathematical Background

### Alexander polynomial from Burau

For an n-braid beta, the Alexander polynomial of its closure is:

```
Delta(t) = det(I - Burau(beta))
```

where Burau(beta) is the reduced Burau representation evaluated at
q = t. The reduced Burau is computed via the Long-Moody construction
with trivial inputs: rho(g_i) = q, rho(sigma_j) = 1.

### Alexander norm and Thurston norm

The **Alexander norm** of a class in H^1(M; Z) is defined via the
Newton polytope of the Alexander polynomial. For a knot (1-component
braid closure), it equals the breadth (max_deg - min_deg) of Delta(t).

The **Thurston norm** measures the minimal complexity of embedded
surfaces representing a homology class. McMullen (2002) showed:

```
Alexander norm <= Thurston norm
```

For alternating braids in B_3 and B_4 that use all generators (hence
pseudo-Anosov), equality holds (author's theorem):

```
Alexander norm = Thurston norm  (alternating braids in B_3, B_4)
```

### Fibered face

The mapping torus of a pseudo-Anosov braid is fibered. The **fibered
face** is the top-dimensional face of the Thurston norm unit ball dual
to the fiber class. In the knot case (H^1 = R), the unit ball is
[-1/c, 1/c] where c = Thurston norm of the fiber, and the fibered
face consists of the two endpoints.

For links (multi-component closures), the fibered face lives in
H^1 = R^k and is determined by the Newton polytope of the
multivariable Alexander polynomial.

## Files

| File | Description |
|------|-------------|
| `alexander_thurston.py` | Pure-SymPy fibered face computation from Alexander polynomial |
| `veering_approach.py` | Reference implementation using SnapPy/Regina/veering (requires SageMath) |
| `test_alternating_braids.py` | Test suite: 35 alternating braids in B_3 and B_4 |

## Dependencies

- **alexander_thurston.py**: Python 3, SymPy (works anywhere)
- **veering_approach.py**: SageMath, SnapPy, Regina, veering, (optionally flipper)
- **test_alternating_braids.py**: Python 3, SymPy

## References

- McMullen, C. T. (2002). The Alexander polynomial of a 3-manifold and
  the Thurston norm on cohomology. *Annales scientifiques de l'ENS*.
- Landry, M., Minsky, Y., & Taylor, S. (2021). A polynomial invariant
  for veering triangulations. *Journal of Topology*.
- Agol, I. (2011). Ideal triangulations of pseudo-Anosov mapping tori.
  *Topology and Geometry in Dimension Three*, AMS.
- Lawrence, R. J. (1990). Homological representations of the Hecke algebra.
