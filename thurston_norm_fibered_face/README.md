# Thurston Norm Fibered Face Computation

## Overview

Compute the obvious fibered face of the Thurston norm unit ball for
pseudo-Anosov braid closures, using layered veering triangulations
(Landry-Minsky-Taylor).

## Files

| File | Description |
|------|-------------|
| `compute_fibered_face.py` | **Main script.** Computes the fibered face via SnapPy + Regina + veering. Requires SageMath environment. |
| `alexander_thurston.py` | Computes the Alexander polynomial and Alexander norm from the reduced Burau representation. Pure SymPy, runs anywhere. Used to validate the veering computation. |
| `test_alternating_braids.py` | Computes Alexander norms for 35 alternating braids in B_3 and B_4. These are ground-truth values for validating `compute_fibered_face.py`. |

## How It Works

### Thurston norm computation (`compute_fibered_face.py`)

For a pseudo-Anosov n-braid beta:

1. **SnapPy** builds the mapping torus M_beta from the braid word.
2. **Regina** enumerates taut angle structures on the triangulation.
3. **veering** identifies layered taut structures and computes:
   - The cone over the fibered face via `cone_in_homology()`
   - The taut polynomial via `taut_polynomial_via_fox_calculus()`

Landry-Minsky-Taylor proved that the cone carried by a layered veering
triangulation is exactly the cone over the corresponding fibered face of
the Thurston norm ball.

### Validation strategy

For alternating braids in B_3 and B_4 that use all generators (hence
pseudo-Anosov), the Alexander norm equals the Thurston norm. So:

1. Run `test_alternating_braids.py` to get Alexander norms (pure SymPy).
2. Run `compute_fibered_face.py` on the same braids (SageMath + veering).
3. Verify: Thurston norm from veering == Alexander norm from Burau.

## Installation

`compute_fibered_face.py` requires a SageMath environment with topology
packages. Here's how to set it up:

### Option 1: SageMath + pip (recommended)

```bash
# Install SageMath
# Ubuntu/Debian:
sudo apt install sagemath
# Or via conda:
conda create -n sage sage -c conda-forge && conda activate sage

# Install topology packages inside Sage's Python
sage -pip install snappy
sage -pip install regina
sage -pip install veering
sage -pip install flipper   # optional, for direct veering triangulations
```

### Option 2: Docker

```bash
docker run -it computop/sage bash
sage -pip install veering flipper
```

## Usage

### Thurston norm (requires SageMath + veering)

```bash
# Run on default examples
sage -python compute_fibered_face.py

# Run on a specific braid
sage -python compute_fibered_face.py 1 -2 3
```

```python
from compute_fibered_face import obvious_fibered_face

data = obvious_fibered_face([1, -2, 3])
print(data['taut_isosig'])       # isosig of layered veering triangulation
print(data['face_rays'])         # rays of cone over fibered face
print(data['taut_polynomial'])   # = Teichmüller polynomial (layered case)
```

### Alexander norm validation (runs anywhere with SymPy)

```bash
python test_alternating_braids.py
```

```python
from alexander_thurston import compute_alexander_data

data = compute_alexander_data([1, -2, 3], n=4)
print(data['alexander_norm'])    # should match Thurston norm for alternating braids
```

## Caveats

- **Homology basis**: `face_rays` from veering are in the veering homology
  basis, dual to the basis used for taut/veering polynomials. You may need
  a change-of-basis to convert to geometric meridian/longitude coordinates.

- **Layered structures**: The code searches for layered taut angle structures
  on the SnapPy triangulation. If none are found (rare for pseudo-Anosov
  braids), use `flipper.bundle()` to construct Agol's veering triangulation
  directly — see `obvious_fibered_face_flipper()`.

## References

- Landry, M., Minsky, Y., & Taylor, S. (2021). A polynomial invariant
  for veering triangulations. *Journal of Topology*.
- Agol, I. (2011). Ideal triangulations of pseudo-Anosov mapping tori.
  *Topology and Geometry in Dimension Three*, AMS.
- McMullen, C. T. (2002). The Alexander polynomial of a 3-manifold and
  the Thurston norm on cohomology. *Annales scientifiques de l'ENS*.
