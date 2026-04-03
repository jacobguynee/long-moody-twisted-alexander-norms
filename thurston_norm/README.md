# Thurston vs Alexander Norm: Fibered Faces of Braid Mapping Tori

## What this does

Given a braid word in Artin generators (e.g. `[1, -2, 1, -2]`), this code:

1. Computes the **obvious fibered face** of the Thurston norm ball for the
   braid's mapping torus.
2. Computes the **Alexander polynomial** via Morton's colored Burau formula.
3. Checks whether the **Alexander face containing the Thurston fibered face
   equals it**, or is strictly larger.

Both norms are expressed in the **meridian homology basis**, so they are
directly comparable.

## Quick start

### Prerequisites

Everything runs inside **SageMath** (>= 9.x).  Install dependencies:

```bash
sage -pip install tnorm snappy
```

TNorm pulls in Regina automatically.

### First thing to check: basis alignment

Before trusting any comparison results, **verify that the Alexander and
Thurston sides are using the same basis** by running the 3-braid tests:

```bash
cd thurston_norm/
sage -python braid_meridian_face.py --test3 --maxlen 8
```

For 3-braids, the Alexander face is **known** to equal the Thurston fibered
face.  If all tests print PASS, the bases are aligned and the pipeline is
working.  If any test prints FAIL, something is wrong with the basis
conversion -- do not trust the 4-braid results until this is resolved.

### Then: search for 4-braid counterexamples

Once 3-braids pass, search for 4-braids where the Alexander face is strictly
larger:

```bash
sage -python braid_meridian_face.py --test4 --maxlen 6
```

Or run the full test suite (3-braid verification + 4-braid search):

```bash
sage -python test_alternating.py --maxlen3 8 --maxlen4 6
```

### Analyze a single braid

```bash
sage -python braid_meridian_face.py --word 1 -2 1 -2 --verbose
```

This prints the Alexander polynomial, Thurston face vertices, and whether the
faces coincide.

### Exhaustive 4-braid search

If the curated 4-braid list doesn't find counterexamples, try brute force:

```bash
sage -python braid_meridian_face.py --search4 --maxlen 6
```

Warning: this is slow (exponential in word length).

## What the output means

For each braid, the code computes the Thurston fibered face vertices v_1, ...,
v_k.  Each vertex lies on the Thurston unit ball, so ||v_i||_T = 1.  It then
computes ||v_i||_A (the Alexander norm) at each vertex.

- **||v_i||_A = 1 for all i**: The Alexander and Thurston faces **coincide**.
  McMullen's inequality ||.||_A <= ||.||_T is tight on this face.

- **||v_i||_A < 1 for some i**: The Alexander face **strictly contains** the
  Thurston face.  This proves the Alexander and Thurston norms differ for
  this manifold.

## Mathematical background

### Setup

An n-strand braid word w determines a mapping torus.  SnapPy's
`Manifold('Braid:[...]')` builds this as the complement of the closed braid
union its axis in S^3.  This is a link complement with one component per
permutation cycle of the braid, plus one for the axis.

We recover a link diagram via `exterior_to_link(check_answer=True)` and then
use `L.exterior()` so that SnapPy's peripheral structure (meridians) is locked
to the link components.  This is essential: it ensures TNorm's Thurston norm
coordinates and the Alexander polynomial exponents live in the **same basis**.

### Alexander side

Morton's formula: the multivariable Alexander polynomial of the closed braid +
axis is det(I - s B) where B is the colored reduced Burau matrix.  Strand
variables t_0, ..., t_{n-1} are identified along permutation cycles to give
meridian variables m_0, ..., m_{k-1} for the k closure components, and s is
the axis meridian.  This is the Python/Sage translation of the MATLAB
`reduced_long_moody.m` specialized to the Burau case (rho(g_i) = t_i,
rho(sigma_i) = 1).

The Alexander norm of a cohomology class phi is the width of the Newton
polytope of Delta in direction phi.

### Thurston side

TNorm computes the full Thurston norm ball of the link exterior.  When the
input manifold comes from `L.exterior()`, TNorm works in the meridian-dual
H_2 basis -- the same coordinate system as the Alexander polynomial exponents.

The obvious fibered face is the top-dimensional fibered face of the Thurston
norm ball containing the fiber class of the open book decomposition.

### Comparison criterion

McMullen proved ||phi||_A <= ||phi||_T for all phi.  Equality at every vertex
of a Thurston face means the Alexander face (the face of the Alexander norm
ball in the same direction) coincides with the Thurston face.

For 3-braids, equality always holds on fibered faces.  For 4-braids, there may
exist examples where the inequality is strict, meaning the Alexander norm
underestimates the Thurston norm and the Alexander face is too large.

## Files

| File | What it does |
|------|-------------|
| `colored_burau.py` | Colored reduced Burau matrix, Morton's Alexander polynomial, Newton polytope, Alexander norm computation |
| `braid_meridian_face.py` | Main pipeline: braid word -> Thurston face -> Alexander comparison. Includes `--test3`, `--test4`, `--search4` modes |
| `test_alternating.py` | Combined test harness: 3-braid verification + 4-braid counterexample search |
| `README.md` | This file |

## Potential issues

- **TNorm fails on some manifolds**: If the triangulation is too complex or
  non-ideal, TNorm may error.  These braids are skipped and reported.

- **`exterior_to_link` fails**: Some mapping tori may not be recognizable as
  link complements by SnapPy's heuristics.  This is rare for small braids.

- **Basis mismatch**: The entire comparison hinges on both sides using the
  meridian basis.  The 3-braid sanity check (`--test3`) exists precisely to
  catch this.  If 3-braids fail, the comparison is meaningless.

- **Pseudo-Anosov assumption**: The pipeline assumes the braid is
  pseudo-Anosov (so the mapping torus is hyperbolic and has a well-defined
  fibered face).  Non-pseudo-Anosov braids (reducible or periodic) will
  likely produce errors or meaningless results.
