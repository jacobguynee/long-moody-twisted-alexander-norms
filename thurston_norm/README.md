# Thurston vs Alexander Norm: Fibered Faces of Braid Mapping Tori

Compute the obvious fibered face of the Thurston norm ball for mapping tori of
pseudo-Anosov braids in the **meridian homology basis** -- the same basis used
for the multivariable Alexander polynomial -- and compare with the Alexander
norm.

## Goal

**When does the Alexander face containing the obvious Thurston fibered face
equal it?**

- For **3-braids**, we verify that the faces always coincide (known theorem).
- For **4-braids**, we search for examples where the Alexander face is
  **strictly larger** than the Thurston face, proving the norms differ.

## Mathematical setup

Given an n-strand braid word w, `Manifold('Braid:[...]')` builds the mapping
torus.  We recover the closed-braid + axis link via `exterior_to_link` and work
with `L.exterior()` so that the peripheral structure (meridians) is canonical.

**Alexander side**: Morton's formula gives the multivariable Alexander polynomial
as det(I - s B) where B is the colored reduced Burau matrix, with variables
identified along permutation cycles.  This is the same Burau representation
computed in the existing MATLAB code (`reduced_long_moody.m`).  The Alexander
norm ball is dual to the Newton polytope of this polynomial.

**Thurston side**: TNorm computes the Thurston norm ball of the link exterior in
the meridian-dual H_2 basis.  The obvious fibered face is the top-dimensional
fibered face containing the fiber class.

**Comparison**: McMullen's inequality says ||phi||_A <= ||phi||_T for all
phi in H^1.  A vertex v of the Thurston unit ball has ||v||_T = 1.  If
||v||_A = 1 for all vertices of the fibered face, the faces coincide.  If
||v||_A < 1 for some vertex, the Alexander face strictly contains the Thurston
face.

## Requirements

Run inside **SageMath** (>= 9.x):

```bash
sage -pip install tnorm snappy
```

TNorm pulls in Regina.  SnapPy handles the braid manifold construction.

## Usage

**Analyze a single braid:**
```bash
sage -python braid_meridian_face.py --word 1 -2 1 -2 --verbose
```

**Verify 3-braids (faces should coincide):**
```bash
sage -python braid_meridian_face.py --test3 --maxlen 8
```

**Search 4-braids for non-coincidence:**
```bash
sage -python braid_meridian_face.py --test4 --maxlen 6
sage -python braid_meridian_face.py --search4 --maxlen 8
```

**Run full test suite:**
```bash
sage -python test_alternating.py --maxlen3 8 --maxlen4 6
```

## Files

- `colored_burau.py` -- Colored reduced Burau matrix, Morton's Alexander
  polynomial, Newton polytope, Alexander norm computation
- `braid_meridian_face.py` -- Main pipeline: braid word -> Thurston face ->
  Alexander comparison.  Includes test harnesses for 3- and 4-braids.
- `test_alternating.py` -- Verification script: confirms 3-braid coincidence,
  searches for 4-braid counterexamples
- `README.md` -- This file

## Connection to existing code

The `colored_burau.py` module is the Python/Sage translation of the MATLAB
`reduced_long_moody.m` specialized to the Burau case (rho(g_i) = t_i,
rho(sigma_i) = 1).  Morton's formula says det(I - B) is the Alexander
polynomial in meridian variables, which is exactly the basis needed for
comparison with the Thurston norm computed by TNorm.
