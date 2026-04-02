# Newton Polytope Conjecture: Burau vs Lawrence k=2 for B_4

## Conjecture

For any braid `beta` in `B_n`, the Newton polytopes of the characteristic
polynomials of the reduced Burau representation and the Lawrence k=2
representation are related by:

```
NP(charpoly(Lawrence_k=2(beta))) / n  =  NP(charpoly(Burau(beta)))
```

where the Newton polytope is computed in (eigenvalue x, parameter q)
coordinates, and `/n` denotes scaling both coordinates by `1/n`.

For B_4 specifically: `NP(Lawrence) / 4 = NP(Burau)`, where Burau is
3x3 and Lawrence is 12x12.

## Background

### Representations

**Reduced Burau** for B_n is an (n-1)-dimensional representation with
entries in Z[q, q^{-1}]. It arises as the reduced Long-Moody construction
with trivial inputs:

- Free group labels: rho(g_i) = q for all i
- Braid images: rho(sigma_j) = 1 for all j

**Lawrence k=2** is an n(n-1)-dimensional representation acting on H_1
of the ordered 2-point configuration space in the n-punctured disk
(Lawrence 1990, Bigelow 2001). It can be computed as an *iterated*
reduced Long-Moody construction (Bigelow-Tian 2008, Section 6):

1. Start with the mixed braid group B_{n,1} (partition [1,...,1,2]).
   Apply reduced Long-Moody with labels [Q,...,Q,a] and trivial braid
   images, producing (n-1)x(n-1) matrices for n mixed braid generators.

2. Recover free group generators g_1,...,g_n via the Birman exact
   sequence: g_n = sigma_{n-1} (at partition boundary), then
   g_i = sigma_i^{-1} g_{i+1} sigma_i.

3. Apply reduced Long-Moody again with rho(g_i) = q * g_i (Step 2
   output scaled by q) and rho(sigma_j) = Step 1 output. Produces
   (n-1)^2 x (n-1)^2 = n(n-1) x n(n-1) matrices (after simplification
   with the partition structure).

### Newton polytope

The **characteristic polynomial** of an (n-1)x(n-1) matrix with entries
in Z[q, q^{-1}] is a polynomial in eigenvalue x with Laurent polynomial
coefficients in q. The **Newton polytope** is the convex hull of all
(x-degree, q-degree) pairs with nonzero coefficients.

### Generic substitution

To speed up symbolic computation of the 12x12 Lawrence charpoly, we
substitute generic numerical values Q=2, a=3 (verified to be generic
across 5 different (Q,a) pairs — see `../verify_generic_substitution.py`
in the repo root). This reduces entries from rational functions in
(Q, q, a) to rational functions in q alone.

## Evidence

Tested on 700+ braids in B_4 with **zero counterexamples**:

| Category | Count | Result |
|----------|-------|--------|
| Alternating braids (monoid s1, s2^-1, s3) | 26 | all match |
| f^a g^b, 16 almost-reducible families | 150+ | all match |
| f^a g^b, 20 additional diverse families | 71+ | all match |
| Reducible braids | 49 | all match |
| Zero-writhe, conjugates, random | 47 | all match |
| General words f^{a1}g^{b1}...f^{ak}g^{bk} | 30+ | all match |
| Random braid words | 20+ | all match |
| Positive braids (word length 10-20) | 20 | all match |
| **B_5 (scaling 1/5, 20x20 Lawrence)** | 24 | all match |

### Structural observations

1. **Vertex count bounded**: Newton polytope hull has at most 6 vertices
   for B_4 (theoretical bound 2(n+1) = 10), at most 4 for B_3.

2. **Monomial containment**: The Lawrence charpoly always has *extra*
   monomials at fractional x-coordinates (1/4, 1/2, 3/4, ...) compared
   to Burau. These always lie strictly inside the convex hull. Burau
   never has monomials that Lawrence (scaled) doesn't.

3. **No exterior cancellations observed**: No braid found where a hull
   vertex monomial has coefficient zero in one representation but not
   the other.

## Files

| File | Description |
|------|-------------|
| `representations.py` | Core: reduced Long-Moody, free group gens, Burau, Lawrence builders |
| `newton_polytope.py` | Newton polytope computation, convex hull, comparison utilities |
| `verify_conjecture.py` | Main test script: runs diverse braids and checks the conjecture |
| `verify_conjugacy.py` | Verifies iterated Long-Moody matches Lawrence (1990) for B_4 |

## Usage

```bash
# Verify the conjecture on diverse braids (~5-10 min)
python verify_conjecture.py

# Verify that iterated Long-Moody = Lawrence for B_4 (~2 min, fully symbolic)
python verify_conjugacy.py
```

Requires Python 3 with SymPy.

## Relationship to repo root files

The MATLAB originals of the representation constructions are in the repo
root: `reduced_lm_colored.m`, `free_group_gens.m`, `lawrence_generators.m`.

Additional test scripts in the repo root cover:
- B_3 tests: `compare_newton_polytopes.py`
- B_4 conjugated braids: `compare_newton_polytopes_B4_conjugated.py`,
  `compare_newton_polytopes_B4_high.py`
- Full twist family: `compare_newton_polytopes_B4_fulltwist.py`
- Positive braids: `compare_newton_polytopes_B4_positive.py`
- Polytope side analysis: `analyze_polytope_sides.py`
- Counterexample hunts: `hunt_counterexample_B4*.py`, `hunt_exterior_cancellation.py`
- B_5 extension: `test_B5_newton_polytopes.py`
- Generic substitution verification: `verify_generic_substitution.py`
- Monomial-level comparison: `hunt_monomial_mismatch.py`

## References

- Lawrence, R. J. (1990). Homological representations of the Hecke algebra.
  *Communications in Mathematical Physics*, 135(1), 141-191.
- Bigelow, S. (2001). Braid groups are linear. *Journal of the AMS*, 14(2), 471-486.
- Bigelow, S. & Tian, Y. (2008). Generalized Long-Moody representations of
  the braid group. Preprint.
