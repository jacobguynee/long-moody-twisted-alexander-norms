"""
representations.py

Core representation constructions for B_n:
  - Reduced Long-Moody (colored/mixed braids)
  - Free group generator recovery
  - Reduced Burau as Long-Moody special case
  - Lawrence k=2 as iterated Long-Moody

Exact port of the MATLAB code in the repo root (reduced_lm_colored.m,
free_group_gens.m). See the repo root files for the original MATLAB
implementations and detailed comments.
"""

from sympy import symbols, eye, zeros, Matrix, cancel, Rational

q, x = symbols('q x')


# ---- Block matrix helpers ----

def build_Ni(i, gi, n, d):
    """Build the (n-1)x(n-1) block matrix N_i for generator sigma_i.
    gi is a d x d matrix (image of the free group generator g_i)."""
    nb = n - 1
    Id = eye(d); Zd = zeros(d)
    blocks = [[Id.copy() if r == c else Zd.copy() for c in range(nb)] for r in range(nb)]
    if i == 0:
        blocks[0][0] = -gi; blocks[0][1] = Zd
        blocks[1][0] = Id;  blocks[1][1] = Id
    elif i == n - 2:
        blocks[nb-2][nb-2] = Id;  blocks[nb-2][nb-1] = gi
        blocks[nb-1][nb-2] = Zd;  blocks[nb-1][nb-1] = -gi
    else:
        p = i - 1
        blocks[p][p] = Id;    blocks[p][p+1] = gi;    blocks[p][p+2] = Zd
        blocks[p+1][p] = Zd;  blocks[p+1][p+1] = -gi; blocks[p+1][p+2] = Zd
        blocks[p+2][p] = Zd;  blocks[p+2][p+1] = Id;  blocks[p+2][p+2] = Id
    return blocks


def left_multiply_blocks(blocks, M, nb):
    return [[M * blocks[r][c] for c in range(nb)] for r in range(nb)]


def multiply_block_matrices(A, B, nb, d):
    C = [[zeros(d) for _ in range(nb)] for _ in range(nb)]
    for r in range(nb):
        for c in range(nb):
            acc = zeros(d)
            for k in range(nb):
                acc = acc + A[r][k] * B[k][c]
            C[r][c] = acc
    return C


def assemble(blocks, nb, d):
    N = nb * d
    M = zeros(N)
    for r in range(nb):
        for c in range(nb):
            M[r*d:(r+1)*d, c*d:(c+1)*d] = blocks[r][c]
    return M


# ---- Core constructions ----

def reduced_lm_colored(rho_g, rho_sigma, partition=None):
    """Reduced Long-Moody construction for colored/mixed braids.

    Inputs:
        rho_g:     list of n  d x d matrices [rho(g_1), ..., rho(g_n)]
        rho_sigma: list of n-1 d x d matrices [rho(sigma_1), ..., rho(sigma_{n-1})]
        partition: length-n list assigning each strand to a part, e.g. [1,1,1,2].
                   At boundary i (partition[i] != partition[i+1]), output is
                   rho^+(sigma_i^2) instead of rho^+(sigma_i).

    Output:
        list of n-1 matrices, each (n-1)*d x (n-1)*d.
    """
    n = len(rho_g); num_sigmas = len(rho_sigma); d = rho_g[0].shape[0]
    if partition is None:
        partition = [1] * n
    is_boundary = [partition[i] != partition[i+1] for i in range(num_sigmas)]
    nb = n - 1
    rho_plus = []
    for i in range(num_sigmas):
        if not is_boundary[i]:
            blocks = build_Ni(i, rho_g[i], n, d)
            blocks = left_multiply_blocks(blocks, rho_sigma[i], nb)
            rho_plus.append(cancel(assemble(blocks, nb, d)))
        else:
            gi_conj = cancel(rho_g[i] * rho_g[i+1] * rho_g[i].inv())
            blocks_N  = build_Ni(i, rho_g[i], n, d)
            blocks_Np = build_Ni(i, gi_conj,  n, d)
            blocks_prod = multiply_block_matrices(blocks_Np, blocks_N, nb, d)
            si_sq = rho_sigma[i] * rho_sigma[i]
            blocks_prod = left_multiply_blocks(blocks_prod, si_sq, nb)
            rho_plus.append(cancel(assemble(blocks_prod, nb, d)))
    return rho_plus


def free_group_gens(rho_sigma, partition=None):
    """Recover free group generator images from braid generator images.

    Uses the Birman exact sequence:
        g_n = sigma_{n-1}^2  (or sigma_{n-1} at a partition boundary)
        g_i = sigma_i^{-1} g_{i+1} sigma_i

    Input:  list of n braid generator matrices (n = len(rho_sigma))
    Output: list of n free group generator matrices [rho(g_1), ..., rho(g_n)]
    """
    n = len(rho_sigma)
    if partition is None:
        partition = [1] * (n + 1)
    last_is_boundary = (partition[n-1] != partition[n])
    gn = rho_sigma[n-1] if last_is_boundary else rho_sigma[n-1] * rho_sigma[n-1]
    rho_g = [None] * n
    rho_g[n-1] = cancel(gn)
    for i in range(n-2, -1, -1):
        rho_g[i] = cancel(rho_sigma[i].inv() * rho_g[i+1] * rho_sigma[i])
    return rho_g


# ---- B_4 representations ----

def build_burau_B4():
    """Reduced Burau representation of B_4.

    This is the Long-Moody construction with:
        rho(g_i) = q  (scalar, for all i = 1..4)
        rho(sigma_j) = 1  (scalar, for all j = 1..3)

    Output: 3 matrices of size 3x3, entries in Z[q, q^{-1}].
    """
    return reduced_lm_colored([Matrix([[q]])] * 4, [Matrix([[1]])] * 3)


def build_lawrence_B4(Q_val=Rational(2), a_val=Rational(3)):
    """Lawrence k=2 representation of B_4 via iterated Long-Moody.

    Step 1: Apply reduced Long-Moody to the mixed braid group B_{4,1}
            with partition [1,1,1,1,2], free group labels [Q,Q,Q,Q,a],
            and trivial braid images. Produces 4x4 matrices for 4 mixed
            braid generators.

    Step 2: Recover free group generators g_1,...,g_4 from the mixed braid
            generators via the Birman exact sequence.

    Step 3: Apply reduced Long-Moody again with:
            rho(g_i) = q * g_i  (4x4 matrices scaled by q)
            rho(sigma_j) = Step 1 output for j = 1..3
            This produces 12x12 matrices.

    Parameters Q_val, a_val are generic numerical substitutions to speed
    up computation. Default (2, 3) verified generic across 5 different
    (Q, a) pairs.

    Output: 3 matrices of size 12x12, entries in Q(q).
    """
    # Step 1
    L4Pre = reduced_lm_colored(
        [Matrix([[Q_val]])] * 4 + [Matrix([[a_val]])],
        [Matrix([[1]])] * 4,
        [1, 1, 1, 1, 2]
    )
    # Step 2
    L4PreFree = free_group_gens(L4Pre, [1, 1, 1, 1, 2])
    # Step 3
    return reduced_lm_colored([g * q for g in L4PreFree], L4Pre[0:3])
