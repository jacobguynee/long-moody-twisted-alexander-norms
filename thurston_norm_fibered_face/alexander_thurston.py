"""
alexander_thurston.py

Compute the Alexander polynomial and Alexander norm for braid closures,
using the reduced Burau representation (via Long-Moody construction).

This does NOT compute the Thurston norm. For the actual Thurston norm
computation via veering triangulations, use compute_fibered_face.py.

The Alexander polynomial of an n-braid closure is:
    Delta(q) = det(I - Burau(beta))
where Burau(beta) is the reduced Burau matrix in Z[q, q^{-1}].
"""

from sympy import symbols, eye, zeros, Matrix, cancel, Poly, Rational, fraction

q = symbols('q')


# ---- Reduced Long-Moody / Burau ----

def build_Ni(i, gi, n, d):
    nb = n - 1; Id = eye(d); Zd = zeros(d)
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


def assemble(blocks, nb, d):
    N = nb * d; M = zeros(N)
    for r in range(nb):
        for c in range(nb):
            M[r*d:(r+1)*d, c*d:(c+1)*d] = blocks[r][c]
    return M


def reduced_lm_colored(rho_g, rho_sigma, partition=None):
    n = len(rho_g); num_sigmas = len(rho_sigma); d = rho_g[0].shape[0]
    if partition is None:
        partition = [1] * n
    is_boundary = [partition[i] != partition[i+1] for i in range(num_sigmas)]
    nb = n - 1; rho_plus = []
    for i in range(num_sigmas):
        if not is_boundary[i]:
            blocks = build_Ni(i, rho_g[i], n, d)
            blocks = left_multiply_blocks(blocks, rho_sigma[i], nb)
            rho_plus.append(cancel(assemble(blocks, nb, d)))
        else:
            gi_conj = cancel(rho_g[i] * rho_g[i+1] * rho_g[i].inv())
            blocks_N  = build_Ni(i, rho_g[i], n, d)
            blocks_Np = build_Ni(i, gi_conj,  n, d)
            from sympy import zeros as _z
            def _mul(A, B, nb, d):
                C = [[_z(d) for _ in range(nb)] for _ in range(nb)]
                for r in range(nb):
                    for c in range(nb):
                        acc = _z(d)
                        for k in range(nb): acc = acc + A[r][k] * B[k][c]
                        C[r][c] = acc
                return C
            blocks_prod = _mul(blocks_Np, blocks_N, nb, d)
            si_sq = rho_sigma[i] * rho_sigma[i]
            blocks_prod = left_multiply_blocks(blocks_prod, si_sq, nb)
            rho_plus.append(cancel(assemble(blocks_prod, nb, d)))
    return rho_plus


def build_burau(n):
    """Reduced Burau generators for B_n. Returns list of (n-1)x(n-1) matrices."""
    return reduced_lm_colored([Matrix([[q]])] * n, [Matrix([[1]])] * (n - 1))


def eval_word(word, gens, inv_gens):
    M = eye(gens[0].shape[0])
    for w in word:
        if w > 0: M = M * gens[w - 1]
        else: M = M * inv_gens[-w - 1]
    return cancel(M)


# ---- Alexander polynomial ----

def alexander_polynomial(word, n):
    """Alexander polynomial of the closure of an n-braid: det(I - Burau(beta)).

    Returns a Laurent polynomial in q (up to unit q^k).
    """
    bg = build_burau(n)
    bgi = [cancel(g.inv()) for g in bg]
    M = eval_word(word, bg, bgi)
    dim = M.shape[0]
    return cancel((eye(dim) - M).det())


def newton_polytope_1d(poly_expr, var=q):
    """Newton polytope of a Laurent polynomial in one variable.

    Returns (min_degree, max_degree).
    """
    num, den = fraction(cancel(poly_expr))
    p_num = Poly(num, var)
    q_shift = Poly(den, var).degree() if den.has(var) else 0
    degs = [d[0] - q_shift for d in p_num.monoms()]
    return (min(degs), max(degs))


def alexander_norm(word, n):
    """Alexander norm of a braid closure = breadth of Alexander polynomial."""
    delta = alexander_polynomial(word, n)
    lo, hi = newton_polytope_1d(delta)
    return hi - lo
