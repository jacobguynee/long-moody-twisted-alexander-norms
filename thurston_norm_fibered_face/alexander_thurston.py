"""
alexander_thurston.py

Compute the Alexander polynomial and Alexander norm for braid closures,
using the reduced Burau representation (via Long-Moody construction).

This does NOT compute the Thurston norm. For the actual Thurston norm
computation via veering triangulations, use compute_fibered_face.py.

Purpose: provide Alexander norm values to VALIDATE the veering computation.
For alternating braids in B_3 and B_4 (using all generators), the Alexander
norm equals the Thurston norm. So we can verify that compute_fibered_face.py
gives the correct Thurston norm by checking agreement with this module on
alternating braids.

The Alexander polynomial of an n-braid closure is:
    Delta(q) = det(I - Burau(beta))
where Burau(beta) is the reduced Burau matrix in Z[q, q^{-1}].
"""

from sympy import symbols, eye, zeros, Matrix, cancel, Poly, Rational, fraction, degree

q, x = symbols('q x')
t = symbols('t')


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

def alexander_polynomial_1var(word, n):
    """1-variable Alexander polynomial: det(I - Burau(beta)).

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


def alexander_norm_1d(poly_expr, var=q):
    """Alexander norm = breadth of the Alexander polynomial.

    For a knot, this is max_deg - min_deg of Delta(t).
    """
    lo, hi = newton_polytope_1d(poly_expr, var)
    return hi - lo


def newton_polytope_2d(word, n):
    """Newton polytope of charpoly of Burau(beta) in (x, q) coordinates.

    Returns set of (x_deg, q_deg) monomials.
    """
    bg = build_burau(n)
    bgi = [cancel(g.inv()) for g in bg]
    M = eval_word(word, bg, bgi)
    cp_expr = cancel(M.charpoly(x).as_expr())
    num, den = fraction(cp_expr)
    p = Poly(num, x, q)
    monoms_num = p.monoms()
    q_shift = Poly(den, q).degree() if den.has(q) else 0
    return set((xd, qd - q_shift) for (xd, qd) in monoms_num)


def convex_hull_2d(points):
    points = sorted(set(points))
    if len(points) <= 2: return points
    def cross(O, A, B):
        return (A[0]-O[0])*(B[1]-O[1]) - (A[1]-O[1])*(B[0]-O[0])
    lower = []
    for p in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0: lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0: upper.pop()
        upper.append(p)
    return sorted(set(lower[:-1] + upper[:-1]))


# ---- Main computation ----

def compute_alexander_data(word, n):
    """Compute Alexander norm data for a braid closure.

    Returns a dict with Alexander polynomial, Newton polytope,
    Alexander norm, and genus estimate.

    To validate veering Thurston norm computations, compare
    alexander_norm from this function against the Thurston norm
    from compute_fibered_face.obvious_fibered_face() on
    alternating braids (where equality is known to hold).
    """
    bg = build_burau(n)
    bgi = [cancel(g.inv()) for g in bg]
    M = eval_word(word, bg, bgi)

    # 1-variable Alexander polynomial
    dim = M.shape[0]
    delta = cancel((eye(dim) - M).det())
    np_1d = newton_polytope_1d(delta)
    alex_norm = np_1d[1] - np_1d[0]

    # 2-variable Newton polytope of charpoly
    cp_expr = cancel(M.charpoly(x).as_expr())
    num, den = fraction(cp_expr)
    p = Poly(num, x, q)
    monoms_num = p.monoms()
    q_shift = Poly(den, q).degree() if den.has(q) else 0
    np_2d = set((xd, qd - q_shift) for (xd, qd) in monoms_num)
    hull_2d = convex_hull_2d(list(np_2d))

    writhe = sum(1 if w > 0 else -1 for w in word)
    genus = Rational(alex_norm, 2)

    return {
        'braid_word': word,
        'n_strands': n,
        'writhe': writhe,
        'burau_matrix': M,
        'alexander_poly': delta,
        'alexander_poly_newton_1d': np_1d,
        'alexander_norm': alex_norm,
        'genus': genus,
        'charpoly_newton_2d': hull_2d,
    }


# Keep old name as alias for backwards compatibility with test script
fibered_face_from_alexander = compute_alexander_data


def print_alexander_data(data):
    """Pretty-print Alexander norm computation."""
    print(f"Braid: {data['braid_word']} in B_{data['n_strands']}")
    print(f"Writhe: {data['writhe']}")
    print(f"Alexander polynomial: {data['alexander_poly']}")
    print(f"  Newton polytope (1D): [{data['alexander_poly_newton_1d'][0]}, "
          f"{data['alexander_poly_newton_1d'][1]}]")
    print(f"  Alexander norm: {data['alexander_norm']}")
    print(f"  Genus (= alex_norm/2): {data['genus']}")
    print(f"Charpoly Newton polytope (x, q): {data['charpoly_newton_2d']}")


# Keep old name as alias
print_fibered_face = print_alexander_data
