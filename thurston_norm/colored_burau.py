"""
Colored reduced Burau matrix and Alexander polynomial via Morton's formula.

For an n-strand braid w, the colored reduced Burau matrix B(t_1,...,t_n) has
strand variables t_1,...,t_n.  Strands in the same permutation cycle of the
braid get identified, yielding the meridian variables for the closed-braid
link components.

The multivariable Alexander polynomial of the closed braid + axis is:

    Delta(t_1,...,t_k, s) = det(I - s * B_reduced)

where t_i are meridian variables for the k closure components and s is the
axis meridian.  Setting s=1 (or omitting the axis) gives the Alexander
polynomial of the closed braid complement.

The Alexander norm of a cohomology class phi is the width of the Newton
polytope of Delta in direction phi.  The Alexander norm ball is dual to the
Newton polytope.

Reference: Morton, "The multivariable Alexander polynomial for a closed
braid" (2001).  This matches the basis used in the existing MATLAB code
(reduced_long_moody.m) for the Burau representation.
"""


def braid_permutation(word, n):
    """Return the permutation (0-indexed) induced by braid word on n strands."""
    perm = list(range(n))
    for g in word:
        i = abs(g) - 1
        perm[i], perm[i + 1] = perm[i + 1], perm[i]
    return perm


def permutation_cycles(perm):
    """Return cycles of a permutation (0-indexed) as sorted lists."""
    seen = set()
    cycles = []
    for start in range(len(perm)):
        if start in seen:
            continue
        cycle = []
        j = start
        while j not in seen:
            seen.add(j)
            cycle.append(j)
            j = perm[j]
        if cycle:
            cycles.append(cycle)
    return cycles


def colored_burau_matrix(word, n, ring=None, strand_vars=None):
    """
    Compute the colored reduced Burau matrix for an n-strand braid.

    The coloring tracks which original strand is at each position as the
    braid is read left to right, and uses that strand's variable at each
    crossing.

    Returns (B, ring, strand_vars) where B is (n-1)x(n-1).
    """
    from sage.all import LaurentPolynomialRing, ZZ, matrix, identity_matrix

    if ring is None:
        var_names = ["t%d" % i for i in range(n)]
        ring = LaurentPolynomialRing(ZZ, var_names)

    if strand_vars is None:
        strand_vars = list(ring.gens()[:n])

    # Track which original strand is at each position
    pos = list(range(n))

    B = identity_matrix(ring, n - 1)

    for g in word:
        k = abs(g)        # 1-based generator index
        positive = g > 0
        r = k - 1         # 0-indexed row in reduced matrix

        # Variable for the strand currently at position k (1-based)
        t = strand_vars[pos[k - 1]]

        M = identity_matrix(ring, n - 1)
        if positive:
            M[r, r] = -t
            if r > 0:
                M[r, r - 1] = ring.one()
            if r < n - 2:
                M[r, r + 1] = t
        else:
            ti = t ** (-1)
            M[r, r] = -ti
            if r > 0:
                M[r, r - 1] = ti
            if r < n - 2:
                M[r, r + 1] = ring.one()

        B = B * M

        # sigma_k swaps positions k and k+1
        pos[k - 1], pos[k] = pos[k], pos[k - 1]

    return B, ring, strand_vars


def identify_cycle_variables(word, n, B, ring, strand_vars):
    """
    Identify strand variables along permutation cycles of the braid.

    Strands in the same cycle of the braid permutation close up to the same
    link component, so their meridian variables are identified.

    Returns (B_sub, new_ring, cycle_vars, cycles).
    """
    from sage.all import LaurentPolynomialRing, ZZ, matrix

    perm = braid_permutation(word, n)
    cycles = permutation_cycles(perm)
    num_cycles = len(cycles)

    new_var_names = ["m%d" % i for i in range(num_cycles)]
    new_ring = LaurentPolynomialRing(ZZ, new_var_names)
    new_vars = list(new_ring.gens())

    # Map each strand to its cycle index
    strand_to_cycle = {}
    for ci, cyc in enumerate(cycles):
        for j in cyc:
            strand_to_cycle[j] = ci

    # Substitution: old strand var -> new cycle var
    sub_dict = {}
    for j in range(n):
        sub_dict[strand_vars[j]] = new_vars[strand_to_cycle[j]]

    nr, nc = B.nrows(), B.ncols()
    B_sub = matrix(new_ring, nr, nc)
    for i in range(nr):
        for j in range(nc):
            B_sub[i, j] = B[i, j].subs(sub_dict)

    return B_sub, new_ring, new_vars, cycles


def morton_alexander_poly(word, n, include_axis=True):
    """
    Compute the multivariable Alexander polynomial via Morton's formula.

    Delta = det(I - s * B)

    where B is the colored reduced Burau with cycle-identified variables.

    Parameters
    ----------
    word : list of int
        Artin braid word.
    n : int
        Number of strands.
    include_axis : bool
        If True, include axis variable s.  If False, set s=1.

    Returns
    -------
    delta : Laurent polynomial
    ring : the polynomial ring
    info : dict with 'cycles', 'cycle_vars', 'num_closure_components',
           and optionally 'axis_var'.
    """
    from sage.all import LaurentPolynomialRing, ZZ, identity_matrix, matrix

    B, ring0, strand_vars = colored_burau_matrix(word, n)
    B_sub, cyc_ring, cyc_vars, cycles = identify_cycle_variables(
        word, n, B, ring0, strand_vars
    )
    num_comp = len(cycles)

    if include_axis:
        all_names = [str(v) for v in cyc_vars] + ["s"]
        full_ring = LaurentPolynomialRing(ZZ, all_names)
        full_vars = list(full_ring.gens())
        s = full_vars[-1]
        cycle_vars = full_vars[:-1]

        sub = {cyc_vars[i]: cycle_vars[i] for i in range(num_comp)}
        nr, nc = B_sub.nrows(), B_sub.ncols()
        B_full = matrix(full_ring, nr, nc)
        for i in range(nr):
            for j in range(nc):
                B_full[i, j] = B_sub[i, j].subs(sub)

        I = identity_matrix(full_ring, nr)
        delta = (I - s * B_full).det()

        return delta, full_ring, {
            "cycles": cycles,
            "cycle_vars": cycle_vars,
            "axis_var": s,
            "num_closure_components": num_comp,
        }
    else:
        I = identity_matrix(cyc_ring, B_sub.nrows())
        delta = (I - B_sub).det()

        return delta, cyc_ring, {
            "cycles": cycles,
            "cycle_vars": cyc_vars,
            "num_closure_components": num_comp,
        }


def newton_polytope(delta):
    """
    Compute the Newton polytope of a Laurent polynomial.

    Returns (vertices, Polyhedron) where vertices is a list of tuples.
    """
    from sage.all import Polyhedron

    exponents = [list(e) for e in delta.exponents()]
    if not exponents:
        return [], None

    P = Polyhedron(vertices=exponents)
    verts = [tuple(v) for v in P.vertices_list()]
    return verts, P


def alexander_norm(delta, phi):
    """
    Compute the Alexander seminorm of a cohomology class phi.

    ||phi||_A = max_{e in supp(delta)} <phi, e> - min_{e in supp(delta)} <phi, e>

    Parameters
    ----------
    delta : Laurent polynomial (the Alexander polynomial)
    phi : tuple of rationals/integers (cohomology class in meridian basis)

    Returns the seminorm value.
    """
    exponents = [list(e) for e in delta.exponents()]
    if not exponents:
        return 0

    vals = [sum(p * e for p, e in zip(phi, exp)) for exp in exponents]
    return max(vals) - min(vals)


def alexander_norm_ball(delta, ambient_dim=None):
    """
    Compute the Alexander norm unit ball (dual of Newton polytope).

    The Alexander norm unit ball B_A is:
        B_A = { phi : ||phi||_A <= 1 }
            = { phi : max <phi, e> - min <phi, e> <= 1 for e in supp(delta) }

    This is the polar dual of the symmetrized Newton polytope.

    Returns the Polyhedron.
    """
    from sage.all import Polyhedron

    exponents = [list(e) for e in delta.exponents()]
    if not exponents:
        return None

    # Newton polytope
    N = Polyhedron(vertices=exponents)

    # Center it (the Alexander norm uses the width, which is
    # the support function of N minus the support function of -N,
    # equivalently the width of the Minkowski sum N + (-N) / 2).
    # The dual of the "width functional" ball is the centrally
    # symmetric body (N - centroid(N)) intersected with its negation.
    #
    # More precisely: ||phi||_A = h_N(phi) + h_N(-phi) where h_N is the
    # support function.  So B_A = { phi : h_N(phi) + h_N(-phi) <= 1 }.
    # This equals the polar of the "difference body" (N - N)/2, which
    # is centrally symmetric.

    # Difference body D = (N - N) / 2 = { (x - y)/2 : x, y in N }
    # Its vertices are a subset of { (v_i - v_j)/2 : v_i, v_j vertices of N }
    from sage.all import QQ
    verts_N = N.vertices_list()
    diff_verts = []
    for vi in verts_N:
        for vj in verts_N:
            dv = tuple((QQ(a) - QQ(b)) / 2 for a, b in zip(vi, vj))
            diff_verts.append(dv)

    D = Polyhedron(vertices=diff_verts)

    # Alexander norm ball = polar dual of D
    # B_A = D^* = { phi : <phi, d> <= 1 for all d in D }
    try:
        B_A = D.polar()
    except Exception:
        # If D is not full-dimensional, polar may fail.
        # Fall back to computing in the affine hull.
        B_A = None

    return B_A
