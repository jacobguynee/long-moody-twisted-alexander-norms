#!/usr/bin/env python3
"""
Analyze vertex counts of Newton polytopes from our symbolic results.
Question: is the number of hull vertices uniformly bounded?
"""

# Collected from all symbolic runs. Format: (name, hull_vertices)
# B3 Burau (2x2, charpoly degree 2 in x, so x in {0,1,2})

b3_data = [
    ('s1.s2^-1',           [(0, 0), (1, -1), (1, 1), (2, 0)]),
    ('s1^2.s2^-1',         [(0, 1), (1, -1), (1, 2), (2, 0)]),
    ('s1.s2^-2',           [(0, -1), (1, -2), (1, 1), (2, 0)]),
    ('s1^2.s2^-2',         [(0, 0), (1, -2), (1, 2), (2, 0)]),
    ('s1^3.s2^-1',         [(0, 2), (1, -1), (1, 3), (2, 0)]),
    ('s1.s2^-3',           [(0, -2), (1, -3), (1, 1), (2, 0)]),
    ('s1.s2^-1.s1.s2^-1',  [(0, 0), (1, -2), (1, 2), (2, 0)]),
]

# B4 Burau (3x3, charpoly degree 3 in x, so x in {0,1,2,3})

b4_alternating = [
    ('s1.s2^-1.s3',                [(0, 1), (1, 2), (2, -1), (3, 0)]),
    ('s3.s2^-1.s1',                [(0, 1), (1, 2), (2, -1), (3, 0)]),
    ('s1.s2^-1.s3.s2^-1',          [(0, 0), (1, 2), (2, -2), (3, 0)]),
    ('s1^2.s2^-1.s3',              [(0, 2), (1, 0), (1, 3), (2, -1), (2, 2), (3, 0)]),
    ('s1.s2^-2.s3',                [(0, 0), (1, 2), (2, -2), (3, 0)]),
    ('s1.s2^-1.s3^2',              [(0, 2), (1, 0), (1, 3), (2, -1), (2, 2), (3, 0)]),
    ('s1^2.s2^-2.s3^2',            [(0, 2), (1, 4), (2, -2), (3, 0)]),
    ('s1.s2^-1.s3.s1.s2^-1.s3',    [(0, 2), (1, 4), (2, -2), (3, 0)]),
    ('s1^2.s2^-1.s3.s2^-1.s1',     [(0, 2), (1, -1), (1, 4), (2, -2), (2, 3), (3, 0)]),
    ('s1.s2^-2.s3^2.s2^-1.s1',     [(0, 1), (1, 4), (2, -3), (3, 0)]),
    ('s1^3.s2^-2.s3',              [(0, 2), (1, -1), (1, 4), (2, -2), (2, 3), (3, 0)]),
    ('s1.s2^-3.s3^2',              [(0, 0), (1, -2), (1, 3), (2, -3), (2, 2), (3, 0)]),
    ('s1.s2^-1.s3.s1.s2^-1.s3.s1', [(0, 3), (1, 0), (1, 5), (2, -2), (2, 3), (3, 0)]),
    ('s1^2.s2^-1.s3^2.s2^-1.s1^2', [(0, 4), (1, 0), (1, 6), (2, -2), (2, 4), (3, 0)]),
    ('s3.s2^-1.s1^2.s2^-1.s3.s2^-1', [(0, 1), (1, 4), (2, -3), (3, 0)]),
    ('s1.s2^-1.s3.s2^-1.s1.s2^-1.s3', [(0, 1), (1, 4), (2, -3), (3, 0)]),
]

b4_fg_family1 = [  # f=s1, g=s2^2.s3.s2^-2
    # same sign positive
    ('f^5.g^5',   [(0, 10), (1, 0), (2, 10), (3, 0)]),
    ('f^6.g^4',   [(0, 10), (1, 0), (2, 10), (3, 0)]),
    ('f^8.g^7',   [(0, 15), (1, 0), (2, 15), (3, 0)]),
    ('f^10.g^7',  [(0, 17), (1, 0), (2, 17), (3, 0)]),
    ('f^15.g^15', [(0, 30), (1, 0), (2, 30), (3, 0)]),
    ('f^25.g^5',  [(0, 30), (1, 0), (2, 30), (3, 0)]),
    # same sign negative
    ('f^-8.g^-7',   [(0, -15), (1, 0), (2, -15), (3, 0)]),
    ('f^-10.g^-7',  [(0, -17), (1, 0), (2, -17), (3, 0)]),
    # mixed sign
    ('f^5.g^-5',  [(0, 0), (1, -5), (1, 5), (2, -5), (2, 5), (3, 0)]),
    ('f^-5.g^5',  [(0, 0), (1, -5), (1, 5), (2, -5), (2, 5), (3, 0)]),
    ('f^8.g^-7',  [(0, 1), (1, -7), (1, 8), (2, -7), (2, 8), (3, 0)]),
    ('f^-7.g^8',  [(0, 1), (1, -7), (1, 8), (2, -7), (2, 8), (3, 0)]),
    ('f^12.g^-3', [(0, 9), (1, -3), (1, 12), (2, -3), (2, 12), (3, 0)]),
    ('f^20.g^-10',[(0, 10), (1, -10), (1, 20), (2, -10), (2, 20), (3, 0)]),
    ('f^-15.g^15',[(0, 0), (1, -15), (1, 15), (2, -15), (2, 15), (3, 0)]),
    ('f^25.g^-10',[(0, 15), (1, -10), (1, 25), (2, -10), (2, 25), (3, 0)]),
]

b4_fg_family2 = [  # f=(s1.s2.s1)^2, g=s3
    # same sign positive
    ('f^8.g^7',   [(0, 55), (1, 25), (1, 54), (2, 1), (2, 30), (3, 0)]),
    ('f^10.g^5',  [(0, 65), (1, 31), (1, 64), (2, 1), (2, 34), (3, 0)]),
    ('f^1.g^14',  [(0, 20), (1, 4), (1, 19), (2, 1), (2, 16), (3, 0)]),
    ('f^15.g^15', [(0, 105), (1, 46), (1, 104), (2, 1), (2, 59), (3, 0)]),
    ('f^25.g^25', [(0, 175), (1, 76), (1, 174), (2, 1), (2, 99), (3, 0)]),
    # same sign negative
    ('f^-8.g^-7',  [(0, -55), (1, -54), (1, -25), (2, -30), (2, -1), (3, 0)]),
    ('f^-10.g^-5', [(0, -65), (1, -64), (1, -31), (2, -34), (2, -1), (3, 0)]),
    # mixed sign
    ('f^10.g^-5',  [(0, 55), (1, 60), (2, -5), (3, 0)]),
    ('f^-5.g^10',  [(0, -20), (1, -30), (2, 10), (3, 0)]),
    ('f^8.g^-7',   [(0, 41), (1, 48), (2, -7), (3, 0)]),
    ('f^-7.g^8',   [(0, -34), (1, -42), (2, 8), (3, 0)]),
    ('f^20.g^-10', [(0, 110), (1, 120), (2, -10), (3, 0)]),
    ('f^-15.g^15', [(0, -75), (1, -90), (2, 15), (3, 0)]),
    ('f^40.g^-10', [(0, 230), (1, 240), (2, -10), (3, 0)]),
    ('f^-25.g^25', [(0, -125), (1, -150), (2, 25), (3, 0)]),
]

b4_positive = [  # positive braids (from symbolic run)
    ('(123)^3.1',  [(0, 10), (1, 6), (2, 4), (3, 0)]),
    ('(321)^3.3',  [(0, 10), (1, 6), (2, 4), (3, 0)]),
    ('1.3.2...',   [(0, 10), (1, 6), (2, 4), (3, 0)]),
    ('1.2.1...',   [(0, 10), (1, 6), (2, 4), (3, 0)]),
    ('(123)^4',    [(0, 12), (3, 0)]),
    ('(1213)^3',   [(0, 12), (3, 0)]),
    ('(2132)^3',   [(0, 12), (1, 6), (2, 6), (3, 0)]),
    ('(1323)^3',   [(0, 12), (3, 0)]),
    ('(12321)^3',  [(0, 15), (1, 12), (2, 3), (3, 0)]),
    ('(13213)^3',  [(0, 15), (1, 12), (2, 3), (3, 0)]),
    ('(21321)^3',  [(0, 15), (1, 9), (2, 6), (3, 0)]),
    ('(123)^6',    [(0, 18), (3, 0)]),
    ('(123212)^3', [(0, 18), (3, 0)]),
    ('(132132)^3', [(0, 18), (3, 0)]),
    ('(213231)^3', [(0, 18), (3, 0)]),
    ('(12132)^4',  [(0, 20), (1, 12), (2, 8), (3, 0)]),
    ('(12321)^4',  [(0, 20), (1, 16), (2, 4), (3, 0)]),
    ('(31213)^4',  [(0, 20), (1, 16), (2, 4), (3, 0)]),
    ('(23123)^4',  [(0, 20), (1, 12), (2, 8), (3, 0)]),
    ('(123)^8',    [(0, 24), (3, 0)]),
]

print('=== Vertex count analysis ===\n')

all_datasets = [
    ('B3 alternating', b3_data),
    ('B4 alternating', b4_alternating),
    ('B4 f=s1,g=s2^2.s3.s2^-2', b4_fg_family1),
    ('B4 f=(s1s2s1)^2,g=s3', b4_fg_family2),
    ('B4 positive braids', b4_positive),
]

global_max = 0
for dataset_name, dataset in all_datasets:
    counts = [len(h) for _, h in dataset]
    print(f'{dataset_name}:')
    print(f'  {len(dataset)} braids, vertex counts: {sorted(set(counts))}')
    print(f'  min={min(counts)}, max={max(counts)}')
    global_max = max(global_max, max(counts))

    # Show distribution
    from collections import Counter
    dist = Counter(counts)
    for k in sorted(dist):
        print(f'    {k} vertices: {dist[k]} braids')
    print()

print(f'GLOBAL MAXIMUM vertex count: {global_max}')
print()

# Theoretical analysis
print('=== Theoretical bound ===')
print()
print('Charpoly det(xI - M) has degree n in x (n = matrix size).')
print('Newton polytope in (x, q) has x-coordinates in {0, 1, ..., n}.')
print('At each x-level, at most 2 hull vertices (min and max q-exponent).')
print('Total: at most 2(n+1) hull vertices.')
print()
print('B3: n=2, bound = 2*3 = 6. Observed max = ', max(len(h) for _,h in b3_data))
print('B4: n=3, bound = 2*4 = 8. Observed max = ', max(max(len(h) for _,h in d) for _,d in all_datasets if 'B4' in _))
print()

# Check: are there ever vertices at EVERY x-level?
print('=== x-levels with hull vertices ===')
for dataset_name, dataset in all_datasets:
    if 'B3' in dataset_name:
        n = 2
    else:
        n = 3
    for name, hull in dataset:
        x_levels = set(v[0] for v in hull)
        verts_per_level = {}
        for v in hull:
            verts_per_level.setdefault(v[0], []).append(v[1])
        # Check if any level has 2 vertices
        multi = {xl: vs for xl, vs in verts_per_level.items() if len(vs) == 2}
        if len(hull) == 2*(n+1):  # maximum possible
            print(f'  MAX VERTICES: {name}: {hull}')
