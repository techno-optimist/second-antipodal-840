#!/usr/bin/env python3
"""Stage 6: (a) t=1,t=2 both-nonempty exact; (b) randomized maximal even-family check
against classification {slice-15, e7+d4 (8 cells), e8+d4 (15 cells)}."""
import itertools, random, time
from uw import (N, FULL, CELLS, popcount, mask_of, lines_of_cell, conflict,
                build_union_graph, mis_bnb)

def fs(*xs): return frozenset(xs)

def induced(adjl, keep):
    idxs = [j for j in range(len(adjl)) if (keep >> j) & 1]
    pos = {v: k for k, v in enumerate(idxs)}
    sub = [0] * len(idxs)
    for a, v in enumerate(idxs):
        m = adjl[v]
        for w in idxs[a + 1:]:
            if (m >> w) & 1:
                sub[a] |= 1 << pos[w]
                sub[pos[w]] |= 1 << a
    return idxs, sub

A = fs(0, 1, 2, 3)
for B, t in [(fs(3, 4, 5, 6), 1), (fs(2, 3, 4, 5), 2)]:
    L, adj = build_union_graph([A, B])
    nv = len(L)
    ALLM = (1 << nv) - 1
    bidx = next(i for i in range(nv) if L[i][0] == 1)
    best = 0
    for a in range(nv):
        if L[a][0] != 0 or ((adj[bidx] >> a) & 1):
            continue
        rem = ALLM & ~(adj[bidx] | (1 << bidx)) & ~(adj[a] | (1 << a))
        idxs, sub = induced(adj, rem)
        m, _ = mis_bnb(sub)
        best = max(best, 2 + m)
    print(f"t={t} both-nonempty max = {best}")

# (b) randomized maximal even-family check
def even_ok(X, fam):
    return all(len(X & Y) in (0, 2) for Y in fam)

def classify(fam):
    """Return type of maximal even family via span decomposition signature."""
    # compute span over GF(2) as masks
    gens = [mask_of(A) for A in fam]
    basis = []
    for g in gens:
        x = g
        for b in basis:
            x = min(x, x ^ b)
        if x:
            basis.append(x)
    span = {0}
    for b in basis:
        span |= {s ^ b for s in span}
    w4 = [s for s in span if popcount(s) == 4]
    dim = len(basis)
    # support sizes of components: build graph on w4 words by intersecting supports
    return (len(fam), dim, len(w4))

random.seed(42)
sigs = {}
t0 = time.time()
for trial in range(20000):
    fam = []
    order = list(range(len(CELLS)))
    random.shuffle(order)
    for i in order:
        X = CELLS[i]
        if even_ok(X, fam):
            fam.append(X)
    sig = classify(fam)
    sigs[sig] = sigs.get(sig, 0) + 1
print("maximal even family signatures (|fam|, span dim, #wt4 words in span):")
for k, v in sorted(sigs.items()):
    print(" ", k, "x", v)
print(f"({time.time()-t0:.0f}s)")
