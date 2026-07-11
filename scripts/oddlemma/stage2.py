#!/usr/bin/env python3
"""Stage 2: color structure, t=3 rigidity, corner collapse, e7/e8 worlds."""
import itertools, json, time
from uw import (N, FULL, CELLS, popcount, mask_of, lines_of_cell, conflict,
                build_union_graph, mis_bnb, check_independent)

def fs(*xs):
    return frozenset(xs)

# ---------- 1. color structure verification ----------
def color_structure(A, B):
    """Return dict color -> (list A-lines, list B-lines) using restriction to R mod 1_R."""
    Rmask = (FULL & ~mask_of(A)) & (FULL & ~mask_of(B))
    def colors(C):
        Hmask = FULL & ~mask_of(C)
        out = {}
        for u in lines_of_cell(C):
            c = min(u & Rmask, (u & Rmask) ^ Rmask)
            out.setdefault(c, []).append(u)
        return out
    return colors(A), colors(B), Rmask

A = fs(0, 1, 2, 3)
print("== t=1 color check ==")
B = fs(3, 4, 5, 6)
ca, cb, Rm = color_structure(A, B)
assert len(ca) == 16 and all(len(v) == 4 for v in ca.values())
assert len(cb) == 16 and all(len(v) == 4 for v in cb.values())
# conflicts exactly = same color, all 16 pairs
bad = 0
for c1, ulist in ca.items():
    for c2, vlist in cb.items():
        for u in ulist:
            for v in vlist:
                conf = conflict(A, u, B, v)
                if (c1 == c2) != conf:
                    bad += 1
print("t=1: 16 colors x (4,4); conflict==samecolor mismatches:", bad)

print("== t=2 color check ==")
B = fs(2, 3, 4, 5)
ca, cb, Rm = color_structure(A, B)
assert len(ca) == 32 and all(len(v) == 2 for v in ca.values())
bad = 0
for c1, ulist in ca.items():
    for c2, vlist in cb.items():
        for u in ulist:
            for v in vlist:
                if (c1 == c2) != conflict(A, u, B, v):
                    bad += 1
print("t=2: 32 colors x (2,2); mismatches:", bad)

print("== t=3 folded-cube check ==")
B = fs(1, 2, 3, 4)
Rmask = (FULL & ~mask_of(A)) & (FULL & ~mask_of(B))
assert popcount(Rmask) == 7
def fold_class(u):
    x = u & Rmask
    return min(x, x ^ Rmask)
# each cell's 64 lines biject onto 64 fold classes
la, lb = lines_of_cell(A), lines_of_cell(B)
fa = {fold_class(u): u for u in la}
fb = {fold_class(v): v for v in lb}
assert len(fa) == 64 and len(fb) == 64
def fold_dist(x, y):
    d = popcount((x ^ y) & Rmask)
    return min(d, 7 - d)
bad = 0
for x, u in fa.items():
    for y, v in fb.items():
        if (fold_dist(x, y) <= 1) != conflict(A, u, B, v):
            bad += 1
print("t=3: bijection to F_2^7/1; conflict==folded-dist<=1 mismatches:", bad)

# ---------- 2. t=3 rigidity: both nonempty => <= 63 ----------
print("== t=3 rigidity ==")
cells = [A, B]
L, adj = build_union_graph(cells)
nA = sum(1 for (ci, u) in L if ci == 0)
best_both = 0
t0 = time.time()
for i in range(len(L)):
    if L[i][0] != 0:
        continue
    # force line i (in A); then find MIS of remaining graph restricted to non-neighbors,
    # but require >=1 in B: since we maximize total, just compute and check B-count.
    P = ((1 << len(L)) - 1) & ~(adj[i] | (1 << i))
    # simple approach: MIS on induced subgraph P with vertex relabeling
    idxs = [j for j in range(len(L)) if (P >> j) & 1]
    sub = [0] * len(idxs)
    pos = {v: k for k, v in enumerate(idxs)}
    for a2, v in enumerate(idxs):
        m = adj[v]
        for w in idxs[a2 + 1:]:
            if (m >> w) & 1:
                sub[a2] |= 1 << pos[w]
                sub[pos[w]] |= 1 << a2
    m2, s2 = mis_bnb(sub)
    # check whether optimum uses any B line; the MIS includes all of A-side always?
    tot = 1 + m2
    # count B lines in s2
    bcount = sum(1 for k in range(len(idxs)) if (s2 >> k) & 1 and L[idxs[k]][0] == 1)
    if bcount >= 1:
        best_both = max(best_both, tot)
    else:
        # optimum for this branch avoided B entirely; also try forcing a B vertex:
        # find max over j in B compatible with i
        for k, v in enumerate(idxs):
            if L[v][0] != 1:
                continue
            P2 = [x & ~((sub[k]) | (1 << k)) for x in sub]
            # restrict to non-neighbors of k
            keep = [(x, xi) for xi, x in enumerate(idxs)]
            mask2 = ((1 << len(idxs)) - 1) & ~(sub[k] | (1 << k))
            idxs3 = [j for j in range(len(idxs)) if (mask2 >> j) & 1]
            sub3 = [0] * len(idxs3)
            pos3 = {v2: k2 for k2, v2 in enumerate(idxs3)}
            for a3, v3 in enumerate(idxs3):
                mrow = sub[v3]
                for w3 in idxs3[a3 + 1:]:
                    if (mrow >> w3) & 1:
                        sub3[a3] |= 1 << pos3[w3]
                        sub3[pos3[w3]] |= 1 << a3
            m3, _ = mis_bnb(sub3)
            best_both = max(best_both, 2 + m3)
        break  # symmetry: one i is enough if we did the forced-B sweep... keep safe: continue loop
print(f"t=3 max with both cells nonempty = {best_both} (elapsed {time.time()-t0:.1f}s)")

# ---------- 3. corner collapse: five 4-subsets of {0..4} ----------
print("== corner collapse (5 cells in a 5-set) ==")
corner = [frozenset(c) for c in itertools.combinations(range(5), 4)]
t0 = time.time()
L, adj = build_union_graph(corner)
m, s = mis_bnb(adj)
print(f"corner world: n={len(L)}, MIS={m} (elapsed {time.time()-t0:.1f}s), indep={check_independent(L, corner, s)}")
