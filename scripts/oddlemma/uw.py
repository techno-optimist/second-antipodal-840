#!/usr/bin/env python3
r"""U-world exact toolkit (stdlib only).

U-world: cells = 4-subsets A of [12] (complement of the 8-support H = [12]\A).
Lines per cell: sign vectors eps in {+-1}^H, even # of -1, mod global flip: 64 lines.
Conflict((A,u),(B,v)) iff |sum_{i in H_A cap H_B} eps_i eps'_i| >= 5  (|cos| > 1/2).

Line encoding: (A frozenset, u int bitmask over [12] with bits only on H=[12]\A,
popcount(u) even), canonical rep = min(u, u ^ maskH).
"""
import itertools, sys
from functools import lru_cache

N = 12
FULL = (1 << N) - 1

def popcount(x):  # py3.8+ int.bit_count in 3.10
    return x.bit_count()

CELLS = [frozenset(c) for c in itertools.combinations(range(N), 4)]  # cell = 4-set A
CELL_INDEX = {c: i for i, c in enumerate(CELLS)}

def mask_of(s):
    m = 0
    for i in s:
        m |= 1 << i
    return m

def lines_of_cell(A):
    """Return list of canonical line masks u (subsets of H=[12]\\A, even weight, mod ^maskH)."""
    Hmask = FULL & ~mask_of(A)
    Hbits = [i for i in range(N) if (Hmask >> i) & 1]
    out = []
    seen = set()
    for bits in range(256):
        u = 0
        for j in range(8):
            if (bits >> j) & 1:
                u |= 1 << Hbits[j]
        if popcount(u) % 2:
            continue
        cu = min(u, u ^ Hmask)
        if cu in seen:
            continue
        seen.add(cu)
        out.append(cu)
    assert len(out) == 64
    return out

def conflict(A, u, B, v):
    """Exact conflict predicate. A,B cells (frozensets); u,v canonical masks."""
    if A == B:
        return False  # same cell always compatible (u==v handled by caller as same line)
    R = (FULL & ~mask_of(A)) & (FULL & ~mask_of(B))  # common support mask
    r = popcount(R)
    if r <= 4:
        return False
    # sum of eps_i eps'_i over R = r - 2*wt((u^v)&R)
    d = popcount((u ^ v) & R)
    return abs(r - 2 * d) >= 5

def build_union_graph(cells):
    """Given a list of cells (frozensets), build vertex list [(cellidx_local,u)] and
    adjacency as bitsets. Returns (verts, adj) with verts[i]=(ci,u)."""
    L = []
    for ci, A in enumerate(cells):
        for u in lines_of_cell(A):
            L.append((ci, u))
    n = len(L)
    adj = [0] * n
    for i in range(n):
        ci, u = L[i]
        A = cells[ci]
        for j in range(i + 1, n):
            cj, v = L[j]
            if ci == cj:
                continue
            if conflict(A, u, cells[cj], v):
                adj[i] |= 1 << j
                adj[j] |= 1 << i
    return L, adj

# ---------- exact MIS (branch & bound with greedy-coloring bound), bitset sets ----------
def mis_bnb(adj, time_budget=None):
    """Exact max independent set on graph given as list of int bitmasks.
    Works on complement: MIS(G) = max clique(complement). We do direct MIS B&B
    with greedy coloring bound on the complement-clique formulation? Simpler:
    classic MIS B&B: bound = greedy clique-cover size of remaining vertices."""
    n = len(adj)
    ALL = (1 << n) - 1
    best = [0]
    best_set = [0]

    order = sorted(range(n), key=lambda i: -popcount(adj[i]))

    def clique_cover_bound(P):
        """Greedy clique cover of P: number of cliques >= MIS upper bound...
        Actually MIS(P) <= #cliques in any clique cover. Greedy: repeatedly grow cliques."""
        cnt = 0
        rem = P
        while rem:
            # start new clique with lowest set bit
            v = (rem & -rem).bit_length() - 1
            clique = 1 << v
            cand = rem & adj[v]
            rem &= ~(1 << v)
            while cand:
                w = (cand & -cand).bit_length() - 1
                clique |= 1 << w
                cand &= adj[w]
                rem &= ~(1 << w)
                cand &= rem
            cnt += 1
        return cnt

    def expand(P, size, cur):
        if size + clique_cover_bound(P) <= best[0]:
            return
        if not P:
            if size > best[0]:
                best[0] = size
                best_set[0] = cur
            return
        # pick vertex of max degree within P
        # (branch: v in MIS -> remove N[v]; v not in MIS -> remove v)
        bestv, bd = -1, -1
        Q = P
        while Q:
            v = (Q & -Q).bit_length() - 1
            Q &= Q - 1
            d = popcount(adj[v] & P)
            if d > bd:
                bd, bestv = d, v
        v = bestv
        if bd == 0:
            # P is independent: take all
            s = size + popcount(P)
            if s > best[0]:
                best[0] = s
                best_set[0] = cur | P
            return
        # branch 1: v in
        expand(P & ~(adj[v] | (1 << v)), size + 1, cur | (1 << v))
        # branch 2: v out
        expand(P & ~(1 << v), size, cur)

    expand(ALL, 0, 0)
    return best[0], best_set[0]

def check_independent(Lverts, cells, chosen_mask):
    """Verify chosen set is independent (exact)."""
    idxs = [i for i in range(len(Lverts)) if (chosen_mask >> i) & 1]
    for a in range(len(idxs)):
        i = idxs[a]
        ci, u = Lverts[i]
        for b in range(a + 1, len(idxs)):
            j = idxs[b]
            cj, v = Lverts[j]
            if ci == cj:
                continue
            if conflict(cells[ci], u, cells[cj], v):
                return False
    return True

if __name__ == "__main__":
    # structural self-tests
    A = frozenset({0, 1, 2, 3})
    # t=1
    B1 = frozenset({3, 4, 5, 6})
    # t=2
    B2 = frozenset({2, 3, 4, 5})
    # t=3
    B3 = frozenset({1, 2, 3, 4})
    # t=0
    B0 = frozenset({4, 5, 6, 7})
    for B, t, expected_edges in [(B1, 1, 256), (B2, 2, 128), (B3, 3, 512), (B0, 0, 0)]:
        L, adj = build_union_graph([A, B])
        e = sum(popcount(x) for x in adj) // 2
        assert e == expected_edges, (t, e)
        print(f"t={t}: edges={e} OK")
    # kappa via exact MIS on pairs
    for B, t in [(B1, 1), (B2, 2), (B3, 3)]:
        L, adj = build_union_graph([A, B])
        m, s = mis_bnb(adj)
        ok = check_independent(L, [A, B], s)
        print(f"t={t}: kappa={m} (independent verified: {ok})")
    print("selftests done")
