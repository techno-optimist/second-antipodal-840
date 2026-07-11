#!/usr/bin/env python3
"""Stronger exact MIS: degree-0/1 reduction + component decomposition, then
incumbent-pruned B&B (clique-cover bound) per component. Exact by construction."""

def popcount(x):
    return x.bit_count()

def _bits(m):
    while m:
        b = m & -m
        yield b.bit_length() - 1
        m ^= b

def _bnb_component(adj, P0):
    """Exact MIS of induced subgraph P0 (assumed nonempty), global-incumbent B&B."""
    best = [0, 0]

    def clique_cover_bound(P):
        cnt = 0
        rem = P
        while rem:
            v = (rem & -rem).bit_length() - 1
            cand = rem & adj[v]
            rem &= ~(1 << v)
            while cand:
                w = (cand & -cand).bit_length() - 1
                cand &= adj[w] & rem
                rem &= ~(1 << w)
            cnt += 1
        return cnt

    def expand(P, size, cur):
        if not P:
            if size > best[0]:
                best[0], best[1] = size, cur
            return
        if size + clique_cover_bound(P) <= best[0]:
            return
        # reductions inside: take deg-0 and deg-1 vertices greedily
        changed = True
        while changed:
            changed = False
            for v in _bits(P):
                if not ((P >> v) & 1):
                    continue
                d = adj[v] & P
                pc = popcount(d)
                if pc == 0:
                    size += 1
                    cur |= 1 << v
                    P &= ~(1 << v)
                    changed = True
                elif pc == 1:
                    size += 1
                    cur |= 1 << v
                    P &= ~((1 << v) | d)
                    changed = True
            if not P:
                if size > best[0]:
                    best[0], best[1] = size, cur
                return
        bestv, bd = -1, -1
        for v in _bits(P):
            dd = popcount(adj[v] & P)
            if dd > bd:
                bd, bestv = dd, v
        v = bestv
        expand(P & ~(adj[v] | (1 << v)), size + 1, cur | (1 << v))
        expand(P & ~(1 << v), size, cur)

    expand(P0, 0, 0)
    return best[0], best[1]

def _components(adj, P):
    comps = []
    rem = P
    while rem:
        seed = rem & -rem
        comp = seed
        frontier = seed
        while frontier:
            nxt = 0
            for v in _bits(frontier):
                nxt |= adj[v] & P
            nxt &= ~comp
            comp |= nxt
            frontier = nxt
        comps.append(comp)
        rem &= ~comp
    return comps

def mis_exact(adj):
    n = len(adj)
    P = (1 << n) - 1
    total, chosen = 0, 0
    # global degree-0/1 reduction first
    changed = True
    while changed:
        changed = False
        for v in _bits(P):
            if not ((P >> v) & 1):
                continue
            d = adj[v] & P
            pc = popcount(d)
            if pc == 0:
                total += 1
                chosen |= 1 << v
                P &= ~(1 << v)
                changed = True
            elif pc == 1:
                total += 1
                chosen |= 1 << v
                P &= ~((1 << v) | d)
                changed = True
    for comp in _components(adj, P):
        m, s = _bnb_component(adj, comp)
        total += m
        chosen |= s
    # verify
    for i in _bits(chosen):
        assert not (adj[i] & chosen & ~(1 << i)), "not independent"
    assert popcount(chosen) == total
    return total, chosen

if __name__ == "__main__":
    import random, sys
    sys.path.insert(0, ".")
    from uw import mis_bnb
    random.seed(1)
    for trial in range(300):
        n = random.randint(1, 42)
        adj = [0] * n
        for i in range(n):
            for j in range(i + 1, n):
                if random.random() < random.choice([0.05, 0.15, 0.4, 0.7]):
                    adj[i] |= 1 << j
                    adj[j] |= 1 << i
        m1, _ = mis_bnb(adj)
        m2, _ = mis_exact(adj)
        assert m1 == m2, (trial, n, m1, m2)
    print("mis2 cross-check vs mis_bnb: 300 random graphs OK")
