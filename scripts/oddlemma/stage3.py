#!/usr/bin/env python3
"""Stage 3: t=3 exact rigidity (clean), e7 world, e8 world with delta-split."""
import itertools, time, json
from uw import (N, FULL, CELLS, popcount, mask_of, lines_of_cell, conflict,
                build_union_graph, mis_bnb, check_independent)

def fs(*xs): return frozenset(xs)

def induced(adj, keep_mask):
    idxs = [j for j in range(len(adj)) if (keep_mask >> j) & 1]
    pos = {v: k for k, v in enumerate(idxs)}
    sub = [0] * len(idxs)
    for a, v in enumerate(idxs):
        m = adj[v]
        for w in idxs[a + 1:]:
            if (m >> w) & 1:
                sub[a] |= 1 << pos[w]
                sub[pos[w]] |= 1 << a
    return idxs, sub

# ---------- 1. t=3 exact both-nonempty maximum, by double forcing ----------
print("== t=3 both-nonempty exact ==")
A, B = fs(0, 1, 2, 3), fs(1, 2, 3, 4)
cells = [A, B]
L, adj = build_union_graph(cells)
nv = len(L)
ALLM = (1 << nv) - 1
t0 = time.time()
best = 0
# fix b = first B-line (transitivity of stabilizer on B-lines), sweep all compatible a
bidx = next(i for i in range(nv) if L[i][0] == 1)
forced_removed = adj[bidx] | (1 << bidx)
for a in range(nv):
    if L[a][0] != 0 or ((adj[bidx] >> a) & 1):
        continue
    rem = ALLM & ~forced_removed & ~(adj[a] | (1 << a))
    idxs, sub = induced(adj, rem)
    m, _ = mis_bnb(sub)
    best = max(best, 2 + m)
print(f"t=3 both-nonempty max = {best} (elapsed {time.time()-t0:.1f}s)")

# verify folded-7-cube connectivity = 7 by exact max-flow (Menger), stdlib
print("== folded 7-cube connectivity check ==")
# vertices: classes of F_2^7 mod 1; adjacency: xor e_i
verts7 = []
seen = set()
for x in range(128):
    c = min(x, x ^ 127)
    if c not in seen:
        seen.add(c)
        verts7.append(c)
vid = {c: i for i, c in enumerate(verts7)}
adj7 = [[] for _ in range(64)]
for c in verts7:
    for i in range(7):
        d = min(c ^ (1 << i), (c ^ (1 << i)) ^ 127)
        if vid[d] not in adj7[vid[c]]:
            adj7[vid[c]].append(vid[d])
assert all(len(a) == 7 for a in adj7)

def vertex_connectivity_pair(s, t, adjlist, n):
    """max # vertex-disjoint s-t paths via unit-capacity node-split max flow (BFS)."""
    # node i -> in=2i, out=2i+1, cap 1 edge in->out (except s,t: inf)
    import collections
    SZ = 2 * n
    cap = {}
    def add(u, v, c):
        cap[(u, v)] = cap.get((u, v), 0) + c
        cap.setdefault((v, u), 0)
    for i in range(n):
        add(2 * i, 2 * i + 1, 10**6 if i in (s, t) else 1)
    for u in range(n):
        for v in adjlist[u]:
            add(2 * u + 1, 2 * v, 10**6)
    src, snk = 2 * s + 1, 2 * t
    flow = 0
    while True:
        # BFS augment
        par = {src: None}
        dq = collections.deque([src])
        while dq:
            u = dq.popleft()
            if u == snk:
                break
            for v in [w for (x, w) in cap if x == u]:
                pass
            # slow; build adjacency once instead
            break
        break
    return None  # replaced below

# faster: dinic-lite
def vconn(s, t, adjlist, n):
    graph = {}
    def add(u, v, c):
        graph.setdefault(u, {})[v] = graph.get(u, {}).get(v, 0) + c
        graph.setdefault(v, {}).setdefault(u, 0)
    for i in range(n):
        add((i, 0), (i, 1), 10**6 if i in (s, t) else 1)
    for u in range(n):
        for v in adjlist[u]:
            add((u, 1), (v, 0), 10**6)
    src, snk = (s, 1), (t, 0)
    import collections
    flow = 0
    while True:
        par = {src: None}
        dq = collections.deque([src])
        found = False
        while dq and not found:
            u = dq.popleft()
            for v, c in graph[u].items():
                if c > 0 and v not in par:
                    par[v] = u
                    if v == snk:
                        found = True
                        break
                    dq.append(v)
        if not found:
            return flow
        v = snk
        while par[v] is not None:
            u = par[v]
            graph[u][v] -= 1
            graph[v][u] += 1
            v = u
        flow += 1
        if flow > 8:
            return flow

t0 = time.time()
minconn = 99
v0 = 0
for tv in range(1, 64):
    f = vconn(v0, tv, adj7, 64)
    if tv not in [x for x in adj7[v0]]:
        minconn = min(minconn, f)
    else:
        # adjacent pairs: connectivity defined via non-adjacent pairs; still record
        pass
# also need pairs not involving v0? vertex-transitive => v0 sweep suffices
print(f"folded-7-cube kappa(v0,t) min over non-neighbors = {minconn} (elapsed {time.time()-t0:.1f}s)")

# ---------- 2. e7 world ----------
print("== e7 world ==")
S7 = list(range(7))
fano_lines = [frozenset({i % 7, (i + 1) % 7, (i + 3) % 7}) for i in range(7)]
e7_cells = [frozenset(set(S7) - set(b)) for b in fano_lines]
for X, Y in itertools.combinations(e7_cells, 2):
    assert len(X & Y) == 2
t0 = time.time()
L7, adj_e7 = build_union_graph(e7_cells)
print(f"e7 world: n={len(L7)}, edges={sum(popcount(x) for x in adj_e7)//2}")
# delta-split on D' = {7,8,9,10,11}
Dmask = mask_of({7, 8, 9, 10, 11})
def dclass(u):
    x = u & Dmask
    return min(x, x ^ Dmask)
classes = {}
for i, (ci, u) in enumerate(L7):
    classes.setdefault(dclass(u), []).append(i)
print(f"delta classes: {len(classes)} sizes {sorted(set(len(v) for v in classes.values()))}")
# verify no cross-class edges
cross = 0
clsof = {}
for c, idxs in classes.items():
    for i in idxs:
        clsof[i] = c
for i in range(len(L7)):
    m = adj_e7[i]
    while m:
        j = (m & -m).bit_length() - 1
        m &= m - 1
        if clsof[i] != clsof[j]:
            cross += 1
print(f"cross-class conflict edges: {cross}")
tot = 0
per = []
for c, idxs in sorted(classes.items()):
    keep = 0
    for i in idxs:
        keep |= 1 << i
    id2, sub = induced(adj_e7, keep)
    m, s = mis_bnb(sub)
    per.append(m)
    tot += m
print(f"e7 per-class MIS: {per}")
print(f"kappa7 = MIS(e7 world) = {tot} (elapsed {time.time()-t0:.1f}s)")

# ---------- 3. e8 world ----------
print("== e8 world ==")
# RM(1,3) = extended Hamming [8,4,4] on {0..7}: evaluate affine functions on F_2^3
pts = list(itertools.product([0, 1], repeat=3))
words = set()
for a1, a2, a3, b in itertools.product([0, 1], repeat=4):
    w = tuple((a1 * x + a2 * y + a3 * z + b) % 2 for (x, y, z) in pts)
    words.add(w)
w4 = [frozenset(i for i in range(8) if w[i]) for w in words if sum(w) == 4]
assert len(w4) == 14
for X, Y in itertools.combinations(w4, 2):
    assert len(X & Y) in (0, 2)
e8_cells = w4
t0 = time.time()
L8, adj_e8 = build_union_graph(e8_cells)
print(f"e8 world: n={len(L8)}, edges={sum(popcount(x) for x in adj_e8)//2}")
Dmask = mask_of({8, 9, 10, 11})
classes = {}
for i, (ci, u) in enumerate(L8):
    x = u & Dmask
    c = min(x, x ^ Dmask)
    classes.setdefault(c, []).append(i)
print(f"delta classes: {len(classes)} sizes {sorted(set(len(v) for v in classes.values()))}")
clsof = {}
for c, idxs in classes.items():
    for i in idxs:
        clsof[i] = c
cross = 0
for i in range(len(L8)):
    m = adj_e8[i]
    while m:
        j = (m & -m).bit_length() - 1
        m &= m - 1
        if clsof[i] != clsof[j]:
            cross += 1
print(f"cross-class conflict edges: {cross}")
tot = 0
per = {}
for c, idxs in sorted(classes.items()):
    keep = 0
    for i in idxs:
        keep |= 1 << i
    id2, sub = induced(adj_e8, keep)
    tt = time.time()
    m, s = mis_bnb(sub)
    per[c] = m
    tot += m
    print(f"  class {c:04b}-ish size {len(idxs)}: MIS={m} ({time.time()-tt:.1f}s)")
print(f"kappa8 = MIS(e8 world) = {tot} (elapsed {time.time()-t0:.1f}s)")
json.dump({"kappa7": None, "e8_per_class": {str(k): v for k, v in per.items()},
           "kappa8": tot}, open("e8_result.json", "w"), indent=1)
