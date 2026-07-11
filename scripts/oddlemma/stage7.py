#!/usr/bin/env python3
"""Stage 7: W7 localization upper bound via per-delta-class exact MIS.
W7 = all 35 cells inside {0..6}; D = {7..11} (5 common support coords).
Valid bound: MIS(W7) <= sum over 16 delta-classes of per-class MIS
(dropping cross-class edges only relaxes). If every class MIS = 4 -> W7 = 64.
"""
import itertools, time, json
from uw import (N, FULL, popcount, mask_of, lines_of_cell, conflict,
                build_union_graph, mis_bnb, check_independent)

W7 = list(range(7))
cells = [frozenset(c) for c in itertools.combinations(W7, 4)]
assert len(cells) == 35
t0 = time.time()
L, adj = build_union_graph(cells)
print(f"W7 world: n={len(L)}, edges={sum(popcount(x) for x in adj)//2}, build {time.time()-t0:.0f}s", flush=True)

Dmask = mask_of({7, 8, 9, 10, 11})
classes = {}
for i, (ci, u) in enumerate(L):
    x = u & Dmask
    c = min(x, x ^ Dmask)
    classes.setdefault(c, []).append(i)
print(f"classes: {len(classes)}, sizes {sorted(set(len(v) for v in classes.values()))}", flush=True)

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

tot = 0
vals = []
for c, idxs in sorted(classes.items()):
    keep = 0
    for i in idxs:
        keep |= 1 << i
    _, sub = induced(adj, keep)
    tt = time.time()
    m, _ = mis_bnb(sub)
    vals.append(m)
    tot += m
    print(f"class {bin(c)}: size {len(idxs)} MIS={m} ({time.time()-tt:.1f}s)", flush=True)
print(f"per-class values: {vals}")
print(f"sum of per-class MIS = {tot}  => MIS(W7) <= {tot}", flush=True)
json.dump({"W7_per_class": vals, "W7_upper": tot}, open("w7_result.json", "w"))
