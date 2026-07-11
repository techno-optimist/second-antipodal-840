#!/usr/bin/env python3
"""Stage 8: W8 per-class MIS (70 cells inside {0..7}; D = {8..11}).
Valid upper bound: MIS(W8) <= sum of 8 per-class MIS. Conjecture per-class 16 -> W8 = 128.
"""
import itertools, time, json
from uw import (FULL, popcount, mask_of, build_union_graph, mis_bnb)

W8 = list(range(8))
cells = [frozenset(c) for c in itertools.combinations(W8, 4)]
assert len(cells) == 70
t0 = time.time()
L, adj = build_union_graph(cells)
print(f"W8 world: n={len(L)}, edges={sum(popcount(x) for x in adj)//2}, build {time.time()-t0:.0f}s", flush=True)

Dmask = mask_of({8, 9, 10, 11})
classes = {}
for i, (ci, u) in enumerate(L):
    x = u & Dmask
    classes.setdefault(min(x, x ^ Dmask), []).append(i)
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
print(f"per-class: {vals}; sum = {tot} => MIS(W8) <= {tot}", flush=True)
json.dump({"W8_per_class": vals, "W8_upper": tot}, open("w8_result.json", "w"))
