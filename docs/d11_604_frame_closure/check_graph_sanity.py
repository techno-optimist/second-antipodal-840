#!/usr/bin/env python3
"""Stdlib exact sanity checker for the d11 604-frame restricted graph.

Recomputes (never echoes) every structural fact the addendum README states
about the weight-4 restricted conflict graph, from the shipped METIS file
and mapping alone:

  1. n = 1328 vertices, m = 25,088 edges; adjacency lists symmetric.
  2. The edge-set SHA-256 equals the constant pinned below AND the value
     pinned inside theta_dual_certificate.json (so the theta certificate
     provably talks about this exact graph).
  3. The 1,328 vertices partition into 166 support cells of exactly 8
     sign-lines each (cell = old_line_index // 8 from the mapping), with
     ZERO intra-cell edges.
  4. All 25,088 edges lie in exactly 1,568 conflicting cell pairs with
     exactly 16 edges each; within every such pair the conflict graph is
     2-regular bipartite and decomposes into FOUR DISJOINT C4s.
  5. The shipped independent-set witness (restricted_redumis.is) has
     exactly 240 vertices, ZERO conflict edges, and consists of exactly
     30 whole cells (8 lines each)  =>  alpha(G) >= 240.

Usage: python3 check_graph_sanity.py [graph.metis mapping.json witness.is cert.json]
(defaults: the copies shipped next to this script).  Exit 0 iff ALL checks pass.
"""
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
EDGES_SHA256 = "237a7f3f031b7be64dbe6bf3a5b42c8a06a39f08e249054fa33fc469a9369809"


def fail(msg):
    print(f"SANITY FAILED: {msg}")
    sys.exit(1)


def main():
    metis = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "restricted_graph.metis"
    mapping = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE / "restricted_graph_mapping.json"
    witness = Path(sys.argv[3]) if len(sys.argv) > 3 else HERE / "restricted_redumis.is"
    cert = Path(sys.argv[4]) if len(sys.argv) > 4 else HERE / "theta_dual_certificate.json"

    # -- 1. parse + symmetry ------------------------------------------------
    lines = metis.read_text().strip().split("\n")
    n, m = map(int, lines[0].split())
    if len(lines) != 1 + n:
        fail(f"METIS row count {len(lines) - 1} != n = {n}")
    adj = [set() for _ in range(n)]
    edges = set()
    for i, ln in enumerate(lines[1 : 1 + n]):
        for tok in ln.split():
            j = int(tok) - 1
            if not (0 <= j < n) or j == i:
                fail(f"bad neighbor {j} in row {i}")
            adj[i].add(j)
            edges.add((min(i, j), max(i, j)))
    if any(i not in adj[j] for i in range(n) for j in adj[i]):
        fail("adjacency lists are not symmetric")
    edges = sorted(edges)
    if (n, m, len(edges)) != (1328, 25088, 25088):
        fail(f"(n, m, |E|) = {(n, m, len(edges))} != (1328, 25088, 25088)")
    print(f"PASS  graph: n = {n}, m = {m}, adjacency symmetric")

    # -- 2. pinned hash (script constant + theta certificate agree) --------
    h = hashlib.sha256(json.dumps(edges).encode()).hexdigest()
    if h != EDGES_SHA256:
        fail(f"edge-set sha256 {h} != pinned {EDGES_SHA256}")
    cert_h = json.load(open(cert))["graph"]["edges_sha256"]
    if cert_h != EDGES_SHA256:
        fail(f"theta certificate pins a DIFFERENT graph ({cert_h})")
    print(f"PASS  edge-set sha256 pinned: {h[:16]}... (matches theta certificate)")

    # -- 3. cell partition, zero intra-cell edges ---------------------------
    mp = json.load(open(mapping))
    old = mp["old_line_indices"]
    if len(old) != n:
        fail(f"mapping has {len(old)} old_line_indices != n")
    cell = [o // 8 for o in old]
    sizes = defaultdict(int)
    for c in cell:
        sizes[c] += 1
    if len(sizes) != 166 or set(sizes.values()) != {8}:
        fail(f"cell partition is {len(sizes)} cells, sizes {set(sizes.values())} (want 166 x 8)")
    intra = sum(1 for a, b in edges if cell[a] == cell[b])
    if intra != 0:
        fail(f"{intra} intra-cell edges (want 0)")
    print("PASS  166 support cells x 8 sign-lines, ZERO intra-cell edges")

    # -- 4. cell-pair structure: 1568 pairs x 16 edges, four disjoint C4s --
    pairs = defaultdict(list)
    for a, b in edges:
        pairs[(min(cell[a], cell[b]), max(cell[a], cell[b]))].append((a, b))
    if len(pairs) != 1568 or set(len(v) for v in pairs.values()) != {16}:
        fail(f"{len(pairs)} conflicting cell pairs, edge counts "
             f"{set(len(v) for v in pairs.values())} (want 1568 x {{16}})")
    for key, es in pairs.items():
        loc = defaultdict(list)
        for a, b in es:
            loc[a].append(b)
            loc[b].append(a)
        if any(len(v) != 2 for v in loc.values()):
            fail(f"cell pair {key} is not 2-regular")
        seen, cyc = set(), []
        for v in loc:
            if v in seen:
                continue
            length, prev, cur = 0, None, v
            while True:
                seen.add(cur)
                length += 1
                nxt = [w for w in loc[cur] if w != prev]
                prev, cur = cur, (nxt[0] if nxt else loc[cur][0])
                if cur == v:
                    break
            cyc.append(length)
        if sorted(cyc) != [4, 4, 4, 4]:
            fail(f"cell pair {key} decomposes as {sorted(cyc)} (want four C4s)")
    print("PASS  1,568 conflicting cell pairs x 16 edges, each = four disjoint C4s")

    # -- 5. the 240-line witness => alpha >= 240 ----------------------------
    bits = witness.read_text().split()
    if len(bits) != n or set(bits) - {"0", "1"}:
        fail("witness file is not n lines of 0/1")
    W = {i for i, x in enumerate(bits) if x == "1"}
    if len(W) != 240:
        fail(f"witness has {len(W)} vertices (want 240)")
    viol = sum(1 for a, b in edges if a in W and b in W)
    if viol != 0:
        fail(f"witness has {viol} conflict edges (want 0)")
    wc = defaultdict(int)
    for v in W:
        wc[cell[v]] += 1
    if len(wc) != 30 or set(wc.values()) != {8}:
        fail(f"witness cell census {len(wc)} cells, sizes {set(wc.values())} (want 30 whole cells)")
    print("PASS  witness: 240 lines, 0 conflicts, exactly 30 whole cells  =>  alpha(G) >= 240")

    print("ALL GRAPH SANITY CHECKS PASSED")


if __name__ == "__main__":
    main()
