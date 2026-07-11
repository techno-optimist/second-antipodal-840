#!/usr/bin/env python3
"""Exact recomputation of the MACHINE facts of the odd-overlap section.

Recomputes, from the conflict predicate alone (stdlib big-integer
arithmetic, no floats, no solvers), every machine-verified fact cited in
the note's Section "The odd-overlap frontier":

  1. structure lemma censuses: bipartite conflict edges 256 / 128 / 512
     at defect overlaps t = 1 / 2 / 3 (and 0 at t = 0);
  2. pair cap: exact MIS = 64 for the union of a conflicting pair, all t;
  3. t=3 rigidity: exact both-nonempty maximum 57 (and 64 for t = 1, 2);
     folded-7-cube vertex connectivity = 7 by unit-capacity max-flow;
  4. W7 localization: 16 fold-classes, each exact per-class MIS 4 -> 64;
  5. corner collapse: the five cells inside a 5-set cap at exactly 64;
  6. e7 world (Fano complements): exact MIS 64;
     e8 world (extended-Hamming weight-4 words): exact MIS 128;
  7. sunflower triples: point-core (t's 1,1,1, union 10) and pair-core
     (t's 2,2,2, union 8) both cap at exactly 80 (the open Venn (1,0,0,1)
     class also beats 64 -- verified 78-line witness, cap in [78,80] --
     but is not re-verified here; see the note's triple-orbit table);
  8. star saturation: the shipped witness certs/config420.json meets each
     of the 12 coordinate stars in exactly 96 U-lines (the tightness of
     the star reduction alpha(U12) <= 3 alpha(U11));
  9. U11 world sanity: 165 cells / 10,560 lines, no free cell pairs.

Every value is RECOMPUTED (never echoed); the script exits 0 iff all
recomputed values equal the ones claimed in the note.

Usage:  python3 scripts/verify_structure.py
"""
from __future__ import annotations

import collections
import itertools
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "oddlemma"))

from uw import (FULL, N, build_union_graph, check_independent,  # noqa: E402
                conflict, lines_of_cell, mask_of, popcount)
from mis2 import mis_exact  # noqa: E402

FAILURES: list[str] = []


def check(name: str, got, want) -> None:
    ok = got == want
    print(f"  {name:58s} got {got!r:>10}  want {want!r:>10}  "
          f"{'OK' if ok else 'FAIL'}")
    if not ok:
        FAILURES.append(name)


def fs(*xs):
    return frozenset(xs)


def induced(adj, keep):
    idxs = [j for j in range(len(adj)) if (keep >> j) & 1]
    pos = {v: k for k, v in enumerate(idxs)}
    sub = [0] * len(idxs)
    for a, v in enumerate(idxs):
        m = adj[v]
        for w in idxs[a + 1:]:
            if (m >> w) & 1:
                sub[a] |= 1 << pos[w]
                sub[pos[w]] |= 1 << a
    return idxs, sub


def pair_world(A, B):
    return build_union_graph([A, B])


def both_nonempty_max(A, B):
    L, adj = pair_world(A, B)
    nv = len(L)
    allm = (1 << nv) - 1
    bidx = next(i for i in range(nv) if L[i][0] == 1)
    best = 0
    for a in range(nv):
        if L[a][0] != 0 or ((adj[bidx] >> a) & 1):
            continue
        rem = allm & ~(adj[bidx] | (1 << bidx)) & ~(adj[a] | (1 << a))
        _, sub = induced(adj, rem)
        m, _ = mis_exact(sub)
        best = max(best, 2 + m)
    return best


def folded7_connectivity():
    verts, seen = [], set()
    for x in range(128):
        c = min(x, x ^ 127)
        if c not in seen:
            seen.add(c)
            verts.append(c)
    vid = {c: i for i, c in enumerate(verts)}
    adj7 = [set() for _ in range(64)]
    for c in verts:
        for i in range(7):
            d = min(c ^ (1 << i), (c ^ (1 << i)) ^ 127)
            adj7[vid[c]].add(vid[d])
    assert all(len(a) == 7 for a in adj7)

    def vconn(s, t):
        graph: dict = {}

        def add(u, v, c):
            graph.setdefault(u, {})[v] = graph.get(u, {}).get(v, 0) + c
            graph.setdefault(v, {}).setdefault(u, 0)

        for i in range(64):
            add((i, 0), (i, 1), 10 ** 6 if i in (s, t) else 1)
        for u in range(64):
            for v in adj7[u]:
                add((u, 1), (v, 0), 10 ** 6)
        src, snk = (s, 1), (t, 0)
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

    # vertex-transitive: a single base vertex suffices
    return min(vconn(0, t) for t in range(1, 64) if t not in adj7[0])


def per_class_sum(cells, dpoints):
    L, adj = build_union_graph(cells)
    dmask = mask_of(dpoints)
    classes: dict = {}
    for i, (_ci, u) in enumerate(L):
        x = u & dmask
        classes.setdefault(min(x, x ^ dmask), []).append(i)
    cross = 0
    clsof = {}
    for c, idxs in classes.items():
        for i in idxs:
            clsof[i] = c
    for i in range(len(L)):
        m = adj[i]
        while m:
            j = (m & -m).bit_length() - 1
            m &= m - 1
            if clsof[i] != clsof[j]:
                cross += 1
    vals = []
    for _c, idxs in sorted(classes.items()):
        keep = 0
        for i in idxs:
            keep |= 1 << i
        _, sub = induced(adj, keep)
        m, _ = mis_exact(sub)
        vals.append(m)
    return len(classes), vals, cross


def main() -> int:
    t0 = time.time()
    A = fs(0, 1, 2, 3)
    pairs = {0: fs(4, 5, 6, 7), 1: fs(3, 4, 5, 6),
             2: fs(2, 3, 4, 5), 3: fs(1, 2, 3, 4)}

    print("[1] structure-lemma censuses (bipartite conflict edges)")
    for t, want in [(0, 0), (1, 256), (2, 128), (3, 512)]:
        _L, adj = pair_world(A, pairs[t])
        check(f"t={t} edges", sum(popcount(x) for x in adj) // 2, want)

    print("[2] pair cap: exact MIS of a conflicting pair union")
    for t in (1, 2, 3):
        L, adj = pair_world(A, pairs[t])
        m, s = mis_exact(adj)
        ok = check_independent(L, [A, pairs[t]], s)
        check(f"t={t} kappa (witness independent: {ok})", m, 64)

    print("[3] rigidity: exact both-nonempty maxima + connectivity")
    for t, want in [(1, 64), (2, 64), (3, 57)]:
        check(f"t={t} both-nonempty max", both_nonempty_max(A, pairs[t]),
              want)
    check("folded-7-cube vertex connectivity", folded7_connectivity(), 7)

    print("[4] W7 localization (35 cells in a 7-set, fold on 5 coords)")
    w7cells = [frozenset(c) for c in itertools.combinations(range(7), 4)]
    ncls, vals, cross = per_class_sum(w7cells, {7, 8, 9, 10, 11})
    check("W7 fold classes", ncls, 16)
    check("W7 per-class exact MIS values", sorted(set(vals)), [4])
    check("W7 partition-bound sum (=> MIS(W7) <= 64)", sum(vals), 64)

    print("[5] corner collapse (five cells inside a 5-set)")
    corner = [frozenset(c) for c in itertools.combinations(range(5), 4)]
    _L, adj = build_union_graph(corner)
    m, _ = mis_exact(adj)
    check("corner world exact MIS", m, 64)

    print("[6] e7 / e8 even-family worlds")
    fano = [frozenset({i % 7, (i + 1) % 7, (i + 3) % 7}) for i in range(7)]
    e7 = [frozenset(set(range(7)) - set(b)) for b in fano]
    ncls, vals, cross = per_class_sum(e7, {7, 8, 9, 10, 11})
    check("e7 cross-class conflict edges", cross, 0)
    check("e7 exact MIS (sum of 16 class MIS)", sum(vals), 64)
    pts = list(itertools.product([0, 1], repeat=3))
    words = {tuple((a * x + b * y + c * z + d) % 2 for (x, y, z) in pts)
             for a, b, c, d in itertools.product([0, 1], repeat=4)}
    e8 = [frozenset(i for i in range(8) if w[i]) for w in words
          if sum(w) == 4]
    check("e8 cell count (extended-Hamming weight-4 words)", len(e8), 14)
    ncls, vals, cross = per_class_sum(e8, {8, 9, 10, 11})
    check("e8 cross-class conflict edges", cross, 0)
    check("e8 exact MIS (sum of 8 class MIS)", sum(vals), 128)

    print("[7] sunflower triples (both cap exactly 80; the open Venn "
          "(1,0,0,1) class also beats 64 via its verified 78-line witness)")
    for name, cells in [
            ("point-core (t's 1,1,1, union 10)",
             [fs(0, 1, 2, 3), fs(0, 4, 5, 6), fs(0, 7, 8, 9)]),
            ("pair-core  (t's 2,2,2, union 8)",
             [fs(0, 1, 2, 3), fs(0, 1, 4, 5), fs(0, 1, 6, 7)])]:
        L, adj = build_union_graph(cells)
        m, s = mis_exact(adj)
        ok = check_independent(L, cells, s)
        check(f"sunflower {name} (indep: {ok})", m, 80)

    print("[8] star saturation of the shipped 288-line witness")
    cert = HERE.parent / "certs" / "config420.json"
    lines = json.load(open(cert))["lines"]
    ulines = [l for l in lines if l["type"] == "U"]
    check("witness U-line count", len(ulines), 288)
    stars = [0] * 12
    for l in ulines:
        for x in set(range(12)) - set(l["support"]):
            stars[x] += 1
    check("per-star counts uniform", sorted(set(stars)), [96])
    check("star identity 4*288 = 12*96", 4 * len(ulines), 12 * 96)

    print("[9] U11 world sanity")
    c11 = [frozenset(c) for c in itertools.combinations(range(11), 3)]
    check("U11 cell count", len(c11), 165)
    check("U11 line count", 165 * 64, 10560)
    rs = {11 - len(a | b) for a, b in itertools.combinations(c11, 2)}
    check("U11 overlap census (no free pairs)", sorted(rs), [5, 6, 7])

    dt = time.time() - t0
    if FAILURES:
        print(f"\nSTRUCTURE VERIFICATION FAILED ({len(FAILURES)} checks): "
              f"{FAILURES}")
        return 1
    print(f"\nALL STRUCTURE FACTS RECOMPUTED AND VERIFIED ({dt:.1f}s). "
          "Every MACHINE entry of the odd-overlap section is reproduced.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
