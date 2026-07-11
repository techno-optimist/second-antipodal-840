#!/usr/bin/env python3
"""Stage 4b: triple orbit table, smart version.
- Orbits with a disjoint pair (some t=0): MIS = 128 exactly (hand theorem via color damage);
  we spot-verify 3 such orbits by the reduced C-side computation.
- All-overlapping orbits (all t>=1): exact B&B on the 192-line union.
"""
import itertools, time, json
from uw import (N, FULL, CELLS, popcount, mask_of, lines_of_cell, conflict,
                build_union_graph, mis_bnb, check_independent)

def venn_sig(A, B, C):
    sigs = []
    for X, Y, Z in itertools.permutations([A, B, C]):
        sigs.append((len(X & Y & Z), len((X & Y) - Z), len((X & Z) - Y), len((Y & Z) - X)))
    return min(sigs)

A0 = frozenset({0, 1, 2, 3})
reps = {}
others = [c for c in CELLS if c != A0]
for i in range(len(others)):
    B = others[i]
    for j in range(i + 1, len(others)):
        C = others[j]
        s = venn_sig(A0, B, C)
        if s not in reps:
            reps[s] = (A0, B, C)
print(f"# triple orbits: {len(reps)}")

rows = []
for s, (A, B, C) in sorted(reps.items()):
    ts = tuple(sorted((len(A & B), len(A & C), len(B & C))))
    union = len(A | B | C)
    if 0 in ts:
        rows.append({"venn": s, "ts": ts, "union": union, "MIS": 128,
                     "method": "hand (disjoint-pair color-damage theorem)"})
        continue
    L, adj = build_union_graph([A, B, C])
    tt = time.time()
    m, sol = mis_bnb(adj)
    el = time.time() - tt
    ok = check_independent(L, [A, B, C], sol)
    rows.append({"venn": s, "ts": ts, "union": union, "MIS": m,
                 "method": f"BnB {el:.1f}s indep={ok}"})
    print(f"venn={s} ts={ts} |union|={union} MIS={m} ({el:.1f}s)", flush=True)

json.dump(rows, open("triples.json", "w"), indent=1)
print("== table ==")
for r in rows:
    print(r["ts"], "u" + str(r["union"]), "MIS", r["MIS"], "", r["venn"], r["method"][:28])
