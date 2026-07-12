#!/usr/bin/env python3
"""Recompute alpha(G) = 240 for the d11 604-frame restricted graph.

This is a RERUN script, not an echo: it rebuilds the exact ILP from the
shipped METIS graph + mapping and solves it to proven optimality.  The model
is identical to the run of record (receipts/scip_ilp_result.json + full SCIP
log receipts/scip_ilp.log — primal 240 = dual 240, gap 0.0, status optimal):

    max  sum_v x_v                        x binary, one var per sign-line
    s.t. x_a + x_b <= 1                   for each of the 25,088 conflict edges
         sum_{v in cell_i} x_v
       + sum_{v in cell_j} x_v <= 8       for each of the 1,568 conflicting
                                          cell pairs (implied cuts: the pair
                                          conflict graph is four disjoint C4s)
         sum_{v in cell_i} x_v <= 8       per-cell cap (166 cells)

Engines (auto-selected, or forced with --engine):

  scip   pyscipopt branch-and-cut — the engine of record.  Closure is
         PROVEN when status == "optimal" (primal = dual = 240, gap 0).
  cpsat  ortools CP-SAT fallback — an independent exact engine.  Closure is
         PROVEN when the solver returns OPTIMAL; the script then ASSERTS
         status == OPTIMAL and objective == 240 (a 241-line set would
         instead surface as objective >= 241 = an immediate d11 record).

Either engine is warm-started/hinted with the shipped 240-line witness so
the primal side is never the bottleneck; all remaining work is the proof
that 241 is infeasible.  Run of record: 2,810 s wall on a Grace-Blackwell
DGX Spark (SCIP); budget accordingly.

Usage:
    python3 rerun_alpha240.py [--engine auto|scip|cpsat]
                              [--time-limit SECONDS]     (default 21600)
                              [--out RECEIPT.json]

Exit 0 with "ALPHA_240_PROVEN" iff optimality was proven at 240.
Exit 1 (honest failure) on timeout / unavailable engine — never a silent pass.
"""
import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load_instance(metis_path, mapping_path, witness_path):
    lines = open(metis_path).read().strip().split("\n")
    n, m = map(int, lines[0].split())
    edges = []
    for i, ln in enumerate(lines[1 : 1 + n]):
        for tok in ln.split():
            j = int(tok) - 1
            if i < j:
                edges.append((i, j))
    assert len(edges) == m, f"parsed {len(edges)} edges, header says {m}"

    mapping = json.load(open(mapping_path))
    cell = [o // 8 for o in mapping["old_line_indices"]]
    cells = defaultdict(list)
    for v, c in enumerate(cell):
        cells[c].append(v)
    pairs = set()
    for a, b in edges:
        pairs.add((min(cell[a], cell[b]), max(cell[a], cell[b])))

    known = [i for i, x in enumerate(open(witness_path).read().split()) if x == "1"]
    assert len(known) == 240, f"witness has {len(known)} lines, want 240"
    return n, m, edges, cells, sorted(pairs), known


def solve_scip(n, edges, cells, pairs, known, time_limit):
    from pyscipopt import Model, quicksum

    model = Model("alpha_604frame_rerun")
    x = {v: model.addVar(vtype="B", name=f"x{v}") for v in range(n)}
    for a, b in edges:
        model.addCons(x[a] + x[b] <= 1)
    for (i, j) in pairs:
        model.addCons(quicksum(x[v] for v in cells[i]) +
                      quicksum(x[v] for v in cells[j]) <= 8)
    for c, vs in cells.items():
        model.addCons(quicksum(x[v] for v in vs) <= 8)
    model.setObjective(quicksum(x[v] for v in range(n)), "maximize")

    sol = model.createSol()  # warm start with the shipped 240 witness
    kn = set(known)
    for v in range(n):
        model.setSolVal(sol, x[v], 1.0 if v in kn else 0.0)
    model.addSol(sol, free=True)

    model.setRealParam("limits/time", time_limit)
    model.optimize()

    primal = model.getObjVal() if model.getNSols() > 0 else None
    dual = model.getDualbound()
    return {
        "engine": "scip",
        "status": model.getStatus(),
        "primal_best": primal,
        "dual_bound": dual,
        "gap": model.getGap(),
        "proven_240": (model.getStatus() == "optimal"
                       and primal is not None and round(primal) == 240
                       and dual < 241.0),
        "found_241": (primal is not None and primal >= 240.5),
    }


def solve_cpsat(n, edges, cells, pairs, known, time_limit):
    from ortools.sat.python import cp_model

    model = cp_model.CpModel()
    x = [model.NewBoolVar(f"x{v}") for v in range(n)]
    for a, b in edges:
        model.Add(x[a] + x[b] <= 1)
    for (i, j) in pairs:
        model.Add(sum(x[v] for v in cells[i]) + sum(x[v] for v in cells[j]) <= 8)
    for c, vs in cells.items():
        model.Add(sum(x[v] for v in vs) <= 8)
    model.Maximize(sum(x))

    kn = set(known)  # hint with the shipped 240 witness
    for v in range(n):
        model.AddHint(x[v], 1 if v in kn else 0)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    status = solver.Solve(model)
    status_name = solver.StatusName(status)
    obj = int(solver.ObjectiveValue()) if status in (cp_model.OPTIMAL, cp_model.FEASIBLE) else None
    proven = status == cp_model.OPTIMAL and obj == 240  # assert-optimal
    return {
        "engine": "cpsat",
        "status": status_name,
        "primal_best": obj,
        "dual_bound": solver.BestObjectiveBound(),
        "gap": None,
        "proven_240": proven,
        "found_241": (obj is not None and obj >= 241),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine", choices=["auto", "scip", "cpsat"], default="auto")
    ap.add_argument("--time-limit", type=float, default=21600.0)
    ap.add_argument("--graph", default=str(HERE / "restricted_graph.metis"))
    ap.add_argument("--mapping", default=str(HERE / "restricted_graph_mapping.json"))
    ap.add_argument("--witness", default=str(HERE / "restricted_redumis.is"))
    ap.add_argument("--out", default=None, help="optional receipt JSON path")
    args = ap.parse_args()

    n, m, edges, cells, pairs, known = load_instance(args.graph, args.mapping, args.witness)
    print(f"instance: {n} vertices, {m} edges, {len(cells)} cells, "
          f"{len(pairs)} conflicting cell pairs, warm start {len(known)}")

    engine = args.engine
    if engine == "auto":
        try:
            import pyscipopt  # noqa: F401
            engine = "scip"
        except ImportError:
            try:
                import ortools  # noqa: F401
                engine = "cpsat"
            except ImportError:
                print("NO ENGINE: install pyscipopt (preferred) or ortools", file=sys.stderr)
                sys.exit(1)
    print(f"engine: {engine}, time limit {args.time_limit:.0f}s")

    solve = solve_scip if engine == "scip" else solve_cpsat
    res = solve(n, edges, cells, pairs, known, args.time_limit)
    res.update({"n": n, "m": m, "time_limit": args.time_limit})
    print(json.dumps(res, indent=2))
    if args.out:
        json.dump(res, open(args.out, "w"), indent=2)

    if res["found_241"]:
        print("FOUND >= 241: this would be a NEW d11 kissing record — verify independently!")
        sys.exit(1)  # contradicts the run of record; do not treat as a pass
    if res["proven_240"]:
        print("ALPHA_240_PROVEN: alpha(G) = 240 exactly; the 604-frame caps at "
              "8 + 54 + 240 = 302 lines = 604 points.")
        sys.exit(0)
    print("NOT CLOSED within the time limit (honest failure; raise --time-limit).")
    sys.exit(1)


if __name__ == "__main__":
    main()
