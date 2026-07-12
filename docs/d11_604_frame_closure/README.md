# Addendum: the d11 604-point frame is closed — no 605th point in the Bianchi family

Receipts for the companion d11 result referenced by this note's campaign: the
weight-4 sector of the published 604-point kissing record in R^11 has
independence number **exactly 240**, so that family caps at
**8 + 54 + 240 = 302 antipodal lines = 604 points**.

## The claim, exactly

1. **The record is antipodal.** The Bianchi et al. 604-point kissing
   configuration in R^11 (arXiv:2606.10402) is antipodal: 302 lines. The
   antipodality observation is not stated in arXiv:2606.10402 — it is ours,
   documented (with an exact Q(√2) reconstruction) in our antipodal-bounds
   note, DOI [10.5281/zenodo.21285878](https://doi.org/10.5281/zenodo.21285878).
   The frame decomposes as 8 coordinate axes + 54 Q(√2) extension lines +
   240 weight-4 lines on 166 integer supports.

2. **The weight-4 sector reduces to a 1,328-vertex restricted graph** G:
   166 support cells × 8 sign-lines = 1,328 vertices, 25,088 conflict edges,
   ZERO intra-cell edges; all edges lie in exactly 1,568 conflicting cell
   pairs with exactly 16 edges each, and every such pair decomposes into four
   disjoint C4s. Every one of these facts is **recomputed** (not echoed) by
   `check_graph_sanity.py` from the shipped METIS file + mapping, with the
   edge-set SHA-256 pinned to the theta certificate's value.

3. **alpha(G) = 240 EXACTLY** — SCIP branch-and-cut on the exact ILP (edge
   constraints + the 1,568 cell-pair cuts + per-cell caps, warm-started at
   the known 240): **primal 240 = dual 240, gap 0.0, status optimal**
   (2,810 s, 2,697 nodes, 261 solutions; `receipts/scip_ilp_result.json`,
   full solver log `receipts/scip_ilp.log`). `rerun_alpha240.py` reproduces
   the solve from the shipped graph — recompute, never echo.

4. **Hence the frame caps at 8 + 54 + 240 = 302 lines = 604 points.** There
   is no 605th point in this family.

5. **The rigorous bracket 240 ≤ alpha(G) ≤ 251 also holds independently**, by
   the rationally-certified Lovász theta dual: Y rounded to dyadic rationals
   (scale 2^16), exact-dyadic t_cert = 251.91494750976562, and
   M = t·I − J − Y proven PSD by verified floating Cholesky (shift 2^−16
   against a rigorous dpotrf error bound of 4.01·10^−10).
   `check_theta_certificate.py` re-derives PASS from the METIS file +
   `theta_dual_certificate.json` alone; the lower bound is the shipped
   240-line witness `restricted_redumis.is`, whose independence (and
   30-whole-cell structure) `check_graph_sanity.py` recomputes.

6. **All 30 neighboring pair-system designs cap at 302 by exact counting**
   (the B1 census): in the 1/2-hybrid alphabet the exposed weight-4 sector is
   all-or-nothing per cell, so W = 8·P with P the max packing of surviving
   supports, and TOTAL_upper = 8 + |E| + 8·P′ (P′ over all ≥1-survivor
   supports) is a rigorous upper bound per variant. Over all 30
   non-isomorphic pair systems (matchings, paths, stars, triangles, cycles,
   K4, augmented matchings, near-1-factors):
   **max TOTAL_upper = 302 < 303, all_closed_rigorous = true** — the
   published perfect matching is optimal in its neighborhood
   (`B1/B1_census.json`, `B1/B1_upper_bounds.json`).

## What is honestly NOT claimed

- **No pure-Boolean UNSAT certificate.** One was attempted — cube-and-conquer
  (598 cubes) on a 192-core farm, plus monolithic kissat/cadical runs of
  17+ hours — and it did **not** complete: 1 cube proven UNSAT, 718
  first-pass timeouts, 0 SAT (so also: no 241-line witness anywhere).
  alpha(G) = 240 therefore rests on the SCIP branch-and-cut optimality proof
  (gap 0), the rationally-certified theta bracket, and the B1 exact counting
  — not on a DRAT-checked Boolean certificate.
- **This does NOT prove k(11) ≤ 604.** It closes the published frame family
  {8 axes + 54 extension lines + weight-4 lines on the 166 exposed supports}
  and its 30 neighborhood pair-system variants. Other frames, other
  architectures, and non-antipodal configurations are untouched.
- The B1 JSONs are exact-counting **receipts** of the producing campaign
  scripts, shipped for the record; they are not stdlib-recomputable here
  (the producer needs the full Q(√2) frame machinery).

## Verify

```
make verify-d11-closure     # from the repo root:
                            #   check_graph_sanity.py       (stdlib, exact)
                            #   check_theta_certificate.py  (needs numpy)

python3 docs/d11_604_frame_closure/rerun_alpha240.py
                            # reproduce alpha = 240 to proven optimality;
                            # needs pyscipopt (engine of record) or ortools
                            # CP-SAT (independent fallback, assert-optimal);
                            # run of record took 2,810 s.
```

## Layout

```
restricted_graph.metis        the 1,328-vertex / 25,088-edge restricted graph
restricted_graph_mapping.json vertex -> original line / support-cell mapping
restricted_redumis.is         the 240-line independent-set witness (KaMIS)
theta_dual_certificate.json   rational Lovasz-theta dual certificate
check_graph_sanity.py         stdlib exact recomputation of facts 2 and the
                              witness half of fact 5 (alpha >= 240)
check_theta_certificate.py    standalone theta-dual checker (fact 5 upper end;
                              adapted from workspace/kissing_theta_cert_check.py,
                              producer workspace/kissing_theta_dual_certify.py)
rerun_alpha240.py             exact ILP re-solve of fact 3 (SCIP / CP-SAT)
receipts/scip_ilp_result.json the run of record: primal=dual=240, gap 0.0
receipts/scip_ilp.log         full SCIP log of that run
receipts/theta_result.json    primal SCS theta = 251.9155933509097 (numeric)
receipts/certify_run.log      theta dual certification run log
B1/B1_census.json             the 30-variant neighborhood census (fact 6)
B1/B1_upper_bounds.json       per-variant rigorous TOTAL_upper; max = 302
```
