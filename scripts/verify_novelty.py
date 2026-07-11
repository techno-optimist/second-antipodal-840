#!/usr/bin/env python3
"""Exact verification of the NOVELTY theorem for the 420-line witness.

THEOREM CERTIFIED.  No line configuration of the form

    "12 axis lines + 51 whole W4 cells"          (the A+W4 family)

— a family that contains the KNOWN modular 420-line configuration in R^12 —
realizes the pairwise |cos| histogram of the witness in certs/config420.json.
Since the multiset of pairwise |cos| values is an isometry invariant, the
witness is NOT isometric to the known modular 420 (nor to any other member
of the family).  This is a packing-independent COUNTING argument: no solver,
no isometry search.

Vocabulary: an axis line is +-e_k; a WHOLE W4 CELL on a 4-coordinate support
S is the set of ALL 8 lines with entries +-1/2 on S (16 sign patterns mod
antipodality).

WHAT THIS SCRIPT CHECKS, all in exact integer/Fraction arithmetic:

  STEP 1 (recompute, never echo): the witness |cos| histogram, rebuilt from
    certs/config420.json — followed by a hard PACKING PRECONDITION: any
    histogram bin with |cos| > 1/2 rejects the input outright (it is not a
    60-degree configuration, so no novelty verdict applies).  This checker
    is therefore sound standalone, not only jointly with verify_420.

  STEP 2 (first-principles constants): every interaction constant of the
    A+W4 family is DERIVED by explicit finite computation on concrete cells —
    within-cell pairs, cross-cell pairs for supports sharing s = 0,1,2,3
    points, axis-cell pairs, axis-axis pairs.  In particular the computation
    itself shows that two cells sharing 3 points produce |cos| = 3/4 > 1/2,
    so share-3 pairs are FORBIDDEN in any 60-degree family (this is where
    "blocks pairwise share <= 2 points" comes from — it is derived, not
    assumed).

  STEP 3 (linear forcing): with P_s = #(cell pairs sharing exactly s points)
    and r_i = #(cells containing point i), the derived constants force
        n_{1/4} = 64 * P1          and   n_{1/2} = 2448 + 32 * P2,
    so matching the witness histogram forces P1 = 506 and P2 = 563
    (any divisibility/negativity failure would itself prove the theorem).

  STEP 4 (convexity forcing): double counting gives
        sum_i C(r_i, 2) = P1 + 2*P2 = 1632,
    and with sum_i r_i = 4*51 = 204 over 12 points the exact identity
        sum_i C(r_i,2) = 1632 + (1/2) * sum_i (r_i - 17)^2
    forces r_i = 17 for every point.

  STEP 5 (pair-degree contradiction): the pair degrees lambda_p = #(cells on
    coordinate pair p) satisfy lambda_p <= 5 (the two extra points of each
    cell through p are pairwise disjoint among the remaining 10 coordinates,
    because share-3 is forbidden by STEP 2).  For each point i,
    sum_{j != i} lambda_ij = 3 * r_i = 51 over 11 values, and an exhaustive
    DP minimization proves sum_j C(lambda_ij, 2) >= 94.  Summing over the 12
    points counts each pair twice:  2 * P2 >= 12 * 94 = 1128, so P2 >= 564.
    But STEP 3 forced P2 = 563.  CONTRADICTION — the theorem holds.

If at any stage the contradiction chain fails to close, the script reports
NOVELTY NOT ESTABLISHED and exits nonzero.  Exit 0 iff the theorem verifies.

Usage:  python3 scripts/verify_novelty.py [certs/novelty.json [certs/config420.json]]
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from fractions import Fraction as Fr
from itertools import combinations, product
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from witness import CERTS, CertificateError, dot, exact_abs_cos, fail, load_witness

HALF = Fr(1, 2)
QUARTER = Fr(1, 4)


# ------------------------- STEP 1: witness histogram -----------------------

def witness_histogram(wit_path):
    d, lines = load_witness(wit_path)
    hist: Counter = Counter()
    n = len(lines)
    for i in range(n):
        _, gi, Ni = lines[i]
        for j in range(i + 1, n):
            _, gj, Nj = lines[j]
            hist[exact_abs_cos(dot(gi, gj), Ni, Nj)] += 1
    return d, n, hist


# ---------------- STEP 2: first-principles family constants ----------------

def whole_cell(support, d):
    """All 8 lines of the whole W4 cell on `support`, as integer reps
    (entries +-1 on support, first support entry +1 canonical), norm^2 = 4."""
    reps = []
    for signs in product((1, -1), repeat=3):
        g = [0] * d
        g[support[0]] = 1
        for i, s in zip(support[1:], signs):
            g[i] = s
        reps.append(tuple(g))
    return reps


def pair_hist(A, B):
    """|cos| histogram between two disjoint line families of W4-type reps
    (norm^2 = 4 each => |cos| = |<g,g'>| / 4)."""
    h: Counter = Counter()
    for a in A:
        for b in B:
            h[Fr(abs(dot(a, b)), 4)] += 1
    return dict(h)


def derive_constants(d):
    """Derive every pairwise-interaction constant of the A+W4 family by
    explicit computation on concrete cells.  Coordinate relabeling maps any
    support pair sharing s points to the representatives used here, and a
    whole cell is sign-symmetric, so these constants are fully general."""
    c = {}
    base = whole_cell((0, 1, 2, 3), d)
    if len(set(base)) != 8:
        fail("internal: whole cell does not contain 8 distinct lines")

    # within one whole cell: C(8,2) = 28 pairs
    h: Counter = Counter()
    for a, b in combinations(base, 2):
        h[Fr(abs(dot(a, b)), 4)] += 1
    c["within"] = dict(h)

    # two whole cells whose supports share s = 0, 1, 2, 3 points
    for s in range(4):
        other = whole_cell(tuple(range(4 - s, 8 - s)), d)
        c[f"share{s}"] = pair_hist(base, other)

    # axis vs whole cell (axis coordinate in / not in the support)
    axis_in = [0] * d
    axis_in[0] = 1
    axis_out = [0] * d
    axis_out[11] = 1
    for name, ax in (("axis_in", axis_in), ("axis_out", axis_out)):
        h = Counter()
        for b in base:
            h[Fr(abs(dot(ax, b)), 2)] += 1  # norms 1 * 4 => sqrt = 2
        c[name] = dict(h)

    # two distinct axes
    c["axis_axis"] = {Fr(0): 1}
    return c


# ----------------------- STEPS 3-5: the counting proof ---------------------

def min_sum_choose2(slots: int, cap: int, total: int) -> int:
    """Exact minimum of sum C(v_k,2) over integer vectors of length `slots`
    with 0 <= v_k <= cap and sum v_k = total, by exhaustive DP."""
    INF = float("inf")
    best = {0: 0}
    for _ in range(slots):
        nxt = {}
        for ssum, val in best.items():
            for v in range(cap + 1):
                key = ssum + v
                if key > total:
                    break
                cand = val + v * (v - 1) // 2
                if cand < nxt.get(key, INF):
                    nxt[key] = cand
        best = nxt
    if total not in best:
        fail(f"no vector of {slots} values in [0,{cap}] sums to {total}")
    return best[total]


def counting_proof(n_quarter: int, n_half: int, fam: dict, constants: dict, log=print):
    """Run the counting argument against a target histogram with n_quarter
    pairs at |cos| = 1/4 and n_half pairs at |cos| = 1/2.

    Returns True iff NO member of the A+W4 family (fam) can realize the
    target histogram (theorem PROVEN); False if the chain fails to close.
    """
    npts, nax, nb, bs = (fam["n_points"], fam["n_axes"],
                         fam["n_blocks"], fam["block_size"])
    if not (npts == nax == 12 and nb == 51 and bs == 4):
        fail(f"unsupported family parameters {fam!r}; this checker certifies "
             "the 12-axes + 51-whole-W4-cells family")
    if nax + 8 * nb != 420:
        fail(f"family size {nax} + 8*{nb} != 420; cannot be isometric to a 420-line witness")

    # legality: share-3 cells produce |cos| = 3/4 > 1/2 (derived in STEP 2)
    illegal = [v for v in constants["share3"] if v > HALF]
    if not illegal:
        fail("internal: share-3 computation found no |cos| > 1/2; "
             "the 'blocks share <= 2 points' premise would be underived")
    log(f"  legality (derived)      : share-3 cells realize |cos| = "
        f"{max(illegal)} > 1/2 => any 60-degree family has P3 = 0")

    sum_r = bs * nb                      # each block covers bs points
    incidences = sum_r                   # all 12 axes present (nax == npts)
    axis_pairs = nax * (nax - 1) // 2

    # assemble each histogram bin from the DERIVED constants
    def bin_equation(value):
        known = (constants["within"].get(value, 0) * nb
                 + constants["axis_in"].get(value, 0) * incidences
                 + constants["axis_out"].get(value, 0) * 0  # count unknown, coeff checked below
                 + constants["axis_axis"].get(value, 0) * axis_pairs)
        if constants["axis_out"].get(value, 0):
            fail(f"axis-out cells contribute to bin {value}; formula would be incomplete")
        coeffs = {s: constants[f"share{s}"].get(value, 0) for s in range(3)}
        return known, coeffs

    known14, co14 = bin_equation(QUARTER)
    if known14 != 0 or co14[0] != 0 or co14[2] != 0 or co14[1] == 0:
        fail(f"derived 1/4-bin structure unexpected: known={known14}, coeffs={co14}")
    if n_quarter % co14[1] != 0 or n_quarter < 0:
        log(f"  STEP 3: n_1/4 = {n_quarter} is not a multiple of {co14[1]} "
            "=> no family member realizes the histogram.  PROVEN outright.")
        return True
    P1 = n_quarter // co14[1]
    log(f"  STEP 3: n_1/4 = {co14[1]}*P1        => P1 = {P1}")

    known12, co12 = bin_equation(HALF)
    if co12[0] != 0 or co12[1] != 0 or co12[2] == 0:
        fail(f"derived 1/2-bin structure unexpected: coeffs={co12}")
    rem = n_half - known12
    if rem < 0 or rem % co12[2] != 0:
        log(f"  STEP 3: n_1/2 - {known12} = {rem} is not a nonnegative multiple "
            f"of {co12[2]} => no family member realizes the histogram.  PROVEN outright.")
        return True
    P2 = rem // co12[2]
    log(f"  STEP 3: n_1/2 = {known12} + {co12[2]}*P2 => P2 = {P2}")

    total_pairs = nb * (nb - 1) // 2
    if P1 + P2 > total_pairs:
        log(f"  STEP 3: P1 + P2 = {P1+P2} > C({nb},2) = {total_pairs}.  PROVEN outright.")
        return True

    # STEP 4: sum_i C(r_i,2) = P1 + 2*P2 (each share-1 pair meets in 1 point,
    # each share-2 pair in 2, share-3 forbidden); convexity identity over Q
    S = P1 + 2 * P2
    mean = Fr(sum_r, npts)
    # sum C(r_i,2) = (sum r_i^2 - sum r_i)/2 and
    # sum r_i^2 = sum (r_i - mean)^2 + mean^2 * npts  (since sum r_i = mean*npts)
    base_min = (mean * mean * npts - sum_r) / 2
    if base_min.denominator != 1:
        fail("internal: convexity base minimum is not an integer")
    base_min = int(base_min)
    log(f"  STEP 4: sum_i C(r_i,2) = P1 + 2*P2 = {S}; convexity minimum over "
        f"sum r_i = {sum_r} is {base_min} (at r_i = {mean})")
    if S < base_min:
        log(f"  STEP 4: {S} < {base_min} is impossible.  PROVEN outright.")
        return True
    if S > base_min or mean.denominator != 1:
        log(f"  STEP 4: forcing fails ({S} > {base_min}: r_i not forced equal). "
            "NOVELTY NOT ESTABLISHED by this argument.")
        return False
    r = int(mean)
    log(f"  STEP 4: equality => r_i = {r} FORCED for every point")

    # STEP 5: lambda_p <= cap; per point i the 11 values lambda_ij sum to 3*r
    cap = (npts - 2) // (bs - 2)  # extras pairwise disjoint among npts-2 points
    per_point_sum = (bs - 1) * r  # each block on i holds bs-1 pairs through i
    m = min_sum_choose2(npts - 1, cap, per_point_sum)
    log(f"  STEP 5: lambda_p <= {cap}; per point sum_j lambda_ij = {per_point_sum} "
        f"over {npts-1} values; DP minimum of sum_j C(lambda_ij,2) = {m}")
    lower = -(-npts * m // 2)  # ceil
    log(f"  STEP 5: 2*P2 = sum_i sum_j C(lambda_ij,2) >= {npts}*{m} = {npts*m} "
        f"=> P2 >= {lower}")
    if P2 < lower:
        log(f"  CONTRADICTION: STEP 3 forced P2 = {P2} < {lower}.  "
            "No family member realizes the histogram.")
        return True
    log(f"  STEP 5: P2 = {P2} >= {lower}; no contradiction. "
        "NOVELTY NOT ESTABLISHED by this argument.")
    return False


def verify(nov_path, wit_path) -> int:
    fam_cert = json.loads(Path(nov_path).read_text())
    fam = fam_cert.get("family")
    if not isinstance(fam, dict) or set(fam) != {"n_points", "n_axes", "n_blocks", "block_size"}:
        fail(f"novelty certificate must store the four family parameters, got {fam!r}")

    print(f"witness : {wit_path}")
    print(f"novelty : {nov_path}  (family: {fam['n_axes']} axes + "
          f"{fam['n_blocks']} whole W4 cells over {fam['n_points']} points)")

    d, n, hist = witness_histogram(wit_path)
    if n != 420:
        fail(f"witness has {n} lines, not 420; the counting theorem targets 420")
    hist_str = ", ".join(f"|cos|={k}: {v}" for k, v in sorted(hist.items()))
    print(f"  STEP 1: recomputed witness histogram {{{hist_str}}}")

    # packing precondition (hard gate): a 60-degree witness realizes only
    # |cos| <= 1/2.  Any bin above 1/2 means the input is not a kissing
    # configuration at all, so no novelty verdict applies -- reject outright
    # rather than letting an invalid input reach the counting argument.
    over = {k: v for k, v in hist.items() if k > HALF and v}
    if over:
        fail(f"witness has {sum(over.values())} pair(s) at |cos| > 1/2 "
             f"({sorted(str(k) for k in over)}); it violates the 60-degree "
             "packing bound, so it is not a kissing configuration and no "
             "novelty verdict applies")

    constants = derive_constants(d)
    print(f"  STEP 2: derived constants — within-cell {fmt(constants['within'])}, "
          f"share0 {fmt(constants['share0'])}, share1 {fmt(constants['share1'])}, "
          f"share2 {fmt(constants['share2'])}, share3 {fmt(constants['share3'])}, "
          f"axis-in {fmt(constants['axis_in'])}, axis-out {fmt(constants['axis_out'])}")

    # quick invariant: LEGAL members of the A+W4 family only realize |cos|
    # values <= 1/2 on the {0, 1/4, 1/2} grid.  The > 1/2 entries of the
    # derived tables (e.g. share-3's 3/4) belong to ILLEGAL family members
    # and must not count as realizable values (they are also excluded by the
    # packing precondition above).
    legal_bins = set()
    for tbl in constants.values():
        legal_bins |= {v for v in tbl if v <= HALF}
    off_grid = {k: v for k, v in hist.items() if k not in legal_bins and v}
    if off_grid:
        print(f"  witness has {sum(off_grid.values())} pair(s) at |cos| values "
              f"{sorted(str(k) for k in off_grid)} that NO A+W4 family member "
              "realizes.  PROVEN outright (quick invariant).")
        print("  VERIFIED: witness is NOT isometric to any member of the family.")
        return 0
    print("  quick invariant silent   : witness |cos| support lies inside the "
          "family's {0, 1/4, 1/2} grid — the counting argument is required")

    proven = counting_proof(hist.get(QUARTER, 0), hist.get(HALF, 0), fam, constants)
    if not proven:
        fail("the counting argument did not close; novelty is NOT established")
    print("  VERIFIED: no 12-axes + 51-whole-W4-cell configuration (including the "
          "known modular 420) realizes the witness histogram; since the pairwise "
          "|cos| multiset is an isometry invariant, the witness is PROVEN "
          "NON-ISOMETRIC to the known modular 420.")
    return 0


def fmt(tbl) -> str:
    return "{" + ", ".join(f"{k}:{v}" for k, v in sorted(tbl.items())) + "}"


def main(argv: list[str]) -> int:
    nov = argv[1] if len(argv) > 1 else str(CERTS / "novelty.json")
    wit = argv[2] if len(argv) > 2 else str(CERTS / "config420.json")
    try:
        return verify(nov, wit)
    except CertificateError as e:
        print(e)
        return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
