#!/usr/bin/env python3
"""Exact verification of the 420-line antipodal kissing witness in R^12.

THEOREM CERTIFIED.  There exist 420 distinct lines through the origin of R^12
whose unit representatives u satisfy |<u_i, u_j>| <= 1/2 for all i != j.
Doubling each line to +-u gives an 840-point ANTIPODAL kissing configuration
(pairwise angles >= 60 degrees), i.e. N60(12) >= 420 lines / 840 points.

WHAT THIS SCRIPT CHECKS (all arithmetic exact integers/Fractions; no floats
on any decision path).  The certificate certs/config420.json stores ONLY the
line descriptions — no counts, no histogram, no verdicts.  This script:

  1. rebuilds all 420 integer representatives g (W2: e_a +- e_b, norm^2 = 2;
     U: +-1 pattern on an 8-coordinate support, even sign class, norm^2 = 8)
     and validates canonical form + pairwise distinctness as lines;
  2. verifies the packing inequality 4*<g,g'>^2 <= N*N' for ALL C(420,2)
     = 87,990 unordered pairs — exactly equivalent to |cos angle| <= 1/2;
  3. RECOMPUTES (never echoes): the line count and composition, the number
     of zero-margin pairs (equality 4*<g,g'>^2 = N*N', i.e. touching pairs
     at exactly 60 degrees), the exact |<u,u'>| histogram, and the U-cell
     structure (distinct 8-supports and lines per cell).

Displayed decimals (none are load-bearing) are rounded UP via directed
rounding.  Exit status 0 iff the witness verifies.

Usage:  python3 scripts/verify_420.py [certs/config420.json]
"""
from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from witness import CERTS, CertificateError, dot, exact_abs_cos, fail, load_witness

N_CLAIMED = 420          # the theorem target: this many lines must be present
COMPOSITION = {"W2": 132, "U": 288}


def verify(path) -> int:
    d, lines = load_witness(path)
    n = len(lines)
    comp = Counter(k for k, _, _ in lines)
    print(f"witness: {path}")
    print(f"  recomputed line count   : {n}  (composition {dict(comp)})")
    if n != N_CLAIMED or dict(comp) != COMPOSITION:
        fail(f"expected exactly {N_CLAIMED} lines with composition {COMPOSITION}, "
             f"recomputed {n} with {dict(comp)}")

    # --- the packing check: every pair at |cos| <= 1/2, exactly ---
    pairs = 0
    violations = 0
    zero_margin = 0
    hist: Counter = Counter()
    for i in range(n):
        _, gi, Ni = lines[i]
        for j in range(i + 1, n):
            _, gj, Nj = lines[j]
            dp = dot(gi, gj)
            pairs += 1
            lhs, rhs = 4 * dp * dp, Ni * Nj
            if lhs > rhs:
                violations += 1
                if violations <= 3:
                    print(f"  VIOLATION: pair ({i},{j}) has 4<g,g'>^2 = {lhs} > "
                          f"N*N' = {rhs}  (|cos| = {exact_abs_cos(dp, Ni, Nj)} > 1/2)")
            elif lhs == rhs:
                zero_margin += 1
            hist[exact_abs_cos(dp, Ni, Nj)] += 1
    print(f"  pairs checked           : {pairs}")
    if pairs != n * (n - 1) // 2:
        fail(f"internal pair-count mismatch: {pairs} != C({n},2)")
    if violations:
        fail(f"{violations} pair(s) violate the 60-degree condition "
             "4<g,g'>^2 <= N*N'")
    print(f"  packing inequality      : 4<g,g'>^2 <= N*N' holds for all "
          f"{pairs} pairs (0 violations)")
    print(f"  zero-margin pairs       : {zero_margin}  (touching at exactly 60 degrees)")
    hist_str = ", ".join(f"|cos|={k}: {v}" for k, v in sorted(hist.items()))
    print(f"  |cos| histogram         : {{{hist_str}}}")
    if sum(hist.values()) != pairs:
        fail("histogram total does not match the pair count")

    # --- recomputed U-cell structure ---
    cells = Counter()
    for kind, g, _ in lines:
        if kind == "U":
            cells[tuple(i for i, v in enumerate(g) if v)] += 1
    sizes = Counter(cells.values())
    print(f"  U cells (8-supports)    : {len(cells)} distinct; "
          f"lines-per-cell histogram {dict(sorted(sizes.items()))}")

    print(f"  VERIFIED: {n} distinct lines in R^{d}, pairwise |cos| <= 1/2 "
          f"=> antipodal kissing configuration with {2*n} points; N60(12) >= {n}.")
    return 0


def main(argv: list[str]) -> int:
    path = argv[1] if len(argv) > 1 else str(CERTS / "config420.json")
    try:
        return verify(path)
    except CertificateError as e:
        print(e)
        return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
