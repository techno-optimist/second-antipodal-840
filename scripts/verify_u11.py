#!/usr/bin/env python3
"""Exact verification of the U11 witnesses (v1.1): alpha(U11) >= 98.

CLAIM CERTIFIED (the v1.1 postscript's answer to Problem 26).  The U11 world
has ground set [0..10]; a cell is a 3-subset A (defect), its support is the
8-set H = [0..10] \\ A; a line of cell A is a vector g in Z^11 with entries
+-1 on H and 0 on A, taken in the even sign class (even number of -1) modulo
global flip (unit representative g/sqrt(8), norm^2 = 8).  Two lines conflict
iff |cos| > 1/2; since <g,g'> is the SIGNED SHARED SUM over the support
overlap R = H cap H' (each shared coordinate contributes s_i * s'_i), the
note's exact predicate becomes

    conflict  <=>  4*<g,g'>^2 > 8*8 = 64  <=>  |<g,g'>| >= 5.

A 98-line independent set therefore proves alpha(U11) >= 98 and answers the
note's Problem 26 (is alpha(U11) <= 96?) IN THE NEGATIVE.

WHAT THIS SCRIPT CHECKS (stdlib only; exact integer arithmetic; no floats on
any decision path).  The certificates store ONLY cell + sign descriptions —
no counts, no sizes, no verdicts.  For every witness in every given file:

  1. rebuilds each line's integer representative g in Z^11 from its cell
     3-subset and canonical even sign pattern (signs[0] = +1), recomputing
     norm^2 = 8, and validates pairwise distinctness as lines;
  2. verifies the packing inequality 4*<g,g'>^2 <= 64 for ALL unordered
     pairs — exactly |cos| <= 1/2, i.e. zero conflicts, an independent set;
  3. RECOMPUTES (never echoes): the witness size, the cell census (distinct
     cells, lines-per-cell histogram) and the support-overlap census of the
     conflict-checked pairs (r = |H cap H'| in {5,6,7} across distinct cells,
     8 within a cell — no free cell pairs).

Exit status 0 iff every witness verifies AND the largest recomputed witness
has at least 98 lines (the v1.1 claim; >= 97 already refutes Problem 26).

Usage:  python3 scripts/verify_u11.py [cert.json ...]
        (default: certs/u11_witness98.json certs/u11_witnesses97.json)
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from witness import CERTS, CertificateError, dot, fail  # noqa: E402

D = 11
NORM2 = 8
CLAIM = 98  # the v1.1 postscript claims alpha(U11) >= 98


def load_u11_witness(name: str, recs) -> list:
    """Structurally validate one witness; return lines as (g, cell) tuples."""
    if not isinstance(recs, list):
        fail(f"witness {name}: 'lines' must be a list, got {type(recs).__name__}")
    lines = []
    for idx, rec in enumerate(recs):
        cell, signs = rec.get("cell"), rec.get("signs")
        if (not isinstance(cell, list) or len(cell) != 3
                or any(not isinstance(i, int) or not 0 <= i < D for i in cell)
                or sorted(set(cell)) != cell):
            fail(f"witness {name}, line {idx}: cell must be 3 distinct sorted "
                 f"coordinates in [0,{D}), got {cell!r}")
        if (not isinstance(signs, list) or len(signs) != 8
                or any(s not in (1, -1) for s in signs)):
            fail(f"witness {name}, line {idx}: signs must be 8 values +-1, "
                 f"got {signs!r}")
        if signs[0] != 1:
            fail(f"witness {name}, line {idx}: non-canonical U11 line "
                 "(signs[0] must be +1)")
        if sum(1 for s in signs if s == -1) % 2 != 0:
            fail(f"witness {name}, line {idx}: U11 line outside the even sign "
                 f"class ({signs.count(-1)} minus signs)")
        support = [i for i in range(D) if i not in cell]  # ascending, |H| = 8
        g = [0] * D
        for i, s in zip(support, signs):
            g[i] = s
        gt = tuple(g)
        if sum(x * x for x in gt) != NORM2:
            fail(f"witness {name}, line {idx}: recomputed norm^2 != {NORM2}")
        lines.append((gt, tuple(cell)))
    seen = {}
    for idx, (gt, _) in enumerate(lines):
        if gt in seen:
            fail(f"witness {name}: duplicate line, index {idx} repeats "
                 f"index {seen[gt]}")
        seen[gt] = idx
    return lines


def verify_witness(name: str, recs) -> int:
    """Verify one witness end to end; return its recomputed size."""
    lines = load_u11_witness(name, recs)
    n = len(lines)
    pairs = violations = zero_margin = 0
    overlap: Counter = Counter()
    for i in range(n):
        gi, ci = lines[i]
        for j in range(i + 1, n):
            gj, cj = lines[j]
            pairs += 1
            r = sum(1 for a, b in zip(gi, gj) if a and b)
            overlap[r] += 1
            dp = dot(gi, gj)  # the signed shared sum over the r-overlap
            lhs = 4 * dp * dp
            if lhs > NORM2 * NORM2:
                violations += 1
                if violations <= 3:
                    print(f"  VIOLATION: witness {name} pair ({i},{j}) has "
                          f"4<g,g'>^2 = {lhs} > 64  (|signed shared sum| = "
                          f"{abs(dp)} >= 5 on the r={r} overlap)")
            elif lhs == NORM2 * NORM2:
                zero_margin += 1
    if pairs != n * (n - 1) // 2:
        fail(f"witness {name}: internal pair-count mismatch")
    if violations:
        fail(f"witness {name}: {violations} pair(s) violate the 60-degree "
             f"condition 4<g,g'>^2 <= 64")
    cells = Counter(c for _, c in lines)
    sizes = Counter(cells.values())
    bad_r = sorted(set(overlap) - {5, 6, 7, 8})
    if bad_r:
        fail(f"witness {name}: impossible support overlaps {bad_r}")
    print(f"  {name}: recomputed size {n}; pairs checked {pairs} "
          f"(0 violations, {zero_margin} zero-margin); "
          f"cells {len(cells)} distinct, lines-per-cell "
          f"{dict(sorted(sizes.items()))}; overlap census "
          f"{ {r: overlap[r] for r in sorted(overlap)} }")
    return n


def main(argv: list[str]) -> int:
    paths = argv[1:] or [str(CERTS / "u11_witness98.json"),
                         str(CERTS / "u11_witnesses97.json")]
    best = 0
    try:
        for path in paths:
            cert = json.loads(Path(path).read_text())
            if not isinstance(cert, dict) or cert.get("world") != "U11":
                fail(f"{path}: not a U11 witness certificate "
                     "(missing 'world': 'U11')")
            print(f"witness file: {path}")
            if "lines" in cert:
                best = max(best, verify_witness(Path(path).stem, cert["lines"]))
            elif "witnesses" in cert:
                if not isinstance(cert["witnesses"], dict) or not cert["witnesses"]:
                    fail(f"{path}: 'witnesses' must be a nonempty dict")
                for wname in sorted(cert["witnesses"]):
                    best = max(best, verify_witness(wname, cert["witnesses"][wname]))
            else:
                fail(f"{path}: certificate has neither 'lines' nor 'witnesses'")
        if best < CLAIM:
            fail(f"largest verified witness has {best} lines; the v1.1 claim "
                 f"needs a {CLAIM}-line witness"
                 + (" (alpha(U11) <= 96 is still refuted)" if best >= 97 else ""))
    except CertificateError as e:
        print(e)
        return 1
    print(f"  VERIFIED: independent sets of up to {best} distinct U11 lines, "
          f"all pairs at |cos| <= 1/2 exactly => alpha(U11) >= {best}.")
    print("  Problem 26 of the note (is alpha(U11) <= 96?) is ANSWERED IN "
          "THE NEGATIVE by explicit witness.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
