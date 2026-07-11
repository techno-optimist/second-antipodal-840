#!/usr/bin/env python3
"""Exact Q(sqrt2) verification of the rotation certificate.

THEOREM CERTIFIED.  Let P be the perfect pairing of the 12 coordinates stored
in certs/rotation.json and T_P the blockwise 45-degree rotation that acts on
each pair (i, j) of P by

    (x_i, x_j)  |->  ( (x_i - x_j)/sqrt(2),  (x_i + x_j)/sqrt(2) ).

T_P is orthogonal, and it maps the 420-line witness of certs/config420.json
ISOMETRICALLY into the INTEGER "A union W4" world of R^12:

    A  (axis lines)  +-e_k,
    W4 (quad lines)  entries +-1/2 on a 4-coordinate support.

Consequence for the note: up to isometry the witness lives inside the integer
A+W4 world, whose maximum 60-degree line family is exactly 420 (Hanani cap),
so the witness is capped in its slice and cannot be extended to 421 there.

EXACTNESS.  sqrt(2) never appears on a decision path: the script works with
h := sqrt(2) * T_P(g), which is the INTEGER vector with h_i = g_i - g_j and
h_j = g_i + g_j on each pair (i, j).  Then

    <T_P g, T_P g'> = <h, h'> / 2.

NOTE that <h,h'> = 2<g,g'> (and <h,h> = 2<g,g>) is an ALGEBRAIC IDENTITY of
the integer map, valid for EVERY input and pairing — the isometry re-check
below is arithmetic self-validation, not a falsifiable gate.  The
falsifiable gates are:

  (i)   the exact 60-degree predicate 4<h,h'>^2 <= (2N)(2N') on ALL image
        pairs — by the identity this is equivalent to the WITNESS itself
        being a valid kissing configuration, so this checker does not
        certify anything about geometry-violating inputs;
  (ii)  image membership in the A union W4 world, checked on the exact
        rational unit vector h / sqrt(2N)  (sqrt(2N) in {2,4}, an integer);
  (iii) the whole-cell (8 lines) vs half-cell (4 lines) census;
  (iv)  image-line injectivity.

RECOMPUTED, never echoed: the image composition (axis / W4 counts), the
distinct W4 image supports, and the whole/half decomposition.  Exit status
0 iff everything verifies.

Usage:  python3 scripts/verify_rotation.py [certs/rotation.json [certs/config420.json]]
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from fractions import Fraction as Fr
from math import isqrt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from witness import CERTS, CertificateError, dot, fail, load_witness


def load_pairing(path, d: int):
    cert = json.loads(Path(path).read_text())
    pairing = cert.get("pairing")
    if (not isinstance(pairing, list) or len(pairing) != d // 2
            or any(not isinstance(p, list) or len(p) != 2 for p in pairing)):
        fail(f"rotation certificate must store {d//2} coordinate pairs, got {pairing!r}")
    flat = [i for p in pairing for i in p]
    if sorted(flat) != list(range(d)):
        fail(f"pairing is not a perfect pairing of 0..{d-1}: {pairing!r}")
    return [tuple(p) for p in pairing]


def rotate_scaled(g: tuple, pairing) -> tuple:
    """Return h = sqrt(2) * T_P(g) as an integer vector."""
    h = [0] * len(g)
    for i, j in pairing:
        h[i] = g[i] - g[j]
        h[j] = g[i] + g[j]
    return tuple(h)


def classify_image(h: tuple, N: int):
    """Classify the unit image T_P(g) = h / sqrt(2N) as 'axis' or 'W4'.

    Requires every unit entry to lie in {0, +-1} (axis) or {0, +-1/2} (W4)
    with the right support size; anything else is outside the A+W4 world.
    """
    s = isqrt(2 * N)
    if s * s != 2 * N:
        fail(f"2N = {2*N} is not a perfect square; cannot form the exact unit image")
    entries = [Fr(x, s) for x in h]
    support = [k for k, e in enumerate(entries) if e != 0]
    vals = {abs(e) for e in entries if e != 0}
    if len(support) == 1 and vals == {Fr(1)}:
        return "axis", tuple(support)
    if len(support) == 4 and vals == {Fr(1, 2)}:
        return "W4", tuple(support)
    fail(f"rotated image h/{s} = {[str(e) for e in entries]} is not in the "
         f"A-union-W4 world (support size {len(support)}, |entries| "
         f"{sorted(str(v) for v in vals)})")


def verify(rot_path, wit_path) -> int:
    d, lines = load_witness(wit_path)
    pairing = load_pairing(rot_path, d)
    print(f"witness : {wit_path}")
    print(f"rotation: {rot_path}")
    print(f"  pairing                 : {pairing}")

    scaled = [(rotate_scaled(g, pairing), N) for _, g, N in lines]

    # --- exact isometry: norms and ALL pairwise Grams preserved ---
    for idx, (h, N) in enumerate(scaled):
        if sum(x * x for x in h) != 2 * N:
            fail(f"image {idx}: <h,h> = {sum(x*x for x in h)} != 2N = {2*N}; "
                 "T_P failed to preserve the norm")
    n = len(lines)
    checked = 0
    for i in range(n):
        gi, (hi, Ni) = lines[i][1], scaled[i]
        for j in range(i + 1, n):
            gj, (hj, Nj) = lines[j][1], scaled[j]
            dh = dot(hi, hj)
            if dh != 2 * dot(gi, gj):
                fail(f"pair ({i},{j}): <h,h'> = {dh} != 2<g,g'> = "
                     f"{2*dot(gi, gj)}; T_P is not an isometry on the witness")
            # REAL GATE: the images (norms 2N, 2N') must satisfy the exact
            # 60-degree bound 4<h,h'>^2 <= (2N)(2N'), i.e. <h,h'>^2 <= N*N'.
            # By the identity <h,h'> = 2<g,g'> this is equivalent to the
            # witness itself being a valid kissing configuration.
            if dh * dh > Ni * Nj:
                fail(f"pair ({i},{j}): rotated images violate the 60-degree "
                     f"bound (<h,h'>^2 = {dh*dh} > N*N' = {Ni*Nj}); the "
                     "input is not a valid kissing configuration, so no "
                     "slice-cap conclusion applies")
            checked += 1
    print(f"  isometry (identity)     : all {n} norms and {checked} pairwise "
          f"Grams preserved exactly (<h,h'> = 2<g,g'>, arithmetic re-check "
          f"of the map identity)")
    print(f"  packing (real gate)     : all {checked} image pairs satisfy "
          f"the exact 60-degree bound <h,h'>^2 <= N*N'")

    # --- image membership + recomputed composition ---
    comp = Counter()
    w4_cells = Counter()
    images = set()
    for h, N in scaled:
        kind, support = classify_image(h, N)
        comp[kind] += 1
        if kind == "W4":
            w4_cells[support] += 1
        # canonical image line rep: first nonzero positive
        first = next(x for x in h if x)
        canon = h if first > 0 else tuple(-x for x in h)
        if canon in images:
            fail("two witness lines map to the same image line; T_P not injective on lines")
        images.add(canon)
    sizes = Counter(w4_cells.values())
    whole = sizes.get(8, 0)
    half = sizes.get(4, 0)
    print(f"  image composition       : {dict(comp)} (all 420 in the integer A+W4 world)")
    print(f"  distinct image supports : {len(w4_cells)} W4 supports + "
          f"{comp['axis']} axes")
    print(f"  W4 cell decomposition   : {whole} whole cells (8 lines) + "
          f"{half} half cells (4 lines); size histogram {dict(sorted(sizes.items()))}")
    if whole * 8 + half * 4 != comp["W4"] or set(sizes) - {4, 8}:
        fail(f"W4 images do not decompose into whole/half cells: {dict(sizes)}")
    print(f"  VERIFIED: T_P is an exact isometry carrying the witness into the "
          f"integer A+W4 world ({comp['axis']} axes + {comp['W4']} W4 lines) "
          f"=> the witness is Hanani-capped in its slice.")
    return 0


def main(argv: list[str]) -> int:
    rot = argv[1] if len(argv) > 1 else str(CERTS / "rotation.json")
    wit = argv[2] if len(argv) > 2 else str(CERTS / "config420.json")
    try:
        return verify(rot, wit)
    except CertificateError as e:
        print(e)
        return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
