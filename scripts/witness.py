#!/usr/bin/env python3
"""Shared exact-witness loader for the d12 420-line note (stdlib only).

Rebuilds every line of certs/config420.json as an INTEGER representative
g in Z^12 with squared norm N:

    W2 {a,b,sign} :  g = e_a + sign*e_b,              N = 2  (unit g/sqrt(2))
    U  {support,signs}:  g_i = signs[k] on support,   N = 8  (unit g/sqrt(8))

so that for two lines the cosine of their angle is

    <u,u'> = <g,g'> / sqrt(N*N'),   sqrt(N*N') in {2, 4, 8} (always an integer),

and every geometric decision below reduces to EXACT integer arithmetic:
    |<u,u'>| <= 1/2   <=>   4*<g,g'>^2 <= N*N'.

No floats anywhere. This module performs only STRUCTURAL validation and the
rebuild; all theorem-level checks live in the verify_* scripts.
"""
from __future__ import annotations

import json
from fractions import Fraction as Fr
from math import isqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent
CERTS = HERE.parent / "certs"


class CertificateError(Exception):
    """A certificate failed validation. The message is the diagnostic."""


def fail(msg: str) -> None:
    raise CertificateError("FAIL: " + msg)


def load_witness(path: Path | str):
    """Load and structurally validate the witness; return (d, lines) where
    lines is a list of (kind, g, N) with g a tuple of ints and N = <g,g>.

    Canonical form enforced (so tuple-distinctness == line-distinctness):
    the first nonzero entry of every g is +1.
    """
    cert = json.loads(Path(path).read_text())
    if not isinstance(cert, dict) or "lines" not in cert or "dimension" not in cert:
        fail("witness missing 'dimension' or 'lines'")
    d = cert["dimension"]
    if d != 12:
        fail(f"witness dimension is {d!r}, this note certifies d = 12")
    lines = []
    for idx, rec in enumerate(cert["lines"]):
        kind = rec.get("type")
        g = [0] * d
        if kind == "W2":
            a, b, s = rec.get("a"), rec.get("b"), rec.get("sign")
            if not (isinstance(a, int) and isinstance(b, int) and 0 <= a < b < d):
                fail(f"line {idx}: W2 needs integer 0 <= a < b < {d}, got a={a!r} b={b!r}")
            if s not in (1, -1):
                fail(f"line {idx}: W2 sign must be +-1, got {s!r}")
            g[a], g[b] = 1, s
        elif kind == "U":
            sup, signs = rec.get("support"), rec.get("signs")
            if (not isinstance(sup, list) or len(sup) != 8
                    or any(not isinstance(i, int) or not 0 <= i < d for i in sup)
                    or sorted(set(sup)) != sup):
                fail(f"line {idx}: U support must be 8 distinct sorted coords in [0,{d}), got {sup!r}")
            if (not isinstance(signs, list) or len(signs) != 8
                    or any(s not in (1, -1) for s in signs)):
                fail(f"line {idx}: U signs must be 8 values +-1, got {signs!r}")
            if signs[0] != 1:
                fail(f"line {idx}: non-canonical U line (signs[0] must be +1)")
            if sum(1 for s in signs if s == -1) % 2 != 0:
                fail(f"line {idx}: U line outside the even sign class "
                     f"({signs.count(-1)} minus signs)")
            for i, s in zip(sup, signs):
                g[i] = s
        else:
            fail(f"line {idx}: unknown line type {kind!r}")
        gt = tuple(g)
        N = sum(x * x for x in gt)
        lines.append((kind, gt, N))
    seen = {}
    for idx, (_, gt, _) in enumerate(lines):
        if gt in seen:
            fail(f"duplicate line: index {idx} repeats index {seen[gt]}")
        seen[gt] = idx
    return d, lines


def dot(g: tuple, h: tuple) -> int:
    return sum(a * b for a, b in zip(g, h))


def exact_abs_cos(dp: int, N1: int, N2: int) -> Fr:
    """|<u,u'>| as an exact Fraction; requires N1*N2 to be a perfect square."""
    s = isqrt(N1 * N2)
    if s * s != N1 * N2:
        fail(f"norm product {N1}*{N2} is not a perfect square; "
             "exact cosine bookkeeping does not apply")
    return Fr(abs(dp), s)


def ceil_str(x: Fr, places: int = 6) -> str:
    """Directed (upward) decimal display of a nonnegative rational — the
    conservative direction when displaying angles/cosines that must stay
    below a threshold is context-dependent, so callers choose; this rounds UP."""
    scale = 10 ** places
    n = -((-x.numerator * scale) // x.denominator)  # ceil
    return f"{n // scale}.{n % scale:0{places}d}"
