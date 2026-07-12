#!/usr/bin/env python3
"""Adversarial self-test for the d12 420-line verification suite.

A verifier that accepts corrupted certificates is a notary for lies, so this
script feeds each checker deliberately-broken variants of the shipped
certificates and PASSES only if the checker REJECTS each one with the
expected diagnostic substring (and still accepts the originals).

Battery (checker <- corruption -> required diagnostic):

  verify_420.py
    A1 pristine witness                        -> accepted, "VERIFIED"
    A2 two U signs flipped (parity-preserving,
       geometry-breaking: |<g,g'>| = 6)        -> "VIOLATION"/"violate"
    A3 one U sign flipped (parity-breaking)    -> "even sign class"
    A4 duplicated line                         -> "duplicate line"
    A5 deleted line (419 remain)               -> "expected exactly 420"
    A6 W2 with a == b                          -> "W2 needs integer"
    A7 W2 sign = 0                             -> "sign must be +-1"
    A8 U support with a repeated coordinate    -> "8 distinct sorted"
    A9 non-canonical U line (all signs flipped;
       same line, signs[0] = -1)               -> "non-canonical"
    A10 W2 sign flipped onto an existing line  -> "duplicate line"

  verify_rotation.py
    B1 pristine certificate                    -> accepted, "VERIFIED"
    B2 repeated coordinate in the pairing      -> "not a perfect pairing"
    B3 WRONG-but-perfect pairing ([[0,1],...])
       (soundness: isometry alone must not
       suffice; image must land in A+W4)       -> "not in the A-union-W4 world"
    B4 five pairs instead of six               -> "coordinate pairs"
    B5 geometry-breaking witness under a valid
       pairing (designed coverage: the image
       packing gate must fire on the exact
       60-degree predicate)                    -> "violate"

  verify_novelty.py
    C1 pristine certificate                    -> accepted, "PROVEN NON-ISOMETRIC"
    C2 SOUNDNESS (load-bearing): the KNOWN modular histogram
       (n_1/4 = 32256, n_1/2 = 20496) must NOT be "proven" non-realizable —
       counting_proof must return False (P2 = 564 meets the bound exactly)
    C3 family with n_blocks = 50               -> "unsupported family parameters"
    C4 family parameters missing               -> "four family parameters"
    C5 corrupted witness (parity-breaking)     -> witness diagnostic propagates
    C6 geometry-breaking witness (parity-preserving, |cos| = 3/4 pairs):
       the packing precondition must reject it
       BEFORE any novelty verdict              -> "not a kissing configuration"

  verify_u11.py (v1.1)
    D1 pristine u11 witnesses                  -> accepted, "VERIFIED" + the
                                                  Problem-26 negative verdict
    D2 one U11 sign flipped (parity-breaking)  -> "even sign class"
    D3 two U11 signs flipped (parity-preserving,
       geometry-breaking; found by determinis-
       tic search, must conflict with another
       witness line)                           -> "VIOLATION"/"violate"
    D4 duplicated U11 line                     -> "duplicate line"
    D5 deleted U11 line (97 remain: alpha <= 96
       still refuted, but the v1.1 claim gate
       alpha(U11) >= 98 must fire)             -> "needs a 98-line witness"

C2 exists specifically to catch "always shout contradiction" regressions of
the counting chain: a checker that proves novelty of the modular histogram
itself is unsound.  Every corruption must exit nonzero AND print the expected
diagnostic substring.

Usage: python3 scripts/selftest.py
"""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
CERTS = HERE.parent / "certs"
PY = sys.executable

sys.path.insert(0, str(HERE))
import verify_novelty  # for the C2 soundness check (module-level call)

WITNESS = json.loads((CERTS / "config420.json").read_text())
ROTATION = json.loads((CERTS / "rotation.json").read_text())
NOVELTY = json.loads((CERTS / "novelty.json").read_text())
U11_98 = json.loads((CERTS / "u11_witness98.json").read_text())
FIRST_U = next(i for i, l in enumerate(WITNESS["lines"]) if l["type"] == "U")

passed = failed = 0


def run(script: str, *args: str) -> tuple[int, str]:
    r = subprocess.run([PY, str(HERE / script), *map(str, args)],
                       capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def tmp_json(obj) -> str:
    f = tempfile.NamedTemporaryFile("w", suffix=".json", delete=False)
    json.dump(obj, f)
    f.close()
    return f.name


def check(label: str, ok: bool, detail: str = ""):
    global passed, failed
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}")
    if not ok:
        failed += 1
        if detail:
            print("        " + detail.replace("\n", "\n        ")[:800])
    else:
        passed += 1


def expect_reject(label: str, script: str, args: list[str], needles: tuple[str, ...]):
    code, out = run(script, *args)
    ok = code != 0 and any(n in out for n in needles)
    check(label, ok, f"exit={code}, expected one of {needles!r} in:\n{out}")


def mutated_witness(mutate) -> str:
    w = json.loads(json.dumps(WITNESS))
    mutate(w)
    return tmp_json(w)


def mutated_u11(mutate) -> str:
    w = json.loads(json.dumps(U11_98))
    mutate(w)
    return tmp_json(w)


def u11_vec(rec) -> tuple:
    """Integer representative of a U11 cert line (mirrors verify_u11)."""
    g = [0] * 11
    support = [i for i in range(11) if i not in rec["cell"]]
    for i, s in zip(support, rec["signs"]):
        g[i] = s
    return tuple(g)


def u11_conflicting_double_flip():
    """Deterministic search: first (line index, sign pos pair) whose double
    flip keeps the line canonical + even but conflicts with another witness
    line (and duplicates none), so ONLY the packing gate can catch it."""
    lines = U11_98["lines"]
    vecs = [u11_vec(r) for r in lines]
    for i, rec in enumerate(lines):
        for p in range(1, 8):          # never touch signs[0]: stay canonical
            for q in range(p + 1, 8):
                cand = json.loads(json.dumps(rec))
                cand["signs"][p] *= -1
                cand["signs"][q] *= -1
                v = u11_vec(cand)
                if any(v == w for w in vecs):
                    continue           # would trip the duplicate gate instead
                for j, w in enumerate(vecs):
                    if j == i:
                        continue
                    dp = sum(a * b for a, b in zip(v, w))
                    if 4 * dp * dp > 64:
                        return i, p, q
    raise AssertionError("no conflicting double flip found in the 98-witness")


def main() -> int:
    print("verify_420.py battery")
    code, out = run("verify_420.py")
    check("A1 pristine witness accepted", code == 0 and "VERIFIED" in out, out)

    def a2(w):
        s = w["lines"][FIRST_U]["signs"]
        s[1] *= -1
        s[3] *= -1  # parity preserved; creates |<g,g'>| = 6 with another U line
    expect_reject("A2 double sign flip -> packing violation", "verify_420.py",
                  [mutated_witness(a2)], ("VIOLATION", "violate"))

    def a3(w):
        w["lines"][FIRST_U]["signs"][1] *= -1
    expect_reject("A3 single sign flip -> parity gate", "verify_420.py",
                  [mutated_witness(a3)], ("even sign class",))

    def a4(w):
        w["lines"][-1] = json.loads(json.dumps(w["lines"][0]))
    expect_reject("A4 duplicated line", "verify_420.py",
                  [mutated_witness(a4)], ("duplicate line",))

    def a5(w):
        w["lines"].pop()
    expect_reject("A5 deleted line", "verify_420.py",
                  [mutated_witness(a5)], ("expected exactly 420",))

    def a6(w):
        w["lines"][0]["b"] = w["lines"][0]["a"]
    expect_reject("A6 W2 with a == b", "verify_420.py",
                  [mutated_witness(a6)], ("W2 needs integer",))

    def a7(w):
        w["lines"][0]["sign"] = 0
    expect_reject("A7 W2 sign = 0", "verify_420.py",
                  [mutated_witness(a7)], ("sign must be +-1",))

    def a8(w):
        sup = w["lines"][FIRST_U]["support"]
        sup[1] = sup[0]
    expect_reject("A8 repeated support coordinate", "verify_420.py",
                  [mutated_witness(a8)], ("8 distinct sorted",))

    def a9(w):
        w["lines"][FIRST_U]["signs"] = [-s for s in w["lines"][FIRST_U]["signs"]]
    expect_reject("A9 non-canonical U line", "verify_420.py",
                  [mutated_witness(a9)], ("non-canonical",))

    def a10(w):
        # every W2 sign variant exists, so a sign flip collides with a line
        w["lines"][0]["sign"] *= -1
    expect_reject("A10 W2 sign flip -> duplicate", "verify_420.py",
                  [mutated_witness(a10)], ("duplicate line",))

    print("verify_rotation.py battery")
    code, out = run("verify_rotation.py")
    check("B1 pristine rotation accepted", code == 0 and "VERIFIED" in out, out)

    bad = json.loads(json.dumps(ROTATION))
    bad["pairing"][0] = [0, 0]
    expect_reject("B2 repeated coordinate in pairing", "verify_rotation.py",
                  [tmp_json(bad)], ("not a perfect pairing",))

    bad = json.loads(json.dumps(ROTATION))
    bad["pairing"] = [[0, 1], [2, 3], [4, 5], [6, 7], [8, 9], [10, 11]]
    expect_reject("B3 wrong-but-perfect pairing (soundness)", "verify_rotation.py",
                  [tmp_json(bad)], ("not in the A-union-W4 world",))

    bad = json.loads(json.dumps(ROTATION))
    bad["pairing"] = bad["pairing"][:5]
    expect_reject("B4 five pairs instead of six", "verify_rotation.py",
                  [tmp_json(bad)], ("coordinate pairs",))

    # a2 preserves sign parity and the cell census but breaks the geometry;
    # the image packing gate (<h,h'>^2 <= N*N') is the check that must fire.
    expect_reject("B5 geometry-breaking witness -> image packing gate",
                  "verify_rotation.py",
                  [str(CERTS / "rotation.json"), mutated_witness(a2)],
                  ("violate the 60-degree bound",
                   "violate", "VIOLATION"))

    print("verify_novelty.py battery")
    code, out = run("verify_novelty.py")
    check("C1 pristine novelty accepted",
          code == 0 and "PROVEN" in out and "NON-ISOMETRIC" in out, out)

    # C2: the modular histogram must NOT be proven non-realizable.
    constants = verify_novelty.derive_constants(12)
    transcript: list[str] = []
    proven = verify_novelty.counting_proof(
        32256, 20496, NOVELTY["family"], constants, log=transcript.append)
    check("C2 SOUNDNESS: modular histogram not 'proven' novel",
          proven is False and any("NOT ESTABLISHED" in t for t in transcript),
          "\n".join(transcript))

    bad = json.loads(json.dumps(NOVELTY))
    bad["family"]["n_blocks"] = 50
    expect_reject("C3 wrong family (n_blocks = 50)", "verify_novelty.py",
                  [tmp_json(bad)], ("unsupported family parameters",))

    bad = json.loads(json.dumps(NOVELTY))
    del bad["family"]["n_blocks"]
    expect_reject("C4 missing family parameter", "verify_novelty.py",
                  [tmp_json(bad)], ("four family parameters",))

    expect_reject("C5 corrupted witness propagates", "verify_novelty.py",
                  [str(CERTS / "novelty.json"), mutated_witness(a3)],
                  ("even sign class", "FAIL:"))

    # C6: a parity-preserving but geometry-breaking witness (|cos| = 3/4
    # pairs) passes the structural loader, so the packing precondition
    # inside verify_novelty itself must reject it before any verdict.
    expect_reject("C6 geometry-breaking witness -> packing precondition",
                  "verify_novelty.py",
                  [str(CERTS / "novelty.json"), mutated_witness(a2)],
                  ("not a kissing configuration",))

    print("verify_u11.py battery (v1.1)")
    code, out = run("verify_u11.py")
    check("D1 pristine u11 witnesses accepted",
          code == 0 and "VERIFIED" in out and "NEGATIVE" in out, out)

    def d2(w):
        w["lines"][0]["signs"][1] *= -1
    expect_reject("D2 single sign flip -> parity gate", "verify_u11.py",
                  [mutated_u11(d2)], ("even sign class",))

    li, p, q = u11_conflicting_double_flip()

    def d3(w):
        w["lines"][li]["signs"][p] *= -1
        w["lines"][li]["signs"][q] *= -1
    expect_reject("D3 double sign flip -> packing violation", "verify_u11.py",
                  [mutated_u11(d3)], ("VIOLATION", "violate"))

    def d4(w):
        w["lines"][-1] = json.loads(json.dumps(w["lines"][0]))
    expect_reject("D4 duplicated line", "verify_u11.py",
                  [mutated_u11(d4)], ("duplicate line",))

    def d5(w):
        w["lines"].pop()
    expect_reject("D5 deleted line -> the >= 98 claim gate", "verify_u11.py",
                  [mutated_u11(d5)], ("needs a 98-line witness",))

    print(f"selftest: {passed} passed, {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
