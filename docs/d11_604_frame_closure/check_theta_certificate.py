#!/usr/bin/env python3
"""Standalone checker for theta_dual_certificate.json (Lovasz theta dual).

Adapted (paths only) from the campaign checker
`arena/workspace/kissing_theta_cert_check.py`; the certificate producer is
`arena/workspace/kissing_theta_dual_certify.py`.  Independent of the producer:
it rebuilds M = t*I - J - Y exactly from the certificate's scaled integers and
the METIS edge list, re-verifies the pinned graph hash, and re-runs the
verified-Cholesky PSD criterion (dpotrf success on M - shift*I with the shift
exceeding 100x a rigorous floating-point error bound).  Prints the certified
conclusion or FAILS loudly.

M PSD certifies theta(G) <= t_cert (an exact dyadic rational), hence
alpha(G) <= floor(t_cert) = 251.  Together with the shipped 240-line witness
(see check_graph_sanity.py) this gives the rigorous bracket
240 <= alpha(G) <= 251, independent of the SCIP branch-and-cut proof.

Requires numpy (the only third-party dependency of this addendum's checks).

Usage: python3 check_theta_certificate.py [certificate.json graph.metis]
(defaults: the copies shipped next to this script).
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def main():
    cert_path = sys.argv[1] if len(sys.argv) > 1 else HERE / "theta_dual_certificate.json"
    metis_path = sys.argv[2] if len(sys.argv) > 2 else HERE / "restricted_graph.metis"

    cert = json.load(open(cert_path))
    lines = open(metis_path).read().strip().split("\n")
    n, m = map(int, lines[0].split())
    edges = set()
    for i, ln in enumerate(lines[1 : 1 + n]):
        for tok in ln.split():
            j = int(tok) - 1
            if j != i:
                edges.add((min(i, j), max(i, j)))
    edges = sorted(edges)
    assert len(edges) == m == cert["graph"]["m"] and n == cert["graph"]["n"]
    h = hashlib.sha256(json.dumps(edges).encode()).hexdigest()
    assert h == cert["graph"]["edges_sha256"], "edge-set hash mismatch"

    scale = 1 << cert["scale_bits"]
    y_int = cert["y_scaled_ints"]
    assert len(y_int) == m
    t_cert = cert["t_cert_scaled_int"] / scale

    M = np.full((n, n), -1.0)
    np.fill_diagonal(M, t_cert - 1.0)
    for (a, b), yi in zip(edges, y_int):
        M[a, b] = M[b, a] = -1.0 - yi / scale

    u = 2.0 ** -53
    gamma = (n + 1) * u / (1.0 - (n + 1) * u)
    bound = gamma * float(np.max(np.sum(np.abs(M), axis=1)))
    shift = 2.0 ** -16
    assert shift > 100.0 * bound, (shift, bound)
    np.linalg.cholesky(M - shift * np.eye(n))  # raises LinAlgError on failure

    alpha_upper = int(np.floor(t_cert))
    print(f"CHECK PASSED: M = t*I - J - Y is PSD (verified, err bound {bound:.2e})")
    print(f"theta(G) <= t = {t_cert!r}  =>  alpha(G) <= {alpha_upper}")


if __name__ == "__main__":
    main()
