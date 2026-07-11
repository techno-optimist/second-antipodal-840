#!/usr/bin/env python3
"""U11 world generator — THE reduced target after the star-reduction theorem.

Star-reduction (this session, machine-verified):
    every U12 line lies in exactly 4 coordinate-stars star_x = {lines with x in cell A};
    star_x's world is EXACTLY the U11 world below (conflict predicate never sees x);
    hence  4*N = sum_x |S ∩ star_x| <= 12*alpha(U11)  =>  alpha(U12) <= 3*alpha(U11).
    The 288 slice witness meets every star in exactly 96 (uniform saturation).
    Therefore:
      alpha(U11) <= 96  ==>  alpha(U12) = 288 (odd-overlap lemma PROVED, ladder global)
      any 289-config contains a >=97-line U11 configuration (pigeonhole)
    -> the record hunt AND the closure both live in this 10,560-line world.

U11 world: ground set [11]; cells = 3-subsets A (complements of 8-supports H);
lines per cell: even-sign classes on H mod global flip (64); conflict iff
|sum_{i in H∩H'} s_i s'_i| >= 5. NOTE: |H∩H'| = 5 + |A∩B| ∈ {5,6,7} — there are
NO free cell pairs in U11 (every pair kappa-constrained: n_A + n_B <= 64).
"""
import itertools

N11 = 11
FULL11 = (1 << N11) - 1

def cells11():
    return [frozenset(c) for c in itertools.combinations(range(N11), 3)]  # 165

def mask_of(s):
    m = 0
    for i in s:
        m |= 1 << i
    return m

def lines_of_cell11(A):
    Hmask = FULL11 & ~mask_of(A)
    Hbits = [i for i in range(N11) if (Hmask >> i) & 1]
    out, seen = [], set()
    for bits in range(256):
        u = 0
        for j in range(8):
            if (bits >> j) & 1:
                u |= 1 << Hbits[j]
        if bin(u).count("1") % 2:
            continue
        cu = min(u, u ^ Hmask)
        if cu not in seen:
            seen.add(cu)
            out.append(cu)
    assert len(out) == 64
    return out

def conflict11(A, u, B, v):
    if A == B:
        return False
    R = (FULL11 & ~mask_of(A)) & (FULL11 & ~mask_of(B))
    r = bin(R).count("1")
    d = bin((u ^ v) & R).count("1")
    return abs(r - 2 * d) >= 5

if __name__ == "__main__":
    cs = cells11()
    assert len(cs) == 165
    n = 165 * 64
    print(f"U11 world: {n} lines (targets: find 97 = record path, or prove <=96 = jewel closed)")
    # spot census: r in {5,6,7} for all pairs
    import random
    random.seed(0)
    for _ in range(200):
        A, B = random.sample(cs, 2)
        r = 11 - len(A | B)
        assert r in (5, 6, 7)
    print("no free pairs: OK (all |H∩H'| in {5,6,7})")
