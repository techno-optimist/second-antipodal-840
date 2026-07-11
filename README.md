# A new 420-line antipodal kissing configuration in R^12 — machine verification

Machine-verification layer for paper #6 of the ProjectForty2 exact-certificate
note series. Every claim below is recomputed from minimal certificates by
stdlib-only exact checkers — the certificates store **no counts, no scores,
no verdicts** (recompute, never echo).

## Certified statements

1. **Witness** (`make verify-witness`): 420 distinct lines in R^12 — 132 W2
   lines (e_a ± e_b)/√2 and 288 U lines (±1/(2√2) on 8-coordinate supports,
   even sign class, 15 partial cells) — with all C(420,2) = 87,990 pairs at
   |cos| ≤ 1/2, checked as the exact integer inequality 4⟨g,g′⟩² ≤ N·N′.
   Doubling to ±u gives an 840-point antipodal kissing configuration:
   **N60(12) ≥ 420 lines / 840 points.** Recomputed along the way: the
   composition, the 20,464 zero-margin (touching) pairs, the exact |cos|
   histogram {0: 35142, 1/4: 32384, 1/2: 20464}, and the U-cell structure.

2. **Rotation** (`make verify-rotation`): the blockwise 45° rotation on the
   stored coordinate pairing (exact over Q(√2), realized as the integer map
   h_i = g_i − g_j, h_j = g_i + g_j) carries the witness **isometrically**
   (all 420 norms + 87,990 Grams preserved exactly; the exact 60° bound is
   re-checked on every image pair as a falsifiable packing gate) into the
   integer A ∪ W4 world: 12 axes + 408 W4 lines (47 whole cells + 8 half
   cells on 55 supports, recomputed). Hence the witness is Hanani-capped in
   its slice.

3. **Novelty** (`make verify-novelty`): a packing-independent counting
   theorem — **no** configuration of the form 12 axes + 51 whole W4 cells
   (the family containing the known modular 420) realizes the witness
   histogram. Every interaction constant (64/P1-pair, 2448 + 32/P2-pair,
   share-3 illegality, λ_p ≤ 5, the DP minimum 94) is derived by explicit
   finite computation inside the checker, not assumed. The chain forces
   P1 = 506, P2 = 563, r_i ≡ 17, then proves P2 ≥ 564: contradiction.
   Since the pairwise |cos| multiset is an isometry invariant, the witness
   is **proven non-isometric to the known modular 420**.

4. **Odd-overlap structure facts** (`make verify-structure`): every MACHINE
   entry of the note's Section "The odd-overlap frontier" is recomputed
   exactly from the conflict predicate (stdlib big-integer branch-and-bound,
   ~5 s): the structure-lemma censuses 256/128/512 at defect overlaps
   t = 1/2/3, the pair cap κ = 64 for all t, the t=3 rigidity 57 (with
   folded-7-cube connectivity 7 by exact max-flow) versus 64 at t = 1/2,
   the W7 localization 16 × 4 = 64, the corner collapse 64, the e7/e8
   even-family caps 64/128, the two sunflower-triple caps of 80, the
   uniform star saturation 96 × 12 of the shipped witness (tightness of the
   star reduction α(U12) ≤ 3·α(U11)), and U11-world sanity (165 cells,
   10,560 lines, no free cell pairs).

Honest scope: these checkers certify the witness, its slice-cap, its
novelty, and the local structure facts. They do NOT certify α(U-world) = 288
or the global optimality ladder (CP-SAT levels 0–3 exact, level ≥ 4
searched); see the note's text for the exact ladder statement, the LP
barrier theorem, and the star reduction to the open problem α(U11) ≤ 96.
Solver receipts (CP-SAT ladder, θ/ILP sweep verdicts) ship under
`certs/receipts/` and are labeled as solver results, not stdlib-recomputable
certificates.

## Usage

```
make verify            # the three core exact checkers (stdlib only, ~30 s)
make verify-structure  # the odd-overlap-section MACHINE facts (~5 s)
make selftest          # adversarial battery: 21 corruption/soundness cases
make pdf               # typeset the note (needs tectonic)
```

## Layout

```
certs/config420.json   the witness: 420 line descriptions, nothing else
certs/rotation.json    the coordinate pairing, nothing else
certs/novelty.json     the four family parameters, nothing else
certs/receipts/        solver receipts (NOT stdlib-recomputable):
                         result_slice_only.json, level2.json, level3.json
                           — CP-SAT OPTIMAL receipts, ladder levels 0/2/3
                             (level 1: V2_DECISION_SUMMARY.json →
                             exact_results.LEVEL1_complete — all 480
                             foreign cells at 288, all OPTIMAL)
                         randk.json, twoslice.json — search beyond level 3
                         witness288.json — the 288 U-line warm start
                         V2_DECISION_SUMMARY.json — the v2 decision record
                         K2_stats.json — mixed-graph edge census (below)
                         SWEEP_VERDICTS.json — Round-88 sweep verdicts
scripts/witness.py     shared exact loader (structural validation only)
scripts/verify_420.py  packing + counts checker
scripts/verify_rotation.py  Q(sqrt2) isometry + A∪W4 membership checker
scripts/verify_novelty.py   counting-theorem checker
scripts/verify_structure.py odd-overlap-section MACHINE-facts checker
scripts/selftest.py    corruption battery (must all be REJECTED)
scripts/oddlemma/      the lemma-hunt exact toolkit (verbatim session
                       scripts): uw.py (predicate + exact MIS B&B),
                       mis2.py, u11_world.py, stage*.py, RESULTS.json
```

## Mixed-graph edge census (recomputed in `certs/receipts/K2_stats.json`)

The full mixed conflict graph on 12 A + 132 W2 + 3,960 W4 + 31,680 U =
35,784 vertices has 40,130,904 edges, distributed as:

| pair | edges | pair | edges |
|---|---|---|---|
| A–A | 0 | W4–W4 | 126,720 |
| A–W2 | 264 | W4–W2 | 23,760 |
| A–W4 | 0 | W4–U | 16,410,240 |
| A–U | 0 | W2–W2 | 0 |
| W2–U | 0 | U–U | 23,569,920 |

(predicate-validation mismatches: 0; same-alphabet entries counted once.)

No third-party dependencies (see `requirements.txt`). No floats on any
decision path; displayed decimals use directed rounding.

## DOI

- This version: [10.5281/zenodo.21306850](https://doi.org/10.5281/zenodo.21306850)
- All versions (concept): [10.5281/zenodo.21306849](https://doi.org/10.5281/zenodo.21306849)
