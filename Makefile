# Makefile for the d12 420-line / 840-point antipodal kissing note (paper #6).
#
# `make verify` is the theorem-grade check: three stdlib-only exact checkers
# rebuild everything from the minimal certificates in certs/ (which store NO
# counts, scores, or verdicts) and RECOMPUTE every claim:
#
#   verify_420.py      the witness: 420 distinct lines in R^12, all 87,990
#                      pairs at |cos| <= 1/2 via the exact integer test
#                      4<g,g'>^2 <= N*N'  =>  N60(12) >= 420 lines / 840 points
#   verify_rotation.py exact Q(sqrt2) blockwise 45-degree rotation carrying
#                      the witness isometrically into the integer A+W4 world
#                      (all Grams preserved; image composition recomputed)
#   verify_novelty.py  the counting theorem: no 12-axes + 51-whole-W4-cell
#                      configuration (including the known modular 420)
#                      realizes the witness |cos| histogram => PROVEN
#                      non-isometric (constants derived, not assumed)
#
# No floats appear on any decision path; displayed decimals (none load-
# bearing) use directed rounding.

PYTHON ?= python3
TEX     = second_antipodal_840.tex

.PHONY: all verify selftest verify-witness verify-rotation verify-novelty \
        verify-structure pdf clean

all: verify

## verify: exact certification of witness + rotation + novelty (stdlib only).
verify: verify-witness verify-rotation verify-novelty

verify-witness:
	$(PYTHON) scripts/verify_420.py certs/config420.json

verify-rotation:
	$(PYTHON) scripts/verify_rotation.py certs/rotation.json certs/config420.json

verify-novelty:
	$(PYTHON) scripts/verify_novelty.py certs/novelty.json certs/config420.json

## verify-structure: recompute every MACHINE fact of the note's odd-overlap
## section (structure-lemma censuses, pair/triple caps, t=3 rigidity,
## W7/e7/e8 localization, sunflower 80s, star saturation of the witness,
## U11 sanity).  Stdlib only, exact, ~5 s.  Additive: not required for the
## three core certificates above.
verify-structure:
	$(PYTHON) scripts/verify_structure.py

## selftest: adversarial check that corrupted certificates are REJECTED
## (and that the novelty chain does NOT "prove" the modular histogram novel).
selftest:
	$(PYTHON) scripts/selftest.py

## pdf: typeset the note.
pdf:
	tectonic $(TEX)

clean:
	rm -f $(TEX:.tex=.aux) $(TEX:.tex=.log) $(TEX:.tex=.out)
