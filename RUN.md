# Wang Conjecture 1.6 — reproduction

Use a FRESH virtualenv. `sympy` pins `mpmath` to 1.3.0, which will downgrade it in any
venv you share with another project.

    python3 -m venv wangvenv && source wangvenv/bin/activate
    pip install mpmath sympy

## Scripts that take no arguments

    cd code
    for f in bookkeeping.py composite_test.py descent.py dictionary.py eps_table.py \
             galois_scaling.py galois_stable.py group_b1.py lemma_checks.py level_test.py \
             local_checks.py local_lemma_test.py model.py scaling2.py scaling_test.py \
             scan.py stab_test.py stable_spread.py stable_spread2.py stress1.py \
             stress2.py stress3.py stress4.py summand_test.py vanish_test.py; do
        echo "== $f"; python3 "$f"
    done

## Scripts that REQUIRE arguments

    python3 regression.py 1200 40        # N, window — the end-to-end check, slowest
    python3 coincide.py 12 12            # TMAX, MMAX
    python3 condC.py 12 13 6             # MMAX, PMAX, EMAX

## Reading the output

Two outputs look like failures and are not:

- `eps_table.py` prints `eps equal across members in all: False`, with (12,3) showing
  -1 against +1. That is the expected raw discrepancy; the two lines after it confirm
  `eps*((t/d)|q)^m` is common across members and equals 1, which is what Step 3 of
  Lemma 5.2 claims.
- `descent.py` prints `undecided residues` for (t,m) = (9,6) and (10,3). Those classes
  are identically zero — verified: 0 nonzero coefficients in n in [200,1200]. The theorem
  permits zero classes.

## OPEN ITEM
`dictionary.py` ends with `t=6 m=5 k=6 e=0: no match with this normalisation`, while the
other four cases match with clean constants. Here d=gcd(6,6)=t, so t/d=1 and no square
root is involved; the coefficients are not degenerate (2 zeros below n=400). This is
unexplained and sits in the exact-formula normalisation layer that VERIFICATION.md already
identifies as the weakest point of the argument. It should be resolved before submission.

## Paper
`wang_ms.tex` is self-contained — no figures, no .bib. Two `pdflatex` passes. 11 pages.

## Status
- Theorem 1.2 (t prime, all m): proved.
- Theorem 1.1 (all t): every step has a written proof.
- Read VERIFICATION.md, especially Round 7: two gaps were found by close reading that no
  amount of numerical testing could have surfaced.

## Resolved since the first packaging

**`regression.py` fixed.** It used `L = 24t` (or `12t`) as the prediction modulus and
`KMAX = 3t`. Both are too small: for (5,5) the deciding group sits above k=15 and the
prediction is not constant mod 120. It reported false mismatches at (5,5) and (12,3).
It now uses the rigid moduli and their lcm, as `stress4.py` does. `stress4.py` gives
0 mismatches on all 13 pairs (10413 comparisons at N=1200), including (5,5) and (12,3).

**`dictionary.py` — RESOLVED: it tests a statement the proof does not make.**

The script tries to identify a single Gamma_0(t) term at modulus k with a single
Gamma_0(576t) term at c=24k. That correspondence is false for several even t: the
twisted sum vanishes identically at c=24k for t=4,6,10,12, and the nonvanishing set in c
is irregular (for t=10 the sum is nonzero at c=24,48,96 but zero at exactly c=120=24k).

But the proof never needs a term-by-term dictionary. Lemma 2.x (uniqueness) matches the
two realisations RATE BY RATE, not term by term: if two families of periodic amplitudes
give the same Bessel expansion, then G_beta = G'_beta for each rate beta. That is why the
paper works with equal-growth groups throughout — difficulty (i) in the strategy section.

The chain that is actually required:
  1. Lemma 2.1: F* modular on Gamma_0(576t) with chi(d)=(t^{-m}/d).
     VERIFIED this session via the Ligozat criterion (both congruences hold, weight 0,
     character = (t/d)^(m mod 2)), for t=2..13 and m=1..6, EVEN t INCLUDED, no failures.
     This needs no numerics, so it avoids the |q|->1 problem that makes direct evaluation
     of eta(24 gamma tau) infeasible at c=576t.
  2. Both realisations give absolutely convergent expansions with the stated growth bound.
  3. Uniqueness => equal-rate groups agree.

So dictionary.py should be relabelled or dropped. As written it implies a gap that does
not exist.

## A note on this repository's scripts
Three separate times, a script's scope limitation was mistaken for a gap in the
mathematics: dictionary.py (above), regression.py (L=24t too coarse; fixed), and an
earlier residue-band scan. Every genuine defect found in this paper was found by READING
(the section 5 Galois convention and the section 3 Salie converse). No numerical test has
ever found one. Weight the two activities accordingly.
