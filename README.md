# Sign periodicity for integer powers of the infinite Borwein product

This repository contains a manuscript proving **Wang's Conjecture 1.6** (Adv. Appl. Math. 141, 2022): for all integers t ≥ 2 and m ≥ 1, the signs of the coefficients of (q;q)^m_∞/(q^t;q^t)^m_∞ are eventually periodic, with least period divisible by t. It also contains the code used to verify every step.

**Status: the proof has not been independently reviewed by anyone with the relevant expertise.**
- The prime-level case (Theorem 1.2) is the most robust part.
- The composite case (Theorem 1.1) relies on §5 — the Galois scaling law and the stable-range proposition. Reviewers should focus there.
- §2's setup, Corollary "decide" and Lemma 5.1 have never been read closely.

`VERIFICATION.md` records nine rounds of self-review: every flaw found and every check run.

**What the verification record shows, and it is worth stating plainly:** both defects found in this paper (the §5 Galois convention and the §3 Salié converse) were found by *reading*. Every quantity involved in each was correct, so roughly 150,000 numerical checks and 15 machine-verified goals passed straight over both. Numerical and SMT verification here is evidence about the stated quantities, not about whether each step follows from its predecessors — and the latter is where the errors actually were.

## Layout
- `wang_ms.tex`, `wang_ms.pdf`: the manuscript (11 pages).
- `VERIFICATION.md`: the end-to-end verification record, rounds 1–9.
- `RUN.md`: how to run everything, and how to read two outputs that look like failures but are not.
- `code/`: verification scripts (Python 3 with `mpmath` and `sympy`).
  - `model.py`: explicit Rademacher-type model of the coefficients.
  - `lemma_checks.py`, `bookkeeping.py`: eta-multiplier and omega formulas.
  - `level_test.py`: the prime-level Level Lemma.
  - `galois_scaling.py`, `summand_test.py`, `eps_table.py`: the Galois scaling law (Lemma 5.2, Prop. 5.3).
  - `galois_stable.py`, `local_lemma_test.py`, `stress2.py`: the stable range (Prop. 5.4).
  - `regression.py`, `stress1.py`, `stress3.py`, `stress4.py`: end-to-end stress tests against brute force.
  - `code/machine/`: machine verification — Z3 (7 goals) and exhaustive residue enumeration (8 goals, complete proofs rather than bounded checks). See its README for what is and is not covered.
- `archive/`: earlier drafts and development notes.

## Reproduce

    python3 -m venv wangvenv && source wangvenv/bin/activate
    pip install mpmath sympy z3-solver
    cd code
    python3 stress4.py 1800        # rigid-only sign prediction vs brute force
    python3 local_lemma_test.py    # local stationary-phase lemma
    python3 stress1.py             # exact summand identity, random terms
    cd machine
    python3 mv1_z3.py              # SMT goals
    python3 mv2_residues.py        # exhaustive residue proofs

Several scripts take command-line arguments; `RUN.md` lists them. Use a fresh virtualenv: `sympy` pins `mpmath` below 1.4 and will downgrade it in a shared environment.

## License
MIT (see `LICENSE`).
