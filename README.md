# Sign periodicity for integer powers of the infinite Borwein product

This repository contains a manuscript proving **Wang's Conjecture 1.6** (Adv. Appl. Math. 141, 2022): for all integers t ≥ 2 and m ≥ 1, the signs of the coefficients of (q;q)^m_∞/(q^t;q^t)^m_∞ are eventually periodic, with least period divisible by t. It also contains the code used to verify every step.

**Status: the proof has not yet been independently reviewed.**
- The prime-level case (Theorem 1.2) is the most robust part.
- The composite case (Theorem 1.1) relies on §5, the Galois scaling law and the stable-range proposition. Reviewers should focus there.

`VERIFICATION.md` records five rounds of self-review: every flaw found, and every numerical check.

## Layout
- `paper/wang_ms.tex`, `paper/wang_ms.pdf`: the manuscript.
- `VERIFICATION.md`: the end-to-end verification record.
- `code/`: verification scripts (Python 3 with `mpmath` and `sympy`).
  - `model.py`: explicit Rademacher-type model of the coefficients.
  - `lemma_checks.py`, `bookkeeping.py`: eta-multiplier and ω formulas.
  - `level_test.py`, `prime_zeros.py`: the prime-level Level Lemma.
  - `galois_scaling.py`, `summand_test.py`, `eps_table.py`: the Galois scaling law (Lemma 5.2 and Proposition 5.3).
  - `galois_stable.py`, `local_lemma_test.py`, `stress2.py`: the stable range (Proposition 5.4).
  - `regression.py`, `stress1.py`, `stress3.py`, `stress4.py`: end-to-end stress tests against brute force.
- `archive/`: earlier drafts and development notes.

## Reproduce
```
pip install mpmath sympy
cd code
python3 stress4.py 1800        # rigid-only sign prediction vs brute force
python3 local_lemma_test.py    # local stationary-phase lemma
python3 stress1.py             # exact summand identity, random terms
```

## License
MIT (see `LICENSE`).
