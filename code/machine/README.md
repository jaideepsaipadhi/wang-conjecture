# Machine verification

Two scripts. Run both; each prints PROVED or a counterexample per goal.

    python3 mv1_z3.py         # Z3 SMT
    python3 mv2_residues.py   # exhaustive residue enumeration

## What is actually PROVED (not sampled)

`mv1_z3.py` — 7 goals discharged by Z3 as UNSAT of the negation, over bounded integers:
  - Step 1: (kq-1) - q(k-1) = q-1
  - Step 2: H = (t/d)h mod q
  - Step 2: q l' = l mod 24  =>  Kql' = Kl mod 24K
  - Step 3: h-coefficient collapses to mu
  - Step 3: hbar-coefficient collapses to gamma + Lambda
  - Ligozat (i) and (ii)

Three further goals returned UNKNOWN (nonlinear modular arithmetic is hard for SMT).
All three depend ONLY on residues, so `mv2_residues.py` closes them COMPLETELY by
enumerating every residue class — that is a full proof of the universally quantified
statement, not a bounded check:
  - gcd(q,6)=1 => q^2 = 1 mod 24   (8 classes; gives 24 | k(q^2-1) for EVERY k)
  - q*qbar = 1 mod 8 => qbar = q mod 8   (4 classes)
  - q odd => 8 | (q-1)(q+1)   (4 classes)

`mv2_residues.py` also discharges, exhaustively over the stated range:
  - Prop 3.x: t/g never divides a unit mod t when g<t   (t<=200, all divisors, all units)
  - the character identity (t^{-m}/d) = (t/d)^(m mod 2)   (full period in d)
  - Salie (c): sum of the Legendre symbol over units = 0   (all odd primes p<200)
  - Salie converse: cos(2 pi x/p) != 0 for integer x, odd p   (odd p<300, all x)

Total: 15 goals, 15 discharged.

## What is NOT machine-verified, and cannot be with these tools

Everything analytic. Specifically:
  - the Rademacher-Zuckerman expansion itself and its convergence;
  - the uniqueness principle (Lemma 2.x) — an asymptotic argument;
  - Prop 5.4's stationary-phase evaluation (the local lemma is checked numerically at
    thousands of parameter values, which is evidence, not proof);
  - Prop 5.3's Kloosterman non-vanishing (the congruence argument mod (1-zeta_q) is a
    proof; the 20477-case check only confirms it);
  - the finiteness/termination argument of section 6;
  - and, most importantly, whether each step FOLLOWS from its predecessors. Every
    individual claim can be true while a deduction between them fails. No tool here
    detects that.

Both defects found in this paper (the section 5 Galois convention and the section 3
Salie converse) were of exactly that kind: every quantity involved was correct, and they
were found by reading. Roughly 150,000 numerical checks never surfaced either one.

A full formalisation would need Lean/mathlib with Dedekind eta multipliers, Kloosterman
sums and the circle method. mathlib does not currently have that machinery.
