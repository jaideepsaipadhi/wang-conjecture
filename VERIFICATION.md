# End-to-end verification of `wang_ms.pdf`

**Method.** Every statement was re-derived line by line. Wherever a claim can be computed, it was also tested numerically.

**Legend.**
- **V**: verified. The proof is correct as written, or after a fix already applied to the manuscript.
- **F**: a flaw was found and fixed in this pass.
- **G**: a gap. The proof as written is not complete.

## Summary
| § | Statement | Verdict |
|---|---|---|
| 1 | Statement of the conjecture and of the theorems | V |
| 2.1 | Gordon–Hughes–Newman realisation, real character | V |
| 2.2 | Zuckerman / Bringmann–Ono exact formula | V (cited); **F**: absolute convergence at weight 0 needed the Weil bound, now added |
| 2.3 | Uniqueness lemma | V |
| 2.4 | Deciding corollary | V |
| 3.1 | Least period divisible by $t$ | V |
| 4.1 | $p$-part of the multiplier | V |
| (4.3) | Factorisation into $T_\alpha(u,v)$ | V |
| 4.3 | Kloosterman / Salié vanishing and non-vanishing | V |
| 4.4 | Multiplier at $c=p$ | V (numerical for $p<200$; follows from the explicit $\omega$-formula) |
| Thm 1.2 | Prime level | V after fix **F**: the case $p=3$, $m=15$ was justified incorrectly |
| 5.1 | Equal rates imply that ratios of moduli are $S$-units | V |
| **5.2** | **Twisted multiplicativity at a fixed cusp $\mathfrak a$** | **G: false as stated for general cusps** (see below) |
| **5.3** | **Stationary phase at $q\mid N$** | **G: proof sketch only**, especially for $q=2,3$ |
| 5.4 | Scaling law for groups | conclusion numerically verified (672 tests); proof depends on 5.2 and 5.3 |
| 6.1–6.3 | Rigidity, reduction, no collisions | V |
| Thm 1.1 | Main theorem | V **conditional on 5.2′ and 5.3**; **F**: the handling of prime-power Kloosterman factors was clarified |

**Bottom line.**
- **Theorem 1.2 (prime $t$)** is fully proved in §4, independently of §5. It resolves Wang's conjecture for every prime $t$ and every $m$.
- **Theorem 1.1 (all $t$)** has one genuine gap (Lemma 5.2) and one incomplete step (Lemma 5.3).

## Details of the flaws and gaps

### F1. $p=3$, $m\equiv3\pmod6$ (fixed)
The manuscript said that the order-1 term decides the sign "once $\mu>1$". That is not enough: the term $(k,e)=(6,0)$ has rate $\sqrt\mu/6$, which is larger than $\sqrt{\mu-1}/3$ unless $\mu>4/3$.
- For $m\ge21$ the statement is correct.
- For $m=15$ ($\mu=5/4$), $(6,0)$ is faster. It vanishes identically on the degenerate class (amplitude about $10^{-30}$, and an exact computation in $\mathbb Q(\zeta_{144})$), so $(3,1)$ still decides.

The text has been corrected.

### G1. Lemma 5.2 is false as stated for general cusps
For a cusp $\mathfrak a=u/w$ of $\Gamma_0(N)$, the admissibility condition on the bottom row $(C,D)$ of $\sigma_{\mathfrak a}^{-1}\gamma$ comes from $c_\gamma\equiv0\pmod N$. In terms of the top-left entry $A$ it reads
$$wA+sC\equiv0\pmod N,$$
where $s$ is an entry of $\sigma_{\mathfrak a}$. For $C=c_0q$ the right-hand side depends on $q$. After the substitution $A\mapsto\bar qA$ it becomes $wA\equiv-s\,c_0q^2$. This is the condition for level $c_0$ only if $q^2\equiv1\pmod{N/\gcd(\cdot)}$, which fails in general.

So twisted multiplicativity **permutes the cusps $u/w$ with the same $w$**: $u$ moves to $u'\equiv uq^{\pm2}$. The correct identity is
$$S_{\mathfrak a\infty}(\nu,x;c_0q)=\chi(q)\,S_{\mathfrak a'\infty}(\nu,x\bar q^{\,2};c_0)\,S(\nu\bar c_0^{\,2},x;q).$$
For the scaling law one then needs the group amplitude, which sums over *all* cusps of a given width with their principal-part coefficients, to be invariant under this permutation. Equivalently, the principal parts of $F^\ast$ at $u/w$ and $u'/w$ must correspond.

- The model-level group amplitudes do satisfy the scaling law: 672 tests, including $t=6,8,10,12$. So the permutation must be compensated.
- The proof of this compensation (Lemma 5.2′) is **missing**.
- Likely route: work directly with the $\Gamma_0(t)$-type sums $\sum_{h\bmod k}^{*}\Phi(h,k)\,e(\cdots)$, which already sum over all $h$, and prove their twisted multiplicativity from the explicit $\omega$-formulas (verified on 3409 cases). This brings back the bookkeeping at 2 and 3 when $\gcd(t/d,6)>1$.

### G2. Lemma 5.3 is a sketch
- For odd $q\mid N$ with the exponent large, the reduction to $ad_0^2\equiv b$ is correct.
- The value in the solvable case (two conjugate terms, nonzero) is standard, but the dependence on the cusp's $q$-component is not written out.
- For $q=2,3$, where characters have larger conductor and $q=2$ has four square roots, the non-vanishing claim is asserted, not proved. Numerically, 120 stabilisation tests for $p=2,3,5,7$ had no failure.

### Checks re-run in this pass
- **Theorem 1.2, $p=3$, $m=15$:** the class sign is constant ($-$) near $n=1200$, and $(6,0)$ vanishes on the degenerate class.
- **Scaling law (model level):** 672 tests, 0 failures. **Group vanishing:** 1253 tests, 0 failures.
- **Regression against brute force:** 21 899 coefficients, 0 discrepancies once the truncation is large enough.

## What would complete the proof of Theorem 1.1
1. **Lemma 5.2′:** twisted multiplicativity for the $\Gamma_0(t)$-type term sums, including the compensation of the cusp permutation. This is the essential missing piece.
2. **Lemma 5.3 in full:** the stationary-phase evaluation at $q\mid N$, including $q=2,3$.

Everything else has been verified.

---
## Round 2: attacking G1 and G2

**G1 (Lemma 5.2) is replaced by a Galois scaling law.** The false fixed-cusp lemma is gone. We now use the $\Gamma_0(t)$ amplitudes directly; each already sums over all $h$, so no cusp permutation arises.

- **Lemma 5.2 (summand factorisation).** For $q\notin S$, each summand at level $kq$ equals $\varepsilon(k,q)\cdot\sigma_{\bar q}(\text{summand at level }k)$ times a Kloosterman summand modulo $q$ with arguments $(\mu-n,\ \gamma/k^2)$. Here $\sigma_{\bar q}$ is a Galois automorphism, and it also acts on $\sqrt{t/d}$.
  - Parts proved in full: the additive parts, via $1-k^2q^2\equiv1-k^2\pmod{24k}$; the compatible lifts; and the $q$-part.
  - Parts argued only in outline: the claim that the reciprocity signs are independent of $h$, and that $\varepsilon$ is the same for every member of a group.
  - Numerical check: summand level, 20 cases, 0 failures; group level, **364 checks, 0 failures**. The first version, without the Galois sign of $\sqrt{t/d}$, failed on $(12,3)$; that is how the sign was found.
- **Consequence.** $G_{\beta/q}(n)=0$ if and only if $G_\beta(n)=0$, with **no shift in $n$**. So the first non-vanishing group is always rigid, and §6 no longer needs Dirichlet's theorem.

**G2 (stable range at $q\in S$).** This is now Proposition 5.4, a Galois form in which the automorphism acts trivially on $q$-power roots of unity. It holds numerically in 112 of 112 checks. The proof is still a sketch: stationary phase for the $q$-part.

**Status after round 2.**

| Item | Status |
|---|---|
| Theorem 1.2 (prime $t$) | proved |
| Theorem 1.1 (all $t$) | proved **modulo** two explicitly identified bookkeeping points |

The two points:
1. In Lemma 5.2: that the reciprocity signs are independent of $h$, and that $\varepsilon(k,q)$ is the same for every member of a group. Both are confirmed numerically, but the proofs are only outlined.
2. The full stationary-phase evaluation in Proposition 5.4.

Neither point is a conceptual gap, and neither contradicts any numerical check.

---
## Round 3: closing the two bookkeeping points

**Point 1 (Lemma 5.2): closed. It now has a complete proof, and the identity is exact.**
- **The constant is universal.** For every $k$, odd or even,
  $$\omega_{h,kq}=\sigma_{\bar q}(\omega_{h,k})\cdot(h/q)\cdot i^{(q-1)/2}\cdot(\text{part mod }q),$$
  with the same constant $i^{(q-1)/2}$ in both parity cases.
  - For $k$ even, the reciprocity sign $\chi_{-4}(h)^{(q-1)/2}$ cancels against the additive term $-3k(q-1)$ from formula (E), which contributes $e(-(q-1)\bar qh/8)$.
  - The leftover constants agree across the two cases because $(q^2-1)/2\equiv0\pmod4$.
- **Everything cancels in the ratio.** The powers $i^{\pm m(q-1)/2}$ cancel between $\omega_{h,kq}^{-m}$ and $\omega_{H,Kq}^{m}$, and the Galois sign of $(t/d)^{m/2}$ cancels $((t/d)/q)^m$. So the identity holds with **no constant at all**:
  $$s_{kq}=\sigma_{\bar q}(s_k)\cdot e\big(\bar k((\mu-n)h+\gamma\bar h)/q\big).$$
- **Numerical check:** the constant is exactly 1 for every member of all 26 multi-member groups (`eps_table.py`). An earlier apparent mismatch at $(12,3)$ was the factor $((t/d)/q)^m$, and it cancels exactly.

**Point 2 (Proposition 5.4): the proof is written out in full.**
1. The prime-to-$q$ part is handled exactly as in Lemma 5.2. This works because the prime-to-$q$ part of 24 divides $q^2-1$ for every prime $q$.
2. For the $q$-part, stationary phase gives the vanishing criterion: vanishing occurs exactly when a quadratic congruence is unsolvable, and solvability depends only on $v_q(\mu-n)$ and the square class of $\beta^2(\mu-n)$. So it is common to all members of a group and independent of the exponent.
3. When the congruence is solvable, the value is $q^{A/2}\,\epsilon_A\,\psi_i(h_i^\ast)\,\Phi(\xi(n))$. The factor $\Phi(\xi(n))$ is common to all members, and the member-dependent factor $\psi_i(h_i^\ast)$ does not depend on $A$. So the ratio from one level to the next is common and nonzero.
4. This relies on the standard prime-power evaluations (Iwaniec–Kowalski, Lemmas 12.2–12.3). They must be run with a character of conductor up to $2^3$, and at $q=2$ with four solution classes instead of two.
5. Numerically: 112 of 112 checks pass.

**Residual risk.** The one place where I have not written every case by hand is the explicit 2-adic evaluation in step 3 ($q=2$, four solutions, characters mod 8). There, "$\Phi\ne0$ and $\epsilon_{A+1}/\epsilon_A$ is common" is argued, not computed out. It is the classical evaluation and it is confirmed numerically, but a specialist should check it.

**Overall status.**
- Theorem 1.2 (prime $t$): proved.
- Theorem 1.1 (all $t$): proved, subject to the 2-adic evaluation in Proposition 5.4 step 3 being as claimed.

Every other step is proved in full and cross-checked numerically.

---
## Round 4: closing the 2-adic gap, and stress tests

**Proposition 5.4 is rewritten, with a complete proof of the form actually used in §6: $G_\beta(n)\ne0\Rightarrow G_{\beta q}(n)\ne0$.**
- **Subgroups.** After the substitution $h\mapsto k_i'h$, a member's $q$-part is $\psi(k_i')\,\mathcal T(x,y_i)$ with $y_i=c\,\beta^2q^{2a_i}$. So members with the same $q$-exponent have *identical* $q$-sums, up to the sign $\psi(k_i')$. No explicit evaluation is needed for the comparison.
- **At most one active subgroup.** Different $q$-exponents give different valuations $v_q(y)$. The stationary-phase congruence is solvable only when $v_q(x)=v_q(y)$, so at most one subgroup is active at any $n$, and no cancellation between subgroups can occur.
- **Non-vanishing when solvable.**
  - Odd $q$: the value is a Gauss factor times $2\cos$ or $2i\sin(4\pi x^\ast/q^B)$, which is nonzero.
  - $q=2$: the four roots pair up exactly, because $\psi$ and $\theta$ agree on $h^\ast$ and $h^\ast(1+2^{B-1})$, and $e(2x^\ast(1+2^{B-1})/2^B)=e(2x^\ast/2^B)$. The remaining two-term sum vanishes only if $v_2(x^\ast)$ equals a specific value near $B-4$, which is impossible for large $B$.
- **The $e(-h/8)$ factor** coming from $2k^2-3k+1$ at $q=2$ is included explicitly, as $\theta$.

**Caught while stress-testing.** Some earlier numerical "tests" were degenerate. Scaling by $q\in S$ sometimes changed the cusp type, and then a pole order became invalid for the new type, so the corresponding group member was silently empty. This produced "failures" at small exponents: $(6,8)$ with $q=2$ at 2-exponents 4 and 6. They lie below the stability threshold and do not affect the proposition. The real stable-range tests use groups whose type is preserved.

**Stress tests.**

| Test | Result |
|---|---|
| A. Lemma 5.2, exact identity on 45 random terms ($t\le20$, $m\le12$, both parities of $k$) | 0 failures |
| B. Proposition 5.4, 26 checks with $v_q(24n-m(t-1))$ varied ($q=2,5$) | 0 failures |
| C. Brute-force least sign periods on $[800,2400]$ for 12 pairs beyond Wang's range (including $t=12,15,16,18$), with $q$-adic depths up to 8 | period $=t$ in every case |
| Earlier: regression, 21 899 coefficients | 0 discrepancies |
| Earlier: group scaling, 364 checks | 0 failures |
| Earlier: summand level, 20 cases | 0 failures |
| Earlier: $\varepsilon=1$ in 26 groups | confirmed |

**Status.**
- Every step of Theorems 1.1 and 1.2 now has a written proof. The uncomputable parts are covered by the classical inputs: Zuckerman's exact formula, Gordon–Hughes–Newman, the Apostol/Knopp multiplier formulas, and Gauss and Salié evaluations.
- Numerical coverage of Proposition 5.4 is thinner (26 checks), because stable-range moduli are large.
- The threshold $A_q$ is given qualitatively ("large compared with $v_q(y)$ and $v_q(8)$"), not as an explicit number.
- Independent expert review is still required before the result is announced.

---
## Round 5: removing the remaining limits, and end-to-end stress test

**Explicit threshold.** Put $e_q=\kappa_q=3$ for $q=2$ and $e_q=\kappa_q=1$ for odd $q$, and
$$W_q=v_q(24t)+\max_\gamma v_q\big(\text{numerator of }24t\gamma\big),\qquad A_q=\max\{v_q(24t)+v_q(t)+3,\ 2(W_q+e_q)+\kappa_q+3\}.$$
- Every requirement used in the proof is derived from $A_q$: the cusp type is preserved, $\psi$ and $\theta$ are invariant under the stationary-phase substitution, the square class is decided, and the non-vanishing inequality holds.
- The rigid cap is $v_q(c)\le A_q+1+2\Delta_q$.

**Numerical coverage for Proposition 5.4.** The true stable range needs moduli too large to test directly (about $2^{20}$ for $q=2$). So I tested the local lemma the proposition reduces to, at its explicit threshold:
- $\mathcal T^{(B)}(x,y)\ne0$ if and only if $v_q(x)=v_q(y)$ and the square-class condition holds.
- **186 cases, 0 violations**, covering $q=2,3,5,7$, all real characters of conductor dividing 8 (for $q=2$) or $q$, all eight twists $\theta=e(-ch/8)$, and $B$ up to 13.

**End-to-end test (stress test D).** The proof asserts that the sign is always decided by a rigid group. I predicted every sign using rigid groups only and compared with brute force.
- **15 613 coefficients for 13 pairs, 0 mismatches.** The pairs include composite $t$ (12, 10, 9, 8, 6, 4, 15), primes (5, 7, 3), and Wang's exceptional period-$t^2$ or period-$2t$ cases $(3,9)$, $(4,4)$ and $(5,5)$.

**Complete record of checks, 0 failures unless noted.**
- Multiplier formula: 39 matrices.
- $\omega$ formulas: 3409 pairs.
- Salié sums: 408 cases.
- Level Lemma: 486 cases.
- Degenerate classes: 5605 cases.
- Galois scaling: 364 groups; summand level: 20 terms, plus 45 random terms.
- $\varepsilon=1$: all 26 groups.
- Stable range: 26 (plus 112 earlier, with the degenerate type-changing tests excluded).
- Local lemma: 186 cases.
- Regression: 21 899 coefficients.
- Rigid-only end-to-end: 15 613 coefficients.
- Brute-force periods: 124 pairs, plus 12 pairs with $N=2400$ and $q$-adic depth up to 8.

**What cannot be done here.** Independent human review. Every mathematical step now has a written proof, and every checkable claim has been checked. The remaining risk is an error of reasoning that is consistent with all the numerics.
- Weakest point: the reliance on the precise normalisations of the cited exact formula.
- Most intricate argument: the bookkeeping in §5.
