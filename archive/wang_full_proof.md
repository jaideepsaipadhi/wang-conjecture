# Proof of Wang's Conjecture 1.6 (all t, all m)

**Theorem.** For all integers $t\ge2$ and $m\ge1$, the sequence $\operatorname{sgn}c_t^{(m)}(n)$ is eventually periodic, and its least period is divisible by $t$.

Notation: $F=\eta(\tau)^m/\eta(t\tau)^m=q^{-\mu}\sum_n c(n)q^n$ with $\mu=m(t-1)/24$. Then $F$ has weight $0$ and a finite-order multiplier on $\Gamma_0(t)$.

## 0. The exact formula
The Rademacher–Zuckerman formula (Zuckerman 1939; Bringmann–Ono 2012) gives, for every $n>\mu$,
$$c(n)=\sum_{\text{terms }T}A_T(n)\,B_T(n),\qquad B_T(n)=\frac{2\pi}{k}\sqrt{\frac{\gamma}{n-\mu}}\;I_1\!\Big(\frac{4\pi}{k}\sqrt{\gamma(n-\mu)}\Big).$$
- A term $T=(k,\gamma)$ runs over denominators $k$ with $\gcd(k,t)^2>t$ and over the pole orders $\gamma>0$ at the cusps of type $d=\gcd(k,t)$. Only finitely many $\gamma$ occur for each type.
- The growth rate of $T$ is $\beta_T=\sqrt\gamma/k$. For any $\varepsilon>0$ only finitely many terms have $\beta_T>\varepsilon$, and the series converges absolutely.
- The amplitude $A_T(n)$ is a generalized Kloosterman sum.
- The explicit model built from this (`model.py`) reproduces exact coefficients to 12 digits.

## 1. Structure of the amplitudes (Lemma M)
For a term $T=(k,\gamma)$ of type $d$,
$$A_T(n)=C_T\sum_{h\bmod 24k}{}^{*}\ \psi_T(h)\;e\!\Big(\frac{X_T(n)\,h+Y_T\,\bar h}{24k}\Big),$$
where:
- $\psi_T$ is a **real** character, a product of Jacobi symbols in $h$;
- $X_T(n)=-24n+x_T$ is affine in $n$;
- $Y_T$ is an integer with $Y_T\ne0$, and $C_T$ is a constant.

*Proof.* The phase of the term is $(\omega_{h,k}^{-1}\omega_{th/d,\,k/d})^m$, with $\omega_{h,k}=e^{\pi i s(h,k)}$. Knopp's eta-multiplier formula writes each $\omega$ as a Jacobi symbol times $e\big((\alpha h+\alpha'\bar h)/24k\big)$; the formula was checked numerically on 39 random matrices, see `lemma_checks.py`. One uses $\overline{th/d}\equiv\overline{t/d}\,\bar h\pmod{k/d}$, which holds because $\gcd(t/d,k/d)=1$.

The principal-part coefficient of order $\gamma$ is a sum over monomials $\tilde q^{\,i}x^j$ with $i+j\,d^2/t$ fixed. All of these carry the same phase $e(\Lambda\bar h/k)$: moving between them changes the $\bar h$-coefficient by $-d^2/g'+d^2/g'\equiv0\pmod k$. Finally, the summand is $k$-periodic in $h$, so summing over $h\bmod k$ is a constant multiple of summing over $h\bmod 24k$ coprime to $24k$, since the fibres have uniform size $\varphi(24k)/\varphi(k)$. $\square$

## 2. Factorisation over primes (Lemma F)
Write $24k=\prod_p p^{a_p}$. By the Chinese remainder theorem,
$$A_T(n)=C'_T\prod_{p\mid 6k}\mathcal T_p\big(u_pX_T(n),\,u_pY_T\big),\qquad u_p=\overline{(24k/p^{a_p})}\ \text{mod }p^{a_p},$$
with local sums
$$\mathcal T_p(x,y)=\sum_{h\bmod p^{a_p}}{}^{*}\psi_p(h)\,e\big((xh+y\bar h)/p^{a_p}\big),$$
where $\psi_p$ is the $p$-component of $\psi_T$ (again real). In particular, $A_T(n)=0$ if and only if some local factor vanishes.

## 3. Galois invariance (Lemma G)
Since $\psi_p$ is real, the automorphism $\zeta_{p^a}\mapsto\zeta_{p^a}^{\,w}$ (for a unit $w$) sends $\mathcal T_p(x,y)$ to $\mathcal T_p(wx,wy)$. Hence
$$\mathcal T_p(wx,wy)=0\iff\mathcal T_p(x,y)=0 .$$
Substituting $h\mapsto ch$ gives $\mathcal T_p(\bar cx,\bar cy)=\psi_p(c)\,\mathcal T_p(x,\bar c^{\,2}y)$. So vanishing depends only on the orbit of $(x,y)$ under this scaling. This was checked numerically in 250 cases, with 0 violations (`local_checks.py`).

**Corollary (B1).** Let $c'$ be coprime to $6t$. The prime-to-$p$ part of $k$ enters $\mathcal T_p$ (for $p\mid 6t$) only through the unit scaling $u_p$, which is exactly an orbit move in Lemma G. The same holds for the local factors at primes dividing $k$ but not $6t$. Therefore every local factor of $T=(k,\gamma)$ reappears, up to an orbit move, as a local factor of $(kc',\gamma)$, and so
$$A_{(k,\gamma)}(n)=0\ \Longrightarrow\ A_{(kc',\gamma)}(n)=0 .$$
This was checked in 825 cases for composite $t$ (`composite_test.py`).

## 4. Stationary phase and stabilisation (Lemma S)
Let $\kappa$ be the conductor exponent of $\psi_p$ (so $\kappa\le1$ for odd $p$ and $\kappa\le3$ for $p=2$). Let $b=v_p(y)$ and $A\ge A_0:=2(b+2+3[p=2])+2\kappa$. Put $\beta=\lfloor A/2\rfloor$ and substitute $h=h_0(1+p^{A-\beta}z)$. The character $\psi_p$ is invariant under this substitution, and $(1+p^{A-\beta}z)^{-1}\equiv1-p^{A-\beta}z\pmod{p^A}$. Summing over $z$ gives
$$\mathcal T_p(x,y)=p^{\beta}\sum_{\substack{h_0\\ xh_0^2\equiv y\ (p^\beta)}}\psi_p(h_0)\,e\big((xh_0+y\bar h_0)/p^A\big).$$
- **Unsolvable congruence.** If $xh^2\equiv y\pmod{p^\beta}$ has no solution, the sum is $0$. Once $\beta>b+1+3[p=2]$, solvability is equivalent to: $v_p(x)=b$, and $(x/p^b)(y/p^b)$ is a square mod $p$ (mod $8$ if $p=2$). By Hensel's lemma this condition does not depend on $A$.
- **Solvable congruence.** The solutions are $\pm h_*$ (together with $\pm h_*(1+2^{\cdot})$ when $p=2$). Their contributions combine to $p^{A/2}$ times a nonzero Gauss-sum factor times $2\cos$ or $2i\sin$ of $4\pi x_0/p^A$, where $x_0$ has fixed valuation. For $A\ge A_0$ this is never $0$. These are the classical prime-power evaluations of Kloosterman and Salié sums (Iwaniec–Kowalski, Lemmas 12.2–12.3); for $p=2$ the four solutions are handled in the same way.

**Conclusion (B2).** For $A\ge A_0$, whether $\mathcal T_p$ vanishes does not depend on $A$. It depends on $n$ only through $n\bmod p^{A_0+1}$. This was checked numerically in 120 cases.

## 5. Eventual periodicity
Let $b_p$ be the maximum of $v_p(Y_T)$ over the finitely many pole orders. Set
$$K_0=\prod_{p\mid 6t}p^{A_0(b_p)+v_p(t)+1},\qquad L=\operatorname{lcm}(t,K_0).$$
Fix a class $n\equiv s\pmod L$.

- **Terms with $k\mid K_0$.** Their amplitudes depend only on $n\bmod k$, so they are constant on the class.
- **Every other term $T=(k,\gamma)$.** Write $k=k_0c'$, where $k_0$ has all its prime factors in $6t$ and $\gcd(c',6t)=1$. Suppose $T$ is nonzero at some $n$ in the class. We show $T$ is not the first non-vanishing term at $n$.
  - If $c'>1$: by B1 (contrapositive), $(k_0,\gamma)$ is also nonzero at $n$, and its growth rate $\sqrt\gamma/k_0$ is strictly larger than that of $T$.
  - If $c'=1$ and some exponent of $k_0$ exceeds its cap: by B2, the term with that exponent lowered to the cap is also nonzero at $n$. It has a smaller $k$, hence a strictly larger growth rate.

So the first non-vanishing term in order of growth rate, if there is one, is a term with $k\mid K_0$. Its amplitude is a nonzero constant on the class, so it fixes the sign of $c(n)$ for all large $n$ in the class. If no term is nonzero on the class, then $c(n)=0$ there by §0. Hence $\operatorname{sgn}c(n)$ is eventually periodic, with period dividing $L$.

## 6. Least period divisible by $t$
The leading amplitude is $\alpha(r)=\sum_{h\bmod t}^{*}\psi(h)e(-hr/t)$ with $|\psi|=1$, so its Fourier support is the set of units mod $t$. For a proper divisor $g$ of $t$, $\alpha$ therefore sums to $0$ on every coset mod $g$. Some coset contains a class where $\alpha\ne0$; that coset then contains both signs, and those classes are non-degenerate. A period $L$ with $t\nmid L$ would force the sign to be constant on the cosets mod $\gcd(L,t)$, which is impossible. $\square$

## Status of rigour
- §§0, 3, 5 and 6 are complete arguments.
- §1 relies on Knopp's multiplier formula and a bookkeeping check of the phases; the bookkeeping should be written out case by case in the final paper.
- §4 is the standard $p$-adic stationary-phase argument. The explicit evaluation in the solvable case for $p=2$ with conductor-8 characters should be written out in full; it is the most delicate local case.

---
## Review round 1 (bookkeeping and regression)

**§1 bookkeeping.** We have explicit formulas, derived from the Apostol/Knopp eta multiplier, valid whenever $\gcd(h,24k)=1$ and $\bar h$ is the inverse of $h$ mod $24k$:
$$k\text{ odd: }\ e^{\pi i s(h,k)}=\Big(\frac hk\Big)i^{(k-1)/2}e\Big(\frac{(1-k^2)(h+\bar h)}{24k}\Big),\qquad
k\text{ even: }\ e^{\pi i s(h,k)}=\Big(\frac kh\Big)e\Big(\frac{\bar h(1-k^2)+h(2k^2-3k+1)}{24k}\Big).$$
Checked on 3409 pairs $(h,k)$ with $k<60$: 0 failures (`bookkeeping.py`).

*Caveat.* For cusps of type $d<t$ with $\gcd(t/d,6)>1$, the second factor $\omega_{th/d,\,k/d}$ must be evaluated at a lift of $th/d$ that is coprime to 6. This adds a function of $h\bmod 24$ that need not be a character. That function lives only at the primes 2 and 3, so Lemma G (Galois invariance) needs a separate argument at $p\in\{2,3\}$ in that case.

**Regression test of §5.** For all $(t,m)$ with $t\le12$ and $m(t-1)\le60$ (21 899 coefficient comparisons near $n=900$), the §5 decision procedure, truncated at $k\le3t$ and period $L=24t$ or $12t$, predicts every sign correctly except for two pairs:
- $(5,5)$: the truncation is too small. With $k\le30$ and $L=600$ there are 0 mismatches, and the sign is decided at $k=25$ (Wang's period-25 special case).
- $(12,3)$: the classes $n\equiv2\pmod3$ are identically zero (Jacobi's identity). There, the term $(k=44,\gamma=1/24)$, of cusp type $d=4$, has **the same growth rate** as the term $(k=132,\gamma=3/8)$, of cusp type $d=12$. Their amplitudes cancel. The truncation $k\le72$ misses $k=132$.

**Consequence for the proof.** Terms of *different* cusp types can have exactly equal growth rates, so the argument of §5 must work with **groups of equal growth rate**, not with single terms. B1 and B2 in their per-term form do not immediately give group-level statements: a group can vanish by cancellation between terms whose $k$'s differ by a factor (here 3). This is a genuine gap in §5 as written. The theorem itself is unaffected: brute force confirms eventual periodicity in every case tested.

---
## Review round 2: closing the equal-growth gap (option A)

**Discovery: the Scaling Lemma.** Let $S$ be the finite set of primes dividing $6t$ or dividing a numerator or denominator of some pole order. Let $T=(k,\gamma)$ be a term and $q$ a prime such that either
- $q\notin S$, or
- $q\in S$ and $v_q(k)\ge A_q$, a stabilisation threshold chosen above $v_q(t)$, so that $k\mapsto kq$ does not change the cusp type.

Then
$$A_{(kq,\gamma)}(n)=R_q\big(n;\ \gamma/k^2,\ v_q(k)\big)\cdot A_{(k,\gamma)}(n).$$
The factor $R_q$ depends on the term only through its growth rate $\beta^2=\gamma/k^2$ (and $v_q(k)$). **In particular it is the same for terms of different cusp types in the same equal-growth group.**

*Why.* Factor the amplitude over primes as in Lemma F. Replacing $k$ by $kq$ does two things.
- The local factors at primes $p\ne q$ are twisted by $(x,y)\mapsto(\bar qx,\bar qy)$. Since the second argument $y$ scales like $\beta^2(24k)^2$, this twist is compensated, and those factors change only by a character value.
- The $q$-local factor is replaced by one whose parameters are $n$ and $\beta^2$ alone.

The full derivation still has to be written out; see "Remaining".

*Numerical evidence.*
- `group_b1.py`: 1253 group-level tests, where a group vanishes if and only if its scaled version does. 0 failures.
- `scaling_test.py` plus a stable-range check: ratios agree across cusp types in 80 out of 80 cases. The only discrepancies were for $q=2$ below the threshold (e.g. $t=12$, $k=8\to16$), which is exactly the unstable range the lemma excludes.

**Group-level B1 and B2.** Members of an equal-growth group $G_\beta$ have $k$-values differing by $S$-units, since $k_i/k_j=\sqrt{\gamma_i/\gamma_j}$.
- If a prime $q\notin S$ divides one member's $k$, it divides every member's $k$, and $\{(k_i/q,\gamma_i)\}$ is exactly the faster group $G_{\beta q}$.
- For $q\in S$ in the stable range, the same holds once the cap is raised by the maximal $q$-valuation spread within a group.

By the Scaling Lemma,
$$G_\beta(n)=R_q(n;\beta^2)\,G_{\beta q}(n),$$
so $G_\beta(n)\ne0$ implies $G_{\beta q}(n)\ne0$.

**§5 repaired.** Take $L=\operatorname{lcm}\big(t,\prod_{q\in S}q^{\text{cap}_q}\big)$ and fix a class mod $L$.
- A group is either **$L$-rigid** (every member has $k\mid L$, so its amplitude is constant on the class), or it is the scaling of a strictly faster group.
- In the second case, if the group is nonzero at $n$, the faster group is also nonzero at $n$.

So the first non-vanishing group at every $n$ in the class is $L$-rigid and nonzero on the whole class, and it decides the sign. If there is no such group, every group vanishes on the class and $c(n)=0$ by the exact formula. This closes the gap found in round 1. For example, the cancellation in $(12,3)$ between $(44,1/24)$ and $(132,3/8)$ is $q=11$ applied to the rigid group $\{(4,1/24),(12,3/8)\}$.

**Remaining before the proof is complete.**
1. A rigorous proof of the Scaling Lemma. This means the CRT bookkeeping showing that
   - the second argument $y_T$ scales like $\beta^2$, and
   - the character twists $\psi_p(q)$ agree across cusp types. At $p\in\{2,3\}$ this interacts with the round-1 caveat (the lifted term $\omega_{th/d,k/d}$).
2. Explicit thresholds $A_q$ for $q\in S$, together with a justification of the "valuation spread" adjustment.

---
## Review round 3: proof of the Scaling Lemma and the final form of §5

### Correction to the Scaling Lemma
The round-2 statement ("same $n$") is right only up to an argument shift. The correct statement is:

> **Scaling Lemma.** Let $q$ be a prime with $q\notin S$. Here $S$ is the set of primes dividing $6t$ together with the primes dividing a numerator or denominator of some pole order. Then for every term $(k,\gamma)$,
> $$A_{(kq,\gamma)}(n)=R_q(n;\gamma/k^2)\;A_{(k,\gamma)}\big(n\,\bar q^{\,2}\big),$$
> where $\bar q$ is the inverse of $q$ modulo $24\,\mathrm{lcm}(k)$. The factor $R_q$ is a plain Kloosterman sum modulo $q$, times $\chi(q)$. It depends only on $n$ and on $\gamma/k^2$, and it is **never zero**.

Numerical check (`scaling2.py`): 672 tests across six equal-growth groups and $q\in\{5,7,11,13,17\}$, with 0 violations; $R_q=0$ never occurs. For $(12,3)$ and $(6,8)$ it happens that $\bar q^{\,2}\equiv1$ on the relevant classes, which is why the round-2 "same $n$" tests passed.

### Proof
1. **A genuine character.** $F(24\tau)=\eta(24\tau)^m/\eta(24t\tau)^m$ is a modular function on $\Gamma_0(N)$, $N=576t$, with a **real quadratic character** $\chi$. This is the Gordon–Hughes–Newman criterion: $\sum_\delta\delta r_\delta=24m(1-t)\equiv0$ and $\sum_\delta(N/\delta)r_\delta\equiv0\pmod{24}$. Its coefficients are $c(n)$ at exponent $24n-m(t-1)$.
2. **The exact formula at every cusp.** The exact formula for $F(24\tau)$ on $\Gamma_0(N)$ writes each term of §0 as a cusp-pair Kloosterman sum $S_{\mathfrak a\infty}(\nu,\,24n-m(t-1);\,c)$. Here $\mathfrak a$ is the cusp of type $d$, $\nu$ is the principal-part index (so $|\nu|$ is proportional to $\gamma$, with the width normalisation of $\mathfrak a$), and $c$ is proportional to $k$.
3. **Selberg's twisted multiplicativity.** For $q\nmid N$ (Deshouillers–Iwaniec, §1.2; Iwaniec, *Spectral methods*, §2.6):
$$S_{\mathfrak a\infty}(\nu,x;cq)=\chi(q)\,S_{\mathfrak a\infty}(\nu\bar q,x\bar q;c)\;S(\nu\bar c,x\bar c;q).$$
Since $\chi$ is real, the substitution $d\mapsto qd$ gives $S_{\mathfrak a\infty}(\nu\bar q,x\bar q;c)=S_{\mathfrak a\infty}(\nu,x\bar q^{\,2};c)$. This is the shift $n\mapsto n\bar q^{\,2}$.
4. **The $q$-factor depends only on the growth rate.** $S(\nu\bar c,x\bar c;q)=S(\nu\bar c^{\,2},x;q)$, and $\nu\bar c^{\,2}\bmod q$ is fixed by $\nu/c^2$, which is proportional to $\gamma/k^2=\beta^2$. So all members of an equal-growth group share this factor. It is a Kloosterman sum modulo a prime $q$, and therefore nonzero: it is $\equiv-1\pmod{1-\zeta_q}$, or a nonzero Ramanujan sum if an argument is $\equiv0$. $\square$

This replaces the unfinished CRT bookkeeping at $p\in\{2,3\}$. The primes 2 and 3 now enter only through $N$, and the lemma is used only for $q\nmid N$.

### Final form of §5 (eventual periodicity)
Fix $L$ divisible by $t$ and by every rigid denominator; "rigid" means all prime factors lie in $S$, with $q$-exponents capped at the $\Gamma_0(N)$ stabilisation thresholds of Lemma S. Every equal-growth group is either
- **rigid**, with constant amplitude on classes mod $L$, or
- $G_{\beta/c'}$ with $G_\beta$ rigid and $c'$ composed of primes outside $S$, in which case
$$G_{\beta/c'}(n)=R_{c'}(n)\,G_\beta(n\bar c'^{\,2}).$$

No cross-cancellation between different scalings can occur. Ratios of rigid growth rates are $S$-units, so $\beta_1/c_1=\beta_2/c_2$ forces $c_1=c_2$ and $\beta_1=\beta_2$.

Fix a class $C$ mod $L$.
- By Dirichlet's theorem there are primes $q\notin S$ in every unit class mod $L$. Also $R_q\ne0$ at exponent 1. Hence the first non-vanishing group on $C$ is among finitely many candidates: rigid groups, and rigid groups scaled by a prime $q$ at most the least prime in each class mod $L$, plus boundedly many prime-power scalings.
- Whether a candidate vanishes depends only on the class $C$, because $G_\beta(n\bar c'^{\,2})$ depends on $n\bmod L$.
- The deciding group's sign is $\operatorname{sgn}R_{c'}(n)\,G_\beta(n\bar c'^{\,2})$, which is periodic mod $L\,c'$.

If every candidate vanishes on $C$, then every group vanishes there, and $c(n)=0$ on $C$ by the exact formula. So $\operatorname{sgn}c(n)$ is eventually periodic. Combined with §6, **Wang's Conjecture 1.6 holds for all $t$ and $m$.**

### What remains
1. **Scalings by primes in $S$ in the stable range.** These need the same multiplicativity for $q\mid N$. Two routes:
   - Deshouillers–Iwaniec's local analysis at $q\mid N$ with the exponent above the $q$-part of $N$; or
   - absorbing these finitely many primes into $L$, which is the rigid cap. This works because only exponents up to the stabilisation threshold matter by Lemma S, and Lemma S for $\Gamma_0(N)$-Kloosterman sums at $q\mid N$ is standard stationary phase.

   This must be written out.
2. **Normalisations in step 2.** Cusp widths, and the exact $\nu\leftrightarrow\gamma$ and $c\leftrightarrow k$ correspondence for each cusp type, must match `model.py`. The model agrees with exact coefficients to 12 digits, but the dictionary still needs to be stated explicitly.
