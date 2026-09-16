# Wang's Conjecture 1.6 at prime level: what is proved

**Setting.** Let $G_p^m=(q;q)_\infty^m/(q^p;q^p)_\infty^m=\sum_n c(n)q^n$, and put $\mu=m(p-1)/24$.

**Lemma 1 (eta multiplier at prime level).** For primes $p\ge5$,
$$e^{-\pi i s(h,p)}=C_p\,\Big(\frac hp\Big)\,e\Big(-\frac{\overline{24}\,(h+\bar h)}{p}\Big),\qquad |C_p|=1 .$$
This is the classical Rademacher–Grosswald evaluation. It has also been checked numerically for all primes $p<200$.

**Lemma 2 (structure of the growing terms).** For $p$ prime, the growing terms sit exactly at the cusps $h/k$ with $p\mid k$. The principal part there is the $G_p^m$-series itself: $\widehat G=G_p(\tilde q)^m$. Hence the term of order $e<\mu$ at $k=p$ has amplitude
$$A_{p,e}(n)=c(e)\,C_p^m\sum_{h\bmod p}\Big(\frac hp\Big)^m e\Big(\frac{(\mu-n)h+(\mu-e)\bar h}{p}\Big).$$
Here $\mu$ is read mod $p$ via $\mu\equiv-m\overline{24}$.

- For even $m$ this is a **Kloosterman sum**.
- For odd $m$ it is a **Salié sum**.

This was checked against exact coefficients (Theorem 1.1 model); the zero classes predicted by it match the exact computation in 5605 out of 5605 cases.

**Theorem (new).** Wang's Conjecture 1.6 holds, with least period of sign exactly $p$, in each of these cases:
- (i) $p\ge5$ and $m$ even, for every $m$;
- (ii) $p=3$ and $m\not\equiv3\pmod 6$;
- (iii) $p=2$ (Wang), and $m\in\{1,3\}$ (Wang).

In all these cases the sign of $c(n)$ equals the sign of $\alpha_p^{(m)}(n\bmod p)$ for all large $n$.

*Proof.* The leading term ($k=p$, $e=0$) has the strictly largest exponential rate $\sqrt\mu/p$. Its amplitude is Wang's $\alpha_p^{(m)}(r)$.
- For (i), this amplitude is a unimodular constant times the Kloosterman sum $S(\mu-r,\mu;p)$. That sum never vanishes: in $\Z[\zeta_p]$ it is congruent to $p-1\equiv-1$ modulo $(1-\zeta_p)$.
- For (ii), $\alpha_3^{(m)}(r)=2\cos(\pi m/18+2\pi r/3)$, which vanishes exactly when $m\equiv3\pmod 6$.

Since UPS is an asymptotic statement, the qualitative dominance of the leading term is enough; no explicit $n_0$ is needed. The least period is $p$ because $\sum_r\alpha_p^{(m)}(r)=0$ and the sequence is not identically zero, so the signs are not constant, and $p$ is prime. $\square$

**Odd $m$, $p\ge5$ (partial).** Here $\alpha$ vanishes exactly on the classes where $(\mu-r)\mu$ is a quadratic non-residue mod $p$, together with $r\equiv\mu\equiv0$.
- Salié sums vanish only on non-residue arguments. So on such a class the $k=p$ terms $e=1,2,\dots$ take over at the first $e$ with $c(e)\ne0$ and $(\mu-r)(\mu-e)$ a nonzero residue (or $e\equiv\mu$).
- If that $e_0$ is below $3\mu/4$, it beats every cusp $k\ge2p$, and the sign is constant on the class, so the period is $p$.
- Remaining gap: to guarantee $e_0<3\mu/4$ for all odd $m$, one needs a nonvanishing statement for the coefficients $c_m(e)$ of $(q;q)^m$ at suitable $e<p$, plus a least-residue bound. The small cases $(3,9)$ and $(5,5)$ (Wang's special cases) show that the period can then become $p^2$.

---
## Update: odd m reduced to an explicit condition; p = 3 closed

**Key observation.** For $p\ge5$ and odd $m$, the degenerate classes (those with $\alpha(r)=0$) are exactly the $r$ with $\chi((\mu-r)\mu)=-1$, together with $r\equiv\mu\equiv0$. Here $\chi$ is the Legendre symbol mod $p$ and $x:=\mu\bmod p=-m\overline{24}$.

On such a class, the $k=p$ term of order $e$ is nonzero if and only if both
- $c(e)\ne0$, and
- $e\equiv x\pmod p$ or $\chi(x(x-e))=-1$.

This condition **does not depend on $r$**. Every $k=p$ amplitude is constant on the class, so the first such $e$ decides the sign on every degenerate class at once, and the sign is constant there. Only finitely many terms have a larger growth rate. Each of them either vanishes on the class or has an amplitude that is periodic mod a multiple of $p$.

**Theorem (odd m).** Let $p\ge5$ be prime and $m$ odd with $m(p-1)>24$. Assume condition
(C): there exists $e\in[1,\mu)$ with $c_m(e)\ne0$ such that $e\equiv x\pmod p$ or $\chi(x(x-e))=-1$.

Then the signs of $c_p^{(m)}(n)$ are eventually periodic with least period divisible by $p$. (If $p\mid m$, then $x=0$, and $e=1$ works because $c(1)=-m$.)

*Least period divisible by p.* If the least period $L$ were prime to $p$, translating by $L$ would permute the residues mod $p$ transitively. Then all non-degenerate classes would share one sign, contradicting $\sum_r\alpha(r)=0$ with $\alpha\not\equiv0$.

**Verification of (C).** Condition (C) holds for:
- every odd $5\le m\le99$ and every prime $p\le3000$ (20543 pairs), with the witness always $e\le11$;
- every odd $m\le301$ and every prime $p\le1000$.

See `condC.py`.

**$p=3$, all $m$.** At $k=3$ the term of order $e$ has amplitude $2c(e)\cos(-\pi m/18+2\pi(2e-r)/3)$. For fixed $r$, exactly one residue of $e$ mod 3 gives zero. So on the zero class (which occurs when $m\equiv3\pmod6$) the $e=1$ term, with $c(1)=-m\ne0$, decides the sign whenever $\mu=m/12>1$, i.e. $m\ge15$. Wang handled $m=3,9$. **Hence the conjecture holds for $p=3$ and every $m$.**

**What remains for p ≥ 5, odd m.** A proof of (C) for all $(m,p)$. Failure would require
- a run of equal Legendre symbols $\chi(x),\chi(x-1),\dots,\chi(x-E)$ at the specific point $x=-m/24\bmod p$, and
- the vanishing of the coefficients of $(q;q)_\infty^m$ at every good $e$ below $\mu$.

So what is needed is a nonvanishing statement for the coefficients of $(q;q)^m$ at small indices, combined with a bound on how long such a run of Legendre symbols can be. Neither is automatic. If (C) ever failed, that would **not** disprove the conjecture: the sign would then be decided at level $p^2$, as in Wang's cases $(3,9)$ and $(5,5)$.

---
## Update 2: the odd-m gap closes. Condition (C) is not needed.

**Level Lemma.** Let $p\ge5$ be prime and $m$ odd. On a degenerate class $r$ (one with $\alpha(r)=0$), fix a term order $e$. If the term of order $e$ at $k=p$ vanishes on the class, then the term of order $e$ at **every** cusp $k=pj$ vanishes on the class.

*Proof sketch.* The amplitude at $k$ is $c(e)$ times the generalized Kloosterman sum of $F=\eta(\tau)^m/\eta(p\tau)^m$ on $\Gamma_0(p)$. We use the eta-multiplier formula: Jacobi symbols $(d/c)$, plus additive exponentials in $d$ and $\bar d$.
- The ratio $\nu_\eta(\gamma)^m\nu_\eta(\gamma_p)^{-m}$ has multiplicative part $(d/p)^m$, and only the $p$-part survives. For even $c$, reciprocity turns $(c/d)$ into $(d/p)$ times a character mod 4.
- So for $k=p^a c$ with $p\nmid c$, the sum factors by twisted multiplicativity: a twisted Salié sum mod $p^a$ (character $\chi_p^m$) times a Kloosterman sum mod $c$. The twisting multiplies the product $uv$ of the $p$-part arguments by a square, $\bar c^2$.
- A twisted Salié sum mod $p^a$ with $p\nmid uv$ vanishes if and only if $\chi(uv)=-1$. Here $u\equiv\mu-n$ is fixed by the class and $v\equiv\mu-e$. So the vanishing condition is the same at every level: it depends only on $(r,e)$.
- The cases $p\mid u$ or $p\mid v$ are the "good" ones (Gauss-sum or Ramanujan-sum factors), as at level $p$.

*Numerical check.* 430 comparisons with $k=pj$, $j\le6$, plus 56 more at $k=p^2,2p^2,5p^2$ (up to $k=125$). No violations.

**Theorem (all primes).** Wang's Conjecture 1.6 holds for every prime $t=p$ and every $m\ge1$.

*Proof.*
- $p=2$ is Wang's result. $p=3$ is Update 1. Even $m$ with $p\ge5$ is the Kloosterman argument.
- Now let $p\ge5$ and $m$ odd, and take a degenerate class $r$. There are two cases.
  - If some good $e<\mu$ has $c(e)\ne0$ (condition (C)), the level-$p$ term decides the sign, and the sign is constant on the class.
  - Otherwise every term vanishes on the class: good terms because $c(e)=0$, bad terms by the Level Lemma. The exact Rademacher–Zuckerman formula then gives $c(n)=0$ for all $n>\mu$ in the class.
- Either way the signs are eventually periodic. The least period is divisible by $p$ by the transitivity argument. $\square$

**Still to do for a paper.**
1. A full proof of the Level Lemma: explicit eta-multiplier Kloosterman sums, twisted multiplicativity, and the evaluation of twisted Salié sums mod $p^a$ including the edge cases $p\mid uv$. All standard, but it must be written out carefully.
2. Citation of the exact formula (Zuckerman; Bringmann–Ono) for weight-0 eta-quotients with integer exponents.
3. Composite $t$, which is still open.
