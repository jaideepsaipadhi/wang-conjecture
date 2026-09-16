# Wang's Conjecture 1.6 for prime level: proof of the Level Lemma

## Setting
Let $p\ge5$ be prime and let $m\ge1$ be odd. Put
$$F(\tau)=\frac{\eta(\tau)^m}{\eta(p\tau)^m}=q^{-\mu}\sum_{n\ge0}c(n)q^n,\qquad \mu=\frac{m(p-1)}{24}.$$
Then $F$ has weight $0$ on $\Gamma_0(p)$. For $M=\begin{psmallmatrix}a&b\\c&d\end{psmallmatrix}\in\Gamma_0(p)$ its multiplier is
$$\nu_F(M)=\Big(\frac{\nu_\eta(M)}{\nu_\eta(M_p)}\Big)^{m},\qquad M_p=\begin{pmatrix}a&pb\\ c/p&d\end{pmatrix}.$$
This holds because $\eta(pM\tau)=\eta(M_p(p\tau))$ and the automorphy factors $(c\tau+d)^{m/2}$ cancel.

The only pole of $F$ is at the cusp $\infty$ of $\Gamma_0(p)$, and its principal part is $\sum_{0\le e<\mu}c(e)q^{e-\mu}$. The Rademacher–Zuckerman exact formula (Zuckerman 1939; Bringmann–Ono) gives, for $n>\mu$,
$$c(n)=\sum_{0\le e<\mu}c(e)\sum_{\substack{k\ge1\\ p\mid k}}\frac{2\pi}{k}\sqrt{\frac{\mu-e}{n-\mu}}\;K_k(e,n)\;I_1\Big(\frac{4\pi}{k}\sqrt{(\mu-e)(n-\mu)}\Big).$$
Here
$$K_k(e,n)=\sum_{d\bmod k}{}^{*}\ \overline{\nu_F(M_d)}\;e\Big(\frac{(e-\mu)a+(n-\mu)d}{k}\Big),\qquad M_d=\begin{pmatrix}a&b\\k&d\end{pmatrix},\ ad\equiv1\pmod k .$$
The explicit model `model.py` reproduces the exact coefficients to 12 digits for composite and prime $t$, so the normalisations agree.

## Lemma 1 (the $p$-part of the multiplier is a Legendre symbol)
We use Knopp's formula for the eta multiplier (Knopp, *Modular Functions in Analytic Number Theory*, Ch. 4, Thm. 2). It was checked numerically on 39 random matrices, with 0 failures (`lemma_checks.py`). For $c>0$:
$$\nu_\eta(M)=\begin{cases}\big(\tfrac dc\big)\,e\!\left(\frac{(a+d)c-bd(c^2-1)-3c}{24}\right), & c\text{ odd},\\[4pt]
\big(\tfrac cd\big)\,e\!\left(\frac{(a+d)c-bd(c^2-1)+3d-3-3cd}{24}\right), & c\text{ even}.\end{cases}$$

**Claim.** $\nu_F(M)=\big(\tfrac dp\big)^m\,\rho(M)$, where $\rho(M)$ has trivial $p$-part. Precisely: write $k=c=p^\alpha c'$ with $p\nmid c'$. Then $24k\cdot\frac{1}{2\pi i}\log\rho(M)$ is an integer divisible by $p^\alpha$.

*Proof.*
- **Character part.** For odd $c$, $\big(\tfrac dc\big)\big/\big(\tfrac d{c/p}\big)=\big(\tfrac dp\big)$. For even $c$, reciprocity turns $\big(\tfrac cd\big)\big/\big(\tfrac{c/p}d\big)=\big(\tfrac pd\big)$ into $\big(\tfrac dp\big)$ times a character mod $4$, which has trivial $p$-part.
- **Additive part.** Put $j=c/p$. The two additive exponents differ by
  $$\frac{p-1}{24}\Big[(a+d)j-bd(pj^2+1)-3j\Big]\qquad(\text{plus } -3dj(p-1)/24 \text{ when } c \text{ is even}).$$
  Multiply by $24k=24c$. The terms become $c^2(a+d)/p$, $c\cdot bd(pj^2+1)$, $3cj$ and $3dc^2/p$ (the last only for even $c$). Each is divisible by $p^{2\alpha-1}$ or by $c$, hence by $p^\alpha$. $\square$

## Lemma 2 (factorisation)
Every term of the exponent of $K_k$, times $24k$, is an integer function of $d\bmod 24k$. Split $e(X/24k)$ by the Chinese remainder theorem into its parts mod $p^\alpha$ and mod $24c'$; this is possible because $p\ge5$. By Lemma 1, the multiplier contributes only $\big(\tfrac dp\big)^m$ to the $p^\alpha$-part. Hence
$$K_k(e,n)=T_\alpha(u,v)\cdot R,\qquad T_\alpha(u,v)=\sum_{d\bmod p^\alpha}{}^{*}\Big(\frac dp\Big)^m e\Big(\frac{ud+v\bar d}{p^\alpha}\Big).$$
Here
$$u\equiv (n-\mu)\,\overline{c'},\qquad v\equiv (e-\mu)\,\overline{c'}\pmod{p^\alpha},$$
where $n-\mu$ and $e-\mu$ are read in $\mathbb Z_{(p)}$ (their denominators divide 24), and $R$ depends only on the part mod $24c'$. In particular
$$uv\equiv (n-\mu)(e-\mu)\,\overline{c'}^{\,2}\pmod{p^\alpha},$$
so the Legendre symbol $\chi(uv)=\chi\big((n-\mu)(e-\mu)\big)$ does not depend on the level $k$.

## Lemma 3 (twisted Salié sums)
Let $m$ be odd and suppose $p\nmid uv$ and $\chi(uv)=-1$. Then $T_\alpha(u,v)=0$ for every $\alpha\ge1$.

*Proof.*
- **$\alpha=1$.** Salié's evaluation gives $T_1(u,v)=\varepsilon_p\sqrt p\,\chi(v)\sum_{x^2\equiv4uv}e(x/p)$, which is an empty sum.
- **$\alpha\ge2$.** Let $\beta=\lfloor\alpha/2\rfloor\ge1$. Substitute $d\mapsto d(1+p^{\alpha-\beta}y)$ and sum over $y\bmod p^\beta$. The character $\chi$ is invariant under this substitution. The sum over $y$ vanishes unless $ud\equiv v\bar d\pmod{p^\beta}$, i.e. $ud^2\equiv v$, which forces $\chi(uv)=1$. $\square$

This was checked numerically on 408 cases with $p\le13$ and $\alpha\le3$: the sum vanishes exactly when $\chi(uv)=-1$.

## Level Lemma
Let $r$ be a degenerate class and $e$ an order with $p\nmid(r-\mu)(e-\mu)$. If the level-$p$ amplitude of order $e$ vanishes on the class (equivalently, $\chi((r-\mu)(e-\mu))=-1$), then $K_k(e,n)=0$ for every $k$ with $p\mid k$ and every $n\equiv r\pmod p$.

*Proof.* Combine Lemmas 2 and 3. $\square$

**Excluded case $p\mid(r-\mu)(e-\mu)$.** Here the level-$p$ amplitude is a Gauss or Ramanujan sum, nonzero, unless $r\equiv e\equiv\mu\equiv0\pmod p$. That can only happen when $p\mid m$. In that situation $e=1$ is good, since $c(1)=-m\ne0$ and $\chi$ is nonzero at the relevant arguments, so the class is decided at level $p$ without the Level Lemma.

## Theorem (all prime levels)
For every prime $t$ and every integer $m\ge1$, the signs of the coefficients of $(q;q)_\infty^m/(q^t;q^t)_\infty^m$ are eventually periodic, with least period divisible by $t$.

*Proof.* The case $t=2$ is Wang's; $t=3$ and $t\ge5$ with $m$ even are in `prime_case_note.md`. For $t=p\ge5$ and $m$ odd, take a degenerate class. There are two cases.
- **Some good order $e$ has $c(e)\ne0$.** Its level-$p$ term is constant on the class. It is the first non-vanishing term in growth order: every term with a larger growth rate is either a level-$p$ term that vanishes on the class, or, by the Level Lemma, a higher-level term that vanishes there too. So it decides the sign.
- **No good order has $c(e)\ne0$.** Then every term vanishes on the class, and the exact formula gives $c(n)=0$ there.

Finally, suppose the least period $L$ were prime to $p$. Then all non-degenerate classes would have the same sign, contradicting $\sum_r\alpha(r)=0$. $\square$

## Remaining points for the paper
- State the exact formula (Zuckerman; Bringmann–Ono) and its normalisation precisely.
- Write out the CRT factorisation in Lemma 2 explicitly for even $c$ and for $3\mid c$. It is routine, but it must be done.
