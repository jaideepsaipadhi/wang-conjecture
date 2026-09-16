# Wang's Conjecture 1.6 for composite t: proof architecture

**Setting.** Let $F=\eta(\tau)^m/\eta(t\tau)^m$ on $\Gamma_0(t)$. Its growing terms sit at the cusps $h/k$ with $\gcd(k,t)^2>t$. Each term is labelled by $(k,\gamma)$; on a residue class $r\bmod t$ its amplitude is a generalized Kloosterman sum of the eta-quotient multiplier.

**Step A (least period divisible by t) — proved for every t.**
The leading term sits at $k=t$. Its amplitude
$$\alpha(r)=\sum_{h\bmod t}{}^{*}\psi(h)\,e(-hr/t),\qquad |\psi(h)|=1,$$
has its Fourier support on the units mod $t$.

Let $g$ be a proper divisor of $t$. Averaging $e(-hr/t)$ over a coset $r+g\mathbb Z/t\mathbb Z$ gives $0$ unless $(t/g)\mid h$, and that is impossible for a unit $h$. So $\alpha$ sums to $0$ on every coset mod $g$. Since $\alpha\not\equiv0$, some coset contains a nonzero value. That coset then contains both signs, and those classes are non-degenerate, so their eventual sign is $\operatorname{sgn}\alpha$.

Now suppose the eventual sign sequence has a period $L$ with $t\nmid L$. Put $g=\gcd(L,t)$. Shifting by $L$ acts on the non-degenerate classes as the subgroup $g\mathbb Z/t\mathbb Z$, so the sign would be constant on each coset mod $g$. That contradicts the previous paragraph. Hence $t\mid L$. This supersedes the prime-only argument.

**Step B (eventual periodicity) — reduced to two local lemmas; both verified numerically.**
- **B1 (coprime propagation).** If the term $(k,\gamma)$ vanishes on the class $r$, then so does $(kc',\gamma)$ for every $c'$ with $\gcd(c',t)=1$. This follows from CRT factorization of the Kloosterman sum. The twisting by $\overline{c'}$ changes local arguments only by squares and units coprime to the local modulus. Verified: 825 checks, $t\in\{10,12,14,15\}$, no violations.
- **B2 (exponent stabilization).** For $p\mid t$, whether the $p$-local factor vanishes depends on $v_p(k)$ only through $\min(v_p(k),A_p)$, for some bound $A_p$. This is $p$-adic stationary phase: for large exponent the sum localizes on the solutions of a quadratic congruence, and their solvability stabilizes by Hensel's lemma. Verified: $k=t\cdot p^j$ for $j\le3$ with $t\in\{6,10,12\}$, $p\in\{2,3\}$. The vanishing pattern is constant for $j\ge1$ in every case.

**Consequence.** Whether a term vanishes on a class depends only on a finite amount of data, the "shape" of $k$. Within a shape, the smallest $k$ has the largest growth rate. So on each class the sign is decided by one term, taken from a finite list of $k$ whose exponents are bounded, or else every term vanishes and $c(n)=0$ by the exact formula. Either way the sign is eventually periodic. Combined with Step A, this gives Wang's Conjecture 1.6 for every $t$.

**What remains to write rigorously.**
1. The eta-quotient multiplier on $\Gamma_0(t)$ at every cusp type, and the CRT factorization of its Kloosterman sums, including the primes 2 and 3 (where 24 is not invertible) and cusps with $\gcd(k,t)\neq t$. It is routine, but long.
2. B2 in general: the $p$-adic stationary-phase evaluation of the local sums, with the explicit bound $A_p$.
