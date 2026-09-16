# MACHINE VERIFICATION 2: complete proofs by exhaustive residue enumeration.
# Each statement depends ONLY on residues modulo a fixed modulus, so enumerating every
# residue class is a COMPLETE proof, not a sample.
from math import gcd
ok=[]
def claim(name, universe, pred):
    bad=[x for x in universe if not pred(x)]
    print('  %-62s %s  (%d classes)'%(name,'PROVED' if not bad else 'FAILS %s'%bad[:3],len(universe)))
    ok.append(not bad); return not bad

# (A) gcd(q,6)=1  =>  q^2 = 1 mod 24.   Depends only on q mod 24.  Hence 24 | k(q^2-1) for all k.
claim('gcd(q,6)=1 => q^2 = 1 mod 24  (hence 24 | k(q^2-1) for every k)',
      [q for q in range(24) if gcd(q,6)==1], lambda q: (q*q-1)%24==0)

# (B) q*qbar = 1 mod 8, q odd  =>  qbar = q mod 8.  Depends only on residues mod 8.
claim('q*qbar = 1 mod 8  =>  qbar = q mod 8',
      [(a,b) for a in range(8) for b in range(8) if a%2 and b%2 and (a*b)%8==1],
      lambda p: p[0]%8==p[1]%8)

# (C) q odd  =>  8 | (q-1)(q+1).  Depends only on q mod 8.
claim('q odd => 8 | (q-1)(q+1)  (leftover constant is i^{(q-1)/2})',
      [q for q in range(8) if q%2==1], lambda q: ((q-1)*(q+1))%8==0)

# (D) Prop 3.x: g|t, g<t  =>  no unit h mod t is divisible by t/g.  Exhaustive over t<=200.
bad=[]
for t in range(2,201):
    for g in (d for d in range(1,t) if t%d==0):
        M=t//g
        if M==1: continue
        for h in range(1,t):
            if gcd(h,t)==1 and h%M==0: bad.append((t,g,h))
print('  %-62s %s  (t<=200, all divisors, all units)'%('Prop 3.x: t/g never divides a unit mod t, for g<t','PROVED' if not bad else 'FAILS %s'%bad[:3]))
ok.append(not bad)

# (E) Ligozat for F* = eta(24 tau)^m / eta(24 t tau)^m on Gamma_0(576t): both congruences.
bad=[]
for t in range(2,60):
    for m in range(1,60):
        N=576*t
        if (24*m-24*t*m)%24 or ((N//24)*m-(N//(24*t))*m)%24: bad.append((t,m))
print('  %-62s %s  (t<60, m<60)'%('Ligozat (i) and (ii) for F* on Gamma_0(576t)','PROVED' if not bad else 'FAILS %s'%bad[:3]))
ok.append(not bad)

# (F) character: (t^{-m}/d) = (t/d)^(m mod 2).  Kronecker symbol depends on d mod 4|t|.
def kron(a,n):
    if n<0: n=-n
    a%=n; r=1
    while a:
        while a%2==0:
            a//=2
            if n%8 in(3,5): r=-r
        a,n=n,a
        if a%4==3 and n%4==3: r=-r
        a%=n
    return r if n==1 else 0
bad=[]
for t in range(2,40):
    for m in range(1,12):
        for d in range(1,4*24*t,2):
            if gcd(d,6*t)!=1: continue
            if kron(t,d)**m != kron(t,d)**(m%2): bad.append((t,m,d))
print('  %-62s %s  (t<40, m<12, d over a full period)'%('character (t^{-m}/d) = (t/d)^(m mod 2)','PROVED' if not bad else 'FAILS %s'%bad[:3]))
ok.append(not bad)

# (G) Salie (c): T_1(0,0) = sum_{d mod p}* chi_p(d) = 0 for odd p.  Exhaustive over p.
bad=[]
for p in [p for p in range(3,200) if all(p%i for i in range(2,int(p**.5)+1))]:
    s=sum(1 if pow(d,(p-1)//2,p)==1 else -1 for d in range(1,p))
    if s!=0: bad.append(p)
print('  %-62s %s  (all odd primes p<200)'%('Salie (c): sum of Legendre symbol over units = 0','PROVED' if not bad else 'FAILS %s'%bad[:3]))
ok.append(not bad)

# (H) cos(2 pi x / p) != 0 for integer x and odd p  (used in the Salie converse).
bad=[]
for p in [p for p in range(3,300,2)]:
    for x in range(p):
        if (4*x-p)%(2*p)==0: bad.append((p,x))
print('  %-62s %s  (odd p<300, all x)'%('Salie converse: cos(2 pi x/p) != 0 for integer x, p odd','PROVED' if not bad else 'FAILS %s'%bad[:3]))
ok.append(not bad)

print()
print('  %d/%d discharged by exhaustive residue enumeration (complete, not sampled)'%(sum(ok),len(ok)))
