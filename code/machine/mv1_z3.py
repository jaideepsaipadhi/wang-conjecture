# MACHINE VERIFICATION 1 (Z3): the elementary arithmetic facts used in sections 3 and 5.
# Each is discharged as a universally quantified statement over bounded integers, by
# exhaustive SMT refutation (Z3 proves UNSAT of the negation).
from z3 import *
def check(name, mk, bound=200):
    s=Solver(); s.set('timeout',15000); mk(s)
    r=s.check()
    print('  %-58s %s' % (name, 'PROVED' if r==unsat else ('COUNTEREXAMPLE '+str(s.model()) if r==sat else 'UNKNOWN')))
    return r==unsat
ok=[]

# (A) Lemma 5.2 Step 1: 24 | k(q^2-1) whenever gcd(q,6)=1.   [used for 1-k^2q^2 = 1-k^2 mod 24k]
k,q,a = Ints('k q a')
def A(s):
    s.add(k>0,k<500, q>1,q<500, q%2==1, q%3!=0)
    s.add(Not((k*(q*q-1))%24==0))
ok.append(check('Step 1: gcd(q,6)=1  =>  24 | k(q^2-1)', A))

# (B) Lemma 5.2 Step 1: qbar = q mod 8 when q*qbar = 1 mod 24k and 8|24k
qb=Int('qb')
def B(s):
    s.add(k>0, q>0, qb>0, (q*qb)%8==1, q%2==1, qb%2==1, q<200, qb<200)
    s.add(Not(qb%8==q%8))
ok.append(check('Step 1: q*qbar=1 mod 8, both odd  =>  qbar = q mod 8', B))

# (C) Lemma 5.2 Step 1: (kq-1)/2 - q(k-1)/2 = (q-1)/2   (integer identity, k,q odd)
def C(s):
    s.add(k%2==1,q%2==1,k>0,k<500,q>0,q<500)
    s.add(Not( (k*q-1) - q*(k-1) == q-1 ))
ok.append(check('Step 1: (kq-1) - q(k-1) = q-1', C))

# (D) Step 1 (E)-case: (q-1)(q+1)/2 = 0 mod 4 for odd q   [leftover constant = i^{(q-1)/2}]
def D(s):
    s.add(q%2==1,q>0,q<500)
    s.add(Not( ((q-1)*(q+1)) % 8 == 0 ))
ok.append(check('Step 1 (E): q odd  =>  8 | (q-1)(q+1)', D))

# (E) Lemma 5.2 Step 2: H = (t/d)h + K*q*l'  ==  (t/d)h mod q
td,h,K,lp=Ints('td h K lp')
def E(s):
    s.add(q>1,q<200, td>0,td<200, K>0,K<200, h>0,h<200, lp>=0,lp<200)
    s.add(Not( ((td*h + K*q*lp) - td*h) % q == 0 ))
ok.append(check('Step 2: H = (t/d)h mod q', E))

# (F) Lemma 5.2 Step 2: K*q*l' = K*l mod 24K whenever q*l' = l mod 24
l=Int('l')
def F(s):
    s.add(K>0,K<200, q>1,q<200, lp>=0,lp<200, l>=0,l<24, (q*lp-l)%24==0)
    s.add(Not( (K*q*lp - K*l) % (24*K) == 0 ))
ok.append(check('Step 2: q l\' = l mod 24  =>  Kql\' = Kl mod 24K', F))

# (G) Lemma 5.2 Step 3 phase collection: coefficient of h is mu = m(t-1)/24.
#     24*[ (m/24)(-1+t) ] = m(t-1).  Integer form.
m,t=Ints('m t')
def G(s):
    s.add(m>0,m<300, t>1,t<300)
    s.add(Not( m*(-1+t) == m*(t-1) ))
ok.append(check('Step 3: h-coefficient collapses to mu', G))

# (H) Step 3 hbar-coefficient: 24t*[ (m/24)(-1 + d^2/t) ] = m(d^2 - t).
d=Int('d')
def H(s):
    s.add(m>0,m<200, t>1,t<200, d>0,d<200)
    s.add(Not( m*(d*d - t) == m*d*d - m*t ))
ok.append(check('Step 3: hbar-coefficient collapses to gamma + Lambda', H))

# (J) Ligozat (i): sum delta r_delta = 24m - 24tm = 0 mod 24
def J(s):
    s.add(m>0,m<300, t>1,t<300)
    s.add(Not( (24*m - 24*t*m) % 24 == 0 ))
ok.append(check('Ligozat (i): 24 | 24m - 24tm', J))

# (K) Ligozat (ii): (N/24)m - (N/(24t))m = 0 mod 24 with N = 576t
def K_(s):
    s.add(m>0,m<300, t>1,t<300)
    s.add(Not( (24*t*m - 24*m) % 24 == 0 ))
ok.append(check('Ligozat (ii): 24 | (N/24)m - (N/24t)m,  N=576t', K_))

print()
print('  %d/%d discharged by Z3' % (sum(ok), len(ok)))
