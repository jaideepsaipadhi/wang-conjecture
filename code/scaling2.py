# Corrected Scaling Lemma: for a prime q outside S (q coprime to 6t and to all pole-order numerators/denominators),
#   A_{(kq,gamma)}(n) = R_q(n; gamma/k^2) * A_{(k,gamma)}(n * qbar^2),   qbar = q^{-1} mod 24*lcm(k's),
# with R_q common to all members of an equal-growth group and R_q(n) != 0.
import mpmath as mp, math
from fractions import Fraction as F
from group_b1 import term_amp
mp.mp.dps=30
cases=[(6,8,[(3,F(1,6)),(6,F(2,3))]),(8,8,[(4,F(1,3)),(8,F(4,3))]),(10,8,[(5,F(1,2)),(10,F(2))]),(12,3,[(4,F(1,24)),(12,F(3,8))]),(12,6,[(8,F(1,12)),(24,F(3,4))]),(10,16,[(5,F(1)),(10,F(4))])]
tot=0; bad=0; zeroR=0
for t,m,G in cases:
    Lk=24*math.lcm(*[k for k,_ in G])
    for q in (5,7,11,13,17):
        if math.gcd(q,6*t)!=1: continue
        qb=pow(q,-1,Lk)
        for n in range(0,70,3):
            rs=[]; zs=[]
            for k,g in G:
                a=term_amp(t,m,k,g,n*qb*qb); b=term_amp(t,m,k*q,g,n)
                zs.append((abs(a)<1e-12,abs(b)<1e-12))
                if abs(a)>1e-12: rs.append(b/a)
            ok=all(abs(r-rs[0])<1e-8 for r in rs) and all(zb for za,zb in zs if za)
            if rs and abs(rs[0])<1e-12: zeroR+=1
            tot+=1; bad+=(not ok)
print('corrected scaling lemma: %d checks, %d violations, %d cases with R=0'%(tot,bad,zeroR))
