# Stable-range Galois scaling for q in S: A_{(kq,gamma)}(n) = K^st(n;beta,alpha) * sigma_g(A_{(k,gamma)}(n)),
# g == q^{-1} on the prime-to-q part of the modulus, g == 1 on the q-part.  Test: ratio common across group members.
import mpmath as mp, math, sys
from fractions import Fraction as F
sys.argv=['x']
exec(open('galois_scaling.py').read().split('cases=[')[0])
mp.mp.dps=30
def crt_g(q,Mbig):
    qa=1
    while Mbig%(qa*q)==0: qa*=q
    rest=Mbig//qa
    # g = qbar mod rest, 1 mod qa
    g1=pow(q,-1,rest)
    return (g1*qa*pow(qa,-1,rest)+1*rest*pow(rest,-1,qa))%Mbig
cases=[(12,3,[(4,F(1,24)),(12,F(3,8))],2,3),(12,6,[(8,F(1,12)),(24,F(3,4))],2,2),(6,8,[(3,F(1,6)),(6,F(2,3))],3,1),(10,8,[(5,F(1,2)),(10,F(2))],5,1)]
tot=0; bad=0
for t,m,G,q,j0 in cases:
    for j in (j0,j0+1):
        Gs=[(k*q**j,gam) for k,gam in G]
        Mbig=24*24*t*8*math.lcm(*[k for k,_ in Gs])*q
        g=crt_g(q,Mbig)
        for n in range(0,40,3):
            rs=[]; ok=True
            for k,gam in Gs:
                a=amp(t,m,k,gam,n,g); b=amp(t,m,k*q,gam,n)
                if abs(a)>1e-12: rs.append(b/a)
                elif abs(b)>1e-12: ok=False
            if rs and any(abs(r-rs[0])>1e-8 for r in rs): ok=False
            tot+=1; bad+=(not ok)
print('stable-range Galois scaling (q in S): %d checks, %d violations'%(tot,bad))
