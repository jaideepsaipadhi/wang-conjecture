# Stable-range scaling at q in S when group members have DIFFERENT q-exponents
# (e.g. t=12, m=3, group {(4,1/24),(12,3/8)} at q=3: 3-exponents differ by 1).
import mpmath as mp, math, sys
from fractions import Fraction as F
sys.argv=['x']
exec(open('galois_stable.py').read().split('\ncases=[')[0])
mp.mp.dps=30
cases=[(12,3,[(4,F(1,24)),(12,F(3,8))],3,(1,2)),(6,8,[(3,F(1,6)),(6,F(2,3))],2,(2,3)),(12,6,[(8,F(1,12)),(24,F(3,4))],3,(1,2))]
tot=0; bad=0; rep=[]
for t,m,G,q,js in cases:
    for j in js:
        Gs=[(k*q**j,gam) for k,gam in G]
        if max(k for k,_ in Gs)*q>400: continue
        Mbig=24*24*t*8*math.lcm(*[k for k,_ in Gs])*q
        g=crt_g(q,Mbig)
        for n in range(0,40,3):
            rs=[]; ok=True
            for k,gam in Gs:
                a=amp(t,m,k,gam,n,g); b=amp(t,m,k*q,gam,n)
                if abs(a)>1e-12: rs.append(b/a)
                elif abs(b)>1e-12: ok=False
            if rs and any(abs(r-rs[0])>1e-8 for r in rs): ok=False; rep.append((t,m,q,j,n,[mp.nstr(r,4) for r in rs]))
            # group-level implication: G_{beta/q}(n)!=0 => G_beta(n)!=0
            Ga=sum(amp(t,m,k,gam,n) for k,gam in Gs); Gb=sum(amp(t,m,k*q,gam,n) for k,gam in Gs)
            if abs(Gb)>1e-10 and abs(Ga)<1e-10: ok=False; rep.append(('IMPLICATION FAILS',t,m,q,j,n))
            tot+=1; bad+=(not ok)
print('spread-exponent stable tests: %d, failures %d'%(tot,bad)); print(rep[:6])
