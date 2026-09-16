import mpmath as mp, math, sys
from fractions import Fraction as F
SA=list(sys.argv); sys.argv=["x"]
exec(open('galois_stable.py').read().split('\ncases=[')[0])
mp.mp.dps=25
t,m,G,q=6,8,[(3,F(1,6)),(6,F(2,3))],2
for j in range(int(SA[1]),int(SA[2])):
    Gs=[(k*q**j,gam) for k,gam in G]
    Mbig=24*24*t*8*math.lcm(*[k for k,_ in Gs])*q; g=crt_g(q,Mbig)
    fails=0; ratio_bad=0; cnt=0
    for n in range(1,49,12):
        A=[amp(t,m,k,gam,n) for k,gam in Gs]; B=[amp(t,m,k*q,gam,n) for k,gam in Gs]
        Ag=[amp(t,m,k,gam,n,g) for k,gam in Gs]
        rs=[b/a for a,b in zip(Ag,B) if abs(a)>1e-12]
        if rs and any(abs(r-rs[0])>1e-8 for r in rs): ratio_bad+=1
        if abs(sum(B))>1e-10 and abs(sum(A))<1e-10: fails+=1
        cnt+=1
    print('j=%d (2-exponents %d,%d): n tested %d, ratio not common %d, implication failures %d'%(j,j,j+1,cnt,ratio_bad,fails),flush=True)
