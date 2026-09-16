# Stress test B: stable-range behaviour at q in S with n chosen to vary v_q(24n - m(t-1)).
import mpmath as mp, math, sys
from fractions import Fraction as F
sys.argv=['x']
exec(open('galois_stable.py').read().split('\ncases=[')[0])
mp.mp.dps=25
def vq(x,q):
    if x==0: return 99
    v=0
    while x%q==0: x//=q; v+=1
    return v
res=[]
for t,m,G,q,js in ((10,8,[(5,F(1,2)),(10,F(2))],5,(1,2)),(12,3,[(4,F(1,24)),(12,F(3,8))],2,(2,3,4)),(12,6,[(8,F(1,12)),(24,F(3,4))],2,(2,3))):
    c=m*(t-1)
    ns=[]
    for s in range(0,4):
        for n in range(0,400):
            if vq(24*n-c,q)==s+vq(24,q)*0 and len([x for x in ns if vq(24*x-c,q)==s])<2: ns.append(n)
    for j in js:
        Gs=[(k*q**j,gam) for k,gam in G]
        if max(k for k,_ in Gs)*q>1300: continue
        Mbig=24*24*t*8*math.lcm(*[k for k,_ in Gs])*q; g=crt_g(q,Mbig)
        for n in ns:
            A=[amp(t,m,k,gam,n) for k,gam in Gs]; B=[amp(t,m,k*q,gam,n) for k,gam in Gs]; Ag=[amp(t,m,k,gam,n,g) for k,gam in Gs]
            rs=[b/a for a,b in zip(Ag,B) if abs(a)>1e-12]
            common=(not rs) or all(abs(r-rs[0])<1e-8 for r in rs)
            member_ok=all(not(abs(a)<1e-12 and abs(b)>1e-12) for a,b in zip(Ag,B))
            impl=not(abs(sum(B))>1e-10 and abs(sum(A))<1e-10)
            res.append((t,m,q,j,n,vq(24*n-c,q),common,member_ok,impl))
bad=[r for r in res if not (r[6] and r[7] and r[8])]
print('stable-range stress: %d (group,level,n) checks; failures: %d'%(len(res),len(bad)))
for r in bad[:10]: print(r)
