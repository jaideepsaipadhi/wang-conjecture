# Group-level B1 test: for cross-type equal-growth groups G (members (k_i,gamma_i)), and c' coprime to 6t,
# does "sum_i A_{(k_i,gamma_i)}(n) = 0"  imply  "sum_i A_{(k_i c',gamma_i)}(n) = 0"  (and conversely)?
import mpmath as mp, math
from model import growing_terms
mp.mp.dps=30
def term_amp(t,m,k,gam,n):
    s=0
    for h in range(k):
        if math.gcd(h,k)!=1: continue
        for g,cf in growing_terms(t,m,h,k):
            if g==gam: s+=cf*mp.expjpi(-2*mp.mpf(h*n)/k)
    return s
cases=[(6,8,[(3,'1/6'),(6,'2/3')]),(8,8,[(4,'1/3'),(8,'4/3')]),(10,8,[(5,'1/2'),(10,'2')]),(12,3,[(4,'1/24'),(12,'3/8')]),(6,16,[(3,'1/3'),(6,'4/3')])]
from fractions import Fraction as F
fwd=0; bwd=0; tot=0; ex=[]
for t,m,G in cases:
    G=[(k,F(g)) for k,g in G]
    for cp in (5,7,11,13):
        if math.gcd(cp,6*t)!=1: continue
        mod=math.lcm(*[k*cp for k,_ in G])
        for n in range(0,mod,max(1,mod//60)):
            a=sum(term_amp(t,m,k,g,n) for k,g in G)
            b=sum(term_amp(t,m,k*cp,g,n) for k,g in G)
            za=abs(a)<1e-12; zb=abs(b)<1e-12; tot+=1
            if za and not zb: fwd+=1; ex.append((t,m,cp,n,'G=0, G*c!=0'))
            if zb and not za: bwd+=1; ex.append((t,m,cp,n,'G*c=0, G!=0'))
print('group-level tests: %d;  G vanishes but G*c does not: %d;  G*c vanishes but G does not: %d'%(tot,fwd,bwd)); print(ex[:8])
