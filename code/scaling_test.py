# Scaling Lemma test: for a prime q and terms (k,gamma) whose cusp type is unchanged by k -> kq,
# A_{(kq,gamma)}(n) = R_q(n; gamma/k^2, v_q(k)) * A_{(k,gamma)}(n), with R independent of the cusp type.
# Test: within each cross-type equal-growth group, the ratios A(kq)/A(k) agree (where defined), and
# a member vanishes at level kq only if it vanishes at level k or all members' ratios vanish.
import mpmath as mp, math
from fractions import Fraction as F
from group_b1 import term_amp
import sys
sys.argv=["x","4","1"]
from coincide import poles, is_sq
mp.mp.dps=30
def groups(t,m,KMAX):
    ds=[d for d in range(1,t+1) if t%d==0 and d*d>t]
    terms=[(k,g) for k in range(1,KMAX+1) if math.gcd(k,t) in ds for g in poles(t,m,math.gcd(k,t))]
    byb={}
    for k,g in terms: byb.setdefault(g/(k*k),[]).append((k,g))
    return [v for v in byb.values() if len({math.gcd(k,t) for k,_ in v})>1]
bad=0; tot=0; ex=[]
for t,m in ((6,8),(8,8),(10,8),(12,3),(6,16),(10,16),(12,6)):
    for G in groups(t,m,2*t):
        for q in (2,3,5,7):
            if any(math.gcd(k*q,t)!=math.gcd(k,t) for k,_ in G): continue
            for n in range(0,60,7):
                rs=[]; 
                for k,g in G:
                    a=term_amp(t,m,k,g,n); b=term_amp(t,m,k*q,g,n)
                    rs.append((a,b))
                nz=[b/a for a,b in rs if abs(a)>1e-12]
                ok=all(abs(r-nz[0])<1e-8 for r in nz) if nz else True
                ok=ok and all(abs(b)<1e-12 for a,b in rs if abs(a)<1e-12)
                tot+=1
                if not ok: bad+=1; ex.append((t,m,q,n,G))
print('scaling-lemma checks: %d, violations: %d'%(tot,bad), ex[:3])
