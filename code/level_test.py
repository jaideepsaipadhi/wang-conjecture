# Claim to test (prime p>=5, odd m): on a degenerate class r (alpha(r)=0), for every term order e,
# the amplitude at EVERY cusp k = p*j vanishes on the class iff the level-p amplitude (k=p) vanishes, i.e.
# vanishing is decided by e alone ("bad e") -- not by k.  Then either (C) holds or the class is identically 0.
import mpmath as mp, math
from model import growing_terms
mp.mp.dps=30
def amps(p,m,k):
    d={}
    for h in range(k):
        if math.gcd(h,k)!=1: continue
        for gam,cf in growing_terms(p,m,h,k):
            d.setdefault(gam,[]).append((h,cf))
    return d
def van(lst,k,p,r):
    for j in range(k):
        n=r+p*j
        if abs(sum(cf*mp.expjpi(-2*mp.mpf(h*n)/k) for h,cf in lst))>1e-12: return False
    return True
viol=[]; checks=0
for p,m in ((5,7),(5,9),(7,5),(7,9),(11,5),(11,7),(13,5),(5,13),(7,11)):
    A={j:amps(p,m,p*j) for j in range(1,7)}
    base=A[1]
    degen=[r for r in range(p) if van(base.get(max(base),[]),p,p,r)]      # gamma max  <-> e=0
    for r in degen:
        for gam,lst in base.items():
            v1=van(lst,p,p,r)
            for j in range(2,7):
                lj=A[j].get(gam,[])
                if not lj: continue
                checks+=1
                vj=van(lj,p*j,p,r)
                if v1 and not vj: viol.append((p,m,r,str(gam),p*j))
print('checks:',checks,' cases where level-p amplitude vanishes on a degenerate class but some k=pj amplitude does not:',viol[:10])
