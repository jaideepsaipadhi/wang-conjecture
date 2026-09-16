# Test the candidate Key Lemma: for each cusp type d (gcd(k,t)=d, d^2>t) and each term gamma, if the amplitude at the
# minimal cusp k=d vanishes for all n = r (mod t), do the amplitudes at every k with gcd(k,t)=d also vanish on n = r (mod t)?
import mpmath as mp, math
from model import growing_terms
mp.mp.dps=30
def amp_by_gamma(t,m,k):
    terms={}
    for h in range(k):
        if math.gcd(h,k)!=1: continue
        for gam,cf in growing_terms(t,m,h,k):
            terms.setdefault(gam,[]).append((h,cf))
    return terms
def vanishes_on_class(terms_list,k,t,r):
    # amplitude A(n)=sum cf e(-h n/k) for all n = r mod t  (check n = r + t*j, j over one period k)
    for j in range(k):
        n=r+t*j
        if abs(sum(cf*mp.expjpi(-2*mp.mpf(h*n)/k) for h,cf in terms_list))>1e-15: return False
    return True
bad=[]; tested=0
for t,m in ((4,12),(6,8),(8,4),(8,6),(9,6),(10,3),(12,3),(12,4),(9,4),(6,5)):
    for d in [d for d in range(1,t+1) if t%d==0 and d*d>t]:
        base=amp_by_gamma(t,m,d)
        ks=[k for k in range(d,4*t+1) if math.gcd(k,t)==d]
        others={k:amp_by_gamma(t,m,k) for k in ks}
        for r in range(t):
            van0=all(vanishes_on_class(v,d,t,r) for v in base.values())
            if not van0: continue
            for k in ks:
                for gam,v in others[k].items():
                    tested+=1
                    if not vanishes_on_class(v,k,t,r): bad.append((t,m,d,r,k,str(gam)))
print('classes/terms tested:',tested,' counterexamples to "vanishing at k=d on class r => vanishing at all k with gcd(k,t)=d":',bad[:10])
