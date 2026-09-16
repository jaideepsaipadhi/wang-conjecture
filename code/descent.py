# For each (t,m): order all growing terms by beta = sqrt(gamma)/k; for each residue class r mod L, find the first beta-group
# whose amplitude (sum over terms in the group of coeff*e(-h n/k)) is nonzero at n=r. Report deciding k's and predicted period.
import mpmath as mp, math, sys
from model import growing_terms
mp.mp.dps=30
def analyse(t,m,K,L):
    groups={}
    for k in range(1,K+1):
        if math.gcd(t,k)**2<=t: continue
        for h in range(k):
            if math.gcd(h,k)!=1: continue
            for gam,cf in growing_terms(t,m,h,k):
                beta=mp.sqrt(mp.mpf(gam.numerator)/gam.denominator)/k
                key=mp.nstr(beta,20)
                groups.setdefault(key,[]).append((beta,k,h,cf))
    order=sorted(groups.values(),key=lambda g:-g[0][0])
    decide={}; signs={}
    for r in range(L):
        for g in order:
            A=sum(cf*mp.expjpi(-2*mp.mpf(h*r)/k) for beta,k,h,cf in g)
            if abs(A)>1e-15:
                decide[r]=max(k for _,k,_,_ in g); signs[r]=1 if mp.re(A)>0 else -1; break
        else:
            decide[r]=None; signs[r]=0
    s=[signs[r] for r in range(L)]
    per=next(p for p in range(1,L+1) if L%p==0 and all(s[i]==s[(i+p)%L] for i in range(L)))
    return per, sorted(set(v for v in decide.values() if v)), [r for r in range(L) if decide[r] is None]
for t in range(2,11):
    for m in range(1,40):
        if m*(t-1)<=24: continue
        if m*(t-1)>48: break
        L=t*12
        per,ks,undec=analyse(t,m,4*t,L)
        print('t=%d m=%d  predicted sign period %d (divisible by t: %s)  deciding k values %s  undecided residues mod %d: %s'%(t,m,per,per%t==0,ks,L,undec[:6]),flush=True)
