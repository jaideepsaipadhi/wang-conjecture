# Stress test D (end-to-end): the proof says the deciding group is always RIGID (moduli built from primes of S with
# bounded exponents). Predict every sign from rigid groups only and compare with brute force on n in [N/3, N].
import mpmath as mp, math, sys
from model import growing_terms, exact
mp.mp.dps=30
def rigid_ks(t,cap_extra=2,KMAX=600):
    ps=[p for p in (2,3,5,7,11,13,17,19) if (6*t)%p==0]
    ks=[1]
    for p in ps:
        ks=[k*p**e for k in ks for e in range(0, (t%p==0 and (lambda v: v)( max(e2 for e2 in range(10) if t%(p**e2)==0)) or 0)+cap_extra+1)]
    return sorted(k for k in set(ks) if k<=KMAX and math.gcd(k,t)**2>t)
def predict(t,m,L):
    groups={}
    for k in rigid_ks(t):
        for h in range(k):
            if math.gcd(h,k)!=1: continue
            for gam,cf in growing_terms(t,m,h,k):
                beta=mp.sqrt(mp.mpf(gam.numerator)/gam.denominator)/k
                groups.setdefault(mp.nstr(beta,18),[]).append((beta,k,h,cf))
    order=sorted(groups.values(),key=lambda g:-g[0][0])
    pred={}
    for r in range(L):
        s=0
        for g in order:
            A=sum(cf*mp.expjpi(-2*mp.mpf(h*r)/k) for beta,k,h,cf in g)
            if abs(A)>1e-12: s=1 if mp.re(A)>0 else -1; break
        pred[r]=s
    return pred
N=int(sys.argv[1]); tot=0; mism={}
for t,m in ((12,3),(12,4),(10,3),(9,4),(8,5),(6,7),(5,7),(7,5),(4,11),(15,2),(5,5),(3,9),(4,4)):
    L=math.lcm(t,*rigid_ks(t))
    if L>5000: L=math.lcm(t,*[k for k in rigid_ks(t) if k<=4*t*t])
    pred=predict(t,m,L); f=exact(t,m,N)
    bad=[n for n in range(N//3,N+1) if ((f[n]>0)-(f[n]<0))!=pred[n%L]]
    tot+=N-N//3+1; mism[(t,m)]=len(bad)
    print('t=%2d m=%2d  rigid moduli %s  L=%d  mismatches %d'%(t,m,rigid_ks(t)[:8],L,len(bad)),flush=True)
print('end-to-end rigid-only prediction: %d coefficient comparisons, total mismatches %d'%(tot,sum(mism.values())))
