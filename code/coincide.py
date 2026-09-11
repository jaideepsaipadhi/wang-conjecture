# Enumerate equal-growth-rate coincidences between terms of DIFFERENT cusp types for G_t^m.
# Pole orders at cusp type d: gamma = g0(d) - e, e in the exponent semigroup N + (d^2/t) N, g0(d) = m(d^2-t)/(24t).
# Terms (k,gamma),(k',gamma') coincide iff sqrt(gamma)/k = sqrt(gamma')/k' with gcd(k,t)=d, gcd(k',t)=d'.
from fractions import Fraction as F
import math, sys
def poles(t,m,d):
    g0=F(m*(d*d-t),24*t); a=F(d*d,t); out=set()
    i=0
    while i<g0:
        j=0
        while i+j*a<g0:
            out.add(g0-i-j*a); j+=1
        i+=1
    return out
def is_sq(q):
    q=F(q); n,dn=q.numerator,q.denominator
    return math.isqrt(n)**2==n and math.isqrt(dn)**2==dn
TMAX,MMAX=int(sys.argv[1]),int(sys.argv[2]); hits=[]
for t in range(4,TMAX+1):
    ds=[d for d in range(1,t+1) if t%d==0 and d*d>t]
    if len(ds)<2: continue
    for m in range(1,MMAX+1):
        P={d:poles(t,m,d) for d in ds}
        found=False
        for i,d in enumerate(ds):
            for d2 in ds[i+1:]:
                for g in P[d]:
                    for g2 in P[d2]:
                        r=g2/g
                        if is_sq(r):
                            s=F(math.isqrt(r.numerator),math.isqrt(r.denominator))   # k2/k = s
                            # need k with gcd(k,t)=d and k2=s*k integer with gcd(k2,t)=d2
                            for k in range(d,40*t,d):
                                if math.gcd(k,t)!=d: continue
                                k2=s*k
                                if k2.denominator==1 and math.gcd(int(k2),t)==d2:
                                    hits.append((t,m,d,str(g),k,d2,str(g2),int(k2))); found=True; break
print('pairs (t,m) with t<=%d, m<=%d having a cross-type equal-growth coincidence: %d'%(TMAX,MMAX,len(set((h[0],h[1]) for h in hits))))
print('examples:',hits[:10])
