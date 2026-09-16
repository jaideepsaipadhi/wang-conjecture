# (1) Explicit Dedekind-sum phase (from Apostol/Knopp eta-multiplier), for gcd(h,24k)=1, hb = inverse of h mod 24k:
#   k odd : exp(pi i s(h,k)) = (h/k) i^{(k-1)/2} e((1-k^2)(h+hb)/(24k))
#   k even: exp(pi i s(h,k)) = (k/h) e((hb(1-k^2) + h(2k^2-3k+1))/(24k))
# (2) Does the term phase (w_{h,k}^{-1} w_{th/d,k/d})^m have the form  (real character in h) * e((a h + b hb)/(24k))  ?
import mpmath as mp, math, itertools
from model import ded
import sys
sys.argv=['x']
exec(open('lemma_checks.py').read().split('# (1) multiplier formula check')[0])
mp.mp.dps=30
def kron(a,n):     # Jacobi symbol (a/n) for odd n>0
    return jacobi(a,n)
def omega_formula(h,k):
    hb=pow(h,-1,24*k)
    if k%2: return kron(h,k)*(1j)**((k-1)//2)*mp.expjpi(2*mp.mpf((1-k*k)*(h+hb))/(24*k))
    # (k/h) for odd h>0 : Jacobi symbol with h as modulus
    return kron(k,h)*mp.expjpi(2*mp.mpf(hb*(1-k*k)+h*(2*k*k-3*k+1))/(24*k))
bad=0; tot=0
for k in range(1,60):
    for h in range(1,24*k):
        if math.gcd(h,24*k)!=1: continue
        if h>200: break
        sd=ded(h%k,k) if k>1 else 0
        w=mp.expjpi(mp.mpf(sd.numerator)/sd.denominator) if k>1 else mp.mpf(1)
        tot+=1; bad+=abs(w-omega_formula(h,k))>1e-20
print('(1) explicit omega formula: %d (h,k) pairs, %d failures'%(tot,bad))
