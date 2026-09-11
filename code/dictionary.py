# Normalisation dictionary check for the main cusp type (d = t): the model amplitude of term (k, gamma = mu - e)
# should equal a constant times the Gamma_0(N) Kloosterman sum of F(24 tau) (N = 576 t, real character chi(d) = (t/d)^m):
#   S(n) = sum_{d mod 24k, (d,24k)=1} chi(d) e( ( (24e - m(t-1)) * dbar + (24n - m(t-1)) * d ) / (24k) ).
import mpmath as mp, math
from fractions import Fraction as F
from group_b1 import term_amp
mp.mp.dps=30
def kron(a,n):   # Kronecker symbol (a/n), n odd positive
    a%=n; r=1
    while a:
        while a%2==0:
            a//=2
            if n%8 in (3,5): r=-r
        a,n=n,a
        if a%4==3 and n%4==3: r=-r
        a%=n
    return r if n==1 else 0
def ghn(t,m,d):   # (t/d)^m for d coprime to 6t (d odd): use Jacobi with d as modulus
    return kron(t,d)**m
def S(t,m,k,e,n,sign):
    M=24*k; s=0
    for d in range(M):
        if math.gcd(d,M)!=1: continue
        db=pow(d,-1,M)
        s+=ghn(t,m,d)*mp.expjpi(2*mp.mpf(((24*e-m*(t-1))*db*sign+(24*n-m*(t-1))*d*sign)%M)/M)
    return s
for t,m,k,e in ((5,2,5,0),(5,3,10,0),(7,4,7,0),(7,5,14,1),(6,5,6,0)):
    mu=F(m*(t-1),24); gam=mu-e
    for sign in (1,-1):
        rat=[]
        for n in range(3,40,4):
            a=term_amp(t,m,k,gam,n); b=S(t,m,k,e,n,sign)
            if abs(b)>1e-10 and abs(a)>1e-10: rat.append(a/b)
            elif (abs(a)<1e-10)!=(abs(b)<1e-10): rat.append(None)
        good=rat and None not in rat and all(abs(r-rat[0])<1e-8 for r in rat)
        if good: print('t=%d m=%d k=%d e=%d sign=%+d: model amplitude = const * Gamma_0(N) Kloosterman sum, const = %s'%(t,m,k,e,sign,mp.nstr(rat[0],6))); break
    else: print('t=%d m=%d k=%d e=%d: no match with this normalisation'%(t,m,k,e))
