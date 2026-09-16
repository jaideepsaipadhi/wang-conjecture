# Galois form of the scaling law: for a prime q outside S,
#   A_{(kq,gamma)}(n) = K_q(n;beta) * sigma_g( A_{(k,gamma)}(n) ),   g = q^{-1} mod M,
# where sigma_g is the Galois automorphism zeta_M -> zeta_M^g applied to the exact cyclotomic value of the amplitude.
# Implementation: recompute the model amplitude with every root of unity e(x) replaced by e(g x) (positive real
# constants (t/d)^{m/2} untouched).
import mpmath as mp, math
from fractions import Fraction as F
import model
mp.mp.dps=30
from sympy import factorint
from sympy.ntheory import jacobi_symbol
def kron(a,n):
    # Kronecker symbol (a/n), n>0
    if n==1: return 1
    r=1
    while n%2==0:
        n//=2
        if a%2==0: return 0
        r*= 1 if a%8 in (1,7) else -1
    return r*(jacobi_symbol(a%n,n) if n>1 else 1)
def gsign(D,m,g):
    # Galois action of sigma_g on sqrt(D)^m (real positive):  sqrt(D0) -> kron(disc,g) sqrt(D0)
    if m%2==0 or g==1: return 1
    f=factorint(D); D0=1
    for p,e in f.items():
        if e%2: D0*=p
    if D0==1: return 1
    disc=D0 if D0%4==1 else 4*D0
    return kron(disc,g)
def amp(t,m,k,gam,n,g=1):
    ej=mp.expjpi
    s=0
    for h in range(k):
        if math.gcd(h,k)!=1: continue
        # replicate growing_terms with Galois twist g
        d=math.gcd(t,k); g0=F(m*(d*d-t),24*t)
        if g0<=0: return 0
        hp=(-pow(h,-1,k))%k; th=(t*h)//d; kd=k//d
        hpd=(-pow(th,-1,kd))%kd if kd>1 else 0
        zq=ej(2*g*mp.mpf(hp)/k); zx=ej(2*g*mp.mpf(d*hpd)/k)
        ser=model.series_prod(m,F(d*d,t),zq,zx,g0)
        sd=-model.ded(h,k)+model.ded(th,kd)
        phase=mp.mpf(t/d)**(mp.mpf(m)/2)*ej(g*m*mp.mpf(sd.numerator)/sd.denominator)*gsign(t//d,m,g)
        for e,c in ser.items():
            if g0-e==gam: s+=phase*c*ej(-2*g*mp.mpf(h*n)/k)
    return s
cases=[(6,8,[(3,F(1,6)),(6,F(2,3))]),(12,3,[(4,F(1,24)),(12,F(3,8))]),(10,8,[(5,F(1,2)),(10,F(2))]),(8,8,[(4,F(1,3)),(8,F(4,3))]),(12,6,[(8,F(1,12)),(24,F(3,4))]),(9,4,[(9,F(1))]),(15,2,[(15,F(7,6))])]
tot=0; bad=0; zero=0
for t,m,G in cases:
    Mbig=24*24*math.lcm(*[k for k,_ in G])*t
    for q in (5,7,11,13):
        if math.gcd(q,6*t)!=1: continue
        g=pow(q,-1,Mbig)
        for n in range(0,40,3):
            rs=[]; ok=True
            for k,gam in G:
                a=amp(t,m,k,gam,n,g); b=amp(t,m,k*q,gam,n)
                if abs(a)>1e-12: rs.append(b/a)
                elif abs(b)>1e-12: ok=False
            if rs and any(abs(r-rs[0])>1e-8 for r in rs): ok=False
            if rs and abs(rs[0])<1e-12: zero+=1
            tot+=1; bad+=(not ok)
print('Galois scaling law: %d checks, %d violations, %d with zero factor'%(tot,bad,zero))

# which cases fail?
fails={}
for t,m,G in cases:
    Mbig=24*24*math.lcm(*[k for k,_ in G])*t
    for q in (5,7,11,13):
        if math.gcd(q,6*t)!=1: continue
        g=pow(q,-1,Mbig)
        for n in range(0,40,3):
            rs=[]
            for k,gam in G:
                a=amp(t,m,k,gam,n,g); b=amp(t,m,k*q,gam,n)
                if abs(a)>1e-12: rs.append(b/a)
            if rs and any(abs(r-rs[0])>1e-8 for r in rs):
                fails.setdefault((t,m,q),[]).append([mp.nstr(r,4) for r in rs])
for key,v in fails.items(): print(key, v[:2])
