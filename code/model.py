# Explicit Rademacher-type model for G_t^m = (q;q)^m/(q^t;q^t)^m (integer m): all growing principal-part terms
# at all cusps h/k (k<=K), using the eta transformation (hh'=-1 mod k convention).
import mpmath as mp, math
from fractions import Fraction as F
mp.mp.dps=50
def ded(h,k):
    t=F(0)
    for j in range(1,k):
        x=F(j,k); y=F(j*h,k)
        fx=x-math.floor(x)-F(1,2) if x.denominator!=1 else 0
        fy=y-math.floor(y)-F(1,2) if y.denominator!=1 else 0
        t+=fx*fy
    return t
def series_prod(m, a_mult, zeta_q, zeta_x, EMAX):
    # expansion of prod_n (1 - zeta_q^n Y^n)^m / prod_n (1 - zeta_x^n Y^{n a})^m  as dict exponent->coeff, exponents < EMAX
    # exponents are rationals (a = d^2/t)
    terms={F(0):mp.mpc(1)}
    def mul(A,B):
        C={}
        for ea,ca in A.items():
            for eb,cb in B.items():
                e=ea+eb
                if e<EMAX: C[e]=C.get(e,0)+ca*cb
        return C
    n=1
    while n<EMAX:
        # (1 - zq^n Y^n)^m
        fac={F(0):mp.mpc(1)}
        for j in range(1,m+1):
            e=F(n*j)
            if e<EMAX: fac[e]=mp.binomial(m,j)*(-zeta_q**n)**j
        terms=mul(terms,fac); n+=1
    n=1
    while n*a_mult<EMAX:
        # 1/(1 - zx^n Y^{na})^m = sum_j binom(m+j-1,j) (zx^n Y^{na})^j
        fac={F(0):mp.mpc(1)}; j=1
        while n*a_mult*j<EMAX:
            fac[n*a_mult*j]=mp.binomial(m+j-1,j)*(zeta_x**n)**j; j+=1
        terms=mul(terms,fac); n+=1
    return terms
def growing_terms(t,m,h,k):
    d=math.gcd(t,k); g0=F(m*(d*d-t),24*t)
    if g0<=0: return []
    hp=(-pow(h,-1,k))%k                          # h h' = -1 mod k
    th=(t*h)//d; kd=k//d
    hpd=(-pow(th,-1,kd))%kd if kd>1 else 0      # (th/d) h'_d = -1 mod k/d
    zq=mp.expjpi(mp.mpf(2*hp)/k); zx=mp.expjpi(mp.mpf(2*d*hpd)/k)
    ser=series_prod(m,F(d*d,t),zq,zx,g0)
    phase=mp.mpf(t/d)**(mp.mpf(m)/2)*mp.expjpi(m*(-ded(h,k)+ded(th,kd)))
    return [(g0-e, phase*c) for e,c in ser.items() if g0-e>0 and abs(c)>1e-30]   # (gamma, coeff): exp(2 pi gamma / z)
def model(t,m,n,K):
    C=2*n-mp.mpf(m*(t-1))/12; tot=0
    for k in range(1,K+1):
        if math.gcd(t,k)**2<=t: continue
        for h in range(k):
            if math.gcd(h,k)!=1: continue
            for gam,cf in growing_terms(t,m,h,k):
                B=2*mp.mpf(gam.numerator)/gam.denominator
                tot+=cf*mp.expjpi(-2*mp.mpf(h*n)/k)*(2*mp.pi/k)*mp.sqrt(B/C)*mp.besseli(1,(2*mp.pi/k)*mp.sqrt(B*C))
    return tot
def exact(t,m,N):
    f=[0]*(N+1); f[0]=1
    for k in range(1,N+1):
        e=m-(m if k%t==0 else 0)
        for _ in range(e):
            for n in range(N,k-1,-1): f[n]-=f[n-k]
    return f
if __name__=='__main__':
    import sys
    for t,m in ((5,10),(6,8),(4,12),(9,5)):
        f=exact(t,m,400)
        for n in (390,391,392,393):
            M=model(t,m,n,4*t)
            print('t=%d m=%d n=%d exact %s  model %s  imag %s'%(t,m,n,mp.nstr(f[n],12),mp.nstr(mp.re(M),12),mp.nstr(mp.im(M),3)))
