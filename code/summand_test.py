# Summand-level check of the Galois/CRT structure: for h coprime to 24kq,
#   s_{kq}(h) = eps(k,q) * sigma_{qbar}( s_k(h) ) * e( kbar*((mu-n) h + gamma hbar)/q ) * ((t/d)|q)^m
# with eps independent of h and n.  (s_k(h) = model summand of term (k,gamma) at h, taken with h mod k.)
import mpmath as mp, math, sys
from fractions import Fraction as F
import model
sys.argv=['x']
exec(open('galois_scaling.py').read().split('cases=[')[0])
mp.mp.dps=30
def summand(t,m,k,gam,n,h,g=1):
    ej=mp.expjpi; d=math.gcd(t,k); g0=F(m*(d*d-t),24*t)
    hh=h%k; hp=(-pow(hh,-1,k))%k; th=(t*hh)//d; kd=k//d
    hpd=(-pow(th,-1,kd))%kd if kd>1 else 0
    zq=ej(2*g*mp.mpf(hp)/k); zx=ej(2*g*mp.mpf(d*hpd)/k)
    ser=model.series_prod(m,F(d*d,t),zq,zx,g0)
    sd=-model.ded(hh,k)+model.ded(th,kd)
    ph=mp.mpf(t/d)**(mp.mpf(m)/2)*ej(g*m*mp.mpf(sd.numerator)/sd.denominator)*gsign(t//d,m,g)
    return sum(ph*c*ej(-2*g*mp.mpf(hh*n)/k) for e,c in ser.items() if g0-e==gam)
bad=0; tot=0
for t,m,k,gam in ((6,8,3,F(1,6)),(6,8,6,F(2,3)),(12,3,4,F(1,24)),(12,3,12,F(3,8)),(10,8,10,F(2)),(8,8,8,F(4,3)),(12,6,8,F(1,12))):
    d=math.gcd(t,k); mu=F(m*(t-1),24)
    for q in (5,7,11):
        if math.gcd(q,6*t*k)!=1: continue
        M=24*k*q; g=pow(q,-1,24*24*k*t*8)
        kb=pow(k,-1,q)
        ratios=[]
        for n in (1,4):
            for h in range(1,M):
                if math.gcd(h,M)!=1: continue
                if len(ratios)>60: break
                a=summand(t,m,k*q,gam,n,h); b=summand(t,m,k,gam,n,h,g)
                # q-part exponent: kbar*((mu-n) h + gamma hbar)/q ; mu,gamma rational with denominators prime to q
                muq=(mu.numerator*pow(mu.denominator,-1,q))%q; gq=(gam.numerator*pow(gam.denominator,-1,q))%q
                qp=mp.expjpi(2*mp.mpf((kb*((muq-n)*h+gq*pow(h,-1,q)))%q)/q)*(kron(t//d,q)**m)
                if abs(b)>1e-15: ratios.append(a/(b*qp))
        tot+=1
        if max(abs(r-ratios[0]) for r in ratios)>1e-10: bad+=1; print('FAIL',t,m,k,q,[mp.nstr(r,4) for r in ratios[:6]])
print('summand-level structure: %d (term,q) cases, %d failures'%(tot,bad))
