# Condition (C) for odd m, prime p>=5 (mu = m(p-1)/24 > 1):  exists e in [1, mu) with c_m(e) != 0 and
#   e == x (mod p)  or  chi(x(x-e)) = -1,   where x = mu mod p = -m/24 mod p, chi = Legendre symbol mod p.
# (C) implies the degenerate classes are decided by a k=p term with constant sign, hence UPS with least period p.
import sys
from sympy import primerange
MMAX=int(sys.argv[1]); PMAX=int(sys.argv[2]); EMAX=int(sys.argv[3])
def etapow(m,N):
    f=[0]*(N+1); f[0]=1
    for k in range(1,N+1):
        for _ in range(m):
            for n in range(N,k-1,-1): f[n]-=f[n-k]
    return f
fails=[]; worst=0; count=0
for m in range(5,MMAX+1,2):
    c=etapow(m,EMAX)
    for p in primerange(5,PMAX+1):
        mu_num=m*(p-1)                       # mu = mu_num/24
        if mu_num<=24: continue
        count+=1
        x=(-m*pow(24,-1,p))%p
        if x==0:                             # degenerate class r=0 only; e=1 works since c(1)=-m
            continue
        found=None
        for e in range(1,min(p,EMAX+1)):
            if 24*e>=mu_num: break           # need e < mu
            if c[e]==0: continue
            if e==x or pow((x*(x-e))%p,(p-1)//2,p)==p-1: found=e; break
        if found is None: fails.append((m,p))
        else: worst=max(worst,found)
print('odd m in [5,%d], primes p<=%d: %d pairs checked; condition (C) fails for %d pairs %s; largest e needed = %d'%(MMAX,PMAX,count,len(fails),fails[:20],worst))
