# Stress test C: brute-force eventual sign periodicity, including n with high q-adic valuation of 24n-m(t-1).
import sys, time
from model import exact
def sgn(x): return (x>0)-(x<0)
def vq(x,q):
    if x==0: return 99
    v=0
    while x%q==0: x//=q; v+=1
    return v
N=int(sys.argv[1]); out=[]
for t,m in ((12,3),(12,4),(10,3),(9,4),(8,5),(6,7),(5,7),(7,5),(4,11),(15,2),(16,2),(18,2)):
    f=exact(t,m,N); lo=N//3
    s=[sgn(f[n]) for n in range(lo,N+1)]
    L=next((p for p in range(1,len(s)//4) if all(s[i]==s[i+p] for i in range(len(s)-p))),None)
    primes=[q for q in (2,3,5,7) if t%q==0]
    deep=max(max(vq(24*n-m*(t-1),q) for n in range(lo,N+1)) for q in primes) if primes else 0
    out.append((t,m,L,L is not None and L%t==0,deep))
    print('t=%2d m=%2d  least sign period on n in [%d,%d]: %s  divisible by t: %s  max q-adic depth covered: %d'%(t,m,lo,N,L,L is not None and L%t==0,deep),flush=True)
