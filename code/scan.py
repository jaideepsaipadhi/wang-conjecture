# Numerical scan of Wang's Conjecture 1.6: for integers t,m, is sgn c_t^(m)(n) eventually periodic
# with least period divisible by t?  Exact integer coefficients of (q;q)^m/(q^t;q^t)^m.
import sys
def coeffs(t,m,N):
    f=[0]*(N+1); f[0]=1
    for k in range(1,N+1):
        e = m - (m if k%t==0 else 0)          # exponent of (1-q^k): m for all k, minus m if t|k
        for _ in range(e):
            for n in range(N,k-1,-1): f[n]-=f[n-k]
    return f
def sgn(x): return (x>0)-(x<0)
def least_period(s):
    L=len(s)
    for p in range(1,L//3+1):
        if all(s[i]==s[i+p] for i in range(L-p)): return p
    return None
N=int(sys.argv[1]); W=int(sys.argv[2])
res=[]
for t in range(2,13):
    for m in range(1,40):
        if m*(t-1)<=24: continue
        if m*(t-1)>72: break
        f=coeffs(t,m,N); s=[sgn(x) for x in f[N-W:N+1]]
        p=least_period(s)
        res.append((t,m,p))
        flag='' if (p is not None and p%t==0) else '   <-- CHECK'
        print('t=%2d m=%2d  m(t-1)=%3d  least sign period in last %d terms: %s%s'%(t,m,m*(t-1),W,p,flag),flush=True)
