# Regression test of the decision procedure (section 5): predict the eventual sign of c_t^(m)(n) on every class
# n mod L from the first non-vanishing term (in growth order) among cusps k <= KMAX, and compare with brute force.
import mpmath as mp, math, sys
from model import growing_terms, exact
mp.mp.dps=30
def predict(t,m,KMAX,L):
    groups={}
    for k in range(1,KMAX+1):
        if math.gcd(t,k)**2<=t: continue
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
            if abs(A)>1e-15: s=1 if mp.re(A)>0 else -1; break
        pred[r]=s
    return pred
if __name__=='__main__':
    N=int(sys.argv[1]); W=int(sys.argv[2]); tot=0; bad=[]
    for t in range(2,13):
        for m in range(1,30):
            if m*(t-1)>60: break
            # L=24t/12t is TOO COARSE: the deciding group can sit at a modulus far above 3t,
            # and the prediction is then not constant mod 24t.  Use the rigid moduli, as stress4.py does.
            from stress4 import rigid_ks
            RK=rigid_ks(t)
            L=math.lcm(t,*RK)
            if L>5000: L=math.lcm(t,*[k for k in RK if k<=4*t*t])
            pred=predict(t,m,max(RK),L)
            f=exact(t,m,N)
            for n in range(N-W,N+1):
                s=(f[n]>0)-(f[n]<0)
                tot+=1
                if s!=pred[n%L]: bad.append((t,m,n)); break
    print('regression: (t,m) pairs with m(t-1)<=60, t<=12; %d coefficient comparisons near n=%d; mismatching pairs: %s'%(tot,N,sorted(set((a,b) for a,b,_ in bad))[:20]))
