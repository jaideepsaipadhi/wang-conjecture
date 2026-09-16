# Stress test A: Lemma 5.2 exact summand identity (no constant) on random terms.
import mpmath as mp, math, sys, random
from fractions import Fraction as F
sys.argv=['x']
exec(open('summand_test.py').read().split('bad=0; tot=0')[0])
sys.argv=["x","4","1"]
from coincide import poles
random.seed(7); mp.mp.dps=30
tot=0; bad=0; ex=[]
tried=0
while tot<45 and tried<400:
    tried+=1
    t=random.randint(2,20); m=random.randint(1,12)
    ds=[d for d in range(1,t+1) if t%d==0 and d*d>t]
    k=random.randint(1,3*t)
    d=math.gcd(k,t)
    if d not in ds: continue
    P=sorted(poles(t,m,d))
    if not P or k*t>400: continue
    gam=random.choice(P)
    q=random.choice([p for p in (5,7,11,13) if math.gcd(p,6*t*k)==1])
    mu=F(m*(t-1),24); M=24*k*q; g=pow(q,-1,24*24*k*t*8); kb=pow(k,-1,q)
    muq=(mu.numerator*pow(mu.denominator,-1,q))%q; gq=(gam.numerator*pow(gam.denominator,-1,q))%q
    worst=0; cnt=0
    for n in (random.randint(0,50),random.randint(0,50)):
        for h in random.sample([h for h in range(1,M) if math.gcd(h,M)==1],12):
            a=summand(t,m,k*q,gam,n,h); b=summand(t,m,k,gam,n,h,g)
            qp=mp.expjpi(2*mp.mpf((kb*((muq-n)*h+gq*pow(h,-1,q)))%q)/q)
            worst=max(worst,abs(a-b*qp)); cnt+=1
    tot+=1
    if worst>1e-12: bad+=1; ex.append((t,m,k,str(gam),q,float(worst)))
print('Lemma 5.2 exact identity, random terms: %d terms, %d failures %s'%(tot,bad,ex[:5]))
