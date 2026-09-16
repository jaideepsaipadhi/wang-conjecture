# Compute eps(k,q) (summand-level constant of Lemma 5.2) for many terms, to find what it depends on,
# and check equality across equal-growth group members.
import mpmath as mp, math, sys, itertools
from fractions import Fraction as F
sys.argv=['x']
exec(open('summand_test.py').read().split('bad=0; tot=0')[0])
def eps(t,m,k,gam,q):
    d=math.gcd(t,k); mu=F(m*(t-1),24); M=24*k*q; g=pow(q,-1,24*24*k*t*8); kb=pow(k,-1,q); n=1
    muq=(mu.numerator*pow(mu.denominator,-1,q))%q; gq=(gam.numerator*pow(gam.denominator,-1,q))%q
    for h in range(1,M):
        if math.gcd(h,M)!=1: continue
        b=summand(t,m,k,gam,n,h,g)
        if abs(b)<1e-15: continue
        a=summand(t,m,k*q,gam,n,h)
        qp=mp.expjpi(2*mp.mpf((kb*((muq-n)*h+gq*pow(h,-1,q)))%q)/q)*(kron(t//d,q)**m)
        return a/(b*qp)
sys.argv=["x","4","1"]
from coincide import poles
rows=[]
for t,m in ((6,8),(12,3),(12,6),(8,8),(10,8),(9,4),(15,2),(12,5),(20,3)):
    ds=[d for d in range(1,t+1) if t%d==0 and d*d>t]
    terms=[(k,g) for k in range(1,2*t+1) if math.gcd(k,t) in ds for g in poles(t,m,math.gcd(k,t))]
    byb={}
    for k,g in terms: byb.setdefault(g/(k*k),[]).append((k,g))
    for grp in byb.values():
        if len(grp)<2: continue
        for q in (5,7,11):
            if math.gcd(q,6*t)!=1: continue
            es=[eps(t,m,k,g,q) for k,g in grp]
            same=all(abs(e-es[0])<1e-9 for e in es)
            rows.append((t,m,q,[(k,str(g)) for k,g in grp],[mp.nstr(e,4) for e in es],same))
print('multi-member groups tested:',len(rows),' eps equal across members in all:',all(r[-1] for r in rows))
for r in rows[:6]: print(r)

ok=True; cnt=0
for t,m,q,grp,es,same in rows:
    vals=[]
    for (k,g),e in zip(grp,es):
        d=math.gcd(t,k); vals.append(complex(mp.mpmathify(e.replace('(','').replace(')','').replace(' ','')))*kron(t//d,q)**m)
    cnt+=1
    if max(abs(v-vals[0]) for v in vals)>1e-3: ok=False; print('differs:',t,m,q,grp,vals)
print('groups:',cnt,' eps*((t/d)|q)^m common across members in all:',ok)

allone=True
for t,m,q,grp,es,same in rows:
    for (k,g),e in zip(grp,es):
        d=math.gcd(t,k); v=complex(mp.mpmathify(e.replace('(','').replace(')','').replace(' ','')))*kron(t//d,q)**m
        if abs(v-1)>1e-3: allone=False; print('not 1:',t,m,q,k,v)
print('eps * ((t/d)|q)^m == 1 for every member of every tested group:',allone)
