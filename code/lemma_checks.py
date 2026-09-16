# Numerical checks of the two analytic ingredients of the Level Lemma.
import mpmath as mp, math, random
mp.mp.dps=40
def eta(t):
    q=mp.expjpi(2*t); return mp.expjpi(t/12)*mp.nprod(lambda n:1-q**n,[1,mp.inf])
def jacobi(a,n):
    a%=n; r=1
    while a:
        while a%2==0:
            a//=2
            if n%8 in (3,5): r=-r
        a,n=n,a
        if a%4==3 and n%4==3: r=-r
        a%=n
    return r if n==1 else 0
def nu_knopp(a,b,c,d):          # eta multiplier (Knopp, Thm 4.2), c>0
    if c%2==1: return jacobi(d,c)*mp.expjpi(mp.mpf((a+d)*c-b*d*(c*c-1)-3*c)/12)
    return jacobi(c,d)*mp.expjpi(mp.mpf((a+d)*c-b*d*(c*c-1)+3*d-3-3*c*d)/12) if d>0 else None
# (1) multiplier formula check
bad=0; tot=0
tau=mp.mpc('0.123','0.87')
for _ in range(60):
    c=random.randint(1,30); d=random.randint(1,40)
    if math.gcd(c,d)!=1: continue
    # solve ad-bc=1
    a=pow(d,-1,c) if c>1 else 1; b=(a*d-1)//c
    nu=nu_knopp(a,b,c,d)
    if nu is None: continue
    lhs=eta((a*tau+b)/(c*tau+d)); rhs=nu*mp.sqrt(c*tau+d)*eta(tau)
    tot+=1; bad+= abs(lhs-rhs)>1e-20
print('(1) eta multiplier formula: %d random matrices, %d failures'%(tot,bad))
# (2) twisted Salie sums mod p^a vanish when p does not divide uv and chi(uv)=-1
viol=0; tot=0
for p in (5,7,11,13):
    for a in (1,2,3):
        M=p**a
        if M>2200: continue
        for _ in range(40):
            u=random.randrange(1,M); v=random.randrange(1,M)
            if u%p==0 or v%p==0: continue
            S=sum(jacobi(dd,p)*mp.expjpi(2*mp.mpf(u*dd+v*pow(dd,-1,M))/M) for dd in range(1,M) if dd%p)
            chi=jacobi(u*v,p); tot+=1
            if chi==-1 and abs(S)>1e-15: viol+=1
            if chi==1 and abs(S)<1e-10: viol+=1
print('(2) twisted Salie sums mod p^a (a<=3): %d cases, vanish exactly when chi(uv)=-1: violations %d'%(tot,viol))
