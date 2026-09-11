# Local sums T(X,Y;p^A) = sum_{h mod p^A, p not| h} psi(h) e((X h + Y hbar)/p^A), psi a real character.
# Check (G) Galois invariance of vanishing:  T(wX,wY)=0 <=> T(X,Y)=0 for units w;
#       (S) stabilisation: for fixed integers X,Y (Y != 0) the vanishing of T is independent of A for A >= 2 v_p(Y)+4 (+3 if p=2).
import cmath, math, random
def chars(p):
    if p==2: return {'1':lambda h:1,'-4':lambda h:(1 if h%4==1 else -1),'8':lambda h:(1 if h%8 in (1,7) else -1),'-8':lambda h:(1 if h%8 in (1,3) else -1)}
    leg=lambda h:(1 if pow(h,(p-1)//2,p)==1 else -1)
    return {'1':lambda h:1,'leg':leg}
def T(psi,X,Y,p,A):
    M=p**A; s=0
    for h in range(1,M):
        if h%p==0: continue
        s+=psi(h)*cmath.exp(2j*math.pi*((X*h+Y*pow(h,-1,M))%M)/M)
    return s
def vp(x,p):
    v=0
    while x%p==0: x//=p; v+=1
    return v
random.seed(1); gv=0; gt=0; sv=0; st=0
for p in (2,3,5,7):
    for name,psi in chars(p).items():
        for _ in range(25):
            A=random.randint(1,5 if p<5 else 3); M=p**A
            X=random.randrange(M); Y=random.randrange(1,M); w=random.choice([u for u in range(1,M) if u%p])
            z1=abs(T(psi,X,Y,p,A))<1e-8; z2=abs(T(psi,w*X%M,w*Y%M,p,A))<1e-8
            gt+=1; gv+=(z1!=z2)
        for _ in range(12):
            b=random.randint(0,1); Y=p**b*random.choice([u for u in range(1,50) if u%p]); X=random.randrange(-300,300)
            A0=2*b+4+(3 if p==2 else 0)
            Amax=A0+2 if p<5 else A0+1
            if p**Amax>5000: Amax=A0
            pats=[abs(T(psi,X,Y,p,A))<1e-8 for A in range(A0,Amax+1)]
            st+=1; sv+=(len(set(pats))>1)
print('(G) Galois invariance of vanishing: %d cases, %d violations'%(gt,gv))
print('(S) stabilisation for A >= 2v_p(Y)+4(+3 at p=2): %d cases, %d violations'%(st,sv))
