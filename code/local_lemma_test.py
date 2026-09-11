# Direct test of the local lemma used in Prop 5.4 (stationary phase with explicit threshold):
# T^(B)(x,y) = sum_{h mod q^B, q not| h} psi(h) theta(h) e((x h + y hbar)/q^B),
# psi real char of conductor | 8 (q=2) or | q (q odd); theta(h)=e(-c h/8) for q=2.
# Claim: for B >= 2(v_q(y)+e_q)+kappa_q+2  (e_q=3,kappa_q=3 for q=2; e_q=1,kappa_q=1 odd q):
#   T != 0  <=>  v_q(x)=v_q(y) and (x/q^v)(y/q^v) is a square mod q (mod 8 if q=2).
import cmath, math, random
def vq(x,q):
    if x==0: return 99
    v=0
    while x%q==0: x//=q; v+=1
    return v
def issq(u,q):
    if q==2: return u%8==1
    return pow(u%q,(q-1)//2,q)==1
def chars(q):
    if q==2: return [lambda h:1,lambda h:(1 if h%4==1 else -1),lambda h:(1 if h%8 in(1,7) else -1),lambda h:(1 if h%8 in(1,3) else -1)]
    return [lambda h:1,lambda h:(1 if pow(h,(q-1)//2,q)==1 else -1)]
def T(q,B,x,y,psi,c):
    M=q**B; s=0
    for h in range(1,M):
        if h%q==0: continue
        s+=psi(h)*cmath.exp(2j*math.pi*(((x*h+y*pow(h,-1,M))%M)/M - (c*h/8 if q==2 else 0)))
    return s
random.seed(3); tot=0; bad=[]
for q,Bs in ((2,(11,12,13)),(3,(7,8)),(5,(5,6)),(7,(4,5))):
    for psi in chars(q):
        for _ in range(12):
            vy=random.randint(0,1); y=q**vy*random.choice([u for u in range(1,60) if u%q])
            e_q,k_q=(3,3) if q==2 else (1,1)
            thr=2*(vy+e_q)+k_q+2
            c=random.randrange(8) if q==2 else 0
            vx=random.choice([vy,vy,vy+1,max(0,vy-1)]); x=q**vx*random.choice([u for u in range(1,60) if u%q])
            for B in Bs:
                if B<thr: continue
                val=T(q,B,x,y,psi,c)
                pred=(vq(x,q)==vq(y,q)) and issq((x//q**vx)*(y//q**vy),q)
                tot+=1
                if (abs(val)>1e-6)!=pred: bad.append((q,B,x,y,c,abs(val),pred))
print('local lemma (explicit threshold): %d cases, %d violations'%(tot,len(bad))); print(bad[:5])
