# Stabilisation in the exponent: for p | t, is the vanishing (on class r mod t) of term gamma at k = t*p^j
# the same for all j >= 1 ?  (Eventual constancy is what the composite-level descent needs.)
import mpmath as mp, math
from level_test import amps, van
mp.mp.dps=25
report=[]
for t,m,p,J in ((6,5,2,3),(6,5,3,2),(10,3,2,3),(12,3,2,3),(6,7,2,3)):
    A={j:amps(t,m,t*p**j) for j in range(0,J+1)}
    changes=0; tot=0
    for r in range(t):
        for gam in A[0]:
            pat=[van(A[j][gam],t*p**j,t,r) if gam in A[j] else None for j in range(J+1)]
            pat=[x for x in pat if x is not None]
            tot+=1
            if len(set(pat[1:]))>1: changes+=1; report.append((t,m,p,r,str(gam),pat))
    print('t=%d m=%d p=%d: %d (class,term) pairs; vanishing pattern not constant for j>=1 in %d'%(t,m,p,tot,changes),flush=True)
print(report[:6])
