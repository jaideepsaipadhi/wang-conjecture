# Composite t: does vanishing of an amplitude (cusp k, term gamma) on a class r mod t propagate to k*c'
# for c' coprime to t (the analogue of the Level Lemma's CRT step)?
import mpmath as mp, math
from level_test import amps, van
mp.mp.dps=30
viol=[]; checks=0
for t,m in ((12,4),(12,3),(10,3),(10,4),(15,2),(14,2)):
    ks=[k for k in range(1,t+1) if t%k==0 and math.gcd(k,t)**2>t]
    for k in ks:
        A=amps(t,m,k)
        for cp in (5,7,11,13):
            if math.gcd(cp,t)!=1: continue
            B=amps(t,m,k*cp)
            for r in range(t):
                for gam,lst in A.items():
                    if gam not in B: continue
                    checks+=1
                    if van(lst,k,t,r) and not van(B[gam],k*cp,t,r): viol.append((t,m,k,cp,r,str(gam)))
print('composite t: checks',checks,' violations of "vanishing at k => vanishing at k*c\' (c\' coprime to t)":',viol[:8])
