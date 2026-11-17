import sys; sys.path.insert(0,"/tmp/claude-0/-home-user-correlations/0ab97fcf-94a2-5aec-8bad-586400481685/scratchpad/agent_bfkl")
from tworho import *
np.seterr(all="ignore")
B=blocks("9ea88190-BFKL.tex"); R=lambda *ls:[rhs(B[l]) for l in ls]
for seed in (1,2):
    env=setup(3,seed); rng=np.random.default_rng(100+seed); u=rng.normal(size=2); v=rng.normal(size=2)
    L1=Finv(env,R(47),u,v); L2=Finv(env,R(69),u,v)
    for cc in (False,True):
        a=Finv(env,R(34),u,v,cc=cc); b=Finv(env,R(40),u,v,cc=cc)
        m=lambda q: np.max(np.abs(q))
        print(seed,"cc",cc,"|Na|=%.3f |Nb|=%.3f |L47|=%.3f |L69|=%.3f"%(m(a),m(b),m(L1),m(L2)),
              "Na-L47 %.2e Nb-L69 %.2e Na-L69 %.2e Nb-L47 %.2e sum-sum %.2e"%(m(a-L1),m(b-L2),m(a-L2),m(b-L1),m(a+b-L1-L2)))
print("sign combinations  L47 vs sa*Na+sb*Nb (no cc):")
for seed in (1,2):
    env=setup(3,seed); rng=np.random.default_rng(100+seed); u=rng.normal(size=2); v=rng.normal(size=2)
    L1=Finv(env,R(47),u,v); a=Finv(env,R(34),u,v); b=Finv(env,R(40),u,v)
    print(seed,{(sa,sb):float("%.2e"%np.max(np.abs(sa*a+sb*b-L1))) for sa in (1,-1) for sb in (1,-1)})
