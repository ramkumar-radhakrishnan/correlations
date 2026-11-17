import sys; sys.path.insert(0,"/tmp/claude-0/-home-user-correlations/0ab97fcf-94a2-5aec-8bad-586400481685/scratchpad/agent_bfkl")
from tworho import *
np.seterr(all="ignore")
B=blocks("9ea88190-BFKL.tex"); t=rhs(B[106])
i=t.index(r"& -\,\frac{1}{(2\pi)^{3}}\,\frac{i\,g^{4}\,f^{hbc}}{32\pi^{4}}"); t1,t2=t[:i],t[i+1:]
for seed in (1,2):
    env=setup(3,seed); rng=np.random.default_rng(100+seed); u=rng.normal(size=2); v=rng.normal(size=2)
    L=Finv(env,[rhs(B[113])],u,v)
    for nm,tx in (("z-term",[t1]),("log-term",[t2])):
        for cc in (True,False):
            a=Finv(env,tx,u,v,cc=cc); fin=np.all(np.isfinite(a))
            if not fin: a=Finv(env,tx,u,v,cc=cc,eps=1e-6); a2=Finv(env,tx,u,v,cc=cc,eps=1e-4)
            print(seed,nm,"cc",cc,"finite at y=x:",fin,"|2rho|=%.4f |2rho-list|=%.3e"%(np.max(np.abs(a)),np.max(np.abs(a-L))), "" if fin else "eps-dependence %.2e"%np.max(np.abs(a-a2)))
