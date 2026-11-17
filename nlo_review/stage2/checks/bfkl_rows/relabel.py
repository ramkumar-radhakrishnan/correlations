import sys; sys.path.insert(0,"/tmp/claude-0/-home-user-correlations/0ab97fcf-94a2-5aec-8bad-586400481685/scratchpad/agent_bfkl")
from tworho import *
np.seterr(all="ignore")
B=blocks("9ea88190-BFKL.tex"); R=lambda *ls:[rhs(B[l]) for l in ls]
rows=[("GI RowIV",(187,),(195,229),True),("GIII RowI",(356,),(363,393),True),("GIII RowVI",(546,),(553,),True),("GIII RowII",(405,),(412,445),True),("GII RowI p1",(246,252),(259,),False)]
for name,br,ls,cc in rows:
    env=setup(3,1); rng=np.random.default_rng(7); u=rng.normal(size=2); v=rng.normal(size=2)
    L=F(env,R(*ls),u,v)/N0            # list read literally: rho(x') at u, rho(x) at v
    types=[(o,c) for c in ((False,"h") if cc else (False,)) for o in ("x","x'","y")]
    opt={}
    for o,c in types:
        opt[(o,c)]=(F(env,R(*br),u,v,cc=c,odd_only=o)/N0, F(env,R(*br),v,u,cc=c,odd_only=o).T/N0)   # (pair->x', odd->x) , (pair->x, odd->x')
    best=[]
    for ch in itertools.product((0,1),repeat=len(types)):
        tot=sum(opt[t][c] for t,c in zip(types,ch)); best.append((np.max(np.abs(tot-L)),ch))
    best.sort(); 
    print(name,"best pointwise relabelling: diff %.1e"%best[0][0],{("%s%s"%(t[0],"+hc" if t[1] else "")):("pair->x',other->x" if c==0 else "pair->x,other->x'") for t,c in zip(types,best[0][1])},"| #choices with diff<1e-9:",sum(1 for d,_ in best if d<1e-9),"| next best %.1e"%best[1][0])
