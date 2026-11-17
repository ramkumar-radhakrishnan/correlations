import sys,time; sys.path.insert(0,"/tmp/claude-0/-home-user-correlations/0ab97fcf-94a2-5aec-8bad-586400481685/scratchpad/agent_bfkl")
from symdec import *
np.seterr(all="ignore")
B4=blocks("2d955658-4rho.tex"); BB=blocks("9ea88190-BFKL.tex"); R=lambda *ls:[rhs(BB[l]) for l in ls]
t=rhs(B4[64]); parts=t.split(r"+\rm{c.c}"); assert len(parts)==3
P1=parts[0]; Bm=parts[1]; k=Bm.rindex(r"\Bigg\{\frac{1}{(2\pi)^{3}}"); P2=Bm[:k]; P3=Bm[k+len(r"\Bigg\{"):]
def dlog(env,fun):
    C=env["const"]; out=[]
    for v in (1e8,1e12):
        C["vee"]=v; out.append(fun())
    C["vee"]=1e9
    return (out[1]-out[0])/np.log(1e4)
def four(env,pos):
    return dlog(env,lambda: inv(env,[P1,P3],pos,cc=True)+inv(env,[P2],pos,cc=False))
if __name__=="__main__":
    for seed in (1,2):
        env=setup(3,seed); rng=np.random.default_rng(100+seed); p=rng.normal(size=(3,2)); t0=time.time()
        op=add(ev(env,[P1,P3],cc=True),ev(env,[P2])); print("4rho.tex eq.(4rho): rho-point orderings present:",sorted(set(len(k) for k in op)),len(op))
        a3=four(env,list(p)); b3=dlog(env,lambda: inv(env,R(34,40),list(p)))
        print(seed,"deg-3 part: |4rho|=%.4f |Na+Nb (BFKL l.34,40)|=%.4f diff=%.2e"%(np.max(np.abs(a3)),np.max(np.abs(b3)),np.max(np.abs(a3-b3))))
        a2=four(env,list(p[:2])); L1=dlog(env,lambda: inv(env,R(47),list(p[:2]))); L2=dlog(env,lambda: inv(env,R(69),list(p[:2])))
        m=lambda q: np.max(np.abs(q))
        print(seed,"deg-2 part: |4rho|=%.4f |L47|=%.4f |L69|=%.4f ; 4rho-(L47+L69)=%.2e ; 4rho-L69=%.2e ; 4rho-L47=%.2e  (%.0fs)"%(m(a2),m(L1),m(L2),m(a2-L1-L2),m(a2-L2),m(a2-L1),time.time()-t0))
