import sys; sys.path.insert(0,"/tmp/claude-0/-home-user-correlations/0ab97fcf-94a2-5aec-8bad-586400481685/scratchpad/agent_bfkl")
from tworho import *
np.seterr(all="ignore")
def split_list(blk):
    r=rhs(blk); a=r.index(r"\Bigg["); b=r.rindex(r"\Bigg]"); head=r[:a]; body=r[a+6:b]
    parts=[p.strip() for p in re.split(r"\\notag\s*\\\\\s*&?",body) if p.strip()]
    out=[]
    for p in parts:   # a few lines hold two monomials separated by '\notag \\ &' already handled; keep as is
        out.append(head+r"\Bigg["+p+r"\Bigg]")
    return out
def fit(bracket_texts,list_block,cc,seeds=(1,2,3),eps=0.0,verbose=True):
    B=blocks("9ea88190-BFKL.tex"); mons=split_list(B[list_block])
    A=[];b=[]
    for seed in seeds:
        env=setup(3,seed); rng=np.random.default_rng(100+seed); u=rng.normal(size=2); v=rng.normal(size=2)
        A.append(np.array([Finv(env,[m],u,v).ravel() for m in mons]).T); b.append(Finv(env,bracket_texts,u,v,cc=cc,eps=eps).ravel())
    A=np.vstack(A); b=np.concatenate(b)
    A2=np.vstack([A.real,A.imag]); b2=np.concatenate([b.real,b.imag])
    c,res,rk,sv=np.linalg.lstsq(A2,b2,rcond=None)
    if verbose:
        print("rank",rk,"of",len(mons),"; residual %.2e ; |b|=%.3f"%(np.linalg.norm(A2@c-b2),np.linalg.norm(b2)))
        print("fitted multipliers of the listed monomials:",np.round(c,4))
    return c
if __name__=="__main__":
    B=blocks("9ea88190-BFKL.tex")
    from row2b import eq
