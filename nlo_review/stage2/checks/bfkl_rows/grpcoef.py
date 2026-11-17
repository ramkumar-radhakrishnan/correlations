import sys; sys.path.insert(0,"/tmp/claude-0/-home-user-correlations/0ab97fcf-94a2-5aec-8bad-586400481685/scratchpad/agent_bfkl")
from tworho import *
np.seterr(all="ignore")
def coef(env,texts,v1=1e9,v2=1e13):
    """coefficient of log(vee) of the operator, pointwise, from two large values of vee"""
    C=env["const"]; o=[]
    for v in (v1,v2):
        C["vee"]=v; tot={}
        for t in texts: tot=add(tot,evaluate(env,t))
        o.append(tot)
    C["vee"]=1e9
    return {k:(o[1].get(k,0)-o[0].get(k,0))/np.log(v2/v1)/N0*env["const"]["kp"] for k in set(o[0])|set(o[1])}
def opdiff(a,b): return norm(add(a,b,-1))
if __name__=="__main__":
    for fn in sys.argv[1:]:
        B=blocks(fn); env=setup(3,1); print("==",fn)
        for l in sorted(B):
            try:
                c=coef(env,[rhs(B[l])]); print(l,len(B[l]),{k:float("%.4g"%norm({k:v})) for k,v in c.items()})
            except Exception as e: print(l,len(B[l]),"ERR",type(e).__name__,str(e)[:150].replace("\n"," "))
