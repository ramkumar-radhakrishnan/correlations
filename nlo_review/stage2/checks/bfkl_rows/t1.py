import sys; sys.path.insert(0,"/tmp/claude-0/-home-user-correlations/0ab97fcf-94a2-5aec-8bad-586400481685/scratchpad/agent_bfkl")
from texeval import *
B=blocks("9ea88190-BFKL.tex")
env=setup(3,1)
for l in sorted(B):
    try:
        r=evaluate(env,rhs(B[l])); print(l,{k:round(norm({k:v}),4) for k,v in r.items()})
    except Exception as e: print(l,"ERR",type(e).__name__,str(e)[:230])
print(set(env.get("typos",[])))
