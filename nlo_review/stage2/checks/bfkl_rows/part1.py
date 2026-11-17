import sys; sys.path.insert(0,"/tmp/claude-0/-home-user-correlations/0ab97fcf-94a2-5aec-8bad-586400481685/scratchpad/agent_bfkl")
from grpcoef import *
G1=blocks("16e62d79-Group_I.tex"); G2=blocks("24c35375-Group_II.tex"); G3=blocks("96b5b029-Group_III.tex"); BB=blocks("9ea88190-BFKL.tex")
def T(B,l,fix=""): return rhs(B[l])+fix
cases=[("GI RowIII N1",[T(G1,85,r"\Bigg\}")],[T(BB,125)]),("GI RowIII N2",[T(G1,88)],[T(BB,132)]),
       ("GI RowIV",[T(G1,132)],[T(BB,187)]),
       ("GII RowI part1",[T(G2,48)],[T(BB,246),T(BB,252)]),
       ("GIII RowI",[T(G3,46)],[T(BB,356)]),("GIII RowII",[T(G3,57)],[T(BB,405)])]
for l,b in ((64,294),(67,299),(70,303),(73,308),(76,312)): cases.append(("GII RowI part2 G%d vs B%d"%(l,b),[T(G2,l)],[T(BB,b)]))
cases.append(("GII RowI part2 sum(64..76) vs sum(294..312)",[T(G2,l) for l in (64,67,70,73,76)],[T(BB,b) for b in (294,299,303,308,312)]))
if __name__=="__main__":
    for nm,g,b in cases:
        out=[]
        for seed in (1,2):
            env=setup(3,seed)
            try:
                cg=coef(env,g); cb=coef(env,b); cg2=coef(env,g,1e11,1e15)
                out.append("|G|=%.4f |B|=%.4f |G-B|=%.1e (vee-stab %.0e) keys %s"%(norm(cg),norm(cb),opdiff(cg,cb),opdiff(cg,cg2),"same" if set(k for k in cg if norm({k:cg[k]})>1e-9)==set(k for k in cb if norm({k:cb[k]})>1e-9) else "DIFF"))
            except Exception as e: out.append("ERR "+str(e)[:150].replace("\n"," "))
        print("%-45s"%nm," || ".join(out))
