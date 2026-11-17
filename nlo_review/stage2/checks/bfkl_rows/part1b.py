import sys; sys.path.insert(0,"/tmp/claude-0/-home-user-correlations/0ab97fcf-94a2-5aec-8bad-586400481685/scratchpad/agent_bfkl")
from grpcoef import *
G2=blocks("24c35375-Group_II.tex"); BB=blocks("9ea88190-BFKL.tex")
fix=lambda t: t.replace(r"U^{ac}(w)\,^{c'b}(z)",r"U^{ac}(w)\,U^{c'b}(z)")
for nm,g,b in [("2a (G64) alone",64,None),("2b G67 vs B294",67,294),("2c G70 vs B299",70,299),("2d G73 vs B303",73,303),("2e G76 vs B308",76,308),("2f G79 vs B312",79,312)]:
    out=[]
    for seed in (1,2):
        env=setup(3,seed); cg=coef(env,[fix(rhs(G2[g]))]); cb=coef(env,[rhs(BB[b])]) if b else {}
        tot=norm({k:v for k,v in evaluate(env,fix(rhs(G2[g]))).items()})/N0
        out.append("|G|=%.4g |B|=%.4g |G-B|=%.1e (size of full G expr %.2g)"%(norm(cg),norm(cb),opdiff(cg,cb),tot))
    print("%-18s"%nm," || ".join(out))
