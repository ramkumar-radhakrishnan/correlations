import sys; sys.path.insert(0,"/tmp/claude-0/-home-user-correlations/0ab97fcf-94a2-5aec-8bad-586400481685/scratchpad/agent_bfkl")
from grpcoef import *
G1=blocks("16e62d79-Group_I.tex"); G2=blocks("24c35375-Group_II.tex"); BB=blocks("9ea88190-BFKL.tex")
# Group II Row IV : BFKL term of Group_II.tex l.262 (first equality) vs BFKL.tex l.345
src=open(UP+"24c35375-Group_II.tex").read().split("\n")
t=src[261]+src[262]; t=t[t.index("&")+1:]; t=t.replace(r"\notag \\","")
for seed in (1,2):
    env=setup(3,seed); a=coef(env,[t]); b=coef(env,[rhs(BB[344])]); print("GII RowIV seed",seed,"|G|=%.4f |B|=%.4f diff %.1e"%(norm(a),norm(b),opdiff(a,b)))
# Group I Row II : what is the BFKL.tex bracket?
tb=rhs(BB[106]); i=tb.index(r"& -\,\frac{1}{(2\pi)^{3}}\,\frac{i\,g^{4}\,f^{hbc}}{32\pi^{4}}"); zterm,logterm=tb[:i],tb[i+1:]
for seed in (1,2):
    env=setup(3,seed); C=env["const"]
    c1=coef(env,[rhs(G1[44])]); c1b=coef(env,[rhs(G1[48])]); c2a=coef(env,[rhs(G1[60])])
    # coefficient of log(1/Lambda) of N_2a, two-Wilson-line colour structure only
    env["ufilter"]=2; o=[]
    for lam in (1e-3,1e-5):
        C["Lam"]=lam; o.append(evaluate(env,rhs(G1[60])))
    cl={k:(o[1][k]-o[0][k])/np.log(1e2)/N0*C["kp"] for k in o[0]}
    cb=coef(env,[logterm]); env["ufilter"]=None
    print("GI RowII seed",seed,"coef of log(vee): N1(l.45) %.1e, N1(l.49) %.1e, N2a(l.61) %.1e ; BFKL log-term %.4f ; [coef of log(k+/Lambda) of N2a, U(w)f{rho,rho} part] %.4f ; sum %.1e ; difference %.1e"%(
          norm(c1),norm(c1b),norm(c2a),norm(cb),norm(cl),norm(add(cb,cl)),opdiff(cb,cl)))
