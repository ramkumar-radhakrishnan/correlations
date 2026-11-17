# Group III rows left as p+ integrals: coefficient of log(vee/k+) from large p+ after y' -> w' = y' + (p+/k+)(y'-z)
import sys; sys.path.insert(0,"/tmp/claude-0/-home-user-correlations/0ab97fcf-94a2-5aec-8bad-586400481685/scratchpad/agent_bfkl")
from grpcoef import *
G3=blocks("96b5b029-Group_III.tex"); BB=blocks("9ea88190-BFKL.tex")
PRE=re.compile(r"[+-]?\s*(?:\\,)?\s*\\frac\{1\}\{\(2\\pi\)\^\{3\}\}")
def pieces(t):
    st=[m.start() for m in PRE.finditer(t)]
    return [t[a:b] for a,b in zip(st,st[1:]+[len(t)])]
def integrand(t):
    t=re.sub(r"\\int_\{0\}\^\{\\vee\s*-\s*k\^\{\+\}\}","",t)
    t=t.replace(r"\frac{dp^{+}}{2\pi}",r"\frac{1}{2\pi}")
    t,n1=re.subn(r"e\^\{-ik\\cdot\s*\\frac\{p\^\{\+\}\}\{k\^\{\+\}\}(?:\\,|\s|\\left)*\(y'-z(?:\\right)?\)\}","",t)
    t,n2=re.subn(r"\\left\[\\frac\{1\}\{p\^\{\+\}\}\\right\]_\{\+\}",r"\\frac{1}{p^{+}}",t)
    assert n1>=1 and n2>=1,(n1,n2)
    return t
def run(name,gl,bl,nint=1):
    for seed in (1,2):
        env=setup(3,seed); C=env["const"]; X=env["X"]; kp=C["kp"]
        g=pieces(rhs(G3[gl])); b=pieces(rhs(BB[bl]))
        gint=[integrand(x) for x in g[:nint]]; gblue=g[nint:]
        sc=lambda op,f: {k:v*f/N0*kp for k,v in op.items()}
        # (1) explicit blue log terms
        c_blue=coef(env,gblue)
        # (2) subtraction term of the plus prescription:  -lim_{p->0} p * integrand
        C["pp"]=1e-10; r0={}
        for x in gint: r0=add(r0,evaluate(env,x))
        c_sub=sc(r0,-C["pp"])
        # (3) large p+ after the change of variable (Jacobian (k+/(p+ + k+))^2), y' = (k+ w' + p+ z)/(p+ + k+)
        ysave=X["y'"].copy(); res=[]
        for pp in (1e6,1e8):
            C["pp"]=pp; X["y'"]=(kp*X["w'"]+pp*X["z"])/(pp+kp); env["Ucache"].clear(); o={}
            for x in gint: o=add(o,evaluate(env,x))
            res.append(sc(o,pp*(kp/(pp+kp))**2))
        X["y'"]=ysave; C["pp"]=3.0
        cb=[coef(env,[x]) for x in b]
        nb=len(b)
        msg="seed %d: BFKL pieces %d | blue: |G|=%.4f |B1|=%.4f diff %.1e | plus-subtraction: |G|=%.4f |B2|=%.4f diff %.1e | B1+B2=%.1e | large p+: |G|=%.4f |B3..|=%.4f diff %.1e (p+ stab %.0e)"%(
            seed,nb,norm(c_blue),norm(cb[0]),opdiff(c_blue,cb[0]),norm(c_sub),norm(cb[1]),opdiff(c_sub,cb[1]),norm(add(cb[0],cb[1])),
            norm(res[1]),norm(cb[2]),opdiff(res[1],cb[2]),opdiff(res[0],res[1]))
        print(name,msg)
if __name__=="__main__":
    run("GIII RowIII part1",74,459); run("GIII RowIII part2",88,486); run("GIII RowV",111,520); run("GIII RowVI",126,546)
