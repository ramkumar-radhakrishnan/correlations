import sympy as sp
from sympy import I
z,x,mm = sp.symbols('zeta xi m', positive=True); e=1-z-x; zb=1-z; eb=1-e
d=sp.eye(2); ep=sp.Matrix([[0,1],[-1,0]])
# amplitude-side and conjugate-side transverse vectors (generic names)
A1,A2,B1,B2 = sp.symbols('A_1 A_2 B_1 B_2', real=True)     # amp: A = R or Rtilde ; B = Y or X
C1,C2,D1,D2 = sp.symbols('C_1 C_2 D_1 D_2', real=True)     # conj: C = Rbar ; D = Ybar or Xbar
A=sp.Matrix([A1,A2]); B=sp.Matrix([B1,B2]); C=sp.Matrix([C1,C2]); D=sp.Matrix([D1,D2])
def eh(l): return sp.Matrix([1,-I*l])
def fh(l): return sp.Matrix([1, I*l])
def legs(a,b,c,Pf,lam1,V,W):
    Gn = sp.Matrix(2,2, lambda l,i: a*d[l,i]+I*lam1*ep[l,i])
    Gm = -lam1*mm*eh(lam1)
    En = lambda lp: sp.Matrix(2,2, lambda r,mI: b*d[r,mI]-I*lp*c*ep[r,mI])
    Em = lambda lp: -lp*mm*(x**2/Pf)*fh(lp)
    gv = sp.Matrix([sum(Gn[l,i]*V[i] for i in range(2)) for l in range(2)])
    return {'P1':(gv, sp.Matrix([sum(En(-lam1)[r,mI]*W[mI] for mI in range(2)) for r in range(2)]), -lam1),
            'P2':(Gm, sp.Matrix([sum(En(+lam1)[r,mI]*W[mI] for mI in range(2)) for r in range(2)]), +lam1),
            'P3':(gv, Em(-lam1), +lam1),
            'P4':(Gm, Em(+lam1), -lam1)}
def cj(M): return M.applyfunc(lambda t: sp.conjugate(t).subs({sp.conjugate(s):s for s in (z,x,mm,A1,A2,B1,B2,C1,C2,D1,D2)}))

def run(name, amp, con):
    acc={}
    for lam1 in (+1,-1):
        LA = legs(*amp, lam1, A, B); LB = legs(*con, lam1, C, D)
        for kA,(GA,EA,lA) in LA.items():
            for kB,(GB,EB,lB) in LB.items():
                if lA!=lB: continue
                val = sum(GA[l]*cj(GB)[l] for l in range(2))*sum(EA[r]*cj(EB)[r] for r in range(2))
                acc[(kA,kB)] = sp.expand(acc.get((kA,kB),0)+val)
    print(f"\n### {name}")
    for k in sorted(acc):
        v = sp.factor(sp.simplify(acc[k]))
        if v!=0: print(f"   {k[0]}x{k[1]}: {v}")
    print("   --- soft limit xi->0 ---")
    for k in sorted(acc):
        v0 = sp.factor(sp.simplify(sp.limit(sp.simplify(acc[k]), x, 0)))
        if v0!=0: print(f"   {k[0]}x{k[1]}: {v0}")
    return acc

Aqb = (2*z-1,   2*e+x, +x, zb)      # antiquark emits : a=2 zeta-1, b=2eta+xi, c=+xi, parent 1-zeta
Aq  = (1-2*e,   2*z+x, -x, eb)      # quark emits     : a=2(1-eta)-1, b=2zeta+xi, c=-xi, parent 1-eta
run("QUARK-emission channel  (their 1.10)", Aq, Aq)
run("INTERFERENCE  (amp: quark emits ; conj: antiquark emits)  (their 1.18)", Aq, Aqb)
