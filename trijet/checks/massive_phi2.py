import sympy as sp
from sympy import I, simplify, expand, factor
z,x,mm = sp.symbols('zeta xi m', positive=True); e=1-z-x; zb=1-z
d=sp.eye(2); ep=sp.Matrix([[0,1],[-1,0]])
R  = sp.Matrix(sp.symbols('R_1 R_2', real=True));   Rb = sp.Matrix(sp.symbols('Rb_1 Rb_2', real=True))
Y  = sp.Matrix(sp.symbols('Y_1 Y_2', real=True));   Yb = sp.Matrix(sp.symbols('Yb_1 Yb_2', real=True))
def ehat(l): return sp.Matrix([1,-I*l])
def fhat(l): return sp.Matrix([1, I*l])

def paths(a,b,c,Pf,lam1):
    G_nf = sp.Matrix(2,2, lambda l,i: a*d[l,i]+I*lam1*ep[l,i])
    G_m  = -lam1*mm*ehat(lam1)
    En   = lambda lp: sp.Matrix(2,2, lambda r,m: b*d[r,m]-I*lp*c*ep[r,m])
    Em   = lambda lp: -lp*mm*(x**2/Pf)*fhat(lp)
    # each path -> (vector over l, vector over r, lam2)  after contracting free idx with R / Y
    return {
     'P1': (sp.Matrix([sum(G_nf[l,i]*R[i] for i in range(2)) for l in range(2)]),
            sp.Matrix([sum(En(-lam1)[r,m]*Y[m] for m in range(2)) for r in range(2)]), -lam1),
     'P2': (G_m, sp.Matrix([sum(En(+lam1)[r,m]*Y[m] for m in range(2)) for r in range(2)]), +lam1),
     'P3': (sp.Matrix([sum(G_nf[l,i]*R[i] for i in range(2)) for l in range(2)]), Em(-lam1), +lam1),
     'P4': (G_m, Em(+lam1), -lam1),
    }
def barred(M):    # same object built from the conjugate-amplitude vectors
    return M.subs({R[0]:Rb[0],R[1]:Rb[1],Y[0]:Yb[0],Y[1]:Yb[1]})
def cj(M): return M.applyfunc(lambda t: sp.conjugate(t).subs({sp.conjugate(s):s for s in (z,x,mm,*R,*Rb,*Y,*Yb)}))

print("=== antiquark channel: each path-pair, fully contracted with R,Rbar,Y,Ybar ===")
a_,b_,c_,P_ = 2*z-1, 2*e+x, x, zb
acc = {}
for lam1 in (+1,-1):
    A = paths(a_,b_,c_,P_,lam1)
    for kA,(GA,EA,l2A) in A.items():
        for kB,(GB,EB,l2B) in A.items():
            if l2A != l2B: continue
            GB2, EB2 = cj(barred(GB)), cj(barred(EB))
            val = sum(GA[l]*GB2[l] for l in range(2))*sum(EA[r]*EB2[r] for r in range(2))
            acc[(kA,kB)] = sp.expand(acc.get((kA,kB),0)+val)
for k in sorted(acc):
    v = sp.simplify(sp.factor(sp.simplify(acc[k])))
    if v!=0: print(f"  {k[0]}x{k[1]}: {v}")

print("\n=== LO-dijet normalisation check (no emission): photon vertex only ===")
lo = 0
for lam1 in (+1,-1):
    G_nf = sp.Matrix(2,2, lambda l,i: (2*z-1)*d[l,i]+I*lam1*ep[l,i])
    Gm   = -lam1*mm*ehat(lam1)
    v1 = sp.Matrix([sum(G_nf[l,i]*R[i] for i in range(2)) for l in range(2)])
    v1b= cj(barred(v1)); lo += sum(v1[l]*v1b[l] for l in range(2))
    lo += sum(Gm[l]*cj(barred(Gm))[l] for l in range(2))
print("  sum_lam1,l |V_gamma|^2 =", sp.factor(sp.simplify(sp.expand(lo))))
print("  expect 4[zeta^2+zetabar^2] (R.Rbar) + 4 m^2  ->",
      sp.simplify(sp.expand(lo) - (4*(z**2+zb**2)*(R.dot(Rb)) + 4*mm**2)))

print("\n=== soft limit xi -> 0 of each path-pair ===")
for k in sorted(acc):
    v0 = sp.simplify(sp.limit(sp.simplify(acc[k]), x, 0))
    print(f"  {k[0]}x{k[1]}:  {sp.factor(v0)}")
