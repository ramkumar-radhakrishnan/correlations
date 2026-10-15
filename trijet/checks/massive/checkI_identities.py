# Check I: the per-helicity contractions quoted in Sec. 3 (Regular x Instantaneous) of Trijet_simplified_massive.tex
import sympy as sp
I=sp.I; z,x=sp.symbols('zeta xi'); e=1-z-x
def eps2(a,b): return [[0,1],[-1,0]][a][b]
def d(a,b): return 1 if a==b else 0
def varphi(i,j,l,zq): return (2*zq-1)*d(i,j) + 2*I*eps2(i,j)*l
def tau(mm,n,l,xi_,par_d): return (2*par_d+xi_)*d(mm,n) + 2*I*xi_*eps2(mm,n)*l
for lam in (sp.Rational(1,2),sp.Rational(-1,2)):
    # (a) qbar reg x qbar inst* : indices (i,m) vector
    A=sp.Matrix(2,2,lambda i,m: sp.factor(sum(varphi(l,i,lam,z)*tau(r,m,lam,x,e)*sp.conjugate(d(l,r)+2*I*eps2(l,r)*lam) for l in (0,1) for r in (0,1))))
    # (b) q reg x q inst*: tau^{(vec m)(pol r)}(xi,zeta) phi^{(pol l)(vec i)}(zeta+xi)
    B=sp.Matrix(2,2,lambda i,m: sp.factor(sum(tau(m,r,lam,x,z)*varphi(l,i,lam,z+x)*sp.conjugate(d(l,r)-2*I*eps2(l,r)*lam) for l in (0,1) for r in (0,1))))
    # (c) q reg x qbar inst*
    C=sp.Matrix(2,2,lambda i,m: sp.factor(sum(tau(m,r,lam,x,z)*varphi(l,i,lam,z+x)*sp.conjugate(d(l,r)+2*I*eps2(l,r)*lam) for l in (0,1) for r in (0,1))))
    # (d) q inst x qbar reg*
    Dm=sp.Matrix(2,2,lambda j,n: sp.factor(sum((d(l,r)-2*I*eps2(l,r)*lam)*sp.conjugate(varphi(l,j,lam,z)*tau(r,n,lam,x,e)) for l in (0,1) for r in (0,1))))
    print("lam=",lam); 
    for nm,Mx in (("a",A),("b",B),("c",C),("d",Dm)): print("  ",nm, sp.simplify(Mx))
