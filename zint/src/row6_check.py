"""Row 6 from scratch: prefactor, colour, transverse scan."""
import sympy as sp, numpy as np
g, pi = sp.symbols('g pi', positive=True); kp, pp = sp.symbols('kplus pplus', positive=True); I = sp.I

print("="*70); print("1. PREFACTOR, built from the three coefficient functions")
print("="*70)
print("   A^dag_i(k+,-w')      : -i g /(sqrt2 pi sqrt(k+))")
print("   Abar^dag_j(p+,-z)    : -i g /(sqrt2 pi sqrt(p+))")
print("   B_2(k+,-w; p+,-z)    : +i g^2 sqrt(p+ k+) / [4 pi^2 (...) (p+ + k+)]")
print("   out front            : 2/(2pi)^3")
Adk  = -I*g/(sp.sqrt(2)*pi*sp.sqrt(kp))
Adp  = -I*g/(sp.sqrt(2)*pi*sp.sqrt(pp))
Bpre =  I*g**2*sp.sqrt(pp*kp)/(4*pi**2*(pp+kp))
tot  = sp.simplify(sp.Integer(2)/(2*pi)**3 * Adk * Adp * Bpre)
print(f"\n   product = {tot}")
print(f"   i-count : (-i)(-i)(+i) = {sp.simplify((-I)*(-I)*I)}")
numeric_correct = sp.simplify(sp.Integer(2)/(2*pi)**3/(sp.sqrt(2)*pi)**2/(4*pi**2))
print(f"   numeric = {numeric_correct}   =  (1/(2pi)^3) x {sp.simplify(numeric_correct*(2*pi)**3)}")
print(f"   QUOTED  = 1/(4 pi^4)  (no 1/(2pi)^3)  = {sp.simplify(1/(4*pi**4))}")
print(f"   ratio quoted/correct = {sp.simplify((1/(4*pi**4))/numeric_correct)}")
print(f"   k+ powers -> {sp.simplify(sp.sqrt(pp*kp)/(sp.sqrt(kp)*sp.sqrt(pp)*(pp+kp)))}"
      "   (divides B's bracket: gives 1/(2(z-w)^2(p++k+)), 1/k+, 1/p+ )")

print()
print("="*70); print("2. COLOUR, term by term"); print("="*70)
print("   BRA  [ -Abar^dag,a_i(k+,-w') + A^dag,b'_i(k+,-w') U^{b'a}(w') ]")
print("        -Abar^dag,a  ->  -U^{ab'}(x') rho^{b'}(x')")
print("        +A^dag,b' U^{b'a}(w')  ->  + U^{b'a}(w') rho^{b'}(x')")
print("        so the bra colour is  [ -U^{ab'}(x') + U^{b'a}(w') ] rho^{b'}(x')")
print("        pulling out the minus:  -[ U^{ab'}(x') - U^{b'a}(w') ] rho^{b'}(x')")
print("        ** your starting line has U^{b'a}(w'); your answer has U^{ab'}(w') -- transposes **")
print()
print("   KET term 1: U^{ac}(w) U^{bd}(z) Abar^dag,b_j(p+,-z) B^{dc}_{2ji}")
print("        Abar^dag,b -> U^{be}(y) rho^e(y) ;  B^{dc} carries f^{c d a} rho^a(x) -> f^{cdc'} rho^{c'}(x)")
print("        => f^{cdc'} U^{ac}(w) U^{bd}(z) U^{be}(y) rho^e(y) rho^{c'}(x)      MATCHES")
print("   KET term 2: -Abar^dag,b_j(p+,-z) Bbar^{ba}_{2ji}")
print("        B^{ba} carries f^{a b a'} rho^{a'}(x); barred -> f^{abc} U^{cd}(x) rho^d(x)")
print("        => -f^{abc} U^{be}(y) U^{cd}(x) rho^e(y) rho^d(x)                   MATCHES")

print()
print("="*70); print("3. TRANSVERSE SCAN (fresh, for this row's integrand)"); print("="*70)
def T(i,j,x,z,w,xi):
    A=(x-z)@(x-z); B=(x-w)@(x-w); zw=z-w; zw2=zw@zw; xz,xw=x-z,x-w
    d = 1.0 if i==j else 0.0
    return ( d*(A-B)/(2*zw2*(xi+1))
           + (xz[j]*zw[i]/zw2 - xz[j]*xw[i]/(2*B))          # the 1/k+ bracket (k+ = 1)
           + (xw[i]*zw[j]/zw2 + xz[j]*xw[i]/(2*A))/xi )     # the 1/p+ bracket
def F(x,y,z,w,V,xi):
    A=(x-z)@(x-z); B=(x-w)@(x-w); yz=y-z
    return sum(V[i]*yz[j]/(yz@yz)*T(i,j,x,z,w,xi) for i in range(2) for j in range(2))/(xi*A+B)
x0=np.array([0.3,-0.2]); y0=np.array([-0.6,0.9]); z0=np.array([0.8,0.5])
w0=np.array([-0.4,-0.7]); V0=np.array([0.4,0.7]); xi=0.4
th=np.linspace(0,2*np.pi,8192,endpoint=False); u=np.stack([np.cos(th),np.sin(th)],-1)
for nm,mv,c in (("x -> z",'x',z0),("x -> w",'x',w0),("z -> w",'z',w0),("y -> z",'y',z0)):
    vals=[]
    for rho in (1e-2,1e-4,1e-6):
        a=dict(x=x0,y=y0,z=z0,w=w0)
        s=[]
        for uu in u:
            a2=dict(a); a2[mv]=c+uu*rho
            s.append(F(a2['x'],a2['y'],a2['z'],a2['w'],V0,xi))
        vals.append(np.mean(s))
    print(f"   {nm:8s}  <F> = " + "  ".join(f"{v:+.6e}" for v in vals) + "   -> O(1), NO UV")
print("   large |z| (no phase on z):", end=" ")
for R in (1e3,1e5):
    print(f"R^2<F>({R:.0e}) = {np.mean([F(x0,y0,uu*R,w0,V0,xi) for uu in u])*R**2:+.8f}", end="  ")
print("  -> IR log")
