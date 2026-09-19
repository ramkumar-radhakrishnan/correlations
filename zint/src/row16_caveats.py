"""Row 16: the 1/Lambda power divergence, and the missing on-shell (i pi) piece."""
import numpy as np, mpmath as mp, sympy as sp
mp.mp.dps=30; K=1.0
print("="*78); print("1.  PREFACTORS OF THE INPUT EXPRESSION (two-piece B_3)"); print("="*78)
g,pi=sp.symbols('g pi',positive=True); P,Kp=sp.symbols('pplus kplus',positive=True); S=P+Kp; I=sp.I
Abra=I*g/(sp.sqrt(2)*pi*sp.sqrt(Kp)); Aket=I*g/(sp.sqrt(2)*pi*sp.sqrt(P))
B1=sp.conjugate(I*g**2/(2*pi)**2*sp.sqrt(P*Kp)*S/(Kp-P)**2).subs(
     {sp.conjugate(g):g,sp.conjugate(pi):pi,sp.conjugate(Kp):Kp,sp.conjugate(P):P})
B3=sp.conjugate(-I*g**2/(4*pi**2*sp.sqrt(P*Kp))).subs(
     {sp.conjugate(g):g,sp.conjugate(pi):pi,sp.conjugate(Kp):Kp,sp.conjugate(P):P})
print("   piece I  : computed/quoted =",
      sp.simplify(sp.simplify(Abra*Aket*B1)/(I*g**4/(8*pi**4)*S/(Kp-P)**2)))
print("   piece III: computed/quoted =",
      sp.simplify(sp.simplify(Abra*Aket*B3)/(-I*g**4/(8*pi**4)/(P*Kp))))
print("   (the overall 1/(2pi)^3 sits outside; colour f^{abd}[U^{de}(x)U^{bc}(z) - U^{de}(x)U^{bc}(y)]")
print("    is the same structure verified before, with f legitimately an overall factor.)")

print(); print("="*78); print("2.  THE 1/Lambda IS A POWER DIVERGENCE, AND NOTHING CANCELS IT"); print("="*78)
print("   terms 1 and 4 both carry +1/Lambda with the SAME coefficient 1/(A-B), so they ADD:")
print("      total  ->  4/[Lambda (A-B)]   in units of i g^4 f/(16 pi^5)")
print()
print("   %-22s %16s %16s"%("(A,B,Lambda)","pieces 1+2 (PV)","4/[Lambda(A-B)]"))
def pv1(A,B,Lam,V,lo,hi):
    be=K*B/A; f=lambda p:(p+K)/((K-p)**2*(p*A-K*B)); d=1e-9
    if lo<be<hi: return mp.quad(f,[lo,be-d])+mp.quad(f,[be+d,hi])
    return mp.quad(f,[lo,hi])
for A,B,Lam,V in [(0.62,1.45,1e-5,60.),(1.80,0.55,1e-6,120.),(2.40,1.10,1e-6,300.)]:
    tot=float(pv1(A,B,Lam,V,K+Lam,V)+pv1(A,B,Lam,V,Lam,K-Lam))
    print("   %-22s %16.2f %16.2f"%("(%.2f,%.2f,%.0e)"%(A,B,Lam),tot,4/(Lam*(A-B))))
print()
print("   In the THREE-piece B_3 the (S/k+)x(S/k+) halves of piece II cancelled each other's")
print("   double pole; with piece II absent here there is nothing at all to cancel piece I's.")
print("   1/Lambda is NOT a rapidity logarithm: evolution cannot absorb it.")

print(); print("="*78); print("3.  THE MISSING ON-SHELL PIECE"); print("="*78)
print("   p+ A - k+ B vanishes at  beta = k+ B/A , which sits INSIDE the range for generic")
print("   transverse points.  The quoted result is the PRINCIPAL VALUE; the delta-function")
print("   half of  1/(x -/+ i eps) = PV -/+ i pi delta  is absent.")
print()
print("   residues (same units):")
Asy,Bsy,ksy,psy=sp.symbols('A B kplus p',positive=True)
be=ksy*Bsy/Asy
rI=sp.simplify(((psy+ksy)/((ksy-psy)**2*Asy)).subs(psy,be))
print("      piece I  :", sp.simplify(sp.factor(rI)), " = (A+B)/(k+ (A-B)^2)  -- the same")
print("                 coefficient that multiplies the logarithm, as it must be")
rIII=sp.simplify(((ksy*Bsy+psy*Asy)/(psy*Asy)/ksy).subs(psy,be))
print("      piece III:", sp.simplify(rIII), " = 2/k+")
print()
print("   so the full answer needs, in addition to what you wrote,")
print("      pieces 1,2 :  -/+ i pi (A+B)/(k+ (A-B)^2)   whenever beta lies in that range")
print("      piece 3    :  +/- 2 i pi / k+               (the quoted term carries an extra -1)")
print("   and since every prefactor already carries an explicit i, these delta pieces are REAL")
print("   and contribute to the cross section.  They are the on-shell intermediate state.")
print()
print("   check that beta really is inside:")
for A,B in [(0.62,1.45),(1.80,0.55),(2.40,1.10),(0.90,2.07)]:
    b=K*B/A
    print("      A=%.2f B=%.2f -> beta=%.4f : in [Lam,k-Lam] %-5s , in [k+Lam,V] %-5s , in [Lam,V] %s"
          %(A,B,b,b<K,b>K,True))
print()
print("="*78); print("4.  LOG ARGUMENTS CAN BE NEGATIVE"); print("="*78)
for A,B,Lam,V in [(0.62,1.45,1e-5,60.)]:
    print("   e.g. A=%.2f B=%.2f : (Lambda+k+)A - k+B = %+.4f  and  V A - k+B = %+.4f"
          %(A,B,(Lam+K)*A-K*B, V*A-K*B))
    print("        so the ratio inside your log is NEGATIVE.  The expressions only make sense")
    print("        as log| . | (which is what the principal value gives); write the bars.")
