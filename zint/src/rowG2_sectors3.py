"""Group II Row IV: the clean BFKL / DGLAP / finite split of the five-line form.

Five lines as written by the user:
  L1  P_N [log(V/Lam)+log(k+/V)] e^{-ikr} Int X [1/epsb - L_r] U          (no dxi, xi=1)
  L2  -(g^2/4pi^2)(1/epsb)  Int_1^{xi0} dxi/xi^2 (Pgg/2) dN_LO(k/xi)
  L3  -(g^2/4pi^2)          Int_1^{xi0} dxi/xi^2 (Pgg/2) [1 - L_r] dN_LO(k/xi)
  L4  +P_N                  Int_1^{xi0} dxi/xi^2 Int [ d_perp xi(1-xi) R
                                                      + (xi/[1-xi]_+ + (1-xi)/xi) X ] U
"""
import sympy as sp

Nc, epsb, Lr, x, dp = sp.symbols('N_c epsbarinv L_r xi d_perp')
XX, RR, plus = sp.symbols('XX RR plusterm')
xb = 1 - x
Chat = x*xb + plus + xb/x            # = Pgg/(2 Nc)
M = epsb + 1 - Lr

print("="*78)
print("1.  THE THREE-WAY SPLIT  (non-delta sector, measure restored to Int_{xi0}^1)")
print("="*78)
# L2+L3+L4 as written carry Int_1^{xi0} = -Int_{xi0}^1 ; flip the sign of each integrand
L23 = -( -Chat*epsb*XX - Chat*(1-Lr)*XX )          # -(L2+L3) integrands, /P_N
L4  = -(  dp*x*xb*RR + (plus + xb/x)*XX )          # -(L4) integrand,     /P_N
nondelta = sp.expand(L23 + L4)
print("  L2+L3 (flipped) :", sp.expand(L23))
print("  L4    (flipped) :", sp.expand(L4))
dglap  = Chat*(epsb - Lr)*XX
finite = x*xb*(XX - dp*RR)
print()
print("  proposed:  DGLAP = Chat [1/epsbar - L_r] X   and   FINITE = xi(1-xi)[X - d_perp R]")
print("  residual  (L2+L3+L4) - DGLAP - FINITE =", sp.simplify(nondelta - dglap - finite))

print()
print("="*78)
print("2.  WHERE THE '+1' BELONGS")
print("="*78)
print("  As written, L3 puts the +1 in the DGLAP sector: Pgg/2 [1 - L_r].")
print("  But L4 keeps -(xi/[1-xi]_+ + xibar/xi) X, and")
print("      (xi/[1-xi]_+ + xibar/xi) * (+1)  -  (xi/[1-xi]_+ + xibar/xi) =",
      sp.simplify((plus + xb/x)*1 - (plus + xb/x)))
print("  so two of the three terms of the +1 cancel L4's X block exactly.")
print("  Only  xi(1-xi) * (+1) * X  survives, and it pairs with -d_perp xi(1-xi) R.")
print("  => the +1 is NOT a DGLAP term.  It belongs to the finite sector.")
print("  leftover check:", sp.simplify(sp.expand(Chat*1*XX - (plus + xb/x)*XX - x*xb*XX)))

print()
print("="*78)
print("3.  WHY THE FINITE SECTOR IS FINITE: the tensor is TRACELESS")
print("="*78)
m = sp.symbols('m')
print("  X - d_perp R  is  X^m X'^{m'} [ delta^{mm'} - d_perp rhat^m rhat^{m'} ] / (...)")
print("  trace :  delta^{mm'} delta_{mm'} - d_perp rhat.rhat  =  d_perp - d_perp =", 0)
print("  The divergent part of G^{mm'} is  (pi/2)(log R^2/r^2 + 1) delta^{mm'} : PURE TRACE.")
print("  A traceless projector annihilates it, so this block cannot carry the log. Finite by structure.")

print()
print("="*78)
print("4.  L3 CANNOT BE WRITTEN AS A CONVOLUTION WITH dN_LO")
print("="*78)
print("  DGLAP factorisation needs the coefficient of dN_LO(k/xi) to be independent of the")
print("  transverse integration variables.  L2's coefficient is 1/epsbar : fine.")
print("  L3's coefficient is [1 - log(4 pi^2 mu^2 (y'-y)^2)] , and (y'-y) is integrated INSIDE")
print("  dN_LO.  So L3 as written is ill-defined; the log must sit inside Int_{x,x',y,y'} ,")
print("  next to X = (x-y).(x'-y')/[(x-y)^2 (x'-y')^2] .")

print()
print("="*78)
print("5.  NORMALISATION CONDITION ON L2/L3")
print("="*78)
g, pi, kp = sp.symbols('g pi kplus', positive=True)
PN = g**4/(8*pi**4)*Nc/kp
need = sp.simplify(PN*Chat/((g**2/(4*pi**2))*Nc*Chat))
print("  with the prefactor g^2/(4 pi^2) now used, L2/L3 match L1/L4 iff")
print("      dN_LO/d^3k = [1/(2pi)^3] *", need, "* Int_{x,x',y,y'} X U")
print("      i.e. dN_LO/d^3k = 1/(2pi)^3 * (g^2 / 2 pi^2) * (1/k+) * Int X U .")
print("  (the earlier g^2/(8 pi^2) version needed g^2/pi^2 -- a factor 2 apart; check your LO row)")
