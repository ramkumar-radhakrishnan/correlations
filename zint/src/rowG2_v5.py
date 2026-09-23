"""Group II Row IV, v5: check (0.25)->(0.26)->(0.27) and the further simplification."""
import sympy as sp
Nc,g,pi,K,eps,Lr=sp.symbols('N_c g pi kplus epsbarinv L_r')
XX,RR,dp,x=sp.symbols('XX RR d_perp xi')
xb=1-x
plus=sp.Symbol('plusterm')                      # stands for xi/[1-xi]_+
Chat=x*xb + plus + xb/x                         # = Pgg/(2 Nc)
Pgg=2*Nc*Chat
P0 =1/(2*pi)**3*g**4/(8*pi**4)*Nc/K             # with N_c
P0p=1/(2*pi)**3*g**4/(8*pi**4)/K                # without N_c
M=eps+1-Lr
# verified ascending form (evolution term aside)
mine = P0*( Chat*M*XX - dp*x*xb*RR - (plus+xb/x)*XX )
print("="*74); print("1.  (0.25), ascending"); print("="*74)
e25 = P0*Chat*M*XX - P0*( dp*x*xb*RR + (plus+xb/x)*XX )
print("   (0.25) - verified =", sp.simplify(sp.expand(e25-mine)))
print()
print("="*74); print("2.  (0.26), ascending  (Nc Chat -> Pgg/2, prefactor loses Nc)"); print("="*74)
e26 = P0p*(Pgg/2)*M*XX - P0*( dp*x*xb*RR + (plus+xb/x)*XX )
print("   (0.26) - (0.25)   =", sp.simplify(sp.expand(e26-e25)))
print()
print("="*74); print("3.  (0.27), ascending  (split M = 1/epsbar + (1-Lr))"); block=None; print("="*74)
LOfac=1/(2*pi)**3*g**2/(pi**2*K)
e27 = g**2/(8*pi**2)*eps*(Pgg/2)*LOfac \
      - P0*( dp*x*xb*RR + (plus+xb/x)*XX ) \
      + P0p*(Pgg/2)*(1-Lr)*XX
print("   note the DGLAP block is written per unit N_LO ; restoring XX :")
e27x = P0p*(Pgg/2)*eps*XX - P0*( dp*x*xb*RR + (plus+xb/x)*XX ) + P0p*(Pgg/2)*(1-Lr)*XX
print("   (g^2/8pi^2) x LOfac - P0' =", sp.simplify(g**2/(8*pi**2)*LOfac - P0p), "  (0 means the DGLAP prefactor is right)")
print("   (0.27) - (0.26)   =", sp.simplify(sp.expand(e27x-e26)))
print()
print("="*74); print("4.  FURTHER SIMPLIFICATION: rows 3 and 4 of (0.27) collapse"); print("="*74)
rows34 = - P0*( dp*x*xb*RR + (plus+xb/x)*XX ) + P0p*(Pgg/2)*(1-Lr)*XX
compact = P0*( -Lr*Chat*XX + x*xb*(XX - dp*RR) )
print("   rows 3+4 - [ -Lr Chat XX + xi xibar (XX - d_perp RR) ] =",
      sp.simplify(sp.expand(rows34-compact)))
print("   => the xi/[1-xi]_+ and xibar/xi entries disappear from the explicit bracket;")
print("      they are already inside Chat = Pgg/(2 Nc).")
