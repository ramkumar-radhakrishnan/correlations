"""Group II Row IV: check the advisor's rearranged final result against the derived one.

Advisor's form (verbatim structure):
  T1 = +P_N [log(V/Lam)+log(k+/V)] Int_{x..}  X  [1/epsb - L_r]  U          (no dxi)
  T2 = -(g^2/8pi^2)(1/epsb) Int_1^{xi0} dxi/xi^2 (Pgg/2) dN_LO(k/xi)
  T3 = +P_N Int_1^{xi0} dxi/xi^2 [ d_perp xi(1-xi) R + (xi/[1-xi]_+ + (1-xi)/xi) X ] U
  T4 = -P_1 Int_1^{xi0} dxi/xi^2 X (Pgg/2) [1 - L_r] U
with P_N = (1/(2pi)^3)(g^4/8pi^4)(N_c/k+), P_1 = same with 1 instead of N_c,
and xi0 = k+/V < 1, so Int_1^{xi0} = - Int_{xi0}^1 .

Derived (verified) form, all xi-integrals written Int_{xi0}^1 :
  D  = P_N { delta(1-xi) log(k+/Lam) Chat_delta + [ Chat M X - d_perp xi xib R
                                                    - (xi/[1-xi]_+ + xib/xi) X ] } U
  with Chat = xi xib + xi/[1-xi]_+ + xib/xi = Pgg/(2Nc),  M = 1/epsb + 1 - L_r .
"""
import sympy as sp

Nc, epsb, Lr, x, dp = sp.symbols('N_c epsbarinv L_r xi d_perp')
XX, RR, plus = sp.symbols('XX RR plusterm')
xb = 1 - x
Chat = x*xb + plus + xb/x           # = Pgg/(2 Nc)
Pgg2 = Nc*Chat                      # = Pgg/2
M = epsb + 1 - Lr

print("="*78)
print("1.  THE NON-DELTA SECTOR, TERM BY TERM  (common factor P_N = ... N_c/k+ stripped)")
print("="*78)

# derived integrand, with the measure written Int_{xi0}^{1}
derived = Chat*M*XX - dp*x*xb*RR - (plus + xb/x)*XX

# advisor integrand as literally written, i.e. with the measure Int_{1}^{xi0}
adv_as_written = ( -Pgg2/Nc*epsb*XX                      # T2, after /N_c to match P_N
                   + dp*x*xb*RR + (plus + xb/x)*XX       # T3
                   - Pgg2/Nc*(1 - Lr)*XX )               # T4, P_1*Pgg/2 = P_N*Chat
print("  advisor (as written, measure Int_1^{xi0}) :", sp.simplify(sp.expand(adv_as_written)))
print("  derived  (measure Int_{xi0}^1)            :", sp.simplify(sp.expand(derived)))
print("  sum (advisor + derived)                   :", sp.simplify(sp.expand(adv_as_written + derived)))
print()
print("  -> advisor_integrand = - derived_integrand, EXACTLY.")
print("     And Int_1^{xi0} = - Int_{xi0}^1 because xi0 = k+/V < 1.")
print("     So the two expressions are THE SAME. The whole sign sits in the reversed limits.")
print()
print("  difference after flipping the limits:",
      sp.simplify(sp.expand((-adv_as_written) - derived)))

print()
print("="*78)
print("2.  THE NORMALISATION OF THE ADVISOR'S COUNTERTERM T2")
print("="*78)
g, pi, kp = sp.symbols('g pi kplus', positive=True)
PN = g**4/(8*pi**4)*Nc/kp                       # /(2pi)^3 stripped from both sides
# T2 = (g^2/8pi^2) * (Pgg/2) * dN_LO ; want it to equal PN * Chat
dNLO_needed = sp.simplify(PN*Chat/((g**2/(8*pi**2))*Pgg2))
print("  T2 matches the rest iff  dN_LO/d^3k = [1/(2pi)^3] *", dNLO_needed, "* Int_{x,x',y,y'} X U")
print("  i.e.  dN_LO/d^3k = 1/(2pi)^3 * (g^2/pi^2) * (1/k+) * Int X U .  Check this against your LO row.")

print()
print("="*78)
print("3.  CAN THE POLE BE SPLIT?  (advisor: absorb the Pgg half, keep the delta half)")
print("="*78)
Lk = sp.Symbol('logkLam')            # log(k+/Lambda)
d1 = sp.Symbol('delta1')             # delta(1-xi)
full_row = (d1*Lk + Chat)*(epsb - Lr)*XX + x*xb*(XX - dp*RR)
# the same row written out the long way (delta sector + non-delta sector)
long_way = d1*Lk*(M-1)*XX + derived   # the delta term inherits M-1: the '+1' of the 1/p+ block
                                      # already cancelled against the -(1/p+) X block
print("  compact - long  =", sp.simplify(sp.expand(full_row - long_way)))
print("  (the delta term carries 1/epsbar - L_r, no +1: advisor's T1 has exactly that. Good.)")
print("  -> ONE bracket [1/epsbar - L_r] multiplies  delta(1-xi) log(k+/Lam) + Pgg(xi)/(2N_c).")
print("     The coefficient of 1/epsbar is  d1*Lk + Chat  =", sp.collect(sp.expand(sp.diff(full_row, epsb)), XX))
print("     An MSbar FF counterterm is proportional to Pgg(xi) ALONE - it has no delta(1-xi)log(k+/Lam).")
print("     So subtracting one cannot remove this pole; it removes only part of it.")

print()
print("="*78)
print("4.  WHERE THE POLE ACTUALLY COMES FROM: short distance vs long distance")
print("="*78)
rho = sp.Symbol('rho', positive=True)
# trace of G near z -> y : (y-z).(y'-z) / [(y-z)^2 (y'-z)^2] with |y-z| = rho -> 0
print("  z -> y :  numerator ~ rho*|r| , denominator ~ rho^2 r^2 , measure d^2z ~ rho drho")
print("            integrand*measure ~ drho / |r|          -> CONVERGENT, no collinear pole.")
print("  z -> oo:  integrand ~ 1/z^2 , measure ~ z dz       -> LOG DIVERGENT.")
print("  Cutoff evaluation gave  Tr G = pi log(R^2/r^2) , R the LARGE-distance cutoff.")
print("  => the single 1/epsbar is an infrared (large transverse distance) divergence,")
print("     NOT a collinear one.  A fragmentation function cannot absorb it.")
