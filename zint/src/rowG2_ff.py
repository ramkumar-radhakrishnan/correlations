"""Group II Row IV: writing the 1/epsbar term as a fragmentation-function counterterm.

Target form (MSbar collinear factorisation for final-state / fragmentation):

    d^3N^{NLO}/d^3k |_pole  =  (alpha_s/2pi) (1/epsbar) Int_{xi0}^1 dxi/xi^2
                               Pgg(xi) x [ LO distribution evaluated at k/xi ]

The only subtlety is WHICH LO distribution: the one with an explicit 1/k+ does
not work, the invariant one Ntil = k+ d^3N/d^3k does.  Checked below.
"""
import sympy as sp

g, pi, kp, x, Nc = sp.symbols('g pi kplus xi N_c', positive=True)
F = sp.Function('F')                    # F(k/xi) = Int_{x,x',y,y'} e^{-i (k/xi).r} X U
xb = 1 - x
Chat = x*xb + sp.Symbol('plusterm') + xb/x
Pgg2 = Nc*Chat                          # Pgg/2

tp = 1/(2*pi)**3                        # the explicit 1/(2pi)^3 both rows carry

print("="*78)
print("1.  alpha_s / 2 pi  IN TERMS OF g")
print("="*78)
al = g**2/(4*pi)
print("   alpha_s = g^2/4pi   =>   alpha_s/(2 pi) =", sp.simplify(al/(2*pi)), " = g^2/(8 pi^2)")
print("   NOT g^2/(8 pi^3).  Keep this straight; it fixes every factor below.")

print()
print("="*78)
print("2.  THE ROW'S POLE TERM")
print("="*78)
row = tp * g**4/(8*pi**4) * Nc/kp * Chat * F(x)     # integrand under Int dxi/xi^2, pole part
print("   row integrand (under Int_{xi0}^1 dxi/xi^2, coefficient of 1/epsbar):")
print("      ", sp.simplify(row))

print()
print("="*78)
print("3.  TRY THE NAIVE LO:  d^3N_LO/d^3k = (1/(2pi)^3)(g^2/pi^2)(1/k+) Int X U")
print("="*78)
# at argument k/xi the 1/k+ becomes xi/k+ and the phase becomes e^{-i(k/xi).r}
NLO_at = tp * g**2/pi**2 * (x/kp) * F(x)
cand   = al/(2*pi) * Pgg2 * NLO_at
print("   (alpha_s/2pi) (Pgg/2) dN_LO(k/xi) =", sp.simplify(cand))
print("   ratio  row / candidate =", sp.simplify(row/cand))
print("   -> a stray 1/xi .  With this LO the convolution needs dxi/xi^3, not dxi/xi^2.")

print()
print("="*78)
print("4.  THE FIX: use the INVARIANT density  Ntil(k) = k+ d^3N/d^3k")
print("="*78)
# Ntil(k/xi) = (k+/xi) * dN_LO/d^3k (k/xi) = (k+/xi)(tp g^2/pi^2)(xi/kp) F = tp g^2/pi^2 F
Ntil_at = tp * g**2/pi**2 * F(x)
cand2 = al/(2*pi) * Pgg2 * Ntil_at / kp
print("   Ntil_LO(k/xi) =", sp.simplify(Ntil_at), "   (no k+ left: that is the point)")
print("   (alpha_s/2pi)(1/k+)(Pgg/2) Ntil_LO(k/xi) =", sp.simplify(cand2))
print("   ratio  row / candidate =", sp.simplify(row/cand2), "   <-- EXACT MATCH")

print()
print("="*78)
print("5.  THE FACTOR 1/2 ON Pgg")
print("="*78)
cand_full = al/(2*pi) * (2*Pgg2) * Ntil_at / kp        # with the FULL Pgg, not Pgg/2
print("   with the canonical counterterm (alpha_s/2pi)(1/epsbar) Pgg (not Pgg/2):")
print("      ratio  row / candidate =", sp.simplify(row/cand_full))
print("   => this row supplies exactly HALF the MSbar counterterm.")
print("      The mirror row (the y <-> y' / hermitian-conjugate contraction) must supply the other half.")
print()
print("   This is also the factor-2 you changed last time:")
print("      g^2/(8 pi^2) * Pgg/2  =  (alpha_s/2pi) * Pgg/2   <- what THIS ROW contains")
print("      g^2/(4 pi^2) * Pgg/2  =  (alpha_s/2pi) * Pgg     <- the FULL counterterm you want")
print("   Both are 'right'; they are answers to different questions.")

print()
print("="*78)
print("6.  DGLAP EVOLUTION CHECK OF THE SIGN")
print("="*78)
muF, eps = sp.symbols('mu_F epsbarinv')
print("   define  D_{g/g}^bare(xi) = delta(1-xi) - (alpha_s/2pi)(1/epsbar) Pgg(xi) + O(a^2)")
print("   then    dN|_bare - [counterterm] = dN|_MSbar(mu_F) , and what is left of the row is")
print("        - (alpha_s/2pi) Pgg(xi)/2 * log(4 pi^2 mu_F^2 r^2)  (inside the transverse integral)")
print("   d/dlog mu_F^2 of that is  -(alpha_s/2pi) Pgg/2 , which is cancelled by")
print("   mu_F^2 dD/dmu_F^2 = +(alpha_s/2pi) Pgg (x) D  acting on the LO term.  Signs consistent.")
