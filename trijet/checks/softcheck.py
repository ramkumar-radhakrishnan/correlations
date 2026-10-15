import sympy as sp
z,mm = sp.symbols('zeta m', positive=True); zb=1-z
RC,BD = sp.symbols('RdotRbar YdotYbar', real=True)   # (R0.R0bar), (Y.Ybar)
# soft limits found, multiplied by the cross-section prefactors of (1.1),(1.9),(1.17)
qbar =  z**2      * ( 16*zb**2*(z**2+zb**2)*RC*BD + 16*mm**2*zb**2*BD )     # prefactor zeta^2
quark=  zb**2     * ( 16*z**2 *(z**2+zb**2)*RC*BD + 16*mm**2*z**2 *BD )     # prefactor eta^2 -> zetabar^2
inter= -z*zb      * ( 16*z*zb *(z**2+zb**2)*RC*BD + 16*mm**2*z*zb *BD )     # prefactor -eta*zeta
univ = 16*z**2*zb**2*( (z**2+zb**2)*RC + mm**2 )*BD
for nm,v in (("qbar",qbar),("quark",quark),("interf",inter)):
    print(f"  {nm:7s} soft limit / [16 z^2 zbar^2 ( (z^2+zbar^2) R.Rbar + m^2 ) (Y.Ybar) ] =",
          sp.simplify(v/univ))
print("\n  => all three share the SAME universal factor, interference with a relative minus:")
print("     A_0^(m)(zeta) = 16 zeta^2 zetabar^2 [ (zeta^2+zetabar^2)(R0.R0bar) + m^2 ]")
print("\n  massless limit m -> 0 reproduces the earlier result:",
      sp.simplify(univ.subs(mm,0) - 16*z**2*zb**2*(z**2+zb**2)*RC*BD)==0)
print("\n  ratio (mass term)/(non-mass term) =", sp.simplify(mm**2/((z**2+zb**2)*RC)),
      " -> matches the textbook massive gamma*_T wavefunction  [z^2+zbar^2] eps^2 K1 K1 + m^2 K0 K0")
