"""Group II Row IV: derive (0.1) from the quoted C_1 and A^(1) definitions."""
import sympy as sp, numpy as np, itertools
g,pi=sp.symbols('g pi',positive=True); P,K=sp.symbols('pplus kplus',positive=True); I=sp.I
S=P+K                                        # the parent + momentum in this row

print("="*78); print("1.  THE OVERALL PREFACTOR OF (0.1), FROM YOUR DEFINITIONS"); print("="*78)
# A^(1)_i(k+_arg, .) = + i g /(sqrt2 pi sqrt(k+_arg))  x  int_x (x-z)^i/(x-z)^2 rho
A  =  I*g/(sp.sqrt(2)*pi*sp.sqrt(S))
Ad = -I*g/(sp.sqrt(2)*pi*sp.sqrt(S))
# C_1(k+_arg = S ; p+_arg = p+) :  -g f /(2 pi sqrt(2 S)) * sqrt(p+ (S-p+))/S ,  S-p+ = k+
C  = -g/(2*pi*sp.sqrt(2*S))*sp.sqrt(P*(S-P))/S
Cd = -g/(2*pi*sp.sqrt(2*S))*sp.sqrt(P*(S-P))/S
print("   A^(1)  at k+_arg = p+ + k+ :  %s"%sp.simplify(A))
print("   Cbar   at k+_arg = p+ + k+ :  %s     (sqrt(p+(S-p+))/S = sqrt(p+ k+)/S)"%sp.simplify(C))
one=sp.simplify(A*C*Ad*Cd)
print("   |A C|^2 = %s"%sp.simplify(sp.powsimp(one,force=True)))
tot=sp.simplify(sp.Integer(4)/(2*pi)**3*one)
quoted=1/(2*pi)**3*g**4/(4*pi**4)*P*K/S**4
print("   (4/(2pi)^3) x that   = %s"%sp.simplify(sp.powsimp(tot,force=True)))
print("   eq (0.1) quotes      = %s"%sp.simplify(quoted))
print("   ratio quoted/derived = %s      <-- 1 means EXACT (magnitude and sign)"
      %sp.simplify(sp.powsimp(quoted/tot,force=True)))
print("   phases: amplitude (i)(-1) = -i ,  conjugate (-i)(-1) = +i ,  product = +1  -> overall + sign OK")

print(); print("="*78); print("2.  THE TENSOR BRACKET:  C_1jki  CALLED AS  C_1jik"); print("="*78)
print("   your definition   C_1jki :  d_im d_jk  - (S/p+) d_jm d_ik  - (S/k+) d_km d_ij")
print("   relabel slot2<->slot3 (i<->k), i.e. call it as C_1jik :")
print("                             :  d_km d_ji  - (S/p+) d_jm d_ki  - (S/k+) d_im d_jk")
print("   eq (0.1) writes           :  d_km d_ij  - (S/p+) d_jm d_ik  - (S/k+) d_im d_kj")
print("   -> IDENTICAL (delta is symmetric).  So (0.1) is consistent with the definition,")
print("      PROVIDED the index contracted with A^(1)_k really is the definition's 3rd slot.")

print(); print("="*78); print("3.  THE DELTA AND THE WW KERNEL"); print("="*78)
print("   C_1(k+_arg, z_arg ; p+_arg, x_arg) carries  delta^(2)[ z_arg + (p+_arg/k+_arg) x_arg ]  and  x_arg^m/x_arg^2")
print("   row calls it with  k+_arg = p+ + k+ ,  z_arg = y - w ,  p+_arg = p+ ,  x_arg = w - z :")
print("     delta^(2)[ (y-w) + p+/(p+ + k+) (w-z) ]     <-- matches (0.1) exactly")
print("     (w-z)^m/(w-z)^2                             <-- matches (0.1) exactly")

print(); print("="*78); print("4.  THE THETA FUNCTION AND THE LOWER LIMIT Lambda"); print("="*78)
print("   Theta( k+_arg - p+_arg - Lambda )  with k+_arg = p+ + k+ , p+_arg = p+  gives")
print("     Theta( k+ - Lambda )   -- a condition on the MEASURED gluon, always satisfied.")
print("   It does NOT enforce p+ > Lambda, yet the row integrates from Lambda.")
print("   The 1/p+ pole of the bracket is the one that needs a cutoff, so the vertex needs")
print("   a Theta(p+ - Lambda) as well (or the row's lower limit must be justified elsewhere).")

print(); print("="*78); print("5.  (0.20):  THE d_perp IN THE DENOMINATOR IS SPURIOUS"); print("="*78)
e=sp.symbols('epsilon'); dperp=2-2*e
print("   from (0.19) the delta coefficient is  (1/d_perp)(1/eps - 1) = %s"%sp.simplify(dperp**-1*(1/e-1)))
print("   so  G^{mm'} = pi^{1-e}(r^2)^{-e} mu^{2e} N(e) [ delta/(2 eps) - r^m r^m'/r^2 ]")
print("   and expanding gives   pi delta^{mm'}/2 [ ... ]  -  pi r^m r^m'/r^2 .")
print("   Your (0.20) writes  pi delta^{mm'}/(2 d_perp) : one factor of d_perp too many.")
c=sp.Symbol('c'); L=1/e+c; XX=sp.Symbol('XX')
print()
print("   consequence for (0.22) part 1.  Contracting delta^{mm'} into the bracket gives a")
print("   COMMON d_perp (first term keeps its explicit one, third and fourth get one from the")
print("   trace), so:")
print("     with your (0.20):   (pi L / (2 d_perp)) * d_perp [..] = (pi/2) L [..]        -> g^4/(8 pi^3) * L")
print("     correct         :   (pi L / 2)          * d_perp [..] = pi (L-1) [..]        -> g^4/(4 pi^3) * (L-1)")
print("   ratio correct/yours = %s"%sp.simplify(sp.pi*(L-1)/(sp.pi*L/2)))
