"""Group II Row IV: is the 1/epsbar a running-coupling (beta_0) pole?

Three tests.
  (1) Order counting: what an O(g^4) running-coupling term must look like.
  (2) The beta_0 coefficient test: where 11N_c/6 comes from, and what this row
      puts in its place.
  (3) The n_f test: beta_0 needs a quark loop, which this row does not have.
"""
import sympy as sp

x, Nc, nf, TR = sp.symbols('xi N_c n_f T_R', positive=True)
xb = 1 - x

print("="*78)
print("1.  WHERE 11 N_c / 6 ACTUALLY COMES FROM  (momentum sum rule)")
print("="*78)
# regulated splitting function: the plus distribution acting on a test function f
def plus_int(f):
    """int_0^1 dxi f(xi)/[1-xi]_+ = int_0^1 dxi (f(xi)-f(1))/(1-xi)"""
    return sp.integrate(sp.simplify((f - f.subs(x, 1))/(1 - x)), (x, 0, 1))

I1 = plus_int(x*x)                      # int z * z/[1-z]_+
I2 = sp.integrate(x*xb/x, (x, 0, 1))    # int z * (1-z)/z
I3 = sp.integrate(x*x*xb, (x, 0, 1))    # int z * z(1-z)
print("  int_0^1 dz z * z/[1-z]_+      =", I1)
print("  int_0^1 dz z * (1-z)/z        =", I2)
print("  int_0^1 dz z * z(1-z)         =", I3)
tot = sp.simplify(2*Nc*(I1 + I2 + I3))
print("  sum x 2N_c                    =", tot)
a = sp.symbols('a')
sol = sp.solve(sp.Eq(tot + a + 2*nf*sp.integrate(x*TR*(x**2+xb**2), (x, 0, 1)), 0), a)[0]
print("  momentum sum rule  int z Pgg + 2 n_f int z Pqg = 0  =>  delta(1-z) coefficient a =",
      sp.simplify(sol))
print("  with T_R = 1/2 :", sp.simplify(sol.subs(TR, sp.Rational(1,2))), " = 11N_c/6 - n_f/3 = beta_0/2 .")
print()
print("  NOTE the -3/2 above: the whole 11/6 hangs off the SOFT endpoint xi -> 1")
print("  of xi/[1-xi]_+ .  Whatever cuts off that endpoint decides what sits on delta(1-xi).")

print()
print("="*78)
print("2.  WHAT THIS ROW PUTS ON delta(1-xi) INSTEAD")
print("="*78)
print("  In Row IV the p+ integral runs Lambda < p+ < V - k+ , i.e. the soft endpoint is cut")
print("  in RAPIDITY, not by transverse kinematics.  The plus prescription therefore gives")
print("      int_Lambda dp+/p+  ->  delta(1-xi) log(k+/Lambda) + [ ]_+ ,")
print("  so the delta(1-xi) coefficient is  log(k+/Lambda)  -- a rapidity logarithm --")
print("  and NOT the constant 11N_c/6 .  No beta_0 can be read off this row.")

print()
print("="*78)
print("3.  ORDER COUNTING:  what an O(g^4) running-coupling term must look like")
print("="*78)
g, mu, r, kp, Lam = sp.symbols('g mu r kplus Lambda', positive=True)
print("  LO   ~ g^2 * (LO structure)")
print("  NLO  ~ g^4 * (...)                    <- this row")
print("  g^2(mu0) = g^2(mu) [ 1 + (alpha_s beta_0/4pi) log(mu^2/mu0^2) + ... ]")
print("  => absorbing a log(mu^2 r^2) into the LO coupling produces, at O(g^4),")
print("        g^4 * beta_0 * log(mu^2 r^2) * (LO structure)    with NO rapidity log,")
print("        NO Pgg(xi), and a pure delta(1-xi).")
print("  Row IV's log(4 pi^2 mu^2 r^2) instead comes multiplied by")
print("        delta(1-xi) log(k+/Lambda) + Pgg(xi)/(2N_c) ,")
print("  neither piece of which has that shape:")
print("     - the delta piece carries an extra log(k+/Lambda)  -> that is evolution, and")
print("       making ITS coupling run is an O(g^6) effect, one order beyond this calculation;")
print("     - the Pgg piece is xi-dependent, so it is a splitting kernel, not a beta function.")

print()
print("="*78)
print("4.  THE n_f TEST")
print("="*78)
print("  beta_0 = 11N_c/3 - 2n_f/3 .  A pure-gluon row cannot generate the -2n_f/3 .")
print("  The n_f piece lives in the g -> q qbar row, whose transverse integral is UV")
print("  divergent at short distance.  Row IV has no quark loop and, by the z -> y power")
print("  counting, no short-distance divergence at all:")
print("     z -> y :  u^m (u+r)^{m'} / [u^2 (u+r)^2] ~ (1/rho)(1/|r|) , d^2 z = rho drho dphi")
print("               -> int drho / |r|   CONVERGENT.")
print("  A coupling is renormalised by a UV pole.  This pole is at |z| -> infinity.")
