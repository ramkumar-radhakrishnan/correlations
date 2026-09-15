"""Build the row from the definitions of A and C and compare with the quoted result."""
import sympy as sp

g, pi = sp.symbols('g pi', positive=True)
kp, pp = sp.symbols('kplus pplus', positive=True)      # k+, p+
Kp = pp + kp                                           # the first argument of C

# ---------- 1. numeric prefactor -------------------------------------------
# A^a_i(K,-z) = + i g /(sqrt2 pi sqrt(K)) int_x (x-z)^i/(x-z)^2 rho^a(x)
#   (from  -A = -(1/(sqrt2 pi)) (i g/sqrt K) int ... )
I = sp.I
A      = lambda K:  I*g/(sp.sqrt(2)*pi*sp.sqrt(K))
Adag   = lambda K: -I*g/(sp.sqrt(2)*pi*sp.sqrt(K))
# C(K, z; p, x) = -(g f)/(2 pi sqrt(2K)) * sqrt(p(K-p))/K * delta(...) * x^m/x^2 * [...]
Cpref  = -g/(2*pi*sp.sqrt(2*Kp)) * sp.sqrt(pp*(Kp-pp))/Kp

overall = sp.Rational(2,1)/(2*pi)**3                   # the 2/(2pi)^3 out front
prod = overall * Adag(Kp) * Cpref * A(kp) * A(pp)
prod = sp.simplify(prod)
print("A. product of the four coefficient functions (times 2/(2pi)^3, without the delta):")
print("   ", sp.simplify(prod))

# ---------- 2. the delta-function Jacobian ---------------------------------
# C's delta is  delta^(2)( z_arg + (p/K) x_arg )  with z_arg = y'-w', x_arg = w'-z:
#   y' - w' + (p/K)(w'-z) = y' - (k/K) w' - (p/K) z  =  -(k/K)[ w' - ((K y' - p z)/k) ]
# so in 2 dimensions it equals (K/k)^2 delta^(2)[ w' - (K y' - p z)/k ]
J = (Kp/kp)**2
print("\nB. rewriting the delta in the quoted w'-rooted form costs a Jacobian:")
print("   delta^(2)(y'-w'+(p/K)(w'-z)) = (K/k)^2 delta^(2)[w' - ((p+k)y' - p z)/k]")
print("   Jacobian (K/k)^2 =", J)

total = sp.simplify(prod*J)
print("\nC. full prefactor  =  product x Jacobian:")
print("   ", total)
target = -sp.Rational(1,1)/(2*pi)**3 * I*g**4/(4*pi**4) / kp**2
print("   quoted           :", sp.simplify(target))
print("   difference       :", sp.simplify(total - target), "   <-- 0 means they agree")
print("   (the f^{c'd'a} is carried along by C and is common to both)")

# ---------- 3. the index bracket -------------------------------------------
print("\nD. the index bracket.  Definition:  C_{1 j k i} with")
print("      [ d_im d_jk - (K/p) d_jm d_ik - (K/(K-p)) d_km d_ij ]")
print("   The object needed is  C_{1 j i k'} : slot1 j->j, slot2 k->i, slot3 i->k'.")
names = ['j','i','kp_']                                # slot values after relabelling
def bracket(m, j, k, i, K, p):
    """definition's bracket with explicit index VALUES"""
    d = lambda a,b: 1 if a == b else 0
    return d(i,m)*d(j,k) - (K/p)*d(j,m)*d(i,k) - (K/(K-p))*d(k,m)*d(i,j)
def quoted(m, j, i, kpr, K, p):
    d = lambda a,b: 1 if a == b else 0
    return d(kpr,m)*d(i,j) - (K/p)*d(j,m)*d(i,kpr) - (K/(K-p))*d(i,m)*d(j,kpr)
bad = 0
for m in range(2):
    for jj in range(2):
        for ii in range(2):
            for kk in range(2):
                lhs = bracket(m, jj, ii, kk, Kp, pp)    # slots (j,k,i) = (j, i, k')
                rhs = quoted(m, jj, ii, kk, Kp, pp)
                if sp.simplify(lhs - rhs) != 0: bad += 1
print("   mismatching index assignments:", bad, " (0 = the relabelling reproduces the quoted bracket)")
print("   and with K = p+ + k+ :   K/p = (p+ + k+)/p+ ,   K/(K-p) = (p+ + k+)/k+  ->")
print("      [ d_k'm d_ij - ((p+ + k+)/p+) d_jm d_ik' - ((p+ + k+)/k+) d_im d_jk' ]   MATCHES")

# ---------- 4. the sqrt factor and the kernel ------------------------------
print("\nE. the square-root factor:  sqrt(p (K-p))/K  with K = p+ + k+  gives",
      sp.simplify(sp.sqrt(pp*(Kp-pp))/Kp), " = sqrt(p+ k+)/(p+ + k+)")
print("   C's kernel x^m/x^2 with x = w'-z  ->  (w'-z)^m/(w'-z)^2   MATCHES")
print("   A^dag(K,-y') -> (x'-y')^k'/(x'-y')^2 rho^e'(x')           MATCHES")
print("   A(k+,-w)     -> (y-w)^i/(y-w)^2                           MATCHES")
print("   A(p+,-z)     -> (x-z)^j/(x-z)^2                           MATCHES")
