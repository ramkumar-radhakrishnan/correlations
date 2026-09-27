"""Row V, part b: the xi -> 0 residue and the soft-gluon region at large |x|.

Notation of the reference note:  r = z-x , s = y-x , P = x'-w' , xi = p+/k+ , q = xibar k .
Integrand (colour stripped):  T(r)/(xibar D) ,  D = (r-s)^2 + (xi/xibar) s^2 .
The xi -> 0 pole sits in T4 = -(xibar/xi) (P^i r^m/P^2 r^2)[ (s-r)^i(-r)^m/r^2 + s^m (s-r)^i/(2 s^2) ].
"""
import numpy as np
def sq(v): return float(v@v)

def T4_residue(x, y, z, xp, wp):
    """lim_{xi->0} xi * [T4 block / (xibar D)]  at fixed x, y, z."""
    P = xp-wp; r = z-x; s = y-x
    br = np.outer(s-r, -r)/sq(r) + np.outer(s-r, s)/(2*sq(s))       # [i,m]
    val = -float(P @ br @ r)/(sq(P)*sq(r))                          # the -(xibar/xi) -> -1/xi
    D0 = sq(r-s)                                                      # D at xi = 0
    return val/D0

y0 = np.array([-0.7, 1.1]); z0 = np.array([0.9, 0.35])
xp0 = np.array([1.6, 0.4]); wp0 = np.array([-0.5, -1.3])
P = xp0-wp0; u = y0-z0
pred = +(P@u)/(2*sq(P)*sq(u))        # (P.(y-z))/(2 P^2 (y-z)^2)
print("="*78)
print("THE xi -> 0 RESIDUE AT LARGE |x|  (y, z fixed; phase at xi = 0 is e^{ik.z}, x-FREE)")
print("="*78)
print("   |x|^2 x < residue >_xhat   ->   predicted  +(P.(y-z)) / [2 P^2 (y-z)^2]")
m = 4096; ph = (np.arange(m)+0.5)*(2*np.pi/m)
for X in (1e1, 1e2, 1e3, 1e4):
    v = np.mean([X**2*T4_residue(X*np.array([np.cos(a),np.sin(a)]), y0, z0, xp0, wp0) for a in ph])
    print("   |x| = %-8.0e  %16.10f      predicted %16.10f" % (X, v, pred))
print()
print("   -> a genuine 1/x^2 tail.  At xi = 0 the total phase exp[-ik.(w' - xibar z - xi x)]")
print("      has no x in it, so Int d^2x of this residue diverges logarithmically at large |x|:")
print("      the soft gluon (p+ -> 0) sitting far from the others.  It is the mirror image,")
print("      under (x <-> z, xi <-> xibar), of the xi -> 1 log that shows up inside the r integral.")
print()
print("   The colour factor there does not switch it off:")
print("      part b :  N_c [U(y) - U(x)]  ->  N_c [U(y) - 1]      as |x| -> oo")
print("      part a :  -f^{adc} f^{fbe} U^{db}(x)[U^{ce}(z) - U^{ce}(x)]")
print("                ->  -f^{adc} f^{fde} [U^{ce}(z) - delta^{ce}]  as |x| -> oo")
