"""Row V 2b: can the 2a-style plus prescription be used? Check the xi -> 0 residue (delta(p+) coefficient)
of (i) the full 2b bracket, (ii) full minus the user's UV piece (with 1/(xibar D)), (iii) full minus the
coordinate counterterm C_UV W s^2/(r^2 (r^2+s^2)).  k+ = 1, overall e^{-ik.(w'-x)} stripped.
Residue at xi = 0:  full -> e^{ik.r} (-B4)/(y-z)^2 ;  user UV piece -> e^{ik.r} W s^2/(r^2 (y-z)^2) ;
coordinate counterterm -> W s^2/(r^2 (r^2+s^2)).
"""
import numpy as np
sq = lambda v: float(v @ v)
x = np.zeros(2); y = np.array([1.1, 0.4]); P = np.array([-0.5, 0.9]) - np.array([0.7, 1.3]); k = np.array([0.8, -0.3])
s = y - x; W = (P@s)/(sq(P)*sq(s))
def B4(z):
    r = z - x; Pv = P/sq(P); K = r/sq(r)
    return Pv@(np.outer(y-z, x-z)/sq(x-z) + np.outer(y-z, y-x)/(2*sq(y-x)))@K
full = lambda z: np.exp(1j*(k@(z-x)))*(-B4(z))/sq(y-z)
user = lambda z: full(z) - np.exp(1j*(k@(z-x)))*W*sq(s)/(sq(z-x)*sq(y-z))
coord = lambda z: full(z) - W*sq(s)/(sq(z-x)*(sq(z-x)+sq(s)))
ph = np.linspace(0, 2*np.pi, 720, endpoint=False)
circ = lambda c, rho: [c + rho*np.array([np.cos(a), np.sin(a)]) for a in ph]
print("rho^2 x angular average (a nonzero constant = log divergence of the z integral there)")
for name, fn in (("(i)   full 2b", full), ("(ii)  full - UV piece with 1/D", user), ("(iii) full - coordinate counterterm", coord)):
    zx = [abs(r*r*np.mean([fn(z) for z in circ(x, r)])) for r in (1e-2, 1e-3, 1e-4)]
    zy = [abs(r*r*np.mean([fn(z) for z in circ(y, r)])) for r in (1e-2, 1e-3, 1e-4)]
    print(f"  {name:38s} z->x: " + " ".join(f"{v:.2e}" for v in zx) + "   z->y: " + " ".join(f"{v:.2e}" for v in zy))
print(f"  predicted constants: z->x in (i): |W| = {abs(W):.4f} ;  z->y in (ii): |W| = {abs(W):.4f}")
