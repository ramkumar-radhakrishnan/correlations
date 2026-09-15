"""Row 3 (gluon B_2).  Is there a transverse UV divergence?
Integrand (constants dropped, k+ = 1, xi = p+, xib = 1-p+):

  F = [xi xib / (xi (y-x)^2 + xib (y-z)^2)] * (x'-w')^i/(x'-w')^2 * (z-x)^m/(z-x)^2 * T^{im}
"""
import numpy as np
dperp = 2.0

def T(i, m, x, y, z, xi):
    xb = 1-xi
    yx, yz, xz = y-x, y-z, x-z
    yx2, yz2, xz2 = yx@yx, yz@yz, xz@xz
    d = 1.0 if i == m else 0.0
    t1 = d*dperp/(2*xz2)*(yx2 - yz2)
    t2 = -1/xb**2*( yx[i]*xz[m]/xz2 - yx[i]*yz[m]/(2*yz2) )
    t3 = d/xi*( (yz@xz)/xz2 + (yx@yz)/(2*yx2) - (yx2-yz2)/(2*xz2) )
    t4 = -1/xi**2*( yz[i]*xz[m]/xz2 + yx[m]*yz[i]/(2*yx2) )
    t5 = d/xb*( (yx@xz)/xz2 - (yx@yz)/(2*yz2) - (yx2-yz2)/(2*xz2) )
    t6 = -1/(xi*xb)*( yx[m]*xz[i]/xz2 - yx[m]*yz[i]/(2*yz2)
                      + yz[m]*xz[i]/xz2 + yx[i]*yz[m]/(2*yx2) )
    return t1+t2+t3+t4+t5+t6, (t1,t2,t3,t4,t5,t6)

def integrand(x, y, z, V, xi):
    """V = (x'-w')^i/(x'-w')^2 , a fixed spectator vector"""
    xb = 1-xi
    den = xi*((y-x)@(y-x)) + xb*((y-z)@(y-z))
    zx = z - x; zx2 = zx@zx
    tot = 0.0; parts = np.zeros(6)
    for i in range(2):
        for m in range(2):
            t, ps = T(i, m, x, y, z, xi)
            w = V[i]*zx[m]/zx2
            tot += w*t
            parts += w*np.array(ps)
    return xi*xb/den*tot, xi*xb/den*parts

rng = np.random.default_rng(2)
x = np.array([0.3, -0.2]); y = np.array([-0.6, 0.9])
Vv = np.array([0.4, 0.7]); xi = 0.35
th = np.linspace(0, 2*np.pi, 4096, endpoint=False)
u = np.stack([np.cos(th), np.sin(th)], -1)

print("A. z -> x   (rho = |z-x|).  rho^2 x <integrand>_theta")
print("   a nonzero limit  =>  integrand ~ 1/rho^2  =>  LOGARITHMIC UV DIVERGENCE")
print(f"   {'rho':>9} {'rho^2<F>':>14}   contributions of the six bracket structures")
for rho in (1e-1, 1e-2, 1e-3, 1e-4, 1e-5):
    vals = []; pars = []
    for uu in u:
        v, p = integrand(x, y, x + uu*rho, Vv, xi)
        vals.append(v); pars.append(p)
    m = np.mean(vals)*rho**2
    mp = np.mean(pars, axis=0)*rho**2
    print(f"   {rho:9.0e} {m:+14.8f}   " + " ".join(f"{q:+8.4f}" for q in mp))
print("   (columns: t1 d_perp, t2 1/xib^2, t3 1/xi, t4 1/xi^2, t5 1/xib, t6 1/(xi xib))")
