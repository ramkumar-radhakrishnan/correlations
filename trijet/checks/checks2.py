import numpy as np
from scipy import integrate

# ---- (5) the lifetime "triangle" transverse integral -------------------------
# J = int d^2z  r_bb'^2/(r_zb^2 r_zb'^2) * ln(min(r_zb^2,r_zb'^2) Qf^2) * Theta(min*Qf^2 - 1)
# claim (leading, r_bb' Qf >> 1):  J -> pi * ln^2(r_bb'^2 Qf^2)
# set Qf = 1, b = (0,0), b' = (r,0)
def J_num(r):
    def integrand(u, phi):          # rho = e^u  ->  rho drho = e^{2u} du
        rho = np.exp(u)
        zx, zy = rho*np.cos(phi), rho*np.sin(phi)
        d2 = (zx-r)**2 + zy**2      # r_zb'^2
        m  = min(rho**2, d2)
        if m <= 1.0: return 0.0
        return (r**2/(rho**2*d2))*np.log(m)*rho**2   # *rho^2 from measure
    val,_ = integrate.dblquad(integrand, 0.0, 2*np.pi,
                              lambda p: 0.0, lambda p: np.log(60*r),
                              epsabs=1e-9, epsrel=1e-9)
    return val

print("=== (5) triangle integral:  J(r)  vs  pi ln^2(r^2)  [Qf=1] ===")
for r in [20.0, 50.0, 150.0, 400.0]:
    num = J_num(r); ana = np.pi*np.log(r**2)**2
    print(f"  r={r:6.1f}  num={num:12.5f}  pi*ln^2(r^2)={ana:12.5f}   ratio={num/ana:.5f}"
          f"   (num-ana)/ln(r^2)={ (num-ana)/np.log(r**2):8.4f}")
