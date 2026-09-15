"""The p+ integral of row 4, in closed form, checked against quadrature."""
import numpy as np
from scipy.integrate import quad

def Tij_integrated(i, j, x, y, z, w, kp, V, Lam):
    """the p+-integrated bracket, closed form (overall 1/k+ factored out)"""
    A = (x-z)@(x-z); B = (x-w)@(x-w); zw = z-w; zw2 = zw@zw
    xz, xw = x-z, x-w
    D = V*A + kp*B                       # V here is (vee - k+)
    vee = V + kp
    d = 1.0 if i == j else 0.0
    t1 = -d/(2*zw2)*np.log(vee*B/D)
    t2 = xw[i]/B*( zw[j]/zw2 + xz[j]/(2*A) )*( np.log(V/Lam) + np.log(kp*B/D) )
    t3 = xz[j]/A*( zw[i]/zw2 - xw[i]/(2*B) )*np.log(D/(kp*B))
    return t1 + t2 + t3

def Tij_quad(i, j, x, y, z, w, kp, V, Lam):
    """direct p+ quadrature of  [bracket]/(p+ A + k+ B),  k+ factored out"""
    A = (x-z)@(x-z); B = (x-w)@(x-w); zw = z-w; zw2 = zw@zw
    xz, xw = x-z, x-w
    d = 1.0 if i == j else 0.0
    def f(p):
        br = ( d*(A-B)/(2*(p+kp)*zw2)
             + (xw[i]*zw[j]/zw2 + xz[j]*xw[i]/(2*A))/p
             + (xz[j]*zw[i]/zw2 - xz[j]*xw[i]/(2*B))/kp )
        return br/(p*A + kp*B)
    ed = np.geomspace(Lam, V, 400)
    return kp*sum(quad(f, a, b, limit=100)[0] for a, b in zip(ed[:-1], ed[1:]))

rng = np.random.default_rng(21)
print("p+-INTEGRATED BRACKET: closed form vs direct quadrature   (k+ = 1, vee-k+ = 37, Lambda = 1e-6)")
kp, V, Lam = 1.0, 37.0, 1e-6
print(f"   {'trial':>5} {'ij':>4} {'closed form':>18} {'quadrature':>18} {'rel':>10}")
for t in range(4):
    x, y, z, w = rng.normal(size=(4,2))
    for (i,j) in [(0,0),(0,1),(1,1)]:
        c = Tij_integrated(i,j,x,y,z,w,kp,V,Lam)
        q = Tij_quad(i,j,x,y,z,w,kp,V,Lam)
        print(f"   {t:5d} {i}{j:>3} {c:18.9f} {q:18.9f} {abs(c-q)/abs(q):10.2e}")
print()
print("the three pieces, with A = (x-z)^2, B = (x-w)^2, D = (vee-k+)A + k+ B :")
print("   I1 = int dp+/[(p+ + k+)(p+A+k+B)] = log(vee B / D)/[k+(B-A)]")
print("      x its coefficient d^ij (A-B)/(2(z-w)^2)  ->  -d^ij log(vee B/D) / (2 k+ (z-w)^2)")
print("      (the (A-B) cancels the (B-A): the A = B singularity is removable and GONE)")
print("   I2 = int dp+/[p+(p+A+k+B)] = [ log((vee-k+)/Lambda) + log(k+B/D) ]/(k+ B)")
print("      the 1/B combines with the (x-w)^i in the numerator into a WW kernel")
print("   I3 = int dp+/[p+A+k+B] = log(D/(k+B))/A")
print("      the 1/A combines with the (x-z)^j into a WW kernel")
