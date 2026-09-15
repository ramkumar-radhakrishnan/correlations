"""Check the user's three-logarithm regrouping against the exact p+ integrals."""
import numpy as np
from scipy.integrate import quad

def exact_quad(i, j, x, y, z, w, kp, P, Lam):
    """direct p+ quadrature of the full bracket / (p+ A + k+ B), times k+ (k+ factored out front)"""
    A=(x-z)@(x-z); B=(x-w)@(x-w); zw=z-w; zw2=zw@zw; xz,xw=x-z,x-w
    d = 1.0 if i==j else 0.0
    def f(p):
        br = ( d*(A-B)/(2*(p+kp)*zw2)
             + (xw[i]*zw[j]/zw2 + xz[j]*xw[i]/(2*A))/p
             + (xz[j]*zw[i]/zw2 - xz[j]*xw[i]/(2*B))/kp )
        return br/(p*A + kp*B)
    ed = np.geomspace(Lam, P, 500)
    return kp*sum(quad(f,a,b,limit=120)[0] for a,b in zip(ed[:-1],ed[1:]))

def regrouped(i, j, x, y, z, w, kp, P, Lam):
    """the user's form: three logarithms"""
    A=(x-z)@(x-z); B=(x-w)@(x-w); zw=z-w; zw2=zw@zw; xz,xw=x-z,x-w
    d = 1.0 if i==j else 0.0
    D  = P*A + kp*B                 # (vee - k+)(x-z)^2 + k+ (x-w)^2
    DL = kp*B + Lam*A               # k+ (x-w)^2 + Lambda (x-z)^2
    vee = P + kp
    # log(D/DL) bracket -- the user's five terms
    br1 = ( d/(2*zw2)
          + xz[j]*zw[i]/(zw2*A)
          - xz[j]*xw[i]/(2*B*A)
          - xw[i]*zw[j]/(zw2*B)
          - xz[j]*xw[i]/(2*A*B) )
    # log(vee/(k+ + Lambda)) : coefficient -d^{ij}/(2 (z-w)^2)  -> the user writes it with +i g^4/16pi^5
    br2 = -d/(2*zw2)
    # log(P/Lambda) bracket
    br3 = ( xw[i]*zw[j]/(zw2*B) + xz[j]*xw[i]/(2*A*B) )
    return br1*np.log(D/DL) + br2*np.log(vee/(kp+Lam)) + br3*np.log(P/Lam)

rng = np.random.default_rng(31)
kp, P, Lam = 1.0, 37.0, 1e-5
print("USER'S THREE-LOG REGROUPING  vs  direct p+ quadrature")
print("   (k+ = 1, vee - k+ = 37, Lambda = 1e-5; the overall 1/k+ and -i g^4/8pi^5 factored out)")
print(f"   {'trial':>5} {'ij':>4} {'regrouped':>18} {'quadrature':>18} {'rel':>10}")
bad = 0
for t in range(5):
    x,y,z,w = rng.normal(size=(4,2))
    for (i,j) in [(0,0),(0,1),(1,0),(1,1)]:
        r = regrouped(i,j,x,y,z,w,kp,P,Lam)
        q = exact_quad(i,j,x,y,z,w,kp,P,Lam)
        rel = abs(r-q)/max(abs(q),1e-12)
        if rel > 1e-5: bad += 1
        print(f"   {t:5d} {i}{j:>3} {r:18.9f} {q:18.9f} {rel:10.2e}")
print(f"   mismatches above 1e-5: {bad}")
print()
print("   NOTE the regrouping is EXACT in Lambda -- it uses")
print("      I1 = [ log(vee/(k+ + Lam)) - log(D/D_Lam) ] / [k+(B-A)]")
print("      I2 = [ log((vee-k+)/Lam) - log(D/D_Lam) ] / (k+ B)")
print("      I3 = log(D/D_Lam) / A")
print("   with D = (vee-k+)(x-z)^2 + k+(x-w)^2 ,  D_Lam = k+(x-w)^2 + Lambda (x-z)^2 .")
