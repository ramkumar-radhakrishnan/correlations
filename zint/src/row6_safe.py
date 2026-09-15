"""The manifestly-finite regrouping: put the delta^{ij} piece on ONE combined log."""
import numpy as np
V, Lam, kp = 37.0, 1e-6, 1.0; vee = V + kp

def T_quoted(i,j,x,z,w,Vv,yz):
    A=(x-z)@(x-z); B=(x-w)@(x-w); zw=z-w; zw2=zw@zw; xz,xw=x-z,x-w
    d = 1.0 if i==j else 0.0
    D, DL = V*A+kp*B, Lam*A+kp*B
    br1 = ( xz[j]*zw[i]/(zw2*A) - xz[j]*xw[i]/(2*B*A)
          - xw[i]*zw[j]/(zw2*B) - xz[j]*xw[i]/(2*A*B) + d/(2*zw2) )
    return ( br1*np.log(D/DL) - (d/(2*zw2))*np.log(vee/(Lam+kp))
           + (xw[i]*zw[j]/(zw2*B) + xz[j]*xw[i]/(2*A*B))*np.log(V/Lam) )

def T_safe(i,j,x,z,w,Vv,yz):
    A=(x-z)@(x-z); B=(x-w)@(x-w); zw=z-w; zw2=zw@zw; xz,xw=x-z,x-w
    d = 1.0 if i==j else 0.0
    D, DL = V*A+kp*B, Lam*A+kp*B
    br_reg = ( xz[j]*zw[i]/(zw2*A) - xw[i]*zw[j]/(zw2*B) - xz[j]*xw[i]/(A*B) )   # 3 terms, no delta
    return ( br_reg*np.log(D/DL)
           + (d/(2*zw2))*np.log(D*(Lam+kp)/(DL*vee))        # <- combined log, vanishes at z=w
           + (xw[i]*zw[j]/(zw2*B) + xz[j]*xw[i]/(2*A*B))*np.log(V/Lam) )

rng = np.random.default_rng(77)
print("1. the safe regrouping is IDENTICAL to the quoted one")
m = 0.0
for _ in range(4):
    x,y,z,w,Vv = rng.normal(size=(5,2))
    for i in range(2):
        for j in range(2):
            m = max(m, abs(T_quoted(i,j,x,z,w,Vv,y-z) - T_safe(i,j,x,z,w,Vv,y-z)))
print(f"   max |quoted - safe| over 4 geometries x 4 index pairs = {m:.2e}")
print()
print("2. and it is manifestly finite at z -> w, term by term")
x0=np.array([0.3,-0.2]); y0=np.array([-0.6,0.9]); w0=np.array([-0.4,-0.7]); Vv=np.array([0.4,0.7])
th=np.linspace(0,2*np.pi,8192,endpoint=False); u=np.stack([np.cos(th),np.sin(th)],-1)
def scanned(f, whichdelta):
    out=[]
    for rho in (1e-2,1e-4,1e-6):
        vals=[]
        for uu in u:
            z = w0 + uu*rho
            s=0.0
            for i in range(2):
                for j in range(2):
                    yz=y0-z; K=Vv[i]*yz[j]/(yz@yz)
                    s += K*f(i,j,x0,z,w0,Vv,yz)
            vals.append(s)
        out.append(np.mean(vals)*rho**2)
    return out
print("   quoted form, delta-piece split across two logs:")
print("      rho^2<F> =", "  ".join(f"{v:+.4e}" for v in scanned(T_quoted, True)))
print("   safe form, delta-piece on one combined log:")
print("      rho^2<F> =", "  ".join(f"{v:+.4e}" for v in scanned(T_safe, True)))
print("   both -> 0 (the totals agree); the difference is that in the safe form")
print("   NO INDIVIDUAL TERM diverges, because log[D(Lam+k+)/(DL vee)] -> 0 at z = w.")
print()
print("3. check that combined log vanishes at z = w:")
for rho in (1e-1,1e-2,1e-3,1e-4):
    z = w0 + np.array([rho,0.0])
    A=(x0-z)@(x0-z); B=(x0-w0)@(x0-w0); D, DL = V*A+kp*B, Lam*A+kp*B
    print(f"   |z-w|={rho:.0e}:  log[D(Lam+k+)/(DL vee)] = {np.log(D*(Lam+kp)/(DL*vee)):+.6e}")
