"""UV (short-distance) test, focused on y' -> z, and on the multi-point collapses.
r = y'-z.  Only the kernel (y'-z)^m/(y'-z)^2 is singular there."""
import numpy as np
from scipy.integrate import quad
rng = np.random.default_rng(3)
zz, xx, xpp, yy, ww = rng.normal(size=(5,2))
kk = np.array([0.9,0.6]); xi = 0.3
d, e = xx-zz, yy-ww

def integ(r):                       # r = y'-z ;  y' = z + r
    r2 = np.sum(r*r,-1); C = xpp-zz-r; C2 = np.sum(C*C,-1)
    ph = np.exp(-1j*xi*(r@kk))
    T1 = np.sum(r*C,-1)/(r2*C2)                 # d_k'm d_ij
    T2 = (r@d)/r2 * (C@e)/C2                    # d_jm d_ik'   <- the pole structure
    T3 = (r@e)/r2 * (C@d)/C2                    # d_im d_jk'
    return {"T1 x ph":ph*T1, "T2 x ph":ph*T2, "T3 x ph":ph*T3,
            "T2 x (ph-1)":(ph-1)*T2, "T2 bare (g(0))":T2+0j}

def disc(eps, nth=2001):
    th = np.linspace(0,2*np.pi,nth,endpoint=False)
    u = np.stack([np.cos(th),np.sin(th)],-1)
    keys = list(integ(u[:1]*1e-3).keys()); out={}
    for k in keys:
        g = lambda rho,k=k: rho*np.mean(integ(u*rho)[k])*2*np.pi
        tot=0j
        ed=np.concatenate(([0.0],np.geomspace(eps*1e-4,eps,25)))
        for p,q in zip(ed[:-1],ed[1:]):
            tot += quad(lambda x:g(x).real,p,q,limit=80)[0]+1j*quad(lambda x:g(x).imag,p,q,limit=80)[0]
        out[k]=tot
    return out

print("A. int_{|y'-z| < eps} d^2y' of each structure   (UV test at y' -> z)")
eps=[1e-1,1e-2,1e-3]; res=[disc(x) for x in eps]
print(f"   {'structure':<18}" + "".join(f"{'eps=%.0e'%x:>13}" for x in eps) + "   ratios/decade")
for k in res[0]:
    v="".join(f"{abs(r[k]):13.3e}" for r in res)
    r1=abs(res[1][k])/abs(res[0][k]); r2=abs(res[2][k])/abs(res[1][k])
    print(f"   {k:<18}{v}   {r1:.4f} {r2:.4f}")
print("   ratio ~0.01/decade  =>  F ~ eps^2 .  The leading 1/rho integrates to ZERO by")
print("   parity (int d^2r r^m/r^2 = 0), so the disc integral is even better than the")
print("   1/rho power counting suggests.  NO ultraviolet divergence at y' -> z.")
print()
print("B. the subtraction improves it further:  ph - 1 = -i xi k.r + O(r^2)")
for x in (1e-2,1e-3,1e-4):
    r=disc(x)
    print(f"   eps={x:.0e}:  |T2 x ph| {abs(r['T2 x ph']):.3e}   |T2 x (ph-1)| {abs(r['T2 x (ph-1)']):.3e}")
print()
print("C. SIMULTANEOUS COLLAPSES (all separations scaling as rho, 2d each)")
print("   points   relative coords   measure        kernels singular   integrand   net")
for nm,npts,nker in [("y'~z",2,1),("y'~z~x'",3,2),("y'~z~x",3,2),("y'~z~x'~x",4,3)]:
    nrel=npts-1; print(f"   {nm:<9} {nrel:^15} rho^{2*nrel} drho/rho {nker:^18} rho^-{nker}"
                       f"     rho^{2*nrel-nker} drho/rho  -> converges")
