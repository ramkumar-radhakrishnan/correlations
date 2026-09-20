import numpy as np, itertools
from scipy import integrate

def dot(a,b): return a[0]*b[0]+a[1]*b[1]
def n2(a):    return a[0]**2+a[1]**2

def integrate_kernel(kern, L=4000.0, N=1):
    """int d^2 z kern(z) over the plane, polar with log radial grid."""
    def f(u,phi):
        rho=np.exp(u); z=(rho*np.cos(phi), rho*np.sin(phi))
        return kern(z)*rho**2
    v,_=integrate.dblquad(f,0.0,2*np.pi, lambda p:-14.0, lambda p:np.log(L),
                          epsabs=1e-10, epsrel=1e-10)
    return v

# ---- CSSV identity (3.27) ---------------------------------------------------
print("=== (6) identity (3.27):  (1/pi) int d^2z [ r_zx.r_zx'/(r_zx^2 r_zx'^2) - r_zx.r_zy/(r_zx^2 r_zy^2) ]  =  ln(r_xy^2/r_xx'^2) ===")
cases = [((0,0),(0.3,0.1),(1.0,0.0)),      # x, x', y
         ((0,0),(1.7,-0.4),(0.5,0.9)),
         ((0.2,0.2),(-0.6,0.5),(1.3,-1.1))]
for x,xp,y in cases:
    x,xp,y = map(np.array,(x,xp,y))
    def k(z, x=x,xp=xp,y=y):
        zx,zxp,zy = z-x, z-xp, z-y
        return dot(zx,zxp)/(n2(zx)*n2(zxp)) - dot(zx,zy)/(n2(zx)*n2(zy))
    num = integrate_kernel(k)/np.pi
    ana = np.log(n2(x-y)/n2(x-xp))
    print(f"   num={num: .8f}   ana={ana: .8f}   diff={num-ana: .2e}")

# ---- the full dipole kernel K_dip of (3.28) ---------------------------------
print("\n=== (7) (1/pi) int d^2z K_dip  =  ln( r_xy^2 r_x'y'^2 / (r_xx'^2 r_yy'^2) ) ===")
cases2 = [((0,0),(1.0,0.0),(0.25,0.15),(1.1,-0.2)),        # x, y, x', y'
          ((0.1,-0.3),(0.9,0.7),(-0.4,0.2),(1.4,0.5)),
          ((0,0),(2.0,0.0),(0.05,0.02),(2.03,-0.04))]
for x,y,xp,yp in cases2:
    x,y,xp,yp = map(np.array,(x,y,xp,yp))
    def kd(z, x=x,y=y,xp=xp,yp=yp):
        zx,zy,zxp,zyp = z-x, z-y, z-xp, z-yp
        return (dot(zx,zxp)/(n2(zx)*n2(zxp)) - dot(zx,zy)/(n2(zx)*n2(zy))
              + dot(zy,zyp)/(n2(zy)*n2(zyp)) - dot(zxp,zyp)/(n2(zxp)*n2(zyp)))
    num = integrate_kernel(kd)/np.pi
    ana = np.log(n2(x-y)*n2(xp-yp)/(n2(x-xp)*n2(y-yp)))
    print(f"   num={num: .8f}   ana={ana: .8f}   diff={num-ana: .2e}")
