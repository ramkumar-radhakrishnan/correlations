"""int d^2z  K(w',w;z) K(x,y;z),   K(a,b;z) = (a-z).(b-z)/[(a-z)^2 (b-z)^2].

Complex coordinates:  a vector v^i/v^2 is the complex number 1/vbar, and
u.v/(u^2 v^2) = Re[1/(ubar v)].  Hence

   K(a,b;z) = Re 1/[(zbar - abar)(z - b)]

so the product of two kernels is (1/2) Re [ (a) + (b) ] with

   (a) = 1/[(zbar-w'bar)(zbar-xbar)(z-w)(z-y)]
   (b) = 1/[(zbar-w'bar)(zbar-ybar)(z-w)(z-x)]

and each is the master

   M(a1,a2;b1,b2) = int d^2z / [(zbar-a1bar)(zbar-a2bar)(z-b1)(z-b2)]
     = -pi/[(a1bar-a2bar)(b1-b2)] * log[ |a1-b1|^2 |a2-b2|^2 / (|a1-b2|^2 |a2-b1|^2) ]

a pure cross-ratio, so the log R^2 of the individual pieces cancels.
"""
import numpy as np
from scipy.integrate import quad

def Kvec(a, b, z):                       # real-vector form
    A = a - z; B = b - z
    return np.sum(A*B, -1)/(np.sum(A*A, -1)*np.sum(B*B, -1))

def Kcpx(a, b, z):                       # complex form, must agree
    return (1/((np.conj(z)-np.conj(a))*(z-b))).real

def L_cut(a, b, R):                      # int_{|z|<R} d^2z / [(zbar-abar)(z-b)]
    c = b - a
    return np.pi*np.log(R**2/abs(c)**2)  # claim

def M(a1, a2, b1, b2):
    num = abs(a1-b1)**2*abs(a2-b2)**2
    den = abs(a1-b2)**2*abs(a2-b1)**2
    return -np.pi/((np.conj(a1)-np.conj(a2))*(b1-b2))*np.log(num/den)

def closed_form(wp, w, x, y):
    """the z-integral, as complex points"""
    return 0.5*(M(wp, x, w, y) + M(wp, y, w, x)).real

def closed_form_vec(wp, w, x, y):
    """same thing written with real 2-vectors"""
    def K(u, v):                          # u.v/(u^2 v^2) for 2-vectors
        return (u@v)/((u@u)*(v@v))
    def s2(u): return u@u
    t1 = K(wp-x, w-y)*np.log(s2(w-wp)*s2(x-y)/(s2(wp-y)*s2(w-x)))
    t2 = K(wp-y, w-x)*np.log(s2(w-wp)*s2(x-y)/(s2(wp-x)*s2(w-y)))
    return -np.pi/2*(t1 + t2)

# ---------------------------------------------------------------- numerics
def integrate_partition(f, pts, ell=1.0, ntheta=3001, Rmax=1e7):
    """int d^2z f(z) with integrable 1/rho singularities at pts.
    Partition of unity chi_i = g_i / sum(g), g_i = 1/|z-p_i|, g_0 = 1/ell."""
    pts = [np.asarray(p, float) for p in pts]
    th = np.linspace(0, 2*np.pi, ntheta, endpoint=False)
    u = np.stack([np.cos(th), np.sin(th)], -1)
    def gsum(z):
        return 1.0/ell + sum(1.0/np.linalg.norm(z-p, axis=-1) for p in pts)
    def radial(center, chi_i):
        def g(rho):
            z = center + u*rho
            return rho*np.mean(chi_i(z)*f(z))*2*np.pi
        tot = 0.0
        edges = np.concatenate(([1e-12], np.geomspace(1e-10, Rmax, 400)))
        for a, b in zip(edges[:-1], edges[1:]):
            tot += quad(g, a, b, limit=200, epsabs=1e-14, epsrel=1e-11)[0]
        return tot
    tot = 0.0
    for p in pts:
        tot += radial(p, lambda z, p=p: (1.0/np.linalg.norm(z-p, axis=-1))/gsum(z))
    cen = sum(pts)/len(pts)
    tot += radial(cen, lambda z: (1.0/ell)/gsum(z))
    return tot

if __name__ == "__main__":
    rng = np.random.default_rng(7)
    print("0. the complex form of the kernel equals the vector form")
    for _ in range(3):
        a, b, z = rng.normal(size=(3, 2))
        ac, bc, zc = [complex(*v) for v in (a, b, z)]
        print(f"   vector {Kvec(a,b,z):+.12f}   complex {Kcpx(ac,bc,zc):+.12f}")
    print()

    print("1. L(a,b) = int_{|z|<R} d^2z/[(zbar-abar)(z-b)]  =  pi log(R^2/|a-b|^2)")
    a, b = complex(0.3,-0.2), complex(-0.5,0.8)
    for R in (20.0, 100.0, 500.0):
        av = np.array([a.real, a.imag]); bv = np.array([b.real, b.imag])
        def f(z):
            zc = z[...,0] + 1j*z[...,1]
            return (1/((np.conj(zc)-np.conj(a))*(zc-b))).real
        th = np.linspace(0, 2*np.pi, 4001, endpoint=False)
        uu = np.stack([np.cos(th), np.sin(th)], -1)
        def g(rho):
            zz = uu*rho
            v = 1/((zz[...,0]-1j*zz[...,1]-np.conj(a))*(zz[...,0]+1j*zz[...,1]-b))
            return rho*np.mean(v.real)*2*np.pi
        edges = np.concatenate(([1e-9], np.geomspace(1e-8, R, 500)))
        tot = sum(quad(g, p, q, limit=200)[0] for p, q in zip(edges[:-1], edges[1:]))
        print(f"   R={R:7.1f}: quad {tot:+.9f}   pi log(R^2/|a-b|^2) {np.pi*np.log(R**2/abs(a-b)**2):+.9f}")
    print()

    print("2. THE FULL z-INTEGRAL  vs  the closed form")
    for trial in range(4):
        P = rng.normal(size=(4, 2))*1.2
        wp, w, x, y = P
        cx = [complex(*v) for v in P]
        f = lambda z: Kvec(wp, z, z)*0  # placeholder
        def integrand(z):
            return Kvec(wp, w, z)*Kvec(x, y, z)
        num = integrate_partition(integrand, [wp, w, x, y], ell=1.0)
        an1 = closed_form(*cx)
        an2 = closed_form_vec(wp, w, x, y)
        print(f"   trial {trial}: quad {num:+.9f}   closed(complex) {an1:+.9f}   "
              f"closed(vector) {an2:+.9f}   rel {abs(num-an2)/abs(an2):.2e}")
