# One-dimensional finite-volume check of Sec. 2.2 (referee point 2): for each (L, b) finds the
# ground state of the dimensionless functional (6) and reports phi_0(0) mod 2pi, the number of
# solitons, and whether the gradient has a max or min at the center. Requires numpy, scipy.
import numpy as np
from scipy.optimize import minimize
def E(phi, h, b):
    d = np.diff(phi)/h; pot = 1-np.cos(phi); w=np.ones_like(pot); w[0]=w[-1]=0.5
    return h*np.sum(0.5*d*d - b*d) + h*np.sum(w*pot)
def grad(phi,h,b):
    d=np.diff(phi)/h; g=np.zeros_like(phi); g[:-1]-= (d-b); g[1:]+= (d-b)
    w=np.ones_like(phi); w[0]=w[-1]=0.5; return g + h*w*np.sin(phi)
def ground(L,b,n=201):
    z=np.linspace(-L/2,L/2,n); h=z[1]-z[0]; best=None
    for N in range(0,8):
        for shift in (0,):
            for slope in (0.3,1,2):
                # N-soliton guess: staircase of N kinks, equally spaced
                if N==0: x0=0.01*z
                else:
                    c=np.linspace(-L/2,L/2,N+2)[1:-1]; x0=sum(4*np.arctan(np.exp(slope*(z-ci))) for ci in c)
                r=minimize(E,x0,args=(h,b),jac=grad,method='L-BFGS-B',options={'maxiter':20000,'gtol':1e-10})
                if best is None or r.fun<best[0]: best=(r.fun,r.x)
    return z,best
for L,b in [(5,0.3),(5,1.0),(5,2.0),(10,1.0),(10,1.5),(10,3),(10,5),(2,2),(2,4)]:
    z,(e,phi)=ground(L,b)
    phi0=phi[len(z)//2] % (2*np.pi); g=np.diff(phi)
    k=np.floor((phi-np.pi)/(2*np.pi)); nsol=int(k.max()-k.min())
    print(f"L={L} b={b}: E={e:.4f}  phi(0) mod 2pi={phi0:.3f}  solitons(crossings of pi mod 2pi)={nsol}  grad at center {'max' if g[len(g)//2]>g[len(g)//2-5] and g[len(g)//2]>=g[len(g)//2+5] else 'min'}")
