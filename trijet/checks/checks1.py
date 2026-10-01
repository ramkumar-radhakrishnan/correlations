import numpy as np
from scipy import integrate, special
gam = 0.5772156649015329
c0  = 2*np.exp(-gam)

print("=== (1) int_a^inf J0(x)/x dx  vs  ln(2/a)-gamma ===")
for a in [1e-1,1e-2,1e-3,1e-4]:
    num,_ = integrate.quad(lambda x: special.j0(x)/x, a, 200, limit=4000)
    tail,_= integrate.quad(lambda x: special.j0(x)/x, 200, np.inf, limit=4000)
    num += tail
    ana = -np.log(a/2)-gam
    print(f"  a={a:8.1e}  num={num: .10f}  ana={ana: .10f}  diff={num-ana: .3e}")

print("\n=== (4) 2D:  I(Lam,D) = int d^2C/(2pi)^2 e^{-iC.D}/C^2 Theta(C^2>Lam^2)  vs (1/2pi)ln(c0/(Lam*D)) ===")
def I2d(Lam,D,Cmax=4000.0):
    # angular average of e^{-iC.D} is J0(C D)
    f = lambda C: special.j0(C*D)/C
    v,_ = integrate.quad(f, Lam, Cmax, limit=8000)
    v2,_= integrate.quad(f, Cmax, np.inf, limit=8000)
    return (v+v2)/(2*np.pi)
for (Lam,D) in [(1e-2,1.0),(1e-3,0.7),(1e-2,0.05),(1e-4,2.3)]:
    ana = np.log(c0/(Lam*D))/(2*np.pi)
    print(f"  Lam={Lam:8.1e} D={D:5.2f}  num={I2d(Lam,D): .8f}  ana={ana: .8f}  diff={I2d(Lam,D)-ana: .2e}")

print("\n=== (2) xi-integral:  int_eps^1 dxi/xi (A + 2 ln xi) = A ln(1/eps) - ln^2(eps) ===")
for eps,A in [(1e-3,2.5),(1e-5,-1.3)]:
    num,_ = integrate.quad(lambda x:(A+2*np.log(x))/x, eps,1, limit=2000)
    ana = A*np.log(1/eps) - np.log(eps)**2
    print(f"  eps={eps:.0e} A={A:5.2f}  num={num: .8f}  ana={ana: .8f}  diff={num-ana: .2e}")

print("\n=== (3) DGLAP:  int_eps^1 dxi (1+(1-xi)^2)/xi = 2 ln(1/eps) - 3/2 ===")
for eps in [1e-3,1e-6]:
    num,_ = integrate.quad(lambda x:(1+(1-x)**2)/x, eps,1, limit=2000)
    ana = 2*np.log(1/eps)-1.5
    print(f"  eps={eps:.0e}  num={num: .8f}  ana={ana: .8f}  diff={num-ana: .2e}")
