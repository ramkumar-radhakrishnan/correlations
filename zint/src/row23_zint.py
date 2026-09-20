"""Row 23: the master z-integral T^{mm'}, derived and verified two independent ways."""
import numpy as np, sympy as sp
from scipy import integrate
print("="*80); print("PART I  --  d = 2 WITH A CUTOFF, BY COMPLEX COORDINATES"); print("="*80)
print(" I.1  the two complex master integrals, checked numerically")
R=60.0
def num2d(f,R,n=4000,m=4000):
    # polar grid, u = rho e^{i phi}
    rho=(np.arange(n)+0.5)*(R/n); ph=(np.arange(m)+0.5)*(2*np.pi/m)
    RR,PH=np.meshgrid(rho,ph,indexing='ij'); U=RR*np.exp(1j*PH)
    return np.sum(f(U)*RR)*(R/n)*(2*np.pi/m)
for aval in (0.7+0.3j, -1.2+0.8j):
    got=num2d(lambda u: 1.0/(u-aval),R); print("     int_{|u|<R} d^2u /(u-a)          a=%-12s  numeric %10.5f%+10.5fi   exact %10.5f%+10.5fi"
          %(aval,got.real,got.imag,(-np.pi*np.conj(aval)).real,(-np.pi*np.conj(aval)).imag))
for aval,bval in ((-1.3+0.4j,0.0),(0.0,-1.3+0.4j)):
    got=num2d(lambda u: 1.0/((np.conj(u)-np.conj(aval))*(u-bval)),R)
    ex=np.pi*np.log(R**2/abs(aval-bval)**2)
    print("     int d^2u /[(ubar-abar)(u-b)]      numeric %12.6f   exact %12.6f"%(got.real,ex))
print()
print(" I.2  the assembled tensor, checked against  (pi/2)(log R^2/r^2 + 1) delta - pi rhat rhat")
for rv in (np.array([0.9,0.0]), np.array([-0.6,1.1])):
    r2=rv@rv; rh=rv/np.sqrt(r2)
    T=np.zeros((2,2))
    for m in range(2):
        for mp in range(2):
            f=lambda u: (np.real(u) if m==0 else np.imag(u))*( (np.real(u)+rv[0]) if mp==0 else (np.imag(u)+rv[1]) )\
                        /(np.abs(u)**2*np.abs(u+rv[0]+1j*rv[1])**2)
            T[m,mp]=num2d(f,R,6000,6000).real
    ex=np.pi/2*(np.log(R*R/r2)+1)*np.eye(2)-np.pi*np.outer(rh,rh)
    print("     r=%s :  numeric [[%8.4f %8.4f],[%8.4f %8.4f]]"%(rv,T[0,0],T[0,1],T[1,0],T[1,1]))
    print("                 exact [[%8.4f %8.4f],[%8.4f %8.4f]]   max diff %.3f"
          %(ex[0,0],ex[0,1],ex[1,0],ex[1,1],np.abs(T-ex).max()))

print(); print("="*80); print("PART II  --  d = 2-2eps BY FEYNMAN PARAMETERS"); print("="*80)
d,eps,De,n=sp.symbols('d epsilon Delta n',positive=True)
print(" II.1  the Feynman trick   1/(AB) = int_0^1 dx /[xA+(1-x)B]^2")
A,B,x=sp.symbols('A B x',positive=True)
print("       int_0^1 dx/[xA+(1-x)B]^2 - 1/(AB)  =  %s"%sp.simplify(sp.integrate(1/(x*A+(1-x)*B)**2,(x,0,1))-1/(A*B)))
print()
print(" II.2  the Euclidean master  I_n = int d^d l /(l^2+Delta)^n = pi^{d/2} Gamma(n-d/2)/Gamma(n) Delta^{d/2-n}")
ll=sp.Symbol('l',positive=True)
for dd,nn in ((1,2),(3,2),(3,3),(5,4)):
    D2=sp.Rational(dd,2)
    Om=2*sp.pi**D2/sp.gamma(D2)
    direct=sp.simplify(Om*sp.integrate(ll**(dd-1)/(ll**2+De)**nn,(ll,0,sp.oo)))
    form=sp.simplify(sp.pi**D2*sp.gamma(nn-D2)/sp.gamma(nn)*De**(D2-nn))
    print("       d=%d n=%d :  radial %-24s  formula %-24s  diff %s"%(dd,nn,direct,form,sp.simplify(direct-form)))
print()
print(" II.3  int d^d l  l^2/(l^2+Delta)^2 = I_1 - Delta I_2 = pi^{d/2} (d/2) Gamma(1-d/2) Delta^{d/2-1}")
dsym=sp.Symbol('d')
lhs=sp.pi**(dsym/2)*sp.gamma(1-dsym/2)*De**(dsym/2-1) - De*sp.pi**(dsym/2)*sp.gamma(2-dsym/2)*De**(dsym/2-2)
print("       I_1 - Delta I_2 - pi^{d/2}(d/2)Gamma(1-d/2)Delta^{d/2-1}  =  %s"
      %sp.simplify(lhs - sp.pi**(dsym/2)*(dsym/2)*sp.gamma(1-dsym/2)*De**(dsym/2-1)))
print()
print(" II.4  the x integral")
e2=sp.Symbol('epsilon',positive=True)
print("       int_0^1 dx [x(1-x)]^{-eps} = %s   (= Beta(1-eps,1-eps))"
      %sp.simplify(sp.integrate(x**(-e2)*(1-x)**(-e2),(x,0,1),conds='none')))
print("       = Gamma(1-eps)^2/Gamma(2-2eps) = 1/(1-2eps) exactly?  %s"
      %sp.simplify(sp.gamma(1-eps)**2/sp.gamma(2-2*eps)-sp.series(sp.gamma(1-eps)**2/sp.gamma(2-2*eps),eps,0,4).removeO()))
print("       series: %s"%sp.series(sp.gamma(1-eps)**2/sp.gamma(2-2*eps),eps,0,3).simplify())
print()
print(" II.5  the normalisation N(eps) = Gamma(1+eps) Gamma(1-eps)^2/Gamma(2-2eps)")
N=sp.gamma(1+eps)*sp.gamma(1-eps)**2/sp.gamma(2-2*eps)
print("       N(eps) = %s"%sp.series(N,eps,0,2))
mu2,r2=sp.symbols('mu2 r2',positive=True)
Ce=sp.pi*N*(sp.pi*mu2*r2)**(-eps)
print("       C_eps/eps  = %s"%sp.simplify(sp.series(Ce/eps,eps,0,1).removeO()))
print("       i.e.  pi[ 1/eps + 2 - gamma_E - log(pi mu^2 r^2) ]")

print(); print("="*80); print("PART III  --  MATCHING, MS-BAR, AND THE '+2 GETS EATEN' TRAP"); print("="*80)
L=sp.Symbol('Lcal')
print(" III.1 matching:  log(R^2/r^2) + 1  =  1/eps + 2 - gamma_E - log(pi mu^2 r^2)")
print("        => log R^2 = 1/eps + 1 - gamma_E - log(pi mu^2)   and the rhat rhat parts agree EXACTLY")
print(" III.2 MS-bar:  1/eps_MSbar = 1/eps - gamma_E + log 4pi  =>")
print("        1/eps + 2 - gamma_E - log(pi mu^2 r^2) = 1/eps_MSbar + 2 - log(4 pi^2 mu^2 r^2) = Lcal + 2")
print("        with  Lcal = 1/eps_MSbar - log(4 pi^2 mu^2 r^2)")
print()
print(" III.3 THE TRAP.  T^{mm'} = (pi/2)(Lcal+2) delta^{mm'} - pi rhat^m rhat^{m'}  is ONLY safe")
print("       when contracted with eps-INDEPENDENT tensors.  Contract with delta_{mm'} in d=2-2eps:")
tr=sp.simplify(sp.pi/2*(L+2)*(2-2*eps)-sp.pi)
print("         tr T = (pi/2)(Lcal+2) d_perp - pi   with d_perp = 2-2eps :  %s"%tr)
print("         but Lcal contains 1/eps, so -2 eps x (pi/2) x (1/eps) = -pi is a FINITE leftover:")
print("         tr T = pi(Lcal+2) - pi - pi = pi Lcal     <-- the +2 is eaten")
print("       whereas  X_m X'_m' T^{mm'} = (pi/2)(Lcal+2) Xcal - pi Rcal   keeps its +2,")
print("       because X, X' carry no eps.  Getting this wrong shifts the finite part.")
