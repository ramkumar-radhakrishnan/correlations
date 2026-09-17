"""The z integral of two WW kernels: every step, verified."""
import numpy as np, sympy as sp
def disc(R,nr=6000,nth=1440,rmin=1e-8):
    rr=np.geomspace(rmin,R,nr); th=np.linspace(0,2*np.pi,nth,endpoint=False)
    RR,TH=np.meshgrid(rr,th,indexing='ij')
    Z=np.stack([RR*np.cos(TH),RR*np.sin(TH)],-1)
    W=RR*np.gradient(rr)[:,None]*(2*np.pi/nth)
    return Z,W
cx=lambda v: v[...,0]+1j*v[...,1]

print("="*76); print("STEP 1.  WHERE IT DIVERGES"); print("="*76)
y=np.array([0.8,-0.5]); yp=np.array([-0.2,0.7]); r=y-yp; r2=r@r
print("   z -> y   : one kernel ~ 1/rho against measure rho drho   -> convergent")
print("   z -> y'  : same                                          -> convergent")
print("   |z| -> oo: (yhat^m zhat^m')/|z|^2, <zhat^m zhat^m'> = delta/2 != 0,")
print("              measure |z|^2 d|z|/|z|   -> LOGARITHMIC, and INFRARED.")
Z,W=disc(3000.)
for Rt in (1e2,1e3,3e3):
    m=(Z**2).sum(-1)<=Rt*Rt
    U=y-Z; V=yp-Z
    F=(U*V).sum(-1)/((U**2).sum(-1)*(V**2).sum(-1))
    print("      trace, |z| < %-6.0f : %12.6f   vs  pi log(R^2/r^2) = %12.6f"
          %(Rt,(F*W*m).sum(),np.pi*np.log(Rt**2/r2)))

print(); print("="*76); print("STEP 2.  TWO MASTER INTEGRALS (complex coordinates)"); print("="*76)
print("   write z = z1 + i z2.  For a 2D vector a,  a^m/a^2  has  (a^1-ia^2)/|a|^2 = 1/a.")
print()
print("   (a)  int_{|z|<R} d^2z/(z-y) = -pi ybar      [for |y| < R]")
print("        proof: split at |z|=|y|.  For |z|>|y| expand 1/(z-y)=sum y^n/z^{n+1}; every")
print("        angular integral vanishes.  For |z|<|y|, -1/y sum (z/y)^n keeps only n=0,")
print("        giving -(1/y) pi|y|^2 = -pi ybar.")
Z,W=disc(50.); zz=cx(Z); yc=cx(y)
print("        numeric %s   exact %s"%(np.round((W/(zz-yc)).sum(),6), np.round(-np.pi*np.conj(yc),6)))
print()
print("   (b)  L(a,b) = int_{|z|<R} d^2z /[(zbar-abar)(z-b)] = pi log(R^2/|a-b|^2)")
print("        proof: 1/(zbar-abar) = d/dzbar log|z-a|^2 and d/dzbar 1/(z-b) = pi delta^2(z-b).")
print("        Integrate by parts: the bulk term gives -pi log|a-b|^2, the boundary")
print("        |z|=R gives +pi log R^2.")
for Rt in (20.,50.):
    Z,W=disc(Rt); zz=cx(Z); a=cx(y); b=cx(yp)
    num=(W/((np.conj(zz)-np.conj(a))*(zz-b))).sum()
    print("        R=%-5.0f numeric %s   exact %12.6f"%(Rt,np.round(num,5),np.pi*np.log(Rt**2/r2)))

print(); print("="*76); print("STEP 3.  THE TWO CONTRACTIONS"); print("="*76)
print("   trace  T^11+T^22 = Re int d^2z /[(z-y)(zbar-ybar')] = L  = pi log(R^2/r^2)")
print("   and    T^11-T^22 - i(T^12+T^21) = int d^2z/[(z-y)(z-y')]")
print("        = 1/(y-y') [ -pi ybar + pi ybar' ] = -pi rbar/r = -pi e^{-2i phi}")
Z,W=disc(50.); zz=cx(Z); a=cx(y); b=cx(yp)
print("        numeric %s   exact %s"%(np.round((W/((zz-a)*(zz-b))).sum(),6),
                                       np.round(-np.pi*np.conj(a-b)/(a-b),6)))
print()
print("   with T^{mm'} = A delta^{mm'} + B rhat^m rhat^{m'} :")
print("      2A + B = pi log(R^2/r^2)   and   B = -pi   =>   A = (pi/2)[log(R^2/r^2) + 1]")

print(); print("="*76); print("STEP 4.  THE CUTOFF ANSWER, CHECKED"); print("="*76)
R=3000.; Z,W=disc(R)
U=y-Z; V=yp-Z; u2=(U**2).sum(-1); v2=(V**2).sum(-1)
Tn=(U[...,:,None]*V[...,None,:]/(u2*v2)[...,None,None]*W[...,None,None]).sum((0,1))
Tc=np.pi/2*(np.log(R**2/r2)+1)*np.eye(2) - np.pi*np.outer(r,r)/r2
print("   numeric :",np.round(Tn,5).tolist())
print("   closed  :",np.round(Tc,5).tolist())
print("   antisymmetric part of numeric: %.2e  -> T IS SYMMETRIC (this matters in step 6)"
      %np.abs(Tn-Tn.T).max())

print(); print("="*76); print("STEP 5.  DIM REG,  d = 2-2eps,  eps > 0"); print("="*76)
eps,r2s,mu=sp.symbols('epsilon r2 mu',positive=True); d=2-2*eps
C=sp.pi**(1-eps)*sp.gamma(1+eps)*sp.gamma(1-eps)**2/sp.gamma(2-2*eps)
print("   Feynman: 1/(u^2 v^2) = int_0^1 dx/[x u^2 + xbar v^2]^2 ,  shift z' = z - (x y + xbar y'),")
print("            Delta = x xbar r^2 ,  numerator -> z'^m z'^m'/d x z'^2  -  x xbar r^m r^m'")
print("   int d^dl/[l^2+D]^2   = pi^{d/2} Gamma(2-d/2) D^{d/2-2}")
print("   int d^dl l^2/[l^2+D]^2 = pi^{d/2}(d/2)Gamma(1-d/2) D^{d/2-1}")
print("   both x integrals give the same Beta:  int_0^1 (x xbar)^{-eps} = Gamma(1-eps)^2/Gamma(2-2eps)")
print()
print("      T^{mm'} = C(eps) (r^2 mu^2)^{-eps} [ delta^{mm'}/(2 eps) - rhat^m rhat^{m'} ]")
print("      C(eps) = pi^{1-eps} Gamma(1+eps)Gamma(1-eps)^2/Gamma(2-2eps) =", sp.simplify(C))
tr=sp.simplify(C*(d/(2*eps)-1)*r2s**(-eps))
trd=sp.pi**(1-eps)*sp.gamma(1+eps)*sp.gamma(1-eps)**2/sp.gamma(1-2*eps)/eps*r2s**(-eps)
print("   trace cross-check against the scalar bubble (u.v = (u^2+v^2-r^2)/2, scaleless -> 0):",
      sp.simplify(tr-trd)==0)
print("   C(eps) = pi exp[-eps(log pi + gamma_E - 2)] :",
      sp.simplify(sp.series(sp.log(C),eps,0,2).removeO()
                  - (sp.log(sp.pi)-eps*(sp.log(sp.pi)+sp.EulerGamma-2))).simplify())
