"""Row 8: the z tensor integral in d = 2-2eps, the contraction, and the p+ integrals."""
import numpy as np, sympy as sp

# ---------------------------------------------------------------- 1. T symmetric?
rng = np.random.default_rng(3)
def Tnum(y,yp,R=3000.0,nr=3000,nth=512):
    rr=np.geomspace(1e-7,R,nr); th=np.linspace(0,2*np.pi,nth,endpoint=False)
    RR,TH=np.meshgrid(rr,th,indexing='ij')
    Z=np.stack([RR*np.cos(TH),RR*np.sin(TH)],-1)
    U=y-Z; V=yp-Z; u2=(U**2).sum(-1); v2=(V**2).sum(-1)
    F=U[...,:,None]*V[...,None,:]/(u2*v2)[...,None,None]
    w=(RR*np.gradient(rr)[:,None]*(2*np.pi/nth))[...,None,None]
    return (F*w).sum((0,1))
y=np.array([0.8,-0.5]); yp=np.array([-0.2,0.7]); r=y-yp; r2=r@r
Tn=Tnum(y,yp)
pred=np.pi/2*(np.log(3000.0**2/r2)+1)*np.eye(2)-np.pi*np.outer(r,r)/r2
print("="*74); print("1.  T^{mm'} = int d^2z (y-z)^m (y'-z)^{m'}/[(y-z)^2 (y'-z)^2]"); print("="*74)
print("   numeric :", np.round(Tn,5).tolist())
print("   closed  :", np.round(pred,5).tolist())
print("   antisymmetric part of numeric T : %.3e   ->  T IS SYMMETRIC" %
      np.abs(Tn-Tn.T).max())

# ------------------------------------------------- 2. dim reg, d = 2-2eps
eps,r2s,mu=sp.symbols('epsilon r2 mu',positive=True)
d=2-2*eps
ceps=sp.pi**(1-eps)*sp.gamma(1-eps)**2/sp.gamma(2-2*eps)
Tpref=ceps*sp.gamma(1+eps)*(r2s)**(-eps)                # T = Tpref*[delta/(2 eps) - rhat rhat]
trace=sp.simplify(Tpref*(d/(2*eps)-1))
trace_direct=sp.pi**(1-eps)*sp.gamma(1+eps)*sp.gamma(1-eps)**2/sp.gamma(1-2*eps)/eps*(r2s)**(-eps)
print(); print("="*74); print("2.  DIM REG  (d = 2-2eps, eps>0 regulates the LARGE-|z| region)"); print("="*74)
print("   T^{mm'} = C(eps) (r^2)^(-eps) [ delta^{mm'}/(2 eps) - rhat^m rhat^{m'} ],")
print("   C(eps) = pi^(1-eps) Gamma(1+eps)Gamma(1-eps)^2/Gamma(2-2eps)")
print("   trace cross-check against the scalar bubble :",
      sp.simplify(trace-trace_direct)==0)
ser=sp.series(sp.simplify(Tpref/(2*eps)).subs(mu,1),eps,0,1).removeO()
print("   delta-part expansion  :", sp.nsimplify(sp.expand(sp.simplify(ser))))
print("   (cutoff gave  (pi/2)(log R^2/r^2 + 1);  1/eps <-> log(R^2 mubar^2)  => IR pole)")

# ------------------------------------------------- 3. contraction with the bracket
print(); print("="*74); print("3.  CONTRACTING  T^{mm'}  WITH THE MOMENTUM BRACKET"); print("="*74)
print("   T is symmetric  =>  delta_{mk'}delta_{m'k} and delta_{mk}delta_{m'k'} give the")
print("   SAME contraction  =>  the whole  (2/k+)(delta_mk' delta_m'k - delta_m'k' delta_km)")
print("   term CANCELS after the z integration.")
# numerical confirmation
Xk=np.array([0.4,1.1]); Xp=np.array([-0.9,0.3])
c2=np.einsum('mn,m,n',Tn,Xp,Xk)-np.einsum('mn,m,n',Tn,Xk,Xp)
print("   numerical check of that cancellation : %.3e" % abs(c2))
P,Kp,dp=sp.symbols('pplus kplus d_perp',positive=True); S=P+Kp
pole=sp.simplify(P/S**2*(d/(2*eps))*(2*eps)/1)   # bookkeeping below instead
XX=sp.Symbol('XdotX')      # (x-y).(x'-y')/[(x-y)^2 (x'-y')^2]
RX=sp.Symbol('RXX')        # (rhat.(x-y)/(x-y)^2)(rhat.(x'-y')/(x'-y')^2)
term1=dp.subs(dp,d)*P/S**2*(XX/(2*eps)-RX)
term4=(P/Kp**2+1/P)*((d/(2*eps))-1)*XX
tot=sp.simplify(sp.expand(term1+term4))
print()
print("   total (in units of C(eps)(r^2)^-eps):")
print("     ", sp.simplify(sp.expand(tot)))
polec=sp.simplify(sp.limit(eps*tot,eps,0))
print("   1/eps POLE coefficient :", sp.factor(sp.simplify(polec/XX)), "* XdotX")
fin=sp.simplify(sp.expand(tot-polec/eps))
print("   finite remainder       :", sp.simplify(fin))
