"""Row 8: (i) is the 1/eps UV or IR?  (ii) the MS-bar form."""
import numpy as np, sympy as sp

print("="*78); print("1.  UV/IR SPLIT OF THE z INTEGRAL  (direct, numerical)"); print("="*78)
def Tsplit(y, yp, a, R, nr=6000, nth=900):
    """int d^2z over  a < |z-y| , a < |z-y'| , |z| < R."""
    rr = np.geomspace(1e-8, R, nr); th = np.linspace(0,2*np.pi,nth,endpoint=False)
    RR,TH = np.meshgrid(rr,th,indexing='ij')
    Z = np.stack([RR*np.cos(TH), RR*np.sin(TH)], -1)
    U = y-Z; V = yp-Z; u2=(U**2).sum(-1); v2=(V**2).sum(-1)
    mask = (u2 > a*a) & (v2 > a*a)
    F = (U*V).sum(-1)/(u2*v2)*mask                      # the TRACE, delta^{mm'} part
    w = RR*np.gradient(rr)[:,None]*(2*np.pi/nth)
    return (F*w).sum()
y = np.array([0.6,0.0]); yp = np.array([-0.6,0.0]); r2 = ((y-yp)**2).sum()
print("  (a) shrink the excluded discs of radius  a  around y and y'   (R = 400 fixed)")
print("      %-12s %14s %14s" % ("a","trace integral","change"))
prev=None
for a in (1e-1,1e-2,1e-3,1e-4,1e-5):
    v = Tsplit(y,yp,a,400.0)
    print("      %-12.0e %14.6f %14s" % (a, v, "-" if prev is None else "%+.6f"%(v-prev)))
    prev=v
print("      -> a perfectly FLAT limit as a->0 :  NO short-distance (UV) divergence.")
print()
print("  (b) grow the outer radius R   (a = 1e-6 fixed).  pi*log(R^2/r^2) predicted")
print("      %-12s %14s %14s %12s" % ("R","trace integral","pi log(R^2/r^2)","difference"))
for R in (50.,200.,800.,3200.):
    v = Tsplit(y,yp,1e-6,R); p = np.pi*np.log(R**2/r2)
    print("      %-12.0f %14.6f %14.6f %12.2e" % (R, v, p, v-p))
print("      -> a clean log in R :  the divergence is INFRARED.")

print()
print("="*78); print("2.  WHERE A UV POLE *DOES* APPEAR -- and why it is spurious"); print("="*78)
print("""  The trace is usually evaluated by writing  u.v = (u^2 + v^2 - r^2)/2 :

      int d^dz (u.v)/(u^2 v^2) = (1/2)[ int d^dz/v^2 + int d^dz/u^2 ] - (r^2/2) G(1,1)

  Taken apart, BOTH pieces are UV divergent at z -> y (and z -> y'):
    * the scaleless tadpole  int d^dz/u^2  is  pi[ 1/eps_IR - 1/eps_UV ] = 0  only
      because its UV and IR poles are set equal -- it is not actually finite;
    * the bubble  G(1,1) = int d^dz/(u^2 v^2)  goes as  rho^{d-3} drho  at z -> y,
      i.e. it IS log UV divergent at d = 2.
  The two spurious 1/eps_UV's cancel between them, leaving the net 1/eps_IR that part 1
  measured directly.  The physical integrand (u.v)/(u^2 v^2) ~ 1/rho at z -> y is
  integrable, exactly as the numbers in 1(a) show.""")
eps=sp.symbols('epsilon'); r2s=sp.symbols('r2',positive=True)
G11 = -2*sp.pi**(1-eps)/eps*sp.gamma(1+eps)*sp.gamma(1-eps)**2/sp.gamma(1-2*eps)*r2s**(-1-eps)
print("  bubble  G(1,1) = ", sp.simplify(G11))
print("  -(r^2/2) G(1,1) =", sp.simplify(-r2s/2*G11), "   <- this 1/eps is 1/eps_UV")

print()
print("="*78); print("3.  MS-BAR"); print("="*78)
mu=sp.symbols('mu',positive=True)
C   = sp.pi**(1-eps)*sp.gamma(1+eps)*sp.gamma(1-eps)**2/sp.gamma(2-2*eps)
Sms = (sp.exp(sp.EulerGamma)/(4*sp.pi))**eps          # the standard MS-bar factor on g^2 mu^{2eps}
expr = sp.simplify(C*Sms*(r2s*mu**2)**(-eps)/eps)
ser  = sp.expand(sp.series(expr, eps, 0, 1).removeO())
print("   C(eps) * (e^gE/4pi)^eps * (r^2 mu^2)^{-eps} / eps  =")
print("      ", sp.simplify(ser))
print("      = pi [ 1/eps + 2 - log(4 pi^2 mu^2 r^2) ]  :",
      sp.simplify(ser - sp.pi*(1/eps + 2 - sp.log(4*sp.pi**2*mu**2*r2s))) == 0)
L = sp.Symbol('Lcal')     # 1/eps - log(4 pi^2 mu^2 r^2)
X,R_ = sp.symbols('X R')
P,Kp = sp.symbols('pplus kplus',positive=True); S=Kp+P
# CDR / MS-bar:  d_perp = 2-2eps   ->  d_perp/(2 eps) = 1/eps - 1
cdr  = P/S**2*((L+2-1)*X - 2*R_) + (P/Kp**2+1/P)*((L+2)-2)*X
dred = P/S**2*((L+2  )*X - 2*R_) + (P/Kp**2+1/P)*((L+2)-2)*X   # d_perp = 2 kept
print()
print("   CDR / MS-bar (d_perp = 2-2eps):")
print("      ", sp.simplify(sp.expand(cdr)), "   [note the +2's cancelled in the delta_mm' term]")
print("   DRED (d_perp = 2 kept):")
print("      ", sp.simplify(sp.expand(dred)))
print("   scheme difference  (DRED - CDR) =", sp.simplify(dred-cdr), "  <- finite, scheme dependent")
print()
print("   cutoff dictionary:  pi/2 (log R^2/r^2 + 1)  ==  pi/2 [1/eps + 2 - log(4 pi^2 mu^2 r^2)]")
Rc=sp.symbols('R',positive=True)
sol=sp.solve(sp.Eq(sp.log(Rc**2/r2s)+1, 2-sp.log(4*sp.pi**2*mu**2*r2s)), Rc)
print("      =>  R =", [sp.simplify(s) for s in sol if s.is_positive is not False])
