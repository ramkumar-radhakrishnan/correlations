"""Group II Row IV: check the tensor algebra, the w,w' integration, and the contraction."""
import numpy as np, sympy as sp, itertools

print("="*78); print("1.  THE 9-TERM TENSOR PRODUCT  (sum over i,j), in D dimensions"); print("="*78)
def prod(D,a,b):
    """sum_ij [d^km d^ij - a d^jm d^ik - b d^im d^kj][d^k'm' d^ij - a d^jm' d^ik' - b d^im' d^k'j]"""
    d=np.eye(D); out=np.zeros((D,)*4)                      # indices (k,m,k',m')
    for k,m,kp,mp in itertools.product(range(D),repeat=4):
        s=0.0
        for i in range(D):
            for j in range(D):
                A=d[k,m]*d[i,j]-a*d[j,m]*d[i,k]-b*d[i,m]*d[k,j]
                B=d[kp,mp]*d[i,j]-a*d[j,mp]*d[i,kp]-b*d[i,mp]*d[kp,j]
                s+=A*B
        out[k,m,kp,mp]=s
    return out
def model(D,c1,c2,c3):
    d=np.eye(D); out=np.zeros((D,)*4)
    for k,m,kp,mp in itertools.product(range(D),repeat=4):
        out[k,m,kp,mp]=c1*d[k,m]*d[kp,mp]+c2*d[m,mp]*d[k,kp]+c3*d[k,mp]*d[m,kp]
    return out
for D in (3,4,5):
    a,b=1.7,2.3
    c1=D-2*(a+b); c2=a**2+b**2; c3=2*a*b
    err=np.abs(prod(D,a,b)-model(D,c1,c2,c3)).max()
    print("   D=%d  a=%.1f b=%.1f :  coeffs (d_perp-2a-2b, a^2+b^2, 2ab)  max err = %.2e"%(D,a,b,err))
print()
print("   with a=S/p+, b=S/k+, S=p++k+ :  a+b = ab = S^2/(p+ k+)")
S,P,K,dp=sp.symbols('S pplus kplus d_perp',positive=True)
a_,b_=S/P,S/K
print("     a+b =",sp.simplify((a_+b_).subs(S,P+K)),"   ab =",sp.simplify((a_*b_).subs(S,P+K)))
pre=P*K/S**4
c1=sp.simplify(pre*(dp-2*(a_+b_))); c2=sp.simplify(pre*(a_**2+b_**2)); c3=sp.simplify(pre*2*a_*b_)
print("   times the p+k+/S^4 of eq (0.1):")
print("     delta^km delta^k'm' :",sp.simplify(sp.expand(c1)))
print("     delta^mm' delta^kk' :",sp.simplify(sp.expand(c2)))
print("     delta^km' delta^mk' :",sp.simplify(sp.expand(c3)))

print(); print("="*78); print("2.  THE w , w' INTEGRATION"); print("="*78)
beta=P/S
print("   delta^(2)[ y - w + (p+/S)(w-z) ] = delta^(2)[ y - (1-b)w - b z ] ,  b = p+/S")
print("     => w = (y - b z)/(1-b) ,  Jacobian 1/(1-b)^2 = S^2/k+^2   per delta")
print("     => w - z = (y-z)/(1-b)   so   (w-z)^m/(w-z)^2 = (1-b) (y-z)^m/(y-z)^2   per kernel")
print("     => w' - w = (y'-y)/(1-b) = (S/k+)(y'-y)      <-- THE PHASE PICKS UP S/k+")
tot=sp.simplify((1/(1-beta))**4 * (1-beta)**2)
print("   net from the two deltas and the two kernels :", sp.simplify(tot.subs(S,P+K)), " = S^2/k+^2")
print("   so the scalar in front of the 3-term tensor is  (p+k+/S^4)(S^2/k+^2) = %s"%sp.simplify(P*K/S**4*S**2/K**2))
print()
print("   YOUR (0.2) has  (1/k+)*(S/k+)*[bracket].  Comparing term by term:")
yours={'km,k\'m\'': dp*P/S**2 - 2/K, 'mm\',kk\'': P/K**2+1/P, 'km\',mk\'': 2/K}
J=S**2/K**2   # from the two deltas and the two kernels
mine ={'km,k\'m\'': c1*J, 'mm\',kk\'': c2*J, 'km\',mk\'': c3*J}
for k in yours:
    ratio=sp.simplify(sp.expand(sp.simplify((S/K**2*yours[k])/mine[k]).subs(S,P+K)))
    print("     %-10s  yours/mine = %s"%(k,ratio))
print("   => the relative weights are RIGHT; the overall (k+ + p+)/k+ should not be there.")
print("      correct prefactor for your bracket is simply 1/k+ .")

print(); print("="*78); print("3.  WHERE d_perp GOES"); print("="*78)
e=sp.symbols('epsilon'); L=sp.Symbol('L')     # L = 1/eps + 2 - log(4 pi^2 mu^2 r^2)
print("   (a) INSIDE G: the 1/d_perp of eq (0.12) meets the d/2 of master (0.13):")
print("       (1/d_perp)*(d/2)*Gamma[1-d/2] = (1/2)Gamma[eps] = Gamma[1+eps]/(2 eps)")
print("       equivalently  (1/d_perp)(1/eps - 1) =", sp.simplify((1/(2-2*e))*(1/e-1)),"  -- EXACT, d_perp gone")
print("   (b) AT THE CONTRACTION it comes back, twice, and does NOT cancel:")
print("       delta^mm' G^mm' = (pi d_perp/2) L - pi = pi(1-eps)L - pi = pi(L - 2)   [eps*L -> 1]")
print("       d_perp * G^{mm'}X_m X'_m' = 2(1-eps)[...] -> pi(L-1)(X.X') - 2pi (X.r)(X'.r)/r^2")

print(); print("="*78); print("4.  THE CONTRACTION"); print("="*78)
XX,Xr,rr=sp.symbols('XdotXp XrXpr one',positive=True)   # XX=(X.X'), Xr=(X.r)(X'.r)/r^2
b1=dp*P/(K*S**2)-2/K**2; b2=S**2/(P*K**3)-2/K**2; b3=2/K**2
print("   bracket coefficients after w,w' :")
print("     b1 (km,k'm')  =", sp.simplify(b1))
print("     b2 (mm',kk')  =", sp.simplify(sp.expand(b2.subs(S,P+K))), " =  p+/k+^3 + 1/(p+ k+)")
print("     b3 (km',mk')  =", b3)
print("   G^{mm'} is symmetric, so  delta^km delta^k'm' and delta^km' delta^mk'  give the SAME")
print("   contraction G_XX' ; hence b1's -2/k+^2 and b3's +2/k+^2 CANCEL identically.")
print("   (that is your '2/k+ (d^mk' d^m'k - d^m'k' d^km)' term: it contributes nothing.)")
print()
print("   C = (pi p+/(k+ S^2)) [ (L-1)(X.X') - 2 (X.r)(X'.r)/r^2 ]")
print("       + pi [ p+/k+^3 + 1/(p+ k+) ] (L-2) (X.X')")

print(); print("="*78); print("5.  THE p+ INTEGRAL   (sigma = S/k+ , kappa = k.r , lam = Lambda/k+)"); print("="*78)
sig=sp.Symbol('sigma',positive=True)
f1=sp.simplify((P/(K*S**2)).subs(P,K*(sig-1)).subs(S,K*sig))
f2=sp.simplify((P/K**3+1/(P*K)).subs(P,K*(sig-1)))
print("   p+/(k+ S^2)      = %s"%f1)
print("   p+/k+^3+1/(p+k+) = %s  =  [ (sigma-1) + 1/(sigma-1) ]/k+^2"%sp.simplify(f2))
print("   dp+/(2 pi) = k+ dsigma/(2 pi) ,  so the overall factor is (k+/2pi)(pi/k+^2) = 1/(2 k+)")
print()
import mpmath as mp
mp.mp.dps=20
def J1(kap,lo,hi=None):
    f=lambda s:(s-1)/s**2*mp.e**(-1j*kap*s)
    T0=lo+8*2*mp.pi/kap
    return mp.quad(f,[lo,T0])+mp.quadosc(f,[T0,mp.inf],period=2*mp.pi/kap)
def J1c(kap):
    z=1j*kap; return (1+z)*mp.e1(z)-mp.e**(-z)
print("   J1 = int_1^inf dsigma (sigma-1)/sigma^2 e^{-i sigma kappa} = (1+i kappa) E1(i kappa) - e^{-i kappa}")
for kap in (0.3,1.0,2.5):
    num=J1(kap,1); print("      kappa=%.1f   numeric %s    closed %s"%(kap,mp.nstr(num,8),mp.nstr(J1c(kap),8)))
print()
print("   J2 = int_{1+lam} dsigma [ (sigma-1) + 1/(sigma-1) ] e^{-i sigma kappa}")
print("      = e^{-i kappa} { [e^{-i kappa T}(1+i kappa T) - 1]/kappa^2  +  E1(i kappa lam) - E1(i kappa T) }")
print("      E1(i kappa lam) -> log(1/lam) - gammaE - log(i kappa)   <-- the rapidity log log(k+/Lambda)")
def J2(kap,lam,T): return mp.quad(lambda t:(t+1/t)*mp.e**(-1j*kap*(1+t)),[lam,T])
def J2c(kap,lam,T):
    z=1j*kap
    return mp.e**(-z)*((mp.e**(-z*T)*(1+z*T)-mp.e**(-z*lam)*(1+z*lam))/kap**2 + mp.e1(z*lam)-mp.e1(z*T))
for kap,lam,T in ((0.7,1e-3,60.0),(1.8,1e-4,40.0)):
    print("      kappa=%.1f lam=%.0e T=%.0f  numeric %s   closed %s"%(kap,lam,T,mp.nstr(J2(kap,lam,T),8),mp.nstr(J2c(kap,lam,T),8)))
print()
print("   NOTE: drop the sigma in the phase (as eq (0.2)/(0.3) do) and the (sigma-1) piece gives")
print("   int^{V} dp+ p+/k+^3 = V^2/(2 k+^3) -- a QUADRATIC divergence in the p+ cutoff.")
print("   The phase's p+ dependence is what makes that term oscillate and drop.")
