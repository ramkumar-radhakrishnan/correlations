"""Group II Row IV: the two 1/epsbar's and where each can go."""
import numpy as np, sympy as sp, mpmath as mp
print("="*76); print("1.  THE BK IDENTITY: real + virtual is large-z convergent"); print("="*76)
a1,a2,b1,b2=sp.symbols('a1 a2 b1 b2',real=True)
a=sp.Matrix([a1,a2]); b=sp.Matrix([b1,b2]); A=(a.T*a)[0]; B=(b.T*b)[0]; AB=(a.T*b)[0]
lhs=AB/(A*B)-sp.Rational(1,2)/A-sp.Rational(1,2)/B
rhs=-sp.Rational(1,2)*((a-b).T*(a-b))[0]/(A*B)
print("   (y-z).(y'-z)/[(y-z)^2(y'-z)^2] - 1/2(y-z)^2 - 1/2(y'-z)^2  +  (y-y')^2/2(y-z)^2(y'-z)^2 =",
      sp.simplify(lhs-rhs))
print("   -> the sum is  -(y-y')^2 / [2 (y-z)^2 (y'-z)^2] ~ 1/z^4 : the log at |z| -> infinity CANCELS.")
print()
def trace_cut(r,R,n=2500,m=768,virtual=False):
    lo,hi=np.log(1e-8),np.log(R)
    xg,wg=np.polynomial.legendre.leggauss(n)
    lr=0.5*(hi-lo)*xg+0.5*(hi+lo); wr=0.5*(hi-lo)*wg; rho=np.exp(lr)
    ph=(np.arange(m)+0.5)*(2*np.pi/m)
    ux=np.outer(rho,np.cos(ph)); uy=np.outer(rho,np.sin(ph))     # u = y-z
    vx=ux+r; vy=uy                                                # v = y'-z , r = y'-y... sign irrelevant
    A=ux**2+uy**2; B=vx**2+vy**2
    f=(ux*vx+uy*vy)/(A*B)
    if virtual: f=f-0.5/A-0.5/B
    w=(wr*rho*rho)[:,None]*(2*np.pi/m)
    return float((f*w).sum())
print("   %-10s %16s %16s"%("R","real only","real + virtual"))
for R in (20.,60.,180.,540.):
    print("   %-10.0f %16.5f %16.5f"%(R,trace_cut(1.0,R),trace_cut(1.0,R,virtual=True)))
print("   real alone grows like pi log R^2 ; the BK combination saturates.  Confirmed.")

print(); print("="*76); print("2.  WHAT THE FF SUBTRACTION LEAVES"); print("="*76)
e,Lr,P=sp.symbols('epsbarinv L_r P_gg')
bare = P/2*(e+1-Lr)
ct   = P/2*e                                  # MSbar FF counterterm, -alpha_s/2pi (1/epsbar) Pgg (x) LO
print("   rows 2+4 of the quoted result, together :  Pgg/2 [ 1/epsbar + 1 - log(4 pi^2 mu^2 r^2) ]")
print("   minus the MSbar FF counterterm Pgg/(2 epsbar) :", sp.simplify(bare-ct))
print("   = Pgg/2 [ 1 - log(4 pi^2 mu_F^2 r^2) ]   -- finite, and carries log mu_F .  <-- what the advisor wants")

print(); print("="*76); print("3.  THE WHOLE ROW HAS ONE 1/epsbar, AND IT FACTORISES"); print("="*76)
Nc,g,pi,K,x=sp.symbols('N_c g pi kplus xi'); XX,RR,dp=sp.symbols('XX RR d_perp')
epsb,Lr=sp.symbols('epsbarinv L_r'); plus=sp.Symbol('plusterm')
xb=1-x; Chat=x*xb+plus+xb/x; M=epsb+1-Lr
full = Chat*M*XX - dp*x*xb*RR - (plus+xb/x)*XX          # verified non-delta integrand
cpt  = Chat*(epsb-Lr)*XX + x*xb*(XX-dp*RR)              # the compact claim
print("   verified integrand - [ Chat (1/epsbar - L_r) XX + xi xibar (XX - d_perp RR) ] =",
      sp.simplify(sp.expand(full-cpt)))
print("   the '+1' of the bracket cancels against the -(xi/[1-xi]_+ + xibar/xi) block.")
print()
print("   and the evolution term carries the SAME bracket:  delta(1-xi) log(k+/Lam) (1/epsbar - L_r) XX")
print("   so the whole row is")
print("     { delta(1-xi) log(k+/Lambda) + Pgg(xi)/(2 Nc) } x [ 1/epsbar - log(4 pi^2 mu^2 r^2) ] (X.X')")
print("     + xi(1-xi) [ (X.X') - d_perp (X.rhat)(X'.rhat) ]")
print("   and  1/epsbar - log(4 pi^2 mu^2 r^2)  =  (1/pi) int d^2z (y-z).(y'-z)/[(y-z)^2 (y'-z)^2] .")
