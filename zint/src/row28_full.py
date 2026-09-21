"""Row V (running-coupling row): full check of the note + the r integration.

x at the origin  =>  r = z-x = z ,  s = y-x = y ,  P = x'-w' ,  xi = p+/k+ , k+ = 1.
"""
import numpy as np, sympy as sp
np.set_printoptions(precision=10, suppress=False)
s_=lambda v: v@v
dd=np.eye(2)

# ----------------------------------------------------------------- the exact bracket
def bracket(y,xv,zv,xi,dperp=2.0):
    """The six groups of eq (1.1), EXACTLY as written there (indices i,m)."""
    xb=1-xi
    A=s_(y-xv); B=s_(y-zv); C=s_(xv-zv)
    yx,yz,xz = y-xv, y-zv, xv-zv
    G=[]
    G.append( dd*dperp/(2*C)*(A-B)*xi*xb )                                        # (1.4)
    G.append( -(xi/xb)*( np.outer(yx,xz)/C - np.outer(yx,yz)/(2*B) ) )            # (1.5)
    G.append( xb*dd*( (yz@xz)/C + (yx@yz)/(2*A) - (A-B)/(2*C) ) )                 # (1.6)
    G.append( -(xb/xi)*( np.outer(yz,xz)/C + np.outer(yz,yx)/(2*A) ) )            # (1.8)
    G.append( xi*dd*( (yx@xz)/C - (yx@yz)/(2*B) - (A-B)/(2*C) ) )                 # (1.9)
    G.append( -( np.outer(xz,yx)/C - np.outer(yz,yx)/(2*B)
                +np.outer(xz,yz)/C + np.outer(yx,yz)/(2*A) ) )                    # (1.10)
    return G

def groups_T(r,s,P,xi,dperp=2.0):
    """contract each group with P^i r^m /(P^2 r^2)"""
    xv=np.zeros(2); zv=r; y=s
    G=bracket(y,xv,zv,xi,dperp)
    w=np.outer(P/s_(P), r/s_(r))
    return np.array([np.einsum('im,im->',w,g) for g in G])

# ----------------------------------------------------------------- the note's rewritings
def note_terms(r,s,P,xi,dperp=2.0,fix18=False,fix19=False):
    xb=1-xi; P2=s_(P); r2=s_(r); s2=s_(s); sr=s@r; Pr=P@r; Ps=P@s; smr=s-r; smr2=s_(smr)
    t4 = dperp*xi*xb/2*( 2*Pr*sr/(P2*r2**2) - Pr/(P2*r2) )
    t5 = (xi/xb)*( Ps/(P2*r2) + Ps*(sr-r2)/(2*P2*r2*smr2) )
    t6 = xb*(Pr/(P2*r2))*( 2 - 2*sr/r2 - sr/(2*s2) )
    half = 0.5 if fix18 else 1.0
    t8 = (xb/xi)*((Ps-Pr)/(P2*r2))*( 1 - half*sr/s2 )
    mid = (s2-sr)/(2*r2*smr2) if fix19 else (s2-sr)/(2*smr2)
    t9 = -xi*(Pr/P2)*( 2*sr/r2**2 + mid - 1/(2*r2) )
    t10= ( 2*Pr*sr/(P2*r2**2) - Pr/(P2*r2) + Ps/(2*P2*s2) - Ps*sr/(2*P2*r2*s2)
          + sr/(2*P2*r2*smr2)*(Ps-Pr) )
    return np.array([t4,t5,t6,t8,t9,t10])

print("="*78); print("1.  GROUP-BY-GROUP:  your (1.4)-(1.10)  vs  the exact contraction"); print("="*78)
rng=np.random.default_rng(3)
lab=["(1.4)","(1.5)","(1.6)/(1.7)","(1.8)","(1.9)","(1.10)"]
bad=np.zeros(6); badfix=np.zeros(6)
for t in range(400):
    r,s,P=[rng.normal(size=2)*1.1 for _ in range(3)]; xi=rng.uniform(.15,.85)
    ex=groups_T(r,s,P,xi)
    bad   =np.maximum(bad,   np.abs(note_terms(r,s,P,xi)-ex))
    badfix=np.maximum(badfix,np.abs(note_terms(r,s,P,xi,fix18=True,fix19=True)-ex))
for l,b,bf in zip(lab,bad,badfix):
    print("   %-12s  as written: %10.3e     with my two fixes: %10.3e"%(l,b,bf))
print()
print("   => (1.8) and (1.9) are the only wrong lines.")
print("      (1.8):  [1 - (r.s)/s^2]          should be   [1 - (r.s)/(2 s^2)]")
print("      (1.9):  (s^2-s.r)/[2 (s-r)^2]    should be   (s^2-s.r)/[2 r^2 (s-r)^2]")

# =================================================================================
print(); print("="*78); print("2.  THE r -> 0 ASYMPTOTICS OF  T(r)/(xibar D)"); print("="*78)
def T_exact(r,s,P,xi,dperp=2.0):   return groups_T(r,s,P,xi,dperp).sum()
def Dfun(r,s,xi):                  return s_(r-s) + (xi/(1-xi))*s_(s)
def F_exact(r,s,P,xi):             return T_exact(r,s,P,xi)/((1-xi)*Dfun(r,s,xi))
CUV=lambda x: x*(1-x)+x/(1-x)+(1-x)/x
def T_uv(r,s,P,xi,dperp=2.0):
    P2=s_(P); r2=s_(r)
    return ( dperp*xi*(1-xi)*(P@r)*(s@r)/(P2*r2**2)
            +(xi/(1-xi)+(1-xi)/xi)*(P@s)/(P2*r2) )
s0=np.array([0.7,-0.4]); P0=np.array([-0.3,0.9])
print("   angular average of  rho^2 * F(rho)  ,  compared with  C_UV(xi) (P.s)/(P^2 s^2)")
print("   %-6s %14s %14s %14s %16s"%("xi","rho=1e-2","1e-3","1e-4","predicted"))
for xi in (0.1,0.25,0.5,0.75,0.9):
    row=[]
    for rho in (1e-2,1e-3,1e-4):
        th=np.linspace(0,2*np.pi,4000,endpoint=False)
        rr=rho*np.stack([np.cos(th),np.sin(th)],1)
        row.append(np.mean([F_exact(v,s0,P0,xi) for v in rr])*rho*rho)
    pred=CUV(xi)*(P0@s0)/(s_(P0)*s_(s0))
    print("   %-6.2f %14.8f %14.8f %14.8f %16.8f"%(xi,*row,pred))
print()
print("   and the SAME check for the subtracted integrand rho^2*[F - T_uv/s^2]  (must -> 0):")
for xi in (0.25,0.5,0.75):
    row=[]
    for rho in (1e-2,1e-3,1e-4):
        th=np.linspace(0,2*np.pi,4000,endpoint=False)
        rr=rho*np.stack([np.cos(th),np.sin(th)],1)
        row.append(np.mean([F_exact(v,s0,P0,xi)-T_uv(v,s0,P0,xi)/s_(s0) for v in rr])*rho*rho)
    print("   xi=%-6.2f %14.3e %14.3e %14.3e"%(xi,*row))

# =================================================================================
print(); print("="*78); print("3.  THE TWO FOURIER MASTERS IN d = 2-2eps  (cutoff cross-check)"); print("="*78)
import mpmath as mp
mp.mp.dps=25
def osc(f,q,lo):     # log-singular near lo, oscillatory tail: split at a few periods
    T0=lo+8*2*mp.pi/q
    return mp.quad(f,[lo,T0])+mp.quadosc(f,[T0,mp.inf],period=2*mp.pi/q)
def M0(q,rho0):      # int_{|r|>rho0} d^2r e^{iqr}/r^2  = 2pi int dr J0(qr)/r
    return 2*mp.pi*osc(lambda t: mp.j0(q*t)/t, q, rho0)
def M2(q,rho0):      # int d^2r e^{iqr} r^i r^m /r^4 : coefficients of delta^{im} and qhat qhat
    a=osc(lambda t: mp.j1(q*t)/(q*t*t), q, rho0)
    b=osc(lambda t:(mp.j0(q*t)-2*mp.j1(q*t)/(q*t))/t, q, rho0)
    return 2*mp.pi*a, 2*mp.pi*b
q=mp.mpf('1.3'); rho0=mp.mpf('1e-6')
L=mp.log(2/(q*rho0))-mp.euler
print("   cutoff:   M0            = %s     2*pi*L = %s"%(mp.nstr(M0(q,rho0),12), mp.nstr(2*mp.pi*L,12)))
A,B=M2(q,rho0)
print("   cutoff:   M2 delta-part = %s     pi*(L+1/2) = %s"%(mp.nstr(A,12), mp.nstr(mp.pi*(L+mp.mpf(1)/2),12)))
print("   cutoff:   M2 qhat-part  = %s     -pi        = %s"%(mp.nstr(B,12), mp.nstr(-mp.pi,12)))
print()
print("   dim reg (MSbar), with  Lq = 1/eps + log(q^2/mubar^2) ,  d = 2-2eps :")
print("       int d^d r e^{iqr}/r^2        = -pi Lq")
print("       int d^d r e^{iqr} r^i r^m/r^4 = -(pi/2) Lq delta^{im} - pi qhat^i qhat^m")
print("   trace check in d dims:  -(pi/2)Lq*d - pi = -pi Lq + pi(eps Lq) - pi = -pi Lq   (eps Lq -> 1)  OK")
print("   cutoff <-> dimreg :  -Lq  <->  2L ;  then the delta-part is pi L , not pi(L+1/2):")
print("   the extra +pi/2 delta^{im} of the cutoff is exactly the pi from  d*(-pi/2 Lq) = -pi Lq + pi .")
print("   Both schemes give the SAME finite leftover below; only the bookkeeping differs.")

# the assembled UV integral
print()
print("   ASSEMBLY:   int d^d r e^{iqr} T_uv(r)")
print("     = d*xi*xibar*(P^i s^m/P^2)*[-(pi/2)Lq delta - pi qhat qhat]  +  A*(P.s/P^2)*(-pi Lq)")
print("     = -pi Lq (P.s/P^2)*[ xi xibar + xi/xibar + xibar/xi ]   +  pi xi xibar (P.s - 2(P.qhat)(s.qhat))/P^2")
print("     = pi C_UV(xi) (P.s/P^2) [ 1/eps_UV + log(mubar^2/q^2) ] + pi xi xibar [P.s-2(P.khat)(s.khat)]/P^2")
print("   (1/eps_UV = -1/eps ;  qhat = khat since q = xibar k)")

# =================================================================================
print(); print("="*78); print("4.  THE xi INTEGRALS  (lambda = Lambda/k+)"); print("="*78)
xi,lam=sp.symbols('xi lamda',positive=True)
Cuv=xi*(1-xi)+xi/(1-xi)+(1-xi)/xi
I1=sp.integrate(Cuv,(xi,lam,1-lam))
I1s=sp.simplify(sp.series(sp.simplify(I1),lam,0,1).removeO())
print("   int_lam^{1-lam} C_UV dxi            = %s"%sp.simplify(I1s))
print("                                       = 2 log(k+/Lambda) - 11/6")
pieces=[(-2*xi/(1-xi)*sp.log(1-xi),"xi/xibar"),(-2*(1-xi)/xi*sp.log(1-xi),"xibar/xi"),
        (-2*xi*(1-xi)*sp.log(1-xi),"xi xibar")]
tot=0
for e,nm in pieces:
    v=sp.integrate(e,(xi,lam,1-lam))
    v=sp.simplify(sp.series(sp.expand(sp.simplify(v)),lam,0,1).removeO()); tot+=v
    print("      %-10s -> %s"%(nm,sp.simplify(v)))
I2s=sp.simplify(tot)
print("   int_lam^{1-lam} C_UV*(-2 log xibar) = %s"%I2s)
print("                                       = log^2(k+/Lambda) + pi^2/3 - 67/18")
print("   numeric check at lam=1e-6:")
f1=mp.quad(lambda x: x*(1-x)+x/(1-x)+(1-x)/x,[1e-6,0.5,1-1e-6])
f2=mp.quad(lambda x:(x*(1-x)+x/(1-x)+(1-x)/x)*(-2*mp.log(1-x)),[1e-6,0.5,1-1e-6])
print("      int C_UV            = %-18s   2 log(1e6)-11/6      = %s"%(mp.nstr(f1,10),mp.nstr(2*mp.log(1e6)-mp.mpf(11)/6,10)))
print("      int C_UV(-2 log xb) = %-18s   log^2(1e6)+pi^2/3-67/18 = %s"%(mp.nstr(f2,10),mp.nstr(mp.log(1e6)**2+mp.pi**2/3-mp.mpf(67)/18,10)))
print("      int_0^1 xi xibar    = 1/6   (the coefficient of the new finite tensor term)")
