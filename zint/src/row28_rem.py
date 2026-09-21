"""Row V: the UV-subtracted remainder  I(xi) = int d^2r e^{iqr} [ T/(xibar D) - T_uv/s^2 ]
   and its behaviour at the two endpoints xi -> 0 and xibar -> 0."""
import numpy as np, mpmath as mp
mp.mp.dps=20
s_=lambda v: v@v; dd=np.eye(2)

def T_of(r,s,P,xi,dperp=2.0):
    """corrected T(r) : the six groups of (1.1) contracted with P^i r^m/(P^2 r^2)"""
    xb=1-xi; P2=s_(P); r2=(r*r).sum(-1); s2=s_(s)
    sr=r@s; Pr=r@P; Ps=P@s; smr=s-r; smr2=(smr*smr).sum(-1)
    t4 = dperp*xi*xb/2*( 2*Pr*sr/(P2*r2**2) - Pr/(P2*r2) )
    t5 = (xi/xb)*( Ps/(P2*r2) + Ps*(sr-r2)/(2*P2*r2*smr2) )
    t6 = xb*(Pr/(P2*r2))*( 2 - 2*sr/r2 - sr/(2*s2) )
    t8 = (xb/xi)*((Ps-Pr)/(P2*r2))*( 1 - sr/(2*s2) )
    t9 = -xi*(Pr/P2)*( 2*sr/r2**2 + (s2-sr)/(2*r2*smr2) - 1/(2*r2) )
    t10= ( 2*Pr*sr/(P2*r2**2) - Pr/(P2*r2) + Ps/(2*P2*s2) - Ps*sr/(2*P2*r2*s2)
          + sr/(2*P2*r2*smr2)*(Ps-Pr) )
    return t4+t5+t6+t8+t9+t10
def Tuv_of(r,s,P,xi,dperp=2.0):
    xb=1-xi; P2=s_(P); r2=(r*r).sum(-1)
    return dperp*xi*xb*(r@P)*(r@s)/(P2*r2**2) + (xi/xb+xb/xi)*(P@s)/(P2*r2)
def integrand(r,s,P,k,xi):
    xb=1-xi; D=((r-s)**2).sum(-1)+(xi/xb)*s_(s)
    return T_of(r,s,P,xi)/(xb*D) - Tuv_of(r,s,P,xi)/s_(s)

def inner(s,P,k,xi,R0,nr=900,nth=512):
    """int_{|r|<R0} d^2r e^{iq.r} [ ... ] by polar quadrature (log-graded in r)."""
    q=(1-xi)*k
    u,wu=np.polynomial.legendre.leggauss(nr)
    lo,hi=np.log(1e-9*R0),np.log(R0)
    lr=0.5*(hi-lo)*u+0.5*(hi+lo); wr=0.5*(hi-lo)*wu
    rr=np.exp(lr)
    th=np.arange(nth)*2*np.pi/nth; wt=2*np.pi/nth
    c,sn=np.cos(th),np.sin(th)
    tot=0j
    for R,w in zip(rr,wr):
        pts=np.stack([R*c,R*sn],1)
        f=integrand(pts,s,P,k,xi)
        ph=np.exp(1j*(pts@q))
        tot+= w*R*R*wt*np.sum(f*ph)      # w is dlog r, so r*dr = R^2 dlog r
    return tot

def tail_exact(s,P,k,xi,R0):
    """int_{|r|>R0} d^2r e^{iqr}[F - Tuv/s^2].  For r >> s,M the bracket collapses:
         1/(xb r^2+xi s^2) * [ (xi/xb)(P.s/P^2)(1/2r^2) + (P.s)/(2 P^2 s^2) ]
         - (xi/xb)(P.s)/(P^2 r^2 s^2)   =   -(P.s)/(2 P^2 s^2 xb r^2) ,
       exactly (the 1/(1+u) pieces cancel).   -> -(pi/xb)(P.s/P^2 s^2) int_R0^inf J0(qr) dr/r"""
    xb=1-xi; qn=float(np.sqrt(s_((1-xi)*k))); Ps=P@s
    T0=R0+8*2*mp.pi/qn
    J=mp.quad(lambda t: mp.j0(qn*t)/t,[R0,T0])+mp.quadosc(lambda t: mp.j0(qn*t)/t,[T0,mp.inf],period=2*mp.pi/qn)
    return -(np.pi/xb)*Ps/(s_(P)*s_(s))*float(J)

def tail_uv(s,P,k,xi,R0):
    """ -1/s^2 * int_{|r|>R0} d^2r e^{iq.r} T_uv(r) , exactly, via Bessel integrals."""
    xb=1-xi; q=(1-xi)*k; qn=float(np.sqrt(s_(q))); P2=s_(P); Ps=P@s
    qh=q/qn; Pq=P@qh; sq=s@qh
    # int_{|r|>R0} d^2r e^{iqr}/r^2 = 2pi int_R0^inf J0(q r) dr/r
    A0=2*mp.pi*mp.quadosc(lambda t: mp.j0(qn*t)/t,[R0,mp.inf],period=2*mp.pi/qn)
    # int_{|r|>R0} d^2r e^{iqr} r^i r^m/r^4 = 2pi int [ J1/(qr) delta + (J0-2J1/(qr)) qh qh ] dr/r
    A1=2*mp.pi*mp.quadosc(lambda t: mp.j1(qn*t)/(qn*t*t),[R0,mp.inf],period=2*mp.pi/qn)
    A2=2*mp.pi*mp.quadosc(lambda t:(mp.j0(qn*t)-2*mp.j1(qn*t)/(qn*t))/t,[R0,mp.inf],period=2*mp.pi/qn)
    val = 2*xi*xb*( float(A1)*Ps + float(A2)*Pq*sq )/P2 + (xi/xb+xb/xi)*Ps/P2*float(A0)
    return -val/s_(s)

s0=np.array([0.7,-0.4]); P0=np.array([-0.3,0.9]); k0=np.array([1.0,0.0])
Ps_=(P0@s0)/(s_(P0)*s_(s0))
print("="*78); print("THE REMAINDER  I(xi)  NEAR  xibar -> 0    (s,P,k fixed; k=1)"); print("="*78)
print("   prediction:  xibar*I  ->  -pi*(P.s/P^2 s^2)*log(1/xibar) + const")
print("   (P.s)/(P^2 s^2) = %.8f ,  -pi*that = %.6f per log"%(Ps_,-np.pi*Ps_))
print()
print("   %-10s %-10s %14s %14s %14s"%("xibar","R0","Re I","xibar*Re I","d/dlog(1/xb)"))
prev=None
for xb in (3e-2,1e-2,3e-3,1e-3,3e-4,1e-4):
    xi=1-xb
    R0=float(np.sqrt(s_(s0)))*xb**-0.75
    val=inner(s0,P0,k0,xi,R0)+tail_exact(s0,P0,k0,xi,R0)
    cur=xb*val.real
    d="" if prev is None else "%14.6f"%((cur-prev[1])/np.log(prev[0]/xb))
    print("   %-10.1e %-10.2f %14.6f %14.6f %s"%(xb,R0,val.real,cur,d))
    prev=(xb,cur)

print()
print("="*78); print("THE REMAINDER  I(xi)  NEAR  xi -> 0"); print("="*78)
print("   here q = xibar k -> k , so the EXPLICIT log(mubar^2/q^2) stays finite;")
print("   the only question is whether I(xi) grows faster than 1/xi.")
print("   (tail beyond R0 is killed by the phase: q R0 = 60 >> 1)")
print("   %-10s %14s %14s %14s"%("xi","Re I","xi*Re I","d/dlog(1/xi)"))
prev=None
for xi in (1e-2,3e-3,1e-3,3e-4,1e-4):
    R0=60.0
    val=inner(s0,P0,k0,xi,R0,1600,2048)
    cur=xi*val.real
    d="" if prev is None else "%14.6f"%((cur-prev[1])/np.log(prev[0]/xi))
    print("   %-10.1e %14.6f %14.6f %s"%(xi,val.real,cur,d))
    prev=(xi,cur)
print()
print("   => xi*I(xi) tends to a CONSTANT: no 1/xi * log(1/xi) at this endpoint.")
