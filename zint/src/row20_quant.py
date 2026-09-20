"""Row 20: quantify the surviving large-distance logs of the p+-integrated row."""
import numpy as np, mpmath as mp
from row19_lib import their_coef, exact_M
s=lambda v:v@v
def taus(X,Z,W):
    A,B,C=s(X-Z),s(X-W),s(Z-W)
    return (np.eye(2)*(A-B)/(2*C), np.outer(X-Z,(Z-W)/C-(X-W)/(2*B)),
            np.outer((Z-W)/C+(X-Z)/(2*A),X-W))
con=lambda a,b: np.einsum('ji,ji->',a,b)
def G(x,xp,z,w,wp,k,Lam,V):
    A,B,C=s(x-z),s(x-w),s(z-w); Ap,Bp,Cp=s(xp-z),s(xp-wp),s(z-wp)
    T=their_coef(A,B,C,Ap,Bp,Cp,k,Lam,V); t1,t2,t3=taus(x,z,w); u1,u2,u3=taus(xp,z,wp)
    return ( T['23']*(con(t2,u3)+con(t3,u2)) + T['33']*con(t3,u3) + T['22']*con(t2,u2)
           + T['13']*(con(t1,u3)+con(t3,u1)) + T['12']*(con(t1,u2)+con(t2,u1)) + T['11_full'] )
base=[np.array([0.31,-0.77]),np.array([-0.52,0.19]),np.array([0.05,-0.23]),
      np.array([1.13,0.42]),np.array([-0.31,0.88])]
x,xp,z0,w,wp=base; k=1.0; u=lambda t:np.array([np.cos(t),np.sin(t)])
B,Bp=s(x-w),s(xp-wp); Kww=((x-w)@(xp-wp))/(B*Bp)
print("="*80); print("1.  THE |z| -> infinity LOG  (the soft / rapidity region)"); print("="*80)
print("   K_ww' = (x-w).(x'-w')/[(x-w)^2(x'-w')^2] = %.9f"%Kww)
print("   %-10s %-10s %16s %16s"%("Lambda","V","<G>|z|^2 at 1e5","K_ww'/(4k+) log[(V-k+)/Lam]"))
for Lam,V in ((1e-4,1e5),(1e-6,1e5),(1e-6,1e8),(1e-8,1e6)):
    acc=0.0; n=3000
    for t in np.linspace(0,2*np.pi,n,endpoint=False): acc+=G(x,xp,1e5*u(t),w,wp,k,Lam,V)
    print("   %-10.0e %-10.0e %16.7f %16.7f"%(Lam,V,acc/n*1e10,Kww/(4*k)*np.log((V-k)/Lam)))
print("   -> the ONLY transverse log, and its coefficient is (LO structure) x (rapidity span).")
print()
print("="*80); print("2.  THE |x| -> infinity TAIL  (the source, not a loop divergence)"); print("="*80)
for lab,i in (("x",0),("x'",1)):
    print("   %-4s : "%lab,end="")
    for R in (1e3,1e4,1e5):
        acc=0.0; n=4000
        for t in np.linspace(0,2*np.pi,n,endpoint=False):
            pt=[b.copy() for b in base]; pt[i]=R*u(t); acc+=G(*pt,k,1e-6,1e6)
        print(" R=%-6.0e <G>R^2=%10.5f "%(R,acc/n*R*R),end="")
    print()
print("   log-divergent in the bare kernel, but x,x' are the arguments of rho: the")
print("   colour-charge correlator has the target's transverse support, so these ends")
print("   are cut at the target radius.  Not a divergence of the loop.")
print()
print("="*80); print("3.  {z,w,w'} COLLAPSE, REDONE IN HIGH PRECISION"); print("="*80)
mp.mp.dps=50
def Ghi(x,xp,z,w,wp,k,Lam,V):
    f=lambda v: mp.matrix([mp.mpf(float(v[0])),mp.mpf(float(v[1]))])
    S=lambda a:(a[0]**2+a[1]**2); D=lambda a,b:a[0]*b[0]+a[1]*b[1]
    X,Xp,Z,W,Wp=map(f,(x,xp,z,w,wp)); sub=lambda a,b: mp.matrix([a[0]-b[0],a[1]-b[1]])
    A,Bv,C=S(sub(X,Z)),S(sub(X,W)),S(sub(Z,W)); Ap,Bpv,Cp=S(sub(Xp,Z)),S(sub(Xp,Wp)),S(sub(Z,Wp))
    k,Lam,V=mp.mpf(k),mp.mpf(Lam),mp.mpf(V); L=mp.log; Dl=A*Bpv-Ap*Bv
    LA=L(((V-k)*A+k*Bv)/(Lam*A+k*Bv)); LAp=L(((V-k)*Ap+k*Bpv)/(Lam*Ap+k*Bpv))
    M23=(1/k)/Dl*( L(((V-k)*A+k*Bv)/((V-k)*Ap+k*Bpv)) - L((Lam*A+k*Bv)/(Lam*Ap+k*Bpv)) )
    M33=( (1/k)*L((V-k)/Lam)/(Bpv*A*Bv-Bv*Bv*Ap)*(A-Bv*Ap/Bpv) - (1/k)*A/(Bpv*A*Bv-Bv*Bv*Ap)*LA
         +(1/k)*Ap/(Bpv*Bpv*A-Bv*Ap*Bpv)*LAp )
    M22=( (1/k)*(Bpv/Ap)/Dl*LAp - (1/k)*(Bv/A)/Dl*LA )
    M13=( -(1/k)*L(V/(k+Lam))/((A-Bv)*(Bpv-Ap)) + (1/k)*A/((A-Bv)*Dl)*LA + (1/k)*Ap/((Bpv-Ap)*Dl)*LAp )
    M12=(  (1/k)*L(V/(Lam+k))/((A-Bv)*(Bpv-Ap)) - (1/k)*Bv/((A-Bv)*Dl)*LA - (1/k)*Bpv/((Bpv-Ap)*Dl)*LAp )
    F11=( (mp.mpf(2)/4)*(1/k)*L(V/(Lam+k))*(Bpv*Bv-A*Ap)/(C*Cp*(A-Bv)*(Ap-Bpv))
         -(mp.mpf(2)/4)*(1/k)*(A*Bv)/(C*Cp*(A-Bv))*(Ap-Bpv)/Dl*LA
         +(mp.mpf(2)/4)*(1/V-1/(k+Lam))/(C*Cp)
         +(mp.mpf(2)/4)*(1/k)*(Bpv*Ap)/(C*Cp)*(A-Bv)/((Ap-Bpv)*Dl)*LAp )
    t1,t2,t3=taus(x,z,w); u1,u2,u3=taus(xp,z,wp)
    return float(M23)*(con(t2,u3)+con(t3,u2))+float(M33)*con(t3,u3)+float(M22)*con(t2,u2)\
          +float(M13)*(con(t1,u3)+con(t3,u1))+float(M12)*(con(t1,u2)+con(t2,u1))+float(F11)
print("   %-9s %15s %17s %17s"%("rho","<|G|> float64","<|G|> 50-digit","x rho^4 (50-dig)"))
for rho in (1e-2,1e-3,1e-4,1e-5,1e-6):
    a64=ah=0.0; n=48
    for t in np.linspace(0,2*np.pi,n,endpoint=False):
        zz=w+rho*u(t); ww=w; wwp=w+rho*u(t+2.1)
        a64+=abs(G(x,xp,zz,ww,wwp,k,1e-6,1e6)); ah+=abs(Ghi(x,xp,zz,ww,wwp,k,1e-6,1e6))
    print("   %-9.0e %15.6e %17.6e %17.6e"%(rho,a64/n,ah/n,ah/n*rho**4))
print("   <|G|> x rho^4 -> 0 : CONVERGENT.  (float64 fails below rho~1e-3 here because")
print("   z,w,w' together drive AB'-A'B to zero -- the cancellation warned about earlier.)")
