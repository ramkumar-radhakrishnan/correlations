"""Row 25: quantify the large-|z| log, the phase cutoff, and the p+ -> 0 behaviour."""
import numpy as np
s=lambda v:v@v; dd=np.eye(2); KT=np.array([0.6,-0.35]); KP=1.0
def groups(x,z,w,P,K,dperp=2.0):
    """the six groups of the bracket, separately, as 2x2 matrices R[m,k']"""
    S=P+K; A,B,C=s(x-z),s(x-w),s(z-w); xz,xw,zw=x-z,x-w,z-w
    return {
    'd_perp'    : dd*dperp*P/(2*C*S**2)*(A-B),
    'S/k+'      : dd*(P/(K*S))*((xz@zw)/C-(xz@xw)/(2*B)-(A-B)/(2*C)),
    'S/p+'      : dd*(1/S)*((xw@zw)/C+(xz@xw)/(2*A)-(A-B)/(2*C)),
    '-1/p+'     : -(1/P)*( np.outer(zw,xw)/C + np.outer(xz,xw)/(2*A) ),
    '-p+/k+^2'  : -(P/K**2)*( np.outer(zw,xz)/C - np.outer(xw,xz)/(2*B) ),
    '-1/k+'     : -(1/K)*( np.outer(xz,zw)/C - np.outer(xz,xw)/(2*B)
                          +np.outer(xw,zw)/C + np.outer(xw,xz)/(2*A) )}
x =np.array([0.41,-0.63]); xp=np.array([-0.58,0.22]); yp=np.array([0.09,0.51]); w=np.array([-0.25,0.86])
u=lambda t:np.array([np.cos(t),np.sin(t)])
def piece(nm,zz,P,K=KP,phase=False):
    N=(yp-zz)/s(yp-zz); X=(xp-yp)/s(xp-yp); D=P*s(x-zz)+K*s(x-w)
    v=np.einsum('m,k,mk->',N,X,groups(x,zz,w,P,K)[nm])/D
    return v*(np.exp(-1j*(P/K)*(KT@(yp-zz))) if phase else 1.0)
print("="*80); print("1.  WHICH GROUP CARRIES THE LARGE-|z| LOG  (phase off)"); print("="*80)
print("   <group> x |z|^2 , angle averaged:")
print("   %-12s %14s %14s %14s   %s"%("group","|z|=1e2","1e3","1e4","behaviour"))
for nm in ('d_perp','S/k+','S/p+','-1/p+','-p+/k+^2','-1/k+'):
    v=[]
    for R in (1e2,1e3,1e4):
        acc=0.0
        for t in np.linspace(0,2*np.pi,4096,endpoint=False): acc+=piece(nm,R*u(t),0.37)
        v.append(acc/4096*R*R)
    flat=abs(v[2])>1e-10 and abs(v[2])>0.25*abs(v[0])
    print("   %-12s %14.6e %14.6e %14.6e   %s"%(nm,*v,"1/|z|^2  -> LOG" if flat else "falls faster"))
print("   => the log comes from the -p+/k+^2 and -1/k+ groups only.  Crucially the -1/p+ group,")
print("      the one that is soft-enhanced, is NOT log divergent in z.  So no double log.")
print()
print("   p+ dependence of that log's coefficient (sum of all groups, phase off, |z|=1e4):")
print("   %-8s %16s %16s"%("p+","coefficient","x k+/1"))
for P in (0.05,0.2,0.5,1.0,3.0):
    acc=0.0
    for t in np.linspace(0,2*np.pi,4096,endpoint=False):
        acc+=sum(piece(nm,1e4*u(t),P) for nm in ('d_perp','S/k+','S/p+','-1/p+','-p+/k+^2','-1/k+'))
    print("   %-8.2f %16.9f %16.9f"%(P,acc/4096*1e8,acc/4096*1e8*P))
print("   -> tends to a CONSTANT as p+ -> 0 (not 1/p+): the z-log is not soft-enhanced.")

print(); print("="*80); print("2.  THE PHASE CUTS THE z INTEGRAL AT |z| ~ k+/(p+ |k|)"); print("="*80)
def zint(P,R,nr=900,nth=1200):
    """radial+angular integral of the full integrand over |z| < R, phase ON"""
    tot=0.0+0j
    lo=np.log(1e-2); hi=np.log(R)
    for lr in (np.arange(nr)+0.5)*(hi-lo)/nr+lo:
        rr=np.exp(lr); acc=0j
        for t in (np.arange(nth)+0.5)*(2*np.pi/nth):
            acc+=sum(piece(nm,rr*u(t),P,phase=True) for nm in
                     ('d_perp','S/k+','S/p+','-1/p+','-p+/k+^2','-1/k+'))
        tot+=acc/nth*2*np.pi*rr*rr*((hi-lo)/nr)
    return tot
print("   int_{|z|<R} d^2z (integrand) , phase ON :")
print("   %-8s %12s %14s %14s %14s   %s"%("p+","k+/(p+|k|)","R=1e2","R=1e3","R=1e4","saturates?"))
for P in (0.2,1.0,4.0):
    v=[abs(zint(P,R)) for R in (1e2,1e3,1e4)]
    print("   %-8.2f %12.1f %14.7f %14.7f %14.7f   %s"%(P,KP/(P*np.sqrt(s(KT))),*v,
          "yes" if abs(v[2]-v[1])<0.06*abs(v[1]) else "still growing"))
print("   -> the integral saturates once R exceeds ~ k+/(p+|k_perp|): the phase")
print("      e^{+i (p+/k+) k.z} is what stops the would-be collinear log.")

print(); print("="*80); print("3.  THE p+ -> 0 BEHAVIOUR AT FIXED TRANSVERSE POINTS"); print("="*80)
z0=np.array([1.07,-0.34])
print("   %-10s %18s %18s"%("p+","integrand","x p+"))
for P in (1e-1,1e-2,1e-3,1e-4,1e-5):
    v=sum(piece(nm,z0,P,phase=True) for nm in ('d_perp','S/k+','S/p+','-1/p+','-p+/k+^2','-1/k+'))
    print("   %-10.0e %18.9f %18.9f"%(P,v.real,v.real*P))
print("   -> integrand x p+ tends to a constant: a simple 1/p+ pole, i.e. ONE rapidity log.")
