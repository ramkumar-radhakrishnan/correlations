"""Row 25 (vectorised): the large-|z| log, the phase cutoff, and p+ -> 0."""
import numpy as np
s=lambda v:v@v; KT=np.array([0.6,-0.35]); KP=1.0; kabs=np.sqrt(KT@KT)
x =np.array([0.41,-0.63]); xp=np.array([-0.58,0.22]); yp=np.array([0.09,0.51]); w=np.array([-0.25,0.86])
X =(xp-yp)/s(xp-yp)
def val(Z,P,K=KP,dperp=2.0,phase=True,only=None):
    """Z : (n,2) array of z points.  Returns the integrand at each."""
    zx,zy=Z[:,0],Z[:,1]
    xz=np.stack([x[0]-zx,x[1]-zy],1); xw=(x-w)[None,:]+0*xz
    zw=np.stack([zx-w[0],zy-w[1]],1); yz=np.stack([yp[0]-zx,yp[1]-zy],1)
    A=(xz**2).sum(1); B=s(x-w); C=(zw**2).sum(1); Dv=P*A+K*B
    N=yz/((yz**2).sum(1))[:,None]
    S=P+K
    NX=(N*X).sum(1)
    def out(a,b): return (a*N).sum(1)*(b*X).sum(1)     # contracts R=outer(a,b) with N^m X^k'
    g={}
    g['d_perp']  = NX*dperp*P/(2*C*S**2)*(A-B)
    g['S/k+']    = NX*(P/(K*S))*(((xz*zw).sum(1))/C-((xz*xw).sum(1))/(2*B)-(A-B)/(2*C))
    g['S/p+']    = NX*(1/S)*(((xw*zw).sum(1))/C+((xz*xw).sum(1))/(2*A)-(A-B)/(2*C))
    g['-1/p+']   = -(1/P)*( out(zw,xw)/C + out(xz,xw)/(2*A) )
    g['-p+/k+^2']= -(P/K**2)*( out(zw,xz)/C - out(xw,xz)/(2*B) )
    g['-1/k+']   = -(1/K)*( out(xz,zw)/C - out(xz,xw)/(2*B) + out(xw,zw)/C + out(xw,xz)/(2*A) )
    tot=g[only] if only else sum(g.values())
    tot=tot/Dv
    if phase: tot=tot*np.exp(-1j*(P/K)*(yz@KT))
    return tot
ring=lambda R,n: R*np.stack([np.cos((np.arange(n)+.5)*2*np.pi/n),np.sin((np.arange(n)+.5)*2*np.pi/n)],1)
print("="*80); print("1.  WHICH GROUP CARRIES THE LARGE-|z| LOG  (phase off)"); print("="*80,flush=True)
print("   <group> x |z|^2 , angle averaged")
print("   %-12s %14s %14s %14s   %s"%("group","|z|=1e2","1e3","1e4","behaviour"))
for nm in ('d_perp','S/k+','S/p+','-1/p+','-p+/k+^2','-1/k+'):
    v=[val(ring(R,20000),0.37,phase=False,only=nm).mean().real*R*R for R in (1e2,1e3,1e4)]
    flat=abs(v[2])>1e-10 and abs(v[2])>0.25*abs(v[0])
    print("   %-12s %14.6e %14.6e %14.6e   %s"%(nm,*v,"~ 1/|z|^2  -> LOG" if flat else "falls faster"))
print("   => the S/k+, -p+/k+^2 and -1/k+ groups. (The -1/p+ group does not, but the 1/D")
print("      factor supplies the 1/p+ anyway -- see row25_soft.py.)")
print()
print("   p+ dependence of the total coefficient (phase off, |z|=1e4):")
print("   %-8s %18s %18s"%("p+","coefficient","x p+"))
for P in (0.02,0.05,0.2,1.0,3.0):
    c=val(ring(1e4,20000),P,phase=False).mean().real*1e8
    print("   %-8.2f %18.9f %18.9f"%(P,c,c*P))
print("   -> coefficient x p+ is constant at small p+, i.e. the z-log coefficient goes like 1/p+.",flush=True)

print(); print("="*80); print("2.  THE PHASE CUTS THE z INTEGRAL AT |z| ~ k+/(p+ |k_perp|)"); print("="*80,flush=True)
def zint(P,R,nr=600,nth=8000):
    lo,hi=np.log(1e-2),np.log(R); lr=(np.arange(nr)+.5)*(hi-lo)/nr+lo; rr=np.exp(lr)
    tot=0j
    for r0 in rr:
        tot+=val(ring(r0,nth),P).mean()*2*np.pi*r0*r0*((hi-lo)/nr)
    return tot
print("   %-7s %13s %13s %13s %13s   %s"%("p+","k+/(p+|k|)","R=1e2","R=1e3","R=1e4","saturated?"))
for P in (0.2,1.0,4.0):
    v=[abs(zint(P,R)) for R in (1e2,1e3,1e4)]
    print("   %-7.2f %13.1f %13.7f %13.7f %13.7f   %s"%(P,KP/(P*kabs),*v,
          "yes" if abs(v[2]-v[1])<0.06*max(abs(v[1]),1e-12) else "still growing"),flush=True)
print("   -> saturates once R exceeds ~ k+/(p+|k_perp|).  The phase e^{+i(p+/k+)k.z} is")
print("      what stops the would-be collinear log: it is cut off kinematically, not by a pole.")

print(); print("="*80); print("3.  p+ -> 0 AT FIXED TRANSVERSE POINTS"); print("="*80)
z0=np.array([[1.07,-0.34]])
print("   %-10s %18s %18s"%("p+","integrand","x p+"))
for P in (1e-1,1e-2,1e-3,1e-4,1e-5):
    v=val(z0,P)[0].real
    print("   %-10.0e %18.9f %18.9f"%(P,v,v*P))
print("   -> a simple 1/p+ pole: ONE rapidity logarithm, no 1/p+^2, no log^2.",flush=True)
