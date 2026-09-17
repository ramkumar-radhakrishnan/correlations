"""Row 14: large-|w| tail -- cumulative radial integral, vectorised, fully resolved."""
import numpy as np
K=1.0; P=0.37; kt=np.array([0.7,-0.4]); kn=np.sqrt(kt@kt)
x=np.array([0.31,-0.77]); xp=np.array([-0.52,0.19]); yp=np.array([0.62,0.41]); z=np.array([0.05,-0.23])
s=lambda v:(v*v).sum(-1); Ss=K+P
A=s(x-z); kap=(kt@(yp-z))/K
Ki=(xp-yp)/s(xp-yp); Kj=(yp-z)/s(yp-z)
def Favg(R,nth):
    th=np.linspace(0,2*np.pi,nth,endpoint=False)
    W=R*np.stack([np.cos(th),np.sin(th)],-1)
    xw=x-W; zw=z-W; B=s(xw); C=s(zw)
    o=lambda a,b: a[...,:,None]*b[...,None,:]
    xz=np.broadcast_to(x-z,W.shape)
    Km=o(zw,xz)/C[:,None,None] - o(xw,xz)/(2*B[:,None,None])
    Pm=o(xw,zw)/C[:,None,None] + o(xw,xz)/(2*A)
    M=np.eye(2)*((A-B)/(2*Ss*C))[:,None,None] + Km/K + Pm/P
    T=(P/Ss)*np.eye(2)*np.trace(M,axis1=1,axis2=2)[:,None,None] - M - (P/K)*M.transpose(0,2,1)
    ph=np.exp(-1j*(W@kt*-1 + (kt@yp)))*np.exp(-1j*kap*P)   # e^{-ik(y'-w)}
    return (ph*np.einsum('nkm,k,m->n',T,Ki,Kj)/(P*A+K*B)).mean()*2*np.pi
print("="*76); print("CUMULATIVE  int_{10<|w|<R} d^2w F   (vectorised, >=400 pts per oscillation)")
print("="*76)
Rs=np.geomspace(10.0,3000.0,600); tot=0j
marks={99,199,299,399,499,599}
print("   %-9s %30s %14s"%("R","cumulative","|increment|"))
prev=0j
for i in range(len(Rs)-1):
    a,b=Rs[i],Rs[i+1]; Rm=0.5*(a+b)
    nth=int(min(max(4096,400*kn*Rm/(2*np.pi)),4_000_000))
    tot+=Favg(Rm,nth)*Rm*(b-a)
    if i in marks:
        print("   %-9.0f %30s %14.3e"%(b,"%.6f%+.6fj"%(tot.real,tot.imag),abs(tot-prev)))
        prev=tot
print()
print("   radial integrand ~ R^{-1/2} x cos(|k|R - phase)  (Fresnel-type):")
print("   %-9s %16s"%("R","R x |<F>|"))
for R in (100.,300.,1000.,3000.):
    nth=int(min(max(4096,400*kn*R/(2*np.pi)),4_000_000))
    print("   %-9.0f %16.6e"%(R,abs(Favg(R,nth))*R))
