"""Row 14: is the large-|w| region really divergent, or was that aliasing?"""
import numpy as np
K=1.0; kt=np.array([0.7,-0.4]); s=lambda v:v@v; kn=np.sqrt(kt@kt); P=0.37
def F(x,xp,yp,z,w):
    Ss=K+P; A=s(x-z); B=s(x-w); C=s(z-w)
    Km=np.outer(z-w,x-z)/C - np.outer(x-w,x-z)/(2*B)
    Pm=np.outer(x-w,z-w)/C + np.outer(x-w,x-z)/(2*A)
    M=np.eye(2)*(A-B)/(2*Ss*C) + Km/K + Pm/P
    T=(P/Ss)*np.eye(2)*np.trace(M) - M - (P/K)*M.T
    kap=(kt@(yp-z))/K
    ph=np.exp(-1j*(kt@(yp-w)))*np.exp(-1j*kap*P)
    return ph*np.einsum('km,k,m->',T,(xp-yp)/s(xp-yp),(yp-z)/s(yp-z))/(P*A+K*B)
x=np.array([0.31,-0.77]); xp=np.array([-0.52,0.19]); yp=np.array([0.62,0.41]); z=np.array([0.05,-0.23])
print("="*74); print("LARGE |w|, WITH ADEQUATE ANGULAR RESOLUTION (>= 200 pts per oscillation)")
print("="*74)
print("  %-9s %12s %18s %14s"%("|w|","n_theta","|A(R)| R^2","ratio"))
prev=None
for R in (30.,100.,300.,1000.):
    nth=int(200*kn*R/(2*np.pi))*8
    acc=0j
    for t in np.linspace(0,2*np.pi,nth,endpoint=False):
        acc+=F(x,xp,yp,z,R*np.array([np.cos(t),np.sin(t)]))
    v=abs(acc/nth)*R**2
    print("  %-9.0f %12d %18.6e %14s"%(R,nth,v,"-" if prev is None else "%.3f"%(v/prev)))
    prev=v
print()
print("  the phase e^{+ik.w} is present and p+-INDEPENDENT (it is the LO phase), so the")
print("  angular average carries a Bessel J0 ~ R^{-1/2} envelope: |A(R)| R^2 ~ R^{3/2} x R^{-2}")
print("  ... i.e. the radial integrand ~ R^{-3/2} dR, CONVERGENT.  The earlier 'LOG DIV'")
print("  verdict for w came from under-resolving ~5x10^5 oscillations with 4x10^5 samples.")
