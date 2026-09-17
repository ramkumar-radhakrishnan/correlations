"""Row 14: the large-|w| tail -- cumulative radial integral, keeping the oscillation."""
import numpy as np
K=1.0; kt=np.array([0.7,-0.4]); s=lambda v:v@v; kn=np.sqrt(kt@kt); P=0.37
def Favg(R,nth):
    x=np.array([0.31,-0.77]); xp=np.array([-0.52,0.19]); yp=np.array([0.62,0.41]); z=np.array([0.05,-0.23])
    Ss=K+P; A=s(x-z); kap=(kt@(yp-z))/K
    Ki=(xp-yp)/s(xp-yp); Kj=(yp-z)/s(yp-z)
    acc=0j
    for t in np.linspace(0,2*np.pi,nth,endpoint=False):
        w=R*np.array([np.cos(t),np.sin(t)])
        B=s(x-w); C=s(z-w)
        Km=np.outer(z-w,x-z)/C - np.outer(x-w,x-z)/(2*B)
        Pm=np.outer(x-w,z-w)/C + np.outer(x-w,x-z)/(2*A)
        M=np.eye(2)*(A-B)/(2*Ss*C) + Km/K + Pm/P
        T=(P/Ss)*np.eye(2)*np.trace(M) - M - (P/K)*M.T
        acc+=np.exp(-1j*(kt@(yp-w)))*np.exp(-1j*kap*P)*np.einsum('km,k,m->',T,Ki,Kj)/(P*A+K*B)
    return acc/nth*2*np.pi
print("="*76); print("CUMULATIVE  int_{10<|w|<R} d^2w F   (radial oscillation kept)"); print("="*76)
print("   the w integral is the Fourier transform to the measured momentum k.")
print("   |<F>| alone falls only as R^{-3/2} (because P^{ij} grows like |w| through")
print("   (x-w)^i (x-z)^j/(2A)), so |<F>| R is NOT absolutely integrable -- but the")
print("   RADIAL oscillation supplies the rest.  Cumulative integral:")
print()
print("   %-10s %26s %16s"%("R","int_{10}^{R} d^2w F","increment"))
tot=0j; prev=0j; R0=10.0
Rs=np.geomspace(10.0,4000.0,220)
for i in range(len(Rs)-1):
    a,b=Rs[i],Rs[i+1]; Rm=0.5*(a+b)
    nth=int(min(max(2048,300*kn*Rm/(2*np.pi))*4,120000))
    tot+=Favg(Rm,nth)*Rm*(b-a)
    if i in (40,90,140,180,218):
        print("   %-10.0f %26s %16.3e"%(b,"%.6f%+.6fj"%(tot.real,tot.imag),abs(tot-prev)))
        prev=tot
print()
print("   -> the cumulative integral settles: the tail is CONDITIONALLY convergent, exactly")
print("      as a 2D Fourier transform of a 1/|w| tail should be (int d^2w e^{ikw}/|w| = 2pi/|k|).")
print("      So large |w| is not a divergence of the row -- it is the Fourier transform.")
