"""Row 15: the large-|x| logarithm -- measure its zeta dependence and test for P_gg."""
import numpy as np
K=1.0; kt=np.array([0.7,-0.4]); s=lambda v:v@v
xp=np.array([-0.52,0.19]); yp=np.array([0.62,0.41]); z=np.array([0.05,-0.23]); w=np.array([1.13,0.42])
def F(x,P):
    Ss=K+P; A=s(x-z); B=s(x-w); C=s(z-w)
    Km=np.outer(z-w,x-z)/C - np.outer(x-w,x-z)/(2*B)
    Pm=np.outer(x-w,z-w)/C + np.outer(x-w,x-z)/(2*A)
    M=np.eye(2)*(A-B)/(2*Ss*C) + Km/K + Pm/P
    T=(P/Ss)*np.eye(2)*np.trace(M) - M - (P/K)*M.T
    kap=(kt@(yp-z))/K
    ph=np.exp(-1j*(kt@(yp-w)))*np.exp(-1j*kap*P)
    return ph*np.einsum('km,k,m->',T,(xp-yp)/s(xp-yp),(yp-z)/s(yp-z))/(P*A+K*B)
def plateau(P,R=3e4,nth=8192):
    acc=0j
    for t in np.linspace(0,2*np.pi,nth,endpoint=False):
        acc+=F(R*np.array([np.cos(t),np.sin(t)]),P)
    return (acc/nth*R**2)
Cuv=lambda zt:(zt/(1-zt) + (1-zt)/zt + zt*(1-zt))
print("="*76); print("LARGE-|x| PLATEAU vs THE SPLITTING VARIABLE  zeta = k+/(p+ + k+)"); print("="*76)
print("   (phase factored out; the transverse structure is p+-independent, so any")
print("    zeta dependence of the plateau is purely longitudinal)")
print()
ph0=np.exp(-1j*(kt@(yp-w)))
print("   %-8s %-9s %16s %14s %16s"%("p+","zeta","plateau/phase","C_UV(zeta)","ratio"))
vals=[]
for Pp in (0.15,0.37,0.7,1.5,3.0,6.0):
    zt=K/(K+Pp)
    kap=(kt@(yp-z))/K
    v=(plateau(Pp)/(ph0*np.exp(-1j*kap*Pp))).real
    c=Cuv(zt); vals.append((zt,v,c))
    print("   %-8.2f %-9.5f %16.8f %14.6f %16.6f"%(Pp,zt,v,c,v/c))
print()
print("   the ratio is NOT constant -> the large-|x| log is not proportional to P_gg(zeta).")
r=[v/c for _,v,c in vals]
print("   ratio spread: min %.5f  max %.5f  (a P_gg coefficient would give a flat column)"
      %(min(r),max(r)))
print()
print("   for comparison, what IS flat: divide instead by the pure soft pole 1/zetabar:")
print("   %-9s %16s %16s"%("zeta","plateau x zetabar","plateau x zeta zetabar"))
for zt,v,c in vals:
    print("   %-9.5f %16.8f %16.8f"%(zt,v*(1-zt),v*zt*(1-zt)))
