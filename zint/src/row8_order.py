"""Row 8: can the + prescription be done BEFORE any transverse integration?
Test 1: does p+ commute with z?     Test 2: is the g(0) subtraction still transverse-safe?
Test 3: what does keeping z un-integrated cost you?"""
import numpy as np, mpmath as mp
mp.mp.dps = 20
Kp=1.0; Pv=30.0; Lam=1e-9; kt=np.array([0.7,-0.4]); Rcut=300.0
def K(a,b): d=a-b; return d/(d@d)

x =np.array([0.31,-0.77]); xp=np.array([-0.52,0.19])
y =np.array([1.13, 0.42]); yp=np.array([-0.31,0.88])
r = yp-y; kap = (kt@r)/Kp
X, Xp = K(x,y), K(xp,yp)

# ---------------- the four p+ integrals, exact -------------------------------
def I0():   return complex(mp.quad(lambda t: mp.e**(-1j*kap*t), [Lam,Pv]))
def I1():   # PLUS PRESCRIPTION form
    return np.log(Pv/Lam) + complex(mp.quad(lambda t:(mp.e**(-1j*kap*t)-1)/t,[0,Pv]))
def I1_direct(): return complex(mp.quad(lambda t: mp.e**(-1j*kap*t)/t, [Lam,Pv]))
def I2():   return complex(mp.quad(lambda t: t*mp.e**(-1j*kap*t), [Lam,Pv]))
def I3():   return complex(mp.quad(lambda t: t*mp.e**(-1j*kap*t)/(Kp+t)**2, [Lam,Pv]))
i0,i1,i1d,i2,i3 = I0(),I1(),I1_direct(),I2(),I3()
print("="*76); print("0.  the + prescription identity itself, at FIXED transverse points"); print("="*76)
print("   direct   int_Lam^Pv dp+ e^{-i kap p+}/p+ = %s" % np.round(i1d,10))
print("   plus     log(Pv/Lam) + int_0^Pv[..-1]/p+ = %s" % np.round(i1,10))
print("   difference = %.3e   (the O(kappa*Lambda) remainder)" % abs(i1-i1d))

# ---------------- z quadrature ------------------------------------------------
nr,nth = 4000, 720
rr=np.geomspace(1e-7,Rcut,nr); th=np.linspace(0,2*np.pi,nth,endpoint=False)
RR,TH=np.meshgrid(rr,th,indexing='ij')
Z=np.stack([RR*np.cos(TH),RR*np.sin(TH)],-1)
W=RR*np.gradient(rr)[:,None]*(2*np.pi/nth)
U=y-Z; V=yp-Z; u2=(U**2).sum(-1); v2=(V**2).sum(-1)
M=U/u2[...,None]; N=V/v2[...,None]
dot=lambda a,b:(a*b).sum(-1)
MN, MX, MXp, NX, NXp = dot(M,N), M@X, M@Xp, N@X, N@Xp
XXp = X@Xp

# ---------------- ROUTE A : p+ FIRST (with + prescription), then z ------------
integrand_A = (MN*XXp*(i1 + i2/Kp**2)          # delta_mm' delta_kk'
               + 2*MX*NXp*i3                    # d_perp = 2 term
               + 2/Kp*(MXp*NX - MX*NXp)*i0)     # the 2/k+ term
A = (integrand_A*W).sum()

# ---------------- ROUTE B : z FIRST (closed form T), then p+ -----------------
r2=r@r
T = np.pi/2*(np.log(Rcut**2/r2)+1)*np.eye(2) - np.pi*np.outer(r,r)/r2
trT = np.trace(T)
B = (trT*XXp*(i1 + i2/Kp**2)
     + 2*(X@T@Xp)*i3
     + 2/Kp*(Xp@T@X - X@T@Xp)*i0)

print(); print("="*76); print("1.  DO THE TWO INTEGRATIONS COMMUTE?"); print("="*76)
print("   route A : p+ first (with the + prescription), then z  = %s" % np.round(A,8))
print("   route B : z first (closed form), then p+              = %s" % np.round(B,8))
print("   relative difference = %.3e   -> THEY COMMUTE" % (abs(A-B)/abs(B)))

print(); print("="*76); print("2.  WHY -- the p+/transverse coupling"); print("="*76)
print("   the ONLY place p+ meets a transverse variable is the phase")
print("       e^{-i k.(y'-y) p+/k+}      kappa = k.(y'-y)/k+ = %.6f" % kap)
print("   it involves y and y' ONLY.  z, x, x' never see p+ at all, so the whole")
print("   z, x, x' structure is a SPECTATOR factor that pulls out of the p+ integral.")
print("   contrast the earlier rows, where p+ sat inside the z phase and the")
print("   p+ -> 0 subtraction stripped it, exposing a bare 1/z^2 tail.")

print(); print("="*76); print("3.  IS THE g(0) SUBTRACTION TERM STILL TRANSVERSE-SAFE?"); print("="*76)
print("   at p+ = 0 the phase becomes e^{-i k.(y'-y)}, NOT 1: the LO oscillation survives.")
print("   large-|y| angular average x R^2, subtraction term (p+ = 0):")
for R in (1e2,1e3,1e4):
    nn=int(max(20000,60*np.sqrt(kt@kt)*R)); t=np.linspace(0,2*np.pi,nn,endpoint=False)
    Y=R*np.stack([np.cos(t),np.sin(t)],-1)
    Mx=Y-np.array([0.05,-0.23]); Mx=Mx/ (Mx**2).sum(-1)[:,None]
    Xx=(np.array([0.31,-0.77])-Y); Xx=Xx/(Xx**2).sum(-1)[:,None]
    ph=np.exp(1j*(Y@kt))                      # p+ = 0 : only the LO phase
    val=abs(((Mx*Xx).sum(-1)*ph).mean())*R**2
    print("      R = %-8.0e   |A| R^2 = %.5e" % (R,val))
print("   falls as R^{-1/2} (Bessel) -> convergent.  The subtraction is safe.")

print(); print("="*76); print("4.  THE PRICE OF LEAVING z UN-INTEGRATED"); print("="*76)
pointwise = 2/Kp*(MXp*NX - MX*NXp)
print("   the 2/k+ term, pointwise in z:  max |integrand| = %.5f  (NOT zero)" % np.abs(pointwise).max())
print("   after int d^2z                                  = %.3e  (zero)"
      % abs((pointwise*W).sum()))
print("   so if you keep z un-integrated you must carry a term that is only")
print("   zero under the integral -- harmless, but do not drop or approximate it")
print("   before integrating.")
print()
print("   and the z integral is still IR log divergent:  tr T = pi log(R^2/r^2) = %.4f at R = %.0f"
      % (trT, Rcut))
print("   -> the + prescription does NOT remove the need to regulate z.")
