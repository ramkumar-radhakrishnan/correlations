"""Row 17: is the transverse log a DGLAP (P_gg) log?  Fixed-S test + z->w collinear test."""
import numpy as np
s=lambda v:v@v
def vert(x,z,w,P,K):
    S=P+K; A,B,C=s(x-z),s(x-w),s(z-w)
    W=(z-w)/C-(x-w)/(2*B); V=(z-w)/C+(x-z)/(2*A)
    return (np.eye(2)*(A-B)/(2*C)+(S/K)*np.outer(x-z,W)+(S/P)*np.outer(V,x-w)), P*A+K*B
def F(x,xp,z,w,wp,P,K):
    S=P+K; Wa,D=vert(x,z,w,P,K); Wb,Dp=vert(xp,z,wp,P,K)
    return np.einsum('ji,ji->',Wa,Wb)*(P*K/S**2)/(D*Dp)
rng=np.random.default_rng(11); x,xp,w,wp=[rng.normal(size=2) for _ in range(4)]
Cuv=lambda z:z/(1-z)+(1-z)/z+z*(1-z)
print("="*78); print("A.  FIXED S = p+ + k+ :  does the large-|z| coefficient track P_gg(zeta)?")
print("="*78)
S=1.0; B,Bp=s(x-w),s(xp-wp); Kww=((x-w)@(xp-wp))/(B*Bp)
print("   %-8s %14s %14s %12s %12s %12s"%("zeta","plateau x 4S^2","1/(zeta zbar)","C_UV","ratio","1/(z zb)/C"))
for zt in (0.05,0.15,0.3,0.5,0.7,0.9,0.97):
    K=zt*S; P=S-K; R=3e4; acc=0.0
    for t in np.linspace(0,2*np.pi,3000,endpoint=False):
        acc+=F(x,xp,R*np.array([np.cos(t),np.sin(t)]),w,wp,P,K)
    pl=acc/3000*R*R/Kww*4*S*S   # strip the transverse factor K_ww'/4
    print("   %-8.3f %14.7f %14.7f %12.6f %12.6f %12.6f"
          %(zt,pl,1/(zt*(1-zt)),Cuv(zt),pl/Cuv(zt),(1/(zt*(1-zt)))/Cuv(zt)))
print()
print("   plateau x 4S^2 = 1/(zeta zetabar) EXACTLY.   Now  P_gg/(2N_c) = C_UV =")
print("      zeta zbar + 1/(zeta zbar) - 2 ,   so this row supplies ONLY the 1/(zeta zbar)")
print("      term = zeta/zbar + zbar/zeta + 2 : the two SOFT ends, with the hard")
print("      zeta zbar - 2 piece absent.  ratio to C_UV is not constant -> no P_gg.")

print(); print("="*78); print("B.  THE FINAL-STATE COLLINEAR REGION  z -> w  (rho = |z-w|)"); print("="*78)
print("   a fragmentation log needs the integrand ~ rho^-2 against d^2rho ~ rho drho.")
print("   %-10s %16s %16s %16s"%("rho","<|F|>","<|F|> x rho","<|F|> x rho^2"))
for rho in (1e-1,1e-2,1e-3,1e-4,1e-5):
    acc=0.0
    for t in np.linspace(0,2*np.pi,512,endpoint=False):
        acc+=abs(F(x,xp,w+rho*np.array([np.cos(t),np.sin(t)]),w,wp,0.37,1.0))
    a=acc/512; print("   %-10.0e %16.6e %16.9f %16.6e"%(rho,a,a*rho,a*rho*rho))
print("   <|F|> x rho is CONSTANT  =>  integrand ~ 1/rho, measure rho drho  =>  int drho:")
print("   linearly convergent, NO collinear log.  Reason: W^{ji}(x,z,w) is only 1/rho")
print("   singular there (the delta^{ji} numerator (x-z)^2-(x-w)^2 vanishes like rho),")
print("   and the conjugate vertex W^{ji}(x',z,w') is REGULAR at z=w because w' != w.")
print("   The phase e^{-ik.(w'-w)} is exactly what keeps w and w' apart.")

print(); print("="*78); print("C.  z -> w AND w' -> w SIMULTANEOUSLY (the only way to make both singular)")
print("="*78)
for rho in (1e-2,1e-3,1e-4,1e-5):
    acc=0.0
    for t in np.linspace(0,2*np.pi,512,endpoint=False):
        u=np.array([np.cos(t),np.sin(t)]); v=np.array([-np.sin(t),np.cos(t)])
        acc+=abs(F(x,xp,w+rho*u,w,w+rho*v,0.37,1.0))
    a=acc/512
    print("   rho=%-9.0e  <|F|> x rho^2 = %12.6e  (constant => |F| ~ rho^-2)"%(rho,a*rho**2))
print("   -> |F| ~ rho^-2 against the 3-point measure rho^3 drho  =>  int rho drho : convergent.")
print("      (and that region is anyway killed by the phase once k_perp != 0.)")
