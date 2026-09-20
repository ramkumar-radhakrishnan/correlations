"""Row 21: is 2b the only rapidity term, what does its coefficient simplify to,
   and does it build the full JIMWLK kernel?"""
import numpy as np, sympy as sp
from row19_lib import their_coef
s=lambda v:v@v
print("="*80); print("1.  THE BLUE REWRITING"); print("="*80)
V,k,Lam=sp.symbols('vee kplus Lambda',positive=True)
lhs=sp.log((V-k)/Lam); rhs=sp.log(V/Lam)+sp.log(1-k/V)
print("   log[(V-k+)/Lam] - { log[V/Lam] + log[1-k+/V] }  =  %s   EXACT"%sp.simplify(sp.expand_log(lhs-rhs,force=True)))
print("   log(1-k+/V) is finite (-> 0 for V >> k+): the rapidity log is log[V/Lambda].")

print(); print("="*80); print("2.  ONLY 2b CARRIES log(1/Lambda) -- AND ONLY ITS FIRST PIECE"); print("="*80)
rng=np.random.default_rng(55); x,xp,z,w,wp=[rng.normal(size=2)*1.3 for _ in range(5)]
A,B,C=s(x-z),s(x-w),s(z-w); Ap,Bp,Cp=s(xp-z),s(xp-wp),s(z-wp); k_=1.0
print("   %-8s %16s %16s %14s  [predicted]"%("part","Lam=1e-5","Lam=1e-9","d/dlog(1/Lam)"))
for nm,pr in (('23',0),('33',1/(k_*B*Bp)),('22',0),('13',0),('12',0),('11_full',0)):
    v=[their_coef(A,B,C,Ap,Bp,Cp,k_,l,1e6)[nm] for l in (1e-5,1e-9)]
    print("   %-8s %16.9f %16.9f %14.9f  [%.9f]"%(nm,*v,(v[1]-v[0])/np.log(1e4),pr))
L=np.log
p2=lambda Lam: -(1/k_)*A/(Bp*A*B-B*B*Ap)*L(((1e6-k_)*A+k_*B)/(Lam*A+k_*B))
p3=lambda Lam: (1/k_)*Ap/(Bp*Bp*A-B*Ap*Bp)*L(((1e6-k_)*Ap+k_*Bp)/(Lam*Ap+k_*Bp))
print("   2b piece 2 : Lam=1e-5 -> %.9f   Lam=1e-9 -> %.9f   (Lambda-FINITE)"%(p2(1e-5),p2(1e-9)))
print("   2b piece 3 : Lam=1e-5 -> %.9f   Lam=1e-9 -> %.9f   (Lambda-FINITE)"%(p3(1e-5),p3(1e-9)))
print("   -> yes: within THIS row, the blue piece is the only source of log(1/Lambda).")

print(); print("="*80); print("3.  YOUR PREFACTOR COLLAPSES TO 1/[(x-w)^2 (x'-w')^2]"); print("="*80)
Asy,Bsy,Apsy,Bpsy=sp.symbols("A B Ap Bp",positive=True)
pref=1/(Bpsy*Asy*Bsy - Bsy**2*Apsy)*(Asy - Bsy*Apsy/Bpsy)
print("   1/[B'AB - B^2A'] x [A - B A'/B']  =  %s"%sp.simplify(pref))
print("   numerically: %.12f   vs   1/(B B') = %.12f"%(
      1/(Bp*A*B-B*B*Ap)*(A-B*Ap/Bp), 1/(B*Bp)))

print(); print("="*80); print("4.  AND THE BRACKET FACTORISES"); print("="*80)
Vz =(z-w)/C  + (x-z)/(2*A)
Vzp=(z-wp)/Cp+ (xp-z)/(2*Ap)
brk=( ((xp-wp)@(x-w))*((z-w)@(z-wp))/(C*Cp) + ((xp-wp)@(x-w))*((xp-z)@(z-w))/(2*C*Ap)
     +((x-w)@(xp-wp))*((x-z)@(z-wp))/(2*A*Cp) + ((x-z)@(xp-z))*((xp-wp)@(x-w))/(4*A*Ap) )
print("   your 4-term bracket                      = %.12f"%brk)
print("   [(x-w).(x'-w')] x [V_z . V_z']           = %.12f"%(((x-w)@(xp-wp))*(Vz@Vzp)))
print("   with  V_z = (z-w)/(z-w)^2 + (x-z)/2(x-z)^2 = K_zw - (1/2) K_zx")
print()
print("   => the whole rapidity term is")
print("      (1/(2pi)^3)(g^4/8pi^5)(1/k+) log[V/Lam] int e^{-ik(w'-w)}")
print("        x  (x-w).(x'-w')/[(x-w)^2 (x'-w')^2]  x  [K_zw - K_zx/2].[K_zw' - K_zx'/2]  x  colour")

print(); print("="*80); print("5.  IS IT THE FULL JIMWLK KERNEL?  NO -- THE 1/2 GIVES IT AWAY"); print("="*80)
K=lambda a,b:(a-b)/s(a-b)
u=lambda t:np.array([np.cos(t),np.sin(t)])
print("   large-|z| behaviour of the z-kernel product, angle-averaged:")
print("   %-8s %20s %20s"%("|z|","yours: (K-K/2).(K-K/2)","full BK: (K-K).(K-K)"))
for R in (1e2,1e3,1e4,1e5):
    a1=a2=0.0; n=4000
    for t in np.linspace(0,2*np.pi,n,endpoint=False):
        zz=R*u(t)
        v1=(K(zz,w)-K(zz,x)/2)@(K(zz,wp)-K(zz,xp)/2)
        v2=(K(zz,w)-K(zz,x))  @(K(zz,wp)-K(zz,xp))
        a1+=v1; a2+=v2
    print("   %-8.0e %20.9f %20.9f"%(R,a1/n*R*R,a2/n*R*R))
print("   yours x |z|^2 -> 1/4 (constant): LOG DIVERGENT.")
print("   full BK x |z|^2 -> 0          : CONVERGENT (it falls like |z|^-4).")
print()
print("   That is exactly the |z| -> infinity IR log of this row.  The full dipole kernel")
print("   K_zw - K_zx has NO large-|z| log; the half-weighted one does.  So this row")
print("   supplies only part of the evolution kernel -- the missing halves come from the")
print("   partner rows, and the z-log closes at the same time as the kernel completes.")
