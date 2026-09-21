"""Row 28: A^(1)dag Cbar^dag x B_2 with C splitting the MEASURED gluon. Prefactor + tensor."""
import numpy as np, sympy as sp, itertools
g,pi=sp.symbols('g pi',positive=True); P,K=sp.symbols('pplus kplus',positive=True); I=sp.I
Dt=sp.Symbol('Dtilde',positive=True)          # Dt = p+(y-x)^2 + (k+-p+)(y-z)^2
print("="*80); print("1.  PREFACTOR"); print("="*80)
A1d=-I*g/(sp.sqrt(2)*pi*sp.sqrt(K))                       # A^(1)dag at k+
C1d=-g/(2*pi*sp.sqrt(2*K))*sp.sqrt(P*(K-P))/K             # Cbar^dag at k+_arg = k+
B2 = I*g**2*sp.sqrt(P*(K-P))/(4*pi**2*Dt*K)               # B_2(k+-p+,-z ; p+,-x), source y
print("   A^(1)dag = %s"%A1d)
print("   Cbar^dag = %s   (sqrt(p+(k+-p+))/k+)"%sp.simplify(C1d))
print("   B_2      = %s   (the two momenta are k+-p+ and p+, so their sum is k+)"%sp.simplify(B2))
prod=sp.simplify(A1d*C1d*B2)
print("   phase (-i)(-1)(+i) = %s"%sp.simplify((-I)*(-1)*I))
tot=sp.simplify(sp.Integer(2)/(2*pi)**3*prod)
print("   (2/(2pi)^3) x product = %s"%tot)
quoted=1/(2*pi)**3*g**4/(8*pi**4)/K**3*P*(K-P)/Dt
print("   your quoted           = %s"%sp.simplify(quoted))
print("   ratio quoted/derived  = %s   <-- magnitude exact, sign opposite"%sp.simplify(quoted/tot))
print("   (the sign lives in your bar-A / bar-B conventions inside the two brackets; check it there.)")

print(); print("="*80); print("2.  THE C DELTA ELIMINATES w"); print("="*80)
w,z,x=sp.symbols('w z x'); xi=sp.Symbol('xi',positive=True)
sol=sp.solve(sp.Eq((w-z)+xi*(z-x),0),w)[0]
print("   delta^(2)[ (w-z) + (p+/k+)(z-x) ]  =>  w = %s = z - xi r ,  r = z-x"%sp.simplify(sol))
print("   Jacobian in w is 1 (coefficient of w is +1).")
wp=sp.Symbol('wp')
print("   w' - w = %s   =>  e^{-ik(w'-w)} = e^{-ik(w'-z)} e^{-i xi k.r}"%sp.simplify(wp-sol))
print("   and with z = x + r :  = e^{-ik.w'} e^{ik.x} e^{i(1-xi) k.r}   =>   q = (1-xi) k = xibar k")
print("   *** your note writes the phase without the xi (eq 1.3); the momentum is q = xibar k. ***")

print(); print("="*80); print("3.  THE TENSOR CONTRACTION"); print("="*80)
print("   B_2 bracket with x_def->y , z_def->x , w_def->z , S->k+ , k+_def->k+-p+ , p+_def->p+ :")
print("     W^{jk} = delta^{jk}[(y-x)^2-(y-z)^2]/2(x-z)^2")
print("            + (k+/(k+-p+))[ (y-x)^j(x-z)^k/(x-z)^2 - (y-x)^j(y-z)^k/2(y-z)^2 ]")
print("            + (k+/p+)    [ (y-z)^k(x-z)^j/(x-z)^2 + (y-x)^j(y-z)^k/2(y-x)^2 ]")
print("   C bracket is called as C_{1 j k i} -- the SAME order as the definition, no relabelling:")
print("     delta_im delta_jk - (k+/p+) delta_jm delta_ik - (k+/(k+-p+)) delta_km delta_ij")
print("   contracting over j,k :")
print("     R^{im} = delta^{im} tr W - (k+/p+) W^{mi} - (k+/(k+-p+)) W^{im}")
print()
s_=lambda v:v@v; dd=np.eye(2)
def Wjk(y,xv,zv,Pv,Kv):
    A=s_(y-xv); Bv=s_(y-zv); C=s_(xv-zv)
    return ( dd*(A-Bv)/(2*C)
            +(Kv/(Kv-Pv))*( np.outer(y-xv,xv-zv)/C - np.outer(y-xv,y-zv)/(2*Bv) )
            +(Kv/Pv)     *( np.outer(xv-zv,y-zv)/C + np.outer(y-xv,y-zv)/(2*A) ) )
def yours(y,xv,zv,Pv,Kv,dperp=2.0):
    A=s_(y-xv); Bv=s_(y-zv); C=s_(xv-zv); Kv2=Kv*Kv
    yx,yz,xz=y-xv,y-zv,xv-zv
    R = dd*dperp/(2*C)*(A-Bv)
    R+= -Kv2/(Kv-Pv)**2*( np.outer(yx,xz)/C - np.outer(yx,yz)/(2*Bv) )
    R+= dd*(Kv/Pv)*( (yz@xz)/C + (yx@yz)/(2*A) - (A-Bv)/(2*C) )
    R+= -Kv2/Pv**2*( np.outer(yz,xz)/C + np.outer(yz,yx)/(2*A) )
    R+= dd*(Kv/(Kv-Pv))*( (yx@xz)/C - (yx@yz)/(2*Bv) - (A-Bv)/(2*C) )
    R+= -Kv2/(Pv*(Kv-Pv))*( np.outer(xz,yx)/C - np.outer(yz,yx)/(2*Bv)
                           +np.outer(xz,yz)/C + np.outer(yx,yz)/(2*A) )
    return R
rng=np.random.default_rng(5)
print("   %-28s %14s"%("trial","||yours - exact||"))
for t in range(4):
    y,xv,zv=[rng.normal(size=2)*1.2 for _ in range(3)]
    Kv=1.0; Pv=rng.uniform(0.15,0.85)
    W=Wjk(y,xv,zv,Pv,Kv)
    ex=dd*np.trace(W) - (Kv/Pv)*W.T - (Kv/(Kv-Pv))*W
    print("   p+=%.3f k+=%.1f              %14.2e"%(Pv,Kv,np.abs(yours(y,xv,zv,Pv,Kv)-ex).max()))
