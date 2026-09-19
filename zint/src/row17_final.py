"""Row 17 (N_2): factorised vertex, the six p+ masters in closed form,
   their Lambda->0 / V->infty asymptotics, and the transverse log."""
import numpy as np, sympy as sp
from scipy.integrate import quad

# ---------------------------------------------------------------- 1. algebra
print("="*78); print("1.  THE VERTEX FACTORISES INTO TWO VECTORS"); print("="*78)
s=lambda v:v@v
def Wji(x,z,w,P,K=1.0):
    S=P+K; xz,xw,zw=x-z,x-w,z-w
    return (np.eye(2)*(s(xz)-s(xw))/(2*s(zw))
            + (S/K)*(np.outer(xz,zw)/s(zw) - np.outer(xz,xw)/(2*s(xw)))
            + (S/P)*(np.outer(zw,xw)/s(zw) + np.outer(xz,xw)/(2*s(xz))))
def taus(x,z,w):
    A,B,C=s(x-z),s(x-w),s(z-w)
    t1=np.eye(2)*(A-B)/(2*C)
    W=(z-w)/C-(x-w)/(2*B)          # tau_2^{ji} = (x-z)^j W^i
    V=(z-w)/C+(x-z)/(2*A)          # tau_3^{ji} = V^j (x-w)^i
    return t1,np.outer(x-z,W),np.outer(V,W*0+ (x-w)),W,V
rng=np.random.default_rng(7)
x,z,w=rng.normal(size=2),rng.normal(size=2),rng.normal(size=2); P,K=0.37,1.0; S=P+K
t1,t2,t3,W,V=taus(x,z,w)
print("   W^{ji} = d^{ji}(A-B)/2C + (S/k+)(x-z)^j W^i + (S/p+) V^j (x-w)^i ,")
print("     W^i = (z-w)^i/(z-w)^2 - (x-w)^i/2(x-w)^2 ,   V^j = (z-w)^j/(z-w)^2 + (x-z)^j/2(x-z)^2")
print("   max|reconstructed - coded|  =  %.3e"%np.abs(t1+(S/K)*t2+(S/P)*t3-Wji(x,z,w,P)).max())

# the seven independent contractions
xp,wp=rng.normal(size=2),rng.normal(size=2)
t1p,t2p,t3p,Wp,Vp=taus(xp,z,wp)
con=lambda a,b: np.einsum('ji,ji->',a,b)
A,B,C=s(x-z),s(x-w),s(z-w); Ap,Bp,Cp=s(xp-z),s(xp-wp),s(z-wp)
pred={
 "t1.t1'":(A-B)*(Ap-Bp)/(2*C*Cp),
 "t1.t2'":(A-B)/(2*C)*((xp-z)@Wp),
 "t1.t3'":(A-B)/(2*C)*(Vp@(xp-wp)),
 "t2.t2'":((x-z)@(xp-z))*(W@Wp),
 "t2.t3'":((x-z)@Vp)*(W@(xp-wp)),
 "t3.t3'":(V@Vp)*((x-w)@(xp-wp)),
 "t3.t2'":(V@(xp-z))*((x-w)@Wp)}
got={"t1.t1'":con(t1,t1p),"t1.t2'":con(t1,t2p),"t1.t3'":con(t1,t3p),
     "t2.t2'":con(t2,t2p),"t2.t3'":con(t2,t3p),"t3.t3'":con(t3,t3p),"t3.t2'":con(t3,t2p)}
print("\n   %-8s %16s %16s %9s"%("contraction","closed form","direct","diff"))
for kk in pred: print("   %-8s %16.9f %16.9f %9.1e"%(kk,pred[kk],got[kk],abs(pred[kk]-got[kk])))

# ------------------------------------------------- 2. masters, symbolic check
print(); print("="*78); print("2.  THE SIX p+ MASTERS: SYMBOLIC ANTIDERIVATIVE CHECK"); print("="*78)
p,k,a,b=sp.symbols('p k alpha beta',positive=True)   # beta = alpha'
L=sp.log
F={
 "h=1          (T2T3', T3T2')": L((p+a)/(p+b))/(b-a),
 "h=p/k+       (T2T2')       ": (-a*L(p+a)+b*L(p+b))/(k*(b-a)),
 "h=k+/p       (T3T3')       ": k*(L(p)/(a*b)+L(p+a)/(a*(a-b))+L(p+b)/(b*(b-a))),
 "h=k+/(p+k+)  (T1T3', T3T1')": k*(L(p+k)/((a-k)*(b-k))+L(p+a)/((k-a)*(b-a))+L(p+b)/((k-b)*(a-b))),
 "h=p/(p+k+)   (T1T2', T2T1')": (-k*L(p+k)/((a-k)*(b-k))-a*L(p+a)/((k-a)*(b-a))-b*L(p+b)/((k-b)*(a-b))),
}
A2=-k/((a-k)*(b-k)); A1=(a*b-k**2)/((a-k)**2*(b-k)**2)
Bc=-a/((k-a)**2*(b-a)); Bp_=-b/((k-b)**2*(a-b))
F["h=p k+/(p+k+)^2 (T1T1')  "]=k*(-A2/(p+k)+A1*L(p+k)+Bc*L(p+a)+Bp_*L(p+b))
H={"h=1          (T2T3', T3T2')":sp.Integer(1),
   "h=p/k+       (T2T2')       ":p/k,
   "h=k+/p       (T3T3')       ":k/p,
   "h=k+/(p+k+)  (T1T3', T3T1')":k/(p+k),
   "h=p/(p+k+)   (T1T2', T2T1')":p/(p+k),
   "h=p k+/(p+k+)^2 (T1T1')  ":p*k/(p+k)**2}
for nm in F:
    d=sp.simplify(sp.diff(F[nm],p)-H[nm]/((p+a)*(p+b)))
    print("   %-30s  dF/dp - h/[(p+a)(p+a')] = %s"%(nm,d))

# ------------------------------------------------- 3. asymptotics
print(); print("="*78); print("3.  LIMITS  p+ in [Lambda, V-k+]   (Lambda->0, V->infty)"); print("="*78)
for nm in F:
    lo=sp.simplify(sp.series(F[nm].subs(p,sp.Symbol('Lam',positive=True)),sp.Symbol('Lam',positive=True),0,1).removeO())
    hi=sp.limit(sp.simplify(F[nm]-sp.Symbol('c')),p,sp.oo) if False else None
    ser=sp.simplify(F[nm].rewrite(sp.log))
    # large-p behaviour
    big=sp.simplify(sp.limit(F[nm]/sp.log(p),p,sp.oo))
    print("   %-30s  F(Lambda->0) = %-38s   F ~ (%s) log V"%(nm,sp.simplify(lo),big))

# ------------------------------------------------- 4. numbers
print(); print("="*78); print("4.  NUMERICAL CHECK OF THE ASSEMBLED p+ INTEGRAL"); print("="*78)
kk_,Av,Bv,Apv,Bpv,Lam,Vv=1.0,0.8,1.7,1.3,0.6,1e-8,1e7
al,be=kk_*Bv/Av,kk_*Bpv/Apv
sub={k:kk_,a:al,b:be}
print("   k+=1  A=%.1f B=%.1f A'=%.1f B'=%.1f  alpha=%.6f alpha'=%.6f  Lam=%g V=%g"%(Av,Bv,Apv,Bpv,al,be,Lam,Vv))
print("   %-30s %16s %16s %9s"%("master","quadrature","closed form","diff"))
for nm in F:
    hf=sp.lambdify(p,H[nm].subs(sub),'numpy'); Ff=sp.lambdify(p,F[nm].subs(sub),'numpy')
    q=quad(lambda u:hf(np.exp(u))*np.exp(u)/((np.exp(u)+al)*(np.exp(u)+be)),
           np.log(Lam),np.log(Vv-kk_),limit=900)[0]/(Av*Apv)
    cf=(Ff(Vv-kk_)-Ff(Lam))/(Av*Apv)
    print("   %-30s %16.9f %16.9f %9.1e"%(nm,q,cf,abs(q-cf)))

# ------------------------------------------------- 5. the rapidity-log residue
print(); print("="*78); print("5.  COEFFICIENT OF log(1/Lambda):  ONLY T3T3'"); print("="*78)
print("   M_33 -> (1/(A A')) k+ log(1/Lam)/(alpha alpha') = log(1/Lam)/(k+ B B')")
print("   times  tau_3.tau_3' = (V.V')[(x-w).(x'-w')]   =>")
print("     log(1/Lam) * (1/k+) * [(x-w).(x'-w')/((x-w)^2 (x'-w')^2)] * (V.V')")
num=(1/(Av*Apv))*kk_/(al*be); print("   check residue: 1/(k+ B B') = %.9f   vs  %.9f"%(1/(kk_*Bv*Bpv),num))

print(); print("   large-|z| behaviour of V.V' :   V -> z/(2|z|^2)  =>  V.V' -> 1/(4|z|^2)")
for R in (1e2,1e3,1e4,1e5):
    acc=0.0
    for t in np.linspace(0,2*np.pi,4096,endpoint=False):
        zz=R*np.array([np.cos(t),np.sin(t)])
        _,_,_,_,Vz=taus(x,zz,w); _,_,_,_,Vzp=taus(xp,zz,wp)
        acc+=Vz@Vzp
    print("     |z|=%-8.0e  <V.V'> x 4|z|^2 = %.9f"%(R,acc/4096*4*R*R))
print("   -> the log(1/Lambda) coefficient still has  int d^2z/(4|z|^2) : a log R^2.")
