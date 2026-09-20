"""Row 22: the two-C row -- prefactor, w,w' integration, tensor contraction, colour."""
import numpy as np, sympy as sp, itertools
g,pi=sp.symbols('g pi',positive=True); P,K=sp.symbols('pplus kplus',positive=True); S=P+K; I=sp.I
print("="*80); print("1.  PREFACTOR"); print("="*80)
A1 = I*g/(sp.sqrt(2)*pi*sp.sqrt(S))                      # A^(1)(p+ + k+) coefficient
C1 = -g/(2*pi*sp.sqrt(2*S))*sp.sqrt(P*(S-P))/S           # C_1 coefficient at k+ -> S
mod=lambda e: sp.simplify(e*sp.conjugate(e).subs({sp.conjugate(g):g,sp.conjugate(pi):pi,
          sp.conjugate(P):P,sp.conjugate(K):K,sp.conjugate(S):S}))
got=sp.simplify(sp.Integer(4)/(2*pi)**3*mod(A1)*mod(C1))
want=1/(2*pi)**3*g**4/(4*pi**4)*P*K/S**4
print("   |A^(1)|^2 = %s"%sp.simplify(mod(A1)))
print("   |C_1|^2   = %s   (uses sqrt(p+(S-p+))/S = sqrt(p+k+)/S)"%sp.simplify(mod(C1)))
print("   (4/(2pi)^3)|A|^2|C|^2 = %s"%sp.simplify(got))
print("   your      (1/(2pi)^3)(g^4/4pi^4) p+k+/S^4  ->  ratio = %s"%sp.simplify(got/want))

print(); print("="*80); print("2.  THE w, w' INTEGRATION"); print("="*80)
print("   delta^(2)[ y - w + (p+/S)(w-z) ] : coefficient of w is -(1 - p+/S) = -k+/S")
print("     => solves to  w = (S y - p+ z)/k+ ,  Jacobian (S/k+)^2 PER delta, (S/k+)^4 total")
print("   but the C kernel rescales:  w - z = (S/k+)(y - z)  =>")
print("     (w-z)^m/(w-z)^2 = (k+/S) (y-z)^m/(y-z)^2 ,  two of them: (k+/S)^2")
print("   NET factor = (S/k+)^4 x (k+/S)^2 = (S/k+)^2 = S^2/k+^2    <-- what you used")
w,z,y,p,k=sp.symbols('w z y pplus kplus'); Ssy=p+k
wsol=sp.solve(sp.Eq(y-w+p/(p+k)*(w-z),0),w)[0]
print("   check:  w solved = %s ;  (w-z) = %s"%(sp.simplify(wsol),sp.simplify(sp.factor(wsol-z))))
wp=wsol.subs(y,sp.Symbol('yp'))
print("   check:  w' - w = %s   =>  phase e^{-ik(w'-w)} = e^{-ik(y'-y)S/k+}"%sp.simplify(sp.factor(wp-wsol)))

print(); print("="*80); print("3.  THE TENSOR CONTRACTION"); print("="*80)
d=2; dd=np.eye(d)
def Mten(Pv,Kv):
    Sv=Pv+Kv; M=np.zeros((d,d,d,d))          # M[m,kk,i,j]
    for m,kk,i,j in itertools.product(range(d),repeat=4):
        M[m,kk,i,j]=(dd[kk,m]*dd[i,j] - (Sv/Pv)*dd[j,m]*dd[i,kk] - (Sv/Kv)*dd[i,m]*dd[kk,j])
    return M
for Pv,Kv in ((0.37,1.0),(2.4,0.6),(1.0,1.0)):
    Sv=Pv+Kv; M=Mten(Pv,Kv); Con=np.einsum('akij,bcij->akbc',M,M)   # sum over i,j
    c1=d-2*Sv/Pv-2*Sv/Kv; c2=Sv**2*(1/Pv**2+1/Kv**2); c3=2*Sv**2/(Pv*Kv)
    pred=np.zeros((d,d,d,d))
    for m,kk,mp,kp in itertools.product(range(d),repeat=4):
        pred[m,kk,mp,kp]=c1*dd[kk,m]*dd[kp,mp]+c2*dd[m,mp]*dd[kk,kp]+c3*dd[m,kp]*dd[kk,mp]
    print("   p+=%.2f k+=%.2f : ||contraction - [c1 d_km d_k'm' + c2 d_mm' d_kk' + c3 d_mk' d_km']|| = %.2e"
          %(Pv,Kv,np.abs(Con-pred).max()))
print("     c1 = d_perp - 2S/p+ - 2S/k+ ,  c2 = S^2(1/p+^2 + 1/k+^2) ,  c3 = 2S^2/(p+ k+)")
print()
print("   now multiply by  prefactor x net-w-factor x k+  =  (p+k+/S^4)(S^2/k+^2)(k+) = p+/S^2 :")
dperp=sp.Symbol('d_perp')
b1=sp.simplify(P/S**2*(dperp-2*S/P-2*S/K)); b2=sp.simplify(P/S**2*S**2*(1/P**2+1/K**2))
b3=sp.simplify(P/S**2*2*S**2/(P*K))
print("     d_km d_k'm'  : %s"%sp.simplify(b1))
print("        yours     : d_perp p+/S^2 - 2/k+          -> difference %s"%sp.simplify(b1-(dperp*P/S**2-2/K)))
print("     d_mm' d_kk'  : %s"%sp.simplify(b2))
print("        yours     : p+/k+^2 + 1/p+                -> difference %s"%sp.simplify(b2-(P/K**2+1/P)))
print("     d_mk' d_km'  : %s"%sp.simplify(b3))
print("        yours     : 2/k+                          -> difference %s"%sp.simplify(b3-2/K))
