"""Row 14: A^(1)dag Cbar^dag x B_2.  Prefactor, delta, and the tensor contraction."""
import numpy as np, sympy as sp
g,pi=sp.symbols('g pi',positive=True); P,Kp=sp.symbols('pplus kplus',positive=True); S=P+Kp; I=sp.I
print("="*78); print("1.  PREFACTOR AND THE C DELTA"); print("="*78)
Ad  = -I*g/(sp.sqrt(2)*pi*sp.sqrt(S))                       # A^(1)dag(S, -y')
Cd  = -g/(2*pi*sp.sqrt(2))/sp.sqrt(S)*sp.sqrt(P*Kp)/S        # Cbar^dag(S,.;p+,.)  (real)
B2  =  I*g**2*sp.sqrt(P*Kp)/(4*pi**2)                        # B_2 one-rho, 1/(D S) pulled into M
over= sp.Integer(4)/(2*pi)**3
print("   phases (-i)(-1)(+i) =", sp.simplify((-I)*(-1)*I), "   -> overall MINUS, matching your -g^4/4pi^4")
num = sp.simplify(over*Ad*Cd*B2/( (-I)*(-1)*I ))
print("   numeric x momentum  =", sp.simplify(over/(sp.sqrt(2)*pi)/(2*pi*sp.sqrt(2))/(4*pi**2)),
      " = (1/(2pi)^3) x", sp.simplify(over/(sp.sqrt(2)*pi)/(2*pi*sp.sqrt(2))/(4*pi**2)*(2*pi)**3))
print("   momentum powers     : (1/sqrt S)(1/sqrt S)(sqrt(p+k+)/S)(sqrt(p+k+)) = p+k+/S^2")
print("   delta^(2)[(y'-w') + (p+/S)(w'-z)]  ->  w' = (S y' - p+ z)/k+ ,  w'-z = (S/k+)(y'-z)")
print("     C kernel gives k+/S ; Jacobian (S/k+)^2  ->  net  (p+k+/S^2)(S^2/k+^2)(k+/S) = p+/S")
print("   => overall  -(1/(2pi)^3)(g^4/4pi^4) x (p+/S)   MATCHES your quoted prefactor")
print()
print("   phase: w'-w = (y'-w) + (p+/k+)(y'-z)   =>")
print("      e^{-ik(w'-w)} = e^{-ik(y'-w)} e^{-ik (p+/k+)(y'-z)}      MATCHES exactly")

print(); print("="*78); print("2.  THE TENSOR CONTRACTION"); print("="*78)
print("   C bracket (slots j,i,k'):  [ d_{k'm} d_{ji} - (S/p+) d_{jm} d_{k'i} - (S/k+) d_{im} d_{k'j} ]")
print("   B_2 gives   M^{ij} = d^{ij}(A-B)/(2 S C) + K^{ij}/k+ + P^{ij}/p+ ,")
print("      K^{ij} = (x-z)^j (z-w)^i/C - (x-z)^j (x-w)^i/(2B)")
print("      P^{ij} = (x-w)^i (z-w)^j/C + (x-w)^i (x-z)^j/(2A)")
print("   contracting i,j:   T^{k'm} = d_{k'm} tr M - (S/p+) M^{k'm} - (S/k+) M^{m k'}")
print("   times the overall p+/S :")
print("      (p+/S) d_{k'm} tr M   -   M^{k'm}   -   (p+/k+) M^{m k'}")
print()
rng=np.random.default_rng(9); s=lambda v:v@v
Kn,Pn = 1.0, 0.37; Sn=Kn+Pn; dperp=2.0
def build(x,z,w):
    A=s(x-z); B=s(x-w); C=s(z-w)
    Km=np.outer(z-w,x-z)/C - np.outer(x-w,x-z)/(2*B)      # K^{ij}, index [i,j]
    Pm=np.outer(x-w,z-w)/C + np.outer(x-w,x-z)/(2*A)      # P^{ij}
    M = np.eye(2)*(A-B)/(2*Sn*C) + Km/Kn + Pm/Pn
    return A,B,C,Km,Pm,M
def Tref(x,z,w):
    A,B,C,Km,Pm,M=build(x,z,w)
    return (Pn/Sn)*np.eye(2)*np.trace(M) - M - (Pn/Kn)*M.T
def Tuser(x,z,w,fix=False):
    A,B,C,Km,Pm,M=build(x,z,w)
    d=np.eye(2)
    t1 = d*dperp*Pn*(A-B)/(2*C*Sn**2)
    t2 = Pn*d/(Kn*Sn)*( ((x-z)@(z-w))/C - ((x-z)@(x-w))/(2*B) - (A-B)/(2*C) )
    t3 = d/Sn*( ((x-w)@(z-w))/C + ((x-z)@(x-w))/(2*A) - (A-B)/(2*C) )
    second = np.outer(x-w,x-z)/(2*A) if fix else np.outer(z-w,x-z)/(2*A)
    t4 = -1/Pn*( np.outer(x-w,z-w)/C + second )            # index [k',m]
    t5 = -Pn/Kn**2*( np.outer(x-z,z-w)/C - np.outer(x-z,x-w)/(2*B) )
    t6 = -1/Kn*( np.outer(z-w,x-z)/C - np.outer(x-w,x-z)/(2*B)
                 + np.outer(z-w,x-w)/C + np.outer(x-z,x-w)/(2*A) )
    return t1+t2+t3+t4+t5+t6
print("   %-6s %16s %16s %16s"%("trial","||T_ref||","||as written||","||with the fix||"))
for t in range(4):
    x,z,w=[rng.normal(size=2)*1.1 for _ in range(3)]
    R=Tref(x,z,w); U=Tuser(x,z,w,False); F=Tuser(x,z,w,True)
    print("   %-6d %16.8f  diff %10.3e   diff %10.3e"
          %(t,np.linalg.norm(R),np.linalg.norm(R-U),np.linalg.norm(R-F)))
print()
print("   -> every term matches EXCEPT the second half of your 1/p+ bracket.")
print("      you wrote   (x-z)^m (z-w)^{k'} / (2 (x-z)^2)")
print("      it should be (x-z)^m (x-w)^{k'} / (2 (x-z)^2)")
print("      because P^{ij} = (x-w)^i[ (z-w)^j/C + (x-z)^j/(2A) ] -- BOTH halves carry (x-w)^i.")
