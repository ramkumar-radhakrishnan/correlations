"""Row 9: the p+ integrals for N1 and N2."""
import sympy as sp
p, K, V, Lam = sp.symbols('p kplus V Lambda', positive=True)
A, B, Ap, Bp = sp.symbols('A B Aprime Bprime', positive=True)   # (x-z)^2,(x-w)^2,(x'-z)^2,(x'-w')^2
S = p + K
Pv = V - K

def do(expr, name):
    prim = sp.integrate(expr, p)
    val  = sp.simplify(prim.subs(p, Pv) - prim.subs(p, Lam))
    return prim, val

print("="*78); print("N1 :  three masters,  D = p A + k B"); print("="*78)
J = {}
J['1/p+ term  J1'] = 1/(p*(p*A + K*B))
J['1/k+ term  J2'] = 1/(p*A + K*B)
J['1/S  term  J3'] = 1/((p+K)*(p*A + K*B))
for nm, e in J.items():
    prim = sp.simplify(sp.integrate(e, p))
    lo   = sp.simplify(sp.limit(prim.subs(p, Lam), Lam, 0, '+'))   # divergent part kept symbolic
    print("   %-16s  primitive = %s" % (nm, sp.simplify(prim)))
print()
print("   evaluated on [Lambda, V-k+], small Lambda:")
print("     J1 = (1/(k+ B)) log[ (V-k+) k+ B / (Lambda ((V-k+)A + k+ B)) ]")
print("        --> (1/(k+ B)) log( k+ B /(Lambda A) )            V-STABLE, carries the rapidity log")
print("     J2 = (1/A) log[ ((V-k+)A + k+ B) / (Lambda A + k+ B) ]")
print("        --> (1/A) log( V A /(k+ B) )                      grows like log V")
print("     J3 = (1/(k+(B-A))) log[ V B / ((V-k+)A + k+ B) ]")
print("        --> (1/(k+(B-A))) log( B/A )                      V-STABLE")
# numeric confirmation
import numpy as np
from scipy.integrate import quad
print()
print("   numerical check (k+=1, A=0.8, B=1.7, Lambda=1e-8, V=1e7):")
Kn,An,Bn,Ln,Vn = 1.0,0.8,1.7,1e-8,1e7
q1 = quad(lambda t: 1/(t*(t*An+Kn*Bn)), Ln, Vn-Kn, limit=500)[0]
q2 = quad(lambda t: 1/(t*An+Kn*Bn), Ln, Vn-Kn, limit=500)[0]
q3 = quad(lambda t: 1/((t+Kn)*(t*An+Kn*Bn)), Ln, Vn-Kn, limit=500)[0]
f1 = 1/(Kn*Bn)*np.log((Vn-Kn)*Kn*Bn/(Ln*((Vn-Kn)*An+Kn*Bn)))
f2 = 1/An*np.log(((Vn-Kn)*An+Kn*Bn)/(Ln*An+Kn*Bn))
f3 = 1/(Kn*(Bn-An))*np.log(Vn*Bn/((Vn-Kn)*An+Kn*Bn))
for nm,a,b in [("J1",q1,f1),("J2",q2,f2),("J3",q3,f3)]:
    print("      %s: quad %14.8f   closed %14.8f   diff %.2e" % (nm,a,b,abs(a-b)))

print()
print("="*78); print("N1 :  the delta^ij term collapses"); print("="*78)
print("   coefficient is  (A-B)/(2 C) * J3  with J3 = log(B/A)/(k+(B-A)) :")
print("      (A-B)/(2C) * 1/(k+(B-A)) log(B/A)  =  -(1/(2 k+ C)) log(B/A)")
print("   i.e.  - delta^ij log[(x-w)^2/(x-z)^2] / ( 2 k+ (z-w)^2 )     -- exact, no V, no Lambda")

print()
print("="*78); print("N2 :  the six longitudinal structures"); print("="*78)
P_, Kk = sp.symbols('pplus kplus', positive=True); Ss = P_+Kk
rows = [("T1 T1'",            sp.simplify(P_*Kk/Ss**2)),
        ("T1 T2' + T2 T1'",   sp.simplify(P_*Kk/Ss**2*(Ss/Kk))),
        ("T1 T3' + T3 T1'",   sp.simplify(P_*Kk/Ss**2*(Ss/P_))),
        ("T2 T2'",            sp.simplify(P_*Kk/Ss**2*(Ss/Kk)**2)),
        ("T3 T3'",            sp.simplify(P_*Kk/Ss**2*(Ss/P_)**2)),
        ("T2 T3' + T3 T2'",   sp.simplify(P_*Kk/Ss**2*(Ss**2/(P_*Kk))))]
for nm,e in rows: print("   %-20s ->  %s" % (nm, e))
print("   (T1 = delta^ij piece, T2 = the S/k+ piece, T3 = the S/p+ piece, all p+-independent)")
print("   the ONLY soft pole is  k+/p+  from T3 T3' : the rapidity log lives there.")
