"""Row 9: the three N1 p+ masters, closed form + numerics (log substitution)."""
import numpy as np
from scipy.integrate import quad
k,A,B,Lam,V = 1.0, 0.8, 1.7, 1e-8, 1e7
Pv = V-k
Q = lambda f,l,h: quad(lambda u: f(np.exp(u))*np.exp(u), np.log(l), np.log(h), limit=900)[0]
print("="*76); print("N1: three p+ masters,  D = p+ A + k+ B"); print("="*76)
print("   k+=1, A=(x-z)^2=0.8, B=(x-w)^2=1.7, Lambda=1e-8, V=1e7\n")
rows = [
 ("J1 = int dp/[p D]        (the 1/p+ term)",
  lambda t: 1/(t*(t*A+k*B)),
  1/(k*B)*np.log(Pv*k*B/(Lam*(Pv*A+k*B))),
  "(1/(k+B)) log[ (V-k+) k+B / (Lambda((V-k+)A+k+B)) ]  ->  (1/(k+B)) log(k+B/(Lambda A))"),
 ("J2 = int dp/D            (the 1/k+ term)",
  lambda t: 1/(t*A+k*B),
  1/A*np.log((Pv*A+k*B)/(Lam*A+k*B)),
  "(1/A) log[ ((V-k+)A+k+B)/(Lambda A+k+B) ]            ->  (1/A) log( V A/(k+B) )   [log V]"),
 ("J3 = int dp/[(p+k) D]    (the delta^ij term)",
  lambda t: 1/((t+k)*(t*A+k*B)),
  1/(k*(B-A))*np.log(V*B/(Pv*A+k*B)),
  "(1/(k+(B-A))) log[ V B/((V-k+)A+k+B) ]               ->  (1/(k+(B-A))) log(B/A)"),
]
for nm,f,cf,form in rows:
    num = Q(f,Lam,Pv)
    print("  %s" % nm)
    print("     %s" % form)
    print("     quad %16.9f   closed %16.9f   diff %.2e\n" % (num,cf,abs(num-cf)))
print("="*76); print("THE delta^ij TERM COLLAPSES EXACTLY"); print("="*76)
C = 0.55   # (z-w)^2
lhs = (A-B)/(2*C) * 1/(k*(B-A))*np.log(B/A)
rhs = -1/(2*k*C)*np.log(B/A)
print("   coefficient (A-B)/(2C) x J3  =  %.12f" % lhs)
print("   -(1/(2 k+ C)) log(B/A)       =  %.12f" % rhs)
print("   -> exactly  - delta^ij log[(x-w)^2/(x-z)^2] / [2 k+ (z-w)^2] :  no V, no Lambda")
