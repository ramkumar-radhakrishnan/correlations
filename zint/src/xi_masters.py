"""The six elementary xi-integrals  J_i = \int_lam^{1-lam} dxi w_i(xi)/Delta(xi),
   Delta = (1-xi) A + xi B ,  A=(s-r)^2=(Y')^2 , B=s^2=Y^2 , C=B-A=2 s.r-r^2 ,
   l_k = ln(k^+/Lambda) = ln(1/lam) ,  L = ln(B/A).
   w = (xi*xib, xi/xib, xib, xib/xi, xi, 1)."""
import numpy as np
from scipy.integrate import quad

def J(A,B,lk):
    C=B-A; L=np.log(B/A)
    return {1: 1/(2*C) + A/C**2 - A*B/C**3*L,
            2: lk/B - A/(B*C)*L,
            3: B/C**2*L - 1/C,
            4: lk/A - B/(A*C)*L,
            5: 1/C - A/C**2*L,
            6: L/C}

if __name__ == "__main__":
    for (A,B) in [(1.0,4.0),(2.3,0.7),(0.9,1.1)]:
        lam=1e-9; lk=np.log(1/lam)
        w=[lambda x:x*(1-x), lambda x:x/(1-x), lambda x:1-x,
           lambda x:(1-x)/x, lambda x:x, lambda x:1.0]
        num=[quad(lambda x,f=f:f(x)/((1-x)*A+x*B),lam,1-lam,limit=400,
                  epsabs=1e-14,epsrel=1e-13)[0] for f in w]
        an=J(A,B,lk)
        print(f"A={A} B={B}")
        for i in range(6):
            print(f"   J{i+1}: analytic {an[i+1]:+.12f}  quad {num[i]:+.12f}  "
                  f"rel {abs(an[i+1]-num[i])/abs(num[i]):.1e}")
    # the paper's C.5 / C.6
    A,B=1.0,4.0; lam=1e-9; lk=np.log(1/lam); L=np.log(B/A)
    q5=quad(lambda x:1/(x*((1-x)*A+x*B)),lam,1-lam,limit=400)[0]
    q6=quad(lambda x:1/((1-x)*((1-x)*A+x*B)),lam,1-lam,limit=400)[0]
    print(f"\n(C.5): quad {q5:.10f}   paper -(1/A)(ln(Lam/k)+ln(B/A)) = {(lk-L)/A:.10f}  OK")
    print(f"(C.6): quad {q6:.10f}   paper -(1/B)(ln(Lam/k)+ln(B/A)) = {(lk-L)/B:.10f}"
          f"   correct value (1/B)(l_k+L) = {(lk+L)/B:.10f}")
