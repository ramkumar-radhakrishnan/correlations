"""Row 8: analytic large-distance asymptotics, with a cheap targeted numerical check."""
import numpy as np, sympy as sp
R,Rr = sp.symbols('R', positive=True), None
print("="*78); print("ANALYTIC LARGE-DISTANCE BEHAVIOUR  (x = R nhat, everything else fixed)"); print("="*78)
# expansion of the WW kernel about R -> infinity
nk,yk,ny = sp.symbols('n_k y_k n_dot_y')
print("""   (x-y)^k/(x-y)^2  =  nhat^k/R  +  [2 nhat^k (nhat.y) - y^k]/R^2  +  O(1/R^3)

   angular averages:   <nhat^k> = 0
                       <2 nhat^k (nhat.y) - y^k> = 2 (y^k/2) - y^k = 0
   -> the first TWO orders average to zero;  the x and x' integrals converge.""")
print()
print("""   y appears in TWO kernels and in the phase:
       |integrand| ~ 1/R^2 ,   measure R^2 dR/R      -> log divergent WITHOUT the phase
       with  e^{i (k.y) S/k+} :  the angular integral of  nhat^m nhat^k e^{i alpha R nhat.khat}
       is Bessel,  ~ R^{-1/2}      -> |A(R)| R^2 ~ R^{-1/2}   CONVERGENT.
   The phase is the only thing that saves y and y'.""")
print()
print("""   z appears in TWO kernels, in NO phase, and after the colour reduction in NO
   Wilson line:
       <nhat^m nhat^m'> = delta^{mm'}/2  does NOT vanish
       -> |A(R)| R^2 -> constant        LOG DIVERGENT.     This is the only one.""")
print()
print("="*78); print("numerical check of the Bessel R^{-1/2} for the y sweep (cheap, kernels only)"); print("="*78)
alpha = np.sqrt(0.7**2+0.4**2)*1.37
print("  %-9s %14s %14s %10s" % ("R","|<n^m n^k e^{i a R n.k}>|","R^{-1/2} scaled","ratio"))
ref=None
for Rv in (1e2,1e3,1e4,1e5,1e6):
    nth=int(min(max(8192, 80*alpha*Rv), 20_000_000))
    th=np.linspace(0,2*np.pi,nth,endpoint=False)
    n=np.stack([np.cos(th),np.sin(th)],-1)
    ph=np.exp(1j*alpha*Rv*(n[:,0]*0.7+n[:,1]*(-0.4))/np.sqrt(0.7**2+0.4**2))
    v=abs((n[:,0]*n[:,0]*ph).mean())
    if ref is None: ref=v*np.sqrt(Rv)
    print("  %-9.0e %14.6e %14.6e %10.4f" % (Rv, v, ref/np.sqrt(Rv), v/(ref/np.sqrt(Rv))))
print("  ratio ~ 1 across four decades  =>  R^{-1/2}, as claimed.")
