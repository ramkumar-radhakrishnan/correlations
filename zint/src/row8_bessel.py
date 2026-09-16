import numpy as np
from scipy.special import j0,j1
print("EXACT:  <n^m n^k e^{i x n.khat}> = (J1(x)/x) d^{mk} + (J0(x) - 2 J1(x)/x) khat^m khat^k")
print("%-9s %14s %14s %14s %14s" % ("x","numeric a","J1(x)/x","numeric b","J0-2J1/x"))
for x in (10.,100.,1000.,1e4,1e5):
    nth=int(max(20000,60*x)); th=np.linspace(0,2*np.pi,nth,endpoint=False)
    c,s=np.cos(th),np.sin(th); ph=np.exp(1j*x*c)
    A11=(c*c*ph).mean(); A22=(s*s*ph).mean()
    a=A22.real; b=(A11-A22).real
    print("%-9.0e %14.6e %14.6e %14.6e %14.6e" % (x,a,j1(x)/x,b,j0(x)-2*j1(x)/x))
print()
print("J0(x) ~ sqrt(2/(pi x)) cos(x-pi/4)  => envelope x^{-1/2};   J1(x)/x ~ x^{-3/2}.")
print("y, y' carry the phase  x = |k| R S/k+  ->  |A(R)| R^2 ~ R^{-1/2} : CONVERGENT.")
print("z carries NO phase     x = 0          ->  J1(x)/x -> 1/2, <n^m n^m'> = d^{mm'}/2,")
print("                                          |A(R)| R^2 -> const : LOG DIVERGENT.")
print("x, x' carry no phase either, but their single kernel's first two angular")
print("moments vanish identically, so they converge as R^{-1}.")
