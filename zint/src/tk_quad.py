import numpy as np, sys
sys.path.insert(0,'/home/user/correlations/zint/src')
from twokernel import Kvec, closed_form_vec, integrate_partition
rng = np.random.default_rng(7)
print("2. THE FULL z-INTEGRAL  vs  the closed form  (random configurations)")
for trial in range(4):
    P = rng.normal(size=(4,2))*1.2
    wp,w,x,y = P
    num = integrate_partition(lambda z: Kvec(wp,w,z)*Kvec(x,y,z), [wp,w,x,y],
                              ell=1.0, ntheta=1501, Rmax=1e4)
    an = closed_form_vec(wp,w,x,y)
    print(f"   trial {trial}: quad {num:+.9f}   closed form {an:+.9f}   rel {abs(num-an)/abs(an):.2e}",
          flush=True)
