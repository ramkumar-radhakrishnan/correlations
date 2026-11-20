import numpy as np
from probe_g2 import total
rng = np.random.default_rng(9)
mx = 0
for _ in range(2000):
    x, xp, z, w, wp = rng.normal(size=(5, 2))
    a = total((x, xp, z, w, wp)); b = total((xp, x, z, wp, w))
    mx = max(mx, abs(a - b) / (abs(a) + abs(b)))
print('max relative asymmetry of sum (7.3)-(7.8) under (x,w)<->(x\',w\'):', mx)
