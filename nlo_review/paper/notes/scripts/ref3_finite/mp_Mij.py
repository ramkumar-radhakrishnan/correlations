# High-precision check of the I.III kernel K^III_z at z -> w (the double-precision scan loses digits there).
import mpmath as mp
mp.mp.dps = 60
def v(a): return [mp.mpf(a[0]), mp.mpf(a[1])]
def dot(a, b): return a[0]*b[0] + a[1]*b[1]
def sub(a, b): return [a[0]-b[0], a[1]-b[1]]
def outer(a, b): return [[a[i]*b[j] for j in range(2)] for i in range(2)]
def madd(*ms): return [[sum(m[i][j] for m in ms) for j in range(2)] for i in range(2)]
def msc(m, s): return [[m[i][j]*s for j in range(2)] for i in range(2)]
def Kv(a): n = dot(a, a); return [a[0]/n, a[1]/n]
def Mij(x, z, w):
    a = sub(x, w); b = sub(x, z); c = sub(z, w); A = dot(a, a); B = dot(b, b); C = dot(c, c)
    L = mp.log(B/A); d = [[1, 0], [0, 1]]
    t1 = msc(madd(msc(d, 1/(2*C)), msc(outer(a, b), -1/(A*B)),
                  msc(madd(msc(outer(a, c), B/A), msc(outer(c, b), A/B), msc(outer(c, a), -1), msc(outer(b, c), -1)), 1/(C*(A-B)))), L)
    A0 = mp.acos(dot(a, c)/mp.sqrt(A*C))/mp.sqrt(A*B - dot(a, b)**2)
    t2 = msc(madd(outer(c, a), outer(b, c), msc(madd(msc(outer(c, b), 1/B), msc(outer(a, c), 1/A)), -dot(a, b))), A0/C)
    t3 = msc(outer(c, b), -mp.log(C/A)/(2*C)/B)
    return madd(t1, t2, t3)
def Kz(x, y, z, w, xp, wp):
    P = Kv(sub(xp, wp)); M = Mij(x, z, w); k = Kv(sub(y, z))
    m = sum(P[i]*M[i][j]*k[j] for i in range(2) for j in range(2))
    ell = -mp.log(dot(sub(z, w), sub(z, w))/dot(sub(x, w), sub(x, w)))/2
    return m + ell*dot(Kv(sub(y, z)), Kv(sub(w, z)))*dot(P, Kv(sub(x, w)))
import random
random.seed(3)
x = v([0.3, -0.2]); y = v([-0.7, 0.5]); w = v([0.9, 0.8]); wp = v([-1.1, -0.4])
for label, args in (("Kz|y->x, z->w", lambda z: (x, x, z, w, x, wp)), ("Kz|x'->y, z->w", lambda z: (x, y, z, w, y, wp))):
    out = []
    for lam in (1e-3, 1e-5, 1e-7, 1e-9):
        mx = 0; mean = 0; n = 64
        for t in range(n):
            th = 2*mp.pi*t/n + 0.1
            z = [w[0] + lam*mp.cos(th), w[1] + lam*mp.sin(th)]
            val = Kz(*args(z)); mx = max(mx, abs(val)); mean += val/n
        out.append((lam, float(mx*lam), float(mx*lam/abs(mp.log(lam))), float(abs(mean)*lam**2)))
    print(label); [print('  lam=%.0e  lam*max=%.4g  lam*max/|ln lam|=%.4g  lam^2*|mean|=%.3g' % o) for o in out]
# triple z,w,w' for Kz|y->x
print("Kz|y->x, {z,w,w'} -> scale")
for lam in (1e-3, 1e-5, 1e-7):
    mx = 0
    for t in range(32):
        th = 2*mp.pi*t/32 + 0.3
        z = [w[0] + lam*mp.cos(th), w[1] + lam*mp.sin(th)]; wpp = [w[0] + lam*0.7*mp.cos(2*th+1), w[1] + lam*0.7*mp.sin(2*th+1)]
        mx = max(mx, abs(Kz(x, x, z, w, x, wpp)))
    print('  lam=%.0e  lam^2*max=%.4g  lam^2*max/|ln lam|=%.4g' % (lam, float(mx*lam**2), float(mx*lam**2/abs(mp.log(lam)))))
