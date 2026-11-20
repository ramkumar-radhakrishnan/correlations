# Large-distance (soft gluon z -> infinity) behaviour of the (2rho) contact lines that carry soft remainders,
# paired as in Table tab:pairs: pair A = I.III K_w lines + III.V Psi lines (y'=w'); pair D = E-part of the ell term of
# the I.III K_z lines + III.VI Psi lines (x<->y relabelled, y'=w').  Smooth, non-uniform U (plane waves: U(z) does not
# approach a limit at large z).  Printed: max over 8 directions of |z|^2 |line|, single member and pair sum.
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, lib, lines
field = lib.SmoothField(seed=7); col = lib.Colour(field)
def find(row, name, drop, keep):
    for u in lines.U:
        if u.row == row and u.name == name and u.drop == drop and u.keep == keep: return u
    raise KeyError((row, name, drop, keep))
def KD_E(P):   # E-part of the ell term of K^III_z (Group I member of pair D)
    x, y, z, w, xp, wp = (P[k] for k in ('x', 'y', 'z', 'w', "x'", "w'"))
    return lib.E_fun(lib.KVEC @ (z - w)) * lib.Kf(y, w, z) * (lib.Kv(xp - wp) @ lib.Kv(x - w))
rng = np.random.default_rng(5)
base = {p: rng.uniform(-1, 1, 2) for p in ['x', 'y', "x'", 'w', "w'"]}
dirs = [np.array([np.cos(a), np.sin(a)]) for a in np.linspace(0.13, 2 * np.pi + 0.13, 8, endpoint=False)]
def run(name, mk_I, mk_III):
    print(name)
    for R in (1e1, 1e2, 1e3, 1e4):
        mI = mS = 0
        for d in dirs:
            P = dict(base); P['z'] = R * d
            vI = mk_I(P); vIII = mk_III(P)
            mI = max(mI, R * R * abs(vI)); mS = max(mS, R * R * abs(vI + vIII))
        print('   |z|=%.0e   max z^2|Group I member| = %.4f   max z^2|pair sum| = %.3e' % (R, mI, mS))
# pair A
for (d_, k_) in (("x'", 'x'), ('y', 'x'), ("x'", 'y')):
    uI = find('I.III', 'Kw', d_, k_); uIII = find('III.V', 'Psi', d_, k_)
    run('pair A, contact %s->%s' % (d_, k_), lambda P: uI.value(col, P),
        lambda P: uIII.value(col, dict(P, **{"y'": P["w'"]})))
# pair D
import copy
for (d_, k_), (d3, k3) in (((("x'", 'x'), ("x'", 'y'))), ((('y', 'x'), ('y', 'x'))), ((("x'", 'y'), ("x'", 'x')))):
    uz = find('I.III', 'Kz', d_, k_); u6 = find('III.VI', 'Psi', d3, k3)
    def mkI(P, uz=uz):
        Q = uz.full(P); return uz.coef * KD_E(Q) * uz.colour(col, Q)
    def mkIII(P, u6=u6):
        Q = dict(P); Q['x'], Q['y'] = P['y'], P['x']; Q["y'"] = P["w'"]
        if u6.drop == 'y' and u6.keep == 'x':     # relabelled contact {x,y}: keep the same physical point
            Q['x'] = Q['y'] = P['x']
        return u6.value(col, Q)
    run('pair D, contact %s->%s (III.VI %s->%s, x<->y)' % (d_, k_, d3, k3), mkI, mkIII)
