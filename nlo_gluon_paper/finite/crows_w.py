"""Group III rows with the splitting after the target (III part 1, V, VI; III part 2), in the w' variables.
At fixed w' the delta function fixes y' = (k w' + p z)/(p + k), with Jacobian (k/P)^2.  The p+ integrand f(p) behaves
as A/p for p -> 0 (y' -> w') and as P/p for p -> infinity (y' -> z).  So, exactly up to O(Lambda/k+, k+/vee),
   int_Lambda^{vee-k} dp/2pi f = [A log(vee/Lambda) + (P - A) L_v]/2pi + int_0^inf dp/2pi [f - theta(k-p) A/p - theta(p-k) P/p].
In the y' variables of the file the region P sits at y' -> z, where K(y'-z) is singular: there the limit vee -> infinity
cannot be taken pointwise."""
import numpy as np
import env  # noqa: F401
import fp
import fixed23 as fx
from fp import G, PI, K, sq
from scipy.integrate import quad


def lquad(f, a, b):
    """int_a^b f(p) dp with p = e^t, panels of half an e-fold"""
    g = lambda t: f(np.exp(t)) * np.exp(t)
    ta, tb = np.log(a), np.log(b)
    edges = np.linspace(ta, tb, int(np.ceil((tb - ta) / 0.5)) + 1)
    tot = 0
    for lo, hi in zip(edges[:-1], edges[1:]):
        tot += (quad(lambda t: np.real(g(t)), lo, hi, limit=200, epsabs=0, epsrel=1e-12)[0]
                + 1j * quad(lambda t: np.imag(g(t)), lo, hi, limit=200, epsabs=0, epsrel=1e-12)[0])
    return tot

PI3 = (2 * PI)**3
ROWS = {'III(1)': fx.g3_3_part1_line68, 'V': fx.g3_5_line101, 'VI': fx.g3_6_line120, 'III(2)': fx.g3_3_part2_line80}


def f_w(c, row, k, p, yp0=None):
    """integrand at fixed w' (incl. Jacobian (k/P)^2 and the phase e^{-ik(w'-w)})"""
    P = p + k
    c.pt['yp'] = (k * c.pt['wp'] + p * c.pt['z']) / P
    return (k / P)**2 * ROWS[row](c, k, p) * np.exp(-1j * c.kv @ (c.pt['wp'] - c.pt['w']))


def colour_at(c, row, ypt):
    """colour factor of the row (corrected sign) with y' at the given point, divided out of the kernel"""
    c.pt['yp'] = ypt
    import ft
    return {'III(1)': ft.col_G3_III1, 'V': ft.col_G3_V, 'VI': ft.col_G3_VI, 'III(2)': ft.col_G3_III2}[row](c)


def A_P(c, row, k):
    """closed forms of the limits A (p -> 0, y' = w') and P (p -> inf, y' = z)"""
    pt = c.pt
    ph = np.exp(-1j * c.kv @ (pt['wp'] - pt['w']))
    Kwz, Kyw, Kxz = K(pt['wp'] - pt['z']), K(pt['y'] - pt['w']), K(pt['x'] - pt['z'])
    if row != 'III(2)':
        pref = 1j * G**4 / (8 * PI**4) if row == 'III(1)' else 1j * G**4 / (4 * PI**4)
        A = -pref / PI3 / k * (Kwz @ Kxz) * (Kyw @ K(pt['xp'] - pt['wp'])) * colour_at(c, row, pt['wp'].copy())
        Pl = -pref / PI3 / k * (Kwz @ Kyw) * (Kxz @ K(pt['xp'] - pt['z'])) * colour_at(c, row, pt['z'].copy())
        return A * ph, Pl * ph
    x, z, w = pt['x'], pt['z'], pt['w']
    xz, xw, zw = x - z, x - w, z - w
    A_, B_ = sq(xz), sq(xw)
    pre = G**4 / (4 * PI**4) / PI3
    b0 = K(pt['xp'] - pt['wp'])
    A = -pre / (k * B_) * ((b0 @ xw) * (Kwz @ zw) / sq(zw) + (b0 @ xw) * (Kwz @ xz) / (2 * A_)) * colour_at(c, row, pt['wp'].copy())
    bz = K(pt['xp'] - pt['z'])
    Pl = -pre / (k * A_) * ((bz @ xz) * (Kwz @ zw) / sq(zw) - (bz @ xz) * (Kwz @ xw) / (2 * B_)) * colour_at(c, row, pt['z'].copy())
    return A * ph, Pl * ph


def run(seeds=(1, 2, 3)):
    out = []
    for seed in seeds:
        c = fp.Cfg(seed, sym=True, names=('x', 'y', 'xp', 'yp', 'z', 'w', 'wp'))
        k = 1.0
        for row in ROWS:
            A, Pl = A_P(c, row, k)
            A_num = 1e-8 * f_w(c, row, k, 1e-8)
            P_num = 1e8 * f_w(c, row, k, 1e8)
            # exact int_Lambda^{vee-k} vs the subtracted form
            errs = []
            sub = lambda p: f_w(c, row, k, p) - (A / p if p < k else Pl / p)
            R = (lquad(sub, 1e-12, k) + lquad(sub, k, 1e12)) / (2 * PI)
            for V, L in ((1e3, 1e-3), (1e5, 1e-5)):
                ex = lquad(lambda p: f_w(c, row, k, p), L, V - k) / (2 * PI)
                pred = (A * np.log(V / L) + (Pl - A) * np.log(V / k)) / (2 * PI) + R
                errs.append(float(abs(ex - pred) / abs(ex)))
            out.append(dict(seed=seed, row=row, A=float(abs(A_num - A) / abs(A)), P=float(abs(P_num - Pl) / abs(Pl)),
                            sub_1e3=errs[0], sub_1e5=errs[1]))
            print(seed, '%-7s A rel %.1e  P rel %.1e  exact vs (A logV/L + (P-A) L_v + R): %.1e (vee=1e3), %.1e (vee=1e5)' % (
                row, out[-1]['A'], out[-1]['P'], errs[0], errs[1]))
    return out


if __name__ == '__main__':
    import json
    json.dump(run(), open('crows_w.json', 'w'), indent=1)
