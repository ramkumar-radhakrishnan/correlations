"""All checks of the note Groups_II_III_check.pdf -> checks23.json.
Uses fp.py (first principles from Calculation_main.tex), user23*.py (the equations as written), fixed23.py (corrected placements)."""
import json
import numpy as np
import sympy as sp
from scipy import integrate
import fp
import user23 as us
import user23b as ub
import user23c as uc
import fixed23 as fx
from fp import F, NC, G, PI, K, sq
from run1 import fp_rows, user_rows

k, Lam, vee = 1.7, 0.11, 23.0
out = {}
FIX = {'G2-I(1) A': fx.g2_1_termA, 'G2-I(1) B': fx.g2_1_termB, 'G2-I(2)': fx.g2_2_line54, 'G3-I': fx.g3_1_line44,
       'G3-II': fx.g3_2_line55, 'G3-III(1)': fx.g3_3_part1_line68, 'G3-III(2)': fx.g3_3_part2_line80,
       'G3-V': fx.g3_5_line101, 'G3-VI': fx.g3_6_line120, 'G2-IV': fx.g2_4_line92}

# 1. pre-integration lines
pre = {}
for key in FIX:
    pre[key] = {'written_sym': [], 'written_gen': [], 'corrected_gen': []}
for seed in (11, 12, 13):
    for sym in (True, False):
        c = fp.Cfg(seed, sym)
        p = 0.37 + 0.3 * seed % 1
        a, b = fp_rows(c, p), user_rows(c, p)
        for key in a:
            if sym:
                pre[key]['written_sym'].append(complex(b[key] / a[key]))
            else:
                pre[key]['written_gen'].append(complex(b[key] / a[key]))
                pre[key]['corrected_gen'].append(complex(FIX[key](c, k, p) / a[key]))
out['pre'] = {kk: {t: [[round(z.real, 8), round(z.imag, 8)] for z in v] for t, v in d.items()} for kk, d in pre.items()}


def pint(f, a=Lam, b=vee - k):
    re = integrate.quad(lambda p: np.real(f(p)), a, b, limit=800, epsabs=0, epsrel=1e-11)[0]
    im = integrate.quad(lambda p: np.imag(f(p)), a, b, limit=800, epsabs=0, epsrel=1e-11)[0]
    return re + 1j * im


# 2. integrated lines of Group II Row I and Group III Rows I, II (general U, corrected colour)
integ = {}
for seed in (21, 22):
    c = fp.Cfg(seed, False)
    A_ex = pint(lambda p: fx.g2_1_termA(c, k, p)) / (2 * PI); B_ex = pint(lambda p: fx.g2_1_termB(c, k, p)) / (2 * PI)
    A_w, B_w = ub.g2_1_line50(c, k, Lam, vee)
    integ.setdefault('G2-I(1) line 50', []).append(max(abs(A_w / A_ex - 1), abs(B_w / B_ex - 1)))
    wts = {'2a': lambda p: (p + k)**2 / (p * k), '2b': lambda p: (p + k)**2 / p**2, '2c': lambda p: (p + k)**2 / k**2,
           '2d': lambda p: (p + k) / p, '2e': lambda p: (p + k) / k, '2f': lambda p: 1.0}
    gname = {'2a': 'pk', '2b': 'pp', '2c': 'kk', '2d': 'p', '2e': 'k', '2f': 'd'}
    written = ub.g2_2_integrated(c, k, Lam, vee)
    for key in wts:
        ex = pint(lambda p: us.g2_2_pref(c, k, p) * wts[key](p) * us.g2_2_groups62(c, k, p)[gname[key]]) / (2 * PI)
        integ.setdefault('G2-I(2) ' + key, []).append(abs(written[key] / ex - 1))
    pp = 0.63
    kern62 = us.g2_2_line62(c, k, pp) / us.g2_2_colour58(c)
    kern54 = us.g2_2_line54(c, k, pp) / us.g2_2_colour54(c)
    integ.setdefault('G2-I(2) line 62 kernel / line 54 kernel', []).append(abs(kern62 / kern54 - 1))
    integ.setdefault('G2-I(2) colour 58 / colour 54 (corrected placements)', []).append(abs(fx.g2_2_colour58(c) / fx.g2_2_colour(c) - 1))
    L1 = np.sum((np.einsum('aB,CD,Cb,B,D->ab', c.U['xp'], c.U['y'], c.U['z'], c.rho['xp'], c.rho['y'])
                 - np.einsum('Db,DE,aC,C,E->ab', c.U['z'], c.U['y'], c.U['wp'], c.rho['xp'], c.rho['y'])) * fx.Rone_fix(c, 'w', 'z', 'x'))
    ex1 = pint(lambda p: fx.g3_1_line44(c, k, p)) / (2 * PI)
    integ.setdefault('G3-I line 48 (G3-II line 59 has the same kernel)', []).append(abs(ub.g3_integrated(c, k, Lam, vee) * L1 / ex1 - 1))
out['integrated'] = {kk: max(v) for kk, v in integ.items()}

# 3. C rows: w' integration and plus prescription (U symmetric = your convention; kernels are colour independent)
crow = {}
for seed in (31, 32, 33):
    c = fp.Cfg(seed, True)
    p0 = 0.63
    col68 = np.einsum('CDa,CDb,ab->', F, fx.Lc_fix(c), fx.two_fix(c, 'w', 'y', 'x'))
    col80 = np.einsum('CDa,CDb,ab->', F, fx.Lc_fix(c), fx.Rone_fix(c, 'w', 'z', 'x'))
    col101 = np.einsum('CDa,CDb,ab->', F, fx.Lc_fix(c), fx.X1amp_fix(c))
    col120 = np.einsum('CDa,CDb,ab->', F, fx.Lc_fix(c), fx.X2amp_fix(c))
    c.pt['wp'] = ((p0 + k) * c.pt['yp'] - p0 * c.pt['z']) / k
    rows = {'III(1) line 72': (fx.g3_3_part1_line68(c, k, p0, sign=1), uc.ker72(c, k, p0, -1j * G**4 / (8 * PI**4)) * col68),
            'III(2) line 85': (fx.g3_3_part2_line80(c, k, p0, sign=1),
                               (uc.ker85_pieces(c, k, p0)[0] + uc.ker85_pieces(c, k, p0)[1] / p0) / uc.phase(c, k, p0) * col80),
            'V line 105': (fx.g3_5_line101(c, k, p0, sign=1), uc.ker72(c, k, p0, -1j * G**4 / (4 * PI**4)) * col101),
            'VI line 124': (fx.g3_6_line120(c, k, p0, sign=1), uc.ker72(c, k, p0, -1j * G**4 / (4 * PI**4)) * col120)}
    for key, (a, b) in rows.items():
        crow.setdefault('w\' integration ' + key, []).append(abs(b / a - 1))
    for key, pref, plog in (('III(1) line 76', -1j * G**4 / (8 * PI**4), 1j * G**4 / (16 * PI**5)),
                            ('V line 109', -1j * G**4 / (4 * PI**4), 1j * G**4 / (8 * PI**5)),
                            ('VI line 128', -1j * G**4 / (4 * PI**4), 1j * G**4 / (8 * PI**5))):
        pieces = lambda p: uc.ker72_pieces(c, k, p, pref)
        exact = pint(lambda p: pieces(p)[0] + pieces(p)[1] / p) / (2 * PI)
        g10 = pieces(0.0)[1]
        sub = lambda p: pieces(p)[0] + (pieces(p)[1] - g10) / p
        logt = uc.logterm72(c, k, Lam, vee, plog)
        w0 = pint(sub, 1e-14) / (2 * PI) + logt
        wL = pint(sub, Lam) / (2 * PI) + logt
        crow.setdefault('plus prescription ' + key + ': written (from 0)', []).append(abs(w0 / exact - 1))
        crow.setdefault('plus prescription ' + key + ': from Lambda', []).append(abs(wL / exact - 1))
    exact = pint(lambda p: uc.ker85_pieces(c, k, p)[0] + uc.ker85_pieces(c, k, p)[1] / p) / (2 * PI)
    g10 = uc.ker85_pieces(c, k, 0.0)[1]
    sub = lambda p: uc.ker85_pieces(c, k, p)[0] + (uc.ker85_pieces(c, k, p)[1] - g10) / p
    logt = uc.logterm90(c, k, Lam, vee)
    crow.setdefault('plus prescription III(2) line 90: written (from 0)', []).append(abs((pint(sub, 1e-14) / (2 * PI) + logt) / exact - 1))
    crow.setdefault('plus prescription III(2) line 90: from Lambda', []).append(abs((pint(sub, Lam) / (2 * PI) + logt) / exact - 1))
out['crows'] = {kk: [min(v), max(v)] for kk, v in crow.items()}

# 4. Row V Wilson-line simplification and Row IV lines 98, 104 (general U, corrected placements)
r4 = {}
for seed in (41, 42):
    c = fp.Cfg(seed, False)
    U, r = c.U, c.rho
    col109f = np.einsum('CDa,CDb,ab->', F, fx.Lc_fix(c), fx.X1amp_fix(c))
    L113f = np.einsum('CE,E->C', U['yp'], r['xp']) - np.einsum('CE,E->C', U['xp'], r['xp'])
    R113 = np.einsum('ce,ad,e,d->ca', U['x'], U['y'], r['x'], r['y']) - np.einsum('ad,ce,e,d->ca', U['w'], U['x'], r['x'], r['y'])
    r4.setdefault('Row V: colour of line 113 / line 109', []).append(abs(np.einsum('Cca,C,ca->', F, L113f, R113) / col109f - 1))
    p = 0.63; P = p + k
    c.pt['w'] = (P * c.pt['y'] - p * c.pt['z']) / k; c.pt['wp'] = (P * c.pt['yp'] - p * c.pt['z']) / k
    l92 = fx.g2_4_line92(c, k, p) * (P / k)**4
    I2 = np.eye(2)
    br = (2 * p / P**2 * np.einsum('km,KM->kmKM', I2, I2) + 2 / k * (np.einsum('mK,Mk->kmKM', I2, I2) - np.einsum('MK,km->kmKM', I2, I2))
          + (p / k**2 + 1 / p) * np.einsum('mM,kK->kmKM', I2, I2))
    kern98 = np.einsum('m,k,M,K,kmKM->', K(c.pt['y'] - c.pt['z']), K(c.pt['x'] - c.pt['y']), K(c.pt['yp'] - c.pt['z']), K(c.pt['xp'] - c.pt['yp']), br)
    Lp = np.einsum('CDa,CDb->ab', F, fx.Lc_fix(c))
    Ra = (np.einsum('cda,ce,db,e->ab', F, U['y'], U['z'], r['x']) - np.einsum('cda,db,ce,e->ab', F, U['z'], U['x'], r['x']))
    col92 = np.sum(Lp * Ra)
    r4.setdefault('Row IV: line 98 / (line 92 at w*, w\'* times (P/k)^4)', []).append(abs(1 / us.pi3 * G**4 / (4 * PI**4) / k * kern98 * col92 / l92 - 1))
    col104f = NC * np.einsum('cE,E,ce,e->', U['yp'] - U['xp'], r['xp'], U['y'] - U['x'], r['x'])
    r4.setdefault('Row IV: colour of line 104 / line 98', []).append(abs(col104f / col92 - 1))
out['row4'] = {kk: max(v) for kk, v in r4.items()}

# 5. z integral: traceless part (2 dims) and assembly of the corrected Row IV (d = 2 - eps)
rng = np.random.default_rng(5)
y, yp = rng.normal(size=2), rng.normal(size=2)
s = y - yp
cen = (y + yp) / 2
num = np.zeros((2, 2))
for a in range(2):
    for b in range(2):
        def f(th, rr):
            z = cen + rr * np.array([np.cos(th), np.sin(th)])
            u, v = y - z, yp - z
            T = np.outer(u, v) / (sq(u) * sq(v))
            return (T - np.trace(T) / 2 * np.eye(2))[a, b] * rr
        edges = [0, 0.5 * np.sqrt(sq(s)), np.sqrt(sq(s)), 3, 30, 300, 3000]
        num[a, b] = sum(integrate.dblquad(f, lo, hi, 0, 2 * np.pi, epsabs=1e-10, epsrel=1e-10)[0] for lo, hi in zip(edges[:-1], edges[1:]))
pred = -np.pi * (np.outer(s, s) / sq(s) - np.eye(2) / 2)
out['J_traceless_maxdiff'] = float(np.abs(num - pred).max())
out['J_symmetric'] = float(abs(num[0, 1] - num[1, 0]))
# assembly: sum over polarizations with J^{ab} = A delta + B s s, A = trace/2, trace = -pi(2/eps + C), B = -pi
asm = []
for trial in range(4):
    eps, C, p = 1e-5, rng.normal(), rng.uniform(0.2, 3)
    P = p + k
    d = 2 - eps
    Kx, Kxp, sh = rng.normal(size=2), rng.normal(size=2), rng.normal(size=2)
    sh = sh / np.sqrt(sq(sh))
    tr = -np.pi * (2 / eps + C)
    Bt = -np.pi * (np.pi**(eps / 2)) * 0 - np.pi
    A_ = (tr - Bt) / d
    J = A_ * np.eye(2) + Bt * np.outer(sh, sh)
    lhs = d * p / P**2 * (Kx @ J @ Kxp) + 2 / k * (Kxp @ J @ Kx - Kx @ J @ Kxp) + (p / k**2 + 1 / p) * (Kx @ Kxp) * tr
    rhs = -np.pi * ((Kx @ Kxp) * ((p / P**2 + p / k**2 + 1 / p) * (2 / eps + C) - p / P**2) + 2 * p / P**2 * (Kx @ sh) * (Kxp @ sh))
    asm.append(abs(lhs - rhs) / abs(rhs))
out['row4_assembly_reldiff'] = max(asm)
xi, pp_, kk_ = sp.symbols('xi p k', positive=True)
phih = xi * (1 - xi) + xi / (1 - xi) + (1 - xi) / xi
out['row4_phihat_identity'] = str(sp.simplify((pp_ / (pp_ + kk_)**2 + pp_ / kk_**2 + 1 / pp_) * kk_ - phih.subs(xi, kk_ / (pp_ + kk_))))
lam_, v_, k_ = sp.symbols('Lambda vee k', positive=True)
x0, x1 = k_ / v_, k_ / (k_ + lam_)
out['row4_endpoint_log'] = str(sp.simplify(sp.expand_log(sp.log((1 - x0) / (1 - x1)) - (sp.log(v_ / lam_) - sp.log(v_ / k_) + sp.log(1 - k_ / v_) + sp.log(1 + lam_ / k_)), force=True)))
# LO normalisation of the pole: (g^4 Nc/(8pi^4))(1/(k xi^2)) Phihat = (alpha_s/2pi) 2 Nc Phihat (1/xi^3) (g^2 xi/(2 pi^2 k))
g_, Nc_, k3, x3 = sp.symbols('g N_c k xi', positive=True)
out['row4_LO_normalisation'] = str(sp.simplify(g_**4 * Nc_ / (8 * sp.pi**4) / (k3 * x3**2)
                                               - g_**2 / (8 * sp.pi**2) * 2 * Nc_ / x3**3 * g_**2 * x3 / (2 * sp.pi**2 * k3)))
json.dump(out, open('checks23.json', 'w'), indent=1)
print(json.dumps({kk: v for kk, v in out.items() if kk != 'pre'}, indent=1))
for kk, d in out['pre'].items():
    print(kk, {t: v[0] for t, v in d.items()})
