"""All checks of Group_I.tex -> checks1.json.

 1. B_3 in momentum space: instantaneous and two-rho entries vs the light-cone wave functions (wf3.py)
 2. dimensional regularization: (1/pi) int_z [K(x,w,z) - K(w,w,z)] in d = 2 - eps (MS-bar) = the Row II bracket
 3. Row II: part 1 against first principles; the rewriting; the bracket of part 2 against A^(3) of main.tex;
    the split into 2a and 2b
 4. Row III: the integrand (Eq. 90) against first principles, with B_3 as written and with E3;
    every p+ integral (N1, N2, N3a, N3b) against the numerical principal-value integral; the 1/Lambda terms
 5. Row IV: the integrand (Eq. 137) and the p+ integral (Eq. 142)
 6. Row V: part 1 (Eq. 151) and part 2 (Eq. 163) against first principles; the plus prescription of part 1
 7. the four-rho parts of Rows III and IV (not in Group_I.tex) against first principles"""
import json
import numpy as np
import sympy as sp
import fp1
import user1
import wf3
from fp1 import Cfg, PI, G
from integ import pv, plain

NAMES = ('x', 'y', 'xp', 'yp', 'z', 'w', 'wp')
out = {}


def rel(a, b):
    return float(abs(a - b) / max(abs(b), 1e-300))


# 1. B_3 entries
out['B3_vs_wavefunction'] = wf3.run()

# 2. dim reg
eps = sp.symbols('epsilon', positive=True)
r2, mu, gE = sp.symbols('r2 mu gamma_E', positive=True)
d = 2 - eps
I = sp.pi**(2 - d / 2) * sp.gamma(d / 2 - 1) * r2**(1 - d / 2) * (mu**2 * sp.exp(sp.EulerGamma) / (4 * sp.pi))**(eps / 2)
ser = sp.series(I / sp.pi, eps, 0, 1).removeO()
target = -2 / eps - 2 * sp.EulerGamma - sp.log(r2 * mu**2 / 4)
out['dimreg_int_K_x_w_z_over_pi_minus_bracket'] = str(sp.simplify(sp.expand_log(ser - target, force=True)))

# 3. Row II
r = []
for sym in (True, False):
    for seed in (1, 2, 3):
        c = Cfg(seed, sym=sym, names=NAMES)
        b = 2.7 - 0.4j
        r.append(rel(user1.r2_part1(c, 1.3, b), fp1.row2_part1(c, 1.3, b)))
out['RowII_part1_vs_first_principles'] = max(r)
# rewriting: 11/3 + log(L/k) - log(v/L) + log(v/k - 1)  =  11/3 + log(v/k) - 2 log(v/L) + log(v/k - 1)
rng = np.random.default_rng(5)
rw = []
for t in range(5):
    kp, Lam = rng.uniform(0.5, 2), rng.uniform(1e-3, 0.4)
    vee = kp + rng.uniform(0.1, 50)
    a = 11 / 3 + np.log(Lam / kp) - np.log(vee / Lam) + np.log(vee / kp - 1)
    b = 11 / 3 + np.log(vee / kp) - 2 * np.log(vee / Lam) + np.log(vee / kp - 1)
    rw.append(abs(a - b))
out['RowII_rewriting'] = max(rw)


def brackets(kv, pv_, kp, Lam, vee):
    qv = kv - pv_
    K2, P2, Q2, KP = kv @ kv, pv_ @ pv_, qv @ qv, kv @ pv_
    E = (K2 * pv_ - KP * kv) / Q2
    C = ((2 * P2 - 2 * KP + K2) * kv + K2 * pv_) / Q2
    D = (P2 * kv + K2 * pv_) / Q2
    LL, LK = np.log(vee / Lam), np.log(vee / kp)
    main = (np.log(Lam / kp) * E - 2 * kv * LL - 2 * kv * np.log(kp * P2 / (vee * K2 + kp * P2))
            + C * np.log(1 + (vee / kp - 1) * K2 / Q2) - D * np.log((vee - kp) / Lam) - kv * LK)
    notes = ((-LL + LK) * E - 2 * kv * LL - 2 * kv * np.log(kp * P2 / (vee * K2 + kp * P2))
             + C * np.log(1 + (vee / kp - 1) * K2 / Q2) - D * (LL + np.log(1 - kp / vee)))
    Aa = -E - 2 * kv - D
    R = (-2 * kv * np.log(P2 / (K2 + kp / vee * P2)) + C * np.log(K2 / Q2 + kp / vee * (1 - K2 / Q2))
         - D * np.log(1 - kp / vee))
    split = Aa * (LL - LK) + R
    return main, notes, split, kv * LK


bm, bs = [], []
for t in range(6):
    kv, pv_ = rng.normal(size=2), rng.normal(size=2)
    kp, Lam = rng.uniform(0.5, 2), rng.uniform(1e-3, 0.4)
    vee = kp + rng.uniform(0.1, 50)
    main, notes, split, extra = brackets(kv, pv_, kp, Lam, vee)
    bm.append(float(np.abs(notes - main - extra).max() / np.abs(main).max()))
    bs.append(float(np.abs(split - main).max() / np.abs(main).max()))
out['RowII_part2_notes_minus_main_minus_k_logvk'] = max(bm)
out['RowII_part2_split_2a_2b_vs_main'] = max(bs)
# exact bracket (Calculation_main corrections X4, X7; Lambda < k+ < vee - Lambda) and its split
bx = []
for t in range(6):
    kv, pv_ = rng.normal(size=2), rng.normal(size=2)
    kp, Lam = rng.uniform(0.5, 2), rng.uniform(1e-3, 0.4)
    vee = kp + Lam + rng.uniform(0.1, 50)
    qv = kv - pv_
    K2, P2, Q2, KP = kv @ kv, pv_ @ pv_, qv @ qv, kv @ pv_
    E = (K2 * pv_ - KP * kv) / Q2
    C = ((2 * P2 - 2 * KP + K2) * kv + K2 * pv_) / Q2
    D = (P2 * kv + K2 * pv_) / Q2
    LL, LK = np.log(vee / Lam), np.log(vee / kp)
    ex = (np.log(Lam / (kp - Lam)) * E - 2 * kv * LL - 2 * kv * np.log((Lam * K2 + kp * P2) / (vee * K2 + kp * P2))
          + C * np.log(1 + (vee / kp - 1) * K2 / Q2) - C * np.log(1 + Lam / kp * K2 / Q2) - D * np.log((vee - kp) / Lam)
          - kv * np.log(vee / (kp + Lam)))
    Rx = (E * np.log(kp / (kp - Lam)) - 2 * kv * np.log((P2 + Lam / kp * K2) / (K2 + kp / vee * P2))
          + C * np.log(K2 / Q2 + kp / vee * (1 - K2 / Q2)) - C * np.log(1 + Lam / kp * K2 / Q2) - D * np.log(1 - kp / vee)
          + kv * np.log(1 + Lam / kp))
    bx.append(float(np.abs((-E - 2 * kv - D) * (LL - LK) + Rx - ex).max() / np.abs(ex).max()))
out['RowII_part2_exact_split'] = max(bx)
rr = []
for sym in (True, False):
    for seed in (1, 2, 3):
        c = Cfg(seed, sym=sym, names=NAMES)
        kv, pv_, B = rng.normal(size=2), rng.normal(size=2), rng.normal(size=2) + 1j * rng.normal(size=2)
        rr.append(rel(user1.r2_part2(c, 1.3, kv, pv_, B), fp1.row2_part2(c, 1.3, kv, pv_, B)))
out['RowII_part2_prefactor_colour'] = max(rr)

# 4. Row III
pts = []
for sym in (True, False):
    for seed in (1, 2, 3):
        c = Cfg(seed, sym=sym, names=NAMES)
        kp, Lam = 1.3, 0.07
        for pp in (0.4, 0.9, 2.1, 3.7):
            u = sum(user1.r3_lines(c, kp, pp, Lam).values())
            pts.append(dict(sym=sym, seed=seed, pp=pp, written=rel(u, fp1.row3(c, kp, pp, Lam, sgt=+1)),
                            E3=rel(u, fp1.row3(c, kp, pp, Lam, sgt=-1))))
out['RowIII_integrand'] = pts

ints = []
for sym, seed, kp, Lam, vee in ((True, 1, 1.3, 0.07, 9.0), (True, 2, 1.3, 0.07, 9.0), (True, 3, 1.3, 0.07, 9.0),
                                (True, 5, 1.3, 0.07, 9.0), (False, 1, 0.8, 0.2, 3.0), (False, 6, 2.0, 0.01, 40.0)):
    c = Cfg(seed, sym=sym, names=NAMES)
    q = user1.r3_kin(c)
    A_, B_ = q['A'], q['B']
    ps = kp * B_ / A_
    N = user1.r3_int(c, kp, Lam, vee)
    Nc = user1.r3_int(c, kp, Lam, vee, absl=False)
    ranges = {'L1': (kp + Lam, vee), 'L2': (Lam, kp - Lam), 'L3': (Lam, vee), 'L4': (kp + Lam, vee), 'L5': (Lam, kp - Lam)}
    groups = {'L1': ['N1.1', 'N1.2', 'N1.3'], 'L2': ['N1.4', 'N1.5', 'N1.6'], 'L3': ['N1.7', 'N1.8'],
              'L4': [k for k in N if k.startswith('N2')], 'L5': [k for k in N if k.startswith('N3')]}
    row = dict(sym=sym, seed=seed, kp=kp, Lam=Lam, vee=vee, pole=float(ps))
    for L, (a, b) in ranges.items():
        g = lambda p, L=L: user1.r3_lines(c, kp, p, Lam, strip=True)[L]
        num = pv(g, a, b, ps, A_) / (2 * PI)
        row[L] = rel(sum(N[k] for k in groups[L]), num)
        row[L + '_complex_log'] = rel(sum(Nc[k] for k in groups[L]), num)
        row[L + '_pole_inside'] = bool(a < ps < b)
    # with E3 the non-instantaneous p+ > k+ entry changes sign: compare -N2 with the first-principles integral
    gfp = lambda p: fp1.row3(c, kp, p, Lam, sgt=-1, parts=('ni',)) * (p * A_ - kp * B_)
    row['N2_E3'] = rel(-sum(N[k] for k in groups['L4']), pv(gfp, kp + Lam, vee, ps, A_) / (2 * PI))
    # rational terms (1/Lambda, 1/(vee-k), 1/(k-Lambda))
    rat_w = N['N1.1'] + N['N1.4'] + N['N2.1'] + N['N3a.1']
    rat_E3 = N['N1.1'] + N['N1.4'] - N['N2.1'] + N['N3a.1']
    row['rational_written'] = complex(rat_w)
    row['rational_written_expected'] = complex(2j * G**4 / (8 * PI**5) * q['KK1'] / (A_ - B_) * (1 / Lam - 1 / (vee - kp))
                                               * user1.col3(c) * fp1.phase(c) / fp1.PI3)
    row['rational_E3'] = abs(rat_E3) / abs(rat_w)
    ints.append(row)
out['RowIII_p_integrals'] = ints

# 5. Row IV
r4 = []
for sym in (True, False):
    for seed in (1, 2, 3):
        c = Cfg(seed, sym=sym, names=NAMES)
        kp, Lam, vee = 1.3, 0.07, 9.0
        pw = max(rel(user1.r4_line(c, kp, pp), fp1.row4(c, kp, pp)) for pp in (0.2, 1.0, 3.5, 7.0))
        num = plain(lambda p: user1.r4_line(c, kp, p), Lam, vee - kp) / (2 * PI)
        r4.append(dict(sym=sym, seed=seed, integrand=pw, integral=rel(user1.r4_int(c, kp, Lam, vee)[0], num)))
out['RowIV'] = r4

# 6. Row V
r5 = []
for sym in (True, False):
    for seed in (1, 2, 3):
        c = Cfg(seed, sym=sym, names=NAMES)
        kp = 1.3
        a1, a2 = [], []
        for pp in (0.15, 0.6, 1.1):
            c.pt['w'] = c.pt['z'] - pp / kp * (c.pt['z'] - c.pt['x'])
            a1.append(rel(user1.r5a_line(c, kp, pp), fp1.row5_part1(c, kp, pp)))
            a2.append(rel(user1.r5b_line(c, kp, pp), fp1.row5_part2(c, kp, pp)))
        r5.append(dict(sym=sym, seed=seed, part1=max(a1), part2=max(a2)))
out['RowV_integrands'] = r5
pl = []
for seed in (1, 2, 3):
    c = Cfg(seed, sym=True, names=NAMES)
    kp = 1.3
    for Lam in (0.1, 0.01, 0.001):
        ex = user1.r5a_exact(c, kp, Lam, plain)
        pl.append(dict(seed=seed, Lam=Lam, written_0_k=rel(user1.r5a_plus(c, kp, Lam, 0.0, kp, plain), ex),
                       limits_Lam=rel(user1.r5a_plus(c, kp, Lam, Lam, kp - Lam, plain), ex)))
out['RowV_part1_plus'] = pl
pl2 = []
for seed in (1, 2, 3):
    c = Cfg(seed, sym=True, names=NAMES)
    kp = 1.3
    for Lam in (0.1, 0.01, 0.001):
        ex = user1.r5_2a_exact(c, kp, Lam, plain)
        pl2.append(dict(seed=seed, Lam=Lam, written_0_k=rel(user1.r5_2a_plus(c, kp, Lam, 0.0, kp, plain), ex),
                        limits_Lam=rel(user1.r5_2a_plus(c, kp, Lam, Lam, kp - Lam, plain), ex)))
out['RowV_2a_plus'] = pl2

# 7. four-rho parts
f4 = []
for sym in (True, False):
    c = Cfg(4, sym=sym, names=NAMES)
    for pp in (0.3, 2.2):
        f4.append(dict(sym=sym, pp=pp, RowIII=rel(user1.r3_4rho(c, 1.3, pp), fp1.row3_4rho(c, 1.3, pp)),
                       RowIV=rel(user1.r4_4rho(c, 1.3, pp), fp1.row4_4rho(c, 1.3, pp))))
out['four_rho'] = f4


def tojson(o):
    if isinstance(o, dict):
        return {k: tojson(v) for k, v in o.items()}
    if isinstance(o, list):
        return [tojson(v) for v in o]
    if isinstance(o, complex):
        return [o.real, o.imag]
    if isinstance(o, (np.floating, np.integer, np.bool_)):
        return o.item()
    return o


if __name__ == '__main__':
    json.dump(tojson(out), open('checks1.json', 'w'), indent=1)
    print(json.dumps(tojson({k: v for k, v in out.items() if k not in ('RowIII_integrand', 'B3_vs_wavefunction')}), indent=1)[:6000])
