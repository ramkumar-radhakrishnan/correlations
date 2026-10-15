"""Runs every check of Main_corrections.tex and writes main_checks.json (read by make_note.py).
Takes about two minutes (the symbolic dim-reg integral of onerho_dimreg.psi3_onerho_xi)."""
import json
import os
import subprocess
import sys
import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lcpt1
import lcpt_vac as LV
import b3_check
import b2mom_check
import tworho2_check
import onerho_dimreg as OD

R = {}


def stat(rs):
    rs = np.array(rs)
    return {'mean_real': float(np.mean(rs.real)), 'mean_imag': float(np.mean(rs.imag)), 'spread': float(np.max(np.abs(rs - rs.mean())))}


def ratio_all(a, b):
    m = np.abs(b) > 1e-8 * np.abs(b).max()
    return a[m] / b[m]


# 1. one gluon -> one gluon, non-instantaneous (onegrho12, onegrho22 and their xi forms)
rng = np.random.default_rng(3)
acc = {k: [] for k in ('onegrho12', 'onegrho12 xi-form', 'onegrho22 (typo read as product)', 'onegrho22 xi-form')}
for t in range(3):
    k, p, v = rng.normal(size=2), rng.normal(size=2), rng.normal(size=8)
    kp = 1.3
    pp = kp * rng.uniform(1.2, 3.0)
    a = lcpt1.psi1(kp, k, pp, p, v)
    acc['onegrho12'] += list(ratio_all(lcpt1.user_onegrho12(kp, k, pp, p, v), a))
    acc['onegrho12 xi-form'] += list(ratio_all(lcpt1.user_onegrho1_xi(kp, k, pp, p, v), a))
    pp = kp * rng.uniform(0.1, 0.9)
    a = lcpt1.psi2(kp, k, pp, p, v)
    acc['onegrho22 (typo read as product)'] += list(ratio_all(lcpt1.user_onegrho22(kp, k, pp, p, v), a))
    acc['onegrho22 xi-form'] += list(ratio_all(lcpt1.user_onegrho2_xi(kp, k, pp, p, v), a))
R['written / first principles (one gluon -> one gluon)'] = {k: stat(v) for k, v in acc.items()}

# 2. B3 momentum space, non-instantaneous entries
rng = np.random.default_rng(4)
gt, lt = [], []
for t in range(3):
    k, p, v = rng.normal(size=2), rng.normal(size=2), rng.normal(size=8)
    kp = 1.3
    pp = kp * rng.uniform(1.2, 3)
    gt += list(ratio_all(b3_check.b3_noninst(kp, k, pp, p, v, 'gt'), lcpt1.psi1(kp, k, pp, p, v)))
    pp = kp * rng.uniform(0.1, 0.9)
    lt += list(ratio_all(b3_check.b3_noninst(kp, k, pp, p, v, 'lt'), lcpt1.psi2(kp, k, pp, p, v)))
R['B3 / wave function (8 pi^3 = %.8f)' % (8 * np.pi**3)] = {'Theta(p+-k+-Lambda) entry': stat(gt), 'Theta(k+-Lambda-p+) entry': stat(lt)}

# 3. vacuum one-gluon terms at order g^3
rng = np.random.default_rng(11)
acc = {k: [] for k in ('instgg', 'eq:gg-uu', '3grhorho', 'psif', 'psiff', '4+5tworho')}
for t in range(3):
    kp = 1.7
    k, p = rng.normal(size=2), rng.normal(size=2)
    P = (kp * rng.uniform(0.1, 0.9), p)
    K = (kp, k)
    Q = LV.sub(K, P)
    acc['instgg'] += list(ratio_all(LV.w_instgg(K, P), LV.fp_psi2a(K, P)))
    acc['eq:gg-uu'] += list(ratio_all(LV.w_gguu(K, P), LV.fp_psi2b(K, P)))
    tot = LV.fp_psi3(K, P)[('p', 'q')] + LV.fp_psi3(K, Q)[('q', 'p')]
    acc['3grhorho'] += list(ratio_all(LV.w_3grhorho(K, P), tot))
    f4, f5 = LV.fp_psi4(K, P), LV.fp_psi5(K, P)
    acc['psif'] += list(ratio_all(LV.w_psif(K, P), f4))
    acc['psiff'] += list(ratio_all(LV.w_psiff(K, P), f5))
    S, C = LV.split(f4 + f5)
    acc['4+5tworho'] += list(ratio_all(LV.w_45tworho(K, P), S))
R['written / first principles (vacuum, one gluon, g^3)'] = {k: stat(v) for k, v in acc.items()}

# 4. B2 momentum space
rng = np.random.default_rng(2)
acc = {'one-rho': [], 'two-rho': []}
for t in range(3):
    kp, pp = rng.uniform(0.3, 2, size=2)
    k, p = rng.normal(size=2), rng.normal(size=2)
    one_w, two_w, _ = b2mom_check.psi(kp, k, pp, p)
    one_b, two_b, _ = b2mom_check.B2(kp, k, pp, p)
    acc['one-rho'] += list(ratio_all(one_b, one_w))
    acc['two-rho'] += list(ratio_all(two_b, two_w))
R['B2 (momentum) / wave function (8 pi^3)'] = {k: stat(v) for k, v in acc.items()}

# 5. one-rho parts in dim reg
xi, eps, L2 = OD.xi, OD.eps, OD.L2
Fa = OD.psi3_onerho_xi()
Fa_claim = ((2 - xi) * sp.log(xi) + (1 + xi) * sp.log(1 - xi) + (1 + xi) * (-2 / eps + L2)) / (xi * (1 - xi))
diff = [complex(sp.N((Fa - Fa_claim).subs({xi: xv, eps: sp.Rational(3, 7), L2: sp.Rational(5, 9)}))) for xv in (sp.Rational(1, 5), sp.Rational(1, 2), sp.Rational(7, 10))]
R['Psi^3 one-rho xi-integrand (dim reg) minus closed form'] = {'max_abs': max(abs(d) for d in diff)}
rows = []
for xv in (0.2, 0.5, 0.7):
    pol, lm, fi = OD.psi45_onerho_x(xv)
    nt = (xv - 2) * np.log(xv) + (xv + 1) * np.log(1 - xv)
    rows.append({'xi': xv, 'pole': pol, 'pole_notes': 2 * (xv - 2), 'ln_mu': lm, 'ln_mu_notes': -(xv - 2), 'finite': fi, 'finite_notes': nt})
R['4+5onerho (dim reg) vs notes'] = rows

# 6. tworho2 (Lambda -> 0 form) vs exact xi integral, antisymmetrised
rng = np.random.default_rng(5)
rows = []
for d in [1e-2, 1e-4, 1e-6]:
    k, p = rng.normal(size=2), rng.normal(size=2)
    q = k - p
    v = tworho2_check.V_num(k, p, d) / (q @ q) - tworho2_check.V_num(k, q, d) / (p @ p)
    n = tworho2_check.V_notes(k, p, d) / (q @ q) - tworho2_check.V_notes(k, q, d) / (p @ p)
    rows.append({'delta': d, 'max_abs_diff': float(np.abs(v - n).max()), 'max_abs_exact': float(np.abs(v).max())})
R['tworho2 (as written) vs exact'] = rows

# 7. exact forms
subprocess.run([sys.executable, os.path.join(HERE, 'exact_forms.py')], check=True, stdout=subprocess.DEVNULL)
R['exact forms'] = json.load(open(os.path.join(HERE, 'exact_forms.json')))

json.dump(R, open(os.path.join(HERE, 'main_checks.json'), 'w'), indent=1)
print(json.dumps(R, indent=1)[:4000])
