#!/usr/bin/env python3
"""dilute2 / e3: assembly of the dilute Group I constant C(q_a,q_b) (point quark; same procedure as ref3_dilute/v9b_final.py)
at the corners of the grid {q_a1,q_a2,q_a3} x {q_b1,q_b2}, and the F_P(q_a)+F_T(q_b) test against the one-loop
Lipatov vertex (BSV hep-ph/0608154 eq:kernelv = Fadin-Lipatov virtual kernel; cross-checked here against
Del Duca-Schmidt hep-ph/9810215 eq:looplip at tau = k^2).

Units g = 1, N_c = 3:  abar = alpha_s N_c/pi = N_c/(4 pi^2).  mu = |k| at every corner (C and X are scale invariant).
D(q_a,q_b) = C - X + (abar/2) ln(k^2) ln(q_a^2/q_b^2),   X = (abar/4){pi^2 + ... } finite part of K^(v)/K^B
(the double pole -2/eps^2 (k^2/mu^2)^eps of K^(v) is cancelled, logarithms included, by the real soft-collinear term,
see the note).  Test: double difference of D over rectangles in (q_a,q_b) = 0.
usage: e3_assemble.py <dir with new outputs> <dir ref3_dilute>"""
import sys, os, numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import curve_fit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from e0_points import PTS
S, S0 = sys.argv[1], sys.argv[2]
Nc = 3; abar = Nc / (4 * np.pi ** 2); c8 = Nc / (8 * np.pi ** 2)
OLD = {'a1b1': ('pt1', [(1e-5, -2.46802238), (1e-6, -3.57719450), (1e-7, -4.89472033), (1e-8, -6.42073748)], -0.04492649),
       'a2b2': ('pt2', [(1e-5, 5.60561281), (1e-6, 8.19570618), (1e-7, 11.27714490), (1e-8, 14.84983799)], -0.00743644)}

def load(name):
    if name in OLD:
        tag, III, Ainf = OLD[name]
        return (np.loadtxt(f'{S0}/v4b_sgrid_{tag}.out'), np.loadtxt(f'{S0}/v5_endpoints_{tag}.out'), np.loadtxt(f'{S0}/v5_nodes_{tag}.out'), III, Ainf)
    III = []; Ainf = []
    for Lam in [L_ for L_ in ('1e-5', '1e-6', '1e-7', '1e-8') if os.path.exists(f'{S}/h2_{name}_{L_}.out')]:
        import re
        m = re.search(r'R_III\(vee=1\)=(\S+)\s+coeff ln vee=(\S+)', open(f'{S}/h2_{name}_{Lam}.out').read())
        III.append((float(Lam), float(m.group(1)))); Ainf.append(float(m.group(2)))
    return (np.loadtxt(f'{S}/v4b_sgrid_{name}.out'), np.loadtxt(f'{S}/v5_endpoints_{name}.out'), np.loadtxt(f'{S}/v5_nodes_{name}.out'), III, float(np.mean(Ainf)))

def X_bsv(qa, qb, with_const=True):
    """finite part of K^(v)/K^B (BSV eq:kernelv, n_f=0, double-pole term and beta_0 pole removed), units abar/4."""
    D2 = (qa + qb) @ (qa + qb); a2, b2 = qa @ qa, qb @ qb; ab = qa @ qb; dqq = (qa - qb) @ (qa - qb)
    def X(a2, b2):
        Lab = np.log(a2 / b2); den = a2 - b2
        br = D2 / den * (1 - D2 * (a2 + b2 - 4 * ab) / (3 * den ** 2)) * Lab - D2 * dqq / (6 * a2 * b2) + D2 ** 2 * (a2 + b2) * (a2 + b2 - 4 * ab) / (6 * a2 * b2 * den ** 2)
        return (np.pi ** 2 if with_const else 0.0) + 3 * D2 / den * Lab - Lab ** 2 + br
    if abs(a2 - b2) < 1e-6 * a2:
        xs = [X(a2 * (1 + h), b2) for h in (1e-3, 2e-3, -1e-3, -2e-3)]
        return (4 * (xs[0] + xs[2]) / 2 - (xs[1] + xs[3]) / 2) / 3
    return X(a2, b2)

def X_dds(qa, qb):
    """2 Re of DDS eq:looplip / N_c at tau = k^2, n_f = 0, without poles and pi^2 constants. DDS: k = q_a - q_b^DDS, q_b^DDS = -q_b."""
    t1, t2 = -(qa @ qa), -(qb @ qb); k2 = (qa + qb) @ (qa + qb); re = -(t1 + t2 + k2) / 2
    x = t1 / t2; L0 = np.log(x) / (1 - x); L2 = (np.log(x) - (x - 1 / x) / 2) / (1 - x) ** 3
    v = -0.5 * np.log(x) ** 2 - (11 / 6) * (t1 + t2 + 2 * re) * L0 / t2 + (1 / 3) * k2 * (-(2 * t1 * t2 + (t1 + t2 + 2 * k2) * re) * L2 / t2 ** 3 - re / (2 * t1 * t2))
    return 2 * v

def groupI(name, verbose=True):
    P = PTS[name]; k = np.array([1.0, 0.0]); l = P['l']; q = k + l
    d4, e, b, III, Ainf = load(name)
    v0 = k / (k @ k) - q / (q @ q); vv = v0 @ v0; Pk = (v0 @ k) / vv; Pq = (v0 @ q) / (q @ q) / vv
    out = {}
    out['I.II comm'] = dict(c2=abar, cL=abar * Pq * np.log(q @ q), cV=0, C=Nc / 24)
    u, sG = d4[:20, 0], d4[:20, 1]; g0, gi = d4[20, 1], d4[21, 1]
    cs = CubicSpline(u, sG); out['I.IV'] = dict(c2=0, cL=g0, cV=gi, C=cs.integrate(-12, 0) - 12 * g0 + cs.integrate(0, 12) - 12 * gi)
    xe, te = e[:, 0], e[:, 1]
    f = lambda x, A, B, C1, D1: A + B * np.log(x) + C1 * x * np.log(x) + D1 * x
    (A, B, _, _), _ = curve_fit(f, xe, xe * te, p0=[-0.01, -0.02, 0, 0])
    xb, tb = b[:, 0], b[:, 1]; xg, wg = np.polynomial.legendre.leggauss(len(xb)); tg = 0.5 * (xg + 1)
    Ifin = np.sum(wg * 0.5 * tg * (tb - (A + B * np.log(xb)) / xb))
    out['I.V res'] = dict(c2=-B, cL=2 * A, cV=0, C=2 * (Ifin - A * np.log(2) + 0.5 * B * np.log(2) ** 2))
    La = np.array([x for x, _ in III]); R = np.array([y for _, y in III]); Ls = np.log(1 / La)
    # (a) free quadratic fit; (b) L^2 coefficient fixed to its analytic value -(Nc/8pi^2) P_q; (c) subsets of three with (b)
    M = np.vstack([Ls ** 2, Ls, np.ones_like(Ls)]).T; (a2, a1, a0), *_ = np.linalg.lstsq(M, R, rcond=None)
    a2x = -c8 * Pq
    Mb = np.vstack([Ls, np.ones_like(Ls)]).T; (b1, b0), *_ = np.linalg.lstsq(Mb, R - a2x * Ls ** 2, rcond=None)
    sub = []
    for drop in range(len(Ls)):
        ii = [i for i in range(len(Ls)) if i != drop]; (s1, s0), *_ = np.linalg.lstsq(Mb[ii], (R - a2x * Ls ** 2)[ii], rcond=None); sub.append((s1, s0))
    rest = {kk: sum(o[kk] for o in out.values()) for kk in ('c2', 'cL', 'cV', 'C')}
    two_wb = -abar * np.log(l @ l); two_wa = -abar * np.log(q @ q)
    # (d) constrained: L^2 analytic and total L coefficient = 2 omega(q_b) imposed -> constant from the mean
    a1c = two_wb - rest['cL']; a0c = np.mean(R - a2x * Ls ** 2 - a1c * Ls)
    res = dict(name=name, Pq=Pq, c2_free=rest['c2'] + a2, c2_fit=a2, c2_ana=a2x,
               cL_a=rest['cL'] + a1, cL_b=rest['cL'] + b1, two_wb=two_wb, cV=rest['cV'] + Ainf, two_wa=two_wa,
               C_a=rest['C'] + a0, C_b=rest['C'] + b0, C_sub=[rest['C'] + s0 for _, s0 in sub], C_d=rest['C'] + a0c, rows=out, III=(a2, a1, a0, b1, b0))
    res['C'] = res['C_b']
    res['dC'] = max(abs(res['C_a'] - res['C_b']), max(abs(x - res['C_b']) for x in res['C_sub']), abs(res['C_d'] - res['C_b']))
    if verbose:
        print(f"== {name}: qa={P['qa']}, qb={P['qb']}, |k|={P['kn']:.5f}; scaled l={l}; P_q={Pq:.5f}")
        for r, o in out.items():
            print(f"   {r:10s}: L^2 {o['c2']: .6f}  L {o['cL']: .6f}  V {o['cV']: .6f}  C {o['C']: .6f}")
        print(f"   I.III     : L^2 {a2: .6f} (analytic {a2x: .6f})  L {a1: .6f} [{b1: .6f} with L^2 fixed]  V {Ainf: .6f}  C {a0: .6f} [{b0: .6f}]")
        print(f"   Group I   : L^2 {res['c2_free']: .6f} (abar={abar:.6f})  L {res['cL_a']: .6f} [{res['cL_b']: .6f}] vs 2w(q_b)={two_wb: .6f}   V {res['cV']: .6f} vs 2w(q_a)={two_wa: .6f}")
        print(f"   C: free fit {res['C_a']:.5f} | L^2 fixed {res['C_b']:.5f} | 3-point subsets {['%.5f' % x for x in res['C_sub']]} | L^2 and L fixed {res['C_d']:.5f}  => C = {res['C']:.4f} +- {res['dC']:.4f}")
    return res

if __name__ == '__main__':
    print("cross-check of the printed one-loop vertex: BSV eq:kernelv braces (without pi^2) against 2Re DDS eq:looplip/N_c (tau=k^2), n_f=0")
    rng = np.random.default_rng(5)
    for _ in range(4):
        qa, qb = rng.normal(size=2), rng.normal(size=2)
        print(f"   BSV {X_bsv(qa, qb, False): .10f}   DDS {X_dds(qa, qb): .10f}")
    names = [n for n in PTS if n in OLD or os.path.exists(f'{S}/h2_{n}_1e-7.out') and os.path.getsize(f'{S}/h2_{n}_1e-7.out') > 0
             and (not os.path.exists(f'{S}/h2_{n}_1e-8.out') or os.path.getsize(f'{S}/h2_{n}_1e-8.out') > 0)
             and os.path.exists(f'{S}/v5_nodes_{n}.out') and len(open(f'{S}/v5_nodes_{n}.out').readlines()) == 12
             and os.path.exists(f'{S}/v4b_sgrid_{n}.out') and len(open(f'{S}/v4b_sgrid_{n}.out').readlines()) == 22
             and all(os.path.getsize(f'{S}/h2_{n}_{L}.out') > 0 for L in ('1e-5', '1e-6', '1e-7'))]
    R = {n: groupI(n) for n in names}
    print("\npoint   |k|      C(CGC)   +-      X(FL,(abar/4)X)  conv=(abar/2)ln k^2 ln(qa^2/qb^2)   C-X       D=C-X+conv   D+abar pi^2/12")
    T = {}
    for n in names:
        P = PTS[n]; X = abar / 4 * X_bsv(P['qa'], P['qb']); conv = abar / 2 * np.log(P['kn'] ** 2) * np.log((P['qa'] @ P['qa']) / (P['qb'] @ P['qb']))
        T[n] = dict(C=R[n]['C'], dC=R[n]['dC'], X=X, conv=conv, lnk2=np.log(P['kn'] ** 2), Cd=R[n]['C_d'], Ca=R[n]['C_a'])
        print(f"{n}  {P['kn']:.5f}  {R[n]['C']: .5f} {R[n]['dC']:.5f}  {X: .6f}  {conv: .6f}   {R[n]['C'] - X: .5f}  {R[n]['C'] - X + conv: .5f}  {R[n]['C'] - X + conv + abar * np.pi ** 2 / 12: .5f}")
    print("\ndouble differences over rectangles  DD[f] = f(a,b) + f(a',b') - f(a,b') - f(a',b):")
    for (a, a_) in (('a1', 'a2'), ('a1', 'a3'), ('a2', 'a3')):
        cor = [a + 'b1', a_ + 'b2', a + 'b2', a_ + 'b1']
        if not all(c in T for c in cor): continue
        sg = [1, 1, -1, -1]
        dd = lambda key: sum(s * T[c][key] for s, c in zip(sg, cor))
        err = np.sqrt(sum(T[c]['dC'] ** 2 for c in cor)); errl = sum(T[c]['dC'] for c in cor)
        print(f" ({a},{a_})x(b1,b2): DD[C] = {dd('C'): .5f} +- {err:.5f} (quadrature; linear sum {errl:.5f})   DD[X] = {dd('X'): .5f}   DD[conv] = {dd('conv'): .5f}   DD[abar ln k^2] = {abar * dd('lnk2'): .5f}")
        print(f"      DD[C - X] = {dd('C') - dd('X'): .5f}   prediction -DD[conv] = {-dd('conv'): .5f}   =>  DD[D] = {dd('C') - dd('X') + dd('conv'): .5f} +- {err:.5f}"
              f"      (variants of the I.III fit: free {dd('Ca') - dd('X') + dd('conv'): .5f}, constrained {dd('Cd') - dd('X') + dd('conv'): .5f})")
    # --- stronger form of the test.  C - (abar/4)X is invariant under rotations and rescaling, so D = F_P(q_a)+F_T(q_b) with
    # D(lam qa, lam qb) - D(qa,qb) = (abar/2) ln lam^2 ln(qa^2/qb^2) forces
    #   F_P = (abar/4) ln^2 qa^2 + c ln qa^2 + const,  F_T = -(abar/4) ln^2 qb^2 - c ln qb^2 + const,  i.e.
    #   E := C - (abar/4) X + (abar/4) ln(qa^2/qb^2) ln(k^4/(qa^2 qb^2)) = c ln(qa^2/qb^2) + e0   (two numbers for all points).
    sym = [n for n in names if abs(PTS[n]['qa'] @ PTS[n]['qa'] - PTS[n]['qb'] @ PTS[n]['qb']) < 1e-9]
    if len(sym) > 1:
        print("\npoints with |q_a| = |q_b| (all conversion terms and F_P + F_T reduce to a constant):  C - (abar/4) X must be the same number")
        for n in sym:
            P = PTS[n]; ang = np.degrees(np.arccos((P['qa'] @ P['qb']) / (P['qa'] @ P['qa'])))
            print(f"   {n}: angle(q_a,q_b) = {ang:6.1f} deg,  q^2/k^2 = {P['qa'] @ P['qa']:.4f}:  C = {T[n]['C']:.5f} +- {T[n]['dC']:.5f}   (abar/4)X = {T[n]['X']:.5f}   C - (abar/4)X = {T[n]['C'] - T[n]['X']:.5f}")
    print("\nstronger form: E = C - (abar/4)X + (abar/4) ln(qa^2/qb^2) ln(k^4/(qa^2 qb^2))  must equal  c*ln(qa^2/qb^2) + e0")
    xs, Es, dEs = [], [], []
    for n in names:
        P = PTS[n]; a2, b2 = P['qa'] @ P['qa'], P['qb'] @ P['qb']; lab = np.log(a2 / b2)
        E = T[n]['C'] - T[n]['X'] + abar / 4 * lab * np.log(P['kn'] ** 4 / (a2 * b2)); xs.append(lab); Es.append(E); dEs.append(T[n]['dC'])
    xs, Es, dEs = map(np.array, (xs, Es, dEs))
    if len(names) >= 3:
        A = np.vstack([xs, np.ones_like(xs)]).T; (c, e0), *_ = np.linalg.lstsq(A, Es, rcond=None)
        print(f"   fit: c = {c:.5f} (abar = {abar:.5f}, c/abar = {c / abar:.4f}),  e0 = {e0:.5f}  (e0 + abar pi^2/12 = {e0 + abar * np.pi ** 2 / 12:.5f})")
        for n, x, E, dE in zip(names, xs, Es, dEs):
            print(f"   {n}: ln(qa^2/qb^2) = {x: .5f}   E = {E: .5f} +- {dE:.5f}   fit residual = {E - c * x - e0: .5f}   [E + abar*x = {E + abar * x: .5f}]")
        print(f"   rms residual {np.sqrt(np.mean((Es - c * xs - e0) ** 2)):.5f} with {len(names)} points and 2 parameters;  hypothesis c = -abar: E + abar*ln(qa^2/qb^2) = {np.mean(Es + abar * xs):.5f} +- {np.std(Es + abar * xs):.5f} (spread)")
    else:
        for n, x, E in zip(names, xs, Es): print(f"   {n}: ln(qa^2/qb^2) = {x: .5f}   E = {E: .5f}")
