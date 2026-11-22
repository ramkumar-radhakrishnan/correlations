# s10_uv_sec6.py -- UV (short-distance) power counting of EVERY term of Sec. 6 ("Finite terms: after cancellation
# and symmetrization", sec:D) of Final_terms_explicit_1_author.tex, at the operator level (colour tensors C^{abc}
# with Wilson lines, before the projectile average), at every coincidence (cluster) of integration points.
# Method of paper/notes/ref3_finiteness.tex: smooth SU(3) adjoint Wilson-line field (plane waves), cluster S scaled to
# size lam = 1e-2..1e-5 (3 base configurations x NOR orientations), |integrand| ~ lam^alpha (max over samples);
# absolutely integrable at S iff alpha > -d_S, d_S = 2(|S|-1) (minus 1 when the delta(s) constraint ties y' to z).
# Kernels: lib.py of paper/notes/scripts/ref3_finite (closed forms of the p+ integrals: E = int[1/p+]_+ e^{..},
# J, Ehat, Psi), and the author-form functions of s2_kernels.py (Row I.V 2nd part, Row III.III 2nd part integrand,
# Row II.I 2nd part), copied verbatim from s2_kernels.py at run time.
# Run: python3 -I s10_uv_sec6.py [quick]   ->  s10_uv_sec6.log, s10_uv_sec6.json
import sys, os, re, itertools, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
R3 = os.path.join(HERE, '..', '..', 'paper', 'notes', 'scripts', 'ref3_finite')
sys.path.insert(0, R3)
import numpy as np
from scipy import integrate, special
import lib
E = lambda *a: np.einsum(*a, optimize='greedy')
f = lib.FABC; NC = lib.NC
# --- author-form kernels copied from s2_kernels.py (exec of the def blocks only) ---
src = open(os.path.join(HERE, 's2_kernels.py')).read()
ns = dict(np=np, integrate=integrate, special=special, gE=np.euler_gamma)
for line in ['dot = lambda a, b: a @ b', 'sq = lambda a: a @ a', 'Kv = lambda v: v/sq(v)']: exec(line, ns)
for name in ['Ecal', 'Ehat', 'Psi', 'Jf', 'quadc', 'rowV2_author', 'III3_author_integrand', 'II1_author']:
    m = re.search(r'^def %s\(.*?(?=^\S)' % name, src, re.S | re.M); exec(m.group(0), ns)
Kv = lib.Kv; KVEC = lib.KVEC
QUICK = len(sys.argv) > 1 and sys.argv[1] == 'quick'
LAMS = [1e-2, 1e-3, 1e-4, 1e-5]
NOR = 4 if QUICK else 10; NCONF = 2 if QUICK else 3

# ---------------------------------------------------------------- colour tensors (sec D, as printed)
def L_(U): return U["x'"] - U["w'"]
def C_I2(U):  return E('ga,gbc->abc', L_(U), E('hbc,gh->gbc', f, U['w']) - E('gde,db,ec->gbc', f, U['x'], U['y']))
def C_I3(u):  return lambda U: E('ga,ghd,db,hc->abc', L_(U), f, U['x'], U[u] - U['y'])
def C_I5(U):  return E('ga,gde,dc,eb->abc', L_(U), f, U['x'], U['z']) - E('ga,gde,eb,dc->abc', L_(U), f, U['y'], U["y'"])
def C_II1a(U):
    X = E('mh,mb->hb', U['z'], U['y'])
    return (E('ga,kbc,gk->abc', U["w'"], f, U['w']) - E('ga,gkl,kb,lc->abc', U["w'"], f, U['z'], U['x'])
            - E('ga,hb,khc,gk->abc', U["x'"], X, f, U['w']) + E('ga,hb,gkl,kh,lc->abc', U["x'"], X, f, U['z'], U['x']))
def C_II1b(U):
    X = E('mh,mc->hc', U['z'], U['y'])
    return (E('gb,kca,gk->abc', U['w'], f, U["w'"]) - E('gb,gkl,kc,la->abc', U['w'], f, U['z'], U['x'])
            - E('gb,hc,kha,gk->abc', U["x'"], X, f, U["w'"]) + E('gb,hc,gkl,kh,la->abc', U["x'"], X, f, U['z'], U['x']))
def C_III2(U):
    return (E('ma,gb,mh,khc,gk->abc', U['y'], U["x'"], U['z'], f, U['w']) - E('ma,gb,mh,gkl,kh,lc->abc', U['y'], U["x'"], U['z'], f, U['z'], U['x'])
            - E('gb,kac,gk->abc', U["x'"], f, U['w']) + E('gb,gkl,ka,lc->abc', U["x'"], f, U['z'], U['x']))
def C_III3(v): return lambda U: E('mnk,ma,nckb->abc', f, U[v] - U["x'"], E('nc,kb->nckb', U['z'], U['w']) - E('kb,nc->nckb', U['y'], U['x']))
def C_III5(v): return lambda U: E('mkl,ma,kb,lc->abc', f, U[v] - U["x'"], U['x'], U['y'] - U['w'])
def C_III6(v): return lambda U: E('mnk,ma,kb,nc->abc', f, U[v] - U["x'"], U['y'], U['x'] - U['z'])
CONTRACT = {'ab': 'abc,abr->rc', 'ac': 'abc,acr->rb', 'bc': 'abc,bcr->ar'}   # f^{abr}, f^{acr}, f^{bcr}

# ---------------------------------------------------------------- kernels (sec D, as printed; c-numbers)
def Kf(a, b, z): return Kv(a - z) @ Kv(b - z)
def kI3z(P): return lib.K_III_z(P)                       # 24-term kernel + (l-E) + E  (= l K(y,w,z) P.K(x-w))
def kI3w(P): return lib.K_III_w(P)                       # -E(k.(z-w)) K(x,w,z) K(x'-w').K(y-w)
def kV1(P):  return lib.K_V1(P)                          # J, E, e^{-is}E(-s) lines
def kT(P):   return lib.K_T(P, False)
def kTs(P):  return lib.K_T(P, True)
def g3(P):
    x, y, z, w, xp, yp = (P[k] for k in ('x', 'y', 'z', 'w', "x'", "y'"))
    A1 = (Kv(yp - z) @ Kv(xp - yp)) * (Kv(y - w) @ Kv(x - z)); B1 = (Kv(yp - z) @ Kv(x - z)) * (Kv(xp - yp) @ Kv(y - w))
    return A1, B1
def kG3reg(P):   # Ehat (A1-B1) - Psi B1, tensor C(y')
    A1, B1 = g3(P); s = lib.s3(P); return lib.Ehat(s) * (A1 - B1) - lib.Psi_fun(s) * B1

# ---------------------------------------------------------------- terms
class T:
    def __init__(s, row, name, pts, lines, dist=None, contact=None, fourier=()):
        s.row, s.name, s.pts, s.lines, s.dist, s.contact = row, name, pts, lines, dist, contact
        s.drop, s.keep = contact[1] if contact else (None, None)
    def full(s, P):
        Q = dict(P)
        if s.drop: Q[s.drop] = Q[s.keep]
        return Q
    def value(s, Uf, P):
        Q = s.full(P); U = {k: Uf(Q[k]) for k in Q}
        tot = 0; k0 = 0
        for kern, col in s.lines:
            k = kern(Q); C = col(U)
            if s.contact: C = E(CONTRACT[s.contact[0]], C, f)
            tot = tot + k * C; k0 = k0 + k
        if s.dist == 'pv':
            sv = lib.s3(Q); tot = tot / sv; k0 = k0 / sv
        return np.sqrt(np.sum(np.abs(tot) ** 2)), abs(k0)
    def label(s):
        c = '' if not s.contact else ' [contact f^{%sr}: %s->%s]' % (s.contact[0], s.drop, s.keep)
        return '%s %s%s%s' % (s.row, s.name, ' (%s)' % s.dist if s.dist else '', c)

G1 = ['x', 'y', 'z', 'w', "x'", "w'"]; G15 = ['x', 'z', 'y', "y'", "x'", "w'"]; G3 = ['w', "x'", "y'", 'z', 'y', 'x']
TERMS = []
def add3(row, name, pts, lines, contacts, dist=None):
    TERMS.append(T(row, name + ' 3rho', pts, lines, dist))
    for c in contacts: TERMS.append(T(row, name + ' contact', pts, lines, dist, contact=c))
# charge order (a,b,c) and contact sets as printed in sec D
add3('I.III', 'Kz.Cz', G1, [(kI3z, C_I3('z'))], [('ab', ("x'", 'x')), ('bc', ('y', 'x')), ('ac', ("x'", 'y'))])
add3('I.III', 'Kw.Cw', G1, [(kI3w, C_I3('w'))], [('ab', ("x'", 'x')), ('bc', ('y', 'x')), ('ac', ("x'", 'y'))])
add3('I.V', 'KV1.C', G15, [(kV1, C_I5)], [('ab', ("x'", 'y')), ('ac', ("x'", "y'"))])
add3('II.I', 'KT(w,wp).C1a', G1, [(kT, C_II1a)], [('ac', ("x'", 'x')), ('bc', ('y', 'x'))])
add3('II.I', 'KT(wp,w).C1b', G1, [(kTs, C_II1b)], [('ab', ("x'", 'x')), ('ac', ('y', 'x'))])
add3('III.II', 'KT.C', G1, [(kT, C_III2)], [('ab', ("x'", 'y')), ('ac', ('y', 'x')), ('bc', ("x'", 'x'))])
G3C = {'III.III': (C_III3, [('ab', ("x'", 'y')), ('ac', ("x'", 'x'))]),
       'III.V': (C_III5, [('ab', ("x'", 'x')), ('bc', ('y', 'x')), ('ac', ("x'", 'y'))]),
       'III.VI': (C_III6, [('ab', ("x'", 'y')), ('bc', ('y', 'x')), ('ac', ("x'", 'x'))])}
for row, (Cf, cts) in G3C.items():
    add3(row, '(Ehat,Psi).C(yp)', G3, [(kG3reg, Cf("y'"))], cts)
    for dist in ('delta', 'pv'):
        add3(row, 'inst[gyp.C(yp)+gz.C(z)]', G3, [(lib.g_delta_yp, Cf("y'")), (lib.g_delta_z, Cf('z'))], cts, dist=dist)

# ---- original (genuine) two-rho terms; 2-index colour tensors (p: first charge, q: second) ----
def c_I2o(U): return E('ap,aq->pq', U["x'"] - U["w'"], U['w'] - U['x'])
def k_I2o(P): return Kv(P["x'"] - P["w'"]) @ Kv(P['x'] - P['w'])
TERMS.append(T('I.II', '(-67/9+pi^2/3) 2rho', ['x', 'w', "x'", "w'"], [(k_I2o, c_I2o)]))
def c_V2a(U): return E('ap,aq->pq', U["x'"] - U["w'"], -E('adc,fbe,db,ce->af', f, f, U['x'], U['z']) + NC * U['y'])
def c_V2b(U): return E('ap,aq->pq', U["x'"] - U["w'"], U['y'] - U['x'])
def k_V2a(P):
    s = KVEC @ (P['x'] - P['z']); return ns['rowV2_author'](s, P['x'], P['y'], P['z'], P["x'"], P["w'"])
def k_V2b(P):
    x, y, z, xp, wp = P['x'], P['y'], P['z'], P["x'"], P["w'"]; s = KVEC @ (x - z)
    return NC * 11 / 6 * np.exp(-1j * s) * (Kv(xp - wp) @ Kv(y - x)) * ((y - x) @ Kv(y - z)) / ((x - z) @ (x - z))
TERMS.append(T('I.V', '2nd part 2rho (reg+soft+11/6)', ['x', 'y', 'z', "x'", "w'"], [(k_V2a, c_V2a), (k_V2b, c_V2b)]))
TERMS.append(T('I.V', '2nd part 2rho without 11/6', ['x', 'y', 'z', "x'", "w'"], [(k_V2a, c_V2a)]))
def c_II2(U):
    t1 = E('xbp,ybq,ax,ay->pq', f, f, U["w'"], U['w'])
    t2 = -E('xbp,acd,ax,cb,dq->pq', f, f, U["w'"], U['z'], U['x'])
    t3 = -E('axd,cbq,ac,xb,dp->pq', f, f, U['w'], U['z'], U["x'"])
    return t1 + t2 + t3 + NC * E('dp,dq->pq', U["x'"], U['x'])
def k_II2(P): return ns['II1_author'](P['x'], P["x'"], P['z'], P['w'], P["w'"])[0]
TERMS.append(T('II.I', '2nd part 2rho (N2a..N2e)', ['x', "x'", 'z', 'w', "w'"], [(k_II2, c_II2)]))
def I_IV(s):   # int_0^inf dt t/(1+t)^2 e^{-is(1+t)}
    return np.exp(-1j * s) * ((1 + 1j * s) * lib.Ehat(s) - 1)
def k_IV(P):
    x, y, xp, yp = P['x'], P['y'], P["x'"], P["y'"]; s = KVEC @ (y - yp)
    return I_IV(s) * (Kv(x - y) @ Kv(xp - yp) - 2 * (Kv(x - y) @ (yp - y)) * (Kv(xp - yp) @ Kv(yp - y)))
def c_IV(U): return E('cp,cq->pq', U["y'"] - U["x'"], U['y'] - U['x'])
TERMS.append(T('II.IV', 'tensor term 2rho', ['x', 'y', "x'", "y'"], [(k_IV, c_IV)]))
def k_IVr(P): return Kf(P['x'], P["x'"], P['z'])
def c_IVr(U): return E('cp,cq->pq', U['z'] - U["x'"], U['z'] - U['x'])
TERMS.append(T('II.IV', 'remnant 2rho', ['x', "x'", 'z'], [(k_IVr, c_IVr)]))
# Group III Row III, second part (genuine 2rho); colour C(v)^{pq} = L^{abp}(v) R^{abq}
def c_33(v):
    def c(U):
        Lt = E('cda,db,cp->abp', f, U['z'], U[v] - U["x'"])
        Rt = E('cbq,ac->abq', f, U['w']) - E('acd,cb,dq->abq', f, U['z'], U['x'])
        return E('abp,abq->pq', Lt, Rt)
    return c
def R33(P):    # int_0^inf dt (author's integrand: regular + soft lines), oscillatory tail with weights
    x, xp, yp, z, w = P['x'], P["x'"], P["y'"], P['z'], P['w']; s = KVEC @ (z - yp)
    g = lambda t: ns['III3_author_integrand'](t, s, x, xp, yp, z, w)
    T0 = 40.0
    I = ns['quadc'](g, 0, T0)
    a = x - w; b = x - z; A, B = a @ a, b @ b
    # tail: g(t) = e^{-ist} G(t) + I_/(A t(1+t)); G evaluated from g
    I_ = (a @ Kv(xp - yp)) * (Kv(z - w) @ Kv(yp - z)) + 0.5 * (Kv(x - z) @ Kv(yp - z)) * (a @ Kv(xp - yp))
    G = lambda t: (g(t) - I_ / (A * t * (1 + t))) * np.exp(1j * s * t)
    wv = abs(s)
    if wv < 1e-12: return I + ns['quadc'](g, T0, np.inf)
    Gr = lambda t: G(t).real; Gi = lambda t: G(t).imag
    c_r = integrate.quad(Gr, T0, np.inf, weight='cos', wvar=wv)[0]; s_r = integrate.quad(Gr, T0, np.inf, weight='sin', wvar=wv)[0]
    c_i = integrate.quad(Gi, T0, np.inf, weight='cos', wvar=wv)[0]; s_i = integrate.quad(Gi, T0, np.inf, weight='sin', wvar=wv)[0]
    sg = np.sign(s)   # e^{-ist} = cos(|s|t) - i sg sin(|s|t)
    tail = (c_r + 1j * c_i) - 1j * sg * (s_r + 1j * s_i) + I_ / A * np.log((1 + T0) / T0)
    return I + tail
def X33(P, v):
    x, z, w, xp, yp = P['x'], P['z'], P['w'], P["x'"], P["y'"]; q = P[v]
    return ((x - z) @ Kv(xp - q)) * (Kv(z - w) @ Kv(yp - z)) - 0.5 * ((x - z) @ Kv(xp - q)) * (Kv(x - w) @ Kv(yp - z))
def k33_reg_y(P):
    a = P['x'] - P['w']; b = P['x'] - P['z']; A, B = a @ a, b @ b; s = KVEC @ (P['z'] - P["y'"])
    return R33(P) + A / B ** 2 * lib.Ehat(s * A / B) * X33(P, "y'")
def k33_reg_z(P):
    a = P['x'] - P['w']; b = P['x'] - P['z']; A, B = a @ a, b @ b; s = KVEC @ (P['z'] - P["y'"])
    return -A / B ** 2 * lib.Ehat(s * A / B) * X33(P, 'z')
def k33_d_y(P): b = P['x'] - P['z']; return -X33(P, "y'") / (b @ b)
def k33_d_z(P): b = P['x'] - P['z']; return +X33(P, 'z') / (b @ b)
P33 = ['w', "x'", "y'", 'z', 'x']
TERMS.append(T('III.III', '2nd part 2rho regular+soft (+Ehat(sA/B) part of inst.)', P33, [(k33_reg_y, c_33("y'")), (k33_reg_z, c_33('z'))]))
for dist in ('delta', 'pv'):
    TERMS.append(T('III.III', '2nd part 2rho instantaneous 1/B part', P33, [(k33_d_y, c_33("y'")), (k33_d_z, c_33('z'))], dist=dist))
def k33_log(P):
    x, z, w, xp, wp = P['x'], P['z'], P['w'], P["x'"], P["w'"]; A = (x - w) @ (x - w); B = (x - z) @ (x - z)
    return A / (B - A) * np.log(B / A) * (Kf(x, xp, z) * (Kv(z - w) @ Kv(wp - z)) - 0.5 * Kf(x, xp, z) * (Kv(x - w) @ Kv(wp - z)))
def c_33log(U): return c_33('z')(U)
TERMS.append(T('III.III', '2nd part 2rho last (log) term', ['w', "x'", "w'", 'z', 'x'], [(k33_log, c_33log)]))

# ---------------------------------------------------------------- scan
def make_P(t, base, S, lam, dirs, ubase, ucl):
    P = {p: base[p].copy() for p in t.pts if p != t.drop}
    c = base[t.keep] if (t.keep in S) else base[S[0]]
    for q in S: P[q] = c + lam * dirs[q]
    P[t.keep if t.keep in S else S[0]] = c.copy()
    if t.dist == 'delta':
        if "y'" in S and 'z' in S: P["y'"] = P['z'] + lam * ucl * lib.EPERP
        elif "y'" in S: P['z'] = np.array([P["y'"][0], P["y'"][1] - ubase])
        else: P["y'"] = P['z'] + ubase * lib.EPERP
    return P
def clusters(t):
    pts = [p for p in t.pts if p != t.drop]
    out = []
    for n in range(2, len(pts) + 1):
        for comb in itertools.combinations(pts, n):
            S = list(comb)
            if t.keep in S: S = [t.keep] + [p for p in S if p != t.keep]
            out.append(S)
    return out
def dimS(t, S):
    d = 2 * (len(S) - 1)
    if t.dist == 'delta' and "y'" in S and 'z' in S: d -= 1
    return d
def scan(t, field, seed):
    rng = np.random.default_rng(seed); Uf = lib.Colour(field).u
    bases = [{p: rng.uniform(-1.5, 1.5, size=2) for p in t.pts} for _ in range(NCONF)]
    res = []
    for S in clusters(t):
        m1 = np.zeros(len(LAMS)); m0 = np.zeros(len(LAMS))
        for base in bases:
            ubase = rng.choice([-1, 1]) * rng.uniform(0.5, 1.5)
            for o in range(NOR):
                dirs = {q: rng.uniform(0.5, 1.5) * np.array([np.cos(a), np.sin(a)]) for q, a in zip(S, rng.uniform(0, 2 * np.pi, len(S)))}
                ucl = rng.choice([-1, 1]) * rng.uniform(0.5, 1.5)
                for i, lam in enumerate(LAMS):
                    v1, v0 = t.value(Uf, make_P(t, base, S, lam, dirs, ubase, ucl))
                    m1[i] = max(m1[i], v1); m0[i] = max(m0[i], v0)
        sl = lambda m, i: np.log10(m[i] / m[i - 1]) / np.log10(LAMS[i] / LAMS[i - 1])
        res.append(dict(S=S, d=dimS(t, S), a=sl(m1, 3), a_prev=sl(m1, 2), a0=sl(m0, 3)))
    return res
def circle_mean(t, field, S2, nang=64, seed=3):
    """|angular mean| of the term on circles around S2[0] for the pair cluster S2, times lam^2 (-> 0 iff the mean
    is weaker than 1/r^2) and the max of lam^2|term|."""
    rng = np.random.default_rng(seed); Uf = lib.Colour(field).u; out = []
    base = {p: rng.uniform(-1.5, 1.5, size=2) for p in t.pts}
    for lam in [1e-2, 1e-3, 1e-4]:
        acc = 0; mx = 0
        for a in np.linspace(0, 2 * np.pi, nang, endpoint=False):
            P = make_P(t, base, S2, lam, {S2[0]: np.zeros(2), S2[1]: np.array([np.cos(a), np.sin(a)])}, 0.9, 1.0)
            Q = t.full(P); U = {k: Uf(Q[k]) for k in Q}
            tot = sum(k(Q) * c(U) for k, c in t.lines); acc = acc + tot; mx = max(mx, np.sqrt(np.sum(np.abs(tot) ** 2)))
        out.append((lam, lam ** 2 * np.sqrt(np.sum(np.abs(acc / nang) ** 2)), lam ** 2 * mx))
    return out

if __name__ == '__main__':
    field = lib.SmoothField(seed=7)
    which = [a for a in sys.argv[1:] if a != 'quick']
    out = []; log = open(os.path.join(HERE, 's10_uv_sec6%s.log' % ('_quick' if QUICK else '')), 'w')
    def P_(s): print(s, flush=True); log.write(s + '\n'); log.flush()
    P_('# s10_uv_sec6: |term| ~ lam^alpha at cluster S (d_S = integration dimension); OK: alpha > -d_S + 0.2; MARG: within 0.2; FAIL')
    t0 = time.time()
    for it, t in enumerate(TERMS):
        if which and not any(w in t.label() for w in which): continue
        r = scan(t, field, it)
        worst = min(r, key=lambda e: e['a'] + e['d']); worst0 = min(r, key=lambda e: e['a0'] + e['d'])
        bad = [e for e in r if e['a'] <= -e['d'] + 0.2]
        P_('%-70s clusters=%3d  worst(with colour): S={%s} d=%d alpha=%.2f (prev %.2f) | kernel alone: S={%s} d=%d alpha0=%.2f | %s  [%.0fs]' % (
            t.label(), len(r), ','.join(worst['S']), worst['d'], worst['a'], worst['a_prev'], ','.join(worst0['S']), worst0['d'], worst0['a0'],
            'OK' if not bad else 'NOT-ABS: ' + '; '.join('{%s} d=%d a=%.2f' % (','.join(e['S']), e['d'], e['a']) for e in bad), time.time() - t0))
        out.append(dict(term=t.label(), res=r))
    json.dump(out, open(os.path.join(HERE, 's10_uv_sec6%s.json' % ('_quick' if QUICK else '')), 'w'), indent=0, default=float)
