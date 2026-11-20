# Verification of the row-by-row symmetrised three-rho equations ((3rho) equations) as PRINTED in a TeX file,
# against the direct contraction of the colour structures of the rows with d^{abc}, including the row factors A
# read off the row equations (app_finite.tex / sec_finite.tex), and of their sum against eq:sodd.
# usage: python3 check_rows.py <sec_finite.tex> <app_finite.tex> [style]   style = paper | fixed
import sys, os, re, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
DABC = os.path.dirname(HERE); SING = os.path.join(os.path.dirname(DABC), 'singlet')
sys.path.insert(0, DABC); sys.path.insert(0, SING)
import numpy as np, sympy as sp
import fasteinsum
from rows import CUB, POINTS
from engine import parse, eval_terms
from sun import SUN
import invariants as I

SEC, APP = sys.argv[1], sys.argv[2]
STYLE = sys.argv[3] if len(sys.argv) > 3 else 'paper'
sec = open(SEC).read(); app = open(APP).read()
EQS = open(sys.argv[4]).read() if STYLE == 'fixed' else sec   # file with the printed (3rho) equations

def eqbody(src, label):
    i = src.index(r'\label{eq:%s}' % label); j = src.index(r'\end{align}', i)
    k = src.rindex(r'\begin{align}', 0, i)
    return src[k:j]

# ---------- row factors A (units of C_f), read from the row equations ----------
g, pi_, I_ = sp.symbols('g'), sp.pi, sp.I
def tex2sym(t):
    t = re.sub(r"f\^\{[^}]*\}", "", t).replace(' ', '')
    t = t.replace(r'\pi', 'pi').replace('g^4', 'g**4').replace('^', '**')
    t = re.sub(r'(\d)(pi|g|i)', r'\1*\2', t); t = re.sub(r'^i', 'I*', t); t = t.replace('ig', 'I*g')
    return sp.sympify(t, locals={'pi': sp.pi, 'I': sp.I, 'g': g})
def appfactor(label):
    b = eqbody(APP, label) if False else eqbody(app, label)
    m = re.search(r"&=([+-]?)\\frac\{1\}\{\(2\\pi\)\^3\}\\,\\frac\{((?:[^{}]|\{[^{}]*\})*)\}\{([^}]*)\}\\,\\frac\{1\}\{k\^\+\}", b)
    assert m, label
    sign = -1 if m.group(1) == '-' else 1
    val = sign * tex2sym(m.group(2)) / tex2sym(m.group(3)) / (g**4 / (16 * sp.pi**5))
    return complex(sp.nsimplify(sp.simplify(val))), m.group(0)
A = {}
A['I.II.N2b'], s1 = appfactor('app_IIb')
assert r'\mathcal L^a\Big[f^{hbc}\,U^{ah}(w)\,\{\rho^b(x),\rho^c(y)\}' in app
# I.III: Eq. GIIII has 2i C_f for the brace with M^{ij} and ell (= K^III_z with C_D); the C_A member is in Eq. Sblock:
assert r"\frac{d^3N^{\rm I.III}}{d^3k}&=2i\,\mathcal C_{\rm f}\,f^{abd}" in sec
assert r"-i\,P\cdot\Kv(y-w)\,\big[2C_A+C_B\big]\,\mathfrak F\big[K(x,\cdot,z)\big]" in sec
assert r"\mathfrak K^{\rm III}_w&\equiv-\mathcal E\big(\v k\cdot(\v z-\v w)\big)\,K(x,w,z)\,P\cdot\Kv(y-w)" in sec
A['I.III(z)'] = 2j
A['I.III(w)'] = (-1j) * 2 * (-1)        # -i*2 (coefficient of C_A in Sblock) times (-1) from the sign in K^III_w
A['I.V.N1'], _ = appfactor('app_V1')    # coefficient of the J-term; K^V_1 starts with +J(s) P.K(z-x) K(y-z).K(y'-x)
assert r"\mathfrak K^{\rm V}_1&\equiv\mathcal J(s)\,P\cdot\Kv(z-x)\,\Kv(y-z)\cdot\Kv(y'-x)" in sec
assert r"\mathcal J(s)\,P\cdot\Kv(z-x)\;\Kv(y-z)\cdot\Kv(y'-x)" in eqbody(app, 'app_V1')
A['II.I.N1a'], _ = appfactor('app_III1a')
b = app[app.index(r'\label{eq:app_III1b}'):app.index(r'\end{align}', app.index(r'\label{eq:app_III1b}'))]
m = re.search(r"&=([+-])\\frac\{1\}\{\(2\\pi\)\^3\}\\,\\frac\{([^}]*)\}\{([^}]*)\}", b)
A['II.I.N1b'] = complex(sp.nsimplify((1 if m.group(1) == '+' else -1) * tex2sym(m.group(2)) / tex2sym(m.group(3)) / (g**4 / (16 * sp.pi**5))))
A['III.II'], _ = appfactor('app_IIIII')
# Group III: brace {Ehat[..]L(y') - [pi delta - iP/s]..[K(x'-y')L(y') - K(x'-z)L(z)]}: L(y') coefficient = P_app*(K_Ehat - K_delta^(y'))
#  = -P_app*(-K_Ehat + K_delta^(y')); L(z): P_app*(+[..]K(x'-z)) = -P_app*K_delta^(z) (K_delta^(z) carries the minus sign).
# The Psi member (block S, pairs A, B, D with Eqs. colid) has the same coefficient -P_app (checked in the text of the script header).
for row, lab in (('III.III.1', 'app_III31'), ('III.V', 'app_III5')):
    P, _ = appfactor(lab)
    A[row + "(y')"] = -P; A[row + "(z)"] = -P
assert "Row~III.VI is Eq.~\\eqref{eq:app_III5} with" in app and "with the same prefactor" in app
A["III.VI(y')"] = A["III.V(y')"]; A["III.VI(z)"] = A["III.V(z)"]
assert r"\mathfrak K_\delta^{(z)}&\equiv-\Big[\pi\delta(s)-i\mathcal P\frac1s\Big]" in sec

# ---------- parse the (3rho) equations ----------
KER = {r"\mathfrak K^{\rm II}_{2b}": 'K2b', r"\mathfrak K^{\rm III}_z": 'Kz3', r"\mathfrak K^{\rm III}_w": 'Kw3',
       r"\mathfrak K^{\rm V}_1": 'KV1', r"\mathfrak K_{\mathcal T}(w,w')": 'KT', r"\mathfrak K_{\mathcal T}(w',w)": 'KTsw',
       r"\Big[-\mathfrak K_{\widehat{\mathcal E}}+\mathfrak K_\Psi+\mathfrak K_\delta^{(y')}\Big]": 'Ky', r"\mathfrak K_\delta^{(z)}": 'Kzd'}
COEF = {'+2i': 2j, '2i': 2j, '-2i': -2j, '-i': -1j, '+i': 1j, r'\frac{1}{2\pi}': 1 / (2 * np.pi), r'+\frac{1}{2\pi}': 1 / (2 * np.pi)}
def pt(s): return s.strip()
def parse_bracket(t):
    terms = []
    t = t.replace(r'\mf{G}', r'\mf G').replace(r'\mf{Q}_{d}', r'\mf Q_d')
    for m in re.finditer(r"([+-]?)\s*(\d*)\\mf (G|Q_d)\(([^)]*)\)", t):
        c = (-1 if m.group(1) == '-' else 1) * (int(m.group(2)) if m.group(2) else 1)
        args = tuple(pt(a) for a in re.split(r'[,;]', m.group(4)))
        terms.append((c, m.group(3), args))
    chk = re.sub(r"([+-]?)\s*(\d*)\\mf (G|Q_d)\(([^)]*)\)", "", t).strip()
    assert chk == '', ('unparsed bracket rest', chk)
    return terms
def eq_terms(label):
    b = eqbody(EQS, label)
    out = []
    pat = re.compile(r"&\\quad([^\n]*?)\\notag\\\\\n&\\quad\\times\\Big\[(.*?)\\Big\]\\,o\(([^)]*)\)")
    for m in pat.finditer(b):
        head = m.group(1).strip()
        for k in sorted(KER, key=len, reverse=True):
            if head.endswith(k):
                coef = head[:-len(k)].rstrip(r'\,').strip()
                break
        else: raise ValueError(head)
        out.append(dict(coef=COEF[coef], ker=KER[k], br=parse_bracket(m.group(2)), o=tuple(pt(a) for a in m.group(3).split(','))))
    return b, out
def evalbr(br, G, U, sub=None):
    v = 0
    for c, kind, a in br:
        a = tuple(sub.get(q, q) for q in a) if sub else a
        if kind == 'G': v += c * I.Gfrak(G, U[a[0]], U[a[1]], U[a[2]])
        else: v += c * I.Qd(G, U[a[0]], U[a[1]], U[a[2]], U[a[3]]).real
    return v

ROWEQ_P = {'oIIb': ('sIIb', {'K2b': 'I.II.N2b'}), 'oIIII': ('sIIII', {'Kz3': 'I.III(z)', 'Kw3': 'I.III(w)'}),
         'oIV1': ('sIV1', {'KV1': 'I.V.N1'}), 'oIII1': ('sIII1', {'KT': 'II.I.N1a', 'KTsw': 'II.I.N1b'}),
         'oIIIII': ('sIIIII', {'KT': 'III.II'}), 'oIIIIII1': ('sIIIIII1', {'Ky': "III.III.1(y')", 'Kzd': 'III.III.1(z)'}),
         'oIIIV': ('sIIIV', {'Ky': "III.V(y')", 'Kzd': 'III.V(z)'}), 'oIIIVI': ('sIIIVI', {'Ky': "III.VI(y')", 'Kzd': 'III.VI(z)'})}
ROWEQ_F = {'s-o-I-II-2b': ('s-I-II-2b', {'K2b': 'I.II.N2b'}), 's-o-I-III': ('s-I-III', {'Kz3': 'I.III(z)', 'Kw3': 'I.III(w)'}),
         's-o-I-V-1': ('s-I-V-1', {'KV1': 'I.V.N1'}), 's-o-II-I-1': ('s-II-I-1', {'KT': 'II.I.N1a', 'KTsw': 'II.I.N1b'}),
         's-o-III-II': ('s-III-II', {'KT': 'III.II'}), 's-o-III-III-1': ('s-III-III-1', {'Ky': "III.III.1(y')", 'Kzd': 'III.III.1(z)'}),
         's-o-III-V': ('s-III-V', {'Ky': "III.V(y')", 'Kzd': 'III.V(z)'}), 's-o-III-VI': ('s-III-VI', {'Ky': "III.VI(y')", 'Kzd': 'III.VI(z)'})}
ROWEQ = ROWEQ_P if STYLE == 'paper' else ROWEQ_F
if STYLE == 'fixed':
    KER.update({r"\mathfrak{K}^{\rm II}_{2b}": 'K2b', r"\mathfrak{K}^{\rm III}_{z}": 'Kz3', r"\mathfrak{K}^{\rm III}_{w}": 'Kw3',
       r"\mathfrak{K}^{\rm V}_{1}": 'KV1', r"\mathfrak{K}_{\mathcal{T}}(w,w')": 'KT', r"\mathfrak{K}_{\mathcal{T}}(w',w)": 'KTsw',
       r"\Big[-\mathfrak{K}_{\widehat{\mathcal{E}}}+\mathfrak{K}_{\Psi}+\mathfrak{K}_{\delta}^{(y')}\Big]": 'Ky', r"\mathfrak{K}_{\delta}^{(z)}": 'Kzd'})
PREF = r"\frac{\mathcal C_{\rm f}\,N_c}{(N_c^2-4)(N_c^2-1)}" if STYLE == 'paper' else r"\frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{16\pi^{5}}\,\frac{1}{k^{+}}\,\frac{N_{c}}{(N_{c}^{2}-4)(N_{c}^{2}-1)}"
parsed = {}
seen = set()
for lab, (slab, kmap) in ROWEQ.items():
    b, terms = eq_terms(lab)
    sb = eqbody(EQS, slab)
    assert PREF in b, lab
    ints = lambda t: re.findall(r"\\int_\{[^}]*\}(?:e\^\{-ik\\cdot\([^)]*\)\})?", t)
    ph = lambda t: re.findall(r"e\^\{-ik\\cdot\([^)]*\)\}", t)
    var = lambda t: set(','.join(re.findall(r"\\int_\{([^}]*)\}", t)).split(','))
    assert ph(b) == ph(sb) and len(ph(b)) == 1, (lab, ph(b), ph(sb))
    assert var(b) >= var(sb), (lab, var(b), var(sb))
    assert ('{\\rm c.c.}' in b) == ('{\\rm c.c.}' in sb), lab
    assert sorted(t['ker'] for t in terms) == sorted(kmap), (lab, [t['ker'] for t in terms])
    for t in terms:
        name = kmap[t['ker']]; seen.add(name)
        assert t['o'] == CUB[name][2], (lab, name, t['o'], CUB[name][2])
        parsed[name] = t
    print('%-9s parsed: %s; integrals/phase and c.c. identical to Eq. %s' % (lab, ', '.join('%s:%s' % (t['ker'], kmap[t['ker']]) for t in terms), slab))
assert seen == set(CUB), set(CUB) - seen
worst = 0
for N in (3, 4):
    G = SUN(N)
    for seed in range(4):
        rng = np.random.default_rng(1000 * N + seed); U = {p: G.randU(rng) for p in POINTS}
        for name, (tex, typ, pts, phi) in CUB.items():
            ev = eval_terms(G, parse(tex), U)
            direct = sum(np.einsum('abc,abc->', T, G.d) for T in ev.values())   # sum over orderings of C_seq d
            expect = A[name] * direct
            got = parsed[name]['coef'] * evalbr(parsed[name]['br'], G, U)
            sc = max(1.0, abs(expect)); r = abs(got - expect) / sc; worst = max(worst, r)
            if seed == 0 and N == 3: print('  SU(3) %-15s A=%-8s t*A*sum C.d = %+.10f%+.10fi  printed %+.10f%+.10fi' % (name, A[name], expect.real, expect.imag, got.real, got.imag))
print('ROW EQUATIONS vs direct contraction (incl. row factor A): worst rel. deviation %.1e' % worst)
RELROW = worst

# ---------- sum vs eq:sodd ----------
if STYLE == 'paper':
    db = eqbody(sec, 'Deltas')
    DEL = {}
    for m in re.finditer(r"\\Delta_(\{\\rm [^}]*\}|\{2b\})(\(u\))?&\\equiv(.*?)\\,,", db):
        DEL[m.group(1).replace(r'{\rm ', '').strip('{}')] = parse_bracket(m.group(3))
    assert sorted(DEL) == sorted(['2b', 'I.III', 'I.V', 'II.I', 'III.II', 'III.III', 'III.V', 'III.VI']), DEL.keys()
    sb = eqbody(sec, 'sodd')
    assert PREF.replace(r'\,N_c}{(N_c^2-4)(N_c^2-1)}', '') or True
    assert r"\frac{\mathcal C_{\rm f}}{N_c^2-1}\,\frac{N_c}{N_c^2-4}" in sb
    SOD = [  # (printed fragment, kernel-term list [(struct, coef, Delta key, u)])
     (r"\frac{1}{2\pi}\,\mathfrak K^{\rm II}_{2b}\,\Delta_{2b}", [('I.II.N2b', 1 / (2 * np.pi), '2b', None)]),
     (r"2i\Big[\mathfrak K^{\rm III}_z\,\Delta_{\rm I.III}(z)+\mathfrak K^{\rm III}_w\,\Delta_{\rm I.III}(w)\Big]", [('I.III(z)', 2j, 'I.III', 'z'), ('I.III(w)', 2j, 'I.III', 'w')]),
     (r"-\int_{x,z,y,y',x',w'}e^{-ik\cdot(z-w')}\,i\,\mathfrak K^{\rm V}_1\,\Delta_{\rm I.V}", [('I.V.N1', -1j, 'I.V', None)]),
     (r"-\int_{x,y,x',z,w,w'}e^{-ik\cdot(w-w')}\,2i\,\mathfrak K_{\mathcal T}(w,w')\,\Delta_{\rm III.II}", [('III.II', -2j, 'III.II', None)]),
     (r"-\int_{w,x',y',z,y,x}e^{-ik\cdot(w-y')}\,2i", [(r + v, -2j, 'III.' + R, u) for R, r in (('III', 'III.III.1'), ('V', 'III.V'), ('VI', 'III.VI')) for v, u in (("(y')", "y'"), ('(z)', 'z'))]),
    ]
    assert r"\sum_{R={\rm III,V,VI}}\Big[\mathfrak K_{y'}\,\Delta_{{\rm III}.R}(y')+\mathfrak K^{(z)}_\delta\,\Delta_{{\rm III}.R}(z)\Big]" in sb
    assert r"\mathfrak K_{y'}\equiv-\mathfrak K_{\widehat{\mathcal E}}+\mathfrak K_\Psi+\mathfrak K^{(y')}_\delta" in sec
    worst = 0
    for N in (3, 4):
        G = SUN(N)
        for seed in range(3):
            rng = np.random.default_rng(77 * N + seed); U = {p: G.randU(rng) for p in POINTS}
            for frag, lst in SOD:
                assert frag in sb, frag
                for name, c, dk, u in lst:
                    key = {'2b': '2b'}.get(dk, r'\rm ' + dk) if dk != '2b' else '2b'
                    br = DEL[dk]
                    v1 = c * evalbr(br, G, U, {'u': u} if u else None)
                    v2 = parsed[name]['coef'] * evalbr(parsed[name]['br'], G, U)
                    r = abs(v1 - v2) / max(1, abs(v2)); worst = max(worst, r)
    print('SUM of the row equations vs eq:sodd (per kernel and structure): worst rel. deviation %.1e' % worst)
    RELSUM = worst

# ---------- II.I: two lines with exp(-ik(w-w')) vs the sin form (discretised w, w' integral) ----------
if True:
    b, terms = eq_terms('oIII1' if STYLE == 'paper' else 's-o-II-I-1'); bs = eqbody(EQS, 'oIII1sin' if STYLE == 'paper' else 's-o-II-I-1sin')
    KTn = r"\mathfrak K_{\mathcal T}(w,w')" if STYLE == 'paper' else r"\mathfrak{K}_{\mathcal{T}}(w,w')"
    SIN = r"4\sin\big(\v k\cdot(\v w-\v w')\big)" if STYLE == 'paper' else r"4\sin\big(k\cdot(w-w')\big)"
    m = re.search(r"&=-" + re.escape(PREF) + r"\\notag\\\\\n&\\quad\\times\\int_\{x,y,z,w,x',w'\}" + re.escape(SIN + r"\," + KTn) + r"\\notag\\\\\n&\\quad\\times\\Big\[(.*?)\\Big\]\\,o\(", bs, re.S)
    assert m, 'sin form not parsed'
    brs = parse_bracket(m.group(1))
    worst = 0
    for N in (3, 4):
        G = SUN(N); rng = np.random.default_rng(5 + N)
        W = ['w%d' % i for i in range(5)]
        U = {p: G.randU(rng) for p in ['x', 'y', "x'", 'z'] + W}
        pos = {p: rng.normal(size=2) for p in W}; k = rng.normal(size=2) * 1.7
        KT = {(a, c): rng.normal() for a in W for c in W}
        lhs = rhs = 0
        for a in W:
            for c in W:
                ph = np.exp(-1j * k @ (pos[a] - pos[c])); sub = {'w': a, "w'": c}
                for t in terms:
                    kv = KT[(a, c)] if t['ker'] == 'KT' else KT[(c, a)]
                    lhs += ph * t['coef'] * kv * evalbr(t['br'], G, U, sub)
                rhs += -4 * np.sin(k @ (pos[a] - pos[c])) * KT[(a, c)] * evalbr(brs, G, U, sub)
        r = abs(lhs - rhs) / abs(rhs); worst = max(worst, r)
        print('  SU(%d) II.I: exp-form %+.12f%+.2ei   sin-form %+.12f' % (N, lhs.real, lhs.imag, rhs.real))
    print('II.I sin form vs two-line form: worst rel. deviation %.1e' % worst)
if STYLE == 'paper':
    sb = eqbody(sec, 'sodd')
    assert r"-\frac{\mathcal C_{\rm f}}{N_c^2-1}\,\frac{N_c}{N_c^2-4}\int_{x,y,z,w,x',w'}4\sin\big(\v k\cdot(\v w-\v w')\big)\,\mathfrak K_{\mathcal T}(w,w')\,\Delta_{\rm II.I}" in sb
    w2 = 0
    for N in (3, 4):
        G = SUN(N); rng = np.random.default_rng(99 + N); U = {p: G.randU(rng) for p in POINTS}
        a, c = evalbr(brs, G, U), evalbr(DEL['II.I'], G, U); w2 = max(w2, abs(a - c) / max(1, abs(a)))
    print('eq:sodd II.I line vs sin form: worst rel. deviation %.1e' % w2)
