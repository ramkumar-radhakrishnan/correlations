# s16c_sec8_colour.py -- colour structures of the new Sec. 8 ("Final result with factorised Wilson-line structures"):
#  (1) every factorised colour tensor AS PRINTED (TeX parsed by s16_lib) == the verified colour tensor of the row
#      (s3_colour.STRUCT; random SU(3), SU(4) adjoint Wilson lines) -- exact identities, no large-N_c;
#  (2) the 'leg' form of s16_sec8data (from which the singlet traces are generated) == the same tensor;
#  (3) the singlet forms are GENERATED here from the legs:
#        (3rho)    t C^{abc} d^{abc} = (t i) sum coef Tr{T^h O^T L_3 D^h L_2^T O}      (G-type, O a Wilson-line leg)
#        (contact) C^{abc} f^{abc}|_{p_i=p_j} = sum coef sgn Tr{T^e O^T L_3 T^e L_2^T O} (F-type)
#      written to s16c_singlet.json, and the PRINTED TeX is parsed and compared with direct contraction;
#  (4) the genuine two-rho colour factors (QUAD8): printed tensor == verified, printed singlet trace == delta-contraction;
#  (5) the generic identities used: F3(O,B,C) = Tr{T^e O^T C T^e B^T O}, G3(O,B,C) = i Tr{T^h O^T C D^h B^T O}.
# Run: python3 -I s16c_sec8_colour.py -> s16c_sec8_colour.log, s16c_singlet.json
import os, sys, json, itertools, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from s16_lib import SUN, ev, evs, randUs
from s3_colour import STRUCT, E
from s16_sec8data import CUBIC8, QUAD8, TEXPT
from singlet_data import CUBIC
out = []
def pr(*a): s = ' '.join(str(x) for x in a); print(s); out.append(s)
SF = {s[0]: s for s in STRUCT}
TT = {'ord': 1, 'r{rr}': 2, '{rr}r': 2}
GEN = {'ord': {'12': 0.5j, '23': 0.5j, '13': -0.5j}, 'r{rr}': {'12': 1j, '13': -1j}, '{rr}r': {'23': 1j, '13': -1j}}
# contacts (pair, from, to, n-args) in the author's Sec. 7 conventions, from singlet_data.CUBIC
CONT = {}
for row, R in CUBIC.items():
    for lab, contacts, *_ in R['structs']:
        CONT[lab] = [(c[0], c[1], c[2], c[5]) for c in contacts]
# ---------------- leg algebra ----------------
def legmat(leg, U, n):
    M = np.eye(n)
    for t in leg:
        if t[0] == 'U': A = U[t[1]]
        elif t[0] == 'Ut': A = U[t[1]].T
        elif t[0] == 'D': A = U[t[1]] - U[t[2]]
        elif t[0] == 'Dt': A = (U[t[1]] - U[t[2]]).T
        M = M @ A
    return M
def legs_tensor(S, U, legs):
    T = 0
    for coef, (L1, L2, L3), q in legs:
        M = [legmat(L, U, S.n) for L in (L1, L2, L3)]
        sub = 'ijk,i%s,j%s,k%s->abc' % tuple('abc'[x - 1] for x in q)
        T = T + coef * np.einsum(sub, S.f, *M, optimize=True)
    return T
def tr_leg(L):
    o = []
    for t in reversed(L):
        o.append((({'U': 'Ut', 'Ut': 'U', 'D': 'Dt', 'Dt': 'D'})[t[0]],) + tuple(t[1:]))
    return o
def ortho(L): return all(t[0] in ('U', 'Ut') for t in L)
def subst(L, fr, to):
    o = []
    for t in L:
        t = tuple(to if (k > 0 and x == fr) else x for k, x in enumerate(t)); o.append(t)
    return o
def zero(L): return any(t[0] in ('D', 'Dt') and t[1] == t[2] for t in L)
def simplify(ch):
    ch = list(ch); changed = True
    while changed and len(ch) > 1:
        changed = False
        for i in range(len(ch)):
            j = (i + 1) % len(ch)
            if i == j: break
            a, b = ch[i], ch[j]
            if a[0] in ('U', 'Ut') and b[0] in ('U', 'Ut') and a[0] != b[0] and a[1] == b[1]:
                ch = [c for k, c in enumerate(ch) if k not in (i, j)]; changed = True; break
    return ch
def canon(ch):
    """rotate so that the trace starts with the first generator (or with the first difference)"""
    if not ch: return ch
    for k, t in enumerate(ch):
        if t[0] in ('T', 'Dg'): return ch[k:] + ch[:k]
    for k, t in enumerate(ch):
        if t[0] in ('Dt', 'Ut'): return ch[k:] + ch[:k]
    return ch
def tex_chain(ch):
    s = ''
    for t in ch:
        if t[0] == 'T': s += r'T^{%s}' % t[1]
        elif t[0] == 'Dg': s += r'D^{%s}' % t[1]
        elif t[0] == 'U': s += r'U(%s)' % TEXPT[t[1]]
        elif t[0] == 'Ut': s += r'U^{\dagger}(%s)' % TEXPT[t[1]]
        elif t[0] == 'D': s += r'\big(U(%s)-U(%s)\big)' % (TEXPT[t[1]], TEXPT[t[2]])
        elif t[0] == 'Dt': s += r'\big(U^{\dagger}(%s)-U^{\dagger}(%s)\big)' % (TEXPT[t[1]], TEXPT[t[2]])
    return r'{\rm Tr}\Big\{' + s + r'\Big\}'
def trace_term(L1, L2, L3, kind):
    """G3 (kind 'G', generator letter h, D insertion) or F3 (kind 'F', T T) of three legs -> (factor, chain, ncpow)"""
    legs = [L1, L2, L3]
    k = next((n for n, L in enumerate(legs) if ortho(L)), None)
    if k is None: raise ValueError('no Wilson-line leg')
    if kind == 'G':
        O, B, C = legs[k], legs[(k + 1) % 3], legs[(k + 2) % 3]     # cyclic rotation: G3 invariant
        ch = [('T', 'h')] + tr_leg(O) + C + [('Dg', 'h')] + tr_leg(B) + O
        fac = 1j
    else:
        rest = [L for n, L in enumerate(legs) if n != k]; O = legs[k]
        B, C = rest[0], rest[1]                                      # F3 totally symmetric
        ch = [('T', 'e')] + tr_leg(O) + C + [('T', 'e')] + tr_leg(B) + O
        fac = 1
    ch = simplify(ch); ncp = 0
    # T^e T^e adjacent (cyclically) -> N_c
    gi = [n for n, t in enumerate(ch) if t[0] in ('T', 'Dg')]
    if len(gi) == 2 and (gi[1] - gi[0] == 1 or (gi[0] == 0 and gi[1] == len(ch) - 1)):
        if all(ch[n][0] == 'T' for n in gi):
            ch = [t for n, t in enumerate(ch) if n not in gi]; ncp = 1
        else:
            return None
    return fac, canon(ch), ncp
def perm_sign(q):
    s = 1; q = list(q)
    for i in range(3):
        for j in range(i + 1, 3):
            if q[i] > q[j]: s = -s
    return s
def combine(terms):
    acc = {}
    for c, ch, ncp in terms:
        key = (tuple(ch), ncp); acc[key] = acc.get(key, 0) + c
    return [(c, list(k[0]), k[1]) for k, c in acc.items() if abs(c) > 1e-12]
def tex_sum(terms):
    s = ''
    for c, ch, ncp in terms:
        c = complex(c); assert abs(c.imag) < 1e-12; c = c.real
        sg = '-' if c < 0 else '+'; a = abs(c)
        num = '' if abs(a - 1) < 1e-12 else ('%d' % round(a) if abs(a - round(a)) < 1e-12 else repr(a))
        s += sg + num + ('N_{c}' if ncp else '') + tex_chain(ch)
    return s[1:] if s.startswith('+') else s
def singlet_3rho(legs):
    terms = []
    for coef, (L1, L2, L3), q in legs:
        r = trace_term(L1, L2, L3, 'G')
        if r is None: continue
        fac, ch, ncp = r
        terms.append((coef, ch, ncp))   # value of the term = coef * i * Tr
    return combine(terms)
def singlet_contact(legs, fr, to):
    terms = []
    for coef, (L1, L2, L3), q in legs:
        Ls = [subst(L, fr, to) for L in (L1, L2, L3)]
        if any(zero(L) for L in Ls): continue
        r = trace_term(*Ls, 'F')
        if r is None: continue
        fac, ch, ncp = r
        terms.append((coef * perm_sign(q), ch, ncp))
    return combine(terms)
# ---------------- checks ----------------
REFQ = {}
def refq(name, S, U):
    f = S.f; L = U['xp'] - U['wp']
    if name == 'I.II N1': return E('ab,ac->bc', L, U['w'] - U['x'])
    if name == 'I.V N2': return E('ab,af->bf', L, S.N * U['y']) - E('aB,adc,fbe,db,ce->Bf', L, f, f, U['x'], U['z'])
    if name == 'I.V 11/6': return E('ab,af->bf', L, U['y'] - U['x'])
    if name == 'II.I second':
        return (E('kbB,lbd,ak,al->Bd', f, f, U['wp'], U['w']) - E('kbB,acd,ak,cb,dE->BE', f, f, U['wp'], U['z'], U['x'])
                - E('akd,cbq,ac,kb,dB->Bq', f, f, U['w'], U['z'], U['xp']) + S.N * E('dB,dE->BE', U['xp'], U['x']))
    if name == 'II.IV': return E('ca,cb->ab', U['yp'] - U['xp'], U['y'] - U['x'])
    if name == 'II.IV rem': return E('ca,cb->ab', U['z'] - U['xp'], U['z'] - U['x'])
    if name.startswith('III.III second'):
        v = 'yp' if "y'" in name else 'z'
        Lv = E('db,cB->cdbB', U['z'], U[v] - U['xp'])
        R2 = E('cbd,ac->abd', f, U['w']) - E('acd,cb,de->abe', f, U['z'], U['x'])
        return E('cda,cdbB,abh->Bh', f, Lv, R2)
    raise KeyError(name)
SING = {'cubic': {}, 'quad': {}}
worst = {}
def upd(k, v): worst[k] = max(worst.get(k, 0), v)
for N in (3, 4):
    S = SUN(N); rng = np.random.default_rng(1700 + N)
    for _ in range(4):
        U = randUs(S, rng)
        # (5) generic identities
        A, B, C = [U[p] - U[q] for p, q in (('x', 'y'), ('z', 'w'), ('xp', 'wp'))]; O = U['yp']
        F3 = lambda P, Q, R: np.einsum('ijk,lmn,il,jm,kn->', S.f, S.f, P, Q, R, optimize=True)
        G3 = lambda P, Q, R: np.einsum('ijk,lmn,il,jm,kn->', S.f, S.d, P, Q, R, optimize=True)
        trF = np.einsum('eij,jk,kl,elm,mn,ni->', S.T, O.T @ C, np.eye(S.n), S.T, B.T @ O, np.eye(S.n), optimize=True)
        trG = 1j * np.einsum('hij,jk,hkl,li->', S.T, O.T @ C, S.D, B.T @ O, optimize=True)
        upd('(5) F3(O,B,C) = Tr{T O^T C T B^T O}', abs(F3(O, B, C) - trF))
        upd('(5) G3(O,B,C) = i Tr{T O^T C D B^T O}', abs(G3(O, B, C) - trG))
        upd('(5) F3 totally symmetric', abs(F3(O, B, C) - F3(B, C, O)) + abs(F3(O, B, C) - F3(O, C, B)))
        upd('(5) G3 cyclic, antisymmetric', abs(G3(O, B, C) - G3(B, C, O)) + abs(G3(O, B, C) + G3(O, C, B)))
        for lab, D8 in CUBIC8.items():
            name, fn, typ, pts = SF[lab]
            ref = fn(S, U); nrm = max(1, np.abs(ref).max())
            v, idx = ev(D8['tex'], S, U)
            assert sorted(idx) == ['a', 'b', 'c'], (lab, idx)
            v = np.transpose(v, [idx.index(q) for q in 'abc'])
            upd('(1) printed factorised C^{abc} == verified: ' + lab, np.abs(v - ref).max() / nrm)
            upd('(2) leg form == verified: ' + lab, np.abs(legs_tensor(S, U, D8['legs']) - ref).max() / nrm)
            # (3) singlets
            s3 = singlet_3rho(D8['legs']); t = TT[typ]
            tex3 = tex_sum(s3)
            exp = t * E('abc,abc->', ref, S.d)
            got = t * 1j * evs(tex3, S, U)
            upd('(3) 3rho singlet: (t i) x printed bracket == t C.d: ' + lab, abs(got - exp) / max(1, abs(exp)))
            SING['cubic'].setdefault(lab, {})['3rho'] = dict(coef=('%+di' % t).replace('+1i', '+i'), tex=tex3, args=[TEXPT[p] for p in pts])
            for pair, fr, to, nargs in CONT[lab]:
                sc = singlet_contact(D8['legs'], fr, to); texc = tex_sum(sc)
                U2 = dict(U); U2[fr] = U[to]
                Phi = E('abc,abc->', fn(S, U2), S.f)
                got = evs(texc, S, U) if texc else 0.0
                upd('(3) contact singlet (printed bracket == C^{abc}f^{abc} at p_i=p_j): %s %s' % (lab, pair), abs(got - Phi) / max(1, abs(Phi)))
                g = GEN[typ][pair]
                SING['cubic'][lab].setdefault('contacts', {})[pair] = dict(gen=repr(g), fr=TEXPT[fr], to=TEXPT[to], n=[TEXPT[p] for p in nargs], tex=texc)
        for name, Q in QUAD8.items():
            ref = refq(name, S, U); nrm = max(1, np.abs(ref).max())
            v, idx = ev(Q['tex'], S, U)
            assert sorted(idx) == sorted(Q['idx']), (name, idx)
            v = np.transpose(v, [idx.index(q) for q in Q['idx']])
            upd('(4) printed two-rho colour factor == verified: ' + name, np.abs(v - ref).max() / nrm)
            upd('(4) printed singlet trace == delta-contraction: ' + name, abs(evs(Q['sing'], S, U) - np.trace(ref)) / max(1, abs(np.trace(ref))))
for k, v in worst.items(): pr('%-110s max rel dev %.1e' % (k, v))
# controls: a wrong sign in one printed difference must be detected
S = SUN(3); rng = np.random.default_rng(5); U = randUs(S, rng)
bad = CUBIC8['II.I 1a']['tex'].replace(r'f^{kbc}-f^{gkl}', r'f^{kbc}+f^{gkl}', 1)
v, idx = ev(bad, S, U); v = np.transpose(v, [idx.index(q) for q in 'abc']); ref = SF['II.I 1a'][1](S, U)
pr('control: II.I 1a with one sign flipped: rel dev %.2f (must be O(1))' % (np.abs(v - ref).max() / np.abs(ref).max()))
json.dump(SING, open(os.path.join(HERE, 's16c_singlet.json'), 'w'), indent=1)
open(os.path.join(HERE, 's16c_sec8_colour.log'), 'w').write('\n'.join(out) + '\n')
