r"""The three-rho terms proportional to log(k+/Lambda) and to log(vee/k+), taken directly from the J and B parts of your
rows (gen.py) and from 4rho.tex, cancel?

Pointwise test (the charges are integration variables): the symmetric part of a product of three charges is evaluated
with classical charges; each term is summed over the 6 placements of its three charge points x, y, x' and, for Groups I
and III (written 'with c.c.'), its c.c. is added: conj(term) with w <-> w' (that restores the phase e^{-ik(w'-w)}).
Group II rows are self-conjugate as written.  The large logarithm is set to 1 (Lambda = k+/e, vee = e k+).
Which terms are three-rho is decided numerically (rho -> 2 rho multiplies the term by 8).

Variants: Group I Row III with its edge-kernel term (as written) or with D (fix F3); the 4rho.tex three-rho lists with
or without c.c.  Output: three_check.json"""
import itertools
import json
import os
import re
import numpy as np
import blocks as BL
import gen as GN
import texeval as TE
import warnings
warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
E = np.e
CH = ('x', 'y', 'xp')
LKL, LVK = GN.LKL, GN.LVK

# fix F3: the edge kernel of Group I Row III (first log term of d^3N_3a) -> D
D_TERM = (r"-\frac{ig^{4}f^{abd}}{8\pi^{5}}\,\frac{1}{k^{+}}\,%s\,\int_{x,y,z,x'}"
          r"\frac{(x'-w')\cdot(y-w)\,(x-z)\cdot(w-z)}{(x'-w')^{2}\,(y-w)^{2}\,(x-z)^{2}\,(w-z)^{2}}\notag \\ & "
          r"\left[U^{ab'}(x')\,\rho^{b'}(x')-U^{ab'}(w')\, \rho^{b'}(x')\right]\times "
          r"\left[U^{de}(x)\,U^{bc}(w)\, \rho^{e}(x)\,\rho^{c}(y)-U^{de}(x)\,U^{bc}(y)\,\rho^{e}(x)\,\rho^{c}(y)\right]")


# Group I Row II part 2a, lines 2 and 3: the z integral is done there (1/pi^4).  For the pointwise test it is written
# back:  log[(a-y)^2/(b-y)^2] = -(1/pi) int_z [K(a,y,z) - K(b,y,z)],  K(a,b,z) = (a-z).(b-z)/((a-z)^2 (b-z)^2)
Z_FORMS = [(r"\log\frac{(w-y)^{2}}{(x-y)^{2}}",
            r"\left(-\frac{1}{\pi}\right)\left[\frac{(w-z)\cdot(y-z)}{(w-z)^{2}\,(y-z)^{2}}-\frac{(x-z)\cdot(y-z)}{(x-z)^{2}\,(y-z)^{2}}\right]"),
           (r"\log\frac{(w-x)^{2}}{(x-y)^{2}}",
            r"\left(-\frac{1}{\pi}\right)\left[\frac{(w-z)\cdot(x-z)}{(w-z)^{2}\,(x-z)^{2}}-\frac{(y-z)\cdot(x-z)}{(y-z)^{2}\,(x-z)^{2}}\right]")]


def z_form(t):
    for a, b in Z_FORMS:
        t = t.replace(a, b)
    return t


def row_terms(part, F3=True):
    """[(label, full LaTeX, add_cc)] of the J ('J') or B ('B') parts of all rows"""
    BLK, SPL = GN.all_splits()
    out = []
    for name, B in BLK.items():
        for n, s in enumerate(SPL[name]):
            t = s.get(part)
            if t is None:
                continue
            if F3 and part == 'J' and (name, n) == ('G1 III N3a', 1):
                t = D_TERM % LKL
            if name == 'G1 II(2a)':
                t = z_form(t)
            out.append(('%s #%d' % (name, n), B.prefix + ' (' + t + ') ' + B.suffix, not name.startswith('G2')))
    return out


def hidden_B():
    """the terms (B) of BFKL_final.pdf (region y' -> z of Group III Rows III(1), III(2), V, VI), y' renamed w'"""
    S = open(os.path.join(HERE, '..', 'bfkl_final', 'static.tex')).read()
    i = S.index(r'\paragraph{The new terms}')
    blk = S[S.index(r'\begin{align}', i) + len(r'\begin{align}'):S.index(r'\end{align}', i)]
    parts = re.split(r'\\text\{Row ([IV()12]+):\}\\quad&', blk)[1:]
    out = []
    for k in range(0, len(parts), 2):
        t = parts[k + 1].replace(r'\Lv', LVK).replace(r'\nc', '')
        t = re.sub(r'[,.]\s*(\\\\)?\s*$', '', t.strip())
        t = t.replace("y'", "w'")
        out.append(('(B) Row ' + parts[k], t, True))
    return out


def four_rho_three(part):
    """three-rho lists of 4rho.tex (the final, rewritten form of 'Two Wilson lines' and 'Three Wilson lines')"""
    L = open(os.path.join(HERE, 'input', '4rho.tex')).read()
    big = LKL if part == 'J' else LVK
    out = []
    for head in (r'We rewrite them (by some simplification)', r'Rewriting with some simplification, we get'):
        i = L.index(head)
        t = L[L.index(r'\begin{align}', i) + len(r'\begin{align}'):L.index(r'\end{align}', i)]
        t = t[re.search(r'&\s*=', t).end():]
        t = t.replace(r'\log{\left(\frac{\vee}{\Lambda}\right)}', big)
        out.append(('4rho.tex ' + head.split()[0] + ' ' + head.split()[1], t))
    return out


def configs(seed):
    c = TE.cfg(seed)
    return c


def permuted(c, perm):
    """the configuration with the charge points (with their U and rho) placed according to perm"""
    new = type(c).__new__(type(c))
    new.__dict__.update(c.__dict__)
    new.pt, new.U, new.rho = dict(c.pt), dict(c.U), dict(c.rho)
    for a, b in zip(CH, perm):
        new.pt[a], new.U[a], new.rho[a] = c.pt[b], c.U[b], c.rho[b]
    return new


def swapped_w(c):
    new = permuted(c, CH)
    new.pt['w'], new.pt['wp'] = c.pt['wp'], c.pt['w']
    new.U['w'], new.U['wp'] = c.U['wp'], c.U['w']
    new.rho['w'], new.rho['wp'] = c.rho['wp'], c.rho['w']
    return new


def scaled(c, lam):
    new = permuted(c, CH)
    new.rho = {k: lam * v for k, v in c.rho.items()}
    return new


def env(c):
    return TE.Env(c, KP=1.0, V=E, L=1.0 / E)


def degree(tex, c):
    a = TE.term_value(tex, env(c))
    b = TE.term_value(tex, env(scaled(c, 2.0)))
    if abs(a) < 1e-14:
        return None
    r = abs(b / a)
    for d in (2, 3, 4):
        if abs(r - 2**d) < 1e-6 * 2**d:
            return d
    raise ValueError('not homogeneous in rho: ratio %s' % r)


def frame(tex):
    """relabelling of the integration variables of a term onto w, w' (phase e^{-ik(w'-w)}), z and the charges x, y, x'"""
    code = TE.translate(tex)
    ph = re.findall(r'KD\("([^"]*)","([^"]*)"\)', code)
    assert len(set(ph)) == 1, ph
    a, p = ph[0]
    charges = []
    for q in re.findall(r'CR\("[^"]*","([^"]*)"\)', code):
        if q not in charges:
            charges.append(q)
    pts = set(re.findall(r'D\("([^"]*)","([^"]*)"\)', code))
    pts = {q for pair in pts for q in pair} | set(re.findall(r'CU\("[^"]*","[^"]*","([^"]*)"\)', code)) | {a, p} | set(charges)
    rest = sorted(pts - {a, p} - set(charges))
    assert len(charges) == 3 and len(rest) <= 1, (charges, rest)
    m = {a: 'wp', p: 'w'}
    m.update(dict(zip(charges, CH)))
    if rest:
        m[rest[0]] = 'z'
    return m


def relabel(c, m):
    new = permuted(c, CH)
    for n, cn in m.items():
        new.pt[n], new.U[n], new.rho[n] = c.pt[cn], c.U[cn], c.rho[cn]
    return new


def sym_value(tex, c, cc):
    m = frame(tex)
    v = 0.0
    for perm in itertools.permutations(CH):
        cp = permuted(c, perm)
        v += TE.term_value(tex, env(relabel(cp, m)))
        if cc:
            v += np.conj(TE.term_value(tex, env(relabel(swapped_w(cp), m))))
    return v


def values(part, seeds):
    """per term and seed: (direct, c.c. part) of the symmetrized value, in units of N e^{-ik(w'-w)}, N = g^4/((2pi)^3 16pi^5 k+)"""
    terms = row_terms(part, True)
    if part == 'J':
        terms += [('G1 III N3a #1 (edge kernel, as written)', t, cc) for lab, t, cc in row_terms('J', False)
                  if lab == 'G1 III N3a #1']
    else:
        terms += hidden_B()
    c0 = configs(seeds[0])
    three = [(lab, t, cc) for lab, t, cc in terms if degree(t, c0) == 3]
    three += [(lab, t, None) for lab, t in four_rho_three(part)]
    out = {}
    for seed in seeds:
        c = configs(seed)
        unit = TE.fp.G**4 / ((2 * np.pi)**3 * 16 * np.pi**5) * np.exp(-1j * (c.kv @ (c.pt['wp'] - c.pt['w'])))
        for lab, t, cc in three:
            a = sym_value(t, c, False) / unit
            b = sym_value(t, c, True) / unit - a
            out.setdefault(lab, {})[seed] = (complex(a), complex(b), cc)
    return out


def total(vals, seed, labels, four_cc=False):
    tot, big = 0.0, 0.0
    for lab in labels:
        a, b, cc = vals[lab][seed]
        v = a + (b if (cc or (cc is None and four_cc)) else 0.0)
        tot += v
        big = max(big, abs(v))
    return abs(tot), big


if __name__ == '__main__':
    seeds = (11, 12, 13)
    out = {}
    for part in ('J', 'B'):
        vals = values(part, seeds)
        base = [lab for lab in vals if 'as written' not in lab]
        variants = {'as in the note': base}
        if part == 'J':
            variants['edge kernel instead of D'] = [l for l in base if l != 'G1 III N3a #1'] + \
                ['G1 III N3a #1 (edge kernel, as written)']
        else:
            variants['without the terms (B)'] = [l for l in base if not l.startswith('(B)')]
        res = {}
        for vname, labels in variants.items():
            for four_cc in (False, True):
                key = vname + (', 4rho.tex lists with an extra c.c.' if four_cc else '')
                r = [total(vals, s, labels, four_cc) for s in seeds]
                res[key] = dict(relative=max(t / b for t, b in r), residual=[t for t, b in r], largest=[b for t, b in r])
                print('%s | %-60s worst rel. residual %.2e' % (part, key, res[key]['relative']))
        table = {lab: dict(direct=v[seeds[0]][0].imag, cc=v[seeds[0]][1].imag, cc_added=v[seeds[0]][2]) for lab, v in vals.items()}
        out[part] = dict(variants=res, seed11_values_imag=table, terms=base)
    json.dump(out, open(os.path.join(HERE, 'three_check.json'), 'w'), indent=1)
