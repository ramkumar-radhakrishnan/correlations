#!/usr/bin/env python3
# gen_rows.py -- generates ../rows/*.tex (the finite terms, row by row) mechanically from Sec. 8 of
# nlo_review/final_explicit/Final_terms_explicit_2_fixed.tex:
#   * rho-explicit form (labels eq:H-...), colour-singlet form (labels eq:S-...), kernel definitions (eq:K-...)
#   * only notation is changed: \mathcal{K} -> \Kv, names of the left-hand sides, the transverse/longitudinal
#     kernel blocks that are repeated in several equations are replaced by named kernels, lines are re-broken.
#   * a round-trip check expands the kernels again and compares every equation with the source (normal form).
# Run: python3 -I scripts/gen_rows.py   (needs pdflatex for the line breaking)
import os, re, sys, subprocess, json, tempfile, shutil
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from rowlib import *
OUT = os.path.join(HERE, '..', 'rows'); os.makedirs(OUT, exist_ok=True)
F1, F2 = load_sec8()

# ---------------------------------------------------------------- read the equations
def read(T, prefix):
    eqs = []
    for a in aligns(T):
        for e in split_equations(a):
            s, lab = flat(e); s, _ = strip_wrappers(s)
            if '=' not in s or r'\equiv' in s.split('=')[0]: continue      # the definition C_{I.V} (not an equation of the result)
            L, R = top_eq_split(s)
            if lab is None:
                m = re.search(r'II\.I\}_\{2([a-e])\}', L); assert m, L
                lab = 'eq:%s-II.I-2%s' % (prefix, m.group(1))
            eqs.append([lab, L, R])
    return eqs
E1 = read(F1, 'H'); E2 = read(F2, 'S')
assert len(E1) == 25 and len(E2) == 25, (len(E1), len(E2))
SRC_RHS = {lab: R for lab, L, R in E1 + E2}

def cosmetic(s):
    s = re.sub(r'(\\Bigg\[|\\bigg\[)\s*\\times', r'\1', s)       # "[ x {kernel}" -> "[ {kernel}"
    s = re.sub(r'\\times\s*\\,\s*', r'\\times ', s)
    s = re.sub(r'\\times\s*\\times', r'\\times', s)
    return s
def cnorm(s):
    s = norm(cosmetic(flat(s)[0]))
    s = s.replace(r'\times', '').replace(r'\,', '')          # explicit multiplication signs and thin spaces are markup
    return s

# ---------------------------------------------------------------- kernels
KER = {}      # name -> content (inside the outer braces, verbatim)
KDELIM = {}   # name -> (open, close) removed around the content
def defk(name, eqlab, opener, closer, nth):
    R = SRC_RHS[eqlab]; i, j = block(R, opener, closer, nth)
    c = R[i+len(opener):j-len(closer)].strip()
    if name in KER: assert cnorm(KER[name]) == cnorm(c), name
    else: KER[name] = c; KDELIM[name] = (opener, closer)
B = (r'\bigg\{', r'\bigg\}'); BB = (r'\Bigg\{', r'\Bigg\}')
defk(r'\mathfrak{K}_{\rm I.II}', 'eq:H-I.II-3', *B, 0)
defk(r'\mathfrak{K}^{(z)}_{\rm I.III}', 'eq:H-I.III-3', *B, 0)
defk(r'\mathfrak{K}^{(w)}_{\rm I.III}', 'eq:H-I.III-3', *B, 1)
defk(r'\mathfrak{K}_{\rm I.V}', 'eq:H-I.V-3', *B, 0)
defk(r'\mathfrak{K}^{\rm quad}_{\rm I.V}', 'eq:H-I.V-2o', *BB, 0)
defk(r'\mathfrak{K}^{(1)}_{\rm II.I}', 'eq:H-II.I-3', *B, 0)
defk(r'\mathfrak{K}^{(2)}_{\rm II.I}', 'eq:H-II.I-3', *B, 1)
defk(r'\mathfrak{K}_{\rm III.II}', 'eq:H-III.II-3', *B, 0)
for r in ('III', 'V', 'VI'):
    defk(r"\mathfrak{K}^{(y')}_{\rm III}", 'eq:H-III.%s-3' % r, *B, 0)
    defk(r'\mathfrak{K}^{(z)}_{\rm III}', 'eq:H-III.%s-3' % r, *B, 1)
# genuine two-rho kernels of Row II.I (a-d): everything between the phase and the colour operator
for x in 'abcd':
    R = SRC_RHS['eq:H-II.I-2' + x]
    a = R.index(r"e^{-ik\cdot(w-w')}\,") + len(r"e^{-ik\cdot(w-w')}\,")
    b = R.index(r"\times\bigg[\Big[U^{ac'}(w')")
    KER[r'\mathfrak{K}^{{\rm quad},%s}_{\rm II.I}' % x] = R[a:b].strip(); KDELIM[r'\mathfrak{K}^{{\rm quad},%s}_{\rm II.I}' % x] = ('', '')
KUSE = {   # which kernels are substituted in which rows
 'I.II': [r'\mathfrak{K}_{\rm I.II}'], 'I.III': [r'\mathfrak{K}^{(z)}_{\rm I.III}', r'\mathfrak{K}^{(w)}_{\rm I.III}'],
 'I.V': [r'\mathfrak{K}_{\rm I.V}', r'\mathfrak{K}^{\rm quad}_{\rm I.V}'],
 'II.I': [r'\mathfrak{K}^{(1)}_{\rm II.I}', r'\mathfrak{K}^{(2)}_{\rm II.I}'] + [r'\mathfrak{K}^{{\rm quad},%s}_{\rm II.I}' % x for x in 'abcd'],
 'II.IV': [], 'III.II': [r'\mathfrak{K}_{\rm III.II}'],
 'III.III': [r"\mathfrak{K}^{(y')}_{\rm III}", r'\mathfrak{K}^{(z)}_{\rm III}'],
 'III.V': [r"\mathfrak{K}^{(y')}_{\rm III}", r'\mathfrak{K}^{(z)}_{\rm III}'],
 'III.VI': [r"\mathfrak{K}^{(y')}_{\rm III}", r'\mathfrak{K}^{(z)}_{\rm III}'],
}
def row_of(lab): return re.match(r'eq:[HS]-([IV]+\.[IV]+)', lab).group(1)

def substitute(lab, R):
    """replace kernel blocks by their names (the blocks are found as strings, after removal of spaces)"""
    used = []
    for name in KUSE[row_of(lab)]:
        o, c = KDELIM[name]; full = o + KER[name] + c
        # locate by a whitespace-insensitive search
        pat = re.escape(re.sub(r'\s+', '', full))
        pat = ''.join(ch if ch == '\\' else ch for ch in pat)
        # build regex allowing optional spaces between characters of the compact form
        comp = re.sub(r'\s+', '', full)
        rx = r'\s*'.join(re.escape(ch) for ch in comp)
        R, n = re.subn(rx, lambda m: ' ' + name + ' ', R)
        if n: used.append((name, n))
    return R, used

def expand(R):
    for name in sorted(KER, key=len, reverse=True):
        o, c = KDELIM[name]
        R = R.replace(name, o + KER[name] + c)
    return R

# ---------------------------------------------------------------- left-hand sides
def new_lhs(lab):
    row = row_of(lab); kind = lab.split('-')[-1]; sing = lab.startswith('eq:S')
    sub = {'3': r'{\rm sym}', '2': r'{\rm con}', '2o': r'{\rm quad}'}.get(kind)
    if row == 'II.IV': sub = r'{\rm fin}'
    if sub is None:
        m = re.match(r'2([a-e])$', kind); sub = r'{\rm quad},%s' % m.group(1)
    if row == 'II.IV': sub = r'{\rm fin}'
    L = r'\frac{d^{3}N^{\rm %s}_{%s}}{d^{3}k}' % (row, sub)
    return (r'\Big\langle' + L + r'\Big\rangle') if sing else L
SRC_LHS_KIND = {'3': r'3\rho', '2': r'2\rho'}

# ---------------------------------------------------------------- line breaking
TOKRX = re.compile(r'\\[A-Za-z]+|\\.|.', re.S)
def candidates(s):
    r"""break candidates: list of (position, kind, bracket depth); kind 1 = before binary +/-, \times; 2 = after \, (fallback)"""
    out = []; bd = 0; lr = 0; rd = 0; sq = 0; prev = None; pos = 0; opens = []
    for m in TOKRX.finditer(s):
        t = m.group(0); pos = m.start()
        if t.isspace(): continue
        if bd == 0 and lr == 0:
            if t in '+-' and rd == 0 and prev not in (None, '[', r'\{', '(', '=', r'\times', ',', '+', '-', r'\cdot', '&'):
                inUdiff = bool(opens) and s[opens[-1]+1:opens[-1]+3] == 'U^' and s[m.end():].lstrip().startswith('U^')
                if not inUdiff: out.append((pos, 1, sq))
            elif t == r'\times' and rd == 0 and prev is not None:
                out.append((pos, 1, sq))
            elif prev == r'\,' and rd == 0 and t not in ('+', '-', r'\times'):
                out.append((pos, 2, sq))
            elif t == r'\int' and rd == 0 and prev not in (None, r'\,'):
                out.append((pos, 2, sq))
        if t == '{': bd += 1
        elif t == '}': bd -= 1
        elif t == r'\left': lr += 1
        elif t == r'\right': lr -= 1
        elif bd == 0:
            if t == '(': rd += 1
            elif t == ')': rd -= 1
            elif t in ('[', r'\{'): sq += 1; opens.append(m.start() + len(t) - 1)
            elif t in (']', r'\}'): sq -= 1; opens.pop()
        prev = t
    return out

class Measurer:
    def __init__(self):
        self.cache = {}
        self.cachefile = os.path.join(HERE, '.widths.json')
        if os.path.exists(self.cachefile): self.cache = json.load(open(self.cachefile))
        main = open(os.path.join(HERE, '..', 'main.tex')).read()
        self.pre = main[:main.index(r'\title{')]
    def measure(self, atoms):
        need = [a for a in dict.fromkeys(atoms) if a not in self.cache]
        if need:
            d = tempfile.mkdtemp()
            for f in ('jheppub.sty', 'tikzit.sty', 'stiles.tikzstyles'):
                shutil.copy(os.path.join(HERE, '..', f), d)
            body = [self.pre, r'\begin{document}', r'\typeout{TW:\the\textwidth}']
            for k, a in enumerate(need):
                body.append(r'\setbox0\hbox{$\displaystyle %s$}\typeout{W:%d:\the\wd0}' % (a, k))
            body.append(r'\end{document}')
            open(os.path.join(d, 'm.tex'), 'w').write('\n'.join(body))
            subprocess.run(['pdflatex', '-interaction=nonstopmode', 'm.tex'], cwd=d, stdout=subprocess.DEVNULL)
            log = open(os.path.join(d, 'm.log'), errors='ignore').read()
            log = re.sub(r'\n', '', log) if False else log
            tw = re.search(r'TW:([\d.]+)pt', log); self.cache['__tw__'] = float(tw.group(1))
            got = dict((int(a), float(b)) for a, b in re.findall(r'W:(\d+):([\d.]+)pt', log))
            assert len(got) == len(need), (len(got), len(need), d)
            for k, a in enumerate(need): self.cache[a] = got[k]
            json.dump(self.cache, open(self.cachefile, 'w'))
            shutil.rmtree(d)
        return [self.cache[a] for a in atoms]
M = Measurer()
QUAD = 10.95

def break_lines(s, first_extra=0.0):
    """returns list of (indent_depth, text)"""
    s = re.sub(r'\s+', ' ', s).strip()
    cands = candidates(s)
    cuts = [0] + [p for p, k, d in cands] + [len(s)]
    kinds = [0] + [k for p, k, d in cands]; depths = [0] + [d for p, k, d in cands]
    atoms = [s[cuts[i]:cuts[i+1]].strip() for i in range(len(cuts) - 1)]
    ws = M.measure(atoms); TW = M.cache['__tw__'] - 6.0
    ws = [w + (5.0 if kinds[i] == 1 else 1.5) for i, w in enumerate(ws)]
    lines = []; i = 0; n = len(atoms); times_next = False
    while i < n:
        depth = min(depths[i], 3) if lines else 0
        lim = TW - (first_extra if not lines else (1 + depth) * QUAD) - (14.0 if times_next else 0.0)
        w = ws[i]; j = i + 1; opts = [(i + 1, kinds[i + 1], depths[i + 1], w)] if i + 1 < n else []
        while j < n and w + ws[j] <= lim:
            w += ws[j]; j += 1
            if j < n: opts.append((j, kinds[j], depths[j], w))
        if j < n:
            prim = [o for o in opts if o[1] == 1 and o[3] >= 0.5 * lim]
            if prim:
                dmin = min(o[2] for o in prim); j = max(o[0] for o in prim if o[2] == dmin)
            # else: fallback break at the farthest candidate (inside a product)
        if ws[i] > lim and j == i + 1: print('  WIDE atom %.0fpt > %.0f: %s' % (ws[i], lim, atoms[i][:80]))
        t = ' '.join(atoms[i:j]).strip()
        t = re.sub(r'(\\,\s*)+$', '', t)
        if times_next: t = r'\times ' + t
        times_next = (j < n and kinds[j] == 2)
        lines.append((depth, t)); i = j
    return lines

def write_eq(fname, label, lhs, rhs, tail='', eqsign='='):
    lines = break_lines(lhs + eqsign + rhs)
    txt = [r'\begin{align}\label{%s}' % label]
    for k, (d, t) in enumerate(lines):
        ind = '' if k == 0 else r'\quad' * (1 + d) + ' '
        end = (tail if k == len(lines) - 1 else r'\notag\\')
        txt.append('&' + ind + t + end)
    txt.append(r'\end{align}')
    open(os.path.join(OUT, fname), 'w').write('\n'.join(txt) + '\n')

def paper_macros(s):
    return s.replace(r'\mathcal{K}', r'\Kv')

report = []
for lab, L, R in E1 + E2:
    R2, used = substitute(lab, cosmetic(R))
    # round trip
    assert cnorm(expand(R2)) == cnorm(SRC_RHS[lab]), lab
    need = {'I.II': 1, 'I.III': 2, 'I.V': 1, 'II.I': 2, 'III.II': 1, 'III.III': 2, 'III.V': 2, 'III.VI': 2}
    kind = lab.split('-')[-1]
    if kind in ('3', '2'): assert len(used) == need[row_of(lab)] and all(n == 1 for _, n in used), (lab, used)
    if re.match(r'2[a-d]$', kind) or (kind == '2o' and row_of(lab) == 'I.V'): assert len(used) == 1, (lab, used)
    fname = lab[3:].replace('.', '_') + '.tex'
    print(lab, used)
    write_eq(fname, lab, new_lhs(lab), paper_macros(R2))
    report.append((lab, L, new_lhs(lab), [u[0] for u in used]))

# kernel definitions
knames = list(KER)
for name in knames:
    c = KER[name]
    if c.startswith('+'): c = c[1:]
    fn = 'K-' + re.sub(r'[^A-Za-z0-9]+', '_', name.replace(r'\mathfrak{K}', '')).strip('_') + '.tex'
    write_eq(fn, 'eq:' + fn[:-4], name, paper_macros(c), eqsign=r'\equiv ')
    report.append(('eq:' + fn[:-4], name, name, []))
json.dump({'kernels': KER, 'delims': KDELIM, 'report': report}, open(os.path.join(HERE, 'gen_rows.json'), 'w'), indent=1)
print('equations written:', len(report))
