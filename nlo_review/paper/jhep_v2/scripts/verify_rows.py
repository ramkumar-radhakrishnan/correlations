#!/usr/bin/env python3
# verify_rows.py -- checks that the equations printed in the paper (../rows/*.tex) are those of Sec. 8 of
# nlo_review/final_explicit/Final_terms_explicit_2_fixed.tex:
#  (1) every row equation, with the named kernels expanded again, has the same normal form as the source equation
#      (normal form: spacing, row breaks, explicit \times and \, , redundant braces removed; \Kv = \mathcal{K});
#  (2) the kernel definitions, the symmetrisation identities and R^i likewise;
#  (3) a copy of the source file is rebuilt with every equation of Sec. 8 replaced by the paper's equation, and the
#      numerical check of the printed colour factors (final_explicit/scripts/s16e_final_checks.py: random SU(3), SU(4)
#      Wilson lines, TeX parsed) is run on it.
# Run: python3 -I scripts/verify_rows.py [--colour]
import os, re, sys, json, shutil, subprocess, tempfile
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from rowlib import *
ROWS = os.path.join(HERE, '..', 'rows')
J = json.load(open(os.path.join(HERE, 'gen_rows.json'))); KER = J['kernels']; KD = J['delims']
def cosmetic(s):
    s = re.sub(r'(\\Bigg\[|\\bigg\[)\s*\\times', r'\1', s)
    return s
def cnorm(s):
    s = norm(cosmetic(flat(s)[0])); return s.replace(r'\times', '').replace(r'\,', '')
def expand(R):
    for name in sorted(KER, key=len, reverse=True):
        o, c = KD[name]; R = R.replace(name, o + KER[name] + c)
    return R
def paper_eq(fn):
    t = open(os.path.join(ROWS, fn)).read()
    t = t.replace(r'\begingroup\small', '').replace(r'\endgroup', '').strip() + '\n'
    m = re.match(r'\\begin\{align\}\\label\{([^}]*)\}\n(.*)\\end\{align\}\s*$', t, flags=re.S); assert m, fn
    body, _ = flat(m.group(2).replace(r'\notag\\', ' ')); assert '\\\\' not in body; return m.group(1), body
# independent reading of the source
F1, F2 = load_sec8(); SRC = {}
for T, pre in ((F1, 'H'), (F2, 'S')):
    for a in aligns(T):
        for e in split_equations(a):
            s, lab = flat(e)
            if r'\equiv' in s and 'C^{abc}_{{\\rm I.V}}' in s: continue
            L, R = top_eq_split(s)
            if lab is None:
                lab = 'eq:%s-II.I-2%s' % (pre, re.search(r'II\.I\}_\{2([a-e])\}', L).group(1))
            SRC[lab] = (L, R)
assert len(SRC) == 50
bad = 0; n = 0; PAPER = {}
for lab, (L, R) in SRC.items():
    plab, body = paper_eq(lab[3:].replace('.', '_') + '.tex'); assert plab == lab
    pL, pR = top_eq_split(body); PAPER[lab] = pR
    row = re.match(r'eq:[HS]-([IV]+\.[IV]+)', lab).group(1)
    assert ('N^{\\rm %s}' % row) in pL and ('N^{\\rm %s}' % row) in L, (lab, pL)
    assert (r'\langle' in pL) == lab.startswith('eq:S'), lab
    ok = cnorm(expand(pR)) == cnorm(R); n += 1
    if not ok: bad += 1; print('DIFFERENT:', lab)
print('(1) row equations compared with the source: %d, different: %d' % (n, bad))
# kernels: definitions printed == kernel strings used in the expansion
nk = 0
for name, c in KER.items():
    fn = 'K-' + re.sub(r'[^A-Za-z0-9]+', '_', name.replace(r'\mathfrak{K}', '')).strip('_') + '.tex'
    plab, body = paper_eq(fn); i = body.index(r'\equiv'); assert cnorm(body[:i]) == cnorm(name)
    rhs = body[i + len(r'\equiv'):].rstrip(' .,').rstrip('\\,').rstrip()
    c2 = c[1:] if c.startswith('+') else c
    assert cnorm(rhs) == cnorm(c2), name; nk += 1
    # and the kernel string occurs in the source
    srcall = cnorm(' '.join(R for L, R in SRC.values())); assert cnorm(c) in srcall, name
doc = open(SRC_PATH if False else os.path.join(HERE, '..', '..', '..', 'final_explicit', 'Final_terms_explicit_2_fixed.tex')).read()
body = aligns(doc[doc.index(r'\label{sec:D}'):])[0].replace(r'\color{red}', '')
ids = split_equations(body)
for e, fn in zip(ids, ('symmetrise', 'symmetrise12', 'symmetrise21')):
    s, lab = flat(e); plab, pb = paper_eq(fn + '.tex')
    assert cnorm(pb.rstrip(' ,.').rstrip('\\,')) .rstrip(',.') == cnorm(s.rstrip(' ,.').rstrip('\\,')).rstrip(',.'), fn
m = re.search(r'where: \$\$(R\^\{i\}.*?)\$\$', doc[doc.index(r'\label{sec:F2}'):], flags=re.S)
plab, pb = paper_eq('Rdef.tex'); assert cnorm(pb.rstrip(' .,').rstrip('\\,')) == cnorm(m.group(1))
print('(2) kernel definitions: %d, symmetrisation identities: 3, R^i: 1 -- all equal to the source' % nk)
def canon(x):
    x = x.replace(r'\times', ' ').replace(r'\,', ' ')
    x = re.sub(r'\\rm\s+', r'\\rm@', x); x = re.sub(r'\s+', '', x).replace(r'\rm@', r'\rm ')
    return x
if '--colour' in sys.argv:
    # (3) rebuild the source file with the paper's equations and run the printed-colour check on it
    d = tempfile.mkdtemp(); top = os.path.join(d, 'final_explicit'); FE = os.path.join(HERE, '..', '..', '..', 'final_explicit')
    shutil.copytree(os.path.join(FE, 'scripts'), os.path.join(top, 'scripts'), ignore=shutil.ignore_patterns('__pycache__', '*.log'))
    shutil.copy(os.path.join(FE, 'Final_terms_explicit_2_author.tex'), top)
    new = doc; pos = {'H': new.index(r'\label{sec:F1}'), 'S': None}
    for lab, (L, R) in SRC.items():
        pe = canon(expand(PAPER[lab]).replace(r'\Kv', r'\mathcal{K}'))
        k = lab[3]
        if pos[k] is None: pos[k] = new.index(r'\label{sec:F2}')
        a = new.index(L, pos[k])                               # the equation, located by its left-hand side (in order)
        m = re.compile(r'\s*&?\s*=\s*&?\s*').match(new, a + len(L)); assert m, lab
        b = m.end()
        if lab[-1] in 'abcde' and 'II.I-2' in lab:
            if lab[-1] == 'e': end = new.index(r'\end{align}', b); tail = r'\notag' + '\n'
            else:
                e = re.compile(r'(?<!\\notag)\\\\').search(new, b)
                while new[:e.start()].rstrip().endswith(r'\notag'): e = re.compile(r'\\\\').search(new, e.end())
                end = e.start(); tail = ''
        else:
            end = new.index(r'\label{%s}' % lab, b); tail = ''
        new = new[:b] + pe + tail + new[end:]; pos[k] = b + len(pe)
    open(os.path.join(top, 'Final_terms_explicit_2_fixed.tex'), 'w').write(new)
    # s16e locates the coefficients with regular expressions that contain the thin spaces "\," of the source;
    # the paper's line breaking moves these, so they are removed from the rebuilt equations (canon) and from the patterns
    q = os.path.join(top, 'scripts', 's16e_final_checks.py'); t = open(q).read(); lines = t.split('\n'); npat = 0
    for i, l in enumerate(lines):
        if l.startswith('pat = re.compile(') or l.lstrip().startswith(('r"|\\\\Big', 'r"|\\\\frac', 'Ls = re.findall(', 'Rs = re.findall(')):
            lines[i] = l.replace('\\\\,', ''); npat += 1
    assert npat == 5, npat
    open(q, 'w').write('\n'.join(lines))
    r = subprocess.run([sys.executable, '-I', 's16e_final_checks.py'], cwd=os.path.join(top, 'scripts'), capture_output=True, text=True)
    print('(3) s16e_final_checks.py on the file rebuilt from the paper equations (exit %d):' % r.returncode)
    print('\n'.join('    ' + l for l in (r.stdout + r.stderr).strip().split('\n')[-14:]))
    shutil.rmtree(d)
