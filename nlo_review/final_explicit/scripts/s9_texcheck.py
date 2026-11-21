# s9_texcheck.py OLD.tex NEW.tex -- TeX-level checks of the plus-prescription rewrite (2026-11-21):
# no closed-form function left in the result, colour markings and NOTEs kept, the exactly cancelling sets
# and the parts without special functions unchanged.  OLD = the version before the rewrite
# (git show 311b225:nlo_review/final_explicit/Final_terms_explicit.tex).  Run: python3 -I s9_texcheck.py OLD NEW
import sys, re
old, new = open(sys.argv[1]).read(), open(sys.argv[2]).read()
def body(s):  # compiled part only
    return s[:re.search(r'^\\end\{document\}', s, re.M).start()]
o, n = body(old), body(new)
res = []
def rep(name, ok, info=''):
    res.append((name, ok, info)); print('%-4s %-70s %s' % ('ok' if ok else 'FAIL', name, info))
# 1. special functions: allowed only in the quoted request (top note) and in the 'Instantaneous terms' / 'equivalently' remarks of Sec. conv
FUN = [r'\mathcal{E}', r'\mathcal{J}', r'\widehat{\mathcal{E}}', r'\Psi', r'\pi\delta(s)', r'\mathcal{P}\frac', r'\mathcal{P}\tfrac', r'{\rm Ci}', r'{\rm Si}', r'E_{1}']
lines = n.split('\n')
allowed = []
for i, l in enumerate(lines):
    if any(f in l for f in FUN):
        allowed.append(i+1)
conv_a = n.index(r'\section{Conventions of this file}'); conv_b = n.index(r'\section{Check of the author')
bad = []
for f in FUN:
    for m in re.finditer(re.escape(f), n):
        p = m.start()
        in_quote = 'Rewrite all the' in n[max(0, p-120):p]
        in_conv = conv_a < p < conv_b and ('Instantaneous terms' in n[n.rfind(r'\item', 0, p):p] or 'Equivalently, in closed form' in n[max(0, p-1500):p])
        if not (in_quote or in_conv): bad.append((f, n.count('\n', 0, p)+1))
rep('no closed-form function outside the allowed remarks', not bad, str(bad[:5]))
for f in FUN: print('     %-24s old %3d  new %3d' % (f, o.count(f), n.count(f)))
# 2. colour markings and NOTEs
for c in ['cI4', 'cUV', 'pA', 'pB', 'pC', 'pD', 'cPr', 'cJt', 'gray', 'red', 'orange']:
    a, b = o.count(r'\color{%s}' % c), n.count(r'\color{%s}' % c)
    rep(r'\color{%s} markings kept (more where a coloured line was split)' % c, b >= a and (b == a or c in ['pA', 'pB', 'pC', 'pD', 'red']), 'old %d new %d' % (a, b))
# every soft-remainder plus prescription in Secs. C3, D, E carries its pair (or blue) colour
SOFT = [r'\Big[\frac{1}{p^{+}}\Big]_{+}', r'\Big[\frac{k^{+}}{p^{+}}\Big]_{+}', r'\Big[\frac{k^{+}}{k^{+}-p^{+}}\Big]_{+}', r'\Big[\frac{k^{+}}{p^{+}\,(p^{+}+k^{+})}\Big]_{+}']
a_, b_ = n.index(r'\subsection{The soft remainders'), n.index(r'\section{Checks}')
unc = [i+1 for i, l in enumerate(n[a_:b_].split('\n')) if any(t in l for t in SOFT) and r'\color{' not in l and not l.startswith(r'\kappa') and 'kappa' not in l and r'\text{' not in l]
# lines inside \begingroup\color{..} blocks are coloured by the group
blk = n[a_:b_].split('\n'); ingroup = False; unc2 = []
for i, l in enumerate(blk):
    if l.startswith(r'\begingroup\color'): ingroup = True
    if l.startswith(r'\endgroup'): ingroup = False
    if (i+1) in unc and not ingroup: unc2.append(i+1)
rep('soft remainders in Secs. C3-E all coloured (pairs A-D / blue)', not unc2, str(unc2[:8]))
rep('NOTEs kept (one added on 2026-11-21)', n.count('NOTE') == o.count('NOTE')+1, 'old %d new %d' % (o.count('NOTE'), n.count('NOTE')))
rep('red phases: none removed', n.count(r'\ph{') >= o.count(r'\ph{'), 'old %d new %d' % (o.count(r'\ph{'), n.count(r'\ph{')))
rep('no phase in the author\'s convention added', len(re.findall(r"e\^\{-ik\\cdot\((w'-w|y'-w|w'-z|y'-y)\)", n)) == len(re.findall(r"e\^\{-ik\\cdot\((w'-w|y'-w|w'-z|y'-y)\)", o)))
# 3. untouched parts
def sec(s, a, b): return s[s.index(a):s.index(b)]
for nm, a, b in [('Secs. 1-4 (BFKL, JIMWLK, DGLAP, RC)', r'\section{BFKL contribution}', r'\section{Finite terms}'),
                 ('C1 Group I Row IV / Group III Row I (exact)', r'\subsection{Group I Row IV and Group III Row I: exact}', r'\subsection{The $\frac{11}{3}$'),
                 ('C2 11/3 UV terms (exact)', r'\subsection{The $\frac{11}{3}$', r'\subsection{The soft remainders'),
                 ('C5-C7 (in-row, joint tail, Coulomb log: exact)', r'\subsection{Cancellations inside single rows', r'\section{What remains'),
                 ('Delta N of Group I Row III (Eq. B-DN)', r'\frac{d^{3}\Delta N^{\rm I.III}}{d^{3}k}', r'\subsubsection*{Row IV}'),
                 ('Group I Row II, Row IV; Group II; Group III Rows I, II (Sec. B)', r'\subsubsection*{Row IV}', r'\subsubsection*{Row V}')]:
    rep('unchanged: '+nm, sec(o, a, b) == sec(n, a, b))
# 4. legend still lists all sets, labels unchanged
lo, ln = set(re.findall(r'\\label\{([^}]*)\}', o)), set(re.findall(r'\\label\{([^}]*)\}', n))
rep('equation labels: none removed except the special-function definitions', lo-ln <= {'eq:Ecal', 'eq:Jcal', 'eq:Ehat', 'eq:Psi', 'eq:distr'}, 'removed %s added %s' % (sorted(lo-ln), sorted(ln-lo)))
rep('no undefined \\eqref targets', all(r in ln for r in re.findall(r'\\eqref\{([^}]*)\}', n)))
print('all ok' if all(r[1] for r in res) else 'FAILURES')
open('s9_texcheck.log', 'w').write('\n'.join('%-4s %s %s' % ('ok' if ok else 'FAIL', nm, info) for nm, ok, info in res)+'\n')
