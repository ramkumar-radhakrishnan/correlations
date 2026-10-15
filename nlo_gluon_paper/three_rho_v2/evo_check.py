"""Re-parse the equations of Evolution_from_final.tex and check them pointwise:
gen + 4rho + (two-rho of the three-rho equation) = JIMWLK x LO;  Plist + 4rho = BFKL x LO;  written J, B = J, B."""
import re
import json
import env
from core import term, symmetrize, Config, pointwise_total
import two_from_final as TF
import final_eq as FE
import target as TG
import refP

TEX = open('Evolution_from_final.tex').read()


def block(label):
    i = TEX.index(r'\label{%s}' % label)
    s = TEX.rindex(r'\begin{align}', 0, i)
    e = TEX.index(r'\end{align}', i)
    return TEX[s:e]


def kern_parse(ktx):
    kern = []
    ks = [tuple(s.split('-')) for s in TF.KRE.findall(ktx)]
    for i in range(0, len(ks), 2):
        kern.append(('KK', ks[i][0], ks[i][1], ks[i + 1][0], ks[i + 1][1]))
    for a in TF.INV.findall(ktx):
        kern.append(('KK', a, 'z', a, 'z'))
    return kern


def entries(label):
    out = []
    for ln in block(label).split('\n'):
        m = re.match(r'^&\\quad (.*?)\\;(.*?)\\;(.*?)(?:\\notag\\\\|\\Bigg\]\\label\{[^}]*\})$', ln)
        if not m:
            continue
        col = m.group(3)
        out.append(term(float(TF.coef_parse(m.group(1))), TF.URE.findall(col), TF.FRE.findall(col), TF.RRE.findall(col),
                        kern_parse(m.group(2)), label))
    return out


def target(label):
    out, kern = [], None
    for ln in block(label).split('\n'):
        m = re.match(r"^&\\quad \+\\underbrace\{(.*)\}_\{S'?_\{\d+\}\}\\Big\[\\notag\\\\$", ln)
        if m:
            kern = kern_parse(m.group(1))
            continue
        m = re.match(r'^&\\qquad\\qquad (.*?)\\;(.*?)(?:\\Big\])?(?:\\notag\\\\|\\Big\\\}?\\Bigg\]\\label\{[^}]*\}|\\Bigg\]\\label\{[^}]*\})$', ln)
        if m and kern is not None:
            col = m.group(2).replace(r'\Big]', '')
            out.append(term(float(TF.coef_parse(m.group(1))), TF.URE.findall(col), TF.FRE.findall(col), TF.RRE.findall(col), kern, label))
    return out


final = [r for r in FE.reparse('ThreeRho_final.tex') if r['color'] != 'red']
three = [s for r in final for t in FE.terms_of(r) for s in symmetrize(t, (2,))]
gen, four, plist = entries('eq:gen'), entries('eq:four'), entries('eq:Plist')
Jw, Bw = target('eq:Jexp'), target('eq:Bexp')
J, B = TG.jimwlk(), refP.bfkl()
print('parsed: gen %d, four %d, Plist %d, J written %d, B written %d' % (len(gen), len(four), len(plist), len(Jw), len(Bw)))
res, ok = [], True
for seed in (3, 7, 11, 19):
    cfg = Config(seed, True)
    v = lambda ts: pointwise_total(ts, cfg)[0].real
    a = v(three) + v(gen) + v(four)
    p = v(plist) + v(four)
    r = dict(seed=seed, sumA=a, J=v(J), Jwritten=v(Jw), sumP=p, B=v(B), Bwritten=v(Bw))
    res.append(r)
    good = all(abs(x - y) < 1e-9 * (1 + abs(y)) for x, y in ((a, r['J']), (r['Jwritten'], r['J']), (p, r['B']), (r['Bwritten'], r['B'])))
    ok &= good
    print('seed %2d: sum A %.8f  J %.8f  J(written) %.8f | sum P %.8f  B %.8f  B(written) %.8f  %s' % (
        seed, a, r['J'], r['Jwritten'], p, r['B'], r['Bwritten'], 'ok' if good else 'MISMATCH'))
print('ALL CHECKS PASS' if ok else 'CHECK FAILED')
json.dump(res, open('evo_check.json', 'w'), indent=1)
