# Summary table of the power counting (scan_allclusters.json, scan_nlocal.json) -> LaTeX rows.
import json, collections
A = json.load(open('scan_allclusters.json')); N = json.load(open('scan_nlocal.json'))
def key(e): return e['unit']
def tex_pt(p): return {'x': 'x', 'y': 'y', 'z': 'z', 'w': 'w', "x'": "x'", "w'": "w'", "y'": "y'"}[p]
NUMERIC_ARTEFACT = {("I.III Kz |y->x", ('z', 'w', "w'")), ("I.III Kz |x'->y", ('z', 'w'))}
rows = collections.OrderedDict()
for e in A:
    k = key(e); S = tuple(e['S'])
    a = e['a_col'] if (k, S) not in NUMERIC_ARTEFACT else e['a_col_prev']
    a0 = e['a_nocol'] if (k, S) not in NUMERIC_ARTEFACT else e['a_col_prev']
    rows.setdefault(k, []).append((S, e['d'], a, a0, e['dist']))
nloc = collections.defaultdict(list)
for e in N:
    nloc[key(e)].append((tuple(e['S']), e['d'], e['a_col'], e['a_nocol']))
out = []
for k, lst in rows.items():
    if 'alone' in k: continue
    drop, keep = k.split('|')[1].split('->')
    with_p = [r for r in lst if keep in r[0]]
    worst = min(with_p, key=lambda r: r[2] + r[1])
    worst0 = min(with_p, key=lambda r: r[3] + r[1])
    need = [r for r in with_p if r[3] + r[1] < 0.2]          # kernel alone not integrable
    allok = all(r[2] + r[1] > 0.2 for r in lst)
    nl = nloc[k]; nl_ok = all(r[2] + r[1] > 0.2 for r in nl); nl_need = [r for r in nl if r[3] + r[1] < 0.2]
    cl = lambda S: '\\{' + ','.join(S) + '\\}'
    reason = ('kernel' if not need else 'colour $=0$ at ' + ', '.join(cl(r[0]) for r in need))
    out.append((k, len(with_p), cl(worst[0]), worst[1], worst[2], cl(worst0[0]) + ', ' + str(worst0[1]), worst0[3], allok, reason, nl_ok,
                ('' if not nl_need else 'colour at ' + ', '.join(cl(r[0]) for r in nl_need))))
with open('table_rows.tex', 'w') as f:
    for (k, n, S, d, a, S0, a0, ok, reason, nlok, nlr) in out:
        row, rest = k.split(' ', 1); name, con = rest.split(' |')
        name = {'Kz': r'$\mathfrak K^{\rm III}_z$', 'Kw': r'$\mathfrak K^{\rm III}_w$', 'KV1': r'$\mathfrak K^{\rm V}_1$',
                "KT(w,w')": r"$\mathfrak K_{\mathcal T}(w,w')$", "KT(w',w)": r"$\mathfrak K_{\mathcal T}(w',w)$",
                'Ehat': r'$\mathfrak K_{\widehat{\mathcal E}}$', 'Psi': r'$\mathfrak K_\Psi$',
                'delta': r"$\mathfrak K_\delta^{(y')}{+}\mathfrak K_\delta^{(z)}$, $\delta$", 'deltaPV': r"$\mathfrak K_\delta^{(y')}{+}\mathfrak K_\delta^{(z)}$, $\mathcal P$"}[name]
        d_, k_ = con.split('->')
        f.write('%s & %s & $%s\\to %s$ & $%s$, $%d$, $%.1f$ & $%s$, $%.1f$ & %s & %s & %s\\\\\n' % (
            row, name, d_, k_, S, d, a, S0, a0, 'yes' if ok else 'NO', reason.replace('->', r'$\to$'),
            'yes' if nlok else 'NO'))
print(open('table_rows.tex').read())
