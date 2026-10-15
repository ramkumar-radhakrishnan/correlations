"""Writes Main_corrections.tex from main_static.tex and main_checks.json (run run_checks.py first)."""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, 'main_checks.json')))


def tt(s):
    return s.replace('_', r'\_').replace('^', r'\^{}').replace('->', r'$\to$').replace('&', r'\&')


def ratio(d):
    return r'%.10g (spread %.1e)' % (d['mean_real'], d['spread'])


rows = []
for name, d in R['written / first principles (one gluon -> one gluon)'].items():
    rows.append((r'%s / first principles' % tt(name), ratio(d)))
for name, d in R['B3 / wave function (8 pi^3 = %.8f)' % (8 * 3.141592653589793**3)].items():
    rows.append((r'$\mf B_3$ %s / wave function' % tt(name).replace('Theta', r'$\Theta$').replace('Lambda', r'$\Lambda$'), ratio(d)))
for name, d in R['written / first principles (vacuum, one gluon, g^3)'].items():
    rows.append((r'%s / first principles' % tt(name), ratio(d)))
for name, d in R['B2 (momentum) / wave function (8 pi^3)'].items():
    rows.append((r'$\mf B_2$ (momentum) %s / wave function' % name, ratio(d)))
rows.append((r'$\Psi^3$ commutator, $\xi$-integrand (dim.\ reg.) minus closed form of X3',
             r'%.1e (three $\xi$ values)' % R['Psi^3 one-rho xi-integrand (dim reg) minus closed form']['max_abs']))
for r in R['4+5onerho (dim reg) vs notes']:
    rows.append((r'\eqn{4+5onerho} at $\xi=%.1f$: pole, $\ln\mu$, finite' % r['xi'],
                 r'%.6g, %.6g, %.10g (notes %.6g, %.6g, %.10g)' % (r['pole'], r['ln_mu'], r['finite'], r['pole_notes'], r['ln_mu_notes'], r['finite_notes'])))
for r in R['tworho2 (as written) vs exact']:
    rows.append((r'\eqn{tworho2} minus exact (odd part), $\delta=%g$' % r['delta'],
                 r'%.2e (exact $\sim$ %.2e)' % (r['max_abs_diff'], r['max_abs_exact'])))
for name, e in R['exact forms'].items():
    extra = ''
    if 'notes_Lambda0_form' in e:
        extra = r'; notes ($\Lambda\to0$) %.8g' % e['notes_Lambda0_form']
    if 'max_component_diff' in e:
        extra += r'; max diff %.1e' % e['max_component_diff']
    if 'vs_notes_with_ln(Lambda/(k+-Lambda))' in e:
        extra += r'; notes with $\ln\frac{\Lambda}{k^+-\Lambda}$: diff %.1e' % e['vs_notes_with_ln(Lambda/(k+-Lambda))']
    rows.append((tt(name), r'exact %.12g, quadrature %.12g (rel.\ %.0e)%s' % (e['exact'], e['quadrature'], e['rel_diff'], extra)))

body = '\n'.join(r'%s & %s\\' % (a, b) for a, b in rows)
src = open(os.path.join(HERE, 'main_static.tex')).read().replace('@@CHECKS@@', body)
open(os.path.join(HERE, 'Main_corrections.tex'), 'w').write(src)
print('wrote Main_corrections.tex with %d check rows' % len(rows))
