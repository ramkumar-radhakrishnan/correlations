"""Write Part2b_MSbar.tex from Part2b_static.tex and part2b_checks.json."""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def sci(x):
    if x == 0:
        return '0'
    m, e = ('%.1e' % x).split('e')
    e = int(e)
    return m if e == 0 else r'%s\times10^{%d}' % (m, e)


ch = json.load(open(os.path.join(HERE, 'part2b_checks.json')))
rows = []
dev = max(r['ratio'] for r in ch['short_distance'])
rows.append(r'short distance: $r^2\langle\text{integrand}\rangle\,\frac{(y-z)^2}{(y-z)^i}$ vs.\ $\xi(1-\xi)+\frac{\xi}{1-\xi}+\frac{1-\xi}{\xi}$; $d=2$ and $d=3$ ($d_\perp=d$); $\xi=0.1,0.37,0.8$; $r=10^{-4}$ & max relative difference $%s$\\' % sci(dev))
sv = ch['S_angular_average']
rows.append(r'Eq.~\eqref{eq:angle}: $r^2\langle K(z,z,x)-K(y,z,x)\rangle$ at $|x-z|/|y-z|=%s$ & %s\\' % (
    ', '.join('%.2f' % r['r_over_R'] for r in sv), ', '.join('%.6f' % (0.0 if abs(r['r2_avg']) < 1e-12 else r['r2_avg']) for r in sv)))
for xi in (0.02, 0.37, 0.95):
    rr = [r for r in ch['rest_small_r'] if r['xi'] == xi]
    rows.append(r'$\xi=%.2f$: $r^2|\langle{\rm 2b}\rangle|$ and $r^2|\langle\text{rest}\rangle|$ at $r=10^{-2},10^{-3},10^{-4}$ & %s; %s\\' % (
        xi, ', '.join('%.4f' % r['r2_orig'] for r in rr), ', '.join('$%s$' % sci(r['r2_rest']) for r in rr)))
for xi in (0.2, 0.7):
    cc = [r for r in ch['rest_x_integral'] if r['xi'] == xi]
    rows.append(r'$\xi=%.1f$: ${\rm Re}\int_{r_0<|x-z|<0.3}d^2x$ (component $i=1$), $r_0=10^{-2},\dots,10^{-5}$: rest; 2b & %s; %s\\' % (
        xi, ', '.join('%.6f' % r['rest'][0] for r in cc), '$' + ',\\ '.join('%.2f' % r['orig'] for r in cc) + '$'))
d5 = max(abs(r['orig_minus_log'] - r['rest'][0]) for r in ch['rest_x_integral'])
rows.append(r'same: 2b minus $\pi\log\frac{0.09}{r_0^2}\times$(coefficient of Eq.~\eqref{eq:limit}) minus rest & max $%s$\\' % sci(d5))
for xi in (0.2, 0.7):
    ll = [r for r in ch['large_r'] if r['xi'] == xi]
    rows.append(r'$\xi=%.1f$, $|x-z|=10,10^2,10^3$: $|x-z|^3\max|{\rm CT}|$; $|x-z|^2\max|{\rm 2b}|$; $|x-z|^2|\langle{\rm 2b}\rangle|$ without phase & %s; %s; %s\\' % (
        xi, ', '.join('%.3f' % r['r3_max_ct'] for r in ll), ', '.join('%.3f' % r['r2_max_orig'] for r in ll),
        ', '.join('%.3f' % r['r2_avg_nophase_F1'] for r in ll)))
for xi in (0.2, 0.7):
    tt = [r for r in ch['rest_tail'] if r['xi'] == xi]
    rows.append(r'$\xi=%.1f$: ${\rm Re}\int_{1.5<|x-z|<R_{\max}}d^2x$ rest (phase kept), $R_{\max}=25,50,100,200$ & %s\\' % (
        xi, '$' + ',\\ '.join('%.3f' % r['tail'][0] for r in tt) + '$'))
rows.append(r'$d_\perp=2-\epsilon$ integrals minus $\pi[\frac2\epsilon+2\gamma_E+\log\frac{R^2\mu^2}{4}]$: $K$ form (Eq.~\eqref{eq:dimreg}), cut, smooth form (\texttt{sympy}) & %s; %s; %s\\' % (
    ch['dimreg_K_minus_target'], ch['dimreg_theta_minus_target'], ch['dimreg_smooth_minus_target']))
rows.append(r'$K$-form formula at $d=3$ minus $\frac14\cdot\frac{4\pi}{R}$ (\texttt{sympy}) & %s\\' % ch['K_formula_d3_minus_known'])
rows.append(r'Eq.~\eqref{eq:xi}: exact integral minus the right-hand side at $\lambda=\frac13,\frac1{10},\frac1{1000}$ (\texttt{sympy}, 30 digits) & $%s$\\' % sci(ch['xi_integral_minus_claim']))
rows.append(r'prefactor of the UV term$/\frac{1}{(2\pi)^3}\frac{g^4N_c}{16\pi^4}\frac{1}{k^+}$; $/\big(\frac{\alpha_sN_c}{2\pi}\times$LO prefactor$\big)$ & %s; %s\\' % (ch['pref_UV_over_g4Nc_16pi4'], ch['pref_UV_over_LO']))
rows.append(r'$\frac{\alpha_sN_c}{2\pi}(-\frac{11}{6})+\frac{\alpha_s\beta_0}{4\pi}$; $\beta_0$ coefficient G1-II(1) + 2b & %s; %s\\' % (ch['beta0_part'], ch['beta0_G1II_plus_2b']))
T = open(os.path.join(HERE, 'Part2b_static.tex')).read()
T = T.replace('@@CHECKS@@', '\n'.join(rows))
assert '@@' not in T
open(os.path.join(HERE, 'Part2b_MSbar.tex'), 'w').write(T)
print('written Part2b_MSbar.tex')
