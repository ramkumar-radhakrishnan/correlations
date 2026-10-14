"""Write RowV_UV.tex from RowV_static.tex and rowV_checks.json."""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def sci(x):
    if x == 0:
        return '0'
    m, e = ('%.1e' % x).split('e')
    e = int(e)
    return m if e == 0 else r'%s\times10^{%d}' % (m, e)


ch = json.load(open(os.path.join(HERE, 'rowV_checks.json')))
rows = []
rows.append(r'contraction $\mf C_1$ vertex $\times$ $\mf B_2$ bracket vs.\ your bracket (V.2), $d_\perp=2,3$, random points & max relative difference $%s$\\' % sci(ch['contraction_rel_diff']))
rows.append(r'prefactor $2|\mf A||\mf C_1||\mf B_2|$ / $\frac{g^4}{8\pi^4}\frac{p^+(k^+-p^+)}{k^{+3}D}$ (\texttt{sympy}) & %s\\' % ch['prefactor_ratio'])
rows.append(r'$f^{adc}f^{fbe}U^{db}U^{ce}=N_cU^{af}$ (same point; general SU(3) adjoint $U$) & $%s$\\' % sci(ch['ffUU_same_point']))
rows.append(r'$f^{abc}f^{cbd}=-N_c\delta^{ad}$; $f^{ebf}=-f^{fbe}$ & $%s$; $%s$\\' % (sci(ch['ff_contraction']), sci(ch['f_antisym'])))
for n in (2, 3):
    last = ch['small_r_d%d' % n]['rows'][-1]
    co = ', '.join('%.6f' % c for c in last['coef'])
    ex = ', '.join('%.4f' % c for c in ch['small_r_d%d' % n]['expected'])
    rows.append(r'$r^2\langle-\frac{r^m}{r^2}B_k^{im}\rangle\cdot Y/Y^2$, $k=1..6$, $d=%d$, $r=10^{-4}$ & %s (expected %s)\\' % (n, co, ex))
fk = [r for r in ch['full_kernel_small_r'] if r['r'] == 1e-4]
dev = max(abs(v - r['Phat']) / r['Phat'] for r in fk for v in r['ratio'])
rows.append(r'full kernel: $r^2\langle\VV\rangle\,Y^2/Y$ vs.\ $\Ph(\xi)$, $\xi=0.1,0.37,0.8$, $r=10^{-4}$ & max relative difference $%s$\\' % sci(dev))
o2 = ch['O2a_small_r']
rows.append(r'$\|\OO_{2a}\|$ at $r=10^{-1},10^{-2},10^{-3}$ (smooth random $U$) & %s (linear in $r$)\\' % ', '.join('%.4f' % r['norm'] for r in o2))
rows.append(r'$\int_\lambda^{1-\lambda}\Ph\,d\xi$; $\int_0^1c(\xi)d\xi$ & $%s$; $%s$\\' % (ch['int_Phat'].replace('log(Lambda)', r'\log\lambda').replace('*', ''), ch['int_plus_remainder']))
rows.append(r'$\int_0^1c(\xi)\log\xi\,d\xi$; $\int_0^1c(\xi)\log[\xi(1-\xi)]d\xi$ & $%s$; $%s$\\' % (ch['int_c_logxi'], ch['int_c_log_xi1mxi']))
rows.append(r'$d_\perp$-dim.\ integrals, Eq.~\eqref{eq:Iphase} and the cut version, minus the quoted $\overline{\rm MS}$ results (\texttt{sympy}) & %s; %s\\' % (ch['dimreg_phase_minus_target'], ch['dimreg_cut_minus_target']))
idm = max(abs(r['lhs'] - r['rhs']) for r in ch['phase_vs_cut_identity'])
rows.append(r'$\int d^2r\,[e^{iq\cdot r}-\theta(R-r)]/r^2=-\pi\log\frac{q^2R^2e^{2\gamma_E}}{4}$ (three $(q,R)$) & max difference $%s$\\' % sci(idm))
pu = ch['plus_remainder_uv']
rows.append(r'2b plus-prescribed remainder: $r^2\langle\cdot\rangle\,Y^2/Y$ at $r=10^{-1},10^{-2},10^{-3}$ & %s (UV coefficient $-\frac{11}{6}$)\\' % ', '.join('%.4f' % r['r2_avg_plus'][0] for r in pu))
rows.append(r'same, after subtracting the counterterm: $|r^2\langle\cdot\rangle|$ & %s (UV finite)\\' % ', '.join('$%s$' % sci(r['r2_avg_after_ct']) for r in pu))
pp = ch['plus_prescription_decomposition']
rows.append(r'Eq.~\eqref{eq:plus}: direct $\xi$ integral minus (plus integral + endpoint logs), $\Lambda/k^+=10^{-2},10^{-3},10^{-4}$ & %s (relative), $\propto\Lambda/k^+$\\' % ', '.join('$%s$' % sci(r['diff'] / r['size']) for r in pp))
T = open(os.path.join(HERE, 'RowV_static.tex')).read()
T = T.replace('@@CHECKS@@', '\n'.join(rows))
assert '@@' not in T
open(os.path.join(HERE, 'RowV_UV.tex'), 'w').write(T)
print('written RowV_UV.tex')
