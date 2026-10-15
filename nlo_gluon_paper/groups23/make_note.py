"""Groups_II_III_check.tex from check_static.tex and checks23.json."""
import json


def sci(x):
    if x == 0:
        return '0'
    m, e = ('%.1e' % x).split('e')
    e = int(e)
    return m if e == 0 else r'%s\times10^{%d}' % (m, e)


ch = json.load(open('checks23.json'))
names = {'G2-I(1) A': 'II-I part 1, block A II(1.1)', 'G2-I(1) B': 'II-I part 1, block B II(1.1)', 'G2-I(2)': 'II-I part 2 II(1.3)',
         'G3-I': 'III-I III(1.1)', 'G3-II': 'III-II III(1.3)', 'G3-III(1)': 'III-III part 1 III(1.5)', 'G3-III(2)': 'III-III part 2 III(1.8)',
         'G3-V': 'III-V III(1.11)', 'G3-VI': 'III-VI III(1.15)', 'G2-IV': 'II-IV II(1.12)'}
rows = [r'\multicolumn{2}{@{}l}{\textbf{Before the integrations: written / first principles} (3 random configurations each)}\\']
for key, d in ch['pre'].items():
    ws = sorted(set(round(v[0], 6) for v in d['written_sym']))
    wg = [v[0] for v in d['written_gen']]
    cg = max(abs(v[0] - 1) + abs(v[1]) for v in d['corrected_gen'])
    rows.append(r'%s & $U=U^T$: %s; general $U$: %s; with Sec.~\ref{sec:M4} placements: $1\pm%s$\\' % (
        names[key], ', '.join('$%g$' % v for v in ws), ', '.join('%.2f' % v for v in wg[:3]), sci(cg)))
rows.append(r'\midrule\multicolumn{2}{@{}l}{\textbf{Integrated lines: |written/exact $-1$|} (general $U$, Sec.~\ref{sec:M4} placements)}\\')
for key, v in ch['integrated'].items():
    key2 = (key.replace('G2-I(1) line 50', 'II-I part 1, $p^+$ integrated II(1.2)').replace('G2-I(2) ', 'II-I part 2, ')
            .replace('G3-I line 48 (G3-II line 59 has the same kernel)', 'III-I, $p^+$ integrated III(1.2); III(1.4) has the same kernel')
            .replace('line 62 kernel / line 54 kernel', 'kernel of II(1.5) / kernel of II(1.3)')
            .replace('colour 58 / colour 54 (corrected placements)', 'colour of II(1.4) / colour of II(1.3)'))
    rows.append(r'%s & $%s$\\' % (key2, sci(v)))
rows.append(r'\midrule\multicolumn{2}{@{}l}{\textbf{$\mathfrak C$ rows of Group III: |written/exact $-1$|} (min--max over 3 configurations)}\\')
for key, (a, b) in ch['crows'].items():
    key2 = (key.replace("w' integration ", "$w'$ integration, ").replace('plus prescription ', 'plus prescription, ')
            .replace('III(1) line 72', 'III-III(1) III(1.6)').replace('III(2) line 85', 'III-III(2) III(1.9)')
            .replace('V line 105', 'III-V III(1.12)').replace('VI line 124', 'III-VI III(1.16)')
            .replace('III(1) line 76', 'III-III(1) III(1.7)').replace('III(2) line 90', 'III-III(2) III(1.10)')
            .replace('V line 109', 'III-V III(1.13)').replace('VI line 128', 'III-VI III(1.17)')
            .replace('written (from 0)', r'as written ($\int_0$)').replace('from Lambda', r'lower limit $\Lambda$'))
    rows.append(r'%s & $%s$ -- $%s$\\' % (key2, sci(a), sci(b)))
rows.append(r'\midrule\multicolumn{2}{@{}l}{\textbf{Row V simplification and Group II, Row IV}}\\')
for key, v in ch['row4'].items():
    key2 = (key.replace('Row V: colour of line 113 / line 109', r'III-V: colour of III(1.14) / III(1.13), $|{\rm ratio}-1|$')
            .replace("Row IV: line 98 / (line 92 at w*, w'* times (P/k)^4)", r"II-IV: II(1.13) / [II(1.12) at $w^*,w'^*$ $\times(\frac{p^++k^+}{k^+})^4$], $|{\rm ratio}-1|$")
            .replace('Row IV: colour of line 104 / line 98', r'II-IV: colour of II(1.14) / II(1.13), $|{\rm ratio}-1|$'))
    rows.append(r'%s & $%s$\\' % (key2, sci(v)))
rows.append(r"$J^{mm'}$: traceless part in 2 dims vs.\ $-\pi(\hat s^m\hat s^{m'}-\frac12\delta^{mm'})$; symmetry $|J^{12}-J^{21}|$ & $%s$; $%s$\\" % (sci(ch['J_traceless_maxdiff']), sci(ch['J_symmetric'])))
rows.append(r"II-IV: polarization sum $\times J^{mm'}$ ($d_\perp=2-\epsilon$) vs.\ the corrected bracket (Sec.~\ref{sec:M3}), relative & $%s$\\" % sci(ch['row4_assembly_reldiff']))
rows.append(r"$k^+(\frac{p^+}{(p^++k^+)^2}+\frac{p^+}{k^{+2}}+\frac1{p^+})-\big[\xi(1-\xi)+\frac{\xi}{1-\xi}+\frac{1-\xi}{\xi}\big]$; endpoint logarithm; LO normalization of the pole (\texttt{sympy}) & %s; %s; %s\\" % (
    ch['row4_phihat_identity'], ch['row4_endpoint_log'], ch['row4_LO_normalisation']))
T = open('check_static.tex').read().replace('@@CHECKS@@', '\n'.join(rows))
assert '@@' not in T
open('Groups_II_III_check.tex', 'w').write(T)
print('ok')
