"""Build Evolution_singlet.tex from singlet_static.tex and the data of render_singlet.py, compare_kernels.py,
compare_rows.py, compare_rows_fast.py, validate.py, orient.py."""
import json
import os
import pickle
import re
from collections import OrderedDict
import env  # noqa: F401
from core import Config
import target as TG
from render_singlet import collect, numeric_value
from render_tex import block_tex_f, kernel_tex, check_factoring
from compare_rows_fast import merge

HERE = os.path.dirname(os.path.abspath(__file__))
J = lambda f: json.load(open(os.path.join(HERE, f)))
P = lambda f: pickle.load(open(os.path.join(HERE, f), 'rb'))


def e(v):
    """1.2e-15 -> 1.2\\times10^{-15}"""
    if v == 0:
        return '$0$'
    m, x = ('%.1e' % v).split('e')
    return r'$%s\times10^{%d}$' % (m, int(x))


def kkey_from_str(s):
    return eval(s)


ROWNAME = OrderedDict([
    ('G1-I_Nafter', ('Group I, Row I(a)', 'G1-I', 'eq:G1Ia')),
    ('G1-I_Nbarbefore', ('Group I, Row I(b)', 'G1-I', 'eq:G1Ib')),
    ('G1-Ia', ("Group I, Row I$'$", "G1-I'", 'eq:G1Ip')),
    ('G1-II_A3', ('Group I, Row II(a)', 'G1-II', 'eq:G1IIa')),
    ('G1-II_A3bar', ('Group I, Row II(b)', 'G1-II', 'eq:G1IIb')),
    ('G1-III', ('Group I, Row III', 'G1-III', 'eq:G1III')),
    ('G1-IV', ('Group I, Row IV', 'G1-IV', 'eq:G1IV')),
    ('G1-V', ('Group I, Row V', 'G1-V', 'eq:G1V')),
    ('G2-I', ('Group II, Row I', 'G2-I', 'eq:G2I')),
    ('G2-II', ('Group II, Row II', 'G2-II', 'eq:G2II')),
    ('G2-III', ('Group II, Row III', 'G2-III', 'eq:G2III')),
    ('G2-IV', ('Group II, Row IV', 'G2-IV', 'eq:G2IV')),
    ('G3-I', ('Group III, Row I', 'G3-I', 'eq:G3I')),
    ('G3-II', ('Group III, Row II', 'G3-II', 'eq:G3II')),
    ('G3-III', ('Group III, Row III', 'G3-III', 'eq:G3III')),
    ('G3-IV', ('Group III, Row IV', 'G3-IV', 'eq:G3IV')),
    ('G3-V', ('Group III, Row V', 'G3-V', 'eq:G3V')),
    ('G3-VI', ('Group III, Row VI', 'G3-VI', 'eq:G3VI')),
])

# what your row contains (for Table tab:rows)
USERROW = {
    'G1-I': r'two-$\rho$ part of lines 1, 2 of (4rho)',
    'G1-II': r"genuine + from three $\rho$'s",
    'G1-III': r"from three $\rho$'s ($N_1$, $N_2$) + line 4 of (4rho)",
    'G1-IV': r"from three $\rho$'s ($N_1$, $N_2$) + line 3 of (4rho)",
    'G1-V': r"genuine + from three $\rho$'s",
    "G1-I'": r'--- (inside your Group I)',
    'G2-I': r"genuine + from three $\rho$'s + line 5 of (4rho)",
    'G2-II': r'line 6 of (4rho)',
    'G2-III': r'line 7 of (4rho)',
    'G2-IV': r'genuine',
    'G3-I': r"from three $\rho$'s ($N_1$, $N_2$) + line 8 of (4rho)",
    'G3-II': r"from three $\rho$'s ($N_1$, $N_2$) + line 9 of (4rho)",
    'G3-III': r"genuine + from three $\rho$'s",
    'G3-IV': r'line 10 of (4rho)',
    'G3-V': r"from three $\rho$'s",
    'G3-VI': r"from three $\rho$'s",
}
REFNAME = {'G1-I': 'I(a)+I(b)', 'G1-II': 'II(a)+II(b)', 'G1-III': 'III', 'G1-IV': 'IV', 'G1-V': 'V', "G1-I'": "I$'$",
           'G2-I': 'I', 'G2-II': 'II', 'G2-III': 'III', 'G2-IV': 'IV', 'G3-I': 'I', 'G3-II': 'II', 'G3-III': 'III',
           'G3-IV': 'IV', 'G3-V': 'V', 'G3-VI': 'VI'}


def ok(v):
    return r'\new{%s}' % e(v) if v < 1e-10 else e(v)


def rowtab():
    C = J('compare_rows_fast.json')
    lines = []
    grp = None
    for k, v in C.items():
        if k in ('Group I + G2-IV', 'all rows'):
            continue
        g = k.split('-')[0]
        if grp is not None and g != grp:
            lines.append(r'\midrule')
        grp = g
        gname = {'G1': 'Group I', 'G2': 'Group II', 'G3': 'Group III'}[g]
        lines.append(r'%s, Row %s & %s & %s & %s\\' % (gname, REFNAME[k], USERROW[k], ok(v['sym']), ok(v['general'])))
    lines.append(r'\midrule')
    v = C['Group I + G2-IV']
    lines.append(r"Group I (all, with I$'$) + Group II Row IV & the same, summed & %s & %s\\" % (ok(v['sym']), ok(v['general'])))
    v = C['all rows']
    lines.append(r'all rows & all your lists & %s & %s\\' % (ok(v['sym']), ok(v['general'])))
    return '\n'.join(lines)


def kertab():
    C = J('compare_kernels.json')
    lines = []
    for p in C['per_kernel']:
        kk = kkey_from_str(p['kernel'])
        lines.append(r'$%s$ & %s & %s & %s\\' % (kernel_tex(kk), 'yes' if p['in_JIMWLK'] else 'no', ok(p['sym']),
                                                 ok(p['general'])))
    return '\n'.join(lines)


def nt(body):
    """one equation number per block"""
    return body.replace('\\\\\n', '\\notag\\\\\n')


def env_align(body, label):
    return '\\begin{align}\n' + nt(body) + '\n\\label{%s}\n\\end{align}' % label


def rows_tex(ref):
    out = []
    grp = None
    for k, (name, short, lab) in ROWNAME.items():
        g = name.split(',')[0]
        if g != grp:
            out.append(r'\subsection{%s}' % (g + (' (with c.c.)' if g != 'Group II' else '')))
            grp = g
        coll = ref[k]
        assert check_factoring(coll), k
        out.append(r'\paragraph{%s.}' % name)
        out.append(env_align(block_tex_f(coll), lab))
    return '\n\n'.join(out)


def lists_tex(rows):
    out = []
    sec = None
    for n, ((s, name), coll) in enumerate(rows.items()):
        if s != sec:
            out.append(r'\subsection{\texorpdfstring{%s}{%s}}' % (s, s.replace('$\\rho$', 'rho')))
            sec = s
        assert check_factoring(coll), name
        out.append(r'\paragraph{%s.}' % name)
        if not coll:
            out.append('(zero after the projection)')
            continue
        out.append(env_align(block_tex_f(coll), 'eq:list%d' % n))
    return '\n\n'.join(out)


if __name__ == '__main__':
    S = P('singlet.pkl')
    RC = P('rows_collected.pkl')
    REFP = P('ref_rows.pkl')
    val = J('validate.json')
    chk = J('singlet_checks.json')
    ker = J('compare_kernels.json')
    crs = J('compare_rows.json')
    ori = J('orient.json')

    rep = {}
    rep['VALMAX'] = e(max(v['max_abs_diff'] for v in val.values()))
    rep['VALL'] = e(val['lists']['max_abs_diff']) + ' (abs.)'
    rep['VALT'] = '%s / %s (abs.)' % (e(val['JIMWLK target']['max_abs_diff']), e(val['BFKL target']['max_abs_diff']))
    rj = max(abs(r['J_collected'] - r['appF_corrected']) / abs(r['appF_corrected']) for r in chk['numeric'])
    rb = max(abs(r['B_collected'] - r['appE_times_minusNc']) / abs(r['appE_times_minusNc']) for r in chk['numeric'])
    rep['NUMJ'] = 'relative difference ' + e(rj)
    rep['NUMB'] = 'relative difference ' + e(rb)

    # reference rows summed vs JIMWLK
    allr = merge(*RC['ref'].values())
    Jg = S['J_gen']
    refj = {}
    for sym in (True, False):
        r = []
        for seed in (3, 7, 11):
            cfg = Config(seed, sym)
            a, b = numeric_value(allr, cfg), numeric_value(Jg, cfg)
            r.append(abs(a - b) / abs(b))
        refj[sym] = max(r)
    rep['REFJ'] = 'relative difference %s for $U=U^T$, %s for a general $U$' % (e(refj[True]), e(refj[False]))
    rep['LISTJ'] = '%s / %s' % (e(ker['total']['sym']), e(ker['total']['general']))
    rep['GENTOT'] = '%.1f (relative, at random points)' % ker['total']['general']
    rep['FOURCHK'] = e(crs['four_rho_lines_vs_your_lists']['sym'])
    r1 = ori['R1']['general']
    r2 = ori['R2']['general']
    rep['ORIENT'] = '%.2f--%.1f (relative); with all lines up to the next $f$: %.2f--%.1f' % (min(r1), max(r1),
                                                                                          min(r2), max(r2))

    rep['LO'] = nt(block_tex_f(collect(TG.lo())))
    rep['JIM'] = nt(block_tex_f(Jg))
    rep['BFKL'] = nt(block_tex_f(S['B_gen']))
    rep['KERTAB'] = kertab()
    rep['ROWTAB'] = rowtab()
    rep['ROWS'] = rows_tex(REFP)
    rep['LISTS'] = lists_tex(S['rows'])

    tex = open(os.path.join(HERE, 'singlet_static.tex')).read()
    for k, v in rep.items():
        tex = tex.replace('@@%s@@' % k, v)
    left = re.findall(r'@@[A-Z]+@@', tex)
    assert not left, left
    open(os.path.join(HERE, 'Evolution_singlet.tex'), 'w').write(tex)
    print('written; refj', refj)
