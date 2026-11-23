# s16_sec8data.py -- colour structures of the finite terms written with differences of Wilson lines (new Sec. 8 of
# Final_terms_explicit_2_fixed.tex).  Every structure is given twice: (i) the TeX printed in the rho-explicit form
# (index notation; charge indices a,b,c in the order of the row, Sec. 6), and (ii) a 'leg' form used to generate the
# colour-singlet traces: C^{abc} = sum_terms coef * f^{ijk} L1^{i q1} L2^{j q2} L3^{k q3}, where L1, L2, L3 are
# matrices (products of U(p), U^T(p), U(p)-U(q), U^T(p)-U^T(q); empty = unit matrix) and (q1,q2,q3) the charge slots.
# Tokens: ('U',p) ('Ut',p) ('D',p,q)=U(p)-U(q) ('Dt',p,q)=U^T(p)-U^T(q).   Checked in s16c_sec8_colour.py.
U = lambda p: ('U', p); Ut = lambda p: ('Ut', p); D = lambda p, q: ('D', p, q); Dt = lambda p, q: ('Dt', p, q)
TEXPT = {'xp': "x'", 'wp': "w'", 'yp': "y'", 'x': 'x', 'y': 'y', 'z': 'z', 'w': 'w'}
def Ud(i, p): return r'U^{%s}(%s)' % (i, TEXPT[p])
def Dd(i, p, q): return r'\big[U^{%s}(%s)-U^{%s}(%s)\big]' % (i, TEXPT[p], i, TEXPT[q])
def DD(i, p, q): return r'\Big[U^{%s}(%s)-U^{%s}(%s)\Big]' % (i, TEXPT[p], i, TEXPT[q])
# ---------------------------------------------------------------- cubic structures
CUBIC8 = {}
CUBIC8['I.II 2b'] = dict(
    tex=r"-" + DD('ga', 'xp', 'wp') + r"f^{gde}\bigg[" + DD('db', 'x', 'w') + Ud('ec', 'y') + r"+" + Ud('db', 'w') + DD('ec', 'y', 'w') + r"\bigg]",
    legs=[(-1, [[D('xp', 'wp')], [D('x', 'w')], [U('y')]], (1, 2, 3)),
          (-1, [[D('xp', 'wp')], [U('w')], [D('y', 'w')]], (1, 2, 3))])
for lab, u in [('I.III C_D (u=z)', 'z'), ('I.III C_A (u=w)', 'w')]:
    CUBIC8[lab] = dict(tex=DD('ga', 'xp', 'wp') + r"f^{ghd}\,U^{db}(x)" + DD('hc', u, 'y'),
                       legs=[(1, [[D('xp', 'wp')], [D(u, 'y')], [U('x')]], (1, 3, 2))])
CUBIC8['I.V first'] = dict(
    tex=DD('ga', 'xp', 'wp') + r"f^{gde}\bigg[" + DD('dc', 'x', 'yp') + Ud('eb', 'z') + r"+" + Ud('dc', 'yp') + DD('eb', 'z', 'y') + r"\bigg]",
    legs=[(1, [[D('xp', 'wp')], [D('x', 'yp')], [U('z')]], (1, 3, 2)),
          (1, [[D('xp', 'wp')], [U('yp')], [D('z', 'y')]], (1, 3, 2))])
CUBIC8['II.I 1a'] = dict(
    tex=DD('ga', 'wp', 'xp') + r"\bigg[" + DD('gk', 'w', 'z') + r"f^{kbc}-f^{gkl}\,U^{kb}(z)" + DD('lc', 'x', 'z') + r"\bigg]"
        + r"+U^{ga}(x')\bigg[U^{mh}(z)" + DD('mb', 'z', 'y') + DD('gk', 'w', 'z') + r"f^{khc}-f^{gml}" + DD('mb', 'z', 'y') + DD('lc', 'x', 'z') + r"\bigg]",
    legs=[(1, [[Dt('w', 'z'), D('wp', 'xp')], [], []], (1, 2, 3)),
          (-1, [[D('wp', 'xp')], [U('z')], [D('x', 'z')]], (1, 2, 3)),
          (1, [[Dt('w', 'z'), U('xp')], [Ut('z'), D('z', 'y')], []], (1, 2, 3)),
          (-1, [[U('xp')], [D('z', 'y')], [D('x', 'z')]], (1, 2, 3))])
CUBIC8['II.I 1b'] = dict(
    tex=DD('gb', 'w', 'xp') + r"\bigg[" + DD('gk', 'wp', 'z') + r"f^{kca}-f^{gkl}\,U^{kc}(z)" + DD('la', 'x', 'z') + r"\bigg]"
        + r"+U^{gb}(x')\bigg[U^{mh}(z)" + DD('mc', 'z', 'y') + DD('gk', 'wp', 'z') + r"f^{kha}-f^{gml}" + DD('mc', 'z', 'y') + DD('la', 'x', 'z') + r"\bigg]",
    legs=[(1, [[Dt('wp', 'z'), D('w', 'xp')], [], []], (2, 3, 1)),
          (-1, [[D('w', 'xp')], [U('z')], [D('x', 'z')]], (2, 3, 1)),
          (1, [[Dt('wp', 'z'), U('xp')], [Ut('z'), D('z', 'y')], []], (2, 3, 1)),
          (-1, [[U('xp')], [D('z', 'y')], [D('x', 'z')]], (2, 3, 1))])
CUBIC8['III.II'] = dict(
    tex=r"-U^{gb}(x')\,U^{mh}(z)" + DD('ma', 'z', 'y') + DD('gk', 'w', 'z') + r"f^{khc}+U^{gb}(x')\,f^{gml}" + DD('ma', 'z', 'y') + DD('lc', 'x', 'z'),
    legs=[(-1, [[Dt('w', 'z'), U('xp')], [Ut('z'), D('z', 'y')], []], (2, 1, 3)),
          (1, [[U('xp')], [D('z', 'y')], [D('x', 'z')]], (2, 1, 3))])
for v, tag in [('yp', "(y')"), ('z', '(z)')]:
    CUBIC8['III.III L' + tag] = dict(
        tex=r"f^{mnk}" + DD('ma', v, 'xp') + r"\bigg[" + DD('nc', 'z', 'x') + Ud('kb', 'w') + r"+" + Ud('nc', 'x') + DD('kb', 'w', 'y') + r"\bigg]",
        legs=[(1, [[D(v, 'xp')], [D('z', 'x')], [U('w')]], (1, 3, 2)),
              (1, [[D(v, 'xp')], [U('x')], [D('w', 'y')]], (1, 3, 2))])
    CUBIC8['III.V L' + tag] = dict(tex=r"f^{mkl}" + DD('ma', v, 'xp') + r"U^{kb}(x)" + DD('lc', 'y', 'w'),
                                   legs=[(1, [[D(v, 'xp')], [U('x')], [D('y', 'w')]], (1, 2, 3))])
    CUBIC8['III.VI L' + tag] = dict(tex=r"f^{mnk}" + DD('ma', v, 'xp') + r"U^{kb}(y)" + DD('nc', 'x', 'z'),
                                    legs=[(1, [[D(v, 'xp')], [D('x', 'z')], [U('y')]], (1, 3, 2))])
# ---------------------------------------------------------------- quadratic (genuine two-rho) colour factors
# tex: two-index colour factor (charge indices as in the row); sing: its delta-contraction (singlet bracket)
QUAD8 = {
 'I.II N1': dict(idx=("b'", 'b'),
    tex=DD("ab'", 'xp', 'wp') + DD('ab', 'w', 'x'),
    sing=r"{\rm Tr}\Big\{\big(U^{\dagger}(x')-U^{\dagger}(w')\big)\big(U(w)-U(x)\big)\Big\}"),
 'I.V N2': dict(idx=("b'", 'f'),
    tex=r"-" + DD("ab'", 'xp', 'wp') + r"f^{adc}f^{fbe}\,U^{db}(x)" + DD('ce', 'z', 'x') + r"+N_{c}" + DD("ab'", 'xp', 'wp') + DD('af', 'y', 'x'),
    sing=r"-{\rm Tr}\Big\{T^{e}U^{\dagger}(x)\big(U(x')-U(w')\big)T^{e}\big(U^{\dagger}(z)-U^{\dagger}(x)\big)U(x)\Big\}+N_{c}{\rm Tr}\Big\{\big(U^{\dagger}(x')-U^{\dagger}(w')\big)\big(U(y)-U(x)\big)\Big\}"),
 'I.V 11/6': dict(idx=("b'", 'f'),
    tex=DD("ab'", 'xp', 'wp') + DD('af', 'y', 'x'),
    sing=r"{\rm Tr}\Big\{\big(U^{\dagger}(x')-U^{\dagger}(w')\big)\big(U(y)-U(x)\big)\Big\}"),
 'II.I second': dict(idx=("d'", 'd'),
    tex=r"\bigg[" + DD("ac'", 'wp', 'z') + r"f^{c'bd'}-f^{ac'e'}\,U^{c'b}(z)" + DD("e'd'", 'xp', 'z') + r"\bigg]\bigg[" + DD('ac', 'w', 'z') + r"f^{cbd}-f^{ace}\,U^{cb}(z)" + DD('ed', 'x', 'z') + r"\bigg]",
    sing=r"N_{c}{\rm Tr}\Big\{\big(U^{\dagger}(w')-U^{\dagger}(z)\big)\big(U(w)-U(z)\big)\Big\}+N_{c}{\rm Tr}\Big\{\big(U^{\dagger}(x')-U^{\dagger}(z)\big)\big(U(x)-U(z)\big)\Big\}"
         r"-{\rm Tr}\Big\{T^{e}U^{\dagger}(z)\big(U(x)-U(z)\big)T^{e}\big(U^{\dagger}(w')-U^{\dagger}(z)\big)U(z)\Big\}-{\rm Tr}\Big\{T^{e}U^{\dagger}(z)\big(U(x')-U(z)\big)T^{e}\big(U^{\dagger}(w)-U^{\dagger}(z)\big)U(z)\Big\}"),
 'II.IV': dict(idx=("e'", 'e'),
    tex=DD("ce'", 'yp', 'xp') + DD('ce', 'y', 'x'),
    sing=r"{\rm Tr}\Big\{\big(U^{\dagger}(y')-U^{\dagger}(x')\big)\big(U(y)-U(x)\big)\Big\}"),
 'II.IV rem': dict(idx=("e'", 'e'),
    tex=DD("ce'", 'z', 'xp') + DD('ce', 'z', 'x'),
    sing=r"{\rm Tr}\Big\{\big(U^{\dagger}(z)-U^{\dagger}(x')\big)\big(U(z)-U(x)\big)\Big\}"),
}
for v, tag in [('yp', "(y')"), ('z', '(z)')]:
    QUAD8['III.III second L' + tag] = dict(idx=("e'", 'd'),
        tex=r"f^{c'd'a}\,U^{d'b}(z)" + DD("c'e'", v, 'xp') + r"\bigg[" + DD('ac', 'w', 'z') + r"f^{cbd}-f^{ace}\,U^{cb}(z)" + DD('ed', 'x', 'z') + r"\bigg]",
        sing=r"N_{c}{\rm Tr}\Big\{\big(U^{\dagger}(%s)-U^{\dagger}(x')\big)\big(U(x)-U(z)\big)\Big\}-{\rm Tr}\Big\{T^{e}U^{\dagger}(z)\big(U(w)-U(z)\big)T^{e}\big(U^{\dagger}(%s)-U^{\dagger}(x')\big)U(z)\Big\}" % (TEXPT[v], TEXPT[v]))
