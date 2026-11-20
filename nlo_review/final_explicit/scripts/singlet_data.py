# Printed colour-singlet brackets (transcribed from paper sec_finite.tex, eqs. sIIa ... oIIIVI),
# used both to generate the TeX of the averaged section and to check it against direct contraction.
# term = (coef, kind, args); kind in D, F, Q, G, Qd ; coef may contain 'Nc' as factor via ('Nc', c)
Nc = 'Nc'
def D(a, b, c=1, nc=False): return (c, 'D', (a, b), nc)
def F(a, b, cc, c=1): return (c, 'F', (a, b, cc), False)
def Q(a, b, cc, d, c=1): return (c, 'Q', (a, b, cc, d), False)
def G(a, b, cc, c=1): return (c, 'G', (a, b, cc), False)
def Qd(a, b, cc, e, c=1): return (c, 'Qd', (a, b, cc, e), False)
def ND(a, b, c=1): return (c, 'D', (a, b), True)   # N_c D
# contacts: (pair label 'p1p2' etc, from-point, to-point, printed coefficient (in units of the row's
# 'A-free' prefactor (1/(2pi)^3) g^4/(16 pi^5) 1/k+ /(N^2-1)), n-arguments, bracket)
# CUBIC rows: name -> dict(A, t, typ, pts, structures: [(label, contacts, 3rho coefficient, 3rho bracket, o-args)])
CUBIC = {
 'I.II 2b': dict(A=-1/(4*3.141592653589793), A_tex=r'-\frac{1}{4\pi}', t=2, typ='r{rr}', pts=('xp', 'x', 'y'),
   structs=[('I.II 2b', [
       ('12', 'xp', 'x', r'-\frac{i}{4\pi}', -1j/(4*3.141592653589793), ('x', 'y'), [ND('x', 'y', -1), ND('x', 'w'), ND('w', 'wp', -1), F('x', 'y', 'wp')]),
       ('13', 'xp', 'y', r'+\frac{i}{4\pi}', 1j/(4*3.141592653589793), ('y', 'x'), [ND('x', 'y', -1), ND('y', 'w'), ND('w', 'wp', -1), F('x', 'y', 'wp')])],
     r'+\frac{1}{2\pi}', 1/(2*3.141592653589793), [G('x', 'y', 'xp'), G('x', 'y', 'wp', -1)], ('xp', 'x', 'y'))]),
 'I.III': dict(A=2j, t=1, typ='ord', pts=('xp', 'x', 'y'), structs=[
     (lab, [
       ('12', 'xp', 'x', '-', -1, ('x', 'y'), [ND('x', 'y'), ND('x', u, -1), F('x', 'y', 'wp', -1), F('x', u, 'wp')]),
       ('23', 'y', 'x', '-', -1, ('x', 'xp'), [ND('x', 'xp'), ND('x', 'wp', -1), F('x', 'xp', u, -1), F('x', u, 'wp')]),
       ('13', 'xp', 'y', '+', 1, ('y', 'x'), [ND('x', 'y'), F('x', 'y', u, -1), F('x', 'y', 'wp', -1), F('x', u, 'wp')])],
      '+2i', 2j, [G('x', 'y', 'xp'), G('x', 'y', 'wp', -1), G('x', 'xp', u), G('x', u, 'wp')], ('xp', 'x', 'y'))
     for lab, u in [('I.III C_D (u=z)', 'z'), ('I.III C_A (u=w)', 'w')]]),
 'I.V first': dict(A=0.5j, t=2, typ='r{rr}', pts=('xp', 'y', 'yp'), structs=[('I.V first', [
       ('12', 'xp', 'y', r'-\frac12', -0.5, ('y', 'yp'), [ND('y', 'yp'), F('x', 'y', 'z', -1), F('x', 'z', 'wp'), F('y', 'yp', 'wp', -1)]),
       ('13', 'xp', 'yp', r'+\frac12', 0.5, ('yp', 'y'), [ND('y', 'yp'), F('x', 'yp', 'z', -1), F('x', 'z', 'wp'), F('y', 'yp', 'wp', -1)])],
     '-i', -1j, [G('x', 'xp', 'z'), G('x', 'z', 'wp'), G('y', 'xp', 'yp'), G('y', 'yp', 'wp')], ('xp', 'y', 'yp'))]),
 'II.I 1a': dict(A=-1j, t=2, typ='{rr}r', pts=('xp', 'y', 'x'), structs=[('II.I 1a', [
       ('23', 'y', 'x', '+', 1, ('x', 'xp'), [ND('x', 'xp'), ND('w', 'wp'), F('x', 'z', 'wp', -1), Q('x', 'z', 'w', 'xp', -1)]),
       ('13', 'xp', 'x', '-', -1, ('x', 'y'), [ND('x', 'y'), ND('w', 'wp'), F('x', 'z', 'wp', -1), Q('x', 'w', 'z', 'y', -1)])],
     '-2i', -2j, [G('x', 'z', 'wp'), G('x', 'y', 'xp', -1), Qd('z', 'y', 'xp', 'w', -1)], ('xp', 'y', 'x'))]),
 'II.I 1b': dict(A=1j, t=2, typ='r{rr}', pts=('x', 'xp', 'y'), structs=[('II.I 1b', [
       ('12', 'xp', 'x', '-', -1, ('x', 'y'), [ND('x', 'y'), ND('w', 'wp'), F('x', 'z', 'w', -1), Q('x', 'wp', 'z', 'y', -1)]),
       ('13', 'y', 'x', '+', 1, ('x', 'xp'), [ND('x', 'xp'), ND('w', 'wp'), F('x', 'z', 'w', -1), Q('x', 'z', 'wp', 'xp', -1)])],
     '+2i', 2j, [G('x', 'z', 'w'), G('x', 'y', 'xp', -1), Qd('z', 'y', 'xp', 'wp', -1)], ('x', 'xp', 'y'))]),
 'III.II': dict(A=-2j, t=1, typ='ord', pts=('y', 'xp', 'x'), structs=[('III.II', [
       ('12', 'xp', 'y', '+', 1, ('y', 'x'), [ND('x', 'y'), ND('y', 'w'), F('x', 'y', 'z', -1), F('y', 'z', 'w', -1)]),
       ('23', 'xp', 'x', '+', 1, ('x', 'y'), [ND('x', 'y'), ND('x', 'z', -1), ND('x', 'w'), Q('x', 'w', 'z', 'y', -1)]),
       ('13', 'y', 'x', '-', -1, ('x', 'xp'), [ND('x', 'xp'), ND('xp', 'w'), F('x', 'xp', 'z', -1), Q('x', 'z', 'w', 'xp', -1)])],
     '-2i', -2j, [G('x', 'y', 'xp'), G('x', 'xp', 'z'), Qd('z', 'y', 'xp', 'w')], ('y', 'xp', 'x'))]),
 'III.III first': dict(A=-1j, t=2, typ='r{rr}', pts=('xp', 'y', 'x'), structs=[
     ("III.III L(y')", [
       ('12', 'xp', 'y', '+', 1, ('y', 'x'), [ND('x', 'y', -1), F('x', 'y', 'yp'), F('y', 'z', 'w'), F('yp', 'z', 'w', -1)]),
       ('13', 'xp', 'x', '-', -1, ('x', 'y'), [ND('x', 'y', -1), F('x', 'y', 'yp'), F('x', 'z', 'w'), F('yp', 'z', 'w', -1)])],
      '-2i', -2j, [G('x', 'y', 'xp'), G('x', 'y', 'yp', -1), G('xp', 'z', 'w', -1), G('yp', 'z', 'w')], ('xp', 'y', 'x')),
     ('III.III L(z)', [
       ('12', 'xp', 'y', '+', 1, ('y', 'x'), [ND('x', 'y', -1), ND('z', 'w', -1), F('x', 'y', 'z'), F('y', 'z', 'w')]),
       ('13', 'xp', 'x', '-', -1, ('x', 'y'), [ND('x', 'y', -1), ND('z', 'w', -1), F('x', 'y', 'z'), F('x', 'z', 'w')])],
      '-2i', -2j, [G('x', 'y', 'xp'), G('x', 'y', 'z', -1), G('xp', 'z', 'w', -1)], ('xp', 'y', 'x'))]),
 'III.V': dict(A=-2j, t=1, typ='ord', pts=('xp', 'x', 'y'), structs=[
     ("III.V L(y')", [
       ('12', 'xp', 'x', '+', 1, ('x', 'y'), [ND('x', 'y', -1), ND('x', 'w'), F('x', 'y', 'yp'), F('x', 'yp', 'w', -1)]),
       ('23', 'y', 'x', '+', 1, ('x', 'xp'), [ND('x', 'xp', -1), ND('x', 'yp'), F('x', 'xp', 'w'), F('x', 'yp', 'w', -1)]),
       ('13', 'xp', 'y', '-', -1, ('y', 'x'), [ND('x', 'y', -1), F('x', 'y', 'yp'), F('x', 'y', 'w'), F('x', 'yp', 'w', -1)])],
      '-2i', -2j, [G('x', 'y', 'xp', -1), G('x', 'y', 'yp'), G('x', 'xp', 'w', -1), G('x', 'yp', 'w')], ('xp', 'x', 'y')),
     ('III.V L(z)', [
       ('12', 'xp', 'x', '+', 1, ('x', 'y'), [ND('x', 'y', -1), ND('x', 'w'), F('x', 'y', 'z'), F('x', 'z', 'w', -1)]),
       ('23', 'y', 'x', '+', 1, ('x', 'xp'), [ND('x', 'xp', -1), ND('x', 'z'), F('x', 'xp', 'w'), F('x', 'z', 'w', -1)]),
       ('13', 'xp', 'y', '-', -1, ('y', 'x'), [ND('x', 'y', -1), F('x', 'y', 'z'), F('x', 'y', 'w'), F('x', 'z', 'w', -1)])],
      '-2i', -2j, [G('x', 'y', 'xp', -1), G('x', 'y', 'z'), G('x', 'xp', 'w', -1), G('x', 'z', 'w')], ('xp', 'x', 'y'))]),
 'III.VI': dict(A=-2j, t=1, typ='ord', pts=('xp', 'y', 'x'), structs=[
     ("III.VI L(y')", [
       ('12', 'xp', 'y', '+', 1, ('y', 'x'), [ND('x', 'y'), ND('y', 'z', -1), F('x', 'y', 'yp', -1), F('y', 'yp', 'z')]),
       ('23', 'y', 'x', '+', 1, ('x', 'xp'), [ND('x', 'xp'), ND('x', 'yp', -1), F('x', 'xp', 'z', -1), F('x', 'yp', 'z')]),
       ('13', 'xp', 'x', '-', -1, ('x', 'y'), [ND('x', 'y'), F('x', 'y', 'yp', -1), F('x', 'y', 'z', -1), F('y', 'yp', 'z')])],
      '-2i', -2j, [G('x', 'y', 'xp', -1), G('x', 'y', 'yp'), G('y', 'xp', 'z'), G('y', 'yp', 'z', -1)], ('xp', 'y', 'x')),
     ('III.VI L(z)', [
       ('12', 'xp', 'y', '+', 1, ('y', 'x'), [ND('x', 'y'), F('x', 'y', 'z', -1)]),
       ('23', 'y', 'x', '+', 1, ('x', 'xp'), [ND('x', 'xp'), F('x', 'xp', 'z', -1)]),
       ('13', 'xp', 'x', '-', -1, ('x', 'y'), [ND('x', 'y'), ND('y', 'z'), F('x', 'y', 'z', -2)])],
      '-2i', -2j, [G('x', 'y', 'xp', -1), G('x', 'y', 'z'), G('y', 'xp', 'z')], ('xp', 'y', 'x'))]),
}
# QUADRATIC structures: name -> (charge points (p,q), bracket of C^{aa} (U at p,q separate))
QUAD = {
 'I.II N1': (('xp', 'x'), [D('x', 'xp', -1), D('x', 'wp'), D('xp', 'w'), D('w', 'wp', -1)]),
 'I.V N2a': (('xp', 'y'), [ND('x', 'xp'), ND('x', 'wp', -1), F('x', 'xp', 'z', -1), F('x', 'z', 'wp')]),
 'I.V Nc[U(y)-U(x)]': (('xp', 'y'), [D('x', 'xp', -1), D('x', 'wp'), D('y', 'xp'), D('y', 'wp', -1)]),   # times N_c outside
 'II.I second': (('xp', 'x'), [ND('x', 'xp'), ND('w', 'wp'), F('x', 'z', 'wp', -1), F('xp', 'z', 'w', -1)]),
 "III.III second L(y')": (('xp', 'x'), [ND('x', 'xp', -1), ND('x', 'yp'), F('xp', 'z', 'w'), F('yp', 'z', 'w', -1)]),
 'III.III second L(z)': (('xp', 'x'), [ND('x', 'xp', -1), ND('x', 'z'), ND('z', 'w', -1), F('xp', 'z', 'w')]),
 'pair C colour C_C': (('xp', 'x'), [ND('x', 'xp'), ND('x', 'wp', -1), F('xp', 'z', 'w', -1), F('wp', 'z', 'w')]),
}
