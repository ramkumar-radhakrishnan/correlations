"""The symmetric three-rho part of every row of the direct set, written as kernel x Phi[p,q],
   Phi[p,q] = f^{abc} CL^a (U(p) rho(x))^b (U(q) rho(y))^c,   CL^a = [U(x') - U(w')]^{ab'} rho^{b'}(x'),
checked against the row itself (classical charges), row by row.  Units p8 = i g^4/(8 pi^5 k+) (times 1/(2pi)^3, phase 1)."""
import json
import numpy as np
import env  # noqa: F401
from fp import G, PI, sq, F
import ft
from ft import PI3
import split as SP
import three_rho_test as TT

p8 = 1j * G**4 / (8 * PI**5)
K = lambda v: v / sq(v)


def Phi(c, p, q):
    return np.einsum('abc,a,b,c->', F, ft.CL(c), c.U[p] @ c.rho['x'], c.U[q] @ c.rho['y'])


def kern(c):
    x, y, z, w = (c.pt[n] for n in ('x', 'y', 'z', 'w'))
    Xp, Y, Zw, Xw, Xz, Yw = K(c.pt['xp'] - c.pt['wp']), K(y - z), K(z - w), K(x - w), K(x - z), K(y - w)
    k = dict(Kb=(Xp @ Xw) * (Y @ Xz), Kh=(Xp @ Zw) * (Y @ Xz), D=-(Xp @ Yw) * (Xz @ Zw), E=(Xp @ Xw) * (Y @ Zw))
    k['KIV'] = k['Kh'] - k['Kb'] - k['E']
    k['Fz'] = TT.Fz(c)
    return k


FORMS = {
    'G1 II 2a': (TT.g1_ii_2a, lambda c, k: 0.5 * k['Fz'] * (Phi(c, 'x', 'y') - Phi(c, 'w', 'w'))),
    'G1 III (F3)': (TT.g1_iii_F3, lambda c, k: (k['Kh'] - k['Kb']) * (Phi(c, 'x', 'y') - Phi(c, 'x', 'z'))
                    + k['D'] * (Phi(c, 'x', 'y') - Phi(c, 'x', 'w'))),
    'K_IV family (G1 IV, G2 I(1) A, G3 I, G3 II)': (
        lambda c: TT.row_fn('G1-IV')(c) + TT.g2_i1_A(c) + TT.row_fn('G3-I')(c) + TT.row_fn('G3-II')(c),
        lambda c, k: k['KIV'] * (Phi(c, 'x', 'z') - Phi(c, 'w', 'w'))),
    'G1 V(1)': (lambda c: TT.g1_v1_t1(c) + TT.g1_v1_t2(c),
                lambda c, k: 0.5 * k['D'] * (Phi(c, 'z', 'w') - Phi(c, 'x', 'y')) + 0.5 * k['E'] * (Phi(c, 'w', 'z') - Phi(c, 'x', 'y'))),
    'K_C family c.c. (G3 III(1), V, VI)': (
        lambda c: TT.crow('G3-III(1)')(c) + TT.crow('G3-V')(c) + TT.crow('G3-VI')(c),
        lambda c, k: -k['D'] * (Phi(c, 'z', 'w') - Phi(c, 'z', 'y') - Phi(c, 'x', 'w') + Phi(c, 'x', 'y'))),
}

if __name__ == '__main__':
    out = {}
    for name, (row, form) in FORMS.items():
        rel = []
        for s in range(1, 6):
            c = TT.cfg(s)
            a = complex(row(c))
            b = complex(p8 * form(c, kern(c)) / PI3)
            rel.append(abs(a - b) / abs(a))
        out[name] = max(rel)
        print('%-48s row vs kernel x Phi form: %.1e' % (name, max(rel)))
    json.dump(out, open('phi_table.json', 'w'), indent=1)
