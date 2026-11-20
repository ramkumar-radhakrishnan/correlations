# The (2rho) contact lines of sec_finite.tex, Eqs. (sIIII), (sIV1), (sIII1), (sIIIII), (sIIIIII1), (sIIIV), (sIIIVI):
# kernel (with the contact substitution drop -> keep) times colour factor (as printed).
import re
import numpy as np
import lib

def parse_colour(s):
    """'N D(x,y) - N D(x,z) - F(x,y,w') + F(x,z,w')' -> list of (coef, kind, args)"""
    s = s.replace(' ', '')
    terms = re.findall(r"([+-]?)(\d*)(N?)([DFQ])\(([^)]*)\)", s)
    out = []
    for sign, num, n, kind, args in terms:
        c = (-1.0 if sign == '-' else 1.0) * (float(num) if num else 1.0) * (lib.NC if n else 1.0)
        out.append((c, kind, [a for a in re.split('[,;]', args)]))
    # sanity: reconstruct
    return out

def colour_fun(s):
    terms = parse_colour(s)
    def f(col, P):
        v = 0.0
        for c, kind, args in terms:
            v += c * getattr(col, kind)(*[P[a] for a in args])
        return v
    f.terms = terms; f.text = s
    return f

G1 = ['x', 'y', 'z', 'w', "x'", "w'"]
G15 = ['x', 'z', 'y', "y'", "x'", "w'"]
G2 = ['x', 'y', 'z', 'w', "x'", "w'"]
G3 = ['w', "x'", "y'", 'z', 'y', 'x']
CHARGES = {'I.III': ["x'", 'x', 'y'], 'I.V': ["x'", 'y', "y'"], 'II.I': ["x'", 'y', 'x'], 'III.II': ['y', "x'", 'x'],
           'III.III': ["x'", 'y', 'x'], 'III.V': ["x'", 'x', 'y'], 'III.VI': ["x'", 'y', 'x']}

class Unit:
    def __init__(self, row, name, contact, base, kern, colour, coef=1.0, dist=None, kern2=None, colour2=None,
                 fourier=(), pair=None):
        self.row, self.name, self.coef, self.dist = row, name, coef, dist
        self.drop, self.keep = contact
        self.base = base
        self.points = [p for p in base if p != self.drop]
        self.kern, self.colour = kern, colour_fun(colour)
        self.kern2 = kern2; self.colour2 = colour_fun(colour2) if colour2 else None
        self.fourier = fourier; self.pair = pair
        ch = CHARGES[row]
        self.remaining_charges = [self.keep] + [c for c in ch if c not in (self.drop, self.keep)]
    def full(self, P):
        Q = dict(P); Q[self.drop] = Q[self.keep]; return Q
    def value(self, col, P, colour=True):
        Q = self.full(P)
        v = self.kern(Q) * (self.colour(col, Q) if colour else 1.0)
        if self.kern2 is not None:
            v = v + self.kern2(Q) * (self.colour2(col, Q) if colour else 1.0)
        return self.coef * v
    def label(self): return '%s %s |%s->%s' % (self.row, self.name, self.drop, self.keep)

U = []
# ---------------- Row I.III, Eq. (sIIII)
for kn, kf, u in (('Kz', lib.K_III_z, 'z'), ('Kw', lib.K_III_w, 'w')):
    U.append(Unit('I.III', kn, ("x'", 'x'), G1, kf, "N D(x,y) - N D(x,%s) - F(x,y,w') + F(x,%s,w')" % (u, u), -1, fourier=('w', "w'")))
    U.append(Unit('I.III', kn, ('y', 'x'), G1, kf, "N D(x,x') - N D(x,w') - F(x,x',%s) + F(x,%s,w')" % (u, u), -1, fourier=('w', "w'")))
    U.append(Unit('I.III', kn, ("x'", 'y'), G1, kf, "N D(x,y) - F(x,y,%s) - F(x,y,w') + F(x,%s,w')" % (u, u), +1, fourier=('w', "w'")))
# ---------------- Row I.V first part, Eq. (sIV1)
U.append(Unit('I.V', 'KV1', ("x'", 'y'), G15, lib.K_V1, "N D(y,y') - F(x,y,z) + F(x,z,w') - F(y,y',w')", -0.5, fourier=('z', "w'")))
U.append(Unit('I.V', 'KV1', ("x'", "y'"), G15, lib.K_V1, "N D(y,y') - F(x,y',z) + F(x,z,w') - F(y,y',w')", +0.5, fourier=('z', "w'")))
# ---------------- Row II.I first part, Eq. (sIII1)
KT = lambda P: lib.K_T(P, False); KTs = lambda P: lib.K_T(P, True)
U.append(Unit('II.I', 'KT(w,w\')', ('y', 'x'), G2, KT, "N D(x,x') + N D(w,w') - F(x,z,w') - Q(x,z;w,x')", +1, fourier=('w', "w'")))
U.append(Unit('II.I', 'KT(w,w\')', ("x'", 'x'), G2, KT, "N D(x,y) + N D(w,w') - F(x,z,w') - Q(x,w;z,y)", -1, fourier=('w', "w'")))
U.append(Unit('II.I', 'KT(w\',w)', ("x'", 'x'), G2, KTs, "N D(x,y) + N D(w,w') - F(x,z,w) - Q(x,w';z,y)", -1, fourier=('w', "w'")))
U.append(Unit('II.I', 'KT(w\',w)', ('y', 'x'), G2, KTs, "N D(x,x') + N D(w,w') - F(x,z,w) - Q(x,z;w',x')", +1, fourier=('w', "w'")))
# ---------------- Row III.II, Eq. (sIIIII)
U.append(Unit('III.II', 'KT(w,w\')', ("x'", 'y'), G2, KT, "N D(x,y) + N D(y,w) - F(x,y,z) - F(y,z,w)", +1, fourier=('w', "w'")))
U.append(Unit('III.II', 'KT(w,w\')', ("x'", 'x'), G2, KT, "N D(x,y) - N D(x,z) + N D(x,w) - Q(x,w;z,y)", +1, fourier=('w', "w'")))
U.append(Unit('III.II', 'KT(w,w\')', ('y', 'x'), G2, KT, "N D(x,x') + N D(x',w) - F(x,x',z) - Q(x,z;w,x')", -1, fourier=('w', "w'")))

# ---------------- Group III: (contact, sign, colour of the L(y') line, colour of the L(z) line, sign of the L(z) line)
G3ROWS = {
 'III.III': [(("x'", 'y'), +1, "-N D(x,y) + F(x,y,y') + F(y,z,w) - F(y',z,w)", "-N D(x,y) - N D(z,w) + F(x,y,z) + F(y,z,w)", +1),
             (("x'", 'x'), -1, "-N D(x,y) + F(x,y,y') + F(x,z,w) - F(y',z,w)", "-N D(x,y) - N D(z,w) + F(x,y,z) + F(x,z,w)", -1)],
 'III.V':   [(("x'", 'x'), +1, "-N D(x,y) + N D(x,w) + F(x,y,y') - F(x,y',w)", "-N D(x,y) + N D(x,w) + F(x,y,z) - F(x,z,w)", +1),
             (('y', 'x'), +1, "-N D(x,x') + N D(x,y') + F(x,x',w) - F(x,y',w)", "-N D(x,x') + N D(x,z) + F(x,x',w) - F(x,z,w)", +1),
             (("x'", 'y'), -1, "-N D(x,y) + F(x,y,y') + F(x,y,w) - F(x,y',w)", "-N D(x,y) + F(x,y,z) + F(x,y,w) - F(x,z,w)", -1)],
 'III.VI':  [(("x'", 'y'), +1, "N D(x,y) - N D(y,z) - F(x,y,y') + F(y,y',z)", "N D(x,y) - F(x,y,z)", +1),
             (('y', 'x'), +1, "N D(x,x') - N D(x,y') - F(x,x',z) + F(x,y',z)", "N D(x,x') - F(x,x',z)", +1),
             (("x'", 'x'), -1, "N D(x,y) - F(x,y,y') - F(x,y,z) + F(y,y',z)", "N D(x,y) + N D(y,z) - 2F(x,y,z)", -1)],
}
for row, lst in G3ROWS.items():
    for contact, sg, cy, cz, sz in lst:
        assert sg == sz
        U.append(Unit(row, 'Ehat', contact, G3, lib.K_Ehat, cy, -sg, fourier=('w', "y'")))
        U.append(Unit(row, 'Psi', contact, G3, lib.K_Psi, cy, +sg, fourier=('w', "y'")))
        U.append(Unit(row, 'delta', contact, G3, lib.g_delta_yp, cy, sg, dist='delta', kern2=lib.g_delta_z, colour2=cz, fourier=('w', "y'")))
        U.append(Unit(row, 'deltaPV', contact, G3, lib.g_delta_yp, cy, sg, dist='pv', kern2=lib.g_delta_z, colour2=cz, fourier=('w', "y'")))
        U.append(Unit(row, 'delta(z) alone', contact, G3, lib.g_delta_z, cz, sg, dist='delta', fourier=('w', "y'")))
