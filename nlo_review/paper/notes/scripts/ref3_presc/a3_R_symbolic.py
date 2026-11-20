# Ref.3 / Major 5(a).  Exact identity for the on-shell kernel of Group I Row III (see a2_rowIII_delta.py):
#   R = P^i Q^{ij}(x,z,w) K^j(y-z),   Q^{ij} = 2 * Brace^{ij},
# where Brace^{ij} is the coefficient of ln[(x-z)^2/(x-w)^2] in the combined kernel M^{ij}, Eq. (6.19) of the draft:
#   Brace^{ij} = d^{ij}/(2C) - a^i b^j/(AB) + [a^i c^j B/A + c^i b^j A/B - c^i a^j - b^i c^j] / (C (A-B)),
# a = x-w, b = x-z, c = z-w = a-b.  Consequence: an i*eta*0 prescription changes M^{ij} by ln(B/A) -> ln(B/A) - i*eta*pi.
import sympy as sp
a1, a2, b1, b2 = sp.symbols('a1 a2 b1 b2', real=True)
a = sp.Matrix([a1, a2]); b = sp.Matrix([b1, b2]); I2 = sp.eye(2)
dot = lambda u, v: (u.T * v)[0]
c = a - b
A, B, C = dot(a, a), dot(b, b), dot(c, c)
V = a - (A / B) * b
outer = lambda u, v: u * v.T
# Q^{ij}: R = P^i Q^{ij} Y^j / Y^2 with the three groups of residues (T1; T3; T4=T5 at p0), k+ = 1
Q = I2 * (A + B) / (B - A) ** 2 - 2 * outer(a, b) / (A * B) \
    + 2 * (outer(c, V) - (B / A) * outer(V, c) - (B / (B - A)) * I2 * dot(c, V)) / (C * (B - A))
Br = I2 / (2 * C) - outer(a, b) / (A * B) + (outer(a, c) * B / A + outer(c, b) * A / B - outer(c, a) - outer(b, c)) / (C * (A - B))
for i in range(2):
    for j in range(2):
        d = sp.simplify(sp.cancel(sp.together(Q[i, j] - 2 * Br[i, j])))
        print(f"Q^{i+1}{j+1} - 2 Brace^{i+1}{j+1} =", d)
# structure of Brace: alpha*delta + beta*epsilon (a rotation-dilation matrix)
print("Brace^11 - Brace^22 =", sp.simplify(sp.cancel(Br[0, 0] - Br[1, 1])), ";  Brace^12 + Brace^21 =", sp.simplify(sp.cancel(Br[0, 1] + Br[1, 0])))
al = sp.factor(sp.cancel(Br[0, 0])); be = sp.factor(sp.cancel(Br[0, 1]))
print("alpha = Brace^11 =", al); print("beta  = Brace^12 =", be)
# small-c behaviour (z -> w): Brace = d/(2C) + O(1/|c|)
t = sp.symbols('t', positive=True); ph = sp.symbols('phi', real=True)
sub = {b1: a1 - t * sp.cos(ph), b2: a2 - t * sp.sin(ph)}
print("C * alpha at z->w :", sp.simplify(sp.limit(sp.simplify((C * al).subs(sub)), t, 0)))
print("C * beta  at z->w :", sp.simplify(sp.limit(sp.simplify((C * be).subs(sub)), t, 0)))
