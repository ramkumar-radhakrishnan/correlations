# Ref.3 / Major 5(a).  Symbolic simplification of the on-shell kernel R of Group I Row III (see a2_rowIII_delta.py).
# Coordinates relative to x: a = x-w, b = x-z, Y = y-z, c = z-w = a-b.  k+ = 1 (R scales like 1/k+).
import sympy as sp
a1, a2, b1, b2, y1, y2, P1, P2 = sp.symbols('a1 a2 b1 b2 Y1 Y2 P1 P2', real=True)
a = sp.Matrix([a1, a2]); b = sp.Matrix([b1, b2]); Y = sp.Matrix([y1, y2]); P = sp.Matrix([P1, P2])
dot = lambda u, v: (u.T * v)[0]
c = a - b
A, B, C, Y2 = dot(a, a), dot(b, b), dot(c, c), dot(Y, Y)
V = a - (A / B) * b
Rvec = Y * (A + B) / (Y2 * (B - A) ** 2) - 2 * a * dot(b, Y) / (A * B * Y2) \
    + 2 * (c * dot(Y, V) - (B / A) * dot(Y, c) * V - B / (B - A) * Y * dot(c, V)) / (Y2 * C * (B - A))
R = sp.together(sp.expand(dot(P, Rvec)))
num, den = sp.fraction(R)
num = sp.factor(sp.expand(num))
den = sp.factor(den)
print("R = num/den with\n den =", den)
print(" num =", num)
# is (B-A) a factor of the numerator (pole cancellation)?
BmA = sp.expand(B - A)
q, r = sp.div(sp.Poly(sp.expand(num), a1, a2, b1, b2, y1, y2, P1, P2), sp.Poly(BmA, a1, a2, b1, b2, y1, y2, P1, P2))
print("remainder of num/(B-A):", r.as_expr() == 0)
q2, r2 = sp.div(q, sp.Poly(BmA, a1, a2, b1, b2, y1, y2, P1, P2))
print("remainder of num/(B-A)^2:", r2.as_expr() == 0)
Rs = sp.factor(sp.cancel(R))
print("\nR simplified =", Rs)
# try a decomposition in terms of simple kernels: write R = P.(alpha * a/A + beta * Y/Y2 + gamma * b/B + ...)
Rvec_s = sp.Matrix([sp.factor(sp.cancel(Rvec[i])) for i in range(2)])
print("\nRvec_1 =", Rvec_s[0])
