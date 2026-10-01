import sympy as sp

# ---------------------------------------------------------------------------
# 1) soft (xi -> 0) limit of the Phi.Phi* structures, using eta = 1 - zeta - xi
# ---------------------------------------------------------------------------
z_, e_, x_ = sp.symbols('zeta eta xi', positive=True)     # zeta, eta, xi
sub = {e_: 1 - z_ - x_}
zb  = 1 - z_                                              # zeta-bar

print("=== soft limits of the spinor structures ===")
# 1.1 antiquark emission
RR_aq  = 8*(e_**2 + (1-z_)**2)*(z_**2 + (1-z_)**2)                    # delta delta
eps_aq = -8*x_*(2*z_-1)*(2*(1-z_)-x_)                                 # eps eps
# 1.2 quark emission
RR_q   = 8*(z_**2 + (1-e_)**2)*(e_**2 + (1-e_)**2)
eps_q  = -8*x_*(2*e_-1)*(2*(1-e_)-x_)
# 1.3 interference
RR_x   = 8*(1-x_-2*z_*e_)*(x_*(1-x_)+2*z_*e_)
eps_x  = -8*x_*(1-2*z_-x_)**2

for name, expr, pre in [("1.1 dd", RR_aq, z_**2), ("1.1 ee", eps_aq, z_**2),
                        ("1.2 dd", RR_q , e_**2), ("1.2 ee", eps_q , e_**2),
                        ("1.3 dd", RR_x , -e_*z_), ("1.3 ee", eps_x , -e_*z_)]:
    lim = sp.simplify(sp.limit(sp.expand(expr.subs(sub)), x_, 0))
    tot = sp.factor(sp.simplify(sp.limit(sp.expand((pre*expr).subs(sub)), x_, 0)))
    print(f"  {name}:  Phi.Phi* -> {sp.factor(lim)}   |  prefactor*Phi.Phi* -> {tot}")

target = sp.factor(16*z_**2*zb**2*(z_**2+zb**2))
print("\n  common factor 16 zeta^2 zetabar^2 [zeta^2+zetabar^2] =", target)

# the O(xi) rate of approach (needed for the subtraction to be finite)
print("\n=== F(xi) - F(0) = O(xi) ? leading term in xi of prefactor*Phi.Phi* ===")
for name, expr, pre in [("1.1", RR_aq, z_**2), ("1.2", RR_q, e_**2), ("1.3", RR_x, -e_*z_)]:
    f = sp.expand((pre*expr).subs(sub))
    ser = sp.series(f, x_, 0, 2).removeO()
    print(f"  {name}: ", sp.factor(sp.simplify(sp.expand(ser - ser.subs(x_,0)))))

# reg x inst and inst x inst prefactors
print("\n=== reg x inst and inst x inst vanish as xi -> 0 ===")
print("  1.1 reg-inst prefactor 8 xi eta/(1-zeta) ->", sp.limit((8*x_*e_/(1-z_)).subs(sub), x_, 0))
print("  1.1 inst-inst prefactor 8 (xi zeta/(1-zeta))^2 ->", sp.limit((8*(x_*z_/(1-z_))**2).subs(sub), x_, 0))
print("  1.3 reg-inst second prefactor 8 xi zeta/(zeta+xi) ->", sp.limit(8*x_*z_/(z_+x_), x_, 0))
