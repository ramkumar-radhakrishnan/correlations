"""ref3_beta / b2: exact coefficients (sympy).
 [1] xi integrals of the paper (B.12/E.3/E.6) and of the quark-loop kernel (BLM remark only)
 [2] d-dimensional polarisation/spin sums: gluon loop  xi(1-xi)*Dpol/d ; quark loop 1-4xi(1-xi)/d
     -> constants in CDR (Dpol=d) and DRED-like (Dpol=2): shift of 67/9, of N_2b UV, of Row II.IV
 [3] LCPT identity behind Row I.II + Row I.V (unresolved pair):  1/(a b) + 1/(b D) = 1/(a D), b=a+D
 [4] pole table with general beta_0 (BLM bookkeeping) sums to -l/2 for any N_f
 [5] independent mu_R, mu_F: RG invariance of  g^2(mu_R) D(mu_F) [LO + NLO]
 [6] DRED <-> CDR: momentum sum rule fixes the delta term of the FF conversion = coupling conversion
Run: python3 -I b2_coefficients.py
"""
import sympy as sp
xi, d, eps, L, l, nf, TR, Nc = sp.symbols('xi d epsilon L ell n_f T_R N_c', positive=True)
dl = sp.symbols('delta', positive=True)
A = xi * (1 - xi)
Kgg = A - 2 + 1 / A                     # paper: [-2 + xi(1-xi) + 1/(xi(1-xi))] = P_gg^real/(2N_c)
Kq2 = xi**2 + (1 - xi)**2               # quark loop shape, units n_f T_R/N_c (b1 [2])
I = lambda f: sp.simplify(sp.integrate(sp.expand(f), (xi, 0, 1)))
print('[1] int (xi(1-xi)-2)                 =', I(A - 2), '   (non-soft part of K_gg)')
r = sp.integrate(Kgg, (xi, dl, 1 - dl)); print('    int_delta^{1-delta} K_gg          =', sp.series(sp.simplify(r), dl, 0, 1).removeO())
print('    int (xi(1-xi)-2) ln(xi(1-xi))     =', I((A - 2) * sp.log(A)))
print('    int xi^2+(1-xi)^2                 =', I(Kq2))
print('    int (xi^2+(1-xi)^2) ln(xi(1-xi))  =', I(Kq2 * sp.log(A)))
print('    int xi(1-xi)                      =', I(A))
# double-plus integral N(xi)/(xi_+ (1-xi)_+) with N = xi(1-xi)[xi(1-xi)+xi/(1-xi)+(1-xi)/xi]... = real P_gg/(2N_c) * xi(1-xi)
Nfun = A * Kgg
# int_0^1 dxi [N(xi) - N(0)(1-xi) - N(1) xi]/(xi(1-xi)) with N(0)=N(1)=1
print('    double-plus integral of real P_gg/(2N_c) =', sp.simplify(sp.integrate(sp.simplify((Nfun - (1 - xi) - xi) / A), (xi, 0, 1))))

# [2] polarisation sums.  Self energy:  Psi1 = -2 * norm * int K^{(d)} [-2/eps + ln(k^2 xi(1-xi)/mu^2)], d = 2-eps
def const_SE(Kd):
    K0 = Kd.subs(eps, 0); K1 = sp.diff(Kd, eps).subs(eps, 0)
    # coefficient of norm, eps^0, without the ln(k^2/mu^2) and the soft logs: -2[ int K0 ln A  - 2 int K1 ]
    return -2 * (sp.integrate(sp.expand(K0 * sp.log(A)), (xi, 0, 1)) - 2 * sp.integrate(sp.expand(K1), (xi, 0, 1)))
dd = 2 - eps
Kg_CDR = A * (dd / dd) - 2 + 1 / A
Kg_DR = A * (2 / dd) - 2 + 1 / A
# finite part: only the non-soft piece matters for the shift
cg_CDR = const_SE(A * (dd / dd) - 2); cg_DR = const_SE(A * (2 / dd) - 2)
print('[2] gluon loop constant (non-soft part, units of norm): CDR %s, Dpol=2: %s, shift %s' % (sp.nsimplify(cg_CDR), sp.nsimplify(cg_DR), sp.nsimplify(cg_DR - cg_CDR)))
print('    -> cross section (alpha_s N_c/2pi)(67/9-pi^2/3)LO becomes (alpha_s N_c/2pi)(67/9-1/3-pi^2/3)LO, i.e. K -> K - N_c/6')
Kq_d = 1 - 4 * A / dd
cq = const_SE(Kq_d)
print('    quark loop (CDR = Dpol=2, no gluon polarisation in the loop): constant %s (units n_f T_R/N_c) ; pole coeff -2*int = %s'
      % (sp.nsimplify(cq), -2 * I(Kq2)))
print('      -> 67/9 - pi^2/3 -> 67/9 - pi^2/3 - (20/9) n_f T_R/N_c ;  11/3 -> 11/3 - (4/3) n_f T_R/N_c = beta_0/N_c')
# N_2b UV of Row I.V: coefficient int K^{(d)} * (2/eps + B): xi(1-xi) * dperp * (1/d) ; CDR dperp=d, DRED dperp=2
for name, fac in [('CDR', dd / dd), ('Dpol=2', 2 / dd)]:
    c = sp.series(sp.integrate(sp.expand(A * fac), (xi, 0, 1)) * (2 / eps), eps, 0, 1).removeO()
    print('    N_2b UV, xi(1-xi) part x (2/eps):', name, sp.simplify(c), '  (finite extra: %s)' % sp.simplify(c - 2 / eps * sp.Rational(1, 6)))
# real collinear (Row II.IV): -C (2/eps) int P^{(d)} Pi ;  P^{(d)} = P + 2N_c xi(1-xi)(Dpol/d - 1)
extra = sp.series(-(2 / eps) * 2 * Nc * A * (2 / dd - 1), eps, 0, 1).removeO()
print('    Row II.IV, Dpol=2 minus CDR: -C int dxi [%s] Pi(xi)' % sp.simplify(-extra))

# [3] identity
a, D = sp.symbols('a D', positive=True); b = a + D
print('[3] 1/(a b) + 1/(b D) - 1/(a D) =', sp.simplify(1 / (a * b) + 1 / (b * D) - 1 / (a * D)))

# [4] pole table, units a*2/eps*LO (UV) and C*2/eps[.(x)Pi] ; beta0/(4N_c) replaces 11/12
b0 = sp.symbols('beta_0')
q = b0 / (4 * Nc); P = sp.symbols('PxPi')
rows = {'I.II': 2 * q - L + l / 2, 'I.V': -2 * q + 2 * L, 'II.IV JIMWLK': -L - l, 'II.IV real Pgg+': -P + q,
        'coupling CT': -q, 'FF CT': P}
print('[4] single-eps sum =', sp.simplify(sum(rows.values())), ' (independent of beta_0 and N_f)')
# two-eps: norm UV +q - L, norm IR -q + L
UV = rows['I.II'] + rows['I.V'] + (q - L) + rows['coupling CT']
IR = (-q + L) + rows['II.IV JIMWLK'] + rows['II.IV real Pgg+'] + rows['FF CT']
print('    two-eps: UV', sp.simplify(UV), ' IR', sp.simplify(IR))

# [5] mu_R, mu_F
als, b0s, LO, Lr, Lf = sp.symbols('alpha_s beta_0 LO L_R L_F')  # L_R = ln mu_R^2, L_F = ln mu_F^2
# alpha_s(mu_R): d alpha/d L_R = -beta0 alpha^2/(4 pi); D(mu_F) delta-part: dD/dL_F = +(alpha/2pi)(beta0/2) D  (+ real P_gg+ part, cancelled by DGLAP real logs)
# sigma = alpha(mu_R) D(mu_F) [LO + alpha (beta0/4pi)(L_R - L_F) LO + (mu-independent)]
f = sp.Function('al')(Lr); g = sp.Function('Dg')(Lf)
sig = f * g * (1 + f * b0s / (4 * sp.pi) * (Lr - Lf))
dR = sp.diff(sig, Lr).subs(sp.Derivative(f, Lr), -b0s * f**2 / (4 * sp.pi))
dF = sp.diff(sig, Lf).subs(sp.Derivative(g, Lf), f / (2 * sp.pi) * b0s / 2 * g)
print('[5] d sigma/d ln mu_R^2 =', sp.simplify(sp.expand(dR)), ' ; d sigma/d ln mu_F^2 =', sp.simplify(sp.expand(dF)), ' (O(alpha^3) only)')

# [6] DRED <-> CDR: real-collinear extra -C*2N_c xi(1-xi) per (x)Pi  -> FF conversion c(xi) = -2 N_c xi(1-xi) + c_delta delta(1-xi)
cdel = sp.symbols('c_delta')
mom = sp.integrate(xi * (-2 * Nc * A), (xi, 0, 1)) + cdel
print('[6] momentum sum rule for the conversion: c_delta =', sp.solve(mom, cdel)[0],
      ' ; (alpha_s/2pi) c_delta =', sp.solve(mom, cdel)[0] / (2 * sp.pi), '* alpha_s  vs  alpha_DR/alpha_MS - 1 = alpha_s N_c/(12 pi)')
