# Order of limits: FF counterterm integrated over the real-emission range a=k+/vee < xi < 1.
# Units: (2/eps) * 2Nc*C*Pi(1) = (2/eps) * u * LO  (u = alpha_s Nc/pi).  Test function Pi(xi)=xi^2*(1+c*xi) (Pi->0 at 0).
import sympy as sp
xi,a,lam,c=sp.symbols('xi a lambda c',positive=True)
Pi=xi**2*(1+c*xi); Pi1=Pi.subs(xi,1)
b0_over_2Nc=sp.Rational(11,12)          # (beta0/2)/(2Nc)
Preal=xi*(1-xi)+xi/(1-xi)+(1-xi)/xi     # P_gg+/(2Nc)
L=sp.log(1/lam); ell=sp.log(1-a)
# Row IV real pole (IR), range a<xi<1-lam, sign -:
real=-sp.integrate(sp.expand(Preal*Pi),(xi,a,1-lam))
# FF CT, reading R1: distribution P_gg on [0,1] applied to theta(xi>a)Pi  (plus subtraction over [0,1])
reg=lambda f: sp.integrate(sp.expand(xi*(1-xi)*f+(1-xi)/xi*f),(xi,a,1))
plusR1=sp.integrate(sp.cancel((xi*Pi-Pi1)/(1-xi)),(xi,a,1)) - Pi1*sp.integrate(1/(1-xi),(xi,0,a))
CT1=reg(Pi)+plusR1+b0_over_2Nc*Pi1
# reading R2: plus subtraction restricted to [a,1]
plusR2=sp.integrate(sp.cancel((xi*Pi-Pi1)/(1-xi)),(xi,a,1))
CT2=reg(Pi)+plusR2+b0_over_2Nc*Pi1
# full-range CT (old table)
CT0=reg(Pi).subs(a,0)+sp.integrate(sp.cancel((xi*Pi-Pi1)/(1-xi)),(xi,0,1))+b0_over_2Nc*Pi1
normIR=(-sp.Rational(11,12)+L)*Pi1
UV=(sp.Rational(11,6)-L+ell/2 + (-sp.Rational(11,6)+2*L) + (sp.Rational(11,12)-L) - sp.Rational(11,12))*Pi1
def tot(CT,norm=True):
    s=real+CT+(normIR if norm else 0)+(UV if norm else (UV-(sp.Rational(11,12)-L)*Pi1))
    s=sp.series(sp.simplify(s),lam,0,1).removeO()
    return sp.simplify(s/Pi1)
for name,CT in [('full range (old)',CT0),('range a<xi<1, R1',CT1),('range a<xi<1, R2',CT2)]:
    t2=tot(CT); t1=tot(CT,False)
    print(name,': two-eps total =',sp.simplify(t2),'| single-eps =',sp.simplify(t1),'| series in a:',sp.series(t2,a,0,3))
    irs=sp.simplify(sp.series(sp.simplify((real+CT+normIR)),lam,0,1).removeO()/Pi1)
    print('    IR part only:',irs,' series:',sp.series(irs,a,0,3))
