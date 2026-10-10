"""BFKL target (region P), symbolic:  <Phi| O[rho + j_p] |Phi> - O  at order eps^2."""
from sym import Term
I=1j
def Mterms(pb,p,kern):
    """rho^{b'}(xb) [U(wb)-U(pb)]^{a b'} [U(w)-U(p)]^{a b} rho^b(x) with U's at given points"""
    out=[]
    for s1,q1 in [(1,'wb'),(-1,pb)]:
        for s2,q2 in [(1,'w'),(-1,p)]:
            out.append((s1*s2,[('a','bp',q1),('a','b',q2)]))
    return out
def bfkl_target():
    res=[]
    LOk=((('xb','wb'),('x','w')),)
    Kuv=(((('u','z'),('v','z'))),)
    # (i)  -1/2 K(u,v,z) [rho^h(u) rho^h(v) O - 2 rho^h(u) O rho^h(v) + O rho^h(u) rho^h(v)]
    for c,Us in Mterms('xb','x',None):
        O=[('bp','xb'),('b','x')]
        res.append(Term(-0.5*c,Us,[],[('h','u'),('h','v')]+O,LOk+Kuv,'Bi'))
        res.append(Term(+1.0*c,Us,[],[('h','u')]+O+[('h','v')],LOk+Kuv,'Bi'))
        res.append(Term(-0.5*c,Us,[],O+[('h','u'),('h','v')],LOk+Kuv,'Bi'))
    # (ii) K(u,v,z) rho^{e'}(u) rho^e(v) (T^{b'} T^b)_{e'e} M(z,z),  kernel K(z-wb).K(z-w)
    for c,Us in Mterms('z','z',None):
        res.append(Term(c*(-1),Us,[('bp','ep','g'),('b','g','e')],[('ep','u'),('e','v')],((('z','wb'),('z','w')),)+Kuv,'Bii'))
    # (iii) K(u,v,z) rho^{e'}(u) rho^{b'}(xb) rho^e(v) (T^b)_{e'e} M(xb,z)   kernel K(xb-wb).K(z-w)
    for c,Us in Mterms('xb','z',None):
        res.append(Term(c*(-I),Us,[('b','ep','e')],[('ep','u'),('bp','xb'),('e','v')],((('xb','wb'),('z','w')),)+Kuv,'Biii'))
    #      + K(u,v,z) rho^{e'}(u) rho^{b}(x) rho^e(v) (T^{b'})_{e'e} M(z,x)  kernel K(z-wb).K(x-w)
    for c,Us in Mterms('z','x',None):
        res.append(Term(c*(-I),Us,[('bp','ep','e')],[('ep','u'),('b','x'),('e','v')],((('z','wb'),('x','w')),)+Kuv,'Biii'))
    return res
