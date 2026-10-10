"""Symbolic H_JIMWLK acting on the LO operator O = rho(xb)[U(wb)-U(xb)]^T [U(w)-U(x)] rho(x)."""
from sym import Term
import itertools
_c=itertools.count()
def LO_terms():
    out=[]
    for s1,p1 in [(1,'wb'),(-1,'xb')]:
        for s2,p2 in [(1,'w'),(-1,'x')]:
            out.append(Term(s1*s2,[('a','bp',p1),('a','b',p2)],[],[('bp','xb'),('b','x')],((( 'xb','wb'),('x','w')),)))
    return out
UPOS=('w','wb','x','xb')
def JR(t,h,spos_tag):
    """J_R^h(s) t : Leibniz over U's; returns list of (term, position)"""
    out=[]
    for k,(i,j,p) in enumerate(t.Us):
        e='r%d'%next(_c)
        Us=list(t.Us); Us[k]=(i,e,p)
        out.append((Term(t.c*(-1j),Us,t.fs+((h,e,j),),t.rs,t.kern),p))    # U T^h : U^{ie}(-i f^{h e j})
    return out
def JL(t,m):
    """J_L^m(s) t = U^{mg}(s) J_R^g(s) t  (the U^{mg}(s) factor is placed in front)"""
    out=[]
    for tt,p in JR(t,'g%d'%next(_c),None):
        g=tt.fs[-1][0]
        out.append((Term(tt.c,((m,g,p),)+tt.Us,tt.fs,tt.rs,tt.kern),p))
    return out
def H_on_LO(factor=-0.5):
    """factor * Htilde O ;  Htilde = sum K(s,s',z){J_R^h J_R^h + J_L^m J_L^m - U^{mh}(z)[J_R^h J_L^m + J_L^m J_R^h]}"""
    res=[]
    for t in LO_terms():
        for kind in ['RR','LL','RL','LR']:
            h='h%d'%next(_c); m='m%d'%next(_c)
            iname = h if kind in ('RR','LR') else m
            oname = h if kind in ('RR','RL') else m
            if kind=='LL': iname=oname=m
            inner = JR(t,iname,None) if kind[1]=='R' else JL(t,iname)
            for t1,p1 in inner:
                outer = JR(t1,oname,None) if kind[0]=='R' else JL(t1,oname)
                for t2,p2 in outer:
                    c=t2.c*factor; Us=t2.Us
                    if kind in ('RL','LR'):
                        c=-c; Us=((m,h,'z'),)+Us
                    res.append(Term(c,Us,t2.fs,t2.rs,t2.kern+(((p2,'z'),(p1,'z')),),'H'))
    return res
