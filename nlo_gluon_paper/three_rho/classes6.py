"""Replace part 2 of Row III by its near-diagonal leading log written in the notes' own variables:
   D = -2 i N int KK(x'-w').KK(y-w) KK(x-z).KK(w-z) f^{abd}[U^{ab'}(x')-U^{ab'}(w')] rho^{b'}(x')
           U^{de}(x) rho^e(x) [U^{bc}(w) - U^{bc}(y)] rho^c(y)  + c.c.
(Row III with the gluon position z -> w and the soft gluon at z between the B3 charge x and w.)"""
import sympy as sp
import classes4 as C4
from engine import tex_colour, I as SI
from sources import SOURCES

def row3_diag_terms():
    src = [s for s in SOURCES if s.key == 'GI.III'][0]
    out = []
    monos = [(sp.Integer(1), [('KK', "x'", "w'", 'y', 'w'), ('KK', 'x', 'z', 'w', 'z')])]
    for n, (sign, colour, word) in enumerate(src.subterms):
        col = [('U', u[1], u[2], 'w') if (u[0] == 'U' and u[3] == 'z') else u for u in colour]
        c = -2 * SI * sign
        out.append(C4.ct('Row III near p+=k+', False, n + 1, c, col, word, monos, {"w'": "w'", 'w': 'w', 'z': 'z'}))
        out.append(C4.ct('Row III near p+=k+', True, n + 1, sp.conjugate(c), col, word, monos, {"w'": 'w', 'w': "w'", 'z': 'z'}))
    return out

if __name__ == '__main__':
    terms = C4.notes_terms() + C4.fix_terms() + row3_diag_terms()
    S = C4.summarize(C4.group(C4.classify(terms)))
    ok = True
    for cls, rows in S.items():
        bad = [(n, t) for n, t in rows if t != 0]
        ok &= not bad
        print(f'{cls:7s}: {len(rows)} colour networks, {len(bad)} with nonzero sum')
        for n, t in bad:
            print('      ', t, tex_colour(n['rep']['col']))
    print('ALL THREE-RHO TERMS CANCEL' if ok else 'residual remains')
