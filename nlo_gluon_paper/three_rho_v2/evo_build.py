"""Books for the evolution note.
A side (coefficient of log(vee/Lambda)): '3rho' = two-rho terms of ThreeRho_final.tex, 'gen' = genuine, '4rho' = 4rho.tex lists, '-JIMWLK'.
P side (coefficient of log vee):        '3rho' = region-P two-rho terms of the three-rho rows (+ A,B, which have P = A),
                                         'gen' = genuine region-P terms, '4rho' = 4rho.tex lists (P = A), '-BFKL'.
python3 evo_build.py -> evo.pkl"""
import pickle
import sys
import env
from core import symmetrize
from book import Book
import final_eq as FE
import target as TG
import refP
import reduce as RD
import reduceP as RP
sys.path.insert(0, env.ALL2)
import sets as S2

THREE = ['3:I.II', '3:I.III', 'D', '3:I.IV', '3:I.V', '3:II.I', '3:III.I', '3:III.II', '3:III.III', '3:III.V', '3:III.VI', 'A,B']
GEN = ['g:I.II', 'g:I.V', 'g:II.I', 'g:II.IV', 'g:III.III']


def merged_book(groups, target, tname):
    B = Book()
    for t in target:
        B.add(dict(t, c=-t['c']), tname)
    for g, ts in groups.items():
        for t in ts:
            B.add(t, g)
    B.finalize()
    return B.dump()


def detailed_book(sets, target, tname):
    B = Book()
    for t in target:
        B.add(dict(t, c=-t['c']), tname)
    for g, ts in sets.items():
        for t in ts:
            B.add(t, g)
    B.finalize()
    return B.dump()


if __name__ == '__main__':
    final = [r for r in FE.reparse('ThreeRho_final.tex') if r['color'] != 'red']
    three_A = [s for r in final for t in FE.terms_of(r) for s in symmetrize(t, (2,))]
    A = S2.A_sets()
    P = S2.P_sets()
    out = {}
    gA = {'3rho': three_A, 'gen': [t for k in GEN for t in A[k]], '4rho': A['4rho']}
    gP = {'3rho': [t for k in THREE if k in P for t in P[k]], 'gen': [t for k in GEN if k in P for t in P[k]], '4rho': P['4rho']}
    for side, groups, target, tname, red in (('A', gA, TG.jimwlk(), '-JIMWLK', 'RD'), ('P', gP, refP.bfkl(), '-BFKL', 'RP')):
        bk = merged_book(groups, target, tname)
        R = {}
        for key, nets in bk.items():
            if red == 'RD':
                table, rel, cols = RD.reduce_structure(nets)
            else:
                table, rel, cols = RP.reduce_structure(nets, (tname,))
            R[key] = dict(table=table, relations=rel, cols=cols, nets=nets)
            nb = sum(1 for r in table if r['total'] != 0)
            print(side, len(nets), len(table), 'nonzero', nb, {g: sum(1 for r in table if g in r['contrib']) for g in ('3rho', 'gen', '4rho', tname)})
        out[side] = dict(book=bk, reduced=R)
    # per-source P two-rho terms (for the explicit list of the log(vee) coefficients)
    out['P_sets'] = {k: P[k] for k in P}
    pickle.dump(out, open('evo.pkl', 'wb'))
