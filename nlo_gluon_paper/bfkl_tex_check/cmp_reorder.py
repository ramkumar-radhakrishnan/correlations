r"""Every 'written as' two-rho list of a three-rho row of BFKL.tex against the two-rho remainder (reorder.py) of the
row's first form (or of the corrected first form).  -> cmp_reorder.json"""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import reorder as RO
import cmp_first as C
sys.path.insert(0, os.path.join(HERE, '..', 'evol_check'))
from cmp_three import kdot
from rows_user import rows

ROWS = {'G1 II': ([106], [114], True), 'G1 III': ([126, 133], [146, 177], True), 'G1 IV': ([188], [196, 230], True),
        'G2 I first part': ([247, 253], [260], False), 'G3 I': ([359], [366, 396], True),
        'G3 II': ([408], [415, 448], True), 'G3 III first part': ([462], [469], True), 'G3 V': ([523], [530], True),
        'G3 VI': ([549], [556], True)}


def run(only=None, first_override=None, seeds=(1, 2)):
    U = {r[3]: r[4] for r in rows('BFKL.tex')}
    out = {}
    for name, (first, written, cc) in ROWS.items():
        if only and name not in only:
            continue
        F = [C.prep(C.plus(U[l])) for l in first]
        if first_override and name in first_override:
            F = [C.prep(t) for t in first_override[name]]
        Wr = [kdot(C.prep(U[l])) for l in written]
        r = []
        for s in seeds:
            c = C.T3.configs(s)
            a = sum(RO.remainder_cc(t, c, cc) for t in F)
            b = sum(RO.two_value(t, c) for t in Wr)
            r.append((complex(b / a), abs(a - b) / max(abs(b), 1e-300)))
        out[name] = dict(ratio_written_over_remainder=[str(complex(round(x[0].real, 10), round(x[0].imag, 10))) for x in r],
                         rel=max(x[1] for x in r))
        print(name, out[name], flush=True)
    return out


if __name__ == '__main__':
    o = run(only=sys.argv[1:] or None)
    json.dump(o, open(os.path.join(HERE, 'cmp_reorder.json'), 'w'), indent=1)
