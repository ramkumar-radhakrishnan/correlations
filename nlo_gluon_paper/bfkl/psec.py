"""Step-by-step reordering of the region-P three-rho sources (generator of the three-rho note)."""
import env  # noqa
import make_tex as MT
from engine import Rho, pieces_for, make_piece, simplify, canonical_rename
import userP as UP


def results(src):
    sres = {'src': src, 'subs': []}
    for n, (sign, colour, word) in enumerate(src.subterms):
        rhos = [Rho(i, p) for i, p in word]
        pcs = []
        for coeff, X, Y, Z, tag in pieces_for(src.form, rhos):
            c0, facs0, (old, new) = make_piece(colour, sign * coeff, X, Y, Z)
            c1, facs1, rules = simplify(c0, facs0)
            facs1c = canonical_rename(facs1) if facs1 else facs1
            pcs.append(dict(tag=tag, X=(X.idx, X.pos), Y=(Y.idx, Y.pos), Z=(Z.idx, Z.pos), coeff_before=c0, facs_before=facs0,
                            subst=(old, new), coeff_after=c1, facs_after=facs1c, rules=rules))
        sres['subs'].append(dict(sign=sign, colour=colour, word=word, pieces=pcs))
    return sres


def section(key):
    src = UP.PSRC[key]
    txt, res = MT.gen_source(src, results(src))
    return txt


if __name__ == '__main__':
    t = section('GI.IV.P')
    print(len(t))
    print(t[:3000])
