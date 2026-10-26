"""bfkl.tex + F3 + the y' -> z terms of Group III Rows III(1), III(2), V, VI: three-rho cancellation, two-rho terms vs B - J,
and row by row vs the earlier verified set."""
import json
import env  # noqa: F401
from core import symmetrize
import rows as R
import two_rho_engine as TE
import three_rho_engine as T3
import rowcmp as RC


def all2(sign=-1):
    FT = TE.file_terms(True)
    for k, f in R.THREE_YZ.items():
        FT[k] = [s for t in f() for s in symmetrize(t, keep=(2,))]
    FT['G3 III(2) yz'] = [s for t in R.g3_iii2_P(sign) for s in symmetrize(t, keep=(2,))]
    return FT


if __name__ == '__main__':
    out = {}
    # three rho
    S = T3.sym3(True)
    for k, f in R.THREE_YZ.items():
        S[k] = [s for t in f() for s in symmetrize(t, keep=(3,))]
    r3 = []
    for seed in (1, 2, 3):
        cfg = T3.Cfg3(seed)
        tot, sc = 0.0, 0.0
        for ts in S.values():
            v, s = T3.value3(ts, cfg)
            tot += v
            sc += s
        r3.append(abs(tot) / sc)
    out['three-rho |sum|/sum|terms|'] = r3
    print('three-rho with F3 + yz:', r3)
    J, B4, Bref = TE.targets()
    for sign in (-1, 1):
        FT = all2(sign)
        groups = dict(FT)
        groups['-(B - J)'] = TE.neg(Bref) + J[-1]['terms']
        res = []
        for seed in (1, 2, 3):
            acc = TE.keyed(groups, seed)
            res.append(len(TE.residual(acc)))
        out['two-rho vs B-J, III(2) yz sign %+d (bad keys)' % sign] = res
        print('III(2) yz sign %+d: bad keys vs B - J:' % sign, res, ' total keys', len(acc))
    json.dump(out, open('yz_test.json', 'w'), indent=1)
