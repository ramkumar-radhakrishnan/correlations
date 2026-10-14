"""For every term of the list, find the matching earlier entry (ThreeRho_cancellation, as written; typos T1, T2 fixed)
and compare the written order of the charges."""
import numpy as np
import env
import parse as P
import cancel_data as CD
from sources import SOURCES
from cmp_entries import fp_user, fp_entry

SRC = {s.key: s for s in SOURCES}
KEYS = {'4rho.a': ['A1', 'A2'], '4rho.b': ['B1', 'B2'], 'I.II.1': ['GI.II'], 'I.II.2': ['GI.II'], 'I.III': ['GI.III'],
        'I.IV': ['GI.IV'], 'I.V.1': ['GI.V1', 'GI.V2'], 'I.V.2': ['GI.V1', 'GI.V2'], 'II.I': ['GII.I1', 'GII.I2'],
        'III.I': ['GIII.I'], 'III.II': ['GIII.II'], 'III.III': ['GIII.III'], 'III.V': ['GIII.V'], 'III.VI': ['GIII.VI']}


def fixed_entries():
    E = P.parse()
    for e in E:
        if e['block'] == 'III.I' and e['n'] == 15:
            e['rs'] = [('l', 'x'), ('e', 'y'), ('c', "x'")]
        if e['block'] == 'III.II' and e['n'] == 13:
            e['Us'] = [('d', 'e', 'y'), ('h', 'l', 'x'), ('b', 'd', 'z'), ('h', 'a', "w'")]
    return E


if __name__ == '__main__':
    E = fixed_entries()
    from collections import Counter
    stats = Counter()
    for e in E:
        ents = []
        for k in KEYS[e['block']]:
            ents += CD.source_entries(SRC[k])[0]
        f = fp_user(e)
        hits = [en for en in ents if en['cls'] == e['cls'] and np.allclose(f, fp_entry(en), atol=1e-10)]
        if not hits:
            print('NO MATCH', e['block'], e['n'])
            continue
        en = hits[0]
        ou = [p for _, p in e['rs']]
        oe = [p for _, i, p in [r for r in en['col'] if r[0] == 'rho']]
        same = ou == oe
        stats[(e['block'], en['cc'], same)] += 1
    for k, v in sorted(stats.items()):
        print(k, v)
