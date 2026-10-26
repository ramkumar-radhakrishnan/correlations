"""Numerical check of the split, term by term: original = J + B + F_exact (any vee, Lambda), F -> F_exact for Lambda -> 0,
the 1/Lambda terms of Group I Row III add up to zero, your two red terms add up to zero.  Output: check_split.json"""
import json
import numpy as np
import blocks as BL
import gen as GN
import texeval as TE

BLK, SPL = GN.all_splits()


def val(B, t, env):
    return TE.term_value(B.term_full(0, t) if False else B.prefix + ' (' + t + ') ' + B.suffix, env)


def E1fix(t):
    return t.replace(r'{\color{blue} \log{\left(\frac{k^{+}}{\Lambda}\right)}}', r'2\,\log{\left(\frac{k^{+}}{\Lambda}\right)}')


if __name__ == '__main__':
    res, skipped = {}, []
    power, red = {}, {}
    for seed in (1, 2, 3):
        c = TE.cfg(seed)
        for kp, V, L in ((1.0, 37.0, 0.013), (0.6, 900.0, 3e-4)):
            env = TE.Env(c, KP=kp, V=V, L=L)
            envs = TE.Env(c, KP=kp, V=V, L=1e-9)
            for name, B in BLK.items():
                for n, t in enumerate(B.terms):
                    s = SPL[name][n]
                    try:
                        o = val(B, E1fix(t) if s.get('E1') else t, env)
                    except TE.TexError as e:
                        if 'integral' in str(e):
                            skipped.append((name, n))
                            continue
                        raise
                    if s.get('power'):
                        power.setdefault((seed, kp), []).append(o)
                        continue
                    if s.get('red'):
                        red.setdefault((seed, kp), []).append(o)
                        continue
                    parts = sum(val(B, s[k], env) for k in ('J', 'B', 'Fx') if s.get(k) is not None)
                    rel = abs(o - parts) / max(abs(o), 1e-30)
                    lim = abs(val(B, s['F'], envs) - val(B, s['Fx'], envs)) / max(abs(o), 1e-30) if s.get('F') else 0.0
                    key = '%s #%d' % (name, n)
                    r = res.setdefault(key, dict(identity=0.0, Lambda_to_0=0.0))
                    r['identity'] = max(r['identity'], float(rel))
                    r['Lambda_to_0'] = max(r['Lambda_to_0'], float(lim))
    out = dict(terms=res, worst_identity=max(r['identity'] for r in res.values()),
               worst_Lambda_to_0=max(r['Lambda_to_0'] for r in res.values()),
               power_terms_sum_rel=max(abs(sum(v)) / max(abs(x) for x in v) for v in power.values()),
               red_terms_sum_rel=max(abs(sum(v)) / max(abs(x) for x in v) for v in red.values()),
               not_evaluated=sorted(set('%s #%d' % k for k in skipped)), n_checked=len(res))
    for k in ('worst_identity', 'worst_Lambda_to_0', 'power_terms_sum_rel', 'red_terms_sum_rel', 'n_checked'):
        print(k, out[k])
    bad = {k: v for k, v in res.items() if v['identity'] > 1e-10 or v['Lambda_to_0'] > 1e-5}
    print('bad', bad)
    print('not evaluated (remaining p+/xi/momentum integrals):', out['not_evaluated'])
    json.dump(out, open('check_split.json', 'w'), indent=1)
