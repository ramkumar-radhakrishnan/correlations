r"""Checks of the split of your new rows into log(k+/Lambda) (JIMWLK), log(vee/k+) (BFKL) and the rest.

 split:  every pointwise-evaluable term = J + B + F_exact (three configurations, two (k+, vee, Lambda)); F -> F_exact for
         Lambda -> 0; the 1/Lambda terms of Group I Row III and your two red terms of d^3N_2 add up to zero.
 three:  the three-rho terms of the J list and of the B list (rows + the three-rho lists of 4rho.tex), symmetrized over
         the 6 placements of the charges, with c.c. for Groups I and III (not for Group II), cancel.
         J: with D (fix F3) and, as a variant, with your edge-kernel line.  B: with (B) (exact identity (I3)).
 poles:  the poles x log(k+/Lambda) of Group I Row II, Group I Row V 2b UV and Group II Row IV, ratio and sum.
Outputs: ll_split.json, ll_three.json, ll_poles.json"""
import json
import os
import sys
import numpy as np
import ll  # noqa: F401  (new inputs)
import gen as GN
import texeval as TE
import three_check as T3
import pole_check as PC

HERE = ll.HERE


def val(B, t, env):
    return TE.term_value(B.prefix + ' (' + t + ') ' + B.suffix, env)


def split_check():
    BLK, SPL = GN.all_splits()
    res, skipped, power, red = {}, [], {}, {}
    for seed in (1, 2, 3):
        c = TE.cfg(seed)
        for kp, V, L in ((1.0, 37.0, 0.013), (0.6, 900.0, 3e-4)):
            env = TE.Env(c, KP=kp, V=V, L=L)
            envs = TE.Env(c, KP=kp, V=V, L=1e-9)
            for name, B in BLK.items():
                for n, t in enumerate(B.terms):
                    s = SPL[name][n]
                    try:
                        o = val(B, t, env)
                    except TE.TexError as e:
                        if 'integral' in str(e):
                            skipped.append('%s #%d' % (name, n))
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
                    r = res.setdefault('%s #%d' % (name, n), dict(identity=0.0, Lambda_to_0=0.0))
                    r['identity'] = max(r['identity'], float(rel))
                    r['Lambda_to_0'] = max(r['Lambda_to_0'], float(lim))
    out = dict(terms=res, worst_identity=max(r['identity'] for r in res.values()),
               worst_Lambda_to_0=max(r['Lambda_to_0'] for r in res.values()),
               power_terms_sum_rel=max(abs(sum(v)) / max(abs(x) for x in v) for v in power.values()),
               red_terms_sum_rel=max(abs(sum(v)) / max(abs(x) for x in v) for v in red.values()),
               not_evaluated=sorted(set(skipped)), n_checked=len(res),
               n_J=sum(1 for L in SPL.values() for s in L if s.get('J')), n_B=sum(1 for L in SPL.values() for s in L if s.get('B')))
    for k in ('worst_identity', 'worst_Lambda_to_0', 'power_terms_sum_rel', 'red_terms_sum_rel', 'n_checked', 'n_J', 'n_B'):
        print(k, out[k], flush=True)
    return out


def three_check():
    seeds = (11, 12, 13)
    out = {}
    for part in ('J', 'B'):
        vals = T3.values(part, seeds)
        base = [lab for lab in vals if 'as written' not in lab]
        variants = {'final': base}
        if part == 'J':
            variants['edge kernel (as written) instead of D'] = [l for l in base if l != 'G1 III N3a #1'] + \
                ['G1 III N3a #1 (edge kernel, as written)']
        else:
            variants['without (B)'] = [l for l in base if not l.startswith('(B)')]
        res = {}
        for vname, labels in variants.items():
            r = [T3.total(vals, s, labels, False) for s in seeds]
            res[vname] = dict(relative=max(t / b for t, b in r), residual=[t for t, b in r], largest=[b for t, b in r])
            print('%s | %-45s worst rel. residual %.2e (%d terms)' % (part, vname, res[vname]['relative'], len(labels)), flush=True)
        table = {lab: dict(direct=v[seeds[0]][0].imag, cc=v[seeds[0]][1].imag, cc_added=v[seeds[0]][2]) for lab, v in vals.items()}
        out[part] = dict(variants=res, seed11_values_imag=table, terms=base)
    return out


def pole_check():
    BLK, SPL = GN.all_splits()
    terms = {}
    for name in PC.NAMES:
        B = BLK[name]
        n = 1 if name == 'G2 IV' else 0
        terms[name] = B.prefix + ' (' + SPL[name][n]['J'] + ') ' + B.suffix
    out = {}
    for seed in (21, 22, 23):
        c = TE.cfg(seed)
        P = {name: PC.pole(t, c, not name.startswith('G2')) for name, t in terms.items()}
        scale = max(abs(v) for v in P.values())
        out[str(seed)] = dict(poles_relative_to_G2_IV={k: float(np.real(v / P['G2 IV'])) for k, v in P.items()},
                              sum_relative=float(abs(sum(P.values())) / scale))
        print(seed, out[str(seed)], flush=True)
    return out


if __name__ == '__main__':
    what = sys.argv[1:] or ['split', 'three', 'poles']
    if 'split' in what:
        json.dump(split_check(), open(os.path.join(HERE, 'll_split.json'), 'w'), indent=1)
    if 'poles' in what:
        json.dump(pole_check(), open(os.path.join(HERE, 'll_poles.json'), 'w'), indent=1)
    if 'three' in what:
        json.dump(three_check(), open(os.path.join(HERE, 'll_three.json'), 'w'), indent=1)
