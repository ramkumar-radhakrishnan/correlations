# Orientation rule for adjoint Wilson lines in expressions without f's (eq.(4rho)):
#   U^{ij}(p): the index contracted with a rho (or with the rho-side index of a second U(z)) is the SECOND one;
#   first indices are contracted among U's.  Returns list of flips (index pair, point) and applies them.
def fix_orientation(m):
    rho_idx = set(i for i, _ in m.rhos)
    side = {}                      # index -> 'rho' or 'first'
    for i in rho_idx: side[i] = 'rho'
    Us = list(m.Us); done = [False] * len(Us)
    changed = True
    while changed:
        changed = False
        for n, (i, j, p) in enumerate(Us):
            if done[n]: continue
            si, sj = side.get(i), side.get(j)
            if si is None and sj is None:
                # U(z) U(z) sharing an index -> shared index is rho-side
                for k, (i2, j2, p2) in enumerate(Us):
                    if k != n and p2 == p:
                        sh = set((i, j)) & set((i2, j2))
                        if len(sh) == 1:
                            s = sh.pop(); side[s] = 'rho'; changed = True
                continue
            if si is None: side[i] = 'first' if sj == 'rho' else 'rho'
            if sj is None: side[j] = 'first' if si == 'rho' else 'rho'
            assert side[i] != side[j], ("inconsistent", m.src, Us[n])
            done[n] = True; changed = True
    assert all(done), ("not oriented", m.src)
    flips = []
    newU = []
    for (i, j, p) in Us:
        if side[i] == 'rho':
            flips.append((i + j, p)); newU.append((j, i, p))
        else:
            newU.append((i, j, p))
    m.Us = newU
    return flips
