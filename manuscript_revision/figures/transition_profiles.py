# Finite-volume ground states of the 1D functional (Eq. 6 of the manuscript) just below and just
# above the vacuum-CSL transition of Fig. 1, for fixed system size Lbar. Requires numpy, scipy, matplotlib.
import numpy as np
from scipy.optimize import minimize, brentq

LBAR, N = 10.0, 801
z = np.linspace(-LBAR/2, LBAR/2, N); h = z[1]-z[0]
w = np.ones(N); w[0] = w[-1] = 0.5

def energy(phi, b):
    d = np.diff(phi)/h
    return h*np.sum(0.5*d*d - b*d) + h*np.sum(w*(1-np.cos(phi)))

def grad(phi, b):
    d = np.diff(phi)/h; g = np.zeros(N); g[:-1] -= d-b; g[1:] += d-b
    return g + h*w*np.sin(phi)

def branch(b, solitons):
    """Local minimum continuously connected to a guess with the given number of solitons."""
    x0 = 0.01*z if solitons == 0 else 4*np.arctan(np.exp(z))
    r = minimize(energy, x0, args=(b,), jac=grad, method='L-BFGS-B',
                 options={'maxiter': 50000, 'gtol': 1e-11, 'ftol': 1e-15})
    return r.fun, r.x

bcr = brentq(lambda b: branch(b, 0)[0] - branch(b, 1)[0], 1.0, 2.0, xtol=1e-6)
print(f"Lbar = {LBAR}: bbar_cr = {bcr:.4f}  (4/pi = {4/np.pi:.4f})")

if __name__ == "__main__":
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 9,
                         "axes.linewidth": 0.6, "xtick.direction": "in", "ytick.direction": "in",
                         "xtick.top": True, "ytick.right": True})
    db = 0.05
    cases = [(bcr-db, 0, "-", rf"$\bar b={bcr-db:.2f}$ (below)"),
             (bcr+db, 1, "--", rf"$\bar b={bcr+db:.2f}$ (above)")]
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(3.4, 4.2), sharex=True)
    for b, ns, ls, lab in cases:
        E0, E1 = branch(b, 0)[0], branch(b, 1)[0]
        E, phi = branch(b, ns)
        assert E <= min(E0, E1) + 1e-9, "plotted state is not the ground state"
        phi = phi - 2*np.pi*np.round(phi[N//2]/(2*np.pi))        # fix the 2*pi ambiguity: phi(0) in (-pi, pi]
        print(f"  bbar={b:.3f}: E(0 solitons)={E0:.4f}, E(1 soliton)={E1:.4f}, phi(0)/pi={phi[N//2]/np.pi:.3f}, "
              f"Delta phi/2pi={(phi[-1]-phi[0])/(2*np.pi):.3f}")
        a1.plot(z, phi/np.pi, "k", ls=ls, lw=1.2, label=lab)
        a2.plot(z[:-1]+h/2, np.diff(phi)/h, "k", ls=ls, lw=1.2)
    for y in (0, 1):
        a1.axhline(y, color="0.6", lw=0.5, ls=":", zorder=0)
    a1.set_ylabel(r"$\phi_0/\pi$"); a2.set_ylabel(r"$\partial_{\bar z}\phi_0$"); a2.set_xlabel(r"$\bar z$")
    a1.legend(frameon=False, loc="upper left", fontsize=8)
    a1.text(0.97, 0.06, "(a)", transform=a1.transAxes, ha="right"); a2.text(0.97, 0.88, "(b)", transform=a2.transAxes, ha="right")
    a2.set_xlim(-LBAR/2, LBAR/2); a2.set_ylim(bottom=0)
    fig.tight_layout(h_pad=0.4)
    fig.savefig("fig_transition.pdf"); fig.savefig("fig_transition.png", dpi=200)
