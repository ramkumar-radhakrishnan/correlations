"""Figures for the CSL seminar talk. All curves computed from the analytic
formulae of Son-Stephanov (2007), Brauner-Yamamoto (2017) and
Brauner-Radhakrishnan (2026)."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import ellipj, ellipe, ellipk
from scipy.optimize import brentq

import os
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figs") + os.sep
plt.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "cm",
    "font.size": 15, "axes.labelsize": 17, "legend.fontsize": 13,
    "axes.linewidth": 1.1, "lines.linewidth": 2.2,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "figure.autolayout": True, "savefig.transparent": True,
})
BLUE, RED, GREEN, ORANGE, GRAY = "#1f4e79", "#c0392b", "#1b7a5a", "#e08214", "#7f7f7f"

E = lambda k: ellipe(k*k)          # complete elliptic integral, modulus k
K = lambda k: ellipk(k*k)
dn = lambda u, k: ellipj(u, k*k)[2]

# ---------------------------------------------------------------- Fig: CSL profile
def fig_profile():
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    for k, c, lab in [(0.5, BLUE, "0.5"), (0.7, GREEN, "0.7"),
                      (0.9, ORANGE, "0.9"), (0.9999, RED, r"$\to 1$")]:
        ell = 2*k*K(k)                       # period in units of 1/m_pi
        z = np.linspace(0, 2*ell, 2000)
        ax.plot(z/ell, 2/k*dn(z/k, k), color=c, label=r"$k=$"+lab)
    ax.set_xlabel(r"$z/\ell$"); ax.set_ylabel(r"$\partial_z\varphi\,/\,m_\pi$")
    ax.set_xlim(0, 2); ax.set_ylim(0, 5.6)
    ax.legend(frameon=False, ncol=4, loc="upper center", handlelength=1.4, columnspacing=1.0, fontsize=12)
    ax.text(0.03, 0.06, r"$\ell=2kK(k)/m_\pi$", transform=ax.transAxes, color=GRAY)
    fig.savefig(OUT+"csl_profile.pdf"); plt.close(fig)

# ---------------------------------------------------------------- Fig: k and period vs B
def fig_kofB():
    b = np.linspace(1.0001, 10, 600)         # B/B_CSL = E(k)/k
    ks = np.array([brentq(lambda k: E(k)/k - bb, 1e-9, 1-1e-12) for bb in b])
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.2, 3.9))
    a1.plot(b, ks, color=BLUE)
    a1.set_xlabel(r"$B/B_{\rm CSL}$"); a1.set_ylabel(r"elliptic modulus $k$")
    a1.set_xlim(0, 10); a1.set_ylim(0, 1.05)
    a1.axvline(1, color=GRAY, ls=":")
    a1.text(1.25, 0.12, r"$B_{\rm CSL}$", color=GRAY)
    a2.plot(b, 2*ks*K(ks), color=RED)
    a2.set_xlabel(r"$B/B_{\rm CSL}$"); a2.set_ylabel(r"$m_\pi\ell$")
    a2.set_xlim(0, 10); a2.set_ylim(0, 8.5)
    a2.axvline(1, color=GRAY, ls=":")
    a2.text(4.4, 5.2, "dilute walls", color=GRAY, fontsize=13)
    a2.annotate("", xy=(1.4, 7.4), xytext=(4.2, 5.6),
                arrowprops=dict(arrowstyle="->", color=GRAY))
    a2.text(6.0, 2.2, "sinusoidal", color=GRAY, fontsize=13)
    fig.savefig(OUT+"csl_kofB.pdf"); plt.close(fig)

# ---------------------------------------------------------------- Fig: phase diagram
def fig_phase():
    fpi, mpi = 0.092, 0.140                    # GeV
    Bcsl = lambda mu: 16*np.pi*fpi**2*mpi/mu   # GeV^2
    def Bbec(mu):
        # min omega^2_{n=0} = B - (m^2/k^2)(2 - k^2 + 2 sqrt(1-k^2+k^4)) = 0,
        # with k fixed by E(k)/k = mu B /(16 pi m f^2)
        def g(B):
            r = mu*B/(16*np.pi*mpi*fpi**2)
            if r <= 1: return -1.0
            k = brentq(lambda kk: E(kk)/kk - r, 1e-12, 1-1e-14)
            return B - mpi**2/k**2*(2 - k**2 + 2*np.sqrt(1-k**2+k**4))
        lo = Bcsl(mu)*1.000001
        hi = 5.0
        if g(lo) < 0 or g(hi) > 0: return np.nan
        return brentq(g, lo, hi)
    mu = np.linspace(180, 800, 400)/1000.
    bc = Bcsl(mu)
    bb = np.array([Bbec(m) for m in mu])
    bchi = 16*np.pi**4*fpi**4/mu**2
    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    ax.fill_between(mu*1000, bc, np.minimum(bb, 0.30), where=(bb > bc),
                    color=BLUE, alpha=0.16, lw=0)
    ax.plot(mu*1000, bc, color=BLUE, ls="--", label=r"$B_{\rm CSL}=16\pi f_\pi^2m_\pi/\mu$")
    ax.plot(mu*1000, bb, color=RED, label=r"$B_{\rm BEC}$  (charged pion BEC)")
    ax.plot(mu*1000, bchi, color=RED, ls=":", lw=1.4, label=r"$B_{\rm BEC}$, chiral limit")
    ax.set_xlim(150, 800); ax.set_ylim(0, 0.30)
    ax.set_xlabel(r"$\mu$ [MeV]"); ax.set_ylabel(r"$B$ [GeV$^2$]")
    ax.text(560, 0.125, "CSL", color=BLUE, fontsize=20)
    ax.text(205, 0.055, "QCD vacuum", color=GRAY, fontsize=15)
    ax.text(520, 0.255, r"$\pi^\pm$ BEC", color=RED, fontsize=16)
    ax.text(232, 0.185, r"$B_{\rm CSL}=16\pi f_\pi^2m_\pi/\mu$", color=BLUE, fontsize=13, rotation=-42)
    ax.text(690, 0.225, r"$B_{\rm BEC}$", color=RED, fontsize=13, rotation=-40)
    sec = ax.secondary_yaxis('right', functions=(lambda x: x*1.7e20, lambda x: x/1.7e20))
    sec.set_ylabel(r"$B$ [G]", fontsize=13)
    sec.set_yticks([0, 1e19, 2e19, 3e19, 4e19, 5e19])
    sec.set_yticklabels(["0", r"$10^{19}$", r"$2$", r"$3$", r"$4$", r"$5\times10^{19}$"],
                        fontsize=10)
    fig.savefig(OUT+"csl_phase.pdf"); plt.close(fig)

# ------------------------------------------------- Fig: finite volume solution families
def fig_finitevolume():
    """Gradient profiles obeying the natural BC d_z phi = Hbar at z = +-L/2."""
    Lb = 5.0
    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    # type (14): z0 = 0, gradient maximal at centre
    for Hb, c in [(2.0, RED)]:
        f = lambda k: 2/k*dn(Lb/(2*k), k) - Hb
        ks = np.linspace(0.02, 1-1e-9, 6000); v = [f(k) for k in ks]
        root = None
        for i in range(len(ks)-1):
            if v[i]*v[i+1] < 0: root = brentq(f, ks[i], ks[i+1]); break
        if root is None: continue
        z = np.linspace(-Lb/2, Lb/2, 800)
        ax.plot(z, 2/root*dn(z/root, root), color=c,
                label=r"$\bar H=%.2f$, $k=%.3f$ (type I)" % (Hb, root))
    # type (15): z0 = kK(k), gradient minimal at centre
    for Hb, c in [(0.9, BLUE)]:
        def f(k):
            return 2/k*dn((Lb/2 - k*K(k))/k, k) - Hb
        ks = np.linspace(0.02, 1-1e-9, 6000); v = [f(k) for k in ks]
        root = None
        for i in range(len(ks)-1):
            if v[i]*v[i+1] < 0: root = brentq(f, ks[i], ks[i+1]); break
        if root is not None:
            z = np.linspace(-Lb/2, Lb/2, 800)
            ax.plot(z, 2/root*dn((z - root*K(root))/root, root), color=c,
                    label=r"$\bar H=%.2f$, $k=%.3f$ (type II)" % (Hb, root))
        ax.axhline(Hb, color=BLUE, ls=":", lw=1.2)
    ax.axhline(2.0, color=RED, ls=":", lw=1.2)
    ax.set_xlabel(r"$\bar z=m_\pi z$"); ax.set_ylabel(r"$\partial_{\bar z}\varphi$")
    ax.set_xlim(-Lb/2, Lb/2); ax.set_ylim(0, 2.9)
    ax.legend(frameon=False, loc="lower center", fontsize=11, ncol=1)
    ax.text(0.02, 0.93, r"$\bar L=5$;  dotted lines: $\partial_{\bar z}\varphi=\bar H$ (natural BC)",
            transform=ax.transAxes, color=GRAY, fontsize=12)
    fig.savefig(OUT+"csl_finitevol.pdf"); plt.close(fig)

# ---------------------------------------------------------------- Fig: domain-wall energy
def fig_dwenergy():
    n = np.arange(1, 400, 2)[:, None]
    a = np.linspace(0.02, 5, 500)[None, :]
    ratio = (16/a*np.sum(np.tanh(np.pi*n*a/2)/(np.pi*n)**3, axis=0))[0]
    a = a[0]
    zeta3 = 1.2020569
    lin = 1 - 14*a*zeta3/np.pi**3
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.plot(a, ratio, color=BLUE, label="exact")
    ax.plot(a, lin, color=RED, ls="--",
            label=r"$1-14\alpha\zeta(3)/\pi^3$")
    ax.plot(a, 14*zeta3/(np.pi**3*a), color=GREEN, ls=":",
            label=r"$14\zeta(3)/(\pi^3\alpha)$")
    ax.set_xlim(0, 5); ax.set_ylim(0, 1.05)
    ax.set_xlabel(r"aspect ratio $\alpha=L_z/L_x$")
    ax.set_ylabel(r"$\langle\mathcal{E}_{\rm DW}\rangle/\mathcal{E}_0$")
    ax.plot([1], [0.5], "o", color="k", ms=7, zorder=5)
    ax.annotate(r"square domain: exactly $1/2$", xy=(1, 0.5), xytext=(1.5, 0.72),
                arrowprops=dict(arrowstyle="->", color="k"), fontsize=12)
    ax.legend(frameon=False, loc="lower left", fontsize=12)
    fig.savefig(OUT+"csl_dwenergy.pdf"); plt.close(fig)

# ---------------------------------------------------------------- Fig: three field geometries
def fig_geometries():
    fig, axs = plt.subplots(1, 3, figsize=(11.4, 3.1))
    x = np.linspace(-1, 1, 22); z = np.linspace(-1, 1, 22)
    X, Z = np.meshgrid(x, z)
    # (a) uniform
    axs[0].quiver(X[::3, ::3], Z[::3, ::3], 0*X[::3, ::3], 1+0*Z[::3, ::3],
                  color=BLUE, scale=14, width=0.008)
    axs[0].set_title(r"(a) uniform $\boldsymbol{B}$", fontsize=16)
    # (b) domain wall
    U = np.zeros_like(X); W = np.tanh(X/0.12)
    axs[1].quiver(X[::3, ::3], Z[::3, ::3], U[::3, ::3], W[::3, ::3],
                  color=RED, scale=14, width=0.008)
    axs[1].axvline(0, color="k", ls="--", lw=1.2)
    axs[1].set_title(r"(b) domain wall", fontsize=16)
    # (c) tangential / closed loops
    psi = np.cos(np.pi*X/2)*np.cos(np.pi*Z/2)
    U = np.gradient(psi, z, axis=0); W = -np.gradient(psi, x, axis=1)
    axs[2].streamplot(x, z, U, W, color=GREEN, density=0.8, linewidth=1.2,
                      arrowsize=0.9)
    axs[2].set_title(r"(c) closed loops", fontsize=16)
    for a in axs:
        a.set_xlim(-1, 1); a.set_ylim(-1, 1); a.set_xticks([]); a.set_yticks([])
        for s in a.spines.values(): s.set_linewidth(1.6)
    fig.savefig(OUT+"csl_geometries.pdf"); plt.close(fig)

for f in (fig_profile, fig_kofB, fig_phase, fig_finitevolume, fig_dwenergy, fig_geometries):
    f(); print("ok", f.__name__)
