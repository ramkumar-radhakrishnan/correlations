# Digitized CSL Lattice Hamiltonian — generalization to arbitrary $n_q$ and $N$

Generalizes the $N=2$ site, one-qubit-per-site derivation in `CSL_quantumcomp.pdf`
(eqs. 1.1–1.21) to **$n_q$ qubits per site** ($Q=2^{n_q}$ field values) and
**$N$ sites**. Every formula below is checked numerically in `verify.py`
(57/57 checks pass, including exact reproduction of eq. 1.21).

---

## 1. Setup

Continuum energy in one spatial dimension:

$$\hat E=f_\pi^2\int dx\left[\frac{\Pi^2}{2f_\pi^4}+\tfrac12(\partial_x\phi)^2+m_\pi^2(1-\cos\phi)-H\,\partial_x\phi\right],
\qquad H\equiv\frac{\mu B}{4\pi^2f_\pi^2}.$$

Discretize on $N$ sites of spacing $a$ and compactify the two $\phi$-derivative
terms exactly as in the note (eq. 1.6),
$\tfrac12(\Delta\phi)^2\to 1-\cos\Delta\phi$ and $\Delta\phi\to\sin\Delta\phi$,
with $\Delta\phi_n=\phi_{n+1}-\phi_n$:

$$\boxed{\;\hat H=t\sum_{n=1}^{N}\hat\Pi_n^{2}
\;+\;J\sum_{\langle n,n+1\rangle}\bigl[\mathbb 1-\cos\Delta\phi_n\bigr]
\;+\;h\sum_{n=1}^{N}\bigl[\mathbb 1-\cos\phi_n\bigr]
\;-\;\kappa\sum_{\langle n,n+1\rangle}\sin\Delta\phi_n\;}$$

$$t=\frac{1}{2f_\pi^2a},\qquad J=\frac{f_\pi^2}{a},\qquad
h=f_\pi^2m_\pi^2a,\qquad \kappa=f_\pi^2H=\frac{\mu B}{4\pi^2}.$$

The bond sum runs over $N-1$ bonds (open BC — the case used in eq. 1.20) or
$N$ bonds (periodic BC). Note that $\kappa$ is $a$-independent: the magnetic
term is topological in origin.

## 2. Digitization: clock and shift operators

$\phi\in[0,2\pi)\to\phi_k=2\pi k/Q$, $k=0,\dots,Q-1$, stored in binary
$k=\sum_{j=0}^{n_q-1}2^jb_j$ on the site register. Define

$$\hat Z_Q\,|k\rangle=\omega^k|k\rangle\;(\,=e^{i\hat\phi}\,),\qquad
\hat X_Q\,|k\rangle=|k{+}1 \bmod Q\rangle,\qquad \omega=e^{2\pi i/Q}.$$

Then **all four terms are exact, for any $Q$**:

| term | operator form |
|---|---|
| $\cos\phi_n$ | $\tfrac12(\hat Z_n+\hat Z_n^\dagger)$ |
| $\sin\phi_n$ | $\tfrac{1}{2i}(\hat Z_n-\hat Z_n^\dagger)$ |
| $\cos\Delta\phi_n$ | $\tfrac12\bigl(\hat Z_{n+1}\hat Z_n^\dagger+\hat Z_{n+1}^\dagger\hat Z_n\bigr)$ |
| $\sin\Delta\phi_n$ | $\tfrac{1}{2i}\bigl(\hat Z_{n+1}\hat Z_n^\dagger-\hat Z_{n+1}^\dagger\hat Z_n\bigr)$ |

For $Q=2$ these collapse to $\cos\phi=Z$, $\cos\Delta\phi=Z_nZ_{n+1}$,
$\sin\Delta\phi=0$ — eqs. (1.10), (1.17), (1.18) of the note.

## 3. The kinetic term for arbitrary $Q$ (new closed form)

$\hat\Pi=-i\partial_\phi$ has integer eigenvalues $l$ in the symmetric window
$l\in\{-Q/2,\dots,Q/2-1\}$ (eq. 1.11 generalized), with
$|l\rangle=Q^{-1/2}\sum_k\omega^{kl}|k\rangle$. Writing
$\hat\Pi^2=\sum_l l^2|l\rangle\langle l|$ in the **position** basis and using
$\hat X_Q|l\rangle=\omega^{-l}|l\rangle$:

$$\boxed{\;\hat\Pi^2=\frac{Q^2+2}{12}\,\mathbb 1
\;+\;\sum_{m=1}^{Q-1}\frac{(-1)^m}{2\sin^2(\pi m/Q)}\;\hat X_Q^{\,m}\;}$$

**Derivation of the coefficients.** With $f(x)=\sum_{l=-Q/2}^{Q/2-1}x^l$ one has
$f(\omega^m)=0$ for $m\not\equiv0$, and $c_m=\frac1Q(x\partial_x)^2f\big|_{x=\omega^m}
=-2(-1)^m\omega^m/(\omega^m-1)^2=(-1)^m/\bigl(2\sin^2(\pi m/Q)\bigr)$.
The diagonal follows from $\sum_l l^2=\tfrac{Q(Q-1)(Q-2)}{12}+\tfrac{Q^2}{4}$,
giving $c_0=(Q^2+2)/12$; equivalently $c_0=-\sum_{m\neq0}c_m$, which is just the
statement that the $l=0$ eigenvalue vanishes.

This is exactly the Fornberg pseudo-spectral second-derivative matrix on a
periodic grid, $\langle k|\hat\Pi^2|k\rangle=Q^2/12+1/6$ — a good consistency
check, since $\hat\Pi^2=-\partial_\phi^2$.

* $Q=2$: $\hat\Pi^2=\tfrac12(\mathbb 1-X)$, so $t\hat\Pi^2=\frac{1}{4f_\pi^2a}(\mathbb 1-X)$ — **eq. (1.15) exactly**.
* $Q=4$: $\hat\Pi^2=\tfrac32\mathbb 1-(\hat X+\hat X^\dagger)+\tfrac12\hat X^2$, spectrum $\{0,1,1,4\}=l^2$ for $l\in\{0,\pm1,-2\}$.

## 4. Qubit encoding

$$\hat Z_Q=\bigotimes_{j=0}^{n_q-1}\begin{pmatrix}1&0\\0&\omega^{2^j}\end{pmatrix}_{\!j}
= i\,e^{-i\pi/Q}\,Z_{n_q-1}\prod_{j=0}^{n_q-2}
\Bigl[\cos\tfrac{\pi 2^j}{Q}\,\mathbb 1-i\sin\tfrac{\pi 2^j}{Q}\,Z_j\Bigr]$$

Expanding the product over subsets $S\subseteq\{0,\dots,n_q-2\}$ gives the
**closed-form Pauli series** (both have exactly $Q/2$ strings, each containing
$Z_{n_q-1}$, because $\phi\to\phi+\pi$ flips the MSB and flips both functions):

$$\cos\hat\phi=\sum_{S}\sin\!\Bigl(\tfrac{\pi}{Q}+|S|\tfrac{\pi}{2}\Bigr)A_S\;Z_{n_q-1}\!\!\prod_{j\in S}\!Z_j,
\qquad
\sin\hat\phi=\sum_{S}\cos\!\Bigl(\tfrac{\pi}{Q}+|S|\tfrac{\pi}{2}\Bigr)A_S\;Z_{n_q-1}\!\!\prod_{j\in S}\!Z_j,$$

$$A_S=\prod_{j\in S}\sin\tfrac{\pi 2^j}{Q}\prod_{j\notin S}\cos\tfrac{\pi 2^j}{Q}.$$

$\hat X_Q$ is the mod-$Q$ increment (a cascade of multi-controlled NOTs).
Its Pauli expansion has $3^{n_q-1}+1$ strings, so for $n_q\gtrsim3$ one should
instead implement $e^{-i\tau t\hat\Pi^2}$ as $\mathrm{QFT}\to$ diagonal
$e^{-i\tau t l^2}$ phases $\to\mathrm{QFT}^\dagger$ ($O(n_q^2)$ gates).

### $n_q=1$ ($Q=2$), Pauli form
$$\cos\phi_n=Z_n,\qquad \sin\phi_n=0,\qquad \hat\Pi_n^2=\tfrac12(\mathbb 1-X_n).$$

### $n_q=2$ ($Q=4$), Pauli form
Writing $Z_{n,1}$ / $Z_{n,0}$ for the MSB / LSB qubit of site $n$:

$$\cos\phi_n=\tfrac12\bigl(Z_{n,1}+Z_{n,1}Z_{n,0}\bigr),\qquad
\sin\phi_n=\tfrac12\bigl(Z_{n,1}-Z_{n,1}Z_{n,0}\bigr),$$
$$\hat\Pi_n^2=\tfrac32\mathbb 1-X_{n,0}-X_{n,1}X_{n,0}+\tfrac12X_{n,1},$$
$$\cos\Delta\phi_n=\tfrac12\bigl(Z_{n,1}Z_{n+1,1}+Z_{n,1}Z_{n,0}Z_{n+1,1}Z_{n+1,0}\bigr),$$
$$\sin\Delta\phi_n=\tfrac12\,Z_{n,1}Z_{n+1,1}\bigl(Z_{n,0}-Z_{n+1,0}\bigr).$$

## 5. Physics remark: the magnetic term needs $\ge 2$ qubits

For $Q=2$, $\Delta\phi\in\{0,\pm\pi\}$ so $\sin\Delta\phi\equiv0$ and $\kappa$
**drops out of the Hamiltonian identically** — a one-qubit register cannot see
the magnetic field, at any $N$. The chiral soliton lattice first appears at
$n_q=2$. This is confirmed numerically in `spectra.py`: the $n_q=1$ spectra are
exactly $\kappa$-independent, while at $n_q=2$ the order parameter
$\langle\sin\Delta\phi\rangle$ grows with $\kappa$ and the gap closes.

Also note that with the *linear* term $-\kappa\sum\Delta\phi$ the sum telescopes
to a boundary term; the compactified $\sin\Delta\phi$ does not, which is exactly
what lets the soliton lattice form.

---

## 6. Results

Conventions: sites $n=1..N$, open BC. Coefficients are exact rationals times
$t,J,h,\kappa$. Generated and cross-checked by `cases.py` / `verify.py`.

### Reference: $n_q=1$, $N=2$ (the note's case)
$$\hat H=t\!\sum_{n=1}^{2}\!\tfrac12(\mathbb 1-X_n)+J(\mathbb 1-Z_1Z_2)+h\!\sum_{n=1}^{2}\!(\mathbb 1-Z_n)$$
reproduces eq. (1.21) entry by entry with $t=1/(2f_\pi^2a)$ (verified numerically).

### A. One qubit, four sites — $n_q=1$, $N=4$ (4 qubits, $16\times16$)

$$\hat H = 2t\,\mathbb 1-\frac{t}{2}\sum_{n=1}^{4}X_n
\;+\;3J\,\mathbb 1-J\!\!\sum_{n=1}^{3}\!Z_nZ_{n+1}
\;+\;4h\,\mathbb 1-h\sum_{n=1}^{4}Z_n$$

i.e. $\hat H=(2t+3J+4h)\mathbb 1-\tfrac t2\sum_n X_n-J\sum_n Z_nZ_{n+1}-h\sum_n Z_n$:
an **open transverse-field Ising chain** with longitudinal field.
12 Pauli strings; $\kappa$ absent. Diagonal element for configuration
$(k_1..k_4)$, $k_n\in\{0,1\}$:
$$2t+2J\,\#\{n:k_n\neq k_{n+1}\}+2h\,\#\{n:k_n=1\},$$
off-diagonal $-t/2$ between configurations differing in one site.

### B. Two qubits, two sites — $n_q=2$, $N=2$ (4 qubits, $16\times16$)

$$\hat H=3t\,\mathbb 1+t\sum_{n=1}^{2}\Bigl[-X_{n,0}+\tfrac12X_{n,1}-X_{n,1}X_{n,0}\Bigr]$$
$$+\,J\,\mathbb 1-\tfrac J2\Bigl[Z_{1,1}Z_{2,1}+Z_{1,1}Z_{1,0}Z_{2,1}Z_{2,0}\Bigr]$$
$$+\,2h\,\mathbb 1-\tfrac h2\sum_{n=1}^{2}\Bigl[Z_{n,1}+Z_{n,1}Z_{n,0}\Bigr]$$
$$+\,\tfrac\kappa2\Bigl[Z_{1,1}Z_{2,1}Z_{2,0}-Z_{1,1}Z_{1,0}Z_{2,1}\Bigr]$$

15 Pauli strings (identity + 14). This is the smallest system in which the
magnetic term is non-trivial. Diagonal in $|k_1k_2\rangle$:
$3t+J(1-\cos\tfrac{\pi(k_2-k_1)}{2})+h\sum_n(1-\cos\tfrac{\pi k_n}{2})-\kappa\sin\tfrac{\pi(k_2-k_1)}{2}$;
off-diagonals $-t$ for $\Delta k=\pm1$ and $+t/2$ for $\Delta k=\pm2$ on one site.

### C. Two qubits, four sites — $n_q=2$, $N=4$ (8 qubits, $256\times256$)

$$\hat H=6t\,\mathbb 1+t\sum_{n=1}^{4}\Bigl[-X_{n,0}+\tfrac12X_{n,1}-X_{n,1}X_{n,0}\Bigr]$$
$$+\,3J\,\mathbb 1-\tfrac J2\sum_{n=1}^{3}\Bigl[Z_{n,1}Z_{n+1,1}+Z_{n,1}Z_{n,0}Z_{n+1,1}Z_{n+1,0}\Bigr]$$
$$+\,4h\,\mathbb 1-\tfrac h2\sum_{n=1}^{4}\Bigl[Z_{n,1}+Z_{n,1}Z_{n,0}\Bigr]$$
$$-\,\tfrac\kappa2\sum_{n=1}^{3}Z_{n,1}Z_{n+1,1}\Bigl[Z_{n,0}-Z_{n+1,0}\Bigr]$$

33 Pauli strings (identity + 32): 12 kinetic (1-/2-local $X$), 6 gradient
(2-/4-local $Z$), 8 mass (1-/2-local $Z$), 6 magnetic (3-local $Z$).
Max Pauli weight 4; all terms commute except the kinetic $X$'s, so a Trotter
step needs 2 layers.

### Resource scaling (general $n_q$, $N$; open BC)

| sector | non-identity Pauli strings | max weight | $n_q{=}1$ | $n_q{=}2$ | $n_q{=}3$ |
|---|---|---|---|---|---|
| kinetic ($X$) | $N\cdot 3^{\,n_q-1}$ | $n_q$ | $N$ | $3N$ | $9N$ |
| mass ($Z$) | $N\cdot Q/2$ | $n_q$ | $N$ | $2N$ | $4N$ |
| gradient ($Z$) | $(N{-}1)\,2^{\,2n_q-3}$ | $2n_q$ | $N{-}1$ | $2(N{-}1)$ | $8(N{-}1)$ |
| magnetic ($Z$) | $(N{-}1)\,2^{\,2n_q-3}$ | $2n_q{-}1$ | $0$ | $2(N{-}1)$ | $8(N{-}1)$ |

(The $n_q=1$ column is the degenerate case: $\sin\phi\equiv0$ collapses the
gradient term to a single string per bond and kills the magnetic term.)
The counts for the two bond sectors follow because $\cos\phi$ and $\sin\phi$
each carry $Q/2$ strings, so $\cos\Delta\phi=C_{n+1}C_n+S_{n+1}S_n$ starts from
$2\cdot(Q/2)^2$ products of which exactly half cancel.

Total qubits $N n_q$, Hilbert dimension $2^{N n_q}=Q^N$. The whole kinetic
sector mutually commutes and the whole potential sector ($J,h,\kappa$) is
diagonal, so one Trotter step is 2 layers.

---

## 7. Numerics ($f_\pi=a=m_\pi=1\Rightarrow t=\tfrac12, J=h=1$, open BC)

Ground-state energy and CSL order parameter $\langle\sin\Delta\phi\rangle$:

| $\kappa$ | $n_q{=}1,N{=}4$: $E_0$ | $\langle s\rangle$ | $n_q{=}2,N{=}2$: $E_0$ | $\langle s\rangle$ | $n_q{=}2,N{=}4$: $E_0$ | $\langle s\rangle$ |
|---|---|---|---|---|---|---|
| 0.0 | 0.9480 | 0 | 1.0807 | 0.0000 | 2.2809 | 0.0000 |
| 1.0 | 0.9480 | 0 | 1.0016 | 0.1760 | 2.1834 | 0.0776 |
| 2.0 | 0.9480 | 0 | 0.6666 | 0.5192 | 1.6014 | 0.3616 |
| 3.0 | 0.9480 | 0 | −0.0104 | 0.7967 | 0.0303 | 0.6667 |

The $n_q=1$ column is flat in $\kappa$, as it must be.

## Files
* `csl_lattice.py` — general construction for any $(n_q,N)$, BC, couplings.
* `verify.py` — 57 checks incl. exact reproduction of eq. (1.21).
* `cases.py` — Pauli-string Hamiltonians for any case.
* `matrices.py` — explicit symbolic matrices for small cases.
* `spectra.py` — spectra and $\langle\sin\Delta\phi\rangle$.
