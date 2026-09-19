"""Inline SVG: DGLAP-like resummation factor K^2 = (Q_s^2/k_T^2)^{2 gamma} vs k_T/Q_s.
Log-log frame: the result is a power law, so each curve is a straight line whose
slope is -2 gamma = -alpha_s beta_0^g/(2 pi)."""
import math

B0G = 11.0*3.0/3.0                       # beta_0^g = 11 Nc / 3, Nc = 3
SERIES = [(0.15, "#2a78d6"), (0.25, "#eb6834"), (0.35, "#1baf7a")]

W, H = 640, 340
L, R, T, B = 60, 112, 24, 48
PW, PH = W-L-R, H-T-B
X0, X1 = 0.15, 3.0
Y0, Y1 = 1.0, 12.0

lx0, lx1 = math.log10(X0), math.log10(X1)
ly0, ly1 = math.log10(Y0), math.log10(Y1)
sx = lambda v: L + PW*(math.log10(v)-lx0)/(lx1-lx0)
sy = lambda v: T + PH*(1 - (math.log10(v)-ly0)/(ly1-ly0))

def K2(r, a):
    return 1.0 if r >= 1.0 else (1.0/r**2) ** (2.0*a*B0G/(4.0*math.pi))

o = [f'<svg viewBox="0 0 {W} {H}" width="100%" role="img" '
     f'aria-label="DGLAP-like resummation factor versus transverse momentum, log-log" '
     f'style="max-width:{W}px;display:block;margin:.2em auto .1em;'
     f'font-family:\'Times New Roman\',Times,serif">',
     '<title>DGLAP-like resummation factor for single-inclusive gluon production</title>']

o.append(f'<rect x="{sx(1.0):.1f}" y="{T}" width="{L+PW-sx(1.0):.1f}" height="{PH}" fill="#f1f1ee"/>')
o.append(f'<text x="{(sx(1.0)+L+PW)/2:.1f}" y="{T+15}" text-anchor="middle" '
         f'font-size="11" fill="#6b6a65">no DGLAP log</text>')

for v in [1,1.5,2,3,5,8,12]:
    y = sy(v)
    o.append(f'<line x1="{L}" y1="{y:.1f}" x2="{L+PW}" y2="{y:.1f}" stroke="#dedcd6" stroke-width="1"/>')
    lab = f"{v:g}"
    o.append(f'<text x="{L-9}" y="{y+4:.1f}" text-anchor="end" font-size="12" fill="#52514e">{lab}</text>')
for v,lab in [(0.15,"0.15"),(0.2,"0.2"),(0.3,"0.3"),(0.5,"0.5"),(0.7,"0.7"),(1.0,"1"),(2.0,"2"),(3.0,"3")]:
    x = sx(v)
    o.append(f'<line x1="{x:.1f}" y1="{T+PH}" x2="{x:.1f}" y2="{T+PH+5}" stroke="#b6b4ad" stroke-width="1"/>')
    o.append(f'<text x="{x:.1f}" y="{T+PH+19}" text-anchor="middle" font-size="12" fill="#52514e">{lab}</text>')
o.append(f'<line x1="{L}" y1="{T+PH}" x2="{L+PW}" y2="{T+PH}" stroke="#8d8b85" stroke-width="1"/>')
o.append(f'<line x1="{L}" y1="{T}" x2="{L}" y2="{T+PH}" stroke="#8d8b85" stroke-width="1"/>')
o.append(f'<line x1="{sx(1.0):.1f}" y1="{T}" x2="{sx(1.0):.1f}" y2="{T+PH}" '
         f'stroke="#8d8b85" stroke-width="1" stroke-dasharray="3 3"/>')

for a, col in SERIES:
    pts, r = [], X0
    while r <= X1 + 1e-9:
        pts.append(f"{sx(r):.2f},{sy(max(min(K2(r,a),Y1),Y0)):.2f}")
        r *= 1.0+0.01
    o.append(f'<polyline points="{" ".join(pts)}" fill="none" stroke="{col}" '
             f'stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/>')
    # direct label at the left end of the line (relief rule: visible labels)
    yv = K2(X0,a)
    o.append(f'<circle cx="{sx(X0):.1f}" cy="{sy(yv):.1f}" r="3.2" fill="{col}" '
             f'stroke="#fcfcfb" stroke-width="1.6"/>')

# legend, one row, labels to the right of the plot at each curve's terminal height
for a, col in SERIES:
    yv = K2(X0,a)
    o.append(f'<text x="{L+PW+10:.1f}" y="{sy(yv)+4:.1f}" font-size="12.5" fill="#0b0b0b">'
             f'<tspan fill="{col}" font-size="15">&#9473;</tspan> '
             f'&#945;<tspan font-size="9.5" dy="3">s</tspan>'
             f'<tspan dy="-3"> = {a:.2f}</tspan></text>')
    o.append(f'<line x1="{sx(X0):.1f}" y1="{sy(yv):.1f}" x2="{L+PW+6:.1f}" y2="{sy(yv):.1f}" '
             f'stroke="{col}" stroke-width="1" stroke-dasharray="2 3" opacity=".55"/>')

o.append(f'<text x="{L+PW/2:.1f}" y="{H-8}" text-anchor="middle" font-size="12.5" fill="#0b0b0b">'
         f'k<tspan font-size="9.5" dy="3">T</tspan><tspan dy="-3"> / Q</tspan>'
         f'<tspan font-size="9.5" dy="3">s</tspan></text>')
o.append(f'<text transform="translate(16,{T+PH/2:.1f}) rotate(-90)" text-anchor="middle" '
         f'font-size="12.5" fill="#0b0b0b">resummation factor</text>')
o.append('</svg>')
open("figure_factor.svg","w").write("\n".join(o))

rows = [0.15,0.2,0.3,0.5,0.7,1.0]
print("k_T/Q_s | " + " | ".join(f"as={a:.2f}" for a,_ in SERIES))
for r in rows:
    print(f"  {r:4.2f}  | " + " | ".join(f" {K2(r,a):5.2f}" for a,_ in SERIES))
print()
for a,_ in SERIES:
    print(f"alpha_s={a:.2f}: gamma={a*B0G/(4*math.pi):.4f}  2gamma={2*a*B0G/(4*math.pi):.4f}")
