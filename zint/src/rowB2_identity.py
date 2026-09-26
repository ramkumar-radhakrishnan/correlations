"""Check the quoted identity  f^{abc} f^{def} U^{be}(z) U^{cf}(z) = N_c S^{ad}(z)
and its use in the 2a / 2b split."""
import numpy as np
from scipy.linalg import expm
rng = np.random.default_rng(11)
Nc = 3
l = np.zeros((8,3,3), dtype=complex)
l[0][0,1]=l[0][1,0]=1
l[1][0,1]=-1j; l[1][1,0]=1j
l[2][0,0]=1;  l[2][1,1]=-1
l[3][0,2]=l[3][2,0]=1
l[4][0,2]=-1j; l[4][2,0]=1j
l[5][1,2]=l[5][2,1]=1
l[6][1,2]=-1j; l[6][2,1]=1j
l[7]=np.diag([1,1,-2])/np.sqrt(3)
t = l/2
f = np.real(-2j*np.einsum('aij,bjk,cki->abc',t,t,t) + 2j*np.einsum('bij,ajk,cki->abc',t,t,t))
def adj(V): return np.real(2*np.einsum('aij,jk,bkl,li->ab', t, V, t, V.conj().T))
def randH():
    H = rng.normal(size=(3,3)) + 1j*rng.normal(size=(3,3)); H=(H+H.conj().T)/2
    return H - np.trace(H)*np.eye(3)/3
def randU(s=0.55): return adj(expm(1j*s*randH()))

print("="*78)
print("1.  THE IDENTITY, AND WHAT S^{ad} IS")
print("="*78)
w1=w2=0.0
for _ in range(50):
    U = randU()
    M = np.einsum('abc,def,be,cf->ad', f, f, U, U)
    w1 = max(w1, abs(M - Nc*U).max())          # S = U ?
    w2 = max(w2, abs(M - Nc*U.T).max())        # S = U^T ?
print("   f^{abc} f^{def} U^{be}(z) U^{cf}(z) - N_c U^{ad}(z)   : max dev = %.3e"%w1)
print("   ... - N_c U^{da}(z)  (the transpose)                  : max dev = %.3e"%w2)
print()
print("   So the identity is CORRECT, with  S^{ad}(z) = U^{ad}(z)  --  the adjoint Wilson line")
print("   itself, not a new object.  Derivation in two lines:")
print("     adjointness:      f^{abc} = U^{aa'} U^{bb'} U^{cc'} f^{a'b'c'}")
print("     multiply by U^{ad} and use U^T U = 1 :   f^{def} U^{be} U^{cf} = U^{gd} f^{gbc}")
print("     contract with f^{abc} and use f^{abc} f^{gbc} = N_c delta^{ag} :  = N_c U^{ad} .")
print("   check  f^{abc} f^{gbc} = N_c delta^{ag} : max dev = %.3e"
      %abs(np.einsum('abc,gbc->ag', f, f) - Nc*np.eye(8)).max())
print("   check  adjointness f^{a'b'c'} U^{aa'}U^{bb'}U^{cc'} = f^{abc} : max dev = %.3e"
      %(lambda U: abs(np.einsum('ABC,aA,bB,cC->abc', f, U, U, U) - f).max())(randU()))

print()
print("="*78)
print("2.  IT ONLY WORKS AT COINCIDENT POINTS")
print("="*78)
Ux, Uz = randU(), randU()
M_mixed = np.einsum('abc,def,be,cf->ad', f, f, Ux, Uz)
print("   with U(x) and U(z) at DIFFERENT points:")
print("     ||M - N_c U(x)|| = %.4f    ||M - N_c U(z)|| = %.4f    (for scale ||M|| = %.4f)"
      %(abs(M_mixed-Nc*Ux).max(), abs(M_mixed-Nc*Uz).max(), abs(M_mixed).max()))
print("   -> so the identity is exactly the x -> z statement, and using it IS taking the UV limit.")

print()
print("="*78)
print("3.  MATCHING THE INDEX NAMES OF YOUR ROW")
print("="*78)
print("   identity  (a,b,c) -> (a,d,c) of your first f ,  (d,e,f) -> (f,b,e) of your second:")
print("     f^{adc} f^{fbe} U^{db}(z) U^{ce}(z) = N_c U^{af}(z)")
Uz2 = randU()
print("     check : max dev = %.3e"
      %abs(np.einsum('adc,fbe,db,ce->af', f, f, Uz2, Uz2) - Nc*Uz2).max())
print("   so your 2a colour  -f^{adc} f^{fbe} U^{db}(x) U^{ce}(z) + N_c U^{af}(x)")
print("   goes to  -N_c U^{af}(z) + N_c U^{af}(x)  ->  0  as x -> z .  Consistent.")

print()
print("="*78)
print("4.  THE LINEAR-ORDER FORM OF THE REMAINDER:  it IS N_c/2 [U(x) - U(z)]")
print("="*78)
print("   The EXACT remainder is  -f^{adc} f^{fbe} U^{db}(x) [ U^{ce}(z) - U^{ce}(x) ] .")
print("   For any tangent variation dU of an adjoint Wilson line,")
print("       -f^{adc} f^{fbe} U^{db} dU^{ce}  =  -(N_c/2) dU^{af}   EXACTLY :")
Uc = randU()
w = 0.0
for _ in range(30):
    H = randH()
    dU = (adj(expm(1j*1e-6*H) @ np.linalg.inv(np.eye(3))) if False else
          (adj(expm(1j*1e-6*H)@expm(0j*np.eye(3))) ))
    # build a genuine tangent vector: d/ds adj(exp(i s H) V) at s=0, V the line behind Uc
    eps = 1e-7
    Vb = expm(1j*0.55*randH())
    U0 = adj(Vb); U1 = adj(expm(1j*eps*H) @ Vb)
    dU = (U1 - U0)/eps
    lhs = -np.einsum('adc,fbe,db,ce->af', f, f, U0, dU)
    rhs = -(Nc/2)*dU
    w = max(w, abs(lhs-rhs).max()/max(abs(rhs).max(),1e-30))
print("       max relative deviation over 30 random tangent directions = %.3e"%w)
print()
print("   So to LEADING order in the separation")
print("       C_2a  =  (N_c/2) [ U^{af}(x) - U^{af}(z) ]  +  O((z-x)^2) ,")
print("   a dipole-like difference with strength N_c/2 , not N_c .  Checked against the exact form:")
x0 = np.array([0.3,-0.2]); H1,H2 = randH(), randH()
def Ufield(v,s=0.5): return adj(expm(1j*s*(v[0]*H1+v[1]*H2)))
Uxf = Ufield(x0)
print("   %-10s %16s %18s"%("|z-x|","||C_2a exact||","rel. dev. from N_c/2 form"))
for r in (3e-1,1e-1,1e-2,1e-3,1e-4):
    d = np.array([0.6,0.8]); Uzf = Ufield(x0+r*d)
    ex = -np.einsum('adc,fbe,db,ce->af', f, f, Uxf, Uzf-Uxf)
    ap = (Nc/2)*(Uxf-Uzf)
    print("   %-10.0e %16.6e %18.3e"%(r, abs(ex).max(), abs(ex-ap).max()/abs(ex).max()))
print("   -> the N_c/2 form is the leading term only; the deviation is O(|z-x|) as expected.")
print("      Use the exact f f U [U(z) - U(x)] form in the note, and the N_c/2 form to see at a")
print("      glance that 2a is a dipole of strength N_c/2 and therefore UV finite.")

print("="*78)
print("5.  WHY PUTTING U(x) (NOT U(z)) IN 2b WAS THE RIGHT CHOICE")
print("="*78)
print("   Both splits are exact by the same identity:")
print("     (a)  2b = N_c[U(y) - U(x)] ,  2a = -f f U^{db}(x) [U^{ce}(z) - U^{ce}(x)]   <- yours")
print("     (b)  2b = N_c[U(y) - U(z)] ,  2a = -f f [U^{db}(x) - U^{db}(z)] U^{ce}(z)")
Uy2 = randU(); Uxx = randU(); Uzz = randU()
A = -np.einsum('adc,fbe,db,ce->af', f, f, Uxx, Uzz-Uxx) + Nc*(Uy2-Uxx)
B = -np.einsum('adc,fbe,db,ce->af', f, f, Uxx-Uzz, Uzz) + Nc*(Uy2-Uzz)
orig = np.einsum('adc,ebf,db,ce->af', f, f, Uxx, Uzz) + Nc*Uy2
print("   (a) - original : %.3e      (b) - original : %.3e"
      %(abs(A-orig).max(), abs(B-orig).max()))
print("   But (a) makes 2b's colour z-INDEPENDENT, so 2b is LO-like and its large-|z| logarithm")
print("   has a clean colour structure.  In (b) the colour still depends on z and the infrared")
print("   bookkeeping is messier.  Your choice is the better one.")
