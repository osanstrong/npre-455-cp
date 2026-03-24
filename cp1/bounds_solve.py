# For Problem 2

import sympy as sp
from sympy import sinh, cosh

C1, C2 = sp.symbols("B M")

D1, L1, X1, S, D2, L2 = sp.symbols("D1 L1 X1 S D2 L2")

a, b, x = sp.symbols("a b x")

phi1 = C1*cosh(x/L1) + S/X1
phi2 = C2*sinh((x-(a+b))/L2)
dphi1 = sp.derive_by_array(phi1, x)
dphi2 = sp.derive_by_array(phi2, x)

e1 = sp.simplify((phi1 - phi2).subs(x, a))
e2 = sp.simplify((D1*dphi1 - D2*dphi2).subs(x, a))

sol = sp.solve((e1, e2), (C1, C2))
print(f"C1 = {sp.simplify(sol[C1])}")
print(f"C2 = {sp.simplify(sol[C2])}")

