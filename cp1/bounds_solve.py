import sympy as sp
from sympy import Symbol as Sym, sinh, cosh

# Prob 2

B = Sym("B")
M = Sym("M")
N = Sym("N")

D1 = Sym("D1")
L1 = Sym("L1")
X1 = Sym("X1")
S = Sym("S")
D2 = Sym("D2")
L2 = Sym("L2")

a = Sym("a")
b = Sym("b")


ex1 = M*sinh((a+b)/L2) + N*cosh((a+b)/L2)
ex2 = M*sinh((a)/L2) + N*cosh((a)/L2) - B*cosh(a/L1) - S/X1
ex3 = D2*(M*cosh((a)/L2) + N*sinh((a)/L2))/L2 - D1*(B*sinh(a/L1)/L1)

sol = sp.solve((ex1, ex2, ex3), (B, M, N))
print(f"B = {sp.simplify(sol[B])}")
print(f"M = {sp.simplify(sol[M])}")
print(f"N = {sp.simplify(sol[N])}")