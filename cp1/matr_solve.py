import sympy as sp
from sympy import Symbol as Sym, Matrix

h = Sym('h')
A = Matrix([
    [1, 1, 1, 1],
    [0, h, 2*h, 3*h],
    [0, h**2/2, 4*h**2/2, 9*h**2/2],
    [0, h**3/6, 8*h**3/6, 27*h**3/6],
])
b = Matrix([
    [0],[0],[1],[0]
])
res = A**-1 * b
print(f"x: {sp.simplify(res)}")
print(f"")