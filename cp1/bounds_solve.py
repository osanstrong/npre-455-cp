import sympy as sp
from sympy import Symbol as Sym, sinh, cosh, sin, cos

import numpy as np
import matplotlib.pyplot as plt

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


# ex1 = M*sinh((a+b)/L2) + N*cosh((a+b)/L2)
# ex2 = M*sinh((a)/L2) + N*cosh((a)/L2) - B*cosh(a/L1) - S/X1
# ex3 = D2*(M*cosh((a)/L2) + N*sinh((a)/L2))/L2 - D1*(B*sinh(a/L1)/L1)

# sol = sp.solve((ex1, ex2, ex3), (B, M, N))
# print(f"B = {sp.simplify(sol[B])}")
# print(f"M = {sp.simplify(sol[M])}")
# print(f"N = {sp.simplify(sol[N])}")

# Prob 3, B1 is buckling
N = Sym("N")
P = Sym("P")
Q = Sym("Q")
# Q = P*sp.tanh((a+b)/L2)
B1 = Sym("B1")
# x = Sym("x") # Alias for B1 to avoid making an explicit solution for B1?


e1 = -N*cos(B1*a) + P*sinh(a/L2) + Q*cosh(a/L2)
e2 = D1*(-N*B1*sin(B1*a)) - (D2/L2)*(P*cosh(a/L2) + Q*sinh(a/L2))
e3 = P*sinh((a+b)/L2) + Q*cosh((a+b)/L2)

# A = sp.Matrix([
#     [-cos(B1*a), sinh(a/L2), cosh(a/L2)],
#     []
# ])

# Target: solve for P and Q in terms of N, to solve for N using integral flux


# Target: Solve for N, P, Q in terms of B -> write the fourth equation entirely in terms of B (or one of the og three?)
solB = sp.solve((e1, e3),  (P, Q))
print(f"P = {sp.simplify(solB[P])}")
print(f"Q = {sp.simplify(solB[Q])}")
# print(f"N = {sp.simplify(solB[N])}")
# solB = sp.solve((e1, e2), (N, P))
# print(f"P = {sp.simplify(solB[P])}")
# print(f"N = {sp.simplify(solB[N])}")

# Values that plug into third equation
P = -N*(B1*D1*L2*sin(B1*a)*cosh(a/L2) + D2*cos(B1*a)*sinh(a/L2))/D2
Q = N*(B1*D1*L2*sin(B1*a)*sinh(a/L2) + D2*cos(B1*a)*cosh(a/L2))/D2

# Solve this transcendental hellpit to find B1 for specific geometry
# 0 = -sinh((a+b)/L2)*(B1*D1*L2*sin(B1*a)*cosh(a/L2) + D2*cos(B1*a)*sinh(a/L2)) + cosh((a+b)/L2)*(B1*D1*L2*sin(B1*a)*sinh(a/L2) + D2*cos(B1*a)*cosh(a/L2))
# Or with known B1 (i.e. critical so B1=(nu*Xf1-Xa1)/D1), just... evaluate it to solve for b

# Set hard evals
D1 = 0.79 # cm
Xa1 = 0.066 # 1/cm
Xf1 = 0.02787
nu = 2.4
L1 = (D1/Xa1)**0.5

D2 = 1
Xa2 = 0.000709 # still 1/cm
L2 = (D2/Xa2)**0.5 


# Solve for b that makes critical for given a
a_vals = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
a_vals = [10, 25, 50]
B1 = (nu*Xf1 - Xa1)/D1 # for k to == 1
for a in a_vals:
    # eq_b = -sinh((a+b)/L2)*(B1*D1*L2*sin(B1*a)*cosh(a/L2) + D2*cos(B1*a)*sinh(a/L2)) + cosh((a+b)/L2)*(B1*D1*L2*sin(B1*a)*sinh(a/L2) + D2*cos(B1*a)*cosh(a/L2))

    # P = -N*(B1*D1*L2*sin(B1*a)*cosh(a/L2) + D2*cos(B1*a)*sinh(a/L2))/D2
    # Q = N*(B1*D1*L2*sin(B1*a)*sinh(a/L2) + D2*cos(B1*a)*cosh(a/L2))/D2
    N = 1
    b_vals = np.linspace(0, 100, 1000)

    # Omitting eq3
    # P = -N*(B1*D1*L2*np.sin(B1*a)*np.cosh(a/L2) + D2*np.cos(B1*a)*np.sinh(a/L2))/D2
    # print(P)
    # Q = N*(B1*D1*L2*np.sin(B1*a)*np.sinh(a/L2) + D2*np.cos(B1*a)*np.cosh(a/L2))/D2

    # Omitting eq2
    P = -N*np.cos(B1*a)*np.cosh((a + b_vals)/L2)/np.sinh(b_vals/L2)
    Q = N*np.cos(B1*a)*np.sinh((a + b_vals)/L2)/np.sinh(b_vals/L2)

    # Omitting eq1
    # P = -B1*D1*L2*N*np.sin(B1*a)*np.cosh((a + b_vals)/L2)/(D2*np.cosh(b_vals/L2))
    # Q = B1*D1*L2*N*np.sin(B1*a)*np.sinh((a + b_vals)/L2)/(D2*np.cosh(b_vals/L2))

    e1 = -N*np.cos(B1*a) + P*np.sinh(a/L2) + Q*np.cosh(a/L2)
    e2 = D1*(-N*B1*np.sin(B1*a)) - (D2/L2)*(P*np.cosh(a/L2) + Q*np.sinh(a/L2))
    e3 = P*np.sinh((a+b_vals)/L2) + Q*np.cosh((a+b_vals)/L2)
    
    eq_b = P*np.sinh((a+b_vals)/L2) + Q*np.cosh((a+b_vals)/L2)
    # plt.plot(b_vals, e2, label=f"b for a={a} cm")
    # plt.xlabel("b (cm)")
    # plt.ylabel("validity (0 is valid)")
    # plt.legend()
    # plt.show()
    # try:
    #     b = sp.nsolve(eq_b, (0,100), solver='bisect', verify=False)
    #     print(f"a={a}->b={b}")
    # except:
    #     print(f"a={a}->Cannot crit :(")

a = 50
b = 50
B1 = Sym("B1")
# eq_b = -sinh((a+b)/L2)*(B1*D1*L2*sin(B1*a)*cosh(a/L2) + D2*cos(B1*a)*sinh(a/L2)) + cosh((a+b)/L2)*(B1*D1*L2*sin(B1*a)*sinh(a/L2) + D2*cos(B1*a)*cosh(a/L2))

P = -N*(B1*D1*L2*sin(B1*a)*cosh(a/L2) + D2*cos(B1*a)*sinh(a/L2))/D2
Q = N*(B1*D1*L2*sin(B1*a)*sinh(a/L2) + D2*cos(B1*a)*cosh(a/L2))/D2
# eq_b = P*sinh((a+b)/L2) + Q*cosh((a+b)/L2)
# e1 = -N*cos(B1*a) + P*sinh(a/L2) + Q*cosh(a/L2)
# e2 = D1*(-N*B1*sin(B1*a)) - (D2/L2)*(P*cosh(a/L2) + Q*sinh(a/L2))
e3 = P*sinh((a+b)/L2) + Q*cosh((a+b)/L2)
try:
    B1 = sp.nsolve(e3, (-100,100), solver='bisect', verify=False)
    print(f"a={a}, b={b} -> B1={B1}")
    B12 = B1**2
    iknuXf = B12*D1 + Xa1
    ik = iknuXf/(nu*Xf1)
    k = 1/ik
    print(f"k: {k}")
except:
    print(f"a={a}, b={b} ->Cannot crit :(")

