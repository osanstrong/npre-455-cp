import sympy as sp
from sympy import Symbol as Sym, sinh, cosh, sin, cos

import numpy as np
import matplotlib.pyplot as plt


D1 = Sym("D_1")
L1 = Sym("L_1")
Xf1 = Sym("X_{f1}")
X1 = Sym("X1")
S = Sym("S")
D2 = Sym("D_2")
L2 = Sym("L_2")

a = Sym("a")
b = Sym("b")

# Prob 3, B is buckling
M = Sym("M")
N = Sym("N")
P = Sym("P")
Q = Sym("Q")
# Q = P*sp.tanh((a+b)/L2)
B = Sym("B")

x = Sym("x")
# phi1 = M*sin(B*x) + N*cos(B*x)
# dphi1 = sp.derive_by_array(phi1, x)
# phi2 = P*sinh(x/L2) + Q*cosh(x/L2)
# dphi2 = sp.derive_by_array(phi2, x)

# # Expressions
# # conserv of flux at a
# e0 = -D1*dphi1.subs(x, 0)
# e1 = phi1.subs(x, a) - phi2.subs(x, a)
# # conserv of current at a
# e2 = D1*dphi1.subs(x, a) - D2*dphi2.subs(x, a)
# # vacuum at a+b
# e3 = phi2.subs(x, (a+b))
# e4 = Xf1*N*sin(B*a)/B - 1

# exprs = (e0, e1, e2, e3, e4)

# cc_idx = 3 # Which one to use as the criticality condition

# ncc_exprs = [exprs[i] for i in range(len(exprs)) if not i==cc_idx]
# # print(exprs)
# # print(ncc_exprs)

# PQ_sol = sp.solve(ncc_exprs, (M, P, Q, N))
# PN = sp.simplify(PQ_sol[P])
# QN = sp.simplify(PQ_sol[Q])
# MN = sp.simplify(PQ_sol[M])
# NN = sp.simplify(PQ_sol[N])
# print(PQ_sol)
# # MN, PN, QN, NN = [sp.simplify(v) for v in PQ_sol]
# print(f"P: {PN}")
# print(f"Q: {QN}")
# print(f"M: {MN}")
# print(f"N: {NN}")
# cc_expr = exprs[cc_idx].subs(P, PN).subs(Q, QN).subs(M, MN)
# cc_expr = sp.simplify(cc_expr)
# print(cc_expr)


D1_val = 0.79 # cm
Xa1 = 0.066 # 1/cm
Xf1_val = 0.02787
nu = 2.4
L1_val = (D1_val/Xa1)**0.5

D2_val = 1
Xa2 = 0.000709 # still 1/cm
L2_val = (D2_val/Xa2)**0.5

# # print(dphi1)
# subs_b = 

# cc_b = cc_expr.subs([(k, subs_b[k]) for k in subs_b])
# b = sp.nsolve(cc_b, [0, 100], solver="bisect", verify=False)
# print(b)


phi1 = N*cos(B*x)
dphi1 = sp.derive_by_array(phi1, x)
iphi1 = (N/B)*sin(B*x)
phi2 = P*sinh((x-(a+b))/L2)
dphi2 = sp.derive_by_array(phi2, x)
iphi2 = L2*P*cosh((x-(a+b))/L2)

e1 = phi1.subs(x, a) - phi2.subs(x, a)
e2 = -D1*dphi1.subs(x, a) + D2*dphi2.subs(x, a)

cc_expr = e2
ncc_expr = e2 if cc_expr == e1 else e1

PN = sp.simplify(sp.solve(cc_expr, P)[0])
print(f"P = {PN}")
print(ncc_expr.subs(P, PN))

B_val = (nu*Xf1_val - Xa1)/D1_val
flux_subs = lambda a_val, b_val: { # known that critical, solve for b given a
D1 : D1_val, # cm
D2 : D2_val,
Xf1 : Xf1_val,
B : B_val,
a : a_val,
b : b_val,
L1 : L1_val,
L2 : L2_val
}

Rf = sp.simplify(Xf1*(iphi1.subs(x, a) - iphi1.subs(x, 0)))
Rf_o_N = sp.simplify(Rf/N)
print(f"Rf: {Rf_o_N}")
ab_vals = [(10, 50), (25, 50), (50, 50)]
for a_val, b_val in ab_vals:
    N_val = 1 / (Rf_o_N.evalf(subs=flux_subs(a_val, b_val)))
    print(f"a={a_val}, b={b_val}, N: {N_val}")

