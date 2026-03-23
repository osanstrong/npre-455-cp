import numpy as np
from numpy import cosh, sinh
import matplotlib.pyplot as plt

import scipy.stats as stat

# Finite difference solution 
D1 = 1 # cm
L1 = 10
X1 = 0.01 # 1/cm
S = 1e12 # / cm3/s
D2 = 1.7
L2 = 5
X2 = 0.068

a = 10 # cm
b = 20

NL = [5, 10, 20, 40]

def fin_diff_sol(N):
    h = (a+b)/N
    A_mat = np.eye(N+1)

    bound = int(N*a/(a+b))+1
    

    D = np.zeros(N+1)
    D[:bound] = D1
    D[bound:] = D2
    L = np.zeros(N+1)
    L[:bound] = L1
    L[bound:] = L2
    S_vals = np.zeros(N+1)
    S_vals[:bound] = S

    X = np.zeros(N+1)
    X[:bound] = X1
    X[bound:] = X2

    b_mat = np.zeros(N+1)
    b_mat[1:N] = np.array(S_vals[1:N])

    A_mat[0,0:3] = [-3, 4, -1]
    # A_mat[0,0:4] = [-11, 18, -9, 2]
    for i in range(1,N):
        lD = D[i-1]
        iD = D[i]
        rD = D[i+1]
        iL = L[i]
        lfrac = lD/(iD+lD)
        rfrac = rD/(iD+rD)
        print(F"D: {lD}, {iD}, {rD}")
        hD = iD / (0.5 * h**2)
        # A_mat[i,i-1:i+2] = np.array([
        #     -hD*lfrac,
        #     hD*(lfrac+rfrac) + X[i],
        #     -hD*rfrac
        # ])

        vals = (D[i]/(h**2/2)) * np.array([
            -D[i-1]/(D[i]+D[i-1]),
            D[i-1]/(D[i]+D[i-1]) + D[i+1]/(D[i]+D[i+1]),
            -D[i+1]/(D[i]+D[i+1])
        ]) + np.array([
            0, X[i], 0
        ])
        print(f"Vals: {vals}")
        A_mat[i,i-1:i+2] = vals
        # A_mat[i,i-1:i+2] = (D[i]/(h**2/2)) * np.array([
        #     -D[i-1]/(D[i]+D[i-1]),
        #     D[i-1]/(D[i]+D[i-1]) + D[i+1]/(D[i]+D[i+1]),
        #     -hD*D[i+1]/(D[i]+D[i+1])
        # ]) + np.array([
        #     0, X[i], 0
        # ])
        # A_mat[i,i-1:i+2] = np.array([
        #     -iD/(h**2),
        #     2*iD/(h**2) + X[i],
        #     -iD/(h**2)
        # ])
    # A_mat[N,N-1:N+1] = [
    #     -0.5*D[N]/h,
    #     0.25 + 0.5*D[N]/h
    # ]
    A_mat[N,N] = 1
    print("_"*10)
    print(f"N={N}, bound={bound}")
    print(f"A: {A_mat}")
    print(f"b: {b_mat}")
    print(f"S: {S_vals}")
    phi = np.linalg.solve(A_mat, b_mat)
    # phi = np.invert(A_mat) * b_mat
    return phi

for N in NL:
    x = np.linspace(0, a+b, N+1)
    phi = fin_diff_sol(N)
    print(x)
    print(phi)
    plt.plot(x, phi, label=f"Numerical, N = {N}")


# def gen_fin_diff(D_vals:np.ndarray, L_vals:np.ndarray, S_left):
#     mat = np.zeros([D_vals.size]*2)

#     pass


# Plot vs analytic

B = -D2*L1*S*cosh(b/L2)/(X1*(D1*L2*sinh(a/L1)*sinh(b/L2) + D2*L1*cosh(a/L1)*cosh(b/L2)))
M = -D1*L2*S*sinh(a/L1)*cosh((a + b)/L2)/(X1*(D1*L2*sinh(a/L1)*sinh(b/L2) + D2*L1*cosh(a/L1)*cosh(b/L2)))
N = D1*L2*S*sinh(a/L1)*sinh((a + b)/L2)/(X1*(D1*L2*sinh(a/L1)*sinh(b/L2) + D2*L1*cosh(a/L1)*cosh(b/L2)))

N_PER_CM = 10
xa = np.linspace(0, a+b, (a+b)*N_PER_CM+1)
phi = np.zeros(xa.shape)
xa1 = xa[:a*N_PER_CM]
xa2 = xa[a*N_PER_CM:]
# phi[:a*N_PER_CM] = B*cosh(xa1/L1) + S/X1
# phi[a*N_PER_CM:] = M*sinh(xa2/L2) + N*cosh(xa2/L2)
def phi_a(xv):
    if xv <= a:
        return B*cosh(xv/L1) + S/X1
    else:
        return M*sinh(xv/L2) + N*cosh(xv/L2)
def D_a(xv):
    if xv <= a:
        return D1
    else: 
        return D2
def curr(xv, h):
    iD = D_a(xv)
    lD = D_a(xv-h)
    rD = D_a(xv+h)
    lfrac = lD / (lD+iD)
    rfrac = rD / (rD+iD)
    hD = iD / (0.5 * h**2)

    iP = phi_a(xv)
    lP = phi_a(xv-h)
    rP = phi_a(xv+h)
    return -hD * (lfrac*(lP-iP) + rfrac*(rP-iP))
phi = np.array([phi_a(xv) for xv in xa])

def vbar(x, label, c=None, style="dashed"):
    [ymin, ymax] = plt.ylim()
    plt.vlines(x, ymin, ymax, label=label, colors=c, linestyles=style)
    plt.ylim(ymin, ymax)

plt.plot(xa, phi, label="Analytical")

# for h in [10, 6, 1, 0.5]:
#     # diff_simp = np.zeros(len(xa))
#     # diff_simp[:len(diff_simp)-1] = [phi_a(xv)]
#     xvals = np.arange(0, a+b+h, h)
#     cvals = np.array([curr(x, h) for x in xvals])
#     print(cvals)
#     plt.plot(xvals, cvals, label=f"-dDdPhi, h={h} cm")

vbar(a, label=f"a ({a} cm)", c='C1')
vbar(a+b, label=f"a+b (b = {b} cm)", c="C2")
plt.xlabel("x (cm)")
plt.ylabel("Flux (neutrons/cm2/s)")
plt.legend()
plt.show()


err_N = [25, 50, 100, 200, 400, 800, 1600]
err_xv = [np.linspace(0, a+b, Nv+1) for Nv in err_N]
err_sols = [fin_diff_sol(Nv) for Nv in err_N]
tar_sols = [np.array([phi_a(xv) for xv in xvals]) for xvals in err_xv]
def err(actual, target):
    diff = actual-target
    return np.linalg.norm(diff) / np.linalg.norm(target)

err_errs = [err(err_sols[i], tar_sols[i]) for i in range(len(err_N))]

res = stat.linregress(np.log(err_N), np.log(err_errs))
print(res)
# LinregressResult(slope=-1.0091264964866835, intercept=-0.2527060789232243, rvalue=-0.9994067245685468, pvalue=1.64610185963301e-08, stderr=0.015552388491767849, intercept_stderr=0.08517539132364403)
# Observed order: -1.01
plt.plot(err_N, err_errs, label="Finite Difference solutions")
plt.xlabel("Number of slots N")
plt.ylabel("Normalized error")
plt.legend()
plt.show()

