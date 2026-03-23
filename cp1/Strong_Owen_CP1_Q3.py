import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt 


# Finite difference solution 
D1 = 0.79 # cm
Xa1 = 0.066 # 1/cm
Xf1 = 0.02787
nu = 2.4

D2 = 1
Xa2 = 0.000709 # still 1/cm

# a = 10 # cm
# b = 20


def fin_diff_mats(N, a, b):
    # For given dimensions, with N cells
    h = (a+b)/N
    M_mat = np.eye(N+1)

    bound = int(N*a/(a+b))+1
    

    D = np.zeros(N+1)
    D[:bound] = D1
    D[bound:] = D2

    Xa = np.zeros(N+1)
    Xa[:bound] = Xa1
    Xa[bound:] = Xa2

    Xf = np.zeros(N+1)
    Xf[:bound] = Xf1
    Xf[bound:] = 0 # Yes redundant but here for clarity

    # b_mat = np.zeros(N+1)
    F_mat = np.zeros((N+1, N+1))
    # for i in range(1:N):
    for i in range(1, N):
        F_mat[i,i] = nu*Xf[i]
    # b_mat[1:N] = np.array(S_vals[1:N])

    M_mat[0,0:3] = [-3, 4, -1]
    # M_mat[0,0:4] = [-11, 18, -9, 2]
    for i in range(1,N):
        # lD = D[i-1]
        # iD = D[i]
        # rD = D[i+1]
        # # lfrac = lD/(iD+lD)
        # # rfrac = rD/(iD+rD)
        # # print(F"D: {lD}, {iD}, {rD}")
        # # hD = iD / (0.5 * h**2)

        vals = (D[i]/(h**2/2)) * np.array([
            -D[i-1]/(D[i]+D[i-1]),
            D[i-1]/(D[i]+D[i-1]) + D[i+1]/(D[i]+D[i+1]),
            -D[i+1]/(D[i]+D[i+1])
        ]) + np.array([
            0, Xa[i], 0
        ])
        # print(f"Vals: {vals}")
        M_mat[i,i-1:i+2] = vals
    return M_mat, F_mat
    # A_mat[N,N] = 1
    # print("_"*10)
    # print(f"N={N}, bound={bound}")
    # print(f"A: {A_mat}")
    # print(f"b: {b_mat}")
    # print(f"S: {S_vals}")
    # phi = np.linalg.solve(A_mat, b_mat)
    # # phi = np.invert(A_mat) * b_mat
    # return phi

def err(actual, target):
    diff = actual-target
    return np.linalg.norm(diff) / np.linalg.norm(target)

def powit_kPhi(M:np.ndarray, F:np.ndarray):
    '''For given M and F operators, power iterate to find k values and normalized flux vector'''
    k0 = 1
    phi0 = np.ones(N+1)
    phi0 /= la.norm(phi0)

    lk = k0
    lphi = phi0

    eps = 1e-6

    N_ITER = 1000
    # print(f"lk: {lk}")
    iM = la.inv(M)
    # print(f"iM: {iM.shape}, F: {F.shape}")
    # print(f"phi0: {lphi}, {lphi.shape}")
    # print(f"next phi: {(F @ lphi).shape}")
    for i in range(N_ITER):

        nphi = iM @ (F @ lphi) / lk
        # print(f"nphi: {nphi.shape}")
        nphi /= la.norm(nphi)
        nk = la.norm((F@nphi)[1:N], 1) / la.norm((M@nphi)[1:N], 1)


        conv_k = (nk-lk)**2 / (lk**2)
        conv_phi = err(nphi, lphi)

        if conv_k < eps and conv_phi < eps:
            print(f"k: {nk}")
            return nk, nphi
        else:
            lk = nk
            lphi = nphi
    print(f"Unable to converge in {N_ITER} iterations")
    return None, None
        

ab_vals = [(10, 50), (25, 50), (50, 50)]

NL = [10, 20, 40, 80, 160, 320, 640, 1280] 
NL = [10*2**i for i in range(5)]

# a, b = ab_vals[0]


phis_for_N = []
ks_for_N = []
for i in range(len(NL)):
    N = NL[i]
    phis = []
    ks = []
    for a, b in ab_vals:

        x = np.linspace(0, a+b, N+1)
        M, F = fin_diff_mats(N, a, b)
        k, phi = powit_kPhi(M, F)
        phis.append(phi)
        ks.append(k)
    phis_for_N.append(phis)
    ks_for_N.append(ks)

for i in range(len(ab_vals)):
    a, b = ab_vals[i]
    ks_for_ab = [ks_for_N[j][i] for j in range(len(NL))]
    strs = [f"& {k:.5f}" for k in ks_for_ab]
    print(f"{a} & {b} "+" ".join(strs)+"\\\\")

    # print(f"x: {x}")
    # print(f"phi: {phi}")
    # if phi is None:
    #     print("Did not solve ")
    #     continue
#     plt.plot(x, phi, label=f"Numerical, N = {N}, k = {k}")

# plt.xlabel("x (cm)")
# plt.ylabel("Normalized flux")
# plt.legend()
# plt.show()