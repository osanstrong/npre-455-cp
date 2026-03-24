import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt 


# Finite difference solution 
D1 = 0.79 # cm
Xa1 = 0.066 # 1/cm
Xf1 = 0.02787
nu = 2.4
B_crit = ((nu*Xf1 - Xa1)/D1)**0.5

D2 = 1
Xa2 = 0.000709 # still 1/cm
L2 = (D2/Xa2)**0.5

# a = 10 # cm
# b = 20
def geom(N, a, b):
    h = (a+b)/N
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
    return h, D, Xa, Xf

def fin_diff_mats(N, a, b):
    # For given dimensions, with N cells
    h, D, Xa, Xf = geom(N, a, b)
    M_mat = np.eye(N+1)


    F_mat = np.zeros((N+1, N+1))
    for i in range(1, N):
        F_mat[i,i] = nu*Xf[i]

    M_mat[0,0:3] = [-3, 4, -1]
    for i in range(1,N):
        lD = D[i-1]
        iD = D[i]
        rD = D[i+1]
        lfrac = lD/(iD+lD)
        rfrac = rD/(iD+rD)

        vals = (D[i]/(h**2/2)) * np.array([
            -lfrac,
            lfrac + rfrac,
            -rfrac
        ]) + np.array([
            0, Xa[i], 0
        ])
        M_mat[i,i-1:i+2] = vals
    return M_mat, F_mat

def err(actual, target):
    diff = actual-target
    return np.linalg.norm(diff) / np.linalg.norm(target)

def Rf(phi: np.ndarray, Xf:np.ndarray, h):
    '''
    For given flux and fission cross section find fission reaction rate
    Uses trapezoidal integration across phi and
    '''
    Rf = sum([0.5*(phi[i-1]*Xf[i-1]+phi[i]*Xf[i])*h for i in range(1,len(Xf))])
    # print(f"Rf: {f}")
    return Rf
    

def powit_kPhi(M:np.ndarray, F:np.ndarray, h:float, Xf:np.ndarray|None=None):
    '''
    For given M and F operators, power iterate to find k values and normalized flux vector.
    
    '''
    k0 = 1
    phi0 = np.ones(N+1)
    phi0 /= la.norm(phi0)

    lk = k0
    lphi = phi0

    eps = 1e-6

    N_ITER = 1000
    iM = la.inv(M)
    for i in range(N_ITER):

        nphi = iM @ (F @ lphi) / lk
        # print(f"nphi: {nphi.shape}")
        rrate = Rf(nphi, Xf, h) if not Xf is None else la.norm(nphi, 1)*h
        nphi /= rrate
        nk = la.norm((F@nphi)[1:N], 1) / la.norm((M@nphi)[1:N], 1)


        conv_k = (nk-lk)**2 / (lk**2)
        conv_phi = err(nphi, lphi)

        if conv_k < eps and conv_phi < eps:
            print(f"k: {nk} after {i} iterations")
            return nk, nphi
        else:
            lk = nk
            lphi = nphi
    print(f"Unable to converge in {N_ITER} iterations")
    return None, None
        

ab_vals = [(10, 50), (25, 50), (50, 50)]
k_vals = [0.97356, 0.999608346014, 1.00796]

NL = [10, 20, 40, 80, 160, 320, 640, 1280] 
NL = [10*2**i for i in range(5)]

def xPhi_an(av, bv, k):
    '''
    Since k has to be analyzed with so many extra steps,
    just...provide it here
    '''
    N_POINTS = 1000
    bound = int(N_POINTS*av/(av+bv))+1
    x = np.linspace(0, av+bv, N_POINTS)
    phi = np.zeros(N_POINTS)
    x1 = x[:bound]
    x2 = x[bound:]


    B = ((nu*Xf1/k - Xa1)/D1)**0.5

    unscl_Rf = (Xf1*np.sin(B*av)/B)#*av
    print(f"Rf = {unscl_Rf}")
    C1 = 1  / unscl_Rf
    # C2 = -B*D1*L2*C1*np.sin(B*av)/(D2*np.cosh(bv/L2))
    
    C2 = C1* np.cos(B*av) / (np.sinh(-bv/L2))
    print(f"{C1*np.cos(B*av)} vs {C2*np.sinh(-bv/L2)}")

    phi[:bound] = C1*np.cos(B*x1)
    phi[bound:] = C2*np.sinh((x2-(av+bv))/L2)
    return x, phi
    
def vbar(x, label, c=None, style="dashed"):
    [ymin, ymax] = plt.ylim()
    plt.vlines(x, ymin, ymax, label=label, colors=c, linestyles=style)
    plt.ylim(ymin, ymax)
# a, b = ab_vals[0]


phis_for_N = []
ks_for_N = []
for i in range(len(NL)):
    N = NL[i]
    phis = []
    ks = []
    for a, b in ab_vals:

        x = np.linspace(0, a+b, N+1)
        h, D, Xa, Xf = geom(N, a, b)
        M, F = fin_diff_mats(N, a, b)
        k, phi = powit_kPhi(M, F, h, Xf)
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
for i in range(len(ab_vals)):
    a, b = ab_vals[i]
    

    for j in range(len(NL)):
        N = NL[j]
        k = ks_for_N[j][i]
        phi = phis_for_N[j][i]
        if phi is None:
            print(f"Did not solve ({a}, {b}), N={N}")
            continue
        x = np.linspace(0, a+b, N+1)
        plt.plot(x, phi, label=f"Numerical, N = {N}, k = {k:.5f}")
    x, phi = xPhi_an(a, b, k_vals[i])
    plt.plot(x, phi, label="Analytical")

    
    vbar(a, label=f"a ({a} cm)", c='C1')
    vbar(a+b, label=f"a+b (b = {b} cm)", c="C2")


    plt.xlabel("x (cm)")
    plt.ylabel("Normalized flux")
    plt.grid(which="both")
    plt.legend()
    plt.show()