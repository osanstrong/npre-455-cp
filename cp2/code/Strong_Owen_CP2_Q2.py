import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
import csv
from math import factorial as fac

def get_flux(path):
    '''
    Return x and corresponding flux of the given reference flux file
    '''
    with open(path, "r") as file:
        reader = csv.reader(file)
        x = []
        f = []
        for row in reader:
            x.append(float(row[0]))
            f.append(float(row[1]))
        return np.array(x), np.array(f)
    
def plot_1050():
    xf, ff = get_flux("fast_flux_a10.csv")
    xt, ft = get_flux("thermal_flux_a10.csv")
    plt.plot(xf, ff, label="Fast (reference)")
    plt.plot(xt, ft, label="Thermal (reference)")

def plot_2550():
    xf, ff = get_flux("fast_flux_a25.csv")
    xt, ft = get_flux("thermal_flux_a25.csv")
    plt.plot(xf, ff, label="Fast (reference)")
    plt.plot(xt, ft, label="Thermal (reference)")



# a = 10 # cm
a = 25 #cm
b = 50 # cm

def plot_reference(ax1, ax2, c=None):
    if a == 10:
        paths = ["fast_flux_a10.csv", "thermal_flux_a10.csv"]
    else:
        paths = ["fast_flux_a25.csv", "thermal_flux_a25.csv"]
    
    xf, ff = get_flux(paths[0])
    xt, ft = get_flux(paths[1])
    ax1.plot(xf, ff, label="Fast (reference)", c=c)
    ax2.plot(xt, ft, label="Thermal (reference)", c=c, linestyle="dashed")

def set_params(new_a):
    global a
    a = new_a
    global D_C1, D_C2, D_M1, D_M2, Xs_C12, Xs_M12, nXf_C1, nXf_C2, Xa_C1, Xa_C2, Xa_M1, Xa_M2
    if a == 10:
        D_C1 = 0.981
        D_C2 = 0.284
        D_M1 = 1.06
        D_M2 = 0.179

        Xs_C12 = 5.36e-5
        Xs_M12 = 3.95e-2

        nXf_C1 = 1.58e-2 # nu times Sigma f, cm⁻¹
        nXf_C2 = 1.56

        Xa_C1 = 1.22e-2
        Xa_C2 = 7.82e-1
        Xa_M1 = 3.39e-4
        Xa_M2 = 1.41e-2
    else:
        D_C1 = 9.27e-1
        D_C2 = 2.85e-1
        D_M1 = 1.04
        D_M2 = 1.79e-1

        Xs_C12 = 1.53e-5
        Xs_M12 = 4.07e-2

        nXf_C1 = 1.29e-2
        nXf_C2 = 1.56

        Xa_C1 = 1.11e-2
        Xa_C2 = 7.8e-1
        Xa_M1 = 3.44e-4
        Xa_M2 = 1.41e-2



def geom(N, a, b):
    h = (a+b)/N

    bound = int(N*a/(a+b))+1
    nodes = N+1
    

    D = np.zeros(nodes*2)
    D[:bound] = D_C1
    D[bound:nodes] = D_M1
    D[nodes:nodes+bound] = D_C2
    D[nodes+bound:] = D_M2

    # Absorption cross section
    Xa = np.zeros(nodes*2)
    Xa[:bound] = Xa_C1
    Xa[bound:nodes] = Xa_M1
    Xa[nodes:nodes+bound] = Xa_C2
    Xa[nodes+bound:] = Xa_M2

    # Just scattering cross section
    # Since these are defined for combinations of energy groups,
    # you would *really* need a list of matrices corresponding to each node.
    # Since scattering is simplified to only one combo (from 1 to 2),
    # we will only define one array, only the length of spatial nodes.
    Xs_12 = np.zeros(nodes)
    Xs_12[:bound] = Xs_C12
    Xs_12[bound:] = Xs_M12

    # Just fission cross section (and multiplication); we assume only fission in left region
    nXf = np.zeros(nodes*2)
    nXf[:bound] = nXf_C1
    nXf[nodes:nodes+bound] = nXf_C2

    return h, D, Xa, nXf, Xs_12

def fin_diff_mats(N, a, b):
    # For given dimensions, with N cells
    h, D, Xa, nXf, Xs_12 = geom(N, a, b)
    nodes = N+1
    M_mat = np.eye(nodes*2)


    F_mat = np.zeros((nodes*2, nodes*2))
    # Assume only fission produced in fast neutrons
    for i in range(1, nodes-1):
        F_mat[i,i] = nXf[i]
        F_mat[i,i+nodes] = nXf[i+nodes]

    right_dx = [-3, 4, -1]
    M_mat[0,0:len(right_dx)] = right_dx
    for i in range(1,nodes-1):
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
            0, Xa[i]+Xs_12[i], 0
        ])
        M_mat[i,i-1:i+2] = vals

    M_mat[nodes,nodes:nodes+len(right_dx)] = right_dx
    for i in range(nodes+1,2*nodes-1):
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

        M_mat[i,i-nodes] = -Xs_12[i-nodes]
    return M_mat, F_mat

def err(actual, target):
    diff = actual-target
    return np.linalg.norm(diff) / np.linalg.norm(target)

def Rf(phi: np.ndarray, nXf: np.ndarray, h):
    '''
    For given flux and fission cross section find neutron production rate
    Uses trapezoidal integration across x of phi and cross section
    '''
    Rf = sum([0.5*(phi[i-1]*nXf[i-1] + phi[i]*nXf[i])*h for i in range(0,len(nXf))])
    return Rf
def trapz(dist: np.ndarray, h):
    return sum([0.5*(dist[i-1]+dist[i])*h for i in range(len(dist))])
    

def powit_kPhi(M:np.ndarray, F:np.ndarray, h:float, nXf:np.ndarray, debug=False):
    '''
    For given M and F operators, power iterate to find k values and normalized flux vector.
    
    '''
    len_phi = M.shape[0]
    k0 = 1
    phi0 = np.ones(len_phi)
    phi0 /= la.norm(phi0, 1)

    lk = k0
    lphi = phi0

    eps = 1e-6

    N_ITER = 1000
    iM = la.inv(M)
    print(nXf)
    for i in range(N_ITER):

        nphi = iM @ (F @ lphi) / lk
        # print(f"nphi: {nphi.shape}")
        # rrate = trapz(F@nphi, h)
        rrate = Rf(nphi, nXf, h)
        # rrate = Rf(nphi, Xf, nu, h) if not Xf is None else la.norm(nphi, 1)*h
        nphi /= rrate
        # nk = la.norm((F@nphi)[1:N], 1) / la.norm((M@nphi)[1:N], 1)
        nk = trapz(F@nphi, h) / trapz(M@nphi, h)


        conv_k = (nk-lk)**2 / (lk**2)
        conv_phi = err(nphi, lphi)

        if conv_k < eps and conv_phi < eps:
            if debug: print(f"k: {nk} after {i} iterations")
            rrate = Rf(nphi, nXf, h)
            nphi /= rrate
            return nk, nphi
        else:
            lk = nk
            lphi = nphi
    if debug: print(f"Unable to converge in {N_ITER} iterations")
    return None, None

# N = 100

def vbar(x, label, c=None, style="dashed"):
    [ymin, ymax] = plt.ylim()
    plt.vlines(x, ymin, ymax, label=label, colors=c, linestyles=style)
    plt.ylim(ymin, ymax)

def calc_kPhi(N):
    h = (a+b)/N
    M, F = fin_diff_mats(N, a, b)
    h, D, Xa, nXf, Xs_12 = geom(N, a, b)
    k, phi = powit_kPhi(M, F, h, nXf)
    return k, phi


def plot_findiff(N, ax1, ax2, c=None):
    h = (a+b)/N
    k, phi = calc_kPhi(N)
    if phi is None: return None, None
    nodes = N+1
    x = h*np.array([n for n in range(nodes)])
    ax1.plot(x, phi[:nodes], label=f"Fast (N={N})", c=c)
    ax2.plot(x, phi[nodes:], label=f"Thermal (N={N})", c=c, linestyle="dashed")
    return k, phi

# 2e: plots and k values
k_vals = {
    10:[],
    25:[]
}
A_VALS = [10, 25]
N_VALS = [20, 40, 80, 500]
for a_val in A_VALS:
    set_params(a_val)
    i = 1
    fig, ax1 = plt.subplots()
    ax1.set_xlabel("x (cm)")
    ax1.set_ylabel("Fast Flux (cm⁻²s⁻¹)")
    ax2 = plt.twinx(ax1)
    ax2.set_ylabel("Thermal Flux (cm⁻²s⁻¹)")
    for N in N_VALS:
        k, phi = plot_findiff(N, ax1, ax2,c=f"C{i}")
        k_vals[a_val].append(k)
        i+=1
    plot_reference(ax1, ax2, f"C{i}")
    vbar(a, f"a = {a} cm", c="C0", style="dotted")
    plt.grid(which="both")
    ax1.grid(True)
    ax2.grid(True)
    ax1.set_axisbelow(True)
    ax2.set_axisbelow(True)
    ax1.legend(loc=5)
    ax2.legend(loc=0)
    plt.show()
print(k_vals)
# 2e: k values
print(f"$N$ & " + " & ".join([str(N) for N in N_VALS]) + "\\\\")
print(f"$a={10}$ & " + " & ".join([f"{k:.5f}" for k in k_vals[10]]) + "\\\\")
print(f"$a={25}$ & " + " & ".join([f"{k:.5f}" for k in k_vals[25]]) + "\\\\")

# 2f: 
set_params(new_a=25)
MAX_N = 500
k, phi = calc_kPhi(N=MAX_N)
h, D, Xa, nXf, Xs_12 = geom(MAX_N, a, b)
nu = 2.4
rates = nXf*phi / nu
# rrate = Rf(phi, nXf , h)
u = h*np.array([n for n in range(MAX_N+1)])
plt.plot(u, rates[:MAX_N+1] + rates[MAX_N+1:], label=f"Reaction Rate (N={MAX_N})")
vbar(a, f"a = {a} cm", c="C0", style="dotted")
plt.xlabel("x (cm)")
plt.ylabel("Flux (cm⁻²s⁻¹)")
plt.grid(which="both")
plt.legend()
plt.show()


# 2h: Legendre expansion
def Pnx(n, x:np.ndarray):
    '''
    Evaluate the n'th Legendre Polynomial over the given values of x.
    Uses Rodrigues' formula
    '''
    c = np.zeros(n+1)
    c[n] = 1
    return np.array(np.polynomial.legendre.legval(x, c))

for a_val in [10, 25]:
    set_params(a_val)
    k, phi_fine = calc_kPhi(MAX_N)
    if not phi_fine is None:
        nodes_big = MAX_N+1
        u = np.linspace(-1, 1, nodes_big*2)
        
        # NOTE: The manual asked for 50 nodes, but it is possible one might desire 50 nodes *per group*, in which case 100 nodes would be appropriate.
        num_terms = 50
        # num_terms = 100
        A = np.zeros((nodes_big*2, num_terms))
        for n in range(num_terms):
            P = Pnx(n, u)
            A[:,n] = P
        # A.T @ A @ x = A.T @ b
        L = la.inv(A.T @ A) @ A.T @ phi_fine

        short_terms = 10
        print(f" $i$ & " + " & ".join([str(i) for i in range(short_terms)]) + "\\\\")
        print(f" $C_i$ & " + " & ".join([f"{Ci:.3f}" for Ci in L[:short_terms]]) + "\\\\")


        L_eval = np.array(np.polynomial.legendre.Legendre(L)(u))

        x = np.linspace(0, a+b, nodes_big)
        plt.plot(x, L_eval[:nodes_big], label="Fast")
        plt.plot(x, L_eval[nodes_big:], label="Thermal")
        vbar(a, f"a = {a} cm", c="C2")
        plt.xlabel("x (cm)")
        plt.ylabel("Flux (cm⁻²s⁻¹)")
        plt.grid(which="both")
        plt.legend()
        plt.show()
