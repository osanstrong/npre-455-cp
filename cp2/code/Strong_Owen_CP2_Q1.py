import matplotlib.pyplot as plt
import numpy as np
from analytic_solution_cp2 import flux1, flux2
import analytic_solution_cp2 as an


a = an.a
b = an.b

D_C1 = an.D1C # Region C, group 1 (fast), cm
D_C2 = an.D2C
D_M1 = an.DM
D_M2 = an.D2M

L_C1 = an.L1C # cm
L_C2 = an.L2C
L_M1 = an.LM
L_M2 = an.L2M

Xa_C1 = D_C1 / (L_C1**2)
Xa_C2 = D_C2 / (L_C2**2)
Xa_M1 = D_M1 / (L_M1**2)
Xa_M2 = D_M2 / (L_M2**2)

Xs_C12 = an.Es12C # Region C, scattering cross section from 1 to 2, cm⁻¹
Xs_M12 = an.Es12M

S_C1 = 1e12 #Region C, group 1 source, cm⁻³s⁻¹



def fin_diff_sol(N, debug=False):
    '''
    Returns a finite difference solution vector to the flux distribution.
    In this vector, the first half (N+1) represents the fast flux across these nodes.
    The second half (N+1, offset N+1) represents the thermal flux across these nodes.
    '''


    h = (a+b)/N
    nodes = N+1 # Nodes of solution, one on the side of each cell. Must be repeated for each energy group
    

    
    A_mat = np.eye(nodes*2)

    bound = int(N*a/(a+b))+1
    

    D = np.zeros(nodes*2)
    D[:bound] = D_C1
    D[bound:nodes] = D_M1
    D[nodes:nodes+bound] = D_C2
    D[nodes+bound:] = D_M2

    L = np.zeros(nodes*2)
    L[:bound] = L_C1
    L[bound:nodes] = L_M1
    L[nodes:nodes+bound] = L_C2
    L[nodes+bound] = L_M2

    S_vals = np.zeros(nodes*2)
    S_vals[:bound] = S_C1

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

    b_mat = np.zeros(nodes*2)
    b_mat[:] = S_vals[:]
    # Now, GENERALLY, we should set the endpoints to known values (either Dirichlet or Neumann).
    # However, here we have no source in symmetry (so the derivative at left is just 0 in both groups)
    # And there is a vacuum at right (flux on right approx. 0, OR set left-facing current to 0, value itself being 0 either way)
    b_mat[0] = 0
    b_mat[nodes-1] = 0
    b_mat[nodes] = 0
    b_mat[2*nodes-1] = 0
    
    # Boundary conditions 

    # Derivative on left boundary (for net-zero current)
    # right_dx = [-1, 1] # 2 nodes
    right_dx = [-3, 4, -1] # 3 nodes
    # right_dx = [-11, 18, -9, 2] # 4 nodes
    # Fast flux
    A_mat[0,0:len(right_dx)] = right_dx
    # Thermal flux
    A_mat[nodes,nodes:nodes+len(right_dx)] = right_dx
    
    # Value on right boundary (for vacuum, optionally replace with leftward current operator)
    A_mat[nodes-1, nodes-1] = 1
    A_mat[nodes*2-1, nodes*2-1] = 1
    

    # Fast Flux
    for i in range(1,nodes-1):
        lD = D[i-1]
        iD = D[i]
        rD = D[i+1]
        hD = iD / (0.5 * h**2)
        lfrac = lD/(iD+lD)
        rfrac = rD/(iD+rD)

        # No upscattering, so just diffusion and absorption+downscattering (removal)
        vals = hD * np.array([
            -lfrac,
            lfrac + rfrac,
            -rfrac
        ]) + np.array([
            # 0, Xa[i]+Xs_12[i], 0
            0, Xa[i], 0 
        ])
        A_mat[i,i-1:i+2] = vals

    if (debug):
        print("_"*10)
        print(f"N={N}, bound={bound}")
        print(f"A: {A_mat}") 
    
    # Thermal Flux
    for i in range(nodes+1, nodes*2-1):
        lD = D[i-1]
        iD = D[i]
        rD = D[i+1]
        hD = iD / (0.5 * h**2)
        lfrac = lD/(iD+lD)
        rfrac = rD/(iD+rD)

        # Downscattering, so both diffusion and absorption/removal out + downscattering in
        vals = hD * np.array([
            -lfrac,
            lfrac + rfrac,
            -rfrac
        ]) + np.array([
            0, Xa[i], 0
        ])
        A_mat[i,i-1:i+2] = vals
        A_mat[i,i-nodes] = -Xs_12[i-nodes]
        

    if (debug):
        print("_"*10)
        print(f"N={N}, bound={bound}")
        print(f"A: {A_mat}")
        print(f"b: {b_mat}")
        print(f"S: {S_vals}")
        print(f"Xa: {Xa}")
        print(f"Xs12: {Xs_12}")
    phi = np.linalg.solve(A_mat, b_mat)
    return phi






def plot_anal(n):
    xpts, phi1 = flux1(n)
    xpts, phi2 = flux2(n)

    plt.plot(xpts, phi1, label="Fast Flux (Analytical)")
    plt.plot(xpts, phi2, label="Thermal Flux (Analytical)")


# N = 10
for N in [1000]:
    np.set_printoptions(linewidth=1000, precision=3)
    phi = fin_diff_sol(N, debug=False)
    x = np.array([n for n in range(N+1)]) * (a+b)/N
    plt.plot(x, phi[:N+1], label=f"Fast Flux (N={N})")
    plt.plot(x, phi[N+1:], label=f"Thermal Flux (N={N})")
# plt.plot(x, phi, label="")
plot_anal(1000)
plt.xlabel("x (cm)")
plt.ylabel("Flux (1/cm²s)")
plt.legend()
plt.show()


