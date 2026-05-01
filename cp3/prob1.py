import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rand



# P1a)

N_VALS = [5, 50, 500, 5000]
xi_sets = [rand.random(n) for n in N_VALS]

X_TOT = 0.1 # 1/cm

def inv_Fs(xi, Xt):
    x = -np.log(1-xi)/Xt
    return x


# x_sets = [inv_Fs(xi_set, X_TOT) for xi_set in xi_sets]

# for i, n in enumerate(N_VALS):
#     xi_set = xi_sets[i]
#     x_set = x_sets[i]
#     mean = sum(x_set) / n
#     plt.hist(x_set, density=True, label=f"N = {n}: d ~ {mean:.3f} cm")
#     plt.xlabel("x (cm)")
#     plt.ylabel("Relative Frequency")
#     plt.legend()
#     plt.show()


#### P1d
N_VALS = [10, 50, 500, 5000]

mX_i = [244, 337, 28, 119]
mX_t = sum(mX_i)
f_i = [mX/mX_t for mX in mX_i]
F_i = [sum(f_i[:i+1]) for i in range(len(f_i))]
print(F_i)
def get_i(xi):
    for i, F in enumerate(F_i):
        if xi < F:
            return i

xi_sets = [rand.random(n) for n in N_VALS]
i_sets = [np.array([get_i(xi) for xi in xi_set]) for xi_set in xi_sets]

for i, n in enumerate(N_VALS):
    xi_set = xi_sets[i]
    i_set = i_sets[i]
    plt.hist(i_set, weights=np.zeros_like(i_set) + 1. / n, label=f"N = {n}")
    plt.xlabel("Reaction Type")
    plt.ylabel("Relative Frequency")
    plt.legend()
    plt.show()


