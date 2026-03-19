import matplotlib.pyplot as plt
import numpy as np


h_vals = [1, 0.5, 0.25]
def phi(x):
    return np.sin(x)

def phi2p(x):
    return -np.sin(x)

def phi2p_approx(x, h):
    # Assumes x values are evenly spaced h apart
    phi2p_vals = (2*phi(x) - 5*phi(x+h) + 4*phi(x+2*h) - phi(x+3*h)) / (h**2)

    # Alternative formulation using one set of known points (i.e. not magically going past the boundary)
    # pv = phi(x)
    # phi2p_vals = np.zeros(len(pv))
    # phi2p_vals[:-3] = (2*pv[:-3] - 5*pv[1:-2] + 4*pv[2:-1] - pv[3:]) / (h**2)
    return phi2p_vals

xa = np.linspace(0, 10, 1000)
plt.plot(xa, phi(xa), label="phi(x) = sin(x)", linestyle='--')
for h in h_vals:
    x = np.arange(0, 10, h)
    y = phi2p_approx(x, h)
    plt.plot(x, y, label=f"phi''(x), h={h}")
plt.plot(xa, phi2p(xa), label="phi''(x), analytical")
plt.xlim(0, 2*np.pi)
ylim = plt.ylim()
plt.ylim(ylim[0], ylim[1]+1)
plt.xlabel("x")
plt.legend()
plt.show()