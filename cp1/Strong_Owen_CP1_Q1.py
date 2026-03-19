import matplotlib.pyplot as plt
import numpy as np


h_vals = [1, 0.5, 0.25]
def phi(x):
    return np.sin(x)

def phi2p(x):
    return -np.sin(x)

def phi2p_approx(x, h):
    # Assumes x values are evenly spaced h apart
    pv = phi(x)
    l = len(pv)-3 # Use 4 terms to need to stop 3 before
    phi2p_vals = np.zeros(l+3)

    # phi2p_vals[:-3] = (pv[:-3] - 5*pv[1:-2] + 4*pv[2:-1] - pv[3:]) / (h**2)
    for i in range(l):
        phi2p_vals[i] = (2*pv[i] - 5*pv[i+1] + 4*pv[i+2] - pv[i+3]) / (h**2)
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