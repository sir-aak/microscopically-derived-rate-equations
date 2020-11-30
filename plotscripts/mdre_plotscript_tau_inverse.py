import numpy as np
import matplotlib.pyplot as plt


scalingFactor = np.array([1e-2, 2e-2, 5e-2, 1e-1, 2e-1, 5e-1, 1.0])
tau_e_inverse = np.array([0.149, 0.0775, 0.0332, 0.0173, 0.00904, 0.00388, 0.00203])
tau_h_inverse = np.array([0.087, 0.0442, 0.0182, 0.00930, 0.00473, 0.00195, 0.000991])


plt.figure(figsize=(5.9, 3.0))
plt.rcParams.update({'font.size': 18})
plt.subplots_adjust(top=0.99, bottom=0.24, left=0.24, right=0.99)
plt.loglog(scalingFactor, tau_e_inverse, color="blue", label=r"$\tau_e^{-1}$")
plt.loglog(scalingFactor, tau_h_inverse, color="red", label=r"$\tau_h^{-1}$")
plt.grid(color="lightgray")
plt.xlabel("scaling factor $s$")
plt.ylabel("inverse scattering" + "\n" + r"lifetime $\tau_b^{-1}$ / ns$^{-1}$")
plt.legend(loc="upper right")

plt.show()
