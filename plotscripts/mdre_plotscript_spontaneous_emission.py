import numpy as np
import matplotlib.pyplot as plt


data        = np.loadtxt("time_series_beta=0.txt")
time_0      = data[:, 0]
intensity_0 = data[:, 1]

data        = np.loadtxt("time_series_beta=1e-9.txt")
time_9      = data[:, 0]
intensity_9 = data[:, 1]

data        = np.loadtxt("time_series_beta=1e-8.txt")
time_8      = data[:, 0]
intensity_8 = data[:, 1]

data        = np.loadtxt("time_series_beta=1e-7.txt")
time_7      = data[:, 0]
intensity_7 = data[:, 1]

data        = np.loadtxt("time_series_beta=1e-6.txt")
time_6      = data[:, 0]
intensity_6 = data[:, 1]

data        = np.loadtxt("time_series_beta=1e-5.txt")
time_5      = data[:, 0]
intensity_5 = data[:, 1]

data        = np.loadtxt("time_series_beta=1e-4.txt")
time_4      = data[:, 0]
intensity_4 = data[:, 1]

data        = np.loadtxt("time_series_beta=1e-3.txt")
time_3      = data[:, 0]
intensity_3 = data[:, 1]

data        = np.loadtxt("time_series_beta=1e-2.txt")
time_2      = data[:, 0]
intensity_2 = data[:, 1]

data        = np.loadtxt("time_series_beta=1e-1.txt")
time_1      = data[:, 0]
intensity_1 = data[:, 1]


# pdf looks different than the output - very strange
plt.figure(figsize=(4.9, 2.72))
plt.rcParams.update({"font.size": 9.0})
# ~ plt.suptitle("MDRE: single laser", fontsize=14)
# ~ plt.title(r"d$t=10^{-3}\,$ns, $\tau = 4\,$ns, $K = 0$, $\eta = 1$, $\varphi = 0$, method: Runge-Kutta 4")
plt.subplots_adjust(top=0.98, bottom=0.14, left=0.10, right=0.99)
plt.plot(time_0, intensity_0, label=r"$\beta=0$")
plt.plot(time_9, intensity_9, label=r"$\beta=10^{-9}$")
plt.plot(time_8, intensity_8, label=r"$\beta=10^{-8}$")
plt.plot(time_7, intensity_7, label=r"$\beta=10^{-7}$")
plt.plot(time_6, intensity_6, label=r"$\beta=10^{-6}$")
plt.plot(time_5, intensity_5, label=r"$\beta=10^{-5}$")
plt.plot(time_4, intensity_4, label=r"$\beta=10^{-4}$")
plt.plot(time_3, intensity_3, label=r"$\beta=10^{-3}$")
plt.plot(time_2, intensity_2, label=r"$\beta=10^{-2}$")
plt.plot(time_1, intensity_1, label=r"$\beta=10^{-1}$")
plt.xlabel(r"time $t$ / ns", size=9.0)
plt.ylabel(r"intensity $|E|^2$", size=9.0)
plt.xlim(-0.05, 5.55)
plt.ylim(-0.05, 3.05)
plt.yticks([0.0, 1.0, 2.0, 3.0])
plt.legend(loc="upper right", ncol=1)
plt.grid(color="lightgray")
plt.show() 

