import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt("time_series_J=10.txt")
time_10      = data[:, 0]
intensity_10 = data[:, 1]

data = np.loadtxt("time_series_J=11.txt")
time_11      = data[:, 0]
intensity_11 = data[:, 1]

data = np.loadtxt("time_series_J=12.txt")
time_12      = data[:, 0]
intensity_12 = data[:, 1]

data = np.loadtxt("time_series_J=13.txt")
time_13      = data[:, 0]
intensity_13 = data[:, 1]

data = np.loadtxt("time_series_J=14.txt")
time_14      = data[:, 0]
intensity_14 = data[:, 1]

data = np.loadtxt("time_series_J=15.txt")
time_15      = data[:, 0]
intensity_15 = data[:, 1]

data = np.loadtxt("time_series_J=16.txt")
time_16      = data[:, 0]
intensity_16 = data[:, 1]

data = np.loadtxt("time_series_J=17.txt")
time_17      = data[:, 0]
intensity_17 = data[:, 1]

data = np.loadtxt("time_series_J=18.txt")
time_18      = data[:, 0]
intensity_18 = data[:, 1]

data = np.loadtxt("time_series_J=19.txt")
time_19      = data[:, 0]
intensity_19 = data[:, 1]


plt.figure(figsize=(4.6, 3.6))
# ~ plt.suptitle("MDRE: single laser", fontsize=14)
# ~ plt.title(r"d$t=10^{-3}\,$ns, $\tau = 4\,$ns, $K = 0$, $\eta = 1$, $\varphi = 0$, method: Runge-Kutta 4")
plt.subplots_adjust(top=0.98, bottom=0.13, left=0.14, right=0.97)
plt.plot(time_19, intensity_19, label=r"$J=19\,$ns$^{-1}$")
plt.plot(time_18, intensity_18, label=r"$J=18\,$ns$^{-1}$")
plt.plot(time_17, intensity_17, label=r"$J=17\,$ns$^{-1}$")
plt.plot(time_16, intensity_16, label=r"$J=16\,$ns$^{-1}$")
plt.plot(time_15, intensity_15, label=r"$J=15\,$ns$^{-1}$")
plt.plot(time_14, intensity_14, label=r"$J=14\,$ns$^{-1}$")
plt.plot(time_13, intensity_13, label=r"$J=13\,$ns$^{-1}$")
plt.plot(time_12, intensity_12, label=r"$J=12\,$ns$^{-1}$")
plt.plot(time_11, intensity_11, label=r"$J=11\,$ns$^{-1}$")
plt.plot(time_10, intensity_10, label=r"$J=10\,$ns$^{-1}$")
plt.xlabel(r"time $t$ / ns", size=9.0)
plt.ylabel(r"intensity $|E|^2$")
plt.xlim(-0.25, 20.25)
plt.ylim(-0.05, np.max(intensity_19) + 0.05)
plt.legend(bbox_to_anchor=(1, 1), loc="upper right", ncol=2)
plt.grid(color="lightgray")
plt.show() 

