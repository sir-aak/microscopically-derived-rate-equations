import numpy as np
import matplotlib.pyplot as plt


data          = np.loadtxt("time_series_K_fb=0,35_phi=1,5.txt")
time          = data[:, 0]
intensity_1_1 = data[:, 1]
intensity_1_2 = data[:, 2]
intensity_1_3 = data[:, 3]
intensity_1_4 = data[:, 4]

data          = np.loadtxt("time_series_K_fb=0,328_phi=2.86.txt")
intensity_2_1 = data[:, 1]
intensity_2_2 = data[:, 2]
intensity_2_3 = data[:, 3]
intensity_2_4 = data[:, 4]

data          = np.loadtxt("time_series_K_fb=0,47_phi=4.txt")
intensity_3_1 = data[:, 1]
intensity_3_2 = data[:, 2]
intensity_3_3 = data[:, 3]
intensity_3_4 = data[:, 4]

data          = np.loadtxt("time_series_K_fb=0,36_phi=3,35.txt")
intensity_4_1 = data[:, 1]
intensity_4_2 = data[:, 2]
intensity_4_3 = data[:, 3]
intensity_4_4 = data[:, 4]

data          = np.loadtxt("time_series_K_fb=0,105_phi=3,55.txt")
intensity_5_1 = data[:, 1]
intensity_5_2 = data[:, 2]
intensity_5_3 = data[:, 3]
intensity_5_4 = data[:, 4]

data          = np.loadtxt("time_series_K_fb=0,25_phi=2,65.txt")
intensity_6_1 = data[:, 1]
intensity_6_2 = data[:, 2]
intensity_6_3 = data[:, 3]
intensity_6_4 = data[:, 4]


fig, axes = plt.subplots(2, 3)
fig.set_size_inches(5.9, 4.5)
fig.subplots_adjust(hspace=0.2, top=0.95, bottom=0.06, left=0.05, right=0.99)
plt.rcParams.update({"font.size": 9})

linewidth = 1.5

# time series
axes[0, 0].plot(time, intensity_1_1, linewidth=linewidth, label="laser 1")
axes[0, 0].plot(time, intensity_1_2, linestyle="dashdot", linewidth=linewidth, label="laser 2")
axes[0, 0].plot(time, intensity_1_3, linestyle="dashed", linewidth=linewidth, label="laser 3")
axes[0, 0].plot(time, intensity_1_4, linestyle="dotted", linewidth=linewidth, label="laser 4")

axes[0, 1].plot(time, intensity_2_1, linewidth=linewidth)
axes[0, 1].plot(time, intensity_2_2, linestyle="dashed", linewidth=linewidth)
axes[0, 1].plot(time, intensity_2_3, linewidth=linewidth)
axes[0, 1].plot(time, intensity_2_4, linestyle="dotted", linewidth=linewidth)

axes[0, 2].plot(time, intensity_3_1, linewidth=linewidth)
axes[0, 2].plot(time, intensity_3_2, linewidth=linewidth)
axes[0, 2].plot(time, intensity_3_3, linestyle="dashed", linewidth=linewidth)
axes[0, 2].plot(time, intensity_3_4, linestyle="dashed", linewidth=linewidth)

axes[1, 0].plot(time, intensity_4_1, linewidth=linewidth)
axes[1, 0].plot(time, intensity_4_2, linewidth=linewidth)
axes[1, 0].plot(time, intensity_4_3, linestyle="dashed", linewidth=linewidth)
axes[1, 0].plot(time, intensity_4_4, linewidth=linewidth)

axes[1, 1].plot(time, intensity_5_1, linewidth=linewidth)
axes[1, 1].plot(time, intensity_5_2, linewidth=linewidth)
axes[1, 1].plot(time, intensity_5_3, linewidth=linewidth)
axes[1, 1].plot(time, intensity_5_4, linewidth=linewidth)

axes[1, 2].plot(time, intensity_6_1, linewidth=linewidth)
axes[1, 2].plot(time, intensity_6_2, linewidth=linewidth)
axes[1, 2].plot(time, intensity_6_3, linewidth=linewidth)
axes[1, 2].plot(time, intensity_6_4, linewidth=linewidth)

legend = axes[0, 0].legend(loc="upper right")
for lh in legend.legendHandles: 
	lh.set_linestyle("solid")

# titles
axes[0, 0].set_title("synchronized", fontsize=9.0)
axes[0, 1].set_title("cluster 1-3", fontsize=9.0)
axes[0, 2].set_title("cluster 2-2", fontsize=9.0)
axes[1, 0].set_title("chimera", fontsize=9.0)
axes[1, 1].set_title("splay", fontsize=9.0)
axes[1, 2].set_title("desynchronized", fontsize=9.0)

# xlims
axes[0, 0].set_xlim(998.0, 1000.0)
axes[0, 0].set_ylim(2.2, 2.65)
axes[0, 1].set_xlim(999.25, 999.9)
axes[0, 2].set_xlim(998.05, 1000.0)
axes[1, 0].set_xlim(998.24, 1000.0)
axes[1, 1].set_xlim(998.85, 999.99)
axes[1, 2].set_xlim(998.0, 1000.0)

# numbers
axes[0, 0].text(0.04, 0.9, "1)", transform = axes[0, 0].transAxes)
axes[0, 1].text(0.04, 0.9, "2)", transform = axes[0, 1].transAxes)
axes[0, 2].text(0.04, 0.9, "3)", transform = axes[0, 2].transAxes)
axes[1, 0].text(0.04, 0.9, "4)", transform = axes[1, 0].transAxes)
axes[1, 1].text(0.04, 0.9, "5)", transform = axes[1, 1].transAxes)
axes[1, 2].text(0.04, 0.9, "6)", transform = axes[1, 2].transAxes)


for i in range(3):
	plt.setp(axes[0, i].get_xticklabels(), visible=False)
	plt.setp(axes[1, i].get_xticklabels(), visible=False)
	plt.setp(axes[0, i].get_yticklabels(), visible=False)
	plt.setp(axes[1, i].get_yticklabels(), visible=False)
	axes[0, i].tick_params(axis="both", which="both", length=0.0)
	axes[1, i].tick_params(axis="both", which="both", length=0.0)


axes[1, 1].set_xlabel(r"time $t$", fontsize=9.0, labelpad=8.5)
# ~ axes[1, 2].set_xlabel(r"inversion $\rho^{act}_{GS,e} + \rho^{act}_{GS,h} - 1$", fontsize=9.0, labelpad=5.0)

axes[0, 0].set_ylabel(r"intensity of $i$-th laser $|E|^2$", fontsize=9)
axes[1, 0].set_ylabel(r"intensity of $i$-th laser $|E|^2$", fontsize=9)

plt.show()
