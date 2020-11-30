import numpy as np
import matplotlib.pyplot as plt


data           = np.loadtxt("time_series_K_inj=0,3_DnuInj=2,8.txt")
time           = data[:, 0]
intensity_1    = data[:, 1]
rho_GS_e_act_1 = data[:, 6]
rho_GS_h_act_1 = data[:, 7]
inversion_1    = rho_GS_e_act_1 + rho_GS_h_act_1 - 1.0

data           = np.loadtxt("time_series_K_inj=0,1_DnuInj=-2.txt")
intensity_2    = data[:, 1]
rho_GS_e_act_2 = data[:, 6]
rho_GS_h_act_2 = data[:, 7]
inversion_2    = rho_GS_e_act_2 + rho_GS_h_act_2 - 1.0

data           = np.loadtxt("time_series_K_inj=0,3_DnuInj=-4,25.txt")
intensity_3    = data[:, 1]
rho_GS_e_act_3 = data[:, 6]
rho_GS_h_act_3 = data[:, 7]
inversion_3    = rho_GS_e_act_3 + rho_GS_h_act_3 - 1.0

data           = np.loadtxt("time_series_K_inj=0,37_DnuInj=-4,1.txt")
intensity_4    = data[:, 1]
rho_GS_e_act_4 = data[:, 6]
rho_GS_h_act_4 = data[:, 7]
inversion_4    = rho_GS_e_act_4 + rho_GS_h_act_4 - 1.0

data           = np.loadtxt("time_series_K_inj=0,3_DnuInj=-3,8.txt")
intensity_5    = data[:, 1]
rho_GS_e_act_5 = data[:, 6]
rho_GS_h_act_5 = data[:, 7]
inversion_5    = rho_GS_e_act_5 + rho_GS_h_act_5 - 1.0


fig, axes = plt.subplots(2, 5)
fig.set_size_inches(5.9, 2.8)
fig.subplots_adjust(hspace=0.3, top=0.99, bottom=0.1, left=0.05, right=0.99)
plt.rcParams.update({"font.size": 9})

# time series
axes[0, 0].plot(time, intensity_1, color="crimson", linewidth=0.75)
axes[0, 1].plot(time, intensity_2, color="crimson", linewidth=0.75)
axes[0, 2].plot(time, intensity_3, color="crimson", linewidth=0.75)
axes[0, 3].plot(time, intensity_4, color="crimson", linewidth=0.75)
axes[0, 4].plot(time, intensity_5, color="crimson", linewidth=0.75)

# phase space
axes[1, 0].plot(inversion_1, intensity_1, color="orange", linewidth=0.75)
axes[1, 1].plot(inversion_2, intensity_2, color="orange", linewidth=0.75)
axes[1, 2].plot(inversion_3, intensity_3, color="orange", linewidth=0.75)
axes[1, 3].plot(inversion_4, intensity_4, color="orange", linewidth=0.75)
axes[1, 4].plot(inversion_5, intensity_5, color="orange", linewidth=0.75)

# xlims
axes[0, 0].set_xlim(50.20, 51.10)
axes[0, 1].set_xlim(50.51, 52.28)
axes[0, 2].set_xlim(50.18, 51.41)
axes[0, 3].set_xlim(50.54, 52.74)
axes[0, 4].set_xlim(50.20, 52.19)

# numbers
axes[1, 0].text(0.825, 0.825, "1)", transform = axes[1, 0].transAxes)
axes[1, 1].text(0.825, 0.825, "2)", transform = axes[1, 1].transAxes)
axes[1, 2].text(0.825, 0.825, "3)", transform = axes[1, 2].transAxes)
axes[1, 3].text(0.825, 0.825, "4)", transform = axes[1, 3].transAxes)
axes[1, 4].text(0.825, 0.825, "5)", transform = axes[1, 4].transAxes)


for i in range(5):
	plt.setp(axes[0, i].get_xticklabels(), visible=False)
	plt.setp(axes[1, i].get_xticklabels(), visible=False)
	plt.setp(axes[0, i].get_yticklabels(), visible=False)
	plt.setp(axes[1, i].get_yticklabels(), visible=False)
	axes[0, i].tick_params(axis="both", which="both", length=0.0)
	axes[1, i].tick_params(axis="both", which="both", length=0.0)


axes[0, 2].set_xlabel(r"time $t$", fontsize=9.0, labelpad=4.5)
axes[1, 2].set_xlabel(r"inversion $\rho^{act}_{GS,e} + \rho^{act}_{GS,h} - 1$", fontsize=9.0, labelpad=5.0)

axes[0, 0].set_ylabel(r"intensity $|E|^2$", fontsize=9)
axes[1, 0].set_ylabel(r"intensity $|E|^2$", fontsize=9)

plt.show()
