import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


# plots intensity time series for MDRE model
def plotIntensity ():
	
	# index boundaries for time 3D plot
	nStart = 0
	nEnd   = 10000
	
	with open("time_series.txt", "r") as file:
		lines = file.readlines()
	
	time           = []
	intensity      = []
	rho_GS_e_act   = []
	rho_GS_h_act   = []
	rho_GS_e_inact = []
	rho_GS_h_inact = []
	rho_ES_e       = []
	rho_ES_h       = []
	E_real         = []
	E_imag         = []
	
	for line in lines:
		time.append(float((line.split('	')[0])))
		intensity.append(float((line.split('	')[1])))
		E_real.append(float((line.split('	')[2])))
		E_imag.append(float((line.split('	')[3])))
		rho_GS_e_act.append(float((line.split('	')[6])))
		rho_GS_h_act.append(float((line.split('	')[7])))
		rho_GS_e_inact.append(float((line.split('	')[8])))
		rho_GS_h_inact.append(float((line.split('	')[9])))
		rho_ES_e.append(float((line.split('	')[10])))
		rho_ES_h.append(float((line.split('	')[11])))
	
	time           = np.array(time)
	intensity      = np.array(intensity)
	E_real         = np.array(E_real)
	E_imag         = np.array(E_imag)
	rho_GS_e_act   = np.array(rho_GS_e_act)
	rho_GS_h_act   = np.array(rho_GS_h_act)
	rho_GS_e_inact = np.array(rho_GS_e_inact)
	rho_GS_h_inact = np.array(rho_GS_h_inact)
	rho_ES_e       = np.array(rho_ES_e)
	rho_ES_h       = np.array(rho_ES_h)
	
	# calculation of inversion
	inversion_GS_act   = rho_GS_e_act   + rho_GS_h_act   - 1.0
	inversion_GS_inact = rho_GS_e_inact + rho_GS_h_inact - 1.0
	inversion_ES       = rho_ES_e       + rho_ES_h       - 1.0
	
	fig, (ax1, ax2) = plt.subplots(1, 2) #sharey=True
	ax12 = ax1.twinx()
	fig.set_size_inches(5.9, 3.2)
	plt.rcParams.update({"font.size": 9})
	fig.subplots_adjust(wspace=0.7, top=0.99, bottom=0.22, left=0.08, right=0.99)
	
	fig.text(0.005, 0.93, "a)")
	ax1.plot(time[nStart:nEnd], intensity[nStart:nEnd], color="crimson")
	ax1.set_xlabel(r"time $t$ / ns", size=9.0)
	ax1.set_ylabel(r"intensity $|E|^2$", color="crimson", size=9.0)
	ax1.set_ylim(np.min(intensity) - 0.1, np.max(intensity) + 0.3)
	ax1.set_xticks([0.0, 5.0, 10.0])
	ax1.set_yticks([0.0, 1.0, 2.0, 3.0])
	ax1.tick_params(axis="x", labelsize=9.0)
	ax1.tick_params(axis="y", colors="crimson", labelsize=9.0)
	ax1.set_zorder(1)
	ax1.set_facecolor("none")
	
	ax12.plot(time[nStart:nEnd], inversion_GS_act[nStart:nEnd], color="orange", label="GS act")
	ax12.plot(time[nStart:nEnd], inversion_GS_inact[nStart:nEnd], color="gray", linestyle="--", label="GS inact")
	ax12.plot(time[nStart:nEnd], inversion_ES[nStart:nEnd], color="cornflowerblue", label="ES")
	ax12.set_ylabel(r"population inversion" + "\n" + r"$\rho_{m,e}^{(in)act} + \rho_{m,h}^{(in)act} - 1$", size=9.0)
	ax12.set_ylim(-1.075, 1.075)
	ax12.set_yticks([-1.0, 0.0, 1.0])
	ax12.tick_params(axis="y", labelsize=9.0)
	ax12.set_zorder(2)
	
	ax12.legend(bbox_to_anchor=(0.44, 0.33))
	
	# ~ fig, ax = plt.subplots()
	# ~ fig.set_size_inches(5.9, 4.8)
	# ~ fig.subplots_adjust(top=0.99, bottom=0.15, left=0.10, right=0.99)
	
	fig.text(0.575, 0.93, "b)")
	ax2.plot(inversion_GS_act, intensity, color="orange", label="GS act")
	ax2.plot(inversion_GS_inact, intensity, color="gray", linestyle="--", label="GS inact")
	ax2.plot(inversion_ES, intensity, color="cornflowerblue", label="ES")
	ax2.set_xlabel(r"population inversion" + "\n" + r"$\rho_{m,e}^{(in)act} + \rho_{m,h}^{(in)act} - 1$", size=9.0)
	ax2.set_ylabel(r"intensity $|E|^2$", color="crimson", size=9.0)
	ax2.set_xlim(-1.075, 1.075)
	ax2.set_ylim(-0.15, 3.15)
	ax2.set_xticks([-1.0, 0.0, 1.0])
	ax2.set_yticks([0.0, 1.0, 2.0, 3.0])
	ax2.tick_params(axis="x", labelsize=9.0)
	ax2.tick_params(axis="y", colors="crimson", labelsize=9.0)
	ax2.grid(color="lightgray")
	
	ax2.legend(loc="upper left")
	
	plt.show()


plotIntensity()

