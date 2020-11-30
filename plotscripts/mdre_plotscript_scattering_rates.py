import numpy as np
import matplotlib.pyplot as plt


# plots relaxation scattering rates time series for MDRE model
def plotScatteringRates ():
	
	with open("time_series_mdre.txt", "r") as file:
		lines = file.readlines()
	
	time         = []
	intensity    = []
	w_e          = []
	w_h          = []
	S_e_rel_diff = []
	S_h_rel_diff = []
	
	for line in lines:
		time.append(float((line.split('	')[0])))
		intensity.append(float((line.split('	')[1])))
		w_e.append(float((line.split('	')[10])))
		w_h.append(float((line.split('	')[11])))
		S_e_rel_diff.append(float((line.split('	')[18])))
		S_h_rel_diff.append(float((line.split('	')[19])))
	
	y_shift1 = 5.75
	y_shift2 = 4.0
	y_shift3 = 1.0
	y_shift4 = 290.0
	
	time         = np.array(time)
	intensity    = np.array(intensity) + y_shift1
	w_e          = np.array(w_e) + y_shift2
	w_h          = np.array(w_h) + y_shift3
	S_e_rel_diff = np.array(S_e_rel_diff) + y_shift4
	S_h_rel_diff = np.array(S_h_rel_diff)
	
	
	fig, ax1 = plt.subplots()
	fig.set_size_inches(10.0, 7.0)
	fig.suptitle("MDRE relaxation scattering rates time series", fontsize=14)
	
	ax1.set_title(r"@ $K=0.06$, $\Delta\nu_{inj}=0.95$ GHz")
	p1 = ax1.plot(time, intensity, color="red", alpha=0.5, label=r"intensity $|E|^2$ + " + str(y_shift1))[0]
	p2 = ax1.plot(time, w_e, color="green", alpha=0.5, label=r"2D electron density $w_e$ + " + str(y_shift2))[0]
	p3 = ax1.plot(time, w_h, color="green", alpha=0.5, linestyle="--", label=r"2D hole density $w_h$ + " + str(y_shift3))[0]
	ax1.set_xlabel(r"time $t$ / ns")
	ax1.set_ylabel(r"intensity and 2D charge carrier density")
	
	ax2 = ax1.twinx()
	p4  = ax2.plot(time, S_e_rel_diff, color="red", label=r"$S_e^{rel,in}(w_e) - S_e^{rel,out}(w_e)$ + " + str(y_shift4))[0]
	p5  = ax2.plot(time, S_h_rel_diff, color="blue", label=r"$S_h^{rel,in}(w_h) - S_h^{rel,out}(w_h)$")[0]
	ax2.set_ylabel(r"relaxation scattering rates $S^{rel,in} - S^{rel,out}$")
	
	plots = [p1, p2, p3, p4, p5]
	ax2.legend(plots, [p.get_label() for p in plots], loc="lower right", fancybox=True, framealpha=0.85)
	
	plt.show()


plotScatteringRates()
