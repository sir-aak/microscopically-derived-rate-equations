import numpy as np
import matplotlib.pyplot as plt


# parameters
K_inj  = 0
DnuInj = 0
K_fb   = 0
phi    = 0
tau    = 0


# plots phase space for MDRE model 
# intensity versus electron density and hole density
def plotPhaseSpaceOccupationProbabilities ():
	
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
	
	for line in lines:
		time.append(float((line.split('	')[0])))
		intensity.append(float((line.split('	')[1])))
		rho_GS_e_act.append(float((line.split('	')[6])))
		rho_GS_h_act.append(float((line.split('	')[7])))
		rho_GS_e_inact.append(float((line.split('	')[8])))
		rho_GS_h_inact.append(float((line.split('	')[9])))
		rho_ES_e.append(float((line.split('	')[10])))
		rho_ES_h.append(float((line.split('	')[11])))
	
	time           = np.array(time)
	intensity      = np.array(intensity)
	rho_GS_e_act   = np.array(rho_GS_e_act)
	rho_GS_h_act   = np.array(rho_GS_h_act)
	rho_GS_e_inact = np.array(rho_GS_e_inact)
	rho_GS_h_inact = np.array(rho_GS_h_inact)
	rho_ES_e       = np.array(rho_ES_e)
	rho_ES_h       = np.array(rho_ES_h)
	
	t_start = time[0]
	t_end   = time[-1]
	
	plt.figure(figsize=(7, 5))
	plt.suptitle("MDRE: phase space: occupation probabilities")
	plt.title(r"$K_{inj}=$" + str(K_inj) + r", $\Delta\nu_{inj}=$" + str(DnuInj) + r", $K_{fb}=$" + str(K_fb) + r", $\varphi=$" + str(phi) + r", $\tau=$" + str(tau) + r"$\,$ns")
	plt.plot(rho_GS_e_act, intensity, color="blue", alpha=0.75, label=r"$\rho_{GS,\,e}^{act}$")
	plt.plot(rho_GS_h_act, intensity, color="blue", linestyle="--", alpha=0.75, label=r"$\rho_{GS,\,h}^{act}$")
	plt.plot(rho_GS_e_inact, intensity, color="gray", alpha=0.75, label=r"$\rho_{GS,\,e}^{inact}$")
	plt.plot(rho_GS_h_inact, intensity, color="gray", linestyle="--", alpha=0.75, label=r"$\rho_{GS,\,h}^{inact}$")
	plt.plot(rho_ES_e, intensity, color="orange", alpha=0.75, label=r"$\rho_{ES,\,e}$")
	plt.plot(rho_ES_h, intensity, color="orange", linestyle="--", alpha=0.75, label=r"$\rho_{ES,\,h}$")
	
	# marks for start and end
	plt.plot(rho_GS_e_act[0], intensity[0], color="black", marker="o", markersize=5, zorder=2)
	plt.plot(rho_GS_e_act[-1], intensity[-1], color="blue", marker="o", markersize=5, zorder=2)
	plt.plot(rho_GS_h_act[-1], intensity[-1], color="blue", marker="o", markersize=5, zorder=2)
	plt.plot(rho_GS_e_inact[-1], intensity[-1], color="gray", marker="o", markersize=5, zorder=2)
	plt.plot(rho_GS_h_inact[-1], intensity[-1], color="gray", marker="o", markersize=5, zorder=2)
	plt.plot(rho_ES_e[-1], intensity[-1], color="orange", marker="o", markersize=5, zorder=2)
	plt.plot(rho_ES_h[-1], intensity[-1], color="orange", marker="o", markersize=5, zorder=2)
	
	# labeled marks for start and end
	# ~ plt.plot(rho_GS_e_act[0], intensity[0], color="black", marker="o", label=r"$\rho_{m\,b}^{in(act)}(t=$" + str(t_start) + r"$\,$ns)", markersize=5, zorder=2)
	# ~ plt.plot(rho_GS_e_act[-1], intensity[-1], color="blue", marker="o", label=r"$\rho_{GS,\,e}^{act}(t=$" + str(t_end) + r"$\,$ns)", markersize=5, zorder=2)
	# ~ plt.plot(rho_GS_h_act[-1], intensity[-1], color="blue", marker="o", label=r"$\rho_{GS,\,h}^{act}(t=$" + str(t_end) + r"$\,$ns)", markersize=5, zorder=2)
	# ~ plt.plot(rho_GS_e_inact[-1], intensity[-1], color="gray", marker="o", label=r"$\rho_{GS,\,e}^{inact}(t=$" + str(t_end) + r"$\,$ns)", markersize=5, zorder=2)
	# ~ plt.plot(rho_GS_h_inact[-1], intensity[-1], color="gray", marker="o", label=r"$\rho_{GS,\,h}^{inact}(t=$" + str(t_end) + r"$\,$ns)", markersize=5, zorder=2)
	# ~ plt.plot(rho_ES_e[-1], intensity[-1], color="orange", marker="o", label=r"$\rho_{ES,\,e}(t=$" + str(t_end) + r"$\,$ns)", markersize=5, zorder=2)
	# ~ plt.plot(rho_ES_h[-1], intensity[-1], color="orange", marker="o", label=r"$\rho_{ES,\,h}(t=$" + str(t_end) + r"$\,$ns)", markersize=5, zorder=2)
	
	plt.xlabel(r"occupation probabilities $\rho_{m\,b}^{in(act)}$")
	plt.ylabel(r"intensity $|E|^2$")
	plt.grid(color="lightgray")
	plt.legend(loc="upper left")
	
	plt.show()


# plots phase space for MDRE model 
# intensity versus electron density and hole density
def plotPhaseSpaceCarrierDensities ():
	
	with open("time_series.txt", "r") as file:
		lines = file.readlines()
	
	time      = []
	intensity = []
	w_e       = []
	w_h       = []
	
	for line in lines:
		time.append(float((line.split('	')[0])))
		intensity.append(float((line.split('	')[1])))
		w_e.append(float((line.split('	')[10])))
		w_h.append(float((line.split('	')[11])))
	
	time      = np.array(time)
	intensity = np.array(intensity)
	w_e       = np.array(w_e)
	w_h       = np.array(w_h)
	
	t_start = time[0]
	t_end   = time[-1]
	
	plt.figure(figsize=(7, 5))
	# ~ plt.suptitle("MDRE: phase space: carrier densities")
	# ~ plt.title(r"$K_{inj}=$" + str(K_inj) + r", $\Delta\nu_{inj}=$" + str(DnuInj) + r", $K_{fb}=$" + str(K_fb) + r", $\varphi=$" + str(phi) + r", $\tau=$" + str(tau) + r"$\,$ns")
	plt.plot(w_e, intensity, color="green", alpha=0.75, label=r"$w_e$")
	plt.plot(w_h, intensity, color="green", linestyle="--", alpha=0.75, label=r"$w_h$")
	
	
	# marks for start and end
	plt.plot(w_e[0], intensity[0], color="black", marker="o", markersize=6, zorder=2)
	plt.plot(w_e[-1], intensity[-1], color="green", marker="o", markersize=6, zorder=2)
	plt.plot(w_h[-1], intensity[-1], color="green", marker="o", markersize=6, zorder=2)
	
	# labeled marks for start and end
	# ~ plt.plot(w_e[0], intensity[0], color="black", marker="o", label=r"$w_b(t=$" + str(t_start) + r"$\,$ns)", markersize=6, zorder=2)
	# ~ plt.plot(w_e[-1], intensity[-1], color="green", marker="o", label=r"$w_e(t=$" + str(t_end) + r"$\,$ns)", markersize=6, zorder=2)
	# ~ plt.plot(w_h[-1], intensity[-1], color="green", marker="o", label=r"$w_h(t=550\,$ns)", markersize=6, zorder=2)
	
	plt.xlabel(r"charge carrier densities $w_b$")
	plt.ylabel(r"intensity $|E|^2$")
	plt.grid(color="lightgray")
	plt.legend(loc="upper left")
	
	plt.show()


plotPhaseSpaceOccupationProbabilities()
plotPhaseSpaceCarrierDensities()
